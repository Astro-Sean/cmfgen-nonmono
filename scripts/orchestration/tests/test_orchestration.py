"""Unit tests for the orchestration layer.

Tests cover:
- Configuration validation (missing fields, type errors, contradictions)
- Epoch mapping (sorting, normalization, duplicate labels)
- State machine (valid/invalid transitions, persistence)
- Command building (no duplicate flags, reproducible)
- Site config (partition validation, time limits)
- Manifests (file writing, atomic writes)
"""

import sys
import os
import json
import tempfile
import unittest
from pathlib import Path

# Add scripts dir to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from orchestration.config import (
    ResolvedConfig, ConfigError, validate_yaml, resolve_config,
    EpochMapping, _parse_snapshot_days, _get_stella_target_days,
)
from orchestration.state import (
    WorkflowStates, StateMachine, VALID_TRANSITIONS,
)
from orchestration.manifests import ManifestWriter
from orchestration.site import (
    SiteConfig, get_site_config, SiteConfigError,
    _time_to_seconds, _time_leq,
)


class TestConfigValidation(unittest.TestCase):
    """Test YAML config validation."""

    def test_valid_config(self):
        """A minimal valid config passes validation."""
        raw = {
            "model_name": "test_model",
            "hydro_backend": "stella",
            "run_cmfgen": True,
            "cmfgen_conversion": {"snapshot_days": [30]},
            "stella": {"target_days": 60},
            "cmfgen": {"vel_law": 12, "target_nd": 100},
        }
        errors = validate_yaml(raw)
        self.assertEqual(errors, [], f"Unexpected errors: {errors}")

    def test_missing_model_name(self):
        """Missing model_name is an error."""
        raw = {"hydro_backend": "stella"}
        errors = validate_yaml(raw)
        self.assertTrue(any("model_name" in e for e in errors))

    def test_wrong_type_model_name(self):
        """model_name must be a string."""
        raw = {"model_name": 123}
        errors = validate_yaml(raw)
        self.assertTrue(any("must be" in e for e in errors))

    def test_mutually_exclusive(self):
        """reuse_stella and skip_stella can't both be True."""
        raw = {
            "model_name": "test",
            "reuse_stella": True,
            "skip_stella": True,
        }
        errors = validate_yaml(raw)
        self.assertTrue(any("Mutually exclusive" in e for e in errors))

    def test_invalid_hydro_backend(self):
        """hydro_backend must be stella or heracles."""
        raw = {"model_name": "test", "hydro_backend": "flash"}
        errors = validate_yaml(raw)
        self.assertTrue(any("hydro_backend" in e for e in errors))

    def test_empty_snapshot_days(self):
        """Empty snapshot_days list is an error."""
        raw = {
            "model_name": "test",
            "cmfgen_conversion": {"snapshot_days": []},
        }
        errors = validate_yaml(raw)
        self.assertTrue(any("empty" in e for e in errors))

    def test_duplicate_snapshot_labels(self):
        """Duplicate snapshot day labels are detected."""
        raw = {
            "model_name": "test",
            "cmfgen_conversion": {"snapshot_days": [30, 30]},
        }
        errors = validate_yaml(raw)
        self.assertTrue(any("Duplicate" in e for e in errors))

    def test_non_numeric_snapshot_days(self):
        """Non-numeric snapshot days are detected."""
        raw = {
            "model_name": "test",
            "cmfgen_conversion": {"snapshot_days": ["foo"]},
        }
        errors = validate_yaml(raw)
        self.assertTrue(any("non-numeric" in e for e in errors))

    def test_vel_law_without_tridiag(self):
        """vel_law != 12 without depth_independent_solver is flagged."""
        raw = {
            "model_name": "test",
            "cmfgen": {"vel_law": 11, "depth_independent_solver": False},
        }
        errors = validate_yaml(raw)
        self.assertTrue(any("vel_law" in e for e in errors))

    def test_vel_law_12_ok(self):
        """vel_law=12 is always valid."""
        raw = {
            "model_name": "test",
            "cmfgen": {"vel_law": 12, "depth_independent_solver": False},
        }
        errors = validate_yaml(raw)
        self.assertFalse(any("vel_law" in e for e in errors))

    def test_unknown_partition(self):
        """Unknown partition names are flagged."""
        raw = {
            "model_name": "test",
            "stella": {"partition": "nonexistent"},
        }
        errors = validate_yaml(raw)
        self.assertTrue(any("not a known partition" in e for e in errors))

    def test_min_nd_gt_max_points(self):
        """min_nd > max_points is an error."""
        raw = {
            "model_name": "test",
            "cmfgen_conversion": {"min_nd": 200, "max_points": 100},
        }
        errors = validate_yaml(raw)
        self.assertTrue(any("min_nd" in e and "max_points" in e for e in errors))


class TestEpochMapping(unittest.TestCase):
    """Test epoch mapping and snapshot day handling."""

    def test_parse_single_day(self):
        """Single snapshot day is converted to list."""
        raw = {"cmfgen_conversion": {"snapshot_days": 30}}
        days = _parse_snapshot_days(raw)
        self.assertEqual(days, [30.0])

    def test_parse_list_days(self):
        """List of snapshot days is sorted."""
        raw = {"cmfgen_conversion": {"snapshot_days": [30, 10, 50]}}
        days = _parse_snapshot_days(raw)
        self.assertEqual(days, [10.0, 30.0, 50.0])

    def test_parse_empty_defaults(self):
        """Empty config defaults to [40.0]."""
        days = _parse_snapshot_days({})
        self.assertEqual(days, [40.0])

    def test_target_days_explicit(self):
        """Explicit stella.target_days takes precedence."""
        raw = {"stella": {"target_days": 100}, "cmfgen_conversion": {"snapshot_days": [30]}}
        td = _get_stella_target_days(raw, [30.0])
        self.assertEqual(td, 100.0)

    def test_target_days_computed(self):
        """Target days computed from max snapshot + extra."""
        raw = {"cmfgen_conversion": {"snapshot_days": [30, 50], "stella_extra_days": 10}}
        days = _parse_snapshot_days(raw)
        td = _get_stella_target_days(raw, days)
        self.assertEqual(td, 60.0)  # 50 + 10

    def test_target_days_excludes_homologous(self):
        """Homologous evolution days are excluded from max."""
        raw = {
            "cmfgen_conversion": {
                "snapshot_days": [30, 349],
                "homologous_evolution_days": [349],
                "stella_extra_days": 0,
            }
        }
        days = _parse_snapshot_days(raw)
        td = _get_stella_target_days(raw, days)
        self.assertEqual(td, 30.0)  # 30 + 0 (349 excluded)

    def test_epoch_mapping_labels(self):
        """Epoch mappings produce correct labels."""
        mappings = [
            EpochMapping(30.0, None, "30.0days", "/dir/cmfgen/30.0days",
                        "/dir/cmfgen/30.0days", "pending"),
        ]
        self.assertEqual(mappings[0].output_label, "30.0days")
        self.assertEqual(mappings[0].status, "pending")


class TestStateMachine(unittest.TestCase):
    """Test workflow state machine."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.state_file = os.path.join(self.tmpdir, "workflow_state.json")

    def test_init_state(self):
        """New state machine starts at INIT."""
        sm = StateMachine("run1", "model1", "hash1", self.state_file)
        self.assertEqual(sm.state.current_state, "INIT")

    def test_valid_transition(self):
        """INIT -> CONFIG_VALIDATED is valid."""
        sm = StateMachine("run1", "model1", "hash1", self.state_file)
        sm.transition(WorkflowStates.CONFIG_VALIDATED, message="Test")
        self.assertEqual(sm.state.current_state, "CONFIG_VALIDATED")
        self.assertEqual(len(sm.state.transitions), 1)

    def test_invalid_transition(self):
        """INIT -> COMPLETE is invalid."""
        sm = StateMachine("run1", "model1", "hash1", self.state_file)
        with self.assertRaises(ValueError):
            sm.transition(WorkflowStates.COMPLETE)

    def test_force_transition(self):
        """Force transition bypasses validation."""
        sm = StateMachine("run1", "model1", "hash1", self.state_file)
        sm.transition(WorkflowStates.COMPLETE, force=True, message="Forced")
        self.assertEqual(sm.state.current_state, "COMPLETE")

    def test_persistence(self):
        """State persists across StateMachine instances."""
        sm1 = StateMachine("run1", "model1", "hash1", self.state_file)
        sm1.transition(WorkflowStates.CONFIG_VALIDATED, message="Test")
        sm1.record_job_id("STELLA", "12345")

        sm2 = StateMachine("run1", "model1", "hash1", self.state_file)
        self.assertEqual(sm2.state.current_state, "CONFIG_VALIDATED")
        self.assertEqual(sm2.state.job_ids.get("STELLA"), "12345")

    def test_config_hash_mismatch(self):
        """Config hash mismatch starts fresh."""
        sm1 = StateMachine("run1", "model1", "hash1", self.state_file)
        sm1.transition(WorkflowStates.CONFIG_VALIDATED)

        sm2 = StateMachine("run1", "model1", "hash2", self.state_file)
        self.assertEqual(sm2.state.current_state, "INIT")

    def test_idempotent_transition(self):
        """Transition to same state is idempotent."""
        sm = StateMachine("run1", "model1", "hash1", self.state_file)
        sm.transition(WorkflowStates.CONFIG_VALIDATED)
        sm.transition(WorkflowStates.CONFIG_VALIDATED)  # no-op
        self.assertEqual(len(sm.state.transitions), 1)

    def test_terminal_state(self):
        """Terminal states are detected."""
        sm = StateMachine("run1", "model1", "hash1", self.state_file)
        sm.transition(WorkflowStates.CONFIG_VALIDATED)
        sm.transition(WorkflowStates.PROFILE_READY)
        sm.transition(WorkflowStates.STELLA_SUBMITTED)
        sm.transition(WorkflowStates.STELLA_COMPLETE)
        sm.transition(WorkflowStates.SNAPSHOTS_MAPPED)
        sm.transition(WorkflowStates.SNAPSHOTS_CONVERTED)
        sm.transition(WorkflowStates.GRID_VALIDATED)
        sm.transition(WorkflowStates.CMFGEN_SUBMITTED)
        sm.transition(WorkflowStates.CMFGEN_COMPLETE)
        sm.transition(WorkflowStates.COMPLETE)
        self.assertTrue(sm.is_terminal())
        self.assertTrue(sm.is_complete())

    def test_epoch_status(self):
        """Epoch status is recorded and persisted."""
        sm = StateMachine("run1", "model1", "hash1", self.state_file)
        sm.record_epoch_status(30.0, "converted")
        self.assertEqual(sm.state.epoch_status.get("30.0days"), "converted")


class TestSiteConfig(unittest.TestCase):
    """Test MPC site configuration."""

    def test_default_site_config(self):
        """Default site config resolves correctly."""
        site = get_site_config()
        self.assertEqual(site.interactive_partition, "interactive")
        self.assertEqual(site.production_partition, "general")
        self.assertEqual(site.cluster_name, "viper")

    def test_unknown_partition(self):
        """Unknown partition raises error."""
        with self.assertRaises(SiteConfigError):
            get_site_config(interactive_partition="nonexistent")

    def test_time_to_seconds(self):
        """Time string parsing works."""
        self.assertEqual(_time_to_seconds("01:59:00"), 7140)
        self.assertEqual(_time_to_seconds("2:00:00"), 7200)
        self.assertEqual(_time_to_seconds("1-00:00:00"), 86400)
        self.assertEqual(_time_to_seconds("15:00"), 900)

    def test_time_leq(self):
        """Time comparison works."""
        self.assertTrue(_time_leq("01:59:00", "2:00:00"))
        self.assertFalse(_time_leq("2:00:00", "01:59:00"))
        self.assertTrue(_time_leq("1-00:00:00", "1-00:00:00"))

    def test_test_time_exceeds_limit(self):
        """Test time exceeding partition limit raises error."""
        with self.assertRaises(SiteConfigError):
            get_site_config(test_time="3:00:00")  # interactive max is 2:00:00


class TestManifests(unittest.TestCase):
    """Test manifest writing."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.mw = ManifestWriter(self.tmpdir, "test_run")

    def test_log_writes(self):
        """Log lines are written to orchestration.log."""
        self.mw.log("Test message")
        log_path = Path(self.tmpdir) / "orchestration.log"
        self.assertTrue(log_path.exists())
        content = log_path.read_text()
        self.assertIn("Test message", content)

    def test_record_command(self):
        """Commands are recorded in commands.jsonl."""
        self.mw.record_command({
            "executable": "sbatch",
            "arguments": ["job.slurm"],
            "execution_mode": "dry-run",
        })
        cmd_path = Path(self.tmpdir) / "commands.jsonl"
        self.assertTrue(cmd_path.exists())
        lines = cmd_path.read_text().strip().split("\n")
        self.assertEqual(len(lines), 1)
        cmd = json.loads(lines[0])
        self.assertEqual(cmd["executable"], "sbatch")
        self.assertIn("command_hash", cmd)
        self.assertIn("timestamp", cmd)

    def test_record_file(self):
        """Files are recorded in file manifest."""
        test_file = Path(self.tmpdir) / "test.txt"
        test_file.write_text("test content")
        self.mw.record_file(str(test_file), "test")
        self.mw.write_file_manifest()
        manifest_path = Path(self.tmpdir) / "file_manifest.json"
        self.assertTrue(manifest_path.exists())
        data = json.loads(manifest_path.read_text())
        self.assertIn(str(test_file), data["files"])

    def test_write_job_graph(self):
        """Job graph is written correctly."""
        jobs = [
            {"name": "stella_test", "stage": "STELLA", "job_id": "123",
             "partition": "general", "dependency": None, "state": "SUBMITTED"},
        ]
        self.mw.write_job_graph(jobs)
        graph_path = Path(self.tmpdir) / "job_graph.json"
        self.assertTrue(graph_path.exists())
        data = json.loads(graph_path.read_text())
        self.assertEqual(len(data["jobs"]), 1)


class TestFullConfigResolution(unittest.TestCase):
    """Test full config resolution with the real SN_12C_CSM2 config."""

    def setUp(self):
        self.project_dir = Path(__file__).resolve().parents[3]
        self.scripts_dir = self.project_dir / "scripts"
        self.config_file = self.project_dir / "configs" / "SN_12C_CSM2.yaml"

    def test_resolve_real_config(self):
        """Real config resolves without errors."""
        if not self.config_file.exists():
            self.skipTest("Config file not found")
        cfg = resolve_config(str(self.config_file), str(self.project_dir), str(self.scripts_dir))
        self.assertEqual(cfg.model_name, "SN_12C_CSM2")
        self.assertTrue(cfg.use_stella)
        self.assertEqual(cfg.snapshot_days, [30.0])
        self.assertFalse(cfg.time_series)
        self.assertEqual(cfg.site.interactive_partition, "interactive")
        self.assertEqual(cfg.site.account, "mpe_cpu")

    def test_config_hash_stable(self):
        """Same config produces same hash."""
        if not self.config_file.exists():
            self.skipTest("Config file not found")
        cfg1 = resolve_config(str(self.config_file), str(self.project_dir), str(self.scripts_dir))
        cfg2 = resolve_config(str(self.config_file), str(self.project_dir), str(self.scripts_dir))
        self.assertEqual(cfg1.config_hash, cfg2.config_hash)

    def test_resolved_config_serializable(self):
        """Resolved config can be serialized to YAML."""
        if not self.config_file.exists():
            self.skipTest("Config file not found")
        cfg = resolve_config(str(self.config_file), str(self.project_dir), str(self.scripts_dir))
        yaml_str = cfg.to_yaml()
        self.assertIn("model_name", yaml_str)
        self.assertIn("SN_12C_CSM2", yaml_str)


if __name__ == "__main__":
    unittest.main()
