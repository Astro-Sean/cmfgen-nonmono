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
from unittest.mock import patch

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
from generate_sn_b_dessart_grid import make_observational_timing
from run_model_grid import (
    apply_nonmono_mem_floor,
    estimate_nonmono_mem_gb,
    check_existing_jobs,
    get_job_resources,
    get_time_limit,
    load_yaml_config,
    resolve_submission,
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

    def _dessart_timing_config(self, source_age, timing=None):
        raw = {
            "model_name": "DX26_timing_test",
            "output": {"model_dir": "./SN_B_dessart_grid/DX26_timing_test"},
            "x_model_from_sn": {
                "source_age_days": source_age,
                "t_start_days": 2.0,
            },
        }
        if timing is not None:
            raw["observational_timing"] = timing
        return raw

    def test_future_dessart_accepts_observed_age_without_metadata(self):
        errors = validate_yaml(self._dessart_timing_config(334.0))
        self.assertEqual(errors, [])

    def test_future_dessart_rejects_unmarked_age_mismatch(self):
        errors = validate_yaml(self._dessart_timing_config(350.0))
        self.assertTrue(any("differs from the observed SN_A age" in e for e in errors))

    def test_future_dessart_accepts_observed_timing(self):
        timing = {
            "sn_a_launch_mjd": 60490.0,
            "sn_b_t0_mjd": 60824.0,
            "observed_source_age_days": 334.0,
            "configured_source_age_days": 334.0,
            "profile_age_days": 336.0,
            "source_age_policy": "observed",
        }
        errors = validate_yaml(self._dessart_timing_config(334.0, timing))
        self.assertEqual(errors, [])

    def test_future_dessart_sensitivity_requires_reason(self):
        timing = {
            "sn_a_launch_mjd": 60490.0,
            "sn_b_t0_mjd": 60824.0,
            "observed_source_age_days": 334.0,
            "configured_source_age_days": 350.0,
            "profile_age_days": 352.0,
            "source_age_policy": "sensitivity",
        }
        errors = validate_yaml(self._dessart_timing_config(350.0, timing))
        self.assertTrue(any("sensitivity_reason is required" in e for e in errors))
        timing["sensitivity_reason"] = "Intentional age scan"
        self.assertEqual(validate_yaml(self._dessart_timing_config(350.0, timing)), [])

    def test_future_dessart_rejects_changed_observation_mjd(self):
        timing = {
            "sn_a_launch_mjd": 60491.0,
            "sn_b_t0_mjd": 60824.0,
            "observed_source_age_days": 333.0,
            "configured_source_age_days": 333.0,
            "profile_age_days": 335.0,
            "source_age_policy": "observed",
        }
        errors = validate_yaml(self._dessart_timing_config(333.0, timing))
        self.assertTrue(any("sn_a_launch_mjd" in e for e in errors))

    def test_runtime_yaml_loader_rejects_off_epoch_future_model(self):
        raw = self._dessart_timing_config(350.0)
        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml") as handle:
            handle.write(json.dumps(raw))
            handle.flush()
            with self.assertRaisesRegex(ValueError, "observed SN_A age"):
                load_yaml_config(handle.name)


class TestDessartTimingGeneration(unittest.TestCase):
    def test_observation_age_is_the_default(self):
        timing = make_observational_timing("DX26_default", 334.0, {})
        self.assertEqual(timing["source_age_policy"], "observed")
        self.assertEqual(timing["observed_source_age_days"], 334.0)

    def test_future_off_epoch_point_requires_explicit_policy(self):
        with self.assertRaises(ValueError):
            make_observational_timing("DX26_off_epoch", 350.0, {})

    def test_age_sensitivity_requires_a_reason(self):
        timing = make_observational_timing(
            "DX26_age_scan",
            350.0,
            {
                "source_age_policy": "sensitivity",
                "source_age_reason": "Check sensitivity to the event separation.",
            },
        )
        self.assertEqual(timing["source_age_policy"], "sensitivity")
        self.assertIn("sensitivity_reason", timing)


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


class TestRunModelGridPartitionPolicy(unittest.TestCase):
    def test_interactive_partition_sets_interactive_resources_and_walltime(self):
        resources = get_job_resources(
            False,
            config_partition="interactive",
            config_mem="40G",
            config_cpus=8,
        )
        self.assertEqual(resources, ("interactive", 8, "30G"))
        self.assertEqual(
            get_time_limit(False, partition="interactive"), "01:59:00"
        )

    def test_small_partition_sets_production_walltime(self):
        self.assertEqual(get_time_limit(False, partition="small"), "23:59:00")

    def test_debug_submission_resolves_to_interactive(self):
        config = {"submission": "debug", "cmfgen": {"partition": "small"}}
        debug, partition = resolve_submission(config)
        self.assertTrue(debug)
        self.assertEqual(partition, "interactive")
        self.assertEqual(config["cmfgen"], {})
        self.assertEqual(get_time_limit(debug, partition=partition), "01:59:00")

    def test_sob_skips_knop_memory_estimate(self):
        cmfgen_cfg = {"non_mono_mode": "SOB"}
        estimate, reason = estimate_nonmono_mem_gb(
            Path("/missing"), cmfgen_cfg
        )
        self.assertIsNone(estimate)
        self.assertIn("KNOP", reason.upper())
        self.assertEqual(
            apply_nonmono_mem_floor(
                "30G", Path("/missing"), cmfgen_cfg, debug=True
            ),
            "30G",
        )

    def test_convergence_monitor_matches_model_cancellation(self):
        result = type("SqueueResult", (), {
            "returncode": 0,
            "stdout": "12112714 dessart_2.82days_convergence\n",
        })()
        with patch("run_model_grid.subprocess.run", return_value=result):
            jobs = check_existing_jobs("dessart")
        self.assertEqual(
            jobs, [("12112714", "dessart_2.82days_convergence")]
        )


class TestSiteConfig(unittest.TestCase):
    """Test MPC site configuration."""

    def test_default_site_config(self):
        """Default site config resolves correctly."""
        site = get_site_config()
        self.assertEqual(site.interactive_partition, "interactive")
        self.assertEqual(site.production_partition, "general")
        # cluster_name is detected from scontrol; on Viper it is "viper",
        # on Raven it is "raven".  Just verify it is a known MPCDF cluster.
        self.assertIn(site.cluster_name, ("viper", "raven"))

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
