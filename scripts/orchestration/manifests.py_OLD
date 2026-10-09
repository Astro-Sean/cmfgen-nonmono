"""Manifest and report writers.

Writes all required observability files:
- run_manifest.json
- resolved_config.yaml
- workflow_state.json (via StateMachine)
- epoch_map.json
- job_graph.json
- commands.jsonl
- file_manifest.json
- validation_summary.json
- interactive_monitor.jsonl (via InteractiveTestRunner)
- interactive_run_report.txt (via InteractiveTestRunner)
- orchestration.log
"""

from __future__ import annotations

import hashlib
import json
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any, Optional, Optional

import yaml

from .config import ResolvedConfig, EpochMapping
from .state import WorkflowState, StateMachine


class ManifestWriter:
    """Writes all manifest and report files to a run directory.

    All files are written atomically (write to temp, then rename).
    """

    def __init__(self, run_dir: str, run_id: str):
        self.run_dir = Path(run_dir)
        self.run_dir.mkdir(parents=True, exist_ok=True)
        self.run_id = run_id
        self.commands: list[dict[str, Any]] = []
        self.files: dict[str, dict[str, Any]] = {}
        self.log_lines: list[str] = []

    def _atomic_write(self, path: Path, content: str):
        """Write content atomically."""
        tmp = path.with_suffix(path.suffix + ".tmp")
        tmp.write_text(content)
        tmp.rename(path)

    def _atomic_write_json(self, path: Path, data: dict):
        """Write JSON atomically."""
        self._atomic_write(path, json.dumps(data, indent=2, default=str))

    def write_resolved_config(self, config: ResolvedConfig):
        """Write resolved_config.yaml."""
        path = self.run_dir / "resolved_config.yaml"
        self._atomic_write(path, config.to_yaml())
        self.log(f"Wrote resolved config: {path}")

    def write_run_manifest(self, config: ResolvedConfig, state: StateMachine):
        """Write run_manifest.json with full run metadata."""
        manifest = {
            "run_id": config.run_id,
            "model_name": config.model_name,
            "config_hash": config.config_hash,
            "config_path": config.config_path,
            "created_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "project_dir": config.project_dir,
            "scripts_dir": config.scripts_dir,
            "model_dir": config.model_dir,
            "hydro_backend": config.hydro_backend,
            "workflow_flags": {
                "reuse_stella": config.reuse_stella,
                "skip_stella": config.skip_stella,
                "run_cmfgen": config.run_cmfgen,
                "run_sedona": config.run_sedona,
                "run_sirocco": config.run_sirocco,
                "run_superlite": config.run_superlite,
                "debug": config.debug,
            },
            "snapshot_days": config.snapshot_days,
            "time_series": config.time_series,
            "stella_target_days": config.stella_target_days,
            "site": asdict(config.site),
            "current_state": state.state.current_state,
            "job_ids": state.state.job_ids,
        }
        path = self.run_dir / "run_manifest.json"
        self._atomic_write_json(path, manifest)
        self.log(f"Wrote run manifest: {path}")

    def write_epoch_map(self, config: ResolvedConfig):
        """Write epoch_map.json."""
        epoch_map = {
            "run_id": config.run_id,
            "model_name": config.model_name,
            "epoch_tolerance": config.epoch_tolerance,
            "mappings": [
                {
                    "requested_day": e.requested_day,
                    "stella_epoch": e.stella_epoch,
                    "output_label": e.output_label,
                    "conversion_dir": e.conversion_dir,
                    "cmfgen_dir": e.cmfgen_dir,
                    "status": e.status,
                    "delta": e.delta,
                }
                for e in config.epoch_mappings
            ],
        }
        path = self.run_dir / "epoch_map.json"
        self._atomic_write_json(path, epoch_map)
        self.log(f"Wrote epoch map: {path}")

    def write_job_graph(self, jobs: list[dict[str, Any]]):
        """Write job_graph.json describing the SLURM dependency DAG.

        Each job entry:
            {
                "name": "stella_SN_12C",
                "stage": "STELLA",
                "job_id": "12345",
                "partition": "general",
                "dependency": null | "afterok:12344",
                "state": "SUBMITTED" | "RUNNING" | "COMPLETED" | ...
            }
        """
        graph = {
            "run_id": self.run_id,
            "jobs": jobs,
            "created_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        }
        path = self.run_dir / "job_graph.json"
        self._atomic_write_json(path, graph)
        self.log(f"Wrote job graph: {path}")

    def record_command(self, command: dict[str, Any]):
        """Record an executed command. Appended to commands.jsonl.

        command dict:
            {
                "executable": "sbatch",
                "arguments": ["job.slurm"],
                "working_dir": "/path",
                "environment": {...},
                "modules": [...],
                "input_files": [...],
                "output_files": [...],
                "job_dependency": "afterok:123",
                "config_hash": "abc123",
                "command_hash": "def456",
                "execution_mode": "submit" | "inline" | "dry-run",
                "resource_request": {...},
                "timestamp": "...",
            }
        """
        command["timestamp"] = time.strftime("%Y-%m-%dT%H:%M:%S")
        # Compute command hash
        cmd_str = json.dumps(command, sort_keys=True, default=str)
        command["command_hash"] = hashlib.sha256(cmd_str.encode()).hexdigest()[:16]
        self.commands.append(command)

        # Append to commands.jsonl
        path = self.run_dir / "commands.jsonl"
        with open(path, "a") as f:
            f.write(json.dumps(command, default=str) + "\n")

    def record_file(self, path: str, file_type: str, hash_: Optional[str] = None,
                    size: Optional[int] = None):
        """Record a file in the file manifest."""
        p = Path(path)
        if size is None and p.exists():
            size = p.stat().st_size
        if hash_ is None and p.exists() and size and size < 100 * 1024 * 1024:
            h = hashlib.sha256()
            with open(p, "rb") as f:
                for chunk in iter(lambda: f.read(8192), b""):
                    h.update(chunk)
            hash_ = h.hexdigest()[:16]
        self.files[str(path)] = {
            "type": file_type,
            "size": size,
            "hash": hash_,
            "recorded_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        }

    def write_file_manifest(self):
        """Write file_manifest.json."""
        manifest = {
            "run_id": self.run_id,
            "files": self.files,
            "created_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        }
        path = self.run_dir / "file_manifest.json"
        self._atomic_write_json(path, manifest)
        self.log(f"Wrote file manifest: {path}")

    def write_validation_summary(self, summary: dict[str, Any]):
        """Write validation_summary.json."""
        summary["run_id"] = self.run_id
        summary["created_at"] = time.strftime("%Y-%m-%dT%H:%M:%S")
        path = self.run_dir / "validation_summary.json"
        self._atomic_write_json(path, summary)
        self.log(f"Wrote validation summary: {path}")

    def write_workflow_state(self, state: StateMachine):
        """Write workflow_state.json (delegates to StateMachine._save)."""
        state._save()
        self.log(f"Wrote workflow state: {state.state_file}")

    def log(self, message: str):
        """Append a line to orchestration.log."""
        ts = time.strftime("%Y-%m-%dT%H:%M:%S")
        line = f"[{ts}] {message}"
        self.log_lines.append(line)
        path = self.run_dir / "orchestration.log"
        with open(path, "a") as f:
            f.write(line + "\n")

    def write_final_report(self, success: bool, details: dict[str, Any]):
        """Write the final success/failure report."""
        report = {
            "run_id": self.run_id,
            "success": success,
            "created_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
            **details,
        }
        path = self.run_dir / "final_report.json"
        self._atomic_write_json(path, report)
        self.log(f"Wrote final report: {path} (success={success})")

    def write_all(self, config: ResolvedConfig, state: StateMachine,
                  jobs: Optional[list] = None, validation: Optional[dict] = None):
        """Write all manifests in one call."""
        self.write_resolved_config(config)
        self.write_run_manifest(config, state)
        self.write_epoch_map(config)
        if jobs:
            self.write_job_graph(jobs)
        self.write_workflow_state(state)
        self.write_file_manifest()
        if validation:
            self.write_validation_summary(validation)
