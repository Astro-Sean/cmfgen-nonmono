"""SLURM job monitoring via squeue/sacct.

Tracks job states, exit codes, resource usage, and output file growth.
Writes machine-readable event logs (interactive_monitor.jsonl) and
human-readable summaries.
"""

from __future__ import annotations

import json
import os
import subprocess
import time
from dataclasses import dataclass, field, asdict
from enum import Enum
from pathlib import Path
from typing import Optional, Any


class JobState(str, Enum):
    """SLURM job states we track."""
    PENDING = "PENDING"
    RUNNING = "RUNNING"
    COMPLETED = "COMPLETED"
    FAILED = "FAILED"
    CANCELLED = "CANCELLED"
    TIMEOUT = "TIMEOUT"
    OUT_OF_MEMORY = "OUT_OF_MEMORY"
    NODE_FAIL = "NODE_FAIL"
    BOOT_FAIL = "BOOT_FAIL"
    UNKNOWN = "UNKNOWN"

    @property
    def is_terminal(self) -> bool:
        return self in (
            JobState.COMPLETED, JobState.FAILED, JobState.CANCELLED,
            JobState.TIMEOUT, JobState.OUT_OF_MEMORY, JobState.NODE_FAIL,
            JobState.BOOT_FAIL,
        )

    @property
    def is_success(self) -> bool:
        return self == JobState.COMPLETED


# Map SLURM state codes to our enum
SLURM_STATE_MAP = {
    "PD": JobState.PENDING,
    "R": JobState.RUNNING,
    "CG": JobState.RUNNING,  # completing
    "CD": JobState.COMPLETED,
    "F": JobState.FAILED,
    "CA": JobState.CANCELLED,
    "TO": JobState.TIMEOUT,
    "OOM": JobState.OUT_OF_MEMORY,
    "NF": JobState.NODE_FAIL,
    "BF": JobState.BOOT_FAIL,
}


@dataclass
class JobInfo:
    """Current info about a SLURM job."""
    job_id: str
    state: JobState
    name: str = ""
    partition: str = ""
    elapsed: str = ""
    time_limit: str = ""
    mem_used: str = ""
    cpu_used: str = ""
    exit_code: Optional[int] = None
    signal: Optional[int] = None
    node: str = ""


@dataclass
class MonitorEvent:
    """A single monitoring event."""
    timestamp: str
    job_id: str
    event_type: str  # job_submitted, job_pending, job_running, stage_started, etc.
    state: str
    message: str = ""
    data: dict[str, Any] = field(default_factory=dict)


class SlurmMonitor:
    """Monitor SLURM jobs using squeue and sacct.

    Usage:
        monitor = SlurmMonitor(job_id="12345", log_dir="/path/to/logs")
        info = monitor.get_job_info()
        monitor.wait_for_completion(timeout=7200)
    """

    def __init__(self, job_id: str, log_dir: str, job_name: str = ""):
        self.job_id = job_id
        self.log_dir = Path(log_dir)
        self.log_dir.mkdir(parents=True, exist_ok=True)
        self.job_name = job_name
        self.events: list[MonitorEvent] = []
        self._last_state: Optional[JobState] = None

    def _run_squeue(self) -> Optional[dict[str, str]]:
        """Query squeue for this job."""
        try:
            result = subprocess.run(
                ["squeue", "-h", "-j", self.job_id,
                 "-o", "%i|%j|%T|%P|%M|%l|%N|%R"],
                capture_output=True, text=True, timeout=10,
            )
            if result.returncode == 0 and result.stdout.strip():
                parts = result.stdout.strip().split("|")
                if len(parts) >= 8:
                    return {
                        "job_id": parts[0],
                        "name": parts[1],
                        "state_code": parts[2],
                        "partition": parts[3],
                        "elapsed": parts[4],
                        "time_limit": parts[5],
                        "node": parts[6],
                        "reason": parts[7],
                    }
        except Exception:
            pass
        return None

    def _parse_state(self, state_str: str) -> JobState:
        """Parse a SLURM state string (short code or full name) to JobState."""
        state_str = state_str.strip()
        # Try short code first
        if state_str in SLURM_STATE_MAP:
            return SLURM_STATE_MAP[state_str]
        # Try full name
        try:
            return JobState(state_str)
        except ValueError:
            pass
        # Try prefix match (e.g. "CANCELLED by user" -> CANCELLED)
        for state in JobState:
            if state_str.startswith(state.value):
                return state
        return JobState.UNKNOWN

    def _run_sacct(self) -> Optional[dict[str, str]]:
        """Query sacct for this job (terminal state info)."""
        try:
            result = subprocess.run(
                ["sacct", "-j", self.job_id, "-nP",
                 "--format=JobID,JobName,State,ExitCode,Elapsed,MaxRSS,Partition,NodeList",
                 "--noheader"],
                capture_output=True, text=True, timeout=10,
            )
            if result.returncode == 0 and result.stdout.strip():
                # Take the main job line (not .batch or .extern)
                for line in result.stdout.strip().split("\n"):
                    parts = line.split("|")
                    if len(parts) >= 8 and "." not in parts[0]:
                        return {
                            "job_id": parts[0],
                            "name": parts[1],
                            "state": parts[2],
                            "exit_code": parts[3],
                            "elapsed": parts[4],
                            "max_rss": parts[5],
                            "partition": parts[6],
                            "node": parts[7],
                        }
        except Exception:
            pass
        return None

    def get_job_info(self) -> JobInfo:
        """Get current job info from squeue or sacct."""
        # Try squeue first (for running/pending jobs)
        sq = self._run_squeue()
        if sq:
            state = self._parse_state(sq["state_code"])
            return JobInfo(
                job_id=self.job_id,
                state=state,
                name=sq["name"],
                partition=sq["partition"],
                elapsed=sq["elapsed"],
                time_limit=sq["time_limit"],
                node=sq["node"],
            )

        # Try sacct (for completed/failed jobs)
        sa = self._run_sacct()
        if sa:
            state_str = sa["state"].split()[0]  # e.g. "COMPLETED" or "CANCELLED by ..."
            state = JobState(state_str) if state_str in [s.value for s in JobState] else JobState.UNKNOWN
            exit_code = None
            signal = None
            if ":" in sa["exit_code"]:
                ec, sig = sa["exit_code"].split(":")
                exit_code = int(ec) if ec else None
                signal = int(sig) if sig else None
            return JobInfo(
                job_id=self.job_id,
                state=state,
                name=sa["name"],
                partition=sa["partition"],
                elapsed=sa["elapsed"],
                exit_code=exit_code,
                signal=signal,
                node=sa["node"],
            )

        return JobInfo(job_id=self.job_id, state=JobState.UNKNOWN)

    def record_event(self, event_type: str, state: JobState,
                     message: str = "", data: Optional[dict] = None):
        """Record a monitoring event and append to JSONL log."""
        event = MonitorEvent(
            timestamp=time.strftime("%Y-%m-%dT%H:%M:%S"),
            job_id=self.job_id,
            event_type=event_type,
            state=state.value,
            message=message,
            data=data or {},
        )
        self.events.append(event)

        # Append to JSONL
        path = self.log_dir / "interactive_monitor.jsonl"
        with open(path, "a") as f:
            f.write(json.dumps(asdict(event), default=str) + "\n")

    def check_output_files(self, expected_files: list[str]) -> dict[str, bool]:
        """Check if expected output files exist and are non-empty."""
        results = {}
        for f in expected_files:
            p = Path(f)
            results[f] = p.exists() and p.stat().st_size > 0
        return results

    def wait_for_completion(
        self,
        timeout: int = 7200,
        poll_interval: int = 10,
        expected_files: Optional[list[str]] = None,
        on_update=None,
    ) -> JobInfo:
        """Wait for job to reach terminal state.

        Args:
            timeout: max wait in seconds
            poll_interval: seconds between polls
            expected_files: files to check for existence after completion
            on_update: callback(job_info) called on each state change

        Returns:
            Final JobInfo
        """
        start = time.time()
        last_info = None

        while time.time() - start < timeout:
            info = self.get_job_info()

            # State change detection
            if info.state != self._last_state:
                if info.state == JobState.PENDING:
                    self.record_event("job_pending", info.state, f"Job pending: {info.node}")
                elif info.state == JobState.RUNNING:
                    self.record_event("job_running", info.state,
                                      f"Job running on {info.node}, elapsed={info.elapsed}")
                elif info.state.is_terminal:
                    msg = f"Job {info.state.value}, exit_code={info.exit_code}"
                    if info.state.is_success:
                        self.record_event("job_completed", info.state, msg)
                    else:
                        self.record_event("job_failed", info.state, msg)
                self._last_state = info.state

                if on_update:
                    on_update(info)

            last_info = info

            if info.state.is_terminal:
                # Check output files
                if expected_files:
                    file_status = self.check_output_files(expected_files)
                    self.record_event("output_check", info.state,
                                      f"Output files: {file_status}",
                                      data={"files": file_status})
                return info

            time.sleep(poll_interval)

        # Timeout
        self.record_event("monitor_timeout", JobState.UNKNOWN,
                          f"Monitor timed out after {timeout}s")
        return last_info or JobInfo(job_id=self.job_id, state=JobState.UNKNOWN)

    def write_report(self, output_path: str, final_info: JobInfo,
                     validation_results: Optional[dict] = None):
        """Write human-readable run report."""
        lines = [
            f"SLURM Job Monitor Report",
            f"{'=' * 60}",
            f"Job ID:      {self.job_id}",
            f"Job Name:    {self.job_name}",
            f"Final State: {final_info.state.value}",
            f"Exit Code:   {final_info.exit_code}",
            f"Partition:   {final_info.partition}",
            f"Node:        {final_info.node}",
            f"Elapsed:     {final_info.elapsed}",
            f"",
            f"Events ({len(self.events)} total):",
            f"{'-' * 60}",
        ]
        for e in self.events:
            lines.append(f"  [{e.timestamp}] {e.event_type}: {e.state} - {e.message}")

        if validation_results:
            lines.append(f"")
            lines.append(f"Validation Results:")
            lines.append(f"{'-' * 60}")
            for key, val in validation_results.items():
                status = "PASS" if val else "FAIL"
                lines.append(f"  {status}: {key}")

        lines.append(f"")
        lines.append(f"Report written: {time.strftime('%Y-%m-%dT%H:%M:%S')}")

        Path(output_path).write_text("\n".join(lines) + "\n")
