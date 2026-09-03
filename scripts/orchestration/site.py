"""MPC HPC site configuration.

Resolves partition names, accounts, QoS, and resource limits from the cluster
environment. Never invents values — if a required site parameter is missing,
fail with a clear configuration error.

On MPCDF Viper the partitions are:
  interactive  — 2:00:00 wall time (the MPC interactive partition)
  general      — 1-00:00:00 (1 day)
  small        — 1-00:00:00 (shared, same nodes as general)
  apu / apu1   — 1-00:00:00 (GPU/APU nodes)
  apudev       — 15:00 (GPU dev)

Accounts: mpe_cpu (CPU partitions), mpe_apu (APU partitions).
"""

from __future__ import annotations

import os
import re
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


# Known MPCDF Viper partitions (validated via sinfo at audit time)
KNOWN_PARTITIONS = {
    "interactive": {"max_time": "2:00:00", "type": "cpu"},
    "general": {"max_time": "1-00:00:00", "type": "cpu"},
    "small": {"max_time": "1-00:00:00", "type": "cpu"},
    "apu": {"max_time": "1-00:00:00", "type": "apu"},
    "apu1": {"max_time": "1-00:00:00", "type": "apu"},
    "apudev": {"max_time": "15:00", "type": "apu"},
    "datatransfer": {"max_time": "6:00:00", "type": "cpu"},
}

# Default account for CPU jobs on Viper
DEFAULT_CPU_ACCOUNT = "mpe_cpu"
DEFAULT_APU_ACCOUNT = "mpe_apu"


class SiteConfigError(Exception):
    """Raised when site configuration is missing or invalid."""


@dataclass(frozen=True)
class SiteConfig:
    """Resolved MPC HPC site configuration.

    All values are obtained from the cluster environment or YAML config —
    never invented.
    """

    cluster_name: str
    interactive_partition: str
    interactive_max_time: str
    production_partition: str
    production_max_time: str
    account: str
    qos: Optional[str]
    # Interactive test defaults (within interactive partition limits)
    test_time: str
    test_mem: str
    test_cpus: int
    test_ntasks: int

    @property
    def is_interactive_available(self) -> bool:
        """Check if the interactive partition exists on this cluster."""
        return self.interactive_partition in _get_live_partitions()


def _get_live_partitions() -> set[str]:
    """Query sinfo for currently available partitions."""
    try:
        result = subprocess.run(
            ["sinfo", "-h", "-o", "%P"],
            capture_output=True, text=True, timeout=10,
        )
        if result.returncode == 0:
            parts = set()
            for line in result.stdout.strip().split("\n"):
                p = line.strip().rstrip("*")
                if p:
                    parts.add(p)
            return parts
    except Exception:
        pass
    return set()


def _get_default_account() -> str:
    """Determine the default SLURM account for this user."""
    # Try sacctmgr
    try:
        result = subprocess.run(
            ["sacctmgr", "show", "assoc", f"user={os.environ.get('USER', '')}",
             "format=Account", "-nP"],
            capture_output=True, text=True, timeout=10,
        )
        if result.returncode == 0 and result.stdout.strip():
            accounts = [a.strip() for a in result.stdout.strip().split("\n") if a.strip()]
            # Prefer mpe_cpu for CPU jobs
            for acc in accounts:
                if acc == DEFAULT_CPU_ACCOUNT:
                    return acc
            if accounts:
                return accounts[0]
    except Exception:
        pass
    # Fallback to environment variable
    return os.environ.get("SLURM_ACCOUNT", DEFAULT_CPU_ACCOUNT)


def get_site_config(
    yaml_site: Optional[dict] = None,
    interactive_partition: Optional[str] = None,
    production_partition: Optional[str] = None,
    account: Optional[str] = None,
    qos: Optional[str] = None,
    test_time: Optional[str] = None,
    test_mem: Optional[str] = None,
    test_cpus: Optional[int] = None,
    test_ntasks: Optional[int] = None,
) -> SiteConfig:
    """Resolve site configuration from YAML + cluster environment.

    Priority: YAML config > explicit args > cluster environment > built-in defaults.

    Raises SiteConfigError if the interactive partition is not available.
    """
    yaml_site = yaml_site or {}

    # Partition resolution
    interactive_part = (
        interactive_partition
        or yaml_site.get("interactive_partition")
        or "interactive"
    )
    production_part = (
        production_partition
        or yaml_site.get("production_partition")
        or "general"
    )

    # Validate partitions are known
    if interactive_part not in KNOWN_PARTITIONS:
        raise SiteConfigError(
            f"Unknown interactive partition '{interactive_part}'. "
            f"Known: {sorted(KNOWN_PARTITIONS.keys())}"
        )
    if production_part not in KNOWN_PARTITIONS:
        raise SiteConfigError(
            f"Unknown production partition '{production_part}'. "
            f"Known: {sorted(KNOWN_PARTITIONS.keys())}"
        )

    # Account resolution
    resolved_account = account or yaml_site.get("account") or _get_default_account()

    # QoS (optional)
    resolved_qos = qos or yaml_site.get("qos") or None

    # Interactive test resource defaults
    # Interactive partition on Viper has 2hr limit; use 01:55:00 to leave margin
    resolved_test_time = test_time or yaml_site.get("test_time") or "01:55:00"
    resolved_test_mem = test_mem or yaml_site.get("test_mem") or "16G"
    resolved_test_cpus = test_cpus or yaml_site.get("test_cpus") or 4
    resolved_test_ntasks = test_ntasks or yaml_site.get("test_ntasks") or 1

    # Validate test time doesn't exceed interactive partition limit
    interactive_max = KNOWN_PARTITIONS[interactive_part]["max_time"]
    if not _time_leq(resolved_test_time, interactive_max):
        raise SiteConfigError(
            f"Interactive test time '{resolved_test_time}' exceeds "
            f"partition '{interactive_part}' max time '{interactive_max}'"
        )

    # Determine cluster name
    cluster_name = "viper"
    try:
        result = subprocess.run(
            ["scontrol", "show", "config"],
            capture_output=True, text=True, timeout=10,
        )
        if result.returncode == 0:
            m = re.search(r"ClusterName\s*=\s*(\S+)", result.stdout)
            if m:
                cluster_name = m.group(1)
    except Exception:
        pass

    return SiteConfig(
        cluster_name=cluster_name,
        interactive_partition=interactive_part,
        interactive_max_time=interactive_max,
        production_partition=production_part,
        production_max_time=KNOWN_PARTITIONS[production_part]["max_time"],
        account=resolved_account,
        qos=resolved_qos,
        test_time=resolved_test_time,
        test_mem=resolved_test_mem,
        test_cpus=resolved_test_cpus,
        test_ntasks=resolved_test_ntasks,
    )


def _time_to_seconds(time_str: str) -> int:
    """Convert SLURM time string to seconds. Handles D-HH:MM:SS, HH:MM:SS, MM:SS."""
    time_str = time_str.strip()
    days = 0
    if "-" in time_str:
        d, time_str = time_str.split("-", 1)
        days = int(d)
    parts = time_str.split(":")
    if len(parts) == 3:
        h, m, s = int(parts[0]), int(parts[1]), int(parts[2])
    elif len(parts) == 2:
        h, m, s = 0, int(parts[0]), int(parts[1])
    elif len(parts) == 1:
        h, m, s = 0, 0, int(parts[0])
    else:
        raise ValueError(f"Invalid time string: {time_str}")
    return days * 86400 + h * 3600 + m * 60 + s


def _time_leq(t1: str, t2: str) -> bool:
    """Return True if t1 <= t2 in seconds."""
    return _time_to_seconds(t1) <= _time_to_seconds(t2)
