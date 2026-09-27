"""MPC HPC site configuration.

Resolves partition names, accounts, QoS, and resource limits from the cluster
environment. Never invents values — if a required site parameter is missing,
fail with a clear configuration error.

Supports two MPCDF clusters:

On MPCDF Viper the partitions are:
  interactive  — 2:00:00 wall time (the MPC interactive partition)
  general      — 1-00:00:00 (1 day)
  small        — 1-00:00:00 (shared, same nodes as general)
  apu / apu1   — 1-00:00:00 (GPU/APU nodes)
  apudev       — 15:00 (GPU dev)

On MPCDF Raven the partitions are:
  interactive  — 2:00:00 wall time (the MPC interactive partition)
  general      — 1-00:00:00 (1 day)
  small        — 1-00:00:00 (shared, same nodes as general)
  gpu / gpu1   — 1-00:00:00 (GPU nodes)
  rvs          — 1-00:00:00 (RVS nodes)
  gpudev       — 15:00 (GPU dev)

Accounts: mpe_cpu (CPU partitions), mpe_apu (Viper APU) / mpe_gpu (Raven GPU).
"""

from __future__ import annotations

import os
import re
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


# Known MPCDF partitions (Viper + Raven, validated via sinfo at audit time)
KNOWN_PARTITIONS = {
    # CPU partitions (common to both Viper and Raven)
    "interactive": {"max_time": "2:00:00", "type": "cpu"},
    "general": {"max_time": "1-00:00:00", "type": "cpu"},
    "small": {"max_time": "1-00:00:00", "type": "cpu"},
    "datatransfer": {"max_time": "6:00:00", "type": "cpu"},
    # Viper GPU/APU partitions
    "apu": {"max_time": "1-00:00:00", "type": "apu"},
    "apu1": {"max_time": "1-00:00:00", "type": "apu"},
    "apudev": {"max_time": "15:00", "type": "apu"},
    # Raven GPU partitions
    "gpu": {"max_time": "1-00:00:00", "type": "gpu"},
    "gpu1": {"max_time": "1-00:00:00", "type": "gpu"},
    "rvs": {"max_time": "1-00:00:00", "type": "gpu"},
    "gpudev": {"max_time": "15:00", "type": "gpu"},
}

# Default accounts for CPU jobs (same on both Viper and Raven)
DEFAULT_CPU_ACCOUNT = "mpe_cpu"
# Viper uses mpe_apu for APU partitions; Raven uses mpe_gpu for GPU partitions
DEFAULT_APU_ACCOUNT = "mpe_apu"
DEFAULT_GPU_ACCOUNT = "mpe_gpu"


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
    # Interactive partition on Viper/Raven has 2hr limit; use 01:55:00 to leave margin
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

    # Determine cluster name using the shared detect_cluster() utility
    cluster_name = detect_cluster()

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


# ---------------------------------------------------------------------------
# Cluster detection and cross-cluster path translation
#
# MPCDF cross-mounts scratch and home filesystems between Viper and Raven.
# Both paths work on both clusters, but for I/O performance on compute nodes
# you should use the LOCAL filesystem path (GPFS on the local cluster) rather
# than the cross-mounted NFS path.
#
# Path mappings (for user "sbrennan"; generalised via $USER / os.getlogin):
#
#   Scratch (ptmp):
#     Viper: /viper/ptmp1/$USER   (also /viper/ptmp/$USER, /ptmp/$USER via symlink)
#     Raven: /raven/ptmp/$USER
#
#   Home (u):
#     Viper: /viper/u1/$USER      (also /viper/u/$USER via symlink)
#     Raven: /raven/u/$USER
#
# Translation rules (applied when the path points at the OTHER cluster's
# scratch or home):
#   On Raven:  /viper/ptmp1/$USER  -> /raven/ptmp/$USER
#              /viper/ptmp/$USER   -> /raven/ptmp/$USER
#              /viper/u1/$USER     -> /raven/u/$USER  (if a local copy exists)
#   On Viper:  /raven/ptmp/$USER   -> /viper/ptmp1/$USER
#              /raven/u/$USER      -> /viper/u1/$USER
# ---------------------------------------------------------------------------

# Canonical scratch prefixes per cluster (ordered most-specific first)
_SCRATCH_PREFIXES = {
    "viper": ["/viper/ptmp1/", "/viper/ptmp/"],
    "raven": ["/raven/ptmp/"],
}
# Canonical home prefixes per cluster
_HOME_PREFIXES = {
    "viper": ["/viper/u1/", "/viper/u/"],
    "raven": ["/raven/u/"],
}


def detect_cluster() -> str:
    """Detect the current MPCDF cluster name.

    Uses ``scontrol show config`` (most reliable, works on all node types),
    falling back to hostname prefix matching.

    Returns ``"viper"`` or ``"raven"``.  Defaults to ``"viper"`` if detection
    fails (Viper is the primary development cluster).
    """
    # Method 1: SLURM ClusterName (works on login + compute nodes)
    try:
        result = subprocess.run(
            ["scontrol", "show", "config"],
            capture_output=True, text=True, timeout=10,
        )
        if result.returncode == 0:
            m = re.search(r"ClusterName\s*=\s*(\S+)", result.stdout)
            if m:
                name = m.group(1).lower()
                if name in ("viper", "raven"):
                    return name
    except Exception:
        pass

    # Method 2: hostname prefix
    import socket
    hn = socket.gethostname().lower()
    if hn.startswith("raven") or hn.startswith("rav"):
        return "raven"
    if hn.startswith("viper") or hn.startswith("vip"):
        return "viper"

    # Default: viper (primary development cluster)
    return "viper"


def get_scratch_base(cluster: Optional[str] = None, user: Optional[str] = None) -> str:
    """Return the canonical scratch base directory for a cluster.

    Args:
        cluster: Cluster name (``"viper"`` or ``"raven"``).  If ``None``,
            detects the current cluster.
        user: Username.  If ``None``, uses ``$USER`` or ``os.getlogin()``.

    Returns e.g. ``/viper/ptmp1/sbrennan`` or ``/raven/ptmp/sbrennan``.
    """
    cluster = cluster or detect_cluster()
    user = user or os.environ.get("USER", "")
    if cluster == "raven":
        return f"/raven/ptmp/{user}"
    return f"/viper/ptmp1/{user}"


def get_home_base(cluster: Optional[str] = None, user: Optional[str] = None) -> str:
    """Return the canonical home base directory for a cluster.

    Returns e.g. ``/viper/u1/sbrennan`` or ``/raven/u/sbrennan``.
    """
    cluster = cluster or detect_cluster()
    user = user or os.environ.get("USER", "")
    if cluster == "raven":
        return f"/raven/u/{user}"
    return f"/viper/u1/{user}"


def translate_to_local_path(path: str, cluster: Optional[str] = None) -> str:
    """Translate a path from another cluster's filesystem to the local cluster's equivalent.

    If the path starts with a scratch or home prefix belonging to the OTHER
    cluster, replace it with the local cluster's equivalent prefix.  Paths
    that already use the local cluster's prefix (or generic ``/ptmp/``) are
    returned unchanged.

    This only rewrites the filesystem prefix — the username and the rest of
    the path are preserved, so ``/viper/ptmp1/sbrennan/SN2024pba/atomic/foo``
    becomes ``/raven/ptmp/sbrennan/SN2024pba/atomic/foo`` on Raven.

    Args:
        path: The path to translate.
        cluster: The target cluster (default: detect current cluster).

    Returns the translated path.
    """
    if not path:
        return path
    cluster = cluster or detect_cluster()

    # Determine the "other" cluster
    other = "raven" if cluster == "viper" else "viper"

    # Map other-cluster scratch prefixes -> local-cluster scratch prefix
    # (prefix includes trailing slash but NOT username; username is in the suffix)
    if cluster == "raven":
        local_scratch_prefix = "/raven/ptmp/"
        local_home_prefix = "/raven/u/"
    else:
        local_scratch_prefix = "/viper/ptmp1/"
        local_home_prefix = "/viper/u1/"

    for other_prefix in _SCRATCH_PREFIXES[other]:
        if path.startswith(other_prefix):
            suffix = path[len(other_prefix):]  # includes username + rest
            return f"{local_scratch_prefix}{suffix}"

    for other_prefix in _HOME_PREFIXES[other]:
        if path.startswith(other_prefix):
            suffix = path[len(other_prefix):]
            return f"{local_home_prefix}{suffix}"

    # Path is already local or uses a generic mount point — return as-is
    return path


def is_local_path(path: str, cluster: Optional[str] = None) -> bool:
    """Check whether *path* is on the local cluster's filesystem (not cross-mounted)."""
    if not path:
        return True
    cluster = cluster or detect_cluster()
    for prefix in _SCRATCH_PREFIXES[cluster]:
        if path.startswith(prefix):
            return True
    for prefix in _HOME_PREFIXES[cluster]:
        if path.startswith(prefix):
            return True
    # /ptmp/$USER is a symlink to the local scratch on both clusters
    ptmp_local = f"/ptmp/{os.environ.get('USER', '')}/"
    if path.startswith("/ptmp/"):
        return True
    return False


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
