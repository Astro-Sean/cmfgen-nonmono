"""Orchestration layer for run_model_grid.py.

Provides typed configuration resolution, workflow state machine, manifests,
SLURM monitoring, and MPC interactive-partition test execution.

Modules:
    config           - ResolvedConfig dataclass, strict YAML schema validation
    state            - WorkflowState machine with persistent state transitions
    manifests        - Manifest and report writers (JSON/YAML/TXT)
    monitor          - SLURM job monitoring via squeue/sacct
    interactive_test - MPC interactive-partition test job generator and monitor
    site             - MPC HPC site configuration (partitions, accounts, QoS)
"""

__version__ = "0.1.0"

from .config import ResolvedConfig, ConfigError, resolve_config, validate_yaml
from .state import WorkflowState, StateMachine, WorkflowStates
from .manifests import ManifestWriter
from .monitor import SlurmMonitor, JobState
from .interactive_test import InteractiveTestRunner, InteractiveTestResult
from .site import SiteConfig, get_site_config
