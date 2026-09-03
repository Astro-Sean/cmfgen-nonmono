"""Workflow state machine with persistent state transitions.

Persist state after every transition. A job ID alone never proves completion —
require SLURM terminal-state and content-validation evidence.

States:
    CONFIG_VALIDATED
    MPC_INTERACTIVE_TEST_SUBMITTED
    MPC_INTERACTIVE_TEST_RUNNING
    MPC_INTERACTIVE_TEST_PASSED
    MPC_INTERACTIVE_TEST_FAILED
    PROFILE_READY
    STELLA_SUBMITTED
    STELLA_COMPLETE
    STELLA_REUSED
    SNAPSHOTS_MAPPED
    SNAPSHOTS_CONVERTED
    GRID_VALIDATED
    POPULATIONS_READY
    CMFGEN_SEEDED
    CMFGEN_SUBMITTED
    CMFGEN_COMPLETE
    OPTIONAL_SPECTRA_SUBMITTED
    COMPLETE
    FAILED
"""

from __future__ import annotations

import json
import time
from dataclasses import dataclass, field, asdict
from enum import Enum
from pathlib import Path
from typing import Optional, Any


class WorkflowStates(str, Enum):
    """All possible workflow states."""
    INIT = "INIT"
    CONFIG_VALIDATED = "CONFIG_VALIDATED"
    MPC_INTERACTIVE_TEST_SUBMITTED = "MPC_INTERACTIVE_TEST_SUBMITTED"
    MPC_INTERACTIVE_TEST_RUNNING = "MPC_INTERACTIVE_TEST_RUNNING"
    MPC_INTERACTIVE_TEST_PASSED = "MPC_INTERACTIVE_TEST_PASSED"
    MPC_INTERACTIVE_TEST_FAILED = "MPC_INTERACTIVE_TEST_FAILED"
    PROFILE_READY = "PROFILE_READY"
    STELLA_SUBMITTED = "STELLA_SUBMITTED"
    STELLA_COMPLETE = "STELLA_COMPLETE"
    STELLA_REUSED = "STELLA_REUSED"
    SNAPSHOTS_MAPPED = "SNAPSHOTS_MAPPED"
    SNAPSHOTS_CONVERTED = "SNAPSHOTS_CONVERTED"
    GRID_VALIDATED = "GRID_VALIDATED"
    POPULATIONS_READY = "POPULATIONS_READY"
    CMFGEN_SEEDED = "CMFGEN_SEEDED"
    CMFGEN_SUBMITTED = "CMFGEN_SUBMITTED"
    CMFGEN_COMPLETE = "CMFGEN_COMPLETE"
    OPTIONAL_SPECTRA_SUBMITTED = "OPTIONAL_SPECTRA_SUBMITTED"
    COMPLETE = "COMPLETE"
    FAILED = "FAILED"


# Valid state transitions (from -> set of allowed next states)
VALID_TRANSITIONS: dict[WorkflowStates, set[WorkflowStates]] = {
    WorkflowStates.INIT: {
        WorkflowStates.CONFIG_VALIDATED,
        WorkflowStates.FAILED,
    },
    WorkflowStates.CONFIG_VALIDATED: {
        WorkflowStates.MPC_INTERACTIVE_TEST_SUBMITTED,
        WorkflowStates.PROFILE_READY,  # skip interactive test in --no-submit
        WorkflowStates.FAILED,
    },
    WorkflowStates.MPC_INTERACTIVE_TEST_SUBMITTED: {
        WorkflowStates.MPC_INTERACTIVE_TEST_RUNNING,
        WorkflowStates.MPC_INTERACTIVE_TEST_FAILED,
        WorkflowStates.FAILED,
    },
    WorkflowStates.MPC_INTERACTIVE_TEST_RUNNING: {
        WorkflowStates.MPC_INTERACTIVE_TEST_PASSED,
        WorkflowStates.MPC_INTERACTIVE_TEST_FAILED,
        WorkflowStates.FAILED,
    },
    WorkflowStates.MPC_INTERACTIVE_TEST_PASSED: {
        WorkflowStates.PROFILE_READY,
        WorkflowStates.FAILED,
    },
    WorkflowStates.MPC_INTERACTIVE_TEST_FAILED: {
        WorkflowStates.FAILED,
    },
    WorkflowStates.PROFILE_READY: {
        WorkflowStates.STELLA_SUBMITTED,
        WorkflowStates.STELLA_REUSED,
        WorkflowStates.SNAPSHOTS_MAPPED,  # skip_stella
        WorkflowStates.FAILED,
    },
    WorkflowStates.STELLA_SUBMITTED: {
        WorkflowStates.STELLA_COMPLETE,
        WorkflowStates.FAILED,
    },
    WorkflowStates.STELLA_COMPLETE: {
        WorkflowStates.SNAPSHOTS_MAPPED,
        WorkflowStates.FAILED,
    },
    WorkflowStates.STELLA_REUSED: {
        WorkflowStates.SNAPSHOTS_MAPPED,
        WorkflowStates.FAILED,
    },
    WorkflowStates.SNAPSHOTS_MAPPED: {
        WorkflowStates.SNAPSHOTS_CONVERTED,
        WorkflowStates.FAILED,
    },
    WorkflowStates.SNAPSHOTS_CONVERTED: {
        WorkflowStates.GRID_VALIDATED,
        WorkflowStates.FAILED,
    },
    WorkflowStates.GRID_VALIDATED: {
        WorkflowStates.POPULATIONS_READY,
        WorkflowStates.CMFGEN_SUBMITTED,  # skip populations if not using proxy
        WorkflowStates.FAILED,
    },
    WorkflowStates.POPULATIONS_READY: {
        WorkflowStates.CMFGEN_SEEDED,
        WorkflowStates.CMFGEN_SUBMITTED,
        WorkflowStates.FAILED,
    },
    WorkflowStates.CMFGEN_SEEDED: {
        WorkflowStates.CMFGEN_SUBMITTED,
        WorkflowStates.FAILED,
    },
    WorkflowStates.CMFGEN_SUBMITTED: {
        WorkflowStates.CMFGEN_COMPLETE,
        WorkflowStates.OPTIONAL_SPECTRA_SUBMITTED,
        WorkflowStates.FAILED,
    },
    WorkflowStates.CMFGEN_COMPLETE: {
        WorkflowStates.OPTIONAL_SPECTRA_SUBMITTED,
        WorkflowStates.COMPLETE,
        WorkflowStates.FAILED,
    },
    WorkflowStates.OPTIONAL_SPECTRA_SUBMITTED: {
        WorkflowStates.COMPLETE,
        WorkflowStates.FAILED,
    },
    WorkflowStates.COMPLETE: set(),
    WorkflowStates.FAILED: set(),
}


@dataclass
class StateTransition:
    """Record of a single state transition."""
    from_state: str
    to_state: str
    timestamp: str
    job_id: Optional[str] = None
    evidence: Optional[dict[str, Any]] = None
    message: str = ""


@dataclass
class WorkflowState:
    """Persistent workflow state for one model run.

    Stored as workflow_state.json in the model output directory.
    """
    run_id: str
    model_name: str
    current_state: str
    config_hash: str
    transitions: list[StateTransition] = field(default_factory=list)
    job_ids: dict[str, str] = field(default_factory=dict)
    epoch_status: dict[str, str] = field(default_factory=dict)
    created_at: str = ""
    updated_at: str = ""

    def to_dict(self) -> dict[str, Any]:
        return {
            "run_id": self.run_id,
            "model_name": self.model_name,
            "current_state": self.current_state,
            "config_hash": self.config_hash,
            "transitions": [asdict(t) for t in self.transitions],
            "job_ids": self.job_ids,
            "epoch_status": self.epoch_status,
            "created_at": self.created_at,
            "updated_at": self.updated_at,
        }


class StateMachine:
    """Manages workflow state transitions with persistence.

    Usage:
        sm = StateMachine(run_id, model_name, config_hash, state_file)
        sm.transition(WorkflowStates.CONFIG_VALIDATED, message="Config validated")
        sm.transition(WorkflowStates.PROFILE_READY, message="Profile generated")
    """

    def __init__(
        self,
        run_id: str,
        model_name: str,
        config_hash: str,
        state_file: str,
    ):
        self.state_file = Path(state_file)
        self.state = self._load_or_init(run_id, model_name, config_hash)

    def _load_or_init(self, run_id, model_name, config_hash) -> WorkflowState:
        """Load existing state or initialize new."""
        if self.state_file.exists():
            try:
                with open(self.state_file, "r") as f:
                    data = json.load(f)
                state = WorkflowState(
                    run_id=data["run_id"],
                    model_name=data["model_name"],
                    current_state=data["current_state"],
                    config_hash=data["config_hash"],
                    job_ids=data.get("job_ids", {}),
                    epoch_status=data.get("epoch_status", {}),
                    created_at=data.get("created_at", ""),
                    updated_at=data.get("updated_at", ""),
                )
                state.transitions = [
                    StateTransition(**t) for t in data.get("transitions", [])
                ]
                # Verify config hash matches (reject stale state from different config)
                if state.config_hash != config_hash:
                    # Config changed — start fresh
                    return self._init_fresh(run_id, model_name, config_hash)
                return state
            except (json.JSONDecodeError, KeyError, TypeError):
                pass  # Corrupt state file — start fresh
        return self._init_fresh(run_id, model_name, config_hash)

    def _init_fresh(self, run_id, model_name, config_hash) -> WorkflowState:
        now = time.strftime("%Y-%m-%dT%H:%M:%S")
        return WorkflowState(
            run_id=run_id,
            model_name=model_name,
            current_state="INIT",
            config_hash=config_hash,
            created_at=now,
            updated_at=now,
        )

    def transition(
        self,
        to_state: WorkflowStates,
        job_id: Optional[str] = None,
        evidence: Optional[dict] = None,
        message: str = "",
        force: bool = False,
    ) -> bool:
        """Execute a state transition. Returns True if successful.

        Raises ValueError on invalid transition unless force=True.
        """
        from_state = WorkflowStates(self.state.current_state)

        if to_state == from_state:
            return True  # Idempotent

        if not force:
            allowed = VALID_TRANSITIONS.get(from_state, set())
            if to_state not in allowed:
                raise ValueError(
                    f"Invalid state transition: {from_state.value} -> {to_state.value}. "
                    f"Allowed: {[s.value for s in allowed] or '(terminal)'}"
                )

        now = time.strftime("%Y-%m-%dT%H:%M:%S")
        transition = StateTransition(
            from_state=from_state.value,
            to_state=to_state.value,
            timestamp=now,
            job_id=job_id,
            evidence=evidence,
            message=message,
        )
        self.state.transitions.append(transition)
        self.state.current_state = to_state.value
        self.state.updated_at = now
        if job_id:
            self.state.job_ids[to_state.value] = job_id

        self._save()
        return True

    def record_job_id(self, stage: str, job_id: str):
        """Record a job ID for a stage."""
        self.state.job_ids[stage] = job_id
        self.state.updated_at = time.strftime("%Y-%m-%dT%H:%M:%S")
        self._save()

    def record_epoch_status(self, day: float, status: str):
        """Record status for a specific epoch."""
        self.state.epoch_status[f"{day:.2f}days"] = status
        self.state.updated_at = time.strftime("%Y-%m-%dT%H:%M:%S")
        self._save()

    def is_terminal(self) -> bool:
        """Check if the workflow is in a terminal state."""
        return self.state.current_state in (
            WorkflowStates.COMPLETE.value,
            WorkflowStates.FAILED.value,
        )

    def is_complete(self) -> bool:
        return self.state.current_state == WorkflowStates.COMPLETE.value

    def is_failed(self) -> bool:
        return self.state.current_state == WorkflowStates.FAILED.value

    def _save(self):
        """Persist state to JSON file."""
        self.state_file.parent.mkdir(parents=True, exist_ok=True)
        with open(self.state_file, "w") as f:
            json.dump(self.state.to_dict(), f, indent=2)

    def to_dict(self) -> dict:
        return self.state.to_dict()
