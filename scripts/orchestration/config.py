"""Typed, immutable, resolved configuration from YAML.

This module provides:
- ResolvedConfig: a frozen dataclass capturing all workflow settings
- validate_yaml: strict schema validation with clear error messages
- resolve_config: load YAML, validate, resolve paths, compute derived values

The ResolvedConfig is the single source of truth — every executed or submitted
command must be generated from it. No defaults are scattered across functions.
"""

from __future__ import annotations

import hashlib
import json
import os
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any, Optional, Union

import yaml

from .site import SiteConfig, get_site_config, SiteConfigError


class ConfigError(Exception):
    """Raised when YAML configuration is invalid."""


# ── Epoch mapping ──────────────────────────────────────────────────────────

@dataclass(frozen=True)
class EpochMapping:
    """One entry in the canonical time-series mapping.

    requested_day -> stella_epoch -> output_label -> conversion_dir -> cmfgen_dir
    """
    requested_day: float
    stella_epoch: Optional[float]
    output_label: str
    conversion_dir: str
    cmfgen_dir: str
    status: str  # "exact", "remapped", "skipped", "unavailable"
    delta: float = 0.0  # |stella_epoch - requested_day|


# ── Resolved configuration ─────────────────────────────────────────────────

@dataclass(frozen=True)
class ResolvedConfig:
    """Immutable, fully-resolved configuration for one model run.

    All paths are absolute. All types are validated. No defaults are hidden
    in downstream functions.
    """

    # Identity
    model_name: str
    description: str
    run_id: str
    config_hash: str
    config_path: str

    # Project paths
    project_dir: str
    scripts_dir: str
    model_dir: str

    # Workflow flags
    hydro_backend: str  # "stella" or "heracles"
    use_stella: bool
    reuse_stella: bool
    restart_stella: bool
    skip_stella: bool
    run_cmfgen: bool
    run_sedona: bool
    run_sirocco: bool
    run_superlite: bool
    debug: bool
    debug_days: float
    cancel_jobs: bool
    force_cancel: bool

    # Time series
    snapshot_days: list[float]
    epoch_tolerance: Optional[float]
    epoch_mappings: list[EpochMapping]
    time_series: bool
    stella_target_days: float

    # STELLA settings (validated dict — code-specific params)
    stella: dict[str, Any]
    stella_partition: str
    stella_memory: str
    stella_time: str

    # CMFGEN conversion settings
    cmfgen_conversion: dict[str, Any]

    # CMFGEN run settings
    cmfgen: dict[str, Any]
    cmfgen_partition: str
    cmfgen_memory: str
    cmfgen_time: str

    # Spectrum codes
    sedona: dict[str, Any]
    superlite: dict[str, Any]
    sirocco: dict[str, Any]

    # Special modes
    kepler_presn: Optional[dict[str, Any]]
    sn_b_from_sn_a: Optional[dict[str, Any]]
    x_model_from_sn: Optional[dict[str, Any]]

    # Output paths
    stella_output_dir: str
    cmfgen_output_dir: str
    sedona_output_dir: str
    superlite_output_dir: str
    plots_dir: str
    output_dir: str

    # MPC site config
    site: SiteConfig

    # Raw YAML (for reference, not for execution decisions)
    raw_yaml: dict[str, Any]

    def to_dict(self) -> dict[str, Any]:
        """Serialize to a JSON-safe dict for manifest writing."""
        d = asdict(self)
        # EpochMapping list -> list of dicts
        d["epoch_mappings"] = [
            {"requested_day": e.requested_day, "stella_epoch": e.stella_epoch,
             "output_label": e.output_label, "conversion_dir": e.conversion_dir,
             "cmfgen_dir": e.cmfgen_dir, "status": e.status, "delta": e.delta}
            for e in self.epoch_mappings
        ]
        d["site"] = asdict(self.site)
        return d

    def to_yaml(self) -> str:
        """Serialize to YAML string for resolved_config.yaml."""
        return yaml.dump(self.to_dict(), default_flow_style=False, sort_keys=True)


# ── Validation ─────────────────────────────────────────────────────────────

# Required top-level fields
REQUIRED_FIELDS = {
    "model_name": str,
}

# Fields with known types for validation
FIELD_TYPES = {
    "model_name": str,
    "description": str,
    "run_cmfgen": bool,
    "run_sedona": bool,
    "run_sirocco": bool,
    "run_superlite": bool,
    "reuse_stella": bool,
    "restart_stella": bool,
    "skip_stella": bool,
    "debug": bool,
    "debug_days": (int, float),
    "cancel_jobs": bool,
    "force_cancel": bool,
    "hydro_backend": str,
}

# Mutually exclusive options
MUTUALLY_EXCLUSIVE = [
    ("reuse_stella", "skip_stella"),  # can't both be True (reuse implies skip)
]


def validate_yaml(raw: dict[str, Any]) -> list[str]:
    """Validate raw YAML dict. Returns list of error messages (empty = valid)."""
    errors = []

    # Required fields
    for fname, ftype in REQUIRED_FIELDS.items():
        if fname not in raw:
            errors.append(f"Missing required field: '{fname}'")
        elif not isinstance(raw[fname], ftype):
            errors.append(
                f"Field '{fname}' must be {ftype.__name__}, "
                f"got {type(raw[fname]).__name__}"
            )

    # Type checking
    for fname, ftype in FIELD_TYPES.items():
        if fname in raw:
            val = raw[fname]
            if not isinstance(val, ftype):
                errors.append(
                    f"Field '{fname}' must be {ftype}, "
                    f"got {type(val).__name__}: {val!r}"
                )

    # Mutually exclusive
    for f1, f2 in MUTUALLY_EXCLUSIVE:
        if raw.get(f1, False) and raw.get(f2, False):
            errors.append(f"Mutually exclusive: '{f1}' and '{f2}' are both True")

    # hydro_backend must be stella or heracles
    backend = raw.get("hydro_backend", "heracles")
    if backend.lower() not in ("stella", "heracles"):
        errors.append(
            f"hydro_backend must be 'stella' or 'heracles', got '{backend}'"
        )

    # snapshot_days validation
    cmfgen_conv = raw.get("cmfgen_conversion", {})
    if "snapshot_days" in cmfgen_conv:
        days = cmfgen_conv["snapshot_days"]
        if isinstance(days, (int, float)):
            pass  # single value is OK
        elif isinstance(days, list):
            if not days:
                errors.append("cmfgen_conversion.snapshot_days is empty list")
            for d in days:
                if not isinstance(d, (int, float)):
                    errors.append(
                        f"snapshot_days contains non-numeric value: {d!r}"
                    )
        else:
            errors.append(
                f"cmfgen_conversion.snapshot_days must be number or list, "
                f"got {type(days).__name__}"
            )

    # Duplicate label check (snapshot_days that would produce same label)
    # Only check if all values are numeric (non-numeric already flagged above)
    if "snapshot_days" in cmfgen_conv and isinstance(cmfgen_conv["snapshot_days"], list):
        days_list = cmfgen_conv["snapshot_days"]
        all_numeric = all(isinstance(d, (int, float)) for d in days_list)
        if all_numeric:
            labels = [f"{float(d):.2f}days" for d in days_list]
            seen = set()
            for lbl in labels:
                if lbl in seen:
                    errors.append(f"Duplicate snapshot day label: '{lbl}'")
                seen.add(lbl)

    # stella.target_days must be positive if set
    stella_cfg = raw.get("stella", {})
    if "target_days" in stella_cfg:
        td = stella_cfg["target_days"]
        if not isinstance(td, (int, float)) or td <= 0:
            errors.append(f"stella.target_days must be positive, got {td!r}")

    # min_nd / max_points relationships (min_nd is a floor, max_points is a ceiling)
    if "min_nd" in cmfgen_conv and "max_points" in cmfgen_conv:
        min_nd = cmfgen_conv["min_nd"]
        max_points = cmfgen_conv["max_points"]
        if isinstance(min_nd, int) and isinstance(max_points, int):
            if min_nd > max_points:
                errors.append(
                    f"cmfgen_conversion.min_nd ({min_nd}) > max_points ({max_points})"
                )

    # cmfgen.target_nd must be positive
    cmfgen_cfg = raw.get("cmfgen", {})
    if "target_nd" in cmfgen_cfg:
        tnd = cmfgen_cfg["target_nd"]
        if not isinstance(tnd, int) or tnd <= 0:
            errors.append(f"cmfgen.target_nd must be positive int, got {tnd!r}")

    # vel_law 11 (SN hydrodynamic, SET_RV_HYDRO_MODEL_V3) and 12 (SNIIN)
    # both support non-monotonic velocity via NON_MONO_SOLVER=T.
    # VEL_LAW=11 reads R/V/SIGMA from SN_HYDRO_DATA and preserves SIGMA_TRUE.
    # VEL_LAW=12 uses input_hydro.dat with tau-based regridding.
    if "vel_law" in cmfgen_cfg:
        vl = cmfgen_cfg["vel_law"]
        if vl not in (11, 12) and not cmfgen_cfg.get("depth_independent_solver", False):
            errors.append(
                f"cmfgen.vel_law={vl} is not supported. "
                f"Use 11 (SN hydrodynamic) or 12 (SNIIN) for non-monotonic velocity."
            )

    # KEPLER pre-SN file must exist if kepler_presn is configured
    kepler = raw.get("kepler_presn", {})
    if kepler and "presn_file" in kepler:
        p = Path(kepler["presn_file"])
        if not p.exists():
            errors.append(f"kepler_presn.presn_file not found: {p}")

    # SN_B from SN_A: sn_a_dir must exist
    sn_b = raw.get("sn_b_from_sn_a", {})
    if sn_b and "sn_a_dir" in sn_b:
        p = Path(sn_b["sn_a_dir"])
        if not p.exists():
            errors.append(f"sn_b_from_sn_a.sn_a_dir not found: {p}")

    # Partition names must be known (if explicitly set)
    for section, key in [("stella", "partition"), ("cmfgen", "partition"),
                         ("sedona", "partition"), ("superlite", "partition")]:
        sec = raw.get(section, {})
        if key in sec:
            from .site import KNOWN_PARTITIONS
            p = sec[key]
            if p not in KNOWN_PARTITIONS:
                errors.append(
                    f"{section}.{key}='{p}' is not a known partition. "
                    f"Known: {sorted(KNOWN_PARTITIONS.keys())}"
                )

    return errors


def _compute_config_hash(raw: dict[str, Any], config_path: str) -> str:
    """Compute SHA256 hash of the config content for identity tracking."""
    h = hashlib.sha256()
    h.update(json.dumps(raw, sort_keys=True, default=str).encode())
    h.update(str(config_path).encode())
    return h.hexdigest()[:16]


def _parse_snapshot_days(raw: dict[str, Any]) -> list[float]:
    """Extract and sort snapshot days from config."""
    cmfgen_conv = raw.get("cmfgen_conversion", {})
    days = cmfgen_conv.get("snapshot_days", 40)
    if isinstance(days, (int, float, str)):
        days = [float(days)]
    days = sorted([float(d) for d in days])
    if not days:
        days = [40.0]
    return days


def _get_stella_target_days(raw: dict[str, Any], snapshot_days: list[float]) -> float:
    """Compute STELLA target_days from config or snapshot days."""
    stella_cfg = raw.get("stella", {})
    cmfgen_conv = raw.get("cmfgen_conversion", {})

    # Explicit target_days wins
    if stella_cfg.get("target_days") is not None:
        return float(stella_cfg["target_days"])

    # Homologous evolution days are excluded from the max
    homo_days = set(cmfgen_conv.get("homologous_evolution_days", []))
    non_homo = [d for d in snapshot_days if d not in homo_days]
    extra = float(cmfgen_conv.get("stella_extra_days", 10))
    if non_homo:
        return max(non_homo) + extra
    return extra


def _build_epoch_mappings(
    snapshot_days: list[float],
    model_dir: str,
    epoch_tolerance: Optional[float],
) -> list[EpochMapping]:
    """Build epoch mapping entries. Without a .swd to read, all are 'pending'."""
    mappings = []
    for day in snapshot_days:
        label = f"{day:.2f}days"
        mappings.append(EpochMapping(
            requested_day=day,
            stella_epoch=None,  # resolved later when .swd is available
            output_label=label,
            conversion_dir=str(Path(model_dir) / "cmfgen" / label),
            cmfgen_dir=str(Path(model_dir) / "cmfgen" / label),
            status="pending",
            delta=0.0,
        ))
    return mappings


def resolve_config(
    config_path: str,
    project_dir: str,
    scripts_dir: str,
    run_id: Optional[str] = None,
    site_overrides: Optional[dict] = None,
) -> ResolvedConfig:
    """Load, validate, and resolve a YAML config into a ResolvedConfig.

    Raises ConfigError on any validation failure.
    """
    config_path = str(Path(config_path).resolve())
    project_dir = str(Path(project_dir).resolve())
    scripts_dir = str(Path(scripts_dir).resolve())

    with open(config_path, "r") as f:
        raw = yaml.safe_load(f)

    if raw is None:
        raise ConfigError(f"Empty config file: {config_path}")
    if not isinstance(raw, dict):
        raise ConfigError(f"Config must be a YAML mapping, got {type(raw).__name__}")

    # Validate
    errors = validate_yaml(raw)
    if errors:
        raise ConfigError(
            f"Config validation failed for {config_path}:\n  "
            + "\n  ".join(errors)
        )

    model_name = raw["model_name"]
    config_hash = _compute_config_hash(raw, config_path)
    run_id = run_id or f"{model_name}_{config_hash}"

    # Resolve paths
    output_cfg = raw.get("output", {})
    model_dir_str = output_cfg.get("model_dir", f"../models/{model_name}")
    model_dir_str = model_dir_str.replace("${model_name}", model_name)
    model_dir = str((Path(project_dir) / model_dir_str).resolve())

    # Snapshot days
    snapshot_days = _parse_snapshot_days(raw)
    epoch_tolerance = raw.get("cmfgen_conversion", {}).get("epoch_tolerance")
    epoch_mappings = _build_epoch_mappings(snapshot_days, model_dir, epoch_tolerance)
    time_series = len(snapshot_days) > 1

    # STELLA target days
    stella_target_days = _get_stella_target_days(raw, snapshot_days)

    # Output subdirectories
    stella_out = output_cfg.get("stella_output", "stella").replace("${model_name}", model_name)
    cmfgen_out = output_cfg.get("cmfgen_output", "cmfgen").replace("${model_name}", model_name)
    sedona_out = output_cfg.get("sedona_output", "sedona").replace("${model_name}", model_name)
    superlite_out = output_cfg.get("superlite_output", "superlite").replace("${model_name}", model_name)
    plots = output_cfg.get("plots", "plots").replace("${model_name}", model_name)
    outputs = output_cfg.get("outputs", "output").replace("${model_name}", model_name)

    # Site config
    yaml_site = raw.get("mpc", raw.get("site", {}))
    try:
        site = get_site_config(yaml_site=yaml_site, **(site_overrides or {}))
    except SiteConfigError as e:
        raise ConfigError(f"Site configuration error: {e}")

    # STELLA partition/memory/time
    stella_cfg = raw.get("stella", {})
    debug = raw.get("debug", False)
    stella_partition = "interactive" if debug else stella_cfg.get("partition", "general")
    stella_memory = "8G" if debug else stella_cfg.get("memory", "16G")
    stella_time = "01:30:00" if debug else stella_cfg.get("time", "23:59:00")

    # CMFGEN partition/memory/time
    cmfgen_cfg = raw.get("cmfgen", {})
    cmfgen_partition = "interactive" if debug else cmfgen_cfg.get("partition", "general")
    cmfgen_memory = cmfgen_cfg.get("mem", "16G")
    cmfgen_time = "01:30:00" if debug else cmfgen_cfg.get("time", "23:59:00")

    return ResolvedConfig(
        model_name=model_name,
        description=raw.get("description", ""),
        run_id=run_id,
        config_hash=config_hash,
        config_path=config_path,
        project_dir=project_dir,
        scripts_dir=scripts_dir,
        model_dir=model_dir,
        hydro_backend=raw.get("hydro_backend", "heracles").lower(),
        use_stella=raw.get("hydro_backend", "heracles").lower() == "stella",
        reuse_stella=raw.get("reuse_stella", False),
        restart_stella=raw.get("restart_stella", False),
        skip_stella=raw.get("skip_stella", False),
        run_cmfgen=raw.get("run_cmfgen", True),
        run_sedona=raw.get("run_sedona", True),
        run_sirocco=raw.get("run_sirocco", True),
        run_superlite=raw.get("run_superlite", False),
        debug=debug,
        debug_days=raw.get("debug_days", 1),
        cancel_jobs=raw.get("cancel_jobs", True),
        force_cancel=raw.get("force_cancel", True),
        snapshot_days=snapshot_days,
        epoch_tolerance=epoch_tolerance,
        epoch_mappings=epoch_mappings,
        time_series=time_series,
        stella_target_days=stella_target_days,
        stella=stella_cfg,
        stella_partition=stella_partition,
        stella_memory=stella_memory,
        stella_time=stella_time,
        cmfgen_conversion=raw.get("cmfgen_conversion", {}),
        cmfgen=cmfgen_cfg,
        cmfgen_partition=cmfgen_partition,
        cmfgen_memory=cmfgen_memory,
        cmfgen_time=cmfgen_time,
        sedona=raw.get("sedona", {}),
        superlite=raw.get("superlite", {}),
        sirocco=raw.get("sirocco", {}),
        kepler_presn=raw.get("kepler_presn"),
        sn_b_from_sn_a=raw.get("sn_b_from_sn_a"),
        x_model_from_sn=raw.get("x_model_from_sn"),
        stella_output_dir=str((Path(model_dir) / stella_out).resolve()),
        cmfgen_output_dir=str((Path(model_dir) / cmfgen_out).resolve()),
        sedona_output_dir=str((Path(model_dir) / sedona_out).resolve()),
        superlite_output_dir=str((Path(model_dir) / superlite_out).resolve()),
        plots_dir=str((Path(model_dir) / plots).resolve()),
        output_dir=str((Path(model_dir) / outputs).resolve()),
        site=site,
        raw_yaml=raw,
    )
