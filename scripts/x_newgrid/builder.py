"""One-pass X_newgrid grid builder.

Orchestrates the full hierarchical grid construction:

    STELLA source
        |
    trim (tau / velocity limits)
        |
    detect features (shocks, interfaces, extrema, dv/dr, ion fronts)
        |
    build log-tau skeleton  (~N_decades * points_per_decade)   [Gap B]
        |
    merge mandatory anchor radii
        |
    hierarchical-priority refinement                          [Gap A]
        |
    coarsen redundant points (anchors protected)              [Gap E]
        |
    local interpolation-error check                           [Gap D]
        |
    write SN_HYDRO_DATA (via grid_mapping.cmfgen_writer)

The two-pass CMFGEN ionisation-front refinement (Pass 2) is a separate
standalone script; this builder produces the Pass-1 grid.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np

from grid_mapping.shell_model import (SourceShellGrid, TargetGrid,
                                      compute_shell_volume,
                                      compute_shell_centres,
                                      density_from_mass_volume)
from grid_mapping.feature_detection import detect_all_features
from grid_mapping.conservative_remap import conservative_remap
from grid_mapping.opacity_tau import compute_tau
from grid_mapping.cmfgen_writer import write_sn_hydro_data
from grid_mapping.ionization_metrics import saha_electron_density

from .config import XNewGridConfig
from .skeleton import build_tau_skeleton
from .hierarchy import hierarchical_refine
from .coarsen import coarsen_grid
from .error_metric import local_interpolation_error, ErrorReport


def _trim_source(source: SourceShellGrid, cfg: XNewGridConfig
                 ) -> SourceShellGrid:
    """Trim the source by tau and velocity limits, re-detecting features."""
    s = source
    n = s.n_shells
    if n == 0:
        return s
    keep = np.ones(n, dtype=bool)
    if len(s.tau_cumulative) == n:
        keep &= s.tau_cumulative >= cfg.tau_min
        keep &= s.tau_cumulative <= cfg.tau_max
    if cfg.v_min_cms > 0:
        keep &= np.abs(s.velocity_cms) >= cfg.v_min_cms
    if cfg.v_max_cms > 0:
        keep &= np.abs(s.velocity_cms) <= cfg.v_max_cms
    if keep.sum() < 2:
        keep[:2] = True

    indices = np.where(keep)[0]
    r_bound_new = np.empty(len(indices) + 1)
    r_bound_new[1:-1] = np.sqrt(s.r_bound_cm[:-1][indices[:-1]] *
                                s.r_bound_cm[1:][indices[:-1]])
    r_bound_new[0] = s.r_bound_cm[0]
    r_bound_new[-1] = s.r_bound_cm[-1]
    for i in range(len(r_bound_new) - 1):
        if r_bound_new[i] <= r_bound_new[i + 1]:
            r_bound_new[i] = r_bound_new[i + 1] * 1.0001

    vol = compute_shell_volume(r_bound_new)
    r_centre = compute_shell_centres(r_bound_new)
    mass = s.shell_mass_g[indices]
    rho = density_from_mass_volume(mass, vol)
    mf = {e: f[indices] for e, f in s.mass_fractions.items()}
    iso = {k: f[indices] for k, f in s.isotope_fractions.items()}
    ni56 = s.ni56_initial[indices] if (
        s.ni56_initial is not None and hasattr(s.ni56_initial, '__len__')) else None
    kappa = s.kappa_cgs[indices] if len(s.kappa_cgs) == n else np.full(len(indices), 0.34)
    tau = compute_tau(r_bound_new, rho, kappa)

    trimmed = SourceShellGrid(
        r_bound_cm=r_bound_new, r_centre_cm=r_centre,
        shell_volume_cm3=vol, shell_mass_g=mass, density_gcm3=rho,
        velocity_cms=s.velocity_cms[indices],
        temperature_K=s.temperature_K[indices],
        kappa_cgs=kappa, mass_fractions=mf, isotope_fractions=iso,
        ni56_initial=ni56, tau_cumulative=tau, time_days=s.time_days,
        source_file=s.source_file, n_zones_original=s.n_zones_original,
    )
    # re-detect features on the trimmed grid (sets masks in-place)
    detect_all_features(trimmed)
    return trimmed


def _centres_to_boundaries(r_centre: np.ndarray, cfg: XNewGridConfig
                           ) -> np.ndarray:
    """Convert centre radii (decreasing) to shell boundaries (decreasing)."""
    r = np.asarray(r_centre, dtype=float)
    n = len(r)
    if n < 2:
        return r.copy()
    r_bound = np.empty(n + 1)
    r_bound[1:-1] = np.sqrt(r[:-1] * r[1:])
    r_bound[0] = r[0] ** 2 / r_bound[1] if r_bound[1] > 0 else r[0] * 1.01
    r_bound[-1] = r[-1] ** 2 / r_bound[-2] if r_bound[-2] > 0 else r[-1] * 0.99
    for i in range(len(r_bound) - 1):
        if r_bound[i] <= r_bound[i + 1]:
            r_bound[i] = r_bound[i + 1] * (1.0 + cfg.min_dr_r)
    return r_bound


def _enforce_min_spacing(r_bound: np.ndarray, min_dr: float) -> np.ndarray:
    """Remove boundaries that are closer than ``min_dr`` to their neighbour.

    Keeps the first and last boundaries.  This prevents duplicate / near-
    duplicate cells that produce zero-mass zones in the conservative remap.
    """
    if len(r_bound) < 3:
        return r_bound
    keep = [0]
    for i in range(1, len(r_bound) - 1):
        if abs(r_bound[keep[-1]] - r_bound[i]) >= min_dr:
            keep.append(i)
    keep.append(len(r_bound) - 1)
    # final check: last pair
    if abs(r_bound[keep[-2]] - r_bound[keep[-1]]) < min_dr:
        keep.pop(-2)
    return r_bound[keep]


def _merge_duplicate_cells(target: TargetGrid, cfg: XNewGridConfig
                            ) -> TargetGrid:
    """Merge cells whose centre radii are not strictly decreasing.

    Returns the target unchanged if the radii are already strictly decreasing.
    Otherwise, drops the duplicate cells and interpolates their quantities
    onto the surviving grid.  This is a simple cell-merge, not a full
    conservative re-remap — the duplicate cells are zero-mass artifacts
    from the conservative remap's neighbour fixup, so dropping them does
    not lose physical mass.
    """
    r = target.r_cm
    if len(r) < 2:
        return target
    dr = np.diff(r)
    if np.all(dr < 0):
        return target
    # find strictly-decreasing subset
    keep = [0]
    for i in range(1, len(r)):
        if r[i] < r[keep[-1]] - cfg.min_dr_r * abs(r[keep[-1]]):
            keep.append(i)
    if len(keep) == len(r):
        return target
    keep_arr = np.array(keep)
    from grid_mapping.shell_model import compute_shell_volume
    from grid_mapping.opacity_tau import compute_tau
    from grid_mapping.solver_metrics import sigma_velocity
    r_kept = r[keep_arr]
    r_bound = _centres_to_boundaries(r_kept, cfg)
    vol = compute_shell_volume(r_bound)
    mass_kept = target.shell_mass_g[keep_arr] if target.shell_mass_g is not None else None
    rho_kept = density_from_mass_volume(mass_kept, vol) if mass_kept is not None else target.density_gcm3[keep_arr]
    kappa_kept = target.kappa_cgs[keep_arr] if len(target.kappa_cgs) == len(r) else np.full(len(keep_arr), 0.34)
    tau_kept = compute_tau(r_bound, rho_kept, kappa_kept)
    v_kept = target.velocity_cms[keep_arr]
    sigma_kept = sigma_velocity(r_kept, v_kept, preserve_sign=cfg.preserve_v_sign)
    merged = TargetGrid(
        r_cm=r_kept, velocity_cms=v_kept,
        temperature_K=target.temperature_K[keep_arr],
        density_gcm3=rho_kept, kappa_cgs=kappa_kept,
        sigma=sigma_kept, tau_cumulative=tau_kept,
        mass_fractions={k: v[keep_arr] for k, v in target.mass_fractions.items()},
        isotope_fractions={k: v[keep_arr] for k, v in
                           (getattr(target, 'isotope_fractions', {}).items() or {})},
        shell_mass_g=mass_kept, shell_volume_cm3=vol,
        time_days=target.time_days,
        provenance=target.provenance[keep_arr] if hasattr(target, 'provenance') and len(target.provenance) == len(r) else np.zeros(len(keep_arr), dtype=int),
        is_anchor=target.is_anchor[keep_arr] if hasattr(target, 'is_anchor') and target.is_anchor is not None and len(target.is_anchor) == len(r) else None,
        electron_density=target.electron_density[keep_arr] if target.electron_density is not None else None,
    )
    return merged


def build_grid(source: SourceShellGrid, cfg: XNewGridConfig,
               outpath: Optional[str] = None,
               lum_lsun: float = 0.0) -> tuple[TargetGrid, dict]:
    """Run the full one-pass X_newgrid construction.

    Parameters
    ----------
    source : un-trimmed SourceShellGrid (features need not be pre-detected)
    cfg    : XNewGridConfig
    outpath: if given, write SN_HYDRO_DATA here
    lum_lsun: luminosity to write into the SN_HYDRO_DATA header

    Returns (target_grid, report).
    """
    report: dict = {'stages': []}

    # 1. Trim + detect features
    src = _trim_source(source, cfg)
    detect_all_features(src)
    report['stages'].append({'stage': 'trim', 'n_zones': src.n_shells})
    if cfg.verbose:
        print(f"[X_newgrid] After trim: {src.n_shells} zones, "
              f"r=[{src.r_centre_cm[0]:.3e},{src.r_centre_cm[-1]:.3e}] cm")

    # 2. Build tau-decade skeleton (Gap B) + merge anchors as boundaries
    r_bound = build_tau_skeleton(src, cfg)
    # Enforce strict decrease with a meaningful minimum spacing relative to
    # the radial span, to avoid duplicate / near-duplicate boundaries in
    # models with a tiny tau span (e.g. uniform-density synthetic tests).
    r_span = r_bound[0] - r_bound[-1]
    min_dr = max(cfg.min_dr_r * r_span / max(len(r_bound), 1),
                 cfg.min_dr_r * r_bound[-1])
    r_bound = _enforce_min_spacing(r_bound, min_dr)
    n_skel_cells = len(r_bound) - 1
    report['stages'].append({'stage': 'skeleton', 'n_cells': n_skel_cells})
    if cfg.verbose:
        print(f"[X_newgrid] Skeleton: {n_skel_cells} cells "
              f"(~{cfg.points_per_tau_decade}/decade, "
              f"{len(r_bound)} boundaries incl. anchors)")

    # 3. Hierarchical-priority refinement (Gap A)
    target, refine_rep = hierarchical_refine(src, r_bound, cfg)
    # Post-refinement: merge any duplicate-radius cells that arose from
    # zero-mass fixup in steep-gradient / optically-thin regions.
    target = _merge_duplicate_cells(target, cfg)
    report['stages'].append({'stage': 'refinement',
                             'iterations': refine_rep['iterations'],
                             'splits': refine_rep['splits'],
                             'n_points': target.n_points,
                             'converged': refine_rep['converged'],
                             'stagnant': refine_rep['stagnant']})
    if cfg.verbose:
        print(f"[X_newgrid] Refinement: {refine_rep['splits']} splits -> "
              f"{target.n_points} points "
              f"(converged={refine_rep['converged']})")

    # 4. Coarsening (Gap E)
    target, coarse_rep = coarsen_grid(target, src, cfg)
    target = _merge_duplicate_cells(target, cfg)
    report['stages'].append({'stage': 'coarsening', **coarse_rep})

    # 5. Local interpolation-error metric (Gap D)
    err = local_interpolation_error(target, src, cfg)
    report['stages'].append({'stage': 'error_metric', **err.to_dict()})
    if cfg.verbose:
        print(f"[X_newgrid] Error metric: passed={err.passed} "
              f"(worst {err.worst_variable}={max(err.per_variable.values()):.4f})")

    # 6. Write SN_HYDRO_DATA
    write_meta = None
    if outpath is not None:
        write_meta = write_sn_hydro_data(
            outpath, target, lum_lsun=lum_lsun,
            t_floor=cfg.t_floor, t_cap=cfg.t_cap,
            trace_threshold=cfg.trace_threshold,
            preserve_v_sign=cfg.preserve_v_sign, vel_law=cfg.vel_law)
        report['output'] = str(outpath)

    report['final_n'] = target.n_points
    report['write_meta'] = write_meta
    return target, report


def build_grid_report(report: dict) -> str:
    """Render a human-readable summary of a build_grid report."""
    lines = ["X_newgrid build summary", "=" * 40]
    for st in report.get('stages', []):
        lines.append(f"  {st.get('stage','?')}: " +
                     " ".join(f"{k}={v}" for k, v in st.items()
                              if k != 'stage'))
    lines.append(f"  final_n = {report.get('final_n','?')}")
    if report.get('write_meta'):
        m = report['write_meta']
        lines.append(f"  output  = {report.get('output','?')} "
                     f"({m['n_points']} pts, {m['n_mass_fractions']} species)")
    return "\n".join(lines)
