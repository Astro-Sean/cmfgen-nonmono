"""Gap E: coarsening / minimisation pass.

After hierarchical refinement converges, sweep through all non-anchor
points and try removing each one.  A removal is accepted only if:

  * the point is NOT an anchor (protected feature);
  * all individual hard limits stay satisfied
    (``coarsen_hard_factor * soft target``);
  * the velocity topology (extrema count) does not degrade;
  * no under-resolved feature appears;
  * the local interpolation error stays within tolerance (Gap D).

This finds the *minimum* viable grid, consistent with the project
AGENTS.md minimal-grid principle: a model using 1000 points when 300
suffice is poorly optimised.
"""

from __future__ import annotations

import numpy as np

from grid_mapping.shell_model import SourceShellGrid, TargetGrid
from grid_mapping.conservative_remap import conservative_remap
from grid_mapping.solver_metrics import count_velocity_reversals

from .error_metric import local_interpolation_error
from .hierarchy import _cell_violations, _feature_violation_mask, PriorityLevel


def _anchor_radii_set(source: SourceShellGrid, cfg) -> set[float]:
    """Set of anchor radii that must never be removed during coarsening."""
    r = source.r_centre_cm
    n = source.n_shells
    out = {float(r[0]), float(r[-1])}

    def _add(mask):
        if mask is not None and len(mask) == n:
            for idx in np.where(mask)[0]:
                out.add(float(r[idx]))

    _add(source.composition_interface_mask)
    _add(source.shock_mask)
    _add(source.contact_discontinuity_mask)
    _add(source.velocity_extremum_mask)
    _add(source.velocity_reversal_mask)
    _add(source.ionization_front_mask)
    _add(source.isotope_boundary_mask)
    return out


def _protected_indices(r_bound: np.ndarray, anchor_radii: set[float],
                       cfg) -> set[int]:
    """Indices of boundary points that coincide with an anchor radius.

    A boundary is protected if it lies within 1% of the local neighbour
    spacing of any anchor radius.  The outer/inner boundaries are always
    protected.
    """
    n = len(r_bound)
    prot = {0, n - 1}
    for i in range(1, n - 1):
        local_dr = abs(r_bound[i - 1] - r_bound[i + 1])
        tol = max(local_dr * 0.01, cfg.min_dr_r * abs(r_bound[i]))
        for ar in anchor_radii:
            if abs(r_bound[i] - ar) < tol:
                prot.add(i)
                break
    return prot


def _individual_hard_limits_ok(target: TargetGrid, cfg) -> bool:
    """Check each criterion stays within ``coarsen_hard_factor * soft target``."""
    n = target.n_points
    if n < 2:
        return True
    hf = cfg.coarsen_hard_factor
    tau = target.tau_cumulative
    dlogtau = np.abs(np.log10((tau[1:] + cfg.tau_floor) /
                              (tau[:-1] + cfg.tau_floor)))
    if float(dlogtau.max()) > cfg.target_dlogtau * hf:
        return False
    rho = target.density_gcm3
    dlogrho = np.abs(np.log10(np.maximum(rho[1:], 1e-30) /
                              np.maximum(rho[:-1], 1e-30)))
    if float(dlogrho.max()) > cfg.target_dlogrho * hf:
        return False
    dv = np.abs(np.diff(target.velocity_cms))
    if float(dv.max()) > cfg.target_dv_kms * 1e5 * hf:
        return False
    T = target.temperature_K
    dlogT = np.abs(np.log10(np.maximum(T[1:], 1.0) /
                            np.maximum(T[:-1], 1.0)))
    if float(dlogT.max()) > cfg.target_dlogT * hf:
        return False
    dr = np.abs(np.diff(target.r_cm))
    r_mid = 0.5 * (target.r_cm[:-1] + target.r_cm[1:])
    dror = dr / np.maximum(r_mid, 1e-30)
    if float(dror.max()) > cfg.target_dr_over_r * hf:
        return False
    return True


def coarsen_grid(target: TargetGrid, source: SourceShellGrid,
                 cfg) -> tuple[TargetGrid, dict]:
    """Remove redundant non-anchor points while preserving all constraints.

    Returns the coarsened TargetGrid and a report dict.
    """
    report = {'points_removed': 0, 'passes': 0, 'attempted': 0, 'final_n': target.n_points}
    if not cfg.enable_coarsening or target.n_points <= cfg.min_points:
        return target, report

    anchor_radii = _anchor_radii_set(source, cfg)
    r_bound = _boundaries_from_target(target)
    protected = _protected_indices(r_bound, anchor_radii, cfg)

    # baseline metrics
    base_err = local_interpolation_error(target, source, cfg)
    base_n_reversals = count_velocity_reversals(target.sigma)
    base_feat = _feature_violation_mask(target, source, cfg)
    base_n_feat = int(base_feat.sum())

    n_removed = 0
    for sweep in range(cfg.coarsen_max_passes):
        n_before = len(r_bound)
        # rank candidate points by local error contribution (lowest first)
        violations = _cell_violations(target, cfg)
        point_score = np.full(len(r_bound), np.inf)
        for i in range(1, len(r_bound) - 1):
            if i in protected:
                continue
            s = 0.0
            if i - 1 < len(violations):
                s = max(s, violations[i - 1].score)
            if i < len(violations):
                s = max(s, violations[i].score)
            point_score[i] = s

        for idx in np.argsort(point_score):
            i = int(idx)
            if i < 1 or i >= len(r_bound) - 1:
                continue
            if i in protected:
                continue
            report['attempted'] += 1

            r_try = np.delete(r_bound, i)
            # min spacing guard
            if abs(r_try[i - 1] - r_try[i]) < \
                    cfg.min_dr_r * min(abs(r_try[i - 1]), abs(r_try[i])):
                continue

            target_try = conservative_remap(source, r_try, source.time_days)
            if not _individual_hard_limits_ok(target_try, cfg):
                continue
            if count_velocity_reversals(target_try.sigma) < base_n_reversals:
                continue
            feat = _feature_violation_mask(target_try, source, cfg)
            if int(feat.sum()) > base_n_feat:
                continue
            err = local_interpolation_error(target_try, source, cfg)
            if not err.passed:
                continue

            # accept removal
            r_bound = r_try
            target = target_try
            protected = _protected_indices(r_bound, anchor_radii, cfg)
            n_removed += 1

        report['passes'] = sweep + 1
        if len(r_bound) == n_before:
            break

    report['points_removed'] = n_removed
    report['final_n'] = target.n_points
    if cfg.verbose and n_removed > 0:
        print(f"[X_newgrid] Coarsening: removed {n_removed} points "
              f"({n_before} -> {target.n_points})")
    return target, report


def _boundaries_from_target(target: TargetGrid) -> np.ndarray:
    """Reconstruct shell boundaries from target centre radii.

    TargetGrid stores centre radii (r_cm).  Boundaries are the geometric
    mean of adjacent centres, with the outer/inner edges extrapolated.
    """
    r = target.r_cm
    n = len(r)
    if n < 2:
        return r.copy()
    r_bound = np.empty(n + 1)
    r_bound[1:-1] = np.sqrt(r[:-1] * r[1:])
    # outer edge: extrapolate using the first cell's ratio
    r_bound[0] = r[0] * (r[0] / r_bound[1]) if r_bound[1] > 0 else r[0] * 1.01
    # inner edge
    r_bound[-1] = r[-1] * (r[-1] / r_bound[-2]) if r_bound[-2] > 0 else r[-1] * 0.99
    # enforce strict decrease
    for i in range(len(r_bound) - 1):
        if r_bound[i] <= r_bound[i + 1]:
            r_bound[i] = r_bound[i + 1] * 1.0001
    return r_bound
