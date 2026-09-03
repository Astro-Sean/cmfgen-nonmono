"""Gap A: explicit hierarchical-priority refinement.

Unlike a weighted-sum AMR score, refinement here follows an *ordered
priority hierarchy* (the Dessart/Hillier methodology surveyed for this
experiment):

    TAU  >  VELOCITY  >  DISCONTINUITY  >  ION_FRONT  >  SMOOTH

A higher-priority violation is always addressed before any lower-priority
one, regardless of magnitude.  Within a single priority level the
normalised violation magnitude breaks ties.  This is the ``S = max[...]``
feature-aware score with an explicit ordering, rather than treating tau,
v, T, rho, X_i as equal-weight AMR variables.

Hard-constraint features (shocks, composition interfaces, velocity
extrema) are *mandatory*: a cell overlapping an under-resolved feature is
refined before ANY smooth-region criterion, regardless of priority level.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import IntEnum
from typing import Optional

import numpy as np

from grid_mapping.shell_model import SourceShellGrid, TargetGrid
from grid_mapping.conservative_remap import conservative_remap


class PriorityLevel(IntEnum):
    """Refinement priority levels, in descending priority."""
    HARD_FEATURE = 0   # mandatory: shocks / interfaces / extrema under-resolved
    TAU = 1            # optical-depth spacing (primary transfer coordinate)
    VELOCITY = 2       # |dv|, |dv|/v, dv/dr extrema
    DISCONTINUITY = 3  # density / T jumps, composition gradients
    ION_FRONT = 4      # ionisation/recombination fronts (Saha estimate)
    SMOOTH = 5         # residual T / rho / geometric smoothness


@dataclass
class CellViolation:
    """A single cell's violation, classified by priority level."""
    level: PriorityLevel
    score: float        # normalised violation magnitude (>=1 means violated)
    criterion: str      # which sub-criterion triggered (for diagnostics)


def _cell_violations(target: TargetGrid, cfg) -> np.ndarray:
    """Return a structured array of per-cell (n-1 cells) violations.

    Each cell gets the *worst* (highest-priority, then largest-score)
    violation.  Returns an object array of CellViolation.
    """
    n = target.n_points
    if n < 2:
        return np.array([], dtype=object)
    r = target.r_cm
    v = target.velocity_cms
    T = target.temperature_K
    rho = target.density_gcm3
    tau = target.tau_cumulative

    dr = np.abs(np.diff(r))
    r_mid = 0.5 * (r[:-1] + r[1:])

    # Only score cells in the radiatively-important region (tau > floor).
    # The optically-thin outer region is covered by the skeleton and is not
    # over-refined; tiny zero-mass cells there are numerical artifacts.
    tau_mid = 0.5 * (tau[:-1] + tau[1:])
    radiative = tau_mid >= cfg.tau_metric_min

    # --- TAU (priority 1) ---
    dlogtau = np.abs(np.log10((tau[1:] + cfg.tau_floor) /
                              (tau[:-1] + cfg.tau_floor)))
    s_tau = dlogtau / cfg.target_dlogtau

    # --- VELOCITY (priority 2) ---
    dv = np.abs(np.diff(v))
    v_mid = 0.5 * (np.abs(v[:-1]) + np.abs(v[1:]))
    s_dv = dv / (cfg.target_dv_kms * 1e5)          # km/s -> cm/s
    s_dvov = dv / np.maximum(v_mid, 1.0) / cfg.target_dv_over_v
    s_vel = np.maximum(s_dv, s_dvov)

    # --- DISCONTINUITY (priority 3) ---
    dlogrho = np.abs(np.log10(np.maximum(rho[1:], 1e-30) /
                              np.maximum(rho[:-1], 1e-30)))
    s_rho = dlogrho / cfg.target_dlogrho
    drho_frac = np.abs(np.diff(rho)) / np.maximum(rho[:-1], 1e-30)
    s_djump = drho_frac / cfg.max_density_jump
    # composition gradient (max over species)
    s_comp = np.zeros(n - 1)
    for frac in target.mass_fractions.values():
        x = np.asarray(frac)
        if len(x) == n:
            dx = np.abs(np.diff(x))
            s_comp = np.maximum(s_comp, dx / max(cfg.composition_jump, 1e-30))
    s_disc = np.maximum(np.maximum(s_rho, s_djump), s_comp)

    # --- ION_FRONT (priority 4) ---
    # On the target grid we only have the (Saha) electron density carried
    # through the remap; use it as the front proxy.  Real CMFGEN fronts are
    # handled by the standalone two-pass script.
    s_ion = np.zeros(n - 1)
    if target.electron_density is not None:
        ne = np.maximum(target.electron_density, 1e-30)
        dne = np.abs(np.log10(ne[1:] / ne[:-1]))
        s_ion = dne / cfg.ne_log_threshold

    # --- SMOOTH (priority 5) ---
    dlogT = np.abs(np.log10(np.maximum(T[1:], 1.0) /
                            np.maximum(T[:-1], 1.0)))
    s_T = dlogT / cfg.target_dlogT
    dT_frac = np.abs(np.diff(T)) / np.maximum(T[:-1], 1.0)
    s_Tjump = dT_frac / cfg.max_temp_jump
    dror = dr / np.maximum(r_mid, 1e-30)
    s_r = dror / cfg.target_dr_over_r
    # cell-size ratio
    s_ratio = np.zeros(n - 1)
    if n >= 3:
        ratios = np.maximum(dr[:-1] / np.maximum(dr[1:], 1e-30),
                            dr[1:] / np.maximum(dr[:-1], 1e-30))
        s_ratio[1:] = ratios / cfg.max_cell_ratio
    s_smooth = np.maximum(np.maximum(s_T, s_Tjump), np.maximum(s_r, s_ratio))

    # --- Assemble worst violation per cell (priority order) ---
    out = np.empty(n - 1, dtype=object)
    for i in range(n - 1):
        # Skip non-radiative (optically thin) cells: no refinement needed.
        if not radiative[i]:
            out[i] = CellViolation(PriorityLevel.SMOOTH, 0.0, "thin")
            continue
        candidates = [
            (PriorityLevel.TAU, s_tau[i], "dlogtau"),
            (PriorityLevel.VELOCITY, s_vel[i],
             "dv" if s_dv[i] >= s_dvov[i] else "dv/v"),
            (PriorityLevel.DISCONTINUITY, s_disc[i],
             "rho" if s_rho[i] >= max(s_djump[i], s_comp[i]) else "comp"),
            (PriorityLevel.ION_FRONT, s_ion[i], "dlogNe"),
            (PriorityLevel.SMOOTH, s_smooth[i], "T/r/geom"),
        ]
        # pick highest-priority violated criterion; within level, largest score
        best = None
        for lvl, sc, crit in candidates:
            if sc < 1.0:
                continue
            if best is None or lvl < best.level or \
                    (lvl == best.level and sc > best.score):
                best = CellViolation(lvl, float(sc), crit)
        if best is None:
            # no violation: use the highest-priority sub-1.0 score for tie-break
            best = CellViolation(PriorityLevel.SMOOTH, 0.0, "ok")
        out[i] = best
    return out


def _feature_violation_mask(target: TargetGrid, source: SourceShellGrid,
                            cfg) -> np.ndarray:
    """Boolean mask (n-1 cells) of cells overlapping an under-resolved feature.

    A feature is under-resolved if fewer than its ``min_points`` grid points
    fall within a local radius window around the feature location.
    """
    n = target.n_points
    mask = np.zeros(n - 1, dtype=bool)
    r_mid = 0.5 * (target.r_cm[:-1] + target.r_cm[1:])

    def _check(mask_arr, min_pts, kind):
        nonlocal mask
        if mask_arr is None or len(mask_arr) != source.n_shells:
            return
        for idx in np.where(mask_arr)[0]:
            r_feat = source.r_centre_cm[idx]
            # window in log-radius around the feature
            dlog = np.abs(np.log10(np.maximum(r_mid, 1e-30)) -
                          np.log10(max(r_feat, 1e-30)))
            n_near = int(np.sum(dlog < 0.05))
            if n_near < min_pts:
                # mark cells in the window as needing refinement
                mask = mask | (dlog < 0.05)

    _check(source.shock_mask, cfg.min_points_per_shock, "shock")
    _check(source.composition_interface_mask,
           cfg.min_points_per_interface, "interface")
    _check(source.velocity_extremum_mask,
           cfg.min_points_per_extremum, "extremum")
    _check(source.velocity_reversal_mask,
           cfg.min_points_per_extremum, "reversal")
    _check(source.ionization_front_mask,
           cfg.min_points_per_ionfront, "ionfront")
    return mask


def _worst_cell(violations: np.ndarray, feat_mask: np.ndarray
                ) -> tuple[int, CellViolation]:
    """Pick the cell to split: hard-feature cells first, then by priority."""
    n = len(violations)
    if n == 0:
        return -1, CellViolation(PriorityLevel.SMOOTH, 0.0, "none")

    # Hard-feature cells take absolute precedence.
    if feat_mask.any():
        cand = np.where(feat_mask)[0]
        # within feature cells, pick the one with the worst (lowest-priority-
        # number, then largest-score) violation
        best_i = -1
        best = None
        for i in cand:
            cv = violations[i]
            if best is None or cv.level < best.level or \
                    (cv.level == best.level and cv.score > best.score):
                best, best_i = cv, i
        return best_i, best

    # Otherwise pick the highest-priority, then largest-score violation.
    best_i = -1
    best = None
    for i in range(n):
        cv = violations[i]
        if cv.score < 1.0 and cv.level == PriorityLevel.SMOOTH:
            continue
        if best is None or cv.level < best.level or \
                (cv.level == best.level and cv.score > best.score):
            best, best_i = cv, i
    if best is None:
        return -1, CellViolation(PriorityLevel.SMOOTH, 0.0, "converged")
    return best_i, best


def _ranked_candidates(violations: np.ndarray, feat_mask: np.ndarray):
    """Return violated cells sorted by priority: feature cells first, then
    by (priority level asc, score desc).  Only cells with score >= 1.0 or
    on the feature mask are included."""
    out = []
    for i, cv in enumerate(violations):
        on_feat = bool(feat_mask[i]) if i < len(feat_mask) else False
        if on_feat or cv.score >= 1.0:
            out.append((i, cv))
    # feature cells first, then by level, then by score
    out.sort(key=lambda t: (0 if feat_mask[t[0]] else 1,
                            int(t[1].level), -t[1].score))
    return out


def _insert_point(r_bound: np.ndarray, edge_idx: int) -> np.ndarray:
    """Insert a new boundary at the geometric-mean radius of cell edge_idx.

    ``r_bound`` is decreasing; cell edge_idx spans
    [r_bound[edge_idx+1], r_bound[edge_idx]].  The new boundary splits it.
    Falls back to the arithmetic mean if the geometric mean is degenerate.
    """
    r_lo = r_bound[edge_idx + 1]
    r_hi = r_bound[edge_idx]
    r_new = np.sqrt(r_lo * r_hi)
    if not (r_lo < r_new < r_hi):
        r_new = 0.5 * (r_lo + r_hi)
    return np.insert(r_bound, edge_idx + 1, r_new)


def _cell_too_thin(r_bound: np.ndarray, edge_idx: int,
                   source: SourceShellGrid, cfg) -> bool:
    """Guard: would splitting this cell create a sub-source-resolution cell?

    Splitting a cell thinner than the local source cell width produces
    zero-mass target zones (the conservative remap has no source mass to
    distribute), whose neighbour-interpolated density is numerically
    meaningless.  Refusing such splits prevents the dlogrho blow-up seen
    in steep-gradient / optically-thin regions.
    """
    r_mid = 0.5 * (r_bound[edge_idx] + r_bound[edge_idx + 1])
    # local source cell width at this radius
    idx_src = np.argmin(np.abs(source.r_centre_cm - r_mid))
    if idx_src <= 0 or idx_src >= source.n_shells - 1:
        return False
    dr_src = abs(source.r_centre_cm[idx_src - 1] - source.r_centre_cm[idx_src + 1])
    dr_cell = abs(r_bound[edge_idx] - r_bound[edge_idx + 1])
    # refuse if the cell is already smaller than ~2x the source resolution
    return dr_cell < 0.5 * dr_src


def hierarchical_refine(source: SourceShellGrid,
                        r_bound_init: np.ndarray,
                        cfg) -> tuple[TargetGrid, dict]:
    """Run the hierarchical-priority constrained refinement loop.

    Parameters
    ----------
    source : trimmed SourceShellGrid (features already detected)
    r_bound_init : (n+1,) initial target boundaries (decreasing), e.g. from
        :func:`build_tau_skeleton` converted to boundaries.
    cfg : XNewGridConfig

    Returns
    -------
    target : TargetGrid (conservatively remapped onto the refined grid)
    report : dict with iteration diagnostics
    """
    r_bound = r_bound_init.copy()
    target = conservative_remap(source, r_bound, source.time_days)

    report = {
        'iterations': 0,
        'splits': 0,
        'history': [],
        'final_n': target.n_points,
        'converged': False,
        'stagnant': False,
    }
    max_iter = 4 * cfg.max_points
    best_worst_score = float('inf')
    stagnant = 0

    for it in range(max_iter):
        report['iterations'] = it
        if target.n_points >= cfg.max_points:
            break

        violations = _cell_violations(target, cfg)
        feat_mask = _feature_violation_mask(target, source, cfg)

        # Collect candidate cells in priority order: feature cells first,
        # then by (priority level asc, score desc).  Try each in turn until
        # one is splittable (not below source resolution).  This avoids
        # re-picking a refused cell every iteration.
        candidates = _ranked_candidates(violations, feat_mask)
        if not candidates:
            report['converged'] = True
            break

        # Stop if the best candidate is not actually violated.
        best_edge, best_cv = candidates[0]
        if (not feat_mask.any() and best_cv.score < 1.0 and
                target.n_points >= cfg.min_points):
            report['converged'] = True
            break

        # Stagnation tracking on the best candidate's score.
        if best_cv.score < best_worst_score * 0.999:
            best_worst_score = best_cv.score
            stagnant = 0
        else:
            stagnant += 1

        # Find a splittable candidate.
        chose = None
        for edge, cv in candidates:
            if _cell_too_thin(r_bound, edge, source, cfg):
                continue
            chose = (edge, cv)
            break

        if chose is None:
            # No violated cell can be split further (all at source-res limit).
            if stagnant >= 5 and target.n_points >= cfg.min_points:
                report['stagnant'] = True
                if cfg.verbose:
                    print(f"[X_newgrid] Stagnation at {target.n_points} pts "
                          f"(all violated cells at source-resolution limit)")
                break
            stagnant += 1
            continue

        edge, cv = chose
        worst_on_feature = bool(feat_mask[edge]) if 0 <= edge < len(feat_mask) else False
        if stagnant >= 40 and worst_on_feature and \
                target.n_points >= cfg.min_points:
            report['stagnant'] = True
            if cfg.verbose:
                print(f"[X_newgrid] Stagnation at {target.n_points} pts "
                      f"(worst {cv.criterion} score={cv.score:.3f} on feature)")
            break

        r_bound = _insert_point(r_bound, edge)
        target = conservative_remap(source, r_bound, source.time_days)
        report['splits'] += 1
        report['history'].append({
            'iter': it, 'edge': int(edge),
            'level': cv.level.name, 'criterion': cv.criterion,
            'score': cv.score, 'n': target.n_points,
        })

    report['final_n'] = target.n_points
    return target, report
