"""Gap B: optical-depth-decade-based skeleton construction.

The skeleton is sized from the physics, not a flat point count:

    N_skeleton ~ points_per_tau_decade * n_tau_decades

where ``n_tau_decades`` is the number of optical-depth decades spanned by the
trimmed source grid.  This follows the Dessart/Hillier CMFGEN practice of
~5-7 points per optical-depth decade.

Mandatory anchor radii (features) are inserted as exact cell **boundaries**
— not centres — so that no cell straddles a physical discontinuity
(composition interface, shock, opacity / ionisation front).  This is the
key point from the surveyed strategy: a discontinuity that becomes a
boundary cannot be numerically washed out by conservative remapping,
because each adjacent cell lies entirely on one side of the jump.
"""

from __future__ import annotations

import numpy as np

from grid_mapping.shell_model import SourceShellGrid


def n_tau_decades(tau: np.ndarray, cfg_tau_floor: float = 1.0e-6) -> float:
    """Number of optical-depth decades spanned by ``tau``."""
    tau = np.asarray(tau, dtype=float)
    if tau.size == 0:
        return 0.0
    log_tau = np.log10(tau + cfg_tau_floor)
    return float(log_tau.max() - log_tau.min())


def skeleton_point_count(tau: np.ndarray, cfg) -> int:
    """Number of log-tau skeleton *boundaries* before anchors are added."""
    ndec = n_tau_decades(tau, cfg.tau_floor)
    n = int(round(cfg.points_per_tau_decade * ndec))
    return max(cfg.skeleton_min, min(n, cfg.target_points))


def build_tau_skeleton(source: SourceShellGrid, cfg) -> np.ndarray:
    """Build the initial target grid **boundaries** (decreasing R).

    Strategy:
    1. Place ``points_per_tau_decade * n_decades`` boundaries uniformly in
       ``log10(tau)`` and map them to radius via the source tau(R).
    2. Insert every anchor radius (composition interfaces, shocks, velocity
       extrema, dv/dr sign changes, ionisation fronts, inner/outer edges)
       as an EXACT boundary so no cell straddles a discontinuity.
    3. Enforce strictly decreasing, unique boundaries with ``min_dr_r``.

    Returns the boundary radii (n_bound decreasing), ready for
    ``conservative_remap``.  Anchors are boundaries, not centres.
    """
    n_src = source.n_shells
    if n_src < 2:
        # degenerate: single boundary pair
        rb = source.r_bound_cm.copy()
        return rb

    tau = np.maximum(source.tau_cumulative, 0.0)
    # Source boundaries' tau: integrate to each boundary.  Approximate by
    # interpolating the centre tau onto the boundary radii (outer->inner).
    # tau increases inward; boundaries lie between centres.
    if tau.size == n_src and source.r_bound_cm.size == n_src + 1:
        # tau at boundaries via linear interp in radius (reversed for interp)
        r_inc = source.r_centre_cm[::-1]
        tau_inc = tau[::-1]
        tau_bound = np.interp(source.r_bound_cm, r_inc, tau_inc)
    else:
        n_skel = cfg.skeleton_min
        return np.geomspace(source.r_bound_cm[0],
                            source.r_bound_cm[-1], n_skel + 1)

    # --- Step 1: log-tau-uniform boundaries ---
    n_skel = skeleton_point_count(tau, cfg)
    log_tau = np.log10(tau_bound + cfg.tau_floor)
    lo, hi = log_tau.min(), log_tau.max()
    if hi <= lo:
        r_skel = np.geomspace(source.r_bound_cm[0],
                              source.r_bound_cm[-1], n_skel + 1)
    else:
        targets = np.linspace(lo, hi, n_skel + 1)
        tau_targets = np.power(10.0, targets) - cfg.tau_floor
        tau_targets = np.maximum(tau_targets, 0.0)
        # tau_bound increases inward (r decreasing); interp needs increasing x
        r_skel = np.interp(tau_targets, tau_bound, source.r_bound_cm)

    # --- Step 2: insert anchor radii as exact boundaries ---
    anchor_r = _anchor_radii(source, cfg)
    r_all = np.unique(np.concatenate([r_skel, anchor_r]))[::-1]  # decreasing

    # --- Step 3: enforce strict decrease + min spacing ---
    r_all = _enforce_strict_decrease(r_all, cfg.min_dr_r)

    # Pin exact outer/inner source boundaries
    r_all[0] = source.r_bound_cm[0]
    r_all[-1] = source.r_bound_cm[-1]
    r_all = _enforce_strict_decrease(r_all, cfg.min_dr_r)
    return r_all


def _anchor_radii(source: SourceShellGrid, cfg) -> np.ndarray:
    """Collect every mandatory anchor radius from the source feature masks.

    These become exact cell boundaries and are NEVER removed during
    coarsening.  We include:
      * composition interfaces, shocks, contact discontinuities
      * velocity extrema and dv/dr sign-change (reversal) locations
      * ionisation-front locations (Saha estimate)
      * isotope boundaries
    plus the inner/outer boundaries.
    """
    r = source.r_centre_cm
    n = source.n_shells
    radii = [float(source.r_bound_cm[0]), float(source.r_bound_cm[-1])]

    def _add(mask):
        if mask is not None and len(mask) == n:
            for idx in np.where(mask)[0]:
                radii.append(float(r[idx]))

    _add(source.composition_interface_mask)
    _add(source.shock_mask)
    _add(source.contact_discontinuity_mask)
    _add(source.velocity_extremum_mask)
    _add(source.velocity_reversal_mask)
    _add(source.ionization_front_mask)
    _add(source.isotope_boundary_mask)

    # Every dv/dr sign change directly from V(R), so even shallow extrema the
    # significance-gated detector skips are preserved as exact boundaries.
    radii.extend(_dvdr_sign_change_radii(source, cfg))
    return np.array(sorted(set(radii)), dtype=float)


def _dvdr_sign_change_radii(source: SourceShellGrid, cfg) -> list[float]:
    """Radii where dV/dR changes sign (velocity extrema / reversals).

    A sign change of the discrete velocity gradient marks a local extremum.
    These are mandatory boundaries: without one at the extremum, PCHIP
    interpolation can shift the velocity gradient substantially even when
    |dv| between the two bracketing cells looks acceptable.
    """
    r = source.r_centre_cm
    v = source.velocity_cms
    n = source.n_shells
    if n < 3:
        return []
    out = []
    for i in range(1, n - 1):
        dv_left = v[i] - v[i - 1]
        dv_right = v[i + 1] - v[i]
        if dv_left == 0.0 or dv_right == 0.0:
            continue
        if (dv_left > 0) != (dv_right > 0):
            # local extremum at i; only keep if significant
            if abs(v[i] - 0.5 * (v[i - 1] + v[i + 1])) > \
                    cfg.extremum_significance_kms * 1e5:  # km/s -> cm/s
                out.append(float(r[i]))
    return out


def _enforce_strict_decrease(r: np.ndarray, min_dr_r: float) -> np.ndarray:
    """Ensure strictly decreasing radius with minimum fractional spacing."""
    r = np.asarray(r, dtype=float)
    out = [r[0]]
    for val in r[1:]:
        if val >= out[-1]:
            val = out[-1] * (1.0 - max(min_dr_r, 1.0e-8))
        else:
            if (out[-1] - val) < min_dr_r * abs(out[-1]):
                val = out[-1] * (1.0 - min_dr_r)
        out.append(val)
    return np.array(out)
