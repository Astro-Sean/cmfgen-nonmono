"""
equidistribute_grid.py

A single-pass, monitor-function-based replacement for the ~10-step
add/remove grid-construction chain in convert_stella_swd_to_cmfgen.py
(stretch_critical_regions -> downsample -> adaptive_grid_refinement ->
thin_degenerate_grid -> insert_outer_boundary_points -> fill_grid_to_minimum
-> enforce_grid_spacing -> enforce_sigma_conformance ->
coarsen_negative_sigma_regions -> enforce_velocity_conformance ->
refine_outer_boundary_grid).

WHY THIS EXISTS
---------------
The current pipeline gets its final point count and shape as an emergent
property of many local passes, each blind to the others. Concretely:
`enforce_sigma_conformance` (insert-only, meant to run last) is followed
immediately by `coarsen_negative_sigma_regions`, which deliberately
collapses negative-sigma (compression) zones back down to 1-2 points --
undoing the previous step for exactly the regions that matter most, and
producing the extreme |sigma| values (seen as low as -72 in practice)
that make the region numerically pathological downstream in CMFGEN.

This module instead:
  1. Builds ONE monitor function on the full-resolution STELLA grid that
     combines every resolution requirement (optical depth spacing,
     velocity-gradient curvature/sigma, temperature and density
     gradients, composition interfaces) with EXPLICIT, named weights.
  2. Equidistributes that monitor function to select exactly N target
     points in one shot -- no downstream pass can silently undo the
     choice made here, because there are no more passes.
  3. Treats "how much should compression zones be resolved" as ONE
     named policy switch (`compression_policy`), not an emergent
     property of function call order.
  4. Validates the result (point count, |sigma|, minimum spacing) and
     raises loudly instead of silently writing a bad grid, iterating a
     bounded number of times to try to fix violations before giving up.

This is a draft: the monitor-function weights and the exact interpolation
choices (linear vs log, which fields need mass-conserving remapping vs.
simple point-sampling) should be tuned against real STELLA snapshots
before this replaces the existing pipeline in production. Wire it in by
replacing the block from the initial downsample call through
`refine_outer_boundary_grid()` in `main()` with a single call to
`build_target_grid()`.
"""

from __future__ import annotations

import numpy as np


# ---------------------------------------------------------------------------
# Core sigma formula, kept identical to convert_stella_swd_to_cmfgen.py's
# sigma_velocity() (bit-for-bit CMFGEN match) so validation against it is
# meaningful. If you wire this module into the real converter, delete this
# copy and import sigma_velocity from there instead.
# ---------------------------------------------------------------------------
def sigma_velocity(r, v_kms):
    r = np.asarray(r, dtype=float)
    v = np.asarray(v_kms, dtype=float)
    n = len(r)
    sigma = np.zeros(n)
    eps = 1e-30
    if n < 2:
        return sigma
    sigma[0] = r[0] * (v[0] - v[1]) / (r[0] - r[1]) / max(abs(v[0]), eps) - 1.0
    sigma[-1] = r[-1] * (v[-2] - v[-1]) / (r[-2] - r[-1]) / max(abs(v[-1]), eps) - 1.0
    if n < 3:
        return sigma
    h_left = r[:-2] - r[1:-1]
    h_right = r[1:-1] - r[2:]
    with np.errstate(divide="ignore", invalid="ignore"):
        T1 = (v[:-2] - v[1:-1]) / h_left
        T2 = (v[1:-1] - v[2:]) / h_right
        dVdR = (T1 * h_right + T2 * h_left) / (h_left + h_right)
        v_mid = np.maximum(np.abs(v[1:-1]), eps)
        sigma_inner = r[1:-1] * dVdR / v_mid - 1.0
    degenerate = (np.abs(h_left) < eps) | (np.abs(h_right) < eps)
    sigma[1:-1] = np.where(degenerate, 0.0, sigma_inner)
    return sigma


# ---------------------------------------------------------------------------
# Monitor function
# ---------------------------------------------------------------------------
def _normalize(w, floor_percentile=50.0):
    """Scale a non-negative weight array to O(1) using a robust reference
    (a percentile rather than the max, so a single spike doesn't wash out
    every other criterion's contribution)."""
    w = np.asarray(w, dtype=float)
    ref = np.percentile(w, floor_percentile)
    if not np.isfinite(ref) or ref <= 0:
        ref = np.max(w) if np.max(w) > 0 else 1.0
    return w / ref


def _smooth(w, width):
    """Light moving-average smoothing so a single-cell gradient spike
    doesn't pull the entire point budget onto one cell while leaving its
    neighbors unresolved. width is a half-window in points."""
    if width <= 0:
        return w
    n = len(w)
    out = np.empty_like(w)
    for i in range(n):
        lo = max(0, i - width)
        hi = min(n, i + width + 1)
        out[i] = np.mean(w[lo:hi])
    return out


def compute_monitor_function(
    r, v, T, rho, tau, mass_fractions, kappa,
    weights=None,
    compression_policy="resolve",
    compression_threshold=-0.1,
    compression_boost=8.0,
    smooth_width=3,
    base_floor=0.05,
):
    """Build the combined monitor function M(x) on the full-resolution
    STELLA grid, where x = ln(r).

    Parameters
    ----------
    r, v, T, rho, tau, kappa : 1-D arrays, same length, any radius order
        (ascending or descending both work; this function is order-agnostic
        except that r must be monotonic).
    mass_fractions : dict[str, array] of the same length as r.
    weights : dict, optional
        Named weights for each criterion. Keys: 'tau', 'sigma', 'temperature',
        'density', 'composition'. Missing keys default to 1.0. Set a weight
        to 0.0 to disable that criterion entirely.
    compression_policy : {'resolve', 'collapse', 'neutral'}
        'resolve'  -- boost the monitor function in compression zones
                      (sigma < compression_threshold) so equidistribution
                      gives them MORE points than a smooth region of
                      comparable size would get. Use this if the JFNK/
                      non-monotonic solver is expected to do real physics
                      through the reversal (e.g. if FIX_T is relaxed there).
        'collapse' -- suppress the monitor function in compression zones so
                      they get close to the minimum (2-3 points), matching
                      the current pipeline's Dessart-style behaviour. Only
                      valid physically if FIX_T=T is guaranteed through
                      those zones.
        'neutral'  -- do not treat compression zones specially; resolve
                      them the same way as any other region with the same
                      |dsigma/dx|.
    base_floor : float
        Minimum monitor value everywhere, as a fraction of the median
        monitor value, so no region (however "boring") gets zero points
        and the outer/inner boundaries are never starved.

    Returns
    -------
    x : array, ln(r), same order as input r.
    M : array, the combined monitor function, same order as input r,
        strictly positive everywhere.
    diagnostics : dict of the individual (normalized) components, for
        plotting/debugging.
    """
    r = np.asarray(r, dtype=float)
    v = np.asarray(v, dtype=float)
    T = np.asarray(T, dtype=float)
    rho = np.asarray(rho, dtype=float)
    tau = np.asarray(tau, dtype=float)

    if weights is None:
        weights = {}
    w_tau = weights.get("tau", 1.0)
    w_sigma = weights.get("sigma", 1.0)
    w_T = weights.get("temperature", 0.5)
    w_rho = weights.get("density", 0.5)
    w_comp = weights.get("composition", 1.5)

    n = len(r)
    x = np.log(r)

    # dx between consecutive points (sign-agnostic; we only ever use |dx|
    # for arc length, so ascending or descending r both work).
    dx = np.empty(n)
    dx[1:-1] = 0.5 * np.abs(x[2:] - x[:-2])
    dx[0] = np.abs(x[1] - x[0])
    dx[-1] = np.abs(x[-1] - x[-2])
    dx = np.maximum(dx, 1e-12)

    # --- tau criterion: resolve log(tau) changes ---
    tau_floor = max(np.min(tau[tau > 0]) * 1e-3, 1e-10) if np.any(tau > 0) else 1e-10
    ln_tau = np.log(np.maximum(tau, tau_floor))
    d_ln_tau = np.gradient(ln_tau, x)
    m_tau = _normalize(np.abs(d_ln_tau))

    # --- sigma criterion: resolve curvature of the velocity law ---
    sigma = sigma_velocity(r, v)
    d_sigma = np.gradient(sigma, x)
    m_sigma = _normalize(np.abs(d_sigma))

    compression_mask = sigma < compression_threshold
    if compression_policy == "resolve":
        m_sigma = np.where(compression_mask, m_sigma * compression_boost, m_sigma)
    elif compression_policy == "collapse":
        m_sigma = np.where(compression_mask, m_sigma / compression_boost, m_sigma)
    elif compression_policy == "neutral":
        pass
    else:
        raise ValueError(f"Unknown compression_policy: {compression_policy!r}")

    # --- thermodynamic criteria ---
    d_ln_T = np.gradient(np.log(np.maximum(T, 1.0)), x)
    m_T = _normalize(np.abs(d_ln_T))
    d_ln_rho = np.gradient(np.log(np.maximum(rho, 1e-30)), x)
    m_rho = _normalize(np.abs(d_ln_rho))

    # --- composition-interface criterion ---
    if mass_fractions:
        comp_grad = np.zeros(n)
        for _, mf in mass_fractions.items():
            mf = np.asarray(mf, dtype=float)
            if len(mf) != n:
                continue
            comp_grad = np.maximum(comp_grad, np.abs(np.gradient(mf, x)))
        m_comp = _normalize(comp_grad)
    else:
        m_comp = np.zeros(n)

    M = (
        base_floor
        + w_tau * m_tau
        + w_sigma * m_sigma
        + w_T * m_T
        + w_rho * m_rho
        + w_comp * m_comp
    )
    M = _smooth(M, smooth_width)
    M = np.maximum(M, base_floor)

    diagnostics = {
        "sigma": sigma,
        "m_tau": m_tau,
        "m_sigma": m_sigma,
        "m_T": m_T,
        "m_rho": m_rho,
        "m_comp": m_comp,
        "compression_mask": compression_mask,
    }
    return x, M, diagnostics


# ---------------------------------------------------------------------------
# Equidistribution
# ---------------------------------------------------------------------------
def equidistribute(x, M, n_target):
    """Select n_target points x_k such that the monitor-weighted arc
    length between consecutive x_k is constant. Returns x_k in the same
    order (ascending or descending) as the input x.

    x, M must be arrays of the same length, x monotonic (either
    direction), M strictly positive.
    """
    x = np.asarray(x, dtype=float)
    M = np.asarray(M, dtype=float)
    n = len(x)
    if n_target < 2:
        raise ValueError("n_target must be >= 2")
    if n_target > n:
        raise ValueError(
            f"n_target ({n_target}) exceeds the number of source points ({n}); "
            "cannot equidistribute onto a finer grid than the input."
        )

    descending = x[0] > x[-1]
    if descending:
        x_work = x[::-1]
        M_work = M[::-1]
    else:
        x_work = x
        M_work = M

    # Cumulative monitor-weighted arc length (trapezoidal).
    dx = np.diff(x_work)
    seg = 0.5 * (M_work[:-1] + M_work[1:]) * dx
    W = np.concatenate([[0.0], np.cumsum(seg)])
    W_total = W[-1]
    if not np.isfinite(W_total) or W_total <= 0:
        raise ValueError("Monitor function integral is non-positive or non-finite.")

    targets = np.linspace(0.0, W_total, n_target)
    # W is monotonically increasing (M > 0 everywhere) -> safe to invert
    # with a plain interpolant.
    x_new = np.interp(targets, W, x_work)
    # Guard against any interpolation round-off letting two targets map to
    # the same x (can happen if M has a near-zero stretch); nudge apart.
    x_new = np.maximum.accumulate(x_new + np.arange(n_target) * 1e-14)

    if descending:
        x_new = x_new[::-1]
    return x_new


# ---------------------------------------------------------------------------
# Interpolation of physical fields onto the new grid
# ---------------------------------------------------------------------------
def _interp_ordered(x_src, y_src, x_new, log=False):
    """np.interp requires ascending x; handle either order transparently."""
    if x_src[0] > x_src[-1]:
        x_src = x_src[::-1]
        y_src = y_src[::-1]
    if log:
        y = np.exp(np.interp(x_new, x_src, np.log(np.maximum(y_src, 1e-300))))
    else:
        y = np.interp(x_new, x_src, y_src)
    return y


def interpolate_fields(x_src, x_new, v, T, rho, tau, kappa, n_e, mass_fractions):
    """Map every physical field from the full-resolution grid onto the
    new grid. Velocity, temperature, density, tau, kappa, n_e are
    interpolated in log-space where positivity makes that meaningful;
    velocity and mass fractions are interpolated linearly (velocity can
    be zero/near-zero at a reversal, and log-interpolating mass fractions
    that cross zero is not meaningful).

    NOTE: tau here is a fast approximation (log-interpolated from the
    fine grid). For a final, production-quality grid, recompute tau on
    the new grid directly via compute_kappa_tau(r_new, rho_new, kappa_new)
    instead of interpolating it, since tau is itself a definite integral
    along the (now different) grid.
    """
    v_new = _interp_ordered(x_src, v, x_new, log=False)
    T_new = _interp_ordered(x_src, T, x_new, log=True)
    rho_new = _interp_ordered(x_src, rho, x_new, log=True)
    tau_new = _interp_ordered(x_src, tau, x_new, log=True)
    kappa_new = _interp_ordered(x_src, kappa, x_new, log=True)
    n_e_new = _interp_ordered(x_src, n_e, x_new, log=True) if n_e is not None else None

    mf_new = {}
    for species, mf in mass_fractions.items():
        mf_new[species] = np.clip(_interp_ordered(x_src, mf, x_new, log=False), 0.0, 1.0)
    if mf_new:
        total = sum(mf_new.values())
        total = np.where(total > 0, total, 1.0)
        for species in mf_new:
            mf_new[species] = mf_new[species] / total  # renormalize to sum to 1

    return v_new, T_new, rho_new, tau_new, kappa_new, n_e_new, mf_new


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------
def validate_grid(r_new, v_new, n_target, max_abs_sigma=2.0, min_dr_r=1e-6,
                   compression_threshold=-0.1, compression_policy="resolve"):
    """Check the candidate grid against hard acceptance criteria. Returns
    (ok, report_dict). Does not raise -- the caller decides whether to
    iterate or fail.
    """
    report = {}
    report["n_points"] = len(r_new)
    report["n_target"] = n_target
    report["point_count_ok"] = len(r_new) == n_target

    sigma_new = sigma_velocity(r_new, v_new)
    report["max_abs_sigma"] = float(np.max(np.abs(sigma_new)))
    if compression_policy == "collapse":
        # Under the collapse policy we only bound sigma OUTSIDE
        # deliberately-collapsed zones; a collapsed reversal is expected
        # to have a large |sigma| by construction.
        mask = sigma_new >= compression_threshold
        worst = np.max(np.abs(sigma_new[mask])) if np.any(mask) else 0.0
        report["max_abs_sigma_outside_compression"] = float(worst)
        report["sigma_ok"] = worst < max_abs_sigma
    else:
        report["sigma_ok"] = report["max_abs_sigma"] < max_abs_sigma

    dr_r = np.abs(np.diff(r_new)) / np.maximum(np.abs(r_new[:-1]), 1e-300)
    report["min_dr_r"] = float(np.min(dr_r))
    report["spacing_ok"] = report["min_dr_r"] > min_dr_r

    report["monotonic_ok"] = bool(
        np.all(np.diff(r_new) > 0) or np.all(np.diff(r_new) < 0)
    )

    ok = (
        report["point_count_ok"]
        and report["sigma_ok"]
        and report["spacing_ok"]
        and report["monotonic_ok"]
    )
    return ok, report


# ---------------------------------------------------------------------------
# Top-level driver
# ---------------------------------------------------------------------------
def build_target_grid(
    r, v, T, rho, tau, kappa, mass_fractions, n_e=None,
    n_target=100,
    compression_policy="resolve",
    compression_threshold=-0.1,
    compression_boost=8.0,
    weights=None,
    max_abs_sigma=2.0,
    min_dr_r=1e-6,
    max_iterations=4,
    boost_growth=2.0,
    verbose=True,
):
    """Single-pass replacement for the stretch/downsample/AMR/thin/
    insert/fill/enforce-spacing/enforce-sigma/coarsen/enforce-velocity/
    refine-outer-boundary chain.

    Iterates a bounded number of times, growing `compression_boost` (if
    compression_policy == 'resolve') when the candidate grid still
    violates the sigma bound, since equidistribution with a fixed weight
    can undershoot in a very sharp, very localized reversal. Raises
    ValueError if it still fails validation after max_iterations --
    by design, this refuses to silently emit a bad grid, unlike the
    current pipeline.

    Returns
    -------
    dict with keys: r, v, T, rho, tau, kappa, n_e, mass_fractions,
    sigma, report (the final validate_grid() report), diagnostics
    (the monitor-function components from the final iteration).
    """
    r = np.asarray(r, dtype=float)
    v = np.asarray(v, dtype=float)

    boost = compression_boost
    last_report = None
    last_diag = None

    for iteration in range(1, max_iterations + 1):
        x, M, diagnostics = compute_monitor_function(
            r, v, T, rho, tau, mass_fractions, kappa,
            weights=weights,
            compression_policy=compression_policy,
            compression_threshold=compression_threshold,
            compression_boost=boost,
        )
        x_new = equidistribute(x, M, n_target)
        r_new = np.exp(x_new)

        (v_new, T_new, rho_new, tau_new, kappa_new, n_e_new, mf_new) = interpolate_fields(
            x, x_new, v, T, rho, tau, kappa, n_e, mass_fractions
        )

        ok, report = validate_grid(
            r_new, v_new, n_target,
            max_abs_sigma=max_abs_sigma,
            min_dr_r=min_dr_r,
            compression_threshold=compression_threshold,
            compression_policy=compression_policy,
        )
        last_report, last_diag = report, diagnostics

        if verbose:
            print(
                f"[build_target_grid] iteration {iteration}: "
                f"N={report['n_points']} max|sigma|={report['max_abs_sigma']:.3g} "
                f"min_dr/r={report['min_dr_r']:.3g} ok={ok}"
            )

        if ok:
            return {
                "r": r_new, "v": v_new, "T": T_new, "rho": rho_new,
                "tau": tau_new, "kappa": kappa_new, "n_e": n_e_new,
                "mass_fractions": mf_new,
                "sigma": sigma_velocity(r_new, v_new),
                "report": report,
                "diagnostics": diagnostics,
            }

        if compression_policy != "resolve" or report["sigma_ok"]:
            # Nothing left to try automatically -- e.g. a spacing or
            # point-count failure isn't fixed by boosting compression
            # weight. Stop iterating and fail below with full context.
            break
        boost *= boost_growth

    raise ValueError(
        "build_target_grid: could not produce a grid satisfying the "
        f"acceptance criteria within {max_iterations} iterations. "
        f"Last report: {last_report}. Consider raising n_target, "
        "raising max_abs_sigma if the physics genuinely requires it, "
        "or switching compression_policy."
    )
