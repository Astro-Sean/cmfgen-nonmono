"""Gap D: local interpolation-error metric.

Compares the remapped CMFGEN grid against the *original* STELLA structure
and reports the **maximum local** relative error per physical variable,
not just a global RMS.  This catches the case where a sharp composition
transition or velocity extremum is numerically washed out even though the
global RMS error is tiny.

For each target grid point we interpolate the original STELLA quantity at
that radius and compare:

    epsilon_q = |q_target - q_STELLA(r_target)| / |q_STELLA(r_target)|

The reported metric is ``max(epsilon_q)`` over all grid points, per
variable.  A grid is accepted only if every variable's max-local error
stays below its configured tolerance.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from grid_mapping.shell_model import SourceShellGrid, TargetGrid


@dataclass
class ErrorReport:
    """Per-variable max-local interpolation error."""
    max_dlogtau: float = 0.0
    max_dlogrho: float = 0.0
    max_dv_frac: float = 0.0
    max_dlogT: float = 0.0
    max_dX: float = 0.0
    worst_location_r: float = 0.0
    worst_variable: str = ""
    passed: bool = False
    per_variable: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        return {
            'max_dlogtau': self.max_dlogtau,
            'max_dlogrho': self.max_dlogrho,
            'max_dv_frac': self.max_dv_frac,
            'max_dlogT': self.max_dlogT,
            'max_dX': self.max_dX,
            'worst_location_r_cm': self.worst_location_r,
            'worst_variable': self.worst_variable,
            'passed': self.passed,
            'per_variable': self.per_variable,
        }


def _interp_from_source(r_target: np.ndarray, r_src: np.ndarray,
                        q_src: np.ndarray, log: bool = False) -> np.ndarray:
    """Interpolate a source quantity onto target radii.

    Source arrays are outer->inner (decreasing r).  ``np.interp`` needs
    increasing x, so we reverse.  If ``log``, interpolate in log10 of the
    quantity (with a floor) to handle many-order-of-magnitude ranges.
    """
    r_src_inc = r_src[::-1]
    if log:
        q = np.log10(np.maximum(q_src, 1e-300))
        out = np.interp(r_target, r_src_inc, q[::-1])
        return np.power(10.0, out)
    out = np.interp(r_target, r_src_inc, q_src[::-1])
    return out


def local_interpolation_error(target: TargetGrid, source: SourceShellGrid,
                              cfg) -> ErrorReport:
    """Compute the max-local interpolation error of ``target`` vs ``source``.

    The comparison is made at the target grid *point* radii: we ask "if we
    had kept the full STELLA resolution, what value would each quantity
    have at this radius, and how far is the remapped value from it?"
    """
    rep = ErrorReport()
    r_t = target.r_cm
    r_s = source.r_centre_cm
    n_s = source.n_shells
    if n_s < 2 or target.n_points < 2:
        return rep

    # Only measure error in the radiatively-important region (tau > floor).
    # The optically-thin outer region is covered by the skeleton; zero-mass
    # cells there are numerical artifacts, not physics errors.
    if len(target.tau_cumulative) == len(r_t):
        radiative = target.tau_cumulative >= cfg.tau_metric_min
    else:
        radiative = np.ones(len(r_t), dtype=bool)

    def _masked_max(arr):
        if not np.any(radiative):
            return 0.0
        return float(np.max(np.where(radiative, arr, 0.0)))

    # --- tau (log) ---
    if len(source.tau_cumulative) == n_s and len(target.tau_cumulative) == len(r_t):
        tau_src = _interp_from_source(r_t, r_s, source.tau_cumulative, log=True)
        eps_tau = np.abs(np.log10((target.tau_cumulative + cfg.tau_floor) /
                                  (tau_src + cfg.tau_floor)))
        rep.max_dlogtau = _masked_max(eps_tau)
        rep.per_variable['dlogtau'] = rep.max_dlogtau
    else:
        eps_tau = np.zeros(len(r_t))

    # --- density (log) ---
    rho_src = _interp_from_source(r_t, r_s, source.density_gcm3, log=True)
    eps_rho = np.abs(np.log10(np.maximum(target.density_gcm3, 1e-30) /
                              np.maximum(rho_src, 1e-30)))
    rep.max_dlogrho = _masked_max(eps_rho)
    rep.per_variable['dlogrho'] = rep.max_dlogrho

    # --- velocity (fractional) ---
    v_src = _interp_from_source(r_t, r_s, source.velocity_cms, log=False)
    eps_v = np.abs(target.velocity_cms - v_src) / np.maximum(np.abs(v_src), 1.0)
    rep.max_dv_frac = _masked_max(eps_v)
    rep.per_variable['dv_frac'] = rep.max_dv_frac

    # --- temperature (log) ---
    T_src = _interp_from_source(r_t, r_s, source.temperature_K, log=True)
    eps_T = np.abs(np.log10(np.maximum(target.temperature_K, 1.0) /
                            np.maximum(T_src, 1.0)))
    rep.max_dlogT = _masked_max(eps_T)
    rep.per_variable['dlogT'] = rep.max_dlogT

    # --- composition (max |dX| over species) ---
    eps_X = np.zeros(len(r_t))
    for elem, frac_t in target.mass_fractions.items():
        if elem not in source.mass_fractions:
            continue
        X_src = _interp_from_source(r_t, r_s, source.mass_fractions[elem],
                                    log=False)
        eps_X = np.maximum(eps_X, np.abs(frac_t - X_src))
    rep.max_dX = _masked_max(eps_X)
    rep.per_variable['dX'] = rep.max_dX

    # --- worst location / variable (within radiative region) ---
    all_eps = {
        'dlogtau': eps_tau,
        'dlogrho': eps_rho,
        'dv_frac': eps_v,
        'dlogT': eps_T,
        'dX': eps_X,
    }
    worst_name = max(rep.per_variable, key=rep.per_variable.get)
    worst_arr = all_eps.get(worst_name, np.zeros(len(r_t)))
    if worst_arr.size and np.any(radiative):
        masked_arr = np.where(radiative, worst_arr, -np.inf)
        rep.worst_variable = worst_name
        rep.worst_location_r = float(r_t[int(np.argmax(masked_arr))])

    # --- pass/fail against tolerances ---
    checks = {
        'dlogtau': (rep.max_dlogtau, cfg.error_tol_dlogtau),
        'dlogrho': (rep.max_dlogrho, cfg.error_tol_dlogrho),
        'dv_frac': (rep.max_dv_frac, cfg.error_tol_dv_frac),
        'dlogT': (rep.max_dlogT, cfg.error_tol_dlogT),
        'dX': (rep.max_dX, cfg.error_tol_dX),
    }
    rep.passed = all(val <= tol for val, tol in checks.values())
    return rep
