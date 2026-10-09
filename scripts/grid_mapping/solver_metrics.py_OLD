"""Exact local CMFGEN velocity-gradient diagnostics.

CMFGEN computes the solver SIGMA from the R and V grid, **not** from the
SIGMA column written to SN_HYDRO_DATA.  The exact formula depends on
``VEL_LAW`` (set per-YAML via ``cmfgen.vel_law``, default 12):

* ``VEL_LAW=12`` -> ``RV_SN_MODEL_SNIIN`` reads R/V from ``input_hydro.dat``
  (written by ``prepare_VADAT.py`` from SN_HYDRO_DATA) and computes SIGMA
  with a distance-weighted central finite difference
  (``rv_sn_model_sniin.f`` L262-269).  ``sigma_velocity()`` reproduces this
  formula exactly.  This is the path the package was originally designed
  for and the project default.

* ``VEL_LAW=11`` -> ``SET_RV_HYDRO_MODEL_V3`` reads R/V/SIGMA from
  SN_HYDRO_DATA, builds its own R grid, fits a **monotone cubic** (PCHIP)
  interpolant to V(R) and takes the **analytic derivative** of that cubic
  (``set_rv_hydro_model_v3.f`` L419-433).  ``sigma_cmfgen_v11()``
  reproduces this.  The SIGMA column read from the file is used **only**
  for the density time-evolution correction (L232), NOT by the solver.

IMPORTANT: in *neither* path does the SIGMA column written to
SN_HYDRO_DATA directly feed the radiative-transfer solver.  Reducing
negative sigma therefore requires fixing the V(R) profile and the R
grid, not editing the sigma column.  Use ``sigma_solver(r, v, vel_law)``
to obtain the sigma the solver will actually compute for a given path.

VEL_LAW=12 finite-difference formula (rv_sn_model_sniin.f L262-269):
    T1 = (V[i-1]-V[i])/(R[i-1]-R[i])     ! left dV/dR
    T2 = (V[i]-V[i+1])/(R[i]-R[i+1])     ! right dV/dR
    dVdR = (T1*h_right + T2*h_left) / (h_left + h_right)  ! distance-weighted
    SIGMA = R*dVdR/V - 1                  ! = dlnV/dlnR - 1

Boundaries use the 2-point one-sided difference.
"""

from __future__ import annotations

import numpy as np

# SIGMA cap applied by CMFGEN (rv_sn_model_sniin.f)
SIGMA_CAP = -0.999
# The -1.2 < SIGMA < -1.0 -> -0.999 clamp range (from AGENTS.md)
SIGMA_CLAMP_LOW = -1.2
SIGMA_CLAMP_HIGH = -1.0


def sigma_velocity(r: np.ndarray, v: np.ndarray, preserve_sign: bool = True) -> np.ndarray:
    """Compute sigma = dlnV/dlnR - 1 using the exact CMFGEN formula.

    Matches ``rv_sn_model_sniin.f`` L262-269 exactly:
    distance-weighted average of left and right linear dV/dR, then
    R*dVdR/V - 1.

    Parameters
    ----------
    r : array, outer->inner (decreasing)
    v : array, velocity in any consistent units (km/s or cm/s)
    preserve_sign : bool
        If True, preserve the sign of V (for infall support).
        If False, use abs(V) in the denominator (legacy behaviour).

    Returns
    -------
    sigma : array, same length as r
    """
    n = len(r)
    sigma = np.zeros(n)
    eps = 1e-30
    # Minimum dR/R ratio to avoid degenerate grid points
    min_dr_r = 1e-10

    v_denom = v if preserve_sign else np.abs(v)

    for i in range(1, n - 1):
        h_left = r[i - 1] - r[i]
        h_right = r[i] - r[i + 1]
        # Check for degenerate grid spacing (dR/R too small)
        if abs(h_left) < eps or abs(h_right) < eps:
            sigma[i] = 0.0
            continue
        if abs(h_left) / r[i] < min_dr_r or abs(h_right) / r[i] < min_dr_r:
            # Degenerate: use one-sided difference from the non-degenerate side
            if abs(h_left) / r[i] >= min_dr_r:
                v_val = v_denom[i]
                if abs(v_val) > eps:
                    sigma[i] = r[i] * (v[i - 1] - v[i]) / h_left / v_val - 1.0
            elif abs(h_right) / r[i] >= min_dr_r:
                v_val = v_denom[i]
                if abs(v_val) > eps:
                    sigma[i] = r[i] * (v[i] - v[i + 1]) / h_right / v_val - 1.0
            continue
        T1 = (v[i - 1] - v[i]) / h_left
        T2 = (v[i] - v[i + 1]) / h_right
        dVdR = (T1 * h_right + T2 * h_left) / (h_left + h_right)
        v_val = v_denom[i]
        if abs(v_val) < eps:
            sigma[i] = 0.0
        else:
            sigma[i] = r[i] * dVdR / v_val - 1.0

    # Boundaries: CMFGEN 2-point one-sided
    if n >= 2:
        h = r[0] - r[1]
        if abs(h) > eps:
            v_val = v_denom[0]
            sigma[0] = r[0] * (v[0] - v[1]) / h / max(abs(v_val), eps) - 1.0 if abs(v_val) > eps else 0.0
        h = r[-2] - r[-1]
        if abs(h) > eps:
            v_val = v_denom[-1]
            sigma[-1] = r[-1] * (v[-2] - v[-1]) / h / max(abs(v_val), eps) - 1.0 if abs(v_val) > eps else 0.0

    return sigma


def sigma_cmfgen_v11(r: np.ndarray, v: np.ndarray,
                     preserve_sign: bool = True) -> np.ndarray:
    """Compute solver SIGMA for VEL_LAW=11 (SET_RV_HYDRO_MODEL_V3).

    Reproduces ``set_rv_hydro_model_v3.f`` L419-433: CMFGEN fits a monotone
    cubic (Fritsch-Carlson / PCHIP) interpolant to V(R) and evaluates the
    **analytic derivative** of that cubic at each grid point, then
    ``SIGMA = R * dV/dR / V - 1``.

    This is the sigma the VEL_LAW=11 solver actually uses.  It differs from
    the VEL_LAW=12 finite-difference ``sigma_velocity`` because the cubic
    derivative is not the same as the distance-weighted central difference.

    Parameters
    ----------
    r : array, outer->inner (decreasing)
    v : array, velocity in any consistent units (km/s or cm/s)
    preserve_sign : bool
        If True preserve the sign of V in the denominator (infall support).
    """
    from scipy.interpolate import PchipInterpolator

    r = np.asarray(r, dtype=float)
    v = np.asarray(v, dtype=float)
    n = len(r)
    sigma = np.zeros(n)
    if n < 2:
        return sigma

    # PCHIP requires increasing x.  r is decreasing (outer->inner), so
    # interpolate on the reversed (increasing) grid and map back.
    order = np.argsort(r)  # increasing
    r_inc = r[order]
    v_inc = v[order]

    # Drop duplicate radii (PCHIP needs strictly increasing x).
    keep = np.concatenate(([True], np.diff(r_inc) > 0))
    r_inc = r_inc[keep]
    v_inc = v_inc[keep]
    if len(r_inc) < 2:
        return sigma

    pchip = PchipInterpolator(r_inc, v_inc)
    # Analytic derivative of the cubic at the original (decreasing) nodes.
    dvdr = pchip(r, 1)  # first derivative evaluated at original r
    v_denom = v if preserve_sign else np.abs(v)
    eps = 1e-30
    mask = np.abs(v_denom) > eps
    sigma[mask] = r[mask] * dvdr[mask] / v_denom[mask] - 1.0
    # CMFGEN zeroes sigma where |sigma| < 1e-4 (set_rv_hydro_model_v3 L437).
    if np.max(np.abs(sigma)) < 1.0e-4:
        sigma[:] = 0.0
    return sigma


def sigma_solver(r: np.ndarray, v: np.ndarray, vel_law: int = 12,
                 preserve_sign: bool = True) -> np.ndarray:
    """Return the SIGMA the CMFGEN solver will compute for a given VEL_LAW.

    Parameters
    ----------
    r, v : arrays, outer->inner (decreasing)
    vel_law : 11 or 12 (project default 12).  Selects the formula matching
        the CMFGEN routine actually called for that VEL_LAW.
    preserve_sign : preserve V sign in the denominator.
    """
    if vel_law == 11:
        return sigma_cmfgen_v11(r, v, preserve_sign=preserve_sign)
    if vel_law == 12:
        return sigma_velocity(r, v, preserve_sign=preserve_sign)
    raise ValueError(f"Unsupported VEL_LAW={vel_law}; expected 11 or 12")


def apply_sigma_cap(sigma: np.ndarray) -> np.ndarray:
    """Apply the CMFGEN SIGMA cap: -1.2 < SIGMA < -1.0 -> -0.999.

    This matches the behaviour of rv_sn_model_sniin.f, adjust_r_grid_v4.f,
    and do_ng_accel_band_v4.f.  See AGENTS.md: do NOT remove this clamp.
    """
    s = sigma.copy()
    mask = (s > SIGMA_CLAMP_LOW) & (s < SIGMA_CLAMP_HIGH)
    s[mask] = SIGMA_CAP
    return s


def sigma_smoothness(sigma: np.ndarray) -> np.ndarray:
    """Compute |d(sigma)| between adjacent points.

    Large values indicate SIGMA jumps that may make the dGAMdR
    computation in fg_j_cmf_v13.f inaccurate.
    """
    return np.abs(np.diff(sigma))


def count_velocity_reversals(sigma: np.ndarray) -> int:
    """Count the number of velocity reversals (SIGMA < -1 regions).

    A reversal is a region where dV/dr changes sign, which corresponds
    to SIGMA < -1 (since SIGMA = dlnV/dlnR - 1, and dlnV/dlnR < 0
    means V decreases outward).
    """
    below = sigma < -1.0
    # Count transitions into below regions
    if not np.any(below):
        return 0
    transitions = np.diff(below.astype(int))
    return int(np.sum(transitions == 1)) + (1 if below[0] else 0)


def find_velocity_extrema(r: np.ndarray, v: np.ndarray) -> dict:
    """Find velocity maxima and minima (local extrema of V).

    Returns dict with 'maxima' and 'minima' lists of indices.
    """
    n = len(v)
    maxima = []
    minima = []
    for i in range(1, n - 1):
        if v[i] > v[i - 1] and v[i] > v[i + 1]:
            maxima.append(i)
        elif v[i] < v[i - 1] and v[i] < v[i + 1]:
            minima.append(i)
    return {'maxima': maxima, 'minima': minima}


def find_dvdr_sign_changes(r: np.ndarray, v: np.ndarray) -> list[int]:
    """Find indices where dV/dr changes sign (velocity gradient reversals).

    Returns list of indices *between* zones where the sign change occurs.
    Index i means the sign change is between shell i-1 and shell i.
    """
    n = len(v)
    if n < 3:
        return []
    # dV/dR at zone centres (outer->inner, so dR < 0)
    dvdr = np.diff(v) / np.diff(r)
    sign_changes = []
    for i in range(1, len(dvdr)):
        if dvdr[i - 1] * dvdr[i] < 0:  # sign change
            sign_changes.append(i)
    return sign_changes


def vol_exp_fac(sigma: np.ndarray, v_cms: np.ndarray, r_cm: np.ndarray,
                 t1_factor: float = 3.33564e-6) -> np.ndarray:
    """Compute the CMFGEN volume expansion factor.

    From rd_sn_data.f L255-257:
        VOL_EXP_FAC = 1/(1 - T1*V/R*(SIGMA+1)) / (1 - T1*V/R)^2

    where T1 = 3.33564E-6 (the CMF Doppler factor constant).

    This is a diagnostic only — the solver computes it internally.
    Note: this uses SIGMA directly (no Dessart guard), so it may
    diverge for SIGMA approaching -1.
    """
    eps = 1e-30
    v_safe = np.where(np.abs(v_cms) < eps, eps, v_cms)
    r_safe = np.where(np.abs(r_cm) < eps, eps, r_cm)
    t1_v_r = t1_factor * v_safe / r_safe
    denom1 = 1.0 - t1_v_r * (sigma + 1.0)
    denom2 = (1.0 - t1_v_r) ** 2
    denom1 = np.where(np.abs(denom1) < eps, eps, denom1)
    denom2 = np.where(np.abs(denom2) < eps, eps, denom2)
    return 1.0 / denom1 / denom2


def dessart_guard(sigma: np.ndarray, mu_sq: np.ndarray) -> np.ndarray:
    """Compute the Dessart+2015 regularized denominator.

    T1 = 1 + SIGMA*MU^2
    T1 = MAX(SQRT(T1*T1 + 0.01), -0.2*(1+SIGMA), 0.1)

    This replaces ABS(1+SIGMA*MU^2) in the stock Sobolev routines.
    Smooth for all SIGMA (no singularity at SIGMA=-1).
    """
    t1 = 1.0 + sigma * mu_sq
    return np.maximum(np.sqrt(t1 * t1 + 0.01), np.maximum(-0.2 * (1.0 + sigma), 0.1))
