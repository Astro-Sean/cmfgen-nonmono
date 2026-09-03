"""Opacity validation and optical-depth reintegration.

Provides routines to compute cumulative optical depth from density and
opacity, identify opacity sources, and find optical-depth landmarks
(tau = 0.01, 0.1, 1, 10, 100, photosphere at tau = 2/3).

Opacity sources are explicitly distinguished:
- Rosseland mean (from STELLA NLTE tables)
- Electron-scattering (constant kappa_es = 0.34 cm^2/g)
- Kramers analytical approximation
"""

from __future__ import annotations

from enum import Enum
from typing import Optional

import numpy as np

from .shell_model import AVOGADRO, K_BOLTZ, M_H


class OpacitySource(Enum):
    """Identifies which opacity source was used for tau computation."""
    STELLA_NLTE_ROSSELAND = "stella_nlte_rosseland"
    ELECTRON_SCATTERING = "electron_scattering"
    KRAMERS = "kramers"
    HYBRID = "hybrid"


# Electron scattering opacity (Thomson)
KAPPA_ES = 0.34  # cm^2/g (solar composition, fully ionized)


def compute_tau(r_bound_cm: np.ndarray, density: np.ndarray, kappa: np.ndarray) -> np.ndarray:
    """Compute cumulative optical depth from outside inward.

    Uses trapezoidal rule between shell centres:
        tau[i] = tau[i-1] + 0.5*(chi[i-1] + chi[i]) * |dr|
    where chi = kappa * rho.

    Parameters
    ----------
    r_bound_cm : (n+1,) shell boundaries, decreasing
    density : (n,) shell-centre densities
    kappa : (n,) shell-centre opacities (cm^2/g)

    Returns
    -------
    tau : (n,) cumulative optical depth at each shell centre, from outside
    """
    r = np.asarray(r_bound_cm, dtype=float)
    rho = np.asarray(density, dtype=float)
    k = np.asarray(kappa, dtype=float)
    n = len(rho)

    if n == 0:
        return np.array([])

    # Shell centre radii (for trapezoidal integration)
    # Use boundaries for dr
    if len(r) == n + 1:
        dr = np.abs(r[:-1] - r[1:])  # (n,) shell widths
    elif len(r) == n:
        # r is already shell centres
        dr = np.abs(np.diff(r))
        dr = np.append(dr, dr[-1])  # approximate last shell width
    else:
        raise ValueError(f"r has {len(r)} elements, expected {n} or {n+1}")

    chi = k * rho  # absorption coefficient (cm^-1)
    tau = np.zeros(n)
    for i in range(1, n):
        tau[i] = tau[i - 1] + 0.5 * (chi[i - 1] + chi[i]) * dr[i]
    return tau


def compute_tau_from_centres(r_cm: np.ndarray, density: np.ndarray, kappa: np.ndarray) -> np.ndarray:
    """Compute cumulative tau when r is shell centres (not boundaries).

    Uses trapezoidal rule between adjacent centres.
    """
    r = np.asarray(r_cm, dtype=float)
    rho = np.asarray(density, dtype=float)
    k = np.asarray(kappa, dtype=float)
    n = len(r)
    if n == 0:
        return np.array([])
    tau = np.zeros(n)
    for i in range(1, n):
        dr = abs(r[i - 1] - r[i])
        chi_avg = 0.5 * (k[i - 1] * rho[i - 1] + k[i] * rho[i])
        tau[i] = tau[i - 1] + chi_avg * dr
    return tau


def electron_scattering_opacity(mass_fractions: dict[str, np.ndarray]) -> np.ndarray:
    """Compute electron-scattering opacity kappa_es.

    kappa_es = sigma_T * n_e / rho = sigma_T * N_A * sum(X_i * Z_i / A_i)

    For fully ionized solar composition: ~0.34 cm^2/g.
    """
    elem_A = {'H': 1.0, 'He': 4.0, 'C': 12.0, 'N': 14.0, 'O': 16.0,
              'Ne': 20.0, 'Na': 23.0, 'Mg': 24.0, 'Al': 27.0, 'Si': 28.0,
              'S': 32.0, 'Ar': 40.0, 'Ca': 40.0, 'Fe': 56.0, 'Ni': 56.0,
              'FePeak': 56.0, 'Ni56': 56.0, 'Ni58': 58.0}
    elem_Z = {'H': 1.0, 'He': 2.0, 'C': 6.0, 'N': 7.0, 'O': 8.0,
              'Ne': 10.0, 'Na': 11.0, 'Mg': 12.0, 'Al': 13.0, 'Si': 14.0,
              'S': 16.0, 'Ar': 18.0, 'Ca': 20.0, 'Fe': 26.0, 'Ni': 28.0,
              'FePeak': 26.0, 'Ni56': 28.0, 'Ni58': 28.0}

    sigma_T = 6.6524e-25  # Thomson cross-section (cm^2)
    # kappa_es = sigma_T * N_A * sum(X_i * Z_i / A_i)
    n = len(next(iter(mass_fractions.values())))
    z_over_a = np.zeros(n)
    for elem, frac in mass_fractions.items():
        A = elem_A.get(elem, 28.0)
        Z = elem_Z.get(elem, 14.0)
        z_over_a += frac * Z / A
    return sigma_T * AVOGADRO * z_over_a


def kramers_opacity(T_K: np.ndarray, rho: np.ndarray,
                    X: np.ndarray, Y: np.ndarray) -> np.ndarray:
    """Kramers analytical opacity (free-free + bound-free).

    kappa_ff = 3.7e22 * rho * (1+X) * T^(-3.5)  (cm^2/g)
    kappa_bf = 4.4e22 * rho * Z * T^(-3.5)      (cm^2/g)

    For SN ejecta, we use a simplified form.
    """
    T = np.maximum(np.asarray(T_K, dtype=float), 3000.0)
    rho_arr = np.asarray(rho, dtype=float)
    X_arr = np.asarray(X, dtype=float)
    Y_arr = np.asarray(Y, dtype=float)
    kappa_ff = 3.7e22 * rho_arr * (1.0 + X_arr) * T ** (-3.5)
    # Bound-free is small for SN ejecta; include with Z ~ 0.02
    Z_arr = np.maximum(1.0 - X_arr - Y_arr, 0.0)
    kappa_bf = 4.4e22 * rho_arr * Z_arr * T ** (-3.5)
    kappa = kappa_ff + kappa_bf
    # Floor at electron scattering
    return np.maximum(kappa, KAPPA_ES * 0.1)


def find_tau_landmarks(tau: np.ndarray, r_cm: np.ndarray,
                       landmarks: Optional[list[float]] = None) -> dict[float, dict]:
    """Find radius and index at specific optical depth values.

    Parameters
    ----------
    tau : cumulative optical depth (outer->inner, increasing)
    r_cm : radius array (outer->inner, decreasing)
    landmarks : tau values to find (default: [2/3, 0.01, 0.1, 1, 10, 100])

    Returns
    -------
    dict mapping tau_value -> {'index': int, 'radius_cm': float} or {} if not found
    """
    if landmarks is None:
        landmarks = [2.0 / 3.0, 0.01, 0.1, 1.0, 10.0, 100.0]

    results = {}
    tau_arr = np.asarray(tau, dtype=float)
    r_arr = np.asarray(r_cm, dtype=float)

    for target_tau in landmarks:
        # Find first index where tau >= target_tau (from outside)
        idx = np.searchsorted(tau_arr, target_tau)
        if idx < len(tau_arr):
            # Interpolate radius
            if idx == 0:
                r_val = r_arr[0]
            else:
                # Linear interpolation in tau
                t0, t1 = tau_arr[idx - 1], tau_arr[idx]
                r0, r1 = r_arr[idx - 1], r_arr[idx]
                if t1 > t0:
                    frac = (target_tau - t0) / (t1 - t0)
                    r_val = r0 + frac * (r1 - r0)
                else:
                    r_val = r_arr[idx]
            results[target_tau] = {'index': int(idx), 'radius_cm': float(r_val)}
    return results


def find_photosphere(tau: np.ndarray, r_cm: np.ndarray) -> Optional[dict]:
    """Find the photosphere at tau = 2/3."""
    landmarks = find_tau_landmarks(tau, r_cm, landmarks=[2.0 / 3.0])
    return landmarks.get(2.0 / 3.0)


def max_delta_tau(tau: np.ndarray) -> float:
    """Maximum per-cell optical depth change."""
    if len(tau) < 2:
        return 0.0
    return float(np.max(np.abs(np.diff(tau))))


def points_per_tau_decade(tau: np.ndarray) -> dict:
    """Count grid points per optical depth decade."""
    tau_arr = np.asarray(tau, dtype=float)
    tau_arr = tau_arr[tau_arr > 0]
    if len(tau_arr) < 2:
        return {}
    log_tau = np.log10(tau_arr)
    decades = {}
    for d in range(int(np.floor(log_tau.min())), int(np.ceil(log_tau.max()))):
        count = np.sum((log_tau >= d) & (log_tau < d + 1))
        if count > 0:
            decades[d] = int(count)
    return decades
