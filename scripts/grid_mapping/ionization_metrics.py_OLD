"""Ionization metrics: ion fractions, Ne, Ye, opacity, and front diagnostics.

Provides Saha-equation electron density calculation (fallback when STELLA
Ne is unavailable) and ionization front diagnostics.
"""

from __future__ import annotations

from typing import Optional

import numpy as np

from .shell_model import AVOGADRO, K_BOLTZ, M_H


# Element properties for Saha calculation
ELEM_IONIZATION_EV = {
    'H': 13.6, 'He': 24.6, 'He2': 54.4, 'C': 11.3, 'N': 14.5,
    'O': 13.6, 'Ne': 21.6, 'Na': 5.1, 'Mg': 7.6, 'Al': 5.99,
    'Si': 8.15, 'S': 10.36, 'Ar': 15.76, 'Ca': 6.11, 'Fe': 7.87,
    'Ni': 7.64, 'FePeak': 7.87, 'Ni56': 7.64, 'Ni58': 7.64,
}
ELEM_A = {
    'H': 1.0, 'He': 4.0, 'C': 12.0, 'N': 14.0, 'O': 16.0,
    'Ne': 20.0, 'Na': 23.0, 'Mg': 24.0, 'Al': 27.0, 'Si': 28.0,
    'S': 32.0, 'Ar': 40.0, 'Ca': 40.0, 'Fe': 56.0, 'Ni': 56.0,
    'FePeak': 56.0, 'Ni56': 56.0, 'Ni58': 58.0,
}
ELEM_Z = {
    'H': 1, 'He': 2, 'C': 6, 'N': 7, 'O': 8,
    'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14,
    'S': 16, 'Ar': 18, 'Ca': 20, 'Fe': 26, 'Ni': 28,
    'FePeak': 26, 'Ni56': 28, 'Ni58': 28,
}

EV_TO_ERG = 1.60218e-12  # eV to erg


def saha_electron_density(rho: np.ndarray, T: np.ndarray,
                          mass_fractions: dict[str, np.ndarray],
                          n_e_ref: float = 1e8) -> np.ndarray:
    """Compute electron density using Saha ionization equilibrium.

    Single-pass Saha with fixed reference n_e_ref. This is a fallback
    when STELLA Ne is unavailable.

    Parameters
    ----------
    rho : density (g/cm^3)
    T : temperature (K)
    mass_fractions : dict of element -> mass fraction array
    n_e_ref : reference electron density for Saha ratio (cm^-3)
    """
    n = len(rho)
    T_arr = np.maximum(np.asarray(T, dtype=float), 3000.0)
    rho_arr = np.asarray(rho, dtype=float)

    # Compute atom density per element
    ne = np.zeros(n)
    for elem, frac in mass_fractions.items():
        if elem not in ELEM_IONIZATION_EV:
            continue
        A = ELEM_A.get(elem, 28.0)
        Z = ELEM_Z.get(elem, 14)
        chi = ELEM_IONIZATION_EV[elem] * EV_TO_ERG

        # Atom number density for this element
        n_elem = rho_arr * frac * AVOGADRO / A

        # Saha ratio: n_{i+1}/n_i = (2*pi*m_e*kT/h^2)^1.5 * 2/n_e * exp(-chi/kT)
        # Simplified: use fixed n_e_ref, compute ionization fraction
        kT = K_BOLTZ * T_arr
        saha_prefactor = 2.0 * (2.0 * np.pi * 9.109e-28 * kT / (6.626e-27 ** 2)) ** 1.5
        ratio = saha_prefactor / n_e_ref * np.exp(-chi / kT)

        # Ionization fraction (simplified: x^2/(1-x) = ratio -> x = sqrt(ratio/(1+ratio)))
        x = np.sqrt(ratio / (1.0 + ratio))
        x = np.clip(x, 0.0, 1.0)

        # For He, allow double ionization (simplified)
        if elem == 'He':
            chi2 = ELEM_IONIZATION_EV['He2'] * EV_TO_ERG
            ratio2 = saha_prefactor / n_e_ref * np.exp(-chi2 / kT)
            x2 = np.sqrt(ratio2 / (1.0 + ratio2))
            x2 = np.clip(x2, 0.0, 1.0)
            # He contributes 1 electron for single + 2 for double
            ne += n_elem * (x * (1 - x2) + 2 * x * x2)
        else:
            # Single ionization releases exactly ONE electron per ionized
            # atom (Saha n_{i+1}/n_i gives the singly-ionized fraction).
            # The previous code multiplied by Z (the full atomic number),
            # overestimating ne by a factor of ~Z (6-28x).  Do NOT multiply
            # by Z here.
            ne += n_elem * x

    return np.maximum(ne, 1.0)  # floor at 1 cm^-3


def compute_electron_fraction(ne: np.ndarray, n_baryon: np.ndarray) -> np.ndarray:
    """Compute Ye = Ne / n_baryon."""
    return ne / np.maximum(n_baryon, 1e-30)


def compute_ionization_diagnostics(grid) -> dict:
    """Compute ionization diagnostics for a grid (SourceShellGrid or TargetGrid).

    Returns dict with:
    - ne: electron density array
    - ye: electron fraction array (if baryon density available)
    - ne_range: (min, max) of Ne
    - ye_range: (min, max) of Ye
    - ionization_front_indices: where |dlog10(Ne)| > threshold
    """
    ne = getattr(grid, 'electron_density', None)
    n_baryon = getattr(grid, 'baryon_density', None)

    diagnostics = {
        'ne_available': ne is not None,
        'ye_available': False,
        'ne_range': None,
        'ye_range': None,
        'ionization_front_indices': [],
    }

    if ne is not None:
        diagnostics['ne_range'] = (float(ne.min()), float(ne.max()))

        if n_baryon is not None:
            ye = compute_electron_fraction(ne, n_baryon)
            diagnostics['ye_available'] = True
            diagnostics['ye_range'] = (float(ye.min()), float(ye.max()))

        # Detect ionization fronts
        log_ne = np.log10(np.maximum(ne, 1e-30))
        dlog_ne = np.abs(np.diff(log_ne))
        front_indices = np.where(dlog_ne > 1.0)[0].tolist()
        diagnostics['ionization_front_indices'] = front_indices

    return diagnostics
