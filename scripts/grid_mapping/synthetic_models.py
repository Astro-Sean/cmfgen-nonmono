"""Synthetic test models for grid_mapping validation.

Generates synthetic SourceShellGrid objects with known properties for
testing the grid mapping pipeline:

- Monotonic velocity (homologous expansion)
- Single velocity reversal (CSM interaction)
- Multiple velocity reversals
- Shock (density jump)
- Composition interface
- Ionization front (Ne jump)

These models are used by the unit tests and can be used for end-to-end
solver validation via run_model_grid.py.
"""

from __future__ import annotations

import numpy as np

from .shell_model import (
    SourceShellGrid, compute_shell_volume, compute_shell_centres,
    density_from_mass_volume, MSUN,
)
from .opacity_tau import compute_tau


def _make_base_grid(n: int = 100,
                    r_outer: float = 10e10,
                    r_inner: float = 1e10) -> tuple:
    """Create base geometry: boundaries, volumes, centres."""
    r_bound = np.linspace(r_outer, r_inner, n + 1)  # decreasing
    vol = compute_shell_volume(r_bound)
    r_centre = compute_shell_centres(r_bound)
    return r_bound, vol, r_centre


def synthetic_monotonic(n: int = 100) -> SourceShellGrid:
    """Homologous expansion: V = R * t_exp^-1, sigma = 0 everywhere."""
    r_bound, vol, r_centre = _make_base_grid(n)

    # Homologous: V proportional to R
    t_exp = 30.0 * 86400.0  # 30 days in seconds
    v = r_centre / t_exp  # cm/s

    T = np.ones(n) * 5000.0
    rho = np.ones(n) * 1e-10
    mass = rho * vol
    kappa = np.ones(n) * 0.34
    mf = {'H': np.ones(n) * 0.7, 'He': np.ones(n) * 0.28, 'O': np.ones(n) * 0.02}

    grid = SourceShellGrid(
        r_bound_cm=r_bound, r_centre_cm=r_centre,
        shell_volume_cm3=vol, shell_mass_g=mass,
        density_gcm3=rho, velocity_cms=v, temperature_K=T,
        kappa_cgs=kappa, mass_fractions=mf,
        time_days=30.0, source_file='synthetic_monotonic',
        n_zones_original=n,
    )
    grid.tau_cumulative = compute_tau(r_bound, rho, kappa)
    return grid


def synthetic_single_reversal(n: int = 100) -> SourceShellGrid:
    """Single velocity reversal: V has a local minimum at mid-radius.

    Simulates CSM interaction where a reverse shock slows the outer
    ejecta, creating a velocity minimum.
    """
    r_bound, vol, r_centre = _make_base_grid(n)

    # Non-monotonic: V decreases to a minimum at 40% of radius, then increases
    t_exp = 30.0 * 86400.0
    v_homologous = r_centre / t_exp

    # Add a velocity minimum at index 40
    reversal_idx = 40
    v = v_homologous.copy()
    # Reduce velocity around the reversal
    for i in range(max(0, reversal_idx - 5), min(n, reversal_idx + 5)):
        dist = abs(i - reversal_idx)
        v[i] *= (1.0 - 0.5 * np.exp(-dist**2 / 4.0))

    T = np.ones(n) * 5000.0
    rho = np.ones(n) * 1e-10
    mass = rho * vol
    kappa = np.ones(n) * 0.34
    mf = {'H': np.ones(n) * 0.7, 'He': np.ones(n) * 0.28, 'O': np.ones(n) * 0.02}

    grid = SourceShellGrid(
        r_bound_cm=r_bound, r_centre_cm=r_centre,
        shell_volume_cm3=vol, shell_mass_g=mass,
        density_gcm3=rho, velocity_cms=v, temperature_K=T,
        kappa_cgs=kappa, mass_fractions=mf,
        time_days=30.0, source_file='synthetic_single_reversal',
        n_zones_original=n,
    )
    grid.tau_cumulative = compute_tau(r_bound, rho, kappa)
    return grid


def synthetic_multi_reversal(n: int = 100) -> SourceShellGrid:
    """Multiple velocity reversals: V has several local minima/maxima.

    Simulates multiple shell collisions or CSM interaction episodes.
    """
    r_bound, vol, r_centre = _make_base_grid(n)

    t_exp = 30.0 * 86400.0
    v = r_centre / t_exp

    # Add multiple reversals at different radii
    for reversal_idx, amplitude in [(25, 0.4), (50, 0.3), (75, 0.2)]:
        for i in range(max(0, reversal_idx - 4), min(n, reversal_idx + 4)):
            dist = abs(i - reversal_idx)
            v[i] *= (1.0 - amplitude * np.exp(-dist**2 / 3.0))

    T = np.ones(n) * 5000.0
    rho = np.ones(n) * 1e-10
    mass = rho * vol
    kappa = np.ones(n) * 0.34
    mf = {'H': np.ones(n) * 0.7, 'He': np.ones(n) * 0.28, 'O': np.ones(n) * 0.02}

    grid = SourceShellGrid(
        r_bound_cm=r_bound, r_centre_cm=r_centre,
        shell_volume_cm3=vol, shell_mass_g=mass,
        density_gcm3=rho, velocity_cms=v, temperature_K=T,
        kappa_cgs=kappa, mass_fractions=mf,
        time_days=30.0, source_file='synthetic_multi_reversal',
        n_zones_original=n,
    )
    grid.tau_cumulative = compute_tau(r_bound, rho, kappa)
    return grid


def synthetic_shock(n: int = 100) -> SourceShellGrid:
    """Strong shock: 10x density jump at mid-radius."""
    r_bound, vol, r_centre = _make_base_grid(n)

    t_exp = 30.0 * 86400.0
    v = r_centre / t_exp

    T = np.ones(n) * 5000.0
    rho = np.ones(n) * 1e-10

    # Shock at index 50: density jumps by 10x
    shock_idx = 50
    rho[shock_idx:] *= 10.0

    # Velocity also jumps at the shock
    v[shock_idx:] *= 0.5  # post-shock velocity is lower

    mass = rho * vol
    kappa = np.ones(n) * 0.34
    mf = {'H': np.ones(n) * 0.7, 'He': np.ones(n) * 0.28, 'O': np.ones(n) * 0.02}

    grid = SourceShellGrid(
        r_bound_cm=r_bound, r_centre_cm=r_centre,
        shell_volume_cm3=vol, shell_mass_g=mass,
        density_gcm3=rho, velocity_cms=v, temperature_K=T,
        kappa_cgs=kappa, mass_fractions=mf,
        time_days=30.0, source_file='synthetic_shock',
        n_zones_original=n,
    )
    grid.tau_cumulative = compute_tau(r_bound, rho, kappa)
    return grid


def synthetic_composition_interface(n: int = 100) -> SourceShellGrid:
    """Composition interface: H/He ratio changes sharply at mid-radius."""
    r_bound, vol, r_centre = _make_base_grid(n)

    t_exp = 30.0 * 86400.0
    v = r_centre / t_exp

    T = np.ones(n) * 5000.0
    rho = np.ones(n) * 1e-10
    mass = rho * vol
    kappa = np.ones(n) * 0.34

    # Composition interface at index 50
    # Outer: H-rich (H=0.7, He=0.28)
    # Inner: He-rich (H=0.1, He=0.8, O=0.1)
    iface_idx = 50
    H = np.ones(n) * 0.7
    He = np.ones(n) * 0.28
    O = np.ones(n) * 0.02
    H[iface_idx:] = 0.1
    He[iface_idx:] = 0.8
    O[iface_idx:] = 0.1
    mf = {'H': H, 'He': He, 'O': O}

    # Add Ni56 only in the inner region
    ni56 = np.zeros(n)
    ni56[iface_idx:] = 0.05

    grid = SourceShellGrid(
        r_bound_cm=r_bound, r_centre_cm=r_centre,
        shell_volume_cm3=vol, shell_mass_g=mass,
        density_gcm3=rho, velocity_cms=v, temperature_K=T,
        kappa_cgs=kappa, mass_fractions=mf, ni56_initial=ni56,
        time_days=30.0, source_file='synthetic_composition_interface',
        n_zones_original=n,
    )
    grid.tau_cumulative = compute_tau(r_bound, rho, kappa)
    return grid


def synthetic_ionization_front(n: int = 100) -> SourceShellGrid:
    """Ionization front: electron density jumps by 5 orders of magnitude."""
    r_bound, vol, r_centre = _make_base_grid(n)

    t_exp = 30.0 * 86400.0
    v = r_centre / t_exp

    T = np.ones(n) * 5000.0
    # Temperature also jumps at the front
    front_idx = 50
    T[front_idx:] = 15000.0  # hotter inner region (ionized)

    rho = np.ones(n) * 1e-10
    mass = rho * vol
    kappa = np.ones(n) * 0.34
    # Opacity also changes at the front
    kappa[front_idx:] = 0.1  # lower opacity in ionized region

    mf = {'H': np.ones(n) * 0.7, 'He': np.ones(n) * 0.28, 'O': np.ones(n) * 0.02}

    # Electron density: 5 orders of magnitude jump
    ne = np.ones(n) * 1e5  # neutral outer region
    ne[front_idx:] = 1e10  # ionized inner region
    nbar = np.ones(n) * 1e15

    grid = SourceShellGrid(
        r_bound_cm=r_bound, r_centre_cm=r_centre,
        shell_volume_cm3=vol, shell_mass_g=mass,
        density_gcm3=rho, velocity_cms=v, temperature_K=T,
        kappa_cgs=kappa, mass_fractions=mf,
        electron_density=ne, baryon_density=nbar,
        time_days=30.0, source_file='synthetic_ionization_front',
        n_zones_original=n,
    )
    grid.tau_cumulative = compute_tau(r_bound, rho, kappa)
    return grid


def synthetic_all_features(n: int = 100) -> SourceShellGrid:
    """Combined model with shock, reversal, composition interface, and ionization front."""
    r_bound, vol, r_centre = _make_base_grid(n)

    t_exp = 30.0 * 86400.0
    v = r_centre / t_exp

    # Velocity reversal at index 30
    for i in range(26, 36):
        dist = abs(i - 30)
        v[i] *= (1.0 - 0.4 * np.exp(-dist**2 / 4.0))

    T = np.ones(n) * 5000.0
    rho = np.ones(n) * 1e-10

    # Shock at index 60
    rho[60:] *= 10.0
    v[60:] *= 0.5

    # Ionization front at index 70
    T[70:] = 15000.0
    kappa = np.ones(n) * 0.34
    kappa[70:] = 0.1

    mass = rho * vol

    # Composition interface at index 50
    H = np.ones(n) * 0.7
    He = np.ones(n) * 0.28
    O = np.ones(n) * 0.02
    H[50:] = 0.1
    He[50:] = 0.8
    O[50:] = 0.1
    mf = {'H': H, 'He': He, 'O': O}

    ni56 = np.zeros(n)
    ni56[50:] = 0.05

    ne = np.ones(n) * 1e5
    ne[70:] = 1e10
    nbar = np.ones(n) * 1e15

    grid = SourceShellGrid(
        r_bound_cm=r_bound, r_centre_cm=r_centre,
        shell_volume_cm3=vol, shell_mass_g=mass,
        density_gcm3=rho, velocity_cms=v, temperature_K=T,
        kappa_cgs=kappa, mass_fractions=mf, ni56_initial=ni56,
        electron_density=ne, baryon_density=nbar,
        time_days=30.0, source_file='synthetic_all_features',
        n_zones_original=n,
    )
    grid.tau_cumulative = compute_tau(r_bound, rho, kappa)
    return grid


# ===========================================================================
# Complex composition models for cocoon / H-free shell testing
# ===========================================================================

def synthetic_hfree_cocoon(n: int = 200,
                           cocoon_frac: float = 0.2,
                           cocoon_center: float = 0.5,
                           optically_thick: bool = True) -> SourceShellGrid:
    """Detached H-free, metal-rich cocoon embedded in H-rich ejecta.

    Structure (outer -> inner):
        H-rich ejecta  |  H-free metal-rich cocoon  |  H-rich ejecta

    The cocoon has:
        X_H = 0, X_He = 0.3, X_O = 0.3, X_Si = 0.2, X_Ca = 0.05, X_Fe = 0.15
    The surrounding ejecta has:
        X_H = 0.7, X_He = 0.28, X_O = 0.02

    Parameters
    ----------
    n : int
        Number of source zones.
    cocoon_frac : float
        Fraction of the radial extent occupied by the cocoon.
    cocoon_center : float
        Radial location of the cocoon centre (fraction of total extent).
    optically_thick : bool
        If True, the cocoon has 10x higher density and opacity, making it
        optically thick.  If False, the cocoon has the same density as
        the surroundings (optically thin).
    """
    r_bound, vol, r_centre = _make_base_grid(n)
    t_exp = 30.0 * 86400.0
    v = r_centre / t_exp
    T = np.ones(n) * 5000.0
    rho = np.ones(n) * 1e-10
    kappa = np.ones(n) * 0.34

    # Cocoon boundaries (in terms of zone index)
    cocoon_half = int(n * cocoon_frac / 2)
    cocoon_mid = int(n * cocoon_center)
    c_outer = cocoon_mid - cocoon_half  # outer boundary index
    c_inner = cocoon_mid + cocoon_half  # inner boundary index

    # Default H-rich composition
    H = np.ones(n) * 0.70
    He = np.ones(n) * 0.28
    O = np.ones(n) * 0.02
    C = np.zeros(n)
    Ne = np.zeros(n)
    Mg = np.zeros(n)
    Si = np.zeros(n)
    S = np.zeros(n)
    Ca = np.zeros(n)
    Fe = np.zeros(n)
    Ni = np.zeros(n)

    # Cocoon: H-free, metal-rich
    H[c_outer:c_inner] = 0.0
    He[c_outer:c_inner] = 0.30
    O[c_outer:c_inner] = 0.30
    Si[c_outer:c_inner] = 0.20
    Ca[c_outer:c_inner] = 0.05
    Fe[c_outer:c_inner] = 0.15

    if optically_thick:
        rho[c_outer:c_inner] *= 10.0
        kappa[c_outer:c_inner] = 3.4  # high opacity
        T[c_outer:c_inner] = 10000.0

    mass = rho * vol
    mf = {'H': H, 'He': He, 'C': C, 'N': np.zeros(n), 'O': O,
          'Ne': Ne, 'Na': np.zeros(n), 'Mg': Mg, 'Al': np.zeros(n),
          'Si': Si, 'S': S, 'Ar': np.zeros(n), 'Ca': Ca, 'Fe': Fe,
          'Ni': Ni}

    grid = SourceShellGrid(
        r_bound_cm=r_bound, r_centre_cm=r_centre,
        shell_volume_cm3=vol, shell_mass_g=mass,
        density_gcm3=rho, velocity_cms=v, temperature_K=T,
        kappa_cgs=kappa, mass_fractions=mf, ni56_initial=Ni,
        time_days=30.0, source_file='synthetic_hfree_cocoon',
        n_zones_original=n,
    )
    grid.tau_cumulative = compute_tau(r_bound, rho, kappa)
    grid.cocoon_outer_idx = c_outer
    grid.cocoon_inner_idx = c_inner
    return grid


def synthetic_large_opacity_jump(n: int = 200) -> SourceShellGrid:
    """Large opacity jump (680x) across a single cell, mimicking the
    production pathological grid failure."""
    r_bound, vol, r_centre = _make_base_grid(n)
    t_exp = 30.0 * 86400.0
    v = r_centre / t_exp
    T = np.ones(n) * 5000.0
    rho = np.ones(n) * 1e-10
    kappa = np.ones(n) * 0.34

    # Create a sharp opacity/density jump at index 150
    jump_idx = 150
    rho[jump_idx:] *= 260.0
    kappa[jump_idx:] *= 2.5  # ~680x chi jump (260 * 2.5)
    T[jump_idx:] = 8000.0

    mass = rho * vol
    H = np.ones(n) * 0.7
    He = np.ones(n) * 0.28
    O = np.ones(n) * 0.02
    mf = {'H': H, 'He': He, 'C': np.zeros(n), 'N': np.zeros(n), 'O': O,
          'Ne': np.zeros(n), 'Na': np.zeros(n), 'Mg': np.zeros(n),
          'Al': np.zeros(n), 'Si': np.zeros(n), 'S': np.zeros(n),
          'Ar': np.zeros(n), 'Ca': np.zeros(n), 'Fe': np.zeros(n),
          'Ni': np.zeros(n)}

    grid = SourceShellGrid(
        r_bound_cm=r_bound, r_centre_cm=r_centre,
        shell_volume_cm3=vol, shell_mass_g=mass,
        density_gcm3=rho, velocity_cms=v, temperature_K=T,
        kappa_cgs=kappa, mass_fractions=mf, ni56_initial=np.zeros(n),
        time_days=30.0, source_file='synthetic_large_opacity_jump',
        n_zones_original=n,
    )
    grid.tau_cumulative = compute_tau(r_bound, rho, kappa)
    grid.jump_idx = jump_idx
    return grid


def synthetic_sharp_density_jump(n: int = 200) -> SourceShellGrid:
    """Sharp density jump (260x) with uniform composition."""
    r_bound, vol, r_centre = _make_base_grid(n)
    t_exp = 30.0 * 86400.0
    v = r_centre / t_exp
    T = np.ones(n) * 5000.0
    rho = np.ones(n) * 1e-10

    jump_idx = 100
    rho[jump_idx:] *= 260.0
    kappa = np.ones(n) * 0.34

    mass = rho * vol
    H = np.ones(n) * 0.7
    He = np.ones(n) * 0.28
    O = np.ones(n) * 0.02
    mf = {'H': H, 'He': He, 'C': np.zeros(n), 'N': np.zeros(n), 'O': O,
          'Ne': np.zeros(n), 'Na': np.zeros(n), 'Mg': np.zeros(n),
          'Al': np.zeros(n), 'Si': np.zeros(n), 'S': np.zeros(n),
          'Ar': np.zeros(n), 'Ca': np.zeros(n), 'Fe': np.zeros(n),
          'Ni': np.zeros(n)}

    grid = SourceShellGrid(
        r_bound_cm=r_bound, r_centre_cm=r_centre,
        shell_volume_cm3=vol, shell_mass_g=mass,
        density_gcm3=rho, velocity_cms=v, temperature_K=T,
        kappa_cgs=kappa, mass_fractions=mf, ni56_initial=np.zeros(n),
        time_days=30.0, source_file='synthetic_sharp_density_jump',
        n_zones_original=n,
    )
    grid.tau_cumulative = compute_tau(r_bound, rho, kappa)
    grid.jump_idx = jump_idx
    return grid


def synthetic_multi_interface(n: int = 200) -> SourceShellGrid:
    """Multiple composition interfaces: H-rich / He-rich / O-rich / Si-rich."""
    r_bound, vol, r_centre = _make_base_grid(n)
    t_exp = 30.0 * 86400.0
    v = r_centre / t_exp
    T = np.ones(n) * 5000.0
    rho = np.ones(n) * 1e-10
    kappa = np.ones(n) * 0.34

    # Four layers at indices 50, 100, 150
    H = np.ones(n) * 0.70
    He = np.ones(n) * 0.28
    O = np.ones(n) * 0.02
    Si = np.zeros(n)
    Ca = np.zeros(n)
    Fe = np.zeros(n)

    # Layer 1 (outer, 0-50): H-rich
    # Layer 2 (50-100): He-rich
    H[50:100] = 0.10; He[50:100] = 0.80; O[50:100] = 0.10
    # Layer 3 (100-150): O-rich
    H[100:150] = 0.0; He[100:150] = 0.20; O[100:150] = 0.60; Si[100:150] = 0.20
    # Layer 4 (inner, 150+): Si-rich (Fe-group)
    H[150:] = 0.0; He[150:] = 0.10; O[150:] = 0.10; Si[150:] = 0.40
    Ca[150:] = 0.05; Fe[150:] = 0.35

    mass = rho * vol
    mf = {'H': H, 'He': He, 'C': np.zeros(n), 'N': np.zeros(n), 'O': O,
          'Ne': np.zeros(n), 'Na': np.zeros(n), 'Mg': np.zeros(n),
          'Al': np.zeros(n), 'Si': Si, 'S': np.zeros(n),
          'Ar': np.zeros(n), 'Ca': Ca, 'Fe': Fe, 'Ni': np.zeros(n)}

    grid = SourceShellGrid(
        r_bound_cm=r_bound, r_centre_cm=r_centre,
        shell_volume_cm3=vol, shell_mass_g=mass,
        density_gcm3=rho, velocity_cms=v, temperature_K=T,
        kappa_cgs=kappa, mass_fractions=mf, ni56_initial=np.zeros(n),
        time_days=30.0, source_file='synthetic_multi_interface',
        n_zones_original=n,
    )
    grid.tau_cumulative = compute_tau(r_bound, rho, kappa)
    return grid


# Registry of all synthetic models
SYNTHETIC_MODELS = {
    'monotonic': synthetic_monotonic,
    'single_reversal': synthetic_single_reversal,
    'multi_reversal': synthetic_multi_reversal,
    'shock': synthetic_shock,
    'composition_interface': synthetic_composition_interface,
    'ionization_front': synthetic_ionization_front,
    'all_features': synthetic_all_features,
    'hfree_cocoon': synthetic_hfree_cocoon,
    'hfree_cocoon_thin': lambda n=200: synthetic_hfree_cocoon(n, optically_thick=False),
    'large_opacity_jump': synthetic_large_opacity_jump,
    'sharp_density_jump': synthetic_sharp_density_jump,
    'multi_interface': synthetic_multi_interface,
}
