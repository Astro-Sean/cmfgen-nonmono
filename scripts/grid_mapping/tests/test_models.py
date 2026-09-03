"""Two-event supernova pathological test model for grid-generation validation.

Creates a synthetic STELLA-like model containing all seven structures
required by the optical-depth-first grid audit:

1. Smooth ejecta (inner homologous core)
2. Dense interaction shell (density jump > 10x)
3. Sharp velocity transition (velocity reversal)
4. Sharp density jump (composition-driven)
5. H-rich region (outer shell, X_H ~ 0.74)
6. H-poor metal-rich region (inner ejecta, X_O ~ 0.96, X_Si ~ 0.57, X_Fe ~ 0.42)
7. Extended optically thin region (outer wind)

The model is designed to be physically pathological: the metal-rich
component is optically thin but spectrally important, and the
interaction shell has a large optical-depth contribution in a narrow
radial range.
"""

from __future__ import annotations

import numpy as np
from dataclasses import dataclass, field


@dataclass
class SyntheticSnapshot:
    """STELLA-like snapshot data (outer→inner, decreasing r)."""
    r: np.ndarray              # cm, outer→inner
    v: np.ndarray              # km/s, outer→inner
    T: np.ndarray              # K
    rho: np.ndarray            # g/cm³
    kappa: np.ndarray          # cm²/g (Rosseland mean)
    mass_fractions: dict       # element → array
    time_days: float = 35.0


def two_event_pathological(n_zones: int = 500) -> SyntheticSnapshot:
    """Create a two-event supernova model with all 7 pathological structures.

    Layout (outer→inner):
      - Zone 0..99:    Extended optically thin wind (H-rich, low rho)
      - Zone 100..199: H-rich outer shell (X_H=0.74, X_He=0.25)
      - Zone 200..249: Dense interaction shell (density jump 10x, velocity reversal)
      - Zone 250..349: H-poor metal-rich ejecta (X_O=0.96, X_Si=0.57, X_Fe=0.42)
      - Zone 350..499: Smooth inner core (homologous, high tau)

    Parameters
    ----------
    n_zones : int
        Number of zones (default 500, matching STELLA).

    Returns
    -------
    SyntheticSnapshot
    """
    # Layout fractions
    f_wind = 0.20       # 20% optically thin wind
    f_h_shell = 0.20    # 20% H-rich shell
    f_int = 0.10        # 10% interaction shell (narrow, dense)
    f_metal = 0.20      # 20% metal-rich ejecta
    f_core = 0.30       # 30% smooth inner core

    n_wind = int(n_zones * f_wind)
    n_h_shell = int(n_zones * f_h_shell)
    n_int = int(n_zones * f_int)
    n_metal = int(n_zones * f_metal)
    n_core = n_zones - n_wind - n_h_shell - n_int - n_metal

    # Radius grid (outer→inner, decreasing)
    # Outer boundary: 5e16 cm (large, optically thin)
    # Inner boundary: 1e14 cm (deep, high tau)
    r_outer = 5e16
    r_inner = 1e14
    # Log-uniform radius
    log_r = np.linspace(np.log10(r_outer), np.log10(r_inner), n_zones)
    r = 10.0 ** log_r

    # Velocity profile: mostly homologous (v ∝ r) but with a reversal
    # at the interaction shell and a jump at the metal-rich boundary
    v_homologous = 5000.0 * (r / r_outer)  # homologous, v_max = 5000 km/s at outer

    # Interaction shell: velocity drops to ~500 km/s then jumps back up
    i_int_start = n_wind + n_h_shell
    i_int_end = i_int_start + n_int
    v = v_homologous.copy()
    # Velocity reversal in interaction shell
    for i in range(i_int_start, i_int_end):
        frac = (i - i_int_start) / max(n_int - 1, 1)
        # Drop to 500 km/s then back up
        v[i] = v_homologous[i] * (1.0 - 0.9 * np.sin(np.pi * frac))

    # Sharp velocity jump at metal-rich boundary
    i_metal_start = i_int_end
    v[i_metal_start:i_metal_start + 5] *= 1.5  # jump up

    # Density profile: power law with jumps
    # Base: rho ∝ r^-8 (steep, typical SN), high enough for tau > 10 at inner
    rho_base = 1e-10 * (r / r_inner) ** (-8.0)

    # Wind: much lower density
    rho = rho_base.copy()
    rho[:n_wind] *= 1e-4

    # Interaction shell: 10x density jump
    rho[i_int_start:i_int_end] *= 10.0

    # Metal-rich region: moderate density
    rho[i_metal_start:i_metal_start + n_metal] *= 0.5

    # Temperature: adiabatic-like, with a hot interaction shell
    T_base = 5000.0 * (r / r_outer) ** (-0.5)
    T = T_base.copy()
    # Hot shell
    T[i_int_start:i_int_end] *= 3.0  # 15000 K in the shell
    # Cold metal-rich region
    T[i_metal_start:i_metal_start + n_metal] *= 0.6  # 3000 K

    # Opacity: hybrid (es + Kramers)
    # Wind: low opacity (cold, neutral)
    # Interaction shell: high opacity (hot, ionized)
    # Metal-rich: line opacity (high)
    kappa = np.full(n_zones, 0.34)  # electron scattering baseline
    kappa[:n_wind] = 0.05           # cold wind, low opacity
    kappa[i_int_start:i_int_end] = 1.0  # hot shell, high opacity
    kappa[i_metal_start:i_metal_start + n_metal] = 0.8  # metal lines

    # Composition
    mf = {
        'H': np.zeros(n_zones),
        'He': np.zeros(n_zones),
        'O': np.zeros(n_zones),
        'Si': np.zeros(n_zones),
        'Ca': np.zeros(n_zones),
        'Fe': np.zeros(n_zones),
    }

    # Wind + H-rich shell: H=0.74, He=0.25, metals=0.01
    mf['H'][:n_wind + n_h_shell] = 0.74
    mf['He'][:n_wind + n_h_shell] = 0.25
    mf['O'][:n_wind + n_h_shell] = 0.005
    mf['Si'][:n_wind + n_h_shell] = 0.002
    mf['Ca'][:n_wind + n_h_shell] = 0.001
    mf['Fe'][:n_wind + n_h_shell] = 0.002

    # Interaction shell: mixed composition (H-poor, some metals)
    mf['H'][i_int_start:i_int_end] = 0.10
    mf['He'][i_int_start:i_int_end] = 0.20
    mf['O'][i_int_start:i_int_end] = 0.30
    mf['Si'][i_int_start:i_int_end] = 0.20
    mf['Ca'][i_int_start:i_int_end] = 0.05
    mf['Fe'][i_int_start:i_int_end] = 0.15

    # Metal-rich ejecta: H=0, O=0.50, Si=0.30, Ca=0.02, Fe=0.18 (sums to 1.0)
    mf['H'][i_metal_start:i_metal_start + n_metal] = 0.0
    mf['He'][i_metal_start:i_metal_start + n_metal] = 0.0
    mf['O'][i_metal_start:i_metal_start + n_metal] = 0.50
    mf['Si'][i_metal_start:i_metal_start + n_metal] = 0.30
    mf['Ca'][i_metal_start:i_metal_start + n_metal] = 0.02
    mf['Fe'][i_metal_start:i_metal_start + n_metal] = 0.18

    # Inner core: same as metal-rich (SN_A ejecta)
    mf['H'][i_metal_start + n_metal:] = 0.0
    mf['He'][i_metal_start + n_metal:] = 0.0
    mf['O'][i_metal_start + n_metal:] = 0.50
    mf['Si'][i_metal_start + n_metal:] = 0.30
    mf['Ca'][i_metal_start + n_metal:] = 0.02
    mf['Fe'][i_metal_start + n_metal:] = 0.18

    # Renormalize to sum=1
    total = sum(mf.values())
    for elem in mf:
        mf[elem] = mf[elem] / np.maximum(total, 1e-30)

    return SyntheticSnapshot(
        r=r, v=v, T=T, rho=rho, kappa=kappa,
        mass_fractions=mf, time_days=35.0,
    )


def smooth_homologous(n_zones: int = 200) -> SyntheticSnapshot:
    """Simple smooth homologous model (no structure).

    Used as a baseline for regression tests.
    """
    r_outer = 1e16
    r_inner = 1e14
    log_r = np.linspace(np.log10(r_outer), np.log10(r_inner), n_zones)
    r = 10.0 ** log_r
    v = 5000.0 * (r / r_outer)  # homologous
    rho = 1e-12 * (r / r_inner) ** (-8.0)
    T = 5000.0 * (r / r_outer) ** (-0.5)
    kappa = np.full(n_zones, 0.34)
    mf = {'H': np.full(n_zones, 0.74), 'He': np.full(n_zones, 0.25),
          'O': np.full(n_zones, 0.005), 'Si': np.full(n_zones, 0.002),
          'Ca': np.full(n_zones, 0.001), 'Fe': np.full(n_zones, 0.002)}
    return SyntheticSnapshot(r=r, v=v, T=T, rho=rho, kappa=kappa,
                             mass_fractions=mf, time_days=35.0)


def steep_density_gradient(n_zones: int = 200) -> SyntheticSnapshot:
    """Model with a steep density gradient (power-law index n=12)."""
    r_outer = 1e16
    r_inner = 1e14
    log_r = np.linspace(np.log10(r_outer), np.log10(r_inner), n_zones)
    r = 10.0 ** log_r
    v = 5000.0 * (r / r_outer)
    rho = 1e-12 * (r / r_inner) ** (-12.0)  # very steep
    T = 5000.0 * (r / r_outer) ** (-0.5)
    kappa = np.full(n_zones, 0.34)
    mf = {'H': np.full(n_zones, 0.74), 'He': np.full(n_zones, 0.25),
          'O': np.full(n_zones, 0.005), 'Si': np.full(n_zones, 0.002),
          'Ca': np.full(n_zones, 0.001), 'Fe': np.full(n_zones, 0.002)}
    return SyntheticSnapshot(r=r, v=v, T=T, rho=rho, kappa=kappa,
                             mass_fractions=mf, time_days=35.0)


def thin_shell(n_zones: int = 300) -> SyntheticSnapshot:
    """Model with a narrow dense shell at mid-radius."""
    r_outer = 1e16
    r_inner = 1e14
    log_r = np.linspace(np.log10(r_outer), np.log10(r_inner), n_zones)
    r = 10.0 ** log_r
    v = 5000.0 * (r / r_outer)
    rho = 1e-12 * (r / r_inner) ** (-8.0)
    # Narrow shell at r ~ 1e15
    shell_mask = np.abs(r - 1e15) < 2e14
    rho[shell_mask] *= 50.0
    T = 5000.0 * (r / r_outer) ** (-0.5)
    T[shell_mask] *= 3.0
    kappa = np.full(n_zones, 0.34)
    kappa[shell_mask] = 1.0
    mf = {'H': np.full(n_zones, 0.74), 'He': np.full(n_zones, 0.25),
          'O': np.full(n_zones, 0.005), 'Si': np.full(n_zones, 0.002),
          'Ca': np.full(n_zones, 0.001), 'Fe': np.full(n_zones, 0.002)}
    return SyntheticSnapshot(r=r, v=v, T=T, rho=rho, kappa=kappa,
                             mass_fractions=mf, time_days=35.0)


def composition_discontinuity(n_zones: int = 300) -> SyntheticSnapshot:
    """Model with a sharp composition discontinuity at mid-radius."""
    r_outer = 1e16
    r_inner = 1e14
    log_r = np.linspace(np.log10(r_outer), np.log10(r_inner), n_zones)
    r = 10.0 ** log_r
    v = 5000.0 * (r / r_outer)
    rho = 1e-12 * (r / r_inner) ** (-8.0)
    T = 5000.0 * (r / r_outer) ** (-0.5)
    kappa = np.full(n_zones, 0.34)

    mf = {'H': np.zeros(n_zones), 'He': np.zeros(n_zones),
          'O': np.zeros(n_zones), 'Si': np.zeros(n_zones),
          'Ca': np.zeros(n_zones), 'Fe': np.zeros(n_zones)}

    # Outer: H-rich
    outer = r > 3e15
    mf['H'][outer] = 0.74
    mf['He'][outer] = 0.25
    mf['O'][outer] = 0.005
    mf['Si'][outer] = 0.002
    mf['Ca'][outer] = 0.001
    mf['Fe'][outer] = 0.002

    # Inner: metal-rich
    inner = ~outer
    mf['O'][inner] = 0.96
    mf['Si'][inner] = 0.57
    mf['Ca'][inner] = 0.036
    mf['Fe'][inner] = 0.42

    # Renormalize
    total = sum(mf.values())
    for elem in mf:
        mf[elem] = mf[elem] / np.maximum(total, 1e-30)

    return SyntheticSnapshot(r=r, v=v, T=T, rho=rho, kappa=kappa,
                             mass_fractions=mf, time_days=35.0)


# Registry of test models
TEST_MODELS = {
    'smooth_homologous': smooth_homologous,
    'steep_density_gradient': steep_density_gradient,
    'thin_shell': thin_shell,
    'composition_discontinuity': composition_discontinuity,
    'two_event_pathological': two_event_pathological,
}
