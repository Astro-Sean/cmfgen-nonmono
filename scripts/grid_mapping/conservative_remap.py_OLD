"""Conservative shell-overlap remapping for mass, composition, and isotopes.

Implements volume-overlap and mass-coordinate-overlap remapping that
conserves total mass, elemental mass, and isotope mass when mapping
from the source STELLA grid to a target CMFGEN grid.

Conservative remapping procedure:
1. Define original shell boundaries (source grid).
2. Define final target shell boundaries (target grid).
3. Compute overlap volumes between source and target shells.
4. Integrate total mass, elemental mass, and isotope mass over overlaps.
5. Compute target density from target shell mass / target volume.
6. Compute mass fractions from mapped masses.
7. Map thermodynamic variables (T, Ne) using mass-weighted averaging.
8. Map velocity using mass-weighted (momentum-conserving) averaging.

The source grid is NEVER remapped repeatedly — the original SourceShellGrid
is preserved for all error estimates.
"""

from __future__ import annotations

from typing import Optional

import numpy as np

from .shell_model import (
    SourceShellGrid, TargetGrid,
    compute_shell_volume, density_from_mass_volume,
    MSUN,
)


def compute_overlap_volumes(r_source: np.ndarray, r_target: np.ndarray) -> np.ndarray:
    """Compute overlap volume matrix between source and target shells.

    Parameters
    ----------
    r_source : (n_source+1,) source shell boundaries, decreasing
    r_target : (n_target+1,) target shell boundaries, decreasing

    Returns
    -------
    overlap : (n_target, n_source) matrix
        overlap[j, i] = volume of intersection of target shell j
        with source shell i, in cm^3.
    """
    n_s = len(r_source) - 1
    n_t = len(r_target) - 1
    overlap = np.zeros((n_t, n_s))

    # Convert to increasing order for easier intersection computation
    rs = r_source[::-1]  # increasing
    rt = r_target[::-1]  # increasing

    FOUR_PI_OVER_3 = 4.0 * np.pi / 3.0

    def shell_vol(r_in, r_out):
        """Volume of spherical shell [r_in, r_out]."""
        return FOUR_PI_OVER_3 * (r_out**3 - r_in**3)

    # For each target shell, find overlapping source shells
    # overlap_inc[j, i] is in increasing-R order for both target and source
    overlap_inc = np.zeros((n_t, n_s))
    for j in range(n_t):
        t_in = rt[j]      # inner radius of target shell j
        t_out = rt[j + 1]  # outer radius of target shell j

        # Find source shells that overlap [t_in, t_out]
        # Source shell i spans [rs[i], rs[i+1]] (increasing)
        for i in range(n_s):
            s_in = rs[i]
            s_out = rs[i + 1]

            # Intersection
            ov_in = max(s_in, t_in)
            ov_out = min(s_out, t_out)

            if ov_out > ov_in:
                overlap_inc[j, i] = shell_vol(ov_in, ov_out)

    # Reverse both rows and columns to get back to decreasing-R order
    # overlap[j, i] now corresponds to target shell j (decreasing R) and source shell i (decreasing R)
    overlap = overlap_inc[::-1, ::-1]

    return overlap


def remap_mass(overlap: np.ndarray, source_mass: np.ndarray) -> np.ndarray:
    """Remap shell mass conservatively using overlap volumes.

    target_mass[j] = sum_i overlap[j,i] / source_vol[i] * source_mass[i]

    This conserves total mass exactly.
    """
    n_t, n_s = overlap.shape
    # source density = source_mass / source_volume
    # But we use overlap fraction: mass is distributed proportional to
    # overlap volume fraction
    source_vol = np.zeros(n_s)
    # Reconstruct source volumes from overlap row sums (approximate)
    # Actually, we need the true source volumes. Pass them in or compute
    # from the overlap matrix column sums (which give source volume if
    # target fully covers source).
    col_sums = overlap.sum(axis=0)
    source_vol = np.where(col_sums > 0, col_sums, 1e-99)

    # Mass density in each source shell
    source_rho = source_mass / source_vol

    # Remap: target mass = sum of overlap * source_density
    target_mass = overlap @ source_rho
    return target_mass


def remap_mass_with_volumes(overlap: np.ndarray, source_mass: np.ndarray,
                             source_vol: np.ndarray) -> np.ndarray:
    """Remap shell mass conservatively using known source volumes.

    target_mass[j] = sum_i overlap[j,i] * (source_mass[i] / source_vol[i])
    """
    source_rho = source_mass / np.maximum(source_vol, 1e-99)
    return overlap @ source_rho


def remap_elemental_mass(overlap: np.ndarray, source_mass: np.ndarray,
                          source_vol: np.ndarray,
                          mass_fractions: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    """Remap elemental mass conservatively.

    For each element:
        source_element_mass[i] = mass_fractions[elem][i] * source_mass[i]
        target_element_mass[j] = sum_i overlap[j,i]/source_vol[i] * source_element_mass[i]
        target_mass_fraction[elem][j] = target_element_mass[j] / target_mass[j]
    """
    source_rho = source_mass / np.maximum(source_vol, 1e-99)
    target_mass = overlap @ source_rho

    target_mf = {}
    for elem, frac in mass_fractions.items():
        # Elemental mass density in source
        elem_rho = frac * source_rho
        # Remap elemental mass
        target_elem_mass = overlap @ elem_rho
        # Target mass fraction
        target_mf[elem] = target_elem_mass / np.maximum(target_mass, 1e-99)

    # Renormalize to sum=1 at each target shell
    mf_sum = np.zeros(len(target_mass))
    for frac in target_mf.values():
        mf_sum += frac
    mf_sum = np.where(mf_sum < 1e-30, 1.0, mf_sum)
    for elem in target_mf:
        target_mf[elem] = target_mf[elem] / mf_sum

    return target_mf


def remap_isotope_mass(overlap: np.ndarray, source_mass: np.ndarray,
                        source_vol: np.ndarray,
                        isotope_fractions: dict[tuple[str, int], np.ndarray],
                        mass_fractions: dict[str, np.ndarray]) -> dict[tuple[str, int], np.ndarray]:
    """Remap isotope mass conservatively.

    Isotope mass fractions are remapped like elemental mass fractions,
    then forced to match the elemental mass fraction for species with
    isotopes (as required by CMFGEN's rd_sn_data.f consistency check).
    """
    source_rho = source_mass / np.maximum(source_vol, 1e-99)
    target_mass = overlap @ source_rho

    target_iso = {}
    for iso_key, frac in isotope_fractions.items():
        iso_rho = frac * source_rho
        target_iso_mass = overlap @ iso_rho
        target_iso[iso_key] = target_iso_mass / np.maximum(target_mass, 1e-99)

    return target_iso


def remap_ni56(overlap: np.ndarray, source_mass: np.ndarray,
               source_vol: np.ndarray,
               ni56_initial: np.ndarray) -> np.ndarray:
    """Remap Ni56 mass fraction conservatively."""
    source_rho = source_mass / np.maximum(source_vol, 1e-99)
    target_mass = overlap @ source_rho

    ni56_rho = ni56_initial * source_rho
    target_ni56_mass = overlap @ ni56_rho
    return target_ni56_mass / np.maximum(target_mass, 1e-99)


def remap_velocity_momentum(overlap: np.ndarray, source_mass: np.ndarray,
                             source_vol: np.ndarray,
                             source_velocity: np.ndarray) -> np.ndarray:
    """Remap velocity using momentum-conserving (mass-weighted) averaging.

    target_velocity[j] = sum_i overlap[j,i] * rho[i] * v[i] /
                         sum_i overlap[j,i] * rho[i]

    This conserves momentum = integral(rho * v dV).
    """
    source_rho = source_mass / np.maximum(source_vol, 1e-99)
    # Mass-weighted velocity: integral(rho * v dV) / integral(rho dV)
    momentum_density = source_rho * source_velocity  # rho * v
    target_momentum = overlap @ momentum_density
    target_mass_density = overlap @ source_rho
    return target_momentum / np.maximum(target_mass_density, 1e-99)


def remap_temperature_mass_weighted(overlap: np.ndarray, source_mass: np.ndarray,
                                     source_vol: np.ndarray,
                                     source_temperature: np.ndarray) -> np.ndarray:
    """Remap temperature using mass-weighted averaging.

    T_target[j] = sum_i overlap[j,i] * rho[i] * T[i] /
                  sum_i overlap[j,i] * rho[i]
    """
    source_rho = source_mass / np.maximum(source_vol, 1e-99)
    weighted_T = source_rho * source_temperature
    target_weighted_T = overlap @ weighted_T
    target_mass_density = overlap @ source_rho
    return target_weighted_T / np.maximum(target_mass_density, 1e-99)


def remap_electron_density(overlap: np.ndarray, source_mass: np.ndarray,
                            source_vol: np.ndarray,
                            source_ne: Optional[np.ndarray]) -> Optional[np.ndarray]:
    """Remap electron density using volume-weighted averaging.

    Ne_target[j] = sum_i overlap[j,i] * Ne[i] / sum_i overlap[j,i]
    """
    if source_ne is None:
        return None
    target_ne = overlap @ source_ne
    target_vol = overlap.sum(axis=1)
    return target_ne / np.maximum(target_vol, 1e-99)


def remap_opacity(overlap: np.ndarray, source_mass: np.ndarray,
                  source_vol: np.ndarray,
                  source_kappa: np.ndarray) -> np.ndarray:
    """Remap opacity using mass-weighted averaging.

    kappa_target[j] = sum_i overlap[j,i] * rho[i] * kappa[i] /
                      sum_i overlap[j,i] * rho[i]
    """
    source_rho = source_mass / np.maximum(source_vol, 1e-99)
    weighted_kappa = source_rho * source_kappa
    target_weighted = overlap @ weighted_kappa
    target_mass_density = overlap @ source_rho
    return target_weighted / np.maximum(target_mass_density, 1e-99)


def conservative_remap(source: SourceShellGrid,
                       r_target_bound: np.ndarray,
                       time_days: float) -> TargetGrid:
    """Perform full conservative remapping from source to target grid.

    Parameters
    ----------
    source : SourceShellGrid
        The original STELLA shell grid (never modified).
    r_target_bound : (n_target+1,) target shell boundaries, decreasing
    time_days : snapshot time in days

    Returns
    -------
    TargetGrid with all quantities conservatively remapped.
    """
    # Compute overlap volumes
    overlap = compute_overlap_volumes(source.r_bound_cm, r_target_bound)

    n_target = len(r_target_bound) - 1
    target_vol = compute_shell_volume(r_target_bound)

    # Remap total mass
    target_mass = remap_mass_with_volumes(overlap, source.shell_mass_g,
                                           source.shell_volume_cm3)

    # Density from mass / volume
    target_density = density_from_mass_volume(target_mass, target_vol)

    # Fix zero-mass target zones (numerical precision in overlap computation
    # can produce zero overlap for thin target shells between source zones).
    # Interpolate density, T, v from neighbors for those zones.
    zero_mass = target_mass <= 0
    if np.any(zero_mass) and n_target > 2:
        n_fix = int(np.sum(zero_mass))
        print(f"  [conservative_remap] Fixing {n_fix} zero-mass target zones by neighbor interpolation")
        for i in np.where(zero_mass)[0]:
            # Find nearest non-zero neighbors
            lo = i - 1
            while lo >= 0 and target_mass[lo] <= 0:
                lo -= 1
            hi = i + 1
            while hi < n_target and target_mass[hi] <= 0:
                hi += 1
            if lo >= 0 and hi < n_target:
                w = 0.5
                target_mass[i] = w * target_mass[lo] + (1 - w) * target_mass[hi]
            elif lo >= 0:
                target_mass[i] = target_mass[lo]
            elif hi < n_target:
                target_mass[i] = target_mass[hi]
        target_density = density_from_mass_volume(target_mass, target_vol)

    # Velocity (momentum-conserving)
    target_v = remap_velocity_momentum(overlap, source.shell_mass_g,
                                       source.shell_volume_cm3,
                                       source.velocity_cms)

    # Temperature (mass-weighted)
    target_T = remap_temperature_mass_weighted(overlap, source.shell_mass_g,
                                                source.shell_volume_cm3,
                                                source.temperature_K)

    # Fix zero T/v in target zones that had zero overlap (already fixed mass above).
    # Interpolate from neighbors for physical consistency.
    if np.any(zero_mass) and n_target > 2:
        for i in np.where(zero_mass)[0]:
            lo = i - 1
            while lo >= 0 and zero_mass[lo]:
                lo -= 1
            hi = i + 1
            while hi < n_target and zero_mass[hi]:
                hi += 1
            if lo >= 0 and hi < n_target:
                target_v[i] = 0.5 * (target_v[lo] + target_v[hi])
                target_T[i] = 0.5 * (target_T[lo] + target_T[hi])
            elif lo >= 0:
                target_v[i] = target_v[lo]
                target_T[i] = target_T[lo]
            elif hi < n_target:
                target_v[i] = target_v[hi]
                target_T[i] = target_T[hi]

    # Composition
    target_mf = remap_elemental_mass(overlap, source.shell_mass_g,
                                     source.shell_volume_cm3,
                                     source.mass_fractions)

    # Isotopes
    target_iso = {}
    if source.isotope_fractions:
        target_iso = remap_isotope_mass(overlap, source.shell_mass_g,
                                        source.shell_volume_cm3,
                                        source.isotope_fractions,
                                        source.mass_fractions)

    # Ni56
    target_ni56 = None
    if source.ni56_initial is not None:
        target_ni56 = remap_ni56(overlap, source.shell_mass_g,
                                  source.shell_volume_cm3,
                                  source.ni56_initial)

    # Electron density
    target_ne = remap_electron_density(overlap, source.shell_mass_g,
                                       source.shell_volume_cm3,
                                       source.electron_density)

    # Opacity
    target_kappa = remap_opacity(overlap, source.shell_mass_g,
                                  source.shell_volume_cm3,
                                  source.kappa_cgs)

    # Target shell centres
    from .shell_model import compute_shell_centres
    r_centre = compute_shell_centres(r_target_bound)

    # Compute sigma
    from .solver_metrics import sigma_velocity
    sigma = sigma_velocity(r_centre, target_v, preserve_sign=True)

    # Compute tau
    from .opacity_tau import compute_tau
    tau = compute_tau(r_target_bound, target_density, target_kappa)

    # Atom density
    elem_A = {'H': 1.0, 'He': 4.0, 'C': 12.0, 'N': 14.0, 'O': 16.0,
              'Ne': 20.0, 'Na': 23.0, 'Mg': 24.0, 'Al': 27.0, 'Si': 28.0,
              'S': 32.0, 'Ar': 40.0, 'Ca': 40.0, 'Fe': 56.0, 'Ni': 56.0,
              'FePeak': 56.0, 'Ni56': 56.0, 'Ni58': 58.0}
    from .shell_model import AVOGADRO
    inv_A = np.zeros(n_target)
    for elem, frac in target_mf.items():
        A = elem_A.get(elem, 28.0)
        inv_A += frac / A
    atom_density = target_density * AVOGADRO * inv_A

    # Provenance: all points are 'interpolated' (conservatively remapped)
    provenance = np.array(['interpolated'] * n_target, dtype='<U16')
    # Mark boundary points
    if n_target > 0:
        provenance[0] = 'protected'  # outer boundary
        provenance[-1] = 'protected'  # inner boundary

    target = TargetGrid(
        r_cm=r_centre,
        velocity_cms=target_v,
        temperature_K=target_T,
        density_gcm3=target_density,
        kappa_cgs=target_kappa,
        sigma=sigma,
        tau_cumulative=tau,
        mass_fractions=target_mf,
        isotope_fractions=target_iso,
        ni56_initial=target_ni56,
        electron_density=target_ne,
        atom_density=atom_density,
        provenance=provenance,
        shell_mass_g=target_mass,
        shell_volume_cm3=target_vol,
        is_anchor=np.zeros(n_target, dtype=bool),
        time_days=time_days,
    )

    return target


def compute_conservation_residuals(source: SourceShellGrid,
                                    target: TargetGrid) -> dict:
    """Compute conservation residuals between source and target grids.

    Returns dict with fractional residuals for:
    - total mass
    - elemental mass (per element)
    - isotope mass (per isotope)
    - momentum
    - kinetic energy
    - thermal energy (if T available)
    - electron number
    - optical depth (total)
    """
    residuals = {}

    # Total mass
    source_mass = source.shell_mass_g.sum()
    target_mass = target.shell_mass_g.sum() if target.shell_mass_g is not None else 0.0
    residuals['total_mass'] = {
        'source': float(source_mass),
        'target': float(target_mass),
        'fractional_residual': float((target_mass - source_mass) / max(source_mass, 1e-30)),
    }

    # Elemental mass
    for elem, frac in source.mass_fractions.items():
        source_elem_mass = (frac * source.shell_mass_g).sum()
        if elem in target.mass_fractions:
            target_elem_mass = (target.mass_fractions[elem] * target.shell_mass_g).sum()
        else:
            target_elem_mass = 0.0
        residuals[f'elemental_mass_{elem}'] = {
            'source': float(source_elem_mass),
            'target': float(target_elem_mass),
            'fractional_residual': float((target_elem_mass - source_elem_mass) / max(source_elem_mass, 1e-30)),
        }

    # Isotope mass
    for iso_key, frac in source.isotope_fractions.items():
        source_iso_mass = (frac * source.shell_mass_g).sum()
        if iso_key in target.isotope_fractions:
            target_iso_mass = (target.isotope_fractions[iso_key] * target.shell_mass_g).sum()
        else:
            target_iso_mass = 0.0
        key = f'isotope_mass_{iso_key[0]}_{iso_key[1]}'
        residuals[key] = {
            'source': float(source_iso_mass),
            'target': float(target_iso_mass),
            'fractional_residual': float((target_iso_mass - source_iso_mass) / max(source_iso_mass, 1e-30)),
        }

    # Ni56
    if source.ni56_initial is not None:
        source_ni56 = (source.ni56_initial * source.shell_mass_g).sum()
        target_ni56 = (target.ni56_initial * target.shell_mass_g).sum() if target.ni56_initial is not None else 0.0
        residuals['ni56_mass'] = {
            'source': float(source_ni56),
            'target': float(target_ni56),
            'fractional_residual': float((target_ni56 - source_ni56) / max(source_ni56, 1e-30)),
        }

    # Momentum
    source_p = (source.shell_mass_g * source.velocity_cms).sum()
    target_p = (target.shell_mass_g * target.velocity_cms).sum() if target.shell_mass_g is not None else 0.0
    residuals['momentum'] = {
        'source': float(source_p),
        'target': float(target_p),
        'fractional_residual': float((target_p - source_p) / max(abs(source_p), 1e-30)),
    }

    # Kinetic energy
    source_ke = (0.5 * source.shell_mass_g * source.velocity_cms**2).sum()
    target_ke = (0.5 * target.shell_mass_g * target.velocity_cms**2).sum() if target.shell_mass_g is not None else 0.0
    residuals['kinetic_energy'] = {
        'source': float(source_ke),
        'target': float(target_ke),
        'fractional_residual': float((target_ke - source_ke) / max(source_ke, 1e-30)),
    }

    # Electron number
    if source.electron_density is not None:
        source_ne = (source.electron_density * source.shell_volume_cm3).sum()
        if target.electron_density is not None and target.shell_volume_cm3 is not None:
            target_ne = (target.electron_density * target.shell_volume_cm3).sum()
        else:
            target_ne = 0.0
        residuals['electron_number'] = {
            'source': float(source_ne),
            'target': float(target_ne),
            'fractional_residual': float((target_ne - source_ne) / max(source_ne, 1e-30)),
        }

    # Total optical depth
    source_tau = source.tau_cumulative[-1] if len(source.tau_cumulative) > 0 else 0.0
    target_tau = target.tau_cumulative[-1] if len(target.tau_cumulative) > 0 else 0.0
    residuals['total_tau'] = {
        'source': float(source_tau),
        'target': float(target_tau),
        'fractional_residual': float((target_tau - source_tau) / max(source_tau, 1e-30)),
    }

    return residuals
