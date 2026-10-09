"""CMFGEN SN_HYDRO_DATA writer with atomic file writes.

Writes the validated target grid to SN_HYDRO_DATA in CMFGEN format
(outer-to-inner, decreasing R). Includes radioactive decay chain
(Ni56 -> Co56 -> Fe56) and isotope consistency.

The writer performs atomic writes (write to temp file, then rename)
to prevent partial output files.
"""

from __future__ import annotations

import os
import tempfile
from typing import Optional

import numpy as np

from .shell_model import TargetGrid, AVOGADRO, LSUN, MSUN
from .solver_metrics import apply_sigma_cap, sigma_solver

# Decay timescales (from STELLA source)
T_NI = 8.8    # days (Ni56 e-folding)
T_CO = 111.3  # days (Co56 e-folding)

# CMFGEN species name mapping
CMFGEN_TO_INTERNAL = {
    'HYD': 'H', 'HE': 'He', 'CARB': 'C', 'OXY': 'O',
    'SIL': 'Si', 'CAL': 'Ca', 'IRON': 'Fe', 'NICK': 'Ni', 'COB': 'Co',
}
INTERNAL_TO_CMFGEN = {v: k for k, v in CMFGEN_TO_INTERNAL.items()}

# Solar fallback abundances for missing elements
SOLAR_ABUNDANCES = {
    'CARB': 0.003, 'SIL': 0.0007, 'CAL': 0.0002,
}

ABUND_FLOOR = 1e-30


def _write_block(f, label: str, values: np.ndarray, per_line: int = 6):
    """Write a labelled data block in CMFGEN format."""
    f.write(f"{label}\n")
    for i in range(0, len(values), per_line):
        chunk = values[i:i + per_line]
        f.write(''.join(f"    {v:22.15E}" for v in chunk) + '\n')
    f.write('\n')


def compute_decay_chain(ni56_initial: np.ndarray, time_days: float) -> dict:
    """Compute Ni56 -> Co56 -> Fe56 decay chain.

    Returns dict with 'Ni56', 'Co56', 'Fe56' mass fraction arrays.
    """
    Ni56_t = ni56_initial * np.exp(-time_days / T_NI)
    Co56_t = ni56_initial * (T_CO / (T_CO - T_NI)) * (
        np.exp(-time_days / T_CO) - np.exp(-time_days / T_NI))
    Fe56_t = ni56_initial - Ni56_t - Co56_t
    return {'Ni56': Ni56_t, 'Co56': Co56_t, 'Fe56': Fe56_t}


def prepare_cmfgen_mass_fractions(target: TargetGrid,
                                  time_days: float,
                                  trace_threshold: float = 1e-4,
                                  t_floor: float = 3000.0,
                                  t_cap: float = 0.0) -> tuple[dict, dict, dict]:
    """Prepare mass fractions, isotopes, and metadata for CMFGEN output.

    Returns (cmfgen_mass_fractions, isotopic_data, metadata).
    """
    n = target.n_points
    mf = dict(target.mass_fractions)  # copy

    # Apply temperature floor/cap
    T = target.temperature_K.copy()
    if t_floor > 0:
        T = np.maximum(T, t_floor)
    if t_cap > 0:
        T = np.minimum(T, t_cap)

    # Radioactive decay
    ni56_initial = target.ni56_initial if target.ni56_initial is not None else mf.get('Ni56', np.zeros(n))
    decay = compute_decay_chain(ni56_initial, time_days)

    # Replace Ni56 with decayed value
    if 'Ni56' in mf:
        del mf['Ni56']
    mf['NICK'] = decay['Ni56']
    mf['COB'] = decay['Co56']

    # Update Fe to include decayed Fe56
    stable_fe = mf.pop('FePeak', np.zeros(n)) if 'FePeak' in mf else mf.pop('Fe', np.zeros(n))
    total_fe = stable_fe + decay['Fe56']
    mf['IRON'] = total_fe

    # Build isotopic data
    isotopic_data = {
        ('NICK', 56): mf['NICK'],
        ('COB', 56): mf['COB'],
        ('IRON', 56): mf['IRON'],
    }

    # Add missing elements with solar abundances
    for cmfgen_name, solar_frac in SOLAR_ABUNDANCES.items():
        internal_name = CMFGEN_TO_INTERNAL.get(cmfgen_name, cmfgen_name)
        if internal_name not in mf:
            mf[internal_name] = np.full(n, solar_frac)

    # Filter to CMFGEN-allowed elements
    allowed = {'H', 'He', 'C', 'O', 'Si', 'Ca', 'Fe', 'Ni', 'Co',
               'HYD', 'HE', 'CARB', 'OXY', 'SIL', 'CAL', 'IRON', 'NICK', 'COB'}
    cmfgen_mf = {}
    for name, frac in mf.items():
        if name not in allowed:
            continue
        cmfgen_name = INTERNAL_TO_CMFGEN.get(name, name.upper())
        cmfgen_mf[cmfgen_name] = frac

    # Renormalize to sum=1
    mf_sum = np.zeros(n)
    for frac in cmfgen_mf.values():
        mf_sum += frac
    mf_sum = np.where(mf_sum < 1e-30, 1.0, mf_sum)
    for elem in cmfgen_mf:
        cmfgen_mf[elem] = cmfgen_mf[elem] / mf_sum

    # Trace threshold zeroing
    if trace_threshold > 0:
        for elem in cmfgen_mf:
            cmfgen_mf[elem] = np.where(cmfgen_mf[elem] < trace_threshold, 0.0, cmfgen_mf[elem])
        # Renormalize
        mf_sum2 = np.zeros(n)
        for frac in cmfgen_mf.values():
            mf_sum2 += frac
        mf_sum2 = np.where(mf_sum2 < 1e-30, 1.0, mf_sum2)
        for elem in cmfgen_mf:
            cmfgen_mf[elem] = cmfgen_mf[elem] / mf_sum2

    # Abundance floor
    for elem in cmfgen_mf:
        cmfgen_mf[elem] = np.maximum(cmfgen_mf[elem], ABUND_FLOOR)

    # Final renormalization
    mf_sum3 = np.zeros(n)
    for frac in cmfgen_mf.values():
        mf_sum3 += frac
    mf_sum3 = np.where(mf_sum3 < 1e-30, 1.0, mf_sum3)
    for elem in cmfgen_mf:
        cmfgen_mf[elem] = cmfgen_mf[elem] / mf_sum3

    # Force isotope = elemental for consistency
    for iso_key, elem_name in [(('NICK', 56), 'NICK'), (('COB', 56), 'COB'), (('IRON', 56), 'IRON')]:
        if iso_key in isotopic_data and elem_name in cmfgen_mf:
            isotopic_data[iso_key] = cmfgen_mf[elem_name].copy()

    metadata = {
        'nfrac': len(cmfgen_mf),
        'n_isotopes': len(isotopic_data),
        'T': T,
        'mf_sum_before_norm': float(np.mean(mf_sum)),
    }

    return cmfgen_mf, isotopic_data, metadata


def write_sn_hydro_data(outpath: str, target: TargetGrid,
                        lum_lsun: float = 0.0,
                        t_floor: float = 3000.0,
                        t_cap: float = 0.0,
                        trace_threshold: float = 1e-4,
                        r_inner_cm: Optional[float] = None,
                        T_inner_K: Optional[float] = None,
                        tau_inner: Optional[float] = None,
                        preserve_v_sign: bool = True,
                        vel_law: int = 12) -> dict:
    """Write CMFGEN SN_HYDRO_DATA file (atomic write via temp file).

    Parameters
    ----------
    vel_law : 11 or 12 (project default 12).  Selects the SIGMA formula used
        for the written Sigma column.  Under VEL_LAW=11 this column is read
        by CMFGEN and used for the density time-evolution correction, so it
        must be consistent with the V(R) the solver will interpolate.  Under
        VEL_LAW=12 the solver recomputes sigma itself from input_hydro.dat,
        so the column is diagnostic only.  In neither case does this column
        directly feed the radiative-transfer solver sigma.

    Returns metadata dict with output statistics.
    """
    n = target.n_points
    if n < 2:
        raise ValueError(f"Target grid has only {n} points; need at least 2")

    # Prepare mass fractions and isotopes
    cmfgen_mf, isotopic_data, meta = prepare_cmfgen_mass_fractions(
        target, target.time_days, trace_threshold, t_floor, t_cap)
    T = meta['T']

    # Recompute SIGMA with the formula matching the target VEL_LAW so the
    # written column is consistent with what CMFGEN will derive from V(R).
    # (The solver recomputes sigma from V(R) in both paths; this column is
    # used for density evolution under VEL_LAW=11 and is diagnostic under 12.)
    sigma_raw = sigma_solver(target.r_cm, target.velocity_cms,
                             vel_law=vel_law, preserve_sign=preserve_v_sign)
    sigma = apply_sigma_cap(sigma_raw)

    # Velocity in km/s (preserve sign)
    v_kms = target.velocity_cms / 1e5

    # Radius in 10^10 cm
    r_1e10 = target.r_cm / 1e10

    # Temperature in 10^4 K
    T_10k = T / 1e4

    # Mean atomic mass and electron density
    elem_A = {'H': 1.0, 'He': 4.0, 'C': 12.0, 'N': 14.0, 'O': 16.0,
              'Ne': 20.0, 'Na': 23.0, 'Mg': 24.0, 'Al': 27.0, 'Si': 28.0,
              'S': 32.0, 'Ar': 40.0, 'Ca': 40.0, 'Fe': 56.0, 'Ni': 56.0,
              'Co': 56.0, 'FePeak': 56.0, 'Ni56': 56.0, 'Ni58': 58.0}
    elem_Z = {'H': 1.0, 'He': 2.0, 'C': 6.0, 'N': 7.0, 'O': 8.0,
              'Ne': 10.0, 'Na': 11.0, 'Mg': 12.0, 'Al': 13.0, 'Si': 14.0,
              'S': 16.0, 'Ar': 18.0, 'Ca': 20.0, 'Fe': 26.0, 'Ni': 28.0,
              'Co': 27.0, 'FePeak': 26.0, 'Ni56': 28.0, 'Ni58': 28.0}

    inv_A = np.zeros(n)
    Z_over_A = np.zeros(n)
    for elem, frac in cmfgen_mf.items():
        internal = CMFGEN_TO_INTERNAL.get(elem, elem)
        A = elem_A.get(internal, 28.0)
        Z = elem_Z.get(internal, 14.0)
        inv_A += frac / A
        Z_over_A += frac * Z / A

    atom_density = target.density_gcm3 * AVOGADRO * inv_A

    # Electron density
    if target.electron_density is not None:
        electron_density = np.clip(target.electron_density, 1.0, None)
    else:
        # Saha fallback
        from .ionization_metrics import saha_electron_density
        electron_density = saha_electron_density(target.density_gcm3, T, cmfgen_mf)

    # Rosseland opacity in 10^{-10} cm^{-1}
    rosseland_chi = target.kappa_cgs * target.density_gcm3 * 1e10

    # Time
    time_days = target.time_days
    if time_days < 0:
        time_days = abs(time_days)

    # Atomic write
    dir_path = os.path.dirname(outpath) or '.'
    fd, tmp_path = tempfile.mkstemp(dir=dir_path, suffix='.tmp')
    try:
        with os.fdopen(fd, 'w') as f:
            f.write(f"Number of data points:{n:>10d}\n")
            f.write(f"Number of mass fractions:{meta['nfrac']:>10d}\n")
            f.write(f"Number of isotopes:{meta['n_isotopes']:>10d}\n")
            f.write(f"Time(days) since explosion:{time_days:>15.7f}\n")
            f.write(f"Energy Lum [lsun]:{lum_lsun:>19.6E}\n")
            if r_inner_cm is not None:
                f.write(f"Inner boundary radius [cm]:{r_inner_cm:>19.6E}\n")
            if T_inner_K is not None:
                f.write(f"Inner boundary temperature [K]:{T_inner_K:>19.6E}\n")
            if tau_inner is not None:
                f.write(f"Inner boundary Rosseland tau:{tau_inner:>19.6E}\n")
            f.write('\n')

            _write_block(f, "Radius grid (10^10cm)", r_1e10)
            _write_block(f, "Velocity (km/s)", v_kms)
            _write_block(f, "Sigma (dlnV/dlnr-1)", sigma)
            _write_block(f, "Temperature (10^4 K)", T_10k)
            _write_block(f, "Density (gm/cm^3)", target.density_gcm3)
            _write_block(f, "Atom density (/cm^3)", atom_density)
            _write_block(f, "Electron density (/cm^3)", electron_density)
            _write_block(f, "Rosseland mean opacty (10^{-10} cm^{-1})", rosseland_chi)
            _write_block(f, "Kappa (cm^2/gm)", target.kappa_cgs)

            for elem, frac in cmfgen_mf.items():
                _write_block(f, f"{elem} mass fraction", frac)

            # Isotopic data
            for (iso_name, baryon_num), iso_frac in isotopic_data.items():
                f.write(f"{iso_name} {baryon_num} mass fraction\n")
                for i in range(0, len(iso_frac), 8):
                    chunk = iso_frac[i:i + 8]
                    f.write(''.join(f"    {v:14.7E}" for v in chunk) + '\n')
                f.write('\n')

        os.rename(tmp_path, outpath)
    except Exception:
        if os.path.exists(tmp_path):
            os.unlink(tmp_path)
        raise

    print(f"Written {outpath}  ({n} depth points, {meta['nfrac']} mass fractions)")
    print(f"  t={time_days:.3f} days, L={lum_lsun:.3e} Lsun")
    print(f"  r=[{r_1e10[0]:.3e},{r_1e10[-1]:.3e}] 1e10cm")
    print(f"  v=[{v_kms[0]:.1f},{v_kms[-1]:.1f}] km/s")

    return {
        'n_points': n,
        'n_mass_fractions': meta['nfrac'],
        'n_isotopes': meta['n_isotopes'],
        'time_days': time_days,
        'lum_lsun': lum_lsun,
        'r_range_1e10': (float(r_1e10[0]), float(r_1e10[-1])),
        'v_range_kms': (float(v_kms[0]), float(v_kms[-1])),
        'vel_law': vel_law,
        'sigma_formula': 'cubic_interp_derivative' if vel_law == 11
                         else 'finite_difference_weighted',
        'sigma_solver_note': ('SN_HYDRO_DATA Sigma column does NOT directly '
                              'feed the RT solver; solver recomputes sigma '
                              'from V(R).'),
    }
