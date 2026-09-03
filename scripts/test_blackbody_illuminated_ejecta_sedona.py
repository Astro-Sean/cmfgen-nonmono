#!/usr/bin/env python3
"""Standalone SEDONA test: SN_A ejecta at 400d + 7000K inner blackbody.

This is a standalone diagnostic that does NOT use run_model_grid.py.
It reads the SN_A STELLA .swd, homologously expands to 400 days,
generates a SEDONA .mod file with a 7000K inner source, and runs SEDONA
directly via srun.

The purpose is to compare SEDONA's spectral synthesis (with full NLTE
treatment for all elements) against SuperLite's limited NLTE (H only).
"""

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
PI = np.pi
DAY_S = 86400.0
M_SUN = 1.989e33
SIGMA_SB = 5.670374419e-5
C_LIGHT = 2.998e10

# SN_A source
SN_A_SWD = '/viper/ptmp1/sbrennan/SN2024pba/SN_A/SN_A_g_T20k_Ni06/stella/SN_A_g_T20k_Ni06.swd'

# SEDONA executable
SEDONA_EXE = '/ptmp/sbrennan/sedona/src/sedona6.ex'
SEDONA_HOME = '/ptmp/sbrennan/sedona'

# Test parameters
T_BB = 7000.0       # inner blackbody temperature [K]
T_OUTER = 2000.0    # outer H/He material temperature [K]
T_FLOOR = 2000.0    # minimum temperature [K]
TARGET_AGE_D = 400.0  # homologous expansion target [days]
MAX_CELLS = 250     # downsample target

# Elements to include in the .mod file
# SEDONA .mod format: H He C O Si Ca Fe (matching the 7-element scheme)
ALL_ELEMENTS = ['H', 'He', 'C', 'O', 'Si', 'Ca', 'Fe']
ELEMENT_Z = {'H': 1, 'He': 2, 'C': 6, 'O': 8, 'Si': 14, 'Ca': 20, 'Fe': 26}


# ---------------------------------------------------------------------------
# Step 1: Read SN_A .swd
# ---------------------------------------------------------------------------
def read_sn_a_profile(swd_path):
    """Read the latest epoch from the SN_A .swd file."""
    sys.path.insert(0, str(Path(__file__).parent))
    from convert_stella_swd_to_cmfgen import read_swd

    data = read_swd(swd_path)
    times = sorted([float(k) for k in data.keys()])
    last_t = times[-1]
    epoch = data[last_t]

    print(f"Read SN_A .swd: {swd_path}")
    print(f"  Epochs: {len(times)}, using last: t={last_t:.4f} d")
    print(f"  Zones: {len(epoch['r'])}")
    print(f"  R: {epoch['r'].min():.3e} - {epoch['r'].max():.3e} cm")
    print(f"  V: {epoch['v'].min():.1f} - {epoch['v'].max():.1f} km/s")
    print(f"  T_gas: {epoch['T'].min():.1f} - {epoch['T'].max():.1f} K")
    print(f"  rho: {epoch['rho'].min():.3e} - {epoch['rho'].max():.3e} g/cm^3")

    return epoch, last_t


# ---------------------------------------------------------------------------
# Step 2: Homologous expansion to 400 days
# ---------------------------------------------------------------------------
def homologous_expand(epoch, t_snap_d, t_target_d):
    """Homologously expand the profile from t_snap to t_target."""
    scale = t_target_d / t_snap_d
    r_new = epoch['r'] * scale
    rho_new = epoch['rho'] / scale**3
    v_new = r_new / (t_target_d * DAY_S)  # cm/s

    # Mass
    r_inner = r_new[::-1]
    rho_inner = rho_new[::-1]
    r_edges = np.concatenate([[0], r_inner])
    vol = 4.0/3.0 * PI * (r_edges[1:]**3 - r_edges[:-1]**3)
    mass = vol * rho_inner
    total_mass = mass.sum()

    print(f"\nHomologous expansion to {t_target_d:.0f} d (scale={scale:.3f}):")
    print(f"  R: {r_new.min():.3e} - {r_new.max():.3e} cm")
    print(f"  V: {v_new.min()/1e5:.1f} - {v_new.max()/1e5:.1f} km/s")
    print(f"  rho: {rho_new.min():.3e} - {rho_new.max():.3e} g/cm^3")
    print(f"  Total mass: {total_mass/M_SUN:.4f} Msun")

    return r_new, rho_new, v_new, total_mass


# ---------------------------------------------------------------------------
# Step 3: Detect elements
# ---------------------------------------------------------------------------
def detect_elements(composition):
    """Detect which elements have non-zero abundance."""
    # STELLA .swd composition columns: H, He, C, N, O, Ne, Na, Mg, Al, Si, S, Ar, Ca, Fe, Ni
    stella_elems = ['H', 'He', 'C', 'N', 'O', 'Ne', 'Na', 'Mg', 'Al',
                    'Si', 'S', 'Ar', 'Ca', 'Fe', 'Ni']
    present = []
    for i, name in enumerate(stella_elems):
        if i < composition.shape[1] and composition[:, i].max() > 1e-10:
            if name in ALL_ELEMENTS:
                present.append(name)
    return present


# ---------------------------------------------------------------------------
# Step 4: Downsample grid
# ---------------------------------------------------------------------------
def downsample_grid(r, rho, v, comp, T_gas, T_rad, elems, max_cells):
    """Downsample by taking every Nth zone, preserving structure."""
    n = len(r)
    if n <= max_cells:
        return r[::-1], rho[::-1], v[::-1], comp[::-1], T_gas[::-1], T_rad[::-1], n

    stride = max(1, n // max_cells)
    # .swd is ordered outer-to-inner; reverse to inner-to-outer
    r_io = r[::-1]
    rho_io = rho[::-1]
    v_io = v[::-1]
    comp_io = comp[::-1]
    T_gas_io = T_gas[::-1]
    T_rad_io = T_rad[::-1]

    # Take every stride-th zone
    r_ds = r_io[::stride]
    rho_ds = rho_io[::stride]
    v_ds = v_io[::stride]
    comp_ds = comp_io[::stride]
    T_gas_ds = T_gas_io[::stride]
    T_rad_ds = T_rad_io[::stride]

    print(f"\nDownsampled: {n} -> {len(r_ds)} zones (stride~{stride})")
    return r_ds, rho_ds, v_ds, comp_ds, T_gas_ds, T_rad_ds, n


# ---------------------------------------------------------------------------
# Step 5: Set temperature profile
# ---------------------------------------------------------------------------
def set_temperature_profile(r_io, comp_io, t_bb, t_outer, t_floor):
    """Set gas temperatures: inner 7000K source, outer 2000K absorber."""
    n = len(r_io)
    h_idx = ALL_ELEMENTS.index('H')
    is_metal = comp_io[:, h_idx] < 0.01

    transition_idx = np.where(is_metal)[0]
    if len(transition_idx) == 0:
        t_gas = np.full(n, t_bb)
    else:
        last_metal = transition_idx[-1]
        ramp_width = min(20, n - last_metal - 1)
        t_gas = np.full(n, t_outer)
        t_gas[:last_metal + 1] = t_bb
        if ramp_width > 0:
            for i in range(last_metal + 1, last_metal + 1 + ramp_width):
                if i < n:
                    frac = (i - last_metal) / ramp_width
                    t_gas[i] = t_bb * (1 - frac) + t_outer * frac

    t_gas = np.maximum(t_gas, t_floor)

    n_hot = np.sum(t_gas >= t_bb - 1)
    n_outer = np.sum(t_gas <= t_outer + 1)
    n_ramp = n - n_hot - n_outer

    print(f"\nTemperature profile:")
    print(f"  Hot (7000 K) zones:   {n_hot} (inner metal-rich source)")
    print(f"  Outer ({t_outer:.0f} K) zones: {n_outer} (H/He absorber)")
    print(f"  Ramp zones:            {n_ramp}")
    print(f"  T_gas range: {t_gas.min():.0f} - {t_gas.max():.0f} K")

    return t_gas


# ---------------------------------------------------------------------------
# Step 6: Write SEDONA .mod file
# ---------------------------------------------------------------------------
def write_sedona_mod(mod_path, r_io, v_io, rho_io, t_gas, comp_io, elems,
                     t_day, l_bol):
    """Write SEDONA 1D_sphere .mod file.

    Format:
      Line 1: '1D_sphere standard'
      Line 2: nzones  r_inner  t_exp  n_elements
      Line 3: Z1 Z2 Z3 ... (atomic numbers)
      Lines 4+: r_outer  v  rho  T  X1  X2  X3 ...

    r_inner is the inner boundary radius (center of star).
    Data lines give the OUTER radius of each zone.
    Zone 1 goes from r_inner to r_io[0].
    """
    n = len(r_io)
    nelem = len(elems)
    t_exp = t_day * DAY_S  # seconds

    # Inner boundary: use a small fraction of the first zone radius
    # so zone 1 has non-zero width
    r_inner = r_io[0] * 0.99

    # Map composition to the elements we're including
    elem_indices = [ALL_ELEMENTS.index(e) for e in elems]
    z_list = [ELEMENT_Z[e] for e in elems]

    with open(mod_path, 'w') as f:
        f.write('1D_sphere standard\n')
        f.write(f'{n}\t{r_inner:.8e}\t{t_exp:.1f}\t{nelem}\n')
        f.write(' '.join(str(z) for z in z_list) + '\n')

        for i in range(n):
            # r_outer, v (cm/s), rho, T, abundances
            vals = [f'{r_io[i]:.8e}', f'{v_io[i]:.8e}',
                    f'{rho_io[i]:.8e}', f'{t_gas[i]:.8e}']
            for ei in elem_indices:
                vals.append(f'{comp_io[i, ei]:.6e}')
            f.write(' '.join(vals) + '\n')

    print(f"\nWrote {mod_path}: {n} cells, {nelem} elements ({', '.join(elems)})")
    print(f"  r_inner = {r_inner:.3e} cm, r_outer = {r_io[-1]:.3e} cm")


# ---------------------------------------------------------------------------
# Step 7: Write SEDONA param.lua
# ---------------------------------------------------------------------------
def write_param_lua(lua_path, mod_name, t_day, r_inner, t_inner, l_bol,
                    ntasks=8, n_thermal=5e8, n_iter=100):
    """Write SEDONA parameter file for the blackbody illumination test."""
    t_start = t_day * 3600.0 * 24.0  # seconds

    content = f"""-- SEDONA standalone test: SN_A ejecta at {t_day:.0f}d + 7000K inner BB
-- Generated by test_blackbody_illuminated_ejecta_sedona.py

grid_type    = "grid_1D_sphere"
model_file   = "{mod_name}.mod"
hydro_module = "none"

-- SEDONA home directory
sedona_home        = os.getenv('SEDONA_HOME') or '{SEDONA_HOME}'
defaults_file      = sedona_home.."/defaults/sedona_defaults.lua"
data_atomic_file   = sedona_home.."/data/cmfgen_levelcap100.hdf5"
data_fuzzline_file = sedona_home.."/data/kurucz_cd23_fuzz.hdf5"

-- Transport properties
transport_nu_grid  = {{1.0e14, 2.0e16, 0.005, 1}}
transport_radiative_equilibrium  = 0  -- Hold T fixed (no thermal balance)
transport_steady_iterate         = {n_iter}

-- DDMC disabled
transport_use_ddmc               = 0
transport_ddmc_tau_threshold     = 1000.0

-- Output spectrum: optical range (2500-12000 A)
spectrum_nu_grid   = {{2.5e14, 1.2e15, 0.0005, 1}}

-- Time of spectrum calculation
tstep_time_start = {t_start:.1f}

-- Particles
particles_n_emit_radioactive = 1e4
particles_n_emit_thermal     = {n_thermal:.0e}
particles_max_total = 2e9
particles_last_iter_pump     = 20

-- Inner boundary source: 7000K blackbody
core_n_emit      = 0
core_radius      = {r_inner:.8e}
core_temperature = {t_inner:.1f}
core_luminosity  = {l_bol:.8e}
core_fix_luminosity = 1

-- Opacity settings: NLTE for H and He (the key absorbers in outer ejecta)
-- Other elements (C, O, Si, Ca, Fe) treated in LTE for computational speed
opacity_grey_opacity         = 0
opacity_electron_scattering  = 1
opacity_fuzz_expansion       = 0
opacity_line_expansion       = 1
opacity_bound_bound          = 0
opacity_bound_free           = 1
opacity_free_free            = 1
opacity_epsilon              = 1
opacity_use_nlte             = 1
opacity_nlte_dens_min        = 1e-20
opacity_atoms_in_nlte        = {{1, 2}}  -- H, He in NLTE; C, O, Si, Ca, Fe in LTE

-- Line Doppler width: turbulent broadening
line_velocity_width          = 2.0e7  -- 200 km/s

-- Output
output_write_radiation = 1
output_file = "spectrum.h5"
"""
    with open(lua_path, 'w') as f:
        f.write(content)

    print(f"Wrote {lua_path}")
    print(f"  core_temperature = {t_inner} K")
    print(f"  core_luminosity = {l_bol:.3e} erg/s")
    print(f"  NLTE elements: H, He (C, O, Si, Ca, Fe in LTE)")
    print(f"  thermal particles = {n_thermal:.0e}")
    print(f"  iterations = {n_iter}")


# ---------------------------------------------------------------------------
# Step 8: Run SEDONA
# ---------------------------------------------------------------------------
def run_sedona(work_dir, exe, nomp, ntasks):
    """Run SEDONA directly using srun with proper module setup."""
    work_dir = Path(work_dir)

    # Build a shell script that loads modules then runs SEDONA
    # This is needed because SEDONA requires GSL, HDF5, OpenMPI libraries
    shell_script = f"""#!/bin/bash
set -euo pipefail
source /mpcdf/soft/RHEL_9/packages/x86_64/python-waterboa/2024.06/bin/activate /ptmp/sbrennan/conda-envs/astro
unset PYTHONPATH PYTHONHOME

module use /mpcdf/soft/RHEL_9/sub/intel_2024_0/modules/mpi
module use /mpcdf/soft/RHEL_9/modules/compilers
module use /mpcdf/soft/RHEL_9/modules/third-partycompilers

module purge
module load intel/2024.0
module load openmpi/4.1
module load gsl/2.6
module load hdf5-serial/1.14.1

export SEDONA_HOME={SEDONA_HOME}
export LD_LIBRARY_PATH=$(gsl-config --prefix)/lib:$(dirname $(dirname $(which h5cc)))/lib:$PWD:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=1

mpirun --bind-to none -np {ntasks} {exe} param.lua
"""
    script_path = work_dir / 'run_sedona.sh'
    with open(script_path, 'w') as f:
        f.write(shell_script)
    os.chmod(script_path, 0o755)

    # Use srun to allocate resources and run the script
    cmd = ['srun', '--partition=interactive',
           '--cpus-per-task=1',
           '--time=01:30:00',
           f'--ntasks={ntasks}',
           f'--ntasks-per-node={ntasks}',
           'bash', str(script_path)]

    print(f"\nRunning SEDONA in {work_dir}")
    print(f"  Executable: {exe}")
    print(f"  Tasks: {ntasks}")
    print(f"  Script: {script_path}")

    result = subprocess.run(cmd, cwd=str(work_dir),
                           capture_output=False, timeout=7200)
    return result.returncode


# ---------------------------------------------------------------------------
# Step 9: Analyze and plot spectrum
# ---------------------------------------------------------------------------
def analyze_and_plot_spectrum(work_dir):
    """Read SEDONA spectrum.h5 and generate plot + analysis."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import h5py

    work_dir = Path(work_dir)
    spec_path = work_dir / 'spectrum.h5'

    if not spec_path.exists():
        print("ERROR: No spectrum.h5 found")
        return False

    with h5py.File(spec_path, 'r') as f:
        # SEDONA stores spectrum as nu vs L_nu
        # Try different possible keys
        keys = list(f.keys())
        print(f"  HDF5 keys: {keys}")

        # Find the spectrum dataset
        nu = None
        lum = None
        for key in keys:
            ds = f[key]
            if hasattr(ds, 'shape'):
                print(f"    {key}: shape={ds.shape}, dtype={ds.dtype}")
                if ds.ndim == 1 and 'nu' in key.lower():
                    nu = ds[:]
                elif ds.ndim == 1 and ('lum' in key.lower() or 'flux' in key.lower()
                                        or 'spec' in key.lower() or key == 'data'):
                    lum = ds[:]

        # If not found by name, try by shape
        if nu is None or lum is None:
            for key in keys:
                ds = f[key]
                if ds.ndim == 1 and ds.shape[0] > 100:
                    if nu is None:
                        nu = ds[:]
                    elif lum is None:
                        lum = ds[:]

    if nu is None or lum is None:
        print("ERROR: Could not find spectrum data in HDF5")
        return False

    # Convert to wavelength (Angstrom) and luminosity
    wl = C_LIGHT / nu * 1e8  # Hz to Angstrom
    # Sort by wavelength
    idx = np.argsort(wl)
    wl = wl[idx]
    lum = lum[idx]

    # Clip to optical range
    mask = (wl >= 2500) & (wl <= 12000)
    wl = wl[mask]
    lum = lum[mask]

    print(f"\n  Spectrum: {len(wl)} bins, {wl[0]:.0f} - {wl[-1]:.0f} A")
    print(f"  L_nu range: {lum.min():.3e} - {lum.max():.3e} erg/s/Hz")

    # Convert to L_lambda (erg/s/A)
    dnu = np.diff(nu[idx][mask])
    dnu = np.append(dnu, dnu[-1])
    lum_lambda = lum * nu[idx][mask]**2 / C_LIGHT * 1e8  # erg/s/A

    # Smoothed spectrum
    window = max(1, len(wl) // 100)
    smooth = np.convolve(lum_lambda, np.ones(window)/window, mode='same')
    smooth[:window] = smooth[window]
    smooth[-window:] = smooth[-window-1]

    # Plot
    fig, axes = plt.subplots(2, 1, figsize=(14, 8), gridspec_kw={'height_ratios': [3, 2]})

    # Top: spectrum
    ax = axes[0]
    ax.plot(wl, lum_lambda, 'b-', linewidth=0.5, alpha=0.3, label='Raw')
    ax.plot(wl, smooth, 'r-', linewidth=1.5, alpha=0.9, label='Smoothed')

    lines_to_mark = [
        (3933, 'CaII HK'), (4102, 'H-delta'), (4340, 'H-gamma'),
        (4861, 'H-beta'), (5169, 'FeII'), (5876, 'HeI'),
        (6355, 'SiII'), (6563, 'H-alpha'), (6580, 'CII'),
        (7774, 'OI'), (8662, 'CaII IR'),
    ]
    ymax = np.percentile(lum_lambda, 99) * 1.4
    for w, name in lines_to_mark:
        ax.axvline(w, color='gray', alpha=0.4, linestyle=':', linewidth=0.8)
        ax.text(w, ymax*0.97, name, rotation=90, fontsize=7,
                ha='center', va='top', color='gray')

    ax.set_xlabel('Wavelength (Angstrom)')
    ax.set_ylabel('Luminosity (erg/s/A)')
    ax.set_title(f'SEDONA: SN_A ejecta at 400d + 7000K inner BB (NLTE: H,He,C,O,Si,Ca,Fe)\n'
                 f'5e8 thermal particles, {len(wl)} spectrum bins')
    ax.set_xlim(2500, 12000)
    ax.set_ylim(0, ymax)
    ax.legend(loc='upper right', fontsize=9)

    # Bottom: continuum-normalized
    ax = axes[1]
    cont_window = max(window * 4, 50)
    cont = np.convolve(lum_lambda, np.ones(cont_window)/cont_window, mode='same')
    cont[:cont_window] = cont[cont_window]
    cont[-cont_window:] = cont[-cont_window-1]
    cont = np.maximum(cont, 1e-30)
    norm = smooth / cont

    ax.plot(wl, norm, 'r-', linewidth=1.5)
    ax.axhline(1.0, color='k', linewidth=0.5, linestyle='--')

    for w, name in lines_to_mark:
        ax.axvline(w, color='gray', alpha=0.4, linestyle=':', linewidth=0.8)
        ax.text(w, ax.get_ylim()[1]*0.95, name, rotation=90, fontsize=7,
                ha='center', va='top', color='gray')

    ax.set_xlabel('Wavelength (Angstrom)')
    ax.set_ylabel('Normalized flux')
    ax.set_title('Continuum-normalized (absorption = below 1.0)')
    ax.set_xlim(2500, 12000)

    plt.tight_layout()
    outpath = work_dir / 'spectrum.png'
    plt.savefig(outpath, dpi=150)
    print(f"\n  Plot saved: {outpath}")

    # Print line diagnostics
    print(f"\n  Key diagnostic lines (smoothed):")
    for name, w in lines_to_mark:
        idx = np.argmin(np.abs(wl - w))
        local = smooth[max(0,idx-30):min(len(wl),idx+31)]
        cont_local = np.median(local)
        depth = (smooth[idx] - cont_local) / cont_local * 100 if cont_local > 0 else 0
        tag = 'ABSORPTION' if depth < -2 else 'EMISSION' if depth > 2 else 'flat'
        print(f"    {name:15s} ({w:5d} A): depth={depth:+.2f}% [{tag}]")

    return True


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description='Standalone SEDONA test: SN_A ejecta at 400d + 7000K inner BB')
    parser.add_argument('--nomp', type=int, default=8,
                        help='CPUs per task (default: 8)')
    parser.add_argument('--ntasks', type=int, default=8,
                        help='MPI tasks (default: 8)')
    parser.add_argument('--n-thermal', type=float, default=5e8,
                        help='Thermal particles (default: 5e8)')
    parser.add_argument('--n-iter', type=int, default=100,
                        help='Transport iterations (default: 100)')
    parser.add_argument('--output-dir', default=None,
                        help='Output directory')
    parser.add_argument('--no-run', action='store_true',
                        help='Generate inputs only, do not run SEDONA')
    args = parser.parse_args()

    project_dir = Path(__file__).parent.parent
    output_dir = Path(args.output_dir) if args.output_dir else \
                 project_dir / 'tests' / 'bb_illuminated_ejecta_sedona'
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"{'='*60}")
    print(f"Standalone SEDONA test: SN_A ejecta at 400d + 7000K inner BB")
    print(f"{'='*60}")
    print(f"  Output directory: {output_dir}")
    print(f"  SEDONA exe:       {SEDONA_EXE}")
    print(f"  SEDONA_HOME:      {SEDONA_HOME}")
    print(f"  MPI tasks:        {args.ntasks}")
    print(f"  Thermal particles: {args.n_thermal:.0e}")

    # Step 1: Read SN_A profile
    epoch, t_snap = read_sn_a_profile(SN_A_SWD)

    # Step 2: Homologous expansion to 400 days
    r, rho, v, total_mass = homologous_expand(epoch, t_snap, TARGET_AGE_D)

    # Step 3: Detect elements
    elems = detect_elements(epoch['composition'])
    print(f"\nElements present: {', '.join(elems)}")

    # Step 4: Downsample grid
    r_io, rho_io, v_io, comp_io, T_gas_orig, T_rad_orig, n_orig = \
        downsample_grid(r, rho, v, epoch['composition'],
                        epoch['T'], epoch['T_rad'], elems, MAX_CELLS)

    # Step 5: Set temperature profile
    t_gas = set_temperature_profile(r_io, comp_io, T_BB, T_OUTER, T_FLOOR)

    # Step 6: Compute L_bol and R_inner
    l_bol = 1.0e42  # erg/s, typical SN luminosity
    r_inner = r_io[0]  # innermost radius

    # Compute tau_es for reporting
    h_idx = ALL_ELEMENTS.index('H')
    X = comp_io[:, h_idx]
    kappa_es = 0.2 * (1.0 + X)
    r_edges = np.concatenate([[0], r_io])
    dr = np.diff(r_edges)
    tau = np.cumsum(kappa_es * rho_io * dr)

    print(f"\nLuminosity and geometry:")
    print(f"  L_bol: {l_bol:.3e} erg/s = {l_bol/3.828e33:.2e} Lsun")
    print(f"  R_inner: {r_inner:.3e} cm")
    print(f"  R_outer: {r_io[-1]:.3e} cm")
    print(f"  tau_es (total): {tau[-1]:.3f}")
    print(f"  T_inner: {T_BB} K, T_outer: {T_OUTER} K")

    # Step 7: Write .mod file
    mod_name = 'bb_test_sedona'
    mod_path = output_dir / f'{mod_name}.mod'
    write_sedona_mod(mod_path, r_io, v_io, rho_io, t_gas, comp_io, elems,
                     TARGET_AGE_D, l_bol)

    # Step 8: Write param.lua
    lua_path = output_dir / 'param.lua'
    write_param_lua(lua_path, mod_name, TARGET_AGE_D, r_inner, T_BB, l_bol,
                    ntasks=args.ntasks, n_thermal=args.n_thermal,
                    n_iter=args.n_iter)

    if args.no_run:
        print(f"\n--no-run specified: inputs generated, not running SEDONA")
        print(f"  To run manually: cd {output_dir} && "
              f"SEDONA_HOME={SEDONA_HOME} mpirun -np {args.ntasks} "
              f"{SEDONA_EXE} param.lua")
        return

    # Step 9: Run SEDONA
    rc = run_sedona(output_dir, SEDONA_EXE, args.nomp, args.ntasks)
    if rc != 0:
        print(f"\nSEDONA exited with code {rc}")
        sys.exit(rc)

    # Step 10: Analyze and plot
    analyze_and_plot_spectrum(output_dir)

    print(f"\n{'='*60}")
    print(f"Test complete. Outputs in: {output_dir}")
    print(f"{'='*60}")


if __name__ == '__main__':
    main()
