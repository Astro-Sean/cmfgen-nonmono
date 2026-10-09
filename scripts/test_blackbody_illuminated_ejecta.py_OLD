#!/usr/bin/env python3
"""Standalone SuperLite test: 400-day-old SN_A ejecta illuminated by an
internal 7000 K blackbody.

This test bypasses run_model_grid.py entirely.  It:

1. Reads the SN_A_g_T20k_Ni06 STELLA .swd at the latest epoch (~62 d).
2. Homologously expands the profile to 400 days.
3. Sets the inner region (metal-rich core) to 7000 K to act as the
   illuminating blackbody source.
4. Sets the outer region (H/He envelope + He-rich shell) to 2000 K
   (cold ejecta providing line opacity).
5. Writes input.str and input.par directly.
6. Runs SuperLite directly (srun or mpirun, not sbatch).

Usage:
    python3 scripts/test_blackbody_illuminated_ejecta.py [--nomp N] [--n2s E]

Defaults: --nomp 8 --n2s 25  (8 threads, 2^25 = 33M packets)
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
SIGMA_SB = 5.670374419e-5      # Stefan-Boltzmann constant [erg cm^-2 s^-1 K^-4]
C_LIGHT = 2.99792458e10        # speed of light [cm/s]
PI = np.pi
M_SUN = 1.989e33               # solar mass [g]
DAY_S = 86400.0                # day in seconds [s]

SN_A_SWD = "/viper/ptmp1/sbrennan/SN2024pba/SN_A/SN_A_g_T20k_Ni06/stella/SN_A_g_T20k_Ni06.swd"
SUPERLITE_EXE = "/ptmp/sbrennan/superlite/superlite"
SUPERLITE_DATA = "/ptmp/sbrennan/superlite"

# Elements in the 15-element STELLA composition array
ALL_ELEMENTS = ['H', 'He', 'C', 'N', 'O', 'Ne', 'Na', 'Mg', 'Al',
                'Si', 'S', 'Ar', 'Ca', 'Fe', 'Ni']

# Elements present in SN_A (detected from .swd)
# We'll detect them at runtime, but these are the expected ones:
# H, He, C, O, Si, Ca, Fe

# Test parameters
T_BB = 7000.0       # blackbody temperature of inner source [K]
T_OUTER = 2000.0    # temperature of outer H/He material [K]
                         # (cold: negligible emission, pure absorption)
T_FLOOR = 2000.0    # SuperLite zeta-table floor [K]
TARGET_AGE_D = None  # None = use native epoch (no expansion)
MAX_CELLS = 250     # downsample target


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
    """Homologously expand the profile from t_snap to t_target.

    r_new = r * (t_target / t_snap)
    rho_new = rho * (t_snap / t_target)^3
    v_new = r_new / t_target  (homologous velocity, in cm/s)
    """
    scale = t_target_d / t_snap_d
    r_new = epoch['r'] * scale
    rho_new = epoch['rho'] / scale**3
    v_new = r_new / (t_target_d * DAY_S)  # cm/s

    # Mass (approximate, using spherical shells)
    # .swd is ordered outer-to-inner; reverse for inner-to-outer
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
# Step 3: Detect elements and downsample
# ---------------------------------------------------------------------------
def detect_elements(comp):
    """Detect which elements have non-zero abundance."""
    present = []
    for i, elem in enumerate(ALL_ELEMENTS):
        if comp[:, i].max() > 1e-10:
            present.append(elem)
    print(f"\nElements present: {', '.join(present)}")
    return present


def downsample_grid(r, rho, v, comp, T_gas, T_rad, elems, max_cells):
    """Downsample the 500-zone grid to max_cells using simple stride.

    For a standalone test we use uniform stride downsampling.
    The adaptive grid is overkill for this diagnostic.
    """
    n_orig = len(r)
    if n_orig <= max_cells:
        return r, rho, v, comp, T_gas, T_rad, n_orig

    stride = n_orig // max_cells
    # Take every stride-th zone, preserving the first and last
    indices = np.linspace(0, n_orig - 1, max_cells, dtype=int)

    # .swd is outer-to-inner; we need inner-to-outer for input.str
    r_io = r[::-1]      # inner to outer
    rho_io = rho[::-1]
    v_io = v[::-1]
    comp_io = comp[::-1]
    T_io = T_gas[::-1]
    Trad_io = T_rad[::-1]

    r_ds = r_io[indices]
    rho_ds = rho_io[indices]
    v_ds = v_io[indices]
    comp_ds = comp_io[indices]
    T_ds = T_io[indices]
    Trad_ds = Trad_io[indices]

    print(f"\nDownsampled: {n_orig} -> {len(indices)} zones (stride~{stride})")
    return r_ds, rho_ds, v_ds, comp_ds, T_ds, Trad_ds, n_orig


# ---------------------------------------------------------------------------
# Step 4: Set temperature profile (inner 7000K source, outer 5000K absorber)
# ---------------------------------------------------------------------------
def set_temperature_profile(r_io, comp_io, t_bb, t_outer, t_floor):
    """Set gas and radiation temperatures.

    Inner metal-rich zones (Si+Fe dominant, H~0): T = t_bb (7000 K)
      — acts as the illuminating blackbody source
    Outer H/He zones: T = t_outer (5000 K)
      — warm enough for H/He ionization (H is partially ionized,
        He is neutral), cooler than the 7000K source so absorption
        occurs at line wavelengths
    Transition: smooth ramp over the composition interface.

    The temperature gradient between the hot inner source and the
    cooler outer material is what produces absorption lines: the
    outer material sees the 7000K radiation field but is only at
    5000K, so it absorbs photons at wavelengths where it has
    bound-bound transitions (H, He, and metal lines).
    """
    n = len(r_io)

    # Determine which zones are metal-rich (H < 0.01)
    h_idx = ALL_ELEMENTS.index('H')
    is_metal = comp_io[:, h_idx] < 0.01

    # Find the transition zone
    transition_idx = np.where(is_metal)[0]
    if len(transition_idx) == 0:
        # All metal-rich
        t_gas = np.full(n, t_bb)
    else:
        last_metal = transition_idx[-1]
        # Smooth transition over ~20 zones
        ramp_width = min(20, n - last_metal - 1)
        t_gas = np.full(n, t_outer)
        t_gas[:last_metal + 1] = t_bb
        if ramp_width > 0:
            for i in range(last_metal + 1, last_metal + 1 + ramp_width):
                if i < n:
                    frac = (i - last_metal) / ramp_width
                    t_gas[i] = t_bb * (1 - frac) + t_outer * frac

    # Clamp to floor
    t_gas = np.maximum(t_gas, t_floor)
    # Radiation temperature = gas temperature initially
    t_rad = t_gas.copy()

    n_hot = np.sum(t_gas >= t_bb - 1)
    n_outer = np.sum(t_gas <= t_outer + 1)
    n_ramp = n - n_hot - n_outer

    print(f"\nTemperature profile:")
    print(f"  Hot (7000 K) zones:   {n_hot} (inner metal-rich source)")
    print(f"  Outer ({t_outer:.0f} K) zones: {n_outer} (H/He absorber)")
    print(f"  Ramp zones:            {n_ramp}")
    print(f"  T_gas range: {t_gas.min():.0f} - {t_gas.max():.0f} K")

    return t_gas, t_rad


# ---------------------------------------------------------------------------
# Step 5: Compute L_bol and R_phot
# ---------------------------------------------------------------------------
def compute_l_bol_and_r_phot(r_io, rho_io, v_io, t_gas, comp_io, elems):
    """Compute L_bol and R_phot for the test.

    L_bol: In SuperLite, the thermal source is distributed volumetrically
      with spectral shape set by the gas temperature (Planck at T_gas).
      The total luminosity is normalised to in_L_bol.  So the temperature
      sets the *spectrum* and in_L_bol sets the *normalisation*.

      We use a physically reasonable SN-like luminosity (1e42 erg/s)
      rather than the raw 4*pi*R^2*sigma*T^4 of the hot region (which
      would be ~1e45 erg/s and is not what a real SN would produce).

    R_phot: grey optical-depth = 2/3 radius, using electron scattering
      opacity as a floor (consistent with the converter's hybrid approach).
      If the ejecta is optically thin (tau < 2/3), use the inner radius.
    """
    n = len(r_io)

    # L_bol: use a typical SN luminosity.  The 7000K gas temperature
    # ensures the source spectrum peaks in the optical (~4100 A).
    l_bol = 1.0e42

    # R_phot: compute grey optical depth
    # Use electron scattering opacity as floor: kappa_es = 0.2 * (1+X) cm^2/g
    h_idx = ALL_ELEMENTS.index('H')
    X = comp_io[:, h_idx]
    kappa_es = 0.2 * (1.0 + X)  # cm^2/g
    kappa = kappa_es

    # Compute optical depth from outside in
    r_edges = np.concatenate([[0], r_io])
    dr = np.diff(r_edges)
    tau_contrib = kappa * rho_io * dr
    tau_cumul = np.cumsum(tau_contrib[::-1])[::-1]  # cumulative from outside

    # Find tau = 2/3
    r_phot = r_io[-1]  # default to outermost
    for i in range(n - 1, -1, -1):
        if tau_cumul[i] >= 2.0/3.0:
            r_phot = r_io[i]
            break

    tau_max = tau_cumul[0]
    print(f"\nLuminosity and photosphere:")
    print(f"  L_bol (set): {l_bol:.3e} erg/s = {l_bol/3.828e33:.2e} Lsun")
    print(f"  kappa_es range: {kappa.min():.3e} - {kappa.max():.3e} cm^2/g")
    print(f"  tau_max (grey, es): {tau_max:.3f}")
    print(f"  R_phot (tau=2/3): {r_phot:.3e} cm")

    if tau_max < 2.0/3.0:
        print(f"  NOTE: ejecta is optically thin (tau_max={tau_max:.3f} < 2/3)")
        print(f"  Using inner radius as R_phot for nebular EOS")
        r_phot = r_io[0]

    return l_bol, r_phot


# ---------------------------------------------------------------------------
# Step 6: Write input.str
# ---------------------------------------------------------------------------
def write_input_str(path, r_io, v_io, rho_io, t_gas, t_rad, comp_io, elems,
                    t_day, model_name):
    """Write SuperLite input.str file.

    Format (inner-to-outer ordering):
      # 1D spherical - ...
      #      ndim    1 1 nelem nabund
      #         x_left       x_right       vx_left      vx_right          mass      rad_temp          temp             h            he             o            si            ca            fe
      <data rows>
    """
    n = len(r_io)

    # Compute cell edges, masses, velocities
    r_edges = np.concatenate([[0], r_io])
    r_left = r_edges[:-1]
    r_right = r_edges[1:]

    # Velocity at edges (homologous: v = r / t)
    t_s = t_day * DAY_S
    v_left = r_left / t_s
    v_right = r_right / t_s

    # Cell masses
    vol = 4.0/3.0 * PI * (r_right**3 - r_left**3)
    mass = vol * rho_io

    # Get abundance columns for present elements
    elem_indices = [ALL_ELEMENTS.index(e) for e in elems]
    abund_cols = comp_io[:, elem_indices]

    # Renormalize abundances per cell
    row_sums = abund_cols.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1.0
    abund_cols = abund_cols / row_sums

    nabund = len(elems)

    with open(path, 'w') as f:
        f.write(f"# 1D spherical - SN_A ejecta at {t_day:.1f}d, "
                f"inner 7000K blackbody illumination test\n")
        # ncol = 7 structure columns (x_left, x_right, vx_left, vx_right,
        #         mass, rad_temp, temp) + nabund abundance columns
        ncol = 7 + nabund
        f.write(f"#      {n}          1 1 {ncol} {nabund}\n")
        # Header: element symbols
        header = "#         x_left       x_right       vx_left      vx_right" \
                 "          mass      rad_temp          temp"
        for e in elems:
            header += f"            {e.lower()}"
        f.write(header + "\n")

        for i in range(n):
            parts = [f"  {r_left[i]:.6e}", f"  {r_right[i]:.6e}",
                     f"  {v_left[i]:.6e}", f"  {v_right[i]:.6e}",
                     f"  {mass[i]:.6e}",
                     f"  {t_rad[i]:.6e}", f"  {t_gas[i]:.6e}"]
            for j in range(nabund):
                parts.append(f"  {abund_cols[i, j]:.6e}")
            f.write("  ".join(parts) + "  \n")

    print(f"\nWrote {path}: {n} cells, {nabund} abundance columns "
          f"({', '.join(elems)})")


# ---------------------------------------------------------------------------
# Step 7: Write input.par
# ---------------------------------------------------------------------------
def write_input_par(path, model_name, t_day, ndim, l_bol, r_phot, nomp,
                    n2s=25, grp_ng=6000, flx_ndim=3500,
                    tauddmc='3d0', taulump='3d0',
                    nlte=True, niter=5, nebular=True):
    """Write SuperLite input.par file."""
    nlte_str = 't' if nlte else 'f'
    nebular_str = 't' if nebular else 'f'

    par = f"""&inputpars

in_name = '{model_name}'
in_comment = 'Standalone test: SN_A ejecta at {t_day:.0f}d, inner 7000K BB'

!-- misc
 in_nomp = {nomp}
 in_io_grabstdout = f
 in_io_opacdump = 'one'
 in_io_profdump = t

!-- grid settings
 in_grd_igeom = 11
 in_ndim = {ndim}, 1, 1

!-- group settings (aligned with 5 A flux bins: edges at 0.001, 5.001, ..., 30000.001)
 in_grp_ng = {grp_ng}
 in_grp_wlmin = 1d-11
 in_grp_wlmax = 30d-5

!-- flux settings (3000-20000 Angstrom, 5 A bins aligned with groups)
 in_flx_ndim = {flx_ndim}, 1, 1
 in_flx_wlmin = 3d-5
 in_flx_wlmax = 2d-4

!-- particle numbers
 in_src_n2s = {n2s}
 in_L_bol = {l_bol:.4e}

!-- ddmc
 in_trn_tauddmc = {tauddmc}
 in_taulump = {taulump}

!-- EOS / ionization
 in_nebular = {nebular_str}
 in_R_phot = {r_phot:.4e}

!-- NLTE
 in_nlte = {nlte_str}
 in_nlte_nelem = 1
 in_niter = {niter}

!-- opacity settings (all physical opacities enabled)
 in_nobbopac = f
 in_nobfopac = f
 in_noffopac = f
 in_nothmson = f

/
"""
    with open(path, 'w') as f:
        f.write(par)
    print(f"Wrote {path} (ndim={ndim}, L_bol={l_bol:.3e} erg/s, "
          f"R_phot={r_phot:.3e} cm)")


# ---------------------------------------------------------------------------
# Step 8: Run SuperLite
# ---------------------------------------------------------------------------
def run_superlite(work_dir, exe, data_dir, nomp):
    """Run SuperLite directly using srun (interactive partition).

    Copies required atomic data files, then runs the executable.
    SuperLite reads:
      - data.ion, data.zeta, data.bf_verner, data.ff_sutherland  (from CWD)
      - Atoms/data.atom.<sym><ion>                                (from CWD/Atoms/)
    """
    work_dir = Path(work_dir)

    # Use the known working reference directory for atomic data
    ref_dir = Path('/viper/ptmp1/sbrennan/SN2024pba/SN_B_grid/M1.0_V10000/superlite/30.0days')

    # Copy top-level data files
    data_files = ['data.bf_verner', 'data.ff_sutherland', 'data.ion',
                  'data.zeta']
    for df in data_files:
        src = ref_dir / df
        dst = work_dir / df
        if src.exists() and not dst.exists():
            shutil.copy2(src, dst)
            print(f"  Copied {df}")

    # Copy Atoms directory
    data_subdir = work_dir / 'Atoms'
    data_subdir.mkdir(exist_ok=True)
    ref_atoms = ref_dir / 'Atoms'
    if ref_atoms.exists():
        for af in ref_atoms.iterdir():
            dst = data_subdir / af.name
            if not dst.exists():
                shutil.copy2(af, dst)
        print(f"  Copied {len(list(data_subdir.iterdir()))} atom files from {ref_atoms}")

    print(f"\nRunning SuperLite in {work_dir}")
    print(f"  Executable: {exe}")
    print(f"  Threads: {nomp}")

    # Set OMP_NUM_THREADS
    env = os.environ.copy()
    env['OMP_NUM_THREADS'] = str(nomp)

    # Run with srun (interactive partition)
    cmd = ['srun', '--partition=interactive', '--cpus-per-task={}'.format(nomp),
           '--time=01:30:00', '--ntasks=1', exe]

    print(f"  Command: {' '.join(cmd)}")
    print(f"  OMP_NUM_THREADS={nomp}")

    result = subprocess.run(cmd, cwd=str(work_dir), env=env,
                           capture_output=False, timeout=7200)
    return result.returncode


# ---------------------------------------------------------------------------
# Step 9: Analyze spectrum
# ---------------------------------------------------------------------------
def analyze_spectrum(work_dir):
    """Read and analyze the SuperLite output spectrum."""
    work_dir = Path(work_dir)

    flx_path = work_dir / 'output.flx_luminos'
    grid_path = work_dir / 'output.flx_grid'
    energy_path = work_dir / 'output.tot_energy'

    if not flx_path.exists():
        print("ERROR: No output.flx_luminos found — SuperLite may have failed")
        return False

    # Read spectrum
    with open(flx_path) as f:
        lines = f.readlines()
    flux = np.array([float(x) for x in lines[1].split()])

    with open(grid_path) as f:
        wl_lines = f.readlines()
    wl = np.array([float(x) for x in wl_lines[1].split()]) * 1e8  # cm to Angstrom

    n = min(len(flux), len(wl))
    flux = flux[:n]
    wl = wl[:n]

    print(f"\n{'='*60}")
    print(f"Spectrum analysis: {n} bins, {wl[0]:.0f}-{wl[-1]:.0f} A")
    print(f"{'='*60}")
    print(f"  Flux range: {flux.min():.3e} - {flux.max():.3e} erg/s")
    print(f"  Total flux: {flux.sum():.3e} erg/s")

    # Energy conservation
    if energy_path.exists():
        with open(energy_path) as f:
            elines = f.readlines()
        # First data line: eout, evelo, sflux, sthermal
        evals = [float(x) for x in elines[2].split()]
        print(f"\n  Energy conservation:")
        print(f"    sthermal (input):  {evals[3]:.3e} erg/s")
        print(f"    eout (emergent):   {evals[0]:.3e} erg/s")
        print(f"    Ratio eout/stherm: {evals[0]/evals[3]:.4f}")

    # Spectrum shape
    print(f"\n  Spectrum shape (1000 A bins):")
    for w0 in range(3000, 20000, 1000):
        mask = (wl >= w0) & (wl < w0 + 1000)
        if mask.sum() > 0:
            f_mean = flux[mask].mean()
            f_std = flux[mask].std()
            cv = f_std / f_mean if f_mean > 0 else 0
            print(f"    {w0:5d}-{w0+1000:5d} A: mean={f_mean:.3e} "
                  f"CV={cv:.4f} {'(lines detected)' if cv > 0.15 else ''}")

    # Look for absorption features at key wavelengths
    print(f"\n  Key diagnostic lines:")
    diagnostics = [
        ('Ca II HK', 3933), ('H delta', 4102), ('H gamma', 4340),
        ('H beta', 4861), ('Fe II', 5169), ('Si II 6355', 6355),
        ('H alpha', 6563), ('Ca II IR', 8662), ('O I 7774', 7774),
        ('He I 5876', 5876), ('C II 6580', 6580),
    ]
    for name, w in diagnostics:
        idx = np.argmin(np.abs(wl - w))
        # Local continuum: average over +/- 100 A, excluding +/- 30 A
        mask = ((wl >= w - 100) & (wl < w - 30)) | \
               ((wl > w + 30) & (wl <= w + 100))
        if mask.sum() > 0 and flux[idx] > 0:
            cont = flux[mask].mean()
            depth = (cont - flux[idx]) / cont if cont > 0 else 0
            tag = 'ABSORPTION' if depth > 0.02 else \
                  'EMISSION' if depth < -0.02 else 'flat'
            print(f"    {name:15s} ({w:5d} A): flux={flux[idx]:.3e} "
                  f"cont={cont:.3e} depth={depth:+.3f} [{tag}]")
        else:
            print(f"    {name:15s} ({w:5d} A): flux={flux[idx]:.3e} (no continuum)")

    return True


# ---------------------------------------------------------------------------
# Step 11: Plot spectrum
# ---------------------------------------------------------------------------
def plot_spectrum(work_dir):
    """Generate a publication-quality spectrum plot."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    work_dir = Path(work_dir)

    flx_path = work_dir / 'output.flx_luminos'
    grid_path = work_dir / 'output.flx_grid'
    energy_path = work_dir / 'output.tot_energy'

    if not flx_path.exists():
        print("ERROR: No output.flx_luminos found")
        return

    # Read spectrum (last iteration)
    with open(flx_path) as f:
        lines = f.readlines()
    flux = np.array([float(x) for x in lines[-1].split()])

    with open(grid_path) as f:
        wl_lines = f.readlines()
    wl = np.array([float(x) for x in wl_lines[1].split()]) * 1e8  # cm to A

    n = min(len(flux), len(wl))
    flux = flux[:n]
    wl = wl[:n]

    # Clip edge artifacts (first/last 3 bins accumulate boundary photons)
    flux_clip = flux.copy()
    flux_clip[:3] = np.median(flux[3:50])
    flux_clip[-3:] = np.median(flux[-50:-3])

    # Smoothed spectrum (50-bin ~ 240 A window for continuum)
    window = 50
    smooth = np.convolve(flux_clip, np.ones(window)/window, mode='same')
    smooth[:window] = smooth[window]
    smooth[-window:] = smooth[-window-1]

    # Energy conservation
    e_in = e_out = 0
    if energy_path.exists():
        with open(energy_path) as f:
            elines = f.readlines()
        evals = [float(x) for x in elines[2].split()]
        e_in = evals[3]
        e_out = evals[0]

    # 7000K Planck for comparison
    h_planck = 6.626e-27
    c_light = 2.998e10
    k_boltz = 1.381e-16
    T_bb = 7000.0
    wl_cm = wl * 1e-8
    planck = (2*h_planck*c_light**2 / wl_cm**5) / \
             (np.exp(h_planck*c_light/(wl_cm*k_boltz*T_bb)) - 1)
    planck_norm = planck / planck[wl > 4000][:100].mean() * \
                  smooth[(wl > 4000) & (wl < 6000)].mean()

    fig, axes = plt.subplots(2, 1, figsize=(14, 8), gridspec_kw={'height_ratios': [3, 2]})

    # Top: Raw + smoothed + Planck
    ax = axes[0]
    ax.plot(wl, flux_clip, 'b-', linewidth=0.3, alpha=0.2, label='Raw (binned)')
    ax.plot(wl, smooth, 'r-', linewidth=1.5, alpha=0.9, label='Smoothed (50-bin)')
    ax.plot(wl, planck_norm, 'k--', linewidth=1, alpha=0.5, label='7000K Planck (normalized)')

    # Mark key lines
    lines_to_mark = [
        (3933, 'CaII HK'), (4102, 'H-delta'), (4340, 'H-gamma'),
        (4861, 'H-beta'), (5169, 'FeII'), (5876, 'HeI'),
        (6355, 'SiII'), (6563, 'H-alpha'), (6580, 'CII'),
        (7774, 'OI'), (8662, 'CaII IR'),
    ]
    ymax = np.percentile(flux_clip[10:-10], 99) * 1.4
    for w, name in lines_to_mark:
        ax.axvline(w, color='gray', alpha=0.4, linestyle=':', linewidth=0.8)
        ax.text(w, ymax*0.97, name, rotation=90, fontsize=7,
                ha='center', va='top', color='gray')

    ax.set_xlabel('Wavelength (Angstrom)')
    ax.set_ylabel('Luminosity (erg/s/bin)')
    title = f'SN_A ejecta + 7000K inner BB (native epoch)\n'
    if e_in > 0:
        title += f'L_in={e_in:.2e}, L_out={e_out:.2e} erg/s, '
    title += f'33M packets, 6000 groups'
    ax.set_title(title)
    ax.set_xlim(3000, 20000)
    ax.set_ylim(0, ymax)
    ax.legend(loc='upper right', fontsize=9)

    # Bottom: Continuum-normalized (smoothed)
    ax = axes[1]
    # Use larger window for continuum
    cont_window = 200
    cont = np.convolve(flux_clip, np.ones(cont_window)/cont_window, mode='same')
    cont[:cont_window] = cont[cont_window]
    cont[-cont_window:] = cont[-cont_window-1]
    # Avoid division by zero
    cont = np.maximum(cont, 1e-30)
    norm = smooth / cont

    ax.plot(wl, norm, 'r-', linewidth=1.5)
    ax.axhline(1.0, color='k', linewidth=0.5, linestyle='--')

    for w, name in lines_to_mark:
        ax.axvline(w, color='gray', alpha=0.4, linestyle=':', linewidth=0.8)
        ax.text(w, ax.get_ylim()[1]*0.95, name, rotation=90, fontsize=7,
                ha='center', va='top', color='gray')

    ax.set_xlabel('Wavelength (Angstrom)')
    ax.set_ylabel('Normalized flux (F / F_continuum)')
    ax.set_title('Continuum-normalized spectrum (absorption = below 1.0)')
    ax.set_xlim(3000, 20000)

    plt.tight_layout()
    outpath = work_dir / 'spectrum.png'
    plt.savefig(outpath, dpi=150)
    print(f"\nPlot saved: {outpath}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description='Standalone SuperLite test: SN_A ejecta + 7000K inner BB')
    parser.add_argument('--nomp', type=int, default=8,
                        help='OpenMP threads (default: 8)')
    parser.add_argument('--n2s', type=int, default=25,
                        help='2^n source particles (default: 25 = 33M)')
    parser.add_argument('--grp-ng', type=int, default=6000,
                        help='Wavelength groups (default: 6000)')
    parser.add_argument('--flx-ndim', type=int, default=3400,
                        help='Flux output bins (default: 3400, aligned with groups)')
    parser.add_argument('--output-dir', default=None,
                        help='Output directory (default: tests/bb_illuminated_ejecta)')
    parser.add_argument('--no-run', action='store_true',
                        help='Generate inputs only, do not run SuperLite')
    args = parser.parse_args()

    project_dir = Path(__file__).parent.parent
    output_dir = Path(args.output_dir) if args.output_dir else \
                 project_dir / 'tests' / 'bb_illuminated_ejecta'
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"{'='*60}")
    print(f"Standalone SuperLite test: SN_A ejecta + 7000K inner BB")
    print(f"{'='*60}")
    print(f"  Output directory: {output_dir}")
    print(f"  SuperLite exe:    {SUPERLITE_EXE}")
    print(f"  Threads:          {args.nomp}")
    print(f"  Particles:        2^{args.n2s} = {2**args.n2s:,}")

    # Step 1: Read SN_A profile
    epoch, t_snap = read_sn_a_profile(SN_A_SWD)

    # Step 2: Use native epoch or homologously expand
    if TARGET_AGE_D is not None and TARGET_AGE_D > t_snap:
        r, rho, v, total_mass = homologous_expand(epoch, t_snap, TARGET_AGE_D)
        t_day = TARGET_AGE_D
    else:
        # Use native epoch — much higher density, optically thick
        r = epoch['r']
        rho = epoch['rho']
        v_km = epoch['v']  # km/s
        v = v_km * 1e5     # convert to cm/s
        t_day = t_snap
        # Compute mass for reporting
        r_io_tmp = r[::-1]
        rho_io_tmp = rho[::-1]
        r_edges_tmp = np.concatenate([[0], r_io_tmp])
        vol_tmp = 4.0/3.0 * PI * (r_edges_tmp[1:]**3 - r_edges_tmp[:-1]**3)
        total_mass = (vol_tmp * rho_io_tmp).sum()
        print(f"\nUsing native epoch (t={t_snap:.1f}d, no expansion):")
        print(f"  R: {r.min():.3e} - {r.max():.3e} cm")
        print(f"  rho: {rho.min():.3e} - {rho.max():.3e} g/cm^3")
        print(f"  Total mass: {total_mass/M_SUN:.4f} Msun")

    # Step 3: Detect elements
    elems = detect_elements(epoch['composition'])

    # Step 4: Downsample grid
    r_io, rho_io, v_io, comp_io, T_gas_orig, T_rad_orig, n_orig = \
        downsample_grid(r, rho, v, epoch['composition'],
                        epoch['T'], epoch['T_rad'], elems, MAX_CELLS)

    # Step 5: Set temperature profile (inner 7000K source, outer 5000K absorber)
    t_gas, t_rad = set_temperature_profile(r_io, comp_io, T_BB, T_OUTER, T_FLOOR)

    # Step 6: Compute L_bol and R_phot
    l_bol, r_phot = compute_l_bol_and_r_phot(r_io, rho_io, v_io, t_gas,
                                              comp_io, elems)

    # Step 7: Write input.str
    str_path = output_dir / 'input.str'
    write_input_str(str_path, r_io, v_io, rho_io, t_gas, t_rad, comp_io, elems,
                    t_day, 'bb_test')

    # Step 8: Write input.par
    par_path = output_dir / 'input.par'
    write_input_par(par_path, 'bb_test', t_day, len(r_io),
                    l_bol, r_phot, args.nomp,
                    n2s=args.n2s, grp_ng=args.grp_ng,
                    flx_ndim=args.flx_ndim)

    if args.no_run:
        print(f"\n--no-run specified: inputs generated, not running SuperLite")
        print(f"  To run manually: cd {output_dir} && OMP_NUM_THREADS={args.nomp} "
              f"{SUPERLITE_EXE}")
        return

    # Step 9: Run SuperLite
    rc = run_superlite(output_dir, SUPERLITE_EXE, SUPERLITE_DATA, args.nomp)
    if rc != 0:
        print(f"\nSuperLite exited with code {rc}")
        sys.exit(rc)

    # Step 10: Analyze spectrum
    analyze_spectrum(output_dir)

    # Step 11: Plot spectrum
    plot_spectrum(output_dir)

    print(f"\n{'='*60}")
    print(f"Test complete. Outputs in: {output_dir}")
    print(f"{'='*60}")


if __name__ == '__main__':
    main()
