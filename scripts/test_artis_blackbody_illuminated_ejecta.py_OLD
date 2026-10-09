#!/usr/bin/env python3
"""Standalone ARTIS test: 400-day-old SN_A ejecta illuminated by an
internal 7000 K blackbody source.

This test bypasses run_model_grid.py entirely.  It:

1. Reads the SN_A_g_T20k_Ni06 STELLA .swd at the latest epoch (~62 d).
2. Homologously expands the profile to 400 days.
3. Sets the inner region (metal-rich core) to 7000 K with a thermal
   energy density (q column) that acts as the illuminating blackbody
   source via ARTIS's USE_MODEL_INITIAL_ENERGY mechanism.
4. Sets the outer region (H/He envelope + He-rich shell) to 2000 K
   (cold ejecta providing line opacity).
5. Writes ARTIS input files: model.txt, abundances.txt, input.txt,
   compositiondata.txt.
6. Runs ARTIS directly (mpirun, not sbatch).

The ARTIS "q" column provides initial thermal energy per gram at t_model.
With USE_MODEL_INITIAL_ENERGY=true and INITIAL_PACKETS_ON=true (both set
in the classic options preset), ARTIS deposits this energy as k-packets
at tmin.  These k-packets thermalize and emerge as r-packets with a
Planck spectrum at the local electron temperature T_e.  This is exactly
the "blackbody within the ejecta" scenario we need for the SN_B test.

Usage:
    python3 scripts/test_artis_blackbody_illuminated_ejecta.py [--nrank N] [--mpkts M]

Defaults: --nrank 8 --mpkts 100000
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
KB = 1.380649e-16              # Boltzmann constant [erg/K]
H_PLANCK = 6.62607015e-27      # Planck constant [erg*s]
EV_TO_ERG = 1.602176634e-12    # eV to erg

SN_A_SWD = "/viper/ptmp1/sbrennan/SN2024pba/SN_A/SN_A_g_T20k_Ni06/stella/SN_A_g_T20k_Ni06.swd"
ARTIS_REPO = "/viper/ptmp1/sbrennan/artis_repo"
ARTIS_EXE = f"{ARTIS_REPO}/sn3d"
ARTIS_ATOMIC = "/viper/ptmp1/sbrennan/artisatomic/artis_files"

# Elements in the 15-element STELLA composition array
ALL_ELEMENTS = ['H', 'He', 'C', 'N', 'O', 'Ne', 'Na', 'Mg', 'Al',
                'Si', 'S', 'Ar', 'Ca', 'Fe', 'Ni']

# ARTIS uses Z=1..30 for abundances.txt.  Map our elements to Z.
ELEMENT_Z = {'H': 1, 'He': 2, 'C': 6, 'N': 7, 'O': 8, 'Ne': 10,
             'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'S': 16,
             'Ar': 18, 'Ca': 20, 'Fe': 26, 'Ni': 28}

# Test parameters
T_BB = 12000.0      # blackbody temperature of inner source [K]
T_OUTER = 4000.0    # temperature of outer H/He material [K]
TARGET_AGE_D = 400  # expand to 400 days
MAX_CELLS = 400     # downsample target (ARTIS handles large grids, but keep it manageable)


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


def downsample_grid(r, rho, v, comp, T_gas, T_rad, max_cells):
    """Downsample the 500-zone grid to max_cells using uniform stride.

    .swd is outer-to-inner; ARTIS model.txt is inner-to-outer.
    """
    n_orig = len(r)
    if n_orig <= max_cells:
        indices = np.arange(n_orig)
    else:
        indices = np.linspace(0, n_orig - 1, max_cells, dtype=int)

    # .swd is outer-to-inner; we need inner-to-outer for ARTIS
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

    print(f"\nDownsampled: {n_orig} -> {len(indices)} zones")
    return r_ds, rho_ds, v_ds, comp_ds, T_ds, Trad_ds, n_orig


# ---------------------------------------------------------------------------
# Step 4: Set temperature profile and compute initial energy (q)
# ---------------------------------------------------------------------------
def set_temperature_and_energy(r_io, comp_io, rho_io, t_bb, t_outer, t_model_d):
    """Set gas temperature and initial thermal energy density.

    Inner metal-rich zones (Si+Fe dominant, H~0): T = t_bb (7000 K)
      — acts as the illuminating blackbody source
    Outer H/He zones: T = t_outer (2000 K)
      — cold ejecta providing line opacity

    The q column (initial energy per gram at t_model) is set for the
    inner hot zones so that ARTIS deposits thermal k-packets there.
    The total energy is chosen to give a luminosity of ~1e42 erg/s
    (typical SN-like luminosity) over the simulation timescale.

    q has units of erg/g at t_model.  ARTIS scales it to tmin:
      q_tmin = q_tmodel * (t_model / tmin)
    (see grid.cc line 763: set_initenergyq(valuein * t_model / tmin))

    The k-packets are emitted with a Planck spectrum at the local T_e,
    so setting T_e = 7000 K in the inner zones gives a 7000 K blackbody.
    """
    n = len(r_io)

    # Determine which zones are metal-rich (H < 0.01)
    h_idx = ALL_ELEMENTS.index('H')
    is_metal = comp_io[:, h_idx] < 0.01

    # Find the transition zone
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

    # Compute initial thermal energy density (q) for the hot inner zones
    # Target luminosity: L = 1e42 erg/s
    # Simulation duration: ~t_model to t_model + ~10 days
    # Total energy: E = L * delta_t ~ 1e42 * 10*86400 ~ 8.6e47 erg
    # But q is energy per gram, so q = E / M_hot
    target_l_bol = 1.0e42  # erg/s
    sim_duration_s = 10.0 * DAY_S  # 10 days
    total_energy = target_l_bol * sim_duration_s  # erg

    # Mass of hot zones
    r_edges = np.concatenate([[0], r_io])
    vol = 4.0/3.0 * PI * (r_edges[1:]**3 - r_edges[:-1]**3)
    mass = vol * rho_io
    mass_hot = mass[is_metal].sum()

    if mass_hot > 0:
        q_hot = total_energy / mass_hot  # erg/g at t_model
    else:
        q_hot = 0.0

    # q array: only hot zones get energy, cold zones get zero
    q = np.zeros(n)
    q[is_metal] = q_hot

    n_hot = np.sum(is_metal)
    n_outer = np.sum(~is_metal)

    print(f"\nTemperature and energy profile:")
    print(f"  Hot ({t_bb:.0f} K) zones:   {n_hot} (inner metal-rich source)")
    print(f"  Outer ({t_outer:.0f} K) zones: {n_outer} (H/He absorber)")
    print(f"  T_gas range: {t_gas.min():.0f} - {t_gas.max():.0f} K")
    print(f"  Mass of hot zones: {mass_hot/M_SUN:.4f} Msun")
    print(f"  q_hot: {q_hot:.3e} erg/g")
    print(f"  Total initial energy: {total_energy:.3e} erg")
    print(f"  Target L_bol: {target_l_bol:.3e} erg/s")

    return t_gas, q, target_l_bol


# ---------------------------------------------------------------------------
# Step 5: Write ARTIS model.txt
# ---------------------------------------------------------------------------
def write_model_txt(path, r_io, v_io, rho_io, t_gas, q, t_model_d):
    """Write ARTIS model.txt file (1D spherical format).

    Format (1D spherical, inner-to-outer):
      Line 1: ncells
      Line 2: t_model (days)
      Then per cell (one line per cell with header):
      #inputcellid velocity_outer_kmps logrho X_Fegroup X_Ni56 X_Co56 X_Fe52 X_Cr48 q

    The velocity is the outer boundary velocity of each cell.
    logrho is log10(rho [g/cm^3]) at t_model.
    q is the initial thermal energy per gram at t_model [erg/g].

    For this test we set all radioactive abundances to zero (no Ni56,
    no Co56, etc.) since the energy source is the thermal q column,
    not radioactive decay.
    """
    n = len(r_io)

    with open(path, 'w') as f:
        f.write(f"{n}\n")
        f.write(f"{t_model_d:.6f}\n")
        # Header specifying columns
        f.write("#inputcellid velocity_outer logrho X_Fegroup X_Ni56 X_Co56 X_Fe52 X_Cr48 q\n")

        for i in range(n):
            cellid = i + 1
            v_outer = v_io[i] / 1e5  # km/s
            log_rho = np.log10(rho_io[i]) if rho_io[i] > 0 else -99.0
            # No radioactive isotopes except a negligible trace of Ni56
            # to satisfy ARTIS's requirement that etot_tinf > 0
            # (the decay energy from this trace is << initial q energy)
            x_fegrp = 0.0
            x_ni56 = 1e-10
            x_co56 = 0.0
            x_fe52 = 0.0
            x_cr48 = 0.0
            q_val = q[i]
            f.write(f"   {cellid} {v_outer:.6e} {log_rho:.6f} "
                    f"{x_fegrp:.6e} {x_ni56:.6e} {x_co56:.6e} "
                    f"{x_fe52:.6e} {x_cr48:.6e} {q_val:.6e}\n")

    print(f"\nWrote {path}: {n} cells, t_model={t_model_d:.1f} d")


# ---------------------------------------------------------------------------
# Step 6: Write ARTIS abundances.txt
# ---------------------------------------------------------------------------
def write_abundances_txt(path, comp_io, elems):
    """Write ARTIS abundances.txt file.

    Format: one line per cell, starting with cell number (1-indexed),
    followed by mass fractions for Z=1 through Z=30 (30 values).
    Cells with zero density should be skipped.
    """
    n = len(comp_io)
    max_z = 30

    with open(path, 'w') as f:
        for i in range(n):
            cellid = i + 1
            # Build abundance array for Z=1..30
            abund = np.zeros(max_z)
            for elem in elems:
                z = ELEMENT_Z.get(elem, 0)
                if z > 0 and z <= max_z:
                    # comp_io has 15 columns indexed by ALL_ELEMENTS
                    elem_idx = ALL_ELEMENTS.index(elem)
                    abund[z - 1] = comp_io[i, elem_idx]

            # Renormalize to sum to 1 (if any non-zero)
            total = abund.sum()
            if total > 0:
                abund = abund / total

            # Write: cellid followed by 30 mass fractions
            vals = ' '.join(f'{a:.8e}' for a in abund)
            f.write(f"{cellid} {vals}\n")

    print(f"Wrote {path}: {n} cells, {len(elems)} elements, Z=1..{max_z}")


# ---------------------------------------------------------------------------
# Step 7: Write ARTIS input.txt
# ---------------------------------------------------------------------------
def write_input_txt(path, t_model_d, t_max_d, n_timesteps, n_lte_ts, mpkts):
    """Write ARTIS input.txt file.

    Key parameters:
    - tmin = t_model (start at the model time)
    - tmax = t_model + t_max_d (end time)
    - num_lte_timesteps: first N timesteps in LTE (for stability)
    - optical_depth_is_thick: cells with tau > this use grey opacity
    - num_grey_timesteps: number of timesteps to use grey opacity

    The simulation starts at tmin = t_model and runs to tmax.
    Initial packets (from the q column) are deposited at tmin.
    """
    tmin = t_model_d
    tmax = t_model_d + t_max_d

    lines = [
        "1281360349               # pre_zseed",
        f"{n_timesteps}                       # ntimesteps",
        f"000 {n_timesteps:03d}                  # timestep_start timestep_finish",
        f"{tmin:.6f} {tmax:.6f}        # tmin_days tmax_days",
        "1.33  1.330000001        # nusyn_min_mev nusyn_max_mev (UNUSED)",
        "80                       # nsyn_time (UNUSED)",
        "3. 0.037                 # start and end times for synthesis (UNUSED)",
        "1                        # model_type (UNUSED, autodetected)",
        "4                        # compute r-light curve (UNUSED)",
        "1                        # n_out_it (UNUSED)",
        "1.0                      # CLIGHT/CLIGHT (UNUSED)",
        "-1                       # gamma_kappagrey (UNUSED)",
        "0 0 1                    # syn_dir (UNUSED)",
        "4                        # opacity_case (UNUSED)",
        "1.0e-10                  # rho_crit_para (UNUSED)",
        "-1                       # debug_packet (UNUSED)",
        "0                        # simulation_continued_from_saved",
        "1e-6                     # rfcut (UNUSED)",
        f"{n_lte_ts}                        # num_lte_timesteps",
        "8.0 999                  # optical_depth_is_thick num_grey_timesteps",
        "-1                       # max_bf_continua (UNUSED)",
        "1                        # nprocs_exspec",
        "1                        # do_emission_res (UNUSED)",
        "0.001 1000               # kpktdiffusion (UNUSED)",
    ]

    with open(path, 'w') as f:
        f.write('\n'.join(lines) + '\n')

    print(f"Wrote {path}: tmin={tmin:.1f}d, tmax={tmax:.1f}d, "
          f"{n_timesteps} timesteps, {n_lte_ts} LTE")


# ---------------------------------------------------------------------------
# Step 8: Write vpkt.txt (virtual packet settings)
# ---------------------------------------------------------------------------
def write_vpkt_txt(path):
    """Write ARTIS vpkt.txt file (virtual packet configuration).

    Virtual packets improve signal-to-noise in 1D but are optional.
    We disable them for this test to keep it simple.
    """
    # Minimal vpkt.txt: no virtual packets
    lines = [
        "0",           # number of virtual packet directions
        "0 10 30",     # UNUSED when n_vpkt_dirs=0
        "0 100",       # UNUSED
        "10",          # UNUSED
        "0",           # UNUSED
        "1",           # UNUSED
        "11.5 21.5",   # UNUSED
    ]
    with open(path, 'w') as f:
        f.write('\n'.join(lines) + '\n')
    print(f"Wrote {path} (virtual packets disabled)")


# ---------------------------------------------------------------------------
# Step 9: Copy atomic data and run ARTIS
# ---------------------------------------------------------------------------
def setup_and_run_artis(work_dir, n_rank, mpkts, t_model_d, t_max_d):
    """Copy atomic data files and run ARTIS.

    ARTIS reads from the current working directory:
    - adata.txt, transitiondata.txt, compositiondata.txt, phixsdata_v2.txt
    - model.txt, abundances.txt, input.txt, vpkt.txt
    - data/ directory (gamma decay files, etc.)
    """
    work_dir = Path(work_dir)

    # Copy atomic data files
    atomic_files = ['adata.txt', 'transitiondata.txt', 'compositiondata.txt',
                    'phixsdata_v2.txt']
    for af in atomic_files:
        src = Path(ARTIS_ATOMIC) / af
        dst = work_dir / af
        if src.exists() and not dst.exists():
            shutil.copy2(src, dst)
            print(f"  Copied {af} ({src.stat().st_size / 1e6:.1f} MB)")

    # Copy data/ directory (gamma decay spectra, etc.)
    data_src = Path(ARTIS_REPO) / "data"
    data_dst = work_dir / "data"
    if data_src.exists() and not data_dst.exists():
        shutil.copytree(data_src, data_dst, dirs_exist_ok=True)
        print(f"  Copied data/ directory")

    # Copy decay data files (betaminusdecays.txt, etc.)
    for df in ['betaminusdecays.txt', 'alphadecays.txt', 'fissiondecays.txt']:
        src = Path(ARTIS_REPO) / "data" / df
        dst = work_dir / df
        if src.exists() and not dst.exists():
            shutil.copy2(src, dst)

    # Set up environment
    env = os.environ.copy()
    env['OMPI_CXX'] = 'g++'
    env['LD_LIBRARY_PATH'] = (
        '/mpcdf/soft/RHEL_9/packages/znver4/gsl/gcc_15-15.1.0/2.7.1/lib:'
        + env.get('LD_LIBRARY_PATH', '')
    )

    # Build mpirun command
    cmd = ['mpirun', '--bind-to', 'none', '-np', str(n_rank), ARTIS_EXE]

    print(f"\nRunning ARTIS in {work_dir}")
    print(f"  Executable: {ARTIS_EXE}")
    print(f"  MPI ranks: {n_rank}")
    print(f"  Command: {' '.join(cmd)}")

    result = subprocess.run(cmd, cwd=str(work_dir), env=env,
                           capture_output=True, text=True, timeout=7200)

    if result.returncode != 0:
        print(f"\nARTIS FAILED (exit code {result.returncode})")
        print("STDERR (last 2000 chars):")
        print(result.stderr[-2000:])
        # Also check log files
        for log_file in sorted(work_dir.glob("output_*.txt")):
            print(f"\n--- {log_file.name} (last 1000 chars) ---")
            content = log_file.read_text(errors='replace')
            print(content[-1000:])
    else:
        print(f"\nARTIS completed successfully")

    return result.returncode


# ---------------------------------------------------------------------------
# Step 10: Analyze spectrum
# ---------------------------------------------------------------------------
def read_spec_out(spec_path):
    """Read ARTIS spec.out file.

    Format:
      Line 1:   "0 t1 t2 ... tN"  (leading 0, then timestep mid-times in days)
      Lines 2+: "freq flux_ts1 flux_ts2 ... flux_tsN"  (one line per freq bin)

    Returns:
      freq:       np.array of frequency bin centers [Hz]
      wl:         np.array of wavelengths [Angstrom]
      ts_times:   np.array of timestep mid-times [days]
      flux_ts:    np.array shape (n_freq, n_timesteps) of flux values
    """
    with open(spec_path) as f:
        lines = f.readlines()

    # First line: leading 0 + timestep mid-times
    header = np.array([float(x) for x in lines[0].split()])
    ts_times = header[1:]  # skip leading 0

    # Subsequent lines: freq + per-timestep flux
    freq_list = []
    flux_list = []
    for line in lines[1:]:
        vals = [float(x) for x in line.split()]
        freq_list.append(vals[0])
        flux_list.append(vals[1:])

    freq = np.array(freq_list)
    flux_ts = np.array(flux_list)  # shape (n_freq, n_timesteps)
    wl = C_LIGHT / freq * 1e8  # Angstrom

    return freq, wl, ts_times, flux_ts


def analyze_spectrum(work_dir):
    """Read and analyze the ARTIS output spectrum."""
    work_dir = Path(work_dir)

    spec_path = work_dir / "spec.out"
    if not spec_path.exists():
        print("ERROR: No spec.out found — ARTIS may have failed")
        return False

    freq, wl, ts_times, flux_ts = read_spec_out(spec_path)
    n_timesteps = len(ts_times)

    print(f"\n{'='*60}")
    print(f"ARTIS Spectrum: {len(freq)} freq bins, {n_timesteps} timesteps")
    print(f"  Wavelength range: {wl.min():.0f} - {wl.max():.0f} A")
    print(f"  Time range: {ts_times[0]:.1f} - {ts_times[-1]:.1f} d")
    print(f"{'='*60}")

    # Time-integrated flux (sum over timesteps)
    flux_sum = flux_ts.sum(axis=1)
    # Last timestep flux
    flux_last = flux_ts[:, -1]

    print(f"  Flux range (sum): {flux_sum.min():.3e} - {flux_sum.max():.3e}")
    print(f"  Total flux (sum): {flux_sum.sum():.3e}")

    # Spectrum shape (1000 A bins)
    print(f"\n  Spectrum shape (1000 A bins):")
    for w0 in range(3000, 12000, 1000):
        mask = (wl >= w0) & (wl < w0 + 1000)
        if mask.sum() > 0:
            f_mean = flux_sum[mask].mean()
            f_std = flux_sum[mask].std()
            cv = f_std / f_mean if f_mean > 0 else 0
            tag = '(lines detected)' if cv > 0.15 else ''
            print(f"    {w0:5d}-{w0+1000:5d} A: mean={f_mean:.3e} "
                  f"CV={cv:.4f} {tag}")

    # Key diagnostic lines
    print(f"\n  Key diagnostic lines:")
    diagnostics = [
        ('Ca II HK', 3933), ('H delta', 4102), ('H gamma', 4340),
        ('H beta', 4861), ('Fe II', 5169), ('Si II 6355', 6355),
        ('H alpha', 6563), ('Ca II IR', 8662), ('O I 7774', 7774),
        ('He I 5876', 5876), ('C II 6580', 6580),
    ]
    for name, w in diagnostics:
        idx = np.argmin(np.abs(wl - w))
        mask = ((wl >= w - 100) & (wl < w - 30)) | \
               ((wl > w + 30) & (wl <= w + 100))
        if mask.sum() > 0 and flux_sum[idx] > 0:
            cont = flux_sum[mask].mean()
            depth = (cont - flux_sum[idx]) / cont if cont > 0 else 0
            tag = 'ABSORPTION' if depth > 0.02 else \
                  'EMISSION' if depth < -0.02 else 'flat'
            print(f"    {name:15s} ({w:5d} A): flux={flux_sum[idx]:.3e} "
                  f"cont={cont:.3e} depth={depth:+.3f} [{tag}]")
        else:
            print(f"    {name:15s} ({w:5d} A): flux={flux_sum[idx]:.3e} (no continuum)")

    return True


# ---------------------------------------------------------------------------
# Step 10b: Plot input structure
# ---------------------------------------------------------------------------
def plot_input_structure(work_dir, r_io, v_io, rho_io, comp_io, t_gas, q,
                          t_model_d, elems):
    """Generate diagnostic plots of the ARTIS input model structure.

    Multi-panel figure:
      1. Density vs velocity
      2. Temperature vs velocity
      3. Initial energy density (q) vs velocity
      4. Elemental abundance profiles vs velocity
      5. Mass profile (cumulative mass vs velocity)
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec

    work_dir = Path(work_dir)

    # Compute cell mid-radius and velocity
    r_edges = np.concatenate([[0], r_io])
    r_mid = 0.5 * (r_edges[:-1] + r_edges[1:])
    v_mid = v_io / 1e5  # km/s
    # Use outer velocity of each cell for plotting (ARTIS model.txt convention)
    v_outer = v_io / 1e5

    # Cumulative mass
    vol = 4.0 / 3.0 * PI * (r_edges[1:]**3 - r_edges[:-1]**3)
    mass = vol * rho_io
    mass_cumul = np.cumsum(mass) / M_SUN

    # Colors for elements
    elem_colors = {
        'H': 'royalblue', 'He': 'skyblue', 'C': 'darkgreen',
        'O': 'forestgreen', 'Si': 'crimson', 'Ca': 'darkorange',
        'Fe': 'firebrick', 'Ni': 'gold', 'S': 'mediumpurple',
        'Ar': 'plum', 'Ne': 'teal', 'Mg': 'olive',
    }

    fig = plt.figure(figsize=(16, 18))
    gs = GridSpec(4, 2, figure=fig, hspace=0.35, wspace=0.25,
                  height_ratios=[1, 1, 1, 1])

    # Panel 1: Density vs velocity
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.step(v_outer, rho_io, where='mid', color='steelblue', linewidth=1.2)
    ax1.set_yscale('log')
    ax1.set_xlabel('Velocity [km/s]')
    ax1.set_ylabel('Density [g/cm^3]')
    ax1.set_title(f'Density profile at t={t_model_d:.0f} d', fontsize=10)
    ax1.grid(True, alpha=0.2)
    ax1.text(0.02, 0.95, f'Max: {rho_io.max():.2e}\nMin: {rho_io.min():.2e}',
             transform=ax1.transAxes, fontsize=8, va='top',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    # Panel 2: Temperature vs velocity
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.step(v_outer, t_gas, where='mid', color='crimson', linewidth=1.2)
    ax2.axhline(T_BB, color='red', linestyle='--', alpha=0.4, label=f'T_BB={T_BB:.0f}K')
    ax2.axhline(T_OUTER, color='blue', linestyle='--', alpha=0.4, label=f'T_outer={T_OUTER:.0f}K')
    ax2.set_xlabel('Velocity [km/s]')
    ax2.set_ylabel('Gas temperature [K]')
    ax2.set_title('Temperature profile', fontsize=10)
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.2)

    # Panel 3: Initial energy density (q) vs velocity
    ax3 = fig.add_subplot(gs[1, 0])
    ax3.step(v_outer, q, where='mid', color='darkorange', linewidth=1.2)
    ax3.set_yscale('log')
    ax3.set_xlabel('Velocity [km/s]')
    ax3.set_ylabel('q [erg/g]')
    ax3.set_title('Initial thermal energy density (q column)', fontsize=10)
    ax3.grid(True, alpha=0.2)
    q_max = q[q > 0].max() if (q > 0).any() else 1
    ax3.text(0.02, 0.95, f'q_hot: {q_max:.3e}\nTotal E: {(q * mass).sum():.3e} erg',
             transform=ax3.transAxes, fontsize=8, va='top',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    if (q == 0).any():
        ax3.set_ylim(q_max * 1e-10, q_max * 10)

    # Panel 4: Elemental abundances vs velocity
    ax4 = fig.add_subplot(gs[1, 1])
    for elem in elems:
        elem_idx = ALL_ELEMENTS.index(elem)
        frac = comp_io[:, elem_idx]
        if frac.max() > 1e-6:
            color = elem_colors.get(elem, 'gray')
            ax4.step(v_outer, frac, where='mid', label=elem,
                     color=color, linewidth=1)
    ax4.set_yscale('log')
    ax4.set_xlabel('Velocity [km/s]')
    ax4.set_ylabel('Mass fraction')
    ax4.set_title('Elemental abundance profiles', fontsize=10)
    ax4.legend(fontsize=7, ncol=2, loc='upper right')
    ax4.set_ylim(1e-4, 1)
    ax4.grid(True, alpha=0.2)

    # Panel 5: Cumulative mass vs velocity
    ax5 = fig.add_subplot(gs[2, 0])
    ax5.step(v_outer, mass_cumul, where='mid', color='darkgreen', linewidth=1.2)
    ax5.set_xlabel('Velocity [km/s]')
    ax5.set_ylabel('Cumulative mass [Msun]')
    ax5.set_title(f'Cumulative mass (total: {mass_cumul[-1]:.4f} Msun)', fontsize=10)
    ax5.grid(True, alpha=0.2)
    # Mark H/He boundary
    h_idx = ALL_ELEMENTS.index('H')
    h_boundary = np.where(comp_io[:, h_idx] > 0.01)[0]
    if len(h_boundary) > 0:
        v_h_start = v_outer[h_boundary[0]]
        ax5.axvline(v_h_start, color='blue', linestyle=':', alpha=0.5,
                    label=f'H/He boundary ({v_h_start:.0f} km/s)')
        ax5.legend(fontsize=8)

    # Panel 6: Radius vs velocity (check homologous expansion)
    ax6 = fig.add_subplot(gs[2, 1])
    ax6.step(v_outer, r_io / 1e15, where='mid', color='purple', linewidth=1.2)
    ax6.set_xlabel('Velocity [km/s]')
    ax6.set_ylabel('Radius [1e15 cm]')
    ax6.set_title(f'Radius profile (R_max: {r_io.max()/1e15:.2f}e15 cm)', fontsize=10)
    ax6.grid(True, alpha=0.2)

    # Panel 7: Composition stacked area vs velocity
    ax7 = fig.add_subplot(gs[3, :])
    # Build stacked arrays
    stack_data = []
    stack_labels = []
    stack_colors = []
    for elem in elems:
        elem_idx = ALL_ELEMENTS.index(elem)
        frac = comp_io[:, elem_idx]
        if frac.max() > 1e-6:
            stack_data.append(frac)
            stack_labels.append(elem)
            stack_colors.append(elem_colors.get(elem, 'gray'))
    if stack_data:
        stack_arr = np.array(stack_data)  # shape (n_elem, n_cells)
        ax7.stackplot(v_outer, stack_arr, labels=stack_labels,
                      colors=stack_colors, alpha=0.7)
    ax7.set_xlabel('Velocity [km/s]')
    ax7.set_ylabel('Mass fraction')
    ax7.set_title('Composition stratification (stacked)', fontsize=10)
    ax7.legend(fontsize=7, ncol=7, loc='upper right')
    ax7.set_xlim(v_outer.min(), v_outer.max())
    ax7.set_ylim(0, 1)
    ax7.grid(True, alpha=0.2)

    fig.suptitle(f'ARTIS input structure: SN_A_g_T20k_Ni06 expanded to {t_model_d:.0f} d\n'
                 f'({len(r_io)} cells, {len(elems)} elements, '
                 f'total mass {mass_cumul[-1]:.4f} Msun)',
                 fontsize=12, y=0.98)

    plot_path = work_dir / "input_structure.png"
    fig.savefig(plot_path, dpi=150, bbox_inches='tight')
    print(f"\nSaved input structure plot to {plot_path}")
    plt.close()



def plot_spectrum(work_dir):
    """Generate diagnostic spectrum plots (3000-10000 A).

    Produces a multi-panel figure:
      1. Full spectrum with step plot and Planck comparison
      2. Smoothed spectrum with line-depth diagnostics
      3. Residual (spectrum / Planck) showing absorption features
      4. Per-timestep light curve evolution
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec

    work_dir = Path(work_dir)
    spec_path = work_dir / "spec.out"

    if not spec_path.exists():
        print("ERROR: No spec.out found for plotting")
        return

    freq, wl_arr, ts_times, flux_ts = read_spec_out(spec_path)

    # Time-integrated flux (sum over all timesteps)
    flux_total = flux_ts.sum(axis=1)
    # Last timestep flux
    flux_last_arr = flux_ts[:, -1]

    # Clip to 3000-10000 A
    mask = (wl_arr >= 3000) & (wl_arr <= 10000)
    wl_plot = wl_arr[mask]
    flux_plot = flux_total[mask]
    flux_last_plot = flux_last_arr[mask]

    # Smoothed spectrum (running mean)
    window = 15
    if len(flux_plot) > window:
        smooth = np.convolve(flux_plot, np.ones(window) / window, mode='same')
        smooth[:window // 2] = smooth[window // 2]
        smooth[-window // 2:] = smooth[-window // 2 - 1]
    else:
        smooth = flux_plot

    # 7000K Planck for comparison
    T_bb = T_BB
    wl_cm = wl_plot * 1e-8
    planck = (2 * H_PLANCK * C_LIGHT**2 / wl_cm**5) / \
             (np.exp(H_PLANCK * C_LIGHT / (wl_cm * KB * T_bb)) - 1)
    cont_mask = (wl_plot > 4000) & (wl_plot < 5500)
    if cont_mask.any() and smooth[cont_mask].mean() > 0:
        planck_norm = planck / np.mean(planck[cont_mask]) * np.mean(smooth[cont_mask])
    else:
        planck_norm = planck * 0

    # Residual: spectrum / planck (shows absorption relative to continuum)
    residual = np.zeros_like(smooth)
    valid = planck_norm > 0
    residual[valid] = smooth[valid] / planck_norm[valid]

    # Diagnostic lines
    lines_to_mark = [
        (3933, 'Ca II H&K'), (4102, 'H-delta'), (4340, 'H-gamma'),
        (4861, 'H-beta'), (5169, 'Fe II 5169'), (5876, 'He I 5876'),
        (6355, 'Si II 6355'), (6563, 'H-alpha'), (7065, 'He I 7065'),
        (7774, 'O I 7774'), (8498, 'Ca II IR'), (8542, 'Ca II IR'),
        (8662, 'Ca II IR'),
    ]

    # Compute line depths for diagnostic panel
    line_depths = []
    for w, name in lines_to_mark:
        idx = np.argmin(np.abs(wl_plot - w))
        cont_lo = ((wl_plot >= w - 80) & (wl_plot < w - 25))
        cont_hi = ((wl_plot > w + 25) & (wl_plot <= w + 80))
        cont_mask_l = cont_lo | cont_hi
        if cont_mask_l.any() and smooth[idx] > 0:
            cont = np.median(smooth[cont_mask_l])
            depth = (cont - smooth[idx]) / cont if cont > 0 else 0
        else:
            cont = 0
            depth = 0
        line_depths.append((w, name, depth, smooth[idx], cont))

    # Read light curve for panel 4
    lc_path = work_dir / "light_curve.out"
    lc_t, lc_lum = [], []
    if lc_path.exists():
        with open(lc_path) as f:
            for line in f:
                vals = line.split()
                if len(vals) >= 2:
                    lc_t.append(float(vals[0]))
                    lc_lum.append(float(vals[1]))

    # ---- Figure ----
    fig = plt.figure(figsize=(16, 14))
    gs = GridSpec(3, 2, figure=fig, hspace=0.35, wspace=0.25,
                  height_ratios=[1.2, 1, 1])

    # Panel 1: Full spectrum (step plot) + Planck
    ax1 = fig.add_subplot(gs[0, :])
    ax1.step(wl_plot, flux_plot, where='mid', color='steelblue',
             linewidth=0.5, alpha=0.5, label='ARTIS (time-integrated, raw)')
    ax1.step(wl_plot, smooth, where='mid', color='crimson',
             linewidth=1.2, alpha=0.9, label='Smoothed (15-bin)')
    ax1.step(wl_plot, planck_norm, where='mid', color='black',
             linewidth=1, alpha=0.5, linestyle='--', label=f'{T_bb:.0f}K Planck (normalized)')

    ymax = np.percentile(flux_plot[10:-10], 99) * 1.5 if len(flux_plot) > 20 else 1
    for w, name in lines_to_mark:
        ax1.axvline(w, color='gray', alpha=0.35, linestyle=':', linewidth=0.8)
        ax1.text(w, ymax * 0.97, name, rotation=90, fontsize=6,
                 ha='center', va='top', color='dimgray')
    ax1.set_xlim(3000, 10000)
    ax1.set_ylim(0, ymax)
    ax1.set_xlabel('Wavelength [Angstrom]')
    ax1.set_ylabel('Flux [arb. units]')
    ax1.set_title('ARTIS SN_B test: SN_A ejecta at 400d + 7000K inner BB source\n'
                  '(SN_A_g_T20k_Ni06, homologously expanded, 400 cells, 500k packets/rank)',
                  fontsize=11)
    ax1.legend(loc='upper right', fontsize=8)
    ax1.grid(True, alpha=0.2)

    # Panel 2: Residual (spectrum / planck)
    ax2 = fig.add_subplot(gs[1, 0])
    ax2.step(wl_plot, residual, where='mid', color='darkgreen', linewidth=0.8)
    ax2.axhline(1.0, color='black', linestyle='--', linewidth=0.5, alpha=0.5)
    for w, name in lines_to_mark:
        ax2.axvline(w, color='gray', alpha=0.3, linestyle=':', linewidth=0.6)
    ax2.set_xlim(3000, 10000)
    ax2.set_ylim(0, 2)
    ax2.set_xlabel('Wavelength [Angstrom]')
    ax2.set_ylabel('Spectrum / Planck')
    ax2.set_title('Residual: deviations from 7000K Planck', fontsize=10)
    ax2.grid(True, alpha=0.2)

    # Panel 3: Line depth bar chart
    ax3 = fig.add_subplot(gs[1, 1])
    names = [ld[1] for ld in line_depths]
    depths = [ld[2] for ld in line_depths]
    colors = ['firebrick' if d > 0.02 else 'steelblue' if d < -0.02
              else 'lightgray' for d in depths]
    y_pos = np.arange(len(names))
    ax3.barh(y_pos, depths, color=colors, height=0.6)
    ax3.set_yticks(y_pos)
    ax3.set_yticklabels(names, fontsize=7)
    ax3.axvline(0, color='black', linewidth=0.5)
    ax3.set_xlabel('Line depth (positive = absorption)')
    ax3.set_title('Diagnostic line depths', fontsize=10)
    ax3.grid(True, alpha=0.2, axis='x')

    # Panel 4: Light curve
    ax4 = fig.add_subplot(gs[2, :])
    if lc_t:
        ax4.plot(lc_t, lc_lum, 'b-', linewidth=1.5)
        ax4.axvline(400, color='gray', linestyle=':', alpha=0.5, label='t_model = 400d')
        ax4.set_xlabel('Time [days]')
        ax4.set_ylabel('Bolometric luminosity [arb. units]')
        ax4.set_title('Bolometric light curve (packets escaping per timestep)', fontsize=10)
        ax4.legend(fontsize=8)
        ax4.grid(True, alpha=0.2)
    else:
        ax4.text(0.5, 0.5, 'No light_curve.out found', ha='center', va='center',
                 transform=ax4.transAxes)

    # Save
    plot_path = work_dir / "spectrum_artis.png"
    fig.savefig(plot_path, dpi=150, bbox_inches='tight')
    print(f"\nSaved spectrum plot to {plot_path}")
    plt.close()

    # Print line depth table
    print(f"\n  Line depth diagnostics (3000-10000 A):")
    print(f"  {'Line':<18s} {'WL(A)':>7s} {'Flux':>12s} {'Continuum':>12s} {'Depth':>8s} {'Type':>10s}")
    print(f"  {'-'*18} {'-'*7} {'-'*12} {'-'*12} {'-'*8} {'-'*10}")
    for w, name, depth, flux_val, cont in line_depths:
        tag = 'ABSORP' if depth > 0.02 else 'EMIT' if depth < -0.02 else 'flat'
        print(f"  {name:<18s} {w:7d} {flux_val:12.3e} {cont:12.3e} {depth:+8.3f} {tag:>10s}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="ARTIS standalone test: SN_A ejecta + 7000K inner blackbody")
    parser.add_argument('--nrank', type=int, default=8,
                        help='Number of MPI ranks (default: 8)')
    parser.add_argument('--mpkts', type=int, default=500000,
                        help='Packets per rank (default: 500000)')
    parser.add_argument('--target-age', type=float, default=TARGET_AGE_D,
                        help=f'Target age in days (default: {TARGET_AGE_D})')
    parser.add_argument('--t-max', type=float, default=50.0,
                        help='Simulation duration in days (default: 50)')
    parser.add_argument('--n-timesteps', type=int, default=50,
                        help='Number of timesteps (default: 50)')
    parser.add_argument('--n-lte', type=int, default=5,
                        help='Number of LTE timesteps (default: 5)')
    parser.add_argument('--max-cells', type=int, default=MAX_CELLS,
                        help=f'Max grid cells (default: {MAX_CELLS})')
    parser.add_argument('--work-dir', type=str,
                        default='/viper/ptmp1/sbrennan/SN2024pba/tests/artis_sn_b_test',
                        help='Working directory')
    args = parser.parse_args()

    work_dir = Path(args.work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("ARTIS SN_B Test: SN_A ejecta + 7000K inner blackbody")
    print("=" * 60)
    print(f"  SN_A model: SN_A_g_T20k_Ni06")
    print(f"  Target age: {args.target_age:.0f} d")
    print(f"  T_BB: {T_BB} K, T_outer: {T_OUTER} K")
    print(f"  MPI ranks: {args.nrank}")
    print(f"  Packets/rank: {args.mpkts}")
    print(f"  Work dir: {work_dir}")

    # Step 1: Read SN_A profile
    epoch, t_snap = read_sn_a_profile(SN_A_SWD)

    # Step 2: Homologous expansion
    r, rho, v, total_mass = homologous_expand(epoch, t_snap, args.target_age)

    # Step 3: Detect elements
    elems = detect_elements(epoch['composition'])

    # Step 4: Downsample grid
    r_io, rho_io, v_io, comp_io, T_gas_io, T_rad_io, n_orig = downsample_grid(
        r, rho, v, epoch['composition'], epoch['T'], epoch['T_rad'],
        args.max_cells)

    # Step 5: Set temperature and energy
    t_gas, q, l_bol = set_temperature_and_energy(
        r_io, comp_io, rho_io, T_BB, T_OUTER, args.target_age)

    # Step 6: Write ARTIS input files
    print(f"\nWriting ARTIS input files to {work_dir}")
    write_model_txt(work_dir / "model.txt", r_io, v_io, rho_io,
                    t_gas, q, args.target_age)
    write_abundances_txt(work_dir / "abundances.txt", comp_io, elems)
    write_input_txt(work_dir / "input.txt", args.target_age, args.t_max,
                    args.n_timesteps, args.n_lte, args.mpkts)
    write_vpkt_txt(work_dir / "vpkt.txt")

    # Step 6b: Plot input structure
    print(f"\nPlotting input structure...")
    plot_input_structure(work_dir, r_io, v_io, rho_io, comp_io, t_gas, q,
                         args.target_age, elems)

    # Step 7: Run ARTIS
    print(f"\nSetting up and running ARTIS...")
    retcode = setup_and_run_artis(work_dir, args.nrank, args.mpkts,
                                  args.target_age, args.t_max)

    if retcode != 0:
        print(f"\nARTIS failed with exit code {retcode}")
        return retcode

    # Step 8: Analyze spectrum
    print(f"\nAnalyzing spectrum...")
    success = analyze_spectrum(work_dir)

    # Step 9: Plot
    if success:
        plot_spectrum(work_dir)

    return 0 if success else 1


if __name__ == '__main__':
    sys.exit(main())
