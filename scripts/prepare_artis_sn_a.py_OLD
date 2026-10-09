#!/usr/bin/env python3
"""Prepare ARTIS input files from a SN_A STELLA .swd snapshot at a given epoch.

This script extracts the STELLA epoch closest to the requested snapshot day
(e.g. 30d post-SBO) and writes ARTIS model.txt, abundances.txt, input.txt,
vpkt.txt, and a SLURM script.  No homologous evolution or blackbody
illumination is applied — the pure SN_A ejecta structure is used directly.

Usage:
    python3 scripts/prepare_artis_sn_a.py --swd SWD_PATH --snapshot-day 30
        --work-dir DIR --nrank 2 --mpkts 500000
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
SIGMA_SB = 5.670374419e-5
C_LIGHT = 2.99792458e10
PI = np.pi
M_SUN = 1.989e33
DAY_S = 86400.0
KB = 1.380649e-16
H_PLANCK = 6.62607015e-27

DEFAULT_SN_A_SWD = "/viper/ptmp1/sbrennan/SN2024pba/SN_A/SN_A_g_T20k_Ni06/stella/SN_A_g_T20k_Ni06.swd"
ARTIS_REPO = "/viper/ptmp1/sbrennan/artis_repo"
ARTIS_EXE = f"{ARTIS_REPO}/sn3d"
ARTIS_ATOMIC = "/viper/ptmp1/sbrennan/artisatomic/artis_files"

STELLA_ELEMENTS = ['H', 'He', 'C', 'N', 'O', 'Ne', 'Na', 'Mg', 'Al',
                   'Si', 'S', 'Ar', 'Ca', 'Fe', 'Ni']

ELEMENT_Z = {'H': 1, 'He': 2, 'C': 6, 'N': 7, 'O': 8, 'Ne': 10,
             'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'S': 16,
             'Ar': 18, 'Ca': 20, 'Fe': 26, 'Ni': 28}

ARTIS_ELEMENTS = {'H', 'He', 'C', 'O', 'Si', 'Ca', 'Fe'}

ELEMENT_REMAP = {
    'N': 'C', 'Ne': 'O', 'Na': 'Si', 'Mg': 'Si', 'Al': 'Si',
    'S': 'Si', 'Ar': 'Ca', 'Ni': 'Fe',
}

DEFAULT_T_MAX_D = 50.0
DEFAULT_N_TIMESTEPS = 50
DEFAULT_N_LTE = 5
DEFAULT_MAX_CELLS = 400
DEFAULT_N_RANK = 2
DEFAULT_MPKTS = 500000


def read_swd_epoch(swd_path, snapshot_day):
    """Read the STELLA .swd and extract the epoch closest to snapshot_day."""
    sys.path.insert(0, str(Path(__file__).parent))
    from convert_stella_swd_to_cmfgen import read_swd, normalize_swd_composition

    data = read_swd(swd_path)
    times = sorted([float(k) for k in data.keys()])

    idx = int(np.argmin([abs(t - snapshot_day) for t in times]))
    t_snap = times[idx]
    epoch = data[t_snap]

    comp_dict = normalize_swd_composition(epoch['composition'])
    comp = np.zeros((len(epoch['r']), len(STELLA_ELEMENTS)))
    for i, elem in enumerate(STELLA_ELEMENTS):
        comp[:, i] = comp_dict[elem]

    v_snap = epoch['r'] / (t_snap * DAY_S)

    print(f"Read SN_A .swd: {swd_path}")
    print(f"  Snapshot epoch: t={t_snap:.4f} d (requested {snapshot_day:.1f}d)")
    print(f"  Zones: {len(epoch['r'])}")
    print(f"  R: {epoch['r'].min():.3e} - {epoch['r'].max():.3e} cm")
    print(f"  V: {v_snap.min()/1e5:.1f} - {v_snap.max()/1e5:.1f} km/s")
    print(f"  T_gas: {epoch['T'].min():.1f} - {epoch['T'].max():.1f} K")
    print(f"  T_rad: {epoch['T_rad'].min():.1f} - {epoch['T_rad'].max():.1f} K")
    print(f"  rho: {epoch['rho'].min():.3e} - {epoch['rho'].max():.3e} g/cm^3")
    print(f"  L_max: {epoch['flux'].max():.3e} erg/s")

    result = {
        'r': epoch['r'],
        'v_cms': v_snap,
        'rho': epoch['rho'],
        'T': epoch['T'],
        'T_rad': epoch['T_rad'],
        'flux': epoch['flux'],
        'composition_normalized': comp,
    }
    return result, t_snap


def map_composition_to_artis(comp_stella):
    """Map 15-element STELLA composition to 7 ARTIS-available elements."""
    n = comp_stella.shape[0]
    comp_artis = np.zeros((n, len(STELLA_ELEMENTS)))

    for i, elem in enumerate(STELLA_ELEMENTS):
        if elem in ARTIS_ELEMENTS:
            comp_artis[:, i] = comp_stella[:, i]
        elif elem in ELEMENT_REMAP:
            target = ELEMENT_REMAP[elem]
            target_idx = STELLA_ELEMENTS.index(target)
            comp_artis[:, target_idx] += comp_stella[:, i]

    totals = comp_artis.sum(axis=1, keepdims=True)
    totals = np.where(totals < 1e-30, 1.0, totals)
    comp_artis = comp_artis / totals

    print("\nComposition mapping (STELLA -> ARTIS):")
    for elem in STELLA_ELEMENTS:
        if elem in ARTIS_ELEMENTS:
            print(f"  {elem}: kept (max={comp_stella[:, STELLA_ELEMENTS.index(elem)].max():.4f})")
        elif elem in ELEMENT_REMAP:
            print(f"  {elem} -> {ELEMENT_REMAP[elem]}")

    return comp_artis


def adaptive_downsample(r, rho, v, comp, T_gas, T_rad, flux, max_cells):
    """Downsample the grid using adaptive sampling. SWD is outer-to-inner;
    ARTIS model.txt is inner-to-outer."""
    n_orig = len(r)

    r_io = r[::-1]
    rho_io_all = rho[::-1]
    v_io_all = v[::-1]
    comp_io_all = comp[::-1]
    T_io_all = T_gas[::-1]
    Trad_io_all = T_rad[::-1]
    flux_io_all = flux[::-1]

    if n_orig <= max_cells:
        print(f"\nNo downsampling needed: {n_orig} zones")
        return (r_io, rho_io_all, v_io_all, comp_io_all,
                T_io_all, Trad_io_all, flux_io_all)

    log_T = np.log10(np.maximum(T_io_all, 1.0))
    dlog_T = np.abs(np.diff(log_T, prepend=log_T[0]))

    log_L = np.log10(np.maximum(np.abs(flux_io_all), 1.0))
    dlog_L = np.abs(np.diff(log_L, prepend=log_L[0]))

    log_rho = np.log10(np.maximum(rho_io_all, 1e-30))
    dlog_rho = np.abs(np.diff(log_rho, prepend=log_rho[0]))

    comp_grad = np.zeros(n_orig)
    for i in range(1, n_orig):
        comp_grad[i] = np.sum(np.abs(comp_io_all[i] - comp_io_all[i-1]))

    log_v = np.log10(np.maximum(v_io_all, 1.0))
    dlog_v = np.abs(np.diff(log_v, prepend=log_v[0]))

    def norm(arr):
        mx = arr.max()
        return arr / mx if mx > 0 else arr

    weights = (0.30 * norm(dlog_T) +
               0.25 * norm(dlog_L) +
               0.15 * norm(dlog_rho) +
               0.15 * norm(comp_grad) +
               0.15 * norm(dlog_v))

    weights[0] = 1.0
    weights[-1] = 1.0
    weights[T_io_all > 5000] = np.maximum(weights[T_io_all > 5000], 0.8)
    pos_flux = flux_io_all > 1e38
    weights[pos_flux] = np.maximum(weights[pos_flux], 0.5)

    prob = weights / weights.sum()
    cumprob = np.cumsum(prob)
    uniform = np.linspace(0, cumprob[-1], max_cells, endpoint=False)
    indices = np.searchsorted(cumprob, uniform)
    indices = np.clip(indices, 0, n_orig - 1)
    indices = np.unique(indices)

    if len(indices) > max_cells:
        indices = np.linspace(0, n_orig - 1, max_cells, dtype=int)

    r_ds = r_io[indices]
    rho_ds = rho_io_all[indices]
    v_ds = v_io_all[indices]
    comp_ds = comp_io_all[indices]
    T_ds = T_io_all[indices]
    Trad_ds = Trad_io_all[indices]
    flux_ds = flux_io_all[indices]

    print(f"\nAdaptive downsampling: {n_orig} -> {len(indices)} zones")
    print(f"  R range: {r_ds.min():.3e} - {r_ds.max():.3e} cm")
    print(f"  V range: {v_ds.min()/1e5:.1f} - {v_ds.max()/1e5:.1f} km/s")
    print(f"  T range: {T_ds.min():.1f} - {T_ds.max():.1f} K")

    return r_ds, rho_ds, v_ds, comp_ds, T_ds, Trad_ds, flux_ds


def write_model_txt(path, r_io, v_io, rho_io, t_gas, q, t_model_d,
                    comp_io=None, ni56_mass_msun=0.06, total_mass_msun=1.72):
    """Write ARTIS model.txt file (1D spherical, inner-to-outer)."""
    n = len(r_io)
    tau_ni = 8.77
    tau_co = 111.5

    frac_ni = np.exp(-t_model_d / tau_ni)
    frac_co_from_ni = (1.0 - frac_ni)
    frac_co_remaining = np.exp(-t_model_d / tau_co)
    ni56_msun = ni56_mass_msun * frac_ni
    co56_msun = ni56_mass_msun * frac_co_from_ni * frac_co_remaining
    fe56_msun = ni56_mass_msun - ni56_msun - co56_msun

    x_ni56_global = ni56_msun / total_mass_msun
    x_co56_global = co56_msun / total_mass_msun
    x_fegrp_global = fe56_msun / total_mass_msun

    print(f"\nRadioactive decay at t={t_model_d:.0f}d (initial Ni56={ni56_mass_msun} Msun):")
    print(f"  Ni56: {ni56_msun:.6e} Msun (X={x_ni56_global:.6e})")
    print(f"  Co56: {co56_msun:.6e} Msun (X={x_co56_global:.6e})")
    print(f"  Fe56 (decay product): {fe56_msun:.6e} Msun (X={x_fegrp_global:.6e})")

    r_edges = np.concatenate([[0], r_io])
    vol = 4.0/3.0 * PI * (r_edges[1:]**3 - r_edges[:-1]**3)
    mass = vol * rho_io

    fe_idx = STELLA_ELEMENTS.index('Fe')
    ni_idx = STELLA_ELEMENTS.index('Ni')

    with open(path, 'w') as f:
        f.write(f"{n}\n")
        f.write(f"{t_model_d:.6f}\n")
        f.write("#inputcellid velocity_outer logrho X_Fegroup X_Ni56 X_Co56 X_Fe52 X_Cr48 q\n")

        for i in range(n):
            cellid = i + 1
            v_outer = v_io[i] / 1e5
            log_rho = np.log10(rho_io[i]) if rho_io[i] > 0 else -99.0

            fe_frac = comp_io[i, fe_idx] if comp_io is not None else 0.0
            if fe_frac > 0.01 and mass[i] > 0:
                scale = fe_frac / 0.5
                x_ni56 = max(x_ni56_global * scale, 1e-10)
                x_co56 = x_co56_global * scale
                x_fegrp = x_fegrp_global * scale
                x_ni56 = min(x_ni56, fe_frac * 0.5)
                x_co56 = min(x_co56, fe_frac * 0.5)
                x_fegrp = min(x_fegrp, fe_frac)
            else:
                x_ni56 = 1e-10
                x_co56 = 0.0
                x_fegrp = 0.0

            x_fe52 = 0.0
            x_cr48 = 0.0
            q_val = q[i]
            f.write(f"   {cellid} {v_outer:.6e} {log_rho:.6f} "
                    f"{x_fegrp:.6e} {x_ni56:.6e} {x_co56:.6e} "
                    f"{x_fe52:.6e} {x_cr48:.6e} {q_val:.6e}\n")

    print(f"\nWrote {path}: {n} cells, t_model={t_model_d:.1f} d")


def write_abundances_txt(path, comp_io, elems):
    """Write ARTIS abundances.txt file."""
    n = len(comp_io)
    max_z = 30

    with open(path, 'w') as f:
        for i in range(n):
            cellid = i + 1
            abund = np.zeros(max_z)
            for elem in elems:
                z = ELEMENT_Z.get(elem, 0)
                if z > 0 and z <= max_z and elem in ARTIS_ELEMENTS:
                    elem_idx = STELLA_ELEMENTS.index(elem)
                    abund[z - 1] = comp_io[i, elem_idx]

            total = abund.sum()
            if total > 0:
                abund = abund / total

            vals = ' '.join(f'{a:.8e}' for a in abund)
            f.write(f"{cellid} {vals}\n")

    print(f"Wrote {path}: {n} cells, Z=1..{max_z}")


def write_input_txt(path, t_model_d, t_max_d, n_timesteps, n_lte_ts, mpkts):
    """Write ARTIS input.txt file."""
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
        f"{mpkts}                  # nprocs_exspec",
        "1                        # do_emission_res (UNUSED)",
        "0.001 1000               # kpktdiffusion (UNUSED)",
    ]

    with open(path, 'w') as f:
        f.write('\n'.join(lines) + '\n')

    print(f"Wrote {path}: tmin={tmin:.1f}d, tmax={tmax:.1f}d, "
          f"{n_timesteps} timesteps, {n_lte_ts} LTE")


def write_vpkt_txt(path):
    """Write ARTIS vpkt.txt (virtual packets disabled)."""
    lines = ["0", "0 10 30", "0 100", "10", "0", "1", "11.5 21.5"]
    with open(path, 'w') as f:
        f.write('\n'.join(lines) + '\n')
    print(f"Wrote {path} (virtual packets disabled)")


def write_slurm_script(path, work_dir, n_rank, mpkts, t_model_d, t_max_d,
                       artis_exe=ARTIS_EXE, job_name=None):
    """Write a SLURM batch script for the ARTIS run."""
    if job_name is None:
        job_name = f"artis_sn_a_{t_model_d:.0f}d"

    script = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --partition=small
#SBATCH --qos=n0001
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks={n_rank}
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output={work_dir}/slurm_%j.out
#SBATCH --error={work_dir}/slurm_%j.err

set -e

# Load GCC 15 and OpenMPI 4.1 for ARTIS
module purge
module load gcc/15
module load openmpi/4.1
export OMPI_CXX=g++
export LD_LIBRARY_PATH=/mpcdf/soft/RHEL_9/packages/x86_64/gcc/15.1.0/lib64:/mpcdf/soft/RHEL_9/packages/znver4/gsl/gcc_15-15.1.0/2.7.1/lib:${{LD_LIBRARY_PATH}}
export LD_PRELOAD=/mpcdf/soft/RHEL_9/packages/x86_64/gcc/15.1.0/lib64/libstdc++.so.6

# Run ARTIS
cd {work_dir}
echo "Starting ARTIS at $(date)"
echo "Executable: {artis_exe}"
echo "MPI ranks: {n_rank}"
echo "Packets/rank: {mpkts}"
echo "t_model: {t_model_d:.1f} d, t_max: {t_max_d:.1f} d"

mpirun --bind-to none -np {n_rank} {artis_exe}

echo "ARTIS finished at $(date) with exit code $?"

# Check for spec.out
if [ -f spec.out ]; then
    echo "spec.out exists — spectrum produced"
else
    echo "ERROR: spec.out not found — ARTIS may have failed"
fi
"""
    with open(path, 'w') as f:
        f.write(script)
    print(f"Wrote {path}")


def copy_atomic_data(work_dir):
    """Copy atomic data files to the working directory."""
    work_dir = Path(work_dir)

    atomic_files = ['adata.txt', 'transitiondata.txt', 'compositiondata.txt',
                    'phixsdata_v2.txt']
    for af in atomic_files:
        src = Path(ARTIS_ATOMIC) / af
        dst = work_dir / af
        if src.exists():
            shutil.copy2(src, dst)
            print(f"  Copied {af} ({src.stat().st_size / 1e6:.1f} MB)")

    data_src = Path(ARTIS_REPO) / "data"
    data_dst = work_dir / "data"
    if data_src.exists() and not data_dst.exists():
        shutil.copytree(data_src, data_dst, dirs_exist_ok=True)
        print(f"  Copied data/ directory")

    for df in ['betaminusdecays.txt', 'alphadecays.txt', 'fissiondecays.txt']:
        src = Path(ARTIS_REPO) / "data" / df
        dst = work_dir / df
        if src.exists() and not dst.exists():
            shutil.copy2(src, dst)


def main():
    parser = argparse.ArgumentParser(
        description="Prepare ARTIS input files from SN_A STELLA .swd at a given epoch")
    parser.add_argument('--swd', type=str, default=DEFAULT_SN_A_SWD,
                        help=f'Path to SN_A .swd file')
    parser.add_argument('--snapshot-day', type=float, default=30.0,
                        help='Snapshot epoch in days (default: 30)')
    parser.add_argument('--t-max', type=float, default=DEFAULT_T_MAX_D,
                        help=f'Simulation duration in days (default: {DEFAULT_T_MAX_D})')
    parser.add_argument('--n-timesteps', type=int, default=DEFAULT_N_TIMESTEPS,
                        help=f'Number of timesteps (default: {DEFAULT_N_TIMESTEPS})')
    parser.add_argument('--n-lte', type=int, default=DEFAULT_N_LTE,
                        help=f'Number of LTE timesteps (default: {DEFAULT_N_LTE})')
    parser.add_argument('--max-cells', type=int, default=DEFAULT_MAX_CELLS,
                        help=f'Max grid cells (default: {DEFAULT_MAX_CELLS})')
    parser.add_argument('--work-dir', type=str, required=True,
                        help='Working directory')
    parser.add_argument('--nrank', type=int, default=DEFAULT_N_RANK,
                        help=f'MPI ranks (default: {DEFAULT_N_RANK})')
    parser.add_argument('--mpkts', type=int, default=DEFAULT_MPKTS,
                        help=f'Packets per rank (default: {DEFAULT_MPKTS})')
    parser.add_argument('--ni56-mass', type=float, default=0.06,
                        help='Initial Ni56 mass in Msun (default: 0.06)')
    parser.add_argument('--total-mass', type=float, default=1.72,
                        help='Total ejecta mass in Msun (default: 1.72)')
    parser.add_argument('--job-name', type=str, default=None,
                        help='SLURM job name')
    args = parser.parse_args()

    work_dir = Path(args.work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print(f"ARTIS Setup: SN_A ejecta at {args.snapshot_day:.0f}d post-SBO")
    print("=" * 70)
    print(f"  SN_A .swd: {args.swd}")
    print(f"  Snapshot day: {args.snapshot_day:.1f} d")
    print(f"  Sim duration: {args.t_max:.0f} d")
    print(f"  Max cells: {args.max_cells}")
    print(f"  MPI ranks: {args.nrank}")
    print(f"  Packets/rank: {args.mpkts}")
    print(f"  Work dir: {work_dir}")

    # Step 1: Read SN_A profile at snapshot epoch
    epoch, t_snap = read_swd_epoch(args.swd, args.snapshot_day)

    # Step 2: Map composition to ARTIS elements
    comp_artis = map_composition_to_artis(epoch['composition_normalized'])

    # Step 3: Adaptive downsampling
    r_io, rho_io, v_io, comp_io, T_io, Trad_io, flux_io = adaptive_downsample(
        epoch['r'], epoch['rho'], epoch['v_cms'], comp_artis,
        epoch['T'], epoch['T_rad'], epoch['flux'], args.max_cells)

    # Step 4: Use STELLA gas temperature directly (no blackbody illumination)
    t_gas = T_io
    q = np.zeros(len(r_io))  # No initial energy deposition — Ni56 decay handles heating

    print(f"\nTemperature profile (direct from STELLA):")
    print(f"  T_gas range: {t_gas.min():.0f} - {t_gas.max():.0f} K")
    print(f"  q: 0 (Ni56 decay heating via ARTIS)")

    # Step 5: Write ARTIS input files
    print(f"\nWriting ARTIS input files to {work_dir}")
    write_model_txt(work_dir / "model.txt", r_io, v_io, rho_io,
                    t_gas, q, t_snap, comp_io=comp_io,
                    ni56_mass_msun=args.ni56_mass, total_mass_msun=args.total_mass)
    write_abundances_txt(work_dir / "abundances.txt", comp_io, STELLA_ELEMENTS)
    write_input_txt(work_dir / "input.txt", t_snap, args.t_max,
                    args.n_timesteps, args.n_lte, args.mpkts)
    write_vpkt_txt(work_dir / "vpkt.txt")

    # Step 6: Copy atomic data
    print(f"\nCopying atomic data...")
    copy_atomic_data(work_dir)

    # Step 7: Write SLURM script
    write_slurm_script(work_dir / "run_artis.slurm", work_dir,
                       args.nrank, args.mpkts, t_snap, args.t_max,
                       job_name=args.job_name)

    # Summary
    print(f"\n{'=' * 70}")
    print(f"Setup complete:")
    print(f"  t_model: {t_snap:.4f} d")
    print(f"  Cells: {len(r_io)}")
    print(f"  R range: {r_io.min():.3e} - {r_io.max():.3e} cm")
    print(f"  V range: {v_io.min()/1e5:.1f} - {v_io.max()/1e5:.1f} km/s")
    print(f"  T range: {t_gas.min():.0f} - {t_gas.max():.0f} K")
    print(f"  Work dir: {work_dir}")
    print(f"  To run: cd {work_dir} && sbatch run_artis.slurm")
    print(f"{'=' * 70}")

    return 0


if __name__ == '__main__':
    sys.exit(main())
