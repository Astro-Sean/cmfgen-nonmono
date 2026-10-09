#!/usr/bin/env python3
"""Prepare and submit a series of ARTIS models for SN_A from 300 to 400 days.

For each epoch in the series, this script:
1. Reads the SN_A STELLA .swd at the latest epoch (~62 d).
2. Homologously expands the profile to the target age.
3. Sets the inner metal-rich zones to T_BB with thermal energy (q column).
4. Writes ARTIS input files: model.txt, abundances.txt, input.txt, vpkt.txt.
5. Copies atomic data files.
6. Creates a SLURM job script.
7. Submits the job.

Usage:
    python3 scripts/prepare_artis_series.py [--t-start 300] [--t-end 400] [--t-step 10]
           [--nrank 8] [--mpkts 500000] [--submit]
"""

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

# Reuse functions from the test script
sys.path.insert(0, str(Path(__file__).parent))
from test_artis_blackbody_illuminated_ejecta import (
    SN_A_SWD,
    ARTIS_REPO,
    ARTIS_EXE,
    ARTIS_ATOMIC,
    ALL_ELEMENTS,
    ELEMENT_Z,
    T_BB,
    T_OUTER,
    DAY_S,
    C_LIGHT,
    PI,
    M_SUN,
    read_sn_a_profile,
    homologous_expand,
    detect_elements,
    downsample_grid,
    set_temperature_and_energy,
    write_model_txt,
    write_abundances_txt,
    write_input_txt,
    write_vpkt_txt,
)

BASE_OUTPUT_DIR = "/viper/ptmp1/sbrennan/SN2024pba/tests/artis_sn_a_series"


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

    # Copy data/ directory
    data_src = Path(ARTIS_REPO) / "data"
    data_dst = work_dir / "data"
    if data_src.exists():
        shutil.copytree(data_src, data_dst, dirs_exist_ok=True)

    # Copy decay data files
    for df in ['betaminusdecays.txt', 'alphadecays.txt', 'fissiondecays.txt']:
        src = Path(ARTIS_REPO) / "data" / df
        dst = work_dir / df
        if src.exists():
            shutil.copy2(src, dst)


def write_slurm_script(work_dir, t_model_d, t_max_d, n_rank, mpkts, partition,
                       time_limit, mem):
    """Write a SLURM batch script for an ARTIS run."""
    work_dir = Path(work_dir)
    slurm_path = work_dir / "run_artis.slurm"

    tmin = t_model_d
    tmax = t_model_d + t_max_d

    lines = [
        "#!/bin/bash",
        f"#SBATCH --job-name=artis_{t_model_d:.0f}d",
        f"#SBATCH --partition={partition}",
        f"#SBATCH --time={time_limit}",
        f"#SBATCH --nodes=1",
        f"#SBATCH --ntasks-per-node={n_rank}",
        f"#SBATCH --mem={mem}",
        f"#SBATCH --output=artis_{t_model_d:.0f}d_%j.out",
        f"#SBATCH --error=artis_{t_model_d:.0f}d_%j.err",
        "",
        "set -e",
        "",
        f"echo 'ARTIS run: t_model={t_model_d:.0f}d, tmax={tmax:.0f}d, "
        f"ranks={n_rank}, mpkts={mpkts}'",
        f"echo 'Working directory: {work_dir}'",
        "echo 'Start: '$(date)",
        "",
        "# Load GCC 15 and OpenMPI 4.1 for ARTIS",
        "module purge",
        "module load gcc/15",
        "module load openmpi/4.1",
        "export OMPI_CXX=g++",
        "export LD_LIBRARY_PATH=/mpcdf/soft/RHEL_9/packages/x86_64/gcc/15.1.0/lib64:/mpcdf/soft/RHEL_9/packages/znver4/gsl/gcc_15-15.1.0/2.7.1/lib:$LD_LIBRARY_PATH",
        "",
        f"cd {work_dir}",
        "",
        f"mpirun --bind-to none -np {n_rank} {ARTIS_EXE}",
        "",
        "echo 'ARTIS finished with exit code: '$?",
        "echo 'End: '$(date)",
        "",
        "# Check for spec.out",
        "if [ -f spec.out ]; then",
        "    echo 'spec.out exists — spectrum produced'",
        "else",
        "    echo 'ERROR: spec.out not found — ARTIS may have failed'",
        "fi",
    ]

    with open(slurm_path, 'w') as f:
        f.write('\n'.join(lines) + '\n')

    os.chmod(slurm_path, 0o755)
    return slurm_path


def prepare_epoch(epoch_day, epoch_data, t_max_d, n_timesteps, n_lte,
                  n_rank, mpkts, max_cells, base_dir, partition, time_limit,
                  mem):
    """Prepare ARTIS input files for a single epoch."""
    work_dir = Path(base_dir) / f"{epoch_day:.0f}d"
    work_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n{'='*60}")
    print(f"Preparing epoch: {epoch_day:.0f} days")
    print(f"  Output dir: {work_dir}")
    print(f"{'='*60}")

    # Unpack shared epoch data
    epoch, t_snap = epoch_data['epoch'], epoch_data['t_snap']
    r, rho, v, total_mass = epoch_data['r'], epoch_data['rho'], epoch_data['v'], epoch_data['total_mass']
    elems = epoch_data['elems']

    # Downsample grid
    r_io, rho_io, v_io, comp_io, T_gas_io, T_rad_io, n_orig = downsample_grid(
        r, rho, v, epoch['composition'], epoch['T'], epoch['T_rad'],
        max_cells)

    # Set temperature and energy
    t_gas, q, l_bol = set_temperature_and_energy(
        r_io, comp_io, rho_io, T_BB, T_OUTER, epoch_day)

    # Write ARTIS input files
    write_model_txt(work_dir / "model.txt", r_io, v_io, rho_io,
                    t_gas, q, epoch_day)
    write_abundances_txt(work_dir / "abundances.txt", comp_io, elems)
    write_input_txt(work_dir / "input.txt", epoch_day, t_max_d,
                    n_timesteps, n_lte, mpkts)
    write_vpkt_txt(work_dir / "vpkt.txt")

    # Copy atomic data
    copy_atomic_data(work_dir)

    # Write SLURM script
    slurm_path = write_slurm_script(work_dir, epoch_day, t_max_d, n_rank,
                                    mpkts, partition, time_limit, mem)

    print(f"  SLURM script: {slurm_path}")
    return work_dir, slurm_path


def submit_job(slurm_path):
    """Submit a SLURM job and return the job ID."""
    result = subprocess.run(['sbatch', str(slurm_path)],
                           capture_output=True, text=True,
                           cwd=str(slurm_path.parent))
    if result.returncode != 0:
        print(f"  ERROR submitting {slurm_path}: {result.stderr}")
        return None

    # Parse job ID from "Submitted batch job 123456"
    job_id = result.stdout.strip().split()[-1]
    print(f"  Submitted job {job_id}")
    return job_id


def main():
    parser = argparse.ArgumentParser(
        description="Prepare and submit ARTIS model series for SN_A (300-400 d)")
    parser.add_argument('--t-start', type=float, default=300.0,
                        help='Start epoch in days (default: 300)')
    parser.add_argument('--t-end', type=float, default=400.0,
                        help='End epoch in days (default: 400)')
    parser.add_argument('--t-step', type=float, default=10.0,
                        help='Epoch step in days (default: 10)')
    parser.add_argument('--nrank', type=int, default=8,
                        help='MPI ranks per job (default: 8)')
    parser.add_argument('--mpkts', type=int, default=500000,
                        help='Packets per rank (default: 500000)')
    parser.add_argument('--t-max', type=float, default=50.0,
                        help='Simulation duration per epoch in days (default: 50)')
    parser.add_argument('--n-timesteps', type=int, default=50,
                        help='Number of timesteps (default: 50)')
    parser.add_argument('--n-lte', type=int, default=5,
                        help='Number of LTE timesteps (default: 5)')
    parser.add_argument('--max-cells', type=int, default=400,
                        help='Max grid cells (default: 400)')
    parser.add_argument('--partition', type=str, default='general',
                        help='SLURM partition (default: general)')
    parser.add_argument('--time', type=str, default='02:00:00',
                        help='Time limit per job (default: 02:00:00)')
    parser.add_argument('--mem', type=str, default='16G',
                        help='Memory per job (default: 16G)')
    parser.add_argument('--submit', action='store_true',
                        help='Submit SLURM jobs after preparing')
    parser.add_argument('--output-dir', type=str, default=BASE_OUTPUT_DIR,
                        help=f'Base output directory (default: {BASE_OUTPUT_DIR})')
    args = parser.parse_args()

    # Generate epoch list
    epochs = np.arange(args.t_start, args.t_end + args.t_step / 2, args.t_step)
    print(f"ARTIS SN_A series: {len(epochs)} epochs from "
          f"{epochs[0]:.0f} to {epochs[-1]:.0f} d (step {args.t_step:.0f} d)")
    print(f"  MPI ranks: {args.nrank}, packets/rank: {args.mpkts}")
    print(f"  Sim duration: {args.t_max:.0f} d, timesteps: {args.n_timesteps}")
    print(f"  Output: {args.output_dir}")
    print(f"  Submit: {args.submit}")

    # Read SN_A profile once (shared across all epochs)
    print(f"\nReading SN_A .swd: {SN_A_SWD}")
    epoch, t_snap = read_sn_a_profile(SN_A_SWD)

    # Detect elements once
    elems = detect_elements(epoch['composition'])

    # Prepare each epoch
    base_dir = Path(args.output_dir)
    base_dir.mkdir(parents=True, exist_ok=True)

    job_ids = []
    for epoch_day in epochs:
        # Homologous expansion for this epoch
        r, rho, v, total_mass = homologous_expand(epoch, t_snap, epoch_day)

        epoch_data = {
            'epoch': epoch,
            't_snap': t_snap,
            'r': r,
            'rho': rho,
            'v': v,
            'total_mass': total_mass,
            'elems': elems,
        }

        work_dir, slurm_path = prepare_epoch(
            epoch_day, epoch_data, args.t_max, args.n_timesteps, args.n_lte,
            args.nrank, args.mpkts, args.max_cells, base_dir,
            args.partition, args.time, args.mem)

        if args.submit:
            job_id = submit_job(slurm_path)
            if job_id:
                job_ids.append((epoch_day, job_id))

    # Summary
    print(f"\n{'='*60}")
    print(f"Summary: prepared {len(epochs)} ARTIS runs")
    if job_ids:
        print(f"Submitted {len(job_ids)} SLURM jobs:")
        for epoch_day, jid in job_ids:
            print(f"  {epoch_day:.0f}d -> job {jid}")
    else:
        print("No jobs submitted (use --submit to submit)")
        print("To submit manually:")
        for epoch_day in epochs:
            sd = base_dir / f"{epoch_day:.0f}d"
            print(f"  sbatch {sd}/run_artis.slurm")

    # Write summary file
    summary_path = base_dir / "series_summary.txt"
    with open(summary_path, 'w') as f:
        f.write(f"ARTIS SN_A series: {len(epochs)} epochs\n")
        f.write(f"Range: {epochs[0]:.0f} - {epochs[-1]:.0f} d (step {args.t_step:.0f} d)\n")
        f.write(f"SN_A model: {SN_A_SWD}\n")
        f.write(f"T_BB: {T_BB} K, T_outer: {T_OUTER} K\n")
        f.write(f"MPI ranks: {args.nrank}, packets/rank: {args.mpkts}\n")
        f.write(f"Sim duration: {args.t_max:.0f} d, timesteps: {args.n_timesteps}\n\n")
        if job_ids:
            f.write("Submitted jobs:\n")
            for epoch_day, jid in job_ids:
                f.write(f"  {epoch_day:.0f}d -> job {jid}\n")
        else:
            f.write("Not submitted. To submit:\n")
            for epoch_day in epochs:
                sd = base_dir / f"{epoch_day:.0f}d"
                f.write(f"  sbatch {sd}/run_artis.slurm\n")

    print(f"\nSummary written to {summary_path}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
