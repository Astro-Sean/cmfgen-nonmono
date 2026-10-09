#!/usr/bin/env python3
"""Extract a snapshot from a STELLA .swd file and create STELLA .hyd + .abn input files.

This allows restarting STELLA from a previous run's snapshot without needing
the original .prf/.crv restart files. The extracted .hyd/.abn files can be
used as input for a new STELLA run that continues evolution for a few more days.

Usage:
  python extract_swd_snapshot_to_stella.py SWD_FILE OUTPUT_DIR [OPTIONS]
    --target-days T    Extract epoch closest to T days (default: last epoch)
    --model-name NAME  Model name for output files (default: from swd filename)
"""

import sys
import os
import numpy as np
from pathlib import Path

# Use the existing SWD reader from the conversion script
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from convert_stella_swd_to_cmfgen import read_swd


def write_stella_hyd(output_path: str, t_start_s: float, nzon: int,
                     km: np.ndarray, r: np.ndarray, rho: np.ndarray,
                     T: np.ndarray, v: np.ndarray, mass_cut: float = 0.0):
    """Write STELLA .hyd file from snapshot data.

    STELLA .hyd format:
      Line 1: timeStart(s)  Nzon  mass_cut(Msun)  Rcen(cm)  rhoCen(g/cc)
      Lines 2..Nzon+1: km  dMr(Msun)  r(cm)  rho(g/cc)  T(K)  u(cm/s)  Mr(Msun)  0

    Input arrays are outer-to-inner (decreasing r) from read_swd.
    We reverse to inner-to-outer (increasing r) for STELLA .hyd.
    """
    # Reverse to inner-to-outer (increasing r)
    # read_swd returns v in km/s; STELLA .hyd needs cm/s
    km_arr = km[::-1]
    r_arr = r[::-1]
    rho_arr = rho[::-1]
    T_arr = T[::-1]
    v_arr = v[::-1] * 1.0e5  # km/s → cm/s

    # Zone masses from cell-edge volumes: geometric-midpoint edges,
    # inner/outer edges extrapolated at the local spacing.  A plain
    # trapezoid leaves zone 1 at zero mass, which strad rejects.
    r_edge = np.empty(nzon + 1)
    r_edge[1:-1] = np.sqrt(r_arr[:-1] * r_arr[1:])
    r_edge[0] = r_arr[0] * (r_arr[0] / r_arr[1]) ** 0.5
    r_edge[-1] = r_arr[-1] * (r_arr[-1] / r_arr[-2]) ** 0.5
    dm_g = (4.0 / 3.0) * np.pi * (r_edge[1:]**3 - r_edge[:-1]**3) * rho_arr

    dMr = dm_g / 1.989e33
    Mr_msun = np.cumsum(dMr)

    # Central values
    r_cen = r_arr[0]
    rho_cen = rho_arr[0]

    with open(output_path, 'w') as f:
        f.write(f"{t_start_s:.6e}  {nzon}  {mass_cut:.6f}  {r_cen:.6e}  {rho_cen:.6e}\n")
        for i in range(nzon):
            f.write(f"{km_arr[i]:.6e}  {dMr[i]:.6e}  {r_arr[i]:.6e}  "
                    f"{rho_arr[i]:.6e}  {T_arr[i]:.6e}  {v_arr[i]:.6e}  "
                    f"{Mr_msun[i]:.6e}  0\n")


def write_stella_abn(output_path: str, nzon: int, km: np.ndarray,
                      composition: np.ndarray):
    """Write STELLA .abn file from snapshot composition.

    STELLA .abn format:
      Lines 1..Nzon: km  0  0  0  H  He  C  N  O  Ne  Na  Mg  Al  Si  S  Ar  Ca  FePeak  Ni58  Ni56

    Input arrays are outer-to-inner (decreasing r) from read_swd.
    We reverse to inner-to-outer to match .hyd file.
    """
    # Reverse to inner-to-outer (increasing r)
    km_arr = km[::-1]
    comp_arr = composition[::-1]

    with open(output_path, 'w') as f:
        for i in range(nzon):
            vals = [f"{km_arr[i]:.6e}", "0", "0", "0"]
            for j in range(comp_arr.shape[1]):
                vals.append(f"{comp_arr[i, j]:.6e}")
            # .swd composition (15 cols) ends with stable Ni, which occupies
            # the .abn Ni58 slot. Only Ni56 (radioactive) is missing.
            vals.append("0.000000e+00")  # Ni56
            f.write("  ".join(vals) + "\n")


def main():
    if len(sys.argv) < 3:
        print("Usage: python extract_swd_snapshot_to_stella.py SWD_FILE OUTPUT_DIR [OPTIONS]")
        print("Options:")
        print("  --target-days T    Extract epoch closest to T days (default: last)")
        print("  --model-name NAME  Model name (default: from swd filename)")
        sys.exit(1)

    swd_file = sys.argv[1]
    output_dir = Path(sys.argv[2])
    output_dir.mkdir(parents=True, exist_ok=True)

    target_days = None
    model_name = None

    for arg in sys.argv[3:]:
        if arg.startswith('--target-days='):
            target_days = float(arg.split('=')[1])
        elif arg.startswith('--model-name='):
            model_name = arg.split('=')[1]

    if model_name is None:
        model_name = Path(swd_file).stem

    print(f"Reading SWD file: {swd_file}")
    epochs = read_swd(swd_file)
    times = sorted(epochs.keys())
    print(f"Found {len(times)} epochs")
    print(f"Time range: {times[0]:.4f} - {times[-1]:.4f} days")

    # Select epoch
    if target_days is not None:
        t_arr = np.array(times)
        idx = np.argmin(np.abs(t_arr - target_days))
        t_selected = times[idx]
        print(f"Selected epoch: t={t_selected:.4f} d (target={target_days} d)")
    else:
        t_selected = times[-1]
        print(f"Selected last epoch: t={t_selected:.4f} d")

    epoch = epochs[t_selected]
    nzon = len(epoch['r'])

    # Write .hyd and .abn
    hyd_path = output_dir / f"{model_name}.hyd"
    abn_path = output_dir / f"{model_name}.abn"

    # .swd epoch times are observer-frame (model time minus the light-
    # travel time across r_max); .hyd timeStart is model-frame.  Add the
    # light-travel time back so the restart begins at the right epoch -
    # the initial-condition epoch (-14 d obs) maps to ~0 model time,
    # matching the +5 s timeStart the original runs used.
    t_lt_s = float(epoch['r'].max()) / 2.998e10
    t_start_s = t_selected * 86400.0 + t_lt_s

    print(f"Writing {hyd_path}...")
    write_stella_hyd(str(hyd_path), t_start_s, nzon,
                     epoch['km'], epoch['r'], epoch['rho'],
                     epoch['T'], epoch['v'])

    print(f"Writing {abn_path}...")
    write_stella_abn(str(abn_path), nzon, epoch['km'], epoch['composition'])

    print(f"Done. NZON={nzon}, t={t_selected:.4f} d")
    print(f"  R range: {epoch['r'].min():.3e} - {epoch['r'].max():.3e} cm")
    print(f"  V range: {epoch['v'].min():.3e} - {epoch['v'].max():.3e} cm/s")
    print(f"  T range: {epoch['T'].min():.3e} - {epoch['T'].max():.3e} K")
    print(f"  rho range: {epoch['rho'].min():.3e} - {epoch['rho'].max():.3e} g/cc")


if __name__ == '__main__':
    main()
