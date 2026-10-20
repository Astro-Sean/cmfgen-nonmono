#!/usr/bin/env python3
"""
Homologously evolve a STELLA .swd snapshot to a later epoch.

In homologous expansion (free expansion phase), each fluid element moves at
constant velocity, so:
    r(t) = r(t0) * (t / t0)
    v(t) = v(t0)                (velocity is constant)
    rho(t) = rho(t0) * (t0 / t)^3   (mass conservation in spherical geometry)
    T(t) = T(t0) * (t0 / t)     (adiabatic expansion, radiation-dominated: T ∝ R^-1)

This script reads a STELLA .swd file, finds the epoch closest to t_source_days,
scales all quantities to t_target_days, and writes a new .swd file (or appends
to an existing one) with the evolved snapshot.

This is useful for extending STELLA snapshots to late times (e.g. 349 days)
without running the full radiation hydrodynamics simulation to that epoch.

Usage:
    python homologous_evolution.py <input.swd> <output.swd> \
        --source-day 60 --target-day 349
"""

import argparse
import sys
from pathlib import Path
import numpy as np

# Import the read_swd function from the conversion script
sys.path.insert(0, str(Path(__file__).parent))
from convert_stella_swd_to_cmfgen import read_swd


def find_closest_epoch(epochs: dict, target_day: float, tolerance: float = 5.0):
    """Find the epoch closest to target_day in the .swd file.

    Returns (t_day, epoch_data) or raises ValueError if no epoch within tolerance.
    """
    available = sorted(epochs.keys())
    closest_t = min(available, key=lambda t: abs(t - target_day))
    diff = abs(closest_t - target_day)
    if diff > tolerance:
        raise ValueError(
            f"No epoch within {tolerance} days of {target_day}. "
            f"Closest: {closest_t} (diff={diff:.2f} days). "
            f"Available epochs: {available[:5]}...{available[-5:]}"
        )
    print(f"  Source epoch: {closest_t} days (requested {target_day}, diff={diff:.2f} days)")
    return closest_t, epochs[closest_t]


def homologous_evolve(epoch_data: dict, t_source: float, t_target: float) -> dict:
    """Scale a STELLA snapshot from t_source to t_target using homologous expansion.

    Parameters:
    -----------
    epoch_data : dict
        Dictionary with keys 'r', 'v', 'rho', 'T', 'flux', 'km', etc.
    t_source : float
        Source epoch in days
    t_target : float
        Target epoch in days

    Returns:
    --------
    dict : Scaled epoch data with same keys
    """
    if t_target <= 0:
        raise ValueError(f"Target day ({t_target}) must be > 0")
    if t_target == t_source:
        return dict(epoch_data)

    # Scale factor for homologous expansion
    scale = t_target / t_source
    scale_rho = scale ** (-3)  # density drops as t^-3 (volume ∝ r^3 ∝ t^3)
    scale_T = scale ** (-1)    # adiabatic: T ∝ R^-1 for radiation-dominated

    print(f"  Homologous expansion: {t_source} -> {t_target} days (scale={scale:.4f})")
    print(f"  rho scale: {scale_rho:.6e}, T scale: {scale_T:.6f}")

    evolved = {}
    for key in epoch_data:
        val = epoch_data[key]
        if val is None:
            evolved[key] = None
            continue

        if key == 'r':
            evolved[key] = val * scale
        elif key == 'v':
            # Velocity is constant in homologous expansion
            evolved[key] = val.copy() if hasattr(val, 'copy') else val
        elif key == 'rho':
            evolved[key] = val * scale_rho
        elif key == 'T':
            evolved[key] = val * scale_T
        elif key == 'flux':
            # Luminosity decreases in homologous expansion
            # For radioactive decay: L ∝ M_Ni * exp(-t/tau_Ni) * (1 - exp(-t/tau_Co))
            # For adiabatic: L ∝ T^4 * R^2 ∝ (R^-1)^4 * R^2 = R^-2 ∝ t^-2
            # Use t^-2 scaling as approximation
            evolved[key] = val * (scale ** (-2))
        elif key == 'km':
            evolved[key] = val.copy() if hasattr(val, 'copy') else val
        elif key == 'kappa':
            # Opacity is approximately constant (electron scattering dominated)
            evolved[key] = val.copy() if hasattr(val, 'copy') else val
        elif key == 'nbar':
            # Baryon number density ∝ rho
            evolved[key] = val * scale_rho
        elif key == 'ne':
            # Electron number density ∝ rho (assuming ionization fraction ~const)
            evolved[key] = val * scale_rho
        elif key == 'composition':
            # Composition is unchanged
            evolved[key] = val.copy() if hasattr(val, 'copy') else val
        else:
            # Default: copy unchanged
            evolved[key] = val.copy() if hasattr(val, 'copy') else val

    return evolved


def write_swd_epoch(fh, t_day: float, epoch_data: dict):
    """Write a single epoch to a .swd file in STELLA Format 181.

    The .swd format (from lbalsw.f) for each zone line:
        t_day  km  ...  log_r  v_1e8  log_T  ...  log_rho  ...  flux_1e40  kappa  nbar  ne  [composition]

    Column layout (0-indexed):
        0: t_day (nonzero for first line of epoch, 0 for subsequent)
        1: km (zone number, 1-based)
        2: (unused)
        3: log10(r_cm)
        4: v in 1e8 cm/s
        5: log10(T_K)
        6: (unused)
        7: log10(rho / 1e-6 g/cc)
        8-10: (unused)
        11: flux in 1e40 erg/s
        12: kappa (Rosseland opacity, cm^2/g)
        13: nbar (baryon number density, cm^-3)
        14: ne (electron number density, cm^-3)
        15-29: composition (15 elements: H, He, C, N, O, Ne, Na, Mg, Al, Si, S, Ar, Ca, Fe, Ni)
    """
    r = epoch_data['r']
    v = epoch_data['v']           # km/s
    T = epoch_data['T']
    rho = epoch_data['rho']
    flux = epoch_data['flux']
    km_arr = epoch_data['km']
    kappa = epoch_data.get('kappa')
    nbar = epoch_data.get('nbar')
    ne = epoch_data.get('ne')
    comp = epoch_data.get('composition')

    n = len(r)
    # STELLA .swd is written inner->outer (km=1 to km=N)
    # But our data is outer->inner (sorted by decreasing r in read_swd)
    # Reverse to inner->outer for writing
    order = np.argsort(r)  # inner->outer

    for i in order:
        km = int(km_arr[i])
        t_val = t_day if km == 1 else 0.0

        log_r = np.log10(r[i])
        v_1e8 = v[i] * 1e5 / 1e8  # km/s -> 1e8 cm/s
        log_T = np.log10(T[i])
        # STELLA SWD: rho is log10(rho / 1e-6 g/cc)
        log_rho = np.log10(rho[i] / 1e-6)
        flux_1e40 = flux[i] / 1e40

        # Build the line
        parts = [
            f"{t_val:12.6f}",
            f"{km:6d}",
            f"{0.0:12.6f}",       # col 2 (unused)
            f"{log_r:12.7f}",     # col 3
            f"{v_1e8:10.4f}",     # col 4
            f"{log_T:10.4f}",     # col 5
            f"{0.0:8.3f}",        # col 6 (unused)
            f"{log_rho:10.4f}",   # col 7
            f"{0.0:8.3f}",        # col 8 (unused)
            f"{0.0:8.3f}",        # col 9 (unused)
            f"{0.0:8.3f}",        # col 10 (unused)
            f"{flux_1e40:12.5E}", # col 11
        ]

        if kappa is not None:
            parts.append(f"{kappa[i]:12.5E}")  # col 12
        if nbar is not None:
            parts.append(f"{nbar[i]:12.5E}")   # col 13
        if ne is not None:
            parts.append(f"{ne[i]:12.5E}")     # col 14
        if comp is not None:
            for j in range(15):
                parts.append(f"{comp[i][j]:16.9E}")

        fh.write(' '.join(parts) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                 formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('input_swd', help='Input STELLA .swd file')
    parser.add_argument('output_swd', help='Output .swd file (will be created/overwritten)')
    parser.add_argument('--source-day', type=float, required=True,
                        help='Source epoch in days (closest available epoch will be used)')
    parser.add_argument('--target-day', type=float, required=True,
                        help='Target epoch in days (homologously evolved)')
    parser.add_argument('--tolerance', type=float, default=10.0,
                        help='Max |source-day - closest_epoch| tolerance (default: 10 days)')
    parser.add_argument('--append', action='store_true',
                        help='Append to existing output .swd file instead of overwriting')
    args = parser.parse_args()

    print(f"\n=== Homologous Evolution: {args.source_day} -> {args.target_day} days ===")
    print(f"  Input: {args.input_swd}")
    print(f"  Output: {args.output_swd}")

    # Read the source .swd file
    print(f"\n  Reading {args.input_swd}...")
    epochs = read_swd(args.input_swd)
    print(f"  Found {len(epochs)} epochs: {sorted(epochs.keys())[:5]}...{sorted(epochs.keys())[-5:]}")

    # Find the closest epoch to source-day
    t_source, epoch_data = find_closest_epoch(epochs, args.source_day, args.tolerance)

    # Homologously evolve
    print(f"\n  Evolving...")
    evolved = homologous_evolve(epoch_data, t_source, args.target_day)

    # Write the output .swd file
    mode = 'a' if args.append else 'w'
    print(f"\n  Writing evolved epoch to {args.output_swd} (mode={mode})...")
    with open(args.output_swd, mode) as fh:
        write_swd_epoch(fh, args.target_day, evolved)

    # Print summary
    r = evolved['r']
    v = evolved['v']
    rho = evolved['rho']
    T = evolved['T']
    print(f"\n  Evolved snapshot at t={args.target_day} days:")
    print(f"    r:   {r.min():.3e} - {r.max():.3e} cm")
    print(f"    v:   {v.min():.1f} - {v.max():.1f} km/s")
    print(f"    rho: {rho.min():.3e} - {rho.max():.3e} g/cc")
    print(f"    T:   {T.min():.3e} - {T.max():.3e} K")
    print(f"\n=== Done ===")


if __name__ == '__main__':
    main()
