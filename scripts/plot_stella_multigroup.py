#!/usr/bin/env python3
"""
Plot STELLA multi-group luminosity evolution at a specific timestep.

STELLA .ph file contains frequency-dependent data including luminosity per frequency group.
This script plots the luminosity for each frequency group at a given timestep.

Usage:
    python plot_stella_multigroup.py <ph_file> --time <days> --output <png>
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
from pathlib import Path
import sys
import argparse


def read_ph_file(ph_path: str, target_time: float = None):
    """Read STELLA .ph file and extract data for a specific epoch or all epochs.
    
    .ph file format:
    - First line: frequency grid (in log frequency)
    - Subsequent lines: time, zone, value, then frequency group data (log luminosity)
    
    Returns tuple: (frequency_grid, epochs_dict)
    """
    with open(ph_path) as f:
        lines = f.readlines()
    
    # First line is frequency grid
    freq_grid = np.array([float(x) for x in lines[0].split()])
    
    epochs = {}
    current_t = None
    buf_data = []
    
    # Process subsequent lines
    for line in lines[1:]:
        cols = line.split()
        if len(cols) < 4:
            continue
        
        try:
            t = float(cols[0])
            zone = int(float(cols[1]))
        except ValueError:
            continue
        
        # New epoch: first column changes
        if t != current_t:
            if current_t is not None and buf_data:
                epochs[current_t] = np.array(buf_data)
            current_t = t
            buf_data = []
        
        # Add data for this zone
        buf_data.append([float(c) for c in cols])
    
    # Don't forget last epoch
    if current_t is not None and buf_data:
        epochs[current_t] = np.array(buf_data)
    
    return freq_grid, epochs


def read_swd_file(swd_path: str, target_time: float = None):
    """Read STELLA .swd file to get radius information."""
    epochs = {}
    current_t = None
    buf_data = []
    
    with open(swd_path) as f:
        for line in f:
            cols = line.split()
            if len(cols) < 13:
                continue
            
            try:
                t = float(cols[0])
                zone = int(float(cols[1]))
            except ValueError:
                continue
            
            # New epoch starts with zone 1 having t != 0
            if t != 0 and zone == 1:
                if current_t is not None and buf_data:
                    epochs[current_t] = np.array(buf_data)
                current_t = t
                buf_data = [[float(c) for c in cols[:13]]]
            elif t != 0 and zone == 500 and current_t is not None:
                buf_data.append([float(c) for c in cols[:13]])
            elif t == 0 and current_t is not None:
                buf_data.append([float(c) for c in cols[:13]])
    
    if current_t is not None and buf_data:
        epochs[current_t] = np.array(buf_data)
    
    return epochs


def plot_multigroup_luminosity(ph_path: str, swd_path: str, target_time: float, output_file: str):
    """Plot multi-group luminosity vs radius at a specific timestep."""
    freq_grid, ph_epochs = read_ph_file(ph_path)
    swd_epochs = read_swd_file(swd_path)
    
    times = sorted(ph_epochs.keys())
    
    # Find closest epoch to target time
    if target_time is not None:
        idx = np.abs(np.array(times) - target_time).argmin()
        selected_time = times[idx]
        print(f"Selected epoch at t = {selected_time:.2f} days (target was {target_time:.2f} days)")
    else:
        selected_time = times[-1]  # Use last epoch
        print(f"Selected epoch at t = {selected_time:.2f} days (last epoch)")
    
    ph_data = ph_epochs[selected_time]
    
    # Find closest time in .swd file (may have different epoch times)
    swd_times = sorted(swd_epochs.keys())
    swd_idx = np.abs(np.array(swd_times) - selected_time).argmin()
    swd_time = swd_times[swd_idx]
    print(f"Using .swd epoch at t = {swd_time:.2f} days (closest to {selected_time:.2f} days)")
    swd_data = swd_epochs[swd_time]
    
    # Extract radius from .swd file (column 3 is logR in cm)
    # STELLA units: V in 1e8 cm/s, rho in log10(rho/1e-6 g/cc)
    r_cm = 10**swd_data[:, 3]  # Radius in cm
    r_15 = r_cm / 1e15  # Radius in 10^15 cm
    
    # Extract frequency group data from .ph file
    # Column 0: time, Column 1: zone, Column 2: something, Columns 3+: frequency group data
    n_groups = ph_data.shape[1] - 3
    n_zones_ph = ph_data.shape[0]
    print(f"Number of frequency groups: {n_groups}")
    print(f"Number of zones in .ph file: {n_zones_ph}")
    print(f"Number of zones in .swd file: {len(r_15)}")
    
    # The .ph file appears to have data for only one zone (likely the photosphere)
    # Extract that zone's data
    ph_zone = int(ph_data[0, 1])
    print(f".ph file zone: {ph_zone}")
    
    # Extract frequency group data for the single zone
    group_data = ph_data[0, 3:]  # Single zone, all frequency groups
    
    # Convert from log to linear (assuming log10)
    # Filter out negative values that might be placeholders
    group_data_linear = np.where(group_data > 0, 10**group_data, np.nan)
    
    # Get radius for this zone from .swd file
    if ph_zone <= len(r_15):
        r_zone_15 = r_15[ph_zone - 1]  # Zone numbers are 1-indexed
        print(f"Radius for zone {ph_zone}: {r_zone_15:.4f} 10^15 cm")
    else:
        print(f"Warning: Zone {ph_zone} not in .swd file (has {len(r_15)} zones)")
        r_zone_15 = r_15[len(r_15) // 2]  # Use middle zone as fallback
    
    # Plot luminosity vs frequency (since we only have one zone)
    # Convert frequency grid from log to linear
    freq_linear = 10**freq_grid
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Plot frequency spectrum
    mask = ~np.isnan(group_data_linear)
    if np.sum(mask) > 0:
        ax.loglog(freq_linear[mask], group_data_linear[mask], 'b-', linewidth=2, marker='o', markersize=4)
    else:
        print("ERROR: No valid data points to plot")
        return
    
    ax.set_xlabel('Frequency [Hz]', fontsize=12)
    ax.set_ylabel('Luminosity [erg/s]', fontsize=12)
    ax.set_title(f'STELLA Multi-Group Luminosity at Zone {ph_zone} (r = {r_zone_15:.2f} 10$^{{15}}$ cm)\nt = {selected_time:.2f} days', 
                 fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    
    # Add wavelength axis on top
    ax2 = ax.twiny()
    wavelength = 3e18 / freq_linear[mask]  # c / f in Angstroms
    ax2.set_xscale('log')
    ax2.set_xlim(wavelength[0], wavelength[-1])
    ax2.set_xlabel('Wavelength [Å]', fontsize=12)
    ax2.invert_xaxis()  # Wavelength decreases with increasing frequency
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Multi-group luminosity plot saved to: {output_file}")
    plt.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('ph_file', help='STELLA .ph file path')
    parser.add_argument('swd_file', help='STELLA .swd file path')
    parser.add_argument('--time', type=float, required=True, help='Target time in days')
    parser.add_argument('--output', default=None, help='Output figure path')
    args = parser.parse_args()
    
    ph_path = Path(args.ph_file)
    swd_path = Path(args.swd_file)
    if not ph_path.exists():
        print(f"ERROR: {ph_path} not found")
        sys.exit(1)
    if not swd_path.exists():
        print(f"ERROR: {swd_path} not found")
        sys.exit(1)
    
    # Default output path
    if args.output:
        output_file = args.output
    else:
        output_file = ph_path.parent / f"{ph_path.stem}_multigroup_t{args.time:.1f}d.png"
    
    plot_multigroup_luminosity(str(ph_path), str(swd_path), args.time, str(output_file))


if __name__ == '__main__':
    main()
