#!/usr/bin/env python3
"""
Plot STELLA model evolution timeseries from .swd file.

This script creates a multi-panel plot showing the evolution of key
physical quantities over time from the STELLA .swd output.

Usage:
    python plot_stella_evolution.py <swd_file> [output_dir]
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
from pathlib import Path
import sys
import argparse


def read_swd_epochs(swd_path: str) -> dict:
    """Read STELLA .swd file and extract epoch data.
    
    STELLA .swd format:
    - Epoch header line: time (days) in first column, zone number in second column
    - Zone data lines: 0 in first column, zone number in second column
    
    Note: Ni56 evolution is not included in .swd output. Ni56 is handled
    internally by STELLA via the .xni file and is not output to .swd.
    
    Returns dict with keys:
        - times: array of epoch times (days)
        - luminosity: array of luminosities (Lsun)
        - T_gas: array of gas temperatures (K) at photosphere
        - T_rad: array of radiation temperatures (K) at photosphere
        - velocity: array of velocities (km/s) at photosphere
        - radius: array of radii (cm) at photosphere
        - density: array of densities (g/cm^3) at photosphere
    """
    epochs = {}
    current_t = None
    buf_data = []
    line_count = 0
    
    with open(swd_path) as f:
        for line in f:
            line_count += 1
            cols = line.split()
            
            if len(cols) < 13:
                continue
            
            # First column is time (days) or 0 for zone data
            try:
                t = float(cols[0])
            except ValueError:
                continue
            
            # New epoch: first column is not 0 (epoch header line)
            if t != 0:
                if current_t is not None and buf_data:
                    # Process previous epoch
                    arr = np.array(buf_data)
                    epochs[current_t] = arr
                current_t = t
                buf_data = []
                # Add epoch header line data (zone 1)
                buf_data.append([float(c) for c in cols[:13]])
            # Zone data line (first column is 0)
            elif t == 0 and current_t is not None:
                # Store data: t, km, log10(M), log10(R), V, log10(Tg), log10(Tr), log10(rho), ...
                buf_data.append([float(c) for c in cols[:13]])
    
    # Don't forget last epoch
    if current_t is not None and buf_data:
        epochs[current_t] = np.array(buf_data)
    
    # Extract photospheric values for each epoch
    times = sorted(t for t in epochs.keys())
    luminosity = []
    T_gas = []
    T_rad = []
    velocity = []
    radius = []
    density = []
    
    for t in times:
        data = epochs[t]
        # Columns: t, km, log10(M), log10(R), V, log10(Tg), log10(Tr), log10(rho), log10(P), log10(Qv), log10(eng), Flum
        # Use outer zone (index 0) for photospheric values
        # Flum is in 1e40 erg/s, convert to Lsun
        Flum = data[0, 11]  # Flum in 1e40 erg/s (column 11)
        L_lsun = Flum * 1e40 / 3.828e33
        
        luminosity.append(L_lsun)
        T_gas.append(10**data[0, 5])  # log10(T_gas) - column 5
        T_rad.append(10**data[0, 6])  # log10(T_rad) - column 6
        velocity.append(data[0, 4]*1e3)   # V in 1000 km/s - column 4
        radius.append(10**data[0, 3])  # log10(R) in cm - column 3
        density.append(10**data[0, 7] * 1e-6)  # log10(rho) in 1e-6 g/cc - column 7
    
    return {
        'times': np.array(times),
        'luminosity': np.array(luminosity),
        'T_gas': np.array(T_gas),
        'T_rad': np.array(T_rad),
        'velocity': np.array(velocity),
        'radius': np.array(radius),
        'density': np.array(density),
    }


def read_swd_epochs_full(swd_path: str) -> dict:
    """Read STELLA .swd file and return full epoch data for radial profiles.
    
    Returns dict keyed by epoch time -> numpy array of zone data.
    """
    epochs = {}
    current_t = None
    buf_data = []
    line_count = 0
    
    with open(swd_path) as f:
        for line in f:
            line_count += 1
            cols = line.split()
            
            if len(cols) < 13:
                continue
            
            # First column is time (days) or 0 for zone data
            try:
                t = float(cols[0])
            except ValueError:
                continue
            
            # New epoch: first column is not 0 (epoch header line)
            if t != 0:
                if current_t is not None and buf_data:
                    # Process previous epoch
                    arr = np.array(buf_data)
                    # Only keep if we don't have this epoch yet, or if this has more zones
                    if current_t not in epochs or len(arr) > len(epochs[current_t]):
                        epochs[current_t] = arr
                current_t = t
                buf_data = []
                # Add epoch header line data (zone 1)
                buf_data.append([float(c) for c in cols[:13]])
            # Zone data line (first column is 0)
            elif t == 0 and current_t is not None:
                # Store data: t, km, log10(M), log10(R), V, log10(Tg), log10(Tr), log10(rho), ...
                buf_data.append([float(c) for c in cols[:13]])
    
    # Don't forget last epoch
    if current_t is not None and buf_data:
        arr = np.array(buf_data)
        if current_t not in epochs or len(arr) > len(epochs[current_t]):
            epochs[current_t] = arr
    
    print(f"DEBUG: Total epochs read: {len(epochs)}")
    for t in sorted(epochs.keys())[:5]:
        print(f"DEBUG: Epoch {t}: {len(epochs[t])} zones")
    return epochs


def plot_evolution(data, output_file: str, pre_sbo_only: bool = False):
    """Create multi-panel evolution plot."""
    # Use original times (no offset) to preserve negative values
    times = data['times']
    
    # Separate pre-SBO (negative) and post-SBO (positive) times
    pre_sbo_mask = times < 0
    post_sbo_mask = times >= 0
    
    # If only plotting pre-SBO, filter data
    if pre_sbo_only:
        # Filter to only negative times (pre-explosion)
        pre_sbo_mask = data['times'] < 0
        post_sbo_mask = np.zeros_like(data['times'], dtype=bool)
        # Filter all data arrays
        times = data['times'][pre_sbo_mask]
        data['luminosity'] = data['luminosity'][pre_sbo_mask]
        data['T_gas'] = data['T_gas'][pre_sbo_mask]
        data['T_rad'] = data['T_rad'][pre_sbo_mask]
        data['velocity'] = data['velocity'][pre_sbo_mask]
        data['radius'] = data['radius'][pre_sbo_mask]
        data['density'] = data['density'][pre_sbo_mask]
        # Update masks
        pre_sbo_mask = np.ones_like(times, dtype=bool)
        post_sbo_mask = np.zeros_like(times, dtype=bool)
        title = 'STELLA Pre-Shock-Breakout Evolution'
    else:
        # Normal mode: use all times, keep negative values
        title = 'STELLA Model Evolution'
    
    fig, axes = plt.subplots(3, 2, figsize=(14, 10))
    fig.suptitle(title, fontsize=16, fontweight='bold')
    
    # Luminosity
    ax = axes[0, 0]
    if np.any(post_sbo_mask):
        ax.plot(times[post_sbo_mask], data['luminosity'][post_sbo_mask], 'b-', linewidth=2, label='post-SBO')
    if np.any(pre_sbo_mask):
        ax.plot(times[pre_sbo_mask], data['luminosity'][pre_sbo_mask], 'b--', linewidth=2, alpha=0.7, label='pre-SBO')
    ax.set_xlabel('Time since explosion (days)')
    ax.set_ylabel(r'Luminosity (L$_\odot$)')
    ax.set_yscale('log')
    ax.grid(True, alpha=0.3)
    ax.set_title('Luminosity')
    ax.legend()
    
    # Temperature
    ax = axes[0, 1]
    if np.any(post_sbo_mask):
        ax.plot(times[post_sbo_mask], data['T_gas'][post_sbo_mask], 'r-', label='T_gas (post-SBO)', linewidth=2)
        ax.plot(times[post_sbo_mask], data['T_rad'][post_sbo_mask], 'b-', label='T_rad (post-SBO)', linewidth=2)
    if np.any(pre_sbo_mask):
        ax.plot(times[pre_sbo_mask], data['T_gas'][pre_sbo_mask], 'r--', label='T_gas (pre-SBO)', linewidth=2, alpha=0.7)
        ax.plot(times[pre_sbo_mask], data['T_rad'][pre_sbo_mask], 'b--', label='T_rad (pre-SBO)', linewidth=2, alpha=0.7)
    ax.set_xlabel('Time since explosion (days)')
    ax.set_ylabel('Temperature (K)')
    ax.set_yscale('log')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_title('Temperature')
    
    # Velocity
    ax = axes[1, 0]
    if np.any(post_sbo_mask):
        ax.plot(times[post_sbo_mask], data['velocity'][post_sbo_mask], 'g-', linewidth=2, label='post-SBO')
    if np.any(pre_sbo_mask):
        ax.plot(times[pre_sbo_mask], data['velocity'][pre_sbo_mask], 'g--', linewidth=2, alpha=0.7, label='pre-SBO')
    ax.set_xlabel('Time since explosion (days)')
    ax.set_ylabel(r'Velocity ($10^4$ km s$^{-1}$)')
    ax.grid(True, alpha=0.3)
    ax.set_title('Photospheric Velocity')
    ax.legend()
    
    # Radius
    ax = axes[1, 1]
    if np.any(post_sbo_mask):
        ax.plot(times[post_sbo_mask], data['radius'][post_sbo_mask] / 1e15, 'm-', linewidth=2, label='post-SBO')
    if np.any(pre_sbo_mask):
        ax.plot(times[pre_sbo_mask], data['radius'][pre_sbo_mask] / 1e15, 'm--', linewidth=2, alpha=0.7, label='pre-SBO')
    ax.set_xlabel('Time since explosion (days)')
    ax.set_ylabel('Radius (10$^{15}$ cm)')
    ax.grid(True, alpha=0.3)
    ax.set_title('Photospheric Radius')
    ax.legend()
    
    # Density
    ax = axes[2, 0]
    if np.any(post_sbo_mask):
        ax.plot(times[post_sbo_mask], data['density'][post_sbo_mask], 'c-', linewidth=2, label='post-SBO')
    if np.any(pre_sbo_mask):
        ax.plot(times[pre_sbo_mask], data['density'][pre_sbo_mask], 'c--', linewidth=2, alpha=0.7, label='pre-SBO')
    ax.set_xlabel('Time since explosion (days)')
    ax.set_ylabel('Density (g/cm$^3$)')
    ax.set_yscale('log')
    ax.grid(True, alpha=0.3)
    ax.set_title('Photospheric Density')
    ax.legend()
    
    # Luminosity vs Velocity (Hertzsprung-Russell-like)
    ax = axes[2, 1]
    sc = ax.scatter(data['velocity'], data['luminosity'], c=times, 
                    cmap='viridis', s=50, alpha=0.7)
    ax.set_xlabel(r'Velocity ($10^4$ km s$^{-1}$)')
    ax.set_ylabel(r'Luminosity (L$_\odot$)')
    # Use log scale only if data has positive values
    if np.any(data['luminosity'] > 0):
        ax.set_yscale('log')
    if np.any(data['velocity'] > 0):
        ax.set_xscale('log')
    ax.grid(True, alpha=0.3)
    ax.set_title('L-V Diagram')
    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label('Time since explosion (days)')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Evolution plot saved to {output_file}")
    plt.close()


def plot_radial_profiles(epochs, output_file: str, n_epochs_max: int = 10, pre_sbo_only: bool = False):
    """Create radial profile plots for each epoch.
    
    Parameters:
    -----------
    epochs : dict
        Dictionary of epoch data from read_swd_epochs
    output_file : str
        Output file path
    n_epochs_max : int
        Maximum number of epochs to plot (evenly spaced in time)
    pre_sbo_only : bool
        Plot only pre-shock-breakout epochs
    """
    # Include all epochs
    times = sorted(t for t in epochs.keys())
    t_offset = times[0]  # Offset to make t=0 at energy injection
    
    # Filter for pre-SBO only if requested
    if pre_sbo_only:
        # Filter to only negative times (pre-explosion)
        times = [t for t in times if t < 0]

    # Select epochs to plot (evenly spaced)
    if len(times) > n_epochs_max:
        indices = np.linspace(0, len(times) - 1, n_epochs_max, dtype=int)
        selected_times = [times[i] for i in indices]
    else:
        selected_times = times
    
    # Create colormap for epochs
    cmap = plt.cm.viridis
    colors = cmap(np.linspace(0, 1, len(selected_times)))
    
    # Update title for pre-SBO only
    if pre_sbo_only:
        title = 'STELLA Pre-Shock-Breakout Radial Profiles'
    else:
        title = 'STELLA Radial Profiles by Epoch'
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(title, fontsize=16, fontweight='bold')
    
    # Temperature vs radius
    ax = axes[0, 0]
    for i, t in enumerate(selected_times):
        data = epochs[t]
        # Columns: t, zone, log10(M), log10(R), V, log10(Tg), log10(Tr), log10(rho), log10(P), log10(Qv), log10(eng), Flum, WRKX
        r_cm = 10**data[:, 3]  # log10(R) in cm
        T_gas = 10**data[:, 5]  # log10(T_gas)
        # Use dashed line for pre-SBO (negative time), solid for post-SBO
        # If pre_sbo_only, all lines are solid with full opacity
        if pre_sbo_only:
            ls = '-'
            alpha = 0.8
            label_suffix = ''
        else:
            ls = '--' if t < 0 else '-'
            alpha = 0.7 if t < 0 else 0.8
            label_suffix = ' (pre-SBO)' if t < 0 else ''
        ax.plot(r_cm / 1e15, T_gas, color=colors[i], 
                label=f't={t:.2f}d{label_suffix}', linewidth=1.5, alpha=alpha, ls=ls)
    ax.set_xlabel('Radius (10$^{15}$ cm)')
    ax.set_ylabel('T_gas (K)')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.grid(True, alpha=0.3)
    ax.set_title('Temperature Profile')
    ax.legend(ncol=2, fontsize=8)
    
    # Density vs radius
    ax = axes[0, 1]
    for i, t in enumerate(selected_times):
        data = epochs[t]
        r_cm = 10**data[:, 3]  # log10(R) in cm
        rho = 10**data[:, 7] * 1e-6  # log10(rho) in 1e-6 g/cc
        # Use same line style logic as temperature
        if pre_sbo_only:
            ls = '-'
            alpha = 0.8
        else:
            ls = '--' if t < 0 else '-'
            alpha = 0.7 if t < 0 else 0.8
        ax.plot(r_cm / 1e15, rho, color=colors[i], 
                label=f't={t:.2f}d', linewidth=1.5, alpha=alpha, ls=ls)
    ax.set_xlabel('Radius (10$^{15}$ cm)')
    ax.set_ylabel('Density (g/cm$^3$)')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.grid(True, alpha=0.3)
    ax.set_title('Density Profile')
    
    # Velocity vs radius
    ax = axes[1, 0]
    for i, t in enumerate(selected_times):
        data = epochs[t]
        r_cm = 10**data[:, 3]  # log10(R) in cm
        v_kms = data[:, 4] * 1e3  # V in 1000 km/s, convert to km/s
        # Use same line style logic as temperature
        if pre_sbo_only:
            ls = '-'
            alpha = 0.8
        else:
            ls = '--' if t < 0 else '-'
            alpha = 0.7 if t < 0 else 0.8
        ax.plot(r_cm / 1e15, v_kms, color=colors[i], 
                label=f't={t:.2f}d', linewidth=1.5, alpha=alpha, ls=ls)
    ax.set_xlabel('Radius (10$^{15}$ cm)')
    ax.set_ylabel('Velocity (km/s)')
    ax.set_xscale('log')
    ax.grid(True, alpha=0.3)
    ax.set_title('Velocity Profile')
    
    # Flum (luminosity flux) vs radius
    ax = axes[1, 1]
    for i, t in enumerate(selected_times):
        data = epochs[t]
        r_cm = 10**data[:, 3]  # log10(R) in cm
        flum = data[:, 11]  # Flum in 1e40 erg/s
        # Use same line style logic as temperature
        if pre_sbo_only:
            ls = '-'
            alpha = 0.8
        else:
            ls = '--' if t < 0 else '-'
            alpha = 0.7 if t < 0 else 0.8
        ax.plot(r_cm / 1e15, flum, color=colors[i], 
                label=f't={t:.2f}d', linewidth=1.5, alpha=alpha, ls=ls)
    ax.set_xlabel('Radius (10$^{15}$ cm)')
    ax.set_ylabel('Flum (10$^{40}$ erg s$^{-1}$)')
    ax.set_xscale('log')
    ax.grid(True, alpha=0.3)
    ax.set_title('Luminosity Flux Profile')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Radial profile plot saved to {output_file}")
    plt.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('swd_file', help='STELLA .swd file path')
    parser.add_argument('output_dir', nargs='?', default='.',
                        help='Output directory (default: current directory)')
    parser.add_argument('--n-epochs', type=int, default=10,
                        help='Number of epochs to plot in radial profiles (default: 10)')
    parser.add_argument('--pre-sbo-only', action='store_true',
                        help='Plot only pre-shock-breakout (negative time) epochs')
    args = parser.parse_args()
    
    swd_path = Path(args.swd_file)
    if not swd_path.exists():
        print(f"ERROR: {swd_path} not found")
        sys.exit(1)
    
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Reading {swd_path}...")
    epochs = read_swd_epochs_full(str(swd_path))
    times = sorted(epochs.keys())
    print(f"Found {len(times)} epochs from {times[0]:.2f} to {times[-1]:.2f} days")
    
    # Plot evolution timeseries
    data = read_swd_epochs(str(swd_path))
    output_file = output_dir / f"{swd_path.stem}_evolution.png"
    if args.pre_sbo_only:
        output_file = output_dir / f"{swd_path.stem}_pre_sbo_evolution.png"
    plot_evolution(data, str(output_file), pre_sbo_only=args.pre_sbo_only)
    
    # Plot radial profiles for each epoch
    output_file_radial = output_dir / f"{swd_path.stem}_radial_profiles.png"
    if args.pre_sbo_only:
        output_file_radial = output_dir / f"{swd_path.stem}_pre_sbo_radial_profiles.png"
    plot_radial_profiles(epochs, str(output_file_radial), n_epochs_max=args.n_epochs, pre_sbo_only=args.pre_sbo_only)


if __name__ == '__main__':
    main()
