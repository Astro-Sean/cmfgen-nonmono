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
    - Zone 1 line: epoch time (days) in first column, zone number (1) in second column
    - Zones 2..N-1: 0 in first column, zone number in second column
    - Zone N (last): epoch time repeated in first column, zone number in second column
    Both zone 1 and zone N carry the epoch time; only zone=1 marks the start of a new epoch.
    
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
            
            try:
                t = float(cols[0])
                zone = int(cols[1])
                vals = [float(c) for c in cols[:13]]
            except ValueError:
                # FORTRAN FORMAT 181/281 can produce concatenated fields when
                # F14.7 fills all 14 chars and F9.4 starts immediately after.
                # Fall back to fixed-width parsing using the FORMAT field positions.
                try:
                    t_str = line[0:11].strip()
                    t = float(t_str) if t_str else 0.0
                    zone = int(line[11:16])
                    vals = [
                        float(line[0:11].strip() or '0'),
                        float(line[11:16].strip()),
                        float(line[16:28].strip()),
                        float(line[28:42].strip()),
                        float(line[42:51].strip()),
                        float(line[51:60].strip()),
                        float(line[60:69].strip()),
                        float(line[69:78].strip()),
                        float(line[78:87].strip()),
                        float(line[87:96].strip()),
                        float(line[96:105].strip()),
                        float(line[105:116].strip()),
                        float(line[116:127].strip()),
                    ]
                except (ValueError, IndexError):
                    continue
            
            # New epoch starts when zone == 1 (first zone always carries the epoch time)
            if zone == 1:
                if current_t is not None and buf_data:
                    # Process previous epoch
                    arr = np.array(buf_data)
                    epochs[current_t] = arr
                current_t = t
                buf_data = []
                buf_data.append(vals)
            elif current_t is not None:
                # All other zones (including the last zone which repeats the epoch time)
                buf_data.append(vals)
    
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
        # Columns: t, km, log10(M), log10(R), V(1e8 cm/s), log10(Tg), log10(Tr), log10(rho/1e-6), log10(P), log10(Qv), log10(eng), Flum, tau
        # Find photosphere: zone where tau (column 12) is closest to 2/3
        tau = data[:, 12]
        i_phot = np.argmin(np.abs(tau - 2.0/3.0))
        
        # Use photospheric zone for values
        # Flum is in 1e40 erg/s, convert to Lsun
        Flum = data[i_phot, 11]  # Flum in 1e40 erg/s (column 11)
        L_lsun = Flum * 1e40 / 3.828e33
        
        luminosity.append(L_lsun)
        T_gas.append(10**data[i_phot, 5])  # log10(T_gas) - column 5
        T_rad.append(10**data[i_phot, 6])  # log10(T_rad) - column 6
        velocity.append(data[i_phot, 4] * 1000.0)  # Convert from 1e8 cm/s to km/s - column 4
        radius.append(10**data[i_phot, 3])  # log10(R) in cm - column 3
        density.append((10**data[i_phot, 7]) * 1e-6)  # Convert from log10(rho/1e-6) to g/cm^3 - column 7
    
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

            try:
                t = float(cols[0])
                zone = int(cols[1])
                vals = [float(c) for c in cols[:13]]
            except ValueError:
                # FORTRAN FORMAT 181/281 can produce concatenated fields when
                # F14.7 fills all 14 chars and F9.4 starts immediately after.
                # Fall back to fixed-width parsing using the FORMAT field positions.
                try:
                    t_str = line[0:11].strip()
                    t = float(t_str) if t_str else 0.0
                    zone = int(line[11:16])
                    vals = [
                        float(line[0:11].strip() or '0'),
                        float(line[11:16].strip()),
                        float(line[16:28].strip()),
                        float(line[28:42].strip()),
                        float(line[42:51].strip()),
                        float(line[51:60].strip()),
                        float(line[60:69].strip()),
                        float(line[69:78].strip()),
                        float(line[78:87].strip()),
                        float(line[87:96].strip()),
                        float(line[96:105].strip()),
                        float(line[105:116].strip()),
                        float(line[116:127].strip()),
                    ]
                except (ValueError, IndexError):
                    continue

            # New epoch starts when zone == 1 (first zone always carries the epoch time)
            if zone == 1:
                if current_t is not None and buf_data:
                    arr = np.array(buf_data)
                    if current_t not in epochs or len(arr) > len(epochs[current_t]):
                        epochs[current_t] = arr
                current_t = t
                buf_data = []
                buf_data.append(vals)
            elif current_t is not None:
                # All other zones (including the last zone which repeats the epoch time)
                buf_data.append(vals)
    
    # Don't forget last epoch
    if current_t is not None and buf_data:
        arr = np.array(buf_data)
        if current_t not in epochs or len(arr) > len(epochs[current_t]):
            epochs[current_t] = arr
    
    print(f"Found {len(epochs)} epochs from {min(epochs):.2f} to {max(epochs):.2f} days")
    for t in sorted(epochs.keys())[:3]:
        print(f"  Epoch {t:.3f}d: {len(epochs[t])} zones")
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
    ax.set_ylabel(r'Velocity (km s$^{-1}$)')
    ax.grid(True, alpha=0.3)
    ax.set_title('Photospheric Velocity')
    ax.legend()
    
    # Radius
    ax = axes[1, 1]
    if np.any(post_sbo_mask):
        ax.plot(times[post_sbo_mask], data['radius'][post_sbo_mask], 'm-', linewidth=2, label='post-SBO')
    if np.any(pre_sbo_mask):
        ax.plot(times[pre_sbo_mask], data['radius'][pre_sbo_mask], 'm--', linewidth=2, alpha=0.7, label='pre-SBO')
    ax.set_xlabel('Time since explosion (days)')
    ax.set_ylabel('Radius (cm)')
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
    # Normalise ylim to post-SBO levels if available — pre-SBO density is
    # ~10^6 g/cc (unexploded stellar interior) and would compress the
    # post-SBO range (~10^-10 to 10^-15) into an invisible band.
    post_sbo_rho = data['density'][post_sbo_mask]
    if len(post_sbo_rho) > 0 and np.any(post_sbo_rho > 0):
        ax.set_ylim(post_sbo_rho.min() * 0.1, post_sbo_rho.max() * 10)
    ax.grid(True, alpha=0.3)
    ax.set_title('Photospheric Density')
    ax.legend()
    
    # Luminosity vs Velocity (Hertzsprung-Russell-like)
    ax = axes[2, 1]
    sc = ax.scatter(data['velocity'], data['luminosity'], c=times, 
                    cmap='viridis', s=50, alpha=0.7, edgecolors='black', linewidth=0.5)
    ax.set_xlabel(r'Velocity (km s$^{-1}$)')
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
    
    # Add time labels to avoid overlap and show evolution direction
    if len(times) > 0:
        # Label every few points to avoid overcrowding
        label_interval = max(1, len(times) // 8)  # Show ~8 labels max
        for i in range(0, len(times), label_interval):
            ax.annotate(f'{times[i]:.1f}d', 
                       (data['velocity'][i], data['luminosity'][i]),
                       xytext=(5, 5), textcoords='offset points',
                       fontsize=8, alpha=0.7,
                       bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.7))
    
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
    
    # Get initial (pre-explosion) radius from first epoch for top axis
    first_epoch = epochs[times[0]]
    r_initial = 10**first_epoch[:, 3]  # cm

    # Collect y-range across all selected epochs for proper ylims.
    # Use post-SBO epochs (t >= 0) only for ALL y-axis ranges — pre-SBO
    # values represent the unexploded stellar interior (T ~ 10^9 K,
    # rho ~ 10^6 g/cc, v ~ 0) and would compress the post-SBO ejecta
    # profiles into an invisible band.
    T_all, rho_all, v_all, flum_all = [], [], [], []
    r_post_sbo = []
    for t in selected_times:
        d = epochs[t]
        if t >= 0:  # post-SBO only for all ylims
            T_all.append(10**d[:, 5])
            rho_all.append(10**d[:, 7] * 1e-6)
            v_all.append(d[:, 4] * 1e3)
            flum_all.append(d[:, 11])
            r_post_sbo.append(10**d[:, 3])
    if T_all:
        T_all = np.concatenate(T_all)
        rho_all = np.concatenate(rho_all)
        v_all = np.concatenate(v_all)
        flum_all = np.concatenate(flum_all)
        r_post_sbo = np.concatenate(r_post_sbo)
    else:
        # Fallback: no post-SBO epochs, use all
        T_all = np.concatenate([10**epochs[t][:, 5] for t in selected_times])
        rho_all = np.concatenate([10**epochs[t][:, 7] * 1e-6 for t in selected_times])
        v_all = np.concatenate([epochs[t][:, 4] * 1e3 for t in selected_times])
        flum_all = np.concatenate([epochs[t][:, 11] for t in selected_times])
        r_post_sbo = np.concatenate([10**epochs[t][:, 3] for t in selected_times])

    # Shared x-range from post-SBO epochs only (pre-SBO star is compact,
    # ~10^10 cm; post-SBO ejecta expands to ~10^15 cm)
    x_min = r_post_sbo.min() * 0.8
    x_max = r_post_sbo.max() * 1.2

    # Helper: add top x-axis showing initial radius (independent xlim from bottom)
    def add_initial_axis(ax_bot, r_init_cm):
        ax_top = ax_bot.twiny()
        ax_top.set_xscale('log')
        ax_top.set_xlabel('Initial R (cm)', fontsize=9, color='gray')
        ax_top.set_xlim(r_init_cm.min() * 0.8, r_init_cm.max() * 1.2)
        ax_top.tick_params(axis='x', colors='gray', labelsize=8)
        from matplotlib.ticker import LogLocator, LogFormatterSciNotation
        ax_top.xaxis.set_major_locator(LogLocator(numticks=6))
        ax_top.xaxis.set_major_formatter(LogFormatterSciNotation())
        return ax_top

    # Temperature vs radius
    ax = axes[0, 0]
    for i, t in enumerate(selected_times):
        data = epochs[t]
        r_cm = 10**data[:, 3]  # log10(R) in cm
        T_gas = 10**data[:, 5]  # log10(T_gas)
        if pre_sbo_only:
            ls = '-'
            alpha = 0.8
            label_suffix = ''
        else:
            ls = '--' if t < 0 else '-'
            alpha = 0.7 if t < 0 else 0.8
            label_suffix = ' (pre-SBO)' if t < 0 else ''
        ax.plot(r_cm, T_gas, color=colors[i], 
                label=f't={t:.2f}d{label_suffix}', linewidth=1.5, alpha=alpha, ls=ls)
        # Mark photosphere (tau ~ 2/3)
        tau_vals = data[:, 12]
        valid = tau_vals > -49
        if np.any(valid):
            i_phot = np.argmin(np.abs(tau_vals[valid] - 2.0/3.0))
            idx_phot = np.where(valid)[0][i_phot]
            ax.plot(r_cm[idx_phot], T_gas[idx_phot], 'o', color=colors[i], markersize=6, 
                    markeredgecolor='white', markeredgewidth=0.5, zorder=5)
    ax.set_xlabel('Radius (cm)')
    ax.set_ylabel('T_gas (K)')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(T_all.min() * 0.5, T_all.max() * 2.0)
    ax.grid(True, alpha=0.3)
    ax.set_title('Temperature Profile')
    ax.legend(ncol=2, fontsize=8)

    # Density vs radius
    ax = axes[0, 1]
    for i, t in enumerate(selected_times):
        data = epochs[t]
        r_cm = 10**data[:, 3]  # log10(R) in cm
        rho = 10**data[:, 7] * 1e-6  # Convert from log10(rho/1e-6) to g/cm^3
        if pre_sbo_only:
            ls = '-'
            alpha = 0.8
        else:
            ls = '--' if t < 0 else '-'
            alpha = 0.7 if t < 0 else 0.8
        ax.plot(r_cm, rho, color=colors[i], 
                label=f't={t:.2f}d', linewidth=1.5, alpha=alpha, ls=ls)
        # Mark photosphere (tau ~ 2/3)
        tau_vals = data[:, 12]
        valid = tau_vals > -49
        if np.any(valid):
            i_phot = np.argmin(np.abs(tau_vals[valid] - 2.0/3.0))
            idx_phot = np.where(valid)[0][i_phot]
            ax.plot(r_cm[idx_phot], rho[idx_phot], 'o', color=colors[i], markersize=6, 
                    markeredgecolor='white', markeredgewidth=0.5, zorder=5)
    ax.set_xlabel('Radius (cm)')
    ax.set_ylabel('Density (g/cm$^3$)')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(rho_all.min() * 0.5, rho_all.max() * 2.0)
    ax.grid(True, alpha=0.3)
    ax.set_title('Density Profile')

    # Velocity vs radius
    ax = axes[1, 0]
    for i, t in enumerate(selected_times):
        data = epochs[t]
        r_cm = 10**data[:, 3]  # log10(R) in cm
        v_kms = data[:, 4] * 1e3  # V in 1000 km/s, convert to km/s
        if pre_sbo_only:
            ls = '-'
            alpha = 0.8
        else:
            ls = '--' if t < 0 else '-'
            alpha = 0.7 if t < 0 else 0.8
        ax.plot(r_cm, v_kms, color=colors[i], 
                label=f't={t:.2f}d', linewidth=1.5, alpha=alpha, ls=ls)
        # Mark photosphere (tau ~ 2/3)
        tau_vals = data[:, 12]
        valid = tau_vals > -49
        if np.any(valid):
            i_phot = np.argmin(np.abs(tau_vals[valid] - 2.0/3.0))
            idx_phot = np.where(valid)[0][i_phot]
            ax.plot(r_cm[idx_phot], v_kms[idx_phot], 'o', color=colors[i], markersize=6, 
                    markeredgecolor='white', markeredgewidth=0.5, zorder=5)
    ax.set_xlabel('Radius (cm)')
    ax.set_ylabel('Velocity (km/s)')
    ax.set_xscale('log')
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(v_all.min() * 0.5, v_all.max() * 1.2)
    ax.grid(True, alpha=0.3)
    ax.set_title('Velocity Profile')

    # Flum (luminosity flux) vs radius
    ax = axes[1, 1]
    for i, t in enumerate(selected_times):
        data = epochs[t]
        r_cm = 10**data[:, 3]  # log10(R) in cm
        flum = data[:, 11]  # Flum in 1e40 erg/s
        if pre_sbo_only:
            ls = '-'
            alpha = 0.8
        else:
            ls = '--' if t < 0 else '-'
            alpha = 0.7 if t < 0 else 0.8
        ax.plot(r_cm, flum, color=colors[i], 
                label=f't={t:.2f}d', linewidth=1.5, alpha=alpha, ls=ls)
        # Mark photosphere (tau ~ 2/3)
        tau_vals = data[:, 12]
        valid = tau_vals > -49
        if np.any(valid):
            i_phot = np.argmin(np.abs(tau_vals[valid] - 2.0/3.0))
            idx_phot = np.where(valid)[0][i_phot]
            ax.plot(r_cm[idx_phot], flum[idx_phot], 'o', color=colors[i], markersize=6, 
                    markeredgecolor='white', markeredgewidth=0.5, zorder=5)
    ax.set_xlabel('Radius (cm)')
    ax.set_ylabel('Flum (10$^{40}$ erg s$^{-1}$)')
    ax.set_xscale('log')
    ax.set_xlim(x_min, x_max)
    flum_pos = flum_all[flum_all > 0]
    flum_neg = flum_all[flum_all < 0]
    ymin = flum_neg.min() * 2.0 if len(flum_neg) > 0 else 0
    ymax = flum_pos.max() * 1.5 if len(flum_pos) > 0 else 1
    ax.set_ylim(ymin, ymax)
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
