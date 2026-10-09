#!/usr/bin/env python3
"""
Plot STELLA radial profiles as a time series (multi-epoch overlay or animation).

This creates either:
1. Multi-epoch overlay: All epochs plotted on 8-panel grid with color-coded time
2. Animation: Frame-by-frame evolution of radial profiles

Usage:
    python plot_stella_radial_timeseries.py <swd_file> [OPTIONS]
    
Options:
    --mode {overlay,animation}  Plot mode (default: overlay)
    --n-epochs N                Number of epochs for overlay (default: 10)
    --fps FPS                   Frames per second for animation (default: 5)
    --output PATH               Output file path
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
from matplotlib.animation import FuncAnimation, PillowWriter
from pathlib import Path
import sys
import argparse


def read_swd_all_epochs(swd_path: str) -> dict:
    """Read all epochs from STELLA .swd file.
    
    Returns dict with epoch times as keys, each containing zone-by-zone data.
    Columns: t(days), zone, logM, logR, V(1e8 cm/s), logTg, logTr, logrho, logP, logQv, logeng, Flum(1e40), WRKX
    """
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
            # Zone 500 with t != 0 is part of current epoch
            elif t != 0 and zone == 500 and current_t is not None:
                buf_data.append([float(c) for c in cols[:13]])
            # Zone data line (t = 0)
            elif t == 0 and current_t is not None:
                buf_data.append([float(c) for c in cols[:13]])
    
    if current_t is not None and buf_data:
        epochs[current_t] = np.array(buf_data)
    
    return epochs


def compute_derived_quantities(data: np.ndarray):
    """Compute derived physical quantities from STELLA data."""
    # Extract columns with proper unit conversions
    # STELLA units: V in 1e8 cm/s, rho in log10(rho/1e-6 g/cc)
    r_cm = 10**data[:, 3]  # Radius in cm
    T_gas = 10**data[:, 5]  # Gas temperature in K
    T_rad = 10**data[:, 6]  # Radiation temperature in K
    v_kms = data[:, 4] * 1000.0  # Convert from 1e8 cm/s to km/s
    rho = (10**data[:, 8]) * 1e-6  # Convert from log10(rho/1e-6) to g/cm^3
    log_eng = data[:, 10]   # log10(internal energy density)
    Flum = data[:, 11] * 1e40  # Luminosity in erg/s
    
    # Avoid T_rad = 0 issues
    T_rad_safe = np.where(T_rad > 0, T_rad, T_gas)
    
    # E_int/E_rad ratio
    a_rad = 7.5657e-15  # erg/cm^3/K^4
    E_rad = a_rad * T_rad_safe**4
    k_B = 1.3806e-16  # erg/K
    m_p = 1.6726e-24  # g
    mu = 0.6
    E_int = (3/2) * rho * k_B * T_gas / (mu * m_p)
    E_ratio = np.where(E_rad > 0, E_int / E_rad, 1e-10)
    
    # Opacity (electron scattering approximation)
    kappa_es = 0.34  # cm^2/g
    T = 10**data[:, 5]
    kappa = np.where(T > 1e4, kappa_es, kappa_es * (T/1e4)**(-3.5))
    opacity = kappa * rho
    
    # Optical depth
    tau = np.zeros_like(r_cm)
    if len(r_cm) > 1:
        dr = np.abs(np.diff(r_cm, prepend=r_cm[0] + (r_cm[0] - r_cm[1])))
        tau[0] = opacity[0] * dr[0]
        for i in range(1, len(r_cm)):
            tau[i] = tau[i-1] + opacity[i] * dr[i]
    
    # F_rad / c E_rad
    c = 2.998e10
    F_rad = Flum / (4 * np.pi * r_cm**2)
    F_over_cE = F_rad / (c * E_rad)
    F_over_cE = np.clip(F_over_cE, -1, 1)
    
    return {
        'r_cm': r_cm,
        'r_15': r_cm / 1e15,
        'T_gas': T_gas,
        'T_rad': T_rad_safe,
        'v_kms': v_kms,
        'rho': rho,
        'E_ratio': E_ratio,
        'opacity': opacity,
        'tau': tau,
        'F_over_cE': F_over_cE,
        'Flum': Flum
    }


def plot_overlay(epochs: dict, n_epochs: int, output_file: str):
    """Create multi-epoch overlay plot with color-coded time."""
    times = sorted(epochs.keys())
    
    # Select epochs at 3-day intervals
    t_min = times[0]
    t_max = times[-1]
    time_grid = np.arange(t_min, t_max + 3, 3)  # 3-day intervals
    
    # Find closest epochs to the 3-day grid points
    selected_times = []
    for t_target in time_grid:
        # Find closest epoch
        idx = np.abs(np.array(times) - t_target).argmin()
        if times[idx] not in selected_times:
            selected_times.append(times[idx])
    
    # Limit to n_epochs if specified
    if len(selected_times) > n_epochs:
        selected_times = selected_times[:n_epochs]
    
    # Filter out negative timesteps
    selected_times = [t for t in selected_times if t >= 0]
    
    if not selected_times:
        print("ERROR: No non-negative timesteps found")
        return
    
    # Create colormap
    cmap = plt.cm.viridis
    norm = plt.Normalize(vmin=min(selected_times), vmax=max(selected_times))
    
    fig, axes = plt.subplots(4, 2, figsize=(14, 16), sharex=True)
    fig.suptitle(f'STELLA Radial Profiles Time Series ({len(selected_times)} epochs)', 
                 fontsize=14, fontweight='bold')
    
    # Common x-axis label
    fig.text(0.5, 0.02, 'Radius (10$^{15}$ cm)', ha='center', fontsize=12)
    
    # Calculate x-axis range and all y-ranges from all selected (post-SBO)
    # epochs together, so no single epoch dominates the scale.  Pre-SBO
    # epochs are already excluded by the t >= 0 filter above.
    x_min_all = 1e10
    x_max_all = 0
    rho_min_all, rho_max_all = np.inf, 0
    T_min_all, T_max_all = np.inf, 0
    v_max_all = 0
    flum_min_pos, flum_max = np.inf, 0
    opacity_min_all, opacity_max_all = np.inf, 0
    tau_min_all, tau_max_all = np.inf, 0
    for t in selected_times:
        data = epochs[t]
        q = compute_derived_quantities(data)
        x_min_all = min(x_min_all, q['r_15'].min())
        x_max_all = max(x_max_all, q['r_15'].max())
        rho_pos = q['rho'][q['rho'] > 0]
        if len(rho_pos) > 0:
            rho_min_all = min(rho_min_all, rho_pos.min())
            rho_max_all = max(rho_max_all, rho_pos.max())
        T_pos = q['T_gas'][q['T_gas'] > 0]
        if len(T_pos) > 0:
            T_min_all = min(T_min_all, T_pos.min())
            T_max_all = max(T_max_all, T_pos.max())
        v_max_all = max(v_max_all, q['v_kms'].max())
        flum_pos = q['Flum'][q['Flum'] > 0]
        if len(flum_pos) > 0:
            flum_min_pos = min(flum_min_pos, flum_pos.min())
            flum_max = max(flum_max, flum_pos.max())
        op_pos = q['opacity'][q['opacity'] > 0]
        if len(op_pos) > 0:
            opacity_min_all = min(opacity_min_all, op_pos.min())
            opacity_max_all = max(opacity_max_all, op_pos.max())
        tau_pos = q['tau'][q['tau'] > 0]
        if len(tau_pos) > 0:
            tau_min_all = min(tau_min_all, tau_pos.min())
            tau_max_all = max(tau_max_all, tau_pos.max())

    x_min = x_min_all * 0.9
    x_max = x_max_all * 1.1
    rho_ylim = (rho_min_all * 0.1, rho_max_all * 10) if rho_max_all > 0 else None
    T_ylim = (T_min_all * 0.5, T_max_all * 2.0) if T_max_all > 0 else (1e3, 1e5)
    v_ylim = (0, v_max_all * 1.1) if v_max_all > 0 else (0, 1)
    flum_ylim = (flum_min_pos * 0.1, flum_max * 2) if flum_max > 0 else (1e40, 1e44)
    opacity_ylim = (opacity_min_all * 0.1, opacity_max_all * 10) if opacity_max_all > 0 else None
    tau_ylim = (tau_min_all * 0.1, tau_max_all * 2) if tau_max_all > 0 else None
    
    # Store legend handles and labels
    legend_handles = []
    legend_labels = []
    
    for t in selected_times:
        data = epochs[t]
        q = compute_derived_quantities(data)
        color = cmap(norm(t))
        r_15 = q['r_15']  # Use 10^15 cm units
        
        # 1. Temperature
        ax = axes[0, 0]
        ax.semilogy(r_15, q['T_gas'], color=color, linewidth=1.5, alpha=0.8)
        ax.set_ylabel('Temperature [K]', fontsize=10)
        ax.set_ylim(T_ylim)
        ax.set_xlim(x_min, x_max)
        ax.grid(True, alpha=0.3)
        ax.set_title('Temperature', fontsize=11)
        
        # 2. E_int/E_rad
        ax = axes[0, 1]
        ax.semilogy(r_15, q['E_ratio'], color=color, linewidth=1.5, alpha=0.8)
        ax.set_ylabel(r'E$_{\mathrm{int}}$/E$_{\mathrm{rad}}$', fontsize=10)
        ax.set_ylim(1e-5, 1e5)
        ax.set_xlim(x_min, x_max)
        ax.grid(True, alpha=0.3)
        ax.set_title('Internal/Radiation Energy Ratio', fontsize=11)
        
        # 3. Velocity
        ax = axes[1, 0]
        ax.plot(r_15, q['v_kms'], color=color, linewidth=1.5, alpha=0.8)
        ax.set_ylabel('Velocity [km s$^{-1}$]', fontsize=10)
        ax.set_ylim(v_ylim)
        ax.set_xlim(x_min, x_max)
        ax.grid(True, alpha=0.3)
        ax.set_title('Velocity', fontsize=11)
        
        # 4. Luminosity
        ax = axes[1, 1]
        ax.semilogy(r_15, q['Flum'], color=color, linewidth=1.5, alpha=0.8)
        ax.set_ylabel(r'Luminosity [erg s$^{-1}$]', fontsize=10)
        ax.set_ylim(flum_ylim)
        ax.set_xlim(x_min, x_max)
        ax.grid(True, alpha=0.3)
        ax.set_title('Luminosity', fontsize=11)
        
        # 5. Density
        ax = axes[2, 0]
        ax.loglog(r_15, q['rho'], color=color, linewidth=1.5, alpha=0.8)
        ax.set_ylabel(r'Density [g cm$^{-3}$]', fontsize=10)
        if rho_ylim is not None:
            ax.set_ylim(rho_ylim)
        ax.set_xlim(x_min, x_max)
        ax.grid(True, alpha=0.3)
        ax.set_title('Density', fontsize=11)
        
        # 6. F_rad/cE_rad
        ax = axes[2, 1]
        ax.plot(r_15, q['F_over_cE'], color=color, linewidth=1.5, alpha=0.8)
        ax.set_ylabel(r'F$_{\mathrm{rad}}$/cE$_{\mathrm{rad}}$', fontsize=10)
        ax.set_ylim(-0.1, 1.1)
        ax.set_xlim(x_min, x_max)
        ax.grid(True, alpha=0.3)
        ax.set_title('Radiation Flux Factor', fontsize=11)
        
        # 7. Opacity
        ax = axes[3, 0]
        ax.loglog(r_15, q['opacity'], color=color, linewidth=1.5, alpha=0.8)
        ax.set_ylabel(r'Opacity [cm$^{-1}$]', fontsize=10)
        if opacity_ylim is not None:
            ax.set_ylim(opacity_ylim)
        ax.set_xlim(x_min, x_max)
        ax.grid(True, alpha=0.3)
        ax.set_title('Opacity', fontsize=11)
        
        # 8. Optical depth
        ax = axes[3, 1]
        ax.semilogy(r_15, q['tau'], color=color, linewidth=1.5, alpha=0.8)
        ax.set_ylabel('Optical depth', fontsize=10)
        if tau_ylim is not None:
            ax.set_ylim(tau_ylim)
        ax.set_xlim(x_min, x_max)
        ax.grid(True, alpha=0.3)
        ax.set_title('Optical Depth', fontsize=11)
        
        # Add to legend
        legend_handles.append(plt.Line2D([0], [0], color=color, linewidth=1.5))
        legend_labels.append(f'{t:.1f} d')
    
    # Add reference lines to last panels
    axes[2, 1].axhline(y=1/3, color='k', linestyle='--', alpha=0.3, linewidth=0.8)
    axes[3, 1].axhline(y=1.0, color='k', linestyle='--', alpha=0.3, linewidth=0.8)
    axes[3, 1].axhline(y=2/3, color='k', linestyle=':', alpha=0.3, linewidth=0.8)
    
    # Add legend to the first panel
    axes[0, 0].legend(legend_handles, legend_labels, fontsize=8, ncol=2, loc='upper right')
    
    plt.tight_layout(rect=[0, 0.03, 1, 0.96])
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Time series overlay saved to: {output_file}")
    plt.close()


def plot_animation(epochs: dict, fps: int, output_file: str):
    """Create animated GIF showing radial profile evolution."""
    times = sorted(epochs.keys())
    
    # Filter out negative timesteps
    times = [t for t in times if t >= 0]
    
    if not times:
        print("ERROR: No non-negative timesteps found")
        return
    
    # Use every Nth epoch for smoother animation
    # Target 30 seconds at 5 FPS = 150 frames
    target_frames = 150
    frame_skip = max(1, len(times) // target_frames)
    selected_times = times[::frame_skip]
    
    fig, axes = plt.subplots(4, 2, figsize=(14, 16), sharex=True)
    fig.suptitle('STELLA Radial Profile Evolution', fontsize=14, fontweight='bold')
    fig.text(0.5, 0.02, 'Radius (10$^{15}$ cm)', ha='center', fontsize=12)
    
    # Time in title - will be updated each frame
    
    x_min = 1e-10
    x_max = 20
    
    # Initialize empty lines for each panel
    lines = []
    for ax_row in axes:
        for ax in ax_row:
            line, = ax.plot([], [], 'b-', linewidth=2)
            lines.append(line)
    
    # Set up axes properties (static with sensible defaults)
    axes[0, 0].set_ylabel('Temperature [K]', fontsize=10)
    axes[0, 0].set_ylim(1e3, 1e5)
    axes[0, 0].set_yscale('log')
    axes[0, 0].set_xlim(x_min, x_max)
    axes[0, 0].grid(True, alpha=0.3)
    axes[0, 0].set_title('Temperature', fontsize=11)
    
    axes[0, 1].set_ylabel(r'E$_{\mathrm{int}}$/E$_{\mathrm{rad}}$', fontsize=10)
    axes[0, 1].set_ylim(1e-10, 1e10)
    axes[0, 1].set_yscale('log')
    axes[0, 1].set_xlim(x_min, x_max)
    axes[0, 1].grid(True, alpha=0.3)
    axes[0, 1].set_title('Internal/Radiation Energy Ratio', fontsize=11)
    
    axes[1, 0].set_ylabel('Velocity [km s$^{-1}$]', fontsize=10)
    axes[1, 0].set_ylim(0, 20000)
    axes[1, 0].set_xlim(x_min, x_max)
    axes[1, 0].grid(True, alpha=0.3)
    axes[1, 0].set_title('Velocity', fontsize=11)
    
    axes[1, 1].set_ylabel(r'Luminosity [erg s$^{-1}$]', fontsize=10)
    axes[1, 1].set_ylim(1e35, 1e45)
    axes[1, 1].set_yscale('log')
    axes[1, 1].set_xlim(x_min, x_max)
    axes[1, 1].grid(True, alpha=0.3)
    axes[1, 1].set_title('Luminosity', fontsize=11)
    
    axes[2, 0].set_ylabel(r'Density [g cm$^{-3}$]', fontsize=10)
    axes[2, 0].set_ylim(1e-20, 1e-8)
    axes[2, 0].set_yscale('log')
    axes[2, 0].set_xlim(x_min, x_max)
    axes[2, 0].grid(True, alpha=0.3)
    axes[2, 0].set_title('Density', fontsize=11)
    
    axes[2, 1].set_ylabel(r'F$_{\mathrm{rad}}$/cE$_{\mathrm{rad}}$', fontsize=10)
    axes[2, 1].set_ylim(-0.1, 1.1)
    axes[2, 1].set_xlim(x_min, x_max)
    axes[2, 1].grid(True, alpha=0.3)
    axes[2, 1].set_title('Radiation Flux Factor', fontsize=11)
    axes[2, 1].axhline(y=1/3, color='k', linestyle='--', alpha=0.3)
    
    axes[3, 0].set_ylabel(r'Opacity [cm$^{-1}$]', fontsize=10)
    axes[3, 0].set_ylim(1e-25, 1e-10)
    axes[3, 0].set_yscale('log')
    axes[3, 0].set_xlim(x_min, x_max)
    axes[3, 0].grid(True, alpha=0.3)
    axes[3, 0].set_title('Opacity', fontsize=11)
    
    axes[3, 1].set_ylabel('Optical depth', fontsize=10)
    axes[3, 1].set_ylim(1e-10, 1e10)
    axes[3, 1].set_yscale('log')
    axes[3, 1].set_xlim(x_min, x_max)
    axes[3, 1].grid(True, alpha=0.3)
    axes[3, 1].set_title('Optical Depth', fontsize=11)
    axes[3, 1].axhline(y=1.0, color='k', linestyle='--', alpha=0.3)
    axes[3, 1].axhline(y=2/3, color='k', linestyle=':', alpha=0.3)
    
    def init():
        for line in lines:
            line.set_data([], [])
        fig.suptitle('STELLA Radial Profile Evolution - t = 0.00 days', fontsize=14, fontweight='bold')
        return lines
    
    def update(frame):
        t = selected_times[frame]
        data = epochs[t]
        q = compute_derived_quantities(data)
        
        # Filter positive values for log scales
        r_T = q['r_15'][q['T_gas'] > 0]
        T_pos = q['T_gas'][q['T_gas'] > 0]
        
        r_E = q['r_15'][q['E_ratio'] > 0]
        E_pos = q['E_ratio'][q['E_ratio'] > 0]
        
        r_Flum = q['r_15'][q['Flum'] > 0]
        Flum_pos = q['Flum'][q['Flum'] > 0]
        
        r_rho = q['r_15'][q['rho'] > 0]
        rho_pos = q['rho'][q['rho'] > 0]
        
        r_op = q['r_15'][q['opacity'] > 0]
        op_pos = q['opacity'][q['opacity'] > 0]
        
        r_tau = q['r_15'][q['tau'] > 0]
        tau_pos = q['tau'][q['tau'] > 0]
        
        lines[0].set_data(r_T, T_pos)
        lines[1].set_data(r_E, E_pos)
        lines[2].set_data(q['r_15'], q['v_kms'])
        lines[3].set_data(r_Flum, Flum_pos)
        lines[4].set_data(r_rho, rho_pos)
        lines[5].set_data(q['r_15'], q['F_over_cE'])
        lines[6].set_data(r_op, op_pos)
        lines[7].set_data(r_tau, tau_pos)
        
        # Update y-limits dynamically
        if len(q['v_kms']) > 0:
            v_max = max(q['v_kms'])
            axes[1, 0].set_ylim(0, max(1, v_max * 1.1))
        if len(Flum_pos) > 0:
            axes[1, 1].set_ylim(Flum_pos.min() * 0.1, Flum_pos.max() * 2)
        if len(rho_pos) > 0:
            axes[2, 0].set_ylim(rho_pos.min() * 0.1, rho_pos.max() * 10)
        if len(op_pos) > 0:
            axes[3, 0].set_ylim(op_pos.min() * 0.1, op_pos.max() * 10)
        if len(tau_pos) > 0:
            axes[3, 1].set_ylim(tau_pos.min() * 0.1, tau_pos.max() * 2)
        
        fig.suptitle(f'STELLA Radial Profile Evolution - t = {t:.2f} days', fontsize=14, fontweight='bold')
        
        return lines
    
    anim = FuncAnimation(fig, update, init_func=init, frames=len(selected_times),
                         interval=1000//fps, blit=False, repeat=True)
    
    # Save with looping
    writer = PillowWriter(fps=fps)
    anim.save(output_file, writer=writer)
    
    # Make the GIF loop infinitely by modifying the file
    import PIL.Image as Image
    from PIL import Image
    
    # Open and re-save with loop=0
    with Image.open(output_file) as img:
        frames = []
        try:
            while True:
                frame = img.copy()
                frames.append(frame)
                img.seek(img.tell() + 1)
        except EOFError:
            pass
        
        if frames:
            frames[0].save(output_file, save_all=True, append_images=frames[1:],
                          duration=1000//fps, loop=0)
    print(f"Animation saved to: {output_file}")
    plt.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('swd_file', help='STELLA .swd file path')
    parser.add_argument('--mode', choices=['overlay', 'animation'], default='overlay',
                        help='Plot mode (default: overlay)')
    parser.add_argument('--n-epochs', type=int, default=10,
                        help='Number of epochs for overlay mode (default: 10)')
    parser.add_argument('--fps', type=int, default=5,
                        help='Frames per second for animation (default: 5)')
    parser.add_argument('--output', default=None, help='Output file path')
    args = parser.parse_args()
    
    swd_path = Path(args.swd_file)
    if not swd_path.exists():
        print(f"ERROR: {swd_path} not found")
        sys.exit(1)
    
    # Default output path
    if args.output:
        output_file = args.output
    else:
        if args.mode == 'overlay':
            output_file = swd_path.parent / f"{swd_path.stem}_radial_timeseries.png"
        else:
            output_file = swd_path.parent / f"{swd_path.stem}_radial_animation.gif"
    
    print(f"Reading {swd_path}...")
    epochs = read_swd_all_epochs(str(swd_path))
    print(f"Found {len(epochs)} epochs")
    
    if args.mode == 'overlay':
        plot_overlay(epochs, args.n_epochs, str(output_file))
    else:
        plot_animation(epochs, args.fps, str(output_file))


if __name__ == '__main__':
    main()
