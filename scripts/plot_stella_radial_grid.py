#!/usr/bin/env python3
"""
Plot STELLA radial profiles in an 8-panel grid similar to Kegerreis et al. (2015).

This creates a multi-panel figure showing radial profiles at a selected epoch:
- Temperature, E_int/E_rad ratio, Velocity, Luminosity (multi-group)
- Density, F_rad/cE_rad, Opacity, Optical depth

Usage:
    python plot_stella_radial_grid.py <swd_file> --time <days> [--output <path>]
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
from pathlib import Path
import sys
import argparse


def read_swd_epoch(swd_path: str, target_time: float = None):
    """Read STELLA .swd file and extract data for a specific epoch or all epochs.
    
    Returns dict with epoch times as keys, each containing zone-by-zone data.
    Columns in .swd: t(days), zone, logM, logR, V(km/s), logTg, logTr, logrho, logP, logQv, logeng, Flum(1e40), WRKX
    
    Note: Zone 1 is always the first line of an epoch (with t != 0), zones 2-500 follow (with t = 0).
    Sometimes zone 500 also has t != 0 if it's a single-zone epoch or file structure quirk.
    """
    epochs = {}
    current_t = None
    buf_data = []
    last_zone = 0
    
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
            # Only start new epoch if zone == 1 AND we have data from previous epoch
            if t != 0 and zone == 1:
                if current_t is not None and buf_data:
                    epochs[current_t] = np.array(buf_data)
                current_t = t
                buf_data = [[float(c) for c in cols[:13]]]
                last_zone = zone
            # Zone 500 with t != 0 is part of current epoch (file structure quirk)
            elif t != 0 and zone > 1 and current_t is not None and zone == 500:
                # This is zone 500 appearing as a header - add to current epoch
                buf_data.append([float(c) for c in cols[:13]])
                last_zone = zone
            # Zone data line (t = 0)
            elif t == 0 and current_t is not None:
                buf_data.append([float(c) for c in cols[:13]])
                last_zone = zone
    
    # Don't forget last epoch
    if current_t is not None and buf_data:
        epochs[current_t] = np.array(buf_data)
    
    return epochs


def find_nearest_epoch(epochs: dict, target_time: float):
    """Find the epoch closest to target_time."""
    times = np.array(sorted(epochs.keys()))
    idx = np.argmin(np.abs(times - target_time))
    return times[idx], epochs[times[idx]]


def compute_opacity(data: np.ndarray) -> np.ndarray:
    """Compute opacity from STELLA data (simplified approximation)."""
    # Extract columns
    log_rho = data[:, 8]  # log10(rho)
    log_Tg = data[:, 5]   # log10(T_gas)
    
    # Approximate electron scattering opacity: kappa = 0.2(1+X) cm^2/g for fully ionized
    # Assume solar composition: X ~ 0.7, so kappa ~ 0.34 cm^2/g
    kappa_es = 0.34  # cm^2/g
    
    rho = 10**log_rho
    T = 10**log_Tg
    
    # Opacity coefficient per unit mass (cm^2/g)
    # At high T (>10^4 K), electron scattering dominates
    # At lower T, Kramers opacity becomes important
    kappa = np.where(T > 1e4, kappa_es, kappa_es * (T/1e4)**(-3.5))
    
    # Opacity per unit length (cm^-1) = kappa * rho
    opacity = kappa * rho
    
    return opacity


def compute_optical_depth(r_cm: np.ndarray, opacity: np.ndarray) -> np.ndarray:
    """Compute optical depth from outside in (assuming spherical symmetry)."""
    # tau = integral of opacity dr from surface to radius
    # Since zones are ordered from outside (zone 1) to inside, we integrate inward
    tau = np.zeros_like(r_cm)
    
    if len(r_cm) < 2:
        return tau
    
    # Compute dr (zone widths)
    dr = np.abs(np.diff(r_cm, prepend=r_cm[0] + (r_cm[0] - r_cm[1]) if len(r_cm) > 1 else r_cm[0]))
    
    # Optical depth from outside: integrate opacity * dr from outer boundary inward
    tau[0] = opacity[0] * dr[0]  # Surface approximation
    
    for i in range(1, len(r_cm)):
        tau[i] = tau[i-1] + opacity[i] * dr[i]
    
    return tau


def plot_radial_grid(epochs: dict, target_time: float, output_file: str):
    """Create 8-panel radial profile plot similar to Kegerreis et al. (2015)."""
    
    # Find nearest epoch
    epoch_time, data = find_nearest_epoch(epochs, target_time)
    print(f"Plotting epoch at t = {epoch_time:.2f} days (target was {target_time:.2f} days)")
    
    # Extract columns
    # Format: t, zone, logM, logR, V, logTg, logTr, logrho, logP, logQv, logeng, Flum, WRKX
    r_cm = 10**data[:, 3]  # Radius in cm
    T_gas = 10**data[:, 5]  # Gas temperature in K
    T_rad = 10**data[:, 6]  # Radiation temperature in K (may be 0 in some zones)
    v_kms = data[:, 4]      # Velocity in km/s
    rho = 10**data[:, 8]    # Density in g/cm^3
    log_eng = data[:, 10]   # log10(internal energy density)
    Flum = data[:, 11] * 1e40  # Luminosity in erg/s (Flum is in 1e40)
    
    # Avoid T_rad = 0 issues
    T_rad_safe = np.where(T_rad > 0, T_rad, T_gas)
    
    # Compute derived quantities
    # E_int/E_rad ratio: internal energy / radiation energy density
    # For ideal gas: E_int = 3/2 P = 3/2 n k T = (3/2) (rho/mu mp) k T
    # For radiation: E_rad = a T^4 where a is radiation constant
    a_rad = 7.5657e-15  # erg/cm^3/K^4
    E_rad = a_rad * T_rad_safe**4
    
    # Approximate internal energy (assuming ideal gas with mu ~ 0.6 for ionized gas)
    k_B = 1.3806e-16  # erg/K
    m_p = 1.6726e-24  # g
    mu = 0.6  # Mean molecular weight for ionized gas
    E_int = (3/2) * rho * k_B * T_gas / (mu * m_p)
    
    E_ratio = E_int / E_rad
    E_ratio = np.where(E_rad > 0, E_ratio, 1e-10)  # Avoid division by zero
    
    # Opacity
    opacity = compute_opacity(data)
    
    # Optical depth
    tau = compute_optical_depth(r_cm, opacity)
    
    # F_rad / c E_rad ratio (flux factor)
    # For now, approximate using luminosity
    c = 2.998e10  # cm/s
    # L = 4 pi r^2 F, so F = L / (4 pi r^2)
    F_rad = Flum / (4 * np.pi * r_cm**2)
    F_over_cE = F_rad / (c * E_rad)
    F_over_cE = np.clip(F_over_cE, -1, 1)  # Clip to physical range
    
    # Convert radius to 10^15 cm for plotting
    r_15 = r_cm / 1e15
    
    # Create figure with 8 panels (4x2 layout) with shared x-axis
    fig, axes = plt.subplots(4, 2, figsize=(12, 16), sharex=True)
    fig.suptitle(f'STELLA Radial Profiles at t = {epoch_time:.2f} days', fontsize=14, fontweight='bold')
    
    # Set common x-axis properties
    x_min = r_15.min()
    x_max = r_15.max()
    
    # Add a common x-axis label at the bottom
    fig.text(0.5, 0.02, 'Radius (10$^{15}$ cm)', ha='center', fontsize=12)
    
    # 1. Temperature [K] (top left)
    ax = axes[0, 0]
    ax.semilogy(r_15, T_gas, 'b-', linewidth=1.5, label='T_gas')
    ax.set_ylabel('Temperature [K]', fontsize=11)
    ax.set_ylim(1e3, 1e5)
    ax.set_xlim(x_min, x_max)
    ax.grid(True, alpha=0.3)
    ax.set_title('Temperature', fontsize=12)
    
    # 2. E_int/E_rad ratio (top right)
    ax = axes[0, 1]
    ax.semilogy(r_15, E_ratio, 'b-', linewidth=1.5)
    ax.set_ylabel(r'E$_{\mathrm{int}}$/E$_{\mathrm{rad}}$', fontsize=11)
    ax.set_ylim(1e-5, 1e5)
    ax.set_xlim(x_min, x_max)
    ax.grid(True, alpha=0.3)
    ax.set_title('Internal/Radiation Energy Ratio', fontsize=12)
    
    # 3. Velocity [km/s] (2nd row left)
    ax = axes[1, 0]
    ax.plot(r_15, v_kms, 'b-', linewidth=1.5)
    ax.set_ylabel('Velocity [km s$^{-1}$]', fontsize=11)
    ax.set_ylim(0, v_kms.max() * 1.1)
    ax.set_xlim(x_min, x_max)
    ax.grid(True, alpha=0.3)
    ax.set_title('Velocity', fontsize=12)
    
    # 4. Luminosity(r) [erg/s] - Multi-group representation using color gradient (2nd row right)
    ax = axes[1, 1]
    # For now, plot bolometric luminosity with color representing radial zones
    # In a full multi-group treatment, this would show multiple frequency groups
    colors = plt.cm.rainbow(np.linspace(0, 1, len(r_15)))
    for i in range(len(r_15)-1):
        ax.plot(r_15[i:i+2], Flum[i:i+2], color=colors[i], linewidth=2)
    ax.set_ylabel(r'Luminosity(r) [erg s$^{-1}$]', fontsize=11)
    ax.set_yscale('log')
    ax.set_ylim(Flum[Flum > 0].min() * 0.1, Flum.max() * 2)
    ax.set_xlim(x_min, x_max)
    ax.grid(True, alpha=0.3)
    ax.set_title('Luminosity (Bolometric)', fontsize=12)
    
    # 5. Density [g/cm^3] (3rd row left)
    ax = axes[2, 0]
    ax.loglog(r_15, rho, 'b-', linewidth=1.5)
    ax.set_ylabel(r'Density [g cm$^{-3}$]', fontsize=11)
    ax.set_ylim(rho.min() * 0.1, rho.max() * 10)
    ax.set_xlim(x_min, x_max)
    ax.grid(True, alpha=0.3)
    ax.set_title('Density', fontsize=12)
    
    # 6. F_rad/cE_rad (3rd row right)
    ax = axes[2, 1]
    ax.plot(r_15, F_over_cE, 'b-', linewidth=1.5)
    ax.set_ylabel(r'F$_{\mathrm{rad}}$/cE$_{\mathrm{rad}}$', fontsize=11)
    ax.set_ylim(-0.1, 1.1)
    ax.set_xlim(x_min, x_max)
    ax.grid(True, alpha=0.3)
    ax.set_title('Radiation Flux Factor', fontsize=12)
    ax.axhline(y=1/3, color='r', linestyle='--', alpha=0.5, label='Diffusion limit (1/3)')
    ax.legend(fontsize=9)
    
    # 7. Opacity [cm^-1] (bottom left)
    ax = axes[3, 0]
    ax.loglog(r_15, opacity, 'b-', linewidth=1.5)
    ax.set_ylabel(r'Opacity [cm$^{-1}$]', fontsize=11)
    ax.set_ylim(opacity.min() * 0.1, opacity.max() * 10)
    ax.set_xlim(x_min, x_max)
    ax.grid(True, alpha=0.3)
    ax.set_title('Opacity', fontsize=12)
    
    # 8. Optical depth (bottom right)
    ax = axes[3, 1]
    ax.semilogy(r_15, tau, 'b-', linewidth=1.5)
    ax.set_ylabel('Optical depth', fontsize=11)
    ax.set_ylim(tau.min() * 0.1, tau.max() * 2)
    ax.set_xlim(x_min, x_max)
    ax.grid(True, alpha=0.3)
    ax.set_title('Optical Depth', fontsize=12)
    ax.axhline(y=1.0, color='r', linestyle='--', alpha=0.5, label=r'$\tau$ = 1')
    ax.axhline(y=2/3, color='g', linestyle='--', alpha=0.5, label=r'$\tau$ = 2/3 (photosphere)')
    ax.legend(fontsize=9)
    
    plt.tight_layout(rect=[0, 0, 1, 0.98])
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Radial grid plot saved to: {output_file}")
    plt.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('swd_file', help='STELLA .swd file path')
    parser.add_argument('--time', type=float, required=True, help='Target time in days')
    parser.add_argument('--output', default=None, help='Output figure path')
    args = parser.parse_args()
    
    swd_path = Path(args.swd_file)
    if not swd_path.exists():
        print(f"ERROR: {swd_path} not found")
        sys.exit(1)
    
    # Default output path
    if args.output:
        output_file = args.output
    else:
        output_file = swd_path.parent / f"{swd_path.stem}_radial_grid_t{args.time:.1f}d.png"
    
    print(f"Reading {swd_path}...")
    epochs = read_swd_epoch(str(swd_path))
    print(f"Found {len(epochs)} epochs")
    
    plot_radial_grid(epochs, args.time, str(output_file))


if __name__ == '__main__':
    main()
