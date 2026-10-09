#!/usr/bin/env python3
"""Plot STELLA input files (.hyd + .abn) to visualize density, temperature, velocity, and composition.

Usage:
    python plot_stella_input.py <hyd_file> <abn_file> [--output <png_file>]

Example:
    python plot_stella_input.py models/X_STELLA/stella/modmake/X_STELLA.hyd \
                              models/X_STELLA/stella/modmake/X_STELLA.abn \
                              --output plots/stella_input_structure.png
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator


def read_hyd(fname: str) -> dict:
    """Read STELLA .hyd file and return physical quantities."""
    with open(fname, 'r') as f:
        lines = f.readlines()
    
    # Skip comment lines (starting with #)
    data_lines = [l for l in lines if not l.strip().startswith('#')]
    
    # Parse header (first data line)
    header = data_lines[0].strip().split()
    t_start = float(header[0])
    nzone = int(header[1])
    
    # Parse zone data
    zones = []
    for line in data_lines[1:]:
        cols = line.strip().split()
        if len(cols) >= 7:
            zones.append({
                'km': int(cols[0]),
                'dMr': float(cols[1]),
                'r': float(cols[2]),
                'rho': float(cols[3]),
                'T': float(cols[4]),
                'u': float(cols[5]),
                'Mr': float(cols[6])
            })
    
    return {
        't_start': t_start,
        'nzone': nzone,
        'km': np.array([z['km'] for z in zones]),
        'r': np.array([z['r'] for z in zones]),
        'rho': np.array([z['rho'] for z in zones]),
        'T': np.array([z['T'] for z in zones]),
        'u': np.array([z['u'] for z in zones]),
        'Mr': np.array([z['Mr'] for z in zones]),
        'dMr': np.array([z['dMr'] for z in zones])
    }


def read_abn(fname: str) -> dict:
    """Read STELLA .abn file and return mass fractions."""
    with open(fname, 'r') as f:
        lines = f.readlines()
    
    # Skip comment lines
    data_lines = [l for l in lines if not l.strip().startswith('#')]
    
    # Columns: km, dum, dum, dum, H, He, C, N, O, Ne, Na, Mg, Al, Si, S, Ar, Ca, FePeak, Ni58, Ni56
    elements = ['H', 'He', 'C', 'N', 'O', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'S', 'Ar', 'Ca', 'FePeak', 'Ni58', 'Ni56']
    mass_fracs = {el: [] for el in elements}
    
    for line in data_lines:
        cols = line.strip().split()
        if len(cols) >= 20:
            mass_fracs['H'].append(float(cols[4]))
            mass_fracs['He'].append(float(cols[5]))
            mass_fracs['C'].append(float(cols[6]))
            mass_fracs['N'].append(float(cols[7]))
            mass_fracs['O'].append(float(cols[8]))
            mass_fracs['Ne'].append(float(cols[9]))
            mass_fracs['Na'].append(float(cols[10]))
            mass_fracs['Mg'].append(float(cols[11]))
            mass_fracs['Al'].append(float(cols[12]))
            mass_fracs['Si'].append(float(cols[13]))
            mass_fracs['S'].append(float(cols[14]))
            mass_fracs['Ar'].append(float(cols[15]))
            mass_fracs['Ca'].append(float(cols[16]))
            mass_fracs['FePeak'].append(float(cols[17]))
            mass_fracs['Ni58'].append(float(cols[18]))
            mass_fracs['Ni56'].append(float(cols[19]))
    
    return {el: np.array(vals) for el, vals in mass_fracs.items()}


def plot_stella_input(hyd_file: str, abn_file: str, output_file: str = None):
    """Create comprehensive plots of STELLA input."""
    
    # Read data
    print(f"Reading {hyd_file}...")
    hyd = read_hyd(hyd_file)
    
    print(f"Reading {abn_file}...")
    abn = read_abn(abn_file)
    
    print(f"Model: {hyd['nzone']} zones")
    print(f"Radius range: {hyd['r'][0]:.2e} to {hyd['r'][-1]:.2e} cm")
    print(f"Density range: {hyd['rho'].min():.2e} to {hyd['rho'].max():.2e} g/cc")
    print(f"Velocity range: {hyd['u'].min()/1e5:.0f} to {hyd['u'].max()/1e5:.0f} km/s")
    
    # Create figure with 4 subplots
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Plot 1: Density vs Radius
    ax1 = axes[0, 0]
    ax1.loglog(hyd['r'], hyd['rho'], 'b-', linewidth=1.5)
    ax1.set_xlabel('Radius (cm)', fontsize=11)
    ax1.set_ylabel('Density (g/cm³)', fontsize=11)
    ax1.set_title('Density Profile', fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.axhline(1e-14, color='r', linestyle='--', alpha=0.5, label='ρ = 1e-14 g/cc')
    ax1.legend()
    
    # Add text annotations for regions
    ax1.text(hyd['r'][0]*2, hyd['rho'][0]*0.5, 'Inner\n(Ejecta A)', 
             fontsize=9, ha='center', alpha=0.7)
    ax1.text(np.sqrt(hyd['r'][0]*hyd['r'][-1]), 1e-12, 'CSM Wind', 
             fontsize=9, ha='center', alpha=0.7)
    ax1.text(hyd['r'][-1]/3, hyd['rho'][-1]*10, 'Outer\n(Ejecta B)', 
             fontsize=9, ha='center', alpha=0.7)
    
    # Plot 2: Temperature vs Radius
    ax2 = axes[0, 1]
    ax2.loglog(hyd['r'], hyd['T'], 'r-', linewidth=1.5)
    ax2.set_xlabel('Radius (cm)', fontsize=11)
    ax2.set_ylabel('Temperature (K)', fontsize=11)
    ax2.set_title('Temperature Profile', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Velocity vs Radius
    ax3 = axes[1, 0]
    ax3.semilogx(hyd['r'], hyd['u']/1e5, 'g-', linewidth=1.5)
    ax3.set_xlabel('Radius (cm)', fontsize=11)
    ax3.set_ylabel('Velocity (km/s)', fontsize=11)
    ax3.set_title('Velocity Profile', fontsize=12, fontweight='bold')
    ax3.grid(True, alpha=0.3)
    ax3.axhline(1000, color='k', linestyle='--', alpha=0.3, label='1,000 km/s')
    ax3.axhline(10000, color='k', linestyle='--', alpha=0.5, label='10,000 km/s')
    ax3.legend()
    
    # Plot 4: Stacked Mass Fractions (Composition)
    ax4 = axes[1, 1]
    
    # Select major elements for plotting (include Ni56 if present)
    major_elements = ['H', 'He', 'C', 'O', 'Si', 'FePeak']
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']
    
    # Add Ni56 if it has significant mass fraction (> 1e-4)
    if 'Ni56' in abn and np.max(abn['Ni56']) > 1e-4:
        major_elements.append('Ni56')
        colors.append('#e377c2')  # Magenta for Ni56
        print(f"Adding Ni56 to composition plot (max fraction: {np.max(abn['Ni56']):.4f})")
    
    # Create stacked area plot
    bottom = np.zeros(len(hyd['r']))
    for el, color in zip(major_elements, colors):
        if el in abn and np.any(abn[el] > 0):
            ax4.fill_between(hyd['r'], bottom, bottom + abn[el], 
                           label=el, color=color, alpha=0.8)
            bottom += abn[el]
    
    ax4.set_xscale('log')
    ax4.set_xlabel('Radius (cm)', fontsize=11)
    ax4.set_ylabel('Mass Fraction', fontsize=11)
    ax4.set_title('Composition Profile (Mass Fractions)', fontsize=12, fontweight='bold')
    ax4.set_ylim(0, 1)
    ax4.legend(loc='upper right', ncol=2, fontsize=9)
    ax4.grid(True, alpha=0.3)

    # Secondary top x-axis: enclosed mass vs radius
    r_arr = hyd['r']
    m_arr = hyd['Mr']
    m_ticks = MaxNLocator(nbins=10, prune='both').tick_values(
        m_arr.min(), m_arr.max()
    )
    r_ticks = np.interp(m_ticks, m_arr, r_arr)
    keep = np.concatenate(([True], np.diff(np.log10(r_ticks)) > 0.12))
    m_ticks = m_ticks[keep]
    for a in (ax1, ax2):
        sec = a.secondary_xaxis(
            'top',
            functions=(
                lambda rr: np.interp(rr, r_arr, m_arr),
                lambda mm: np.interp(mm, m_arr, r_arr),
            ),
        )
        sec.set_xticks(m_ticks)
        sec.set_xticklabels([f"{t:g}" for t in m_ticks])
        sec.set_xlabel('Enclosed Mass (M$_\\odot$)', fontsize=10)

    # Add overall title
    fig.suptitle(f'STELLA Input Structure: {Path(hyd_file).stem}',
                 fontsize=14, fontweight='bold', y=0.98)

    plt.tight_layout()
    
    # Save or show
    if output_file:
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"Saved plot to: {output_file}")
    else:
        plt.show()


def main():
    parser = argparse.ArgumentParser(description='Plot STELLA input files')
    parser.add_argument('hyd_file', help='Path to .hyd file')
    parser.add_argument('abn_file', help='Path to .abn file')
    parser.add_argument('--output', '-o', help='Output PNG file (optional)')
    args = parser.parse_args()
    
    if not Path(args.hyd_file).exists():
        print(f"Error: {args.hyd_file} not found")
        sys.exit(1)
    
    if not Path(args.abn_file).exists():
        print(f"Error: {args.abn_file} not found")
        sys.exit(1)
    
    plot_stella_input(args.hyd_file, args.abn_file, args.output)


if __name__ == '__main__':
    main()
