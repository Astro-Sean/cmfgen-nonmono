#!/usr/bin/env python3
"""Plot input structure (R, V, T, rho) and mass fractions for SEDONA .mod and SuperLite input.str.

Both codes read from the same SN_HYDRO_DATA (via CMFGEN conversion), so this script
verifies that the inputs are consistent and produces diagnostic plots for each run.

Usage:
  python plot_input_structure.py --sedona-mod <path.mod> --superlite-str <path/input.str> \
    --sn-hydro <path/SN_HYDRO_DATA> -o <output_dir>
"""

import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path


def read_sn_hydro_data(filepath):
    """Read CMFGEN SN_HYDRO_DATA file."""
    with open(filepath, 'r') as f:
        lines = f.readlines()
    header = {}
    for line in lines:
        if line.startswith('Number of data points'):
            header['nd'] = int(line.split(':')[1].strip())
        elif line.startswith('Time(days)'):
            header['t_days'] = float(line.split(':')[1].strip())
        elif line.startswith('Energy Lum'):
            header['lum_lsun'] = float(line.split(':')[1].strip())
    nd = header['nd']
    def read_section(keyword):
        for i, line in enumerate(lines):
            if line.strip().startswith(keyword):
                data = []
                j = i + 1
                while len(data) < nd:
                    parts = lines[j].split()
                    for p in parts:
                        data.append(float(p))
                    j += 1
                return np.array(data[:nd])
        raise ValueError(f"Section '{keyword}' not found")
    r = read_section('Radius grid') * 1e10  # cm
    v = read_section('Velocity') * 1e5      # cm/s
    T = read_section('Temperature') * 1e4   # K
    rho = read_section('Density')           # g/cm^3
    mf_labels = ['HYD', 'HE', 'CARB', 'OXY', 'SIL', 'CAL', 'NICK', 'COB', 'IRON']
    mass_fracs = {}
    for label in mf_labels:
        mass_fracs[label] = read_section(f'{label} mass fraction')
    return header, r, v, T, rho, mass_fracs


def read_sedona_mod(filepath):
    """Read SEDONA .mod file (1D_sphere standard format)."""
    with open(filepath, 'r') as f:
        lines = f.readlines()
    # Line 1: nzones rmin t_exp nsources
    parts = lines[1].split()
    nzones = int(parts[0])
    # Line 2: species Z list
    # Data lines: r v rho T X(H) X(He) ...
    data = []
    for line in lines[3:]:
        vals = line.split()
        if len(vals) >= 4:
            data.append([float(x) for x in vals])
    data = np.array(data)
    r = data[:, 0]
    v = data[:, 1]
    rho = data[:, 2]
    T = data[:, 3]
    # Mass fractions start at column 4
    species_line = lines[2].split()
    z_list = [int(z) for z in species_line]
    elem_names = {1: 'H', 2: 'He', 6: 'C', 8: 'O', 14: 'Si', 20: 'Ca', 26: 'Fe', 27: 'Co', 28: 'Ni'}
    mf = {}
    for i, z in enumerate(z_list):
        name = elem_names.get(z, f'Z{z}')
        mf[name] = data[:, 4 + i] if 4 + i < data.shape[1] else np.zeros(nzones)
    return r, v, rho, T, mf


def read_superlite_str(filepath):
    """Read SuperLite input.str file."""
    with open(filepath, 'r') as f:
        lines = f.readlines()
    # Line 2: # nx ny nz ncol nabund
    parts = lines[1].split()
    nx = int(parts[1])
    nabund = int(parts[5])
    # Line 3: column labels
    labels = lines[2].split()
    # Data lines
    data = []
    for line in lines[3:]:
        vals = line.split()
        if len(vals) >= 7:
            data.append([float(x) for x in vals])
    data = np.array(data)
    # Columns: x_left x_right vx_left vx_right mass rad_temp temp h he c o si ca ni co fe
    r_left = data[:, 0]
    r_right = data[:, 1]
    r_center = np.sqrt(r_left * r_right)
    v_left = data[:, 2]
    v_right = data[:, 3]
    v_center = 0.5 * (v_left + v_right)
    mass = data[:, 4]
    T = data[:, 6]
    # Compute density from mass and volume
    vol = (4.0 / 3.0) * np.pi * (r_right**3 - r_left**3)
    rho = mass / vol
    # Mass fractions start at column 7
    mf_labels_str = labels[7:7 + nabund]
    mf = {}
    elem_map = {'h': 'H', 'he': 'He', 'c': 'C', 'o': 'O', 'si': 'Si', 'ca': 'Ca',
                'ni': 'Ni', 'co': 'Co', 'fe': 'Fe'}
    for i, lbl in enumerate(mf_labels_str):
        name = elem_map.get(lbl.lower(), lbl)
        mf[name] = data[:, 7 + i]
    return r_center, v_center, rho, T, mf


def plot_structure(r, v, T, rho, label, color, ax_r, ax_v, ax_T, ax_rho, t_days=None):
    """Plot structure on the given axes."""
    # Convert to useful units
    r_solar = r / 6.96e10  # solar radii
    v_kms = v / 1e5
    # Radius
    ax_r.plot(r_solar, label=label, color=color, linewidth=1.5)
    # Velocity
    ax_v.plot(r_solar, v_kms, label=label, color=color, linewidth=1.5)
    # Temperature
    ax_T.plot(r_solar, T, label=label, color=color, linewidth=1.5)
    # Density
    ax_rho.plot(r_solar, rho, label=label, color=color, linewidth=1.5)


def plot_mass_fractions(r, mf, label_prefix, color_map, ax):
    """Plot mass fractions vs radius."""
    for name, vals in mf.items():
        if vals.max() > 1e-10:  # Only plot non-zero elements
            ax.plot(r / 6.96e10, vals, label=f'{label_prefix}{name}',
                    color=color_map.get(name, 'gray'), linewidth=1.2, linestyle='-')


def main():
    parser = argparse.ArgumentParser(description='Plot input structure for SEDONA and SuperLite')
    parser.add_argument('--sedona-mod', type=str, default=None, help='SEDONA .mod file')
    parser.add_argument('--superlite-str', type=str, default=None, help='SuperLite input.str file')
    parser.add_argument('--sn-hydro', type=str, default=None, help='CMFGEN SN_HYDRO_DATA file')
    parser.add_argument('-o', '--output-dir', type=str, required=True, help='Output directory for plots')
    parser.add_argument('--model-name', type=str, default='', help='Model name for titles')
    args = parser.parse_args()

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Color scheme
    colors = {'SN_HYDRO_DATA': '#2c3e50', 'SEDONA': '#e74c3c', 'SuperLite': '#2980b9'}
    elem_colors = {
        'H': '#e74c3c', 'He': '#3498db', 'C': '#2ecc71', 'O': '#9b59b6',
        'Si': '#f39c12', 'Ca': '#1abc9c', 'Fe': '#e67e22', 'Ni': '#95a5a6', 'Co': '#7f8c8d'
    }

    # --- Structure comparison plot ---
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    ax_r, ax_v = axes[0]
    ax_T, ax_rho = axes[1]

    datasets = []

    # Read SN_HYDRO_DATA (CMFGEN source)
    if args.sn_hydro and Path(args.sn_hydro).exists():
        header, r_h, v_h, T_h, rho_h, mf_h = read_sn_hydro_data(args.sn_hydro)
        t_days = header['t_days']
        datasets.append(('SN_HYDRO_DATA', r_h, v_h, T_h, rho_h, mf_h, colors['SN_HYDRO_DATA']))
        print(f"SN_HYDRO_DATA: {len(r_h)} zones, t={t_days:.1f}d")

    # Read SEDONA .mod
    if args.sedona_mod and Path(args.sedona_mod).exists():
        r_s, v_s, rho_s, T_s, mf_s = read_sedona_mod(args.sedona_mod)
        datasets.append(('SEDONA', r_s, v_s, T_s, rho_s, mf_s, colors['SEDONA']))
        print(f"SEDONA .mod: {len(r_s)} zones")

    # Read SuperLite input.str
    if args.superlite_str and Path(args.superlite_str).exists():
        r_sl, v_sl, rho_sl, T_sl, mf_sl = read_superlite_str(args.superlite_str)
        datasets.append(('SuperLite', r_sl, v_sl, T_sl, rho_sl, mf_sl, colors['SuperLite']))
        print(f"SuperLite input.str: {len(r_sl)} cells")

    if not datasets:
        print("ERROR: No input files found")
        return

    # Plot structure
    for name, r, v, T, rho, mf, color in datasets:
        plot_structure(r, v, T, rho, name, color, ax_r, ax_v, ax_T, ax_rho)

    # Format structure plots
    ax_r.set_ylabel('Radius (R$_\\odot$)', fontsize=12)
    ax_r.set_yscale('log')
    ax_r.legend(fontsize=10)
    ax_r.set_title('Radius Profile', fontsize=13)
    ax_r.grid(True, alpha=0.3)

    ax_v.set_ylabel('Velocity (km/s)', fontsize=12)
    ax_v.legend(fontsize=10)
    ax_v.set_title('Velocity Profile', fontsize=13)
    ax_v.grid(True, alpha=0.3)

    ax_T.set_ylabel('Temperature (K)', fontsize=12)
    ax_T.set_yscale('log')
    ax_T.legend(fontsize=10)
    ax_T.set_title('Temperature Profile', fontsize=13)
    ax_T.grid(True, alpha=0.3)

    ax_rho.set_ylabel('Density (g cm$^{-3}$)', fontsize=12)
    ax_rho.set_yscale('log')
    ax_rho.set_xlabel('Radius (R$_\\odot$)', fontsize=12)
    ax_rho.legend(fontsize=10)
    ax_rho.set_title('Density Profile', fontsize=13)
    ax_rho.grid(True, alpha=0.3)

    title = f'Input Structure Comparison - {args.model_name}' if args.model_name else 'Input Structure Comparison'
    fig.suptitle(title, fontsize=15, y=1.01)
    plt.tight_layout()
    outpath = outdir / 'input_structure_comparison.png'
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {outpath}")

    # --- Mass fractions comparison plot ---
    fig, axes = plt.subplots(1, len(datasets), figsize=(7 * len(datasets), 8), sharey=True)
    if len(datasets) == 1:
        axes = [axes]

    for ax, (name, r, v, T, rho, mf, color) in zip(axes, datasets):
        r_solar = r / 6.96e10
        for elem, vals in mf.items():
            if np.any(vals > 1e-10):
                ax.plot(r_solar, vals, label=elem,
                        color=elem_colors.get(elem, 'gray'), linewidth=1.5)
        ax.set_xlabel('Radius (R$_\\odot$)', fontsize=12)
        ax.set_title(f'{name} Mass Fractions', fontsize=13)
        ax.set_yscale('log')
        ax.set_ylim(1e-6, 1)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=9, ncol=2)

    axes[0].set_ylabel('Mass Fraction', fontsize=12)
    title = f'Mass Fractions - {args.model_name}' if args.model_name else 'Mass Fractions'
    fig.suptitle(title, fontsize=15, y=1.01)
    plt.tight_layout()
    outpath = outdir / 'input_mass_fractions_comparison.png'
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {outpath}")

    # --- Per-code structure + mass fractions combined plot ---
    for name, r, v, T, rho, mf, color in datasets:
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        ax_r, ax_v, ax_rho = axes[0]
        ax_T, ax_mf, ax_mfr = axes[1]

        r_solar = r / 6.96e10
        v_kms = v / 1e5
        zone_idx = np.arange(len(r))

        # Radius vs zone
        ax_r.plot(zone_idx, r_solar, color=color, linewidth=1.5)
        ax_r.set_ylabel('Radius (R$_\\odot$)', fontsize=12)
        ax_r.set_xlabel('Zone Index', fontsize=12)
        ax_r.set_title('Radius vs Zone', fontsize=13)
        ax_r.set_yscale('log')
        ax_r.grid(True, alpha=0.3)

        # Velocity vs zone
        ax_v.plot(zone_idx, v_kms, color=color, linewidth=1.5)
        ax_v.set_ylabel('Velocity (km/s)', fontsize=12)
        ax_v.set_xlabel('Zone Index', fontsize=12)
        ax_v.set_title('Velocity vs Zone', fontsize=13)
        ax_v.grid(True, alpha=0.3)

        # Density vs zone
        ax_rho.plot(zone_idx, rho, color=color, linewidth=1.5)
        ax_rho.set_ylabel('Density (g cm$^{-3}$)', fontsize=12)
        ax_rho.set_xlabel('Zone Index', fontsize=12)
        ax_rho.set_title('Density vs Zone', fontsize=13)
        ax_rho.set_yscale('log')
        ax_rho.grid(True, alpha=0.3)

        # Temperature vs zone
        ax_T.plot(zone_idx, T, color=color, linewidth=1.5)
        ax_T.set_ylabel('Temperature (K)', fontsize=12)
        ax_T.set_xlabel('Zone Index', fontsize=12)
        ax_T.set_title('Temperature vs Zone', fontsize=13)
        ax_T.set_yscale('log')
        ax_T.grid(True, alpha=0.3)

        # Mass fractions vs zone
        for elem, vals in mf.items():
            if np.any(vals > 1e-10):
                ax_mf.plot(zone_idx, vals, label=elem,
                           color=elem_colors.get(elem, 'gray'), linewidth=1.5)
        ax_mf.set_ylabel('Mass Fraction', fontsize=12)
        ax_mf.set_xlabel('Zone Index', fontsize=12)
        ax_mf.set_title('Mass Fractions vs Zone', fontsize=13)
        ax_mf.set_yscale('log')
        ax_mf.set_ylim(1e-6, 1)
        ax_mf.legend(fontsize=9, ncol=3)
        ax_mf.grid(True, alpha=0.3)

        # Mass fractions vs radius (additional view)
        for elem, vals in mf.items():
            if np.any(vals > 1e-10):
                ax_mfr.plot(r_solar, vals, label=elem,
                            color=elem_colors.get(elem, 'gray'), linewidth=1.5)
        ax_mfr.set_ylabel('Mass Fraction', fontsize=12)
        ax_mfr.set_xlabel('Radius (R$_\\odot$)', fontsize=12)
        ax_mfr.set_title('Mass Fractions vs Radius', fontsize=13)
        ax_mfr.set_yscale('log')
        ax_mfr.set_ylim(1e-6, 1)
        ax_mfr.legend(fontsize=9, ncol=3)
        ax_mfr.grid(True, alpha=0.3)

        fig.suptitle(f'{name} Input Structure - {args.model_name}', fontsize=15, y=1.01)
        plt.tight_layout()
        safe_name = name.replace(' ', '_').replace('.', '_')
        outpath = outdir / f'{safe_name.lower()}_input_structure.png'
        plt.savefig(outpath, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"Saved: {outpath}")


if __name__ == '__main__':
    main()
