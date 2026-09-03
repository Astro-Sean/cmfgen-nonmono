#!/usr/bin/env python3
"""Generate overview plots from SN_HYDRO_DATA for CMFGEN model directory.

Called by submit_cmfgen.sh after CMFGEN directory preparation:
    python3 plot_sn_hydro_cmfgen_inputs.py <MODEL_DIR>

Reads <MODEL_DIR>/SN_HYDRO_DATA and writes PNGs to <MODEL_DIR>/plots/:
  - structure_overview.png: R, V, T, rho vs radius
  - mass_fractions.png: stacked mass fractions vs radius
"""
import sys
import os
import re
import numpy as np
from pathlib import Path


def parse_sn_hydro_data(filepath):
    with open(filepath) as f:
        lines = f.readlines()

    def find_section(name):
        for i, line in enumerate(lines):
            if line.strip().lower().startswith(name.lower()):
                return i
        return -1

    def read_array(start_idx, nd):
        vals = []
        i = start_idx + 1
        while len(vals) < nd:
            parts = lines[i].split()
            for p in parts:
                vals.append(float(p))
                if len(vals) >= nd:
                    break
            i += 1
        return vals

    nd = int(re.search(r'\d+', lines[0]).group())
    time_days = float(lines[3].split()[-1])
    lum_lsun = float(lines[4].split()[-1])

    r_idx = find_section('Radius grid')
    v_idx = find_section('Velocity')
    t_idx = find_section('Temperature')
    rho_idx = find_section('Density')
    ed_idx = find_section('Electron density')

    radius = read_array(r_idx, nd)
    velocity = read_array(v_idx, nd)
    temperature = read_array(t_idx, nd)
    density = read_array(rho_idx, nd)
    ed = read_array(ed_idx, nd) if ed_idx >= 0 else None

    mass_fracs = {}
    elem_labels = ['HYD', 'HE', 'OXY', 'SIL', 'CAL', 'NICK', 'IRON', 'NICK 56', 'COB 56']
    for label in elem_labels:
        idx = find_section(f'{label} mass fraction')
        if idx >= 0:
            mass_fracs[label] = read_array(idx, nd)

    return {
        'nd': nd, 'time_days': time_days, 'lum_lsun': lum_lsun,
        'radius': radius, 'velocity': velocity,
        'temperature': temperature, 'density': density,
        'ed': ed, 'mass_fracs': mass_fracs,
    }


def plot_structure_overview(hydro, output_file):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("WARNING: matplotlib not available; skipping plots")
        return

    RSUN_CM = 6.96e10
    r = np.array(hydro['radius'])[::-1] * 1e10 / RSUN_CM  # R_sun, inner->outer
    v = np.array(hydro['velocity'])[::-1]  # km/s
    t = np.array(hydro['temperature'])[::-1] * 1e4  # K
    rho = np.array(hydro['density'])[::-1]  # g/cm^3

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    fig.suptitle(f"SN_HYDRO_DATA Structure (t={hydro['time_days']:.1f}d, L={hydro['lum_lsun']:.2e} Lsun)",
                 fontsize=12)

    ax = axes[0, 0]
    ax.plot(r, v, 'b-', lw=1.5)
    ax.set_xlabel('Radius (R_sun)')
    ax.set_ylabel('Velocity (km/s)')
    ax.set_title('Velocity Profile')
    ax.set_xscale('log')

    ax = axes[0, 1]
    ax.plot(r, t, 'r-', lw=1.5)
    ax.set_xlabel('Radius (R_sun)')
    ax.set_ylabel('Temperature (K)')
    ax.set_title('Temperature Profile')
    ax.set_xscale('log')
    ax.set_yscale('log')

    ax = axes[1, 0]
    ax.plot(r, rho, 'g-', lw=1.5)
    ax.set_xlabel('Radius (R_sun)')
    ax.set_ylabel('Density (g/cm^3)')
    ax.set_title('Density Profile')
    ax.set_xscale('log')
    ax.set_yscale('log')

    ax = axes[1, 1]
    if hydro['ed'] is not None:
        ed = np.array(hydro['ed'])[::-1]
        ax.plot(r, ed, 'm-', lw=1.5)
        ax.set_ylabel('Electron Density (cm^-3)')
    else:
        ax.text(0.5, 0.5, 'No electron density data', transform=ax.transAxes, ha='center')
    ax.set_xlabel('Radius (R_sun)')
    ax.set_title('Electron Density')
    ax.set_xscale('log')
    ax.set_yscale('log')

    plt.tight_layout()
    fig.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Wrote {output_file}")


def plot_mass_fractions(hydro, output_file):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        return

    RSUN_CM = 6.96e10
    r = np.array(hydro['radius'])[::-1] * 1e10 / RSUN_CM

    elem_names = {
        'HYD': 'H', 'HE': 'He', 'OXY': 'O', 'SIL': 'Si',
        'CAL': 'Ca', 'NICK': 'Ni', 'IRON': 'Fe',
        'NICK 56': 'Ni56', 'COB 56': 'Co56',
    }

    fig, ax = plt.subplots(figsize=(10, 6))
    fig.suptitle(f"Mass Fractions (t={hydro['time_days']:.1f}d)", fontsize=12)

    colors = plt.cm.tab10(np.linspace(0, 1, len(hydro['mass_fracs'])))
    bottom = np.zeros(len(r))
    for (label, fracs), color in zip(hydro['mass_fracs'].items(), colors):
        f = np.array(fracs)[::-1]
        name = elem_names.get(label, label)
        ax.fill_between(r, bottom, bottom + f, label=name, alpha=0.7, color=color)
        bottom += f

    ax.set_xlabel('Radius (R_sun)')
    ax.set_ylabel('Mass Fraction')
    ax.set_xscale('log')
    ax.set_ylim(0, max(1.1, np.max(bottom)))
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=9)

    plt.tight_layout()
    fig.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Wrote {output_file}")


def main():
    if len(sys.argv) < 2:
        print("Usage: plot_sn_hydro_cmfgen_inputs.py <MODEL_DIR>")
        sys.exit(1)

    model_dir = Path(sys.argv[1])
    hydro_path = model_dir / 'SN_HYDRO_DATA'

    if not hydro_path.exists():
        print(f"WARNING: SN_HYDRO_DATA not found at {hydro_path}")
        sys.exit(0)

    hydro = parse_sn_hydro_data(hydro_path)
    plots_dir = model_dir / 'plots'
    plots_dir.mkdir(exist_ok=True)

    print(f"Writing SN_HYDRO_DATA overview plots -> {plots_dir}")
    plot_structure_overview(hydro, plots_dir / 'structure_overview.png')
    plot_mass_fractions(hydro, plots_dir / 'mass_fractions.png')


if __name__ == '__main__':
    main()
