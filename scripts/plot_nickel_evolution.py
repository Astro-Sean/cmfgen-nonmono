#!/usr/bin/env python3
"""
Plot Nickel56 evolution from STELLA simulation.

Shows initial Ni56 distribution (grey dashed) from .xni file and
Ni56 evolution at multiple timesteps from STELLA .abn files.

Usage:
    python plot_nickel_evolution.py MODEL_DIR [OPTIONS]
    
Options:
    --xni-file PATH          Path to .xni file (default: stella/eve/MODEL.xni)
    --abn-dir PATH           Directory containing .abn files (default: stella/modmake)
    --n-epochs N             Number of epochs to plot (default: 10)
    --output PATH            Output plot path (default: nickel_evolution.png)
"""

import struct
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
import sys
import argparse
import glob


def read_xni(xni_path: str) -> tuple:
    """Read .xni binary file and return Ni56 mass fractions and metadata.
    
    Returns:
        tuple: (ni56_fractions, header_metadata)
    """
    with open(xni_path, 'rb') as f:
        data = f.read()
    
    # .xni format: 2 header floats + 500 zone pairs (value + placeholder)
    n_floats = len(data) // 4
    floats = struct.unpack('f' * n_floats, data)
    
    header = floats[:2]
    nzon = (n_floats - 2) // 2
    
    # Extract zone values (every other float starting at index 2)
    ni56_fractions = np.array([floats[2 + 2*z] for z in range(nzon)])
    
    return ni56_fractions, header


def read_abn(abn_path: str) -> dict:
    """Read STELLA .abn file for mass fractions.
    
    Format: km  dum  dum  dum  H  He  C  N  O  Ne  Na  Mg  Al  Si  S  Ar  Ca  FePeak  Ni58  Ni56
    
    Returns:
        dict with element mass fractions
    """
    data = []
    with open(abn_path) as f:
        for line in f:
            stripped = line.strip()
            if not stripped:
                continue
            cols = [float(x) for x in stripped.split()]
            if len(cols) >= 20:
                data.append(cols)
    
    arr = np.array(data)
    return {
        'km': arr[:, 0].astype(int),
        'H': arr[:, 4],
        'He': arr[:, 5],
        'C': arr[:, 6],
        'N': arr[:, 7],
        'O': arr[:, 8],
        'Ne': arr[:, 9],
        'Na': arr[:, 10],
        'Mg': arr[:, 11],
        'Al': arr[:, 12],
        'Si': arr[:, 13],
        'S': arr[:, 14],
        'Ar': arr[:, 15],
        'Ca': arr[:, 16],
        'Fe': arr[:, 17],
        'Ni58': arr[:, 18],
        'Ni56': arr[:, 19],
    }


def read_swd_for_mass_and_radius(swd_path: str, epoch_days: float) -> tuple:
    """Read STELLA .swd file and extract cumulative mass and radius data for a specific epoch.
    
    Returns (cumulative_mass_Msun, radius_cm) for each zone.
    """
    with open(swd_path) as f:
        for line in f:
            cols = line.split()
            if len(cols) < 13:
                continue
            
            t = float(cols[0])
            km = int(cols[1])
            
            # Find the requested epoch
            if km == 1 and t == epoch_days:
                # Read all zones for this epoch
                # Columns: t(days), zone, logM, logR, V(1e8 cm/s), logTg, logTr, logrho, logP, logQv, logeng, Flum, WRKX
                logM = []
                logR = []
                logM.append(float(cols[2]))  # log10(M_ext_Msun)
                logR.append(float(cols[3]))  # log10(R_cm)
                
                for line in f:
                    cols = line.split()
                    if len(cols) < 13:
                        continue
                    t_next = float(cols[0])
                    km_next = int(cols[1])
                    
                    if t_next == 0:  # Zone data
                        logM.append(float(cols[2]))
                        logR.append(float(cols[3]))
                    elif km_next == 1:  # Next epoch
                        break
                
                return 10**np.array(logM), 10**np.array(logR)  # Convert from log10 to linear
    
    # If exact epoch not found, return None
    return None, None


def find_abn_files(abn_dir: str) -> list:
    """Find all .abn files in directory sorted by time."""
    abn_files = sorted(Path(abn_dir).glob("*.abn"))
    return abn_files


def extract_time_from_abn_filename(abn_path: str) -> float:
    """Extract time in days from .abn filename."""
    # Format: X1_lam2-debug_T0.2000.abn or similar
    stem = Path(abn_path).stem
    try:
        # Try to extract time from common patterns
        if '_T' in stem:
            time_str = stem.split('_T')[-1].replace('.abn', '')
            return float(time_str)
        elif '_t' in stem:
            time_str = stem.split('_t')[-1].replace('.abn', '')
            return float(time_str)
        else:
            # Fallback: try to parse as float
            return float(stem.split('_')[-1])
    except (ValueError, IndexError):
        return 0.0


def plot_nickel_evolution(xni_path: str, abn_dir: str, swd_path: str,
                          n_epochs: int, output_path: str, model_name: str):
    """Create nickel evolution plot: Ni56 mass fraction vs radius over time.

    Reads the initial Ni56 distribution from the .xni file and the time-evolved
    Ni56 mass fractions from STELLA .abn output files, then plots them as a
    function of radius (from the .swd file) at several epochs.
    """

    # Read initial Ni56 distribution from .xni
    print(f"Reading initial Ni56 from {xni_path}")
    ni56_initial, header = read_xni(xni_path)
    nzon = len(ni56_initial)

    # Find .abn files
    abn_files = find_abn_files(abn_dir)
    print(f"Found {len(abn_files)} .abn files")

    if len(abn_files) == 0:
        print("ERROR: No .abn files found")
        return

    # Select epochs to plot
    if len(abn_files) > n_epochs:
        indices = np.linspace(0, len(abn_files) - 1, n_epochs, dtype=int)
        selected_abn = [abn_files[i] for i in indices]
    else:
        selected_abn = abn_files

    # Read mass and radius data from .swd for each epoch
    print(f"Reading mass and radius data from {swd_path}")
    mass_profiles = []
    radius_profiles = []
    ni56_profiles = []
    epoch_times = []

    for abn_file in selected_abn:
        time_days = extract_time_from_abn_filename(str(abn_file))
        epoch_times.append(time_days)
        mass, radius = read_swd_for_mass_and_radius(swd_path, time_days)
        if mass is not None:
            mass_profiles.append(mass)
            radius_profiles.append(radius)
        else:
            # Fallback to zone index
            mass_profiles.append(np.arange(nzon))
            radius_profiles.append(np.arange(nzon))

        # Read Ni56 mass fraction from .abn file
        abn_data = read_abn(str(abn_file))
        ni56_profiles.append(abn_data['Ni56'])

    # Create plot: Ni56 mass fraction vs radius
    fig, ax = plt.subplots(figsize=(10, 6))

    colors = plt.cm.viridis(np.linspace(0, 1, len(selected_abn)))

    # Plot Ni56 evolution from .abn files
    for i, (ni56, radius) in enumerate(zip(ni56_profiles, radius_profiles)):
        time_days = epoch_times[i]
        ax.plot(radius, ni56, color=colors[i], linewidth=2, alpha=0.8,
                label=f't = {time_days:.2f} d')

    # Plot initial Ni56 distribution (grey dashed) using first epoch's radius
    if radius_profiles:
        ax.plot(radius_profiles[0], ni56_initial, 'k--', linewidth=2, alpha=0.7,
                label='Initial (from .xni)')

    # Formatting
    ax.set_xlabel('Radius (cm)', fontsize=12)
    ax.set_ylabel('Ni56 Mass Fraction', fontsize=12)
    ax.set_title(f'Ni56 Evolution: {model_name}', fontsize=14)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(bottom=1e-10)
    ax.legend(loc='upper left', fontsize=10)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    print(f"Saved plot to {output_path}")


def main():
    parser = argparse.ArgumentParser(description='Plot Ni56 evolution from STELLA')
    parser.add_argument('model_dir', help='Model directory path')
    parser.add_argument('--xni-file', help='Path to .xni file')
    parser.add_argument('--abn-dir', help='Directory containing .abn files')
    parser.add_argument('--swd-file', help='Path to .swd file for radius data')
    parser.add_argument('--n-epochs', type=int, default=10, 
                       help='Number of epochs to plot (default: 10)')
    parser.add_argument('--output', default='nickel_evolution.png',
                       help='Output plot path (default: nickel_evolution.png)')
    
    args = parser.parse_args()
    
    model_dir = Path(args.model_dir)
    model_name = model_dir.name
    
    # Set default paths
    xni_path = args.xni_file or str(model_dir / 'stella' / 'eve' / f'{model_name}.xni')
    abn_dir = args.abn_dir or str(model_dir / 'stella' / 'modmake')
    swd_path = args.swd_file or str(model_dir / 'stella' / f'{model_name}.swd')
    
    # Check if files exist
    if not Path(xni_path).exists():
        print(f"ERROR: .xni file not found: {xni_path}")
        sys.exit(1)
    
    if not Path(abn_dir).exists():
        print(f"ERROR: .abn directory not found: {abn_dir}")
        sys.exit(1)
    
    if not Path(swd_path).exists():
        print(f"WARNING: .swd file not found: {swd_path}")
        print("Will use zone index instead of radius")
    
    plot_nickel_evolution(xni_path, abn_dir, swd_path, args.n_epochs, 
                          args.output, model_name)


if __name__ == '__main__':
    main()
