#!/usr/bin/env python3
"""
Plot composition profile from STELLA .abn file and evolution from .swd file.

This script creates a multi-panel plot showing the radial profile of mass fractions
from the STELLA .abn file (initial stratified composition) and optionally shows
time evolution from the .swd file.

Usage:
    python plot_composition_evolution.py <abn_file> <output_dir> [OPTIONS]
    
Options:
    --swd_file FILE            Optional .swd file for time evolution
    --elements LIST             Comma-separated list of elements to plot (default: H,He,O,Si,Ca,Fe)
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
from pathlib import Path
import sys
import argparse


def read_swd_all_epochs(swd_path: str) -> dict:
    """Read all epochs from STELLA .swd file with composition data.
    
    Returns dict with epoch times as keys, each containing zone-by-zone data.
    Columns: t(days), zone, logM, logR, V(1e8 cm/s), logTg, logTr, logrho, logP, logQv, logeng, Flum(1e40), WRKX, n_bar, n_e, composition(15 mass fractions)
    Composition columns (15-29): H, He, C, N, O, Ne, Na, Mg, Al, Si, S, Ar, Ca, Fe, Ni
    NOTE: The .swd file writes Y_k*A_k (number fraction × atomic mass), not true mass fractions.
    These are normalized to true mass fractions (X_k = Y_k*A_k / <A>) when plotting.
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
            
            # Zone 1 marks the start of a new epoch
            if zone == 1:
                if current_t is not None and len(buf_data) > 0:
                    epochs[current_t] = np.array(buf_data)
                current_t = t
                buf_data = []
            
            # Parse zone data
            try:
                data = [float(cols[i]) for i in range(len(cols))]
                buf_data.append(data)
            except ValueError:
                continue
    
    # Don't forget the last epoch
    if current_t is not None and len(buf_data) > 0:
        epochs[current_t] = np.array(buf_data)
    
    return epochs


def read_abn(abn_path: str) -> dict:
    """Read STELLA *.abn file with composition data.
    
    Returns dict with element mass fractions per zone.
    Columns: km, dum, dum, dum, H, He, C, N, O, Ne, Na, Mg, Al, Si, S, Ar, Ca, FePeak, Ni58, Ni56
    """
    data = []
    with open(abn_path) as f:
        for line in f:
            if line.strip().startswith('#'):
                continue
            cols = [float(x) for x in line.split()]
            if len(cols) >= 20:
                data.append(cols)
    
    arr = np.array(data)
    # .abn file columns: km, dum, dum, dum, H, He, C, N, O, Ne, Na, Mg, Al, Si, S, Ar, Ca, FePeak, Ni58, Ni56
    # We need columns 4-19 (0-indexed: 4-19)
    return {
        'H':     arr[:, 4],
        'He':    arr[:, 5],
        'C':     arr[:, 6],
        'N':     arr[:, 7],
        'O':     arr[:, 8],
        'Ne':    arr[:, 9],
        'Na':    arr[:, 10],
        'Mg':    arr[:, 11],
        'Al':    arr[:, 12],
        'Si':    arr[:, 13],
        'S':     arr[:, 14],
        'Ar':    arr[:, 15],
        'Ca':    arr[:, 16],
        'FePeak': arr[:, 17],
        'Ni58':  arr[:, 18],
        'Ni56':  arr[:, 19],
        'zone':  arr[:, 0].astype(int),
    }


def plot_composition_profile(abn_data: dict, output_file: str, elements: list = None,
                              swd_data: dict = None, n_epochs_max: int = 10):
    """Create multi-panel plot showing radial profile of mass fractions.

    Uses dual x-axes:
      - Bottom: post-explosion radius from .swd snapshots (10^15 cm)
      - Top: initial pre-explosion radius from .swd t=0 epoch (10^11 cm)
    The .abn initial composition is plotted on the top axis, .swd epoch
    compositions on the bottom axis.  Both share the y-axis (mass fraction).
    """

    # Element names in STELLA YABUN array (order matches .swd columns 15-29)
    swd_element_names = ['H', 'He', 'C', 'N', 'O', 'Ne', 'Na', 'Mg', 'Al',
                         'Si', 'S', 'Ar', 'Ca', 'Fe', 'Ni']
    element_indices = {elem: i for i, elem in enumerate(swd_element_names)}

    # Use only the elements specified
    if elements is None:
        elements = ['H', 'He', 'C', 'O', 'Si', 'Ca', 'Fe']

    # Map element names for plotting (.abn uses 'FePeak', .swd uses 'Fe')
    element_map = {'Fe': 'FePeak'}

    # Filter to only elements that exist in the .abn data
    valid_elements = []
    display_names = []
    for element in elements:
        data_element = element_map.get(element, element)
        if data_element in abn_data:
            valid_elements.append(data_element)
            display_names.append(element)
        else:
            print(f"WARNING: Element {element} not found in .abn data, skipping")

    # --- Determine radii -------------------------------------------------
    if swd_data:
        times_sorted = sorted(swd_data.keys())
        # Initial (pre-explosion) radius from the first epoch (t≈0)
        first_epoch = swd_data[times_sorted[0]]
        r_initial = 10**first_epoch[:, 3]  # cm
        # Post-explosion radius from the last epoch
        last_epoch = swd_data[times_sorted[-1]]
        r_final = 10**last_epoch[:, 3]  # cm
        # Zone numbers for matching
        zones_initial = first_epoch[:, 1].astype(int)
        zones_final = last_epoch[:, 1].astype(int)
    else:
        r_initial = None
        r_final = None

    # --- Create figure ---------------------------------------------------
    n_elements = len(valid_elements)
    n_cols = 3
    n_rows = (n_elements + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(16, 4.5 * n_rows))
    if n_rows == 1:
        axes = axes.reshape(1, -1)

    if swd_data:
        title = 'Composition: Initial (top axis) vs Post-Explosion Snapshots (bottom axis)'
    else:
        title = 'Radial Profile of Mass Fractions (Initial Composition)'
    fig.suptitle(title, fontsize=14, fontweight='bold')

    # --- Plot each element ----------------------------------------------
    for i, (data_element, display_name) in enumerate(zip(valid_elements, display_names)):
        row = i // n_cols
        col = i % n_cols
        ax_bot = axes[row, col]

        # Determine the .swd element index
        swd_elem_name = 'Fe' if data_element == 'FePeak' else data_element
        elem_idx = element_indices.get(swd_elem_name, None)

        if swd_data and elem_idx is not None:
            # --- Bottom axis: .swd epoch compositions vs post-explosion radius
            times = times_sorted
            if len(times) > n_epochs_max:
                idx_sel = np.linspace(0, len(times) - 1, n_epochs_max, dtype=int)
                times = [times[j] for j in idx_sel]
            colors = plt.cm.viridis(np.linspace(0, 1, len(times)))

            for j, t in enumerate(times):
                epoch_data = swd_data[t]
                r = 10**epoch_data[:, 3]  # cm
                raw_comp = epoch_data[:, 15:30]
                comp_sum = raw_comp.sum(axis=1, keepdims=True)
                comp_sum = np.where(comp_sum < 1e-30, 1.0, comp_sum)
                mf = (raw_comp / comp_sum)[:, elem_idx]
                ax_bot.plot(r, mf, color=colors[j], linewidth=1.5, alpha=0.8,
                            linestyle='-', label=f'{t:.1f} d')

            ax_bot.set_xscale('log')
            ax_bot.set_xlabel('Post-explosion R (cm)', fontsize=9)
            r_final_val = 10**last_epoch[:, 3]
            ax_bot.set_xlim(r_final_val.min() * 0.7, r_final_val.max() * 1.3)

            # --- Top axis: initial .abn composition vs pre-explosion radius
            ax_top = ax_bot.twiny()
            # Match .abn zones to initial .swd zones by zone number
            abn_zones = abn_data['zone']
            # Build radius array matching .abn zone order
            r_init_by_zone = {}
            for kz, rz in zip(zones_initial, r_initial):
                r_init_by_zone[kz] = rz
            r_abn = np.array([r_init_by_zone.get(z, np.nan) for z in abn_zones])
            # Filter out any zones without radius match
            mask = ~np.isnan(r_abn)
            r_abn_plot = r_abn[mask]  # cm (pre-explosion scale)
            mf_abn = abn_data[data_element][mask]

            ax_top.plot(r_abn_plot, mf_abn, linewidth=2, color='gray',
                        linestyle='--', alpha=0.8, label='Initial (.abn)')
            ax_top.set_xscale('log')
            ax_top.set_xlabel('Initial R (cm)', fontsize=9, color='gray')
            r_init_val = r_initial
            ax_top.set_xlim(r_init_val.min() * 0.7, r_init_val.max() * 1.3)
            ax_top.tick_params(axis='x', colors='gray', labelsize=8)

            # Legend on first panel (combine both axes)
            if i == 0:
                lines_bot, labels_bot = ax_bot.get_legend_handles_labels()
                lines_top, labels_top = ax_top.get_legend_handles_labels()
                ax_bot.legend(lines_bot + lines_top, labels_bot + labels_top,
                              fontsize=7, loc='upper right')
        else:
            # No .swd data — just plot .abn vs zone number
            ax_bot.plot(abn_data['zone'], abn_data[data_element],
                        linewidth=2, color='gray', linestyle='--', alpha=0.8,
                        label='Initial (.abn)')
            ax_bot.set_xlabel('Zone Number (inner → outer)', fontsize=9)
            if i == 0:
                ax_bot.legend(fontsize=8, loc='upper right')

        ax_bot.set_ylabel(f'{display_name} mass fraction', fontsize=10)
        ax_bot.set_ylim(0, 1.1)
        ax_bot.grid(True, alpha=0.3)
        ax_bot.set_title(display_name, fontsize=12, fontweight='bold')
        ax_bot.tick_params(labelsize=8)

    # Hide unused subplots
    for i in range(n_elements, n_rows * n_cols):
        row = i // n_cols
        col = i % n_cols
        axes[row, col].set_visible(False)

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Composition plot saved to: {output_file}")
    plt.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('abn_file', help='STELLA .abn file path')
    parser.add_argument('output_dir', nargs='?', default='.', help='Output directory for plots')
    parser.add_argument('--swd_file', help='Optional STELLA .swd file for time evolution')
    parser.add_argument('--elements', help='Comma-separated list of elements to plot (e.g., H,He,O,Si,Ca,Fe)')
    parser.add_argument('--n-epochs', type=int, default=10, help='Max number of .swd epochs to plot (default: 10)')
    
    args = parser.parse_args()
    
    # Parse elements if provided
    elements = None
    if args.elements:
        elements = [e.strip() for e in args.elements.split(',')]
    
    # Read .abn file
    print(f"Reading .abn file: {args.abn_file}")
    abn_data = read_abn(args.abn_file)
    print(f"  Found {len(abn_data['zone'])} zones")
    print(f"  Available elements: {list(abn_data.keys())}")
    
    # Read .swd file if provided
    swd_data = None
    if args.swd_file:
        print(f"Reading .swd file: {args.swd_file}")
        swd_data = read_swd_all_epochs(args.swd_file)
        print(f"  Found {len(swd_data)} epochs")
        print(f"  Epoch times: {sorted(swd_data.keys())}")
    
    # Create output directory if needed
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate output filename
    abn_name = Path(args.abn_file).stem
    if swd_data:
        output_file = output_dir / f"{abn_name}_composition_evolution.png"
    else:
        output_file = output_dir / f"{abn_name}_composition_profile.png"
    
    # Create plot
    plot_composition_profile(abn_data, str(output_file), elements, swd_data, n_epochs_max=args.n_epochs)


if __name__ == '__main__':
    main()
