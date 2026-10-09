#!/usr/bin/env python3
"""
Plot all STELLA lightcurves found in the models' output folders.

By default, plots only the R band. To plot additional bands,
modify the BANDS_TO_PLOT list below.

Usage:
    python plot_all_stella_lightcurves.py
    python plot_all_stella_lightcurves.py --output lightcurves.png
"""

import argparse
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ============================================
# CONFIGURATION: Modify these to change behavior
# ============================================

# Bands to plot (add more bands here to plot multiple bands)
# Available bands: U, B, V, R, I, u, g, r, i, z
BANDS_TO_PLOT = ['R']

# Band colors for plotting (used internally for reading all bands)
BAND_COLORS = {
    'U': '#7B00D4', 'B': '#0000FF', 'V': '#00AA00',
    'R': '#FF4400', 'I': '#AA0000',
    'u': '#7B00D4', 'g': '#228B22', 'r': '#FF4400',
    'i': '#AA0000', 'z': '#555555',
}

# Line styles for different bands (add more if plotting multiple bands)
BAND_LINESTYLES = {
    'R': '-',
    'V': '--',
    'B': ':',
    'I': '-.',
}

# Marker style for data points
MARKER = 'o'
MARKER_SIZE = 3
MARKER_ALPHA = 0.5
# Show markers every N points to avoid overcrowding
MARKEVERY = 10

# ============================================


def find_all_lightcurves(models_dir):
    """
    Find all STELLA lightcurve files in the models directory.
    Prefers native STELLA output (*_lightcurve.csv) over pysynphot (*_lightcurve_complete.csv).

    Parameters:
    -----------
    models_dir : Path
        Path to the models directory

    Returns:
    --------
    list of tuples
        (model_name, csv_path) for each found lightcurve
    """
    lightcurves = []
    models_path = Path(models_dir)

    if not models_path.exists():
        print(f"ERROR: Models directory not found: {models_dir}", file=sys.stderr)
        return lightcurves

    # Scan all subdirectories in models/
    for model_dir in sorted(models_path.iterdir()):
        if not model_dir.is_dir():
            continue

        output_dir = model_dir / 'output'
        if not output_dir.exists():
            continue

        # First, look for native STELLA output: *_lightcurve.csv
        candidates = sorted(output_dir.glob('*_lightcurve.csv'))
        # Filter out the _complete versions
        native_candidates = [c for c in candidates if '_complete' not in c.name]
        if native_candidates:
            lightcurves.append((model_dir.name, native_candidates[0]))
            continue

        # Fallback to pysynphot output: *_lightcurve_complete.csv
        pysyn_candidates = sorted(output_dir.glob('*_lightcurve_complete.csv'))
        if pysyn_candidates:
            lightcurves.append((model_dir.name, pysyn_candidates[0]))
            continue

        # Fallback: lightcurve_complete.csv (old naming)
        fallback = output_dir / 'lightcurve_complete.csv'
        if fallback.exists():
            lightcurves.append((model_dir.name, fallback))

    return lightcurves


def read_lightcurve(csv_path):
    """
    Read a lightcurve CSV and return dict of band -> (times, mags).
    Handles both native STELLA format and pysynphot format.

    Parameters:
    -----------
    csv_path : Path
        Path to the lightcurve CSV file

    Returns:
    --------
    dict
        band -> (times, mags) dictionary
    """
    times = []
    rows = []
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            times.append(float(row['time']))
            rows.append(row)

    # Detect format from column names
    # Native STELLA: time,U,B,V,R,I,bol_mag,bol_ergs
    # Pysynphot: time,U,B,V,R,I,u,g,r,i,z,bol_mag,bol_ergs
    columns = list(rows[0].keys()) if rows else []
    is_native_format = 'bol_ergs' in columns and 'u' not in columns

    data = {}
    # Use appropriate band set based on format
    if is_native_format:
        # Native STELLA bands
        available_bands = ['U', 'B', 'V', 'R', 'I']
    else:
        # Pysynphot bands
        available_bands = list(BAND_COLORS.keys())

    for band in available_bands:
        mags = []
        for row in rows:
            v = row.get(band, '')
            if v == '':
                mags.append(float('nan'))
            else:
                mags.append(float(v))
        data[band] = (times, mags)
    return data


def plot_lightcurves(lightcurves, bands, output_file=None, ylim=None):
    """
    Plot lightcurves for all models.

    Parameters:
    -----------
    lightcurves : list of tuples
        (model_name, csv_path) for each lightcurve
    bands : list
        List of bands to plot
    output_file : Path or None
        Output file path, or None to display
    ylim : tuple or None
        Y-axis limits (ymin, ymax), or None for auto
    """
    fig, ax = plt.subplots(figsize=(12, 8))

    # Use colormap to assign different colors to different models
    colormap = matplotlib.colormaps['tab20']
    n_models = len(lightcurves)

    for i, (model_name, csv_path) in enumerate(lightcurves):
        data = read_lightcurve(csv_path)

        # Assign color based on model index
        color = colormap(i / max(n_models - 1, 1))

        for band in bands:
            if band not in data:
                continue

            times, mags = data[band]
            linestyle = BAND_LINESTYLES.get(band, '-')

            # Use model name as label
            if len(bands) == 1:
                label = model_name
            else:
                label = f'{model_name} {band}'

            ax.plot(times, mags, color=color, linestyle=linestyle,
                    marker=MARKER, markersize=MARKER_SIZE, markevery=MARKEVERY,
                    label=label, linewidth=1.5, alpha=0.8)

    ax.invert_yaxis()
    ax.set_xlabel('Time since explosion (days)', fontsize=12)
    ax.set_ylabel('Absolute magnitude', fontsize=12)
    ax.set_xlim(left=-5)

    # Set y-axis limits (only if specified)
    if ylim is not None:
        ax.set_ylim(ylim)

    # Set title based on bands being plotted
    if len(bands) == 1:
        title = f'STELLA Lightcurves ({bands[0]} band)'
    else:
        title = f'STELLA Lightcurves ({", ".join(bands)} bands)'
    ax.set_title(title, fontsize=14)

    ax.legend(loc='upper right', fontsize=8, ncol=2)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=200)
        print(f"Saved figure: {output_file}")
    else:
        plt.show()
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description='Plot all STELLA lightcurves from models directory'
    )
    parser.add_argument('--models-dir', type=str, default='/ptmp/sbrennan/SN2024pba/models',
                        help='Path to models directory (default: /ptmp/sbrennan/SN2024pba/models)')
    parser.add_argument('--output', '-o', type=str, default=None,
                        help='Output PNG file (default: display)')
    parser.add_argument('--ylim', type=str, default=None,
                        help='Y-axis limits as "ymin,ymax" (default: -20,-17)')
    args = parser.parse_args()

    models_dir = Path(args.models_dir)
    print(f"Scanning {models_dir} for lightcurve files...")

    lightcurves = find_all_lightcurves(models_dir)

    if not lightcurves:
        print("No lightcurve files found.", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(lightcurves)} lightcurve(s):")
    for model_name, csv_path in lightcurves:
        print(f"  {model_name}: {csv_path}")

    print(f"\nPlotting bands: {', '.join(BANDS_TO_PLOT)}")
    print("To plot additional bands, modify BANDS_TO_PLOT in the script.")

    # Parse ylim if provided
    ylim = None
    if args.ylim:
        try:
            ylim = tuple(float(x) for x in args.ylim.split(','))
            if len(ylim) != 2:
                raise ValueError("ylim must have exactly 2 values")
            print(f"Y-axis limits: {ylim}")
        except Exception as e:
            print(f"ERROR: Invalid ylim format: {e}", file=sys.stderr)
            print("Use format: --ylim ymin,ymax (e.g., --ylim -20,-17)", file=sys.stderr)
            sys.exit(1)

    plot_lightcurves(lightcurves, BANDS_TO_PLOT, args.output, ylim)


if __name__ == '__main__':
    main()
