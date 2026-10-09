#!/usr/bin/env python3
"""Plot STELLA lightcurves for all models in a directory.

Searches a given directory (e.g., grid/) for all available STELLA lightcurves
and creates separate plots for each selected band. Saves plots to the same directory.

Usage:
    python plot_model_lightcurves.py grid/ --bands V R
    python plot_model_lightcurves.py models/ --bands U B V R I --days-max 60

Options:
    --bands BAND [BAND ...]  Bands to plot (default: V R)
    --days-max DAYS          Clip x-axis at DAYS
    --title TITLE            Figure title
"""

import argparse
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BAND_COLORS = {
    'U': '#7B00D4', 'B': '#0000FF', 'V': '#00AA00',
    'R': '#FF4400', 'I': '#AA0000',
    'u': '#7B00D4', 'g': '#228B22', 'r': '#FF4400',
    'i': '#AA0000', 'z': '#555555',
}


def find_lightcurve(model_dir: Path):
    """Find the lightcurve_complete CSV in a model directory."""
    output_dir = model_dir / 'output'
    if not output_dir.exists():
        return None
    # Look for *_lightcurve_complete.csv
    candidates = sorted(output_dir.glob('*_lightcurve_complete.csv'))
    if candidates:
        return candidates[0]
    # Fallback: lightcurve_complete.csv (old naming)
    fallback = output_dir / 'lightcurve_complete.csv'
    if fallback.exists():
        return fallback
    return None


def read_lightcurve(csv_path: Path):
    """Read a lightcurve CSV and return dict of band -> (times, mags)."""
    times = []
    rows = []
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            times.append(float(row['time']))
            rows.append(row)

    data = {}
    for band in BAND_COLORS:
        mags = []
        for row in rows:
            v = row.get(band, '')
            if v == '':
                mags.append(float('nan'))
            else:
                mags.append(float(v))
        data[band] = (times, mags)
    return data


def find_all_lightcurves(search_dir: Path):
    """Find all STELLA lightcurves in subdirectories of search_dir."""
    lightcurves = []
    
    if not search_dir.exists():
        print(f"ERROR: Directory {search_dir} does not exist", file=sys.stderr)
        return lightcurves
    
    # Search for subdirectories containing lightcurve files
    for model_dir in sorted(search_dir.iterdir()):
        if not model_dir.is_dir():
            continue
        
        csv_path = find_lightcurve(model_dir)
        if csv_path:
            lightcurves.append((model_dir, csv_path))
            print(f"Found lightcurve: {model_dir.name}")
    
    return lightcurves


def plot_band(lightcurves, band, output_dir, days_max=None, title=None):
    """Create a plot for a single band."""
    fig, ax = plt.subplots(figsize=(10, 6))
    
    for model_dir, csv_path in lightcurves:
        data = read_lightcurve(csv_path)
        
        if band not in data:
            print(f"WARNING: Band {band} not found in {model_dir.name}", file=sys.stderr)
            continue
        
        times, mags = data[band]
        color = BAND_COLORS.get(band, '#333333')
        label = model_dir.name
        
        # Plot this model's lightcurve
        ax.plot(times, mags, color=color, linestyle='-',
                label=label, linewidth=1.5, alpha=0.7)
    
    ax.invert_yaxis()
    ax.set_xlabel('Time since explosion (days)', fontsize=12)
    ax.set_ylabel(f'{band} absolute magnitude', fontsize=12)
    
    if days_max:
        ax.set_xlim(left=-5, right=days_max)
    else:
        ax.set_xlim(left=-5)
    
    plot_title = title or f'STELLA Lightcurve Comparison - {band} band'
    ax.set_title(plot_title, fontsize=13)
    ax.legend(loc='upper right', fontsize=8, ncol=2)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    
    # Save the plot
    output_path = output_dir / f'lightcurve_{band.lower()}.png'
    plt.savefig(output_path, dpi=200)
    print(f"Saved {band} band plot: {output_path}")
    plt.close()


def main():
    parser = argparse.ArgumentParser(description='Plot STELLA model lightcurves from directory')
    parser.add_argument('search_dir', type=Path,
                        help='Directory to search for models (e.g., grid/)')
    parser.add_argument('--bands', nargs='+', default=['V', 'R'],
                        help='Bands to plot (default: V R)')
    parser.add_argument('--days-max', type=float, default=None,
                        help='Max time in days to show')
    parser.add_argument('--title', type=str, default=None,
                        help='Figure title prefix')
    args = parser.parse_args()

    # Deduplicate requested bands while preserving order
    bands = []
    seen = set()
    for b in args.bands:
        b = b.strip()
        if b not in seen:
            bands.append(b)
            seen.add(b)

    # Find all lightcurves in the search directory
    print(f"Searching for lightcurves in {args.search_dir}...")
    lightcurves = find_all_lightcurves(args.search_dir)
    
    if not lightcurves:
        print(f"ERROR: No lightcurves found in {args.search_dir}", file=sys.stderr)
        sys.exit(1)
    
    print(f"Found {len(lightcurves)} lightcurve(s)")
    
    # Create a separate plot for each band
    for band in bands:
        plot_band(lightcurves, band, args.search_dir, args.days_max, args.title)
    
    print(f"\nSaved {len(bands)} plot(s) to {args.search_dir}/")


if __name__ == '__main__':
    main()