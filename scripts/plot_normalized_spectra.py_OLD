#!/usr/bin/env python3
"""
Plot multiple CMFGEN spectra normalized to 1 with vertical offsets.
Creates presentation-ready plots for PowerPoint.
"""

import argparse
import sys
from pathlib import Path

try:
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib import rcParams
except ImportError as e:
    print(f"ERROR: Missing required Python package: {e}", file=sys.stderr)
    print("Install with: pip install numpy matplotlib", file=sys.stderr)
    sys.exit(1)

# Set publication-quality style for PowerPoint
rcParams.update({
    'font.size': 14,
    'axes.labelsize': 16,
    'axes.titlesize': 18,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'legend.fontsize': 14,
    'figure.figsize': (14, 10),
    'lines.linewidth': 1.5,
})


def read_obs_fin_10(filename):
    """
    Read CMFGEN obs_fin_10 (OBSFRAME) ASCII format.
    
    Returns:
        wavelength (array): Wavelength in Angstroms
        flux (array): Flux in erg/s/cm^2/A
    """
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # Skip leading blank lines to find header
    header_idx = 0
    for i, line in enumerate(lines):
        if 'Continuum Frequencies' in line:
            header_idx = i
            break
    
    # Parse header to get number of points
    header_line = lines[header_idx].strip()
    nfreq = int(header_line.split('(')[1].split(')')[0].strip())
    
    # Find the flux header line
    flux_header_idx = None
    for i, line in enumerate(lines):
        if 'Observed intensity' in line:
            flux_header_idx = i
            break
    
    if flux_header_idx is None:
        raise ValueError("Could not find 'Observed intensity' header line")
    
    # Read frequencies
    freq = []
    for i in range(header_idx + 1, flux_header_idx):
        line = lines[i].strip()
        if line:
            values = line.split()
            freq.extend([float(v) for v in values])
    freq = np.array(freq[:nfreq])
    
    # Read fluxes
    flux_jy = []
    for i in range(flux_header_idx + 1, len(lines)):
        line = lines[i].strip()
        if line:
            values = line.split()
            flux_jy.extend([float(v) for v in values])
    flux_jy = np.array(flux_jy[:nfreq])
    
    # Convert frequency (10^15 Hz) to wavelength (Angstroms)
    c = 2.99792458e10  # speed of light in cm/s
    freq_hz = freq * 1e15
    wavelength_cm = c / freq_hz
    wavelength = wavelength_cm * 1e8  # convert to Angstroms
    
    # Convert flux from Janskys to erg/s/cm^2/A
    flux_cgs = flux_jy * 1e-23
    flux = flux_cgs * (c / wavelength_cm**2)
    
    return wavelength, flux


def normalize_flux(flux):
    """Normalize flux to maximum of 1."""
    if np.max(flux) > 0:
        return flux / np.max(flux)
    return flux


def offset_spectra(wavelengths, fluxes, offset_step=1.2):
    """
    Offset spectra vertically for clear visibility.
    
    Parameters:
        wavelengths (list): List of wavelength arrays
        fluxes (list): List of flux arrays
        offset_step (float): Vertical offset between spectra
    
    Returns:
        offset_fluxes (list): Flux arrays with vertical offsets
    """
    offset_fluxes = []
    for i, flux in enumerate(flux):
        # Normalize first
        norm_flux = normalize_flux(flux)
        # Add vertical offset
        offset_flux = norm_flux + i * offset_step
        offset_fluxes.append(offset_flux)
    
    return offset_fluxes


def plot_normalized_spectra(wavelengths, fluxes, labels, output_file, 
                          wavelength_range=None, offset_step=1.2):
    """
    Plot multiple spectra normalized to 1 with vertical offsets.
    
    Parameters:
        wavelengths (list): List of wavelength arrays
        fluxes (list): List of flux arrays
        labels (list): List of labels for legend
        output_file (str): Output PNG file path
        wavelength_range (tuple): (min_wl, max_wl) in Angstroms
        offset_step (float): Vertical offset between spectra
    """
    fig, ax = plt.subplots(figsize=(14, 10))
    
    # Offset and normalize spectra
    offset_fluxes = offset_spectra(fluxes, offset_step)
    
    # Plot each spectrum
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', 
              '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
    
    for i, (wl, fl, label) in enumerate(zip(wavelengths, offset_fluxes, labels)):
        color = colors[i % len(colors)]
        ax.plot(wl, fl, linewidth=1.5, color=color, label=label)
    
    # Set wavelength range if specified
    if wavelength_range:
        ax.set_xlim(wavelength_range[0], wavelength_range[1])
    else:
        # Use full range from all spectra
        all_wl = np.concatenate(wavelengths)
        ax.set_xlim(all_wl.min(), all_wl.max())
    
    # Set y-axis to accommodate all offsets
    max_offset = len(fluxes) * offset_step
    ax.set_ylim(-0.1, max_offset + 1.0)
    
    # Remove y-axis ticks (normalized flux)
    ax.set_yticks([])
    ax.set_ylabel('Normalized Flux (offset)', fontsize=16)
    
    # Add wavelength axis
    ax.set_xlabel('Wavelength (Å)', fontsize=16)
    
    # Add legend
    ax.legend(loc='upper right', framealpha=0.95, shadow=True)
    
    # Add grid for better readability
    ax.grid(True, alpha=0.3, which='both', linestyle='--')
    
    # Set background to white for PowerPoint
    ax.set_facecolor('white')
    fig.patch.set_facecolor('white')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    
    print(f"Plot saved to {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description='Plot multiple CMFGEN spectra normalized to 1 with vertical offsets'
    )
    parser.add_argument('input_files', nargs='+', help='Input obs_fin_10 files')
    parser.add_argument('--labels', nargs='+', help='Labels for each spectrum (must match number of input files)')
    parser.add_argument('--out', '-o', default='normalized_spectra.png', 
                       help='Output PNG file (default: normalized_spectra.png)')
    parser.add_argument('--wavelength-range', nargs=2, type=float, metavar=('MIN', 'MAX'),
                       help='Wavelength range in Angstroms (e.g., --wavelength-range 3000 7000)')
    parser.add_argument('--offset-step', type=float, default=1.2,
                       help='Vertical offset between spectra (default: 1.2)')
    
    args = parser.parse_args()
    
    # Validate inputs
    if args.labels and len(args.labels) != len(args.input_files):
        print(f"ERROR: Number of labels ({len(args.labels)}) must match number of input files ({len(args.input_files)})", 
              file=sys.stderr)
        sys.exit(1)
    
    # Read all spectra
    wavelengths = []
    fluxes = []
    labels = args.labels if args.labels else [Path(f).stem for f in args.input_files]
    
    for input_file in args.input_files:
        input_path = Path(input_file)
        if not input_path.exists():
            print(f"ERROR: Input file not found: {input_file}", file=sys.stderr)
            sys.exit(1)
        
        print(f"Reading spectrum from {input_file}...")
        try:
            wl, fl = read_obs_fin_10(input_file)
            wavelengths.append(wl)
            fluxes.append(fl)
            print(f"  Read {len(wl)} wavelength points")
            print(f"  Wavelength range: {wl.min():.1f} - {wl.max():.1f} Å")
        except Exception as e:
            print(f"ERROR: Failed to read {input_file}: {e}", file=sys.stderr)
            sys.exit(1)
    
    # Plot
    wavelength_range = tuple(args.wavelength_range) if args.wavelength_range else None
    plot_normalized_spectra(wavelengths, fluxes, labels, args.out, wavelength_range, args.offset_step)
    
    print("Done.")


if __name__ == '__main__':
    main()
