#!/usr/bin/env python3
"""
Plot CMFGEN observed spectrum (obs_fin_10 format).
Plots F_lambda in cgs units (erg/s/cm^2/A) vs wavelength in Angstroms.
"""

import argparse
import sys
from pathlib import Path

try:
    import numpy as np
    import matplotlib.pyplot as plt
except ImportError as e:
    print(f"ERROR: Missing required Python package: {e}", file=sys.stderr)
    print("Install with: pip install numpy matplotlib", file=sys.stderr)
    sys.exit(1)


def read_obs_fin_10(filename):
    """
    Read CMFGEN obs_fin_10 (OBSFRAME) ASCII format.
    
    The format contains:
    - Header line with number of frequencies
    - Frequency array (in Hz) - multiple values per line
    - Flux header line
    - Flux array (in Janskys) - multiple values per line
    
    Converts to wavelength (Angstroms) and F_lambda (erg/s/cm^2/A).
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
    # Format: " Continuum Frequencies ( N )"
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
    
    # Read frequencies (lines after header_idx to flux_header_idx-1)
    freq = []
    for i in range(header_idx + 1, flux_header_idx):
        line = lines[i].strip()
        if line:
            values = line.split()
            freq.extend([float(v) for v in values])
    freq = np.array(freq[:nfreq])
    
    # Read fluxes (lines after flux_header_idx)
    flux_jy = []
    for i in range(flux_header_idx + 1, len(lines)):
        line = lines[i].strip()
        if line:
            values = line.split()
            flux_jy.extend([float(v) for v in values])
    flux_jy = np.array(flux_jy[:nfreq])
    
    # Convert frequency (10^15 Hz) to wavelength (Angstroms)
    # CMFGEN stores frequencies in units of 10^15 Hz
    # lambda = c / nu
    c = 2.99792458e10  # speed of light in cm/s
    freq_hz = freq * 1e15  # Convert from 10^15 Hz to Hz
    wavelength_cm = c / freq_hz  # wavelength in cm
    wavelength = wavelength_cm * 1e8  # convert to Angstroms
    
    # Convert flux from Janskys to erg/s/cm^2/A
    # 1 Jansky = 1e-26 W/m^2/Hz = 1e-23 erg/s/cm^2/Hz
    # F_lambda = F_nu * |dnu/dlambda| = F_nu * c / lambda^2
    flux_cgs = flux_jy * 1e-23  # Convert Jy to erg/s/cm^2/Hz
    flux = flux_cgs * (c / wavelength_cm**2)  # Convert to erg/s/cm^2/A
    
    return wavelength, flux


def plot_spectrum(wavelength, flux, output_file):
    """
    Plot F_lambda vs wavelength.
    
    Parameters:
    -----------
    wavelength : array
        Wavelength in Angstroms
    flux : array
        Flux in erg/s/cm^2/A
    output_file : str
        Output PNG file path
    """
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Plot F_lambda vs wavelength
    ax.plot(wavelength, flux, linewidth=0.5, color='black')
    
    ax.set_xlabel('Wavelength (Å)', fontsize=12)
    ax.set_ylabel(r'$F_\lambda$ (erg s$^{-1}$ cm$^{-2}$ Å$^{-1}$)', fontsize=12)
    ax.set_title('CMFGEN Observed Spectrum', fontsize=14)
    
    # Use log scale for y-axis (flux typically spans many orders of magnitude)
    ax.set_yscale('log')
    
    # Set x-axis limits to full wavelength range
    ax.set_xlim(wavelength.min(), wavelength.max())
    
    # Add grid for better readability
    ax.grid(True, alpha=0.3, which='both')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"Plot saved to {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description='Plot CMFGEN observed spectrum (obs_fin_10 format)'
    )
    parser.add_argument('input_file', help='Input obs_fin_10 file')
    parser.add_argument('--out', '-o', default='spectrum.png', 
                       help='Output PNG file (default: spectrum.png)')
    
    args = parser.parse_args()
    
    input_path = Path(args.input_file)
    if not input_path.exists():
        print(f"ERROR: Input file not found: {args.input_file}", file=sys.stderr)
        sys.exit(1)
    
    print(f"Reading spectrum from {args.input_file}...")
    try:
        wavelength, flux = read_obs_fin_10(args.input_file)
    except Exception as e:
        print(f"ERROR: Failed to read obs_fin_10 file: {e}", file=sys.stderr)
        sys.exit(1)
    
    print(f"Read {len(wavelength)} wavelength points")
    print(f"  Wavelength range: {wavelength.min():.1f} - {wavelength.max():.1f} Å")
    print(f"  Flux range: {flux.min():.3e} - {flux.max():.3e} erg/s/cm^2/Å")
    
    print(f"Plotting spectrum...")
    plot_spectrum(wavelength, flux, args.out)
    
    print("Done.")


if __name__ == '__main__':
    main()
