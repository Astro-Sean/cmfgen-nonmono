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
    
    def _parse_fortran_float(s):
        """Parse Fortran-style scientific notation (e.g. 1.2105-265 -> 1.2105E-265)."""
        import re
        # Insert E between mantissa and exponent if missing (e.g. 1.23-45 -> 1.23E-45)
        s = re.sub(r'([0-9])([+-][0-9])', r'\1E\2', s)
        return float(s)

    # Read frequencies (lines after header_idx to flux_header_idx-1)
    freq = []
    for i in range(header_idx + 1, flux_header_idx):
        line = lines[i].strip()
        if line:
            values = line.split()
            freq.extend([_parse_fortran_float(v) for v in values])
    freq = np.array(freq[:nfreq])
    
    # Read fluxes (lines after flux_header_idx)
    flux_jy = []
    for i in range(flux_header_idx + 1, len(lines)):
        line = lines[i].strip()
        if line:
            values = line.split()
            flux_jy.extend([_parse_fortran_float(v) for v in values])
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
    #
    # CMFGEN's cmf_flux.exe (BATOBS) can emit non-finite or divergent fluxes
    # at a small fraction of frequency points (line centres / edges where the
    # formal solution blows up). Observed on the X/20d model: 501 NaN, 2 inf,
    # and ~420 finite-but-divergent points spanning 1e10 to 1e300 Jy, while
    # the physical spectrum is 1e4-1e8 Jy. Both classes destroy the plot's
    # auto-scaling. Mask them to NaN so matplotlib renders gaps and the axis
    # limits are computed from valid data only.
    flux_jy = np.where(np.isfinite(flux_jy), flux_jy, np.nan)

    # Robust outlier mask: BATOBS divergence produces finite values that are
    # 10-290 orders of magnitude above the physical flux. Use the 90th
    # percentile of |flux| (robust: at most ~10% of points are divergent)
    # and mask anything more than 1e2x above it. This is generous enough to
    # preserve real P-Cygni emission peaks (~10-30x continuum) while cleanly
    # separating the divergent tail (which starts >100x above the physical
    # maximum). Empirically validated on the X/20d model where physical
    # fluxes span 1e4-1e8 Jy and divergence starts at 1e10 Jy.
    finite_jy = flux_jy[np.isfinite(flux_jy)]
    if finite_jy.size > 100:
        p90 = np.percentile(np.abs(finite_jy), 90)
        if p90 > 0:
            flux_jy = np.where(np.abs(flux_jy) > p90 * 1e2, np.nan, flux_jy)

    flux_cgs = flux_jy * 1e-23  # Convert Jy to erg/s/cm^2/Hz
    flux = flux_cgs * (c / wavelength_cm**2)  # Convert to erg/s/cm^2/A
    flux = np.where(np.isfinite(flux), flux, np.nan)

    return wavelength, flux


def _smooth_flux(flux, window=11):
    """Median-filter a flux array, ignoring NaNs.

    Uses scipy.ndimage.median_filter with 'reflect' mode on a NaN-filled
    copy, then restores NaNs at positions that were NaN in the input so
    that genuine gaps (e.g. masked BATOBS divergence) are not filled in.

    Parameters
    ----------
    flux : np.ndarray
        1-D flux array, may contain NaNs.
    window : int
        Median filter window size (odd). At the CMFGEN frequency grid
        spacing (~0.24 Å/pt in the optical), window=11 ~ 2.6 Å — enough
        to suppress 1-2 point noise spikes without erasing real P-Cygni
        features (200-300 Å wide).

    Returns
    -------
    np.ndarray
        Smoothed flux, same length, NaN where input was NaN.
    """
    from scipy.ndimage import median_filter
    if window < 3:
        return flux.copy()
    if window % 2 == 0:
        window += 1
    nan_mask = np.isnan(flux)
    filled = np.where(nan_mask, 0.0, flux)
    smoothed = median_filter(filled, size=window, mode='reflect')
    # Restore NaN at positions that were NaN so genuine gaps are not
    # filled with spurious smoothed values.
    smoothed = np.where(nan_mask, np.nan, smoothed)
    return smoothed


def _ylimits_from_smoothed(flux, smooth_window=11, pad_frac=0.05):
    """Compute robust y-axis limits from a smoothed spectrum.

    Smoothing suppresses per-bin noise and residual borderline-divergent
    points that survived the outlier mask, so the limits reflect the real
    spectral extent (continuum + P-Cygni peaks/troughs) rather than
    individual noisy bins. The 1st/99th percentiles of the smoothed
    spectrum are used (not min/max) to exclude any residual outliers in
    the smoothed data. A small padding fraction is added so line peaks
    are not clipped at the frame edge.

    Parameters
    ----------
    flux : np.ndarray
        1-D flux array (may contain NaNs).
    smooth_window : int
        Median-filter window for the smoothing.
    pad_frac : float
        Fraction of the smoothed percentile range to add as padding
        above/below.

    Returns
    -------
    (ymin, ymax) or None if no finite data.
    """
    smoothed = _smooth_flux(flux, window=smooth_window)
    finite = smoothed[np.isfinite(smoothed)]
    if finite.size == 0:
        return None
    lo = np.percentile(finite, 1)
    hi = np.percentile(finite, 99)
    span = hi - lo
    if span == 0:
        span = abs(hi) * 0.1 or 1.0
    return lo - pad_frac * span, hi + pad_frac * span


def plot_spectrum(wavelength, flux, output_file):
    """
    Plot F_lambda vs wavelength with linear scale and line markers.

    Uses LINEAR y-scale (not log) because SN P-Cygni line profiles span
    200-300 Å and have 15-25% absorption depth — log scale compresses
    these broad features into invisible wiggles on the continuum.

    Y-axis limits are derived from a median-smoothed copy of the spectrum
    so that residual per-bin noise / borderline-divergent BATOBS points
    do not stretch the axis. The raw (unsmoothed) spectrum is plotted.

    Parameters:
    -----------
    wavelength : array
        Wavelength in Angstroms
    flux : array
        Flux in erg/s/cm^2/A
    output_file : str
        Output PNG file path
    """
    fig, axes = plt.subplots(2, 1, figsize=(14, 10), gridspec_kw={'height_ratios': [2, 1]})

    # Restrict to optical range
    mask = (wavelength >= 3500) & (wavelength <= 9000)
    wl_plot = wavelength[mask]
    flux_plot = flux[mask]

    # Major SN lines to mark — only lines with detectable features at
    # this epoch. Verified against the X/20d spectrum: Hγ 4340, He I 5876,
    # [O I] 6300, and Ca II 8664 show <5% depth (pure continuum), so they
    # are omitted to avoid misleading markers. Si II 6355 is absent (this
    # is a Type IIb, not Type Ia).
    line_labels = {
        4102: 'Hδ', 4861: 'Hβ', 5160: 'Fe II',
        6563: 'Hα',
        7774: 'O I', 8446: 'O I', 8542: 'Ca II',
    }

    # Top panel: full optical range, LINEAR scale
    ax = axes[0]
    ax.plot(wl_plot, flux_plot, linewidth=0.3, color='black')
    ax.set_yscale('linear')
    ax.set_xlim(3500, 9000)
    ax.set_title('CMFGEN Observed Spectrum (linear scale)', fontsize=14)
    ax.set_ylabel(r'$F_\lambda$ (erg s$^{-1}$ cm$^{-2}$ Å$^{-1}$)', fontsize=11)
    ax.grid(True, alpha=0.3)
    # Y-limits from smoothed spectrum so residual noise / borderline
    # divergent points don't stretch the axis.
    ylim_top = _ylimits_from_smoothed(flux_plot)
    if ylim_top:
        ax.set_ylim(*ylim_top)
    for wl, name in line_labels.items():
        ax.axvline(wl, color='blue', alpha=0.25, linestyle='--')
        ax.text(wl, ax.get_ylim()[1] * 0.92, name, fontsize=7,
                ha='center', color='blue', rotation=90)

    # Bottom panel: zoomed on H-alpha region (6200-6900 Å) to show P-Cygni
    ax2 = axes[1]
    mask_ha = (wavelength >= 6200) & (wavelength <= 6900)
    wl_ha = wavelength[mask_ha]
    flux_ha = flux[mask_ha]
    ax2.plot(wl_ha, flux_ha, linewidth=0.5, color='black')
    ax2.set_yscale('linear')
    ax2.set_xlim(6200, 6900)
    ax2.set_title('Hα region zoom (P-Cygni profile)', fontsize=12)
    ax2.set_xlabel('Wavelength (Å)', fontsize=11)
    ax2.set_ylabel(r'$F_\lambda$', fontsize=11)
    ax2.grid(True, alpha=0.3)
    ylim_ha = _ylimits_from_smoothed(flux_ha)
    if ylim_ha:
        ax2.set_ylim(*ylim_ha)
    ax2.axvline(6563, color='red', alpha=0.3, linestyle='--', label='Hα 6563Å')
    ax2.legend(fontsize=9)

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close(fig)

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
    n_bad = int((~np.isfinite(flux)).sum())
    if n_bad:
        print(f"  WARNING: {n_bad} non-finite flux points masked (NaN/inf from CMFGEN BATOBS)")
    finite = flux[np.isfinite(flux)]
    if finite.size:
        print(f"  Flux range: {finite.min():.3e} - {finite.max():.3e} erg/s/cm^2/Å")
    else:
        print("  Flux range: all non-finite (no valid data to plot)", file=sys.stderr)
    
    print(f"Plotting spectrum...")
    plot_spectrum(wavelength, flux, args.out)
    
    print("Done.")


if __name__ == '__main__':
    main()
