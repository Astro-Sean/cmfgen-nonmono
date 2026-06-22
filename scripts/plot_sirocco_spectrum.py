#!/usr/bin/env python3
"""
Plot SIROCCO output spectrum.

Reads the .spec file produced by SIROCCO and creates a publication-quality
plot of the synthetic spectrum compared to observed data (if available).

Usage:
    python plot_sirocco_spectrum.py <spec_file> [observed_file]

Example:
    python plot_sirocco_spectrum.py X1.spec
    python plot_sirocco_spectrum.py X1.spec SN2024pba_observed.txt
"""

import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

# Distance to SN2024pba for scaling
D_MPC = 83.58
D_CM = D_MPC * 1e6 * 3.085677581e18  # cm
FLUX_TO_LUM = 4 * np.pi * D_CM**2      # multiply observed flux by this to get L_λ

def read_sirocco_spec(spec_file):
    """
    Read SIROCCO .spec file.

    Returns dictionaries with:
    - freq: frequency (Hz)
    - wave: wavelength (Angstrom)
    - total: total spectrum (erg/s/cm^2/A at reference distance)
    - components: dict of individual components (Created, Emitted, Wind, etc.)
    - angles: dict of angle-dependent spectra
    """
    with open(spec_file, 'r') as f:
        lines = f.readlines()

    # Parse header to find column names and data start
    header_lines = []
    data_start = 0
    col_names = []

    for i, line in enumerate(lines):
        if line.startswith('#'):
            header_lines.append(line)
        elif data_start == 0 and not line.startswith('#'):
            # First non-comment line with data
            data_start = i
            parts = line.split()
            if len(parts) >= 2:
                # Try to parse the previous line as column headers
                if i > 0:
                    prev_line = lines[i-1].strip()
                    if prev_line.startswith('#'):
                        col_names = prev_line[1:].split()
                    else:
                        col_names = ['Freq.', 'Lambda', 'Created', 'WCreated',
                                     'Emitted', 'CenSrc', 'Disk', 'Wind',
                                     'HitSurf', 'Scattered']
                break

    if data_start == 0:
        raise ValueError(f"Could not find data in {spec_file}")

    # Parse data
    data = []
    for line in lines[data_start:]:
        if line.strip():
            try:
                vals = [float(x) for x in line.split()]
                data.append(vals)
            except ValueError:
                pass

    data = np.array(data)

    if len(data) == 0:
        raise ValueError(f"No data found in {spec_file}")

    # Extract columns
    ncols = data.shape[1]
    result = {
        'freq': data[:, 0],      # Hz
        'wave': data[:, 1],      # Angstrom
    }

    # Standard columns (up to 10)
    standard_cols = ['Created', 'WCreated', 'Emitted', 'CenSrc',
                     'Disk', 'Wind', 'HitSurf', 'Scattered']
    for i, name in enumerate(standard_cols):
        if i + 2 < ncols:
            result[name] = data[:, i + 2]

    # Angle-dependent columns (from column 10 onwards)
    result['angles'] = {}
    if ncols > 10:
        for i in range(10, ncols):
            angle_name = f"A{i-10}"
            if i - 10 < len(col_names) - 10:
                angle_name = col_names[i] if i < len(col_names) else f"A{i-10}"
            result['angles'][angle_name] = data[:, i]

    # Total spectrum (sum of all components at the reference angle)
    if 'Emitted' in result:
        result['total'] = result['Emitted']
    elif 'Created' in result:
        result['total'] = result['Created']
    else:
        result['total'] = np.zeros_like(result['freq'])

    return result


def load_observed(obs_file):
    """Load observed spectrum: wavelength(A), flux(erg/s/cm2/A), error."""
    wave, flux, err = np.loadtxt(obs_file, unpack=True)
    idx = np.argsort(wave)
    return wave[idx], flux[idx], err[idx]


def plot_spectrum(spec_file, obs_file=None, output_file=None):
    """Create spectrum plot."""
    spec_file = Path(spec_file)
    if not spec_file.exists():
        print(f"ERROR: Spectrum file not found: {spec_file}")
        return False

    # Read SIROCCO spectrum
    spec = read_sirocco_spec(spec_file)
    wave = spec['wave']
    total = spec['total']

    # Scale from flux at reference distance to intrinsic luminosity
    # SIROCCO spectra are at 100 pc by default
    ref_distance_pc = 100.0
    ref_distance_cm = ref_distance_pc * 3.085677581e18
    scale = (ref_distance_cm / D_CM) ** 2
    total_scaled = total / scale  # Convert to intrinsic luminosity

    fig, ax = plt.subplots(figsize=(12, 6))

    # Plot SIROCCO synthetic spectrum
    ax.plot(wave, total_scaled, color='C0', linewidth=1.0, label='SIROCCO', zorder=3)

    # Plot components if available
    for comp_name, color in [('Wind', 'C1'), ('CenSrc', 'C2'), ('Disk', 'C3')]:
        if comp_name in spec and np.any(spec[comp_name] > 0):
            comp_scaled = spec[comp_name] / scale
            ax.plot(wave, comp_scaled, color=color, linewidth=0.7,
                    alpha=0.7, label=comp_name, zorder=2)

    # Plot angle-dependent spectra if available
    if spec['angles']:
        for i, (angle_name, angle_spec) in enumerate(sorted(spec['angles'].items())):
            if i < 3:  # Only plot first 3 angles to avoid clutter
                angle_scaled = angle_spec / scale
                ax.plot(wave, angle_scaled, linewidth=0.5, alpha=0.5,
                        label=f'{angle_name}', zorder=1)

    # Overlay observed spectrum if provided
    if obs_file and Path(obs_file).exists():
        try:
            obs_wave, obs_flux, obs_err = load_observed(obs_file)
            mask = (obs_wave >= 3000) & (obs_wave <= 10000)
            obs_wave = obs_wave[mask]
            obs_flux = obs_flux[mask]
            obs_err = obs_err[mask]

            obs_flux_scaled = obs_flux * FLUX_TO_LUM
            obs_err_scaled = obs_err * FLUX_TO_LUM

            ax.step(obs_wave, obs_flux_scaled, where='mid', color='black',
                    linewidth=0.8, label=f'SN2024pba (scaled, d={D_MPC} Mpc)', zorder=4)
            ax.fill_between(obs_wave, obs_flux_scaled - obs_err_scaled,
                            obs_flux_scaled + obs_err_scaled,
                            color='black', alpha=0.2, step='mid', zorder=2)
        except Exception as e:
            print(f"Warning: Could not read observed spectrum: {e}")

    ax.set_xlabel(r'Wavelength ($\AA$)')
    ax.set_ylabel(r'Luminosity (erg s$^{-1}$ $\AA^{-1}$)')
    ax.set_xlim(1000, 25000)
    ax.set_title(f'SIROCCO Spectrum: {spec_file.stem}')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper right', fontsize='small')

    plt.tight_layout()

    if output_file is None:
        output_file = spec_file.with_suffix('_spectrum.png')
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Spectrum plot saved to {output_file}")
    plt.close()

    return True


def main():
    if len(sys.argv) < 2:
        print("Usage: python plot_sirocco_spectrum.py <spec_file> [observed_file]")
        sys.exit(1)

    spec_file = sys.argv[1]
    obs_file = sys.argv[2] if len(sys.argv) > 2 else None

    try:
        plot_spectrum(spec_file, obs_file)
    except Exception as e:
        print(f"ERROR: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()
