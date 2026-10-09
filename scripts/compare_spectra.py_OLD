#!/usr/bin/env python3
"""Compare SuperLite model spectra to observed SN2024pba spectrum.

Reads SuperLite output.flx_luminos + output.flx_grid, converts to
observed flux at SN2024pba distance (z=0.0186, D_L=77.5 Mpc),
redshifts the wavelengths, and plots against the observed spectrum.

Usage:
    python compare_spectra.py <obs_file> <model_dir1> [model_dir2 ...]
"""
import numpy as np
import sys
from pathlib import Path

# SN2024pba parameters
Z = 0.0186
try:
    from astropy.cosmology import FlatLambdaCDM
    cosmo = FlatLambdaCDM(H0=73, Om0=0.3)
    D_L_Mpc = cosmo.luminosity_distance(Z).value
except ImportError:
    D_L_Mpc = 73.0 * Z / 73.0 * (1 + Z) * 3e5 / 73.0  # rough
D_L_cm = D_L_Mpc * 3.086e24
FLUX_TO_LUM = 4 * np.pi * D_L_cm**2

def read_superlite_spectrum(model_dir):
    """Read SuperLite spectrum from output.flx_luminos + output.flx_grid.

    output.flx_grid: wavelength bin edges in cm (first line is header '# ndim ...')
    output.flx_luminos: luminosity per bin in erg/s (same format, may have header)
    """
    model_dir = Path(model_dir)
    flx_file = model_dir / 'output.flx_luminos'
    grid_file = model_dir / 'output.flx_grid'

    if not flx_file.exists():
        return None

    # Read wavelength grid (wavelength bin edges in cm)
    with open(grid_file) as f:
        lines = f.readlines()
    # Skip comment/header lines
    wl_data = []
    for line in lines:
        if line.strip().startswith('#'):
            continue
        wl_data.extend([float(x) for x in line.split()])
    wl_edges = np.array(wl_data)

    # Read luminosities (erg/s per bin)
    with open(flx_file) as f:
        lines = f.readlines()
    lum_data = []
    for line in lines:
        if line.strip().startswith('#'):
            continue
        lum_data.extend([float(x) for x in line.split()])
    lum_per_bin = np.array(lum_data)

    # Wavelength bin centers (cm -> Angstrom)
    # Grid has N+1 edges for N bins; luminosities has N values
    n_bins = len(lum_per_bin)
    wl_edges = wl_edges[:n_bins + 1]  # trim to match
    wl_centers = 0.5 * (wl_edges[:-1] + wl_edges[1:]) * 1e8  # Angstrom
    dwl = (wl_edges[1:] - wl_edges[:-1]) * 1e8  # Angstrom

    # Luminosity per bin -> luminosity per Angstrom
    lum_per_A = lum_per_bin / np.maximum(dwl, 1e-30)

    # Sort by wavelength
    order = np.argsort(wl_centers)
    wl = wl_centers[order]
    lum = lum_per_A[order]

    return wl, lum

def read_observed_spectrum(obs_file):
    """Read observed spectrum (wl_A, flux_erg_s_cm2_A, err)."""
    data = np.loadtxt(obs_file)
    wl = data[:, 0]
    flux = data[:, 1]
    # Filter NaN
    mask = ~np.isnan(flux) & ~np.isnan(wl)
    return wl[mask], flux[mask]

def main():
    if len(sys.argv) < 3:
        print("Usage: python compare_spectra.py <obs_file> <model_dir1> [model_dir2 ...]")
        sys.exit(1)

    obs_file = sys.argv[1]
    model_dirs = sys.argv[2:]

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    # Read observed spectrum
    obs_wl, obs_flux = read_observed_spectrum(obs_file)
    obs_name = Path(obs_file).stem

    # Rest-frame wavelengths
    obs_wl_rest = obs_wl / (1 + Z)

    fig, axes = plt.subplots(2, 1, figsize=(16, 12), gridspec_kw={'height_ratios': [2, 1]})

    # Plot observed spectrum (redshifted to rest frame)
    axes[0].plot(obs_wl_rest, obs_flux, 'k-', alpha=0.5, linewidth=0.5,
                 label=f'Observed ({obs_name})')

    colors = plt.cm.tab10(np.linspace(0, 1, len(model_dirs)))

    for i, model_dir in enumerate(model_dirs):
        result = read_superlite_spectrum(model_dir)
        if result is None:
            print(f"  {model_dir}: no output.flx_luminos found, skipping")
            continue

        wl_model, lum_model = result

        # Convert model luminosity to observed flux
        # F = L / (4 * pi * D_L^2)
        flux_model = lum_model / FLUX_TO_LUM

        # Redshift model wavelengths to observed frame
        wl_model_obs = wl_model * (1 + Z)

        # Plot in rest frame
        wl_model_rest = wl_model
        model_name = Path(model_dir).parent.parent.name + '/' + Path(model_dir).name

        axes[0].plot(wl_model_rest, flux_model, color=colors[i], linewidth=1.0,
                     label=model_name, alpha=0.8)

        # Ratio plot (model / observed) in overlapping wavelength range
        from scipy.interpolate import interp1d
        mask = (wl_model_rest > obs_wl_rest.min()) & (wl_model_rest < obs_wl_rest.max())
        if mask.sum() > 10:
            ip = interp1d(wl_model_rest[mask], flux_model[mask], kind='linear',
                          bounds_error=False, fill_value=np.nan)
            model_at_obs = ip(obs_wl_rest)
            # Smooth observed for ratio
            from scipy.ndimage import median_filter
            obs_smooth = median_filter(obs_flux, size=15)
            ratio = model_at_obs / np.maximum(np.abs(obs_smooth), 1e-30)
            valid = np.abs(ratio) < 10  # filter extremes
            axes[1].plot(obs_wl_rest[valid], ratio[valid], color=colors[i],
                        linewidth=0.5, alpha=0.7, label=model_name)

        print(f"  {model_name}: wl=[{wl_model.min():.0f}, {wl_model.max():.0f}] A, "
              f"max flux={flux_model.max():.2e}")

    axes[0].set_xlabel('Rest-frame Wavelength (Å)')
    axes[0].set_ylabel('Flux (erg/s/cm²/Å)')
    axes[0].set_title(f'SN2024pba: Model vs Observed (z={Z}, D_L={D_L_Mpc:.1f} Mpc)')
    axes[0].legend(fontsize=8, loc='upper right')
    axes[0].set_xlim(3000, 10000)
    axes[0].set_ylim(-1e-15, 2e-15)
    axes[0].grid(True, alpha=0.3)

    axes[1].set_xlabel('Rest-frame Wavelength (Å)')
    axes[1].set_ylabel('Model / Observed')
    axes[1].axhline(1.0, color='k', linestyle='--', alpha=0.5)
    axes[1].set_xlim(3000, 10000)
    axes[1].set_ylim(0, 3)
    axes[1].legend(fontsize=7)
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    out_file = '/viper/ptmp1/sbrennan/SN2024pba/results/spectrum_comparison.png'
    Path(out_file).parent.mkdir(exist_ok=True)
    plt.savefig(out_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\nSaved: {out_file}")

if __name__ == '__main__':
    main()
