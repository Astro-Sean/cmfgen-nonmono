#!/usr/bin/env python3
"""
Compare CMFGEN observed spectrum (obs_fin_10) with an observed SN spectrum.

The CMFGEN obs_fin_10 file contains flux in Janskys, normalized to a
source at 1 kpc (see obs_frame_sub_v5.f line ~851:
    OBS_FLUX = OBS_FLUX * 6.59934 * R(1)^2   !Janskys (1kpc)
where 6.59934 = 1e23 * 2*pi * 1e20 / (3.0856e21^2)).

To compare with an observed spectrum at redshift z:

1. Scale CMFGEN flux from 1 kpc to the luminosity distance:
   F_obs = F_1kpc * (1 kpc / D_L[kpc])^2

2. Convert CMFGEN F_nu (Jy) to F_lambda (erg/s/cm^2/A):
   F_lambda = F_nu[Jy] * 1e-23 * c[cm/s] / lambda[cm]^2 * 1e-8
   (the 1e-8 converts per-cm to per-Angstrom)

3. De-redshift the observed spectrum:
   - Wavelength: lambda_rest = lambda_obs / (1+z)
   - Flux: F_lambda_rest = F_lambda_obs * (1+z)
     (accounts for bandwidth compression by redshift)

Usage:
    python compare_cmfgen_to_observation.py \
        --cmfgen /path/to/obs_fin_10 \
        --obs /path/to/observed_spectrum.txt \
        --redshift 0.0186 \
        --out comparison.png
"""

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# Cosmology
# ---------------------------------------------------------------------------

def luminosity_distance(z, H0=70.0, Omega_M=0.3, Omega_L=0.7):
    """
    Compute luminosity distance in Mpc for a flat LCDM cosmology.

    Uses the 2nd-order expansion in z (accurate to <0.1% for z < 0.1):
        D_L = (c/H0) * [z + 0.5*(1-q0)*z^2] * (1+z)
    where q0 = 0.5*(Omega_M - 2*Omega_L).
    """
    c_km_s = 2.99792458e5  # km/s
    q0 = 0.5 * (Omega_M - 2.0 * Omega_L)
    D_L = (c_km_s / H0) * (z + 0.5 * (1.0 - q0) * z**2) * (1.0 + z)
    return D_L  # Mpc


# ---------------------------------------------------------------------------
# File readers
# ---------------------------------------------------------------------------

def read_obs_fin_10(filename):
    """
    Read CMFGEN obs_fin_10 (OBSFRAME) ASCII format.

    Returns (wavelength_A, flux_Jy_at_1kpc) where:
      - wavelength_A is in Angstroms
      - flux_Jy_at_1kpc is in Janskys, normalized to 1 kpc
    """
    with open(filename, 'r') as f:
        lines = f.readlines()

    # Find frequency header
    header_idx = None
    for i, line in enumerate(lines):
        if 'Continuum Frequencies' in line:
            header_idx = i
            break
    if header_idx is None:
        raise ValueError("Could not find 'Continuum Frequencies' header")

    nfreq = int(lines[header_idx].split('(')[1].split(')')[0].strip())

    # Find flux header
    flux_header_idx = None
    for i, line in enumerate(lines):
        if 'Observed intensity' in line:
            flux_header_idx = i
            break
    if flux_header_idx is None:
        raise ValueError("Could not find 'Observed intensity' header")

    def _parse_fortran_float(s):
        s = re.sub(r'([0-9])([+-][0-9])', r'\1E\2', s)
        return float(s)

    # Read frequencies (in units of 10^15 Hz)
    freq = []
    for i in range(header_idx + 1, flux_header_idx):
        line = lines[i].strip()
        if line:
            freq.extend([_parse_fortran_float(v) for v in line.split()])
    freq = np.array(freq[:nfreq])

    # Read fluxes (in Jy at 1 kpc)
    flux_jy = []
    for i in range(flux_header_idx + 1, len(lines)):
        line = lines[i].strip()
        if line:
            flux_jy.extend([_parse_fortran_float(v) for v in line.split()])
    flux_jy = np.array(flux_jy[:nfreq])

    # Convert frequency to wavelength
    c_cm_s = 2.99792458e10  # cm/s
    freq_hz = freq * 1e15  # 10^15 Hz -> Hz
    wavelength_cm = c_cm_s / freq_hz
    wavelength_A = wavelength_cm * 1e8  # cm -> Angstroms

    # Sort by wavelength
    idx = np.argsort(wavelength_A)
    return wavelength_A[idx], flux_jy[idx]


def read_observed_spectrum(filename):
    """
    Read observed spectrum in whitespace-delimited ASCII format.

    Expected columns: wavelength(A)  F_lambda(erg/s/cm^2/A)  uncertainty

    Handles NaN values and negative fluxes (sets to NaN).
    """
    wl, fl, er = [], [], []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                w = float(parts[0])
                f_val = float(parts[1])
                e = float(parts[2]) if len(parts) > 2 else 0.0
            except ValueError:
                continue
            wl.append(w)
            fl.append(f_val)
            er.append(e)

    wl = np.array(wl)
    fl = np.array(fl)
    er = np.array(er)

    # Mark NaN and negative fluxes as invalid
    fl[np.isnan(fl)] = 0.0
    fl[fl < 0] = 0.0

    # Sort by wavelength
    idx = np.argsort(wl)
    return wl[idx], fl[idx], er[idx]


# ---------------------------------------------------------------------------
# Flux conversion and scaling
# ---------------------------------------------------------------------------

def cmfgen_jy_to_flambda(wavelength_A, flux_jy_1kpc, D_L_Mpc):
    """
    Convert CMFGEN flux from Jy@1kpc to F_lambda at the observer.

    Parameters:
        wavelength_A : array of wavelengths in Angstroms (rest frame)
        flux_jy_1kpc : array of flux in Janskys, normalized to 1 kpc
        D_L_Mpc : luminosity distance in Mpc

    Returns:
        F_lambda in erg/s/cm^2/A at the observer distance
    """
    c_cm_s = 2.99792458e10  # cm/s
    D_L_kpc = D_L_Mpc * 1e3  # Mpc -> kpc

    # Scale from 1 kpc to D_L
    distance_factor = (1.0 / D_L_kpc) ** 2
    flux_jy_obs = flux_jy_1kpc * distance_factor

    # Convert Jy -> erg/s/cm^2/Hz
    flux_cgs_nu = flux_jy_obs * 1e-23  # Jy -> erg/s/cm^2/Hz

    # Convert F_nu -> F_lambda (per Angstrom)
    wavelength_cm = wavelength_A * 1e-8  # A -> cm
    flux_cgs_lambda = flux_cgs_nu * (c_cm_s / wavelength_cm ** 2)  # erg/s/cm^2/cm
    flux_cgs_lambda *= 1e-8  # per-cm -> per-Angstrom

    return flux_cgs_lambda


def deredshift_observed(wavelength_obs, flux_obs, z):
    """
    Convert observed-frame spectrum to rest-frame.

    Parameters:
        wavelength_obs : observed wavelengths in Angstroms
        flux_obs : F_lambda in erg/s/cm^2/A (observed frame)
        z : redshift

    Returns:
        (wavelength_rest, flux_rest) in rest frame
    """
    wavelength_rest = wavelength_obs / (1.0 + z)
    flux_rest = flux_obs * (1.0 + z)  # bandwidth compression correction
    return wavelength_rest, flux_rest


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

# Major SN lines (rest-frame wavelengths)
SN_LINES = {
    3889: 'He I', 4102: 'Hδ', 4340: 'Hγ', 4861: 'Hβ',
    4924: 'Fe II', 5018: 'Fe II', 5160: 'Fe II',
    5876: 'He I', 6300: '[O I]', 6364: '[O I]',
    6563: 'Hα', 7065: 'He I', 7774: 'O I',
    8446: 'O I', 8542: 'Ca II', 8662: 'Ca II',
    8498: 'Ca II',
}


def wl_to_velocity(wavelength, rest_wl):
    """Convert wavelength to velocity (km/s) relative to rest wavelength."""
    c_kms = 2.99792458e5
    return c_kms * (wavelength - rest_wl) / rest_wl


def smooth_spectrum(flux, window=15):
    """
    Light moving-average smoothing of a 1-D flux array.

    Preserves broad P-Cygni profiles (200-300 A) while suppressing
    pixel-to-pixel noise.  Edges are handled with a shrinking window
    so no NaNs are introduced.  window should be odd; if even, +1.

    Parameters:
        flux : 1-D array
        window : smoothing window size in pixels (default 15)

    Returns:
        smoothed 1-D array (same length, no NaNs)
    """
    if window < 3:
        return flux.copy()
    if window % 2 == 0:
        window += 1
    n = len(flux)
    if n < window:
        return flux.copy()
    half = window // 2
    # Use cumulative-sum trick for O(N) moving average
    f = np.asarray(flux, dtype=float)
    csum = np.concatenate(([0.0], np.cumsum(f)))
    smoothed = np.empty(n)
    for i in range(n):
        lo = max(0, i - half)
        hi = min(n, i + half + 1)
        smoothed[i] = (csum[hi] - csum[lo]) / (hi - lo)
    return smoothed


def plot_comparison(wl_cmf, fl_cmf, wl_obs, fl_obs, model_label, obs_label,
                    z, D_L_Mpc, output_file, wl_range=None):
    """
    Plot CMFGEN vs observed spectrum comparison.

    Panel 1: absolute F_lambda, full optical range
    Panel 2: continuum-normalized, full optical range
    Panel 3: Hα region in velocity space (km/s relative to 6563 A)
    Panel 4: Hβ region in velocity space (km/s relative to 4861 A)
    """
    if wl_range is None:
        wl_range = (3500, 9000)

    fig, axes = plt.subplots(4, 1, figsize=(16, 18),
                             gridspec_kw={'height_ratios': [3, 2, 2, 2]})

    mask_c = (wl_cmf >= wl_range[0]) & (wl_cmf <= wl_range[1])
    mask_o = (wl_obs >= wl_range[0]) & (wl_obs <= wl_range[1])

    # Compute continuum levels for normalization
    fl_cmf_masked = fl_cmf[mask_c].copy()
    fl_obs_masked = fl_obs[mask_o].copy()
    cmf_cont = np.median(fl_cmf_masked[fl_cmf_masked > 0]) if np.any(fl_cmf_masked > 0) else 1.0
    obs_cont = np.median(fl_obs_masked[fl_obs_masked > 0]) if np.any(fl_obs_masked > 0) else 1.0

    # ---- Panel 1: absolute flux comparison ----
    ax = axes[0]
    ax.plot(wl_cmf[mask_c], fl_cmf[mask_c], linewidth=0.4, color='blue',
            alpha=0.8, label=f'{model_label} (CMFGEN, scaled to D_L={D_L_Mpc:.1f} Mpc)')
    fl_obs_smooth_p1 = smooth_spectrum(fl_obs[mask_o], window=15)
    ax.plot(wl_obs[mask_o], fl_obs[mask_o], linewidth=0.5, color='red',
            alpha=0.25, label=f'{obs_label} (raw)')
    ax.plot(wl_obs[mask_o], fl_obs_smooth_p1, linewidth=1.2, color='red',
            alpha=0.9, label=f'{obs_label} (smoothed)')

    ax.set_yscale('linear')
    ax.set_xlim(wl_range)
    ax.set_title(f'CMFGEN vs Observed Spectrum — {model_label} vs {obs_label}  '
                 f'(z={z}, D_L={D_L_Mpc:.1f} Mpc, mu={5*np.log10(D_L_Mpc*1e6)-5:.2f})',
                 fontsize=13)
    ax.set_ylabel(r'$F_\lambda$ (erg s$^{-1}$ cm$^{-2}$ Å$^{-1}$)', fontsize=11)
    ax.legend(fontsize=10, loc='upper right')
    ax.grid(True, alpha=0.3)
    for wl, name in SN_LINES.items():
        if wl_range[0] <= wl <= wl_range[1]:
            ax.axvline(wl, color='gray', alpha=0.2, linestyle='--')
            ax.text(wl, ax.get_ylim()[1] * 0.95, name, fontsize=6,
                    ha='center', color='gray', rotation=90)

    # ---- Panel 2: continuum-normalized ----
    ax2 = axes[1]
    fl_cmf_norm = fl_cmf_masked / cmf_cont
    fl_obs_norm = fl_obs_masked / obs_cont
    fl_obs_norm_smooth = smooth_spectrum(fl_obs_norm, window=15)

    ax2.plot(wl_cmf[mask_c], fl_cmf_norm, linewidth=0.4, color='blue',
             alpha=0.8, label=f'{model_label} (norm. to continuum)')
    ax2.plot(wl_obs[mask_o], fl_obs_norm, linewidth=0.5, color='red',
             alpha=0.25, label=f'{obs_label} (raw, norm.)')
    ax2.plot(wl_obs[mask_o], fl_obs_norm_smooth, linewidth=1.2, color='red',
             alpha=0.9, label=f'{obs_label} (smoothed, norm.)')

    ax2.set_yscale('linear')
    ax2.set_xlim(wl_range)
    ax2.set_ylim(-0.2, 2.5)
    ax2.set_title('Continuum-normalized comparison (each divided by own median)',
                  fontsize=12)
    ax2.set_ylabel(r'$F_\lambda / F_{\rm cont}$', fontsize=11)
    ax2.axhline(1.0, color='black', alpha=0.3, linestyle=':')
    ax2.legend(fontsize=10, loc='upper right')
    ax2.grid(True, alpha=0.3)
    for wl, name in SN_LINES.items():
        if wl_range[0] <= wl <= wl_range[1]:
            ax2.axvline(wl, color='gray', alpha=0.2, linestyle='--')
            ax2.text(wl, 2.4, name, fontsize=6,
                     ha='center', color='gray', rotation=90)

    # ---- Panel 3: Hα in velocity space ----
    ax3 = axes[2]
    ha_rest = 6563.0
    ha_wl_range = (6100, 7100)  # ~21,000 km/s range
    mask_c_ha = (wl_cmf >= ha_wl_range[0]) & (wl_cmf <= ha_wl_range[1])
    mask_o_ha = (wl_obs >= ha_wl_range[0]) & (wl_obs <= ha_wl_range[1])

    if np.any(mask_c_ha):
        v_cmf_ha = wl_to_velocity(wl_cmf[mask_c_ha], ha_rest)
        fl_cmf_ha = fl_cmf[mask_c_ha] / cmf_cont
        ax3.plot(v_cmf_ha, fl_cmf_ha, linewidth=0.5, color='blue',
                 alpha=0.8, label=f'{model_label}')

    if np.any(mask_o_ha):
        v_obs_ha = wl_to_velocity(wl_obs[mask_o_ha], ha_rest)
        fl_obs_ha = fl_obs[mask_o_ha] / obs_cont
        fl_obs_ha_smooth = smooth_spectrum(fl_obs_ha, window=15)
        ax3.plot(v_obs_ha, fl_obs_ha, linewidth=0.6, color='red',
                 alpha=0.25, label=f'{obs_label} (raw)')
        ax3.plot(v_obs_ha, fl_obs_ha_smooth, linewidth=1.4, color='red',
                 alpha=0.9, label=f'{obs_label} (smoothed)')

    ax3.set_yscale('linear')
    v_lim = 20000  # km/s
    ax3.set_xlim(-v_lim, v_lim)
    ax3.set_title('Hα 6563 Å — velocity space (continuum-normalized)', fontsize=12)
    ax3.set_xlabel('Velocity (km/s)', fontsize=11)
    ax3.set_ylabel(r'$F_\lambda / F_{\rm cont}$', fontsize=11)
    ax3.axhline(1.0, color='black', alpha=0.3, linestyle=':')
    ax3.axvline(0, color='red', alpha=0.3, linestyle='--', label='Hα rest')
    # Mark typical velocities
    for v_mark, label in [(850, 'CSM ~850'), (-850, None), (6400, 'SN ~6400'), (-6400, None)]:
        ax3.axvline(v_mark, color='green', alpha=0.2, linestyle=':')
        if label:
            ax3.text(v_mark, ax3.get_ylim()[1] * 0.85, label, fontsize=7,
                     ha='center', color='green', rotation=90)
    ax3.legend(fontsize=9, loc='upper right')
    ax3.grid(True, alpha=0.3)

    # ---- Panel 4: Hβ in velocity space ----
    ax4 = axes[3]
    hb_rest = 4861.0
    hb_wl_range = (4600, 5100)  # ~30,000 km/s range
    mask_c_hb = (wl_cmf >= hb_wl_range[0]) & (wl_cmf <= hb_wl_range[1])
    mask_o_hb = (wl_obs >= hb_wl_range[0]) & (wl_obs <= hb_wl_range[1])

    if np.any(mask_c_hb):
        v_cmf_hb = wl_to_velocity(wl_cmf[mask_c_hb], hb_rest)
        fl_cmf_hb = fl_cmf[mask_c_hb] / cmf_cont
        ax4.plot(v_cmf_hb, fl_cmf_hb, linewidth=0.5, color='blue',
                 alpha=0.8, label=f'{model_label}')

    if np.any(mask_o_hb):
        v_obs_hb = wl_to_velocity(wl_obs[mask_o_hb], hb_rest)
        fl_obs_hb = fl_obs[mask_o_hb] / obs_cont
        fl_obs_hb_smooth = smooth_spectrum(fl_obs_hb, window=15)
        ax4.plot(v_obs_hb, fl_obs_hb, linewidth=0.6, color='red',
                 alpha=0.25, label=f'{obs_label} (raw)')
        ax4.plot(v_obs_hb, fl_obs_hb_smooth, linewidth=1.4, color='red',
                 alpha=0.9, label=f'{obs_label} (smoothed)')

    ax4.set_yscale('linear')
    ax4.set_xlim(-20000, 20000)
    ax4.set_title('Hβ 4861 Å — velocity space (continuum-normalized)', fontsize=12)
    ax4.set_xlabel('Velocity (km/s)', fontsize=11)
    ax4.set_ylabel(r'$F_\lambda / F_{\rm cont}$', fontsize=11)
    ax4.axhline(1.0, color='black', alpha=0.3, linestyle=':')
    ax4.axvline(0, color='red', alpha=0.3, linestyle='--', label='Hβ rest')
    ax4.legend(fontsize=9, loc='upper right')
    ax4.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Comparison plot saved to {output_file}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Compare CMFGEN obs_fin_10 spectrum with an observed SN spectrum.'
    )
    parser.add_argument('--cmfgen', required=True,
                        help='Path to CMFGEN obs_fin_10 file')
    parser.add_argument('--obs', required=True,
                        help='Path to observed spectrum (ASCII: wl[A] F_lambda err)')
    parser.add_argument('--redshift', type=float, default=0.0186,
                        help='Redshift of the SN (default: 0.0186 for SN2024pba)')
    parser.add_argument('--H0', type=float, default=70.0,
                        help='Hubble constant in km/s/Mpc (default: 70)')
    parser.add_argument('--model-label', default='CMFGEN model',
                        help='Label for the CMFGEN model in the plot')
    parser.add_argument('--obs-label', default='Observed',
                        help='Label for the observed spectrum in the plot')
    parser.add_argument('--out', '-o', default='comparison.png',
                        help='Output PNG file (default: comparison.png)')
    parser.add_argument('--wl-min', type=float, default=3500,
                        help='Minimum wavelength to plot (default: 3500 A)')
    parser.add_argument('--wl-max', type=float, default=9000,
                        help='Maximum wavelength to plot (default: 9000 A)')
    parser.add_argument('--no-deredshift', action='store_true',
                        help='Do not de-redshift the observed spectrum '
                             '(use if already in rest frame)')

    args = parser.parse_args()

    # --- Read CMFGEN spectrum ---
    cmfgen_path = Path(args.cmfgen)
    if not cmfgen_path.exists():
        print(f"ERROR: CMFGEN file not found: {args.cmfgen}", file=sys.stderr)
        sys.exit(1)

    print(f"Reading CMFGEN spectrum from {args.cmfgen}...")
    wl_cmf, fl_jy_1kpc = read_obs_fin_10(args.cmfgen)
    print(f"  {len(wl_cmf)} points, wavelength range: "
          f"{wl_cmf.min():.1f} - {wl_cmf.max():.1f} A")

    # --- Read observed spectrum ---
    obs_path = Path(args.obs)
    if not obs_path.exists():
        print(f"ERROR: Observed spectrum not found: {args.obs}", file=sys.stderr)
        sys.exit(1)

    print(f"Reading observed spectrum from {args.obs}...")
    wl_obs_raw, fl_obs_raw, er_obs = read_observed_spectrum(args.obs)
    print(f"  {len(wl_obs_raw)} points, wavelength range: "
          f"{wl_obs_raw.min():.1f} - {wl_obs_raw.max():.1f} A")

    # --- De-redshift observed spectrum ---
    if args.no_deredshift:
        wl_obs = wl_obs_raw
        fl_obs = fl_obs_raw
        print("  (skipping de-redshift -- spectrum assumed to be in rest frame)")
    else:
        wl_obs, fl_obs = deredshift_observed(wl_obs_raw, fl_obs_raw, args.redshift)
        print(f"  De-redshifted by z={args.redshift}: "
              f"wavelength range {wl_obs.min():.1f} - {wl_obs.max():.1f} A")

    # --- Compute luminosity distance ---
    D_L_Mpc = luminosity_distance(args.redshift, H0=args.H0)
    D_L_kpc = D_L_Mpc * 1e3
    mu = 5 * np.log10(D_L_Mpc * 1e6) - 5
    print(f"\nLuminosity distance: D_L = {D_L_Mpc:.2f} Mpc ({D_L_kpc:.1f} kpc)")
    print(f"Distance modulus: mu = {mu:.3f} mag")
    print(f"Flux scaling factor (1 kpc -> D_L): {(1.0/D_L_kpc)**2:.6e}")

    # --- Convert CMFGEN flux to observed F_lambda ---
    fl_cmf = cmfgen_jy_to_flambda(wl_cmf, fl_jy_1kpc, D_L_Mpc)

    # --- Print flux ranges for sanity check ---
    mask_opt = (wl_cmf >= args.wl_min) & (wl_cmf <= args.wl_max)
    if np.any(mask_opt):
        fl_cmf_opt = fl_cmf[mask_opt]
        fl_cmf_pos = fl_cmf_opt[fl_cmf_opt > 0]
        if len(fl_cmf_pos) > 0:
            print(f"\nCMFGEN F_lambda (optical, at D_L={D_L_Mpc:.1f} Mpc):")
            print(f"  Range: {fl_cmf_pos.min():.3e} - {fl_cmf_pos.max():.3e} "
                  f"erg/s/cm^2/A")

    mask_obs_opt = (wl_obs >= args.wl_min) & (wl_obs <= args.wl_max)
    if np.any(mask_obs_opt):
        fl_obs_opt = fl_obs[mask_obs_opt]
        fl_obs_pos = fl_obs_opt[fl_obs_opt > 0]
        if len(fl_obs_pos) > 0:
            print(f"\nObserved F_lambda (optical, rest frame):")
            print(f"  Range: {fl_obs_pos.min():.3e} - {fl_obs_pos.max():.3e} "
                  f"erg/s/cm^2/A")

    # --- Plot comparison ---
    print(f"\nGenerating comparison plot...")
    plot_comparison(
        wl_cmf, fl_cmf,
        wl_obs, fl_obs,
        model_label=args.model_label,
        obs_label=args.obs_label,
        z=args.redshift,
        D_L_Mpc=D_L_Mpc,
        output_file=args.out,
        wl_range=(args.wl_min, args.wl_max),
    )

    print("Done.")


if __name__ == '__main__':
    main()
