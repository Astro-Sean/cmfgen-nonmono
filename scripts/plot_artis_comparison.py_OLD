#!/usr/bin/env python3
"""Compare ARTIS model spectra with observed SN2024pba spectra.

ARTIS spec.out is in F_nu [erg/s/cm2/Hz] at 1 Mpc.
We scale to the observed distance and compare with the two observed spectra.

Usage:
    python plot_artis_comparison.py [--distance-mpc D] [--timestep TS]
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse
import os
from pathlib import Path

LSUN = 3.828e33
C = 2.998e10  # cm/s
MPC_CM = 3.086e24


def read_artis_spec(spec_path, timestep=0):
    """Read ARTIS spec.out and return wavelength [A], F_lambda [erg/s/cm2/A] at 1 Mpc."""
    data = np.loadtxt(spec_path)
    freqs = data[1:, 0]
    F_nu = data[1:, 1 + timestep]  # erg/s/cm2/Hz at 1 Mpc

    wave_A = C / freqs * 1e8
    F_lambda = F_nu * C / wave_A**2 * 1e8  # erg/s/cm2/A at 1 Mpc

    # Sort by wavelength
    idx = np.argsort(wave_A)
    return wave_A[idx], F_lambda[idx]


def read_observed_spec(path):
    """Read observed spectrum: wavelength [A], flux [erg/s/cm2/A]."""
    data = np.loadtxt(path)
    return data[:, 0], data[:, 1]


def main():
    parser = argparse.ArgumentParser(description="Compare ARTIS spectra with observations")
    parser.add_argument('--series-dir', type=str,
                        default='/viper/ptmp1/sbrennan/SN2024pba/tests/artis_sn_a_series',
                        help='Directory with ARTIS model series')
    parser.add_argument('--obs-dir', type=str,
                        default='/ptmp/sbrennan/SN2024pba/observational_data',
                        help='Directory with observed spectra')
    parser.add_argument('--distance-mpc', type=float, default=100.0,
                        help='Distance to SN in Mpc (for scaling model flux)')
    parser.add_argument('--timestep', type=int, default=0,
                        help='ARTIS timestep to plot (0 = first = peak)')
    parser.add_argument('--output', type=str,
                        default='/ptmp/sbrennan/SN2024pba/tests/artis_sn_a_series/comparison.png',
                        help='Output plot path')
    args = parser.parse_args()

    # Read observed spectra
    obs1_path = os.path.join(args.obs_dir,
                             'SN2024pba_60523.465_PHASE-326.5_Palomar-5.1m-Hale_DBSP.txt')
    obs2_path = os.path.join(args.obs_dir,
                             'SN2024pba_60875.145_PHASE25.1_Nordic-Optical-Telescope_ALFOSC.txt')

    obs1_wave, obs1_flux = read_observed_spec(obs1_path)
    obs2_wave, obs2_flux = read_observed_spec(obs2_path)

    # Scale factor: model is at 1 Mpc, scale to observed distance
    scale = 1.0 / (args.distance_mpc**2)

    # Read all available ARTIS spectra
    days = [300, 310, 320, 330, 340, 350, 360, 370, 380, 390, 400]
    artis_specs = {}
    for day in days:
        spec_path = os.path.join(args.series_dir, f'{day}d', 'spec.out')
        if os.path.isfile(spec_path):
            wave, flux_1mpc = read_artis_spec(spec_path, args.timestep)
            flux_scaled = flux_1mpc * scale
            artis_specs[day] = (wave, flux_scaled)
            print(f"Read {day}d: {len(wave)} bins, peak F_lambda={flux_scaled[flux_scaled>0].max():.3e}")

    if not artis_specs:
        print("No ARTIS spectra found!")
        return

    # Create plot
    fig, axes = plt.subplots(2, 1, figsize=(14, 10), sharex=True)

    # Top panel: all ARTIS models + observations
    ax = axes[0]
    colors = plt.cm.viridis(np.linspace(0, 1, len(artis_specs)))
    for (day, (wave, flux)), color in zip(sorted(artis_specs.items()), colors):
        mask = (wave > 3000) & (wave < 11000)
        ax.plot(wave[mask], flux[mask], alpha=0.7, color=color, label=f'ARTIS {day}d')

    ax.plot(obs1_wave, obs1_flux, 'k--', alpha=0.8, label='Obs (phase -326.5d)')
    ax.plot(obs2_wave, obs2_flux, 'r--', alpha=0.8, label='Obs (phase +25.1d)')

    ax.set_ylabel(r'$F_\lambda$ [erg/s/cm$^2$/\AA]')
    ax.set_title(f'ARTIS SN_A models vs SN2024pba observations (d={args.distance_mpc} Mpc, ts={args.timestep})')
    ax.legend(fontsize=8, ncol=2)
    ax.set_yscale('log')
    ax.set_ylim(1e-19, 1e-13)
    ax.set_xlim(3000, 11000)

    # Bottom panel: ratio of best-matching model to observations
    ax2 = axes[1]
    # Find best matching model (by chi-squared in optical)
    best_day = None
    best_chi2 = np.inf
    for day, (wave, flux) in artis_specs.items():
        # Interpolate model to observed wavelengths
        mask = (wave > 3000) & (wave < 10000) & (flux > 0)
        if not mask.any():
            continue
        model_interp = np.interp(obs2_wave, wave[mask], flux[mask], left=0, right=0)
        obs_mask = obs2_flux > 0
        if obs_mask.any():
            # Normalize by matching total flux
            norm = np.trapezoid(obs2_flux[obs_mask], obs2_wave[obs_mask]) / \
                   max(np.trapezoid(model_interp[obs_mask], obs2_wave[obs_mask]), 1e-30)
            model_normed = model_interp * norm
            chi2 = np.sum(((model_normed[obs_mask] - obs2_flux[obs_mask]) /
                          (obs2_flux[obs_mask] + 1e-20))**2)
            if chi2 < best_chi2:
                best_chi2 = chi2
                best_day = day

    if best_day is not None:
        wave, flux = artis_specs[best_day]
        mask = (wave > 3000) & (wave < 10000) & (flux > 0)
        model_interp = np.interp(obs2_wave, wave[mask], flux[mask], left=0, right=0)
        obs_mask = obs2_flux > 0
        norm = np.trapezoid(obs2_flux[obs_mask], obs2_wave[obs_mask]) / \
               max(np.trapezoid(model_interp[obs_mask], obs2_wave[obs_mask]), 1e-30)
        model_normed = model_interp * norm

        ax2.plot(obs2_wave, model_normed, 'b-', alpha=0.7, label=f'ARTIS {best_day}d (normed)')
        ax2.plot(obs2_wave, obs2_flux, 'r--', alpha=0.7, label='Obs (phase +25.1d)')
        ax2.set_ylabel(r'$F_\lambda$ [erg/s/cm$^2$/\AA]')
        ax2.set_xlabel(r'Wavelength [\AA]')
        ax2.set_title(f'Best match: ARTIS {best_day}d (chi2={best_chi2:.1f}, norm={norm:.2e})')
        ax2.legend()
        ax2.set_yscale('log')
        ax2.set_ylim(1e-18, 1e-14)
        ax2.set_xlim(3000, 11000)
    else:
        ax2.text(0.5, 0.5, 'No valid model found', transform=ax2.transAxes, ha='center')

    plt.tight_layout()
    plt.savefig(args.output, dpi=150)
    print(f"\nSaved comparison plot to {args.output}")
    if best_day is not None:
        print(f"Best matching model: {best_day}d (chi2={best_chi2:.1f})")


if __name__ == '__main__':
    main()
