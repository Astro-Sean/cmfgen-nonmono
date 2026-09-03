#!/usr/bin/env python3
"""
Compare ARTIS spec.out files with observed SN2024pba spectra.

Reads spec.out from each ARTIS run directory, converts to flux at the
observed distance, and compares with the observed spectrum using
continuum shape and line absorption depths.

Usage:
    python compare_artis_to_obs.py [--distance MPC] [--series-dir DIR]
"""
import argparse
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

# Constants
C_LIGHT_ANG = 2.998e18   # Angstrom/s
C_LIGHT_CM = 2.998e10    # cm/s
PC_CM = 3.0857e18        # cm/pc
MPC_CM = 1e6 * PC_CM
DAY_S = 86400.0

# Observational data
OBS_FILES = {
    'phase_p25': '/ptmp/sbrennan/SN2024pba/observational_data/'
                 'SN2024pba_60875.145_PHASE25.1_Nordic-Optical-Telescope_ALFOSC.txt',
    'phase_m327': '/ptmp/sbrennan/SN2024pba/observational_data/'
                  'SN2024pba_60523.465_PHASE-326.5_Palomar-5.1m-Hale_DBSP.txt',
}

# Lines of interest (same as analyze_grid.py)
LINES = [
    {'name': 'Si II 6355', 'center': 6355.0, 'abs': (6280, 6420),
     'cont_b': (6450, 6550), 'cont_r': (6600, 6700)},
    {'name': 'O I 7774', 'center': 7774.0, 'abs': (7720, 7790),
     'cont_b': (7650, 7710), 'cont_r': (7820, 7880)},
    {'name': 'Ca II 8662', 'center': 8662.0, 'abs': (8620, 8690),
     'cont_b': (8720, 8780), 'cont_r': (8570, 8610)},
    {'name': 'Ca II 8542', 'center': 8542.0, 'abs': (8500, 8570),
     'cont_b': (8400, 8470), 'cont_r': (8610, 8660)},
    {'name': 'Halpha', 'center': 6563.0, 'abs': (6500, 6620),
     'cont_b': (6450, 6500), 'cont_r': (6700, 6750)},
    {'name': 'He I 5876', 'center': 5876.0, 'abs': (5820, 5920),
     'cont_b': (5750, 5800), 'cont_r': (5950, 6000)},
]


def read_observed_spectrum(path):
    """Read observed spectrum: wavelength (A), flux (erg/s/cm^2/A), error."""
    data = np.loadtxt(path)
    wave = data[:, 0]
    flux = data[:, 1]
    err = data[:, 2] if data.shape[1] > 2 else None
    return wave, flux, err


def read_artis_spec(spec_path, distance_mpc, timestep=-1):
    """Read ARTIS spec.out and convert to wavelength (A) vs flux (erg/s/cm^2/A).

    spec.out format: frequency (Hz) + luminosity per timestep (erg/s/Hz)
    Convert: F_nu = L_nu / (4*pi*d^2)
             F_lambda = F_nu * c / lambda^2
    """
    data = np.loadtxt(spec_path)
    nu = data[:, 0]          # Hz
    lum = data[:, 1:]        # erg/s/Hz, one column per timestep

    # Filter out zero-frequency rows
    valid = nu > 0
    nu = nu[valid]
    lum = lum[valid]

    # Use the specified timestep (default: last)
    if timestep < 0:
        timestep = lum.shape[1] - 1
    lum_ts = lum[:, timestep]

    # Convert frequency to wavelength
    wave = C_LIGHT_ANG / nu  # Angstrom

    # Convert luminosity to flux
    d_cm = distance_mpc * MPC_CM
    f_nu = lum_ts / (4.0 * np.pi * d_cm**2)  # erg/s/cm^2/Hz

    # Convert F_nu to F_lambda
    f_lambda = f_nu * C_LIGHT_ANG / wave**2   # erg/s/cm^2/A

    # Sort by wavelength
    idx = np.argsort(wave)
    wave = wave[idx]
    f_lambda = f_lambda[idx]

    return wave, f_lambda


def measure_line_depths(wave, flux):
    """Measure absorption line depths relative to local continuum."""
    results = {}
    for line in LINES:
        center = line['center']
        abs_lo, abs_hi = line['abs']
        cb_lo, cb_hi = line['cont_b']
        cr_lo, cr_hi = line['cont_r']

        # Continuum from blue and red windows
        blue_mask = (wave >= cb_lo) & (wave <= cb_hi)
        red_mask = (wave >= cr_lo) & (wave <= cr_hi)

        if blue_mask.sum() < 2 or red_mask.sum() < 2:
            results[line['name']] = None
            continue

        # Linear continuum interpolation
        blue_w = wave[blue_mask].mean()
        blue_f = np.median(flux[blue_mask])
        red_w = wave[red_mask].mean()
        red_f = np.median(flux[red_mask])

        # Absorption region
        abs_mask = (wave >= abs_lo) & (wave <= abs_hi)
        if abs_mask.sum() < 2:
            results[line['name']] = None
            continue

        abs_w = wave[abs_mask]
        abs_f = flux[abs_mask]

        # Continuum at absorption wavelengths
        if red_w != blue_w:
            cont = blue_f + (red_f - blue_f) * (abs_w - blue_w) / (red_w - blue_w)
        else:
            cont = np.full_like(abs_w, (blue_f + red_f) / 2)

        # Line depth = 1 - min(F_abs / F_cont)
        ratio = abs_f / np.maximum(cont, 1e-30)
        depth = 1.0 - np.min(ratio)
        depth = max(0.0, min(1.0, depth))

        results[line['name']] = depth

    return results


def compute_chi2(model_wave, model_flux, obs_wave, obs_flux,
                 wave_range=(4000, 9000)):
    """Compute normalized chi-squared between model and observed spectrum."""
    # Interpolate model onto observed wavelength grid
    mask = (obs_wave >= wave_range[0]) & (obs_wave <= wave_range[1])
    if mask.sum() < 10:
        return np.inf

    obs_w = obs_wave[mask]
    obs_f = obs_flux[mask]

    # Interpolate model
    model_interp = np.interp(obs_w, model_wave, model_flux, left=0, right=0)

    # Normalize both to their median in the comparison window
    med_obs = np.median(obs_f[obs_f > 0]) if np.any(obs_f > 0) else 1
    med_mod = np.median(model_interp[model_interp > 0]) if np.any(model_interp > 0) else 1

    obs_norm = obs_f / med_obs
    mod_norm = model_interp / med_mod

    # Chi-squared (shape comparison, no error weighting)
    valid = (obs_norm > 0) & (mod_norm > 0)
    if valid.sum() < 10:
        return np.inf

    chi2 = np.mean((obs_norm[valid] - mod_norm[valid])**2)
    return chi2


def main():
    parser = argparse.ArgumentParser(
        description='Compare ARTIS spectra with observed SN2024pba')
    parser.add_argument('--series-dir', type=str,
                        default='/viper/ptmp1/sbrennan/SN2024pba/tests/artis_sn_a_series',
                        help='Directory containing ARTIS run subdirectories')
    parser.add_argument('--distance', type=float, default=83.5,
                        help='Distance to SN2024pba in Mpc (default: 83.5)')
    parser.add_argument('--timestep', type=int, default=-1,
                        help='ARTIS timestep to use (default: last)')
    parser.add_argument('--obs', type=str, default='phase_p25',
                        choices=['phase_p25', 'phase_m327', 'both'],
                        help='Which observed spectrum to compare with')
    parser.add_argument('--output', type=str,
                        default='/ptmp/sbrennan/SN2024pba/tests/artis_comparison.png',
                        help='Output plot path')
    args = parser.parse_args()

    # Read observed spectra
    obs_data = {}
    for key, path in OBS_FILES.items():
        if os.path.isfile(path):
            w, f, e = read_observed_spectrum(path)
            obs_data[key] = {'wave': w, 'flux': f, 'err': e}
            print(f"Read observed spectrum: {key} ({len(w)} points, "
                  f"{w.min():.0f}-{w.max():.0f} A)")
        else:
            print(f"WARNING: Observed spectrum not found: {path}")

    if not obs_data:
        print("ERROR: No observed spectra found")
        sys.exit(1)

    # Find all ARTIS run directories with spec.out
    series_dir = Path(args.series_dir)
    run_dirs = []
    for d in sorted(series_dir.iterdir()):
        spec = d / 'spec.out'
        if spec.is_file():
            # Check if it's from the current run (not stale)
            run_dirs.append(d)

    if not run_dirs:
        print(f"ERROR: No spec.out files found in {series_dir}")
        sys.exit(1)

    print(f"\nFound {len(run_dirs)} ARTIS runs with spec.out")

    # Compare each model with observations
    results = []
    for d in run_dirs:
        epoch = d.name  # e.g., "300d"
        spec_path = d / 'spec.out'

        try:
            mw, mf = read_artis_spec(str(spec_path), args.distance,
                                     timestep=args.timestep)
        except Exception as e:
            print(f"  ERROR reading {spec_path}: {e}")
            continue

        print(f"\n  {epoch}: {len(mw)} points, "
              f"{mw.min():.0f}-{mw.max():.0f} A")

        # Measure line depths
        model_lines = measure_line_depths(mw, mf)

        # Compare with each observed spectrum
        for obs_key, obs in obs_data.items():
            if args.obs != 'both' and args.obs != obs_key:
                continue

            obs_lines = measure_line_depths(obs['wave'], obs['flux'])
            chi2 = compute_chi2(mw, mf, obs['wave'], obs['flux'])

            # Line depth differences
            line_diffs = {}
            for name in model_lines:
                if model_lines[name] is not None and obs_lines.get(name) is not None:
                    line_diffs[name] = abs(model_lines[name] - obs_lines[name])

            mean_line_diff = np.mean(list(line_diffs.values())) if line_diffs else np.inf

            results.append({
                'epoch': epoch,
                'obs': obs_key,
                'chi2': chi2,
                'mean_line_diff': mean_line_diff,
                'model_lines': model_lines,
                'obs_lines': obs_lines,
                'model_wave': mw,
                'model_flux': mf,
            })

            print(f"    vs {obs_key}: chi2={chi2:.4f}, "
                  f"mean_line_diff={mean_line_diff:.4f}")

    if not results:
        print("\nERROR: No valid comparisons made")
        sys.exit(1)

    # Rank by combined metric (lower is better)
    for r in results:
        r['score'] = r['chi2'] + 0.5 * r['mean_line_diff']

    results.sort(key=lambda x: x['score'])

    print("\n" + "=" * 70)
    print("RANKING (best to worst):")
    print(f"{'Rank':<5} {'Epoch':<8} {'Obs':<12} {'Chi2':<10} {'LineDiff':<10} {'Score':<10}")
    print("-" * 70)
    for i, r in enumerate(results):
        print(f"{i+1:<5} {r['epoch']:<8} {r['obs']:<12} "
              f"{r['chi2']:<10.4f} {r['mean_line_diff']:<10.4f} "
              f"{r['score']:<10.4f}")

    # Plot best few models against observations
    best_epochs = list(dict.fromkeys(r['epoch'] for r in results[:4]))

    fig, axes = plt.subplots(len(best_epochs) + 1, 1,
                             figsize=(14, 4 * (len(best_epochs) + 1)),
                             sharex=True)
    if len(best_epochs) == 0:
        axes = [axes]

    # Top panel: observed spectra
    ax = axes[0]
    for key, obs in obs_data.items():
        if args.obs != 'both' and args.obs != key:
            continue
        ax.plot(obs['wave'], obs['flux'], label=f'Obs: {key}', alpha=0.8)
    ax.set_ylabel('Flux (erg/s/cm$^2$/\\AA)')
    ax.set_title('Observed SN2024pba spectra')
    ax.legend()
    ax.set_xlim(3500, 9500)

    # Lower panels: best models vs best-matching obs
    best_obs = results[0]['obs']
    obs = obs_data[best_obs]

    for i, epoch in enumerate(best_epochs):
        ax = axes[i + 1]
        # Find the result for this epoch
        r = next(r for r in results if r['epoch'] == epoch)

        # Normalize both to continuum near 7000 A
        mask = (r['model_wave'] > 6900) & (r['model_wave'] < 7100)
        mod_norm = np.median(r['model_flux'][mask]) if mask.any() else 1
        obs_mask = (obs['wave'] > 6900) & (obs['wave'] < 7100)
        obs_norm = np.median(obs['flux'][obs_mask]) if obs_mask.any() else 1

        ax.plot(obs['wave'], obs['flux'] / obs_norm, 'k-', alpha=0.5,
                label=f'Obs ({best_obs})')
        ax.plot(r['model_wave'], r['model_flux'] / mod_norm, 'r-',
                label=f'ARTIS {epoch} (chi2={r["chi2"]:.3f})')
        ax.set_ylabel('Normalized flux')
        ax.set_title(f'ARTIS {epoch} vs observed ({best_obs})')
        ax.legend(fontsize=8)
        ax.set_xlim(3500, 9500)

    axes[-1].set_xlabel('Wavelength (\\AA)')
    plt.tight_layout()
    plt.savefig(args.output, dpi=150)
    print(f"\nSaved comparison plot to {args.output}")

    # Print line depth table
    print("\n" + "=" * 70)
    print("LINE DEPTHS (model vs observed):")
    print(f"{'Epoch':<8} {'Obs':<12}", end='')
    for line in LINES:
        print(f" {line['name']:<14}", end='')
    print()
    print("-" * 70)

    for r in results[:6]:
        print(f"{r['epoch']:<8} {r['obs']:<12}", end='')
        for line in LINES:
            ml = r['model_lines'].get(line['name'])
            ol = r['obs_lines'].get(line['name'])
            if ml is not None and ol is not None:
                print(f" {ml:.2f}/{ol:.2f}    ", end='')
            else:
                print(f" {'--':<14}", end='')
        print()


if __name__ == '__main__':
    main()
