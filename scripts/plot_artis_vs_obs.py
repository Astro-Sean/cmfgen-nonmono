#!/usr/bin/env python3
"""Compare ARTIS SN_A series spectra with SN2024pba observations.

The observational data is at phase ~25d while ARTIS epochs are 300-400d,
so we normalize each spectrum for shape comparison.

Usage:
    python3 scripts/plot_artis_vs_obs.py [--base-dir <dir>] [--obs-file <path>]
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

C_LIGHT = 2.99792458e10  # cm/s

OBS_FILE = "/ptmp/sbrennan/SN2024pba/observational_data/SN2024pba_60875.145_PHASE25.1_Nordic-Optical-Telescope_ALFOSC.txt"
BASE_DIR = "/viper/ptmp1/sbrennan/SN2024pba/tests/artis_sn_a_series"


def read_observation(obs_path):
    """Read observational spectrum: wavelength (A), flux, error."""
    wl, flux, err = [], [], []
    with open(obs_path) as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 2:
                wl.append(float(parts[0]))
                flux.append(float(parts[1]))
                if len(parts) >= 3:
                    err.append(float(parts[2]))
                else:
                    err.append(0.0)
    return np.array(wl), np.array(flux), np.array(err)


def read_artis_spec(spec_path):
    """Read ARTIS spec.out: header=timestep times, then freq + flux per ts."""
    with open(spec_path) as f:
        lines = f.readlines()
    header = np.array([float(x) for x in lines[0].split()])
    ts_times = header[1:]
    freq_list = []
    flux_list = []
    for line in lines[1:]:
        vals = [float(x) for x in line.split()]
        freq_list.append(vals[0])
        flux_list.append(vals[1:])
    freq = np.array(freq_list)
    flux_ts = np.array(flux_list)
    wl = C_LIGHT / freq * 1e8  # Angstrom
    return wl, ts_times, flux_ts


def normalize(wl, flux, wl_min=3500, wl_max=9500, smooth_window=15):
    """Normalize spectrum to unit max in the optical range, with optional smoothing."""
    mask = (wl >= wl_min) & (wl <= wl_max)
    f = flux[mask].copy()
    if smooth_window > 0 and len(f) > smooth_window:
        kernel = np.ones(smooth_window) / smooth_window
        f = np.convolve(f, kernel, mode='same')
        f[:smooth_window // 2] = f[smooth_window // 2]
        f[-smooth_window // 2:] = f[-smooth_window // 2 - 1]
    f_max = f.max()
    if f_max > 0:
        f_norm = f / f_max
    else:
        f_norm = f
    wl_plot = wl[mask]
    return wl_plot, f_norm


def main():
    parser = argparse.ArgumentParser(
        description="Compare ARTIS spectra with SN2024pba observations")
    parser.add_argument('--base-dir', type=str, default=BASE_DIR)
    parser.add_argument('--obs-file', type=str, default=OBS_FILE)
    args = parser.parse_args()

    base_dir = Path(args.base_dir)
    obs_path = Path(args.obs_file)

    # Read observation
    obs_wl, obs_flux, obs_err = read_observation(obs_path)
    obs_wl_n, obs_flux_n = normalize(obs_wl, obs_flux, smooth_window=5)
    print(f"Observation: {len(obs_wl)} points, {obs_wl[0]:.0f}-{obs_wl[-1]:.0f} A, "
          f"phase 25.1d")

    # Find completed ARTIS epochs
    epoch_dirs = []
    for d in sorted(base_dir.iterdir()):
        if d.is_dir() and d.name.endswith('d'):
            spec = d / "spec.out"
            if spec.exists() and spec.stat().st_size > 1000:
                try:
                    epoch_day = float(d.name[:-1])
                    epoch_dirs.append((epoch_day, d))
                except ValueError:
                    pass

    if not epoch_dirs:
        print("No completed ARTIS epochs found!")
        return 1

    print(f"Found {len(epoch_dirs)} completed ARTIS epochs: "
          f"{epoch_dirs[0][0]:.0f}d to {epoch_dirs[-1][0]:.0f}d")

    # --- Plot 1: All ARTIS epochs + observation (normalized) ---
    fig, ax = plt.subplots(figsize=(16, 8))

    # Observation
    ax.step(obs_wl_n, obs_flux_n, where='mid', color='black',
            linewidth=1.5, alpha=0.8, label='SN2024pba obs (phase 25.1d)')

    # ARTIS epochs
    colors = plt.cm.plasma(np.linspace(0.1, 0.9, len(epoch_dirs)))
    for idx, (epoch_day, work_dir) in enumerate(epoch_dirs):
        spec_path = work_dir / "spec.out"
        wl, ts_times, flux_ts = read_artis_spec(spec_path)
        flux_total = flux_ts.sum(axis=1)
        wl_n, f_n = normalize(wl, flux_total)
        ax.step(wl_n, f_n, where='mid', color=colors[idx],
                linewidth=1.2, alpha=0.75, label=f'ARTIS {epoch_day:.0f}d')

    # Line markers
    lines_to_mark = [
        (3933, 'Ca II H&K'), (4102, 'H-delta'), (4340, 'H-gamma'),
        (4861, 'H-beta'), (5169, 'Fe II'), (5876, 'He I'),
        (6355, 'Si II'), (6563, 'H-alpha'), (7065, 'He I'),
        (7774, 'O I'), (8498, 'Ca II IR'), (8662, 'Ca II IR'),
    ]
    ymax = 1.6
    for w, name in lines_to_mark:
        ax.axvline(w, color='gray', alpha=0.3, linestyle=':', linewidth=0.8)
        ax.text(w, ymax * 0.97, name, rotation=90, fontsize=7,
                ha='center', va='top', color='dimgray')

    ax.set_xlim(3400, 10000)
    ax.set_ylim(0, ymax)
    ax.set_xlabel('Wavelength [Angstrom]', fontsize=12)
    ax.set_ylabel('Normalized flux (peak=1)', fontsize=12)
    ax.set_title('ARTIS SN_A series vs SN2024pba observation\n'
                 '(normalized — observation at phase 25.1d, ARTIS at 300-400d)',
                 fontsize=13)
    ax.legend(fontsize=9, ncol=2, loc='upper right')
    ax.grid(True, alpha=0.2)

    plt.tight_layout()
    plot_path = base_dir / "artis_vs_obs_comparison.png"
    fig.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {plot_path}")

    # --- Plot 2: Each ARTIS epoch in its own panel vs observation ---
    n_epochs = len(epoch_dirs)
    n_cols = min(3, n_epochs)
    n_rows = (n_epochs + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(7 * n_cols, 5 * n_rows),
                             squeeze=False)

    for idx, (epoch_day, work_dir) in enumerate(epoch_dirs):
        row, col = idx // n_cols, idx % n_cols
        ax = axes[row][col]

        spec_path = work_dir / "spec.out"
        wl, ts_times, flux_ts = read_artis_spec(spec_path)
        flux_total = flux_ts.sum(axis=1)
        wl_n, f_n = normalize(wl, flux_total)

        ax.step(obs_wl_n, obs_flux_n, where='mid', color='black',
                linewidth=1.2, alpha=0.7, label='Obs (25.1d)')
        ax.step(wl_n, f_n, where='mid', color=colors[idx],
                linewidth=1.0, alpha=0.8, label=f'ARTIS {epoch_day:.0f}d')

        for w, name in lines_to_mark:
            ax.axvline(w, color='gray', alpha=0.2, linestyle=':', linewidth=0.6)

        ax.set_xlim(3400, 10000)
        ax.set_ylim(0, 1.5)
        ax.set_xlabel('Wavelength [A]', fontsize=10)
        ax.set_ylabel('Normalized flux', fontsize=10)
        ax.set_title(f'ARTIS {epoch_day:.0f}d vs SN2024pba (phase 25.1d)',
                     fontsize=11)
        ax.legend(fontsize=8, loc='upper right')
        ax.grid(True, alpha=0.15)

    # Hide unused subplots
    for idx in range(n_epochs, n_rows * n_cols):
        row, col = idx // n_cols, idx % n_cols
        axes[row][col].set_visible(False)

    plt.tight_layout()
    plot_path2 = base_dir / "artis_vs_obs_panels.png"
    fig.savefig(plot_path2, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {plot_path2}")

    # --- Plot 3: Last timestep vs time-integrated for each epoch ---
    fig, axes2 = plt.subplots(n_rows, n_cols, figsize=(7 * n_cols, 5 * n_rows),
                              squeeze=False)

    for idx, (epoch_day, work_dir) in enumerate(epoch_dirs):
        row, col = idx // n_cols, idx % n_cols
        ax = axes2[row][col]

        spec_path = work_dir / "spec.out"
        wl, ts_times, flux_ts = read_artis_spec(spec_path)
        flux_total = flux_ts.sum(axis=1)
        flux_last = flux_ts[:, -1]

        wl_n1, f_n1 = normalize(wl, flux_total)
        wl_n2, f_n2 = normalize(wl, flux_last)

        ax.step(obs_wl_n, obs_flux_n, where='mid', color='black',
                linewidth=1.2, alpha=0.6, label='Obs (25.1d)')
        ax.step(wl_n1, f_n1, where='mid', color='steelblue',
                linewidth=0.8, alpha=0.5, label='Time-integrated')
        ax.step(wl_n2, f_n2, where='mid', color='crimson',
                linewidth=1.0, alpha=0.8, label=f'Last ts ({ts_times[-1]:.1f}d)')

        for w, name in lines_to_mark:
            ax.axvline(w, color='gray', alpha=0.2, linestyle=':', linewidth=0.6)

        ax.set_xlim(3400, 10000)
        ax.set_ylim(0, 1.5)
        ax.set_xlabel('Wavelength [A]', fontsize=10)
        ax.set_ylabel('Normalized flux', fontsize=10)
        ax.set_title(f'ARTIS {epoch_day:.0f}d: integrated vs last timestep',
                     fontsize=11)
        ax.legend(fontsize=8, loc='upper right')
        ax.grid(True, alpha=0.15)

    for idx in range(n_epochs, n_rows * n_cols):
        row, col = idx // n_cols, idx % n_cols
        axes2[row][col].set_visible(False)

    plt.tight_layout()
    plot_path3 = base_dir / "artis_integrated_vs_lastts.png"
    fig.savefig(plot_path3, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {plot_path3}")

    print(f"\nDone! {len(epoch_dirs)} ARTIS epochs compared with observation.")
    return 0


if __name__ == '__main__':
    sys.exit(main())
