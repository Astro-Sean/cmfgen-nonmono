#!/usr/bin/env python3
"""Monitor and plot ARTIS SN_A series results.

For each completed epoch:
1. Plot the optical spectrum (spec.out) — 3000-10000 A
2. Plot the temperature profile before (input model.txt) and after (estimators) ARTIS
3. Generate a multi-epoch comparison plot

Usage:
    python3 scripts/plot_artis_series.py [--base-dir <dir>] [--monitor]
"""

import argparse
import os
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

C_LIGHT = 2.99792458e10        # cm/s
DAY_S = 86400.0
KB = 1.380649e-16              # erg/K
H_PLANCK = 6.62607015e-27      # erg*s
SIGMA_SB = 5.670374419e-5      # erg cm^-2 s^-1 K^-4
PI = np.pi
M_SUN = 1.989e33

ALL_ELEMENTS = ['H', 'He', 'C', 'N', 'O', 'Ne', 'Na', 'Mg', 'Al',
                'Si', 'S', 'Ar', 'Ca', 'Fe', 'Ni']
ELEMENT_Z = {'H': 1, 'He': 2, 'C': 6, 'N': 7, 'O': 8, 'Ne': 10,
             'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'S': 16,
             'Ar': 18, 'Ca': 20, 'Fe': 26, 'Ni': 28}

BASE_DIR = "/viper/ptmp1/sbrennan/SN2024pba/tests/artis_sn_a_series"


def read_spec_out(spec_path):
    """Read ARTIS spec.out file."""
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
    return freq, wl, ts_times, flux_ts


def read_model_txt(model_path):
    """Read ARTIS model.txt to get input structure.

    Returns: v_outer (km/s), logrho, q (erg/g), t_model (d), n_cells
    """
    with open(model_path) as f:
        lines = f.readlines()
    n_cells = int(lines[0].strip())
    t_model = float(lines[1].strip())
    # Skip header line (line 2)
    v_outer = []
    logrho = []
    q = []
    for i in range(3, 3 + n_cells):
        parts = lines[i].split()
        v_outer.append(float(parts[1]))
        logrho.append(float(parts[2]))
        q.append(float(parts[7]))
    return np.array(v_outer), np.array(logrho), np.array(q), t_model, n_cells


def read_abundances_txt(abund_path, n_cells):
    """Read abundances.txt to get composition per cell."""
    max_z = 30
    comp = np.zeros((n_cells, max_z))
    with open(abund_path) as f:
        for i, line in enumerate(f):
            if i >= n_cells:
                break
            parts = line.split()
            cellid = int(parts[0])
            for j in range(max_z):
                comp[cellid - 1, j] = float(parts[1 + j])
    return comp


def read_estimators_temp(est_path, n_cells):
    """Read ARTIS estimators file and extract T_e, T_R per cell for each timestep.

    Returns: dict with 'timesteps', 'Te' (n_ts, n_cells), 'TR' (n_ts, n_cells),
             'W' (n_ts, n_cells), 'nne' (n_ts, n_cells), 'tdays' (n_ts,)
    """
    ts_data = {}  # ts -> list of (mgi, TR, Te, W, nne, tdays)

    with open(est_path) as f:
        for line in f:
            if line.startswith('timestep '):
                parts = line.split()
                ts = int(parts[1])
                mgi = int(parts[3].split('=')[-1] if '=' in parts[3] else parts[3])
                # Parse: timestep N modelgridindex M titeration K TR ... Te ... W ... TJ ... grey_depth ... thick ... nne ... Ye ... tdays ...
                tr = None
                te = None
                w = None
                nne = None
                tdays = None
                for j, p in enumerate(parts):
                    if p == 'TR':
                        tr = float(parts[j + 1])
                    elif p == 'Te':
                        te = float(parts[j + 1])
                    elif p == 'W':
                        w = float(parts[j + 1])
                    elif p == 'nne':
                        nne = float(parts[j + 1])
                    elif p == 'tdays':
                        tdays = float(parts[j + 1])
                if ts not in ts_data:
                    ts_data[ts] = {}
                ts_data[ts][mgi] = (tr, te, w, nne, tdays)

    if not ts_data:
        return None

    timesteps = sorted(ts_data.keys())
    n_ts = len(timesteps)

    Te = np.full((n_ts, n_cells), np.nan)
    TR = np.full((n_ts, n_cells), np.nan)
    W = np.full((n_ts, n_cells), np.nan)
    NNe = np.full((n_ts, n_cells), np.nan)
    tdays = np.full(n_ts, np.nan)

    for i, ts in enumerate(timesteps):
        for mgi, (tr, te, w, nne, td) in ts_data[ts].items():
            if mgi < n_cells:
                Te[i, mgi] = te
                TR[i, mgi] = tr
                W[i, mgi] = w
                NNe[i, mgi] = nne
                tdays[i] = td

    return {
        'timesteps': timesteps,
        'Te': Te,
        'TR': TR,
        'W': W,
        'nne': NNe,
        'tdays': tdays,
    }


def get_input_temperature(model_path, abund_path):
    """Determine the input temperature profile from model.txt and abundances.

    The input script sets T_BB=7000K for metal-rich zones (H < 0.01)
    and T_OUTER=2000K for H/He zones, with a ramp transition.
    We reconstruct this from the abundances.
    """
    v_outer, logrho, q, t_model, n_cells = read_model_txt(model_path)
    comp = read_abundances_txt(abund_path, n_cells)

    # H is Z=1, index 0
    h_frac = comp[:, 0]
    is_metal = h_frac < 0.01

    T_BB = 7000.0
    T_OUTER = 2000.0

    t_input = np.full(n_cells, T_OUTER)
    t_input[is_metal] = T_BB

    # Ramp transition
    metal_indices = np.where(is_metal)[0]
    if len(metal_indices) > 0:
        last_metal = metal_indices[-1]
        ramp_width = min(20, n_cells - last_metal - 1)
        for i in range(last_metal + 1, last_metal + 1 + ramp_width):
            if i < n_cells:
                frac = (i - last_metal) / ramp_width
                t_input[i] = T_BB * (1 - frac) + T_OUTER * frac

    return v_outer, t_input, q, logrho, t_model, n_cells, comp


def plot_spectrum(work_dir, epoch_day):
    """Plot the ARTIS spectrum for a single epoch."""
    work_dir = Path(work_dir)
    spec_path = work_dir / "spec.out"
    if not spec_path.exists():
        return False

    freq, wl, ts_times, flux_ts = read_spec_out(spec_path)
    flux_total = flux_ts.sum(axis=1)
    flux_last = flux_ts[:, -1]

    # Clip to optical
    mask = (wl >= 3000) & (wl <= 10000)
    wl_plot = wl[mask]
    flux_plot = flux_total[mask]
    flux_last_plot = flux_last[mask]

    # Smoothed
    window = 15
    if len(flux_plot) > window:
        smooth = np.convolve(flux_plot, np.ones(window) / window, mode='same')
        smooth[:window // 2] = smooth[window // 2]
        smooth[-window // 2:] = smooth[-window // 2 - 1]
    else:
        smooth = flux_plot

    # Diagnostic lines
    lines_to_mark = [
        (3933, 'Ca II H&K'), (4102, 'H-delta'), (4340, 'H-gamma'),
        (4861, 'H-beta'), (5169, 'Fe II 5169'), (5876, 'He I 5876'),
        (6355, 'Si II 6355'), (6563, 'H-alpha'), (7065, 'He I 7065'),
        (7774, 'O I 7774'), (8498, 'Ca II IR'), (8662, 'Ca II IR'),
    ]

    fig, axes = plt.subplots(2, 1, figsize=(14, 10), gridspec_kw={'height_ratios': [2, 1]})

    # Top: Spectrum
    ax1 = axes[0]
    ax1.step(wl_plot, flux_plot, where='mid', color='steelblue',
             linewidth=0.5, alpha=0.5, label='Time-integrated (raw)')
    ax1.step(wl_plot, smooth, where='mid', color='crimson',
             linewidth=1.2, alpha=0.9, label='Smoothed (15-bin)')
    if len(ts_times) > 1:
        ax1.step(wl_plot, flux_last_plot, where='mid', color='darkgreen',
                 linewidth=0.8, alpha=0.6, label=f'Last timestep ({ts_times[-1]:.1f}d)')

    ymax = np.percentile(flux_plot[10:-10], 99) * 1.5 if len(flux_plot) > 20 else 1
    for w, name in lines_to_mark:
        ax1.axvline(w, color='gray', alpha=0.35, linestyle=':', linewidth=0.8)
        ax1.text(w, ymax * 0.97, name, rotation=90, fontsize=6,
                 ha='center', va='top', color='dimgray')
    ax1.set_xlim(3000, 10000)
    ax1.set_ylim(0, ymax)
    ax1.set_xlabel('Wavelength [Angstrom]')
    ax1.set_ylabel('Flux [arb. units]')
    ax1.set_title(f'ARTIS spectrum: SN_A at {epoch_day:.0f} d '
                  f'({len(freq)} freq bins, {len(ts_times)} timesteps)',
                  fontsize=11)
    ax1.legend(loc='upper right', fontsize=8)
    ax1.grid(True, alpha=0.2)

    # Bottom: Light curve
    ax2 = axes[1]
    lc_path = work_dir / "light_curve.out"
    if lc_path.exists():
        lc_t, lc_lum = [], []
        with open(lc_path) as f:
            for line in f:
                vals = line.split()
                if len(vals) >= 2:
                    lc_t.append(float(vals[0]))
                    lc_lum.append(float(vals[1]))
        ax2.plot(lc_t, lc_lum, 'b-', linewidth=1.5)
        ax2.axvline(epoch_day, color='gray', linestyle=':', alpha=0.5,
                    label=f't_model = {epoch_day:.0f}d')
        ax2.set_xlabel('Time [days]')
        ax2.set_ylabel('Bolometric luminosity [arb. units]')
        ax2.set_title('Bolometric light curve', fontsize=10)
        ax2.legend(fontsize=8)
        ax2.grid(True, alpha=0.2)
    else:
        ax2.text(0.5, 0.5, 'No light_curve.out', ha='center', va='center',
                 transform=ax2.transAxes)

    plt.tight_layout()
    plot_path = work_dir / f"spectrum_{epoch_day:.0f}d.png"
    fig.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved spectrum plot: {plot_path}")
    return True


def plot_temperature_profile(work_dir, epoch_day):
    """Plot temperature profile: input (model.txt) vs ARTIS output (estimators)."""
    work_dir = Path(work_dir)

    model_path = work_dir / "model.txt"
    abund_path = work_dir / "abundances.txt"
    est_path = work_dir / "estimators_0000.out"

    if not model_path.exists() or not abund_path.exists():
        print(f"  Missing model.txt or abundances.txt in {work_dir}")
        return False

    # Get input structure
    v_outer, t_input, q, logrho, t_model, n_cells, comp = get_input_temperature(
        model_path, abund_path)

    # Read ARTIS output temperatures
    artis_data = None
    if est_path.exists():
        artis_data = read_estimators_temp(est_path, n_cells)

    # Create figure
    if artis_data is not None:
        fig = plt.figure(figsize=(16, 14))
        gs = GridSpec(3, 2, figure=fig, hspace=0.35, wspace=0.25,
                      height_ratios=[1, 1, 1])
    else:
        fig = plt.figure(figsize=(16, 10))
        gs = GridSpec(2, 2, figure=fig, hspace=0.35, wspace=0.25)

    # Panel 1: Input temperature profile
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.step(v_outer, t_input, where='mid', color='crimson', linewidth=1.2,
             label='Input T_gas')
    ax1.axhline(7000, color='red', linestyle='--', alpha=0.4, label='T_BB=7000K')
    ax1.axhline(2000, color='blue', linestyle='--', alpha=0.4, label='T_outer=2000K')
    ax1.set_xlabel('Velocity [km/s]')
    ax1.set_ylabel('Temperature [K]')
    ax1.set_title(f'Input temperature profile (t_model={t_model:.0f} d)', fontsize=10)
    ax1.legend(fontsize=8)
    ax1.grid(True, alpha=0.2)

    # Panel 2: Input energy density (q)
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.step(v_outer, q, where='mid', color='darkorange', linewidth=1.2)
    ax2.set_yscale('log')
    ax2.set_xlabel('Velocity [km/s]')
    ax2.set_ylabel('q [erg/g]')
    ax2.set_title('Initial thermal energy density (q column)', fontsize=10)
    ax2.grid(True, alpha=0.2)
    q_max = q[q > 0].max() if (q > 0).any() else 1
    ax2.set_ylim(q_max * 1e-10, q_max * 10)

    # Panel 3: H abundance (to show metal/H boundary)
    ax3 = fig.add_subplot(gs[1, 0])
    h_frac = comp[:, 0]  # Z=1
    he_frac = comp[:, 1]  # Z=2
    ax3.step(v_outer, h_frac, where='mid', color='royalblue', linewidth=1, label='H')
    ax3.step(v_outer, he_frac, where='mid', color='skyblue', linewidth=1, label='He')
    si_idx = ELEMENT_Z['Si'] - 1
    fe_idx = ELEMENT_Z['Fe'] - 1
    ax3.step(v_outer, comp[:, si_idx], where='mid', color='crimson', linewidth=1, label='Si')
    ax3.step(v_outer, comp[:, fe_idx], where='mid', color='firebrick', linewidth=1, label='Fe')
    ax3.set_xlabel('Velocity [km/s]')
    ax3.set_ylabel('Mass fraction')
    ax3.set_title('Composition (showing H/metal boundary)', fontsize=10)
    ax3.legend(fontsize=8)
    ax3.set_ylim(-0.05, 1.0)
    ax3.grid(True, alpha=0.2)

    # Panel 4: Density
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.step(v_outer, 10**logrho, where='mid', color='steelblue', linewidth=1.2)
    ax4.set_yscale('log')
    ax4.set_xlabel('Velocity [km/s]')
    ax4.set_ylabel('Density [g/cm^3]')
    ax4.set_title('Density profile', fontsize=10)
    ax4.grid(True, alpha=0.2)

    if artis_data is not None:
        n_ts = len(artis_data['timesteps'])
        ts_list = artis_data['timesteps']
        Te = artis_data['Te']
        TR = artis_data['TR']
        W = artis_data['W']
        tdays = artis_data['tdays']

        # Panel 5: T_e and T_R evolution (selected timesteps)
        ax5 = fig.add_subplot(gs[2, 0])
        ax5.step(v_outer, t_input, where='mid', color='black', linewidth=1.5,
                 alpha=0.5, label='Input T_gas', linestyle='--')

        # Plot a few selected timesteps
        ts_to_plot = [0, min(5, n_ts - 1), min(n_ts // 2, n_ts - 1), n_ts - 1]
        ts_to_plot = sorted(set(ts_to_plot))
        colors_ts = plt.cm.viridis(np.linspace(0.1, 0.9, len(ts_to_plot)))

        for idx, ts_i in enumerate(ts_to_plot):
            if ts_i < n_ts:
                td = tdays[ts_i] if not np.isnan(tdays[ts_i]) else ts_i
                ax5.step(v_outer, Te[ts_i], where='mid', color=colors_ts[idx],
                         linewidth=1, alpha=0.8,
                         label=f'T_e ts{ts_i} ({td:.1f}d)')
                ax5.step(v_outer, TR[ts_i], where='mid', color=colors_ts[idx],
                         linewidth=1, alpha=0.4, linestyle=':')

        ax5.set_xlabel('Velocity [km/s]')
        ax5.set_ylabel('Temperature [K]')
        ax5.set_title('T_e (solid) & T_R (dotted) evolution through ARTIS', fontsize=10)
        ax5.legend(fontsize=7, loc='upper right')
        ax5.grid(True, alpha=0.2)
        ax5.set_ylim(0, 10000)

        # Panel 6: Dilution factor W
        ax6 = fig.add_subplot(gs[2, 1])
        for idx, ts_i in enumerate(ts_to_plot):
            if ts_i < n_ts:
                td = tdays[ts_i] if not np.isnan(tdays[ts_i]) else ts_i
                ax6.step(v_outer, W[ts_i], where='mid', color=colors_ts[idx],
                         linewidth=1, alpha=0.8,
                         label=f'W ts{ts_i} ({td:.1f}d)')
        ax6.set_xlabel('Velocity [km/s]')
        ax6.set_ylabel('Dilution factor W')
        ax6.set_title('Radiation dilution factor W', fontsize=10)
        ax6.legend(fontsize=7)
        ax6.set_ylim(0, 1.1)
        ax6.grid(True, alpha=0.2)

    fig.suptitle(f'ARTIS temperature profile: SN_A at {epoch_day:.0f} d\n'
                 f'Input vs ARTIS output ({n_cells} cells, t_model={t_model:.0f}d)',
                 fontsize=12, y=0.98)

    plot_path = work_dir / f"temperature_{epoch_day:.0f}d.png"
    fig.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved temperature plot: {plot_path}")
    return True


def plot_multi_epoch_comparison(base_dir, completed_epochs):
    """Create a multi-epoch comparison plot showing spectral evolution."""
    if len(completed_epochs) < 2:
        return False

    fig, axes = plt.subplots(2, 1, figsize=(16, 12), gridspec_kw={'height_ratios': [2, 1]})

    ax1 = axes[0]
    ax2 = axes[1]

    colors = plt.cm.plasma(np.linspace(0.1, 0.9, len(completed_epochs)))

    all_spectra = []
    for idx, (epoch_day, work_dir) in enumerate(completed_epochs):
        spec_path = Path(work_dir) / "spec.out"
        if not spec_path.exists():
            continue

        freq, wl, ts_times, flux_ts = read_spec_out(spec_path)
        flux_total = flux_ts.sum(axis=1)

        mask = (wl >= 3000) & (wl <= 10000)
        wl_plot = wl[mask]
        flux_plot = flux_total[mask]

        # Normalize by max for comparison
        f_max = flux_plot.max() if flux_plot.max() > 0 else 1
        flux_norm = flux_plot / f_max

        # Smoothed
        window = 15
        if len(flux_plot) > window:
            smooth = np.convolve(flux_norm, np.ones(window) / window, mode='same')
            smooth[:window // 2] = smooth[window // 2]
            smooth[-window // 2:] = smooth[-window // 2 - 1]
        else:
            smooth = flux_norm

        ax1.step(wl_plot, smooth, where='mid', color=colors[idx],
                 linewidth=1.2, alpha=0.85, label=f'{epoch_day:.0f} d')

        # Light curve
        lc_path = Path(work_dir) / "light_curve.out"
        if lc_path.exists():
            lc_t, lc_lum = [], []
            with open(lc_path) as f:
                for line in f:
                    vals = line.split()
                    if len(vals) >= 2:
                        lc_t.append(float(vals[0]))
                        lc_lum.append(float(vals[1]))
            ax2.plot(lc_t, lc_lum, '-', color=colors[idx], linewidth=1.2,
                     alpha=0.85, label=f'{epoch_day:.0f} d')

        all_spectra.append((epoch_day, wl_plot, flux_plot, f_max))

    # Line markers on spectrum
    lines_to_mark = [
        (3933, 'Ca II'), (4340, 'H-g'), (4861, 'H-b'),
        (5876, 'He I'), (6563, 'H-a'), (7774, 'O I'), (8662, 'Ca IR'),
    ]
    ymax = 1.5
    for w, name in lines_to_mark:
        ax1.axvline(w, color='gray', alpha=0.3, linestyle=':', linewidth=0.8)
        ax1.text(w, ymax * 0.97, name, rotation=90, fontsize=7,
                 ha='center', va='top', color='dimgray')

    ax1.set_xlim(3000, 10000)
    ax1.set_ylim(0, ymax)
    ax1.set_xlabel('Wavelength [Angstrom]')
    ax1.set_ylabel('Normalized flux')
    ax1.set_title(f'ARTIS spectral evolution: SN_A {completed_epochs[0][0]:.0f}–'
                  f'{completed_epochs[-1][0]:.0f} d (normalized, smoothed)',
                  fontsize=11)
    ax1.legend(fontsize=8, ncol=3, loc='upper right')
    ax1.grid(True, alpha=0.2)

    ax2.set_xlabel('Time [days]')
    ax2.set_ylabel('Bolometric luminosity [arb. units]')
    ax2.set_title('Bolometric light curves', fontsize=10)
    ax2.legend(fontsize=8, ncol=3)
    ax2.grid(True, alpha=0.2)

    plt.tight_layout()
    plot_path = Path(base_dir) / "multi_epoch_comparison.png"
    fig.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved multi-epoch comparison: {plot_path}")
    return True


def plot_multi_epoch_temperature(base_dir, completed_epochs):
    """Create a multi-epoch temperature comparison plot."""
    if len(completed_epochs) < 2:
        return False

    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    ax1 = axes[0]  # Input T
    ax2 = axes[1]  # ARTIS output T_e (last timestep)

    colors = plt.cm.plasma(np.linspace(0.1, 0.9, len(completed_epochs)))

    for idx, (epoch_day, work_dir) in enumerate(completed_epochs):
        model_path = Path(work_dir) / "model.txt"
        abund_path = Path(work_dir) / "abundances.txt"
        est_path = Path(work_dir) / "estimators_0000.out"

        if not model_path.exists() or not abund_path.exists():
            continue

        v_outer, t_input, q, logrho, t_model, n_cells, comp = get_input_temperature(
            model_path, abund_path)

        ax1.step(v_outer, t_input, where='mid', color=colors[idx],
                 linewidth=1.2, alpha=0.8, label=f'{epoch_day:.0f} d')

        if est_path.exists():
            artis_data = read_estimators_temp(est_path, n_cells)
            if artis_data is not None:
                n_ts = len(artis_data['timesteps'])
                Te_last = artis_data['Te'][-1]
                ax2.step(v_outer, Te_last, where='mid', color=colors[idx],
                         linewidth=1.2, alpha=0.8, label=f'{epoch_day:.0f} d (last ts)')

    ax1.set_xlabel('Velocity [km/s]')
    ax1.set_ylabel('Temperature [K]')
    ax1.set_title('Input temperature profiles (before ARTIS)', fontsize=11)
    ax1.legend(fontsize=8, ncol=2)
    ax1.set_ylim(0, 10000)
    ax1.grid(True, alpha=0.2)

    ax2.set_xlabel('Velocity [km/s]')
    ax2.set_ylabel('Temperature [K]')
    ax2.set_title('ARTIS output T_e (last timestep)', fontsize=11)
    ax2.legend(fontsize=8, ncol=2)
    ax2.set_ylim(0, 10000)
    ax2.grid(True, alpha=0.2)

    fig.suptitle(f'Temperature profile evolution: SN_A {completed_epochs[0][0]:.0f}–'
                 f'{completed_epochs[-1][0]:.0f} d',
                 fontsize=12, y=0.98)
    plt.tight_layout()
    plot_path = Path(base_dir) / "multi_epoch_temperature.png"
    fig.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved multi-epoch temperature comparison: {plot_path}")
    return True


def check_job_complete(work_dir):
    """Check if an ARTIS run has completed (spec.out exists and is non-trivial)."""
    spec_path = Path(work_dir) / "spec.out"
    if not spec_path.exists():
        return False
    # Check file is not empty
    if spec_path.stat().st_size < 1000:
        return False
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Monitor and plot ARTIS SN_A series results")
    parser.add_argument('--base-dir', type=str, default=BASE_DIR,
                        help=f'Base output directory (default: {BASE_DIR})')
    parser.add_argument('--monitor', action='store_true',
                        help='Keep running until all epochs complete')
    parser.add_argument('--poll-interval', type=int, default=60,
                        help='Poll interval in seconds for --monitor (default: 60)')
    args = parser.parse_args()

    base_dir = Path(args.base_dir)

    # Find all epoch directories
    epoch_dirs = []
    for d in sorted(base_dir.iterdir()):
        if d.is_dir() and d.name.endswith('d'):
            try:
                epoch_day = float(d.name[:-1])
                epoch_dirs.append((epoch_day, d))
            except ValueError:
                pass

    if not epoch_dirs:
        print(f"No epoch directories found in {base_dir}")
        return 1

    print(f"Found {len(epoch_dirs)} epoch directories: "
          f"{epoch_dirs[0][0]:.0f}d to {epoch_dirs[-1][0]:.0f}d")

    import time

    completed = []
    plotted = set()

    while True:
        # Check completion status
        newly_completed = []
        for epoch_day, work_dir in epoch_dirs:
            if epoch_day in plotted:
                continue
            if check_job_complete(work_dir):
                print(f"\n{'='*50}")
                print(f"Epoch {epoch_day:.0f}d: COMPLETE")
                # Plot spectrum
                ok1 = plot_spectrum(work_dir, epoch_day)
                # Plot temperature profile
                ok2 = plot_temperature_profile(work_dir, epoch_day)
                if ok1 and ok2:
                    plotted.add(epoch_day)
                    completed.append((epoch_day, str(work_dir)))
                    newly_completed.append((epoch_day, str(work_dir)))

        if newly_completed:
            # Update multi-epoch plots
            all_completed = [(ed, wd) for ed, wd in epoch_dirs if ed in plotted]
            if len(all_completed) >= 2:
                plot_multi_epoch_comparison(base_dir, all_completed)
                plot_multi_epoch_temperature(base_dir, all_completed)

        if not args.monitor:
            break

        if len(plotted) == len(epoch_dirs):
            print(f"\nAll {len(epoch_dirs)} epochs plotted!")
            break

        remaining = len(epoch_dirs) - len(plotted)
        print(f"  [{time.strftime('%H:%M:%S')}] {len(plotted)}/{len(epoch_dirs)} "
              f"complete, {remaining} pending. Waiting {args.poll_interval}s...")
        time.sleep(args.poll_interval)

    # Final summary
    print(f"\n{'='*60}")
    print(f"Summary: {len(plotted)}/{len(epoch_dirs)} epochs plotted")
    if len(plotted) < len(epoch_dirs):
        print("Pending epochs:")
        for epoch_day, work_dir in epoch_dirs:
            if epoch_day not in plotted:
                spec_exists = (Path(work_dir) / "spec.out").exists()
                print(f"  {epoch_day:.0f}d: spec.out={'yes' if spec_exists else 'no'}")

    # Always try multi-epoch plots with whatever we have
    all_completed = [(ed, wd) for ed, wd in epoch_dirs if ed in plotted]
    if len(all_completed) >= 2:
        plot_multi_epoch_comparison(base_dir, all_completed)
        plot_multi_epoch_temperature(base_dir, all_completed)

    return 0


if __name__ == '__main__':
    sys.exit(main())
