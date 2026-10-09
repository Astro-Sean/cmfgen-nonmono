#!/usr/bin/env python3
"""Raw observations vs models reddened by MW + host extinction.

The observed photometry is plotted with NO extinction correction
(apparent mags minus the distance modulus).  Each model is shifted to
its best-fit time offset under MW-only extinction, then progressively
reddened into the observed frame:

    Mabs_reddened(band) = Mabs_model(band) + E(B-V)_tot * coeff_G23(band)

One column per model; four curves per panel for MW only (E=0.15) and
MW+host with E(B-V)_host = 0.10 / 0.20 / 0.30, all on the Gordon+2023
Rv=3.1 law.  This isolates what extra extinction would do to the
lightcurve and the g-r / o-g / o-r colours at fixed model physics.

Usage:
    python scripts/plot_extinction_test.py [MODEL ...]
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import plot_dessart_grid_vs_obs as base

EBV_MW = 0.15
# G23 Rv=3.1 band-averaged coefficients (same as the main plot)
COEFF = {'g': 3.776, 'r': 2.695, 'o': 2.350}
OBS2MOD = {'g': 'g', 'r': 'R', 'o': 'o'}

HOST_LEVELS = [0.0, 0.10, 0.20, 0.30]
LEVEL_COLORS = ['#1f77b4', '#2ca02c', '#ff7f0e', '#d62728']


def redden(df, band, ebv_tot):
    """Model Mabs with extinction added (into the observed frame)."""
    return df[OBS2MOD[band]] + ebv_tot * COEFF[band]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('models', nargs='*',
                    default=['DX14_V8000_sa330',
                             'DX14_V8000_w3e17_sa350'],
                    help='model directory names to plot')
    ap.add_argument('--grid-dir', default=base.GRID_DIR)
    ap.add_argument('--fit-end', type=float, default=100.0)
    args = ap.parse_args()

    # Raw observations: apparent mags minus DM, zero extinction
    base.set_ebv(0.0)
    obs = base.load_obs()

    # MW-only obs for the shift/score evaluation
    base.set_ebv(EBV_MW)
    obs_mw = base.load_obs()

    nmod = len(args.models)
    fig, axes = plt.subplots(4, nmod, figsize=(5.5 * nmod, 13),
                             sharex=True, sharey='row',
                             gridspec_kw={'hspace': 0.12, 'wspace': 0.08})
    if nmod == 1:
        axes = axes[:, None]

    for col, name in enumerate(args.models):
        mdir = os.path.join(args.grid_dir, name)
        tt = base.find_tt(mdir)
        if tt is None:
            print(f'{name}: no .tt output, skipping')
            continue
        df = base.load_tt(tt)
        tpk = df.loc[df['R'].idxmin(), 'time']
        s = base.best_phase_shift(df, obs_mw, tpk, args.fit_end)
        score, dt = s[0], s[1]
        x = df['time'] - tpk + dt
        print(f'{name}: dt={dt:+.2f} d, MW-only score={score:.1f}')

        for eh, colr in zip(HOST_LEVELS, LEVEL_COLORS):
            et = EBV_MW + eh
            lab = (f'E(B-V)={et:.2f} MW only' if eh == 0
                   else f'E(B-V)={et:.2f} (+{eh:.2f} host)')
            kw = dict(color=colr, lw=1.5, alpha=0.9)
            axes[0, col].plot(x, redden(df, 'r', et), label=lab, **kw)
            axes[1, col].plot(
                x, redden(df, 'g', et) - redden(df, 'r', et),
                label=lab, **kw)
            axes[2, col].plot(
                x, redden(df, 'o', et) - redden(df, 'g', et),
                label=lab, **kw)
            axes[3, col].plot(
                x, redden(df, 'o', et) - redden(df, 'r', et),
                label=lab, **kw)

        axes[0, col].set_title(f'{name}\nMW-only score {score:.1f}',
                               fontsize=10)
        axes[0, col].legend(fontsize=8)

    for col in range(nmod):
        axes[0, col].scatter(obs['r']['phase'], obs['r']['Mabs'], s=16,
                             c='k', zorder=10, label='obs r raw')
        axes[0, col].scatter(obs['o']['phase'], obs['o']['Mabs'], s=12,
                             facecolors='none', edgecolors='gray',
                             zorder=9, label='obs o raw')
        for row, (b1, b2) in enumerate((('g', 'r'), ('o', 'g'),
                                        ('o', 'r')), start=1):
            pairs = base.obs_color_pairs(obs, b1, b2, max_sep=1.5)
            axes[row, col].errorbar(
                pairs['phase'], pairs['col'], yerr=pairs['err'],
                fmt='ko', ms=3, lw=0.8, capsize=2, zorder=10,
                label='obs raw')

    axes[0, 0].set_ylabel('M abs (m - DM, raw)')
    axes[1, 0].set_ylabel('g - r')
    axes[2, 0].set_ylabel('o - g')
    axes[3, 0].set_ylabel('o - r')
    # sharey='row': inverting one axis per row flips both columns
    for row in range(4):
        axes[row, 0].invert_yaxis()  # bright / blue at top
    for ax in axes.flat:
        ax.grid(True, alpha=0.3)
        ax.set_xlim(-25, 135)
    for ax in axes[-1]:
        ax.set_xlabel('phase (d from r-band peak)')

    out = os.path.join(args.grid_dir, 'extinction_test.png')
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Saved: {out}')


if __name__ == '__main__':
    main()
