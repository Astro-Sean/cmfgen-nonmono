#!/usr/bin/env python3
"""Bolometric comparison of the two-component CSM X-grid models.

Reads the STELLA .swd snapshots (via extract_lightcurve_from_swd) and
plots L_bol, photospheric T_bb and R_bb for X, X_M4 and the three
two-component CSM variants.  A model still running (partial .swd) is
drawn with a thin line and its endpoint marked.

Usage:
    python scripts/plot_x_csm2_bolometric.py [--grid-dir SN_B_X_grid]
"""

import argparse
import glob
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import AutoMinorLocator
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from extract_lightcurve_from_swd import extract_lightcurve
from plot_x_grid_vs_obs import MODEL_STYLE, style_for

GRID_DIR = '/ptmp/sbrennan/SN2024pba/SN_B_X_grid'

MODELS = ['X', 'X_M4', 'X_M4_rs5w3', 'X_M4_rs6w3', 'X_M4_rs5n04']

# labels tuned for this CSM study
LABELS = {
    'X': 'fiducial X',
    'X_M4': r'$M_{\rm sh}=4\,M_\odot$ (single CSM)',
    'X_M4_rs5w3': (r'rs$=5{\times}10^{15}$, $w_{\rm out}'
                   r'=3{\times}10^{17}$'),
    'X_M4_rs6w3': (r'rs$=6{\times}10^{15}$, $w_{\rm out}'
                   r'=3{\times}10^{17}$'),
    'X_M4_rs5n04': (r'rs$=5{\times}10^{15}$, $n_{\rm out}=0.4$'),
}

plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['DejaVu Serif'],
    'mathtext.fontset': 'dejavuserif',
    'font.size': 8,
    'axes.labelsize': 8,
    'legend.fontsize': 6.5,
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    'axes.linewidth': 0.7,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.top': True,
    'ytick.right': True,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})


def find_swd(mdir):
    name = os.path.basename(mdir)
    cands = [
        os.path.join(mdir, 'stella', f'{name}.swd'),
        os.path.join(mdir, 'stella', 'strad', f'{name}.swd'),
    ]
    for c in cands:
        if os.path.exists(c):
            return c
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--grid-dir', default=GRID_DIR)
    ap.add_argument('--out', default=None)
    ap.add_argument('--kappa', type=float, default=0.34)
    args = ap.parse_args()

    out = args.out or os.path.join(args.grid_dir, 'x_csm2_bolometric.png')

    fig, axes = plt.subplots(3, 1, figsize=(3.54, 5.2), sharex=True)
    handles = []

    for name in MODELS:
        mdir = os.path.join(args.grid_dir, name)
        swd = find_swd(mdir)
        if swd is None:
            print(f'  {name}: no .swd, skipped')
            continue
        lc = extract_lightcurve(swd, kappa=args.kappa)
        t = lc['time']
        # .tt exists only once STELLA finishes; its absence means the
        # run is still in progress and the .swd curve is partial
        partial = not os.path.exists(
            os.path.join(mdir, 'stella', f'{name}.tt'))
        st = style_for(name)
        lab = LABELS.get(name, st['label'])
        if partial:
            lab += ' (running)'
        print(f'  {name}: {len(t)} epochs, {t[0]:.1f}-{t[-1]:.1f} d')

        for ax, key in zip(axes, ['L_bol', 'Tbb', 'rbb']):
            ax.plot(t, lc[key], color=st['color'], lw=1.0, ls=st['ls'],
                    zorder=st['zorder'])
            if partial:
                ax.plot(t[-1], lc[key][-1], 'o', ms=3, color=st['color'],
                        zorder=st['zorder'] + 1)
        handles.append(Line2D([], [], color=st['color'], ls=st['ls'],
                              lw=1.0, label=lab))

    axes[0].set_yscale('log')
    axes[0].set_ylabel(r'$L_{\rm bol}$ (erg s$^{-1}$)')
    axes[1].set_yscale('log')
    axes[1].set_ylabel(r'$T_{\rm bb}$ (K)')
    axes[2].set_ylabel(r'$R_{\rm bb}$ (cm)')
    axes[2].set_xlabel('days since shell launch')
    axes[2].set_xlim(-15, 170)

    for ax in axes:
        ax.grid(True, alpha=0.3)
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        if ax.get_yscale() != 'log':
            ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.grid(True, which='minor', alpha=0.12)

    fig.legend(handles=handles, fontsize=6, ncol=2, frameon=False,
               loc='upper center', bbox_to_anchor=(0.5, 1.0),
               columnspacing=0.8, handlelength=1.8,
               handletextpad=0.5, labelspacing=0.35)
    fig.align_ylabels(axes)
    plt.tight_layout(rect=[0, 0, 1, 0.88])

    fig.savefig(out, dpi=300, bbox_inches='tight')
    pdf = os.path.splitext(out)[0] + '.pdf'
    fig.savefig(pdf, bbox_inches='tight')
    print(f'Saved: {out}')
    print(f'Saved: {pdf}')


if __name__ == '__main__':
    main()
