#!/usr/bin/env python3
"""A&A-style plot of model Y vs SN 2024pba first peak.

Model Y = SN_A_g_T20k_Ni075, the best SN_A fit under the focused
metric (i-band non-detection rise gate + peak/<15 d weighting).
Alongside it, a selection of models that each differ from model Y by
one or two parameters only (Ni56 mass, envelope mass/radius/T_base,
explosion energy), read as a one-parameter-at-a-time study.

Same formatting as plot_x_grid_vs_obs.py: 2x2 panel layout, serif
fonts, inward ticks, model legend above the figure.  Panels: g, r, i
absolute mags plus an i-band zoom on the rise window showing the
non-detection limits used by the metric.  Phase axis: days since the
observed i-band peak (MJD 60507.493); each model's own i peak is
placed at 0 plus its best bounded shift (+-5 d).

Usage:
    python scripts/plot_sn_a_modelY_vs_obs.py [--out FIG.png]
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
from plot_sn_a_pba_grid_vs_obs import (BASE, OBS_FILE, OBS_PEAK_MJD,
                                     PEAK_BAND, EXP_MJD,
                                     load_model, load_observations,
                                     score_model, score_model_shift)
from plot_model import read_2011dh_photometry

MODEL_DIRS = [os.path.join(BASE, 'SN_A_Y_grid'),
              os.path.join(BASE, 'SN_A_Y_grid')]
OUT_FILE = os.path.join(BASE, 'SN_A_Y_grid', 'sn_a_modelY_vs_obs.png')
OBS_2011DH_DIR = os.path.join(BASE, 'observational_data', 'SN2011dh')
# SN 2011dh Johnson band -> panel (2024pba g on V, r on R, i on I)
DH_BAND_MAP = {'g': 'V', 'r': 'R', 'i': 'I'}
DH_COLORS = {'V': 'yellowgreen', 'R': 'crimson', 'I': 'chocolate'}

# A&A single column is ~9 cm (3.54 in) wide.
AA_COL_IN = 3.54
plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['DejaVu Serif'],
    'mathtext.fontset': 'dejavuserif',
    'font.size': 8,
    'axes.labelsize': 8,
    'axes.titlesize': 8,
    'legend.fontsize': 6,
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

# Model selection: model Y plus variants differing by <= 2 parameters.
# Model Y = E 9x10^50 erg, M_env 0.05 Msun at 200 Rsun (X_H=0.7,
# T_base 20 kK), Ni56 0.075 Msun.  (name, grid dir, label, style)
MODELS = [
    ('SN_A_g_T20k_Ni075', 'SN_A_Y_grid',
     r'model Y',
     dict(color='k', lw=1.9, ls='-', zorder=19)),
    # Ni56 axis (same envelope as Y)
    ('SN_A_g_T20k_Ni06', 'SN_A_Y_grid',
     r'$M_{\rm Ni}=0.06\,M_\odot$',
     dict(color='tab:blue', lw=1.1, ls='-', zorder=12)),
    ('SN_A_g_T20k_Ni05', 'SN_A_Y_grid',
     r'$M_{\rm Ni}=0.05\,M_\odot$',
     dict(color='tab:blue', lw=1.1, ls='--', zorder=12)),
    ('SN_A_g_T20k_Ni04', 'SN_A_Y_grid',
     r'$M_{\rm Ni}=0.04\,M_\odot$',
     dict(color='tab:blue', lw=1.1, ls=':', zorder=12)),
    # Envelope temperature
    ('SN_A_g_T15k_Ni06', 'SN_A_Y_grid',
     r'$T_{\rm base}=15$ kK, $M_{\rm Ni}=0.06$',
     dict(color='tab:red', lw=1.1, ls='-', zorder=12)),
    # Envelope mass
    ('SN_A_g_M10', 'SN_A_Y_grid',
     r'$M_{\rm env}=0.10\,M_\odot$, $M_{\rm Ni}=0$',
     dict(color='tab:green', lw=1.1, ls='-', zorder=12)),
    # Envelope radius
    ('SN_A_g_R300', 'SN_A_Y_grid',
     r'$R_{\rm env}=300\,R_\odot$, $M_{\rm Ni}=0$',
     dict(color='tab:purple', lw=1.1, ls='-', zorder=12)),
    # Compact envelope (R 200 -> 10 Rsun, Ni -> 0.05)
    ('Y_C05', 'SN_A_Y_grid',
     r'compact: $R_{\rm env}=10\,R_\odot$, $M_{\rm Ni}=0.05$',
     dict(color='tab:orange', lw=1.1, ls='-', zorder=12)),
    # Explosion energy (E 9 -> 15, Ni -> 0.05)
    ('Y_E15', 'SN_A_Y_grid',
     r'$E=15\times10^{50}$ erg, $M_{\rm Ni}=0.05$',
     dict(color='tab:cyan', lw=1.1, ls='-', zorder=12)),
]


def bin_1day(phase, mag, width=1.0):
    edges = np.arange(np.floor(phase.min()),
                      np.ceil(phase.max()) + width, width)
    idx = np.digitize(phase, edges)
    bx, by = [], []
    for k in range(1, len(edges)):
        sel = idx == k
        if sel.any():
            bx.append(np.mean(phase[sel]))
            by.append(np.median(mag[sel]))
    return np.array(bx), np.array(by)


def save_both(fig, out_png):
    fig.savefig(out_png, dpi=300)
    pdf = os.path.splitext(out_png)[0] + '.pdf'
    fig.savefig(pdf)
    print(f'saved {out_png} and {pdf}')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default=OUT_FILE)
    args = ap.parse_args()

    obs = load_observations(OBS_FILE)
    for band in ('g', 'r', 'i'):
        for k in ('det', 'nd'):
            obs[band][k] = obs[band][k].copy()
            obs[band][k]['phase'] = \
                obs[band][k]['MJD'] - OBS_PEAK_MJD

    # SN 2011dh photometry, shifted to its own I-band peak and
    # median-binned into 1-day bins (multiple nightly points -> 1/day)
    dh = None
    dh_shift = 0.0
    if os.path.exists(os.path.join(OBS_2011DH_DIR, 'tablec4.dat')):
        dh = read_2011dh_photometry(OBS_2011DH_DIR)
        valid_i = ~np.isnan(dh['I'])
        if valid_i.any():
            dh_shift = dh['phase'][valid_i][np.argmin(dh['I'][valid_i])]

    scored = []
    for name, grid, label, style in MODELS:
        mdir = os.path.join(BASE, grid, name)
        df, src = load_model(mdir)
        if df is None or not len(df):
            print(f'  {name:<34} no .tt yet - skipped')
            continue
        tpk = df.loc[df[PEAK_BAND].idxmin(), 'time']
        # Bounded extra shift (+-5 d, 0.25 d step): the observed i peak
        # is only bracketed to a few days, so a small per-model shift is
        # physical.  score_model_shift returns the best score and delta.
        sc, delta = score_model_shift(df, obs, -tpk)
        scored.append(dict(name=name, df=df, label=label, style=style,
                           tpk=tpk, dt=delta, score=sc))
        print(f"  {name:<34} dt={delta:+5.2f} comb={sc['combined']:.2f} "
              f"rise={sc.get('rise_rms', np.nan):.2f}")
    print(f'plotting {len(scored)} models')

    BANDS = [
        ('g', 'o', 'tab:green', 'SN 2024pba $g$'),
        ('r', 'o', 'tab:red', 'SN 2024pba $r$'),
        ('i', 's', 'goldenrod', 'SN 2024pba $i$'),
    ]
    x_lo, x_hi = -40.0, 60.0

    fig, axes = plt.subplots(2, 2, figsize=(2.0 * AA_COL_IN, 4.6))
    band_ax = dict(zip('gri', axes.flat[:3]))
    ax_zoom = axes.flat[3]

    for band, fmt, col, lab in BANDS:
        ax = band_ax[band]
        det = obs[band]['det']
        ax.errorbar(det['phase'], det['Mabs'], yerr=det['Mabs_err'],
                    fmt=fmt, color=col, ms=3, mec='k', mew=0.5,
                    capsize=1.0, elinewidth=0.6, alpha=0.85, zorder=20)
        ax.set_ylim(top=-18.2, bottom=-13.9)
        ax.set_xlim(x_lo, x_hi)
        ax.text(0.03, 0.06, f'${band}$', transform=ax.transAxes,
                fontsize=8)
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.grid(True, alpha=0.3)
        ax.grid(True, which='minor', alpha=0.12)
    band_ax['g'].set_ylabel('Absolute magnitude')
    ax_zoom.set_ylabel('Absolute magnitude')

    # Non-detection limits: white-faced inverted triangles edged in the
    # band colour; the i-band ones before first detection are the
    # metric's rise gate and also go on the zoom panel.
    for band, _, col, _ in BANDS:
        nd = obs[band]['nd']
        if not len(nd):
            continue
        band_ax[band].plot(nd['phase'], nd['Mabs_lim'], 'v', ms=3.5,
                           mfc='white', mec=col, mew=0.8, alpha=0.9,
                           zorder=19)
    nd_i = obs['i']['nd']
    nd_pre = nd_i[nd_i['phase'] < obs['i']['det']['phase'].min()]
    ax_zoom.plot(nd_pre['phase'], nd_pre['Mabs_lim'], 'v', ms=3.5,
                 mfc='white', mec='goldenrod', mew=0.8, alpha=0.9,
                 zorder=19)

    # SN 2011dh reference (open circles, 1-day median bins)
    if dh is not None:
        for band, _, _, _ in BANDS:
            dh_band = DH_BAND_MAP[band]
            valid = ~np.isnan(dh[dh_band])
            bx, by = bin_1day(dh['phase'][valid] - dh_shift,
                              dh[dh_band][valid])
            band_ax[band].plot(bx, by, 'o', mfc='none',
                               mec=DH_COLORS[dh_band], ms=3,
                               alpha=0.7, zorder=8)
        valid = ~np.isnan(dh['I'])
        bx, by = bin_1day(dh['phase'][valid] - dh_shift,
                          dh['I'][valid])
        ax_zoom.plot(bx, by, 'o', mfc='none', mec=DH_COLORS['I'],
                     ms=3, alpha=0.7, zorder=8)

    # 4th panel: i-band zoom on the rise/peak window that the metric
    # prioritises.
    ax_zoom.errorbar(obs['i']['det']['phase'], obs['i']['det']['Mabs'],
                     yerr=obs['i']['det']['Mabs_err'], fmt='s',
                     color='goldenrod', ms=3, mec='k', mew=0.5,
                     capsize=1.0, elinewidth=0.6, alpha=0.85, zorder=20)
    ax_zoom.set_ylim(top=-17.6, bottom=-13.2)
    ax_zoom.set_xlim(-35.0, 20.0)
    ax_zoom.text(0.03, 0.06, r'$i$ (rise zoom)',
                 transform=ax_zoom.transAxes, fontsize=8)
    ax_zoom.xaxis.set_minor_locator(AutoMinorLocator())
    ax_zoom.yaxis.set_minor_locator(AutoMinorLocator())
    ax_zoom.grid(True, alpha=0.3)
    ax_zoom.grid(True, which='minor', alpha=0.12)

    handles = []
    for rank, s in enumerate(scored):
        st = s['style']
        ph = s['df']['time'] - s['tpk'] + s['dt']
        for band, _, _, _ in BANDS:
            band_ax[band].plot(ph, s['df'][band], color=st['color'],
                               lw=st['lw'], ls=st['ls'],
                               zorder=st['zorder'])
        ax_zoom.plot(ph, s['df']['i'], color=st['color'], lw=st['lw'],
                     ls=st['ls'], zorder=st['zorder'])
        handles.append(Line2D([], [], color=st['color'], lw=st['lw'],
                              ls=st['ls'], label=s['label']))

    obs_handles = [
        Line2D([], [], marker='o', color='tab:green', ls='', ms=3,
               mec='k', mew=0.5, label='$g$'),
        Line2D([], [], marker='o', color='tab:red', ls='', ms=3,
               mec='k', mew=0.5, label='$r$'),
        Line2D([], [], marker='s', color='goldenrod', ls='', ms=3,
               mec='k', mew=0.5, label='$i$'),
        Line2D([], [], marker='v', color='goldenrod', ls='', ms=3.5,
               mfc='white', mec='goldenrod', mew=0.8,
               label='5$\\sigma$ lim'),
        Line2D([], [], marker='o', color='0.35', ls='', ms=3,
               mfc='none', mec='0.35', label='SN 2011dh'),
    ]
    band_ax['r'].legend(handles=obs_handles, fontsize=6, ncol=1,
                        frameon=False, loc='upper right',
                        handletextpad=0.4, borderaxespad=0.3,
                        labelspacing=0.3)
    fig.legend(handles=handles, fontsize=5.4, ncol=3, frameon=False,
               loc='upper center', bbox_to_anchor=(0.5, 1.0),
               columnspacing=0.6, handlelength=1.6,
               handletextpad=0.4)

    for ax in axes[0]:
        ax.tick_params(labelbottom=False)
    for ax in axes[-1]:
        ax.set_xlabel('Days since SN 2024pba $i$-band peak')
    fig.align_ylabels(axes)
    n_leg_rows = max(1, -(-len(handles) // 3))
    rect_top = min(0.92, 1.0 - 0.045 * n_leg_rows)
    plt.tight_layout(rect=[0, 0, 1, rect_top])
    save_both(fig, args.out)


if __name__ == '__main__':
    main()
