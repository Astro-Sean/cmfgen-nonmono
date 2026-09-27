#!/usr/bin/env python3
"""Plot ALL SN_A-family STELLA lightcurves vs SN 2024pba first-peak obs.

Discovers every model under SN_A_grid/ that has a STELLA .tt
lightcurve (the only source used - the pysynphot *_lightcurve*.csv
products are never read) and overlays their g/r/i lightcurves on the
dereddened SN 2024pba photometry.

Same conventions as plot_sn_a_pba_grid_vs_obs.py:
  DM = 34.52 (80.4 Mpc), E(B-V) = 0.15, A_g/A_r/A_i = 0.495/0.343/0.255
  x-axis = days since i-band peak (each model's i peak at 0, observed
  i peak at MJD 60507.493) unless --no-align-peak.

Observed points follow the superbol band colours (g=green, r=red,
i=goldenrod), filled with black edges.  Model lines are coloured per
model name (tab20 cycle).

Usage:
    python scripts/plot_sn_a_all_models.py [--no-align-peak] [--out FIG.png]
"""
import argparse
import glob
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_sn_a_pba_grid_vs_obs import (BASE, OBS_FILE, EXP_MJD,
                                     OBS_PEAK_MJD, PEAK_BAND, BAND_COLORS,
                                     load_model, load_observations,
                                     score_model)
from plot_model import read_2011dh_photometry

GRID_NEW_DIR = os.path.join(BASE, 'SN_A_grid')
OUT_FILE = os.path.join(GRID_NEW_DIR, 'sn_a_all_models_vs_obs.png')
OBS_2011DH_DIR = os.path.join(BASE, 'observational_data', 'SN2011dh')
# SN 2011dh Johnson band -> panel (same convention as plot_model.py:
# 2024pba g on V, r on R, i on I)
DH_BAND_MAP = {'g': 'V', 'r': 'R', 'i': 'I'}
DH_COLORS = {'V': 'yellowgreen', 'R': 'crimson', 'I': 'chocolate'}
CO_LINE_YMID = -16.4  # M_r of the 56Co-decay line midpoint (at day 40)
LEAK_LINE_YMID = -16.5  # M_r of the leakage-line midpoint (at day 80)
LEAK_SLOPE = 0.021  # mag/day, 2011dh R-band tail rate at ~100d (Ergon+15, Table 1)


def discover_all_sn_a_models():
    """All model dirs under SN_A_grid/ with lightcurves."""
    dirs = [d for d in sorted(glob.glob(os.path.join(GRID_NEW_DIR, '*')))
            if os.path.isdir(d)]
    models = []
    for d in dirs:
        df, src = load_model(d)
        if df is not None and len(df):
            models.append({'name': os.path.basename(d), 'dir': d,
                           'df': df, 'src': src})
    return models


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--exp-mjd', type=float, default=EXP_MJD,
                    help='Explosion MJD for the obs phase axis when '
                         f'--no-align-peak is used (default {EXP_MJD})')
    ap.add_argument('--align-peak', action=argparse.BooleanOptionalAction,
                    default=True,
                    help='Peak-align all lightcurves at the r-band peak '
                         '(default). --no-align-peak uses explosion date.')
    ap.add_argument('--out', default=OUT_FILE)
    args = ap.parse_args()

    obs = load_observations(OBS_FILE)
    for band in ('g', 'r', 'i'):
        for k in ('det', 'nd'):
            obs[band][k] = obs[band][k].copy()
            ref = OBS_PEAK_MJD if args.align_peak else args.exp_mjd
            obs[band][k]['phase'] = obs[band][k]['MJD'] - ref

    models = discover_all_sn_a_models()
    if not models:
        print('No SN_A-family models with lightcurves found')
        return
    print(f'Loaded {len(models)} SN_A-family models:\n')
    print(f"{'Model':<22} {'src':>6} {'i_peak':>7} {'pk_day':>6} "
          f"{'rms_g':>6} {'rms_r':>6} {'rms_i':>6} {'rms':>6}")
    print('-' * 80)
    for m in models:
        df = m['df']
        t_ipeak = df.loc[df[PEAK_BAND].idxmin(), 'time']
        m['shift'] = -t_ipeak if args.align_peak else 0.0
        sc = score_model(df, obs, m['shift'])
        m['score'] = sc
        print(f"{m['name']:<22} {m['src']:>6} {df[PEAK_BAND].min():7.2f} "
              f"{t_ipeak:6.1f} {sc.get('g', float('nan')):6.2f} "
              f"{sc.get('r', float('nan')):6.2f} "
              f"{sc.get('i', float('nan')):6.2f} {sc['combined']:6.2f}")

    # SN 2011dh photometry, shifted to its own I-band peak and median-binned
    # into 1-day bins (multiple nightly points -> one per day)
    dh = None
    dh_shift = 0.0
    if os.path.exists(os.path.join(OBS_2011DH_DIR, 'tablec4.dat')):
        dh = read_2011dh_photometry(OBS_2011DH_DIR)
        valid_i = ~np.isnan(dh['I'])
        if valid_i.any():
            dh_shift = dh['phase'][valid_i][np.argmin(dh['I'][valid_i])]

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

    # Neutral muted-pastel model colours so the band-coloured data stands out
    model_colors = sns.color_palette('muted', n_colors=max(len(models), 1))
    fig, axes = plt.subplots(3, 1, figsize=(11, 12), sharex=True)
    for ax, band, ylabel in ((axes[0], 'g', 'M_g (AB)'),
                             (axes[1], 'r', 'M_r (AB)'),
                             (axes[2], 'i', 'M_i (AB)')):
        det, nd = obs[band]['det'], obs[band]['nd']
        bc = BAND_COLORS[band]
        ax.errorbar(det['phase'], det['Mabs'], yerr=det['Mabs_err'],
                    fmt='o', color=bc, mfc=bc, mec='k', ms=8, capsize=4,
                    zorder=10)
        if len(nd):
            ax.scatter(nd['phase'], nd['Mabs_lim'], marker='v',
                       facecolors=bc, edgecolors='k', s=50,
                       zorder=9)
        if dh is not None:
            dh_band = DH_BAND_MAP[band]
            valid = ~np.isnan(dh[dh_band])
            bx, by = bin_1day(dh['phase'][valid] - dh_shift,
                              dh[dh_band][valid])
            ax.plot(bx, by, 'o', mfc='none', mec=DH_COLORS[dh_band],
                    ms=3.5, alpha=0.7, zorder=8)
        compact_h, ext_h = [], []
        for i, m in enumerate(models):
            compact = (m['name'] == 'SN_A' or '_S_' in m['name']
                       or 'pba' in m['name'])
            line, = ax.plot(m['df']['time'] + m['shift'], m['df'][band],
                            color=model_colors[i], lw=1.4, alpha=0.9,
                            ls='-' if compact else '--', label=m['name'])
            (compact_h if compact else ext_h).append(line)
        # 56Co->56Fe decay slope reference on the r panel only
        # (full gamma trapping, T=111.3 d -> 0.0098 mag/day)
        if band == 'r':
            co_slope = 2.5 * np.log10(np.e) / 111.26  # mag/day
            xs = np.array([30.0, 50.0])
            ax.plot(xs, CO_LINE_YMID + co_slope * (xs - 40.0), '--',
                    color='grey', lw=1.2, zorder=5)
            # Partial gamma-ray trapping (measured 2011dh tail slope)
            xs2 = np.array([50.0, 110.0])
            ax.plot(xs2, LEAK_LINE_YMID + LEAK_SLOPE * (xs2 - 80.0), ':',
                    color='grey', lw=1.6, zorder=5)
        ax.set_ylabel(ylabel)
        ax.set_ylim(-13.8, -17.6)
        ax.grid(alpha=0.3)
        lc = ax.legend(handles=compact_h, loc='upper right',
                       bbox_to_anchor=(0.78, 1.0), frameon=False,
                       fontsize=6.5, title='Compact', title_fontsize=7)
        ax.add_artist(lc)
        ax.legend(handles=ext_h, loc='upper right',
                  bbox_to_anchor=(1.0, 1.0), frameon=False,
                  fontsize=6.5, title='Extended', title_fontsize=7)
        if ax is axes[0]:
            g_compact_h, g_ext_h = compact_h, ext_h

    # Figure-level legend: which points are which supernova
    sn_handles = [
        Line2D([], [], marker='o', color='grey', mfc='grey', mec='k',
               ls='', ms=6, label='SN 2024pba (SN_A)'),
        Line2D([], [], marker='v', color='grey', mfc='grey', mec='k',
               ls='', ms=6, label='SN 2024pba 5-sig lim'),
        Line2D([], [], marker='o', mfc='none', mec='grey',
               ls='', ms=6, label='SN 2011dh'),
        Line2D([], [], color='grey', lw=1.4, label='STELLA model'),
    ]
    # One row, directly above the top panel, flush with its upper spine
    leg_sn = axes[0].legend(handles=sn_handles, loc='lower center', ncol=4,
                            frameon=False, fontsize=9,
                            bbox_to_anchor=(0.5, 1.0))
    axes[0].add_artist(leg_sn)
    # Re-add the g-panel model legends (the call above replaced them)
    lc = axes[0].legend(handles=g_compact_h, loc='upper right',
                        bbox_to_anchor=(0.78, 1.0), frameon=False,
                        fontsize=6.5, title='Compact', title_fontsize=7)
    axes[0].add_artist(lc)
    axes[0].legend(handles=g_ext_h, loc='upper right',
                   bbox_to_anchor=(1.0, 1.0), frameon=False,
                   fontsize=6.5, title='Extended', title_fontsize=7)

    if args.align_peak:
        axes[2].set_xlabel(f'Days since {PEAK_BAND}-band peak')
        axes[2].set_xlim(-30, 180)
    else:
        axes[2].set_xlabel(f'Days since explosion (MJD {args.exp_mjd:.1f})')
        axes[2].set_xlim(-5, 90)

    plt.tight_layout(rect=[0, 0, 1, 0.94])
    # Angled labels along the two reference slopes on the r panel
    fig.canvas.draw()
    ax = axes[1]

    def angled_label(slope, ymid, xmid, xa, xb, text):
        ya = ymid + slope * (xa - xmid)
        yb = ymid + slope * (xb - xmid)
        pa = ax.transData.transform((xa, ya))
        pb = ax.transData.transform((xb, yb))
        ang = np.degrees(np.arctan2(pb[1] - pa[1], pb[0] - pa[0]))
        n = np.array([-(pb[1] - pa[1]), pb[0] - pa[0]])
        n /= np.hypot(*n)
        if n[1] < 0:
            n = -n
        mid = ax.transData.inverted().transform((pa + pb) / 2 + 10 * n)
        ax.text(mid[0], mid[1], text, rotation=ang, rotation_mode='anchor',
                ha='center', va='bottom', fontsize=8, color='grey')

    co_slope = 2.5 * np.log10(np.e) / 111.26
    angled_label(co_slope, CO_LINE_YMID, 40.0, 30.0, 50.0,
                 r'$^{56}$Co radioactive decay')

    angled_label(LEAK_SLOPE, LEAK_LINE_YMID, 80.0, 50.0, 110.0,
                 r'$\gamma$-ray leakage (Ergon+15)')
    plt.savefig(args.out, dpi=150)
    print(f'\nPlot saved: {args.out}')


if __name__ == '__main__':
    main()
