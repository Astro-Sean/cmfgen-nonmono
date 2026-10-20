#!/usr/bin/env python3
"""A&A-style plot of the best SN_A models vs SN 2024pba first peak.

Same formatting as plot_x_grid_vs_obs.py: 2x2 panel layout, serif
fonts, inward ticks, model legend above the figure.  All models under
SN_A_grid/ (and SN_A_Y_grid/ once their .tt files land) are scored
with the focused metric from plot_sn_a_pba_grid_vs_obs.py (i-band
non-detection rise gate + peak/<15 d weighting) and only the
N_BEST best-fitting are drawn.

Panels: g, r, i absolute mags plus an i-band zoom on the rise window
(-35..+15 d) showing the two non-detection limits used by the metric.
Phase axis: days since the observed i-band peak (MJD 60507.493); each
model's own i peak is placed at 0.

Usage:
    python scripts/plot_sn_a_best_vs_obs.py [--n-best 5] [--out FIG.png]
"""

import argparse
import glob
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import ConnectionPatch, Rectangle
from matplotlib.ticker import AutoMinorLocator
import numpy as np
import pandas as pd


sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_sn_a_pba_grid_vs_obs import (BASE, OBS_FILE, OBS_PEAK_MJD,
                                     PEAK_BAND, EXP_MJD, A_R, DM,
                                     load_model, load_observations,
                                     score_model, score_model_shift)
from plot_model import read_2011dh_photometry

MODEL_DIRS = [os.path.join(BASE, 'SN_A_Y_grid')]
OUT_FILE = os.path.join(BASE, 'SN_A_Y_grid', 'sn_a_best_vs_obs.png')
OBS_2011DH_DIR = os.path.join(BASE, 'observational_data', 'SN2011dh')
# GOTO forced photometry in the broad L filter (~400-700 nm); closest
# STELLA band is r, so r-band extinction applies.
GOTO_FILE = os.path.join(BASE, 'observational_data',
                         'SN2024pba_GOTO_forced_photometry.csv')
GOTO_COLOR = 'purple'
PLOT_GOTO = False
ATLAS_FILE = os.path.join(BASE, 'observational_data',
                          'SN2024pba_atclean_cleaned.csv')
A_O = 2.350 * 0.15  # same coefficient as the SN_B figure
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

# Short legend labels for the known best families; anything else falls
# back to the directory name.
MODEL_LABEL = {
    'SN_A_g_T20k_Ni075': r'$E=9$, $M_{\rm env}=0.05$, $R=200$, Ni 0.075',
    'SN_A_g_T20k_Ni06': r'$E=9$, $M_{\rm env}=0.05$, $R=200$, Ni 0.06',
    'SN_A_g_T15k_Ni06': r'$E=9$, $T_{\rm base}=15$k, Ni 0.06',
    'SN_A_g_M10': r'$E=9$, $M_{\rm env}=0.10$, Ni 0',
    'SN_A_g_T20k_Ni05': r'$E=9$, $M_{\rm env}=0.05$, $R=200$, Ni 0.05',
    'SN_A_g_T15k_Ni05': r'$E=9$, $T_{\rm base}=15$k, Ni 0.05',
    'SN_A_g_R300': r'$E=9$, $R=300$, Ni 0',
    'SN_A_S_11dh': r'2011dh stripped, Ni 0.075',
}
RANK_COLOR = ['k', 'tab:red', 'tab:blue', 'tab:green', 'tab:purple']
RANK_LS = ['-', '--', '-.', ':', '-']
# Score-equivalent duplicate runs that should not take two legend slots.
# S_11dh is excluded on request: it does not match the lightcurve well.
DUP_NAMES = {'SN_A_g_T20k_Ni06_r2', 'SN_A_S_11dh'}


def discover_models():
    models = []
    for root in MODEL_DIRS:
        for mdir in sorted(glob.glob(os.path.join(root, '*'))):
            name = os.path.basename(mdir)
            if not os.path.isdir(mdir) or name in DUP_NAMES:
                continue
            df, src = load_model(mdir)
            if df is not None and len(df):
                models.append({'name': name, 'df': df, 'src': src,
                               'grid': os.path.basename(root)})
    return models


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
    ap.add_argument('--n-best', type=int, default=5)
    ap.add_argument('--out', default=OUT_FILE)
    args = ap.parse_args()

    obs = load_observations(OBS_FILE)
    if os.path.exists(ATLAS_FILE):
        atl = pd.read_csv(ATLAS_FILE)
        sub = atl[atl['F'] == 'o'].copy()
        sub['sn'] = sub['uJy'] / sub['duJy_new']
        good = (sub['Mask'].isin(['0x400800', '0x400900'])
                & (sub['chi/N'].isna() | (sub['chi/N'] < 5.0)))
        o_det = sub[(sub['sn'] >= 3.0) & good].copy()
        o_det['Mabs_err'] = o_det['dm'].where(o_det['dm'].notna(),
                                            1.086 / o_det['sn'])
        o_det['phase'] = o_det['MJD'] - OBS_PEAK_MJD
        o_det['Mabs'] = o_det['m'] - DM - A_O
        o_det['night'] = np.floor(o_det['phase'])
        rows = []
        for night, grp in o_det.groupby('night'):
            w = 1.0 / grp['Mabs_err']**2
            rows.append({
                'phase': np.average(grp['phase'], weights=w),
                'Mabs': np.average(grp['Mabs'], weights=w),
                'Mabs_err': 1.0 / np.sqrt(w.sum()),
            })
        obs['o'] = {'det': pd.DataFrame(rows), 'nd': pd.DataFrame()}
    else:
        obs['o'] = {'det': pd.DataFrame(), 'nd': pd.DataFrame()}
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
    for m in discover_models():
        df = m['df']
        tpk = df.loc[df[PEAK_BAND].idxmin(), 'time']
        # Bounded extra shift (+-5 d, 0.25 d step): the observed i peak
        # is only bracketed to a few days, so a small per-model shift is
        # physical.  score_model_shift returns the best score and delta.
        sc, delta = score_model_shift(df, obs, -tpk)
        if np.isnan(sc['combined']):
            continue
        scored.append(dict(m, tpk=tpk, dt=delta, score=sc))
        print(f"  {m['name']:<34} dt={delta:+5.2f} comb={sc['combined']:.2f} "
              f"rise={sc.get('rise_rms', np.nan):.2f}")
    scored.sort(key=lambda s: s['score']['combined'])
    # Ni075 is the adopted reference model: its i-band match is better
    # and its stripped envelope is the physically preferred structure,
    # so it takes the best slot even though S_11dh edges it on g/r.
    for i, s in enumerate(scored):
        if s['name'] == 'SN_A_g_T20k_Ni075':
            scored.insert(0, scored.pop(i))
            break
    # Remaining SN_A_Y_grid models are drawn thin and grey behind the
    # highlighted best set so the whole Y grid is visible at a glance.
    y_extra = [s for s in scored[args.n_best:]
               if s['grid'] == 'SN_A_Y_grid']
    scored = scored[:args.n_best]
    print(f'plotting best {len(scored)} of the scored models')

    BANDS = [
        ('g', 'o', 'tab:green', 'SN 2024pba $g$'),
        ('r', 'o', 'tab:red', 'SN 2024pba $r$'),
        ('i', 's', 'goldenrod', 'SN 2024pba $i$'),
        ('o', '^', 'darkorange', 'SN 2024pba $o$ (ATLAS)'),
    ]
    x_lo, x_hi = -27.0, 60.0

    fig, axes = plt.subplots(2, 2, figsize=(2.0 * AA_COL_IN, 4.6),
                             sharex=True)
    band_ax = dict(zip('grio', axes.flat))
    ax_zoom = band_ax['i'].inset_axes([0.35, 0.10, 0.32, 0.45])
    ax_zoom.set_facecolor('white')
    ax_zoom.patch.set_alpha(1.0)

    for band, fmt, col, lab in BANDS:
        ax = band_ax[band]
        det = obs[band]['det']
        ax.errorbar(det['phase'], det['Mabs'], yerr=det['Mabs_err'],
                    fmt=fmt, color=col, ms=3, mec='k', mew=0.5,
                    capsize=1.0, elinewidth=0.6, alpha=0.85, zorder=20)
        ax.set_ylim(top=-17.6, bottom=-13.9)
        ax.set_xlim(x_lo, x_hi)
        ax.text(0.03, 0.06, f'${band}$', transform=ax.transAxes,
                fontsize=8, color=col)
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.grid(True, alpha=0.3)
        ax.grid(True, which='minor', alpha=0.12)
    for ax in axes[:, 0]:
        ax.set_ylabel('Absolute magnitude')

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

    # GOTO forced photometry (broad L band, ~r) on the r panel.
    if PLOT_GOTO and os.path.exists(GOTO_FILE):
        goto = pd.read_csv(GOTO_FILE)
        gdet = goto[goto['forced_mag'].notna()].copy()
        gdet['phase'] = gdet['mjd'] - OBS_PEAK_MJD
        gdet['Mabs'] = gdet['forced_mag'] - DM - A_R
        band_ax['r'].errorbar(gdet['phase'], gdet['Mabs'],
                              yerr=gdet['forced_mag_uncert'], fmt='D',
                              color=GOTO_COLOR, ms=2.8, mec='k',
                              mew=0.4, capsize=0.8, elinewidth=0.5,
                              alpha=0.85, zorder=20)
        glim = goto[goto['forced_mag'].isna()].copy()
        glim['phase'] = glim['mjd'] - OBS_PEAK_MJD
        glim['Mabs'] = glim['forced_mag_sigma_limit'] - DM - A_R
        glim = glim[np.isfinite(glim['Mabs'])
                    & (glim['phase'] > x_lo) & (glim['phase'] < x_hi)]
        band_ax['r'].plot(glim['phase'], glim['Mabs'], 'v', ms=3.5,
                          mfc='white', mec=GOTO_COLOR, mew=0.8,
                          alpha=0.9, zorder=19)

    # i-band zoom on the last non-detection and first detection.
    ax_zoom.errorbar(obs['i']['det']['phase'], obs['i']['det']['Mabs'],
                     yerr=obs['i']['det']['Mabs_err'], fmt='s',
                     color='goldenrod', ms=3, mec='k', mew=0.5,
                     capsize=1.0, elinewidth=0.6, alpha=0.85, zorder=20)
    zx_lo, zx_hi = -21.0, -16.5
    ax_zoom.set_ylim(top=-15.4, bottom=-14.7)
    ax_zoom.set_xlim(zx_lo, zx_hi)
    ax_zoom.xaxis.set_minor_locator(AutoMinorLocator())
    ax_zoom.yaxis.set_minor_locator(AutoMinorLocator())




    ax_zoom.grid(True, alpha=0.3)
    ax_zoom.tick_params(labelsize=5.5)
    for tick_lab in list(ax_zoom.get_yticklabels()):
        tick_lab.set_bbox({'facecolor': 'white', 'edgecolor': 'none',
                           'pad': 0.8, 'alpha': 0.95})
    # Zoom window on the i panel, connectors to the inset (same
    # convention as the r-band insets on the SN_B figure).
    main_i = band_ax['i']
    main_i.add_patch(Rectangle(
        (zx_lo, -15.5), zx_hi - zx_lo, 0.8, fill=False, edgecolor='k',
        lw=0.5, zorder=30))
    for data_xy, inset_xy in (((zx_hi, -14.7), (0.0, 0.0)),
                              ((zx_hi, -15.5), (0.0, 1.0))):
        connector = ConnectionPatch(
            xyA=data_xy, coordsA='data', axesA=main_i,
            xyB=inset_xy, coordsB='axes fraction', axesB=ax_zoom,
            color='k', lw=0.5, zorder=29,
            clip_on=False)
        fig.add_artist(connector)

    handles = []
    for rank, s in enumerate(scored):
        col = RANK_COLOR[rank % len(RANK_COLOR)]
        ls = RANK_LS[rank % len(RANK_LS)]
        lw = 1.8 if rank == 0 else 1.1
        lab = MODEL_LABEL.get(s['name'], s['name'])
        if rank == 0:
            lab += ' (best)'
        ph = s['df']['time'] - s['tpk'] + s['dt']
        for band, _, _, _ in BANDS:
            if band not in s['df']:
                continue
            band_ax[band].plot(ph, s['df'][band], color=col, lw=lw,
                               ls=ls, zorder=10 + rank)
        ax_zoom.plot(ph, s['df']['i'], color=col, lw=lw, ls=ls,
                     zorder=10 + rank)
        handles.append(Line2D([], [], color=col, lw=lw, ls=ls,
                              label=lab))

    for s in y_extra:
        ph = s['df']['time'] - s['tpk'] + s['dt']
        for band, _, _, _ in BANDS:
            if band not in s['df']:
                continue
            band_ax[band].plot(ph, s['df'][band], color='0.7', lw=0.7,
                               alpha=0.6, zorder=5)
        ax_zoom.plot(ph, s['df']['i'], color='0.7', lw=0.7,
                     alpha=0.6, zorder=5)
    if y_extra:
        handles.append(Line2D([], [], color='0.7', lw=0.7,
                              label='other $Y$-grid'))

    fig.legend(handles=handles, fontsize=5.4, ncol=2, frameon=False,
               loc='upper center', bbox_to_anchor=(0.5, 0.94),
               columnspacing=0.6, handlelength=1.6,
               handletextpad=0.4)

    for ax in axes[1]:
        ax.set_xlabel('Days since SN 2024pba $i$-band peak')
    fig.align_ylabels(axes)
    n_leg_rows = max(1, -(-len(handles) // 2))
    rect_top = min(0.92, 1.0 - 0.045 * n_leg_rows)
    plt.tight_layout(rect=[0, 0, 1, rect_top])
    save_both(fig, args.out)


if __name__ == '__main__':
    main()
