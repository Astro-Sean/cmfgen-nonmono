#!/usr/bin/env python3
"""Plot the observed SN 2024pba SN_B lightcurve vs SN 1994W
and SN 1998S.

Photometry, distance moduli, extinction and phasing conventions are
imported unchanged from plot_sn_b_grid_vs_obs.py:

  - SN 2024pba: PSF r/g (SN_B_lightcurve.dat) and ATLAS o, absolute
    mag = m - DM(34.52) - A_band, phase = MJD - 60850 (the observed
    second-peak date).  Non-detections shown as 5-sigma upper limits.
  - SN 1994W: V and R bands (OSC/astrocatalogs), phased to its own
    V-band peak.
  - SN 1998S: R band (OSC/astrocatalogs), phased to its own V-band
    peak.
  - Detections with photometric error > ERR_MAX are dropped.

Usage:
    python scripts/plot_sn_b_obs_vs_sn94w.py [--out PATH]
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_sn_b_grid_vs_obs import (ATLAS_FILE, OBS_FILE, SN94W_AV,
                                   SN94W_DM, SN94W_FILE, SNB_EXP_PHASE,
                                   load_atlas, load_observations)

SN94W_OSC_FILE = '/ptmp/sbrennan/SN2024pba/observational_data/SN1994W_osc.csv'
SN98S_OSC_FILE = '/ptmp/sbrennan/SN2024pba/observational_data/SN1998S_osc.csv'
OUT_DEFAULT = '/ptmp/sbrennan/SN2024pba/results/sn_b_obs_vs_sn94w.png'

# SN 1994W per-band extinction, scaled from A_V = 0.05 (SN94W_AV in
# plot_sn_b_grid_vs_obs) with standard R_V = 3.1 ratios.
SN94W_ABAND = {'V': SN94W_AV, 'B': 1.30 * SN94W_AV, 'R': 0.80 * SN94W_AV}

# SN 1998S (NGC 3877): D ~ 17 Mpc -> DM ~ 31.1; Galactic E(B-V) ~ 0.02
# -> A_V ~ 0.06 (Fassia et al. 2001 find little host reddening).
SN98S_DM = 31.1
SN98S_AV = 0.06
SN98S_ABAND = {'V': SN98S_AV, 'B': 1.30 * SN98S_AV, 'R': 0.80 * SN98S_AV}

# Detections with larger photometric errors are dropped.
ERR_MAX = 0.5


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', default=OUT_DEFAULT,
                    help=f'output figure path (default: {OUT_DEFAULT})')
    ap.add_argument('--xlim', nargs=2, type=float, default=None,
                    metavar=('XLO', 'XHI'), help='phase window (days)')
    args = ap.parse_args()

    obs = load_observations(OBS_FILE)
    obs['o'] = load_atlas(ATLAS_FILE, band='o')

    # SN 1994W: prefer the per-source OSC photometry (mjd, band, mag,
    # e_mag, source); fall back to the plain V-band CSV.  Phase zero =
    # the V-band peak of the homogeneous IAUC visual-estimate series.
    if os.path.exists(SN94W_OSC_FILE):
        w94raw = pd.read_csv(SN94W_OSC_FILE)
        vmain = w94raw[w94raw['source'] == 'IAUC_vis']
        vpk = w94raw[w94raw['band'] == 'V']
        ref = (vmain if len(vmain) else vpk)
        t0 = ref.loc[ref['mag'].idxmin(), 'mjd']
        w94raw['phase'] = w94raw['mjd'] - t0
        w94raw['Mabs'] = [m - SN94W_DM - SN94W_ABAND.get(b, 0.0)
                          for m, b in zip(w94raw['mag'], w94raw['band'])]
        w94 = w94raw
    else:
        w94 = pd.read_csv(SN94W_FILE)
        w94 = w94[w94['V'].notna()].copy()
        w94['band'] = 'V'
        w94['source'] = 'IAUC_vis'
        w94['e_mag'] = 0.2
        w94['phase'] = w94['mjd'] - w94.loc[w94['V'].idxmin(), 'mjd']
        w94['Mabs'] = w94['V'] - SN94W_DM - SN94W_AV
        w94 = (w94.assign(night=np.floor(w94['phase']))
                  .groupby('night')
                  .agg(Mabs=('Mabs', 'median'), e_mag=('e_mag', 'median'))
                  .reset_index()
                  .rename(columns={'night': 'phase'}))
        w94['band'] = 'V'
        w94['source'] = 'IAUC_vis'

    # SN 1998S (mjd, band, mag, e_mag, source).  Phase zero = the
    # Fassia et al. (2001) V-band peak.
    w98 = None
    if os.path.exists(SN98S_OSC_FILE):
        w98 = pd.read_csv(SN98S_OSC_FILE)
        vpk = w98[(w98['source'] == 'Fassia01') & (w98['band'] == 'V')]
        if not len(vpk):
            vpk = w98[w98['band'] == 'V']
        t0 = vpk.loc[vpk['mag'].idxmin(), 'mjd']
        w98['phase'] = w98['mjd'] - t0
        w98['Mabs'] = [m - SN98S_DM - SN98S_ABAND.get(b, 0.0)
                       for m, b in zip(w98['mag'], w98['band'])]

    fig, ax = plt.subplots(figsize=(10, 6.5))

    band_style = {'r': dict(color='tab:red', marker='o',
                            label='SN 2024pba r (PSF)'),
                  'g': dict(color='tab:green', marker='s',
                            label='SN 2024pba g (PSF)'),
                  'o': dict(color='darkorange', marker='^',
                            label='SN 2024pba o (ATLAS)')}
    for b, st in band_style.items():
        det = obs[b]['det']
        det = det[det['Mabs_err'] <= ERR_MAX]
        ax.errorbar(det['phase'], det['Mabs'], yerr=det['Mabs_err'],
                    fmt=st['marker'], color=st['color'], ms=5,
                    mec='k', mew=0.7, capsize=2, elinewidth=0.7,
                    alpha=0.85, zorder=20, label=st['label'])
        nd = obs[b]['nd']
        if len(nd):
            ax.scatter(nd['phase'], nd['Mabs_lim'], marker='v',
                       facecolors='none', edgecolors=st['color'],
                       s=30, alpha=0.5, zorder=15)

    # SN 1994W R band: the only published R photometry is the late-time
    # JKT/NOT CCD set from Sollerman, Cumming & Lundqvist (1998) - there
    # is no R-band plateau photometry, so the 94W comparison here covers
    # only the post-plateau decline.  Points with error > ERR_MAX are
    # dropped.
    SN94W_SOURCES = [
        ('IAUC_vis_early', 'V', 'o',
         'Early visual est. (IAUCs, 1994)', 0.30),
        ('IAUC_vis', 'V', 's',
         'Visual est. (IAUC 6072+, 1994)', 0.95),
        ('IAUC6044', 'V', '^',
         'Mikuz & Bragaglia (IAUC 6044, 1994)', 0.95),
        ('Tsvetkov95', 'V', 'D',
         'Tsvetkov (IBVS 4253, 1995)', 0.95),
    ]
    w94_handles = []
    w94_arrow_y = []
    arrow_lbl = None
    for src_key, band, marker, lbl, alpha in SN94W_SOURCES:
        sub = w94[(w94['source'] == src_key) & (w94['band'] == band)
                  & (w94['e_mag'] <= ERR_MAX)]
        if not len(sub):
            continue
        h = ax.errorbar(sub['phase'], sub['Mabs'], yerr=sub['e_mag'],
                        fmt=marker, color='0.35', mec='0.35',
                        mfc='white', ms=7 if alpha > 0.5 else 5,
                        mew=1.3, capsize=2.5, elinewidth=0.9,
                        ls='none', alpha=alpha, zorder=10, label=lbl)
        w94_handles.append(h)

    # The post-plateau collapse (SCL98 late V/R: ~3.5 mag in ~12 d)
    # is shown as a short arrow from the end of the plateau rather
    # than as individual points.
    v_end = w94[(w94['source'].isin(('IAUC_vis', 'Tsvetkov95')))
                & (w94['band'] == 'V')]
    if len(v_end):
        p_end = v_end.loc[v_end['phase'].idxmax()]
        late = w94[(w94['source'] == 'SCL98')]
        if len(late):
            p_late = late.loc[late['phase'].idxmin()]
            dx = 0.45 * (p_late['phase'] - p_end['phase'])
            # Start the arrow just below the last plateau marker so it
            # does not overlap the marker error bars.  On this inverted
            # axis the tip must stay <= -14.9 to remain visible.
            x0 = p_end['phase'] + 2.0
            y0 = p_end['Mabs'] + 0.30
            end_mabs = min(y0 + 0.9, -15.1)
            dy = end_mabs - y0
            tip = (x0 + dx, end_mabs)
            ax.annotate('', xy=tip, xytext=(x0, y0),
                        arrowprops=dict(arrowstyle='->', color='0.35',
                                        lw=2.0, shrinkA=4))
            arrow_lbl = (x0, y0, tip)
            w94_arrow_y = [y0, end_mabs]
        h = Line2D([], [], color='0.35', lw=2.0,
                   label='post-plateau drop (SCL98)')
        w94_handles.append(h)

    # SN 1998S R band: homogeneous CCD photometry (Fassia et al. 2000,
    # 2001) drawn as a single line, lighter grey than SN 1994W.
    w98_handles = []
    if w98 is not None:
        sub = w98[(w98['source'].isin(('Fassia00', 'Fassia01')))
                  & (w98['band'] == 'R')
                  & (w98['e_mag'] <= ERR_MAX)].sort_values('phase')
        if len(sub):
            h, = ax.plot(sub['phase'], sub['Mabs'], color='0.55',
                         lw=1.8, alpha=0.9, zorder=10,
                         label='Fassia et al. (2000, 2001) [R]')
            w98_handles.append(h)

    ax.axvline(SNB_EXP_PHASE, color='k', ls=':', lw=1.2, alpha=0.6,
               zorder=1, label='SN_B explosion (est.)')

    if args.xlim:
        x_lo, x_hi = args.xlim
    else:
        det_phases = np.concatenate(
            [obs[b]['det']['phase'].values for b in band_style
             if len(obs[b]['det'])])
        x_lo = max(-80.0, det_phases.min() - 5.0)
        x_hi = 150.0
    ax.set_xlim(x_lo, x_hi)

    in_win = np.concatenate(
        [obs[b]['det'].loc[
            (obs[b]['det']['phase'] >= x_lo)
            & (obs[b]['det']['phase'] <= x_hi)
            & (obs[b]['det']['Mabs_err'] <= ERR_MAX), 'Mabs'].values
         for b in band_style]
        + [w94.loc[(w94['phase'] >= x_lo) & (w94['phase'] <= x_hi)
                   & (w94['band'] == 'V'), 'Mabs'].values,
           np.asarray(w94_arrow_y)]
        + ([w98.loc[(w98['phase'] >= x_lo) & (w98['phase'] <= x_hi)
                   & (w98['band'] == 'R'), 'Mabs'].values]
           if w98 is not None else []))
    ax.set_ylim(bottom=-14.9, top=np.nanmin(in_win) - 0.4)

    ax.set_xlabel('Phase (days since peak)')
    ax.set_ylabel('Absolute magnitude')
    ax.set_title('SN 2024pba SN_B observed lightcurve '
                 'vs SN 1994W and SN 1998S')
    ax.grid(True, alpha=0.3)
    # Two legends: SN 2024pba photometry/limits, and a separate
    # SN 1994W legend keyed by data source.
    lim_handle = Line2D([], [], marker='v', color='0.5', ls='',
                        markerfacecolor='none', ms=6,
                        label='5-sigma upper limit')
    handles, labels = ax.get_legend_handles_labels()
    pba = [(h, l) for h, l in zip(handles, labels)
           if l.startswith(('SN 2024pba', 'SN_B'))]
    leg1 = ax.legend([h for h, _ in pba] + [lim_handle],
                     [l for _, l in pba] + ['5-sigma upper limit'],
                     fontsize=9, frameon=False, loc='upper right')
    ax.add_artist(leg1)
    leg94 = ax.legend(handles=w94_handles, title='SN 1994W (V and R)',
                      fontsize=9, frameon=False, loc='lower left')
    ax.add_artist(leg94)
    if w98_handles:
        ax.legend(handles=w98_handles, title='SN 1998S (R)',
                  fontsize=9, frameon=False, loc='lower right')

    fig.tight_layout()

    # Label parallel to the drop arrow.  The angle must be computed
    # after the layout is final so the data->display transform uses the
    # real pixel scaling of the inverted magnitude axis.
    if arrow_lbl is not None:
        fig.canvas.draw()
        (x0, y0), tip = (arrow_lbl[:2], arrow_lbl[2])
        p0 = ax.transData.transform((x0, y0))
        p1 = ax.transData.transform(tip)
        d = p1 - p0
        ang = np.degrees(np.arctan2(d[1], d[0]))
        mid = ((x0 + tip[0]) / 2.0 + 4.5, (y0 + tip[1]) / 2.0)
        ax.text(mid[0], mid[1], 'SN 1994W drop-off', rotation=ang,
                rotation_mode='anchor', ha='center', va='center',
                fontsize=7.5, color='0.35')

    fig.savefig(args.out, dpi=150, bbox_inches='tight')
    print(f'Saved: {args.out}')


if __name__ == '__main__':
    main()
