#!/usr/bin/env python3
"""Plot the grid's best COLOUR-matching models vs SN 2024pba.

Ranks every completed model by chi2_col alone (g-r + o-g + o-r) and
colours the top --n-top of them; the rest stay grey.  The r-band
panel sits on top so the lightcurve cost of a good colour match is
visible in the same figure.

Usage:
    python scripts/plot_dessart_grid_best_color.py [--n-top 5]
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import AutoMinorLocator
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import plot_dessart_grid_vs_obs as base


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--grid-dir', default=base.GRID_DIR)
    ap.add_argument('--out', default=None)
    ap.add_argument('--fit-end', type=float, default=130.0)
    ap.add_argument('--ebv', type=float, default=None,
                    help='total E(B-V) incl. host; default keeps '
                         'EBV=0.15 (MW foreground only)')
    ap.add_argument('--n-top', type=int, default=5)
    args = ap.parse_args()
    if args.ebv is not None:
        base.set_ebv(args.ebv)
        print(f'  E(B-V) total = {base.EBV:.2f} '
              f'(A_g={base.A_G:.2f} A_r={base.A_R:.2f} '
              f'A_o={base.A_O:.2f})')
    out_file = args.out or os.path.join(
        args.grid_dir, 'dessart_grid_vs_obs_colorbest.png')

    import glob
    models = []
    for mdir in sorted(d for d in
                       glob.glob(os.path.join(args.grid_dir, '*'))
                       if os.path.isdir(d)):
        tt = base.find_tt(mdir)
        if tt is None:
            continue
        df = base.load_tt(tt)
        if df is not None and len(df):
            models.append((os.path.basename(mdir), df))
    if not models:
        raise SystemExit(f'no completed models in {args.grid_dir}')

    obs = base.load_obs()
    scored = []
    for name, df in models:
        tpk = df.loc[df['R'].idxmin(), 'time']
        score, dt, cr, ct, cc, pen = base.best_phase_shift(
            df, obs, tpk, args.fit_end)
        scored.append({
            'name': name, 'df': df, 'tpk': tpk, 'dt': dt,
            'score': score, 'chi_col': cc,
            'ok': np.isfinite(pen),
        })
    # rank by the colour term only
    scored.sort(key=lambda s: (not s['ok'], s['chi_col']))
    n_top = args.n_top
    print(f'loaded {len(scored)} models')
    print('  --- best colour matches:',
          ', '.join(s['name'] for s in scored[:n_top]))

    colors = plt.cm.tab10(np.arange(10))
    fig, axes = plt.subplots(4, 1, figsize=(10, 12), sharex=True,
                             gridspec_kw={
                                 'height_ratios': [2, 1.4, 1.4, 1.4]})
    axr, ax_gr, ax_og, ax_or = axes
    colour_ax = {('g', 'r'): ax_gr, ('o', 'g'): ax_og,
                 ('o', 'r'): ax_or}

    det = obs['r']
    axr.errorbar(det['phase'], det['Mabs'], yerr=det['Error'],
                 fmt='o', color='tab:red', ms=5, mec='k', mew=0.7,
                 capsize=2, elinewidth=0.7, alpha=0.85, zorder=20,
                 label='2024pba r')
    axr.set_ylim(top=-20.6, bottom=-13.2)
    axr.set_ylabel('r-band magnitude')
    axr.axvspan(-45.0, -19.0, ymin=0.26, ymax=0.37,
                color='0.85', alpha=0.5, zorder=0)

    for (b1, b2), ax in colour_ax.items():
        sep = 0.75 if b1 == 'g' else 1.5
        pairs = base.obs_color_pairs(obs, b1, b2, max_sep=sep)
        if len(pairs):
            ax.errorbar(pairs['phase'], pairs['col'],
                        yerr=pairs['err'], fmt='o', color='k', ms=5,
                        mec='k', mew=0.7, capsize=2, elinewidth=0.7,
                        alpha=0.85, zorder=20,
                        label=f'2024pba {b1} - {b2}')
            c_lo = min(pairs['col'].min(), -0.4) - 0.3
            c_hi = max(pairs['col'].max(), 0.4) + 0.3
            ax.set_ylim(c_hi, c_lo)  # blue at the top
        ax.set_ylabel(f'{b1} - {b2} (mag)')

    for ax in axes:
        ax.axvline(base.SNB_EXP_PHASE, color='k', ls=':', lw=1.2,
                   alpha=0.6)
        ax.axvline(-19.0, color='0.4', ls=':', lw=1.0, alpha=0.8)
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.grid(True, alpha=0.3)
        ax.grid(True, which='minor', alpha=0.12)

    handles = []
    for rank, s in enumerate(scored):
        df, tpk, dt = s['df'], s['tpk'], s['dt']
        ph = df['time'] - tpk + dt
        ls = base.family_ls(s['name'])
        if rank < n_top:
            c = colors[rank]
            lab = (f"{base.short_label(s['name'])} "
                   f"({dt:+.0f} d, col {s['chi_col']:.0f})")
            axr.plot(ph, df['R'], color=c, lw=1.8, ls=ls, zorder=10)
            for (b1, b2), ax in colour_ax.items():
                ax.plot(ph, df[base.OBS2MOD[b1]] - df[base.OBS2MOD[b2]],
                        color=c, lw=1.8, ls=ls, zorder=10)
            handles.append(Line2D([], [], color=c, lw=1.8, ls=ls,
                                  label=lab))
        else:
            axr.plot(ph, df['R'], color='0.6', lw=0.9, ls=ls,
                     alpha=0.45, zorder=2)
            for (b1, b2), ax in colour_ax.items():
                ax.plot(ph,
                        df[base.OBS2MOD[b1]] - df[base.OBS2MOD[b2]],
                        color='0.6', lw=0.9, ls=ls, alpha=0.45,
                        zorder=2)

    axr.legend(handles=handles, fontsize=8, ncol=2, frameon=False,
               loc='lower center', bbox_to_anchor=(0.5, 1.0),
               title='best colour matches', title_fontsize=8)
    for ax in colour_ax.values():
        ax.legend(fontsize=7, ncol=4, frameon=False, loc='upper left',
                  bbox_to_anchor=(0.0, 1.0))
    ax_or.set_xlim(-60.0, 130.0)
    ax_or.set_xlabel('Days until SN 2024pba peak')

    plt.tight_layout()
    plt.savefig(out_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Saved: {out_file}')


if __name__ == '__main__':
    main()
