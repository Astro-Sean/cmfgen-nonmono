#!/usr/bin/env python3
"""Per-family model-vs-obs plots for the Dessart-inspired SN_B grid.

One figure per DX family (r, o and g-r panels).  Every model in the
family is drawn with its own best time shift and labelled only by the
name tokens that vary inside the family, so each panel isolates the
parameter that family scans (n for DX9, M/V for DX11, sa for DX12,
and so on).

Usage:
    python scripts/plot_dessart_families_vs_obs.py \
        [--grid-dir SN_B_dessart_grid] [--out-dir family_plots] \
        [--families DX9 DX11] [--fit-end 130] [--ebv 0.15]
"""

import argparse
import os
import re
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import AutoMinorLocator
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import plot_dessart_grid_vs_obs as base


def family_key(name):
    """DX family prefix of a model name (DX7_, DX11_, ... or DX_)."""
    m = re.match(r'(DX\d+_|DX_)', name)
    return m.group(1) if m else None


def family_sort_key(fam):
    m = re.match(r'DX(\d+)_', fam)
    return int(m.group(1)) if m else -1


def varying_label(name, members):
    """Name tokens that differ between the family's members.

    Tokens present in every member (family prefix included) are
    dropped, leaving just the scanned parameter(s).  Falls back to
    the full name for single-member families.
    """
    toks = name.split('_')
    if len(members) == 1:
        return name
    common = set(toks)
    for other in members:
        common &= set(other[0].split('_'))
    keep = [t for t in toks if t not in common]
    return ' '.join(keep) if keep else name


def plot_family(fam, members, obs, args, out_dir):
    scored = []
    for name, df in members:
        tpk = df.loc[df['R'].idxmin(), 'time']
        score, dt, chi_rise, chi_tail, chi_col, pen = (
            base.best_phase_shift(df, obs, tpk, args.fit_end))
        scored.append({
            'name': name, 'df': df, 'tpk': tpk, 'dt': dt,
            'score': score, 'ok': np.isfinite(pen),
            'label': varying_label(name, members),
        })
    scored.sort(key=lambda s: s['score'])

    cmap = plt.cm.viridis
    cols = cmap(np.linspace(0.05, 0.9, max(len(scored), 2)))

    fig, axes = plt.subplots(3, 1, figsize=(9, 9.5), sharex=True,
                             gridspec_kw={'height_ratios': [2, 2, 1.4]})
    axr, axo, axc = axes

    det = obs['r']
    axr.errorbar(det['phase'], det['Mabs'], yerr=det['Error'],
                 fmt='o', color='tab:red', ms=5, mec='k', mew=0.7,
                 capsize=2, elinewidth=0.7, alpha=0.85, zorder=20,
                 label='2024pba r')
    det = obs['o']
    axo.errorbar(det['phase'], det['Mabs'], yerr=det['Mabs_err'],
                 fmt='^', color='darkorange', ms=5, mec='k', mew=0.7,
                 capsize=2, elinewidth=0.7, alpha=0.85, zorder=20,
                 label='2024pba o (ATLAS)')
    pairs = base.obs_color_pairs(obs, 'g', 'r', max_sep=0.75)
    if len(pairs):
        axc.errorbar(pairs['phase'], pairs['col'], yerr=pairs['err'],
                     fmt='o', color='k', ms=5, mec='k', mew=0.7,
                     capsize=2, elinewidth=0.7, alpha=0.85, zorder=20,
                     label='2024pba g - r')

    handles = []
    for i, s in enumerate(scored):
        ph = s['df']['time'] - s['tpk'] + s['dt']
        c = cols[i]
        mark = '' if s['ok'] else '  [foot FAIL]'
        lab = f"{s['label']} ({s['dt']:+.0f} d, {s['score']:.0f}){mark}"
        axr.plot(ph, s['df']['R'], color=c, lw=1.7, zorder=10)
        axo.plot(ph, s['df']['o'], color=c, lw=1.7, zorder=10)
        axc.plot(ph, s['df']['g'] - s['df']['R'], color=c, lw=1.7,
                 zorder=10)
        handles.append(Line2D([], [], color=c, lw=1.7, label=lab))

    for ax in (axr, axo, axc):
        ax.axvline(base.SNB_EXP_PHASE, color='k', ls=':', lw=1.2,
                   alpha=0.6)
        ax.axvline(-19.0, color='0.4', ls=':', lw=1.0, alpha=0.8)
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.grid(True, alpha=0.3)
        ax.grid(True, which='minor', alpha=0.12)
    axr.axvspan(-45.0, -19.0, ymin=0.26, ymax=0.37,
                color='0.85', alpha=0.5, zorder=0)

    axr.set_ylim(top=-20.6, bottom=-13.2)
    axo.set_ylim(top=-20.6, bottom=-13.2)
    if len(pairs):
        c_lo = min(pairs['col'].min(), -0.4) - 0.3
        c_hi = max(pairs['col'].max(), 0.4) + 0.3
        # blue at the top, matching the magnitude panels
        axc.set_ylim(c_hi, c_lo)
    axr.set_ylabel('r-band magnitude')
    axo.set_ylabel('o-band magnitude')
    axc.set_ylabel('g - r (mag)')
    axc.set_xlim(-60.0, 130.0)
    axc.set_xlabel('Days until SN 2024pba peak')
    # family tag inside the panel - a suptitle collides with the
    # legend anchored above the axes
    axr.text(0.02, 0.95, fam.rstrip('_'), transform=axr.transAxes,
             fontsize=11, fontweight='bold', va='top')
    axr.legend(handles=handles, fontsize=7.5, ncol=2, frameon=False,
               loc='lower center', bbox_to_anchor=(0.5, 1.0))

    out = os.path.join(out_dir, f'dessart_{fam.rstrip("_")}_vs_obs.png')
    plt.tight_layout()
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  {fam:7s} {len(scored):2d} models -> {out}')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--grid-dir', default=base.GRID_DIR)
    ap.add_argument('--out-dir', default=None,
                    help='default: <grid-dir>/family_plots')
    ap.add_argument('--families', nargs='*', default=None,
                    help='restrict to these families, e.g. DX9 DX11')
    ap.add_argument('--fit-end', type=float, default=130.0,
                    help='last obs phase (d) included in the best-dt fit')
    ap.add_argument('--ebv', type=float, default=None,
                    help='total E(B-V) incl. host; default keeps '
                         'EBV=0.15 (MW foreground only)')
    args = ap.parse_args()
    if args.ebv is not None:
        base.set_ebv(args.ebv)
        print(f'  E(B-V) total = {base.EBV:.2f} '
              f'(A_g={base.A_G:.2f} A_r={base.A_R:.2f} '
              f'A_o={base.A_O:.2f})')
    out_dir = args.out_dir or os.path.join(args.grid_dir, 'family_plots')
    os.makedirs(out_dir, exist_ok=True)

    families = {}
    import glob
    for mdir in sorted(d for d in
                       glob.glob(os.path.join(args.grid_dir, '*'))
                       if os.path.isdir(d)):
        name = os.path.basename(mdir)
        fam = family_key(name)
        if fam is None:
            continue
        tt = base.find_tt(mdir)
        if tt is None:
            continue
        df = base.load_tt(tt)
        if df is not None and len(df):
            families.setdefault(fam, []).append((name, df))
    if not families:
        raise SystemExit(f'no completed models found in {args.grid_dir}')

    obs = base.load_obs()
    wanted = (None if args.families is None else
              {f if f.endswith('_') else f + '_' for f in args.families})
    for fam in sorted(families, key=family_sort_key):
        if wanted is not None and fam not in wanted:
            continue
        plot_family(fam, families[fam], obs, args, out_dir)


if __name__ == '__main__':
    main()
