#!/usr/bin/env python3
"""Compare two STELLA runs of the same model at different frequency grids.

Panels:
  1-3: g, r, o absolute-magnitude lightcurves vs observations.
  4-6: g-r, o-g, o-r colour curves vs observed colour pairs.
  7:   SEDs (L_lambda) at matched epochs for both runs.

Usage:
    python scripts/plot_nfreq_comparison.py \
        --old-dir SN_B_dessart_grid/DX15_V7000_sa350 \
        --new-dir SN_B_dessart_grid/DX15_V7000_sa350_nfreq200 \
        --out SN_B_dessart_grid/DX15_nfreq_comparison.png
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import AutoMinorLocator

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import plot_dessart_grid_vs_obs as base
import verify_stella_photometry as vp


def load_tt_cols(path):
    rows = []
    for line in open(path):
        s = line.strip()
        if not s or s[0].isalpha():
            continue
        try:
            r = [float(x) for x in s.split()]
        except ValueError:
            continue
        if len(r) == 24:
            rows.append(r)
    a = np.array(rows)
    return dict(time=a[:, 0], g=a[:, 13], r=a[:, 14], o=a[:, 21])


def read_flx_safe(path):
    """Read a .flx with its own grid; returns times, nfrus, L_lambda cube."""
    mfreq = vp.detect_mfreq(path)
    freqob = vp.build_freqob(mfreq - 30, mfreq)
    times, nfrus, flx = vp.read_flx(path, mfreq)
    lam_A = 1e8 * vp.CS / (freqob * vp.UFREQ)
    return times, nfrus, flx, freqob, lam_A


def sed_at(times, nfrus, flx, freqob, t_req):
    idx = int(np.argmin(np.abs(times - t_req)))
    lnu = flx[:, idx] * freqob**3 * 10**vp.CLUMNU
    lam_A = 1e8 * vp.CS / (freqob * vp.UFREQ)
    llam = lnu * vp.C_AA / lam_A**2
    llam[nfrus[idx]:] = np.nan
    llam[flx[:, idx] <= 0] = np.nan
    return times[idx], llam


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--old-dir', required=True)
    ap.add_argument('--new-dir', required=True)
    ap.add_argument('--old-label', default='100 freq bins')
    ap.add_argument('--new-label', default='200 freq bins')
    ap.add_argument('--out', default=None)
    ap.add_argument('--sed-epochs', nargs='*', type=float,
                    default=[20.0, 40.0, 60.0])
    args = ap.parse_args()

    name_old = os.path.basename(args.old_dir.rstrip('/'))
    name_new = os.path.basename(args.new_dir.rstrip('/'))
    tt_old = base.find_tt(args.old_dir)
    tt_new = base.find_tt(args.new_dir)
    if tt_new is None:
        sys.exit(f"no .tt yet for {name_new}")

    m_old = load_tt_cols(tt_old)
    m_new = load_tt_cols(tt_new)
    obs = base.load_obs()

    # score both with the shared scorer
    df_old = base.load_tt(tt_old)
    df_new = base.load_tt(tt_new)
    tpk_old = float(df_old.loc[df_old['R'].idxmin(), 'time'])
    tpk_new = float(df_new.loc[df_new['R'].idxmin(), 'time'])
    sc_old = base.best_phase_shift(df_old, obs, tpk_old)
    sc_new = base.best_phase_shift(df_new, obs, tpk_new)
    dt_old, dt_new = sc_old[1], sc_new[1]
    print(f"{name_old}: score={sc_old[0]:.1f} dt={dt_old:+.2f}")
    print(f"{name_new}: score={sc_new[0]:.1f} dt={dt_new:+.2f}")

    fig, axes = plt.subplots(3, 3, figsize=(17, 13), sharex=False)
    plt.subplots_adjust(hspace=0.32, wspace=0.26)

    # .tt columns are native SDSS/ATLAS abmag; ZTF g/r and ATLAS o
    # observations are also AB - like-for-like comparison.
    band_panels = [('g', obs['g'], 'g', 'g band'),
                   ('r', obs['r'], 'r', 'r band'),
                   ('o', obs['o'], 'o', 'o band (ATLAS)')]
    for col, (band, od, mc, band_lab) in enumerate(band_panels):
        ax = axes[0, col]
        ph_old = m_old['time'] - tpk_old + dt_old
        ph_new = m_new['time'] - tpk_new + dt_new
        ax.plot(ph_old, m_old[mc], '-', color='tab:blue', lw=1.4,
                label=args.old_label)
        ax.plot(ph_new, m_new[mc], '-', color='tab:red', lw=1.4,
                label=args.new_label)
        ax.errorbar(od['phase'], od['Mabs'],
                    yerr=od.get('Error', od.get('Mabs_err')),
                    fmt='o', ms=4, color='k', label='Observed')
        ax.invert_yaxis()
        ax.set_xlabel('Days since observed peak')
        ax.set_ylabel('Absolute magnitude')
        ax.set_title(band_lab)
        ax.set_xlim(-50, 140)
        if col == 0:
            ax.legend(fontsize=8)
        ax.xaxis.set_minor_locator(AutoMinorLocator())

    colour_pairs = [('g', 'r', 'g - r'), ('o', 'g', 'o - g'),
                    ('o', 'r', 'o - r')]
    for col, (b1, b2, lab) in enumerate(colour_pairs):
        ax = axes[1, col]
        col_old = m_old[b1] - m_old[b2]
        col_new = m_new[b1] - m_new[b2]
        ax.plot(m_old['time'] - tpk_old + dt_old, col_old, '-',
                color='tab:blue', lw=1.4)
        ax.plot(m_new['time'] - tpk_new + dt_new, col_new, '-',
                color='tab:red', lw=1.4)
        pairs = base.obs_color_pairs(obs, b1, b2,
                                     max_sep=1.5 if 'o' in (b1, b2) else 0.75)
        if len(pairs):
            ax.errorbar(pairs['phase'], pairs['col'], yerr=pairs['err'],
                        fmt='o', ms=4, color='k')
        ax.set_xlabel('Days since observed peak')
        ax.set_ylabel(f'{lab} (mag)')
        ax.set_title(f'{lab} colour')
        ax.set_xlim(-50, 140)
        ax.xaxis.set_minor_locator(AutoMinorLocator())

    # SED panel: merged into bottom row across all three columns
    for c in range(3):
        axes[2, c].remove()
    ax = fig.add_subplot(3, 1, 3)
    flx_old = os.path.join(args.old_dir, 'stella', f'{name_old}.flx')
    flx_new = os.path.join(args.new_dir, 'stella', f'{name_new}.flx')
    if os.path.exists(flx_old) and os.path.exists(flx_new):
        to_, nf_, f_, fo_o, lam_o = read_flx_safe(flx_old)
        tn_, nfn_, fn_, fo_n, lam_n = read_flx_safe(flx_new)
        cmap = plt.cm.viridis
        for k, te in enumerate(args.sed_epochs):
            t_act, ll_o = sed_at(to_, nf_, f_, fo_o, te)
            t_actn, ll_n = sed_at(tn_, nfn_, fn_, fo_n, te)
            col = cmap(k / max(len(args.sed_epochs) - 1, 1))
            ax.plot(lam_o, ll_o, '-', color=col, lw=1.3,
                    label=f'{args.old_label}, {t_act:.0f} d')
            ax.plot(lam_n, ll_n, '--', color=col, lw=1.3,
                    label=f'{args.new_label}, {t_actn:.0f} d')
        ax.set_xlim(2500, 9500)
        ax.set_ylim(1e37, 5e40)
        ax.set_yscale('log')
        ax.set_xlabel('Wavelength (Angstrom)')
        ax.set_ylabel(r'Luminosity per wavelength (erg/s/A)')
        ax.set_title('Spectra at matched epochs')
        ax.legend(fontsize=7, ncol=2)
    else:
        ax.text(0.5, 0.5, 'spectra not available yet',
                transform=ax.transAxes, ha='center')

    fig.suptitle(f'{name_old}: {args.old_label} vs {args.new_label}'
                 f'  (fit scores {sc_old[0]:.1f} / {sc_new[0]:.1f})',
                 fontsize=12)
    out = args.out or os.path.join(args.new_dir, 'nfreq_comparison.png')
    fig.savefig(out, dpi=140)
    print('saved', out)


if __name__ == '__main__':
    main()
