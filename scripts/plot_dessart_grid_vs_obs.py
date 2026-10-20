#!/usr/bin/env python3
"""Plot the Dessart-2010jl analog grid vs SN 2024pba observations.

Three band panels sharing a phase axis (days since each model's own
peak): model g/r/o absolute magnitudes vs observed SDSS g/r (PSF)
and ATLAS o.

Usage:
    python scripts/plot_dessart_grid_vs_obs.py \
        [--grid-dir SN_B_dessart_grid] [--out PATH]
"""

import argparse
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import AutoMinorLocator
import numpy as np
import pandas as pd

GRID_DIR = '/ptmp/sbrennan/SN2024pba/SN_B_dessart_grid'
OBS_FILE = '/ptmp/sbrennan/SN2024pba/observational_data/SN_B_lightcurve.dat'
ATLAS_FILE = '/ptmp/sbrennan/SN2024pba/observational_data/SN2024pba_atclean_cleaned.csv'
DL_MPC = 80.4
EBV = 0.15  # MW foreground only (IRSA dust maps: S&F 0.149 at the
            # SN position); host extinction is NOT included - use
            # --ebv to test a total value
PEAK_MJD = 60850.0
SNB_EXP_MJD = 60805.0
SNB_EXP_PHASE = SNB_EXP_MJD - PEAK_MJD
# Extinction coefficients: Gordon et al. 2023 MW curve (Rv=3.1),
# band-averaged over the SVO filter transmissions (g 4719 A,
# r 6185 A, i 7480 A, ATLAS o 6866 A).  These replace the old Schlafly &
# Finkbeiner 2011 values (3.303/2.285); the E(B-V) map value is
# unchanged - only the per-band curve is updated.
A_G = 3.776 * EBV
A_R = 2.695 * EBV
A_I = 2.04 * EBV
A_O = 2.350 * EBV


def set_ebv(ebv):
    """Rescale the band extinctions for a different total E(B-V)."""
    global EBV, A_G, A_R, A_I, A_O
    EBV = ebv
    A_G = 3.776 * EBV
    A_R = 2.695 * EBV
    A_I = 2.04 * EBV
    A_O = 2.350 * EBV


def find_tt(mdir):
    name = os.path.basename(mdir)
    for path in (
            os.path.join(mdir, 'stella', f'{name}.tt'),
            os.path.join(mdir, 'output', f'{name}.tt'),
            os.path.join(mdir, 'stella', 'strad', f'{name}.tt')):
        if os.path.exists(path):
            return path
    return None


def load_tt(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith(('time', 'XNianm')):
                continue
            try:
                row = [float(x) for x in s.split()]
            except ValueError:
                continue
            if len(row) >= 14:
                rows.append(row)
    if not rows:
        return None
    arr = np.array(rows)
    # tt columns: 6=Mbol, 9=MV, 11=MR(Johnson); 13/14/15=Mg/Mr/Mi(SDSS),
    # 20/21=Mc/Mo(ATLAS).  The obs r/g are SDSS-like, so score on Mr.
    df = pd.DataFrame({
        'time': arr[:, 0], 'V': arr[:, 9],
        'R': arr[:, 14], 'g': arr[:, 13], 'i': arr[:, 15],
        'o': arr[:, 21],
    })
    return df[df['time'] >= 0].copy()


def load_obs():
    raw = pd.read_csv(OBS_FILE)
    obs = {}
    for band, A in (('g', A_G), ('r', A_R), ('i', A_I)):
        sub = raw[raw['Filter'] == band].copy()
        det = sub[sub['Mag'].apply(lambda x: x != '-' and not pd.isna(x))]
        det = det.copy()
        det['Mag'] = det['Mag'].astype(float)
        det['Error'] = det['Error'].astype(float)
        det['phase'] = det['MJD'] - PEAK_MJD
        det['Mabs'] = det['Mag'] - 5 * np.log10(DL_MPC * 1e5) - A
        obs[band] = det
        nd = sub[sub['Mag'].apply(lambda x: x == '-')]
        if len(nd):
            nd = nd.copy()
            lim = pd.to_numeric(nd['Limit_5p0S2N'], errors='coerce')
            lim = lim.where(lim.isna(),
                            lim - 5 * np.log10(DL_MPC * 1e5) - A)
            obs.setdefault('lim', {})[band] = pd.DataFrame({
                'phase': nd['MJD'] - PEAK_MJD,
                'Mabs': lim,
            }).dropna()
    atl = pd.read_csv(ATLAS_FILE)
    sub = atl[atl['F'] == 'o'].copy()
    sub['sn'] = sub['uJy'] / sub['duJy_new']
    good = (sub['Mask'].isin(['0x400800', '0x400900'])
            & (sub['chi/N'].isna() | (sub['chi/N'] < 5.0)))
    det = sub[(sub['sn'] >= 3.0) & good].copy()
    det['mag_err'] = det['dm'].where(det['dm'].notna(), 1.086 / det['sn'])
    det['phase'] = det['MJD'] - PEAK_MJD
    det['Mabs'] = det['m'] - 5 * np.log10(DL_MPC * 1e5) - A_O
    det['Mabs_err'] = det['mag_err']
    if len(det):
        det['night'] = np.floor(det['phase'])
        rows = []
        for night, g in det.groupby('night'):
            w = 1.0 / g['Mabs_err']**2
            rows.append({
                'phase': np.average(g['phase'], weights=w),
                'Mabs': np.average(g['Mabs'], weights=w),
                'Mabs_err': 1.0 / np.sqrt(w.sum()),
            })
        det = pd.DataFrame(rows)
    obs['o'] = det
    return obs


def short_label(name):
    return name.replace('DX3_', 'DX3 ').replace('DX2_', 'DX2 ').replace(
        'DX_', '').replace('_', ' ')


# Linestyle per model family so the curve type is readable even among
# the grey background set.
FAMILY_LS = {
    'DX_': (0, (1, 1)),
    'DX2_': (0, (5, 2)),
    'DX3_': (0, (3, 1, 1, 1)),
    'DX4_': '-.',
    'DX5_': (0, (5, 2, 1, 2)),
    'DX6_': (0, (3, 1, 1, 1, 1, 1)),
    'DX7_': '-',
    'DX8_': '--',
    'DX9_': (0, (2, 1)),
    'DX10_': ':',
    'DX11_': (0, (4, 1, 1, 1)),
    'DX12_': (0, (6, 2, 2, 2)),
    'DX13_': (0, (1, 1, 3, 1)),
    'DX14_': (0, (2, 2, 4, 2)),
    'DX15_': (0, (4, 1, 4, 1, 1, 1)),
    'DX16_': (0, (6, 1, 1, 1, 2, 1)),
    'DX17_': (0, (1, 1, 4, 1, 4, 1)),
    'DX18_': (0, (7, 2, 1, 2, 1, 2, 1, 2)),
    'DX19_': (0, (3, 1, 1, 1, 4, 1, 1, 1)),
    'DX20_': (0, (8, 2, 2, 2, 4, 2, 2, 2)),
    'DX21_': (0, (2, 1, 6, 1, 3, 1)),
    'DX22_': (0, (5, 1, 2, 1, 2, 1, 5, 1)),
}


def family_ls(name):
    """Dash pattern keyed on the DX family prefix of the model name."""
    for fam, ls in FAMILY_LS.items():
        if name.startswith(fam):
            return ls
    return '-'


def model_mag_at(df, t, col='R'):
    """Model magnitude at model-clock t; NaN outside coverage.

    A small tolerance clamps queries landing a hair outside the .tt
    time range - a shift anchored so the first recorded epoch sits
    exactly at -19 d must not be discarded by float roundoff.
    """
    t0, t1 = df['time'].iloc[0], df['time'].iloc[-1]
    if t < t0 - 0.02 or t > t1 + 0.02:
        return np.nan
    t = min(max(t, t0), t1)
    return float(np.interp(t, df['time'], df[col]))


def model_R_at(df, t):
    return model_mag_at(df, t, 'R')


# Obs band letter -> model .tt column (SDSS Mg/Mr, ATLAS Mo).
OBS2MOD = {'g': 'g', 'r': 'R', 'i': 'i', 'o': 'o'}


def obs_color_pairs(obs, b1, b2, max_sep=0.75):
    """Observed b1-b2 colours from detections paired within max_sep d.

    One point per b1 epoch, matched to the nearest b2 detection; the
    same pairing rule feeds the score and the colour panels so both
    show the same quantity.
    """
    err1 = 'Error' if b1 in 'gri' else 'Mabs_err'
    err2 = 'Error' if b2 in 'gri' else 'Mabs_err'
    rows = []
    for _, row1 in obs[b1].iterrows():
        ref = obs[b2]
        sep = (ref['phase'] - row1['phase']).abs()
        near = ref[sep <= max_sep]
        if len(near):
            row2 = near.loc[sep[near.index].idxmin()]
            rows.append({
                'phase': row1['phase'],
                'col': row1['Mabs'] - row2['Mabs'],
                'err': np.hypot(row1[err1], row2[err2]),
            })
    return pd.DataFrame(rows)


def best_phase_shift(df, obs, t_pk, fit_end=None, foot_phase=-19.0,
                     foot_slack=10.0, early_phase=80.0, early_wt=5.0,
                     tail_wt=0.5, col_wt=0.0, peak_gate=0.4):
    """Best time offset (days added to the model's own-peak phase).

    The turn-on is anchored to the observed -19 d r-band point by a
    single hard gate: at foot_phase and at every earlier phase the
    model must be fainter than or equal to the point
    (M_model >= M_obs(-19); the model never exceeds its brightness).
    Missing coverage counts as dark.  The foot crossing may therefore
    land at foot_phase or later (a faster, more compressed rise) but
    never earlier - that would put the model above the point at -19 d.
    If no dt is feasible the best-scoring dt overall is kept and the
    model is flagged via pen=inf.

    The model should also stay fainter than the observed plateau
    points at -45..-19 d (soft penalty).  The score also carries
    colour chi2 terms - g-r, o-g and o-r - over the fit window:
    colour is the SED slope, so it discriminates between models that
    match the band lightcurves with different temperatures (the
    known too-blue failure mode).  ATLAS o pairs at 1.5 d cadence vs
    0.75 d for g-r (g only starts at +29 d anyway).
    Returns (score, dt, chi2_rise, chi2_tail, chi2_col, pen).
    """
    # Default: score over every fitted-band detection, so the whole
    # observed tail contributes rather than just the first ~130 d.
    if fit_end is None:
        fit_end = max(obs[b]['phase'].max() for b in ('r', 'o', 'i')
                      if len(obs[b]))
    # ATLAS o carries the rise/peak (dense cadence); r and g carry the
    # tail, which is what each band actually covers.  i joins the early
    # window since its detections are pre-+100 d only.
    plat = obs['o'][(obs['o']['phase'] >= -45) & (obs['o']['phase'] <= -19)]
    fit_r = obs['r'][(obs['r']['phase'] > -19) &
                     (obs['r']['phase'] <= fit_end)]
    fit_g = obs['g'][(obs['g']['phase'] > -19) &
                     (obs['g']['phase'] <= fit_end)]
    fit_o = obs['o'][(obs['o']['phase'] > -19) &
                     (obs['o']['phase'] <= fit_end)]
    fit_i = obs['i'][(obs['i']['phase'] > -19) &
                     (obs['i']['phase'] <= fit_end)]
    fit_cols = []
    for b1, b2, sep in (('g', 'r', 0.75), ('o', 'g', 1.5),
                        ('o', 'r', 1.5), ('r', 'i', 0.75)):
        pairs = obs_color_pairs(obs, b1, b2, max_sep=sep)
        pairs = pairs[(pairs['phase'] > foot_phase) &
                      (pairs['phase'] <= fit_end)]
        fit_cols.append((b1, b2, pairs))

    # Level of the -19 d r point; the model time where its rise first
    # reaches it defines the foot.  A model that never gets that bright
    # is assigned its last time, so it stays faint through the foot and
    # the dark-floor chi2 sinks it anyway.  A model already brighter at
    # its first recorded epoch cannot satisfy the gate under any shift:
    # anchoring its first epoch to foot_phase still leaves it brighter
    # than the point at -19 d.
    # The foot anchor uses the r-band point at -19 d: the model's r rise
    # must cross that magnitude at or after the point.  The peak gate
    # stays on o (densest cadence around maximum).
    r_foot = float(np.interp(foot_phase,
                             obs['r']['phase'], obs['r']['Mabs']))
    rise = df[df['R'] <= r_foot]
    t_cross = (float(rise['time'].iloc[0]) if len(rise)
               else float(df['time'].iloc[-1]))
    dt_anchor = t_pk - t_cross + foot_phase
    # Observed o peak: magnitude and phase.  A model must reproduce it
    # within peak_gate (in either direction) or pay a penalty - total
    # chi2 over the long tail otherwise lets a model win while missing
    # the peak and pre-80 d decline entirely.
    i_pk = fit_o['Mabs'].idxmin()
    obs_r_pk = float(fit_o.loc[i_pk, 'Mabs'])
    pk_ph = float(fit_o.loc[i_pk, 'phase'])
    # Phases up to and including the anchor point itself.
    early_ph = np.arange(-60.0, foot_phase + 1e-9, 0.1)

    def violates_foot(dt):
        for ph in early_ph:
            rm = model_mag_at(df, ph - dt + t_pk, 'R')
            if np.isfinite(rm) and rm < r_foot - 0.05:
                return True  # brighter than the point at or before -19 d
        return False

    best = (np.inf, 0.0, np.inf, np.inf, np.inf, np.inf)
    fallback = best
    for dt in np.arange(dt_anchor - 1.0,
                        dt_anchor + foot_slack + 1e-9, 0.25):
        pen = 0.0
        for _, row in plat.iterrows():
            rm = model_mag_at(df, row['phase'] - dt + t_pk, 'o')
            if np.isfinite(rm):
                excess = row['Mabs'] - rm
                if excess > 0:
                    pen += (excess / 0.15) ** 2
        rm_pk = model_mag_at(df, pk_ph - dt + t_pk, 'o')
        if not np.isfinite(rm_pk):
            pen += (obs_r_pk - (-12.0)) ** 2
        else:
            d_pk = abs(rm_pk - obs_r_pk)
            if d_pk > peak_gate:
                pen += ((d_pk - peak_gate) / 0.15) ** 2
        # Total (not mean) chi2 so large shifts can't win by shrinking
        # coverage; a model absent at a detection is effectively dark
        # (Rm = -12 floor) and underpredicts it.  Detections up to
        # early_phase are upweighted so the rise/peak/early-decline
        # fidelity is a prerequisite rather than a minority term among
        # the much more numerous tail points.
        chi2_rise = 0.0
        chi2_tail = 0.0
        # Every band contributes to the early window (the dense ATLAS o
        # cadence dominates it anyway); past early_phase only r and g
        # are scored - o coverage ends near +250 d and the systematic
        # o-band offset in the models would otherwise contaminate the
        # tail fit.
        for det, err_col, mcol in ((fit_r, 'Error', 'R'),
                                   (fit_g, 'Error', 'g'),
                                   (fit_o, 'Mabs_err', 'o'),
                                   (fit_i, 'Error', 'i')):
            for _, row in det.iterrows():
                if row['phase'] > early_phase and mcol not in ('R', 'g'):
                    continue
                rm = model_mag_at(df, row['phase'] - dt + t_pk, mcol)
                if not np.isfinite(rm):
                    rm = -12.0
                err = np.hypot(row[err_col], 0.25)
                c2 = ((rm - row['Mabs']) / err) ** 2
                if row['phase'] <= early_phase:
                    chi2_rise += early_wt * c2
                else:
                    chi2_tail += tail_wt * c2
        # Colour chi2 over every paired epoch; 0.15 mag floor on the
        # colour error, missing coverage counted as a 0.75 mag residual
        # so a model that ends early is penalised rather than excused.
        chi2_col = 0.0
        for b1, b2, pairs in fit_cols:
            m1, m2 = OBS2MOD[b1], OBS2MOD[b2]
            for _, row in pairs.iterrows():
                t_m = row['phase'] - dt + t_pk
                resid = (model_mag_at(df, t_m, m1)
                         - model_mag_at(df, t_m, m2)) - row['col']
                if not np.isfinite(resid):
                    resid = 0.75
                chi2_col += col_wt * (
                    resid / np.hypot(row['err'], 0.15)) ** 2
        score = chi2_rise + chi2_tail + chi2_col + pen
        if score < fallback[0]:
            fallback = (score, dt, chi2_rise, chi2_tail, chi2_col, np.inf)
        if violates_foot(dt):
            continue
        if score < best[0]:
            best = (score, dt, chi2_rise, chi2_tail, chi2_col, pen)
    if not np.isfinite(best[0]):
        best = fallback  # no dt satisfies the foot; flag via pen=inf
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--grid-dir', default=GRID_DIR)
    ap.add_argument('--out', default=None)
    ap.add_argument('--fit-end', type=float, default=None,
                    help='last obs phase (d) in the best-dt fit; '
                         'default uses the latest fitted-band detection')
    ap.add_argument('--ebv', type=float, default=None,
                    help='total E(B-V) incl. host; default keeps EBV=0.15 '
                         '(MW foreground only)')
    args = ap.parse_args()
    if args.ebv is not None:
        set_ebv(args.ebv)
        print(f'  E(B-V) total = {EBV:.2f} '
              f'(A_g={A_G:.2f} A_r={A_R:.2f} A_o={A_O:.2f})')
    out_file = args.out or os.path.join(
        args.grid_dir, 'dessart_grid_vs_obs.png')

    models = []
    for mdir in sorted(d for d in
                       __import__('glob').glob(os.path.join(args.grid_dir, '*'))
                       if os.path.isdir(d)):
        tt = find_tt(mdir)
        if tt is None:
            continue
        df = load_tt(tt)
        if df is not None and len(df):
            models.append((os.path.basename(mdir), df))
    if not models:
        raise SystemExit(f'no completed models found in {args.grid_dir}')
    print(f'loaded {len(models)} models')

    obs = load_obs()

    colors = plt.cm.tab10(np.arange(10))
    x_lo, x_hi = -60.0, 130.0

    # --- three band panels ---
    BANDS = [
        ('g', 'g', 'Error', 'o', 'tab:green', '2024pba g'),
        ('r', 'R', 'Error', 'o', 'tab:red', '2024pba r'),
        ('o', 'o', 'Mabs_err', '^', 'darkorange', '2024pba o (ATLAS)'),
    ]
    fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True,
                             gridspec_kw={'height_ratios': [2, 2, 2]})
    axg, axr, axo = axes
    band_ax = {'g': axg, 'r': axr, 'o': axo}
    band_modcol = {'g': 'g', 'r': 'R', 'o': 'o'}

    for band, mcol, ecol, fmt, col, lab in BANDS:
        ax = band_ax[band]
        det = obs[band]
        ax.errorbar(det['phase'], det['Mabs'], yerr=det[ecol],
                    fmt=fmt, color=col, ms=5, mec='k', mew=0.7,
                    capsize=2, elinewidth=0.7, alpha=0.85, zorder=20,
                    label=lab if band != 'g' else None)
        ax.axvline(SNB_EXP_PHASE, color='k', ls=':', lw=1.2, alpha=0.6)
        ax.set_ylim(top=-20.6, bottom=-13.2)
        ax.set_ylabel(f'{band}-band magnitude')
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.grid(True, alpha=0.3)
        ax.grid(True, which='minor', alpha=0.12)
    # Early r-band plateau: at -19 d the model may be at the -19 d
    # point's magnitude or fainter, never brighter (see the foot
    # anchor in best_phase_shift)
    axr.axvspan(-45.0, -19.0, ymin=0.26, ymax=0.37,
                color='0.85', alpha=0.5, zorder=0)
    axr.axvline(-19.0, color='0.4', ls=':', lw=1.0, alpha=0.8)

    r_foot = float(np.interp(-19.0, obs['r']['phase'],
                             obs['r']['Mabs']))

    scored = []
    for name, df in models:
        tpk = df.loc[df['R'].idxmin(), 'time']
        score, dt, chi_rise, chi_tail, chi_gr, pen = best_phase_shift(
            df, obs, tpk, args.fit_end)
        r19 = model_R_at(df, -19.0 - dt + tpk)
        # pen=inf means no dt satisfied the -19 d foot gates: the
        # model either stays fainter than the observed point at -19 d
        # (rise too slow / no coverage) or is brighter than it at some
        # earlier phase.
        ok = np.isfinite(pen)
        scored.append(
            (score, dt, ok, r19, chi_rise, chi_tail, chi_gr, name, df))
        print(f'  {name:36s} dt={dt:+6.2f}  score={score:8.1f}  '
              f'(rise {chi_rise:7.1f} tail {chi_tail:8.1f} '
              f'col {chi_gr:7.1f})  '
              f'R(-19d)={r19:6.2f} vs {r_foot:5.2f}  '
              f'pen={pen:8.1f}  {"OK" if ok else "VIOLATES"}')

    # plateau constraint is a gate: violators rank below all passers
    scored.sort(key=lambda s: (not s[2], s[0]))
    n_top = 5
    print('  --- top 5:', ', '.join(s[7] for s in scored[:n_top]))
    for rank, (score, dt, ok, r19, chi_rise, chi_tail, chi_gr, name,
               df) in enumerate(scored):
        tpk = df.loc[df['R'].idxmin(), 'time']
        ph = df['time'] - tpk + dt
        ls = family_ls(name)
        if rank < n_top:
            c = colors[rank]
            lab = f'{short_label(name)} ({dt:+.0f}d)'
            for band, mcol, ecol, fmt, col, olab in BANDS:
                band_ax[band].plot(
                    ph, df[mcol], color=c, lw=1.8, ls=ls, zorder=10,
                    label=lab if band == 'g' else None)
        else:
            for band in band_ax:
                band_ax[band].plot(ph, df[band_modcol[band]], color='0.6',
                                   lw=0.9, ls=ls, alpha=0.45, zorder=2)

    obs_handles = [
        Line2D([], [], marker='o', color='tab:green', ls='', ms=5,
               mec='k', mew=0.7, label='2024pba g'),
        Line2D([], [], marker='o', color='tab:red', ls='', ms=5,
               mec='k', mew=0.7, label='2024pba r'),
        Line2D([], [], marker='^', color='darkorange', ls='', ms=5,
               mec='k', mew=0.7, label='2024pba o (ATLAS)'),
    ]
    h, l = axg.get_legend_handles_labels()
    axg.legend(handles=h + obs_handles, fontsize=8, ncol=4,
               frameon=False, loc='lower center',
               bbox_to_anchor=(0.5, 1.0))

    fam_seen = [fam for fam in FAMILY_LS
                if any(s[7].startswith(fam) for s in scored)]
    fam_handles = [
        Line2D([], [], color='0.3', lw=1.4, ls=FAMILY_LS[fam],
               label=fam.rstrip('_'))
        for fam in fam_seen]
    if fam_handles:
        axo.legend(handles=fam_handles, fontsize=7, ncol=2,
                   frameon=False, loc='lower left',
                   title='model family', title_fontsize=7)


    axo.set_xlim(x_lo, x_hi)
    axo.set_xlabel('Days until SN 2024pba peak')

    plt.tight_layout()
    plt.savefig(out_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Saved: {out_file}')

    # --- colour plot: g-r, o-g, o-r panels (same shifts/styles) ---
    # Same pairing as the score: 0.75 d for g-r, 1.5 d for the
    # sparser ATLAS o epochs.
    colour_defs = [('g', 'r', 0.75), ('o', 'g', 1.5), ('o', 'r', 1.5)]
    figc, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True)
    for ax, (b1, b2, sep) in zip(axes, colour_defs):
        pairs = obs_color_pairs(obs, b1, b2, max_sep=sep)
        if len(pairs):
            ax.errorbar(pairs['phase'], pairs['col'], yerr=pairs['err'],
                        fmt='o', color='k', ms=5, mec='k', mew=0.7,
                        capsize=2, elinewidth=0.7, alpha=0.85,
                        zorder=20, label=f'2024pba {b1} - {b2}')
        for rank, (score, dt, ok, r19, chi_rise, chi_tail, chi_gr,
                   name, df) in enumerate(scored):
            tpk = df.loc[df['R'].idxmin(), 'time']
            ph = df['time'] - tpk + dt
            mcol = df[OBS2MOD[b1]] - df[OBS2MOD[b2]]
            if rank < n_top:
                ax.plot(ph, mcol, color=colors[rank], lw=1.8,
                        ls=family_ls(name), zorder=10,
                        label=(f'{short_label(name)} ({dt:+.0f}d)'
                               if b1 == 'g' else None))
            else:
                ax.plot(ph, mcol, color='0.6', lw=0.9,
                        ls=family_ls(name), alpha=0.45, zorder=2)
        ax.axvline(SNB_EXP_PHASE, color='k', ls=':', lw=1.2, alpha=0.6)
        ax.axvline(-19.0, color='0.4', ls=':', lw=1.0, alpha=0.8)
        ax.set_ylabel(f'{b1} - {b2} (mag)')
        ax.grid(True, alpha=0.3)
        ax.grid(True, which='minor', alpha=0.12)
        if len(pairs):
            c_lo = min(pairs['col'].min(), -0.4) - 0.3
            c_hi = max(pairs['col'].max(), 0.4) + 0.3
            # blue (smaller colour) at the top, matching the
            # magnitude panels' inverted axis convention
            ax.set_ylim(c_hi, c_lo)
        ax.legend(fontsize=7, ncol=4, frameon=False,
                  loc='upper left', bbox_to_anchor=(0.0, 1.0))
    axes[-1].set_xlim(x_lo, x_hi)
    axes[-1].set_xlabel('Days until SN 2024pba peak')
    plt.tight_layout()
    gr_file = out_file.replace('.png', '_gr.png')
    plt.savefig(gr_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Saved: {gr_file}')


if __name__ == '__main__':
    main()
