#!/usr/bin/env python3
"""Plot SN_A_pba grid STELLA lightcurves vs SN 2024pba first-peak observations.

The STELLA model lightcurves give absolute magnitudes vs days since
explosion (t=0).  The only source used is the STELLA .tt file.
New-build .tt files (24 cols) carry real SDSS ugriz in AB and those
are used directly; older 14-col files only carry Johnson UBVRI and
fall back to the V ~ g, R ~ r, I ~ i proxy (~0.1-0.3 mag band-pass
systematic).  The *_lightcurve*.csv products are never read:
the 'complete' csv is a pysynphot re-extraction from .flx spectra,
which is a different and partly broken photometric reduction.

The observations (observational_data/SN_A_lightcurve.dat) are apparent
SDSS g/r/i magnitudes vs MJD, with non-detections as upper limits.

Conventions:
  - Redshift z = 0.0186 (TNS UCSC classification report, 2025TNSCR2900)
  - Distance modulus DM = 34.52  (~80 Mpc, H0=70, Omega_m=0.3, Omega_L=0.7)
  - Galactic reddening E(B-V) = 0.15 (IRSA S&F2011, RA=25.3077 Dec=49.0635)
    A_g = 3.303*E(B-V) = 0.495,  A_r = 2.285*E(B-V) = 0.343,
    A_i = 1.698*E(B-V) = 0.255
  - Explosion epoch is bracketed by the i-band: non-detection at
    MJD 60487.4 (i > 19.8), first detection at 60490.4 (i = 19.7).
    Default EXP_MJD = 60488.0 (--exp-mjd to change; --align-peak instead
    shifts each model so its r-band peak falls on the observed r peak).

Each model is scored by RMS magnitude deviation against the g/r/i
detections at MJD < 60544 (upper limits excluded; later epochs belong
to the second event and are not part of the SN_A fit).

By default lightcurves are peak-aligned: the x-axis is days since the
i-band peak, with each model's I/i peak placed at zero together with the
observed i peak (MJD 60507.493).  Use --no-align-peak to switch to
explosion-date alignment (phase zero = --exp-mjd, default MJD 60488,
the i-band explosion bracket).

Usage:
    python scripts/plot_sn_a_pba_grid_vs_obs.py [--no-align-peak]
                                              [--exp-mjd 60488.0]
                                              [--out FIG.png]
"""
import argparse
import glob
import os
import re

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE = '/viper/ptmp1/sbrennan/SN2024pba'
GRID_DIR = os.path.join(BASE, 'SN_A_pba_grid')
OBS_FILE = os.path.join(BASE, 'observational_data', 'SN_A_lightcurve.dat')
OUT_FILE = os.path.join(GRID_DIR, 'sn_a_pba_grid_vs_obs.png')

Z = 0.0186
DL_MPC = 80.4
DM = 5.0 * np.log10(DL_MPC * 1e6 / 10.0)          # ~34.52
EBV = 0.15
A_G = 3.303 * EBV                                  # ~0.495
A_R = 2.285 * EBV                                  # ~0.343
A_I = 1.698 * EBV                                  # ~0.255
EXTINCTION = {'g': A_G, 'r': A_R, 'i': A_I}

# Explosion epoch bracketed by i-band: non-det MJD 60487.389 (i > 19.8),
# first det MJD 60490.428 (i = 19.68).  Midpoint ~60488.
EXP_MJD = 60488.0
# Observed i-band peak: brightest secure i detection, i=17.586+-0.038 at
# MJD 60507.493.  (The i=17.5 point at MJD 60698 is formally brighter but
# is a late-time subtraction artifact.)  Used as the peak-alignment
# reference (x = 0).
OBS_PEAK_MJD = 60507.493
PEAK_BAND = 'i'

# Best-fit judgements for SN_A use only detections before this MJD.
# Later photometry belongs to the second-event / SN_B epoch and must
# not count against an SN_A-only model.
FIT_MAX_MJD = 60544.0

# superbol (mnicholl/superbol) band colours for the observed points
BAND_COLORS = {'g': 'g', 'r': 'r', 'i': 'goldenrod'}

# Focused-fit weights:
# - the rise is anchored in the i band: the two non-detections before the
#   first i detection are one-sided upper limits, the first i detection is
#   matched directly.  Each i constraint weighs W_RISE.
# - detections at phase <= EARLY_PHASE_CUT (days since i peak; covers the
#   rise, the peak and the first two weeks post-peak) weigh W_EARLY,
#   later detections weigh W_LATE.
EARLY_PHASE_CUT = 15.0
W_EARLY = 1.0
W_LATE = 0.25
W_RISE = 2.0
N_I_NONDET = 2
DARK_MAG = -12.0   # model magnitude when the phase is outside coverage


def parse_grid_name(name):
    """'E05_Ni0p075' -> (5.0, 0.075, ''); envelope variants keep a suffix
    ('E05_Ni0p04_M0p05' -> (5.0, 0.04, '_M0p05'))."""
    m = re.match(r'E(\d+(?:p\d+)?)_Ni(\d+(?:p\d+)?)(_.*)?$', name)
    if not m:
        raise ValueError(name)
    return (float(m.group(1).replace('p', '.')),
            float(m.group(2).replace('p', '.')),
            m.group(3) or '')


def find_tt(mdir):
    """Locate the model's STELLA .tt lightcurve file.

    Search order: canonical stella/ output, the output/ copy, the
    in-progress strad/ run dir, then any stella/<name>*.tt variant
    (e.g. _sigmoid).  backup_* dirs are excluded - stale data.
    """
    name = os.path.basename(mdir)
    candidates = [
        os.path.join(mdir, 'stella', f'{name}.tt'),
        os.path.join(mdir, 'output', f'{name}.tt'),
        os.path.join(mdir, 'stella', 'strad', f'{name}.tt'),
    ]
    candidates += sorted(
        glob.glob(os.path.join(mdir, 'stella', f'{name}*.tt')))
    # Fallback for symlinked/renamed model dirs where the .tt keeps the
    # original basename: any lightcurve under stella/ or output/.
    # Skip debug runs and backup dirs (stale data).
    candidates += sorted(
        p for p in glob.glob(os.path.join(mdir, 'stella', '*.tt'))
        + glob.glob(os.path.join(mdir, 'output', '*.tt'))
        if 'debug' not in os.path.basename(p))
    for path in candidates:
        if os.path.exists(path) and 'backup' not in path:
            return path
    return None


def load_model(mdir):
    """Return (df, source) for a model dir, or (None, None).

    Always reads the STELLA .tt file.  New-build .tt files (>=24 cols)
    carry real SDSS ugriz in AB (cols 12-16) and those are used for
    g/r/i; old 14-col files fall back to the Johnson V/R/I proxy.
    source is 'tt_sdss' or 'tt' accordingly.  The *_lightcurve*.csv
    products are never used - the 'complete' csv is a pysynphot
    re-extraction from .flx spectra, which is a different (and partly
    broken) photometric reduction.
    """
    tt = find_tt(mdir)
    if tt is None:
        return None, None
    rows = []
    with open(tt) as fh:
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
        return None, None
    arr = np.array(rows)
    df = pd.DataFrame({
        'time': arr[:, 0], 'U': arr[:, 7], 'B': arr[:, 8],
        'V': arr[:, 9], 'I': arr[:, 10], 'R': arr[:, 11],
        'bol_mag': arr[:, 6]})
    df = df[df['time'] >= 0].copy()
    if arr.shape[1] >= 24:
        # New STELLA .tt: real SDSS ugriz in AB at cols 12-16
        # (Mu Mg Mr Mi Mz).  Prefer them over the Johnson proxy.
        # Cols 20/21 are ATLAS c/o (AB); exposed for comparison with
        # ATLAS forced photometry, not part of the g/r/i score.
        df['g'] = arr[df.index, 13]
        df['r'] = arr[df.index, 14]
        df['i'] = arr[df.index, 15]
        df['c'] = arr[df.index, 20]
        df['o'] = arr[df.index, 21]
        src = 'tt_sdss'
    else:
        # Old 14-col .tt: Johnson only - V~g, R~r, I~i proxy
        df['g'] = df['V']
        df['r'] = df['R']
        df['i'] = df['I']
        src = 'tt'
    return df, src


def discover_models(grid_dir):
    models = []
    for mdir in sorted(glob.glob(os.path.join(grid_dir, 'E*_Ni*'))):
        name = os.path.basename(mdir)
        if not os.path.isdir(mdir):
            continue
        try:
            E, Ni, suffix = parse_grid_name(name)
        except ValueError:
            continue
        df, src = load_model(mdir)
        if df is not None and len(df):
            models.append({'name': name, 'E': E, 'Ni': Ni,
                           'suffix': suffix, 'df': df, 'src': src})
    models.sort(key=lambda m: (m['E'], m['Ni'], m['suffix']))
    return models


def load_observations(path):
    raw = pd.read_csv(path)
    obs = {}
    for band, A in EXTINCTION.items():
        sub = raw[raw['Filter'] == band].copy()
        det = sub[sub['Mag'].apply(lambda x: x != '-' and not pd.isna(x))].copy()
        det['Mag'] = det['Mag'].astype(float)
        det['Error'] = det['Error'].astype(float)
        det['Mabs'] = det['Mag'] - DM - A
        det['Mabs_err'] = det['Error']
        nd = sub[sub['Mag'].apply(lambda x: x == '-' or pd.isna(x))].copy()
        if len(nd):
            lim5 = pd.to_numeric(nd['Limit_5p0S2N'], errors='coerce')
            nd['Mabs_lim'] = lim5 - DM - A
        obs[band] = {'det': det, 'nd': nd}
    return obs


def score_model(df, obs, shift):
    """Focused weighted-RMS fit vs SN 2024pba first peak.

    Three ingredients, all in magnitudes:
      1. i-band rise gate: at the last N_I_NONDET i non-detections before
         the first i detection the model must be fainter than the 5-sigma
         limit (one-sided residual), and the first i detection is matched
         directly.  Each constraint carries weight W_RISE.
      2. g/r/i detections at phase <= EARLY_PHASE_CUT (rise, peak and
         <15 d post-peak) weigh W_EARLY; later detections up to
         FIT_MAX_MJD weigh W_LATE.
      3. combined = weighted RMS over all residuals.

    model time + shift = phase (days since observed i peak).  Phases
    outside the model coverage count as DARK_MAG for the first i
    detection (a model that has no flux that early is dark and satisfies
    a non-detection by construction).
    """
    scores = {}
    num = 0.0
    den = 0.0
    for band in ('g', 'r', 'i'):
        det = obs[band]['det']
        t_mod = df['time'].values + shift
        m_mod = df[band].values
        ok = (det['phase'].values >= t_mod.min()) & \
             (det['phase'].values <= t_mod.max()) & \
             (det['MJD'].values < FIT_MAX_MJD)
        if ok.sum() < 3:
            scores[band] = np.nan
            continue
        phases = det['phase'].values[ok]
        interp = np.interp(phases, t_mod, m_mod)
        resid = interp - det['Mabs'].values[ok]
        w = np.where(phases <= EARLY_PHASE_CUT, W_EARLY, W_LATE)
        scores[band] = float(np.sqrt(np.sum(w * resid**2) / np.sum(w)))
        scores[f'n_{band}'] = int(ok.sum())
        num += float(np.sum(w * resid**2))
        den += float(np.sum(w))

    # i-band rise gate: non-detection limits + first detection
    iobs = obs['i']
    det_i = iobs['det'].sort_values('phase')
    nd_i = iobs['nd'].sort_values('phase')
    t_i = df['time'].values + shift
    m_i = df['i'].values
    rise_pen = 0.0
    n_rise = 0
    if len(det_i):
        first_ph = det_i['phase'].values[0]
        nd_use = nd_i[nd_i['phase'] < first_ph].tail(N_I_NONDET)
        for ph, mlim in zip(nd_use['phase'].values,
                            nd_use['Mabs_lim'].values):
            if ph > t_i.max():
                continue
            # phases before model t=0 are dark -> satisfy the limit
            m_mod = (np.interp(ph, t_i, m_i) if ph >= t_i.min()
                     else np.inf)
            # one-sided: only penalise a model brighter than the limit
            pen = max(0.0, mlim - m_mod)
            rise_pen += pen * pen
            n_rise += 1
        # first i detection: direct match (dark if outside coverage)
        m_det = det_i['Mabs'].values[0]
        if first_ph <= t_i.max():
            m_mod = (np.interp(first_ph, t_i, m_i)
                     if first_ph >= t_i.min() else DARK_MAG)
            pen = m_mod - m_det
            rise_pen += pen * pen
            n_rise += 1
    scores['rise_rms'] = float(np.sqrt(rise_pen / n_rise)) if n_rise \
        else np.nan
    num += W_RISE * rise_pen
    den += W_RISE * n_rise

    parts = [s for s in (scores.get('g'), scores.get('r'), scores.get('i'))
             if s is not None and not np.isnan(s)]
    scores['combined'] = float(np.sqrt(num / den)) \
        if den > 0 and parts else np.nan
    return scores


def score_model_shift(df, obs, base_shift, dmax=5.0, dstep=0.25):
    """Scan a small extra time offset around base_shift and return the
    best (scores, delta).  The observed i peak is only bracketed to a
    few days by sparse sampling, so a bounded shift is physical - it is
    NOT a free fit parameter; dmax keeps it honest."""
    best_sc, best_d = None, 0.0
    for d in np.arange(-dmax, dmax + 0.5 * dstep, dstep):
        sc = score_model(df, obs, base_shift + d)
        c = sc['combined']
        if not np.isnan(c) and (best_sc is None or c < best_sc['combined']):
            best_sc, best_d = sc, float(d)
    return best_sc, best_d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--exp-mjd', type=float, default=EXP_MJD,
                    help='Adopted explosion MJD for the observed phase axis '
                         f'when --no-align-peak is used (default {EXP_MJD})')
    ap.add_argument('--align-peak', action=argparse.BooleanOptionalAction,
                    default=True,
                    help='Align by peak time (default): x-axis is days since '
                         'the r-band peak, model and observed peaks at 0. '
                         'Use --no-align-peak for explosion-date alignment.')
    ap.add_argument('--out', default=OUT_FILE)
    args = ap.parse_args()

    obs = load_observations(OBS_FILE)
    for band in ('g', 'r', 'i'):
        for k in ('det', 'nd'):
            obs[band][k] = obs[band][k].copy()
            if args.align_peak:
                # Days since the observed i-band peak (peak = 0)
                obs[band][k]['phase'] = obs[band][k]['MJD'] - OBS_PEAK_MJD
            else:
                obs[band][k]['phase'] = obs[band][k]['MJD'] - args.exp_mjd

    models = discover_models(GRID_DIR)
    if not models:
        print(f'No completed models found in {GRID_DIR}')
        return
    ref_desc = (f'days since {PEAK_BAND}-band peak' if args.align_peak
                else f'days since explosion (MJD {args.exp_mjd})')
    print(f'Loaded {len(models)} grid models ({ref_desc}):\n')

    print(f"{'Model':<14} {'E':>5} {'Ni56':>6} {'src':>6} "
          f"{'i_peak':>7} {'pk_day':>6} {'rms_g':>6} {'rms_r':>6} "
          f"{'rms_i':>6} {'rise':>6} {'rms':>6}")
    print('-' * 88)
    for m in models:
        df = m['df']
        t_ipeak = df.loc[df[PEAK_BAND].idxmin(), 'time']
        # Peak alignment: model x = time - t_ipeak (i peak at 0, same as obs)
        shift = -t_ipeak if args.align_peak else 0.0
        m['shift'] = shift
        sc = score_model(df, obs, shift)
        m['score'] = sc
        # Report the model's i-peak time since its own explosion
        t_pk = t_ipeak
        print(f"{m['name']:<14} {m['E']:5.1f} {m['Ni']:6.3f} {m['src']:>6} "
              f"{df[PEAK_BAND].min():7.2f} {t_pk:6.1f} "
              f"{sc.get('g', float('nan')):6.2f} "
              f"{sc.get('r', float('nan')):6.2f} "
              f"{sc.get('i', float('nan')):6.2f} "
              f"{sc.get('rise_rms', float('nan')):6.2f} "
              f"{sc['combined']:6.2f}")

    best = min(models, key=lambda m: m['score']['combined']
               if not np.isnan(m['score']['combined']) else 9e9)
    print(f"\nBest match: {best['name']} "
          f"(combined RMS = {best['score']['combined']:.2f} mag)")
    for band in ('g', 'r', 'i'):
        # Restrict peak search to the first-peak window; later points can
        # include subtraction artifacts (e.g. i=17.5 at MJD 60698).
        det = obs[band]['det']
        near = det[det['phase'] <= 60]
        ref = near if len(near) else det
        print(f"  {band} peak: M_{band} ~ {ref['Mabs'].min():.2f} at phase "
              f"{ref.loc[ref['Mabs'].idxmin(), 'phase']:.1f} d "
              f"({len(det)} det)")

    # ---- Figure: g, r, i — absolute mag vs days since r peak ----
    fig, axes = plt.subplots(3, 1, figsize=(11, 12), sharex=True)
    e_vals = sorted({m['E'] for m in models})
    e_palette = plt.get_cmap('viridis')(
        np.linspace(0.05, 0.95, max(len(e_vals), 2)))
    e_colors = dict(zip(e_vals, e_palette))
    ni_ls = {0.02: '-.', 0.04: ':', 0.075: '-', 0.15: '--'}
    # Envelope-variant suffixes get point markers to distinguish them from
    # the plain E/Ni members at the same (E, Ni) corner.
    suffix_marker = {'_M': 's', '_R': '^'}

    for ax, band, ylabel in ((axes[0], 'g', 'M_g (AB)'),
                             (axes[1], 'r', 'M_r (AB)'),
                             (axes[2], 'i', 'M_i (AB)')):
        det = obs[band]['det']
        nd = obs[band]['nd']
        bc = BAND_COLORS[band]
        ax.errorbar(det['phase'], det['Mabs'], yerr=det['Mabs_err'],
                    fmt='o', color=bc, mfc=bc, mec='k', ms=4, capsize=2,
                    zorder=10, label=f'SN 2024pba {band} det')
        if len(nd):
            ax.scatter(nd['phase'], nd['Mabs_lim'], marker='v',
                       facecolors=bc, edgecolors='k', s=25,
                       zorder=9, label=f'{band} 5-sig lim')
        for m in models:
            df = m['df']
            mark = next((v for k, v in suffix_marker.items()
                         if m['suffix'].startswith(k)), None)
            mevery = max(1, len(df) // 14) if mark else None
            ax.plot(df['time'] + m['shift'], df[band],
                    color=e_colors.get(m['E'], '0.5'),
                    ls=ni_ls.get(m['Ni'], '-'), lw=1.4, alpha=0.9,
                    marker=mark, ms=3.5, markevery=mevery,
                    label=m['name'])
        ax.set_ylabel(ylabel)
        ax.invert_yaxis()
        ax.grid(alpha=0.3)
        ax.legend(fontsize=6.5, ncol=2, loc='lower left')

    if args.align_peak:
        axes[2].set_xlabel(f'Days since {PEAK_BAND}-band peak')
        axes[2].set_xlim(-30, 60)
    else:
        axes[2].set_xlabel(f'Days since explosion (MJD {args.exp_mjd:.1f})')
        axes[2].set_xlim(-5, 90)
    note = ('colour=E_burst; style=Ni56 (dash-dot 0.02, dotted 0.04, '
            'solid 0.075, dashed 0.15); markers = envelope variants '
            '(square=_M, triangle=_R)')
    axes[0].set_title('SN_A_pba grid vs SN 2024pba first peak\n' + note,
                      fontsize=10)
    plt.tight_layout()
    plt.savefig(args.out, dpi=150)
    print(f'\nPlot saved: {args.out}')


if __name__ == '__main__':
    main()
