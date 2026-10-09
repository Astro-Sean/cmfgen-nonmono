#!/usr/bin/env python3
"""Plot SN_B_v2 T6k grid STELLA lightcurves vs SN 2024pba SN_B observations.

The STELLA model lightcurves (the .tt files - the only source used;
the pysynphot *_lightcurve*.csv products are never read) give absolute
Johnson-Cousins UBVRI magnitudes vs time (days since shell launch, t=0).

The observations (observational_data/SN_B_lightcurve.dat) are PSF
photometry giving apparent SDSS g/r magnitudes vs MJD, with
non-detections as upper limits.

Conventions:
  - Redshift z = 0.0186 (TNS UCSC classification report, 2025TNSCR2900)
  - Distance modulus DM = 34.52  (~80 Mpc, H0=70, Omega_m=0.3, Omega_L=0.7)
  - Galactic reddening E(B-V) = 0.15  (IRSA Schlafly & Finkbeiner 2011,
    at RA=25.3077, Dec=49.0635)
  - Observed peak date = MJD 60850.0  (user-specified; confirmed by TNS
    astronote 2025-224: second peak ~2025-06-24)

Both model and observed lightcurves are plotted as "days since peak"
(x=0 = peak), so the peaks align.  Each model is shifted by its own
R-band (or V-band) peak time.

Caveat: the model bands are Johnson-Cousins UBVRI while the observations
are SDSS g and r.  Johnson R is compared to SDSS r and Johnson V to SDSS g;
these are close but not identical passbands, so expect ~0.1-0.3 mag
systematic offsets (especially V vs g).
"""
import argparse
import os
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
GRID_DIR = '/ptmp/sbrennan/SN2024pba/SN_B_v2_T6k_grid'
OBS_FILE = '/ptmp/sbrennan/SN2024pba/observational_data/SN_B_lightcurve.dat'
ATLAS_FILE = '/ptmp/sbrennan/SN2024pba/observational_data/SN2024pba_atclean_cleaned.csv'
SN09IP_FILE = '/ptmp/sbrennan/SN2024pba/observational_data/SN2009ip.csv'
SN10JL_FILE = '/ptmp/sbrennan/SN2024pba/observational_data/SN2010jl.csv'
SN94W_FILE = '/ptmp/sbrennan/SN2024pba/observational_data/SN1994W.csv'
OUT_FILE = os.path.join(GRID_DIR, 'sn_b_t6k_grid_vs_obs.png')

Z = 0.0186                    # TNS redshift (UCSC report)
# Luminosity distance for z=0.0186 (H0=70, Omega_m=0.3, Omega_L=0.7)
# d_L ≈ cz/H0 * (1 + z/2) = 0.0186*299792/70 * 1.0093 ≈ 80.4 Mpc
DL_MPC = 80.4
DM = 5.0 * np.log10(DL_MPC * 1e6 / 10.0)   # ≈ 34.52
EBV = 0.15                    # Galactic E(B-V) from IRSA S&F2011
PEAK_MJD = 60850.0
SNB_EXP_MJD = 60805.0          # adopted SN_B explosion epoch
SNB_EXP_PHASE = SNB_EXP_MJD - PEAK_MJD   # -45 d            # observed peak date (user-specified)
# SN 2009ip comparison lightcurve (already absolute mags in the CSV).
# Phase zero = R-band peak of Event B (2012b), MJD 56207.4.
SN09IP_PEAK_MJD = 56207.4
# SN 2010jl (Ofek+ 2014, J/ApJ/781/42): Time col is days since the
# paper's reference epoch; phase-aligned to its own R peak instead.
# DM = 33.4 (D ~ 48 Mpc), Galactic A_R ~ 0.07.
SN10JL_DM = 33.4
SN10JL_AR = 0.07
# SN 1994W (OSC/astrocatalogs, Sollerman+98 photometry): V band.
# NGC 4041, D ~ 25 Mpc -> DM = 32.0, Galactic A_V ~ 0.05.
SN94W_DM = 32.0
SN94W_AV = 0.05

# SDSS-band extinction coefficients from Schlafly & Finkbeiner (2011), R_V=3.1:
#   R_g = 3.303, R_r = 2.285
A_G = 3.303 * EBV            # ~0.495
A_R = 2.285 * EBV            # ~0.343
# ATLAS o band (lambda_eff ~ 663 nm): A_o/A_V ~ 0.78 (Tonry et al. 2018)
A_O = 2.40 * EBV             # ~0.36


def parse_grid_name(name):
    """Parse a grid model name -> (M Msun, V km/s, W label or '').

    Accepts all conventions:
      bare: M3.0_V10000[_NW]          (SN_B_grid shell mass/velocity)
      old:  T6k_M2.0_V4000[_W3]       (W in units of 1e-3 Msun/yr)
      new:  X_M7p0_V1e4[_W3e-2]       (literal values, p=., e-notation)
    The W label is normalised to units of 1e-3 Msun/yr ('W3', 'W30').
    """
    M = V = None
    W = ''
    for p in name.split('_'):
        if p.startswith('M') and len(p) > 1:
            try:
                M = float(p[1:].replace('p', '.'))
            except ValueError:
                continue
        elif p.startswith('V') and len(p) > 1:
            try:
                V = float(p[1:].replace('p', '.'))
            except ValueError:
                continue
        elif p.startswith('W') and len(p) > 1:
            wtok = p[1:]
            if 'e' in wtok.lower() or 'p' in wtok:
                mdot = float(wtok.replace('p', '.'))
            else:
                mdot = float(wtok) * 1e-3
            W = f'W{mdot * 1e3:g}'
    if M is None or V is None:
        raise ValueError(f'unrecognised grid model name: {name}')
    return M, V, W


def find_tt(mdir):
    """Locate the model's STELLA .tt lightcurve file.

    Search order: canonical stella/ output, the output/ copy, the
    in-progress strad/ run dir, then any stella/<name>*.tt variant.
    backup_* dirs are excluded - stale data.
    """
    name = os.path.basename(mdir)
    candidates = [
        os.path.join(mdir, 'stella', f'{name}.tt'),
        os.path.join(mdir, 'output', f'{name}.tt'),
        os.path.join(mdir, 'stella', 'strad', f'{name}.tt'),
    ]
    candidates += sorted(
        glob.glob(os.path.join(mdir, 'stella', f'{name}*.tt')))
    for path in candidates:
        if os.path.exists(path) and 'backup' not in path:
            return path
    return None


def load_model_tt(tt_path):
    """Load lightcurve from a STELLA .tt file (cols: time Tbb rbb Teff
    Rlast_sc Rtau23 Mbol MU MB MV MI MR Mbolavg gdepos)."""
    rows = []
    with open(tt_path) as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith('time') or s.startswith('XNianm'):
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
    df = pd.DataFrame({
        'time': arr[:, 0], 'U': arr[:, 7], 'B': arr[:, 8], 'V': arr[:, 9],
        'I': arr[:, 10], 'R': arr[:, 11], 'bol_mag': arr[:, 6],
    })
    return df[df['time'] >= 0].copy()


def discover_models(grid_dir):
    """Find one lightcurve per grid model dir.

    Only models with a completed STELLA run are used, i.e. those with
    a .tt file (UBVRI) or the csv lightcurve extracted from it:
      1. stella/<name>.tt               (UBVRI)
      2. output/<name>.tt               (copy of the same)
    Partial .res lightcurves from in-progress runs are excluded.

    Returns list of (name, M, V, W, df, source_tag).
    Diagnostic dirs (*_DIAG, *_TRI) and names that do not parse as
    M<mass>_V<velocity> variants are excluded.
    """
    models = []
    dirs = sorted(d for d in glob.glob(os.path.join(grid_dir, '*'))
                  if os.path.isdir(d))
    for mdir in dirs:
        name = os.path.basename(mdir)
        if not os.path.isdir(mdir):
            continue
        if name.endswith(('_DIAG', '_TRI')):
            continue
        try:
            M, V, W = parse_grid_name(name)
        except (ValueError, IndexError):
            continue

        # Always use the STELLA .tt lightcurve - the *_lightcurve*.csv
        # products (pysynphot re-extraction from .flx) are never used.
        df = None
        src = None
        tt_path = find_tt(mdir)
        if tt_path is not None:
            df = load_model_tt(tt_path)
            src = 'tt'
        if df is not None and len(df):
            models.append((name, M, V, W, df, src))
    models.sort(key=lambda x: (x[1], x[2], x[3], x[0]))
    return models


def load_observations(path):
    raw = pd.read_csv(path)
    obs = {}
    for band, A in (('g', A_G), ('r', A_R)):
        sub = raw[raw['Filter'] == band].copy()
        # Detections: Mag is numeric
        det = sub[sub['Mag'].apply(lambda x: x != '-' and not pd.isna(x))].copy()
        det['Mag'] = det['Mag'].astype(float)
        det['Error'] = det['Error'].astype(float)
        det['phase'] = det['MJD'] - PEAK_MJD       # days since observed peak
        det['Mabs'] = det['Mag'] - DM - A           # absolute magnitude
        det['Mabs_err'] = det['Error']
        # Non-detections: Mag == '-', use 5-sigma limit
        nd = sub[sub['Mag'].apply(lambda x: x == '-' or pd.isna(x))].copy()
        if len(nd):
            nd['phase'] = nd['MJD'] - PEAK_MJD
            nd['Mabs_lim'] = nd['Limit_5p0S2N'].astype(float) - DM - A
        else:
            nd['Mabs_lim'] = pd.Series([], dtype=float)
        obs[band] = {'det': det, 'nd': nd}
    return obs


def load_atlas(path, band='o', sn_min=3.0, chi2_max=5.0):
    """Load ATLAS forced-photometry (atclean) lightcurve.

    Columns used: MJD, m (mag), uJy, duJy_new (flux err), F (filter),
    chi/N, mag5sig.
    Detections: S/N = uJy/duJy_new >= sn_min, chi/N < chi2_max, and a
    clean chip Mask (0x400800/0x400900).  Other mask flags (e.g.
    0x1400802, 0x400a01) mark subtraction artifacts with bogus fluxes.
    Non-detections: S/N < sn_min, plotted at mag5sig.
    Magnitude errors are taken from dm where present, else 1.086/SN.
    Returns {'det': df, 'nd': df} with phase/Mabs columns, like
    load_observations().
    """
    raw = pd.read_csv(path)
    sub = raw[raw['F'] == band].copy()
    sub['sn'] = sub['uJy'] / sub['duJy_new']
    good = (sub['Mask'].isin(['0x400800', '0x400900'])
            & (sub['chi/N'].isna() | (sub['chi/N'] < chi2_max)))

    det = sub[(sub['sn'] >= sn_min) & good].copy()
    det['mag_err'] = det['dm'].where(det['dm'].notna(), 1.086 / det['sn'])
    det['phase'] = det['MJD'] - PEAK_MJD
    det['Mabs'] = det['m'] - DM - A_O
    det['Mabs_err'] = det['mag_err']

    # Rebin detections into nightly bins (inverse-variance weighted)
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

    nd = sub[(sub['sn'] < sn_min) & good].copy()
    nd['phase'] = nd['MJD'] - PEAK_MJD
    nd['Mabs_lim'] = nd['mag5sig'] - DM - A_O

    return {'det': det, 'nd': nd}


def main():
    ap = argparse.ArgumentParser(
        description='Plot SN_B grid STELLA lightcurves vs SN 2024pba '
                    'SN_B observations.')
    ap.add_argument('--grid-dir', default=GRID_DIR,
                    help=f'Grid directory to scan for models '
                         f'(default: {GRID_DIR})')
    ap.add_argument('--out', default=None,
                    help='Output figure path (default: <grid_dir>/'
                         'sn_b_grid_vs_obs.png)')
    args = ap.parse_args()
    grid_dir = args.grid_dir
    out_file = args.out or os.path.join(grid_dir, 'sn_b_grid_vs_obs.png')

    # Discover all model lightcurves in the grid
    models = discover_models(grid_dir)
    print(f"Loaded {len(models)} grid models:")
    for name, M, V, W, df, src in models:
        has_r = df['R'].notna().any()
        r_info = (f"R peak M={df['R'].min():.2f} at "
                  f"t={df.loc[df['R'].idxmin(), 'time']:.1f}d" if has_r
                  else "no R band")
        tpeak_v = df.loc[df['V'].idxmin(), 'time']
        print(f"  {name}: M={M} Msun, V={V} km/s"
              f"{', ' + W if W else ''} [{src}], "
              f"{r_info}, "
              f"V peak M={df['V'].min():.2f} at t={tpeak_v:.1f}d")

    obs = load_observations(OBS_FILE)
    obs['o'] = load_atlas(ATLAS_FILE, band='o')

    # SN 2009ip comparison R-band lightcurve (already absolute mags);
    # phase relative to the 2012b (Event B) R-band peak.  The file
    # interleaves photometry blocks, so sort and nightly-median bin.
    ip09 = pd.read_csv(SN09IP_FILE)
    ip09 = ip09[ip09['R'].notna()].copy()
    ip09['phase'] = ip09['mjd'] - SN09IP_PEAK_MJD
    rows = []
    for night, g in ip09.groupby(np.floor(ip09['phase'])):
        e = g['R_err'].dropna()
        err = np.sqrt(1.0 / (1.0 / e**2).sum()) if len(e) else 0.0
        rows.append({'phase': night, 'R': g['R'].median(),
                     'R_err': err})
    ip09 = pd.DataFrame(rows)

    # SN 2010jl R-band (Ofek+ 2014), apparent -> absolute, phased to its
    # own R peak; nightly-median bin like 09ip.
    jl10 = pd.read_csv(SN10JL_FILE)
    jl10 = jl10[jl10['R'].notna()].copy()
    jl10['phase'] = jl10['mjd'] - jl10.loc[jl10['R'].idxmin(), 'mjd']
    jl10['Mabs'] = jl10['R'] - SN10JL_DM - SN10JL_AR
    rows = []
    for night, g in jl10.groupby(np.floor(jl10['phase'])):
        e = g['R_err'].dropna()
        err = np.sqrt(1.0 / (1.0 / e**2).sum()) if len(e) else 0.0
        rows.append({'phase': night, 'Mabs': g['Mabs'].median(),
                     'err': err})
    jl10 = pd.DataFrame(rows)

    # SN 1994W V-band (OSC/Sollerman+98), absolute mags, phased to its
    # own V peak; nightly-median bin.
    w94 = pd.read_csv(SN94W_FILE)
    w94 = w94[w94['V'].notna()].copy()
    w94['phase'] = w94['mjd'] - w94.loc[w94['V'].idxmin(), 'mjd']
    w94['Mabs'] = w94['V'] - SN94W_DM - SN94W_AV
    w94 = (w94.assign(night=np.floor(w94['phase']))
              .groupby('night')['Mabs'].median().reset_index()
              .rename(columns={'night': 'phase'}))
    print(f"\nObservations: "
          f"{len(obs['r']['det'])} r-band det, "
          f"{len(obs['g']['det'])} g-band det, "
          f"{len(obs['o']['det'])} ATLAS o-band det")
    r_det = obs['r']['det']
    g_det = obs['g']['det']
    o_det = obs['o']['det']
    print(f"  r peak: Mabs={r_det['Mabs'].min():.2f} at phase="
          f"{r_det.loc[r_det['Mabs'].idxmin(), 'phase']:.1f}d")
    print(f"  g peak: Mabs={g_det['Mabs'].min():.2f} at phase="
          f"{g_det.loc[g_det['Mabs'].idxmin(), 'phase']:.1f}d")
    print(f"  o peak: Mabs={o_det['Mabs'].min():.2f} at phase="
          f"{o_det.loc[o_det['Mabs'].idxmin(), 'phase']:.1f}d")

    # Tag for non-M/V/W suffixes (e.g. _we0p15 geometry variants)
    def suffix_of(name):
        toks = name.split('_')[1:]
        extra = [t for t in toks if not t.startswith(('M', 'V', 'W'))]
        return '_'.join(extra)

    models_ext = [(n, M, V, W, suffix_of(n), df, src)
                  for n, M, V, W, df, src in models]

    def lbl(m):
        n, M, V, W, suf, df, src = m
        return (f'M={M:g}, V={V:g}' + (f', {W}' if W else '')
                + (f', {suf}' if suf else ''))

    # ---- Panel slices: vary one parameter, hold the rest fixed ----
    # 1) fixed (M, W) -> vary V      2) fixed (V, W) -> vary M
    # 3) fixed (M, V) -> vary W/geometry
    def mdot_str(W):
        """'W10' -> 'Ṁ_wind = 1e-2 M_sun/yr'; '' -> base wind Ṁ=1e-3."""
        mdot = float(W[1:]) * 1e-3 if W else 1e-3
        mant, exp = f'{mdot:.0e}'.split('e')
        sci = f'${float(mant):g}\\times10^{{{int(exp)}}}$'
        return (f'$\\dot{{M}}_{{\\rm wind}}$ = {sci} '
                f'M$_\\odot$ yr$^{{-1}}$')

    slices = []
    for keyfunc, varied, tit in (
            (lambda m: (m[1], m[3], m[4]), 'V', 'vary V'),
            (lambda m: (m[2], m[3], m[4]), 'M', 'vary M'),
            (lambda m: (m[1], m[2]), 'W', 'vary wind')):
        groups = {}
        for m in models_ext:
            groups.setdefault(keyfunc(m), []).append(m)
        for k, ms in sorted(groups.items(), key=lambda kv: str(kv[0])):
            if len(ms) < 2:
                continue
            if varied == 'V':
                title = (f'M$_{{\\rm shell}}$ = {k[0]:g} M$_\\odot$, '
                         f'{mdot_str(k[1])}'
                         + (f', {k[2]}' if k[2] else '')
                         + ' — vary v$_{\\rm shell}$')
            elif varied == 'M':
                title = (f'v$_{{\\rm shell}}$ = {k[0]:g} km s$^{{-1}}$, '
                         f'{mdot_str(k[1])}'
                         + (f', {k[2]}' if k[2] else '')
                         + ' — vary M$_{\\rm shell}$')
            else:
                title = (f'M$_{{\\rm shell}}$ = {k[0]:g} M$_\\odot$, '
                         f'v$_{{\\rm shell}}$ = {k[1]:g} km s$^{{-1}}$'
                         ' — vary wind/geometry')
            slices.append((title, varied, ms))

    fig, axes = plt.subplots(max(1, (len(slices) + 1) // 2), 2,
                             figsize=(17, 5.2 * ((len(slices) + 1) // 2)),
                             sharex=True, sharey=True)
    axes = np.atleast_1d(axes).ravel()

    # x-range highlights the SN_B event only
    x_lo, x_hi = -80.0, 100.0

    # Brightest point in the plotted window (obs or model) sets the top
    # of the inverted magnitude axis, with a small pad.
    brightest = []
    for obs_band in ('r', 'g', 'o'):
        det = obs[obs_band]['det']
        in_win = det[(det['phase'] >= x_lo) & (det['phase'] <= x_hi)]
        if len(in_win):
            brightest.append(in_win['Mabs'].min())
    for name, M, V, W, df, src in models:
        for mod_band in ('R', 'V'):
            band = df[mod_band]
            if not band.notna().any():
                continue
            t_peak = df.loc[band.idxmin(), 'time']
            in_win = band[(df['time'] - t_peak >= x_lo)
                          & (df['time'] - t_peak <= x_hi)]
            if in_win.notna().any():
                brightest.append(in_win.min())
    y_top = min(brightest) - 0.3 if brightest else -18.8

    # Observations: all three bands on every panel, coloured markers
    # with black outlines.
    band_style = {'r': dict(color='tab:red', marker='o', label='r (PSF)'),
                  'g': dict(color='tab:green', marker='s', label='g (PSF)'),
                  'o': dict(color='darkorange', marker='^',
                            label='o (ATLAS)')}

    def plot_obs(ax):
        for b, st in band_style.items():
            det = obs[b]['det']
            ax.errorbar(det['phase'], det['Mabs'], yerr=det['Mabs_err'],
                        fmt=st['marker'], color=st['color'], ms=5,
                        mec='k', mew=0.7, capsize=2, elinewidth=0.7,
                        alpha=0.85, zorder=20)

    # Obs band legend: single row across the top of the figure
    from matplotlib.lines import Line2D
    obs_handles = [Line2D([], [], marker=st['marker'], color=st['color'],
                          mec='k', mew=0.7, ms=7, ls='',
                          label=f"{st['label']}")
                   for b, st in band_style.items()]
    obs_handles.append(Line2D([], [], marker='o', color='0.55',
                              mec='k', mew=0.7, ms=7, ls='',
                              label='SN 2009ip R'))
    obs_handles.append(Line2D([], [], color='lightblue', lw=1.4,
                              label='SN 2010jl R'))
    obs_handles.append(Line2D([], [], color='plum', lw=1.4,
                              label='SN 1994W V'))
    obs_handles.append(Line2D([], [], color='k', ls=':', lw=1.2,
                              label='SN_B explosion (est.)'))
    fig.legend(handles=obs_handles, ncol=7, loc='upper center',
               frameon=False, fontsize=10, bbox_to_anchor=(0.5, 1.0))

    # Model curves are coloured by band (model R red, model V green,
    # matching the obs markers); linestyle distinguishes models within
    # a panel.
    model_band_color = {'R': 'tab:red', 'V': 'tab:green'}
    model_ls = ['-', '--', '-.', ':', (0, (5, 1)), (0, (1, 1)),
                (0, (3, 1, 1, 1)), (0, (3, 2, 1, 2, 1, 2))]

    for ax, (title, varied, ms_) in zip(axes, slices):
        plot_obs(ax)
        ax.errorbar(ip09['phase'], ip09['R'], yerr=ip09['R_err'],
                    fmt='o', color='0.55', ms=4, capsize=2,
                    elinewidth=0.6, alpha=0.6, zorder=0)
        jl10w = jl10[(jl10['phase'] >= x_lo) & (jl10['phase'] <= x_hi)]
        ax.plot(jl10w['phase'], jl10w['Mabs'], color='lightblue',
                lw=2.0, alpha=0.9, zorder=0)
        w94w = w94[(w94['phase'] >= x_lo) & (w94['phase'] <= x_hi)]
        ax.plot(w94w['phase'], w94w['Mabs'], color='plum',
                lw=2.0, alpha=0.9, zorder=0)
        ax.axvline(SNB_EXP_PHASE, color='k', ls=':', lw=1.2,
                   alpha=0.6, zorder=1)
        for i, m in enumerate(ms_):
            n, M, V, W, suf, df, src = m
            lab = ({'V': f'v$_{{\\rm shell}}$ = {V:g} km s$^{{-1}}$',
                    'M': f'M$_{{\\rm shell}}$ = {M:g} M$_\\odot$',
                    'W': mdot_str(W) + (f' {suf}' if suf else '')}[varied])
            ls = model_ls[i % len(model_ls)]
            for mod_band in ('R', 'V'):
                band = df[mod_band]
                if not band.notna().any():
                    continue
                t_peak = df.loc[band.idxmin(), 'time']
                ax.plot(df['time'] - t_peak, band,
                        color=model_band_color[mod_band],
                        linestyle=ls, linewidth=1.8, alpha=0.9,
                        label=lab if mod_band == 'R' else None)
        ax.set_title(title, fontsize=11)
        ax.set_xlabel('Phase (days since model peak)')
        ax.grid(True, alpha=0.3)
        ax.set_xlim(left=x_lo, right=x_hi)
        ax.set_ylim(bottom=-14.2, top=y_top)
        leg1 = ax.legend(fontsize=7, ncol=1, frameon=False, loc='best')
        ax.add_artist(leg1)
        name_handles = [
            Line2D([], [], color='0.3', ls=model_ls[i % len(model_ls)],
                   lw=1.5, label=m[0])
            for i, m in enumerate(ms_)]
        ax.legend(handles=name_handles, fontsize=5.5, ncol=1,
                  frameon=False, loc='lower left')
    for ax in axes[len(slices):]:
        ax.set_visible(False)
    for i in range(0, len(slices), 2):
        axes[i].set_ylabel('Absolute magnitude')

    # Filter-mismatch caveat
    fig.text(0.5, 0.015,
             'Caveat: model curves are Johnson-Cousins R (red) and V '
             '(green); linestyle distinguishes models within a panel.  '
             'Observations are SDSS r,g and ATLAS o.  Expect ~0.1-0.3 mag '
             'systematic offsets (V vs g and R vs o are the largest).  '
             'Only models with completed STELLA runs (.tt photometry) '
             'are shown.',
             ha='center', fontsize=9, style='italic', color='0.3')

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(out_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\nSaved: {out_file}")


if __name__ == '__main__':
    main()
