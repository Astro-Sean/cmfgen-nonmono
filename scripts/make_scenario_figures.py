#!/usr/bin/env python3
"""Generate all figures for latex/scenario.tex (SN 2024pba two-event scenario).

Reads:
  - observational_data/SN2024pba_atclean_cleaned.csv   (ATLAS o/c photometry)
  - observational_data/SN_A_lightcurve.dat             (PSF g/r around first peak)
  - observational_data/SN_B_lightcurve.dat             (PSF g/r around second peak)
  - observational_data/SN2024pba_*.txt                 (DBSP + ALFOSC spectra)
  - SN_A/SN_A_g_T20k_Ni06 output lightcurve / input profile
  - SN_B_v2_T6k_grid model lightcurves, input profiles and STELLA .swd snapshots

Writes PDF figures into latex/figures/.

Distance / extinction adopted (scenario.md):
  D_L = 80.4 Mpc -> DM = 34.53,  z = 0.0186,  E(B-V) = 0.15
  A_g = 3.303 E, A_r = 2.285 E, A_o ~ 2.40 E, A_c ~ 3.3 E (approx g)
"""
import os
import re
import sys
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, Circle, FancyArrowPatch

ROOT = '/viper/ptmp1/sbrennan/SN2024pba'
OBS = os.path.join(ROOT, 'observational_data')
GRID = os.path.join(ROOT, 'SN_B_v2_T6k_grid')
SN_A_DIR = os.path.join(ROOT, 'SN_A', 'SN_A_g_T20k_Ni06')
FIGDIR = os.path.join(ROOT, 'latex', 'figures')
os.makedirs(FIGDIR, exist_ok=True)

sys.path.insert(0, os.path.join(ROOT, 'scripts'))
from convert_stella_swd_to_cmfgen import read_swd, normalize_swd_composition  # noqa: E402

Z = 0.0186
DL_MPC = 80.4
DM = 5.0 * np.log10(DL_MPC * 1e6 / 10.0)        # 34.53
EBV = 0.15
A_G, A_R = 3.303 * EBV, 2.285 * EBV
A_O, A_C = 2.40 * EBV, 3.30 * EBV
PEAK_MJD = 60850.0
SNA_PEAK_MJD = 60513.0
SNA_EXP_MJD = 60490.0
BUMP_LO, BUMP_HI = 60810.0, 60830.0

plt.rcParams.update({
    'font.size': 9, 'axes.labelsize': 9, 'axes.titlesize': 9,
    'legend.fontsize': 7.5, 'xtick.labelsize': 8, 'ytick.labelsize': 8,
    'axes.linewidth': 0.8, 'lines.linewidth': 1.1,
    'figure.dpi': 150, 'savefig.bbox': 'tight',
})

C_O = '#d95f02'          # ATLAS o (orange)
C_C = '#1b9e77'          # ATLAS c (cyan)
C_R = '#7570b3'          # PSF r
C_G = '#66a61e'          # PSF g


def save(fig, name):
    out = os.path.join(FIGDIR, name)
    fig.savefig(out)
    plt.close(fig)
    print('wrote', out)


def load_atlas():
    df = pd.read_csv(os.path.join(OBS, 'SN2024pba_atclean_cleaned.csv'))
    # drop rows with NaN dm or extreme residual flux (obvious outliers)
    df = df[np.isfinite(df['dm'])].copy()
    df = df[df['duJy_new'] < 60].copy()
    out = {}
    for f, A in (('o', A_O), ('c', A_C)):
        s = df[df['F'] == f].copy()
        s['Mabs'] = s['m'] - DM - A
        out[f] = s
    return out


def load_psf(path):
    raw = pd.read_csv(path)
    obs = {}
    for band, A in (('g', A_G), ('r', A_R)):
        sub = raw[raw['Filter'] == band].copy()
        det = sub[sub['Mag'].apply(lambda x: x != '-' and not pd.isna(x))].copy()
        det['Mag'] = det['Mag'].astype(float)
        det['Error'] = det['Error'].astype(float)
        det['Mabs'] = det['Mag'] - DM - A
        obs[band] = det
    return obs


# ---------------------------------------------------------------------------
# Fig. 1: observed ATLAS lightcurve with event anchors
# ---------------------------------------------------------------------------
def fig_observed_lc():
    at = load_atlas()
    fig, ax = plt.subplots(figsize=(6.6, 3.4))
    for f, c, lab in (('c', C_C, 'ATLAS c'), ('o', C_O, 'ATLAS o')):
        s = at[f]
        ax.errorbar(s['MJD'], s['Mabs'], yerr=s['dm'], fmt='o', ms=1.6,
                    color=c, alpha=0.55, elinewidth=0.4, capsize=0,
                    label=lab, zorder=2)
    anchors = [
        (SNA_EXP_MJD, 'SN$_{\\rm A}$ explosion'),
        (SNA_PEAK_MJD, 'SN$_{\\rm A}$ peak'),
        (60850.0, 'SN$_{\\rm B}$ main peak'),
    ]
    for mjd, lab in anchors:
        ax.axvline(mjd, color='0.4', lw=0.7, ls='--', zorder=1)
        ax.text(mjd - 3, -20.3, lab, rotation=90, fontsize=7, va='top',
                ha='right', color='0.25')
    ax.axvspan(BUMP_LO, BUMP_HI, color='0.85', zorder=0)
    ax.text(BUMP_LO - 4, -14.2, 'precursor\nbump', fontsize=7,
            ha='right', va='top', color='0.25')
    ax.annotate('$o\\approx-18.2$', xy=(60850, -18.2),
                xytext=(60880, -18.8), fontsize=7,
                arrowprops=dict(arrowstyle='->', lw=0.6, color='0.3'),
                color='0.2')
    ax.annotate('$M\\approx-17.2$', xy=(60515, -17.2),
                xytext=(60530, -16.9), fontsize=7,
                arrowprops=dict(arrowstyle='->', lw=0.6, color='0.3'),
                color='0.2')
    ax.set_xlim(60450, 61160)
    ax.set_ylim(-12.5, -20.5)
    ax.set_xlabel('MJD')
    ax.set_ylabel('absolute magnitude')
    ax.legend(loc='lower left', frameon=False, markerscale=2)
    save(fig, 'fig_obs_lc.pdf')


# ---------------------------------------------------------------------------
# Fig. 2: schematic geometry
# ---------------------------------------------------------------------------
def fig_schematic():
    fig, ax = plt.subplots(figsize=(6.6, 6.4))
    ax.set_aspect('equal')
    ax.axis('off')

    # concentric regions (radii chosen for legibility, not to scale):
    # shell | wind gap | SN_A ejecta   (buffer omitted for clarity)
    r_shell, r_gap, r_ej = 0.60, 1.75, 2.30
    ax.add_patch(Wedge((0, 0), r_ej, 0, 360, width=r_ej - r_gap,
                       fc='#c2a5cf', ec='none', alpha=0.7))
    ax.add_patch(Wedge((0, 0), r_gap, 0, 360, width=r_gap - r_shell,
                       fc='#e0f3f8', ec='none', alpha=0.95))
    ax.add_patch(Wedge((0, 0), r_shell, 0, 360, fc='#80cdc1', ec='none'))
    ax.add_patch(Wedge((0, 0), r_ej, 0, 360, fill=False, ec='0.3', lw=0.8))
    for rr in (r_shell, r_gap):
        ax.add_patch(Circle((0, 0), rr, fill=False, ec='0.3', lw=0.6))

    # forward shock at the shell/wind interface (outer edge of the shell)
    r_sh = r_shell
    ax.add_patch(Circle((0, 0), r_sh, fill=False, ec='#b2182b', lw=2.0))

    # radiative heating ahead of the shock: hatching that gets sparser
    # with distance from the shell/wind interface
    heat_bands = [
        (r_shell, 1.00, '//////'),
        (1.00, 1.40, '////'),
        (1.40, r_gap, '///'),
        (r_gap, 2.05, '//'),
        (2.05, r_ej, '/'),
    ]
    for r0, r1, h in heat_bands:
        ax.add_patch(Wedge((0, 0), r1, 0, 360, width=r1 - r0,
                           fc='none', ec='#d05a1e', hatch=h, lw=0,
                           alpha=0.55))
    # shell expansion arrows
    for ang in (210, 250, 290, 330):
        a = np.deg2rad(ang)
        ax.add_patch(FancyArrowPatch(
            (0.30 * np.cos(a), 0.30 * np.sin(a)),
            (0.52 * np.cos(a), 0.52 * np.sin(a)),
            arrowstyle='-|>', mutation_scale=7, color='0.25', lw=0.8))

    def leader(xy, xytext, text, ha='left'):
        ax.annotate(text, xy=xy, xytext=xytext, fontsize=8.5, ha=ha,
                    va='center', color='0.1',
                    arrowprops=dict(arrowstyle='-', lw=0.5, color='0.4'))

    a45 = np.deg2rad(45)
    leader((0.3 * np.cos(np.deg2rad(80)), 0.3 * np.sin(np.deg2rad(80))),
           (0.0, 3.15),
           'SN$_{\\rm B}$ shell: $M_{\\rm sh}\\approx5$-$10$ M$_\\odot$,\n'
           '$v_{\\rm sh}\\approx10^{4}$ km s$^{-1}$, solar composition',
           ha='center')
    leader((1.2 * np.cos(np.deg2rad(-25)), 1.2 * np.sin(np.deg2rad(-25))),
           (2.5, 0.35),
           'wind gap: $\\dot{M}\\approx10^{-3}$-$3\\times10^{-2}$\n'
           'M$_\\odot$ yr$^{-1}$, $v_{\\rm w}=100$ km s$^{-1}$,\n'
           '$\\rho\\propto r^{-2}$', ha='left')
    leader((2.02 * np.cos(np.deg2rad(-135)), 2.02 * np.sin(np.deg2rad(-135))),
           (-0.6, -2.6),
           'SN$_{\\rm A}$ ejecta at 400 d: metal rich\n'
           '(He,C,O,Si,Ca,Fe), homologous,\n'
           '$v\\approx1800$-$10\\,000$ km s$^{-1}$', ha='right')
    leader((r_sh * np.cos(-a45), r_sh * np.sin(-a45)), (2.5, -1.0),
           'forward shock at the\nshell/wind interface', ha='left')
    leader((1.55 * np.cos(np.deg2rad(30)), 1.55 * np.sin(np.deg2rad(30))),
           (2.5, 1.85),
           'radiatively heated material\n(heating declines outward)', ha='left')
    ax.annotate('observer', xy=(-2.37, 0.25), xytext=(-2.37, 0.45),
                fontsize=8.5, color='0.2', va='bottom', ha='left')
    ax.add_patch(FancyArrowPatch((-2.37, 0.18), (-2.18, 0.18),
                                 arrowstyle='-|>', mutation_scale=10,
                                 color='0.2', lw=1.1))

    ax.set_xlim(-2.4, 2.4)
    ax.set_ylim(-3.15, 3.3)
    save(fig, 'fig_schematic.pdf')


# ---------------------------------------------------------------------------
# Fig. 3: initial model structure (density, velocity, composition)
# ---------------------------------------------------------------------------
def _load_input_profile(path):
    arr = np.loadtxt(path, comments='#')
    return pd.DataFrame(arr, columns=[
        'r', 'rho', 'p', 'u', 'T', 'fH', 'fHe', 'fO', 'fSi', 'fCa', 'fFe', 'fNi'])


def fig_input_structure():
    prof = _load_input_profile(os.path.join(
        GRID, 'T6k_M10.0_V10000', 'stella', 'input_profile'))
    r, rho, v = prof['r'], prof['rho'], prof['u'] / 1e5
    fH, fHe, fO = prof['fH'], prof['fHe'], prof['fO']
    fSi, fCa, fFe = prof['fSi'], prof['fCa'], prof['fFe']

    # region boundaries: shell inner edge r_min; shell->wind where drho slope
    # breaks to r^-2; buffer->ejecta where metal fractions rise (fO > 0.02)
    r_min = r.min()
    i_ej = np.where(fO.values > 0.02)[0]
    r_ej_in = r[i_ej].min() if len(i_ej) else np.nan
    # wind->buffer: strong density upturn before ejecta; take radius of max
    # second-derivative / simplest: where rho rises above smooth r^-2 near edge
    i_w = (r < r_ej_in)
    rho_w = rho[i_w]
    r_w = r[i_w]
    # shell outer edge: max of rho plateau inside wind (highest-density zone)
    i_sh = np.argmax(rho.values * (r.values < r_w.values[len(r_w) // 3]
                                   if len(r_w) else 1))
    r_shell_out = r.iloc[i_sh]
    # buffer outer edge: where H mass fraction drops below 0.5 outside wind
    i_buf = np.where((r > r_w.max()) & (fH > 0.5))[0]
    r_buf_out = r[i_buf].max() if len(i_buf) else r_ej_in

    regions = [
        (r_min, r_shell_out, '#80cdc1', 'SN$_{\\rm B}$ shell'),
        (r_shell_out, r_w.max(), '#e0f3f8', 'wind'),
        (r_w.max(), r_ej_in, '#fee0b6', 'buffer'),
        (r_ej_in, r.max(), '#c2a5cf', 'SN$_{\\rm A}$ ejecta'),
    ]

    fig, axes = plt.subplots(3, 1, figsize=(6.4, 6.4), sharex=True)
    for ax in axes:
        for r0, r1, c, _ in regions:
            ax.axvspan(r0, r1, color=c, alpha=0.45, lw=0)
    axes[0].loglog(r, rho, 'k-', lw=1.0)
    axes[0].set_ylabel(r'$\rho$ [g cm$^{-3}$]')
    axes[1].semilogx(r, v, 'k-', lw=1.0)
    axes[1].set_ylabel(r'$v$ [km s$^{-1}$]')
    for f, lab, c in ((fH, 'H', '#1b9e77'), (fHe, 'He', '#7570b3'),
                      (fO, 'O', '#d95f02'), (fSi, 'Si', '#e7298a'),
                      (fCa, 'Ca', '#66a61e'), (fFe, 'Fe+Ni', '#a6761d')):
        fsum = f + (prof['fNi'] if lab == 'Fe+Ni' else 0)
        axes[2].semilogx(r, np.maximum(fsum, 1e-6), lw=1.0, label=lab, color=c)
    axes[2].set_ylabel('mass fraction')
    axes[2].set_xlabel('radius [cm]')
    axes[2].set_ylim(8e-4, 1.6)
    axes[2].set_yscale('log')
    axes[2].legend(loc='center left', frameon=False, ncol=6, columnspacing=0.9,
                   handlelength=1.2)
    ylab = [3e-9, 3e-8, 3e-9, 3e-8]
    for (r0, r1, c, lab), yy in zip(regions, ylab):
        axes[0].text(np.sqrt(r0 * r1), yy, lab, rotation=90, fontsize=6.5,
                     ha='center', va='bottom', color='0.2')
    axes[0].set_ylim(5e-18, 8e-8)
    fig.align_ylabels(axes)
    save(fig, 'fig_input_structure.pdf')


# ---------------------------------------------------------------------------
# Fig. 4: SN_A lightcurve vs observed first peak
# ---------------------------------------------------------------------------
def fig_sn_a_lc():
    lc = pd.read_csv(os.path.join(
        SN_A_DIR, 'output', 'SN_A_g_T20k_Ni06_lightcurve.csv'))
    lc = lc.drop_duplicates(subset='time')
    lc = lc[np.isfinite(lc['R']) & (lc['R'] < 0)]  # drop glitch rows
    # model time is days since explosion; align explosion (t=0) to the
    # adopted SN_A explosion epoch MJD 60490 (= -23 d relative to SN_A peak)
    ph = lc['time'] - (SNA_PEAK_MJD - SNA_EXP_MJD)

    obs = load_psf(os.path.join(OBS, 'SN_A_lightcurve.dat'))
    at = load_atlas()

    fig, ax = plt.subplots(figsize=(6.4, 3.2))
    ax.plot(ph, lc['R'], color='k', lw=1.2, label='SN\\_A model $R$')
    ax.plot(ph, lc['bol_mag'], color='0.55', lw=1.0, ls='--',
            label='SN\\_A model bol.')
    for band, c, lab in (('r', C_R, 'obs $r$'), ('g', C_G, 'obs $g$')):
        d = obs[band]
        ax.errorbar(d['MJD'] - SNA_PEAK_MJD, d['Mabs'], yerr=d['Error'],
                    fmt='o', ms=2.5, color=c, elinewidth=0.5, capsize=0,
                    label=lab, zorder=3)
    s = at['o']
    sel = (s['MJD'] > 60495) & (s['MJD'] < 60560)
    ax.errorbar(s.loc[sel, 'MJD'] - SNA_PEAK_MJD, s.loc[sel, 'Mabs'],
                yerr=s.loc[sel, 'dm'], fmt='s', ms=2.0, color=C_O, alpha=0.7,
                elinewidth=0.5, capsize=0, label='ATLAS $o$')
    ax.axvline(0, color='0.4', lw=0.7, ls='--')
    ax.axvline(SNA_EXP_MJD - SNA_PEAK_MJD, color='0.4', lw=0.7, ls=':')
    ax.text(1.5, -11.4, 'SN$_{\\rm A}$ peak (MJD 60513)', fontsize=7,
            color='0.3')
    ax.text(SNA_EXP_MJD - SNA_PEAK_MJD + 1, -11.9,
            'explosion (MJD $\\sim$60490)', fontsize=7, color='0.3',
            ha='left')
    ax.set_xlabel('days since SN$_{\\rm A}$ peak')
    ax.set_ylabel('absolute magnitude')
    ax.set_xlim(-30, 60)
    ax.set_ylim(-11, -19)
    ax.legend(frameon=False, loc='lower right')
    save(fig, 'fig_sn_a_lc.pdf')


# ---------------------------------------------------------------------------
# Fig. 5: T6k_M10.0_V10000 lightcurve vs observed SN_B lightcurve
# ---------------------------------------------------------------------------
def fig_model_lc():
    name = 'T6k_M10.0_V10000'
    lc = pd.read_csv(os.path.join(
        GRID, name, 'output', f'{name}_lightcurve.csv'))
    lc = lc[lc['time'] >= -5].drop_duplicates(subset='time')
    t_peak = lc.loc[lc['R'].idxmin(), 'time']
    ph = lc['time'] - t_peak

    obs = load_psf(os.path.join(OBS, 'SN_B_lightcurve.dat'))
    at = load_atlas()

    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    ax.plot(ph, lc['R'], color='k', lw=1.3,
            label=f'{name.replace("_", "\\_")} $R$')
    ax.plot(ph, lc['bol_mag'], color='0.55', lw=1.0, ls='--',
            label='model bolometric')
    d = obs['r']
    ax.errorbar(d['MJD'] - PEAK_MJD, d['Mabs'], yerr=d['Error'], fmt='o',
                ms=2.4, color=C_R, elinewidth=0.5, capsize=0,
                label='obs $r$', zorder=3)
    s = at['o']
    sel = (s['MJD'] > 60780) & (s['MJD'] < 61160)
    ax.errorbar(s.loc[sel, 'MJD'] - PEAK_MJD, s.loc[sel, 'Mabs'],
                yerr=s.loc[sel, 'dm'], fmt='s', ms=2.0, color=C_O, alpha=0.7,
                elinewidth=0.5, capsize=0, label='ATLAS $o$', zorder=2)
    ax.axvspan(BUMP_LO - PEAK_MJD, BUMP_HI - PEAK_MJD, color='0.88', zorder=0)
    ax.text(-30, -18.6, 'precursor\nbump', fontsize=7, ha='center',
            va='top', color='0.3')
    ax.annotate('shock-cooling /\nwind diffusion bump', xy=(-38, -16.6),
                xytext=(-82, -16.3), fontsize=7,
                arrowprops=dict(arrowstyle='->', lw=0.6, color='0.3'),
                color='0.25')
    ax.annotate('main peak\n$M\\approx-19.4$', xy=(0, -19.4),
                xytext=(25, -18.3),
                fontsize=7, arrowprops=dict(arrowstyle='->', lw=0.6,
                                          color='0.3'), color='0.25')
    ax.axvline(0, color='0.4', lw=0.7, ls='--')
    ax.set_xlabel('days since SN$_{\\rm B}$ peak (MJD 60850)')
    ax.set_ylabel('absolute magnitude')
    ax.set_xlim(-90, 150)
    ax.set_ylim(-13, -20)
    ax.legend(frameon=False, loc='lower left')
    save(fig, 'fig_model_lc.pdf')


# ---------------------------------------------------------------------------
# Fig. 6: radiative heating of SN_A ejecta (swd snapshots)
# ---------------------------------------------------------------------------
def fig_heating():
    swd = os.path.join(
        GRID, 'T6k_M10.0_V10000', 'stella', 'backup_20260915_114237',
        'T6k_M10.0_V10000.swd')
    epochs = [20, 40, 57, 75]
    cols = plt.cm.viridis(np.linspace(0.05, 0.85, len(epochs)))

    fig, axes = plt.subplots(3, 1, figsize=(6.4, 6.6), sharex=True)
    for targ, c in zip(epochs, cols):
        ep = read_swd(swd, target_days=targ)
        tt = [k for k, v in ep.items() if v is not None][0]
        s = ep[tt]
        r, v, T, Tr, rho = s['r'], s['v'], s['T'], s['T_rad'], s['rho']
        mf = normalize_swd_composition(s['composition'])
        O = mf['O']
        i_ej = np.where(O > 0.05)[0]
        r_ej = r[i_ej].min()
        # shock: largest |dv| inside the ejecta inner edge
        order = np.argsort(r)
        rs, vs = r[order], v[order]
        dv = np.abs(np.diff(vs))
        j = np.argmax(dv[rs[:-1] < r_ej])
        r_sh = rs[j]
        face = i_ej[np.argmin(r[i_ej])]

        axes[0].loglog(r, rho, color=c, lw=0.9)
        axes[1].semilogx(r, v, color=c, lw=0.9)
        axes[2].loglog(r, T, color=c, lw=1.1,
                       label=f'$t={tt:.0f}$ d')
        axes[2].loglog(r, Tr, color=c, lw=0.7, ls=':')
        for ax in axes:
            ax.axvline(r_ej, color=c, lw=0.6, ls='--', alpha=0.8)
            ax.axvline(r_sh, color=c, lw=0.6, ls=':', alpha=0.8)
        axes[2].plot([r_ej], [T[face]], 'o', ms=3, color=c, zorder=5)

    axes[0].set_ylabel(r'$\rho$ [g cm$^{-3}$]')
    axes[1].set_ylabel(r'$v$ [km s$^{-1}$]')
    axes[2].set_ylabel(r'$T$ [K]')
    axes[2].set_xlabel('radius [cm]')
    axes[2].set_xlim(2e14, 5e16)

    # annotate line styles on middle panel
    axes[1].plot([], [], color='0.3', ls='--', lw=0.7,
                 label='SN$_{\\rm A}$ ejecta inner edge')
    axes[1].plot([], [], color='0.3', ls=':', lw=0.7,
                 label='velocity jump (shock)')
    axes[1].legend(frameon=False, loc='lower right', fontsize=7)
    axes[2].plot([], [], 'k-', lw=1.0, label='$T_{\\rm gas}$')
    axes[2].plot([], [], 'k:', lw=0.9, label='$T_{\\rm rad}$')
    axes[2].plot([], [], 'ko', ms=3,
                 label='$T_{\\rm gas}$ at ejecta inner edge')
    axes[2].legend(frameon=False, loc='center left', fontsize=7)
    fig.align_ylabels(axes)
    save(fig, 'fig_heating.pdf')


# ---------------------------------------------------------------------------
# Fig. 7: grid comparison (all X_*/T6k_* lightcurves vs obs)
# ---------------------------------------------------------------------------
_RES_LC_RE = re.compile(
    r'^\s*(-?\d+\.\d+)\s+E\s+'
    r'(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+'
    r'(\d+\.\d+E[+-]\d+)\s+'
    r'(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+'
    r'(-?\d+\.\d+)\s+(-?\d+\.\d+)')


def _load_model_res(res_path):
    rows = []
    with open(res_path) as fh:
        for line in fh:
            m = _RES_LC_RE.match(line.strip())
            if m:
                rows.append([float(m.group(k)) for k in range(1, 13)])
    if not rows:
        return None
    arr = np.array(rows)
    _, idx = np.unique(arr[:, 0], return_index=True)
    arr = arr[np.sort(idx)]
    df = pd.DataFrame({'time': arr[:, 0], 'bol_mag': arr[:, 6],
                       'V': arr[:, 9], 'R': np.nan})
    return df[df['time'] >= -5].copy()


def _load_model_tt(tt_path):
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
    return pd.DataFrame({'time': arr[:, 0], 'bol_mag': arr[:, 6],
                         'V': arr[:, 9], 'R': arr[:, 11]})


def fig_grid():
    models = []
    for mdir in sorted(glob.glob(os.path.join(GRID, 'X_*'))
                       + glob.glob(os.path.join(GRID, 'T6k_*'))):
        name = os.path.basename(mdir)
        if name.endswith(('_DIAG', '_TRI')):
            continue
        csv_p = os.path.join(mdir, 'output', f'{name}_lightcurve.csv')
        tt_p = os.path.join(mdir, 'stella', f'{name}.tt')
        res_p = os.path.join(mdir, 'stella', 'strad', f'{name}.res')
        df = src = None
        if os.path.exists(csv_p):
            df = pd.read_csv(csv_p).drop_duplicates(subset='time')
            df = df[df['time'] >= -5]
            src = 'csv'
        elif os.path.exists(tt_p):
            df, src = _load_model_tt(tt_p), 'tt'
        elif os.path.exists(res_p):
            df, src = _load_model_res(res_p), 'res'
        if df is not None and len(df):
            models.append((name, df, src))

    obs = load_psf(os.path.join(OBS, 'SN_B_lightcurve.dat'))
    at = load_atlas()

    fig, ax = plt.subplots(figsize=(6.6, 3.6))
    cmap = plt.cm.tab10
    for i, (name, df, src) in enumerate(models):
        col = 'V' if np.isfinite(df.get('V', pd.Series(dtype=float))).any() \
            else 'bol_mag'
        tpk = df.loc[df[col].idxmin(), 'time']
        ax.plot(df['time'] - tpk, df[col], lw=0.9, color=cmap(i % 10),
                label=name.replace('_', '\\_'))
    d = obs['r']
    ax.errorbar(d['MJD'] - PEAK_MJD, d['Mabs'], yerr=d['Error'], fmt='o',
                ms=2.4, color='k', elinewidth=0.5, capsize=0,
                label='obs $r$', zorder=10)
    s = at['o']
    sel = (s['MJD'] > 60790) & (s['MJD'] < 61120)
    ax.errorbar(s.loc[sel, 'MJD'] - PEAK_MJD, s.loc[sel, 'Mabs'],
                yerr=s.loc[sel, 'dm'], fmt='s', ms=1.8, color=C_O, alpha=0.6,
                elinewidth=0.4, capsize=0, label='ATLAS $o$', zorder=9)
    ax.axvline(0, color='0.4', lw=0.7, ls='--')
    ax.set_xlabel('days since own peak')
    ax.set_ylabel('absolute magnitude ($V$/bol.)')
    ax.set_xlim(-60, 120)
    ax.set_ylim(-13, -20)
    ax.legend(frameon=False, fontsize=6, ncol=2, loc='lower left')
    save(fig, 'fig_grid.pdf')


# ---------------------------------------------------------------------------
# Fig. 8: observed spectra with line identifications
# ---------------------------------------------------------------------------
def fig_spectra():
    specs = [
        ('SN2024pba_60523.465_PHASE-326.5_Palomar-5.1m-Hale_DBSP.txt',
         'DBSP, MJD 60523.5 ($+10$ d after SN$_{\\rm A}$ peak)'),
        ('SN2024pba_60875.145_PHASE25.1_Nordic-Optical-Telescope_ALFOSC.txt',
         'ALFOSC, MJD 60875.1 ($+25$ d after SN$_{\\rm B}$ peak)'),
    ]
    lines = [3934, 3968, 4102, 4340, 4861, 4924, 5018, 5169, 5876, 6355,
             6563, 6678, 7065, 7774, 8498, 8542, 8662]
    labs = {3934: 'Ca II K', 3968: 'Ca II H', 4861: 'H$\\beta$',
            5169: 'Fe II', 5876: 'He I', 6355: 'Si II', 6563: 'H$\\alpha$',
            7065: 'He I', 7774: 'O I', 8542: 'Ca II'}
    fig, axes = plt.subplots(2, 1, figsize=(6.6, 4.6))
    for ax, (fn, lab) in zip(axes, specs):
        d = np.loadtxt(os.path.join(OBS, fn))
        lam, fl = d[:, 0] / (1 + Z), d[:, 1]
        med = np.nanmedian(fl[(lam > 5000) & (lam < 7000)])
        fl = fl / med
        ax.plot(lam, fl, color='k', lw=0.6)
        lo, hi = np.nanpercentile(fl[(lam > 3600) & (lam < 9000)], [1, 99])
        ax.set_ylim(lo - 0.15 * (hi - lo), hi + 0.45 * (hi - lo))
        ax.set_xlim(3400, 9200)
        ax.set_ylabel('normalised flux')
        ax.text(0.015, 0.06, lab, transform=ax.transAxes, fontsize=8,
                va='bottom')
        ymax = ax.get_ylim()[1]
        for wl in lines:
            ax.axvline(wl, color='0.78', lw=0.5, ls='-', zorder=0)
            if wl in labs:
                ax.text(wl, ymax, labs[wl], rotation=90, fontsize=6,
                        ha='center', va='top', color='0.35')
    axes[1].set_xlabel('rest wavelength [\\AA]')
    fig.align_ylabels(axes)
    save(fig, 'fig_spectra.pdf')


if __name__ == '__main__':
    fig_observed_lc()
    fig_schematic()
    fig_input_structure()
    fig_sn_a_lc()
    fig_model_lc()
    fig_heating()
    fig_grid()
    fig_spectra()
    print('done')
