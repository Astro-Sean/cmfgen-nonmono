#!/usr/bin/env python3
"""Estimate the SN_B explosion epoch from the ~M=-15.5 shock-cooling bump.

Two wind-embedded-shock constraints are combined:

1) Bump-slope fit.  The r-band bump (Mabs ~ -15.5, MJD 60810-60831)
   slowly brightens.  A shock dissipating in an optically thick wind
   leaks radiation through a declining column; to first order
       L(t) = L_p * (1 - exp(-(t - t0)/tau))
   where tau is the diffusion time of the shocked shell+wind mass
   (Arnett t_m for the assumed shell mass/velocity) and t0 the launch.
   With tau fixed per (M, v) we fit (t0, L_p).

2) Crossing-time inversion.  The steep main-peak rise starts when the
   forward shock reaches the wind-ejecta buffer at r = v_inner*t_age.
   With the shock propagating at ~v_shell:
       t_rise = t0 + v_inner*(t0 - t_SNA)/v_shell
       =>  t0 = (t_rise + k*t_SNA)/(1 + k),   k = v_inner/v_shell

Both are evaluated across the assumed (M_shell, v_shell) grid to give
the span of explosion times and the implied SN_A ejecta age at launch
(= source_age_days for the X_* models).
"""

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

PEAK_MJD = 60850.0
SN_A_EXP = 60490.0            # SN_A explosion epoch
T_RISE   = 60833.0            # onset of steep main-peak rise (~ +/-3 d)
DM = 34.526                   # distance modulus, DL = 80.4 Mpc
A_R = 2.285 * 0.15
V_INNER = 2370.0              # km/s, SN_A ejecta inner edge = r_ej_inner/t_age
BETA, C, MSUN = 13.7, 2.998e10, 1.989e33

OBS = '/ptmp/sbrennan/SN2024pba/observational_data/SN_B_lightcurve.dat'

# ---- bump photometry -------------------------------------------------
raw = pd.read_csv(OBS)
r = raw[raw['Filter'] == 'r'].copy()
det = r[r['Mag'].apply(lambda x: x != '-' and not pd.isna(x))].copy()
det['MJD'] = det['MJD'].astype(float)
det['Mag'] = det['Mag'].astype(float)
det['Error'] = det['Error'].astype(float)
bump = det[(det['MJD'] >= 60805) & (det['MJD'] <= 60832)].copy()
bump['Mabs'] = bump['Mag'] - DM - A_R
# fit in linear luminosity proxy (10^(-0.4 M))
bump['Lrel'] = 10.0**(-0.4 * bump['Mabs'])
bump['Lrel_err'] = bump['Lrel'] * 0.4 * np.log(10) * bump['Error']

def t_arnett(M_msun, v_kms, kap=0.34):
    return np.sqrt(2 * kap * M_msun * MSUN / (BETA * C * v_kms * 1e5)) / 86400.0

def bump_model(t, t0, Lp, tau):
    return Lp * (1.0 - np.exp(-(t - t0) / tau)) * (t > t0)

print('=== Wind-embedded shock bump fit ===')
print(f'bump points: {len(bump)} r-band det, Mabs '
      f'{bump["Mabs"].min():.2f}..{bump["Mabs"].max():.2f}, '
      f'MJD {bump["MJD"].min():.1f}..{bump["MJD"].max():.1f}\n')

results = []
print(f'{"M":>5} {"v":>6} {"tau_d":>6} {"t0_fit":>8} {"age":>5}   '
      f'{"t0_cross":>8} {"age":>5}')
for M in (5.0, 7.0, 10.0):
    for v in (6000.0, 8000.0, 10000.0):
        tau = t_arnett(M, v)
        # slope fit: t0 and plateau amplitude free, tau fixed
        try:
            p, _ = curve_fit(
                lambda t, t0, Lp: bump_model(t, t0, Lp, tau),
                bump['MJD'], bump['Lrel'], sigma=bump['Lrel_err'],
                p0=[60790.0, 2.0e6], bounds=([60600, 0], [60820, 1e8]))
            t0_fit = p[0]
        except Exception:
            t0_fit = np.nan
        # crossing inversion
        k = V_INNER / v
        t0_x = (T_RISE + k * SN_A_EXP) / (1.0 + k)
        results.append((M, v, tau, t0_fit, t0_x))
        print(f'{M:5.1f} {v:6.0f} {tau:6.1f} {t0_fit:8.1f} '
              f'{t0_fit - SN_A_EXP:5.0f}   {t0_x:8.1f} {t0_x - SN_A_EXP:5.0f}')

t0s = np.array([r[3] for r in results] + [r[4] for r in results])
t0s = t0s[np.isfinite(t0s)]
print(f'\nCombined explosion-time span: MJD {t0s.min():.0f} - '
      f'{t0s.max():.0f}  (peak {-60850.0 + t0s.min():.0f} to '
      f'{-60850.0 + t0s.max():.0f} d before SN_B max)')
print(f'Implied SN_A age at launch: {t0s.min()-SN_A_EXP:.0f} - '
      f'{t0s.max()-SN_A_EXP:.0f} d  (source_age_days)')

# ---- figure ----------------------------------------------------------
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

A_G, A_O = 3.303 * 0.15, 2.40 * 0.15
ATLAS = '/ptmp/sbrennan/SN2024pba/observational_data/SN2024pba_atclean_cleaned.csv'
OUT = '/ptmp/sbrennan/SN2024pba/SN_B_v2_T6k_grid/snb_explosion_estimate.png'

# context photometry: r/g detections + nightly-binned ATLAS o
ctx = []
for b, A, mk, cl in (('r', A_R, 'o', 'tab:red'), ('g', A_G, 's', 'tab:green')):
    s = raw[raw['Filter'] == b].copy()
    d = s[s['Mag'].apply(lambda x: x != '-' and not pd.isna(x))].copy()
    d['Mag'] = d['Mag'].astype(float); d['Error'] = d['Error'].astype(float)
    d['Mabs'] = d['Mag'] - DM - A
    d['Mabs_err'] = d['Error']
    ctx.append((d, mk, cl, f'{b} (PSF)'))
at = pd.read_csv(ATLAS)
o = at[at['F'] == 'o'].copy()
o['sn'] = o['uJy'] / o['duJy_new']
good = o['Mask'].isin(['0x400800', '0x400900']) & (o['chi/N'].isna() | (o['chi/N'] < 5.0))
o = o[(o['sn'] >= 3) & good].copy()
o['Mabs'] = o['m'] - DM - A_O
o['Mabs_err'] = o['dm'].where(o['dm'].notna(), 1.086 / o['sn'])
rows = []
for night, g in o.groupby(np.floor(o['MJD'])):
    w = 1.0 / g['Mabs_err']**2
    rows.append({'MJD': np.average(g['MJD'], weights=w),
                 'Mabs': np.average(g['Mabs'], weights=w),
                 'Mabs_err': 1.0 / np.sqrt(w.sum())})
ctx.append((pd.DataFrame(rows), '^', 'darkorange', 'o (ATLAS)'))

t0_lo, t0_hi = t0s.min(), t0s.max()
colors = {5.0: 'tab:blue', 7.0: 'tab:purple', 10.0: 'tab:brown'}
fig, (ax, ax2) = plt.subplots(1, 2, figsize=(15.5, 6.2),
                              gridspec_kw={'width_ratios': [2.1, 1]})

for d, mk, cl, lb in ctx:
    ax.errorbar(d['MJD'], d['Mabs'], yerr=d['Mabs_err'], fmt=mk,
                color=cl, ms=4.5, mec='k', mew=0.6, capsize=2,
                elinewidth=0.6, alpha=0.85, zorder=10, label=lb)

# fitted bump rise curves, drawn from t0 to the main-peak rise onset
tt = np.linspace(60600, 60840, 800)
for M, v, tau, t0_fit, t0_x in results:
    Lp = np.median(bump['Lrel']) / (1 - np.exp(-(T_RISE - t0_fit) / tau))
    mabs = -2.5 * np.log10(bump_model(tt, t0_fit, Lp, tau))
    mabs[tt < t0_fit] = np.nan
    mabs[tt > T_RISE] = np.nan
    ax.plot(tt, mabs, color=colors[M], lw=1.3, alpha=0.85,
            ls={6000.0: ':', 8000.0: '--', 10000.0: '-'}[v])

ax.axvspan(t0_lo, t0_hi, color='0.6', alpha=0.3, zorder=0,
           label=f'explosion window\n{t0_lo:.0f}-{t0_hi:.0f}')
for t0_x in {r[4] for r in results}:
    ax.axvline(t0_x, color='0.55', lw=0.8, ls=':', alpha=0.8)
ax.axvline(PEAK_MJD, color='k', lw=1.0, ls='-', alpha=0.6)
ax.text(PEAK_MJD + 0.6, -19.0, 'SN_B peak', rotation=90, fontsize=8)
ax.axvline(60805.0, color='k', lw=1.2, ls='--', alpha=0.8)
ax.text(60805.6, -14.4, 'adopted t0 = 60805', rotation=90, fontsize=8)
ax.annotate('main-peak rise\nbegins ~60833', xy=(T_RISE, -16.2),
            xytext=(60790, -17.1), fontsize=8,
            arrowprops=dict(arrowstyle='->', lw=0.8))
ax.set_xlim(60680, 60870)
ax.set_ylim(-19.4, -14.2)
ax.invert_yaxis()
ax.set_xlabel('MJD')
ax.set_ylabel('Absolute magnitude')
ax.grid(alpha=0.3)
mass_h = [Line2D([], [], color=colors[M], lw=1.5,
                 label=f'M$_{{\\rm shell}}$ = {M:g} M$_\\odot$')
          for M in colors]
vel_h = [Line2D([], [], color='0.3', lw=1.5, ls=s,
                label=f'v$_{{\\rm shell}}$ = {v:g} km s$^{{-1}}$')
         for v, s in ((10000.0, '-'), (8000.0, '--'), (6000.0, ':'))]
leg1 = ax.legend(loc='lower left', fontsize=8, frameon=False)
ax.add_artist(leg1)
ax.legend(handles=mass_h + vel_h, loc='upper left', fontsize=7.5,
          ncol=2, frameon=False)

# right panel: explosion epoch vs v_shell for both methods
for M, v, tau, t0_fit, t0_x in results:
    ax2.plot(v, t0_fit, 'o', color=colors[M], ms=6, mec='k', mew=0.5)
    ax2.plot(v, t0_x, 'v', color=colors[M], ms=6, mec='k', mew=0.5)
ax2.axhspan(t0_lo, t0_hi, color='0.6', alpha=0.3)
ax2.axhline(SN_A_EXP + 320, color='0.3', ls='--', lw=1)
ax2.text(6050, SN_A_EXP + 320 + 1.5, 'sa320 run', fontsize=8)
ax2.set_xlabel('v$_{\\rm shell}$ (km s$^{-1}$)')
ax2.set_ylabel('Explosion MJD')
ax2.set_xticks([6000, 8000, 10000])
ax2.grid(alpha=0.3)
meth_h = [Line2D([], [], marker='o', ls='', color='0.3', mec='k',
                 label='bump-slope fit'),
          Line2D([], [], marker='v', ls='', color='0.3', mec='k',
                 label='crossing inversion')]
ax2.legend(handles=meth_h, loc='lower right', fontsize=8, frameon=False)
ax2.set_title('SN_B explosion estimates', fontsize=10)

fig.suptitle('SN 2024pba — SN_B explosion epoch from the wind-embedded '
             'shock bump', fontsize=11)
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig(OUT, dpi=150, bbox_inches='tight')
print(f'Saved: {OUT}')
