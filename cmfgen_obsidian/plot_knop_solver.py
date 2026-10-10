#!/usr/bin/env python3
"""
Why the Knop et al. (2009) formal solver is needed for Type IIn
velocity fields.

Left column  - homologous expansion (monotonic): V = V0*r so
               sigma = d ln V / d ln r = 1 everywhere.  The comoving-
               frame frequency shifts monotonically along any ray, so
               every (frequency, line) pair has exactly ONE resonance
               point - the Sobolev approximation the stock solver uses.
Right column - interacting ejecta (Type IIn): the fast shell is
               decelerated into a cold dense shell at r = 5, where
               sigma < -1.  The comoving frequency along a ray folds
               back, so one line resonates at THREE points and the sign
               of gamma = (V/r)(1 + sigma*mu^2) flips inside the
               compression zone.  The Sobolev / stock assumptions fail.
Bottom row   - the Knop fix: chi_hat = chi + |gamma|/dln_nu stays
               positive, the gamma*dI/dln_nu term moves into a source
               term S' = S + Q*I_upwind, and sign(gamma) picks the
               upwind frequency.  The "next" frequency at compressing
               points is taken from the previous great iteration
               (Gauss-Seidel lagged I_NEXT, checkpointed in
               FGINEXT_GS).

Output: attachments/knop_solver.png
"""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

C_KMS = 2.99792458e5  # speed of light, km/s

R_SHOCK = 5.0

r = np.linspace(0.1, 10.0, 8000)
p = 3.5                                   # ray impact parameter
s = np.linspace(-10.0, 10.0, 24000)
r_ray = np.sqrt(p ** 2 + s ** 2)
mu = s / r_ray


def ray_quantities(v_prof):
    """Projected Doppler factor along the ray: nu_cmf = nu_lab*(1-dop)."""
    v_ray = np.interp(r_ray, r, v_prof)
    return v_ray * mu / C_KMS


# ------------------------------------------------------------------
# Left model: homologous expansion, V = V0*r, sigma = 1 everywhere.
# ------------------------------------------------------------------
v_mono = 2.0e3 * r
sigma_mono = np.gradient(np.log(v_mono), np.log(r))
dop_mono = ray_quantities(v_mono)
delta_mono = 0.008
cross_mono = np.where(np.diff(np.sign(dop_mono - delta_mono)) != 0)[0]

# ------------------------------------------------------------------
# Right model: same homologous ejecta, but the shell is decelerated
# into a cold dense shell at r = 5 (one compression region), then
# joins the faster shocked CSM beyond it.
# ------------------------------------------------------------------
sw1 = 0.5 * (1.0 + np.tanh((r - 5.0) / 0.15))   # into the CDS
sw2 = 0.5 * (1.0 + np.tanh((r - 6.2) / 0.30))   # CDS -> shocked CSM
v_out = (1.0 - sw2) * 1.5e3 + sw2 * 4.0e3
v_int = (1.0 - sw1) * v_mono + sw1 * v_out
sigma_int = np.gradient(np.log(v_int), np.log(r))
dop_int = ray_quantities(v_int)

# Detuning inside the fold: between the CDS dip minimum and the
# shocked-CSM level, the line is crossed three times on the receding
# side.
pos = s > 0
g_pos = np.gradient(dop_int[pos], s[pos])
turn = np.where(np.diff(np.sign(g_pos)) != 0)[0]
dop_turn = dop_int[pos][turn]
delta_int = 0.5 * (dop_turn[1] + dop_int[s > 9.0][-1])
cross_int = np.where(np.diff(np.sign(dop_int - delta_int)) != 0)[0]

# ------------------------------------------------------------------
fig = plt.figure(figsize=(11, 12))
gs = fig.add_gridspec(4, 2, height_ratios=[1.0, 0.65, 1.0, 0.75],
                      hspace=0.48, wspace=0.30)
fig.suptitle('Sobolev vs non-monotonic transfer: why the Knop solver '
             'is needed', fontsize=13)

col_titles = ['Homologous expansion', 'Interacting ejecta (Type IIn)']

# ---- row 1: V(r) ---------------------------------------------------
ax = fig.add_subplot(gs[0, 0])
ax.plot(r, v_mono, color='k', lw=1.6)
ax.set_title(col_titles[0], fontsize=11)
ax.set_ylabel(r'$V(r)$ [km s$^{-1}$]')
ax.set_xlim(0, 10)
ax.set_ylim(0, 21000)
ax.text(0.4, 17500, r'$V \propto r$', fontsize=10)

ax = fig.add_subplot(gs[0, 1])
ax.plot(r, v_int, color='k', lw=1.6)
ax.fill_between(r, 0, 10500, where=sigma_int < -1.0, color='tab:red',
                alpha=0.15)
ax.axvline(R_SHOCK, color='gray', lw=0.8, ls=':')
ax.annotate('reverse shock:\nshell stalls into the\ncold dense shell',
            xy=(5.05, 5800), xytext=(0.6, 7800), fontsize=9,
            arrowprops=dict(arrowstyle='->', color='k', lw=0.8))
ax.text(6.6, 4700, 'shocked CSM', fontsize=9)
ax.text(0.4, 3800, 'fast ejecta\nshell', fontsize=9)
ax.set_title(col_titles[1], fontsize=11)
ax.set_xlim(0, 10)
ax.set_ylim(0, 10500)

# ---- row 2: sigma(r) ----------------------------------------------
ax = fig.add_subplot(gs[1, 0])
ax.plot(r, sigma_mono, color='k', lw=1.6)
ax.axhline(-1.0, color='tab:red', lw=0.8, ls=':')
ax.set_ylabel(r'$\sigma = d\ln V/d\ln r$')
ax.set_xlim(0, 10)
ax.set_ylim(-32, 12)
ax.text(0.4, -20, r'$\sigma = 1$: monotonic'
        '\n' r'$\sigma > -1$ always', fontsize=9)

ax = fig.add_subplot(gs[1, 1])
ax.plot(r, sigma_int, color='k', lw=1.6)
ax.axhline(-1.0, color='tab:red', lw=0.8, ls=':')
ax.fill_between(r, -32, 12, where=sigma_int < -1.0, color='tab:red',
                alpha=0.15)
ax.set_xlim(0, 10)
ax.set_ylim(-32, 12)
ax.annotate(r'$\sigma < -1$: compression,'
            '\nsolver breaks here', xy=(5.1, -20), xytext=(6.4, -16),
            fontsize=9, color='tab:red',
            arrowprops=dict(arrowstyle='->', color='tab:red', lw=0.8))

# ---- row 3: v*mu along the ray -------------------------------------
ax = fig.add_subplot(gs[2, 0])
ax.plot(s, dop_mono, color='k', lw=1.6)
ax.axhline(delta_mono, color='tab:green', lw=1.0, ls='--')
for k in cross_mono:
    ax.plot(s[k], dop_mono[k], 'o', ms=9, color='tab:green', zorder=5)
    ax.annotate('one resonance\n(Sobolev point)', (s[k], dop_mono[k]),
                textcoords='offset points', xytext=(-110, -25),
                fontsize=9, color='tab:green',
                arrowprops=dict(arrowstyle='->', color='tab:green',
                                lw=0.8))
ax.set_ylabel(r'$v\mu$  [units of $c$]')
ax.set_xlabel('ray coordinate $s$')
ax.set_xlim(-10, 10)
ax.set_ylim(-0.025, 0.025)
ax.set_title('Along a ray: one line, one resonance\n'
             '-> Sobolev works', fontsize=10)

ax = fig.add_subplot(gs[2, 1])
ax.plot(s, dop_int, color='k', lw=1.6)
ax.axhline(delta_int, color='tab:green', lw=1.0, ls='--')
lbl = [(-30, 14), (-4, -28), (22, 14)]
for i, k in enumerate(cross_int):
    j = min(i, len(lbl) - 1)
    ax.plot(s[k], dop_int[k], 'o', ms=9, color='tab:green', zorder=5)
    ax.annotate(f'root {i + 1}', (s[k], dop_int[k]),
                textcoords='offset points', xytext=lbl[j],
                ha='center', fontsize=9, color='tab:green',
                arrowprops=dict(arrowstyle='-', color='tab:green',
                                lw=0.8))
ax.set_xlabel('ray coordinate $s$')
ax.set_xlim(-10, 10)
ax.set_ylim(-0.025, 0.025)
ax.set_title('Same ray, folded profile: one line, '
             f'{len(cross_int)} resonances\n+ flipped coupling '
             '-> Sobolev fails', fontsize=10)

# ---- row 4: the Knop fix -------------------------------------------
ax = fig.add_subplot(gs[3, :])
ax.axis('off')
ax.set_xlim(0, 10)
ax.set_ylim(0, 1.2)

ax.text(0.1, 1.02, 'The Knop et al. (2009) fix', fontsize=11,
        fontweight='bold')
ax.text(0.1, 0.86, r'$\chi_{\rm hat} = \chi + |\gamma|/d\ln\nu$'
        '  (generalized opacity, always positive)', fontsize=10)
ax.text(0.1, 0.68, r"$S' = S + Q\,I_{\rm upwind}$"
        r'$,  Q = (\gamma/d\ln\nu)/\chi_{\rm hat}$'
        '  (velocity coupling moved into the source)', fontsize=10)
ax.text(0.1, 0.50, 'sign(' r'$\gamma$' ') picks the upwind frequency:',
        fontsize=10)
ax.text(0.35, 0.32, r'$\gamma > 0$ expanding:  use previous $\nu$',
        fontsize=10, color='tab:blue')
ax.text(0.35, 0.16, r'$\gamma < 0$ compressing:  use next $\nu$'
        '  (from previous iteration - lagged I_NEXT)',
        fontsize=10, color='tab:red')

# Upwinding sketch on the right half: frequency grid with arrows.
x0, x1, x2 = 6.6, 7.6, 8.6
y_line = 0.52
ax.plot([x0 - 0.5, x2 + 0.5], [y_line, y_line], color='k', lw=1.2)
for x, lab in [(x0, r'$\nu_{i-1}$'), (x1, r'$\nu_i$'),
               (x2, r'$\nu_{i+1}$')]:
    ax.plot([x, x], [y_line - 0.05, y_line + 0.05], color='k', lw=1.2)
    ax.text(x, y_line - 0.14, lab, ha='center', fontsize=10)
ax.annotate('', xy=(x1, y_line + 0.10), xytext=(x0, y_line + 0.10),
            arrowprops=dict(arrowstyle='->', color='tab:blue', lw=1.6))
ax.text((x0 + x1) / 2, y_line + 0.16, r'$\gamma > 0$',
        ha='center', fontsize=9, color='tab:blue')
ax.annotate('', xy=(x1, y_line + 0.34), xytext=(x2, y_line + 0.34),
            arrowprops=dict(arrowstyle='->', color='tab:red', lw=1.6))
ax.text((x1 + x2) / 2, y_line + 0.40,
        r'$\gamma < 0$: lagged $I_{\rm next}$',
        ha='center', fontsize=9, color='tab:red')
ax.text(7.6, 0.02, 'intensity always upwinds from the correct side',
        ha='center', fontsize=9, style='italic')

out = Path(__file__).resolve().parent / 'attachments' / 'knop_solver.png'
plt.savefig(out, dpi=160, bbox_inches='tight')
print(f'wrote {out}  (roots: monotonic={len(cross_mono)}, '
      f'interacting={len(cross_int)})')
