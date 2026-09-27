#!/usr/bin/env python3
"""Diagrams for latex/cmfgen_knop_solution.tex.

Pedagogical figures explaining (a) why the original CMFGEN formal solver
required a monotonic velocity field, and (b) how the Knop et al. (2009)
non-monotonic formal solver works.  Pure schematic figures -- no model data.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ROOT = '/viper/ptmp1/sbrennan/SN2024pba'
FIGDIR = os.path.join(ROOT, 'latex', 'figures')
os.makedirs(FIGDIR, exist_ok=True)

plt.rcParams.update({
    'font.size': 9, 'axes.labelsize': 9, 'axes.titlesize': 9,
    'legend.fontsize': 8, 'xtick.labelsize': 8, 'ytick.labelsize': 8,
    'axes.linewidth': 0.8, 'lines.linewidth': 1.4,
    'figure.dpi': 150, 'savefig.bbox': 'tight',
})


def save(fig, name):
    out = os.path.join(FIGDIR, name)
    fig.savefig(out)
    plt.close(fig)
    print('wrote', out)


# ---------------------------------------------------------------------------
# Fig 1: velocity profile + comoving-frequency march along a ray,
#        monotonic vs non-monotonic (2x2)
# ---------------------------------------------------------------------------
def fig_vel_crossings():
    r = np.linspace(0.0, 1.0, 400)

    v_mono = 0.15 + 0.85 * r                              # monotonic outflow
    # non-monotonic: outward rise, compression dip, outward again
    v_nm = v_mono - 0.28 * np.exp(-((r - 0.55) / 0.07) ** 2)
    v_nm = np.clip(v_nm, 0.0, None)

    nu_obs = 1.0
    nu_cmf_mono = nu_obs * (1.0 - 0.35 * v_mono)          # ~1 - mu v/c
    nu_cmf_nm = nu_obs * (1.0 - 0.35 * v_nm)
    nu0 = 0.85                                            # line rest freq

    fig, axes = plt.subplots(2, 2, figsize=(6.8, 5.2), sharex='col')

    for ax in axes[0]:
        ax.set_ylabel('$v$ (arb.)')
    for ax in axes[1]:
        ax.set_ylabel(r'$\nu_{\rm CMF}$ along ray')
        ax.set_xlabel('path coordinate $s$')
    for ax in axes.flat:
        ax.set_xticks([]); ax.set_yticks([])

    # --- monotonic column ---
    ax = axes[0, 0]
    ax.plot(r, v_mono, 'k-')
    ax.annotate('$dv/dr>0$ everywhere\n(monotonic expansion)',
                xy=(0.55, 0.30), fontsize=8, ha='center')
    ax.set_ylim(0, 1.1)
    ax.set_title('Original solver regime', fontsize=9)

    ax = axes[1, 0]
    ax.plot(r, nu_cmf_mono, 'k-')
    ax.axhline(nu0, color='#b2182b', ls='--', lw=1.0)
    ix = np.argmin(np.abs(nu_cmf_mono - nu0))
    ax.plot(r[ix], nu_cmf_mono[ix], 'o', ms=6, mfc='none', mec='#b2182b',
            mew=1.5)
    ax.annotate('one resonance\n(one Sobolev point)',
                xy=(r[ix], nu0), xytext=(0.10, 0.62), fontsize=8,
                arrowprops=dict(arrowstyle='->', lw=0.6, color='0.3'))
    ax.text(0.97, nu0 + 0.012, r'$\nu_0$ (line rest frame)', fontsize=7.5,
            color='#b2182b', ha='right', va='bottom')
    ax.annotate('', xy=(0.86, nu_cmf_mono[344]),
                xytext=(0.70, nu_cmf_mono[280]),
                arrowprops=dict(arrowstyle='->', lw=0.7, color='0.3'))
    ax.text(0.75, 0.585, 'decreasing\n(redshift)', fontsize=7,
            ha='center', va='bottom', color='0.3')
    ax.set_ylim(0.55, 1.02)

    # --- non-monotonic column ---
    ax = axes[0, 1]
    ax.plot(r, v_nm, 'k-')
    ax.axvspan(0.45, 0.65, color='#b2182b', alpha=0.12)
    ax.annotate('compression zone\n$dv/dr<0$ ($\\sigma<-1$)',
                xy=(0.55, v_nm.min() + 0.03), xytext=(0.14, 0.14),
                fontsize=8,
                arrowprops=dict(arrowstyle='->', lw=0.6, color='0.3'))
    ax.set_ylim(0, 1.1)
    ax.set_title('Interacting-SN regime', fontsize=9)

    ax = axes[1, 1]
    ax.plot(r, nu_cmf_nm, 'k-')
    ax.axhline(nu0, color='#b2182b', ls='--', lw=1.0)
    ax.axvspan(0.45, 0.65, color='#b2182b', alpha=0.12)
    cross = np.where(np.diff(np.sign(nu_cmf_nm - nu0)))[0]
    offsets = [(-0.055, 0.05), (-0.015, -0.115), (0.025, 0.05)]
    for i, j in enumerate(cross):
        ax.plot(r[j], nu_cmf_nm[j], 'o', ms=6, mfc='none', mec='#b2182b',
                mew=1.5)
        dx, dy = offsets[i]
        ax.annotate(f'{i+1}', xy=(r[j], nu_cmf_nm[j]),
                    xytext=(r[j] + dx, nu_cmf_nm[j] + dy),
                    fontsize=8, color='#b2182b')
    ax.text(0.97, nu0 - 0.018, r'$\nu_0$', fontsize=7.5, color='#b2182b',
            ha='right', va='top')
    ax.annotate('blueshift here\n(wrong direction!)',
                xy=(0.575, nu_cmf_nm[230] + 0.005), xytext=(0.75, 0.955),
                fontsize=8, ha='left',
                arrowprops=dict(arrowstyle='->', lw=0.6, color='0.3'))
    ax.text(0.28, 0.575, 'three resonances\n(three Sobolev points)',
            fontsize=8, ha='center', color='0.25')
    ax.set_ylim(0.55, 1.02)

    fig.align_ylabels(axes)
    fig.tight_layout(h_pad=0.4)
    save(fig, 'fig_knop_crossings.pdf')


# ---------------------------------------------------------------------------
# Fig 2: upwind-frequency choice on the frequency grid
# ---------------------------------------------------------------------------
def fig_upwind():
    fig, axes = plt.subplots(1, 2, figsize=(6.8, 2.5), sharey=True)

    nus = np.linspace(0, 1, 9)
    for ax, gamma_neg, title in (
            (axes[0], False, 'expanding: $\\Gamma>0$'),
            (axes[1], True, 'compressing: $\\Gamma<0$')):
        ax.plot(nus, np.zeros_like(nus), 'k-', lw=1.2)
        for i, x in enumerate(nus):
            ax.plot([x, x], [-0.06, 0.06], 'k-', lw=1.0)
            ax.text(x, 0.11, f'$\\nu_{{{9-i}}}$', ha='center', fontsize=7)
        ax.text(0.0, 0.32, 'high $\\nu$', fontsize=7.5, ha='left',
                color='0.4')
        ax.text(1.0, 0.32, 'low $\\nu$', fontsize=7.5, ha='right',
                color='0.4')
        # sweep order = decreasing frequency = left to right on this axis
        ax.add_patch(FancyArrowPatch((0.15, 0.23), (0.85, 0.23),
                                     arrowstyle='-|>', mutation_scale=9,
                                     color='0.5', lw=0.9))
        ax.text(0.5, 0.25, 'sweep order', fontsize=7, ha='center',
                color='0.5')

        i_c = 4   # current point
        ax.plot(nus[i_c], 0, 'o', ms=8, mfc='#1f78b4', mec='k')
        ax.text(nus[i_c] - 0.07, -0.20, 'current\npoint', ha='center',
                fontsize=8, color='#1f78b4')
        if not gamma_neg:
            i_up = i_c - 1    # previous (higher nu, already computed)
            lab = 'upwind = previous\n(already computed)'
            xtxt = nus[i_up] - 0.02
        else:
            i_up = i_c + 1    # next (lower nu, NOT yet computed)
            lab = 'upwind = next\n(not yet computed!)'
            xtxt = nus[i_up] + 0.04
        ax.plot(nus[i_up], 0, 's', ms=8, mfc='#e66101', mec='k')
        ax.annotate(lab, xy=(nus[i_up], -0.09),
                    xytext=(xtxt, -0.45), fontsize=8, ha='center',
                    color='#e66101',
                    arrowprops=dict(arrowstyle='->', lw=0.7,
                                    color='#e66101'))
        ax.set_title(title, fontsize=9)
        ax.set_xlim(-0.12, 1.12)
        ax.set_ylim(-0.62, 0.44)
        ax.axis('off')
    fig.tight_layout()
    save(fig, 'fig_knop_upwind.pdf')


# ---------------------------------------------------------------------------
# Fig 3: solver flowchart — original vs Knop non-monotonic path
# ---------------------------------------------------------------------------
def fig_flowchart():
    fig, ax = plt.subplots(figsize=(7.0, 4.8))
    ax.axis('off')
    ax.set_xlim(-3.8, 10.2)
    ax.set_ylim(0, 10)

    def box(x, y, w, h, text, fc='#e0f3f8', ec='0.3', fs=7.5):
        ax.add_patch(FancyBboxPatch((x, y), w, h,
                                    boxstyle='round,pad=0.08',
                                    fc=fc, ec=ec, lw=0.8))
        ax.text(x + w / 2, y + h / 2, text, ha='center', va='center',
                fontsize=fs)

    def arrow(x0, y0, x1, y1, color='0.3'):
        ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1),
                                     arrowstyle='-|>', mutation_scale=9,
                                     color=color, lw=0.9))

    ax.text(2.6, 9.65, 'Original CMFGEN solver', fontsize=8.5,
            fontweight='bold', ha='center')
    ax.text(7.4, 9.65, 'Knop et al. (2009) solver', fontsize=8.5,
            fontweight='bold', ha='center')

    # left column
    box(0.6, 8.0, 4.0, 1.0,
        'velocity field $v(r)$\nassume $dv/dr \\geq 0$')
    box(0.6, 6.2, 4.0, 1.2,
        'generalised opacity\n$\\hat\\chi = \\chi + '
        '\\Gamma/\\Delta\\ln\\nu$')
    box(0.6, 4.4, 4.0, 1.2,
        'upwind = previous $\\nu$\n(sweep high $\\rightarrow$ low $\\nu$)')
    box(0.6, 2.6, 4.0, 1.2,
        'short characteristics\n$dI/d\\hat\\tau = I - S\'$')
    box(0.6, 0.8, 4.0, 1.2,
        'moments $J_\\nu$ $\\rightarrow$ rates\n$\\rightarrow$ next '
        'iteration')

    arrow(2.6, 8.0, 2.6, 7.4)
    arrow(2.6, 6.2, 2.6, 5.6)
    arrow(2.6, 4.4, 2.6, 3.8)
    arrow(2.6, 2.6, 2.6, 2.0)

    # failure branch (red, to the left)
    box(-3.7, 5.9, 3.0, 1.8,
        'if $\\Gamma<0$ (compression):\n'
        '$\\hat\\chi<0$ $\\rightarrow$ $e^{+\\tau}$\n'
        'blow-up, wrong upwind,\n'
        'negative depths, NaNs',
        fc='#fde0dd', ec='#b2182b', fs=7)
    arrow(0.6, 6.8, -0.65, 6.8, color='#b2182b')

    # right column
    box(5.4, 8.0, 4.0, 1.0,
        'same $v(r)$, compression allowed\n'
        '(true $\\sigma$ kept in SIGMA_TRUE)')
    box(5.4, 6.2, 4.0, 1.2,
        'opacity uses magnitude only\n$\\hat\\chi = \\chi + '
        '|\\Gamma|/\\Delta\\ln\\nu > 0$')
    box(5.4, 4.4, 4.0, 1.2,
        'upwind by sign of $\\Gamma$:\n'
        '$\\Gamma\\geq0$: previous $\\nu$; $\\Gamma<0$: next $\\nu$')
    box(5.4, 2.6, 4.0, 1.2,
        'piecewise short characteristics;\n$I_{\\rm NEXT}$ stored from '
        'previous\niteration (Gauss--Seidel)')
    box(5.4, 0.8, 4.0, 1.2,
        'moments $J_\\nu$ $\\rightarrow$ rates\n$\\rightarrow$ next '
        'great iteration\n(converges over several rounds)', fs=7)

    arrow(7.4, 8.0, 7.4, 7.4)
    arrow(7.4, 6.2, 7.4, 5.6)
    arrow(7.4, 4.4, 7.4, 3.8)
    arrow(7.4, 2.6, 7.4, 2.0)
    ax.add_patch(FancyArrowPatch((9.55, 1.4), (9.55, 3.2),
                                 arrowstyle='-|>', mutation_scale=9,
                                 color='#b2182b', lw=0.9))
    ax.text(9.62, 2.3, '$I_{\\rm NEXT}$\nstored', fontsize=6.5,
            color='#b2182b', va='center')

    save(fig, 'fig_knop_flowchart.pdf')


# ---------------------------------------------------------------------------
# Fig 4: the river analogy — three panels
# ---------------------------------------------------------------------------
def fig_river():
    from matplotlib.patches import Rectangle, Circle

    fig, axes = plt.subplots(3, 1, figsize=(6.8, 6.8))

    town_x = [1.3, 3.3, 5.3, 7.3, 9.0]
    town_lab = [r'$\nu_{i+2}$', r'$\nu_{i+1}$', r'$\nu_i$',
                r'$\nu_{i-1}$', r'$\nu_{i-2}$']
    gorge = (4.3, 6.3)          # reversed-current stretch
    rng = np.random.default_rng(7)
    leaf_x = rng.uniform(0.6, 9.6, 26)
    leaf_y = rng.uniform(0.62, 1.10, 26)

    def draw_river(ax, title, reversed_flow):
        ax.set_xlim(0, 10); ax.set_ylim(0, 3.0); ax.axis('off')
        # water
        ax.add_patch(Rectangle((0.3, 0.55), 9.4, 1.0, fc='#d6eaf8',
                               ec='0.4', lw=0.8, zorder=0))
        if reversed_flow:
            ax.add_patch(Rectangle((gorge[0], 0.55), gorge[1] - gorge[0],
                                   1.0, fc='#f5b7b1', ec='none',
                                   alpha=0.55, zorder=1))
            ax.text((gorge[0] + gorge[1]) / 2, 1.62,
                    'reversed current', ha='center', va='bottom',
                    fontsize=7, color='#922b21')
        # leaves (photons)
        for x, y in zip(leaf_x, leaf_y):
            ax.add_patch(Circle((x, y), 0.045, fc='#27ae60',
                                ec='none', zorder=2))
        # current arrows inside the river (high nu -> low nu is L->R)
        for x0 in np.arange(0.7, 9.6, 0.9):
            mid = x0 + 0.35
            rev = reversed_flow and gorge[0] < mid < gorge[1]
            dx = -0.55 if rev else 0.55
            ax.add_patch(FancyArrowPatch((x0, 1.05), (x0 + dx, 1.05),
                                         arrowstyle='-|>',
                                         mutation_scale=7,
                                         color='#21618c', lw=0.8,
                                         zorder=3))
        # towns on the river
        for j, (x, lab) in enumerate(zip(town_x, town_lab)):
            hl = (j == 2)
            ax.add_patch(Rectangle((x - 0.22, 0.87), 0.44, 0.36,
                                   fc='#fdebd0' if hl else 'white',
                                   ec='#b2182b' if hl else 'k',
                                   lw=1.2 if hl else 0.8, zorder=4))
            ax.text(x, 0.83, lab, ha='center', va='top', fontsize=6.5,
                    zorder=5)
        ax.text(0.3, 2.72, title, fontsize=8.5, ha='left', va='center',
                fontweight='bold')
        return ax

    def ask_arrow(ax, x_from, x_to, color, ls='-'):
        ax.add_patch(FancyArrowPatch((x_from, 1.85), (x_to, 1.85),
                                     arrowstyle='-|>', mutation_scale=9,
                                     color=color, lw=1.0, ls=ls,
                                     zorder=6))

    xi = town_x[2]
    x_prev = town_x[1]     # higher nu, left
    x_next = town_x[3]     # lower nu, right

    # --- panel 1: expanding ---
    ax = axes[0]
    draw_river(ax, 'Expanding ($\\Gamma>0$): normal stretch', False)
    ask_arrow(ax, xi - 0.15, x_prev + 0.25, '#1f78b4')
    ax.text((xi + x_prev) / 2, 1.97, 'asks upstream\n(already visited)',
            ha='center', fontsize=7, color='#1f78b4')
    ax.text(5.0, 0.22, 'leaves always arrive from upstream (high $\\nu$)',
            ha='center', fontsize=7, color='0.35')

    # --- panel 2: compression, old solver fails ---
    ax = axes[1]
    draw_river(ax, 'Compressing ($\\Gamma<0$): the old solver breaks',
               True)
    ask_arrow(ax, xi - 0.15, x_prev + 0.25, '#b2182b', ls=':')
    ax.plot((xi + x_prev) / 2, 1.85, 'x', ms=7, color='#b2182b', mew=2)
    ax.text((xi + x_prev) / 2, 1.97, 'old solver still\nasks upstream',
            ha='center', fontsize=7, color='#b2182b')
    ask_arrow(ax, x_next - 0.25, xi + 0.15, '#27ae60')
    ax.text((xi + x_next) / 2, 1.97, 'but leaves actually\n'
            'come from downstream', ha='center', fontsize=7,
            color='#1e8449')
    ax.text(5.0, 0.22, 'signed current also gives negative muddiness: '
            'leaves multiply, intensities blow up',
            ha='center', fontsize=7, color='0.35')

    # --- panel 3: Knop fix ---
    ax = axes[2]
    draw_river(ax, 'Compressing ($\\Gamma<0$): the Knop fix', True)
    ask_arrow(ax, xi + 0.15, x_next - 0.25, '#e66101')
    ax.text((xi + x_next) / 2, 1.97, 'asks downstream: $I_{\\rm NEXT}$\n'
            'remembered from last round', ha='center', fontsize=7,
            color='#e66101')
    ax.text(5.0, 0.22, 'muddiness uses $|\\Gamma|$ only: always positive, '
            'always attenuates', ha='center', fontsize=7, color='0.35')

    fig.text(0.5, 0.012, 'towns = frequency-grid points,  '
             'leaves = photons ($I_\\nu$),  current = velocity coupling '
             '$\\Gamma$,  muddiness = opacity $\\hat\\chi$',
             ha='center', fontsize=7.5, color='0.3')
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    save(fig, 'fig_knop_river.pdf')


if __name__ == '__main__':
    fig_vel_crossings()
    fig_upwind()
    fig_flowchart()
    fig_river()
    print('done')
