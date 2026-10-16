#!/usr/bin/env python3
"""A&A-style plot of the SN_B_X_grid parameter study vs SN 2024pba.

Four band panels (g, r, i, ATLAS o) in a 2x2 layout sharing the phase
axis used for the Dessart-analog grid figure, plus a companion
g-r / o-g / o-r colour figure.  Every model in the grid is a
single-parameter variation of the fiducial model X, so each curve is
labelled by its changed parameter and grouped visually by parameter
family.  The same -19 d r-band foot anchor as
plot_dessart_grid_vs_obs.py sets the per-model time shift.

Usage:
    python scripts/plot_x_grid_vs_obs.py [--grid-dir SN_B_X_grid]
"""

import argparse
import glob
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.ticker import AutoMinorLocator
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import plot_dessart_grid_vs_obs as base

GRID_DIR = '/ptmp/sbrennan/SN2024pba/SN_B_X_grid'
SN2009IP_FILE = '/ptmp/sbrennan/SN2024pba/observational_data/SN2009ip.csv'
# 2012B g/r maximum; on this phase axis the 2012A precursor sits near
# -40 d, the same window as the observed 2024pba r-band foot.
SN2009IP_PEAK_MJD = 56206.0
SN2009IP_BIN_D = 2.0

# A&A single column is ~9 cm (3.54 in) wide.
AA_COL_IN = 3.54
plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['DejaVu Serif'],
    'mathtext.fontset': 'dejavuserif',
    'font.size': 8,
    'axes.labelsize': 8,
    'axes.titlesize': 8,
    'legend.fontsize': 6,
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    'axes.linewidth': 0.7,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.top': True,
    'ytick.right': True,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})

# Curve style and legend label per model.  Families share a colour so
# the panel reads as a one-parameter-at-a-time study; the fiducial X
# is the thick black reference.  Models absent from this table fall
# back to the directory name so new grid members still plot.
MODEL_STYLE = {
    'X': dict(color='k', lw=1.8, ls='-', zorder=5,
              label=r'optimum ($M_{\rm sh}=4$, $v_{\rm sh}=7000$, '
                   r'wall $6\times10^{17}$)'),
    'X_M4': dict(color='tab:blue', lw=1.1, ls='-', zorder=12,
                 label=r'$M_{\rm sh}=4\,M_\odot$'),
    'X_M10': dict(color='tab:blue', lw=1.1, ls='--', zorder=12,
                  label=r'$M_{\rm sh}=10\,M_\odot$'),
    'X_V6000': dict(color='tab:red', lw=1.1, ls='-', zorder=12,
                    label=r'$v_{\rm sh}=6000$ km s$^{-1}$'),
    'X_V8000': dict(color='tab:red', lw=1.1, ls='--', zorder=12,
                    label=r'$v_{\rm sh}=8000$ km s$^{-1}$'),
    'X_V10000': dict(color='tab:red', lw=1.1, ls='-.', zorder=12,
                     label=r'$v_{\rm sh}=10\,000$ km s$^{-1}$'),
    'X_V12000': dict(color='tab:red', lw=1.1, ls=':', zorder=12,
                     label=r'$v_{\rm sh}=12\,000$ km s$^{-1}$'),
    'X_pl0': dict(color='tab:green', lw=1.1, ls='-', zorder=12,
                  label=r'$n=0.0$'),
    'X_pl05': dict(color='tab:green', lw=1.1, ls=':', zorder=12,
                   label=r'$n=0.5$'),
    'X_pl1': dict(color='tab:green', lw=1.1, ls='--', zorder=12,
                  label=r'$n=1.0$'),
    'X_pl2': dict(color='tab:green', lw=1.1, ls='-.', zorder=12,
                  label=r'$n=2.0$'),
    'X_w3e17': dict(color='tab:purple', lw=1.1, ls='-', zorder=12,
                    label=r'$w=3\times10^{17}$ g cm$^{-1}$'),
    'X_wind': dict(color='tab:orange', lw=1.1, ls='-', zorder=12,
                   label=(r'no ejecta, CSM trunc.\ at '
                          r'$6\times10^{15}$ cm')),
    'X_windext': dict(color='tab:orange', lw=1.1, ls='--', zorder=12,
                      label=(r'no ejecta, CSM extended to '
                             r'$3.2\times10^{16}$ cm')),
    'X_M4_rs5w3': dict(color='tab:cyan', lw=1.1, ls='-', zorder=12,
                       label=(r'$M_{\rm sh}=4$, $w_{\rm out}'
                              r'=3\times10^{17}$ at $r>5\times10^{15}$')),
    'X_M4_rs6w3': dict(color='tab:brown', lw=1.1, ls='-', zorder=12,
                       label=(r'$M_{\rm sh}=4$, $w_{\rm out}'
                              r'=3\times10^{17}$ at $r>6\times10^{15}$')),
    'X_M4_rs5n04': dict(color='tab:olive', lw=1.1, ls='-', zorder=12,
                        label=(r'$M_{\rm sh}=4$, $n_{\rm out}=0.4$ '
                               r'at $r>5\times10^{15}$')),
    'X_M4_sh3w6': dict(color='tab:cyan', lw=1.1, ls='--', zorder=12,
                       label=(r'$M_{\rm sh}=4$, shell $w=6\times10^{17}$ '
                              r'at $3.2-4\times10^{15}$')),
    'X_M4_sh3w20': dict(color='tab:cyan', lw=1.1, ls=':', zorder=12,
                        label=(r'$M_{\rm sh}=4$, shell $w=2\times10^{18}$ '
                               r'at $3.2-3.6\times10^{15}$')),
    'X_M4_sh3w60': dict(color='teal', lw=1.1, ls='-.', zorder=12,
                        label=(r'$M_{\rm sh}=4$, shell $w=6\times10^{18}$ '
                               r'at $3-3.5\times10^{15}$')),
    'X_M4_rs4w6ej8c7': dict(color='tab:blue', lw=1.1, ls='-.', zorder=12,
                           label=(r'$M_{\rm sh}=4$, $w_{\rm out}=6\times10^{17}$ '
                                  r'at $r>4\times10^{15}$, '
                                  r'$r_{\rm ej}=8\times10^{15}$')),
    'X_M4_rs4w20ej8c7': dict(color='tab:purple', lw=1.1, ls='--', zorder=12,
                            label=(r'$M_{\rm sh}=4$, $w_{\rm out}=2\times10^{18}$ '
                                   r'at $r>4\times10^{15}$, '
                                   r'$r_{\rm ej}=8\times10^{15}$')),
    'X_M4_sh3w6V85': dict(color='darkturquoise', lw=1.1, ls='-', zorder=12,
                          label=(r'$M_{\rm sh}=4$, shell $w=6\times10^{17}$, '
                                 r'$v_{\rm sh}=8500$ km s$^{-1}$')),
    'X_M4_sh3w6r24': dict(color='darkturquoise', lw=1.1, ls='--', zorder=12,
                          label=(r'$M_{\rm sh}=4$, shell $w=6\times10^{17}$ '
                                 r'at $2.4-3.3\times10^{15}$')),
    'X_M4_sh3w6w25': dict(color='darkturquoise', lw=1.1, ls=':', zorder=12,
                          label=(r'$M_{\rm sh}=4$, inner $w=2.5\times10^{17}$, '
                                 r'shell $w=6\times10^{17}$')),
    'X_M4_sh3w6fr': dict(color='darkturquoise', lw=1.1, ls='-.', zorder=12,
                         label=(r'$M_{\rm sh}=4$, fast rise '
                                r'($v_{\rm sh}=8500$, wall $2.4-3.3\times10^{15}$, '
                                r'inner $w=2.5\times10^{17}$)')),
    'X_M4_sh3w6w40': dict(color='cadetblue', lw=1.1, ls='-', zorder=12,
                          label=(r'$M_{\rm sh}=4$, inner $w=4\times10^{17}$, '
                                 r'shell $w=6\times10^{17}$')),
    'X_M4_sh3w6w4n05': dict(color='cadetblue', lw=1.1, ls='--', zorder=12,
                           label=(r'$M_{\rm sh}=4$, inner $w=4\times10^{17}$ '
                                  r'$n=0.5$, shell $w=6\times10^{17}$')),
    'X_M4_sh3w6w40c55': dict(color='cadetblue', lw=1.1, ls='-.', zorder=12,
                             label=(r'$M_{\rm sh}=4$, inner $w=4\times10^{17}$, '
                                    r'wall to $5.5\times10^{15}$')),
    'X_M4_sh3w6nr4': dict(color='tab:pink', lw=1.1, ls='-', zorder=12,
                          label=(r'$M_{\rm sh}=4$, wall $6\times10^{17}$, '
                                 r'shell $n_\rho=4$')),
    'X_M4_sh3w6nr12': dict(color='tab:pink', lw=1.1, ls='--', zorder=12,
                           label=(r'$M_{\rm sh}=4$, wall $6\times10^{17}$, '
                                  r'shell $n_\rho=12$')),
    'X_M4_sh3w6v0hi': dict(color='tab:pink', lw=1.1, ls='-.', zorder=12,
                           label=(r'$M_{\rm sh}=4$, wall $6\times10^{17}$, '
                                  r'$v_0=5000$ km s$^{-1}$')),
    'X_M4_sh3w6M25': dict(color='tab:pink', lw=1.1, ls=':', zorder=12,
                          label=(r'$M_{\rm sh}=2.5$, wall $6\times10^{17}$')),
    'X_M4_sh3w6redo': dict(color='tab:cyan', lw=1.4, ls='-', zorder=13,
                           label=r'$M_{\rm sh}=4$, wall $6\times10^{17}$ (redo)'),
    'X_M4_sh3w6bq10': dict(color='mediumturquoise', lw=1.4, ls='-', zorder=13,
                           label=(r'$M_{\rm sh}=4$, wall $6\times10^{17}$ '
                                  r'(BQ=10)')),
    'X_M4_cut4': dict(color='tab:gray', lw=1.1, ls='-', zorder=12,
                      label=(r'$M_{\rm sh}=4$, sharp CSM edge at '
                             r'$4\times10^{15}$')),
    # v2 grid: X is the bq=10 fiducial; every other entry is a
    # single-parameter perturbation of it.
    'X_M2': dict(color='tab:blue', lw=1.1, ls=':', zorder=12,
                 label=r'$M_{\rm sh}=2\,M_\odot$'),
    'X_nr4': dict(color='tab:pink', lw=1.1, ls='-', zorder=12,
                  label=r'shell $n_\rho=4$'),
    'X_nr12': dict(color='tab:pink', lw=1.1, ls='--', zorder=12,
                   label=r'shell $n_\rho=12$'),
    'X_v0lo': dict(color='tab:olive', lw=1.1, ls='-', zorder=12,
                   label=r'$v_0=1500$ km s$^{-1}$'),
    'X_v0hi': dict(color='tab:olive', lw=1.1, ls='--', zorder=12,
                   label=r'$v_0=5000$ km s$^{-1}$'),
    'X_w1e17': dict(color='tab:purple', lw=1.1, ls='-', zorder=12,
                    label=r'inner $w=1\times10^{17}$'),
    'X_w4e17': dict(color='tab:purple', lw=1.1, ls='--', zorder=12,
                    label=r'inner $w=4\times10^{17}$'),
    'X_wall3e17': dict(color='tab:cyan', lw=1.1, ls='-', zorder=12,
                       label=r'wall $w=3\times10^{17}$'),
    'X_wall2e18': dict(color='tab:cyan', lw=1.1, ls='--', zorder=12,
                       label=r'wall $w=2\times10^{18}$'),
    'X_nowall': dict(color='tab:cyan', lw=1.1, ls=':', zorder=12,
                     label=r'no wall (continuous CSM)'),
    'X_rs24': dict(color='tab:brown', lw=1.1, ls='-', zorder=12,
                   label=r'$r_{\rm split}=2.4\times10^{15}$'),
    'X_cut33': dict(color='teal', lw=1.1, ls='-', zorder=12,
                    label=r'cutoff $3.3\times10^{15}$'),
    'X_cut55': dict(color='teal', lw=1.1, ls='--', zorder=12,
                    label=r'cutoff $5.5\times10^{15}$'),
    'X_V10000_cut70': dict(color='teal', lw=1.1, ls='-.', zorder=12,
                           label=(r'$v_{\rm sh}=10\,000$ km s$^{-1}$, '
                                  r'cutoff $7\times10^{15}$')),
    'X_V10000_cut70_ss': dict(color='tab:olive', lw=1.4, ls='-', zorder=13,
                              label=(r'cut70 + slow shell $w=10^{17}$, '
                                     r'$u=150$ km s$^{-1}$')),
    'X_V10000_cut250w5': dict(color='teal', lw=1.4, ls='-', zorder=13,
                              label=(r'$v_{\rm sh}=10\,000$ km s$^{-1}$, '
                                     r'wall $5\times10^{17}$')),
    'X_V10000_w5ss25': dict(color='tab:olive', lw=1.4, ls='--', zorder=13,
                            label=(r'wall $5\times10^{17}$ + slow wind '
                                   r'$w=2.5\times10^{17}$')),
    'X_V10000_w5ss50': dict(color='tab:olive', lw=1.4, ls='-.', zorder=13,
                            label=(r'slow wind $w=5\times10^{17}$')),
    'X_V10000_w5ss50x': dict(color='darkgoldenrod', lw=1.4, ls='-', zorder=13,
                             label=(r'slow wind $w=5\times10^{17}$ '
                                    r'to $3\times10^{16}$')),
    'X_V10000_w5ss1e18x': dict(color='tab:brown', lw=1.4, ls='-.', zorder=13,
                               label=(r'slow wind $w=10^{18}$ '
                                      r'to $3\times10^{16}$')),
    'X_V10000_wind': dict(color='tab:purple', lw=1.4, ls='-', zorder=13,
                          label=(r'wind $\dot M=0.12$, '
                                 r'SN-A blanket')),
    'X_V10000_wind03': dict(color='tab:purple', lw=1.0, ls=':', zorder=13,
                            label=r'wind $\dot M=0.03$'),
    'X_V10000_wind06': dict(color='tab:purple', lw=1.0, ls='--', zorder=13,
                            label=r'wind $\dot M=0.06$'),
    'X_V10000_wind09': dict(color='tab:purple', lw=1.0, ls='-.', zorder=13,
                            label=r'wind $\dot M=0.09$'),
    'X_V10000_wind035': dict(color='tab:purple', lw=1.0, ls='--', zorder=13,
                             label=r'wind $\dot M=0.035$'),
    'X_V10000_wind04': dict(color='tab:purple', lw=1.0, ls='-.', zorder=13,
                            label=r'wind $\dot M=0.04$'),
    'X_V10000_wind045': dict(color='tab:purple', lw=1.6, ls='-', zorder=14,
                             label=(r'wind $\dot M=0.045$, '
                                    r'SN-A blanket')),
    'X_V10000_wind045b': dict(color='tab:purple', lw=1.6, ls='-.', zorder=14,
                              label=(r'wind $\dot M=0.045$, '
                                     r'blanket $\rho\times5$')),
    'X_V10000_wind045s': dict(color='tab:purple', lw=1.6, ls='--', zorder=14,
                              label=(r'wind $\dot M=0.045$ + shell '
                                     r'$6{\times}10^{17}$ at $1.8\times10^{16}$')),
    'X_V10000_wind02i': dict(color='tab:purple', lw=1.0, ls=':', zorder=13,
                             label=r'wind $\dot M=0.02$, intact SN-A'),
    'X_V10000_wind03i': dict(color='tab:purple', lw=1.0, ls=(0, (1, 1)),
                             zorder=13,
                             label=r'wind $\dot M=0.03$, intact SN-A'),
    'X_V10000_wind045i': dict(color='tab:purple', lw=1.6, ls='-', zorder=14,
                              label=(r'wind $\dot M=0.045$, '
                                     r'intact SN-A')),
    'X_V10000_wind06i': dict(color='tab:purple', lw=1.0, ls='--', zorder=13,
                             label=r'wind $\dot M=0.06$, intact SN-A'),
    'X_V10000_wind09i': dict(color='tab:purple', lw=1.0, ls='-.', zorder=13,
                             label=r'wind $\dot M=0.09$, intact SN-A'),
    'X_V10000_w2m1i': dict(color='blueviolet', lw=1.4, ls='-', zorder=14,
                           label=(r'2-$\dot M$ wind $2{\to}5{\times}10^{17}$ '
                                  r'at $3.5{\times}10^{15}$, intact SN-A')),
    'X_V10000_w2m2i': dict(color='blueviolet', lw=1.4, ls='--', zorder=14,
                           label=(r'2-$\dot M$ wind $2{\to}8{\times}10^{17}$ '
                                  r'at $4.5{\times}10^{15}$, intact SN-A')),
    'X_V10000_w2m3i': dict(color='mediumvioletred', lw=1.4, ls='-',
                           zorder=14,
                           label=(r'2-$\dot M$ wind $2{\to}3.5{\times}10^{17}$'
                                  r' at $3{\times}10^{15}$, intact SN-A')),
    'X_V10000_w2m4i': dict(color='mediumvioletred', lw=1.4, ls='--',
                           zorder=14,
                           label=(r'2-$\dot M$ wind $2{\to}4{\times}10^{17}$'
                                  r' at $3.5{\times}10^{15}$, $n=3$ outer,'
                                  ' intact SN-A')),
    'X_V10000_wind12i': dict(color='tab:pink', lw=1.0, ls='-', zorder=13,
                             label=r'wind $\dot M=0.12$, intact SN-A'),
    'X_V10000_w2m5i': dict(color='crimson', lw=1.4, ls='-', zorder=14,
                           label=(r'2-$\dot M$ wind $2{\to}6{\times}10^{17}$ '
                                  r'at $3{\times}10^{15}$, intact SN-A')),
    'X_V10000_w2m6i': dict(color='crimson', lw=1.4, ls='--', zorder=14,
                           label=(r'2-$\dot M$ wind $2{\to}6{\times}10^{17}$'
                                  r' at $3{\times}10^{15}$, $n=3$ outer,'
                                  ' intact SN-A')),
    'X_V10000_w2m7i': dict(color='tab:orange', lw=1.4, ls='-', zorder=14,
                           label=(r'2-$\dot M$ wind $2{\to}5{\times}10^{17}$ '
                                  r'at $2{\times}10^{15}$, intact SN-A')),
    'X_V10000_w2m8i': dict(color='tab:orange', lw=1.4, ls='--', zorder=14,
                           label=(r'2-$\dot M$ wind $1{\to}5{\times}10^{17}$ '
                                  r'at $3{\times}10^{15}$, intact SN-A')),
    'X_V10000_wind15i': dict(color='tab:pink', lw=1.0, ls='--', zorder=13,
                             label=r'wind $\dot M=0.15$, intact SN-A'),
    'X_V10000_wind20i': dict(color='tab:pink', lw=1.0, ls='-.', zorder=13,
                             label=r'wind $\dot M=0.20$, intact SN-A'),
    'X_V10000_wind13i': dict(color='tab:pink', lw=1.2, ls=':', zorder=13,
                             label=r'wind $\dot M=0.13$, intact SN-A'),
    'X_V10000_wind14i': dict(color='tab:pink', lw=1.2,
                             ls=(0, (1, 1)), zorder=13,
                             label=r'wind $\dot M=0.14$, intact SN-A'),
    'X_V10000_wind12li': dict(color='tab:gray', lw=1.4, ls='-', zorder=14,
                              label=(r'wind $\dot M=0.12$, $v_\infty=100$, '
                                     'intact SN-A')),
    'X_V10000_wind09li': dict(color='tab:gray', lw=1.4, ls='--', zorder=14,
                              label=(r'wind $\dot M=0.09$, $v_\infty=100$, '
                                     'intact SN-A')),
    'X_V10000_w2m9i': dict(color='tab:brown', lw=1.4, ls='-', zorder=14,
                           label=(r'2-$\dot M$ wind $2{\to}1{\times}10^{18}$ '
                                  r'at $3{\times}10^{15}$, intact SN-A')),
    'X_V10000_w2m10i': dict(color='tab:brown', lw=1.4, ls='--', zorder=14,
                            label=(r'2-$\dot M$ wind $2{\to}2{\times}10^{18}$ '
                                   r'at $3{\times}10^{15}$, intact SN-A')),
    'X_V10000_w2m11i': dict(color='tab:brown', lw=1.4, ls='-.', zorder=14,
                            label=(r'2-$\dot M$ wind $2{\to}1.5{\times}10^{18}$'
                                   r' at $4{\times}10^{15}$, intact SN-A')),
    'X_V10000_w2m12i': dict(color='darkred', lw=1.4, ls='-', zorder=14,
                            label=(r'2-$\dot M$ wind $2{\to}6{\times}10^{17}$'
                                   r' at $3{\times}10^{15}$, steep CD,'
                                   ' intact SN-A')),
    'X_V10000_w2m13i': dict(color='darkgreen', lw=1.4, ls='-', zorder=14,
                            label=(r'2-$\dot M$ wind $2{\to}5{\times}10^{17}$'
                                   r' at $4{\times}10^{15}$, steep CD,'
                                   ' intact SN-A')),
    'X_V10000_w2m14i': dict(color='darkgreen', lw=1.4, ls='--', zorder=14,
                            label=(r'2-$\dot M$ wind $2{\to}4{\times}10^{17}$'
                                   r' at $4.5{\times}10^{15}$, steep CD,'
                                   ' intact SN-A')),
    'X_V10000_w2m15i': dict(color='teal', lw=1.6, ls='-', zorder=14,
                            label=(r'eruption $t=-180$ d: '
                                   r'$2{\to}5{\times}10^{17}$ at '
                                   r'$4{\times}10^{15}$, intact SN-A')),
    'X_V8000_w2m13i_sh6': dict(color='darkorange', lw=1.4, ls='-', zorder=14,
                              label=(r'6 $M_\odot$ shell at 8000 km s$^{-1}$,'
                                     ' w2m13 wind, intact SN-A')),
    'X_V14000_w2m13i_sh2': dict(color='darkorange', lw=1.4, ls=':', zorder=14,
                               label=(r'2 $M_\odot$ shell at 14000 km s$^{-1}$,'
                                      ' w2m13 wind, intact SN-A')),
    'X_V10000_w2m13i_sh6': dict(color='darkorange', lw=1.4, ls='--', zorder=14,
                               label=(r'6 $M_\odot$ shell at 10000 km s$^{-1}$,'
                                      ' w2m13 wind, intact SN-A')),
    'X_V12000_w2m13i_sh4': dict(color='darkorange', lw=1.4, ls='-.', zorder=14,
                               label=(r'4 $M_\odot$ shell at 12000 km s$^{-1}$,'
                                      ' w2m13 wind, intact SN-A')),
    'X_V10000_wind2m_inv': dict(color='tab:purple', lw=1.6, ls=':', zorder=14,
                              label=(r'wind $\dot{M}$: $1\to3$'
                                     r'$\times 10^{17}$ at $5{\times}10^{15}$,'
                                     ' SN-A blanket')),
    'X_V10000_wind2m': dict(color='tab:purple', lw=1.6, ls='--', zorder=14,
                            label=(r'wind $\dot M$: $2.5\to1.3\times10^{17}$ '
                                   r'at $10^{16}$, SN-A blanket')),
    'X_V10000_wallwind': dict(color='tab:purple', lw=1.4, ls='--', zorder=13,
                              label=(r'two-$\dot M$ wind '
                                     r'$1.5\to0.5\times10^{18}$, '
                                     r'SN-A blanket')),
    'X_V12000_cut250w5': dict(color='tab:red', lw=1.4, ls='-', zorder=13,
                              label=(r'$v_{\rm sh}=12\,000$ km s$^{-1}$, '
                                     r'wall $5\times10^{17}$')),
    'X_V10000_cut250w5_hhe': dict(color='mediumturquoise', lw=1.4, ls='-',
                                 zorder=13,
                                 label=r'cut250w5 + H/He CSM'),
    'X_V10000_cut250w5_ss2e16': dict(color='tab:cyan', lw=1.4, ls='-',
                                    zorder=13,
                                    label=(r'cut250w5 + thin slow wind '
                                           r'$w=2\times10^{16}$')),
    'X_V10000_cut250w5_ss5e16': dict(color='tab:cyan', lw=1.4, ls='--',
                                    zorder=13,
                                    label=(r'cut250w5 + thin slow wind '
                                           r'$w=5\times10^{16}$')),
    'X_V10000_cut250w5_ni03': dict(color='tab:brown', lw=1.4, ls='-',
                                  zorder=13,
                                  label=r'cut250w5 + Ni $0.3\,M_\odot$'),
    'X_V10000_cut250w5_ni06': dict(color='tab:brown', lw=1.4, ls='--',
                                  zorder=13,
                                  label=r'cut250w5 + Ni $0.6\,M_\odot$'),
    'X_n05': dict(color='tab:green', lw=1.1, ls='-', zorder=12,
                  label=r'CSM $n=0.5$'),
    'X_n10': dict(color='tab:green', lw=1.1, ls='--', zorder=12,
                  label=r'CSM $n=1.0$'),
    'X_t1': dict(color='tab:orange', lw=1.1, ls='-', zorder=12,
                 label=r'$t_{\rm start}=1$ d'),
}


# Legend order: fiducial, then each parameter family in turn.
MODEL_ORDER = ['X', 'X_M2', 'X_M10', 'X_V6000', 'X_V8000',
               'X_V10000', 'X_V12000',
               'X_nr4', 'X_nr12', 'X_v0lo', 'X_v0hi',
               'X_w1e17', 'X_w4e17', 'X_wall3e17', 'X_wall2e18',
               'X_nowall',
               'X_rs24', 'X_cut33', 'X_cut55', 'X_V10000_cut70',
               'X_V10000_cut70_ss', 'X_V10000_cut250w5',
               'X_V10000_w5ss25', 'X_V10000_w5ss50', 'X_V10000_w5ss50x',
               'X_V10000_w5ss1e18x', 'X_V10000_wind', 'X_V10000_wind03',
               'X_V10000_wind06', 'X_V10000_wind09', 'X_V10000_wind035',
               'X_V10000_wind04', 'X_V10000_wind045',
               'X_V10000_wind045b', 'X_V10000_wind045s',
               'X_V10000_wind02i', 'X_V10000_wind03i',
               'X_V10000_wind045i', 'X_V10000_wind06i',
               'X_V10000_wind09i',
               'X_V10000_w2m1i', 'X_V10000_w2m2i', 'X_V10000_w2m3i',
               'X_V10000_w2m4i', 'X_V10000_w2m5i', 'X_V10000_w2m6i',
               'X_V10000_w2m7i', 'X_V10000_w2m8i',
               'X_V10000_w2m9i', 'X_V10000_w2m10i', 'X_V10000_w2m11i',
               'X_V10000_w2m12i', 'X_V10000_w2m13i', 'X_V10000_w2m14i',
               'X_V10000_w2m15i',
               'X_V8000_w2m13i_sh6', 'X_V14000_w2m13i_sh2',
               'X_V10000_w2m13i_sh6', 'X_V12000_w2m13i_sh4',
               'X_V10000_wind12i', 'X_V10000_wind13i', 'X_V10000_wind14i',
               'X_V10000_wind15i', 'X_V10000_wind20i',
               'X_V10000_wind12li', 'X_V10000_wind09li',
               'X_V10000_wind2m', 'X_V10000_wind2m_inv',
               'X_V10000_wallwind',
               'X_n05', 'X_n10', 'X_t1',
               'X_V12000_cut250w5', 'X_V10000_cut250w5_hhe',
               'X_V10000_cut250w5_ss2e16', 'X_V10000_cut250w5_ss5e16',
               'X_V10000_cut250w5_ni03', 'X_V10000_cut250w5_ni06']


def style_for(name):
    return MODEL_STYLE.get(name, dict(color='0.45', lw=0.5, ls='-',
                                      zorder=10, label=name))


# Reference runs outside the live grid (e.g. archived bq=1 best model).
# Each entry maps a model directory to the name used for style/legend.
REFERENCE_DIRS = []


def load_models(grid_dir, min_tmax=0.0):
    models = []
    skipped = []
    for mdir in sorted(d for d in glob.glob(os.path.join(grid_dir, '*'))
                       if os.path.isdir(d)):
        tt = base.find_tt(mdir)
        if tt is None:
            continue
        df = base.load_tt(tt)
        if df is not None and len(df):
            if df['time'].max() < min_tmax:
                skipped.append((os.path.basename(mdir),
                                df['time'].max()))
                continue
            models.append((os.path.basename(mdir), df))
    for rel_dir, name in REFERENCE_DIRS:
        mdir = os.path.normpath(os.path.join(grid_dir, rel_dir))
        tt = base.find_tt(mdir)
        if tt is None:
            continue
        df = base.load_tt(tt)
        if df is not None and len(df):
            if df['time'].max() < min_tmax:
                skipped.append((name, df['time'].max()))
                continue
            models.append((name, df))
    if skipped:
        txt = ', '.join(f'{n} ({t:.0f}d)' for n, t in skipped)
        print(f'  skipped {len(skipped)} models below '
              f'{min_tmax:.0f}d: {txt}')
    return models


def save_both(fig, out_png):
    fig.savefig(out_png, dpi=300, bbox_inches='tight')
    pdf = os.path.splitext(out_png)[0] + '.pdf'
    fig.savefig(pdf, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved: {out_png}')
    print(f'Saved: {pdf}')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--grid-dir', default=GRID_DIR)
    ap.add_argument('--out', default=None)
    ap.add_argument('--fit-end', type=float, default=None,
                    help='last obs phase (d) included in the best-dt fit')
    ap.add_argument('--ebv', type=float, default=None,
                    help='total E(B-V) incl. host; default keeps '
                         'EBV=0.15 (MW foreground only)')
    ap.add_argument('--all-models', action='store_true',
                    help='plot every model; default keeps only the '
                         'intact-SN_A family (names ending in i)')
    ap.add_argument('--min-tmax', type=float, default=200.0,
                    help='drop models whose .tt ends before this '
                         'model time (d); 0 disables')
    args = ap.parse_args()
    if args.ebv is not None:
        base.set_ebv(args.ebv)
        print(f'  E(B-V) total = {base.EBV:.2f} '
              f'(A_g={base.A_G:.2f} A_r={base.A_R:.2f} '
              f'A_o={base.A_O:.2f})')
    out_file = args.out or os.path.join(
        args.grid_dir, 'x_grid_vs_obs.png')

    models = load_models(args.grid_dir, min_tmax=args.min_tmax)
    if not args.all_models:
        # The intact-SN_A family is the current baseline; the older
        # truncated-ejecta models misrepresent SN_A composition.
        models = [(n, d) for n, d in models if n.endswith('i')]
    if not models:
        raise SystemExit(f'no completed models found in {args.grid_dir}')
    print(f'loaded {len(models)} models')

    obs = base.load_obs()

    scored = []
    for name, df in models:
        tpk = df.loc[df['R'].idxmin(), 'time']
        score, dt, cr, ct, cc, pen = base.best_phase_shift(
            df, obs, tpk, args.fit_end)
        scored.append(dict(name=name, df=df, tpk=tpk, dt=dt,
                           score=score, ok=np.isfinite(pen)))
        print(f'  {name:12s} dt={dt:+6.2f}  score={score:8.1f}  '
              f'(early={cr:7.1f} tail={ct:7.1f} col={cc:7.1f} '
              f'pen={pen:7.1f})  '
              f'{"OK" if np.isfinite(pen) else "VIOLATES"}')
    # fiducial first in the legend, then the parameter variations
    order = {name: i for i, name in enumerate(MODEL_ORDER)}
    scored.sort(key=lambda s: order.get(s['name'], len(order)))
    # Highlight the best-fitting model in black so it stands out
    # regardless of which family it belongs to.
    ok = [s for s in scored if s['ok']]
    best_name = (min(ok, key=lambda s: s['score'])['name']
                 if ok else None)

    def style_for_scored(s):
        st = style_for(s['name'])
        if s['name'] == best_name:
            st = dict(st, color='k', lw=2.2, zorder=19,
                      label=st.get('label', s['name']) + ' (best)')
        return st

    BANDS = [
        ('g', 'g', 'Error', 'o', 'tab:green', 'SN 2024pba $g$'),
        ('r', 'R', 'Error', 'o', 'tab:red', 'SN 2024pba $r$'),
        ('i', 'i', 'Error', 's', 'goldenrod', 'SN 2024pba $i$'),
        ('o', 'o', 'Mabs_err', '^', 'darkorange',
         'SN 2024pba $o$ (ATLAS)'),
    ]
    x_lo, x_hi = -60.0, 400.0

    fig, axes = plt.subplots(2, 2, figsize=(2.0 * AA_COL_IN, 4.6),
                           sharex=True, sharey=True)
    band_ax = dict(zip('grio', axes.flat))

    for band, mcol, ecol, fmt, col, lab in BANDS:
        ax = band_ax[band]
        det = obs[band]
        ax.errorbar(det['phase'], det['Mabs'], yerr=det[ecol],
                    fmt=fmt, color=col, ms=3, mec='k', mew=0.5,
                    capsize=1.0, elinewidth=0.6, alpha=0.85, zorder=20)
        ax.set_ylim(top=-19.2, bottom=-12.9)
        ax.text(0.03, 0.06, f'${band}$', transform=ax.transAxes,
                fontsize=8)
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.grid(True, alpha=0.3)
        ax.grid(True, which='minor', alpha=0.12)
    band_ax['g'].set_ylabel('Absolute magnitude')
    band_ax['i'].set_ylabel('Absolute magnitude')
    # grey marker behind the early r-band plateau points, sized from
    # the detections themselves
    axr = band_ax['r']
    plat = obs['r'][(obs['r']['phase'] >= -45) &
                    (obs['r']['phase'] <= -19)]
    if len(plat):
        axr.add_patch(Rectangle(
            (-45.0, plat['Mabs'].min() - 0.05), 26.0,
            plat['Mabs'].max() - plat['Mabs'].min() + 0.10,
            fc='0.85', ec='none', alpha=0.5, zorder=0))

    handles = []
    for s in scored:
        st = style_for_scored(s)
        ph = s['df']['time'] - s['tpk'] + s['dt']
        lab = st['label']
        for band, mcol, _, _, _, _ in BANDS:
            band_ax[band].plot(ph, s['df'][mcol], color=st['color'],
                               lw=st['lw'], ls=st['ls'],
                               zorder=st['zorder'])
        handles.append(Line2D([], [], color=st['color'], lw=st['lw'],
                              ls=st['ls'], label=lab))

    # SN 2009ip comparison (absolute mags, phased to the 2012B peak).
    # Johnson V/R cover both the 2012A precursor and 2012B; the SDSS
    # columns only cover 2012B. Rebin to 2-day bins, dropping NaNs.
    sn09 = pd.read_csv(SN2009IP_FILE)
    win = (sn09['mjd'] > 56140.0) & (sn09['mjd'] < 56460.0)
    for band, mcol in (('g', 'V'), ('r', 'R')):
        sub = sn09.loc[win, ['mjd', mcol]].dropna()
        sub = sub[np.isfinite(sub['mjd']) & np.isfinite(sub[mcol])]
        sub = sub.copy()
        sub['bin'] = np.floor(sub['mjd'] / SN2009IP_BIN_D)
        day = (sub.groupby('bin')
                  .agg(mjd=('mjd', 'mean'), mag=(mcol, 'median'))
                  .dropna()
                  .sort_values('mjd'))
        ph = day['mjd'].to_numpy() - SN2009IP_PEAK_MJD
        mag = day['mag'].to_numpy()
        # keep the connecting line off long observing gaps
        seg = np.split(np.arange(len(day)),
                       np.where(np.diff(ph) > 10.0)[0] + 1)
        for s in seg:
            band_ax[band].plot(ph[s], mag[s], '-', color='0.45',
                               lw=0.7, marker='o', ms=1.5, mec='none',
                               alpha=0.8, zorder=8)
    obs_handles = [
        Line2D([], [], marker='o', color='tab:green', ls='', ms=3,
               mec='k', mew=0.5, label='$g$'),
        Line2D([], [], marker='o', color='tab:red', ls='', ms=3,
               mec='k', mew=0.5, label='$r$'),
        Line2D([], [], marker='s', color='goldenrod', ls='', ms=3,
               mec='k', mew=0.5, label='$i$'),
        Line2D([], [], marker='^', color='darkorange', ls='', ms=3,
               mec='k', mew=0.5, label='$o$'),
        Line2D([], [], color='0.45', lw=0.7, marker='o', ms=2,
               mec='none', label='SN 2009ip'),
    ]
    # filters key inside the upper-right panel; model key above the
    # figure
    band_ax['r'].legend(handles=obs_handles, fontsize=6, ncol=1,
                        frameon=False, loc='upper right',
                        handletextpad=0.4, borderaxespad=0.3,
                        labelspacing=0.3)
    fig.legend(handles=handles, fontsize=5.4, ncol=4, frameon=False,
               loc='upper center', bbox_to_anchor=(0.5, 1.0),
               columnspacing=0.6, handlelength=1.6,
               handletextpad=0.4)

    for ax in axes[-1]:
        ax.set_xlabel('Days until SN 2024pba peak')
    axes[-1, 0].set_xlim(x_lo, x_hi)
    fig.align_ylabels(axes)
    # Reserve enough top margin for the multi-row model legend so it
    # does not overlap the top panel.
    n_leg_rows = max(1, -(-len(handles) // 4))
    rect_top = min(0.92, 1.0 - 0.030 * n_leg_rows)
    plt.tight_layout(rect=[0, 0, 1, rect_top])
    save_both(fig, out_file)

    # --- colour companion figure: g-r, o-g, o-r, r-i ---
    colour_defs = [('g', 'r', 0.75), ('o', 'g', 1.5), ('o', 'r', 1.5),
                   ('r', 'i', 0.75)]
    figc, axesc = plt.subplots(4, 1, figsize=(AA_COL_IN, 6.4),
                               sharex=True)
    for ax, (b1, b2, sep) in zip(axesc, colour_defs):
        pairs = base.obs_color_pairs(obs, b1, b2, max_sep=sep)
        if len(pairs):
            ax.errorbar(pairs['phase'], pairs['col'], yerr=pairs['err'],
                        fmt='o', color='k', ms=3, mec='k', mew=0.5,
                        capsize=1.0, elinewidth=0.6, alpha=0.85,
                        zorder=20,
                        label=f'SN 2024pba ${b1} - {b2}$')
        for s in scored:
            st = style_for_scored(s)
            ph = s['df']['time'] - s['tpk'] + s['dt']
            ax.plot(ph,
                    s['df'][base.OBS2MOD[b1]] - s['df'][base.OBS2MOD[b2]],
                    color=st['color'], lw=st['lw'], ls=st['ls'],
                    zorder=st['zorder'])
        ax.set_ylabel(f'${b1} - {b2}$ (mag)')
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.grid(True, alpha=0.3)
        ax.grid(True, which='minor', alpha=0.12)
        if len(pairs):
            c_lo = min(pairs['col'].min(), -0.4) - 0.3
            c_hi = max(pairs['col'].max(), 0.4) + 0.3
            ax.set_ylim(c_hi, c_lo)  # blue at the top
            ax.legend(fontsize=6, frameon=False, loc='lower right',
                      handletextpad=0.4, borderaxespad=0.3)
    figc.legend(handles=handles, fontsize=5.4, ncol=3, frameon=False,
                loc='upper center', bbox_to_anchor=(0.5, 1.0),
                columnspacing=0.6, handlelength=1.6,
                handletextpad=0.4)
    axesc[-1].set_xlim(x_lo, x_hi)
    axesc[-1].set_xlabel('Days until SN 2024pba peak')
    figc.align_ylabels(axesc)
    n_leg_rows_c = max(1, -(-len(handles) // 3))
    rect_top_c = min(0.90, 1.0 - 0.030 * n_leg_rows_c)
    plt.tight_layout(rect=[0, 0, 1, rect_top_c])
    save_both(figc, out_file.replace('.png', '_gr.png'))


if __name__ == '__main__':
    main()
