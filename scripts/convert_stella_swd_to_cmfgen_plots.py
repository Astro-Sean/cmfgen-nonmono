"""Diagnostic plotting functions for the STELLA .swd -> CMFGEN SN_HYDRO_DATA converter.

Extracted from convert_stella_swd_to_cmfgen.py to separate visualization from
conversion logic.  All functions accept numpy arrays in outer->inner (decreasing
radius) order, matching the converter's internal convention.

Imported by the main converter via:
    from convert_stella_swd_to_cmfgen_plots import (
        plot_comparison, plot_mass_fractions_stacked, plot_zone_density,
        generate_grid_diagnostic_report, plot_dtau_resolution_panels,
        plot_sigma_vs_index, plot_structure_overview,
    )
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# Physical constants (duplicated from main converter to keep module standalone)
# ---------------------------------------------------------------------------
AVOGADRO = 6.02214076e23

# Atomic mass and charge dictionaries used for density computations in plots.
# Includes both internal STELLA names and CMFGEN output aliases (e.g. HYD, HE, ...)
_ELEM_A = {
    'H': 1.0, 'HYD': 1.0, 'He': 4.0, 'HE': 4.0,
    'C': 12.0, 'CARB': 12.0,
    'N': 14.0, 'NIT': 14.0,
    'O': 16.0, 'OXY': 16.0,
    'Ne': 20.0, 'NEON': 20.0,
    'Na': 23.0, 'SOD': 23.0,
    'Mg': 24.0, 'MAG': 24.0,
    'Al': 27.0, 'ALUM': 27.0,
    'Si': 28.0, 'SIL': 28.0,
    'S': 32.0, 'SUL': 32.0,
    'Ar': 40.0, 'ARG': 40.0,
    'Ca': 40.0, 'CAL': 40.0,
    'FePeak': 56.0, 'Fe': 56.0, 'Fe56': 56.0, 'IRON': 56.0,
    'Ni58': 58.0, 'Ni56': 56.0, 'NICK': 56.0, 'COB': 56.0,
}

_ELEM_Z = {
    'H': 1.0, 'HYD': 1.0, 'He': 2.0, 'HE': 2.0,
    'C': 6.0, 'CARB': 6.0,
    'N': 7.0, 'NIT': 7.0,
    'O': 8.0, 'OXY': 8.0,
    'Ne': 10.0, 'NEON': 10.0,
    'Na': 11.0, 'SOD': 11.0,
    'Mg': 12.0, 'MAG': 12.0,
    'Al': 13.0, 'ALUM': 13.0,
    'Si': 14.0, 'SIL': 14.0,
    'S': 16.0, 'SUL': 16.0,
    'Ar': 18.0, 'ARG': 18.0,
    'Ca': 20.0, 'CAL': 20.0,
    'FePeak': 26.0, 'Fe': 26.0, 'Fe56': 26.0, 'IRON': 26.0,
    'Ni58': 28.0, 'Ni56': 28.0, 'NICK': 28.0, 'COB': 27.0,
}

_MF_COLORS = {
    'H': 'tab:blue', 'He': 'tab:cyan', 'C': 'tab:green', 'N': 'olive',
    'O': 'tab:orange', 'Ne': 'tab:pink', 'Na': 'wheat', 'Mg': 'lime',
    'Al': 'silver', 'Si': 'tab:purple', 'S': 'gold', 'Ar': 'tan',
    'Ca': 'tab:red', 'FePeak': 'tab:brown', 'IRON': 'tab:brown',
    'Fe': 'tab:brown', 'NICK': 'darkgray', 'Ni56': 'darkgray',
    'Ni58': 'dimgray', 'COB': 'teal',
}

_PLOT_ORDER = [
    'H', 'He', 'O', 'C', 'Si', 'S', 'Ca', 'Fe', 'IRON', 'FePeak',
    'NICK', 'Ni56', 'Ni58', 'COB', 'Mg', 'Na', 'Ne', 'N', 'Al', 'Ar',
]


def _compute_atom_density(rho_arr, mf_dict):
    """Compute atom number density n_atom = rho * N_A * sum(X_i / A_i)."""
    inv_A = np.zeros_like(rho_arr)
    if mf_dict is not None:
        for elem, frac in mf_dict.items():
            A = _ELEM_A.get(elem, 28.0)
            inv_A += frac / A
    inv_A = np.where(inv_A > 1e-30, inv_A, 1.0 / 28.0)
    return rho_arr * AVOGADRO * inv_A


def _compute_ne_max(rho_arr, mf_dict):
    """Compute full-ionization electron density n_e,max = rho * N_A * sum(X_i * Z_i / A_i)."""
    Z_over_A = np.zeros_like(rho_arr)
    if mf_dict is not None:
        for elem, frac in mf_dict.items():
            A = _ELEM_A.get(elem, 28.0)
            Z = _ELEM_Z.get(elem, 14.0)
            Z_over_A += frac * Z / A
    Z_over_A = np.where(Z_over_A > 1e-30, Z_over_A, 1.0 / 28.0)
    return rho_arr * AVOGADRO * Z_over_A


def _save_fig(fig, output_file, suffix, dpi=150):
    """Save figure to the same directory as output_file with given suffix."""
    plot_file = str(Path(output_file).parent / f'SN_HYDRO_DATA_{suffix}.png')
    os.makedirs(os.path.dirname(plot_file) or '.', exist_ok=True)
    fig.savefig(plot_file, dpi=dpi, bbox_inches='tight')
    print(f"Plot saved to {plot_file}")
    plt.close(fig)


def _add_depth_axis(ax, r_desc, xlog=True):
    """Secondary top axis with the CMFGEN depth index (1 = outer boundary).

    ``r_desc`` is the radius grid in outer->inner (decreasing) order.
    Tick depths are dropped where adjacent labels would overprint --
    the grid compresses near shocks and composition interfaces, so
    evenly spaced depth ticks pile up there in a log-radius axis.
    """
    r_desc = np.asarray(r_desc, dtype=float)
    nd = len(r_desc)
    if nd < 2:
        return
    lo, hi = ax.get_xlim()
    candidates = np.arange(1, nd + 1, 10, dtype=int)
    if candidates[-1] != nd:
        candidates = np.append(candidates, nd)
    pos = r_desc[np.clip(candidates - 1, 0, nd - 1)]
    # Drop ticks that fall outside the drawn radius range.
    in_range = (pos >= min(lo, hi)) & (pos <= max(lo, hi))
    candidates, pos = candidates[in_range], pos[in_range]
    if len(candidates) == 0:
        return
    if xlog:
        axis_pos = np.log(pos)
        span = abs(np.log(max(lo, hi)) - np.log(min(lo, hi)))
    else:
        axis_pos = pos.astype(float)
        span = abs(hi - lo)
    keep = [0]
    for i in range(1, len(candidates)):
        if abs(axis_pos[i] - axis_pos[keep[-1]]) > 0.07 * span:
            keep.append(i)
    sec = ax.secondary_xaxis("top")
    sec.set_xticks(pos[keep])
    sec.set_xticklabels([str(t) for t in candidates[keep]])
    sec.set_xlabel("CMFGEN depth index (1 = outer)", fontsize=9)


# ---------------------------------------------------------------------------
# Plot functions
# ---------------------------------------------------------------------------

def plot_comparison(orig_r, orig_rho, orig_T, orig_v, orig_tau, orig_sigma,
                   new_r, new_rho, new_T, new_v, new_tau, new_sigma,
                   output_file, time_days, tau_min=1e-5, tau_max=30.0,
                   orig_mass_fractions=None, new_mass_fractions=None,
                   orig_ne=None, new_ne=None,
                   ne_from_stella=False, orig_tau_nlte=None, orig_tau_const=None, kappa=0.34):
    """Plot before/after comparison of STELLA to CMFGEN conversion.

    Uses drawstyle='steps-mid' for all profile plots to reflect the discrete
    zone-centered nature of the data.
    """
    orig_r_rsun = orig_r  # cm
    new_r_rsun = new_r  # cm

    # Determine consistent x-axis limits (radius)
    r_min = min(np.min(orig_r_rsun), np.min(new_r_rsun)) * 0.9
    r_max = max(np.max(orig_r_rsun), np.max(new_r_rsun)) * 1.1

    # 3x2 layout: density | temperature, velocity | sigma, tau | densities.
    # Mass fractions are covered by the dedicated SN_HYDRO_DATA_mass_fraction
    # plot and are no longer repeated here.
    fig, axes = plt.subplots(3, 2, figsize=(16, 18))
    fig.suptitle(f'STELLA to CMFGEN Conversion: {time_days:.1f}-day snapshot', fontsize=16)

    # ---- Row 1: Density | Temperature ----
    # Plot 1: Density
    ax = axes[0, 0]
    ax.step(orig_r_rsun, orig_rho, 'bs', alpha=0.7, markersize=3, where='mid', label=f'STELLA ({len(orig_r)} points)')
    ax.step(new_r_rsun, new_rho, 'ro', alpha=0.7, markersize=2, where='mid', label=f'CMFGEN ({len(new_r)} points)')
    ax.set_ylabel('Density (g/cm³)', fontsize=12)
    ax.set_title('Density Profile', fontsize=14)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_xlim(r_min, r_max)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Radius (cm)', fontsize=10)

    # Plot 2: Temperature
    ax = axes[0, 1]
    ax.step(orig_r_rsun, orig_T, 'bs', alpha=0.7, markersize=3, where='mid', label=f'STELLA ({len(orig_r)} points)')
    ax.step(new_r_rsun, new_T, 'ro', alpha=0.7, markersize=2, where='mid', label=f'CMFGEN ({len(new_r)} points)')
    ax.set_ylabel('Temperature (K)', fontsize=12)
    ax.set_title('Temperature Profile', fontsize=14)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_xlim(r_min, r_max)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Radius (cm)', fontsize=10)

    # ---- Row 2: Velocity | Sigma ----
    # Plot 3: Velocity
    ax = axes[1, 0]
    ax.step(orig_r_rsun, orig_v, 'bs', alpha=0.7, markersize=4, where='mid', label=f'STELLA ({len(orig_r)} points)')
    ax.step(new_r_rsun, new_v, 'ro', alpha=0.8, markersize=2, where='mid', label=f'CMFGEN ({len(new_r)} points)')
    ax.set_ylabel('Velocity (km/s)', fontsize=12)
    ax.set_title('Velocity Profile', fontsize=14)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_xlim(r_min, r_max)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Radius (cm)', fontsize=10)

    # Plot 4: Sigma
    ax = axes[1, 1]
    ax.step(orig_r_rsun, orig_sigma, 'bs', alpha=0.7, markersize=3, where='mid', label=f'STELLA ({len(orig_r)} points)')
    ax.step(new_r_rsun, new_sigma, 'ro', alpha=0.7, markersize=2, where='mid', label=f'CMFGEN ({len(new_r)} points)')
    ax.axhline(y=-0.99, color='k', linestyle=':', alpha=0.5, label='σ=-0.99')
    ax.axhline(y=-1.0, color='gray', linestyle='--', alpha=0.5, label='σ=-1.0')
    ax.set_xlabel('Radius (cm)', fontsize=10)
    ax.set_ylabel('Sigma (σ)', fontsize=12)
    ax.set_title('Sigma Profile', fontsize=14)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_yscale('symlog', linthresh=1.0)
    ax.set_xlim(r_min, r_max)
    ax.set_xscale('log')

    # ---- Row 3: Optical depth | Electron density ----
    # Plot 5: Optical depth
    ax = axes[2, 0]
    orig_tau_plot = np.clip(orig_tau, 1e-10, None)
    new_tau_plot = np.clip(new_tau, 1e-10, None)
    ax.step(orig_r_rsun, orig_tau_plot, 'bs', alpha=0.5, markersize=3, where='mid', label=f'STELLA analytical ({len(orig_r)} pts)')
    if orig_tau_nlte is not None:
        orig_tau_nlte_plot = np.clip(orig_tau_nlte, 1e-10, None)
        ax.step(orig_r_rsun, orig_tau_nlte_plot, 'g^', alpha=0.5, markersize=2, where='mid', label=f'STELLA NLTE ({len(orig_r)} pts)')
    if orig_tau_const is not None:
        orig_tau_const_plot = np.clip(orig_tau_const, 1e-10, None)
        ax.step(orig_r_rsun, orig_tau_const_plot, 'b--', alpha=0.5, markersize=2, where='mid', label=f'STELLA κ={kappa} ({len(orig_r)} pts)')
    ax.step(new_r_rsun, new_tau_plot, 'ro', alpha=0.7, markersize=2, where='mid', label=f'CMFGEN κ={kappa} ({len(new_r)} points)')
    ax.axhline(y=2/3, color='orange', linestyle=':', alpha=0.8, linewidth=2, label='τ = 2/3')
    ax.axhline(y=tau_min, color='gray', linestyle='--', alpha=0.5, linewidth=1, label=f'τ_min={tau_min:.0e}')
    ax.axhline(y=tau_max, color='gray', linestyle='--', alpha=0.5, linewidth=1, label=f'τ_max={tau_max:.0f}')
    ax.set_xlabel('Radius (cm)', fontsize=10)
    ax.set_ylabel('Optical Depth τ', fontsize=12)
    ax.set_title('Optical Depth Profile', fontsize=14)
    ax.legend(fontsize=7, ncol=2)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_yscale('log')
    ax.set_xlim(r_min, r_max)
    ax.set_xscale('log')
    all_tau = [new_tau, orig_tau]
    if orig_tau_nlte is not None:
        all_tau.append(orig_tau_nlte)
    if orig_tau_const is not None:
        all_tau.append(orig_tau_const)
    tau_lo = max(min(t.min() for t in all_tau), 1e-10)
    tau_hi = max(t.max() for t in all_tau)
    ax.set_ylim(tau_lo * 0.3, tau_hi * 2.0)

    # Plot 6: Electron density + Atom density
    ax = axes[2, 1]
    has_ne = orig_ne is not None and new_ne is not None
    has_mf = new_mass_fractions is not None

    if has_ne or has_mf:
        ne_source = 'STELLA' if ne_from_stella else 'Saha'
        if has_ne:
            ax.step(orig_r_rsun, orig_ne, 'bs', alpha=0.7, markersize=3, where='mid', label=f'{ne_source} n_e ({len(orig_r)})')
            ax.step(new_r_rsun, new_ne, 'ro', alpha=0.7, markersize=2, where='mid', label=f'{ne_source} n_e interp ({len(new_r)})')
        if has_mf:
            new_natom = _compute_atom_density(new_rho, new_mass_fractions)
            ax.step(new_r_rsun, new_natom, 'g^', alpha=0.7, markersize=2, where='mid', label=f'n_atom CMFGEN ({len(new_r)})')
            if orig_mass_fractions is not None:
                orig_natom = _compute_atom_density(orig_rho, orig_mass_fractions)
                ax.step(orig_r_rsun, orig_natom, 'gD', alpha=0.5, markersize=3, where='mid', label=f'n_atom STELLA ({len(orig_r)})')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('Radius (cm)', fontsize=10)
        ax.set_ylabel('Number Density (cm⁻³)', fontsize=12)
        ax.set_title('Electron & Atom Density', fontsize=14)
        ax.legend(fontsize=8, ncol=2)
        ax.grid(True, alpha=0.3, which='both')
        ax.set_xlim(r_min, r_max)
    else:
        ax.text(0.5, 0.5, 'Density data not available', ha='center', va='center', transform=ax.transAxes)
        ax.set_title('Electron & Atom Density', fontsize=14)
        ax.set_xlabel('Radius (cm)', fontsize=10)
        ax.set_xlim(r_min, r_max)
        ax.set_xscale('log')

    # Depth-index top axis on every panel (all share the same radius range).
    for ax in axes.flat:
        _add_depth_axis(ax, new_r)

    plt.tight_layout()
    _save_fig(fig, output_file, 'comparison')


def plot_mass_fractions_stacked(radius_cm, density, mass_fractions, time_days, output_file,
                                trimmed_radius_cm=None, trimmed_mass_fractions=None,
                                tau=None, tau_min=1e-5, tau_max=30.0):
    """Create a mass fraction plot vs radius.
    
    Plots original STELLA data as dashed lines and trimmed CMFGEN data as solid lines.
    Isotopes (NICK, COB, IRON) are shown with distinct dotted/dash-dot styles.
    """
    # Original data (STELLA - full model)
    r_orig = np.asarray(radius_cm)
    
    # Trimmed data (CMFGEN - selected tau range)
    if trimmed_radius_cm is not None and trimmed_mass_fractions is not None:
        r_trim = np.asarray(trimmed_radius_cm)
        mf_trim = trimmed_mass_fractions
        has_trimmed = True
    else:
        has_trimmed = False
    
    # Elemental species: (key_in_orig, key_in_trim, display_label, color)
    ELEM_SPECS = [
        ('H',      'H',    'H',    '#FF6B6B'),
        ('He',     'He',   'He',   '#4ECDC4'),
        ('C',      'C',    'C',    '#45B7D1'),
        ('O',      'O',    'O',    '#96CEB4'),
        ('Si',     'Si',   'Si',   '#9B59B6'),
        ('Ca',     'Ca',   'Ca',   '#E74C3C'),
        ('FePeak', 'IRON', 'Fe',   '#7F8C8D'),
    ]
    # Isotope species shown only in trimmed (CMFGEN) data (radioactive parents only)
    ISO_SPECS = [
        ('NICK', 'Ni56 (isotope)', '#e377c2', '-.'),
        ('COB',  'Co56 (isotope)', '#ff7f0e', '-.'),
    ]
    
    fig, ax = plt.subplots(figsize=(12, 6))

    x_orig = r_orig  # cm

    # Shade the CMFGEN selection region so the user can see which part of the
    # STELLA profile was selected (solid) vs excluded (dashed only).
    if has_trimmed:
        x_trim = r_trim  # cm
        r_sel_min, r_sel_max = x_trim.min(), x_trim.max()
        ax.axvspan(r_sel_min, r_sel_max, color='#E8E8E8', alpha=0.5, zorder=0)
        ax.axvline(r_sel_min, color='gray', linewidth=0.8, linestyle=':', alpha=0.7, zorder=1)
        ax.axvline(r_sel_max, color='gray', linewidth=0.8, linestyle=':', alpha=0.7, zorder=1)

    # Plot original STELLA data as dashed lines (FULL profile, all zones)
    for key_orig, key_trim, label, color in ELEM_SPECS:
        if key_orig in mass_fractions and np.any(mass_fractions[key_orig] > 1e-30):
            ax.step(x_orig, mass_fractions[key_orig], color=color, where='mid',
                    linewidth=1.8, linestyle='--', alpha=0.85, zorder=5)

    # Plot trimmed CMFGEN elemental data as solid lines ON TOP of the dashed
    # (only within the selection region)
    if has_trimmed:
        for key_orig, key_trim, label, color in ELEM_SPECS:
            if key_trim in mf_trim and np.any(mf_trim[key_trim] > 1e-30):
                ax.step(x_trim, mf_trim[key_trim], color=color, where='mid',
                        linewidth=2.5, linestyle='-', alpha=0.9, label=label, zorder=6)
        # Plot isotopes as dash-dot
        for key, label, color, ls in ISO_SPECS:
            if key in mf_trim and np.any(mf_trim[key] > 1e-30):
                ax.step(x_trim, mf_trim[key], color=color, where='mid',
                        linewidth=2, linestyle=ls, alpha=0.8, label=label, zorder=6)

    # Build legend with explicit STELLA / CMFGEN proxy entries
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    stella_proxy = Line2D([0], [0], color='gray', linewidth=1.8, linestyle='--', alpha=0.85, label='STELLA (dashed)')
    cmfgen_proxy = Line2D([0], [0], color='gray', linewidth=2.5, linestyle='-', alpha=0.9, label='CMFGEN (solid)')
    sel_proxy = Patch(facecolor='#E8E8E8', edgecolor='gray', linewidth=0.8, linestyle=':', alpha=0.5, label='CMFGEN selection')
    handles, labels = ax.get_legend_handles_labels()
    # Deduplicate element labels (solid lines already labelled)
    seen = set()
    deduped_handles, deduped_labels = [], []
    for h, l in zip(handles, labels):
        if l not in seen:
            seen.add(l)
            deduped_handles.append(h)
            deduped_labels.append(l)
    handles = deduped_handles + [stella_proxy, cmfgen_proxy, sel_proxy]
    labels = deduped_labels + ['STELLA (dashed)', 'CMFGEN (solid)', 'CMFGEN selection']

    ax.set_xlabel(r'Radius (cm)', fontsize=12)
    ax.set_ylabel('Mass Fraction', fontsize=12)
    ax.set_title(f'Mass Fractions: STELLA (dashed, full profile) vs CMFGEN (solid, selection)\nt = {time_days:.1f} days', fontsize=13)
    ax.set_xlim(x_orig.min() * 0.9, x_orig.max() * 1.1)
    ax.set_ylim(0, 1.05)
    ax.set_xscale('log')
    ax.legend(handles=handles, labels=labels, loc='center left', bbox_to_anchor=(1.01, 0.5),
              fontsize=9, framealpha=0.9)
    ax.grid(True, alpha=0.3, which='both')

    # CMFGEN depth index on the trimmed grid (1 = outer boundary).
    if has_trimmed:
        _add_depth_axis(ax, x_trim)

    # Reserve space for the legend on the right
    fig.subplots_adjust(right=0.72, left=0.07, top=0.88, bottom=0.12)
    _save_fig(fig, output_file, 'mass_fraction')


def plot_zone_density(r_cm, rho, kappa, tau, output_file, time_days,
                      tau_breakpoints=None, thin_tau_threshold=1.0):
    """Plot zone density (zones per equal-log-radius bin) vs radius.

    Shows how grid points are distributed in radius, with tau overlaid
    to highlight the tau-based zoning strategy.  Vertical lines mark
    tau breakpoints and the thin-region threshold.
    """
    n = len(r_cm)
    r_rsun = r_cm  # cm

    # Compute zone density via histogram in log(radius)
    log_r = np.log10(r_cm)
    n_bins = 60
    counts, bin_edges = np.histogram(log_r, bins=n_bins)
    bin_centers = 10 ** (0.5 * (bin_edges[:-1] + bin_edges[1:]))
    bin_widths_log = bin_edges[1:] - bin_edges[:-1]
    # zones per dex in log(r)
    density = counts / bin_widths_log

    fig, ax1 = plt.subplots(figsize=(14, 7))

    # Zone density bar plot
    ax1.bar(bin_centers, density, width=np.diff(10**bin_edges),
            color='steelblue', alpha=0.7, label='Zones per dex in log(r)')
    ax1.set_xlabel('Radius (10$^{10}$ cm)', fontsize=12)
    ax1.set_ylabel('Zone density (zones / dex)', fontsize=12, color='steelblue')
    ax1.tick_params(axis='y', labelcolor='steelblue')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.grid(True, alpha=0.3, which='both')

    # Tau overlay on right y-axis
    ax2 = ax1.twinx()
    ax2.step(r_cm, tau, 'r-', linewidth=2, where='mid', alpha=0.8, label='Cumulative τ')
    ax2.set_ylabel('Cumulative optical depth τ', fontsize=12, color='r')
    ax2.tick_params(axis='y', labelcolor='r')
    ax2.set_yscale('log')

    # Mark tau breakpoints
    if tau_breakpoints is not None:
        for bp in tau_breakpoints:
            # Find radius where tau = bp
            idx = np.searchsorted(tau, bp)
            if 0 < idx < n:
                r_bp = r_cm[idx]
                ax1.axvline(r_bp, color='green', linestyle='--', alpha=0.7, linewidth=1.5)
                ax1.annotate(f'τ={bp}', xy=(r_bp, ax1.get_ylim()[1]),
                             fontsize=9, color='green', ha='center', va='top')

    # Mark thin-region threshold
    idx_thin = np.searchsorted(tau, thin_tau_threshold)
    if 0 < idx_thin < n:
        r_thin = r_cm[idx_thin]
        ax1.axvline(r_thin, color='orange', linestyle=':', alpha=0.8, linewidth=2)
        ax1.annotate(f'τ={thin_tau_threshold} (thin threshold)',
                     xy=(r_thin, ax1.get_ylim()[1] * 0.5),
                     fontsize=9, color='orange', ha='center', va='top')

    # Title and legend
    ax1.set_title(f'Zone Density vs Radius ({n} zones)\nt = {time_days:.1f} days', fontsize=14)

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right', fontsize=10)

    _add_depth_axis(ax1, r_cm)

    plt.tight_layout()
    _save_fig(fig, output_file, 'zone_density')


def generate_grid_diagnostic_report(
    r_cm: np.ndarray,
    v_kms: np.ndarray,
    rho: np.ndarray,
    T_K: np.ndarray,
    kappa: np.ndarray,
    tau: np.ndarray,
    sigma: np.ndarray,
    mass_fractions: dict,
    time_days: float,
    output_file: str,
    n_e: np.ndarray | None = None,
) -> dict:
    """Generate a machine-readable grid diagnostic report.

    Produces a JSON file with all quantities specified in the task:
    ND, bounds, max deltas, worst-resolved cells, composition interfaces, etc.

    Parameters
    ----------
    r_cm : array, outer→inner (decreasing)
    v_kms, rho, T_K, kappa, tau, sigma : arrays, same orientation
    mass_fractions : dict of element → array
    time_days : epoch time
    output_file : path to SN_HYDRO_DATA (report written alongside as .json)

    Returns
    -------
    dict : the report data (also written as JSON)
    """
    n = len(r_cm)
    eps = 1e-30

    # Per-cell quantities
    dtau = np.diff(tau)
    dtau_abs = np.abs(dtau)
    r_mid = 0.5 * (r_cm[:-1] + r_cm[1:])

    # dlog(tau) — guard against tau=0
    tau_safe = np.maximum(tau, 1e-10)
    dlogtau = np.abs(np.diff(np.log10(tau_safe)))

    # dtau/dr and dtau/dln(r)
    dr = np.abs(np.diff(r_cm))
    dr_safe = np.maximum(dr, eps)
    dtau_dr = dtau_abs / dr_safe
    dtau_dlnr = dtau_dr * r_mid

    # Physical gradients
    dlogrho = np.abs(np.diff(np.log10(np.maximum(rho, eps))))
    dlogT = np.abs(np.diff(np.log10(np.maximum(T_K, eps))))
    dv = np.abs(np.diff(v_kms))

    # Electron density gradient
    dlogne = None
    if n_e is not None:
        dlogne = np.abs(np.diff(np.log10(np.maximum(n_e, eps))))

    # Composition jumps (max over all elements)
    max_comp_jump = 0.0
    worst_comp_elem = ""
    worst_comp_idx = -1
    for elem, mf in mass_fractions.items():
        if np.all(mf < 1e-10):
            continue
        jumps = np.abs(np.diff(mf))
        # Relative jump: |dX| / max(X, floor)
        mf_safe = np.maximum(mf[:-1], 1e-6)
        rel_jumps = jumps / mf_safe
        idx = int(np.argmax(rel_jumps))
        if rel_jumps[idx] > max_comp_jump:
            max_comp_jump = float(rel_jumps[idx])
            worst_comp_elem = elem
            worst_comp_idx = idx

    # Fraction of tau in largest cell
    total_tau = float(np.sum(dtau_abs))
    frac_largest_cell = float(np.max(dtau_abs) / total_tau) if total_tau > 0 else 0.0

    # Composition interfaces (count)
    n_interfaces = 0
    interface_indices = []
    for elem, mf in mass_fractions.items():
        if np.all(mf < 1e-10):
            continue
        jumps = np.abs(np.diff(mf))
        mf_safe = np.maximum(mf[:-1], 1e-6)
        rel_jumps = jumps / mf_safe
        for i in np.where(rel_jumps > 0.3)[0]:
            if i not in interface_indices:
                interface_indices.append(int(i))
                n_interfaces += 1

    # Velocity turning points
    n_turns = 0
    for i in range(1, n - 1):
        if (v_kms[i] > v_kms[i - 1] and v_kms[i] > v_kms[i + 1]) or \
           (v_kms[i] < v_kms[i - 1] and v_kms[i] < v_kms[i + 1]):
            n_turns += 1

    # Negative sigma regions
    neg_sigma = sigma < 0
    n_neg_sigma = int(np.sum(neg_sigma))

    # Worst cells
    worst_dtau_idx = int(np.argmax(dtau_abs))
    worst_dlogtau_idx = int(np.argmax(dlogtau))
    worst_drho_idx = int(np.argmax(dlogrho))
    worst_dT_idx = int(np.argmax(dlogT))
    worst_dv_idx = int(np.argmax(dv))

    report = {
        "time_days": float(time_days),
        "ND": int(n),
        "r_inner_cm": float(r_cm[-1]),
        "r_outer_cm": float(r_cm[0]),
        "v_inner_kms": float(v_kms[-1]),
        "v_outer_kms": float(v_kms[0]),
        "tau_inner": float(tau[-1]),
        "tau_outer": float(tau[0]),
        "max_delta_log_tau": float(np.max(dlogtau)),
        "median_delta_log_tau": float(np.median(dlogtau)),
        "max_delta_rho": float(np.max(dlogrho)),
        "max_delta_T": float(np.max(dlogT)),
        "max_delta_v": float(np.max(dv)),
        "max_delta_ne": float(np.max(dlogne)) if dlogne is not None else None,
        "max_composition_jump": float(max_comp_jump),
        "maximum_cell_delta_tau": float(np.max(dtau_abs)),
        "fraction_of_tau_in_largest_cell": float(frac_largest_cell),
        "number_of_refinement_regions": int(n_interfaces),  # approx
        "number_of_composition_interfaces": int(n_interfaces),
        "number_of_interaction_regions": int(n_turns),  # velocity turns as proxy
        "number_of_velocity_turns": int(n_turns),
        "negative_sigma_zones": int(n_neg_sigma),
        "sigma_min": float(np.min(sigma)),
        "sigma_max": float(np.max(sigma)),
        "worst_optical_depth_cell": worst_dtau_idx,
        "worst_dlogtau_cell": worst_dlogtau_idx,
        "worst_density_gradient_cell": worst_drho_idx,
        "worst_temperature_gradient_cell": worst_dT_idx,
        "worst_velocity_gradient_cell": worst_dv_idx,
        "worst_composition_interface": {
            "element": worst_comp_elem,
            "cell_index": worst_comp_idx,
            "relative_jump": float(max_comp_jump),
        },
        "composition_interface_indices": interface_indices,
        "dtau_dr_max": float(np.max(dtau_dr)),
        "dtau_dlnr_max": float(np.max(dtau_dlnr)),
        "cells_dtau_gt_1": int(np.sum(dtau_abs > 1.0)),
        "cells_dtau_gt_0.3": int(np.sum(dtau_abs > 0.3)),
        "cells_dlogtau_gt_0.1": int(np.sum(dlogtau > 0.1)),
        "mass_fractions_sum_min": float(np.min(sum(mass_fractions.values()))),
        "mass_fractions_sum_max": float(np.max(sum(mass_fractions.values()))),
    }

    # Write JSON report
    report_path = str(Path(output_file).parent / 'SN_HYDRO_DATA_grid_report.json')
    os.makedirs(os.path.dirname(report_path) or '.', exist_ok=True)
    with open(report_path, 'w') as f:
        json.dump(report, f, indent=2)
    print(f"Grid diagnostic report saved to {report_path}")

    return report


def plot_dtau_resolution_panels(r_cm, tau, output_file, time_days,
                                rho=None, T=None, v_kms=None, kappa=None):
    """Multi-panel optical-depth resolution diagnostic.

    Panels:
      1. dtau/dr vs r  — optical-depth gradient
      2. dtau/dln(r) vs r  — logarithmic optical-depth gradient
      3. Per-cell dtau (bar) — which cells dominate the tau integral
      4. Cumulative tau vs r — where optical depth accumulates

    These plots make it immediately obvious whether a single cell
    dominates the optical-depth integral due to inadequate sampling.
    """
    n = len(r_cm)
    if n < 3:
        return

    eps = 1e-30
    dtau = np.diff(tau)
    dtau_abs = np.abs(dtau)
    r_mid = 0.5 * (r_cm[:-1] + r_cm[1:])
    dr = np.abs(np.diff(r_cm))
    dr_safe = np.maximum(dr, eps)

    dtau_dr = dtau_abs / dr_safe
    dtau_dlnr = dtau_dr * r_mid

    # Per-cell fraction of total tau
    total_tau = np.sum(dtau_abs)
    frac_per_cell = dtau_abs / total_tau if total_tau > 0 else np.zeros_like(dtau_abs)

    fig, axes = plt.subplots(2, 2, figsize=(18, 14))
    fig.suptitle(f'Optical-Depth Resolution Diagnostics — t = {time_days:.1f} d, ND = {n}',
                 fontsize=14, fontweight='bold')

    # Panel 1: dtau/dr vs r
    ax = axes[0, 0]
    ax.step(r_mid, dtau_dr, 'b-o', markersize=2, linewidth=1.0, where='mid', label='dτ/dr')
    ax.set_xlabel('Radius (cm)', fontsize=11)
    ax.set_ylabel('dτ/dr (cm⁻¹)', fontsize=11, color='b')
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.grid(True, alpha=0.3, which='both')
    ax.set_title('Optical-Depth Gradient: dτ/dr', fontsize=13)
    # Mark the cell with max dtau/dr
    idx_max = int(np.argmax(dtau_dr))
    ax.annotate(f'max: cell {idx_max}\n{dtau_dr[idx_max]:.2e}',
                xy=(r_mid[idx_max], dtau_dr[idx_max]),
                fontsize=8, color='red',
                arrowprops=dict(arrowstyle='->', color='red'),
                xytext=(r_mid[idx_max] * 1.5, dtau_dr[idx_max] * 0.3))

    # Panel 2: dtau/dln(r) vs r
    ax = axes[0, 1]
    ax.step(r_mid, dtau_dlnr, 'g-o', markersize=2, linewidth=1.0, where='mid', label='dτ/dlnr')
    ax.set_xlabel('Radius (cm)', fontsize=11)
    ax.set_ylabel('dτ/dln(r)', fontsize=11, color='g')
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.grid(True, alpha=0.3, which='both')
    ax.set_title('Logarithmic Optical-Depth Gradient: dτ/dln(r)', fontsize=13)
    idx_max2 = int(np.argmax(dtau_dlnr))
    ax.annotate(f'max: cell {idx_max2}\n{dtau_dlnr[idx_max2]:.2e}',
                xy=(r_mid[idx_max2], dtau_dlnr[idx_max2]),
                fontsize=8, color='red',
                arrowprops=dict(arrowstyle='->', color='red'),
                xytext=(r_mid[idx_max2] * 1.5, dtau_dlnr[idx_max2] * 0.3))

    # Panel 3: Per-cell tau contribution (bar chart)
    ax = axes[1, 0]
    colors = ['red' if f > 0.1 else 'orange' if f > 0.05 else 'steelblue'
              for f in frac_per_cell]
    ax.bar(np.arange(len(frac_per_cell)), frac_per_cell * 100, color=colors, width=1.0)
    ax.set_xlabel('Cell Index (0=outer)', fontsize=11)
    ax.set_ylabel('% of total τ per cell', fontsize=11)
    ax.set_title('Per-Cell Optical-Depth Contribution', fontsize=13)
    ax.axhline(100.0 / (n - 1), color='gray', linestyle='--', alpha=0.5,
               label=f'Uniform ({100.0/(n-1):.2f}%)')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, axis='y')
    # Mark cells > 10%
    for i in np.where(frac_per_cell > 0.1)[0]:
        ax.annotate(f'{frac_per_cell[i]*100:.1f}%', xy=(i, frac_per_cell[i] * 100),
                    fontsize=7, color='red', ha='center', va='bottom')

    # Panel 4: Cumulative tau vs r
    ax = axes[1, 1]
    ax.step(r_cm, tau, 'k-', linewidth=1.5, where='mid', label='τ(r)')
    ax.set_xlabel('Radius (cm)', fontsize=11)
    ax.set_ylabel('Cumulative τ', fontsize=11)
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.grid(True, alpha=0.3, which='both')
    ax.set_title('Cumulative Optical Depth', fontsize=13)
    ax.axhline(2.0/3.0, color='purple', linestyle=':', alpha=0.6, label='τ = 2/3')
    ax.legend(fontsize=9)
    # Mark photosphere
    idx_phot = int(np.argmin(np.abs(tau - 2.0/3.0)))
    ax.annotate(f'τ=2/3 at r={r_cm[idx_phot]:.2e}',
                xy=(r_cm[idx_phot], 2.0/3.0),
                fontsize=8, color='purple',
                arrowprops=dict(arrowstyle='->', color='purple'),
                xytext=(r_cm[idx_phot] * 3, 2.0/3.0 * 5))

    # Depth-index top axis on the radius panels (panel 3 already uses index).
    for a in (axes[0, 0], axes[0, 1], axes[1, 1]):
        _add_depth_axis(a, r_cm)

    plt.tight_layout()
    _save_fig(fig, output_file, 'dtau_resolution')


def plot_sigma_vs_index(r_cm, v_kms, output_file, time_days, sigma_fn=None):
    """Plot sigma (dlnV/dlnR - 1) vs grid index.

    CMFGEN indexes depth from outside (index 1) to inside (index ND).
    Sigma = -1 marks the velocity minimum where dV/dR = 0 (Sobolev breakdown).
    The T1 floor at 0.1 caps the Sobolev optical depth amplification at 10x.

    Parameters
    ----------
    sigma_fn : callable, optional
        Function (r, v) -> sigma.  If None, a simple local sigma is computed.
    """
    if sigma_fn is not None:
        sigma = sigma_fn(r_cm, v_kms)
    else:
        # Fallback: simple local sigma computation
        eps = 1e-30
        dlnv = np.diff(np.log(np.maximum(v_kms, eps)))
        dlnr = np.diff(np.log(np.maximum(r_cm, eps)))
        # dlnr is negative (r decreases outer->inner).  Preserve the sign —
        # clipping negatives to +eps produces sigma ~1e28 (wrong units).
        dlnr_safe = np.where(np.abs(dlnr) < eps, eps, dlnr)
        sigma = dlnv / dlnr_safe - 1.0

    n = len(sigma)
    indices = np.arange(1, n + 1)  # 1-based to match CMFGEN depth indexing

    fig, ax = plt.subplots(figsize=(14, 7))

    ax.step(indices, sigma, 'b-', linewidth=1.2, where='mid', label='σ = dlnV/dlnR − 1')
    ax.axhline(0.0, color='gray', linestyle='-', alpha=0.5, linewidth=0.8)
    ax.axhline(-1.0, color='r', linestyle='--', alpha=0.7, linewidth=1.0,
               label='σ = −1 (velocity minimum, Sobolev breakdown)')

    # Highlight zones where |sigma + 1| < 0.1 (T1 floor active region)
    near_m1 = np.abs(sigma + 1.0) < 0.1
    if np.any(near_m1):
        ax.fill_between(indices, -1.5, 2.5, where=near_m1, alpha=0.15, color='red',
                        label=f'|σ+1| < 0.1 (T1 floor active, {np.sum(near_m1)} zones)')

    ax.set_xlabel('Grid Index (CMFGEN depth, 1=outermost)', fontsize=12)
    ax.set_ylabel('σ = dlnV/dlnR − 1', fontsize=12)
    ax.set_title(f'Sigma vs Grid Index\nt = {time_days:.1f} days, ND = {n}', fontsize=14)
    ax.legend(loc='best', fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(1, n)
    ax.invert_xaxis()  # depth ND (inner) on left, depth 1 (outer) on right

    plt.tight_layout()
    _save_fig(fig, output_file, 'sigma_vs_index')


def plot_structure_overview(r_cm, v_kms, rho, T_K, mass_fractions, kappa, time_days, output_file,
                           n_e=None, n_atom=None, ne_from_stella=False,
                           sigma_fn=None, tau_fn=None, suptitle=None,
                           ne_source_label=None):
    """Create 6-panel structure overview plot with consistent x-axis formatting.

    Input arrays are in outer->inner order (from STELLA .swd).
    For plotting, convert to inner->outer order so x-axis increases left-to-right.

    Parameters
    ----------
    sigma_fn : callable, optional
        Function (r, v) -> sigma.  If None, a simple local sigma is computed.
    tau_fn : callable, optional
        Function (r, rho, kappa) -> tau.  If None, a simple cumulative tau is computed.
    """
    # Convert to inner->outer order for plotting (x-axis increases left-to-right)
    # Note: This reverses the arrays in-place for plotting only
    r_plt = r_cm[::-1]
    v_plt = v_kms[::-1] / 1e3  # STELLA .swd velocity is naturally in 10^3 km/s
    rho_plt = rho[::-1]
    T_plt = T_K[::-1]
    kappa_plt = kappa[::-1]
    mf_plt = {elem: mass_fractions[elem][::-1] for elem in mass_fractions}

    r_plt = r_plt  # cm (no conversion)

    # Determine consistent x-axis limits
    r_min = np.min(r_plt) * 0.9
    r_max = np.max(r_plt) * 1.1

    # Create 2x4 grid
    fig, axes = plt.subplots(2, 4, figsize=(24, 12))
    fig.suptitle(suptitle or f'STELLA Structure Overview — t = {time_days:.1f} days',
                 fontsize=16, fontweight='bold')
    
    # Plot 1: Density
    ax = axes[0, 0]
    ax.step(r_plt, rho_plt, 'b-', linewidth=2, where='mid')
    ax.set_ylabel('Density (g/cm³)', fontsize=12)
    ax.set_title('Density Profile', fontsize=14)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_xlim(r_min, r_max)
    ax.set_xscale('log')
    ax.set_yscale('log')

    # Plot 2: Temperature
    ax = axes[0, 1]
    ax.step(r_plt, T_plt, 'r-', linewidth=2, where='mid')
    ax.set_ylabel('Temperature (K)', fontsize=12)
    ax.set_title('Temperature Profile', fontsize=14)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_xlim(r_min, r_max)
    ax.set_xscale('log')
    ax.set_yscale('log')

    # Plot 3: Velocity
    ax = axes[0, 2]
    ax.step(r_plt, v_plt, 'g-', linewidth=2, where='mid')
    ax.set_ylabel('Velocity (10^3 km/s)', fontsize=12)
    ax.set_title('Velocity Profile', fontsize=14)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_xlim(r_min, r_max)
    ax.set_xscale('log')
    ax.set_yscale('log')

    # Plot 4: Opacity
    ax = axes[0, 3]
    ax.step(r_plt, kappa_plt, 'm-', linewidth=2, where='mid')
    ax.set_ylabel('Opacity (cm²/g)', fontsize=12)
    ax.set_title('Rosseland Opacity', fontsize=14)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_xlim(r_min, r_max)
    ax.set_xscale('log')
    ax.set_yscale('log')
    
    # Plot 5: Mass fractions (stacked)
    ax = axes[1, 0]
    # Use all elements present in mass_fractions, skipping zero-valued ones
    # Map internal names to display names for clarity
    display_names = {
        'H': 'H', 'He': 'He', 'C': 'C', 'N': 'N', 'O': 'O', 'Ne': 'Ne',
        'Si': 'Si', 'S': 'S', 'Ca': 'Ca', 'FePeak': 'Fe', 'Fe': 'Fe',
        'NICK': 'Ni', 'COB': 'Co', 'IRON': 'Fe', 'Ni58': 'Ni', 'Ni56': 'Ni'
    }
    # Extended color palette for more elements
    colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4', '#FFEAA7', '#DDA0DD',
              '#F8B500', '#6C5CE7', '#A8E6CF', '#FFD93D', '#FF8E53', '#6A0572',
              '#AB83A1', '#C44569', '#F8B195', '#F67280', '#C06C84', '#6C5B7B']

    # Filter out elements that are all zero
    sorted_elements = sorted([elem for elem in mf_plt.keys() if np.any(mf_plt[elem] > 1e-10)])
    mf_array = np.array([mf_plt[elem] for elem in sorted_elements])
    display_labels = [display_names.get(elem, elem) for elem in sorted_elements]

    if len(sorted_elements) > 0:
        ax.stackplot(r_plt, mf_array, labels=display_labels,
                    colors=colors[:len(sorted_elements)],
                    alpha=0.8)
        ax.legend(fontsize=8, ncol=2, loc='upper right')
    ax.set_ylabel('Mass Fraction', fontsize=12)
    ax.set_title('Mass Fractions', fontsize=14)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(r_min, r_max)
    ax.set_xscale('log')

    # Plot 6: Electron and atom densities
    ax = axes[1, 1]
    if ne_source_label is None:
        ne_source_label = 'STELLA' if ne_from_stella else 'Saha (analytical)'
    if n_e is not None:
        n_e_plt = n_e[::-1]  # Convert to inner->outer order
        ax.step(r_plt, n_e_plt, 'b-', linewidth=2, where='mid', label=f'n_e [{ne_source_label}]')
    # Compute n_atom and the full-ionization electron ceiling from the same
    # composition and atomic data so the two curves are mutually consistent.
    if mass_fractions:
        n_atom_plt = _compute_atom_density(rho_plt, mf_plt)
        ne_max_plt = _compute_ne_max(rho_plt, mf_plt)
        ax.step(r_plt, n_atom_plt, 'r-', linewidth=2, where='mid', label='n_atom (cm⁻³)')
        ax.step(r_plt, ne_max_plt, 'b--', linewidth=1.5, where='mid', alpha=0.7, label='n_e (full ionization)')
    elif n_atom is not None:
        n_atom_plt = n_atom[::-1]
        ax.step(r_plt, n_atom_plt, 'r-', linewidth=2, where='mid', label='n_atom (cm⁻³)')
    ax.set_xscale('log')
    ax.set_yscale('log')
    if n_e is not None or n_atom is not None or mass_fractions:
        ax.legend(fontsize=9)
    ax.set_ylabel('Number Density (cm⁻³)', fontsize=12)
    ax.set_title(f'Electron & Atom Densities (n_e: {ne_source_label})', fontsize=13)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_xlim(r_min, r_max)

    # Plot 7: Sigma (velocity gradient) — placed directly below velocity panel
    ax = axes[1, 2]
    if sigma_fn is not None:
        # sigma_fn (e.g. sigma_velocity) returns one value per zone (n points),
        # defined at zone centers using the CMFGEN distance-weighted formula.
        sigma_plt = sigma_fn(r_cm, v_kms)[::-1]  # Convert to inner->outer order
        r_sigma = r_plt  # n points, one per zone
    else:
        # Fallback: simple logarithmic central difference (n-1 interface points).
        eps = 1e-30
        dlnv = np.diff(np.log(np.maximum(v_kms, eps)))
        dlnr = np.diff(np.log(np.maximum(r_cm, eps)))
        # dlnr is negative (r decreases outer->inner).  Preserve the sign —
        # clipping negatives to +eps produces sigma ~1e28 (wrong units).
        dlnr_safe = np.where(np.abs(dlnr) < eps, eps, dlnr)
        sigma_plt = (dlnv / dlnr_safe - 1.0)[::-1]
        r_sigma = r_plt[:-1]  # n-1 interface points
    ax.step(r_sigma, sigma_plt, 'c-', linewidth=2, where='mid')
    ax.set_ylabel('σ = dlnV/dlnr - 1', fontsize=12)
    ax.set_title('Velocity Gradient (σ)', fontsize=14)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(r_min, r_max)
    ax.set_xscale('log')

    # Plot 8: Optical depth
    ax = axes[1, 3]
    if tau_fn is not None:
        tau = tau_fn(r_cm, rho, kappa)
    else:
        # Simple cumulative tau: integrate kappa * rho * dr
        dr = np.abs(np.diff(r_cm))
        dtau = kappa[:-1] * rho[:-1] * dr
        tau = np.cumsum(dtau)
        tau = np.concatenate([[0], tau])
    tau_plt = tau[::-1]  # Convert to inner->outer order
    ax.step(r_plt, tau_plt, 'k-', linewidth=2, where='mid')
    ax.axhline(y=2/3, color='orange', linestyle='--', alpha=0.8, linewidth=2, label='τ = 2/3')
    ax.set_ylabel('Optical Depth τ', fontsize=12)
    ax.set_title('Optical Depth', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_xlim(r_min, r_max)
    ax.set_xscale('log')
    ax.set_yscale('log')

    # Set common x-axis label on bottom row
    for ax in axes[1, :]:
        ax.set_xlabel(r'Radius (cm)')

    # Depth-index top axis on every panel (all share the same radius range).
    for ax in axes.flat:
        _add_depth_axis(ax, r_cm)

    plt.tight_layout()
    _save_fig(fig, output_file, 'structure_overview')
