#!/usr/bin/env python3
"""
Test/evaluate the density-matched wind-ejecta stitching option.

Builds two profiles using the SAME source ejecta and shell parameters
(as in X_nomono/config.yaml):

  A) CURRENT  : combine_shell_wind_ejecta  (log-linear density buffer,
                 fixed Mdot, wind composition held across the whole buffer)
  B) DENSITY-MATCHED : new combine_shell_wind_ejecta_density_matched
                 (Mdot derived so rho_wind(r_interface) = rho_ejecta(r_interface),
                  thin sigmoid transition for v/T/composition, no extended buffer)

Then produces comparison plots so the user can evaluate before any change
is made to the production pipeline.

This is a DIAGNOSTIC script only -- it does NOT write .hyd/.abn/.dat and does
NOT submit any STELLA job.  It only reads the source .swd and prints/plots.

Usage:
    python scripts/test_density_matched_wind.py
"""
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# Reuse the production helpers
sys.path.insert(0, str(Path(__file__).parent))
from prepare_SN_B_model import (
    read_source_ejecta,
    create_hhe_shell,
    create_wind,
    combine_shell_wind_ejecta,
    _interface_blend, _interface_ramp_weights,
    MSUN, YR_TO_S, SOLAR_COMPOSITION,
)


# ----------------------------------------------------------------------------
# New: density-matched wind-ejecta combiner
# ----------------------------------------------------------------------------
def combine_shell_wind_ejecta_density_matched(
    shell, ejecta, nzon_total=500,
    wind_velocity_kms=100.0, wind_T_K=2000.0,
    composition='solar', ejecta_T_floor=2000.0,
    transition_zones=5, sw_transition_dex=0.5,
):
    """Density-matched contact-discontinuity stitching at the wind-ejecta
    interface, with a log-linear shell-wind buffer preserved.

    The shell-wind boundary is a genuine physical discontinuity (dense shell
    running into rarefied wind), so a log-linear density buffer is kept there
    (same as the current combiner).

    At the wind-ejecta interface, Mdot is derived so that
    rho_wind(r_interface) = rho_ejecta(r_interface), i.e. the wind and ejecta
    are in pressure/density balance at the contact radius.  This eliminates
    the density cliff that motivates the extended log-linear buffer in the
    current combiner.  A thin sigmoid (transition_zones, default 5) blends
    velocity, temperature and composition across the contact discontinuity.

    Returns (combined_grid, wind_dict, diagnostics).
    """
    from scipy.interpolate import interp1d

    # Sort ejecta inner->outer
    ej_order = np.argsort(ejecta['r'])
    ej_r_sorted = ejecta['r'][ej_order]
    ej_rho_sorted = ejecta['rho'][ej_order]
    ej_v_sorted = ejecta['v'][ej_order]
    ej_T_sorted = ejecta['T'][ej_order]

    # Region boundaries
    r_min = shell['r_inner']
    r_shell_outer = shell['r_outer']
    r_ej_inner = ej_r_sorted[0] * 0.99  # small overlap, as in production
    r_max = ej_r_sorted[-1]

    log_r_min = np.log10(r_min)
    log_r_shell = np.log10(r_shell_outer)
    log_r_ej = np.log10(r_ej_inner)
    log_r_max = np.log10(r_max)

    v_wind = wind_velocity_kms * 1e5  # cm/s
    T_wind = wind_T_K

    # --- DERIVE Mdot from density matching at the contact radius ---
    log_r_ej_src = np.log10(ej_r_sorted)
    log_rho_ej_src = np.log10(np.maximum(ej_rho_sorted, 1e-300))
    rho_ej_at_interface = 10.0 ** np.interp(
        np.log10(r_ej_inner), log_r_ej_src, log_rho_ej_src)

    mdot_matched = 4.0 * np.pi * r_ej_inner**2 * v_wind * rho_ej_at_interface
    mdot_msun_yr = mdot_matched / (MSUN / YR_TO_S)

    print(f"\n  DENSITY-MATCHED wind-ejecta combiner:")
    print(f"    r_shell_outer = {r_shell_outer:.3e} cm")
    print(f"    r_ej_inner    = {r_ej_inner:.3e} cm")
    print(f"    rho_ejecta(r_interface) = {rho_ej_at_interface:.3e} g/cc")
    print(f"    => derived Mdot = {mdot_msun_yr:.4e} Msun/yr "
          f"(v_wind={wind_velocity_kms} km/s)")

    # --- Shell-wind buffer (same logic as current combiner) ---
    # The shell edge density is much higher than the wind density at that
    # radius, so a log-linear buffer is needed to keep gradients < 16.
    total_gap_dex = log_r_ej - log_r_shell
    max_buffer_dex = total_gap_dex / 3.0
    sw_transition_dex = min(sw_transition_dex, max_buffer_dex)

    r_sw_buffer_end = r_shell_outer * (10 ** sw_transition_dex)
    rho_sw_start = shell['rho'][-1]
    rho_sw_end = mdot_matched / (4.0 * np.pi * r_sw_buffer_end**2 * v_wind)
    sw_contrast = np.log10(rho_sw_start / max(rho_sw_end, 1e-300))
    sw_grad = sw_contrast / sw_transition_dex
    nzon_sw_buffer = max(20, int(sw_contrast / 16 * sw_transition_dex / 0.02))
    nzon_sw_buffer = min(nzon_sw_buffer, 60)

    print(f"    Shell-wind buffer: {nzon_sw_buffer} zones, "
          f"r=[{r_shell_outer:.3e}, {r_sw_buffer_end:.3e}] cm")
    print(f"      rho: {rho_sw_start:.3e} -> {rho_sw_end:.3e} g/cc "
          f"(contrast={sw_contrast:.1f} dex, est. gradient={sw_grad:.1f})")

    # Zone allocation
    nzon_shell = max(50, int(nzon_total * 0.15))
    nzon_ejecta = max(100, int(nzon_total * 0.55))
    nzon_wind = nzon_total - nzon_shell - nzon_ejecta - nzon_sw_buffer - transition_zones
    nzon_wind = max(40, nzon_wind)
    n_total = (nzon_shell + nzon_sw_buffer + nzon_wind + transition_zones
               + nzon_ejecta)
    print(f"    Zone allocation: {nzon_shell} shell + {nzon_sw_buffer} SW buffer "
          f"+ {nzon_wind} wind + {transition_zones} CD-sigmoid + "
          f"{nzon_ejecta} ejecta = {n_total}")

    # --- Shell grid ---
    shell_r = np.logspace(log_r_min, log_r_shell, nzon_shell)
    log_r_src = np.log10(shell['r'])
    log_rho_src = np.log10(np.maximum(shell['rho'], 1e-300))
    shell_rho = 10.0 ** np.interp(np.log10(shell_r), log_r_src, log_rho_src)
    shell_v = np.interp(np.log10(shell_r), log_r_src, shell['v'])
    shell_T = np.interp(np.log10(shell_r), log_r_src, shell['T'])

    # Shell mass conservation
    if 'rho_0' in shell and 'r0' in shell:
        n_rho_val = 8.0
        frac_v0 = shell['r0'] / shell['r_outer']
        factor = 1 + (3.0 / (n_rho_val - 3)) * (1 - frac_v0**(n_rho_val - 3))
        M_shell_target = (4/3) * np.pi * shell['r0']**3 * shell['rho_0'] * factor
    else:
        r_lo_src = np.concatenate([[shell['r_inner']], shell['r'][:-1]])
        M_shell_target = ((4/3) * np.pi *
                          (shell['r']**3 - r_lo_src**3) * shell['rho']).sum()
    r_lo_shell = np.concatenate([[shell_r[0] * 0.5], shell_r[:-1]])
    M_shell_interp = ((4/3) * np.pi *
                      (shell_r**3 - r_lo_shell**3) * shell_rho).sum()
    if M_shell_interp > 0:
        shell_rho *= M_shell_target / M_shell_interp

    # --- Shell-wind buffer (log-linear density ramp) ---
    sw_buf_r = np.logspace(log_r_shell, np.log10(r_sw_buffer_end),
                           nzon_sw_buffer + 1)[1:]
    sw_log_rho = np.linspace(np.log10(rho_sw_start),
                             np.log10(max(rho_sw_end, 1e-300)),
                             nzon_sw_buffer)
    sw_buf_rho = 10.0 ** sw_log_rho
    sw_buf_v = np.linspace(shell['v'][-1], v_wind, nzon_sw_buffer)
    sw_buf_T = np.linspace(shell['T'][-1], T_wind, nzon_sw_buffer)

    # --- Wind grid (1/r^2 with derived Mdot) ---
    # Wind spans from r_sw_buffer_end to r_ej_inner, with the last
    # transition_zones points being the CD sigmoid.
    wind_r = np.logspace(np.log10(r_sw_buffer_end), log_r_ej,
                         nzon_wind + transition_zones + 1)[1:]
    wind_r_pure = wind_r[:nzon_wind]
    cd_r = wind_r[nzon_wind:]

    wind_rho_pure = mdot_matched / (4.0 * np.pi * wind_r_pure**2 * v_wind)
    wind_v_pure = np.full(nzon_wind, v_wind)
    wind_T_pure = np.full(nzon_wind, T_wind)

    # --- Contact-discontinuity sigmoid (thin, transition_zones points) ---
    # Density: matched by construction -- wind 1/r^2 law continues through
    # the CD points, and at r_ej_inner rho_wind = rho_ejecta by Mdot choice.
    cd_rho = mdot_matched / (4.0 * np.pi * cd_r**2 * v_wind)

    v_ej_inner = ej_v_sorted[0]
    cd_v = _interface_blend(transition_zones, v_wind, v_ej_inner,
                            transition_zones=transition_zones)
    cd_T = _interface_blend(transition_zones, T_wind,
                            max(ej_T_sorted[0], ejecta_T_floor),
                            transition_zones=transition_zones)

    # --- Ejecta grid (log-spaced, interpolated from source) ---
    ej_r_new = np.logspace(log_r_ej, log_r_max, nzon_ejecta + 1)[1:]

    def log_interp(y, r_target):
        log_y = np.log10(np.maximum(y, 1e-300))
        ip = interp1d(log_r_ej_src, log_y, kind='linear',
                      bounds_error=False,
                      fill_value=(log_y[0], log_y[-1]))
        return 10.0 ** ip(np.log10(r_target))

    def lin_interp(y, r_target):
        ip = interp1d(log_r_ej_src, y, kind='linear',
                      bounds_error=False,
                      fill_value=(y[0], y[-1]))
        return ip(np.log10(r_target))

    ej_rho = log_interp(ej_rho_sorted, ej_r_new)
    ej_v = lin_interp(ej_v_sorted, ej_r_new)
    ej_T = np.maximum(log_interp(ej_T_sorted, ej_r_new), ejecta_T_floor)

    # Light T smoothing (same as production)
    if len(ej_T) > 5:
        kernel = np.ones(5) / 5.0
        ej_T_smooth = np.convolve(ej_T, kernel, mode='same')
        ej_T_smooth[:2] = ej_T[:2]
        ej_T_smooth[-2:] = ej_T[-2:]
        ej_T = np.maximum(ej_T_smooth, ejecta_T_floor)

    # Light rho smoothing (same as production)
    for _ in range(3):
        log_ej_rho = np.log10(np.maximum(ej_rho, 1e-300))
        log_ej_r = np.log10(ej_r_new)
        grads = np.abs(np.diff(log_ej_rho) / np.diff(log_ej_r))
        smoothed = ej_rho.copy()
        for i in range(1, len(ej_rho) - 1):
            g_left = grads[i - 1]
            g_right = grads[i]
            if g_left > 16 or g_right > 16:
                smoothed[i] = np.sqrt(max(ej_rho[i - 1] * ej_rho[i + 1], 1e-300))
        ej_rho = smoothed

    # --- Combine ---
    r_combined = np.concatenate([shell_r, sw_buf_r, wind_r_pure, cd_r, ej_r_new])
    rho_combined = np.concatenate([shell_rho, sw_buf_rho, wind_rho_pure, cd_rho, ej_rho])
    v_combined = np.concatenate([shell_v, sw_buf_v, wind_v_pure, cd_v, ej_v])
    T_combined = np.concatenate([shell_T, sw_buf_T, wind_T_pure, cd_T, ej_T])

    # --- Composition ---
    all_elem_keys = ['fH', 'fHe', 'fO', 'fSi', 'fCa', 'fFe']
    ej_comp = {}
    for key in all_elem_keys:
        if key in ejecta:
            ej_comp[key] = lin_interp(ejecta[key][ej_order], ej_r_new)
        else:
            ej_comp[key] = np.zeros(nzon_ejecta)

    # Shell composition
    shell_comp_vals = {k: float(shell.get(k, [0.0])[0]) for k in all_elem_keys}

    # Wind composition (solar or H/He)
    if composition == 'solar':
        sc = SOLAR_COMPOSITION
        wind_comp_vals = {k: sc[k] for k in all_elem_keys}
    else:
        wind_comp_vals = {'fH': 0.7, 'fHe': 0.3,
                          'fO': 0.0, 'fSi': 0.0, 'fCa': 0.0, 'fFe': 0.0}

    # CD sigmoid composition: blend wind -> ejecta inner-edge composition
    ej_comp_inner = {}
    for key in all_elem_keys:
        if key in ejecta:
            ej_comp_inner[key] = float(ejecta[key][ej_order][0])
        else:
            ej_comp_inner[key] = 0.0

    cd_w = _interface_ramp_weights(transition_zones, transition_zones=transition_zones)

    comp_combined = {}
    for key in all_elem_keys:
        sv = shell_comp_vals[key]
        wv = wind_comp_vals[key]
        ev = ej_comp_inner[key]
        cd_comp = wv * (1.0 - cd_w) + ev * cd_w
        comp_combined[key] = np.concatenate([
            np.full(nzon_shell, sv),
            np.full(nzon_sw_buffer, sv),   # SW buffer = shell composition
            np.full(nzon_wind, wv),
            cd_comp,
            ej_comp[key],
        ])

    # Renormalize
    comp_sum = sum(comp_combined.values())
    comp_sum = np.maximum(comp_sum, 1e-10)
    for key in comp_combined:
        comp_combined[key] /= comp_sum

    # --- Diagnostics ---
    r_lo = np.concatenate([[r_combined[0] * 0.5], r_combined[:-1]])
    dm = (4/3) * np.pi * (r_combined**3 - r_lo**3) * rho_combined
    dm = np.maximum(dm, 0.0)
    M_total = np.sum(dm)

    idx_shell = nzon_shell
    idx_sw_buf = idx_shell + nzon_sw_buffer
    idx_wind = idx_sw_buf + nzon_wind
    idx_cd = idx_wind + transition_zones

    print(f"    Mass: shell={np.sum(dm[:idx_shell])/MSUN:.3f} Msun, "
          f"SW buf={np.sum(dm[idx_shell:idx_sw_buf])/MSUN:.4f} Msun, "
          f"wind={np.sum(dm[idx_sw_buf:idx_wind])/MSUN:.4f} Msun, "
          f"CD={np.sum(dm[idx_wind:idx_cd])/MSUN:.6f} Msun, "
          f"ejecta={np.sum(dm[idx_cd:])/MSUN:.3f} Msun, "
          f"total={M_total/MSUN:.3f} Msun")

    # Density gradients
    log_r = np.log10(r_combined)
    log_rho_arr = np.log10(np.maximum(rho_combined, 1e-300))
    if len(log_r) > 1:
        gradients = np.abs(np.diff(log_rho_arr) / np.diff(log_r))
        max_grad = np.max(gradients)
        grad_interior = gradients[:-5]
        max_grad_interior = np.max(grad_interior) if len(grad_interior) > 0 else 0
        print(f"    Max |dlog rho/dlog r|: {max_grad:.1f} "
              f"(interior: {max_grad_interior:.1f})")
        # Gradients across the CD sigmoid
        cd_grads = gradients[idx_wind-1:idx_cd]
        if len(cd_grads) > 0:
            print(f"    CD sigmoid gradients: max={np.max(cd_grads):.2f}, "
                  f"mean={np.mean(cd_grads):.2f}")
        # Gradients across the shell-wind buffer
        sw_grads = gradients[idx_shell-1:idx_sw_buf]
        if len(sw_grads) > 0:
            print(f"    SW buffer gradients: max={np.max(sw_grads):.2f}, "
                  f"mean={np.mean(sw_grads):.2f}")

    wind_dict = {
        'r': wind_r_pure, 'rho': wind_rho_pure, 'v': wind_v_pure, 'T': wind_T_pure,
        'r_inner': r_sw_buffer_end, 'r_outer': r_ej_inner,
        'M_wind': np.sum(dm[idx_sw_buf:idx_wind]),
        'mdot_msun_yr': mdot_msun_yr,
        'mdot_cgs': mdot_matched,
    }

    diagnostics = {
        'mdot_msun_yr': mdot_msun_yr,
        'rho_ej_at_interface': rho_ej_at_interface,
        'rho_wind_at_interface': mdot_matched / (4.0 * np.pi * r_ej_inner**2 * v_wind),
        'r_interface': r_ej_inner,
        'v_wind_cms': v_wind,
        'T_wind': T_wind,
        'M_wind_msun': np.sum(dm[idx_sw_buf:idx_wind]) / MSUN,
        'M_total_msun': M_total / MSUN,
        'nzon_shell': nzon_shell,
        'nzon_sw_buffer': nzon_sw_buffer,
        'nzon_wind': nzon_wind,
        'nzon_cd': transition_zones,
        'nzon_ejecta': nzon_ejecta,
    }

    combined = {
        'r': r_combined, 'rho': rho_combined,
        'v': v_combined, 'T': T_combined,
        **comp_combined,
    }
    return combined, wind_dict, diagnostics


# ----------------------------------------------------------------------------
# Plotting
# ----------------------------------------------------------------------------
def plot_comparison(grid_current, grid_matched, diag_current, diag_matched,
                    out_path):
    """Side-by-side comparison of the two stitching options."""
    fig, axes = plt.subplots(4, 2, figsize=(14, 16), sharex='col')
    fig.suptitle('Wind-Ejecta stitching: current log-linear buffer vs '
                 'density-matched contact discontinuity',
                 fontsize=13, y=0.995)

    colors = {'current': 'C0', 'matched': 'C3'}
    labels = {
        'current': f"Current (Mdot={diag_current['mdot_msun_yr']:.2e} Msun/yr, "
                   f"buffer zones)",
        'matched': f"Density-matched (Mdot={diag_matched['mdot_msun_yr']:.2e} "
                   f"Msun/yr, {diag_matched['nzon_cd']}-zone CD sigmoid)",
    }

    for col, (grid, diag, key) in enumerate([
        (grid_current, diag_current, 'current'),
        (grid_matched, diag_matched, 'matched'),
    ]):
        r = grid['r']
        rho = grid['rho']
        v = grid['v'] / 1e5  # km/s
        T = grid['T']
        c = colors[key]
        lab = labels[key]

        # Row 0: density
        ax = axes[0, col]
        ax.plot(r, rho, '-', color=c, lw=0.8, label=lab)
        ax.set_xscale('log'); ax.set_yscale('log')
        ax.set_ylabel('rho  [g/cc]')
        ax.set_title(lab, fontsize=10)
        # Mark the interface
        r_int = diag.get('r_interface', None)
        if r_int is not None:
            ax.axvline(r_int, color='k', ls='--', alpha=0.4, lw=0.8)
            ax.text(r_int, ax.get_ylim()[1]*0.3, ' r_interface',
                    fontsize=8, color='k', alpha=0.6, va='top')
        ax.grid(True, which='both', alpha=0.3)

        # Row 1: velocity
        ax = axes[1, col]
        ax.plot(r, v, '-', color=c, lw=0.8)
        ax.set_xscale('log'); ax.set_yscale('log')
        ax.set_ylabel('v  [km/s]')
        if r_int is not None:
            ax.axvline(r_int, color='k', ls='--', alpha=0.4, lw=0.8)
        ax.grid(True, which='both', alpha=0.3)

        # Row 2: temperature
        ax = axes[2, col]
        ax.plot(r, T, '-', color=c, lw=0.8)
        ax.set_xscale('log'); ax.set_yscale('log')
        ax.set_ylabel('T  [K]')
        if r_int is not None:
            ax.axvline(r_int, color='k', ls='--', alpha=0.4, lw=0.8)
        ax.grid(True, which='both', alpha=0.3)

        # Row 3: composition
        ax = axes[3, col]
        for elem, c_elem in [('fH', 'C2'), ('fHe', 'C4'),
                             ('fO', 'C1'), ('fSi', 'C8'),
                             ('fCa', 'C9'), ('fFe', 'C5')]:
            ax.plot(r, grid[elem], '-', color=c_elem, lw=0.8, label=elem)
        ax.set_xscale('log'); ax.set_yscale('log')
        ax.set_ylabel('mass fraction')
        ax.set_xlabel('r  [cm]')
        if col == 0:
            ax.legend(fontsize=8, loc='lower left', ncol=2)
        if r_int is not None:
            ax.axvline(r_int, color='k', ls='--', alpha=0.4, lw=0.8)
        ax.grid(True, which='both', alpha=0.3)
        ax.set_ylim(1e-6, 2)

    plt.tight_layout(rect=[0, 0, 1, 0.985])
    plt.savefig(out_path, dpi=140, bbox_inches='tight')
    print(f"\n  Wrote comparison plot: {out_path}")
    plt.close(fig)


def plot_overlay(grid_current, grid_matched, diag_current, diag_matched,
                 out_path):
    """Overlay both profiles on the same axes for direct comparison."""
    fig, axes = plt.subplots(4, 1, figsize=(11, 16), sharex=True)
    fig.suptitle('Wind-Ejecta stitching: overlay '
                 '(blue=current log-linear buffer, red=density-matched CD)',
                 fontsize=12, y=0.995)

    r_int = diag_matched['r_interface']

    # Density
    ax = axes[0]
    ax.plot(grid_current['r'], grid_current['rho'], '-', color='C0',
            lw=0.9, label=f"current (Mdot={diag_current['mdot_msun_yr']:.2e})")
    ax.plot(grid_matched['r'], grid_matched['rho'], '-', color='C3',
            lw=0.9, label=f"density-matched (Mdot={diag_matched['mdot_msun_yr']:.2e})")
    ax.axvline(r_int, color='k', ls='--', alpha=0.5, lw=0.8,
               label=f'r_interface={r_int:.2e} cm')
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_ylabel('rho  [g/cc]')
    ax.legend(fontsize=9)
    ax.grid(True, which='both', alpha=0.3)

    # Velocity
    ax = axes[1]
    ax.plot(grid_current['r'], grid_current['v']/1e5, '-', color='C0', lw=0.9)
    ax.plot(grid_matched['r'], grid_matched['v']/1e5, '-', color='C3', lw=0.9)
    ax.axvline(r_int, color='k', ls='--', alpha=0.5, lw=0.8)
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_ylabel('v  [km/s]')
    ax.grid(True, which='both', alpha=0.3)

    # Temperature
    ax = axes[2]
    ax.plot(grid_current['r'], grid_current['T'], '-', color='C0', lw=0.9)
    ax.plot(grid_matched['r'], grid_matched['T'], '-', color='C3', lw=0.9)
    ax.axvline(r_int, color='k', ls='--', alpha=0.5, lw=0.8)
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_ylabel('T  [K]')
    ax.grid(True, which='both', alpha=0.3)

    # Composition (overlay both -- H/He/Fe only to keep it readable)
    ax = axes[3]
    for elem, c_elem in [('fH', 'C2'), ('fHe', 'C4'),
                         ('fO', 'C1'), ('fFe', 'C5')]:
        ax.plot(grid_current['r'], grid_current[elem], '-', color=c_elem,
                lw=0.7, alpha=0.55)
        ax.plot(grid_matched['r'], grid_matched[elem], '-', color=c_elem,
                lw=1.1, label=f"{elem} (matched)")
    ax.axvline(r_int, color='k', ls='--', alpha=0.5, lw=0.8)
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_ylabel('mass fraction')
    ax.set_xlabel('r  [cm]')
    ax.legend(fontsize=8, loc='lower left', ncol=2)
    ax.grid(True, which='both', alpha=0.3)
    ax.set_ylim(1e-6, 2)

    plt.tight_layout(rect=[0, 0, 1, 0.985])
    plt.savefig(out_path, dpi=140, bbox_inches='tight')
    print(f"  Wrote overlay plot: {out_path}")
    plt.close(fig)


def plot_zoom_interface(grid_current, grid_matched, diag_current, diag_matched,
                        out_path):
    """Zoom on the wind-ejecta interface region."""
    r_int = diag_matched['r_interface']
    # Zoom window: 1 dex below to 1 dex above r_interface
    r_lo = r_int / 10**1.0
    r_hi = r_int * 10**1.0

    fig, axes = plt.subplots(3, 1, figsize=(11, 12), sharex=True)
    fig.suptitle(f'Zoom on wind-ejecta interface (r_interface={r_int:.2e} cm)\n'
                 f'blue=current log-linear buffer, red=density-matched CD',
                 fontsize=11, y=0.995)

    for grid, diag, c, lab in [
        (grid_current, diag_current, 'C0',
         f"current (Mdot={diag_current['mdot_msun_yr']:.2e}, "
         f"M_wind={diag_current.get('M_wind_msun', float('nan')):.4f} Msun)"),
        (grid_matched, diag_matched, 'C3',
         f"density-matched (Mdot={diag_matched['mdot_msun_yr']:.2e}, "
         f"M_wind={diag_matched['M_wind_msun']:.4f} Msun)"),
    ]:
        r = grid['r']
        m = (r >= r_lo) & (r <= r_hi)
        ax = axes[0]
        ax.plot(r[m], grid['rho'][m], '-', color=c, lw=1.4, label=lab)
        ax = axes[1]
        ax.plot(r[m], grid['v'][m]/1e5, '-', color=c, lw=1.4)
        ax = axes[2]
        ax.plot(r[m], grid['T'][m], '-', color=c, lw=1.4)

    for ax, ylabel in [(axes[0], 'rho  [g/cc]'),
                       (axes[1], 'v  [km/s]'),
                       (axes[2], 'T  [K]')]:
        ax.axvline(r_int, color='k', ls='--', alpha=0.5, lw=0.8)
        ax.set_xscale('log'); ax.set_yscale('log')
        ax.set_ylabel(ylabel)
        ax.grid(True, which='both', alpha=0.3)
    axes[0].legend(fontsize=9)
    axes[2].set_xlabel('r  [cm]')

    plt.tight_layout(rect=[0, 0, 1, 0.985])
    plt.savefig(out_path, dpi=140, bbox_inches='tight')
    print(f"  Wrote interface zoom plot: {out_path}")
    plt.close(fig)


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main():
    # Use the SAME parameters as X_nomono/config.yaml
    source_dir = '/viper/ptmp1/sbrennan/SN2024pba/SN_A/SN_A_g_T20k_Ni06'
    source_age = 400.0
    shell_mass = 3.0
    shell_velocity = 7000.0
    shell_xh = 0.7381
    shell_yhe = 0.2485
    shell_composition = 'solar'
    n_rho = 8.0
    v0_kms = 3000.0
    n_T = 0.4
    T_interface = 2000.0
    T_max = 65000.0
    t_start = 5.0
    ejecta_T_floor = 2000.0
    ejecta_rho_min = 1.0e-20
    nzon = 500

    # Wind parameters (current X_nomono config)
    wind_mdot_current = 0.01        # Msun/yr (current config)
    wind_velocity_kms = 100.0
    wind_T_K = 2000
    sw_transition_dex = 0.5
    we_transition_dex = 0.5

    out_dir = Path('/viper/ptmp1/sbrennan/SN2024pba/X_nomono/plots')
    out_dir.mkdir(parents=True, exist_ok=True)

    print("\n" + "="*70)
    print("TEST: density-matched wind-ejecta stitching")
    print("="*70)

    # ---- Step 1: read source ejecta (same for both) ----
    print("\n--- Step 1: Reading source ejecta ---")
    ejecta = read_source_ejecta(source_dir, source_age=source_age,
                                t_start_days=t_start)
    print(f"  Ejecta: {len(ejecta['r'])} zones, "
          f"r=[{ejecta['r'].min():.3e}, {ejecta['r'].max():.3e}] cm")
    print(f"  Ejecta: rho=[{ejecta['rho'].min():.3e}, {ejecta['rho'].max():.3e}] g/cc")
    ej_inner_idx = int(np.argmin(ejecta['r']))
    print(f"  Ejecta inner edge: r={ejecta['r'].min():.3e} cm, "
          f"rho={ejecta['rho'][ej_inner_idx]:.3e} g/cc")

    # Convert ejecta velocities to cm/s (production does this too)
    ejecta['v'] = ejecta['v'] * 1e5

    # Apply rho_min truncation (same as production)
    if ejecta_rho_min > 0:
        mask = ejecta['rho'] >= ejecta_rho_min
        n_before = len(ejecta['r'])
        n_truncated = n_before - mask.sum()
        if n_truncated > 0:
            for key in list(ejecta.keys()):
                if isinstance(ejecta[key], np.ndarray) and len(ejecta[key]) == n_before:
                    ejecta[key] = ejecta[key][mask]
            print(f"  Ejecta truncation (rho_min={ejecta_rho_min:.1e}): "
                  f"{n_before} -> {mask.sum()} zones")

    # ---- Step 2: build the shell (same for both, NOT immediate contact) ----
    print("\n--- Step 2: Creating H/He shell (Dessart 2010jl structure) ---")
    shell = create_hhe_shell(shell_mass, shell_velocity,
                             shell_xh, shell_yhe,
                             n_rho=n_rho, v0_kms=v0_kms,
                             n_T=n_T, T_interface_K=T_interface,
                             T_max_K=T_max, t_start_days=t_start,
                             composition=shell_composition)
    print(f"  Shell: r=[{shell['r_inner']:.3e}, {shell['r_outer']:.3e}] cm, "
          f"rho_edge={shell['rho'][-1]:.3e} g/cc")

    # ---- Step 3a: CURRENT profile (fixed Mdot, log-linear buffer) ----
    print("\n--- Step 3a: CURRENT profile (fixed Mdot, log-linear buffer) ---")
    wind_r_inner = shell['r_outer'] * 1.01
    wind_r_outer = ejecta['r'].min() * 0.99
    wind_current = create_wind(wind_mdot_current, wind_velocity_kms, wind_T_K,
                               wind_r_inner, wind_r_outer,
                               composition=shell_composition)
    grid_current = combine_shell_wind_ejecta(
        shell, wind_current, ejecta, nzon_total=nzon,
        sw_transition_dex=sw_transition_dex,
        we_transition_dex=we_transition_dex,
        ejecta_T_floor=ejecta_T_floor)

    # Diagnostics for current
    r_lo = np.concatenate([[grid_current['r'][0]*0.5], grid_current['r'][:-1]])
    dm = (4/3) * np.pi * (grid_current['r']**3 - r_lo**3) * grid_current['rho']
    dm = np.maximum(dm, 0.0)
    diag_current = {
        'mdot_msun_yr': wind_mdot_current,
        'r_interface': wind_r_outer,
        'M_wind_msun': wind_current['M_wind'] / MSUN,
        'M_total_msun': np.sum(dm) / MSUN,
    }
    print(f"  Current: M_wind={diag_current['M_wind_msun']:.4f} Msun, "
          f"M_total={diag_current['M_total_msun']:.3f} Msun")

    # ---- Step 3b: DENSITY-MATCHED profile ----
    print("\n--- Step 3b: DENSITY-MATCHED profile ---")
    grid_matched, wind_matched, diag_matched = \
        combine_shell_wind_ejecta_density_matched(
            shell, ejecta, nzon_total=nzon,
            wind_velocity_kms=wind_velocity_kms,
            wind_T_K=wind_T_K,
            composition=shell_composition,
            ejecta_T_floor=ejecta_T_floor,
            transition_zones=5)

    # ---- Step 4: Plots ----
    print("\n--- Step 4: Generating comparison plots ---")
    plot_comparison(grid_current, grid_matched, diag_current, diag_matched,
                    out_dir / 'wind_ejecta_stitching_comparison.png')
    plot_overlay(grid_current, grid_matched, diag_current, diag_matched,
                 out_dir / 'wind_ejecta_stitching_overlay.png')
    plot_zoom_interface(grid_current, grid_matched, diag_current, diag_matched,
                        out_dir / 'wind_ejecta_stitching_zoom.png')

    # ---- Summary ----
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    ej_inner_idx = int(np.argmin(ejecta['r']))
    print(f"  Source ejecta inner edge: r={ejecta['r'].min():.3e} cm, "
          f"rho={ejecta['rho'][ej_inner_idx]:.3e} g/cc")
    print(f"  Shell outer edge:         r={shell['r_outer']:.3e} cm, "
          f"rho={shell['rho'][-1]:.3e} g/cc")
    print()
    print(f"  CURRENT (fixed Mdot, log-linear buffer):")
    print(f"    Mdot                = {diag_current['mdot_msun_yr']:.4e} Msun/yr")
    print(f"    M_wind              = {diag_current['M_wind_msun']:.4f} Msun")
    print(f"    M_total             = {diag_current['M_total_msun']:.3f} Msun")
    print()
    print(f"  DENSITY-MATCHED (derived Mdot, thin CD sigmoid):")
    print(f"    Mdot                = {diag_matched['mdot_msun_yr']:.4e} Msun/yr")
    print(f"    rho_ej(interface)  = {diag_matched['rho_ej_at_interface']:.3e} g/cc")
    print(f"    rho_wind(interface)= {diag_matched['rho_wind_at_interface']:.3e} g/cc")
    print(f"    M_wind              = {diag_matched['M_wind_msun']:.4f} Msun")
    print(f"    M_total             = {diag_matched['M_total_msun']:.3f} Msun")
    print(f"    Zones: shell={diag_matched['nzon_shell']}, "
          f"wind={diag_matched['nzon_wind']}, "
          f"CD={diag_matched['nzon_cd']}, "
          f"ejecta={diag_matched['nzon_ejecta']}")
    print()
    print(f"  Plots written to {out_dir}/")
    print(f"    wind_ejecta_stitching_comparison.png  (side-by-side)")
    print(f"    wind_ejecta_stitching_overlay.png     (overlay)")
    print(f"    wind_ejecta_stitching_zoom.png        (interface zoom)")
    print("="*70)


if __name__ == '__main__':
    main()
