#!/usr/bin/env python3
"""Gradient diagnostic for the wind-ejecta stitching comparison.

Plots dlog(rho)/dlog(r) and dlog(v)/dlog(r) zone-by-zone for both stitching
options, with the solver-stall threshold (|grad| > 16) marked.  This makes
the numerical-stability comparison explicit.
"""
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
from test_density_matched_wind import (
    read_source_ejecta, create_hhe_shell, create_wind,
    combine_shell_wind_ejecta, combine_shell_wind_ejecta_density_matched,
    MSUN,
)


def gradients(grid, key):
    r = grid['r']
    y = np.maximum(grid[key], 1e-300)
    log_r = np.log10(r)
    log_y = np.log10(y)
    return np.diff(log_y) / np.diff(log_r), 0.5*(r[:-1] + r[1:])


def main():
    # Same params as test_density_matched_wind.py
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
    wind_mdot_current = 0.01
    wind_velocity_kms = 100.0
    wind_T_K = 2000
    sw_transition_dex = 0.5
    we_transition_dex = 0.5

    out_dir = Path('/viper/ptmp1/sbrennan/SN2024pba/X_nomono/plots')

    print("Reading source ejecta...")
    ejecta = read_source_ejecta(source_dir, source_age=source_age,
                                t_start_days=t_start)
    ejecta['v'] = ejecta['v'] * 1e5
    if ejecta_rho_min > 0:
        mask = ejecta['rho'] >= ejecta_rho_min
        n_before = len(ejecta['r'])
        for key in list(ejecta.keys()):
            if isinstance(ejecta[key], np.ndarray) and len(ejecta[key]) == n_before:
                ejecta[key] = ejecta[key][mask]

    print("Creating shell...")
    shell = create_hhe_shell(shell_mass, shell_velocity,
                             shell_xh, shell_yhe,
                             n_rho=n_rho, v0_kms=v0_kms,
                             n_T=n_T, T_interface_K=T_interface,
                             T_max_K=T_max, t_start_days=t_start,
                             composition=shell_composition)

    print("Building CURRENT profile...")
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

    print("Building DENSITY-MATCHED profile...")
    grid_matched, _, diag_matched = combine_shell_wind_ejecta_density_matched(
        shell, ejecta, nzon_total=nzon,
        wind_velocity_kms=wind_velocity_kms,
        wind_T_K=wind_T_K,
        composition=shell_composition,
        ejecta_T_floor=ejecta_T_floor,
        transition_zones=5)

    r_int = diag_matched['r_interface']

    # --- Plot ---
    fig, axes = plt.subplots(2, 2, figsize=(14, 10), sharex=True)
    fig.suptitle('Wind-Ejecta stitching: gradient diagnostic\n'
                 '|dlog(rho)/dlog(r)| and |dlog(v)/dlog(r)| '
                 '(solver-stall threshold = 16)',
                 fontsize=12, y=0.995)

    for col, (grid, lab, c) in enumerate([
        (grid_current, 'CURRENT (log-linear buffer)', 'C0'),
        (grid_matched, 'DENSITY-MATCHED (CD sigmoid)', 'C3'),
    ]):
        # rho gradient
        grad_rho, r_mid = gradients(grid, 'rho')
        ax = axes[0, col]
        ax.plot(r_mid, np.abs(grad_rho), '-', color=c, lw=0.7)
        ax.axhline(16, color='k', ls='--', alpha=0.5, lw=1.0, label='stall threshold (16)')
        ax.axvline(r_int, color='gray', ls=':', alpha=0.6, lw=1.0)
        ax.set_xscale('log'); ax.set_yscale('log')
        ax.set_ylabel('|dlog rho / dlog r|')
        ax.set_title(lab, fontsize=10)
        ax.set_ylim(1e-1, 2e2)
        ax.grid(True, which='both', alpha=0.3)
        if col == 0:
            ax.legend(fontsize=8, loc='upper left')

        # v gradient
        grad_v, r_mid = gradients(grid, 'v')
        ax = axes[1, col]
        ax.plot(r_mid, np.abs(grad_v), '-', color=c, lw=0.7)
        ax.axhline(16, color='k', ls='--', alpha=0.5, lw=1.0)
        ax.axvline(r_int, color='gray', ls=':', alpha=0.6, lw=1.0)
        ax.set_xscale('log'); ax.set_yscale('log')
        ax.set_ylabel('|dlog v / dlog r|')
        ax.set_xlabel('r  [cm]')
        ax.set_ylim(1e-2, 1e2)
        ax.grid(True, which='both', alpha=0.3)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    out = out_dir / 'wind_ejecta_stitching_gradients.png'
    plt.savefig(out, dpi=140, bbox_inches='tight')
    print(f"Wrote gradient plot: {out}")
    plt.close(fig)

    # --- Print worst-gradient locations ---
    print("\nWorst |dlog rho/dlog r| locations:")
    for grid, lab in [(grid_current, 'CURRENT'), (grid_matched, 'MATCHED')]:
        g, r_mid = gradients(grid, 'rho')
        idx = np.argmax(np.abs(g))
        print(f"  {lab}: |grad|={abs(g[idx]):.1f} at r={r_mid[idx]:.3e} cm "
              f"(zone ~{idx})")
        # Worst in the wind/interface region only (between shell outer and ejecta inner+1dex)
        r_shell_outer = shell['r_outer']
        r_ej_inner = ejecta['r'].min() * 0.99
        m = (r_mid > r_shell_outer * 0.5) & (r_mid < r_ej_inner * 3)
        if m.any():
            g_wind = np.abs(g[m])
            r_wind = r_mid[m]
            idx_w = np.argmax(g_wind)
            print(f"    wind/interface region: |grad|={g_wind[idx_w]:.2f} "
                  f"at r={r_wind[idx_w]:.3e} cm")


if __name__ == '__main__':
    main()
