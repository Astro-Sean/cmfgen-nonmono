#!/usr/bin/env python3
"""Resolution convergence study for the STELLA → CMFGEN grid.

Generates grids at progressively stricter resolution levels:
  coarse   (~150 points)
  medium   (~250 points)
  fine     (~400 points)
  very_fine (~600 points)

For each level, compares:
  - tau(r)
  - dtau/dr
  - delta_log_tau
  - density
  - temperature
  - velocity
  - composition
  - location of tau=1
  - interaction-front position

Outputs a comparison report and plots.
"""

import sys
import os
import json
import numpy as np
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from convert_stella_swd_to_cmfgen import (
    read_swd, normalize_swd_composition, compute_kappa_tau,
    approximate_rosseland_opacity, sigma_velocity,
)
from tau_amr_optimiser import (
    AMRConfig, detect_all_features, recursive_amr, escalate_grid,
    pathological_cell_report,
)


def run_convergence_study(swd_file, target_days, output_dir, kappa=0.34):
    """Run the resolution convergence study."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Read STELLA data
    data = read_swd(swd_file)
    epochs = sorted(data.keys())
    best = min(epochs, key=lambda t: abs(t - target_days))
    ep = data[best]
    print(f"Epoch: {best:.2f}d, {len(ep['r'])} zones")

    r = ep['r']; v = ep['v']; T = ep['T']; rho = ep['rho']
    mf = normalize_swd_composition(ep['composition'])

    # Compute kappa (hybrid: STELLA NLTE floored at e-scattering)
    kappa_stella = ep.get('kappa', None)
    if kappa_stella is not None:
        kappa_arr = np.maximum(kappa_stella, kappa)
    else:
        kappa_arr = compute_kappa_total(rho, T, mf)
        kappa_arr = np.maximum(kappa_arr, kappa)

    tau = compute_kappa_tau(r, rho, kappa_arr)

    # Tau trim (same as converter)
    tau_const = compute_kappa_tau(r, rho, np.full(len(r), kappa))
    tau_max = 100.0; tau_min = 0.01
    i_outer = 0
    while i_outer < len(r) and tau[i_outer] < tau_min: i_outer += 1
    i_inner = len(r) - 1
    while i_inner > 0 and tau_const[i_inner] > tau_max: i_inner -= 1

    sl = slice(i_outer, i_inner + 1)
    r = r[sl]; v = v[sl]; T = T[sl]; rho = rho[sl]
    kappa_arr = kappa_arr[sl]; tau = tau[sl]
    mf = {k: val[sl] for k, val in mf.items()}

    print(f"Tau-trimmed: {len(r)} zones, r=[{r[0]:.2e}, {r[-1]:.2e}]")

    # Resolution levels
    levels = {
        'coarse':    {'max_points': 150, 'min_points': 100, 'target_dlogtau': 0.15},
        'medium':    {'max_points': 250, 'min_points': 150, 'target_dlogtau': 0.10},
        'fine':      {'max_points': 400, 'min_points': 200, 'target_dlogtau': 0.06},
        'very_fine': {'max_points': 600, 'min_points': 300, 'target_dlogtau': 0.04},
    }

    results = {}

    for level_name, params in levels.items():
        print(f"\n=== {level_name} (max={params['max_points']}, dlogtau={params['target_dlogtau']}) ===")
        cfg = AMRConfig(
            max_points=params['max_points'],
            min_points=params['min_points'],
            base_grid_points=max(40, params['min_points'] // 2),
            target_dlogtau=params['target_dlogtau'],
            kappa_source='hybrid',
            verbose=False,
        )
        feats, comp_f, shock_f, ext_f = detect_all_features(r, v, rho, T, mf, tau, cfg)
        interface_indices = np.array([f.index for f in comp_f], dtype=int)

        r_grid, v_grid, T_grid, rho_grid, kappa_grid, tau_grid, mf_grid, amr_report = \
            recursive_amr(r, v, T, rho, kappa_arr, tau, mf, feats, cfg, interface_indices)

        # Escalation
        r_src = r.copy(); v_src = v.copy(); T_src = T.copy()
        rho_src = rho.copy(); kappa_src = kappa_arr.copy(); tau_src = tau.copy()
        mf_src = {k: val.copy() for k, val in mf.items()}
        r_grid, v_grid, T_grid, rho_grid, kappa_grid, tau_grid, mf_grid, esc_report = \
            escalate_grid(r_grid, v_grid, T_grid, rho_grid, kappa_grid, tau_grid, mf_grid,
                          r_src, v_src, T_src, rho_src, kappa_src, tau_src, mf_src,
                          feats, cfg, interface_indices, max_rounds=5)

        nd = len(r_grid)
        sigma = sigma_velocity(r_grid, v_grid)
        dtau = np.abs(np.diff(tau_grid))
        dlogtau = np.abs(np.diff(np.log10(np.maximum(tau_grid, 1e-10))))

        # Find tau=1 location
        idx_tau1 = int(np.argmin(np.abs(tau_grid - 1.0)))
        r_tau1 = float(r_grid[idx_tau1])

        # Find interaction front (max |sigma| or density jump)
        idx_int = int(np.argmax(np.abs(sigma)))
        r_int = float(r_grid[idx_int])

        result = {
            'level': level_name,
            'ND': nd,
            'max_points': params['max_points'],
            'target_dlogtau': params['target_dlogtau'],
            'r_outer': float(r_grid[0]),
            'r_inner': float(r_grid[-1]),
            'tau_inner': float(tau_grid[-1]),
            'r_tau1': r_tau1,
            'r_interaction': r_int,
            'max_dlogtau': float(np.max(dlogtau)),
            'median_dlogtau': float(np.median(dlogtau)),
            'max_dtau': float(np.max(dtau)),
            'frac_largest_cell': float(np.max(dtau) / np.sum(dtau)),
            'n_velocity_turns': sum(1 for i in range(1, nd-1)
                                    if (v_grid[i] > v_grid[i-1] and v_grid[i] > v_grid[i+1]) or
                                       (v_grid[i] < v_grid[i-1] and v_grid[i] < v_grid[i+1])),
            'n_neg_sigma': int(np.sum(sigma < 0)),
            'sigma_min': float(np.min(sigma)),
            'r_grid': r_grid.tolist(),
            'tau_grid': tau_grid.tolist(),
            'v_grid': v_grid.tolist(),
            'rho_grid': rho_grid.tolist(),
            'T_grid': T_grid.tolist(),
        }
        results[level_name] = result
        print(f"  ND={nd}, max_dlogtau={result['max_dlogtau']:.3f}, "
              f"r_tau1={r_tau1:.2e}, frac_largest={result['frac_largest_cell']:.3f}")

    # Convergence analysis
    print("\n=== Convergence Analysis ===")
    levels_list = list(results.keys())
    for i in range(1, len(levels_list)):
        prev = results[levels_list[i-1]]
        curr = results[levels_list[i]]
        d_r_tau1 = abs(curr['r_tau1'] - prev['r_tau1']) / prev['r_tau1']
        d_r_int = abs(curr['r_interaction'] - prev['r_interaction']) / prev['r_interaction']
        d_max_dlogtau = (prev['max_dlogtau'] - curr['max_dlogtau']) / prev['max_dlogtau']
        d_frac = (prev['frac_largest_cell'] - curr['frac_largest_cell']) / prev['frac_largest_cell']
        print(f"  {prev['level']} → {curr['level']}: "
              f"Δr(τ=1)={d_r_tau1:.3f}, Δr(int)={d_r_int:.3f}, "
              f"Δmax_dlogtau={d_max_dlogtau:.3f}, Δfrac_largest={d_frac:.3f}")

    # Save report
    report_path = output_dir / 'convergence_report.json'
    # Strip large arrays for the summary
    summary = {k: {k2: v2 for k2, v2 in v.items() if k2 not in
                   ['r_grid', 'tau_grid', 'v_grid', 'rho_grid', 'T_grid']}
               for k, v in results.items()}
    with open(report_path, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"\nReport saved to {report_path}")

    # Generate comparison plot
    try:
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(2, 3, figsize=(20, 12))
        fig.suptitle(f'Resolution Convergence Study — {best:.1f}d', fontsize=14)

        for level_name in levels_list:
            r_grid = np.array(results[level_name]['r_grid'])
            tau_grid = np.array(results[level_name]['tau_grid'])
            v_grid = np.array(results[level_name]['v_grid'])
            rho_grid = np.array(results[level_name]['rho_grid'])
            T_grid = np.array(results[level_name]['T_grid'])
            nd = results[level_name]['ND']

            # tau vs r
            axes[0,0].plot(r_grid, tau_grid, '-o', markersize=2, label=f'{level_name} (ND={nd})')
            # v vs r
            axes[0,1].plot(r_grid, v_grid, '-o', markersize=2, label=f'{level_name}')
            # rho vs r
            axes[0,2].plot(r_grid, rho_grid, '-o', markersize=2, label=f'{level_name}')
            # T vs r
            axes[1,0].plot(r_grid, T_grid, '-o', markersize=2, label=f'{level_name}')
            # dlogtau vs index
            dlogtau = np.abs(np.diff(np.log10(np.maximum(tau_grid, 1e-10))))
            axes[1,1].plot(np.arange(len(dlogtau)), dlogtau, '-o', markersize=2, label=f'{level_name}')
            # dtau per cell
            dtau = np.abs(np.diff(tau_grid))
            axes[1,2].plot(np.arange(len(dtau)), dtau, '-o', markersize=2, label=f'{level_name}')

        for ax, ylabel, title in [
            (axes[0,0], 'τ', 'Optical Depth'),
            (axes[0,1], 'V (km/s)', 'Velocity'),
            (axes[0,2], 'ρ (g/cm³)', 'Density'),
            (axes[1,0], 'T (K)', 'Temperature'),
            (axes[1,1], 'Δlog τ', 'dlog(tau) per cell'),
            (axes[1,2], 'Δτ', 'dtau per cell'),
        ]:
            ax.set_ylabel(ylabel, fontsize=11)
            ax.set_title(title, fontsize=13)
            ax.legend(fontsize=8)
            ax.grid(True, alpha=0.3)
            if ax in [axes[0,0], axes[0,1], axes[0,2], axes[1,0]]:
                ax.set_xscale('log')
            if ax in [axes[0,0], axes[0,2], axes[1,0]]:
                ax.set_yscale('log')

        axes[1,1].set_xlabel('Cell Index')
        axes[1,2].set_xlabel('Cell Index')

        plt.tight_layout()
        plot_path = output_dir / 'convergence_comparison.png'
        plt.savefig(plot_path, dpi=150, bbox_inches='tight')
        print(f"Comparison plot saved to {plot_path}")
        plt.close()
    except ImportError:
        print("matplotlib not available, skipping plot")


if __name__ == '__main__':
    swd = sys.argv[1] if len(sys.argv) > 1 else '/viper/ptmp1/sbrennan/SN2024pba/X/stella/strad/M3.0_V7000.swd'
    days = float(sys.argv[2]) if len(sys.argv) > 2 else 35.57
    out = sys.argv[3] if len(sys.argv) > 3 else '/tmp/convergence_study'
    run_convergence_study(swd, days, out)
