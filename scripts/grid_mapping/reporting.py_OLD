"""Reporting: JSON, text, plots, and provenance diagnostics.

Writes the required diagnostic files:
- grid_audit.json
- grid_features.json
- grid_validation.json
- grid_validation.txt
- original and final profile plots
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Optional

import numpy as np

from .shell_model import SourceShellGrid, TargetGrid
from .feature_detection import FeatureReport
from .validation import ValidationReport


def write_json(filepath: str, data: dict):
    """Write JSON with atomic rename."""
    dir_path = os.path.dirname(filepath) or '.'
    fd, tmp_path = tempfile_mkstemp(dir_path)
    try:
        with os.fdopen(fd, 'w') as f:
            json.dump(data, f, indent=2, default=str)
        os.rename(tmp_path, filepath)
    except Exception:
        if os.path.exists(tmp_path):
            os.unlink(tmp_path)
        raise


def tempfile_mkstemp(dir_path: str):
    import tempfile
    return tempfile.mkstemp(dir=dir_path, suffix='.tmp')


def write_grid_audit(filepath: str, source: SourceShellGrid,
                     target: TargetGrid, config: dict):
    """Write grid_audit.json with source and target grid summary."""
    audit = {
        'source': {
            'n_zones': source.n_shells,
            'n_zones_original': source.n_zones_original,
            'time_days': source.time_days,
            'source_file': source.source_file,
            'r_range_cm': (float(source.r_bound_cm[0]), float(source.r_bound_cm[-1])),
            'v_range_cms': (float(source.velocity_cms.min()), float(source.velocity_cms.max())),
            'T_range_K': (float(source.temperature_K.min()), float(source.temperature_K.max())),
            'rho_range': (float(source.density_gcm3.min()), float(source.density_gcm3.max())),
            'total_mass_msun': float(source.shell_mass_g.sum() / 1.989e33),
            'tau_range': (float(source.tau_cumulative[0]), float(source.tau_cumulative[-1])),
            'has_ne': source.electron_density is not None,
            'has_ni56': source.ni56_initial is not None,
            'n_mass_fractions': len(source.mass_fractions),
            'n_isotopes': len(source.isotope_fractions),
        },
        'target': {
            'n_points': target.n_points,
            'time_days': target.time_days,
            'r_range_cm': (float(target.r_cm[0]), float(target.r_cm[-1])),
            'v_range_cms': (float(target.velocity_cms.min()), float(target.velocity_cms.max())),
            'T_range_K': (float(target.temperature_K.min()), float(target.temperature_K.max())),
            'rho_range': (float(target.density_gcm3.min()), float(target.density_gcm3.max())),
            'sigma_range': (float(target.sigma.min()), float(target.sigma.max())),
            'tau_range': (float(target.tau_cumulative[0]), float(target.tau_cumulative[-1])),
            'has_ne': target.electron_density is not None,
            'has_ni56': target.ni56_initial is not None,
            'n_mass_fractions': len(target.mass_fractions),
            'n_isotopes': len(target.isotope_fractions),
            'n_anchors': int(target.is_anchor.sum()) if target.is_anchor is not None else 0,
            'provenance_counts': _provenance_counts(target),
        },
        'config': config,
    }
    write_json(filepath, audit)


def _provenance_counts(target: TargetGrid) -> dict:
    """Count provenance types."""
    if len(target.provenance) == 0:
        return {}
    unique, counts = np.unique(target.provenance, return_counts=True)
    return {str(u): int(c) for u, c in zip(unique, counts)}


def write_grid_features(filepath: str, features: FeatureReport):
    """Write grid_features.json with detected feature summary."""
    write_json(filepath, features.to_dict())


def write_grid_validation_json(filepath: str, report: ValidationReport):
    """Write grid_validation.json."""
    write_json(filepath, report.to_dict())


def write_grid_validation_txt(filepath: str, report: ValidationReport,
                              source: SourceShellGrid, target: TargetGrid):
    """Write grid_validation.txt (human-readable validation report)."""
    lines = []
    lines.append("=" * 60)
    lines.append("GRID VALIDATION REPORT")
    lines.append("=" * 60)
    lines.append("")
    lines.append(f"Source grid: {source.n_shells} shells ({source.n_zones_original} original)")
    lines.append(f"Target grid: {target.n_points} points")
    lines.append(f"Time: {target.time_days:.3f} days")
    lines.append("")

    # Overall status
    status = "PASS" if report.passed else "FAIL"
    lines.append(f"Overall status: {status}")
    lines.append(f"Errors: {len(report.errors)}")
    lines.append(f"Warnings: {len(report.warnings)}")
    lines.append("")

    # Errors
    if report.errors:
        lines.append("ERRORS:")
        for e in report.errors:
            lines.append(f"  - {e}")
        lines.append("")

    # Warnings
    if report.warnings:
        lines.append("WARNINGS:")
        for w in report.warnings:
            lines.append(f"  - {w}")
        lines.append("")

    # Ordering
    lines.append("Ordering check:")
    for k, v in report.ordering_check.items():
        lines.append(f"  {k}: {v}")
    lines.append("")

    # Finite check
    lines.append("Finite check:")
    for k, v in report.finite_check.items():
        lines.append(f"  {k}: {v}")
    lines.append("")

    # Solver check
    lines.append("Solver check:")
    for k, v in report.solver_check.items():
        lines.append(f"  {k}: {v}")
    lines.append("")

    # Feature check
    lines.append("Feature check:")
    for k, v in report.feature_check.items():
        lines.append(f"  {k}: {v}")
    lines.append("")

    # Conservation residuals
    lines.append("Conservation residuals:")
    for key, val in report.conservation_residuals.items():
        if isinstance(val, dict) and 'fractional_residual' in val:
            lines.append(f"  {key}: source={val['source']:.6e}, target={val['target']:.6e}, "
                        f"residual={val['fractional_residual']:.2e}")
    lines.append("")

    with open(filepath, 'w') as f:
        f.write('\n'.join(lines) + '\n')


def write_profile_plots(plot_dir: str, source: SourceShellGrid,
                        target: TargetGrid, features: FeatureReport):
    """Write diagnostic plots comparing source and target grids."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("[Reporting] matplotlib not available, skipping plots")
        return

    plot_path = Path(plot_dir)
    plot_path.mkdir(parents=True, exist_ok=True)

    # Profile comparison plot
    fig, axes = plt.subplots(4, 2, figsize=(14, 16))

    r_src = source.r_centre_cm / 1e10
    r_tgt = target.r_cm / 1e10

    # Radius
    axes[0, 0].semilogy(r_src, 'b-', label='Source', alpha=0.5)
    axes[0, 0].semilogy(r_tgt, 'r.', label='Target', markersize=2)
    axes[0, 0].set_ylabel('R (1e10 cm)')
    axes[0, 0].set_title('Radius')
    axes[0, 0].legend()

    # Velocity
    v_src = source.velocity_cms / 1e5
    v_tgt = target.velocity_cms / 1e5
    axes[0, 1].plot(r_src, v_src, 'b-', label='Source', alpha=0.5)
    axes[0, 1].plot(r_tgt, v_tgt, 'r.', label='Target', markersize=2)
    axes[0, 1].set_ylabel('V (km/s)')
    axes[0, 1].set_title('Velocity (non-monotonic preserved)')
    axes[0, 1].legend()

    # Density
    axes[1, 0].loglog(r_src, source.density_gcm3, 'b-', label='Source', alpha=0.5)
    axes[1, 0].loglog(r_tgt, target.density_gcm3, 'r.', label='Target', markersize=2)
    axes[1, 0].set_ylabel('rho (g/cm^3)')
    axes[1, 0].set_title('Density')
    axes[1, 0].legend()

    # Temperature
    axes[1, 1].semilogy(r_src, source.temperature_K, 'b-', label='Source', alpha=0.5)
    axes[1, 1].semilogy(r_tgt, target.temperature_K, 'r.', label='Target', markersize=2)
    axes[1, 1].set_ylabel('T (K)')
    axes[1, 1].set_title('Temperature')
    axes[1, 1].legend()

    # Sigma
    from .solver_metrics import sigma_velocity
    sigma_src = sigma_velocity(source.r_centre_cm, source.velocity_cms, preserve_sign=True)
    axes[2, 0].plot(r_src, sigma_src, 'b-', label='Source', alpha=0.5)
    axes[2, 0].plot(r_tgt, target.sigma, 'r.', label='Target', markersize=2)
    axes[2, 0].axhline(y=-1.0, color='k', linestyle='--', alpha=0.3, label='sigma=-1')
    axes[2, 0].set_ylabel('sigma')
    axes[2, 0].set_title('SIGMA (dlnV/dlnR - 1)')
    axes[2, 0].legend()

    # Tau
    axes[2, 1].semilogy(r_src, source.tau_cumulative, 'b-', label='Source', alpha=0.5)
    axes[2, 1].semilogy(r_tgt, target.tau_cumulative, 'r.', label='Target', markersize=2)
    axes[2, 1].set_ylabel('tau')
    axes[2, 1].set_title('Cumulative optical depth')
    axes[2, 1].legend()

    # Mass fractions (H, He)
    for elem in ['H', 'He']:
        if elem in source.mass_fractions and elem in target.mass_fractions:
            axes[3, 0].semilogy(r_src, source.mass_fractions[elem], 'b-', alpha=0.5, label=f'{elem} src')
            axes[3, 0].semilogy(r_tgt, target.mass_fractions[elem], 'r.', markersize=2, label=f'{elem} tgt')
    axes[3, 0].set_ylabel('mass fraction')
    axes[3, 0].set_title('Composition (H, He)')
    axes[3, 0].legend()

    # Feature markers
    axes[3, 1].plot(r_src, np.zeros(len(r_src)), 'b-', alpha=0.3)
    for idx in features.shock_indices:
        if idx < len(r_src):
            axes[3, 1].axvline(x=r_src[idx], color='r', alpha=0.5, label='shock' if idx == features.shock_indices[0] else '')
    for idx in features.velocity_reversal_indices:
        if idx < len(r_src):
            axes[3, 1].axvline(x=r_src[idx], color='g', alpha=0.5, label='reversal' if idx == features.velocity_reversal_indices[0] else '')
    axes[3, 1].set_ylabel('features')
    axes[3, 1].set_title('Feature markers')
    axes[3, 1].legend()

    for ax_row in axes:
        for ax in ax_row:
            ax.set_xlabel('zone index' if ax == axes[0, 0] else 'R (1e10 cm)')

    plt.tight_layout()
    plt.savefig(plot_path / 'grid_comparison.png', dpi=150)
    plt.close()

    print(f"[Reporting] Plots written to {plot_dir}")
