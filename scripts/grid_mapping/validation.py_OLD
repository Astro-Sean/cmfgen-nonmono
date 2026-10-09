"""Validation: conservation, ordering, feature, and solver checks.

Validates the target CMFGEN grid before writing production files.
All checks return a list of issues (empty = pass).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from .shell_model import SourceShellGrid, TargetGrid
from .conservative_remap import compute_conservation_residuals
from .solver_metrics import sigma_velocity, apply_sigma_cap, count_velocity_reversals
from .feature_detection import FeatureReport


@dataclass
class ValidationReport:
    """Comprehensive validation report."""
    passed: bool = True
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    conservation_residuals: dict = field(default_factory=dict)
    ordering_check: dict = field(default_factory=dict)
    feature_check: dict = field(default_factory=dict)
    solver_check: dict = field(default_factory=dict)
    finite_check: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        return {
            'passed': self.passed,
            'errors': self.errors,
            'warnings': self.warnings,
            'conservation_residuals': self.conservation_residuals,
            'ordering_check': self.ordering_check,
            'feature_check': self.feature_check,
            'solver_check': self.solver_check,
            'finite_check': self.finite_check,
        }


def validate_target_grid(target: TargetGrid,
                         source: SourceShellGrid,
                         source_features: FeatureReport,
                         conservation_tolerance: float = 1e-3,
                         sigma_tolerance: float = 1e-6) -> ValidationReport:
    """Full validation of the target grid.

    Checks:
    1. Radius ordering (strictly decreasing, no duplicates)
    2. All arrays finite and equal length
    3. Physical signs (density > 0, temperature > 0)
    4. SIGMA consistency (recomputed from R, V matches stored sigma)
    5. Velocity reversal count preserved
    6. Feature anchors retained
    7. Conservation residuals within tolerance
    8. Mass fractions sum to 1
    """
    report = ValidationReport()
    n = target.n_points

    # --- 1. Radius ordering ---
    r = target.r_cm
    if n < 2:
        report.errors.append("Target grid has fewer than 2 points")
        report.passed = False
        return report

    diffs = np.diff(r)
    n_nondecreasing = int(np.sum(diffs >= 0))
    n_duplicate = int(np.sum(np.abs(diffs) < 1e-30))
    report.ordering_check = {
        'n_points': n,
        'strictly_decreasing': n_nondecreasing == 0,
        'n_non_decreasing': n_nondecreasing,
        'n_duplicate': n_duplicate,
    }
    if n_nondecreasing > 0:
        report.errors.append(f"Radius not strictly decreasing ({n_nondecreasing} violations)")
        report.passed = False

    # --- 2. Array lengths and finiteness ---
    array_names = ['r_cm', 'velocity_cms', 'temperature_K', 'density_gcm3',
                   'kappa_cgs', 'sigma', 'tau_cumulative']
    length_issues = []
    finite_issues = []
    for name in array_names:
        arr = getattr(target, name)
        if len(arr) != n:
            length_issues.append(f"{name}: {len(arr)} != {n}")
        if not np.all(np.isfinite(arr)):
            n_bad = int(np.sum(~np.isfinite(arr)))
            finite_issues.append(f"{name}: {n_bad} non-finite")
    report.finite_check = {
        'length_issues': length_issues,
        'finite_issues': finite_issues,
    }
    if length_issues:
        report.errors.extend(length_issues)
        report.passed = False
    if finite_issues:
        report.errors.extend(finite_issues)
        report.passed = False

    # --- 3. Physical signs ---
    if np.any(target.density_gcm3 <= 0):
        n_bad = int(np.sum(target.density_gcm3 <= 0))
        report.errors.append(f"density has {n_bad} non-positive values")
        report.passed = False
    if np.any(target.temperature_K <= 0):
        n_bad = int(np.sum(target.temperature_K <= 0))
        report.errors.append(f"temperature has {n_bad} non-positive values")
        report.passed = False

    # --- 4. SIGMA consistency ---
    sigma_recomputed = sigma_velocity(target.r_cm, target.velocity_cms, preserve_sign=True)
    sigma_diff = np.max(np.abs(sigma_recomputed - target.sigma))
    report.solver_check['sigma_max_diff'] = float(sigma_diff)
    if sigma_diff > sigma_tolerance:
        report.warnings.append(f"SIGMA recomputed differs from stored (max diff={sigma_diff:.2e})")

    # SIGMA cap check
    sigma_capped = apply_sigma_cap(target.sigma)
    n_capped = int(np.sum(sigma_capped != target.sigma))
    report.solver_check['n_sigma_capped'] = n_capped

    # --- 5. Velocity reversal count ---
    source_reversals = count_velocity_reversals(
        sigma_velocity(source.r_centre_cm, source.velocity_cms, preserve_sign=True))
    target_reversals = count_velocity_reversals(sigma_recomputed)
    report.solver_check['source_velocity_reversals'] = source_reversals
    report.solver_check['target_velocity_reversals'] = target_reversals
    if target_reversals < source_reversals:
        report.warnings.append(
            f"Velocity reversal count changed: {source_reversals} -> {target_reversals}")

    # --- 6. Feature anchors ---
    if target.is_anchor is not None:
        n_anchors = int(np.sum(target.is_anchor))
    else:
        n_anchors = 0
    report.feature_check = {
        'n_source_anchors': len(source_features.anchor_indices),
        'n_target_anchors': n_anchors,
    }
    if n_anchors < len(source_features.anchor_indices):
        report.warnings.append(
            f"Anchor count reduced: {len(source_features.anchor_indices)} -> {n_anchors}")

    # --- 7. Conservation residuals ---
    residuals = compute_conservation_residuals(source, target)
    report.conservation_residuals = residuals

    for key, val in residuals.items():
        if isinstance(val, dict) and 'fractional_residual' in val:
            if abs(val['fractional_residual']) > conservation_tolerance:
                report.warnings.append(
                    f"Conservation residual {key}: {val['fractional_residual']:.2e} "
                    f"(tolerance={conservation_tolerance:.0e})")

    # --- 8. Mass fraction sums ---
    if target.mass_fractions:
        mf_sum = np.zeros(n)
        for frac in target.mass_fractions.values():
            if len(frac) == n:
                mf_sum += frac
        max_dev = float(np.max(np.abs(mf_sum - 1.0)))
        report.finite_check['mass_fraction_max_deviation'] = max_dev
        if max_dev > 0.01:
            report.warnings.append(f"Mass fractions deviate from 1.0 (max={max_dev:.4f})")

    return report
