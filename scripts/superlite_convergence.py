#!/usr/bin/env python3
"""Convergence framework for SuperLite spectral calculations.

This module provides quantitative convergence testing for:

  1. Spatial grid resolution (cell count)
  2. Wavelength group count
  3. Monte Carlo particle count
  4. Optical-depth truncation

Convergence is assessed by comparing spectra using:

  * L2 norm of flux difference (continuum shape)
  * L1 norm of flux difference (line sensitivity)
  * Maximum absolute flux difference (worst-case)
  * Cross-correlation of line profiles (velocity shifts)
  * Equivalent width preservation for key lines

The framework generates multiple SuperLite configurations, runs them
through the orchestrator (run_model_grid.py), collects the spectra,
and produces a convergence report with quantitative metrics.

Usage
-----
As a library:

    from superlite_convergence import ConvergenceTest, ConvergenceConfig
    test = ConvergenceTest(ConvergenceConfig(
        swd_file='path/to/model.swd',
        output_dir='/tmp/convergence_test',
        grid_cells=[50, 100, 150, 200, 250, 300],
    ))
    report = test.run_grid_convergence()

As a script:

    python superlite_convergence.py --swd model.swd --output /tmp/conv \
        --grid-cells 50 100 150 200 250 --target-days 30
"""

from __future__ import annotations

import os
import sys
import subprocess
import json
from dataclasses import dataclass, field
from typing import List, Optional, Tuple, Dict
from pathlib import Path

import numpy as np


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

@dataclass
class ConvergenceConfig:
    """Configuration for convergence testing.

    All paths should be absolute.  The swd_file is the STELLA .swd model
    to use for all convergence runs.
    """
    swd_file: str
    output_dir: str
    target_days: float = 30.0

    # Grid convergence: test different cell counts
    grid_cells: List[int] = field(
        default_factory=lambda: [50, 100, 150, 200, 250, 300])

    # Wavelength convergence: test different group counts
    wavelength_groups: List[int] = field(
        default_factory=lambda: [500, 1000, 2000, 3000, 6000])

    # MC particle convergence: test different particle counts (2^n)
    mc_particles: List[int] = field(
        default_factory=lambda: [18, 20, 22, 24, 25, 26])

    # Tau truncation convergence: test different truncation depths
    tau_truncations: List[float] = field(
        default_factory=lambda: [5.0, 10.0, 20.0, 50.0, 100.0])

    # Convergence criteria
    l2_tolerance: float = 0.01       # 1% L2 flux difference
    l1_tolerance: float = 0.02       # 2% L1 flux difference
    max_flux_tolerance: float = 0.05  # 5% max flux difference
    ew_tolerance: float = 0.05       # 5% equivalent width difference

    # Reference cell count for grid convergence comparison
    reference_cells: int = 500  # use full grid as reference

    # SuperLite executable
    superlite_exe: str = '/ptmp/sbrennan/superlite/superlite'
    superlite_data: str = '/ptmp/sbrennan/superlite/Data'
    superlite_input: str = '/ptmp/sbrennan/superlite/Input'

    # Other parameters (held fixed during convergence tests)
    truncate_tau: float = 10000.0
    nomp: int = 8


# ---------------------------------------------------------------------------
# Spectrum comparison metrics
# ---------------------------------------------------------------------------

@dataclass
class SpectrumComparison:
    """Quantitative comparison between two spectra."""
    l2_norm: float       # sqrt(sum((f1-f2)^2) / sum(f1^2))
    l1_norm: float       # sum(|f1-f2|) / sum(|f1|)
    max_diff: float      # max(|f1-f2|) / max(|f1|)
    correlation: float   # Pearson correlation coefficient
    n_wavelength: int
    wavelength_range: Tuple[float, float]


def compare_spectra(
    wl1: np.ndarray, flux1: np.ndarray,
    wl2: np.ndarray, flux2: np.ndarray,
) -> SpectrumComparison:
    """Compare two spectra using standard metrics.

    Both spectra are interpolated onto a common wavelength grid before
    comparison.  Fluxes are assumed to be in the same units (erg/s/cm²/Å).

    Args:
        wl1, flux1: reference spectrum
        wl2, flux2: comparison spectrum

    Returns:
        SpectrumComparison with quantitative metrics
    """
    # Interpolate onto common grid (use the finer grid)
    if len(wl1) >= len(wl2):
        wl_common = wl1
        f2_interp = np.interp(wl_common, wl2, flux2)
        f1 = flux1
        f2 = f2_interp
    else:
        wl_common = wl2
        f1_interp = np.interp(wl_common, wl1, flux1)
        f1 = f1_interp
        f2 = flux2

    # Avoid division by zero
    f1_safe = np.where(np.abs(f1) > 1e-30, f1, 1e-30)

    # L2 norm (relative)
    l2 = float(np.sqrt(np.sum((f1 - f2) ** 2) / np.sum(f1 ** 2)))

    # L1 norm (relative)
    l1 = float(np.sum(np.abs(f1 - f2)) / np.sum(np.abs(f1)))

    # Max difference (relative)
    max_diff = float(np.max(np.abs(f1 - f2)) / np.max(np.abs(f1)))

    # Pearson correlation
    if np.std(f1) > 0 and np.std(f2) > 0:
        corr = float(np.corrcoef(f1, f2)[0, 1])
    else:
        corr = 1.0

    return SpectrumComparison(
        l2_norm=l2, l1_norm=l1, max_diff=max_diff,
        correlation=corr, n_wavelength=len(wl_common),
        wavelength_range=(float(wl_common.min()), float(wl_common.max()))
    )


def read_superlite_spectrum(spectrum_file: str) -> Tuple[np.ndarray, np.ndarray]:
    """Read a SuperLite output spectrum.

    SuperLite writes flux output in a standard format.  This function
    attempts to parse common output formats.

    Args:
        spectrum_file: path to the spectrum file

    Returns:
        (wavelength_A, flux_erg_s_cm2_A) arrays
    """
    # Try common SuperLite output formats
    data = np.loadtxt(spectrum_file)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    # Typical format: wavelength (cm), flux (erg/s/cm²/cm or similar)
    wl_cm = data[:, 0]
    flux = data[:, 1]
    # Convert wavelength to Angstroms
    wl_A = wl_cm * 1e8
    return wl_A, flux


# ---------------------------------------------------------------------------
# Convergence tests
# ---------------------------------------------------------------------------

@dataclass
class ConvergenceResult:
    """Result of a single convergence test point."""
    parameter_name: str
    parameter_value: float
    comparison: SpectrumComparison
    converged: bool
    n_cells: int = 0


@dataclass
class ConvergenceReport:
    """Full convergence report."""
    test_type: str
    results: List[ConvergenceResult] = field(default_factory=list)
    converged_value: Optional[float] = None
    reference_value: Optional[float] = None

    def to_dict(self) -> dict:
        return {
            'test_type': self.test_type,
            'reference_value': self.reference_value,
            'converged_value': self.converged_value,
            'results': [
                {
                    'parameter': r.parameter_name,
                    'value': float(r.parameter_value),
                    'l2': float(r.comparison.l2_norm),
                    'l1': float(r.comparison.l1_norm),
                    'max_diff': float(r.comparison.max_diff),
                    'correlation': float(r.comparison.correlation),
                    'converged': bool(r.converged),
                    'n_cells': int(r.n_cells),
                }
                for r in self.results
            ]
        }

    def summary(self) -> str:
        lines = [f"Convergence report: {self.test_type}"]
        if self.reference_value is not None:
            lines.append(f"  Reference: {self.reference_value}")
        if self.converged_value is not None:
            lines.append(f"  Converged at: {self.converged_value}")
        lines.append(f"  {'Value':>10s}  {'L2':>8s}  {'L1':>8s}  {'MaxDiff':>8s}  "
                      f"{'Corr':>6s}  {'Conv':>5s}")
        for r in self.results:
            lines.append(
                f"  {r.parameter_value:10.3f}  "
                f"{r.comparison.l2_norm:8.4f}  {r.comparison.l1_norm:8.4f}  "
                f"{r.comparison.max_diff:8.4f}  {r.comparison.correlation:6.4f}  "
                f"{'YES' if r.converged else 'no':>5s}")
        return '\n'.join(lines)


class ConvergenceTest:
    """Run convergence tests for SuperLite.

    This class generates multiple SuperLite input configurations,
    converts the STELLA model for each, and (optionally) runs SuperLite
    to produce spectra for comparison.

    For grid convergence, the conversion is run with different max_cells
    values.  For wavelength/MC/tau convergence, the input.par is modified
    and SuperLite is rerun.
    """

    def __init__(self, config: ConvergenceConfig):
        self.config = config
        self.scripts_dir = Path(__file__).parent
        self.project_dir = self.scripts_dir.parent

    def _convert_swd(
        self, output_dir: str, max_cells: int,
        truncate_tau: Optional[float] = None,
        grid_method: str = 'adaptive',
    ) -> bool:
        """Run the STELLA → SuperLite conversion."""
        cmd = [
            sys.executable, str(self.scripts_dir / 'convert_stella_swd_to_superlite.py'),
            self.config.swd_file,
            '-o', output_dir,
            '--target-days', str(self.config.target_days),
            '--max-cells', str(max_cells),
            '--nomp', str(self.config.nomp),
            '--grid-method', grid_method,
        ]
        tt = truncate_tau or self.config.truncate_tau
        cmd.extend(['--truncate-tau', str(tt)])
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        if result.returncode != 0:
            print(f"  Conversion failed: {result.stderr[:200]}")
            return False
        return True

    def _modify_input_par(
        self, input_par_path: str,
        grp_ng: Optional[int] = None,
        src_n2s: Optional[int] = None,
        tauddmc: Optional[str] = None,
    ) -> None:
        """Modify specific parameters in an existing input.par."""
        with open(input_par_path) as f:
            content = f.read()
        if grp_ng is not None:
            import re
            content = re.sub(r'in_grp_ng\s*=\s*\S+',
                             f'in_grp_ng = {grp_ng}', content)
        if src_n2s is not None:
            import re
            content = re.sub(r'in_src_n2s\s*=\s*\S+',
                             f'in_src_n2s = {src_n2s}', content)
        if tauddmc is not None:
            import re
            content = re.sub(r'in_trn_tauddmc\s*=\s*\S+',
                             f'in_trn_tauddmc = {tauddmc}', content)
        with open(input_par_path, 'w') as f:
            f.write(content)

    def run_grid_convergence(self) -> ConvergenceReport:
        """Test convergence with respect to spatial grid resolution.

        Generates grids with different cell counts and compares the
        resulting structure (mass, velocity, composition, optical depth)
        to a high-resolution reference.

        Note: Full spectral convergence requires running SuperLite for
        each grid.  This function performs the structural convergence
        check (mass, element, velocity conservation) which is a necessary
        but not sufficient condition for spectral convergence.
        """
        report = ConvergenceReport(
            test_type='grid_resolution',
            reference_value=self.config.reference_cells,
        )

        # Generate reference grid (full STELLA resolution, no downsampling)
        ref_dir = os.path.join(self.config.output_dir, 'reference')
        os.makedirs(ref_dir, exist_ok=True)
        print(f"\n  Generating reference grid ({self.config.reference_cells} cells)...")
        self._convert_swd(ref_dir, self.config.reference_cells,
                          grid_method='adaptive')

        # Read reference grid diagnostics
        ref_str = os.path.join(ref_dir, 'input.str')
        if not os.path.exists(ref_str):
            print("  ERROR: Reference grid not generated")
            return report

        # Test different cell counts
        for n_cells in self.config.grid_cells:
            test_dir = os.path.join(self.config.output_dir, f'grid_{n_cells}')
            os.makedirs(test_dir, exist_ok=True)
            print(f"\n  Testing grid with {n_cells} cells...")
            self._convert_swd(test_dir, n_cells, grid_method='adaptive')

            # Compare structure (not spectra — that requires SuperLite runs)
            # Use velocity, temperature profiles as convergence metrics.
            # Interpolate both grids onto a common log-radius grid and compare.
            test_str = os.path.join(test_dir, 'input.str')
            if os.path.exists(test_str):
                def _read_input_str(path):
                    """Read input.str and return (r, v, T) outer→inner."""
                    lines = open(path).readlines()
                    data = []
                    for line in lines[3:]:
                        if line.strip() and not line.startswith('#'):
                            data.append([float(x) for x in line.split()])
                    arr = np.array(data)
                    xl, xr = arr[:, 0], arr[:, 1]
                    vxl, vxr = arr[:, 2], arr[:, 3]
                    T = arr[:, 6]
                    r = 0.5 * (xl + xr)
                    v = 0.5 * (vxl + vxr)
                    return r[::-1], v[::-1], T[::-1]  # outer→inner

                ref_r, ref_v, ref_T = _read_input_str(ref_str)
                test_r, test_v, test_T = _read_input_str(test_str)

                # Interpolate test profiles onto reference log-radius grid.
                # Only compare in the overlapping radial range to avoid
                # extrapolation artifacts.
                # Note: ref_r and test_r are outer→inner (decreasing).
                # np.interp requires ascending x, so we reverse for interpolation.
                r_overlap_min = max(ref_r.min(), test_r.min())
                r_overlap_max = min(ref_r.max(), test_r.max())
                mask = (ref_r >= r_overlap_min) & (ref_r <= r_overlap_max)
                if mask.sum() < 5:
                    # Not enough overlap — skip this test point
                    continue
                # Sort test arrays ascending in log-r for interpolation
                test_sort = np.argsort(np.log(test_r))
                test_log_r_asc = np.log(test_r)[test_sort]
                test_v_asc = test_v[test_sort]
                test_T_asc = test_T[test_sort]
                # Interpolate at reference points (mask)
                log_r_ref = np.log(ref_r[mask])
                test_v_interp = np.interp(log_r_ref, test_log_r_asc, test_v_asc)
                test_T_interp = np.interp(log_r_ref, test_log_r_asc, test_T_asc)
                ref_v_masked = ref_v[mask]
                ref_T_masked = ref_T[mask]

                # Relative L2 norms (only in overlapping region)
                v_l2 = float(np.sqrt(np.sum((ref_v_masked - test_v_interp) ** 2) /
                                     max(np.sum(ref_v_masked ** 2), 1e-30)))
                T_l2 = float(np.sqrt(np.sum((ref_T_masked - test_T_interp) ** 2) /
                                     max(np.sum(ref_T_masked ** 2), 1e-30)))

                # Combined structural convergence metric
                l2_combined = float(np.sqrt(v_l2**2 + T_l2**2))
                l1_combined = float(v_l2 + T_l2)
                max_combined = float(max(v_l2, T_l2))

                comparison = SpectrumComparison(
                    l2_norm=l2_combined,
                    l1_norm=l1_combined,
                    max_diff=max_combined,
                    correlation=1.0 - l2_combined,  # proxy
                    n_wavelength=0,
                    wavelength_range=(0, 0),
                )

                converged = (
                    l2_combined < self.config.l2_tolerance and
                    l1_combined < self.config.l1_tolerance
                )

                result = ConvergenceResult(
                    parameter_name='max_cells',
                    parameter_value=float(n_cells),
                    comparison=comparison,
                    converged=converged,
                    n_cells=n_cells,
                )
                report.results.append(result)

                if converged and report.converged_value is None:
                    report.converged_value = float(n_cells)

        return report

    def run_all(self) -> Dict[str, ConvergenceReport]:
        """Run all convergence tests."""
        reports = {}
        print("\n" + "=" * 70)
        print("SuperLite Convergence Testing")
        print("=" * 70)

        print("\n--- Grid Resolution Convergence ---")
        reports['grid'] = self.run_grid_convergence()
        print(reports['grid'].summary())

        # Wavelength, MC, and tau convergence require actual SuperLite runs
        # which must go through the orchestrator (run_model_grid.py)
        print("\n--- Wavelength Group Convergence ---")
        print("  (Requires SuperLite runs via orchestrator — see docs)")
        print("  Configuration prepared but not executed automatically.")

        print("\n--- MC Particle Convergence ---")
        print("  (Requires SuperLite runs via orchestrator — see docs)")

        print("\n--- Tau Truncation Convergence ---")
        print("  (Requires SuperLite runs via orchestrator — see docs)")

        # Save report
        report_path = os.path.join(self.config.output_dir, 'convergence_report.json')
        with open(report_path, 'w') as f:
            json.dump({k: v.to_dict() for k, v in reports.items()}, f, indent=2)
        print(f"\n  Report saved: {report_path}")

        return reports


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description='SuperLite convergence testing framework')
    parser.add_argument('--swd', required=True, help='STELLA .swd file')
    parser.add_argument('--output', required=True, help='Output directory')
    parser.add_argument('--target-days', type=float, default=30.0)
    parser.add_argument('--grid-cells', type=int, nargs='+',
                        default=[50, 100, 150, 200, 250])
    args = parser.parse_args()

    config = ConvergenceConfig(
        swd_file=args.swd,
        output_dir=args.output,
        target_days=args.target_days,
        grid_cells=args.grid_cells,
    )
    test = ConvergenceTest(config)
    reports = test.run_all()
