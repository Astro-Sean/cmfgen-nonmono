#!/usr/bin/env python3
"""Tests for cmfgen_convergence.py — synthetic convergence histories.

Tests that the classifier correctly distinguishes:
1. Strictly decreasing convergence
2. Convergence with occasional increases
3. Damped oscillation
4. Persistent oscillation
5. Stagnation
6. Clear divergence
7. Numerical noise around a very small correction
8. A negligible-abundance species with huge fractional changes
9. A genuinely important species with large absolute changes
"""

import sys
import math
import os
from pathlib import Path

# Add scripts dir to path
sys.path.insert(0, str(Path(__file__).parent))

from cmfgen_convergence import (
    ConvergenceHistory, IterationRecord, ConvergenceAssessment,
    assess_convergence, parse_outgen, format_assessment,
)


def make_history(corrections, lum_inner=1e9, lum_outer=1e8,
                 max_inc_depths=None, crashed=False, completed=False):
    """Build a ConvergenceHistory from a list of maxch values."""
    its = []
    for i, c in enumerate(corrections, 1):
        d = max_inc_depths[i - 1] if max_inc_depths else 1
        its.append(IterationRecord(
            iteration=i,
            max_inc=c,
            max_inc_depth=d,
            max_dec=min(c * 0.9, 110.0),
            max_dec_depth=d,
            maxch=c,
            lum_inner=lum_inner,
            lum_outer=lum_outer,
        ))
    return ConvergenceHistory(iterations=its, crashed=crashed, completed=completed)


# ---------------------------------------------------------------------------
# Test cases
# ---------------------------------------------------------------------------

def test_strictly_decreasing():
    """1. Strictly decreasing convergence → CONVERGING or CONVERGED."""
    corr = [1e20, 1e15, 1e10, 1e5, 1e2, 50, 20, 10, 5, 2]
    h = make_history(corr)
    a = assess_convergence(h, window=10)
    assert a.classification in ("CONVERGING", "CONVERGED"), \
        f"Expected CONVERGING/CONVERGED, got {a.classification}: {a.message}"
    assert a.trend == "decreasing", f"Expected decreasing trend, got {a.trend}"
    print(f"  test_strictly_decreasing: {a.classification} (slope={a.slope_log10:.3f}) ✓")


def test_convergence_with_occasional_increases():
    """2. Convergence with occasional increases → CONVERGING.

    The key test: occasional increases must NOT trigger divergence.
    """
    corr = [0.45, 0.21, 0.16, 0.19, 0.11, 0.08, 0.09, 0.04]
    # Scale up to % (CMFGEN reports in %)
    corr = [c * 100 for c in corr]
    h = make_history(corr)
    a = assess_convergence(h, window=8)
    assert a.classification in ("CONVERGING", "CONVERGED", "DAMPED_OSCILLATION"), \
        f"Expected CONVERGING/CONVERGED/DAMPED_OSCILLATION, got {a.classification}: {a.message}"
    # Must NOT be classified as diverging
    assert a.classification != "DIVERGING", \
        "Occasional increases incorrectly classified as DIVERGING"
    print(f"  test_convergence_with_occasional_increases: {a.classification} (slope={a.slope_log10:.3f}) ✓")


def test_damped_oscillation():
    """3. Damped oscillation → DAMPED_OSCILLATION or CONVERGING."""
    # Oscillating but amplitude decreasing
    corr = [30, 18, 22, 14, 17, 9, 12, 7, 8, 4]
    h = make_history(corr)
    a = assess_convergence(h, window=10)
    assert a.classification in ("DAMPED_OSCILLATION", "CONVERGING", "CONVERGED"), \
        f"Expected DAMPED_OSCILLATION/CONVERGING, got {a.classification}: {a.message}"
    print(f"  test_damped_oscillation: {a.classification} (osc={a.oscillation_indicator}, slope={a.slope_log10:.3f}) ✓")


def test_persistent_oscillation():
    """4. Persistent oscillation → PERSISTENT_OSCILLATION."""
    # Oscillating without amplitude decrease, above tolerance
    corr = [10, 18, 9, 20, 11, 19, 10, 20, 11, 19]
    h = make_history(corr)
    a = assess_convergence(h, window=10)
    assert a.classification == "PERSISTENT_OSCILLATION", \
        f"Expected PERSISTENT_OSCILLATION, got {a.classification}: {a.message}"
    print(f"  test_persistent_oscillation: {a.classification} (osc={a.oscillation_indicator}) ✓")


def test_stagnation():
    """5. Stagnation → STAGNATING."""
    # Plateau above tolerance, no trend
    corr = [12, 11, 12, 11, 12, 11, 12, 11, 12, 11]
    h = make_history(corr)
    a = assess_convergence(h, window=10)
    assert a.classification == "STAGNATING", \
        f"Expected STAGNATING, got {a.classification}: {a.message}"
    print(f"  test_stagnation: {a.classification} (stag={a.stagnation_indicator}) ✓")


def test_clear_divergence():
    """6. Clear divergence → DIVERGING."""
    corr = [10, 16, 25, 38, 57, 85, 120, 180, 270, 400]
    h = make_history(corr)
    a = assess_convergence(h, window=10)
    assert a.classification == "DIVERGING", \
        f"Expected DIVERGING, got {a.classification}: {a.message}"
    assert a.trend == "increasing", f"Expected increasing trend, got {a.trend}"
    print(f"  test_clear_divergence: {a.classification} (slope={a.slope_log10:.3f}) ✓")


def test_numerical_noise_small_correction():
    """7. Numerical noise around a very small correction → CONVERGED."""
    # Corrections below tolerance, with tiny noise
    corr = [3.1, 2.9, 3.0, 2.8, 3.1, 2.9, 3.0, 2.8, 3.1, 2.9]
    h = make_history(corr)
    a = assess_convergence(h, window=10, converged_tol=10.0)
    assert a.classification == "CONVERGED", \
        f"Expected CONVERGED, got {a.classification}: {a.message}"
    print(f"  test_numerical_noise_small_correction: {a.classification} (max={a.max_correction:.2f}%) ✓")


def test_negligible_species_huge_fractional():
    """8. Negligible-abundance species with huge fractional changes.

    In CMFGEN, the 'Maximm changes' is the maximum unscaled STEQ residual.
    A trace species can have enormous fractional corrections while the
    major species are well-converged. The major_scale metric helps
    distinguish this.

    Here we simulate: maxch is huge (trace species) but major_scale
    is small (major species converging). The classification should
    not be DIVERGING if the major species are converging.

    However, since we only use maxch for the trend, this test verifies
    that the classifier doesn't crash and produces a reasonable output.
    The real mitigation is in the assessment message noting that
    maxch may be dominated by trace species.
    """
    corr = [1e15, 1e12, 1e10, 1e8, 1e6, 1e4, 1e2, 50, 20, 10]
    h = make_history(corr)
    # Add major_scale info showing major species are converging
    for i, it in enumerate(h.iterations):
        it.major_scale = 1e-6 * (1 + i * 0.1)  # increasing (less damping)
    a = assess_convergence(h, window=10)
    # The trend is decreasing, so should be CONVERGING
    assert a.classification in ("CONVERGING", "CONVERGED"), \
        f"Expected CONVERGING/CONVERGED, got {a.classification}: {a.message}"
    print(f"  test_negligible_species_huge_fractional: {a.classification} "
          f"(slope={a.slope_log10:.3f}) ✓")


def test_important_species_large_absolute():
    """9. Genuinely important species with large absolute changes → DIVERGING."""
    # Large corrections that are increasing — a real problem
    corr = [1e5, 1e8, 1e12, 1e15, 1e20, 1e25, 1e30, 1e35, 1e40, 1e45]
    h = make_history(corr)
    a = assess_convergence(h, window=10)
    assert a.classification == "DIVERGING", \
        f"Expected DIVERGING, got {a.classification}: {a.message}"
    print(f"  test_important_species_large_absolute: {a.classification} "
          f"(slope={a.slope_log10:.3f}) ✓")


def test_negative_luminosity_warning():
    """Test that negative outer luminosity is flagged."""
    corr = [1e5, 1e4, 1e3, 1e2, 50, 20, 10, 5, 3, 2]
    h = make_history(corr, lum_inner=1e9, lum_outer=-5e8)
    a = assess_convergence(h, window=10)
    assert a.negative_luminosity, "Expected negative_luminosity flag"
    assert "negative" in a.message.lower(), "Expected negative luminosity warning in message"
    print(f"  test_negative_luminosity_warning: negative_lum={a.negative_luminosity} ✓")


def test_crashed():
    """Test that crash is detected."""
    h = make_history([1e5, 1e4], crashed=True)
    a = assess_convergence(h, window=10)
    assert a.classification == "CRASHED", f"Expected CRASHED, got {a.classification}"
    print(f"  test_crashed: {a.classification} ✓")


def test_completed_converged():
    """Test that completed + small corrections = CONVERGED."""
    corr = [1e5, 1e4, 1e3, 1e2, 50, 20, 10, 5, 3, 2]
    h = make_history(corr, completed=True)
    a = assess_convergence(h, window=10, converged_tol=10.0)
    assert a.classification == "CONVERGED", f"Expected CONVERGED, got {a.classification}"
    print(f"  test_completed_converged: {a.classification} (max={a.max_correction:.2f}%) ✓")


def test_completed_not_converged():
    """Test that completed + large corrections = STAGNATING."""
    corr = [1e5, 1e4, 1e3, 500, 400, 350, 300, 280, 250, 230]
    h = make_history(corr, completed=True)
    a = assess_convergence(h, window=10, converged_tol=10.0)
    # Completed but corrections still high — should be STAGNATING or CONVERGING
    # depending on whether the trend is flat or decreasing
    assert a.classification in ("STAGNATING", "CONVERGING"), \
        f"Expected STAGNATING/CONVERGING, got {a.classification}: {a.message}"
    print(f"  test_completed_not_converged: {a.classification} "
          f"(max={a.max_correction:.2e}%, slope={a.slope_log10:.3f}) ✓")


def test_too_few_iterations():
    """Test that < 3 iterations = STARTING."""
    h = make_history([1e5, 1e4])
    a = assess_convergence(h, window=10)
    assert a.classification == "STARTING", f"Expected STARTING, got {a.classification}"
    print(f"  test_too_few_iterations: {a.classification} ✓")


def test_outgen_parsing():
    """Test parsing a real OUTGEN file if available."""
    # Try to find an OUTGEN file
    candidates = [
        Path('/viper/ptmp1/sbrennan/SN2024pba/X/cmfgen/35.57days/OUTGEN'),
        Path('/ptmp/sbrennan/SN2024pba/X/cmfgen/35.57days/OUTGEN'),
    ]
    outgen = None
    for c in candidates:
        if c.exists():
            outgen = c
            break
    if outgen is None:
        print("  test_outgen_parsing: SKIPPED (no OUTGEN found)")
        return

    history = parse_outgen(outgen)
    assert history.n_iterations > 0, "Parsed 0 iterations from OUTGEN"
    print(f"  test_outgen_parsing: parsed {history.n_iterations} iterations ✓")

    # Assess convergence
    a = assess_convergence(history, window=10)
    print(f"  test_outgen_parsing: classification={a.classification}, "
          f"max={a.max_correction:.2e}%, slope={a.slope_log10:.3f}")


def test_rolling_window():
    """Test that the rolling window correctly uses only recent iterations."""
    # First 20 iterations: huge, then 10 iterations: small
    corr = [1e20] * 20 + [5, 4, 3, 2, 1, 5, 4, 3, 2, 1]
    h = make_history(corr)
    a = assess_convergence(h, window=10)
    # The window should only see the small corrections
    assert a.max_correction < 100, \
        f"Window should only see recent small corrections, got max={a.max_correction:.2e}"
    assert a.classification == "CONVERGED", \
        f"Expected CONVERGED, got {a.classification}: {a.message}"
    print(f"  test_rolling_window: max={a.max_correction:.2f}%, "
          f"window={a.window_size} ✓")


def run_all_tests():
    """Run all tests and report results."""
    print("=" * 60)
    print("CMFGEN Convergence Diagnostics Tests")
    print("=" * 60)
    print()

    tests = [
        test_strictly_decreasing,
        test_convergence_with_occasional_increases,
        test_damped_oscillation,
        test_persistent_oscillation,
        test_stagnation,
        test_clear_divergence,
        test_numerical_noise_small_correction,
        test_negligible_species_huge_fractional,
        test_important_species_large_absolute,
        test_negative_luminosity_warning,
        test_crashed,
        test_completed_converged,
        test_completed_not_converged,
        test_too_few_iterations,
        test_rolling_window,
        test_outgen_parsing,
    ]

    passed = 0
    failed = 0
    for test in tests:
        try:
            test()
            passed += 1
        except AssertionError as e:
            print(f"  {test.__name__}: FAILED - {e}")
            failed += 1
        except Exception as e:
            print(f"  {test.__name__}: ERROR - {e}")
            failed += 1

    print()
    print("=" * 60)
    print(f"Results: {passed} passed, {failed} failed, {len(tests)} total")
    print("=" * 60)
    return failed == 0


if __name__ == '__main__':
    success = run_all_tests()
    sys.exit(0 if success else 1)
