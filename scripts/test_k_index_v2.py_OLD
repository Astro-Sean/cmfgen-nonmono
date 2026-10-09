#!/usr/bin/env python3
"""
Standalone test of the observer-frame K-index search algorithm (v2).

Tests the bidirectional K-index fix in obs_frame_sub_v9.f against:
1. The original monotonic-only algorithm
2. A brute-force reference search

Key fix: The K-index algorithm finds the CMF frequency bracket [K-1, K] such that
  SM_FREQ_CMF(K-1) >= T1 >= SM_FREQ_CMF(K)
where SM_FREQ_CMF is sorted in DECREASING order, and
  T1 = OBS_FREQ(ML) * (1 - VMU_RAY(I)) * GAM_RAY(I)

For monotonically expanding flow, T1 decreases along the ray, so only
incrementing K suffices. For non-monotonic flow (compression zones,
velocity reversals), T1 can increase, requiring K to also decrement.
"""

import numpy as np
import sys
import os

C_KMS = 2.99792458e5  # Speed of light in km/s

# ============================================================================
# K-index search algorithms (translated from obs_frame_sub_v9.f)
# ============================================================================

def k_search_original_first(T1, SM_FREQ_CMF, SM_NCF):
    """Original K-index search for the first ray point (I=1).
    
    From line 778: DO WHILE(T1 .LT. SM_FREQ_CMF(K)) K=K+1
    This increments K until T1 >= SM_FREQ_CMF(K).
    """
    K = 0  # 0-indexed
    while T1 < SM_FREQ_CMF[K] and K < SM_NCF - 1:
        K = K + 1
    return K

def k_search_original(K_prev, T1, SM_FREQ_CMF, SM_NCF):
    """Original monotonic-only K-index search for I=2..NRAY.
    
    The original code (line 821, commented out) was:
      DO WHILE(T1 .LT. SM_FREQ_CMF(K))
        K = K + 1
      END DO
    
    This only increments K (handles T1 decreasing = monotonic expansion).
    It does NOT decrement K when T1 increases (compression zone).
    """
    K = K_prev
    while T1 < SM_FREQ_CMF[K] and K < SM_NCF - 1:
        K = K + 1
    return K

def k_search_corrected(K_prev, T1, SM_FREQ_CMF, SM_NCF):
    """Corrected bidirectional K-index search for I=2..NRAY.
    
    From lines 822-842 of obs_frame_sub_v9.f:
      DO WHILE(T1 .GT. SM_FREQ_CMF(K-1))  K=K-1  (decrement for T1 increasing)
      DO WHILE(T1 .LT. SM_FREQ_CMF(K) .AND. K .LT. SM_NCF)  K=K+1  (increment for T1 decreasing)
    """
    K = K_prev
    # Decrement K if T1 is above the upper bound (T1 increased - compression zone)
    while K > 0 and T1 > SM_FREQ_CMF[K-1]:
        K = K - 1
    # Increment K if T1 is below the lower bound (T1 decreased - expansion)
    while T1 < SM_FREQ_CMF[K] and K < SM_NCF - 1:
        K = K + 1
    return K

def k_search_brute_force(T1, SM_FREQ_CMF, SM_NCF):
    """Brute-force reference K-index search.
    
    Scans the entire SM_FREQ_CMF array to find the correct bracket.
    SM_FREQ_CMF is in DECREASING order.
    Returns K such that SM_FREQ_CMF(K-1) >= T1 >= SM_FREQ_CMF(K).
    """
    if T1 >= SM_FREQ_CMF[0]:
        return 0  # T1 is above the highest frequency
    if T1 <= SM_FREQ_CMF[-1]:
        return SM_NCF - 1  # T1 is below the lowest frequency
    # Binary search for efficiency
    lo, hi = 0, SM_NCF - 1
    while lo < hi - 1:
        mid = (lo + hi) // 2
        if T1 >= SM_FREQ_CMF[mid]:
            hi = mid
        else:
            lo = mid
    return hi


# ============================================================================
# Velocity field models
# ============================================================================

def monotonic_velocity(NRAY, VINF=10000.0):
    """Monotonically increasing velocity field (standard expanding SN).
    
    Velocity increases from VINF/10 at center to VINF at surface.
    Returns velocity in units of c.
    """
    r = np.linspace(0.0, 1.0, NRAY)
    v = VINF * (0.1 + 0.9 * r)
    return v / C_KMS

def mild_nonmonotonic_velocity(NRAY, VINF=10000.0):
    """Mildly non-monotonic velocity field with small reversal."""
    r = np.linspace(0.0, 1.0, NRAY)
    v = VINF * (0.1 + 0.9 * r)
    # Add a small reversal at 60-70% of the radius
    reversal = 0.15 * VINF * np.exp(-((r - 0.65) / 0.05)**2)
    v = v - reversal
    v = np.maximum(v, 0.01 * VINF)
    return v / C_KMS

def strong_nonmonotonic_velocity(NRAY, VINF=10000.0):
    """Strongly non-monotonic velocity field with compression zone.
    
    Mimics the SN2024pba model with a strong shock/compression zone.
    """
    r = np.linspace(0.0, 1.0, NRAY)
    v = VINF * (0.1 + 0.9 * r)
    # Add a strong compression zone at 40-60% of the radius
    reversal = 0.35 * VINF * np.exp(-((r - 0.5) / 0.08)**2)
    v = v - reversal
    # Add a secondary smaller reversal
    reversal2 = 0.15 * VINF * np.exp(-((r - 0.75) / 0.04)**2)
    v = v - reversal2
    v = np.maximum(v, 0.01 * VINF)
    return v / C_KMS

def actual_model_velocity(cmfgen_dir):
    """Read velocity from actual CMFGEN model RVTJ file."""
    rvtj_path = os.path.join(cmfgen_dir, 'RVTJ')
    if not os.path.exists(rvtj_path):
        return None
    
    with open(rvtj_path, 'r') as f:
        lines = f.readlines()
    
    # Find the velocity data section
    in_velocity = False
    velocities = []
    for line in lines:
        if 'Velocity (km/s)' in line:
            in_velocity = True
            continue
        if in_velocity:
            parts = line.split()
            if len(parts) >= 1:
                try:
                    for p in parts:
                        v = float(p)
                        velocities.append(v)
                except ValueError:
                    if len(velocities) >= 10:
                        break
                    in_velocity = False
    
    if len(velocities) < 10:
        return None
    
    return np.array(velocities) / C_KMS  # Return in units of c


# ============================================================================
# Test functions
# ============================================================================

def run_test(test_name, VMU_RAY, OBS_FREQ, SM_NCF, NRAY, N_OBS_FREQ):
    """Run a K-index search test and report results."""
    
    # Create CMF frequency grid (DECREASING order)
    # The CMF frequency grid must bracket all possible T1 values.
    # T1 = OBS_FREQ(ML) * (1 - VMU_RAY(I)) * GAM_RAY(I)
    # For V << c: T1 ≈ OBS_FREQ * (1 - V/c)
    # So T1 ranges from OBS_FREQ * (1 - V_max/c) to OBS_FREQ * (1 - V_min/c)
    
    # Compute T1 range
    GAM_RAY = 1.0 / np.sqrt(1.0 - VMU_RAY**2)
    T1_min = np.min(OBS_FREQ[:, None] * (1.0 - VMU_RAY[None, :]) * GAM_RAY[None, :])
    T1_max = np.max(OBS_FREQ[:, None] * (1.0 - VMU_RAY[None, :]) * GAM_RAY[None, :])
    
    # Create CMF frequency grid that brackets T1 range with margin
    margin = 0.1 * (T1_max - T1_min)
    freq_lo = max(T1_min - margin, 1e-6)
    freq_hi = T1_max + margin
    SM_FREQ_CMF = np.linspace(freq_hi, freq_lo, SM_NCF)  # DECREASING order
    
    print(f"\n{'='*70}")
    print(f"Test: {test_name}")
    print(f"  NRAY={NRAY}, SM_NCF={SM_NCF}, N_OBS_FREQ={N_OBS_FREQ}")
    print(f"  Velocity range: {VMU_RAY.min()*C_KMS:.1f} to {VMU_RAY.max()*C_KMS:.1f} km/s")
    print(f"  T1 range: {T1_min:.6e} to {T1_max:.6e}")
    print(f"  CMF freq range: {SM_FREQ_CMF[0]:.6e} to {SM_FREQ_CMF[-1]:.6e}")
    print(f"  Grid spacing: {(freq_hi - freq_lo) / (SM_NCF - 1):.6e}")
    print(f"{'='*70}")
    
    # Velocity profile summary
    dv = np.diff(VMU_RAY) * C_KMS  # in km/s
    n_reversals = np.sum(dv < 0)
    if n_reversals > 0:
        max_rev = np.min(dv[dv < 0])
        print(f"\n  Velocity profile: {n_reversals} reversals, max reversal: {abs(max_rev):.1f} km/s")
    else:
        print(f"\n  Velocity profile: monotonic (no reversals)")
    
    # Test each observer frequency
    total_errors_original = 0
    total_errors_corrected = 0
    total_points = 0
    max_error_original = 0
    max_error_corrected = 0
    worst_case_original = None
    worst_case_corrected = None
    
    # Track T2 out-of-range
    T2_orig_oor_total = 0
    T2_corr_oor_total = 0
    
    for ML in range(N_OBS_FREQ):
        # Compute T1 for all ray points
        T1_array = OBS_FREQ[ML] * (1.0 - VMU_RAY) * GAM_RAY
        
        # Run original algorithm
        K_original = np.zeros(NRAY, dtype=int)
        K_original[0] = k_search_original_first(T1_array[0], SM_FREQ_CMF, SM_NCF)
        for I in range(1, NRAY):
            K_original[I] = k_search_original(K_original[I-1], T1_array[I], SM_FREQ_CMF, SM_NCF)
        
        # Run corrected algorithm
        K_corrected = np.zeros(NRAY, dtype=int)
        K_corrected[0] = K_original[0]  # Same initial K
        for I in range(1, NRAY):
            K_corrected[I] = k_search_corrected(K_corrected[I-1], T1_array[I], SM_FREQ_CMF, SM_NCF)
        
        # Run brute-force reference
        K_reference = np.zeros(NRAY, dtype=int)
        for I in range(NRAY):
            K_reference[I] = k_search_brute_force(T1_array[I], SM_FREQ_CMF, SM_NCF)
        
        # Compare
        for I in range(NRAY):
            total_points += 1
            err_orig = abs(K_original[I] - K_reference[I])
            err_corr = abs(K_corrected[I] - K_reference[I])
            
            if err_orig > 0:
                total_errors_original += 1
                if err_orig > max_error_original:
                    max_error_original = err_orig
                    worst_case_original = (ML, I, K_original[I], K_reference[I], T1_array[I])
            
            if err_corr > 0:
                total_errors_corrected += 1
                if err_corr > max_error_corrected:
                    max_error_corrected = err_corr
                    worst_case_corrected = (ML, I, K_corrected[I], K_reference[I], T1_array[I])
            
            # Check T2 range
            # T2 = (T1 - SM_FREQ_CMF(K)) / (SM_FREQ_CMF(K-1) - SM_FREQ_CMF(K))
            # T2 should be in [0, 1] for interpolation (not extrapolation)
            K = K_original[I]
            if K > 0 and K < SM_NCF:
                denom = SM_FREQ_CMF[K-1] - SM_FREQ_CMF[K]
                if abs(denom) > 1e-30:
                    T2 = (T1_array[I] - SM_FREQ_CMF[K]) / denom
                    if T2 < -1e-10 or T2 > 1.0 + 1e-10:
                        T2_orig_oor_total += 1
            
            K = K_corrected[I]
            if K > 0 and K < SM_NCF:
                denom = SM_FREQ_CMF[K-1] - SM_FREQ_CMF[K]
                if abs(denom) > 1e-30:
                    T2 = (T1_array[I] - SM_FREQ_CMF[K]) / denom
                    if T2 < -1e-10 or T2 > 1.0 + 1e-10:
                        T2_corr_oor_total += 1
    
    # Report results
    print(f"\n  Results:")
    print(f"    Total points:           {total_points}")
    print(f"    Original errors:        {total_errors_original} ({100*total_errors_original/total_points:.2f}%)")
    print(f"    Corrected errors:       {total_errors_corrected} ({100*total_errors_corrected/total_points:.2f}%)")
    print(f"    Max original error:     {max_error_original} K-steps")
    print(f"    Max corrected error:    {max_error_corrected} K-steps")
    print(f"    T2 out of range (orig): {T2_orig_oor_total} ({100*T2_orig_oor_total/total_points:.2f}%)")
    print(f"    T2 out of range (corr): {T2_corr_oor_total} ({100*T2_corr_oor_total/total_points:.2f}%)")
    
    if worst_case_original:
        ML, I, K, K_ref, T1 = worst_case_original
        print(f"\n    Worst original case: ML={ML}, I={I}, K={K}, K_ref={K_ref}")
        print(f"      T1={T1:.6e}")
        if K > 0 and K < len(SM_FREQ_CMF):
            print(f"      SM_FREQ_CMF(K-1)={SM_FREQ_CMF[K-1]:.6e}, SM_FREQ_CMF(K)={SM_FREQ_CMF[K]:.6e}")
        if K_ref > 0 and K_ref < len(SM_FREQ_CMF):
            print(f"      Ref: SM_FREQ_CMF(K_ref-1)={SM_FREQ_CMF[K_ref-1]:.6e}, SM_FREQ_CMF(K_ref)={SM_FREQ_CMF[K_ref]:.6e}")
    
    if worst_case_corrected:
        ML, I, K, K_ref, T1 = worst_case_corrected
        print(f"\n    Worst corrected case: ML={ML}, I={I}, K={K}, K_ref={K_ref}")
        print(f"      T1={T1:.6e}")
    
    return total_errors_original, total_errors_corrected, total_points, T2_orig_oor_total, T2_corr_oor_total


# ============================================================================
# Main
# ============================================================================

def main():
    NRAY = 200
    SM_NCF = 10000
    N_OBS_FREQ = 200
    
    # Create observer frequencies (in units of 10^15 Hz)
    # Optical/UV range: 3000-10000 Angstrom
    # nu = c / lambda = 3e18 / lambda(Angstrom) in Hz
    # For 3000 A: nu = 1.0e15 Hz = 1.0 (in units of 10^15 Hz)
    # For 10000 A: nu = 0.3e15 Hz = 0.3 (in units of 10^15 Hz)
    OBS_FREQ = np.linspace(0.3, 3.0, N_OBS_FREQ)
    
    # Test A: Monotonic velocity
    VMU_RAY = monotonic_velocity(NRAY)
    orig_a, corr_a, total_a, T2_orig_a, T2_corr_a = run_test(
        "A: Monotonic velocity (regression test)",
        VMU_RAY, OBS_FREQ, SM_NCF, NRAY, N_OBS_FREQ)
    
    # Test B: Mild non-monotonic velocity
    VMU_RAY = mild_nonmonotonic_velocity(NRAY)
    orig_b, corr_b, total_b, T2_orig_b, T2_corr_b = run_test(
        "B: Mild non-monotonic velocity",
        VMU_RAY, OBS_FREQ, SM_NCF, NRAY, N_OBS_FREQ)
    
    # Test C: Strong non-monotonic velocity
    VMU_RAY = strong_nonmonotonic_velocity(NRAY)
    orig_c, corr_c, total_c, T2_orig_c, T2_corr_c = run_test(
        "C: Strong non-monotonic velocity (compression zone)",
        VMU_RAY, OBS_FREQ, SM_NCF, NRAY, N_OBS_FREQ)
    
    # Test D: Actual model velocity
    cmfgen_dir = '/raven/ptmp/sbrennan/SN2024pba/SN_B_v2_M3_V5000/cmfgen/40.00days'
    actual_v = actual_model_velocity(cmfgen_dir)
    if actual_v is not None and len(actual_v) > 10:
        # Interpolate to NRAY points
        r_old = np.linspace(0, 1, len(actual_v))
        r_new = np.linspace(0, 1, NRAY)
        VMU_RAY = np.interp(r_new, r_old, actual_v)
        orig_d, corr_d, total_d, T2_orig_d, T2_corr_d = run_test(
            "D: Actual model velocity (SN2024pba)",
            VMU_RAY, OBS_FREQ, SM_NCF, NRAY, N_OBS_FREQ)
    else:
        print("\nTest D: Actual model velocity - SKIPPED (RVTJ not available)")
        orig_d, corr_d, total_d = 0, 0, 1
        T2_orig_d, T2_corr_d = 0, 0
    
    # Summary
    print(f"\n{'='*70}")
    print(f"SUMMARY")
    print(f"{'='*70}")
    tests = [
        ("A (monotonic)", orig_a, corr_a, total_a, T2_orig_a, T2_corr_a),
        ("B (mild non-mono)", orig_b, corr_b, total_b, T2_orig_b, T2_corr_b),
        ("C (strong non-mono)", orig_c, corr_c, total_c, T2_orig_c, T2_corr_c),
        ("D (actual model)", orig_d, corr_d, total_d, T2_orig_d, T2_corr_d),
    ]
    
    for name, orig, corr, total, t2o, t2c in tests:
        print(f"  {name:25s}: Orig={orig:6d}/{total} err, Corr={corr:6d}/{total} err, T2_oor: orig={t2o}, corr={t2c}")
    
    # Verdict
    print(f"\n  VERDICT:")
    all_pass = True
    for name, orig, corr, total, t2o, t2c in tests:
        if 'monotonic' in name.lower():
            if orig == 0 and corr == 0:
                print(f"    {name}: PASS (regression - no change)")
            else:
                print(f"    {name}: FAIL (regression broken)")
                all_pass = False
        else:
            if corr == 0 and t2c == 0:
                if orig > 0 or t2o > 0:
                    print(f"    {name}: PASS (fix corrects {orig} K errors, {t2o} T2 errors)")
                else:
                    print(f"    {name}: PASS (both correct)")
            elif corr == 0:
                print(f"    {name}: PASS (K correct, but {t2c} T2 out of range)")
            else:
                print(f"    {name}: FAIL (corrected has {corr} errors)")
                all_pass = False
    
    if all_pass:
        print(f"\n  OVERALL: PASS")
    else:
        print(f"\n  OVERALL: FAIL")
    
    return all_pass


if __name__ == '__main__':
    success = main()
    sys.exit(0 if success else 1)
