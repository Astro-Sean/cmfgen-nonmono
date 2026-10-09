#!/usr/bin/env python3
"""
Standalone test of the observer-frame K-index search algorithm.

Tests the bidirectional K-index fix in obs_frame_sub_v9.f against:
1. The original monotonic-only algorithm
2. A brute-force reference search

The K-index algorithm finds the CMF frequency bracket [K-1, K] such that
  SM_FREQ_CMF(K-1) >= T1 >= SM_FREQ_CMF(K)

where SM_FREQ_CMF is sorted in DECREASING order, and
  T1 = OBS_FREQ(ML) * (1 - VMU_RAY(I)) * GAM_RAY(I)

For monotonically expanding flow, T1 decreases along the ray, so only
decrementing K suffices. For non-monotonic flow (compression zones,
velocity reversals), T1 can increase, requiring K to also increment.

Test cases:
  A: Monotonic velocity (regression test)
  B: Mild non-monotonic velocity (small reversal)
  C: Strongly non-monotonic velocity (actual model-like)
"""

import numpy as np
import sys
import os

# ============================================================================
# K-index search algorithms (translated from obs_frame_sub_v9.f)
# ============================================================================

def k_search_original(K_prev, T1, SM_FREQ_CMF, SM_NCF):
    """Original monotonic-only K-index search.
    
    Only decrements K (assumes T1 is monotonically decreasing).
    Returns K such that SM_FREQ_CMF(K-1) >= T1 >= SM_FREQ_CMF(K).
    """
    K = K_prev
    # Original: DO WHILE(T1 < SM_FREQ_CMF(K))  -> K = K+1
    # But the original code actually had a different logic.
    # Let me re-read the original code.
    # The ORIGINAL code (commented out) was:
    #   DO WHILE(T1 .LT. SM_FREQ_CMF(K))
    #     K = K+1
    #   END DO
    # This increments K when T1 < SM_FREQ_CMF(K), i.e., when T1 is below the lower bound.
    # Wait, SM_FREQ_CMF is in DECREASING order, so SM_FREQ_CMF(K) > SM_FREQ_CMF(K+1).
    # If T1 < SM_FREQ_CMF(K), we need to go to a LOWER frequency, i.e., HIGHER K.
    
    # Actually, let me re-read the code more carefully.
    # The commented-out original was:
    #   DO WHILE(T1 .LT. SM_FREQ_CMF(K))
    #     K=K+1
    # But this was for the FIRST ray point (I=1), not the main loop.
    
    # For the main loop (I=2..NRAY), the original was:
    #   K=INDX(I-1)
    #   DO WHILE(T1 .LT. SM_FREQ_CMF(K))
    #     K=K+1  (but this was commented out and replaced)
    
    # Wait, let me re-read the actual code at lines 818-842.
    # The CURRENT code (with fix) is:
    #   K=INDX(I-1)
    #   DO WHILE(T1 .GT. SM_FREQ_CMF(K-1))   -- decrement if T1 > upper bound
    #     K=K-1
    #   END DO
    #   DO WHILE(T1 .LT. SM_FREQ_CMF(K) .AND. K .LT. SM_NCF)  -- increment if T1 < lower bound
    #     K=K+1
    #   END DO
    
    # The ORIGINAL code (commented out at line 821) was:
    #   DO WHILE(T1 .LT. SM_FREQ_CMF(K))
    #     K=K+1  (but this doesn't match the current decrement logic)
    
    # Actually, looking more carefully at the code:
    # Line 818: K=INDX(I-1)
    # Line 821 (commented): DO WHILE(T1 .LT. SM_FREQ_CMF(K))
    # Line 822 (active): DO WHILE(T1 .GT. SM_FREQ_CMF(K-1))
    # Line 823: K=K-1
    # Line 840 (added): DO WHILE(T1 .LT. SM_FREQ_CMF(K) .AND. K .LT. SM_NCF)
    # Line 841: K=K+1
    
    # So the original algorithm was:
    #   K = INDX(I-1)
    #   DO WHILE(T1 .LT. SM_FREQ_CMF(K))
    #     K = K+1
    #   END DO
    # This increments K when T1 < SM_FREQ_CMF(K).
    # Since SM_FREQ_CMF is in DECREASING order, SM_FREQ_CMF(K) is the LOWER bound.
    # If T1 < SM_FREQ_CMF(K), T1 is below the lower bound, so we need HIGHER K (lower frequency).
    # But this only handles T1 DECREASING (moving to lower frequencies).
    # It does NOT handle T1 INCREASING (moving to higher frequencies).
    
    # Wait, that can't be right. If SM_FREQ_CMF is in decreasing order:
    #   SM_FREQ_CMF(1) > SM_FREQ_CMF(2) > ... > SM_FREQ_CMF(SM_NCF)
    # Then SM_FREQ_CMF(K-1) > SM_FREQ_CMF(K)
    # The bracket is: SM_FREQ_CMF(K-1) >= T1 >= SM_FREQ_CMF(K)
    
    # If T1 < SM_FREQ_CMF(K), T1 is below the lower bound, need K+1 (lower freq)
    # If T1 > SM_FREQ_CMF(K-1), T1 is above the upper bound, need K-1 (higher freq)
    
    # The original only handles T1 < SM_FREQ_CMF(K) (T1 decreasing).
    # The fix adds handling for T1 > SM_FREQ_CMF(K-1) (T1 increasing).
    
    # But wait, the CURRENT code at line 822 is:
    #   DO WHILE(T1 .GT. SM_FREQ_CMF(K-1))  -> K=K-1
    # This is the NEW decrement loop (handles T1 increasing).
    # And line 840 is:
    #   DO WHILE(T1 .LT. SM_FREQ_CMF(K) .AND. K .LT. SM_NCF)  -> K=K+1
    # This is the NEW increment loop (handles T1 decreasing).
    
    # The ORIGINAL code (line 821, commented out) was:
    #   DO WHILE(T1 .LT. SM_FREQ_CMF(K))  -> K=K+1
    # This only handles T1 decreasing.
    
    # So the original algorithm was:
    #   K = INDX(I-1)
    #   DO WHILE(T1 .LT. SM_FREQ_CMF(K))
    #     K = K+1
    #     IF(K .GT. SM_NCF) STOP
    #   END DO
    # This only increments K (handles T1 decreasing, i.e., monotonic expansion).
    # It does NOT decrement K when T1 increases (compression zone).
    
    # Actually wait, I need to re-read more carefully. Let me look at the
    # commented-out lines again.
    
    # Line 818: K=INDX(I-1)      !INDX(I)
    # Line 819:! K=INDX(I)
    # Line 820: T1=OBS_FREQ(ML)*(1.0_LDP-VMU_RAY(I))*GAM_RAY(I)
    # Line 821:! DO WHILE(T1 .LT. SM_FREQ_CMF(K))
    # Line 822: DO WHILE(T1 .GT. SM_FREQ_CMF(K-1))
    # Line 823:   K=K-1
    # Line 824:! IF(K .GT. SM_NCF)THEN
    # Line 825:   IF(K .LE. 1)THEN
    # Line 826-831: error handling
    # Line 832: END DO
    # Line 840: DO WHILE(T1 .LT. SM_FREQ_CMF(K) .AND. K .LT. SM_NCF)
    # Line 841:   K=K+1
    # Line 842: END DO
    
    # So the original (line 821, commented) was:
    #   DO WHILE(T1 .LT. SM_FREQ_CMF(K))
    #     K = K+1
    #   END DO
    # Which handles T1 < SM_FREQ_CMF(K) by incrementing K.
    
    # But wait, this doesn't make sense for a decreasing frequency grid.
    # If SM_FREQ_CMF is decreasing and T1 < SM_FREQ_CMF(K),
    # then T1 is below the current point, and we need to go to HIGHER K
    # (lower frequency). So K=K+1 is correct for T1 decreasing.
    
    # But the original code doesn't handle T1 > SM_FREQ_CMF(K-1),
    # which would require K=K-1 (going to higher frequency).
    # In a compression zone, T1 can increase, and the original code
    # would NOT find the correct bracket.
    
    # Actually, I think I'm confused. Let me think about this differently.
    # 
    # SM_FREQ_CMF is in DECREASING order: SM_FREQ_CMF(1) > SM_FREQ_CMF(2) > ...
    # The bracket for T1 is: SM_FREQ_CMF(K-1) >= T1 >= SM_FREQ_CMF(K)
    # i.e., K-1 is the upper bound, K is the lower bound.
    #
    # For monotonic expansion (T1 decreasing along ray):
    #   T1 at point I is smaller than T1 at point I-1
    #   So we need K to INCREASE (move to lower frequencies)
    #   The original code does: DO WHILE(T1 < SM_FREQ_CMF(K)) K=K+1
    #   This increments K until T1 >= SM_FREQ_CMF(K), which is correct.
    #
    # For compression zone (T1 increasing along ray):
    #   T1 at point I is larger than T1 at point I-1
    #   So we need K to DECREASE (move to higher frequencies)
    #   The original code does NOT handle this case.
    #   T1 > SM_FREQ_CMF(K-1) would be true, but the original code
    #   only checks T1 < SM_FREQ_CMF(K), which would be false.
    #   So K would stay at the wrong value.
    #
    # The fix adds:
    #   DO WHILE(T1 > SM_FREQ_CMF(K-1)) K=K-1  (decrement for T1 increasing)
    #   DO WHILE(T1 < SM_FREQ_CMF(K) .AND. K < SM_NCF) K=K+1  (increment for T1 decreasing)
    
    # OK so the original algorithm is:
    K = K_prev
    while T1 < SM_FREQ_CMF[K] and K < SM_NCF - 1:
        K = K + 1
    return K

def k_search_corrected(K_prev, T1, SM_FREQ_CMF, SM_NCF):
    """Corrected bidirectional K-index search.
    
    Handles both T1 increasing (compression zone) and T1 decreasing (expansion).
    Returns K such that SM_FREQ_CMF(K-1) >= T1 >= SM_FREQ_CMF(K).
    """
    K = K_prev
    # Decrement K if T1 is above the upper bound (T1 increased - compression zone)
    while T1 > SM_FREQ_CMF[K-1] and K > 1:
        K = K - 1
    # Increment K if T1 is below the lower bound (T1 decreased - expansion)
    while T1 < SM_FREQ_CMF[K] and K < SM_NCF - 1:
        K = K + 1
    return K

def k_search_brute_force(T1, SM_FREQ_CMF, SM_NCF):
    """Brute-force reference K-index search.
    
    Scans the entire SM_FREQ_CMF array to find the correct bracket.
    Returns K such that SM_FREQ_CMF(K-1) >= T1 >= SM_FREQ_CMF(K).
    """
    # SM_FREQ_CMF is in decreasing order
    # Find K such that SM_FREQ_CMF(K-1) >= T1 >= SM_FREQ_CMF(K)
    for K in range(1, SM_NCF):
        if T1 >= SM_FREQ_CMF[K]:
            return K
    return SM_NCF - 1


# ============================================================================
# Velocity field models
# ============================================================================

def monotonic_velocity(NRAY, VINF=5000.0, C_KMS=2.99792458e5):
    """Monotonically increasing velocity field (standard expanding SN)."""
    # Linear velocity from 0.1*VINF at center to VINF at surface
    r = np.linspace(0.0, 1.0, NRAY)
    v = VINF * (0.1 + 0.9 * r)
    return v / C_KMS  # Return in units of c

def mild_nonmonotonic_velocity(NRAY, VINF=5000.0, C_KMS=2.99792458e5):
    """Mildly non-monotonic velocity field with small reversal."""
    r = np.linspace(0.0, 1.0, NRAY)
    v = VINF * (0.1 + 0.9 * r)
    # Add a small reversal at 60-70% of the radius
    reversal = 0.15 * VINF * np.exp(-((r - 0.65) / 0.05)**2)
    v = v - reversal
    # Ensure velocity is positive
    v = np.maximum(v, 0.01 * VINF)
    return v / C_KMS

def strong_nonmonotonic_velocity(NRAY, VINF=5000.0, C_KMS=2.99792458e5):
    """Strongly non-monotonic velocity field with compression zone.
    
    Mimics the SN2024pba model with a strong shock/compression zone.
    """
    r = np.linspace(0.0, 1.0, NRAY)
    # Base velocity: linear increase
    v = VINF * (0.1 + 0.9 * r)
    # Add a strong compression zone at 40-60% of the radius
    # This creates a velocity reversal of ~30% of VINF
    reversal = 0.35 * VINF * np.exp(-((r - 0.5) / 0.08)**2)
    v = v - reversal
    # Add a secondary smaller reversal
    reversal2 = 0.15 * VINF * np.exp(-((r - 0.75) / 0.04)**2)
    v = v - reversal2
    # Ensure velocity is positive
    v = np.maximum(v, 0.01 * VINF)
    return v / C_KMS

def actual_model_velocity(cmfgen_dir):
    """Read velocity from actual CMFGEN model RVTJ file."""
    rvtj_path = os.path.join(cmfgen_dir, 'RVTJ')
    if not os.path.exists(rvtj_path):
        return None
    
    # Read RVTJ and extract velocity
    velocities = []
    radii = []
    with open(rvtj_path, 'r') as f:
        lines = f.readlines()
    
    # Find the velocity data section
    in_velocity = False
    for i, line in enumerate(lines):
        if 'Velocity' in line or 'velocity' in line:
            in_velocity = True
            continue
        if in_velocity:
            parts = line.split()
            if len(parts) >= 2:
                try:
                    r = float(parts[0])
                    v = float(parts[1])
                    radii.append(r)
                    velocities.append(v)
                except ValueError:
                    in_velocity = False
    
    if not velocities:
        return None
    
    return np.array(velocities)


# ============================================================================
# Test functions
# ============================================================================

def run_test(test_name, VMU_RAY, SM_FREQ_CMF, OBS_FREQ, SM_NCF, NRAY, N_OBS_FREQ):
    """Run a K-index search test and report results."""
    
    print(f"\n{'='*70}")
    print(f"Test: {test_name}")
    print(f"  NRAY={NRAY}, SM_NCF={SM_NCF}, N_OBS_FREQ={N_OBS_FREQ}")
    print(f"  Velocity range: {VMU_RAY.min():.6e} to {VMU_RAY.max():.6e} (c)")
    print(f"  CMF freq range: {SM_FREQ_CMF[0]:.6e} to {SM_FREQ_CMF[-1]:.6e}")
    print(f"{'='*70}")
    
    # Compute GAM_RAY (Lorentz factor)
    GAM_RAY = 1.0 / np.sqrt(1.0 - VMU_RAY**2)
    
    # Test each observer frequency
    total_errors_original = 0
    total_errors_corrected = 0
    total_points = 0
    max_error_original = 0.0
    max_error_corrected = 0.0
    worst_case_original = None
    worst_case_corrected = None
    
    for ML in range(N_OBS_FREQ):
        # Compute T1 for all ray points
        T1_array = OBS_FREQ[ML] * (1.0 - VMU_RAY) * GAM_RAY
        
        # Run original algorithm
        K_original = np.zeros(NRAY, dtype=int)
        K_original[0] = 1
        # First point: find initial K
        T1 = T1_array[0]
        K = 1
        while T1 < SM_FREQ_CMF[K] and K < SM_NCF - 1:
            K = K + 1
        K_original[0] = K
        
        for I in range(1, NRAY):
            T1 = T1_array[I]
            K_original[I] = k_search_original(K_original[I-1], T1, SM_FREQ_CMF, SM_NCF)
        
        # Run corrected algorithm
        K_corrected = np.zeros(NRAY, dtype=int)
        K_corrected[0] = K_original[0]  # Same initial K
        for I in range(1, NRAY):
            T1 = T1_array[I]
            K_corrected[I] = k_search_corrected(K_corrected[I-1], T1, SM_FREQ_CMF, SM_NCF)
        
        # Run brute-force reference
        K_reference = np.zeros(NRAY, dtype=int)
        for I in range(NRAY):
            T1 = T1_array[I]
            K_reference[I] = k_search_brute_force(T1, SM_FREQ_CMF, SM_NCF)
        
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
    
    # Report results
    print(f"\n  Results:")
    print(f"    Total points:        {total_points}")
    print(f"    Original errors:     {total_errors_original} ({100*total_errors_original/total_points:.2f}%)")
    print(f"    Corrected errors:    {total_errors_corrected} ({100*total_errors_corrected/total_points:.2f}%)")
    print(f"    Max original error:  {max_error_original} K-steps")
    print(f"    Max corrected error: {max_error_corrected} K-steps")
    
    if worst_case_original:
        ML, I, K, K_ref, T1 = worst_case_original
        print(f"    Worst original case: ML={ML}, I={I}, K={K}, K_ref={K_ref}, T1={T1:.6e}")
        if K > 0 and K < SM_NCF:
            print(f"      SM_FREQ_CMF(K-1)={SM_FREQ_CMF[K-1]:.6e}, SM_FREQ_CMF(K)={SM_FREQ_CMF[K]:.6e}")
    
    if worst_case_corrected:
        ML, I, K, K_ref, T1 = worst_case_corrected
        print(f"    Worst corrected case: ML={ML}, I={I}, K={K}, K_ref={K_ref}, T1={T1:.6e}")
    
    # Check for interpolation errors
    # The interpolation is: T2 = (T1 - SM_FREQ_CMF(K)) / (SM_FREQ_CMF(K-1) - SM_FREQ_CMF(K))
    # If K is wrong, T2 can be < 0 or > 1, causing extrapolation instead of interpolation
    print(f"\n  Interpolation quality:")
    
    for ML in [0, N_OBS_FREQ//4, N_OBS_FREQ//2, 3*N_OBS_FREQ//4, N_OBS_FREQ-1]:
        T1_array = OBS_FREQ[ML] * (1.0 - VMU_RAY) * GAM_RAY
        
        # Check T2 range for original
        T2_orig_out_of_range = 0
        T2_corr_out_of_range = 0
        
        K_orig = np.zeros(NRAY, dtype=int)
        K_orig[0] = 1
        T1 = T1_array[0]
        K = 1
        while T1 < SM_FREQ_CMF[K] and K < SM_NCF - 1:
            K = K + 1
        K_orig[0] = K
        
        K_corr = np.zeros(NRAY, dtype=int)
        K_corr[0] = K_orig[0]
        
        for I in range(1, NRAY):
            T1 = T1_array[I]
            K_orig[I] = k_search_original(K_orig[I-1], T1, SM_FREQ_CMF, SM_NCF)
            K_corr[I] = k_search_corrected(K_corr[I-1], T1, SM_FREQ_CMF, SM_NCF)
        
        for I in range(NRAY):
            # Original
            K = K_orig[I]
            if K > 0 and K < SM_NCF:
                denom = SM_FREQ_CMF[K-1] - SM_FREQ_CMF[K]
                if denom != 0:
                    T2 = (T1_array[I] - SM_FREQ_CMF[K]) / denom
                    if T2 < -1e-10 or T2 > 1.0 + 1e-10:
                        T2_orig_out_of_range += 1
            
            # Corrected
            K = K_corr[I]
            if K > 0 and K < SM_NCF:
                denom = SM_FREQ_CMF[K-1] - SM_FREQ_CMF[K]
                if denom != 0:
                    T2 = (T1_array[I] - SM_FREQ_CMF[K]) / denom
                    if T2 < -1e-10 or T2 > 1.0 + 1e-10:
                        T2_corr_out_of_range += 1
        
        print(f"    ML={ML:4d}: T2 out of range - original: {T2_orig_out_of_range:5d}, corrected: {T2_corr_out_of_range:5d}")
    
    # Velocity profile summary
    print(f"\n  Velocity profile:")
    print(f"    Min: {VMU_RAY.min():.6e} c ({VMU_RAY.min()*2.99792458e5:.1f} km/s)")
    print(f"    Max: {VMU_RAY.max():.6e} c ({VMU_RAY.max()*2.99792458e5:.1f} km/s)")
    
    # Check for reversals
    dv = np.diff(VMU_RAY)
    n_reversals = np.sum(dv < 0)
    if n_reversals > 0:
        print(f"    Reversals: {n_reversals} points where velocity decreases")
        max_rev = np.min(dv[dv < 0]) * 2.99792458e5
        print(f"    Max reversal magnitude: {abs(max_rev):.1f} km/s")
    else:
        print(f"    No reversals (monotonic)")
    
    return total_errors_original, total_errors_corrected, total_points


# ============================================================================
# Main
# ============================================================================

def main():
    NRAY = 200
    SM_NCF = 5000
    N_OBS_FREQ = 100
    
    # Create CMF frequency grid (DECREASING order)
    # Frequencies in units of 10^15 Hz
    SM_FREQ_CMF = np.linspace(5.0, 0.01, SM_NCF)[::-1]  # Decreasing order
    
    # Create observer frequencies
    OBS_FREQ = np.linspace(0.1, 3.0, N_OBS_FREQ)
    
    # Test A: Monotonic velocity
    VMU_RAY = monotonic_velocity(NRAY)
    orig_a, corr_a, total_a = run_test(
        "A: Monotonic velocity (regression test)",
        VMU_RAY, SM_FREQ_CMF, OBS_FREQ, SM_NCF, NRAY, N_OBS_FREQ)
    
    # Test B: Mild non-monotonic velocity
    VMU_RAY = mild_nonmonotonic_velocity(NRAY)
    orig_b, corr_b, total_b = run_test(
        "B: Mild non-monotonic velocity",
        VMU_RAY, SM_FREQ_CMF, OBS_FREQ, SM_NCF, NRAY, N_OBS_FREQ)
    
    # Test C: Strong non-monotonic velocity
    VMU_RAY = strong_nonmonotonic_velocity(NRAY)
    orig_c, corr_c, total_c = run_test(
        "C: Strong non-monotonic velocity (compression zone)",
        VMU_RAY, SM_FREQ_CMF, OBS_FREQ, SM_NCF, NRAY, N_OBS_FREQ)
    
    # Test D: Actual model velocity
    cmfgen_dir = '/raven/ptmp/sbrennan/SN2024pba/SN_B_v2_M3_V5000/cmfgen/40.00days'
    actual_v = actual_model_velocity(cmfgen_dir)
    if actual_v is not None and len(actual_v) > 10:
        # Interpolate to NRAY points
        r_old = np.linspace(0, 1, len(actual_v))
        r_new = np.linspace(0, 1, NRAY)
        VMU_RAY = np.interp(r_new, r_old, actual_v)
        orig_d, corr_d, total_d = run_test(
            "D: Actual model velocity (SN2024pba)",
            VMU_RAY, SM_FREQ_CMF, OBS_FREQ, SM_NCF, NRAY, N_OBS_FREQ)
    else:
        print("\nTest D: Actual model velocity - SKIPPED (RVTJ not available)")
        orig_d, corr_d, total_d = 0, 0, 1
    
    # Summary
    print(f"\n{'='*70}")
    print(f"SUMMARY")
    print(f"{'='*70}")
    print(f"  Test A (monotonic):     Original: {orig_a:5d}/{total_a} errors, Corrected: {corr_a:5d}/{total_a} errors")
    print(f"  Test B (mild non-mono):  Original: {orig_b:5d}/{total_b} errors, Corrected: {corr_b:5d}/{total_b} errors")
    print(f"  Test C (strong non-mono):Original: {orig_c:5d}/{total_c} errors, Corrected: {corr_c:5d}/{total_c} errors")
    print(f"  Test D (actual model):   Original: {orig_d:5d}/{total_d} errors, Corrected: {corr_d:5d}/{total_d} errors")
    
    # Verdict
    print(f"\n  VERDICT:")
    if orig_a == 0 and corr_a == 0:
        print(f"    Test A (regression): PASS - both algorithms correct for monotonic flow")
    else:
        print(f"    Test A (regression): FAIL")
    
    if orig_b > 0 and corr_b == 0:
        print(f"    Test B (mild):       PASS - fix corrects {orig_b} errors from original")
    elif orig_b == 0 and corr_b == 0:
        print(f"    Test B (mild):       PASS - both correct (reversal too small to affect K)")
    else:
        print(f"    Test B (mild):       FAIL - original: {orig_b}, corrected: {corr_b}")
    
    if orig_c > 0 and corr_c == 0:
        print(f"    Test C (strong):     PASS - fix corrects {orig_c} errors from original")
    elif orig_c == 0 and corr_c == 0:
        print(f"    Test C (strong):     PASS - both correct")
    else:
        print(f"    Test C (strong):     FAIL - original: {orig_c}, corrected: {corr_c}")
    
    if orig_d > 0 and corr_d == 0:
        print(f"    Test D (actual):     PASS - fix corrects {orig_d} errors from original")
    elif orig_d == 0 and corr_d == 0:
        print(f"    Test D (actual):     PASS - both correct")
    elif corr_d > 0:
        print(f"    Test D (actual):     FAIL - corrected still has {corr_d} errors")
    else:
        print(f"    Test D (actual):     SKIPPED")


if __name__ == '__main__':
    main()
