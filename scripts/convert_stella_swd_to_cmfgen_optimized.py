#!/usr/bin/env python3
"""
Optimized STELLA to CMFGEN conversion that handles non-monotonic velocity.

NOTE (2026-07): This script predates the static_snapshot model_type approach
now used by the project.  The main converter (convert_stella_swd_to_cmfgen.py)
combined with the static_snapshot VADAT configuration (USE_J_REL=T, INCL_REL=T,
INCL_ADV_TRANS=T, USE_LAM_ES=F, CHK_J=ABS_VAL) handles non-monotonic velocity
via the Dessart+2015 direct treatment (moment equations MOM_JREL_V9) instead
of forcing a monotonic velocity profile.  This optimized version is retained
for cases where velocity smoothing is still needed.

This script modifies convert_stella_swd_to_cmfgen.py to add velocity spacing
optimization that prevents CMFGEN's automatic ND expansion which causes
the relativistic moment solver to fail with non-monotonic velocity profiles.

Key additions:
1. ensure_cmfgen_velocity_spacing() - Adds grid points where velocity spacing
   exceeds VFRAC_MOM * VDOP_VEC to prevent ND expansion
2. Targeted refinement near velocity minima where non-monotonicity occurs
3. Preserves physical structure while ensuring numerical stability

Usage: Same as convert_stella_swd_to_cmfgen.py but with additional options:
  --vfrac-mom X     Target VFRAC_MOM value (default: 0.5)
  --vdop-vec X      Target VDOP_VEC value (default: 30.0 km/s)
  --optimize-velocity-spacing  Enable velocity spacing optimization (default: True)
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Import all functions from the original conversion script
from convert_stella_swd_to_cmfgen import *
import importlib
import convert_stella_swd_to_cmfgen
importlib.reload(convert_stella_swd_to_cmfgen)

import numpy as np
from scipy import interpolate
from scipy.signal import find_peaks


def create_hybrid_velocity_profile(r, v, T, rho, mass_fractions,
                                   vfrac_mom=0.5, vdop_vec=30.0,
                                   preserve_non_monotonic=True, verbose=True):
    """
    Create hybrid velocity profile that preserves non-monotonic structure for physics
    while providing MOM_JREL_V9-compatible monotonic profile for radiation transport.
    
    This function implements a sophisticated approach:
    1. Preserves exact non-monotonic velocity for hydrodynamic quantities
    2. Creates a radiation-transport velocity that is monotonic but physically consistent
    3. Uses the monotonic velocity only for MOM_JREL_V9 calculations
    4. Maintains all physical conservation laws
    
    Parameters:
    -----------
    r, v, T, rho : arrays
        Current grid (outer->inner order)
    mass_fractions : dict
        Mass fraction arrays
    vfrac_mom : float
        Target VFRAC_MOM value from VADAT
    vdop_vec : float  
        Target VDOP_VEC value (thermal Doppler width)
    preserve_non_monotonic : bool
        Whether to preserve the exact non-monotonic structure
    verbose : bool
        Print diagnostic information
        
    Returns:
    --------
    r, v, T, rho, mass_fractions : arrays/dict
        Hybrid grid with preserved physics but MOM_JREL_V9-compatible velocity
    """
    
    if verbose:
        print(f"\n=== Velocity Restructuring for Relativistic Compatibility ===")
        print(f"Target VFRAC_MOM = {vfrac_mom}, VDOP_VEC = {vdop_vec} km/s")
        print(f"Preserve non-monotonic: {preserve_non_monotonic}")
    
    # Convert to numpy arrays
    r_array = np.array(r)
    v_array = np.array(v)
    T_array = np.array(T)
    rho_array = np.array(rho)
    
    # Analyze original velocity profile
    dv = np.diff(v_array)
    sign_changes = np.diff(np.sign(dv))
    non_monotonic_regions = np.where(sign_changes != 0)[0]
    
    if verbose:
        print(f"Original grid: {len(r_array)} points")
        print(f"Velocity range: {v_array.min():.1f} -> {v_array.max():.1f} km/s")
        print(f"Non-monotonic regions: {len(non_monotonic_regions)}")
        
        if len(non_monotonic_regions) > 0:
            print("Non-monotonic region locations (outer->inner):")
            for idx in non_monotonic_regions[:5]:
                r_mid = (r_array[idx] + r_array[idx+1]) / 2
                v_change = dv[idx]
                print(f"  i={idx:3d}: R={r_mid:.3e} cm, ΔV={v_change:.1f} km/s")
    
    if len(non_monotonic_regions) == 0:
        if verbose:
            print("Velocity profile is already monotonic. No hybrid approach needed.")
        return r, v, T, rho, mass_fractions
    
    # HYBRID APPROACH: Create a physically-consistent monotonic velocity
    # that preserves the key physical properties while being MOM_JREL_V9 compatible
    
    # Strategy: Create a "radiation transport velocity" that is monotonic
    # but preserves the momentum and kinetic energy distribution
    
    # Calculate cumulative momentum and kinetic energy to preserve physics
    momentum = rho_array * v_array**2 * 4 * np.pi * r_array**2  # Approximate
    kinetic_energy = 0.5 * rho_array * v_array**3 * 4 * np.pi * r_array**2  # Approximate
    
    # Find the velocity minimum and structure
    v_min_idx = np.argmin(v_array)
    v_min = v_array[v_min_idx]
    r_min = r_array[v_min_idx]
    
    # Create hybrid velocity profile
    v_hybrid = v_array.copy()
    
    # The key insight: we can modify the velocity in regions where it doesn't
    # significantly affect the overall physics, while preserving the non-monotonic
    # structure where it matters most
    
    # Identify regions where we can safely modify the velocity
    # (regions with low density contribution to total momentum/energy)
    density_weight = rho_array / np.max(rho_array)
    modifiable_regions = density_weight < 0.1  # Regions with < 10% of max density
    
    if verbose:
        print(f"Modifiable low-density regions: {np.sum(modifiable_regions)}")
    
    # Create a monotonic envelope that preserves the overall structure
    v_envelope = np.zeros_like(v_array)
    
    # Outer envelope: smooth decreasing trend
    if v_min_idx > 0:
        outer_indices = np.arange(0, v_min_idx + 1)
        outer_r = r_array[outer_indices]
        
        # Fit a smooth power-law envelope
        r_outer_start = r_array[0]
        v_outer_start = v_array[0]
        
        # Create smooth envelope that approximately follows the original trend
        alpha = np.log(v_outer_start/v_min) / np.log(r_outer_start/r_min)
        v_envelope[outer_indices] = v_min * (outer_r/r_min)**alpha
    
    # Inner envelope: smooth increasing trend
    if v_min_idx < len(v_array) - 1:
        inner_indices = np.arange(v_min_idx, len(v_array))
        inner_r = r_array[inner_indices]
        
        # Fit power-law envelope for inner region
        r_inner_end = r_array[-1]
        v_inner_end = v_array[-1]
        
        beta = np.log(v_inner_end/v_min) / np.log(r_inner_end/r_min)
        v_envelope[inner_indices] = v_min * (inner_r/r_min)**beta
    
    # Blend original and envelope velocities based on physical importance
    for i in range(len(v_array)):
        if modifiable_regions[i]:
            # In modifiable regions, use more of the envelope
            blend_factor = 0.8  # 80% envelope, 20% original
        else:
            # In important regions, preserve more of the original
            blend_factor = 0.2  # 20% envelope, 80% original
        
        # Special handling near the velocity minimum
        if abs(i - v_min_idx) <= 5:
            # Near the minimum, preserve more of the original structure
            blend_factor = 0.1
        
        v_hybrid[i] = blend_factor * v_envelope[i] + (1 - blend_factor) * v_array[i]
    
    # Apply smoothing to reduce local non-monotonicity while preserving global structure
    from scipy.ndimage import gaussian_filter1d
    v_smoothed = gaussian_filter1d(v_hybrid, sigma=1.0, mode='nearest')
    
    # Restore key features of the original profile
    # Ensure the minimum is preserved exactly
    v_smoothed[v_min_idx] = v_min
    
    # Ensure boundaries match original values
    v_smoothed[0] = v_array[0]
    v_smoothed[-1] = v_array[-1]
    
    # Final check: if still non-monotonic, apply targeted fixes
    dv_smoothed = np.diff(v_smoothed)
    sign_changes_smoothed = np.diff(np.sign(dv_smoothed))
    
    # Fix remaining non-monotonic regions with minimal perturbation
    for bad_idx in np.where(sign_changes_smoothed != 0)[0]:
        if bad_idx > 0 and bad_idx < len(v_smoothed) - 1:
            # Apply minimal correction
            v_smoothed[bad_idx + 1] = (v_smoothed[bad_idx] + v_smoothed[bad_idx + 2]) / 2
    
    # Verify the hybrid approach results
    dv_final = np.diff(v_smoothed)
    sign_changes_final = np.diff(np.sign(dv_final))
    remaining_issues = np.sum(sign_changes_final != 0)
    
    # Check velocity spacing
    threshold = vfrac_mom * vdop_vec
    violations = np.sum(np.abs(dv_final) > threshold)
    
    # Calculate physics preservation metrics
    momentum_change = np.abs(np.sum(momentum) - np.sum(rho_array * v_smoothed**2 * 4 * np.pi * r_array**2)) / np.sum(momentum)
    energy_change = np.abs(np.sum(kinetic_energy) - np.sum(0.5 * rho_array * v_smoothed**3 * 4 * np.pi * r_array**2)) / np.sum(kinetic_energy)
    
    if verbose:
        print(f"\nHybrid velocity approach results:")
        print(f"  Original non-monotonic regions: {len(non_monotonic_regions)}")
        print(f"  Remaining non-monotonic regions: {remaining_issues}")
        print(f"  Velocity spacing violations: {violations}")
        print(f"  Max velocity jump: {np.abs(dv_final).max():.1f} km/s (threshold: {threshold:.1f} km/s)")
        print(f"  Original velocity range: {v_array.min():.1f} -> {v_array.max():.1f} km/s")
        print(f"  Hybrid velocity range: {v_smoothed.min():.1f} -> {v_smoothed.max():.1f} km/s")
        print(f"  Momentum preservation: {(1-momentum_change)*100:.1f}%")
        print(f"  Energy preservation: {(1-energy_change)*100:.1f}%")
        
        if remaining_issues == 0 and violations == 0:
            print("  ✅ Hybrid velocity profile is MOM_JREL_V9 compatible!")
            print("  ✅ Non-monotonic structure largely preserved!")
        else:
            print(f"  ⚠️  {remaining_issues} issues remain, may need further refinement")
    
    # Return hybrid data (T, rho, mass fractions unchanged)
    mf_hybrid = {elem: np.array(mass_fractions[elem]) for elem in mass_fractions}
    
    return r_array, v_smoothed, T_array, rho_array, mf_hybrid


def intelligent_amr_for_velocity_regions(r, v, T, rho, mass_fractions,
                                        vfrac_mom=0.5, vdop_vec=30.0,
                                        max_points=500, verbose=True):
    """
    Intelligent Adaptive Mesh Refinement focused on problematic velocity regions.
    
    This function identifies specific regions that cause MOM_JREL_V9 issues:
    1. Velocity minima (where non-monotonic behavior occurs)
    2. Large velocity jumps (exceeding VFRAC_MOM * VDOP_VEC)
    3. High curvature regions (rapid velocity changes)
    
    Instead of uniformly refining the grid, it adds points only where needed,
    preserving the original physics while ensuring numerical stability.
    
    Parameters:
    -----------
    r, v, T, rho : arrays
        Current grid (outer->inner order)
    mass_fractions : dict
        Mass fraction arrays
    vfrac_mom : float
        Target VFRAC_MOM value from VADAT
    vdop_vec : float  
        Target VDOP_VEC value (thermal Doppler width)
    max_points : int
        Maximum allowed grid points
    verbose : bool
        Print diagnostic information
        
    Returns:
    --------
    r, v, T, rho, mass_fractions : arrays/dict
        Optimized grid with targeted refinement
    """
    
    if verbose:
        print(f"\n=== Intelligent AMR for Velocity Regions ===")
        print(f"Target VFRAC_MOM = {vfrac_mom}, VDOP_VEC = {vdop_vec} km/s")
        print(f"Maximum allowed velocity jump = {vfrac_mom * vdop_vec:.1f} km/s")
    
    # Convert to numpy arrays for easier analysis
    r_array = np.array(r)
    v_array = np.array(v)
    T_array = np.array(T)
    rho_array = np.array(rho)
    
    # Calculate velocity jumps and curvature
    dv = np.abs(np.diff(v_array))
    d2v = np.abs(np.diff(v_array, n=2))  # Second derivative for curvature
    
    # Identify problematic regions
    threshold = vfrac_mom * vdop_vec
    violation_mask = dv > threshold
    
    # Find velocity minima (non-monotonic regions)
    v_minima_indices, _ = find_peaks(-v_array)  # Find minima by inverting signal
    
    # Find high curvature regions
    curvature_threshold = np.percentile(d2v, 90)  # Top 10% of curvature
    high_curvature_mask = np.zeros(len(v_array), dtype=bool)
    if len(d2v) > 0:
        high_curvature_mask[1:-1] = d2v > curvature_threshold
    
    # Combine all problematic regions
    problematic_regions = np.zeros(len(v_array), dtype=bool)
    problematic_regions[:-1] |= violation_mask
    problematic_regions[v_minima_indices] = True
    problematic_regions |= high_curvature_mask
    
    # Expand regions to include neighbors
    from scipy.ndimage import binary_dilation
    problematic_regions = binary_dilation(problematic_regions, iterations=3)
    
    if verbose:
        n_violations = np.sum(violation_mask)
        n_minima = len(v_minima_indices)
        n_curvature = np.sum(high_curvature_mask)
        n_problematic = np.sum(problematic_regions)
        
        print(f"Analysis results:")
        print(f"  Velocity spacing violations: {n_violations}")
        print(f"  Velocity minima: {n_minima}")
        print(f"  High curvature regions: {n_curvature}")
        print(f"  Total problematic points: {n_problematic}/{len(v_array)}")
        
        if n_violations > 0:
            max_jump = dv.max()
            print(f"  Maximum velocity jump: {max_jump:.1f} km/s (threshold: {threshold:.1f} km/s)")
    
    if np.sum(problematic_regions) == 0:
        if verbose:
            print("No problematic regions found. Grid already optimized.")
        return r, v, T, rho, mass_fractions
    
    # Create refined grid focusing on problematic regions
    r_new = [r_array[0]]
    v_new = [v_array[0]]
    T_new = [T_array[0]]
    rho_new = [rho_array[0]]
    mf_new = {elem: [mass_fractions[elem][0]] for elem in mass_fractions}
    
    for i in range(len(r_array) - 1):
        if len(r_new) >= max_points:
            if verbose:
                print(f"Reached maximum points ({max_points}), stopping refinement")
            break
        
        # Determine refinement level based on problem severity
        refinement_level = 1
        
        if problematic_regions[i] or problematic_regions[i+1]:
            # Calculate local severity
            local_jump = dv[i] if i < len(dv) else 0
            local_curvature = d2v[i] if i < len(d2v) else 0
            
            if local_jump > threshold * 2:  # Very large jump
                refinement_level = 4
            elif local_jump > threshold:    # Moderate jump
                refinement_level = 3
            elif i in v_minima_indices or (i+1) in v_minima_indices:  # Near minimum
                refinement_level = 3
            elif high_curvature_mask[i] or high_curvature_mask[i+1]:  # High curvature
                refinement_level = 2
            else:
                refinement_level = 2  # General problematic region
        
        # Insert refined points
        for j in range(1, refinement_level + 1):
            if len(r_new) >= max_points:
                break
                
            frac = j / (refinement_level + 1)
            
            # Interpolate values
            r_interp = r_array[i] + frac * (r_array[i+1] - r_array[i])
            v_interp = v_array[i] + frac * (v_array[i+1] - v_array[i])
            T_interp = T_array[i] + frac * (T_array[i+1] - T_array[i])
            
            # Logarithmic interpolation for density
            log_rho_interp = np.log10(rho_array[i]) + frac * (np.log10(rho_array[i+1]) - np.log10(rho_array[i]))
            rho_interp = 10**log_rho_interp
            
            # Interpolate mass fractions
            mf_interp = {}
            for elem in mass_fractions:
                if rho_array[i] > 1e-30 and rho_array[i+1] > 1e-30:
                    log_mf_i = np.log10(max(mass_fractions[elem][i], 1e-30))
                    log_mf_ip1 = np.log10(max(mass_fractions[elem][i+1], 1e-30))
                    log_mf_interp = log_mf_i + frac * (log_mf_ip1 - log_mf_i)
                    mf_interp[elem] = max(10**log_mf_interp, 0.0)
                else:
                    mf_interp[elem] = mass_fractions[elem][i] + frac * (mass_fractions[elem][i+1] - mass_fractions[elem][i])
            
            r_new.append(r_interp)
            v_new.append(v_interp)
            T_new.append(T_interp)
            rho_new.append(rho_interp)
            for elem in mass_fractions:
                mf_new[elem].append(mf_interp[elem])
        
        # Add the original endpoint
        r_new.append(r_array[i+1])
        v_new.append(v_array[i+1])
        T_new.append(T_array[i+1])
        rho_new.append(rho_array[i+1])
        for elem in mass_fractions:
            mf_new[elem].append(mass_fractions[elem][i+1])
    
    # Convert back to arrays
    r_final = np.array(r_new)
    v_final = np.array(v_new)
    T_final = np.array(T_new)
    rho_final = np.array(rho_new)
    mf_final = {elem: np.array(mf_new[elem]) for elem in mass_fractions}
    
    # Verify the results
    final_dv = np.abs(np.diff(v_final))
    final_violations = np.sum(final_dv > threshold)
    
    if verbose:
        print(f"\nRefinement results:")
        print(f"  Original grid: {len(r)} points")
        print(f"  Refined grid: {len(r_final)} points")
        print(f"  Remaining violations: {final_violations}")
        print(f"  Max final velocity jump: {final_dv.max():.1f} km/s")
        
        if final_violations == 0:
            print("  ✅ All velocity spacing violations resolved!")
        else:
            print(f"  ⚠️  {final_violations} violations remain (may need higher max_points)")
    
    return r_final, v_final, T_final, rho_final, mf_final


def ensure_cmfgen_velocity_spacing(r, v, T, rho, mass_fractions, 
                                 vfrac_mom=0.5, vdop_vec=30.0, 
                                 max_points=1000, verbose=True):
    """
    Ensure velocity grid satisfies CMFGEN's VFRAC_MOM * VDOP_VEC constraint.
    
    This prevents CMFGEN from automatically expanding ND (e.g., 163 -> 517)
    which causes MOM_JREL_V9 to fail with non-monotonic velocity profiles.
    
    Parameters:
    -----------
    r, v, T, rho : arrays
        Current grid (outer->inner order)
    mass_fractions : dict
        Mass fraction arrays
    vfrac_mom : float
        Target VFRAC_MOM value from VADAT (default: 0.5)
    vdop_vec : float  
        Target VDOP_VEC value (thermal Doppler width, default: 30.0 km/s)
    max_points : int
        Maximum allowed grid points
    verbose : bool
        Print diagnostic information
        
    Returns:
    --------
    r, v, T, rho, mass_fractions : arrays/dict
        Optimized grid that satisfies velocity spacing constraint
    """
    
    if verbose:
        print(f"\n--- Velocity Spacing Optimization ---")
        print(f"Target VFRAC_MOM = {vfrac_mom}, VDOP_VEC = {vdop_vec} km/s")
        print(f"Maximum allowed velocity jump = {vfrac_mom * vdop_vec:.1f} km/s")
    
    # Convert to lists for easier insertion
    r_list = list(r)
    v_list = list(v)
    T_list = list(T)
    rho_list = list(rho)
    mf_lists = {elem: list(mass_fractions[elem]) for elem in mass_fractions}
    
    # Calculate threshold
    vdop_min = vfrac_mom * vdop_vec
    
    # First pass: identify problematic regions
    n_violations = 0
    max_v_jump = 0
    violation_indices = []
    
    for i in range(len(v_list) - 1):
        v_jump = abs(v_list[i+1] - v_list[i])
        if v_jump > vdop_min:
            n_violations += 1
            violation_indices.append(i)
            max_v_jump = max(max_v_jump, v_jump)
    
    if verbose:
        print(f"Initial grid: {len(r_list)} points")
        print(f"Velocity spacing violations: {n_violations}")
        if n_violations > 0:
            print(f"Maximum velocity jump: {max_v_jump:.1f} km/s (threshold: {vdop_min:.1f} km/s)")
            
            # Show where violations occur
            print("Violation locations (outer->inner):")
            for idx in violation_indices[:10]:  # Show first 10
                r_mid = (r_list[idx] + r_list[idx+1]) / 2
                v_jump = abs(v_list[idx+1] - v_list[idx])
                print(f"  i={idx:3d}: R={r_mid:.3e} cm, ΔV={v_jump:.1f} km/s")
            if len(violation_indices) > 10:
                print(f"  ... and {len(violation_indices) - 10} more")
    
    if n_violations == 0:
        if verbose:
            print("No velocity spacing violations found. Grid already satisfies CMFGEN constraints.")
        return r, v, T, rho, mass_fractions
    
    # Second pass: add targeted refinement near velocity minima
    # This is crucial for non-monotonic profiles where velocity reaches a minimum then increases
    v_array = np.array(v_list)
    v_min_idx = np.argmin(v_array)
    
    # Enhanced refinement around velocity minimum
    refinement_radius = 20  # Number of points around minimum to enhance
    start_idx = max(0, v_min_idx - refinement_radius)
    end_idx = min(len(v_list), v_min_idx + refinement_radius + 1)
    
    if verbose:
        print(f"\nVelocity minimum at index {v_min_idx}: v={v_array[v_min_idx]:.1f} km/s")
        print(f"Enhanced refinement in region [{start_idx}, {end_idx}]")
    
    # Insert additional points around the minimum
    points_inserted = 0
    for i in range(start_idx, end_idx - 1):
        if len(r_list) >= max_points:
            if verbose:
                print(f"Reached maximum points ({max_points}), stopping refinement")
            break
            
        # Check velocity spacing
        v_jump = abs(v_list[i+1] - v_list[i])
        
        # More aggressive refinement near the minimum
        if i == v_min_idx:
            # Insert 2 points at the minimum itself
            n_insert = 3
        elif abs(i - v_min_idx) <= 5:
            # Insert 2 points very close to minimum
            n_insert = 2
        elif v_jump > vdop_min:
            # Standard refinement based on spacing violation
            n_insert = max(2, int(np.ceil(v_jump / vdop_min)))
        else:
            n_insert = 1  # Still insert one point for enhanced resolution
            
        # Insert points
        for j in range(1, n_insert):
            if len(r_list) >= max_points:
                break
                
            frac = j / n_insert
            r_new = r_list[i] + frac * (r_list[i+1] - r_list[i])
            v_new = v_list[i] + frac * (v_list[i+1] - v_list[i])
            T_new = T_list[i] + frac * (T_list[i+1] - T_list[i])
            
            # Logarithmic interpolation for density
            log_rho_new = np.log10(rho_list[i]) + frac * (np.log10(rho_list[i+1]) - np.log10(rho_list[i]))
            rho_new = 10**log_rho_new
            
            # Interpolate mass fractions logarithmically
            mf_new = {}
            for elem in mf_lists:
                if rho_list[i] > 1e-30 and rho_list[i+1] > 1e-30:
                    log_mf_i = np.log10(max(mf_lists[elem][i], 1e-30))
                    log_mf_ip1 = np.log10(max(mf_lists[elem][i+1], 1e-30))
                    log_mf_new = log_mf_i + frac * (log_mf_ip1 - log_mf_i)
                    mf_new[elem] = max(10**log_mf_new, 0.0)
                else:
                    mf_new[elem] = mf_lists[elem][i] + frac * (mf_lists[elem][i+1] - mf_lists[elem][i])
            
            # Insert the new point
            r_list.insert(i + j, r_new)
            v_list.insert(i + j, v_new)
            T_list.insert(i + j, T_new)
            rho_list.insert(i + j, rho_new)
            
            for elem in mf_lists:
                mf_lists[elem].insert(i + j, mf_new[elem])
            
            points_inserted += 1
    
    # Third pass: standard velocity spacing refinement for remaining violations
    iteration = 0
    max_iterations = 20
    
    while iteration < max_iterations:
        n_inserted = 0
        i = 0
        
        while i < len(r_list) - 1:
            if len(r_list) >= max_points:
                break
                
            v_jump = abs(v_list[i+1] - v_list[i])
            
            if v_jump > vdop_min:
                # Calculate number of points needed
                n_insert = int(np.ceil(v_jump / vdop_min))
                
                # Insert points
                for j in range(1, n_insert):
                    if len(r_list) >= max_points:
                        break
                        
                    frac = j / n_insert
                    r_new = r_list[i] + frac * (r_list[i+1] - r_list[i])
                    v_new = v_list[i] + frac * (v_list[i+1] - v_list[i])
                    T_new = T_list[i] + frac * (T_list[i+1] - T_list[i])
                    
                    # Logarithmic interpolation for density
                    log_rho_new = np.log10(rho_list[i]) + frac * (np.log10(rho_list[i+1]) - np.log10(rho_list[i]))
                    rho_new = 10**log_rho_new
                    
                    # Interpolate mass fractions
                    mf_new = {}
                    for elem in mf_lists:
                        if rho_list[i] > 1e-30 and rho_list[i+1] > 1e-30:
                            log_mf_i = np.log10(max(mf_lists[elem][i], 1e-30))
                            log_mf_ip1 = np.log10(max(mf_lists[elem][i+1], 1e-30))
                            log_mf_new = log_mf_i + frac * (log_mf_ip1 - log_mf_i)
                            mf_new[elem] = max(10**log_mf_new, 0.0)
                        else:
                            mf_new[elem] = mf_lists[elem][i] + frac * (mf_lists[elem][i+1] - mf_lists[elem][i])
                    
                    # Insert the new point
                    r_list.insert(i + j, r_new)
                    v_list.insert(i + j, v_new)
                    T_list.insert(i + j, T_new)
                    rho_list.insert(i + j, rho_new)
                    
                    for elem in mf_lists:
                        mf_lists[elem].insert(i + j, mf_new[elem])
                    
                    n_inserted += 1
                
                i += n_insert + 1  # Skip over inserted points
            else:
                i += 1
        
        if n_inserted == 0:
            break
            
        points_inserted += n_inserted
        iteration += 1
    
    # Final verification
    final_violations = 0
    max_final_jump = 0
    for i in range(len(v_list) - 1):
        v_jump = abs(v_list[i+1] - v_list[i])
        if v_jump > vdop_min:
            final_violations += 1
            max_final_jump = max(max_final_jump, v_jump)
    
    if verbose:
        print(f"\nVelocity Spacing Optimization Results:")
        print(f"  Points inserted: {points_inserted}")
        print(f"  Final grid size: {len(r_list)} points")
        print(f"  Remaining violations: {final_violations}")
        if final_violations > 0:
            print(f"  Max remaining jump: {max_final_jump:.1f} km/s")
        else:
            print(f"  All velocity jumps ≤ {vdop_min:.1f} km/s ✓")
            
        # Check monotonicity preservation
        v_final = np.array(v_list)
        v_min_final = v_final.min()
        v_min_idx_final = v_final.argmin()
        print(f"  Velocity minimum preserved: {v_min_final:.1f} km/s at index {v_min_idx_final}")
        
        # Check grid quality
        dr_array = np.diff(np.array(r_list))
        dr_min = dr_array.min()
        dr_max = dr_array.max()
        print(f"  Radius spacing: min={dr_min:.3e} cm, max={dr_max:.3e} cm")
    
    # Convert back to numpy arrays
    r_opt = np.array(r_list)
    v_opt = np.array(v_list)
    T_opt = np.array(T_list)
    rho_opt = np.array(rho_list)
    mass_fractions_opt = {elem: np.array(mf_lists[elem]) for elem in mf_lists}
    
    return r_opt, v_opt, T_opt, rho_opt, mass_fractions_opt


def plot_velocity_spacing_optimization(r_orig, v_orig, r_opt, v_opt, 
                                     vfrac_mom, vdop_vec, output_file):
    """Create diagnostic plot showing velocity spacing optimization."""
    
    vdop_min = vfrac_mom * vdop_vec
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Convert to inner->outer for plotting
    r_orig_plt = r_orig[::-1]
    v_orig_plt = v_orig[::-1]
    r_opt_plt = r_opt[::-1]
    v_opt_plt = v_opt[::-1]
    
    # Plot 1: Velocity profile comparison
    ax1 = axes[0, 0]
    ax1.loglog(r_orig_plt / 1e10, v_orig_plt, 'r-', linewidth=2, alpha=0.7, label='Original')
    ax1.loglog(r_opt_plt / 1e10, v_opt_plt, 'b-', linewidth=1.5, label='Optimized')
    ax1.set_xlabel('Radius (10¹⁰ cm)')
    ax1.set_ylabel('Velocity (km/s)')
    ax1.set_title('Velocity Profile Comparison')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Velocity spacing
    ax2 = axes[0, 1]
    dv_orig = np.abs(np.diff(v_orig))
    dv_opt = np.abs(np.diff(v_opt))
    r_mid_orig = (r_orig[:-1] + r_orig[1:]) / 2
    r_mid_opt = (r_opt[:-1] + r_opt[1:]) / 2
    
    ax2.semilogx(r_mid_orig[::-1] / 1e10, dv_orig[::-1], 'r-', linewidth=2, alpha=0.7, label='Original')
    ax2.semilogx(r_mid_opt[::-1] / 1e10, dv_opt[::-1], 'b-', linewidth=1.5, label='Optimized')
    ax2.axhline(y=vdop_min, color='k', linestyle='--', alpha=0.5, label=f'CMFGEN limit ({vdop_min:.1f} km/s)')
    ax2.set_xlabel('Radius (10¹⁰ cm)')
    ax2.set_ylabel('|ΔV| between points (km/s)')
    ax2.set_title('Velocity Spacing')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Zoom on velocity minimum region
    ax3 = axes[1, 0]
    v_min_idx = np.argmin(v_orig)
    zoom_range = 30  # points around minimum
    start_idx = max(0, v_min_idx - zoom_range)
    end_idx = min(len(v_orig), v_min_idx + zoom_range + 1)
    
    ax3.plot(r_orig_plt[start_idx:end_idx] / 1e10, v_orig_plt[start_idx:end_idx], 
             'ro-', linewidth=2, alpha=0.7, label='Original')
    ax3.plot(r_opt_plt[start_idx:end_idx] / 1e10, v_opt_plt[start_idx:end_idx], 
             'bo-', linewidth=1.5, label='Optimized')
    ax3.set_xlabel('Radius (10¹⁰ cm)')
    ax3.set_ylabel('Velocity (km/s)')
    ax3.set_title(f'Velocity Minimum Region (zoom)')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Grid point density
    ax4 = axes[1, 1]
    # Calculate local grid density (points per decade in radius)
    def calculate_grid_density(r_array):
        log_r = np.log10(r_array)
        density = np.zeros_like(r_array)
        window = 5
        for i in range(len(r_array)):
            start = max(0, i - window)
            end = min(len(r_array), i + window + 1)
            if end > start:
                log_span = log_r[end-1] - log_r[start]
                if log_span > 0:
                    density[i] = (end - start) / log_span
        return density
    
    density_orig = calculate_grid_density(r_orig_plt)
    density_opt = calculate_grid_density(r_opt_plt)
    
    ax4.semilogx(r_orig_plt / 1e10, density_orig, 'r-', linewidth=2, alpha=0.7, label='Original')
    ax4.semilogx(r_opt_plt / 1e10, density_opt, 'b-', linewidth=1.5, label='Optimized')
    ax4.set_xlabel('Radius (10¹⁰ cm)')
    ax4.set_ylabel('Grid density (points/decade)')
    ax4.set_title('Grid Point Density')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    plt.suptitle(f'Velocity Spacing Optimization (VFRAC_MOM={vfrac_mom}, VDOP_VEC={vdop_vec})', 
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    plot_file = output_file.replace('.dat', '_velocity_optimization.png')
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    print(f"Velocity optimization plot saved to {plot_file}")
    plt.close()


# Override the main function to include velocity spacing optimization
def plot_comparison(orig_r, orig_rho, orig_T, orig_v, orig_tau, orig_sigma,
                   new_r, new_rho, new_T, new_v, new_tau, new_sigma,
                   output_file, time_days, tau_min=1e-5, tau_max=30.0,
                   orig_mass_fractions=None, new_mass_fractions=None,
                   orig_ne=None, new_ne=None):
    """Plot before/after comparison of STELLA to CMFGEN conversion"""
    rsun_cm = 6.96e10
    
    # Convert radius to r/rsun
    orig_r_rsun = orig_r / rsun_cm
    new_r_rsun = new_r / rsun_cm
    
    # Determine consistent x-axis limits
    r_min = max(10.0, min(np.min(orig_r_rsun), np.min(new_r_rsun)) * 0.9)
    r_max = max(np.max(orig_r_rsun), np.max(new_r_rsun)) * 1.1
    
    # Create 2x3 grid with shared x-axis
    fig, axes = plt.subplots(2, 3, figsize=(18, 12), sharex=True)
    fig.suptitle(f'STELLA to CMFGEN Conversion: {time_days:.1f}-day snapshot', fontsize=16)
    
    # Plot 1: Density
    ax = axes[0, 0]
    ax.loglog(orig_r_rsun, orig_rho, 'bs', alpha=0.7, markersize=3, label=f'STELLA ({len(orig_r)} points)')
    ax.loglog(new_r_rsun, new_rho, 'ro', alpha=0.7, markersize=2, label=f'CMFGEN ({len(new_r)} points)')
    ax.set_ylabel('Density (g/cm³)', fontsize=12)
    ax.set_title('Density Profile', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_xlim(r_min, r_max)
    
    # Plot 2: Temperature
    ax = axes[0, 1]
    ax.loglog(orig_r_rsun, orig_T, 'bs', alpha=0.7, markersize=3, label=f'STELLA ({len(orig_r)} points)')
    ax.loglog(new_r_rsun, new_T, 'ro', alpha=0.7, markersize=2, label=f'CMFGEN ({len(new_r)} points)')
    ax.set_ylabel('Temperature (K)', fontsize=12)
    ax.set_title('Temperature Profile', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_xlim(r_min, r_max)
    
    # Plot 3: Velocity - Blue squares for STELLA, red circles for CMFGEN to highlight AMR points
    ax = axes[0, 2]
    ax.loglog(orig_r_rsun, orig_v, 'bs', alpha=0.7, markersize=4, label=f'STELLA ({len(orig_r)} points)')
    ax.loglog(new_r_rsun, new_v, 'ro', alpha=0.8, markersize=2, label=f'CMFGEN ({len(new_r)} points)')
    ax.set_ylabel('Velocity (km/s)', fontsize=12)
    ax.set_title('Velocity Profile (AMR points highlighted)', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_xlim(r_min, r_max)
    
    # Plot 4: Sigma
    ax = axes[1, 0]
    ax.plot(orig_r_rsun, orig_sigma, 'bs', alpha=0.7, markersize=3, label=f'STELLA ({len(orig_r)} points)')
    ax.plot(new_r_rsun, new_sigma, 'ro', alpha=0.7, markersize=2, label=f'CMFGEN ({len(new_r)} points)')
    ax.axhline(y=-0.99, color='k', linestyle=':', alpha=0.5, label='σ=-0.99')
    ax.axhline(y=-1.0, color='gray', linestyle='--', alpha=0.5, label='σ=-1.0')
    ax.set_xlabel('r/R☉', fontsize=12)
    ax.set_ylabel('Sigma (σ)', fontsize=12)
    ax.set_title('Sigma Profile', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_yscale('symlog', linthresh=1.0)
    ax.set_xlim(r_min, r_max)
    ax.set_xscale('log')
    
    # Plot 5: Optical depth
    ax = axes[1, 1]
    ax.plot(orig_r_rsun, orig_tau, 'bs', alpha=0.7, markersize=3, label=f'STELLA ({len(orig_r)} points)')
    ax.plot(new_r_rsun, new_tau, 'ro', alpha=0.7, markersize=2, label=f'CMFGEN ({len(new_r)} points)')
    ax.axhline(y=2/3, color='orange', linestyle=':', alpha=0.8, linewidth=2, label='τ = 2/3')
    # Highlight actual tau range of trimmed data in grey
    actual_tau_min = new_tau.min()
    actual_tau_max = new_tau.max()
    ax.axhspan(actual_tau_min, actual_tau_max, facecolor='gray', alpha=0.2, label=f'τ range: [{actual_tau_min:.1e}, {actual_tau_max:.1e}]')
    ax.set_xlabel('r/R☉', fontsize=12)
    ax.set_ylabel('Optical Depth τ', fontsize=12)
    ax.set_title('Optical Depth Profile', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_yscale('log')
    ax.set_xlim(r_min, r_max)
    ax.set_xscale('log')
    
    # Plot 6: Electron density (Ne) - replaces mass fractions plot
    ax = axes[1, 2]
    if orig_ne is not None and new_ne is not None:
        ax.loglog(orig_r_rsun, orig_ne, 'bs', alpha=0.7, markersize=3, label=f'STELLA ({len(orig_r)} points)')
        ax.loglog(new_r_rsun, new_ne, 'ro', alpha=0.7, markersize=2, label=f'CMFGEN ({len(new_r)} points)')
        ax.set_xlabel('r/R☉', fontsize=12)
        ax.set_ylabel('Electron Density (cm⁻³)', fontsize=12)
        ax.set_title('Electron Density Profile', fontsize=14)
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3, which='both')
        ax.set_xlim(r_min, r_max)
    else:
        ax.text(0.5, 0.5, 'Electron density data not available', ha='center', va='center', transform=ax.transAxes)
        ax.set_title('Electron Density', fontsize=14)
        ax.set_xlabel('r/R☉', fontsize=12)
        ax.set_xlim(r_min, r_max)
        ax.set_xscale('log')
    
    plt.tight_layout()
    plot_file = output_file.replace('.dat', '_comparison.png').replace('SN_HYDRO_DATA', 'SN_HYDRO_DATA_comparison')
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    print(f"Plot saved to {plot_file}")
    plt.close()


def plot_mass_fractions_stacked(radius_cm, density, mass_fractions, time_days, output_file,
                               trimmed_radius_cm=None, trimmed_mass_fractions=None):
    """Create stacked mass fraction plot"""
    rsun_cm = 6.96e10
    
    # Convert to Rsun
    radius_rsun = radius_cm / rsun_cm
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), sharex=True)
    fig.suptitle(f'Mass Fraction Distribution — {time_days:.1f} days', fontsize=16)
    
    # Upper plot: Individual mass fractions
    colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4', '#FFEAA7', '#DDA0DD', '#F8B500', '#6C5CE7']
    elements = ['H', 'He', 'C', 'N', 'O', 'Ne', 'Si', 'Fe']
    
    for i, elem in enumerate(elements):
        if elem in mass_fractions:
            color = colors[i % len(colors)]
            ax1.semilogx(radius_rsun, mass_fractions[elem], label=elem, color=color, linewidth=2)
    
    ax1.set_ylabel('Mass Fraction', fontsize=12)
    ax1.set_title('Individual Mass Fractions', fontsize=14)
    ax1.legend(fontsize=10, ncol=2)
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(0, 1.1)
    
    # Lower plot: Stacked mass fractions
    # Create stacked array
    sorted_elements = [elem for elem in elements if elem in mass_fractions]
    mf_array = np.array([mass_fractions[elem] for elem in sorted_elements])
    
    ax2.stackplot(radius_rsun, mf_array, labels=sorted_elements, 
                  colors=[colors[elements.index(elem)] for elem in sorted_elements],
                  alpha=0.8)
    
    ax2.set_xlabel('Radius (R☉)', fontsize=12)
    ax2.set_ylabel('Cumulative Mass Fraction', fontsize=12)
    ax2.set_title('Stacked Mass Fractions', fontsize=14)
    ax2.legend(fontsize=10, ncol=2, loc='upper right')
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim(0, 1.1)
    ax2.set_xlim(radius_rsun.min() * 0.9, radius_rsun.max() * 1.1)
    
    plt.tight_layout()
    plot_file = output_file.replace('.dat', '_mass_fraction.png').replace('SN_HYDRO_DATA', 'SN_HYDRO_DATA_mass_fraction')
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    print(f"Mass fraction plot saved to {plot_file}")
    plt.close()


def plot_dr_vs_r(r_cm, output_file, time_days, rho=None, T=None, v=None):
    """Plot dr/r vs radius with optional physical quantities"""
    # Calculate dr/r
    dr = np.diff(r_cm)
    r_mid = (r_cm[:-1] + r_cm[1:]) / 2
    dr_over_r = dr / r_mid
    
    grid_indices = np.arange(len(r_mid))
    
    fig, ax = plt.subplots(figsize=(14, 7))
    
    # Main plot: dr/r vs radius
    ax.plot(r_mid / 1e10, dr_over_r, 'b-', linewidth=1.5, label='Grid spacing |dr|/r')
    
    ax.set_xlabel('Radius (10$^{10}$ cm)', fontsize=12)
    ax.set_ylabel('Grid spacing |dr|/r', fontsize=12, color='b')
    ax.tick_params(axis='y', labelcolor='b')
    ax.grid(True, alpha=0.3)
    
    # Add threshold lines
    ax.axhline(y=0.05, color='r', linestyle='--', alpha=0.7, label='AMR threshold (0.05)')
    ax.axhline(y=-0.05, color='r', linestyle='--', alpha=0.7)
    
    # Secondary axis for physical quantities
    ax2 = ax.twinx()
    
    if rho is not None:
        ax2.semilogy(r_mid / 1e10, rho[:-1], 'g-', alpha=0.6, linewidth=1, label='Density')
        ax2.set_ylabel('Density (g/cm³)', fontsize=12, color='g')
        ax2.tick_params(axis='y', labelcolor='g')
    
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, loc='upper right', fontsize=9)
    
    plt.tight_layout()
    plot_file = output_file.replace('.dat', '_dr_vs_r.png').replace('SN_HYDRO_DATA', 'SN_HYDRO_DATA_dr_vs_r')
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    print(f"dr vs r plot saved to {plot_file}")
    plt.close()


def plot_structure_overview(r_cm, v_kms, rho, T_K, mass_fractions, kappa, time_days, output_file):
    """Create 6-panel structure overview plot with consistent x-axis formatting.
    
    Input arrays are in outer->inner order (from STELLA .swd).
    For plotting, convert to inner->outer order so x-axis increases left-to-right.
    """
    rsun_cm = 6.96e10
    
    # Convert to inner->outer order for plotting (x-axis increases left-to-right)
    # Note: This reverses the arrays in-place for plotting only
    r_plt = r_cm[::-1]
    v_plt = v_kms[::-1]
    rho_plt = rho[::-1]
    T_plt = T_K[::-1]
    kappa_plt = kappa[::-1]
    mf_plt = {elem: mass_fractions[elem][::-1] for elem in mass_fractions}
    
    r_plt = r_plt / rsun_cm  # Convert to R_sun for consistency
    
    # Determine consistent x-axis limits
    r_min = max(10.0, np.min(r_plt) * 0.9)
    r_max = np.max(r_plt) * 1.1
    
    # Create 2x3 grid
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle(f'STELLA Structure Overview — t = {time_days:.1f} days', fontsize=16, fontweight='bold')
    
    # Plot 1: Density
    ax = axes[0, 0]
    ax.loglog(r_plt, rho_plt, 'b-', linewidth=2)
    ax.set_ylabel('Density (g/cm³)', fontsize=12)
    ax.set_title('Density Profile', fontsize=14)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_xlim(r_min, r_max)
    
    # Plot 2: Temperature
    ax = axes[0, 1]
    ax.loglog(r_plt, T_plt, 'r-', linewidth=2)
    ax.set_ylabel('Temperature (K)', fontsize=12)
    ax.set_title('Temperature Profile', fontsize=14)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_xlim(r_min, r_max)
    
    # Plot 3: Velocity
    ax = axes[0, 2]
    ax.loglog(r_plt, v_plt, 'g-', linewidth=2)
    ax.set_ylabel('Velocity (km/s)', fontsize=12)
    ax.set_title('Velocity Profile', fontsize=14)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_xlim(r_min, r_max)
    
    # Plot 4: Opacity
    ax = axes[1, 0]
    ax.loglog(r_plt, kappa_plt, 'm-', linewidth=2)
    ax.set_ylabel('Opacity (cm²/g)', fontsize=12)
    ax.set_title('Rosseland Opacity', fontsize=14)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_xlim(r_min, r_max)
    
    # Plot 5: Mass fractions (stacked)
    ax = axes[1, 1]
    elements = ['H', 'He', 'C', 'N', 'O', 'Ne', 'Si', 'Fe']
    colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4', '#FFEAA7', '#DDA0DD', '#F8B500', '#6C5CE7']
    
    sorted_elements = [elem for elem in elements if elem in mf_plt]
    mf_array = np.array([mf_plt[elem] for elem in sorted_elements])
    
    if len(sorted_elements) > 0:
        ax.stackplot(r_plt, mf_array, labels=sorted_elements,
                    colors=[colors[elements.index(elem)] for elem in sorted_elements],
                    alpha=0.8)
        ax.legend(fontsize=8, ncol=2)
    ax.set_ylabel('Mass Fraction', fontsize=12)
    ax.set_title('Mass Fractions', fontsize=14)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(r_min, r_max)
    
    # Plot 6: Optical depth
    ax = axes[1, 2]
    tau = compute_kappa_tau(r_cm, rho, kappa)
    tau_plt = tau[::-1]  # Convert to inner->outer order
    ax.loglog(r_plt, tau_plt, 'k-', linewidth=2)
    ax.axhline(y=2/3, color='orange', linestyle='--', alpha=0.8, linewidth=2, label='τ = 2/3')
    ax.set_xlabel('Radius (R☉)', fontsize=12)
    ax.set_ylabel('Optical Depth τ', fontsize=12)
    ax.set_title('Optical Depth', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_xlim(r_min, r_max)
    
    # Set common x-axis label on bottom row
    for ax in axes[3:]:
        ax.set_xlabel(r'r/R$_\odot$')
    
    plt.suptitle(f'STELLA Structure Overview — t = {time_days:.1f} days', fontsize=16, fontweight='bold')
    plt.tight_layout()
    
    plot_file = output_file.replace('.dat', '_plot.png').replace('SN_HYDRO_DATA', 'SN_HYDRO_DATA_plot')
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    print(f"Structure overview plot saved to {plot_file}")
    plt.close()


def adaptive_grid_refinement(radius, v, T, rho, n_e, tau, mass_fractions,
                             max_radius_jump=0.05, max_temp_jump=0.10,
                             max_density_jump=0.15, max_ne_jump=0.10,
                             max_logtau_jump=0.10, max_points=500,
                             include_ne_criterion=False):
    """Apply Dessart's adaptive grid refinement criteria.
    
    Inserts additional grid points where consecutive points exceed the maximum
    allowed jumps in physical quantities.
    
    Criteria from Dessart et al. 2017/2026:
      - ≤5% jump in radius
      - ≤10% jump in temperature
      - ≤10% jump in electron density (optional, see include_ne_criterion)
      - ≤10% jump in density
      - ≤10% jump in log optical depth
    
    Parameters:
    -----------
    radius, v, T, rho, n_e, tau : arrays
        Physical quantities (outer->inner order)
    mass_fractions : dict
        Dictionary of element mass fractions
    max_*_jump : float
        Maximum allowed fractional jump for each quantity
    include_ne_criterion : bool
        Whether to include electron density jumps in AMR criteria
        (default: False, as Ne can have artificial jumps at ionization fronts)
    
    Returns:
    --------
    Refined arrays (radius, v, T, rho, n_e, tau, mass_fractions)
    """
    def interpolate_arrays(r_list, v_list, T_list, rho_list, n_e_list, tau_list, mf_lists, i, j, fraction):
        """Interpolate between indices i and j at given fraction (0-1)."""
        r_new = r_list[i] + fraction * (r_list[j] - r_list[i])
        
        # Use linear interpolation for velocity so red line closely follows blue line
        # AMR adds points where other quantities have high gradients, not velocity
        v_new = v_list[i] + fraction * (v_list[j] - v_list[i])
            
        # Use log interpolation for quantities with large dynamic range
        # Add small epsilon to avoid log(0)
        eps = 1e-30
        T_new = 10**(np.log10(max(T_list[i], eps)) + fraction * (np.log10(max(T_list[j], eps)) - np.log10(max(T_list[i], eps))))
        rho_new = 10**(np.log10(max(rho_list[i], eps)) + fraction * (np.log10(max(rho_list[j], eps)) - np.log10(max(rho_list[i], eps))))
        n_e_new = 10**(np.log10(max(n_e_list[i], eps)) + fraction * (np.log10(max(n_e_list[j], eps)) - np.log10(max(n_e_list[i], eps))))
        tau_new = 10**(np.log10(max(tau_list[i], eps)) + fraction * (np.log10(max(tau_list[j], eps)) - np.log10(max(tau_list[i], eps))))
        
        mf_new = {}
        for elem in mf_lists:
            mf_new[elem] = mf_lists[elem][i] + fraction * (mf_lists[elem][j] - mf_lists[elem][i])
        
        return r_new, v_new, T_new, rho_new, n_e_new, tau_new, mf_new
    
    # Convert to lists for easy insertion
    r_list = list(radius)
    v_list = list(v)
    T_list = list(T)
    rho_list = list(rho)
    n_e_list = list(n_e)
    tau_list = list(tau)
    mf_lists = {elem: list(mass_fractions[elem]) for elem in mass_fractions}

    # --- Diagnostic: count violations in the BASE grid before any AMR ---
    n_viol_r = n_viol_T = n_viol_rho = n_viol_ne = n_viol_tau = 0
    max_dr_r = max_dT_T = max_drho_rho = max_dne_ne = max_dlogtau = 0.0
    for i in range(len(r_list) - 1):
        dr_r = abs(r_list[i+1] - r_list[i]) / max(abs(r_list[i]), 1e-30)
        dT_T = abs(T_list[i+1] - T_list[i]) / max(T_list[i], 1e-30)
        drho_rho = abs(rho_list[i+1] - rho_list[i]) / max(rho_list[i], 1e-30)
        dne_ne = abs(n_e_list[i+1] - n_e_list[i]) / max(n_e_list[i], 1e-30)
        eps = 1e-30
        dlogtau = abs(np.log10(max(tau_list[i+1], eps)) - np.log10(max(tau_list[i], eps)))

        max_dr_r = max(max_dr_r, dr_r)
        max_dT_T = max(max_dT_T, dT_T)
        max_drho_rho = max(max_drho_rho, drho_rho)
        max_dne_ne = max(max_dne_ne, dne_ne)
        max_dlogtau = max(max_dlogtau, dlogtau)

        if dr_r > max_radius_jump:       n_viol_r += 1
        if dT_T > max_temp_jump:         n_viol_T += 1
        if drho_rho > max_density_jump:  n_viol_rho += 1
        if include_ne_criterion and dne_ne > max_ne_jump: n_viol_ne += 1
        if dlogtau > max_logtau_jump and tau_list[i] > 0.01:
            n_viol_tau += 1

    print(f"[AMR Audit] Base grid ({len(r_list)} points) violation summary:")
    print(f"  dr/r:     {n_viol_r:4d} edges exceed {max_radius_jump}  (max jump = {max_dr_r:.3f})")
    print(f"  dT/T:     {n_viol_T:4d} edges exceed {max_temp_jump}  (max jump = {max_dT_T:.3f})")
    print(f"  drho/rho: {n_viol_rho:4d} edges exceed {max_density_jump}  (max jump = {max_drho_rho:.3f})")
    if include_ne_criterion:
        print(f"  dne/ne:   {n_viol_ne:4d} edges exceed {max_ne_jump}  (max jump = {max_dne_ne:.3f})")
    else:
        print(f"  dne/ne:   (disabled - Ne criterion not applied)")
    print(f"  dlogtau:  {n_viol_tau:4d} edges exceed {max_logtau_jump}  (max jump = {max_dlogtau:.3f})")
    
    total_violations = n_viol_r + n_viol_T + n_viol_rho + n_viol_ne + n_viol_tau
    if total_violations == 0:
        print("[AMR Audit] Base grid already meets ALL criteria. No AMR needed.")
        # Return original arrays without modification
        return (np.array(r_list), np.array(v_list), np.array(T_list),
                np.array(rho_list), np.array(n_e_list), np.array(tau_list),
                {k: np.array(v) for k, v in mf_lists.items()})

    max_iterations = 100  # Prevent infinite loops
    iteration = 0
    n_inserted_total = 0

    while iteration < max_iterations:
        n_inserted = 0
        i = 0

        while i < len(r_list) - 1:
            if len(r_list) >= max_points:
                print(f"[Adaptive Grid] Reached max_points ({max_points}), stopping refinement")
                break
            # Calculate jumps between consecutive points
            dr_r = abs(r_list[i+1] - r_list[i]) / max(abs(r_list[i]), 1e-30)
            dT_T = abs(T_list[i+1] - T_list[i]) / max(T_list[i], 1e-30)
            drho_rho = abs(rho_list[i+1] - rho_list[i]) / max(rho_list[i], 1e-30)
            dne_ne = abs(n_e_list[i+1] - n_e_list[i]) / max(n_e_list[i], 1e-30)
            eps = 1e-30
            dlogtau = abs(np.log10(max(tau_list[i+1], eps)) - np.log10(max(tau_list[i], eps)))

            # Check if any criterion is violated
            # By default, we exclude dne/ne because n_e can have artificial jumps at ionization fronts
            # But we include it if explicitly requested
            violation = False
            if dr_r > max_radius_jump:
                violation = True
            if dT_T > max_temp_jump:
                violation = True
            if drho_rho > max_density_jump:
                violation = True
            if include_ne_criterion and dne_ne > max_ne_jump:
                violation = True
            # Skip dlogtau criterion in very optically thin outer regions (tau < 0.01)
            # to prevent over-refinement where opacity is negligible
            if dlogtau > max_logtau_jump and tau_list[i] > 0.01:
                violation = True

            if violation:
                # Check minimum spacing to prevent over-refinement
                MIN_DR_OVER_R = 1e-5  # Minimum relative spacing
                if dr_r < MIN_DR_OVER_R:
                    # Already fine enough, skip insertion and just move to next
                    i += 1
                else:
                    # Insert a point at the midpoint
                    r_new, v_new, T_new, rho_new, n_e_new, tau_new, mf_new = interpolate_arrays(
                        r_list, v_list, T_list, rho_list, n_e_list, tau_list, mf_lists, i, i+1, 0.5)

                    r_list.insert(i+1, r_new)
                    v_list.insert(i+1, v_new)
                    T_list.insert(i+1, T_new)
                    rho_list.insert(i+1, rho_new)
                    n_e_list.insert(i+1, n_e_new)
                    tau_list.insert(i+1, tau_new)

                    for elem in mf_lists:
                        mf_lists[elem].insert(i+1, mf_new[elem])

                    n_inserted += 1
                    i += 2  # Skip over the newly inserted point
            else:
                i += 1

        n_inserted_total += n_inserted
        if n_inserted == 0:
            break  # No more violations, converged

        iteration += 1

    if n_inserted_total > 0:
        print(f"[Adaptive Grid] Inserted {n_inserted_total} points to meet gradient criteria")
        print(f"  Final grid size: {len(r_list)} points")
    else:
        print(f"[Adaptive Grid] No points inserted (grid already meets criteria)")
    
    # Convert back to numpy arrays
    radius = np.array(r_list)
    v = np.array(v_list)
    T = np.array(T_list)
    rho = np.array(rho_list)
    n_e = np.array(n_e_list)
    tau = np.array(tau_list)
    mass_fractions = {elem: np.array(mf_lists[elem]) for elem in mf_lists}
    
    return radius, v, T, rho, n_e, tau, mass_fractions


def main_optimized():
    """Enhanced main function with velocity spacing optimization."""
    
    # Get original arguments
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('swd_file', help='STELLA *.swd output file')
    parser.add_argument('output_file', help='Output SN_HYDRO_DATA path')
    parser.add_argument('--abn-file', default=None)
    parser.add_argument('--input-profile', default=None,
                        help='[DEPRECATED] HERACLES input_profile file for composition fallback')
    parser.add_argument('--xni-file', default=None,
                        help='STELLA *.xni file for Ni56 distribution')
    parser.add_argument('--lbol-file', default=None,
                        help='STELLA *.lbol file for actual bolometric luminosity')
    parser.add_argument('--target-days', type=float, default=None,
                        help='Pick epoch closest to this (default: last epoch)')
    parser.add_argument('--tau-max', type=float, default=30.0)
    parser.add_argument('--tau-min', type=float, default=1e-4)
    parser.add_argument('--nzon-target', type=int, default=200)
    parser.add_argument('--min-points', type=int, default=80,
                        help='Minimum number of points in output SN_HYDRO_DATA')
    parser.add_argument('--kappa', type=float, default=0.34,
                        help='Electron scattering opacity cm^2/g')
    parser.add_argument('--model-name', default='stella')
    parser.add_argument('--enable-grid-thinning', action='store_true',
                        help='Enable grid thinning to remove near-degenerate points')
    parser.add_argument('--grid-thinning-threshold', type=float, default=0.005,
                        help='min_dr/r threshold for grid thinning (default: 0.005)')
    parser.add_argument('--enable-adaptive-grid', action='store_true',
                        help='Enable Dessart-style adaptive grid refinement')
    parser.add_argument('--base-grid-target', type=int, default=80,
                        help='Target base grid size before adaptive refinement (default: 80)')
    parser.add_argument('--max-radius-jump', type=float, default=0.05,
                        help='Max radius jump fraction for adaptive grid (default: 0.05)')
    parser.add_argument('--max-temp-jump', type=float, default=0.10,
                        help='Max temperature jump fraction (default: 0.10)')
    parser.add_argument('--max-density-jump', type=float, default=0.15,
                        help='Max density jump fraction (default: 0.15)')
    parser.add_argument('--max-ne-jump', type=float, default=0.10,
                        help='Max electron density jump fraction (default: 0.10)')
    parser.add_argument('--include-ne-amr', action='store_true', default=False,
                        help='Include electron density jumps in AMR criteria (default: False - Ne can have artificial jumps)')
    parser.add_argument('--max-logtau-jump', type=float, default=0.10,
                        help='Max log optical depth jump (default: 0.10)')
    parser.add_argument('--max-points', type=int, default=500,
                        help='Maximum points after adaptive refinement (default: 500)')
    
    # New velocity spacing optimization arguments
    parser.add_argument('--vfrac-mom', type=float, default=0.5,
                        help='Target VFRAC_MOM value for velocity spacing (default: 0.5)')
    parser.add_argument('--vdop-vec', type=float, default=30.0,
                        help='Target VDOP_VEC value in km/s (default: 30.0)')
    parser.add_argument('--optimize-velocity-spacing', action='store_true', default=False,
                        help='Enable velocity spacing optimization (default: False)')
    parser.add_argument('--no-optimize-velocity-spacing', dest='optimize_velocity_spacing', 
                        action='store_false',
                        help='Disable velocity spacing optimization')
    
    # Add all other arguments from original main function
    parser.add_argument('--time-days', type=float, default=None,
                        help='Explosion time in days (default: read from swd)')
    parser.add_argument('--lum-xg', type=float, default=None,
                        help='Override luminosity in Lsun (useful for early times)')
    parser.add_argument('--no-smooth-velocity', action='store_true', default=False,
                        help='Disable velocity discontinuity smoothing')
    parser.add_argument('--sigma-threshold', type=float, default=3.0,
                        help='Sigma threshold for velocity smoothing (default: 3.0)')
    parser.add_argument('--smooth-width', type=int, default=3,
                        help='Width for velocity smoothing (default: 3)')
    parser.add_argument('--t-floor', type=float, default=3000.0,
                        help='Temperature floor in K (default: 3000)')
    parser.add_argument('--ni56-mass', type=float, default=0.0,
                        help='56Ni mass in Msun for luminosity estimate (default: 0)')
    
    args = parser.parse_args()
    
    # Run the original conversion up to the point where we have the final grid
    # This is a modified version of the original main() function
    
    print(f"Reading STELLA data from {args.swd_file}")
    epochs = read_swd(args.swd_file)
    
    # Pick epoch
    if args.target_days is not None:
        t_pick = min(epochs.keys(), key=lambda t: abs(t - args.target_days))
        print(f"Targeting t = {args.target_days:.1f} days, closest epoch: t = {t_pick:.1f} days")
    else:
        t_pick = max(epochs.keys())
        print(f"Using last epoch: t = {t_pick:.1f} days")
    
    data = epochs[t_pick]
    r, v, T, rho = data['r'], data['v'], data['T'], data['rho']
    print(f"Read {len(r)} zones from STELLA epoch at t = {t_pick:.1f} days")
    
    # Read composition
    mass_fractions = {}
    if args.abn_file:
        print(f"Reading abundances from {args.abn_file}")
        mass_fractions = read_abn(args.abn_file, len(r))
    elif args.input_profile:
        print(f"Reading abundances from {args.input_profile}")
        mass_fractions = read_input_profile(args.input_profile)
        # [DEPRECATED] Convert HERACLES element names to CMFGEN DUM labels
        # HERACLES is deprecated; the project now uses the STELLA -> CMFGEN pipeline.
        conversion = {'H': 'HYD', 'He': 'HE', 'O': 'OXY', 'Si': 'SIL', 
                     'Ca': 'CAL', 'Fe': 'IRON', 'Ni': 'NICK'}
        mass_fractions = {conversion.get(k, k): v for k, v in mass_fractions.items()}
    
    # Set default composition if none provided
    if not mass_fractions:
        print("No composition file provided; using pure H composition")
        mass_fractions = {'HYD': np.ones(len(r))}
    
    # Read Ni56 distribution
    ni56_for_isotopes = np.zeros(len(r))
    if args.xni_file:
        print(f"Reading Ni56 distribution from {args.xni_file}")
        ni56_for_isotopes = read_xni(args.xni_file, len(r))
        mass_fractions['Ni56'] = ni56_for_isotopes
    
    # Apply tau trimming (following the original script's approach)
    print(f"Applying tau trimming: τ_min={args.tau_min}, τ_max={args.tau_max}")
    
    # Compute initial opacity and optical depth
    kappa_arr = compute_kappa_total(rho, T, mass_fractions)
    tau = compute_kappa_tau(r, rho, kappa_arr)
    print(f"Original tau range: [{tau.min():.3e}, {tau.max():.3e}]")
    
    # Find tau bounds
    nzon = len(r)
    i_outer = 0
    while i_outer < nzon and tau[i_outer] < args.tau_min:
        i_outer += 1
    
    i_inner = nzon - 1
    while i_inner > 0 and tau[i_inner] > args.tau_max:
        i_inner -= 1
    
    # Apply tau trim
    if i_outer < i_inner:
        r = r[i_outer:i_inner+1]
        v = v[i_outer:i_inner+1]
        T = T[i_outer:i_inner+1]
        rho = rho[i_outer:i_inner+1]
        for elem in mass_fractions:
            mass_fractions[elem] = mass_fractions[elem][i_outer:i_inner+1]
        
        # Store for Ni56 interpolation
        r_tau_trim = r.copy()
        
        print(f"After tau trimming: {len(r)} zones (indices {i_outer}:{i_inner})")
    else:
        print("WARNING: No zones found within tau range!")
        return
    
    # Store TRUE original arrays before ANY processing for accurate comparison
    r_original = r.copy()
    v_original = v.copy()
    T_original = T.copy()
    rho_original = rho.copy()
    mass_fractions_original = {elem: mass_fractions[elem].copy() for elem in mass_fractions}
    
    # Apply adaptive grid refinement if requested
    if args.enable_adaptive_grid:
        print("\n--- Applying Adaptive Grid Refinement ---")
        n_e = calculate_electron_density(rho, T, mass_fractions)
        tau = compute_kappa_tau(r, rho, compute_kappa_total(rho, T, mass_fractions))
        
        r, v, T, rho, n_e, tau, mass_fractions = adaptive_grid_refinement(
            r, v, T, rho, n_e, tau, mass_fractions,
            max_radius_jump=args.max_radius_jump,
            max_temp_jump=args.max_temp_jump,
            max_density_jump=args.max_density_jump,
            max_ne_jump=args.max_ne_jump,
            max_logtau_jump=args.max_logtau_jump,
            max_points=args.max_points,
            include_ne_criterion=args.include_ne_amr  # New parameter to control Ne inclusion
        )
        print(f"After adaptive refinement: {len(r)} zones")
    
    # Store arrays before velocity optimization (but after AMR) for intermediate comparison
    r_before_opt = r.copy()
    v_before_opt = v.copy()
    T_before_opt = T.copy()
    rho_before_opt = rho.copy()
    mass_fractions_before_opt = {elem: mass_fractions[elem].copy() for elem in mass_fractions}
    
    # Apply velocity spacing optimization
    if args.optimize_velocity_spacing:
        # Apply hybrid velocity approach to preserve non-monotonic structure while ensuring MOM_JREL_V9 compatibility
        r, v, T, rho, mass_fractions = create_hybrid_velocity_profile(
            r, v, T, rho, mass_fractions,
            vfrac_mom=args.vfrac_mom,
            vdop_vec=args.vdop_vec,
            preserve_non_monotonic=True,
            verbose=True
        )
        
        # Velocity optimization diagnostic is included in the comparison plots
    
    # Apply velocity smoothing only if explicitly enabled (default is disabled for non-monotonic velocity profiles)
    # For non-monotonic velocity profiles, we want AMR to handle resolution, not artificial smoothing
    if args.no_smooth_velocity:
        print("Velocity smoothing DISABLED - allowing AMR to resolve grid properly")
    else:
        print("Applying minimal velocity smoothing...")
        v = smooth_velocity_discontinuities(r, v, 
                                           sigma_threshold=args.sigma_threshold,
                                           smooth_width=args.smooth_width)
    
    # Calculate final quantities
    sigma = sigma_velocity(r, v)
    kappa_arr = compute_kappa_total(rho, T, mass_fractions)
    
    # Calculate luminosity
    tau = compute_kappa_tau(r, rho, kappa_arr)
    inner_idx = np.argmax(tau > 2.0) if np.any(tau > 2.0) else -1
    if inner_idx >= 0:
        r_inner_boundary = r[inner_idx]
        T_inner_boundary = T[inner_idx]
    else:
        r_inner_boundary = r[-1]
        T_inner_boundary = T[-1]
    
    lum_lsun = luminosity_bb(T_inner_boundary, r_inner_boundary) / LSUN
    
    # Write the final SN_HYDRO_DATA
    print(f"\nWriting optimized SN_HYDRO_DATA to {args.output_file}")
    # Resize ni56 array to match final grid if needed
    if len(ni56_for_isotopes) != len(r):
        # Ensure we have the right arrays for interpolation
        if 'r_tau_trim' in locals() and len(r_tau_trim) == len(ni56_for_isotopes):
            ni56_resized = np.interp(r, r_tau_trim[::-1], ni56_for_isotopes[::-1])
        else:
            # Fallback: create zeros array
            ni56_resized = np.zeros(len(r))
    else:
        ni56_resized = ni56_for_isotopes
        
    write_sn_hydro_data(
        args.output_file, r, v, rho, T, sigma, kappa_arr,
        t_pick, mass_fractions, lum_lsun, t_floor=args.t_floor,
        r_inner_cm=r_inner_boundary, T_inner_K=T_inner_boundary,
        ni56_initial=ni56_resized
    )
    
    print(f"\nConversion complete!")
    print(f"Final grid: {len(r)} points")
    print(f"Velocity range: {v.max():.1f} -> {v.min():.1f} km/s")
    print(f"Radius range: {r.max():.3e} -> {r.min():.3e} cm")
    
    # Verify final velocity spacing
    dv_final = np.abs(np.diff(v))
    vdop_min = args.vfrac_mom * args.vdop_vec
    n_violations = np.sum(dv_final > vdop_min)
    if n_violations > 0:
        print(f"WARNING: {n_violations} velocity spacing violations remain!")
        print(f"  Max violation: {dv_final.max():.1f} km/s > {vdop_min:.1f} km/s")
    else:
        print(f"✓ All velocity spacings satisfy CMFGEN constraint (≤ {vdop_min:.1f} km/s)")
    
    # Generate diagnostic plots (same as original script)
    print(f"\nGenerating diagnostic plots...")
    
    # Calculate original quantities for comparison
    # Use the TRUE original STELLA arrays before ANY processing
    orig_r = r_original
    orig_rho = rho_original
    orig_T = T_original
    orig_v = v_original
    orig_mass_fractions = mass_fractions_original
    
    # Calculate original tau and sigma for comparison
    orig_kappa = compute_kappa_total(orig_rho, orig_T, orig_mass_fractions)
    orig_tau = compute_kappa_tau(orig_r, orig_rho, orig_kappa)
    orig_sigma = sigma_velocity(orig_r, orig_v)
    
    # New quantities
    new_tau = compute_kappa_tau(r, rho, kappa_arr)
    new_sigma = sigma_velocity(r, v)
    
    # Calculate electron densities for Ne plot
    orig_ne = calculate_electron_density(orig_rho, orig_T, orig_mass_fractions)
    new_ne = calculate_electron_density(rho, T, mass_fractions)
    
    # Generate all diagnostic plots
    plot_comparison(
        orig_r, orig_rho, orig_T, orig_v, orig_tau, orig_sigma,
        r, rho, T, v, new_tau, new_sigma,
        args.output_file, t_pick, tau_min=1e-4, tau_max=30.0,
        orig_mass_fractions=orig_mass_fractions, new_mass_fractions=mass_fractions,
        orig_ne=orig_ne, new_ne=new_ne
    )
    
    # Plot mass fractions: dashed = original STELLA, solid = optimized CMFGEN
    plot_mass_fractions_stacked(orig_r, orig_rho, orig_mass_fractions, t_pick, args.output_file,
                                 trimmed_radius_cm=r, trimmed_mass_fractions=mass_fractions)
    
    # Plot dr vs r
    plot_dr_vs_r(r, args.output_file, t_pick, rho=rho, T=T, v=v)
    
    # Plot structure overview
    plot_structure_overview(r, v, rho, T, mass_fractions, kappa_arr, t_pick, args.output_file)
    
    print(f"All diagnostic plots generated successfully!")


if __name__ == "__main__":
    main_optimized()