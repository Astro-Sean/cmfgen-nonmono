#!/usr/bin/env python3
"""
Convert HERACLES snapshot to SN_HYDRO_DATA format for CMFGEN
Uses global multi-physics gradient-based adaptive mesh refinement (AMR)

AMR Strategy:
1. Temperature gradient AMR (Step 6.5): Global |dlnT/dlnr| > 0.5 for shock boundaries
2. Multi-physics gradient AMR (Step 6.5b): Global coverage with priority weighting
   - Applied across entire grid (all depths)
   - Density: |dlnrho/dlnr| > 0.25 (20% tighter in problematic region)
   - Velocity: |dlnv/dlnr| > 0.15 (20% tighter in problematic region)
   - Temperature: |dlnT/dlnr| > 0.30 (20% tighter in problematic region)
   - Problematic region (log tau ~ -6.5 to -5.5) gets priority refinement
3. Outer boundary spacing (Step 6.6): Fixes abrupt dr changes
4. Optical depth enforcement (Step 7): Ensures dtau > 0.01 in outer 10%

This ensures density, velocity, and temperature gradients are well resolved
globally, with extra refinement where CMFGEN encounters numerical instabilities.
"""

import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.interpolate import CubicSpline
from scipy.signal import savgol_filter

# Set up HERACLES paths
heracles_root = os.environ.get("HERACLES_ROOT", "/ptmp/sbrennan/heracles")
heracles_dir = '/ptmp/sbrennan/heracles/explosion_hydro_1e51'
target_days = 10

# Parse command-line arguments
if len(sys.argv) > 1:
    heracles_dir = sys.argv[1]
if len(sys.argv) > 2:
    try:
        target_days = float(sys.argv[2])
    except ValueError:
        target_days = 10
if len(sys.argv) > 3:
    output_file = sys.argv[3]
else:
    output_file = 'SN_HYDRO_DATA'
if len(sys.argv) > 4:
    try:
        tau_min = float(sys.argv[4])
    except ValueError:
        tau_min = None
else:
    tau_min = None
if len(sys.argv) > 5:
    try:
        tau_max = float(sys.argv[5])
    except ValueError:
        tau_max = None
else:
    tau_max = None

# Add HERACLES visualization tools to path
visu1 = os.path.join(heracles_root, "heracles-2.0", "visu", "python", "data")
visu2 = os.path.join(heracles_root, "heracles-2.0", "visu", "python")
if str(visu1) not in sys.path:
    sys.path.insert(0, str(visu1))
if str(visu2) not in sys.path:
    sys.path.insert(0, str(visu2))

import rd_heracles as rdh

def find_closest_heracles_snapshot(target_days, heracles_dir):
    """Find HERACLES snapshot closest to target time using rd_heracles to read actual simulation time"""
    snapshot_times = []
    for snap_dir in os.listdir(heracles_dir):
        snap_path = os.path.join(heracles_dir, snap_dir)
        if os.path.isdir(snap_path) and snap_dir.isdigit():
            snap_id = int(snap_dir)
            # Read actual simulation time using rd_heracles
            try:
                d = rdh.HeraclesData(idump=snap_id, form='bin', verbose=False,
                                     dir=heracles_dir)
                time_s = float(d.time)
                time_days = time_s / 86400.0
                snapshot_times.append((snap_path, time_days))
            except Exception as e:
                print(f"Warning: Could not read time for snapshot {snap_id}: {e}")
                pass
    
    # Find closest snapshot
    closest = None
    min_diff = float('inf')
    for snap_dir, time_days in snapshot_times:
        diff = abs(time_days - target_days)
        if diff < min_diff:
            min_diff = diff
            closest = (snap_dir, time_days, diff)
    
    if closest:
        print(f"Found closest snapshot: {os.path.basename(closest[0])} at {closest[1]:.2f} days (diff: {closest[2]:.2f} days)")
    else:
        print(f"Warning: Could not find any snapshots with readable time in {heracles_dir}")
    
    return closest

def calculate_temperature_with_radiation(density, velocity_cms, e_tot, mu_loc):
    """
    Solve for temperature T [K] given:
      density [g/cm³], velocity_cms [cm/s], e_tot [erg/cm³], mu_loc [amu]
    Uses Newton-Raphson to enforce:
      e_tot = e_kin + e_int = 0.5*rho*v² + e_gas + e_rad
      e_gas = P_gas/(gamma-1) = (rho*k_B*T)/(mu*m_H)/(gamma-1)
      e_rad = a*T⁴ (volumetric energy density, erg/cm³)
    """
    KB = 1.38062e-16    # erg/K
    UMA = 1.660531e-24  # g (atomic mass unit)
    a_rad = 7.5657e-15  # erg cm⁻³ K⁻⁴
    GAMMA = 5.0/3.0
    m_H = UMA

    # Add safety floors to prevent division by zero
    density_safe = np.maximum(density, 1e-30)
    mu_safe = np.maximum(mu_loc, 1e-30)

    # Kinetic energy (erg/cm³)
    e_kin = 0.5 * density_safe * velocity_cms**2
    # Internal energy (erg/cm³) - ensure positive
    e_int = np.maximum(e_tot - e_kin, 1e-30)

    # Initial guess: ideal gas only (no radiation)
    T_guess = e_int * (GAMMA-1) * mu_safe * UMA / (density_safe * KB)
    T_guess = np.maximum(T_guess, 3000.0)  # Avoid too-low initial guess (global T_MIN = 3000K for NLTE)

    # Vectorized Newton-Raphson (20 iterations max)
    T = T_guess.copy()
    for _ in range(20):
        # Radiation energy density (volumetric, erg/cm³) - FIXED: removed /density
        e_rad = a_rad * T**4
        # Gas energy density (from total internal energy)
        e_gas = e_int - e_rad
        # Gas pressure (from energy)
        P_gas = (GAMMA-1) * e_gas
        # Ideal gas law pressure
        P_gas_ideal = density_safe * KB * T / (mu_safe * UMA)

        # Residual: P_gas (from energy) - P_gas (from ideal gas)
        residual = P_gas - P_gas_ideal

        # Derivative of residual w.r.t. T
        de_rad_dT = 4 * a_rad * T**3  # FIXED: removed /density
        de_gas_dT = -de_rad_dT
        dP_gas_dT = (GAMMA-1) * de_gas_dT
        dP_gas_ideal_dT = density_safe * KB / (mu_safe * UMA)
        d_residual_dT = dP_gas_dT - dP_gas_ideal_dT

        # Newton step (avoid division by zero)
        step = np.where(np.abs(d_residual_dT) > 1e-30,
                       residual / d_residual_dT,
                       0.0)
        T_new = T - step

        # Check convergence
        if np.allclose(T_new, T, rtol=1e-4, atol=1.0):
            break
        T = T_new

    # Add debug diagnostics before clipping
    n_negative = np.sum(T < 0)
    n_nonfinite = np.sum(~np.isfinite(T))
    n_above_max = np.sum(T > 5.0e5)
    n_below_min = np.sum(T < 3000.0)  # Match VADAT T_MIN=0.3 (3000 K for NLTE)
    if n_negative > 0 or n_nonfinite > 0 or n_below_min > 0 or n_above_max > 0:
        print(f"  Temperature solver diagnostics:")
        print(f"    Negative values: {n_negative}")
        print(f"    Non-finite values: {n_nonfinite}")
        print(f"    Values above 5e5 K: {n_above_max}")
        print(f"    Values below 3000 K: {n_below_min}")
        print(f"    Temperature range before clipping: {T.min():.2e} to {T.max():.2e} K")

    # Clip to physical range (3000 K to 5e5 K).
    # Upper limit: CMFGEN LTE opacity tables / GFF table become unreliable at very high T;
    # the post-shock reverse-shock spike (up to ~1.5e6 K) is deeply buried and CMFGEN's
    # LTE cold-start crashes if T exceeds ~5e5 K at those depths.
    T = np.clip(T, 3000.0, 5.0e5)
    return T


def calculate_sigma_with_smoothing(r, v, min_dlnr=1e-5, sigma_min=-0.99, sigma_max=10.0, smooth_window=3):
    """
    Calculate sigma = dlnV/dlnr - 1 with physical constraints and smoothing.
    
    Parameters:
    -----------
    r : array
        Radius in 10^10 cm
    v : array
        Velocity in km/s
    min_dlnr : float
        Minimum dlnr threshold to avoid numerical instability
    sigma_min : float
        Minimum sigma value (physical lower bound)
    sigma_max : float
        Maximum sigma value (physical upper bound)
    smooth_window : int
        Window size for moving average smoothing
    
    Returns:
    --------
    sigma : array
        Sigma values with physical constraints
    """
    sigma = np.zeros_like(r)
    n = len(r)
    
    for i in range(1, n-1):
        dlnv = np.log(v[i+1]) - np.log(v[i-1])
        dlnr = np.log(r[i+1]) - np.log(r[i-1])
        
        if abs(dlnr) > min_dlnr:
            sigma[i] = dlnv / dlnr - 1
        else:
            # Fall back to one-sided differences
            if i > 0:
                dlnv = np.log(v[i]) - np.log(v[i-1])
                dlnr = np.log(r[i]) - np.log(r[i-1])
                if abs(dlnr) > min_dlnr:
                    sigma[i] = dlnv / dlnr - 1
                else:
                    sigma[i] = sigma[i-1]
            else:
                sigma[i] = 0.0
    
    # Boundary values
    sigma[0] = sigma[1]
    sigma[-1] = sigma[-2]
    
    # Apply conservative clipping to prevent numerical extremes
    # CMFGEN requires sigma > -1 to keep density scaling positive (see set_rv_hydro_model_v3.f)
    # sigma = -1 corresponds to constant velocity; values below -1 are numerical artifacts
    sigma_min = -0.99  # Just above -1 to ensure (sigma + 1) > 0
    sigma_max = 20.0   # Captures all realistic supernova velocity laws
    sigma = np.clip(sigma, sigma_min, sigma_max)
    
    # Smooth sigma to avoid discontinuities
    if smooth_window > 1:
        sigma_smoothed = np.copy(sigma)
        half_window = smooth_window // 2
        for i in range(len(sigma)):
            start = max(0, i - half_window)
            end = min(len(sigma), i + half_window + 1)
            sigma_smoothed[i] = np.mean(sigma[start:end])
        sigma = sigma_smoothed
    
    return sigma


def calculate_sigma(r, v, min_dlnr=1e-5):
    """
    Calculate sigma = dlnV/dlnr - 1 with minimum dlnr protection.
    
    Parameters:
    -----------
    r : array
        Radius in 10^10 cm
    v : array
        Velocity in km/s
    min_dlnr : float
        Minimum dlnr threshold to avoid numerical instability
    
    Returns:
    --------
    sigma : array
        Sigma values
    """
    return calculate_sigma_with_smoothing(r, v, min_dlnr=min_dlnr, sigma_min=-0.99, sigma_max=20.0, smooth_window=9)


def smooth_sigma_savgol(r, velocity, window_length=9, polyorder=3, threshold=2.0):
    """
    Apply Savitzky-Golay smoothing to velocity profile and recalculate Sigma.
    
    This function applies Savitzky-Golay filtering to velocity in regions where Sigma
    exceeds a threshold, then recalculates Sigma from the smoothed velocity. This
    prevents CMFGEN segmentation faults by smoothing sharp velocity gradients while
    preserving physical consistency.
    
    Parameters:
    -----------
    r : array
        Radius values (10^10 cm)
    velocity : array
        Velocity values (km/s)
    window_length : int
        Length of the filter window (must be odd, default: 9)
    polyorder : int
        Order of the polynomial used to fit the samples (default: 3)
    threshold : float
        Sigma threshold above which velocity smoothing is applied (default: 2.0)
    
    Returns:
    --------
    sigma_recalculated : array
        Recalculated Sigma values from smoothed velocity
    velocity_smoothed : array
        Smoothed velocity values (km/s)
    """
    # Convert velocity to cm/s for Sigma calculation
    v_cgs = velocity * 1e5  # km/s to cm/s
    r_cgs = r * 1e10  # 10^10 cm to cm
    
    # Calculate original Sigma
    sigma_original = np.zeros_like(r)
    for i in range(1, len(r)-1):
        dlnv = np.log(v_cgs[i+1]) - np.log(v_cgs[i-1])
        dlnr = np.log(r_cgs[i+1]) - np.log(r_cgs[i-1])
        if abs(dlnr) > 1e-10:
            sigma_original[i] = dlnv / dlnr - 1
    sigma_original[0] = sigma_original[1]
    sigma_original[-1] = sigma_original[-2]
    
    velocity_smoothed = v_cgs.copy()
    
    # Ensure window_length is odd
    if window_length % 2 == 0:
        window_length += 1
    
    # Ensure window_length is not larger than array
    if window_length >= len(velocity_smoothed):
        window_length = len(velocity_smoothed) - 1
        if window_length % 2 == 0:
            window_length -= 1
    
    # Identify shock regions where Sigma > threshold
    shock_mask = np.abs(sigma_original) > threshold
    
    if not np.any(shock_mask):
        print("[Sigma Smoothing] No shock regions detected (Sigma < threshold)")
        return sigma_original, velocity
    
    # Find contiguous shock regions
    shock_regions = []
    in_shock = False
    shock_start = 0
    
    for i in range(len(sigma_original)):
        if shock_mask[i] and not in_shock:
            in_shock = True
            shock_start = i
        elif not shock_mask[i] and in_shock:
            in_shock = False
            shock_regions.append((shock_start, i))
    
    if in_shock:
        shock_regions.append((shock_start, len(sigma_original)))
    
    # Apply Savitzky-Golay smoothing to velocity in each shock region with padding
    for start, end in shock_regions:
        # Add padding around shock region for smooth transitions
        pad = window_length // 2
        padded_start = max(0, start - pad)
        padded_end = min(len(velocity_smoothed), end + pad)
        
        # Apply Savitzky-Golay filter to padded region
        region_velocity = velocity_smoothed[padded_start:padded_end]
        if len(region_velocity) >= window_length:
            smoothed_region = savgol_filter(region_velocity, window_length, polyorder)
            velocity_smoothed[padded_start:padded_end] = smoothed_region
    
    # Print statistics
    print(f"[Sigma Smoothing] Applied Savitzky-Golay smoothing to velocity")
    print(f"[Sigma Smoothing] Window length: {window_length}, Polynomial order: {polyorder}")
    print(f"[Sigma Smoothing] Threshold: {threshold}")
    print(f"[Sigma Smoothing] Found {len(shock_regions)} shock regions")
    print(f"[Sigma Smoothing] Original Sigma range: {sigma_original.min():.3f} to {sigma_original.max():.3f}")
    print(f"[Sigma Smoothing] Velocity range before: {v_cgs.min()/1e5:.3f} to {v_cgs.max()/1e5:.3f} km/s")
    print(f"[Sigma Smoothing] Velocity range after:  {velocity_smoothed.min()/1e5:.3f} to {velocity_smoothed.max()/1e5:.3f} km/s")
    
    # Recalculate Sigma from smoothed velocity
    sigma_recalculated = np.zeros_like(r)
    for i in range(1, len(r)-1):
        dlnv = np.log(velocity_smoothed[i+1]) - np.log(velocity_smoothed[i-1])
        dlnr = np.log(r_cgs[i+1]) - np.log(r_cgs[i-1])
        if abs(dlnr) > 1e-10:
            sigma_recalculated[i] = dlnv / dlnr - 1
    sigma_recalculated[0] = sigma_recalculated[1]
    sigma_recalculated[-1] = sigma_recalculated[-2]
    
    print(f"[Sigma Smoothing] Recalculated Sigma range:  {sigma_recalculated.min():.3f} to {sigma_recalculated.max():.3f}")
    
    # Apply additional Savitzky-Golay smoothing to the recalculated Sigma values
    # This further reduces sharp gradients that remain after velocity smoothing
    shock_mask = np.abs(sigma_recalculated) > threshold
    if np.any(shock_mask):
        for start, end in shock_regions:
            pad = window_length // 2
            padded_start = max(0, start - pad)
            padded_end = min(len(sigma_recalculated), end + pad)
            
            region_sigma = sigma_recalculated[padded_start:padded_end]
            if len(region_sigma) >= window_length:
                smoothed_region = savgol_filter(region_sigma, window_length, polyorder)
                sigma_recalculated[padded_start:padded_end] = smoothed_region
        
        print(f"[Sigma Smoothing] Applied additional Sigma smoothing")
        print(f"[Sigma Smoothing] Sigma range after additional smoothing:  {sigma_recalculated.min():.3f} to {sigma_recalculated.max():.3f}")
    
    return sigma_recalculated, velocity_smoothed / 1e5  # Return velocity in km/s
MAX_DRHO_PCT = 15.0  # Max 15% jump in density between adjacent grid points
MAX_DT_PCT = 15.0  # Max 15% jump in temperature between adjacent grid points
MAX_DV_PCT = 10.0  # Max 10% jump in velocity between adjacent grid points
MAX_DR_R_PCT = 5.0  # Max 5% radial grid spacing (dr/r) for radiation transfer stability (increased from 2% per Lechat review)
# NOTE: Increased grid points to prioritize problematic regions (e.g., depths 49, 51-55)
# where NaN values occur in STEQ due to matrix ill-conditioning. These regions correspond
# to sharp ionization fronts that need finer resolution for CMFGEN stability.
MAX_GRID_POINTS = 800  # Increased from 600 for better problematic region coverage
MAX_GRID_POINTS_OUTER = 1200  # Allow more points for outer boundary dtau/dr enforcement
N_BASE = 100  # Number of base points for initial grid
N_ITERATIONS = 20  # Number of refinement iterations (increased to 20 for stricter enforcement)
MIN_DRHO_PCT = 2.0  # Min 2% jump in density for coarsening
MIN_DT_PCT = 1.0  # Min 1% jump in temperature for coarsening
MIN_DV_PCT = 1.0  # Min 1% jump in velocity for coarsening
MIN_DR_FRACTION = 0.001  # Minimum dr as fraction of local radius (0.1% per Lechat review)
MIN_DR_R_THIN = 0.005    # Aggressive thinning threshold: remove points with dr/r < 0.5% before writing

def thin_degenerate_grid(radius, *arrays):
    """Remove near-degenerate grid points to prevent singular BA matrices in CMFGEN.

    CMFGEN's CMF_BLKBAND_V3 solver becomes singular when adjacent depth points have
    dr/r < ~0.5%, causing IFAIL errors and SIGSEGV crashes.  This function greedily
    walks the grid and merges consecutive points that are too close together: a point
    is kept only if the accumulated radial displacement since the last kept point is
    >= MIN_DR_R_THIN of the local radius.  The first and last points are always kept.

    Args:
        radius: 1-D array of radii in ANY consistent unit (outer->inner = decreasing).
        *arrays: Additional arrays (velocity, density, temperature, ...) filtered in sync.

    Returns:
        Tuple (radius_thinned, *arrays_thinned).
    """
    n = len(radius)
    keep = np.zeros(n, dtype=bool)
    keep[0] = True   # always keep outer boundary
    keep[n-1] = True # always keep inner boundary

    last_kept = 0
    for i in range(1, n - 1):
        dr = abs(radius[last_kept] - radius[i])
        r_ref = max(abs(radius[last_kept]), abs(radius[i]), 1e-30)
        if dr / r_ref >= MIN_DR_R_THIN:
            keep[i] = True
            last_kept = i

    n_removed = int(np.sum(~keep))
    n_kept = int(np.sum(keep))
    if n_removed > 0:
        print(f"[Grid Thinning] Removed {n_removed} near-degenerate points (dr/r < {MIN_DR_R_THIN*100:.1f}%)")
        print(f"[Grid Thinning] Grid size: {n} -> {n_kept} points")
        radius = radius[keep]
        arrays = tuple(arr[keep] for arr in arrays)
    else:
        print(f"[Grid Thinning] No degenerate points found (all dr/r >= {MIN_DR_R_THIN*100:.1f}%)")
    return (radius,) + arrays


# Optical depth-based grid refinement parameters
MAX_DTAU_PCT = 30.0  # Maximum percentage change in optical depth
MIN_DTAU = 0.01  # Minimum optical depth change between points
DTAU_OUTER_THRESHOLD = 0.1  # Target dtau for outer boundary points


def calculate_optical_depth(r, rho, kappa):
    """Calculate optical depth from outer boundary inward.
    
    Parameters:
    -----------
    r : array (outer->inner order)
    rho : array (outer->inner order)
    kappa : array (outer->inner order)
    
    Returns:
    --------
    tau : array (outer->inner order, tau[0]=0 at outer boundary)
    """
    tau = np.zeros_like(r)
    tau[0] = 0  # Outer boundary
    for i in range(1, len(r)):
        dr = abs(r[i-1] - r[i])
        kappa_avg = 0.5 * (kappa[i-1] + kappa[i])
        rho_avg = 0.5 * (rho[i-1] + rho[i])
        tau[i] = tau[i-1] + kappa_avg * rho_avg * dr
    return tau


def apply_targeted_temperature_gradient_fix(temperature, dtdr_threshold=1e-3, min_change_pct=1e-3, temp_min=3000.0):
    """Apply targeted temperature gradient fix to prevent CMFGEN instability.
    
    Contained-region algorithm: Identifies flat regions (consecutive cells where dT/T < threshold)
    and applies minimal cooling within each region only. Cooling cascade is contained — it stops
    at the region boundary. All HERACLES temperature structure outside flat regions is untouched.
    
    Args:
        temperature: Temperature array in K (outer->inner order, index 0 = outermost)
        dtdr_threshold: dT/T threshold in percent (default: 1e-3%)
        min_change_pct: Minimum change to apply when threshold exceeded (default: 1e-3%)
        temp_min: Minimum temperature floor (default: 3000K). System may go slightly below
                  this if needed to avoid flat regions and maintain gradients.
    
    Returns:
        temperature_fixed: Temperature array with targeted gradient fix applied
    """
    temperature_fixed = temperature.copy()
    adjustments_made = 0
    n = len(temperature)
    
    # Step 1: Identify flat regions (consecutive cells where dT/T < dtdr_threshold)
    flat_regions = []
    in_flat = False
    flat_start = 0
    
    for i in range(n - 1):
        t_outer = temperature[i]      # Outer cell (index i)
        t_inner = temperature[i+1]    # Inner cell (index i+1)
        
        if t_inner > 1e-6:
            dtdr_pct = abs(t_outer - t_inner) / t_inner * 100
        else:
            dtdr_pct = 0
        
        if dtdr_pct < dtdr_threshold and not in_flat:
            in_flat = True
            flat_start = i
        elif dtdr_pct >= dtdr_threshold and in_flat:
            in_flat = False
            flat_regions.append((flat_start, i))
    
    if in_flat:
        flat_regions.append((flat_start, n-1))
    
    # Step 2: Fix each flat region — inner to outer, contained within region
    for start, end in flat_regions:
        # Process from inner (end) to outer (start)
        # The innermost cell of the region (end) stays at its original temperature
        # Each successive outer cell is cooled relative to its inner neighbor
        # No re-checking: the region was identified as flat, so all outer cells need cooling
        for i in range(end - 1, start - 1, -1):
            t_inner = temperature_fixed[i+1]  # Already-fixed inner neighbor
            t_new = t_inner * (1.0 - min_change_pct / 100.0)
            temperature_fixed[i] = t_new
            adjustments_made += 1
    
    if adjustments_made > 0:
        print(f"[Targeted Temperature Gradient Fix] Applied {adjustments_made} adjustments")
        print(f"[Targeted Temperature Gradient Fix] Found {len(flat_regions)} flat regions")
        print(f"[Targeted Temperature Gradient Fix] dT/dr threshold: {dtdr_threshold:.4f}%")
        print(f"[Targeted Temperature Gradient Fix] Minimum change: {min_change_pct:.4f}%")
        print(f"[Targeted Temperature Gradient Fix] Temperature range: {temperature_fixed.min():.0f}K to {temperature_fixed.max():.0f}K")
    else:
        print(f"[Targeted Temperature Gradient Fix] No adjustments needed (all gradients above {dtdr_threshold:.4f}%)")
    
    return temperature_fixed


def break_isothermal_regions(temperature, t_floor=3000.0, min_rel_change=5e-4):
    """
    Break flat/isothermal temperature regions by adding a tiny monotonic gradient.

    SEDONA's NLTE solver can crash with a singular matrix when many consecutive
    zones have exactly the same temperature (especially at low densities). This
    function identifies contiguous blocks of cells at the floor temperature and
    applies a small linear perturbation so each cell differs from its neighbour,
    keeping the rate matrix well-conditioned.

    Parameters:
    -----------
    temperature : array (K), outer->inner order
    t_floor : float
        Temperature floor value to identify flat regions (default: 3000.0)
    min_rel_change : float
        Minimum fractional dT/T to enforce (default: 5e-4 = 0.05%)

    Returns:
    --------
    temperature : array with flat regions broken
    """
    temperature = np.asarray(temperature, dtype=float).copy()
    n = len(temperature)
    if n < 2:
        return temperature

    tol = 1e-9
    at_floor = np.abs(temperature - t_floor) < tol
    if not np.any(at_floor):
        return temperature

    n_broken = 0
    i = 0
    while i < n:
        if not at_floor[i]:
            i += 1
            continue
        # Find contiguous floored block
        j = i
        while j < n and at_floor[j]:
            j += 1
        block_len = j - i
        if block_len >= 2:
            # Apply tiny perturbation: outermost cell gets t_floor*(1+delta),
            # next gets t_floor*(1+2*delta), etc., so every pair has
            # |dT/T| = delta > min_rel_change.  This is a ~0.05% effect.
            delta = max(min_rel_change * 1.1, 1e-6)
            for k in range(i, j):
                temperature[k] = t_floor * (1.0 + (k - i + 1) * delta)
            n_broken += block_len
        i = j

    if n_broken > 0:
        print(f"[Anti-isothermal] Broken {n_broken} floored cells "
              f"(T == {t_floor:.1f} K) with {min_rel_change:.2e} gradient")
    return temperature


def calculate_dr_over_r_log(r1, r2):
    """Calculate dr/r in log-space (consistent with logarithmic grids).
    
    Parameters:
    -----------
    r1, r2 : float
        Radii in cm (can be in any order)
    
    Returns:
    --------
    dr_r_pct : float
        dr/r as percentage calculated in log-space (always positive)
    """
    dlogr = abs(np.log10(r2) - np.log10(r1))  # Take absolute value
    return (10**dlogr - 1) * 100  # Convert to %


def calculate_electron_density(density, temperature, mass_fractions):
    """
    Calculate electron density (cm^-3) accounting for ionization state.
    
    Args:
        density: Density in g/cm³ (array)
        temperature: Temperature in K (array)
        mass_fractions: Dict of mass fractions (e.g., {'HYD': X_H, 'HE': X_He, ...})
    
    Returns:
        n_e: Electron density in cm^-3 (array)
    """
    N_A = 6.02214076e23  # Avogadro's number (mol^-1)
    
    # Atomic data: (atomic mass A, atomic number Z)
    # Map CMFGEN element names to atomic data
    atomic_data = {
        'HYD': (1.0,  1),   # Hydrogen
        'HE':  (4.0,  2),   # Helium
        'CARB': (12.0, 6),  # Carbon
        'NITR': (14.0, 7),  # Nitrogen
        'OXY': (16.0, 8),   # Oxygen
        'NEON': (20.0, 10), # Neon
        'MAGN': (24.0, 12), # Magnesium
        'SIL': (28.0, 14),  # Silicon
        'SULF': (32.0, 16), # Sulfur
        'IRON': (56.0, 26), # Iron
        'SI':  (28.0, 14),  # Silicon (alternate name)
        'FE':  (56.0, 26),  # Iron (alternate name)
        'C':   (12.0, 6),   # Carbon (alternate name)
        'O':   (16.0, 8),   # Oxygen (alternate name)
    }
    
    n_e = np.zeros_like(density)
    
    for elem, (A, Z) in atomic_data.items():
        if elem in mass_fractions:
            X = mass_fractions[elem]  # Mass fraction
            if np.any(X > 0):
                # Mean ionization state ⟨Z⟩ for this element
                if elem in ['HYD', 'H']:  # Hydrogen
                    # Fit to Saha equation for H I <-> H II
                    # ⟨Z⟩_H ≈ 1 / (1 + exp(157800/T - 10))
                    Z_mean = 1.0 / (1.0 + np.exp(157800.0 / temperature - 10.0))
                elif elem in ['HE', 'HELIUM']:  # Helium
                    # Fit to Saha equation for He I <-> He II and He II <-> He III
                    # ⟨Z⟩_He ≈ 2 / (1 + exp(285000/T - 15))
                    Z_mean = 2.0 / (1.0 + np.exp(285000.0 / temperature - 15.0))
                else:  # Metals (C, O, Ne, Mg, Si, S, Fe)
                    # Assume fully ionized for T > 10,000 K (typical for SNe)
                    # For lower T, use a smooth transition
                    T_transition = 20000.0  # K (adjust based on element)
                    Z_mean = Z * (1.0 / (1.0 + np.exp(-(temperature - T_transition) / 5000.0)))
                
                # Contribution to electron density: ρ * X * ⟨Z⟩ / A * N_A
                n_e += density * X * Z_mean / A * N_A
    
    return n_e


def approximate_rosseland_opacity(T, rho, X, Y, Z, n_e=None):
    """
    Approximate Rosseland mean opacity (cm²/g) for supernova ejecta.
    
    Parameters:
    -----------
    T : array
        Temperature in K
    rho : array
        Density in g/cm³
    X : array
        Hydrogen mass fraction
    Y : array
        Helium mass fraction
    Z : array
        Metal mass fraction
    n_e : array, optional
        Electron density in cm^-3 (not used in current implementation)
    
    Returns:
    --------
    kappa_rosseland : array
        Rosseland mean opacity in cm²/g
    """
    # Avoid division by zero and prevent blow-up at low T
    T_floor = 1000.0  # Minimum temperature for opacity calculation (K)
    T_safe = np.maximum(T, T_floor)
    
    # Electron scattering with composition dependence
    kappa_es = 0.2 * (1.0 + X) / (1.0 + 3.0*X + 0.25*Y)
    
    # Bound-free (photoionization)
    Z_over_A = 0.5  # Approximate for solar-like metals
    kappa_bf = 4.34e25 * Z_over_A * Z * rho * (T_safe ** -3.5)
    
    # Free-free (bremsstrahlung) - use composition-dependent formula
    kappa_ff = 1.7e20 * (1.0 + X) * rho * (T_safe ** -3.5)
    
    # Total Rosseland opacity
    kappa_rosseland = kappa_es + kappa_bf + kappa_ff
    
    return kappa_rosseland


def refine_grid_forward(radius_cm, vel_cms, rho_cgs, temp_104k, ed_cgs=None, output_dir=None):
    """
    Simplified robust grid refinement using logarithmic spacing with adaptive refinement.
    
    Algorithm:
    1. Start with logarithmic spacing from inner to outer boundary
    2. Iteratively refine regions with large gradients
    3. Add explicit points at outer boundary for CMFGEN stability
    
    This replaces the complex binary search approach that was getting stuck.
    """
    # Convert to inner->outer order
    r_io = radius_cm[::-1]
    v_io = np.abs(vel_cms[::-1])
    rho_io = rho_cgs[::-1]
    t_io = temp_104k[::-1]
    ne_io = ed_cgs[::-1] if ed_cgs is not None else np.zeros_like(radius_cm)
    tau_io = np.zeros_like(radius_cm)  # Placeholder
    
    # Start with logarithmic grid
    log_r_min = np.log10(r_io[0])
    log_r_max = np.log10(r_io[-1])
    
    # Initial grid: 100 points logarithmically spaced
    n_base = 100
    log_r_grid = np.linspace(log_r_min, log_r_max, n_base)
    r_grid = 10**log_r_grid
    
    # Interpolate values onto logarithmic grid
    log_r_io = np.log10(r_io)
    v_grid = 10**np.interp(log_r_grid, log_r_io, np.log10(v_io))
    rho_grid = 10**np.interp(log_r_grid, log_r_io, np.log10(rho_io))
    t_grid = 10**np.interp(log_r_grid, log_r_io, np.log10(t_io))
    ne_grid = 10**np.interp(log_r_grid, log_r_io, np.log10(np.maximum(ne_io, 1e-10)))
    tau_grid = np.zeros_like(r_grid)  # Placeholder
    
    # Iterative refinement: add points at large gradients
    for iteration in range(N_ITERATIONS):
        if len(r_grid) >= MAX_GRID_POINTS:
            break
        
        # Calculate gradients
        dr = np.diff(r_grid)
        dv = np.abs(np.diff(v_grid))
        drho = np.abs(np.diff(rho_grid))
        dt = np.abs(np.diff(t_grid))
        
        # Calculate percentage jumps
        dv_pct = dv / v_grid[:-1] * 100
        drho_pct = drho / rho_grid[:-1] * 100
        dt_pct = dt / t_grid[:-1] * 100
        
        # Find violations
        violations = (dv_pct > MAX_DV_PCT) | (drho_pct > MAX_DRHO_PCT) | (dt_pct > MAX_DT_PCT)
        
        if not np.any(violations):
            break
        
        # Add points at violation locations (geometric mean)
        insert_indices = np.where(violations)[0]
        
        for idx in insert_indices:
            if len(r_grid) >= MAX_GRID_POINTS:
                break
            
            # Insert geometric mean
            r_new = np.sqrt(r_grid[idx] * r_grid[idx + 1])
            log_r_new = np.log10(r_new)
            
            # Interpolate values
            v_new = 10**np.interp(log_r_new, np.log10(r_grid), np.log10(v_grid))
            rho_new = 10**np.interp(log_r_new, np.log10(r_grid), np.log10(rho_grid))
            t_new = 10**np.interp(log_r_new, np.log10(r_grid), np.log10(t_grid))
            ne_new = 10**np.interp(log_r_new, np.log10(r_grid), np.log10(np.maximum(ne_grid, 1e-10)))
            tau_new = 0.0  # Placeholder
            
            # Insert
            r_grid = np.insert(r_grid, idx + 1, r_new)
            v_grid = np.insert(v_grid, idx + 1, v_new)
            rho_grid = np.insert(rho_grid, idx + 1, rho_new)
            t_grid = np.insert(t_grid, idx + 1, t_new)
            ne_grid = np.insert(ne_grid, idx + 1, ne_new)
            tau_grid = np.insert(tau_grid, idx + 1, tau_new)
    
    print(f"[Grid Refinement] After adaptive refinement: {len(r_grid)} points")
    
    # Add explicit points at outer boundary for CMFGEN radiation transport
    outer_boundary_spacing_pct = 2.0
    n_outer_points = 20
    
    r_outer = r_grid[-1]
    r_inner_outer = r_grid[-2] if len(r_grid) > 1 else r_grid[-1] * 0.9
    
    if len(r_grid) < MAX_GRID_POINTS - n_outer_points:
        # Calculate expansion time from outer boundary
        t_expansion = r_outer / v_grid[-1] if v_grid[-1] > 0 else 400 * 86400
        rho_norm = rho_grid[-1] * r_outer**2
        T_outer = t_grid[-1]
        
        # Collect outer boundary points
        outer_points_r = []
        outer_points_v = []
        outer_points_rho = []
        outer_points_t = []
        
        for i in range(1, n_outer_points + 1):
            spacing_fraction = outer_boundary_spacing_pct / 100.0
            r_new = r_outer * (1 - i * spacing_fraction)
            
            if r_new > r_inner_outer and r_new < r_outer:
                v_new = r_new / t_expansion
                rho_new = rho_norm / r_new**2
                t_new = T_outer
                
                outer_points_r.append(r_new)
                outer_points_v.append(v_new)
                outer_points_rho.append(rho_new)
                outer_points_t.append(t_new)
        
        # Sort and insert
        if outer_points_r:
            sort_indices = np.argsort(outer_points_r)
            outer_points_r = [outer_points_r[i] for i in sort_indices]
            outer_points_v = [outer_points_v[i] for i in sort_indices]
            outer_points_rho = [outer_points_rho[i] for i in sort_indices]
            outer_points_t = [outer_points_t[i] for i in sort_indices]
            
            for r_new, v_new, rho_new, t_new in zip(outer_points_r, outer_points_v, outer_points_rho, outer_points_t):
                r_grid = np.insert(r_grid, -1, r_new)
                v_grid = np.insert(v_grid, -1, v_new)
                rho_grid = np.insert(rho_grid, -1, rho_new)
                t_grid = np.insert(t_grid, -1, t_new)
            
            print(f"[Grid Refinement] Added {len(outer_points_r)} points at outer boundary")
    
    print(f"[Grid Refinement] Final grid: {len(r_grid)} points")
    
    # Generate diagnostic plot
    if output_dir:
        try:
            # Create point_iterations array (all zeros since we don't track iterations in this method)
            point_iterations = np.zeros(len(r_grid), dtype=int)
            
            plot_grid_refinement_diagnostics(
                r_io, v_io, rho_io, t_io,
                r_grid, v_grid, rho_grid, t_grid,
                point_iterations,
                output_dir  # Pass directory, not full path
            )
        except ImportError:
            print("[Grid Refinement] matplotlib not available, skipping diagnostic plots")
        except Exception as e:
            print(f"[Grid Refinement] WARNING: Failed to generate diagnostic plot: {e}")
    
    # Return in inner->outer order
    return r_grid, v_grid, rho_grid, t_grid, ne_grid


def detect_gradient_violations(radius_cm, rho, v, T, max_drho_pct=MAX_DRHO_PCT, max_dt_pct=MAX_DT_PCT, max_dv_pct=MAX_DV_PCT):
    """
    Detect regions where gradient constraints are violated.

    Args:
        radius_cm: Radius array (cm, outer->inner order).
        rho: Density array (g/cm³).
        v: Velocity array (cm/s).
        T: Temperature array (K).
        max_drho_pct: Max allowed % jump in density.
        max_dt_pct: Max allowed % jump in temperature.
        max_dv_pct: Max allowed % jump in velocity.

    Returns:
        violation_mask: Boolean array (True = violation at this cell boundary).
    """
    n = len(radius_cm)
    violation_mask = np.zeros(n - 1, dtype=bool)  # True if violation between i and i+1

    for i in range(n - 1):
        # Density jump
        if rho[i] > 1e-30:
            drho_pct = abs(rho[i+1] - rho[i]) / rho[i] * 100
        else:
            drho_pct = 0

        # Temperature jump
        if T[i] > 1e-6:
            dt_pct = abs(T[i+1] - T[i]) / T[i] * 100
        else:
            dt_pct = 0

        # Velocity jump
        if v[i] > 1e-10:
            dv_pct = abs(v[i+1] - v[i]) / v[i] * 100
        else:
            dv_pct = 0

        # Check violations
        if (drho_pct > max_drho_pct or
            dt_pct > max_dt_pct or
            dv_pct > max_dv_pct):
            violation_mask[i] = True

    return violation_mask


def smooth_violations(radius_cm, rho, v, T, violation_mask, window_size=3):
    """
    Strictly enforce gradient limits by monotonically blending across violating intervals.
    For velocity, apply Gaussian smoothing in the shock region to avoid artificial kinks.

    Args:
        radius_cm: Radius array (cm, outer->inner order).
        rho: Density array (g/cm³).
        v: Velocity array (cm/s).
        T: Temperature array (K).
        violation_mask: Boolean array (True = violation at boundary i).
        window_size: Unused (kept for API compatibility); blending is log-space monotone.

    Returns:
        rho_smoothed, v_smoothed, T_smoothed: Arrays with gradients strictly within limits.
    """
    n = len(rho)
    rho_smoothed = rho.copy()
    v_smoothed = v.copy()
    T_smoothed = T.copy()

    # --- Pass 1: find contiguous violating regions and re-interpolate monotonically ---
    # Group consecutive violating boundaries into spans [i_start, i_end] (boundary indices).
    # For each span, anchor the endpoints and redistribute intermediate values in log-space
    # so that every step stays within the per-boundary budget.
    violating_boundaries = np.where(violation_mask)[0]
    if len(violating_boundaries) == 0:
        return rho_smoothed, v_smoothed, T_smoothed

    # Build contiguous groups
    groups = []
    group = [violating_boundaries[0]]
    for b in violating_boundaries[1:]:
        if b == group[-1] + 1:
            group.append(b)
        else:
            groups.append(group)
            group = [b]
    groups.append(group)

    for group in groups:
        # Cell indices that need resampling: from group[0] to group[-1]+1 inclusive
        i0 = group[0]       # first violating boundary → between cell i0 and i0+1
        i1 = group[-1] + 1  # last cell in span
        # Widen by 1 on each side to give smooth transitions at the edges
        i0 = max(0, i0 - 1)
        i1 = min(n - 1, i1 + 1)
        span = i1 - i0 + 1  # number of cells in the span
        if span < 2:
            continue

        # Anchor values at both ends; redistribute log-linearly inside span
        for arr, key in [(rho_smoothed, 'rho'), (T_smoothed, 'T'), (v_smoothed, 'v')]:
            v0, v1 = arr[i0], arr[i1]
            if v0 <= 0 or v1 <= 0:
                # Fall back to linear for non-positive arrays
                arr[i0:i1+1] = np.linspace(v0, v1, span)
            else:
                arr[i0:i1+1] = np.exp(np.linspace(np.log(v0), np.log(v1), span))

    # --- Pass 2: Gaussian velocity smoothing in regions with >5% velocity jump ---
    # This blurs the velocity profile slightly as permitted, reducing sharp kinks
    # that arise from shock-front grid resolution.
    from scipy.ndimage import gaussian_filter1d
    dv_pct = np.zeros(n - 1)
    for i in range(n - 1):
        if v_smoothed[i] > 1e-10:
            dv_pct[i] = abs(v_smoothed[i+1] - v_smoothed[i]) / v_smoothed[i] * 100.0
    shock_mask = dv_pct > 5.0  # flag boundaries with >5% velocity jump
    if np.any(shock_mask):
        # Build a per-cell weight: 1 where in shock region, 0 elsewhere
        cell_weight = np.zeros(n)
        for i in np.where(shock_mask)[0]:
            cell_weight[max(0, i-2):min(n, i+3)] = 1.0
        v_gauss = gaussian_filter1d(v_smoothed, sigma=2.0)
        v_smoothed = np.where(cell_weight > 0, v_gauss, v_smoothed)
        # Ensure velocity stays positive
        v_smoothed = np.maximum(v_smoothed, 1e2)  # floor at 1 cm/s

    return rho_smoothed, v_smoothed, T_smoothed


def plot_smoothing_evolution(smoothing_history, output_dir, r_input=None, rho_input=None, v_input=None, T_input=None):
    """
    Generate diagnostic plots showing evolution of density, temperature, and velocity during smoothing.

    Args:
        smoothing_history: Dictionary containing iteration data (r, rho, v, T, violations)
        output_dir: Directory to save plots
        r_input: Original input radius array (optional, for comparison)
        rho_input: Original input density array (optional, for comparison)
        v_input: Original input velocity array (optional, for comparison)
        T_input: Original input temperature array (optional, for comparison)
    """
    import matplotlib.pyplot as plt

    n_iterations = len(smoothing_history['iteration'])
    if n_iterations < 2:
        print("[plot_smoothing_evolution] No smoothing iterations performed, skipping plots")
        return

    # Create figure with 3 subplots (one for each quantity) + gradient jump panel
    fig, axes = plt.subplots(3, 2, figsize=(14, 14))
    fig.suptitle(f'Gradient Smoothing Evolution ({n_iterations - 1} iterations)', fontsize=14)

    colors = plt.cm.viridis(np.linspace(0, 1, n_iterations))

    # Plot density evolution
    ax = axes[0, 0]
    # Plot input if available
    if r_input is not None and rho_input is not None:
        ax.loglog(r_input / 1e10, rho_input, color='black', linestyle='--', linewidth=2,
                 alpha=0.8, label='Input (pre-refinement)')
    for i in range(n_iterations):
        r = smoothing_history['r'][i]
        rho = smoothing_history['rho'][i]
        ax.loglog(r / 1e10, rho, color=colors[i], alpha=0.7, linewidth=1.5,
                 label=f'Iter {smoothing_history["iteration"][i]}')
    ax.set_ylabel('Density (g/cm³)', fontsize=12)
    ax.set_xlabel('Radius (10$^{10}$ cm)', fontsize=12)
    ax.set_title('Density Evolution', fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.legend(loc='best', fontsize=7)

    # Plot density gradient jumps
    ax = axes[0, 1]
    # Plot input if available
    if r_input is not None and rho_input is not None:
        drho_pct_input = np.zeros(len(r_input) - 1)
        for j in range(len(r_input) - 1):
            if rho_input[j] > 1e-30:
                drho_pct_input[j] = abs(rho_input[j+1] - rho_input[j]) / rho_input[j] * 100
        ax.semilogx(r_input[:-1] / 1e10, drho_pct_input, color='black', linestyle='--', linewidth=2,
               alpha=0.8, label='Input (pre-refinement)')
    for i in range(n_iterations):
        r = smoothing_history['r'][i]
        rho = smoothing_history['rho'][i]
        # Calculate percentage jumps
        drho_pct = np.zeros(len(r) - 1)
        for j in range(len(r) - 1):
            if rho[j] > 1e-30:
                drho_pct[j] = abs(rho[j+1] - rho[j]) / rho[j] * 100
        ax.semilogx(r[:-1] / 1e10, drho_pct, color=colors[i], alpha=0.7, linewidth=1.5,
               label=f'Iter {smoothing_history["iteration"][i]}')
    ax.axhline(y=15, color='red', linestyle='--', alpha=0.5, label='MAX_DRHO_PCT')
    ax.set_ylabel('Density Jump (%)', fontsize=12)
    ax.set_xlabel('Radius (10$^{10}$ cm)', fontsize=12)
    ax.set_title('Density Gradient Jumps', fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.legend(loc='best', fontsize=7)

    # Plot temperature evolution
    ax = axes[1, 0]
    # Plot input if available
    if r_input is not None and T_input is not None:
        ax.loglog(r_input / 1e10, T_input, color='black', linestyle='--', linewidth=2,
                 alpha=0.8, label='Input (pre-refinement)')
    for i in range(n_iterations):
        r = smoothing_history['r'][i]
        T = smoothing_history['T'][i]
        ax.loglog(r / 1e10, T, color=colors[i], alpha=0.7, linewidth=1.5,
                 label=f'Iter {smoothing_history["iteration"][i]}')
    ax.set_ylabel('Temperature (K)', fontsize=12)
    ax.set_xlabel('Radius (10$^{10}$ cm)', fontsize=12)
    ax.set_title('Temperature Evolution', fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.legend(loc='best', fontsize=7)

    # Plot temperature gradient jumps
    ax = axes[1, 1]
    # Plot input if available
    if r_input is not None and T_input is not None:
        dt_pct_input = np.zeros(len(r_input) - 1)
        for j in range(len(r_input) - 1):
            if T_input[j] > 1e-6:
                dt_pct_input[j] = abs(T_input[j+1] - T_input[j]) / T_input[j] * 100
        ax.semilogx(r_input[:-1] / 1e10, dt_pct_input, color='black', linestyle='--', linewidth=2,
               alpha=0.8, label='Input (pre-refinement)')
    for i in range(n_iterations):
        r = smoothing_history['r'][i]
        T = smoothing_history['T'][i]
        # Calculate percentage jumps
        dt_pct = np.zeros(len(r) - 1)
        for j in range(len(r) - 1):
            if T[j] > 1e-6:
                dt_pct[j] = abs(T[j+1] - T[j]) / T[j] * 100
        ax.semilogx(r[:-1] / 1e10, dt_pct, color=colors[i], alpha=0.7, linewidth=1.5,
               label=f'Iter {smoothing_history["iteration"][i]}')
    ax.axhline(y=20, color='red', linestyle='--', alpha=0.5, label='MAX_DT_PCT')
    ax.set_ylabel('Temperature Jump (%)', fontsize=12)
    ax.set_xlabel('Radius (10$^{10}$ cm)', fontsize=12)
    ax.set_title('Temperature Gradient Jumps', fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.legend(loc='best', fontsize=7)

    # Plot velocity evolution
    ax = axes[2, 0]
    # Plot input if available
    if r_input is not None and v_input is not None:
        ax.loglog(r_input / 1e10, v_input, color='black', linestyle='--', linewidth=2,
                 alpha=0.8, label='Input (pre-refinement)')
    for i in range(n_iterations):
        r = smoothing_history['r'][i]
        v = smoothing_history['v'][i]
        ax.loglog(r / 1e10, v, color=colors[i], alpha=0.7, linewidth=1.5,
                 label=f'Iter {smoothing_history["iteration"][i]}')
    ax.set_ylabel('Velocity (cm/s)', fontsize=12)
    ax.set_xlabel('Radius (10$^{10}$ cm)', fontsize=12)
    ax.set_title('Velocity Evolution', fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.legend(loc='best', fontsize=7)

    # Plot velocity gradient jumps
    ax = axes[2, 1]
    # Plot input if available
    if r_input is not None and v_input is not None:
        dv_pct_input = np.zeros(len(r_input) - 1)
        for j in range(len(r_input) - 1):
            if v_input[j] > 1e-10:
                dv_pct_input[j] = abs(v_input[j+1] - v_input[j]) / v_input[j] * 100
        ax.semilogx(r_input[:-1] / 1e10, dv_pct_input, color='black', linestyle='--', linewidth=2,
               alpha=0.8, label='Input (pre-refinement)')
    for i in range(n_iterations):
        r = smoothing_history['r'][i]
        v = smoothing_history['v'][i]
        # Calculate percentage jumps
        dv_pct = np.zeros(len(r) - 1)
        for j in range(len(r) - 1):
            if v[j] > 1e-10:
                dv_pct[j] = abs(v[j+1] - v[j]) / v[j] * 100
        ax.semilogx(r[:-1] / 1e10, dv_pct, color=colors[i], alpha=0.7, linewidth=1.5,
               label=f'Iter {smoothing_history["iteration"][i]}')
    ax.axhline(y=8, color='red', linestyle='--', alpha=0.5, label='MAX_DV_PCT')
    ax.set_ylabel('Velocity Jump (%)', fontsize=12)
    ax.set_xlabel('Radius (10$^{10}$ cm)', fontsize=12)
    ax.set_title('Velocity Gradient Jumps', fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.legend(loc='best', fontsize=7)

    plt.tight_layout()

    # Save plot
    output_file = output_dir / 'SN_HYDRO_DATA_gradient_smoothing.png'
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"[plot_smoothing_evolution] Saved diagnostic plot to {output_file}")

    # Also plot violation count vs iteration
    if smoothing_history['violations']:
        fig, ax = plt.subplots(figsize=(8, 5))
        iterations = range(1, len(smoothing_history['violations']) + 1)
        violations = smoothing_history['violations']
        ax.plot(iterations, violations, 'o-', linewidth=2, markersize=8, color='red')
        ax.set_xlabel('Smoothing Iteration', fontsize=12)
        ax.set_ylabel('Number of Violations', fontsize=12)
        ax.set_title('Gradient Violations vs Smoothing Iteration', fontsize=12)
        ax.grid(True, alpha=0.3)
        ax.set_xticks(iterations)

        output_file2 = output_dir / 'SN_HYDRO_DATA_gradient_violations_history.png'
        plt.savefig(output_file2, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"[plot_smoothing_evolution] Saved violations plot to {output_file2}")


def ensure_monotonic_grid(r, *arrays, tolerance=1e-10):
    """
    Ensure radius grid is strictly monotonic by removing duplicate or near-duplicate values.
    CMFGEN requires strictly monotonic radius grid (no duplicate values).

    Parameters:
    -----------
    r : array (outer->inner order or inner->outer order)
    *arrays : variable number of arrays corresponding to r (v, rho, t, ne, sigma, kappa, etc.)
    tolerance : float, tolerance for considering radii as duplicate (relative to r value)

    Returns:
    --------
    tuple: (r_monotonic, *arrays_monotonic)
    """
    r = np.asarray(r)
    arrays = [np.asarray(arr) if arr is not None else np.zeros_like(r) for arr in arrays]

    # Find indices where radius is not strictly decreasing (for outer->inner order)
    # or strictly increasing (for inner->outer order)
    is_outer_to_inner = r[0] > r[-1]

    if is_outer_to_inner:
        # Outer->inner: should be strictly decreasing
        mask = np.diff(r) < -tolerance * np.abs(r[:-1])
    else:
        # Inner->outer: should be strictly increasing
        mask = np.diff(r) > tolerance * np.abs(r[:-1])

    # Keep first point, then keep points where the difference is significant
    keep_indices = [0]
    for i in range(1, len(r)):
        if mask[i-1]:
            keep_indices.append(i)
        else:
            # Skip this duplicate/near-duplicate point
            print(f"[ensure_monotonic_grid] Removing duplicate radius at index {i}: r={r[i]:.3e}")

    if len(keep_indices) < len(r):
        print(f"[ensure_monotonic_grid] Removed {len(r) - len(keep_indices)} duplicate radius values")
        print(f"[ensure_monotonic_grid] Grid points: {len(r)} -> {len(keep_indices)}")

    r_monotonic = r[keep_indices]
    arrays_monotonic = [arr[keep_indices] for arr in arrays]

    return r_monotonic, *arrays_monotonic


def refine_grid(radius_cm, vel_cms, rho_cgs, temp_104k, ed_cgs=None, kappa_cgs=None, output_dir=None):
    """
    Optimized adaptive grid refinement: minimal points with strict gradient enforcement.
    Algorithm:
    1. Start with 50 logarithmic base points (reduced from 100).
    2. For up to 10 iterations, add points at geometric mean where gradients exceed thresholds.
    3. Coarsen grid by removing points in low-gradient regions (if removal doesn't violate thresholds).
    4. Add outer boundary points based on optical depth spacing.

    Parameters:
    -----------
    radius_cm : array (outer->inner order, CMFGEN format)
    vel_cms : array (outer->inner)
    rho_cgs : array (outer->inner)
    temp_104k : array (outer->inner)
    ed_cgs : array (optional, outer->inner)
    kappa_cgs : array (optional, outer->inner) - Opacity in cm^2/g
    output_dir : Path (optional)

    Returns:
    --------
    tuple: (r_refined, v_refined, rho_refined, t_refined, ne_refined) in outer->inner order.
    """
    # Convert to outer->inner order for consistency
    r_oi = radius_cm[::-1]  # outer->inner
    v_oi = np.abs(vel_cms[::-1])
    rho_oi = rho_cgs[::-1]
    t_oi = temp_104k[::-1]
    ne_oi = ed_cgs[::-1] if ed_cgs is not None else np.zeros_like(radius_cm)
    kappa_oi = kappa_cgs[::-1] if kappa_cgs is not None else np.full_like(radius_cm, 0.34)

    # --- Step 1: Start with a minimal base grid (50 points) ---
    n_base = 50  # Reduced from 100 to minimize initial points
    log_r_oi = np.log10(r_oi)
    log_r_base = np.linspace(log_r_oi[0], log_r_oi[-1], n_base)
    r_refined = 10 ** log_r_base

    # Interpolate to base grid (log-log space for smoothness)
    v_refined = 10 ** np.interp(log_r_base, log_r_oi, np.log10(v_oi))
    rho_refined = 10 ** np.interp(log_r_base, log_r_oi, np.log10(rho_oi))
    t_refined = 10 ** np.interp(log_r_base, log_r_oi, np.log10(t_oi))
    ne_refined = 10 ** np.interp(log_r_base, log_r_oi, np.log10(np.maximum(ne_oi, 1e-10)))
    kappa_refined = np.interp(log_r_base, log_r_oi, kappa_oi)

    # Calculate optical depth from outer boundary inward
    tau_refined = calculate_optical_depth(r_refined, rho_refined, kappa_refined)

    print(f"[Grid Refinement] Starting with {n_base} base points")
    print(f"[Grid Refinement] Optical depth range: {tau_refined[0]:.2e} to {tau_refined[-1]:.2e}")

    # --- Step 2: Iteratively add points where gradients exceed thresholds ---
    for iteration in range(N_ITERATIONS):
        if len(r_refined) >= MAX_GRID_POINTS:
            break

        # Find all consecutive pairs with gradient violations
        violations = []
        for i in range(len(r_refined) - 1):
            r1, r2 = r_refined[i], r_refined[i+1]
            v1, v2 = v_refined[i], v_refined[i+1]
            rho1, rho2 = rho_refined[i], rho_refined[i+1]
            t1, t2 = t_refined[i], t_refined[i+1]
            tau1, tau2 = tau_refined[i], tau_refined[i+1]

            # Calculate percentage jumps
            drho_pct = abs(rho2 - rho1) / rho1 * 100 if rho1 > 1e-30 else 0
            dt_pct = abs(t2 - t1) / t1 * 100 if t1 > 1e-6 else 0
            dv_pct = abs(v2 - v1) / v1 * 100 if v1 > 1e-10 else 0
            dr_r_pct = calculate_dr_over_r_log(r1, r2) if r1 > 1e-10 else 0  # Log-space calculation
            dtau_pct = abs(tau2 - tau1) / tau1 * 100 if tau1 > 1e-10 else 0  # Optical depth percentage jump

            # Check for violations (including dr/r constraint for radiation transfer stability)
            if (drho_pct > MAX_DRHO_PCT or
                dt_pct > MAX_DT_PCT or
                dv_pct > MAX_DV_PCT or
                dr_r_pct > MAX_DR_R_PCT or
                dtau_pct > MAX_DTAU_PCT):
                violations.append(i)

        if not violations:
            break  # No violations left

        # For each violation, insert a point at the geometric mean
        new_points = []
        for i in violations:
            if len(r_refined) >= MAX_GRID_POINTS:
                break
            r_mid = np.sqrt(r_refined[i] * r_refined[i+1])  # Geometric mean
            log_r_mid = np.log10(r_mid)

            # Interpolate from ORIGINAL data (not refined grid)
            v_mid = 10 ** np.interp(log_r_mid, log_r_oi, np.log10(v_oi))
            rho_mid = 10 ** np.interp(log_r_mid, log_r_oi, np.log10(rho_oi))
            t_mid = 10 ** np.interp(log_r_mid, log_r_oi, np.log10(t_oi))
            ne_mid = 10 ** np.interp(log_r_mid, log_r_oi, np.log10(np.maximum(ne_oi, 1e-10)))
            kappa_mid = np.interp(log_r_mid, log_r_oi, kappa_oi)

            new_points.append((i+1, r_mid, v_mid, rho_mid, t_mid, ne_mid, kappa_mid))

        # Insert points in reverse order to maintain indices
        for idx, r, v, rho, t, ne, kappa in sorted(new_points, key=lambda x: -x[0]):
            r_refined = np.insert(r_refined, idx, r)
            v_refined = np.insert(v_refined, idx, v)
            rho_refined = np.insert(rho_refined, idx, rho)
            t_refined = np.insert(t_refined, idx, t)
            ne_refined = np.insert(ne_refined, idx, ne)
            kappa_refined = np.insert(kappa_refined, idx, kappa)

        # Recalculate optical depth after inserting points
        tau_refined = calculate_optical_depth(r_refined, rho_refined, kappa_refined)

        print(f"[Grid Refinement] Iteration {iteration+1}: {len(r_refined)} points")

    # Check if dr/r constraint was satisfied (log-space calculation)
    max_dr_r = 0.0
    for i in range(len(r_refined) - 1):
        r1, r2 = r_refined[i], r_refined[i+1]
        dr_r = calculate_dr_over_r_log(r1, r2) if r1 > 1e-10 else 0  # Log-space
        max_dr_r = max(max_dr_r, dr_r)
    
    if max_dr_r > MAX_DR_R_PCT:
        print(f"[Grid Refinement] WARNING: Could not satisfy dr/r constraint ({max_dr_r:.2f}% > {MAX_DR_R_PCT}%)")
        print(f"[Grid Refinement] WARNING: Hit MAX_GRID_POINTS limit ({MAX_GRID_POINTS})")
        print(f"[Grid Refinement] WARNING: Consider increasing MAX_GRID_POINTS for better stability")
    else:
        print(f"[Grid Refinement] dr/r constraint satisfied: max dr/r = {max_dr_r:.2f}%")

    # --- Step 2.5: Post-refinement gradient smoothing ---
    # Detect gradient violations and apply localized smoothing to satisfy constraints
    max_smoothing_iterations = 30
    smoothing_window = 3

    # Save iteration data for plotting
    smoothing_history = {
        'iteration': [],
        'r': [],
        'rho': [],
        'v': [],
        'T': [],
        'violations': []
    }

    # Save initial state (iteration 0)
    T_k = t_refined * 1e4
    smoothing_history['iteration'].append(0)
    smoothing_history['r'].append(r_refined.copy())
    smoothing_history['rho'].append(rho_refined.copy())
    smoothing_history['v'].append(v_refined.copy())
    smoothing_history['T'].append(T_k.copy())

    # Record initial violations
    initial_violation_mask = detect_gradient_violations(
        r_refined, rho_refined, v_refined, T_k,
        max_drho_pct=MAX_DRHO_PCT,
        max_dt_pct=MAX_DT_PCT,
        max_dv_pct=MAX_DV_PCT
    )
    smoothing_history['violations'].append(np.sum(initial_violation_mask))
    print(f"[Gradient Smoothing] Initial state: {np.sum(initial_violation_mask)} violations")

    for smooth_iter in range(max_smoothing_iterations):
        # Detect violations
        violation_mask = detect_gradient_violations(
            r_refined, rho_refined, v_refined, T_k,
            max_drho_pct=MAX_DRHO_PCT,
            max_dt_pct=MAX_DT_PCT,
            max_dv_pct=MAX_DV_PCT
        )

        smoothing_history['violations'].append(np.sum(violation_mask))

        if not np.any(violation_mask):
            print(f"[Gradient Smoothing] All constraints satisfied after {smooth_iter} iterations")
            # Save final state
            smoothing_history['iteration'].append(smooth_iter + 1)
            smoothing_history['r'].append(r_refined.copy())
            smoothing_history['rho'].append(rho_refined.copy())
            smoothing_history['v'].append(v_refined.copy())
            smoothing_history['T'].append(T_k.copy())
            break

        print(f"[Gradient Smoothing] Iteration {smooth_iter + 1}: {np.sum(violation_mask)} violations")

        # Smooth violations
        rho_refined, v_refined, T_k = smooth_violations(
            r_refined, rho_refined, v_refined, T_k, violation_mask, window_size=smoothing_window
        )
        t_refined = T_k / 1e4  # Convert back to 10^4 K

        # Save state after this iteration
        smoothing_history['iteration'].append(smooth_iter + 1)
        smoothing_history['r'].append(r_refined.copy())
        smoothing_history['rho'].append(rho_refined.copy())
        smoothing_history['v'].append(v_refined.copy())
        smoothing_history['T'].append(T_k.copy())
    else:
        print(f"[Gradient Smoothing] WARNING: Could not satisfy all constraints after {max_smoothing_iterations} iterations")

    # Generate diagnostic plots if output_dir is provided
    if output_dir is not None:
        try:
            # Convert original temperature to K for plotting
            T_input_k = t_oi * 1e4
            plot_smoothing_evolution(smoothing_history, output_dir,
                                    r_input=r_oi, rho_input=rho_oi,
                                    v_input=v_oi, T_input=T_input_k)
        except ImportError:
            print("[Gradient Smoothing] matplotlib not available, skipping diagnostic plots")
        except Exception as e:
            print(f"[Gradient Smoothing] WARNING: Failed to generate smoothing plots: {e}")

    # --- Step 3: Coarsen grid by removing unnecessary points ---
    # Only remove points if:
    # 1. All adjacent gradients are below MIN_*_PCT (low-gradient region).
    # 2. Removing the point doesn't create a gradient violation.
    points_to_remove = []
    for i in range(1, len(r_refined) - 1):
        # Get neighbors
        r_prev, r_curr, r_next = r_refined[i-1], r_refined[i], r_refined[i+1]
        rho_prev, rho_curr, rho_next = rho_refined[i-1], rho_refined[i], rho_refined[i+1]
        t_prev, t_curr, t_next = t_refined[i-1], t_refined[i], t_refined[i+1]
        v_prev, v_curr, v_next = v_refined[i-1], v_refined[i], v_refined[i+1]

        # Compute gradients with neighbors
        drho_prev = abs(rho_curr - rho_prev) / rho_prev * 100 if rho_prev > 1e-30 else 0
        drho_next = abs(rho_next - rho_curr) / rho_curr * 100 if rho_curr > 1e-30 else 0
        dt_prev = abs(t_curr - t_prev) / t_prev * 100 if t_prev > 1e-6 else 0
        dt_next = abs(t_next - t_curr) / t_curr * 100 if t_curr > 1e-6 else 0
        dv_prev = abs(v_curr - v_prev) / v_prev * 100 if v_prev > 1e-10 else 0
        dv_next = abs(v_next - v_curr) / v_curr * 100 if v_curr > 1e-10 else 0
        dr_r_prev = calculate_dr_over_r_log(r_prev, r_curr) if r_prev > 1e-10 else 0  # Log-space
        dr_r_next = calculate_dr_over_r_log(r_curr, r_next) if r_curr > 1e-10 else 0  # Log-space

        # Check if this is a low-gradient region (including dr/r)
        if (max(drho_prev, drho_next) < MIN_DRHO_PCT and
            max(dt_prev, dt_next) < MIN_DT_PCT and
            max(dv_prev, dv_next) < MIN_DV_PCT and
            max(dr_r_prev, dr_r_next) < 1.0):  # dr/r already small (< 1%)
            # Check if removing this point would violate thresholds
            drho_combined = abs(rho_next - rho_prev) / rho_prev * 100 if rho_prev > 1e-30 else 0
            dt_combined = abs(t_next - t_prev) / t_prev * 100 if t_prev > 1e-6 else 0
            dv_combined = abs(v_next - v_prev) / v_prev * 100 if v_prev > 1e-10 else 0
            dr_r_combined = calculate_dr_over_r_log(r_prev, r_next) if r_prev > 1e-10 else 0  # Log-space

            if (drho_combined <= MAX_DRHO_PCT and
                dt_combined <= MAX_DT_PCT and
                dv_combined <= MAX_DV_PCT and
                dr_r_combined <= MAX_DR_R_PCT):
                points_to_remove.append(i)

    if points_to_remove:
        print(f"[Grid Refinement] Removing {len(points_to_remove)} unnecessary points")
        r_refined = np.delete(r_refined, points_to_remove)
        v_refined = np.delete(v_refined, points_to_remove)
        rho_refined = np.delete(rho_refined, points_to_remove)
        t_refined = np.delete(t_refined, points_to_remove)
        ne_refined = np.delete(ne_refined, points_to_remove)
        kappa_refined = np.delete(kappa_refined, points_to_remove)

        # Recalculate optical depth after removing points
        tau_refined = calculate_optical_depth(r_refined, rho_refined, kappa_refined)

    # --- Step 4: Enforce minimum dr to prevent numerical instability ---
    # Only remove a point if it is too close to its neighbour AND doing so does not
    # create a gradient violation (checked for rho, T, v).  This prevents Step 4 from
    # undoing the gradient enforcement performed in Step 2.5.
    print(f"[Grid Refinement] Enforcing minimum dr (0.1% of local radius) - dtau pruning disabled")
    i = 1
    points_removed = 0
    while i < len(r_refined) - 1:
        dr = r_refined[i] - r_refined[i-1]  # outer->inner: r decreasing, so dr > 0
        min_dr = MIN_DR_FRACTION * r_refined[i-1]
        if dr < min_dr:
            # Check that removing point i won't violate gradient limits with its neighbours
            if i + 1 < len(r_refined):
                rho_prev, rho_next = rho_refined[i-1], rho_refined[i+1]
                t_prev, t_next = t_refined[i-1], t_refined[i+1]
                v_prev, v_next = v_refined[i-1], v_refined[i+1]
                drho_combined = abs(rho_next - rho_prev) / max(rho_prev, 1e-30) * 100
                dt_combined   = abs(t_next - t_prev)   / max(t_prev, 1e-6) * 100
                dv_combined   = abs(v_next - v_prev)   / max(v_prev, 1e-10) * 100
                if drho_combined > MAX_DRHO_PCT or dt_combined > MAX_DT_PCT or dv_combined > MAX_DV_PCT:
                    # Removal would create a violation — skip this point
                    i += 1
                    continue
            r_refined   = np.delete(r_refined,   i)
            v_refined   = np.delete(v_refined,   i)
            rho_refined = np.delete(rho_refined, i)
            t_refined   = np.delete(t_refined,   i)
            ne_refined  = np.delete(ne_refined,  i)
            kappa_refined = np.delete(kappa_refined, i)
            tau_refined = np.delete(tau_refined, i)
            points_removed += 1
        else:
            i += 1
    if points_removed > 0:
        print(f"[Grid Refinement] Removed {points_removed} points to enforce minimum dr")

    # --- Step 5: Add outer boundary points with interpolated values ---
    # Re-enabled with interpolated density/temperature (not r^-2/constant T) and fewer points
    # to reduce grid size while avoiding artificial structure
    print(f"[Grid Refinement] Adding outer boundary points with interpolated values")
    n_outer_points = 20  # Reduced from 50 to keep grid size manageable
    spacing_fraction = 0.02  # 2% spacing
    
    r_outer = r_refined[-1]
    r_inner_outer = r_refined[-2] if len(r_refined) > 1 else r_refined[-1] * 0.9
    
    # Add points moving inward with 2% spacing
    for i in range(1, n_outer_points + 1):
        if len(r_refined) >= MAX_GRID_POINTS:
            break
        
        r_new = r_outer * (1 - i * spacing_fraction)
        
        if r_new <= r_inner_outer or r_new >= r_outer:
            break  # Reached inner boundary or invalid radius
        
        # Interpolate ALL quantities from ORIGINAL data (not refined grid) to avoid artificial structure
        log_r_new = np.log10(r_new)
        v_new = 10 ** np.interp(log_r_new, log_r_oi, np.log10(v_oi))
        rho_new = 10 ** np.interp(log_r_new, log_r_oi, np.log10(rho_oi))  # Interpolated density (not r^-2)
        t_new = 10 ** np.interp(log_r_new, log_r_oi, np.log10(t_oi))  # Interpolated temperature (not constant)
        kappa_new = np.interp(log_r_new, log_r_oi, kappa_oi)
        
        r_refined = np.insert(r_refined, -1, r_new)
        v_refined = np.insert(v_refined, -1, v_new)
        rho_refined = np.insert(rho_refined, -1, rho_new)
        t_refined = np.insert(t_refined, -1, t_new)
        ne_refined = np.insert(ne_refined, -1, 0.0)  # Placeholder
        kappa_refined = np.insert(kappa_refined, -1, kappa_new)

    # Recalculate optical depth after adding outer boundary points
    tau_refined = calculate_optical_depth(r_refined, rho_refined, kappa_refined)

    print(f"[Grid Refinement] Added {n_outer_points} outer boundary points with interpolated values")

    print(f"[Grid Refinement] Final grid: {len(r_refined)} points")

    # --- Step 6: Ensure grid is strictly monotonic (remove duplicate radius values) ---
    # CMFGEN requires strictly monotonic radius grid (no duplicate values)
    r_refined, v_refined, rho_refined, t_refined, ne_refined, kappa_refined, tau_refined = ensure_monotonic_grid(
        r_refined, v_refined, rho_refined, t_refined, ne_refined, kappa_refined, tau_refined
    )

    # --- Step 6.5: AMR for steep temperature gradients (Type IIn shocks) ---
    print(f"\n[Grid Refinement] Step 6.5: AMR for steep temperature gradients...")
    
    # Refine regions with steep temperature gradients (e.g., shock boundaries)
    TEMP_GRAD_THRESHOLD = 0.5  # Refine if |dlnT/dlnr| > 0.5
    MAX_TEMP_REFINE_POINTS = 50
    
    temp_grad_points_added = 0
    for iteration in range(5):
        dlnT = np.diff(np.log(t_refined))
        dlnR = np.diff(np.log(r_refined))
        temp_grad = np.abs(dlnT / dlnR)
        steep_indices = np.where(temp_grad > TEMP_GRAD_THRESHOLD)[0]
        
        if len(steep_indices) == 0 or temp_grad_points_added >= MAX_TEMP_REFINE_POINTS:
            break
        
        for idx in sorted(steep_indices, reverse=True):
            if len(r_refined) >= 10000 - 10:
                print(f"  WARNING: Approaching max grid points")
                break
            if temp_grad_points_added >= MAX_TEMP_REFINE_POINTS:
                break
            
            r_new = 0.5 * (r_refined[idx] + r_refined[idx+1])
            log_r_new = np.log10(r_new)
            
            v_new = 10 ** np.interp(log_r_new, np.log10(r_refined), np.log10(v_refined))
            rho_new = 10 ** np.interp(log_r_new, np.log10(r_refined), np.log10(rho_refined))
            t_new = 10 ** np.interp(log_r_new, np.log10(r_refined), np.log10(t_refined))
            kappa_new = np.interp(log_r_new, np.log10(r_refined), kappa_refined)
            
            r_refined = np.insert(r_refined, idx+1, r_new)
            v_refined = np.insert(v_refined, idx+1, v_new)
            rho_refined = np.insert(rho_refined, idx+1, rho_new)
            t_refined = np.insert(t_refined, idx+1, t_new)
            ne_refined = np.insert(ne_refined, idx+1, 0.0)
            kappa_refined = np.insert(kappa_refined, idx+1, kappa_new)
            
            temp_grad_points_added += 1
        
        tau_refined = calculate_optical_depth(r_refined, rho_refined, kappa_refined)
    
    print(f"[Grid Refinement] Added {temp_grad_points_added} points for temperature gradient AMR")
    
    # --- Step 6.5b: Global multi-physics gradient-based AMR ---
    # Applies across entire grid to ensure density, velocity, and temperature gradients
    # are well resolved everywhere. Prioritizes problematic regions but covers all depths.
    print(f"\n[Grid Refinement] Step 6.5b: Global multi-physics gradient AMR...")
    
    # Gradient thresholds for refinement (log-space gradients)
    # These are tighter than standard to ensure good resolution of all physics
    DRHO_GRAD_THRESHOLD = 0.25  # Refine if |dlnrho/dlnr| > 0.25
    DV_GRAD_THRESHOLD = 0.15    # Refine if |dlnv/dlnr| > 0.15
    DT_GRAD_THRESHOLD = 0.30     # Refine if |dlnT/dlnr| > 0.30
    
    # Prioritize problematic optical depth region (log tau ~ -6.5 to -5.5)
    # with tighter thresholds, but apply globally too
    FG_MIN_TAU = 10**(-6.5)  # ~3.16e-7
    FG_MAX_TAU = 10**(-5.5)  # ~3.16e-6
    
    MAX_PHYSICS_REFINE_POINTS = 200
    
    physics_points_added = 0
    for iteration in range(20):
        # Calculate all gradients in log-space across entire grid
        dlnr = np.abs(np.diff(np.log(r_refined)))
        dlnrho = np.abs(np.diff(np.log(rho_refined)))
        dlnv = np.abs(np.diff(np.log(v_refined)))
        dlnT = np.abs(np.diff(np.log(t_refined)))
        
        # Avoid division by zero
        dlnr = np.maximum(dlnr, 1e-30)
        
        # Compute normalized gradients
        rho_grad = dlnrho / dlnr
        v_grad = dlnv / dlnr
        T_grad = dlnT / dlnr
        
        # Check which regions we're in for priority weighting
        in_problematic_region = (tau_refined[:-1] >= FG_MIN_TAU) & (tau_refined[:-1] <= FG_MAX_TAU)
        
        # Find all gradient violations across entire grid
        violation_indices = []
        for idx in range(len(rho_grad)):
            # Use tighter thresholds in problematic region
            if in_problematic_region[idx]:
                thresh_rho = DRHO_GRAD_THRESHOLD * 0.8  # 20% tighter
                thresh_v = DV_GRAD_THRESHOLD * 0.8
                thresh_T = DT_GRAD_THRESHOLD * 0.8
            else:
                thresh_rho = DRHO_GRAD_THRESHOLD
                thresh_v = DV_GRAD_THRESHOLD
                thresh_T = DT_GRAD_THRESHOLD
            
            if (rho_grad[idx] > thresh_rho or 
                v_grad[idx] > thresh_v or 
                T_grad[idx] > thresh_T):
                violation_indices.append(idx)
        
        if len(violation_indices) == 0 or physics_points_added >= MAX_PHYSICS_REFINE_POINTS:
            if iteration == 0:
                print(f"  No gradient violations found - grid already well resolved")
            else:
                print(f"  No more violations after {iteration} iterations")
            break
        
        # Add points at violation locations (work from outside in to maintain indices)
        for idx in sorted(violation_indices, reverse=True):
            if len(r_refined) >= MAX_GRID_POINTS or physics_points_added >= MAX_PHYSICS_REFINE_POINTS:
                break
            
            # Geometric mean for new point
            r_new = np.sqrt(r_refined[idx] * r_refined[idx+1])
            log_r_new = np.log10(r_new)
            
            # Interpolate from original HERACLES data (not refined grid)
            v_new = 10 ** np.interp(log_r_new, log_r_oi, np.log10(v_oi))
            rho_new = 10 ** np.interp(log_r_new, log_r_oi, np.log10(rho_oi))
            t_new = 10 ** np.interp(log_r_new, log_r_oi, np.log10(t_oi))
            kappa_new = np.interp(log_r_new, log_r_oi, kappa_oi)
            
            # Insert new point
            r_refined = np.insert(r_refined, idx+1, r_new)
            v_refined = np.insert(v_refined, idx+1, v_new)
            rho_refined = np.insert(rho_refined, idx+1, rho_new)
            t_refined = np.insert(t_refined, idx+1, t_new)
            ne_refined = np.insert(ne_refined, idx+1, 0.0)
            kappa_refined = np.insert(kappa_refined, idx+1, kappa_new)
            
            physics_points_added += 1
        
        # Recalculate optical depth after grid modification
        tau_refined = calculate_optical_depth(r_refined, rho_refined, kappa_refined)
    
    print(f"[Grid Refinement] Added {physics_points_added} points for global multi-physics gradient AMR")
    print(f"  Thresholds: drho={DRHO_GRAD_THRESHOLD}, dv={DV_GRAD_THRESHOLD}, dT={DT_GRAD_THRESHOLD}")
    print(f"  Problematic region (log tau=-6.5 to -5.5) uses 20% tighter thresholds")
    
    # --- Step 6.6: Refine coarse outer boundary spacing ---
    print(f"\n[Grid Refinement] Step 6.6: Refining coarse outer boundary spacing...")
    
    # Target: outer 15% of grid where spacing becomes coarse
    n_outer_15 = max(1, int(0.15 * len(r_refined)))
    dr_ratio_threshold = 3.0  # Add point if dr[i]/dr[i+1] > 3 or dr[i+1]/dr[i] > 3
    max_outer_refine = 30
    
    outer_points_added = 0
    for iteration in range(3):
        if outer_points_added >= max_outer_refine:
            break
            
        # Check spacing in outer region
        dr_outer = np.diff(r_refined[:n_outer_15+1])
        if len(dr_outer) < 2:
            break
            
        # Find abrupt spacing changes
        spacing_ratios = np.abs(dr_outer[:-1] / np.maximum(dr_outer[1:], 1e-30))
        reverse_ratios = np.abs(dr_outer[1:] / np.maximum(dr_outer[:-1], 1e-30))
        
        coarse_indices = np.where((spacing_ratios > dr_ratio_threshold) | 
                                  (reverse_ratios > dr_ratio_threshold))[0]
        
        if len(coarse_indices) == 0:
            break
        
        # Add points at coarse spacing locations (work from outside in)
        for idx in sorted(coarse_indices, reverse=True):
            if len(r_refined) >= 10000 - 5 or outer_points_added >= max_outer_refine:
                break
                
            # Add point at problematic interval
            r_new = 0.5 * (r_refined[idx] + r_refined[idx+1])
            log_r_new = np.log10(r_new)
            
            v_new = 10 ** np.interp(log_r_new, np.log10(r_refined), np.log10(v_refined))
            rho_new = 10 ** np.interp(log_r_new, np.log10(r_refined), np.log10(rho_refined))
            t_new = 10 ** np.interp(log_r_new, np.log10(r_refined), np.log10(t_refined))
            kappa_new = np.interp(log_r_new, np.log10(r_refined), kappa_refined)
            
            r_refined = np.insert(r_refined, idx+1, r_new)
            v_refined = np.insert(v_refined, idx+1, v_new)
            rho_refined = np.insert(rho_refined, idx+1, rho_new)
            t_refined = np.insert(t_refined, idx+1, t_new)
            ne_refined = np.insert(ne_refined, idx+1, 0.0)
            kappa_refined = np.insert(kappa_refined, idx+1, kappa_new)
            
            outer_points_added += 1
        
        tau_refined = calculate_optical_depth(r_refined, rho_refined, kappa_refined)
    
    print(f"[Grid Refinement] Added {outer_points_added} points for outer boundary spacing")
    
    # --- Step 7: Enforce minimum dtau in outer grid for CMFGEN stability ---
    # Relaxed constraint: only enforce in regions where tau > 0.1 (not the very outer boundary)
    # The very outer optically thin regions (tau < 0.1) don't need fine dtau resolution
    # CMFGEN stability is only affected when dtau becomes extremely small across many zones
    print(f"[Grid Refinement] Enforcing dtau > 0.01 in outer regions where tau > 0.1")
    
    # Identify outer regions where tau > 0.1 (not the very thin outer boundary)
    n_points = len(r_refined)
    # Find where tau exceeds 0.1 in the outer portion
    tau_threshold_zone = 0
    for i in range(n_points):
        if tau_refined[i] > 0.1:
            tau_threshold_zone = i
            break
    
    # Only apply dtau enforcement from outer boundary up to tau = 0.1
    # Skip the very optically thin outermost regions where over-refinement occurs
    n_outer_for_dtau = max(1, tau_threshold_zone) if tau_threshold_zone > 0 else max(1, int(0.1 * n_points))
    
    print(f"  Applying dtau enforcement to outer {n_outer_for_dtau} zones (tau < 0.1)")
    
    # Calculate dtau for all intervals
    dtau = np.zeros(n_points - 1)
    for i in range(1, n_points):
        dtau[i-1] = abs(tau_refined[i-1] - tau_refined[i])
    
    # Check dtau in selected outer region
    outer_dtau = dtau[:n_outer_for_dtau]
    print(f"[Grid Refinement] Target region dtau range: {outer_dtau.min():.6f} to {outer_dtau.max():.6f}")
    print(f"[Grid Refinement] Target region tau range: {tau_refined[:n_outer_for_dtau].min():.6e} to {tau_refined[:n_outer_for_dtau].max():.6e}")
    
    # Add points where dtau <= 0.01 in target region
    points_added = 0
    max_iterations = 5  # Reduced from 10 to prevent over-refinement
    iteration = 0
    
    while iteration < max_iterations:
        # Recalculate dtau
        dtau = np.zeros(len(r_refined) - 1)
        for i in range(1, len(r_refined)):
            dtau[i-1] = abs(tau_refined[i-1] - tau_refined[i])
        
        # Check target outer region
        n_points_current = len(r_refined)
        tau_threshold_zone_current = 0
        for i in range(n_points_current):
            if tau_refined[i] > 0.1:
                tau_threshold_zone_current = i
                break
        n_outer_for_dtau_current = max(1, tau_threshold_zone_current) if tau_threshold_zone_current > 0 else max(1, int(0.1 * n_points_current))
        
        # Find intervals where dtau <= 0.01
        intervals_to_split = []
        for i in range(min(n_outer_for_dtau_current, len(dtau))):
            if dtau[i] <= 0.01:
                intervals_to_split.append(i)
        
        if len(intervals_to_split) == 0:
            print(f"[Grid Refinement] All target intervals have dtau > 0.01")
            break
        
        # Split intervals by adding intermediate points
        # Work from outermost inward to maintain monotonicity
        MIN_DR_OVER_R = 1e-5  # Minimum relative spacing to prevent over-refinement
        
        for i in sorted(intervals_to_split, reverse=True):
            if len(r_refined) >= MAX_GRID_POINTS_OUTER:
                print(f"[Grid Refinement] Reached MAX_GRID_POINTS_OUTER ({MAX_GRID_POINTS_OUTER}), stopping dtau enforcement")
                break
            
            # Check if spacing is already very small - skip if so
            dr_current = abs(r_refined[i+1] - r_refined[i])
            r_avg = 0.5 * (r_refined[i] + r_refined[i+1])
            if dr_current / r_avg < MIN_DR_OVER_R:
                continue  # Skip this interval - already fine enough
            
            # Add point at midpoint of interval
            r_new = 0.5 * (r_refined[i] + r_refined[i+1])
            
            # Interpolate quantities
            log_r_new = np.log10(r_new)
            v_new = 10 ** np.interp(log_r_new, np.log10(r_refined), np.log10(v_refined))
            rho_new = 10 ** np.interp(log_r_new, np.log10(r_refined), np.log10(rho_refined))
            t_new = 10 ** np.interp(log_r_new, np.log10(r_refined), np.log10(t_refined))
            kappa_new = np.interp(log_r_new, np.log10(r_refined), kappa_refined)
            
            # Insert new point
            r_refined = np.insert(r_refined, i+1, r_new)
            v_refined = np.insert(v_refined, i+1, v_new)
            rho_refined = np.insert(rho_refined, i+1, rho_new)
            t_refined = np.insert(t_refined, i+1, t_new)
            ne_refined = np.insert(ne_refined, i+1, 0.0)  # Placeholder
            kappa_refined = np.insert(kappa_refined, i+1, kappa_new)
            
            points_added += 1
        
        # Recalculate optical depth after adding points
        tau_refined = calculate_optical_depth(r_refined, rho_refined, kappa_refined)
        
        iteration += 1
        if iteration >= max_iterations:
            print(f"[Grid Refinement] Reached max iterations ({max_iterations}) for dtau enforcement")
            break
    
    print(f"[Grid Refinement] Added {points_added} points to enforce dtau > 0.01 (tau > 0.1 region only, min_spacing enforced)")
    
    # Final verification
    n_points_final = len(r_refined)
    n_outer_10_final = max(1, int(0.1 * n_points_final))
    dtau_final = np.zeros(n_points_final - 1)
    for i in range(1, n_points_final):
        dtau_final[i-1] = abs(tau_refined[i-1] - tau_refined[i])
    
    outer_dtau_final = dtau_final[:n_outer_10_final]
    print(f"[Grid Refinement] Final outer 10% dtau range: {outer_dtau_final.min():.6f} to {outer_dtau_final.max():.6f}")
    print(f"[Grid Refinement] Minimum dr in outer 10%: {np.min(np.abs(np.diff(r_refined[:n_outer_10_final]))):.6e}")
    
    # If we exceeded MAX_GRID_POINTS, trim back to MAX_GRID_POINTS
    if len(r_refined) > MAX_GRID_POINTS:
        print(f"[Grid Refinement] Exceeded MAX_GRID_POINTS ({MAX_GRID_POINTS}), trimming back")
        # Keep outer boundary points, trim from inner side
        n_remove = len(r_refined) - MAX_GRID_POINTS
        r_refined = r_refined[n_remove:]
        v_refined = v_refined[n_remove:]
        rho_refined = rho_refined[n_remove:]
        t_refined = t_refined[n_remove:]
        ne_refined = ne_refined[n_remove:]
        kappa_refined = kappa_refined[n_remove:]
        tau_refined = tau_refined[n_remove:]
        print(f"[Grid Refinement] Trimmed {n_remove} inner points to meet MAX_GRID_POINTS")



    # --- Diagnostic plot ---
    if output_dir is not None:
        try:
            import matplotlib.pyplot as plt
            point_iterations = np.zeros(len(r_refined), dtype=int)
            plot_grid_refinement_diagnostics(
                r_oi, v_oi, rho_oi, t_oi,
                r_refined, v_refined, rho_refined, t_refined,
                point_iterations, output_dir
            )
        except ImportError:
            print("[Grid Refinement] matplotlib not available, skipping diagnostic plots")
        except Exception as e:
            print(f"[Grid Refinement] WARNING: Failed to generate diagnostic plot: {e}")

    # Return in outer->inner order (CMFGEN format)
    return r_refined[::-1], v_refined[::-1], rho_refined[::-1], t_refined[::-1], ne_refined[::-1]


def plot_grid_refinement_diagnostics(r_orig, v_orig, rho_orig, t_orig, 
                                   r_refined, v_refined, rho_refined, t_refined,
                                   point_iterations, output_dir):
    """Generate diagnostic plots for grid refinement process with R on x-axis.
    
    All inputs should be in inner->outer order for consistency.
    """
    import matplotlib.pyplot as plt
    
    # Ensure all data is in inner->outer order
    # If original is in outer->inner, reverse it
    if r_orig[0] > r_orig[-1]:
        r_orig = r_orig[::-1]
        v_orig = np.abs(v_orig[::-1])
        rho_orig = rho_orig[::-1]
        t_orig = t_orig[::-1]
    
    # Ensure refined is in inner->outer order
    if r_refined[0] > r_refined[-1]:
        r_refined = r_refined[::-1]
        v_refined = v_refined[::-1]
        rho_refined = rho_refined[::-1]
        t_refined = t_refined[::-1]
    
    # Compute spacing constraints for refined grid
    n_refined = len(r_refined)
    drho_pct = np.zeros(n_refined - 1)
    dt_pct = np.zeros(n_refined - 1)
    dv_pct = np.zeros(n_refined - 1)
    dr_r_pct = np.zeros(n_refined - 1)
    
    for i in range(n_refined - 1):
        drho_pct[i] = (rho_refined[i+1] - rho_refined[i]) / rho_refined[i] * 100 if rho_refined[i] > 1e-30 else 0
        dt_pct[i] = (t_refined[i+1] - t_refined[i]) / t_refined[i] * 100 if t_refined[i] > 1e-6 else 0
        dv_pct[i] = (v_refined[i+1] - v_refined[i]) / v_refined[i] * 100 if v_refined[i] > 1e-10 else 0
        dr_r_pct[i] = (r_refined[i+1] - r_refined[i]) / r_refined[i] * 100 if r_refined[i] > 1e-10 else 0
    
    # Create figure with subplots (2 rows x 4 columns)
    fig, axes = plt.subplots(2, 4, figsize=(24, 12))

    
    # Upper row: Density, Temperature, Velocity profiles
    # Plot 1: Density profile (R on x-axis)
    ax = axes[0, 0]
    ax.step(r_orig / 1e10, rho_orig, 'o-', label='Original HERACLES', markersize=4, alpha=0.7, where='post')
    ax.step(r_refined / 1e10, rho_refined, 's-', label=f'Refined Grid ({len(r_refined)} pts)', markersize=3, alpha=0.7, where='post')
    ax.set_xlabel('Radius (10$^{10}$ cm)')
    ax.set_ylabel('Density (g/cm$^3$)')
    ax.set_title('Density Profile')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_xscale('log')
    ax.set_yscale('log')
    
    # Plot 2: Temperature profile (R on x-axis)
    ax = axes[0, 1]
    ax.step(r_orig / 1e10, t_orig * 1e4, 'o-', label='Original HERACLES', markersize=4, alpha=0.7, where='post')
    ax.step(r_refined / 1e10, t_refined * 1e4, 's-', label='Refined Grid', markersize=3, alpha=0.7, where='post')
    ax.set_xlabel('Radius (10$^{10}$ cm)')
    ax.set_ylabel('Temperature (K)')
    ax.set_title('Temperature Profile')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_xscale('log')
    ax.set_yscale('log')
    
    # Plot 3: Velocity profile (R on x-axis)
    ax = axes[0, 2]
    ax.step(r_orig / 1e10, v_orig / 1e5, 'o-', label='Original HERACLES', markersize=4, alpha=0.7, where='post')
    ax.step(r_refined / 1e10, v_refined / 1e5, 's-', label='Refined Grid', markersize=3, alpha=0.7, where='post')
    ax.set_xlabel('Radius (10$^{10}$ cm)')
    ax.set_ylabel('Velocity (km/s)')
    ax.set_title('Velocity Profile')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_xscale('log')
    ax.set_yscale('log')
    
    # Bottom row: Jumps
    # Plot 4: Density jumps vs R
    ax = axes[1, 0]
    r_centers = (r_refined[:-1] + r_refined[1:]) / 2
    ax.scatter(r_centers / 1e10, drho_pct, s=20, alpha=0.7, edgecolors='none')
    ax.axhline(y=MAX_DRHO_PCT, color='red', linestyle='--', label=f'Max (+{MAX_DRHO_PCT}%)')
    ax.axhline(y=-MAX_DRHO_PCT, color='red', linestyle='--', label=f'Min (-{MAX_DRHO_PCT}%)')
    ax.set_xlabel('Radius (10$^{10}$ cm)')
    ax.set_ylabel('Density Jumps ($\\frac{\\rho_{i+1} - \\rho_i}{\\rho_i} \\times 100$)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_xscale('log')
    
    # Plot 5: Temperature jumps vs R
    ax = axes[1, 1]
    ax.scatter(r_centers / 1e10, dt_pct, s=20, alpha=0.7, edgecolors='none')
    ax.axhline(y=MAX_DT_PCT, color='red', linestyle='--', label=f'Max (+{MAX_DT_PCT}%)')
    ax.axhline(y=-MAX_DT_PCT, color='red', linestyle='--', label=f'Min (-{MAX_DT_PCT}%)')
    ax.set_xlabel('Radius (10$^{10}$ cm)')
    ax.set_ylabel('Temperature Jumps ($\\frac{T_{i+1} - T_i}{T_i} \\times 100$)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_xscale('log')
    
    # Plot 6: Velocity jumps vs R
    ax = axes[1, 2]
    ax.scatter(r_centers / 1e10, dv_pct, s=20, alpha=0.7, edgecolors='none')
    ax.axhline(y=MAX_DV_PCT, color='red', linestyle='--', label=f'Max (+{MAX_DV_PCT}%)')
    ax.axhline(y=-MAX_DV_PCT, color='red', linestyle='--', label=f'Min (-{MAX_DV_PCT}%)')
    ax.set_xlabel('Radius (10$^{10}$ cm)')
    ax.set_ylabel('Velocity Jumps ($\\frac{v_{i+1} - v_i}{v_i} \\times 100$)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_xscale('log')
    
    # Plot 7: dr/r vs R (new column)
    ax = axes[1, 3]
    ax.scatter(r_centers / 1e10, dr_r_pct, s=20, alpha=0.7, edgecolors='none')
    ax.axhline(y=MAX_DR_R_PCT, color='red', linestyle='--', label=f'Max (+{MAX_DR_R_PCT}%)')
    ax.set_xlabel('Radius (10$^{10}$ cm)')
    ax.set_ylabel('dr/r ($\\frac{r_{i+1} - r_i}{r_i} \\times 100$)')
    ax.set_title('Radial Grid Spacing')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_xscale('log')
    
    plt.tight_layout()
    
    # Save plot
    plot_path = Path(output_dir) / "SN_HYDRO_DATA_rgrid.png"
    try:
        plt.savefig(plot_path, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"[Grid Refinement] Diagnostic plot saved to {plot_path}")
    except Exception as e:
        print(f"[Grid Refinement] WARNING: Failed to save diagnostic plot: {e}")
        plt.close()
    
    print(f"[Grid Refinement] Max jumps: rho={np.max(drho_pct):.2f}%, T={np.max(dt_pct):.2f}%, v={np.max(dv_pct):.2f}%")

def compute_gradient_weights(radius_cm, velocity_kms=None, temperature=None, density=None, gradient_threshold=0.5, weight_scale=2.0):
    """
    Compute a weighting function based on velocity, temperature, and density gradients.
    Regions with higher gradients will have higher weights.

    Parameters:
    -----------
    radius_cm : array
        Radius array in cm (outer to inner)
    velocity_kms : array (optional)
        Velocity array in km/s
    temperature : array (optional)
        Temperature array in K
    density : array (optional)
        Density array in g/cm³
    gradient_threshold : float
        Threshold for gradient normalization
    weight_scale : float
        Scale factor for weight range (final range = [1, 1 + weight_scale])
        Higher values concentrate more grid points in high-gradient regions

    Returns:
    --------
    weights : array
        Weighting function (higher where gradients are larger)
    """
    log_r = np.log10(radius_cm)
    weights = np.ones_like(radius_cm)

    # Velocity gradient contribution
    if velocity_kms is not None:
        log_velocity = np.log10(np.abs(velocity_kms) + 1e-300)
        velocity_gradient = np.abs(np.gradient(log_velocity, log_r))
        velocity_gradient_norm = velocity_gradient / (velocity_gradient.max() + 1e-300)
        weights += velocity_gradient_norm

    # Temperature gradient contribution
    if temperature is not None:
        log_temperature = np.log10(temperature + 1e-300)
        temperature_gradient = np.abs(np.gradient(log_temperature, log_r))
        temperature_gradient_norm = temperature_gradient / (temperature_gradient.max() + 1e-300)
        weights += temperature_gradient_norm
    
    # Density gradient contribution
    if density is not None:
        log_density = np.log10(density + 1e-300)
        density_gradient = np.abs(np.gradient(log_density, log_r))
        density_gradient_norm = density_gradient / (density_gradient.max() + 1e-300)
        weights += density_gradient_norm

    # Normalize weights to [1, 1 + weight_scale] for better resolution in high-gradient regions
    # Higher weight_scale concentrates more grid points where gradients are large
    weights = 1.0 + weight_scale * (weights - weights.min()) / (weights.max() - weights.min() + 1e-300)

    return weights

def generate_weighted_grid(
    radius_cm,
    density,
    velocity_kms=None,
    temperature=None,
    n_grid=100,
    gradient_threshold=0.5,
):
    """
    Generate a new radius grid weighted by velocity and temperature gradients,
    then interpolate density onto this grid.

    Parameters:
    -----------
    radius_cm : array
        Radius array in cm (outer to inner)
    density : array
        Density array in g/cm³
    velocity_kms : array (optional)
        Velocity array in km/s
    temperature : array (optional)
        Temperature array in K
    n_grid : int
        Target number of grid points
    gradient_threshold : float
        Threshold for gradient normalization

    Returns:
    --------
    r_new : array
        New radius grid (outer to inner)
    density_new : array
        Density interpolated onto the new grid
    """
    # Compute weights based on velocity, temperature, and density gradients
    weights = compute_gradient_weights(
        radius_cm,
        velocity_kms=velocity_kms,
        temperature=temperature,
        density=density,
        gradient_threshold=gradient_threshold,
    )

    # Normalize weights to create a probability distribution
    weights_normalized = weights / weights.sum()
    cumulative_weights = np.cumsum(weights_normalized)

    # Generate target cumulative weights for n_grid points
    target_cumulative = np.linspace(0, 1, n_grid)

    # Find the indices where the cumulative weights cross the target values
    indices = np.searchsorted(cumulative_weights, target_cumulative)
    indices = np.clip(indices, 0, len(radius_cm) - 1)

    # Get the new radius grid
    r_new = radius_cm[indices]

    # Interpolate density onto the new grid (in log-space for better accuracy)
    # Use quadratic interpolation to balance structure preservation and gradient control
    log_density = np.log10(density + 1e-300)
    from scipy import interpolate
    interp_func = interpolate.interp1d(radius_cm, log_density, kind='quadratic', fill_value='extrapolate')
    log_density_new = interp_func(r_new)
    density_new = 10**log_density_new

    print(f"  Weighted grid: {len(r_new)} points based on velocity/temperature gradients")
    print(f"  Weight range: {weights.min():.3f} to {weights.max():.3f}")

    return r_new, density_new

def calculate_tau_from_outer(r_cm, rho_cgs, kappa_cgs=None):
    """
    Calculate optical depth from the outside in (tau=0 at outer boundary).
    Arrays are in order outer->inner (r_cm[0]=outermost, r_cm[-1]=innermost).
    Optical depth increases with decreasing radius (as we go inward).
    """
    KAPPA_ES_CGS = 0.34  # Electron scattering opacity
    n = len(r_cm)
    
    if kappa_cgs is None:
        kappa = np.full_like(rho_cgs, KAPPA_ES_CGS, dtype=float)
    else:
        kappa = np.asarray(kappa_cgs, float)
    
    tau = np.zeros(n, dtype=float)
    if n < 2:
        return tau
    
    # Integrate from outer surface inward
    # tau[0] = 0 at outer boundary
    # tau increases as we go inward (to smaller radii)
    for i in range(1, n):
        dr_cm = abs(r_cm[i] - r_cm[i - 1])
        chimid = 0.5 * (kappa[i] * rho_cgs[i] + kappa[i - 1] * rho_cgs[i - 1])
        tau[i] = tau[i - 1] + max(chimid * dr_cm, 0.0)
    
    return tau

def trim_data_by_tau(r_cm, rho_cgs, velocity_cms, temperature_k, mass_fractions_dict,
                   kappa_cgs=None, tau_min=0.01, tau_max=30.0):
    """
    Trim hydro data to keep only zones within tau range [tau_min, tau_max].
    
    Args:
        r_cm: Radius in cm
        rho_cgs: Density in g/cm^3
        velocity_cms: Velocity in cm/s
        temperature_k: Temperature in K
        mass_fractions_dict: Dictionary of mass fraction arrays
        kappa_cgs: Opacity in cm^2/g (optional)
        tau_min: Minimum optical depth
        tau_max: Maximum optical depth
    
    Returns:
        Tuple of trimmed arrays (r_cm, rho_cgs, velocity_cms, temperature_k, mass_fractions_dict, kappa_cgs)
    """
    tau = calculate_tau_from_outer(r_cm, rho_cgs, kappa_cgs)
    n = len(tau)
    
    # Explicitly check and enforce outer->inner order for tau
    # tau should increase from outer (tau=0) to inner (tau=max)
    if tau[0] > tau[-1]:  # tau is in inner->outer order, reverse to outer->inner
        tau = tau[::-1]
        r_cm = r_cm[::-1]
        rho_cgs = rho_cgs[::-1]
        velocity_cms = velocity_cms[::-1]
        temperature_k = temperature_k[::-1]
        for species in mass_fractions_dict:
            mass_fractions_dict[species] = mass_fractions_dict[species][::-1]
        if kappa_cgs is not None:
            kappa_cgs = kappa_cgs[::-1]
    
    # Find indices where tau is within range
    # tau is calculated from outer->inner: tau[0] = 0 (outer), tau[-1] = max (inner)
    # Optical depth increases with decreasing radius (inward direction)
    # We want to keep zones where tau_min <= tau <= tau_max
    mask = (tau >= tau_min) & (tau <= tau_max)
    indices = np.where(mask)[0]
    
    if len(indices) == 0:
        print(f"WARNING: No zones found within tau range [{tau_min}, {tau_max}].")
        print(f"  tau range in data: [{tau.min():.3e}, {tau.max():.3e}]")
        # Find zones closest to the desired range
        if tau.max() < tau_min:
            # All zones have tau < tau_min, keep all
            print(f"  All zones have tau < tau_min. Keeping all zones.")
            i_min, i_max = 0, n - 1
        elif tau.min() > tau_max:
            # All zones have tau > tau_max, keep all
            print(f"  All zones have tau > tau_max. Keeping all zones.")
            i_min, i_max = 0, n - 1
        else:
            # Find the zone where tau crosses tau_min and tau_max
            i_min = np.where(tau >= tau_min)[0][0] if np.any(tau >= tau_min) else 0
            i_max = np.where(tau <= tau_max)[0][-1] if np.any(tau <= tau_max) else n - 1
            print(f"  Using closest zones: {i_min} to {i_max}")
    else:
        i_min, i_max = indices[0], indices[-1]
    print(f"Trimming hydro data: original ND={n}, tau range [{tau.min():.3e}, {tau.max():.3e}]")
    print(f"  Keeping zones {i_min} to {i_max} (ND_new={i_max - i_min + 1})")
    print(f"  New tau range: [{tau[i_max]:.3e}, {tau[i_min]:.3e}]")
    
    # Trim all arrays
    r_cm_trimmed = r_cm[i_min:i_max+1].copy()
    rho_cgs_trimmed = rho_cgs[i_min:i_max+1].copy()
    velocity_cms_trimmed = velocity_cms[i_min:i_max+1].copy()
    temperature_k_trimmed = temperature_k[i_min:i_max+1].copy()
    
    mass_fractions_trimmed = {}
    for species, mf in mass_fractions_dict.items():
        mass_fractions_trimmed[species] = mf[i_min:i_max+1].copy()
    
    if kappa_cgs is not None:
        kappa_cgs_trimmed = kappa_cgs[i_min:i_max+1].copy()
    else:
        kappa_cgs_trimmed = None
    
    # Post-trimming dr/r check (per Lechat review)
    dr_r_pct_trimmed = np.zeros(len(r_cm_trimmed) - 1)
    for i in range(len(r_cm_trimmed) - 1):
        dr_r_pct_trimmed[i] = calculate_dr_over_r_log(r_cm_trimmed[i], r_cm_trimmed[i+1]) if r_cm_trimmed[i] > 1e-10 else 0
    
    max_dr_r_trimmed = np.max(dr_r_pct_trimmed) if len(dr_r_pct_trimmed) > 0 else 0
    if max_dr_r_trimmed > MAX_DR_R_PCT:
        print(f"[trim_data_by_tau] WARNING: After trimming, max dr/r = {max_dr_r_trimmed:.2f}% > {MAX_DR_R_PCT}%")
        print(f"[trim_data_by_tau] WARNING: Consider re-running refinement on trimmed data")
    else:
        print(f"[trim_data_by_tau] dr/r after trimming: max = {max_dr_r_trimmed:.2f}%")
    
    return (r_cm_trimmed, rho_cgs_trimmed, velocity_cms_trimmed, temperature_k_trimmed,
            mass_fractions_trimmed, kappa_cgs_trimmed)

def interp_outer_to_inner(radius_old, values_old, radius_new, method="linear"):
    """
    Interpolate from outer-to-inner radius grid to new outer-to-inner radius grid.
    
    Args:
        radius_old: Original radius array (outer->inner)
        values_old: Original values array (outer->inner)
        radius_new: New radius array (outer->inner)
        method: Interpolation method ('linear' or 'cubic')
    
    Returns:
        Interpolated values array (outer->inner)
    """
    radius_old = np.asarray(radius_old)
    radius_new = np.asarray(radius_new)
    values_old = np.asarray(values_old)
    
    # Validate inputs
    assert np.all(radius_old > 0), "ERROR: radius_old must be strictly positive"
    assert np.all(radius_new > 0), "ERROR: radius_new must be strictly positive"
    assert np.all(np.isfinite(radius_old)), "ERROR: radius_old contains non-finite values"
    assert np.all(np.isfinite(radius_new)), "ERROR: radius_new contains non-finite values"
    assert np.all(np.isfinite(values_old)), "ERROR: values_old contains non-finite values"
    # Check monotonicity (either direction is fine since we reverse)
    old_diff = np.diff(radius_old)
    new_diff = np.diff(radius_new)
    assert np.all(old_diff != 0), "ERROR: radius_old contains duplicate points"
    assert np.all(new_diff != 0), "ERROR: radius_new contains duplicate points"
    
    x_old = np.log(radius_old[::-1])
    y_old = values_old[::-1]
    x_new = np.log(radius_new[::-1])

    if method == "cubic" and len(x_old) >= 4:
        from scipy.interpolate import CubicSpline
        y_new = CubicSpline(x_old, y_old, extrapolate=False)(x_new)
    else:
        y_new = np.interp(x_new, x_old, y_old)

    return y_new[::-1]

def validate_final_model(radius_cm, density, velocity_kms, temperature, tau, sigma, kappa, mass_fractions):
    """
    Validate that all final model arrays have the same length.
    
    Args:
        radius_cm: Radius array (cm)
        density: Density array (g/cm^3)
        velocity_kms: Velocity array (km/s)
        temperature: Temperature array (K)
        tau: Optical depth array
        sigma: Sigma array
        kappa: Opacity array (cm^2/g)
        mass_fractions: Dictionary of mass fraction arrays
    
    Returns:
        n: Number of grid points
    
    Raises:
        ValueError: If any array has a different length than radius_cm
    """
    n = len(radius_cm)

    required = {
        "density": density,
        "velocity_kms": velocity_kms,
        "temperature": temperature,
        "tau": tau,
        "sigma": sigma,
        "kappa": kappa,
    }

    for name, arr in required.items():
        if len(arr) != n:
            raise ValueError(f"{name} has length {len(arr)} but radius has length {n}")

    for name, arr in mass_fractions.items():
        if len(arr) != n:
            raise ValueError(f"mass fraction {name} has length {len(arr)} but radius has length {n}")

    # dr/r validation (per Lechat review)
    dr_r_pct = np.zeros(n - 1)
    for i in range(n - 1):
        dr_r_pct[i] = calculate_dr_over_r_log(radius_cm[i], radius_cm[i+1]) if radius_cm[i] > 1e-10 else 0
    
    max_dr_r = np.max(dr_r_pct) if len(dr_r_pct) > 0 else 0
    if max_dr_r > MAX_DR_R_PCT:
        print(f"[validate_final_model] WARNING: Final grid violates dr/r constraint! Max dr/r = {max_dr_r:.2f}%")
        print(f"[validate_final_model] WARNING: This may cause CMFGEN radiation transport instability")
    else:
        print(f"[validate_final_model] dr/r constraint satisfied: max dr/r = {max_dr_r:.2f}%")

    # Check gradient constraints
    velocity_cms = velocity_kms * 1e5  # Convert km/s to cm/s
    violation_mask = detect_gradient_violations(
        radius_cm, density, velocity_cms, temperature,
        max_drho_pct=MAX_DRHO_PCT,
        max_dt_pct=MAX_DT_PCT,
        max_dv_pct=MAX_DV_PCT
    )

    if np.any(violation_mask):
        print(f"[validate_final_model] WARNING: {np.sum(violation_mask)} gradient violations remain!")
        print(f"[validate_final_model] WARNING: Consider increasing smoothing iterations or window size")
    else:
        print(f"[validate_final_model] All gradient constraints satisfied")

    return n

def convert_heracles_to_cmfgen(heracles_dir, target_days=30, output_file='SN_HYDRO_DATA', 
                               tau_min=None, tau_max=None, n_grid=100, adaptive_grid=True, 
                               gradient_threshold=0.5, initial_profile=None,
                               smooth_sigma=False, sigma_window=9, sigma_polyorder=3, sigma_threshold=2.0,
                               velocity_max_kms=None):
    """Convert HERACLES snapshot to CMFGEN SN_HYDRO_DATA format
    
    Args:
        heracles_dir: Directory containing HERACLES snapshots
        target_days: Target snapshot time in days
        output_file: Output filename for SN_HYDRO_DATA
        tau_min: Minimum optical depth for trimming (optional)
        tau_max: Maximum optical depth for trimming (optional)
        n_grid: Number of grid points for CMFGEN (regridding target)
        adaptive_grid: Use adaptive grid based on density gradients
        gradient_threshold: Threshold for density gradient to trigger adaptive sampling
        initial_profile: Path to initial profile file with temperature data (optional)
        smooth_sigma: Enable Savitzky-Golay smoothing of Sigma values (default: False)
        sigma_window: Window length for Savitzky-Golay filter (default: 9)
        sigma_polyorder: Polynomial order for Savitzky-Golay filter (default: 3)
        sigma_threshold: Sigma threshold above which smoothing is applied (default: 2.0)
        velocity_max_kms: Maximum velocity in km/s for trimming outer layers (optional)
    """
    
    print(f"Finding HERACLES snapshot closest to {target_days} days...")
    snap_dir, actual_days, diff = find_closest_heracles_snapshot(target_days, heracles_dir)
    print(f"Found snapshot: {snap_dir} at {actual_days} days (diff: {diff:.1f} days)")
    
    # Flag to track if we're using initial profile temperature
    use_initial_temperature = False
    
    # Read initial profile if provided
    if initial_profile is not None and os.path.exists(initial_profile):
        print(f"Reading initial profile from {initial_profile}...")
        try:
            profile_data = np.loadtxt(initial_profile)
            # Profile format: r_cm rho_cgs p_cgs u_cgs fx_H fx_He fx_C fx_O fx_Si fx_Fe
            initial_r = profile_data[:, 0]  # cm
            initial_rho = profile_data[:, 1]  # g/cm³
            initial_p = profile_data[:, 2]  # dyn/cm²
            # Calculate temperature from ideal gas law: T = p * mu * m_H / (rho * k_B)
            # Use mean molecular weight mu = 0.6 for typical H/He mixture
            k_B = 1.38062e-16  # erg/K
            m_H = 1.673e-24  # g
            mu = 0.6
            initial_T = initial_p * mu * m_H / (initial_rho * k_B)
            initial_T = np.maximum(initial_T, 1000.0)  # Floor at 1000 K
            print(f"  Initial profile radius range: {initial_r.min():.2e} to {initial_r.max():.2e} cm")
            print(f"  Initial profile temperature range: {initial_T.min():.2e} to {initial_T.max():.2e} K")
            use_initial_temperature = True
        except Exception as e:
            print(f"  Warning: Could not read initial profile: {e}")
            print(f"  Will use HERACLES-calculated temperature instead")
    
    # Read HERACLES snapshot
    print(f"Reading HERACLES snapshot from {snap_dir}...")
    
    # Check what files are available in the snapshot directory
    snap_files = os.listdir(snap_dir)
    print(f"Files in snapshot directory: {snap_files}")
    
    # Extract the snapshot number from the directory path
    snap_num = os.path.basename(snap_dir)
    
    # Try to read using rd_heracles HeraclesData
    try:
        print(f"Reading HERACLES data using HeraclesData (snapshot {snap_num})...")
        d = rdh.HeraclesData(idump=int(snap_num), form='bin', verbose=False,
                             dir=heracles_dir)
        
        # Extract physical quantities from HERACLES data structure
        radius_cm = np.asarray(d.xcg).ravel()          # cm
        density   = np.asarray(d.rho).ravel()          # g/cm^3
        velocity_cms = np.asarray(d.u)[..., 0].ravel() # cm/s
        velocity_kms = velocity_cms / 1e5              # km/s
        
        # Check if HERACLES provides optical depth data
        if hasattr(d, 'tau'):
            tau_heracles = np.asarray(d.tau).ravel()
            print(f"Found HERACLES optical depth data: min={tau_heracles.min():.2e}, max={tau_heracles.max():.2e}")
        else:
            tau_heracles = None
            print("HERACLES does not provide optical depth data - will calculate from opacity")
        
        # Derive temperature: prefer T_rad from radiation field (hydrorad runs),
        # fall back to Newton-Raphson inversion of E_tot (pure hydro runs).
        A_VALS = [1.008, 4.003, 12.011, 15.999, 28.085, 55.845]
        Z_VALS = [1.0,   2.0,   6.0,    8.0,   14.0,   26.0  ]

        # Mass fractions (needed for mu and later composition)
        fx_arr = np.asarray(getattr(d, "fx", np.zeros((len(radius_cm), 6)))).reshape(len(radius_cm), -1)
        nfx    = fx_arr.shape[1]
        mf_raw = [fx_arr[:, i] if i < nfx else np.zeros(len(radius_cm)) for i in range(6)]
        mf_sum = np.maximum(np.sum(mf_raw, axis=0), 1e-300)
        mf_n   = [mf_raw[i] / mf_sum for i in range(6)]
        inv_mu = np.zeros(len(radius_cm))
        for i in range(6):
            inv_mu += mf_n[i] * (1.0 + Z_VALS[i]) / A_VALS[i]
        mu_loc = 1.0 / np.maximum(inv_mu, 1e-30)

        er_field = getattr(d, 'er', None)
        if er_field is not None:
            # Hydrorad run: use T_rad = (E_rad / a_rad)^(1/4) directly.
            # This is the self-consistent radiation temperature solved by HERACLES M1.
            a_rad = 7.5657e-15  # erg cm^-3 K^-4
            er = np.asarray(er_field).ravel()
            er_safe = np.maximum(er, 0.0)
            temperature = (er_safe / a_rad) ** 0.25
            temperature = np.where(np.isfinite(temperature) & (temperature > 0), temperature, 1000.0)
            print(f"Temperature from HERACLES radiation field T_rad = (E_rad/a)^(1/4)")
            print(f"  E_rad range: {er_safe.min():.2e} to {er_safe.max():.2e} erg/cm³")
            print(f"  Temperature range: {temperature.min():.2e} to {temperature.max():.2e} K")
        else:
            # Pure hydro run: invert total energy via Newton-Raphson
            E_tot = np.asarray(d.E).ravel()
            temperature = calculate_temperature_with_radiation(density, velocity_cms, E_tot, mu_loc)
            temperature = np.where(np.isfinite(temperature) & (temperature > 0), temperature, 1000.0)
            print(f"Temperature derived from energy field via Newton-Raphson (no radiation field available)")
            print(f"  Mean molecular weight mu: min={mu_loc.min():.3e}, max={mu_loc.max():.3e}")
        print(f"  Density: min={density.min():.3e}, max={density.max():.3e} g/cm³")
        print(f"  Temperature (from HERACLES): min={temperature.min():.3e}, max={temperature.max():.3e} K")
        
        # Smooth temperature to remove extreme spikes that cause NRAY_MAX errors in CMFGEN
        # Use a simple moving average to spread spikes over adjacent cells
        def smooth_temperature(temp, window_size=5):
            """Apply moving average smoothing to temperature array."""
            temp_smoothed = np.copy(temp)
            half_window = window_size // 2
            for i in range(len(temp)):
                start = max(0, i - half_window)
                end = min(len(temp), i + half_window + 1)
                temp_smoothed[i] = np.mean(temp[start:end])
            return temp_smoothed
        
        # Disable temperature smoothing - grid refinement already handles gradients
        # temperature = smooth_temperature(temperature, window_size=3)
        # print(f"  Temperature (after smoothing): min={temperature.min():.3e}, max={temperature.max():.3e} K")
        print(f"  Temperature (no smoothing): min={temperature.min():.3e}, max={temperature.max():.3e} K")
        
        # Smooth velocity to reduce dV/dr gradient that causes NRAY_MAX errors in CMFGEN
        # Use a simple moving average to spread steep velocity transitions over adjacent cells
        def smooth_velocity(velocity, window_size=5):
            """Apply moving average smoothing to velocity profile."""
            if window_size < 1:
                window_size = 1
            if window_size % 2 == 0:
                window_size += 1  # Ensure odd number
            
            # Pad with edge values
            pad_width = window_size // 2
            velocity_padded = np.pad(velocity, (pad_width, pad_width), mode='edge')
            
            # Apply moving average
            smoothed = np.convolve(velocity_padded, np.ones(window_size)/window_size, mode='valid')
            
            return smoothed
        
        # No sigma clipping - pass values directly to CMFGEN

        # Apply sigmoid smoothing to transition between ejecta A (H-rich) and ejecta B (metal-rich)
        # This creates a smooth temperature transition instead of a sharp cutoff
        # DISABLE sigmoid smoothing when using initial profile temperature to preserve intended gradient
        if not use_initial_temperature:
            H_fraction = mf_raw[0]  # Hydrogen mass fraction
            T_original = temperature.copy()
            
            # Sigmoid temperature overwrite DISABLED - use raw solver output
            # This was forcing metal-rich zones toward 5000 K, which is arbitrary and not physics-based
            print(f"  Temperature (sigmoid smoothing DISABLED - using raw solver output)")
            print(f"    Temperature range: {temperature.min():.2e} to {temperature.max():.2e} K")
        else:
            print(f"  Temperature (sigmoid smoothing DISABLED - using initial profile temperature)")
            print(f"    Temperature range: {temperature.min():.2e} to {temperature.max():.2e} K")
        
        # Get mass fractions
        fx = np.asarray(getattr(d, "fx", np.zeros((len(radius_cm), 6))))
        if len(fx.shape) == 4:
            fx = fx.reshape(-1, 6)
        
        # --- DEBUG: Validate and fix mass fractions ---
        print("\n=== MASS FRACTION DEBUG ===")
        print(f"fx shape: {fx.shape}")
        print(f"fx min/max: {fx.min():.3e} / {fx.max():.3e}")
        
        # Ensure fx is 2D
        if len(fx.shape) != 2:
            fx = fx.reshape(-1, fx.shape[-1])
        
        # Fill missing columns with zeros if nfx < 6
        nfx = fx.shape[1] if len(fx.shape) > 1 else 0
        if nfx < 6:
            print(f"WARNING: HERACLES provides only {nfx} mass fractions (expected 6). Filling missing with zeros.")
            fx_padded = np.zeros((fx.shape[0], 6))
            fx_padded[:, :nfx] = fx
            fx = fx_padded
        
        # Extract mass fractions and canonicalize to CMFGEN names immediately
        mass_fractions = {}
        # HERACLES column order: H, He, C, O, Si, Fe
        # CMFGEN canonical names: HYD, HE, CARB, OXY, SIL, IRON
        elem_name_map = {
            'H': 'HYD',
            'He': 'HE',
            'C': 'CARB',
            'O': 'OXY',
            'Si': 'SIL',
            'Fe': 'IRON'
        }
        heracles_names = ['H', 'He', 'C', 'O', 'Si', 'Fe']
        for i, elem in enumerate(heracles_names):
            cmfgen_name = elem_name_map[elem]
            mass_fractions[cmfgen_name] = fx[:, i]
        
        # Clip to [0, 1] and ensure sum is ~1
        for elem in mass_fractions:
            mass_fractions[elem] = np.clip(mass_fractions[elem], 0.0, 1.0)
        
        # Apply floor: set mass fractions < 1e-8 to 0.0 to avoid numerical noise
        MASS_FRACTION_FLOOR = 1e-8
        for elem in mass_fractions:
            mass_fractions[elem] = np.where(mass_fractions[elem] < MASS_FRACTION_FLOOR, 0.0, mass_fractions[elem])
        
        mf_sum = np.sum([mass_fractions[elem] for elem in mass_fractions], axis=0)
        if np.any(mf_sum < 1e-10):
            print(f"WARNING: Mass fraction sum is near zero in some zones. Preserving dominant species.")
            # Find dominant species at each zone and set it to 1.0 where sum is near zero
            zero_mask = mf_sum < 1e-10
            # Stack all mass fractions to find argmax
            all_mf = np.stack([mass_fractions[elem] for elem in mass_fractions], axis=0)
            dominant_idx = np.argmax(all_mf, axis=0)
            elem_list = list(mass_fractions.keys())
            for i, elem in enumerate(elem_list):
                is_dominant = zero_mask & (dominant_idx == i)
                mass_fractions[elem] = np.where(is_dominant, 1.0, mass_fractions[elem])
                mass_fractions[elem] = np.where(zero_mask & (dominant_idx != i), 0.0, mass_fractions[elem])
            mf_sum = np.sum([mass_fractions[elem] for elem in mass_fractions], axis=0)
        
        # Renormalize
        mf_sum_safe = np.where(mf_sum > 1e-10, mf_sum, 1.0)
        for elem in mass_fractions:
            mass_fractions[elem] /= mf_sum_safe
        
        # Verify
        mf_sum_final = np.sum([mass_fractions[elem] for elem in mass_fractions], axis=0)
        print(f"Mass fraction sum after fixes: min={mf_sum_final.min():.6f}, max={mf_sum_final.max():.6f}")
        print("=== END MASS FRACTION DEBUG ===\n")
        
        # Print element statistics (using canonical CMFGEN names)
        canonical_elems = ['HYD', 'HE', 'CARB', 'OXY', 'SIL', 'IRON']
        for elem in canonical_elems:
            if elem in mass_fractions:
                mf = mass_fractions[elem]
                print(f"{elem}: min={mf.min():.3e}, max={mf.max():.3e}, std={mf.std():.3e}")
        
        # Add additional elements with zeros if needed (canonical CMFGEN names)
        if 'NITR' not in mass_fractions:
            mass_fractions['NITR'] = np.zeros(len(radius_cm))
        if 'NEON' not in mass_fractions:
            mass_fractions['NEON'] = np.zeros(len(radius_cm))
        if 'MAGN' not in mass_fractions:
            mass_fractions['MAGN'] = np.zeros(len(radius_cm))
        if 'SULF' not in mass_fractions:
            mass_fractions['SULF'] = np.zeros(len(radius_cm))
        
        print(f"\nBEFORE REGRIDDING (n_points={len(radius_cm)}):")
        for elem in canonical_elems:
            if elem in mass_fractions:
                mf = mass_fractions[elem]
                print(f"{elem}: min={mf.min():.3e}, max={mf.max():.3e}, std={mf.std():.3e}")
                print(f"  First 5: {mf[:5]}")
                print(f"  Last 5: {mf[-5:]}")
        
        n_points = len(radius_cm)
        print(f"Successfully read HERACLES data: {n_points} points")
        print(f"  Radius range: {radius_cm.min():.2e} to {radius_cm.max():.2e} cm")
        print(f"  Velocity range: {velocity_kms.min():.2f} to {velocity_kms.max():.2f} km/s")
        print(f"  Density range: {density.min():.2e} to {density.max():.2e} g/cm^3")
        print(f"  Temperature range: {temperature.min():.2e} to {temperature.max():.2e} K")
        
    except Exception as e:
        print(f"ERROR reading HERACLES data: {e}")
        print("HERACLES data is required for temperature calculation")
        return
    
    # STEP 1: Load original model (already done above)
    # STEP 2: Reverse once so arrays are outer-to-inner
    # HERACLES grid is inner->outer (radius increasing).
    # CMFGEN SN_HYDRO_DATA expects outer->inner (radius decreasing).
    # Reverse all arrays now so index 0 = outermost point throughout.
    radius_cm    = radius_cm[::-1].copy()
    density      = density[::-1].copy()
    velocity_cms = velocity_cms[::-1].copy()
    temperature  = temperature[::-1].copy()
    for elem in mass_fractions:
        mass_fractions[elem] = mass_fractions[elem][::-1].copy()
    
    print(f"After reversal — radius[0]={radius_cm[0]:.3e} cm (outer), radius[-1]={radius_cm[-1]:.3e} cm (inner)")
    print(f"Velocity range: min={velocity_cms.min():.2f} to max={velocity_cms.max():.2f} cm/s")
    assert radius_cm[0] > radius_cm[-1], "ERROR: After reversal, radius should be outer->inner (decreasing)"
    
    # Save original data for comparison plotting
    original_radius_cm = radius_cm.copy()
    original_density = density.copy()
    original_temperature = temperature.copy()
    original_velocity_cms = velocity_cms.copy()
    original_velocity_kms = velocity_cms / 1e5
    original_mass_fractions = {elem: arr.copy() for elem, arr in mass_fractions.items()}
    
    # Calculate original sigma for comparison plotting
    original_radius = original_radius_cm / 1e10  # 10^10 cm
    original_sigma_unclipped = calculate_sigma(original_radius, original_velocity_kms)
    
    # STEP 3: Compute original tau
    # Calculate Rosseland mean opacity including electron scattering + bound-free + free-free
    X = mass_fractions.get('HYD', np.zeros(len(radius_cm)))
    Y = mass_fractions.get('HE', np.zeros(len(radius_cm)))
    
    # Calculate metal mass fraction Z (sum of all metals except H and He)
    Z = np.zeros(len(radius_cm))
    for elem in mass_fractions:
        if elem not in ['HYD', 'HE']:
            Z += mass_fractions[elem]
    
    # Calculate electron density with temperature-dependent ionization
    electron_density = calculate_electron_density(density, temperature, mass_fractions)
    
    # Calculate Rosseland mean opacity using actual electron density
    kappa_cm2_per_g_array = approximate_rosseland_opacity(temperature, density, X, Y, Z)
    
    print(f"\nCalculating optical depth using Rosseland mean opacity (ES + bf + ff)...")
    print(f"  Electron density range: {electron_density.min():.2e} to {electron_density.max():.2e} cm^-3")
    print(f"  κ range: {kappa_cm2_per_g_array.min():.3e} to {kappa_cm2_per_g_array.max():.3e} cm²/g")
    
    tau = calculate_tau_from_outer(radius_cm, density, kappa_cm2_per_g_array)
    print(f"  Optical depth range: {tau.min():.2e} to {tau.max():.2e}")
    
    # Save original tau for comparison plotting
    original_tau = tau.copy()
    
    # STEP 4: Trim by tau_min and tau_max
    # Use CLI values if provided, otherwise use defaults
    if tau_min is None:
        tau_min = 1e-4
    if tau_max is None:
        tau_max = 30.0
    print(f"Using tau range: [{tau_min}, {tau_max}]")
    
    # Find zone just behind tau_max for inner boundary (NOT in SN_HYDRO_DATA)
    i_boundary = None
    r_inner_boundary_cm = None
    T_inner_boundary_K = None
    if np.any(tau > tau_max):
        i_max = np.where(tau <= tau_max)[0][-1] if np.any(tau <= tau_max) else 0
        i_boundary = min(i_max + 1, len(tau) - 1)
        r_inner_boundary_cm = radius_cm[i_boundary]
        T_inner_boundary_K = temperature[i_boundary]
        print(f"Inner boundary zone (not in SN_HYDRO_DATA): depth {i_boundary}, "
              f"T={T_inner_boundary_K:.1f} K, R={r_inner_boundary_cm:.2e} cm")
    else:
        # All zones below tau_max, use innermost zone
        r_inner_boundary_cm = radius_cm[-1]
        T_inner_boundary_K = temperature[-1]
    
    print(f"\nTrimming HERACLES snapshot to optical depth range [{tau_min}, {tau_max}]")
    
    (radius_cm, density, velocity_cms, temperature, mass_fractions, kappa_cm2_per_g_array) = \
        trim_data_by_tau(radius_cm, density, velocity_cms, temperature, mass_fractions,
                        kappa_cgs=kappa_cm2_per_g_array, tau_min=tau_min, tau_max=tau_max)
    
    # Validate mass fractions after trimming
    mf_sum = np.zeros(len(radius_cm))
    for elem in mass_fractions:
        mf_sum += mass_fractions[elem]
    if np.any(np.abs(mf_sum - 1.0) > 1e-6):
        print(f"WARNING: Mass fractions not normalized! Sum range: {mf_sum.min():.6f} to {mf_sum.max():.6f}")
        # Renormalize if needed
        for elem in mass_fractions:
            mass_fractions[elem] /= mf_sum
    
    print(f"After trimming: {len(radius_cm)} points")
    print(f"  Radius range: {radius_cm.min():.2e} to {radius_cm.max():.2e} cm")
    assert radius_cm[0] > radius_cm[-1], "ERROR: After trimming, radius should be outer->inner (decreasing)"
    
    # STEP 5: Trim by velocity_max_kms if specified
    velocity_trimmed_indices = None
    if velocity_max_kms is not None:
        print(f"\nTrimming outer layers where velocity > {velocity_max_kms} km/s")
        velocity_kms = velocity_cms / 1e5
        # Keep zones where velocity <= velocity_max_kms
        # Outer boundary has highest velocity, so we trim from the outer end
        mask = velocity_kms <= velocity_max_kms
        velocity_trimmed_indices = np.where(~mask)[0]  # Indices that will be trimmed
        
        if len(velocity_trimmed_indices) > 0:
            print(f"  Trimming {len(velocity_trimmed_indices)} zones with velocity > {velocity_max_kms} km/s")
            print(f"  Velocity range before: {velocity_kms.min():.2f} to {velocity_kms.max():.2f} km/s")
            
            # Trim all arrays
            radius_cm = radius_cm[mask]
            density = density[mask]
            velocity_cms = velocity_cms[mask]
            temperature = temperature[mask]
            for elem in mass_fractions:
                mass_fractions[elem] = mass_fractions[elem][mask]
            kappa_cm2_per_g_array = kappa_cm2_per_g_array[mask]
            
            velocity_kms = velocity_cms / 1e5
            print(f"  Velocity range after: {velocity_kms.min():.2f} to {velocity_kms.max():.2f} km/s")
            print(f"  New ND: {len(radius_cm)}")
        else:
            print(f"  No zones exceed velocity limit of {velocity_max_kms} km/s")
            velocity_trimmed_indices = np.array([])
    
    # STEP 6: Recompute trimmed tau, velocity_kms, sigma
    velocity_kms = velocity_cms / 1e5
    tau = calculate_tau_from_outer(radius_cm, density, kappa_cm2_per_g_array)
    
    radius = radius_cm / 1e10  # 10^10 cm
    sigma_unclipped = calculate_sigma(radius, velocity_kms)
    sigma = sigma_unclipped  # Pass sigma values directly to CMFGEN
    
    print(f"  Velocity range: {velocity_kms.min():.2f} to {velocity_kms.max():.2f} km/s")
    print(f"  Optical depth range: {tau.min():.2e} to {tau.max():.2e}")
    
    # Save trimmed data for comparison plotting
    trimmed_radius_cm = radius_cm.copy()
    trimmed_density = density.copy()
    trimmed_temperature = temperature.copy()
    trimmed_velocity_cms = velocity_cms.copy()
    trimmed_velocity_kms = velocity_kms.copy()
    trimmed_mass_fractions = {elem: arr.copy() for elem, arr in mass_fractions.items()}
    trimmed_tau = tau.copy()
    trimmed_sigma = sigma.copy()
    
    # STEP 6: Build refined grid using adaptive grid refinement (replaces AMR)
    print(f"\nApplying adaptive grid refinement (replaces AMR)...")
    
    # Get output directory for diagnostic plots
    output_dir = Path(output_file).parent if output_file else Path(".")
    
    # Apply adaptive grid refinement
    # Note: electron_density not computed yet, will be None (refinement handles this)
    r_new, velocity_cms_new, density_new, temperature_new, ed_new = refine_grid(
        radius_cm=radius_cm,
        vel_cms=velocity_cms,
        rho_cgs=density,
        temp_104k=temperature / 1e4,  # Convert K to 10^4 K
        ed_cgs=None,  # Will be computed later
        kappa_cgs=kappa_cm2_per_g_array,
        output_dir=output_dir
    )
    
    # Convert back to original units
    temperature_new = temperature_new * 1e4  # Convert 10^4 K back to K
    velocity_kms_new = velocity_cms_new / 1e5
    
    print(f"  New radius range: {r_new[0]:.2e} to {r_new[-1]:.2e} cm")
    assert r_new[0] > r_new[-1], "ERROR: After refine_grid, radius should be outer->inner (decreasing)"
    
    # STEP 7: Adaptive grid refinement already interpolated all quantities
    # Just update the main arrays with refined values
    # refine_grid returns same order as input (outer->inner)
    velocity_cms = velocity_cms_new
    velocity_kms = velocity_cms / 1e5
    density = density_new
    temperature = temperature_new
    # electron_density will be computed later from refined density
    
    # Apply velocity smoothing to reduce dV/dr gradient
    velocity_cms = smooth_velocity(velocity_cms, window_size=5)
    velocity_kms = velocity_cms / 1e5
    print(f"  Velocity (after smoothing): min={velocity_kms.min():.2f} to {velocity_kms.max():.2f} km/s")


# STEP 8: Interpolate remaining arrays (opacity, mass fractions) onto refined grid
    # Adaptive grid refinement already handled R, V, rho, T - only interpolate opacity and mass fractions
    # Note: r_new is in outer->inner order from adaptive grid refinement (CMFGEN format)
    # radius_cm is in outer->inner order, so no reversal needed for interpolation
    kappa_cm2_per_g_array = interp_outer_to_inner(radius_cm, kappa_cm2_per_g_array, r_new, method="linear")
    
    # Interpolate mass fractions
    print("\nInterpolating mass fractions onto refined grid...")
    for elem in mass_fractions:
        mass_fractions[elem] = interp_outer_to_inner(radius_cm, mass_fractions[elem], r_new, method="linear")
        # Clip to valid range and replace NaN/inf with 0
        mass_fractions[elem] = np.nan_to_num(mass_fractions[elem], nan=0.0, posinf=0.0, neginf=0.0)
        mass_fractions[elem] = np.clip(mass_fractions[elem], 0.0, 1.0)
        # Apply floor: set mass fractions < 1e-8 to 0.0 to avoid numerical noise
        MASS_FRACTION_FLOOR = 1e-8
        mass_fractions[elem] = np.where(mass_fractions[elem] < MASS_FRACTION_FLOOR, 0.0, mass_fractions[elem])
    
    # Renormalize mass fractions to sum to exactly 1.0 at each depth
    mf_sum = np.zeros(len(r_new))
    for elem in mass_fractions:
        mf_sum += mass_fractions[elem]
    
    print(f"Mass fraction sum after interpolation: min={mf_sum.min():.6f}, max={mf_sum.max():.6f}")
    
    # Avoid division by zero - if sum is too small, use uniform distribution
    mf_sum_safe = np.where(mf_sum > 1e-10, mf_sum, 1.0)
    for elem in mass_fractions:
        mass_fractions[elem] /= mf_sum_safe
    
    # Handle case where sum was too small - preserve dominant species instead of uniform
    zero_mask = mf_sum < 1e-10
    if np.any(zero_mask):
        print(f"WARNING: Preserving dominant species at {np.sum(zero_mask)} zones where mass fraction sum was near zero")
        all_mf = np.stack([mass_fractions[elem] for elem in mass_fractions], axis=0)
        dominant_idx = np.argmax(all_mf, axis=0)
        elem_list = list(mass_fractions.keys())
        for i, elem in enumerate(elem_list):
            is_dominant = zero_mask & (dominant_idx == i)
            mass_fractions[elem] = np.where(is_dominant, 1.0, mass_fractions[elem])
            mass_fractions[elem] = np.where(zero_mask & (dominant_idx != i), 0.0, mass_fractions[elem])
    
    # Verify final normalization
    mf_sum_final = np.zeros(len(r_new))
    for elem in mass_fractions:
        mf_sum_final += mass_fractions[elem]
    print(f"Mass fraction sum after renormalization: min={mf_sum_final.min():.6f}, max={mf_sum_final.max():.6f}")
    
    # Recalculate kappa after regridding using updated composition and temperature
    X = mass_fractions.get('HYD', np.zeros(len(r_new)))
    Y = mass_fractions.get('HE', np.zeros(len(r_new)))
    
    # Calculate metal mass fraction Z (sum of all metals except H and He)
    Z = np.zeros(len(r_new))
    for elem in mass_fractions:
        if elem not in ['HYD', 'HE']:
            Z += mass_fractions[elem]
    
    # Calculate electron density with temperature-dependent ionization
    electron_density = calculate_electron_density(density_new, temperature_new, mass_fractions)
    
    # Calculate Rosseland mean opacity
    kappa_cm2_per_g_array = approximate_rosseland_opacity(temperature_new, density_new, X, Y, Z)
    
    # Update radius_cm to new grid in outer->inner order (consistent with interpolated arrays)
    radius_cm = r_new
    
    print(f"  After regridding:")
    print(f"    Velocity range: {velocity_kms.min():.2f} to {velocity_kms.max():.2f} km/s")
    print(f"    Density range: {density.min():.2e} to {density.max():.2e} g/cm^3")
    
    # Generate dr vs r plot comparing original and new grids
    import matplotlib.pyplot as plt
    # Original grid
    dr_orig = np.abs(np.diff(original_radius_cm))
    r_mid_orig = (original_radius_cm[:-1] + original_radius_cm[1:]) / 2
    # New grid
    dr_new = np.abs(np.diff(radius_cm))
    r_mid_new = (radius_cm[:-1] + radius_cm[1:]) / 2
    
    # Get density and temperature at mid-points for plotting
    rho_mid_new = (density[:-1] + density[1:]) / 2
    temp_mid_new = (temperature[:-1] + temperature[1:]) / 2
    
    fig, ax1 = plt.subplots(figsize=(12, 7))
    
    # Plot dr on left axis
    ax1.loglog(r_mid_orig, dr_orig, 'b--', alpha=0.5, linewidth=1, markersize=3, label=f'Original HERACLES ({len(original_radius_cm)} points)')
    ax1.loglog(r_mid_new, dr_new, 'r-', linewidth=1.5, markersize=4, label=f'AMR Grid ({len(radius_cm)} points)')
    ax1.set_xlabel('Radius (cm)', fontsize=12)
    ax1.set_ylabel('dr (cm)', fontsize=12, color='r')
    ax1.tick_params(axis='y', labelcolor='r')
    ax1.grid(True, which='both', linestyle='--', alpha=0.3)
    
    # Create first right axis for density
    ax2 = ax1.twinx()
    ax2.loglog(r_mid_new, rho_mid_new, 'g-', linewidth=1.5, alpha=0.7, label='Density')
    ax2.set_ylabel('Density (g/cm³)', fontsize=12, color='g')
    ax2.tick_params(axis='y', labelcolor='g')
    
    # Create second right axis for temperature (offset from first)
    ax3 = ax1.twinx()
    ax3.spines['right'].set_position(('axes', 1.15))
    ax3.loglog(r_mid_new, temp_mid_new, 'm-', linewidth=1.5, alpha=0.7, label='Temperature')
    ax3.set_ylabel('Temperature (K)', fontsize=12, color='m')
    ax3.tick_params(axis='y', labelcolor='m')
    
    # Create third right axis for radius (offset from second)
    ax4 = ax1.twinx()
    ax4.spines['right'].set_position(('axes', 1.3))
    ax4.loglog(r_mid_new, r_mid_new, 'c-', linewidth=1.5, alpha=0.7, label='Radius')
    ax4.set_ylabel('Radius (cm)', fontsize=12, color='c')
    ax4.tick_params(axis='y', labelcolor='c')
    
    # Combine legends from all axes
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    lines3, labels3 = ax3.get_legend_handles_labels()
    lines4, labels4 = ax4.get_legend_handles_labels()
    ax1.legend(lines1 + lines2 + lines3 + lines4, labels1 + labels2 + labels3 + labels4, 
              loc='upper right', fontsize=10)
    
    plt.title('Grid Spacing Comparison: Original HERACLES vs AMR\n(with Density, Temperature, and Radius profiles)', fontsize=14)
    # Set xlim to range of new AMR grid
    ax1.set_xlim(r_mid_new.min(), r_mid_new.max())
    plt.tight_layout()
    dr_plot_file = output_file.replace('.dat', '_dr_vs_r.png').replace('SN_HYDRO_DATA', 'SN_HYDRO_DATA_dr_vs_r')
    plt.savefig(dr_plot_file, dpi=150, bbox_inches='tight')
    print(f"  dr vs r plot saved to {dr_plot_file}")
    plt.close()
    
    # STEP 8: Convert velocity_cms -> velocity_kms (already done above)
    radius = radius_cm / 1e10  # 10^10 cm
    sigma_unclipped = calculate_sigma(radius, velocity_kms)
    sigma = sigma_unclipped  # Pass sigma values directly to CMFGEN
    
    # Apply Savitzky-Golay smoothing to velocity if requested
    if smooth_sigma:
        print(f"\n{'='*60}")
        print(f"Applying Sigma smoothing (Savitzky-Golay to velocity)...")
        print(f"{'='*60}")
        sigma, velocity_kms = smooth_sigma_savgol(radius, velocity_kms, window_length=sigma_window, 
                                                  polyorder=sigma_polyorder, threshold=sigma_threshold)
    
    # Apply density floor to prevent CMFGEN zero RGU error
    # CMFGEN fails when density is too low (E-89 or lower)
    DENSITY_MIN = 1e-22  # Minimum density floor in g/cm³
    density_before_floor = density.copy()
    density = np.maximum(density, DENSITY_MIN)
    
    # DISABLED: Distance-dependent temperature floor was crushing the temperature structure
    # The HERACLES temperature data should be preserved as-is
    # Original HERACLES already has reasonable temperature gradients
    print(f"  Distance-dependent temperature floor: DISABLED (preserving HERACLES structure)")
    print(f"  Temperature range: {temperature.min():.0f}K to {temperature.max():.0f}K")

    # STEP 8.5: Apply targeted temperature gradient fix
    # Preserves all HERACLES temperature structure outside flat regions
    # T_MIN set to 3000K, but system may go slightly below if needed to avoid flat regions
    temperature_before_fix = temperature.copy()
    temperature = apply_targeted_temperature_gradient_fix(temperature,
                                                          dtdr_threshold=1e-3,
                                                          min_change_pct=1e-3, temp_min=3000.0)
    if not np.array_equal(temperature, temperature_before_fix):
        print(f"  Temperature before fix: min={temperature_before_fix.min():.0f}K, max={temperature_before_fix.max():.0f}K")
        print(f"  Temperature after fix: min={temperature.min():.0f}K, max={temperature.max():.0f}K")

    # STEP 9: Recompute final tau and final sigma from the regridded arrays
    # Recalculate tau AFTER applying density floor for consistency
    tau = calculate_tau_from_outer(radius_cm, density, kappa_cm2_per_g_array)
    
    print(f"  Final optical depth range: {tau.min():.2e} to {tau.max():.2e}")
    
    # STEP 9.5: Final strict gradient enforcement pass
    # Runs after all post-processing (sigma fix, T fix, outer points) to guarantee
    # no violation survives into the output file.
    print("\n[Final Gradient Enforcement] Checking for residual violations...")
    velocity_cms_final = velocity_kms * 1e5
    final_violation_mask = detect_gradient_violations(
        radius_cm, density, velocity_cms_final, temperature,
        max_drho_pct=MAX_DRHO_PCT, max_dt_pct=MAX_DT_PCT, max_dv_pct=MAX_DV_PCT
    )
    n_final_violations = int(np.sum(final_violation_mask))
    if n_final_violations > 0:
        print(f"[Final Gradient Enforcement] {n_final_violations} violations — applying monotone log-linear blending...")
        for _pass in range(50):
            final_violation_mask = detect_gradient_violations(
                radius_cm, density, velocity_cms_final, temperature,
                max_drho_pct=MAX_DRHO_PCT, max_dt_pct=MAX_DT_PCT, max_dv_pct=MAX_DV_PCT
            )
            if not np.any(final_violation_mask):
                print(f"[Final Gradient Enforcement] All violations resolved after {_pass} passes")
                break
            density, velocity_cms_final, temperature = smooth_violations(
                radius_cm, density, velocity_cms_final, temperature, final_violation_mask
            )
        else:
            n_remaining = int(np.sum(detect_gradient_violations(
                radius_cm, density, velocity_cms_final, temperature,
                max_drho_pct=MAX_DRHO_PCT, max_dt_pct=MAX_DT_PCT, max_dv_pct=MAX_DV_PCT
            )))
            print(f"[Final Gradient Enforcement] WARNING: {n_remaining} violations remain after 50 passes")
        velocity_kms = velocity_cms_final / 1e5
        # Recompute tau and sigma after final enforcement
        tau = calculate_tau_from_outer(radius_cm, density, kappa_cm2_per_g_array)
    else:
        print("[Final Gradient Enforcement] No violations — output is clean.")

    # STEP 10: Validate lengths
    n_points = len(radius_cm)
    print(f"\nValidating final model (n_points={n_points})...")
    validate_final_model(radius_cm, density, velocity_kms, temperature, tau, sigma, 
                       kappa_cm2_per_g_array, mass_fractions)
    print("  All arrays have consistent length")
    
    # STEP 11: Write output (using only final regridded arrays)
    # Set up output parameters
    n_mass_fractions = 8  # HYD, HE, CARB, OXY, SIL, IRON, COB, NICK
    n_isotopes = 0
    
    # Ensure positive time for CMFGEN/SEDONA compatibility
    # HERACLES may have negative times for pre-explosion phases
    if actual_days < 0:
        print(f"[Time Correction] Negative time detected: {actual_days:.2f} days")
        print(f"  Converting to positive time since explosion bomb introduction")
        time_days = abs(actual_days)
    else:
        time_days = actual_days
    
    energy_lum = 6.071478e6
    # Cap temperature at CMFGEN-safe maximum before unit conversion.
    # Post-shock reverse-shock gas can reach ~1.5e6 K in HERACLES, which is beyond
    # the range of CMFGEN's LTE opacity tables and causes GFF_VEC to receive
    # extreme values that lead to a SIGSEGV during the LTE cold-start.
    T_MAX_CMFGEN = 5.0e5  # 50 x 10^4 K
    n_capped = int(np.sum(temperature > T_MAX_CMFGEN))
    if n_capped > 0:
        print(f"[T cap] Capping {n_capped} depths from max={temperature.max():.2e} K to T_MAX={T_MAX_CMFGEN:.2e} K")
        temperature = np.minimum(temperature, T_MAX_CMFGEN)
    # Ensure floor before breaking isothermal regions
    t_floor = 3000.0
    n_floored = int(np.sum(temperature < t_floor))
    if n_floored > 0:
        print(f"[T floor] Flooring {n_floored} depths to T_FLOOR={t_floor:.0f} K")
        temperature = np.maximum(temperature, t_floor)
    # Break flat isothermal regions after floor to prevent SEDONA/CMFGEN singular-matrix crashes
    temperature = break_isothermal_regions(temperature, t_floor=t_floor)
    temperature_10k = temperature / 1e4
    
    # Calculate atom density and electron density on regridded arrays
    # Using mean atomic mass from mass fractions (fully ionised approximation)
    atomic_masses = {'HYD': 1.0, 'HE': 4.0, 'CARB': 12.0, 'N': 14.0, 'OXY': 16.0,
                     'NE': 20.0, 'MG': 24.0, 'SIL': 28.0, 'S': 32.0, 'IRON': 56.0}
    atomic_Z      = {'HYD': 1,   'HE': 2,   'CARB': 6,    'N': 7,    'OXY': 8,
                     'NE': 10,  'MG': 12,  'SIL': 14,  'S': 16,   'IRON': 26}
    
    mean_A  = np.zeros(n_points)
    mean_Z  = np.zeros(n_points)
    for elem, A in atomic_masses.items():
        if elem in mass_fractions:
            Xi = mass_fractions[elem]
            mean_A += Xi / A          # number fraction contribution
            mean_Z += Xi / A * atomic_Z[elem]
    
    # mean_A is now sum(Xi/Ai) — the mean number of particles per amu
    # atom_density = rho * sum(Xi/Ai) * N_Avogadro
    N_A = 6.02214076e23  # mol^-1
    atom_density     = density * mean_A * N_A
    # Calculate electron density with temperature-dependent ionization
    electron_density = calculate_electron_density(density, temperature, mass_fractions)
    
    # Write output file in CMFGEN format
    print(f"Writing CMFGEN format to {output_file}...")
    
    # Verify arrays are in outer->inner order (CMFGEN expects this)
    # CMFGEN expects outer->inner order (decreasing radius)
    # Arrays should already be in outer->inner order from interpolation step
    if radius[0] < radius[-1]:
        print("WARNING: radius array is in inner->outer order, reversing to outer->inner for CMFGEN")
        radius = radius[::-1]
        velocity_kms = velocity_kms[::-1]
        sigma = sigma[::-1]
        temperature_10k = temperature_10k[::-1]
        density = density[::-1]
        electron_density = electron_density[::-1]
        atom_density = atom_density[::-1]
        kappa_cm2_per_g_array = kappa_cm2_per_g_array[::-1]
        for elem in mass_fractions:
            mass_fractions[elem] = mass_fractions[elem][::-1]
    
    # Verify orientation is correct
    assert radius[0] > radius[-1], "ERROR: radius array should be in outer->inner order (decreasing)"
    print(f"Verified: radius[0]={radius[0]:.2e} > radius[-1]={radius[-1]:.2e} (outer->inner order)")
    
    # Note: Velocity trimming is already applied earlier in the conversion function (STEP 5)
    # No need to apply it again here

    # Remove near-degenerate depth points (dr/r < 0.5%) that cause singular BA matrices
    # in CMF_BLKBAND_V3, leading to IFAIL errors and SIGSEGV. This targets flat plateau
    # zones (e.g. post-shock ejecta with near-constant velocity/density over many cells).
    print(f"\nApplying grid thinning to remove near-degenerate points...")
    mf_arrays = [mass_fractions[elem] for elem in mass_fractions]
    result = thin_degenerate_grid(
        radius, velocity_kms, density, temperature_10k, electron_density,
        atom_density, sigma, kappa_cm2_per_g_array, *mf_arrays
    )
    n_fixed = 8  # number of named arrays before mass fractions
    radius, velocity_kms, density, temperature_10k, electron_density, atom_density, sigma, kappa_cm2_per_g_array = result[:n_fixed]
    for i, elem in enumerate(mass_fractions):
        mass_fractions[elem] = result[n_fixed + i]
    n_points = len(radius)

    # Ensure radius grid is strictly monotonic before writing SN_HYDRO_DATA
    # CMFGEN requires strictly monotonic radius grid (no duplicate values)
    radius, velocity_kms, density, temperature_10k, electron_density, atom_density, sigma, kappa_cm2_per_g_array = ensure_monotonic_grid(
        radius, velocity_kms, density, temperature_10k, electron_density, atom_density, sigma, kappa_cm2_per_g_array
    )

    with open(output_file, 'w') as f:
        # Header (exact format matching working CMFGEN SN_HYDRO_DATA)
        f.write(f"Number of data points:               {n_points:>3d}\n")
        f.write(f"Number of mass fractions:            {n_mass_fractions:>2d}\n")
        f.write(f"Number of isotopes:                  {n_isotopes:>3d}\n")
        f.write(f"Time(days) since explosion:          {time_days:10.5E}\n")
        f.write(f"Energy Lum [lsun]:                   {energy_lum:11.6E}\n")
        f.write(f"Inner boundary radius [cm]:          {r_inner_boundary_cm:15.6E}\n")
        f.write(f"Inner boundary temperature [K]:      {T_inner_boundary_K:15.6E}\n")
        f.write("\n")
        
        # Radius grid (10^10 cm)
        f.write("Radius grid (10^10cm)\n")
        write_array_8_per_line(f, radius)
        
        # Velocity (km/s)
        f.write("Velocity (km/s)\n")
        write_array_8_per_line(f, velocity_kms)
        
        # Sigma
        f.write("Sigma (dlnV/dlnr-1)\n")
        write_array_8_per_line(f, sigma)
        
        # Temperature (10^4 K)
        f.write("Temperature (10^4 K)\n")
        write_array_8_per_line(f, temperature_10k)
        
        # Density (gm/cm^3)
        f.write("Density (gm/cm^3)\n")
        write_array_8_per_line(f, density)
        
        # Electron density (1/cm^3)
        f.write("Electron density (1/cm^3)\n")
        write_array_8_per_line(f, electron_density)
        
        # Atom density (1/cm^3)
        f.write("Atom density (1/cm^3)\n")
        write_array_8_per_line(f, atom_density)
        
        # Kappa (cm^2/gm)
        f.write("Kappa (cm^2/gm)\n")
        write_array_8_per_line(f, kappa_cm2_per_g_array)
        
        # Mass fractions - already in CMFGEN canonical names from ingest
        # Output 8 mass fractions: HYD, HE, CARB, OXY, SIL, IRON, COB, NICK
        all_mass_fractions = {}
        canonical_elems = ['HYD', 'HE', 'CARB', 'OXY', 'SIL', 'IRON']
        for elem in canonical_elems:
            if elem in mass_fractions:
                all_mass_fractions[elem] = mass_fractions[elem]
            else:
                all_mass_fractions[elem] = np.zeros(n_points)
        
        # COB and NICK are zero (not tracked in HERACLES)
        all_mass_fractions['COB'] = np.zeros(n_points)
        all_mass_fractions['NICK'] = np.zeros(n_points)
        
        # Normalize all 8 mass fractions to sum to exactly 1 at each depth
        
        all_elements = ['HYD', 'HE', 'CARB', 'OXY', 'SIL', 'IRON', 'COB', 'NICK']
        mf_sum = np.zeros(n_points)
        for elem in all_elements:
            if elem in all_mass_fractions:
                mf_sum += all_mass_fractions[elem]
        
        print(f"Mass fraction sum before normalization (all elements): min={mf_sum.min():.6f}, max={mf_sum.max():.6f}")
        
        # Avoid division by zero
        mf_sum_safe = np.where(mf_sum > 1e-10, mf_sum, 1.0)
        for elem in all_elements:
            if elem in all_mass_fractions:
                all_mass_fractions[elem] = all_mass_fractions[elem] / mf_sum_safe
        
        # Verify normalization
        mf_sum_after = np.zeros(n_points)
        for elem in all_elements:
            if elem in all_mass_fractions:
                mf_sum_after += all_mass_fractions[elem]
        print(f"Mass fraction sum after normalization (all elements): min={mf_sum_after.min():.6f}, max={mf_sum_after.max():.6f}")
        
        # If total mass fraction sum < 1e-8 at any depth, set to zero (negligible material)
        # Use post-normalization sum to check for negligible material (density-based check would be better)
        zero_mask = mf_sum_after < 0.5  # If sum is significantly less than 1 after normalization, set to zero
        if np.any(zero_mask):
            print(f"Setting mass fractions to zero at {np.sum(zero_mask)} depths where total < 0.5")
        
        # Mass fractions are already in outer->inner order from interpolation step
        # No reversal needed - they are consistent with radius array
        
        # Write element mass fractions (normalized without artificial floor)
        # MINIMUM_ISO_POP in CMFGEN's rd_sn_data.f floors isotope populations to 1e-20.
        # The element consistency check then fails if the element mass fraction is < 1e-20
        # while the isotope is floored to 1e-20.  Apply the same floor to isotope-bearing
        # element species so element and isotope stay consistent after CMFGEN's flooring.
        CMFGEN_MIN_ISO_POP = 1.0e-20
        ISO_BEARING = {'IRON', 'COB', 'NICK'}
        normalized_mf = {}
        for elem in ['HYD', 'HE', 'CARB', 'OXY', 'SIL', 'IRON', 'COB', 'NICK']:
            if elem in all_mass_fractions:
                f.write(f"{elem} mass fraction\n")
                # Normalize mass fractions to sum to 1
                mf_normalized = all_mass_fractions[elem]
                # Set to zero where total mass fraction sum < 1e-8
                mf_normalized = np.where(zero_mask, 0.0, mf_normalized)
                # For isotope-bearing species, floor to CMFGEN's MINIMUM_ISO_POP so that
                # the element mass fraction never falls below the isotope floor
                if elem in ISO_BEARING:
                    mf_normalized = np.where(
                        (mf_normalized > 0) & (mf_normalized < CMFGEN_MIN_ISO_POP),
                        CMFGEN_MIN_ISO_POP,
                        mf_normalized,
                    )
                write_array_8_per_line(f, mf_normalized)
                normalized_mf[elem] = mf_normalized
        
        # Isotope mass fractions (only write if n_isotopes > 0)
        # Set isotope mass fractions to match element mass fractions exactly
        # This prevents CMFGEN from renormalizing and causing mismatches
        if n_isotopes > 0:
            f.write("NICK 56 mass fraction\n")
            write_array_8_per_line(f, normalized_mf['NICK'][::-1])
            f.write("COB 56 mass fraction\n")
            write_array_8_per_line(f, normalized_mf['COB'][::-1])
            f.write("IRON 56 mass fraction\n")
            # Use exact same normalized value as IRON element
            write_array_8_per_line(f, normalized_mf['IRON'][::-1])
    
    # Create stacked mass fraction plot vs enclosed mass
    print("Creating mass fraction stacked plot...")
    # Use original full profile data to show complete 13 Msun profile
    plot_mass_fractions_stacked(original_radius_cm, original_density, original_mass_fractions, time_days, output_file)
    
    print(f"Conversion complete: {output_file}")
    print(f"  Data points: {n_points}")
    print(f"  Time: {time_days:.1f} days")
    print(f"  Radius range: {radius.min():.2e} to {radius.max():.2e} (10^10 cm)")
    print(f"  Velocity range: {velocity_kms.min():.2f} to {velocity_kms.max():.2f} km/s")
    print(f"  Sigma range: {sigma.min():.2e} to {sigma.max():.2e}")
    print(f"  Sigma clipping: DISABLED (passing values directly to CMFGEN)")
    print(f"  Density range: {density.min():.2e} to {density.max():.2e} g/cm^3")
    print(f"  Temperature range: {temperature.min():.2e} to {temperature.max():.2e} K")
    print(f"  Optical depth range: {tau.min():.2e} to {tau.max():.2e}")
    
    # Plot before/after comparison
    # Read the SN_HYDRO_DATA file directly to show exactly what CMFGEN will read
    print("Reading SN_HYDRO_DATA file directly for comparison plot...")
    cmfgen_data = read_cmfgen_file(output_file)
    
    # mass_fractions already uses canonical CMFGEN names from ingest
    cmfgen_mass_fractions = {elem: mass_fractions[elem] for elem in mass_fractions}
    
    # Plot comparison using data directly from SN_HYDRO_DATA file
    # SN_HYDRO_DATA is now written in inner->outer order (increasing radius)
    # Original HERACLES data is in outer->inner order, so reverse it to inner->outer for comparison
    plot_comparison(original_radius_cm[::-1], original_density[::-1], original_temperature[::-1], 
                   original_velocity_kms[::-1], original_tau[::-1], original_sigma_unclipped[::-1],
                   cmfgen_data['radius'] * 1e10,  # Convert back to cm
                   cmfgen_data['density'],
                   cmfgen_data['temperature'] * 1e4,  # Convert back to K
                   cmfgen_data['velocity'],
                   cmfgen_data.get('tau', np.zeros_like(cmfgen_data['radius'])),  # Tau not in file
                   cmfgen_data['sigma'],
                   cmfgen_data['sigma'],  # sigma_unclipped same as sigma
                   output_file, time_days, 
                   orig_mass_fractions={k: v[::-1] for k, v in original_mass_fractions.items()}, 
                   new_mass_fractions=cmfgen_mass_fractions,
                   velocity_trimmed_indices=velocity_trimmed_indices)

def plot_comparison(orig_r, orig_rho, orig_T, orig_v, orig_tau, orig_sigma,
                   new_r, new_rho, new_T, new_v, new_tau, new_sigma, new_sigma_unclipped,
                   output_file, time_days, orig_mass_fractions=None, new_mass_fractions=None,
                   velocity_trimmed_indices=None):
    """Plot before/after comparison of HERACLES to CMFGEN conversion"""
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon
    
    # Solar radius in cm
    rsun_cm = 6.96e10
    
    # Convert radius to r/rsun for all subplots (actual values, not log)
    orig_r_rsun = orig_r / rsun_cm
    new_r_rsun = new_r / rsun_cm

    # Add grey vertical span for velocity-trimmed material if indices provided
    velocity_trimmed_rsun = None
    if velocity_trimmed_indices is not None and len(velocity_trimmed_indices) > 0:
        # Convert trimmed indices to radius range
        velocity_trimmed_rsun = orig_r_rsun[velocity_trimmed_indices]
        if len(velocity_trimmed_rsun) > 0:
            trim_min = velocity_trimmed_rsun.min()
            trim_max = velocity_trimmed_rsun.max()
            print(f"[Plot] Adding grey vertical span for velocity-trimmed material: {trim_min:.2f} to {trim_max:.2f} R☉")

    # Note: Radius arrays may be in different orders; plotting handles this correctly
    
    # Create a 2x3 grid for the main comparison plots
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle(f'HERACLES to CMFGEN Conversion: {time_days}-day snapshot', fontsize=16)
    
    # Plot 1: Density vs radius (both axes log scale)
    ax = axes[0, 0]
    ax.loglog(orig_r_rsun, orig_rho, 'b--', alpha=0.7, linewidth=2, label='Original HERACLES (input grid)')
    ax.loglog(new_r_rsun, new_rho, 'r-', alpha=0.7, linewidth=2, label='CMFGEN (outline)')
    # Add grey vertical span for velocity-trimmed material
    if velocity_trimmed_rsun is not None:
        ax.axvspan(trim_min, trim_max, facecolor='grey', alpha=0.3, label='Velocity-trimmed')
    ax.set_xscale('log')
    ax.set_xlabel('r/R☉', fontsize=12)
    ax.set_ylabel('Density (g/cm³)', fontsize=12)
    ax.set_title('Density Profile', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(left=10.0)  # Minimum at 10 Rsun
    
    # Plot 2: Temperature vs radius (both axes log scale)
    ax = axes[0, 1]
    ax.loglog(orig_r_rsun, orig_T, 'b--', alpha=0.7, linewidth=2, label='Original HERACLES (input grid)')
    ax.loglog(new_r_rsun, new_T, 'r-', alpha=0.7, linewidth=2, label='CMFGEN (outline)')
    # Add grey vertical span for velocity-trimmed material
    if velocity_trimmed_rsun is not None:
        ax.axvspan(trim_min, trim_max, facecolor='grey', alpha=0.3)
    ax.set_xscale('log')
    ax.set_xlabel('r/R☉', fontsize=12)
    ax.set_ylabel('Temperature (K)', fontsize=12)
    ax.set_title('Temperature Profile', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(left=10.0)  # Minimum at 10 Rsun
    
    # Plot 3: Velocity vs radius (both axes log scale)
    ax = axes[0, 2]
    ax.loglog(orig_r_rsun, orig_v, 'b--', alpha=0.7, linewidth=2, label='Original HERACLES (input grid)')
    ax.loglog(new_r_rsun, new_v, 'r-', alpha=0.7, linewidth=2, label='CMFGEN (outline)')
    # Add grey vertical span for velocity-trimmed material
    if velocity_trimmed_rsun is not None:
        ax.axvspan(trim_min, trim_max, facecolor='grey', alpha=0.3)
    ax.set_xscale('log')
    ax.set_xlabel('r/R☉', fontsize=12)
    ax.set_ylabel('Velocity (km/s)', fontsize=12)
    ax.set_title('Velocity Profile', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(left=10.0)  # Minimum at 10 Rsun
    
    # Plot 4: Sigma (unclipped) vs radius - use symlog for negative values
    ax = axes[1, 0]
    ax.step(orig_r_rsun, orig_sigma, 'b--', alpha=0.7, linewidth=2, label='Original HERACLES (input grid)', where='post')
    ax.step(new_r_rsun, new_sigma_unclipped, 'r-', alpha=0.7, linewidth=2, label='CMFGEN (outline)', where='post')
    # Add grey vertical span for velocity-trimmed material
    if velocity_trimmed_rsun is not None:
        ax.axvspan(trim_min, trim_max, facecolor='grey', alpha=0.3)
    ax.axhline(y=-0.99, color='k', linestyle=':', alpha=0.5, label='σ=-0.99 (clipping limit)')
    ax.axhline(y=-1.0, color='gray', linestyle='--', alpha=0.5, label='σ=-1.0 (CMFGEN limit)')
    ax.set_xscale('log')
    ax.set_xlabel('r/R☉', fontsize=12)
    ax.set_ylabel('Sigma (σ)', fontsize=12)
    ax.set_title('Sigma Profile', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_yscale('symlog', linthresh=1.0)
    ax.set_xlim(left=10.0)  # Minimum at 10 Rsun
    
    # Plot 5: Optical depth vs radius
    ax = axes[1, 1]
    ax.step(orig_r_rsun, orig_tau, 'b--', alpha=0.7, linewidth=2, label='Original HERACLES (input grid)', where='post')
    ax.step(new_r_rsun, new_tau, 'r-', alpha=0.7, linewidth=2, label='CMFGEN (outline)', where='post')
    ax.axhline(y=2/3, color='orange', linestyle=':', alpha=0.8, linewidth=2, label='τ = 2/3 (photosphere)')
    # Add horizontal shaded region for tau_min,max (red)
    tau_min, tau_max = np.min(new_tau), np.max(new_tau)
    ax.axhspan(tau_min, tau_max, facecolor='red', alpha=0.2, label=f'τ range: {tau_min:.2e} to {tau_max:.2e}')
    # Add grey vertical span for velocity-trimmed material
    if velocity_trimmed_rsun is not None:
        ax.axvspan(trim_min, trim_max, facecolor='grey', alpha=0.3)
    ax.set_xscale('log')
    ax.set_xlabel('r/R☉', fontsize=12)
    ax.set_ylabel('Optical Depth τ', fontsize=12)
    ax.set_title('Optical Depth Profile', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_yscale('log')
    ax.set_xlim(left=10.0)  # Minimum at 10 Rsun
    
    # Plot 6: Mass fractions vs radius (if available)
    ax = axes[1, 2]
    if orig_mass_fractions is not None and new_mass_fractions is not None:
        # Both dictionaries now use canonical CMFGEN names
        canonical_elems = ['HYD', 'HE', 'CARB', 'OXY', 'SIL', 'IRON']
        colors = ['blue', 'green', 'red', 'orange', 'purple', 'brown']
        
        # Plot each element with single label
        for i, elem in enumerate(canonical_elems):
            if elem in orig_mass_fractions and elem in new_mass_fractions:
                # Plot input grid (dashed, no label)
                ax.step(orig_r_rsun, orig_mass_fractions[elem], '--', color=colors[i], alpha=0.5, linewidth=1, where='post')
                # Plot CMFGEN (solid, with label)
                ax.step(new_r_rsun, new_mass_fractions[elem], '-', color=colors[i], alpha=0.8, linewidth=2, label=elem, where='post')
        
        # Add grey vertical span for velocity-trimmed material
        if velocity_trimmed_rsun is not None:
            ax.axvspan(trim_min, trim_max, facecolor='grey', alpha=0.3)
        ax.set_xscale('log')
        ax.set_xlabel('r/R☉', fontsize=12)
        ax.set_ylabel('Mass Fraction', fontsize=12)
        ax.set_title('Mass Fractions', fontsize=14)
        ax.legend(fontsize=8, loc='upper right')
        ax.grid(True, alpha=0.3)
        ax.set_xlim(left=10.0)  # Minimum at 10 Rsun
    else:
        ax.text(0.5, 0.5, 'Mass fractions not available', 
                ha='center', va='center', transform=ax.transAxes, fontsize=14)
        ax.set_title('Mass Fractions', fontsize=14)
        ax.set_title('Mass Fraction Profile', fontsize=14)
    
    plt.tight_layout()
    
    # Save plot
    plot_file = output_file.replace('.dat', '_comparison.png').replace('SN_HYDRO_DATA', 'SN_HYDRO_DATA_comparison')
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    print(f"Plot saved to {plot_file}")
    plt.close()

def plot_mass_fractions_stacked(radius_cm, density, mass_fractions, time_days, output_file):
    """
    Create a stacked area plot of mass fractions vs enclosed mass.
    
    Args:
        radius_cm: Radius array in cm (can be either order, will be converted to inner->outer)
        density: Density array in g/cm^3
        mass_fractions: Dict of element mass fractions (e.g., {'HYD': array, 'HE': array, ...})
        time_days: Time since explosion in days
        output_file: Path to SN_HYDRO_DATA file (used to determine output plot path)
    """
    import matplotlib.pyplot as plt
    
    # Convert to inner->outer order (center to surface) for enclosed mass calculation
    r = np.asarray(radius_cm)
    rho = np.asarray(density)
    
    # Check if array is in outer->inner order (decreasing radius)
    if r[0] > r[-1]:
        # Reverse to inner->outer order
        r = r[::-1]
        rho = rho[::-1]
        for elem in mass_fractions:
            mass_fractions[elem] = mass_fractions[elem][::-1]
    
    # Calculate enclosed mass from center to each radius
    # M(<r) = ∫_0^r 4π r'^2 ρ(r') dr'
    
    # Integrate from center outward to get enclosed mass at each radius
    integrand = 4.0 * np.pi * (r ** 2) * rho  # g/cm
    if hasattr(np, "trapezoid"):
        enclosed_mass = np.trapezoid(integrand, r, axis=0)
    else:
        enclosed_mass = np.trapz(integrand, r)
    
    # Cumulative sum to get mass at each shell (more accurate for discrete data)
    dr = np.diff(r, prepend=r[0])
    dr[0] = r[0]  # First shell has radius r[0]
    enclosed_mass_cum = np.cumsum(4.0 * np.pi * (r ** 2) * rho * dr)
    
    # Convert to solar masses for log scale plot
    msun_g = 1.98847e33
    enclosed_mass_msun = enclosed_mass_cum / msun_g
    
    # Prepare mass fraction data - order elements by atomic mass (roughly)
    element_order = ['HYD', 'HE', 'CARB', 'OXY', 'SIL', 'IRON', 'COB', 'NICK']
    element_colors = {
        'HYD': '#FF6B6B',   # Light red
        'HE': '#4ECDC4',     # Teal
        'CARB': '#45B7D1',   # Blue
        'OXY': '#96CEB4',    # Green
        'SIL': '#FFEAA7',    # Yellow
        'IRON': '#DDA0DD',   # Plum
        'COB': '#98D8C8',    # Mint
        'NICK': '#F7DC6F'    # Gold
    }
    element_labels = {
        'HYD': 'H',
        'HE': 'He',
        'CARB': 'C',
        'OXY': 'O',
        'SIL': 'Si',
        'IRON': 'Fe',
        'COB': 'Co',
        'NICK': 'Ni'
    }
    
    # Stack mass fractions
    mf_arrays = []
    labels = []
    colors = []
    
    for elem in element_order:
        if elem in mass_fractions and np.any(mass_fractions[elem] > 1e-30):
            mf_arrays.append(mass_fractions[elem])
            labels.append(element_labels.get(elem, elem))
            colors.append(element_colors.get(elem, 'gray'))
    
    if not mf_arrays:
        print("No mass fraction data to plot")
        return
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Create stacked area plot
    x = enclosed_mass_msun
    y = np.vstack(mf_arrays)
    
    ax.stackplot(x, y, labels=labels, colors=colors, alpha=0.8)
    
    # Labels and formatting
    ax.set_xlabel('Enclosed Mass (M$_\odot$)', fontsize=12)
    ax.set_ylabel('Mass Fraction', fontsize=12)
    ax.set_title(f'Mass Fraction Composition vs Enclosed Mass\nt = {time_days:.1f} days', fontsize=14)
    
    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(0, 1)
    
    # Legend
    ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), fontsize=10)
    
    ax.grid(True, alpha=0.3, linestyle='--')
    
    plt.tight_layout()
    
    # Save plot
    plot_file = output_file.replace('.dat', '_mass_fraction.png').replace('SN_HYDRO_DATA', 'SN_HYDRO_DATA_mass_fraction')
    if plot_file == output_file:
        plot_file = str(Path(output_file).with_suffix('')) + '_mass_fraction.png'
    
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    print(f"Mass fraction plot saved to {plot_file}")
    plt.close()

def write_array_8_per_line(f, array, format_str='%.15E'):
    """Write array with 8 values per line in scientific notation"""
    n = len(array)
    for i in range(0, n, 8):
        chunk = array[i:i+8]
        line = '    '.join([format_str % val for val in chunk])
        f.write(line + '\n')
    f.write('\n')

def read_cmfgen_file(filename):
    """Read the CMFGEN SN_HYDRO_DATA file and return data dictionary"""
    # Read the file
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # Parse header
    n_points = int(lines[0].split()[-1])
    time_days_raw = float(lines[3].split()[-1])
    energy_lum = float(lines[4].split()[-1])
    
    # Ensure positive time for consistency
    if time_days_raw < 0:
        time_days = abs(time_days_raw)
    else:
        time_days = time_days_raw
    
    # Find data sections
    data_sections = {}
    current_section = None
    data_lines = []
    
    for line in lines[6:]:  # Skip header
        line = line.strip()
        if not line:
            if current_section and data_lines:
                data_sections[current_section] = data_lines
                data_lines = []
            continue
        
        # Check if this is a section header
        # Section headers either have units in parentheses, or contain "mass fraction"
        is_section_header = False
        if '(' in line and ')' in line and any(unit in line for unit in ['cm', 'km/s', 'K', 'gm', '1/cm', 'cm^2']):
            is_section_header = True
        elif 'mass fraction' in line.lower():
            is_section_header = True
        
        if is_section_header:
            if current_section and data_lines:
                data_sections[current_section] = data_lines
                data_lines = []
            current_section = line
            continue
        
        if current_section:
            data_lines.append(line)
    
    # Parse data from each section
    data = {}
    for section, lines_list in data_sections.items():
        values = []
        for line in lines_list:
            values.extend([float(x) for x in line.split()])
        data[section] = np.array(values)
    
    # Check for size mismatches
    sizes = {section: len(values) for section, values in data.items()}
    unique_sizes = set(sizes.values())
    if len(unique_sizes) > 1:
        # Use the most common size (likely the correct one)
        most_common_size = max(set(sizes.values()), key=list(sizes.values()).count)
        # Truncate or pad arrays to match
        for section in data:
            if len(data[section]) != most_common_size:
                if len(data[section]) > most_common_size:
                    data[section] = data[section][:most_common_size]
                else:
                    # Pad with zeros
                    pad_length = most_common_size - len(data[section])
                    data[section] = np.pad(data[section], (0, pad_length), 'constant')
    
    # Extract key arrays
    result = {
        'radius': data.get('Radius grid (10^10cm)', np.zeros(n_points)),
        'velocity': data.get('Velocity (km/s)', np.zeros(n_points)),
        'sigma': data.get('Sigma (dlnV/dlnr-1)', np.zeros(n_points)),
        'temperature': data.get('Temperature (10^4 K)', np.zeros(n_points)),
        'density': data.get('Density (gm/cm^3)', np.zeros(n_points)),
        'electron_density': data.get('Electron density (1/cm^3)', np.zeros(n_points)),
        'atom_density': data.get('Atom density (1/cm^3)', np.zeros(n_points)),
        'kappa': data.get('Kappa (cm^2/gm)', np.zeros(n_points)),
    }
    
    # Extract mass fractions
    for section in data:
        if 'mass fraction' in section.lower():
            elem_name = section.split()[0]  # First word is element name
            result[elem_name] = data[section]
    
    return result

def plot_cmfgen_file(filename):
    """Read and plot the CMFGEN SN_HYDRO_DATA file"""
    import matplotlib.pyplot as plt
    
    # Read the file
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # Parse header
    n_points = int(lines[0].split()[-1])
    time_days = float(lines[3].split()[-1])
    energy_lum = float(lines[4].split()[-1])
    
    print(f"File info:")
    print(f"  Data points: {n_points}")
    print(f"  Time: {time_days} days")
    print(f"  Energy Lum: {energy_lum} lsun")
    
    # Find data sections
    data_sections = {}
    current_section = None
    data_lines = []
    
    for line in lines[6:]:  # Skip header
        line = line.strip()
        if not line:
            if current_section and data_lines:
                data_sections[current_section] = data_lines
                data_lines = []
            continue
        
        # Check if this is a section header
        # Section headers either have units in parentheses, or contain "mass fraction"
        is_section_header = False
        if '(' in line and ')' in line and any(unit in line for unit in ['cm', 'km/s', 'K', 'gm', '1/cm', 'cm^2']):
            is_section_header = True
        elif 'mass fraction' in line.lower():
            is_section_header = True
        
        if is_section_header:
            if current_section and data_lines:
                data_sections[current_section] = data_lines
                data_lines = []
            current_section = line
            continue
        
        if current_section:
            data_lines.append(line)
    
    # Parse data from each section
    data = {}
    for section, lines_list in data_sections.items():
        values = []
        for line in lines_list:
            values.extend([float(x) for x in line.split()])
        data[section] = np.array(values)
        print(f"  Section '{section}': {len(values)} points")
    
    # Check for size mismatches
    sizes = {section: len(values) for section, values in data.items()}
    unique_sizes = set(sizes.values())
    if len(unique_sizes) > 1:
        print(f"WARNING: Size mismatch detected!")
        for section, size in sizes.items():
            print(f"  {section}: {size} points")
        # Use the most common size (likely the correct one)
        most_common_size = max(set(sizes.values()), key=list(sizes.values()).count)
        print(f"Using most common size: {most_common_size}")
        # Truncate or pad arrays to match
        for section in data:
            if len(data[section]) != most_common_size:
                if len(data[section]) > most_common_size:
                    data[section] = data[section][:most_common_size]
                    print(f"  Truncated {section} to {most_common_size} points")
                else:
                    # Pad with zeros
                    pad_length = most_common_size - len(data[section])
                    data[section] = np.pad(data[section], (0, pad_length), 'constant')
                    print(f"  Padded {section} to {most_common_size} points")
    
    # Create plots
    fig, axes = plt.subplots(3, 3, figsize=(18, 18))
    fig.suptitle(f'CMFGEN SN_HYDRO_DATA: {time_days} days', fontsize=16)
    
    # SN_HYDRO_DATA is written in outer->inner order for CMFGEN
    # Reverse to inner->outer for plotting
    if 'Radius grid (10^10cm)' in data:
        radius = data['Radius grid (10^10cm)'][::-1] * 1e10  # Convert back to cm and reverse
    else:
        radius = None
    
    # Plot 1: Radius vs Velocity
    if 'Velocity (km/s)' in data and radius is not None:
        ax = axes[0, 0]
        velocity = data['Velocity (km/s)'][::-1]  # Reverse to match radius
        ax.step(radius, velocity, 'b-', linewidth=2, where='post')
        ax.set_xlabel('Radius (cm)')
        ax.set_ylabel('Velocity (km/s)')
        ax.set_title('Radius vs Velocity')
        ax.set_xscale('log')
        ax.grid(True, alpha=0.3)
    
    # Plot 2: Density profile
    if 'Density (gm/cm^3)' in data and radius is not None:
        ax = axes[0, 1]
        density = data['Density (gm/cm^3)'][::-1]  # Reverse to match radius
        ax.step(radius, density, 'r-', linewidth=2, where='post')
        ax.set_xlabel('Radius (cm)')
        ax.set_ylabel('Density (g/cm³)')
        ax.set_title('Density Profile')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.grid(True, alpha=0.3)
    
    # Plot 3: Temperature profile
    if 'Temperature (10^4 K)' in data and radius is not None:
        ax = axes[0, 2]
        temp = data['Temperature (10^4 K)'][::-1] * 1e4  # Reverse and convert back to K
        ax.step(radius, temp, 'g-', linewidth=2, where='post')
        ax.set_xlabel('Radius (cm)')
        ax.set_ylabel('Temperature (K)')
        ax.set_title('Temperature Profile')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.grid(True, alpha=0.3)
    
    # Plot 4: Sigma
    if 'Sigma (dlnV/dlnr-1)' in data and radius is not None:
        ax = axes[1, 0]
        sigma = data['Sigma (dlnV/dlnr-1)'][::-1]  # Reverse to match radius
        ax.step(radius, sigma, 'm-', linewidth=2, where='post')
        ax.set_xlabel('Radius (cm)')
        ax.set_ylabel('Sigma')
        ax.set_title('Sigma Profile')
        ax.set_xscale('log')
        ax.axhline(y=-1, color='r', linestyle='--', label='Sigma = -1')
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    # Plot 5: Electron density
    if 'Electron density (1/cm^3)' in data and radius is not None:
        ax = axes[1, 1]
        ne = data['Electron density (1/cm^3)'][::-1]  # Reverse to match radius
        ax.step(radius, ne, 'c-', linewidth=2, where='post')
        ax.set_xlabel('Radius (cm)')
        ax.set_ylabel('Electron density (cm$^{-3}$)')
        ax.set_title('Electron Density Profile')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.grid(True, alpha=0.3)
    
    # Plot 6: Mass fractions (linear scale)
    if radius is not None:
        ax = axes[1, 2]
        
        # Find mass fraction sections
        for section in data:
            if 'mass fraction' in section.lower():
                elem = section.split()[0]
                ax.step(radius, data[section][::-1], label=elem, linewidth=2, where='post')  # Reverse to match radius
    
    ax.set_xlabel('Radius (cm)')
    ax.set_ylabel('Mass Fraction')
    ax.set_title('Mass Fractions')
    ax.set_xscale('log')
    ax.set_ylim(-0.05, 1.05)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Plot 7: Acceleration
    if 'Velocity (km/s)' in data and radius is not None:
        ax = axes[2, 0]
        velocity = data['Velocity (km/s)'][::-1] * 1e5  # Reverse and convert km/s to cm/s
        
        # Calculate acceleration: a = dv/dt = (dv/dr) * (dr/dt) = (dv/dr) * v
        # Use central difference for derivative
        accel = np.zeros_like(radius)
        for i in range(1, len(radius)-1):
            dv = velocity[i+1] - velocity[i-1]
            dr = radius[i+1] - radius[i-1]
            accel[i] = (dv / dr) * velocity[i]
        
        # Boundary values
        accel[0] = accel[1]
        accel[-1] = accel[-2]
        
        # Convert to km/s²
        accel_km = accel / 1e5
        
        ax.step(radius, accel_km, 'k-', linewidth=2, where='post')
        ax.set_xlabel('Radius (cm)')
        ax.set_ylabel('Acceleration (km/s²)')
        ax.set_title('Acceleration Profile')
        ax.set_xscale('log')
        ax.grid(True, alpha=0.3)
    
    # Plot 8: Kappa
    if 'Kappa (cm^2/gm)' in data and radius is not None:
        ax = axes[2, 1]
        kappa = data['Kappa (cm^2/gm)'][::-1]  # Reverse to match radius
        ax.step(radius, kappa, 'orange', linewidth=2, where='post')
        ax.set_xlabel('Radius (cm)')
        ax.set_ylabel('Kappa (cm²/g)')
        ax.set_title('Opacity Profile')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.grid(True, alpha=0.3)
    
    # Plot 9: Atom density
    if 'Atom density (1/cm^3)' in data and radius is not None:
        ax = axes[2, 2]
        na = data['Atom density (1/cm^3)'][::-1]  # Reverse to match radius
        ax.step(radius, na, 'purple', linewidth=2, where='post')
        ax.set_xlabel('Radius (cm)')
        ax.set_ylabel('Atom Density (cm⁻³)')
        ax.set_title('Atom Density Profile')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plot_file = filename.replace('.dat', '_final_plot.png').replace('_DATA', '_DATA_plot')
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    print(f"Plot saved to {plot_file}")
    plt.close()
    
    # Generate separate dr vs r plot for AMR grid visualization
    # This plot compares original HERACLES grid vs new AMR grid
    # Remove this from here - will be generated in main conversion function where both grids are available

if __name__ == '__main__':
    # Use command-line arguments if provided, otherwise use defaults
    # Default paths are relative to SN2024pba base directory
    if len(sys.argv) > 1:
        heracles_dir = sys.argv[1]
    else:
        heracles_dir = '../heracles/hydro'
    
    if len(sys.argv) > 2:
        try:
            target_days = float(sys.argv[2])
        except ValueError:
            target_days = 15
    else:
        target_days = 15
    
    if len(sys.argv) > 3:
        output_file = sys.argv[3]
    else:
        output_file = '../cmfgen/cmfgen_supernova/SN_HYDRO_DATA'

    # Parse optional arguments
    tau_min = None
    tau_max = None
    n_grid = 100
    adaptive_grid = False  # Disabled for NRAY_MAX stability
    gradient_threshold = 0.5
    initial_profile = None
    smooth_sigma = False
    sigma_window = 9
    sigma_polyorder = 3
    sigma_threshold = 2.0
    velocity_max_kms = None
    
    # Parse optional flags first (they can be anywhere in the argument list)
    positional_args = []
    for arg in sys.argv[1:]:
        if arg.startswith('--'):
            if arg == '--adaptive-grid':
                adaptive_grid = True
            elif arg.startswith('--gradient-threshold='):
                gradient_threshold = float(arg.split('=')[1])
            elif arg.startswith('--initial-profile='):
                initial_profile = arg.split('=')[1]
            elif arg == '--smooth-sigma':
                smooth_sigma = True
            elif arg.startswith('--sigma-window='):
                sigma_window = int(arg.split('=')[1])
            elif arg.startswith('--sigma-polyorder='):
                sigma_polyorder = int(arg.split('=')[1])
            elif arg.startswith('--sigma-threshold='):
                sigma_threshold = float(arg.split('=')[1])
            elif arg.startswith('--velocity-max-kms='):
                velocity_max_kms = float(arg.split('=')[1])
        else:
            positional_args.append(arg)
    
    # Parse positional arguments
    if len(positional_args) >= 1:
        heracles_dir = positional_args[0]
    if len(positional_args) >= 2:
        target_days = float(positional_args[1])
    if len(positional_args) >= 3:
        output_file = positional_args[2]
    if len(positional_args) >= 4:
        tau_min_str = positional_args[3]
        tau_min = float(tau_min_str) if tau_min_str.lower() != 'none' else None
    if len(positional_args) >= 5:
        tau_max_str = positional_args[4]
        tau_max = float(tau_max_str) if tau_max_str.lower() != 'none' else None
    if len(positional_args) >= 6:
        n_grid = int(positional_args[5])
    
    convert_heracles_to_cmfgen(heracles_dir, target_days=target_days, output_file=output_file,
                               tau_min=tau_min, tau_max=tau_max, n_grid=n_grid,
                               adaptive_grid=adaptive_grid, gradient_threshold=gradient_threshold,
                               initial_profile=initial_profile,
                               smooth_sigma=smooth_sigma, sigma_window=sigma_window,
                               sigma_polyorder=sigma_polyorder, sigma_threshold=sigma_threshold,
                               velocity_max_kms=velocity_max_kms)
    
    # Read and plot the created file
    print("\n" + "="*60)
    print("Reading and plotting the created CMFGEN file...")
    print("="*60)
    
    try:
        plot_cmfgen_file(output_file)
    except Exception as e:
        print(f"ERROR plotting file: {e}")
