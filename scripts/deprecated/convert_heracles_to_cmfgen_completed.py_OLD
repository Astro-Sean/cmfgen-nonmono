#!/usr/bin/env python3
"""
[DEPRECATED] Convert HERACLES snapshot to SN_HYDRO_DATA format for CMFGEN

DEPRECATED (2026-07-09): HERACLES does not fit the science case for this project
(interacting Type IIn supernovae). The project now uses the STELLA -> CMFGEN
pipeline exclusively. This file is retained for reference only.

For the active conversion path, see convert_stella_swd_to_cmfgen.py.

---

Convert HERACLES snapshot to SN_HYDRO_DATA format for CMFGEN
Uses adaptive grid refinement method to replace AMR
"""

import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.interpolate import CubicSpline

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
      e_rad = a*T⁴/rho
    """
    KB = 1.38062e-16    # erg/K
    UMA = 1.660531e-24  # g (atomic mass unit)
    a_rad = 7.5657e-15  # erg cm⁻³ K⁻⁴
    GAMMA = 5.0/3.0
    m_H = UMA

    # Kinetic energy (erg/cm³)
    e_kin = 0.5 * density * velocity_cms**2
    # Internal energy (erg/cm³) - ensure positive
    e_int = np.maximum(e_tot - e_kin, 1e-30)

    # Initial guess: ideal gas only (no radiation)
    T_guess = e_int * (GAMMA-1) * mu_loc * UMA / (density * KB)
    T_guess = np.maximum(T_guess, 2000.0)  # Avoid too-low initial guess

    # Vectorized Newton-Raphson (20 iterations max)
    T = T_guess.copy()
    for _ in range(20):
        # Radiation energy density
        e_rad = a_rad * T**4 / density
        # Gas energy density (from total internal energy)
        e_gas = e_int - e_rad
        # Gas pressure (from energy)
        P_gas = (GAMMA-1) * e_gas
        # Ideal gas law pressure
        P_gas_ideal = density * KB * T / (mu_loc * UMA)

        # Residual: P_gas (from energy) - P_gas (from ideal gas)
        residual = P_gas - P_gas_ideal

        # Derivative of residual w.r.t. T
        de_rad_dT = 4 * a_rad * T**3 / density
        de_gas_dT = -de_rad_dT
        dP_gas_dT = (GAMMA-1) * de_gas_dT
        dP_gas_ideal_dT = density * KB / (mu_loc * UMA)
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

    # Clip to physical range (2000 K to 1e8 K)
    T = np.clip(T, 2000.0, 1e8)
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
    
    # Apply physical constraints
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
    return calculate_sigma_with_smoothing(r, v, min_dlnr=min_dlnr, sigma_min=-0.99, sigma_max=10.0, smooth_window=3)
MAX_DRHO_PCT = 15.0  # Max 15% jump in density (lowered for better shock resolution)
MAX_DT_PCT = 20.0  # Max 20% jump in temperature
MAX_DV_PCT = 8.0  # Max 8% jump in velocity (lowered for better shock resolution)
MAX_GRID_POINTS = 400  # Maximum number of grid points for CMFGEN (safe limit for NRAY_MAX)
N_BASE = 100  # Number of base points for initial grid
N_ITERATIONS = 15  # Number of refinement iterations (increased for better shock resolution)
MIN_DRHO_PCT = 2.0  # Min 2% jump in density for coarsening
MIN_DT_PCT = 1.0  # Min 1% jump in temperature for coarsening
MIN_DV_PCT = 1.0  # Min 1% jump in velocity for coarsening


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
    # Avoid division by zero
    T_safe = T
    
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


def compute_tau_from_outer(radius_cm, rho_cgs, kappa_cgs=None):
    """
    Compute optical depth from outer boundary inward.
    
    Parameters:
    -----------
    radius_cm : array
        Radius in cm (inner->outer)
    rho_cgs : array
        Density in g/cm³ (inner->outer)
    kappa_cgs : array, optional
        Opacity in cm²/g (inner->outer). If None, use electron scattering fallback.
    
    Returns:
    --------
    tau : array
        Optical depth (inner->outer, with tau[-1] = 0 at outer boundary)
    """
    n = len(radius_cm)
    tau = np.zeros(n)
    
    # Use electron scattering opacity if kappa not provided
    if kappa_cgs is None or len(kappa_cgs) != n:
        # Approximate electron scattering: kappa_es = 0.2 * (1 + X) cm²/g
        # For H-rich ejecta, kappa_es ≈ 0.34 cm²/g
        kappa = np.ones(n) * 0.34
    else:
        kappa = np.asarray(kappa_cgs, float)
    
    # Integrate from outer in (tau[-1] = 0)
    for j in range(n - 2, -1, -1):
        dr_cm = abs(radius_cm[j + 1] - radius_cm[j])
        if dr_cm <= 0.0:
            continue
        chimid = kappa[j] * rho_cgs[j]
        tau[j] = tau[j + 1] + max(chimid * dr_cm, 0.0)
    
    return tau


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


def refine_grid(radius_cm, vel_cms, rho_cgs, temp_104k, ed_cgs=None, output_dir=None):
    """
    Optimized adaptive grid refinement: minimal points with strict gradient enforcement.
    Algorithm:
    1. Start with 50 logarithmic base points (reduced from 100).
    2. For up to 10 iterations, add points at geometric mean where gradients exceed thresholds.
    3. Coarsen grid by removing points in low-gradient regions (if removal doesn't violate thresholds).
    4. Add 10 outer boundary points (reduced from 20) for CMFGEN stability.

    Parameters:
    -----------
    radius_cm : array (outer->inner order, CMFGEN format)
    vel_cms : array (outer->inner)
    rho_cgs : array (outer->inner)
    temp_104k : array (outer->inner)
    ed_cgs : array (optional, outer->inner)
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

    print(f"[Grid Refinement] Starting with {n_base} base points")

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

            # Calculate percentage jumps
            drho_pct = abs(rho2 - rho1) / rho1 * 100 if rho1 > 1e-30 else 0
            dt_pct = abs(t2 - t1) / t1 * 100 if t1 > 1e-6 else 0
            dv_pct = abs(v2 - v1) / v1 * 100 if v1 > 1e-10 else 0

            # Check for violations
            if (drho_pct > MAX_DRHO_PCT or
                dt_pct > MAX_DT_PCT or
                dv_pct > MAX_DV_PCT):
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

            new_points.append((i+1, r_mid, v_mid, rho_mid, t_mid, ne_mid))

        # Insert points in reverse order to maintain indices
        for idx, r, v, rho, t, ne in sorted(new_points, key=lambda x: -x[0]):
            r_refined = np.insert(r_refined, idx, r)
            v_refined = np.insert(v_refined, idx, v)
            rho_refined = np.insert(rho_refined, idx, rho)
            t_refined = np.insert(t_refined, idx, t)
            ne_refined = np.insert(ne_refined, idx, ne)

        print(f"[Grid Refinement] Iteration {iteration+1}: {len(r_refined)} points")

    # --- Step 3: Coarsen grid by removing unnecessary points ---
    # Only remove points if:
    # 1. All adjacent gradients are below MIN_*_PCT (low-gradient region).
    # 2. Removing the point doesn't create a gradient violation.
    points_to_remove = []
    for i in range(1, len(r_refined) - 1):
        # Get neighbors
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

        # Check if this is a low-gradient region
        if (max(drho_prev, drho_next) < MIN_DRHO_PCT and
            max(dt_prev, dt_next) < MIN_DT_PCT and
            max(dv_prev, dv_next) < MIN_DV_PCT):
            # Check if removing this point would violate thresholds
            drho_combined = abs(rho_next - rho_prev) / rho_prev * 100 if rho_prev > 1e-30 else 0
            dt_combined = abs(t_next - t_prev) / t_prev * 100 if t_prev > 1e-6 else 0
            dv_combined = abs(v_next - v_prev) / v_prev * 100 if v_prev > 1e-10 else 0

            if (drho_combined <= MAX_DRHO_PCT and
                dt_combined <= MAX_DT_PCT and
                dv_combined <= MAX_DV_PCT):
                points_to_remove.append(i)

    if points_to_remove:
        print(f"[Grid Refinement] Removing {len(points_to_remove)} unnecessary points")
        r_refined = np.delete(r_refined, points_to_remove)
        v_refined = np.delete(v_refined, points_to_remove)
        rho_refined = np.delete(rho_refined, points_to_remove)
        t_refined = np.delete(t_refined, points_to_remove)
        ne_refined = np.delete(ne_refined, points_to_remove)

    # --- Step 4: Enforce minimum dr to prevent numerical instability (DISABLED) ---
    # Disabled - too aggressive for this data, causes extreme gradient violations
    # Robust sigma calculation (with min_dlnr=1e-5) handles small dr gracefully
    # MIN_DR_FRACTION = 0.001  # 0.1% of total radius range
    # min_dr = MIN_DR_FRACTION * (r_refined[-1] - r_refined[0])
    # print(f"[Grid Refinement] Enforcing min_dr = {min_dr:.2e} cm")
    # i = 1
    # points_removed = 0
    # while i < len(r_refined) - 1:
    #     if r_refined[i+1] - r_refined[i] < min_dr:
    #         r_refined = np.delete(r_refined, i+1)
    #         v_refined = np.delete(v_refined, i+1)
    #         rho_refined = np.delete(rho_refined, i+1)
    #         t_refined = np.delete(t_refined, i+1)
    #         ne_refined = np.delete(ne_refined, i+1)
    #         points_removed += 1
    #     else:
    #         i += 1
    # if points_removed > 0:
    #     print(f"[Grid Refinement] Removed {points_removed} points to enforce min_dr")
    # print(f"[Grid Refinement] Final points after min_dr enforcement: {len(r_refined)}")

    # --- Step 5: Add outer boundary points (DISABLED to avoid velocity jumps) ---
    n_outer_points = 0  # DISABLED - uses only HERACLES data, no artificial points
    if len(r_refined) < MAX_GRID_POINTS - n_outer_points:
        r_outer = r_refined[-1]
        r_inner_outer = r_refined[-2] if len(r_refined) > 1 else r_refined[-1] * 0.9
        rho_norm = rho_refined[-1] * r_outer**2
        T_outer = t_refined[-1]

        outer_points = []
        for i in range(1, n_outer_points + 1):
            spacing_fraction = 2.0 / 100.0  # 2% spacing
            r_new = r_outer * (1 - i * spacing_fraction)
            if r_new > r_inner_outer and r_new < r_outer:
                # Interpolate from ORIGINAL data (not refined grid) to avoid velocity jumps
                log_r_new = np.log10(r_new)
                v_new = 10 ** np.interp(log_r_new, log_r_oi, np.log10(v_oi))
                rho_new = rho_norm / r_new**2  # r^-2 density profile
                t_new = T_outer  # Constant temperature
                outer_points.append((r_new, v_new, rho_new, t_new))

        if outer_points:
            outer_points.sort(key=lambda x: x[0])  # Sort by radius
            for r, v, rho, t in outer_points:
                if len(r_refined) >= MAX_GRID_POINTS:
                    break
                r_refined = np.insert(r_refined, -1, r)
                v_refined = np.insert(v_refined, -1, v)
                rho_refined = np.insert(rho_refined, -1, rho)
                t_refined = np.insert(t_refined, -1, t)
                ne_refined = np.insert(ne_refined, -1, 0.0)  # Placeholder

            print(f"[Grid Refinement] Added {len(outer_points)} outer boundary points")

    print(f"[Grid Refinement] Final grid: {len(r_refined)} points")

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
    
    for i in range(n_refined - 1):
        drho_pct[i] = (rho_refined[i+1] - rho_refined[i]) / rho_refined[i] * 100 if rho_refined[i] > 1e-30 else 0
        dt_pct[i] = (t_refined[i+1] - t_refined[i]) / t_refined[i] * 100 if t_refined[i] > 1e-6 else 0
        dv_pct[i] = (v_refined[i+1] - v_refined[i]) / v_refined[i] * 100 if v_refined[i] > 1e-10 else 0
    
    # Create figure with subplots (2 rows x 3 columns)
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))

    
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
        raise ValueError(f"No zones found within tau range [{tau_min}, {tau_max}]. "
                        f"tau range in data: [{tau.min():.3e}, {tau.max():.3e}]")
    
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
    x_old = np.log(np.asarray(radius_old)[::-1])
    y_old = np.asarray(values_old)[::-1]
    x_new = np.log(np.asarray(radius_new)[::-1])

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

    return n

def convert_heracles_to_cmfgen(heracles_dir, target_days=30, output_file='SN_HYDRO_DATA', 
                               tau_min=None, tau_max=None, n_grid=100, adaptive_grid=True, 
                               gradient_threshold=0.5, initial_profile=None):
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
            # New profile format: r_cm rho_cgs p_cgs u_cgs T_K fx_H fx_He fx_C fx_O fx_Si fx_Fe
            initial_r = profile_data[:, 0]  # cm
            initial_T = profile_data[:, 4]  # K
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
        
        # Derive temperature from total energy with radiation pressure
        KB    = 1.38062e-16   # erg/K
        UMA   = 1.660531e-24  # g (atomic mass unit)
        GAMMA = 5.0 / 3.0
        A_VALS = [1.008, 4.003, 12.011, 15.999, 28.085, 55.845]
        Z_VALS = [1.0,   2.0,   6.0,    8.0,   14.0,   26.0  ]
        
        # Only calculate temperature from HERACLES energy if not using initial profile temperature
        if not use_initial_temperature:
            E_tot = np.asarray(d.E).ravel()
            
            # Mean molecular weight from mass fractions
            fx_arr = np.asarray(getattr(d, "fx", np.zeros((len(radius_cm), 6)))).reshape(len(radius_cm), -1)
            nfx    = fx_arr.shape[1]
            mf_raw = [fx_arr[:, i] if i < nfx else np.zeros(len(radius_cm)) for i in range(6)]
            mf_sum = np.maximum(np.sum(mf_raw, axis=0), 1e-300)
            mf_n   = [mf_raw[i] / mf_sum for i in range(6)]
            inv_mu = np.zeros(len(radius_cm))
            for i in range(6):
                inv_mu += mf_n[i] * (1.0 + Z_VALS[i]) / A_VALS[i]
            mu_loc = 1.0 / np.maximum(inv_mu, 1e-30)
            
            # Calculate temperature with radiation pressure
            temperature = calculate_temperature_with_radiation(density, velocity_cms, E_tot, mu_loc)
            temperature = np.where(np.isfinite(temperature) & (temperature > 0), temperature, 1000.0)
            print(f"Temperature derived from energy field (gas + radiation pressure)")
            print(f"  Mean molecular weight mu: min={mu_loc.min():.3e}, max={mu_loc.max():.3e}")
            print(f"  Density: min={density.min():.3e}, max={density.max():.3e} g/cm³")
            print(f"  Temperature (from HERACLES): min={temperature.min():.3e}, max={temperature.max():.3e} K")
        else:
            # Interpolate initial profile temperature onto HERACLES grid
            print(f"Interpolating initial profile temperature onto HERACLES grid...")
            from scipy.interpolate import interp1d
            # Log-log interpolation for better accuracy over large dynamic range
            log_r_initial = np.log10(initial_r)
            log_T_initial = np.log10(initial_T)
            log_r_heracles = np.log10(radius_cm)
            
            # Create interpolation function
            interp_func = interp1d(log_r_initial, log_T_initial, kind='linear', 
                                   bounds_error=False, fill_value='extrapolate')
            log_T_interp = interp_func(log_r_heracles)
            temperature = 10**log_T_interp
            
            print(f"  Interpolated temperature range: {temperature.min():.2e} to {temperature.max():.2e} K")
            print(f"  HERACLES radius range: {radius_cm.min():.2e} to {radius_cm.max():.2e} cm")
            print(f"  Initial profile radius range: {initial_r.min():.2e} to {initial_r.max():.2e} cm")
            
            # Check if HERACLES grid includes the innermost high-temperature region
            if radius_cm.min() > initial_r.min():
                print(f"  WARNING: HERACLES grid does not include innermost region")
                print(f"  HERACLES r_min = {radius_cm.min():.2e} cm, Initial profile r_min = {initial_r.min():.2e} cm")
                print(f"  Highest temperature in HERACLES grid will be at r = {radius_cm.min():.2e} cm")
                print(f"  This corresponds to initial profile T ≈ {temperature.max():.2e} K")
                print(f"  The intended T_inner_max = {initial_T.max():.2e} K is at r = {initial_r[np.argmax(initial_T)]:.2e} cm")
            
            print(f"Using initial profile temperature - skipping HERACLES energy calculation")
            
            # Still need to calculate mu_loc for later use (electron density calculation)
            fx_arr = np.asarray(getattr(d, "fx", np.zeros((len(radius_cm), 6)))).reshape(len(radius_cm), -1)
            nfx    = fx_arr.shape[1]
            mf_raw = [fx_arr[:, i] if i < nfx else np.zeros(len(radius_cm)) for i in range(6)]
            mf_sum = np.maximum(np.sum(mf_raw, axis=0), 1e-300)
            mf_n   = [mf_raw[i] / mf_sum for i in range(6)]
            inv_mu = np.zeros(len(radius_cm))
            for i in range(6):
                inv_mu += mf_n[i] * (1.0 + Z_VALS[i]) / A_VALS[i]
            mu_loc = 1.0 / np.maximum(inv_mu, 1e-30)
        
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
            
            # Define sigmoid transition parameters
            H_transition = 0.5  # Transition point at H = 0.5
            transition_width = 0.1  # Width of transition region
            T_target = 5000.0  # Target temperature for metal-rich regions (K)
            
            # Apply sigmoid function: f(H) = 1 / (1 + exp(-(H - H_center) / width))
            # This gives weight to original temperature (near 1 for H-rich) vs target (near 0 for metal-rich)
            sigmoid_weight = 1.0 / (1.0 + np.exp(-(H_fraction - H_transition) / transition_width))
            temperature = sigmoid_weight * T_original + (1.0 - sigmoid_weight) * T_target
            
            print(f"  Temperature (sigmoid smoothing applied)")
            print(f"    Transition: H={H_transition} with width={transition_width}")
            print(f"    Target temperature: {T_target} K")
            print(f"    Temperature range: {temperature.min():.2e} to {temperature.max():.2e} K")
        else:
            print(f"  Temperature (sigmoid smoothing DISABLED - using initial profile temperature)")
            print(f"    Temperature range: {temperature.min():.2e} to {temperature.max():.2e} K")
        
        # Get mass fractions
        fx = np.asarray(getattr(d, "fx", np.zeros((len(radius_cm), 6))))
        if len(fx.shape) == 4:
            fx = fx.reshape(-1, 6)
        
        print(f"fx shape: {fx.shape}")
        print(f"fx min/max: {fx.min():.3e} / {fx.max():.3e}")
        print(f"fx first 5 rows:\n{fx[:5]}")
        print(f"fx last 5 rows:\n{fx[-5:]}")
        
        # Build element arrays (H, He, C, O, Si, Fe) from HERACLES
        nfx = fx.shape[1] if len(fx.shape) > 1 else 0
        mass_fractions = {}
        elem_names = ['H', 'He', 'C', 'O', 'Si', 'Fe']
        for i, elem in enumerate(elem_names):
            if i < nfx:
                mass_fractions[elem.upper()] = fx[:, i]
                print(f"{elem.upper()}: min={mass_fractions[elem.upper()].min():.3e}, max={mass_fractions[elem.upper()].max():.3e}, std={mass_fractions[elem.upper()].std():.3e}")
            else:
                mass_fractions[elem.upper()] = np.zeros(len(radius_cm))
        
        # Add additional elements with zeros if needed
        if 'N' not in mass_fractions:
            mass_fractions['N'] = np.zeros(len(radius_cm))
        if 'NE' not in mass_fractions:
            mass_fractions['NE'] = np.zeros(len(radius_cm))
        if 'MG' not in mass_fractions:
            mass_fractions['MG'] = np.zeros(len(radius_cm))
        if 'S' not in mass_fractions:
            mass_fractions['S'] = np.zeros(len(radius_cm))
        
        print(f"\nBEFORE REGRIDDING (n_points={len(radius_cm)}):")
        for elem in ['H', 'HE', 'C', 'O', 'SI', 'FE']:
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
        print("Attempting to read profile files directly...")
        profile_files = [f for f in os.listdir(snap_dir) if 'profile' in f.lower()]
        if profile_files:
            profile_file = os.path.join(snap_dir, profile_files[0])
            print(f"Reading profile from {profile_file}")
            profile_data = np.loadtxt(profile_file)
            
            # New profile format: r_cm rho_cgs p_cgs u_cgs T_K fx_H fx_He fx_C fx_O fx_Si fx_Fe
            radius_cm = profile_data[:, 0]  # cm
            density = profile_data[:, 1]  # g/cm³
            pressure = profile_data[:, 2]  # g/cm/s²
            velocity_cms = profile_data[:, 3]  # cm/s
            temperature = profile_data[:, 4]  # K
            velocity_kms = velocity_cms / 1e5  # Convert to km/s
            
            # Mass fractions
            n_points = len(radius_cm)
            mass_fractions = {
                'H': profile_data[:, 5],
                'HE': profile_data[:, 6],
                'C': profile_data[:, 7],
                'O': profile_data[:, 8],
                'SI': profile_data[:, 9],
                'FE': profile_data[:, 10]
            }
            
            print(f"Using initial profile temperature (not HERACLES energy)")
            print(f"  Temperature range: {temperature.min():.2e} to {temperature.max():.2e} K")
            use_initial_temperature = True
        else:
            print(f"ERROR: No profile file found in {snap_dir}")
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
    # Force fixed tau range for CMFGEN stability
    tau_min = 1e-4  # Force minimum tau
    tau_max = 30.0  # Force maximum tau
    print(f"Using fixed tau range: [{tau_min}, {tau_max}]")
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
    
    # STEP 5: Recompute trimmed tau, velocity_kms, sigma
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
        output_dir=output_dir
    )
    
    # Convert back to original units
    temperature_new = temperature_new * 1e4  # Convert 10^4 K back to K
    velocity_kms_new = velocity_cms_new / 1e5
    
    print(f"  New radius range: {r_new[0]:.2e} to {r_new[-1]:.2e} cm")
    
    # STEP 7: Adaptive grid refinement already interpolated all quantities
    # Just update the main arrays with refined values
    # refine_grid returns outer->inner order (CMFGEN format), so no reversal needed
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
    for elem in mass_fractions:
        mass_fractions[elem] = interp_outer_to_inner(radius_cm, mass_fractions[elem], r_new, method="linear")
        mass_fractions[elem] = np.maximum(mass_fractions[elem], 0.0)
    
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
    
    # Apply density floor to prevent CMFGEN zero RGU error
    # CMFGEN fails when density is too low (E-89 or lower)
    DENSITY_MIN = 1e-22  # Minimum density floor in g/cm³
    density_before_floor = density.copy()
    density = np.maximum(density, DENSITY_MIN)
    
    # Apply distance-dependent temperature floor to maintain reasonable DTDR
    # Without a floor, outer regions have very flat temperature gradients (DTDR ~ 1e-23)
    # causing CMFGEN to fail. Floor scales with radius to preserve inner gradient.
    temperature_before_floor = temperature.copy()
    
    # Distance-dependent floor: higher at outer boundary, decreases inward
    # T_floor(r) = T_outer * (r / r_max)^beta
    # where beta controls how quickly floor decreases inward
    r_norm = radius_cm / radius_cm.max()  # Normalized radius (0 to 1)
    T_outer_floor = 1000.0  # Floor temperature at outer boundary (K)
    beta = 0.5  # Floor decreases with radius (higher beta = steeper decrease)
    T_floor = T_outer_floor * (r_norm ** beta)
    
    # Apply floor
    temperature = np.maximum(temperature, T_floor)
    
    n_floored = np.sum(temperature_before_floor < T_floor)
    if n_floored > 0:
        print(f"  Applied distance-dependent temperature floor: {n_floored} zones raised")
        print(f"    Temperature before floor: min={temperature_before_floor.min():.2e}, max={temperature_before_floor.max():.2e} K")
        print(f"    Temperature after floor: min={temperature.min():.2e}, max={temperature.max():.2e} K")
        print(f"    Floor range: {T_floor.min():.2e} to {T_floor.max():.2e} K")

    # STEP 9: Recompute final tau and final sigma from the regridded arrays
    # Recalculate tau AFTER applying density floor for consistency
    # ... (rest of the code remains the same)
    tau = calculate_tau_from_outer(radius_cm, density, kappa_cm2_per_g_array)
    
    print(f"  Final optical depth range: {tau.min():.2e} to {tau.max():.2e}")
    
    # STEP 10: Validate lengths
    n_points = len(radius_cm)
    print(f"\nValidating final model (n_points={n_points})...")
    validate_final_model(radius_cm, density, velocity_kms, temperature, tau, sigma, 
                       kappa_cm2_per_g_array, mass_fractions)
    print("  All arrays have consistent length")
    
    # STEP 11: Write output (using only final regridded arrays)
    # Set up output parameters
    n_mass_fractions = 8  # HYD, HE, CARB, OXY, SIL, IRON, COB, NICK
    n_isotopes = 3
    time_days = actual_days
    energy_lum = 6.071478e6
    temperature_10k = temperature / 1e4
    
    # Calculate atom density and electron density on regridded arrays
    # Using mean atomic mass from mass fractions (fully ionised approximation)
    atomic_masses = {'H': 1.0, 'HE': 4.0, 'C': 12.0, 'N': 14.0, 'O': 16.0,
                     'NE': 20.0, 'MG': 24.0, 'SI': 28.0, 'S': 32.0, 'FE': 56.0}
    atomic_Z      = {'H': 1,   'HE': 2,   'C': 6,    'N': 7,    'O': 8,
                     'NE': 10,  'MG': 12,  'SI': 14,  'S': 16,   'FE': 26}
    
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
    
    # Verify arrays are in inner->outer order (CMFGEN default)
    # CMFGEN expects inner->outer order (increasing radius)
    # Arrays are in outer->inner order from interpolation step, so reverse them
    if radius[0] > radius[-1]:
        print("WARNING: radius array is in outer->inner order, reversing to inner->outer for CMFGEN")
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
    assert radius[0] < radius[-1], "ERROR: radius array should be in inner->outer order (increasing)"
    print(f"Verified: radius[0]={radius[0]:.2e} < radius[-1]={radius[-1]:.2e} (inner->outer order)")
    
    with open(output_file, 'w') as f:
        # Header (exact format matching working CMFGEN SN_HYDRO_DATA)
        f.write(f"Number of data points:               {n_points:>3d}\n")
        f.write(f"Number of mass fractions:            {n_mass_fractions:>2d}\n")
        f.write(f"Number of isotopes:                  {n_isotopes:>3d}\n")
        f.write(f"Time(days) since explosion:          {time_days:10.5E}\n")
        f.write(f"Energy Lum [lsun]:                   {energy_lum:11.6E}\n")
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
        
        # Mass fractions (rename to match CMFGEN expected names)
        # Output 8 mass fractions to match working SN_HYDRO_DATA: HYD, HE, CARB, OXY, SIL, IRON, COB, NICK
        elem_name_map = {
            'H': 'HYD',
            'HE': 'HE',
            'C': 'CARB',
            'O': 'OXY',
            'SI': 'SIL',
            'FE': 'IRON'
        }
        # Normalize mass fractions to sum to 1 at each depth
        # Include COB and NICK with placeholder values in the sum
        cob_placeholder = np.full(n_points, 1.0e-8)
        nick_placeholder = np.full(n_points, 1.0e-8)
        
        # Create a dictionary with all mass fractions including placeholders
        # Map HERACLES element names to CMFGEN names
        all_mass_fractions = {}
        for heracles_elem, cmfgen_elem in elem_name_map.items():
            if heracles_elem in mass_fractions:
                all_mass_fractions[cmfgen_elem] = mass_fractions[heracles_elem]
        all_mass_fractions['COB'] = cob_placeholder
        all_mass_fractions['NICK'] = nick_placeholder
        
        # Calculate sum of real mass fractions (exclude placeholders COB, NICK)
        real_elements = ['HYD', 'HE', 'CARB', 'OXY', 'SIL', 'IRON']
        mf_sum = np.zeros(n_points)
        for elem in real_elements:
            if elem in all_mass_fractions:
                mf_sum += all_mass_fractions[elem]
        
        print(f"Mass fraction sum before normalization (real elements only): min={mf_sum.min():.6f}, max={mf_sum.max():.6f}")
        
        # Normalize real mass fractions to sum to exactly 1 at each depth
        # This prevents CMFGEN from having to renormalize internally
        # Avoid division by zero
        mf_sum_safe = np.where(mf_sum > 1e-10, mf_sum, 1.0)
        for elem in real_elements:
            if elem in all_mass_fractions:
                all_mass_fractions[elem] = all_mass_fractions[elem] / mf_sum_safe
        
        # Set placeholder elements to zero (or tiny floor if needed)
        all_mass_fractions['COB'] = np.zeros(n_points)
        all_mass_fractions['NICK'] = np.zeros(n_points)
        
        # Verify normalization
        mf_sum_after = np.zeros(n_points)
        for elem in real_elements:
            if elem in all_mass_fractions:
                mf_sum_after += all_mass_fractions[elem]
        print(f"Mass fraction sum after normalization (real elements only): min={mf_sum_after.min():.6f}, max={mf_sum_after.max():.6f}")
        
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
            write_array_8_per_line(f, normalized_mf['NICK'])
            f.write("COB 56 mass fraction\n")
            write_array_8_per_line(f, normalized_mf['COB'])
            f.write("IRON 56 mass fraction\n")
            write_array_8_per_line(f, normalized_mf['IRON'])
    
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
    # Convert mass fractions to CMFGEN naming for plotting
    cmfgen_mass_fractions = {}
    elem_map = {'H': 'HYD', 'HE': 'HE', 'C': 'CARB', 'O': 'OXY', 'SI': 'SIL', 'FE': 'IRON'}
    for heracles_elem, cmfgen_elem in elem_map.items():
        if heracles_elem in mass_fractions:
            cmfgen_mass_fractions[cmfgen_elem] = mass_fractions[heracles_elem]
    
    # Plot comparison using refined grid data in inner->outer order (natural order from adaptive grid refinement)
    # SN_HYDRO_DATA is written in outer->inner order for CMFGEN, but comparison plot shows the actual refined grid
    # Original data is in outer->inner order, so reverse it to inner->outer to match refined grid
    plot_comparison(original_radius_cm[::-1], original_density[::-1], original_temperature[::-1], 
                   original_velocity_kms[::-1], original_tau[::-1], original_sigma_unclipped[::-1],
                   radius_cm, density, temperature, velocity_kms, tau, sigma, sigma_unclipped,
                   output_file, time_days, orig_mass_fractions={k: v[::-1] for k, v in original_mass_fractions.items()}, new_mass_fractions=cmfgen_mass_fractions)

def plot_comparison(orig_r, orig_rho, orig_T, orig_v, orig_tau, orig_sigma,
                   new_r, new_rho, new_T, new_v, new_tau, new_sigma, new_sigma_unclipped,
                   output_file, time_days, orig_mass_fractions=None, new_mass_fractions=None):
    """Plot before/after comparison of HERACLES to CMFGEN conversion"""
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon
    
    # Solar radius in cm
    rsun_cm = 6.96e10
    
    # Convert radius to r/rsun for all subplots (actual values, not log)
    orig_r_rsun = orig_r / rsun_cm
    new_r_rsun = new_r / rsun_cm
    
    # Create a 2x3 grid for the main comparison plots
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle(f'HERACLES to CMFGEN Conversion: {time_days}-day snapshot', fontsize=16)
    
    # Plot 1: Density vs radius (both axes log scale)
    ax = axes[0, 0]
    ax.loglog(orig_r_rsun, orig_rho, 'b--', alpha=0.7, linewidth=2, label='Original HERACLES (input grid)')
    # Create filled polygon for CMFGEN data
    poly_x = np.concatenate([new_r_rsun, new_r_rsun[::-1]])
    poly_y = np.concatenate([new_rho, np.zeros_like(new_rho)])
    ax.add_patch(Polygon(np.column_stack([poly_x, poly_y]), closed=True, 
                       facecolor='red', alpha=0.3, label='CMFGEN (filled)'))
    ax.loglog(new_r_rsun, new_rho, 'r-', alpha=0.7, linewidth=2, label='CMFGEN (outline)')
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
    poly_x = np.concatenate([new_r_rsun, new_r_rsun[::-1]])
    poly_y = np.concatenate([new_T, np.zeros_like(new_T)])
    ax.add_patch(Polygon(np.column_stack([poly_x, poly_y]), closed=True, 
                       facecolor='red', alpha=0.3, label='CMFGEN (filled)'))
    ax.loglog(new_r_rsun, new_T, 'r-', alpha=0.7, linewidth=2, label='CMFGEN (outline)')
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
    poly_x = np.concatenate([new_r_rsun, new_r_rsun[::-1]])
    poly_y = np.concatenate([new_v, np.zeros_like(new_v)])
    ax.add_patch(Polygon(np.column_stack([poly_x, poly_y]), closed=True, 
                       facecolor='red', alpha=0.3, label='CMFGEN (filled)'))
    ax.loglog(new_r_rsun, new_v, 'r-', alpha=0.7, linewidth=2, label='CMFGEN (outline)')
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
    ax.step(new_r_rsun, new_sigma_unclipped, 'r-', alpha=0.7, linewidth=2, label='CMFGEN (unclipped)', where='post')
    ax.axhline(y=-0.99, color='k', linestyle=':', alpha=0.5, label='σ=-0.99 (clipping limit)')
    ax.axhline(y=-1.0, color='gray', linestyle='--', alpha=0.5, label='σ=-1.0 (CMFGEN limit)')
    ax.set_xscale('log')
    ax.set_xlabel('r/R☉', fontsize=12)
    ax.set_ylabel('Sigma (unclipped)', fontsize=12)
    ax.set_title('Sigma Profile (Unclipped) - Symlog', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_yscale('symlog', linthresh=1.0)
    ax.set_xlim(left=10.0)  # Minimum at 10 Rsun
    
    # Plot 5: Optical depth vs radius
    ax = axes[1, 1]
    ax.step(orig_r_rsun, orig_tau, 'b--', alpha=0.7, linewidth=2, label='Original HERACLES (input grid)', where='post')
    poly_x = np.concatenate([new_r_rsun, new_r_rsun[::-1]])
    poly_y = np.concatenate([new_tau, np.zeros_like(new_tau)])
    ax.add_patch(Polygon(np.column_stack([poly_x, poly_y]), closed=True, 
                       facecolor='red', alpha=0.3, label='CMFGEN (filled)'))
    ax.step(new_r_rsun, new_tau, 'r-', alpha=0.7, linewidth=2, label='CMFGEN (outline)', where='post')
    ax.axhline(y=2/3, color='orange', linestyle=':', alpha=0.8, linewidth=2, label='τ = 2/3 (photosphere)')
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
        # Element mapping: HERACLES names -> CMFGEN names
        elem_map = {'H': 'HYD', 'HE': 'HE', 'C': 'CARB', 'O': 'OXY', 'SI': 'SIL', 'FE': 'IRON'}
        colors = ['blue', 'green', 'red', 'orange', 'purple', 'brown']
        
        # Plot each element
        for i, (heracles_elem, cmfgen_elem) in enumerate(elem_map.items()):
            if heracles_elem in orig_mass_fractions and cmfgen_elem in new_mass_fractions:
                # Plot input grid (dashed)
                ax.step(orig_r_rsun, orig_mass_fractions[heracles_elem], 
                       linestyle='--', alpha=0.7, linewidth=2, 
                       color=colors[i], label=f'{heracles_elem} (input)', where='post')
                # Plot CMFGEN (solid)
                ax.step(new_r_rsun, new_mass_fractions[cmfgen_elem], 
                       linestyle='-', alpha=0.7, linewidth=2, 
                       color=colors[i], label=f'{cmfgen_elem} (CMFGEN)', where='post')
        
        ax.set_xscale('log')
        ax.set_xlabel('r/R☉', fontsize=12)
        ax.set_ylabel('Mass Fraction', fontsize=12)
        ax.set_title('Mass Fraction Profile', fontsize=14)
        ax.legend(fontsize=7, ncol=3)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(left=10.0)  # Minimum at 10 Rsun
    else:
        ax.text(0.5, 0.5, 'Mass fractions not available', 
                ha='center', va='center', transform=ax.transAxes, fontsize=14)
        ax.set_title('Mass Fraction Profile', fontsize=14)
    
    plt.tight_layout()
    
    # Save plot
    plot_file = output_file.replace('.dat', '_comparison.png').replace('SN_HYDRO_DATA', 'SN_HYDRO_DATA_comparison')
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    print(f"Plot saved to {plot_file}")
    plt.close()

def write_array_8_per_line(f, array, format_str='%.15E'):
    """Write array with 8 values per line in scientific notation"""
    n = len(array)
    for i in range(0, n, 8):
        chunk = array[i:i+8]
        line = '    '.join([format_str % val for val in chunk])
        f.write(line + '\n')
    f.write('\n')

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
    plot_file = filename.replace('.dat', '_plot.png').replace('_DATA', '_DATA_plot')
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
                               initial_profile=initial_profile)
    
    # Read and plot the created file
    print("\n" + "="*60)
    print("Reading and plotting the created CMFGEN file...")
    print("="*60)
    
    try:
        plot_cmfgen_file(output_file)
    except Exception as e:
        print(f"ERROR plotting file: {e}")
