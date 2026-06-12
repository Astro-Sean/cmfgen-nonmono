"""Convert a HERACLES-style initial profile to STELLA *.hyd + *.abn input files.

STELLA input format (see OpenStella/readme.stella):

*.hyd  -- hydrodynamic profile
  Line 1: timeStart(s)  Nzon  mass_cut(Msun)  Rcen(cm)  rhoCen(g/cc)
  Lines 2..Nzon+1: km  dMr(Msun)  r(cm)  rho(g/cc)  T(K)  u(cm/s)  Mr(Msun)  0

*.abn  -- abundance profile (mass fractions per zone)
  Lines 1..Nzon: km  0  0  0  H  He  C  N  O  Ne  Na  Mg  Al  Si  S  Ar  Ca  FePeak  Ni58  Ni56

HERACLES profile columns (from header comment):
  r_cm  rho_cgs  p_cgs  u_cgs  fx_H  fx_He  fx_C  fx_O  fx_Si  fx_Fe

Usage:
  python convert_profile_to_stella.py INPUT_PROFILE MODEL_NAME OUTPUT_DIR [OPTIONS]
    --nzon=N          Target number of Lagrangian zones (default: 200)
    --nzon-max=N      Hard maximum zones (default: 500)
    --t-start=T       Start time in seconds (default: 0)
    --mass-cut=M      Core mass in Msun to exclude (default: 0)
    --xni=X           56Ni mass fraction for opacity table (default: from abn)
    --rcen=R          Central (inner boundary) radius in cm (default: first r in profile)
    --rho-cen=R       Central density g/cc (default: innermost rho)
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np


def parse_heracles_profile(fname: str) -> dict:
    """Read HERACLES input_profile (r, rho, p, u, fH, fHe, fC, fO, fSi, fFe)."""
    import re
    data = []
    explosion_age_days = None
    with open(fname) as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith('#'):
                m = re.search(r'SN explosion age\s*=\s*([\d.eE+\-]+)', stripped)
                if m:
                    explosion_age_days = float(m.group(1))
                continue
            cols = [float(x) for x in stripped.split()]
            if len(cols) < 10:
                continue
            data.append(cols)
    arr = np.array(data)
    return {
        'r_cm':               arr[:, 0],
        'rho':                arr[:, 1],
        'p':                  arr[:, 2],
        'u_cms':              arr[:, 3],
        'fH':                 arr[:, 4],
        'fHe':                arr[:, 5],
        'fC':                 arr[:, 6],
        'fO':                 arr[:, 7],
        'fSi':                arr[:, 8],
        'fFe':                arr[:, 9],
        'explosion_age_days': explosion_age_days,
    }


def temperature_from_pressure(p_cgs, rho_cgs: np.ndarray,
                               fH=None, fHe=None, fC=None, fO=None,
                               fSi=None, fFe=None,
                               mu: float = 0.62) -> np.ndarray:
    """Estimate temperature from pressure.

    Uses the larger of radiation temperature (P_rad = a*T^4/3) and ideal-gas
    temperature, with a floor of 3000 K (STELLA needs reasonable initial T).

    If mass-fraction arrays are provided, mu is computed per-zone from the
    composition (neutral-gas approximation).  Otherwise the fallback mu is used.
    """
    p = np.asarray(p_cgs, dtype=float)
    rho = np.asarray(rho_cgs, dtype=float)
    k_B = 1.380649e-16    # erg/K
    m_H = 1.6735575e-24   # g
    a_rad = 7.566e-15     # erg/cm^3/K^4

    # Compute mean molecular weight from composition if available
    if fH is not None:
        fH  = np.asarray(fH,  dtype=float)
        fHe = np.asarray(fHe, dtype=float) if fHe is not None else np.zeros_like(fH)
        fC  = np.asarray(fC,  dtype=float) if fC  is not None else np.zeros_like(fH)
        fO  = np.asarray(fO,  dtype=float) if fO  is not None else np.zeros_like(fH)
        fSi = np.asarray(fSi, dtype=float) if fSi is not None else np.zeros_like(fH)
        fFe = np.asarray(fFe, dtype=float) if fFe is not None else np.zeros_like(fH)
        mu_arr = 1.0 / (fH/1.0 + fHe/4.0 + fC/12.0 + fO/16.0 + fSi/28.0 + fFe/56.0)
        # Guard against pathological all-zero compositions
        mu_arr = np.where(np.isfinite(mu_arr), mu_arr, mu)
    else:
        mu_arr = mu

    T_gas = p * mu_arr * m_H / (rho * k_B)
    T_rad = (3.0 * p / a_rad) ** 0.25
    T = np.maximum(T_gas, T_rad)
    T = np.clip(T, 3000.0, 1e8)
    return T


def detect_composition_interfaces(r, fH, fHe, fC, fO, fSi, fFe, threshold=0.1):
    """Detect composition interfaces where mass fractions change rapidly.
    
    Parameters:
    -----------
    r: array
        Radius array (cm)
    fH, fHe, fC, fO, fSi, fFe: arrays
        Mass fraction arrays
    threshold: float
        Minimum relative change to trigger refinement (default: 0.1 = 10%)
    
    Returns:
    --------
    interface_indices: array
        Indices where composition interfaces are detected
    """
    interface_indices = []
    
    for f, name in zip([fH, fHe, fC, fO, fSi, fFe], ['H', 'He', 'C', 'O', 'Si', 'Fe']):
        # Calculate relative change in mass fraction
        df = np.abs(np.diff(f))
        df_max = np.maximum(f[:-1], f[1:])
        df_rel = df / np.maximum(df_max, 1e-10)
        
        # Find where relative change exceeds threshold
        interfaces = np.where(df_rel > threshold)[0]
        interface_indices.extend(interfaces)
    
    return np.unique(interface_indices)


def detect_shocks(r, rho, p, u, pressure_threshold=0.5, velocity_threshold=0.3):
    """Detect shock fronts based on pressure and velocity jumps.
    
    Parameters:
    -----------
    r: array
        Radius array (cm)
    rho: array
        Density array (g/cc)
    p: array
        Pressure array (dyne/cm^2)
    u: array
        Velocity array (cm/s)
    pressure_threshold: float
        Minimum relative pressure jump for shock detection (default: 0.5 = 50%)
    velocity_threshold: float
        Minimum relative velocity jump for shock detection (default: 0.3 = 30%)
    
    Returns:
    --------
    shock_indices: array
        Indices where shocks are detected
    """
    shock_indices = []
    
    # Pressure jumps
    dp = np.abs(np.diff(p))
    dp_max = np.maximum(p[:-1], p[1:])
    dp_rel = dp / np.maximum(dp_max, 1e-10)
    pressure_shocks = np.where(dp_rel > pressure_threshold)[0]
    shock_indices.extend(pressure_shocks)
    
    # Velocity jumps
    du = np.abs(np.diff(u))
    du_max = np.maximum(np.abs(u[:-1]), np.abs(u[1:]))
    du_rel = du / np.maximum(du_max, 1e-10)
    velocity_shocks = np.where(du_rel > velocity_threshold)[0]
    shock_indices.extend(velocity_shocks)
    
    return np.unique(shock_indices)


def detect_photosphere(r, rho, T, kappa_guess=0.3):
    """Detect photosphere region where optical depth tau ≈ 1-2/3.
    
    Parameters:
    -----------
    r: array
        Radius array (cm)
    rho: array
        Density array (g/cc)
    T: array
        Temperature array (K)
    kappa_guess: float
        Guess for opacity (cm^2/g), default 0.3 for electron scattering
    
    Returns:
    --------
    photosphere_indices: array
        Indices near photosphere
    """
    # Calculate approximate optical depth from outside in
    # tau = integral(kappa * rho * dr) from infinity to r
    dr = np.diff(r)
    dtau = kappa_guess * rho[1:] * dr
    tau = np.cumsum(dtau[::-1])[::-1]  # Integrate from outside in
    
    # Find where tau is between 0.5 and 2.0 (photosphere region)
    photosphere_mask = (tau > 0.5) & (tau < 2.0)
    photosphere_indices = np.where(photosphere_mask)[0]
    
    return photosphere_indices


def detect_high_curvature(r, quantity, threshold=0.5):
    """Detect regions with high second derivative (high curvature).
    
    Parameters:
    -----------
    r: array
        Radius array (cm)
    quantity: array
        Physical quantity (e.g., density, temperature)
    threshold: float
        Threshold for normalized second derivative
    
    Returns:
    --------
    curvature_indices: array
        Indices where curvature is high
    """
    # Calculate second derivative using central differences
    dq_dr = np.gradient(quantity, r)
    d2q_dr2 = np.gradient(dq_dr, r)
    
    # Normalize by quantity to get relative curvature
    curvature = np.abs(d2q_dr2 * r**2 / np.maximum(quantity, 1e-10))
    
    # Find where curvature exceeds threshold
    curvature_indices = np.where(curvature > threshold)[0]
    
    return curvature_indices


def rezone(r, rho, p, u, fH, fHe, fC, fO, fSi, fFe, nzon_target: int,
           outer_fraction=0.3, outer_bias=2.0, transition_radius=1e13,
           check_gradients=True, max_dr_r=0.15, max_dT_T=0.20, max_drho_rho=0.20,
           max_du_u=0.20, adaptive_refinement=True, max_iterations=10, nzon_max=500,
           detect_composition=True, detect_shock_features=True, detect_photo_features=True,
           detect_curvature=True):
    """Resample profile to nzon_target zones using logarithmic radius spacing.

    Log-r spacing is essential for profiles spanning many decades in radius
    (e.g. RSG shell at 1e12 cm through ejecta at 3e16 cm = 4.5 dex).
    Pure mass-based spacing fails when regions differ greatly in mass: a
    2 Msun shell overwhelms a 0.009 Msun wind, leaving the wind with 0 zones.

    The new grid is spaced uniformly in log(r), which gives dr/r = constant —
    appropriate for a self-similar power-law structure.  All quantities are
    then interpolated from the original fine grid onto the new radii.

    If adaptive_refinement is enabled, the grid is iteratively refined in
    regions with large gradients to ensure smooth profiles, up to nzon_max.
    Additional refinement criteria can be enabled:
    - detect_composition: Add zones at composition interfaces
    - detect_shock_features: Add zones at shock fronts
    - detect_photo_features: Add zones near photosphere (tau ≈ 1)
    - detect_curvature: Add zones in high-curvature regions

    Parameters:
    -----------
    check_gradients: bool
        If True, check for large jumps and report problematic regions
    max_dr_r: float
        Maximum allowed relative jump in radius (default: 0.15 = 15%)
    max_dT_T: float
        Maximum allowed relative jump in temperature (default: 0.20 = 20%)
    max_drho_rho: float
        Maximum allowed relative jump in density (default: 0.20 = 20%)
    max_du_u: float
        Maximum allowed relative jump in velocity (default: 0.20 = 20%)
    adaptive_refinement: bool
        If True, iteratively add zones in high-gradient regions (default: True)
    max_iterations: int
        Maximum number of adaptive refinement iterations (default: 10)
    nzon_max: int
        Maximum number of zones allowed (hard limit, default: 500)
    detect_composition: bool
        If True, detect and refine composition interfaces (default: True)
    detect_shock_features: bool
        If True, detect and refine shock fronts (default: True)
    detect_photo_features: bool
        If True, detect and refine photosphere region (default: True)
    detect_curvature: bool
        If True, detect and refine high-curvature regions (default: True)
    """
    n = len(r)
    if n == nzon_target:
        return r, rho, p, u, fH, fHe, fC, fO, fSi, fFe

    # Start with uniform log-r spacing
    r_new = np.logspace(np.log10(r[0]), np.log10(r[-1]), nzon_target)

    def interp_logr(x, r_target):
        return np.interp(r_target, r, x)

    # Adaptive refinement loop
    if adaptive_refinement:
        print(f"\n[Adaptive Refinement] Starting with {nzon_target} zones (max: {nzon_max})...")
        
        for iteration in range(max_iterations):
            # Check if we've reached the maximum
            if len(r_new) >= nzon_max:
                print(f"  Reached maximum zone limit ({nzon_max}), stopping refinement")
                break

            # Interpolate quantities onto current grid
            rho_new = interp_logr(rho, r_new)
            p_new   = interp_logr(p, r_new)
            u_new   = interp_logr(u, r_new)
            fH_new  = np.clip(interp_logr(fH, r_new),  0, 1)
            fHe_new = np.clip(interp_logr(fHe, r_new), 0, 1)
            fC_new  = np.clip(interp_logr(fC, r_new),  0, 1)
            fO_new  = np.clip(interp_logr(fO, r_new),  0, 1)
            fSi_new = np.clip(interp_logr(fSi, r_new), 0, 1)
            fFe_new = np.clip(interp_logr(fFe, r_new), 0, 1)

            # Renormalise mass fractions
            mf_sum = fH_new + fHe_new + fC_new + fO_new + fSi_new + fFe_new
            mf_sum = np.where(mf_sum < 1e-10, 1.0, mf_sum)
            fH_new  /= mf_sum
            fHe_new /= mf_sum
            fC_new  /= mf_sum
            fO_new  /= mf_sum
            fSi_new /= mf_sum
            fFe_new /= mf_sum

            # Compute temperature
            T_new = temperature_from_pressure(p_new, rho_new, fH=fH_new, fHe=fHe_new,
                                              fC=fC_new, fO=fO_new, fSi=fSi_new, fFe=fFe_new)

            # Check for gradient violations
            dr_r = np.abs(np.diff(r_new) / r_new[:-1])
            dT_T = np.abs(np.diff(T_new) / np.maximum(T_new[:-1], 1e-10))
            drho_rho = np.abs(np.diff(rho_new) / np.maximum(rho_new[:-1], 1e-30))
            du_u = np.abs(np.diff(u_new) / np.maximum(np.abs(u_new[:-1]), 1e-10))

            # Find violations
            violations = (dr_r > max_dr_r) | (dT_T > max_dT_T) | (drho_rho > max_drho_rho) | (du_u > max_du_u)
            n_violations = np.sum(violations)

            if n_violations == 0:
                print(f"  ✓ All gradients within limits after {iteration} iteration(s)")
                break

            if iteration == 0:
                print(f"  Iteration {iteration+1}: {n_violations} gradient violations detected")

            # Add zones at violation locations
            # Insert midpoint between zones with violations
            violation_indices = np.where(violations)[0]
            new_radii = []
            for idx in violation_indices:
                # Add midpoint between zone idx and idx+1
                r_mid = np.sqrt(r_new[idx] * r_new[idx+1])  # Geometric mean
                new_radii.append(r_mid)

            if new_radii:
                # Check if adding these would exceed max
                if len(r_new) + len(new_radii) > nzon_max:
                    # Only add enough to reach the limit
                    n_add = nzon_max - len(r_new)
                    if n_add > 0:
                        new_radii = new_radii[:n_add]
                        r_new = np.sort(np.concatenate([r_new, new_radii]))
                        print(f"    Added {len(new_radii)} zones (reached max {nzon_max})")
                    else:
                        print(f"    Cannot add more zones (at max {nzon_max})")
                    break
                else:
                    # Insert new radii and sort
                    r_new = np.sort(np.concatenate([r_new, new_radii]))
                    print(f"    Added {len(new_radii)} zones at high-gradient regions")
            else:
                break

        print(f"  Final grid size: {len(r_new)} zones")
        
        # Additional feature-based refinement
        if len(r_new) < nzon_max:
            print(f"\n[Feature-Based Refinement] Adding zones at special features...")
            
            # Interpolate quantities for feature detection
            rho_new = interp_logr(rho, r_new)
            p_new   = interp_logr(p, r_new)
            u_new   = interp_logr(u, r_new)
            fH_new  = np.clip(interp_logr(fH, r_new),  0, 1)
            fHe_new = np.clip(interp_logr(fHe, r_new), 0, 1)
            fC_new  = np.clip(interp_logr(fC, r_new),  0, 1)
            fO_new  = np.clip(interp_logr(fO, r_new),  0, 1)
            fSi_new = np.clip(interp_logr(fSi, r_new), 0, 1)
            fFe_new = np.clip(interp_logr(fFe, r_new), 0, 1)
            
            # Renormalise mass fractions
            mf_sum = fH_new + fHe_new + fC_new + fO_new + fSi_new + fFe_new
            mf_sum = np.where(mf_sum < 1e-10, 1.0, mf_sum)
            fH_new  /= mf_sum
            fHe_new /= mf_sum
            fC_new  /= mf_sum
            fO_new  /= mf_sum
            fSi_new /= mf_sum
            fFe_new /= mf_sum
            
            # Compute temperature
            T_new = temperature_from_pressure(p_new, rho_new, fH=fH_new, fHe=fHe_new,
                                              fC=fC_new, fO=fO_new, fSi=fSi_new, fFe=fFe_new)
            
            feature_radii = []
            
            # Composition interfaces
            if detect_composition:
                comp_indices = detect_composition_interfaces(r_new, fH_new, fHe_new, fC_new, 
                                                             fO_new, fSi_new, fFe_new, threshold=0.1)
                if len(comp_indices) > 0:
                    print(f"  Composition interfaces: {len(comp_indices)} detected")
                    for idx in comp_indices:
                        r_mid = np.sqrt(r_new[idx] * r_new[idx+1])
                        feature_radii.append(r_mid)
            
            # Shock fronts
            if detect_shock_features:
                shock_indices = detect_shocks(r_new, rho_new, p_new, u_new, 
                                             pressure_threshold=0.5, velocity_threshold=0.3)
                if len(shock_indices) > 0:
                    print(f"  Shock fronts: {len(shock_indices)} detected")
                    for idx in shock_indices:
                        r_mid = np.sqrt(r_new[idx] * r_new[idx+1])
                        feature_radii.append(r_mid)
            
            # Photosphere
            if detect_photo_features:
                photo_indices = detect_photosphere(r_new, rho_new, T_new, kappa_guess=0.3)
                if len(photo_indices) > 0:
                    print(f"  Photosphere region: {len(photo_indices)} zones")
                    # Add zones at photosphere boundaries
                    if len(photo_indices) > 0:
                        r_mid = np.sqrt(r_new[photo_indices[0]] * r_new[photo_indices[-1]])
                        feature_radii.append(r_mid)
            
            # High curvature regions
            if detect_curvature:
                curv_indices = detect_high_curvature(r_new, rho_new, threshold=0.5)
                if len(curv_indices) > 0:
                    print(f"  High curvature regions: {len(curv_indices)} detected")
                    for idx in curv_indices:
                        r_mid = np.sqrt(r_new[idx] * r_new[idx+1])
                        feature_radii.append(r_mid)
            
            # Add feature zones if space available
            if feature_radii:
                feature_radii = np.unique(feature_radii)
                if len(r_new) + len(feature_radii) > nzon_max:
                    n_add = nzon_max - len(r_new)
                    if n_add > 0:
                        feature_radii = feature_radii[:n_add]
                        r_new = np.sort(np.concatenate([r_new, feature_radii]))
                        print(f"  Added {len(feature_radii)} feature zones (reached max {nzon_max})")
                    else:
                        print(f"  Cannot add feature zones (at max {nzon_max})")
                else:
                    r_new = np.sort(np.concatenate([r_new, feature_radii]))
                    print(f"  Added {len(feature_radii)} feature zones")
            
            print(f"  Final grid size after feature refinement: {len(r_new)} zones")

    # Final interpolation
    rho_new = interp_logr(rho, r_new)
    p_new   = interp_logr(p, r_new)
    u_new   = interp_logr(u, r_new)
    fH_new  = np.clip(interp_logr(fH, r_new),  0, 1)
    fHe_new = np.clip(interp_logr(fHe, r_new), 0, 1)
    fC_new  = np.clip(interp_logr(fC, r_new),  0, 1)
    fO_new  = np.clip(interp_logr(fO, r_new),  0, 1)
    fSi_new = np.clip(interp_logr(fSi, r_new), 0, 1)
    fFe_new = np.clip(interp_logr(fFe, r_new), 0, 1)

    # Renormalise mass fractions
    mf_sum = fH_new + fHe_new + fC_new + fO_new + fSi_new + fFe_new
    mf_sum = np.where(mf_sum < 1e-10, 1.0, mf_sum)
    fH_new  /= mf_sum
    fHe_new /= mf_sum
    fC_new  /= mf_sum
    fO_new  /= mf_sum
    fSi_new /= mf_sum
    fFe_new /= mf_sum

    # Compute temperature for final gradient check
    T_new = temperature_from_pressure(p_new, rho_new, fH=fH_new, fHe=fHe_new,
                                      fC=fC_new, fO=fO_new, fSi=fSi_new, fFe=fFe_new)

    # Final gradient report
    if check_gradients:
        print(f"\n[Final Gradient Check] Analyzing {len(r_new)} zones...")
        
        dr_r = np.abs(np.diff(r_new) / r_new[:-1])
        dT_T = np.abs(np.diff(T_new) / np.maximum(T_new[:-1], 1e-10))
        drho_rho = np.abs(np.diff(rho_new) / np.maximum(rho_new[:-1], 1e-30))
        du_u = np.abs(np.diff(u_new) / np.maximum(np.abs(u_new[:-1]), 1e-10))
        
        n_dr_r = np.sum(dr_r > max_dr_r)
        n_dT_T = np.sum(dT_T > max_dT_T)
        n_drho_rho = np.sum(drho_rho > max_drho_rho)
        n_du_u = np.sum(du_u > max_du_u)
        
        print(f"  Radius jumps > {max_dr_r*100:.0f}%: {n_dr_r} zones")
        print(f"  Temperature jumps > {max_dT_T*100:.0f}%: {n_dT_T} zones")
        print(f"  Density jumps > {max_drho_rho*100:.0f}%: {n_drho_rho} zones")
        print(f"  Velocity jumps > {max_du_u*100:.0f}%: {n_du_u} zones")
        
        total_violations = n_dr_r + n_dT_T + n_drho_rho + n_du_u
        if total_violations == 0:
            print(f"  ✓ All gradients within acceptable limits")
        else:
            print(f"  ⚠ {total_violations} gradient violations remain")
            print(f"    Consider increasing max_iterations or nzon for better resolution")

    return r_new, rho_new, p_new, u_new, fH_new, fHe_new, fC_new, fO_new, fSi_new, fFe_new


def write_hyd(fname: str, r_cm, rho, T, u_cms, t_start_s: float,
              mass_cut_msun: float, rcen_cm: float, rho_cen: float):
    """Write STELLA *.hyd file.

    Note: u_cms is expected to be in cm/s (from input_profile u_cgs).
    """
    Msun = 1.989e33
    n = len(r_cm)

    # Cumulative Lagrangian mass
    dr = np.concatenate([[r_cm[1] - r_cm[0]], np.diff(r_cm)])
    r_lo = np.concatenate([[rcen_cm], r_cm[:-1]])
    r_hi = r_cm
    dm_cgs = (4.0 / 3.0) * np.pi * (r_hi**3 - r_lo**3) * rho
    dm_cgs = np.maximum(dm_cgs, 0.0)
    M_r = np.cumsum(dm_cgs)

    with open(fname, 'w') as f:
        # First line: timetobegin  Nzone  BM1  Rcen  rhoCen
        # Format: (1p,e12.3,i6,3e13.5)
        # timetobegin: time from explosion in seconds (0 = pre-explosion, bomb will trigger)
        # BM1: NS core mass (may be zero)
        # Rcen: core radius (must match rcen_cm used for zone masses)
        # rhoCen: density to left of first rho on grid (must be nonzero)
        t_start_use = max(t_start_s, 0.0)  # Allow t_start=0 for bomb mode
        f.write(f"{t_start_use:12.3E}{n:6d}{mass_cut_msun:13.5E}{rcen_cm:13.5E}{rho_cen:13.5E}\n")
        for km in range(n):
            dMr = dm_cgs[km] / Msun
            Mr  = M_r[km] / Msun + mass_cut_msun
            # Format: (1x,i4,1p,e12.4,e15.7,5e12.4)
            # km  dMr  rHyd  rhoHyd  TpHyd  uHyd  aMr  dum
            f.write(f" {km+1:4d}{dMr:12.4E}{r_cm[km]:15.7E}{rho[km]:12.4E}{T[km]:12.4E}{u_cms[km]:12.4E}{Mr:12.4E}{0.0:12.4E}\n")
    print(f"Written {fname}  ({n} zones, M_total={M_r[-1]/Msun:.3f} Msun)")


def write_abn(fname: str, fH, fHe, fC, fO, fSi, fFe, r_cm):
    """Write STELLA *.abn file.
    Columns: km  0  0  0  H  He  C  N  O  Ne  Na  Mg  Al  Si  S  Ar  Ca  FePeak  Ni58  Ni56
    We set N=Na=Ne=Mg=Al=S=Ar=Ca=Ni58=Ni56=0 and put FePeak=fFe.
    """
    n = len(fH)
    with open(fname, 'w') as f:
        for km in range(n):
            # Format: (i4,1p,19e10.3)
            # km  dum  dum  dum  H  He  C  N  O  Ne  Na  Mg  Al  Si  S  Ar  Ca  FePeak  Ni58  Ni56
            f.write(f"{km+1:4d}{0.0:10.3E}{0.0:10.3E}{0.0:10.3E}"
                    f"{fH[km]:10.3E}{fHe[km]:10.3E}{fC[km]:10.3E}"
                    f"{0.0:10.3E}{fO[km]:10.3E}{0.0:10.3E}{0.0:10.3E}{0.0:10.3E}{0.0:10.3E}{fSi[km]:10.3E}"
                    f"{0.0:10.3E}{0.0:10.3E}{0.0:10.3E}{fFe[km]:10.3E}{0.0:10.3E}{0.0:10.3E}\n")
    print(f"Written {fname}  ({n} zones)")


def write_dat(fname: str, model_name: str, target_days: float,
              explosion_age_days: float = 0.0, nout: int = 30,
              snapshot_interval_days: float = 5.0,
              explosion_energy: float = 0.0, explosion_mass: float = 0.0,
              explosion_time: float = 1.0, explosion_delay: float = 0.0,
              initial_kinetic_energy: float = 0.0,
              nstmax: int = 900000, disable_scattering: bool = False,
              eps: float = 0.003, hmax: float = 5e5,
              enable_convection: bool = False, fit_tau: float = 15.0,
              givdtl: bool = True, custom_output_times: list = None):
    """Write STELLA run/strad/*.dat control file for model_name.

    explosion_age_days: age of the SN at start of the .hyd profile (days since explosion).
    target_days: additional days of evolution beyond the explosion age.
    snapshot_interval_days: output snapshots every N days (default: 5.0).
    explosion_energy: EBurst in units of 1e50 erg (0 = no bomb). e.g. 10 = 1e51 erg.
    explosion_mass: AMHT — heated core mass in Msun (used with EBurst).
    explosion_time: tBurst — energy injection duration in seconds.
    initial_kinetic_energy: EKO in units of 1e50 erg (pre-existing kinetic energy of profile).
    nstmax: maximum number of timesteps before STELLA halts (default 900000).
    TcurA/TcurB in STELLA .dat are in real days; STELLA stops when t*UTIME >= TcurB*86400.
    """
    t_start = explosion_age_days  # days at start of .hyd (= timetobegin/86400)
    t_end = explosion_age_days + target_days  # absolute day to stop

    # Generate output epochs in absolute days, spaced from t_start to t_end
    # STELLA TO array has KOMAX=80 slots — must not exceed this
    KOMAX = 80
    
    # Use custom output times if provided
    if custom_output_times is not None and len(custom_output_times) > 0:
        epochs = custom_output_times
        if len(epochs) > KOMAX:
            print(f"WARNING: {len(epochs)} custom output times exceed KOMAX={KOMAX}, truncating")
            epochs = epochs[:KOMAX]
    # For 1-day simulations, use specific epochs
    elif abs(target_days - 1.0) < 0.1:
        epochs = [0.2, 0.4, 0.6, 0.8, 1.0]
    # For 15-day simulations (2 days past SBO at 13 days), use specific epochs
    elif abs(target_days - 15.0) < 0.1:
        epochs = [0.4, 0.8, 1.2, 1.6, 2.0, 4.0, 8.0, 15.0]
    else:
        n_epochs = min(KOMAX, max(5, int(target_days / snapshot_interval_days + 1)))
        epochs = [round(t_start + (i + 1) * target_days / n_epochs, 4)
                  for i in range(n_epochs)]
        epochs = sorted(set(epochs))

    with open(fname, 'w') as f:
        f.write(f"   EPS   (STIFF ACCURACY)   HMIN (MIN.STEP(SEC)) HMAX (MAX.STEP(SEC))\n")
        f.write(f"   {eps:.4g}                    1.D-14               {hmax:.4g}\n")
        f.write(f"   METH  ADAMS(1)-GEAR(2)-GEAR+BGH(3)   JSTART  (0-START_ELSE_ORDER) MAXORD KNadap\n")
        f.write(f"    2                        0                            4        -4\n")
        f.write(f"   abs(NSTA)>abs(NSTB)=>reads last   TcurA     TcurB\n")
        f.write(f"   -9000     -1                      {t_start:.4f}   {t_end:.4f}\n")
        f.write(f"   NSTMAX  NDebug  NOUT    IOUT (0 -full <0 print zones)   MBATCH\n")
        f.write(f"   {nstmax}   100000  {len(epochs)}      10                             900\n")
        f.write(f"   AMNI (CONTAM. MASS)     XMNI (56NI MASS)\n")
        f.write(f"   8.D-01                  0.0d0\n")
        f.write(f"   AMHT(Solar)  EBurst(1e50) tBurst(s)  tstart -- Heated Core, Energy & time\n")
        f.write(f"   {explosion_mass:.4g}d0      {explosion_energy:.4g}d0      {explosion_time:.4g}d0      {explosion_delay:.4g}d0\n")
        f.write(f"   EKO (kin.energy,1e50) mass fract.tri. u profile  us (out +1)\n")
        f.write(f"   {initial_kinetic_energy:.4g}d0                      0.d0  9.d-1  1.d0          +1.\n")
        conv_flag = "T" if enable_convection else "F"
        givdtl_flag = "T" if givdtl else "F"
        f.write(f"   THRMAT   CRAP     CONV     EDTM     CHNCND  Givdtl\n")
        f.write(f"   1.D-30   1.D0     {conv_flag}        T        T       {givdtl_flag}\n")
        # Temperature floor: 3000 K for all runs (physically motivated minimum)
        temp_floor = "3.D-02"  # 3000 K
        
        f.write(f"   FLOOR(R) FLOOR(V) FLOOR(T) FLOOR(RADIAT)    FLOOR(UCONV)\n")
        f.write(f"   1.D-12   1.D-03   {temp_floor}   1.D-02           1.D-04\n")
        f.write(f"   Wacc(R)  Wacc(V)  Wacc(T)  Wacc(RADIAT)\n")
        f.write(f"   1.d-0    1.d-0    1.D-01   1.D-00\n")
        # Set scattering based on disable_scattering flag
        scat_flag = "F" if disable_scattering else "T"
        f.write(f"   FitTau   TAUTOL   Rvis   AQ    BQ       DRT     NRT     SCAT\n")
        f.write(f"   {fit_tau:.1f}D0   1.3D0     1.D0   2.D0  1.D+00   1.      1       {scat_flag}\n")
        f.write(f"    NTO then    DATA TO\n")
        f.write(f"    {len(epochs)}\n")
        for t_out in epochs:
            f.write(f"    {t_out:.4f}\n")
    print(f"Written {fname}  ({len(epochs)} output epochs, last={epochs[-1]:.2f} days)")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input_profile', help='HERACLES input_profile path')
    parser.add_argument('model_name', help='STELLA model name (used for output filenames)')
    parser.add_argument('output_dir', help='Directory to write *.hyd, *.abn, *.dat')
    parser.add_argument('--nzon', type=int, default=500, help='Target Lagrangian zone count (default: 500)')
    parser.add_argument('--nzon-max', type=int, default=500)
    parser.add_argument('--t-start', type=float, default=0.0, help='Start time in seconds')
    parser.add_argument('--mass-cut', type=float, default=0.0, help='Core mass in Msun')
    parser.add_argument('--rcen', type=float, default=None, help='Inner boundary radius (cm)')
    parser.add_argument('--rho-cen', type=float, default=None, help='Central density (g/cc)')
    parser.add_argument('--r-min', type=float, default=1e7, help='Minimum radius to include (cm) - excludes inner wind (default: 1e7)')
    parser.add_argument('--target-days', type=float, default=15.0,
                        help='Target simulation endpoint (days)')
    parser.add_argument('--snapshot-interval', type=float, default=5.0,
                        help='Output snapshot interval in days (default: 5.0)')
    parser.add_argument('--explosion-energy', type=float, default=0.0,
                        help='EBurst in units of 1e50 erg (0=no bomb, 10=1e51 erg)')
    parser.add_argument('--explosion-mass', type=float, default=0.0,
                        help='AMHT: heated core mass in Msun (default: 0.0)')
    parser.add_argument('--explosion-time', type=float, default=1.0,
                        help='tBurst: energy injection duration in seconds (default: 1.0)')
    parser.add_argument('--explosion-delay', type=float, default=0.0,
                        help='Delay before thermal bomb triggers in seconds (default: 0.0)')
    parser.add_argument('--initial-kinetic-energy', type=float, default=0.0,
                        help='EKO: pre-existing kinetic energy in units of 1e50 erg (default: 0.0)')
    parser.add_argument('--nstmax', type=int, default=900000,
                        help='NSTMAX: maximum number of timesteps (default: 900000)')
    parser.add_argument('--disable-scattering', action='store_true',
                        help='Disable scattering for faster execution')
    parser.add_argument('--eps', type=float, default=0.003,
                        help='STIFF accuracy parameter EPS (default: 0.003)')
    parser.add_argument('--hmax', type=float, default=5e5,
                        help='Maximum timestep HMAX in seconds (default: 5e5 = ~5.8 days)')
    parser.add_argument('--enable-convection', action='store_true',
                        help='Enable convection (CONV=T) for extended progenitors')
    parser.add_argument('--fit-tau', type=float, default=15.0,
                        help='FitTau parameter for radiation transport (default: 15.0)')
    parser.add_argument('--givdtl', action='store_true', default=True,
                        help='Enable detailed output (Givdtl=T, default: True)')
    parser.add_argument('--custom-output-times', type=str, default=None,
                        help='Comma-separated list of custom output times in days (e.g., "3.1,6.2,9.3")')
    parser.add_argument('--outer-fraction', type=float, default=0.3,
                        help='Fraction of zones in outer (log-spaced) region (default: 0.3)')
    parser.add_argument('--transition-radius', type=float, default=1e13,
                        help='Radius (cm) where grid transitions from linear to log spacing (default: 1e13)')
    parser.add_argument('--check-gradients', action='store_true', default=True,
                        help='Enable gradient checking for smooth profiles (default: True)')
    parser.add_argument('--max-dr-r', type=float, default=0.15,
                        help='Max allowed relative radius jump (default: 0.15 = 15%)')
    parser.add_argument('--max-dT-T', type=float, default=0.20,
                        help='Max allowed relative temperature jump (default: 0.20 = 20%)')
    parser.add_argument('--max-drho-rho', type=float, default=0.20,
                        help='Max allowed relative density jump (default: 0.20 = 20%)')
    parser.add_argument('--max-du-u', type=float, default=0.20,
                        help='Max allowed relative velocity jump (default: 0.20 = 20%)')
    parser.add_argument('--adaptive-refinement', action='store_true', default=True,
                        help='Enable adaptive grid refinement in high-gradient regions (default: True)')
    parser.add_argument('--max-iterations', type=int, default=10,
                        help='Maximum adaptive refinement iterations (default: 10)')
    parser.add_argument('--detect-composition', action='store_true', default=True,
                        help='Detect and refine composition interfaces (default: True)')
    parser.add_argument('--detect-shock-features', action='store_true', default=True,
                        help='Detect and refine shock fronts (default: True)')
    parser.add_argument('--detect-photo-features', action='store_true', default=True,
                        help='Detect and refine photosphere region (default: True)')
    parser.add_argument('--detect-curvature', action='store_true', default=True,
                        help='Detect and refine high-curvature regions (default: True)')
    args = parser.parse_args()

    nzon = min(args.nzon, args.nzon_max)
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Reading profile: {args.input_profile}")
    prof = parse_heracles_profile(args.input_profile)

    r   = prof['r_cm']
    rho = prof['rho']
    p   = prof['p']
    u   = prof['u_cms']

    # Filter out inner wind region (low density, causes stiffness in STELLA)
    if args.r_min > 0:
        mask = r >= args.r_min
        n_total = len(r)
        n_keep = np.sum(mask)
        if n_keep < n_total:
            print(f"  Excluding inner region: r < {args.r_min:.2e} cm ({n_total - n_keep} zones -> {n_keep} zones)")
            r = r[mask]
            rho = rho[mask]
            p = p[mask]
            u = u[mask]
            prof['fH'] = prof['fH'][mask]
            prof['fHe'] = prof['fHe'][mask]
            prof['fC'] = prof['fC'][mask]
            prof['fO'] = prof['fO'][mask]
            prof['fSi'] = prof['fSi'][mask]
            prof['fFe'] = prof['fFe'][mask]

    # Thermal bomb mode: force t_start = 0 (stationary profile, STELLA does the explosion)
    if args.explosion_energy > 0:
        args.t_start = 0.0
        print(f"  Thermal bomb mode: t_start forced to 0 (EBurst={args.explosion_energy:.4g} x1e50 erg)")
    # Auto-detect start time from profile header if not explicitly provided
    elif args.t_start == 0.0 and prof.get('explosion_age_days') is not None:
        args.t_start = prof['explosion_age_days'] * 86400.0
        print(f"  Auto-detected explosion age: {prof['explosion_age_days']:.1f} days "
              f"-> t_start = {args.t_start:.3e} s")
    # For EKO mode with stationary profiles, set negative start time to include pre-explosion phase
    elif args.initial_kinetic_energy > 0 and args.t_start == 0.0:
        args.t_start = -20.0 * 86400.0  # Start at -20 days
        print(f"  EKO mode with stationary profile: t_start set to -20 days for pre-explosion phase")
    fH  = prof['fH']
    fHe = prof['fHe']
    fC  = prof['fC']
    fO  = prof['fO']
    fSi = prof['fSi']
    fFe = prof['fFe']

    print(f"Input profile: {len(r)} zones, r=[{r[0]:.2e},{r[-1]:.2e}] cm, "
          f"v=[{u[0]/1e5:.0f},{u[-1]/1e5:.0f}] km/s")

    # If rcen was not explicitly provided and inner zones were filtered out,
    # use the last excluded radius as the inner boundary so zone 1 has a
    # non-zero volume.  Otherwise fall back to r[0] (or 0.0 if that is r[0]).
    if args.rcen is not None:
        rcen = args.rcen
    else:
        # r was already masked above; the original full array is no longer
        # available, but we can infer a sensible rcen from the first kept zone.
        # Use a small offset (1%) to guarantee non-zero volume for zone 1.
        rcen = r[0] * 0.95
    rho_cen = args.rho_cen if args.rho_cen is not None else rho[0]

    if len(r) != nzon:
        print(f"Rezoning {len(r)} -> {nzon} zones (uniform log-r spacing with adaptive refinement)...")
        r, rho, p, u, fH, fHe, fC, fO, fSi, fFe = rezone(
            r, rho, p, u, fH, fHe, fC, fO, fSi, fFe, nzon,
            check_gradients=args.check_gradients,
            max_dr_r=args.max_dr_r,
            max_dT_T=args.max_dT_T,
            max_drho_rho=args.max_drho_rho,
            max_du_u=args.max_du_u,
            adaptive_refinement=args.adaptive_refinement,
            max_iterations=args.max_iterations,
            nzon_max=args.nzon_max,
            detect_composition=args.detect_composition,
            detect_shock_features=args.detect_shock_features,
            detect_photo_features=args.detect_photo_features,
            detect_curvature=args.detect_curvature
        )

    # Inject kinetic energy if specified (EKO mode)
    if args.initial_kinetic_energy > 0:
        Msun = 1.989e33
        E_target = args.initial_kinetic_energy * 1e50  # erg
        dr = np.concatenate([[r[1] - r[0]], np.diff(r)])
        r_lo = np.concatenate([[0.0], r[:-1]])
        r_hi = r
        dm_cgs = (4.0 / 3.0) * np.pi * (r_hi**3 - r_lo**3) * rho
        dm_cgs = np.maximum(dm_cgs, 0.0)
        M_total = np.sum(dm_cgs)
        E_current = 0.5 * np.sum(dm_cgs * u**2)
        if E_current > 0:
            scale_factor = np.sqrt(E_target / E_current)
            u = u * scale_factor
            print(f"  Injecting kinetic energy: EKO={args.initial_kinetic_energy:.4g} x1e50 erg")
            print(f"    Scaling velocities by {scale_factor:.3f} (E_current={E_current/1e50:.3f} x1e50 -> E_target={E_target/1e50:.3f} x1e50)")
        else:
            # Profile is stationary, assign homologous velocities based on radius
            v_scale = np.sqrt(2 * E_target / M_total)
            u = v_scale * (r / r[-1])
            print(f"  Injecting kinetic energy into stationary profile: EKO={args.initial_kinetic_energy:.4g} x1e50 erg")
            print(f"    Assigning homologous velocities: v_max={v_scale/1e5:.0f} km/s")

    # Compute temperature from pressure using actual composition
    T = temperature_from_pressure(p, rho, fH=fH, fHe=fHe, fC=fC, fO=fO, fSi=fSi, fFe=fFe)

    # Parse custom output times if provided
    custom_output_times = None
    if args.custom_output_times:
        try:
            custom_output_times = [float(t.strip()) for t in args.custom_output_times.split(',')]
            print(f"  Using custom output times: {custom_output_times}")
        except ValueError as e:
            print(f"  ERROR: Invalid custom_output_times format: {e}")
            print(f"  Expected comma-separated list of numbers, e.g., '3.1,6.2,9.3'")
            custom_output_times = None

    hyd_file = outdir / f"{args.model_name}.hyd"
    abn_file = outdir / f"{args.model_name}.abn"
    dat_file = outdir / f"{args.model_name}.dat"

    write_hyd(str(hyd_file), r, rho, T, u,
              t_start_s=args.t_start, mass_cut_msun=args.mass_cut,
              rcen_cm=rcen, rho_cen=rho_cen)
    write_abn(str(abn_file), fH, fHe, fC, fO, fSi, fFe, r)
    write_dat(str(dat_file), args.model_name, target_days=args.target_days,
              explosion_age_days=args.t_start / 86400.0,
              snapshot_interval_days=args.snapshot_interval,
              explosion_energy=args.explosion_energy,
              explosion_mass=args.explosion_mass,
              explosion_time=args.explosion_time,
              explosion_delay=args.explosion_delay,
              initial_kinetic_energy=args.initial_kinetic_energy,
              nstmax=args.nstmax, disable_scattering=args.disable_scattering,
              eps=args.eps, hmax=args.hmax,
              enable_convection=args.enable_convection, fit_tau=args.fit_tau, givdtl=args.givdtl,
              custom_output_times=custom_output_times)

    # -------------------------------------------------------
    # Final step: plot the generated .hyd profile for verification
    # -------------------------------------------------------
    plot_hyd_profile(str(hyd_file), str(outdir / f"{args.model_name}_hyd_profile.png"), input_profile=args.input_profile)

    print(f"\nSTELLA input files written to {outdir}/")
    print(f"  {hyd_file.name}  {abn_file.name}  {dat_file.name}")
    print(f"  NZON={nzon}  (compile STELLA with Mzon={nzon} in zone.inc)")


def plot_hyd_profile(hyd_file: str, output_png: str, input_profile: str = None):
    """Read the generated .hyd file and plot all columns for verification.

    Produces a 2x3 multi-panel figure showing density, temperature, velocity,
    cumulative mass, zone mass increment, and temperature–density phase space.
    If input_profile is provided, overlays the original profile for comparison.
    """
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("WARNING: matplotlib not available; skipping .hyd profile plot")
        return

    # Parse .hyd file
    with open(hyd_file) as f:
        lines = f.readlines()

    header = lines[0].strip()
    print(f"  Plotting .hyd profile: header = {header[:60]}...")

    data = []
    for line in lines[1:]:
        cols = line.strip().split()
        if len(cols) >= 8:
            data.append([float(c) for c in cols[:8]])
    arr = np.array(data)
    km   = arr[:, 0].astype(int)
    dMr  = arr[:, 1]
    r    = arr[:, 2]
    rho  = arr[:, 3]
    T    = arr[:, 4]
    u    = arr[:, 5]
    Mr   = arr[:, 6]

    # Parse input profile if provided
    r_in, rho_in, p_in, u_in, fH_in, fHe_in, fC_in, fO_in, fSi_in, fFe_in = None, None, None, None, None, None, None, None, None, None
    if input_profile and Path(input_profile).exists():
        prof = parse_heracles_profile(input_profile)
        r_in = prof['r_cm']
        rho_in = prof['rho']
        p_in = prof['p']
        u_in = prof['u_cms']
        fH_in = prof['fH']
        fHe_in = prof['fHe']
        fC_in = prof['fC']
        fO_in = prof['fO']
        fSi_in = prof['fSi']
        fFe_in = prof['fFe']
        # Compute temperature for input profile
        T_in = temperature_from_pressure(p_in, rho_in, fH=fH_in, fHe=fHe_in,
                                       fC=fC_in, fO=fO_in, fSi=fSi_in, fFe=fFe_in)

    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    if input_profile and r_in is not None:
        fig.suptitle(f"STELLA .hyd profile vs Input: {Path(hyd_file).name}", fontsize=12)
    else:
        fig.suptitle(f"STELLA .hyd profile: {Path(hyd_file).name}", fontsize=12)

    def _plot(ax, x, y, xlabel, ylabel, logx=True, logy=True, color='blue', label=None, linestyle='-'):
        ax.plot(x, y, lw=1.2, color=color, linestyle=linestyle, label=label)
        if logx:
            ax.set_xscale('log')
        if logy:
            ax.set_yscale('log')
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.grid(True, ls='--', alpha=0.4)
        if label:
            ax.legend()

    # Panel 1: density vs radius
    _plot(axes[0, 0], r, rho, 'Radius [cm]', r'Density [g/cm$^3$]', color='blue', label='STELLA')
    if r_in is not None:
        _plot(axes[0, 0], r_in, rho_in, 'Radius [cm]', r'Density [g/cm$^3$]', color='red', linestyle='--', label='Input')
    
    # Panel 2: temperature vs radius
    _plot(axes[0, 1], r, T, 'Radius [cm]', 'Temperature [K]', color='blue', label='STELLA')
    if r_in is not None:
        _plot(axes[0, 1], r_in, T_in, 'Radius [cm]', 'Temperature [K]', color='red', linestyle='--', label='Input')
    
    # Panel 3: velocity vs radius
    _plot(axes[0, 2], r, np.abs(u) + 1e-10, 'Radius [cm]', '|Velocity| [cm/s]', color='blue', label='STELLA')
    if r_in is not None:
        _plot(axes[0, 2], r_in, np.abs(u_in) + 1e-10, 'Radius [cm]', '|Velocity| [cm/s]', color='red', linestyle='--', label='Input')
    
    # Panel 4: cumulative mass vs radius
    _plot(axes[1, 0], r, Mr, 'Radius [cm]', r'Cumulative Mass [$M_\odot$]', logy=False, color='blue', label='STELLA')
    if r_in is not None:
        # Compute cumulative mass for input profile
        dr_in = np.concatenate([[r_in[1] - r_in[0]], np.diff(r_in)])
        r_lo_in = np.concatenate([[0.0], r_in[:-1]])
        r_hi_in = r_in
        dm_in = (4.0 / 3.0) * np.pi * (r_hi_in**3 - r_lo_in**3) * rho_in
        Mr_in = np.cumsum(dm_in) / 1.989e33  # Convert to Msun
        _plot(axes[1, 0], r_in, Mr_in, 'Radius [cm]', r'Cumulative Mass [$M_\odot$]', logy=False, color='red', linestyle='--', label='Input')
    
    # Panel 5: dMr vs radius
    _plot(axes[1, 1], r, dMr, 'Radius [cm]', r'Zone Mass [$M_\odot$]', color='blue', label='STELLA')
    if r_in is not None:
        _plot(axes[1, 1], r_in, dm_in / 1.989e33, 'Radius [cm]', r'Zone Mass [$M_\odot$]', color='red', linestyle='--', label='Input')
    
    # Panel 6: T vs rho (phase space)
    axes[1, 2].scatter(rho, T, c=km, cmap='viridis', s=10, label='STELLA', alpha=0.7)
    if r_in is not None:
        axes[1, 2].scatter(rho_in, T_in, c='red', s=5, label='Input', alpha=0.3)
    axes[1, 2].set_xscale('log')
    axes[1, 2].set_yscale('log')
    axes[1, 2].set_xlabel(r'Density [g/cm$^3$]')
    axes[1, 2].set_ylabel('Temperature [K]')
    axes[1, 2].legend()
    axes[1, 2].grid(True, ls='--', alpha=0.4)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(output_png, dpi=150)
    plt.close()
    print(f"  Profile plot saved: {output_png}")


if __name__ == '__main__':
    main()
