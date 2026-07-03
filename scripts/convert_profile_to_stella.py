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
    """Read HERACLES input_profile (r, rho, p, u, T, fH, fHe, fO, fSi, fCa, fFe, fNi)."""
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
            if len(cols) < 12:
                continue
            data.append(cols)
    arr = np.array(data)
    return {
        'r_cm':               arr[:, 0],
        'rho':                arr[:, 1],
        'p':                  arr[:, 2],
        'u_cms':              arr[:, 3],
        'T':                  arr[:, 4],  # Temperature from profile
        'fH':                 arr[:, 5],
        'fHe':                arr[:, 6],
        'fO':                 arr[:, 7],
        'fSi':                arr[:, 8],
        'fCa':                arr[:, 9],
        'fFe':                arr[:, 10],
        'fNi':                arr[:, 11],
        'explosion_age_days': explosion_age_days,
    }


def temperature_from_pressure(p_cgs, rho_cgs: np.ndarray,
                               fH=None, fHe=None, fO=None,
                               fSi=None, fCa=None, fFe=None, fNi=None,
                               mu: float = 0.62, T_from_profile=None) -> np.ndarray:
    """Estimate temperature from pressure.

    Uses the larger of radiation temperature (P_rad = a*T^4/3) and ideal-gas
    temperature, with a floor of 6000 K (STELLA needs reasonable initial T).

    If T_from_profile is provided, use it directly (preferred for LTE stability).
    Otherwise compute from pressure.

    If mass-fraction arrays are provided, mu is computed per-zone from the
    composition (neutral-gas approximation).  Otherwise the fallback mu is used.
    """
    # If temperature is provided in profile, use it directly
    if T_from_profile is not None:
        T = np.asarray(T_from_profile, dtype=float)
        T = np.maximum(T, 6000.0)  # Force minimum 6000K for LTE stability
        T = np.minimum(T, 1e8)  # Cap at reasonable maximum
        return T
    
    p = np.asarray(p_cgs, dtype=float)
    rho = np.asarray(rho_cgs, dtype=float)
    k_B = 1.380649e-16    # erg/K
    m_H = 1.6735575e-24   # g
    a_rad = 7.566e-15     # erg/cm^3/K^4

    # Compute mean molecular weight from composition if available
    if fH is not None:
        fH  = np.asarray(fH,  dtype=float)
        fHe = np.asarray(fHe, dtype=float) if fHe is not None else np.zeros_like(fH)
        fO  = np.asarray(fO,  dtype=float) if fO  is not None else np.zeros_like(fH)
        fSi = np.asarray(fSi, dtype=float) if fSi is not None else np.zeros_like(fH)
        fCa = np.asarray(fCa, dtype=float) if fCa is not None else np.zeros_like(fH)
        fFe = np.asarray(fFe, dtype=float) if fFe is not None else np.zeros_like(fH)
        fNi = np.asarray(fNi, dtype=float) if fNi is not None else np.zeros_like(fH)
        mu_arr = 1.0 / (fH/1.0 + fHe/4.0 + fO/16.0 + fSi/28.0 + fCa/40.0 + fFe/56.0 + fNi/58.0)
        # Guard against pathological all-zero compositions
        mu_arr = np.where(np.isfinite(mu_arr), mu_arr, mu)
    else:
        mu_arr = mu

    T_gas = p * mu_arr * m_H / (rho * k_B)
    T_rad = (3.0 * p / a_rad) ** 0.25
    T = np.maximum(T_gas, T_rad)
    T = np.maximum(T, 6000.0)  # Force minimum 6000K for LTE stability
    T = np.minimum(T, 1e8)  # Cap at reasonable maximum
    return T


def detect_composition_interfaces(r, fH, fHe, fO, fSi, fCa, fFe, fNi, threshold=0.1):
    """Detect composition interfaces where mass fractions change rapidly.
    
    Parameters:
    -----------
    r: array
        Radius array (cm)
    fH, fHe, fO, fSi, fCa, fFe, fNi: arrays
        Mass fraction arrays
    threshold: float
        Minimum relative change to trigger refinement (default: 0.1 = 10%)
    
    Returns:
    --------
    interface_indices: array
        Indices where composition interfaces are detected
    """
    interface_indices = []
    
    for f, name in zip([fH, fHe, fO, fSi, fCa, fFe, fNi], ['H', 'He', 'O', 'Si', 'Ca', 'Fe', 'Ni']):
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


def rezone(r, rho, p, u, T, fH, fHe, fO, fSi, fCa, fFe, fNi, nzon_target: int,
           outer_fraction=0.3, outer_bias=2.0, transition_radius=1e13,
           check_gradients=True, max_dr_r=0.15, max_dT_T=0.20, max_drho_rho=0.20,
           max_du_u=0.20, adaptive_refinement=True, max_iterations=10, nzon_max=500,
           detect_composition=True, detect_shock_features=True, detect_photo_features=True,
           detect_curvature=True, track_refinement=True):
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
        return r, rho, p, u, T, fH, fHe, fO, fSi, fCa, fFe, fNi

    # Start with uniform log-r spacing
    # If nzon_target == nzon_max, start with fewer zones to allow AMR to work
    if adaptive_refinement and nzon_target == nzon_max:
        # Start with 80% of target, let AMR add the remaining 20%
        initial_nzon = int(nzon_target * 0.8)
        print(f"  [AMR] Starting with {initial_nzon} zones (will refine to {nzon_max})")
        r_new = np.logspace(np.log10(r[0]), np.log10(r[-1]), initial_nzon)
    else:
        initial_nzon = nzon_target
        r_new = np.logspace(np.log10(r[0]), np.log10(r[-1]), nzon_target)

    # Track refinement locations for visualization
    refinement_locations = []  # Radii where zones were added via AMR
    gradient_violations = []  # Radii where gradient violations were detected

    def interp_logr(x, r_target):
        return np.interp(r_target, r, x)

    # Adaptive refinement loop
    if adaptive_refinement:
        print(f"\n[Adaptive Refinement] Starting with {len(r_new)} zones (max: {nzon_max})...")
        
        for iteration in range(max_iterations):
            # Check if we've reached the maximum
            if len(r_new) >= nzon_max:
                print(f"  Reached maximum zone limit ({nzon_max}), stopping refinement")
                break

            # Interpolate quantities onto current grid
            rho_new = interp_logr(rho, r_new)
            p_new   = interp_logr(p, r_new)
            u_new   = interp_logr(u, r_new)
            T_new   = interp_logr(T, r_new)
            T_new   = np.maximum(T_new, 1000.0)  # Floor for numerical stability
            fH_new  = np.clip(interp_logr(fH, r_new),  0, 1)
            fHe_new = np.clip(interp_logr(fHe, r_new), 0, 1)
            fO_new  = np.clip(interp_logr(fO, r_new),  0, 1)
            fSi_new = np.clip(interp_logr(fSi, r_new), 0, 1)
            fCa_new = np.clip(interp_logr(fCa, r_new), 0, 1)
            fFe_new = np.clip(interp_logr(fFe, r_new), 0, 1)
            fNi_new = np.clip(interp_logr(fNi, r_new), 0, 1)

            # Renormalise mass fractions
            mf_sum = fH_new + fHe_new + fO_new + fSi_new + fCa_new + fFe_new + fNi_new
            mf_sum = np.where(mf_sum < 1e-10, 1.0, mf_sum)
            fH_new  /= mf_sum
            fHe_new /= mf_sum
            fO_new  /= mf_sum
            fSi_new /= mf_sum
            fCa_new /= mf_sum
            fFe_new /= mf_sum
            fNi_new /= mf_sum

            # Use interpolated temperature from profile (preserves original structure)
            # Fallback to pressure-based estimate only if T is not available

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

            # Track violation locations for visualization
            if track_refinement:
                violation_indices = np.where(violations)[0]
                for idx in violation_indices:
                    r_violation = (r_new[idx] + r_new[idx+1]) / 2
                    gradient_violations.append(r_violation)

            # Add zones at violation locations
            # Insert midpoint between zones with violations
            violation_indices = np.where(violations)[0]
            new_radii = []
            for idx in violation_indices:
                # Add midpoint between zone idx and idx+1
                r_mid = np.sqrt(r_new[idx] * r_new[idx+1])  # Geometric mean
                new_radii.append(r_mid)
                if track_refinement:
                    refinement_locations.append(r_mid)

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
            T_new   = interp_logr(T, r_new)
            T_new   = np.maximum(T_new, 1000.0)
            fH_new  = np.clip(interp_logr(fH, r_new),  0, 1)
            fHe_new = np.clip(interp_logr(fHe, r_new), 0, 1)
            fO_new  = np.clip(interp_logr(fO, r_new),  0, 1)
            fSi_new = np.clip(interp_logr(fSi, r_new), 0, 1)
            fCa_new = np.clip(interp_logr(fCa, r_new), 0, 1)
            fFe_new = np.clip(interp_logr(fFe, r_new), 0, 1)
            fNi_new = np.clip(interp_logr(fNi, r_new), 0, 1)

            # Renormalise mass fractions
            mf_sum = fH_new + fHe_new + fO_new + fSi_new + fCa_new + fFe_new + fNi_new
            mf_sum = np.where(mf_sum < 1e-10, 1.0, mf_sum)
            fH_new  /= mf_sum
            fHe_new /= mf_sum
            fO_new  /= mf_sum
            fSi_new /= mf_sum
            fCa_new /= mf_sum
            fFe_new /= mf_sum
            fNi_new /= mf_sum
            
            feature_radii = []
            
            # Composition interfaces
            if detect_composition:
                comp_indices = detect_composition_interfaces(r_new, fH_new, fHe_new,
                                                             fO_new, fSi_new, fCa_new, fFe_new, fNi_new, threshold=0.1)
                if len(comp_indices) > 0:
                    print(f"  Composition interfaces: {len(comp_indices)} detected")
                    for idx in comp_indices:
                        if idx < len(r_new) - 1:  # Guard against out-of-bounds
                            r_mid = np.sqrt(r_new[idx] * r_new[idx+1])
                            feature_radii.append(r_mid)
                            if track_refinement:
                                refinement_locations.append(r_mid)
            
            # Shock fronts
            if detect_shock_features:
                shock_indices = detect_shocks(r_new, rho_new, p_new, u_new,
                                             pressure_threshold=0.5, velocity_threshold=0.3)
                if len(shock_indices) > 0:
                    print(f"  Shock fronts: {len(shock_indices)} detected")
                    for idx in shock_indices:
                        if idx < len(r_new) - 1:  # Guard against out-of-bounds
                            r_mid = np.sqrt(r_new[idx] * r_new[idx+1])
                            feature_radii.append(r_mid)
                            if track_refinement:
                                refinement_locations.append(r_mid)
            
            # Photosphere
            if detect_photo_features:
                photo_indices = detect_photosphere(r_new, rho_new, T_new, kappa_guess=0.3)
                if len(photo_indices) > 0:
                    print(f"  Photosphere region: {len(photo_indices)} zones")
                    # Add zones at photosphere boundaries
                    if len(photo_indices) > 0 and photo_indices[-1] < len(r_new) - 1:
                        r_mid = np.sqrt(r_new[photo_indices[0]] * r_new[photo_indices[-1]])
                        feature_radii.append(r_mid)
                        if track_refinement:
                            refinement_locations.append(r_mid)
            
            # High curvature regions
            if detect_curvature:
                curv_indices = detect_high_curvature(r_new, rho_new, threshold=0.5)
                if len(curv_indices) > 0:
                    print(f"  High curvature regions: {len(curv_indices)} detected")
                    for idx in curv_indices:
                        if idx < len(r_new) - 1:  # Guard against out-of-bounds
                            r_mid = np.sqrt(r_new[idx] * r_new[idx+1])
                            feature_radii.append(r_mid)
                            if track_refinement:
                                refinement_locations.append(r_mid)
            
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

    # Pad to exact nzon_max if needed: fill largest gaps first
    if len(r_new) < nzon_max:
        n_pad = nzon_max - len(r_new)
        print(f"\n[Gap Filling] Padding from {len(r_new)} to {nzon_max} zones ({n_pad} added)")
        for _ in range(n_pad):
            gaps = np.diff(np.log10(r_new))
            max_gap_idx = np.argmax(gaps)
            r_mid = np.sqrt(r_new[max_gap_idx] * r_new[max_gap_idx + 1])
            r_new = np.sort(np.append(r_new, r_mid))
            if track_refinement:
                refinement_locations.append(r_mid)

    # Final interpolation
    rho_new = interp_logr(rho, r_new)
    p_new   = interp_logr(p, r_new)
    u_new   = interp_logr(u, r_new)
    T_new   = interp_logr(T, r_new)
    T_new   = np.maximum(T_new, 1000.0)
    fH_new  = np.clip(interp_logr(fH, r_new),  0, 1)
    fHe_new = np.clip(interp_logr(fHe, r_new), 0, 1)
    fO_new  = np.clip(interp_logr(fO, r_new),  0, 1)
    fSi_new = np.clip(interp_logr(fSi, r_new), 0, 1)
    fCa_new = np.clip(interp_logr(fCa, r_new), 0, 1)
    fFe_new = np.clip(interp_logr(fFe, r_new), 0, 1)
    fNi_new = np.clip(interp_logr(fNi, r_new), 0, 1)

    # Renormalise mass fractions
    mf_sum = fH_new + fHe_new + fO_new + fSi_new + fCa_new + fFe_new + fNi_new
    mf_sum = np.where(mf_sum < 1e-10, 1.0, mf_sum)
    fH_new  /= mf_sum
    fHe_new /= mf_sum
    fO_new  /= mf_sum
    fSi_new /= mf_sum
    fCa_new /= mf_sum
    fFe_new /= mf_sum
    fNi_new /= mf_sum

    # Temperature already interpolated from profile above
    # Ensure it's sane for gradient checking
    T_new = np.maximum(T_new, 1000.0)

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

        # Track final gradient violation locations by type for visualization
        if track_refinement and total_violations > 0:
            # Track violations by type with dictionaries
            gradient_violations_by_type = {
                'radius': [],
                'temperature': [],
                'density': [],
                'velocity': []
            }

            # Radius violations
            radius_violations = dr_r > max_dr_r
            radius_indices = np.where(radius_violations)[0]
            for idx in radius_indices:
                r_viol = (r_new[idx] + r_new[idx+1]) / 2
                gradient_violations_by_type['radius'].append(r_viol)

            # Temperature violations
            temp_violations = dT_T > max_dT_T
            temp_indices = np.where(temp_violations)[0]
            for idx in temp_indices:
                r_viol = (r_new[idx] + r_new[idx+1]) / 2
                gradient_violations_by_type['temperature'].append(r_viol)

            # Density violations
            density_violations = drho_rho > max_drho_rho
            density_indices = np.where(density_violations)[0]
            for idx in density_indices:
                r_viol = (r_new[idx] + r_new[idx+1]) / 2
                gradient_violations_by_type['density'].append(r_viol)

            # Velocity violations
            velocity_violations = du_u > max_du_u
            velocity_indices = np.where(velocity_violations)[0]
            for idx in velocity_indices:
                r_viol = (r_new[idx] + r_new[idx+1]) / 2
                gradient_violations_by_type['velocity'].append(r_viol)

            # Store as part of gradient_violations (will be unpacked in plotting)
            gradient_violations = gradient_violations_by_type
        else:
            gradient_violations = None

        # If violations remain and we have room for more zones, continue AMR
        if total_violations > 0 and len(r_new) < nzon_max and adaptive_refinement:
            print(f"  [AMR] Continuing refinement to resolve {total_violations} violations...")
            remaining_iterations = max_iterations - (max_iterations if iteration else 0)
            # Continue AMR loop with remaining iterations
            for extra_iter in range(remaining_iterations):
                if len(r_new) >= nzon_max:
                    print(f"    Reached maximum zone limit ({nzon_max})")
                    break

                # Interpolate quantities
                rho_new = interp_logr(rho, r_new)
                p_new   = interp_logr(p, r_new)
                u_new   = interp_logr(u, r_new)
                fH_new  = np.clip(interp_logr(fH, r_new),  0, 1)
                fHe_new = np.clip(interp_logr(fHe, r_new), 0, 1)
                fO_new  = np.clip(interp_logr(fO, r_new),  0, 1)
                fSi_new = np.clip(interp_logr(fSi, r_new), 0, 1)
                fCa_new = np.clip(interp_logr(fCa, r_new), 0, 1)
                fFe_new = np.clip(interp_logr(fFe, r_new), 0, 1)
                fNi_new = np.clip(interp_logr(fNi, r_new), 0, 1)

                # Renormalise mass fractions
                mf_sum = fH_new + fHe_new + fO_new + fSi_new + fCa_new + fFe_new + fNi_new
                mf_sum = np.where(mf_sum < 1e-10, 1.0, mf_sum)
                fH_new  /= mf_sum
                fHe_new /= mf_sum
                fO_new  /= mf_sum
                fSi_new /= mf_sum
                fCa_new /= mf_sum
                fFe_new /= mf_sum
                fNi_new /= mf_sum

                # Temperature already interpolated from profile above
                T_new = np.maximum(T_new, 1000.0)

                # Check for gradient violations
                dr_r = np.abs(np.diff(r_new) / r_new[:-1])
                dT_T = np.abs(np.diff(T_new) / np.maximum(T_new[:-1], 1e-10))
                drho_rho = np.abs(np.diff(rho_new) / np.maximum(rho_new[:-1], 1e-30))
                du_u = np.abs(np.diff(u_new) / np.maximum(np.abs(u_new[:-1]), 1e-10))

                violations = (dr_r > max_dr_r) | (dT_T > max_dT_T) | (drho_rho > max_drho_rho) | (du_u > max_du_u)
                n_violations = np.sum(violations)

                if n_violations == 0:
                    print(f"    ✓ All gradients resolved after adding {len(r_new) - initial_nzon} zones")
                    break

                # Add zones at violation locations
                violation_indices = np.where(violations)[0]
                new_radii = []
                for idx in violation_indices:
                    r_mid = np.sqrt(r_new[idx] * r_new[idx+1])
                    new_radii.append(r_mid)
                    if track_refinement:
                        refinement_locations.append(r_mid)

                if new_radii:
                    if len(r_new) + len(new_radii) > nzon_max:
                        n_add = nzon_max - len(r_new)
                        if n_add > 0:
                            new_radii = new_radii[:n_add]
                            r_new = np.sort(np.concatenate([r_new, new_radii]))
                            print(f"    Added {len(new_radii)} zones (reached max {nzon_max})")
                        else:
                            print(f"    Cannot add more zones (at max {nzon_max})")
                        break
                    else:
                        r_new = np.sort(np.concatenate([r_new, new_radii]))
                        print(f"    Added {len(new_radii)} zones (total: {len(r_new)})")
                else:
                    break

            # Recompute final gradient check
            rho_new = interp_logr(rho, r_new)
            p_new   = interp_logr(p, r_new)
            u_new   = interp_logr(u, r_new)
            fH_new  = np.clip(interp_logr(fH, r_new),  0, 1)
            fHe_new = np.clip(interp_logr(fHe, r_new), 0, 1)
            fO_new  = np.clip(interp_logr(fO, r_new),  0, 1)
            fSi_new = np.clip(interp_logr(fSi, r_new), 0, 1)
            fCa_new = np.clip(interp_logr(fCa, r_new), 0, 1)
            fFe_new = np.clip(interp_logr(fFe, r_new), 0, 1)
            fNi_new = np.clip(interp_logr(fNi, r_new), 0, 1)

            mf_sum = fH_new + fHe_new + fO_new + fSi_new + fCa_new + fFe_new + fNi_new
            mf_sum = np.where(mf_sum < 1e-10, 1.0, mf_sum)
            fH_new  /= mf_sum
            fHe_new /= mf_sum
            fO_new  /= mf_sum
            fSi_new /= mf_sum
            fCa_new /= mf_sum
            fFe_new /= mf_sum
            fNi_new /= mf_sum

            # Temperature already interpolated from profile above
            T_new = np.maximum(T_new, 1000.0)

            dr_r = np.abs(np.diff(r_new) / r_new[:-1])
            dT_T = np.abs(np.diff(T_new) / np.maximum(T_new[:-1], 1e-10))
            drho_rho = np.abs(np.diff(rho_new) / np.maximum(rho_new[:-1], 1e-30))
            du_u = np.abs(np.diff(u_new) / np.maximum(np.abs(u_new[:-1]), 1e-10))

            n_dr_r = np.sum(dr_r > max_dr_r)
            n_dT_T = np.sum(dT_T > max_dT_T)
            n_drho_rho = np.sum(drho_rho > max_drho_rho)
            n_du_u = np.sum(du_u > max_du_u)

            print(f"  [Final Check after extended AMR]")
            print(f"    Radius jumps > {max_dr_r*100:.0f}%: {n_dr_r} zones")
            print(f"    Temperature jumps > {max_dT_T*100:.0f}%: {n_dT_T} zones")
            print(f"    Density jumps > {max_drho_rho*100:.0f}%: {n_drho_rho} zones")
            print(f"    Velocity jumps > {max_du_u*100:.0f}%: {n_du_u} zones")

            # Update gradient violations for visualization
            if track_refinement and (n_dr_r + n_dT_T + n_drho_rho + n_du_u) > 0:
                gradient_violations_by_type = {
                    'radius': [],
                    'temperature': [],
                    'density': [],
                    'velocity': []
                }

                radius_violations = dr_r > max_dr_r
                radius_indices = np.where(radius_violations)[0]
                for idx in radius_indices:
                    r_viol = (r_new[idx] + r_new[idx+1]) / 2
                    gradient_violations_by_type['radius'].append(r_viol)

                temp_violations = dT_T > max_dT_T
                temp_indices = np.where(temp_violations)[0]
                for idx in temp_indices:
                    r_viol = (r_new[idx] + r_new[idx+1]) / 2
                    gradient_violations_by_type['temperature'].append(r_viol)

                density_violations = drho_rho > max_drho_rho
                density_indices = np.where(density_violations)[0]
                for idx in density_indices:
                    r_viol = (r_new[idx] + r_new[idx+1]) / 2
                    gradient_violations_by_type['density'].append(r_viol)

                velocity_violations = du_u > max_du_u
                velocity_indices = np.where(velocity_violations)[0]
                for idx in velocity_indices:
                    r_viol = (r_new[idx] + r_new[idx+1]) / 2
                    gradient_violations_by_type['velocity'].append(r_viol)

                gradient_violations = gradient_violations_by_type

    # Convert to numpy arrays for return
    refinement_locations = np.array(refinement_locations) if refinement_locations else np.array([])

    return r_new, rho_new, p_new, u_new, T_new, fH_new, fHe_new, fO_new, fSi_new, fCa_new, fFe_new, fNi_new, refinement_locations, gradient_violations


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


def write_abn(fname: str, fH, fHe, fO, fSi, fCa, fFe, fNi, r_cm):
    """Write STELLA *.abn file.
    Columns: km  0  0  0  H  He  C  N  O  Ne  Na  Mg  Al  Si  S  Ar  Ca  FePeak  Ni58  Ni56
    We set C=N=Na=Ne=Mg=Al=S=Ar=Ni58=0 and put Ca=fCa, FePeak=fFe.
    Ni56 is set to 0 here - radioactive Ni56 comes from the .xni file only.
    """
    n = len(fH)
    with open(fname, 'w') as f:
        for km in range(n):
            # Format: (i4,1p,19e10.3)
            # km  dum  dum  dum  H  He  C  N  O  Ne  Na  Mg  Al  Si  S  Ar  Ca  FePeak  Ni58  Ni56
            f.write(f"{km+1:4d}{0.0:10.3E}{0.0:10.3E}{0.0:10.3E}"
                    f"{fH[km]:10.3E}{fHe[km]:10.3E}{0.0:10.3E}"
                    f"{0.0:10.3E}{fO[km]:10.3E}{0.0:10.3E}{0.0:10.3E}{0.0:10.3E}{0.0:10.3E}{fSi[km]:10.3E}"
                    f"{0.0:10.3E}{0.0:10.3E}{fCa[km]:10.3E}{fFe[km]:10.3E}{0.0:10.3E}{0.0:10.3E}\n")
    print(f"Written {fname}  ({n} zones, Ni56=0 - radioactive Ni56 from .xni file only)")


def write_dat(fname: str, model_name: str, target_days: float,
              explosion_age_days: float = 0.0, nout: int = 30,
              snapshot_interval_days: float = 5.0,
              explosion_energy: float = 0.0, explosion_mass: float = 0.0,
              explosion_time: float = 1.0, explosion_delay: float = 0.0,
              initial_kinetic_energy: float = 0.0,
              nstmax: int = 900000, disable_scattering: bool = False,
              eps: float = 0.003, hmin: float = 1e-14, hmax: float = 5e5,
              meth: int = 2, enable_convection: bool = False, fit_tau: float = 15.0,
              tau_tol: float = 1.3,
              givdtl: bool = True, custom_output_times: list = None,
              wacc_r: float = 1.0, wacc_v: float = 1.0, wacc_t: float = 0.5, wacc_rad: float = 1.0,
              temp_floor: float = 0.1, iout: int = 10):
    """Write STELLA run/strad/*.dat control file for model_name.

    explosion_age_days: age of the SN at start of the .hyd profile (days since explosion).
    target_days: additional days of evolution beyond the explosion age.
    snapshot_interval_days: output snapshots every N days (default: 5.0).
    explosion_energy: EBurst in units of 1e50 erg (0 = no bomb). e.g. 10 = 1e51 erg.
    explosion_mass: AMHT — heated core mass in Msun (used with EBurst).
    explosion_time: tBurst — energy injection duration in seconds.
    initial_kinetic_energy: EKO in units of 1e50 erg (pre-existing kinetic energy of profile).
    nstmax: maximum number of timesteps before STELLA halts (default 900000).
    wacc_r: Wacc(R) - radius accuracy weight (default: 1.0).
    wacc_v: Wacc(V) - velocity accuracy weight (default: 1.0).
    wacc_t: Wacc(T) - temperature accuracy weight (default: 0.5 for stability).
    wacc_rad: Wacc(RADIAT) - radiation accuracy weight (default: 1.0).
    temp_floor: FLOOR(T) - minimum temperature in 10^4 K units (default: 0.1 = 1000 K).
    iout: output frequency (1=every timestep, 10=every 10 timesteps, 0=full zone data).
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
        f.write(f"   {eps:.4g}                    {hmin:.4g}               {hmax:.4g}\n")
        f.write(f"   METH  ADAMS(1)-GEAR(2)-GEAR+BGH(3)   JSTART  (0-START_ELSE_ORDER) MAXORD KNadap\n")
        f.write(f"    {meth}                        0                            4        -4\n")
        f.write(f"   abs(NSTA)>abs(NSTB)=>reads last   TcurA     TcurB\n")
        f.write(f"   -9000     -1                      {t_start:.4f}   {t_end:.4f}\n")
        f.write(f"   NSTMAX  NDebug  NOUT    IOUT (0 -full <0 print zones)   MBATCH\n")
        f.write(f"   {nstmax}   100000  {len(epochs)}      {iout}                             900\n")
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
        # Temperature floor in 10^4 K units, written as Fortran D-format
        temp_floor_str = f"{temp_floor:.2g}D+00"

        f.write(f"   FLOOR(R) FLOOR(V) FLOOR(T) FLOOR(RADIAT)    FLOOR(UCONV)\n")
        f.write(f"   1.D-12   1.D-03   {temp_floor_str}   1.D-02           1.D-04\n")
        f.write(f"   Wacc(R)  Wacc(V)  Wacc(T)  Wacc(RADIAT)\n")
        f.write(f"   {wacc_r:.2g}d-0    {wacc_v:.2g}d-0    {wacc_t:.2g}D-00   {wacc_rad:.2g}D-00\n")
        # Set scattering based on disable_scattering flag
        scat_flag = "F" if disable_scattering else "T"
        f.write(f"   FitTau   TAUTOL   Rvis   AQ    BQ       DRT     NRT     SCAT\n")
        f.write(f"   {fit_tau:.1f}D0   {tau_tol:.1f}D0     1.D0   2.D0  1.D+00   1.      1       {scat_flag}\n")
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
    parser.add_argument('--tau-min', type=float, default=1e-4,
                        help='Trim outer zones where cumulative optical depth from surface < tau_min (default: 1e-4)')
    parser.add_argument('--rho-min', type=float, default=0.0,
                        help='Trim outer zones where density < rho_min (g/cm^3) (default: 0 = no trimming)')
    parser.add_argument('--kappa', type=float, default=0.34,
                        help='Opacity (cm^2/g) for tau-min trimming and photosphere detection (default: 0.34)')
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
    parser.add_argument('--hmin', type=float, default=1e-14,
                        help='Minimum timestep HMIN in seconds (default: 1e-14)')
    parser.add_argument('--hmax', type=float, default=5e5,
                        help='Maximum timestep HMAX in seconds (default: 5e5 = ~5.8 days)')
    parser.add_argument('--meth', type=int, default=2, choices=[1, 2, 3],
                        help='Integrator method: 1=Adams, 2=GEAR, 3=GEAR+BGH (default: 2)')
    parser.add_argument('--enable-convection', action='store_true',
                        help='Enable convection (CONV=T) for extended progenitors')
    parser.add_argument('--fit-tau', type=float, default=15.0,
                        help='FitTau parameter for radiation transport (default: 15.0)')
    parser.add_argument('--tau-tol', type=float, default=1.3,
                        help='TauTol: Ncnd changes when tau outside [FitTau/TauTol, FitTau*TauTol] '
                             '(default: 1.3; use 1.1 for scattering runs to avoid large Ncnd jumps)')
    parser.add_argument('--givdtl', action='store_true', default=True,
                        help='Enable detailed output (Givdtl=T, default: True)')
    parser.add_argument('--custom-output-times', type=str, default=None,
                        help='Comma-separated list of custom output times in days (e.g., "3.1,6.2,9.3")')
    parser.add_argument('--iout', type=int, default=10,
                        help='Output frequency: 1=every timestep (full lightcurve), 10=every 10 timesteps (default: 10)')
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
    parser.add_argument('--wacc-r', type=float, default=1.0,
                        help='Wacc(R) - radius accuracy weight (default: 1.0)')
    parser.add_argument('--wacc-v', type=float, default=1.0,
                        help='Wacc(V) - velocity accuracy weight (default: 1.0)')
    parser.add_argument('--wacc-t', type=float, default=0.5,
                        help='Wacc(T) - temperature accuracy weight (default: 0.5 for stability)')
    parser.add_argument('--wacc-rad', type=float, default=1.0,
                        help='Wacc(RADIAT) - radiation accuracy weight (default: 1.0)')
    parser.add_argument('--temp-floor', type=float, default=0.1,
                        help='FLOOR(T) - minimum temperature in 10^4 K units (default: 0.1 = 1000 K)')
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
    fH  = prof['fH']
    fHe = prof['fHe']
    fO  = prof['fO']
    fSi = prof['fSi']
    fCa = prof['fCa']
    fFe = prof['fFe']
    fNi = prof['fNi']

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
            if 'T' in prof:
                prof['T'] = prof['T'][mask]
            fH  = fH[mask]
            fHe = fHe[mask]
            fO  = fO[mask]
            fSi = fSi[mask]
            fCa = fCa[mask]
            fFe = fFe[mask]
            fNi = fNi[mask]

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
    fO  = prof['fO']
    fSi = prof['fSi']
    fCa = prof['fCa']
    fFe = prof['fFe']
    fNi = prof['fNi']

    # Trim outer regions where cumulative optical depth from surface is below threshold
    if args.tau_min is not None and args.tau_min > 0:
        kappa = args.kappa
        dr = np.diff(r)
        dtau_zone = kappa * rho[1:] * dr
        # tau_zone_cum[i] = optical depth from r[i] to the surface (integrated from outside in)
        tau_zone_cum = np.cumsum(dtau_zone[::-1])[::-1]
        keep_zone = tau_zone_cum >= args.tau_min
        if np.any(keep_zone):
            last_keep = np.where(keep_zone)[0][-1]  # outermost zone with enough tau
            n_keep = last_keep + 1  # keep everything from inner boundary to this zone
            if n_keep < len(r):
                n_trim = len(r) - n_keep
                print(f"  Trimming outer {n_trim} zones (tau < {args.tau_min:.2e} at r={r[n_keep-1]:.2e} cm)")
                print(f"    Kept {n_keep} zones, r=[{r[0]:.2e},{r[n_keep-1]:.2e}] cm")
                r = r[:n_keep]
                rho = rho[:n_keep]
                p = p[:n_keep]
                u = u[:n_keep]
                if 'T' in prof:
                    prof['T'] = prof['T'][:n_keep]
                prof['fH'] = prof['fH'][:n_keep]
                prof['fHe'] = prof['fHe'][:n_keep]
                prof['fCa'] = prof['fCa'][:n_keep]
                prof['fO'] = prof['fO'][:n_keep]
                prof['fSi'] = prof['fSi'][:n_keep]
                prof['fFe'] = prof['fFe'][:n_keep]
                prof['fNi'] = prof['fNi'][:n_keep]
            else:
                print(f"  No outer trimming needed (all zones tau >= {args.tau_min:.2e})")
        else:
            print(f"  Warning: entire profile is optically thin (tau < {args.tau_min:.2e}); keeping full grid")

    # Trim outer regions where density falls below threshold
    if args.rho_min is not None and args.rho_min > 0:
        keep_rho = rho >= args.rho_min
        if np.any(keep_rho):
            last_keep = np.where(keep_rho)[0][-1]  # outermost zone dense enough
            n_keep = last_keep + 1
            if n_keep < len(r):
                n_trim = len(r) - n_keep
                print(f"  Trimming outer {n_trim} zones (rho < {args.rho_min:.2e} g/cm^3 at r={r[n_keep]:.2e} cm)")
                print(f"    Kept {n_keep} zones, r=[{r[0]:.2e},{r[n_keep-1]:.2e}] cm")
                r = r[:n_keep]
                rho = rho[:n_keep]
                p = p[:n_keep]
                u = u[:n_keep]
                if 'T' in prof:
                    prof['T'] = prof['T'][:n_keep]
                prof['fH'] = prof['fH'][:n_keep]
                prof['fHe'] = prof['fHe'][:n_keep]
                prof['fCa'] = prof['fCa'][:n_keep]
                prof['fO'] = prof['fO'][:n_keep]
                prof['fSi'] = prof['fSi'][:n_keep]
                prof['fFe'] = prof['fFe'][:n_keep]
                prof['fNi'] = prof['fNi'][:n_keep]
            else:
                print(f"  No outer trimming needed (all zones rho >= {args.rho_min:.2e})")
        else:
            print(f"  Warning: entire profile below rho_min={args.rho_min:.2e}; keeping full grid")

    # Re-assign abundance variables after any trimming so lengths match r
    fH  = prof['fH']
    fHe = prof['fHe']
    fO  = prof['fO']
    fSi = prof['fSi']
    fCa = prof['fCa']
    fFe = prof['fFe']
    fNi = prof['fNi']

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

    # Get temperature from profile BEFORE rezone so it is preserved through grid refinement
    T = prof.get('T', None)
    if T is None:
        # Fallback: compute from pressure if not in profile
        T = temperature_from_pressure(p, rho, fH=fH, fHe=fHe, fO=fO, fSi=fSi, fCa=fCa, fFe=fFe, fNi=fNi)

    if len(r) != nzon:
        print(f"Rezoning {len(r)} -> {nzon} zones (uniform log-r spacing with adaptive refinement)...")
        r, rho, p, u, T, fH, fHe, fO, fSi, fCa, fFe, fNi, refinement_locations, gradient_violations = rezone(
            r, rho, p, u, T, fH, fHe, fO, fSi, fCa, fFe, fNi, nzon,
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
            detect_curvature=args.detect_curvature,
            track_refinement=True
        )
    else:
        refinement_locations = np.array([])
        gradient_violations = None

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

    # Temperature already interpolated by rezone() from profile; ensure it's sane
    if T is None or len(T) != len(r):
        # Fallback: compute from pressure if not available
        T_from_profile = prof.get('T', None)
        T = temperature_from_pressure(p, rho, fH=fH, fHe=fHe, fO=fO, fSi=fSi, fCa=fCa, fFe=fFe, fNi=fNi, T_from_profile=T_from_profile)
    T = np.maximum(T, 1000.0)  # Absolute floor for numerical stability

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
    write_abn(str(abn_file), fH, fHe, fO, fSi, fCa, fFe, fNi, r)
    write_dat(str(dat_file), args.model_name, target_days=args.target_days,
              explosion_age_days=args.t_start / 86400.0,
              snapshot_interval_days=args.snapshot_interval,
              explosion_energy=args.explosion_energy,
              explosion_mass=args.explosion_mass,
              explosion_time=args.explosion_time,
              explosion_delay=args.explosion_delay,
              initial_kinetic_energy=args.initial_kinetic_energy,
              nstmax=args.nstmax, disable_scattering=args.disable_scattering,
              eps=args.eps, hmin=args.hmin, hmax=args.hmax, meth=args.meth,
              enable_convection=args.enable_convection, fit_tau=args.fit_tau,
              tau_tol=args.tau_tol, givdtl=args.givdtl,
              custom_output_times=custom_output_times,
              wacc_r=args.wacc_r, wacc_v=args.wacc_v, wacc_t=args.wacc_t, wacc_rad=args.wacc_rad,
              temp_floor=args.temp_floor, iout=args.iout)

    # -------------------------------------------------------
    # Final step: plot the generated .hyd profile for verification
    # -------------------------------------------------------
    plot_hyd_profile(str(hyd_file), str(outdir / f"{args.model_name}_hyd_profile.png"),
                     input_profile=args.input_profile,
                     refinement_locations=refinement_locations,
                     gradient_violations=gradient_violations)

    print(f"\nSTELLA input files written to {outdir}/")
    print(f"  {hyd_file.name}  {abn_file.name}  {dat_file.name}")
    print(f"  NZON={nzon}  (compile STELLA with Mzon={nzon} in zone.inc)")


def plot_hyd_profile(hyd_file: str, output_png: str, input_profile: str = None,
                     refinement_locations=None, gradient_violations=None):
    """Read the generated .hyd file and plot all columns for verification.

    Produces a 2x3 multi-panel figure showing density, temperature, velocity,
    cumulative mass, zone mass increment, and temperature–density phase space.
    If input_profile is provided, overlays the original profile for comparison.
    If refinement_locations is provided, marks where AMR added zones.
    If gradient_violations is provided, marks where gradient violations were detected.
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
    r_in, rho_in, p_in, u_in, fH_in, fHe_in, fO_in, fSi_in, fCa_in, fFe_in, fNi_in = None, None, None, None, None, None, None, None, None, None, None
    if input_profile and Path(input_profile).exists():
        prof = parse_heracles_profile(input_profile)
        r_in = prof['r_cm']
        rho_in = prof['rho']
        p_in = prof['p']
        u_in = prof['u_cms']
        fH_in = prof['fH']
        fHe_in = prof['fHe']
        fO_in = prof['fO']
        fSi_in = prof['fSi']
        fCa_in = prof['fCa']
        fFe_in = prof['fFe']
        fNi_in = prof['fNi']
        # Compute temperature for input profile
        T_from_profile_in = prof.get('T', None)
        T_in = temperature_from_pressure(p_in, rho_in, fH=fH_in, fHe=fHe_in,
                                       fO=fO_in, fSi=fSi_in, fCa=fCa_in, fFe=fFe_in, fNi=fNi_in,
                                       T_from_profile=T_from_profile_in)

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

    # Mark AMR refinement locations and gradient violations on radius-based panels only
    # (panels 1-5 have radius on x-axis, panel 6 is density vs temperature phase space)
    radius_panels = [axes[0, 0], axes[0, 1], axes[0, 2], axes[1, 0], axes[1, 1]]

    if refinement_locations is not None and len(refinement_locations) > 0:
        for ax in radius_panels:
            for r_ref in refinement_locations:
                ax.axvline(x=r_ref, color='green', linestyle=':', alpha=0.5, linewidth=1)

    # Mark gradient violations by type with different colors
    if gradient_violations is not None and isinstance(gradient_violations, dict):
        # New format: dictionary by type
        violation_colors = {
            'radius': 'orange',
            'temperature': 'red',
            'density': 'purple',
            'velocity': 'cyan'
        }
        violation_labels = {
            'radius': 'Radius',
            'temperature': 'Temperature',
            'density': 'Density',
            'velocity': 'Velocity'
        }

        for ax in radius_panels:
            for vtype, locations in gradient_violations.items():
                if len(locations) > 0:
                    for r_viol in locations:
                        ax.axvline(x=r_viol, color=violation_colors[vtype],
                                  linestyle='--', alpha=0.4, linewidth=1.5)
    elif gradient_violations is not None and len(gradient_violations) > 0:
        # Old format: single array (backward compatibility)
        for ax in radius_panels:
            for r_viol in gradient_violations:
                ax.axvline(x=r_viol, color='orange', linestyle='--', alpha=0.3, linewidth=1)

    # Add legend for AMR markers and gradient violations if any were used
    if (refinement_locations is not None and len(refinement_locations) > 0) or \
       (gradient_violations is not None and (
           (isinstance(gradient_violations, dict) and any(len(v) > 0 for v in gradient_violations.values())) or
           (not isinstance(gradient_violations, dict) and len(gradient_violations) > 0)
       )):
        from matplotlib.lines import Line2D
        legend_elements = []
        if refinement_locations is not None and len(refinement_locations) > 0:
            legend_elements.append(Line2D([0], [0], color='green', linestyle=':', alpha=0.5,
                                         label=f'AMR zones ({len(refinement_locations)})'))

        if gradient_violations is not None:
            if isinstance(gradient_violations, dict):
                # New format: add each violation type
                violation_colors = {
                    'radius': 'orange',
                    'temperature': 'red',
                    'density': 'purple',
                    'velocity': 'cyan'
                }
                violation_labels = {
                    'radius': 'Radius',
                    'temperature': 'Temperature',
                    'density': 'Density',
                    'velocity': 'Velocity'
                }
                for vtype, locations in gradient_violations.items():
                    if len(locations) > 0:
                        legend_elements.append(Line2D([0], [0], color=violation_colors[vtype],
                                                     linestyle='--', alpha=0.4,
                                                     label=f'{violation_labels[vtype]} violations ({len(locations)})'))
            elif len(gradient_violations) > 0:
                # Old format
                legend_elements.append(Line2D([0], [0], color='orange', linestyle='--', alpha=0.3,
                                             label=f'Gradient violations ({len(gradient_violations)})'))

        if legend_elements:
            axes[0, 0].legend(handles=legend_elements, loc='upper right', fontsize=7)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(output_png, dpi=150)
    plt.close()
    print(f"  Profile plot saved: {output_png}")


if __name__ == '__main__':
    main()
