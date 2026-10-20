"""Convert a [DEPRECATED] HERACLES-style initial profile to STELLA *.hyd + *.abn input files.

HERACLES is deprecated and does not fit the science case (interacting Type IIn SNe).
The project now uses the STELLA -> CMFGEN pipeline exclusively. STELLA .hyd is the
active input format; HERACLES profile parsing is retained for backward compatibility.

STELLA input format (see OpenStella/readme.stella):

*.hyd  -- hydrodynamic profile
  Line 1: timeStart(s)  Nzon  mass_cut(Msun)  Rcen(cm)  rhoCen(g/cc)
  Lines 2..Nzon+1: km  dMr(Msun)  r(cm)  rho(g/cc)  T(K)  u(cm/s)  Mr(Msun)  0

*.abn  -- abundance profile (mass fractions per zone)
  Lines 1..Nzon: km  0  0  0  H  He  C  N  O  Ne  Na  Mg  Al  Si  S  Ar  Ca  FePeak  Ni58  Ni56

[DEPRECATED] HERACLES profile columns (from header comment):
  r_cm  rho_cgs  p_cgs  u_cgs  fx_H  fx_He  fx_C  fx_O  fx_Si  fx_Fe
  Note: STELLA .hyd format is the active path for the STELLA -> CMFGEN pipeline.

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
    """Read input_profile (r, rho, p, u, T, fH, fHe, [fC], fO, fSi, fCa, fFe, fNi).

    Supports both 12-column (legacy, C merged into O) and 13-column (C separate)
    formats.  When 13 columns are detected, fC is read from column 7 and fO from
    column 8.  When 12 columns are detected (legacy), fC is set to zeros and fO
    retains the merged O+C value.

    Note: STELLA .hyd is the active input format; this parsing is retained for
    backward compatibility with the STELLA -> CMFGEN pipeline.
    """
    import re
    data = []
    csm_age_days = None  # CSM spatial structure age (not simulation start time)
    with open(fname) as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith('#'):
                # Read CSM age for spatial structure (not simulation start time)
                m = re.search(r'Ejecta B age\s*=\s*([\d.eE+\-]+)', stripped)
                if m:
                    csm_age_days = float(m.group(1))
                continue
            cols = [float(x) for x in stripped.split()]
            if len(cols) < 12:
                continue
            data.append(cols)
    arr = np.array(data)
    ncols = arr.shape[1]
    if ncols >= 13:
        # 13-column format: C is separate (col 7 = fC, col 8 = fO)
        return {
            'r_cm':               arr[:, 0],
            'rho':                arr[:, 1],
            'p':                  arr[:, 2],
            'u_cms':              arr[:, 3],
            'T':                  arr[:, 4],
            'fH':                 arr[:, 5],
            'fHe':                arr[:, 6],
            'fC':                 arr[:, 7],
            'fO':                 arr[:, 8],
            'fSi':                arr[:, 9],
            'fCa':                arr[:, 10],
            'fFe':                arr[:, 11],
            'fNi':                arr[:, 12],
            'csm_age_days':       csm_age_days,
        }
    else:
        # 12-column legacy format: C merged into O (col 7 = fO+C)
        return {
            'r_cm':               arr[:, 0],
            'rho':                arr[:, 1],
            'p':                  arr[:, 2],
            'u_cms':              arr[:, 3],
            'T':                  arr[:, 4],
            'fH':                 arr[:, 5],
            'fHe':                arr[:, 6],
            'fC':                 np.zeros(len(arr)),
            'fO':                 arr[:, 7],
            'fSi':                arr[:, 8],
            'fCa':                arr[:, 9],
            'fFe':                arr[:, 10],
            'fNi':                arr[:, 11],
            'csm_age_days':       csm_age_days,
        }


def temperature_from_pressure(p_cgs, rho_cgs: np.ndarray,
                               fH=None, fHe=None, fC=None, fO=None,
                               fSi=None, fCa=None, fFe=None, fNi=None,
                               mu: float = 0.62, T_from_profile=None) -> np.ndarray:
    """Estimate temperature from pressure.

    Uses the larger of radiation temperature (P_rad = a*T^4/3) and ideal-gas
    temperature. No floor is applied here — the FLOOR(T) parameter in the
    STELLA .dat file controls the minimum temperature at runtime.

    If T_from_profile is provided, use it directly (preferred for LTE stability).
    Otherwise compute from pressure.

    If mass-fraction arrays are provided, mu is computed per-zone from the
    composition (neutral-gas approximation).  Otherwise the fallback mu is used.
    """
    # If temperature is provided in profile, use it directly
    if T_from_profile is not None:
        T = np.asarray(T_from_profile, dtype=float)
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
        fC  = np.asarray(fC,  dtype=float) if fC  is not None else np.zeros_like(fH)
        fO  = np.asarray(fO,  dtype=float) if fO  is not None else np.zeros_like(fH)
        fSi = np.asarray(fSi, dtype=float) if fSi is not None else np.zeros_like(fH)
        fCa = np.asarray(fCa, dtype=float) if fCa is not None else np.zeros_like(fH)
        fFe = np.asarray(fFe, dtype=float) if fFe is not None else np.zeros_like(fH)
        fNi = np.asarray(fNi, dtype=float) if fNi is not None else np.zeros_like(fH)
        mu_arr = 1.0 / (fH/1.0 + fHe/4.0 + fC/12.0 + fO/16.0 + fSi/28.0 + fCa/40.0 + fFe/56.0 + fNi/58.0)
        # Guard against pathological all-zero compositions
        mu_arr = np.where(np.isfinite(mu_arr), mu_arr, mu)
    else:
        mu_arr = mu

    T_gas = p * mu_arr * m_H / (rho * k_B)
    T_rad = (3.0 * p / a_rad) ** 0.25
    T = np.maximum(T_gas, T_rad)
    T = np.minimum(T, 1e8)  # Cap at reasonable maximum
    return T


def detect_composition_interfaces(r, fH, fHe, fO, fSi, fCa, fFe, fNi, threshold=0.1, fC=None):
    """Detect composition interfaces where mass fractions change rapidly.

    Parameters:
    -----------
    r: array
        Radius array (cm)
    fH, fHe, fO, fSi, fCa, fFe, fNi: arrays
        Mass fraction arrays
    fC: array or None
        Carbon mass fraction (optional, for 13-column format)
    threshold: float
        Minimum relative change to trigger refinement (default: 0.1 = 10%)

    Returns:
    --------
    interface_indices: array
        Indices where composition interfaces are detected
    """
    interface_indices = []

    species = [fH, fHe, fO, fSi, fCa, fFe, fNi]
    names = ['H', 'He', 'O', 'Si', 'Ca', 'Fe', 'Ni']
    if fC is not None:
        species.insert(2, fC)
        names.insert(2, 'C')

    for f, name in zip(species, names):
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
           fC=None,
           outer_fraction=0.3, outer_bias=2.0, transition_radius=1e13,
           check_gradients=True, max_dr_r=0.15, max_dT_T=0.20, max_drho_rho=0.20,
           max_du_u=0.20, adaptive_refinement=True, max_iterations=10, nzon_max=500,
           detect_composition=True, detect_shock_features=True, detect_photo_features=True,
           detect_curvature=True, track_refinement=True, min_dr_r=0.005):
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
    min_dr_r: float
        Minimum relative zone spacing dr/r threshold for warnings (default: 0.005 = 0.5%).
        Zones closer than this trigger a warning but are NOT removed, because
        removing zones would make Nzon != Mzon, causing a fatal STELLA crash
        (EVE uses Mzon for state-vector indexing, STRAD uses Nzon).
    """
    n = len(r)
    if fC is None:
        fC = np.zeros(n)
    if n == nzon_target:
        return r, rho, p, u, T, fH, fHe, fO, fSi, fCa, fFe, fNi, fC

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
            T_new   = np.maximum(T_new, 10.0)  # Floor for numerical stability
            fH_new  = np.clip(interp_logr(fH, r_new),  0, 1)
            fHe_new = np.clip(interp_logr(fHe, r_new), 0, 1)
            fC_new  = np.clip(interp_logr(fC, r_new),  0, 1)
            fO_new  = np.clip(interp_logr(fO, r_new),  0, 1)
            fSi_new = np.clip(interp_logr(fSi, r_new), 0, 1)
            fCa_new = np.clip(interp_logr(fCa, r_new), 0, 1)
            fFe_new = np.clip(interp_logr(fFe, r_new), 0, 1)
            fNi_new = np.clip(interp_logr(fNi, r_new), 0, 1)

            # Renormalise mass fractions
            mf_sum = fH_new + fHe_new + fC_new + fO_new + fSi_new + fCa_new + fFe_new + fNi_new
            mf_sum = np.where(mf_sum < 1e-10, 1.0, mf_sum)
            fH_new  /= mf_sum
            fHe_new /= mf_sum
            fC_new  /= mf_sum
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
                # Enforce minimum dr/r: skip if either resulting pair would be too small
                dr_left = (r_mid - r_new[idx]) / r_new[idx]
                dr_right = (r_new[idx+1] - r_mid) / r_mid
                if min_dr_r > 0 and (dr_left < min_dr_r or dr_right < min_dr_r):
                    continue
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
            T_new   = np.maximum(T_new, 10.0)
            fH_new  = np.clip(interp_logr(fH, r_new),  0, 1)
            fHe_new = np.clip(interp_logr(fHe, r_new), 0, 1)
            fC_new  = np.clip(interp_logr(fC, r_new),  0, 1)
            fO_new  = np.clip(interp_logr(fO, r_new),  0, 1)
            fSi_new = np.clip(interp_logr(fSi, r_new), 0, 1)
            fCa_new = np.clip(interp_logr(fCa, r_new), 0, 1)
            fFe_new = np.clip(interp_logr(fFe, r_new), 0, 1)
            fNi_new = np.clip(interp_logr(fNi, r_new), 0, 1)

            # Renormalise mass fractions
            mf_sum = fH_new + fHe_new + fC_new + fO_new + fSi_new + fCa_new + fFe_new + fNi_new
            mf_sum = np.where(mf_sum < 1e-10, 1.0, mf_sum)
            fH_new  /= mf_sum
            fHe_new /= mf_sum
            fC_new  /= mf_sum
            fO_new  /= mf_sum
            fSi_new /= mf_sum
            fCa_new /= mf_sum
            fFe_new /= mf_sum
            fNi_new /= mf_sum

            feature_radii = []

            # Helper: check if inserting r_mid between r_new[idx] and r_new[idx+1]
            # would violate min_dr_r. Returns True if insertion is allowed.
            def _can_insert(idx, r_mid):
                if min_dr_r <= 0:
                    return True
                dr_left = (r_mid - r_new[idx]) / r_new[idx]
                dr_right = (r_new[idx+1] - r_mid) / r_mid
                return dr_left >= min_dr_r and dr_right >= min_dr_r

            # Composition interfaces
            if detect_composition:
                comp_indices = detect_composition_interfaces(r_new, fH_new, fHe_new,
                                                             fO_new, fSi_new, fCa_new, fFe_new, fNi_new, threshold=0.1,
                                                             fC=fC_new)
                if len(comp_indices) > 0:
                    print(f"  Composition interfaces: {len(comp_indices)} detected")
                    for idx in comp_indices:
                        if idx < len(r_new) - 1:  # Guard against out-of-bounds
                            r_mid = np.sqrt(r_new[idx] * r_new[idx+1])
                            if not _can_insert(idx, r_mid):
                                continue
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
                            if not _can_insert(idx, r_mid):
                                continue
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
                        # Check min_dr_r against the actual spanning zones
                        if min_dr_r <= 0 or (
                            (r_mid - r_new[photo_indices[0]]) / r_new[photo_indices[0]] >= min_dr_r
                            and (r_new[photo_indices[-1]] - r_mid) / r_mid >= min_dr_r):
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
                            if not _can_insert(idx, r_mid):
                                continue
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
    # Skip gaps where splitting would violate min_dr_r (too narrow)
    if len(r_new) < nzon_max:
        n_pad = nzon_max - len(r_new)
        print(f"\n[Gap Filling] Padding from {len(r_new)} to {nzon_max} zones ({n_pad} to add)")
        n_actually_added = 0
        for _ in range(n_pad):
            gaps = np.diff(np.log10(r_new))
            max_gap_idx = np.argmax(gaps)
            r_mid = np.sqrt(r_new[max_gap_idx] * r_new[max_gap_idx + 1])
            # Enforce min_dr_r: skip if either resulting pair would be too small
            if min_dr_r > 0:
                dr_left = (r_mid - r_new[max_gap_idx]) / r_new[max_gap_idx]
                dr_right = (r_new[max_gap_idx + 1] - r_mid) / r_mid
                if dr_left < min_dr_r or dr_right < min_dr_r:
                    # This gap is too narrow to split; zero it out and try next
                    gaps[max_gap_idx] = -1
                    # Find next largest gap
                    max_gap_idx = np.argmax(gaps)
                    if gaps[max_gap_idx] < 0:
                        # No more splittable gaps
                        print(f"  No more gaps splittable with min_dr_r={min_dr_r:.4f}; "
                              f"stopped at {len(r_new)} zones")
                        break
                    r_mid = np.sqrt(r_new[max_gap_idx] * r_new[max_gap_idx + 1])
            r_new = np.sort(np.append(r_new, r_mid))
            n_actually_added += 1
            if track_refinement:
                refinement_locations.append(r_mid)
        print(f"  Added {n_actually_added} zones via gap filling")

    # Report minimum dr/r (enforced during AMR — should be >= min_dr_r now).
    # The initial uniform log-r grid may have pairs below min_dr_r if the
    # radius range is very large and nzon is small; these are not inserted
    # by AMR so they persist from the initial grid.
    if min_dr_r > 0:
        dr_over_r = np.diff(r_new) / r_new[:-1]
        n_tight = np.sum(dr_over_r < min_dr_r)
        min_val = dr_over_r.min()
        if n_tight > 0:
            print(f"\n[Min dr/r] INFO: {n_tight} zone pairs have dr/r < {min_dr_r:.4f} "
                  f"(smallest dr/r = {min_val:.6f}). These are from the initial "
                  f"uniform grid; AMR-inserted zones respect min_dr_r={min_dr_r:.4f}.")
        else:
            print(f"\n[Min dr/r] OK: All zone pairs have dr/r >= {min_dr_r:.4f} "
                  f"(smallest dr/r = {min_val:.6f})")

    # CRITICAL: Nzon in the .hyd MUST equal compiled Mzon in zone.inc.
    # EVE indexes the state vector Y as Y(K + 2*Mzon) for temperature,
    # but STRAD reads it as Y(K + 2*Nzon). If Nzon != Mzon, STRAD reads
    # garbage memory and the run crashes immediately (Tp=0, ifail 33).
    # If min_dr_r prevented us from reaching nzon_max, relax the constraint
    # and fill remaining gaps (better to have a few small zones than a crash).
    if len(r_new) < nzon_max:
        n_remaining = nzon_max - len(r_new)
        print(f"\n[Gap Filling] Relaxing min_dr_r to fill {n_remaining} remaining zones")
        for _ in range(n_remaining):
            gaps = np.diff(np.log10(r_new))
            max_gap_idx = np.argmax(gaps)
            r_mid = np.sqrt(r_new[max_gap_idx] * r_new[max_gap_idx + 1])
            r_new = np.sort(np.append(r_new, r_mid))
            if track_refinement:
                refinement_locations.append(r_mid)
    if len(r_new) != nzon_max:
        raise RuntimeError(
            f"FATAL: Final grid has {len(r_new)} zones but compiled Mzon={nzon_max}. "
            f"Nzon must equal Mzon or STELLA will crash due to state-vector index mismatch "
            f"between EVE (uses Mzon) and STRAD (uses Nzon)."
        )

    # Final interpolation
    rho_new = interp_logr(rho, r_new)
    p_new   = interp_logr(p, r_new)
    u_new   = interp_logr(u, r_new)
    T_new   = interp_logr(T, r_new)
    T_new   = np.maximum(T_new, 10.0)
    fH_new  = np.clip(interp_logr(fH, r_new),  0, 1)
    fHe_new = np.clip(interp_logr(fHe, r_new), 0, 1)
    fC_new  = np.clip(interp_logr(fC, r_new),  0, 1)
    fO_new  = np.clip(interp_logr(fO, r_new),  0, 1)
    fSi_new = np.clip(interp_logr(fSi, r_new), 0, 1)
    fCa_new = np.clip(interp_logr(fCa, r_new), 0, 1)
    fFe_new = np.clip(interp_logr(fFe, r_new), 0, 1)
    fNi_new = np.clip(interp_logr(fNi, r_new), 0, 1)

    # Renormalise mass fractions
    mf_sum = fH_new + fHe_new + fC_new + fO_new + fSi_new + fCa_new + fFe_new + fNi_new
    mf_sum = np.where(mf_sum < 1e-10, 1.0, mf_sum)
    fH_new  /= mf_sum
    fHe_new /= mf_sum
    fC_new  /= mf_sum
    fO_new  /= mf_sum
    fSi_new /= mf_sum
    fCa_new /= mf_sum
    fFe_new /= mf_sum
    fNi_new /= mf_sum

    # Temperature already interpolated from profile above
    # Ensure it's sane for gradient checking
    T_new = np.maximum(T_new, 10.0)

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
                T_new = np.maximum(T_new, 10.0)

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
            T_new = np.maximum(T_new, 10.0)

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

    return r_new, rho_new, p_new, u_new, T_new, fH_new, fHe_new, fO_new, fSi_new, fCa_new, fFe_new, fNi_new, fC_new, refinement_locations, gradient_violations


def write_hyd(fname: str, r_cm, rho, T, u_cms, t_start_s: float,
              mass_cut_msun: float, rcen_cm: float, rho_cen: float):
    """Write STELLA *.hyd file.

    Note: u_cms is expected to be in cm/s (from input_profile u_cgs).
    """
    Msun = 1.989e33

    # Fix duplicate or near-duplicate radii that create zero-mass zones
    # (these cause STELLA's Jacobian to become structurally singular)
    # Instead of removing zones (which mismatches STELLA's compiled Mzon),
    # perturb the duplicate radius slightly to give a non-zero mass.
    # Must use >1e-7 relative spacing to be visible in e15.7 format (7 sig digits).
    if len(r_cm) > 1:
        min_rel_spacing = 1e-6  # Minimum relative spacing (visible in e15.7 format)
        for i in range(1, len(r_cm)):
            if r_cm[i] - r_cm[i-1] < r_cm[i-1] * min_rel_spacing:
                # Perturb this zone's radius to be slightly larger than previous
                r_cm[i] = r_cm[i-1] * (1.0 + min_rel_spacing)
                print(f"  WARNING: perturbed zone {i+1} radius to {r_cm[i]:.6e} (was duplicate of zone {i})")

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
        t_start_use = 0.0  # Always start at t=0 (explosion time), CSM structure is pre-expansion
        f.write(f"{t_start_use:12.3E}{n:6d}{mass_cut_msun:13.5E}{rcen_cm:13.5E}{rho_cen:13.5E}\n")
        for km in range(n):
            dMr = dm_cgs[km] / Msun
            Mr  = M_r[km] / Msun + mass_cut_msun
            # Format: (1x,i4,1p,e12.4,e15.7,5e12.4)
            # km  dMr  rHyd  rhoHyd  TpHyd  uHyd  aMr  dum
            f.write(f" {km+1:4d}{dMr:12.4E}{r_cm[km]:15.7E}{rho[km]:12.4E}{T[km]:12.4E}{u_cms[km]:12.4E}{Mr:12.4E}{0.0:12.4E}\n")
    print(f"Written {fname}  ({n} zones, M_total={M_r[-1]/Msun:.3f} Msun)")


def write_abn(fname: str, fH, fHe, fO, fSi, fCa, fFe, fNi, r_cm, fC=None):
    """Write STELLA *.abn file.
    Columns: km  0  0  0  H  He  C  N  O  Ne  Na  Mg  Al  Si  S  Ar  Ca  FePeak  Ni58  Ni56
    C is written as a separate element (fC).  If fC is not provided, C=0
    (legacy behaviour with C merged into O).

    The LAST column (Ni56 / XNifor) is the radioactive Ni56 mass fraction that
    eve2 reads (evesn2.f: read(23,*)...,XNifor(j)) and writes to the .xni file.
    eve2 does NOT read a pre-existing .xni -- it only writes one -- so Ni56
    MUST be placed here in the .abn last column.  We write fNi (from the
    input_profile's Ni column) into this column.
    """
    n = len(fH)
    fNi_arr = np.asarray(fNi, dtype=float) if fNi is not None else np.zeros(n)
    fC_arr  = np.asarray(fC, dtype=float) if fC is not None else np.zeros(n)
    n_ni = int(np.sum(fNi_arr > 1e-10))
    with open(fname, 'w') as f:
        for km in range(n):
            # Format: (i4,1p,19e10.3)
            # km  dum  dum  dum  H  He  C  N  O  Ne  Na  Mg  Al  Si  S  Ar  Ca  FePeak  Ni58  Ni56
            stable = np.array([fH[km], fHe[km], fC_arr[km], fO[km], fSi[km], fCa[km], fFe[km]], dtype=float)
            stable *= (1.0 - fNi_arr[km]) / max(stable.sum(), 1e-30)
            f.write(f"{km+1:4d}{0.0:10.3E}{0.0:10.3E}{0.0:10.3E}"
                    f"{stable[0]:10.3E}{stable[1]:10.3E}{stable[2]:10.3E}"
                    f"{0.0:10.3E}{stable[3]:10.3E}{0.0:10.3E}{0.0:10.3E}{0.0:10.3E}{0.0:10.3E}{stable[4]:10.3E}"
                    f"{0.0:10.3E}{0.0:10.3E}{stable[5]:10.3E}{stable[6]:10.3E}{0.0:10.3E}{fNi_arr[km]:10.3E}\n")
    print(f"Written {fname}  ({n} zones, Ni56/XNifor from input_profile fNi column, {n_ni} zones with Ni56>0)")


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
              temp_floor: float = 0.03, floor_rad: float = 0.01, iout: int = 10,
              aq: float = 2.0, bq: float = 1.0,
              amni: float = 0.0, xmni: float = 0.0,
              tcura_days: float = 0.0,
              r_max_cm: float = 0.0, v_max_cms: float = 0.0):
    """Write STELLA run/strad/*.dat control file for model_name.

    explosion_age_days: age of the SN at start of the .hyd profile (days since explosion).
                       NOTE: This parameter is deprecated - simulation always starts at t=0.
                       The CSM spatial structure represents pre-explosion wind, but timing starts at explosion.
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
    temp_floor: FLOOR(T) - minimum gas temperature in 10^5 K units (default: 0.03 = 3000 K, STELLA UTP=1e5).
    floor_rad: FLOOR(RADIAT) - minimum radiation temperature in 10^5 K units (default: 0.01 = 1000 K).
    iout: output frequency (1=every timestep, 10=every 10 timesteps, 0=full zone data).
    TcurA/TcurB in STELLA .dat are in real days; STELLA stops when t*UTIME >= TcurB*86400.
    """
    t_start = 0.0  # Epoch axis stays anchored to explosion time
    tcura = tcura_days  # TcurA: STELLA clock value at run start (days)

    # TcurB is a PROPER_T stop, but target_days is an OBS.TIME target and
    # STELLA reports t_ob = PROPER_T - r_max(t)/c.  r_max grows with the
    # outermost velocity, so solve TcurB = target + light_time(TcurB)
    # self-consistently; otherwise the run ends early by ~r_max/c.
    c_light = 2.998e10
    light_time_0 = r_max_cm / c_light / 86400.0 if r_max_cm > 0 else 0.0
    if 0.0 < v_max_cms < c_light and light_time_0 > 0.0:
        denom = 1.0 - v_max_cms / c_light
        numer = target_days + light_time_0 - v_max_cms * tcura / c_light
        t_end = numer / denom
        print(f"  Light-travel-time correction: initial={light_time_0:.2f} d "
              f"(r_max={r_max_cm:.3e} cm, v_max={v_max_cms/1e5:.0f} km/s); "
              f"TcurB={t_end:.2f} d for OBS.TIME target={target_days} d")
    else:
        t_end = target_days + light_time_0

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
        f.write(f"   -9000     -1                      {tcura:.4f}   {t_end:.4f}\n")
        f.write(f"   NSTMAX  NDebug  NOUT    IOUT (0 -full <0 print zones)   MBATCH\n")
        # NOUT controls how often PRIBAL writes .swd in the main loop (every NOUT steps).
        # Must be < 200 to avoid per-step PRIBAL calls (stradsep5tt.f:640 IF(NOUT>=200) block),
        # which causes buffer overflow/SIGSEGV.
        # With NOUT=10 and ~23000 timesteps, the .swd fills up (~445MB / 2378 epochs) at ~5d.
        # NOUT=100 gives 10x fewer writes, covering the full 0-9d range with ~238 epochs.
        nout_val = nout if nout is not None else len(epochs)
        f.write(f"   {nstmax}   100000  {nout_val}      {iout}                             900\n")
        f.write(f"   AMNI (CONTAM. MASS)     XMNI (56NI MASS)\n")
        f.write(f"   {amni:.4g}D0                  {xmni:.4g}D0\n")
        f.write(f"   AMHT(Solar)  EBurst(1e50) tBurst(s)  tstart -- Heated Core, Energy & time\n")
        f.write(f"   {explosion_mass:.4g}d0      {explosion_energy:.4g}d0      {explosion_time:.4g}d0      {explosion_delay:.4g}d0\n")
        f.write(f"   EKO (kin.energy,1e50) mass fract.tri. u profile  us (out +1)\n")
        # EKO>0: STELLA OVERWRITES all .hyd velocities with a triangular (linear-in-mass)
        #         profile scaled to EKO. mass_fract_tri (AI1) only controls WHERE the
        #         triangular starts (0=innermost zone), NOT whether it overwrites.
        #         See begtt.f / begradsep.f lines 218-258.
        # EKO=0: STELLA preserves the .hyd velocity profile as-is.
        f.write(f"   {initial_kinetic_energy:.4g}d0                      0.d0  9.d-1  1.d0          +1.\n")
        conv_flag = "T" if enable_convection else "F"
        givdtl_flag = "T" if givdtl else "F"
        f.write(f"   THRMAT   CRAP     CONV     EDTM     CHNCND  Givdtl\n")
        f.write(f"   1.D-30   1.D0     {conv_flag}        T        T       {givdtl_flag}\n")
        # Temperature floor in 10^5 K units (STELLA UTP=1e5 K), written as Fortran D-format
        temp_floor_str = f"{temp_floor:.2g}D+00"

        f.write(f"   FLOOR(R) FLOOR(V) FLOOR(T) FLOOR(RADIAT)    FLOOR(UCONV)\n")
        floor_rad_str = f"{floor_rad:.2g}D+00"
        f.write(f"   1.D-12   1.D-03   {temp_floor_str}   {floor_rad_str}           1.D-04\n")
        f.write(f"   Wacc(R)  Wacc(V)  Wacc(T)  Wacc(RADIAT)\n")
        f.write(f"   {wacc_r:.2g}d-0    {wacc_v:.2g}d-0    {wacc_t:.2g}D-00   {wacc_rad:.2g}D-00\n")
        # Set scattering based on disable_scattering flag
        scat_flag = "F" if disable_scattering else "T"
        f.write(f"   FitTau   TAUTOL   Rvis   AQ    BQ       DRT     NRT     SCAT\n")
        f.write(f"   {fit_tau:.2f}D0   {tau_tol:.2f}D0     1.D0   {aq:.1f}D0  {bq:.1f}D+00   1.      1       {scat_flag}\n")
        f.write(f"    NTO then    DATA TO\n")
        f.write(f"    {len(epochs)}\n")
        for t_out in epochs:
            f.write(f"    {t_out:.4f}\n")
    print(f"Written {fname}  ({len(epochs)} output epochs, last={epochs[-1]:.2f} days)")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input_profile', help='[DEPRECATED] HERACLES input_profile path (STELLA .hyd is the active format)')
    parser.add_argument('model_name', help='STELLA model name (used for output filenames)')
    parser.add_argument('output_dir', help='Directory to write *.hyd, *.abn, *.dat')
    parser.add_argument('--nzon', type=int, default=500, help='Target Lagrangian zone count (default: 500)')
    parser.add_argument('--nzon-max', type=int, default=500)
    parser.add_argument('--t-start', type=float, default=0.0, help='Start time in seconds')
    parser.add_argument('--tcura-days', type=float, default=0.0,
                        help='TcurA: STELLA clock start time in days (default: 0.0)')
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
    parser.add_argument('--aq', type=float, default=2.0,
                        help='Artificial viscosity coefficient AQ for shock capturing (default: 2.0). '
                             'Increase for stronger shocks (e.g., 5.0-10.0 for non-monotonic velocity).')
    parser.add_argument('--bq', type=float, default=1.0,
                        help='Pseudo-viscosity coefficient BQ for Rayleigh-Taylor (default: 1.0)')
    parser.add_argument('--amni', type=float, default=0.0,
                        help='AMNI: contaminated (Ni56 core) mass in Msun (default: 0.0 = auto). '
                             'Set to match nickel56 distribution_mass_msun for KEPLER models.')
    parser.add_argument('--xmni', type=float, default=0.0,
                        help='XMNI: 56Ni mass in Msun (default: 0.0 = auto). '
                             'Set to match nickel56 total_mass_msun for KEPLER models.')
    parser.add_argument('--custom-output-times', type=str, default=None,
                        help='Comma-separated list of custom output times in days (e.g., "3.1,6.2,9.3")')
    parser.add_argument('--iout', type=int, default=10,
                        help='Output frequency: 1=every timestep (full lightcurve), 10=every 10 timesteps (default: 10)')
    parser.add_argument('--nout', type=int, default=100,
                        help='NOUT: .swd snapshot write frequency in timesteps (default: 100). '
                             'Must be < 200. NOUT=10 fills .swd at ~5d; NOUT=100 covers ~50d.')
    parser.add_argument('--outer-fraction', type=float, default=0.3,
                        help='Fraction of zones in outer (log-spaced) region (default: 0.3)')
    parser.add_argument('--transition-radius', type=float, default=1e13,
                        help='Radius (cm) where grid transitions from linear to log spacing (default: 1e13)')
    parser.add_argument('--check-gradients', action='store_true', default=True,
                        help='Enable gradient checking for smooth profiles (default: True)')
    parser.add_argument('--max-dr-r', type=float, default=0.15,
                        help='Max allowed relative radius jump (default: 0.15 = 15%%)')
    parser.add_argument('--min-dr-r', type=float, default=0.01,
                        help='Minimum relative zone spacing dr/r (default: 0.01 = 1%%). '
                             'Zones closer than this are NOT inserted during AMR, preventing '
                             'tiny zones that force small timesteps.')
    parser.add_argument('--max-dT-T', type=float, default=0.20,
                        help='Max allowed relative temperature jump (default: 0.20 = 20%%)')
    parser.add_argument('--max-drho-rho', type=float, default=0.20,
                        help='Max allowed relative density jump (default: 0.20 = 20%%)')
    parser.add_argument('--max-du-u', type=float, default=0.20,
                        help='Max allowed relative velocity jump (default: 0.20 = 20%%)')
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
    parser.add_argument('--temp-floor', type=float, default=0.03,
                        help='FLOOR(T) - minimum temperature in 10^5 K units (default: 0.03 = 3000 K, STELLA UTP=1e5)')
    parser.add_argument('--floor-rad', type=float, default=0.01,
                        help='FLOOR(RADIAT) - minimum radiation temperature in 10^5 K units (default: 0.01 = 1000 K)')
    args = parser.parse_args()

    nzon = min(args.nzon, args.nzon_max)
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Reading profile: {args.input_profile}")
    # [DEPRECATED] HERACLES profile parsing; STELLA .hyd is the active input format
    prof = parse_heracles_profile(args.input_profile)

    r   = prof['r_cm']
    rho = prof['rho']
    p   = prof['p']
    u   = prof['u_cms']
    fH  = prof['fH']
    fHe = prof['fHe']
    fC  = prof.get('fC', np.zeros_like(fH))
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
            prof['fH']  = prof['fH'][mask]
            prof['fHe'] = prof['fHe'][mask]
            prof['fC']  = prof['fC'][mask]
            prof['fO']  = prof['fO'][mask]
            prof['fSi'] = prof['fSi'][mask]
            prof['fCa'] = prof['fCa'][mask]
            prof['fFe'] = prof['fFe'][mask]
            prof['fNi'] = prof['fNi'][mask]
            fH  = prof['fH']
            fHe = prof['fHe']
            fC  = prof['fC']
            fO  = prof['fO']
            fSi = prof['fSi']
            fCa = prof['fCa']
            fFe = prof['fFe']
            fNi = prof['fNi']

    # Thermal bomb mode: force t_start = 0 (stationary profile, STELLA does the explosion)
    if args.explosion_energy > 0:
        args.t_start = 0.0
        print(f"  Thermal bomb mode: t_start forced to 0 (EBurst={args.explosion_energy:.4g} x1e50 erg)")
    # EKO mode: start at t=0 (explosion time); CSM structure is pre-explosion
    # NOTE: EKO>0 causes STELLA to overwrite velocities with a triangular profile.
    # For profiles with correct velocities, use EKO=0 instead.
    elif args.initial_kinetic_energy > 0 and args.t_start == 0.0:
        args.t_start = 0.0  # Start at explosion time (CSM structure is pre-explosion)
        print(f"  EKO mode: t_start=0 (KE={args.initial_kinetic_energy:.4g} x1e50 erg)")
        print(f"  WARNING: STELLA will overwrite .hyd velocities with triangular profile!")
        if prof.get('csm_age_days') is not None:
            print(f"  CSM spatial structure represents {prof['csm_age_days']:.1f} day old wind")
    # CSM structure already represents pre-explosion wind, so simulation should start at t=0
    elif args.t_start == 0.0:
        args.t_start = 0.0
        if prof.get('csm_age_days') is not None:
            print(f"  CSM spatial structure represents {prof['csm_age_days']:.1f} day old wind")
        print(f"  Simulation starts at t=0 (explosion time, EKO=0: velocities from .hyd preserved)")
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
                prof['fC'] = prof['fC'][:n_keep]
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
                prof['fC'] = prof['fC'][:n_keep]
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
    fC  = prof['fC']
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
        T = temperature_from_pressure(p, rho, fH=fH, fHe=fHe, fC=fC, fO=fO, fSi=fSi, fCa=fCa, fFe=fFe, fNi=fNi)

    if len(r) != nzon:
        print(f"Rezoning {len(r)} -> {nzon} zones (uniform log-r spacing with adaptive refinement)...")
        r, rho, p, u, T, fH, fHe, fO, fSi, fCa, fFe, fNi, fC, refinement_locations, gradient_violations = rezone(
            r, rho, p, u, T, fH, fHe, fO, fSi, fCa, fFe, fNi, nzon,
            fC=fC,
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
            track_refinement=True,
            min_dr_r=args.min_dr_r
        )
    else:
        refinement_locations = np.array([])
        gradient_violations = None

    # Inject kinetic energy if specified (EKO mode)
    # WARNING: EKO>0 in the .dat file causes STELLA to OVERWRITE all .hyd velocities
    # with a triangular (linear-in-mass) profile. The scaling below adjusts the .hyd
    # velocities, but STELLA will ignore them and use its own triangular profile.
    # This .hyd scaling is only useful for energy accounting before STELLA overwrites.
    # For interacting systems with physically motivated velocities, use EKO=0.
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
        print(f"  WARNING: EKO={args.initial_kinetic_energy:.4g} x1e50 erg > 0")
        print(f"    STELLA will OVERWRITE all .hyd velocities with a triangular profile!")
        print(f"    The input velocity structure will be LOST. Use EKO=0 to preserve it.")
        if E_current > 0:
            scale_factor = np.sqrt(E_target / E_current)
            u = u * scale_factor
            print(f"  Scaling .hyd velocities by {scale_factor:.3f} for energy accounting")
            print(f"    (E_current={E_current/1e50:.3f} x1e50 -> E_target={E_target/1e50:.3f} x1e50)")
        else:
            # Profile is stationary, assign homologous velocities based on radius
            v_scale = np.sqrt(2 * E_target / M_total)
            u = v_scale * (r / r[-1])
            print(f"  Assigning homologous velocities to stationary profile: v_max={v_scale/1e5:.0f} km/s")

    # Temperature already interpolated by rezone() from profile; ensure it's sane
    if T is None or len(T) != len(r):
        # Fallback: compute from pressure if not available
        T_from_profile = prof.get('T', None)
        T = temperature_from_pressure(p, rho, fH=fH, fHe=fHe, fC=fC, fO=fO, fSi=fSi, fCa=fCa, fFe=fFe, fNi=fNi, T_from_profile=T_from_profile)
    T = np.maximum(T, 10.0)  # Absolute floor for numerical stability

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
    write_abn(str(abn_file), fH, fHe, fO, fSi, fCa, fFe, fNi, r, fC=fC)
    write_dat(str(dat_file), args.model_name, target_days=args.target_days,
              explosion_age_days=0.0,  # Always start at t=0 (explosion time)
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
              temp_floor=args.temp_floor, floor_rad=args.floor_rad, iout=args.iout,
              nout=args.nout,
              aq=args.aq, bq=args.bq,
              amni=args.amni, xmni=args.xmni,
              tcura_days=args.tcura_days,
              r_max_cm=float(r[-1]), v_max_cms=float(u[-1]))

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
    r_in, rho_in, p_in, u_in, fH_in, fHe_in, fC_in, fO_in, fSi_in, fCa_in, fFe_in, fNi_in = None, None, None, None, None, None, None, None, None, None, None, None
    if input_profile and Path(input_profile).exists():
        # [DEPRECATED] HERACLES profile parsing; STELLA .hyd is the active input format
        prof = parse_heracles_profile(input_profile)
        r_in = prof['r_cm']
        rho_in = prof['rho']
        p_in = prof['p']
        u_in = prof['u_cms']
        fH_in = prof['fH']
        fHe_in = prof['fHe']
        fC_in = prof.get('fC', None)
        fO_in = prof['fO']
        fSi_in = prof['fSi']
        fCa_in = prof['fCa']
        fFe_in = prof['fFe']
        fNi_in = prof['fNi']
        # Compute temperature for input profile
        T_from_profile_in = prof.get('T', None)
        T_in = temperature_from_pressure(p_in, rho_in, fH=fH_in, fHe=fHe_in, fC=fC_in,
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
