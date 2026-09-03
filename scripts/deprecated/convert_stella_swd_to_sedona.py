#!/usr/bin/env python3
"""
Convert a STELLA *.swd snapshot to SEDONA .mod format.

This script reads a STELLA .swd file and converts it to the SEDONA .mod format
for direct use with SEDONA radiative transfer, bypassing CMFGEN.

STELLA *.swd file format (ASCII, written by PRIBAL in lbalsw.trf):
  Each epoch block: one line per zone, columns:
    t_day   km   log10(M_ext_Msun)   log10(R_cm)   V_kms
    log10(T_gas_K)   log10(T_rad_K)   log10(rho_gcc)
    log10(P_cgs)   log10(Qv)   log10(eng)   Flum_1e40   WRKX

  The km=1 (outer boundary) line has t_day as first field;
  all interior lines start with 0 or blank in that column.

Usage:
  python convert_stella_swd_to_sedona.py SWD_FILE OUTPUT_FILE [OPTIONS]
    --target-days T    Pick the epoch closest to T days (default: last epoch)
    --tau-max X        Outer optical depth cutoff (default: 30)
    --tau-min X        Inner optical depth cutoff (default: 1e-4)
    --nzon-target N    Target number of output grid points (default: 200)
    --time-days T      Explosion time in days (default: read from swd)
    --kappa K          Electron-scattering opacity cm^2/g (default: 0.34)
    --species LIST     Comma-separated atomic numbers (default: 1,2,6,8,14,26,27,28)
    --homologous       Write correct time and radius for hydro_module="homologous"
    --no-feature-aware-downsample  Disable feature-aware downsampling (use simple log-spacing)
    --velocity-threshold X    Velocity gradient threshold for feature detection (default: 0.05)
    --density-threshold X     Density gradient threshold for feature detection (default: 0.15)
    --temperature-threshold X Temperature gradient threshold for feature detection (default: 0.10)

The output .mod file is in SEDONA format (inner->outer, increasing r).
"""

from __future__ import annotations

import argparse
import struct
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter


# Physical constants
C_LIGHT   = 2.99792458e10   # cm/s
K_BOLTZ   = 1.380649e-16    # erg/K
M_H       = 1.6735575e-24   # g
SIGMA_SB  = 5.670374419e-5  # erg/cm^2/s/K^4
LSUN      = 3.828e33        # erg/s
MSUN      = 1.989e33        # g
AVOGADRO  = 6.02214076e23

# STELLA Natom=15 element order (from azdat.f / YABUN array)
STELLA_ELEM_NAMES = ['H', 'He', 'C', 'N', 'O', 'Ne', 'Na', 'Mg', 'Al', 'Si',
                     'S', 'Ar', 'Ca', 'Fe', 'Ni']
# Atomic masses AZ from STELLA azdat.f (Common/AZZn/)
STELLA_AZ = np.array([1.0, 4.0, 12.0, 14.0, 16.0, 20.0, 23.0, 24.0, 27.0,
                      28.0, 32.0, 40.0, 40.0, 56.0, 58.0])


def normalize_swd_composition(composition: np.ndarray) -> dict[str, np.ndarray]:
    """Normalize STELLA .swd composition from Y_k*A_k to true mass fractions.

    The .swd file writes Yat(ii) = YABUN(ii,Km) * AZ(ii), i.e. number fraction
    times atomic mass.  These values sum to <A> = sum(Y_k * A_k), not 1.
    The true mass fraction is X_k = Y_k * A_k / <A>.

    This function divides each element column by the per-zone sum of all 15
    elements, producing normalized mass fractions that sum to 1 per zone.

    Args:
        composition: array of shape (nzon, 15) with Y_k*A_k values from .swd

    Returns:
        dict mapping element name -> mass fraction array (nzon,), normalized to sum=1
    """
    comp = np.asarray(composition, dtype=float)
    if comp.ndim != 2 or comp.shape[1] != 15:
        raise ValueError(f"Expected composition shape (nzon, 15), got {comp.shape}")
    # Per-zone sum of all 15 Y_k*A_k values = <A>
    total = comp.sum(axis=1, keepdims=True)
    total = np.where(total < 1e-30, 1.0, total)
    mass_fractions = comp / total  # X_k = Y_k * A_k / <A>
    return {name: mass_fractions[:, i] for i, name in enumerate(STELLA_ELEM_NAMES)}


def read_swd(swd_path: str) -> dict[float, dict]:
    """Parse a STELLA *.swd file.

    Returns a dict keyed by t_day -> {'r': [], 'v_kms': [], 'rho': [], 'T': [], 'flux': [], 'km': [], 'composition': []}.
    All values are numpy arrays sorted outer->inner (decreasing r).
    
    STELLA .swd format (from convert_stella_swd_to_cmfgen.py):
    - Column 0: time (non-zero for header lines with km=1 or km=Nzon)
    - Column 1: km (zone index)
    - Column 3: log10(R_cm)
    - Column 4: velocity in 1e8 cm/s (multiply by 1000 to get km/s)
    - Column 5: log10(T_gas_K)
    - Column 7: log10(rho) where rho is in units of 1e-6 g/cc
    - Column 11: Flum*1e-40 (luminosity in 1e40 erg/s)
    - Columns 15-29: composition (15 elements: H, He, C, N, O, Ne, Na, Mg, Al, Si, S, Ar, Ca, Fe, Ni)
    """
    epochs: dict[float, dict] = {}
    current_t: float | None = None
    buf_km     = []
    buf_r      = []
    buf_v      = []
    buf_T      = []
    buf_rho    = []
    buf_flux   = []
    buf_composition = []

    def _flush():
        nonlocal current_t, buf_km, buf_r, buf_v, buf_T, buf_rho, buf_flux, buf_composition
        if current_t is not None and buf_r:
            order = np.argsort(buf_r)[::-1]  # outer->inner
            epochs[current_t] = {
                'km':    np.array(buf_km)[order],
                'r':     np.array(buf_r)[order],
                'v_kms': np.array(buf_v)[order],
                'T':     np.array(buf_T)[order],
                'rho':   np.array(buf_rho)[order],
                'flux':  np.array(buf_flux)[order],
                'composition': np.array(buf_composition)[order] if buf_composition else None,
            }
        buf_km    = []
        buf_r     = []
        buf_v     = []
        buf_T     = []
        buf_rho   = []
        buf_flux  = []
        buf_composition = []

    with open(swd_path) as fh:
        for raw_line in fh:
            line = raw_line.rstrip()
            if not line.strip():
                continue
            cols = line.split()
            if len(cols) < 13:
                continue
            try:
                # Detect epoch header line: non-zero first field when km==1 or km==Nzon
                t_candidate = float(cols[0])
                km          = int(cols[1])
                log_r       = float(cols[3])
                v_1e8_cms   = float(cols[4])  # velocity in 1e8 cm/s
                log_T       = float(cols[5])
                log_rho     = float(cols[7])
                flux_1e40   = float(cols[11])
            except (ValueError, IndexError):
                continue

            if t_candidate != 0.0:
                # Format 181 line (km=1 or km=Nzon): nonzero time in col 0
                # Only start a new epoch when km==1
                if km == 1:
                    _flush()
                    current_t = t_candidate
            else:
                km = int(cols[1])

            r_cm  = 10.0 ** log_r  # log_r is log10(radius in cm)
            T_K   = 10.0 ** log_T
            # STELLA SWD uses internal units: rho is log10(rho / 1e-6 g/cc)
            # Convert to cgs: rho_cgs = 10^log_rho * 1e-6
            rho   = (10.0 ** log_rho) * 1e-6
            flux  = flux_1e40 * 1e40  # Convert to cgs: erg/s
            v_kms = v_1e8_cms * 1000.0  # Convert 1e8 cm/s to km/s
            
            # Composition data (15 values) starts at column 15 (0-indexed)
            composition = []
            if len(cols) >= 30:  # 15 structure columns + 15 composition columns
                for i in range(15):
                    composition.append(float(cols[15 + i]))

            buf_km.append(km)
            buf_r.append(r_cm)
            buf_v.append(v_kms)
            buf_T.append(T_K)
            buf_rho.append(rho)
            buf_flux.append(flux)
            if composition:
                buf_composition.append(composition)

    _flush()
    if not epochs:
        raise ValueError(f"No epochs found in {swd_path}")
    return epochs


def read_swd_composition(swd_data: dict) -> dict[str, np.ndarray]:
    """Read composition from STELLA .swd data dictionary.
    
    Args:
        swd_data: Dictionary from read_swd() containing 'composition' array
        
    Returns dict mapping element name -> mass fraction array (nzon,).
    """
    if 'composition' not in swd_data or swd_data['composition'] is None:
        print(f"WARNING: No composition data found in .swd data")
        return {name: np.zeros(1) for name in STELLA_ELEM_NAMES}
    
    composition = swd_data['composition']
    return normalize_swd_composition(composition)


def calculate_optical_depth(r: np.ndarray, rho: np.ndarray, kappa: float = 0.34) -> np.ndarray:
    """Calculate electron scattering optical depth from outer boundary.
    
    tau = integral from r to infinity of kappa * rho dr
    Approximate as cumulative sum from outer boundary inward.
    """
    # Sort by radius (outer->inner)
    sort_idx = np.argsort(r)[::-1]
    r_sorted = r[sort_idx]
    rho_sorted = rho[sort_idx]
    
    # Calculate shell widths
    dr = np.zeros_like(r_sorted)
    dr[1:] = np.abs(r_sorted[1:] - r_sorted[:-1])
    dr[0] = r_sorted[0] * 0.01  # guess for outermost shell
    
    # Optical depth (integrate from outer boundary inward)
    dtau = kappa * rho_sorted * dr
    tau = np.cumsum(dtau)
    
    # Map back to original order
    tau_original = np.zeros_like(r)
    tau_original[sort_idx] = tau
    
    return tau_original


def downsample_grid(r: np.ndarray, v: np.ndarray, rho: np.ndarray, T: np.ndarray,
                    composition: dict[str, np.ndarray], target_n: int = 200) -> tuple:
    """Downsample grid to target number of points using log spacing in radius.
    
    Returns (r_new, v_new, rho_new, T_new, composition_new)
    """
    n = len(r)
    if n <= target_n:
        return r, v, rho, T, composition
    
    # Log-spaced indices
    log_r = np.log10(r)
    idx = np.linspace(0, n-1, target_n, dtype=int)
    
    r_new = r[idx]
    v_new = v[idx]
    rho_new = rho[idx]
    T_new = T[idx]
    
    composition_new = {}
    for elem, comp in composition.items():
        composition_new[elem] = comp[idx]
    
    return r_new, v_new, rho_new, T_new, composition_new


def downsample_grid_feature_aware(r: np.ndarray, v: np.ndarray, rho: np.ndarray, T: np.ndarray,
                                  composition: dict[str, np.ndarray], target_n: int = 200,
                                  velocity_threshold: float = 0.05, 
                                  density_threshold: float = 0.15,
                                  temperature_threshold: float = 0.10) -> tuple:
    """Downsample grid preserving high-gradient regions (shock fronts, composition boundaries).
    
    Identifies regions with large gradients in velocity, density, or temperature
    and forces inclusion of points from these regions. Fills remaining budget
    with log-spaced points.
    
    Parameters
    ----------
    r : np.ndarray
        Radius array (cm)
    v : np.ndarray  
        Velocity array (km/s)
    rho : np.ndarray
        Density array (g/cm^3)
    T : np.ndarray
        Temperature array (K)
    composition : dict
        Dictionary mapping element names to mass fraction arrays
    target_n : int
        Target number of output grid points
    velocity_threshold : float
        Relative velocity gradient threshold for feature detection
    density_threshold : float
        Relative density gradient threshold for feature detection  
    temperature_threshold : float
        Relative temperature gradient threshold for feature detection
    
    Returns (r_new, v_new, rho_new, T_new, composition_new)
    """
    n = len(r)
    if n <= target_n:
        return r, v, rho, T, composition
    
    # Calculate relative gradients
    log_r = np.log10(r)
    log_v = np.log10(v + 1e-10)  # Add small value to avoid log(0)
    log_rho = np.log10(rho + 1e-20)
    log_T = np.log10(T + 1e-10)
    
    # Compute gradients
    dv_dr = np.abs(np.gradient(log_v, log_r))
    drho_dr = np.abs(np.gradient(log_rho, log_r))
    dT_dr = np.abs(np.gradient(log_T, log_r))
    
    # Normalize gradients
    dv_dr_norm = dv_dr / (np.max(dv_dr) + 1e-10)
    drho_dr_norm = drho_dr / (np.max(drho_dr) + 1e-10)
    dT_dr_norm = dT_dr / (np.max(dT_dr) + 1e-10)
    
    # Identify high-gradient regions
    high_grad_mask = (dv_dr_norm > velocity_threshold) | \
                     (drho_dr_norm > density_threshold) | \
                     (dT_dr_norm > temperature_threshold)
    
    high_grad_indices = np.where(high_grad_mask)[0]
    
    # Always include first and last points
    forced_indices = set([0, n-1])
    
    # Add high-gradient region points
    for idx in high_grad_indices:
        # Add the point and its neighbors to capture the feature
        for offset in [-1, 0, 1]:
            j = idx + offset
            if 0 <= j < n:
                forced_indices.add(j)
    
    # Convert to sorted array
    forced_indices = sorted(forced_indices)
    
    # Calculate remaining budget for log-spaced points
    remaining_budget = target_n - len(forced_indices)
    if remaining_budget < 10:
        remaining_budget = 10  # Ensure some log-spaced points
    
    # Fill remaining budget with log-spaced points, avoiding forced indices
    all_indices = set(range(n))
    available_indices = sorted(all_indices - set(forced_indices))
    
    if len(available_indices) > 0:
        # Log-spaced selection from available indices
        log_available = np.log10(np.array(available_indices) + 1)
        log_idx = np.linspace(0, len(log_available)-1, remaining_budget, dtype=int)
        selected_available = [available_indices[i] for i in log_idx]
    else:
        selected_available = []
    
    # Combine forced and selected indices
    final_indices = sorted(forced_indices + selected_available)
    
    # Ensure we don't exceed target_n
    if len(final_indices) > target_n:
        # Trim from the log-spaced section first
        excess = len(final_indices) - target_n
        final_indices = final_indices[:target_n]
    
    # Apply downsampling
    r_new = r[final_indices]
    v_new = v[final_indices]
    rho_new = rho[final_indices]
    T_new = T[final_indices]
    
    composition_new = {}
    for elem, comp in composition.items():
        composition_new[elem] = comp[final_indices]
    
    print(f"  Feature-aware downsampling: {n} -> {len(final_indices)} zones")
    print(f"    High-gradient regions detected: {len(high_grad_indices)} zones")
    print(f"    Forced points: {len(forced_indices)}, Log-spaced points: {len(selected_available)}")
    
    return r_new, v_new, rho_new, T_new, composition_new


def convert_to_sedona(swd_data: dict, output_file: str, species_list: list[int] | None = None,
                      time_days: float | None = None, homologous: bool = False,
                      kappa: float = 0.34, tau_max: float = 30.0, 
                      tau_min: float = 1e-4, target_n: int = 200,
                      feature_aware_downsample: bool = True,
                      velocity_threshold: float = 0.05,
                      density_threshold: float = 0.15,
                      temperature_threshold: float = 0.10) -> None:
    """Convert STELLA snapshot to SEDONA .mod format.
    
    Parameters
    ----------
    swd_data : dict
        Dictionary with 'r', 'v_kms', 'rho', 'T', 'flux', 'composition' arrays
    output_file : str
        Path to output .mod file
    species_list : list[int] or None
        Atomic numbers to include. Default: H, He, C, O, Si, Fe, Co, Ni
    time_days : float or None
        Time since explosion in days
    homologous : bool
        If True, writes correct time and outer radius for hydro_module="homologous"
    kappa : float
        Electron scattering opacity for optical depth calculation
    tau_max : float
        Outer optical depth cutoff
    tau_min : float
        Inner optical depth cutoff
    target_n : int
        Target number of grid points after downsampling
    feature_aware_downsample : bool
        If True, use feature-aware downsampling to preserve shock fronts
    velocity_threshold : float
        Relative velocity gradient threshold for feature detection
    density_threshold : float
        Relative density gradient threshold for feature detection  
    temperature_threshold : float
        Relative temperature gradient threshold for feature detection
    """
    # Extract data
    r = swd_data['r']              # cm
    v_kms = swd_data['v_kms']      # km/s
    rho = swd_data['rho']          # g/cm^3
    T = swd_data['T']              # K
    flux = swd_data['flux']        # erg/s
    
    # Get and normalize composition from swd_data
    if 'composition' in swd_data and swd_data['composition'] is not None:
        composition = normalize_swd_composition(swd_data['composition'])
    else:
        print(f"WARNING: No composition data in swd_data, using zeros")
        composition = {name: np.zeros(len(r)) for name in STELLA_ELEM_NAMES}
    
    # Calculate optical depth
    tau = calculate_optical_depth(r, rho, kappa)
    
    # Apply optical depth cuts
    mask = (tau <= tau_max) & (tau >= tau_min)
    if not np.any(mask):
        print(f"WARNING: No zones within tau range [{tau_min}, {tau_max}]")
        print(f"  tau range: [{tau.min():.3e}, {tau.max():.3e}]")
        print(f"  r range: [{r.min():.3e}, {r.max():.3e}] cm")
        print(f"  rho range: [{rho.min():.3e}, {rho.max():.3e}] g/cm^3")
        # If tau is extremely high, skip tau filtering and use all zones
        if tau.max() > 1e10:
            print(f"  Tau values are extremely high - skipping tau filtering")
            mask = np.ones_like(tau, dtype=bool)
        else:
            mask = (tau <= tau_max * 10) & (tau >= tau_min * 0.1)  # Expand range
    
    r = r[mask]
    v_kms = v_kms[mask]
    rho = rho[mask]
    T = T[mask]
    tau = tau[mask]
    
    # Filter composition (already in swd_data, just apply mask)
    composition_filtered = {}
    for elem, comp in composition.items():
        if len(comp) == len(mask):
            composition_filtered[elem] = comp[mask]
        else:
            # Interpolate composition to match grid
            from scipy.interpolate import interp1d
            orig_r = np.linspace(0, 1, len(comp))
            new_r = np.linspace(0, 1, np.sum(mask))
            comp_interp = interp1d(orig_r, comp, kind='linear', fill_value='extrapolate')
            composition_filtered[elem] = comp_interp(new_r)
    
    # Sort by radius (inner->outer for SEDONA)
    sort_idx = np.argsort(r)
    r = r[sort_idx]
    v_kms = v_kms[sort_idx]
    rho = rho[sort_idx]
    T = T[sort_idx]
    
    for elem in composition_filtered:
        composition_filtered[elem] = composition_filtered[elem][sort_idx]
    
    # Downsample if needed
    n_current = len(r)
    if n_current > target_n:
        if feature_aware_downsample:
            print(f"Using feature-aware downsampling to preserve shock fronts")
            r, v_kms, rho, T, composition_filtered = downsample_grid_feature_aware(
                r, v_kms, rho, T, composition_filtered, target_n,
                velocity_threshold, density_threshold, temperature_threshold
            )
        else:
            print(f"Using simple log-spaced downsampling")
            r, v_kms, rho, T, composition_filtered = downsample_grid(
                r, v_kms, rho, T, composition_filtered, target_n
            )
    
    # Default species list
    if species_list is None:
        species_list = [1, 2, 6, 8, 14, 26, 27, 28]  # H, He, C, O, Si, Fe, Co, Ni
    
    # Map element names to atomic numbers
    elem_to_z = {
        'H': 1, 'He': 2, 'C': 6, 'N': 7, 'O': 8, 'Ne': 10, 'Na': 11,
        'Mg': 12, 'Al': 13, 'Si': 14, 'S': 16, 'Ar': 18, 'Ca': 20, 'Fe': 26, 'Ni': 28
    }
    
    # Build mass fraction array for each species
    n_species = len(species_list)
    n_zones = len(r)
    species_mass_fracs = np.zeros((n_zones, n_species))
    
    for i, z in enumerate(species_list):
        # Find element with this atomic number
        found = False
        for elem, elem_z in elem_to_z.items():
            if elem_z == z and elem in composition_filtered:
                species_mass_fracs[:, i] = composition_filtered[elem]
                found = True
                break
        if not found:
            species_mass_fracs[:, i] = 0.0
    
    # Remove zero-abundance species
    non_zero_mask = np.any(species_mass_fracs > 0, axis=0)
    if not np.all(non_zero_mask):
        removed_z = [species_list[i] for i in range(len(species_list)) if not non_zero_mask[i]]
        print(f"Removing zero-abundance species: {removed_z}")
        species_list = [species_list[i] for i in range(len(species_list)) if non_zero_mask[i]]
        species_mass_fracs = species_mass_fracs[:, non_zero_mask]
        n_species = len(species_list)
    
    # Normalize mass fractions
    for i in range(n_zones):
        total = np.sum(species_mass_fracs[i, :])
        if total > 0:
            species_mass_fracs[i, :] /= total
    
    # Convert units
    v_cms = v_kms * 1e5  # km/s -> cm/s
    
    # Determine header values
    t_model = time_days if time_days is not None else 0.0
    t_sec = t_model * 24.0 * 3600.0
    r_max = r[-1]
    v_max = np.max(v_cms) if len(v_cms) > 0 else 1e9

    # SEDONA expects species in "Z.A" format (e.g., "1.1 2.4 8.16 14.28 26.56")
    # Map atomic number Z to most common isotope mass number A
    Z_TO_A = {
        1: 1,    # H
        2: 4,    # He
        6: 12,   # C
        7: 14,   # N
        8: 16,   # O
        10: 20,  # Ne
        11: 22,  # Na
        12: 24,  # Mg
        14: 28,  # Si
        16: 32,  # S
        18: 36,  # Ar
        20: 40,  # Ca
        22: 44,  # Ti
        24: 48,  # Cr
        26: 56,  # Fe
        27: 59,  # Co
        28: 58,  # Ni
    }
    species_z_a = [f"{z}.{Z_TO_A.get(z, z*2)}" for z in species_list]

    # Write output file
    with open(output_file, 'w') as f:
        # Header: nzon  rmin  t_now  n_species
        # SEDONA reads: rmin (inner radius), texp (time in seconds), n_elems
        # See grid_1D_sphere.cpp:288-293
        # rmin: inner boundary radius. Must NOT be 0 for a trimmed CSM shell,
        # because SEDONA creates zone 0 from rmin to r[0]. If rmin=0 and r[0]=6.7e14,
        # zone 0 gets volume 4/3*pi*r[0]^3 = 1.28e45 cm^3, creating ~16 Msun of
        # fake mass and tau~1200, trapping all photons.
        # Set rmin to just inside r[0] (half a zone width inward).
        if n_zones > 1:
            dr0 = r[1] - r[0]
        else:
            dr0 = r[0] * 0.01
        rmin = max(r[0] - dr0, 0.0)

        f.write("1D_sphere standard\n")
        if homologous:
            f.write(f"{n_zones}\t{t_sec:.1f}\t{t_sec:.1f}\t{n_species}\n")
            print(f"Header (homologous): {n_zones} zones, t_start=t_now={t_sec:.1f}s ({t_model:.2f}d)")
        else:
            # Static snapshot: rmin=inner boundary, t_now=actual time in seconds
            f.write(f"{n_zones}\t{rmin:.4e}\t{t_sec:.1f}\t{n_species}\n")
            print(f"Header (static): {n_zones} zones, rmin={rmin:.3e} cm, t_now={t_sec:.1f}s ({t_model:.2f}d)")
        f.write(" ".join(species_z_a) + "\n")
        
        # Data lines: r, v, rho, T, mass_fracs...
        for i in range(n_zones):
            line = f"{r[i]:.4e} {v_cms[i]:.4e} {rho[i]:.4e} {T[i]:.4e}"
            for j in range(n_species):
                line += f" {species_mass_fracs[i, j]:.4f}"
            line += "\n"
            f.write(line)
    
    print(f"Written {output_file}: {n_zones} zones, {n_species} species")
    print(f"  Time: {t_model:.3f} days")
    print(f"  Radius range: [{r[0]:.3e}, {r[-1]:.3e}] cm")
    print(f"  Velocity range: [{v_cms[0]:.3e}, {v_cms[-1]:.3e}] cm/s")
    print(f"  Density range: [{rho[0]:.3e}, {rho[-1]:.3e}] g/cm^3")
    print(f"  Temperature range: [{T[0]:.1f}, {T[-1]:.1f}] K")
    print(f"  Optical depth range: [{tau[0]:.3e}, {tau[-1]:.3e}]")


def main():
    parser = argparse.ArgumentParser(
        description="Convert STELLA .swd snapshot to SEDONA .mod format"
    )
    parser.add_argument('swd_file', help='Path to STELLA .swd file')
    parser.add_argument('output_file', help='Path to output .mod file')
    parser.add_argument('--target-days', type=float, default=None,
                        help='Target epoch in days (default: last epoch)')
    parser.add_argument('--kappa', type=float, default=0.34,
                        help='Electron scattering opacity (cm^2/g, default: 0.34)')
    parser.add_argument('--tau-max', type=float, default=30.0,
                        help='Outer optical depth cutoff (default: 30)')
    parser.add_argument('--tau-min', type=float, default=1e-4,
                        help='Inner optical depth cutoff (default: 1e-4)')
    parser.add_argument('--nzon-target', type=int, default=200,
                        help='Target number of grid points (default: 200)')
    parser.add_argument('--time-days', type=float, default=None,
                        help='Time since explosion in days')
    parser.add_argument('--species', type=str, default=None,
                        help='Comma-separated atomic numbers (e.g. 1,2,6,8,14,26)')
    parser.add_argument('--homologous', action='store_true',
                        help='Write correct time and radius for hydro_module="homologous"')
    parser.add_argument('--no-feature-aware-downsample', action='store_true',
                        help='Disable feature-aware downsampling (use simple log-spacing)')
    parser.add_argument('--velocity-threshold', type=float, default=0.05,
                        help='Velocity gradient threshold for feature detection (default: 0.05)')
    parser.add_argument('--density-threshold', type=float, default=0.15,
                        help='Density gradient threshold for feature detection (default: 0.15)')
    parser.add_argument('--temperature-threshold', type=float, default=0.10,
                        help='Temperature gradient threshold for feature detection (default: 0.10)')
    
    args = parser.parse_args()
    
    print(f"Reading {args.swd_file}...")
    epochs = read_swd(args.swd_file)
    
    if not epochs:
        print("ERROR: No epochs found in .swd file")
        sys.exit(1)
    
    # Select epoch
    if args.target_days is not None:
        closest_day = min(epochs.keys(), key=lambda x: abs(x - args.target_days))
        delta = abs(closest_day - args.target_days)
        print(f"Target days: {args.target_days}, selected epoch: {closest_day:.3f} days")
        if delta > 5.0:
            print(f"WARNING: Closest available epoch ({closest_day:.2f} d) is {delta:.2f} days "
                  f"from requested ({args.target_days:.2f} d). "
                  f"STELLA may have crashed before reaching this epoch.")
            available = sorted(epochs.keys())
            print(f"  Available epoch range: {available[0]:.2f} ... {available[-1]:.2f} days")
        elif delta > 1.0:
            print(f"  Note: Using epoch {delta:.2f} days from requested (within 5-day tolerance)")
    else:
        closest_day = max(epochs.keys())
        print(f"Using last epoch: {closest_day:.3f} days")
    
    swd_data = epochs[closest_day]
    
    # Parse species list
    species_list = None
    if args.species:
        species_list = [int(x.strip()) for x in args.species.split(',')]
    
    # Convert
    print(f"Converting to SEDONA format...")
    # Use the selected epoch time if time_days not specified
    time_to_use = args.time_days if args.time_days is not None else closest_day
    convert_to_sedona(
        swd_data, args.output_file,
        species_list=species_list,
        time_days=time_to_use,
        homologous=args.homologous,
        kappa=args.kappa,
        tau_max=args.tau_max,
        tau_min=args.tau_min,
        target_n=args.nzon_target,
        feature_aware_downsample=not args.no_feature_aware_downsample,
        velocity_threshold=args.velocity_threshold,
        density_threshold=args.density_threshold,
        temperature_threshold=args.temperature_threshold
    )


if __name__ == '__main__':
    main()