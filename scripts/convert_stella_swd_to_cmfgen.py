"""Convert a STELLA *.swd snapshot to CMFGEN SN_HYDRO_DATA format.

STELLA *.swd file format (ASCII, written by PRIBAL in lbalsw.trf):
  Each epoch block: one line per zone, columns:
    t_day   km   log10(M_ext_Msun)   log10(R_cm)   V_kms
    log10(T_gas_K)   log10(T_rad_K)   log10(rho_gcc)
    log10(P_cgs)   log10(Qv)   log10(eng)   Flum_1e40   WRKX

  The km=1 (outer boundary) line has t_day as first field;
  all interior lines start with 0 or blank in that column.

  Lines with km=1 or km=Nzon are "header" lines (t_day appears);
  interior lines start with a leading zero or spaces in the t field.

Usage:
  python convert_stella_swd_to_cmfgen.py SWD_FILE ABN_FILE OUTPUT_FILE [OPTIONS]
    --target-days T    Pick the epoch closest to T days (default: last epoch)
    --tau-max X        Outer optical depth cutoff (default: 30)
    --tau-min X        Inner optical depth cutoff (default: 1e-4)
    --nzon-target N    Target number of output grid points (default: 200)
    --time-days T      Explosion time in days (default: read from swd)
    --model-name NAME  For header comment (default: stella)
    --kappa K          Electron-scattering opacity cm^2/g (default: 0.34)
    --xni-mass X       56Ni mass in Msun for luminosity estimate (default: 0)
    --abn-file PATH    *.abn file from STELLA for mass fractions (optional)

The output SN_HYDRO_DATA is in CMFGEN outer->inner (decreasing r) order.
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


def read_swd(swd_path: str) -> dict[float, dict]:
    """Parse a STELLA *.swd file.

    Returns a dict keyed by t_day -> {'r': [], 'v_kms': [], 'rho': [], 'T': [], 'flux': [], 'km': []}.
    All values are numpy arrays sorted outer->inner (decreasing r).
    """
    epochs: dict[float, dict] = {}
    current_t: float | None = None
    buf_km   = []
    buf_r    = []
    buf_v    = []
    buf_T    = []
    buf_rho  = []
    buf_flux = []

    def _flush():
        nonlocal current_t, buf_km, buf_r, buf_v, buf_T, buf_rho, buf_flux
        if current_t is not None and buf_r:
            order = np.argsort(buf_r)[::-1]  # outer->inner
            epochs[current_t] = {
                'km':   np.array(buf_km)[order],
                'r':    np.array(buf_r)[order],
                'v':    np.array(buf_v)[order],
                'T':    np.array(buf_T)[order],
                'rho':  np.array(buf_rho)[order],
                'flux': np.array(buf_flux)[order],
            }
        buf_km   = []
        buf_r    = []
        buf_v    = []
        buf_T    = []
        buf_rho  = []
        buf_flux = []

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
                v_kms       = float(cols[4]) * 1000.0  # .swd unit: 1e8 cm/s -> km/s
                log_T       = float(cols[5])
                log_rho     = float(cols[7])
                flux_1e40   = float(cols[12])  # Flux in units of 1e40 erg/s/cm^2
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
            flux  = flux_1e40 * 1e40  # Convert to cgs: erg/s/cm^2

            buf_km.append(km)
            buf_r.append(r_cm)
            buf_v.append(v_kms)
            buf_T.append(T_K)
            buf_rho.append(rho)
            buf_flux.append(flux)

    _flush()
    if not epochs:
        raise ValueError(f"No epochs found in {swd_path}")
    return epochs


def read_abn(abn_path: str, nzon: int) -> dict[str, np.ndarray]:
    """Read STELLA *.abn file -> mass fractions per zone.

    Columns: km  0  0  0  H  He  C  N  O  Ne  Na  Mg  Al  Si  S  Ar  Ca  FePeak  Ni58  Ni56
    Returns dict of element -> array[nzon], sorted outer->inner (reversed).
    """
    elements = ['H', 'He', 'C', 'N', 'O', 'Ne', 'Na', 'Mg', 'Al', 'Si',
                'S', 'Ar', 'Ca', 'FePeak', 'Ni58', 'Ni56']
    data = {e: [] for e in elements}
    with open(abn_path) as fh:
        for line in fh:
            cols = line.split()
            if len(cols) < 20 or not cols[0].strip().lstrip('-').isdigit():
                continue
            try:
                # Map columns 4-19 to elements (16 elements total)
                # Actual ABN format: cols[4]=H, cols[5]=He, cols[6]=C, cols[7]=N, cols[8]=O, cols[9]=Ne,
                # cols[10]=Na, cols[11]=Mg, cols[12]=Al, cols[13]=Si, cols[14]=S,
                # cols[15]=Ar, cols[16]=Ca, cols[17]=FePeak, cols[18]=Ni58, cols[19]=Ni56
                vals = [float(cols[4]), float(cols[5]), float(cols[6]), float(cols[7]), float(cols[8]),
                        float(cols[9]), float(cols[10]), float(cols[11]), float(cols[12]), float(cols[13]),
                        float(cols[14]), float(cols[15]), float(cols[16]), float(cols[17]), float(cols[18]),
                        float(cols[19])]
            except ValueError:
                continue
            for e, v in zip(elements, vals):
                data[e].append(v)
    for e in elements:
        arr = np.array(data[e])
        if len(arr) < nzon:
            raise ValueError(f"ABN file has {len(arr)} zones, expected {nzon}")
        # If ABN has more zones than SWD, truncate to match SWD
        if len(arr) > nzon:
            arr = arr[:nzon]
        data[e] = arr[::-1]  # inner->outer to outer->inner
    return data


def read_input_profile(input_profile_path: str) -> dict[str, np.ndarray]:
    """Read HERACLES input_profile file -> mass fractions per zone.

    Columns: r_cm rho_cgs p_cgs u_cgs fx_H fx_He fx_O fx_Si fx_Ca fx_Fe fx_Ni
    Returns dict of element -> array[nzon], sorted outer->inner (reversed).
    """
    elements = ['H', 'He', 'O', 'Si', 'Ca', 'Fe', 'Ni']
    data = {e: [] for e in elements}
    with open(input_profile_path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            cols = line.split()
            if len(cols) < 11:
                continue
            try:
                # Map columns to elements
                # cols[4]=fx_H, cols[5]=fx_He, cols[6]=fx_O, cols[7]=fx_Si, cols[8]=fx_Ca, cols[9]=fx_Fe, cols[10]=fx_Ni
                vals = [float(cols[4]), float(cols[5]), float(cols[6]), float(cols[7]),
                        float(cols[8]), float(cols[9]), float(cols[10])]
            except ValueError:
                continue
            for e, v in zip(elements, vals):
                data[e].append(v)
    for e in elements:
        arr = np.array(data[e])
        data[e] = arr[::-1]  # inner->outer to outer->inner
    return data


def read_xni(xni_path: str, nzon: int) -> np.ndarray:
    """Read STELLA *.xni file -> Ni56 mass fraction per zone.
    
    The .xni file is binary with 32-bit floats:
    - 2 header floats (metadata)
    - nzon zone values (Ni56 mass fraction, alternating with placeholders)
    - nzon padding values (duplicate pattern)
    
    Returns array[nzon] of Ni56 mass fractions, sorted outer->inner.
    """
    with open(xni_path, 'rb') as f:
        data = f.read()
    
    # Unpack as 32-bit floats
    n_floats = len(data) // 4
    floats = struct.unpack('f' * n_floats, data)
    
    # Extract Ni56 values (indices 2, 4, 6, ... up to 2+nzon)
    # Based on generate_xni_file.py pattern: even indices have mass fraction
    ni56_values = []
    for i in range(nzon):
        idx = 2 + 2 * i  # Skip 2 header, then every other value
        if idx < len(floats):
            ni56_values.append(floats[idx])
    
    arr = np.array(ni56_values)
    if len(arr) < nzon:
        raise ValueError(f"XNI file has {len(arr)} zones, expected {nzon}")
    if len(arr) > nzon:
        arr = arr[:nzon]
    
    # The .xni file is written inner->outer, but we need outer->inner
    arr = arr[::-1]
    
    return arr


def read_lbol(lbol_path):
    """Read STELLA *.lbol file -> bolometric luminosity vs time.
    
    Returns dict keyed by t_day -> L_bol (Lsun).
    """
    data = np.loadtxt(lbol_path, skiprows=1)
    # Column 0: time (days), Column 2: L_bol (Lsun)
    lbol_dict = {}
    for row in data:
        t_day = row[0]
        L_bol = row[2]
        lbol_dict[t_day] = L_bol
    return lbol_dict


# CMFGEN element name abbreviations
CMFGEN_ELEM_NAMES = {
    'H': 'HYD',
    'He': 'HE',
    'C': 'CARB',
    'N': 'NIT',
    'O': 'OXY',
    'Ne': 'NEON',
    'Na': 'SOD',
    'Mg': 'MAG',
    'Al': 'ALUM',
    'Si': 'SIL',
    'S': 'SUL',
    'Ar': 'ARG',
    'Ca': 'CAL',
    'FePeak': 'IRON',
    'Ni58': 'NI58',
    'Ni56': 'NI56',
}

# Element name mapping: STELLA → HERACLES-style (for opacity calculation)
STELLA_TO_HERACLES_ELEM = {
    'H': 'HYD',
    'He': 'HE',
    'C': 'CARB',
    'N': 'NITR',
    'O': 'OXY',
    'Ne': 'NEON',
    'Na': 'SOD',
    'Mg': 'MAGN',
    'Al': 'ALUM',
    'Si': 'SIL',
    'S': 'SULF',
    'Ar': 'ARG',
    'Ca': 'CAL',
    'FePeak': 'IRON',
    'Ni58': 'NI58',
    'Ni56': 'NI56',
}

# HERACLES element name mapping: HERACLES → STELLA (for converting back)
HERACLES_TO_STELLA_ELEM = {v: k for k, v in STELLA_TO_HERACLES_ELEM.items()}


def calculate_electron_density(density, temperature, mass_fractions):
    """Calculate electron density (cm^-3) accounting for ionization state.
    
    Args:
        density: Density in g/cm³ (array)
        temperature: Temperature in K (array)
        mass_fractions: Dict of mass fractions (e.g., {'H': X_H, 'He': X_He, ...})
    
    Returns:
        n_e: Electron density in cm^-3 (array)
    """
    N_A = 6.02214076e23  # Avogadro's number (mol^-1)
    
    # Atomic data: (atomic mass A, atomic number Z)
    # Map STELLA element names to atomic data
    atomic_data = {
        'H': (1.0,  1),   # Hydrogen
        'He': (4.0,  2),  # Helium
        'C':  (12.0, 6),  # Carbon
        'N':  (14.0, 7),  # Nitrogen
        'O':  (16.0, 8),  # Oxygen
        'Ne': (20.0, 10), # Neon
        'Na': (23.0, 11), # Sodium
        'Mg': (24.0, 12), # Magnesium
        'Al': (27.0, 13), # Aluminum
        'Si': (28.0, 14), # Silicon
        'S':  (32.0, 16), # Sulfur
        'Ar': (40.0, 18), # Argon
        'Ca': (40.0, 20), # Calcium
        'FePeak': (56.0, 26), # Iron group
        'Ni58': (58.0, 28), # Nickel-58
        'Ni56': (56.0, 28), # Nickel-56
    }
    
    n_e = np.zeros_like(density)
    
    for elem, (A, Z) in atomic_data.items():
        if elem in mass_fractions:
            X = mass_fractions[elem]  # Mass fraction
            if np.any(X > 0):
                # Mean ionization state ⟨Z⟩ for this element
                if elem in ['H']:  # Hydrogen
                    # Fit to Saha equation for H I <-> H II
                    # ⟨Z⟩_H ≈ 1 / (1 + exp(157800/T - 10))
                    Z_mean = 1.0 / (1.0 + np.exp(157800.0 / temperature - 10.0))
                elif elem in ['He']:  # Helium
                    # Fit to Saha equation for He I <-> He II and He II <-> He III
                    # ⟨Z⟩_He ≈ 2 / (1 + exp(285000/T - 15))
                    Z_mean = 2.0 / (1.0 + np.exp(285000.0 / temperature - 15.0))
                else:  # Metals (C, O, Ne, Mg, Si, S, Ca, Fe)
                    # Assume fully ionized for T > 10,000 K (typical for SNe)
                    # For lower T, use a smooth transition
                    T_transition = 20000.0  # K (adjust based on element)
                    Z_mean = Z * (1.0 / (1.0 + np.exp(-(temperature - T_transition) / 5000.0)))
                
                # Contribution to electron density: ρ * X * ⟨Z⟩ / A * N_A
                n_e += density * X * Z_mean / A * N_A
    
    return n_e


def approximate_rosseland_opacity(T, rho, X, Y, Z, n_e=None, mass_fractions=None):
    """Approximate Rosseland mean opacity (cm²/g) for supernova ejecta.
    
    Includes electron scattering, bound-free, free-free, and line opacity.
    
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
        Electron density in cm^-3
    mass_fractions : dict, optional
        Dictionary of element mass fractions for better Z/A calculation
    
    Returns:
    --------
    kappa_rosseland : array
        Rosseland mean opacity in cm²/g
    """
    # Avoid division by zero
    T_safe = np.maximum(T, 1.0)  # Minimum 1 K
    rho_safe = np.maximum(rho, 1e-30)  # Avoid zero density
    
    # Electron scattering with composition dependence (more accurate formula)
    kappa_es = 0.2 * (1.0 + X) / (1.0 + 3.0*X + 0.25*Y)
    
    # Calculate Z/A from actual composition if available
    if mass_fractions is not None:
        Z_over_A = 0.0
        Z_sum = 0.0
        A_sum = 0.0
        # Atomic masses (approximate)
        atomic_masses = {
            'H': 1.0, 'He': 4.0, 'C': 12.0, 'N': 14.0, 'O': 16.0,
            'Ne': 20.0, 'Na': 23.0, 'Mg': 24.0, 'Al': 27.0, 'Si': 28.0,
            'S': 32.0, 'Ar': 40.0, 'Ca': 40.0, 'Fe': 56.0, 'Ni': 58.0
        }
        for elem, mass_frac in mass_fractions.items():
            if elem in atomic_masses:
                A = atomic_masses[elem]
                Z = 1.0  # Simplified: assume Z ≈ 1 for all (weighted average)
                if elem == 'He': Z = 2.0
                elif elem in ['C', 'N', 'O']: Z = 6.0 + list(atomic_masses.keys()).index(elem) - 2
                elif elem == 'Ne': Z = 10.0
                elif elem == 'Na': Z = 11.0
                elif elem == 'Mg': Z = 12.0
                elif elem == 'Al': Z = 13.0
                elif elem == 'Si': Z = 14.0
                elif elem == 'S': Z = 16.0
                elif elem == 'Ar': Z = 18.0
                elif elem == 'Ca': Z = 20.0
                elif elem == 'Fe': Z = 26.0
                elif elem == 'Ni': Z = 28.0
                Z_sum += Z * mass_frac
                A_sum += A * mass_frac
        if np.any(A_sum > 0):
            Z_over_A = Z_sum / A_sum
        else:
            Z_over_A = 0.5  # Fallback
    else:
        Z_over_A = 0.5  # Approximate for solar-like metals
    
    # Bound-free (photoionization) using Kramers formula with improved Z/A
    kappa_bf = 4.34e25 * Z_over_A * Z * rho_safe * (T_safe ** -3.5)
    
    # Free-free (bremsstrahlung) using Kramers formula
    kappa_ff = 1.7e20 * (1.0 + X) * rho_safe * (T_safe ** -3.5)
    
    # Note: a kappa_line ~ Z*rho*T^{-2} term was previously included here but removed.
    # It created a hard discontinuity at T=10000K (the mask threshold) that produced
    # unphysical jumps of ~2000x in kappa across adjacent zones with identical composition.
    # CMFGEN computes line opacities self-consistently from atomic data; this crude
    # approximation is not needed and actively misleads the tau/grid calculation.

    # Total Rosseland opacity
    kappa_rosseland = kappa_es + kappa_bf + kappa_ff
    
    return kappa_rosseland


def compute_kappa_total(rho, T_K, mass_fractions):
    """Compute total Rosseland opacity using HERACLES-style calculation.
    
    Uses approximate_rosseland_opacity with composition-dependent
    electron scattering, bound-free, free-free, and line opacity components.
    """
    # Extract mass fractions
    X = mass_fractions.get('H', np.ones(len(rho)) * 0.7)
    Y = mass_fractions.get('He', np.ones(len(rho)) * 0.28)
    
    # Calculate metal mass fraction Z (sum of all metals except H and He)
    Z = np.zeros(len(rho))
    for elem in mass_fractions:
        if elem not in ['H', 'He']:
            Z += mass_fractions[elem]
    
    # Calculate Rosseland opacity with improved composition dependence
    kappa_total = approximate_rosseland_opacity(T_K, rho, X, Y, Z, mass_fractions=mass_fractions)
    
    return kappa_total


def compute_kappa_tau(r_cm, rho, kappa_arr):
    """Compute cumulative optical depth tau[i] = integral from r[0] inward."""
    n = len(r_cm)
    tau = np.zeros(n)
    for i in range(1, n):
        dr = abs(r_cm[i-1] - r_cm[i])
        chi = kappa_arr[i] * 0.5 * (rho[i-1] + rho[i])
        tau[i] = tau[i-1] + chi * dr
    return tau


def luminosity_bb(T_eff_K: float, r_phot_cm: float) -> float:
    return 4.0 * np.pi * r_phot_cm**2 * SIGMA_SB * T_eff_K**4


def sigma_velocity(r, v_kms):
    """Compute dlnV/dlnr - 1 = sigma."""
    n = len(r)
    sigma = np.zeros(n)
    for i in range(1, n - 1):
        # dlnV/dlnr: use central differences
        # For radial grid where r decreases inward (r[0] = outermost):
        # dlnv = ln(v[i+1]) - ln(v[i-1])  (v increases inward)
        # dlnr = ln(r[i+1]) - ln(r[i-1])  (r decreases inward, so dlnr < 0)
        dlnv = np.log(max(abs(v_kms[i+1]), 1e-30)) - np.log(max(abs(v_kms[i-1]), 1e-30))
        dlnr = np.log(max(r[i+1], 1e-30)) - np.log(max(r[i-1], 1e-30))
        if abs(dlnr) > 1e-30:
            sigma[i] = dlnv / dlnr - 1.0
    sigma[0]  = sigma[1]
    sigma[-1] = sigma[-2]
    return sigma


def smooth_velocity_discontinuities(r, v_kms, sigma_threshold=2.0, smooth_width=5):
    """Smooth sharp velocity discontinuities to reduce sigma spikes.
    
    Detects regions where sigma exceeds threshold and applies aggressive moving average smoothing
    to create a smooth velocity transition at ejecta interfaces.
    
    Args:
        r: Radius array (cm), outer->inner
        v_kms: Velocity array (km/s), outer->inner
        sigma_threshold: Sigma value above which smoothing is applied (default: 2.0)
        smooth_width: Number of points on each side for smoothing (default: 5)
    
    Returns:
        Smoothed velocity array
    """
    v_smooth = v_kms.copy()
    sigma = sigma_velocity(r, v_kms)
    
    # Find indices where sigma exceeds threshold
    spike_indices = np.where(np.abs(sigma) > sigma_threshold)[0]
    
    if len(spike_indices) == 0:
        return v_smooth
    
    print(f"[Velocity Smoothing] Found {len(spike_indices)} points with |sigma| > {sigma_threshold}")
    print(f"[Velocity Smoothing] WARNING: This is a physical ejecta interface. Consider increasing grid resolution instead.")
    
    # Group contiguous spike indices into regions
    regions = []
    if len(spike_indices) > 0:
        current_region = [spike_indices[0]]
        for idx in spike_indices[1:]:
            if idx == current_region[-1] + 1:
                current_region.append(idx)
            else:
                regions.append(current_region)
                current_region = [idx]
        regions.append(current_region)
    
    # Apply aggressive moving average smoothing to each region
    for region in regions:
        region_start = max(0, region[0] - smooth_width)
        region_end = min(len(v_kms) - 1, region[-1] + smooth_width)
        
        if region_end - region_start < 4:
            continue
        
        # Apply multiple passes of moving average with larger window
        window_size = 5
        for pass_num in range(3):  # 3 passes
            v_temp = v_smooth.copy()
            for i in range(region_start + window_size, region_end - window_size):
                # 5-point moving average
                v_temp[i] = np.mean(v_smooth[i-window_size:i+window_size+1])
            v_smooth = v_temp
    
    # Recompute sigma to verify improvement
    new_sigma = sigma_velocity(r, v_smooth)
    max_sigma_reduction = np.max(np.abs(sigma)) - np.max(np.abs(new_sigma))
    print(f"[Velocity Smoothing] Max |sigma| reduced from {np.max(np.abs(sigma)):.2f} to {np.max(np.abs(new_sigma)):.2f}")
    
    return v_smooth


def thin_degenerate_grid(radius, *arrays, min_dr_r=0.005):
    """Remove near-degenerate grid points (dr/r < min_dr_r)."""
    n = len(radius)
    keep = np.zeros(n, dtype=bool)
    keep[0] = True
    keep[n-1] = True
    last_kept = 0
    for i in range(1, n - 1):
        dr = abs(radius[last_kept] - radius[i])
        r_ref = max(abs(radius[last_kept]), abs(radius[i]), 1e-30)
        if dr / r_ref >= min_dr_r:
            keep[i] = True
            last_kept = i
    n_removed = int(np.sum(~keep))
    if n_removed > 0:
        print(f"[Grid Thinning] Removed {n_removed} near-degenerate points "
              f"(dr/r < {min_dr_r*100:.1f}%), {n} -> {int(np.sum(keep))} points")
        radius = radius[keep]
        arrays = tuple(a[keep] for a in arrays)
    return (radius,) + arrays


def downsample_grid(radius, *arrays, target_points=100):
    """Downsample grid to approximately target_points using logarithmic spacing.
    
    This implements Dessart's approach: start with a coarse base grid (80-100 points)
    before adaptive refinement.
    """
    n = len(radius)
    if n <= target_points:
        print(f"[Downsample] Grid already at {n} points <= target {target_points}, skipping")
        return (radius,) + arrays
    
    # Use logarithmic spacing in radius for downsampling
    # This ensures better resolution in inner regions where gradients are steep
    log_r = np.log10(radius)
    log_r_min, log_r_max = log_r.min(), log_r.max()
    
    # Create target indices using log spacing
    target_indices = np.unique(np.linspace(0, n-1, target_points, dtype=int))
    
    # Always keep first and last points
    if target_indices[0] != 0:
        target_indices = np.insert(target_indices, 0, 0)
    if target_indices[-1] != n-1:
        target_indices = np.append(target_indices, n-1)
    
    n_removed = n - len(target_indices)
    print(f"[Downsample] Reduced from {n} to {len(target_indices)} points (target: {target_points})")
    
    radius = radius[target_indices]
    arrays = tuple(a[target_indices] for a in arrays)
    
    return (radius,) + arrays


def adaptive_grid_refinement(radius, v, T, rho, n_e, tau, mass_fractions,
                             max_radius_jump=0.05, max_temp_jump=0.10,
                             max_density_jump=0.15, max_ne_jump=0.10,
                             max_logtau_jump=0.10, max_points=500):
    """Apply Dessart's adaptive grid refinement criteria.
    
    Inserts additional grid points where consecutive points exceed the maximum
    allowed jumps in physical quantities.
    
    Criteria from Dessart et al. 2017/2026:
      - ≤5% jump in radius
      - ≤10% jump in temperature
      - ≤10% jump in electron density
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
    
    Returns:
    --------
    Refined arrays (radius, v, T, rho, n_e, tau, mass_fractions)
    """
    def interpolate_arrays(r_list, v_list, T_list, rho_list, n_e_list, tau_list, mf_lists, i, j, fraction):
        """Interpolate between indices i and j at given fraction (0-1)."""
        r_new = r_list[i] + fraction * (r_list[j] - r_list[i])
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
    n_viol_r = n_viol_T = n_viol_rho = n_viol_tau = 0
    max_dr_r = max_dT_T = max_drho_rho = max_dlogtau = 0.0
    for i in range(len(r_list) - 1):
        dr_r = abs(r_list[i+1] - r_list[i]) / max(abs(r_list[i]), 1e-30)
        dT_T = abs(T_list[i+1] - T_list[i]) / max(T_list[i], 1e-30)
        drho_rho = abs(rho_list[i+1] - rho_list[i]) / max(rho_list[i], 1e-30)
        eps = 1e-30
        dlogtau = abs(np.log10(max(tau_list[i+1], eps)) - np.log10(max(tau_list[i], eps)))

        max_dr_r = max(max_dr_r, dr_r)
        max_dT_T = max(max_dT_T, dT_T)
        max_drho_rho = max(max_drho_rho, drho_rho)
        max_dlogtau = max(max_dlogtau, dlogtau)

        if dr_r > max_radius_jump:       n_viol_r += 1
        if dT_T > max_temp_jump:         n_viol_T += 1
        if drho_rho > max_density_jump:  n_viol_rho += 1
        if dlogtau > max_logtau_jump and tau_list[i] > 0.01:
            n_viol_tau += 1

    print(f"[AMR Audit] Base grid ({len(r_list)} points) violation summary:")
    print(f"  dr/r:     {n_viol_r:4d} edges exceed {max_radius_jump}  (max jump = {max_dr_r:.3f})")
    print(f"  dT/T:     {n_viol_T:4d} edges exceed {max_temp_jump}  (max jump = {max_dT_T:.3f})")
    print(f"  drho/rho: {n_viol_rho:4d} edges exceed {max_density_jump}  (max jump = {max_drho_rho:.3f})")
    print(f"  dlogtau:  {n_viol_tau:4d} edges exceed {max_logtau_jump}  (max jump = {max_dlogtau:.3f})")
    total_violations = n_viol_r + n_viol_T + n_viol_rho + n_viol_tau
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
            eps = 1e-30
            dlogtau = abs(np.log10(max(tau_list[i+1], eps)) - np.log10(max(tau_list[i], eps)))

            # Check if any criterion is violated
            # NOTE: dne/ne is excluded because n_e is computed from an approximate
            # ionization formula that creates artificial jumps at ionization fronts.
            # The physically meaningful quantities (r, T, rho, tau) are checked instead.
            violation = False
            if dr_r > max_radius_jump:
                violation = True
            if dT_T > max_temp_jump:
                violation = True
            if drho_rho > max_density_jump:
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


def write_sn_hydro_data(outpath: str, r_cm, v_kms, rho, T_K, sigma, kappa,
                        time_days: float, mass_fractions: dict,
                        lum_lsun: float, t_floor: float = 3000.0,
                        r_inner_cm: float = None, T_inner_K: float = None,
                        ni56_initial: np.ndarray = None):
    """Write CMFGEN SN_HYDRO_DATA file (outer->inner order)."""
    n = len(r_cm)
    # nfrac is computed after isotope handling below (mass_fractions may be modified)
    r_1e10 = r_cm / 1e10
    
    # Ensure positive time for CMFGEN/SEDONA compatibility
    # Some hydro simulations may have negative times for pre-explosion phases
    if time_days < 0:
        print(f"[Time Correction] Negative time detected: {time_days:.2f} days")
        print(f"  Converting to positive time since explosion bomb introduction")
        time_days = abs(time_days)
    
    # Apply temperature floor
    n_floored = int(np.sum(T_K < t_floor))
    if n_floored > 0:
        print(f"[T floor] Flooring {n_floored} depths from {T_K.min():.0f} K to T_FLOOR={t_floor:.0f} K")
        T_K = np.maximum(T_K, t_floor)

    # Break flat isothermal regions after floor to prevent SEDONA/CMFGEN singular-matrix crashes
    T_K = break_isothermal_regions(T_K, t_floor=t_floor)

    # Recalculate kappa consistently from the floored/capped temperature so that
    # kappa and T stored in SN_HYDRO_DATA are self-consistent.  The kappa passed
    # in was computed before the T floor, so cold zones had kappa_bf ~ T^{-3.5}
    # blown up by orders of magnitude relative to the floored temperature.
    kappa = compute_kappa_total(rho, T_K, mass_fractions)
    print(f"[kappa] Recalculated after T floor: min={kappa.min():.3e}, max={kappa.max():.3e} cm²/g")

    T_10k  = T_K / 1e4

    # Mean atomic mass and electron density (approximate, assuming full ionisation)
    elem_A = {'H': 1.0, 'He': 4.0, 'C': 12.0, 'N': 14.0, 'O': 16.0,
              'Ne': 20.0, 'Na': 23.0, 'Mg': 24.0, 'Al': 27.0, 'Si': 28.0,
              'S': 32.0, 'Ar': 40.0, 'Ca': 40.0, 'FePeak': 56.0, 'Ni58': 58.0, 'Ni56': 56.0}
    elem_Z = {'H': 1.0, 'He': 2.0, 'C': 6.0, 'N': 7.0, 'O': 8.0,
              'Ne': 10.0, 'Na': 11.0, 'Mg': 12.0, 'Al': 13.0, 'Si': 14.0,
              'S': 16.0, 'Ar': 18.0, 'Ca': 20.0, 'FePeak': 26.0, 'Ni58': 28.0, 'Ni56': 28.0}

    inv_A = np.zeros(n)
    Z_over_A = np.zeros(n)
    for elem, frac in mass_fractions.items():
        A = elem_A.get(elem, 28.0)
        Z = elem_Z.get(elem, 14.0)
        inv_A    += frac / A
        Z_over_A += frac * Z / A

    mean_A = np.where(inv_A > 1e-30, 1.0 / inv_A, 28.0)
    atom_density  = rho * AVOGADRO * inv_A          # atoms/cm^3
    electron_density = rho * AVOGADRO * Z_over_A    # e^-/cm^3

    def write_block(f, label, values, per_line=8):
        f.write(f"{label}\n")
        for i in range(0, len(values), per_line):
            chunk = values[i:i+per_line]
            f.write(''.join(f"  {v:22.15E}" for v in chunk) + '\n')
        f.write('\n')

    # Calculate radioactive decay chain: Ni56 -> Co56 -> Fe56
    # Decay timescales (from STELLA source: volenpumnoint.f line 99)
    t_Ni = 8.8      # days (e-folding time for Ni56)
    t_Co = 111.3    # days (e-folding time for Co56)

    # Get initial Ni56 mass fraction from parameter (from .xni file)
    # If not provided, check mass_fractions (backward compatibility)
    if ni56_initial is not None:
        Ni56_initial = ni56_initial
    else:
        Ni56_initial = mass_fractions.get('Ni56', np.zeros(n))

    # Calculate decayed fractions at snapshot time
    Ni56_t = Ni56_initial * np.exp(-time_days / t_Ni)
    Co56_t = Ni56_initial * (t_Co / (t_Co - t_Ni)) * (np.exp(-time_days / t_Co) - np.exp(-time_days / t_Ni))
    Fe56_t = Ni56_initial - Ni56_t - Co56_t  # Total decayed to Fe56

    # Print decay information
    if np.any(Ni56_initial > 1e-10):
        max_Ni56_init = np.max(Ni56_initial)
        max_Ni56_t = np.max(Ni56_t)
        max_Co56_t = np.max(Co56_t)
        max_Fe56_t = np.max(Fe56_t)
        print(f"[Radioactive decay] t={time_days:.1f} days:")
        print(f"  Initial Ni56 max fraction: {max_Ni56_init:.4e}")
        print(f"  Current Ni56 max fraction: {max_Ni56_t:.4e} ({100*max_Ni56_t/max_Ni56_init:.1f}% remaining)")
        print(f"  Co56 max fraction: {max_Co56_t:.4e}")
        print(f"  Fe56 (from decay) max fraction: {max_Fe56_t:.4e}")

    # Update elemental mass fractions and prepare isotopic data for CMFGEN.
    #
    # The Ni56 decay chain is: Ni56 -> Co56 -> Fe56
    #
    # From rd_sn_data.f (ELSE branch, lines ~448-474):
    #   For every element that has isotopes listed in SN_HYDRO_DATA, CMFGEN:
    #     1. Sums the isotope mass fractions for that element
    #     2. Checks that elemental mass fraction is within 20% of the isotope sum (STOP if not)
    #     3. Replaces the elemental population with the isotope sum
    #
    # Therefore, for EACH species that has an isotope entry, the elemental mass
    # fraction must be present AND must equal the sum of its isotopes.
    #
    #   NICK elemental  = Ni56_t    (all Ni56 is radioactive; no stable Ni56 in progenitor)
    #   NICK 56 isotope = Ni56_t
    #
    #   COB elemental   = Co56_t    (all Co56 comes from Ni56 decay; no stable Co56 in progenitor)
    #   COB 56 isotope  = Co56_t
    #
    #   IRON elemental  = stable_fe + Fe56_t   (progenitor Fe + Fe56 from decay)
    #   IRON 56 isotope = stable_fe + Fe56_t   (must equal elemental)
    #
    # The original Ni56 key in mass_fractions is the initial value (pre-decay), so we
    # replace it with the current value. FePeak gets updated to include decayed Fe56.

    # Replace initial Ni56 with current (decayed) value
    if 'Ni56' in mass_fractions:
        del mass_fractions['Ni56']
    mass_fractions['NICK'] = Ni56_t

    # Add Co56 as a new elemental species (none existed before in .abn)
    mass_fractions['COB'] = Co56_t

    # Update FePeak to include Fe56 produced by decay
    stable_fe = np.zeros(n)
    if 'FePeak' in mass_fractions:
        stable_fe = mass_fractions.pop('FePeak')
    total_fe = stable_fe + Fe56_t
    mass_fractions['IRON'] = total_fe

    # Build isotopic data: (species_name, baryon_number) -> mass_fraction array
    # All three isotopes must be present and match their elemental counterparts.
    isotopic_data = {
        ('NICK', 56): Ni56_t,    # matches NICK elemental
        ('COB',  56): Co56_t,    # matches COB  elemental
        ('IRON', 56): total_fe,  # matches IRON elemental
    }

    n_isotopes = len(isotopic_data)
    
    # Map CMFGEN species names to our internal names
    cmfgen_to_internal = {
        'HYD': 'H',
        'HE': 'He', 
        'OXY': 'O',
        'SIL': 'Si',        # CMFGEN uses 'SIL' not 'SILI'
        'CAL': 'Ca',        # CMFGEN uses 'CAL' not 'CALC'
        'IRON': 'Fe',
        'NICK': 'Ni',
        'COB': 'Co'
    }
    
    # Add missing elements with solar abundances if not present
    # This is needed when reusing STELLA output that doesn't have all elements
    solar_abundances = {
        'SIL': 0.0007,  # Silicon (CMFGEN uses 'SIL' not 'SILI')
        'CAL': 0.0002,  # Calcium (CMFGEN uses 'CAL' not 'CALC')
    }
    
    print(f"  Current elements in mass_fractions: {list(mass_fractions.keys())}")
    
    for cmfgen_name, solar_frac in solar_abundances.items():
        internal_name = cmfgen_to_internal.get(cmfgen_name, cmfgen_name)
        if internal_name not in mass_fractions:
            mass_fractions[internal_name] = np.full(n, solar_frac)
            print(f"  Added {cmfgen_name} with solar abundance {solar_frac}")
        else:
            print(f"  {internal_name} already present, skipping solar fallback")
    
    # Filter to only include the 5 elements we want: H, He, O, Si, Ca, Fe + radioactive Ni, Co
    # All other elements (C, N, Ne, Mg, S, etc.) are excluded
    # Include both internal names and CMFGEN names in the allowed list
    # CMFGEN species names: HYD, HE, OXY, SIL, CAL, IRON, NICK, COB
    allowed_elements = {'H', 'He', 'O', 'Si', 'Ca', 'Fe', 'Ni', 'Co',
                       'HYD', 'HE', 'OXY', 'SIL', 'CAL', 'IRON', 'NICK', 'COB'}
    
    # Convert internal names to CMFGEN names for output
    cmfgen_mass_fractions = {}
    for internal_name, frac_array in mass_fractions.items():
        # Only keep allowed elements
        if internal_name not in allowed_elements:
            print(f"  Skipping {internal_name} (not in allowed element list)")
            continue
        
        # Find corresponding CMFGEN name
        cmfgen_name = None
        for cmfgen_name_key, internal in cmfgen_to_internal.items():
            if internal == internal_name:
                cmfgen_name = cmfgen_name_key
                break
        if cmfgen_name is None:
            # Use the internal name as-is if no mapping found
            cmfgen_name = internal_name.upper()
        cmfgen_mass_fractions[cmfgen_name] = frac_array
    
    # Replace mass_fractions with CMFGEN-named version
    mass_fractions = cmfgen_mass_fractions
    
    nfrac = len(mass_fractions)  # Count after filtering

    # Build plot-only mass fractions for visualisation
    mass_fractions_plot = mass_fractions.copy()  # already has NICK, COB, IRON

    with open(outpath, 'w') as f:
        f.write(f"Number of data points:            {n:8d}\n")
        f.write(f"Number of mass fractions:         {nfrac:8d}\n")
        f.write(f"Number of isotopes:               {n_isotopes:8d}\n")
        f.write(f"Time(days) since explosion:      {time_days:15.5E}\n")
        f.write(f"Energy Lum [lsun]:               {lum_lsun:15.6E}\n")
        if r_inner_cm is not None:
            f.write(f"Inner boundary radius [cm]:      {r_inner_cm:15.6E}\n")
        if T_inner_K is not None:
            f.write(f"Inner boundary temperature [K]:  {T_inner_K:15.6E}\n")
        f.write('\n')
        write_block(f, "Radius grid (10^10cm)", r_1e10)
        write_block(f, "Velocity (km/s)",       v_kms)
        write_block(f, "Sigma (dlnV/dlnr-1)",   sigma)
        write_block(f, "Temperature (10^4 K)",  T_10k)
        write_block(f, "Density (gm/cm^3)",     rho)
        write_block(f, "Electron density (number density)", electron_density)
        write_block(f, "Atom density (number density)",     atom_density)
        write_block(f, "Kappa (cm^2/gm)",       kappa)
        for elem, frac in mass_fractions.items():
            # mass_fractions now already has CMFGEN names from our conversion
            write_block(f, f"{elem} mass fraction", frac)

        # Write isotopic data
        # Format: "SPECIES_NAME BARYON_NUMBER mass fraction"
        per_line = 8
        for (iso_name, baryon_num), iso_frac in isotopic_data.items():
            f.write(f"{iso_name} {baryon_num} mass fraction\n")
            for i in range(0, len(iso_frac), per_line):
                chunk = iso_frac[i:i+per_line]
                f.write(''.join(f"  {v:22.15E}" for v in chunk) + '\n')
            f.write('\n')

    print(f"Written {outpath}  ({n} depth points, {nfrac} mass fractions)")
    print(f"  t={time_days:.3f} days, L={lum_lsun:.3e} Lsun")
    print(f"  r=[{r_1e10[0]:.3e},{r_1e10[-1]:.3e}] 1e10cm")
    print(f"  v=[{v_kms[0]:.1f},{v_kms[-1]:.1f}] km/s")


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


def plot_comparison(orig_r, orig_rho, orig_T, orig_v, orig_tau, orig_sigma,
                   new_r, new_rho, new_T, new_v, new_tau, new_sigma,
                   output_file, time_days, tau_min=1e-5, tau_max=30.0,
                   orig_mass_fractions=None, new_mass_fractions=None):
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
    ax.loglog(orig_r_rsun, orig_rho, 'b--', alpha=0.7, linewidth=2, label='Original STELLA')
    ax.loglog(new_r_rsun, new_rho, 'r-', alpha=0.7, linewidth=2, label='CMFGEN')
    ax.set_ylabel('Density (g/cm³)', fontsize=12)
    ax.set_title('Density Profile', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_xlim(r_min, r_max)
    
    # Plot 2: Temperature
    ax = axes[0, 1]
    ax.loglog(orig_r_rsun, orig_T, 'b--', alpha=0.7, linewidth=2, label='Original STELLA')
    ax.loglog(new_r_rsun, new_T, 'r-', alpha=0.7, linewidth=2, label='CMFGEN')
    ax.set_ylabel('Temperature (K)', fontsize=12)
    ax.set_title('Temperature Profile', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_xlim(r_min, r_max)
    
    # Plot 3: Velocity
    ax = axes[0, 2]
    ax.loglog(orig_r_rsun, orig_v, 'b--', alpha=0.7, linewidth=2, label='Original STELLA')
    ax.loglog(new_r_rsun, new_v, 'r-', alpha=0.7, linewidth=2, label='CMFGEN')
    ax.set_ylabel('Velocity (km/s)', fontsize=12)
    ax.set_title('Velocity Profile', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_xlim(r_min, r_max)
    
    # Plot 4: Sigma
    ax = axes[1, 0]
    ax.step(orig_r_rsun, orig_sigma, 'b--', alpha=0.7, linewidth=2, label='Original STELLA', where='post')
    ax.step(new_r_rsun, new_sigma, 'r-', alpha=0.7, linewidth=2, label='CMFGEN', where='post')
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
    ax.step(orig_r_rsun, orig_tau, 'b--', alpha=0.7, linewidth=2, label='Original STELLA', where='post')
    ax.step(new_r_rsun, new_tau, 'r-', alpha=0.7, linewidth=2, label='CMFGEN', where='post')
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
    
    # Plot 6: Mass fractions (if available)
    ax = axes[1, 2]
    if new_mass_fractions and len(new_mass_fractions) > 0:
        # Determine selected regions based on tau
        selected_mask_new = (new_tau >= tau_min) & (new_tau <= tau_max)
        selected_mask_orig = (orig_tau >= tau_min) & (orig_tau <= tau_max)
        
        for elem, color in zip(['H', 'He', 'C', 'O', 'Si', 'Fe'], 
                               ['#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4', '#FFEAA7', '#DDA0DD']):
            if elem in new_mass_fractions:
                # Plot selected region (solid lines)
                if np.any(selected_mask_new):
                    ax.plot(new_r_rsun[selected_mask_new], new_mass_fractions[elem][selected_mask_new], 
                           label=elem, color=color, linewidth=2, linestyle='-')
                # Plot unselected region (dashed lines)
                if np.any(~selected_mask_new):
                    ax.plot(new_r_rsun[~selected_mask_new], new_mass_fractions[elem][~selected_mask_new], 
                           color=color, linewidth=1.5, linestyle='--', alpha=0.6)
        ax.set_xlabel('r/R☉', fontsize=12)
        ax.set_ylabel('Mass Fraction', fontsize=12)
        ax.set_title('Mass Fractions (CMFGEN)', fontsize=14)
        ax.legend(fontsize=9, ncol=2)
        ax.grid(True, alpha=0.3, which='both')
        ax.set_ylim(0, 1.1)
        ax.set_xlim(r_min, r_max)
        ax.set_xscale('log')
    else:
        ax.text(0.5, 0.5, 'No mass fraction data', ha='center', va='center', transform=ax.transAxes)
        ax.set_title('Mass Fractions', fontsize=14)
        ax.set_xlabel('r/R☉', fontsize=12)
        ax.set_xlim(r_min, r_max)
        ax.set_xscale('log')
    
    plt.tight_layout()
    plot_file = output_file.replace('.dat', '_comparison.png').replace('SN_HYDRO_DATA', 'SN_HYDRO_DATA_comparison')
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    print(f"Plot saved to {plot_file}")
    plt.close()


def plot_adaptive_grid_diagnostics(orig_r, orig_v, orig_T, orig_rho, orig_tau,
                                   downsampled_r, downsampled_v, downsampled_T, downsampled_rho, downsampled_tau,
                                   refined_r, refined_v, refined_T, refined_rho, refined_tau,
                                   output_file, time_days,
                                   max_radius_jump=0.05, max_temp_jump=0.10,
                                   max_density_jump=0.15, max_ne_jump=0.10,
                                   max_logtau_jump=0.10, max_points=500):
    """Create diagnostic plots for adaptive grid refinement process.

    Upper row: Density, Temperature, Velocity profiles (with grid markers)
    Lower row: Histograms of jumps before/after AMR with threshold lines
    """
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 3, figsize=(20, 12))
    fig.suptitle(f'Adaptive Grid Diagnostics - {time_days:.1f} days', fontsize=16)

    RSUN = 6.96e10

    # Convert to Rsun for plotting
    orig_r_rsun = orig_r / RSUN
    downsampled_r_rsun = downsampled_r / RSUN
    refined_r_rsun = refined_r / RSUN

    # --- UPPER ROW: Physical profiles with grid markers ---

    n_orig = len(orig_r)
    n_base = len(downsampled_r)
    n_ref = len(refined_r)

    # Plot 1: Density profile
    ax = axes[0, 0]
    ax.plot(orig_r_rsun, orig_rho, '-', linewidth=1.5, label=f'Tau-trimmed ({n_orig})', alpha=0.5, color='gray')
    ax.plot(downsampled_r_rsun, downsampled_rho, 'o', markersize=4, linestyle='None',
            label=f'Base grid ({n_base})', alpha=0.7, color='blue')
    ax.plot(refined_r_rsun, refined_rho, '.', markersize=2, linestyle='None',
            label=f'AMR refined ({n_ref})', alpha=0.9, color='red')
    ax.set_xlabel('r / R☉', fontsize=12)
    ax.set_ylabel('Density (g/cm³)', fontsize=12)
    ax.set_title('Density Profile', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xscale('log')
    ax.set_yscale('log')

    # Plot 2: Temperature profile
    ax = axes[0, 1]
    ax.plot(orig_r_rsun, orig_T, '-', linewidth=1.5, label=f'Tau-trimmed ({n_orig})', alpha=0.5, color='gray')
    ax.plot(downsampled_r_rsun, downsampled_T, 'o', markersize=4, linestyle='None',
            label=f'Base grid ({n_base})', alpha=0.7, color='blue')
    ax.plot(refined_r_rsun, refined_T, '.', markersize=2, linestyle='None',
            label=f'AMR refined ({n_ref})', alpha=0.9, color='red')
    ax.set_xlabel('r / R☉', fontsize=12)
    ax.set_ylabel('Temperature (K)', fontsize=12)
    ax.set_title('Temperature Profile', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xscale('log')
    ax.set_yscale('log')

    # Plot 3: Velocity profile
    ax = axes[0, 2]
    ax.plot(orig_r_rsun, orig_v / 1e5, '-', linewidth=1.5, label=f'Tau-trimmed ({n_orig})', alpha=0.5, color='gray')
    ax.plot(downsampled_r_rsun, downsampled_v / 1e5, 'o', markersize=4, linestyle='None',
            label=f'Base grid ({n_base})', alpha=0.7, color='blue')
    ax.plot(refined_r_rsun, refined_v / 1e5, '.', markersize=2, linestyle='None',
            label=f'AMR refined ({n_ref})', alpha=0.9, color='red')
    ax.set_xlabel('r / R☉', fontsize=12)
    ax.set_ylabel('Velocity (km/s)', fontsize=12)
    ax.set_title('Velocity Profile', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xscale('log')

    # Calculate jumps for base and refined grids
    def calculate_jumps(r, T, rho, tau):
        eps = 1e-30
        dr_r = np.abs(np.diff(r) / r[:-1])
        dT_T = np.abs(np.diff(T) / T[:-1])
        drho_rho = np.abs(np.diff(rho) / rho[:-1])
        dlogtau = np.abs(np.diff(np.log10(np.maximum(tau, eps))))
        return dr_r, dT_T, drho_rho, dlogtau

    dr_r_base, dT_T_base, drho_rho_base, dlogtau_base = calculate_jumps(
        downsampled_r, downsampled_T, downsampled_rho, downsampled_tau)
    dr_r_ref, dT_T_ref, drho_rho_ref, dlogtau_ref = calculate_jumps(
        refined_r, refined_T, refined_rho, refined_tau)

    # --- LOWER ROW: Jump histograms before/after AMR ---

    # Plot 4: Density jump histogram (linear x-axis)
    ax = axes[1, 0]
    bins = np.linspace(0, 1, 40)
    ax.hist(drho_rho_base, bins=bins, alpha=0.5, label=f'Base grid ({n_base}) max={drho_rho_base.max():.3f}', color='blue')
    ax.hist(drho_rho_ref, bins=bins, alpha=0.5, label=f'AMR refined ({n_ref}) max={drho_rho_ref.max():.3f}', color='red')
    ax.axvline(max_density_jump, color='k', linestyle='--', linewidth=2, label=f'Threshold ({max_density_jump})')
    ax.set_xlabel('drho / rho', fontsize=12)
    ax.set_ylabel('Count', fontsize=12)
    ax.set_title('Density Jump Distribution', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(left=0)
    ax.set_yscale('log')

    # Plot 5: Temperature jump histogram (linear x-axis)
    ax = axes[1, 1]
    ax.hist(dT_T_base, bins=bins, alpha=0.5, label=f'Base grid ({n_base}) max={dT_T_base.max():.3f}', color='blue')
    ax.hist(dT_T_ref, bins=bins, alpha=0.5, label=f'AMR refined ({n_ref}) max={dT_T_ref.max():.3f}', color='red')
    ax.axvline(max_temp_jump, color='k', linestyle='--', linewidth=2, label=f'Threshold ({max_temp_jump})')
    ax.set_xlabel('dT / T', fontsize=12)
    ax.set_ylabel('Count', fontsize=12)
    ax.set_title('Temperature Jump Distribution', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(left=0)
    ax.set_yscale('log')

    # Plot 6: Velocity (dr/r) jump histogram (linear x-axis)
    ax = axes[1, 2]
    ax.hist(dr_r_base, bins=bins, alpha=0.5, label=f'Base grid ({n_base}) max={dr_r_base.max():.3f}', color='blue')
    ax.hist(dr_r_ref, bins=bins, alpha=0.5, label=f'AMR refined ({n_ref}) max={dr_r_ref.max():.3f}', color='red')
    ax.axvline(max_radius_jump, color='k', linestyle='--', linewidth=2, label=f'Threshold ({max_radius_jump})')
    ax.set_xlabel('dr / r', fontsize=12)
    ax.set_ylabel('Count', fontsize=12)
    ax.set_title('Radius Jump Distribution', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(left=0)
    ax.set_yscale('log')

    plt.tight_layout()
    plot_file = output_file.replace('.dat', '_adaptive_grid.png').replace('SN_HYDRO_DATA', 'SN_HYDRO_DATA_adaptive_grid')
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    print(f"Adaptive grid diagnostic plot saved to {plot_file}")
    plt.close()


def plot_mass_fractions_stacked(radius_cm, density, mass_fractions, time_days, output_file,
                                trimmed_radius_cm=None, trimmed_mass_fractions=None,
                                tau=None, tau_min=1e-5, tau_max=30.0):
    """Create a mass fraction plot vs radius.
    
    Plots original STELLA data as dashed lines and trimmed CMFGEN data as solid lines.
    Isotopes (NICK, COB, IRON) are shown with distinct dotted/dash-dot styles.
    """
    rsun_cm = 6.96e10
    
    # Original data (STELLA - full model)
    r_orig = np.asarray(radius_cm)
    
    # Trimmed data (CMFGEN - selected tau range)
    if trimmed_radius_cm is not None and trimmed_mass_fractions is not None:
        r_trim = np.asarray(trimmed_radius_cm)
        mf_trim = trimmed_mass_fractions
        has_trimmed = True
    else:
        has_trimmed = False
    
    # Elemental species: (key_in_orig, key_in_trim, display_label, color)
    ELEM_SPECS = [
        ('H',      'H',    'H',    '#FF6B6B'),
        ('He',     'He',   'He',   '#4ECDC4'),
        ('C',      'C',    'C',    '#45B7D1'),
        ('O',      'O',    'O',    '#96CEB4'),
        ('Si',     'Si',   'Si',   '#FFEAA7'),
        ('FePeak', 'IRON', 'Fe',   '#DDA0DD'),
    ]
    # Isotope species shown only in trimmed (CMFGEN) data (radioactive parents only)
    ISO_SPECS = [
        ('NICK', 'Ni56 (isotope)', '#e377c2', '-.'),
        ('COB',  'Co56 (isotope)', '#ff7f0e', '-.'),
    ]
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Plot original STELLA data as dashed background
    x_orig = r_orig / rsun_cm
    for key_orig, key_trim, label, color in ELEM_SPECS:
        if key_orig in mass_fractions and np.any(mass_fractions[key_orig] > 1e-30):
            ax.plot(x_orig, mass_fractions[key_orig], color=color,
                    linewidth=1.5, linestyle='--', alpha=0.5)
    
    # Plot trimmed CMFGEN elemental data as solid foreground
    if has_trimmed:
        x_trim = r_trim / rsun_cm
        for key_orig, key_trim, label, color in ELEM_SPECS:
            if key_trim in mf_trim and np.any(mf_trim[key_trim] > 1e-30):
                ax.plot(x_trim, mf_trim[key_trim], color=color,
                        linewidth=2.5, linestyle='-', label=label)
        # Plot isotopes as dash-dot
        for key, label, color, ls in ISO_SPECS:
            if key in mf_trim and np.any(mf_trim[key] > 1e-30):
                ax.plot(x_trim, mf_trim[key], color=color,
                        linewidth=2, linestyle=ls, label=label)
    
    ax.set_xlabel(r'r / R$_\odot$', fontsize=12)
    ax.set_ylabel('Mass Fraction', fontsize=12)
    ax.set_title(f'Mass Fractions: STELLA (dashed) vs CMFGEN (solid) + isotopes (dash-dot)\nt = {time_days:.1f} days', fontsize=14)
    ax.set_xlim(x_orig.min() * 0.9, x_orig.max() * 1.1)
    ax.set_ylim(0, 1.05)
    ax.set_xscale('log')
    ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), fontsize=10)
    ax.grid(True, alpha=0.3, which='both')
    
    plt.tight_layout()
    plot_file = output_file.replace('.dat', '_mass_fraction.png').replace('SN_HYDRO_DATA', 'SN_HYDRO_DATA_mass_fraction')
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    print(f"Mass fraction plot saved to {plot_file}")
    plt.close()


def plot_dr_vs_r(r_cm, output_file, time_days, rho=None, T=None, v=None):
    """Plot grid spacing (|dr|/r) vs radius to visualize grid resolution.
    
    Features:
    - Log scale x-axis for radius
    - Upper x-axis showing grid index
    - Right y-axes showing actual physical quantities (density, temperature, velocity)
    """
    dr = np.diff(r_cm)
    r_mid = 0.5 * (r_cm[1:] + r_cm[:-1])
    
    # Calculate relative spacing |dr|/r (absolute value since grid is outer->inner)
    dr_over_r = np.abs(dr) / r_mid
    
    # Grid indices (0 = outermost, n-1 = innermost)
    grid_indices = np.arange(len(r_mid))
    
    fig, ax = plt.subplots(figsize=(14, 7))
    
    # Main plot: dr/r vs radius
    ax.plot(r_mid / 1e10, dr_over_r, 'b-', linewidth=1.5, label='Grid spacing |dr|/r')
    
    ax.set_xlabel('Radius (10$^{10}$ cm)', fontsize=12)
    ax.set_ylabel('Grid spacing |dr|/r', fontsize=12, color='b')
    ax.tick_params(axis='y', labelcolor='b')
    ax.set_title(f'Grid Resolution: |dr|/r vs r\nt = {time_days:.1f} days', fontsize=14)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid(True, alpha=0.3, which='both')
    
    # Add reference line for 0.5% threshold
    ax.axhline(0.005, color='r', linestyle='--', alpha=0.7, label='0.5% threshold')
    
    # Upper x-axis: Grid Index (linear progression mapped to log radius)
    ax_top = ax.twiny()
    ax_top.set_xlim(ax.get_xlim())
    
    # Create tick locations: evenly spaced in grid index
    # r_mid is ordered outer->inner, so index 0 is outermost (largest r)
    # grid index 0 -> r_max, grid index n-1 -> r_min
    n_ticks = 7
    tick_indices = np.linspace(0, len(r_mid)-1, n_ticks, dtype=int)
    
    # Get corresponding radii for these grid indices
    tick_radii = r_mid[tick_indices] / 1e10
    
    ax_top.set_xticks(tick_radii)
    ax_top.set_xticklabels([f'{idx}' for idx in tick_indices])
    ax_top.set_xlabel('Grid Index (0=outer, n=inner) [linear]', fontsize=12)
    
    # Right y-axes: Physical quantities (actual values, not normalized)
    if rho is not None or T is not None or v is not None:
        # Interpolate to r_mid points
        rho_mid = 0.5 * (rho[1:] + rho[:-1]) if rho is not None else None
        T_mid = 0.5 * (T[1:] + T[:-1]) if T is not None else None
        v_mid = 0.5 * (v[1:] + v[:-1]) if v is not None else None
        
        offset = 0
        if rho is not None:
            ax_rho = ax.twinx()
            ax_rho.semilogy(r_mid / 1e10, rho_mid, 'g--', linewidth=1.0, alpha=0.6, label='Density')
            ax_rho.set_ylabel('Density (g/cm$^3$)', fontsize=11, color='g')
            ax_rho.tick_params(axis='y', labelcolor='g')
            offset += 60
        
        if T is not None:
            ax_T = ax.twinx()
            ax_T.spines['right'].set_position(('outward', offset))
            ax_T.semilogy(r_mid / 1e10, T_mid, 'orange', linestyle='--', linewidth=1.0, alpha=0.6, label='Temperature')
            ax_T.set_ylabel('Temperature (K)', fontsize=11, color='orange')
            ax_T.tick_params(axis='y', labelcolor='orange')
            offset += 60
        
        if v is not None:
            ax_v = ax.twinx()
            ax_v.spines['right'].set_position(('outward', offset))
            ax_v.plot(r_mid / 1e10, v_mid, 'm--', linewidth=1.0, alpha=0.6, label='Velocity')
            ax_v.set_ylabel('Velocity (km/s)', fontsize=11, color='m')
            ax_v.tick_params(axis='y', labelcolor='m')
    
    # Combine legends from all axes
    lines, labels = ax.get_legend_handles_labels()
    if rho is not None:
        l, lab = ax_rho.get_legend_handles_labels()
        lines.extend(l)
        labels.extend(lab)
    if T is not None:
        l, lab = ax_T.get_legend_handles_labels()
        lines.extend(l)
        labels.extend(lab)
    if v is not None:
        l, lab = ax_v.get_legend_handles_labels()
        lines.extend(l)
        labels.extend(lab)
    
    ax.legend(lines, labels, loc='upper right', fontsize=9)
    
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
    
    # Create subplots with shared x-axis
    fig, axes = plt.subplots(2, 3, figsize=(18, 10), sharex=True)
    axes = axes.flatten()
    
    # Density
    ax1 = axes[0]
    ax1.loglog(r_plt, rho_plt, 'b-', linewidth=2)
    ax1.set_ylabel('Density (g cm$^{-3}$)')
    ax1.set_title('Density')
    ax1.grid(True, alpha=0.3, which='both')
    ax1.set_xlim(r_min, r_max)
    
    # Velocity
    ax2 = axes[1]
    ax2.loglog(r_plt, v_plt, 'r-', linewidth=2)
    ax2.set_ylabel('Velocity (km s$^{-1}$)')
    ax2.set_title('Velocity')
    ax2.grid(True, alpha=0.3, which='both')
    ax2.set_xlim(r_min, r_max)
    
    # Temperature
    ax3 = axes[2]
    ax3.loglog(r_plt, T_plt, 'g-', linewidth=2)
    ax3.set_ylabel('Temperature (K)')
    ax3.set_title('Temperature')
    ax3.grid(True, alpha=0.3, which='both')
    ax3.set_xlim(r_min, r_max)
    
    # Composition — elemental (solid) + isotopes (dash-dot)
    ax4 = axes[3]
    ELEM_COMP = [
        ('H',    'H',    '#FF6B6B', '-'),
        ('He',   'He',   '#4ECDC4', '-'),
        ('C',    'CARB', '#45B7D1', '-'),
        ('O',    'OXY',  '#96CEB4', '-'),
        ('Si',   'SIL',  '#FFEAA7', '-'),
        ('IRON', 'IRON', '#DDA0DD', '-'),
    ]
    ISO_COMP = [
        ('NICK', 'Ni56', '#e377c2', '-.'),
        ('COB',  'Co56', '#ff7f0e', '-.'),
    ]
    for key, label, color, ls in ELEM_COMP:
        if key in mf_plt and np.any(mf_plt[key] > 1e-10):
            ax4.plot(r_plt, mf_plt[key], label=label, color=color,
                     linewidth=2, linestyle=ls)
    for key, label, color, ls in ISO_COMP:
        if key in mf_plt and np.any(mf_plt[key] > 1e-10):
            ax4.plot(r_plt, mf_plt[key], label=label, color=color,
                     linewidth=1.5, linestyle=ls)
    ax4.set_ylabel('Mass Fraction')
    ax4.set_title('Composition')
    ax4.legend(fontsize=8, ncol=2)
    ax4.grid(True, alpha=0.3, which='both')
    ax4.set_yscale('log')
    ax4.set_ylim(1e-4, 1.1)
    ax4.set_xlim(r_min, r_max)
    ax4.set_xscale('log')
    
    # Opacity (kappa)
    ax5 = axes[4]
    ax5.loglog(r_plt, kappa_plt, 'k-', linewidth=2)
    ax5.set_ylabel(r'Opacity κ (cm$^2$/g)')
    ax5.set_title(rf'Total Opacity (κ$_{{min}}$={kappa_plt.min():.3e}, κ$_{{max}}$={kappa_plt.max():.3e})')
    ax5.grid(True, alpha=0.3, which='both')
    ax5.set_xlim(r_min, r_max)
    
    # Sigma
    ax6 = axes[5]
    sigma = sigma_velocity(r_cm, v_kms)
    sigma_plt = sigma[::-1]  # Reverse sigma for plotting
    ax6.step(r_plt, sigma_plt, 'g-', linewidth=2, where='post')
    ax6.set_ylabel('Sigma (σ)')
    ax6.set_title('Sigma = dlnV/dlnr - 1')
    ax6.grid(True, alpha=0.3, which='both')
    ax6.set_xlim(r_min, r_max)
    ax6.set_xscale('log')
    ax6.axhline(y=0, color='k', linestyle='--', alpha=0.5, label='σ=0 (homologous)')
    ax6.legend(fontsize=8)
    
    # Set common x-axis label on bottom row
    for ax in axes[3:]:
        ax.set_xlabel(r'r/R$_\odot$')
    
    plt.suptitle(f'STELLA Structure Overview — t = {time_days:.1f} days', fontsize=16, fontweight='bold')
    plt.tight_layout()
    
    plot_file = output_file.replace('.dat', '_plot.png').replace('SN_HYDRO_DATA', 'SN_HYDRO_DATA_plot')
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    print(f"Structure overview plot saved to {plot_file}")
    plt.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('swd_file',    help='STELLA *.swd output file')
    parser.add_argument('output_file', help='Output SN_HYDRO_DATA path')
    parser.add_argument('--abn-file',     default=None)
    parser.add_argument('--input-profile', default=None,
                        help='HERACLES input_profile file for composition fallback')
    parser.add_argument('--xni-file',     default=None,
                        help='STELLA *.xni file for Ni56 distribution')
    parser.add_argument('--lbol-file',    default=None,
                        help='STELLA *.lbol file for actual bolometric luminosity')
    parser.add_argument('--target-days',  type=float, default=None,
                        help='Pick epoch closest to this (default: last epoch)')
    parser.add_argument('--tau-max',      type=float, default=30.0)
    parser.add_argument('--tau-min',      type=float, default=1e-4)
    parser.add_argument('--nzon-target',  type=int,   default=200)
    parser.add_argument('--min-points',   type=int,   default=80,
                        help='Minimum number of points in output SN_HYDRO_DATA')
    parser.add_argument('--kappa',        type=float, default=0.34,
                        help='Electron scattering opacity cm^2/g')
    parser.add_argument('--model-name',   default='stella')
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
    parser.add_argument('--max-logtau-jump', type=float, default=0.10,
                        help='Max log optical depth jump (default: 0.10)')
    parser.add_argument('--max-points',   type=int,   default=500,
                        help='Maximum number of grid points for adaptive refinement (default: 500)')
    parser.add_argument('--t-floor',      type=float, default=3000.0,
                        help='Minimum temperature floor in K (default: 3000)')
    parser.add_argument('--smooth-velocity', action='store_true',
                        help='Enable velocity smoothing to reduce sigma spikes at ejecta interfaces')
    parser.add_argument('--sigma-threshold', type=float, default=2.0,
                        help='Sigma threshold for velocity smoothing (default: 2.0)')
    parser.add_argument('--smooth-width', type=int, default=5,
                        help='Number of points for velocity smoothing window (default: 5)')
    args = parser.parse_args()

    print(f"Reading {args.swd_file}...")
    epochs = read_swd(args.swd_file)
    available = sorted(epochs.keys())
    print(f"Available epochs: {available[0]:.2f} ... {available[-1]:.2f} days "
          f"({len(available)} snapshots)")

    if args.target_days is None:
        t_pick = available[-1]
    else:
        # Handle negative times (pre-SBO) vs positive times (post-SBO)
        # If target_days is positive, prefer positive epochs (post-SBO)
        # If target_days is negative, prefer negative epochs (pre-SBO)
        if args.target_days >= 0:
            # Filter to positive epochs only (post-SBO)
            positive_epochs = [t for t in available if t >= 0]
            if positive_epochs:
                t_pick = min(positive_epochs, key=lambda t: abs(t - args.target_days))
                # If target is far from any available positive epoch, use last positive epoch
                if abs(t_pick - args.target_days) > max(1.0, 0.1 * (positive_epochs[-1] - positive_epochs[0])):
                    t_pick = positive_epochs[-1]
                    print(f"Requested {args.target_days:.2f} d far from available post-SBO range; "
                          f"using last post-SBO epoch")
            else:
                # No positive epochs available, use closest epoch (could be pre-SBO)
                t_pick = min(available, key=lambda t: abs(t - args.target_days))
                print(f"WARNING: No post-SBO epochs available; using closest epoch (pre-SBO)")
        else:
            # Target is negative (pre-SBO), prefer negative epochs
            negative_epochs = [t for t in available if t < 0]
            if negative_epochs:
                t_pick = min(negative_epochs, key=lambda t: abs(t - args.target_days))
                # If target is far from any available negative epoch, use first negative epoch
                if abs(t_pick - args.target_days) > max(1.0, 0.1 * (negative_epochs[-1] - negative_epochs[0])):
                    t_pick = negative_epochs[0]
                    print(f"Requested {args.target_days:.2f} d far from available pre-SBO range; "
                          f"using earliest pre-SBO epoch")
            else:
                # No negative epochs available, use closest epoch (could be post-SBO)
                t_pick = min(available, key=lambda t: abs(t - args.target_days))
                print(f"WARNING: No pre-SBO epochs available; using closest epoch (post-SBO)")
    print(f"Using epoch t = {t_pick:.4f} days "
          f"(requested {args.target_days:.2f})" if args.target_days else
          f"Using last epoch t = {t_pick:.4f} days")

    snap = epochs[t_pick]
    r    = snap['r']
    v    = snap['v']
    T    = snap['T']
    rho  = snap['rho']
    nzon = len(r)
    print(f"Snapshot: {nzon} zones")
    
    # Save original data for mass fraction plot (before trimming)
    orig_r_for_plot = r.copy()
    orig_rho_for_plot = rho.copy()
    orig_T_for_plot = T.copy()

    # Load or construct mass fractions
    if args.abn_file:
        print(f"Reading abundances from {args.abn_file}...")
        abn_raw = read_abn(args.abn_file, nzon)
        # Map to CMFGEN-friendly set
        mass_fractions = {
            'H':  abn_raw['H'],
            'He': abn_raw['He'],
            'C':  abn_raw['C'],
            'O':  abn_raw['O'],
            'Si': abn_raw['Si'],
            'S':  abn_raw['S'],
            'Ca': abn_raw['Ca'],
            'FePeak': abn_raw['FePeak'],
        }
        
        # Check if Ca is zero in .abn file and use input_profile as fallback
        if np.all(mass_fractions['Ca'] == 0) and args.input_profile:
            print(f"WARNING: Ca is zero in .abn file, reading from input_profile...")
            try:
                profile_raw = read_input_profile(args.input_profile)
                # Interpolate profile data to match SWD grid
                from scipy.interpolate import interp1d
                # Read profile radii for interpolation
                profile_r = []
                with open(args.input_profile) as fh:
                    for line in fh:
                        if line.startswith('#'):
                            continue
                        cols = line.split()
                        if len(cols) >= 11:
                            profile_r.append(float(cols[0]))
                profile_r = np.array(profile_r)[::-1]  # outer->inner
                # Interpolate Ca to SWD grid
                interp_ca = interp1d(profile_r, profile_raw['Ca'], kind='linear', 
                                     bounds_error=False, fill_value='extrapolate')
                mass_fractions['Ca'] = interp_ca(r)
                print(f"  Ca interpolated from input_profile: min={np.min(mass_fractions['Ca']):.4e}, max={np.max(mass_fractions['Ca']):.4e}")
            except Exception as e:
                print(f"  WARNING: Failed to read Ca from input_profile: {e}")
                print(f"  Using solar abundance fallback instead")
        
        # Save original mass fractions before trimming (before normalization)
        orig_mass_fractions = {k: v.copy() for k, v in mass_fractions.items()}
        
        # Initialize ni56_for_isotopes (will be filled from .xni or .abn)
        ni56_for_isotopes = np.zeros(nzon)
        
        # Read Ni56 from .xni file if provided (overrides .abn Ni56)
        if args.xni_file:
            print(f"Reading Ni56 distribution from {args.xni_file}...")
            try:
                ni56_from_xni = read_xni(args.xni_file, nzon)
                # Store Ni56 separately for isotopic data, don't add to mass_fractions yet
                ni56_for_isotopes = ni56_from_xni
                print(f"  Ni56 max fraction: {np.max(ni56_from_xni):.4e}")
            except Exception as e:
                print(f"  WARNING: Failed to read .xni file: {e}")
                print(f"  Using Ni56 from .abn file instead")
                ni56_for_isotopes = abn_raw['Ni56'].copy()
        else:
            # Use Ni56 from .abn file
            ni56_for_isotopes = abn_raw['Ni56'].copy()
        
        # Remove Ni56 from elemental mass fractions (will be in isotopic_data)
        # This must happen BEFORE normalization
        if 'Ni56' in mass_fractions:
            del mass_fractions['Ni56']
        
        # Normalise elemental mass fractions (stable material only)
        mf_sum = sum(mass_fractions.values())
        mf_sum = np.where(mf_sum < 1e-10, 1.0, mf_sum)
        mass_fractions = {k: v / mf_sum for k, v in mass_fractions.items()}
        
        # Update orig_mass_fractions to match (remove Ni56, normalize)
        if 'Ni56' in orig_mass_fractions:
            del orig_mass_fractions['Ni56']
        orig_mf_sum = sum(orig_mass_fractions.values())
        orig_mf_sum = np.where(orig_mf_sum < 1e-10, 1.0, orig_mf_sum)
        orig_mass_fractions = {k: v / orig_mf_sum for k, v in orig_mass_fractions.items()}
    else:
        # Simple 2-component H+He placeholder
        print("WARNING: no --abn-file supplied; using H=0.7/He=0.3 everywhere")
        mass_fractions = {
            'H':  np.full(nzon, 0.70),
            'He': np.full(nzon, 0.28),
            'C':  np.full(nzon, 0.01),
            'O':  np.full(nzon, 0.01),
        }
        # Save original mass fractions before trimming
        orig_mass_fractions = {k: v.copy() for k, v in mass_fractions.items()}
        # No Ni56 in this case
        ni56_for_isotopes = np.zeros(nzon)

    # Optical depth trimming with adaptive relaxation to meet minimum point requirement
    # tau is ordered from outside (small tau) to inside (large tau)
    # Use composition-dependent opacity with HERACLES-style Rosseland mean
    kappa_arr = compute_kappa_total(rho, T, mass_fractions)
    tau = compute_kappa_tau(r, rho, kappa_arr)

    print(f"Original tau range: [{tau.min():.3e}, {tau.max():.3e}]")

    # Initial tau trimming (Dessart+2015: truncate at tau_Rosseland=30)
    # STRICT: do NOT relax tau_max if too few points. Adaptive grid will add
    # points where gradients require it, but we keep the tau=30 boundary.
    print(f"[Tau Trim] Dessart+2015 style: truncating at tau_Rosseland={args.tau_max:.1f}")
    tau_max_current = args.tau_max
    tau_min_current = args.tau_min

    # Find outer boundary: first zone where tau >= tau_min
    i_outer = 0
    while i_outer < nzon and tau[i_outer] < tau_min_current:
        i_outer += 1
    # Find inner boundary: last zone where tau <= tau_max
    i_inner = nzon - 1
    while i_inner > 0 and tau[i_inner] > tau_max_current:
        i_inner -= 1
    i_outer = max(0, i_outer - 1)

    n_after_tau = i_inner - i_outer + 1

    if n_after_tau < args.min_points:
        print(f"[Tau Trim] WARNING: only {n_after_tau} points within tau range "
              f"(min_points={args.min_points}), but keeping strict tau_max={tau_max_current}")
        print(f"  Adaptive grid will add points where gradients require resolution.")

    # Inner boundary: zone just behind tau_max (NOT in SN_HYDRO_DATA)
    i_boundary = min(i_inner + 1, nzon - 1)
    r_inner_boundary = r[i_boundary]
    T_inner_boundary = T[i_boundary]
    print(f"Inner boundary zone (not in SN_HYDRO_DATA): depth {i_boundary}, "
          f"T={T_inner_boundary:.1f} K, R={r_inner_boundary:.2e} cm")

    # Final tau trim: exclude the boundary zone
    sl = slice(i_outer, i_inner + 1)
    print(f"Tau trim: keep depths {i_outer}..{i_inner}  "
          f"(tau=[{tau[i_outer]:.2e},{tau[i_inner]:.2e}])")

    r    = r[sl];   v    = v[sl];   T    = T[sl]
    rho  = rho[sl]; kappa_arr = kappa_arr[sl]
    mass_fractions = {k: v2[sl] for k, v2 in mass_fractions.items()}
    
    # Trim Ni56 array to match the trimmed grid
    if 'ni56_for_isotopes' in locals():
        ni56_for_isotopes = ni56_for_isotopes[sl]

    # Save post-tau-trim grid for Ni56 interpolation after grid processing
    r_tau_trim = r.copy()

    # Disable artifact filter - it was removing optically thick zones (high density)
    # The high density zones in the inner region are physical, not artifacts
    print("[Artifact Filter] DISABLED - keeping all zones after tau trim")

    # --- Outer boundary refinement ---
    # The photosphere and spectrum formation region are in the outer layers.
    # Add extra points where tau < 1 to improve Sobolev/flux accuracy.
    def refine_outer_boundary(r, v, T, rho, kappa, mass_fractions, tau, max_extra=20, tau_thresh=1.0):
        """Insert additional points in the outer low-tau region (tau < tau_thresh)."""
        n = len(r)
        # Find the index where tau first exceeds tau_thresh (outer->inner order)
        i_thresh = n
        for i in range(n):
            if tau[i] > tau_thresh:
                i_thresh = i
                break
        # The outer region is indices 0 .. i_thresh-1
        n_outer = i_thresh
        if n_outer < 3:
            print(f"[Outer Refine] Outer region too small ({n_outer} pts), skipping")
            return r, v, T, rho, kappa, mass_fractions, tau

        # We insert points between consecutive pairs in the outer region
        # Limit to avoid bloat
        n_insert_target = min(max_extra, n_outer)
        # Choose evenly spaced insertion points
        insert_indices = []
        step = max(1, n_outer // n_insert_target)
        for i in range(0, n_outer - 1, step):
            insert_indices.append(i)

        # Build new arrays by inserting at midpoints
        def interp_linear(a, b, frac):
            return a + frac * (b - a)
        def interp_log(a, b, frac):
            eps = 1e-30
            return 10 ** (np.log10(max(a, eps)) + frac * (np.log10(max(b, eps)) - np.log10(max(a, eps))))

        r_new, v_new, T_new, rho_new, kappa_new, tau_new = list(r), list(v), list(T), list(rho), list(kappa), list(tau)
        mf_new = {k: list(v) for k, v in mass_fractions.items()}

        # Insert from the end backwards to keep indices valid
        offset = 0
        for i in sorted(insert_indices, reverse=True):
            j = i + offset  # current position in growing list
            if j >= len(r_new) - 1:
                continue
            r_mid = interp_linear(r_new[j], r_new[j+1], 0.5)
            v_mid = interp_linear(v_new[j], v_new[j+1], 0.5)
            T_mid = interp_log(T_new[j], T_new[j+1], 0.5)
            rho_mid = interp_log(rho_new[j], rho_new[j+1], 0.5)
            kappa_mid = interp_log(kappa_new[j], kappa_new[j+1], 0.5)
            tau_mid = interp_log(tau_new[j], tau_new[j+1], 0.5)

            r_new.insert(j+1, r_mid)
            v_new.insert(j+1, v_mid)
            T_new.insert(j+1, T_mid)
            rho_new.insert(j+1, rho_mid)
            kappa_new.insert(j+1, kappa_mid)
            tau_new.insert(j+1, tau_mid)
            for elem in mf_new:
                mf_mid = interp_linear(mf_new[elem][j], mf_new[elem][j+1], 0.5)
                mf_new[elem].insert(j+1, mf_mid)
            offset += 1

        n_added = len(r_new) - n
        print(f"[Outer Refine] Inserted {n_added} points in outer region (tau < {tau_thresh})")
        print(f"  Outer region: {n_outer} -> {n_outer + n_added} points")
        return (np.array(r_new), np.array(v_new), np.array(T_new),
                np.array(rho_new), np.array(kappa_new),
                {k: np.array(v) for k, v in mf_new.items()}, np.array(tau_new))

    # Apply outer boundary refinement before AMR
    tau = compute_kappa_tau(r, rho, kappa_arr)
    r, v, T, rho, kappa_arr, mass_fractions, tau = refine_outer_boundary(
        r, v, T, rho, kappa_arr, mass_fractions, tau,
        max_extra=20, tau_thresh=1.0)

    # Dessart-style adaptive grid workflow
    if args.enable_adaptive_grid:
        print("\n--- Dessart-style Adaptive Grid Workflow ---")

        # Save original grid after tau trim
        orig_r_adapt = r.copy()
        orig_v_adapt = v.copy()
        orig_T_adapt = T.copy()
        orig_rho_adapt = rho.copy()
        orig_tau_adapt = tau.copy()

        n_base = len(r)

        # Step 1: Only downsample if tau-trimmed grid exceeds max_points.
        # Otherwise keep the tau-trimmed resolution as the base grid.
        if n_base > args.max_points:
            print(f"[Base Grid] Tau-trimmed grid has {n_base} points > max_points={args.max_points}, downsampling to {args.base_grid_target}")
            mf_list = list(mass_fractions.values())
            result = downsample_grid(r, v, T, rho, kappa_arr, *mf_list,
                                    target_points=args.base_grid_target)
            n_fixed = 4
            r, v, T, rho, kappa_arr = result[:n_fixed + 1]
            for i, key in enumerate(mass_fractions):
                mass_fractions[key] = result[n_fixed + 1 + i]
            n_base = len(r)
        else:
            print(f"[Base Grid] Keeping tau-trimmed grid: {n_base} points (<= max_points={args.max_points})")

        # Save base grid (before adaptive refinement)
        downsampled_r = r.copy()
        downsampled_v = v.copy()
        downsampled_T = T.copy()
        downsampled_rho = rho.copy()
        downsampled_tau = compute_kappa_tau(r, rho, kappa_arr)

        # Step 2: Calculate electron density and optical depth for adaptive refinement
        n_e = calculate_electron_density(rho, T, mass_fractions)
        tau = compute_kappa_tau(r, rho, kappa_arr)

        # Step 3: Apply adaptive grid refinement (adds points only where gradients are steep)
        r, v, T, rho, n_e, tau, mass_fractions = adaptive_grid_refinement(
            r, v, T, rho, n_e, tau, mass_fractions,
            max_radius_jump=args.max_radius_jump,
            max_temp_jump=args.max_temp_jump,
            max_density_jump=args.max_density_jump,
            max_ne_jump=args.max_ne_jump,
            max_logtau_jump=args.max_logtau_jump,
            max_points=args.max_points
        )

        # Recalculate kappa with refined grid
        kappa_arr = compute_kappa_total(rho, T, mass_fractions)

        # Generate adaptive grid diagnostic plots
        print("\n--- Generating Adaptive Grid Diagnostic Plots ---")
        plot_adaptive_grid_diagnostics(
            orig_r_adapt, orig_v_adapt, orig_T_adapt, orig_rho_adapt, orig_tau_adapt,
            downsampled_r, downsampled_v, downsampled_T, downsampled_rho, downsampled_tau,
            r, v, T, rho, tau,
            args.output_file, t_pick,
            max_radius_jump=args.max_radius_jump,
            max_temp_jump=args.max_temp_jump,
            max_density_jump=args.max_density_jump,
            max_ne_jump=args.max_ne_jump,
            max_logtau_jump=args.max_logtau_jump,
            max_points=args.max_points
        )
        
    elif args.enable_grid_thinning:
        # Legacy grid thinning (if adaptive grid not enabled)
        mf_list = list(mass_fractions.values())
        result = thin_degenerate_grid(r, v, T, rho, kappa_arr, *mf_list, 
                                      min_dr_r=args.grid_thinning_threshold)
        n_fixed = 4
        r, v, T, rho, kappa_arr = result[:n_fixed + 1]
        for i, key in enumerate(mass_fractions):
            mass_fractions[key] = result[n_fixed + 1 + i]
    else:
        print(f"[Grid Processing] DISABLED - keeping all {len(r)} zones")

    # Resample Ni56 to the final grid after downsampling/refinement
    # Both r and r_tau_trim are in outer->inner (decreasing) order.
    # np.interp requires xp to be INCREASING, so we reverse both arrays.
    if len(ni56_for_isotopes) != len(r):
        print(f"[Ni56 Resample] Interpolating Ni56 from {len(ni56_for_isotopes)} to {len(r)} grid points")
        ni56_for_isotopes = np.interp(r, r_tau_trim[::-1], ni56_for_isotopes[::-1])

    # Final check: ensure we meet minimum point requirement
    if len(r) < args.min_points:
        print(f"WARNING: Final grid has {len(r)} points, below minimum {args.min_points}")
        print(f"  Consider increasing STELLA nzon or relaxing tau bounds further")
    else:
        print(f"[Grid Check] Final grid has {len(r)} points >= minimum {args.min_points}")

    # Apply velocity smoothing if requested
    if args.smooth_velocity:
        v = smooth_velocity_discontinuities(r, v, 
                                            sigma_threshold=args.sigma_threshold,
                                            smooth_width=args.smooth_width)

    # Sigma
    sigma = sigma_velocity(r, v)

    # Ensure ni56_for_isotopes is defined (even if empty)
    if 'ni56_for_isotopes' not in locals():
        ni56_for_isotopes = np.zeros(len(r))

    # Luminosity: blackbody at inner boundary (zone just behind tau_max)
    lum_lsun = luminosity_bb(T_inner_boundary, r_inner_boundary) / LSUN
    print(f"Inner boundary (zone behind tau_max): T={T_inner_boundary:.1f} K, "
          f"R={r_inner_boundary:.2e} cm -> L={lum_lsun:.3e} Lsun (blackbody)")
    
    # Also print photospheric luminosity for reference
    tau2 = compute_kappa_tau(r, rho, kappa_arr)
    i_phot = np.searchsorted(tau2, 2.0 / 3.0)
    i_phot = min(i_phot, len(r) - 1)
    T_phot = T[i_phot]
    r_phot = r[i_phot]
    lum_phot_lsun = luminosity_bb(T_phot, r_phot) / LSUN
    print(f"Photosphere (tau=2/3): depth {i_phot}, T={T_phot:.1f} K, "
          f"R={r_phot:.2e} cm -> L={lum_phot_lsun:.3e} Lsun (for reference)")

    write_sn_hydro_data(
        args.output_file, r, v, rho, T, sigma, kappa_arr,
        time_days=t_pick, mass_fractions=mass_fractions, lum_lsun=lum_lsun,
        t_floor=args.t_floor,
        r_inner_cm=r_inner_boundary, T_inner_K=T_inner_boundary,
        ni56_initial=ni56_for_isotopes
    )

    # Create mass_fractions_plot with Co56 and Fe56 for visualization
    # The write function already calculated the decay chain, so we need to recalculate for plotting
    mass_fractions_plot = mass_fractions.copy()
    t_Ni = 8.8
    t_Co = 111.3
    Ni56_initial = ni56_for_isotopes.copy()
    Ni56_t = Ni56_initial * np.exp(-abs(t_pick) / t_Ni)
    Co56_t = Ni56_initial * (t_Co / (t_Co - t_Ni)) * (np.exp(-abs(t_pick) / t_Co) - np.exp(-abs(t_pick) / t_Ni))
    Fe56_t = Ni56_initial - Ni56_t - Co56_t
    mass_fractions_plot['Co56'] = Co56_t
    mass_fractions_plot['Fe56'] = Fe56_t
    
    # Generate diagnostic plots
    print("\n--- Generating diagnostic plots ---")
    
    # Plot comparison (original vs final)
    orig_r = snap['r']
    orig_v = snap['v']
    orig_T = snap['T']
    orig_rho = snap['rho']
    orig_kappa = compute_kappa_total(orig_rho, orig_T, orig_mass_fractions)
    orig_tau = compute_kappa_tau(orig_r, orig_rho, orig_kappa)
    orig_sigma = sigma_velocity(orig_r, orig_v)
    
    new_tau = compute_kappa_tau(r, rho, kappa_arr)
    new_sigma = sigma_velocity(r, v)
    
    plot_comparison(
        orig_r, orig_rho, orig_T, orig_v, orig_tau, orig_sigma,
        r, rho, T, v, new_tau, new_sigma,
        args.output_file, t_pick, tau_min=args.tau_min, tau_max=args.tau_max,
        new_mass_fractions=mass_fractions
    )
    
    # Plot mass fractions: dashed = original STELLA, solid = trimmed CMFGEN
    plot_mass_fractions_stacked(orig_r_for_plot, orig_rho_for_plot, orig_mass_fractions, t_pick, args.output_file,
                                 trimmed_radius_cm=r, trimmed_mass_fractions=mass_fractions_plot)
    
    # Plot dr vs r
    plot_dr_vs_r(r, args.output_file, t_pick, rho=rho, T=T, v=v)
    
    # Plot structure overview
    plot_structure_overview(r, v, rho, T, mass_fractions, kappa_arr, t_pick, args.output_file)


if __name__ == '__main__':
    main()
