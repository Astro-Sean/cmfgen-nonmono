#!/usr/bin/env python3
"""
Generate STELLA .xni file (nickel56 distribution) from YAML configuration.

The .xni file is Fortran unformatted (form='unformatted'), containing a
single record of Nzon REAL*8 (double precision) values representing the
Ni56 mass fraction in each zone.  STELLA reads this with:
    read(28) (XNifor(i), i=1, Nzon)
where XNifor is REAL*8 (from IMPLICIT REAL*8(A-H,O-Z)).

Fortran unformatted format: [4-byte rec_len] [data] [4-byte rec_len]

Usage:
    python generate_xni_file.py CONFIG_FILE OUTPUT_FILE
"""

import struct
import numpy as np
import yaml
import sys
from pathlib import Path


def load_hyd_zone_masses(hyd_path, nzon):
    """Read STELLA .hyd zone masses for exact Ni56 placement.

    .hyd row format (write_x_hyd in prepare_SN_B_model.py):
        header: t_start  nzon  mass_cut  rcen  rho_cen
        zones:  km  dMr(Msun)  r(cm)  rho  T  u  Mr(Msun, incl. mass_cut)  0

    Returns (dMr, M_ejecta) in Msun, inner->outer order, or None on failure.
    """
    try:
        with open(hyd_path) as f:
            header = f.readline().split()
            mass_cut = float(header[2])
            dMr, Mr = [], []
            for line in f:
                cols = line.split()
                if len(cols) >= 7:
                    dMr.append(float(cols[1]))
                    Mr.append(float(cols[6]))
        if len(dMr) != nzon:
            print(f"WARNING: .hyd has {len(dMr)} zones, expected {nzon} - ignoring")
            return None
        return np.asarray(dMr), np.asarray(Mr) - mass_cut
    except Exception as e:
        print(f"WARNING: could not read .hyd zone masses from {hyd_path}: {e}")
        return None


def _xni_masks_from_hyd(ni_boundaries_msun, dMr, M_ejecta):
    """Zone masks for mass-bounded Ni56 regions using actual .hyd masses.

    ni_boundaries_msun: cumulative ejecta-mass boundaries, e.g. [0.5] for a
    single inner region or [dist_a, dist_a+dist_b] for the dual case.
    Returns list of bool masks (inner->outer) covering (prev, boundary].
    """
    masks = []
    lo = 0.0
    for hi in ni_boundaries_msun:
        masks.append((M_ejecta > lo) & (M_ejecta <= hi))
        lo = hi
    return masks


def generate_xni(nickel_config, input_profile_path=None, nzon=500, hyd_path=None):
    """
    Generate .xni binary data from nickel56 configuration.
    Supports dual distribution for ejecta A and ejecta B.

    Parameters:
    - nickel_config: dict with keys:
        - ejecta_a: dict for inner ejecta (fresh Ni56)
        - ejecta_b: dict for outer ejecta (decayed Co56)
        OR legacy single distribution:
        - total_mass_msun: Total Ni56 mass in Msun
        - distribution_mass_msun: Mass range over which Ni56 is distributed
        - distribution_type: 'uniform' or 'gaussian'
        - inner_zone_fraction: Fraction of inner zones to place Ni56 (default: 0.4)
    - input_profile_path: Path to input_profile file for mass-based distribution.
    - nzon: Number of zones (default: 500)

    Returns:
    - np.ndarray: Ni56 mass fraction per zone (length nzon, dtype float64)
    """
    hyd_masses = load_hyd_zone_masses(hyd_path, nzon) if hyd_path else None

    # Check for dual distribution (ejecta_a + ejecta_b)
    if 'ejecta_a' in nickel_config and 'ejecta_b' in nickel_config:
        return generate_xni_dual(nickel_config, input_profile_path, nzon,
                                 hyd_masses=hyd_masses)

    # Legacy single distribution
    total_mass = nickel_config.get('total_mass_msun', 0.05)
    dist_mass = nickel_config.get('distribution_mass_msun', 0.05)
    dist_type = nickel_config.get('distribution_type', 'uniform')
    inner_zone_frac = nickel_config.get('inner_zone_fraction', 0.4)

    # Calculate mass fraction within the distribution region
    mass_fraction = total_mass / dist_mass if dist_mass > 0 else 0.0
    mass_fraction = min(mass_fraction, 1.0)

    zone_values = np.zeros(nzon, dtype=np.float64)

    # Preferred path: exact mass-based placement on the real .hyd grid.
    # The input_profile path maps onto an assumed log-spaced radius grid,
    # which mismatches the actual .hyd zoning and over/under-estimates the
    # Ni56 mass.  .hyd gives per-zone dMr/Mr directly.
    if hyd_masses is not None:
        dMr, M_ejecta = hyd_masses
        mask = _xni_masks_from_hyd([dist_mass], dMr, M_ejecta)[0]
        m_masked = dMr[mask].sum()
        if mask.any() and m_masked > 0:
            frac = min(total_mass / m_masked, 1.0)
            zone_values[mask] = frac
            print(f"Using .hyd mass-based distribution: {hyd_path}")
            print(f"  Total Ni56 mass: {total_mass} Msun")
            print(f"  Ni56 mass fraction: {frac:.6f} over {mask.sum()}/{nzon} zones")
            print(f"  Ni56 covers inner {M_ejecta[mask][-1]:.4f} Msun of ejecta")
            print(f"  Achieved Ni56 mass: {frac * m_masked:.5f} Msun")
            return zone_values
        print("WARNING: .hyd mask empty - falling back to input_profile mapping")

    # Determine zone range for Ni56
    if input_profile_path and Path(input_profile_path).exists():
        # Mass-based distribution: read [DEPRECATED] HERACLES input_profile and find zones within mass range
        # Note: STELLA .hyd format is the active path; HERACLES profile format is deprecated
        print(f"Using mass-based distribution from input_profile: {input_profile_path}")
        r_heracles = []
        rho_heracles = []
        with open(input_profile_path) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                cols = line.split()
                if len(cols) >= 2:
                    r_heracles.append(float(cols[0]))
                    rho_heracles.append(float(cols[1]))
        
        r_heracles = np.array(r_heracles)
        rho_heracles = np.array(rho_heracles)
        
        # Calculate cumulative mass from inner boundary outward
        dr = np.concatenate([[r_heracles[0]], np.diff(r_heracles)])
        r_lo = np.concatenate([[0.0], r_heracles[:-1]])
        dm = (4.0/3.0) * np.pi * (r_heracles**3 - r_lo**3) * rho_heracles
        dm = np.maximum(dm, 0.0)
        m_cum = np.cumsum(dm) / 1.989e33  # cumulative mass in Msun
        
        # Find zones within the Ni56 distribution mass range
        ni56_mask = m_cum <= dist_mass
        n_ni_heracles = np.sum(ni56_mask)
        
        # Map [DEPRECATED] HERACLES zones to STELLA zones
        # STELLA uses log-spaced zones, so we need to interpolate
        # Note: project uses STELLA .hyd format as the active path
        r_stella = np.logspace(np.log10(r_heracles[0]), np.log10(r_heracles[-1]), nzon)
        log_r_heracles = np.log10(r_heracles)
        log_r_stella = np.log10(r_stella)
        
        # Find which STELLA zones are within the mass range
        r_ni_max = r_heracles[ni56_mask][-1] if len(ni56_mask) > 0 else r_heracles[0]
        ni56_stella_mask = r_stella <= r_ni_max
        n_ni_zones = np.sum(ni56_stella_mask)
        
        print(f"  Total Ni56 mass: {total_mass} Msun")
        print(f"  Distribution mass range: {dist_mass} Msun")
        print(f"  Ni56 mass fraction (in distribution region): {mass_fraction:.6f}")
        print(f"  Distribution type: {dist_type}")
        print(f"  Ni56 extends to radius: {r_ni_max:.2e} cm")
        print(f"  Ni56 extends to mass: {m_cum[ni56_mask][-1]:.3f} Msun")
        print(f"  STELLA zones with Ni56: {n_ni_zones}/{nzon}")
    else:
        # Zone-based distribution (fallback)
        print(f"Using zone-based distribution (input_profile not provided)")
        n_ni_zones = max(1, int(nzon * inner_zone_frac))
        
        print(f"  Total Ni56 mass: {total_mass} Msun")
        print(f"  Distribution mass range: {dist_mass} Msun")
        print(f"  Ni56 mass fraction (in distribution region): {mass_fraction:.6f}")
        print(f"  Distribution type: {dist_type}")
        print(f"  Inner zone fraction: {inner_zone_frac} ({n_ni_zones}/{nzon} zones)")

    # Build the zone value array (inner->outer order)
    # Zone 0 = innermost, zone nzon-1 = outermost
    zone_values = np.zeros(nzon, dtype=np.float64)
    if input_profile_path and Path(input_profile_path).exists():
        # Use mass-based mask
        zone_values[ni56_stella_mask] = mass_fraction
    else:
        # Use zone-based distribution
        for z in range(n_ni_zones):
            zone_values[z] = mass_fraction

    return zone_values


def generate_xni_dual(nickel_config, input_profile_path=None, nzon=500,
                    hyd_masses=None):
    """
    Generate .xni binary data with dual distribution (ejecta A + ejecta B).

    Ejecta A: Inner H-rich shell with fresh Ni56
    Ejecta B: Outer metal-rich ejecta with decayed Co56 from 400-day old SN
    """
    ejecta_a_cfg = nickel_config.get('ejecta_a', {})
    ejecta_b_cfg = nickel_config.get('ejecta_b', {})

    # Ejecta A parameters
    mass_a = ejecta_a_cfg.get('total_mass_msun', 0.05)
    dist_mass_a = ejecta_a_cfg.get('distribution_mass_msun', 3.0)
    zone_frac_a = ejecta_a_cfg.get('inner_zone_fraction', 0.4)

    # Ejecta B parameters
    mass_b = ejecta_b_cfg.get('total_mass_msun', 0.0017)
    dist_mass_b = ejecta_b_cfg.get('distribution_mass_msun', 2.0)
    zone_frac_b = ejecta_b_cfg.get('inner_zone_fraction', 1.0)

    # Calculate mass fractions
    mass_frac_a = mass_a / dist_mass_a if dist_mass_a > 0 else 0.0
    mass_frac_b = mass_b / dist_mass_b if dist_mass_b > 0 else 0.0
    mass_frac_a = min(mass_frac_a, 1.0)
    mass_frac_b = min(mass_frac_b, 1.0)

    print("Dual Ni56 distribution:")
    print(f"  Ejecta A (inner): {mass_a} Msun over {dist_mass_a} Msun (fraction: {mass_frac_a:.6f})")
    print(f"  Ejecta B (outer): {mass_b} Msun over {dist_mass_b} Msun (fraction: {mass_frac_b:.6f})")

    # Preferred path: exact mass-based placement on the real .hyd grid.
    if hyd_masses is not None:
        dMr, M_ejecta = hyd_masses
        mask_a, mask_b = _xni_masks_from_hyd(
            [dist_mass_a, dist_mass_a + dist_mass_b], dMr, M_ejecta)
        zone_values = np.zeros(nzon, dtype=np.float64)
        m_a = dMr[mask_a].sum()
        m_b = dMr[mask_b].sum()
        if mask_a.any() and m_a > 0:
            zone_values[mask_a] = min(mass_a / m_a, 1.0)
        if mask_b.any() and m_b > 0:
            zone_values[mask_b] = min(mass_b / m_b, 1.0)
        print(f"Using .hyd mass-based distribution ({mask_a.sum()} ejecta-A "
              f"zones, {mask_b.sum()} ejecta-B zones)")
        print(f"  Achieved Ni56 masses: A={(zone_values * dMr)[mask_a].sum():.5f} Msun, "
              f"B={(zone_values * dMr)[mask_b].sum():.5f} Msun")
        return zone_values

    # Determine zone ranges
    if input_profile_path and Path(input_profile_path).exists():
        # [DEPRECATED] HERACLES input_profile format; STELLA .hyd is the active path
        print(f"Using mass-based distribution from input_profile: {input_profile_path}")
        r_heracles = []
        rho_heracles = []
        with open(input_profile_path) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                cols = line.split()
                if len(cols) >= 2:
                    r_heracles.append(float(cols[0]))
                    rho_heracles.append(float(cols[1]))
        
        r_heracles = np.array(r_heracles)
        rho_heracles = np.array(rho_heracles)
        
        # Calculate cumulative mass
        dr = np.concatenate([[r_heracles[0]], np.diff(r_heracles)])
        r_lo = np.concatenate([[0.0], r_heracles[:-1]])
        dm = (4.0/3.0) * np.pi * (r_heracles**3 - r_lo**3) * rho_heracles
        dm = np.maximum(dm, 0.0)
        m_cum = np.cumsum(dm) / 1.989e33  # cumulative mass in Msun
        
        # Map to STELLA zones (from [DEPRECATED] HERACLES profile; STELLA .hyd is active path)
        r_stella = np.logspace(np.log10(r_heracles[0]), np.log10(r_heracles[-1]), nzon)
        
        # Ejecta A: inner regions (up to dist_mass_a)
        idx_a = min(np.searchsorted(m_cum, dist_mass_a), len(r_heracles) - 1)
        mask_a = r_stella <= r_heracles[idx_a]
        n_zones_a = np.sum(mask_a)
        
        # Ejecta B: outer regions (from dist_mass_a to dist_mass_a + dist_mass_b)
        # Clamp to available mass
        total_mass_available = m_cum[-1]
        mass_b_max = min(dist_mass_a + dist_mass_b, total_mass_available)
        idx_b_min = min(np.searchsorted(m_cum, dist_mass_a), len(r_heracles) - 1)
        idx_b_max = min(np.searchsorted(m_cum, mass_b_max), len(r_heracles) - 1)
        r_b_min = r_heracles[idx_b_min]
        r_b_max = r_heracles[idx_b_max]
        mask_b = (r_stella > r_b_min) & (r_stella <= r_b_max)
        n_zones_b = np.sum(mask_b)
        
        print(f"  Total mass available: {total_mass_available:.3f} Msun")
        print(f"  Ejecta B mass range: {dist_mass_a:.3f} to {mass_b_max:.3f} Msun")
        
        print(f"  Ejecta A zones: {n_zones_a}/{nzon} (inner)")
        print(f"  Ejecta B zones: {n_zones_b}/{nzon} (outer)")
    else:
        # Zone-based distribution
        print(f"Using zone-based distribution (input_profile not provided)")
        n_zones_a = max(1, int(nzon * zone_frac_a))
        n_zones_b = max(1, int(nzon * zone_frac_b))
        mask_a = np.zeros(nzon, dtype=bool)
        mask_b = np.zeros(nzon, dtype=bool)
        mask_a[:n_zones_a] = True
        # Ejecta B in outer zones (last n_zones_b zones)
        mask_b[-n_zones_b:] = True
        print(f"  Ejecta A zones: {n_zones_a}/{nzon} (inner)")
        print(f"  Ejecta B zones: {n_zones_b}/{nzon} (outer)")
    
    # Create zone value array
    zone_values = np.zeros(nzon, dtype=np.float64)
    zone_values[mask_a] = mass_frac_a
    zone_values[mask_b] = mass_frac_b

    return zone_values


def main():
    if len(sys.argv) < 3:
        print("Usage: python generate_xni_file.py CONFIG_FILE OUTPUT_FILE [INPUT_PROFILE] [--hyd-file HYD_FILE]")
        print("  CONFIG_FILE: YAML configuration file")
        print("  OUTPUT_FILE: Output .xni file path")
        print("  INPUT_PROFILE: Optional input_profile path for mass-based distribution")
        print("  --hyd-file: Optional STELLA .hyd path; uses actual zone masses")
        print("              (dMr/Mr columns) for exact mass-based Ni56 placement")
        sys.exit(1)

    # Split positional args from --hyd-file option (and its value)
    hyd_file = None
    skip = set()
    for i, a in enumerate(sys.argv[1:], start=1):
        if a == '--hyd-file' and i + 1 < len(sys.argv):
            hyd_file = Path(sys.argv[i + 1])
            skip.update([i, i + 1])
        elif a.startswith('--hyd-file='):
            hyd_file = Path(a.split('=', 1)[1])
            skip.add(i)
    args = [a for i, a in enumerate(sys.argv[1:], start=1) if i not in skip]

    config_file = Path(args[0])
    output_file = Path(args[1])
    input_profile_path = Path(args[2]) if len(args) > 2 else None

    # Load YAML config
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)

    # Extract nickel56 configuration
    stella_cfg = config.get('stella', {})
    nickel_cfg = stella_cfg.get('nickel56', {})

    if not nickel_cfg.get('enabled', False):
        print("WARNING: nickel56.enabled is false, generating zero-nickel .xni file")
        nickel_cfg['total_mass_msun'] = 0.0
        nickel_cfg['distribution_mass_msun'] = 1.0

    # Generate .xni data (returns zone_values array, not raw bytes)
    nzon = stella_cfg.get('nzon', 500)
    zone_values = generate_xni(nickel_cfg, input_profile_path, nzon,
                               hyd_path=hyd_file)

    # Write as Fortran unformatted: [rec_len] [Nzon * REAL*8] [rec_len]
    output_file.parent.mkdir(parents=True, exist_ok=True)
    rec_len = nzon * 8
    with open(output_file, 'wb') as f:
        f.write(struct.pack('i', rec_len))
        f.write(struct.pack(f'{nzon}d', *zone_values))
        f.write(struct.pack('i', rec_len))

    print(f"Written .xni file: {output_file} ({4 + rec_len + 4} bytes, Fortran unformatted REAL*8)")
    print(f"  Non-zero zones: {np.sum(zone_values > 1e-30)}/{nzon}")
    print(f"  Max XNifor: {zone_values.max():.6f}")


if __name__ == '__main__':
    main()
