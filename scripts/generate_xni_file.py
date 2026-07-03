#!/usr/bin/env python3
"""
Generate STELLA .xni file (nickel56 distribution) from YAML configuration.

The .xni file format is binary with 32-bit floats:
- 2 header floats (metadata)
- 500 zone values (Ni56 mass fraction per zone)
- 500 additional floats (padding/duplicate)

Usage:
    python generate_xni_file.py CONFIG_FILE OUTPUT_FILE
"""

import struct
import numpy as np
import yaml
import sys
from pathlib import Path


def generate_xni(nickel_config, nzon=500):
    """
    Generate .xni binary data from nickel56 configuration.

    Parameters:
    - nickel_config: dict with keys:
        - total_mass_msun: Total Ni56 mass in Msun
        - distribution_mass_msun: Mass range over which Ni56 is distributed
        - distribution_type: 'uniform' or 'gaussian'
        - inner_zone_fraction: Fraction of inner zones to place Ni56 (default: 0.4)
          Ni56 is concentrated in the inner metal-rich ejecta, not the outer H/He shell.
    - nzon: Number of zones (default: 500)

    Returns:
    - bytes: Binary data for .xni file
    """
    total_mass = nickel_config.get('total_mass_msun', 0.05)
    dist_mass = nickel_config.get('distribution_mass_msun', 0.05)
    dist_type = nickel_config.get('distribution_type', 'uniform')
    inner_zone_frac = nickel_config.get('inner_zone_fraction', 0.4)

    # Calculate mass fraction within the distribution region
    mass_fraction = total_mass / dist_mass if dist_mass > 0 else 0.0
    mass_fraction = min(mass_fraction, 1.0)

    # Determine zone range for Ni56 (inner zones only)
    n_ni_zones = max(1, int(nzon * inner_zone_frac))
    # Place Ni56 in inner zones (indices 0..n_ni_zones-1 in inner->outer order)
    # STELLA .xni uses inner->outer order (zone 0 = innermost)

    print(f"Generating .xni file:")
    print(f"  Total Ni56 mass: {total_mass} Msun")
    print(f"  Distribution mass range: {dist_mass} Msun")
    print(f"  Ni56 mass fraction (in distribution region): {mass_fraction:.6f}")
    print(f"  Distribution type: {dist_type}")
    print(f"  Inner zone fraction: {inner_zone_frac} ({n_ni_zones}/{nzon} zones)")

    # Create array of (2 + 2*nzon) floats:
    #   2 header floats + nzon zone pairs (value + placeholder)
    # The file stores nzon zones, each zone = 2 floats.
    n_floats = 2 + 2 * nzon
    floats = [0.0] * n_floats

    # Header floats (based on existing file structure)
    floats[0] = 5.605194e-42  # Metadata float 1
    floats[1] = -1.417409e-34  # Metadata float 2

    # Build the zone value array (inner->outer order)
    # Zone 0 = innermost, zone nzon-1 = outermost
    zone_values = np.zeros(nzon)
    for z in range(n_ni_zones):
        zone_values[z] = mass_fraction

    # Write zone pairs: each zone = [value, placeholder]
    for z in range(nzon):
        floats[2 + 2*z] = zone_values[z]
        floats[2 + 2*z + 1] = -1.417409e-34

    # Convert to binary
    data = struct.pack('f' * n_floats, *floats)

    return data


def main():
    if len(sys.argv) != 3:
        print("Usage: python generate_xni_file.py CONFIG_FILE OUTPUT_FILE")
        sys.exit(1)
    
    config_file = Path(sys.argv[1])
    output_file = Path(sys.argv[2])
    
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
    
    # Generate .xni data
    nzon = stella_cfg.get('nzon', 500)
    data = generate_xni(nickel_cfg, nzon)
    
    # Write output file
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, 'wb') as f:
        f.write(data)
    
    print(f"Written .xni file: {output_file} ({len(data)} bytes)")


if __name__ == '__main__':
    main()
