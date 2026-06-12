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
    - nzon: Number of zones (default: 500)
    
    Returns:
    - bytes: Binary data for .xni file
    """
    total_mass = nickel_config.get('total_mass_msun', 0.05)
    dist_mass = nickel_config.get('distribution_mass_msun', 0.05)
    dist_type = nickel_config.get('distribution_type', 'uniform')
    
    # Calculate mass fraction
    # If Ni56 is distributed over X Msun, the mass fraction is total_mass / dist_mass
    mass_fraction = total_mass / dist_mass if dist_mass > 0 else 0.0
    
    # Cap at 1.0 (physically impossible to have >100% Ni56)
    mass_fraction = min(mass_fraction, 1.0)
    
    print(f"Generating .xni file:")
    print(f"  Total Ni56 mass: {total_mass} Msun")
    print(f"  Distribution mass range: {dist_mass} Msun")
    print(f"  Ni56 mass fraction: {mass_fraction:.6f}")
    print(f"  Distribution type: {dist_type}")
    
    # Create array of 1002 floats (2 header + 500 zones + 500 padding)
    floats = [0.0] * 1002
    
    # Header floats (based on existing file structure)
    floats[0] = 5.605194e-42  # Metadata float 1
    floats[1] = -1.417409e-34  # Metadata float 2
    
    # Zone values (indices 2-501)
    # Based on existing file: alternating pattern of mass_fraction and near-zero
    for i in range(nzon):
        idx = 2 + i
        if i % 2 == 0:
            # Even indices: mass fraction
            floats[idx] = mass_fraction
        else:
            # Odd indices: near-zero placeholder
            floats[idx] = -1.417409e-34
    
    # Padding/duplicate values (indices 502-1001)
    # Same pattern as zone values
    for i in range(nzon):
        idx = 502 + i
        if i % 2 == 0:
            floats[idx] = mass_fraction
        else:
            floats[idx] = -1.417409e-34
    
    # Convert to binary
    data = struct.pack('f' * 1002, *floats)
    
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
