#!/usr/bin/env python3
"""
Fix SEDONA param.lua after STELLA job generates it.

Reads the .mod file to detect zero-abundance species, filters them from
opacity_atoms_in_nlte, and ensures core_radius is strictly less than the
first zone radius to prevent zones inside the core.
"""

import sys
import re
from pathlib import Path


def fix_param_lua(sedona_dir: Path):
    """Fix param.lua in the given SEDONA directory."""
    sedona_dir = Path(sedona_dir)
    param_file = sedona_dir / 'param.lua'
    mod_file = next(sedona_dir.glob('*.mod'), None)

    if not param_file.exists():
        print(f"WARNING: param.lua not found at {param_file}")
        return False

    if not mod_file:
        print(f"WARNING: No .mod file found in {sedona_dir}")
        return False

    # Read .mod file to detect zero-abundance species and first zone radius
    with open(mod_file) as f:
        mod_lines = f.readlines()

    mod_species_z = []
    zero_abundance_z = set()
    r_first_zone = None

    if len(mod_lines) >= 4:
        species_line = mod_lines[2].strip()
        try:
            mod_species_z = [int(x) for x in species_line.split()]
        except ValueError:
            mod_species_z = []

        for col, z in enumerate(mod_species_z):
            has_mass = False
            for line in mod_lines[3:]:
                parts = line.strip().split()
                if len(parts) >= 5 + col:
                    if float(parts[4 + col]) > 0:
                        has_mass = True
                        break
            if not has_mass:
                zero_abundance_z.add(z)

        first_data = mod_lines[3].strip().split()
        if len(first_data) >= 1:
            r_first_zone = float(first_data[0])

    # Read current param.lua
    with open(param_file) as f:
        text = f.read()

    # Fix core_radius if needed
    r_match = re.search(r'core_radius\s*=\s*([0-9.eE+-]+)', text)
    if r_match and r_first_zone is not None:
        core_radius = float(r_match.group(1))
        if core_radius >= r_first_zone:
            new_radius = r_first_zone * 0.9999
            text = text.replace(r_match.group(0), f'core_radius      = {new_radius:.6e}')
            print(f"  Adjusted core_radius from {core_radius:.6e} to {new_radius:.6e} cm")
            print(f"    (first zone radius = {r_first_zone:.6e} cm)")

    # Filter species from NLTE list that are missing from .mod or have zero abundance
    missing_from_mod = set()
    if mod_species_z:
        nlte_match = re.search(r'opacity_atoms_in_nlte\s*=\s*\{([^}]*)\}', text)
        if nlte_match:
            all_z = [int(x.strip()) for x in nlte_match.group(1).split(',') if x.strip()]
            missing_from_mod = set(all_z) - set(mod_species_z)

    species_to_remove = zero_abundance_z | missing_from_mod
    if species_to_remove:
        nlte_match = re.search(r'opacity_atoms_in_nlte\s*=\s*\{([^}]*)\}', text)
        if nlte_match:
            original = nlte_match.group(0)
            all_z = [int(x.strip()) for x in nlte_match.group(1).split(',') if x.strip()]
            filtered = [z for z in all_z if z not in species_to_remove]
            new_nlte = f'opacity_atoms_in_nlte = {{{", ".join(str(z) for z in filtered)}}}'
            text = text.replace(original, new_nlte)
            print(f"  Filtered opacity_atoms_in_nlte: removed {sorted(species_to_remove)}")
            print(f"    New list: {filtered}")

    with open(param_file, 'w') as f:
        f.write(text)

    print(f"  Fixed {param_file}")
    return True


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: fix_sedona_param.py <sedona_directory>")
        sys.exit(1)

    sedona_dir = Path(sys.argv[1])
    print(f"Fixing SEDONA param.lua in {sedona_dir}...")
    fix_param_lua(sedona_dir)
