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


def parse_sn_hydro_data_header(sn_hydro_path: Path):
    """Parse key inner-boundary values from a CMFGEN SN_HYDRO_DATA header."""
    values = {
        'time_days': None,
        'lum_lsun': None,
        'r_inner_cm': None,
        't_inner_k': None,
    }
    if not sn_hydro_path.exists():
        return values

    with open(sn_hydro_path) as f:
        for line in f:
            stripped = line.strip()
            if stripped.startswith('Time(days) since explosion:'):
                try:
                    values['time_days'] = float(stripped.split()[-1])
                except ValueError:
                    pass
            elif stripped.startswith('Energy Lum [lsun]:'):
                try:
                    values['lum_lsun'] = float(stripped.split()[-1])
                except ValueError:
                    pass
            elif stripped.startswith('Inner boundary radius [cm]:'):
                try:
                    values['r_inner_cm'] = float(stripped.split()[-1])
                except ValueError:
                    pass
            elif stripped.startswith('Inner boundary temperature [K]:'):
                try:
                    values['t_inner_k'] = float(stripped.split()[-1])
                except ValueError:
                    pass
            # Header values are all before the radius grid section
            if stripped.startswith('Radius grid'):
                break
    return values


def substitute_param_placeholders(text: str, sedona_dir: Path):
    """Fill any remaining @...@ placeholders in param.lua from SN_HYDRO_DATA."""
    # Look for SN_HYDRO_DATA in several candidate locations:
    #   1. sedona_dir/../cmfgen/SN_HYDRO_DATA          (non-timeseries layout)
    #   2. sedona_dir/../cmfgen/<day>/SN_HYDRO_DATA     (timeseries: sedona/5.0days -> cmfgen/5.0days)
    #   3. sedona_dir/SN_HYDRO_DATA                     (alongside param.lua)
    day_label = sedona_dir.name  # e.g. "5.0days"
    model_root = sedona_dir.parent.parent  # e.g. .../X_3
    candidates = [
        model_root / 'cmfgen' / day_label / 'SN_HYDRO_DATA',
        sedona_dir.parent / 'cmfgen' / 'SN_HYDRO_DATA',
        model_root / 'cmfgen' / 'SN_HYDRO_DATA',
        sedona_dir / 'SN_HYDRO_DATA',
    ]
    sn_hydro = None
    for c in candidates:
        if c.exists():
            sn_hydro = c
            break
    if sn_hydro is None:
        # Fall back to the original default for the warning message
        sn_hydro = sedona_dir.parent / 'cmfgen' / 'SN_HYDRO_DATA'
    header = parse_sn_hydro_data_header(sn_hydro)

    # Determine model name from existing .mod file or parent directory
    mod_files = list(sedona_dir.glob('*.mod'))
    model_name = None
    if mod_files:
        model_name = mod_files[0].stem
    if not model_name:
        model_name = sedona_dir.parent.name

    if '@model_name@' in text and model_name:
        text = text.replace('@model_name@', model_name)
        print(f"  Set model_name = {model_name}")

    if '@snapshot_days@' in text and header['time_days'] is not None:
        text = text.replace('@snapshot_days@', f"{header['time_days']:.6f}")
        print(f"  Set snapshot_days = {header['time_days']:.6f} days")

    if header['r_inner_cm'] is not None:
        r_inner = header['r_inner_cm']
        # Ensure core_radius is within the grid bounds:
        #   >= rmin (grid inner boundary, from .mod header)
        #   < r_first_zone (outer radius of first zone, from first data line)
        mod_file = mod_files[0] if mod_files else None
        if mod_file:
            with open(mod_file) as f:
                mod_lines = f.readlines()
            if len(mod_lines) >= 4:
                # Read rmin from .mod header (line 2: "nzon  rmin  texp  n_elems")
                header_parts = mod_lines[1].strip().split()
                if len(header_parts) >= 2:
                    try:
                        rmin = float(header_parts[1])
                        if r_inner < rmin:
                            r_inner = rmin
                    except ValueError:
                        pass
                first_data = mod_lines[3].strip().split()
                if first_data:
                    r_first_zone = float(first_data[0])
                    if r_inner >= r_first_zone:
                        r_inner = r_first_zone * 0.9999
        if '@r_inner@' in text:
            text = text.replace('@r_inner@', f"{r_inner:.6e}")
            print(f"  Set core_radius = {r_inner:.6e} cm")
        else:
            # Also update an already-substituted core_radius line if it is unsafe
            r_match = re.search(r'core_radius\s*=\s*([0-9.eE+-]+)', text)
            if r_match and r_inner != float(r_match.group(1)):
                text = text.replace(r_match.group(0), f'core_radius      = {r_inner:.6e}')
                print(f"  Updated core_radius = {r_inner:.6e} cm")

    if '@temperature@' in text and header['t_inner_k'] is not None:
        text = text.replace('@temperature@', f"{header['t_inner_k']:.2e}")
        print(f"  Set core_temperature = {header['t_inner_k']:.2e} K")

    if '@luminosity@' in text and header['lum_lsun'] is not None:
        L_sun = 3.828e33
        luminosity = header['lum_lsun'] * L_sun
        text = text.replace('@luminosity@', f"{luminosity:.3e}")
        print(f"  Set core_luminosity = {luminosity:.3e} erg/s ({header['lum_lsun']:.3e} Lsun)")

    return text


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
    rmin = None

    if len(mod_lines) >= 4:
        species_line = mod_lines[2].strip()
        try:
            # Species tokens are "Z.A" (e.g. "1.1 2.4 26.56"); take the Z
            mod_species_z = [int(x.split('.')[0]) for x in species_line.split()]
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

        # Read rmin from .mod header (line 2: "nzon  rmin  texp  n_elems")
        header_parts = mod_lines[1].strip().split()
        if len(header_parts) >= 2:
            try:
                rmin = float(header_parts[1])
            except ValueError:
                pass

        first_data = mod_lines[3].strip().split()
        if len(first_data) >= 1:
            r_first_zone = float(first_data[0])

    # Read current param.lua
    with open(param_file) as f:
        text = f.read()

    # Fill any unsubstituted template placeholders using SN_HYDRO_DATA header values
    text = substitute_param_placeholders(text, sedona_dir)

    # Fix core_radius if needed
    # core_radius must be >= rmin (grid inner boundary) so that particles
    # emitted from the core are on the grid (get_zone returns a valid index).
    # core_radius must be < r_first_zone (outer radius of first zone) so that
    # particles are emitted inside the first zone, not outside it.
    r_match = re.search(r'core_radius\s*=\s*([0-9.eE+-]+)', text)
    if r_match and r_first_zone is not None:
        core_radius = float(r_match.group(1))
        if core_radius >= r_first_zone:
            new_radius = r_first_zone * 0.9999
            text = text.replace(r_match.group(0), f'core_radius      = {new_radius:.6e}')
            print(f"  Adjusted core_radius from {core_radius:.6e} to {new_radius:.6e} cm")
            print(f"    (first zone radius = {r_first_zone:.6e} cm)")
        elif rmin is not None and core_radius < rmin:
            # core_radius is inside the grid inner boundary — particles would
            # be emitted off-grid, causing get_zone to return -1 and triggering
            # an assertion failure in lorentz_transform (relativity.cpp:86).
            # Set core_radius to rmin (just inside the first zone, on the grid).
            new_radius = rmin
            text = text.replace(r_match.group(0), f'core_radius      = {new_radius:.6e}')
            print(f"  Adjusted core_radius from {core_radius:.6e} to {new_radius:.6e} cm")
            print(f"    (rmin = {rmin:.6e} cm; was inside grid inner boundary)")

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
