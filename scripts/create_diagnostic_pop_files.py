#!/usr/bin/env python3
"""
Create diagnostic POP files from *_IN population files for observer testing.

This is a DIAGNOSTIC tool, NOT a science product. The POP files created
here use intermediate (non-converged) populations from the *_IN files.
They are only for testing the observer K-index fix with non-monotonic
velocity fields.

Usage:
    python create_diagnostic_pop_files.py <cmfgen_dir>

Where <cmfgen_dir> contains the *_IN files and RVTJ.
"""

import sys
import os
import re
import math
from pathlib import Path

# Species mapping: POP file prefix -> (species_abbr, list of (ion_id, in_file_suffix))
# GEN_ION_ID for observer: '0','I','2','III','IV','V','SIX','SEV','VIII','IX','X',...
# J=1:'0', J=2:'I', J=3:'2', J=4:'III', J=5:'IV', J=6:'V', ...
SPECIES_MAP = {
    'HYD':  {'abbr': 'H',  'ions': [('HI', 'HI_IN'), ('H2', None)]},
    'HE':   {'abbr': 'He', 'ions': [('HeI', 'HeI_IN'), ('He2', None)]},
    'OXY':  {'abbr': 'O',  'ions': [('OI', 'OI_IN'), ('O2', 'O2_IN')]},
    'SIL':  {'abbr': 'Sk', 'ions': [('SkI', 'SkI_IN'), ('Sk2', 'Sk2_IN')]},
    'CAL':  {'abbr': 'Ca', 'ions': [('CaI', 'CaI_IN'), ('Ca2', 'Ca2_IN')]},
    'IRON': {'abbr': 'Fe', 'ions': [('FeI', 'FeI_IN'), ('Fe2', 'Fe2_IN'),
                                     ('FeIII', 'FeIII_IN'), ('FeIV', 'FeIV_IN')]},
}

# Atomic weights for abundance computation
ATOMIC_WEIGHTS = {
    'HYD': 1.008,
    'HE': 4.0026,
    'OXY': 15.999,
    'SIL': 28.085,
    'CAL': 40.078,
    'IRON': 55.845,
}

# Format date that RD_POP_DIM accepts
FORMAT_DATE = '27-JAN-1992'

def parse_in_file(filepath):
    """Parse a *_IN population file.
    
    Returns (nlevels, nd, populations) where populations is a list of nd lists,
    each containing nlevels population values.
    """
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    # Find the line with NLEVELS and ND
    # Format: V_INF  0.0  NLEVELS  ND
    nlevels = None
    nd = None
    data_start = 0
    for i, line in enumerate(lines):
        parts = line.split()
        if len(parts) >= 4:
            try:
                v1 = float(parts[0])
                v2 = float(parts[1])
                nl = int(parts[2])
                nd_val = int(parts[3])
                if nl > 0 and nd_val > 0 and nl < 1000 and nd_val < 10000:
                    nlevels = nl
                    nd = nd_val
                    data_start = i + 1
                    break
            except (ValueError, IndexError):
                continue
    
    if nlevels is None or nd is None:
        raise ValueError(f"Could not parse header from {filepath}")
    
    # Parse population data
    # After the header, for each depth:
    #   - 1 header line (R, density, T, ...)
    #   - ceil(nlevels/5) population lines
    #   - 1 blank line (usually)
    
    populations = []
    idx = data_start
    
    for d in range(nd):
        # Skip blank lines
        while idx < len(lines) and lines[idx].strip() == '':
            idx += 1
        
        # Skip header line (starts with a number - R value)
        if idx < len(lines):
            idx += 1
        
        # Read nlevels population values
        level_pops = []
        while len(level_pops) < nlevels and idx < len(lines):
            line = lines[idx].strip()
            if line == '':
                idx += 1
                continue
            parts = line.split()
            for p in parts:
                try:
                    level_pops.append(float(p))
                except ValueError:
                    break
            idx += 1
        
        if len(level_pops) < nlevels:
            print(f"Warning: depth {d+1} in {filepath}: expected {nlevels} levels, got {len(level_pops)}")
            level_pops.extend([0.0] * (nlevels - len(level_pops)))
        
        populations.append(level_pops[:nlevels])
    
    return nlevels, nd, populations


def compute_ion_populations(level_pops):
    """Compute ion population (sum of level populations) at each depth."""
    return [sum(levels) for levels in level_pops]


def read_mass_fractions(scratch_dir):
    """Read mass fractions from SN_HYDRO_DATA."""
    hydro_file = os.path.join(scratch_dir, 'SN_HYDRO_DATA')
    if not os.path.exists(hydro_file):
        return {}
    
    with open(hydro_file, 'r') as f:
        content = f.read()
    
    # Find mass fraction sections
    mass_fracs = {}
    species_names = ['HYD', 'HE', 'OXY', 'SIL', 'CAL', 'IRON']
    
    for sp in species_names:
        # Look for "<SP> mass fraction" followed by values
        pattern = rf'{sp} mass fraction\s*\n((?:\s*[\d.E+\-]+\s*)+)'
        match = re.search(pattern, content)
        if match:
            values = match.group(1).split()
            mass_fracs[sp] = [float(v) for v in values]
    
    return mass_fracs


def compute_abundance(mass_fracs, species):
    """Compute number abundance from mass fractions.
    
    Returns the abundance at the outermost depth point (depth 1).
    """
    if not mass_fracs or species not in mass_fracs:
        return 1.0e-4  # Default
    
    # Compute number fraction at depth 1
    total_num = 0.0
    for sp, fracs in mass_fracs.items():
        if sp in ATOMIC_WEIGHTS and len(fracs) > 0:
            total_num += fracs[0] / ATOMIC_WEIGHTS[sp]
    
    if total_num == 0:
        return 1.0e-4

    abund = mass_fracs[species][0] / ATOMIC_WEIGHTS[species] / total_num
    # Ensure non-zero abundance to avoid division-by-zero in observer
    if abund == 0.0:
        abund = 1.0e-30
    return abund


def write_pop_file(filepath, species, nd, total_pop, ions_data, abundance, completion_time):
    """Write a POP file in the format expected by the observer.
    
    ions_data: list of (ion_id, nlevels, per_level_pops, ion_pops)
        per_level_pops: list of nd lists, each with nlevels values
        ion_pops: list of nd values
    """
    with open(filepath, 'w') as f:
        # Header - abundance must be at column 30 for Fortran T30 format
        f.write(f' Output format date:         {FORMAT_DATE}         \n')
        f.write(f' Completion of Model:        {completion_time}\n')
        f.write(f' ND:                           {nd}\n')
        # Species name + " abundance:" then pad to column 30, then the number
        label = f'{species} abundance:'
        # Pad to column 29 (0-indexed), so number starts at column 30 (1-indexed)
        f.write(f'{label:<29}{abundance:16.9E}\n')
        
        # Total species population (8 per line)
        f.write('    ')
        for i, p in enumerate(total_pop):
            f.write(f'   {p:14.7E}')
            if (i + 1) % 8 == 0 and i < len(total_pop) - 1:
                f.write('\n    ')
        f.write('\n')
        
        # Per-ion data
        for ion_id, nlevels, per_level_pops, ion_pops in ions_data:
            # Number of levels must be at column 30 for Fortran T30,I8 format
            label = f' Number of {ion_id} levels:'
            f.write(f'{label:<29}{nlevels:8d}\n')
            # Oscillator date must be at column 30 for Fortran T30,A format
            label = f'{ion_id} Oscillator date:'
            f.write(f'{label:<29}12-OCT-1995\n')
            
            # Per-level populations (nlevels * nd values, 8 per line)
            f.write('    ')
            count = 0
            for d in range(nd):
                for l in range(nlevels):
                    f.write(f'   {per_level_pops[d][l]:14.7E}')
                    count += 1
                    if count % 8 == 0:
                        f.write('\n    ')
            if count % 8 != 0:
                f.write('\n')
            
            # Ion populations (nd values, 8 per line)
            f.write('    ')
            for i, p in enumerate(ion_pops):
                f.write(f'   {p:14.7E}')
                if (i + 1) % 8 == 0 and i < len(ion_pops) - 1:
                    f.write('\n    ')
            f.write('\n')


def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <cmfgen_dir>")
        sys.exit(1)
    
    cmfgen_dir = sys.argv[1]
    cmfgen_dir = os.path.abspath(cmfgen_dir)
    
    if not os.path.isdir(cmfgen_dir):
        print(f"Error: {cmfgen_dir} is not a directory")
        sys.exit(1)
    
    # Read mass fractions for abundance computation
    mass_fracs = read_mass_fractions(cmfgen_dir)
    
    # Get completion time from RVTJ
    rvtj_path = os.path.join(cmfgen_dir, 'RVTJ')
    completion_time = '09-Sep-2026 17:19:40'
    if os.path.exists(rvtj_path):
        with open(rvtj_path, 'r') as f:
            for line in f:
                if 'Completion of Model:' in line:
                    # Extract time from column 30
                    completion_time = line[29:49].strip()
                    break
    
    print(f"Creating diagnostic POP files in {cmfgen_dir}")
    print(f"Completion time: {completion_time}")
    print(f"WARNING: These POP files use intermediate (non-converged) populations!")
    print()
    
    for species, info in SPECIES_MAP.items():
        pop_file = os.path.join(cmfgen_dir, f'POP{species}')
        
        ions_data = []
        total_pop = None
        nd = None
        
        for ion_id, in_suffix in info['ions']:
            if in_suffix is None:
                continue
            
            in_file = os.path.join(cmfgen_dir, in_suffix)
            if not os.path.exists(in_file):
                print(f"  Skipping {ion_id}: {in_suffix} not found")
                continue
            
            try:
                nlevels, nd_parsed, level_pops = parse_in_file(in_file)
                if nd is None:
                    nd = nd_parsed
                elif nd != nd_parsed:
                    print(f"  Warning: ND mismatch in {in_suffix}: {nd_parsed} vs {nd}")
                
                ion_pops = compute_ion_populations(level_pops)
                ions_data.append((ion_id, nlevels, level_pops, ion_pops))
                
                # Accumulate total population
                if total_pop is None:
                    total_pop = list(ion_pops)
                else:
                    for i in range(min(len(total_pop), len(ion_pops))):
                        total_pop[i] += ion_pops[i]
                
                print(f"  {ion_id}: {nlevels} levels, {nd_parsed} depths")
                
            except Exception as e:
                print(f"  Error parsing {in_suffix}: {e}")
        
        if not ions_data or total_pop is None or nd is None:
            print(f"  No valid ion data for {species}, skipping")
            continue
        
        # Compute abundance
        abundance = compute_abundance(mass_fracs, species)
        
        # Write POP file
        write_pop_file(pop_file, species, nd, total_pop, ions_data, abundance, completion_time)
        print(f"  Written: {pop_file} (abundance={abundance:.4e}, total_pop[0]={total_pop[0]:.4e})")
        print()
    
    print("Done. These POP files are DIAGNOSTIC ONLY - not science products.")


if __name__ == '__main__':
    main()
