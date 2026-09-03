#!/usr/bin/env python3
"""
Convert CMFGEN SN_HYDRO_DATA to SIROCCO 1D spherical import format.

SIROCCO spherical grid format (per-cell columns):
  i  inwind  r  v_r  rho  [t_e  t_r]  [mass_fractions...]

Where:
  - i        : cell index (0-based)
  - inwind   : 0 = in wind, -1 = not in wind (guard cells)
  - r        : radial coordinate [cm]
  - v_r      : radial velocity [cm/s]
  - rho      : mass density [g/cm3]
  - t_e      : electron temperature [K] (optional, written if available)
  - t_r      : radiation temperature [K] (optional, written if available)
  - mass_fractions : mass fractions for each element (if available)

Guard cells:
  - 1 guard cell at inner edge (inwind = -1)
  - 2 guard cells at outer edge (inwind = -1)

Reference: https://sirocco-rt.readthedocs.io/en/latest/wind_models/importing_models.html
"""

import sys
import re
import argparse
import numpy as np
from pathlib import Path


def read_sn_hydro_data(filepath):
    """
    Parse CMFGEN SN_HYDRO_DATA file.

    Returns dict with keys: r_cm, v_cms, rho, T (temperature), all as 1D arrays.
    Also returns inner_boundary_radius and inner_boundary_temperature from the header.
    Also returns mass_fractions as a dict of element names to arrays.
    Grid is inner->outer (r increasing).
    """
    filepath = Path(filepath)
    text = filepath.read_text()
    lines = text.splitlines()

    def _read_block(marker):
        try:
            idx = lines.index(marker)
        except ValueError:
            return None
        vals = []
        for ln in lines[idx + 1:]:
            if not ln.strip():
                break
            for tok in ln.split():
                try:
                    vals.append(float(tok))
                except ValueError:
                    pass
        return np.array(vals) if vals else None

    def _read_header_value(marker):
        """Read a scalar value from header lines like 'Inner boundary radius [cm]: 1.470459E+15'"""
        for ln in lines:
            if marker in ln:
                parts = ln.split(':')
                if len(parts) > 1:
                    try:
                        return float(parts[-1].strip())
                    except ValueError:
                        pass
        return None

    def _read_header_int(marker):
        """Read an integer value from header lines like 'Number of mass fractions: 8'"""
        for ln in lines:
            if marker in ln:
                parts = ln.split(':')
                if len(parts) > 1:
                    try:
                        return int(parts[-1].strip())
                    except ValueError:
                        pass
        return None

    # Read inner boundary conditions from header
    inner_boundary_radius = _read_header_value("Inner boundary radius")
    inner_boundary_temperature = _read_header_value("Inner boundary temperature")
    n_mass_fractions = _read_header_int("Number of mass fractions")

    # SN_HYDRO_DATA is stored outer->inner (R decreasing), so we reverse
    r_1e10 = _read_block("Radius grid (10^10cm)")
    v_kms = _read_block("Velocity (km/s)")
    rho = _read_block("Density (gm/cm^3)")
    T = _read_block("Temperature (10^4 K)")

    if r_1e10 is None or v_kms is None or rho is None:
        raise ValueError(f"Could not parse required blocks from {filepath}")

    # Reverse so array is inner->outer (r increasing)
    r_cm = r_1e10[::-1] * 1.0e10
    v_cms = v_kms[::-1] * 1.0e5
    rho = rho[::-1]
    T_K = T[::-1] * 1.0e4 if T is not None else None

    # If header values not available, use first grid cell values
    if inner_boundary_radius is None:
        inner_boundary_radius = r_cm[0]
    if inner_boundary_temperature is None:
        inner_boundary_temperature = T_K[0] if T_K is not None else 1000.0

    # Read mass fractions if available
    mass_fractions = {}
    if n_mass_fractions is not None and n_mass_fractions > 0:
        # Look for mass fraction blocks (they have names like "HYDR mass fraction", "HELI mass fraction", etc.)
        for ln in lines:
            if "mass fraction" in ln:
                element_name = ln.split()[0]  # Get the element name (e.g., "HYDR")
                mass_frac_block = _read_block(f"{element_name} mass fraction")
                if mass_frac_block is not None:
                    # Reverse to match inner->outer order
                    mass_fractions[element_name] = mass_frac_block[::-1]

    return {
        'r_cm': r_cm,
        'v_cms': v_cms,
        'rho': rho,
        'T_K': T_K,
        'inner_boundary_radius': inner_boundary_radius,
        'inner_boundary_temperature': inner_boundary_temperature,
        'mass_fractions': mass_fractions,
    }


def write_sirocco_import(data, outfile, n_guard_inner=1, n_guard_outer=2):
    """
    Write SIROCCO 1D spherical import file.

    Format matches SIROCCO's expected format:
      i  r  v_r  rho  t_e  t_r  [mass_fractions...]
    
    Guard cells are prepended/appended with rho=0 and zero mass fractions.
    """
    outfile = Path(outfile)
    r = data['r_cm']
    v = data['v_cms']
    rho = data['rho']
    T = data['T_K']
    mass_fractions = data.get('mass_fractions', {})
    n = len(r)

    # Build guard cell radii by extrapolating
    # Inner guard: place at 90% of first wind cell radius (ensures strictly smaller)
    r_inner_guards = [r[0] * 0.9 ** (n_guard_inner - i) for i in range(n_guard_inner)]
    v_inner_guards = [v[0] * 0.9 ** (n_guard_inner - i) for i in range(n_guard_inner)]

    # Outer guard: extend beyond last wind cell with 1% spacing
    dr_outer = max(r[-1] - r[-2], r[-1] * 0.01)
    r_outer_guards = [r[-1] + (i + 1) * dr_outer for i in range(n_guard_outer)]
    v_outer_guards = [v[-1] + (i + 1) * (v[-1] - v[-2]) for i in range(n_guard_outer)]

    has_T = T is not None
    has_mass_fractions = len(mass_fractions) > 0
    element_names = sorted(mass_fractions.keys()) if has_mass_fractions else []
    
    # Check if composition is uniform
    is_uniform_composition = True
    if has_mass_fractions:
        for elem in element_names:
            if len(mass_fractions[elem]) > 1:
                if not np.allclose(mass_fractions[elem], mass_fractions[elem][0], rtol=1e-10):
                    is_uniform_composition = False
                    break
    
    if is_uniform_composition:
        print(f"  Composition is uniform - will write same mass fractions for all cells")
    else:
        print(f"  Composition varies spatially - will write per-cell mass fractions")

    with open(outfile, 'w') as f:
        # Write header (SIROCCO doesn't parse headers, but useful for humans)
        f.write("# SIROCCO 1D spherical import file\n")
        f.write("# Columns: i r v_r rho t_e t_r")
        if has_mass_fractions:
            for elem in element_names:
                f.write(f" {elem}")
        f.write("\n")

        cell_idx = 0

        # Inner guard cells
        for i in range(n_guard_inner):
            r_g = max(r_inner_guards[i], 1.0)
            v_g = v_inner_guards[i]
            # Use first wind cell values for guard cells to avoid numerical issues
            rho_val = max(rho[0], 1e-20)  # Avoid zero density
            if has_T:
                t_e_val = max(T[0], 100.0)  # Avoid zero temperature
                t_r_val = t_e_val * 1.1
            else:
                t_e_val = 1000.0
                t_r_val = 1100.0
            line = f"{cell_idx} {r_g:.6e} {v_g:.6e} {rho_val:.6e} {t_e_val:.6e} {t_r_val:.6e}"
            if has_mass_fractions:
                for elem in element_names:
                    line += f" {mass_fractions[elem][0]:.6e}"
            f.write(line + "\n")
            cell_idx += 1

        # Wind cells
        for i in range(n):
            rho_val = max(rho[i], 0.0)
            if has_T:
                t_e_val = T[i] if T[i] > 0 else 1000.0
                t_r_val = t_e_val * 1.1  # Default t_r = 1.1 * t_e
            else:
                t_e_val = 1000.0
                t_r_val = 1100.0
            line = f"{cell_idx} {r[i]:.6e} {v[i]:.6e} {rho_val:.6e} {t_e_val:.6e} {t_r_val:.6e}"
            if has_mass_fractions:
                if is_uniform_composition:
                    # Write same mass fractions for all cells
                    for elem in element_names:
                        line += f" {mass_fractions[elem][0]:.6e}"
                else:
                    # Write per-cell mass fractions
                    for elem in element_names:
                        mf_val = mass_fractions[elem][i] if i < len(mass_fractions[elem]) else 0.0
                        line += f" {mf_val:.6e}"
            f.write(line + "\n")
            cell_idx += 1

        # Outer guard cells
        for i in range(n_guard_outer):
            r_g = r_outer_guards[i]
            v_g = v_outer_guards[i]
            # Use last wind cell values for guard cells to avoid numerical issues
            rho_val = max(rho[-1], 1e-25)  # Avoid zero density
            if has_T:
                t_e_val = max(T[-1], 100.0)  # Avoid zero temperature
                t_r_val = t_e_val * 1.1
            else:
                t_e_val = 1000.0
                t_r_val = 1100.0
            line = f"{cell_idx} {r_g:.6e} {v_g:.6e} {rho_val:.6e} {t_e_val:.6e} {t_r_val:.6e}"
            if has_mass_fractions:
                for elem in element_names:
                    line += f" {mass_fractions[elem][-1]:.6e}"
            f.write(line + "\n")
            cell_idx += 1

    total_cells = n_guard_inner + n + n_guard_outer
    return {
        'n_wind': n,
        'n_total': total_cells,
        'r_min': r[0],
        'r_max': r[-1],
        'v_min': v[0],
        'v_max': v[-1],
        'n_mass_fractions': len(element_names),
    }


def convert(sn_hydro_path, output_path):
    """Main conversion: SN_HYDRO_DATA -> SIROCCO .import.txt"""
    print(f"Reading SN_HYDRO_DATA: {sn_hydro_path}")
    data = read_sn_hydro_data(sn_hydro_path)
    n = len(data['r_cm'])
    print(f"  Grid: {n} zones")
    print(f"  r: {data['r_cm'][0]:.3e} - {data['r_cm'][-1]:.3e} cm")
    print(f"  v: {data['v_cms'][0]:.3e} - {data['v_cms'][-1]:.3e} cm/s")
    print(f"  rho: {data['rho'].min():.3e} - {data['rho'].max():.3e} g/cm3")
    if data['T_K'] is not None:
        print(f"  T: {data['T_K'].min():.3e} - {data['T_K'].max():.3e} K")
    print(f"  Inner boundary: R={data['inner_boundary_radius']:.3e} cm, T={data['inner_boundary_temperature']:.3e} K")
    if data['mass_fractions']:
        print(f"  Mass fractions: {len(data['mass_fractions'])} elements")
        for elem in sorted(data['mass_fractions'].keys()):
            print(f"    {elem}: {data['mass_fractions'][elem].min():.3e} - {data['mass_fractions'][elem].max():.3e}")

    print(f"Writing SIROCCO import file: {output_path}")
    stats = write_sirocco_import(data, output_path)
    print(f"  Wind cells: {stats['n_wind']}, Total (with guards): {stats['n_total']}")
    if stats['n_mass_fractions'] > 0:
        print(f"  Mass fractions per cell: {stats['n_mass_fractions']}")
    print(f"Done: {output_path}")
    
    # Write boundary conditions file for SIROCCO parameter file
    boundary_file = Path(output_path).parent / "boundary_conditions.txt"
    with open(boundary_file, 'w') as f:
        f.write(f"inner_boundary_radius_cm = {data['inner_boundary_radius']:.6e}\n")
        f.write(f"inner_boundary_temperature_K = {data['inner_boundary_temperature']:.6e}\n")
    print(f"Boundary conditions written: {boundary_file}")
    
    return stats


def main():
    parser = argparse.ArgumentParser(
        description='Convert CMFGEN SN_HYDRO_DATA to SIROCCO 1D spherical import format'
    )
    parser.add_argument('sn_hydro_data', help='Path to CMFGEN SN_HYDRO_DATA file')
    parser.add_argument('output', help='Output path for SIROCCO .import.txt file')
    args = parser.parse_args()

    convert(args.sn_hydro_data, args.output)


if __name__ == '__main__':
    main()
