#!/usr/bin/env python3
"""
Generate a grid of STELLA models with varying explosion energies and RSG shell masses.

This script:
1. Takes a base YAML configuration as template
2. Varies explosion_energy and M_ejecta_A over specified ranges
3. Creates individual config files for each grid point
4. Saves them to grid/ directory with descriptive names
5. Can optionally submit the grid using run_model_grid.py
6. Can disable CMFGEN to run STELLA only

Usage:
    python generate_stella_grid.py --base-config configs/X1.yaml
    python generate_stella_grid.py --base-config configs/X1.yaml --submit
    python generate_stella_grid.py --base-config configs/X1.yaml --no-cmfgen --submit
"""

import yaml
import argparse
import numpy as np
from pathlib import Path
import subprocess
import sys


def load_base_config(base_config_path):
    """Load the base YAML configuration."""
    with open(base_config_path, 'r') as f:
        config = yaml.safe_load(f)
    return config


def generate_grid_name(base_name, energy, mass):
    """Generate a descriptive model name for the grid point."""
    # Format: X1_E{energy}M{mass}
    # Energy in 1e50 erg, mass in Msun
    energy_str = f"{energy:.1f}".replace('.', 'p')
    mass_str = f"{mass:.1f}".replace('.', 'p')
    return f"{base_name}_E{energy_str}M{mass_str}"


def create_grid_config(base_config, energy, mass, output_dir, run_cmfgen=False):
    """Create a modified config for a specific grid point."""
    config = base_config.copy()
    
    # Modify the specific parameters
    config['stella']['explosion_energy'] = energy
    config['ejecta']['M_ejecta_A'] = mass
    
    # Update model name and description
    model_name = generate_grid_name(base_config['model_name'], energy, mass)
    config['model_name'] = model_name
    config['description'] = f"Grid model: E={energy}e50 erg, M_RSG={mass} Msun"
    
    # Disable CMFGEN for grid models (STELLA only)
    config['run_cmfgen'] = run_cmfgen
    
    # Set output directory to grid/ instead of models/
    # Use absolute paths to ensure everything goes to grid/
    grid_dir = output_dir.resolve()
    if 'output' not in config:
        config['output'] = {}
    config['output']['model_dir'] = str(grid_dir / model_name)
    config['output']['base_dir'] = str(grid_dir)
    config['output']['stella_output'] = str(grid_dir / model_name / 'stella')
    config['output']['plots_dir'] = str(grid_dir / model_name / 'plots')
    
    # Also set cmfgen base_dir to avoid models/ prefix
    if 'cmfgen' not in config:
        config['cmfgen'] = {}
    config['cmfgen']['base_dir'] = str(grid_dir.parent / 'base')
    
    # Save to grid directory
    output_path = output_dir / f"{model_name}.yaml"
    with open(output_path, 'w') as f:
        yaml.dump(config, f, default_flow_style=False)
    
    return output_path


def generate_grid(base_config_path, energies, masses, output_dir='grid', run_cmfgen=False):
    """Generate the full grid of configurations."""
    base_config = load_base_config(base_config_path)
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    grid_configs = []
    
    print(f"Generating STELLA parameter grid:")
    print(f"  Base config: {base_config_path}")
    print(f"  Energy range: {min(energies)}-{max(energies)} x 1e50 erg")
    print(f"  Mass range: {min(masses)}-{max(masses)} Msun")
    print(f"  Output directory: {output_dir}")
    print(f"  Run CMFGEN: {run_cmfgen}")
    print()
    
    for energy in energies:
        for mass in masses:
            config_path = create_grid_config(base_config, energy, mass, output_path, run_cmfgen)
            grid_configs.append(config_path)
            print(f"  Created: {config_path.name}")
    
    print(f"\nTotal grid points: {len(grid_configs)}")
    return grid_configs


def submit_grid(grid_configs, debug=False, debug_days=1, no_cmfgen=False):
    """Submit the grid using run_model_grid.py."""
    print(f"\nSubmitting grid with {len(grid_configs)} models...")
    
    # Build command
    cmd = ['python', 'scripts/run_model_grid.py']
    
    if debug:
        cmd.extend(['--debug', f'--debug-days={int(debug_days)}'])
    
    if no_cmfgen:
        cmd.append('--no-cmfgen')
    
    cmd.extend([str(c) for c in grid_configs])
    
    print(f"Command: {' '.join(cmd)}")
    
    # Submit from project directory (scripts will handle grid directory)
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    print(result.stdout)
    if result.stderr:
        print("STDERR:", result.stderr)
    
    return result.returncode == 0


def main():
    parser = argparse.ArgumentParser(description='Generate STELLA parameter grid')
    parser.add_argument('--base-config', type=str, default='configs/X1.yaml',
                        help='Base YAML configuration file (default: configs/X1.yaml)')
    parser.add_argument('--output-dir', type=str, default='grid',
                        help='Output directory for grid configs (default: grid)')
    parser.add_argument('--energies', type=str, default='0.5,1.0,2.0,3.0',
                        help='Comma-separated explosion energies in 1e50 erg (default: 0.5,1.0,2.0,3.0)')
    parser.add_argument('--masses', type=str, default='0.1,0.5,1.0,2.0,3.0,4.0',
                        help='Comma-separated RSG shell masses in Msun (default: 0.1,0.5,1.0,2.0,3.0,4.0)')
    parser.add_argument('--no-cmfgen', action='store_true',
                        help='Disable CMFGEN (STELLA only)')
    parser.add_argument('--submit', action='store_true',
                        help='Submit the grid after generation')
    parser.add_argument('--debug', action='store_true',
                        help='Run in debug mode (1 day evolution)')
    parser.add_argument('--debug-days', type=float, default=1.0,
                        help='Debug mode evolution time in days (default: 1.0)')
    
    args = parser.parse_args()
    
    # Parse parameter ranges
    energies = [float(x) for x in args.energies.split(',')]
    masses = [float(x) for x in args.masses.split(',')]
    
    # Generate grid
    grid_configs = generate_grid(args.base_config, energies, masses, args.output_dir, not args.no_cmfgen)
    
    # Optionally submit
    if args.submit:
        success = submit_grid(grid_configs, args.debug, args.debug_days, args.no_cmfgen)
        if success:
            print("\n✓ Grid submitted successfully")
        else:
            print("\n✗ Grid submission failed")
            sys.exit(1)
    else:
        print(f"\nGrid configs saved to {args.output_dir}/")
        print("To submit, run:")
        cmfgen_flag = "--no-cmfgen" if args.no_cmfgen else ""
        print(f"  python scripts/generate_stella_grid.py --base-config {args.base_config} --submit {cmfgen_flag}")


if __name__ == '__main__':
    main()