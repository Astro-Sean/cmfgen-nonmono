#!/usr/bin/env python3
"""
Generate single-stage HERACLES namelist from YAML configuration.

Usage:
    python generate_heracles_namelist.py <config_file.yaml> [output_namelist]
"""

import sys
import yaml
import argparse
from pathlib import Path


def load_yaml_config(config_file):
    """Load YAML configuration file."""
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    return config


def generate_namelist_from_config(config, output_file=None):
    """
    Generate HERACLES namelist from YAML configuration.
    
    Parameters:
    -----------
    config : dict
        YAML configuration dictionary
    output_file : str or Path, optional
        Output namelist file path. If None, uses default.
    
    Returns:
    --------
    output_file : Path
        Path to generated namelist file
    """
    # Extract HERACLES parameters
    heracles = config.get('heracles', {})
    mode = heracles.get('mode', 'hydro')
    simulation_days = float(heracles.get('simulation_days', 90))
    Lbox_cm = float(heracles.get('Lbox_cm', 6.0e16))
    r_inner_cm = float(heracles.get('r_inner_cm', 1.0e7))
    r_transition_cm = float(heracles.get('r_transition_cm', 1.0e12))
    nx_inner_linear = int(heracles.get('nx_inner_linear', 3000))
    nx_hybrid_blend = int(heracles.get('nx_hybrid_blend', 0))
    hybrid_sigmoid_k = float(heracles.get('hybrid_sigmoid_k', 10))
    dt_max = float(heracles.get('dt_max', 300))
    nout = int(heracles.get('nout', 100))
    
    # Extract ejecta parameters for grid
    ejecta = config.get('ejecta', {})
    n_cells = ejecta.get('n_cells', 10000)
    
    # Extract model info
    model_name = config.get('model_name', 'unknown')
    output = config.get('output', {})
    base_dir = output.get('base_dir', './models').replace('${model_name}', model_name)
    heracles_output = output.get('heracles_output', 'heracles').replace('${model_name}', model_name)
    
    print(f"[DEBUG] base_dir: {base_dir}")
    print(f"[DEBUG] heracles_output (from config): {heracles_output}")
    
    # Resolve to absolute path from project directory
    project_dir = Path(__file__).parent.parent
    print(f"[DEBUG] project_dir: {project_dir}")
    
    # If heracles_output is a relative path without 'models', assume it's relative to model_dir
    if not str(heracles_output).startswith('models') and not str(heracles_output).startswith('./models') and not str(heracles_output).startswith('../'):
        # It's a relative path from the model directory
        model_dir = project_dir / base_dir / model_name
        heracles_output = model_dir / heracles_output
    else:
        # It's a relative path from project directory
        heracles_output = (project_dir / heracles_output).absolute()
    
    print(f"[DEBUG] heracles_output (resolved): {heracles_output}")
    
    # Create directory
    Path(heracles_output).mkdir(parents=True, exist_ok=True)
    
    # Set output file path
    if output_file is None:
        output_file = Path(heracles_output) / 'liste.nml'
    else:
        output_file = Path(output_file)
    
    # Calculate simulation time in seconds
    tend = simulation_days * 86400.0
    
    # Set solver based on mode
    if mode == 'hydrorad':
        solver = 'exact'
        ray = '.true.'
        diffu = '.true.'
    else:  # hydro
        solver = 'exact'
        ray = '.false.'
        diffu = '.false.'
    
    # Generate namelist content
    namelist_content = f"""! Single-stage HERACLES simulation for model: {model_name}
! Mode: {mode}
! Simulation time: {simulation_days} days ({tend:.1e} s)
! Generated from YAML configuration

&PARAM
reprise        = .false.
i_restart      = 0
tend           = {tend:.6e}
dt_max         = {dt_max}
Tcpu_max       = {tend * 1.1:.6e}
kmax           = 100000000
Nimp           = 100000
nfx            = 6
nout           = {nout}
slope_type     = 3
Solver         = '{solver}'
RAY            = {ray}
DIFFU          = {diffu}
format_sortie  = 'bin'
geometrie      = 'rs'
&END

&PARAM_GRILLE
ndim           = 1
nx_glob(1)     = {n_cells}
nx_glob(2)     = 1
nx_glob(3)     = 1
icl_glob(1)    = 1
icl_glob(2)    = 1
icl_glob(3)    = 1
icl_glob(4)    = 1
icl_glob(5)    = 1
icl_glob(6)    = 1
icl_ray_glob(1) = 1
icl_ray_glob(2) = 1
icl_ray_glob(3) = 1
icl_ray_glob(4) = 1
icl_ray_glob(5) = 1
icl_ray_glob(6) = 1
icl_dif_glob(1) = 1
icl_dif_glob(2) = 1
icl_dif_glob(3) = 1
icl_dif_glob(4) = 1
icl_dif_glob(5) = 1
icl_dif_glob(6) = 1
nsx            = 1
nsy            = 1
nsz            = 1
&END
"""

    # Add radiation and diffusion blocks for hydrorad mode
    if mode == 'hydrorad':
        namelist_content += f"""
&PARAM_RAY
alpha             = 1.0
varrel_ray        = 0.20
nitemax_ray       = 100000
matrix_inv_solver = 'GS'
precond_ray       = .true.
alpha_explicite   = 3.0
cal_valp          = .true.
valp_min          = 0.00
slope_type_ray    = 1
&END

&PARAM_DIF
varrel_dif = 0.05
&END

&DIR
dir_valp='/ptmp/sbrennan/heracles/heracles-2.0/src/rad_transfer/m1'
&END
"""

    namelist_content += f"""

&CONDINIT
Lbox(1)        = {Lbox_cm:.6e}
Lbox(2)        = 1.0
Lbox(3)        = 1.0
&END

&PROFILEINIT
 use_profile           = .true.
 profile_file          = 'input_profile'
 radial_stretch_power  = 1
 use_hybrid_radial_mesh = .false.
 r_inner_face           = {r_inner_cm:.6e}
 r_transition           = {r_transition_cm:.6e}
 nx_inner_linear        = {nx_inner_linear}
 nx_hybrid_blend        = {nx_hybrid_blend}
 hybrid_sigmoid_k       = {hybrid_sigmoid_k:.6e}

&END

&BOMBINIT
use_bomb       = .false.
e_bomb_ergs    = 0.0
n_cells_bomb   = 0
&END

&JOB
&END

&REPR
&END

&OUTPUT
&END
"""
    
    # Write namelist
    with open(output_file, 'w') as f:
        f.write(namelist_content)
    
    print(f"Generated HERACLES namelist: {output_file}")
    print(f"  Mode: {mode}")
    print(f"  Simulation time: {simulation_days} days")
    print(f"  Grid: {n_cells} cells, Lbox = {Lbox_cm:.2e} cm")
    
    return output_file


def main():
    parser = argparse.ArgumentParser(
        description='Generate single-stage HERACLES namelist from YAML configuration'
    )
    parser.add_argument('config_file', type=str, help='YAML configuration file')
    parser.add_argument('output_namelist', type=str, nargs='?', help='Output namelist file (optional)')
    
    args = parser.parse_args()
    
    # Load configuration
    config = load_yaml_config(args.config_file)
    
    # Generate namelist
    output_file = generate_namelist_from_config(config, args.output_namelist)
    
    print(f"\nNamelist generated successfully: {output_file}")


if __name__ == '__main__':
    main()
