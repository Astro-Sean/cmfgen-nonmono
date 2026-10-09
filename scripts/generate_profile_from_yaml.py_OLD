#!/usr/bin/env python3
"""
Generate [DEPRECATED] HERACLES input profile from YAML configuration with diagnostic plots.

HERACLES is deprecated and does not fit the science case (interacting Type IIn SNe).
The project now uses the STELLA -> CMFGEN pipeline exclusively. STELLA is the active
backend; HERACLES options are retained for backward compatibility.

Usage:
    python generate_profile_from_yaml.py <config_file.yaml> [output_profile]
"""

import sys
import yaml
import argparse
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

# [DEPRECATED] write_heracles_profile is for the deprecated HERACLES backend;
# project uses STELLA -> CMFGEN pipeline. STELLA is the active backend.
from create_embedded_hhe_shell import create_embedded_hhe_profile, write_heracles_profile, plot_profile


def load_yaml_config(config_file):
    """Load YAML configuration file."""
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    return config


def generate_profile_from_config(config, output_file=None, hydro_backend='heracles', debug_mode=False):
    """
    Generate [DEPRECATED] HERACLES profile from YAML configuration.
    The project uses STELLA -> CMFGEN pipeline; HERACLES is deprecated.
    
    Parameters:
    -----------
    config : dict
        YAML configuration dictionary
    output_file : str or Path, optional
        Output profile file path. If None, uses default.
    hydro_backend : str, optional
        Hydro backend: 'heracles' [DEPRECATED] or 'stella' (active). Determines output directory.
    debug_mode : bool, optional
        If True, override days_ejecta_B to 0 for debug runs (start from explosion).
    
    Returns:
    --------
    output_file : Path
        Path to generated profile file
    """
    # Extract ejecta parameters
    ejecta = config.get('ejecta', {})
    M_ejecta_A = ejecta.get('M_ejecta_A', 9.8)
    Mdot_wind = float(ejecta.get('Mdot_wind', 1e-3))
    M_ejecta_B = ejecta.get('M_ejecta_B', 0.1)
    E_kin_ejecta_A = float(ejecta.get('E_kin_ejecta_A', 1e51))
    E_kin_ejecta_B = float(ejecta.get('E_kin_ejecta_B', 0))
    days_ejecta_A = ejecta.get('days_ejecta_A', 11.6)
    days_ejecta_B = ejecta.get('days_ejecta_B', 400)
    fx_H_A = ejecta.get('fx_H_A', 0.70)
    fx_He_A = ejecta.get('fx_He_A', 0.28)
    fx_C_A = ejecta.get('fx_C_A', 0.003)
    fx_O_A = ejecta.get('fx_O_A', 0.009)
    fx_Si_A = ejecta.get('fx_Si_A', 0.0007)
    fx_Ca_A = ejecta.get('fx_Ca_A', 0.0002)
    fx_Fe_A = ejecta.get('fx_Fe_A', 0.0013)
    fx_H_wind = ejecta.get('fx_H_wind', 0.70)
    fx_He_wind = ejecta.get('fx_He_wind', 0.28)
    fx_H_B = ejecta.get('fx_H_B', 0.01)
    fx_He_B = ejecta.get('fx_He_B', 0.40)
    fx_C_B = ejecta.get('fx_C_B', 0.03)
    fx_O_B = ejecta.get('fx_O_B', 0.15)
    fx_Si_B = ejecta.get('fx_Si_B', 0.05)
    fx_Ca_B = ejecta.get('fx_Ca_B', 0.02)
    fx_Fe_B = ejecta.get('fx_Fe_B', 0.03)
    fx_Ni_B = ejecta.get('fx_Ni_B', 0.02)
    V_ejecta_B_max_kms = ejecta.get('V_ejecta_B_max_kms', 10000.0)
    N_rho_ejecta_B = ejecta.get('N_rho_ejecta_B', 7)
    n_cells = ejecta.get('n_cells', 10000)
    
    # Extract transition layer parameters
    transition_width_fraction = ejecta.get('transition_width_fraction', 0.2)
    transition_power_law_exponent = ejecta.get('transition_power_law_exponent', 0.5)
    sigmoid_k = ejecta.get('sigmoid_k', 20)
    # A/B specific transition smoothing (H/He -> metal)
    transition_AB_width_fraction = ejecta.get('transition_AB_width_fraction', 0.10)
    transition_AB_sigmoid_k = ejecta.get('transition_AB_sigmoid_k', 10)
    
    # Thermal bomb mode: stationary inner H/He shell
    stationary_ejecta_A = ejecta.get('stationary_ejecta_A', False)
    
    # New approach: give inner ejecta initial velocity for collision
    V_ejecta_A_initial_kms = ejecta.get('V_ejecta_A_initial_kms', 0.0)
    # Shell-CSM gap fraction: wind region between fast shell and slow CSM.
    # Only used when V_ejecta_A_initial_kms > 0. Prevents velocity reversal.
    # If collision_time_days > 0, gap width is set from collision time instead.
    shell_csm_gap_fraction = ejecta.get('shell_csm_gap_fraction', 0.02)
    collision_time_days = ejecta.get('collision_time_days', 0.0)

    # RSG progenitor shell + wind parameters
    R_star_cm  = float(ejecta['R_star_cm']) if 'R_star_cm' in ejecta else None
    n_shell    = ejecta.get('n_shell', 2.0)
    wind_velocity_type = ejecta.get('wind_velocity_type', 'beta_law')
    V_inf_kms  = ejecta.get('V_inf_kms', 100.0)
    beta_wind  = ejecta.get('beta_wind', 5.0)
    wind_time_days = ejecta.get('wind_time_days', 400.0)
    R_CSM_MAX  = ejecta.get('R_CSM_MAX', None)
    
    # Extract temperature profile parameters
    temperature = config.get('temperature', {})
    temp_profile_type = temperature.get('temp_profile_type', 'adiabatic')
    N_T_ejecta = temperature.get('N_T_ejecta', 0.6)
    T_inner_max = temperature.get('T_inner_max', 100000)
    T_csm = temperature.get('T_csm', 6000.0)
    
    # Extract output directories
    output = config.get('output', {})
    model_name = config.get('model_name', 'unknown')
    base_dir = output.get('base_dir', './models').replace('${model_name}', model_name)
    model_dir = output.get('model_dir', './models/${model_name}').replace('${model_name}', model_name)
    
    # Use different output directory based on hydro_backend
    # [DEPRECATED] HERACLES backend; project uses STELLA -> CMFGEN pipeline
    if hydro_backend.lower() == 'stella':
        heracles_output = output.get('stella_output', 'stella').replace('${model_name}', model_name)
    else:
        heracles_output = output.get('heracles_output', 'heracles').replace('${model_name}', model_name)
    plots_dir = output.get('plots', 'plots').replace('${model_name}', model_name)
    
    print(f"[DEBUG] base_dir: {base_dir}")
    print(f"[DEBUG] model_dir: {model_dir}")
    print(f"[DEBUG] heracles_output (from config): {heracles_output}")
    print(f"[DEBUG] plots_dir (from config): {plots_dir}")
    
    # Resolve to absolute path from project directory
    project_dir = Path(__file__).parent.parent
    print(f"[DEBUG] project_dir: {project_dir}")
    
    # Helper function to resolve relative paths
    def resolve_path(path):
        if not str(path).startswith('models') and not str(path).startswith('./models') and not str(path).startswith('../'):
            # It's a relative path from the model directory
            return project_dir / base_dir / model_name / path
        else:
            # It's a relative path from project directory
            return (project_dir / path).absolute()
    
    model_dir = resolve_path(model_dir)
    heracles_output = resolve_path(heracles_output)
    plots_dir = resolve_path(plots_dir)
    
    print(f"[DEBUG] model_dir (resolved): {model_dir}")
    print(f"[DEBUG] heracles_output (resolved): {heracles_output}")
    print(f"[DEBUG] plots_dir (resolved): {plots_dir}")
    
    # Create directories
    Path(heracles_output).mkdir(parents=True, exist_ok=True)
    Path(plots_dir).mkdir(parents=True, exist_ok=True)
    
    # Set output file path
    if output_file is None:
        output_file = Path(heracles_output) / 'input_profile'
    else:
        output_file = Path(output_file)
    
    # Create profile
    print(f"Generating HERACLES profile for model: {model_name}")
    print(f"  M_ejecta_A = {M_ejecta_A} Msun")
    print(f"  Mdot_wind = {Mdot_wind} Msun/yr")
    print(f"  M_ejecta_B = {M_ejecta_B} Msun")
    print(f"  n_cells = {n_cells}")
    
    profile = create_embedded_hhe_profile(
        days_ejecta_A=days_ejecta_A,
        days_ejecta_B=days_ejecta_B,
        n_cells=n_cells,
        fx_Fe_ejecta=fx_Fe_B,
        fx_Si_ejecta=fx_Si_B,
        fx_He_ejecta=fx_He_B,
        fx_H_ejecta=fx_H_B,
        M_ejecta_A=M_ejecta_A,
        Mdot_wind=Mdot_wind,
        M_ejecta_B=M_ejecta_B,
        E_kin_ejecta_A=E_kin_ejecta_A,
        E_kin_ejecta_B=E_kin_ejecta_B,
        fx_H_A=fx_H_A,
        fx_He_A=fx_He_A,
        fx_C_A=fx_C_A,
        fx_O_A=fx_O_A,
        fx_Si_A=fx_Si_A,
        fx_Ca_A=fx_Ca_A,
        fx_Fe_A=fx_Fe_A,
        fx_H_wind=fx_H_wind,
        fx_He_wind=fx_He_wind,
        fx_H_B=fx_H_B,
        fx_He_B=fx_He_B,
        fx_C_B=fx_C_B,
        fx_O_B=fx_O_B,
        fx_Si_B=fx_Si_B,
        fx_Ca_B=fx_Ca_B,
        fx_Fe_B=fx_Fe_B,
        fx_Ni_B=fx_Ni_B,
        temp_profile_type=temp_profile_type,
        N_T_ejecta=N_T_ejecta,
        T_inner_max=T_inner_max,
        T_csm=T_csm,
        transition_width_fraction=transition_width_fraction,
        transition_power_law_exponent=transition_power_law_exponent,
        sigmoid_k=sigmoid_k,
        transition_AB_width_fraction=transition_AB_width_fraction,
        transition_AB_sigmoid_k=transition_AB_sigmoid_k,
        stationary_ejecta_A=stationary_ejecta_A,
        R_star_cm=R_star_cm,
        n_shell=n_shell,
        wind_velocity_type=wind_velocity_type,
        V_inf_kms=V_inf_kms,
        beta_wind=beta_wind,
        wind_time_days=wind_time_days,
        R_CSM_MAX=R_CSM_MAX,
        V_ejecta_B_max_kms=V_ejecta_B_max_kms,
        N_rho_ejecta_B=N_rho_ejecta_B,
        V_ejecta_A_initial_kms=V_ejecta_A_initial_kms,
        shell_csm_gap_fraction=shell_csm_gap_fraction,
        collision_time_days=collision_time_days,
        days_ejecta_B_ref=ejecta.get('days_ejecta_B_ref', 1.0),
        T_inner_max_B=ejecta.get('T_inner_max_B', 2.0e6),
        m_ejecta_B=ejecta.get('m_ejecta_B', 2.0),
        v_H_cut_ejecta_B_kms=ejecta.get('v_H_cut_ejecta_B_kms', 8000.0),
        fx_H_B_outer=ejecta.get('fx_H_B_outer', 0.60),
        fx_He_B_outer=ejecta.get('fx_He_B_outer', 0.38),
        fx_C_B_outer=ejecta.get('fx_C_B_outer', 0.002),
        fx_O_B_outer=ejecta.get('fx_O_B_outer', 0.01),
        fx_Si_B_outer=ejecta.get('fx_Si_B_outer', 0.003),
        fx_Ca_B_outer=ejecta.get('fx_Ca_B_outer', 0.001),
        fx_Fe_B_outer=ejecta.get('fx_Fe_B_outer', 0.001),
        gamma_ejecta_B=ejecta.get('gamma_ejecta_B', 4.0/3.0),
        outer_wind_enabled=ejecta.get('outer_wind_enabled', True),
        T_interface_shell=temperature.get('T_interface_shell', 2000.0),
        buffer_zone_width_rstar=ejecta.get('buffer_zone_width_rstar', 0.0),
        buffer_zone_density_fraction=ejecta.get('buffer_zone_density_fraction', 0.01),
    )
    
    # Write profile
    # [DEPRECATED] HERACLES profile writer; STELLA -> CMFGEN is the active pipeline
    write_heracles_profile(profile, output_file)
    
    # Generate diagnostic plots using the same style as project plots
    print(f"\nGenerating diagnostic plots in {plots_dir}")
    try:
        plot_profile(profile, output_dir=Path(plots_dir))
    except ImportError:
        print("  Warning: matplotlib not available, skipping plots")
    except Exception as e:
        print(f"  Warning: Could not generate plots: {e}")
    
    return output_file


def main():
    parser = argparse.ArgumentParser(
        description='Generate [DEPRECATED] HERACLES input profile from YAML configuration '
                    '(project uses STELLA -> CMFGEN pipeline)'
    )
    parser.add_argument('config_file', type=str, help='YAML configuration file')
    parser.add_argument('output_profile', type=str, nargs='?', help='Output profile file (optional)')
    parser.add_argument('--hydro-backend', type=str, default='heracles',
                        help='Hydro backend: heracles [DEPRECATED] or stella (active). '
                             'Default: heracles [DEPRECATED]. Project uses STELLA -> CMFGEN pipeline.')
    
    args = parser.parse_args()
    
    # Load configuration
    config = load_yaml_config(args.config_file)
    
    # Generate profile
    output_file = generate_profile_from_config(config, args.output_profile, args.hydro_backend)
    
    print(f"\nProfile generated successfully: {output_file}")


if __name__ == '__main__':
    main()
