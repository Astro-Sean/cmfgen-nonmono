#!/usr/bin/env python3
"""
Wrapper script to convert HERACLES snapshot to CMFGEN format using YAML configuration.

This script:
1. Reads YAML configuration
2. Finds HERACLES snapshot closest to target_days
3. Converts snapshot to CMFGEN SN_HYDRO_DATA format
4. Saves output to model-specific directory

Usage:
    python convert_heracles_to_cmfgen_wrapper.py <config_file.yaml>
"""

import sys
import yaml
import argparse
from pathlib import Path
import subprocess
import os


def load_yaml_config(config_file):
    """Load YAML configuration file."""
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    return config


def convert_heracles_to_cmfgen(config):
    """
    Convert HERACLES snapshot to CMFGEN format using YAML configuration.
    
    Parameters:
    -----------
    config : dict
        YAML configuration dictionary
    
    Returns:
    --------
    output_file : Path
        Path to generated CMFGEN input file
    """
    # Extract model info
    model_name = config.get('model_name', 'unknown')
    output = config.get('output', {})
    heracles_output = output.get('heracles_output', '../models/${model_name}/heracles').replace('${model_name}', model_name)
    
    print(f"[DEBUG] model_name: {model_name}")
    print(f"[DEBUG] heracles_output (before replace): {output.get('heracles_output', '../models/${model_name}/heracles')}")
    print(f"[DEBUG] heracles_output (after replace): {heracles_output}")
    
    # CMFGEN conversion output path (from cmfgen_conversion section)
    cmfgen_conversion = config.get('cmfgen_conversion', {})
    cmfgen_output_file = cmfgen_conversion.get('output_file', '../models/${model_name}/cmfgen/SN_HYDRO_DATA').replace('${model_name}', model_name)
    
    print(f"[DEBUG] cmfgen_output_file (before replace): {cmfgen_conversion.get('output_file', '../models/${model_name}/cmfgen/SN_HYDRO_DATA')}")
    print(f"[DEBUG] cmfgen_output_file (after replace): {cmfgen_output_file}")
    
    # Use heracles_dir from cmfgen_conversion section if available, otherwise fall back to heracles_output
    heracles_dir = cmfgen_conversion.get('heracles_dir', heracles_output).replace('${model_name}', model_name)
    
    print(f"[DEBUG] heracles_dir (from cmfgen_conversion): {heracles_dir}")
    
    # Also get the initial profile file path (in the same directory as heracles_dir)
    initial_profile = Path(heracles_dir) / 'input_profile'
    print(f"[DEBUG] initial_profile path: {initial_profile}")
    print(f"[DEBUG] Does initial_profile exist? {initial_profile.exists()}")
    
    # Extract directory from output file
    cmfgen_output = str(Path(cmfgen_output_file).parent)
    
    plots_dir = output.get('plots', '../models/${model_name}/plots').replace('${model_name}', model_name)
    
    print(f"[DEBUG] plots_dir (before replace): {output.get('plots', '../models/${model_name}/plots')}")
    print(f"[DEBUG] plots_dir (after replace): {plots_dir}")
    
    # Extract HERACLES parameters
    heracles = config.get('heracles', {})
    target_days = heracles.get('target_days', 15)
    
    print(f"[DEBUG] target_days: {target_days}")
    
    # Extract CMFGEN conversion parameters
    cmfgen_conversion = config.get('cmfgen_conversion', {})
    tau_max = cmfgen_conversion.get('tau_max', 30)
    tau_min = cmfgen_conversion.get('tau_min', None)
    n_grid = cmfgen_conversion.get('n_grid', 100)
    opacity_model = cmfgen_conversion.get('opacity_model', 'electron_scattering')
    kappa_constant = cmfgen_conversion.get('kappa_constant', 0.34)
    
    # Create directories
    print(f"[DEBUG] Creating directory: {cmfgen_output}")
    Path(cmfgen_output).mkdir(parents=True, exist_ok=True)
    print(f"[DEBUG] Creating directory: {plots_dir}")
    Path(plots_dir).mkdir(parents=True, exist_ok=True)
    
    # Set output file path from YAML config
    output_file = Path(cmfgen_output_file)
    
    # Find closest snapshot
    print(f"Finding HERACLES snapshot closest to {target_days} days...")
    
    # Get HERACLES root from environment or use default
    heracles_root = os.environ.get('HERACLES_ROOT', '/ptmp/sbrennan/heracles')
    print(f"[DEBUG] heracles_root: {heracles_root}")
    print(f"[DEBUG] Current working directory: {os.getcwd()}")
    print(f"[DEBUG] Looking for snapshots in: {heracles_dir}")
    print(f"[DEBUG] Absolute path to heracles_dir: {os.path.abspath(heracles_dir)}")
    print(f"[DEBUG] Does heracles_dir exist? {os.path.exists(heracles_dir)}")
    
    # Call find_closest_snapshot script
    find_snap_script = Path(__file__).parent / 'find_closest_snapshot.py'
    print(f"[DEBUG] find_snap_script: {find_snap_script}")
    print(f"[DEBUG] Does find_snap_script exist? {find_snap_script.exists()}")
    print(f"[DEBUG] About to call find_closest_snapshot...")
    
    print(f"[DEBUG] Running find_closest_snapshot...")
    result = subprocess.run(
        [sys.executable, str(find_snap_script), str(heracles_dir), str(target_days),
         '--heracles-root', heracles_root],
        capture_output=True,
        text=True,
        timeout=300  # 5 minute timeout to prevent hanging
    )
    
    print(f"[DEBUG] find_closest_snapshot returncode: {result.returncode}")
    
    if result.returncode != 0:
        print(f"ERROR: Failed to find snapshot")
        print(result.stderr)
        return None
    
    # Parse snapshot path from output (last line is the actual path, debug lines are [DEBUG])
    lines = result.stdout.strip().split('\n')
    # Find the last line that doesn't start with [DEBUG]
    snapshot_path = None
    for line in reversed(lines):
        if not line.strip().startswith('[DEBUG]'):
            snapshot_path = line.strip()
            break
    
    print(f"[DEBUG] snapshot_path: {snapshot_path}")
    
    if snapshot_path is None:
        print(f"ERROR: Could not parse snapshot path from output")
        print(result.stdout)
        return None
    
    print(f"Found snapshot: {snapshot_path}")
    
    # Call convert_heracles_to_cmfgen script
    print(f"Converting snapshot to CMFGEN format...")
    
    convert_script = Path(__file__).parent / 'convert_heracles_to_cmfgen.py'
    
    # Build command line arguments
    # convert_heracles_to_cmfgen.py expects: heracles_dir target_days output_file [tau_min] [tau_max] [n_grid]
    cmd = [
        sys.executable,
        str(convert_script),
        str(heracles_dir),
        str(target_days),
        str(output_file)
    ]
    
    # Add initial profile path if it exists
    if initial_profile.exists():
        cmd.append(f'--initial-profile={initial_profile}')
        print(f"[DEBUG] Adding initial profile to command: {initial_profile}")
    
    # Add tau_min, tau_max, and n_grid if specified in config
    # IMPORTANT: Order must be: tau_min tau_max n_grid
    # If any are None, we need to pass placeholders or adjust the order
    if 'cmfgen_conversion' in config:
        conv_config = config['cmfgen_conversion']
        tau_min = conv_config.get('tau_min')
        tau_max = conv_config.get('tau_max')
        n_grid = conv_config.get('n_grid')
        
        # Always pass tau_min and tau_max (use None if not specified)
        # This ensures correct argument parsing
        if tau_min is not None and tau_min != 'null':
            cmd.append(str(tau_min))
        else:
            cmd.append('None')  # Placeholder for tau_min
        
        if tau_max is not None and tau_max != 'null':
            cmd.append(str(tau_max))
        else:
            cmd.append('None')  # Placeholder for tau_max
        
        if n_grid is not None:
            cmd.append(str(n_grid))
        adaptive_grid = conv_config.get('adaptive_grid', False)
        if adaptive_grid:
            cmd.append('--adaptive-grid')
        gradient_threshold = conv_config.get('gradient_threshold', 0.5)
        if adaptive_grid:
            cmd.append(f'--gradient-threshold={gradient_threshold}')
        
        # Add sigma smoothing parameters
        smooth_sigma = conv_config.get('smooth_sigma', False)
        if smooth_sigma:
            cmd.append('--smooth-sigma')
            sigma_window = conv_config.get('sigma_window', 9)
            sigma_polyorder = conv_config.get('sigma_polyorder', 3)
            sigma_threshold = conv_config.get('sigma_threshold', 2.0)
            cmd.append(f'--sigma-window={sigma_window}')
            cmd.append(f'--sigma-polyorder={sigma_polyorder}')
            cmd.append(f'--sigma-threshold={sigma_threshold}')
        
        # Add velocity trimming parameter
        velocity_max_kms = conv_config.get('velocity_max_kms')
        if velocity_max_kms is not None:
            cmd.append(f'--velocity-max-kms={velocity_max_kms}')
    
    print(f"[DEBUG] About to call convert_heracles_to_cmfgen...")
    print(f"[DEBUG] Command: {' '.join(cmd)}")
    
    print(f"[DEBUG] Running convert_heracles_to_cmfgen...")
    result = subprocess.run(cmd, timeout=600)  # 10 minute timeout, no capture to avoid pipe deadlock
    
    print(f"[DEBUG] convert_heracles_to_cmfgen returncode: {result.returncode}")
    
    if result.returncode != 0:
        print(f"ERROR: Conversion failed")
        return None
    
    print(f"\nCMFGEN input file created: {output_file}")
    
    return output_file


def main():
    parser = argparse.ArgumentParser(
        description='Convert HERACLES snapshot to CMFGEN format using YAML configuration'
    )
    parser.add_argument('config_file', type=str, help='YAML configuration file')
    
    args = parser.parse_args()
    
    # Load configuration
    config = load_yaml_config(args.config_file)
    
    # Convert snapshot
    output_file = convert_heracles_to_cmfgen(config)
    
    if output_file:
        print(f"\nConversion completed successfully")
        sys.exit(0)
    else:
        print(f"\nConversion failed")
        sys.exit(1)


if __name__ == '__main__':
    main()
