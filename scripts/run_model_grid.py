#!/usr/bin/env python3
"""
Master workflow script for automated supernova model grid.

This script:
1. Scans configs/ for YAML files
2. For each config:
   - Generates HERACLES input profile
   - Generates HERACLES namelist
   - Submits HERACLES job (named: heracles_{model_name})
   - Submits CMFGEN conversion job with dependency (named: convert_{model_name})
   - Submits CMFGEN run job with dependency (named: cmfgen_{model_name})
3. Returns job IDs for monitoring

Usage:
    python run_model_grid.py [config1.yaml config2.yaml ...]
    python run_model_grid.py --all  # Process all configs in configs/
"""

import sys
import yaml
import argparse
import subprocess
from pathlib import Path
import os
import re


def load_yaml_config(config_file):
    """Load YAML configuration file."""
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    return config


def substitute_template(template, variables):
    """Substitute variables in template string."""
    result = template
    for key, value in variables.items():
        placeholder = "{" + key + "}"
        result = result.replace(placeholder, str(value))
    return result


def generate_profile(config, scripts_dir, project_dir, config_path, hydro_backend='heracles'):
    """Generate HERACLES input profile from YAML config."""
    model_name = config.get('model_name', 'unknown')
    print(f"\n{'='*60}")
    print(f"Generating profile for model: {model_name}")
    print(f"{'='*60}")
    
    script = scripts_dir / 'generate_profile_from_yaml.py'
    config_file = config_path
    
    result = subprocess.run(
        [sys.executable, str(script), str(config_file), '--hydro-backend', hydro_backend],
        capture_output=True,
        text=True
    )
    
    if result.returncode != 0:
        print(f"ERROR: Profile generation failed")
        print(result.stdout)
        print(result.stderr)
        return False
    
    print(result.stdout)
    return True


def generate_namelist(config, scripts_dir, project_dir, config_path, debug=False):
    """Generate HERACLES namelist from YAML config."""
    model_name = config.get('model_name', 'unknown')
    print(f"\n{'='*60}")
    print(f"Generating HERACLES namelist for model: {model_name}")
    print(f"{'='*60}")
    
    script = scripts_dir / 'generate_heracles_namelist.py'
    config_file = config_path
    
    result = subprocess.run(
        [sys.executable, str(script), str(config_file)],
        capture_output=True,
        text=True
    )
    
    if result.returncode != 0:
        print(f"ERROR: Namelist generation failed")
        print(result.stdout)
        print(result.stderr)
        return False
    
    print(result.stdout)
    return True


def check_existing_jobs(model_name):
    """Check if there are existing Slurm jobs for this model."""
    result = subprocess.run(
        ['squeue', '-u', os.environ.get('USER', os.environ.get('LOGNAME'))],
        capture_output=True,
        text=True
    )
    
    if result.returncode != 0:
        return []
    
    jobs = []
    for line in result.stdout.strip().split('\n')[1:]:  # Skip header
        parts = line.split()
        if len(parts) >= 3:
            job_id = parts[0]  # Column 0 is the JOBID
            job_name = parts[2]  # Column 2 is the job NAME
            # Match exact model name with underscore prefix/suffix to avoid partial matches
            # e.g., "X" matches "heracles_X", "cmf_X_init", "cmfgen_X_INIT" but not "heracles_X1"
            # Also match CMFGEN jobs that contain the model name in their job name
            if f'_{model_name}_' in job_name or job_name.endswith(f'_{model_name}') or model_name in job_name:
                jobs.append(job_id)
    
    return jobs


def prepare_cmfgen_directory(config, scripts_dir, project_dir):
    """Prepare CMFGEN directory by copying base files."""
    model_name = config.get('model_name', 'unknown')
    print(f"\n{'='*60}")
    print(f"Preparing CMFGEN directory for model: {model_name}")
    print(f"{'='*60}")
    
    # Extract paths
    output = config.get('output', {})
    model_dir = output.get('model_dir', '../models/${model_name}').replace('${model_name}', model_name)
    model_dir = (project_dir / model_dir).absolute()
    
    cmfgen = config.get('cmfgen', {})
    base_dir = cmfgen.get('base_dir', '../base')
    base_dir = (project_dir / base_dir).absolute()

    # model_type selects the CMFGEN executable used in batch.sh
    # Supported types and their exe names (relative to $cmfdist/exe/):
    #   supernova (default) -> cmfgen_dev.exe
    #   lte                 -> main_lte.exe
    MODEL_TYPE_EXE = {
        'supernova': 'cmfgen_dev.exe',
        'lte':       'main_lte.exe',
    }
    model_type = cmfgen.get('model_type', 'supernova')
    exe_name = MODEL_TYPE_EXE.get(model_type)
    if exe_name is None:
        print(f"ERROR: Unknown model_type '{model_type}'. Supported: {list(MODEL_TYPE_EXE)}")
        return False

    cmfgen_dir = model_dir / 'cmfgen'
    
    # Copy base directory to cmfgen directory
    if not base_dir.exists():
        print(f"ERROR: Base directory not found: {base_dir}")
        return False
    
    if cmfgen_dir.exists():
        print(f"Removing existing cmfgen directory: {cmfgen_dir}")
        import shutil
        shutil.rmtree(cmfgen_dir)
    
    print(f"Copying base directory from {base_dir} to {cmfgen_dir}...")
    import shutil
    shutil.copytree(base_dir, cmfgen_dir)

    # Patch batch.sh to use the selected executable
    batch_sh = cmfgen_dir / 'batch.sh'
    if batch_sh.exists():
        text = batch_sh.read_text()
        import re
        patched = re.sub(
            r'(?m)^(setenv CMFGEN_PROG\s+.*?/)[a-z_]+\.exe\s*$',
            r'\g<1>' + exe_name,
            text,
        )
        if patched != text:
            batch_sh.write_text(patched)
            print(f"  Patched batch.sh: CMFGEN_PROG -> {exe_name}  (model_type={model_type})")
        else:
            print(f"  WARNING: Could not patch CMFGEN_PROG in batch.sh (line not found); using as-is")
    else:
        print(f"  WARNING: batch.sh not found in {cmfgen_dir}; skipping exe patch")

    print(f"CMFGEN directory prepared: {cmfgen_dir}  [model_type={model_type}, exe={exe_name}]")
    return True


def prepare_sedona_directory(config, scripts_dir, project_dir):
    """Prepare SEDONA directory by copying base files and creating param.lua."""
    model_name = config.get('model_name', 'unknown')
    print(f"\n{'='*60}")
    print(f"Preparing SEDONA directory for model: {model_name}")
    print(f"{'='*60}")
    
    # Extract paths
    output = config.get('output', {})
    model_dir = output.get('model_dir', '../models/${model_name}').replace('${model_name}', model_name)
    model_dir = (project_dir / model_dir).absolute()
    
    sedona_dir = model_dir / 'sedona'
    sedona_dir.mkdir(parents=True, exist_ok=True)
    
    # Copy sedona_base files
    sedona_base = project_dir / 'sedona_base'
    if sedona_base.exists():
        print(f"Copying SEDONA base files from {sedona_base}...")
        import shutil
        for f in sedona_base.iterdir():
            if f.is_file():
                dest = sedona_dir / f.name
                shutil.copy2(f, dest)
                print(f"  Copied {f.name}")
    else:
        print(f"WARNING: sedona_base directory not found at {sedona_base}")
    
    # Convert CMFGEN SN_HYDRO_DATA to SEDONA .mod format
    sn_hydro = model_dir / 'cmfgen' / 'SN_HYDRO_DATA'
    mod_file = sedona_dir / f'{model_name}.mod'
    
    if sn_hydro.exists():
        print(f"Converting CMFGEN SN_HYDRO_DATA to SEDONA .mod format...")
        # Use the conversion script from sedona
        sedona_convert = project_dir.parent / 'sedona' / 'convert_cmfgen_to_sedona_v2.py'
        if sedona_convert.exists():
            # Remove existing .mod file to ensure fresh conversion
            if mod_file.exists():
                print(f"  Removing existing {mod_file.name} to ensure fresh conversion")
                mod_file.unlink()
            
            result = subprocess.run(
                [sys.executable, str(sedona_convert), str(sn_hydro), str(mod_file)],
                capture_output=True,
                text=True
            )
            if result.returncode == 0:
                print(f"  Successfully created {model_name}.mod")
                print(result.stdout)
            else:
                print(f"  ERROR: Conversion failed")
                print(result.stdout)
                print(result.stderr)
                return False
        else:
            print(f"  ERROR: Conversion script not found at {sedona_convert}")
            return False
    else:
        print(f"  ERROR: SN_HYDRO_DATA not found at {sn_hydro}")
        print(f"    SEDONA requires CMFGEN SN_HYDRO_DATA to create the .mod file")
        return False

    # Analyze .mod file: read innermost zone radius and detect zero-abundance species
    mod_species_z = []
    zero_abundance_z = set()
    r_first_zone = None
    if mod_file.exists():
        with open(mod_file) as f:
            mod_lines = f.readlines()
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
        if zero_abundance_z:
            print(f"  WARNING: Zero-abundance species in .mod: {sorted(zero_abundance_z)}")
            print(f"    These will be removed from opacity_atoms_in_nlte")
    
    # Get snapshot days from cmfgen_conversion config (fallback only)
    cmfgen_conv = config.get('cmfgen_conversion', {})
    snapshot_days_config = cmfgen_conv.get('snapshot_days', 40)
    
    # Read SN_HYDRO_DATA for inner boundary info (tau-trimmed model)
    # and time since energy injection
    sn_hydro = model_dir / 'cmfgen' / 'SN_HYDRO_DATA'
    r_inner = 6.0e14      # default fallback
    T_inner = 5000.0      # default fallback (K)
    luminosity = 8.32e42  # default fallback (erg/s)
    snapshot_days = snapshot_days_config  # fallback
    
    if sn_hydro.exists():
        try:
            import numpy as np
            with open(sn_hydro) as f:
                lines = f.readlines()
            
            nd = 500  # default number of zones
            radius_data = []
            temp_data = []
            luminosity_lsun = None  # Will read from header
            
            # Parse the file to find sections and header values
            section = None
            for i, line in enumerate(lines):
                stripped = line.strip()
                
                # Find ND
                if 'Number of data points' in stripped:
                    parts = stripped.split()
                    for j, p in enumerate(parts):
                        if p == 'points:':
                            nd = int(parts[j+1])
                            break
                
                # Find luminosity from header (Energy Lum [lsun])
                if 'Energy Lum [lsun]' in stripped:
                    parts = stripped.split()
                    for j, p in enumerate(parts):
                        if p == '[lsun]:':
                            luminosity_lsun = float(parts[j+1])
                            break
                
                # Find time since explosion (energy injection) from header
                if 'Time(days) since explosion' in stripped:
                    parts = stripped.split()
                    for j, p in enumerate(parts):
                        if p == 'explosion:':
                            snapshot_days = float(parts[j+1])
                            break
                
                # Identify sections (exact header matches)
                if stripped.startswith('Radius grid'):
                    section = 'radius'
                    continue
                elif stripped.startswith('Velocity'):
                    section = 'velocity'
                    continue
                elif stripped.startswith('Temperature'):
                    section = 'temperature'
                    continue
                elif stripped.startswith('Sigma'):
                    section = 'sigma'
                    continue
                elif stripped == 'Density (gm/cm^3)':
                    section = 'density'
                    continue
                elif stripped.startswith('Electron density'):
                    section = 'electron_density'
                    continue
                elif stripped.startswith('Atom density'):
                    section = 'atom_density'
                    continue
                elif stripped.startswith('Kappa'):
                    section = 'kappa'
                    continue
                elif stripped == '':
                    section = None
                    continue
                
                # Collect data based on section
                if section == 'radius' and stripped and not stripped.startswith('Number'):
                    values = [float(x) for x in stripped.split()]
                    radius_data.extend(values)
                elif section == 'temperature' and stripped:
                    values = [float(x) for x in stripped.split()]
                    temp_data.extend(values)
            
            # Get innermost zone (last entry = innermost in CMFGEN format)
            # This is the boundary of the tau-trimmed region
            if radius_data and temp_data:
                # Innermost is the last zone (deepest, smallest radius)
                r_inner_10e10 = radius_data[nd-1]  # 10^10 cm units
                T_inner_10e4 = temp_data[nd-1]     # 10^4 K units (from stellar sim)
                
                r_inner = r_inner_10e10 * 1e10      # convert to cm
                T_inner = T_inner_10e4 * 1e4        # convert to K (physical T from sim)
                
                # Get luminosity from SN_HYDRO_DATA header, fallback to default
                if luminosity_lsun is not None:
                    luminosity = luminosity_lsun * 3.828e33  # convert L_sun to erg/s
                else:
                    # Calculate from R and T if header L not available
                    sigma_sb = 5.6704e-5
                    luminosity = 4 * np.pi * r_inner**2 * sigma_sb * T_inner**4
                    print(f"  WARNING: Could not read luminosity from SN_HYDRO_DATA, calculating from R,T")

                # Try reading new inner-boundary header lines (written on previous runs)
                r_inner_from_header = None
                T_inner_from_header = None
                L_inner_from_header = None
                for line in lines:
                    stripped = line.strip()
                    if stripped.startswith('Inner boundary radius [cm]:'):
                        try:
                            r_inner_from_header = float(stripped.split()[-1])
                        except (ValueError, IndexError):
                            pass
                    elif stripped.startswith('Inner boundary temperature [K]:'):
                        try:
                            T_inner_from_header = float(stripped.split()[-1])
                        except (ValueError, IndexError):
                            pass
                    elif stripped.startswith('Inner boundary luminosity [erg/s]:'):
                        try:
                            L_inner_from_header = float(stripped.split()[-1])
                        except (ValueError, IndexError):
                            pass

                if r_inner_from_header is not None:
                    r_inner = r_inner_from_header
                if T_inner_from_header is not None:
                    T_inner = T_inner_from_header
                if L_inner_from_header is not None:
                    luminosity = L_inner_from_header

                print(f"  Tau-trimmed model (SN_HYDRO_DATA):")
                print(f"    Time since energy injection: {snapshot_days:.2f} days")
                print(f"    Inner boundary (zone {nd}): R={r_inner:.6e} cm, T={T_inner:.0f} K")
                print(f"    Luminosity: {luminosity:.3e} erg/s")
                print(f"    (Values from stellar simulation at tau_max boundary)")
        except Exception as e:
            print(f"  WARNING: Could not parse SN_HYDRO_DATA: {e}")
            print(f"    Using defaults: R={r_inner:.3e} cm, T={T_inner:.0f} K")
    else:
        print(f"  WARNING: SN_HYDRO_DATA not found at {sn_hydro}")
        print(f"    Using defaults: R={r_inner:.3e} cm, T={T_inner:.0f} K")
    
    # Create param.lua from template
    param_template = sedona_dir / 'param.lua'
    if param_template.exists():
        text = param_template.read_text()
        # Replace placeholders (use @ delimiters to avoid Lua table syntax conflicts)
        text = text.replace('@model_name@', model_name)
        text = text.replace('@snapshot_days@', str(snapshot_days))

        # Ensure core_radius does not exceed first zone radius (avoid zones inside core)
        if r_first_zone is not None and r_inner >= r_first_zone:
            core_radius = r_first_zone * 0.9999
            print(f"  Adjusted core_radius from {r_inner:.6e} to {core_radius:.6e} cm")
            print(f"    (first zone radius = {r_first_zone:.6e} cm)")
        else:
            core_radius = r_inner
        text = text.replace('@r_inner@', f"{core_radius:.6e}")
        text = text.replace('@temperature@', f"{T_inner:.2e}")
        text = text.replace('@luminosity@', f"{luminosity:.3e}")

        # Filter species from NLTE list that are missing from .mod or have zero abundance
        missing_from_mod = set()
        if mod_species_z:
            nlte_match = re.search(r'opacity_atoms_in_nlte\s*=\s*\{([^}]*)\}', text)
            if nlte_match:
                all_z = [int(x.strip()) for x in nlte_match.group(1).split(',') if x.strip()]
                missing_from_mod = set(all_z) - set(mod_species_z)

        species_to_remove = zero_abundance_z | missing_from_mod

        # Temperature-based NLTE filtering: heavy species cause singular matrices at low T
        NLTE_TEMP_THRESHOLD = 5000.0  # K
        COLD_NLTE_ATOMS = {14, 26}    # Si, Fe (and Ni, Co if ever added)
        if T_inner < NLTE_TEMP_THRESHOLD:
            cold_atoms_in_list = COLD_NLTE_ATOMS & set(all_z if 'all_z' in dir() else [])
            if not cold_atoms_in_list and mod_species_z:
                # all_z may not be defined yet if no missing_from_mod; rebuild it
                nlte_match = re.search(r'opacity_atoms_in_nlte\s*=\s*\{([^}]*)\}', text)
                if nlte_match:
                    all_z = [int(x.strip()) for x in nlte_match.group(1).split(',') if x.strip()]
                    cold_atoms_in_list = COLD_NLTE_ATOMS & set(all_z)
            if cold_atoms_in_list:
                species_to_remove = species_to_remove | cold_atoms_in_list
                print(f"  Low T inner boundary ({T_inner:.0f} K < {NLTE_TEMP_THRESHOLD:.0f} K):")
                print(f"    Removing heavy NLTE species {sorted(cold_atoms_in_list)} to avoid singular matrix")

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

        param_template.write_text(text)
        print(f"  Generated param.lua:")
        print(f"    Time since energy injection: {snapshot_days:.2f} days")
        print(f"    Inner boundary: R={core_radius:.6e} cm, T={T_inner:.0f} K, L={luminosity:.3e} erg/s")
    else:
        print(f"  WARNING: param.lua template not found in sedona_base")
    
    print(f"SEDONA directory prepared: {sedona_dir}")
    return True


def submit_heracles_job(config, scripts_dir, project_dir, debug=False, config_file=None):
    """Submit HERACLES Slurm job."""
    model_name = config.get('model_name', 'unknown')
    print(f"\n{'='*60}")
    print(f"Submitting HERACLES job for model: {model_name}")
    print(f"{'='*60}")
    
    # Extract paths
    output = config.get('output', {})
    model_dir = output.get('model_dir', '../models/${model_name}').replace('${model_name}', model_name)
    model_dir = (project_dir / model_dir).absolute()
    
    heracles = config.get('heracles', {})
    mode = heracles.get('mode', 'hydro')
    simulation_days = 2 if debug else heracles.get('simulation_days', 90)
    Lbox_cm = heracles.get('Lbox_cm', 6.0e16)
    r_inner_cm = heracles.get('r_inner_cm', 1.0e7)
    r_transition_cm = heracles.get('r_transition_cm', 1.0e12)
    nx_inner_linear = heracles.get('nx_inner_linear', 3000)
    nx_hybrid_blend = heracles.get('nx_hybrid_blend', 0)
    hybrid_sigmoid_k = heracles.get('hybrid_sigmoid_k', 10)
    n_cells = 10000
    dt_max = heracles.get('dt_max', 300)
    
    # Resource settings
    partition = 'interactive' if debug else heracles.get('partition', 'small')
    cpus = heracles.get('cpus', 2)
    memory = heracles.get('memory', '16G')
    # Enforce partition time limits: interactive = 1:59:00 max, general/small = 1-00:00:00 (1 day)
    if debug:
        time = '01:59:00'
    else:
        time = heracles.get('time', '23:59:00')
        # Override to 1-day limit if partition is general or small and time is > 1 day
        if partition in ('general', 'small') and time != '1-00:00:00':
            time = '1-00:00:00'
    
    if debug:
        print(f"DEBUG MODE: Using partition={partition}, time={time}, simulation_days={simulation_days}")
    
    # HERACLES root
    heracles_root = os.environ.get('HERACLES_ROOT', '/ptmp/sbrennan/heracles')
    # Force correct path if environment variable points to heracles-2.0
    if heracles_root.endswith('/heracles-2.0'):
        heracles_root = heracles_root.replace('/heracles-2.0', '')
    
    # Load template
    template_file = scripts_dir / 'heracles_job_template.slurm'
    with open(template_file, 'r') as f:
        template = f.read()
    
    # Use debug config file if provided, otherwise use original config
    if config_file is None:
        config_file = project_dir / 'configs' / f'model_{model_name}.yaml'
    
    # Substitute variables (convert to proper types for formatting)
    variables = {
        'MODEL_NAME': model_name,
        'MODEL_DIR': str(model_dir),
        'HERACLES_DIR': str(model_dir / 'heracles'),
        'HERACLES_ROOT': heracles_root,
        'MODE': mode,
        'SIMULATION_DAYS': simulation_days,
        'PARTITION': partition,
        'CPUS': cpus,
        'MEMORY': memory,
        'TIME': time,
        'USE_HYBRID_MESH': 'true',
        'LBOX_CM': f"{float(Lbox_cm):.6e}",
        'N_CELLS': n_cells,
        'R_TRANSITION_CM': f"{float(r_transition_cm):.6e}",
        'NX_INNER_LINEAR': nx_inner_linear,
        'HYBRID_SIGMOID_K': hybrid_sigmoid_k,
        'R_INNER_CM': f"{float(r_inner_cm):.6e}",
        'PROJECT_DIR': str(project_dir),
        'CONFIG_FILE': str(config_file)
    }
    
    job_script = substitute_template(template, variables)
    
    # Write job script
    heracles_dir = model_dir / 'heracles'
    heracles_dir.mkdir(parents=True, exist_ok=True)
    job_file = heracles_dir / f'heracles_{model_name}.slurm'
    
    with open(job_file, 'w') as f:
        f.write(job_script)
    
    # Submit job
    result = subprocess.run(
        ['sbatch', str(job_file)],
        capture_output=True,
        text=True
    )
    
    if result.returncode != 0:
        print(f"ERROR: HERACLES job submission failed")
        print(result.stdout)
        print(result.stderr)
        return None
    
    # Parse job ID from output (format: "Submitted batch job 12345")
    job_id = result.stdout.strip().split()[-1]
    print(f"HERACLES job submitted: {job_id}")
    
    return job_id


def submit_conversion_job(config, scripts_dir, project_dir, heracles_job_id, debug=False):
    """Submit CMFGEN conversion job with dependency on HERACLES.

    NOTE: This function is currently unused. The HERACLES workflow embeds the conversion
    step at the end of heracles_job_template.slurm, and the STELLA workflow embeds it
    inside stella_job_template.slurm. Kept for potential future use as a standalone
    conversion job.
    """
    model_name = config.get('model_name', 'unknown')
    print(f"\n{'='*60}")
    print(f"Submitting CMFGEN conversion job for model: {model_name}")
    print(f"{'='*60}")
    
    # Extract paths
    output = config.get('output', {})
    model_dir = output.get('model_dir', '../models/${model_name}').replace('${model_name}', model_name)
    model_dir = (project_dir / model_dir).absolute()
    
    # Resource settings (use same as HERACLES for conversion)
    heracles = config.get('heracles', {})
    partition = 'interactive' if debug else heracles.get('partition', 'small')
    cpus = heracles.get('cpus', 2)
    memory = heracles.get('memory', '16G')
    # Enforce partition time limits: interactive = 1:59:00 max, general/small = 1-00:00:00 (1 day)
    if debug:
        time = '01:59:00'
    else:
        time = '02:00:00'
        # Override to 1-day limit if partition is general or small
        if partition in ('general', 'small'):
            time = '1-00:00:00'
    
    # HERACLES root
    heracles_root = os.environ.get('HERACLES_ROOT', '/ptmp/sbrennan/heracles')
    # Force correct path if environment variable points to heracles-2.0
    if heracles_root.endswith('/heracles-2.0'):
        heracles_root = heracles_root.replace('/heracles-2.0', '')
    
    # Load template
    template_file = scripts_dir / 'conversion_job_template.slurm'
    with open(template_file, 'r') as f:
        template = f.read()
    
    # Substitute variables
    variables = {
        'MODEL_NAME': model_name,
        'MODEL_DIR': str(model_dir),
        'PROJECT_DIR': str(project_dir),
        'HERACLES_ROOT': heracles_root,
        'HERACLES_JOB_ID': heracles_job_id,
        'PARTITION': partition,
        'CPUS': cpus,
        'MEMORY': memory,
        'TIME': time
    }
    
    job_script = substitute_template(template, variables)
    
    # Write job script
    cmfgen_dir = model_dir / 'cmfgen'
    cmfgen_dir.mkdir(parents=True, exist_ok=True)
    job_file = cmfgen_dir / f'convert_{model_name}.slurm'
    
    with open(job_file, 'w') as f:
        f.write(job_script)
    
    # Submit job
    result = subprocess.run(
        ['sbatch', str(job_file)],
        capture_output=True,
        text=True
    )
    
    if result.returncode != 0:
        print(f"ERROR: Conversion job submission failed")
        print(result.stdout)
        print(result.stderr)
        return None
    
    # Parse job ID
    job_id = result.stdout.strip().split()[-1]
    
    print(f"Conversion job submitted: {job_id}")
    return job_id


def submit_cmfgen_jobs(config, scripts_dir, project_dir, dependency_job_id=None, debug=False, config_file=None):
    """Submit CMFGEN jobs (init, final, obs) using submit_cmfgen_jobs.sh script."""
    model_name = config.get('model_name', 'unknown')
    print(f"\n{'='*60}")
    print(f"Submitting CMFGEN jobs for model: {model_name}")
    print(f"{'='*60}")
    
    # Extract paths
    output = config.get('output', {})
    model_dir = output.get('model_dir', '../models/${model_name}').replace('${model_name}', model_name)
    model_dir = (project_dir / model_dir).absolute()
    cmfgen_dir = model_dir / 'cmfgen'
    
    # Read init_num_its and single_final from config
    # Note: fix_temperature is handled internally by submit_cmfgen.sh (FIX_T=T by default)
    cmfgen_cfg = config.get('cmfgen', {})
    init_num_its = cmfgen_cfg.get('init_num_its', 3)
    single_final = cmfgen_cfg.get('single_final', True)
    
    # Submit script
    submit_script = scripts_dir / 'submit_cmfgen_jobs.sh'
    
    # Make script executable
    os.chmod(submit_script, 0o755)
    
    # Build command
    dep_id = dependency_job_id if dependency_job_id else ""
    # config_file should always be set by process_config, but provide fallback for safety
    config_file_arg = str(config_file) if config_file else ""
    # Enforce partition time limits: interactive = 1:59:00 max, general = 1-00:00:00 (1 day)
    cmfgen_partition = 'interactive' if debug else 'general'
    cmfgen_time = '01:59:00' if debug else '1-00:00:00'  # general partition = 1 day max
    cmd = [str(submit_script), model_name, str(model_dir), str(cmfgen_dir), dep_id,
           str(init_num_its),
           cmfgen_partition, cmfgen_time,
           '1' if single_final else '0',
           config_file_arg]
    if debug:
        cmd.insert(1, '--debug')
    
    # Call the script
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        cwd=project_dir
    )
    
    if result.returncode != 0:
        print(f"ERROR: CMFGEN jobs submission failed")
        print(result.stdout)
        print(result.stderr)
        return None
    
    # Parse job IDs from output (format: "INIT:FINAL:OBS")
    # In single-final mode, INIT is "(none)"
    job_ids_str = result.stdout.strip().split('\n')[-1]
    init_job, final_job, obs_job = job_ids_str.split(':')
    init_job = None if init_job.strip() == '(none)' else init_job.strip()
    
    print(f"CMFGEN jobs submitted:")
    print(f"  INIT:  {init_job if init_job else '(none - single-final mode)'}")
    print(f"  FINAL: {final_job}")
    print(f"  OBS:   {obs_job}")
    
    return {
        'init': init_job,
        'final': final_job,
        'obs': obs_job
    }


def submit_sedona_job(config, scripts_dir, project_dir, dependency_job_id=None, debug=False):
    """Submit SEDONA job to produce spectrum from existing hydro data."""
    model_name = config.get('model_name', 'unknown')
    print(f"\n{'='*60}")
    print(f"Submitting SEDONA job for model: {model_name}")
    print(f"{'='*60}")
    
    # Extract paths
    output_cfg = config.get('output', {})
    model_dir = output_cfg.get('model_dir', '../models/${model_name}').replace('${model_name}', model_name)
    model_dir = (project_dir / model_dir).absolute()
    sedona_dir = model_dir / 'sedona'
    
    # Get SEDONA settings from config
    sedona_cfg = config.get('sedona', {})
    heracles_cfg = config.get('heracles', {})
    
    # Resource settings
    partition = 'interactive' if debug else sedona_cfg.get('partition', heracles_cfg.get('partition', 'general'))
    cpus = sedona_cfg.get('cpus', 4)
    memory = sedona_cfg.get('memory', '16G')
    time_limit = '01:59:00' if debug else sedona_cfg.get('time', '24:00:00')
    # Enforce 1-day limit for general or small partition
    if partition in ('general', 'small') and not debug:
        time_limit = '1-00:00:00'
    
    # Check for required files (warn only - .mod will be created by STELLA/HERACLES conversion)
    mod_file = sedona_dir / f'{model_name}.mod'
    if not mod_file.exists():
        # Try to find any .mod file
        mod_files = list(sedona_dir.glob('*.mod'))
        if mod_files:
            mod_file = mod_files[0]
            print(f"  Found alternative .mod file: {mod_file.name}")
        else:
            print(f"  WARNING: No .mod file found in {sedona_dir}")
            print(f"    The .mod file should be created by the hydro conversion step.")
            print(f"    SEDONA job will be submitted with dependency - it will fail if .mod is missing at runtime.")
    
    # Create SEDONA slurm script from template
    job_file = sedona_dir / f'sedona_{model_name}.slurm'
    
    # Create sedona directory if it doesn't exist
    sedona_dir.mkdir(parents=True, exist_ok=True)
    
    # Load template from sedona_base
    template_file = project_dir / 'sedona_base' / 'sedona_job_template.slurm'
    if not template_file.exists():
        print(f"ERROR: SEDONA job template not found at {template_file}")
        return None
    
    with open(template_file, 'r') as f:
        template = f.read()
    
    # Substitute variables
    variables = {
        'MODEL_NAME': model_name,
        'TIME_LIMIT': time_limit,
        'PARTITION': partition,
        'CPUS': cpus,
        'MEMORY': memory,
        'WORK_DIR': str(sedona_dir),
        'CMFGEN_DATA': str(model_dir / 'cmfgen' / 'SN_HYDRO_DATA'),
        'MOD_FILE': f'{model_name}.mod'
    }
    
    job_script = substitute_template(template, variables)
    
    with open(job_file, 'w') as f:
        f.write(job_script)
    
    # Build dependency string if provided
    dep_arg = f"--dependency=afterok:{dependency_job_id}" if dependency_job_id else ""
    
    # Submit job
    cmd = ['sbatch']
    if dep_arg:
        cmd.append(dep_arg)
    cmd.append(str(job_file))
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"ERROR: SEDONA job submission failed")
        print(result.stdout)
        print(result.stderr)
        return None
    
    job_id = result.stdout.strip().split()[-1]
    print(f"SEDONA job submitted: {job_id}")
    return job_id




def prepare_sirocco_directory(config, scripts_dir, project_dir):
    """Prepare SIROCCO directory by copying base files and creating the .pf parameter file."""
    model_name = config.get('model_name', 'unknown')
    print(f"\n{'='*60}")
    print(f"Preparing SIROCCO directory for model: {model_name}")
    print(f"{'='*60}")

    output = config.get('output', {})
    model_dir = output.get('model_dir', '../models/${model_name}').replace('${model_name}', model_name)
    model_dir = (project_dir / model_dir).absolute()

    sirocco_dir = model_dir / 'sirocco'
    sirocco_dir.mkdir(parents=True, exist_ok=True)

    # Copy sirocco_base files
    sirocco_base = project_dir / 'sirocco_base'
    if sirocco_base.exists():
        print(f"Copying SIROCCO base files from {sirocco_base}...")
        import shutil
        for f in sirocco_base.iterdir():
            if f.is_file():
                dest = sirocco_dir / f.name
                shutil.copy2(f, dest)
                print(f"  Copied {f.name}")
    else:
        print(f"WARNING: sirocco_base directory not found at {sirocco_base}")

    # Convert SN_HYDRO_DATA -> SIROCCO import file
    sn_hydro = model_dir / 'cmfgen' / 'SN_HYDRO_DATA'
    import_file = sirocco_dir / f'{model_name}.import.txt'

    n_wind = 0
    r_inner_cm = 1.0e14
    t_inner_k = 10000.0
    central_mass_msol = 1.4

    if sn_hydro.exists():
        print(f"Converting CMFGEN SN_HYDRO_DATA to SIROCCO import format...")
        convert_script = scripts_dir / 'convert_cmfgen_to_sirocco.py'
        if convert_script.exists():
            # Remove existing import file to ensure fresh conversion
            if import_file.exists():
                print(f"  Removing existing {import_file.name} to ensure fresh conversion")
                import_file.unlink()
            
            result = subprocess.run(
                [sys.executable, str(convert_script), str(sn_hydro), str(import_file)],
                capture_output=True, text=True
            )
            if result.returncode == 0:
                print(f"  Successfully created {model_name}.import.txt")
                print(result.stdout)
            else:
                print(f"  ERROR: SIROCCO conversion failed")
                print(result.stdout)
                print(result.stderr)
                return False
        else:
            print(f"  ERROR: Conversion script not found at {convert_script}")
            return False

        # Read boundary conditions from conversion script output
        boundary_file = sirocco_dir / 'boundary_conditions.txt'
        if boundary_file.exists():
            print(f"  Reading boundary conditions from {boundary_file.name}...")
            for line in boundary_file.read_text().splitlines():
                if 'inner_boundary_radius_cm' in line:
                    r_inner_cm = float(line.split('=')[1].strip())
                elif 'inner_boundary_temperature_K' in line:
                    t_inner_k = float(line.split('=')[1].strip())
            print(f"    Using inner boundary: R={r_inner_cm:.3e} cm, T={t_inner_k:.0f} K")

        # Read key scalars from import file to fill .pf template
        if import_file.exists():
            import numpy as np
            lines = import_file.read_text().splitlines()
            data_lines = [l for l in lines if not l.startswith('#') and l.strip()]
            # Count total data lines (all cells in domain 0 for 1D spherical)
            wind_lines = [l for l in data_lines if len(l.split()) >= 2 and l.split()[1] == '0']
            if wind_lines:
                n_wind = len(wind_lines)
                print(f"  Detected {n_wind} wind cells from import file")
        else:
            print(f"  WARNING: SN_HYDRO_DATA not found at {sn_hydro}")
            print(f"    SIROCCO import file will be created at runtime by the job script.")

    # Patch .pf template
    pf_template = sirocco_dir / 'sirocco.pf'
    pf_dest = sirocco_dir / f'{model_name}.pf'
    if pf_template.exists():
        text = pf_template.read_text()
        replacements = {
            '@model_name@': model_name,
            '@r_inner_cm@': f'{r_inner_cm:.6e}',
            '@t_inner_k@': f'{t_inner_k:.1f}',
            '@central_mass_msol@': f'{central_mass_msol:.2f}',
            '@n_wind@': str(n_wind + 3) if n_wind > 0 else '503',  # +3 guard cells
        }
        for k, v in replacements.items():
            text = text.replace(k, v)
        pf_dest.write_text(text)
        print(f"  Wrote {pf_dest.name} (n_wind={n_wind}, r_inner={r_inner_cm:.3e} cm, T={t_inner_k:.0f} K)")
        pf_template.unlink()
    else:
        print(f"  WARNING: sirocco.pf template not found")

    print(f"SIROCCO directory prepared: {sirocco_dir}")
    return True


def submit_sirocco_job(config, scripts_dir, project_dir, dependency_job_id=None, debug=False):
    """Submit SIROCCO job to produce spectrum from existing hydro data."""
    model_name = config.get('model_name', 'unknown')
    print(f"\n{'='*60}")
    print(f"Submitting SIROCCO job for model: {model_name}")
    print(f"{'='*60}")

    output_cfg = config.get('output', {})
    model_dir = output_cfg.get('model_dir', '../models/${model_name}').replace('${model_name}', model_name)
    model_dir = (project_dir / model_dir).absolute()
    sirocco_dir = model_dir / 'sirocco'

    sirocco_cfg = config.get('sirocco', {})
    heracles_cfg = config.get('heracles', {})

    partition = 'interactive' if debug else sirocco_cfg.get('partition', heracles_cfg.get('partition', 'general'))
    cpus = sirocco_cfg.get('cpus', 4)
    memory = sirocco_cfg.get('memory', '16G')
    time_limit = '01:59:00' if debug else sirocco_cfg.get('time', '24:00:00')
    # Enforce 1-day limit for general or small partition
    if partition in ('general', 'small') and not debug:
        time_limit = '1-00:00:00'

    import_file = sirocco_dir / f'{model_name}.import.txt'
    if not import_file.exists():
        print(f"  WARNING: No import file found at {import_file}")
        print(f"    SIROCCO job will be submitted with dependency - it will fail if import file is missing at runtime.")

    job_file = sirocco_dir / f'sirocco_{model_name}.slurm'
    sirocco_dir.mkdir(parents=True, exist_ok=True)

    template_file = project_dir / 'sirocco_base' / 'sirocco_job_template.slurm'
    if not template_file.exists():
        print(f"ERROR: SIROCCO job template not found at {template_file}")
        return None

    with open(template_file, 'r') as f:
        template = f.read()

    variables = {
        'MODEL_NAME': model_name,
        'TIME_LIMIT': time_limit,
        'PARTITION': partition,
        'CPUS': str(cpus),
        'MEMORY': memory,
        'WORK_DIR': str(sirocco_dir),
        'CMFGEN_DATA': str(model_dir / 'cmfgen' / 'SN_HYDRO_DATA'),
        'IMPORT_FILE': str(import_file),
    }

    job_script = substitute_template(template, variables)

    with open(job_file, 'w') as f:
        f.write(job_script)

    dep_arg = f"--dependency=afterok:{dependency_job_id}" if dependency_job_id else ""

    cmd = ['sbatch']
    if dep_arg:
        cmd.append(dep_arg)
    cmd.append(str(job_file))

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"ERROR: SIROCCO job submission failed")
        print(result.stdout)
        print(result.stderr)
        return None

    job_id = result.stdout.strip().split()[-1]
    print(f"SIROCCO job submitted: {job_id}")
    return job_id


def should_use_multistage(target_days):
    """Determine if multi-stage simulation should be used for long runs."""
    return False  # Disabled - always use single-stage execution


def get_multistage_plan(target_days):
    """Get multi-stage simulation plan for long runs."""
    if target_days <= 30:
        return None
    
    stages = []
    if target_days <= 40:
        # 40-day specific: 3 stages
        stages = [
            {"start": 0, "end": 15, "name": "SBO", "nstmax": 900000, "temp_floor": "3.D-02"},
            {"start": 15, "end": 30, "name": "transition", "nstmax": 1500000, "temp_floor": "5.D-02"},
            {"start": 30, "end": target_days, "name": "nebular", "nstmax": 2000000, "temp_floor": "1.D-01"}
        ]
    else:
        # General long runs: 4 stages
        stage_length = target_days / 4
        for i in range(4):
            start = i * stage_length
            end = (i + 1) * stage_length
            nstmax = 900000 + i * 300000
            temp_floor = f"{3 + i * 2}.D-02"
            stages.append({
                "start": start, "end": end, "name": f"stage_{i+1}",
                "nstmax": nstmax, "temp_floor": temp_floor
            })
    
    return stages


def submit_stella_job(config, scripts_dir, project_dir, debug=False, config_file=None):
    """Submit STELLA Slurm job (runs STELLA + converts swd -> SN_HYDRO_DATA)."""
    model_name = config.get('model_name', 'unknown')
    
    output = config.get('output', {})
    model_dir = output.get('model_dir', '../models/${model_name}').replace('${model_name}', model_name)
    model_dir = (project_dir / model_dir).resolve()

    stella_cfg = config.get('stella', {})
    heracles_cfg = config.get('heracles', {})

    # Support arbitrary simulation days with adaptive stability measures
    stella_target = stella_cfg.get('target_days', None)
    heracles_target = heracles_cfg.get('target_days', 15)
    target_days   = stella_target if stella_target is not None else heracles_target
    
    nzon          = stella_cfg.get('nzon', 200)
    nzon_max      = stella_cfg.get('nzon_max', nzon)  # Allow AMR headroom if > nzon
    kappa         = stella_cfg.get('kappa', 0.34)
    tau_max       = stella_cfg.get('tau_max', 30.0)
    tau_min       = stella_cfg.get('tau_min', 0)
    rho_min       = stella_cfg.get('rho_min') or 0.0
    
    # Use explicit snapshot_interval from config if provided; otherwise compute adaptively
    if stella_cfg.get('snapshot_interval') is not None:
        snapshot_interval = float(stella_cfg['snapshot_interval'])
    else:
        # Adaptive output strategy based on simulation length
        if target_days <= 15:
            # Short simulations: standard 8 outputs
            snapshot_interval = target_days / 8.0
        elif target_days <= 30:
            # Medium simulations: 12 outputs
            snapshot_interval = target_days / 12.0
        else:
            # Long simulations: 20 outputs max, more frequent early on
            snapshot_interval = target_days / 20.0
    
    # Ensure minimum interval of 0.5 days for numerical stability
    snapshot_interval = max(snapshot_interval, 0.5)
    r_min         = stella_cfg.get('r_min', 0.0)
    explosion_energy = stella_cfg.get('explosion_energy', 0.0)
    explosion_mass   = stella_cfg.get('explosion_mass', 0.0)
    explosion_time   = stella_cfg.get('explosion_time', 1.0)
    explosion_delay  = stella_cfg.get('explosion_delay', 0.0)
    initial_kinetic_energy = stella_cfg.get('initial_kinetic_energy', 0.0)
    
    # Detect RSG-like models (extended stationary envelope) for robust defaults
    ejecta_cfg = config.get('ejecta', {})
    is_rsg = (ejecta_cfg.get('stationary_ejecta_A', False) and
              float(ejecta_cfg.get('R_star_cm', 0)) > 1e13)

    # Adaptive numerical parameters for long simulations
    if target_days <= 15:
        nstmax = 900000  # Standard for short runs
    elif target_days <= 30:
        nstmax = 1500000  # Increased for medium runs
    else:
        nstmax = 2000000  # Maximum for long runs
    
    # Override with user-specified value if provided
    nstmax = stella_cfg.get('nstmax', nstmax)
    disable_scattering = stella_cfg.get('disable_scattering', False)

    # RSG-robust defaults: looser accuracy, capped timestep, disable convection (causes instability)
    eps_default = 0.003
    hmax_default = 86400.0 if is_rsg else 5e5  # 1 day vs ~5.8 days
    enable_convection_default = False  # Disabled due to numerical instability in innermost zone

    eps = stella_cfg.get('eps', eps_default)
    hmax = stella_cfg.get('hmax', hmax_default)
    enable_convection = stella_cfg.get('enable_convection', enable_convection_default)
    fit_tau = stella_cfg.get('fit_tau', 15.0)
    tau_tol = stella_cfg.get('tau_tol', 1.3)
    abs_fraction = stella_cfg.get('abs_fraction', 1e-2)
    givdtl = stella_cfg.get('givdtl', True)
    custom_output_times = stella_cfg.get('custom_output_times', None)
    outer_fraction = stella_cfg.get('outer_fraction', 0.3)
    transition_radius = stella_cfg.get('transition_radius', 1e13)
    
    print(f"\n{'='*60}")
    print(f"Submitting STELLA job for model: {model_name}")
    print(f"{'='*60}")
    
    # Show simulation parameters
    print(f"Simulation parameters:")
    print(f"  Target days: {target_days}")
    if custom_output_times:
        print(f"  Custom output times: {custom_output_times}")
    else:
        print(f"  Snapshot interval: {snapshot_interval:.2f} days")
    print(f"  NSTMAX: {nstmax}")
    print(f"  Scattering: {'Disabled' if disable_scattering else 'Enabled'}")
    if is_rsg:
        print(f"  [RSG mode] EPS={eps}, HMAX={hmax:.0f}s, CONV={'T' if enable_convection else 'F'}")
    if should_use_multistage(target_days):
        stages = get_multistage_plan(target_days)
        print(f"  Multi-stage simulation plan:")
        for i, stage in enumerate(stages):
            print(f"    Stage {i+1}: {stage['start']:.1f}-{stage['end']:.1f} days ({stage['name']})")
            print(f"      NSTMAX: {stage['nstmax']}, Temp floor: {stage['temp_floor']}")
    
    cmfgen_conv_cfg  = config.get('cmfgen_conversion', {})
    cmfgen_snapshot_days = cmfgen_conv_cfg.get('snapshot_days', target_days)
    enable_grid_thinning = cmfgen_conv_cfg.get('enable_grid_thinning', False)
    grid_thinning_threshold = cmfgen_conv_cfg.get('grid_thinning_threshold', 0.005)
    # Adaptive grid parameters
    enable_adaptive_grid = cmfgen_conv_cfg.get('enable_adaptive_grid', False)
    base_grid_target = cmfgen_conv_cfg.get('base_grid_target', 100)
    max_radius_jump = cmfgen_conv_cfg.get('max_radius_jump', 0.05)
    max_temp_jump = cmfgen_conv_cfg.get('max_temp_jump', 0.10)
    max_density_jump = cmfgen_conv_cfg.get('max_density_jump', 0.10)
    max_ne_jump = cmfgen_conv_cfg.get('max_ne_jump', 0.10)
    max_logtau_jump = cmfgen_conv_cfg.get('max_logtau_jump', 0.10)
    partition     = 'interactive' if debug else stella_cfg.get('partition', heracles_cfg.get('partition', 'general'))
    memory        = stella_cfg.get('memory', heracles_cfg.get('memory', '16G'))
    # Enforce partition time limits: interactive = 1:59:00 max, general/small = 1-00:00:00 (1 day)
    if debug:
        time_limit = '01:59:00'
    else:
        time_limit = stella_cfg.get('time', '12:00:00')
        # Override to 1-day limit if partition is general or small and time is > 1 day
        if partition in ('general', 'small') and time_limit != '1-00:00:00':
            time_limit = '1-00:00:00'

    stella_root   = stella_cfg.get('stella_root', '/ptmp/sbrennan/stella/OpenStella')
    input_profile = stella_cfg.get('input_profile',
                        str(model_dir / 'stella' / 'input_profile'))

    stella_run_dir = model_dir / 'stella'
    stella_run_dir.mkdir(parents=True, exist_ok=True)

    # Generate .xni file if nickel56 is enabled
    nickel_cfg = stella_cfg.get('nickel56', {})
    if nickel_cfg.get('enabled', False):
        print(f"Generating .xni file from nickel56 configuration...")
        xni_script = scripts_dir / 'generate_xni_file.py'
        xni_output = stella_run_dir / 'eve' / f'{model_name}.xni'
        xni_output.parent.mkdir(parents=True, exist_ok=True)
        
        result = subprocess.run(
            [sys.executable, str(xni_script), str(config_file), str(xni_output)],
            capture_output=True,
            text=True
        )
        
        if result.returncode != 0:
            print(f"WARNING: .xni file generation failed")
            print(result.stdout)
            print(result.stderr)
        else:
            print(result.stdout)
    else:
        print("nickel56 not enabled, skipping .xni file generation")

    template_file = scripts_dir / 'stella_job_template.slurm'
    with open(template_file, 'r') as f:
        template = f.read()

    # Convert custom_output_times list to comma-separated string for template
    custom_output_times_str = ','.join(str(t) for t in custom_output_times) if custom_output_times else ''
    
    variables = {
        'MODEL_NAME':    model_name,
        'MODEL_DIR':     str(model_dir),
        'PROJECT_DIR':   str(project_dir),
        'SCRIPTS_DIR':   str(scripts_dir),
        'STELLA_ROOT':   stella_root,
        'STELLA_RUN_DIR': str(stella_run_dir),
        'INPUT_PROFILE': input_profile,
        'TARGET_DAYS':   target_days,
        'NZON':          nzon,
        'NZON_MAX':      nzon_max,
        'KAPPA':         kappa,
        'TAU_MAX':       tau_max,
        'TAU_MIN':       tau_min,
        'RHO_MIN_ARG':   f'--rho-min="{rho_min}"' if rho_min > 0 else '',
        'SNAPSHOT_INTERVAL': snapshot_interval,
        'R_MIN':         r_min,
        'EXPLOSION_ENERGY': explosion_energy,
        'EXPLOSION_MASS':   explosion_mass,
        'EXPLOSION_TIME':   explosion_time,
        'EXPLOSION_DELAY':  explosion_delay,
        'INITIAL_KINETIC_ENERGY': initial_kinetic_energy,
        'NSTMAX':           nstmax,
        'DISABLE_SCATTERING_ARG': '--disable-scattering' if disable_scattering else '',
        'EPS': eps,
        'HMAX': hmax,
        'ENABLE_CONV_ARG': '--enable-convection' if enable_convection else '',
        'FIT_TAU': fit_tau,
        'TAU_TOL': tau_tol,
        'ABS_FRACTION': abs_fraction,
        'GIVDTL_ARG': '--givdtl' if givdtl else '',
        'CMFGEN_SNAPSHOT_DAYS': cmfgen_snapshot_days,
        'ENABLE_GRID_THINNING_ARG': '--enable-grid-thinning' if enable_grid_thinning else '',
        'GRID_THINNING_THRESHOLD': grid_thinning_threshold,
        'ENABLE_ADAPTIVE_GRID_ARG': '--enable-adaptive-grid' if enable_adaptive_grid else '',
        'BASE_GRID_TARGET': base_grid_target,
        'MAX_RADIUS_JUMP': max_radius_jump,
        'MAX_TEMP_JUMP': max_temp_jump,
        'MAX_DENSITY_JUMP': max_density_jump,
        'MAX_NE_JUMP': max_ne_jump,
        'MAX_LOGTAU_JUMP': max_logtau_jump,
        'CUSTOM_OUTPUT_TIMES_ARG': f'--custom-output-times="{custom_output_times_str}"' if custom_output_times else '',
        'OUTER_FRACTION': outer_fraction,
        'TRANSITION_RADIUS': transition_radius,
        'PARTITION':     partition,
        'MEMORY':        memory,
        'TIME':          time_limit,
        'CONFIG_FILE':   str(config_file) if config_file else '',
    }

    job_script = substitute_template(template, variables)
    job_file = stella_run_dir / f'stella_{model_name}.slurm'
    with open(job_file, 'w') as f:
        f.write(job_script)

    result = subprocess.run(
        ['sbatch', str(job_file)],
        capture_output=True,
        text=True
    )
    if result.returncode != 0:
        print(f"ERROR: STELLA job submission failed")
        print(result.stdout)
        print(result.stderr)
        return None

    job_id = result.stdout.strip().split()[-1]
    print(f"STELLA job submitted: {job_id}")
    return job_id


def heracles_output_exists(config, project_dir):
    """Check if SN_HYDRO_DATA (conversion output) exists in the cmfgen directory."""
    model_name = config.get('model_name', 'unknown')
    output = config.get('output', {})
    model_dir = output.get('model_dir', '../models/${model_name}').replace('${model_name}', model_name)
    model_dir = (project_dir / model_dir).absolute()
    sn_hydro = model_dir / 'cmfgen' / 'SN_HYDRO_DATA'
    return sn_hydro.exists(), sn_hydro


def run_stella_conversion(config, scripts_dir, project_dir):
    """Run STELLA .swd -> SN_HYDRO_DATA conversion if output is missing."""
    model_name = config.get('model_name', 'unknown')
    output = config.get('output', {})
    model_dir = output.get('model_dir', '../models/${model_name}').replace('${model_name}', model_name)
    model_dir = (project_dir / model_dir).absolute()
    
    sn_hydro = model_dir / 'cmfgen' / 'SN_HYDRO_DATA'
    if sn_hydro.exists():
        print(f"  SN_HYDRO_DATA already exists: {sn_hydro}")
        return True
    
    # Find .swd file (check strad directory first, then main stella directory)
    swd_file = model_dir / 'stella' / 'strad' / f'{model_name}.swd'
    if not swd_file.exists():
        swd_file = model_dir / 'stella' / f'{model_name}.swd'
    if not swd_file.exists():
        print(f"ERROR: STELLA .swd not found for model {model_name}")
        print(f"  Checked: {model_dir}/stella/strad/{model_name}.swd")
        print(f"  Checked: {model_dir}/stella/{model_name}.swd")
        return False
    
    # Find .abn file
    abn_file = model_dir / 'stella' / 'modmake' / f'{model_name}.abn'
    if not abn_file.exists():
        abn_file = None
    
    # Get conversion parameters from config
    stella_cfg = config.get('stella', {})
    cmfgen_conv_cfg = config.get('cmfgen_conversion', {})
    target_days = cmfgen_conv_cfg.get('snapshot_days', stella_cfg.get('target_days', 40))
    kappa = stella_cfg.get('kappa', 0.34)
    tau_max = stella_cfg.get('tau_max', 30.0)
    tau_min = stella_cfg.get('tau_min', 0)
    enable_grid_thinning = cmfgen_conv_cfg.get('enable_grid_thinning', False)
    enable_adaptive_grid = cmfgen_conv_cfg.get('enable_adaptive_grid', True)
    base_grid_target = cmfgen_conv_cfg.get('base_grid_target', 100)
    max_radius_jump = cmfgen_conv_cfg.get('max_radius_jump', 0.05)
    max_temp_jump = cmfgen_conv_cfg.get('max_temp_jump', 0.10)
    max_density_jump = cmfgen_conv_cfg.get('max_density_jump', 0.10)
    max_ne_jump = cmfgen_conv_cfg.get('max_ne_jump', 0.10)
    max_logtau_jump = cmfgen_conv_cfg.get('max_logtau_jump', 0.10)
    max_points = cmfgen_conv_cfg.get('max_points', 500)
    
    print(f"\nConverting STELLA .swd to SN_HYDRO_DATA...")
    print(f"  .swd file: {swd_file}")
    print(f"  output: {sn_hydro}")
    
    cmd = [
        sys.executable, str(scripts_dir / 'convert_stella_swd_to_cmfgen.py'),
        str(swd_file), str(sn_hydro),
        '--target-days', str(target_days),
        '--kappa', str(kappa),
        '--tau-max', str(tau_max),
        '--tau-min', str(tau_min),
        '--model-name', model_name,
        '--base-grid-target', str(base_grid_target),
        '--max-radius-jump', str(max_radius_jump),
        '--max-temp-jump', str(max_temp_jump),
        '--max-density-jump', str(max_density_jump),
        '--max-ne-jump', str(max_ne_jump),
        '--max-logtau-jump', str(max_logtau_jump),
        '--max-points', str(max_points),
    ]
    if abn_file:
        cmd += ['--abn-file', str(abn_file)]
    if enable_grid_thinning:
        cmd.append('--enable-grid-thinning')
    if enable_adaptive_grid:
        cmd.append('--enable-adaptive-grid')
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"ERROR: STELLA conversion failed")
        print(result.stdout)
        print(result.stderr)
        return False
    
    print(f"  SN_HYDRO_DATA created: {sn_hydro}")
    return True


def stella_output_exists(config, project_dir):
    """Check if STELLA output folder exists (stella/ directory with .tt file)."""
    model_name = config.get('model_name', 'unknown')
    output = config.get('output', {})
    model_dir = output.get('model_dir', '../models/${model_name}').replace('${model_name}', model_name)
    model_dir = (project_dir / model_dir).absolute()
    stella_dir = model_dir / 'stella'
    tt_file = stella_dir / f'{model_name}.tt'
    return stella_dir.exists() and tt_file.exists(), stella_dir


def stella_is_complete(config, project_dir, stella_dir):
    """Check if STELLA run is complete by checking if simulation reached target time."""
    model_name = config.get('model_name', 'unknown')
    tt_file = stella_dir / f'{model_name}.tt'
    
    if not tt_file.exists():
        return False
    
    # Read .tt file to check final simulation time
    try:
        import numpy as np
        data = np.loadtxt(str(tt_file))
        if data.size == 0:
            return False
        # Last column is time in days
        final_time_days = data[-1, -1]
        target_days = config.get('stella', {}).get('target_days',
                         config.get('heracles', {}).get('target_days', 40))
        # Consider complete if within 1% of target or reached target
        return final_time_days >= target_days * 0.99
    except Exception:
        return False


def plot_stella_evolution_if_available(config, scripts_dir, project_dir):
    """Plot STELLA evolution from .swd file if available."""
    model_name = config.get('model_name', 'unknown')
    output = config.get('output', {})
    model_dir = output.get('model_dir', '../models/${model_name}').replace('${model_name}', model_name)
    model_dir = (project_dir / model_dir).absolute()
    
    # Find .swd file (check strad directory first, then main stella directory)
    swd_file = model_dir / 'stella' / 'strad' / f'{model_name}.swd'
    if not swd_file.exists():
        swd_file = model_dir / 'stella' / f'{model_name}.swd'
    if not swd_file.exists():
        return False
    
    plots_dir = model_dir / 'plots'
    plots_dir.mkdir(parents=True, exist_ok=True)
    output_png = plots_dir / 'stella_evolution.png'
    
    plot_script = scripts_dir / 'plot_stella_evolution.py'
    if not plot_script.exists():
        print(f"WARNING: plot_stella_evolution.py not found, skipping evolution plot")
        return False
    
    print(f"\n{'='*60}")
    print(f"Plotting STELLA evolution for {model_name}")
    print(f"{'='*60}")
    result = subprocess.run(
        [sys.executable, str(plot_script), str(swd_file), str(plots_dir)],
        capture_output=True, text=True
    )
    if result.returncode == 0:
        print(f"Evolution plot saved: {output_png}")
        return True
    else:
        print(f"WARNING: Evolution plotting failed")
        print(result.stderr)
        return False


def plot_stella_lightcurve_if_available(config, scripts_dir, project_dir):
    """Plot STELLA lightcurve if .tt file exists (e.g. from a previous run)."""
    model_name = config.get('model_name', 'unknown')
    output = config.get('output', {})
    model_dir = output.get('model_dir', '../models/${model_name}').replace('${model_name}', model_name)
    model_dir = (project_dir / model_dir).absolute()
    
    tt_file = model_dir / 'cmfgen' / f'{model_name}.tt'
    if not tt_file.exists():
        tt_file = model_dir / 'stella' / f'{model_name}.tt'
    if not tt_file.exists():
        return False
    
    plots_dir = model_dir / 'plots'
    plots_dir.mkdir(parents=True, exist_ok=True)
    output_png = plots_dir / 'stella_lightcurve.png'
    
    plot_script = scripts_dir / 'plot_stella_lightcurve.py'
    if not plot_script.exists():
        print(f"WARNING: plot_stella_lightcurve.py not found, skipping lightcurve plot")
        return False
    
    target_days = config.get('stella', {}).get('target_days',
                     config.get('heracles', {}).get('target_days', 40))
    
    print(f"\n{'='*60}")
    print(f"Plotting STELLA lightcurve for {model_name}")
    print(f"{'='*60}")
    result = subprocess.run(
        [sys.executable, str(plot_script), str(tt_file),
         '--output', str(output_png), '--days-max', str(target_days)],
        capture_output=True, text=True
    )
    if result.returncode == 0:
        print(f"Lightcurve plot saved: {output_png}")
        return True
    else:
        print(f"WARNING: Lightcurve plotting failed")
        print(result.stderr)
        return False


def process_config(config_file, scripts_dir, project_dir, debug=False, debug_days=1, redo_heracles=True, plot_lightcurve=True, skip_stella=False, run_cmfgen=None, run_sedona=None, run_sirocco=None):
    """Process a single configuration file.
    
    Args:
        run_cmfgen: If True, submit CMFGEN job to produce spectrum (requires SN_HYDRO_DATA). If None, read from config.
        run_sedona: If True, submit SEDONA job to produce spectrum (requires .mod file). If None, read from config.
        run_sirocco: If True, submit SIROCCO job to produce spectrum (requires SN_HYDRO_DATA). If None, read from config.
    """
    config = load_yaml_config(config_file)
    model_name = config.get('model_name', 'unknown')
    
    # Read spectrum code flags from config, allowing CLI override
    run_cmfgen = config.get('run_cmfgen', True) if run_cmfgen is None else run_cmfgen
    run_sedona = config.get('run_sedona', True) if run_sedona is None else run_sedona
    run_sirocco = config.get('run_sirocco', True) if run_sirocco is None else run_sirocco
    
    # Determine which codes are being used
    use_stella = config.get('hydro_backend', 'heracles').lower() == 'stella'
    hydro_code = 'STELLA' if use_stella else 'HERACLES'
    
    print(f"\n{'='*70}")
    print(f"Processing model: {model_name}")
    print(f"Config file: {config_file}")
    print(f"Hydro code: {hydro_code} (hydrodynamic radiation modeling)")
    active_codes = [c for c, flag in [('CMFGEN', run_cmfgen), ('SEDONA', run_sedona), ('SIROCCO', run_sirocco)] if flag]
    if active_codes:
        print(f"Spectrum codes: {' + '.join(active_codes)} (independent spectrum generation)")
    else:
        print(f"Spectrum code: None (all spectrum codes disabled)")
    if debug:
        print(f"DEBUG MODE ENABLED")
    print(f"{'='*70}")
    
    # Check for existing jobs
    existing_jobs = check_existing_jobs(model_name)
    if existing_jobs:
        print(f"\nWARNING: Found {len(existing_jobs)} existing job(s) for model {model_name}:")
        for job_id in existing_jobs:
            print(f"  - {job_id}")
        print(f"\nAutomatically cancelling existing jobs...")
        import subprocess
        try:
            subprocess.run(['scancel'] + existing_jobs, check=True, capture_output=True, text=True)
            print(f"Cancelled {len(existing_jobs)} job(s)")
        except subprocess.CalledProcessError as e:
            print(f"WARNING: Failed to cancel some jobs: {e.stderr}")
        # Wait a moment for cancellation to take effect
        import time
        time.sleep(2)
    
    # Extract paths
    output = config.get('output', {})
    
    # If debug mode, override simulation days and target days
    debug_config_file = None
    if debug:
        config = config.copy()
        config['model_name'] = f"{model_name}-debug"
        config['heracles'] = config.get('heracles', {}).copy()
        config['heracles']['simulation_days'] = debug_days
        config['heracles']['target_days'] = debug_days  # Target same as simulation for debug
        # Also override STELLA target_days
        config['stella'] = config.get('stella', {}).copy()
        config['stella']['target_days'] = debug_days
        # Update model_dir to include -debug suffix
        output['model_dir'] = output.get('model_dir', '../models/${model_name}').replace('${model_name}', f"{model_name}-debug")
        # Also update heracles_output path to include -debug suffix
        output['heracles_output'] = output.get('heracles_output', 'heracles').replace('${model_name}', f"{model_name}-debug")
        # Also update plots_dir path to include -debug suffix
        output['plots_dir'] = output.get('plots_dir', 'plots').replace('${model_name}', f"{model_name}-debug")
        # Write debug config file to model directory
        model_dir = (project_dir / output['model_dir']).absolute()
        model_dir.mkdir(parents=True, exist_ok=True)
        debug_config_path = model_dir / 'config_debug.yaml'
        with open(debug_config_path, 'w') as f:
            yaml.dump(config, f)
        debug_config_file = str(debug_config_path)
        config_file = debug_config_file
    
    # Step 1: Generate profile (for both HERACLES and STELLA)
    use_stella = config.get('hydro_backend', 'heracles').lower() == 'stella'
    hydro_backend = 'stella' if use_stella else 'heracles'
    if not generate_profile(config, scripts_dir, project_dir, config_file, hydro_backend):
        print(f"ERROR: Profile generation failed for {model_name}")
        return None
    
    # Step 1b: Read profile max radius and update Lbox_cm to match
    if not use_stella:
        output_cfg = config.get('output', {})
        heracles_out = output_cfg.get('heracles_output', 'heracles').replace('${model_name}', model_name)
        base_dir_cfg = output_cfg.get('base_dir', './models').replace('${model_name}', model_name)
        if not str(heracles_out).startswith('models') and not str(heracles_out).startswith('./models') and not str(heracles_out).startswith('../'):
            profile_path = project_dir / base_dir_cfg / model_name / heracles_out / 'input_profile'
        else:
            profile_path = (project_dir / heracles_out / 'input_profile').absolute()
        
        if profile_path.exists():
            import numpy as np
            data = np.loadtxt(str(profile_path), comments='#')
            r_max = data[-1, 0]  # max radius from profile (first column)
            old_lbox = float(config.get('heracles', {}).get('Lbox_cm', 6.0e16))
            config.setdefault('heracles', {})['Lbox_cm'] = float(r_max)
            print(f"\nUpdated Lbox_cm: {old_lbox:.2e} -> {r_max:.2e} (from profile R_max)")
            
            # Write updated config for namelist generator
            model_dir_resolved = output_cfg.get('model_dir', './models/${model_name}').replace('${model_name}', model_name)
            model_dir_resolved = (project_dir / model_dir_resolved).absolute()
            model_dir_resolved.mkdir(parents=True, exist_ok=True)
            config_file = str(model_dir_resolved / 'config_with_lbox.yaml')
            with open(config_file, 'w') as f:
                yaml.dump(config, f)
    
    # Step 2: Generate namelist (only for HERACLES)
    if not use_stella:
        if not generate_namelist(config, scripts_dir, project_dir, config_file, debug):
            print(f"ERROR: Namelist generation failed for {model_name}")
            return None
    
    # Step 3: Prepare CMFGEN directory (copy base files) - only if running CMFGEN
    if run_cmfgen:
        if not prepare_cmfgen_directory(config, scripts_dir, project_dir):
            print(f"ERROR: CMFGEN directory preparation failed for {model_name}")
            return None

    # Create output directory for model artifacts (lightcurve CSV, etc.)
    output_cfg = config.get('output', {})
    model_dir_resolved = output_cfg.get('model_dir', '../models/${model_name}').replace('${model_name}', model_name)
    model_dir_resolved = (project_dir / model_dir_resolved).absolute()
    output_dir = model_dir_resolved / 'output'
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"  Output directory: {output_dir}")

    # Step 4: Submit hydro job (STELLA or HERACLES) + conversion -> SN_HYDRO_DATA
    exists, sn_hydro_path = heracles_output_exists(config, project_dir)
    use_stella = config.get('hydro_backend', 'heracles').lower() == 'stella'
    

    if use_stella:
        # STELLA path: check for existing stella/ folder
        if skip_stella:
            print(f"\n--skip-stella set: skipping STELLA job submission")
            hydro_job_id = ""
        else:
            stella_exists, stella_path = stella_output_exists(config, project_dir)
            if stella_exists and not redo_heracles:
                # Check if STELLA run is complete
                is_complete = stella_is_complete(config, project_dir, stella_path)
                if is_complete:
                    print(f"\nFound complete STELLA output folder: {stella_path}")
                    print(f"  Simulation reached target time, skipping STELLA job submission")
                    print(f"  Proceeding directly to CMFGEN/SEDONA")
                    print(f"  (use --redo-heracles to force resubmission)")
                    hydro_job_id = ""
                else:
                    print(f"\nFound incomplete STELLA output folder: {stella_path}")
                    print(f"  (use --redo-heracles to force resubmission)")
                    response = input("Use existing incomplete stella/ folder? [Y/n]: ").strip().lower()
                    if response != 'n':
                        print(f"Using existing stella/ folder, skipping STELLA job submission")
                        hydro_job_id = ""
                    else:
                        print(f"Proceeding with new STELLA job submission...")
                        hydro_job_id = submit_stella_job(config, scripts_dir, project_dir, debug, config_file)
                        if hydro_job_id is None:
                            print(f"ERROR: STELLA job submission failed for {model_name}")
                            return None
            else:
                if stella_exists and redo_heracles:
                    print(f"\nFound existing STELLA output folder: {stella_path}")
                    print(f"  (--redo-heracles set, resubmitting anyway)")
                hydro_job_id = submit_stella_job(config, scripts_dir, project_dir, debug, config_file)
                if hydro_job_id is None:
                    print(f"ERROR: STELLA job submission failed for {model_name}")
                    return None
    else:
        # HERACLES path
        if not redo_heracles and exists:
            hydro_label = 'HERACLES'
            print(f"\nSkipping {hydro_label}: SN_HYDRO_DATA found at {sn_hydro_path}")
            print(f"  (use --redo-heracles to force resubmission)")
            hydro_job_id = ""
        else:
            if not redo_heracles and not exists:
                hydro_label = 'HERACLES'
                print(f"\nWARNING: --no-redo-heracles set but SN_HYDRO_DATA not found at {sn_hydro_path}")
                print(f"  Submitting {hydro_label} anyway...")
            hydro_job_id = submit_heracles_job(config, scripts_dir, project_dir, debug, config_file)
            if hydro_job_id is None:
                print(f"ERROR: HERACLES job submission failed for {model_name}")
                return None
    # Keep heracles_job_id alias for CMFGEN dependency
    heracles_job_id = hydro_job_id

    # Step 4b: Prepare SEDONA directory if SEDONA will be run (skip for STELLA - done in STELLA job)
    if run_sedona and not use_stella:
        if not prepare_sedona_directory(config, scripts_dir, project_dir):
            print(f"WARNING: SEDONA directory preparation failed for {model_name}")
            # Continue anyway to allow manual setup

    # Step 4c: Prepare SIROCCO directory if SIROCCO will be run
    if run_sirocco:
        if not prepare_sirocco_directory(config, scripts_dir, project_dir):
            print(f"WARNING: SIROCCO directory preparation failed for {model_name}")
            # Continue anyway to allow manual setup

    # Step 5: Submit spectrum generation jobs (CMFGEN, SEDONA and/or SIROCCO)
    cmfgen_job_ids = None
    sedona_job_id = None
    sirocco_job_id = None
    
    if run_cmfgen:
        print(f"\nSubmitting CMFGEN spectrum generation job...")
        cmfgen_job_ids = submit_cmfgen_jobs(config, scripts_dir, project_dir, heracles_job_id, debug, config_file)
        if cmfgen_job_ids is None:
            print(f"ERROR: CMFGEN jobs submission failed for {model_name}")
            # Continue anyway if SEDONA/SIROCCO also requested
    
    if run_sedona:
        print(f"\nSubmitting SEDONA spectrum generation job...")
        sedona_job_id = submit_sedona_job(config, scripts_dir, project_dir, heracles_job_id, debug)
        if sedona_job_id is None:
            print(f"ERROR: SEDONA job submission failed for {model_name}")
            # Continue anyway

    if run_sirocco:
        print(f"\nSubmitting SIROCCO spectrum generation job...")
        sirocco_job_id = submit_sirocco_job(config, scripts_dir, project_dir, heracles_job_id, debug)
        if sirocco_job_id is None:
            print(f"ERROR: SIROCCO job submission failed for {model_name}")
            # Continue anyway
    
    hydro_label = 'STELLA' if use_stella else 'HERACLES'
    print(f"\n{'='*70}")
    print(f"Model {model_name} job chain submitted successfully")
    print(f"  Hydro ({hydro_label}): {heracles_job_id if heracles_job_id else '(skipped - existing SN_HYDRO_DATA)'}")
    if cmfgen_job_ids:
        if cmfgen_job_ids['init'] is None:
            print(f"  CMFGEN:       {cmfgen_job_ids['final']} (single-final mode)")
        else:
            print(f"  CMFGEN INIT:  {cmfgen_job_ids['init']}")
            print(f"  CMFGEN FINAL: {cmfgen_job_ids['final']}")
        print(f"  CMFGEN OBS:   {cmfgen_job_ids['obs']}")
    if sedona_job_id:
        print(f"  SEDONA:       {sedona_job_id}")
    if sirocco_job_id:
        print(f"  SIROCCO:      {sirocco_job_id}")
    print(f"{'='*70}")

    return {
        'model_name': model_name,
        'heracles_job_id': heracles_job_id,
        'cmfgen_job_ids': cmfgen_job_ids,
        'sedona_job_id': sedona_job_id,
        'sirocco_job_id': sirocco_job_id,
    }


def main():
    parser = argparse.ArgumentParser(
        description='Master workflow script for automated supernova model grid'
    )
    parser.add_argument('configs', nargs='*', help='YAML config files (or use --all)')
    parser.add_argument('--all', action='store_true', help='Process all configs in configs/')
    parser.add_argument('--debug', action='store_true', help='Debug mode: use interactive partition, 1hr59m time limit')
    parser.add_argument('--debug-days', type=int, default=1, help='Simulation days in debug mode (default: 1)')
    parser.add_argument('--redo-heracles', dest='redo_heracles', action='store_true', default=True,
                        help='Resubmit HERACLES+conversion jobs (default: True)')
    parser.add_argument('--no-redo-heracles', dest='redo_heracles', action='store_false',
                        help='Skip HERACLES if SN_HYDRO_DATA already exists; resubmit only if missing')
    parser.add_argument('--plot', dest='plot_lightcurve', action='store_true', default=True,
                        help='Plot STELLA lightcurve if .tt file exists (default: True)')
    parser.add_argument('--no-plot', dest='plot_lightcurve', action='store_false',
                        help='Skip lightcurve plotting')
    parser.add_argument('--skip-stella', dest='skip_stella', action='store_true', default=False,
                        help='Skip STELLA run and use existing stella/ output (default: False)')
    parser.add_argument('--run-cmfgen', dest='run_cmfgen', action='store_true', default=None,
                        help='Submit CMFGEN job to produce spectrum (default: from config)')
    parser.add_argument('--no-cmfgen', dest='run_cmfgen', action='store_false',
                        help='Skip CMFGEN job')
    parser.add_argument('--run-sedona', dest='run_sedona', action='store_true', default=None,
                        help='Submit SEDONA job to produce spectrum (default: from config)')
    parser.add_argument('--no-sedona', dest='run_sedona', action='store_false',
                        help='Skip SEDONA job')
    parser.add_argument('--run-sirocco', dest='run_sirocco', action='store_true', default=None,
                        help='Submit SIROCCO job to produce spectrum (default: from config)')
    parser.add_argument('--no-sirocco', dest='run_sirocco', action='store_false',
                        help='Skip SIROCCO job')

    args = parser.parse_args()
    
    # Determine project directory
    scripts_dir = Path(__file__).parent
    project_dir = scripts_dir.parent
    
    # Get config files
    if args.all:
        configs_dir = project_dir / 'configs'
        config_files = list(configs_dir.glob('*.yaml'))
        if not config_files:
            print(f"ERROR: No YAML config files found in {configs_dir}")
            sys.exit(1)
    else:
        if not args.configs:
            print("ERROR: No config files specified. Use --all or specify config files.")
            sys.exit(1)
        config_files = [Path(f) for f in args.configs]
    
    print(f"{'='*70}")
    print(f"SN2024pba Model Grid Workflow")
    print(f"{'='*70}")
    print(f"Project directory: {project_dir}")
    print(f"Scripts directory: {scripts_dir}")
    print(f"Number of configs: {len(config_files)}")
    if args.debug:
        print(f"DEBUG MODE ENABLED")
    if not args.redo_heracles:
        print(f"HERACLES: skipped if SN_HYDRO_DATA exists (--no-redo-heracles)")
    if args.skip_stella:
        print(f"STELLA: skipped (--skip-stella)")
    if args.run_cmfgen is None:
        print(f"Spectrum generation: CMFGEN (from config)")
    elif args.run_cmfgen:
        print(f"Spectrum generation: CMFGEN enabled (--run-cmfgen)")
    else:
        print(f"Spectrum generation: CMFGEN disabled (--no-cmfgen)")
    if args.run_sedona is None:
        print(f"Spectrum generation: SEDONA (from config)")
    elif args.run_sedona:
        print(f"Spectrum generation: SEDONA enabled (--run-sedona)")
    else:
        print(f"Spectrum generation: SEDONA disabled (--no-sedona)")
    if args.run_sirocco is None:
        print(f"Spectrum generation: SIROCCO (from config)")
    elif args.run_sirocco:
        print(f"Spectrum generation: SIROCCO enabled (--run-sirocco)")
    else:
        print(f"Spectrum generation: SIROCCO disabled (--no-sirocco)")
    print(f"{'='*70}")
    
    # Process each config
    results = []
    for config_file in config_files:
        result = process_config(config_file, scripts_dir, project_dir, args.debug, args.debug_days, args.redo_heracles, args.plot_lightcurve, args.skip_stella, args.run_cmfgen, args.run_sedona, args.run_sirocco)
        if result:
            results.append(result)
    
    # Summary
    print(f"\n{'='*70}")
    print(f"Workflow Summary")
    print(f"{'='*70}")
    print(f"Total models: {len(results)}")
    print(f"\nJob IDs:")
    for result in results:
        print(f"  {result['model_name']}:")
        # Show hydro job info
        hydro_job = result['heracles_job_id']
        hydro_label = 'HERACLES'  # Default, could be enhanced to track actual hydro code used
        print(f"    Hydro:      {hydro_job if hydro_job else '(skipped - existing data)'}")
        # Show spectrum generation jobs
        if result.get('cmfgen_job_ids'):
            if result['cmfgen_job_ids']['init'] is None:
                print(f"    CMFGEN:     {result['cmfgen_job_ids']['final']} (single-final)")
            else:
                print(f"    CMFGEN:     INIT={result['cmfgen_job_ids']['init']}, FINAL={result['cmfgen_job_ids']['final']}")
            print(f"    CMFGEN OBS: {result['cmfgen_job_ids']['obs']}")
        if result.get('sedona_job_id'):
            print(f"    SEDONA:     {result['sedona_job_id']}")
        if result.get('sirocco_job_id'):
            print(f"    SIROCCO:    {result['sirocco_job_id']}")
    print(f"\nMonitor with: squeue -u $USER")
    print(f"{'='*70}")


if __name__ == '__main__':
    main()
