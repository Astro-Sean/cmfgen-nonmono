#!/usr/bin/env python3
"""
Master workflow script for automated supernova model grid.

This script:
1. Scans configs/ for YAML files
2. For each config:
   - Generates input profile (STELLA or [DEPRECATED] HERACLES)
   - Generates [DEPRECATED] HERACLES namelist (only if hydro_backend='heracles')
   - Submits STELLA job (named: stella_{model_name}) or [DEPRECATED] HERACLES job
   - Submits CMFGEN conversion job with dependency (named: convert_{model_name})
   - Submits CMFGEN run job with dependency (named: cmfgen_{model_name})
3. Returns job IDs for monitoring

NOTE: HERACLES is [DEPRECATED]. The project now uses the STELLA -> CMFGEN pipeline
exclusively for interacting Type IIn supernovae. HERACLES options still work at runtime
but are not recommended for new science runs.

Config options:
  - reuse_stella: If True, use existing STELLA output from current model directory instead of running STELLA
  - stella_source_model: Model name to copy STELLA output from (deprecated, no longer needed)
  - run_cmfgen: If True, submit CMFGEN job to produce spectrum
  - run_sedona: If True, submit SEDONA job to produce spectrum
  - run_sirocco: If True, submit SIROCCO job to produce spectrum
  - run_superlite: If True, submit SuperLite job to produce spectrum

Usage:
    python run_model_grid.py [config1.yaml config2.yaml ...]
    python run_model_grid.py --all  # Process all configs in configs/
"""

import sys
import yaml
import argparse
import subprocess
import shutil
from pathlib import Path
import os
import re
import numpy as np


DEBUG_TIME_LIMIT = '01:59:00'
DEFAULT_TIME_LIMIT = '24:00:00'

# ── Resource presets ────────────────────────────────────────────────
# Debug mode: interactive partition, 4 CPUs, 30G memory.
# Production mode: general partition, 8 CPUs, 32G memory.
# (The general partition may map to 'small' on some sites — this is fine.)
DEBUG_PARTITION = 'interactive'
DEBUG_CPUS = 4
DEBUG_MEM = '30G'
PROD_PARTITION = 'general'
PROD_CPUS = 8
PROD_MEM = '32G'


def get_job_resources(debug, config_partition=None, config_mem=None, config_cpus=None):
    """Return (partition, cpus, memory) for SLURM job submission.

    Debug mode always uses the interactive preset regardless of config.
    Production mode uses config values if provided, else the production preset.
    """
    if debug:
        return DEBUG_PARTITION, DEBUG_CPUS, DEBUG_MEM
    partition = config_partition if config_partition else PROD_PARTITION
    cpus = config_cpus if config_cpus else PROD_CPUS
    mem = config_mem if config_mem else PROD_MEM
    return partition, cpus, mem


def get_time_limit(debug, config_value=None):
    """Return the SLURM time limit.

    Single source of truth for all job types.
    - Debug mode: always DEBUG_TIME_LIMIT (01:30:00), config value ignored
      (interactive partition requires < 2h)
    - Non-debug: use config_value if provided, else DEFAULT_TIME_LIMIT (24:00:00)
    """
    if debug:
        return DEBUG_TIME_LIMIT
    return config_value if config_value else DEFAULT_TIME_LIMIT


def load_yaml_config(config_file):
    """Load YAML configuration file."""
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    return config


def parse_snapshot_days(config):
    """Return a sorted list of snapshot days from cmfgen_conversion.snapshot_days.

    Supports both a single value (legacy) and a list (time-series workflow).
    """
    cmfgen_conv = config.get('cmfgen_conversion', {})
    days = cmfgen_conv.get('snapshot_days', 40)
    if isinstance(days, (int, float, str)):
        days = [float(days) if isinstance(days, str) else days]
    days = sorted([float(d) for d in days])
    if not days:
        days = [40.0]
    return days


def get_stella_epochs(config, project_dir):
    """Read available epoch times (in days) from the STELLA .swd file.

    Returns a sorted list of epoch times, or None if the .swd file cannot be read.
    """
    model_name = config.get('model_name', 'unknown')
    output = config.get('output', {})
    model_dir = output.get('model_dir', '../models/${model_name}').replace('${model_name}', model_name)
    model_dir = (project_dir / model_dir).absolute()

    swd_source = config.get('cmfgen_conversion', {}).get('stella_source_swd')
    swd_file = Path(swd_source).expanduser() if swd_source else model_dir / 'stella' / 'strad' / f'{model_name}.swd'
    if not swd_file.exists() and not swd_source:
        swd_file = model_dir / 'stella' / f'{model_name}.swd'
    if not swd_file.exists():
        return None

    try:
        # Import the read_swd function from the conversion script
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "convert_stella_swd_to_cmfgen",
            str(project_dir / 'scripts' / 'convert_stella_swd_to_cmfgen.py'),
        )
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        swd_data = mod.read_swd(str(swd_file))
        epochs = sorted(swd_data.keys())
        return epochs
    except Exception as e:
        print(f"  WARNING: Could not read STELLA epochs from {swd_file}: {e}")
        return None


def filter_snapshot_days_by_stella_epochs(config, project_dir, snapshot_days):
    """Filter and remap snapshot_days to closest matching STELLA epochs within dt days.

    Reads the epoch_tolerance (dt) from cmfgen_conversion config.
    If the .swd file cannot be read or epoch_tolerance is not set, returns the
    original list unchanged.

    Accounts for the simulation start time: STELLA simulations often start at
    negative time (pre-SBO). The "days since explosion" is the time after t=0
    (SBO). For positive target days, only positive (post-SBO) epochs are
    considered as candidates.

    If a requested day has no matching epoch within dt tolerance, it is skipped.
    If a requested day has a matching epoch, the day is remapped to the closest
    epoch value (so the CMFGEN directory name reflects the actual STELLA epoch).
    """
    cmfgen_conv = config.get('cmfgen_conversion', {})
    dt = cmfgen_conv.get('epoch_tolerance', None)

    # If no tolerance is configured, don't filter
    if dt is None:
        return snapshot_days

    dt = float(dt)
    epochs = get_stella_epochs(config, project_dir)
    if epochs is None:
        print(f"  WARNING: Cannot read STELLA epochs — skipping epoch tolerance filter")
        return snapshot_days

    # Report simulation start time for context
    sim_start = epochs[0]
    print(f"\n  STELLA .swd has {len(epochs)} epochs: "
          f"{epochs[0]:.2f} ... {epochs[-1]:.2f} days")
    print(f"  Simulation starts at t={sim_start:.2f} d (pre-SBO if negative)")
    print(f"  Epoch tolerance (dt): {dt} days")

    filtered = []
    skipped = []
    for day in snapshot_days:
        # Find closest epoch (prefer positive epochs for positive target days)
        # This accounts for the simulation start time: post-SBO epochs (t>=0)
        # are the relevant ones for "days since explosion" matching.
        if day >= 0:
            positive = [t for t in epochs if t >= 0]
            candidates = positive if positive else epochs
        else:
            negative = [t for t in epochs if t < 0]
            candidates = negative if negative else epochs
        if not candidates:
            skipped.append(day)
            continue
        closest = min(candidates, key=lambda t: abs(t - day))
        delta = abs(closest - day)
        if delta <= dt:
            # Don't remap — the conversion script will select the closest epoch.
            # The directory name uses the requested day (e.g., "20.0days").
            filtered.append(day)
        else:
            skipped.append(day)
            print(f"  SKIP day {day}: closest STELLA epoch is {closest:.2f} d "
                  f"(delta={delta:.2f} > dt={dt})")

    if skipped:
        print(f"  Skipping {len(skipped)} snapshot day(s) with no matching STELLA epoch: {skipped}")
    if filtered:
        print(f"  Submitting jobs for {len(filtered)} snapshot day(s): {filtered}")

    return filtered


def get_stella_target_days(config):
    """Compute STELLA target_days from snapshot days + extra margin.

    STELLA must model the full evolution for all snapshot days.
    """
    snapshot_days = parse_snapshot_days(config)
    cmfgen_conv = config.get('cmfgen_conversion', {})
    extra = cmfgen_conv.get('stella_extra_days', 10)
    if snapshot_days:
        return max(snapshot_days) + float(extra)
    else:
        return float(extra)


def calculate_kinetic_energy_from_profile(profile_path, r_min=0.0):
    """Compute total kinetic energy (in 1e50 erg) from an input_profile.

    The profile is expected to have columns:
        r_cm  rho_cgs  p_cgs  u_cgs  fx_H  fx_He  ...

    Shell masses are computed using the same spherical-shell convention as
    convert_profile_to_stella.py (inner boundary = previous radius or centre).

    If r_min > 0, zones with r < r_min are excluded (matching STELLA's
    --r-min trimming) so the computed KE matches what STELLA actually sees.
    """
    if not Path(profile_path).exists():
        print(f"  WARNING: Profile {profile_path} not found, cannot compute KE")
        return None

    try:
        data = np.loadtxt(str(profile_path), comments='#')
        if data.ndim == 1:
            data = data.reshape(1, -1)
        r = data[:, 0]          # cm
        rho = data[:, 1]        # g/cm^3
        u = data[:, 3]          # cm/s  (column index 3 = u_cgs)

        # Apply r_min trimming to match STELLA's --r-min filtering
        r_min = float(r_min)
        if r_min > 0:
            mask = r >= r_min
            n_excluded = len(r) - np.sum(mask)
            if n_excluded > 0:
                r = r[mask]
                rho = rho[mask]
                u = u[mask]
                print(f"  (Excluded {n_excluded} zones with r < {r_min:.2e} cm)")

        # Match shell-mass prescription in convert_profile_to_stella.py:
        # inner boundary = rcen (defaults to r[0]*0.95 in converter)
        rcen = r[0] * 0.95
        r_lo = np.concatenate([[rcen], r[:-1]])
        r_hi = r
        dm = (4.0 / 3.0) * np.pi * (r_hi**3 - r_lo**3) * rho
        dm = np.maximum(dm, 0.0)

        ekin_cgs = 0.5 * np.sum(dm * u**2)
        ekin_1e50 = ekin_cgs / 1e50
        m_total_msun = np.sum(dm) / 1.989e33

        print(f"  Auto-computed KE from profile:")
        print(f"    M_total = {m_total_msun:.3f} Msun")
        print(f"    E_kin   = {ekin_1e50:.3f} x 1e50 erg")
        print(f"    v_rms   = {np.sqrt(np.average(u**2, weights=dm))/1e5:.1f} km/s")
        return float(ekin_1e50)
    except Exception as exc:
        print(f"  WARNING: Failed to compute KE from profile: {exc}")
        return None


def substitute_template(template, variables):
    """Substitute variables in template string."""
    result = template
    for key, value in variables.items():
        placeholder = "{" + key + "}"
        result = result.replace(placeholder, str(value))
    return result


def generate_profile(config, scripts_dir, project_dir, config_path, hydro_backend='heracles'):
    """Generate input profile from YAML config.

    [DEPRECATED] The default hydro_backend='heracles' is deprecated.
    Use 'stella' as the hydro backend instead (STELLA -> CMFGEN pipeline).

    If a 'kepler_presn' section is present in the config, uses
    convert_kepler_presn_to_stella.py to convert a KEPLER pre-SN model
    instead of generate_profile_from_yaml.py.
    """
    model_name = config.get('model_name', 'unknown')
    print(f"\n{'='*60}")
    print(f"Generating profile for model: {model_name}")
    print(f"{'='*60}")

    # Check for KEPLER pre-SN mode
    kepler_cfg = config.get('kepler_presn', None)
    if kepler_cfg:
        return generate_kepler_presn_profile(config, scripts_dir, project_dir, config_path)

    # Check for SWD snapshot restart mode (extract .hyd/.abn from existing .swd)
    swd_snap_cfg = config.get('stella_from_swd_snapshot', None)
    if swd_snap_cfg:
        return generate_swd_snapshot_profile(config, scripts_dir, project_dir, config_path)

    # Check for SN_B mode (H/He shell inside SN_A ejecta)
    sn_b_cfg = config.get('sn_b_from_sn_a', None)
    if sn_b_cfg:
        return generate_sn_b_profile(config, scripts_dir, project_dir, config_path)

    # Check for X model mode (H/He shell + wind + old SN ejecta)
    x_model_cfg = config.get('x_model_from_sn', None)
    if x_model_cfg:
        return generate_x_model_profile(config, scripts_dir, project_dir, config_path)

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


def generate_kepler_presn_profile(config, scripts_dir, project_dir, config_path):
    """Generate input_profile from a KEPLER pre-SN model.

    Calls convert_kepler_presn_to_stella.py to:
    1. Read the KEPLER pre-SN file
    2. Excise the inner compact remnant (mass_cut)
    3. Remap onto nzon STELLA zones
    4. Add H surface layer if configured
    5. Write .hyd, .abn, .xni, .dat, and input_profile files

    The input_profile is written to {model_dir}/stella/input_profile
    so the rest of the pipeline (stella_job_template.slurm) can find it.
    """
    model_name = config.get('model_name', 'unknown')
    kepler_cfg = config.get('kepler_presn', {})
    stella_cfg = config.get('stella', {})
    output_cfg = config.get('output', {})

    presn_file = kepler_cfg.get('presn_file', '')
    mass_cut = kepler_cfg.get('mass_cut_msun', 1.4)
    h_surface_fraction = kepler_cfg.get('h_surface_fraction', 0.0)
    h_surface_mass_msun = kepler_cfg.get('h_surface_mass_msun', 0.0)
    t_floor = kepler_cfg.get('t_floor_K', 3000.0)

    # Extended H envelope parameters (Bersten 2012 approach for SN 2011dh:
    # compact He core + extended ~0.1 Msun H envelope at R~200 Rsun)
    h_ext_fraction = kepler_cfg.get('h_ext_fraction', 0.0)
    h_ext_mass_msun = kepler_cfg.get('h_ext_mass_msun', 0.0)
    h_ext_r_outer_rsun = kepler_cfg.get('h_ext_r_outer_rsun', 0.0)
    h_ext_n_zones = kepler_cfg.get('h_ext_n_zones', 50)
    h_ext_T_surface_K = kepler_cfg.get('h_ext_T_surface_K', 6000.0)
    h_ext_T_base_K = kepler_cfg.get('h_ext_T_base_K', None)
    h_ext_n_rho = kepler_cfg.get('h_ext_n_rho', 3.5)
    h_ext_n_T = kepler_cfg.get('h_ext_n_T', 2.0)

    # CSM wind parameters (prevents radiation solver stall when photosphere
    # exits the ejecta at ~10 days)
    wind_cfg = kepler_cfg.get('csm_wind', {})
    wind_zones = wind_cfg.get('n_zones', 0)
    wind_mdot = wind_cfg.get('mdot_msun_yr', 0.0)
    wind_velocity = wind_cfg.get('velocity_kms', 100.0)
    wind_duration = wind_cfg.get('duration_yr', 1000.0)
    wind_T = wind_cfg.get('T_K', 5000.0)
    wind_density_profile = wind_cfg.get('density_profile', 'uniform')
    wind_power_law_index = wind_cfg.get('power_law_index', None)

    # Surface density reduction (smoother outer profile)
    reduce_cfg = kepler_cfg.get('reduce_surface_density', {})
    reduce_surf_zones = reduce_cfg.get('n_zones', 0)
    reduce_surf_rho = reduce_cfg.get('rho_target', 0.0)

    nzon = stella_cfg.get('nzon', 500)
    target_days = stella_cfg.get('target_days', 60.0)
    explosion_energy = stella_cfg.get('explosion_energy', 10.0)
    explosion_mass = stella_cfg.get('explosion_mass', 0.5)
    explosion_time = stella_cfg.get('explosion_time', 1.0)
    initial_kinetic_energy = stella_cfg.get('initial_kinetic_energy', 0.0)

    # rho_min: trim outer KEPLER zones with density below this threshold before
    # remapping.  Removes the extended low-density atmosphere that causes STELLA
    # radiation-transport instabilities (BADSTE / KFLAG=-2) when the photosphere
    # reaches the very low density outer region.
    rho_min = stella_cfg.get('rho_min') or 0.0
    if isinstance(rho_min, str):
        rho_min = float(rho_min)

    # Ni56 configuration
    nickel_cfg = stella_cfg.get('nickel56', {})
    ni56_mass = nickel_cfg.get('total_mass_msun', 0.06)
    ni56_dist = nickel_cfg.get('distribution_mass_msun', 0.5)

    # Resolve output directory (same logic as generate_profile_from_yaml.py)
    model_dir = output_cfg.get('model_dir', './models/${model_name}').replace('${model_name}', model_name)
    stella_output = output_cfg.get('stella_output', 'stella').replace('${model_name}', model_name)

    # Resolve to absolute path
    if str(model_dir).startswith('/'):
        model_dir_abs = Path(model_dir)
    else:
        model_dir_abs = (project_dir / model_dir).absolute()

    if str(stella_output).startswith('/'):
        stella_out_abs = Path(stella_output)
    elif not str(stella_output).startswith('models') and not str(stella_output).startswith('./models'):
        stella_out_abs = model_dir_abs / stella_output
    else:
        stella_out_abs = (project_dir / stella_output).absolute()

    stella_out_abs.mkdir(parents=True, exist_ok=True)

    # Also create modmake subdirectory for .hyd/.abn/.dat files
    modmake_dir = stella_out_abs / 'modmake'
    modmake_dir.mkdir(parents=True, exist_ok=True)

    print(f"  KEPLER pre-SN mode:")
    print(f"    presn_file: {presn_file}")
    print(f"    mass_cut: {mass_cut} Msun")
    print(f"    h_surface_fraction: {h_surface_fraction}")
    print(f"    h_surface_mass_msun: {h_surface_mass_msun}")
    if h_ext_mass_msun > 0:
        print(f"    h_ext_fraction: {h_ext_fraction}")
        print(f"    h_ext_mass_msun: {h_ext_mass_msun}")
        print(f"    h_ext_r_outer_rsun: {h_ext_r_outer_rsun}")
        print(f"    h_ext_n_zones: {h_ext_n_zones}")
    print(f"    output_dir: {stella_out_abs}")

    script = scripts_dir / 'convert_kepler_presn_to_stella.py'
    cmd = [
        sys.executable, str(script),
        str(presn_file),
        str(model_name),
        '--mass-cut', str(mass_cut),
        '--nzon', str(nzon),
        '--target-days', str(target_days),
        '--explosion-energy', str(explosion_energy),
        '--explosion-mass', str(explosion_mass),
        '--explosion-time', str(explosion_time),
        '--initial-kinetic-energy', str(initial_kinetic_energy),
        '--ni56-mass', str(ni56_mass),
        '--ni56-dist', str(ni56_dist),
        '--t-floor', str(t_floor),
        '--outdir', str(modmake_dir),
    ]
    if h_surface_fraction > 0 and h_surface_mass_msun > 0:
        cmd.extend(['--h-surface-fraction', str(h_surface_fraction),
                    '--h-surface-mass-msun', str(h_surface_mass_msun)])
    if h_ext_fraction > 0 and h_ext_mass_msun > 0 and h_ext_r_outer_rsun > 0:
        cmd.extend(['--h-ext-fraction', str(h_ext_fraction),
                    '--h-ext-mass-msun', str(h_ext_mass_msun),
                    '--h-ext-r-outer-rsun', str(h_ext_r_outer_rsun),
                    '--h-ext-n-zones', str(h_ext_n_zones),
                    '--h-ext-T-surface-K', str(h_ext_T_surface_K),
                    '--h-ext-T-base-K', str(h_ext_T_base_K) if h_ext_T_base_K is not None else 'None',
                    '--h-ext-n-rho', str(h_ext_n_rho),
                    '--h-ext-n-T', str(h_ext_n_T)])
    if rho_min > 0:
        cmd.extend(['--rho-min', str(rho_min)])

    # CSM wind: append low-density wind/extension beyond ejecta to prevent
    # radiation solver stall when photosphere exits ejecta (~10d for Type IIb).
    # 'extrapolate' mode needs no Mdot -- it derives density purely from the
    # ejecta's own surface value/slope (negligible-mass, gradient-matched).
    if wind_zones > 0 and (wind_mdot > 0 or wind_density_profile in ('extrapolate', 'sigmoid')):
        cmd.extend(['--wind-zones', str(wind_zones),
                    '--wind-mdot', str(wind_mdot),
                    '--wind-velocity', str(wind_velocity),
                    '--wind-duration', str(wind_duration),
                    '--wind-T', str(wind_T),
                    '--wind-density-profile', str(wind_density_profile)])
        if wind_power_law_index is not None:
            cmd.extend(['--wind-power-law-index', str(wind_power_law_index)])
        print(f"    csm_wind: {wind_zones} zones, Mdot={wind_mdot:.1e} Msun/yr, "
              f"v={wind_velocity:.0f} km/s, t={wind_duration:.0f} yr, "
              f"profile={wind_density_profile}", end='')
        if wind_power_law_index is not None:
            print(f", n={wind_power_law_index}", end='')
        print()

    # Surface density reduction
    if reduce_surf_zones > 0 and reduce_surf_rho > 0:
        cmd.extend(['--reduce-surf-zones', str(reduce_surf_zones),
                    '--reduce-surf-rho', str(reduce_surf_rho)])
        print(f"    reduce_surface_density: {reduce_surf_zones} zones -> {reduce_surf_rho:.1e} g/cc")

    # Smooth steep density gradients if configured
    # KEPLER pre-SN models can have |dlog rho/dlog r| > 8 at the He core edge,
    # which causes STELLA GEAR+BGH solver to stall (BADSTE / KFLAG=-2) during
    # the thermal bomb phase.  Mild Gaussian smoothing reduces the steepest
    # gradients without significantly altering the profile.
    smooth_profile = kepler_cfg.get('smooth_profile', False)
    max_log_gradient = kepler_cfg.get('smooth_max_gradient', 5.0)
    smooth_sigma = kepler_cfg.get('smooth_sigma', 3.0)
    if smooth_profile:
        cmd.extend(['--smooth-profile',
                    '--max-log-gradient', str(max_log_gradient),
                    '--smooth-sigma', str(smooth_sigma)])

    # Composition-interface smoothing (runs unconditionally, independent of
    # smooth_profile, inside convert_kepler_presn_to_stella.py). Relax these
    # from the tight built-in defaults to preserve sharp shell boundaries
    # (e.g. the O-shell) when requested via config.
    comp_smooth_sigma = kepler_cfg.get('comp_smooth_sigma', 10.0)
    comp_smooth_max_jump = kepler_cfg.get('comp_smooth_max_jump', 0.02)
    comp_smooth_max_dx = kepler_cfg.get('comp_smooth_max_dx', 0.08)
    cmd.extend(['--comp-smooth-sigma', str(comp_smooth_sigma),
                '--comp-smooth-max-jump', str(comp_smooth_max_jump),
                '--comp-smooth-max-dx', str(comp_smooth_max_dx)])
    print(f"    comp_smoothing: sigma={comp_smooth_sigma}, max_jump={comp_smooth_max_jump}, max_dx={comp_smooth_max_dx}")

    # Pass plot directory for composition smoothing before/after plot
    plots_dir = model_dir_abs / 'plots'
    plots_dir.mkdir(parents=True, exist_ok=True)
    cmd.extend(['--plot-dir', str(plots_dir)])

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"ERROR: KEPLER pre-SN profile generation failed")
        print(result.stdout)
        print(result.stderr)
        return False

    print(result.stdout)

    # Copy input_profile to the stella output directory (where the slurm template expects it)
    import shutil
    src_profile = modmake_dir / 'input_profile'
    dst_profile = stella_out_abs / 'input_profile'
    if src_profile.exists():
        shutil.copy2(str(src_profile), str(dst_profile))
        print(f"  Copied input_profile: {src_profile} -> {dst_profile}")

    return True


def generate_swd_snapshot_profile(config, scripts_dir, project_dir, config_path):
    """Extract a STELLA snapshot from an existing .swd file for restart.

    Calls extract_swd_snapshot_to_stella.py to:
    1. Read the specified .swd file
    2. Find the epoch closest to target_days
    3. Write .hyd and .abn files to the modmake directory
    4. Create a .dat file with EKO=0 (preserve velocities) and TcurA/TcurB
       set for the restart time range
    """
    model_name = config.get('model_name', 'unknown')
    swd_cfg = config.get('stella_from_swd_snapshot', {})
    swd_file = swd_cfg.get('swd_file', '')
    target_days = swd_cfg.get('target_days', 35.0)

    if not swd_file or not Path(swd_file).exists():
        print(f"ERROR: SWD file not found: {swd_file}", file=sys.stderr)
        return False

    output_cfg = config.get('output', {})
    model_dir_str = output_cfg.get('model_dir', '.')
    if str(model_dir_str).startswith('/'):
        model_dir = Path(model_dir_str)
    else:
        model_dir = project_dir / model_dir_str
    stella_output = output_cfg.get('stella_output', 'stella')
    stella_out_abs = model_dir / stella_output
    modmake_dir = stella_out_abs / 'modmake'
    modmake_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n  SWD snapshot restart mode")
    print(f"  Source SWD: {swd_file}")
    print(f"  Target days: {target_days}")

    script = scripts_dir / 'extract_swd_snapshot_to_stella.py'
    result = subprocess.run(
        [sys.executable, str(script), str(swd_file), str(modmake_dir),
         f'--target-days={target_days}', f'--model-name={model_name}'],
        capture_output=True, text=True
    )
    print(result.stdout)
    if result.returncode != 0:
        print(f"ERROR: SWD snapshot extraction failed", file=sys.stderr)
        print(result.stderr, file=sys.stderr)
        return False

    # Create .dat file with EKO=0 (preserve velocities from .hyd)
    # and TcurA=target_days, TcurB=stella.target_days
    stella_cfg = config.get('stella', {})
    tcurb = stella_cfg.get('target_days', 42.0)

    # Read the original .dat file to use as template
    orig_dat = Path(swd_file).parent / f'{model_name}.dat'
    if orig_dat.exists():
        with open(orig_dat, 'r') as f:
            dat_lines = f.readlines()
        # Modify NSTA/NSTB line (line 7, 0-indexed=6) and TcurA/TcurB
        # Format: abs(NSTA)>abs(NSTB)=>reads last   TcurA     TcurB
        #         -9000     -1                      5.0000   84.8589
        # For fresh start from .hyd: NSTA=1, NSTB=-1 (no restart from .prf)
        # TcurA = target_days (start from snapshot time)
        # TcurB = stella.target_days (end time)
        for i, line in enumerate(dat_lines):
            if 'TcurA' in line and 'TcurB' in line:
                # This is the header line, modify the next line
                if i + 1 < len(dat_lines):
                    parts = dat_lines[i + 1].split()
                    if len(parts) >= 4:
                        # NSTA=1 (fresh start, read .hyd), NSTB=-1
                        parts[0] = '1'
                        parts[1] = '-1'
                        parts[2] = f'{target_days:.4f}'
                        parts[3] = f'{tcurb:.4f}'
                        dat_lines[i + 1] = '   ' + '     '.join(parts) + '\n'
                break

        # Also set EKO=0 to preserve velocities
        for i, line in enumerate(dat_lines):
            if 'EKO' in line and 'kin.energy' in line:
                if i + 1 < len(dat_lines):
                    parts = dat_lines[i + 1].split()
                    if len(parts) >= 1:
                        parts[0] = '0.0d0'
                        dat_lines[i + 1] = ' ' + '   '.join(parts) + '\n'
                break

        dat_path = modmake_dir / f'{model_name}.dat'
        with open(dat_path, 'w') as f:
            f.writelines(dat_lines)
        print(f"  Written .dat: {dat_path} (TcurA={target_days}, TcurB={tcurb}, EKO=0)")
    else:
        print(f"  WARNING: Original .dat not found at {orig_dat}, using template default")

    # Copy .hyd and .abn to stella run directory
    stella_run_dir = stella_out_abs
    for ext in ['.hyd', '.abn']:
        src = modmake_dir / f'{model_name}{ext}'
        dst = stella_run_dir / f'{model_name}{ext}'
        if src.exists():
            shutil.copy2(str(src), str(dst))
            print(f"  Copied {ext}: {src} -> {dst}")

    # Copy .dat to stella run directory
    dat_src = modmake_dir / f'{model_name}.dat'
    dat_dst = stella_run_dir / f'{model_name}.dat'
    if dat_src.exists():
        shutil.copy2(str(dat_src), str(dat_dst))
        print(f"  Copied .dat: {dat_src} -> {dat_dst}")

    return True


def generate_sn_b_profile(config, scripts_dir, project_dir, config_path):
    """Generate SN_B initial model from SN_A ejecta + H/He shell.

    Calls prepare_sn_b.py to:
    1. Read the SN_A STELLA .swd output at 349 days
    2. Homologously evolve to the requested age if needed
    3. Create a H/He shell (given mass, velocity, composition)
    4. Combine shell (inner) + SN_A ejecta (outer) into STELLA .hyd format
    5. Write .hyd, .abn, .xni, .dat, and input_profile files
    """
    model_name = config.get('model_name', 'unknown')
    sn_b_cfg = config.get('sn_b_from_sn_a', {})
    stella_cfg = config.get('stella', {})
    output_cfg = config.get('output', {})

    sn_a_dir = sn_b_cfg.get('sn_a_dir', '')
    sn_a_age = sn_b_cfg.get('sn_a_age_days', 349.0)
    shell_mass = sn_b_cfg.get('shell_mass_msun', 1.0)
    shell_velocity = sn_b_cfg.get('shell_velocity_kms', 5000)
    shell_xh = sn_b_cfg.get('shell_xh', 0.7)
    shell_yhe = sn_b_cfg.get('shell_yhe', 0.3)
    explosion_energy = sn_b_cfg.get('explosion_energy', 0.0)
    aq = sn_b_cfg.get('aq', 4.0)

    nzon = stella_cfg.get('nzon', 500)
    target_days = stella_cfg.get('target_days', 50.0)

    # Get snapshot day for .dat generation
    cmfgen_conv = config.get('cmfgen_conversion', {})
    snapshot_days = cmfgen_conv.get('snapshot_days', [40])
    snapshot_day = snapshot_days[0] if isinstance(snapshot_days, list) else snapshot_days

    # Resolve output directory
    model_dir = output_cfg.get('model_dir', './models/${model_name}').replace('${model_name}', model_name)
    if str(model_dir).startswith('/'):
        model_dir_abs = Path(model_dir)
    else:
        model_dir_abs = (project_dir / model_dir).absolute()

    model_dir_abs.mkdir(parents=True, exist_ok=True)

    print(f"  SN_B mode (H/He shell inside SN_A ejecta):")
    print(f"    SN_A directory: {sn_a_dir}")
    print(f"    SN_A age: {sn_a_age} days")
    print(f"    Shell: M={shell_mass} Msun, v={shell_velocity} km/s, X={shell_xh}, Y={shell_yhe}")
    print(f"    Target days: {target_days}, Snapshot day: {snapshot_day}")
    print(f"    Output: {model_dir_abs}")

    script = scripts_dir / 'prepare_sn_b.py'
    cmd = [
        sys.executable, str(script),
        str(sn_a_dir),
        str(model_dir_abs),
        '--shell-mass', str(shell_mass),
        '--shell-velocity', str(shell_velocity),
        '--shell-xh', str(shell_xh),
        '--shell-yhe', str(shell_yhe),
        '--target-days', str(target_days),
        '--snapshot-day', str(snapshot_day),
        '--nzon', str(nzon),
        '--sn-a-age', str(sn_a_age),
        '--explosion-energy', str(explosion_energy),
        '--aq', str(aq),
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"ERROR: SN_B profile generation failed")
        print(result.stdout)
        print(result.stderr)
        return False

    print(result.stdout)
    return True


def generate_x_model_profile(config, scripts_dir, project_dir, config_path):
    """Generate X initial model from H/He shell + wind + old SN ejecta.

    Calls prepare_x_model.py to:
    1. Read the source STELLA .swd output (latest epoch)
    2. Homologously evolve to the requested age (e.g. 400 days)
    3. Create a H/He shell (given mass, velocity, composition)
    4. Create a connecting wind (Mdot, v_wind)
    5. Combine shell (inner) + wind (middle) + ejecta (outer) into STELLA format
    6. Write .hyd, .abn, .xni, .dat, and input_profile files
    """
    model_name = config.get('model_name', 'unknown')
    x_cfg = config.get('x_model_from_sn', {})
    stella_cfg = config.get('stella', {})
    output_cfg = config.get('output', {})

    source_dir = x_cfg.get('source_dir', '')
    source_age = x_cfg.get('source_age_days', 400.0)
    shell_mass = x_cfg.get('shell_mass_msun', 1.0)
    shell_velocity = x_cfg.get('shell_velocity_kms', 5000)
    shell_xh = x_cfg.get('shell_xh', 0.7)
    shell_yhe = x_cfg.get('shell_yhe', 0.3)
    shell_composition = x_cfg.get('shell_composition', 'hhe')
    n_rho = x_cfg.get('n_rho', 8.0)
    v0_kms = x_cfg.get('v0_kms', 3000.0)
    n_T = x_cfg.get('n_T', 0.4)
    T_interface = x_cfg.get('T_interface_K', 5000.0)
    T_max = x_cfg.get('T_max_K', 65000.0)
    t_start = x_cfg.get('t_start_days', 10.0)
    ejecta_T_floor = x_cfg.get('ejecta_T_floor_K', 1000.0)
    ejecta_rho_min = x_cfg.get('ejecta_rho_min', 0.0)
    ejecta_r_max_cm = x_cfg.get('ejecta_r_max_cm', 0.0)
    wind_mdot = x_cfg.get('wind_mdot_msun_yr', 1e-3)
    wind_velocity = x_cfg.get('wind_velocity_kms', 600)
    wind_T = x_cfg.get('wind_T_K', 5000)
    no_wind = x_cfg.get('no_wind', False)
    se_transition_dex = x_cfg.get('se_transition_dex', 0.5)
    density_matched_wind = x_cfg.get('density_matched_wind', False)
    cd_transition_zones = x_cfg.get('cd_transition_zones', 5)
    cd_sigmoid_k = x_cfg.get('cd_sigmoid_k', 5.0)
    v_blend_zones = x_cfg.get('v_blend_zones', 30)
    force_homologous_velocity = x_cfg.get('force_homologous_velocity', False)
    explosion_energy = x_cfg.get('explosion_energy', 0.0)
    aq = x_cfg.get('aq', 10.0)
    fit_tau = x_cfg.get('fit_tau', 5.0)
    meth = x_cfg.get('meth', 2)
    temp_floor = x_cfg.get('temp_floor', 0.03)
    floor_rad = x_cfg.get('floor_rad', 0.01)
    hmin = x_cfg.get('hmin', 1e-14)
    hmax = x_cfg.get('hmax', 3600.0)

    nzon = stella_cfg.get('nzon', 500)
    target_days = stella_cfg.get('target_days', 45.0)

    # Get snapshot day for .dat generation
    cmfgen_conv = config.get('cmfgen_conversion', {})
    snapshot_days = cmfgen_conv.get('snapshot_days', [30])
    snapshot_day = snapshot_days[0] if isinstance(snapshot_days, list) else snapshot_days

    # Resolve output directory
    model_dir = output_cfg.get('model_dir', './${model_name}').replace('${model_name}', model_name)
    if str(model_dir).startswith('/'):
        model_dir_abs = Path(model_dir)
    else:
        model_dir_abs = (project_dir / model_dir).absolute()

    model_dir_abs.mkdir(parents=True, exist_ok=True)

    print(f"  X model mode (H/He shell + wind + old SN ejecta):")
    print(f"    Source directory: {source_dir}")
    print(f"    Source age: {source_age} days")
    print(f"    Shell: M={shell_mass} Msun, v={shell_velocity} km/s, X={shell_xh}, Y={shell_yhe}, "
          f"composition={shell_composition}")
    if no_wind:
        print(f"    Wind: SKIPPED (no_wind=True)")
    elif density_matched_wind:
        print(f"    Wind: DENSITY-MATCHED (Mdot derived, v={wind_velocity} km/s, "
              f"CD sigmoid k={cd_sigmoid_k}, zones={cd_transition_zones})")
    else:
        print(f"    Wind: Mdot={wind_mdot} Msun/yr, v={wind_velocity} km/s")
    print(f"    Target days: {target_days}, Snapshot day: {snapshot_day}")
    print(f"    Output: {model_dir_abs}")

    script = scripts_dir / 'prepare_x_model.py'
    cmd = [
        sys.executable, str(script),
        str(source_dir),
        str(model_dir_abs),
        '--model-name', str(model_name),
        '--shell-mass', str(shell_mass),
        '--shell-velocity', str(shell_velocity),
        '--shell-xh', str(shell_xh),
        '--shell-yhe', str(shell_yhe),
        '--shell-composition', str(shell_composition),
        '--n-rho', str(n_rho),
        '--v0', str(v0_kms),
        '--n-T', str(n_T),
        '--T-interface', str(T_interface),
        '--T-max', str(T_max),
        '--t-start', str(t_start),
        '--ejecta-T-floor', str(ejecta_T_floor),
        '--ejecta-rho-min', str(ejecta_rho_min),
        '--ejecta-r-max-cm', str(ejecta_r_max_cm),
        '--wind-mdot', str(wind_mdot),
        '--wind-velocity', str(wind_velocity),
        '--wind-T', str(wind_T),
        '--se-transition-dex', str(se_transition_dex),
        '--target-days', str(target_days),
        '--snapshot-day', str(snapshot_day),
        '--nzon', str(nzon),
        '--source-age', str(source_age),
        '--explosion-energy', str(explosion_energy),
        '--aq', str(aq),
        '--fit-tau', str(fit_tau),
        '--meth', str(meth),
        '--temp-floor', str(temp_floor),
        '--floor-rad', str(floor_rad),
        '--hmin', str(hmin),
        '--hmax', str(hmax),
    ]
    if no_wind:
        cmd.append('--no-wind')
    if density_matched_wind:
        cmd.append('--density-matched-wind')
        cmd.extend(['--cd-transition-zones', str(cd_transition_zones),
                    '--cd-sigmoid-k', str(cd_sigmoid_k),
                    '--v-blend-zones', str(v_blend_zones)])
    if force_homologous_velocity:
        cmd.append('--force-homologous-velocity')

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"ERROR: X model profile generation failed")
        print(result.stdout)
        print(result.stderr)
        return False

    print(result.stdout)
    return True


def generate_namelist(config, scripts_dir, project_dir, config_path, debug=False):
    """[DEPRECATED] Generate HERACLES namelist from YAML config.

    HERACLES is deprecated; use the STELLA -> CMFGEN pipeline instead.
    Note: generate_heracles_namelist.py has been moved to scripts/deprecated/.
    """
    model_name = config.get('model_name', 'unknown')
    print(f"\n{'='*60}")
    print(f"Generating [DEPRECATED] HERACLES namelist for model: {model_name}")
    print(f"{'='*60}")
    
    script = scripts_dir / 'deprecated' / 'generate_heracles_namelist.py'
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
            # Match model name with underscore or hyphen delimiters to avoid partial matches
            # Matches: "X1_heating_*", "X1-debug_*", "*_X1_heating_*", or "*_X1_heating"
            # Does NOT match: "X1_heating2_*", "X1_heating_scat_*" (when looking for "X1_heating")
            if (job_name.startswith(f'{model_name}_') or
                job_name.startswith(f'{model_name}-') or
                f'_{model_name}_' in job_name or
                f'_{model_name}-' in job_name or
                job_name.endswith(f'_{model_name}') or
                job_name.endswith(f'-{model_name}')):
                jobs.append((job_id, job_name))
    
    return jobs


def cancel_existing_jobs(model_name, force=False, verbose=True, job_filter=None):
    """Cancel existing Slurm jobs for this model.

    Args:
        model_name: Model name to match
        force: If True, cancel without asking. If False, prompt user.
        verbose: If True, print detailed messages. If False, be silent.
        job_filter: Optional substring to filter job names (e.g. 'sedona' to only
            cancel SEDONA jobs, leaving CMFGEN/STELLA jobs running). Case-insensitive.

    Returns:
        Number of jobs cancelled
    """
    existing_jobs = check_existing_jobs(model_name)

    # Apply optional job name filter (e.g. only cancel sedona jobs)
    if job_filter:
        job_filter_lower = job_filter.lower()
        existing_jobs = [(jid, jn) for jid, jn in existing_jobs if job_filter_lower in jn.lower()]

    if not existing_jobs:
        return 0

    if verbose:
        label = f" (filter: '{job_filter}')" if job_filter else ""
        print(f"\nFound {len(existing_jobs)} existing job(s) for model '{model_name}'{label}:")
        for job_id, job_name in existing_jobs:
            print(f"  {job_id}: {job_name}")

    if not force:
        response = input(f"Cancel these jobs? (y/N): ").strip().lower()
        if response != 'y':
            if verbose:
                print("Skipping job cancellation.")
            return 0
    
    cancelled_count = 0
    for job_id, job_name in existing_jobs:
        try:
            result = subprocess.run(
                ['scancel', job_id],
                capture_output=True,
                text=True
            )
            if result.returncode == 0:
                if verbose:
                    print(f"  Cancelled job {job_id} ({job_name})")
                cancelled_count += 1
            else:
                if verbose:
                    print(f"  Failed to cancel job {job_id}: {result.stderr}")
        except Exception as e:
            if verbose:
                print(f"  Error cancelling job {job_id}: {e}")
    
    if cancelled_count > 0 and verbose:
        print(f"Cancelled {cancelled_count} job(s)")
        # Wait a moment for Slurm to process the cancellations
        import time
        time.sleep(2)
    
    return cancelled_count


def clean_previous_sedona_runs(config, project_dir, force_cancel=True, verbose=True):
    """Cancel running SEDONA jobs and delete previous SEDONA output for a model.

    This ONLY affects SEDONA jobs and output — CMFGEN, STELLA, and SIROCCO jobs
    are never cancelled or touched.

    Args:
        config: YAML configuration dictionary
        project_dir: Path to project directory
        force_cancel: If True, cancel without prompting.
        verbose: If True, print detailed messages.

    Returns:
        Number of SEDONA jobs cancelled.
    """
    model_name = config.get('model_name', 'unknown')
    output = config.get('output', {})
    model_dir = output.get('model_dir', '../models/${model_name}').replace('${model_name}', model_name)
    model_dir = (project_dir / model_dir).absolute()
    sedona_root = model_dir / 'sedona'

    cancelled = 0

    # 1. Cancel only SEDONA jobs (leave CMFGEN/STELLA/SIROCCO untouched)
    if force_cancel:
        cancelled = cancel_existing_jobs(model_name, force=True, verbose=verbose, job_filter='sedona')

    # 2. Delete previous SEDONA output directories and files
    import shutil
    import glob as _glob

    if sedona_root.exists():
        if verbose:
            print(f"\n  Cleaning previous SEDONA runs in: {sedona_root}")

        # Delete per-day subdirectories (e.g. 30.0days/, 40.0days/)
        for subdir in list(sedona_root.iterdir()):
            if subdir.is_dir():
                shutil.rmtree(subdir)
                if verbose:
                    print(f"    Removed directory: {subdir.name}")

        # Delete output files in the sedona root (spectra, plots, logs, slurm scripts)
        stale_patterns = ['spectrum_*.dat', 'spectrum_*.h5', 'spectrum_auto_plot.png',
                          'plt_*.dat', 'plt_*.h5', 'integrated_quantities.dat',
                          'sedona_*.out', 'sedona_*.err', 'sedona_*.slurm', 'spectrum.h5',
                          '*.mod', 'param.lua', 'plot_spectrum.py']
        removed = 0
        for pattern in stale_patterns:
            for f in _glob.glob(str(sedona_root / pattern)):
                os.remove(f)
                removed += 1
        if verbose and removed > 0:
            print(f"    Removed {removed} stale file(s) from sedona root")

        # Remove the sedona root directory itself if it's now empty
        try:
            sedona_root.rmdir()
            if verbose:
                print(f"    Removed empty sedona root directory")
        except OSError:
            pass  # Directory not empty (may have other files)

    return cancelled


def _read_hydro_section(hydro_path, section_name):
    """Read a named section from SN_HYDRO_DATA (e.g. 'Temperature', 'Velocity').
    Returns a 1D numpy array or None if the section is not found."""
    try:
        with open(hydro_path) as f:
            lines = f.readlines()
    except Exception:
        return None
    in_section = False
    data = []
    for line in lines:
        s = line.strip()
        if not s:
            continue
        if section_name.lower() in s.lower():
            in_section = True
            continue
        if in_section:
            if s.startswith('!') and 'mass' not in s.lower() and 'fraction' not in s.lower():
                in_section = False
                continue
            try:
                vals = [float(x) for x in s.split()]
                data.extend(vals)
            except ValueError:
                in_section = False
                continue
    if not data:
        return None
    return np.array(data)


def generate_grey_scl_fac_in(sn_hydro_path):
    """Generate GREY_SCL_FAC_IN from SN_HYDRO_DATA.

    Reads the STELLA T profile and Rosseland opacity from SN_HYDRO_DATA,
    computes the grey temperature (Eddington approximation) using T_eff
    derived from the inner boundary condition, and writes the ratio
    T_STELLA/T_grey as a function of log10(Tau_Ross).

    This file is read by CMFGEN's SCALE_GREY routine to scale the
    grey temperature distribution towards the STELLA hydrodynamic T,
    giving better initial temperature estimates for convergence.
    """
    import re
    sn_hydro_path = Path(sn_hydro_path)
    if not sn_hydro_path.exists():
        return False

    lines = sn_hydro_path.read_text().splitlines()
    sections = {}
    current = None
    data = []
    for ln in lines:
        s = ln.strip()
        if not s:
            continue
        if re.match(r'^[A-Z]', s) and not re.match(r'^[0-9.E+-]', s):
            if current and data:
                sections[current] = data
            current = s
            data = []
        elif s.startswith('!'):
            if current and data:
                sections[current] = data
            current = None
            data = []
        else:
            try:
                data.extend([float(x) for x in s.split()])
            except ValueError:
                if current and data:
                    sections[current] = data
                current = None
                data = []
    if current and data:
        sections[current] = data

    # Need Radius, Temperature, and Rosseland opacity
    r_key = next((k for k in sections if 'Radius' in k), None)
    t_key = next((k for k in sections if 'Temperature' in k), None)
    chi_key = next((k for k in sections if 'Rosseland' in k or 'opacty' in k.lower()), None)
    if not (r_key and t_key and chi_key):
        print("  GREY_SCL_FAC_IN: missing required sections in SN_HYDRO_DATA")
        return False

    R_cm = np.array(sections[r_key])
    # Detect unit: header says (10^10cm) but verify
    if R_cm[0] < 1e10:
        R_cm = R_cm * 1e10  # convert from 10^10 cm to cm
    T_K = np.array(sections[t_key])
    if T_K[0] < 1e3:
        T_K = T_K * 1e4  # convert from 10^4 K to K
    chi = np.array(sections[chi_key])
    if chi[0] < 1e-5:
        chi = chi * 1e-10  # convert from 10^{-10} cm^{-1} to cm^{-1}

    nd = len(R_cm)
    if nd < 3:
        print(f"  GREY_SCL_FAC_IN: too few depth points ({nd})")
        return False

    # Compute Tau_Ross (outer to inner, tau increases inward)
    tau = np.zeros(nd)
    for i in range(1, nd):
        dR = R_cm[i - 1] - R_cm[i]
        chi_mid = 0.5 * (chi[i - 1] + chi[i])
        tau[i] = tau[i - 1] + chi_mid * dR
    tau[0] = chi[0] * (R_cm[0] - R_cm[1]) * 0.5

    # Check for non-positive tau
    if tau.min() <= 0:
        print(f"  GREY_SCL_FAC_IN: non-positive tau detected (min={tau.min():.4e}), skipping")
        return False

    log_tau = np.log10(tau)

    # Check monotonicity (required by SCALE_GREY)
    for i in range(1, nd):
        if log_tau[i] <= log_tau[i - 1]:
            print(f"  GREY_SCL_FAC_IN: log_tau not monotonic at depth {i+1}, skipping")
            return False

    # T_eff from inner boundary: T_inner = T_eff * (0.75*(tau_inner + 2/3))^0.25
    T_inner = T_K[-1]
    tau_inner = tau[-1]
    T_eff = T_inner / (0.75 * (tau_inner + 2.0 / 3.0)) ** 0.25

    # Grey temperature (Eddington approximation)
    T_grey = T_eff * (0.75 * (tau + 2.0 / 3.0)) ** 0.25

    # Scale factor = T_STELLA / T_grey
    scale_fac = T_K / T_grey

    # Write GREY_SCL_FAC_IN
    out_path = sn_hydro_path.parent / 'GREY_SCL_FAC_IN'
    with open(out_path, 'w') as f:
        f.write('! GREY_SCL_FAC_IN generated from STELLA SN_HYDRO_DATA\n')
        f.write('! Scale factor = T_STELLA / T_grey (Eddington grey approximation)\n')
        f.write(f'! T_eff = {T_eff:.0f} K (from inner T={T_inner:.0f} K, tau={tau_inner:.1f})\n')
        f.write(f'! ND = {nd}\n')
        f.write(f'{nd}\n')
        for i in range(nd):
            f.write(f'{log_tau[i]:12.6f}  {scale_fac[i]:12.6f}\n')

    print(f"  GREY_SCL_FAC_IN written: {out_path} ({nd} points, "
          f"scale range {scale_fac.min():.3f}-{scale_fac.max():.3f}, "
          f"T_eff={T_eff:.0f} K)")
    return True


def prepare_cmfgen_directory(config, scripts_dir, project_dir, day=None):
    """Prepare CMFGEN directory by copying base files.

    Args:
        day: Optional snapshot day label. If provided, directory is
            cmfgen/<day>days for time-series runs.

    ── Diagnostic: files and scripts involved ───────────────────────────

    Copies the CMFGEN base directory (cmfgen_base/) to the model's run
    directory.  The base directory contains template files that are
    customized at runtime:

      cmfgen_base/MODEL_SPEC    — Atomic model specification (ND, NC, NP,
                                   NUM_BNDS, ISF entries for each ion).
                                   ND/NP are updated by prepare_VADAT.py
                                   at job start to match SN_HYDRO_DATA.
      cmfgen_base/VADAT_INIT    — VADAT template (radiative transfer controls).
      cmfgen_base/RDINR        — Input hydrodynamics reader configuration.
      cmfgen_base/T_IN         — Temperature initialization profile.
      cmfgen_base/IN_ITS       — Iteration control parameters.
      cmfgen_base/batch.sh     — CMFGEN batch script (executes cmfgen_dev.exe
                                   for static_snapshot model_type).
      cmfgen_base/submit_cmfgen.sh — CMFGEN submission wrapper (calls
                                   prepare_VADAT.py, prepare_IN_INTS.py,
                                   then submits SLURM jobs).

    The model_type config key selects the CMFGEN executable:
      'static_snapshot' -> cmfgen_dev.exe  (Dessart+2015 direct treatment)
      'lte'             -> main_lte.exe    (LTE approximation)

    Existing SN_HYDRO_DATA and plots/ are preserved during directory
    refresh to avoid race conditions with running conversion jobs.
    ─────────────────────────────────────────────────────────────────────
    """
    model_name = config.get('model_name', 'unknown')
    day_label = f"{day:.2f}days" if day is not None else None
    print(f"\n{'='*60}")
    if day_label:
        print(f"Preparing CMFGEN directory for model: {model_name} day {day}")
    else:
        print(f"Preparing CMFGEN directory for model: {model_name}")
    print(f"{'='*60}")

    # Extract paths
    output = config.get('output', {})
    model_dir = output.get('model_dir', '../models/${model_name}').replace('${model_name}', model_name)
    model_dir = (project_dir / model_dir).absolute()

    cmfgen = config.get('cmfgen', {})
    base_dir = cmfgen.get('base_dir', '../cmfgen_base')
    base_dir = (project_dir / base_dir).absolute()

    # model_type selects the CMFGEN executable and VADAT configuration used in batch.sh
    MODEL_TYPE_EXE = {
        'static_snapshot': 'cmfgen_dev.exe',
        'lte':            'main_lte.exe',
    }
    model_type = cmfgen.get('model_type', 'static_snapshot')
    exe_name = MODEL_TYPE_EXE.get(model_type)
    if exe_name is None:
        print(f"ERROR: Unknown model_type '{model_type}'. Supported: {list(MODEL_TYPE_EXE)}")
        return False

    cmfgen_dir = model_dir / 'cmfgen' / day_label if day_label else model_dir / 'cmfgen'

    # Copy base directory to cmfgen directory
    if not base_dir.exists():
        print(f"ERROR: Base directory not found: {base_dir}")
        return False

    if cmfgen_dir.exists():
        import shutil, tempfile
        # Preserve existing SN_HYDRO_DATA if it exists (avoid race condition where
        # prepare wipes the file before the CMFGEN job can read it).
        sn_hydro_existing = cmfgen_dir / 'SN_HYDRO_DATA'
        sn_hydro_backup = None
        if sn_hydro_existing.exists():
            sn_hydro_backup = Path(tempfile.mkdtemp()) / 'SN_HYDRO_DATA'
            shutil.copy2(sn_hydro_existing, sn_hydro_backup)
            print(f"  Preserving existing SN_HYDRO_DATA during directory refresh")

        # Preserve existing plots/ subdirectory (SN_HYDRO_DATA diagnostic plots
        # from convert_stella_swd_to_cmfgen.py). These are regenerated when the
        # conversion runs, but if the conversion fails or is skipped (e.g.
        # SN_HYDRO_DATA already exists), the plots from the previous run would
        # be lost permanently without this preservation.
        plots_existing = cmfgen_dir / 'plots'
        plots_backup = None
        if plots_existing.exists() and any(plots_existing.iterdir()):
            plots_backup = Path(tempfile.mkdtemp()) / 'plots'
            shutil.copytree(plots_existing, plots_backup)
            print(f"  Preserving existing plots/ during directory refresh")

        # Preserve EDDFACTOR + EDDFACTOR_INFO for USE_FIXED_J initialization.
        # When use_fixed_j=True in config, CMFGEN requires a valid EDDFACTOR
        # file from a previous run. The orchestrator wipes the directory on
        # resubmission, so we must preserve these files.
        eddfactor_backup = None
        eddfactor_info_backup = None
        cmfgen_cfg = config.get('cmfgen', {})
        preserve_edd = cmfgen_cfg.get('vadat_overrides', {}).get('USE_FIXED_J', 'F').upper().startswith('T')
        if preserve_edd:
            edd_existing = cmfgen_dir / 'EDDFACTOR'
            edd_info_existing = cmfgen_dir / 'EDDFACTOR_INFO'
            if edd_existing.exists() and edd_info_existing.exists():
                tmpdir = Path(tempfile.mkdtemp())
                eddfactor_backup = tmpdir / 'EDDFACTOR'
                eddfactor_info_backup = tmpdir / 'EDDFACTOR_INFO'
                shutil.copy2(edd_existing, eddfactor_backup)
                shutil.copy2(edd_info_existing, eddfactor_info_backup)
                print(f"  Preserving existing EDDFACTOR + EDDFACTOR_INFO for USE_FIXED_J initialization")

        # Preserve restart/output files from previous CMFGEN run.
        # These are expensive to regenerate and needed for restart, analysis,
        # or downstream codes. The directory wipe would otherwise destroy them.
        restart_files = [
            'OUTGEN',              # CMFGEN output populations (restart)
            'RVTJ',                # R, V, T, J structure (restart + OBS input)
            'MODEL',               # Current model structure (restart)
            'CUR_MODEL_DATA',      # Current model data (restart)
            'SPECIES_MASSES',      # Species mass inventory (diagnostic)
            'OBSFLUX',             # Observer-frame flux (OBS output)
            'RS_XRAY_FLUXES',      # X-ray fluxes (diagnostic)
            'SN_HYDRO_FOR_NEXT_MODEL',  # Hydro for next time step
            'input_hydro.dat',     # Input hydro (diagnostic)
            'IT_SPECIFIER',        # Iteration specifier (restart)
        ]
        restart_backups = {}
        for fname in restart_files:
            fpath = cmfgen_dir / fname
            if fpath.exists():
                tmpdir = Path(tempfile.mkdtemp())
                backup = tmpdir / fname
                shutil.copy2(fpath, backup)
                restart_backups[fname] = (backup, tmpdir)
                print(f"  Preserving existing {fname} during directory refresh")

        print(f"Removing existing cmfgen directory: {cmfgen_dir}")
        shutil.rmtree(cmfgen_dir)

        print(f"Copying base directory from {base_dir} to {cmfgen_dir}...")
        shutil.copytree(base_dir, cmfgen_dir)

        # Restore SN_HYDRO_DATA if it was preserved
        if sn_hydro_backup is not None:
            shutil.copy2(sn_hydro_backup, cmfgen_dir / 'SN_HYDRO_DATA')
            shutil.rmtree(sn_hydro_backup.parent)

        # Restore plots/ if it was preserved
        if plots_backup is not None:
            shutil.copytree(plots_backup, cmfgen_dir / 'plots')
            shutil.rmtree(plots_backup.parent)

        # Restore EDDFACTOR if it was preserved
        if eddfactor_backup is not None:
            shutil.copy2(eddfactor_backup, cmfgen_dir / 'EDDFACTOR')
            shutil.copy2(eddfactor_info_backup, cmfgen_dir / 'EDDFACTOR_INFO')
            shutil.rmtree(eddfactor_backup.parent)

        # Restore restart/output files
        for fname, (backup, tmpdir) in restart_backups.items():
            shutil.copy2(backup, cmfgen_dir / fname)
            shutil.rmtree(tmpdir)
    else:
        import shutil
        print(f"Copying base directory from {base_dir} to {cmfgen_dir}...")
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


def prepare_sedona_directory(config, scripts_dir, project_dir, use_stella=False, day=None):
    """Prepare SEDONA directory by copying base files and creating param.lua.

    Args:
        config: YAML configuration dictionary
        scripts_dir: Path to scripts directory
        project_dir: Path to project directory
        use_stella: If True, convert from STELLA .swd file instead of CMFGEN SN_HYDRO_DATA
        day: Optional snapshot day label. If provided, directory is
            sedona/<day>days for time-series runs.
    """
    model_name = config.get('model_name', 'unknown')
    day_label = f"{day:.2f}days" if day is not None else None
    print(f"\n{'='*60}")
    if day_label:
        print(f"Preparing SEDONA directory for model: {model_name} day {day}")
    else:
        print(f"Preparing SEDONA directory for model: {model_name}")
    print(f"{'='*60}")

    # Extract paths
    output = config.get('output', {})
    model_dir = output.get('model_dir', '../models/${model_name}').replace('${model_name}', model_name)
    model_dir = (project_dir / model_dir).absolute()

    sedona_dir = model_dir / 'sedona' / day_label if day_label else model_dir / 'sedona'
    sedona_dir.mkdir(parents=True, exist_ok=True)

    # Clean up stale output files from previous SEDONA runs to ensure a fresh start.
    # SEDONA appends to spectrum_*.dat/plt_*.dat and accumulates integrated_quantities.dat,
    # so stale files from prior runs would corrupt or contaminate new output.
    stale_patterns = ['spectrum_*.dat', 'spectrum_*.h5', 'spectrum_auto_plot.png',
                      'plt_*.dat', 'plt_*.h5', 'integrated_quantities.dat',
                      'sedona_*.out', 'sedona_*.err', 'spectrum.h5']
    removed_count = 0
    import glob as _glob
    for pattern in stale_patterns:
        for stale_file in _glob.glob(str(sedona_dir / pattern)):
            os.remove(stale_file)
            removed_count += 1
    if removed_count > 0:
        print(f"  Cleaned {removed_count} stale output file(s) from previous run")

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
    
    # Convert to SEDONA .mod format (either from STELLA .swd or CMFGEN SN_HYDRO_DATA)
    mod_file = sedona_dir / f'{model_name}.mod'

    # Read SEDONA settings
    sedona_cfg = config.get('sedona', {})
    use_homologous = sedona_cfg.get('use_homologous_expansion', False)
    use_nlte = sedona_cfg.get('use_nlte', not use_homologous)  # default: off for homologous, on for static snapshot
    nzon_target = sedona_cfg.get('nzon_target', 300)
    transport_steady_iterate = sedona_cfg.get('transport_steady_iterate', 25)
    core_n_emit = sedona_cfg.get('core_n_emit', 0)
    particles_n_emit_radioactive = sedona_cfg.get('particles_n_emit_radioactive', 5e4)
    particles_n_emit_thermal = sedona_cfg.get('particles_n_emit_thermal', 5e8)  # 5e8 thermal photons (high count for fine spectrum grid + low escape fraction)
    transport_nu_grid_spacing = float(sedona_cfg.get('transport_nu_grid_spacing', 0.005))
    spectrum_nu_grid_spacing = float(sedona_cfg.get('spectrum_nu_grid_spacing', 0.0005))
    # Spectrum frequency range: narrow to optical for maximum resolution where it matters.
    # Defaults: 2.5e14 Hz (12000 A) to 1.2e15 Hz (2500 A) — covers optical + near-UV.
    # Full range would be {1.0e14, 2.0e16} (150-30000 A) but wastes bins on UV/IR.
    spectrum_nu_min = float(sedona_cfg.get('spectrum_nu_min', 2.5e14))
    spectrum_nu_max = float(sedona_cfg.get('spectrum_nu_max', 1.2e15))
    transport_rad_eq = sedona_cfg.get('transport_radiative_equilibrium', 0)
    particles_max_total = sedona_cfg.get('particles_max_total', 2e9)  # Override SEDONA default max for high-photon runs
    particles_last_iter_pump = sedona_cfg.get('particles_last_iter_pump', 20)  # Radiation field averaging window (more iters = more photons in spectrum)
    transport_use_ddmc = sedona_cfg.get('transport_use_ddmc', 0)
    transport_ddmc_tau_threshold = sedona_cfg.get('transport_ddmc_tau_threshold', 1000.0)
    # Line opacity configuration
    # Sobolev (line_expansion) uses local |dv/dr| from source code fix for non-monotonic velocity.
    # Direct Voigt (bound_bound) causes photon trapping in optically thick ejecta (tau_es > 10).
    # Fuzz expansion provides iron-group line blanketing (uses 1/t_exp Sobolev approximation).
    opacity_line_expansion = sedona_cfg.get('opacity_line_expansion', 1)  # Sobolev ON with local |dv/dr| (default for non-monotonic velocity)
    opacity_bound_bound = sedona_cfg.get('opacity_bound_bound', 0)  # Direct Voigt OFF (causes photon trapping)
    opacity_fuzz_expansion = sedona_cfg.get('opacity_fuzz_expansion', 1 if use_homologous else 0)  # Default: ON for homologous, OFF for static
    line_velocity_width_kms = sedona_cfg.get('line_velocity_width_kms', 200.0)  # Turbulent broadening (km/s); 0=thermal
    line_velocity_width_cgs = line_velocity_width_kms * 1.0e5  # km/s -> cm/s

    # Conversion source selection (mirrors superlite.input_source):
    #   'swd'           : direct STELLA .swd -> .mod via
    #                     convert_stella_swd_to_sedona.py (one step, no CMFGEN
    #                     dependency, native STELLA grid + kappa, full 15-element
    #                     composition).  Skips the CMFGEN SN_HYDRO_DATA step.
    #   'sn_hydro_data' : (default) CMFGEN SN_HYDRO_DATA -> .mod.  When
    #                     use_stella=True this runs the two-step pipeline
    #                     (.swd -> SN_HYDRO_DATA -> .mod).
    input_source = sedona_cfg.get('input_source', 'sn_hydro_data')

    # Prefer the CMFGEN tau-trimmed SN_HYDRO_DATA if it exists (matches the
    # spectrum calculation and avoids the STELLA .swd path issue).
    # In time-series mode, SN_HYDRO_DATA lives at cmfgen/<day>days/SN_HYDRO_DATA.
    # NOTE: skipped when input_source='swd' (direct .swd path requested explicitly).
    if day_label:
        sn_hydro = model_dir / 'cmfgen' / day_label / 'SN_HYDRO_DATA'
    else:
        raise ValueError("day is required — SN_HYDRO_DATA must live in cmfgen/<day>days/, not cmfgen/ root")
    if input_source != 'swd' and sn_hydro.exists() and use_stella:
        print(f"  Found CMFGEN SN_HYDRO_DATA; using it for SEDONA instead of STELLA .swd")
        use_stella = False

    # Read snapshot days from config
    cmfgen_conv = config.get('cmfgen_conversion', {})
    snapshot_days_config = cmfgen_conv.get('snapshot_days', 40)
    stella_cfg = config.get('stella', {})
    target_days = stella_cfg.get('target_days', 60.0)

    if input_source == 'swd':
        # --- Direct STELLA .swd -> SEDONA .mod (one step) ---
        # Uses convert_stella_swd_to_sedona.py: reads the .swd snapshot
        # directly, tau-truncates the optically-thick inner region, preserves
        # the native STELLA grid + Rosseland kappa, and writes a SEDONA .mod
        # with the full element set present in the model.  No CMFGEN
        # SN_HYDRO_DATA is required, so SEDONA can run even when CMFGEN is
        # skipped.
        stella_dir = model_dir / 'stella'
        swd_file = stella_dir / 'strad' / f'{model_name}.swd'
        if not swd_file.exists():
            swd_file = stella_dir / f'{model_name}.swd'

        convert_script = scripts_dir / 'convert_stella_swd_to_sedona.py'
        if not convert_script.exists():
            print(f"  ERROR: SEDONA conversion script not found at {convert_script}")
            return False

        if swd_file.exists():
            print(f"Converting STELLA .swd -> SEDONA .mod (direct, one step)...")
            if use_homologous:
                print(f"  Homologous expansion enabled: correct time + radius in .mod header")
            convert_day = day if day is not None else target_days
            if mod_file.exists():
                print(f"  Removing existing {mod_file.name} to ensure fresh conversion")
                mod_file.unlink()
            cmd = [sys.executable, str(convert_script), str(swd_file), str(mod_file),
                   '--target-days', str(convert_day),
                   '--time-days', str(convert_day),
                   '--model-name', str(model_name)]
            # tau-truncation threshold (default 100; matches CMFGEN tau_max ceiling)
            sedona_tau_max = sedona_cfg.get('tau_max', 100.0)
            cmd.extend(['--tau-max', str(sedona_tau_max)])
            # temperature floor
            sedona_t_floor = sedona_cfg.get('t_floor_K', 2000.0)
            if sedona_t_floor > 0:
                cmd.extend(['--t-floor', str(sedona_t_floor)])
            # optional density floor for near-vacuum outer zones
            sedona_rho_floor = sedona_cfg.get('rho_floor', 0.0)
            if sedona_rho_floor > 0:
                cmd.extend(['--rho-floor', str(sedona_rho_floor)])
            # optional downsampling target
            if nzon_target > 0:
                cmd.extend(['--nzon-target', str(nzon_target)])
            # optional explicit species list
            sedona_species = sedona_cfg.get('species', None)
            if sedona_species:
                cmd.extend(['--species', str(sedona_species)])
            # use radiation temperature instead of gas temperature
            if sedona_cfg.get('use_T_rad', False):
                cmd.append('--use-T-rad')
            if use_homologous:
                cmd.append('--homologous')

            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode == 0:
                print(f"  Successfully created {model_name}.mod")
                print(result.stdout)
            else:
                print(f"  ERROR: STELLA .swd -> SEDONA .mod conversion failed")
                print(result.stdout)
                print(result.stderr)
                return False
            # Core properties are read from the .mod file below (use_stella path)
            use_stella = True
        else:
            print(f"  WARNING: STELLA .swd not found at {swd_file}")
            print(f"    The conversion will be attempted at job runtime as a fallback.")
            use_stella = True
    elif use_stella:
        # Convert from STELLA .swd file via SN_HYDRO_DATA (two-step pipeline):
        #   .swd -> SN_HYDRO_DATA (convert_stella_swd_to_cmfgen.py)
        #   SN_HYDRO_DATA -> .mod (convert_cmfgen_to_sedona_v2.py, correct rmin/texp)
        stella_dir = model_dir / 'stella'
        swd_file = stella_dir / 'strad' / f'{model_name}.swd'
        if not swd_file.exists():
            swd_file = stella_dir / f'{model_name}.swd'

        if swd_file.exists():
            print(f"Converting STELLA .swd -> SN_HYDRO_DATA -> SEDONA .mod format...")
            if use_homologous:
                print(f"  Homologous expansion enabled: correct time + radius in .mod header")

            # Step 1: .swd -> SN_HYDRO_DATA
            convert_day = day if day is not None else target_days
            if day_label:
                sn_hydro = model_dir / 'cmfgen' / day_label / 'SN_HYDRO_DATA'
            else:
                raise ValueError("day is required — cannot write SN_HYDRO_DATA to cmfgen/ root")
            cmfgen_convert = scripts_dir / 'convert_stella_swd_to_cmfgen.py'
            if not cmfgen_convert.exists():
                print(f"  ERROR: Conversion script not found at {cmfgen_convert}")
                return False
            cmfgen_cmd = [
                sys.executable, str(cmfgen_convert), str(swd_file), str(sn_hydro),
                '--target-days', str(convert_day),
                '--model-name', str(model_name),
            ]
            result = subprocess.run(cmfgen_cmd, capture_output=True, text=True)
            if result.returncode != 0:
                print(f"  ERROR: .swd -> SN_HYDRO_DATA conversion failed")
                print(result.stdout)
                print(result.stderr)
                return False

            # Step 2: SN_HYDRO_DATA -> .mod via v2
            sedona_convert = project_dir.parent / 'sedona' / 'convert_cmfgen_to_sedona_v2.py'
            if not sedona_convert.exists():
                print(f"  ERROR: Conversion script not found at {sedona_convert}")
                return False
            if mod_file.exists():
                print(f"  Removing existing {mod_file.name} to ensure fresh conversion")
                mod_file.unlink()
            cmd = [sys.executable, str(sedona_convert), str(sn_hydro), str(mod_file)]
            if use_homologous:
                cmd.extend(['--homologous', '--time-days', str(convert_day)])

            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode == 0:
                print(f"  Successfully created {model_name}.mod")
                print(result.stdout)
            else:
                print(f"  ERROR: SN_HYDRO_DATA -> .mod conversion failed")
                print(result.stdout)
                print(result.stderr)
                return False
        else:
            print(f"  ERROR: STELLA .swd file not found at {swd_file}")
            print(f"    SEDONA requires STELLA .swd file to create the .mod file")
            return False
    else:
        # Convert from CMFGEN SN_HYDRO_DATA (original behavior)
        if day_label:
            sn_hydro = model_dir / 'cmfgen' / day_label / 'SN_HYDRO_DATA'
        else:
            raise ValueError("day is required — SN_HYDRO_DATA must live in cmfgen/<day>days/")

        if sn_hydro.exists():
            print(f"Converting CMFGEN SN_HYDRO_DATA to SEDONA .mod format...")
            if use_homologous:
                print(f"  Homologous expansion enabled: correct time + radius in .mod header")
            # Use the conversion script from sedona
            sedona_convert = project_dir.parent / 'sedona' / 'convert_cmfgen_to_sedona_v2.py'
            if sedona_convert.exists():
                # Remove existing .mod file to ensure fresh conversion
                if mod_file.exists():
                    print(f"  Removing existing {mod_file.name} to ensure fresh conversion")
                    mod_file.unlink()

                cmd = [sys.executable, str(sedona_convert), str(sn_hydro), str(mod_file)]
                convert_day = day if day is not None else target_days
                if use_homologous:
                    cmd.extend(['--homologous', '--time-days', str(convert_day)])

                result = subprocess.run(cmd, capture_output=True, text=True)
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
    rmin_from_mod = None
    # Radioactive species (Z values) that produce decay photons
    RADIOACTIVE_Z = {28, 27, 24, 22, 20, 18, 16}  # Ni, Co, Cr, Ti, Ca, Ar, S
    if mod_file.exists():
        with open(mod_file) as f:
            mod_lines = f.readlines()
        if len(mod_lines) >= 4:
            # Read rmin from .mod header (line 2: "nzon  rmin  texp  n_elems")
            header_parts = mod_lines[1].strip().split()
            if len(header_parts) >= 2:
                try:
                    rmin_from_mod = float(header_parts[1])
                except ValueError:
                    pass
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

    # Disable radioactive photon emission if no radioactive species are present in .mod
    radioactive_in_mod = set(mod_species_z) & RADIOACTIVE_Z
    if not radioactive_in_mod and float(particles_n_emit_radioactive) > 0:
        print(f"  No radioactive species (Z in {sorted(RADIOACTIVE_Z)}) found in .mod")
        print(f"    Disabling particles_n_emit_radioactive (was {float(particles_n_emit_radioactive):.0e})")
        particles_n_emit_radioactive = 0
    
    # Get snapshot days from config (different for STELLA vs CMFGEN)
    # But prefer the per-day 'day' parameter when available (time_series mode)
    if use_stella:
        stella_cfg = config.get('stella', {})
        snapshot_days_config = stella_cfg.get('target_days', 60.0)
    else:
        cmfgen_conv = config.get('cmfgen_conversion', {})
        snapshot_days_config = cmfgen_conv.get('snapshot_days', 40)
    # Override with the actual snapshot day when running time series
    if day is not None:
        snapshot_days_config = day
    
    # Read inner boundary info from hydro model (STELLA .swd or CMFGEN SN_HYDRO_DATA)
    r_inner = 6.0e14      # default fallback
    T_inner = 5000.0      # default fallback (K)
    luminosity = 8.32e42  # default fallback (erg/s)
    snapshot_days = snapshot_days_config  # fallback
    
    if use_stella:
        # Read core properties from the .mod file (already converted from STELLA .swd)
        # The .mod file header has rmin, and the first data line has r[0], v, rho, T[0]
        mod_file_sedona = sedona_dir / f'{model_name}.mod'
        if mod_file_sedona.exists():
            try:
                with open(mod_file_sedona) as f:
                    f.readline()  # "1D_sphere standard"
                    header = f.readline().split()
                    rmin = float(header[1])
                    t_sec = float(header[2])
                    f.readline()  # species line
                    first_data = f.readline().split()
                    r0 = float(first_data[0])
                    T0 = float(first_data[3])

                r_inner = rmin
                T_inner = T0
                # Luminosity from Stefan-Boltzmann: L = 4*pi*r^2*sigma*T^4
                sigma_sb = 5.6704e-5
                luminosity = 4 * 3.14159265 * r0**2 * sigma_sb * T0**4
                snapshot_days = t_sec / 86400.0

                print(f"  Core properties from .mod file:")
                print(f"    rmin (core_radius): {r_inner:.6e} cm")
                print(f"    T[0] (core_temperature): {T_inner:.0f} K")
                print(f"    L (4*pi*r^2*sigma*T^4): {luminosity:.3e} erg/s")
                print(f"    snapshot_days: {snapshot_days:.2f}")
            except Exception as e:
                print(f"  WARNING: Failed to read .mod file for core properties: {e}")
                print(f"    Using defaults: R={r_inner:.3e} cm, T={T_inner:.0f} K")
        else:
            print(f"  WARNING: .mod file not found at {mod_file_sedona}")
            print(f"    Using defaults: R={r_inner:.3e} cm, T={T_inner:.0f} K")
    else:
        # Read from CMFGEN SN_HYDRO_DATA (original behavior)
        if day_label:
            sn_hydro = model_dir / 'cmfgen' / day_label / 'SN_HYDRO_DATA'
        else:
            raise ValueError("day is required — SN_HYDRO_DATA must live in cmfgen/<day>days/")

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

        # Ensure core_radius is within the grid bounds:
        #   >= rmin (grid inner boundary, from .mod header) so particles are on the grid
        #   < r_first_zone (outer radius of first zone) so particles are in the first zone
        if r_first_zone is not None and r_inner >= r_first_zone:
            core_radius = r_first_zone * 0.9999
            print(f"  Adjusted core_radius from {r_inner:.6e} to {core_radius:.6e} cm")
            print(f"    (first zone radius = {r_first_zone:.6e} cm)")
        elif rmin_from_mod is not None and r_inner < rmin_from_mod:
            core_radius = rmin_from_mod
            print(f"  Adjusted core_radius from {r_inner:.6e} to {core_radius:.6e} cm")
            print(f"    (rmin = {rmin_from_mod:.6e} cm; was inside grid inner boundary)")
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

        # Configure homologous expansion if requested
        if use_homologous:
            text = text.replace('hydro_module = "none"', 'hydro_module = "homologous"')
            # Ensure data_fuzzline_file is set (fuzz lines for iron-group blanketing)
            if 'data_fuzzline_file' not in text:
                text = text.replace(
                    'data_atomic_file   = sedona_home.."/data/cmfgen_levelcap100.hdf5"',
                    'data_atomic_file   = sedona_home.."/data/cmfgen_levelcap100.hdf5"\ndata_fuzzline_file = sedona_home.."/data/kurucz_cd23_fuzz.hdf5"  -- Iron-group line blanketing'
                )

        # Toggle SEDONA NLTE solver based on config (robust to template state)
        # GasState::set_atoms_in_nlte() forces use_nlte_=1 when the list is non-empty,
        # and use_nlte_=0 when empty — regardless of opacity_use_nlte. So the atoms
        # list is what actually controls NLTE. See GasState.cpp:109-112.
        if use_nlte:
            text = text.replace('opacity_use_nlte             = 0', 'opacity_use_nlte             = 1')
            # Populate the NLTE atoms list with species present in the .mod file.
            # Use mod_species_z if available (parsed from .mod), otherwise default to
            # the common SN species: H(1), He(2), O(8), Si(14).
            if mod_species_z:
                nlte_z = [z for z in mod_species_z if z in (1, 2, 8, 14)]
            else:
                nlte_z = [1, 2, 8, 14]
            if nlte_z:
                nlte_list_str = ", ".join(str(z) for z in nlte_z)
                text = re.sub(r'opacity_atoms_in_nlte\s*=\s*\{[^}]*\}',
                              f'opacity_atoms_in_nlte        = {{{nlte_list_str}}}  -- NLTE species: H, He, O, Si (from .mod)',
                              text)
        else:
            text = text.replace('opacity_use_nlte             = 1', 'opacity_use_nlte             = 0')
            # Also clear the atoms_in_nlte list — GasState::set_atoms_in_nlte() forces
            # use_nlte_=1 whenever the list is non-empty, ignoring opacity_use_nlte=0.
            # See GasState.cpp:109-112. Empty list {} is the only way to truly disable NLTE.
            text = re.sub(r'opacity_atoms_in_nlte\s*=\s*\{[^}]*\}',
                          'opacity_atoms_in_nlte        = {}  -- empty: NLTE fully disabled (GasState forces use_nlte=1 if list non-empty)',
                          text)

        # Configure line opacity: Sobolev (expansion) vs direct Voigt (bound-bound)
        # Sobolev with local |dv/dr| (source code fix) handles non-monotonic velocity correctly.
        # Direct Voigt causes photon trapping in optically thick ejecta (tau_es > 10).
        # Fuzz expansion provides iron-group line blanketing (uses 1/t_exp Sobolev approx).
        text = re.sub(r'opacity_line_expansion\s*=\s*\d+',
                      f'opacity_line_expansion       = {opacity_line_expansion}  -- {"Sobolev ON with local |dv/dr|" if opacity_line_expansion else "Sobolev OFF"}', text)
        text = re.sub(r'opacity_bound_bound\s*=\s*\d+',
                      f'opacity_bound_bound          = {opacity_bound_bound}  -- {"Direct Voigt ON (P-Cygni features)" if opacity_bound_bound else "Direct Voigt OFF"}', text)
        text = re.sub(r'opacity_fuzz_expansion\s*=\s*\d+',
                      f'opacity_fuzz_expansion       = {opacity_fuzz_expansion}  -- {"Fuzz line blanketing ON" if opacity_fuzz_expansion else "Fuzz OFF"}', text)
        # Line velocity width (turbulent broadening for direct Voigt lines)
        text = re.sub(r'line_velocity_width\s*=\s*[\d.eE+-]+',
                      f'line_velocity_width          = {line_velocity_width_cgs:.1e}  -- {line_velocity_width_kms:.0f} km/s turbulent broadening', text)

        # Configure resolution and convergence settings
        # Use regex for transport_steady_iterate to avoid partial-match bugs
        # (e.g. replacing "10" inside "100" would give "1000" instead of "100")
        text = re.sub(r'transport_steady_iterate\s*=\s*\d+',
                      f'transport_steady_iterate         = {transport_steady_iterate}', text)
        text = re.sub(r'core_n_emit\s*=\s*[\d.eE+-]+',
                      f'core_n_emit      = {float(core_n_emit):.0e}', text)
        text = re.sub(r'particles_n_emit_radioactive\s*=\s*[\d.eE+-]+',
                      f'particles_n_emit_radioactive = {float(particles_n_emit_radioactive):.0e}', text)
        text = re.sub(r'particles_n_emit_thermal\s*=\s*[\d.eE+-]+',
                      f'particles_n_emit_thermal     = {float(particles_n_emit_thermal):.0e}', text)
        # Particles_last_iter_pump: radiation field averaging window (lower = less memory, higher = smoother)
        text = re.sub(r'particles_last_iter_pump\s*=\s*\d+',
                      f'particles_last_iter_pump     = {int(particles_last_iter_pump)}', text)

        # Override particles_max_total if specified (SEDONA default is 1e7, too low for high-photon runs)
        if particles_max_total is not None:
            # Match any existing particles_max_total line and replace its value
            text = re.sub(r'particles_max_total\s*=\s*[\d.eE+-]+',
                          f'particles_max_total = {float(particles_max_total):.0e}', text)
            # If no existing line, add one after particles_n_emit_thermal
            if 'particles_max_total' not in text:
                text = text.replace(
                    f'particles_n_emit_thermal     = {float(particles_n_emit_thermal):.0e}',
                    f'particles_n_emit_thermal     = {float(particles_n_emit_thermal):.0e}\nparticles_max_total = {float(particles_max_total):.0e}  -- Override default max (1e7) for high-photon runs')

        # Configure radiative equilibrium (paper: hold T fixed from hydro)
        text = text.replace('transport_radiative_equilibrium  = 0', f'transport_radiative_equilibrium  = {transport_rad_eq}')

        # Configure DDMC (Discrete Diffusion Monte Carlo) for high optical depth
        # Update DDMC settings in param.lua (template has defaults: use_ddmc=2, tau_threshold=3.0)
        text = re.sub(r'transport_use_ddmc\s*=\s*\d+', f'transport_use_ddmc               = {transport_use_ddmc}', text)
        text = re.sub(r'transport_ddmc_tau_threshold\s*=\s*[\d.]+', f'transport_ddmc_tau_threshold     = {transport_ddmc_tau_threshold}', text)

        # Update frequency grid spacing for transport and spectrum independently.
        # transport_nu_grid: controls opacity binning and DDMC matrix size (coarse = faster)
        #   — range configurable via transport_nu_min/transport_nu_max in YAML
        #   — default: full range (1e14-2e16 Hz = 150-30000 A) for correct UV photon transport
        #   — for validation: can narrow to optical range to reduce frequency bins
        transport_nu_min = float(sedona_cfg.get('transport_nu_min', 1.0e14))
        transport_nu_max = float(sedona_cfg.get('transport_nu_max', 2.0e16))
        # spectrum_nu_grid: controls output spectrum binning (fine = better spectral resolution)
        #   — narrowed to optical range (spectrum_nu_min to spectrum_nu_max) to concentrate
        #     all bins where the science is. With q=0.0005 over 2.5e14-1.2e15 Hz (2500-12000 A),
        #     this gives ~3100 bins at ~3 A/bin (~150 km/s velocity resolution).
        text = text.replace('transport_nu_grid  = {1.0e14, 2.0e16, 0.001, 1}',
                          f'transport_nu_grid  = {{{transport_nu_min:.1e}, {transport_nu_max:.1e}, {transport_nu_grid_spacing}, 1}}')
        # Replace spectrum grid: support both the template default and any pre-existing value
        text = re.sub(r'spectrum_nu_grid\s*=\s*\{[^}]*\}',
                      f'spectrum_nu_grid   = {{{spectrum_nu_min:.1e}, {spectrum_nu_max:.1e}, {spectrum_nu_grid_spacing}, 1}}  -- Spectrum grid: {spectrum_nu_min:.1e}-{spectrum_nu_max:.1e} Hz, q={spectrum_nu_grid_spacing}',
                      text)

        if use_homologous:
            print(f"  Configured param.lua for homologous expansion:")
            print(f"    hydro_module = 'homologous'")
            print(f"    opacity_fuzz_expansion = {opacity_fuzz_expansion}")
            print(f"    opacity_use_nlte = 0")
        else:
            print(f"  Configured param.lua for static snapshot mode (paper Sect. 2.4):")
            print(f"    hydro_module = 'none' (no homologous expansion)")
            print(f"    opacity_use_nlte = {int(use_nlte)} ({'NLTE' if use_nlte else 'LTE (Saha-Boltzmann)'})")
            print(f"    opacity_line_expansion = {opacity_line_expansion} ({'Sobolev ON with local |dv/dr|' if opacity_line_expansion else 'Sobolev OFF'})")
            print(f"    opacity_bound_bound = {opacity_bound_bound} ({'Direct Voigt ON (P-Cygni)' if opacity_bound_bound else 'OFF'})")
            print(f"    opacity_fuzz_expansion = {opacity_fuzz_expansion}")
            print(f"    transport_use_ddmc = {transport_use_ddmc}")
            if opacity_bound_bound:
                print(f"    line_velocity_width = {line_velocity_width_kms:.0f} km/s (turbulent broadening for Voigt lines)")
            print(f"    transport_radiative_equilibrium = {transport_rad_eq} (T held from hydro)")
            print(f"    core_n_emit = {float(core_n_emit):.0e} (no artificial core)")
            print(f"    particles_n_emit_thermal = {float(particles_n_emit_thermal):.0e} (thermal photons from shock-heated zones)")
            print(f"    particles_n_emit_radioactive = {float(particles_n_emit_radioactive):.0e} (radioactive decay photons)")
            print(f"    particles_last_iter_pump = {int(particles_last_iter_pump)} (spectrum averaging window)")
            print(f"    particles_max_total = {float(particles_max_total):.0e} (max photon packets)")
            print(f"    transport_steady_iterate = {transport_steady_iterate} (iterations)")
            # Spectrum grid info
            spec_lam_min = 2.998e10 / spectrum_nu_max * 1e8  # Angstrom
            spec_lam_max = 2.998e10 / spectrum_nu_min * 1e8  # Angstrom
            import math
            n_spec_bins = int(math.log(spectrum_nu_max / spectrum_nu_min) / math.log(1 + spectrum_nu_grid_spacing))
            dlam = (spec_lam_max - spec_lam_min) / n_spec_bins
            dv_kms = spectrum_nu_grid_spacing * 2.998e5  # km/s
            print(f"    spectrum_nu_grid = {{{spectrum_nu_min:.1e}, {spectrum_nu_max:.1e}, {spectrum_nu_grid_spacing}, 1}}")
            print(f"      -> {spec_lam_min:.0f}-{spec_lam_max:.0f} A, {n_spec_bins} bins, ~{dlam:.1f} A/bin, ~{dv_kms:.0f} km/s velocity resolution")
        if transport_use_ddmc > 0:
            ddmc_methods = {1: 'IMD', 2: 'DDMC (Fleck & Canfield 1984)', 3: 'RandomWalk'}
            print(f"    transport_use_ddmc = {transport_use_ddmc} ({ddmc_methods.get(transport_use_ddmc, 'unknown')})")
            print(f"    transport_ddmc_tau_threshold = {transport_ddmc_tau_threshold} (DDMC active in zones with tau > threshold)")

        param_template.write_text(text)
        print(f"  Generated param.lua:")
        print(f"    Time since energy injection: {snapshot_days:.2f} days")
        print(f"    Inner boundary: R={core_radius:.6e} cm, T={T_inner:.0f} K, L={luminosity:.3e} erg/s")
    else:
        print(f"  WARNING: param.lua template not found in sedona_base")
    
    print(f"SEDONA directory prepared: {sedona_dir}")
    return True


def submit_heracles_job(config, scripts_dir, project_dir, debug=False, config_file=None):
    """[DEPRECATED] Submit HERACLES Slurm job.

    HERACLES is deprecated; use the STELLA -> CMFGEN pipeline instead.
    """
    model_name = config.get('model_name', 'unknown')
    print(f"\n{'='*60}")
    print(f"Submitting [DEPRECATED] HERACLES job for model: {model_name}")
    print(f"{'='*60}")
    
    # Extract paths
    output = config.get('output', {})
    model_dir = output.get('model_dir', '../models/${model_name}').replace('${model_name}', model_name)
    model_dir = (project_dir / model_dir).absolute()
    
    heracles = config.get('heracles', {})  # DEPRECATED: HERACLES backend, use STELLA instead
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
    partition, cpus, memory = get_job_resources(
        debug,
        config_partition=heracles.get('partition'),
        config_mem=heracles.get('memory'),
        config_cpus=heracles.get('cpus'),
    )
    time = get_time_limit(debug, heracles.get('time'))
    
    if debug:
        print(f"DEBUG MODE: Using partition={partition}, time={time}, simulation_days={simulation_days}")
    
    # [DEPRECATED] HERACLES root
    heracles_root = os.environ.get('HERACLES_ROOT', '/ptmp/sbrennan/heracles')
    # Force correct path if environment variable points to heracles-2.0
    if heracles_root.endswith('/heracles-2.0'):
        heracles_root = heracles_root.replace('/heracles-2.0', '')
    
    # Load template
    template_file = scripts_dir / 'deprecated' / 'heracles_job_template.slurm'
    with open(template_file, 'r') as f:
        template = f.read()
    
    # Use debug config file if provided, otherwise use original config.
    # Config files are named <model_name>.yaml (not model_<model_name>.yaml).
    if config_file is None:
        config_file = project_dir / 'configs' / f'{model_name}.yaml'
    
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
        print(f"ERROR: [DEPRECATED] HERACLES job submission failed")
        print(result.stdout)
        print(result.stderr)
        return None
    
    # Parse job ID from output (format: "Submitted batch job 12345")
    job_id = result.stdout.strip().split()[-1]
    print(f"[DEPRECATED] HERACLES job submitted: {job_id}")
    
    return job_id


def submit_conversion_job(config, scripts_dir, project_dir, heracles_job_id, debug=False):
    """[DEPRECATED] Submit CMFGEN conversion job with dependency on HERACLES.

    NOTE: This function is currently unused. The HERACLES workflow embeds the conversion
    step at the end of heracles_job_template.slurm, and the STELLA workflow embeds it
    inside stella_job_template.slurm. Kept for potential future use as a standalone
    conversion job. HERACLES is deprecated; use the STELLA -> CMFGEN pipeline instead.
    """
    model_name = config.get('model_name', 'unknown')
    print(f"\n{'='*60}")
    print(f"Submitting CMFGEN conversion job for model: {model_name}")
    print(f"{'='*60}")
    
    # Extract paths
    output = config.get('output', {})
    model_dir = output.get('model_dir', '../models/${model_name}').replace('${model_name}', model_name)
    model_dir = (project_dir / model_dir).absolute()
    
    # Resource settings (use same as [DEPRECATED] HERACLES for conversion)
    heracles = config.get('heracles', {})  # DEPRECATED: HERACLES backend, use STELLA instead
    partition, cpus, memory = get_job_resources(
        debug,
        config_partition=heracles.get('partition'),
        config_mem=heracles.get('memory'),
        config_cpus=heracles.get('cpus'),
    )
    time = get_time_limit(debug, heracles.get('time'))
    
    # [DEPRECATED] HERACLES root
    heracles_root = os.environ.get('HERACLES_ROOT', '/ptmp/sbrennan/heracles')
    # Force correct path if environment variable points to heracles-2.0
    if heracles_root.endswith('/heracles-2.0'):
        heracles_root = heracles_root.replace('/heracles-2.0', '')

    # Load template
    template_file = scripts_dir / 'deprecated' / 'conversion_job_template.slurm'
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


def run_inline_cmfgen_conversion(config, scripts_dir, project_dir, day):
    """Run the STELLA .swd -> CMFGEN SN_HYDRO_DATA conversion inline (not as a SLURM job).

    Used when reuse_stella=True or skip_stella in the non-time-series path,
    where there is no post-STELLA conversion job to create SN_HYDRO_DATA.

    Returns True if SN_HYDRO_DATA was created, False on failure.
    Always regenerates the file even if it already exists.

    ── Diagnostic: scripts and methods invoked ──────────────────────────

    This function calls build_cmfgen_conversion_cmd() which constructs a
    shell command invoking one of two converter scripts:

      converter='v1' (default):
          scripts/convert_stella_swd_to_cmfgen.py
          — Monolithic converter with tau-based AMR, sigma-aware grid
            coarsening, composition-aware hybrid interpolation, and
            ionization/recombination front detection.
          — Base grid methods: 'log_tau', 'tau_amr_recursive',
            'feature_aware', 'log_radius'.

      converter='v2':
          scripts/convert_stella_swd_to_cmfgen_v2.py
          — Modular grid_mapping converter supporting 1000+ grid points
            via conservative shell-overlap remapping.

    The converter reads the STELLA .swd file (via read_swd()), selects the
    requested epoch (closest to --target-days), and writes SN_HYDRO_DATA
    in CMFGEN outer-to-inner depth ordering (R decreasing).

    After conversion, check_sn_hydro_sigma() reports the fraction of zones
    with non-monotonic velocity (sigma < -1) for diagnostic purposes.
    ─────────────────────────────────────────────────────────────────────
    """
    model_name = config.get('model_name', 'unknown')
    output = config.get('output', {})
    model_dir = output.get('model_dir', '../models/${model_name}').replace('${model_name}', model_name)
    model_dir = (project_dir / model_dir).absolute()

    day_label = f"{day:.2f}days" if day is not None else None
    if day_label:
        cmfgen_dir = model_dir / 'cmfgen' / day_label
        sn_hydro = cmfgen_dir / 'SN_HYDRO_DATA'
    else:
        cmfgen_dir = model_dir / 'cmfgen'
        sn_hydro = cmfgen_dir / 'SN_HYDRO_DATA'

    if sn_hydro.exists():
        print(f"  SN_HYDRO_DATA already exists at {sn_hydro} — removing to regenerate with current config")
        sn_hydro.unlink()

    # Build and run the conversion command
    cmd_str = build_cmfgen_conversion_cmd(config, scripts_dir, project_dir, day)
    print(f"  Running inline conversion for day {day}...")
    result = subprocess.run(cmd_str, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        # v2 converter may return non-zero due to validation warnings even when
        # SN_HYDRO_DATA was successfully written. Check if the file exists.
        if sn_hydro.exists():
            print(f"  WARNING: Converter returned non-zero exit code ({result.returncode})")
            print(f"  But SN_HYDRO_DATA exists — treating validation warnings as non-fatal")
            # Print last few lines of output for diagnostics
            if result.stdout:
                for line in result.stdout.split('\n')[-5:]:
                    if line.strip():
                        print(f"    {line}")
        else:
            print(f"  ERROR: Inline conversion failed for day {day}")
            print(result.stdout[-500:] if result.stdout else "")
            print(result.stderr[-500:] if result.stderr else "")
            return False
    if not sn_hydro.exists():
        print(f"  ERROR: SN_HYDRO_DATA not created at {sn_hydro}")
        return False
    print(f"  SN_HYDRO_DATA created: {sn_hydro}")

    # Generate GREY_SCL_FAC_IN from SN_HYDRO_DATA for better initial T estimates
    generate_grey_scl_fac_in(sn_hydro)

    # Post-conversion: check non-monotonic velocity severity and optionally enforce homologous expansion
    cmfgen_conv_cfg = config.get('cmfgen_conversion', {})
    enforce_homologous = cmfgen_conv_cfg.get('enforce_homologous', False)
    sigma_warn_fraction = cmfgen_conv_cfg.get('sigma_warn_fraction', 0.20)
    check_sn_hydro_sigma(sn_hydro, enforce_homologous=enforce_homologous,
                         sigma_warn_fraction=sigma_warn_fraction)

    return True


def build_cmfgen_conversion_cmd(config, scripts_dir, project_dir, day):
    """Return the shell command to convert the STELLA .swd to CMFGEN SN_HYDRO_DATA for a given day.

    ── Diagnostic: converter script selection ───────────────────────────

    The converter script is selected by cmfgen_conversion.converter in YAML:

      'v1' (default) -> scripts/convert_stella_swd_to_cmfgen.py
          Monolithic converter.  Key internal methods (by base_grid_method):
            'log_tau':
                downsample_grid_tau_based()
                — Log-uniform tau spacing with composition-aware hybrid
                  interpolation.  Composition interfaces detected by
                  detect_composition_interfaces(); grid points inserted
                  at interfaces when composition_grid_insertion=True.
                — Post-base-grid AMR: adaptive_grid_refinement() inserts
                  points at edges violating dlogtau, dlogrho, dv, dlogT,
                  dr/r, sigma, and ne thresholds.
                — Optional: coarsen_negative_sigma_regions() removes zones
                  in velocity-reversal regions (DISABLE for non-monotonic
                  models — set coarsen_negative_sigma=false).
                — Optional: enforce_sigma_conformance() inserts points to
                  satisfy sigma criteria on the final grid (DISABLE to
                  preserve STELLA velocity structure — set
                  enable_sigma_conformance=false).
            'tau_amr_recursive':
                tau_amr_optimiser.recursive_amr()
                — Recursive cell-based AMR with HARD constraints:
                  composition interfaces (min 5 forced points each),
                  shocks (min 7 points), velocity extrema (min 5 points).
                — Optical depth reconstructed via hybrid kappa (STELLA
                  NLTE + Kramers + e-scattering floor).
                — Conservative zone-integrated remap onto final grid.
                — Coarsening pass removes redundant points after AMR
                  converges (all cells within target thresholds).
            'feature_aware':
                downsample_grid_feature_aware()
                — Preserves high-gradient regions during downsampling.
            'log_radius':
                Log-uniform radius spacing (legacy, no tau awareness).

      'v2' -> scripts/convert_stella_swd_to_cmfgen_v2.py
          Modular grid_mapping converter.  Supports 1000+ grid points.
          Strips unsupported flags (ionization-front-aware, etc.).

    ── Diagnostic: key YAML config sections ─────────────────────────────

    cmfgen_conversion:
      cmfgen_tau_max      — Inner boundary optical depth cutoff (>=50
                             for DIFFUSION boundary, DIF=T in VADAT).
      base_grid_method    — Selects the grid construction algorithm.
      base_grid_target    — Initial coarse grid point count.
      max_points / min_points — AMR point budget bounds.
      composition_aware   — Detect composition interfaces for hybrid
                             mass-fraction interpolation.
      composition_grid_insertion — Insert grid points AT composition
                             interfaces (true = preserve shell boundaries).
      composition_gradient_threshold — Fractional |dX|/X that triggers
                             interface detection (0.3 = 30%).
      coarsen_negative_sigma — Remove zones in velocity-reversal regions
                             (set FALSE to preserve STELLA velocity).
      enable_sigma_conformance — Insert points to enforce sigma criteria
                             on final grid (set FALSE to avoid modifying
                             the STELLA velocity structure).
      sigma_optimize      — Use exact CMFGEN sigma formula for AMR scoring
                             (set FALSE to avoid sigma-driven modifications).
      use_tau_amr         — Enable tau-based AMR after base grid.
    ─────────────────────────────────────────────────────────────────────
    """
    model_name = config.get('model_name', 'unknown')
    output = config.get('output', {})
    model_dir = output.get('model_dir', '../models/${model_name}').replace('${model_name}', model_name)
    model_dir = (project_dir / model_dir).absolute()

    stella_cfg = config.get('stella', {})
    cmfgen_conv_cfg = config.get('cmfgen_conversion', {})
    cmfgen_cfg = config.get('cmfgen', {})
    kappa = stella_cfg.get('kappa', 0.34)
    # INC_RAD_DECAYS: when False, the converter strips isotope blocks from
    # SN_HYDRO_DATA so CMFGEN does not solve the Ni56->Co56->Fe56 decay chain.
    # This is the default for SN_B interaction models where radioactive heating
    # is not the focus and the implicit Ni/Co species cause convergence issues
    # at low temperatures (T < 4000 K) in the inner ejecta.
    inc_rad_decays = cmfgen_cfg.get('inc_rad_decays', False)
    # Converter selection: 'v2' uses the modular grid_mapping converter
    # (convert_stella_swd_to_cmfgen_v2.py), 'v1' (default) uses the monolithic
    # convert_stella_swd_to_cmfgen.py.  The v2 converter supports up to 1000+
    # grid points via conservative shell-overlap remapping.
    converter_version = cmfgen_conv_cfg.get('converter', 'v1')
    converter_script = 'convert_stella_swd_to_cmfgen_v2.py' if converter_version == 'v2' else 'convert_stella_swd_to_cmfgen.py'
    tau_max = cmfgen_conv_cfg.get('cmfgen_tau_max', stella_cfg.get('stella_tau_max', 100.0))
    tau_min = cmfgen_conv_cfg.get('cmfgen_tau_min', stella_cfg.get('stella_tau_min', 0))
    if isinstance(tau_max, str):
        tau_max = float(tau_max)
    if isinstance(tau_min, str):
        tau_min = float(tau_min)
    if tau_max < 50:
        print(f"  WARNING: cmfgen_tau_max={tau_max} < 50 — inner boundary tau will be < 50.")
        print(f"           DIFFUSION boundary (DIF=T) requires tau>=50 for LTE validity.")
        print(f"           With tau<50, FG_J_CMF_V13 may produce NEGATIVE mean intensities")
        print(f"           and halt immediately (observed on X_v2 with tau_max=30).")
        print(f"           For interaction models with steep tau jumps at shell-ejecta")
        print(f"           interfaces, use tau_max>=100 to push the inner boundary past")
        print(f"           the jump into the smooth inner shell.")
        print(f"           prepare_VADAT.py will auto-detect tau at inner boundary from SN_HYDRO_DATA")
        print(f"           and use ZERO_FLUX (DIF=F) if tau<50, or DIFFUSION (DIF=T) if tau>=50.")
    base_grid_target = cmfgen_conv_cfg.get('base_grid_target', 50)
    base_grid_method = cmfgen_conv_cfg.get('base_grid_method', 'tau_amr_recursive')
    tau_amr_config = cmfgen_conv_cfg.get('tau_amr_config', None)
    tau_breakpoints = cmfgen_conv_cfg.get('tau_breakpoints', [0.01, 1.0])
    tau_zone_fractions = cmfgen_conv_cfg.get('tau_zone_fractions', [0.10, 0.25, 0.65])
    tau_spacing = cmfgen_conv_cfg.get('tau_spacing', 'auto')
    # Composition-aware grid options
    composition_aware = cmfgen_conv_cfg.get('composition_aware', True)
    composition_gradient_threshold = cmfgen_conv_cfg.get('composition_gradient_threshold', 0.5)
    composition_elements = cmfgen_conv_cfg.get('composition_elements', None)
    min_zones_per_shell = cmfgen_conv_cfg.get('min_zones_per_shell', 3)
    composition_trace_threshold = cmfgen_conv_cfg.get('composition_trace_threshold', 1e-4)
    composition_grid_insertion = cmfgen_conv_cfg.get('composition_grid_insertion', True)
    composition_interface_cap = cmfgen_conv_cfg.get('composition_interface_cap', 0.15)
    max_thin_zones = cmfgen_conv_cfg.get('max_thin_zones', 0)
    max_thick_tau = cmfgen_conv_cfg.get('max_thick_tau', 1e10)
    min_points = cmfgen_conv_cfg.get('min_points', 150)
    max_points = cmfgen_conv_cfg.get('max_points', 150)
    min_nd = cmfgen_conv_cfg.get('min_nd', 150)
    # --- target_nd override: enforce a hard grid-point cap in the converter ---
    # When cmfgen.target_nd is set, the converter itself produces a grid with
    # exactly target_nd points, eliminating the need for RDINR downsampling at
    # CMFGEN runtime.  This overrides the cmfgen_conversion grid-size settings
    # so all post-AMR steps (outer boundary insertion, min-nd fill, sigma
    # conformance) are capped at target_nd.
    target_nd = cmfgen_cfg.get('target_nd', None)
    if target_nd is not None:
        target_nd = int(target_nd)
        max_points = target_nd
        min_points = target_nd
        min_nd = target_nd
        base_grid_target = target_nd
        sigma_conformance_headroom = 0
        outer_boundary_extra_points = 0
        print(f"[CMFGEN] target_nd={target_nd}: overriding converter grid size "
              f"(max_points=min_points=min_nd=base_grid_target={target_nd}, "
              f"sigma_conformance_headroom=0, outer_boundary_extra_points=0)")
    else:
        sigma_conformance_headroom = cmfgen_conv_cfg.get('sigma_conformance_headroom', 40)
        outer_boundary_extra_points = cmfgen_conv_cfg.get('outer_boundary_extra_points', 5)
    max_radius_jump = cmfgen_conv_cfg.get('max_radius_jump', 0.05)
    max_temp_jump = cmfgen_conv_cfg.get('max_temp_jump', 0.10)
    max_density_jump = cmfgen_conv_cfg.get('max_density_jump', 0.15)
    max_ne_jump = cmfgen_conv_cfg.get('max_ne_jump', 0.10)
    max_logtau_jump = cmfgen_conv_cfg.get('max_logtau_jump', 0.10)
    max_velocity_jump = cmfgen_conv_cfg.get('max_velocity_jump', 0.10)
    max_sigma_jump = cmfgen_conv_cfg.get('max_sigma_jump', 1.0)
    max_cell_ratio = cmfgen_conv_cfg.get('max_cell_ratio', 1.5)
    enable_adaptive_grid = cmfgen_conv_cfg.get('enable_adaptive_grid', True)
    use_tau_amr = cmfgen_conv_cfg.get('use_tau_amr', False)
    smooth_physical_discontinuities = cmfgen_conv_cfg.get('smooth_physical_discontinuities', False)
    phys_smooth_sigma_zones = cmfgen_conv_cfg.get('phys_smooth_sigma_zones', 10.0)
    stretch_trouble_regions = cmfgen_conv_cfg.get('stretch_trouble_regions', False)
    stretch_factor = cmfgen_conv_cfg.get('stretch_factor', 3.0)
    max_stretch_zones = cmfgen_conv_cfg.get('max_stretch_zones', 200)
    min_zones_to_stretch = cmfgen_conv_cfg.get('min_zones_to_stretch', 2)
    # Critical-region detection and stretching (for sharp discontinuities like
    # reverse shocks where STELLA has very few zones spanning the transition).
    critical_regions = cmfgen_conv_cfg.get('critical_regions', False)
    stretch_critical_regions = cmfgen_conv_cfg.get('stretch_critical_regions', False)
    critical_stretch_factor = cmfgen_conv_cfg.get('critical_stretch_factor', 5.0)
    critical_max_stretch_zones = cmfgen_conv_cfg.get('critical_max_stretch_zones', 300)
    critical_dT_T = cmfgen_conv_cfg.get('critical_dT_T', 0.10)
    critical_drho_rho = cmfgen_conv_cfg.get('critical_drho_rho', 0.10)
    critical_dlogtau = cmfgen_conv_cfg.get('critical_dlogtau', 0.10)
    critical_dv_v = cmfgen_conv_cfg.get('critical_dv_v', 0.30)
    critical_dne_ne = cmfgen_conv_cfg.get('critical_dne_ne', 0.10)
    critical_pad_zones = cmfgen_conv_cfg.get('critical_pad_zones', 2)
    thin_grid_tau = cmfgen_conv_cfg.get('thin_grid_tau', 1.0)
    thin_grid_priority = cmfgen_conv_cfg.get('thin_grid_priority', 1.0)
    if target_nd is None:
        outer_boundary_extra_points = cmfgen_conv_cfg.get('outer_boundary_extra_points', 5)
    outer_boundary_ratio = cmfgen_conv_cfg.get('outer_boundary_ratio', 1.5)
    min_dr_r = cmfgen_conv_cfg.get('min_dr_r', 1e-6)
    max_dr_r = cmfgen_conv_cfg.get('max_dr_r', 0.5)
    sigma_priority_threshold = cmfgen_conv_cfg.get('sigma_priority_threshold', 2.0)
    sigma_priority = cmfgen_conv_cfg.get('sigma_priority', 5.0)
    dsigma_priority_threshold = cmfgen_conv_cfg.get('dsigma_priority_threshold', 0.5)
    dsigma_priority = cmfgen_conv_cfg.get('dsigma_priority', 5.0)
    t_floor = cmfgen_conv_cfg.get('t_floor_K', 3000.0)
    t_cap = cmfgen_conv_cfg.get('t_cap_K', 0.0)
    trace_zero_threshold = cmfgen_conv_cfg.get('trace_zero_threshold', cmfgen_conv_cfg.get('trace_threshold', 0.0))
    velocity_max_kms = cmfgen_conv_cfg.get('velocity_max_kms', 0.0)
    v_min_km_s = cmfgen_conv_cfg.get('v_min_km_s', 0.0)
    auto_boundaries = cmfgen_conv_cfg.get('auto_boundaries', False)
    enable_sigma_conformance = cmfgen_conv_cfg.get('enable_sigma_conformance', True)
    if target_nd is None:
        sigma_conformance_headroom = cmfgen_conv_cfg.get('sigma_conformance_headroom', 40)
    max_abs_sigma = cmfgen_conv_cfg.get('max_abs_sigma', 5.0)
    sigma_optimize = cmfgen_conv_cfg.get('sigma_optimize', True)

    # Ionization-front-aware grid enhancement (v2: ensures fronts resolved)
    ionization_front_aware = cmfgen_conv_cfg.get('ionization_front_aware', True)
    ionization_front_threshold = cmfgen_conv_cfg.get('ionization_front_threshold', 0.5)
    ionization_front_min_points = cmfgen_conv_cfg.get('ionization_front_min_points', 5)
    ionization_front_width = cmfgen_conv_cfg.get('ionization_front_width', 3)
    ionization_front_cap = cmfgen_conv_cfg.get('ionization_front_cap', 0.15)

    # Recombination-front resolution (resolve under-resolved STELLA fronts)
    resolve_recombination_fronts = cmfgen_conv_cfg.get('resolve_recombination_fronts', True)
    recombination_t_jump = cmfgen_conv_cfg.get('recombination_t_jump', 0.3)
    recombination_front_width = cmfgen_conv_cfg.get('recombination_front_width', 5)

    # Ne AMR and sigma smoothing — must match run_stella_conversion() logic
    include_ne_amr = cmfgen_conv_cfg.get('include_ne_amr', False)
    ne_priority_threshold = cmfgen_conv_cfg.get('ne_priority_threshold', 2.0)
    ne_priority = cmfgen_conv_cfg.get('ne_priority', 1.0)
    max_dsigma = cmfgen_conv_cfg.get('max_dsigma', 1.0)
    sigma_smooth_passes = cmfgen_conv_cfg.get('sigma_smooth_passes', 10)
    sigma_smooth_window = cmfgen_conv_cfg.get('sigma_smooth_window', 3)

    # Velocity smoothing flags — must match run_stella_conversion() logic
    smooth_velocity = cmfgen_conv_cfg.get('smooth_velocity', False)
    no_velocity_smoothing = cmfgen_conv_cfg.get('no_velocity_smoothing', True)
    smooth_sigma = cmfgen_conv_cfg.get('smooth_sigma', False)
    smooth_wiggles = cmfgen_conv_cfg.get('smooth_wiggles', False)
    wiggle_window = cmfgen_conv_cfg.get('wiggle_window', 3)
    wiggle_max_change = cmfgen_conv_cfg.get('wiggle_max_change', 0.02)
    velocity_smooth_window_kms = cmfgen_conv_cfg.get('velocity_smooth_window_kms', 0.0)
    smooth_sigma_scatter = cmfgen_conv_cfg.get('smooth_sigma_scatter', False)
    max_dsigma_scatter = cmfgen_conv_cfg.get('max_dsigma_scatter', 0.5)
    sigma_scatter_window = cmfgen_conv_cfg.get('sigma_scatter_window', 7)
    sigma_scatter_dv_threshold = cmfgen_conv_cfg.get('sigma_scatter_dv_threshold', 0.01)
    smooth_sigma_global = cmfgen_conv_cfg.get('smooth_sigma_global', False)
    sigma_global_window = cmfgen_conv_cfg.get('sigma_global_window', 7)
    sigma_global_max_rel_change = cmfgen_conv_cfg.get('sigma_global_max_rel_change', 0.05)
    sigma_global_preserve_threshold = cmfgen_conv_cfg.get('sigma_global_preserve_threshold', 2.0)
    # Final Savitzky-Golay smoothing (default OFF — creates sigma spikes at shock fronts)
    enable_final_smoothing = cmfgen_conv_cfg.get('enable_final_smoothing', False)

    # Sigma-aware grid coarsening (Luc Dessart approach): removes grid points
    # from negative-sigma (velocity reversal) regions.  Default True for
    # backward compatibility, but should be False for non-monotonic models
    # where the velocity reversal structure must be preserved.
    coarsen_negative_sigma = cmfgen_conv_cfg.get('coarsen_negative_sigma', True)
    coarsen_sigma_threshold = cmfgen_conv_cfg.get('coarsen_sigma_threshold', -0.1)
    coarsen_max_zones = cmfgen_conv_cfg.get('coarsen_max_zones', 2)
    coarsen_min_region_zones = cmfgen_conv_cfg.get('coarsen_min_region_zones', 4)

    day_label = f"{day:.2f}days"
    cmfgen_day_dir = model_dir / 'cmfgen' / day_label
    swd_source = cmfgen_conv_cfg.get('stella_source_swd')
    swd_file = Path(swd_source).expanduser() if swd_source else model_dir / 'stella' / 'strad' / f'{model_name}.swd'
    fallback_swd = swd_file if swd_source else model_dir / 'stella' / f'{model_name}.swd'
    cmfgen_out = cmfgen_day_dir / 'SN_HYDRO_DATA'

    # Diagnostic plots go to plots/ subdirectory, not cluttering the run directory
    plots_dir = cmfgen_day_dir / 'plots'

    rho_file = model_dir / 'stella' / 'eve' / f'{model_name}.rho'
    abn_file = model_dir / 'stella' / 'modmake' / f'{model_name}.abn'
    input_profile = model_dir / 'stella' / 'input_profile'
    xni_file = model_dir / 'stella' / 'eve' / f'{model_name}.xni'

    adaptive_arg = " --enable-adaptive-grid" if enable_adaptive_grid else ""
    phys_smooth_arg = ""
    if smooth_physical_discontinuities:
        phys_smooth_arg = f" --smooth-physical-discontinuities --phys-smooth-sigma-zones {phys_smooth_sigma_zones}"
    stretch_arg = ""
    if stretch_trouble_regions:
        stretch_arg = f" --stretch-trouble-regions --stretch-factor {stretch_factor} --max-stretch-zones {max_stretch_zones} --min-zones-to-stretch {min_zones_to_stretch}"

    # Velocity smoothing flags — must match run_stella_conversion() logic.
    # v2 converter defaults: smooth_wiggles=False, smooth_sigma=False,
    # smooth_sigma_scatter=False, smooth_sigma_global=False.
    # Orchestrator defaults match the converter defaults (v2 fidelity-first).
    smooth_args = ""
    if smooth_velocity and not no_velocity_smoothing:
        smooth_args += " --smooth-velocity"
    if smooth_sigma:
        smooth_args += " --smooth-sigma"
    else:
        smooth_args += " --no-smooth-sigma"
    if smooth_sigma:
        smooth_args += f" --max-dsigma {max_dsigma}"
        smooth_args += f" --sigma-smooth-passes {sigma_smooth_passes}"
        smooth_args += f" --sigma-smooth-window {sigma_smooth_window}"
    if not smooth_wiggles:
        smooth_args += " --no-smooth-wiggles"
    else:
        smooth_args += f" --wiggle-window {wiggle_window}"
        smooth_args += f" --wiggle-max-change {wiggle_max_change}"
    if velocity_smooth_window_kms > 0:
        smooth_args += f" --velocity-smooth-window-kms {velocity_smooth_window_kms}"
    if smooth_sigma_scatter:
        smooth_args += " --smooth-sigma-scatter"
        smooth_args += f" --max-dsigma-scatter {max_dsigma_scatter}"
        smooth_args += f" --sigma-scatter-window {sigma_scatter_window}"
        smooth_args += f" --sigma-scatter-dv-threshold {sigma_scatter_dv_threshold}"
    else:
        smooth_args += " --no-smooth-sigma-scatter"
    if smooth_sigma_global:
        smooth_args += " --smooth-sigma-global"
        smooth_args += f" --sigma-global-window {sigma_global_window}"
        smooth_args += f" --sigma-global-max-rel-change {sigma_global_max_rel_change}"
        smooth_args += f" --sigma-global-preserve-threshold {sigma_global_preserve_threshold}"
    else:
        smooth_args += " --no-smooth-sigma-global"
    # Final Savitzky-Golay smoothing (default OFF)
    if enable_final_smoothing:
        smooth_args += " --enable-final-smoothing"
    else:
        smooth_args += " --no-final-smoothing"

    # STELLA models the full evolution for all snapshot days
    homo_prefix = (
        f"SWD='{swd_file}'; "
        f"[[ -f \"$SWD\" ]] || SWD='{fallback_swd}'; "
    )

    # Pre-compute composition-elements argument string (avoids nested f-string quoting issues)
    comp_elements_arg = ''
    if composition_elements:
        comp_elements_arg = f"--composition-elements {','.join(str(x) for x in composition_elements)} "

    if converter_version == 'v2':
        auto_boundary_arg = " --auto-boundaries" if auto_boundaries else ""
        # Non-monotonic mode: preserve negative sigma, enable shock-weighted grid clustering
        non_mono_arg = " --non-monotonic"
        shock_weight = cmfgen_conv_cfg.get('shock_weight', 5.0)
        shock_arg = f" --shock-weight {shock_weight}"
        # No sigma clamp — the formal solver handles extreme sigma via
        # semi-implicit Q-splitting (Baron & Hauschildt 2004).
        return (
            f"{homo_prefix}"
            f"CMFGEN_RHO_ARG=''; CMFGEN_XNI_ARG=''; "
            f"[[ -f '{rho_file}' ]] && CMFGEN_RHO_ARG=\"--rho-file '{rho_file}'\"; "
            f"[[ -f '{xni_file}' ]] && CMFGEN_XNI_ARG=\"--xni-file '{xni_file}'\"; "
            f"mkdir -p '{cmfgen_day_dir}' '{plots_dir}'; "
            f"{sys.executable} {scripts_dir}/{converter_script} \"$SWD\" '{cmfgen_out}' "
            f"--target-days {day} --tau-max {tau_max} --tau-min {tau_min} "
            f"--t-floor {t_floor} --target-points {base_grid_target}"
            f"{auto_boundary_arg}{non_mono_arg}{shock_arg} "
            f"$CMFGEN_RHO_ARG $CMFGEN_XNI_ARG"
        )

    cmd = (
        f"{homo_prefix}"
        f"CMFGEN_RHO_ARG=''; CMFGEN_ABN_ARG=''; CMFGEN_PROFILE_ARG=''; CMFGEN_XNI_ARG=''; "
        f"[[ -f '{rho_file}' ]] && CMFGEN_RHO_ARG=\"--rho-file '{rho_file}'\"; "
        f"[[ -f '{abn_file}' ]] && CMFGEN_ABN_ARG=\"--abn-file '{abn_file}'\"; "
        f"[[ -f '{input_profile}' ]] && CMFGEN_PROFILE_ARG=\"--input-profile '{input_profile}'\"; "
        f"[[ -f '{xni_file}' ]] && CMFGEN_XNI_ARG=\"--xni-file '{xni_file}'\"; "
        f"{sys.executable} {scripts_dir}/{converter_script} \"$SWD\" '{cmfgen_out}' "
        f"--target-days {day} --kappa {kappa} --tau-max {tau_max} --tau-min {tau_min} "
        f"--model-name {model_name} --base-grid-target {base_grid_target} "
        f"--base-grid-method {base_grid_method} "
        f"{f'--tau-amr-config {tau_amr_config} ' if tau_amr_config else ''}"
        f"{'--use-tau-amr ' if use_tau_amr else ''}"
        f"--tau-breakpoints {','.join(str(x) for x in tau_breakpoints)} "
        f"--tau-zone-fractions {','.join(str(x) for x in tau_zone_fractions)} "
        f"--tau-spacing {tau_spacing} "
        f"--composition-gradient-threshold {composition_gradient_threshold} "
        f"--min-zones-per-shell {min_zones_per_shell} "
        f"--trace-threshold {composition_trace_threshold} "
        f"--composition-interface-cap {composition_interface_cap} "
        f"{'--no-composition-grid-insertion' if not composition_grid_insertion else ''} "
        f"{'--no-composition-aware' if not composition_aware else ''} "
        f"{comp_elements_arg}"
        f"--max-thin-zones {max_thin_zones} --max-thick-tau {max_thick_tau} "
        f"--min-points {min_points} "
        f"--max-points {max_points} --min-nd {min_nd} --max-radius-jump {max_radius_jump} --max-temp-jump {max_temp_jump} "
        f"--max-density-jump {max_density_jump} --max-ne-jump {max_ne_jump} --max-logtau-jump {max_logtau_jump} "
        f"--max-velocity-jump {max_velocity_jump} --max-sigma-jump {max_sigma_jump} --max-cell-ratio {max_cell_ratio} "
        f"--thin-grid-tau {thin_grid_tau} --thin-grid-priority {thin_grid_priority} "
        f"--outer-boundary-extra-points {outer_boundary_extra_points} "
        f"--outer-boundary-ratio {outer_boundary_ratio} "
        f"--min-dr-r {min_dr_r} --max-dr-r {max_dr_r} "
        f"--sigma-priority-threshold {sigma_priority_threshold} --sigma-priority {sigma_priority} "
        f"--dsigma-priority-threshold {dsigma_priority_threshold} --dsigma-priority {dsigma_priority} "
        f"--ne-priority-threshold {ne_priority_threshold} --ne-priority {ne_priority} "
        f"--t-floor {t_floor} --t-cap {t_cap} --trace-zero-threshold {trace_zero_threshold} --velocity-max-kms {velocity_max_kms} --v-min-km-s {v_min_km_s} "
        f"--sigma-conformance-headroom {sigma_conformance_headroom} "
        f"--max-abs-sigma {max_abs_sigma} "
        f"{'--enable-sigma-conformance' if enable_sigma_conformance else '--no-sigma-conformance'} "
        f"{'--no-sigma-optimize' if not sigma_optimize else ''} "
        f"{'--ionization-front-aware' if ionization_front_aware else '--no-ionization-front-aware'} "
        f"--ionization-front-threshold {ionization_front_threshold} "
        f"--ionization-front-min-points {ionization_front_min_points} "
        f"--ionization-front-width {ionization_front_width} "
        f"--ionization-front-cap {ionization_front_cap} "
        f"{'--resolve-recombination-fronts' if resolve_recombination_fronts else '--no-resolve-recombination-fronts'} "
        f"--recombination-t-jump {recombination_t_jump} "
        f"--recombination-front-width {recombination_front_width} "
        f"--plot-dir {plots_dir} "
        f"{'--no-inc-rad-decays' if not inc_rad_decays else ''} "
        f"$CMFGEN_RHO_ARG $CMFGEN_ABN_ARG $CMFGEN_PROFILE_ARG $CMFGEN_XNI_ARG{adaptive_arg}{phys_smooth_arg}{stretch_arg}{smooth_args}"
        f"{' --critical-regions' if critical_regions else ''}"
        f"{' --stretch-critical-regions' if stretch_critical_regions else ''}"
        f"{f' --critical-stretch-factor {critical_stretch_factor}' if stretch_critical_regions else ''}"
        f"{f' --critical-max-stretch-zones {critical_max_stretch_zones}' if stretch_critical_regions else ''}"
        f"{f' --critical-dT-T {critical_dT_T}' if critical_regions else ''}"
        f"{f' --critical-drho-rho {critical_drho_rho}' if critical_regions else ''}"
        f"{f' --critical-dlogtau {critical_dlogtau}' if critical_regions else ''}"
        f"{f' --critical-dv-v {critical_dv_v}' if critical_regions else ''}"
        f"{f' --critical-dne-ne {critical_dne_ne}' if critical_regions else ''}"
        f"{f' --critical-pad-zones {critical_pad_zones}' if critical_regions else ''}"
    )
    # Coarsen-negative-sigma flag: v2 converter doesn't support this option
    if converter_version != 'v2':
        cmd += (
            f" {'--no-coarsen-negative-sigma' if not coarsen_negative_sigma else '--coarsen-negative-sigma --coarsen-sigma-threshold ' + str(coarsen_sigma_threshold) + ' --coarsen-max-zones ' + str(coarsen_max_zones) + ' --coarsen-min-region-zones ' + str(coarsen_min_region_zones)}"
        )

    # Handle include_ne_amr: if False, override max_ne_jump to 1e10 to disable Ne AMR
    # (must match run_stella_conversion() logic at line 3196)
    if not include_ne_amr:
        cmd = cmd.replace("--max-ne-jump {}".format(max_ne_jump), "--max-ne-jump 1e10")

    # v2 converter: strip flags it doesn't support
    if converter_version == 'v2':
        import re as _re
        _v2_unsupported = [
            r'--ionization-front-aware\s+',
            r'--no-ionization-front-aware\s+',
            r'--ionization-front-threshold\s+\S+\s*',
            r'--ionization-front-min-points\s+\S+\s*',
            r'--ionization-front-width\s+\S+\s*',
            r'--ionization-front-cap\s+\S+\s*',
            r'--resolve-recombination-fronts\s+',
            r'--no-resolve-recombination-fronts\s+',
            r'--recombination-t-jump\s+\S+\s*',
            r'--recombination-front-width\s+\S+\s*',
            r'--refine-outer-boundary\s+',
            r'--outer-refine-cells\s+\S+\s*',
            r'--outer-refine-max-points\s+\S+\s*',
            r'--outer-refine-subdivisions\s+\S+\s*',
            r'--outer-refine-target-dr-r\s+\S+\s*',
        ]
        for pattern in _v2_unsupported:
            cmd = _re.sub(pattern, '', cmd)
        # Clean up any double spaces left by removals
        cmd = _re.sub(r'  +', ' ', cmd)

    return cmd


def optimize_model_spec_from_hydro(config, scripts_dir, project_dir, day):
    """Run the MODEL_SPEC optimiser against SN_HYDRO_DATA and save an audit report.

    ── Diagnostic: script invoked ───────────────────────────────────────

    Calls scripts/optimize_model_spec_v2.py which:
      1. Parses SN_HYDRO_DATA (HydroReader using proxy_populations_v2.parse_sn_hydro_data)
      2. Computes Saha/LTE ionisation fractions across all depth zones
         (SahaEstimator) as a DIAGNOSTIC candidate generator
      3. Groups ejecta into physical zones by composition/T/ionisation
         (PhysicalZoneDetector)
      4. Selects required ionisation stages per element (IonSelector):
         - ESSENTIAL: max Saha fraction > 1e-6
         - OPTIONAL:  max Saha fraction < 1e-6 (can be omitted)
      5. Selects NLTE/super/full levels per ion (LevelSelector) based on
         LTE population thresholds against available F_OSCDAT data
      6. Generates an optimised MODEL_SPEC (ModelSpecGenerator)
      7. Produces a machine- and human-readable audit report
         (AuditReporter) documenting why each ion/level was selected

    The audit report is saved to <cmfgen_dir>/MODEL_SPEC.audit.md.
    The optimised MODEL_SPEC is NOT automatically written to the run
    directory — the cmfgen_base/MODEL_SPEC template (which has been
    updated with the correct ISF entries) is used by prepare_VADAT.py.
    The audit report serves as the machine-readable justification for
    the ion/level selection.

    This implements the AGENTS.md §13-18 requirement that MODEL_SPEC be
    generated from the physical conditions in SN_HYDRO_DATA rather than
    manually maintained.
    ─────────────────────────────────────────────────────────────────────
    """
    model_name = config.get('model_name', 'unknown')
    output = config.get('output', {})
    model_dir = output.get('model_dir', '../models/${model_name}').replace('${model_name}', model_name)
    model_dir = (project_dir / model_dir).absolute()

    day_label = f"{day:.2f}days"
    cmfgen_dir = model_dir / 'cmfgen' / day_label
    sn_hydro = cmfgen_dir / 'SN_HYDRO_DATA'
    model_spec = cmfgen_dir / 'MODEL_SPEC'

    if not sn_hydro.exists():
        print(f"  MODEL_SPEC optimization: SN_HYDRO_DATA not found at {sn_hydro} — skipping")
        return

    # The optimizer needs a MODEL_SPEC to read current ISF values.
    # If the run-dir copy doesn't exist yet, use the base template.
    if not model_spec.exists():
        base_spec = project_dir / 'cmfgen_base' / 'MODEL_SPEC'
        if base_spec.exists():
            model_spec = base_spec
        else:
            print(f"  MODEL_SPEC optimization: no MODEL_SPEC found — skipping")
            return

    optimizer_script = scripts_dir / 'optimize_model_spec_v2.py'
    if not optimizer_script.exists():
        print(f"  MODEL_SPEC optimization: {optimizer_script.name} not found — skipping")
        return

    print(f"\n  Optimizing MODEL_SPEC from SN_HYDRO_DATA...")
    cmd = f"python3 {optimizer_script} {sn_hydro} {model_spec} --cmfgen-dir {cmfgen_dir} --out-dir {cmfgen_dir}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    # Save the audit report
    audit_path = cmfgen_dir / 'MODEL_SPEC.audit.md'
    if result.stdout:
        audit_path.write_text(result.stdout)
        print(f"  MODEL_SPEC audit report saved: {audit_path}")

    # Print a summary (last 20 lines of stdout for the recommendations)
    if result.stdout:
        lines = result.stdout.strip().split('\n')
        for line in lines[-20:]:
            if line.strip():
                print(f"    {line}")
    else:
        print(f"  WARNING: MODEL_SPEC optimizer produced no output")
        if result.stderr:
            print(f"  stderr: {result.stderr[-300:]}")


def submit_stella_post_conversion_job(config, scripts_dir, project_dir, snapshot_days, dependency_job_id=None, debug=False):
    """Submit a single job that converts the STELLA .swd to all CMFGEN/SEDONA snapshot files.

    This job runs once after STELLA finishes and produces:
      - cmfgen/<day>days/SN_HYDRO_DATA for each day in snapshot_days
      - sedona/<day>days/<model>.mod for each day in snapshot_days

    ── Diagnostic: scripts invoked in the conversion job ────────────────

    For each snapshot day, the job runs:

      1. build_cmfgen_conversion_cmd() — constructs the shell command for
         the STELLA .swd -> SN_HYDRO_DATA conversion.  See the docstring
         of build_cmfgen_conversion_cmd() for details on converter script
         selection (v1 vs v2) and base_grid_method options.

      2. scripts/convert_stella_swd_to_cmfgen.py (or _v2.py)
         — Reads the STELLA .swd, selects the epoch closest to the
           requested day, and constructs the CMFGEN grid using the
           configured base_grid_method (log_tau, tau_amr_recursive, etc.).
         — Writes SN_HYDRO_DATA + diagnostic plots to cmfgen/<day>days/.

      3. scripts/convert_stella_swd_to_superlite.py (if SuperLite enabled)
         — Converts SN_HYDRO_DATA to SuperLite input format.

      4. SEDONA .mod conversion (if SEDONA enabled)
         — Converts SN_HYDRO_DATA to SEDONA .mod format.

    This job is submitted with a SLURM dependency on the STELLA job
    (afterok:STELLA_JOB_ID) so it only runs after STELLA completes.
    When skip_stella/reuse_stella is set, the conversion runs inline
    instead (see run_inline_cmfgen_conversion()).
    ─────────────────────────────────────────────────────────────────────
    """
    model_name = config.get('model_name', 'unknown')
    print(f"\n{'='*60}")
    print(f"Submitting post-STELLA conversion job for model: {model_name}")
    print(f"  Snapshots: {snapshot_days}")
    print(f"{'='*60}")

    output = config.get('output', {})
    model_dir = output.get('model_dir', '../models/${model_name}').replace('${model_name}', model_name)
    model_dir = (project_dir / model_dir).absolute()

    sedona_cfg = config.get('sedona', {})
    nzon_target = sedona_cfg.get('nzon_target', 200)
    use_homologous = sedona_cfg.get('use_homologous_expansion', False)
    cmfgen_conv_cfg = config.get('cmfgen_conversion', {})

    swd_source = cmfgen_conv_cfg.get('stella_source_swd')
    swd_file = Path(swd_source).expanduser() if swd_source else model_dir / 'stella' / 'strad' / f'{model_name}.swd'
    fallback_swd = swd_file if swd_source else model_dir / 'stella' / f'{model_name}.swd'

    partition, _, memory = get_job_resources(debug)
    time = get_time_limit(debug, config_value='01:00:00')

    # Build per-day conversion commands
    cmfgen_cmds = []
    sedona_cmds = []
    for day in snapshot_days:
        day_label = f"{day:.2f}days"
        cmfgen_day_dir = model_dir / 'cmfgen' / day_label
        sedona_day_dir = model_dir / 'sedona' / day_label
        cmfgen_out = cmfgen_day_dir / 'SN_HYDRO_DATA'
        sedona_mod = sedona_day_dir / f'{model_name}.mod'

        cmfgen_cmd = build_cmfgen_conversion_cmd(config, scripts_dir, project_dir, day)
        # SEDONA .mod via v2: convert SN_HYDRO_DATA (already produced by cmfgen_cmd above)
        # to SEDONA format. Uses convert_cmfgen_to_sedona_v2.py which writes correct
        # rmin/texp header (the old convert_stella_swd_to_sedona.py is deprecated).
        sedona_v2 = project_dir.parent / 'sedona' / 'convert_cmfgen_to_sedona_v2.py'
        sedona_cmd = (
            f"mkdir -p '{sedona_day_dir}'; "
            f"if [[ -f '{cmfgen_out}' ]]; then "
            f"{sys.executable} {sedona_v2} '{cmfgen_out}' '{sedona_mod}'"
        )
        if use_homologous:
            sedona_cmd += f" --homologous --time-days {day}"
        sedona_cmd += (
            f"; else echo 'WARNING: SN_HYDRO_DATA not found for day {day} "
            f"(epoch may not exist in .swd)'; fi"
        )

        cmfgen_cmds.append(f"mkdir -p '{cmfgen_day_dir}' && {cmfgen_cmd}")
        sedona_cmds.append(sedona_cmd)

    dep_arg = f"#SBATCH --dependency=afterok:{dependency_job_id}" if dependency_job_id else ""
    # Post-conversion sigma check/enforce for SLURM path
    enforce_homologous = cmfgen_conv_cfg.get('enforce_homologous', False)
    sigma_warn_fraction = cmfgen_conv_cfg.get('sigma_warn_fraction', 0.20)

    sigma_check_block = ""
    if enforce_homologous:
        # Build a Python one-liner to enforce homologous expansion on each SN_HYDRO_DATA
        for day in snapshot_days:
            day_label = f"{day:.2f}days"
            cmfgen_day_dir = model_dir / 'cmfgen' / day_label
            hydro_file = cmfgen_day_dir / 'SN_HYDRO_DATA'
            sigma_check_block += (
                f'  if [[ -f "{hydro_file}" ]]; then '
                f'{sys.executable} -c "'
                f"import sys; sys.path.insert(0, '{scripts_dir}'); "
                f"from run_model_grid import check_sn_hydro_sigma; "
                f"from pathlib import Path; "
                f"check_sn_hydro_sigma(Path('{hydro_file}'), enforce_homologous=True, "
                f"sigma_warn_fraction={sigma_warn_fraction})"
                f'" || echo "WARNING: sigma check failed for day {day}"; fi\n'
            )

    # Generate GREY_SCL_FAC_IN for each snapshot day
    grey_scl_block = ""
    for day in snapshot_days:
        day_label = f"{day:.2f}days"
        cmfgen_day_dir = model_dir / 'cmfgen' / day_label
        hydro_file = cmfgen_day_dir / 'SN_HYDRO_DATA'
        grey_scl_block += (
            f'  if [[ -f "{hydro_file}" ]]; then '
            f'{sys.executable} -c "'
            f"import sys; sys.path.insert(0, '{scripts_dir}'); "
            f"from run_model_grid import generate_grey_scl_fac_in; "
            f"from pathlib import Path; "
            f"generate_grey_scl_fac_in(Path('{hydro_file}'))"
            f'" || echo "WARNING: GREY_SCL_FAC_IN generation failed for day {day}"; fi\n'
        )

    cmfgen_block = "\n".join(
        f'  echo "Converting CMFGEN snapshot {day}..."; '
        f'{cmd} || echo "WARNING: CMFGEN conversion failed for day {day} (epoch may not exist in .swd)"'
        for day, cmd in zip(snapshot_days, cmfgen_cmds)
    )
    sedona_block = "\n".join(
        f'  echo "Converting SEDONA snapshot {day}..."; '
        f'{cmd} || echo "WARNING: SEDONA conversion failed for day {day} (epoch may not exist in .swd)"'
        for day, cmd in zip(snapshot_days, sedona_cmds)
    )

    # Build a list of expected output files for the verification step
    expected_files = []
    for day in snapshot_days:
        day_label = f"{day:.2f}days"
        cmfgen_day_dir = model_dir / 'cmfgen' / day_label
        sedona_day_dir = model_dir / 'sedona' / day_label
        expected_files.append((day, str(cmfgen_day_dir / 'SN_HYDRO_DATA'), str(sedona_day_dir / f'{model_name}.mod')))

    verify_block = "\n".join(
        f'  if [[ -f "{hydro}" ]]; then echo "  OK: CMFGEN day {day} ({hydro})"; '
        f'else echo "  SKIP: CMFGEN day {day} (no SN_HYDRO_DATA)"; fi; '
        f'if [[ -f "{mod}" ]]; then echo "  OK: SEDONA day {day} ({mod})"; '
        f'else echo "  SKIP: SEDONA day {day} (no .mod)"; fi'
        for day, hydro, mod in expected_files
    )

    job_script = f"""#!/bin/bash
#SBATCH --job-name={model_name}_convert_all
{dep_arg}
#SBATCH --partition={partition}
#SBATCH --qos=a0008
#SBATCH --time={time}
#SBATCH --mem={memory}
#SBATCH --cpus-per-task=2
#SBATCH --output={model_dir}/logs/convert_all_{model_name}.out
#SBATCH --error={model_dir}/logs/convert_all_{model_name}.err

# NOTE: No set -e — individual snapshot conversions may fail if STELLA crashed
# before that epoch. We continue and convert whatever epochs are available.
set -uo pipefail

mkdir -p {model_dir}/logs

# Resolve .swd file (strad/ subdirectory or main stella directory)
SWD='{swd_file}'
if [[ ! -f "$SWD" ]]; then
  SWD='{fallback_swd}'
fi

if [[ ! -f "$SWD" ]]; then
  echo "ERROR: STELLA .swd file not found at {swd_file} or {fallback_swd}" >&2
  echo "  STELLA may have crashed before writing any output." >&2
  exit 1
fi

SWD_SIZE=$(stat -c%s "$SWD" 2>/dev/null || stat -f%z "$SWD" 2>/dev/null || echo 0)
echo "Using .swd: $SWD ($SWD_SIZE bytes)"
echo "Converting STELLA snapshots to CMFGEN + SEDONA formats..."
echo "(Individual snapshots that are beyond STELLA's final epoch will be skipped)"

# CMFGEN SN_HYDRO_DATA for all days (non-fatal per day)
{cmfgen_block}

# Post-conversion: enforce homologous expansion if enabled
{sigma_check_block}

# SEDONA .mod for all days (non-fatal per day)
{sedona_block}

# Verify which snapshots were actually produced
echo ""
echo "=== Snapshot verification ==="
{verify_block}

echo ""
echo "Snapshot conversions complete (check above for OK/SKIP status)"

# Exit 0 if at least one CMFGEN or SEDONA snapshot was produced; otherwise exit 1.
# This allows downstream CMFGEN/SEDONA jobs (afterok dependency) to proceed for
# the snapshots that were successfully converted, even if STELLA crashed early.
SUCCESS_COUNT=0
{chr(10).join(f'[ -f "{hydro}" ] && SUCCESS_COUNT=$((SUCCESS_COUNT+1)); [ -f "{mod}" ] && SUCCESS_COUNT=$((SUCCESS_COUNT+1))' for _, hydro, mod in expected_files)}

if [ "$SUCCESS_COUNT" -gt 0 ]; then
  echo "$SUCCESS_COUNT snapshot file(s) produced. Conversion job succeeded."
  exit 0
else
  echo "ERROR: No snapshot files were produced. STELLA may have crashed before any output." >&2
  exit 1
fi
"""

    logs_dir = model_dir / 'logs'
    logs_dir.mkdir(parents=True, exist_ok=True)
    job_file = logs_dir / f'convert_all_{model_name}.slurm'
    with open(job_file, 'w') as f:
        f.write(job_script)

    cmd = ['sbatch', str(job_file)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"ERROR: Post-STELLA conversion job submission failed")
        print(result.stdout)
        print(result.stderr)
        return None

    job_id = result.stdout.strip().split()[-1]
    print(f"Post-STELLA conversion job submitted: {job_id}")
    return job_id


def submit_cmfgen_seed_job(config, scripts_dir, project_dir, day, prev_day, prev_final_job_id, debug=False):
    """Submit a seed job that copies the previous day's CMFGEN restart files and runs out2in.sh."""
    model_name = config.get('model_name', 'unknown')
    day_label = f"{day:.2f}days"
    prev_label = f"{prev_day:.2f}days"
    print(f"\n{'='*60}")
    print(f"Submitting CMFGEN seed job for model: {model_name} day {day} (from {prev_label})")
    print(f"{'='*60}")

    output = config.get('output', {})
    model_dir = output.get('model_dir', '../models/${model_name}').replace('${model_name}', model_name)
    model_dir = (project_dir / model_dir).absolute()

    curr_dir = model_dir / 'cmfgen' / day_label
    prev_dir = model_dir / 'cmfgen' / prev_label

    partition, _, memory = get_job_resources(debug)
    time = '00:30:00' if debug else '01:00:00'

    dep_arg = f"#SBATCH --dependency=afterok:{prev_final_job_id}" if prev_final_job_id else ""

    job_script = f"""#!/bin/bash
#SBATCH --job-name={model_name}_seed_{day_label}
{dep_arg}
#SBATCH --partition={partition}
#SBATCH --qos=a0008
#SBATCH --time={time}
#SBATCH --mem={memory}
#SBATCH --cpus-per-task=1
#SBATCH --output={model_dir}/logs/seed_{model_name}_{day_label}.out
#SBATCH --error={model_dir}/logs/seed_{model_name}_{day_label}.err

set -euo pipefail

mkdir -p {curr_dir}

echo "Seeding CMFGEN directory for day {day} from {prev_label}..."

# Copy restart files from previous day's converged CMFGEN run
copied_out=0
for src in {prev_dir}/*OUT; do
  if [[ -f "$src" ]]; then
    cp -p "$src" {curr_dir}/
    copied_out=$((copied_out + 1))
  fi
done
for src in {prev_dir}/RVTJ {prev_dir}/GAMMAS {prev_dir}/EDDFACTOR {prev_dir}/EDDFACTOR_INFO; do
  if [[ -f "$src" ]]; then
    cp -p "$src" {curr_dir}/
  fi
done
if [[ $copied_out -eq 0 || ! -f {prev_dir}/RVTJ ]]; then
  echo "ERROR: previous epoch is not a converged CMFGEN seed: require *OUT and RVTJ in {prev_dir}" >&2
  exit 1
fi

# Generate *_IN files from the copied restart state
CMFDIST="${{CMFDIST:-/viper/ptmp1/sbrennan/SN2024pba/cmfgen_nonmono_dev}}"
if [[ -f "$CMFDIST/com/out2in.sh" ]]; then
  cd {curr_dir}
  "$CMFDIST/com/out2in.sh"
  if [[ ! -f HI_IN ]]; then
    echo "ERROR: out2in.sh did not create species input files" >&2
    exit 1
  fi
  echo "out2in.sh completed successfully"
else
  echo "ERROR: out2in.sh not found at $CMFDIST/com/out2in.sh" >&2
  exit 1
fi

echo "Seed job complete for day {day}"
"""

    logs_dir = model_dir / 'logs'
    logs_dir.mkdir(parents=True, exist_ok=True)
    job_file = logs_dir / f'seed_{model_name}_{day_label}.slurm'
    with open(job_file, 'w') as f:
        f.write(job_script)

    cmd = ['sbatch', str(job_file)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"ERROR: CMFGEN seed job submission failed for day {day}")
        print(result.stdout)
        print(result.stderr)
        return None

    job_id = result.stdout.strip().split()[-1]
    print(f"CMFGEN seed job submitted: {job_id}")
    return job_id


def _generate_eddfactor_from_stella_edd(config, scripts_dir, project_dir, cmfgen_dir, day=None):
    """Convert STELLA .edd file to CMFGEN EDDFACTOR.

    This uses the STELLA radiation field (J_nu, H_nu, Eddington factors)
    directly as the initial guess for CMFGEN, bypassing the 1-iteration
    CMFGEN run entirely.

    The STELLA .edd file is located at:
        model_dir/stella/strad/{model_name}.prf.edd

    The conversion interpolates J_nu from the STELLA depth grid onto the
    CMFGEN depth grid (defined by SN_HYDRO_DATA), handles frequency grid
    remapping, depth ordering reversal, and unit conversion.

    Args:
        config: YAML config dict
        scripts_dir: Path to scripts directory
        project_dir: Path to project directory
        cmfgen_dir: Path to CMFGEN run directory (must contain SN_HYDRO_DATA)
        day: Optional snapshot day label

    Returns:
        True if EDDFACTOR was generated successfully, False otherwise.
    """
    model_name = config.get('model_name', 'unknown')
    output_cfg = config.get('output', {})
    model_dir_str = output_cfg.get('model_dir', '../models/${model_name}').replace('${model_name}', model_name)
    model_dir = (project_dir / model_dir_str).absolute()

    # Locate STELLA .edd file
    stella_output = output_cfg.get('stella_output', 'stella')
    if str(stella_output).startswith('/'):
        stella_out_abs = Path(stella_output)
    else:
        stella_out_abs = model_dir / stella_output

    edd_file = stella_out_abs / 'strad' / f'{model_name}.prf.edd'
    if not edd_file.exists():
        # Try without strad/ subdirectory
        edd_file = stella_out_abs / f'{model_name}.prf.edd'
    if not edd_file.exists():
        print(f"ERROR: STELLA .edd file not found at {edd_file}", file=sys.stderr)
        print(f"  Looked in: {stella_out_abs / 'strad' / f'{model_name}.prf.edd'}", file=sys.stderr)
        print(f"  And: {stella_out_abs / f'{model_name}.prf.edd'}", file=sys.stderr)
        return False

    print(f"  STELLA .edd file: {edd_file}")
    print(f"  .edd file size: {edd_file.stat().st_size} bytes")

    # SN_HYDRO_DATA must exist in the CMFGEN directory
    sn_hydro = cmfgen_dir / 'SN_HYDRO_DATA'
    if not sn_hydro.exists():
        print(f"ERROR: SN_HYDRO_DATA not found in {cmfgen_dir}", file=sys.stderr)
        return False

    # Convert .edd → EDDFACTOR
    edd_converter = project_dir / 'X_stella_radfield' / 'convert_stella_edd_to_cmfgen.py'
    if not edd_converter.exists():
        # Fallback: check if it's in scripts/
        edd_converter = scripts_dir / 'convert_stella_edd_to_cmfgen.py'
    if not edd_converter.exists():
        print(f"ERROR: convert_stella_edd_to_cmfgen.py not found", file=sys.stderr)
        return False

    edd_output = cmfgen_dir / 'EDDFACTOR'
    edd_info = cmfgen_dir / 'EDDFACTOR_INFO'

    # Remove existing EDDFACTOR
    if edd_output.exists():
        edd_output.unlink()
    if edd_info.exists():
        edd_info.unlink()

    convert_cmd = [
        sys.executable, str(edd_converter),
        str(edd_file),
        str(edd_output),
        '--sn-hydro-data', str(sn_hydro),
    ]
    print(f"  Converting STELLA .edd → CMFGEN EDDFACTOR...")
    convert_result = subprocess.run(convert_cmd, capture_output=True, text=True,
                                    cwd=str(project_dir))
    if convert_result.stdout:
        for line in convert_result.stdout.splitlines():
            print(f"    {line}")
    if convert_result.returncode != 0:
        print(f"ERROR: .edd → EDDFACTOR conversion failed", file=sys.stderr)
        if convert_result.stderr:
            print(convert_result.stderr[-500:], file=sys.stderr)
        return False

    # Verify EDDFACTOR was created
    if not edd_output.exists():
        print(f"ERROR: EDDFACTOR not created at {edd_output}", file=sys.stderr)
        return False

    print(f"  EDDFACTOR created: {edd_output} ({edd_output.stat().st_size} bytes)")

    # Verify EDDFACTOR structure
    verify_script = project_dir / 'X_stella_radfield' / 'verify_eddfactor.py'
    if verify_script.exists():
        # Read ND from SN_HYDRO_DATA
        nd_val = None
        try:
            with open(sn_hydro, 'r') as f:
                for line in f:
                    if 'Number of data points' in line:
                        nd_val = int(line.split()[-1])
                        break
        except Exception:
            pass
        if nd_val:
            verify_cmd = [sys.executable, str(verify_script), str(edd_output), str(nd_val)]
            verify_result = subprocess.run(verify_cmd, capture_output=True, text=True)
            if verify_result.stdout:
                for line in verify_result.stdout.splitlines():
                    print(f"    {line}")
            if verify_result.returncode != 0:
                print(f"  WARNING: EDDFACTOR verification reported issues", file=sys.stderr)

    # Generate CGS diagnostic plots
    plot_cgs_eddfactor(cmfgen_dir, sn_hydro, edd_output, model_name, day)

    return True


def plot_cgs_eddfactor(cmfgen_dir, sn_hydro_path, eddfactor_path, model_name, day=None):
    """Generate CGS diagnostic plots for EDDFACTOR and SN_HYDRO_DATA.

    Plots:
      1. R (cm), V (cm/s), T (K), rho (g/cm³) vs depth index
      2. J_nu (erg/cm²/s/Hz/sr) at selected frequencies vs R (cm)
      3. J_nu contour plot (frequency vs radius) in CGS
      4. Eddington factor f=K/J vs R
    """
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import numpy as np
        import struct
    except ImportError:
        print(f"  WARNING: matplotlib/numpy not available, skipping CGS plots")
        return

    plots_dir = cmfgen_dir / 'plots'
    plots_dir.mkdir(exist_ok=True)

    # ── Read SN_HYDRO_DATA (R, V, T, rho in CMFGEN units) ──
    nd = None
    n_mf = 0
    n_iso = 0
    with open(sn_hydro_path, 'r') as f:
        lines = f.readlines()
    for line in lines:
        s = line.strip()
        if 'Number of data points' in s:
            nd = int(s.split()[-1])
        elif 'Number of mass fractions' in s:
            n_mf = int(s.split()[-1])
        elif 'Number of isotopes' in s:
            n_iso = int(s.split()[-1])

    if nd is None:
        print(f"  WARNING: Could not read ND from SN_HYDRO_DATA")
        return

    vals = []
    in_data = False
    for line in lines:
        s = line.strip()
        if not s:
            in_data = True
            continue
        if in_data:
            for p in s.split():
                try:
                    vals.append(float(p))
                except ValueError:
                    pass

    # SN_HYDRO_DATA blocks: R(ND), V(ND), RHO(ND), T(ND), E(ND), mass_fracs, isotopes
    # R is in 1e10 cm, V in km/s, RHO in g/cm³, T in K
    r_1e10 = np.array(vals[:nd])
    v_kms = np.array(vals[nd:2*nd])
    rho = np.array(vals[2*nd:3*nd])
    T = np.array(vals[3*nd:4*nd])

    r_cm = r_1e10 * 1e10  # Convert to cm
    v_cms = v_kms * 1e5   # Convert to cm/s

    # ── Read EDDFACTOR ──
    # Read EDDFACTOR_INFO
    edd_info_path = cmfgen_dir / 'EDDFACTOR_INFO'
    edd_nd = None
    edd_rec_size = None
    nfrus = None
    if edd_info_path.exists():
        with open(edd_info_path, 'r') as f:
            edd_info_lines = f.readlines()
        # Parse Fortran-format INFO file:
        # Line 1: ' 12-Apr-2017 ... !INFO format date'
        # Line 2: ' 23-Jan-2017 ... !File format date'
        # Line 3: ' ND RECL WORD_SIZE UNIT_SIZE INT_SIZE LIT_END'
        # Line 4: column labels
        for line in edd_info_lines:
            stripped = line.strip()
            if not stripped or stripped.startswith('!'):
                continue
            # Skip lines that contain 'format date' or 'LIT_END' or 'ND' header
            if 'format date' in stripped or 'LIT_END' in stripped or stripped.startswith('ND'):
                continue
            # This should be the data line with ND, RECL, ...
            parts = stripped.split()
            if len(parts) >= 2:
                try:
                    edd_nd = int(parts[0])
                    edd_rec_size = int(parts[1]) * 4  # RECL is in 4-byte units
                except ValueError:
                    pass
                break

    if edd_nd is None:
        edd_nd = nd
        edd_rec_size = 8 * (nd + 1)

    # Read R, V from EDDFACTOR direct-access file
    # Direct-access: each record is edd_rec_size bytes, no record markers.
    # Record N is at offset (N-1) * edd_rec_size (1-indexed).
    with open(eddfactor_path, 'rb') as f:
        # Record 10 (INITIAL_RV_REC): R(1:ND) as float64
        f.seek(9 * edd_rec_size)  # record 10 at offset (10-1)*rec_size
        r_bytes = f.read(edd_rec_size)
        if len(r_bytes) >= 8 * edd_nd:
            edd_r = np.frombuffer(r_bytes, dtype=np.float64)[:edd_nd]
        else:
            edd_r = r_cm

        # Record 11: V(1:ND) as float64
        f.seek(10 * edd_rec_size)  # record 11 at offset (11-1)*rec_size
        v_bytes = f.read(edd_rec_size)
        if len(v_bytes) >= 8 * edd_nd:
            edd_v = np.frombuffer(v_bytes, dtype=np.float64)[:edd_nd]
        else:
            edd_v = v_cms

        # Read FINISH_REC (record 5) to get completion marker
        f.seek(4 * edd_rec_size)  # record 5 at offset (5-1)*rec_size
        finish_bytes = f.read(8)
        finish_val = np.frombuffer(finish_bytes, dtype=np.float64)[0] if len(finish_bytes) >= 8 else 0.0

        # Read EDD_CONT_REC (record 3) to get ACCESS_F (points to Eddington factor data)
        f.seek(2 * edd_rec_size)  # record 3 at offset (3-1)*rec_size
        edd_cont_bytes = f.read(4)
        edd_cont_access = np.frombuffer(edd_cont_bytes, dtype=np.int32)[0] if len(edd_cont_bytes) >= 4 else 0

        # Read RV_REC (record 2) to get IREC (points to R,V data)
        f.seek(1 * edd_rec_size)  # record 2 at offset (2-1)*rec_size
        rv_rec_bytes = f.read(4)
        rv_irec = np.frombuffer(rv_rec_bytes, dtype=np.int32)[0] if len(rv_rec_bytes) >= 4 else 10

        # Determine NFRUS: count records from INITIAL_ACCESS_REC (15) to end of file
        file_size = f.seek(0, 2)  # seek to end
        nfrus = (file_size // edd_rec_size) - 15 + 1
        if nfrus < 1:
            nfrus = 50  # fallback

        # Read J_nu from records 15+ (INITIAL_ACCESS_REC)
        # Each record: (RJ(1:ND), FREQ) — ND float64 + 1 float64
        jnu_data = []
        freqs = []
        for i in range(nfrus):
            rec_num = 15 + i
            f.seek((rec_num - 1) * edd_rec_size)
            rec_bytes = f.read(edd_rec_size)
            if len(rec_bytes) < 8 * (edd_nd + 1):
                continue
            rec = np.frombuffer(rec_bytes, dtype=np.float64)
            freqs.append(rec[edd_nd])  # frequency is the last value
            jnu_data.append(rec[:edd_nd])

        freqs = np.array(freqs)
        jnu = np.array(jnu_data) if jnu_data else None

        # Read Eddington factors from records pointed to by EDD_CONT_REC
        # ACCESS_F from record 3 points to where Eddington factor data starts
        edd_factors = []
        if edd_cont_access > 0 and edd_cont_access < 15:
            for i in range(nfrus):
                rec_num = edd_cont_access + i
                f.seek((rec_num - 1) * edd_rec_size)
                rec_bytes = f.read(edd_rec_size)
                if len(rec_bytes) < 8 * (edd_nd + 1):
                    continue
                rec = np.frombuffer(rec_bytes, dtype=np.float64)
                edd_factors.append(rec[:edd_nd])
        edd_factors = np.array(edd_factors) if edd_factors else None

    # ── Plot 1: Structure (R, V, T, rho) in CGS ──
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f'{model_name} CMFGEN Grid (day={day}) — CGS units', fontsize=14)

    ax = axes[0, 0]
    ax.plot(range(nd), r_cm, 'b-', linewidth=0.8)
    ax.set_xlabel('Depth index (outer→inner)')
    ax.set_ylabel('R (cm)')
    ax.set_title('Radius')
    ax.set_yscale('log')

    ax = axes[0, 1]
    ax.plot(range(nd), v_cms, 'r-', linewidth=0.8)
    ax.set_xlabel('Depth index (outer→inner)')
    ax.set_ylabel('V (cm/s)')
    ax.set_title('Velocity')

    ax = axes[1, 0]
    ax.plot(range(nd), T, 'g-', linewidth=0.8)
    ax.set_xlabel('Depth index (outer→inner)')
    ax.set_ylabel('T (K)')
    ax.set_title('Temperature')
    ax.set_yscale('log')

    ax = axes[1, 1]
    ax.plot(range(nd), rho, 'm-', linewidth=0.8)
    ax.set_xlabel('Depth index (outer→inner)')
    ax.set_ylabel('ρ (g/cm³)')
    ax.set_title('Density')
    ax.set_yscale('log')

    plt.tight_layout()
    plot1 = plots_dir / 'structure_cgs.png'
    fig.savefig(plot1, dpi=150)
    plt.close(fig)
    print(f"  CGS plot: {plot1}")

    # ── Plot 2: J_nu at selected frequencies vs R ──
    if jnu is not None and jnu.shape[0] > 0:
        fig, ax = plt.subplots(figsize=(10, 7))
        # Select ~8 frequencies spanning the range
        n_select = min(8, jnu.shape[0])
        indices = np.linspace(0, jnu.shape[0] - 1, n_select, dtype=int)
        colors = plt.cm.viridis(np.linspace(0, 1, n_select))
        for idx, color in zip(indices, colors):
            freq_Hz = freqs[idx]
            wavelength_A = 3e18 / freq_Hz  # Hz to Å
            j_vals = jnu[idx, :]
            # Filter positive values for log plot
            mask = j_vals > 0
            if np.any(mask):
                ax.plot(r_cm[mask], j_vals[mask], '-', color=color, linewidth=0.8,
                        label=f'ν={freq_Hz:.2e} Hz (λ={wavelength_A:.0f} Å)')
        ax.set_xlabel('R (cm)')
        ax.set_ylabel('J_ν (erg cm⁻² s⁻¹ Hz⁻¹ sr⁻¹)')
        ax.set_title(f'{model_name} J_ν from STELLA .edd (day={day})')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.legend(fontsize=8)
        plt.tight_layout()
        plot2 = plots_dir / 'jnu_vs_r_cgs.png'
        fig.savefig(plot2, dpi=150)
        plt.close(fig)
        print(f"  CGS plot: {plot2}")

        # ── Plot 3: J_nu contour (frequency vs radius) ──
        fig, ax = plt.subplots(figsize=(10, 7))
        # Use log10 of J_nu for the contour
        log_jnu = np.log10(np.maximum(jnu, 1e-300))
        # Frequency axis: convert to wavelength in Å for readability
        wavelengths_A = 3e18 / freqs
        im = ax.pcolormesh(r_cm, wavelengths_A, log_jnu, shading='auto', cmap='inferno')
        ax.set_xlabel('R (cm)')
        ax.set_ylabel('Wavelength (Å)')
        ax.set_title(f'{model_name} log₁₀(J_ν) from STELLA .edd (day={day})')
        ax.set_xscale('log')
        ax.set_yscale('log')
        plt.colorbar(im, ax=ax, label='log₁₀ J_ν [erg cm⁻² s⁻¹ Hz⁻¹ sr⁻¹]')
        plt.tight_layout()
        plot3 = plots_dir / 'jnu_contour_cgs.png'
        fig.savefig(plot3, dpi=150)
        plt.close(fig)
        print(f"  CGS plot: {plot3}")

    # ── Plot 4: Eddington factor vs R ──
    if edd_factors is not None and edd_factors.shape[0] > 0:
        fig, ax = plt.subplots(figsize=(10, 7))
        n_select = min(8, edd_factors.shape[0])
        indices = np.linspace(0, edd_factors.shape[0] - 1, n_select, dtype=int)
        colors = plt.cm.viridis(np.linspace(0, 1, n_select))
        for idx, color in zip(indices, colors):
            freq_Hz = freqs[idx]
            wavelength_A = 3e18 / freq_Hz
            f_vals = edd_factors[idx, :]
            mask = np.isfinite(f_vals) & (f_vals > 0)
            if np.any(mask):
                ax.plot(r_cm[mask], f_vals[mask], '-', color=color, linewidth=0.8,
                        label=f'λ={wavelength_A:.0f} Å')
        ax.axhline(y=1/3, color='gray', linestyle='--', alpha=0.5, label='f=1/3 (isotropic)')
        ax.set_xlabel('R (cm)')
        ax.set_ylabel('f = K_ν / J_ν')
        ax.set_title(f'{model_name} Eddington factor from STELLA .edd (day={day})')
        ax.set_xscale('log')
        ax.legend(fontsize=8)
        plt.tight_layout()
        plot4 = plots_dir / 'eddington_factor_cgs.png'
        fig.savefig(plot4, dpi=150)
        plt.close(fig)
        print(f"  CGS plot: {plot4}")


def generate_eddfactor_for_cmfgen(config, scripts_dir, project_dir, cmfgen_dir, config_file=None, day=None, dependency_job_id=None):
    """Generate a physics-motivated EDDFACTOR for CMFGEN FIXED_J initialization.

    This implements the fixed-J workflow (CMFGEN manual p.36-38):
      1. Run CMFGEN for 1 iteration with USE_FIXED_J=F, COMP_F=T
         -> generates EDDFACTOR with correct frequency grid + J from
            STELLA-derived hydro (T, R, V, density)
      2. Optionally replace J with W*B(T,nu) (j_source="grey")
         or keep CMFGEN-computed J (j_source="cmfgen_1iter")
      3. Generate J comparison plot

    The EDDFACTOR is then used by the fixed-J run (USE_FIXED_J=T, COMP_F=F)
    to break the J-feedback loop.

    Args:
        config: YAML config dict
        scripts_dir: Path to scripts directory
        project_dir: Path to project directory
        cmfgen_dir: Path to CMFGEN run directory
        day: Optional snapshot day label

    Returns:
        True if EDDFACTOR was generated successfully, False otherwise.
    """
    cmfgen_cfg = config.get('cmfgen', {})
    fjw = cmfgen_cfg.get('fixed_j_workflow', {})
    if not fjw.get('enabled', False):
        return True  # not enabled — nothing to do

    j_source = fjw.get('j_source', 'cmfgen_1iter')
    timeout_min = fjw.get('init_timeout_minutes', 15)
    # The SLURM job time limit for the 1-iteration run (short — the run itself
    # takes seconds to minutes).  This is separate from the orchestrator wait
    # timeout, which must be much longer to allow for queue pending time.
    edd_slurm_time_min = fjw.get('init_slurm_time_minutes', 15)
    plot_comparison = fjw.get('plot_j_comparison', True)

    print(f"\n{'='*60}")
    print(f"Fixed-J workflow: EDDFACTOR generation")
    print(f"  CMFGEN dir: {cmfgen_dir}")
    print(f"  J source: {j_source}")
    if j_source != 'stella_edd':
        print(f"  Init timeout: {timeout_min} min (orchestrator wait)")
        print(f"  SLURM time: {edd_slurm_time_min} min (job runtime limit)")
    print(f"{'='*60}")

    # ── STELLA .edd source: convert STELLA radiation field to EDDFACTOR ──
    # This bypasses the 1-iteration CMFGEN run entirely. Instead, we use
    # the STELLA .edd file (containing J_nu, H_nu, Eddington factors) and
    # interpolate J_nu onto the CMFGEN grid defined by SN_HYDRO_DATA.
    if j_source == 'stella_edd':
        return _generate_eddfactor_from_stella_edd(
            config, scripts_dir, project_dir, cmfgen_dir, day
        )

    gen_script = scripts_dir / 'generate_eddfactor.py'
    if not gen_script.exists():
        print(f"ERROR: {gen_script} not found", file=sys.stderr)
        return False

    # Prepare proxy populations and VADAT for the 1-iteration CMFGEN run.
    # The SLURM script normally does this on the compute node, but
    # generate_eddfactor.py runs CMFGEN locally, so we need to prepare
    # everything here first.
    print(f"  Preparing proxy populations and VADAT for EDDFACTOR generation...")

    # Create atomic data symlinks (normally done by batch.sh)
    atomic_dir = os.environ.get('ATOMIC', '/u/sbrennan/atomic')
    atomic_links = {
        'HI_F_OSCDAT': f'{atomic_dir}/HYD/I/5dec96/hi_osc.dat',
        'HeI_F_OSCDAT': f'{atomic_dir}/HE/I/5dec96/heioscdat_a7.dat',
        'OI_F_OSCDAT': f'{atomic_dir}/OXY/I/19apr23/osc_data',
        'SkI_F_OSCDAT': f'{atomic_dir}/SIL/I/19apr23/osc_data',
        'Sk2_F_OSCDAT': f'{atomic_dir}/SIL/II/19apr23/osc_data',
        'CaI_F_OSCDAT': f'{atomic_dir}/CA/I/19apr23/osc_data',
        'FeIII_F_OSCDAT': f'{atomic_dir}/FE/III/19apr23/osc_data',
        'FeI_F_OSCDAT': f'{atomic_dir}/FE/I/19apr23/osc_data',
        'Fe2_F_OSCDAT': f'{atomic_dir}/FE/II/19apr23/osc_data',
        'Ca2_F_OSCDAT': f'{atomic_dir}/CA/II/19apr23/osc_data',
    }
    for link_name, target in atomic_links.items():
        link_path = cmfgen_dir / link_name
        if link_path.is_symlink() or link_path.exists():
            link_path.unlink()
        os.symlink(target, link_path)

    # Run prepare_VADAT.py (same args as SLURM script)
    vadat_script = (scripts_dir / 'prepare_VADAT.py').resolve()
    config_file_path = Path(config_file) if config_file else None
    vel_law = cmfgen_cfg.get('vel_law', 12)
    model_type = cmfgen_cfg.get('model_type', 'static_snapshot')
    inc_rad_decays = cmfgen_cfg.get('inc_rad_decays', False)
    no_rel = not cmfgen_cfg.get('incl_rel', False)
    use_hydro_grid = cmfgen_cfg.get('use_hydro_grid', False)
    smooth_velocity_jumps = cmfgen_cfg.get('smooth_velocity_jumps', False)
    smooth_max_dv_ratio = cmfgen_cfg.get('smooth_max_dv_ratio', 2.0)
    smooth_n_insert = cmfgen_cfg.get('smooth_n_insert', 8)
    smooth_rho_sigma = cmfgen_cfg.get('smooth_rho_sigma', 0.0)
    smooth_temp_sigma = cmfgen_cfg.get('smooth_temp_sigma', 0.0)
    blend_composition = cmfgen_cfg.get('blend_composition', False)
    blend_sigma = cmfgen_cfg.get('blend_sigma', 5.0)
    blend_solar_weight = cmfgen_cfg.get('blend_solar_weight', 0.0)
    min_mass_fraction = cmfgen_cfg.get('min_mass_fraction', 0.0)
    monotonize_velocity = cmfgen_cfg.get('monotonize_velocity', False)
    chemical_n_insert = cmfgen_cfg.get('chemical_n_insert', 0)
    chemical_jump_threshold = cmfgen_cfg.get('chemical_jump_threshold', 0.1)
    reversal_n_insert = cmfgen_cfg.get('reversal_n_insert', 0)
    vadat_cmd = [
        sys.executable, str(vadat_script),
        '--single-final', '--convergence-mode', 'stable',
        '--vel-law', str(vel_law),
        '--model-type', model_type,
        '--inc-rad-decays', '1' if inc_rad_decays else '0',
    ]
    if no_rel:
        vadat_cmd.append('--no-relativistic')
    if use_hydro_grid:
        vadat_cmd.append('--use-hydro-grid')
    if smooth_velocity_jumps:
        vadat_cmd.append('--smooth-velocity-jumps')
        vadat_cmd.extend(['--smooth-max-dv-ratio', str(smooth_max_dv_ratio)])
        vadat_cmd.extend(['--smooth-n-insert', str(smooth_n_insert)])
        if smooth_rho_sigma > 0:
            vadat_cmd.extend(['--smooth-rho-sigma', str(smooth_rho_sigma)])
        if smooth_temp_sigma > 0:
            vadat_cmd.extend(['--smooth-temp-sigma', str(smooth_temp_sigma)])
    if blend_composition:
        vadat_cmd.append('--blend-composition')
        vadat_cmd.extend(['--blend-sigma', str(blend_sigma)])
        vadat_cmd.extend(['--blend-solar-weight', str(blend_solar_weight)])
        if min_mass_fraction > 0:
            vadat_cmd.extend(['--min-mass-fraction', str(min_mass_fraction)])
    if monotonize_velocity:
        vadat_cmd.append('--monotonize-velocity')
    if chemical_n_insert > 0:
        vadat_cmd.extend(['--chemical-n-insert', str(chemical_n_insert)])
        vadat_cmd.extend(['--chemical-jump-threshold', str(chemical_jump_threshold)])
    if reversal_n_insert > 0:
        vadat_cmd.extend(['--reversal-n-insert', str(reversal_n_insert)])
    if config_file_path and config_file_path.exists():
        vadat_cmd.extend(['--vadat-overrides', str(config_file_path.resolve())])
    vadat_result = subprocess.run(vadat_cmd, capture_output=True, text=True,
                                  cwd=str(cmfgen_dir))
    if vadat_result.returncode != 0:
        print(f"ERROR: prepare_VADAT.py failed", file=sys.stderr)
        if vadat_result.stderr:
            print(vadat_result.stderr[-500:], file=sys.stderr)
        return False

    # Run prepare_IN_INTS.py
    inits_script = (scripts_dir / 'prepare_IN_INTS.py').resolve()
    init_num_its = cmfgen_cfg.get('init_num_its', 100000)
    inits_cmd = [sys.executable, str(inits_script), '--single-final',
                 '--num-its', str(init_num_its)]
    inits_result = subprocess.run(inits_cmd, capture_output=True, text=True,
                                  cwd=str(cmfgen_dir))
    if inits_result.returncode != 0:
        print(f"ERROR: prepare_IN_INTS.py failed", file=sys.stderr)
        return False

    # Run proxy_populations_v2.py
    proxy_v2 = cmfgen_cfg.get('proxy_populations_v2', True)
    if proxy_v2:
        proxy_script = (scripts_dir / 'proxy_populations_v2.py').resolve()
        min_ion_frac = cmfgen_cfg.get('proxy_populations_min_ion_frac', '1e-100')
        # t_exc_min: minimum excitation temperature for LTE populations.
        # Should match T_EXC in VADAT (set by prepare_VADAT from photosphere T).
        # Default 0.3 (3000K) is too low — CMFGEN uses T_EXC~1.0 (10000K) for
        # excitation, so proxy populations at 3000K are 30+ orders of magnitude
        # off from SE, causing 1e+90+ convergence corrections.
        t_exc_min = cmfgen_cfg.get('proxy_populations_t_exc_min', '1.0')
        proxy_cmd = [
            sys.executable, str(proxy_script), '.',
            '--hydro-data', './SN_HYDRO_DATA',
            '--output-dir', '.',
            '--strategy', 'GREY',
            '--min-ion-frac', str(min_ion_frac),
            '--t-exc-min', str(t_exc_min),
            '--write-t-ion',
            '--smooth',
            '--t-ion-source', 'stella',
        ]
        # Mass fraction floor: prevent catastrophic corrections for elements
        # with zero abundance at some depths (e.g. He in O/Si-rich core).
        mass_frac_floor = cmfgen_cfg.get('proxy_mass_frac_floor', '1e-8')
        proxy_cmd.extend(['--mass-frac-floor', str(mass_frac_floor)])
        # Complex model: enable velocity-aware tau, front-aware smoothing,
        # temperature-aware smoothing, adaptive ne damping, kramers fallback
        if cmfgen_cfg.get('proxy_complex_model', False):
            proxy_cmd.append('--complex-model')
        print(f"  Proxy command: {' '.join(proxy_cmd)}")
        proxy_result = subprocess.run(proxy_cmd, capture_output=True, text=True,
                                      cwd=str(cmfgen_dir))
        if proxy_result.stdout:
            print(f"  Proxy stdout (last 500): {proxy_result.stdout[-500:]}")
        if proxy_result.stderr:
            print(f"  Proxy stderr (last 500): {proxy_result.stderr[-500:]}")
        if proxy_result.returncode != 0:
            print(f"ERROR: proxy_populations_v2.py failed", file=sys.stderr)
            if proxy_result.stderr:
                print(proxy_result.stderr[-500:], file=sys.stderr)
            return False

        # Patch VADAT for proxy populations (same as SLURM script)
        # DC_METH selection: TR interpolates in T space, but REGRID_LOG_DC_V1
        # has a known issue with non-monotonic T (regrid_log_dc_v1.f:284:
        # "This code may have issues if T is very nonmonotonic"). With multiple
        # T matches per depth, the TR option picks the wrong interval, producing
        # pathological DCs and 10^22+ population corrections on iteration 1.
        # For non-monotonic T models, use DC_METH=R (radius interpolation) instead.
        vadat_path = cmfgen_dir / 'VADAT'
        if vadat_path.exists():
            # Check T monotonicity from SN_HYDRO_DATA
            dc_meth = 'TR'
            hydro_path = cmfgen_dir / 'SN_HYDRO_DATA'
            if hydro_path.exists():
                try:
                    t_arr = _read_hydro_section(hydro_path, 'Temperature')
                    if t_arr is not None and len(t_arr) > 1:
                        n_non_mono = int(np.sum(np.diff(t_arr) > 0))
                        if n_non_mono > len(t_arr) * 0.1:
                            dc_meth = 'R'
                            print(f"  [DC_METH] T is non-monotonic ({n_non_mono}/{len(t_arr)-1} zones), "
                                  f"using DC_METH=R instead of TR to avoid REGRID_LOG_DC_V1 interpolation failure")
                except Exception:
                    pass
            lines = vadat_path.read_text().splitlines()
            out = []
            for ln in lines:
                if '[LTE_EST]' in ln: out.append('F [LTE_EST]')
                elif '[LIN_INT]' in ln: out.append('T [LIN_INT]')
                elif '[DO_DDT]' in ln: out.append('F [DO_DDT]')
                elif '[USE_DJDT_RTE]' in ln: out.append('F [USE_DJDT_RTE]')
                elif '[DC_METH]' in ln: out.append(f'{dc_meth} [DC_METH]')
                else: out.append(ln)
            vadat_path.write_text('\n'.join(out) + '\n')

    print(f"  Proxy populations and VADAT prepared")

    # Phase 1: Run CMFGEN for 1 iteration to generate EDDFACTOR
    # (always needed — generates the correct frequency grid)
    #
    # Submit as a SLURM job to the interactive partition (even 1 iteration
    # can take several minutes for large models). We modify VADAT/IN_ITS
    # in-place, submit the job, wait for completion, then restore originals.
    print(f"  Preparing 1-iteration CMFGEN run for EDDFACTOR generation...")

    # Backup VADAT and IN_ITS so we can restore after the SLURM job
    vadat_path = cmfgen_dir / 'VADAT'
    in_its_path = cmfgen_dir / 'IN_ITS'
    vadat_backup = vadat_path.read_text() if vadat_path.exists() else ""
    in_its_backup = in_its_path.read_text() if in_its_path.exists() else ""

    # Modify IN_ITS: NUM_ITS=1
    in_its_lines = in_its_backup.splitlines()
    modified_its = []
    for line in in_its_lines:
        if '[NUM_ITS]' in line:
            modified_its.append("1 [NUM_ITS] ")
        else:
            modified_its.append(line)
    in_its_path.write_text("\n".join(modified_its) + "\n")

    # Modify VADAT: USE_FIXED_J=F, COMP_F=T
    vadat_lines = vadat_backup.splitlines()
    modified_v = []
    for line in vadat_lines:
        if "[USE_FIXED_J]" in line:
            modified_v.append("F [USE_FIXED_J]")
        elif "[COMP_F]" in line:
            modified_v.append("T [COMP_F]")
        else:
            modified_v.append(line)
    vadat_path.write_text("\n".join(modified_v) + "\n")

    # Remove existing EDDFACTOR (fresh start for COMP_F=T)
    edd_path = cmfgen_dir / 'EDDFACTOR'
    edd_info_path = cmfgen_dir / 'EDDFACTOR_INFO'
    if edd_path.exists():
        edd_path.unlink()
    if edd_info_path.exists():
        edd_info_path.unlink()

    # Create SLURM script for the 1-iteration EDDFACTOR generation run
    edd_slurm = cmfgen_dir / 'run_eddfactor_init.slurm'
    # SLURM job time limit (short — just the 1-iteration run)
    edd_hours = edd_slurm_time_min // 60
    edd_mins = edd_slurm_time_min % 60
    if edd_hours < 100:
        edd_time = f"{edd_hours:02d}:{edd_mins:02d}:00"
    else:
        edd_days = edd_hours // 24
        edd_hours_rem = edd_hours % 24
        edd_time = f"{edd_days:02d}-{edd_hours_rem:02d}:{edd_mins:02d}:00"
    edd_partition, edd_cpus, edd_mem = get_job_resources(
        debug,
        config_partition=cmfgen_cfg.get('init_partition', cmfgen_cfg.get('partition')),
        config_mem=cmfgen_cfg.get('init_mem'),
        config_cpus=cmfgen_cfg.get('cpus_per_task'),
    )
    model_name = config.get('model_name', 'model')
    edd_slurm.write_text(f"""#!/bin/bash
#SBATCH --job-name={model_name}_eddfactor_init
#SBATCH --output=eddfactor_%j.out
#SBATCH --error=eddfactor_%j.err
#SBATCH --time={edd_time}
#SBATCH --partition={edd_partition}
#SBATCH --qos=n0008
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={edd_cpus}
#SBATCH --mem={edd_mem}
#SBATCH --chdir={cmfgen_dir}

set -euo pipefail

echo "EDDFACTOR init job started on $(date)"

# Run CMFGEN via batch.sh (1 iteration, USE_FIXED_J=F, COMP_F=T)
tcsh -f {cmfgen_dir}/batch.sh

echo "EDDFACTOR init job finished on $(date)"
""")

    # Submit the SLURM job (with dependency if provided)
    submit_cmd = ['sbatch', '--parsable']
    if dependency_job_id:
        submit_cmd.extend(['--dependency=afterok:' + str(dependency_job_id)])
        print(f"  Submitting EDDFACTOR init job to partition={edd_partition}, time={edd_time}, dependency=afterok:{dependency_job_id}")
    else:
        print(f"  Submitting EDDFACTOR init job to partition={edd_partition}, time={edd_time}")
    submit_cmd.append(str(edd_slurm))
    submit_result = subprocess.run(submit_cmd, capture_output=True, text=True)
    if submit_result.returncode != 0:
        print(f"ERROR: sbatch failed for EDDFACTOR init job", file=sys.stderr)
        if submit_result.stderr:
            print(submit_result.stderr, file=sys.stderr)
        in_its_path.write_text(in_its_backup)
        vadat_path.write_text(vadat_backup)
        return False

    edd_job_id = submit_result.stdout.strip()
    print(f"  EDDFACTOR init job submitted: {edd_job_id}")

    # Wait for the job to complete (poll sacct every 15 seconds)
    import time as _time
    max_wait = timeout_min * 60 + 300  # timeout + 5 min buffer
    waited = 0
    poll_interval = 15
    job_done = False
    job_state = "UNKNOWN"
    while waited < max_wait:
        _time.sleep(poll_interval)
        waited += poll_interval
        sacct_cmd = ['sacct', '-j', edd_job_id, '--format=State', '--noheader', '--parsable2']
        sacct_result = subprocess.run(sacct_cmd, capture_output=True, text=True)
        states = [s.strip() for s in sacct_result.stdout.strip().split('\n') if s.strip()]
        if states:
            job_state = states[0]
            if job_state in ('COMPLETED', 'FAILED', 'CANCELLED', 'TIMEOUT', 'OUT_OF_MEMORY'):
                job_done = True
                break
        print(f"  Waiting for EDDFACTOR init job {edd_job_id} (state={job_state}, waited={waited}s)...")

    # Restore original IN_ITS and VADAT regardless of job outcome
    in_its_path.write_text(in_its_backup)
    vadat_path.write_text(vadat_backup)
    print(f"  Restored original IN_ITS and VADAT")

    if not job_done:
        print(f"ERROR: EDDFACTOR init job {edd_job_id} did not complete within {max_wait}s", file=sys.stderr)
        subprocess.run(['scancel', edd_job_id], capture_output=True)
        return False

    if job_state != 'COMPLETED':
        print(f"ERROR: EDDFACTOR init job {edd_job_id} ended with state={job_state}", file=sys.stderr)
        edd_out = cmfgen_dir / f'eddfactor_{edd_job_id}.out'
        if edd_out.exists():
            print(f"  Last 10 lines of output:", file=sys.stderr)
            for line in edd_out.read_text().splitlines()[-10:]:
                print(f"    {line}", file=sys.stderr)
        return False

    print(f"  EDDFACTOR init job {edd_job_id} completed successfully")

    # Verify EDDFACTOR was generated
    if not edd_path.exists() or not edd_info_path.exists():
        print(f"ERROR: EDDFACTOR not found after init job", file=sys.stderr)
        return False

    # Check FINISH_REC
    try:
        import struct as _struct
        sys.path.insert(0, str(scripts_dir))
        from generate_eddfactor import read_eddfactor_info
        nd_info, _, _, _ = read_eddfactor_info(edd_info_path)
        rec_size_info = 8 * (nd_info + 1)
        with open(edd_path, 'rb') as f:
            f.seek(4 * rec_size_info)
            finish_val = _struct.unpack('d', f.read(8))[0]
        if finish_val < 0.5:
            print(f"ERROR: EDDFACTOR incomplete (FINISH_REC={finish_val})", file=sys.stderr)
            return False
        print(f"  EDDFACTOR valid: FINISH_REC={finish_val}, ND={nd_info}")
    except Exception as e:
        print(f"ERROR: Failed to verify EDDFACTOR: {e}", file=sys.stderr)
        return False

    # Phase 2: If j_source="grey", replace J with W*B(T,nu)
    if j_source == 'grey':
        print(f"  Replacing J with W*B(T,nu) (j_source=grey)...")
        phase_cmd = [
            sys.executable, str(gen_script),
            str(cmfgen_dir),
            '--phase', 'replace',
        ]
        phase_result = subprocess.run(phase_cmd, capture_output=True, text=True,
                                      cwd=str(project_dir))
        print(phase_result.stdout)
        if phase_result.returncode != 0:
            print(f"ERROR: J replacement failed", file=sys.stderr)
            if phase_result.stderr:
                print(phase_result.stderr, file=sys.stderr)
            return False

    # Generate J comparison plot (always with --grey for comparison)
    if plot_comparison:
        plot_script = scripts_dir / 'plot_eddfactor_j.py'
        if plot_script.exists():
            plot_cmd = [sys.executable, str(plot_script), str(cmfgen_dir / 'EDDFACTOR'), '--grey']
            print(f"  Generating J comparison plot...")
            plot_result = subprocess.run(plot_cmd, capture_output=True, text=True,
                                         cwd=str(project_dir))
            if plot_result.returncode == 0:
                print(f"  J comparison plot generated")
            else:
                print(f"  WARNING: J plot generation failed (non-fatal)")
                if plot_result.stderr:
                    print(f"  {plot_result.stderr.strip()[:200]}")

    print(f"  EDDFACTOR generated successfully: {edd_path}")
    return True


def apply_fixed_j_workflow_overrides(config, j_release=False):
    """Apply VADAT overrides for the fixed-J workflow.

    For 'setup' phase (j_release=False):
      - USE_FIXED_J=T, COMP_F=F (read existing EDDFACTOR, don't compute new)
      - All other convergence settings from config

    For 'j_release' phase (j_release=True):
      - USE_FIXED_J=F, COMP_F=T (compute J self-consistently)
      - DO_LAM_AUTO=T (allow BA/full linearization now that populations are converged)
      - FIX_T=F, FIX_T_AUTO=T (allow T to vary for final convergence)
      - DO_NG=T (enable Ng acceleration for faster final convergence)

    Modifies config in-place.
    """
    cmfgen_cfg = config.get('cmfgen', {})
    fjw = cmfgen_cfg.get('fixed_j_workflow', {})
    if not fjw.get('enabled', False):
        return

    overrides = cmfgen_cfg.get('vadat_overrides', {})

    if j_release:
        # J-release phase: J solved self-consistently with converged populations
        print(f"\n  Fixed-J workflow: J-RELEASE phase")
        print(f"    USE_FIXED_J=F, COMP_F=T (J solved self-consistently)")
        print(f"    DO_LAM_AUTO=T (allow BA/full linearization)")
        print(f"    FIX_T=F, FIX_T_AUTO=T (T can vary for final convergence)")
        print(f"    DO_NG=T (Ng acceleration enabled)")
        overrides['USE_FIXED_J'] = 'F'
        overrides['COMP_F'] = 'T'
        overrides['DO_LAM_AUTO'] = 'T'
        overrides['FIX_T'] = 'F'
        overrides['FIX_T_AUTO'] = 'T'
        overrides['DO_NG'] = 'T'
    else:
        # Setup phase: fixed J from EDDFACTOR
        print(f"\n  Fixed-J workflow: SETUP phase")
        print(f"    USE_FIXED_J=T, COMP_F=F (read existing EDDFACTOR)")
        print(f"    DO_LAM_IT=T, DO_LAM_AUTO=T (Lambda + BA acceleration)")
        print(f"    JC_W_EDD=T, JBAR_W_EDD=T (keep EDDFACTOR open; RD_CONT_J handles freq mismatch)")
        overrides['USE_FIXED_J'] = 'T'
        overrides['COMP_F'] = 'F'
        overrides['DO_LAM_IT'] = 'T'
        overrides['DO_LAM_AUTO'] = 'T'

    cmfgen_cfg['vadat_overrides'] = overrides


def submit_j_release_job(config, scripts_dir, project_dir, cmfgen_dir,
                         fixed_j_job_id, config_file=None, day=None):
    """Submit the J-release CMFGEN job (stage 3 of the fixed-J workflow).

    This runs AFTER the fixed-J job has converged populations. It:
    - Modifies VADAT in-place ON THE COMPUTE NODE: USE_FIXED_J=F, COMP_F=T,
      DO_LAM_AUTO=T, FIX_T=F, FIX_T_AUTO=T, DO_NG=T
    - Preserves converged populations (does NOT regenerate proxy _IN files)
    - Submits with --dependency=afterok:<fixed_j_job_id>
    - Also submits the OBS job depending on the J-release job

    Returns dict with 'final' and 'obs' job IDs, or None on failure.
    """
    model_name = config.get('model_name', 'unknown')
    cmfgen_cfg = config.get('cmfgen', {})
    day_label = f"{day:.2f}days" if day is not None else None

    print(f"\n{'='*60}")
    print(f"Submitting J-RELEASE CMFGEN job (stage 3)")
    print(f"  Model: {model_name}, Day: {day}")
    print(f"  CMFGEN dir: {cmfgen_dir}")
    print(f"  Dependency: afterok:{fixed_j_job_id}")
    print(f"{'='*60}")

    # Create SLURM script for the J-release run.
    # The VADAT modification happens ON THE COMPUTE NODE (after the fixed-J
    # job has completed and written its output VADAT).
    partition, cpus, mem = get_job_resources(
        debug,
        config_partition=cmfgen_cfg.get('partition'),
        config_mem=cmfgen_cfg.get('mem'),
        config_cpus=cmfgen_cfg.get('cpus_per_task'),
    )
    time_limit = cmfgen_cfg.get('time', '01:55:00')
    obs_time = cmfgen_cfg.get('obs_time', time_limit)

    j_release_slurm = cmfgen_dir / 'run_cmfgen_j_release.slurm'
    j_release_slurm.write_text(f"""#!/bin/bash
#SBATCH --job-name={model_name}_j_release
#SBATCH --output=cmf_jrelease_%j.out
#SBATCH --error=cmf_jrelease_%j.err
#SBATCH --time={time_limit}
#SBATCH --partition={partition}
#SBATCH --qos=n0008
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={mem}
#SBATCH --chdir={cmfgen_dir}

set -euo pipefail

echo "J-RELEASE job started on $(date)"

# Copy converged populations from fixed-J run: _OUT -> _IN
# The fixed-J phase writes population output files as *_OUT.
# The J-release phase reads initial populations from *_IN.
# By copying _OUT to _IN, we use the converged populations from fixed-J
# as the starting point for the J-release phase.
echo "Copying converged populations (_OUT -> _IN)..."
for f in *_OUT; do
    if [ -f "$f" ]; then
        base="${{f%_OUT}}"
        cp "$f" "${{base}}_IN"
        echo "  Copied $f -> ${{base}}_IN"
    fi
done

# Also set LTE_EST=F and LIN_INT=T so CMFGEN reads the _IN files
# instead of computing LTE populations.
python3 -c "
import re
with open('VADAT', 'r') as f:
    lines = f.readlines()
out = []
for line in lines:
    if '[USE_FIXED_J]' in line: out.append('F [USE_FIXED_J]')
    elif '[COMP_F]' in line: out.append('T [COMP_F]')
    elif '[DO_LAM_AUTO]' in line: out.append('T [DO_LAM_AUTO]')
    elif '[FIX_T]' in line: out.append('F [FIX_T]')
    elif '[FIX_T_AUTO]' in line: out.append('T [FIX_T_AUTO]')
    elif '[DO_NG]' in line: out.append('T [DO_NG]')
    elif '[MAX_LAM]' in line: out.append('1.0D5 [MAX_LAM]')
    elif '[LAM_SCALE_OPT]' in line: out.append('NONE [LAM_SCALE_OPT]')
    elif '[SCALE_OPT]' in line: out.append('MAJOR [SCALE_OPT]')
    elif '[LTE_EST]' in line: out.append('F [LTE_EST]')
    elif '[LIN_INT]' in line: out.append('T [LIN_INT]')
    else: out.append(line.rstrip('\\n'))
with open('VADAT', 'w') as f:
    f.write('\\n'.join(out) + '\\n')
print('Modified VADAT for J-release: USE_FIXED_J=F, COMP_F=T, DO_LAM_AUTO=T, LTE_EST=F, LIN_INT=T')
"

# Run CMFGEN via batch.sh (J released, full linearization, converged pops)
tcsh -f {cmfgen_dir}/batch.sh

echo "J-RELEASE job finished on $(date)"
""")

    # Submit the J-release job with afterok dependency on the fixed-J job
    submit_cmd = ['sbatch', '--parsable',
                  '--dependency=afterok:' + str(fixed_j_job_id),
                  str(j_release_slurm)]
    print(f"  Submitting J-RELEASE job: partition={partition}, time={time_limit}")
    submit_result = subprocess.run(submit_cmd, capture_output=True, text=True)
    if submit_result.returncode != 0:
        print(f"ERROR: sbatch failed for J-RELEASE job", file=sys.stderr)
        if submit_result.stderr:
            print(submit_result.stderr, file=sys.stderr)
        return None

    j_release_job_id = submit_result.stdout.strip()
    print(f"  J-RELEASE job submitted: {j_release_job_id}")

    # Submit OBS job depending on the J-release job
    obs_job_id = None
    if not cmfgen_cfg.get('skip_obs', False):
        obs_slurm = cmfgen_dir / 'run_cmfgen_j_release_OBS.slurm'
        obs_slurm.write_text(f"""#!/bin/bash
#SBATCH --job-name={model_name}_j_release_OBS
#SBATCH --output=cmf_jrelease_obs_%j.out
#SBATCH --error=cmf_jrelease_obs_%j.err
#SBATCH --time={obs_time}
#SBATCH --partition={partition}
#SBATCH --qos=n0008
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={mem}
#SBATCH --chdir={cmfgen_dir}

set -euo pipefail

echo "J-RELEASE OBS job started on $(date)"

# Run CMFGEN OBS (formal solution for observed-frame spectrum)
# The OBS script is generated by submit_cmfgen_jobs.sh — we reuse it if it exists
if [ -f "{cmfgen_dir}/run_cmfgen_OBS.slurm" ]; then
    # Modify the existing OBS slurm to depend on j_release instead of FINAL
    # Just run batch.sh with OBS settings
    tcsh -f {cmfgen_dir}/batch.sh
else
    echo "WARNING: No OBS slurm template found, running batch.sh directly"
    tcsh -f {cmfgen_dir}/batch.sh
fi

echo "J-RELEASE OBS job finished on $(date)"

# Email guard: only notify on success when obs_fin_10 exists and is non-empty
OBS_FIN="{cmfgen_dir}/obs_fin_10"
if [[ ! -s "$OBS_FIN" ]]; then
    echo "ERROR: BATOBS did not produce obs_fin_10 or it is empty" >&2
    exit 3
fi

# Copy spectrum output if available
SPEC_OUT="{cmfgen_dir}/obs_fin_10"
NOTIFY_EMAIL="{cmfgen_cfg.get('notification_email', 'sbrennan@mpe.mpg.de')}"
printf "Subject: CMFGEN OBS completed for {model_name}\\n\\nCMFGEN OBS job completed successfully.\\nobs_fin_10 is non-empty.\\nSpectrum: %s\\n" "$SPEC_OUT" | mailx -s "CMFGEN OBS completed for {model_name}" "$NOTIFY_EMAIL" 2>/dev/null || true
echo "Notification email sent to $NOTIFY_EMAIL"
""")

        obs_submit_cmd = ['sbatch', '--parsable',
                          '--dependency=afterok:' + str(j_release_job_id),
                          str(obs_slurm)]
        obs_result = subprocess.run(obs_submit_cmd, capture_output=True, text=True)
        if obs_result.returncode == 0:
            obs_job_id = obs_result.stdout.strip()
            print(f"  J-RELEASE OBS job submitted: {obs_job_id}")
        else:
            print(f"  WARNING: OBS job submission failed (non-fatal)")
            if obs_result.stderr:
                print(f"  {obs_result.stderr.strip()[:200]}")

    return {
        'final': j_release_job_id,
        'obs': obs_job_id
    }




def submit_cmfgen_jobs(config, scripts_dir, project_dir, dependency_job_id=None, debug=False, config_file=None, no_relativistic=False, day=None, prev_cmfgen_dir=None, proxy_populations=False, proxy_populations_min_ion_frac=None, proxy_populations_grey_max_iter=None, proxy_populations_grey_tol=None):
    """Submit CMFGEN jobs (init, final, obs) using submit_cmfgen_jobs.sh script.

    Args:
        day: Optional snapshot day label for time-series runs (directory cmfgen/<day>days).
        prev_cmfgen_dir: Optional previous day's CMFGEN directory. If provided, the
            first CMFGEN run is skipped and the final job restarts from the previous
            day's converged populations.

    ── Diagnostic: submission chain and scripts invoked ─────────────────

    This function calls scripts/submit_cmfgen_jobs.sh, which on the compute
    node executes the following chain:

      1. scripts/prepare_VADAT.py
         — Reads SN_HYDRO_DATA to determine ND (number of depth points).
         — Updates MODEL_SPEC: ND = ND_hydro + 2*NBND_INS, NP = ND + NC.
         — Updates VADAT from VADAT_INIT + YAML vadat_overrides.
         — Auto-detects inner boundary tau: uses DIFFUSION (DIF=T) if
           tau>=50, ZERO_FLUX (DIF=F) if tau<50.
         — Sets RSTAR, LSTAR, VINF, VEL_LAW from SN_HYDRO_DATA.
         — When proxy_populations=True: generates dilute-LTE _IN files
           via proxy_populations.py (or proxy_populations_v2.py) and
           overrides VADAT for proxy mode (LTE_EST=F, LIN_INT=T, etc.).

      2. scripts/prepare_IN_INTS.py
         — Generates IN_INTS file for iteration control.

      3. batch.sh (from cmfgen_base/, patched for model_type)
         — Executes the CMFGEN binary (cmfgen_dev.exe for static_snapshot).
         — single_final=True: LTE cold-start + final convergence in one job.
         — single_final=False: INIT (cold start) + STABLE (convergence) stages.

      4. OBS job (dependency on FINAL)
         — Runs CMFGEN observation mode to extract the emergent spectrum.

    The CMFGEN job depends on the STELLA-to-CMFGEN conversion job
    (dependency_job_id) unless SN_HYDRO_DATA already exists, in which
    case the dependency is cleared (no need to wait for conversion).
    ─────────────────────────────────────────────────────────────────────
    """
    model_name = config.get('model_name', 'unknown')
    day_label = f"{day:.2f}days" if day is not None else None
    print(f"\n{'='*60}")
    if day_label:
        print(f"Submitting CMFGEN jobs for model: {model_name} day {day}")
    else:
        print(f"Submitting CMFGEN jobs for model: {model_name}")
    print(f"{'='*60}")

    # Extract paths
    output = config.get('output', {})
    model_dir = output.get('model_dir', '../models/${model_name}').replace('${model_name}', model_name)
    model_dir = (project_dir / model_dir).absolute()
    cmfgen_dir = model_dir / 'cmfgen' / day_label if day_label else model_dir / 'cmfgen'

    # Read all CMFGEN config from YAML (no env vars)
    cmfgen_cfg = config.get('cmfgen', {})
    init_num_its = cmfgen_cfg.get('init_num_its', 3)
    # single_final is always True now — INIT/STABLE stages removed.
    # The only distinction is cold start (single_final) vs warm restart (prev_cmfgen_dir set).
    single_final = True
    inc_rad_decays = cmfgen_cfg.get('inc_rad_decays', False)
    vel_law = cmfgen_cfg.get('vel_law', 12)
    model_type = cmfgen_cfg.get('model_type', 'static_snapshot')
    target_nd = cmfgen_cfg.get('target_nd', None)
    use_hydro_grid = cmfgen_cfg.get('use_hydro_grid', False)
    smooth_velocity_jumps = cmfgen_cfg.get('smooth_velocity_jumps', False)
    smooth_max_dv_ratio = cmfgen_cfg.get('smooth_max_dv_ratio', 2.0)
    smooth_n_insert = cmfgen_cfg.get('smooth_n_insert', 8)
    smooth_rho_sigma = cmfgen_cfg.get('smooth_rho_sigma', 0.0)
    smooth_temp_sigma = cmfgen_cfg.get('smooth_temp_sigma', 0.0)
    blend_composition = cmfgen_cfg.get('blend_composition', False)
    blend_sigma = cmfgen_cfg.get('blend_sigma', 5.0)
    blend_solar_weight = cmfgen_cfg.get('blend_solar_weight', 0.0)
    min_mass_fraction = cmfgen_cfg.get('min_mass_fraction', 0.0)
    monotonize_velocity = cmfgen_cfg.get('monotonize_velocity', False)
    chemical_n_insert = cmfgen_cfg.get('chemical_n_insert', 0)
    chemical_jump_threshold = cmfgen_cfg.get('chemical_jump_threshold', 0.1)
    reversal_n_insert = cmfgen_cfg.get('reversal_n_insert', 0)
    if target_nd is not None:
        print(f"[CMFGEN] target_nd={target_nd} set in config. "
              f"Converter will produce {target_nd}-point SN_HYDRO_DATA; "
              f"prepare_VADAT.py will set MODEL_SPEC ND={target_nd}+2*NBND_INS.")
    # SLURM resources
    # Debug mode: interactive partition, 4 CPUs, 30G memory.
    # Production mode: general partition, 8 CPUs, 32G memory.
    # Config values override production defaults but not debug presets.
    cmfgen_cpus = cmfgen_cfg.get('cpus_per_task', None)
    partition, cmfgen_cpus, slurm_mem = get_job_resources(
        debug,
        config_partition=cmfgen_cfg.get('partition'),
        config_mem=cmfgen_cfg.get('mem'),
        config_cpus=cmfgen_cpus if cmfgen_cpus else None,
    )
    is_interactive = debug or partition == 'interactive'
    if is_interactive:
        partition = 'interactive'
    time_limit = get_time_limit(debug, cmfgen_cfg.get('time'))
    init_slurm_mem = cmfgen_cfg.get('init_mem', slurm_mem)
    # Grid/solver
    n_ob_ins = cmfgen_cfg.get('n_ob_ins', 30)
    # FINAL requeue
    final_requeue_when = cmfgen_cfg.get('final_requeue_when', 'wall')
    final_requeue_max = cmfgen_cfg.get('final_requeue_max', 500)
    # Blackbody override
    auto_set_rstar = cmfgen_cfg.get('auto_set_rstar', False)
    auto_set_lstar = cmfgen_cfg.get('auto_set_lstar', False)
    phot_tau = cmfgen_cfg.get('phot_tau', 0.6666667)
    bb_teff_k = cmfgen_cfg.get('bb_teff_k', 9100)
    lstar_radius_mode = cmfgen_cfg.get('lstar_radius_mode', 'RMAX')

    # Submit script
    submit_script = scripts_dir / 'submit_cmfgen_jobs.sh'

    # Make script executable
    os.chmod(submit_script, 0o755)

    # Build command
    dep_id = dependency_job_id if dependency_job_id else ""
    config_file_arg = str(Path(config_file).resolve()) if config_file else ""
    cmfgen_partition = partition  # already set to 'interactive' if debug or config requests it
    cmfgen_time = time_limit
    cmd = [str(submit_script), model_name, str(model_dir), str(cmfgen_dir), dep_id,
           str(init_num_its),
           cmfgen_partition, cmfgen_time,
           '1' if single_final else '0',
           config_file_arg,
           '1' if inc_rad_decays else '0',
           str(vel_law),
           '1' if no_relativistic else '0',
           model_type]
    # Positional args 14-26 must be in fixed order for submit_cmfgen_jobs.sh:
    # 14=DAY_LABEL, 15=PREV_CMFGEN_DIR, 16=TARGET_ND,
    # 17=SLURM_MEM, 18=INIT_SLURM_MEM, 19=N_OB_INS,
    # 20=FINAL_REQUEUE_WHEN, 21=FINAL_REQUEUE_MAX,
    # 22=AUTO_SET_RSTAR, 23=AUTO_SET_LSTAR, 24=PHOT_TAU, 25=BB_TEFF_K, 26=LSTAR_RADIUS_MODE
    cmd.append(day_label if day_label else '')
    cmd.append(str(prev_cmfgen_dir) if prev_cmfgen_dir is not None else '')
    cmd.append(str(target_nd) if target_nd is not None else '')
    cmd.append(str(slurm_mem))
    cmd.append(str(init_slurm_mem))
    cmd.append(str(n_ob_ins))
    cmd.append(str(final_requeue_when))
    cmd.append(str(final_requeue_max))
    cmd.append('1' if auto_set_rstar else '0')
    cmd.append('1' if auto_set_lstar else '0')
    cmd.append(str(phot_tau))
    cmd.append(str(bb_teff_k))
    cmd.append(str(lstar_radius_mode))
    cmd.append('1' if proxy_populations else '0')
    cmd.append(str(proxy_populations_min_ion_frac) if proxy_populations_min_ion_frac is not None else '')
    # 29=PROXY_IN_NLTE_SOLVER
    proxy_populations_nlte_solver = cmfgen_cfg.get('proxy_populations_nlte_solver', False)
    cmd.append('1' if proxy_populations_nlte_solver else '0')
    # 30=PROXY_POPULATIONS_GREY_MAX_ITER, 31=PROXY_POPULATIONS_GREY_TOL
    cmd.append(str(proxy_populations_grey_max_iter) if proxy_populations_grey_max_iter is not None else "")
    cmd.append(str(proxy_populations_grey_tol) if proxy_populations_grey_tol is not None else "")
    # 32=OBS_TIME (optional override for OBS job walltime)
    obs_time = cmfgen_cfg.get('obs_time', '')
    cmd.append(str(obs_time) if obs_time else '')
    # 33=USE_HYDRO_GRID (bypass tau regridding in RV_SN_MODEL_SNIIN)
    cmd.append('1' if use_hydro_grid else '0')
    # 34=SMOOTH_VELOCITY_JUMPS, 35=SMOOTH_MAX_DV_RATIO, 36=SMOOTH_N_INSERT
    # 37=SMOOTH_RHO_SIGMA, 38=SMOOTH_TEMP_SIGMA
    # 39=BLEND_COMPOSITION, 40=BLEND_SIGMA, 41=BLEND_SOLAR_WEIGHT
    cmd.append('1' if smooth_velocity_jumps else '0')
    cmd.append(str(smooth_max_dv_ratio))
    cmd.append(str(smooth_n_insert))
    cmd.append(str(smooth_rho_sigma))
    cmd.append(str(smooth_temp_sigma))
    cmd.append('1' if blend_composition else '0')
    cmd.append(str(blend_sigma))
    cmd.append(str(blend_solar_weight))
    cmd.append(str(min_mass_fraction))
    cmd.append('1' if monotonize_velocity else '0')
    cmd.append(str(chemical_n_insert))
    cmd.append(str(chemical_jump_threshold))
    cmd.append(str(reversal_n_insert))
    # 47=PROXY_POPULATIONS_T_EXC_MIN (min excitation T for proxy LTE, 10^4 K)
    proxy_t_exc_min = cmfgen_cfg.get('proxy_populations_t_exc_min', '')
    cmd.append(str(proxy_t_exc_min) if proxy_t_exc_min else '')
    # 48=PROXY_POPULATIONS_V2 (use proxy_populations_v2.py instead of V1)
    proxy_v2_enabled = cmfgen_cfg.get('proxy_populations_v2', False)
    cmd.append('1' if proxy_v2_enabled else '0')
    # 49=CPUS_PER_TASK (CPU count for CMFGEN SLURM jobs)
    cmd.append(str(cmfgen_cpus) if cmfgen_cpus else '')
    # Pass complex model and mass frac floor as named CLI args
    if cmfgen_cfg.get('proxy_complex_model', False):
        cmd.append('--proxy-complex-model=1')
    proxy_mass_frac_floor = cmfgen_cfg.get('proxy_mass_frac_floor', '')
    if proxy_mass_frac_floor:
        cmd.append(f'--proxy-populations-mass-frac-floor={proxy_mass_frac_floor}')
    if debug:
        cmd.insert(1, '--debug')

    # Call the script (no env var modifications — all config via CLI args)
    # Pass SKIP_OBS via environment if configured
    env = os.environ.copy()
    cmfgen_install_dir = cmfgen_cfg.get('install_dir')
    if cmfgen_install_dir:
        cmfgen_install_dir = str(Path(cmfgen_install_dir).resolve())
        cmfgen_executable = Path(cmfgen_install_dir) / 'exe' / 'cmfgen_dev.exe'
        if not cmfgen_executable.is_file():
            print(f"ERROR: configured CMFGEN executable not found: {cmfgen_executable}")
            return None
        env['CMFDIST'] = cmfgen_install_dir
    if cmfgen_cfg.get('skip_obs', False):
        env['SKIP_OBS'] = '1'
    if cmfgen_cfg.get('proxy_populations_v2', False):
        env['USE_PROXY_V2'] = '1'
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        cwd=project_dir,
        env=env,
    )

    if result.returncode != 0:
        print(f"ERROR: CMFGEN jobs submission failed")
        print(result.stdout)
        print(result.stderr)
        return None

    # Parse job IDs from output (format: ":FINAL:OBS" — INIT removed)
    job_ids_str = result.stdout.strip().split('\n')[-1]
    parts = job_ids_str.split(':')
    # Format is now ":FINAL:OBS" (empty first field since INIT stage removed)
    final_job = parts[1] if len(parts) > 1 else None
    obs_job = parts[2] if len(parts) > 2 else None

    print(f"CMFGEN jobs submitted:")
    print(f"  FINAL: {final_job}")
    print(f"  OBS:   {obs_job}")

    return {
        'final': final_job,
        'obs': obs_job
    }


def submit_sedona_job(config, scripts_dir, project_dir, dependency_job_id=None, debug=False, day=None):
    """Submit SEDONA job to produce spectrum from existing hydro data.

    Args:
        day: Optional snapshot day label for time-series runs (directory sedona/<day>days).
    """
    model_name = config.get('model_name', 'unknown')
    day_label = f"{day:.2f}days" if day is not None else None
    print(f"\n{'='*60}")
    if day_label:
        print(f"Submitting SEDONA job for model: {model_name} day {day}")
    else:
        print(f"Submitting SEDONA job for model: {model_name}")
    print(f"{'='*60}")

    # Extract paths
    output_cfg = config.get('output', {})
    model_dir = output_cfg.get('model_dir', '../models/${model_name}').replace('${model_name}', model_name)
    model_dir = (project_dir / model_dir).absolute()
    sedona_dir = model_dir / 'sedona' / day_label if day_label else model_dir / 'sedona'
    
    # Get SEDONA settings from config
    sedona_cfg = config.get('sedona', {})
    heracles_cfg = config.get('heracles', {})  # DEPRECATED: HERACLES backend, use STELLA instead

    # Resource settings
    partition, _, memory = get_job_resources(
        debug,
        config_partition=sedona_cfg.get('partition', heracles_cfg.get('partition')),
        config_mem=sedona_cfg.get('memory', '128G'),
    )
    ntasks = sedona_cfg.get('ntasks', 8)
    time_limit = get_time_limit(debug, sedona_cfg.get('time'))
    
    # Check for required files (warn only - .mod will be created by STELLA/[DEPRECATED] HERACLES conversion)
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

    # If SN_HYDRO_DATA and .mod already exist, no dependency is needed
    cmfgen_dir = model_dir / 'cmfgen' / day_label if day_label else model_dir / 'cmfgen'
    sn_hydro = cmfgen_dir / 'SN_HYDRO_DATA'
    if sn_hydro.exists() and mod_file.exists():
        if dependency_job_id:
            print(f"  SN_HYDRO_DATA and .mod already exist — clearing dependency (no need to wait for conversion)")
            dependency_job_id = None
    
    # Create SEDONA slurm script from template
    job_name = f"{model_name}_{day_label}_sedona" if day_label else f"{model_name}_sedona"
    job_file = sedona_dir / f'sedona_{job_name}.slurm'

    # Create sedona directory if it doesn't exist
    sedona_dir.mkdir(parents=True, exist_ok=True)

    # Load template from sedona_base
    template_file = project_dir / 'sedona_base' / 'sedona_job_template.slurm'
    if not template_file.exists():
        print(f"ERROR: SEDONA job template not found at {template_file}")
        return None

    with open(template_file, 'r') as f:
        template = f.read()

    # Per-day CMFGEN SN_HYDRO_DATA for time-series
    cmfgen_data_dir = model_dir / 'cmfgen' / day_label if day_label else model_dir / 'cmfgen'

    # Output directory for centralized spectrum/plot collection
    outputs_rel = config.get('outputs', 'outputs')
    output_dir = (project_dir / outputs_rel).absolute()
    output_dir.mkdir(parents=True, exist_ok=True)

    # Epoch label for output filenames (e.g. "4.9d" for 5.0-day snapshot)
    if day is not None:
        epoch_label = f"{day:g}d"
    else:
        epoch_label = "0d"

    # Substitute variables
    # SEDONA direct-conversion fallback parameters (used by the slurm template
    # if the .mod file is missing at runtime — normally prepare_sedona_directory
    # creates it before submission).
    sedona_cfg_vars = config.get('sedona', {})
    sedona_target_day = day if day is not None else config.get('cmfgen_conversion', {}).get('snapshot_days', [40])[0]
    sedona_tau_max = sedona_cfg_vars.get('tau_max', 100.0)
    sedona_t_floor = sedona_cfg_vars.get('t_floor_K', 2000.0)
    # STELLA_SWD: prefer strad/ then stella/ (matches _find_stella_swd convention)
    swd_candidate = model_dir / 'stella' / 'strad' / f'{model_name}.swd'
    if not swd_candidate.exists():
        swd_candidate = model_dir / 'stella' / f'{model_name}.swd'
    variables = {
        'MODEL_NAME': job_name,
        'MODEL_BASE': model_name,
        'EPOCH_LABEL': epoch_label,
        'OUTPUT_DIR': str(output_dir),
        'PLOTS_DIR': str(model_dir / 'plots'),
        'TIME_LIMIT': time_limit,
        'PARTITION': partition,
        'NTASKS': ntasks,
        'MEMORY': memory,
        'WORK_DIR': str(sedona_dir),
        'STELLA_SWD': str(swd_candidate),
        'SEDONA_TARGET_DAYS': f"{sedona_target_day:g}",
        'SEDONA_TAU_MAX': f"{sedona_tau_max:g}",
        'SEDONA_T_FLOOR': f"{sedona_t_floor:g}",
        'CMFGEN_DATA': str(cmfgen_data_dir / 'SN_HYDRO_DATA'),
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


def clean_previous_superlite_runs(config, project_dir, force_cancel=True, verbose=True):
    """Cancel running SuperLite jobs and delete previous SuperLite output for a model.

    This ONLY affects SuperLite jobs and output — CMFGEN, STELLA, SEDONA, and SIROCCO
    jobs are never cancelled or touched.

    Args:
        config: YAML configuration dictionary
        project_dir: Path to project directory
        force_cancel: If True, cancel without prompting.
        verbose: If True, print detailed messages.

    Returns:
        Number of SuperLite jobs cancelled.
    """
    model_name = config.get('model_name', 'unknown')
    output = config.get('output', {})
    model_dir = output.get('model_dir', '../models/${model_name}').replace('${model_name}', model_name)
    model_dir = (project_dir / model_dir).absolute()
    superlite_root = model_dir / 'superlite'

    cancelled = 0

    if force_cancel:
        cancelled = cancel_existing_jobs(model_name, force=True, verbose=verbose, job_filter='superlite')

    import shutil
    import glob as _glob

    if superlite_root.exists():
        if verbose:
            print(f"\n  Cleaning previous SuperLite runs in: {superlite_root}")

        for subdir in list(superlite_root.iterdir()):
            if subdir.is_dir():
                # Remove output files but preserve input.str and input.par
                # (user may have customized grid/particle settings in input.par)
                for f in subdir.iterdir():
                    if f.name.startswith('output.') or f.name.startswith('superlite_'):
                        f.unlink()
                if verbose:
                    print(f"    Cleaned output files in: {subdir.name}")

        stale_patterns = ['superlite_*.out', 'superlite_*.err', 'superlite_*.slurm',
                          'flux.dat', 'spectrum.dat', 'fluxwl.dat', 'input.str',
                          '*.log']
        removed = 0
        for pattern in stale_patterns:
            for f in _glob.glob(str(superlite_root / pattern)):
                os.remove(f)
                removed += 1
        if verbose and removed > 0:
            print(f"    Removed {removed} stale file(s) from superlite root")

        try:
            superlite_root.rmdir()
            if verbose:
                print(f"    Removed empty superlite root directory")
        except OSError:
            pass

    return cancelled


def prepare_superlite_directory(config, scripts_dir, project_dir, day=None):
    """Prepare SuperLite directory by converting hydro data to input.str format.

    SuperLite reads a 1D spherical structure file (input.str) and a parameter
    file (input.par).  Two conversion paths are supported, selected by the
    ``superlite.input_source`` config key:

      * ``sn_hydro_data`` (default): CMFGEN SN_HYDRO_DATA -> input.str via
        convert_sn_hydro_to_superlite.py.  Requires the CMFGEN conversion to
        have run first.  Limited to 9 grouped abundance columns.
      * ``swd``: STELLA .swd -> input.str directly via
        convert_stella_swd_to_superlite.py.  Skips the CMFGEN step, uses the
        native STELLA Rosseland kappa for accurate tau truncation, passes all
      In debug mode, generates a quick-test input.par with reduced particles
      and wavelength groups for fast interactive iteration.
        15 elements individually, and applies MC-aware grid downsampling to
        reduce SuperLite runtime in optically-thick regions.

    Args:
        config: YAML configuration dictionary
        scripts_dir: Path to scripts directory
        project_dir: Path to project directory
        day: Optional snapshot day label. If provided, directory is
            superlite/<day>days for time-series runs.
    """
    model_name = config.get('model_name', 'unknown')
    debug = config.get('debug', False)
    day_label = f"{day:.2f}days" if day is not None else None
    print(f"\n{'='*60}")
    if day_label:
        print(f"Preparing SuperLite directory for model: {model_name} day {day}")
    else:
        print(f"Preparing SuperLite directory for model: {model_name}")
    print(f"{'='*60}")

    output = config.get('output', {})
    model_dir = output.get('model_dir', '../models/${model_name}').replace('${model_name}', model_name)
    model_dir = (project_dir / model_dir).absolute()

    superlite_dir = model_dir / 'superlite' / day_label if day_label else model_dir / 'superlite'
    superlite_dir.mkdir(parents=True, exist_ok=True)

    # Clean stale output files
    stale_patterns = ['superlite_*.out', 'superlite_*.err', 'flux.dat', 'spectrum.dat', 'fluxwl.dat']
    import glob as _glob
    removed_count = 0
    for pattern in stale_patterns:
        for stale_file in _glob.glob(str(superlite_dir / pattern)):
            os.remove(stale_file)
            removed_count += 1
    if removed_count > 0:
        print(f"  Cleaned {removed_count} stale output file(s) from previous run")

    # Read SuperLite settings from config
    sl_cfg = config.get('superlite', {})
    input_source = sl_cfg.get('input_source', 'sn_hydro_data')
    truncate_tau = sl_cfg.get('truncate_tau', None)
    nomp = sl_cfg.get('cpus_per_task', 4)

    input_str = superlite_dir / 'input.str'
    if input_str.exists():
        print(f"  Removing existing input.str to ensure fresh conversion")
        input_str.unlink()

    if input_source == 'swd':
        # --- Direct STELLA .swd -> SuperLite (MC-aware) ---
        return _prepare_superlite_from_swd(
            config, scripts_dir, model_dir, superlite_dir,
            day, truncate_tau, nomp, sl_cfg, debug)
    else:
        # --- CMFGEN SN_HYDRO_DATA -> SuperLite (legacy) ---
        return _prepare_superlite_from_sn_hydro(
            config, scripts_dir, model_dir, superlite_dir,
            day_label, truncate_tau, nomp)


def _find_stella_swd(model_dir, model_name):
    """Locate the STELLA .swd file for a model (strad/ then stella/)."""
    swd_file = model_dir / 'stella' / 'strad' / f'{model_name}.swd'
    if not swd_file.exists():
        swd_file = model_dir / 'stella' / f'{model_name}.swd'
    return swd_file


def _prepare_superlite_from_swd(config, scripts_dir, model_dir, superlite_dir,
                                day, truncate_tau, nomp, sl_cfg, debug=False):
    """Direct .swd -> input.str path (MC-aware, 15-element, native kappa).

    In debug mode, generates a quick-test input.par with reduced particles
    and wavelength groups for fast interactive iteration.
    """
    model_name = config.get('model_name', 'unknown')
    convert_script = scripts_dir / 'convert_stella_swd_to_superlite.py'
    if not convert_script.exists():
        print(f"  ERROR: Conversion script not found at {convert_script}")
        return False

    swd_file = _find_stella_swd(model_dir, model_name)
    if not swd_file.exists():
        print(f"  WARNING: STELLA .swd not found at {swd_file}")
        print(f"    The conversion will be attempted at job runtime as a fallback.")
        return True

    cmd = [sys.executable, str(convert_script), str(swd_file),
           '-o', str(superlite_dir), '--model-name', model_name]
    if day is not None:
        cmd.extend(['--target-days', str(day)])
    if truncate_tau:
        cmd.extend(['--truncate-tau', str(truncate_tau)])
    # MC-aware grid downsampling (default 150 cells; configurable)
    max_cells = sl_cfg.get('max_cells', 150)
    cmd.extend(['--max-cells', str(max_cells)])
    min_cells = sl_cfg.get('min_cells', 20)
    cmd.extend(['--min-cells', str(min_cells)])
    cmd.extend(['--nomp', str(nomp)])
    if sl_cfg.get('no_downsample', False):
        cmd.append('--no-downsample')
    rho_floor = sl_cfg.get('rho_floor', 0.0)
    if rho_floor > 0:
        cmd.extend(['--rho-floor', str(rho_floor)])
    t_floor = sl_cfg.get('t_floor_K', 2000.0)
    if t_floor > 0:
        cmd.extend(['--t-floor', str(t_floor)])
    # Grid method: adaptive (physics-aware AMR) or tau_weighted (legacy)
    grid_method = sl_cfg.get('grid_method', 'adaptive')
    cmd.extend(['--grid-method', grid_method])
    # AMR parameters (only used when grid_method=adaptive)
    if grid_method == 'adaptive':
        amr_cfg = sl_cfg.get('amr', {})
        if 'rho_grad' in amr_cfg:
            cmd.extend(['--amr-rho-grad', str(amr_cfg['rho_grad'])])
        if 'temp_grad' in amr_cfg:
            cmd.extend(['--amr-temp-grad', str(amr_cfg['temp_grad'])])
        if 'vel_grad' in amr_cfg:
            cmd.extend(['--amr-vel-grad', str(amr_cfg['vel_grad'])])
        if 'comp_grad' in amr_cfg:
            cmd.extend(['--amr-comp-grad', str(amr_cfg['comp_grad'])])
        if 'cannon_interface_cells' in amr_cfg:
            cmd.extend(['--amr-cannon-interface-cells', str(amr_cfg['cannon_interface_cells'])])
        if 'shock_dv_threshold' in amr_cfg:
            cmd.extend(['--amr-shock-dv-threshold', str(amr_cfg['shock_dv_threshold'])])
    # SuperLite numerical parameters (from superlite.numerical YAML section)
    num_cfg = sl_cfg.get('numerical', {})
    if 'grp_ng' in num_cfg:
        cmd.extend(['--sl-grp-ng', str(num_cfg['grp_ng'])])
    if 'grp_wlmin' in num_cfg:
        cmd.extend(['--sl-grp-wlmin', str(num_cfg['grp_wlmin'])])
    if 'grp_wlmax' in num_cfg:
        cmd.extend(['--sl-grp-wlmax', str(num_cfg['grp_wlmax'])])
    if 'flx_ndim' in num_cfg:
        cmd.extend(['--sl-flx-ndim', str(num_cfg['flx_ndim'])])
    if 'flx_wlmin' in num_cfg:
        cmd.extend(['--sl-flx-wlmin', str(num_cfg['flx_wlmin'])])
    if 'flx_wlmax' in num_cfg:
        cmd.extend(['--sl-flx-wlmax', str(num_cfg['flx_wlmax'])])
    if 'src_n2s' in num_cfg:
        cmd.extend(['--sl-src-n2s', str(num_cfg['src_n2s'])])
    if 'trn_tauddmc' in num_cfg:
        cmd.extend(['--sl-trn-tauddmc', str(num_cfg['trn_tauddmc'])])
    if 'taulump' in num_cfg:
        cmd.extend(['--sl-taulump', str(num_cfg['taulump'])])
    if 'nlte' in num_cfg:
        cmd.extend(['--sl-nlte', str(num_cfg['nlte'])])
    if 'nlte_nelem' in num_cfg:
        cmd.extend(['--sl-nlte-nelem', str(num_cfg['nlte_nelem'])])
    if 'niter' in num_cfg:
        cmd.extend(['--sl-niter', str(num_cfg['niter'])])
    if 'nebular' in num_cfg:
        cmd.extend(['--sl-nebular', str(num_cfg['nebular'])])
    # Quick-test mode in debug: reduced particles/groups for fast iteration
    if debug:
        cmd.append('--quick-test')

    print(f"  Converting STELLA .swd -> SuperLite (MC-aware, input_source=swd)")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode == 0:
        print(f"  Successfully created input.str and input.par (direct .swd path)")
        print(result.stdout)
    else:
        print(f"  ERROR: SuperLite .swd conversion failed")
        print(result.stdout)
        print(result.stderr)
        return False

    print(f"SuperLite directory prepared: {superlite_dir}")
    return True


def _prepare_superlite_from_sn_hydro(config, scripts_dir, model_dir, superlite_dir,
                                     day_label, truncate_tau, nomp):
    """Legacy SN_HYDRO_DATA -> input.str path (9 grouped abundances)."""
    # Find SN_HYDRO_DATA (from CMFGEN conversion)
    if day_label:
        sn_hydro = model_dir / 'cmfgen' / day_label / 'SN_HYDRO_DATA'
    else:
        raise ValueError("day is required — SN_HYDRO_DATA must live in cmfgen/<day>days/")

    if not sn_hydro.exists():
        print(f"  WARNING: SN_HYDRO_DATA not found at {sn_hydro}")
        print(f"    SuperLite requires SN_HYDRO_DATA to create input.str")
        print(f"    The conversion will be attempted at job runtime as a fallback.")
        return True

    convert_script = scripts_dir / 'convert_sn_hydro_to_superlite.py'
    if not convert_script.exists():
        print(f"  ERROR: Conversion script not found at {convert_script}")
        return False

    cmd = [sys.executable, str(convert_script), str(sn_hydro), '-o', str(superlite_dir)]
    if truncate_tau:
        cmd.extend(['--truncate-tau', str(truncate_tau)])
    cmd.extend(['--nomp', str(nomp)])

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode == 0:
        print(f"  Successfully created input.str and input.par")
        print(result.stdout)
    else:
        print(f"  ERROR: SuperLite conversion failed")
        print(result.stdout)
        print(result.stderr)
        return False

    print(f"SuperLite directory prepared: {superlite_dir}")
    return True


def submit_superlite_job(config, scripts_dir, project_dir, dependency_job_id=None, debug=False, day=None):
    """Submit SuperLite job to produce spectrum from existing hydro data.

    Args:
        day: Optional snapshot day label for time-series runs (directory superlite/<day>days).
    """
    model_name = config.get('model_name', 'unknown')
    day_label = f"{day:.2f}days" if day is not None else None
    print(f"\n{'='*60}")
    if day_label:
        print(f"Submitting SuperLite job for model: {model_name} day {day}")
    else:
        print(f"Submitting SuperLite job for model: {model_name}")
    print(f"{'='*60}")

    output_cfg = config.get('output', {})
    model_dir = output_cfg.get('model_dir', '../models/${model_name}').replace('${model_name}', model_name)
    model_dir = (project_dir / model_dir).absolute()
    superlite_dir = model_dir / 'superlite' / day_label if day_label else model_dir / 'superlite'

    sl_cfg = config.get('superlite', {})

    # Resource settings
    partition, sl_cpus, memory = get_job_resources(
        debug,
        config_partition=sl_cfg.get('partition'),
        config_mem=sl_cfg.get('memory', '32G'),
        config_cpus=sl_cfg.get('cpus_per_task'),
    )
    ntasks = sl_cfg.get('ntasks', 1)
    cpus_per_task = sl_cpus
    time_limit = get_time_limit(debug, sl_cfg.get('time'))

    # SuperLite paths
    superlite_install = sl_cfg.get('install_dir', '/ptmp/sbrennan/superlite')
    superlite_exe = sl_cfg.get('exe', f'{superlite_install}/superlite')
    superlite_data = f'{superlite_install}/Data'
    superlite_input = f'{superlite_install}/Input'

    # Check for input.str (warn only - conversion may happen at runtime)
    input_str = superlite_dir / 'input.str'
    if not input_str.exists():
        print(f"  WARNING: No input.str found in {superlite_dir}")
        print(f"    The SLURM job will attempt conversion from SN_HYDRO_DATA at runtime.")

    # If SN_HYDRO_DATA and input.str already exist, no dependency is needed
    cmfgen_dir = model_dir / 'cmfgen' / day_label if day_label else model_dir / 'cmfgen'
    sn_hydro = cmfgen_dir / 'SN_HYDRO_DATA'
    if sn_hydro.exists() and input_str.exists():
        if dependency_job_id:
            print(f"  SN_HYDRO_DATA and input.str already exist — clearing dependency (no need to wait for conversion)")
            dependency_job_id = None

    # Create SuperLite slurm script from template
    job_name = f"{model_name}_{day_label}_superlite" if day_label else f"{model_name}_superlite"
    job_file = superlite_dir / f'superlite_{job_name}.slurm'

    superlite_dir.mkdir(parents=True, exist_ok=True)
    (model_dir / 'plots').mkdir(parents=True, exist_ok=True)

    template_file = project_dir / 'superlite_base' / 'superlite_job_template.slurm'
    if not template_file.exists():
        print(f"ERROR: SuperLite job template not found at {template_file}")
        return None

    with open(template_file, 'r') as f:
        template = f.read()

    # Per-day CMFGEN SN_HYDRO_DATA for time-series
    cmfgen_data_dir = model_dir / 'cmfgen' / day_label if day_label else model_dir / 'cmfgen'

    # Output directory for centralized spectrum collection
    outputs_rel = config.get('outputs', 'outputs')
    output_dir = (project_dir / outputs_rel).absolute()
    output_dir.mkdir(parents=True, exist_ok=True)

    if day is not None:
        epoch_label = f"{day:g}d"
    else:
        epoch_label = "0d"

    variables = {
        'MODEL_NAME': job_name,
        'MODEL_BASE': model_name,
        'EPOCH_LABEL': epoch_label,
        'OUTPUT_DIR': str(output_dir),
        'PLOTS_DIR': str(model_dir / 'plots'),
        'TIME_LIMIT': time_limit,
        'PARTITION': partition,
        'NTASKS': ntasks,
        'CPUS_PER_TASK': cpus_per_task,
        'MEMORY': memory,
        'WORK_DIR': str(superlite_dir),
        'CMFGEN_DATA': str(cmfgen_data_dir / 'SN_HYDRO_DATA'),
        'SUPERLITE_EXE': superlite_exe,
        'SUPERLITE_DATA': superlite_data,
        'SUPERLITE_INPUT': superlite_input,
        # MC-aware direct .swd path (used by runtime fallback in slurm template)
        'INPUT_SOURCE': sl_cfg.get('input_source', 'sn_hydro_data'),
        'STELLA_SWD': str(_find_stella_swd(model_dir, model_name)),
        'TARGET_DAYS': epoch_label.rstrip('d'),
        'TRUNCATE_TAU': str(sl_cfg.get('truncate_tau', 20)),
        'MAX_CELLS': str(sl_cfg.get('max_cells', 150)),
        'T_FLOOR': str(sl_cfg.get('t_floor_K', 2000.0)),
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
        print(f"ERROR: SuperLite job submission failed")
        print(result.stdout)
        print(result.stderr)
        return None

    job_id = result.stdout.strip().split()[-1]
    print(f"SuperLite job submitted: {job_id}")
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
    heracles_cfg = config.get('heracles', {})  # DEPRECATED: HERACLES backend, use STELLA instead

    partition, cpus, memory = get_job_resources(
        debug,
        config_partition=sirocco_cfg.get('partition', heracles_cfg.get('partition')),
        config_mem=sirocco_cfg.get('memory', '16G'),
        config_cpus=sirocco_cfg.get('cpus'),
    )
    time_limit = get_time_limit(debug, sirocco_cfg.get('time'))

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


def submit_stella_job(config, scripts_dir, project_dir, debug=False, config_file=None, dependency=None):
    """Submit STELLA Slurm job (runs STELLA + converts swd -> SN_HYDRO_DATA)."""
    model_name = config.get('model_name', 'unknown')
    
    output = config.get('output', {})
    model_dir = output.get('model_dir', '../models/${model_name}').replace('${model_name}', model_name)
    model_dir = (project_dir / model_dir).resolve()

    stella_cfg = config.get('stella', {})
    heracles_cfg = config.get('heracles', {})  # DEPRECATED: HERACLES backend, use STELLA instead

    # Support arbitrary simulation days with adaptive stability measures.
    # For the time-series workflow, target_days is computed from the snapshot list
    # plus an extra margin unless the user explicitly set stella.target_days.
    stella_target = stella_cfg.get('target_days', None)
    heracles_target = heracles_cfg.get('target_days', 15)  # DEPRECATED: HERACLES fallback, use stella target_days
    if stella_target is not None:
        target_days = stella_target
    elif parse_snapshot_days(config) != [40.0]:
        target_days = get_stella_target_days(config)
    else:
        target_days = heracles_target
    
    nzon          = stella_cfg.get('nzon', 200)
    nzon_max      = stella_cfg.get('nzon_max', nzon)  # Allow AMR headroom if > nzon
    kappa         = stella_cfg.get('kappa', 0.34)
    tau_max       = stella_cfg.get('stella_tau_max', 30.0)
    tau_min       = stella_cfg.get('stella_tau_min', 0)
    # Convert to float if string (YAML may return scientific notation as string)
    if isinstance(tau_max, str):
        tau_max = float(tau_max)
    if isinstance(tau_min, str):
        tau_min = float(tau_min)
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
    r_min         = float(stella_cfg.get('r_min', 0.0))

    # Detect RSG-like models (extended stationary envelope) for robust defaults
    ejecta_cfg = config.get('ejecta', {})
    is_rsg = (ejecta_cfg.get('stationary_ejecta_A', False) and
              float(ejecta_cfg.get('R_star_cm', 0)) > 1e13)

    # RSG thermal bomb defaults: spread energy over more mass and longer time
    # to avoid the stiff shock that causes STELLA KFLAG=-2 / singular-Jacobian
    # convergence failures. Long runs (target_days > 30) need more gradual energy
    # injection and a lower timestep ceiling to keep the Gear solver stable.
    # These are only used if the user does not explicitly set the parameters.
    if is_rsg:
        explosion_energy_default = 5.0      # 5e50 erg (typical SN)
        if target_days > 30:
            explosion_mass_default   = 1.5      # 1.5 Msun (more mass than short-run 1.0)
            explosion_time_default   = 5.0      # 5 s (more gradual than 1 s)
        else:
            explosion_mass_default   = 1.0      # 1.0 Msun (more mass than 0.15)
            explosion_time_default   = 1.0      # 1 second (longer than 0.1s)
    else:
        explosion_energy_default = 0.0
        explosion_mass_default   = 0.0
        explosion_time_default   = 1.0

    explosion_energy = stella_cfg.get('explosion_energy', explosion_energy_default)
    explosion_mass   = stella_cfg.get('explosion_mass', explosion_mass_default)
    explosion_time   = stella_cfg.get('explosion_time', explosion_time_default)
    explosion_delay  = stella_cfg.get('explosion_delay', 0.0)
    initial_kinetic_energy = stella_cfg.get('initial_kinetic_energy', 'auto')  # EKO mode: 'auto' = compute from profile velocities, or explicit value (e.g. 4.0 = 4e50 erg)

    # Warn if the thermal bomb energy density is very high (common cause of KFLAG=-2)
    if isinstance(explosion_mass, (int, float)) and explosion_mass > 0 and isinstance(explosion_time, (int, float)) and explosion_time > 0:
        energy_density = explosion_energy / explosion_mass  # 1e50 erg / Msun
        power = explosion_energy / explosion_time  # 1e50 erg / s
        if is_rsg and (energy_density > 1e3 or power > 1e3):
            print(f"  WARNING: High RSG thermal bomb energy density")
            print(f"    Energy density: {energy_density:.1e} x 1e50 erg/Msun")
            print(f"    Power: {power:.1e} x 1e50 erg/s")
            print(f"    Consider increasing explosion_mass or explosion_time")

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

    # Robust solver defaults (validated on SN_12C_60days, 2026-08-06):
    # EPS=0.003, HMIN=1e-14, METH=2 (standard GEAR) are the STELLA reference
    # values (m030307mhomo.dat). Previous runs used EPS=0.03, HMIN=1e-10,
    # METH=3 (GEAR+BGH) which caused KFLAG=-2 crashes at td=3.46d because
    # HMIN=1e-10 was 10000x too large — the solver hit the floor and gave up
    # before reaching the tiny timesteps needed at the shock/reverse-shock
    # interface. METH=2 (standard GEAR) is more stable than METH=3 for
    # non-monotonic velocity profiles with CSM winds.
    eps_default = 0.003
    if is_rsg and target_days > 30:
        hmax_default = 43200.0  # 0.5 day cap for long RSG runs (stiffer late-time evolution)
        meth_default = 2        # Standard GEAR (METH=3 caused KFLAG=-2 crashes)
    elif is_rsg:
        hmax_default = 86400.0  # 1 day cap for short RSG runs
        meth_default = 2        # Standard GEAR for short RSG runs
    else:
        hmax_default = 5e5      # ~5.8 day cap for non-RSG
        meth_default = 2
    enable_convection_default = False  # Disabled due to numerical instability in innermost zone

    meth = int(stella_cfg.get('meth', meth_default))
    if meth not in (1, 2, 3):
        print(f"WARNING: meth={meth} not in {{1,2,3}}, resetting to default {meth_default}")
        meth = meth_default

    eps = float(stella_cfg.get('eps', eps_default))
    hmax = float(stella_cfg.get('hmax', hmax_default))
    hmin = float(stella_cfg.get('hmin', 1e-14))
    enable_convection = stella_cfg.get('convection',
                                       stella_cfg.get('enable_convection',
                                                      enable_convection_default))
    fit_tau = float(stella_cfg.get('fit_tau', 15.0))
    tau_tol = float(stella_cfg.get('tau_tol', 1.3))
    # Artificial viscosity for shock capturing.
    # Non-monotonic velocity profiles (e.g., inner ejecta faster than outer) create
    # shell collisions that require higher AQ to stabilize the GEAR solver.
    aq = float(stella_cfg.get('aq', 2.0))
    bq = float(stella_cfg.get('bq', 1.0))
    # Minimum relative zone spacing dr/r for AMR (prevents tiny zones that
    # force small timesteps). Default 0.01 = 1%.
    min_dr_r = float(stella_cfg.get('min_dr_r', 0.01))
    # Derive He mass fraction for Kramers free-free opacity tables.
    # Use wind composition if available (most representative of CSM),
    # fallback to shell A, then a solar default.
    ejecta_cfg2 = config.get('ejecta', {})
    fHe = float(ejecta_cfg2.get('fx_He_wind')
                or ejecta_cfg2.get('fx_He_A')
                or ejecta_cfg2.get('fx_He_B')
                or 0.28)
    fH  = float(ejecta_cfg2.get('fx_H_wind')
                or ejecta_cfg2.get('fx_H_A')
                or ejecta_cfg2.get('fx_H_B')
                or (1.0 - fHe - 0.02))  # fallback: solar metals ~2%
    abs_fraction_floor = stella_cfg.get('abs_fraction_floor', 1e-4)
    # Legacy abs_fraction: if explicitly set in YAML, pass as constant override;
    # otherwise leave as None so generate_stella_opacity_tables.py uses Kramers.
    abs_fraction_override = stella_cfg.get('abs_fraction', None)
    # OPAL Rosseland mean opacity tables (Iglesias & Rogers 1996).
    # When enabled, per-zone composition is read from the input_profile and
    # the SNEC OPAL tables (codataa..e) are used for hpbanab1.
    use_opal_tables  = stella_cfg.get('use_opal_tables', False)
    opal_table_dir   = stella_cfg.get('opal_table_dir',
                           str(scripts_dir.parent / 'SNEC' / 'SNEC')
                           if (scripts_dir.parent / 'SNEC' / 'SNEC' / 'tables').exists()
                           else str(Path('/ptmp/sbrennan/SNEC/SNEC')))
    input_profile_path = stella_cfg.get('input_profile', None)
    # ronfict: use STELLA's native vladsf opacity preprocessor (xronfict.exe) to generate
    # fully frequency-dependent multigroup opacity tables via Saha + bound-free/free-free/line
    # physics.  This supersedes the grey generate_stella_opacity_tables.py approach.
    # Requires: xronfict.exe, lineatom.dat, yakovlev/ cross-section data, and the eve .mod/.xni files.
    use_ronfict = stella_cfg.get('use_ronfict', False)
    # Resolve stella_root early so defaults for ronfict_exe / yakovlev_dir are correct.
    stella_root_for_ronfict = stella_cfg.get('stella_root', '/ptmp/sbrennan/stella/OpenStella')
    ronfict_exe = stella_cfg.get('ronfict_exe',
                      str(Path(stella_root_for_ronfict) / 'run' / 'vladsf' / 'xronfict.exe'))
    yakovlev_dir = stella_cfg.get('yakovlev_dir',
                      str(Path(stella_root_for_ronfict) / 'run' / 'vladsf' / 'yakovlev'))
    givdtl = stella_cfg.get('givdtl', True)
    custom_output_times = stella_cfg.get('custom_output_times', None)
    iout = stella_cfg.get('iout', 10)  # Output frequency: 1=every timestep, 10=every 10 timesteps
    outer_fraction = stella_cfg.get('outer_fraction', 0.3)
    transition_radius = stella_cfg.get('transition_radius', 1e13)
    
    # Wacc parameters for temperature convergence control
    wacc_r = float(stella_cfg.get('wacc_r', 1.0))
    wacc_v = float(stella_cfg.get('wacc_v', 1.0))
    wacc_t = float(stella_cfg.get('wacc_t', 0.5))  # Default 0.5 for stability with heating
    wacc_rad = float(stella_cfg.get('wacc_rad', 1.0))
    temp_floor = float(stella_cfg.get('temp_floor', 0.03))  # Default 0.03 = 300 K (STELLA reference)
    floor_rad = float(stella_cfg.get('floor_rad', 0.01))  # Default 0.01 = 100 K (STELLA reference)
    
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
        method_name = {1: 'Adams', 2: 'GEAR', 3: 'GEAR+BGH'}.get(meth, 'Unknown')
        print(f"  [RSG mode] EPS={eps}, HMAX={hmax:.0f}s, METH={meth} ({method_name}), CONV={'T' if enable_convection else 'F'}")
    if should_use_multistage(target_days):
        stages = get_multistage_plan(target_days)
        print(f"  Multi-stage simulation plan:")
        for i, stage in enumerate(stages):
            print(f"    Stage {i+1}: {stage['start']:.1f}-{stage['end']:.1f} days ({stage['name']})")
            print(f"      NSTMAX: {stage['nstmax']}, Temp floor: {stage['temp_floor']}")
    
    cmfgen_conv_cfg  = config.get('cmfgen_conversion', {})
    cmfgen_snapshot_days = parse_snapshot_days(config)
    cmfgen_tau_max = cmfgen_conv_cfg.get('cmfgen_tau_max', tau_max)
    cmfgen_tau_min = cmfgen_conv_cfg.get('cmfgen_tau_min', tau_min)
    if isinstance(cmfgen_tau_max, str):
        cmfgen_tau_max = float(cmfgen_tau_max)
    if isinstance(cmfgen_tau_min, str):
        cmfgen_tau_min = float(cmfgen_tau_min)
    enable_grid_thinning = cmfgen_conv_cfg.get('enable_grid_thinning', False)
    grid_thinning_threshold = cmfgen_conv_cfg.get('grid_thinning_threshold', 0.005)
    # Adaptive grid parameters
    enable_adaptive_grid = cmfgen_conv_cfg.get('enable_adaptive_grid', True)
    base_grid_target = cmfgen_conv_cfg.get('base_grid_target', 50)
    min_points = cmfgen_conv_cfg.get('min_points', 150)
    max_radius_jump = cmfgen_conv_cfg.get('max_radius_jump', 0.05)
    max_temp_jump = cmfgen_conv_cfg.get('max_temp_jump', 0.10)
    max_density_jump = cmfgen_conv_cfg.get('max_density_jump', 0.15)
    max_ne_jump = cmfgen_conv_cfg.get('max_ne_jump', 0.10)
    max_logtau_jump = cmfgen_conv_cfg.get('max_logtau_jump', 0.10)
    max_velocity_jump = cmfgen_conv_cfg.get('max_velocity_jump', 0.10)
    max_sigma_jump = cmfgen_conv_cfg.get('max_sigma_jump', 1.0)
    max_cell_ratio = cmfgen_conv_cfg.get('max_cell_ratio', 1.5)
    max_points = cmfgen_conv_cfg.get('max_points', 150)
    min_nd = cmfgen_conv_cfg.get('min_nd', 150)
    smooth_physical_discontinuities = cmfgen_conv_cfg.get('smooth_physical_discontinuities', False)
    phys_smooth_sigma_zones = cmfgen_conv_cfg.get('phys_smooth_sigma_zones', 10.0)
    stretch_trouble_regions = cmfgen_conv_cfg.get('stretch_trouble_regions', False)
    stretch_factor = cmfgen_conv_cfg.get('stretch_factor', 3.0)
    max_stretch_zones = cmfgen_conv_cfg.get('max_stretch_zones', 200)
    min_zones_to_stretch = cmfgen_conv_cfg.get('min_zones_to_stretch', 2)
    # Velocity smoothing parameters
    smooth_velocity = cmfgen_conv_cfg.get('smooth_velocity', False)
    sigma_threshold = cmfgen_conv_cfg.get('sigma_threshold', 2.0)
    smooth_width = cmfgen_conv_cfg.get('smooth_width', 5)
    # Negative-sigma coarsening (Luc Dessart approach)
    coarsen_negative_sigma = cmfgen_conv_cfg.get('coarsen_negative_sigma', True)
    # DEPRECATED: heracles_cfg fallback kept for backward compat; use stella_cfg instead
    partition, _, memory = get_job_resources(
        debug,
        config_partition=stella_cfg.get('partition', heracles_cfg.get('partition')),
        config_mem=stella_cfg.get('memory', heracles_cfg.get('memory')),
    )
    time_limit = get_time_limit(debug, stella_cfg.get('time'))

    stella_root   = stella_cfg.get('stella_root', '/ptmp/sbrennan/stella/OpenStella')
    input_profile = stella_cfg.get('input_profile',
                        str(model_dir / 'stella' / 'input_profile'))

    # Compute kinetic energy from profile and set EKO automatically.
    # EKO (initial_kinetic_energy) in STELLA's .dat file:
    #   EKO > 0: STELLA OVERWRITES all .hyd velocities with a triangular
    #            (linear-in-mass) profile scaled to match EKO. The input
    #            velocity structure is completely lost. The mass_fract_tri
    #            parameter (AI1) only controls WHERE the triangular starts
    #            (0=innermost zone), NOT whether it overwrites.
    #            See begtt.f / begradsep.f lines 218-258.
    #   EKO = 0: STELLA preserves the velocities from the .hyd file as-is.
    #            Use this when the input profile already has correct velocities.
    # If initial_kinetic_energy is "auto", EKO is computed from the input
    # profile velocities (WARNING: this overwrites them!). Set to 0.0 to
    # preserve the .hyd velocity profile.
    if Path(input_profile).exists():
        computed_ke = calculate_kinetic_energy_from_profile(input_profile, r_min=r_min)
        if computed_ke is not None:
            if initial_kinetic_energy == 'auto':
                initial_kinetic_energy = computed_ke
                print(f"  EKO set from profile KE: {initial_kinetic_energy:.3f} x1e50 erg (auto)")
                print(f"  WARNING: EKO>0 causes STELLA to OVERWRITE input velocities with a triangular profile!")
                print(f"           Set initial_kinetic_energy: 0.0 to preserve .hyd velocities.")
            elif isinstance(initial_kinetic_energy, (int, float)) and initial_kinetic_energy > 0.0:
                print(f"  EKO={initial_kinetic_energy:.3f} x1e50 erg (profile KE={computed_ke:.3f})")
                print(f"  WARNING: EKO>0 causes STELLA to OVERWRITE input velocities with a triangular profile!")
                print(f"           The mass_fract_tri=0 does NOT preserve velocities - it just starts triangular from zone 1.")
                if abs(initial_kinetic_energy - computed_ke) > 0.5:
                    print(f"  WARNING: EKO differs from profile KE by {abs(initial_kinetic_energy - computed_ke):.3f} x1e50 erg")
                    print(f"           Consider setting initial_kinetic_energy to 'auto' to match profile")
            else:
                print(f"  INFO: Profile KE = {computed_ke:.3f} x 1e50 erg (EKO=0, velocities from .hyd preserved)")

    stella_run_dir = model_dir / 'stella'
    stella_run_dir.mkdir(parents=True, exist_ok=True)

    # Generate .xni file if nickel56 is enabled
    nickel_cfg = stella_cfg.get('nickel56', {})
    if nickel_cfg.get('enabled', False):
        print(f"Generating .xni file from nickel56 configuration...")
        xni_script = scripts_dir / 'generate_xni_file.py'
        xni_output = stella_run_dir / 'eve' / f'{model_name}.xni'
        xni_output.parent.mkdir(parents=True, exist_ok=True)

        # Pass input_profile path for mass-based Ni56 distribution
        input_profile_path = stella_run_dir / 'input_profile'
        xni_args = [sys.executable, str(xni_script), str(config_file), str(xni_output)]
        if input_profile_path.exists():
            xni_args.append(str(input_profile_path))

        result = subprocess.run(
            xni_args,
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
        # nickel56 not enabled: generate a ZERO .xni file to prevent eve2 from
        # creating one with default Ni56 (which would add spurious radioactive decay).
        print("nickel56 not enabled, generating ZERO .xni file (no radioactive decay)")
        xni_output = stella_run_dir / 'eve' / f'{model_name}.xni'
        xni_output.parent.mkdir(parents=True, exist_ok=True)
        nzon_for_xni = int(stella_cfg.get('nzon', 500))
        # .xni format: Fortran unformatted sequential record of nzon REAL*8 values.
        # STELLA reads: read(28) (XNifor(i), i=1, Nzon) where XNifor is REAL*8.
        # All Ni56 fractions = 0.0 means no radioactive decay.
        import struct as _struct
        rec_len = nzon_for_xni * 8  # nzon REAL*8 (double precision) values
        with open(xni_output, 'wb') as f:
            f.write(_struct.pack('i', rec_len))           # record marker
            f.write(_struct.pack(f'{nzon_for_xni}d', *([0.0] * nzon_for_xni)))  # data
            f.write(_struct.pack('i', rec_len))           # trailing marker
        print(f"  Written zero .xni file: {xni_output} ({4 + rec_len + 4} bytes, {nzon_for_xni} zones, Fortran unformatted REAL*8)")

    template_file = scripts_dir / 'stella_job_template.slurm'
    with open(template_file, 'r') as f:
        template = f.read()

    # Convert custom_output_times list to comma-separated string for template
    custom_output_times_str = ','.join(str(t) for t in custom_output_times) if custom_output_times else ''
    
    plots_dir = model_dir / 'plots'
    plots_dir.mkdir(parents=True, exist_ok=True)

    variables = {
        'MODEL_NAME':    model_name,
        'MODEL_DIR':     str(model_dir),
        'PROJECT_DIR':   str(project_dir),
        'PLOTS_DIR':     str(plots_dir),
        'SCRIPTS_DIR':   str(scripts_dir),
        'STELLA_ROOT':   stella_root,
        'STELLA_RUN_DIR': str(stella_run_dir),
        'INPUT_PROFILE': input_profile,
        'TARGET_DAYS':   target_days,
        'NZON':          nzon,
        'NZON_MAX':      nzon_max,
        'ZERO_NI56':     'true' if not nickel_cfg.get('enabled', False) else 'false',
        'KEPLER_MODE':   'true' if (config.get('kepler_presn', None) or
                                     config.get('sn_b_from_sn_a', None) or
                                     config.get('x_model_from_sn', None) or
                                     config.get('stella_from_swd_snapshot', None)) else 'false',
        'KAPPA':         kappa,
        'TAU_MAX':       tau_max,
        'TAU_MIN':       tau_min,
        'RHO_MIN_ARG':   f'--rho-min="{rho_min}"' if rho_min > 0 else '',
        'SNAPSHOT_INTERVAL': snapshot_interval,
        'NI56_DIST_MASS': nickel_cfg.get('distribution_mass_msun', 0.0),
        'NI56_TOTAL_MASS': nickel_cfg.get('total_mass_msun', 0.0),
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
        'HMIN': hmin,
        'METH': meth,
        'ENABLE_CONV_ARG': '--enable-convection' if enable_convection else '',
        'FIT_TAU': fit_tau,
        'TAU_TOL': tau_tol,
        'AQ': aq,
        'BQ': bq,
        'MIN_DR_R': min_dr_r,
        'FH': fH,
        'FHE': fHe,
        'ABS_FRACTION_FLOOR': abs_fraction_floor,
        'ABS_FRACTION_ARG': f'--abs-fraction {abs_fraction_override}' if abs_fraction_override is not None else '',
        'OPAL_TABLES_ARG': (
            f'--opal-tables '
            f'--opal-table-dir "{opal_table_dir}" '
            f'--input-profile "{input_profile_path or input_profile}"'
        ) if use_opal_tables else '',
        'GIVDTL_ARG': '--givdtl' if givdtl else '',
        'RESTART_STELLA': 'true' if config.get('restart_stella', False) else 'false',
        # In the time-series workflow (snapshot_days is a list), the slurm
        # template's Step 8 conversion is redundant — per-day SN_HYDRO_DATA
        # files are created by submit_stella_post_conversion_job().  Use the
        # first snapshot day here so the template's single --target-days call
        # gets a valid float instead of a Python list repr like "[5, 10, ...]".
        'CMFGEN_SNAPSHOT_DAYS': (cmfgen_snapshot_days[0]
                                 if isinstance(cmfgen_snapshot_days, list)
                                 else cmfgen_snapshot_days),
        # Day label for cmfgen/<day>days/ directory structure.
        # Single-day: "30.0days" -> cmfgen/30.0days/
        # Time-series: use first day as fallback (per-day dirs created by post-conversion job)
        'CMFGEN_DAY_LABEL': f"{cmfgen_snapshot_days[0]:.2f}days"
                            if isinstance(cmfgen_snapshot_days, list)
                            else f"{cmfgen_snapshot_days:.2f}days",
        'CMFGEN_TAU_MAX':       cmfgen_tau_max,
        'CMFGEN_TAU_MIN':       cmfgen_tau_min,
        'ENABLE_GRID_THINNING_ARG': '--enable-grid-thinning' if enable_grid_thinning else '',
        'GRID_THINNING_THRESHOLD': grid_thinning_threshold,
        'ENABLE_ADAPTIVE_GRID_ARG': '--enable-adaptive-grid' if enable_adaptive_grid else '',
        'PHYS_SMOOTH_ARG': '--smooth-physical-discontinuities --phys-smooth-sigma-zones {}'.format(phys_smooth_sigma_zones) if smooth_physical_discontinuities else '',
        'STRETCH_ARG': '--stretch-trouble-regions --stretch-factor {} --max-stretch-zones {} --min-zones-to-stretch {}'.format(stretch_factor, max_stretch_zones, min_zones_to_stretch) if stretch_trouble_regions else '',
        'BASE_GRID_TARGET': base_grid_target,
        'MIN_POINTS': min_points,
        'MAX_RADIUS_JUMP': max_radius_jump,
        'MAX_TEMP_JUMP': max_temp_jump,
        'MAX_DENSITY_JUMP': max_density_jump,
        'MAX_NE_JUMP': max_ne_jump,
        'MAX_LOGTAU_JUMP': max_logtau_jump,
        'MAX_POINTS': max_points,
        'MIN_ND': min_nd,
        'SMOOTH_VELOCITY_ARG': '--smooth-velocity' if smooth_velocity else '',
        'SIGMA_THRESHOLD': sigma_threshold,
        'SMOOTH_WIDTH': smooth_width,
        'COARSEN_NEGATIVE_SIGMA_ARG': '--no-coarsen-negative-sigma' if not coarsen_negative_sigma else '',
        'WACC_R': wacc_r,
        'WACC_V': wacc_v,
        'WACC_T': wacc_t,
        'WACC_RAD': wacc_rad,
        'TEMP_FLOOR': temp_floor,
        'FLOOR_RAD': floor_rad,
        'CUSTOM_OUTPUT_TIMES_ARG': f'--custom-output-times="{custom_output_times_str}"' if custom_output_times else '',
        'IOUT_ARG': f'--iout={iout}',
        'NOUT_ARG': f'--nout=100',  # NOUT=100: .swd write every 100 steps (covers 0-50d)
        'OUTER_FRACTION': outer_fraction,
        'TRANSITION_RADIUS': transition_radius,
        'PARTITION':     partition,
        'MEMORY':        memory,
        'TIME':          time_limit,
        'CONFIG_FILE':   str(config_file) if config_file else '',
        # ronfict multigroup opacity table generation
        'USE_RONFICT':   'true' if use_ronfict else 'false',
        'RONFICT_EXE':   ronfict_exe,
        'YAKOVLEV_DIR':  yakovlev_dir,
    }

    job_script = substitute_template(template, variables)
    job_file = stella_run_dir / f'stella_{model_name}.slurm'
    with open(job_file, 'w') as f:
        f.write(job_script)

    sbatch_cmd = ['sbatch']
    if dependency:
        sbatch_cmd.extend(['--dependency', dependency])
        print(f"  STELLA job dependency: {dependency}")
    sbatch_cmd.append(str(job_file))

    result = subprocess.run(
        sbatch_cmd,
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
    """[DEPRECATED] Check if SN_HYDRO_DATA (conversion output) exists in the cmfgen directory.

    Name retained for backward compat; works for both STELLA and HERACLES pipelines.
    """
    model_name = config.get('model_name', 'unknown')
    output = config.get('output', {})
    model_dir = output.get('model_dir', '../models/${model_name}').replace('${model_name}', model_name)
    model_dir = (project_dir / model_dir).absolute()
    sn_hydro = model_dir / 'cmfgen' / 'SN_HYDRO_DATA'
    return sn_hydro.exists(), sn_hydro


def run_stella_conversion(config, scripts_dir, project_dir, source_model_name=None, day=None):
    """Run STELLA .swd -> SN_HYDRO_DATA conversion if output is missing.

    Args:
        source_model_name: If specified, use .swd file from this model name (for reuse_stella case).
                          When reuse_stella=True, this is the current model name (uses existing output).
        day: Optional snapshot day label. If provided, output goes to
            cmfgen/<day>days/SN_HYDRO_DATA for time-series runs.

    ── Diagnostic: script invoked ───────────────────────────────────────

    Calls build_cmfgen_conversion_cmd() which invokes:
      scripts/convert_stella_swd_to_cmfgen.py (v1, default)
      scripts/convert_stella_swd_to_cmfgen_v2.py (v2)

    See build_cmfgen_conversion_cmd() docstring for full details on
    converter selection, base_grid_method options, and sigma-preserving
    config settings.
    ─────────────────────────────────────────────────────────────────────
    """
    model_name = config.get('model_name', 'unknown')
    output = config.get('output', {})
    model_dir = output.get('model_dir', '../models/${model_name}').replace('${model_name}', model_name)
    model_dir = (project_dir / model_dir).absolute()

    day_label = f"{day:.2f}days" if day is not None else None

    # Ensure cmfgen directory exists for output
    cmfgen_dir = model_dir / 'cmfgen' / day_label if day_label else model_dir / 'cmfgen'
    cmfgen_dir.mkdir(parents=True, exist_ok=True)

    # Ensure plots directory exists for diagnostic plots
    plots_dir = model_dir / 'plots'
    plots_dir.mkdir(parents=True, exist_ok=True)

    sn_hydro = cmfgen_dir / 'SN_HYDRO_DATA'
    if sn_hydro.exists():
        print(f"  SN_HYDRO_DATA already exists: {sn_hydro}")
        print(f"  Removing existing file to regenerate with current tau settings")
        sn_hydro.unlink()

    # Use source model name for .swd file if reusing STELLA output
    swd_model_name = source_model_name if source_model_name else model_name

    # Find .swd file (check strad directory first, then main stella directory)
    swd_file = model_dir / 'stella' / 'strad' / f'{swd_model_name}.swd'
    if not swd_file.exists():
        swd_file = model_dir / 'stella' / f'{swd_model_name}.swd'
    if not swd_file.exists():
        print(f"ERROR: STELLA .swd not found for model {model_name}")
        print(f"  Checked: {model_dir}/stella/strad/{swd_model_name}.swd")
        print(f"  Checked: {model_dir}/stella/{swd_model_name}.swd")
        return False

    # Find .rho file (final evolved composition from STELLA)
    rho_model_name = source_model_name if source_model_name else model_name
    rho_file = model_dir / 'stella' / 'eve' / f'{rho_model_name}.rho'
    if not rho_file.exists():
        rho_file = None

    # Find .abn file (initial composition, fallback if .rho not available)
    abn_model_name = source_model_name if source_model_name else model_name
    abn_file = model_dir / 'stella' / 'modmake' / f'{abn_model_name}.abn'
    if not abn_file.exists():
        abn_file = None

    # Find input_profile for composition fallback (if .rho/.abn has zero Ca)
    input_profile = model_dir / 'stella' / 'input_profile'
    if not input_profile.exists():
        input_profile = None

    # Find .xni file (Ni56 distribution from STELLA eve2)
    xni_model_name = source_model_name if source_model_name else model_name
    xni_file = model_dir / 'stella' / 'eve' / f'{xni_model_name}.xni'
    if not xni_file.exists():
        xni_file = None

    # Get conversion parameters from config
    stella_cfg = config.get('stella', {})
    cmfgen_conv_cfg = config.get('cmfgen_conversion', {})
    cmfgen_cfg_rad = config.get('cmfgen', {})
    inc_rad_decays = cmfgen_cfg_rad.get('inc_rad_decays', False)
    converter_version = cmfgen_conv_cfg.get('converter', 'v1')
    converter_script = 'convert_stella_swd_to_cmfgen_v2.py' if converter_version == 'v2' else 'convert_stella_swd_to_cmfgen.py'
    target_days = day if day is not None else cmfgen_conv_cfg.get('snapshot_days', stella_cfg.get('target_days', 40))

    kappa = stella_cfg.get('kappa', 0.34)
    # cmfgen_tau_max/cmfgen_tau_min for CMFGEN conversion: read from cmfgen_conversion section first,
    # fall back to stella section.  The stella.stella_tau_min controls pre-run .hyd trimming
    # and is typically much larger (0.1) than the value needed for CMFGEN grid selection.
    tau_max = cmfgen_conv_cfg.get('cmfgen_tau_max', stella_cfg.get('stella_tau_max', 100.0))
    tau_min = cmfgen_conv_cfg.get('cmfgen_tau_min', stella_cfg.get('stella_tau_min', 0))
    # Convert to float if string (YAML may return scientific notation as string)
    if isinstance(tau_max, str):
        tau_max = float(tau_max)
    if isinstance(tau_min, str):
        tau_min = float(tau_min)
    if tau_max < 50:
        print(f"  WARNING: cmfgen_tau_max={tau_max} < 50 — inner boundary tau will be < 50.")
        print(f"           DIFFUSION boundary (DIF=T) requires tau>=50 for LTE validity.")
        print(f"           With tau<50, FG_J_CMF_V13 may produce NEGATIVE mean intensities")
        print(f"           and halt immediately (observed on X_v2 with tau_max=30).")
        print(f"           For interaction models with steep tau jumps at shell-ejecta")
        print(f"           interfaces, use tau_max>=100 to push the inner boundary past")
        print(f"           the jump into the smooth inner shell.")
        print(f"           prepare_VADAT.py will auto-detect tau at inner boundary from SN_HYDRO_DATA")
        print(f"           and use ZERO_FLUX (DIF=F) if tau<50, or DIFFUSION (DIF=T) if tau>=50.")
    enable_grid_thinning = cmfgen_conv_cfg.get('enable_grid_thinning', False)
    enable_adaptive_grid = cmfgen_conv_cfg.get('enable_adaptive_grid', True)
    use_tau_amr = cmfgen_conv_cfg.get('use_tau_amr', False)
    smooth_physical_discontinuities = cmfgen_conv_cfg.get('smooth_physical_discontinuities', False)
    phys_smooth_sigma_zones = cmfgen_conv_cfg.get('phys_smooth_sigma_zones', 10.0)
    stretch_trouble_regions = cmfgen_conv_cfg.get('stretch_trouble_regions', False)
    stretch_factor = cmfgen_conv_cfg.get('stretch_factor', 3.0)
    max_stretch_zones = cmfgen_conv_cfg.get('max_stretch_zones', 200)
    min_zones_to_stretch = cmfgen_conv_cfg.get('min_zones_to_stretch', 2)
    # Critical-region detection and stretching (for sharp discontinuities like
    # reverse shocks where STELLA has very few zones spanning the transition).
    critical_regions = cmfgen_conv_cfg.get('critical_regions', False)
    stretch_critical_regions = cmfgen_conv_cfg.get('stretch_critical_regions', False)
    critical_stretch_factor = cmfgen_conv_cfg.get('critical_stretch_factor', 5.0)
    critical_max_stretch_zones = cmfgen_conv_cfg.get('critical_max_stretch_zones', 300)
    critical_dT_T = cmfgen_conv_cfg.get('critical_dT_T', 0.10)
    critical_drho_rho = cmfgen_conv_cfg.get('critical_drho_rho', 0.10)
    critical_dlogtau = cmfgen_conv_cfg.get('critical_dlogtau', 0.10)
    critical_dv_v = cmfgen_conv_cfg.get('critical_dv_v', 0.30)
    critical_pad_zones = cmfgen_conv_cfg.get('critical_pad_zones', 2)
    base_grid_target = cmfgen_conv_cfg.get('base_grid_target', 50)
    base_grid_method = cmfgen_conv_cfg.get('base_grid_method', 'tau_amr_recursive')
    tau_amr_config = cmfgen_conv_cfg.get('tau_amr_config', None)
    tau_breakpoints = cmfgen_conv_cfg.get('tau_breakpoints', [0.01, 1.0])
    tau_zone_fractions = cmfgen_conv_cfg.get('tau_zone_fractions', [0.10, 0.25, 0.65])
    tau_spacing = cmfgen_conv_cfg.get('tau_spacing', 'auto')
    # Composition-aware grid options
    composition_aware = cmfgen_conv_cfg.get('composition_aware', True)
    composition_gradient_threshold = cmfgen_conv_cfg.get('composition_gradient_threshold', 0.5)
    composition_elements = cmfgen_conv_cfg.get('composition_elements', None)
    min_zones_per_shell = cmfgen_conv_cfg.get('min_zones_per_shell', 3)
    composition_trace_threshold = cmfgen_conv_cfg.get('composition_trace_threshold', 1e-4)
    composition_grid_insertion = cmfgen_conv_cfg.get('composition_grid_insertion', True)
    composition_interface_cap = cmfgen_conv_cfg.get('composition_interface_cap', 0.15)
    max_thin_zones = cmfgen_conv_cfg.get('max_thin_zones', 0)
    max_thick_tau = cmfgen_conv_cfg.get('max_thick_tau', 1e10)
    min_points = cmfgen_conv_cfg.get('min_points', 150)
    max_radius_jump = cmfgen_conv_cfg.get('max_radius_jump', 0.05)
    max_temp_jump = cmfgen_conv_cfg.get('max_temp_jump', 0.10)
    max_density_jump = cmfgen_conv_cfg.get('max_density_jump', 0.15)
    max_ne_jump = cmfgen_conv_cfg.get('max_ne_jump', 0.10)
    max_logtau_jump = cmfgen_conv_cfg.get('max_logtau_jump', 0.10)
    max_velocity_jump = cmfgen_conv_cfg.get('max_velocity_jump', 0.10)
    max_sigma_jump = cmfgen_conv_cfg.get('max_sigma_jump', 1.0)
    max_cell_ratio = cmfgen_conv_cfg.get('max_cell_ratio', 1.5)
    max_points = cmfgen_conv_cfg.get('max_points', 150)
    min_nd = cmfgen_conv_cfg.get('min_nd', 150)
    thin_grid_tau = cmfgen_conv_cfg.get('thin_grid_tau', 1.0)
    thin_grid_priority = cmfgen_conv_cfg.get('thin_grid_priority', 1.0)
    outer_boundary_extra_points = cmfgen_conv_cfg.get('outer_boundary_extra_points', 5)
    outer_boundary_ratio = cmfgen_conv_cfg.get('outer_boundary_ratio', 1.5)
    min_dr_r = cmfgen_conv_cfg.get('min_dr_r', 1e-6)
    max_dr_r = cmfgen_conv_cfg.get('max_dr_r', 0.5)
    sigma_priority_threshold = cmfgen_conv_cfg.get('sigma_priority_threshold', 2.0)
    sigma_priority = cmfgen_conv_cfg.get('sigma_priority', 5.0)
    dsigma_priority_threshold = cmfgen_conv_cfg.get('dsigma_priority_threshold', 0.5)
    dsigma_priority = cmfgen_conv_cfg.get('dsigma_priority', 5.0)
    ne_priority_threshold = cmfgen_conv_cfg.get('ne_priority_threshold', 2.0)
    ne_priority = cmfgen_conv_cfg.get('ne_priority', 1.0)
    t_floor = cmfgen_conv_cfg.get('t_floor_K', 3000.0)
    t_cap = cmfgen_conv_cfg.get('t_cap_K', 0.0)
    trace_zero_threshold = cmfgen_conv_cfg.get('trace_zero_threshold', cmfgen_conv_cfg.get('trace_threshold', 0.0))
    velocity_max_kms = cmfgen_conv_cfg.get('velocity_max_kms', 0.0)
    v_min_km_s = cmfgen_conv_cfg.get('v_min_km_s', 0.0)
    enable_sigma_conformance = cmfgen_conv_cfg.get('enable_sigma_conformance', True)
    sigma_conformance_headroom = cmfgen_conv_cfg.get('sigma_conformance_headroom', 40)
    max_abs_sigma = cmfgen_conv_cfg.get('max_abs_sigma', 5.0)
    sigma_optimize = cmfgen_conv_cfg.get('sigma_optimize', True)

    ionization_front_width = cmfgen_conv_cfg.get('ionization_front_width', 3)
    ionization_front_cap = cmfgen_conv_cfg.get('ionization_front_cap', 0.15)

    # Recombination-front resolution (resolve under-resolved STELLA fronts)
    resolve_recombination_fronts = cmfgen_conv_cfg.get('resolve_recombination_fronts', True)
    recombination_t_jump = cmfgen_conv_cfg.get('recombination_t_jump', 0.3)
    recombination_front_width = cmfgen_conv_cfg.get('recombination_front_width', 5)

    # Velocity smoothing parameters
    smooth_velocity = cmfgen_conv_cfg.get('smooth_velocity', False)
    no_velocity_smoothing = cmfgen_conv_cfg.get('no_velocity_smoothing', True)  # Default to True for non-monotonic profiles
    enable_velocity_spacing_optimization = cmfgen_conv_cfg.get('enable_velocity_spacing_optimization', False)
    include_ne_amr = cmfgen_conv_cfg.get('include_ne_amr', False)  # Default to False - Ne can have artificial jumps
    ne_priority_threshold = cmfgen_conv_cfg.get('ne_priority_threshold', 2.0)
    ne_priority = cmfgen_conv_cfg.get('ne_priority', 1.0)
    sigma_threshold = cmfgen_conv_cfg.get('sigma_threshold', 2.0)
    smooth_width = cmfgen_conv_cfg.get('smooth_width', 5)
    smooth_sigma = cmfgen_conv_cfg.get('smooth_sigma', False)
    smooth_wiggles = cmfgen_conv_cfg.get('smooth_wiggles', False)
    wiggle_window = cmfgen_conv_cfg.get('wiggle_window', 3)
    wiggle_max_change = cmfgen_conv_cfg.get('wiggle_max_change', 0.02)
    max_dsigma = cmfgen_conv_cfg.get('max_dsigma', 1.0)
    sigma_smooth_passes = cmfgen_conv_cfg.get('sigma_smooth_passes', 10)
    sigma_smooth_window = cmfgen_conv_cfg.get('sigma_smooth_window', 3)
    velocity_smooth_window_kms = cmfgen_conv_cfg.get('velocity_smooth_window_kms', 0.0)
    smooth_sigma_scatter = cmfgen_conv_cfg.get('smooth_sigma_scatter', False)
    max_dsigma_scatter = cmfgen_conv_cfg.get('max_dsigma_scatter', 0.5)
    sigma_scatter_window = cmfgen_conv_cfg.get('sigma_scatter_window', 7)
    sigma_scatter_dv_threshold = cmfgen_conv_cfg.get('sigma_scatter_dv_threshold', 0.01)
    smooth_sigma_global = cmfgen_conv_cfg.get('smooth_sigma_global', False)
    sigma_global_window = cmfgen_conv_cfg.get('sigma_global_window', 7)
    sigma_global_max_rel_change = cmfgen_conv_cfg.get('sigma_global_max_rel_change', 0.05)
    sigma_global_preserve_threshold = cmfgen_conv_cfg.get('sigma_global_preserve_threshold', 2.0)
    # Final Savitzky-Golay smoothing (default OFF — creates sigma spikes at shock fronts)
    enable_final_smoothing = cmfgen_conv_cfg.get('enable_final_smoothing', False)

    print(f"\nConverting STELLA .swd to SN_HYDRO_DATA...")
    print(f"  .swd file: {swd_file}")
    print(f"  output: {sn_hydro}")
    print(f"  tau_max: {tau_max}")
    print(f"  tau_min: {tau_min}")
    
    cmd = [
        sys.executable, str(scripts_dir / converter_script),
        str(swd_file), str(sn_hydro),
        '--target-days', str(target_days),
        '--kappa', str(kappa),
        '--tau-max', str(tau_max),
        '--tau-min', str(tau_min),
        '--model-name', model_name,
        '--base-grid-target', str(base_grid_target),
        '--base-grid-method', str(base_grid_method),
    ]
    if tau_amr_config:
        cmd.extend(['--tau-amr-config', str(tau_amr_config)])
    if use_tau_amr:
        cmd.append('--use-tau-amr')
    cmd += [
        '--tau-breakpoints', ','.join(str(x) for x in tau_breakpoints),
        '--tau-zone-fractions', ','.join(str(x) for x in tau_zone_fractions),
        '--tau-spacing', str(tau_spacing),
        '--composition-gradient-threshold', str(composition_gradient_threshold),
        '--min-zones-per-shell', str(min_zones_per_shell),
        '--trace-threshold', str(composition_trace_threshold),
        '--composition-interface-cap', str(composition_interface_cap),
    ]
    if not composition_grid_insertion:
        cmd.append('--no-composition-grid-insertion')
    if not composition_aware:
        cmd.append('--no-composition-aware')
    if composition_elements:
        cmd += ['--composition-elements', ','.join(str(x) for x in composition_elements)]
    cmd += [
        '--max-thin-zones', str(max_thin_zones),
        '--max-thick-tau', str(max_thick_tau),
        '--min-points', str(min_points),
        '--max-radius-jump', str(max_radius_jump),
        '--max-temp-jump', str(max_temp_jump),
        '--max-density-jump', str(max_density_jump),
        '--max-ne-jump', str(max_ne_jump),
        '--max-logtau-jump', str(max_logtau_jump),
        '--max-velocity-jump', str(max_velocity_jump),
        '--max-sigma-jump', str(max_sigma_jump),
        '--max-cell-ratio', str(max_cell_ratio),
        '--max-points', str(max_points),
        '--min-nd', str(min_nd),
        '--thin-grid-tau', str(thin_grid_tau),
        '--thin-grid-priority', str(thin_grid_priority),
        '--outer-boundary-extra-points', str(outer_boundary_extra_points),
        '--outer-boundary-ratio', str(outer_boundary_ratio),
        '--min-dr-r', str(min_dr_r),
        '--max-dr-r', str(max_dr_r),
        '--sigma-priority-threshold', str(sigma_priority_threshold),
        '--sigma-priority', str(sigma_priority),
        '--dsigma-priority-threshold', str(dsigma_priority_threshold),
        '--dsigma-priority', str(dsigma_priority),
        '--ne-priority-threshold', str(ne_priority_threshold),
        '--ne-priority', str(ne_priority),
        '--t-floor', str(t_floor),
        '--t-cap', str(t_cap),
        '--trace-zero-threshold', str(trace_zero_threshold),
        '--velocity-max-kms', str(velocity_max_kms),
        '--v-min-km-s', str(v_min_km_s),
        '--sigma-conformance-headroom', str(sigma_conformance_headroom),
        '--max-abs-sigma', str(max_abs_sigma),
    ]
    cmd.append('--enable-sigma-conformance' if enable_sigma_conformance else '--no-sigma-conformance')
    if not sigma_optimize:
        cmd.append('--no-sigma-optimize')

    # Ionization-front-aware grid enhancement (v1 only — v2 doesn't support these)
    if converter_version != 'v2':
        cmd.append('--ionization-front-aware' if ionization_front_aware else '--no-ionization-front-aware')
        cmd += [
            '--ionization-front-threshold', str(ionization_front_threshold),
            '--ionization-front-min-points', str(ionization_front_min_points),
            '--ionization-front-width', str(ionization_front_width),
            '--ionization-front-cap', str(ionization_front_cap),
        ]

        # Recombination-front resolution
        cmd.append('--resolve-recombination-fronts' if resolve_recombination_fronts else '--no-resolve-recombination-fronts')
        cmd += [
            '--recombination-t-jump', str(recombination_t_jump),
            '--recombination-front-width', str(recombination_front_width),
        ]

    # Handle electron density AMR criterion
    if not include_ne_amr:
        # Override max_ne_jump with a huge value to effectively disable Ne-based refinement
        for i, arg in enumerate(cmd):
            if arg == '--max-ne-jump':
                cmd[i+1] = '1e10'
                break
        print("  Electron density AMR criterion DISABLED (max_ne_jump set to 1e10)")
    else:
        print("  Electron density AMR criterion ENABLED")
    if rho_file:
        cmd += ['--rho-file', str(rho_file)]
    elif abn_file:
        cmd += ['--abn-file', str(abn_file)]
    if input_profile:
        cmd += ['--input-profile', str(input_profile)]
    if xni_file:
        cmd += ['--xni-file', str(xni_file)]
        print(f"  Passing Ni56 distribution from .xni file: {xni_file}")
    if enable_grid_thinning:
        cmd.append('--enable-grid-thinning')
    if enable_adaptive_grid:
        cmd.append('--enable-adaptive-grid')
    if smooth_physical_discontinuities:
        cmd.append('--smooth-physical-discontinuities')
        cmd += ['--phys-smooth-sigma-zones', str(phys_smooth_sigma_zones)]
    if stretch_trouble_regions:
        cmd.append('--stretch-trouble-regions')
        cmd += ['--stretch-factor', str(stretch_factor)]
        cmd += ['--max-stretch-zones', str(max_stretch_zones)]
        cmd += ['--min-zones-to-stretch', str(min_zones_to_stretch)]
    
    # Critical-region detection and stretching
    if critical_regions:
        cmd.append('--critical-regions')
        cmd += ['--critical-dT-T', str(critical_dT_T)]
        cmd += ['--critical-drho-rho', str(critical_drho_rho)]
        cmd += ['--critical-dlogtau', str(critical_dlogtau)]
        cmd += ['--critical-dv-v', str(critical_dv_v)]
        cmd += ['--critical-dne-ne', str(critical_dne_ne)]
        cmd += ['--critical-pad-zones', str(critical_pad_zones)]
        if stretch_critical_regions:
            cmd.append('--stretch-critical-regions')
            cmd += ['--critical-stretch-factor', str(critical_stretch_factor)]
            cmd += ['--critical-max-stretch-zones', str(critical_max_stretch_zones)]
            print(f"  Critical region stretching ENABLED (factor={critical_stretch_factor}, max_zones={critical_max_stretch_zones})")
        print(f"  Critical region detection ENABLED (dT/T={critical_dT_T}, drho/rho={critical_drho_rho}, pad={critical_pad_zones})")

    # Handle velocity smoothing based on config
    # Note: velocity smoothing is OFF by default in convert_stella_swd_to_cmfgen.py
    # (opt-in via --smooth-velocity). We only pass the flag when explicitly enabled.
    if smooth_velocity and not no_velocity_smoothing:
        cmd.append('--smooth-velocity')
        cmd += ['--sigma-threshold', str(sigma_threshold)]
        cmd += ['--smooth-width', str(smooth_width)]
        print("  Velocity smoothing ENABLED (smooth_velocity=True)")
    else:
        print("  Velocity smoothing DISABLED - preserving non-monotonic velocity (Dessart+2015)")

    # Sigma-transition smoothing: smooth the sigma gradient (not magnitude)
    if smooth_sigma:
        cmd.append('--smooth-sigma')
        cmd += ['--max-dsigma', str(max_dsigma)]
        cmd += ['--sigma-smooth-passes', str(sigma_smooth_passes)]
        cmd += ['--sigma-smooth-window', str(sigma_smooth_window)]
        print(f"  Sigma-transition smoothing ENABLED (max_dsigma={max_dsigma}, passes={sigma_smooth_passes})")
    else:
        cmd.append('--no-smooth-sigma')
        print("  Sigma-transition smoothing DISABLED")
    if not smooth_wiggles:
        cmd.append('--no-smooth-wiggles')
        print("  Velocity wiggle smoothing DISABLED - preserving the input velocity profile")
    else:
        cmd += ['--wiggle-window', str(wiggle_window)]
        cmd += ['--wiggle-max-change', str(wiggle_max_change)]
        print(f"  Velocity wiggle smoothing ENABLED (window={wiggle_window} pts, max_change={wiggle_max_change*100:.0f}%)")
    if velocity_smooth_window_kms > 0:
        cmd += ['--velocity-smooth-window-kms', str(velocity_smooth_window_kms)]
        print(f"  Velocity smoothing window: {velocity_smooth_window_kms} km/s")
    if smooth_sigma_scatter:
        cmd.append('--smooth-sigma-scatter')
        cmd += ['--max-dsigma-scatter', str(max_dsigma_scatter)]
        cmd += ['--sigma-scatter-window', str(sigma_scatter_window)]
        cmd += ['--sigma-scatter-dv-threshold', str(sigma_scatter_dv_threshold)]
        print(f"  Sigma scatter smoothing ENABLED (max_dsigma={max_dsigma_scatter}, window={sigma_scatter_window} pts)")
    else:
        cmd.append('--no-smooth-sigma-scatter')
        print("  Sigma scatter smoothing DISABLED")
    if smooth_sigma_global:
        cmd.append('--smooth-sigma-global')
        cmd += ['--sigma-global-window', str(sigma_global_window)]
        cmd += ['--sigma-global-max-rel-change', str(sigma_global_max_rel_change)]
        cmd += ['--sigma-global-preserve-threshold', str(sigma_global_preserve_threshold)]
        print(f"  Global sigma smoothing ENABLED (window={sigma_global_window}, max_change={sigma_global_max_rel_change*100:.0f}%, preserve>{sigma_global_preserve_threshold})")
    else:
        cmd.append('--no-smooth-sigma-global')
        print("  Global sigma smoothing DISABLED")

    # Final Savitzky-Golay smoothing (default OFF)
    if enable_final_smoothing:
        cmd.append('--enable-final-smoothing')
        print("  Final S-G smoothing ENABLED (may create sigma spikes at shock fronts)")
    else:
        cmd.append('--no-final-smoothing')
        print("  Final S-G smoothing DISABLED — preserving raw velocities at shock fronts")

    # Sigma-aware negative-sigma coarsening (Luc Dessart approach)
    # v2 converter doesn't support this option
    if converter_version != 'v2':
        if coarsen_negative_sigma:
            cmd.append('--coarsen-negative-sigma')
            cmd += ['--coarsen-sigma-threshold', str(coarsen_sigma_threshold)]
            cmd += ['--coarsen-max-zones', str(coarsen_max_zones)]
            cmd += ['--coarsen-min-region-zones', str(coarsen_min_region_zones)]
            print(f"  Negative-sigma coarsening ENABLED (threshold={coarsen_sigma_threshold}, max_zones={coarsen_max_zones})")
        else:
            cmd.append('--no-coarsen-negative-sigma')
            print("  Negative-sigma coarsening DISABLED — preserving full non-monotonic velocity structure")
    else:
        print(f"  Negative-sigma coarsening: N/A (v2 converter)")

    # INC_RAD_DECAYS: strip isotope blocks from SN_HYDRO_DATA when disabled
    if not inc_rad_decays:
        cmd.append('--no-inc-rad-decays')
        print("  Radioactive decay isotopes DISABLED — SN_HYDRO_DATA will have NUM_ISOTOPES=0")

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"ERROR: STELLA conversion failed")
        print(result.stdout)
        print(result.stderr)
        return False
    
    print(f"  SN_HYDRO_DATA created: {sn_hydro}")

    # Generate GREY_SCL_FAC_IN from SN_HYDRO_DATA for better initial T estimates
    generate_grey_scl_fac_in(sn_hydro)

    # Post-conversion: check non-monotonic velocity severity and optionally enforce homologous expansion
    enforce_homologous = cmfgen_conv_cfg.get('enforce_homologous', False)
    sigma_warn_fraction = cmfgen_conv_cfg.get('sigma_warn_fraction', 0.20)
    check_sn_hydro_sigma(sn_hydro, enforce_homologous=enforce_homologous,
                         sigma_warn_fraction=sigma_warn_fraction)

    return True


def check_sn_hydro_sigma(sn_hydro_path, enforce_homologous=False, sigma_warn_fraction=0.20):
    """Post-conversion diagnostic: check non-monotonic velocity severity in SN_HYDRO_DATA.

    Reads the SN_HYDRO_DATA file, computes sigma = dln(v)/dln(r) - 1 at each zone
    interface, and reports the fraction of zones with sigma < -1 (non-monotonic).

    If enforce_homologous=True, overrides the velocity to V = R / t_sn using the
    SN age from the file header, following the approach of Dessart+2022 (A&A 660,
    L22, arXiv:2204.00446) for interacting SNe with optically thin CSM. This
    eliminates all non-monotonic velocity issues by enforcing homologous expansion,
    allowing the standard CMFGEN solver (rather than the Sobolev nonmonotonic
    solver) to be used.

    References:
      - Dessart, Audit & Hillier 2015, MNRAS 449, 4304 (arXiv:1503.05463):
        Sobolev nonmonotonic solver, 100 grid points, tau_Ross truncation at 30.
      - Dessart & Hillier 2022, A&A 660, L22 (arXiv:2204.00446):
        Homologous expansion enforced for optically thin CSM interaction,
        shock power deposited as additional energy source in CMFGEN.
      - Dessart, Hillier & Audit 2016, MNRAS 458, 3543 (arXiv:1602.02977):
        Nonmonotonic solver with arbitrary velocity fields, Sobolev for lines.
      - Dessart, Hillier & Audit 2017, A&A 605, A83 (arXiv:1704.01697):
        RSG explosions in CSM, nonmonotonic solver with HERACLES snapshots.

    Args:
        sn_hydro_path: Path to SN_HYDRO_DATA file.
        enforce_homologous: If True, rewrite velocity as V = R / t_sn.
        sigma_warn_fraction: Warn if fraction of zones with sigma < -1 exceeds this.
    """
    import numpy as np

    if not sn_hydro_path.exists():
        print(f"  WARNING: Cannot check sigma — SN_HYDRO_DATA not found: {sn_hydro_path}")
        return

    try:
        with open(sn_hydro_path) as f:
            lines = f.readlines()

        # Parse header
        nd = None
        t_sn_days = None
        for line in lines[:10]:
            if 'Number of data points' in line:
                nd = int(line.split(':')[1].strip())
            if 'Time(days)' in line:
                t_sn_days = float(line.split(':')[1].strip())

        if nd is None or t_sn_days is None:
            print(f"  WARNING: Could not parse SN_HYDRO_DATA header for sigma check")
            return

        # Read sections
        sections = {}
        for i, line in enumerate(lines):
            for kw in ['Radius grid', 'Velocity', 'Density', 'Temperature']:
                if kw in line:
                    sections[kw] = i
                    break

        def read_nd_vals(lines, start, nd):
            vals = []
            for i in range(start, len(lines)):
                parts = lines[i].split()
                for p in parts:
                    try:
                        vals.append(float(p))
                    except ValueError:
                        pass
                if len(vals) >= nd:
                    break
            return np.array(vals[:nd])

        r = read_nd_vals(lines, sections['Radius grid'] + 1, nd)
        v = read_nd_vals(lines, sections['Velocity'] + 1, nd)

        if len(r) != nd or len(v) != nd:
            print(f"  WARNING: Could not parse R/V arrays for sigma check (r={len(r)}, v={len(v)}, nd={nd})")
            return

        # Compute sigma = dln(v)/dln(r) - 1
        dln_v = np.diff(np.log(v))
        dln_r = np.diff(np.log(r))
        sigma = dln_v / dln_r - 1

        n_neg = int(np.sum(sigma < -1.0))
        n_total = len(sigma)
        frac = n_neg / n_total if n_total > 0 else 0

        # Find the most extreme sigma
        min_sigma = float(np.min(sigma)) if n_total > 0 else 0.0
        max_dv = float(np.max(np.abs(np.diff(v))))

        print(f"  Non-monotonic velocity diagnostic:")
        print(f"    ND={nd}, sigma < -1 zones: {n_neg}/{n_total} ({frac*100:.1f}%)")
        print(f"    Min sigma: {min_sigma:.2f}, Max |dV|: {max_dv:.1f} km/s")
        print(f"    V range: {v.min():.1f} - {v.max():.1f} km/s")

        if frac > sigma_warn_fraction:
            print(f"  WARNING: {frac*100:.1f}% of zones have sigma < -1 (threshold: {sigma_warn_fraction*100:.0f}%)")
            print(f"    The Sobolev nonmonotonic solver (USE_LAM_ES=T, GLOBAL_LINE=SOB) may struggle")
            print(f"    to converge with this many non-monotonic zones. Consider:")
            print(f"    1. Setting 'enforce_homologous: true' in cmfgen_conversion config")
            print(f"       (Dessart+2022 approach: V=R/t, eliminates non-monotonic velocity)")
            print(f"    2. Reducing grid resolution in the non-monotonic region")
            print(f"    3. Using fewer zones (Dessart+2015 uses 100 grid points)")

        if enforce_homologous:
            # Check if V=R/t would produce super-relativistic velocities
            # (happens for interacting SNe where outer ejecta is much older than t_sn)
            t_sn_s = t_sn_days * 86400.0  # seconds
            v_homologous = r * 1e10 / t_sn_s / 1e5  # R in cm -> km/s
            v_max_homologous = float(np.max(v_homologous))

            if v_max_homologous > 1.0e5:
                # Super-relativistic: use cumulative maximum monotonization instead
                # This preserves the original velocity range while removing non-monotonic decreases.
                # Physically motivated: in interacting SNe, the outer SN_A ejecta has been expanding
                # for much longer than t_sn (the SN_B interaction time), so V=R/t_sn is invalid.
                # Cumulative maximum ensures V increases monotonically with R (sigma >= 0).
                # SN_HYDRO_DATA stores values outer-to-inner (R decreasing), so we reverse,
                # apply cumulative maximum (inner-to-outer), then reverse back.
                print(f"  Monotonizing velocity (cumulative maximum) for interacting SN")
                print(f"    V=R/t_sn would give V_max={v_max_homologous:.1f} km/s (super-relativistic)")
                print(f"    Using cumulative maximum instead to preserve original velocity range")
                print(f"    Original V range: {v.min():.1f} - {v.max():.1f} km/s")
                v_reversed = v[::-1]
                v_mono = np.maximum.accumulate(v_reversed)
                v_new = v_mono[::-1]
                print(f"    Monotonized V range: {v_new.min():.1f} - {v_new.max():.1f} km/s")
            else:
                print(f"  Enforcing homologous expansion (V = R / t_sn) per Dessart+2022 (arXiv:2204.00446)")
                print(f"    t_sn = {t_sn_days:.4f} days = {t_sn_days * 86400:.2e} s")
                print(f"    This overrides the non-monotonic velocity with homologous expansion.")
                print(f"    Use with VEL_LAW=12 and USE_LAM_ES=T in VADAT (sigma will be ~0 everywhere).")
                v_new = v_homologous

            # Rewrite the velocity section in SN_HYDRO_DATA
            v_section_start = sections['Velocity'] + 1
            v_section_end = v_section_start
            v_count = 0
            while v_count < nd and v_section_end < len(lines):
                parts = lines[v_section_end].split()
                v_count += len(parts)
                v_section_end += 1

            # Format velocity values (6 per line, Fortran scientific notation)
            new_v_lines = []
            for i in range(0, nd, 6):
                chunk = v_new[i:i+6]
                line_str = '     ' + '     '.join(f'{x:16.13E}' for x in chunk)
                new_v_lines.append(line_str)

            # Replace velocity section
            lines = lines[:v_section_start] + new_v_lines + lines[v_section_end:]

            with open(sn_hydro_path, 'w') as f:
                f.writelines(lines)

            # Verify
            print(f"    Velocity overridden: V range {v_new.min():.1f} - {v_new.max():.1f} km/s")
            dln_v_new = np.diff(np.log(v_new))
            sigma_new = dln_v_new / dln_r - 1
            n_neg_new = int(np.sum(sigma_new < -1.0))
            print(f"    After override: sigma < -1 zones: {n_neg_new}/{n_total} (was {n_neg})")

    except Exception as e:
        print(f"  WARNING: Sigma check failed: {e}")


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
        # First column is time in days (was incorrectly using last column gdepos)
        final_time_days = data[-1, 0]
        target_days = config.get('stella', {}).get('target_days',
                         config.get('heracles', {}).get('target_days', 40))  # DEPRECATED: HERACLES fallback
        # Consider complete if within 1% of target or reached target
        return final_time_days >= target_days * 0.99
    except Exception:
        return False


def generate_initial_plots(config, scripts_dir, project_dir):
    """Generate initial diagnostic plots (structure + composition) from .hyd/.abn files.

    Calls plot_model.py --step stella_input to produce:
    - Structure plot (r, rho, v, T)
    - Composition plot (mass fractions vs radius)
    - Mass fraction distribution (mass fractions vs enclosed mass)
    """
    model_name = config.get('model_name', 'unknown')
    output = config.get('output', {})
    model_dir_str = output.get('model_dir', '../models/${model_name}').replace('${model_name}', model_name)
    model_dir = (project_dir / model_dir_str).absolute()

    plot_script = scripts_dir / 'plot_model.py'
    if not plot_script.exists():
        print(f"WARNING: plot_model.py not found, skipping initial plots")
        return False

    # Check that .hyd file exists
    hyd_file = model_dir / 'stella' / 'modmake' / f'{model_name}.hyd'
    if not hyd_file.exists():
        print(f"WARNING: .hyd file not found ({hyd_file}), skipping initial plots")
        return False

    print(f"\n{'='*60}")
    print(f"Generating initial diagnostic plots for {model_name}")
    print(f"{'='*60}")

    plots_dir = model_dir / 'plots'
    plots_dir.mkdir(parents=True, exist_ok=True)

    result = subprocess.run(
        [sys.executable, str(plot_script), str(model_dir), '--step', 'stella_input'],
        capture_output=True, text=True
    )

    if result.returncode == 0:
        print(result.stdout)
        print(f"Initial plots saved to: {plots_dir}")
        return True
    else:
        print(f"WARNING: Initial plot generation failed")
        print(result.stderr)
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


def plot_stella_output_plots(config, scripts_dir, project_dir):
    """Generate stella_output plots (structure, composition, lightcurve) from .swd."""
    model_name = config.get('model_name', 'unknown')
    output = config.get('output', {})
    model_dir_str = output.get('model_dir', '../models/${model_name}').replace('${model_name}', model_name)
    model_dir = (project_dir / model_dir_str).absolute()

    plot_script = scripts_dir / 'plot_model.py'
    if not plot_script.exists():
        print(f"WARNING: plot_model.py not found, skipping stella_output plots")
        return False

    # Check that .swd file exists
    swd_file = model_dir / 'stella' / 'strad' / f'{model_name}.swd'
    if not swd_file.exists():
        swd_file = model_dir / 'stella' / f'{model_name}.swd'
    if not swd_file.exists():
        print(f"WARNING: .swd file not found, skipping stella_output plots")
        return False

    print(f"\n  Generating STELLA output plots for {model_name}")
    result = subprocess.run(
        [sys.executable, str(plot_script), str(model_dir), '--step', 'stella_output'],
        capture_output=True, text=True
    )
    if result.returncode == 0:
        print(result.stdout)
        return True
    else:
        print(f"WARNING: stella_output plot generation failed")
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
                     config.get('heracles', {}).get('target_days', 40))  # DEPRECATED: HERACLES fallback
    
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


def extract_stella_lightcurve_pysynphot(config, scripts_dir, project_dir):
    """Run pysynphot lightcurve extraction after STELLA completes.

    Produces a gap-free multi-band lightcurve (UBVRI + SDSS ugriz) from
    STELLA's .flx SED and .swd snapshots. Output is written to
    <model_dir>/output/lightcurve.csv in the same format as the existing
    lightcurve.csv with ugriz bands appended.
    """
    model_name = config.get('model_name', 'unknown')
    output = config.get('output', {})
    model_dir = output.get('model_dir', '../models/${model_name}').replace('${model_name}', model_name)
    model_dir = (project_dir / model_dir).absolute()

    # Locate input files
    swd_file = model_dir / 'stella' / 'strad' / f'{model_name}.swd'
    if not swd_file.exists():
        swd_file = model_dir / 'stella' / f'{model_name}.swd'

    flx_file = model_dir / 'stella' / 'strad' / f'{model_name}.flx'
    if not flx_file.exists():
        flx_file = model_dir / 'stella' / f'{model_name}.flx'

    tt_file = model_dir / 'stella' / f'{model_name}.tt'

    if not swd_file.exists() or not flx_file.exists():
        print(f"  [pysynphot] Missing .swd or .flx for {model_name}, skipping lightcurve extraction")
        return False

    out_dir = model_dir / 'output'
    out_dir.mkdir(parents=True, exist_ok=True)
    out_csv = out_dir / 'lightcurve.csv'
    plot_png = model_dir / 'plots' / 'lightcurve_pysynphot.png'
    plot_png.parent.mkdir(parents=True, exist_ok=True)

    extract_script = scripts_dir / 'extract_lightcurve_pysynphot.py'
    if not extract_script.exists():
        print(f"WARNING: extract_lightcurve_pysynphot.py not found, skipping")
        return False

    print(f"\n{'='*60}")
    print(f"Extracting pysynphot lightcurve for {model_name}")
    print(f"{'='*60}")

    cmd = [
        sys.executable, str(extract_script),
        str(swd_file), str(flx_file),
        '--tt-file', str(tt_file),
        '--output', str(out_csv),
        '--plot', str(plot_png),
        '--title', f'{model_name} Lightcurve (pysynphot)',
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode == 0:
        print(f"Lightcurve CSV saved: {out_csv}")
        return True
    else:
        print(f"WARNING: pysynphot lightcurve extraction failed")
        print(result.stderr)
        return False


def enforce_cmfgen_sn_iin_settings(config):
    """
    [CRITICAL] Mutates the config dictionary to enforce the correct CMFGEN VADAT
    settings for STELLA snapshots of interacting transients (Type IIn SNe).

    Settings are divided into two categories:

    1. **Non-overridable defaults** — fundamental physics settings that must
       always be set for STELLA snapshots.  These are applied unconditionally
       and cannot be overridden via YAML ``vadat_overrides``.

    2. **User-overridable defaults** — settings that have a sensible default
       but may legitimately need to change for specific models (e.g. high-tau
       inner boundary requires ``DIF=T`` / ``IB_METH=DIFFUSION``; non-monotonic
       velocity requires ``NON_MONO_SOLVER=T``).  These are only set if the
       user has NOT already specified them in the YAML ``vadat_overrides``.

    This mirrors the ``user_overridable`` mechanism in
    ``prepare_VADAT.py:enforce_sn_iin_physics()`` so that the two enforcers
    agree on which keys the user is allowed to control.
    """
    if 'cmfgen' not in config:
        config['cmfgen'] = {}

    cmfgen = config['cmfgen']

    # 1. Execution & Model Type (non-overridable)
    cmfgen['model_type'] = 'static_snapshot'
    # VEL_LAW: respect config value for non-monotonic solver selection.
    # VEL_LAW=12 (RV_SN_MODEL_SNIIN) is designed for non-monotonic velocity
    # flows (interacting/IIn SNe). VEL_LAW=11 (SET_RV_HYDRO_MODEL_V3) reads
    # R/V directly from RDINR. Default to 12 for non-monotonic solver.
    if 'vel_law' not in cmfgen:
        cmfgen['vel_law'] = 12
    cmfgen['incl_rel'] = False

    # 2. VADAT Overrides
    if 'vadat_overrides' not in cmfgen:
        cmfgen['vadat_overrides'] = {}
    v = cmfgen['vadat_overrides']

    # ── Non-monotonic solver detection ────────────────────────────────
    # The non-monotonic solver (Knop et al. 2009, FG_J_CMF_NONMONO_V1)
    # is activated by NON_MONO_SOLVER=T in VADAT.  It is also implied by
    # VEL_LAW=11 or VEL_LAW=12 (SN hydrodynamic structure readers that
    # preserve non-monotonic velocity from STELLA).
    #
    # When the non-monotonic solver is active, the following settings are
    # recommended defaults.  They are NOT forced — the user may override
    # them via YAML vadat_overrides if they have a specific reason:
    #
    #   USE_LAM_ES=T   — moment equations (MOM_J_CMF_V11) are unstable
    #                    for sigma < -1; bypassing them is recommended.
    #                    The user may set USE_LAM_ES=F to enable the dual-
    #                    sigma architecture (clamped SIGMA for moments,
    #                    SIGMA_TRUE for formal solver).
    #   N_TYPE=N_ON_J  — G_ONLY (N/H) is singular at velocity reversals
    #                    where H=0; N/J is always well-defined (J>0).
    #   SOL_METH       — TRIDIAG (NUM_BNDS=3) is the CMFGEN-recommended
    #                    setting (manual: "best compromise between
    #                    convergence and memory requirements").  DIAG
    #                    (NUM_BNDS=1) is slower but may be more robust
    #                    for strongly non-monotonic flows.
    #
    # KNOWN ISSUE (BA_AUDIT_REPORT.txt): The BA variation operator
    # (VAR_MOM_J_CMF_V12) uses clamped SIGMA but receives Eddington
    # factors from the non-monotonic formal solve (GET_MOMS_NON_REL).
    # This inconsistency can produce zero BA corrections (NaN silently
    # skipped by comparison logic) in strongly reversing flows,
    # regardless of whether DIAG or TRIDIAG is used.  The proper fix
    # (audit option a) is to use consistent Eddington factors for the
    # BA variation, not to switch from TRIDIAG to DIAG.
    non_mono_active = (
        str(v.get('NON_MONO_SOLVER', '')).strip().upper() == 'T'
        or cmfgen.get('vel_law', 12) in (11, 12)
    )

    if non_mono_active:
        # Warn if USE_LAM_ES=F (moment equations unstable for sigma < -1)
        if str(v.get('USE_LAM_ES', '')).strip().upper() == 'F':
            print(f"  [NON_MONO] WARNING: USE_LAM_ES=F with non-monotonic solver. "
                  f"Moment equations (MOM_J_CMF_V11) are unstable for sigma < -1. "
                  f"Recommend USE_LAM_ES=T to bypass moment solver.")
        # Warn if N_TYPE=G_ONLY (singular at velocity reversals)
        if str(v.get('N_TYPE', '')).strip().upper() == 'G_ONLY':
            print(f"  [NON_MONO] WARNING: N_TYPE=G_ONLY with non-monotonic solver. "
                  f"G_ONLY (N/H) is singular at velocity reversals where H=0. "
                  f"Recommend N_TYPE=N_ON_J (N/J is well-defined since J>0).")

    # Keys that the user may override via YAML vadat_overrides.
    # If already present in v, the user's value is preserved.
    # This list must stay consistent with prepare_VADAT.py's
    # enforce_sn_iin_physics() user_overridable set.
    user_overridable = {
        'DIF', 'IB_METH', 'IB_STAB', 'MAX_LAM', 'MAX_LIN',
        'BEG_NG', 'IBEG_NG', 'BW_NG', 'N_OB_INS', 'N_IB_INS',
        'INC_GRID', 'ALL_FREQ', 'SCL_LN', 'DO_AV', 'ITS/NG',
        'BA_CHK_FAC', 'T_INIT_TAU', 'MIN_CF', 'LAM_VAL', 'NUM_LAM',
        'N_FIX_BA', 'NEW_LINE_BA', 'IND_BA_METH',
        'ST_INT', 'END_INT', 'IT_ON_T', 'T_ON_R', 'SCL_SL_OPAC',
        'NON_MONO_SOLVER', 'FIX_NE', 'FIX_IMP', 'DO_NG',
        'DO_LAM_IT', 'DO_LAM_AUTO', 'USE_LAM_ES',
        'OB_RAT', 'THK_CONT', 'STOP_IF_BP', 'REV_RGRID',
        'LTE_EST', 'FIX_T', 'FIX_T_AUTO',
        'NEG_OPAC_OPT', 'DC_METH',
    }

    # Helper: set a default only if the user hasn't already specified it.
    def _default(key, val):
        if key not in v:
            v[key] = val

    # Grid & Interpolation (Handles non-monotonic T and shocks)
    _default('T_ON_R', 'T')
    _default('DC_METH', 'TR')
    _default('INC_GRID', 'F')
    _default('ALL_FREQ', 'F')

    # Inner Boundary
    # Default: ZERO_FLUX / DIF=F (optically thin core).
    # User may override to DIF=T / IB_METH=DIFFUSION for high-tau (>=50)
    # inner boundaries where the diffusion boundary condition is physically
    # appropriate (required by PERTJFEAU_IBC for Eddington-factor variation
    # matrices).
    _default('IB_METH', 'ZERO_FLUX')
    _default('DIF', 'F')
    _default('USE_LAM_ES', 'T')
    v['N_TYPE'] = 'N_ON_J'       # non-overridable
    # SOL_METH: controlled by depth_independent_solver config flag.
    # depth_independent_solver=True  -> DIAG (NUM_BNDS=1, each depth independent)
    # depth_independent_solver=False -> TRIDIAG (NUM_BNDS=3, couples adjacent depths)
    # CMFGEN manual: "TRIDIAG is the best compromise between convergence and
    # memory requirements. DIAG ... has a slower rate of convergence."
    # TRIDIAG is the default recommendation for both monotonic and non-monotonic
    # models.  DIAG may be used as a fallback if BA convergence issues arise.
    # KNOWN ISSUE: VAR_MOM_J_CMF_V12 uses clamped SIGMA with non-monotonic
    # Eddington factors — this can produce zero BA corrections in strongly
    # reversing flows regardless of DIAG/TRIDIAG (see BA_AUDIT_REPORT.txt).
    depth_indep = cmfgen.get('depth_independent_solver', True)
    if depth_indep:
        v['SOL_METH'] = 'DIAG'       # non-overridable when depth_independent
    else:
        v['SOL_METH'] = 'TRI'        # non-overridable when TRIDIAG (couples depths)
    v['REL_OBS'] = 'F'           # non-overridable
    v['REL_CMF'] = 'F'           # non-overridable

    # AUTO_SMOOTH: forced off — population smoothing destabilizes non-monotonic
    # velocity models (causes divergent SMOOTH_POP_CORRECTIONS). Not overridable.
    v['AUTO_SMOOTH'] = 'F'

    # Static Snapshot (No time derivatives) — non-overridable
    v['DO_DDT'] = 'F'
    v['INC_ADV'] = 'F'
    v['INC_AD'] = 'F'
    v['INCL_DJDT'] = 'F'
    v['USE_DJDT_RTE'] = 'F'

    # Line Formation (Sobolev)
    v['GLOBAL_LINE'] = 'SOB'     # non-overridable
    v['CHK_L_POS'] = 'T'         # non-overridable
    _default('NEG_OPAC_OPT', 'SRCE_CHK')
    _default('SCL_SL_OPAC', 'T')

    # Initialization & Convergence
    v['LIN_INT'] = 'F'           # non-overridable
    _default('LTE_EST', 'T')
    _default('FIX_T', 'T')         # Static snapshot: T fixed to hydro
    _default('FIX_T_AUTO', 'T')    # Signals DO_LAM_AUTO=T (exit Lambda for BA)
    _default('FIX_NE', 'T')        # Static snapshot: n_e fixed to hydro
    _default('FIX_IMP', 'T')       # Static snapshot: impurity pops fixed
    _default('DO_LAM_IT', 'T')
    _default('DO_LAM_AUTO', 'T')   # Auto-switch to BA when MAXCH<50%
    _default('MAX_LAM', '1.0D5')   # Manual default: 1e5 (Knop solver is stable)
    _default('NUM_LAM', '10')
    _default('DO_NG', 'T')
    _default('BEG_NG', '1.0D-01')
    _default('IBEG_NG', '50')
    _default('BW_NG', '60')

    # 3. SEDONA Overrides (If SEDONA is used alongside CMFGEN)
    if 'sedona' in config:
        sed = config['sedona']
        sed['use_homologous_expansion'] = False
        sed['opacity_line_expansion'] = 1
        sed['opacity_bound_bound'] = 0

    return config


def process_config(config_file, scripts_dir, project_dir, j_release=False):
    """Process a single configuration file.

    All parameters are read from the YAML config file — no CLI overrides.
    """
    config = load_yaml_config(config_file)
    config = enforce_cmfgen_sn_iin_settings(config)
    model_name = config.get('model_name', 'unknown')

    # Read all workflow parameters from config
    debug = config.get('debug', False)
    debug_days = config.get('debug_days', 1)
    skip_stella = config.get('skip_stella', False)
    cancel_jobs = config.get('cancel_jobs', True)
    force_cancel = config.get('force_cancel', True)
    plot_lightcurve = config.get('plot_lightcurve', True)
    redo_heracles = config.get('redo_heracles', True)
    stella_dependency = os.environ.get('RUN_MODEL_GRID_STELLA_DEP', None)

    # Spectrum code flags from config
    run_cmfgen = config.get('run_cmfgen', True)
    run_sedona = config.get('run_sedona', True)
    run_sirocco = config.get('run_sirocco', True)
    run_superlite = config.get('run_superlite', False)
    
    # Read reuse_stella flag from config
    reuse_stella = config.get('reuse_stella', False)
    # stella_source_model is deprecated - no longer needed
    
    # Determine which codes are being used
    # DEPRECATED: default 'heracles' backend is deprecated; use 'stella' (STELLA -> CMFGEN pipeline)
    use_stella = config.get('hydro_backend', 'heracles').lower() == 'stella'
    hydro_code = 'STELLA' if use_stella else '[DEPRECATED] HERACLES'
    
    print(f"\n{'='*70}")
    print(f"Processing model: {model_name}")
    print(f"Config file: {config_file}")
    print(f"Hydro code: {hydro_code} (hydrodynamic radiation modeling)")
    active_codes = [c for c, flag in [('CMFGEN', run_cmfgen), ('SEDONA', run_sedona), ('SIROCCO', run_sirocco), ('SuperLite', run_superlite)] if flag]
    if active_codes:
        print(f"Spectrum codes: {' + '.join(active_codes)} (independent spectrum generation)")
    else:
        print(f"Spectrum code: None (all spectrum codes disabled)")
    if debug:
        print(f"DEBUG MODE ENABLED")
    print(f"{'='*70}")
    
    # Check for existing jobs (cancel by default unless --keep-jobs is specified)
    # When only SEDONA is being rerun (run_cmfgen=False), only cancel SEDONA jobs
    # to preserve running CMFGEN/STELLA jobs. The SEDONA-specific cleanup is also
    # handled by clean_previous_sedona_runs() later in the workflow.
    if cancel_jobs:
        if run_cmfgen:
            # Rerunning CMFGEN too — cancel all jobs for this model
            cancelled = cancel_existing_jobs(model_name, force=force_cancel, verbose=True)
            if cancelled == 0:
                print(f"No existing jobs found for model {model_name}")
        elif run_sedona:
            # Only rerunning SEDONA — cancel only SEDONA jobs, preserve CMFGEN/STELLA
            print(f"\nSEDONA-only rerun: cancelling only SEDONA jobs (preserving CMFGEN/STELLA)")
            cancelled = cancel_existing_jobs(model_name, force=force_cancel, verbose=True, job_filter='sedona')
            if cancelled == 0:
                print(f"No existing SEDONA jobs found for model {model_name}")
        elif run_superlite:
            # Only rerunning SuperLite — cancel only SuperLite jobs, preserve CMFGEN/STELLA
            print(f"\nSuperLite-only rerun: cancelling only SuperLite jobs (preserving CMFGEN/STELLA)")
            cancelled = cancel_existing_jobs(model_name, force=force_cancel, verbose=True, job_filter='superlite')
            if cancelled == 0:
                print(f"No existing SuperLite jobs found for model {model_name}")
        else:
            print(f"\nNo spectrum codes enabled — skipping job cancellation")
    else:
        print(f"\n--keep-jobs specified: skipping job cancellation")
    
    # Extract paths
    output = config.get('output', {})
    
    # If debug mode, override simulation days and target days
    debug_config_file = None
    if debug:
        config = config.copy()
        config['model_name'] = f"{model_name}-debug"
        config['heracles'] = config.get('heracles', {}).copy()  # DEPRECATED: HERACLES backend, use STELLA instead
        config['heracles']['simulation_days'] = debug_days
        config['heracles']['target_days'] = debug_days  # Target same as simulation for debug
        # Also override STELLA target_days
        config['stella'] = config.get('stella', {}).copy()
        config['stella']['target_days'] = debug_days
        # Override snapshot_days to a single epoch at debug_days so the time-series
        # workflow produces just one CMFGEN + one SEDONA job (not the full list).
        config['cmfgen_conversion'] = config.get('cmfgen_conversion', {}).copy()
        config['cmfgen_conversion']['snapshot_days'] = [float(debug_days)]
        # single_final is always True now — no INIT/STABLE stages to skip.
        config['cmfgen'] = config.get('cmfgen', {}).copy()
        # Update model_dir to include -debug suffix
        output['model_dir'] = output.get('model_dir', '../models/${model_name}').replace('${model_name}', f"{model_name}-debug")
        # Also update stella_output path to include -debug suffix
        output['stella_output'] = output.get('stella_output', 'stella').replace('${model_name}', f"{model_name}-debug")
        # Also update plots_dir path to include -debug suffix
        output['plots_dir'] = output.get('plots_dir', 'plots').replace('${model_name}', f"{model_name}-debug")
        # Write debug config file to model directory
        model_dir = (project_dir / output['model_dir']).absolute()
        model_dir.mkdir(parents=True, exist_ok=True)
        debug_config_path = model_dir / 'config_debug.yaml'
        # Revert model_name to the original (without -debug suffix) before
        # writing, so repeated runs don't accumulate -debug suffixes.
        config_to_write = config.copy()
        config_to_write['model_name'] = model_name
        with open(debug_config_path, 'w') as f:
            yaml.dump(config_to_write, f)
        debug_config_file = str(debug_config_path)
        config_file = debug_config_file
    
    # Step 1: Generate profile (for both [DEPRECATED] HERACLES and STELLA)
    use_stella = config.get('hydro_backend', 'heracles').lower() == 'stella'  # DEPRECATED: default 'heracles', use 'stella'
    hydro_backend = 'stella' if use_stella else 'heracles'  # DEPRECATED: 'heracles' backend, use 'stella'

    # When skip_stella/reuse_stella is set and SN_HYDRO_DATA already exists
    # for all requested snapshots, skip profile generation — the STELLA
    # output has already been converted and the source .swd may not be
    # available at the original source_dir.
    #
    # Also skip if the .swd file exists in the model's stella directory
    # (reuse_stella=True with remade STELLA output but no SN_HYDRO_DATA yet).
    _skip_profile = False
    if skip_stella or reuse_stella:
        _cmfgen_conv_chk = config.get('cmfgen_conversion', {})
        _snap_days_chk = _cmfgen_conv_chk.get('snapshot_days', [])
        _output_cfg_chk = config.get('output', {})
        _model_dir_str_chk = _output_cfg_chk.get('model_dir', './${model_name}').replace('${model_name}', model_name)
        if str(_model_dir_str_chk).startswith('/'):
            _model_dir_chk = Path(_model_dir_str_chk)
        else:
            _model_dir_chk = (project_dir / _model_dir_str_chk).absolute()
        if _snap_days_chk and all(
            (_model_dir_chk / 'cmfgen' / f'{d:.2f}days' / 'SN_HYDRO_DATA').exists()
            for d in _snap_days_chk
        ):
            print(f"  skip_stella/reuse_stella=True and SN_HYDRO_DATA exists for all snapshots — skipping profile generation")
            _skip_profile = True
        elif reuse_stella:
            # Check if .swd exists in the model's stella directory
            _stella_out_chk = _output_cfg_chk.get('stella_output', 'stella').replace('${model_name}', model_name)
            _stella_dir_chk = _model_dir_chk / _stella_out_chk
            _swd_candidates = list(_stella_dir_chk.glob('*.swd'))
            if _swd_candidates:
                print(f"  reuse_stella=True and .swd found in {_stella_dir_chk} — skipping profile generation, will convert inline")
                _skip_profile = True

    if not _skip_profile:
        if not generate_profile(config, scripts_dir, project_dir, config_file, hydro_backend):
            print(f"ERROR: Profile generation failed for {model_name}")
            return None

    # Step 1a: Verify the generated profile actually exists on disk
    output_cfg = config.get('output', {})
    # Use the updated model_name from config (may have -debug suffix)
    current_model_name = config.get('model_name', model_name)
    base_dir_cfg = output_cfg.get('base_dir', './models').replace('${model_name}', current_model_name)
    if use_stella:
        hydro_out = output_cfg.get('stella_output', 'stella').replace('${model_name}', current_model_name)
    else:
        hydro_out = output_cfg.get('heracles_output', 'heracles').replace('${model_name}', current_model_name)  # DEPRECATED: HERACLES backend
    
    # Check if base_dir is an absolute path (for grid models)
    if str(base_dir_cfg).startswith('/'):
        # Use absolute path directly
        profile_path = Path(base_dir_cfg) / current_model_name / hydro_out / 'input_profile'
    else:
        # Use model_dir from config (handles custom paths like ./SN_A)
        model_dir_str = output_cfg.get('model_dir', '../models/${model_name}').replace('${model_name}', current_model_name)
        model_dir = (project_dir / model_dir_str).absolute()
        profile_path = model_dir / hydro_out / 'input_profile'
    if not profile_path.exists():
        # SWD snapshot mode doesn't generate an input_profile — it creates .hyd/.abn/.dat directly
        if config.get('stella_from_swd_snapshot', None):
            print(f"  SWD snapshot mode: input_profile not required (using .hyd/.abn/.dat directly)")
        else:
            print(f"ERROR: Expected profile not found after generation: {profile_path}")
            return None
    else:
        print(f"  Verified profile exists: {profile_path}")

    # Step 1a: Generate initial diagnostic plots (structure + composition)
    generate_initial_plots(config, scripts_dir, project_dir)

    # Step 1b: Read profile max radius and update Lbox_cm to match  # DEPRECATED: HERACLES-only step
    if not use_stella:
        output_cfg = config.get('output', {})
        heracles_out = output_cfg.get('heracles_output', 'heracles').replace('${model_name}', model_name)  # DEPRECATED: HERACLES backend
        base_dir_cfg = output_cfg.get('base_dir', './models').replace('${model_name}', model_name)
        if not str(heracles_out).startswith('models') and not str(heracles_out).startswith('./models') and not str(heracles_out).startswith('../'):
            profile_path = project_dir / base_dir_cfg / model_name / heracles_out / 'input_profile'
        else:
            profile_path = (project_dir / heracles_out / 'input_profile').absolute()
        
        if profile_path.exists():
            import numpy as np
            data = np.loadtxt(str(profile_path), comments='#')
            r_max = data[-1, 0]  # max radius from profile (first column)
            old_lbox = float(config.get('heracles', {}).get('Lbox_cm', 6.0e16))  # DEPRECATED: HERACLES backend
            config.setdefault('heracles', {})['Lbox_cm'] = float(r_max)  # DEPRECATED: HERACLES backend
            print(f"\nUpdated Lbox_cm: {old_lbox:.2e} -> {r_max:.2e} (from profile R_max)")
            
            # Write updated config for namelist generator
            model_dir_resolved = output_cfg.get('model_dir', './models/${model_name}').replace('${model_name}', model_name)
            model_dir_resolved = (project_dir / model_dir_resolved).absolute()
            model_dir_resolved.mkdir(parents=True, exist_ok=True)
            config_file = str(model_dir_resolved / 'config_with_lbox.yaml')
            with open(config_file, 'w') as f:
                yaml.dump(config, f)
    
    # Step 2: Generate namelist (only for [DEPRECATED] HERACLES)
    if not use_stella:
        if not generate_namelist(config, scripts_dir, project_dir, config_file, debug):
            print(f"ERROR: Namelist generation failed for {model_name}")
            return None
    
    # Determine if we are running the time-series snapshot workflow
    snapshot_days = parse_snapshot_days(config)

    # Filter snapshot days by available STELLA epochs (only when reusing existing
    # STELLA output, since the .swd file already exists). For fresh STELLA runs,
    # the epochs don't exist yet — STELLA will produce them.
    # Days that already have SN_HYDRO_DATA are also filtered against the epoch
    # tolerance — stale data from a previous homologous expansion (which is now
    # disabled) must not be silently reused if no STELLA epoch is close enough.
    if (reuse_stella or skip_stella) and (run_cmfgen or run_sedona or run_superlite):
        # Check which days already have SN_HYDRO_DATA
        output_cfg = config.get('output', {})
        model_dir_snap = output_cfg.get('model_dir', '../models/${model_name}').replace('${model_name}', model_name)
        model_dir_snap = (project_dir / model_dir_snap).absolute()
        days_with_data = []
        days_without_data = []
        for day in snapshot_days:
            day_label = f"{day:.2f}days"
            # Check both per-day subdirectory (time-series) and base cmfgen dir (debug/single)
            sn_hydro_per_day = model_dir_snap / 'cmfgen' / day_label / 'SN_HYDRO_DATA'
            sn_hydro_base = model_dir_snap / 'cmfgen' / 'SN_HYDRO_DATA'
            if sn_hydro_per_day.exists() or sn_hydro_base.exists():
                days_with_data.append(day)
            else:
                days_without_data.append(day)

        # Filter ALL days (with and without existing data) against epoch tolerance.
        # This prevents reusing stale SN_HYDRO_DATA from a previous homologous
        # expansion when no STELLA epoch is close to the requested day.
        all_filtered = filter_snapshot_days_by_stella_epochs(config, project_dir, snapshot_days)

        # Days with data that passed the filter are kept
        days_with_data_filtered = [d for d in days_with_data if d in all_filtered]
        days_with_data_skipped = [d for d in days_with_data if d not in all_filtered]
        # Days without data that passed the filter need conversion
        filtered_no_data = [d for d in days_without_data if d in all_filtered]
        days_without_data_skipped = [d for d in days_without_data if d not in all_filtered]

        snapshot_days = sorted(days_with_data_filtered + filtered_no_data)

        if days_with_data_skipped:
            print(f"  WARNING: Skipping {len(days_with_data_skipped)} day(s) with existing SN_HYDRO_DATA "
                  f"(no matching STELLA epoch — stale data from homologous expansion?): {days_with_data_skipped}")
        if days_without_data_skipped:
            print(f"  Skipped days without SN_HYDRO_DATA (no matching STELLA epoch): {days_without_data_skipped}")

        if not snapshot_days:
            print(f"ERROR: No snapshot days have a matching STELLA epoch within dt tolerance")
            return None
        if days_with_data_filtered:
            print(f"  Keeping {len(days_with_data_filtered)} day(s) with existing SN_HYDRO_DATA: {days_with_data_filtered}")

    time_series = len(snapshot_days) > 1
    if time_series:
        print(f"\nTime-series workflow enabled: snapshots at {snapshot_days} days")
        print(f"STELLA will run to {get_stella_target_days(config)} days")

    # Step 3: Prepare CMFGEN directory (copy base files) - only if running CMFGEN
    # Always use cmfgen/<day>days/ directory structure, even for single-day runs.
    # This ensures CMFGEN runs in the same directory where SN_HYDRO_DATA is written
    # by the post-STELLA conversion job, avoiding stale-file mismatches.
    #
    # Exception: J-release phase of fixed-J workflow.
    # In j_release, we must NOT wipe the directory — it contains converged
    # populations from the fixed-J run. We only update VADAT overrides.
    if run_cmfgen:
        fjw_enabled = config.get('cmfgen', {}).get('fixed_j_workflow', {}).get('enabled', False)
        skip_dir_prep = fjw_enabled and j_release
        if skip_dir_prep:
            print(f"\n  Fixed-J workflow: J-RELEASE phase — skipping directory preparation")
            print(f"    (preserving converged populations from fixed-J run)")
        if time_series and not skip_dir_prep:
            for day in snapshot_days:
                if not prepare_cmfgen_directory(config, scripts_dir, project_dir, day=day):
                    print(f"ERROR: CMFGEN directory preparation failed for {model_name} day {day}")
                    return None
        elif not time_series and not skip_dir_prep:
            # Single-day workflow: use cmfgen/<day>days/ directory
            day = snapshot_days[0] if snapshot_days else None
            if day is not None:
                if not prepare_cmfgen_directory(config, scripts_dir, project_dir, day=day):
                    print(f"ERROR: CMFGEN directory preparation failed for {model_name} day {day}")
                    return None
            else:
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

    # Step 4: Submit hydro job (STELLA or [DEPRECATED] HERACLES)
    exists, sn_hydro_path = heracles_output_exists(config, project_dir)
    use_stella = config.get('hydro_backend', 'heracles').lower() == 'stella'

    if use_stella:
        if reuse_stella:
            print(f"\nreuse_stella: True - using existing STELLA output from current model directory")
            target_stella_dir = model_dir_resolved / 'stella'
            if not target_stella_dir.exists():
                print(f"ERROR: STELLA directory not found in current model: {target_stella_dir}")
                return None
            print(f"  Using existing STELLA output from {target_stella_dir}")
            hydro_job_id = ""
        elif skip_stella:
            print(f"\n--skip-stella set: skipping STELLA job submission")
            hydro_job_id = ""
        else:
            stella_exists, stella_path = stella_output_exists(config, project_dir)
            if stella_exists and not redo_heracles:
                is_complete = stella_is_complete(config, project_dir, stella_path)
                if is_complete:
                    print(f"\nFound complete STELLA output folder: {stella_path}")
                    print(f"  Simulation reached target time, skipping STELLA job submission")
                    print(f"  (use --redo-heracles to force resubmission)")
                    extract_stella_lightcurve_pysynphot(config, scripts_dir, project_dir)
                    # Generate output plots from .swd
                    plot_stella_evolution_if_available(config, scripts_dir, project_dir)
                    plot_stella_output_plots(config, scripts_dir, project_dir)
                    print(f"  Proceeding directly to CMFGEN/SEDONA")
                    hydro_job_id = ""
                else:
                    print(f"\nFound incomplete STELLA output folder: {stella_path}")
                    print(f"  (use --redo-heracles to force resubmission)")
                    response = input("Use existing incomplete stella/ folder? [Y/n]: ").strip().lower()
                    if response != 'n':
                        print(f"Using existing stella/ folder, skipping STELLA job submission")
                        # Generate output plots (including multiband/2011dh comparison)
                        # even for incomplete runs, since .swd/.tt may still be usable
                        extract_stella_lightcurve_pysynphot(config, scripts_dir, project_dir)
                        plot_stella_evolution_if_available(config, scripts_dir, project_dir)
                        plot_stella_output_plots(config, scripts_dir, project_dir)
                        hydro_job_id = ""
                    else:
                        print(f"Proceeding with new STELLA job submission...")
                        hydro_job_id = submit_stella_job(config, scripts_dir, project_dir, debug, config_file, dependency=stella_dependency)
                        if hydro_job_id is None:
                            print(f"ERROR: STELLA job submission failed for {model_name}")
                            return None
            else:
                if stella_exists and redo_heracles:
                    print(f"\nFound existing STELLA output folder: {stella_path}")
                    print(f"  (--redo-heracles set, resubmitting anyway)")
                hydro_job_id = submit_stella_job(config, scripts_dir, project_dir, debug, config_file, dependency=stella_dependency)
                if hydro_job_id is None:
                    print(f"ERROR: STELLA job submission failed for {model_name}")
                    return None
    else:
        # [DEPRECATED] HERACLES path
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
                print(f"ERROR: [DEPRECATED] HERACLES job submission failed for {model_name}")
                return None
    heracles_job_id = hydro_job_id

    # ── Step 4a: Inline STELLA -> CMFGEN conversion ─────────────────────
    # When skip_stella or reuse_stella is set (no STELLA job submitted), the
    # STELLA .swd -> SN_HYDRO_DATA conversion runs inline here, not as a
    # SLURM job.  This calls run_inline_cmfgen_conversion() which invokes
    # build_cmfgen_conversion_cmd() -> scripts/convert_stella_swd_to_cmfgen.py
    # (or _v2.py) to produce SN_HYDRO_DATA for the requested snapshot day.
    #
    # The converter uses the base_grid_method from YAML config:
    #   'tau_amr_recursive' -> tau_amr_optimiser.recursive_amr() with hard
    #     constraints for composition interfaces, shocks, velocity extrema.
    #   'log_tau'           -> downsample_grid_tau_based() with composition-
    #     aware hybrid interpolation + adaptive_grid_refinement() AMR.
    #
    # Key sigma-preserving settings (set in YAML cmfgen_conversion):
    #   coarsen_negative_sigma=false  — preserve velocity reversal zones
    #   enable_sigma_conformance=false — do not insert sigma-driven points
    #   sigma_optimize=false           — do not modify sigma values
    #
    # This must run BEFORE SEDONA/SuperLite prep so all codes use the
    # same SN_HYDRO_DATA.
    if not time_series:
        single_day = snapshot_days[0] if snapshot_days else None
        if (reuse_stella or skip_stella) and not heracles_job_id and (run_cmfgen or run_sedona or run_superlite):
            if single_day is not None:
                # Skip inline conversion if SN_HYDRO_DATA already exists and CMFGEN
                # is not being run (e.g. SuperLite-only rerun using existing data).
                output_cfg = config.get('output', {})
                model_dir_conv = output_cfg.get('model_dir', '../models/${model_name}').replace('${model_name}', model_name)
                model_dir_conv = (project_dir / model_dir_conv).absolute()
                sn_hydro_check = model_dir_conv / 'cmfgen' / f'{single_day:.2f}days' / 'SN_HYDRO_DATA'
                if sn_hydro_check.exists() and not run_cmfgen:
                    print(f"\n  SN_HYDRO_DATA already exists at {sn_hydro_check}")
                    print(f"  CMFGEN not requested — skipping inline conversion (using existing data)")
                else:
                    print(f"\n{'='*60}")
                    print(f"Inline STELLA -> CMFGEN conversion (reuse_stella/skip_stella, no hydro job)")
                    print(f"{'='*60}")
                    if not run_inline_cmfgen_conversion(config, scripts_dir, project_dir, single_day):
                        print(f"ERROR: Inline SN_HYDRO_DATA conversion failed for {model_name}")
                        return None

        # ── Step 4a-bis: Optimize MODEL_SPEC from SN_HYDRO_DATA ───────
        # Runs scripts/optimize_model_spec_v2.py to analyse the physical
        # structure (T, rho, v, composition, tau) in SN_HYDRO_DATA and
        # generate an audit report (MODEL_SPEC.audit.md) documenting which
        # ions are needed and why.
        if run_cmfgen and single_day is not None:
            optimize_model_spec_from_hydro(config, scripts_dir, project_dir, single_day)
    else:
        single_day = None  # Set later in time-series loop

    # Step 4b: Prepare SEDONA directory if SEDONA will be run
    if run_sedona:
        # Clean up previous SEDONA runs: delete all previous SEDONA output for a fresh start.
        # Running SEDONA jobs are cancelled only if cancel_jobs is True (--keep-jobs preserves them).
        # CMFGEN/STELLA jobs are NEVER cancelled or touched.
        print(f"\n{'='*60}")
        print(f"Cleaning previous SEDONA runs for {model_name} (CMFGEN/STELLA untouched)")
        print(f"{'='*60}")
        clean_previous_sedona_runs(config, project_dir, force_cancel=cancel_jobs, verbose=True)

        if time_series:
            for day in snapshot_days:
                if not prepare_sedona_directory(config, scripts_dir, project_dir, use_stella=use_stella, day=day):
                    print(f"WARNING: SEDONA directory preparation failed for {model_name} day {day}")
        else:
            if not prepare_sedona_directory(config, scripts_dir, project_dir, use_stella=use_stella, day=single_day):
                print(f"WARNING: SEDONA directory preparation failed for {model_name}")

    # Step 4c: Prepare SIROCCO directory if SIROCCO will be run
    if run_sirocco:
        if not prepare_sirocco_directory(config, scripts_dir, project_dir):
            print(f"WARNING: SIROCCO directory preparation failed for {model_name}")

    # Step 4d: Prepare SuperLite directory if SuperLite will be run
    if run_superlite:
        print(f"\n{'='*60}")
        print(f"Cleaning previous SuperLite runs for {model_name} (CMFGEN/STELLA untouched)")
        print(f"{'='*60}")
        clean_previous_superlite_runs(config, project_dir, force_cancel=cancel_jobs, verbose=True)

        if time_series:
            for day in snapshot_days:
                if not prepare_superlite_directory(config, scripts_dir, project_dir, day=day):
                    print(f"WARNING: SuperLite directory preparation failed for {model_name} day {day}")
        else:
            single_day = snapshot_days[0] if snapshot_days else None
            if not prepare_superlite_directory(config, scripts_dir, project_dir, day=single_day):
                print(f"WARNING: SuperLite directory preparation failed for {model_name}")

    # ── Step 5: Submit spectrum generation jobs ──────────────────────────
    # CMFGEN:   scripts/submit_cmfgen_jobs.sh -> prepare_VADAT.py +
    #           prepare_IN_INTS.py + batch.sh (cmfgen_dev.exe) + OBS job.
    #           MODEL_SPEC ND is auto-updated by prepare_VADAT.py to match
    #           SN_HYDRO_DATA (ND = N_hydro + 2*NBND_INS, NP = ND + NC).
    # SEDONA:   scripts/sedona_job_template.slurm (MPI radiative transfer).
    # SIROCCO:  scripts/convert_cmfgen_to_sirocco.py + sirocco submission.
    # SuperLite: scripts/superlite_job_template.slurm (Monte Carlo transport).
    cmfgen_job_ids = None
    sedona_job_id = None
    sirocco_job_id = None
    superlite_job_id = None

    # Compute config-level relativistic flag once for CMFGEN
    config_no_relativistic = None
    if run_cmfgen:
        cmfgen_config = config.get('cmfgen', {})
        incl_rel = cmfgen_config.get('incl_rel', False)
        config_no_relativistic = not incl_rel
        # proxy_populations from config
        proxy_populations = cmfgen_config.get('proxy_populations', False)
        proxy_populations_min_ion_frac = cmfgen_config.get('proxy_populations_min_ion_frac', None)
        proxy_populations_grey_max_iter = cmfgen_config.get("proxy_populations_grey_max_iter", None)
        proxy_populations_grey_tol = cmfgen_config.get("proxy_populations_grey_tol", None)
        if config_no_relativistic:
            print(f"CMFGEN: Static snapshot mode - USE_LAM_ES=T (J=TC formal solution), no time evolution")
        else:
            print(f"CMFGEN: Relativistic effects ENABLED (incl_rel=True in config)")

        # ── Fixed-J workflow ───────────────────────────────────────────
        # Apply VADAT overrides for the requested phase (setup or j_release).
        # In 'setup' phase, also generate the EDDFACTOR before submitting.
        fjw = cmfgen_config.get('fixed_j_workflow', {})
        j_source_is_stella_edd = fjw.get('j_source', '') == 'stella_edd'
        if fjw.get('enabled', False):
            apply_fixed_j_workflow_overrides(config, j_release=j_release)
            # Write modified config to a file on the shared filesystem so the
            # compute node gets the fixed-J overrides (USE_FIXED_J, COMP_F, etc.).
            # The in-memory config dict is modified by apply_fixed_j_workflow_overrides,
            # but submit_cmfgen.sh on the compute node reads the config file from disk.
            import yaml as _yaml
            _tmp_cfg_path = project_dir / f'.{model_name}_fixed_j_config.yaml'
            with open(_tmp_cfg_path, 'w') as _tmp_cfg:
                _yaml.dump(config, _tmp_cfg, default_flow_style=False, sort_keys=False)
            config_file = str(_tmp_cfg_path)
            print(f"  Wrote fixed-J config to: {config_file}")

    if not time_series:
        # Single-day workflow: use cmfgen/<day>days/ directory (same as time-series)
        # SN_HYDRO_DATA was already created by inline conversion at Step 4a above.
        if run_cmfgen:
            # ── Fixed-J workflow: generate EDDFACTOR before submitting ──
            fjw = config.get('cmfgen', {}).get('fixed_j_workflow', {})
            if fjw.get('enabled', False) and not j_release:
                # Define cmfgen_dir for both branches (needed for auto-chain later)
                output_cfg = config.get('output', {})
                model_dir_str = output_cfg.get('model_dir', '../models/${model_name}').replace('${model_name}', model_name)
                model_dir_resolved = (project_dir / model_dir_str).absolute()
                day_label = f"{single_day:.2f}days" if single_day is not None else None
                cmfgen_dir = model_dir_resolved / 'cmfgen' / day_label if day_label else model_dir_resolved / 'cmfgen'

                # Check if STELLA is being run (EDDFACTOR from stella_edd needs STELLA output)
                if j_source_is_stella_edd and not skip_stella and hydro_job_id:
                    print(f"\n  Fixed-J workflow with stella_edd: STELLA job {hydro_job_id} submitted.")
                    print(f"  EDDFACTOR generation deferred until STELLA completes.")
                    print(f"  After STELLA finishes, rerun:")
                    print(f"    python scripts/run_model_grid.py {config_file} --skip-stella")
                    print(f"  to generate EDDFACTOR and submit CMFGEN jobs.")
                    # Skip CMFGEN submission — STELLA must complete first
                    run_cmfgen = False
                else:
                    # Setup phase: generate EDDFACTOR from 1-iteration CMFGEN run
                    if not generate_eddfactor_for_cmfgen(config, scripts_dir, project_dir, cmfgen_dir, config_file=config_file, day=single_day, dependency_job_id=heracles_job_id):
                        print(f"ERROR: EDDFACTOR generation failed for {model_name}")
                        print(f"  Cannot proceed with fixed-J workflow. Check generate_eddfactor.py output.")
                        return

            if run_cmfgen:
                print(f"\nSubmitting CMFGEN spectrum generation job...")
                cmfgen_job_ids = submit_cmfgen_jobs(config, scripts_dir, project_dir, heracles_job_id, debug, config_file, config_no_relativistic, day=single_day, proxy_populations=proxy_populations, proxy_populations_min_ion_frac=proxy_populations_min_ion_frac, proxy_populations_grey_max_iter=proxy_populations_grey_max_iter, proxy_populations_grey_tol=proxy_populations_grey_tol)
                if cmfgen_job_ids is None:
                    print(f"ERROR: CMFGEN jobs submission failed for {model_name}")
                else:
                    # ── Auto-chain: submit J-release job (stage 3) ──
                    # If fixed-J workflow is enabled with auto_chain, automatically
                    # submit the J-release run with afterok dependency on the fixed-J job.
                    fjw_check = config.get('cmfgen', {}).get('fixed_j_workflow', {})
                    if (fjw_check.get('enabled', False) and fjw_check.get('auto_chain', False)
                            and not j_release and cmfgen_job_ids.get('final')):
                        fixed_j_job_id = cmfgen_job_ids['final']
                        print(f"\n{'='*60}")
                        print(f"Auto-chain: submitting J-RELEASE job (stage 3)")
                        print(f"  Fixed-J job: {fixed_j_job_id}")
                        print(f"{'='*60}")
                        j_release_ids = submit_j_release_job(
                            config, scripts_dir, project_dir, cmfgen_dir,
                            fixed_j_job_id, config_file=config_file, day=single_day)
                        if j_release_ids is not None:
                            print(f"  J-RELEASE: {j_release_ids['final']}")
                            if j_release_ids.get('obs'):
                                print(f"  J-RELEASE OBS: {j_release_ids['obs']}")
                            cmfgen_job_ids['j_release'] = j_release_ids['final']
                            cmfgen_job_ids['j_release_obs'] = j_release_ids.get('obs')
                        else:
                            print(f"  WARNING: J-RELEASE job submission failed (non-fatal)")

        if run_sedona:
            print(f"\nSubmitting SEDONA spectrum generation job...")
            sedona_job_id = submit_sedona_job(config, scripts_dir, project_dir, heracles_job_id, debug, day=single_day)
            if sedona_job_id is None:
                print(f"ERROR: SEDONA job submission failed for {model_name}")

        if run_sirocco:
            print(f"\nSubmitting SIROCCO spectrum generation job...")
            sirocco_job_id = submit_sirocco_job(config, scripts_dir, project_dir, heracles_job_id, debug)
            if sirocco_job_id is None:
                print(f"ERROR: SIROCCO job submission failed for {model_name}")

        if run_superlite:
            print(f"\nSubmitting SuperLite spectrum generation job...")
            superlite_job_id = submit_superlite_job(config, scripts_dir, project_dir, heracles_job_id, debug, day=single_day)
            if superlite_job_id is None:
                print(f"ERROR: SuperLite job submission failed for {model_name}")
    else:
        # Time-series workflow: single post-STELLA conversion, then per-day CMFGEN/SEDONA
        all_cmfgen_ids = []
        all_sedona_ids = []
        all_superlite_ids = []
        prev_cmfgen_final = None
        prev_cmfgen_dir = None
        prev_day = None

        # One conversion job after STELLA creates all per-day SN_HYDRO_DATA and .mod files.
        # Always regenerate SN_HYDRO_DATA to ensure consistency with current config
        # (tau settings, AMR parameters, grid refinement, etc.).
        convert_job_id = None
        need_conversion = True

        # Check if all requested epochs already have SN_HYDRO_DATA files.
        # If so, skip conversion (reuse existing files from a prior conversion job).
        # This allows resubmitting CMFGEN jobs without redoing the conversion.
        force_reconvert = config.get('cmfgen_conversion', {}).get('force_reconvert', False)
        if not force_reconvert:
            existing_count = 0
            for day in snapshot_days:
                day_label = f"{day:.2f}days"
                sn_hydro = model_dir_resolved / 'cmfgen' / day_label / 'SN_HYDRO_DATA'
                if sn_hydro.exists() and sn_hydro.stat().st_size > 0:
                    existing_count += 1
            if existing_count == len(snapshot_days):
                print(f"\n  All {existing_count} epochs already have SN_HYDRO_DATA — skipping conversion.")
                need_conversion = False

        if need_conversion and (skip_stella or reuse_stella) and not heracles_job_id:
            # Remove existing SN_HYDRO_DATA files so the conversion job regenerates them
            for day in snapshot_days:
                day_label = f"{day:.2f}days"
                sn_hydro = model_dir_resolved / 'cmfgen' / day_label / 'SN_HYDRO_DATA'
                if sn_hydro.exists():
                    print(f"  Removing existing SN_HYDRO_DATA for day {day} to regenerate")
                    sn_hydro.unlink()
            print(f"\n  Regenerating all per-epoch SN_HYDRO_DATA files.")

        if need_conversion and (run_cmfgen or run_sedona or run_superlite):
            convert_job_id = submit_stella_post_conversion_job(
                config, scripts_dir, project_dir, snapshot_days,
                dependency_job_id=heracles_job_id, debug=debug
            )
            if convert_job_id is None:
                print(f"ERROR: Post-STELLA conversion job failed for {model_name}")
                return None

        for day in snapshot_days:
            print(f"\n{'='*60}")
            print(f"Submitting day {day} snapshot workflow")
            print(f"{'='*60}")

            # CMFGEN job for this day (SN_HYDRO_DATA already created by post-STELLA conversion)
            if run_cmfgen:
                if prev_cmfgen_final is None:
                    # First day: single_final (LTE cold start) + OBS
                    cmfgen_job_ids = submit_cmfgen_jobs(
                        config, scripts_dir, project_dir,
                        dependency_job_id=convert_job_id, debug=debug, config_file=config_file,
                        no_relativistic=config_no_relativistic, day=day, proxy_populations=proxy_populations,
                        proxy_populations_min_ion_frac=proxy_populations_min_ion_frac, proxy_populations_grey_max_iter=proxy_populations_grey_max_iter, proxy_populations_grey_tol=proxy_populations_grey_tol
                    )
                else:
                    # Subsequent days: seed from previous day's converged populations, then final-only
                    seed_job_id = submit_cmfgen_seed_job(
                        config, scripts_dir, project_dir, day, prev_day, prev_cmfgen_final, debug=debug
                    )
                    if seed_job_id is None:
                        print(f"ERROR: CMFGEN seed job failed for {model_name} day {day}")
                        continue
                    # CMFGEN waits for conversion (SN_HYDRO_DATA) and seed (populations)
                    # When convert_job_id is None (existing SN_HYDRO_DATA), only depend on seed
                    cmfgen_dep = f"{convert_job_id}:{seed_job_id}" if convert_job_id else str(seed_job_id)
                    cmfgen_job_ids = submit_cmfgen_jobs(
                        config, scripts_dir, project_dir,
                        dependency_job_id=cmfgen_dep, debug=debug, config_file=config_file,
                        no_relativistic=config_no_relativistic, day=day,
                        prev_cmfgen_dir=str(prev_cmfgen_dir), proxy_populations=proxy_populations,
                        proxy_populations_min_ion_frac=proxy_populations_min_ion_frac, proxy_populations_grey_max_iter=proxy_populations_grey_max_iter, proxy_populations_grey_tol=proxy_populations_grey_tol
                    )

                if cmfgen_job_ids is None:
                    print(f"ERROR: CMFGEN jobs submission failed for {model_name} day {day}")
                else:
                    prev_cmfgen_final = cmfgen_job_ids.get('obs') or cmfgen_job_ids.get('final')
                    prev_cmfgen_dir = model_dir_resolved / 'cmfgen' / f"{day:.2f}days"
                    all_cmfgen_ids.append((day, cmfgen_job_ids))

            # SEDONA job for this day (.mod already created by post-STELLA conversion)
            if run_sedona:
                sedona_job_id = submit_sedona_job(
                    config, scripts_dir, project_dir,
                    dependency_job_id=convert_job_id, debug=debug, day=day
                )
                if sedona_job_id is None:
                    print(f"ERROR: SEDONA job submission failed for {model_name} day {day}")
                else:
                    all_sedona_ids.append((day, sedona_job_id))

            # SuperLite job for this day (input.str already created by prepare_superlite_directory)
            if run_superlite:
                sl_job_id = submit_superlite_job(
                    config, scripts_dir, project_dir,
                    dependency_job_id=convert_job_id, debug=debug, day=day
                )
                if sl_job_id is None:
                    print(f"ERROR: SuperLite job submission failed for {model_name} day {day}")
                else:
                    all_superlite_ids.append((day, sl_job_id))

            prev_day = day

        # Use the last day's IDs for the summary
        if all_cmfgen_ids:
            cmfgen_job_ids = all_cmfgen_ids[-1][1]
        if all_sedona_ids:
            sedona_job_id = all_sedona_ids[-1][1]
        if all_superlite_ids:
            superlite_job_id = all_superlite_ids[-1][1]

        if run_sirocco:
            print(f"\nSubmitting SIROCCO spectrum generation job...")
            sirocco_job_id = submit_sirocco_job(config, scripts_dir, project_dir, heracles_job_id, debug)
            if sirocco_job_id is None:
                print(f"ERROR: SIROCCO job submission failed for {model_name}")

    hydro_label = 'STELLA' if use_stella else 'HERACLES'
    print(f"\n{'='*70}")
    print(f"Model {model_name} job chain submitted successfully")
    print(f"  Hydro ({hydro_label}): {heracles_job_id if heracles_job_id else '(skipped - existing data)'}")
    if cmfgen_job_ids:
        print(f"  CMFGEN:       {cmfgen_job_ids['final']}")
        print(f"  CMFGEN OBS:   {cmfgen_job_ids['obs']}")
    if sedona_job_id:
        print(f"  SEDONA:       {sedona_job_id}")
    if sirocco_job_id:
        print(f"  SIROCCO:      {sirocco_job_id}")
    if superlite_job_id:
        print(f"  SuperLite:    {superlite_job_id}")
    print(f"{'='*70}")

    return {
        'model_name': model_name,
        'heracles_job_id': heracles_job_id,
        'cmfgen_job_ids': cmfgen_job_ids,
        'sedona_job_id': sedona_job_id,
        'sirocco_job_id': sirocco_job_id,
        'superlite_job_id': superlite_job_id,
    }


def main():
    parser = argparse.ArgumentParser(
        description='Master workflow script for automated supernova model grid. '
                    'All parameters are controlled by the YAML config file.'
    )
    parser.add_argument('configs', nargs='*', help='YAML config files (or use --all)')
    parser.add_argument('--all', action='store_true', help='Process all configs in configs/')

    # ── Orchestration modes ────────────────────────────────────────────
    parser.add_argument('--dry-run', action='store_true',
                        help='Validate YAML, resolve paths, build job graph, '
                             'record all commands — no SLURM submission or '
                             'destructive operation')
    parser.add_argument('--interactive-test', action='store_true',
                        help='Submit a full validation workflow as a SLURM job '
                             'to the MPC interactive partition. Production is '
                             'blocked unless this passes.')
    parser.add_argument('--validate-only', action='store_true',
                        help='Validate YAML config only — no submission, no '
                             'file changes')
    parser.add_argument('--print-resolved-config', action='store_true',
                        help='Print the fully resolved configuration and exit')
    parser.add_argument('--status', metavar='RUN_ID', default=None,
                        help='Show workflow status for a run ID')
    parser.add_argument('--resume', metavar='RUN_ID', default=None,
                        help='Resume a workflow from its persisted state')
    parser.add_argument('--force-stage', metavar='STAGE', default=None,
                        help='Force transition to a workflow stage (with confirmation)')
    parser.add_argument('--no-submit', action='store_true',
                        help='Prepare everything but do not submit SLURM jobs')
    parser.add_argument('--stella-dependency', default=None,
                        help='SLURM dependency for STELLA job (e.g. afterok:JOBID)')
    parser.add_argument('--epochs', nargs='*', type=float, default=None,
                        help='Override snapshot_days with selected epochs only')
    parser.add_argument('--interactive-test-dir', default=None,
                        help='Directory for interactive test outputs (default: '
                             '<model_dir>/interactive_test)')
    parser.add_argument('--interactive-test-iters', type=int, default=3,
                        help='CMFGEN init iterations for interactive test (default: 3)')
    parser.add_argument('--j-release', action='store_true',
                        help='J-release phase: submit CMFGEN with USE_FIXED_J=F '
                             'using converged populations from a prior fixed-J run. '
                             'Requires fixed_j_workflow.enabled=true in config.')

    args = parser.parse_args()

    # Determine project directory
    scripts_dir = Path(__file__).parent
    project_dir = scripts_dir.parent

    # ── Import orchestration layer ─────────────────────────────────────
    try:
        sys.path.insert(0, str(scripts_dir))
        from orchestration.config import resolve_config, ConfigError
        from orchestration.state import StateMachine, WorkflowStates
        from orchestration.manifests import ManifestWriter
        from orchestration.interactive_test import InteractiveTestRunner
        from orchestration.monitor import SlurmMonitor, JobState
        HAS_ORCHESTRATION = True
    except ImportError as e:
        print(f"WARNING: Orchestration layer not available: {e}")
        HAS_ORCHESTRATION = False

    # ── Handle --status ────────────────────────────────────────────────
    if args.status and HAS_ORCHESTRATION:
        _handle_status(args.status, project_dir)
        return

    # ── Handle --resume ────────────────────────────────────────────────
    if args.resume and HAS_ORCHESTRATION:
        _handle_resume(args.resume, project_dir, scripts_dir)
        return

    # ── Handle --force-stage ───────────────────────────────────────────
    if args.force_stage and HAS_ORCHESTRATION:
        _handle_force_stage(args.force_stage, project_dir)
        return

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
    if args.dry_run:
        print(f"Mode: DRY RUN (no submission)")
    elif args.validate_only:
        print(f"Mode: VALIDATE ONLY")
    elif args.interactive_test:
        print(f"Mode: MPC INTERACTIVE TEST")
    elif args.no_submit:
        print(f"Mode: NO SUBMIT (prepare only)")
    else:
        print(f"Mode: PRODUCTION")
    print(f"{'='*70}")

    # ── Handle --validate-only ─────────────────────────────────────────
    if args.validate_only:
        all_valid = True
        for config_file in config_files:
            print(f"\nValidating: {config_file}")
            if HAS_ORCHESTRATION:
                try:
                    cfg = resolve_config(str(config_file), str(project_dir), str(scripts_dir))
                    print(f"  PASSED: {cfg.model_name} (hash={cfg.config_hash})")
                    print(f"  Snapshot days: {cfg.snapshot_days}")
                    print(f"  STELLA target: {cfg.stella_target_days}d")
                    print(f"  Site: {cfg.site.cluster_name} "
                          f"interactive={cfg.site.interactive_partition} "
                          f"account={cfg.site.account}")
                except ConfigError as e:
                    print(f"  FAILED: {e}")
                    all_valid = False
            else:
                # Fallback: just load YAML
                config = load_yaml_config(config_file)
                print(f"  Loaded: {config.get('model_name', 'unknown')}")
        sys.exit(0 if all_valid else 1)

    # ── Handle --print-resolved-config ─────────────────────────────────
    if args.print_resolved_config:
        for config_file in config_files:
            if HAS_ORCHESTRATION:
                try:
                    cfg = resolve_config(str(config_file), str(project_dir), str(scripts_dir))
                    print(f"\n# Resolved config for {config_file}")
                    print(cfg.to_yaml())
                except ConfigError as e:
                    print(f"ERROR: {e}")
                    sys.exit(1)
            else:
                print("ERROR: Orchestration layer required for --print-resolved-config")
                sys.exit(1)
        return

    # ── Handle --dry-run ───────────────────────────────────────────────
    if args.dry_run:
        for config_file in config_files:
            print(f"\n{'='*60}")
            print(f"Dry run: {config_file}")
            print(f"{'='*60}")
            if HAS_ORCHESTRATION:
                try:
                    cfg = resolve_config(str(config_file), str(project_dir), str(scripts_dir))
                    print(f"  Model: {cfg.model_name}")
                    print(f"  Hash: {cfg.config_hash}")
                    print(f"  Snapshot days: {cfg.snapshot_days}")
                    print(f"  STELLA target: {cfg.stella_target_days}d")
                    print(f"  STELLA partition: {cfg.stella_partition}")
                    print(f"  CMFGEN partition: {cfg.cmfgen_partition}")
                    print(f"  Site: {cfg.site.cluster_name}")

                    # Build job graph (no submission)
                    jobs = _build_job_graph(cfg)
                    print(f"\n  Job graph ({len(jobs)} jobs):")
                    for j in jobs:
                        dep = j.get('dependency', 'none')
                        print(f"    {j['name']:30s} partition={j['partition']:12s} "
                              f"dep={dep}")

                    # Write manifests
                    run_dir = Path(cfg.model_dir) / 'output' / cfg.run_id
                    sm = StateMachine(cfg.run_id, cfg.model_name, cfg.config_hash,
                                     str(run_dir / 'workflow_state.json'))
                    sm.transition(WorkflowStates.CONFIG_VALIDATED,
                                 message="Dry run: config validated")
                    mw = ManifestWriter(str(run_dir), cfg.run_id)
                    mw.write_all(cfg, sm, jobs=jobs)
                    mw.log(f"Dry run complete for {cfg.model_name}")
                    print(f"\n  Manifests written to: {run_dir}")
                except ConfigError as e:
                    print(f"  CONFIG ERROR: {e}")
            else:
                print("  ERROR: Orchestration layer required for --dry-run")
        return

    # ── Handle --interactive-test ──────────────────────────────────────
    if args.interactive_test:
        if not HAS_ORCHESTRATION:
            print("ERROR: Orchestration layer required for --interactive-test")
            sys.exit(1)
        all_passed = True
        for config_file in config_files:
            print(f"\n{'='*60}")
            print(f"MPC Interactive Test: {config_file}")
            print(f"{'='*60}")
            try:
                cfg = resolve_config(str(config_file), str(project_dir), str(scripts_dir))
                # Determine test directory
                if args.interactive_test_dir:
                    test_dir = args.interactive_test_dir
                else:
                    test_dir = str(Path(cfg.model_dir) / 'interactive_test')
                runner = InteractiveTestRunner(
                    cfg, test_dir,
                    init_iterations=args.interactive_test_iters,
                )
                result = runner.run()
                if result.passed:
                    print(f"\n  INTERACTIVE TEST PASSED")
                    print(f"  Job ID: {result.job_id}")
                    print(f"  State: {result.final_state}")
                    print(f"  Elapsed: {result.elapsed}")
                    print(f"  Gates: {sum(result.gates.values())}/{len(result.gates)} passed")
                else:
                    print(f"\n  INTERACTIVE TEST FAILED")
                    print(f"  Job ID: {result.job_id}")
                    print(f"  State: {result.final_state}")
                    print(f"  Exit code: {result.exit_code}")
                    print(f"  Failure reasons:")
                    for r in result.failure_reasons:
                        print(f"    - {r}")
                    print(f"  Test dir preserved: {result.test_dir}")
                    all_passed = False
            except ConfigError as e:
                print(f"  CONFIG ERROR: {e}")
                all_passed = False
            except Exception as e:
                print(f"  ERROR: {e}")
                all_passed = False

        if all_passed:
            print(f"\n{'='*70}")
            print(f"All interactive tests PASSED — production submission allowed")
            print(f"{'='*70}")
        else:
            print(f"\n{'='*70}")
            print(f"Interactive test FAILED — production submission BLOCKED")
            print(f"Fix the issues above before submitting production jobs")
            print(f"{'='*70}")
            sys.exit(1)
        return

    # ── Normal production workflow ─────────────────────────────────────
    # If --no-submit, set a flag that process_config can check
    if args.no_submit:
        os.environ['RUN_MODEL_GRID_NO_SUBMIT'] = '1'

    if args.stella_dependency:
        os.environ['RUN_MODEL_GRID_STELLA_DEP'] = args.stella_dependency
    # Process each config
    results = []
    for config_file in config_files:
        # Write manifests if orchestration is available
        if HAS_ORCHESTRATION:
            try:
                cfg = resolve_config(str(config_file), str(project_dir), str(scripts_dir))
                run_dir = Path(cfg.model_dir) / 'output' / cfg.run_id
                sm = StateMachine(cfg.run_id, cfg.model_name, cfg.config_hash,
                                 str(run_dir / 'workflow_state.json'))
                sm.transition(WorkflowStates.CONFIG_VALIDATED,
                             message="Production workflow started")
                mw = ManifestWriter(str(run_dir), cfg.run_id)
                mw.write_resolved_config(cfg)
                mw.write_epoch_map(cfg)
                mw.log(f"Production workflow started for {cfg.model_name}")
            except Exception as e:
                print(f"  WARNING: Manifest writing failed: {e}")

        result = process_config(config_file, scripts_dir, project_dir)
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
            print(f"    CMFGEN:     {result['cmfgen_job_ids']['final']}")
            print(f"    CMFGEN OBS: {result['cmfgen_job_ids']['obs']}")
        if result.get('sedona_job_id'):
            print(f"    SEDONA:     {result['sedona_job_id']}")
        if result.get('sirocco_job_id'):
            print(f"    SIROCCO:    {result['sirocco_job_id']}")
        if result.get('superlite_job_id'):
            print(f"    SuperLite:  {result['superlite_job_id']}")
    print(f"\nMonitor with: squeue -u $USER")
    print(f"{'='*70}")


def _build_job_graph(cfg):
    """Build the SLURM job dependency graph from resolved config (dry-run)."""
    jobs = []
    dep = None

    # STELLA job
    if not cfg.reuse_stella and not cfg.skip_stella:
        jobs.append({
            'name': f'stella_{cfg.model_name}',
            'stage': 'STELLA',
            'job_id': '(dry-run)',
            'partition': cfg.stella_partition,
            'dependency': None,
            'state': 'PLANNED',
        })
        dep = 'afterok:stella'

    # Conversion job (time-series)
    if cfg.time_series and (cfg.run_cmfgen or cfg.run_sedona or cfg.run_superlite):
        jobs.append({
            'name': f'convert_{cfg.model_name}',
            'stage': 'CONVERSION',
            'job_id': '(dry-run)',
            'partition': cfg.stella_partition,
            'dependency': dep,
            'state': 'PLANNED',
        })
        dep = f'afterok:convert'

    # Per-epoch jobs
    for epoch in cfg.epoch_mappings:
        if cfg.run_cmfgen:
            jobs.append({
                'name': f'cmfgen_{cfg.model_name}_{epoch.output_label}',
                'stage': 'CMFGEN',
                'job_id': '(dry-run)',
                'partition': cfg.cmfgen_partition,
                'dependency': dep,
                'state': 'PLANNED',
            })
        if cfg.run_sedona:
            jobs.append({
                'name': f'sedona_{cfg.model_name}_{epoch.output_label}',
                'stage': 'SEDONA',
                'job_id': '(dry-run)',
                'partition': cfg.sedona.get('partition', 'general'),
                'dependency': dep,
                'state': 'PLANNED',
            })
        if cfg.run_superlite:
            jobs.append({
                'name': f'superlite_{cfg.model_name}_{epoch.output_label}',
                'stage': 'SuperLite',
                'job_id': '(dry-run)',
                'partition': cfg.superlite.get('partition', 'general'),
                'dependency': dep,
                'state': 'PLANNED',
            })

    return jobs


def _handle_status(run_id, project_dir):
    """Handle --status RUN_ID: show workflow state."""
    import json
    # Search for workflow_state.json in model directories
    found = False
    for model_dir in Path(project_dir).iterdir():
        if not model_dir.is_dir():
            continue
        for state_file in model_dir.rglob('workflow_state.json'):
            try:
                data = json.loads(state_file.read_text())
                if data.get('run_id', '').startswith(run_id) or run_id in data.get('run_id', ''):
                    print(f"\nWorkflow State: {data['run_id']}")
                    print(f"  Model: {data['model_name']}")
                    print(f"  Current state: {data['current_state']}")
                    print(f"  Config hash: {data['config_hash']}")
                    print(f"  Created: {data.get('created_at', '?')}")
                    print(f"  Updated: {data.get('updated_at', '?')}")
                    print(f"  Job IDs:")
                    for stage, jid in data.get('job_ids', {}).items():
                        print(f"    {stage}: {jid}")
                    print(f"  Epoch status:")
                    for epoch, status in data.get('epoch_status', {}).items():
                        print(f"    {epoch}: {status}")
                    print(f"  Transitions ({len(data.get('transitions', []))}):")
                    for t in data.get('transitions', []):
                        print(f"    [{t['timestamp']}] {t['from_state']} -> {t['to_state']}: {t.get('message', '')}")
                    print(f"  State file: {state_file}")
                    found = True
            except (json.JSONDecodeError, KeyError):
                pass
    if not found:
        print(f"No workflow state found for run ID: {run_id}")
        print(f"Searched in: {project_dir}")


def _handle_resume(run_id, project_dir, scripts_dir):
    """Handle --resume RUN_ID: resume workflow from persisted state."""
    import json
    # Find the state file
    for model_dir in Path(project_dir).iterdir():
        if not model_dir.is_dir():
            continue
        for state_file in model_dir.rglob('workflow_state.json'):
            try:
                data = json.loads(state_file.read_text())
                if run_id in data.get('run_id', ''):
                    print(f"Resuming workflow: {data['run_id']}")
                    print(f"  Current state: {data['current_state']}")
                    print(f"  Model: {data['model_name']}")
                    # Find the config file
                    config_path = None
                    for cf in (Path(project_dir) / 'configs').glob('*.yaml'):
                        try:
                            from orchestration.config import resolve_config
                            cfg = resolve_config(str(cf), str(project_dir), str(scripts_dir))
                            if cfg.config_hash == data['config_hash']:
                                config_path = cf
                                break
                        except Exception:
                            pass
                    if config_path:
                        print(f"  Config: {config_path}")
                        print(f"  Resuming by re-running process_config...")
                        result = process_config(config_path, scripts_dir, project_dir)
                        if result:
                            print(f"  Resume complete")
                    else:
                        print(f"  ERROR: Could not find matching config file")
                    return
            except (json.JSONDecodeError, KeyError):
                pass
    print(f"No workflow state found for run ID: {run_id}")


def _handle_force_stage(stage, project_dir):
    """Handle --force-stage STAGE: force transition to a stage."""
    import json
    print(f"Force transition to: {stage}")
    response = input(f"Are you sure you want to force stage to {stage}? [y/N]: ").strip().lower()
    if response != 'y':
        print("Aborted")
        return
    # Find and update state files
    for model_dir in Path(project_dir).iterdir():
        if not model_dir.is_dir():
            continue
        for state_file in model_dir.rglob('workflow_state.json'):
            try:
                data = json.loads(state_file.read_text())
                from orchestration.state import StateMachine, WorkflowStates
                sm = StateMachine(data['run_id'], data['model_name'],
                                 data['config_hash'], str(state_file))
                try:
                    target = WorkflowStates(stage)
                except ValueError:
                    print(f"Unknown stage: {stage}")
                    print(f"Valid stages: {[s.value for s in WorkflowStates]}")
                    return
                sm.transition(target, force=True, message="Forced by --force-stage")
                print(f"  Forced {data['run_id']} -> {stage}")
            except Exception as e:
                print(f"  Error updating {state_file}: {e}")


if __name__ == '__main__':
    main()
