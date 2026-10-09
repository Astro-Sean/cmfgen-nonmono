#!/usr/bin/env python3
"""
Setup SN_B model from a config.yaml file.

This is a convenience wrapper around prepare_sn_b.py that reads all parameters
from a YAML config file instead of requiring CLI arguments. It supports two
modes for specifying the source SN_A ejecta:

  1. Direct .swd file path (for experimentation):
       sn_b_from_sn_a:
         sn_a_swd: "/path/to/SN_A.swd"

  2. SN_A model directory (auto-discovers the .swd):
       sn_b_from_sn_a:
         sn_a_dir: "/path/to/SN_A"

This allows you to iterate on the SN_A STELLA run and re-generate the SN_B
profile from any .swd snapshot without needing a full model directory structure.

Usage:
    python setup_SN_B.py <config.yaml>

Example config.yaml:
    model_name: "SN_B"
    sn_b_from_sn_a:
      sn_a_swd: "/ptmp/sbrennan/SN2024pba/SN_A/stella/SN_A.swd"
      sn_a_age_days: 349.0
      shell_mass_msun: 1.0
      shell_velocity_kms: 5000
      shell_xh: 0.7
      shell_yhe: 0.3
    stella:
      target_days: 50.0
      nzon: 500
    cmfgen_conversion:
      snapshot_days: [40]
    output:
      model_dir: "./SN_B"
"""

import argparse
import sys
import subprocess
from pathlib import Path

import yaml

# Import shared utilities from prepare_sn_b.py
sys.path.insert(0, str(Path(__file__).parent))
from prepare_sn_b import (
    create_hhe_shell,
    combine_shell_and_csm,
    write_sn_b_hyd,
    write_sn_b_abn,
    write_sn_b_xni,
    write_sn_b_dat,
    write_sn_b_input_profile,
)
from convert_stella_swd_to_cmfgen import read_swd, normalize_swd_composition
from homologous_evolution import find_closest_epoch, homologous_evolve


def read_sn_a_ejecta_from_swd(swd_path, sn_a_age=349, tolerance=20):
    """Read SN_A ejecta directly from a .swd file path.

    This bypasses the directory-based lookup in prepare_sn_b.read_sn_a_ejecta()
    and reads the .swd file directly, allowing experimentation with arbitrary
    .swd files without a full model directory.
    """
    swd_path = Path(swd_path)
    if not swd_path.exists():
        raise FileNotFoundError(f"SN_A .swd not found: {swd_path}")

    print(f"  Reading SN_A .swd: {swd_path}")
    epochs = read_swd(str(swd_path))
    available = sorted(epochs.keys())
    print(f"  Available epochs: {available[0]:.2f} ... {available[-1]:.2f} days "
          f"({len(available)} total)")

    # If requested age is beyond available data, use the latest epoch
    latest_epoch = available[-1]
    if sn_a_age > latest_epoch + tolerance:
        print(f"  Requested age {sn_a_age} days is beyond available data "
              f"(latest: {latest_epoch:.2f} days). Using latest epoch and "
              f"homologously evolving.")
        t_source = latest_epoch
        epoch_data = epochs[t_source]
    else:
        t_source, epoch_data = find_closest_epoch(epochs, sn_a_age, tolerance)

    # Homologously evolve if needed
    if abs(t_source - sn_a_age) > 1.0:
        print(f"  Homologously evolving from {t_source:.2f} to {sn_a_age} days...")
        epoch_data = homologous_evolve(epoch_data, t_source, sn_a_age)
    else:
        print(f"  Using epoch at {t_source:.2f} days (matches requested age)")

    # Normalize composition
    if epoch_data.get('composition') is not None:
        comp = normalize_swd_composition(epoch_data['composition'])
        epoch_data['fH'] = comp['H']
        epoch_data['fHe'] = comp['He']
        epoch_data['fC'] = comp['C']
        epoch_data['fO'] = comp['O']
        epoch_data['fSi'] = comp['Si']
        epoch_data['fCa'] = comp['Ca']
        epoch_data['fFe'] = comp['Fe']
        epoch_data['fNi'] = comp['Ni']
    else:
        for key in ['fH', 'fHe', 'fC', 'fO', 'fSi', 'fCa', 'fFe', 'fNi']:
            epoch_data[key] = np.zeros(len(epoch_data['r']))
        epoch_data['fHe'] = np.ones(len(epoch_data['r']))

    return epoch_data


def resolve_swd_path(sn_b_cfg, project_dir):
    """Resolve the .swd file path from either sn_a_swd or sn_a_dir.

    Priority:
      1. sn_a_swd: direct path to .swd file
      2. sn_a_dir: directory containing stella/<model>.swd or stella/strad/<model>.swd
    """
    import numpy as np  # needed by read_sn_a_ejecta_from_swd

    sn_a_swd = sn_b_cfg.get('sn_a_swd', None)
    sn_a_dir = sn_b_cfg.get('sn_a_dir', None)

    if sn_a_swd:
        swd_path = Path(sn_a_swd)
        if not swd_path.is_absolute():
            swd_path = (project_dir / swd_path).resolve()
        if not swd_path.exists():
            raise FileNotFoundError(f"sn_a_swd path does not exist: {swd_path}")
        return swd_path

    if sn_a_dir:
        sn_a_path = Path(sn_a_dir)
        if not sn_a_path.is_absolute():
            sn_a_path = (project_dir / sn_a_path).resolve()
        model_name = sn_a_path.name
        # Try stella/strad/<model>.swd first, then stella/<model>.swd
        for sub in ['stella/strad', 'stella']:
            candidate = sn_a_path / sub / f'{model_name}.swd'
            if candidate.exists():
                return candidate
        raise FileNotFoundError(
            f"Could not find .swd in sn_a_dir '{sn_a_path}'. "
            f"Looked for stella/strad/{model_name}.swd and stella/{model_name}.swd. "
            f"Use sn_a_swd to specify the .swd path directly."
        )

    raise ValueError("config must specify either sn_b_from_sn_a.sn_a_swd or sn_b_from_sn_a.sn_a_dir")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                 formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('config', help='Path to config.yaml file')
    parser.add_argument('--force', '-f', action='store_true',
                        help='Overwrite existing output directory')
    args = parser.parse_args()

    config_path = Path(args.config).resolve()
    if not config_path.exists():
        print(f"ERROR: Config file not found: {config_path}")
        sys.exit(1)

    project_dir = config_path.parent

    with open(config_path) as f:
        config = yaml.safe_load(f)

    model_name = config.get('model_name', 'SN_B')
    sn_b_cfg = config.get('sn_b_from_sn_a', {})
    stella_cfg = config.get('stella', {})
    output_cfg = config.get('output', {})

    # --- Extract parameters ---
    sn_a_age = sn_b_cfg.get('sn_a_age_days', 349.0)
    shell_mass = sn_b_cfg.get('shell_mass_msun', 1.0)
    shell_velocity = sn_b_cfg.get('shell_velocity_kms', 5000)
    shell_xh = sn_b_cfg.get('shell_xh', 0.7)
    shell_yhe = sn_b_cfg.get('shell_yhe', 0.3)
    explosion_energy = sn_b_cfg.get('explosion_energy', 0.0)
    aq = sn_b_cfg.get('aq', 4.0)
    fit_tau = sn_b_cfg.get('fit_tau', 5.0)
    meth = sn_b_cfg.get('meth', 2)

    nzon = stella_cfg.get('nzon', 500)
    target_days = stella_cfg.get('target_days', 50.0)

    cmfgen_conv = config.get('cmfgen_conversion', {})
    snapshot_days = cmfgen_conv.get('snapshot_days', [40])
    snapshot_day = snapshot_days[0] if isinstance(snapshot_days, list) else snapshot_days

    # --- Resolve output directory ---
    model_dir = output_cfg.get('model_dir', './${model_name}').replace('${model_name}', model_name)
    if str(model_dir).startswith('/'):
        model_dir_abs = Path(model_dir)
    else:
        model_dir_abs = (project_dir / model_dir).resolve()

    # --- Resolve .swd path ---
    swd_path = resolve_swd_path(sn_b_cfg, project_dir)

    print(f"\n=== Setting up SN_B model: {model_name} ===")
    print(f"  Config: {config_path}")
    print(f"  Source .swd: {swd_path}")
    print(f"  SN_A age: {sn_a_age} days")
    print(f"  Shell: M={shell_mass} Msun, v={shell_velocity} km/s, "
          f"X={shell_xh}, Y={shell_yhe}")
    print(f"  STELLA target: {target_days} days, Snapshot: {snapshot_day} days")
    print(f"  Output: {model_dir_abs}")

    # --- Create output directories ---
    stella_dir = model_dir_abs / 'stella'
    modmake_dir = stella_dir / 'modmake'
    plots_dir = model_dir_abs / 'plots'
    modmake_dir.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(parents=True, exist_ok=True)

    # --- Step 1: Read SN_A ejecta from .swd ---
    print(f"\n--- Step 1: Reading SN_A ejecta at {sn_a_age} days ---")
    csm = read_sn_a_ejecta_from_swd(swd_path, sn_a_age=sn_a_age)

    print(f"  Ejecta: {len(csm['r'])} zones, "
          f"r=[{csm['r'].min():.3e}, {csm['r'].max():.3e}] cm")
    print(f"  Ejecta: v=[{csm['v'].min():.1f}, {csm['v'].max():.1f}] km/s")
    print(f"  Ejecta: rho=[{csm['rho'].min():.3e}, {csm['rho'].max():.3e}] g/cc")

    # Convert velocities from km/s to cm/s
    csm['v'] = csm['v'] * 1e5  # km/s -> cm/s

    # --- Step 2: Create H/He shell ---
    print(f"\n--- Step 2: Creating H/He shell ---")
    shell = create_hhe_shell(shell_mass, shell_velocity, shell_xh, shell_yhe)

    # --- Step 3: Combine shell + CSM ---
    print(f"\n--- Step 3: Combining shell + CSM ---")
    combined = combine_shell_and_csm(shell, csm, nzon_total=nzon)

    print(f"  Combined: {len(combined['r'])} zones")
    print(f"  r: [{combined['r'].min():.3e}, {combined['r'].max():.3e}] cm")
    print(f"  v: [{combined['v'].min():.1f}, {combined['v'].max():.1f}] cm/s "
          f"([{combined['v'].min()/1e5:.1f}, {combined['v'].max()/1e5:.1f}] km/s)")
    print(f"  rho: [{combined['rho'].min():.3e}, {combined['rho'].max():.3e}] g/cc")

    # --- Step 4: Write STELLA input files ---
    print(f"\n--- Step 4: Writing STELLA input files ---")
    hyd_path = modmake_dir / f'{model_name}.hyd'
    abn_path = modmake_dir / f'{model_name}.abn'
    xni_path = modmake_dir / f'{model_name}.xni'
    dat_path = modmake_dir / f'{model_name}.dat'
    ip_path = modmake_dir / 'input_profile'

    write_sn_b_hyd(str(hyd_path), combined)
    write_sn_b_abn(str(abn_path), combined)
    write_sn_b_xni(str(xni_path), nzon)
    write_sn_b_dat(str(dat_path), model_name, target_days, snapshot_day,
                   explosion_energy=explosion_energy, aq=aq)
    write_sn_b_input_profile(str(ip_path), combined)

    # Copy input_profile to stella dir
    import shutil
    shutil.copy2(str(ip_path), str(stella_dir / 'input_profile'))

    # --- Step 5: Generate plots ---
    print(f"\n--- Step 5: Generating plots ---")
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        fig.suptitle(f'{model_name} - Initial Model (Shell + SN_A CSM)',
                     fontsize=14, fontweight='bold')
        r_rsun = combined['r'] / 6.96e10

        axes[0, 0].loglog(r_rsun, combined['rho'], 'b-', linewidth=0.8)
        axes[0, 0].set_xlabel('Radius (R_sun)')
        axes[0, 0].set_ylabel('Density (g/cm³)')
        axes[0, 0].set_title('Density Profile')
        axes[0, 0].grid(True, alpha=0.3)

        axes[0, 1].semilogx(r_rsun, combined['v'] / 1e5, 'r-', linewidth=0.8)
        axes[0, 1].set_xlabel('Radius (R_sun)')
        axes[0, 1].set_ylabel('Velocity (km/s)')
        axes[0, 1].set_title('Velocity Profile')
        axes[0, 1].grid(True, alpha=0.3)

        axes[1, 0].semilogx(r_rsun, combined['T'], 'g-', linewidth=0.8)
        axes[1, 0].set_xlabel('Radius (R_sun)')
        axes[1, 0].set_ylabel('Temperature (K)')
        axes[1, 0].set_title('Temperature Profile')
        axes[1, 0].set_yscale('log')
        axes[1, 0].grid(True, alpha=0.3)

        for key, color, label in [('fH', 'blue', 'H'), ('fHe', 'orange', 'He'),
                                   ('fO', 'green', 'O'), ('fSi', 'red', 'Si'),
                                   ('fFe', 'purple', 'Fe'), ('fCa', 'brown', 'Ca')]:
            if np.any(combined[key] > 1e-10):
                axes[1, 1].semilogx(r_rsun, combined[key], color=color,
                                    label=label, linewidth=0.8)
        axes[1, 1].set_xlabel('Radius (R_sun)')
        axes[1, 1].set_ylabel('Mass Fraction')
        axes[1, 1].set_title('Composition Profile')
        axes[1, 1].set_yscale('log')
        axes[1, 1].set_ylim(1e-10, 1)
        axes[1, 1].legend(fontsize=8)
        axes[1, 1].grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(str(plots_dir / 'sn_b_initial.png'), dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  Saved: {plots_dir}/sn_b_initial.png")
    except Exception as e:
        print(f"  Plot generation failed: {e}")

    print(f"\n=== SN_B setup complete ===")
    print(f"  STELLA input files in: {modmake_dir}")
    print(f"  Run full pipeline with: python scripts/run_model_grid.py {config_path}")


if __name__ == '__main__':
    import numpy as np  # needed for plot generation
    main()
