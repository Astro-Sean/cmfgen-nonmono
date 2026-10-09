#!/usr/bin/env python3
"""Generate SN_A_pba batch-2 configs: lower E/Ni + envelope variations.

Follow-up to generate_sn_a_pba_grid.py.  The first grid (E 5/9/13 x Ni56
0.04/0.075/0.15) showed all models decay too slowly post-peak and are
~0.4-0.6 mag too blue in g.  This batch explores three directions:

  1. Lower E x lower Ni  (redder photosphere + faster-decaying tail):
       E03/E04 x Ni 0.02/0.04  +  E05_Ni0p02
  2. Less H-envelope mass at the E05_Ni0p04 baseline
       (faster post-peak decay):  h_ext_mass_msun 0.10 -> 0.05, 0.03
  3. Larger H-envelope radius at the E05_Ni0p04 baseline
       (larger R_ph -> redder):  h_ext_r_outer_rsun 200 -> 400, 500

Naming: E<ee>_Ni<n>[_M<env mass>|_R<env radius>].

Same base as batch 1: SN_A_v4 physics, parallel OpenStella, STELLA-only.
"""
import os
import copy
import yaml
from pathlib import Path

BASE_DIR = Path("/viper/ptmp1/sbrennan/SN2024pba")
GRID_DIR = BASE_DIR / "SN_A_pba_grid"
TEMPLATE = BASE_DIR / "SN_A_v4" / "config.yaml"
PRESN_SOURCE = BASE_DIR / "SN_A" / "s12.0-presn-stripped"
PARALLEL_STELLA = "/ptmp/sbrennan/stella/OpenStella_parallel"

# (name, E_burst, Ni56, env_mass_msun or None, env_r_out_rsun or None)
MODELS = [
    ("E03_Ni0p02", 3.0, 0.02, None, None),
    ("E03_Ni0p04", 3.0, 0.04, None, None),
    ("E04_Ni0p02", 4.0, 0.02, None, None),
    ("E04_Ni0p04", 4.0, 0.04, None, None),
    ("E05_Ni0p02", 5.0, 0.02, None, None),
    ("E05_Ni0p04_M0p05", 5.0, 0.04, 0.05, None),
    ("E05_Ni0p04_M0p03", 5.0, 0.04, 0.03, None),
    ("E05_Ni0p04_R400", 5.0, 0.04, None, 400.0),
    ("E05_Ni0p04_R500", 5.0, 0.04, None, 500.0),
]


def make_config(base, model_name, e_burst, ni56, env_mass, env_r):
    cfg = copy.deepcopy(base)

    env_desc = "0.10 Msun H env at 200 Rsun"
    if env_mass is not None:
        env_desc = f"{env_mass} Msun H env at 200 Rsun"
    if env_r is not None:
        env_desc = f"0.10 Msun H env at {env_r:g} Rsun"

    cfg["model_name"] = model_name
    cfg["description"] = (
        f"SN_A_pba batch2: SN_A_v4 physics (1.6 Msun He core + {env_desc}) "
        f"with E_burst={e_burst} (~{e_burst/13:.1f} foe) and Ni56={ni56} "
        f"Msun. Target: SN 2024pba first-peak g/r/i lightcurve. "
        f"Parallel STELLA, 65 d."
    )

    cfg["run_cmfgen"] = False
    cfg["run_sirocco"] = False
    cfg["run_sedona"] = False
    cfg["run_superlite"] = False

    cfg["skip_stella"] = False
    cfg["reuse_stella"] = False
    cfg["restart_stella"] = False
    cfg["cancel_jobs"] = True
    cfg["force_cancel"] = True
    cfg["plot_lightcurve"] = True
    cfg["debug"] = False

    cfg["stella"]["explosion_energy"] = e_burst
    cfg["stella"]["nickel56"]["total_mass_msun"] = ni56
    if env_mass is not None:
        cfg["kepler_presn"]["h_ext_mass_msun"] = env_mass
    if env_r is not None:
        cfg["kepler_presn"]["h_ext_r_outer_rsun"] = env_r

    cfg["stella"]["stella_root"] = PARALLEL_STELLA
    cfg["stella"]["cpus_per_task"] = 3

    cfg["kepler_presn"]["presn_file"] = (
        f"/ptmp/sbrennan/SN2024pba/SN_A_pba_grid/{model_name}/s12.0-presn-stripped")
    cfg["output"]["model_dir"] = f"./SN_A_pba_grid/{model_name}"

    return cfg


def main():
    with open(TEMPLATE) as f:
        base = yaml.safe_load(f)

    GRID_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Generating {len(MODELS)} SN_A_pba batch-2 configs in {GRID_DIR}/\n")
    print(f"{'Model':<20} {'E_burst':>8} {'Ni56':>6} {'env_M':>6} {'env_R':>6}")
    print("-" * 50)

    configs = []
    for name, e, ni, env_mass, env_r in MODELS:
        model_dir = GRID_DIR / name
        model_dir.mkdir(parents=True, exist_ok=True)

        kep_link = model_dir / "s12.0-presn-stripped"
        if not kep_link.exists():
            os.symlink(PRESN_SOURCE, kep_link)

        cfg = make_config(base, name, e, ni, env_mass, env_r)
        config_path = model_dir / "config.yaml"

        header = (
            f"# SN_A_pba batch-2 model: {name}\n"
            f"# Base: SN_A_v4 (KEPLER s12.0 stripped + H env)\n"
            f"# Axes: E_burst={e}, Ni56={ni} Msun"
            + (f", env_mass={env_mass} Msun" if env_mass else "")
            + (f", env_r_out={env_r:g} Rsun" if env_r else "") + "\n"
            f"# Target: SN 2024pba first-peak g/r/i lightcurve\n"
            f"# Codes: STELLA only (parallel OpenStella build)\n\n"
        )
        with open(config_path, "w") as f:
            f.write(header)
            yaml.dump(cfg, f, default_flow_style=False, sort_keys=False,
                      width=120)

        configs.append(str(config_path))
        print(f"{name:<20} {e:8.1f} {ni:6.3f} "
              f"{env_mass if env_mass else '-':>6} "
              f"{env_r if env_r else '-':>6}")

    print(f"\nGenerated {len(configs)} configs.")
    print("\nSubmit (chained behind the running batch-1 queue):")
    print("  python scripts/run_model_grid.py "
          "SN_A_pba_grid/<new>/config.yaml ... "
          "--stella-dependency afterany:<last_batch1_jobid>")


if __name__ == "__main__":
    main()
