#!/usr/bin/env python3
"""Generate SN_A variant models: S (surface-layer) and X (extended-env) families.

Two 3-model families, each a single-parameter perturbation of a proven basis:

  SN_A_S_* — SN_A basis (H mixed into outer ejecta surface layers,
             compact ~1 Rsun progenitor, mdot=1e-6 wind, hmax=500k s)
    SN_A_S_E07   E_burst 9 -> 7          (cooler photosphere -> redder)
    SN_A_S_Ni05  Ni56 0.075 -> 0.05      (faster post-peak decay)
    SN_A_S_H15   h_surface 0.10 -> 0.15  (more H opacity -> redder)

  SN_A_X_* — SN_A_v4 basis (detached 0.10 Msun H envelope, 1.2->200 Rsun,
             rho~r^-3.5, T 50kK->6kK, mdot=1e-9 wind, hmax=3600 s)
    SN_A_X_E07   E_burst 9 -> 7          (cooler photosphere -> redder)
    SN_A_X_Ni05  Ni56 0.075 -> 0.05      (faster post-peak decay)
    SN_A_X_R300  env radius 200 -> 300   (larger R_ph -> redder + more SBO)

All: parallel OpenStella, cpus_per_task=3, STELLA-only, 65 d.
Naming: SN_A_{S,X}_{param varied}.  Directories at top level so
plot_sn_a_all_models.py picks them up via the SN_A* glob.
"""
import os
import copy
import yaml
from pathlib import Path

BASE_DIR = Path("/viper/ptmp1/sbrennan/SN2024pba")
PARALLEL_STELLA = "/ptmp/sbrennan/stella/OpenStella_parallel"
PRESN_SOURCE = BASE_DIR / "SN_A" / "s12.0-presn-stripped"

# (name, template_dir, description of change, edits dict)
# edits: {("stella","explosion_energy"): value, ...}
MODELS = [
    # --- S family (SN_A basis) ---
    ("SN_A_S_E07", "SN_A", "E_burst 9->7 (redder)",
     {("stella", "explosion_energy"): 7.0}),
    ("SN_A_S_Ni05", "SN_A", "Ni56 0.075->0.05 (faster decay)",
     {("stella", "nickel56", "total_mass_msun"): 0.05}),
    ("SN_A_S_H15", "SN_A", "h_surface_mass 0.10->0.15 Msun (redder)",
     {("kepler_presn", "h_surface_mass_msun"): 0.15}),
    # --- X family (SN_A_v4 basis) ---
    ("SN_A_X_E07", "SN_A_v4", "E_burst 9->7 (redder)",
     {("stella", "explosion_energy"): 7.0}),
    ("SN_A_X_Ni05", "SN_A_v4", "Ni56 0.075->0.05 (faster decay)",
     {("stella", "nickel56", "total_mass_msun"): 0.05}),
    ("SN_A_X_R300", "SN_A_v4", "h_ext_r_outer 200->300 Rsun (redder+SBO)",
     {("kepler_presn", "h_ext_r_outer_rsun"): 300.0}),
]


def set_nested(cfg, path, value):
    d = cfg
    for k in path[:-1]:
        d = d[k]
    d[path[-1]] = value


def make_config(base, name, change_desc, edits):
    cfg = copy.deepcopy(base)
    family = "surface-layer (compact)" if "_S_" in name else \
             "extended-envelope (200 Rsun)"
    cfg["model_name"] = name
    cfg["description"] = (
        f"SN_A variant {name}: {family} H distribution. "
        f"Change vs basis: {change_desc}. "
        f"Target: SN 2024pba first-peak g/r/i lightcurve. "
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

    for path, val in edits.items():
        set_nested(cfg, path, val)

    cfg["stella"]["stella_root"] = PARALLEL_STELLA
    cfg["stella"]["cpus_per_task"] = 3

    cfg["kepler_presn"]["presn_file"] = (
        f"/ptmp/sbrennan/SN2024pba/{name}/s12.0-presn-stripped")
    cfg["output"]["model_dir"] = f"./{name}"
    return cfg


def main():
    print(f"Generating {len(MODELS)} SN_A variant configs\n")
    print(f"{'Model':<14} {'Basis':<9} {'Change'}")
    print("-" * 60)

    configs = []
    for name, template_dir, change_desc, edits in MODELS:
        with open(BASE_DIR / template_dir / "config.yaml") as f:
            base = yaml.safe_load(f)

        model_dir = BASE_DIR / name
        model_dir.mkdir(parents=True, exist_ok=True)
        kep_link = model_dir / "s12.0-presn-stripped"
        if not kep_link.exists():
            os.symlink(PRESN_SOURCE, kep_link)

        cfg = make_config(base, name, change_desc, edits)
        config_path = model_dir / "config.yaml"
        header = (
            f"# SN_A variant: {name}\n"
            f"# Basis: {template_dir} "
            f"({'surface H layer' if '_S_' in name else 'extended H env @200Rsun'})\n"
            f"# Change: {change_desc}\n"
            f"# Target: SN 2024pba first-peak g/r/i lightcurve\n"
            f"# Codes: STELLA only (parallel OpenStella build)\n\n"
        )
        with open(config_path, "w") as f:
            f.write(header)
            yaml.dump(cfg, f, default_flow_style=False, sort_keys=False,
                      width=120)
        configs.append(str(config_path))
        print(f"{name:<14} {template_dir:<9} {change_desc}")

    print(f"\nGenerated {len(configs)} configs.")


if __name__ == "__main__":
    main()
