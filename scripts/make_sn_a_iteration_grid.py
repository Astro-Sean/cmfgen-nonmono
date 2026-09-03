#!/usr/bin/env python3
"""Generate SN_A iteration grid configs for SN 2011dh matching.

Based on the best grid model (SN_A_g_T20k_Ni06), this script creates 8
focused iterations targeting the remaining discrepancies:

1. SBO flash too bright at day 7-10 (1-2 mag excess bolometric)
2. U-V colour too blue (0.5-1.3 mag)
3. V peak slightly too faint and early (17.7d vs 20d, -17.07 vs -17.17)

Key physics knobs:
  - h_ext_mass_msun:  0.05 -> 0.10  (more H envelope dims SBO, reddens colours)
  - h_ext_r_outer_rsun: 200 -> 300  (larger radius dims/broadens SBO)
  - h_ext_T_base_K:   20k -> 15k    (cooler envelope reddens U-V)
  - kappa:            0.34 -> 0.40  (higher opacity broadens peak, dims SBO)
  - Ni56:             0.06 -> 0.07  (slightly brighter peak)

The original SN_A config used M_env=0.10 Msun (specifically for 2011dh),
but the grid never tested this value — all grid models used 0.05 Msun.

Models are placed in SN_A/ with descriptive names.
"""
import os, yaml, copy

BASE_DIR = "/viper/ptmp1/sbrennan/SN2024pba"
TEMPLATE = f"{BASE_DIR}/SN_A/SN_A_g_T20k_Ni06/config.yaml"
PRESN_SOURCE = f"{BASE_DIR}/SN_A/s12.0-presn-stripped"

with open(TEMPLATE) as f:
    base = yaml.safe_load(f)

# ── Iteration grid ──────────────────────────────────────────────────────────
# Each entry: (name, description, parameter overrides)
# Names are descriptive of what changed from the baseline.
models = [
    ("baseline",
     "Baseline: best grid model (M_env=0.05, R=200, T=20k, kappa=0.34, Ni=0.06, E=9.0)",
     {}),
    ("Menv010",
     "More H envelope: M_env=0.10 (original SN_A value, never tested in grid)",
     {"h_ext_mass_msun": 0.10}),
    ("R300",
     "Larger radius: R=300 Rsun (SN 2011dh YSG progenitor range)",
     {"h_ext_r_outer_rsun": 300}),
    ("Menv010_R300",
     "More H + larger radius: M_env=0.10, R=300 (combined SBO dimming)",
     {"h_ext_mass_msun": 0.10, "h_ext_r_outer_rsun": 300}),
    ("T15k",
     "Cooler envelope: T_base=15kK (redder U-V, less UV excess)",
     {"h_ext_T_base_K": 15000}),
    ("Menv010_R300_T15k",
     "Best-guess combo: M_env=0.10, R=300, T_base=15k (all SBO+colour fixes)",
     {"h_ext_mass_msun": 0.10, "h_ext_r_outer_rsun": 300, "h_ext_T_base_K": 15000}),
    ("Menv010_R300_T15k_k40",
     "Combo + higher kappa: M_env=0.10, R=300, T=15k, kappa=0.40 (broader peak)",
     {"h_ext_mass_msun": 0.10, "h_ext_r_outer_rsun": 300, "h_ext_T_base_K": 15000,
      "kappa": 0.40}),
    ("Menv010_R300_T15k_k40_Ni07",
     "Combo + higher kappa + more Ni: M_env=0.10, R=300, T=15k, kappa=0.40, Ni=0.07",
     {"h_ext_mass_msun": 0.10, "h_ext_r_outer_rsun": 300, "h_ext_T_base_K": 15000,
      "kappa": 0.40, "ni56": 0.07}),
]

print(f"Generating {len(models)} iteration models in SN_A/\n")

for name, desc, overrides in models:
    model_dir = f"{BASE_DIR}/SN_A/{name}"
    os.makedirs(model_dir, exist_ok=True)

    # Symlink KEPLER pre-SN source
    kep_link = f"{model_dir}/s12.0-presn-stripped"
    if not os.path.exists(kep_link):
        os.symlink(PRESN_SOURCE, kep_link)

    cfg = copy.deepcopy(base)
    cfg["model_name"] = name
    cfg["description"] = desc
    cfg["run_cmfgen"] = False
    cfg["run_sedona"] = False
    cfg["run_sirocco"] = False
    cfg["run_superlite"] = False

    # Apply overrides
    if "h_ext_mass_msun" in overrides:
        cfg["kepler_presn"]["h_ext_mass_msun"] = overrides["h_ext_mass_msun"]
    if "h_ext_r_outer_rsun" in overrides:
        cfg["kepler_presn"]["h_ext_r_outer_rsun"] = overrides["h_ext_r_outer_rsun"]
    if "h_ext_T_base_K" in overrides:
        cfg["kepler_presn"]["h_ext_T_base_K"] = overrides["h_ext_T_base_K"]
    if "kappa" in overrides:
        cfg["stella"]["kappa"] = overrides["kappa"]
    if "ni56" in overrides:
        cfg["stella"]["nickel56"]["total_mass_msun"] = overrides["ni56"]

    # Update paths for nested directory
    cfg["output"]["model_dir"] = f"./SN_A/{name}"
    cfg["kepler_presn"]["presn_file"] = (
        f"/ptmp/sbrennan/SN2024pba/SN_A/{name}/s12.0-presn-stripped")

    config_path = f"{model_dir}/config.yaml"
    with open(config_path, "w") as f:
        yaml.dump(cfg, f, default_flow_style=False, sort_keys=False, width=120)

    # Print summary
    m_env = cfg["kepler_presn"]["h_ext_mass_msun"]
    r_env = cfg["kepler_presn"]["h_ext_r_outer_rsun"]
    t_base = cfg["kepler_presn"]["h_ext_T_base_K"]
    kappa = cfg["stella"]["kappa"]
    ni56 = cfg["stella"]["nickel56"]["total_mass_msun"]
    e_burst = cfg["stella"]["explosion_energy"]
    print(f"  {name:<30}: M_env={m_env}, R={r_env}, T={t_base}K, "
          f"kappa={kappa}, Ni56={ni56}, E={e_burst}")

print(f"\nGenerated {len(models)} configs in SN_A/")
