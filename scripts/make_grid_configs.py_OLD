#!/usr/bin/env python3
"""Generate SN_A grid configs for SN 2011dh matching.

Grid explores:
- Envelope mass: 0.02, 0.05, 0.10 Msun
- Envelope radius: 50, 100, 200 Rsun
- T_base: 20000, 35000, 50000 K
- Explosion energy: 6.0, 9.0, 12.0
- n_rho: 2.0, 3.5, 6.0
- H fraction: 0.4, 0.7

Each model gets a unique name encoding its key parameters.
"""
import os, yaml, copy

BASE_DIR = "/viper/ptmp1/sbrennan/SN2024pba"
TEMPLATE = "/viper/ptmp1/sbrennan/SN2024pba/SN_A_v4/config.yaml"

with open(TEMPLATE) as f:
    base = yaml.safe_load(f)

# Grid definitions: (name_suffix, overrides)
grid = [
    # --- Group 1: Envelope mass variation (fix R=200, T=50k, E=9, n_rho=3.5, H=0.7) ---
    ("M02",  {"h_ext_mass_msun": 0.02, "h_ext_r_outer_rsun": 200, "h_ext_T_base_K": 50000, "explosion_energy": 9.0, "h_ext_n_rho": 3.5, "h_ext_fraction": 0.7}),
    ("M05",  {"h_ext_mass_msun": 0.05, "h_ext_r_outer_rsun": 200, "h_ext_T_base_K": 50000, "explosion_energy": 9.0, "h_ext_n_rho": 3.5, "h_ext_fraction": 0.7}),
    # SN_A_v4 is M10 (already done)

    # --- Group 2: Envelope radius variation (fix M=0.05, T=50k, E=9, n_rho=3.5, H=0.7) ---
    ("R050", {"h_ext_mass_msun": 0.05, "h_ext_r_outer_rsun": 50,  "h_ext_T_base_K": 50000, "explosion_energy": 9.0, "h_ext_n_rho": 3.5, "h_ext_fraction": 0.7}),
    ("R100", {"h_ext_mass_msun": 0.05, "h_ext_r_outer_rsun": 100, "h_ext_T_base_K": 50000, "explosion_energy": 9.0, "h_ext_n_rho": 3.5, "h_ext_fraction": 0.7}),
    # R200 = M05 above

    # --- Group 3: T_base variation (fix M=0.05, R=200, E=9, n_rho=3.5, H=0.7) ---
    ("T20k", {"h_ext_mass_msun": 0.05, "h_ext_r_outer_rsun": 200, "h_ext_T_base_K": 20000, "explosion_energy": 9.0, "h_ext_n_rho": 3.5, "h_ext_fraction": 0.7}),
    ("T35k", {"h_ext_mass_msun": 0.05, "h_ext_r_outer_rsun": 200, "h_ext_T_base_K": 35000, "explosion_energy": 9.0, "h_ext_n_rho": 3.5, "h_ext_fraction": 0.7}),
    # T50k = M05 above

    # --- Group 4: Explosion energy variation (fix M=0.05, R=200, T=50k, n_rho=3.5, H=0.7) ---
    ("E06",  {"h_ext_mass_msun": 0.05, "h_ext_r_outer_rsun": 200, "h_ext_T_base_K": 50000, "explosion_energy": 6.0, "h_ext_n_rho": 3.5, "h_ext_fraction": 0.7}),
    ("E12",  {"h_ext_mass_msun": 0.05, "h_ext_r_outer_rsun": 200, "h_ext_T_base_K": 50000, "explosion_energy": 12.0, "h_ext_n_rho": 3.5, "h_ext_fraction": 0.7}),
    # E09 = M05 above

    # --- Group 5: Density profile variation (fix M=0.05, R=200, T=50k, E=9, H=0.7) ---
    ("N02",  {"h_ext_mass_msun": 0.05, "h_ext_r_outer_rsun": 200, "h_ext_T_base_K": 50000, "explosion_energy": 9.0, "h_ext_n_rho": 2.0, "h_ext_fraction": 0.7}),
    ("N06",  {"h_ext_mass_msun": 0.05, "h_ext_r_outer_rsun": 200, "h_ext_T_base_K": 50000, "explosion_energy": 9.0, "h_ext_n_rho": 6.0, "h_ext_fraction": 0.7}),
    # N35 = M05 above

    # --- Group 6: H fraction variation (fix M=0.05, R=200, T=50k, E=9, n_rho=3.5) ---
    ("H04",  {"h_ext_mass_msun": 0.05, "h_ext_r_outer_rsun": 200, "h_ext_T_base_K": 50000, "explosion_energy": 9.0, "h_ext_n_rho": 3.5, "h_ext_fraction": 0.4}),
    # H07 = M05 above

    # --- Group 7: Best-guess combinations ---
    # Ergon's lower H estimate (0.02-0.05) + smaller radius + lower T
    ("Ergon_lo", {"h_ext_mass_msun": 0.03, "h_ext_r_outer_rsun": 150, "h_ext_T_base_K": 35000, "explosion_energy": 8.0, "h_ext_n_rho": 3.5, "h_ext_fraction": 0.6}),
    # Compact + low energy to minimize SBO
    ("Compact",  {"h_ext_mass_msun": 0.03, "h_ext_r_outer_rsun": 100, "h_ext_T_base_K": 30000, "explosion_energy": 7.0, "h_ext_n_rho": 4.0, "h_ext_fraction": 0.5}),
    # Higher energy + low mass for faster rise
    ("E12_M03",  {"h_ext_mass_msun": 0.03, "h_ext_r_outer_rsun": 200, "h_ext_T_base_K": 40000, "explosion_energy": 12.0, "h_ext_n_rho": 3.5, "h_ext_fraction": 0.7}),
]

for suffix, overrides in grid:
    model_name = f"SN_A_g_{suffix}"
    model_dir = f"{BASE_DIR}/{model_name}"
    os.makedirs(model_dir, exist_ok=True)

    # Symlink KEPLER source
    kep_link = f"{model_dir}/s12.0-presn-stripped"
    if not os.path.exists(kep_link):
        os.symlink("/ptmp/sbrennan/SN2024pba/SN_A/s12.0-presn-stripped", kep_link)

    cfg = copy.deepcopy(base)
    cfg["model_name"] = model_name
    cfg["run_cmfgen"] = False  # STELLA only
    cfg["run_sedona"] = False
    cfg["run_sirocco"] = False
    cfg["run_superlite"] = False

    # Apply overrides
    kp = cfg["kepler_presn"]
    kp["h_ext_mass_msun"] = overrides["h_ext_mass_msun"]
    kp["h_ext_r_outer_rsun"] = overrides["h_ext_r_outer_rsun"]
    kp["h_ext_T_base_K"] = overrides["h_ext_T_base_K"]
    kp["h_ext_fraction"] = overrides["h_ext_fraction"]
    kp["h_ext_n_rho"] = overrides["h_ext_n_rho"]
    cfg["stella"]["explosion_energy"] = overrides["explosion_energy"]

    # Update paths
    cfg["output"]["model_dir"] = f"./{model_name}"
    cfg["kepler_presn"]["presn_file"] = f"/ptmp/sbrennan/SN2024pba/{model_name}/s12.0-presn-stripped"

    # Description
    m = overrides["h_ext_mass_msun"]
    r = overrides["h_ext_r_outer_rsun"]
    t = overrides["h_ext_T_base_K"]
    e = overrides["explosion_energy"]
    n = overrides["h_ext_n_rho"]
    h = overrides["h_ext_fraction"]
    cfg["description"] = f"Grid model {suffix}: M_env={m} Msun, R={r} Rsun, T_base={t} K, EBurst={e}, n_rho={n}, X_H={h}"

    config_path = f"{model_dir}/config.yaml"
    with open(config_path, "w") as f:
        yaml.dump(cfg, f, default_flow_style=False, sort_keys=False, width=120)
    print(f"  {model_name}: M={m}, R={r}, T={t}, E={e}, n={n}, H={h}")

print(f"\nGenerated {len(grid)} grid configs")
