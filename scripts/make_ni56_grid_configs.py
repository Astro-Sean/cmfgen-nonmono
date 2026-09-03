#!/usr/bin/env python3
"""Generate SN_A Ni56 grid configs for SN 2011dh matching.

Based on the best-matching model (SN_A_g_T20k: M_env=0.05, R=200, T_base=20k,
E=9.0, n_rho=3.5, X_H=0.7), this grid varies only the Ni56 mass to chase down
the ~0.3 mag over-brightness in V/B at peak.

Ni56 masses: 0.040, 0.050, 0.060, 0.075 (baseline), 0.090

The 2011dh observed V peak is -16.90 at 19.5d; the T20k model peaks at -17.26
(0.36 mag too bright). Arnett's rule suggests L_peak ~ M_Ni, so reducing Ni56
by ~30% (0.075 -> 0.05) should bring the peak down by ~0.3-0.4 mag.
"""
import os, yaml, copy

BASE_DIR = "/viper/ptmp1/sbrennan/SN2024pba"
TEMPLATE = "/viper/ptmp1/sbrennan/SN2024pba/SN_A_g_T20k/config.yaml"

with open(TEMPLATE) as f:
    base = yaml.safe_load(f)

# Ni56 grid: (name_suffix, total_mass_msun)
ni56_grid = [
    ("Ni04", 0.040),
    ("Ni05", 0.050),
    ("Ni06", 0.060),
    ("Ni075", 0.075),  # baseline (same as T20k)
    ("Ni09", 0.090),
]

for suffix, ni56_mass in ni56_grid:
    model_name = f"SN_A_g_T20k_{suffix}"
    model_dir = f"{BASE_DIR}/{model_name}"
    os.makedirs(model_dir, exist_ok=True)

    # Symlink KEPLER source (same as T20k)
    kep_link = f"{model_dir}/s12.0-presn-stripped"
    if not os.path.exists(kep_link):
        os.symlink("/ptmp/sbrennan/SN2024pba/SN_A/s12.0-presn-stripped", kep_link)

    cfg = copy.deepcopy(base)
    cfg["model_name"] = model_name
    cfg["run_cmfgen"] = False
    cfg["run_sedona"] = False
    cfg["run_sirocco"] = False
    cfg["run_superlite"] = False

    # Override Ni56 mass
    cfg["stella"]["nickel56"]["total_mass_msun"] = ni56_mass

    # Update paths
    cfg["output"]["model_dir"] = f"./{model_name}"
    cfg["kepler_presn"]["presn_file"] = f"/ptmp/sbrennan/SN2024pba/{model_name}/s12.0-presn-stripped"

    # Description
    cfg["description"] = (
        f"Ni56 grid model {suffix}: based on T20k (M_env=0.05, R=200, T_base=20k, "
        f"E=9.0, n_rho=3.5, X_H=0.7) with Ni56={ni56_mass} Msun"
    )

    config_path = f"{model_dir}/config.yaml"
    with open(config_path, "w") as f:
        yaml.dump(cfg, f, default_flow_style=False, sort_keys=False, width=120)
    print(f"  {model_name}: Ni56={ni56_mass} Msun")

print(f"\nGenerated {len(ni56_grid)} Ni56 grid configs")
print("\nTo run all:")
for suffix, _ in ni56_grid:
    print(f"  python scripts/run_model_grid.py SN_A_g_T20k_{suffix}/config.yaml --skip-stella")
print("\n(Note: --skip-stella not needed for first run; use it only to rerun CMFGEN/SEDONA)")
print("\nTo run STELLA only (no CMFGEN/SEDONA):")
for suffix, _ in ni56_grid:
    print(f"  python scripts/run_model_grid.py SN_A_g_T20k_{suffix}/config.yaml")
