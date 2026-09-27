#!/usr/bin/env python3
"""Generate the SN_B_v2 T6k grid: (shell mass, velocity) variants of
SN_B_v2_M3_V5000 targeting Tbb ~ 6000 K at 30-40 d post-max.

Grid: 3x3 base + 2x2 high end extension
  Shell mass:     2.0, 3.0, 5.0 Msun           (base grid, 55 d)
                  7.0, 10.0 Msun               (extension, 100 d)
  Shell velocity: 4000, 5000, 6000 km/s        (base grid, 55 d)
                  8000, 10000 km/s             (extension, 100 d)
  Wind:           Mdot=1e-3 Msun/yr, v=100 km/s (as SN_B_v2_M3_V5000)
  Ni56:           1e-4 Msun (as parent)
  Source:         SN_A_g_T20k_Ni06 (400-day-old SN_A ejecta)
  STELLA:         55 days base (covers 30-40 d post-max; parent peaks ~8.3 d);
                  100 days extension (heavier/faster shells peak later)
Codes: STELLA only (no CMFGEN/SEDONA/SuperLite)

Each config is the parent configs/SN_B_v2_M3_V5000.yaml with only the
keys listed above modified — all other physics/numerics are identical.
Existing config.yaml files are not overwritten; re-running only writes
missing grid points.
"""
import yaml
from pathlib import Path

BASE_DIR = Path("/viper/ptmp1/sbrennan/SN2024pba")
GRID_DIR = BASE_DIR / "SN_B_v2_T6k_grid"
PARENT_CONFIG = BASE_DIR / "configs" / "SN_B_v2_M3_V5000.yaml"

MASSES = [2.0, 3.0, 5.0]
VELOCITIES = [4000, 5000, 6000]
TARGET_DAYS = 55.0

# High-end extension: heavier/faster shells, longer evolution.
EXT_MASSES = [7.0, 10.0]
EXT_VELOCITIES = [8000, 10000]
EXT_TARGET_DAYS = 100.0

GRID_POINTS = (
    [(m, v, TARGET_DAYS) for m in MASSES for v in VELOCITIES]
    + [(m, v, EXT_TARGET_DAYS) for m in EXT_MASSES for v in EXT_VELOCITIES]
)


def main():
    with open(PARENT_CONFIG) as f:
        parent = yaml.safe_load(f)

    GRID_DIR.mkdir(parents=True, exist_ok=True)
    configs = []
    for mass, vel, target_days in GRID_POINTS:
        model_name = f"T6k_M{mass:.1f}_V{vel}"
        model_dir = GRID_DIR / model_name
        config_path = model_dir / "config.yaml"
        if config_path.exists():
            print(f"  Exists, skipping: {config_path}")
            continue
        model_dir.mkdir(parents=True, exist_ok=True)

        cfg = dict(parent)
        cfg["model_name"] = model_name
        cfg["description"] = (
            f"T6k grid point of SN_B_v2_M3_V5000: H/He shell "
            f"M={mass} Msun, v={vel} km/s, solar comp; wind "
            f"Mdot=1e-3 Msun/yr, v=100 km/s; 400d-old "
            f"SN_A_g_T20k_Ni06 ejecta; Ni56=1e-4 Msun. "
            f"STELLA to {target_days:.0f} d. Target: Tbb ~ 6000 K "
            f"at 30-40 d post-max."
        )

        # STELLA only
        cfg["run_cmfgen"] = False
        cfg["run_sirocco"] = False
        cfg["run_sedona"] = False
        cfg["run_superlite"] = False
        cfg["skip_stella"] = False
        cfg["reuse_stella"] = False
        cfg["submission"] = "normal"
        cfg["debug"] = False

        x = dict(parent["x_model_from_sn"])
        x["shell_mass_msun"] = mass
        x["shell_velocity_kms"] = vel
        cfg["x_model_from_sn"] = x

        st = dict(parent["stella"])
        st["target_days"] = target_days
        cfg["stella"] = st

        out = dict(parent["output"])
        out["model_dir"] = f"./SN_B_v2_T6k_grid/{model_name}"
        cfg["output"] = out

        header = (
            f"# SN_B_v2 T6k grid model: {model_name}\n"
            f"# Variant of SN_B_v2_M3_V5000 targeting Tbb ~ 6000 K at 30-40 d post-max\n"
            f"# Shell: M={mass} Msun, v={vel} km/s, solar composition\n"
            f"# Wind: Mdot=1e-3 Msun/yr, v=100 km/s; Ni56=1e-4 Msun\n"
            f"# Source: SN_A_g_T20k_Ni06 (400-day-old SN_A ejecta)\n"
            f"# STELLA evolution: {target_days:.0f} days; STELLA only\n\n"
        )
        with open(config_path, "w") as f:
            f.write(header)
            yaml.dump(cfg, f, default_flow_style=False, sort_keys=False)
        configs.append(str(config_path))
        print(f"  Written: {config_path}")

    print(f"\nGenerated {len(configs)} configs in {GRID_DIR}/")
    print("\nSubmit all with:")
    for c in configs:
        print(f"  python scripts/run_model_grid.py {c}")


if __name__ == "__main__":
    main()
