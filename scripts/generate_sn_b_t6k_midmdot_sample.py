#!/usr/bin/env python3
"""Generate the intermediate-Mdot SN_B_v2 T6k sample (parallel STELLA).

Motivation: the base T6k models (Mdot=1e-3) reproduce the observed
double-peak morphology but decay too fast post-peak (~6-10 mag/100d vs
observed ~2).  The high-Mdot X_* models (1e-2..3e-2) match the peak and
slow decay but smear the flash+dip+main-peak into a single broad maximum
because the wind is optically thick from the start (tau ~ 3-8).

This sample probes the intermediate regime (Mdot = 3e-3, 5e-3; tau ~ 1-2)
at the bright family (M=7,10 Msun, v=10000 km/s) to find where the dip
survives while the decay slows.

Models run with the OpenMP-parallel STELLA build
(/ptmp/sbrennan/stella/OpenStella_parallel, cpus_per_task=3), matching
the validated parallel_test_p1/p2/p3 runs.

Sample (4 models): the two brightest no-wind models, each with
  Wind Mdot:      3e-3 (W3e-3), 5e-3 (W5e-3) Msun/yr  (3x, 5x base 1e-3)
  Shell mass:     7.0, 10.0 Msun
  Shell velocity: 10000 km/s
  Wind velocity:  100 km/s (unchanged)
  Ni56:           1e-4 Msun (unchanged)
  Source:         SN_A_g_T20k_Ni06 (400-day-old SN_A ejecta)
  STELLA:         150 days (same as the X_* high-Mdot family)

Each config is the corresponding T6k model config with only the wind
Mdot, target_days, stella_root and cpus_per_task modified.
"""
import yaml
from pathlib import Path

BASE_DIR = Path("/viper/ptmp1/sbrennan/SN2024pba")
GRID_DIR = BASE_DIR / "SN_B_v2_T6k_grid"
PARALLEL_STELLA = "/ptmp/sbrennan/stella/OpenStella_parallel"

# (parent model name, mdot) pairs
PARENTS = ["T6k_M7.0_V10000", "T6k_M10.0_V10000"]
MDOTS = [3e-3, 5e-3]          # Msun/yr
TARGET_DAYS = 150.0
CPUS_PER_TASK = 3
# SN_A ejecta age at SN_B STELLA start: SN_A launch MJD 60490,
# STELLA start MJD 60822 -> 332 days (supersedes the sa320 estimate).
SOURCE_AGE_DAYS = 332.0


def _p(x, nd=1):
    """7.0 -> '7p0', 2.5 -> '2p5', 10.0 -> '10p0' (decimal point -> p)."""
    return f"{x:.{nd}f}".replace(".", "p")


def _e(x):
    """10000 -> '1e4', 0.003 -> '3e-3' (compact e-notation)."""
    s = f"{x:.0e}"
    return s.replace("e+0", "e").replace("e+", "e").replace("e-0", "e-")


def grid_name(mass_msun, vel_kms, mdot_msun_yr=None):
    """X_M7p0_V1e4_W3e-3 : shell mass (Msun), velocity (km/s), wind Mdot."""
    name = f"X_M{_p(mass_msun)}_V{_e(vel_kms)}"
    if mdot_msun_yr is not None and mdot_msun_yr != 1e-3:
        name += f"_W{_e(mdot_msun_yr)}"
    return name


def main():
    GRID_DIR.mkdir(parents=True, exist_ok=True)
    configs = []
    for parent_name in PARENTS:
        parent_path = GRID_DIR / parent_name / "config.yaml"
        with open(parent_path) as f:
            parent = yaml.safe_load(f)
        mass = parent["x_model_from_sn"]["shell_mass_msun"]
        vel = parent["x_model_from_sn"]["shell_velocity_kms"]

        for mdot in MDOTS:
            model_name = grid_name(mass, vel, mdot)
            model_dir = GRID_DIR / model_name
            config_path = model_dir / "config.yaml"
            if config_path.exists():
                print(f"  Exists, skipping: {config_path}")
                continue
            model_dir.mkdir(parents=True, exist_ok=True)

            cfg = dict(parent)
            cfg["model_name"] = model_name
            cfg["description"] = (
                f"T6k intermediate-Mdot sample (parallel STELLA): H/He "
                f"shell M={mass} Msun, v={vel} km/s, solar comp; wind "
                f"Mdot={mdot:.0e} Msun/yr, v=100 km/s; 332d-old "
                f"SN_A_g_T20k_Ni06 ejecta; Ni56=1e-4 Msun. STELLA to "
                f"{TARGET_DAYS:.0f} d. Tests whether the flash+dip "
                f"survives at tau_wind ~ 1-2 while decay slows."
            )

            cfg["run_cmfgen"] = False
            cfg["run_sirocco"] = False
            cfg["run_sedona"] = False
            cfg["run_superlite"] = False
            cfg["skip_stella"] = False
            cfg["reuse_stella"] = False
            cfg["submission"] = "normal"
            cfg["debug"] = False

            x = dict(parent["x_model_from_sn"])
            x["wind_mdot_msun_yr"] = mdot
            x["source_age_days"] = SOURCE_AGE_DAYS
            cfg["x_model_from_sn"] = x

            st = dict(parent["stella"])
            st["target_days"] = TARGET_DAYS
            st["stella_root"] = PARALLEL_STELLA
            st["cpus_per_task"] = CPUS_PER_TASK
            cfg["stella"] = st

            out = dict(parent["output"])
            out["model_dir"] = f"./SN_B_v2_T6k_grid/{model_name}"
            cfg["output"] = out

            header = (
                f"# SN_B_v2 T6k intermediate-Mdot sample: {model_name}\n"
                f"# Mdot between thin-wind (1e-3, double peak, fast decay)\n"
                f"# and thick-wind (>=1e-2, single broad peak) regimes\n"
                f"# Shell: M={mass} Msun, v={vel} km/s, solar composition\n"
                f"# Wind: Mdot={mdot:.0e} Msun/yr, v=100 km/s; "
                f"Ni56=1e-4 Msun\n"
                f"# Source: SN_A_g_T20k_Ni06 (332-day-old SN_A ejecta; "
                f"SN_A MJD 60490 -> STELLA start MJD 60822)\n"
                f"# STELLA evolution: {TARGET_DAYS:.0f} days; STELLA only\n"
                f"# Parallel STELLA: {PARALLEL_STELLA} "
                f"(cpus_per_task={CPUS_PER_TASK})\n\n"
            )
            with open(config_path, "w") as f:
                f.write(header)
                yaml.dump(cfg, f, default_flow_style=False, sort_keys=False)
            configs.append(str(config_path))
            print(f"  Written: {config_path}")

    print(f"\nGenerated {len(configs)} configs in {GRID_DIR}/")
    print("\nSubmit all with:")
    print(f"  python scripts/run_model_grid.py {' '.join(configs)}")


if __name__ == "__main__":
    main()
