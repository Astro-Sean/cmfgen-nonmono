#!/usr/bin/env python3
"""Generate the SN_B_v2 T6k wind grid: higher wind Mdot to match the
observed SN 2024pba SN_B lightcurve peak (M_r = -18.2) and slow decay
(~2.5 mag/100d at +50d, ~1.9 mag/100d at +100d).

The original T6k grid used Mdot=1e-3 Msun/yr and peaked at M_R ~ -16.2
with a decay of 5-8 mag/100d — 2 mag too faint and 2-3x too fast.

Physics:
  L_peak ~ Mdot * v_shell^3 / v_wind   (CSM interaction luminosity)
  t_decay ~ Mdot / (v_wind * v_shell)  (CSM crossing time)

To get ~6x brighter (2 mag) AND ~2.5x slower decay, increase Mdot by
3-10x while keeping moderate shell velocities.

Grid (12 models):
  Wind Mdot:      3e-3, 5e-3, 1e-2 Msun/yr  (3x, 5x, 10x original)
  Shell mass:     3.0, 5.0 Msun
  Shell velocity:  5000, 6000 km/s
  Wind velocity:  100 km/s (unchanged)
  Ni56:           1e-4 Msun (unchanged)
  Source:         SN_A_g_T20k_Ni06 (400-day-old SN_A ejecta)
  STELLA:         120 days (covers +100d post peak)

Codes: STELLA only (no CMFGEN/SEDONA/SuperLite)

Each config is the parent configs/SN_B_v2_M3_V5000.yaml with only the
wind Mdot, shell mass, shell velocity, and target_days modified.
"""
import yaml
from pathlib import Path

BASE_DIR = Path("/viper/ptmp1/sbrennan/SN2024pba")
GRID_DIR = BASE_DIR / "SN_B_v2_T6k_grid"
PARENT_CONFIG = BASE_DIR / "configs" / "SN_B_v2_M3_V5000.yaml"

# Grid parameters
MDOTS = [3e-3, 5e-3, 1e-2]        # Msun/yr
MASSES = [3.0, 5.0]               # Msun
VELOCITIES = [5000, 6000]         # km/s
TARGET_DAYS = 120.0               # covers +100d post peak

GRID_POINTS = [
    (m, v, mdot, TARGET_DAYS)
    for mdot in MDOTS
    for m in MASSES
    for v in VELOCITIES
]


def _p(x, nd=1):
    """7.0 -> '7p0', 2.5 -> '2p5' (decimal point -> p)."""
    return f"{x:.{nd}f}".replace(".", "p")


def _e(x):
    """5000 -> '5e3', 0.003 -> '3e-3' (compact e-notation)."""
    s = f"{x:.0e}"
    return s.replace("e+0", "e").replace("e+", "e").replace("e-0", "e-")


def grid_name(mass_msun, vel_kms, mdot_msun_yr=None):
    """X_M3p0_V5e3_W3e-3 : shell mass, velocity, wind Mdot (literal values)."""
    name = f"X_M{_p(mass_msun)}_V{_e(vel_kms)}"
    if mdot_msun_yr is not None and mdot_msun_yr != 1e-3:
        name += f"_W{_e(mdot_msun_yr)}"
    return name


def main():
    with open(PARENT_CONFIG) as f:
        parent = yaml.safe_load(f)

    GRID_DIR.mkdir(parents=True, exist_ok=True)
    configs = []
    for mass, vel, mdot, target_days in GRID_POINTS:
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
            f"T6k wind grid: H/He shell M={mass} Msun, v={vel} km/s, "
            f"solar comp; wind Mdot={mdot:.0e} Msun/yr, v=100 km/s; "
            f"400d-old SN_A_g_T20k_Ni06 ejecta; Ni56=1e-4 Msun. "
            f"STELLA to {target_days:.0f} d. Target: M_peak ~ -18.2, "
            f"slow decay ~2.5 mag/100d."
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
        x["wind_mdot_msun_yr"] = mdot
        cfg["x_model_from_sn"] = x

        st = dict(parent["stella"])
        st["target_days"] = target_days
        cfg["stella"] = st

        out = dict(parent["output"])
        out["model_dir"] = f"./SN_B_v2_T6k_grid/{model_name}"
        cfg["output"] = out

        header = (
            f"# SN_B_v2 T6k wind grid model: {model_name}\n"
            f"# Higher wind Mdot to match observed SN_B peak M_r=-18.2 "
            f"and slow decay ~2.5 mag/100d\n"
            f"# Shell: M={mass} Msun, v={vel} km/s, solar composition\n"
            f"# Wind: Mdot={mdot:.0e} Msun/yr, v=100 km/s; Ni56=1e-4 Msun\n"
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
