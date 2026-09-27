#!/usr/bin/env python3
"""Generate a small high-Mdot wind sample on the peak-matching T6k models.

Motivation: the bright T6k models (M7-10, V8000-10000) match the observed
SN_B peak (M_o ~ -18.6, M_r ~ -18.2) and the double-peak morphology, but
decay post-peak at ~6-10 mag/100d vs the observed ~2 mag/100d.  A higher
wind mass-loss rate increases the CSM column and diffusion/interaction
timescale, which should delay the decline.

Sample (4 models): the two brightest no-wind models, each with
  Wind Mdot:      1e-2 (W10), 3e-2 (W30) Msun/yr  (10x, 30x original 1e-3)
  Shell mass:     7.0, 10.0 Msun
  Shell velocity: 10000 km/s
  Wind velocity:  100 km/s (unchanged)
  Ni56:           1e-4 Msun (unchanged)
  Source:         SN_A_g_T20k_Ni06 (400-day-old SN_A ejecta)
  STELLA:         150 days (peak is ~60d post-launch; covers +90d post-peak)

Each config is the corresponding no-wind model config with only the wind
Mdot and target_days modified.
"""
import yaml
from pathlib import Path

BASE_DIR = Path("/viper/ptmp1/sbrennan/SN2024pba")
GRID_DIR = BASE_DIR / "SN_B_v2_T6k_grid"

# (parent model name, mdot) pairs
PARENTS = ["T6k_M7.0_V10000", "T6k_M10.0_V10000"]
MDOTS = [1e-2, 3e-2]          # Msun/yr
TARGET_DAYS = 150.0           # peak ~60d + ~90d post-peak coverage


def _p(x, nd=1):
    """7.0 -> '7p0', 2.5 -> '2p5', 10.0 -> '10p0' (decimal point -> p)."""
    return f"{x:.{nd}f}".replace(".", "p")


def _e(x):
    """10000 -> '1e4', 0.03 -> '3e-2' (compact e-notation, no + or 0-pad)."""
    s = f"{x:.0e}"                      # '1e+04', '3e-02'
    s = s.replace("e+0", "e").replace("e+", "e").replace("e-0", "e-")
    return s


def grid_name(mass_msun, vel_kms, mdot_msun_yr=None):
    """X_M7p0_V1e4_W3e-2 : shell mass (Msun), velocity (km/s), wind Mdot.

    The wind suffix is omitted for the base 1e-3 Msun/yr wind.
    """
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
                f"T6k high-Mdot sample: H/He shell M={mass} Msun, "
                f"v={vel} km/s, solar comp; wind Mdot={mdot:.0e} Msun/yr, "
                f"v=100 km/s; 400d-old SN_A_g_T20k_Ni06 ejecta; "
                f"Ni56=1e-4 Msun. STELLA to {TARGET_DAYS:.0f} d. "
                f"Target: M_peak ~ -18.5 with decay ~2 mag/100d."
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
            cfg["x_model_from_sn"] = x

            st = dict(parent["stella"])
            st["target_days"] = TARGET_DAYS
            cfg["stella"] = st

            out = dict(parent["output"])
            out["model_dir"] = f"./SN_B_v2_T6k_grid/{model_name}"
            cfg["output"] = out

            header = (
                f"# SN_B_v2 T6k high-Mdot sample: {model_name}\n"
                f"# Higher wind Mdot to slow post-peak decay "
                f"(observed ~2 mag/100d)\n"
                f"# Shell: M={mass} Msun, v={vel} km/s, solar composition\n"
                f"# Wind: Mdot={mdot:.0e} Msun/yr, v=100 km/s; "
                f"Ni56=1e-4 Msun\n"
                f"# Source: SN_A_g_T20k_Ni06 (400-day-old SN_A ejecta)\n"
                f"# STELLA evolution: {TARGET_DAYS:.0f} days; STELLA only\n\n"
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
