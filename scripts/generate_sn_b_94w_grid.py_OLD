#!/usr/bin/env python3
"""Generate the SN_B_v2_94W grid aimed at a fainter peak + smoother decay.

Round 1 showed the bright narrow spike (M_r ~ -19.65) is the shell
breakout flash: it is insensitive to wind Mdot (5e-2 -> 3e-3) and to
shell velocity (4000 -> 2500 km/s).  The broad interaction peak sits at
M_r ~ -18.6..-19.0, ~1 mag brighter than the observed -17.9.

Round 2 therefore attacks the breakout energetics and geometry:
  * shell mass 8 -> 4 Msun         (halve breakout/deposited energy)
  * shell velocity 4000 -> 3000    (with mass cut: ~4x less KE)
  * n_rho 8 -> 4                   (softer shell edge -> smeared breakout)
  * sw_transition_dex 0.5 -> 1.0   (broader shell->wind transition)
  * Mdot 0.05 -> 0.1               (more tau above breakout: smother test)
  * M4 + V3k + matched r^-1.5 CSM  (combined faint+smooth candidate)

Round 2 confirmed the spike is a startup transient: identical to <0.03
mag across every hydrodynamic variation.  Round 3 explores the shell
mass itself at fixed v=4000 km/s, Mdot=1e-2:
  * M = 2, 6 Msun            (with M4, M8 -> mass leverage 2-8)
  * M4, frac_v0 = 0.5        (more mass in the power-law envelope)
  * M4, t_start = 16 d       (same mass at 8x volume -> lower density)

Each config is the parent configs/SN_B_v2_94W.yaml with only the listed
keys changed.  All use the parallel (PARDISO) OpenStella build inherited
from the parent.  STELLA only, 150 d.
"""
import yaml
from pathlib import Path

BASE_DIR = Path("/viper/ptmp1/sbrennan/SN2024pba")
GRID_DIR = BASE_DIR / "SN_B_v2_94W_grid"
PARENT_CONFIG = BASE_DIR / "configs" / "SN_B_v2_94W.yaml"

# name_suffix -> dict of x_model_from_sn overrides (on top of the parent)
GRID = {
    # ---- round 1: wind Mdot x shell velocity + matched-CSM variants ----
    "V4k_W1e-2": {"shell_velocity_kms": 4000, "wind_mdot_msun_yr": 1e-2},
    "V4k_W3e-3": {"shell_velocity_kms": 4000, "wind_mdot_msun_yr": 3e-3},
    "V2p5k_W1e-2": {"shell_velocity_kms": 2500, "wind_mdot_msun_yr": 1e-2},
    "V2p5k_W3e-3": {"shell_velocity_kms": 2500, "wind_mdot_msun_yr": 3e-3},
    "V4k_Wmatch": {
        "shell_velocity_kms": 4000,
        "csm_density_match": True,
        "wind_mdot_msun_yr": 0.0,
    },
    "V4k_PLn1p5": {
        "shell_velocity_kms": 4000,
        "beta_wind": False,
        "powerlaw_csm": True,
        "csm_powerlaw_n": 1.5,
        "csm_density_match": True,
        "wind_mdot_msun_yr": 0.0,
    },
    # ---- round 2: breakout energetics / geometry ----
    "M4_V4k_W1e-2": {
        "shell_mass_msun": 4.0,
        "shell_velocity_kms": 4000,
        "wind_mdot_msun_yr": 1e-2,
    },
    "M4_V3k_W1e-2": {
        "shell_mass_msun": 4.0,
        "shell_velocity_kms": 3000,
        "wind_mdot_msun_yr": 1e-2,
    },
    "n4_V4k_W1e-2": {
        "n_rho": 4.0,
        "shell_velocity_kms": 4000,
        "wind_mdot_msun_yr": 1e-2,
    },
    "SW1p0_V4k_W1e-2": {
        "sw_transition_dex": 1.0,
        "shell_velocity_kms": 4000,
        "wind_mdot_msun_yr": 1e-2,
    },
    "V4k_W1e-1": {"shell_velocity_kms": 4000, "wind_mdot_msun_yr": 1e-1},
    "M4_V3k_PLn1p5": {
        "shell_mass_msun": 4.0,
        "shell_velocity_kms": 3000,
        "beta_wind": False,
        "powerlaw_csm": True,
        "csm_powerlaw_n": 1.5,
        "csm_density_match": True,
        "wind_mdot_msun_yr": 0.0,
    },
    # ---- round 3: shell mass leverage + mass distribution ----
    "M2_V4k_W1e-2": {
        "shell_mass_msun": 2.0,
        "shell_velocity_kms": 4000,
        "wind_mdot_msun_yr": 1e-2,
    },
    "M6_V4k_W1e-2": {
        "shell_mass_msun": 6.0,
        "shell_velocity_kms": 4000,
        "wind_mdot_msun_yr": 1e-2,
    },
    "M4_fV05_V4k_W1e-2": {
        "shell_mass_msun": 4.0,
        "shell_velocity_kms": 4000,
        "v0_kms": 2000.0,
        "wind_mdot_msun_yr": 1e-2,
    },
    "M4_ts16_V4k_W1e-2": {
        "shell_mass_msun": 4.0,
        "shell_velocity_kms": 4000,
        "t_start_days": 16.0,
        "wind_mdot_msun_yr": 1e-2,
    },
}


def main():
    with open(PARENT_CONFIG) as f:
        parent = yaml.safe_load(f)

    GRID_DIR.mkdir(parents=True, exist_ok=True)
    configs = []
    for suffix, over in GRID.items():
        model_name = f"94W_{suffix}"
        model_dir = GRID_DIR / model_name
        config_path = model_dir / "config.yaml"
        if config_path.exists():
            print(f"  Skipping existing: {config_path}")
            configs.append(str(config_path))
            continue
        model_dir.mkdir(parents=True, exist_ok=True)

        x = dict(parent["x_model_from_sn"])
        x.update(over)
        # v0_kms must stay below shell_velocity_kms; keep the parent ratio
        # unless the variant pins v0_kms itself.
        if "v0_kms" not in over:
            x["v0_kms"] = x["shell_velocity_kms"] * parent["x_model_from_sn"][
                "v0_kms"
            ] / parent["x_model_from_sn"]["shell_velocity_kms"]

        cfg = dict(parent)
        cfg["model_name"] = model_name
        cfg["x_model_from_sn"] = x
        desc = ", ".join(f"{k}={v}" for k, v in over.items())
        cfg["description"] = (
            f"94W grid {suffix}: {desc}; unsmoothed SN_A_g_T20k_Ni06 "
            "ejecta @1000d clipped at 10000 km/s. STELLA 150 d, STELLA only."
        )

        out = dict(parent["output"])
        out["model_dir"] = f"./SN_B_v2_94W_grid/{model_name}"
        cfg["output"] = out

        header = (
            f"# SN_B_v2_94W grid model: {model_name}\n"
            f"# overrides vs parent: {desc}\n"
            f"# Parallel STELLA: /ptmp/sbrennan/stella/OpenStella_parallel\n\n"
        )
        with open(config_path, "w") as f:
            f.write(header)
            yaml.dump(cfg, f, default_flow_style=False, sort_keys=False)
        configs.append(str(config_path))
        print(f"  Written: {config_path}")

    print(f"\nTotal configs in {GRID_DIR}/: {len(configs)}")
    print("\nSubmit with:")
    for c in configs:
        print(f"  python scripts/run_model_grid.py {c}")


if __name__ == "__main__":
    main()
