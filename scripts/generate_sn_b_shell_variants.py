#!/usr/bin/env python3
"""Generate shell-configuration variants of the best-fit SN_B model.

Parent: SN_B_double_peak_grid/SB_M10_ts8_W3e-4 (best double-peak fit,
2026-09-21 scoring).  Goal: dim the first plateau (~-16.5, target -15.5)
by ~1 mag while keeping the second peak (~-19.4, target -18.54) and the
dip morphology.

The first plateau is powered by shell->wind forward-shock dissipation:
photosphere sits at r ~ 2-4e15 cm, v ~ 8000 km/s, rho ~ 1e-15 (wind
densities).  L_plateau ~ 0.5 * (Mdot/v_wind) * v_s**3, so the levers are
wind density (Mdot, v_wind) and shell velocity; n_rho controls how much
mass sits in the fast leading envelope.
"""
import os
import sys

import yaml

BASE = "/viper/ptmp1/sbrennan/SN2024pba"
GRID = os.path.join(BASE, "SN_B_double_peak_grid")
PARENT = os.path.join(GRID, "SB_M10_ts8_W3e-4", "config.yaml")

# (name, {x_model_from_sn key: value}, rationale)
VARIANTS = [
    ("SB_M10_ts8_W1e-4",
     {"wind_mdot_msun_yr": 1e-4},
     "mdot 3e-4 -> 1e-4: plateau ~ -1.2 mag (L ~ Mdot)"),
    ("SB_M10_ts8_vw300",
     {"wind_velocity_kms": 300.0},
     "v_wind 100 -> 300 km/s: D=Mdot/v_w /3 -> plateau ~ -1.2 mag; "
     "toward observed ~850 km/s narrow-P-Cyg CSM velocity"),
    ("SB_M10_ts8_vw850",
     {"wind_velocity_kms": 850.0},
     "v_wind -> 850 km/s = TNS narrow H-alpha P-Cyg absorption velocity; "
     "D /8.5 -> plateau ~ -2.3 mag (dim bracket)"),
    ("SB_M10_ts8_V9k",
     {"shell_velocity_kms": 9000.0},
     "v_shell 10000 -> 9000: plateau x0.73 (~ -0.34 mag), peak 2 dims too"),
    ("SB_M10_ts8_V8k",
     {"shell_velocity_kms": 8000.0},
     "v_shell 10000 -> 8000: plateau x0.51 (~ -0.7 mag), peak 2 dims too"),
    ("SB_M10_ts8_nr12",
     {"n_rho": 12.0},
     "n_rho 8 -> 12: steeper envelope, less fast leading-edge mass, "
     "weaker early shock power"),
    ("SB_M10_ts8_W1e-4vw300",
     {"wind_mdot_msun_yr": 1e-4, "wind_velocity_kms": 300.0},
     "mdot=1e-4 + v_wind=300: D /9 -> plateau ~ -2.4 mag (extreme bracket)"),
]


def main():
    with open(PARENT) as f:
        parent = yaml.safe_load(f)

    for name, patch, why in VARIANTS:
        if len(name) > 28:
            sys.exit(f"name too long: {name} ({len(name)} chars)")
        cfg = dict(parent)
        cfg["model_name"] = name
        cfg["description"] = (
            f"Shell-config variant of SB_M10_ts8_W3e-4 (best double-peak "
            f"fit). {why}. All other physics identical to parent.")
        x = dict(parent["x_model_from_sn"])
        x.update(patch)
        cfg["x_model_from_sn"] = x
        out = dict(parent["output"])
        out["model_dir"] = f"./SN_B_double_peak_grid/{name}"
        cfg["output"] = out

        mdir = os.path.join(GRID, name)
        os.makedirs(mdir, exist_ok=True)
        path = os.path.join(mdir, "config.yaml")
        with open(path, "w") as f:
            f.write(f"# Variant of SB_M10_ts8_W3e-4: {why}\n")
            f.write(f"# Parent config: {PARENT}\n")
            yaml.dump(cfg, f, default_flow_style=False, sort_keys=False)
        print("wrote", path)


if __name__ == "__main__":
    main()
