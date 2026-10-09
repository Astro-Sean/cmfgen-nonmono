#!/usr/bin/env python3
"""Generate the optimized SN_B variant batch (180 d STELLA runs).

Parent: SN_B_double_peak_grid/SB_M10_ts8_W3e-4 (best double-peak fit).

Motivation from the four-criteria trend review (TRENDS.md):

1. Every double-peak model peaks at t ~ 57-65 d, so a 100 d run leaves
   <26 d of post-peak coverage and the decay criterion is untestable.
   All variants here use stella.target_days=180 (nstmax raised to 4e6).
2. The completed W1e-4 / vw300 runs (100 d) went single-peaked and dim
   (-16.4): lowering wind density 3x suppressed the second peak within
   83 d.  At 180 d we can see whether the peak is delayed rather than
   absent, and measure the plateau level properly.
3. Shell KE (v_shell) dims the main peak; wind density dims the plateau.
   Combined variants test both levers at once.

SN_A ejecta source (SN_A_g_T20k_Ni06 @ 400 d) is fixed in all models.
"""
import os
import sys

import yaml

BASE = "/viper/ptmp1/sbrennan/SN2024pba"
GRID = os.path.join(BASE, "SN_B_double_peak_grid")
PARENT = os.path.join(GRID, "SB_M10_ts8_W3e-4", "config.yaml")
PARENT_M8 = os.path.join(GRID, "SB_M8_ts8_W3e-4", "config.yaml")

TARGET_DAYS = 180.0

# (name, parent, {x_model_from_sn patch}, rationale)
VARIANTS = [
    ("SB_M10_ts8_W3e-4_t180", PARENT, {},
     "incumbent winner at 180 d: measures the post-peak decay "
     "(unmeasurable in the 100 d run; peak at t=57 d)"),
    ("SB_M8_ts8_W3e-4_t180", PARENT_M8, {},
     "M8 runner-up at 180 d: best rise time (21 d vs obs 19) in the "
     "100 d run; measures decay for the lighter shell"),
    ("SB_M10_ts8_W1e-4_t180", PARENT, {"wind_mdot_msun_yr": 1e-4},
     "mdot/3 at 180 d: 100 d run was single-peaked -16.4; check if the "
     "second peak is delayed rather than absent, and the decay"),
    ("SB_M10_ts8_vw300_t180", PARENT, {"wind_velocity_kms": 300.0},
     "v_wind x3 at 180 d: same D=Mdot/v_w /3 as W1e-4; also single-peaked "
     "at 100 d; toward observed ~850 km/s CSM velocity"),
    ("SB_M10_ts8_vw850_t180", PARENT, {"wind_velocity_kms": 850.0},
     "v_wind = 850 km/s = observed narrow H-alpha P-Cyg velocity; "
     "D /8.5, thinnest wind in the batch"),
    ("SB_M10_ts8_V9k_t180", PARENT, {"shell_velocity_kms": 9000.0},
     "v_shell 10000 -> 9000: dims main peak via lower shell KE "
     "(plateau ~ -0.34 mag)"),
    ("SB_M10_ts8_V8k_t180", PARENT, {"shell_velocity_kms": 8000.0},
     "v_shell 10000 -> 8000: stronger KE dim (plateau ~ -0.7 mag)"),
    ("SB_M10_ts8_W1e-4vw300_t180", PARENT,
     {"wind_mdot_msun_yr": 1e-4, "wind_velocity_kms": 300.0},
     "thinnest wind (D /9): plateau bracket lower edge"),
    ("SB_M10_ts8_vw300V9k_t180", PARENT,
     {"shell_velocity_kms": 9000.0, "wind_velocity_kms": 300.0},
     "optimized combo: vw300 dims plateau ~ -1.2 mag, V9k dims main "
     "peak via lower KE; targets plateau -15.6 and peak -18.5"),
    # --- single-parameter interpolation of the double/single transition.
    # Wind density D = mdot/v_w is the control parameter: W1e-4 and vw300
    # (both D=1e-6) gave identical lightcurves (-16.44 vs -16.43).
    # D axis: 3e-4/100=3e-6 (double) -> these points -> 1e-6 (single).
    ("SB_M10_ts8_W25e-5_t180", PARENT, {"wind_mdot_msun_yr": 2.5e-4},
     "D=2.5e-6: thin-side transition point between double (3e-6) and "
     "single (1e-6) morphology"),
    ("SB_M10_ts8_W2e-4_t180", PARENT, {"wind_mdot_msun_yr": 2e-4},
     "D=2e-6: thin-side transition point"),
    ("SB_M10_ts8_W15e-5_t180", PARENT, {"wind_mdot_msun_yr": 1.5e-4},
     "D=1.5e-6: thin-side transition point just above the 1e-6 single"),
    # --- dense side of the transition: double (3e-6) -> bright single
    # (1e-5 at mdot=1e-3).
    ("SB_M10_ts8_W5e-4_t180", PARENT, {"wind_mdot_msun_yr": 5e-4},
     "D=5e-6: dense-side transition point"),
    ("SB_M10_ts8_W75e-5_t180", PARENT, {"wind_mdot_msun_yr": 7.5e-4},
     "D=7.5e-6: dense-side transition point"),
]


def main():
    for name, parent_path, patch, why in VARIANTS:
        if len(name) > 28:
            sys.exit(f"name too long: {name} ({len(name)} chars)")
        with open(parent_path) as f:
            parent = yaml.safe_load(f)
        cfg = dict(parent)
        cfg["model_name"] = name
        cfg["description"] = (
            f"Optimized 180 d variant of {os.path.basename(os.path.dirname(parent_path))} "
            f"(best double-peak family). {why}. SN_A source fixed "
            f"(SN_A_g_T20k_Ni06 @ 400 d).")
        x = dict(parent["x_model_from_sn"])
        x.update(patch)
        cfg["x_model_from_sn"] = x
        st = dict(parent["stella"])
        st["target_days"] = TARGET_DAYS
        st["nstmax"] = 4000000
        cfg["stella"] = st
        out = dict(parent["output"])
        out["model_dir"] = f"./SN_B_double_peak_grid/{name}"
        cfg["output"] = out

        mdir = os.path.join(GRID, name)
        os.makedirs(mdir, exist_ok=True)
        path = os.path.join(mdir, "config.yaml")
        with open(path, "w") as f:
            f.write(f"# Optimized 180 d variant: {why}\n")
            f.write(f"# Parent config: {parent_path}\n")
            yaml.dump(cfg, f, default_flow_style=False, sort_keys=False)
        print("wrote", path)


if __name__ == "__main__":
    main()
