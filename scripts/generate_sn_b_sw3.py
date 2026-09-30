#!/usr/bin/env python3
"""Generate the SW3 batch: dense-target (young contact-age) variants.

The SW2 decline (~2.4 mag/20 d post-peak vs obs ~2.3 mag/100 d) fails
because the shock runs through an optically thin target: at the 330 d
contact age the SN_A ejecta column is too small to trap radiation or
decelerate the shock into the self-similar s=2 regime.  A steady wind
cannot supply the missing mass (swept wind mass in 150 d stays <1% of
the shell mass even at mdot=1e-1).

The density lever that remains is the contact age T_c: the SN_A ejecta
is homologous, rho ~ T_c^-3, so a younger target at contact is both
denser (more swept mass -> deceleration -> slow decline) and more
optically thick (diffusion tail).  This batch scans T_c while holding
the SW2 geometry otherwise fixed.

Bookkeeping (same convention as SW2, per contact age T_c):
  t_start_days = T_c * v_in / V_shell   (v_in = 2070 km/s)
  source_age   = T_c - t_start          (SN_A age at SN_B launch)

Members:
  A250/A200/A150/A120: target density x2.4/x4.5/x10.6/x20.8 vs 330 d,
    bracketing the optical depth where the diffusion-limited rise
    matches the observed ~40-50 d (A120 expected to overshoot).
  M2 at A150: light-piston falsification of the "low ejecta mass"
    hypothesis -- a 2 Msun piston decelerates faster, predicted to
    decline FASTER, not slower.
  M7 at A150: heavy piston, mass-matched to the swept target.
  V10k at A150: +KE for peak power at fixed density.

All Ni56=0, all paired_wind_control for the ejecta-attribution test.
"""
import os
import sys

import yaml

BASE = "/viper/ptmp1/sbrennan/SN2024pba"
CONFIGS = os.path.join(BASE, "configs")
PARENT = os.path.join(CONFIGS, "SW2_M5_V07k_W3m4.yaml")

V_IN_KMS = 2070.0  # SN_A ejecta inner-edge velocity (homologous)

# (name, M_shell, V_shell, contact_age_d, rationale)
VARIANTS = [
    ("SW3_M5_V07k_A250", 5.0, 7000.0, 250.0,
     "target density x2.4 vs 330 d: first step off the SW2 baseline"),
    ("SW3_M5_V07k_A200", 5.0, 7000.0, 200.0,
     "target density x4.5: more swept mass, deeper photosphere"),
    ("SW3_M5_V07k_A150", 5.0, 7000.0, 150.0,
     "target density x10.6: expect tau~1-10 through the metal ejecta, "
     "diffusion rise approaching obs ~40-50 d"),
    ("SW3_M5_V07k_A120", 5.0, 7000.0, 120.0,
     "target density x20.8: overshoot bracket, likely too thick/slow "
     "a rise -- bounds the density axis"),
    ("SW3_M2_V07k_A150", 2.0, 7000.0, 150.0,
     "light piston on dense target: tests the low-ejecta-mass "
     "hypothesis; decelerates sooner, predicted to decline FASTER"),
    ("SW3_M7_V07k_A150", 7.0, 7000.0, 150.0,
     "heavy piston, mass-matched to swept target: predicted slowest "
     "decline of the A150 set"),
    ("SW3_M5_V10k_A150", 5.0, 10000.0, 150.0,
     "+KE bracket on the dense target: peak power ~ v^3, watch for "
     "the wind-interface stall seen on RSG V10k"),
]


def timing_for(v_kms, contact_age):
    """Shell reaches ejecta inner edge at t_start = T_c * v_in/v_shell."""
    ts = contact_age * V_IN_KMS / v_kms
    return ts, contact_age - ts  # (t_start_days, source_age_days)


def main():
    with open(PARENT) as f:
        parent = yaml.safe_load(f)

    for name, m, v, tc, why in VARIANTS:
        if len(name) > 28:
            sys.exit(f"name too long: {name} ({len(name)} chars)")
        ts, sa = timing_for(v, tc)
        cfg = dict(parent)
        cfg["model_name"] = name
        cfg["description"] = (
            f"SW3 dense-target: n_rho=8 shell (M={m} Msun, v={v:.0f} "
            f"km/s, cold 2000 K) pre-expanded {ts:.1f} d, mdot=3e-4 "
            f"wind, SN_A ejecta at contact age {tc:.0f} d "
            f"(rho x{(330.0 / tc) ** 3:.1f} vs 330 d).  {why}."
        )
        x = dict(parent["x_model_from_sn"])
        x.update(
            {
                "source_age_days": sa,
                "shell_mass_msun": m,
                "shell_velocity_kms": v,
                "t_start_days": ts,
            }
        )
        cfg["x_model_from_sn"] = x

        out = dict(parent["output"])
        out["model_dir"] = f"./SN_B_cocoon/{name}"
        cfg["output"] = out

        path = os.path.join(CONFIGS, f"{name}.yaml")
        with open(path, "w") as f:
            f.write("# SW3 batch: dense-target (young contact-age) Ni-free variants.\n")
            f.write(f"# {why}\n")
            f.write(f"# Parent config: {PARENT}\n")
            yaml.dump(cfg, f, default_flow_style=False, sort_keys=False)
        print("wrote", path)


if __name__ == "__main__":
    main()
