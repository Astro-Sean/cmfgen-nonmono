#!/usr/bin/env python3
"""Generate the SW5 batch: systematic decline-rate scan.

The completed grid now brackets the observed peak and rise:
  - peak amplitude is set by target density: SW3 A150-A200 reach
    -18.3 to -18.9 vs obs -18.2
  - early sustainment scales with shell KE/mass: SW2 M7 held ~1 mag
    brighter than M5 over the first 30 d post-peak
but every model still decays ~3-8x too fast (obs ~1.5 mag/100 d).

SW5 scans the two levers not yet systematically explored for the
post-peak decline, all with the SN_A ejecta retained:

  Mass axis (V=7000, dense target): M2-M10.  A heavier piston carries
    more momentum and sweeps more target mass; per Khatami & Kasen a
    larger shocked mass powers a stronger shock-cooling tail, and
    eta = M_ejecta/M_shell moves the CI regime.
  Velocity axis (M=5): V05k-V06k.  A slower shell stretches every
    timescale ~V^-1 while the dense target compensates the ~v^3
    drop in peak power -- the combination should slow the decline
    while keeping the peak near obs.
  Age fill: A175 interpolates between A150 and A200.

Bookkeeping (SW2/SW3 convention, per contact age T_c):
  t_start_days = T_c * v_in / V_shell   (v_in = 2070 km/s)
  source_age   = T_c - t_start          (SN_A age at SN_B launch)

Ni56=0 throughout.  paired_wind_control is OFF for this batch --
the SW2/SW3 windctls already provide the attribution baselines and
queue time is better spent on ejecta-bearing runs.
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
    ("SW5_M2_V07k_A200", 2.0, 7000.0, 200.0,
     "light piston on moderate-density target: eta~1 bracket, "
     "decelerates fastest - predicted to decline fastest"),
    ("SW5_M7_V07k_A200", 7.0, 7000.0, 200.0,
     "heavy piston + A200 target: extends the M7 sustainment seen at "
     "the 330 d baseline onto the right-peak density"),
    ("SW5_M10_V07k_A150", 10.0, 7000.0, 150.0,
     "heaviest piston + densest useful target: maximal swept-mass "
     "budget for the shock-cooling tail, may overshoot the peak"),
    ("SW5_M10_V07k_A200", 10.0, 7000.0, 200.0,
     "heavy piston + moderate target: balance point if A150 overshoots"),
    ("SW5_M10_V07k_A250", 10.0, 7000.0, 250.0,
     "heavy piston + thin target: if KE alone does the decline, this "
     "is the least-obscured version"),
    ("SW5_M5_V05k_A150", 5.0, 5000.0, 150.0,
     "slow shell + dense target: stretches all timescales x1.4 while "
     "the dense target keeps the peak near obs"),
    ("SW5_M5_V05k_A200", 5.0, 5000.0, 200.0,
     "slow shell + moderate target: peak ~-17.5 expected, slowest "
     "crossing time of the grid"),
    ("SW5_M5_V06k_A150", 5.0, 6000.0, 150.0,
     "intermediate velocity + dense target: fills the V-axis between "
     "V05k and V07k"),
    ("SW5_M5_V06k_A200", 5.0, 6000.0, 200.0,
     "intermediate velocity + moderate target"),
    ("SW5_M5_V07k_A175", 5.0, 7000.0, 175.0,
     "age interpolation between A150 and A200 where the peak matched"),
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
            f"SW5 decline-scan: n_rho=8 shell (M={m} Msun, v={v:.0f} "
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
                "paired_wind_control": False,
            }
        )
        cfg["x_model_from_sn"] = x

        out = dict(parent["output"])
        out["model_dir"] = f"./SN_B_cocoon/{name}"
        cfg["output"] = out

        path = os.path.join(CONFIGS, f"{name}.yaml")
        with open(path, "w") as f:
            f.write("# SW5 batch: systematic decline-rate scan "
                    "(mass x velocity x contact age).\n")
            f.write(f"# {why}\n")
            f.write(f"# Parent config: {PARENT}\n")
            yaml.dump(cfg, f, default_flow_style=False, sort_keys=False)
        print("wrote", path)


if __name__ == "__main__":
    main()
