#!/usr/bin/env python3
"""Generate the SW2 follow-up batch: Ni-free n_rho=8 shell + wind + ejecta.

The compact RSG-shell batch (n_rho=11, t_start~1.4 d) peaked at only
Mo ~ -13.4 -- ~5 mag below the -18.54 target, and ~2 mag below the old
expanded-shell SW twin (-15.28 -> -18.06 going n_rho 11 -> 8).  This
batch returns to the SW geometry: a stripped-star shell (n_rho=8,
Dessart 2010jl convention) launched early enough to expand across the
cavity before the STELLA start, a continuous rho ~ r^-2 wind filling
the gap, and the SN_A ejecta at its homologous 330 d radius.

Bookkeeping (same convention as the CM_/SW_ batches):
  t_start_days = 330 * v_in / V_shell  with v_in = 2070 km/s the SN_A
  ejecta inner-edge velocity, so the shell just reaches the ejecta at
  STELLA start and effective SN_A age = source_age + t_start = 330 d.

All models are interaction-only: cold shell (T_interface = T_max =
2000 K), Ni56 disabled, and every primary carries
paired_wind_control so the orchestrator auto-generates the no-ejecta
windctl sibling for the component-attribution test.

Axes: peak boost via shell KE (M, V) and tail sustainment via wind
mass-loading (mdot).  The known deficits of the -18.06 SW model are a
too-fast post-peak decline (~2.4 mag/20 d vs obs 2.3/100 d) and no
early -15.6 plateau; more wind ahead of the shock is the Ni-free knob
that extends dissipation past peak.
"""
import os
import sys

import yaml

BASE = "/viper/ptmp1/sbrennan/SN2024pba"
CONFIGS = os.path.join(BASE, "configs")
PARENT = os.path.join(CONFIGS, "RSG_M5_V07k_W3m4.yaml")

V_IN_KMS = 2070.0          # SN_A ejecta inner-edge velocity (~330 d)
EFFECTIVE_AGE = 330.0      # SN_A age at STELLA t=0, all models

# (name, M_shell, V_shell, mdot, rationale)
VARIANTS = [
    ("SW2_M5_V07k_W3m4", 5.0, 7000.0, 3e-4,
     "reproduction of SW_M5_V07k_W3m4 (-18.06): validates the pipeline "
     "after the stradsep5tt T-floor patch"),
    ("SW2_M5_V08k_W3m4", 5.0, 8000.0, 3e-4,
     "KE +31% at fixed mass: tests whether the -0.5 mag peak deficit "
     "closes with faster shock"),
    ("SW2_M7_V07k_W3m4", 7.0, 7000.0, 3e-4,
     "KE +40% via mass: later/broader dissipation, tests tail slowing"),
    ("SW2_M5_V07k_W5m4", 5.0, 7000.0, 5e-4,
     "denser wind: more mass ahead of the shock -> sustained post-peak "
     "dissipation (the tail is the residual mismatch)"),
    ("SW2_M5_V07k_W2m4", 5.0, 7000.0, 2e-4,
     "thinner-wind bracket: checks whether W3m4 is already at the "
     "mass-loading optimum for this gap"),
    ("SW2_M5_V10k_W3m4", 5.0, 10000.0, 3e-4,
     "upper-KE bracket (2x base): RSG V10k stalled at day ~2.8, but "
     "that was the compact-shell geometry; watch for the same stall"),
]


def t_start_for(v_kms):
    """Cavity-crossing time: shell reaches ejecta inner edge at t_start."""
    ts = EFFECTIVE_AGE * V_IN_KMS / v_kms
    return ts, EFFECTIVE_AGE - ts  # (t_start_days, source_age_days)


def main():
    with open(PARENT) as f:
        parent = yaml.safe_load(f)

    for name, m, v, mdot, why in VARIANTS:
        if len(name) > 28:
            sys.exit(f"name too long: {name} ({len(name)} chars)")
        ts, sa = t_start_for(v)
        cfg = dict(parent)
        cfg["model_name"] = name
        cfg["description"] = (
            f"SW2 Ni-free shell+wind+ejecta: n_rho=8 shell "
            f"(M={m} Msun, v={v:.0f} km/s, cold 2000 K) pre-expanded "
            f"for {ts:.1f} d, mdot={mdot} wind to SN_A ejecta at "
            f"330 d.  {why}."
        )
        x = dict(parent["x_model_from_sn"])
        x.update(
            {
                "source_age_days": sa,
                "shell_mass_msun": m,
                "shell_velocity_kms": v,
                "n_rho": 8.0,
                "T_interface_K": 2000.0,
                "T_max_K": 2000.0,
                "t_start_days": ts,
                "wind_T_K": 2000.0,
                "wind_mdot_msun_yr": mdot,
                # Narrow gap: the wind-ejecta contact discontinuity is a
                # physical 2.3-dex jump that cannot fit a <16 gradient ramp
                # inside the 0.05-dex gap.  The original SW batch ran stably
                # with grad~211 (verified on the archived CM twin .hyd).
                "allow_steep_interface": True,
            }
        )
        cfg["x_model_from_sn"] = x

        ni_cfg = dict(parent.get("stella", {}).get("nickel56", {}))
        ni_cfg.update({"enabled": False, "total_mass_msun": 0.0})
        ni_cfg.pop("distribution_mass_msun", None)
        ni_cfg.pop("distribution_type", None)
        ni_cfg.pop("inner_zone_fraction", None)
        st = dict(parent.get("stella", {}))
        st["nickel56"] = ni_cfg
        cfg["stella"] = st

        out = dict(parent["output"])
        out["model_dir"] = f"./SN_B_cocoon/{name}"
        cfg["output"] = out

        path = os.path.join(CONFIGS, f"{name}.yaml")
        with open(path, "w") as f:
            f.write("# SW2 batch: Ni-free n_rho=8 shell + wind + ejecta.\n")
            f.write(f"# {why}\n")
            f.write(f"# Parent config: {PARENT}\n")
            yaml.dump(cfg, f, default_flow_style=False, sort_keys=False)
        print("wrote", path)


if __name__ == "__main__":
    main()
