#!/usr/bin/env python3
"""Generate the SW4 batch: compact-shell wind-CI morphology test.

Motivation (Khatami & Kasen 2024, ApJ 972, 140): the observed SN_B
morphology -- ~20 d plateau, fast ~20 d rise to peak, slow post-peak
decline -- matches the "continued interaction -> emergence/cooling"
sequence of an interior-breakout event.  In the thin (xi << 1) regime
the lightcurve tracks the instantaneous shock power L ~ rho v^3, so
the morphology is set by the density profile the shock traverses:

  * plateau:  shock crossing a wind cavity (rho ~ r^-2, s=2 gives
              L ~ t^-0.4, a slowly declining plateau)
  * fast rise: the density jump at the wind->ejecta edge (x~30 in rho
              = ~3.7 mag in L)
  * decline:  shock in the steep homologous ejecta (the known weak
              point; density ladder covered by SW3)

The SW2/RSG geometries fail the plateau for opposite reasons: SW2
pre-expands the shell so the wind cavity is only ~0.05 dex (no CI
phase at all), while the compact RSG n_rho=11 shell loses ~2.8 mag of
edge KE (too faint).  SW4 keeps the stripped n_rho=8 shell of SW2 but
at the COMPACT t_start of the RSG batch (~1.4 d), so a real ~2e15 cm
wind cavity exists for the shock to dissipate into, and raises mdot to
1e-3 so the wind phase radiates near the observed plateau level
(L ~ 2 pi r^2 rho v^3 ~ 1e42 erg/s at r~7e14 -> M ~ -16).

All Ni56=0, paired_wind_control for the ejecta-attribution test.
Timing convention follows the RSG batch: effective SN_A age 330 d at
STELLA t=0 (source_age + t_start).
"""
import os
import sys

import yaml

BASE = "/viper/ptmp1/sbrennan/SN2024pba"
CONFIGS = os.path.join(BASE, "configs")
PARENT = os.path.join(CONFIGS, "SW2_M5_V07k_W3m4.yaml")

# (name, M_shell, V_shell, mdot, t_start, source_age, rationale)
VARIANTS = [
    ("SW4_M5_V07k_W1m3c", 5.0, 7000.0, 1e-3, 1.4, 328.6,
     "compact n_rho=8 shell (t_start=1.4d) + mdot=1e-3 wind cavity: "
     "wind-CI plateau then ejecta-edge density-jump rise"),
    ("SW4_M5_V07k_W3m3c", 5.0, 7000.0, 3e-3, 1.4, 328.6,
     "brighter/thicker wind bracket: tests whether a denser cavity "
     "lifts the plateau to the obs -15.6 and slows the decline"),
]


def main():
    with open(PARENT) as f:
        parent = yaml.safe_load(f)

    for name, m, v, mdot, ts, sa, why in VARIANTS:
        if len(name) > 28:
            sys.exit(f"name too long: {name} ({len(name)} chars)")
        cfg = dict(parent)
        cfg["model_name"] = name
        cfg["description"] = (
            f"SW4 wind-CI morphology test: compact n_rho=8 shell "
            f"(M={m} Msun, v={v:.0f} km/s, cold 2000 K, t_start={ts}d "
            f"so shell edge ~{v * 1e5 * ts * 86400.0:.1e} cm), "
            f"mdot={mdot:g} wind spanning the cavity to SN_A ejecta "
            f"at ~330 d (source_age={sa} + t_start).  {why}."
        )
        x = dict(parent["x_model_from_sn"])
        x.update(
            {
                "source_age_days": sa,
                "shell_mass_msun": m,
                "shell_velocity_kms": v,
                "t_start_days": ts,
                "wind_mdot_msun_yr": mdot,
            }
        )
        cfg["x_model_from_sn"] = x

        out = dict(parent["output"])
        out["model_dir"] = f"./SN_B_cocoon/{name}"
        cfg["output"] = out

        path = os.path.join(CONFIGS, f"{name}.yaml")
        with open(path, "w") as f:
            f.write("# SW4 batch: compact shell + wind cavity -> K&K "
                    "wind-CI morphology test.\n")
            f.write(f"# {why}\n")
            f.write(f"# Parent config: {PARENT}\n")
            yaml.dump(cfg, f, default_flow_style=False, sort_keys=False)
        print("wrote", path)


if __name__ == "__main__":
    main()
