#!/usr/bin/env python3
"""Generate batch 6 of the SN_B cocoon/internal-shock tests.

Batch 4/5 (cold shell) established that pure-interaction lightcurves
saturate at Mo ~ -14.5 to -15.1 with ~1% radiative coupling: the
shock energy goes into bulk ejecta motion and the photosphere stays
cold (~1500 K) and deep (~5e14 cm).  Even KE=1e51 (5 Msun @ 10k) does
not break the floor.

The observed SN_B needs ~3e49-5e49 erg radiated over ~50 d around a
main peak at ~+40 d after the first faint bump.  Two honest ways to
get there inside the cocoon scenario:

(a) SN_B is a genuine second explosion: the "shell" carries Ni56 in
    its inner zones, so the early plateau is the interaction flash and
    the main -18.5 peak is radioactive.  Ni ~0.03-0.05 Msun peaks at
    Mbol ~ -18.5 to -19 at ~40-60 d -- exactly the observed delay.

(b) Pure interaction with a maximal shell (10 Msun @ 10k = 1e52 erg):
    tests whether the cocoon can reach -18.5 without radioactivity.

Detached thin CSM (2 Msun at 850 km/s) is retained in all models for
the narrow P-Cygni, but it is tau~0.01 (transparent): the narrow line
is kinematic, not powered by the CSM blanket.
"""
import os
import sys

import yaml

BASE = "/viper/ptmp1/sbrennan/SN2024pba"
GRID = os.path.join(BASE, "SN_B_cocoon")
PARENT = os.path.join(BASE, "SN_B_cocoon/CD_M1_V10k_CM2/config.yaml")

# (name, M_shell, V_shell, M_CSM, Ni56, rationale)
VARIANTS = [
    ("CF_M5_V10k_CM2_Ni05", 5.0, 10000.0, 2.0, 0.05,
     "SN-scale shell + Ni56: interaction bump + radioactive main peak"),
    ("CF_M2_V10k_CM2_Ni05", 2.0, 10000.0, 2.0, 0.05,
     "lighter shell, same Ni56: isolates interaction vs Ni share"),
    ("CF_M2_V10k_CM2_Ni03", 2.0, 10000.0, 2.0, 0.03,
     "Ni low end"),
    ("CF_M5_V10k_CM1_Ni05", 5.0, 10000.0, 1.0, 0.05,
     "thinner CSM blanket under same Ni56"),
    ("CF_M10_V10k_CM2", 10.0, 10000.0, 2.0, 0.0,
     "pure-interaction ceiling: 1e52 erg, no Ni56"),
]


def main():
    with open(PARENT) as f:
        parent = yaml.safe_load(f)

    for name, m, v, cm, ni, why in VARIANTS:
        if len(name) > 28:
            sys.exit(f"name too long: {name} ({len(name)} chars)")
        cfg = dict(parent)
        cfg["model_name"] = name
        cfg["description"] = (
            f"Cocoon batch-6 (cold shell, Ni56={ni}): H/He shell "
            f"(M={m} Msun, v={v:.0f} km/s, T=2000 K) in immediate "
            f"contact with SN_A ejecta + detached outer CSM "
            f"(M={cm:.1f} Msun, 850 km/s, r=6e16-1.2e17 cm).  {why}."
        )
        x = dict(parent["x_model_from_sn"])
        x.update(
            {
                "shell_mass_msun": m,
                "shell_velocity_kms": v,
                "v0_kms": 0.3 * v,
                "outer_csm_mass_msun": cm,
                "outer_csm_r_inner_cm": 6.0e16,
                "outer_csm_r_outer_cm": 1.2e17,
            }
        )
        cfg["x_model_from_sn"] = x

        ni_cfg = dict(parent.get("stella", {}).get("nickel56", {}))
        ni_cfg.update(
            {
                "enabled": ni > 0,
                "total_mass_msun": ni,
                "distribution_mass_msun": min(0.5, m),
                "distribution_type": "uniform",
                "inner_zone_fraction": 0.15,
            }
        )
        st = dict(parent.get("stella", {}))
        st["nickel56"] = ni_cfg
        cfg["stella"] = st

        out = dict(parent["output"])
        out["model_dir"] = f"./SN_B_cocoon/{name}"
        cfg["output"] = out

        mdir = os.path.join(GRID, name)
        os.makedirs(mdir, exist_ok=True)
        path = os.path.join(mdir, "config.yaml")
        with open(path, "w") as f:
            f.write("# Cocoon batch-6: cold shell + Ni56 / max-KE tests.\n")
            f.write(f"# {why}\n")
            f.write(f"# Parent config: {PARENT}\n")
            yaml.dump(cfg, f, default_flow_style=False, sort_keys=False)
        print("wrote", path)


if __name__ == "__main__":
    main()
