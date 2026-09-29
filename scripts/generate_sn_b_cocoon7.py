#!/usr/bin/env python3
"""Generate batch 7 of the SN_B cocoon/internal-shock tests.

Batch 6 established:
- Pure cold interaction saturates at Mo ~ -15.4 even at KE=1e52 erg
  (CF_M10): the cocoon geometry cannot reach -18.5 on kinetic power
  alone.  Radiative coupling is ~1%; extra KE becomes bulk motion.
- Ni56 in the shell gives the observed morphology: faint interaction
  floor (~-15) -> delayed radioactive second peak.  Ni05 -> -17.2,
  Ni03 -> -16.9.  Observed -18.5 needs Ni ~0.1-0.15 Msun.
- Peak delay scales with shell mass (25 d at 5 Msun vs 13 d at 2 Msun);
  the observed main peak is ~40 d after the faint bump.

Batch 7 scales Ni56 up and tests a denser cocoon (younger SN_A age,
source_age_days=200 -> ejecta ~(315/200)^3 ~ 3.9x denser -> hotter
photosphere and longer Ni diffusion time).
"""
import os
import sys

import yaml

BASE = "/viper/ptmp1/sbrennan/SN2024pba"
GRID = os.path.join(BASE, "SN_B_cocoon")
PARENT = os.path.join(BASE, "SN_B_cocoon/CF_M5_V10k_CM2_Ni05/config.yaml")

# (name, M_shell, V_shell, Ni56, sn_a_age_days, t_start_days, rationale)
VARIANTS = [
    ("CG_M5_V10k_CM2_Ni15", 5.0, 10000.0, 0.15, 315.0, 69.0,
     "Ni upper bound: expect peak ~-18.2..-18.5, delay ~30 d"),
    ("CG_M5_V10k_CM2_Ni10", 5.0, 10000.0, 0.10, 315.0, 69.0,
     "Ni mid: bracket the -18.5 target"),
    ("CG_M2_V10k_CM2_Ni10", 2.0, 10000.0, 0.10, 315.0, 69.0,
     "same Ni in lighter shell: earlier peak, hotter per-zone heat"),
    ("CG_M10_V10k_CM2_Ni05", 10.0, 10000.0, 0.05, 315.0, 69.0,
     "1e52 erg interaction floor plus Ni05 peak"),
    ("CG_M5_V10k_CM2_Ni10a200", 5.0, 10000.0, 0.10, 200.0, 44.0,
     "younger cocoon: ejecta ~4x denser -> hotter photosphere, "
     "slower Ni diffusion"),
]


def main():
    with open(PARENT) as f:
        parent = yaml.safe_load(f)

    for name, m, v, ni, age, tstart, why in VARIANTS:
        if len(name) > 28:
            sys.exit(f"name too long: {name} ({len(name)} chars)")
        cfg = dict(parent)
        cfg["model_name"] = name
        cfg["description"] = (
            f"Cocoon batch-7 (cold shell + Ni56={ni}): H/He shell "
            f"(M={m} Msun, v={v:.0f} km/s, T=2000 K) in immediate "
            f"contact with SN_A ejecta (source age {age:.0f} d) + "
            f"detached outer CSM (2 Msun, 850 km/s, 6e16-1.2e17 cm).  "
            f"{why}."
        )
        x = dict(parent["x_model_from_sn"])
        x.update(
            {
                "shell_mass_msun": m,
                "shell_velocity_kms": v,
                "v0_kms": 0.3 * v,
                "source_age_days": age,
                "t_start_days": tstart,
                "outer_csm_mass_msun": 2.0,
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
            f.write("# Cocoon batch-7: Ni scaling + denser cocoon.\n")
            f.write(f"# {why}\n")
            f.write(f"# Parent config: {PARENT}\n")
            yaml.dump(cfg, f, default_flow_style=False, sort_keys=False)
        print("wrote", path)


if __name__ == "__main__":
    main()
