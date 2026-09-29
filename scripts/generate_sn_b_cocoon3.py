#!/usr/bin/env python3
"""Generate batch 3 of the SN_B cocoon/internal-shock tests.

Batch 2 (CB_*) showed the detached outer CSM smears the prompt flash
into a broad ~40-60 d peak with the observed delayed-rise morphology,
but the lightcurve sits at Mbol ~ -25 -- ~40x too luminous.  Energy
conservation is the binding constraint: the observed SN_B event
radiates ~3e49 erg (L~6e42 x ~50 d), while a 5 Msun shell at 10k km/s
carries 1e51 erg.  The CSM cannot hide energy, only delay it.

Batch 3 therefore drops the shell kinetic energy by ~10-20x while
keeping the batch-2 geometry (immediate contact with the SN_A cocoon +
2 Msun detached CSM shell at 850 km/s, tau~30).  KE targets
~3e49 - 2e50 erg: the observed peak should then emerge once the
trapped/reprocessed energy diffuses out on the CSM timescale.
"""
import os
import sys

import yaml

BASE = "/viper/ptmp1/sbrennan/SN2024pba"
GRID = os.path.join(BASE, "SN_B_cocoon")
PARENT = os.path.join(BASE, "SN_B_cocoon/CB_M5_V10k_CM2/config.yaml")

# (name, M_shell, V_shell, csm_mass_msun, KE_1e50, rationale)
VARIANTS = [
    ("CC_M05_V10k_CM2", 0.5, 10000.0, 2.0, 1.0,
     "lowest KE consistent with obs E_rad"),
    ("CC_M1_V10k_CM2", 1.0, 10000.0, 2.0, 2.0,
     "KE bracket high end"),
    ("CC_M1_V07k_CM2", 1.0, 7000.0, 2.0, 1.0,
     "same KE, slower piston: weaker flash + later arrival"),
    ("CC_M2_V07k_CM2", 2.0, 7000.0, 2.0, 2.0,
     "heavier slow piston: more momentum into cocoon"),
    ("CC_M1_V10k_CM4", 1.0, 10000.0, 4.0, 2.0,
     "low KE + denser blanket: maximum trapping"),
]


def main():
    with open(PARENT) as f:
        parent = yaml.safe_load(f)

    for name, m, v, cm, ke, why in VARIANTS:
        if len(name) > 28:
            sys.exit(f"name too long: {name} ({len(name)} chars)")
        cfg = dict(parent)
        cfg["model_name"] = name
        cfg["description"] = (
            f"Cocoon batch-3 (low-KE): H/He shell (M={m} Msun, "
            f"v={v:.0f} km/s, KE={ke:.1f}e50 erg) in immediate contact "
            f"with SN_A ejecta + detached outer CSM (M={cm:.1f} Msun, "
            f"850 km/s, r=6e16-1.2e17 cm).  {why}.  Ni56 disabled."
        )
        x = dict(parent["x_model_from_sn"])
        x.update(
            {
                "shell_mass_msun": m,
                "shell_velocity_kms": v,
                "v0_kms": 0.3 * v,
                "outer_csm_mass_msun": cm,
            }
        )
        cfg["x_model_from_sn"] = x

        out = dict(parent["output"])
        out["model_dir"] = f"./SN_B_cocoon/{name}"
        cfg["output"] = out

        mdir = os.path.join(GRID, name)
        os.makedirs(mdir, exist_ok=True)
        path = os.path.join(mdir, "config.yaml")
        with open(path, "w") as f:
            f.write("# Cocoon batch-3: internal shock + outer CSM, low KE.\n")
            f.write(f"# {why}\n")
            f.write(f"# Parent config: {PARENT}\n")
            yaml.dump(cfg, f, default_flow_style=False, sort_keys=False)
        print("wrote", path)


if __name__ == "__main__":
    main()
