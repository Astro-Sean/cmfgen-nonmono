#!/usr/bin/env python3
"""Generate batch 4 of the SN_B cocoon/internal-shock tests.

Post-mortem on batches 2-3: the .tt bolometric lightcurves radiate
~1.8e52-5.7e52 erg over ~150 d, while the available kinetic energy is
only ~2e50-1e51 erg.  The excess comes from the Dessart+2015 shell
initialization: T_max_K=65000 makes the low-density shell radiation
dominated (E_rad ~ aT^4 V ~ 3e52 erg at start), and eve2 relaxes the
radiation field in the ejecta/gap/CSM to Tr ~ 4-20 kK as well
(~1e52 erg there).  In the evacuated-cavity geometry that reservoir
leaks out as a fake interaction lightcurve; in the dense-wind grids it
stayed trapped and went into pdV work, which is why the failure was not
visible before.

Batch 4 keeps the batch-3 geometry (immediate contact + detached CSM at
850 km/s) but starts the shell cold: T_interface_K=T_max_K=2000 K.
All heating then comes from the shocks, so the lightcurve amplitude is
set by the shell/ejecta kinetic energy as intended.
"""
import os
import sys

import yaml

BASE = "/viper/ptmp1/sbrennan/SN2024pba"
GRID = os.path.join(BASE, "SN_B_cocoon")
PARENT = os.path.join(BASE, "SN_B_cocoon/CB_M5_V10k_CM2/config.yaml")

# (name, M_shell, V_shell, csm_mass_msun, KE_1e50, rationale)
VARIANTS = [
    ("CD_M05_V10k_CM2", 0.5, 10000.0, 2.0, 1.0,
     "cold-shell anchor at lowest KE"),
    ("CD_M1_V10k_CM2", 1.0, 10000.0, 2.0, 2.0,
     "cold-shell counterpart of CC_M1_V10k_CM2"),
    ("CD_M2_V10k_CM2", 2.0, 10000.0, 2.0, 4.0,
     "cold shell, KE mid-point: sustains the CSM phase longer"),
    ("CD_M2_V07k_CM2", 2.0, 7000.0, 2.0, 2.0,
     "cold shell, slower piston: weaker flash + later CSM arrival"),
    ("CD_M1_V10k_CM4", 1.0, 10000.0, 4.0, 2.0,
     "cold shell + heavier blanket: maximum trapping at low KE"),
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
            f"Cocoon batch-4 (cold shell): H/He shell (M={m} Msun, "
            f"v={v:.0f} km/s, KE={ke:.1f}e50 erg, T=2000 K) in immediate "
            f"contact with SN_A ejecta + detached outer CSM "
            f"(M={cm:.1f} Msun, 850 km/s, r=6e16-1.2e17 cm).  {why}.  "
            f"Ni56 disabled.  Fixes the 65 kK radiation-reservoir "
            f"artifact of batches 2-3."
        )
        x = dict(parent["x_model_from_sn"])
        x.update(
            {
                "shell_mass_msun": m,
                "shell_velocity_kms": v,
                "v0_kms": 0.3 * v,
                "outer_csm_mass_msun": cm,
                # Cold initial shell: in the evacuated cavity the hot
                # Dessart T-profile is a radiation bomb, not a minor
                # detail.  Shock provides all heating.
                "T_interface_K": 2000.0,
                "T_max_K": 2000.0,
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
            f.write("# Cocoon batch-4: internal shock + outer CSM, cold shell.\n")
            f.write(f"# {why}\n")
            f.write(f"# Parent config: {PARENT}\n")
            yaml.dump(cfg, f, default_flow_style=False, sort_keys=False)
        print("wrote", path)


if __name__ == "__main__":
    main()
