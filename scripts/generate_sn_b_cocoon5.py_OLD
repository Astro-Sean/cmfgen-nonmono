#!/usr/bin/env python3
"""Generate batch 5 of the SN_B cocoon/internal-shock tests.

Batch 4 (CD_*: cold 2000 K shell) removed the radiation-reservoir
artifact: E_rad fell from ~5e52 to ~2.5e48 erg -- the lightcurve is now
honestly powered by the interaction.  But all CD models sit at a flat
Mo ~ -14.5: the shock into the thin detached CSM (tau~30, rho~5e-19)
is too weakly dissipative to feed the photosphere, so the emitting
surface stays deep inside at ~7e14 cm / 5000 K.

To reach the observed -18.5 the photosphere must sit in the CSM at
~1e17 cm / ~5000-6000 K for tens of days.  That needs a CSM that is
optically thick enough to become the photosphere but thin enough that
the shock-deposited energy diffuses out on ~weeks, not years --
tau ~ 5-10 rather than 30-60.  Batch 5 therefore sweeps CSM mass
0.5-2 Msun and placement (detached 6e16-1.2e17 vs attached-dense
5.5e16-9e16), all with the cold shell.
"""
import os
import sys

import yaml

BASE = "/viper/ptmp1/sbrennan/SN2024pba"
GRID = os.path.join(BASE, "SN_B_cocoon")
PARENT = os.path.join(BASE, "SN_B_cocoon/CD_M1_V10k_CM2/config.yaml")

# (name, M_shell, V_shell, M_CSM, r_in, r_out, rationale)
VARIANTS = [
    ("CE_M1_V10k_CM05", 1.0, 10000.0, 0.5, 6.0e16, 1.2e17,
     "thinnest CSM tau~3: shock radiates promptly"),
    ("CE_M1_V10k_CM1", 1.0, 10000.0, 1.0, 6.0e16, 1.2e17,
     "tau~5 detached blanket"),
    ("CE_M2_V10k_CM1", 2.0, 10000.0, 1.0, 5.5e16, 9.0e16,
     "more arriving KE + denser attached CSM tau~10"),
    ("CE_M2_V10k_CM2", 2.0, 10000.0, 2.0, 5.5e16, 9.0e16,
     "attached 2 Msun: thicker photosphere, slower release"),
    ("CE_M5_V10k_CM2", 5.0, 10000.0, 2.0, 6.0e16, 1.2e17,
     "cold version of CB flagship; KE=1e51 reference"),
]


def main():
    with open(PARENT) as f:
        parent = yaml.safe_load(f)

    for name, m, v, cm, rin, rout, why in VARIANTS:
        if len(name) > 28:
            sys.exit(f"name too long: {name} ({len(name)} chars)")
        cfg = dict(parent)
        cfg["model_name"] = name
        cfg["description"] = (
            f"Cocoon batch-5 (cold shell + moderate-tau CSM): "
            f"H/He shell (M={m} Msun, v={v:.0f} km/s, T=2000 K) in "
            f"immediate contact with SN_A ejecta + outer CSM "
            f"(M={cm:.2f} Msun, 850 km/s, r=[{rin:.1e},{rout:.1e}] cm).  "
            f"{why}.  Ni56 disabled."
        )
        x = dict(parent["x_model_from_sn"])
        x.update(
            {
                "shell_mass_msun": m,
                "shell_velocity_kms": v,
                "v0_kms": 0.3 * v,
                "outer_csm_mass_msun": cm,
                "outer_csm_r_inner_cm": rin,
                "outer_csm_r_outer_cm": rout,
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
            f.write("# Cocoon batch-5: cold shell + moderate-tau CSM.\n")
            f.write(f"# {why}\n")
            f.write(f"# Parent config: {PARENT}\n")
            yaml.dump(cfg, f, default_flow_style=False, sort_keys=False)
        print("wrote", path)


if __name__ == "__main__":
    main()
