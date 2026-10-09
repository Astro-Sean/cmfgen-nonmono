#!/usr/bin/env python3
"""Generate batch 2 of the SN_B cocoon/internal-shock tests.

Batch 1 (SN_B_cocoon, CO_*) showed that a shell in immediate contact
with the SN_A ejecta always produces a ~-23..-24 prompt flash: the
retained ejecta above the contact (~1.6 Msun, tau ~ 0.3 at contact age)
cannot trap the shock radiation, and the shock recedes in mass
coordinate as the whole ejecta is accelerated to ~shell speed within
days.

Batch 2 keeps the same cocoon geometry but adds the dense, slow CSM
that the observed narrow P-Cygni (~850 km/s) independently requires.
The CSM is a detached shell ABOVE the SN_A ejecta outer edge (a pre-SN
eruption at ~850 km/s; SN_A's 14500 km/s leading edge reaches only
~4.8e16 cm in 384 d, so a shell at 6e16+ cm is untouched until
SN_B-era material arrives).  The internal shock crosses the metal
cocoon first (mostly hidden phase -> faint early emission), then the
accelerated ejecta runs into the CSM, whose optical depth
(tau ~ 15 per Msun for r_in=6e16, dr=6e16) traps and reprocesses the
prompt flash and powers the observed -18.5 peak with a rise set by
the outflow arrival-time spread.
"""
import os
import sys

import yaml

BASE = "/viper/ptmp1/sbrennan/SN2024pba"
GRID = os.path.join(BASE, "SN_B_cocoon")
PARENT = os.path.join(BASE, "SN_B_cocoon/CO_M5_V10k_sa315/config.yaml")

# Detached CSM shell at ~850 km/s (narrow P-Cygni speed).
CSM_R_IN = 6.0e16
CSM_R_OUT = 1.2e17
CSM_V = 850.0

# (name, M_shell, V_shell, csm_mass_msun, rationale)
VARIANTS = [
    ("CB_M5_V10k_CM2", 5.0, 10000.0, 2.0,
     "base: batch-1 geometry + 2 Msun detached CSM (tau~30)"),
    ("CB_M2_V10k_CM2", 2.0, 10000.0, 2.0,
     "lighter piston: less KE, weaker sustained phase"),
    ("CB_M5_V06k_CM2", 5.0, 6000.0, 2.0,
     "weaker internal shock: dimmer prompt leakage"),
    ("CB_M5_V10k_CM4", 5.0, 10000.0, 4.0,
     "denser CSM (tau~60): more trapping, later/brighter peak"),
    ("CB_M5_V10k_CM1", 5.0, 10000.0, 1.0,
     "thinner CSM (tau~15): flash partially leaks"),
]


def main():
    with open(PARENT) as f:
        parent = yaml.safe_load(f)

    for name, m, v, cm, why in VARIANTS:
        if len(name) > 28:
            sys.exit(f"name too long: {name} ({len(name)} chars)")
        cfg = dict(parent)
        cfg["model_name"] = name
        cfg["description"] = (
            f"Cocoon batch-2: H/He shell (M={m} Msun, v={v:.0f} km/s) in "
            f"immediate contact with SN_A ejecta (internal shock through "
            f"the metal cocoon) + detached outer CSM shell "
            f"(M={cm:.1f} Msun, v={CSM_V:.0f} km/s, "
            f"r={CSM_R_IN:.1e}-{CSM_R_OUT:.1e} cm).  {why}.  "
            f"Ni56 disabled.  Full ejecta retained (no v_max cut)."
        )
        x = dict(parent["x_model_from_sn"])
        x.update(
            {
                "shell_mass_msun": m,
                "shell_velocity_kms": v,
                "v0_kms": 0.3 * v,
                "ejecta_v_max_kms": 0.0,
                "outer_csm_mass_msun": cm,
                "outer_csm_velocity_kms": CSM_V,
                "outer_csm_r_inner_cm": CSM_R_IN,
                "outer_csm_r_outer_cm": CSM_R_OUT,
                "outer_csm_T_K": 2000.0,
                "outer_csm_zones": 60,
            }
        )
        cfg["x_model_from_sn"] = x

        stella = dict(parent["stella"])
        stella["target_days"] = 150.0
        stella["nstmax"] = 4000000
        cfg["stella"] = stella

        out = dict(parent["output"])
        out["model_dir"] = f"./SN_B_cocoon/{name}"
        cfg["output"] = out

        mdir = os.path.join(GRID, name)
        os.makedirs(mdir, exist_ok=True)
        path = os.path.join(mdir, "config.yaml")
        with open(path, "w") as f:
            f.write("# Cocoon batch-2: internal shock + detached outer CSM.\n")
            f.write(f"# {why}\n")
            f.write(f"# Parent config: {PARENT}\n")
            yaml.dump(cfg, f, default_flow_style=False, sort_keys=False)
        print("wrote", path)


if __name__ == "__main__":
    main()
