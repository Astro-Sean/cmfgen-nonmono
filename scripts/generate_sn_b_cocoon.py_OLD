#!/usr/bin/env python3
"""Generate the SN_B "cocoon / internal shock" test batch.

Physics: in every previous SN_B grid the cavity between the launched
shell and the 315-400 d old SN_A ejecta was filled with a steady wind
(mdot 1e-5..1e-1 Msun/yr).  The wind dominated the lightcurve
(L ~ 0.5*(Mdot/v_w)*v_s^3), producing a too-bright first plateau and a
second peak only in a narrow t_start/mdot corridor.  None reproduce the
observed SN_B morphology (r ~ -15.5 precursor, o = -18.54 peak, +2
mag/100d tail).

This batch tests the alternative geometry: an evacuated cavity and the
SN_B shell launched in IMMEDIATE CONTACT with the SN_A ejecta inner
edge (no_wind: true).  The shell piston drives a forward shock directly
into the metal-rich ejecta -- the old SN_A ejecta is the confining
"cocoon" and the shock propagates inside it (internal shock), giving
shock breakout/cooling through the ejecta surface rather than
shell-wind dissipation.  In the Khatami & Kasen (2024, ApJ 972, 140)
framework the SN_A ejecta plays the role of the CSM shell with
eta = M_ejecta/M_shell and breakout parameter xi set by the shock
speed vs the ejecta diffusion depth.

t_start is set to the cavity crossing time so that the simulation
starts at the moment of contact: t_c = v_in * source_age / (V - v_in),
with v_in ~ 1800 km/s the SN_A ejecta inner-edge velocity (homologous).
"""
import os
import sys

import yaml

BASE = "/viper/ptmp1/sbrennan/SN2024pba"
GRID = os.path.join(BASE, "SN_B_cocoon")
PARENT = os.path.join(
    BASE, "SN_B_double_peak_grid/SB_M10_ts8_W3e-4/config.yaml"
)

# SN_A ejecta inner-edge velocity at ~1 yr (homologous; ~1500-2200 km/s
# in the source .swd).  Used only for the self-consistent t_start.
V_IN_KMS = 1800.0


def t_contact(source_age, v_shell):
    return V_IN_KMS * source_age / (v_shell - V_IN_KMS)


# (name, M_shell, V_shell, source_age, rationale)
VARIANTS = [
    ("CO_M5_V10k_sa315", 5.0, 10000.0, 315.0,
     "base cocoon: adopted launch epoch (SN_A age 315 d)"),
    ("CO_M2_V10k_sa315", 2.0, 10000.0, 315.0,
     "light piston: eta=M_ej/M_sh~4, shock may stall inside cocoon"),
    ("CO_M10_V10k_sa315", 10.0, 10000.0, 315.0,
     "heavy piston: eta~0.9, shock plows through ejecta"),
    ("CO_M5_V15k_sa315", 5.0, 15000.0, 315.0,
     "stronger shock: earlier breakout, brighter peak"),
    ("CO_M5_V06k_sa315", 5.0, 6000.0, 315.0,
     "weak shock bracket: longer crossing, dimmer"),
    ("CO_M5_V10k_sa400", 5.0, 10000.0, 400.0,
     "age anchor matching the wind-filled grids (400 d)"),
]


def main():
    with open(PARENT) as f:
        parent = yaml.safe_load(f)

    for name, m, v, sa, why in VARIANTS:
        if len(name) > 28:
            sys.exit(f"name too long: {name} ({len(name)} chars)")
        ts = round(t_contact(sa, v))
        cfg = dict(parent)
        cfg["model_name"] = name
        cfg["description"] = (
            f"Cocoon/internal-shock SN_B test: H/He shell (M={m} Msun, "
            f"v={v:.0f} km/s, solar) in immediate contact with SN_A "
            f"ejecta inner edge (no cavity wind).  {why}.  "
            f"t_start={ts} d = cavity crossing time (contact at sim "
            f"start).  Ni56 disabled."
        )
        x = dict(parent["x_model_from_sn"])
        x.update(
            {
                "shell_mass_msun": m,
                "shell_velocity_kms": v,
                "v0_kms": 0.3 * v,
                "source_age_days": sa,
                "t_start_days": float(ts),
                "no_wind": True,
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
            f.write(f"# Cocoon/internal-shock SN_B test. {why}\n")
            f.write(f"# no_wind=True -> immediate contact: shell outer "
                    f"edge = 0.99 * ejecta inner edge; shock propagates "
                    f"through metal-rich SN_A ejecta (the cocoon).\n")
            f.write(f"# Parent config: {PARENT}\n")
            yaml.dump(cfg, f, default_flow_style=False, sort_keys=False)
        print("wrote", path, f"(t_start={ts}d)")


if __name__ == "__main__":
    main()
