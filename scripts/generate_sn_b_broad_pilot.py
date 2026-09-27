#!/usr/bin/env python3
"""Generate the 10-model pilot subset of SN_B_broad_grid.

New folder SN_B_broad_pilot/ with fixed-width names (all 21 chars):

    SB_M10_V10k_n08_W3e-3
         M{M:02g}   shell mass, Msun
         V{V:02g}k  shell velocity, 1000 km/s
         n{n:02g}   envelope density exponent
         W{mdot}    wind mass-loss rate tag

Pilot selection: OFAT cross on the incumbent anchor
SB_M10_V10k_n08_W3e-3 (D=3.5e-6 ~ parent D=3e-6 at v_w=860), one
endpoint per axis, plus mdot extremes for leverage.

All models: v_wind=860 km/s, t_start=8 d, v0=0.3*V_shell, 180 d STELLA,
SN_A fixed (SN_A_g_T20k_Ni06 @ 400 d).
"""
import math
import os
import sys

import yaml

BASE = "/viper/ptmp1/sbrennan/SN2024pba"
GRID = os.path.join(BASE, "SN_B_broad_pilot")
PARENT = os.path.join(BASE, "SN_B_double_peak_grid",
                      "SB_M10_ts8_W3e-4", "config.yaml")

V_WIND = 860.0
T_START = 8.0
TARGET_DAYS = 180.0
V0_FRAC = 0.3

# (M, V, n_rho, mdot, why)
PILOT = [
    (10.0, 10000.0, 8.0, 3e-3,
     "anchor: D=3.5e-6 matches incumbent parent at v_w=860"),
    (3.0, 10000.0, 8.0, 3e-3, "low shell mass endpoint"),
    (12.0, 10000.0, 8.0, 3e-3, "high shell mass endpoint"),
    (10.0, 5000.0, 8.0, 3e-3, "low shell velocity endpoint"),
    (10.0, 15000.0, 8.0, 3e-3, "high shell velocity endpoint"),
    (10.0, 10000.0, 6.0, 3e-3, "shallow envelope endpoint"),
    (10.0, 10000.0, 12.0, 3e-3, "steep envelope endpoint"),
    (10.0, 10000.0, 8.0, 1e-4, "thin wind (D=1.2e-6)"),
    (10.0, 10000.0, 8.0, 1e-2, "dense wind (D=1.2e-5)"),
    (10.0, 10000.0, 8.0, 1e-1, "extreme dense wind (D=1.2e-4)"),
]


def tag_mdot(m):
    exp = math.floor(math.log10(m) + 1e-9)
    mant = m / 10.0**exp
    if abs(mant - 1.0) < 1e-9:
        return f"W1e{exp}"
    return f"W{mant:g}e{exp}"


def name_of(M, V, n, md):
    return (f"SB_M{M:02g}_V{V / 1000.0:02g}k_n{n:02g}_{tag_mdot(md)}")


def main():
    with open(PARENT) as f:
        parent = yaml.safe_load(f)
    for M, V, n, md, why in PILOT:
        name = name_of(M, V, n, md)
        assert len(name) == 21, (name, len(name))
        cfg = dict(parent)
        cfg["model_name"] = name
        cfg["description"] = (
            f"Broad pilot: M_shell={M:g} Msun, V_shell={V:g} km/s, "
            f"n_rho={n:g}, mdot={md:g} Msun/yr, v_w={V_WIND:g} km/s, "
            f"t_start={T_START:g} d, v0={V0_FRAC * V:g} km/s; {why}. "
            f"SN_A fixed (SN_A_g_T20k_Ni06 @ 400 d), STELLA "
            f"{TARGET_DAYS:g} d.")
        x = dict(parent["x_model_from_sn"])
        x.update(shell_mass_msun=M, shell_velocity_kms=V,
                 v0_kms=V0_FRAC * V, n_rho=n, wind_mdot_msun_yr=md,
                 wind_velocity_kms=V_WIND)
        cfg["x_model_from_sn"] = x
        st = dict(parent["stella"])
        st.update(target_days=TARGET_DAYS, nstmax=4000000)
        cfg["stella"] = st
        out = dict(parent["output"])
        out["model_dir"] = f"./SN_B_broad_pilot/{name}"
        cfg["output"] = out
        mdir = os.path.join(GRID, name)
        os.makedirs(mdir, exist_ok=True)
        path = os.path.join(mdir, "config.yaml")
        with open(path, "w") as f:
            f.write(f"# SN_B_broad_pilot: {why}\n# Parent: {PARENT}\n")
            yaml.dump(cfg, f, default_flow_style=False, sort_keys=False)
        print("wrote", path)


if __name__ == "__main__":
    main()
