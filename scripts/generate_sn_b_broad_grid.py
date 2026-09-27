#!/usr/bin/env python3
"""Generate the broad SN_B parameter grid (SN_B_broad_grid/).

Axes (single model per grid point; everything else inherited from the
SB_M10_ts8_W3e-4 parent, SN_A ejecta fixed at SN_A_g_T20k_Ni06 @ 400 d):

    M_shell   : 3, 5, 7, 10, 12 Msun
    V_shell   : 5000, 7000, 10000, 15000 km/s
    n_rho     : 6, 8, 10, 12          (shell density exponent)
    wind mdot : 1e-5 .. 1e-1          (D-stratified: 5 points, 1 dex steps)
    v_wind    : 860 km/s (fixed; observed narrow P-Cygni CSM velocity)
    t_start   : 8 d (fixed)
    duration  : 180 d STELLA (post-peak decay measurable)

v0_kms is set to 0.3 * V_shell so the flat-core fraction of the profile
is preserved for every shell velocity (a fixed v0=3000 would violate
v0 < v_shell for V=5000 and distort the shape for the others).

Best-fit anchor: the parent has D = mdot/v_w = 3e-4/100 = 3e-6.  At
v_w=860 the nearest grid point is mdot=3e-3 (D=3.5e-6), so
SB_M10_V10k_n08_W3e-3 is the rough incumbent match.

Naming (fixed width, user directive 22-Sep-2026):
    SB_M<M:02>_V<V:02>k_n<n:02>_W<mdot>   e.g. SB_M10_V07k_n08_W3e-3

NOTE: re-running this generator after the padding change creates NEW
zero-padded directories alongside the legacy unpadded ones
(SB_M10_V10k_n8_W3e-3 vs SB_M10_V10k_n08_W3e-3) -- do not regenerate
over the existing grid or the same physics will be submitted twice
under different names.
"""
import os
import sys

import yaml

BASE = "/viper/ptmp1/sbrennan/SN2024pba"
GRID = os.path.join(BASE, "SN_B_broad_grid")
PARENT = os.path.join(BASE, "SN_B_double_peak_grid",
                      "SB_M10_ts8_W3e-4", "config.yaml")

M_SHELL = [3.0, 5.0, 7.0, 10.0, 12.0]
V_SHELL = [5000.0, 7000.0, 10000.0, 15000.0]
N_RHO = [6.0, 8.0, 10.0, 12.0]
MDOT = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1]
V_WIND = 860.0
T_START = 8.0
TARGET_DAYS = 180.0
V0_FRAC = 0.3


def tag_mdot(m):
    # 1e-5 -> W1e-5 ; keep one mantissa digit for non-decade values
    import math
    exp = math.floor(math.log10(m) + 1e-9)
    mant = m / 10.0**exp
    if abs(mant - 1.0) < 1e-9:
        return f"W1e{exp}"
    return f"W{mant:g}e{exp}"


def tag_v(v):
    # fixed width (user directive 22-Sep-2026): V05k/V07k/V10k/V15k
    return f"V{v / 1000.0:02g}k"


def tag_n(n):
    return f"n{n:02g}"


def tag_m(m):
    return f"M{m:02g}"


def main():
    with open(PARENT) as f:
        parent = yaml.safe_load(f)

    # grid points: (M, V, n_rho, mdot)
    points = [(M, V, n, md) for M in M_SHELL for V in V_SHELL
              for n in N_RHO for md in MDOT]
    # incumbent anchor: D=3e-3/860=3.5e-6 ~ parent D=3e-6; not in the
    # decade-step mdot axis, added so the grid contains a best-fit match
    points.append((10.0, 10000.0, 8.0, 3e-3))
    names = []
    for M, V, n, md in points:
        name = (f"SB_{tag_m(M)}_{tag_v(V)}_{tag_n(n)}"
                f"_{tag_mdot(md)}")
        if len(name) > 28:
            sys.exit(f"name too long: {name}")
        names.append(name)
        cfg = dict(parent)
        cfg["model_name"] = name
        cfg["description"] = (
            f"Broad grid: M_shell={M:g} Msun, "
            f"V_shell={V:g} km/s, n_rho={n:g}, "
            f"wind mdot={md:g} Msun/yr at v_w={V_WIND:g} "
            f"km/s, t_start={T_START:g} d, "
            f"v0={V0_FRAC * V:g} km/s (0.3*V). "
            f"SN_A source fixed (SN_A_g_T20k_Ni06 @ 400 d). "
            f"STELLA to {TARGET_DAYS:g} d.")
        x = dict(parent["x_model_from_sn"])
        x["shell_mass_msun"] = M
        x["shell_velocity_kms"] = V
        x["v0_kms"] = V0_FRAC * V
        x["n_rho"] = n
        x["wind_mdot_msun_yr"] = md
        x["wind_velocity_kms"] = V_WIND
        cfg["x_model_from_sn"] = x
        st = dict(parent["stella"])
        st["target_days"] = TARGET_DAYS
        st["nstmax"] = 4000000
        cfg["stella"] = st
        out = dict(parent["output"])
        out["model_dir"] = f"./SN_B_broad_grid/{name}"
        cfg["output"] = out

        mdir = os.path.join(GRID, name)
        os.makedirs(mdir, exist_ok=True)
        path = os.path.join(mdir, "config.yaml")
        with open(path, "w") as f:
            f.write(
                f"# SN_B_broad_grid: M={M:g} V={V:g} "
                f"n_rho={n:g} mdot={md:g} v_w={V_WIND:g}\n")
            f.write(f"# Parent: {PARENT}\n")
            yaml.dump(cfg, f, default_flow_style=False,
                      sort_keys=False)
    print(f"wrote {len(names)} configs to {GRID}")
    print("example:", sorted(names)[len(names) // 2])
    anchor = "SB_M10_V10k_n8_W3e-3"
    print(f"incumbent-matching point: {anchor} (D={3e-3 / 860:.2e} vs "
          f"parent 3e-6)")


if __name__ == "__main__":
    main()
