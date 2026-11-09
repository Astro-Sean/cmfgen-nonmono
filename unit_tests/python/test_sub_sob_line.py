#!/usr/bin/env python3
"""Unit tests for SUB_SOB_LINE_V3 (new_main/mod_subs/sub_sob_line_v3.f).

The 05-Oct-2026 alterations are exercised through the LAMBDA_ITERATION
branch (COMPUTE_BA=T, LAMBDA_ITERATION=T), which reaches every guard
without pulling in the DO_VAR_CONT machinery:

  * SIGMA vs SIGMA_TRUE dispatch under DO_NONMONO_SOB /
    DO_SOB_CLAMPED_LINES: the clamped case must reproduce the flag-off
    baseline exactly, the unclamped case must differ.
  * DO_NONLOC_SOB -> SOB_RH_COUPLE: a velocity dip produces nonzero
    DRC_SIM which shifts the STEQ residual.
  * ETAL_MAT <= 0 guard: VB_SIM/VC_SIM are zeroed at that depth rather
    than carrying a wrong-sign derivative.
  * NEG_OPACITY: VB_SIM zeroed at the flagged depth.
  * POPS(NUP)=0: the expanded dRATE_dUP stays finite (the algebraic
    RATE_FAC*ZNET/POPS(NUP) form would NaN).
  * U_STAR_RATIO=0: the expanded dRATE_dT stays finite.
  * LST_ITERATION SOBEW dispatch (EW output).

DRC_SIM zero-denominator checks on BB_COR/CHIL/ETAL_MAT are pure
defensive code: an exact zero in any of them already produces
non-finite values inside SOBJBAR_SIM before the guard can run, so
those branches are unreachable from a live call and are not claimed
as covered.
"""
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import (
    BINDIR,
    F90DIR,
    PROD_MOD_INC,
    TREE,
    Report,
    compile_objs,
    link_prod,
    prod_obj,
    prod_objs,
    require_fresh_archives,
    run,
)

HELPERS = [
    "tools/set_kind_module.f",
]

# Repeat pass: libsubs members pull module objects from libmod and
# vice versa; the linker resolves archives in order.
LIBS = [
    TREE / "lib" / "libmod.a",
    TREE / "lib" / "libsubs.a",
    TREE / "lib" / "libmod.a",
    TREE / "lib" / "libsubs.a",
]

ND = 8


def build():
    require_fresh_archives()
    md = BINDIR / "obj_sub_sob_line"
    objs = prod_objs([TREE / h for h in HELPERS])
    objs += [prod_obj(TREE / "new_main" / "mod_subs"
                          / "sub_sob_line_v3.f")]
    objs += compile_objs(
        "drv_sub_sob_line",
        [F90DIR / "drv_sub_sob_line.f90"],
        moddir=md,
        extra_flags=PROD_MOD_INC,
    )
    return link_prod("drv_sub_sob_line", objs + LIBS, moddir=md)


def parse(out):
    vals = {"VB": {}, "V2": {}, "ZN": {}, "DR": {}, "ST": {},
            "BAU": {}, "BAT": {}, "NAN": None, "EW": None}
    for ln in out.splitlines():
        t = ln.split()
        if not t:
            continue
        if t[0] == "VB" and len(t) == 4:
            vals["VB"][int(t[1])] = (float(t[2]), float(t[3]))
        elif t[0] == "V2" and len(t) == 4:
            vals["V2"][int(t[1])] = (float(t[2]), float(t[3]))
        elif t[0] == "ZN" and len(t) == 4:
            vals["ZN"][int(t[1])] = (float(t[2]), float(t[3]))
        elif t[0] == "DR" and len(t) == 3:
            vals["DR"][int(t[1])] = float(t[2])
        elif t[0] == "ST" and len(t) == 5:
            vals["ST"][int(t[1])] = (float(t[2]), float(t[3]),
                                     float(t[4]))
        elif t[0] == "BAU" and len(t) == 6:
            vals["BAU"][int(t[1])] = tuple(float(x) for x in t[2:])
        elif t[0] == "BAT" and len(t) == 5:
            vals["BAT"][int(t[1])] = tuple(float(x) for x in t[2:])
        elif t[0] == "NAN":
            vals["NAN"] = int(t[1])
        elif t[0] == "EW":
            vals["EW"] = (float(t[1]), float(t[2]))
    return vals


def finite(x):
    return math.isfinite(x)


def main():
    rep = Report("SUB_SOB_LINE_V3 altered paths")
    binary = build()

    base = parse(run(binary, "base").stdout)
    rep.check("base: all outputs finite", base["NAN"] == 0)
    rep.check(
        "base: VB_2/VC_2 zeroed under LAMBDA_ITERATION",
        all(v == (0.0, 0.0) for v in base["V2"].values()),
    )

    # -- SIGMA dispatch --------------------------------------------------
    st = parse(run(binary, "sig_true").stdout)
    sc = parse(run(binary, "sig_clmp").stdout)
    rep.check("sig_true: all outputs finite", st["NAN"] == 0)
    rep.check(
        "sig_clmp equals flag-off baseline exactly",
        base["ST"] == sc["ST"] and base["VB"] == sc["VB"]
        and base["ZN"] == sc["ZN"],
    )
    differs = any(
        abs(st["ST"][k][c] - base["ST"][k][c]) > 0.0
        for k in st["ST"] for c in range(3)
    )
    rep.check(
        "sig_true routes SIGMA_TRUE into the line operator",
        differs,
    )

    # -- ETAL_MAT < 0 guard ----------------------------------------------
    ne = parse(run(binary, "neg_etal").stdout)
    rep.check("neg_etal: all outputs finite", ne["NAN"] == 0)
    rep.check(
        "neg_etal: VB_SIM/VC_SIM zeroed at ETAL_MAT<0 depth",
        ne["VB"][3] == (0.0, 0.0),
    )
    # Valid depths must satisfy the guard's analytic values:
    #   VB_SIM = -JBAR/ETAL_MAT, VC_SIM = JBAR*CHIL_MAT/ETAL_MAT^2.
    chil_v = (3.0e4, 4.5e4, 6.5e4, 9.5e4, 1.4e5, 2.1e5, 3.2e5, 5.0e5)
    etal_mat1 = [0.91 * c for c in chil_v]      # 0.7*1.3*CHIL
    chil_mat1 = [0.6 * c for c in chil_v]
    ok = all(
        abs(ne["VB"][i][0] + ne["ZN"][i][1] / etal_mat1[i - 1])
        < 1e-10 * max(1.0, abs(ne["VB"][i][0]))
        for i in (1, 2, 4, 5, 6, 7, 8)
    )
    rep.check("neg_etal: VB_SIM=-JBAR/ETAL_MAT at valid depths", ok)
    ok = all(
        abs(ne["VB"][i][1] - ne["ZN"][i][1] * chil_mat1[i - 1]
            / etal_mat1[i - 1] ** 2)
        < 1e-10 * max(1.0, abs(ne["VB"][i][1]))
        for i in (1, 2, 4, 5, 6, 7, 8)
    )
    rep.check("neg_etal: VC_SIM=JBAR*CHIL_MAT/ETAL_MAT^2", ok)

    # -- NEG_OPACITY -------------------------------------------------------
    no = parse(run(binary, "neg_opac").stdout)
    rep.check("neg_opac: all outputs finite", no["NAN"] == 0)
    rep.check(
        "neg_opac: VB_SIM zeroed at flagged depth",
        no["VB"][4][0] == 0.0 and no["VB"][4][1] != 0.0,
    )

    # -- Expanded dRATE forms: no NaN with zero population/ratio ---------
    zp = parse(run(binary, "zeropop").stdout)
    zr = parse(run(binary, "zerorat").stdout)
    rep.check("zeropop: all outputs finite", zp["NAN"] == 0)
    rep.check("zerorat: all outputs finite", zr["NAN"] == 0)
    rep.check(
        "zeropop: BA entries all finite",
        all(finite(x) for v in zp["BAU"].values() for x in v)
        and all(finite(x) for v in zp["BAT"].values() for x in v),
    )
    rep.check(
        "zerorat: BA entries all finite",
        all(finite(x) for v in zr["BAU"].values() for x in v)
        and all(finite(x) for v in zr["BAT"].values() for x in v),
    )

    # -- DO_NONLOC_SOB ---------------------------------------------------
    nl = parse(run(binary, "nonloc").stdout)
    rep.check("nonloc: all outputs finite", nl["NAN"] == 0)
    rep.check(
        "nonloc: DRC_SIM nonzero somewhere",
        any(v != 0.0 for v in nl["DR"].values()),
    )
    rep.check(
        "nonloc: DRC_SIM zero under flag-off baseline",
        all(v == 0.0 for v in base["DR"].values()),
    )
    shifted = any(
        nl["ST"][k][0] != base["ST"][k][0]
        or nl["ST"][k][1] != base["ST"][k][1]
        for k in nl["ST"]
    )
    rep.check("nonloc: remote term shifts the STEQ residual", shifted)

    # -- LST_ITERATION SOBEW dispatch -------------------------------------
    et = parse(run(binary, "ew_true").stdout)
    ec = parse(run(binary, "ew_clmp").stdout)
    rep.check("ew_true: EW output produced", et["EW"] is not None)
    rep.check("ew_clmp: EW output produced", ec["EW"] is not None)
    rep.check(
        "ew: unclamped EW differs from clamped",
        et["EW"] is not None and ec["EW"] is not None
        and et["EW"][0] != ec["EW"][0],
    )

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
