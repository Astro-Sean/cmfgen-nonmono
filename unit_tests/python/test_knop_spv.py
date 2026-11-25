#!/usr/bin/env python3
"""Runtime ABI/purity test for KN_ORD_SPV_IMPL.

Phase-5.5 follow-up (hardening): the OPTIONAL-argument design was
replaced by a module-scope worker (KN_ORD_SPV_IMPL, required
LOGICAL L_CNT) behind two explicit-interface entry points:
KN_ORD_SPV (production, .TRUE.) and KN_ORD_SPV_SHADOW (shadow,
.FALSE.).

Coverage is reported in four distinct tiers:

  RUNTIME WORKER  - drv_knop_spv.f90 USEs FG_J_CMF_NONMONO_MOD_V1
                    and calls the compiled worker with
                    L_CNT=.TRUE./.FALSE. across an 11-case branch
                    matrix: GAM=0, expanding L>1 and L=1,
                    compressing with two / JPR1-only / JPR2-only /
                    no valid donors, compressing L+NNF2>NCF
                    boundary, FG_FREQ_MATCH on (NNF2=2 and NNF2=1),
                    and the unallocated KN_DBG_NDON path.
  INDEPENDENT     - every case's expected SPV is recomputed here
                    ARITHMETIC        from the manufactured inputs
                    using the published branch formulas
                    (including the production NNF2 definition) and
                    compared at 1e-12.
  OPERAND         - every input array is snapshotted before
                    INTEGRITY         EACH worker call and compared
                    elementwise after it (NWS/NWP per case).  GAM
                    is checked in its branch-configured state -
                    the driver sets it, the worker must never
                    write it.  KN_DBG_NDON is the one permitted
                    mutable input: it is diffed elementwise and
                    the single +1 at (2,3) is verified explicitly.
                    Matrix-level shape+SUM checksums remain as a
                    secondary guard.
  WRAPPER SOURCE  - the two thin wrappers are internal to
                    PINS              KNOP_RAY_SOLVE_V1 and cannot
                    be called without the full ray context, so
                    their forwarding is pinned at source level
                    (single IMPL call, correct args, .TRUE./.FALSE.).
                    The wrappers are NOT runtime-covered here.
"""

import math
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import (
    F90DIR,
    TREE,
    Report,
    compile_and_link_prod,
    run,
)

HELPERS = [
    "tools/set_kind_module.f",
    "unix/term_io.f",
    "unix/gen_asci_open.f",
    "unix/verbose_info.f",
    "subs/set_case.f",
    "subs/phys_con.f",
    "subs/mon_int_funs_v2.f",
    "subs/mon_interp.f",
    "subs/mod_trap_derivatives.f",
    "subs/nordtau.f",
    "subs/qfromf.f",
    "subs/integrate.f",
    "subs/write_vec.f",
    "unix/dir_acc_pars_gen.f",
    "plane/mod_ray_mom_store.f",
    "new_main/mod_subs/knop_status_mod.f90",
    "plane/get_moms_non_rel.f",
    "new_main/mod_subs/mod_cmfgen.f",
    "new_main/mod_subs/control_variable_mod.f",
    "new_main/mod_subs/jfnk_cmfgen_int_mod.f",
    "new_main/mod_subs/jfnk_solver_v1.f90",
    "new_main/mod_subs/check_param_consistency.f",
    "newsubs/rd_store_log.f",
    "new_main/mod_subs/iter_ledger_mod.f90",
    "new_main/mod_subs/iter_log_fmt_mod.f90",
    "plane/fg_j_cmf_nonmono_v1.f",
    "plane/mom_j_cmf_nonmono_v1.f",
    "unit_tests/f90/stub_vadat.f",
]

FG = (TREE / "plane/fg_j_cmf_nonmono_v1.f").read_text()
NCFX = 6


def manufactured(r, l):
    """Reproduce the driver's manufactured operand values."""
    sr = 1.5 + 0.01 * r + 0.001 * l
    q = 0.3 + 0.01 * r
    fld = 0.1 * r + 0.01 * l
    drf = 0.4 + 0.1 * r
    winx = 0.25
    jpr1 = max(1, r - 1)
    jpr2 = r
    return sr, q, fld, drf, winx, jpr1, jpr2


def expected(name, r, l):
    """Independent expected SPV from the branch definitions."""
    sr, q, fld, drf, winx, jpr1, jpr2 = manufactured(r, l)
    if name in ("gam0", "gamp1"):
        return sr
    if name in ("gamp", "unall"):
        fld_prev = 0.1 * r + 0.01 * (l - 1)
        return sr + q * fld_prev
    if name == "gamnn":
        return sr
    if name == "gamnb":
        return sr
    # Compressing with donors: NNF2 per the production definition.
    fm = name in ("fmtch", "fmtcb")
    nnf2 = 1
    if fm and l < NCFX:
        kn_dln = 0.2
        t1 = drf / kn_dln
        if t1 >= NCFX - l:
            nnf2 = NCFX - l
        else:
            nnf2 = max(1, math.ceil(t1))
    if l + nnf2 > NCFX:
        return sr
    slot = l + nnf2
    f1 = 0.1 * jpr1 + 0.01 * slot
    f2 = 0.1 * jpr2 + 0.01 * slot
    if name == "gamn2":
        inext = (1.0 - winx) * f1 + winx * f2
    elif name == "gamn1":
        inext = f1
    elif name == "gamnj2":
        inext = f2
    else:  # fmtch, fmtcb: both donors valid
        inext = (1.0 - winx) * f1 + winx * f2
    return sr + abs(q) * inext


EXP_NDP = {
    "gam0": 0,
    "gamp": 1,
    "gamp1": 0,
    "gamn2": 0,
    "gamn1": 0,
    "gamnj2": 0,
    "gamnn": 0,
    "gamnb": 0,
    "fmtch": 0,
    "fmtcb": 0,
    "unall": -1,
}

# Arrays whose PRE/POST checksum must be identical.  GAM is
# covered by the per-invocation elementwise compare (the driver
# deliberately changes it to configure each branch).  NDON is
# handled separately (one expected +1 delta at (2,3)).
OPERANDS = [
    "SR", "Q", "FLD", "DRF", "WINX", "JPR1", "JPR2",
    "IKC1", "IKC2", "DLN",
]


def build():
    return compile_and_link_prod(
        "knspv", [TREE / h for h in HELPERS] + [F90DIR / "drv_knop_spv.f90"]
    )


def main():
    rep = Report("knop_spv")
    exe = build()
    out = run(exe).stdout
    print(out)

    # ---- RUNTIME WORKER: parse the branch matrix ----
    cases = {}
    for m in re.finditer(
        r"SPCASE\s+(\S+)\s+R\s*(\d+)\s+L\s*(\d+)\s+SVT\s*(\S+)\s+"
        r"SVF\s*(\S+)\s+BEQ\s*([TF])\s+NDP\s*(-?\d+)\s+NDS\s*(-?\d+)"
        r"\s+NWS\s*(-?\d+)\s+NWP\s*(-?\d+)",
        out,
    ):
        cases[m.group(1).strip()] = {
            "R": int(m.group(2)),
            "L": int(m.group(3)),
            "SVT": float(m.group(4)),
            "SVF": float(m.group(5)),
            "BEQ": m.group(6) == "T",
            "NDP": int(m.group(7)),
            "NDS": int(m.group(8)),
            "NWS": int(m.group(9)),
            "NWP": int(m.group(10)),
        }

    for name, ndp in EXP_NDP.items():
        c = cases.get(name)
        rep.check(f"{name} parsed", c is not None)
        if c is None:
            continue
        rep.check(f"{name} SPV bit-identical", c["BEQ"])
        rep.check(f"{name} NDP={ndp}", c["NDP"] == ndp, str(c["NDP"]))
        exp_nds = -1 if name == "unall" else 0
        rep.check(f"{name} NDS={exp_nds}", c["NDS"] == exp_nds)
        rep.check(
            f"{name} inputs unchanged shadow call", c["NWS"] == 0,
        )
        rep.check(
            f"{name} inputs unchanged prod call", c["NWP"] == 0,
        )

    # ---- INDEPENDENT ARITHMETIC ----
    for name, c in cases.items():
        want = expected(name, c["R"], c["L"])
        rep.check(
            f"{name} independent value",
            abs(c["SVT"] - want) < 1e-12,
            f"SVT={c['SVT']} want={want}",
        )

    # ---- OPERAND INTEGRITY ----
    pre, post = {}, {}
    for m in re.finditer(
        r"(PRE|POST)\s+(\S+)\s+(\d+)\s+(\d+)\s+SUM\s*(\S+)", out
    ):
        (pre if m.group(1) == "PRE" else post)[m.group(2)] = (
            int(m.group(3)),
            int(m.group(4)),
            m.group(5),
        )
    for op in OPERANDS:
        rep.check(
            f"{op} operand unchanged",
            op in pre and op in post and pre[op] == post[op],
            f"pre={pre.get(op)} post={post.get(op)}",
        )

    m = re.search(r"NDCHG\s+(\d+)", out)
    ndchg = int(m.group(1)) if m else -1
    deltas = re.findall(r"NDELT\s+(\d+)\s+(\d+)\s+(-?\d+)", out)
    rep.check(
        "NDON exactly one change: (2,3) +1",
        ndchg == 1
        and len(deltas) == 1
        and (int(deltas[0][0]), int(deltas[0][1]), int(deltas[0][2]))
        == (2, 3, 1),
        f"ndchg={ndchg} deltas={deltas}",
    )

    # ---- WRAPPER SOURCE PINS (not runtime coverage) ----
    prod_wrap = FG[
        FG.index("SUBROUTINE KN_ORD_SPV(FLD,R,L,SPV)"):
    ].split("\n\tEND")[0]
    shad_wrap = FG[
        FG.index("SUBROUTINE KN_ORD_SPV_SHADOW(FLD,R,L,SPV)"):
    ].split("\n\tEND")[0]
    rep.check(
        "prod wrapper: 4-arg, single IMPL call, .TRUE.",
        "SUBROUTINE KN_ORD_SPV(FLD,R,L,SPV)" in prod_wrap
        and prod_wrap.count("CALL KN_ORD_SPV_IMPL(") == 1
        and "PK_SR,PK_Q,DRF_L,JPR1,JPR2,WINX," in prod_wrap
        and "IKC1,IKC2,LS,NCF,FLD,R,L,SPV,.TRUE.)" in prod_wrap,
    )
    rep.check(
        "shadow wrapper: 4-arg, single IMPL call, .FALSE.",
        "SUBROUTINE KN_ORD_SPV_SHADOW(FLD,R,L,SPV)" in shad_wrap
        and shad_wrap.count("CALL KN_ORD_SPV_IMPL(") == 1
        and "PK_SR,PK_Q,DRF_L,JPR1,JPR2,WINX," in shad_wrap
        and "IKC1,IKC2,LS,NCF,FLD,R,L,SPV,.FALSE.)" in shad_wrap,
    )
    rep.check(
        "no OPTIONAL in any SPV routine",
        "OPTIONAL" not in prod_wrap
        and "OPTIONAL" not in shad_wrap
        and "OPTIONAL" not in FG[
            FG.index("SUBROUTINE KN_ORD_SPV_IMPL"):
            FG.index("END SUBROUTINE KN_ORD_SPV_IMPL")
        ],
    )
    rep.check(
        "call sites: shadow eval uses shadow entry",
        "CALL KN_ORD_SPV_SHADOW(FSF,RS,L,SV(KR))" in FG,
    )
    rep.check(
        "call sites: production keeps 4-arg KN_ORD_SPV",
        "CALL KN_ORD_SPV(FLD,R,L2,T1)" in FG
        and FG.count("CALL KN_ORD_SPV(FLD,J,L,SPJ)") == 2,
    )

    nfail = rep.summary()
    if nfail:
        sys.exit(1)
    print("all KNOP SPV ABI/purity checks passed")


if __name__ == "__main__":
    main()
