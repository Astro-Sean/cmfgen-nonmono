#!/usr/bin/env python3
"""FG_J_CMF_V12 intensity-repair guards
(plane/fg_j_cmf_v12.f vs the frozen reference).

Alterations under test, compared directly against the pristine
reference source built with production flags.  All sit on the INTEGRAL
solve path inside the OpenMP ray loop:

1. Inner-boundary I_CORE floor.  With DIF on, the boundary intensity is
   (ETA_RAY + DBB*mu_z)/CHI_RAY; when the inward diffusion term DBB
   dominates ETA_RAY on an optically thin base the seed is negative.
   Prod floors it at 1e-30; upstream hands the negative seed to the
   outward sweep.  PROVEN: case 'neg' (DBB=-1e6) makes I_CORE about
   -3e13 - ref propagates it (I_M=-3.3e13 in its own diagnostic dump)
   and dies before finishing, while prod completes with JNU>0 and
   0<FEDD<=1 at every depth.

2. Negative I_M/I_P sweep floor (1e-30) with a 'negative I_X repaired
   on ray' diagnostic.  NOT triggered here: on the INIT=T INTEGRAL
   path with positive sources the monotonic-clamped -dS corrections
   can at most reach ~1/6 of the positive source terms, so the sweep
   cannot undershoot in a toy atmosphere.  The guard exists for
   iterated (Q/=0) and model-noise states.  Verified NON-SPURIOUS:
   the steep-source 'jump', opaque-wall 'wall' and source-free 'thin'
   profiles all give prod==ref bit-for-bit - the repair never fires
   where it should not.

3. JNU_STORE noise-floor repair before FEDD=K/J.  Same status as (2):
   requires quadrature noise on an underflow-floor intensity field.

The driver builds a 12-deep expanding atmosphere, IMPAR impact
parameters (6 core + 12 tangent rays) and the production quadrature
weights, then calls the solver once with INIT=T, NEW_FREQ=T, DIF=T,
METHOD='ZERO' (trapezoid dTAU) on the INTEGRAL path.  THK=F keeps the
fixed ray grid.
"""
import re
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import F90DIR, REF, TREE, Report, compile_and_link_prod, run

DRV = F90DIR / "drv_fg12.f90"
# Helpers the solver calls.  All but mon_int_funs_v2 are identical to
# the frozen reference, so production objects serve both binaries.
HELPERS = [
    TREE / "subs" / "genangqw_v2.f",
    TREE / "subs" / "nordangqw_v2.f",
    TREE / "subs" / "jweight_v2.f",
    TREE / "subs" / "hweight_v2.f",
    TREE / "subs" / "kweight_v2.f",
    TREE / "subs" / "nweight_v2.f",
    TREE / "subs" / "tuvghd_rh.f",
    TREE / "subs" / "write_direct_info_v3.f",
    TREE / "unix" / "dir_acc_pars_gen.f",
    TREE / "unix" / "term_io.f",
    TREE / "unix" / "tune.f",
    TREE / "subs" / "set_acc_dmu.f",
    TREE / "tools" / "chk_little_endian.f",
    TREE / "tools" / "get_lu.f",
]
PROD_SOURCES = [
    TREE / "plane" / "mod_ray_mom_store.f",
    TREE / "plane" / "fg_j_cmf_v12.f",
    TREE / "subs" / "mon_int_funs_v2.f",
]
REF_SOURCES = [
    # REF mod first: its GMNR_NMJ_G-free interface must shadow the
    # production .mod when the REF solver body compiles.
    REF / "plane" / "mod_ray_mom_store.f",
    REF / "plane" / "fg_j_cmf_v12.f",
    # mon_int_funs_v2 carries its own alteration; the ref binary gets
    # the ref version so the comparison stays strictly before/after.
    REF / "subs" / "mon_int_funs_v2.f",
]


def floats(tag, out):
    m = re.search(tag + r"=\s*((?:[-+0-9.Ee]+\s*)+)", out)
    return [float(x) for x in m.group(1).split()]


def main():
    rep = Report("FG_J_CMF_V12 intensity repair")
    prod = compile_and_link_prod(
        "drv_fg12", PROD_SOURCES + HELPERS + [DRV])
    ref = compile_and_link_prod(
        "drv_fg12_ref", REF_SOURCES + HELPERS + [DRV])

    with tempfile.TemporaryDirectory() as d:
        p = run(prod, "ctl", cwd=d).stdout
        r = run(ref, "ctl", cwd=d).stdout
        rep.check("ctl: both binaries complete",
                  "DONE_FG12" in p and "DONE_FG12" in r)
        pj, rj = floats("JNU", p), floats("JNU", r)
        pf, rf = floats("FEDD", p), floats("FEDD", r)
        rep.check("ctl: JNU finite and positive (prod)",
                  all(x > 0.0 for x in pj))
        rep.check("ctl: FEDD in (0,1] (prod)",
                  all(0.0 < x <= 1.0 + 1e-8 for x in pf))
        rep.check("ctl: prod JNU identical to ref (guards inert)",
                  pj == rj)
        rep.check("ctl: prod FEDD identical to ref", pf == rf)

        p = run(prod, "neg", cwd=d).stdout
        r = run(ref, "neg", cwd=d).stdout
        rep.check("neg: prod completes", "DONE_FG12" in p)
        pj = floats("JNU", p)
        pf = floats("FEDD", p)
        rep.check("neg: prod JNU stays positive despite I_CORE<0",
                  all(x > 0.0 for x in pj))
        rep.check("neg: prod FEDD in (0,1] after I_CORE floor",
                  all(0.0 < x <= 1.0 + 1e-8 for x in pf))
        rep.check("neg: ref fails - unfloored I_CORE seed propagates",
                  "DONE_FG12" not in r)
        rep.check("neg: ref diagnostic shows negative sweep intensity",
                  "Error for I_P" in r or "Error for I_M" in r)

        for case in ("jump", "wall", "thin"):
            p = run(prod, case, cwd=d).stdout
            r = run(ref, case, cwd=d).stdout
            rep.check(f"{case}: both complete",
                      "DONE_FG12" in p and "DONE_FG12" in r)
            rep.check(f"{case}: JNU identical (repairs not spurious)",
                      floats("JNU", p) == floats("JNU", r))
            rep.check(f"{case}: FEDD identical",
                      floats("FEDD", p) == floats("FEDD", r))

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
