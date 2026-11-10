#!/usr/bin/env python3
"""STEQ_ADVEC_V4 non-positive-V rejection
(new_main/subs/steq_advec_v4.f).

Alteration under test: INC_ADV=T combined with DO_NONMONO_SOB or
DO_NONMONO_KNOP and a non-positive velocity now STOPs with a
labeled error.  The advection stencil is a backward difference
using the inner neighbor as upstream, which is only valid for
radial outflow -- a sign-changing field would need the opposite
stencil, so guessing would be wrong physics.

Cases:
  nonpos  : SOB + INC_ADV, V(3)<0  -> STOP with labeled error
  zerov   : KNOP + INC_ADV, V(2)=0 -> STOP with labeled error
  flagoff : flags off, V(3)<0      -> guard skipped, returns clean
  good    : SOB + INC_ADV, V>0     -> passes guard, returns clean
  noinc   : SOB, INC_ADV=F, V<0    -> early return before guard
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (BINDIR, F90DIR, PROD_MOD_INC, TREE, Report,
                    compile_objs, link_prod, prod_obj, run)

GUARD = "the non-monotonic solver"


def build():
    md = BINDIR / "obj_steq_adv"
    objs = [prod_obj(TREE / "new_main" / "subs" / "steq_advec_v4.f")]
    objs += compile_objs(
        "steq_adv_drv", [F90DIR / "drv_steq_adv.f90"],
        moddir=md, extra_flags=PROD_MOD_INC)
    libs = [TREE / "lib" / a for a in
            ("libsubs.a", "libmod.a", "libsubs.a")]
    return link_prod("steq_adv", objs + libs, moddir=md)


def main():
    rep = Report("steq_adv")
    exe = build()

    for case in ("nonpos", "zerov"):
        r = run(exe, case, check=False)
        rep.check(
            f"{case} STOPs with labeled non-monotonic-advection "
            "error",
            GUARD in r.stdout and "DONE_" not in r.stdout,
            r.stdout.strip().splitlines()[-1]
            if r.stdout.strip() else "no output")

    for case in ("flagoff", "good", "noinc"):
        r = run(exe, case, check=False)
        rep.check(
            f"{case} passes the guard and returns",
            "DONE_" in r.stdout and GUARD not in r.stdout,
            r.stdout.strip().splitlines()[-1]
            if r.stdout.strip() else "no output")

    rep.summary()


if __name__ == "__main__":
    main()
