#!/usr/bin/env python3
"""ADJUST_ATM_R_GRID velocity guard (misc/adjust_atm_r_grid.f).

Alteration under test: the connection search and log
interpolation use OLD_V as a coordinate, so it must be positive
and strictly decreasing.  Upstream silently picked an arbitrary
radius branch on a repeated velocity.

The driver links the production misc/adjust_atm_r_grid.o and
passes a synthetic R/V/TAU deck:

  flat : OLD_V repeats a value -> labelled STOP
  neg  : OLD_V goes negative -> labelled STOP
  good : decreasing positive OLD_V -> returns a full R grid
"""
import re
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (BINDIR, F90DIR, PROD_MOD_INC, TREE,
                    Report, compile_objs, link_prod, prod_obj, run)

GUARD = "OLD_V must be positive and strictly decreasing"


def build():
    md = BINDIR / "obj_adjgrid"
    objs = [prod_obj(TREE / "misc" / "adjust_atm_r_grid.f")]
    objs += compile_objs(
        "adjgrid", [F90DIR / "drv_adj_grid.f90"],
        moddir=md, extra_flags=PROD_MOD_INC)
    libs = [TREE / "lib" / a for a in ("libsubs.a", "libmod.a")]
    return link_prod("adjgrid", objs + libs, moddir=md)


def run_case(exe, case):
    # R_GRID_SELECTION and fort.* land in cwd; keep them in a
    # per-case scratch dir so the test tree stays clean.
    with tempfile.TemporaryDirectory() as d:
        r = run(exe, case, check=False, cwd=d)
    return r.stdout + r.stderr


def main():
    rep = Report("ADJUST_ATM_R_GRID velocity guard")
    exe = build()

    out = run_case(exe, "flat")
    rep.check("flat OLD_V: labelled error", GUARD in out)

    out = run_case(exe, "neg")
    rep.check("negative OLD_V: labelled error", GUARD in out)

    out = run_case(exe, "good")
    rep.check("good deck: no guard error", GUARD not in out)
    rep.check("good deck: completes", "DONE_" in out)
    vals = [float(m.group(1)) for m in re.finditer(
        r"^\s*\d+\s+(\S+)$", out, re.M)]
    rep.check("good deck: R grid returned and decreasing",
              len(vals) == 10
              and all(vals[i] > vals[i + 1] for i in range(9))
              and all(v > 0 for v in vals))

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
