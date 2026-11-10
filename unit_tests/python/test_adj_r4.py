#!/usr/bin/env python3
"""ADJUST_R_GRID_V4 nonmono gates (new_main/mod_subs/adjust_r_grid_v4.f).

Alterations under test:
  - GRID_METH=VEL_SPACE is rejected under DO_NONMONO_SOB /
    DO_NONMONO_KNOP: DO_VEL_REGRID_V2 interpolates R(V), which
    needs a single-valued monotonic velocity field.
  - the interpolated SIGMA < -1 clamp/STOP is skipped under the
    same flags: sigma < -1 is expected across a reversal and is
    handled by RECOMPUTE_SIGMA_TRUE (SOB) or the KNOP solver.

The driver fills MOD_CMFGEN arrays and writes an
ADJUST_R_DEFAULTS deck per case.  OB_ONLY nudges R(2) so the
sigma block is always reached; a V(ND) spike forces the
interpolated sigma below -1.
"""
import re
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (BINDIR, F90DIR, PROD_MOD_INC, TREE,
                    Report, compile_objs, link_prod, prod_obj, run)


def build():
    md = BINDIR / "obj_adjr4"
    objs = [prod_obj(TREE / "new_main" / "mod_subs" /
                     "adjust_r_grid_v4.f")]
    objs += compile_objs(
        "adjr4", [F90DIR / "drv_adj_r4.f90"],
        moddir=md, extra_flags=PROD_MOD_INC)
    libs = [TREE / "lib" / a for a in ("libsubs.a", "libmod.a")]
    return link_prod("adjr4", objs + libs, moddir=md)


def run_case(exe, case):
    with tempfile.TemporaryDirectory() as d:
        r = run(exe, case, check=False, cwd=d)
    return r.stdout + r.stderr


def main():
    rep = Report("ADJUST_R_GRID_V4 nonmono gates")
    exe = build()

    out = run_case(exe, "velsob")
    rep.check("VEL_SPACE + SOB: rejected",
              "VEL_SPACE is incompatible" in out)

    out = run_case(exe, "velknop")
    rep.check("VEL_SPACE + KNOP: rejected",
              "VEL_SPACE is incompatible" in out)

    out = run_case(exe, "sigoff")
    rep.check("reversal V + flags off: upstream STOP kept",
              "SIGMA is not monotonic" in out)

    out = run_case(exe, "sigsob")
    rep.check("reversal V + SOB: completes",
              "SIGMA is not monotonic" not in out
              and "ADJ_ DONE_R_REV=T" in out)
    m = re.search(r"ADJ_ SIGMA_NDm1=\s*(\S+)", out)
    rep.check("reversal V + SOB: sigma left below -1",
              m is not None and float(m.group(1)) < -1.0)

    out = run_case(exe, "sigknop")
    rep.check("reversal V + KNOP: completes",
              "SIGMA is not monotonic" not in out
              and "ADJ_ DONE_R_REV=T" in out)

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
