#!/usr/bin/env python3
"""EXTRAP_SCR velocity guards (misc/extrap_scr.f -> exe/extrap_scr.exe).

Alteration under test: the unequal-grid branch maps populations on
LOG(V) and the R update is Lagrangian expansion, so each model's V
must be positive and strictly decreasing and the expanded R grid
must stay ordered.  Upstream fed bad coordinates into MON_INTERP
or wrote a reordered checkpoint; the fork validates first.

Fixtures are real checkpoints: the mk_scr helper writes SCRTEMP,
POINT1 and POINT2 through the production SCR_RITE_V2, and MODEL /
VADAT are fabricated text files.  The production extrap_scr.exe
is then driven through stdin:

  dir1 path, dir2 path, then the SN_AGE3 answer
  (default = age2*age2/age1; a younger age makes T1 negative so
   R+V*T1 reorders, reaching the expansion guard).
"""
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (BINDIR, F90DIR, PROD_MOD_INC, TREE,
                    Report, compile_objs, link_prod, prod_obj)

EXE = TREE / "exe" / "extrap_scr.exe"
ND = 700
NT = 4


def build():
    md = BINDIR / "obj_mkscr"
    objs = [prod_obj(TREE / "subs" / "scr_read_v2.f")]
    objs += compile_objs(
        "mkscr", [F90DIR / "mk_scr.f90"],
        moddir=md, extra_flags=PROD_MOD_INC)
    libs = [TREE / "lib" / a for a in ("libsubs.a", "libmod.a")]
    return link_prod("mkscr", objs + libs, moddir=md)


def mk_model_dir(parent, name, mkscr, vcase, sn_age):
    d = Path(parent, name)
    d.mkdir()
    (d / "MODEL").write_text(
        "%d !Number of depth points\n%d !Total number of variables\n"
        % (ND, NT))
    (d / "VADAT").write_text("%.1f [SN_AGE]\n" % sn_age)
    r = subprocess.run([str(mkscr), vcase], cwd=d, capture_output=True,
                       text=True, timeout=60)
    ok = (d / "SCRTEMP").exists() and (d / "POINT1").exists()
    return ok


def run_case(mkscr, c1, c2, age1, age2, age3):
    with tempfile.TemporaryDirectory() as d:
        mk_model_dir(d, "m1", mkscr, c1, age1)
        mk_model_dir(d, "m2", mkscr, c2, age2)
        stdin = "%s/m1/\n%s/m2/\n%s\n" % (d, d, age3)
        r = subprocess.run([str(EXE)], cwd=d, input=stdin,
                           capture_output=True, text=True, timeout=120)
        new_scr = Path(d, "NEW_SCRTEMP")
        return r.stdout + r.stderr, new_scr.exists()


def main():
    rep = Report("EXTRAP_SCR velocity guards")
    rep.check("production exe present", EXE.exists())
    mkscr = build()

    out, _ = run_case(mkscr, "flat", "mono", 12.5, 25.0, 50.0)
    rep.check("flat V in model 1: rejected",
              "V must be positive and strictly decreasing" in out)

    out, _ = run_case(mkscr, "mono", "flat", 12.5, 25.0, 50.0)
    rep.check("flat V in model 2: rejected",
              "V must be positive and strictly decreasing" in out)

    out, _ = run_case(mkscr, "mono", "neg", 12.5, 25.0, 50.0)
    rep.check("nonpositive V tail: rejected",
              "V must be positive" in out)

    out, _ = run_case(mkscr, "mono", "mono", 12.5, 25.0, 5.0)
    rep.check("backwards age: expansion reorder rejected",
              "reordered the R grid" in out)

    out, new_scr = run_case(mkscr, "mono", "mono", 12.5, 25.0, 50.0)
    rep.check("good decks: completes", "You should now do" in out)
    rep.check("good decks: NEW_SCRTEMP written", new_scr)

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
