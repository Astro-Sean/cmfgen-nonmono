#!/usr/bin/env python3
"""MERGE_SCR velocity guards (misc/merge_scr.f).

Alteration under test: identical guards to extrap_scr - each
model's V must be positive and strictly decreasing for the LOG(V)
population map, and the Lagrangian-expanded R grid must stay
ordered.  MERGE_SCR has no shipped exe target (it is absent from
the misc EXE list), so the production source is compiled with the
probed production flag set and linked exactly like the other
misc exes (OBJ1 + libs).

Fixtures are real checkpoints written by the mk_scr helper
through production SCR_RITE_V2, plus fabricated MODEL / VADAT
text files.
"""
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (BINDIR, F90DIR, PROD_MOD_INC, TREE,
                    Report, compile_objs, link_prod, prod_obj)

ND = 700
NT = 4


def build(name):
    md = BINDIR / ("obj_" + name)
    src = TREE / "misc" / "merge_scr.f"
    if name == "mkscr":
        objs = [prod_obj(TREE / "subs" / "scr_read_v2.f")]
        objs += compile_objs(
            "mkscr", [F90DIR / "mk_scr.f90"],
            moddir=md, extra_flags=PROD_MOD_INC)
        libs = [TREE / "lib" / a for a in ("libsubs.a", "libmod.a")]
        return link_prod("mkscr", objs + libs, moddir=md)
    # merge_scr: production source, no shipped object -> compile
    # with the probed production flags (same flags the misc
    # Makefile recipes use), then link like a misc exe.
    objs = compile_objs("merge_scr", [src], moddir=md, prod=True,
                        extra_flags=PROD_MOD_INC)
    obj1 = [TREE / "misc" / o for o in (
        "set_hydro_hd_loc.o", "update_pg_filename.o",
        "adjust_atm_r_grid.o", "rd_f_to_s_ids.o", "rd_ba_mat.o",
        "rd_f_to_s_ids_v2.o", "write_osc_v2.o")]
    obj1 = [o for o in obj1 if o.exists()]
    libs = [TREE / "lib" / a for a in (
        "libdisp.a", "libpgplt.a", "libsubs.a", "libspec.a",
        "libsubs.a", "libmod.a")]
    return link_prod("merge_scr", objs + obj1 + libs, moddir=md)


def mk_model_dir(parent, name, mkscr, vcase, sn_age):
    d = Path(parent, name)
    d.mkdir()
    (d / "MODEL").write_text(
        "%d !Number of depth points\n%d !Total number of variables\n"
        % (ND, NT))
    (d / "VADAT").write_text("%.1f [SN_AGE]\n" % sn_age)
    subprocess.run([str(mkscr), vcase], cwd=d, capture_output=True,
                   text=True, timeout=60)
    return (d / "SCRTEMP").exists() and (d / "POINT1").exists()


def run_case(mkscr, exe, c1, c2, age1, age2, age3):
    with tempfile.TemporaryDirectory() as d:
        mk_model_dir(d, "m1", mkscr, c1, age1)
        mk_model_dir(d, "m2", mkscr, c2, age2)
        stdin = "%s/m1/\n%s/m2/\n%s\n" % (d, d, age3)
        r = subprocess.run([str(exe)], cwd=d, input=stdin,
                           capture_output=True, text=True, timeout=120)
        new_scr = Path(d, "NEW_SCRTEMP")
        return r.stdout + r.stderr, new_scr.exists()


def main():
    rep = Report("MERGE_SCR velocity guards")
    mkscr = build("mkscr")
    exe = build("merge_scr")
    rep.check("merge_scr built from production source",
              Path(exe).exists())

    out, _ = run_case(mkscr, exe, "flat", "mono", 12.5, 25.0, 50.0)
    rep.check("flat V in model 1: rejected",
              "V must be positive and strictly decreasing" in out)

    out, _ = run_case(mkscr, exe, "mono", "flat", 12.5, 25.0, 50.0)
    rep.check("flat V in model 2: rejected",
              "V must be positive and strictly decreasing" in out)

    out, _ = run_case(mkscr, exe, "mono", "neg", 12.5, 25.0, 50.0)
    rep.check("nonpositive V tail: rejected",
              "V must be positive" in out)

    out, _ = run_case(mkscr, exe, "mono", "mono", 12.5, 25.0, 5.0)
    rep.check("backwards age: expansion reorder rejected",
              "reordered the R grid" in out)

    out, new_scr = run_case(mkscr, exe, "mono", "mono", 12.5, 25.0,
                            50.0)
    rep.check("good decks: completes", "You should now do" in out
              or "Set new R" in out)
    rep.check("good decks: NEW_SCRTEMP written", new_scr)

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
