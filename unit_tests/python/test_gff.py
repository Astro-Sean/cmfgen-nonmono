#!/usr/bin/env python3
"""GFF_VEC bad-T diagnostic (subs/gff.f).

The alteration flags a non-positive (or NaN) temperature to stdout
before LOG10 turns it into a NaN, so a corrupt temperature can be
traced to its depth index.

  ok   : all-positive T -> no diagnostic, finite positive Gaunts
  bad  : T(3) <= 0 -> exactly one 'bad T' line naming I=3
  nan  : T(2) = NaN -> .NOT.(T.GT.0) catches NaN too
  gff  : scalar GFF agrees with the vector path on valid T
"""
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (BINDIR, F90DIR, PROD_MOD_INC, TREE, WORKDIR,
                    Report, compile_objs, link_prod, prod_obj, run)


def build():
    md = BINDIR / "obj_gff"
    objs = [prod_obj(TREE / "subs" / "gff.f")]
    objs += compile_objs(
        "gff_drv", [F90DIR / "drv_gff.f90"],
        moddir=md, extra_flags=PROD_MOD_INC)
    return link_prod("gff", objs, moddir=md)


def gvs(out):
    return {
        int(m.group(1)): float(m.group(2))
        for m in re.finditer(r"GV\s+(\d+)\s+(\S+)", out)
    }


def main():
    rep = Report("GFF_VEC bad-T diagnostic")
    binary = build()
    d = WORKDIR / "gff"
    d.mkdir(parents=True, exist_ok=True)

    r = run(binary, "ok", cwd=d)
    vals = gvs(r.stdout)
    rep.check("ok: no diagnostic", "bad T" not in r.stdout)
    rep.check("ok: all Gaunts finite+positive",
              len(vals) == 5 and all(0 < v == v < 20 for v in vals.values()))

    r = run(binary, "bad", cwd=d, check=False)
    hits = re.findall(r"bad T: I=\s*(\d+)", r.stdout)
    rep.check("bad: one diagnostic", hits == ["3"], detail=str(hits))
    rep.check("bad: diagnostic names RNU/ZION",
              "RNU=" in r.stdout and "ZION=" in r.stdout)

    r = run(binary, "nan", cwd=d, check=False)
    hits = re.findall(r"bad T: I=\s*(\d+)", r.stdout)
    rep.check("nan: NaN T caught", hits == ["2"], detail=str(hits))

    r = run(binary, "gff", cwd=d)
    vals = gvs(r.stdout)
    sc = float(re.search(r"SC\s+(\S+)", r.stdout).group(1))
    rep.check("gff: scalar vs vector agree",
              abs(sc - vals[2]) / abs(vals[2]) < 1e-12,
              detail=f"scalar={sc} vec={vals[2]}")

    rep.summary()
    sys.exit(0)


if __name__ == "__main__":
    main()
