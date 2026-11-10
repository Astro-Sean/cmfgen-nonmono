#!/usr/bin/env python3
"""FIDDLE_POP_CORRECTIONS_V2 nonmono scale-cap + NaN floor
(new_main/subs/fiddle_pop_corrections_v2.f).

Alterations under test:
1. Nonmono scale-scan cap: each row's contribution to the SCALE
   scan is capped at the depth's apply-clamp bound, so a huge
   demanded correction (1E70 here, the KNOP cold-start mode) no
   longer freezes every other correction at the depth.  Stock
   uses +-HUGE, so the big row still freezes the column.
2. DO_SOB_FG_KNOP floors non-positive/non-finite corrected POPS:
   species to 1e-100, the T row to T_MIN.

Arithmetic (CHANGE_LIM=5 -> BIG_LIM=0.8, LIT_LIM=-4,
POP_LIM=500 -> depth cap 0.998):
  depth1, stock   : SCALE=0.8/1E70 -> POPS(1)=0.2 (own clamp),
                    POPS(2..4) frozen at 1.0
  depth1, nonmono : SCALE=0.8/0.998=0.8016 -> POPS(1)=0.002,
                    POPS(2)=0.992, POPS(3)=1.008, T=2.9976
  depth2, fgknop  : all-NaN STEQ -> POPS(1:3)=1e-100, T=0.4
  depth2, stocknan: NaN propagates into POPS
"""
import math
import re
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (BINDIR, F90DIR, PROD_MOD_INC, TREE, WORKDIR,
                    Report, compile_objs, link_prod, prod_obj, run)


def build():
    md = BINDIR / "obj_fpc"
    objs = [
        prod_obj(TREE / "new_main" / "subs"
                 / "fiddle_pop_corrections_v2.f")]
    objs += compile_objs(
        "fpc_drv", [F90DIR / "drv_fpc.f90"],
        moddir=md, extra_flags=PROD_MOD_INC)
    libs = [TREE / "lib" / a for a in
            ("libsubs.a", "libmod.a", "libsubs.a")]
    return link_prod("fpc", objs + libs, moddir=md)


def casedir(tag):
    d = WORKDIR / f"fpc_{tag}"
    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True)
    return d


def _num(x):
    if x == "NaN":
        return math.nan
    m = re.match(r"^(-?\d+\.\d+)([+-]\d+)$", x)
    return float(f"{m.group(1)}E{m.group(2)}") if m else float(x)


def row(out, key):
    m = re.search(rf"^{key}\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)",
                  out, re.M)
    assert m, f"{key} row missing in:\n{out}"
    return [_num(x) for x in m.groups()]


def close(a, b, tol=2.0e-3):
    return abs(a - b) <= tol * max(1.0, abs(b))


def main():
    rep = Report("fpc")
    exe = build()

    # stock: 1E70 row freezes the column at SCALE~1E-70
    d = casedir("stock")
    out = run(exe, "stock", cwd=d).stdout
    rep.check("stock completes", "DONE_" in out)
    p1 = row(out, "P1")
    rep.check("stock P1(1)=0.2 (own clamp)",
              close(p1[0], 0.2), f"{p1[0]}")
    rep.check("stock P1(2)~1.0 (frozen by SCALE~1E-70)",
              close(p1[1], 1.0, 1e-6), f"{p1[1]}")
    rep.check("stock P1(3)~1.0 (frozen)",
              close(p1[2], 1.0, 1e-6), f"{p1[2]}")
    p2 = row(out, "P2")
    rep.check("stock P2 corrects normally",
              close(p2[0], 0.995), f"{p2[0]}")

    # nonmono: depth cap (=BIG_LIM=0.8) rescues the small
    # corrections; the big row still hits its own apply clamp
    d = casedir("nonmono")
    out = run(exe, "nonmono", cwd=d).stdout
    rep.check("nonmono completes", "DONE_" in out)
    p1 = row(out, "P1")
    rep.check("nonmono P1(1)=0.2 (own apply clamp, same bound)",
              close(p1[0], 0.2), f"{p1[0]}")
    rep.check("nonmono P1(2)=0.995 (full correction applied)",
              close(p1[1], 0.995), f"{p1[1]}")
    rep.check("nonmono P1(3)=1.005 (increase applied)",
              close(p1[2], 1.005), f"{p1[2]}")
    rep.check("nonmono T row corrected",
              close(p1[3], 2.997), f"{p1[3]}")

    # fgknop: NaN STEQ depth floored instead of poisoning POPS
    d = casedir("fgknop")
    out = run(exe, "fgknop", cwd=d).stdout
    rep.check("fgknop completes", "DONE_" in out)
    p2 = row(out, "P2")
    rep.check("fgknop species rows floored to 1e-100",
              all(close(v, 1.0e-100, 1.0) for v in p2[:3]),
              f"{p2[:3]}")
    rep.check("fgknop T row floored to T_MIN=0.4",
              close(p2[3], 0.4), f"{p2[3]}")
    p1 = row(out, "P1")
    rep.check("fgknop depth-1 corrections still applied",
              close(p1[1], 0.995), f"{p1[1]}")

    # stocknan: flag off, NaN propagates (documents old behavior)
    d = casedir("stocknan")
    out = run(exe, "stocknan", cwd=d).stdout
    rep.check("stocknan completes", "DONE_" in out)
    p2 = row(out, "P2")
    rep.check("stocknan POPS rows are NaN (floor not applied)",
              all(math.isnan(v) for v in p2), f"{p2}")

    rep.summary()


if __name__ == "__main__":
    main()
