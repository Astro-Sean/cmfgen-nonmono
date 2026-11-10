#!/usr/bin/env python3
"""COMP_GREY_V2/V4 grey-work sigma clamp
(new_main/mod_subs/comp_grey_v2.f and comp_grey_v4.f; the V4
twin carries the identical block tagged COMP_GREY_V4).

Alteration under test: under DO_NONMONO_SOB/DO_NONMONO_KNOP a
work copy SIGMA_GREY is clamped to >= -0.999 before reaching the
grey solvers -- their monotonic-expansion moment structures
sign-flip or divide near a pole when sigma < -1.  SIGMA and
SIGMA_TRUE stay untouched; stock runs are bit-identical (plain
copy, no message).  The clamp count prints at CONV_DIAG_LVL >= 1.

Cases:
  clamp : KNOP + CONV_DIAG_LVL=1, SIGMA=(-0.5,-1.5,-2.0,-0.9)
          -> "clamped 2 of 4" message, SIGMA still holds -1.5/-2.0
  quiet : same + CONV_DIAG_LVL=0 -> clamp runs, no message
  pos   : KNOP + all sigma >= -0.999 -> nothing to clamp, no msg
  off   : flags off + sub--1 sigmas -> no block, no message
"""
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (BINDIR, F90DIR, PROD_MOD_INC, TREE, Report,
                    compile_objs, link_prod, prod_obj, run)

MSG = "clamped"


def build():
    md = BINDIR / "obj_grey"
    objs = [
        prod_obj(TREE / "new_main" / "mod_subs" / "comp_grey_v2.f"),
        prod_obj(TREE / "new_main" / "mod_subs" / "comp_grey_v4.f"),
    ]
    objs += compile_objs(
        "grey_drv", [F90DIR / "drv_grey.f90"],
        moddir=md, extra_flags=PROD_MOD_INC)
    libs = [TREE / "lib" / a for a in
            ("libsubs.a", "libmod.a", "libsubs.a")]
    return link_prod("grey", objs + libs, moddir=md)


def sig_row(out):
    m = re.search(r"^SIG\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)", out, re.M)
    assert m, f"SIG row missing in:\n{out}"
    return [float(x) for x in m.groups()]


def main():
    rep = Report("grey")
    exe = build()

    out = run(exe, "clamp").stdout
    rep.check("clamp completes", "DONE_" in out)
    rep.check("clamp message reports 2 of 4",
              "clamped" in out and "2" in out.split("clamped")[1][:30])
    rep.check("clamped -0.999 bound in message",
              "-0.999" in out)
    s = sig_row(out)
    rep.check("SIGMA untouched by the work-copy clamp",
              s == [-0.5, -1.5, -2.0, -0.9], f"{s}")

    out = run(exe, "quiet").stdout
    rep.check("quiet completes", "DONE_" in out)
    rep.check("CONV_DIAG_LVL=0 suppresses the clamp message",
              MSG not in out)

    out = run(exe, "pos").stdout
    rep.check("pos completes", "DONE_" in out)
    rep.check("no message when no sigma needs clamping",
              MSG not in out)

    out = run(exe, "off").stdout
    rep.check("off completes", "DONE_" in out)
    rep.check("flag-off: block skipped entirely, no message",
              MSG not in out)
    s = sig_row(out)
    rep.check("off: SIGMA still holds sub--1 values (stock path)",
              s == [-0.5, -1.5, -2.0, -0.9], f"{s}")

    # V4 twin: identical clamp block, message tagged COMP_GREY_V4
    out = run(exe, "v4").stdout
    rep.check("v4 completes", "DONE_" in out)
    rep.check("v4 clamp message reports 2 of 4",
              "COMP_GREY_V4: clamped" in out
              and "2" in out.split("clamped")[1][:30])
    s = sig_row(out)
    rep.check("v4: SIGMA untouched",
              s == [-0.5, -1.5, -2.0, -0.9], f"{s}")

    rep.summary()


if __name__ == "__main__":
    main()
