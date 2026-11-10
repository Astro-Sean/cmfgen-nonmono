#!/usr/bin/env python3
"""REV_RVSIG velocity-coordinate gate (misc/rev_rvsig.f).

Alteration under test: revisions that locate an interval by V
(TAU, NEWG, SCLR, MDOT, STOP) are undefined on a velocity field
that is not positive and strictly decreasing - a repeated V picks
an arbitrary radius branch.  The fork checks OLD_V up front and
rejects those options with a labelled STOP while leaving the
radius-coordinate edits (SCLV and friends) usable.

The shipped exe is driven through stdin with a fabricated
RVSIG_COL_OLD: ND row of 'R V SIGMA' under a '!Number of depth
points' header line.
"""
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import TREE, Report

EXE = TREE / "exe" / "rev_rvsig.exe"
R = [20.0, 15.0, 10.0, 5.0]
V_MONO = [800.0, 600.0, 400.0, 200.0]
V_FLAT = [800.0, 600.0, 600.0, 200.0]
V_NEG = [800.0, 600.0, -1.0, 200.0]
GUARD = "requires a positive"


def write_rvsig(d, v):
    rows = ["%d !Number of depth points" % len(R)]
    rows += ["%.6f %.6f %.6f" % (R[i], v[i], -1.0)
             for i in range(len(R))]
    Path(d, "RVSIG_COL_OLD").write_text("\n".join(rows) + "\n")


def run_case(v, option, extra=""):
    with tempfile.TemporaryDirectory() as d:
        write_rvsig(d, v)
        # GEN_IN prompts: file name (blank = default), RD_DENSITY,
        # OPTION.  Extra answers feed whichever prompt comes next.
        stdin = "\nF\n%s\n%s" % (option, extra)
        r = subprocess.run([str(EXE)], cwd=d, input=stdin,
                           capture_output=True, text=True, timeout=60)
    return r.stdout + r.stderr


def main():
    rep = Report("REV_RVSIG velocity-coordinate gate")
    rep.check("production exe present", EXE.exists())

    for opt in ("TAU", "NEWG", "SCLR", "MDOT", "STOP"):
        out = run_case(V_FLAT, opt)
        rep.check("flat V: %s rejected" % opt, GUARD in out)

    out = run_case(V_NEG, "TAU")
    rep.check("negative V: TAU rejected", GUARD in out)

    out = run_case(V_FLAT, "SCLV")
    rep.check("flat V: radius-coord SCLV passes the gate",
              GUARD not in out)

    out = run_case(V_MONO, "TAU")
    rep.check("monotonic V: TAU passes the gate", GUARD not in out)

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
