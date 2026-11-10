#!/usr/bin/env python3
"""REV_RDINR IV-path velocity guards (misc/rev_rdinr.f).

Alteration under test: the IV option inverts R through VEL and
inserts user-supplied NEW_XV abscissae.  Both are only defined on
a positive, strictly decreasing velocity coordinate.  The fork
validates VEL before the inversion and NEW_XV after it, and
OUT_RGRID refuses to write RVSIG_HOPE when the interpolated V
goes non-positive.

stdin prompts: FILE_IN, FILE_OUT, OPTION, then for IV the
IST/IEND pair and NEW_XV values (0.0 terminates).
"""
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import TREE, Report

EXE = TREE / "exe" / "rev_rdinr.exe"
R = [20.0, 15.0, 10.0, 5.0]
V_MONO = [800.0, 600.0, 400.0, 200.0]
V_FLAT = [800.0, 600.0, 600.0, 200.0]
V_NEG = [800.0, 600.0, -1.0, 200.0]


def write_dc(d, v):
    rows = ["hdr1", "hdr2", "hdr3",
            "5.0  1.0e10  1  %d" % len(R), " "]
    for i in range(len(R)):
        rows.append("%.6f 1.0e8 1.0e12 3.0 1 %.6f 1.0" % (R[i], v[i]))
        rows.append("   1.0")
        rows.append(" ")
    Path(d, "dc_in.dat").write_text("\n".join(rows) + "\n")


def run_case(v, answers):
    with tempfile.TemporaryDirectory() as d:
        write_dc(d, v)
        stdin = "dc_in.dat\ndc_out.dat\nIV\n" + "\n".join(answers) + "\n"
        r = subprocess.run([str(EXE)], cwd=d, input=stdin,
                           capture_output=True, text=True, timeout=60)
        out_file = Path(d, "dc_out.dat")
        ftext = out_file.read_text() if out_file.exists() else ""
    return r.stdout + r.stderr, ftext


def main():
    rep = Report("REV_RDINR IV velocity guards")
    rep.check("production exe present", EXE.exists())

    out, _ = run_case(V_FLAT, ["1", "4"])
    rep.check("flat VEL: IV rejected before inversion",
              "IV requires" in out
              and "strictly decreasing velocity" in out)

    out, _ = run_case(V_NEG, ["1", "4"])
    rep.check("negative VEL: IV rejected", "IV requires" in out)

    out, _ = run_case(V_MONO, ["1", "4", "700", "900", "0"])
    rep.check("increasing NEW_XV: rejected after input",
              "IV velocity sequence" in out)

    out, _ = run_case(V_MONO, ["1", "4", "-5", "0"])
    rep.check("nonpositive NEW_XV: rejected",
              "IV velocity sequence" in out)

    out, ftext = run_case(V_MONO, ["1", "4", "700", "500", "300", "0"])
    rep.check("good IV run: passes both guards",
              "Error in REV_RDINR" not in out
              and "IV velocity sequence" not in out)
    rep.check("good IV run: reaches the R-grid remap",
              "MON_INTERP" in out or "RVSIG_HOPE" in out
              or len(ftext) > 3 * 4)

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
