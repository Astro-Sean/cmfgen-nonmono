#!/usr/bin/env python3
"""AVE_DC velocity guards (misc/ave_dc.f -> exe/ave_dc.exe).

Alteration under test: the population mapping uses V as the
independent coordinate with an intermediate Lagrangian-expanded
R_NEW grid, so both V1 and V2 must be positive and strictly
decreasing.  Upstream fed non-monotonic V into MON_INTERP and
died opaquely; the fork validates first with labelled STOPs.

The production executable (built Nov 9, post-alteration) is run
directly against fabricated fixtures:

  DC_FILE_INFO : two directory names then the DC file name
  <dir>/TEST_OUT : '!Format date' header, RSTAR LSTAR N ND,
                   per-depth rows 'R DION ED T JNK V' + N DC vals
  TEST__F_OSCDAT : copied from the real He2 oscillator file so
                   GENOSC_V9 succeeds before the guards run

  nonpos_v1 : V1 contains 0 -> 'V1 must be positive'
  flat_v1   : V1 repeats a value -> 'V1 and expanded R must be
              decreasing'
  nonpos_v2 : V2 negative -> 'V2 must be positive'
  flat_v2   : V2 rises -> 'V2 must be strictly decreasing'
  good      : identical V grids -> 'Done', output file with
              ED=sqrt(ED1*ED2) per row
"""
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import TREE, Report

EXE = TREE / "exe" / "ave_dc.exe"
OSC_SRC = (TREE.parent / "handoff_smoke" / "jfnk_probe_fz" /
           "He2_F_OSCDAT")

R = [20.0, 15.0, 10.0, 5.0]
V = [800.0, 600.0, 400.0, 200.0]
ED1 = [1.0e12, 8.0e11, 5.0e11, 3.0e11]
ED2 = [4.0e11, 3.2e11, 2.0e11, 1.2e11]
N_LEV = 30


def dc_file(v, ed, dc):
    rows = ["24-FEB-2004  !Format date",
            "20.0  1.0e10  %d  %d" % (N_LEV, len(R))]
    for i in range(len(R)):
        rows.append("%.6f 1.0e8 %.6e 3.0 0 %.6f" %
                    (R[i], ed[i], v[i]))
        rows.append(" ".join("%.6e" % dc for _ in range(N_LEV)))
    return "\n".join(rows) + "\n"


def casedir(d, v1, v2):
    Path(d, "dirone").mkdir()
    Path(d, "dirtwo").mkdir()
    Path(d, "DC_FILE_INFO").write_text("dirone/\ndirtwo/\nTEST_OUT\n")
    Path(d, "dirone", "TEST_OUT").write_text(dc_file(v1, ED1, 1.5))
    Path(d, "dirtwo", "TEST_OUT").write_text(dc_file(v2, ED2, 4.0))
    shutil.copy(OSC_SRC, Path(d, "TEST__F_OSCDAT"))


def run_case(name, v1, v2):
    with tempfile.TemporaryDirectory() as d:
        casedir(d, v1, v2)
        r = subprocess.run([str(EXE)], cwd=d, capture_output=True,
                           text=True, timeout=120)
        out = r.stdout + r.stderr
        out_file = Path(d, "TEST__IN")
        ftext = out_file.read_text() if out_file.exists() else ""
    return out, ftext


def main():
    rep = Report("AVE_DC velocity guards")
    rep.check("production exe present", EXE.exists())
    rep.check("real OSCDAT fixture present", OSC_SRC.exists())

    out, _ = run_case("nonpos_v1", [800, 600, 0.0, 200], V)
    rep.check("nonpos_v1: labelled V1 error",
              "V1 must be positive" in out)

    out, _ = run_case("flat_v1", [800, 600, 600.0, 200], V)
    rep.check("flat_v1: labelled V1/R_NEW error",
              "V1 and expanded R must be decreasing" in out)

    out, _ = run_case("nonpos_v2", V, [800, -1.0, 400, 200])
    rep.check("nonpos_v2: labelled V2 error",
              "V2 must be positive" in out)

    out, _ = run_case("flat_v2", V, [800, 700, 800.0, 200])
    rep.check("flat_v2: labelled V2 error",
              "V2 must be strictly decreasing" in out)

    out, ftext = run_case("good", V, V)
    rep.check("good: completes to 'Done'", "Done" in out)
    rep.check("good: output file written", len(ftext) > 0)
    eds = [float(m.group(3)) for m in re.finditer(
        r"^ *(\S+) +(\S+) +(\S+) +(\S+) +(\S+) +(\S+) +(\S+) +\d+",
        ftext, re.M)]
    exp = [(e1 * e2) ** 0.5 for e1, e2 in zip(ED1, ED2)]
    rep.check("good: ED = sqrt(ED1*ED2) at every depth",
              len(eds) == len(R) and all(
                  abs(a - e) < 1e-4 * e for a, e in zip(eds, exp)))

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
