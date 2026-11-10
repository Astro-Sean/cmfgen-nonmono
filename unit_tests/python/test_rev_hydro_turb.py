#!/usr/bin/env python3
"""REV_HYDRO_TURB INT/TST guards (misc/rev_hydro_turb.f).

Alteration under test: the INT diagnostic walks outward to the
unique V=30 km/s crossing - undefined on a velocity field that is
not positive and strictly decreasing, and undefined when no
crossing exists.  TST takes LOG of a flux-like product that can
go non-positive on shocked profiles.  The fork validates V before
the walk, bounds the loop, and checks the log argument.

The shipped exe is run against a fabricated HYDRO file
(header + ND rows of 'R V E VdVdR dPdR g_TOT g_RAD g_ELEC Gamma'
plus the mass/gravity trailer) and MODEL_SPEC for ND.
"""
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import TREE, Report

EXE = TREE / "exe" / "rev_hydro_turb.exe"
R = [20.0, 15.0, 10.0, 5.0]


def write_hydro(d, v, vdvdr=5.0):
    rows = ["hydro header"]
    for i in range(len(R)):
        rows.append(
            "%.6f %.6f 1.0 %.6f 100.0 90.0 10.0 1.0 1.4"
            % (R[i], v[i], vdvdr))
    rows += [" ",
             "Surface gravity is: 1000.0",
             "Photospheric radius is: 5.0",
             "Stars mass is: 10.0"]
    Path(d, "HYDRO").write_text("\n".join(rows) + "\n")
    Path(d, "MODEL_SPEC").write_text("%d  [ND]\n" % len(R))


def run_case(v, xopt, vdvdr=5.0, tail="EX"):
    with tempfile.TemporaryDirectory() as d:
        write_hydro(d, v, vdvdr)
        # prompts: hydro file, output file, mass, vturb,
        # low/high depth, plane-parallel, then XOPT loop
        stdin = ("HYDRO\nREV_HYDRO\n10\n0\n1\n4\nF\n%s\n%s\n"
                 % (xopt, tail))
        r = subprocess.run([str(EXE)], cwd=d, input=stdin,
                           capture_output=True, text=True, timeout=60)
    return r.stdout + r.stderr


def main():
    rep = Report("REV_HYDRO_TURB INT/TST guards")
    rep.check("production exe present", EXE.exists())

    out = run_case([800, 600, 600, 200], "INT")
    rep.check("flat V: INT rejected with labelled error",
              "INT requires" in out and "strictly decreasing" in out)

    out = run_case([800, 600, -1.0, 200], "INT")
    rep.check("negative V: INT rejected", "INT requires" in out)

    out = run_case([20, 15, 10, 5], "INT")
    rep.check("no V=30 crossing: bounded-walk error",
              "could not locate" in out)

    out = run_case([800, 600, 40, 10], "INT")
    rep.check("crossing present: INT runs the table",
              "could not locate" not in out
              and "INT requires" not in out)

    out = run_case([800, 600, 400, 200], "TST", vdvdr=-5.0)
    rep.check("nonpositive log arg: TST rejected",
              "TST log argument" in out)

    out = run_case([800, 600, 400, 200], "TST")
    rep.check("positive log arg: TST completes",
              "TST log argument" not in out)

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
