#!/usr/bin/env python3
"""INCIDENT_INTENSITY guards in GET_IBOUND (plane/get_ibound.f).

The altered code validates the incident grid at load (NCF >= 2,
finite positive strictly-decreasing NU_GRID), bounds the saved
frequency walk (CUR_NU_INDX / ML clamps), and returns IBOUND=0 for
non-finite FREQ.  The routine loads the file once on its first call,
so each case runs the driver in a scratch dir holding that case's
file.

The valid-file arithmetic is verified against hand-computed values:
with NU=[30,20,10], INT=[100,200,300], dilution 0.5, ANG=[1,0],
AVAR=[2,4], MU=[0.5,0.75] -> DIST=[3.0,2.5]; FREQ=25 gives
T2=75 -> IBOUND=[225,187.5]; FREQ=15 gives T2=125 -> [375,312.5].
"""
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (BINDIR, F90DIR, PROD_MOD_INC, TREE, WORKDIR,
                    Report, compile_objs, link_prod, prod_obj, run)

HEADER = """! manufactured incident file for unit tests
!Format date 2026
      3   !Number of continuum frequencies
      2   !Number of angles
    0.5   !Dilution factor used for scaling
!Angle grid
  1.0 0.0
!Intensity variation with angle
  2.0 4.0
"""

GRID_GOOD = "  30.0 100.0\n  20.0 200.0\n  10.0 300.0\n"
GRID_NONPOS = "  30.0 100.0\n  -5.0 200.0\n  10.0 300.0\n"
GRID_NOTDEC = "  30.0 100.0\n  40.0 200.0\n  10.0 300.0\n"

GOOD_FILE = HEADER + GRID_GOOD


def ncf1_file():
    return GOOD_FILE.replace(
        "      3   !Number of continuum frequencies",
        "      1   !Number of continuum frequencies",
    ).replace(GRID_GOOD, "  30.0 100.0\n")


def build():
    md = BINDIR / "obj_gib"
    objs = [prod_obj(TREE / "plane" / "get_ibound.f")]
    objs += compile_objs(
        "gib_drv", [F90DIR / "drv_get_ibound.f90"],
        moddir=md, extra_flags=PROD_MOD_INC)
    libs = [TREE / "lib" / a for a in ("libsubs.a", "libmod.a")]
    return link_prod("gib", objs + libs, moddir=md)


def casedir(tag, text):
    d = WORKDIR / f"gib_{tag}"
    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True)
    if text is not None:
        (d / "INCIDENT_INTENSITY").write_text(text)
    return d


def ib_values(stdout):
    vals = []
    for ln in stdout.splitlines():
        t = ln.split()
        if len(t) == 3 and t[0] == "IB":
            vals.append((float(t[1]), float(t[2])))
    return vals


def main():
    rep = Report("GET_IBOUND incident-grid guards")
    binary = build()

    d = casedir("good", GOOD_FILE)
    r = run(binary, "good", cwd=d, check=False)
    vals = ib_values(r.stdout)
    rep.check("good: driver completed",
              "DONE" in r.stdout and len(vals) == 5)
    if len(vals) == 5:
        rep.check("good: FREQ above grid -> 0", vals[0] == (0.0, 0.0))
        rep.check("good: FREQ=25 interpolation",
                  abs(vals[1][0] - 225.0) < 1e-9
                  and abs(vals[1][1] - 187.5) < 1e-9,
                  detail=str(vals[1]))
        rep.check("good: FREQ=15 interpolation",
                  abs(vals[2][0] - 375.0) < 1e-9
                  and abs(vals[2][1] - 312.5) < 1e-9,
                  detail=str(vals[2]))
        rep.check("good: FREQ below grid -> 0",
                  vals[3] == (0.0, 0.0))
        rep.check("good: NaN FREQ -> 0", vals[4] == (0.0, 0.0),
                  detail=str(vals[4]))

    cases = [
        ("ncf1", ncf1_file(), "need at least 2"),
        ("nonpos", HEADER + GRID_NONPOS, "non-finite or non-positive"),
        ("notdec", HEADER + GRID_NOTDEC, "not strictly decreasing"),
        ("badmu", GOOD_FILE, "Invalid incident MU range"),
        ("reorder", GOOD_FILE, "Frequencies out of order"),
    ]
    for tag, text, msg in cases:
        d = casedir(tag, text)
        r = run(binary, tag, cwd=d, check=False)
        blob = r.stdout + r.stderr
        rep.check(f"{tag}: {msg.strip()}", msg in blob,
                  detail=blob.strip().splitlines()[-1]
                  if blob.strip() else "")
        rep.check(f"{tag}: no clean finish", "DONE" not in blob)

    rep.summary()
    sys.exit(0)


if __name__ == "__main__":
    main()
