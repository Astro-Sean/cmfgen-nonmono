#!/usr/bin/env python3
"""READ/WRITE_BCD_MAT error-path alterations (new_main/subs/read_bcd_mat.f).

Alterations under test:
  - unrecognized TYPE now STOPs in both routines (upstream printed
    the error and returned, leaving output arrays untouched)
  - WRITE_BCD_MAT's error now names the correct routine
  - valid round-trips must stay bitwise exact (unformatted I/O)

Cases:
  rt BC / rt C / rt D  write+read round-trip of every array TYPE
                       touches, bitwise compared
  rd missing           absent scratch file -> 'Unable to open' STOP
  rd bad type          valid file, TYPE='X' -> 'Unrecognized TYPE'
                       + 'Error in READ_BCD_MAT' + STOP
  wr bad type          TYPE='X' -> 'Unrecognized TYPE'
                       + 'Error in WRITE_BCD_MAT' + STOP
"""
import os
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (BINDIR, F90DIR, PROD_MOD_INC, TREE, WORKDIR,
                    Report, compile_objs, link_prod, prod_obj, run)

DEPTH_INDX = 7


def build():
    md = BINDIR / "obj_bcd"
    objs = [prod_obj(TREE / "new_main" / "subs" / "read_bcd_mat.f")]
    objs += compile_objs(
        "bcd_drv", [F90DIR / "drv_bcd_mat.f90"],
        moddir=md, extra_flags=PROD_MOD_INC)
    libs = [TREE / "lib" / a for a in ("libsubs.a", "libmod.a")]
    return link_prod("bcd", objs + libs, moddir=md)


def casedir(tag):
    d = WORKDIR / f"bcd_{tag}"
    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True)
    return d


def run_case(binary, tag, case, typ, d):
    os.environ["SCRATCH_DIR"] = str(d) + "/"
    return run(binary, case, typ, cwd=d, check=False)


def main():
    rep = Report("READ/WRITE_BCD_MAT alterations")
    binary = build()

    for typ in ("BC", "C", "D"):
        d = casedir(f"rt_{typ.lower()}")
        r = run_case(binary, f"rt_{typ}", "rt", typ, d)
        rep.check(f"rt {typ}: bitwise round-trip",
                  "RT_OK" in r.stdout,
                  detail=r.stdout.strip().splitlines()[-1]
                  if r.stdout.strip() else "")

    d = casedir("rd_missing")
    r = run_case(binary, "rd_missing", "rd", "BC", d)
    blob = r.stdout + r.stderr
    rep.check("rd missing: open error STOP",
              "Unable to open" in blob and "Error in READ_BCD_MAT" in blob)
    rep.check("rd missing: did not return READ",
              "\nREAD\n" not in "\n" + r.stdout)

    d = casedir("rd_badtype")
    (d / f"XSCRATCH{DEPTH_INDX:04d}").write_bytes(b"")
    r = run_case(binary, "rd_badtype", "rd", "X", d)
    blob = r.stdout + r.stderr
    rep.check("rd bad TYPE: STOPs", "Unrecognized TYPE" in blob)
    rep.check("rd bad TYPE: names READ_BCD_MAT",
              "Error in READ_BCD_MAT" in blob)
    rep.check("rd bad TYPE: did not return READ",
              "\nREAD\n" not in "\n" + r.stdout)

    d = casedir("wr_badtype")
    r = run_case(binary, "wr_badtype", "wr", "X", d)
    blob = r.stdout + r.stderr
    rep.check("wr bad TYPE: STOPs", "Unrecognized TYPE" in blob)
    rep.check("wr bad TYPE: names WRITE_BCD_MAT (corrected msg)",
              "Error in WRITE_BCD_MAT" in blob)
    rep.check("wr bad TYPE: no wrong routine name",
              "Error in READ_BCD_MAT" not in blob)

    rep.summary()
    sys.exit(0)


if __name__ == "__main__":
    main()
