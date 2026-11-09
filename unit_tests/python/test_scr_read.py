#!/usr/bin/env python3
"""SCR_RITE_V2/SCR_READ_V2 checkpoint system (subs/scr_read_v2.f).

Alterations under test:
  - WRITFAIL output flag on IOSTAT write errors (torn
    POINT/SCRTEMP pairs used to resume silently)
  - DUMP_FG_INEXT_GS / DUMP_MOM_JNEXT_GS checkpoint hooks for the
    Knop GS state (no-op here, arrays unallocated)
  - INQUIRE/CLOSE unit hygiene on failed opens

ND=700, NT=4 forces multi-record R/V/SIGMA and POPS writes
(N_PER_REC=2047 on this platform).

  rt     : write + read-last round trip, bitwise on all arrays,
           flags and metadata returned
  rt2    : two iterations written; last-read returns write 2,
           explicit IREC_RD=1 returns write 1
  newmod : no files -> NEWMOD=T, IREC_RD=0
  wrfail : unwritable dir -> WRITFAIL=T instead of torn pair
  norv   : WRITE_RVSIG=F -> RVSIG_WRITTEN=F; header R/V recs
           still round-trip (they are written unconditionally),
           only per-iteration R/V copies are skipped
"""
import os
import re
import shutil
import stat
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (BINDIR, F90DIR, PROD_MOD_INC, TREE, WORKDIR,
                    Report, compile_objs, link_prod, prod_obj, run)


def build():
    md = BINDIR / "obj_scr"
    objs = [prod_obj(TREE / "subs" / "scr_read_v2.f")]
    objs += compile_objs(
        "scr_drv", [F90DIR / "drv_scr.f90"],
        moddir=md, extra_flags=PROD_MOD_INC)
    libs = [TREE / "lib" / a for a in ("libsubs.a", "libmod.a")]
    return link_prod("scr", objs + libs, moddir=md)


def casedir(tag):
    d = WORKDIR / f"scr_{tag}"
    if d.exists():
        d.chmod(0o755)
        shutil.rmtree(d)
    d.mkdir(parents=True)
    return d


def main():
    rep = Report("SCR_RITE/READ_V2 checkpoint alterations")
    binary = build()

    d = casedir("rt")
    r = run(binary, "rt", cwd=d)
    m = re.search(r"RD\s+(\w)\s+(\d+)\s+(\w)", r.stdout)
    rep.check("rt: write ok (WRITFAIL=F, IREC=1)",
              re.search(r"WF1\s+F\s+1", r.stdout) is not None)
    rep.check("rt: NEWMOD=F, IREC_RD=1, RVSIG=T",
              m and m.groups() == ("F", "1", "T"),
              detail=r.stdout.strip().splitlines()[-3] if m else "?")
    rep.check("rt: NITSF=17, LST_NG=8 returned",
              re.search(r"META\s+17\s+-1\s+8", r.stdout) is not None)
    rep.check("rt: bitwise round-trip all arrays",
              re.search(r"EQ1\s+T\s*$", r.stdout, re.M) is not None)
    rep.check("rt: SCRTEMP + POINT files written",
              (d / "SCRTEMP").exists() and (d / "POINT1").exists())

    d = casedir("rt2")
    r = run(binary, "rt2", cwd=d)
    rep.check("rt2: second write IREC=2",
              re.search(r"WF2\s+F\s+2", r.stdout) is not None)
    rep.check("rt2: last read -> write 2",
              re.search(r"EQ1\s+T\s*$", r.stdout, re.M) is not None)
    rep.check("rt2: IREC_RD=1 -> write 1",
              re.search(r"EQ2\s+T\s+1", r.stdout) is not None)

    d = casedir("newmod")
    r = run(binary, "newmod", cwd=d, check=False)
    rep.check("newmod: NEWMOD=T, IREC_RD=0",
              re.search(r"NM\s+T\s+0", r.stdout) is not None,
              detail=r.stdout.strip().splitlines()[-1]
              if r.stdout.strip() else "")

    d = casedir("wrfail")
    d.chmod(0o555)
    try:
        r = run(binary, "wrfail", cwd=d, check=False)
    finally:
        d.chmod(0o755)
    rep.check("wrfail: WRITFAIL=T on unwritable dir",
              re.search(r"WF\s+T", r.stdout) is not None,
              detail=(r.stdout + r.stderr).strip().splitlines()[-1]
              if (r.stdout + r.stderr).strip() else "")

    d = casedir("norv")
    r = run(binary, "norv", cwd=d)
    m = re.search(r"RD\s+(\w)\s+(\d+)\s+(\w)", r.stdout)
    rep.check("norv: RVSIG_WRITTEN=F", m and m.group(3) == "F")
    rep.check("norv: header R/V + POPS round-trip",
              re.search(r"EQ1\s+T\s*$", r.stdout, re.M) is not None)

    rep.summary()
    sys.exit(0)


if __name__ == "__main__":
    main()
