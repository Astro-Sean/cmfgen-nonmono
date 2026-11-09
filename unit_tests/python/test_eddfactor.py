#!/usr/bin/env python3
"""OPEN_RW_EDDFACTOR alterations (new_main/subs/open_rw_eddfactor.f).

Two fixes under test, both file-driven:

  * Filename gate: upstream tested ``FILENAME .EQ. FILENAME`` (always
    true), so the FINISH_REC "not computed" recompute fired for
    every file name.  The fix gates on 'EDDFACTOR'/'ES_J_CONV'.
    Verified by creating a file pair under each name and re-opening:
    EDDFACTOR warns + recomputes (ACCESS_F=15, COMPUTE=T); NOT_EDD
    returns silently (ACCESS_F=0, COMPUTE=F).

  * RV_REC IOSTAT guard in OUT_RV_TO_EDDFACTOR: a file truncated
    below the pointer record reports IOSTAT, warns "treating file as
    empty", and writes a fresh pointer instead of aborting.

Plus the fresh-create path: files + _INFO are written and the
RV_REC pointer record holds INITIAL_RV_REC=10.
"""
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (BINDIR, F90DIR, PROD_MOD_INC, TREE, WORKDIR,
                    Report, compile_objs, link_prod, prod_obj, run)


def build():
    md = BINDIR / "obj_edd"
    objs = [prod_obj(TREE / "new_main" / "subs"
                     / "open_rw_eddfactor.f")]
    objs += compile_objs(
        "edd_drv", [F90DIR / "drv_eddfactor.f90"],
        moddir=md, extra_flags=PROD_MOD_INC)
    libs = [TREE / "lib" / a for a in ("libsubs.a", "libmod.a")]
    return link_prod("edd", objs + libs, moddir=md)


def casedir(tag):
    d = WORKDIR / f"edd_{tag}"
    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True)
    return d


def marker(stdout, key):
    for ln in stdout.splitlines():
        t = ln.split()
        if len(t) >= 2 and t[0] == key:
            return t[1:]
    return None


def main():
    rep = Report("OPEN_RW_EDDFACTOR alterations")
    binary = build()

    d = casedir("fresh")
    r = run(binary, "fresh", cwd=d, check=False)
    rep.check("fresh: completed", "DONE" in r.stdout)
    af = marker(r.stdout, "AF")
    rep.check("fresh: ACCESS_F=15, COMPUTE=T",
              af is not None and af[0] == "15" and af[1] == "T",
              detail=str(af))
    rep.check("fresh: EDDFACTOR written",
              marker(r.stdout, "FEX") == ["T"])
    rep.check("fresh: _INFO written",
              marker(r.stdout, "IEX") == ["T"])
    rep.check("fresh: RV_REC pointer = 10",
              marker(r.stdout, "RVP") == ["10"])

    d = casedir("gate")
    r = run(binary, "gate", cwd=d, check=False)
    blob = r.stdout
    rep.check("gate: completed", "DONE" in blob)
    warn_n = blob.count("All values not")
    # EDDFACTOR warns twice (create leaves FINISH_REC=0 -> second call
    # warns); NOT_EDD must produce no warning under the fixed gate.
    rep.check("gate: warnings only for EDDFACTOR name", warn_n == 1,
              detail=f"{warn_n} warnings")
    ed = marker(blob, "EDDAF")
    rep.check("gate: EDDFACTOR recomputes (AF=15,T)",
              ed is not None and ed[0] == "15" and ed[1] == "T",
              detail=str(ed))
    nd = marker(blob, "NOTED")
    rep.check("gate: NOT_EDD passes silently (AF=0,F)",
              nd is not None and nd[0] == "0" and nd[1] == "F",
              detail=str(nd))

    d = casedir("trunc")
    r = run(binary, "trunc", cwd=d, check=False)
    blob = r.stdout
    rep.check("trunc: completed", "DONE" in blob)
    rep.check("trunc: empty-file warning",
              "treating file as empty" in blob)
    af = marker(blob, "AF")
    rep.check("trunc: ACCESS_F bumped by 3",
              af is not None and af[0] == "18", detail=str(af))
    rep.check("trunc: fresh RV pointer written",
              marker(blob, "RVP") == ["15"])

    rep.summary()
    sys.exit(0)


if __name__ == "__main__":
    main()
