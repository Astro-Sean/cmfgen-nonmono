#!/usr/bin/env python3
"""RD_SN_DATA isotope-header sentinel (new_main/subs/rd_sn_data.f).

The alteration: NISO defaults to -1 so an ABSENT 'Number of
isotopes:' header is distinguished from an explicit '0' (upstream
used 0 for both and rejected explicit-zero files at the gate).

  niso_absent  header line missing -> STOP 'NISO or model time
               undefined' (sentinel path)
  niso_zero    'Number of isotopes: 0' explicit -> passes the gate,
               prints 'Read in SN_HYDRO_DATA data headers', dies on
               the truncated body -- proves NISO=0 is accepted
  age_missing  'Time(days) since explosion' absent -> same STOP as
               niso_absent (other sentinel checks still work)

Only the header gate is exercised; full SN_HYDRO_DATA ingestion
needs isotope/decay model state and stays integration-level.
"""
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (BINDIR, F90DIR, PROD_MOD_INC, TREE, WORKDIR,
                    Report, compile_objs, link_prod, prod_obj, run)

UNDEF = "NISO or model time undefined"
BANNER = "Read in SN_HYDRO_DATA data headers"


def sn_file(with_iso=True, iso_n=0, with_age=True):
    lines = [
        "! manufactured SN_HYDRO_DATA header for unit tests",
        "Number of data points: 4",
        "Number of mass fractions: 1",
    ]
    if with_iso:
        lines.append(f"Number of isotopes: {iso_n}")
    if with_age:
        lines.append("Time(days) since explosion: 100.0")
    lines.append("Radius grid")
    lines.append("  4.0 3.0 2.0 1.0")
    return "\n".join(lines) + "\n"


def build():
    md = BINDIR / "obj_rsn"
    objs = [prod_obj(TREE / "new_main" / "subs" / "rd_sn_data.f")]
    objs += compile_objs(
        "rsn_drv", [F90DIR / "drv_rd_sn_data.f90"],
        moddir=md, extra_flags=PROD_MOD_INC)
    libs = [TREE / "lib" / a for a in ("libsubs.a", "libmod.a")]
    return link_prod("rsn", objs + libs, moddir=md)


def casedir(tag, text):
    d = WORKDIR / f"rsn_{tag}"
    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True)
    (d / "SN_HYDRO_DATA").write_text(text)
    return d


def main():
    rep = Report("RD_SN_DATA isotope-header sentinel")
    binary = build()

    d = casedir("niso_absent", sn_file(with_iso=False))
    r = run(binary, "niso_absent", cwd=d, check=False)
    blob = r.stdout + r.stderr
    rep.check("niso_absent: gate STOPs", UNDEF in blob,
              detail=blob.strip().splitlines()[-1]
              if blob.strip() else "")
    rep.check("niso_absent: banner not reached", BANNER not in blob)

    d = casedir("niso_zero", sn_file(with_iso=True, iso_n=0))
    r = run(binary, "niso_zero", cwd=d, check=False)
    blob = r.stdout + r.stderr
    rep.check("niso_zero: gate passed", BANNER in blob,
              detail=blob.strip().splitlines()[-1]
              if blob.strip() else "")
    rep.check("niso_zero: no undefined rejection", UNDEF not in blob)

    d = casedir("age_missing", sn_file(with_age=False))
    r = run(binary, "age_missing", cwd=d, check=False)
    blob = r.stdout + r.stderr
    rep.check("age_missing: gate STOPs", UNDEF in blob)
    rep.check("age_missing: banner not reached", BANNER not in blob)

    rep.summary()
    sys.exit(0)


if __name__ == "__main__":
    main()
