#!/usr/bin/env python3
"""DO_V_PARAMS rejection under non-monotonic velocity flags
(new_main/subs/do_wind_vel_v1.f, altered 05-Oct-2026).

The guard stops only when the optional DO_V_PARAMS file is present,
so each case runs in its own scratch directory with the file written
(or withheld) by the harness:

  sob_file   DO_NONMONO_SOB=T + file present -> STOP with the
             documented incompatibility message
  knop_file  DO_NONMONO_KNOP=T + file present -> STOP, same message
  sob_nofile DO_NONMONO_SOB=T + no file -> clean return (default
             no-file path must be unchanged)
  off_file   flags off + file present (FST_IT > MAIN_COUNTER) ->
             clean early return, DONE_V_REVISION stays .FALSE.

The flag-off file parses through RD_OPTIONS_INTO_STORE before the
FIRST_IT early return, so the check covers both sides of the new
branch -- flags-off semantics are the upstream ones.
"""
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (BINDIR, F90DIR, PROD_MOD_INC, TREE, WORKDIR,
                    Report, compile_objs, link_prod, prod_obj, run)

STOP_MSG = "DO_V_PARAMS is incompatible with non-monotonic"

# Minimal option file; FST_IT=99 keeps MAIN_COUNTER=0 below the
# first-revision threshold so the flag-off path returns before any
# hydrodynamic work.
PARAMS = """      5       [N_GRAD_ITS]
      2       [DST]
     99       [FST_IT]
     10       [MAX_IT]
     20.0     [VSOUND]
"""


def build():
    md = BINDIR / "obj_dwv"
    objs = [prod_obj(TREE / "new_main" / "subs" / "do_wind_vel_v1.f")]
    objs += compile_objs(
        "dwv_drv", [F90DIR / "drv_do_wind_vel.f90"],
        moddir=md, extra_flags=PROD_MOD_INC)
    libs = [TREE / "lib" / a for a in ("libsubs.a", "libmod.a")]
    return link_prod("dwv", objs + libs, moddir=md)


def casedir(tag, with_file):
    d = WORKDIR / f"dwv_{tag}"
    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True)
    if with_file:
        (d / "DO_V_PARAMS").write_text(PARAMS)
    return d


def main():
    rep = Report("DO_WIND_VEL_V1 nonmono rejection")
    binary = build()

    for tag in ("sob_file", "knop_file"):
        d = casedir(tag, with_file=True)
        r = run(binary, tag, cwd=d, check=False)
        blob = r.stdout + r.stderr
        rep.check(f"{tag}: DO_V_PARAMS rejected", STOP_MSG in blob,
                  detail=blob.strip().splitlines()[-1]
                  if blob.strip() else "")
        rep.check(f"{tag}: no clean return",
                  "RETURNED" not in blob)

    d = casedir("sob_nofile", with_file=False)
    r = run(binary, "sob_nofile", cwd=d, check=False)
    blob = r.stdout + r.stderr
    rep.check("sob_nofile: returns without file",
              "RETURNED T" not in blob and "RETURNED" in blob
              and STOP_MSG not in blob,
              detail=blob.strip().splitlines()[-1]
              if blob.strip() else "")

    d = casedir("off_file", with_file=True)
    r = run(binary, "off_file", cwd=d, check=False)
    blob = r.stdout + r.stderr
    rep.check("off_file: flag-off path returns", "RETURNED" in blob,
              detail=blob.strip().splitlines()[-1]
              if blob.strip() else "")
    rep.check("off_file: no rejection text", STOP_MSG not in blob)
    rep.check("off_file: no V revision done",
              "RETURNED F" in blob)
    rep.check("off_file: no RVSIG_COL written",
              not (d / "RVSIG_COL").exists())

    rep.summary()
    sys.exit(0)


if __name__ == "__main__":
    main()
