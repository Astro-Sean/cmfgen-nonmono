#!/usr/bin/env python3
"""GET_JH_AT_PREV_TIME_STEP file/velocity guards
(plane/get_jh_at_prev_time_step.f).

Alterations under test:
  - IOSTAT on header/R-V/J-H record reads: truncated
    JH_AT_OLD_TIME fails loudly instead of feeding garbage to the
    interpolation.
  - OLD_V (file grid) and V (model grid) must be positive and
    strictly decreasing before LOG(V) is used as the MON_INTERP
    coordinate.

Manufactured file: ND_OLD=4, NCF_OLD=100, ST_IREC=4, NUST=100..1,
JST(:,ML)=ML so the interpolated NEW_J is exactly the frequency
weight.  Cases:

  good      NU=50 sits exactly on rec 51 -> NEW_J=51, DT=1e4
  ogrey     grey branch -> NEW_J=1.5 everywhere
  nu_high   NU above grid -> max-freq clamp + warning
  nu_low    NU below grid -> last-rec clamp + warning
  bad_hdr   REC 3 hole reads as zeros -> ST_IREC=0 fails at the
            R/V IOSTAT check
  bad_rv    file truncated at REC 3 -> 'R/V rec' STOP
  oldv_neg  OLD_V<=0 -> 'must be positive' STOP
  oldv_flat OLD_V repeated -> 'strictly decreasing' STOP
  v_neg / v_flat / inner_mm / outer_big / dt_zero
"""
import re
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (BINDIR, F90DIR, PROD_MOD_INC, TREE, WORKDIR,
                    Report, compile_objs, link_prod, prod_obj, run)

STOP_CASES = [
    ("bad_hdr", "Error reading JH_AT_OLD_TIME R/V rec"),
    ("bad_rv", "Error reading JH_AT_OLD_TIME R/V rec"),
    ("oldv_neg", "OLD_V must be positive"),
    ("oldv_flat", "OLD_V must be strictly decreasing"),
    ("v_neg", "V must be positive"),
    ("v_flat", "V must be strictly decreasing"),
    ("inner_mm", "Velocities at inner boundary are not identical"),
    ("outer_big", "Velocity at outer boundary is too large"),
    ("dt_zero", "0 ov -ve time step"),
]


def build():
    md = BINDIR / "obj_jh"
    objs = [
        prod_obj(TREE / "plane" / "get_jh_at_prev_time_step.f"),
        prod_obj(TREE / "subs" / "write_direct_info_v3.f"),
    ]
    objs += compile_objs(
        "jh_drv", [F90DIR / "drv_jh_prev.f90"],
        moddir=md, extra_flags=PROD_MOD_INC)
    libs = [TREE / "lib" / a for a in ("libsubs.a", "libmod.a")]
    return link_prod("jh", objs + libs, moddir=md)


def nj(out):
    return {
        int(m.group(1)): float(m.group(2))
        for m in re.finditer(r"NJ\s+(\d+)\s+(\S+)", out)
    }


def casedir(tag):
    d = WORKDIR / f"jh_{tag}"
    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True)
    return d


def main():
    rep = Report("GET_JH_AT_PREV_TIME_STEP alterations")
    binary = build()

    d = casedir("good")
    r = run(binary, "good", cwd=d)
    vals = nj(r.stdout)
    dt = float(re.search(r"DT\s+(\S+)", r.stdout).group(1))
    rep.check("good: NEW_J=51 on exact rec",
              all(abs(v - 51.0) < 1e-9 for v in vals.values()),
              detail=str(vals))
    rep.check("good: NEW_H finite", "NH" in r.stdout)
    rep.check("good: DT=1e4", abs(dt - 1.0e4) < 1.0, detail=str(dt))

    d = casedir("ogrey")
    r = run(binary, "ogrey", cwd=d)
    vals = nj(r.stdout)
    rep.check("ogrey: NEW_J=1.5 from grey rec",
              all(abs(v - 1.5) < 1e-9 for v in vals.values()),
              detail=str(vals))

    d = casedir("nu_high")
    r = run(binary, "nu_high", cwd=d)
    vals = nj(r.stdout)
    rep.check("nu_high: max-freq warning",
              "Invalid maximum frequency" in r.stdout)
    rep.check("nu_high: clamps to JST(:,1)=1",
              all(abs(v - 1.0) < 1e-9 for v in vals.values()))

    d = casedir("nu_low")
    r = run(binary, "nu_low", cwd=d)
    vals = nj(r.stdout)
    rep.check("nu_low: min-freq warning",
              "Invalid minimum frequency" in r.stdout)
    rep.check("nu_low: clamps to JST(:,100)=100",
              all(abs(v - 100.0) < 1e-6 for v in vals.values()))

    for tag, msg in STOP_CASES:
        d = casedir(tag)
        r = run(binary, tag, cwd=d, check=False)
        rep.check(f"{tag}: STOP '{msg[:28]}'", msg in r.stdout + r.stderr,
                  detail=(r.stdout + r.stderr).strip().splitlines()[-1]
                  if (r.stdout + r.stderr).strip() else "")

    rep.summary()
    sys.exit(0)


if __name__ == "__main__":
    main()
