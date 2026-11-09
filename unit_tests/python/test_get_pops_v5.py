#!/usr/bin/env python3
"""GET_POPS_AT_PREV_TIME_STEP_V5 velocity-coordinate guards
(new_main/subs/get_pops_at_prev_time_step_v5.f).

Alteration under test (05-Oct-2026): before LOG(V) interpolation
and the R(V) remap, both the current V (MOD_CMFGEN) and the
previous-time-step TMP_V (OLD_MODEL_DATA) must be positive and
strictly decreasing. Non-monotonic or non-positive grids used to
feed LOG() and the interpolator silently.

The driver writes a real OLD_MODEL_DATA unformatted file
('21-Jul-2008' format with NUM_SPECIES_RD=0 so the species section
ends at EOF; SPECIES_BEG_ID=0 makes the post-read species loops
no-ops) with ND_OLD=6 against a current ND=5 - different grids, so
the MON_INTERP log-V and linear-V paths genuinely interpolate.

  good         : completes; OLD_R monotone with exact endpoints,
                 T row interpolated (exact at endpoints)
  v_neg        : V(2)<0     -> STOP 'Current V must be positive'
  v_flat       : V(3)=V(2)  -> STOP 'strictly decreasing'
  oldv_neg     : TMP_V(2)<0 -> STOP 'Previous-time-step V ... positive'
  oldv_flat    : equal TMP_V-> STOP 'Previous-time-step ... decreasing'
  v_ib_mismatch: TMP_V(ND_OLD)/=V(ND) -> STOP 'inner boundary unequal'
  outer_small  : TMP_V(1)<V(1) -> warning only (upstream), then
                 MON_INTERP still stops on the out-of-range
                 coordinate - the warning-only path cannot
                 actually complete
"""
import re
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (BINDIR, F90DIR, PROD_MOD_INC, TREE, WORKDIR,
                    Report, compile_objs, link_prod, prod_obj, run)


def build():
    md = BINDIR / "obj_get_pops_v5"
    objs = [
        prod_obj(TREE / "new_main" / "subs"
                 / "get_pops_at_prev_time_step_v5.f"),
        prod_obj(TREE / "new_main" / "subs"
                 / "read_seq_time_file_v1.f"),
    ]
    objs += compile_objs(
        "gp5_drv", [F90DIR / "drv_get_pops_v5.f90"],
        moddir=md, extra_flags=PROD_MOD_INC)
    libs = [TREE / "lib" / a for a in ("libsubs.a", "libmod.a")]
    return link_prod("get_pops_v5", objs + libs, moddir=md)


def casedir(tag):
    d = WORKDIR / f"gp5_{tag}"
    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True)
    return d


def val(out, key):
    m = re.search(rf"{key}=([TF])", out)
    return m.group(1) == "T" if m else None


def fvals(out, key):
    m = re.search(rf"{key}=\s*([-+0-9.Ee]+)\s+([-+0-9.Ee]+)", out)
    return (float(m.group(1)), float(m.group(2))) if m else None


def run_case(binary, tag):
    d = casedir(tag)
    return run(binary, tag, cwd=d)


def main():
    rep = Report("GET_POPS_AT_PREV_TIME_STEP_V5 V guards")
    binary = build()

    r = run_case(binary, "good")
    rep.check("good: completes on valid grids",
              "COMPLETED" in r.stdout)
    e = fvals(r.stdout, "OLD_R_END")
    rep.check("good: OLD_R endpoints exact",
              e is not None and abs(e[0] - 12.0) < 1e-10
              and abs(e[1] - 3.0) < 1e-10)
    rep.check("good: remapped OLD_R strictly decreasing",
              val(r.stdout, "OLD_R_MONO") is True)
    rep.check("good: interpolated T finite/positive",
              val(r.stdout, "T_FIN") is True)
    e = fvals(r.stdout, "T_END")
    rep.check("good: T endpoints exact on node-coincident grid",
              e is not None and abs(e[0] - 5.0) < 1e-10
              and abs(e[1] - 3.0) < 1e-10)
    rep.check("good: electron row zero with no species",
              val(r.stdout, "ED_ZERO") is True)

    r = run_case(binary, "v_neg")
    rep.check("v_neg: STOP on non-positive current V",
              "Current V must be positive" in r.stdout
              and "COMPLETED" not in r.stdout)
    rep.check("v_neg: depth and value reported",
              "Depth, V=" in r.stdout)

    r = run_case(binary, "v_flat")
    rep.check("v_flat: STOP on non-decreasing current V",
              "Current V must be strictly decreasing" in r.stdout
              and "COMPLETED" not in r.stdout)

    r = run_case(binary, "oldv_neg")
    rep.check("oldv_neg: STOP on non-positive old V",
              "Previous-time-step V must be positive" in r.stdout
              and "COMPLETED" not in r.stdout)

    r = run_case(binary, "oldv_flat")
    rep.check("oldv_flat: STOP on non-decreasing old V",
              "Previous-time-step V must be strictly decreasing"
              in r.stdout and "COMPLETED" not in r.stdout)

    r = run_case(binary, "v_ib_mismatch")
    rep.check("v_ib_mismatch: STOP on inner-boundary V mismatch",
              "Velocities at inner boundary are unequal" in r.stdout
              and "COMPLETED" not in r.stdout)

    r = run_case(binary, "outer_small")
    rep.check("outer_small: warning emitted",
              "Old velocity at outer boundary is too small"
              in r.stdout)
    rep.check("outer_small: MON_INTERP then rejects the range",
              "Error in MON_INTERP - values outside range"
              in r.stdout and "COMPLETED" not in r.stdout)

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
