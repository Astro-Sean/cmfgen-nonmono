#!/usr/bin/env python3
"""EVAL_ADIABATIC_V3 SIGMA_TRUE selection
(new_main/subs/eval_adiabatic_v3.f).

Alteration under test: the adiabatic coolant's B coefficient
uses the instantaneous velocity gradient.  Under the nonmono
flags SIGMA can be stale (the forced-monotone LOG(V) cubic
from ADJUST_R_GRID_V4), so SIG_AD selects SIGMA_TRUE when
allocated and falls back to SIGMA otherwise - stock runs are
bit-identical.

B(I) = SCALE*(POP_ATOM+ED)*V*(3+SIG_AD)/R and
AD_CR_V(I) = B(I)*T(I)*4e-10*PI, so the ratio
AD_CR_V_a/AD_CR_V_b isolates (3+SIG_AD) exactly.

  sob_true   : DO_NONMONO_SOB=T, SIGMA_TRUE allocated ->
               AD_CR_V scaled by (3+SIGMA_TRUE)/(3+SIGMA)
  knop_true  : DO_NONMONO_KNOP=T, same selection
  sob_unalloc: flag on but SIGMA_TRUE absent -> identical
               to flag-off (SIGMA fallback)
  off        : both flags F, SIGMA_TRUE present -> SIGMA
               still used (flag gate works)
"""
import re
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (BINDIR, F90DIR, PROD_MOD_INC, TREE, Report,
                    compile_objs, link_prod, prod_obj, run)

SIGMA = [-0.5] * 5
SIGMA_TRUE = [0.3, -0.8, 0.5, -0.2, 0.7]


def vec(out, key):
    m = re.search(rf"{key}=\s*(.*)", out)
    return [float(v) for v in m.group(1).split()] if m else None


def main():
    rep = Report("EVAL_ADIABATIC_V3 SIGMA_TRUE selection")
    md = BINDIR / "obj_adiab"
    objs = [prod_obj(TREE / "new_main" / "subs" /
                     "eval_adiabatic_v3.f")]
    objs += compile_objs(
        "adiab_drv", [F90DIR / "drv_adiab.f90"],
        moddir=md, extra_flags=PROD_MOD_INC)
    libs = [TREE / "lib" / a for a in ("libsubs.a", "libmod.a")]
    binary = link_prod("adiab", objs + libs, moddir=md)

    out = {}
    for case in ("sob_true", "knop_true", "sob_unalloc", "off"):
        with tempfile.TemporaryDirectory() as d:
            r = run(binary, case, cwd=d)
            out[case] = r.stdout
            rep.check(f"{case}: completes",
                      f"DONE_{case}" in r.stdout)

    ref_v = vec(out["off"], "ADCRV")
    sob_v = vec(out["sob_true"], "ADCRV")
    kn_v = vec(out["knop_true"], "ADCRV")
    un_v = vec(out["sob_unalloc"], "ADCRV")

    rep.check("off: AD_CR_V finite and nonzero",
              ref_v is not None and all(
                  v == v and v != 0.0 for v in ref_v))
    rep.check("sob: uses SIGMA_TRUE - ratio (3+ST)/(3+SIGMA)",
              sob_v is not None and all(
                  abs(s / r - (3.0 + st) / (3.0 + sm)) < 1e-5
                  for s, r, st, sm in
                  zip(sob_v, ref_v, SIGMA_TRUE, SIGMA)))
    rep.check("knop: same SIGMA_TRUE selection",
              kn_v is not None and all(
                  abs(k / r - (3.0 + st) / (3.0 + sm)) < 1e-5
                  for k, r, st, sm in
                  zip(kn_v, ref_v, SIGMA_TRUE, SIGMA)))
    rep.check("unalloc: SIGMA fallback identical to flag-off",
              un_v is not None and all(
                  abs(u - r) < 1e-6 * abs(r)
                  for u, r in zip(un_v, ref_v)))

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
