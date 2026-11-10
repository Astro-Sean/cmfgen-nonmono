#!/usr/bin/env python3
"""SET_ABUND_CLUMP non-monotonic VTURB guards
(new_main/mod_subs/set_abund_clump.f).

Alteration under test (05-Oct-2026): UPDATE_VTURB normalizes the
turbulent-width profile by a *speed*, not the signed radial
velocity. With DO_NONMONO_SOB/KNOP the POW law uses
ABS(V)/MAXVAL(ABS(V)) (all-zero V falls back to VTURB_MIN) and
TRUNC_POW uses MIN(ABS(V),VEND)/VEND. Upstream took V/V(1) and
MIN(V,VEND) - a sign reversal or repeated V made the ratio
non-positive and NaN'd the whole vector.

Flags off must preserve the upstream formula exactly: monotone V
gives the same power law, and a non-monotonic V still produces
the upstream NaN (the fix is deliberately flag-gated).

  pow_nm      : DO_NONMONO_SOB + signed/reversing V -> finite
                ABS-based profile, max at |V|max
  pow_nm_knop : same under DO_NONMONO_KNOP
  pow_off     : flags off, monotone V -> upstream formula result
  pow_off_nm  : flags off, signed V -> upstream NaN (unchanged
                legacy behaviour, proving the gate)
  pow_zero    : DO_NONMONO + all-zero V -> VTURB_MIN everywhere
  tpow_nm     : TRUNC_POW nonmono -> MIN(ABS,VEND)/VEND profile
  tpow_off    : TRUNC_POW flag-off -> upstream formula result
  bad_law     : unrecognized law STOPs
"""
import math
import re
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (BINDIR, F90DIR, PROD_MOD_INC, TREE, WORKDIR,
                    Report, compile_objs, link_prod, prod_obj, run)

VMIN, VMAX, VPOW, VEND = 5.0, 25.0, 0.5, 15.0
V_NM = [40.0, 10.0, -5.0, 20.0, 30.0]
V_MO = [40.0, 30.0, 20.0, 15.0, 10.0]


def build():
    md = BINDIR / "obj_vturb"
    objs = [prod_obj(TREE / "new_main" / "mod_subs"
                     / "set_abund_clump.f")]
    objs += compile_objs(
        "vt_drv", [F90DIR / "drv_vturb.f90"],
        moddir=md, extra_flags=PROD_MOD_INC)
    libs = [TREE / "lib" / a for a in
            ("libsubs.a", "libmod.a", "libsubs.a")]
    return link_prod("vturb", objs + libs, moddir=md)


def casedir(tag):
    d = WORKDIR / f"vt_{tag}"
    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True)
    return d


def vt(out):
    m = re.search(r"VT=\s+([-+0-9.EeNaNinf\s]+)", out)
    if not m:
        return None
    return [float(x) for x in m.group(1).split()]


def flag(out, key):
    m = re.search(rf"{key}=([TF])", out)
    return m.group(1) == "T" if m else None


def expect_abs(v):
    mx = max(abs(x) for x in v)
    return [VMIN + (VMAX - VMIN) * (abs(x) / mx) ** VPOW for x in v]


def expect_old(v):
    return [VMIN + (VMAX - VMIN) * (x / v[0]) ** VPOW for x in v]


def expect_tpow(v, use_abs):
    out = []
    for x in v:
        t = (abs(x) if use_abs else x)
        t = min(t, VEND) / VEND
        out.append(VMIN + (VMAX - VMIN) * t ** VPOW)
    return out


def close(a, b, tol=2e-4):
    return a is not None and len(a) == len(b) and all(
        abs(x - y) < tol for x, y in zip(a, b))


def run_case(binary, tag):
    d = casedir(tag)
    return run(binary, tag, cwd=d)


def main():
    rep = Report("SET_ABUND_CLUMP VTURB non-mono guards")
    binary = build()

    for tag in ("pow_nm", "pow_nm_knop"):
        r = run_case(binary, tag)
        rep.check(f"{tag}: completes", f"DONE_{tag}" in r.stdout)
        rep.check(f"{tag}: ABS-profile matches formula",
                  close(vt(r.stdout), expect_abs(V_NM)))
        rep.check(f"{tag}: all finite (no NaN)",
                  flag(r.stdout, "VT_FIN") is True
                  and flag(r.stdout, "VT_NAN") is False)

    r = run_case(binary, "pow_off")
    rep.check("pow_off: upstream formula preserved flag-off",
              close(vt(r.stdout), expect_old(V_MO)))

    r = run_case(binary, "pow_off_nm")
    rep.check("pow_off_nm: upstream NaN preserved flag-off",
              flag(r.stdout, "VT_NAN") is True)

    r = run_case(binary, "pow_zero")
    rep.check("pow_zero: all-zero V gives VTURB_MIN",
              vt(r.stdout) == [VMIN] * 5)

    r = run_case(binary, "tpow_nm")
    rep.check("tpow_nm: MIN(ABS,VEND)/VEND profile",
              close(vt(r.stdout), expect_tpow(V_NM, True)))

    r = run_case(binary, "tpow_off")
    rep.check("tpow_off: upstream MIN(V,VEND) formula",
              close(vt(r.stdout), expect_tpow(V_MO, False)))

    r = run_case(binary, "bad_law")
    rep.check("bad_law: unrecognized VTURB_LAW STOPs",
              "Error VTURB law not recognized" in r.stdout
              and "DONE_" not in r.stdout)

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
