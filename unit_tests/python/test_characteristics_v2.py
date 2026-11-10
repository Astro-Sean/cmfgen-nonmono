#!/usr/bin/env python3
"""CHARACTERISTICS_V2 tangent-clamp + ABS(dV) insertion
(plane/characteristics_v2.f vs the frozen reference).

Two alterations under test, compared directly against the
pristine reference source built with production flags:

1. Impact-parameter clamp: Z(I)=SQRT(R_GRID(I)**2-P(IP)**2)
   goes NaN when a tangent ray's P exceeds R_GRID(I) -- grid
   noise is enough.  The fork clamps the effective impact
   radius to 0.9999*R_GRID(I) and counts the event.  The ref
   binary NaNs on the nudged ray; the prod binary stays finite.

2. ABS(dV) on the fine-grid insertion test.  Along a tangent
   ray Z*V/R falls toward the tangent point even for a smooth
   monotonic V, so the upstream signed test dV > threshold
   skipped every negative-dV interval and under-refined the
   ray.  With the fix, NINS=|dV|/threshold insertions land on
   both signs: prod NZ totals exceed ref on the same inputs.

Grid follows misc/tst_char_v2.f: R(1)=20 outer to R(ND)=2
inner (geometric), core rays fill P < R(ND), each tangent ray
lands exactly on a shell.  The 'nudge' case pushes two P values
2e-4 above their shell so the clamp is exercised.
"""
import re
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (F90DIR, REF, TREE, Report, compile_and_link_prod,
                    run)

SHARED = [
    TREE / "tools" / "set_kind_module.f",
    TREE / "plane" / "mod_space_grid_v2.f",
    TREE / "subs" / "mon_interp.f",
    TREE / "unix" / "term_io.f",
    TREE / "plane" / "velocity_law.f",
    TREE / "plane" / "derivs.f",
    TREE / "plane" / "derivr.f",
    TREE / "plane" / "runge_kutta.f",
]
DRV = F90DIR / "drv_characteristics.f90"
GAM_SOURCES = [
    TREE / "tools" / "set_kind_module.f",
    TREE / "plane" / "mod_space_grid_v2_gam.f",
    TREE / "subs" / "mon_interp.f",
    TREE / "unix" / "term_io.f",
    TREE / "plane" / "velocity_law.f",
    TREE / "plane" / "derivs.f",
    TREE / "plane" / "derivr.f",
    TREE / "plane" / "runge_kutta.f",
]
DRV_GAM = F90DIR / "drv_characteristics_gam.f90"


def rays(out):
    info = {}
    for m in re.finditer(
            r"RAY\s*(\d+) NZ=\s*(\d+) BAD=\s*(\d+) MU=\s*"
            r"([-+0-9.Ee]+)\s+([-+0-9.Ee]+) RMIN=\s*([-+0-9.Ee]+)",
            out):
        info[int(m.group(1))] = (int(m.group(2)), int(m.group(3)),
                                 float(m.group(4)), float(m.group(5)),
                                 float(m.group(6)))
    return info


def main():
    rep = Report("CHARACTERISTICS_V2 clamp + ABS(dV)")
    prod = compile_and_link_prod(
        "drv_characteristics",
        SHARED + [TREE / "plane" / "characteristics_v2.f", DRV])
    ref = compile_and_link_prod(
        "drv_characteristics_ref",
        SHARED + [REF / "plane" / "characteristics_v2.f", DRV])

    with tempfile.TemporaryDirectory() as d:
        p_out = run(prod, "mono", cwd=d).stdout
        r_out = run(ref, "mono", cwd=d).stdout
        rep.check("mono: both binaries complete",
                  "DONE_mono" in p_out and "DONE_mono" in r_out)
        rep.check("mono: insertion flag set", "INS=T" in p_out)
        pi, ri = rays(p_out), rays(r_out)
        rep.check("mono: prod rays all finite",
                  all(v[1] == 0 for v in pi.values()))
        rep.check("mono: ref rays also finite",
                  all(v[1] == 0 for v in ri.values()))
        p_nz = int(re.search(r"NZTOT=\s*(\d+)", p_out).group(1))
        r_nz = int(re.search(r"NZTOT=\s*(\d+)", r_out).group(1))
        rep.check("mono: V~R - prod within 2% of ref (tangent "
                  "clamp perturbs only the last interval)",
                  abs(p_nz - r_nz) <= 0.02 * r_nz)
        rep.check("mono: |MU| <= 1+eps on all rays",
                  all(v[2] >= -1.001 and v[3] <= 1.001
                      for v in pi.values()))
        rep.check("mono: ray radius stays at/above impact",
                  all(v[4] > 0.0 for v in pi.values()))

        p_out = run(prod, "nudge", cwd=d).stdout
        r_out = run(ref, "nudge", cwd=d).stdout
        rep.check("nudge: prod completes", "DONE_nudge" in p_out)
        pi, ri = rays(p_out), rays(r_out)
        rep.check("nudge: clamped ray stays finite",
                  pi.get(11) is not None and pi[11][1] == 0)
        rep.check("nudge: all prod rays finite",
                  all(v[1] == 0 for v in pi.values()))
        rep.check("nudge: clamp keeps Z real - prod refines the "
                  "clamped interval, ref skips it",
                  pi.get(11) is not None and ri.get(11) is not None
                  and pi[11][0] > ri[11][0])
        rep.check("nudge: outer nudged ray also refined",
                  pi.get(25) is not None and ri.get(25) is not None
                  and pi[25][0] > ri[25][0])

        p_out = run(prod, "dip", cwd=d).stdout
        r_out = run(ref, "dip", cwd=d).stdout
        rep.check("dip: completes", "DONE_dip" in p_out)
        pi = rays(p_out)
        rep.check("dip: velocity reversal stays finite",
                  all(v[1] == 0 for v in pi.values()))
        p_nz = int(re.search(r"NZTOT=\s*(\d+)", p_out).group(1))
        r_nz = int(re.search(r"NZTOT=\s*(\d+)", r_out).group(1))
        rep.check("dip: ABS(dV) refines the reversal flanks ref "
                  "misses", p_nz > r_nz)

    gam_prod = compile_and_link_prod(
        "drv_characteristics_gam",
        GAM_SOURCES + [TREE / "plane" / "characteristics_v2_gam.f",
                       DRV_GAM])
    gam_ref = compile_and_link_prod(
        "drv_characteristics_gam_ref",
        GAM_SOURCES + [REF / "plane" / "characteristics_v2_gam.f",
                       DRV_GAM])
    with tempfile.TemporaryDirectory() as d:
        p_out = run(gam_prod, "mono", cwd=d).stdout
        r_out = run(gam_ref, "mono", cwd=d).stdout
        rep.check("gam mono: both complete",
                  "DONE_mono" in p_out and "DONE_mono" in r_out)
        pi = rays(p_out)
        rep.check("gam mono: all rays finite",
                  all(v[1] == 0 for v in pi.values()))
        p_nz = int(re.search(r"NZTOT=\s*(\d+)", p_out).group(1))
        r_nz = int(re.search(r"NZTOT=\s*(\d+)", r_out).group(1))
        rep.check("gam mono: prod within 2% of ref",
                  abs(p_nz - r_nz) <= 0.02 * r_nz)

        p_out = run(gam_prod, "nudge", cwd=d).stdout
        r_out = run(gam_ref, "nudge", cwd=d).stdout
        pi, ri = rays(p_out), rays(r_out)
        rep.check("gam nudge: all prod rays finite",
                  all(v[1] == 0 for v in pi.values()))
        rep.check("gam nudge: clamped ray refined vs ref",
                  pi.get(11) is not None and ri.get(11) is not None
                  and pi[11][0] > ri[11][0])

        p_out = run(gam_prod, "dip", cwd=d).stdout
        r_out = run(gam_ref, "dip", cwd=d).stdout
        pi = rays(p_out)
        rep.check("gam dip: reversal stays finite",
                  all(v[1] == 0 for v in pi.values()))
        p_nz = int(re.search(r"NZTOT=\s*(\d+)", p_out).group(1))
        r_nz = int(re.search(r"NZTOT=\s*(\d+)", r_out).group(1))
        rep.check("gam dip: ABS(dV) out-refines ref", p_nz > r_nz)

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
