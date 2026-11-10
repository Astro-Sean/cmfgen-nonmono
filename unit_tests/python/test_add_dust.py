#!/usr/bin/env python3
"""ADD_ETA_DUST innermost-point bound
(plane/add_dust_eta.f).

Alteration under test: the DO WHILE scan for the first ray point
at or inside R(1) gained an L < RAY(IP)%NZ bound.  A tangent ray
whose impact radius stays above the inner boundary has every
R_RAY > R(1); upstream, L ran to NZ+1 and read R_RAY(NZ+1) out
of bounds.  With the bound, NS=1 and MON_INTERP maps the inner
boundary emissivity onto the single innermost point, then the
R^-4 dilution fills the outer points.

The driver uses R=[2,4,6,8] (R(1)=2 inner), NP=8 rays feeding
the DW(IP) quadrature weights, and two ray geometries (R_RAY
runs inner to outer, matching the R ordering):

  outside  : R_RAY=[2.5,3,4,6,8], all above R(1).  Guard stops
             at L=NZ=5 -> the run completes; NS=1 so only the
             outermost point is interpolated, and points 1..4
             carry interp*(R(1)/R_RAY(5))**4 plus base ETA.
  crossing : R_RAY=[2,3,4,6,8] touches the boundary at L=1 --
             the normal path interpolates at every point; the
             endpoints land exactly on R(1) and R(ND), letting
             the quadrature-weight accumulation be checked to
             the printed precision.
"""
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import F90DIR, TREE, Report, compile_and_link_prod, run

SOURCES = [
    TREE / "tools" / "set_kind_module.f",
    TREE / "plane" / "mod_space_grid_v2.f",
    TREE / "subs" / "mon_interp.f",
    TREE / "unix" / "term_io.f",
    TREE / "plane" / "add_dust_eta.f",
    F90DIR / "drv_add_dust.f90",
]

ETA = [1.0e-4, 2.0e-4, 3.0e-4, 4.0e-4, 5.0e-4]


def vec(out, key):
    m = re.search(rf"{key}=\s*(.*)", out)
    return [float(v) for v in m.group(1).split()] if m else None


def main():
    rep = Report("ADD_ETA_DUST innermost-point bound")
    binary = compile_and_link_prod("drv_add_dust", SOURCES)

    out = run(binary, "outside").stdout
    rep.check("outside: completes (upstream read R_RAY(NZ+1))",
              "DONE_outside" in out)
    ep = vec(out, "EPRAY")
    rep.check("outside: all points written", ep is not None and
              all(v > -1.0 for v in ep))
    if ep:
        scale = (2.0 / 8.0) ** 4
        base = ep[4] - ETA[4]
        rep.check("outside: single interp positive/finite",
                  base == base and base > 0.0)
        rep.check("outside: remaining points carry R^-4 fill",
                  all(abs((ep[i] - ETA[i]) - base * scale) < 1e-10
                      for i in range(4)))
    em = vec(out, "EMRAY")
    rep.check("outside: minus stream filled too",
              em is not None and all(v > -1.0 for v in em))

    out = run(binary, "crossing").stdout
    rep.check("crossing: completes", "DONE_crossing" in out)
    ep = vec(out, "EPRAY")
    em = vec(out, "EMRAY")
    exp_p1 = 9.948e-4
    exp_p4 = 5.48e-4
    exp_m1 = 2.4172e-3
    rep.check("crossing: inner end equals ETA_P(1) from weights",
              ep is not None and
              abs(ep[0] - ETA[0] - exp_p1) < 5e-9)
    rep.check("crossing: outer end equals ETA_P(4) from weights",
              ep is not None and
              abs(ep[4] - ETA[4] - exp_p4) < 5e-9)
    rep.check("crossing: interior interp finite",
              ep is not None and all(
                  (ep[i] - ETA[i]) == (ep[i] - ETA[i])
                  for i in range(5)))
    rep.check("crossing: minus stream inner end matches weights",
              em is not None and
              abs(em[0] - ETA[0] - exp_m1) < 5e-9)

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
