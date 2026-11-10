#!/usr/bin/env python3
"""DO_TAU_REGRID_V2 floored-opacity head fix
(new_main/subs/do_tau_regrid_v2.f vs the frozen reference).

Alteration under test: TAU_OLD(1) now uses the floored opacity
TA(1)=MAX(|FLUX_MEAN(1)|,ESEC(1))*CLUMP_FAC(1) instead of raw
FLUX_MEAN(1).  When FLUX_MEAN(1) is zero (e.g. a floored
|RLUMST| at the outer depth) upstream put LOG10(0)=-Inf into
TAU_OLD(1) and NaN'd the whole new tau grid; the fork already
used TA in the VTAU branch and now does so everywhere, plus a
1e-2*TAU_OLD(2) floor.

Driver: ND=8 grid, GRID_TYPE=UNIFORM via the option store,
FLUX_MEAN(1)=0 -> 'zero_head'.  The ref binary reads the same
inputs and either NaNs the grid or dies in MON_INTERP; prod
completes with a finite, ordered R grid and DONE_R_REV=T.
'normal' (nonzero FLUX_MEAN(1)) shows prod and ref agree.
"""
import re
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (BINDIR, F90DIR, PROD_MOD_INC, REF, TREE, Report,
                    compile_objs, link_prod, prod_obj, run)

OPTS = "UNIFORM   [GRID_TYPE]\n"


def build(tag, tau_src):
    md = BINDIR / f"obj_tau_{tag}"
    if tau_src.is_relative_to(TREE):
        objs = [prod_obj(tau_src)]
        drv_srcs = [F90DIR / "drv_tau_regrid.f90"]
    else:
        objs = []
        drv_srcs = [tau_src, F90DIR / "drv_tau_regrid.f90"]
    objs += [prod_obj(TREE / "newsubs" / "rd_store_log.f"),
             prod_obj(TREE / "subs" / "mon_interp.f"),
             prod_obj(TREE / "subs" / "nordtau.f"),
             prod_obj(TREE / "unix" / "term_io.f")]
    objs += compile_objs(
        f"tau_{tag}_drv", drv_srcs,
        moddir=md, extra_flags=PROD_MOD_INC)
    libs = [TREE / "lib" / a for a in ("libsubs.a", "libmod.a")]
    return link_prod(f"tau_{tag}", objs + libs, moddir=md)


def rvec(out):
    m = re.search(r"RNEW=\s*(.*)", out)
    return [float(v) for v in m.group(1).split()] if m else None


def main():
    rep = Report("DO_TAU_REGRID_V2 floored-opacity head")
    prod = build("prod", TREE / "new_main" / "subs" /
                 "do_tau_regrid_v2.f")
    ref = build("ref", REF / "new_main" / "subs" /
                "do_tau_regrid_v2.f")

    for case in ("normal", "zero_head"):
        with tempfile.TemporaryDirectory() as d:
            Path(d, "tau_opts.txt").write_text(OPTS)
            p = run(prod, case, check=False, cwd=d)
            r = run(ref, case, check=False, cwd=d)
            pr, rr = rvec(p.stdout), rvec(r.stdout)

        if case == "normal":
            rep.check("normal: prod completes",
                      f"DONE_{case}" in p.stdout)
            rep.check("normal: ref completes too",
                      f"DONE_{case}" in r.stdout)
            rep.check("normal: regrid applied",
                      "REV=T" in p.stdout)
            rep.check("normal: prod/ref R grids agree",
                      pr is not None and rr is not None and
                      all(abs(a - b) < 1e-6 * abs(b)
                          for a, b in zip(pr, rr)))
        else:
            rep.check("zero_head: prod completes",
                      f"DONE_{case}" in p.stdout)
            rep.check("zero_head: finite Tau(min) printed",
                      "-Infinity" not in p.stdout and
                      "Tau(min)" in p.stdout)
            rep.check("zero_head: new R grid finite+ordered",
                      pr is not None and
                      all(v == v and v > 0 for v in pr) and
                      all(pr[i] > pr[i + 1]
                          for i in range(len(pr) - 1)))
            rep.check("zero_head: endpoints preserved",
                      pr is not None and
                      abs(pr[0] - 20.0) < 1e-6 and
                      abs(pr[-1] - 2.0) < 1e-6)
            ref_died = ("DONE_zero_head" not in r.stdout or
                        rr is None or
                        any(v != v or v <= 0 for v in rr))
            rep.check("zero_head: ref NaNs/dies on LOG10(0)",
                      ref_died)

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
