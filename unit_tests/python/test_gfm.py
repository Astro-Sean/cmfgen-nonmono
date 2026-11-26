#!/usr/bin/env python3
"""GENERATE_FULL_MATRIX_V3 JFNK_QUIET_IO gating
(new_main/subs/generate_full_matrix_v3.f).

Alteration under test: the BA_ASCI diagnostic dumps at K=5, K=41
and K=ND are now gated on .NOT. JFNK_QUIET_IO.  JFNK's assembly
paths re-enter this routine once per residual eval and Krylov
matvec -- an ungated dump rewrites ~100 MB of scratch every eval.
The production SOLVEBA writes are unchanged.

Cases (each in a fresh case dir):
  loud  : JFNK_QUIET_IO=F -> BA_ASCI_N_D5 and BA_ASCI_N_DND exist
  quiet : JFNK_QUIET_IO=T -> neither file exists
  ncdef : BA_NC_MAXPOP=F -> NC row at EQ_SPECIES (row 4)
  ncmax : BA_NC_MAXPOP=T -> NC row at the most populated level
          (row 3), displaced rate eq relocated to EQ_SPECIES
          (checklist WP-7c, drv_gfm_nc.f90)
"""
import shutil
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (SRC, BINDIR, F90DIR, PROD_MOD_INC, TREE, WORKDIR,
                    Report, compile_objs, link_prod, prod_obj, run)


def build():
    md = BINDIR / "obj_gfm"
    objs = [
        prod_obj(SRC / "new_main" / "subs"
                 / "generate_full_matrix_v3.f")]
    objs += compile_objs(
        "gfm_drv", [F90DIR / "drv_gfm.f90"],
        moddir=md, extra_flags=PROD_MOD_INC)
    libs = [TREE / "lib" / a for a in
            ("libsubs.a", "libmod.a", "libsubs.a")]
    return link_prod("gfm", objs + libs, moddir=md)


def casedir(tag):
    d = WORKDIR / f"gfm_{tag}"
    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True)
    return d


def build_nc():
    md = BINDIR / "obj_gfm"
    objs = [
        prod_obj(SRC / "new_main" / "subs"
                 / "generate_full_matrix_v3.f")]
    objs += compile_objs(
        "gfm_nc_drv", [F90DIR / "drv_gfm_nc.f90"],
        moddir=md, extra_flags=PROD_MOD_INC)
    libs = [TREE / "lib" / a for a in
            ("libsubs.a", "libmod.a", "libsubs.a")]
    return link_prod("gfm_nc", objs + libs, moddir=md)


def parse_block(out, tag, n):
    lines = out.splitlines()
    i = lines.index(tag)
    return np.array(
        [[float(x) for x in ln.split()] for ln in lines[i + 1 : i + 1 + n]]
    )


def check_nc_case(rep, exe, case):
    """Fixture: ion 1 has level rows 1 (g.s.), 2 and 3, NC at eq 4;
    pseudo ion 2 sits on EQ_SPECIES=4; ED/T at rows 5/6. POPS puts
    the argmax on row 3. Column j of C_MAT is scaled by
    POPS(j,K)=(1,5,10,0.5,2,2)."""
    d = casedir(case)
    out = run(exe, case, cwd=d).stdout
    rep.check(f"{case} completes", "DONE_" in out)
    cmat = parse_block(out, "C_MAT", 6)
    steq = parse_block(out, "STEQ_VEC", 1)[0]
    pops_col = np.array([1.0, 5.0, 10.0, 0.5, 2.0, 2.0])
    nc = np.array([92.0, 93.0, 94.0, 95.0, 0.0, 0.0]) * pops_col
    rate3 = np.array([32.0, 33.0, 34.0, 0.0, 0.0, 0.0]) * pops_col
    row1 = np.array([12.0, 13.0, 14.0, 0.0, 0.0, 0.0]) * pops_col
    row5 = np.full(6, 0.01) * pops_col
    row6 = np.full(6, 0.02) * pops_col
    if case == "ncmax":
        exp3, exp4 = nc, rate3
        exp_steq34 = (4.8, 1.3)
    else:
        exp3, exp4 = rate3, nc
        exp_steq34 = (1.3, 4.8)
    for tag, got, exp in (
        ("row1 rate eq", cmat[0], row1),
        ("NC/rate row 3", cmat[2], exp3),
        ("NC/rate row 4", cmat[3], exp4),
        ("ED row", cmat[4], row5),
        ("T row", cmat[5], row6),
    ):
        err = np.abs(got - exp).max()
        rep.check(f"{case}: {tag} placement", err < 1e-10,
                  f"maxerr={err:.2e}")
    err = max(abs(steq[2] - exp_steq34[0]), abs(steq[3] - exp_steq34[1]))
    rep.check(f"{case}: STEQ NC/rate placement", err < 1e-12,
              f"maxerr={err:.2e}")


def main():
    rep = Report("gfm")
    exe = build()

    for case, expect in (("loud", True), ("quiet", False)):
        d = casedir(case)
        out = run(exe, case, cwd=d).stdout
        rep.check(f"{case} completes", "DONE_" in out)
        for name in ("BA_ASCI_N_D5", "BA_ASCI_N_DND"):
            rep.check(
                f"{case}: {name} "
                f"{'written' if expect else 'suppressed'}",
                (d / name).exists() == expect)
        if expect:
            body = (d / "BA_ASCI_N_D5").read_text()
            rep.check("loud dump holds STEQ_VEC block",
                      "STEQ_VEC" in body or "C_MAT" in body)

    exe_nc = build_nc()
    check_nc_case(rep, exe_nc, "ncdef")
    check_nc_case(rep, exe_nc, "ncmax")

    rep.summary()


if __name__ == "__main__":
    main()
