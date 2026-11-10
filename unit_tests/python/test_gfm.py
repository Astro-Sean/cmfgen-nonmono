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
"""
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (BINDIR, F90DIR, PROD_MOD_INC, TREE, WORKDIR,
                    Report, compile_objs, link_prod, prod_obj, run)


def build():
    md = BINDIR / "obj_gfm"
    objs = [
        prod_obj(TREE / "new_main" / "subs"
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

    rep.summary()


if __name__ == "__main__":
    main()
