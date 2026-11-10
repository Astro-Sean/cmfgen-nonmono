#!/usr/bin/env python3
"""PAR_OPAC_SUB non-monotonic V guard (opac/par_opac_sub.f).

Alteration under test: this standalone opacity/flux utility
uses the stock CMF formal solvers, which require a positive
single-valued velocity field.  The fork rejects a non-monotonic
or non-positive V at entry instead of evaluating an unsupported
formal solution; upstream ran into the machinery regardless.

  flat/reversal/neg : STOP with the labelled non-monotonic error
  good              : the guard does NOT fire (the run then dies
                      downstream on missing model state - that is
                      expected; only the guard's absence is
                      asserted)
"""
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (BINDIR, F90DIR, PROD_MOD_INC, TREE, Report,
                    compile_objs, link_prod, prod_obj, run)

MSG = "non-monotonic velocity fields are not supported"


def build():
    md = BINDIR / "obj_par_opac"
    objs = [prod_obj(TREE / "opac" / "mod_par_opac.f"),
            prod_obj(TREE / "opac" / "par_opac_sub.f"),
            prod_obj(TREE / "opac" / "ins_line_v5.f")]
    objs += compile_objs(
        "par_opac_drv",
        [F90DIR / "drv_par_opac.f90",
         F90DIR / "stub_par_opac_deps.f90"],
        moddir=md, extra_flags=PROD_MOD_INC)
    libs = [TREE / "lib" / a for a in ("libsubs.a", "libmod.a")]
    return link_prod("par_opac", objs + libs, moddir=md)


def main():
    rep = Report("PAR_OPAC_SUB non-monotonic V guard")
    binary = build()

    for case in ("flat", "reversal", "neg", "good"):
        with tempfile.TemporaryDirectory() as d:
            r = run(binary, case, check=False, timeout=60, cwd=d)
            out_flux = Path(d, "OUT_FLUX")
            text = (out_flux.read_text()
                    if out_flux.exists() else "")
            text += r.stdout + r.stderr
            rep.check(f"{case}: does not complete",
                      f"DONE_{case}" not in r.stdout)
            if case == "good":
                rep.check("good: guard does not fire",
                          MSG not in text)
            else:
                rep.check(f"{case}: labelled error in OUT_FLUX",
                          MSG in text)

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
