#!/usr/bin/env python3
"""CHECK_IONS_PRESENT depth-index default
(new_main/subs/check_ions_present.f).

Alteration under test: DPTH_INDX and ED_VAL are now initialised to
the inner boundary (ND) before the argmax loop over T. Upstream
left both undefined when no T(K) exceeded the initial T_VAL=0 -
an all-zero or all-negative T array then indexed T/ED/POP_SPECIES
with a garbage depth in the summary block below.

The driver sets every SPECIES_PRES=F so all species loops are
no-ops and the routine runs straight to the T_VAL argmax and the
summary prints on WARNING_LU (fort.4). The printed temperature is
T(DPTH_INDX), which exposes which depth index was used.

  zero_t : T=0 everywhere -> defaults to ND, prints T(ND)=0
  pos_t  : max T at depth 2 -> prints T=3.0 (argmax respected)
  neg_t  : T<0 everywhere -> defaults to ND, prints T(ND)=-1
"""
import re
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (BINDIR, F90DIR, PROD_MOD_INC, TREE, WORKDIR,
                    Report, compile_objs, link_prod, prod_obj, run)


def build():
    md = BINDIR / "obj_check_ions"
    objs = [prod_obj(TREE / "new_main" / "subs"
                     / "check_ions_present.f")]
    objs += compile_objs(
        "ci_drv", [F90DIR / "drv_check_ions.f90"],
        moddir=md, extra_flags=PROD_MOD_INC)
    libs = [TREE / "lib" / a for a in
            ("libsubs.a", "libmod.a", "libsubs.a")]
    return link_prod("check_ions", objs + libs, moddir=md)


def casedir(tag):
    d = WORKDIR / f"ci_{tag}"
    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True)
    return d


def warn_temp(d):
    f = d / "fort.4"
    if not f.exists():
        return None
    m = re.search(
        r"Temperature \(10\^4 K\) is\s+([-+0-9.Ee]+)",
        f.read_text())
    return float(m.group(1)) if m else None


def main():
    rep = Report("CHECK_IONS_PRESENT depth default")
    binary = build()

    for tag, want in (("zero_t", 0.0), ("pos_t", 3.0),
                      ("neg_t", -1.0)):
        d = casedir(tag)
        r = run(binary, tag, cwd=d)
        rep.check(f"{tag}: routine completes",
                  f"DONE_{tag}" in r.stdout)
        got = warn_temp(d)
        rep.check(
            f"{tag}: summary uses depth index T={want}",
            got is not None and abs(got - want) < 1e-6,
            detail=f"got {got}")
        rep.check(f"{tag}: fort.4 written (WARNING_LU live)",
                  (d / "fort.4").exists())

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
