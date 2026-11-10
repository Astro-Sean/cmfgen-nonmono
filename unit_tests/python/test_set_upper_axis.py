#!/usr/bin/env python3
"""SET_UPPER_AXIS ED monotonicity guard (pgplt/set_upper_axis.f).

Alteration under test: the top-axis labels XV by inverting
XV(ED), which requires a positive, strictly increasing ED.
Upstream fed LOG10 of non-positive ED (NaN) or produced a
nonsense axis when ED was flat/decreasing.  The fork validates
ED first and STOPs with a labelled error.

  good : strictly increasing ED -> completes, NXED and the
         label come back through the TOPBORD common block
  flat : ED(I+1)=ED(I) -> STOP, 'strictly increasing'
  drop : ED decreases mid-array -> STOP
  neg  : ED(1)<0 -> STOP
  zero : ED(1)=0 -> STOP
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import F90DIR, TREE, Report, compile_and_link_prod, run

SOURCES = [
    TREE / "tools" / "set_kind_module.f",
    TREE / "subs" / "mon_interp.f",
    TREE / "unix" / "term_io.f",
    TREE / "pgplt" / "set_upper_axis.f",
    F90DIR / "drv_set_upper_axis.f90",
]


def main():
    rep = Report("SET_UPPER_AXIS ED guard")
    binary = compile_and_link_prod("drv_set_upper_axis", SOURCES)

    r = run(binary, "good", check=False)
    rep.check("good: completes", "DONE_good" in r.stdout)
    rep.check("good: NXED populated", "NXED=" in r.stdout)
    rep.check("good: label returned", "LAB=TEST" in r.stdout)

    for case in ("flat", "drop", "neg", "zero"):
        r = run(binary, case, check=False)
        rep.check(f"{case}: does not complete",
                  f"DONE_{case}" not in r.stdout)
        rep.check(f"{case}: labelled error",
                  "strictly increasing" in r.stdout or
                  "strictly increasing" in r.stderr)

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
