#!/usr/bin/env python3
"""Super-level -> full-atom mapping audit (SUP_TO_FULL_V3B).

Compiles the production subs/sup_to_full_v3b.f against a synthetic-atom
driver and checks the invariants the line-opacity star ratios depend
on:

  A. sum over full levels of a super level = super-level population
     (exact conservation, non-LTE input).
  B. members share the super-level population in proportion to
     C2LTE_F_ON_S (same departure coefficient per group, INT_SEQ=0).
  C. LTE input (C2_S = C2LTE_S) returns LTE full-atom populations.
  D. a zero-population super level yields zero members, no NaN
     (regression for the unguarded C2_S/SUM division; upstream had
     the same 0/0 hazard, fixed here).
  E. extreme population contrast still conserves and stays finite.
"""

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from common import F90DIR, TREE, Report, compile_and_link_prod, run


SOURCES = [
    TREE / "tools" / "set_kind_module.f",
    TREE / "new_main" / "mod_subs" / "mod_lev_dis_blk.f",
    TREE / "subs" / "occupation_prob.f",
    TREE / "subs" / "sup_to_full_v3b.f",
    F90DIR / "drv_sup_to_full.f90",
]


def main():
    rep = Report("sup_to_full mapping")
    binary = compile_and_link_prod("drv_sup_to_full", SOURCES)
    out = run(binary).stdout

    for tag in ("A", "B", "C", "D", "E"):
        m = re.search(rf"^\s*{tag}\s.*err=\s*(\S+)", out, re.M)
        rep.check(m is not None, f"case {tag} reported", out)
        if m:
            rep.check(float(m.group(1).replace("D", "E").replace("d", "e"))
                      < 1.0e-8, f"case {tag} invariant err={m.group(1)}")
    rep.check("FAIL_COUNT= 0" in out and "PASS" in out, "driver PASS", out)
    sys.exit(rep.summary())


if __name__ == "__main__":
    main()
