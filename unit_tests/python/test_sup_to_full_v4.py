#!/usr/bin/env python3
"""SUP_TO_FULL_V4 MOD_HOLD_LEV_DIS gate
(new_main/mod_subs/sup_to_full_v4.f).

Alteration under test: the POPION recompute + COMP_LEV_DIS_BLK
call is wrapped in IF(.NOT. MOD_HOLD_LEV_DIS). With the hold flag
set, POPION and the level-dissolution constants keep their
incoming values (the caller has already computed them on a
consistent ion-density). Flag-off preserves upstream: POPION is
the sum of POPS rows whose Z_POP exceeds 0.01, then the
dissolution constants are recomputed.

The driver builds a two-ion synthetic atom (5 full levels in 2
super levels on ion 1, single ground member on ion 2 with
XzV_PRES=F), the same atom the V3B audit uses, so the
super->full conservation invariant is checked alongside the gate.

  hold_f : POPION = rows 1+2+4 summed (row 3 skipped: Z_POP<0.01),
           B_LEV_DIS recomputed (1e30 fill for DO_LEV_DIS=F),
           conservation holds
  hold_t : POPION keeps its 777 sentinel, B/A/X_LEV_DIS keep the
           888 sentinel, and the super->full mapping still runs
"""
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import F90DIR, TREE, Report, compile_and_link_prod, run

SOURCES = [
    TREE / "tools" / "set_kind_module.f",
    TREE / "new_main" / "mod_subs" / "mod_cmfgen.f",
    TREE / "new_main" / "mod_subs" / "mod_lev_dis_blk.f",
    TREE / "subs" / "occupation_prob.f",
    TREE / "subs" / "sup_to_full_v3b.f",
    TREE / "new_main" / "mod_subs" / "sup_to_full_v4.f",
    F90DIR / "drv_sup_to_full_v4.f90",
]

POPS_SUM = [1.0e-3 + 4.0e-5 + 7.0e-2, 2.0e-4 + 8.0e-3 + 3.0e-2]


def f2(out, key):
    m = re.search(rf"{key}=\s*([-+0-9.Ee]+)\s+([-+0-9.Ee]+)", out)
    return (float(m.group(1)), float(m.group(2))) if m else None


def f1(out, key):
    m = re.search(rf"{key}=\s*([-+0-9.Ee]+)", out)
    return float(m.group(1)) if m else None


def main():
    rep = Report("SUP_TO_FULL_V4 MOD_HOLD_LEV_DIS gate")
    binary = compile_and_link_prod("drv_sup_to_full_v4", SOURCES)

    out = run(binary, "hold_f").stdout
    rep.check("hold_f: completes", "DONE_hold_f" in out)
    p = f2(out, "POPION")
    rep.check("hold_f: POPION = sum of Z_POP>0.01 rows",
              p is not None and
              abs(p[0] - POPS_SUM[0]) < 1e-10 and
              abs(p[1] - POPS_SUM[1]) < 1e-10)
    b = f2(out, "BLD")
    rep.check("hold_f: dissolution constants recomputed",
              b is not None and b[0] == 1.0e30)
    rep.check("hold_f: super->full conservation",
              (f1(out, "CONS_ERR") or 1.0) < 1e-11)

    out = run(binary, "hold_t").stdout
    rep.check("hold_t: completes", "DONE_hold_t" in out)
    p = f2(out, "POPION")
    rep.check("hold_t: POPION keeps incoming values",
              p == (777.0, 777.0))
    b = f2(out, "BLD")
    rep.check("hold_t: dissolution constants untouched",
              b == (888.0, 888.0))
    rep.check("hold_t: super->full mapping still runs",
              (f1(out, "CONS_ERR") or 1.0) < 1e-11)

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
