#!/usr/bin/env python3
"""test_mom_aa_state -- audit-D2 GET_MOM_AA_STATE state/activity split.

Links the production MOD_MOM_J_NONMONO_V1 (plane/
mom_j_cmf_nonmono_v1.f) plus its real externals, drives the real
GET_MOM_AA_STATE getter through unit_tests/f90/drv_mom_aa_state.f90,
and asserts:

  A  refresh ran, both AAs enabled        -> flags + counts verbatim
  B  refresh ran, inner AA disabled       -> AA_ON=F, counts still pass
  C  refresh ran, outer AA disabled       -> AA_OUT_ON=F
  D  refresh did NOT run, stale counters  -> all counters 0, RF_RAN=F
  E  formal-J path (never ran)            -> RF_RAN=F, counters 0
  F  no refresh + outer alloc latch       -> AFAIL=T while counters 0
     and enables still reported (alloc fail != 'disabled')
  G  ran, enabled, zero applications      -> enables T with counts 0
     (enabled state must not be inferred from activity counters)

This exercises the production getter directly; nothing reimplements
its logic.

Run: python3 unit_tests/python/test_mom_aa_state.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import common
from common import F90DIR, TREE, WORKDIR, Report, run

REP = common.Report("test_mom_aa_state (D2 GET_MOM_AA_STATE)")

WORK = WORKDIR / "mom_aa_state"
WORK.mkdir(parents=True, exist_ok=True)

# Same production link set as test_knop_ray: resolves every external
# the mom solver object references so the real object links whole.
HELPERS = [
    "tools/set_kind_module.f",
    "unix/term_io.f",
    "unix/gen_asci_open.f",
    "unix/verbose_info.f",
    "subs/set_case.f",
    "subs/phys_con.f",
    "subs/mon_int_funs_v2.f",
    "subs/mon_interp.f",
    "subs/mod_trap_derivatives.f",
    "subs/nordtau.f",
    "subs/qfromf.f",
    "subs/integrate.f",
    "subs/write_vec.f",
    "unix/dir_acc_pars_gen.f",
    "plane/mod_ray_mom_store.f",
    "new_main/mod_subs/knop_status_mod.f90",
    "plane/get_moms_non_rel.f",
    "new_main/mod_subs/mod_cmfgen.f",
    "new_main/mod_subs/control_variable_mod.f",
    "new_main/mod_subs/jfnk_cmfgen_int_mod.f",
    "new_main/mod_subs/jfnk_solver_v1.f90",
    "new_main/mod_subs/check_param_consistency.f",
    "newsubs/rd_store_log.f",
    "new_main/mod_subs/iter_ledger_mod.f90",
    "plane/fg_j_cmf_nonmono_v1.f",
    "plane/mom_j_cmf_nonmono_v1.f",
    "unit_tests/f90/stub_vadat.f",
]


def build():
    return common.compile_and_link_prod(
        "momaa", [TREE / h for h in HELPERS] + [F90DIR / "drv_mom_aa_state.f90"]
    )


def parse(out):
    d = {}
    for ln in out.splitlines():
        if "=" in ln:
            k, v = ln.split("=", 1)
            d[k.strip()] = v.strip()
    return d


BIN = build()

CASES = {
    "a": dict(AA_ON="T", AA_OUT_ON="T", RF_RAN="T",
              AAPP="3", ASKIP="1", ASNG="0", ARST="2",
              AAPP_O="2", AREJ_O="1", ASNG_O="0",
              AFAIL="F", FREQM="F"),
    "b": dict(AA_ON="F", AA_OUT_ON="T", RF_RAN="T",
              AAPP="4", AAPP_O="1"),
    "c": dict(AA_ON="T", AA_OUT_ON="F", RF_RAN="T",
              AAPP="5", AAPP_O="0", AREJ_O="0", ASNG_O="0"),
    "d": dict(AA_ON="T", AA_OUT_ON="T", RF_RAN="F",
              AAPP="0", ASKIP="0", ASNG="0", ARST="0",
              AAPP_O="0", AREJ_O="0", ASNG_O="0"),
    "e": dict(AA_ON="T", AA_OUT_ON="T", RF_RAN="F",
              AAPP="0", AAPP_O="0"),
    "f": dict(AA_ON="T", AA_OUT_ON="T", RF_RAN="F", AFAIL="T",
              AAPP="0", AAPP_O="0"),
    "g": dict(AA_ON="T", AA_OUT_ON="T", RF_RAN="T",
              AAPP="0", ASKIP="0", ASNG="0", ARST="0",
              AAPP_O="0", AREJ_O="0", ASNG_O="0"),
}

for mode, exp in CASES.items():
    r = run(BIN, mode, cwd=WORK)
    got = parse(r.stdout)
    REP.check(f"driver ran ({mode})", r.returncode == 0,
              detail=r.stdout[-200:])
    for k, v in exp.items():
        REP.check(f"{mode}:{k}={v}", got.get(k) == v,
                  detail=f"got {got.get(k)}")

sys.exit(REP.summary())
