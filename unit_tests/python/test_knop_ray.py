#!/usr/bin/env python3
"""Runtime test for KNOP_RAY_SOLVE_V1 (plane/fg_j_cmf_nonmono_v1.f).

Drives the production per-ray refresh solver directly through the
FG_J_CMF_NONMONO_MOD_V1 module state on a manufactured core-ray
problem (NI=8, NCF=6, one compressing band at depths 3-5).  This is
the first runtime coverage of the refresh fixed-point solver - the
machinery previously had static checks only.

Compiles the exact production file (HELPERS recipe identical to
test_mom_sweep_table.py) plus f90/drv_knop_ray.f90 and asserts on
the printed census fields: exit class (NXIT), residual, store
writes, matching census (NMATCH/NMAX), and store finiteness.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import F90DIR, TREE, WORKDIR, Report, compile_objs, link, run

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
    "plane/fg_j_cmf_nonmono_v1.f",
    "plane/mom_j_cmf_nonmono_v1.f",
    "unit_tests/f90/stub_vadat.f",
]


def build():
    objs = compile_objs(
        "knray", [TREE / h for h in HELPERS] + [F90DIR / "drv_knop_ray.f90"]
    )
    return link("knray", objs)


def parse(out):
    d = {}
    store = {}
    import re

    for ln in out.splitlines():
        t = ln.split()
        if not t:
            continue
        if t[0] == "STORE":
            store[(int(t[1]), int(t[2]))] = (float(t[3]), float(t[4]))
            continue
        m = re.match(r"^(\w+)=\s*(\S+)", ln)
        if m:
            k, v = m.group(1), m.group(2)
            try:
                d[k] = float(v)
            except ValueError:
                d[k] = v
        elif t[0] == "RESULT_OK":
            d["RESULT_OK"] = True
    return d, store


def main():
    rep = Report("KNOP_RAY_SOLVE_V1")
    bin_ = build()

    # -- T1: alternating sweep on the manufactured compressing band.
    r = run(bin_, "conv", cwd=WORKDIR, check=False)
    d, store = parse(r.stdout)
    ok = (
        "RESULT_OK" in d
        and d.get("NXIT", -1) in (0, 1, 2)
        and d.get("RESID", 1e9) < 1.0e-1
        and d.get("NCOMP", 0) > 0
        and all(abs(v) < 1e30 for v in store.values() for v in v)
    )
    rep.check(
        "T1 alternating sweep converges on manufactured band",
        ok,
        f"NXIT={d.get('NXIT')} RESID={d.get('RESID'):.3e} "
        f"NPASS={d.get('NPASS')} NCOMP={d.get('NCOMP')} rc={r.returncode}",
    )

    # -- T2: ordered sweep (KN_RAY_ORDERED) solves the same field.
    r = run(bin_, "ord", cwd=WORKDIR, check=False)
    d2, store2 = parse(r.stdout)
    ok = "RESULT_OK" in d2 and d2.get("RESID", 1e9) < 1.0e-1
    rep.check(
        "T2 ordered sweep converges",
        ok,
        f"NXIT={d2.get('NXIT')} RESID={d2.get('RESID'):.3e} "
        f"NPASS={d2.get('NPASS')}",
    )

    # -- T3: ordered vs alternating committed stores agree to the
    #    solve tolerance (both march the same discrete operator).
    common = set(store) & set(store2)
    if common:
        dd = max(
            abs(store[k][j] - store2[k][j]) / max(store2[k][j], 1e-30)
            for k in common
            for j in (0, 1)
        )
    else:
        dd = 9.9e9
    rep.check(
        "T3 ord/alt committed stores agree",
        "RESULT_OK" in d2 and dd < 5.0e-2,
        f"max rel store diff {dd:.2e} over {len(common)} slots",
    )

    # -- T4: DO_NONMONO_FREQ_MATCH with drift > 1 bin must activate
    #    multi-bin matched slots (NMATCH>0 and NMAX>1).
    r = run(bin_, "match", cwd=WORKDIR, check=False)
    d3, _ = parse(r.stdout)
    ok = "RESULT_OK" in d3 and d3.get("NMATCH", 0) > 0 and d3.get(
        "NMAX", 1
    ) > 1
    rep.check(
        "T4 FREQ_MATCH activates multi-bin slots under strong drift",
        ok,
        f"NMATCH={d3.get('NMATCH')} NMAX={d3.get('NMAX')} "
        f"RESID={d3.get('RESID'):.3e}",
    )

    # -- T5: a non-converged solve must exit visibly, never as
    #    'converged': forced stall/cap via KN_RAY_TOL=1e-14, MAXPASS=4
    #    must produce NXIT=4 (cap) - or another failure class - and
    #    a residual that stays above the production tolerance.
    r = run(bin_, "cap", cwd=WORKDIR, check=False)
    d4, _ = parse(r.stdout)
    xit = d4.get("NXIT", -1)
    ok = "RESULT_OK" in d4 and xit in (3, 4, 5, 6) and d4.get(
        "RESID", 0.0
    ) > 1.0e-14
    rep.check(
        "T5 forced non-convergence exits visibly (not as converged)",
        ok,
        f"NXIT={xit} RESID={d4.get('RESID'):.3e} NPASS={d4.get('NPASS')}",
    )

    # -- T6: warm start near the fixed point must show a small first-
    #    pass residual and converge in very few passes.
    r = run(bin_, "fixed", cwd=WORKDIR, check=False)
    d5, _ = parse(r.stdout)
    ok = "RESULT_OK" in d5 and d5.get("NXIT", -1) == 1 and d5.get(
        "NPASS", 99
    ) <= 4
    rep.check(
        "T6 near-fixed-point seed converges immediately",
        ok,
        f"NXIT={d5.get('NXIT')} RESID={d5.get('RESID'):.3e} "
        f"NPASS={d5.get('NPASS')}",
    )

    rep.summary()


if __name__ == "__main__":
    main()
