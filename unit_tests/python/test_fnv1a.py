#!/usr/bin/env python3
"""Unit test for KN_FNV1A_I1 / KN_FNV1A_COEF (plane/fg_j_cmf_nonmono_v1.f).

Compiles the exact production file plus f90/drv_fnv1a.f90 (HELPERS
identical to test_knop_ray.py) and asserts the driver's self-checks:
FNV-1a published vectors, complete byte coverage of REAL(LDP),
local/compensating change detection, array independence, determinism,
dimension mixing, and a production-size timing line for the
incremental-cost report.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import (
    F90DIR,
    TREE,
    Report,
    compile_and_link_prod,
    run,
)

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
    return compile_and_link_prod(
        "fnv1a", [TREE / h for h in HELPERS] + [F90DIR / "drv_fnv1a.f90"]
    )


def main():
    rep = Report("FNV1A checksum")
    exe = build()
    r = run(exe, check=False, timeout=600)
    out = r.stdout
    fails = [ln for ln in out.splitlines() if "FNV1A FAIL" in ln]
    rep.check("driver exit status", r.returncode == 0)
    rep.check("no FNV1A FAIL lines", not fails, "; ".join(fails[:5]))
    rep.check(
        "result line reports zero failures",
        "FNV1A RESULT:" in out and " 0 failed" in out,
    )
    rep.check(
        "production-size timing reported",
        "FNV1A perf two prod-size arrays:" in out,
    )
    nfail = rep.summary()
    if nfail:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
