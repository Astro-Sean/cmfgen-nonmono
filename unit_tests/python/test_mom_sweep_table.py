#!/usr/bin/env python3
"""Runtime check for the MOM_J_CMF_NONMONO_V1 GS-sweep report.

The pre-fix 'matrix'/'eq resid' formats had more arg-consuming
descriptors than arguments, so MOM_MIN_DRA and MOM_MAX_COH (both
real) landed on A descriptors and printed garbage every report.
A static text check cannot see this; the test compiles
plane/mom_j_cmf_nonmono_v1.f, runs DUMP_MOM_JNEXT_GS through
drv_mom_sweep.f90, and requires the seeded values to appear on
their rows.
"""

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import (
    F90DIR,
    TREE,
    WORKDIR,
    Report,
    compile_objs,
    link,
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
    "plane/fg_j_cmf_nonmono_v1.f",
    "plane/mom_j_cmf_nonmono_v1.f",
]


def main():
    rep = Report("MOM_J_CMF_NONMONO_V1 GS-sweep table")
    objs = compile_objs(
        "mom_sweep",
        [TREE / h for h in HELPERS]
        + [F90DIR / "stub_vadat.f", F90DIR / "drv_mom_sweep.f90"],
    )
    exe = link("mom_sweep", objs)
    p = run(exe, cwd=WORKDIR)
    rep.check(
        "driver ran the sweep report",
        "MOM_J_CMF_NONMONO_V1: GS sweep" in p.stdout,
        p.stdout[-400:] if p.returncode != 0 else "",
    )
    matrix = re.search(r"^\s*matrix\s.*$", p.stdout, re.M)
    rep.check(
        "matrix row prints both DR ratios",
        matrix is not None
        and "1.100E-03" in matrix.group(0)
        and "4.500E-06" in matrix.group(0),
        p.stdout[-400:] if matrix is None else "",
    )
    eqres = re.search(r"^\s*eq resid\s.*$", p.stdout, re.M)
    rep.check(
        "eq resid row prints residual and ESEC/CHI-1",
        eqres is not None
        and "9.900E-02" in eqres.group(0)
        and "1.300E+01" in eqres.group(0),
        p.stdout[-400:] if eqres is None else "",
    )
    rep.summary()


if __name__ == "__main__":
    main()
