"""Regression tests for the KN_ORD_P2DAMP state machine (phase 2.29A).

Covers R1-R10, R14 of the phase plan: the stale-state fix (insignificant
cells clear P2_SGN/P2_FLIP/P2_DMP on the ELSE branch at all three march
call sites) plus the KN_P2_TRACK/KN_P2_CLEAR/KN_P2_MID primitive
semantics. R11 (off-state bitwise equivalence) is covered by the
existing test_knop_ray.py suite, which asserts driver outputs against
independent expectations with KN_ORD_P2DAMP off.
"""

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import SRC, F90DIR, Report, compile_and_link_prod, run

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
    "new_main/mod_subs/iter_log_fmt_mod.f90",
    "plane/fg_j_cmf_nonmono_v1.f",
    "plane/mom_j_cmf_nonmono_v1.f",
]


def main():
    rep = Report("knop_p2_state")
    exe = compile_and_link_prod(
        "knp2", [SRC / h for h in HELPERS]
        + [F90DIR / "stub_vadat.f", F90DIR / "drv_knop_p2.f90"]
    )
    out = run(exe).stdout
    
    lines = [l for l in out.splitlines() if l.strip()]
    for l in lines:
        print("   ", l)

    # each row: TAG sgn flip dmp bool -> bool must be T
    for l in lines:
        parts = l.split()
        tag = parts[0]
        ok = parts[-1] == "T"
        rep.check(f"driver {tag}", ok, l)

        src = (SRC / "plane/fg_j_cmf_nonmono_v1.f").read_text()

    # three tracker call sites, each inside a significance block
    n_track = len(re.findall(r"CALL KN_P2_TRACK", src))
    n_clear = len(re.findall(r"CALL KN_P2_CLEAR", src))
    rep.check("R6/R7: three TRACK call sites", n_track == 3,
              f"found {n_track}")
    rep.check("R6/R7: three CLEAR call sites", n_clear == 3,
              f"found {n_clear}")

    # every CLEAR must sit on the ELSE branch of a significance IF
    sig_re = re.compile(
        r"IF\(MAX\(VOLD,ABS\(VNEW\)\) .GT. 1\.0E-4_LDP\*VMAXF\(L\)\)THEN"
        r"[\s\S]{0,3000}?ELSE\s*\n[\s\S]{0,200}?"
        r"CALL KN_P2_CLEAR")
    rep.check("CLEAR on significance ELSE branches",
              len(sig_re.findall(src)) == 3,
              f"found {len(sig_re.findall(src))}")

    # the NI boundary row indexes the NI arrays
    rep.check("boundary row tracks P2_*(NI,L)",
              "KN_P2_TRACK(P2_SGN(NI,L)" in src
              and "KN_P2_CLEAR(P2_SGN(NI,L)" in src)

    # alloc failure keeps DO_P2 off; state is per-solve allocatable
    rep.check("R10: alloc failure disables P2",
              re.search(r"ALLOCATE\(P2_SGN\(NI,NCF\).*?STAT=IOS", src,
                        re.S) is not None
              and "DO_P2=.FALSE." in src)

    # deallocation on exit -> no cross-solve state (R13)
    rep.check("R13: state deallocated on solve exit",
              "IF(DO_P2)DEALLOCATE(P2_SGN,P2_FLIP,P2_DMP)" in src)

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
