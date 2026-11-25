#!/usr/bin/env python3
"""Static regression checks for the per-iteration OUTGEN layout
refactor (ITER_LOG_FMT_MOD + CONV_DIAG_LVL gating).

Reads the production Fortran sources as text and asserts:

- ITER_LOG_FMT_MOD exists, is free-form-compiled into libmod.a, and
  exports the shared helpers (rule/section/kv/banner);
- every Makefile that builds a consumer file carries the module dep
  so the build order cannot regress;
- the parser-visible tokens downstream scripts grep keep their exact
  wording ('Current iteration count is', 'POPULATION SOLVE',
  'Luminosity of star', 'Time taken for iteration', 'Maximum changes
  as returned by SOLVEBA_V13 is', 'Average number of calls to
  FG_J_CMF is', 'Maximum % increase', 'Maximum % decrease',
  'DEC_VEC:'/'INC_VEC:', 'ITER_TYPE', 'Model Convergence:');
- the verbosity contract holds: banner/summary surfaces stay
  unconditional or >=1, per-call trace chatter is >=2;
- the RUNAWAY census is one aligned table (single header + separator
  + top-N rows + suppressed-count line), not a stream of ad-hoc
  record lines.

Static checks only: the writes live inside the iteration loop and
cannot be driven standalone.
"""

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CMF = (ROOT / "new_main/cmfgen_sub.f").read_text()
CJB = (ROOT / "new_main/mod_subs/comp_j_blank.f").read_text()
POPS = (ROOT / "new_main/mod_subs/solve_for_pops.f").read_text()
BA = (ROOT / "new_main/subs/solveba_v13.f").read_text()
ILF = (ROOT / "new_main/mod_subs/iter_log_fmt_mod.f90").read_text()
MK_MOD = (ROOT / "new_main/mod_subs/Makefile").read_text()
MK_MAIN = (ROOT / "new_main/Makefile").read_text()
MK_SUBS = (ROOT / "new_main/subs/Makefile").read_text()
NT_SE = (ROOT / "new_main/subs/non_therm/se_ba_non_therm_v3.f90").read_text()
NT_EL = (
    ROOT / "new_main/subs/non_therm/electron_non_therm_spec.f90"
).read_text()
RDBA = (ROOT / "new_main/subs/read_ba_data_v3.f").read_text()
FXP = (ROOT / "new_main/subs/fixpop_in_ba_v3.f").read_text()
FPC = (ROOT / "new_main/subs/fiddle_pop_corrections_v2.f").read_text()


def check(name, passed):
    if not passed:
        raise AssertionError(name)
    print(f"[PASS] {name}")


def test_module_exists_with_shared_helpers():
    check(
        "ITER_LOG_FMT_MOD declares module + implicit none",
        re.search(r"module iter_log_fmt_mod", ILF) is not None
        and re.search(r"implicit none", ILF) is not None,
    )
    for sym in (
        "ilf_rule",
        "ilf_section",
        "ilf_kv",
        "ilf_kvd",
        "ilf_kvi",
        "ilf_set_banner",
        "ilf_banner",
    ):
        check(f"public helper {sym} exported", f"public :: " in ILF and
              re.search(rf"public ::[^\n]*{sym}\b|\b{sym}\b", ILF)
              is not None and f"subroutine {sym}" in ILF)
    check(
        "section width pinned at 80 columns",
        re.search(r"ilf_lw = 80\b", ILF) is not None,
    )
    check(
        "banner reports mode -> mode, BA status, conv %, elapsed s",
        "' -> '" in ILF and "' | BA '" in ILF and "' | Conv '" in ILF
        and "' s'" in ILF,
    )
    check(
        "banner skipped when no state exported (staleness guard)",
        re.search(r"ilf_it /= it\s*\)\s*return", ILF) is not None,
    )


def test_makefile_wiring():
    check(
        "mod_subs/Makefile builds iter_log_fmt_mod.o free-form",
        "iter_log_fmt_mod.f90" in MK_MOD and "FFREE" in MK_MOD,
    )
    check(
        "mod_subs/Makefile MOD list archives iter_log_fmt_mod",
        "(iter_log_fmt_mod.o)" in MK_MOD,
    )
    check(
        "new_main/Makefile MODS includes iter_log_fmt_mod",
        "libmod.a(iter_log_fmt_mod.o)" in MK_MAIN,
    )
    check(
        "subs/Makefile gives solveba_v13.o the ILF dep",
        re.search(r"solveba_v13\.o\) : \$\(ILF\)", MK_SUBS) is not None,
    )


def test_consumers_use_module():
    for name, src in (
        ("cmfgen_sub.f", CMF),
        ("comp_j_blank.f", CJB),
        ("solve_for_pops.f", POPS),
        ("solveba_v13.f", BA),
    ):
        check(
            f"{name} USEs ITER_LOG_FMT_MOD",
            re.search(r"USE ITER_LOG_FMT_MOD", src) is not None,
        )


def test_parser_tokens_preserved():
    pairs = [
        ("' Current iteration count is'", CMF),
        ("' Time taken for iteration'", CMF),
        ("'Luminosity of star", CMF),
        ("POPULATION SOLVE", POPS),
        ("'Model Convergence:'", POPS),
        ("ITER_TYPE", POPS),
        ("Maximum changes as returned by SOLVEBA_V13 is", BA),
        ("' DEC_VEC: '", BA),
        ("' INC_VEC: '", BA),
        ("' Maximum % increase at depth '", BA),
        ("' Maximum % decrease at depth '", BA),
        ("'Average number of calls to FG_J_CMF is'", CJB),
    ]
    for token, src in pairs:
        check(f"token {token.strip()} preserved", token in src)


def test_round2_helpers_and_moment_table():
    check(
        "fixed-width scientific helper uses two-digit exponent first",
        "write (str, '(es11.4e2)') x" in ILF,
    )
    check(
        "fixed-width scientific helper has three-digit fallback",
        "write (str, '(es11.3e3)') x" in ILF,
    )
    check(
        "tag helper emits fixed level/source grammar",
        "'['//trim(lev)//']['//trim(src)//'] '" in ILF,
    )
    check(
        "issue census prints an empty iteration explicitly",
        "'ISSUES ('" in ILF,
    )
    check(
        "shared consecutive-message compressor exists",
        "subroutine ilf_repeat_write" in ILF
        and "subroutine ilf_repeat_flush" in ILF,
    )
    mom = Path("plane/mom_j_cmf_nonmono_v1.f").read_text()
    check(
        "moment refresh prints one pass-table header",
        mom.count("CALL ILF_TABLE_HEADER(LUER,TBLC,TBLW,INDENT=2)") == 1,
    )
    check(
        "moment refresh row uses the fixed-width numeric helper",
        "RES_W=ILF_ES(MOM_BLK_RESID)" in mom
        and "CHG_W=ILF_ES(PASS_CHG)" in mom,
    )
    check(
        "best-pass detail line is verbose-only",
        "IF(CONV_DIAG_LVL .GE. 3)WRITE(LUER,'(A,I4,A,A,A,I12,A,I12,A,ES12.4)')"
        in mom,
    )
    check(
        "XM rejection warnings go through the shared compressor",
        "CALL ILF_REPEAT_WRITE(LUER,'WARN','MOM',TRIM(XM_MSG))" in mom,
    )
    check(
        "GS sweep no longer repeats the stored-field sentence",
        "stored-field change reported in the moment-refresh summary" in mom,
    )
    check(
        "iteration end prints the issue census before timing",
        CMF.index("CALL ILF_ISSUE_REPORT(LUER)")
        < CMF.index("' Time taken for iteration'"),
    )


def test_banner_call_and_state_export():
    check(
        "CMFGEN_SUB prints the IT banner at the timing site",
        "CALL ILF_BANNER" in CMF,
    )
    check(
        "SOLVE_FOR_POPS exports mode/BA/conv for the banner",
        "CALL ILF_SET_BANNER(MODE_NOW,MODE_NEXT,BA_CH,CONV_NPCT," in POPS,
    )
    check(
        "banner state export sits inside the ITER_SUM block",
        POPS.index("CALL ILF_SET_BANNER")
        > POPS.index("Model Convergence:"),
    )


def test_verbosity_levels():
    check(
        "setup chatter (Zeroing BA) gated at level 2",
        "IF(CONV_DIAG_LVL .GE. 2)WRITE(6,*)'Zeroing BA matrices'" in CMF,
    )
    check(
        "DTDR report gated at level 2",
        re.search(
            r"IF\(CONV_DIAG_LVL \.GE\. 2\)\s*\n\s*1\s*WRITE\(LUER,\*\)"
            r"'The value of DTDR is :',DTDR",
            CMF,
        )
        is not None,
    )
    check(
        "non-thermal species dump gated at level 2",
        "DPTH_INDX .EQ. 1 .AND. CONV_DIAG_LVL .GE. 2" in NT_SE,
    )
    check(
        "electron non-therm enter/exit chatter gated at level 2",
        NT_EL.count("CONV_DIAG_LVL .GE. 2") >= 2,
    )
    check(
        "BA_READ_V3 progress chatter gated at level 2",
        RDBA.count("CONV_DIAG_LVL .GE. 2") >= 3,
    )
    check(
        "held-fixed notices gated at level 2",
        FXP.count("CONV_DIAG_LVL .GE. 2") >= 3,
    )
    check(
        "POPULATION SOLVE section prints at level >=1",
        re.search(
            r"IF\(CONV_DIAG_LVL \.GE\. 1\)THEN\s*\n\s*CALL ILF_SECTION"
            r"\(LUER,'POPULATION SOLVE'\)",
            POPS,
        )
        is not None,
    )
    check(
        "DEC/INC vectors print only at level 2 (trace noise)",
        re.search(
            r"IF\(CONV_DIAG_LVL \.GE\. 2\)THEN\s*\n\s*WRITE\(6,"
            r"'\(10X,10I10\)'\)",
            BA,
        )
        is not None,
    )
    check(
        "Major-species scale line gated at level >=1",
        "CONV_DIAG_LVL .GE. 1" in FPC
        and "minimum value of scale for Major species" in FPC,
    )


def test_runaway_table_layout():
    check(
        "RUNAWAY census keeps the 'RUNAWAY it=' header token",
        "' RUNAWAY it='" in BA,
    )
    check(
        "RUNAWAY header reuses the row field widths (A5/A30/A11/A2/A4)",
        "(1X,A5,2X,A5,1X,A30,4(1X,A11)," in BA,
    )
    check(
        "RUNAWAY row format is indexed fixed-width (I5,A30,4xES11.4)",
        "(1X,I5,2X,I5,1X,A30,4(1X,ES11.4)," in BA,
    )
    check(
        "suppressed-row count line exists and counts RW_NSUP-RW_NDONE",
        "RW_NSUP-RW_NDONE" in BA and "RW_NDONE=0" in BA,
    )
    check(
        "RUNAWAY block prints at level >=1 with topn kept",
        "POP_RUNAWAY_TOPN" in BA and "CONV_DIAG_LVL .GE. 1" in BA,
    )


def test_no_numerics_in_module():
    # Comments name SOLVE_FOR_POPS; the state pins look for actual
    # solver data structures, not routine names.
    check(
        "formatter stores no solver state (banner state only)",
        "STEQ(" not in ILF and "POPS(" not in ILF and "BA(" not in ILF,
    )


def main():
    test_module_exists_with_shared_helpers()
    test_makefile_wiring()
    test_consumers_use_module()
    test_parser_tokens_preserved()
    test_round2_helpers_and_moment_table()
    test_banner_call_and_state_export()
    test_verbosity_levels()
    test_runaway_table_layout()
    test_no_numerics_in_module()
    print("All iteration-output layout checks passed.")


if __name__ == "__main__":
    main()
