#!/usr/bin/env python3
"""Static regression checks for the GS_DRIFT_DIAG field-drift
diagnostic and the FG_RF_DCHG stale-delta fix.

Reads the production Fortran sources as text and asserts:

- the VADAT key is declared with default F and parsed via
  RD_STORE_LOG;
- the diagnostic block in SOLVE_FOR_POPS separates the three drift
  quantities the headline channel conflates (raw sweep drift,
  refresh residual, retained-field delta) and prints the adaptive
  omega, residual history, argmax cell, censuses, and
  freeze/release state;
- the GET_FG_GS_RELAX / GET_MOM_GS_RELAX accessors export the
  module state the line needs;
- FG_RF_DCHG is reset at sweep entry (the stale-delta defect: it
  was only cleared inside the refresh, so refresh-skipped
  iterations reported the previous refresh's retained delta);
- the MOM_SW_* snapshot cannot go stale (MOM_SW_NTOT is zeroed at
  the FREQ_INDX=1 reset).

Static checks only: the metric loops live inside the solver sweeps
and cannot be driven standalone (see docs/audits/FIELD_DRIFT_AUDIT.md
section 9, correction P2).
"""

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
FG = (ROOT / "plane/fg_j_cmf_nonmono_v1.f").read_text()
MOM = (ROOT / "plane/mom_j_cmf_nonmono_v1.f").read_text()
POPS = (ROOT / "new_main/mod_subs/solve_for_pops.f").read_text()
RD_CONTROL = (ROOT / "new_main/mod_subs/rd_control_variables.f").read_text()
CTRL_MOD = (ROOT / "new_main/mod_subs/control_variable_mod.f").read_text()


def section(source, start, end):
    first = source.index(start)
    last = source.index(end, first)
    return source[first:last]


def check(name, passed):
    if not passed:
        raise AssertionError(name)
    print(f"[PASS] {name}")


def test_vadat_key_declared_and_parsed():
    check(
        "GS_DRIFT_DIAG declared LOGICAL in CONTROL_VARIABLE_MOD",
        re.search(r"LOGICAL\s+GS_DRIFT_DIAG", CTRL_MOD) is not None,
    )
    check(
        "GS_DRIFT_DIAG parsed via RD_STORE_LOG",
        re.search(
            r"GS_DRIFT_DIAG=\.FALSE\.\s*\n\s*CALL RD_STORE_LOG\("
            r"GS_DRIFT_DIAG,'GS_DRIFT_DIAG'",
            RD_CONTROL,
        )
        is not None,
    )


def test_diag_block_separates_raw_relaxed_retained():
    block = section(POPS, "GS_DRIFT_DIAG .AND.", "END BLOCK")
    check(
        "diagnostic gated on the lagged-store paths only",
        "DO_NONMONO_KNOP .OR. DO_SOB_FG_KNOP" in block,
    )
    check(
        "ray line prints sweep drift, refresh residual, retained delta",
        "GS_SW_MAXCHG,GS_MAXCHG,GS_RF_DCHG" in block,
    )
    check(
        "mom line prints sweep drift, maxchg, retained delta",
        "MOM_SW_MAXCHG,MOM_MAXCHG,MOM_RF_DCHG" in block,
    )
    check(
        "relaxation factors and residual histories reported",
        "GET_FG_GS_RELAX" in block and "GET_MOM_GS_RELAX" in block,
    )
    check(
        "argmax cells and churn censuses reported",
        "FG_AK,FG_AR,FG_AL" in block and "MOM_AI,MOM_AL" in block
        and "GS_NCONV,GS_NTOT" in block and "MOM_SW_NBAD50,MOM_SW_NTOT" in block,
    )
    check(
        "freeze/release state and demand reported",
        "GS_FREEZE,GS_JUST_REL" in block and "N_MATURE" in block
        and "MAXCH" in block and "LAMBDA_ITERATION" in block,
    )


def test_relax_accessors_export_module_state():
    fg_body = section(
        FG,
        "SUBROUTINE GET_FG_GS_RELAX",
        "RETURN\n\tEND",
    )
    check(
        "GET_FG_GS_RELAX exports INEXT_RELAX and KN_RF_RESID_PRV",
        "OMEGA=INEXT_RELAX" in fg_body and "RESID_PRV=KN_RF_RESID_PRV" in fg_body,
    )
    check(
        "GET_FG_GS_RELAX exports the in-sweep argmax indices",
        "ARGK=FG_GS_ARGK" in fg_body and "ARGR=FG_GS_ARGR" in fg_body
        and "ARGL=FG_GS_ARGL" in fg_body,
    )
    mom_body = section(
        MOM,
        "SUBROUTINE GET_MOM_GS_RELAX",
        "RETURN\n\tEND",
    )
    check(
        "GET_MOM_GS_RELAX exports JNEXT_RELAX and MOM_RF_RESID_PRV",
        "OMEGA=JNEXT_RELAX" in mom_body
        and "RESID_PRV=MOM_RF_RESID_PRV" in mom_body,
    )
    check(
        "GET_MOM_GS_RELAX exports refresh ran/conv and argmax",
        "RF_RAN=MOM_RF_RAN" in mom_body and "RF_CONV=MOM_RF_CONV" in mom_body
        and "ARGI=MOM_GS_ARGI" in mom_body and "ARGL=MOM_GS_ARGL" in mom_body,
    )


def test_fg_rf_dchg_reset_at_sweep_entry():
    init_block = section(FG, "IF(INIT)THEN", "LOG_R_EXT(1:ND_EXT)")
    check(
        "FG_RF_DCHG cleared in the per-iteration INIT reset block",
        "FG_RF_DCHG=0.0_LDP" in init_block,
    )
    check(
        "FG_SW_* sweep snapshot still cleared alongside",
        "FG_SW_MAXCHG=0.0_LDP" in init_block and "FG_SW_NCONV=0" in init_block,
    )
    mom_init = section(MOM, "MOM_GS_NTOT=INT(ND_NM,8)", "MOM_NONMM_CNT=0")
    check(
        "MOM_SW_NTOT and MOM_RF_DCHG cleared per iteration",
        "MOM_SW_NTOT=0" in mom_init and "MOM_RF_DCHG=0.0_LDP" in mom_init
        and "MOM_RF_RAN=.FALSE." in mom_init,
    )


def test_mom_sw_snapshot_not_stale():
    getter = section(MOM, "SUBROUTINE GET_MOM_GS_MATURITY", "END")
    check(
        "maturity getter falls back to live counters when no refresh ran",
        "IF(MOM_SW_NTOT .GT. 0)THEN" in getter
        and "SW_MAXCHG=MOM_GS_MAXCHG" in getter,
    )


def main():
    test_vadat_key_declared_and_parsed()
    test_diag_block_separates_raw_relaxed_retained()
    test_relax_accessors_export_module_state()
    test_fg_rf_dchg_reset_at_sweep_entry()
    test_mom_sw_snapshot_not_stale()
    print("test_gs_drift_diag: all checks passed")


if __name__ == "__main__":
    main()
