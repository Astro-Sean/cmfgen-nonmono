#!/usr/bin/env python3
"""Static regression checks for the KNOP lagged-field drift fixes.

Reads the production Fortran sources as text and asserts the guards
landed: fail-closed refresh commit defaults, the NEW_FREQ sweep-write
gate, the fully-poisoned-row heal, the rollback-restore value audit,
the FG_TRIAL hold exclusion, and the raw-release counter wiring.
"""

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
FG = (ROOT / "plane/fg_j_cmf_nonmono_v1.f").read_text()
MOM = (ROOT / "plane/mom_j_cmf_nonmono_v1.f").read_text()
POPS = (ROOT / "new_main/mod_subs/solve_for_pops.f").read_text()
RD_CONTROL = (ROOT / "new_main/mod_subs/rd_control_variables.f").read_text()
SOLVEBA = (ROOT / "new_main/subs/solveba_v13.f").read_text()


def section(source, start, end):
    first = source.index(start)
    last = source.index(end, first)
    return source[first:last]


def check(name, passed):
    if not passed:
        raise AssertionError(name)
    print(f"[PASS] {name}")


def test_heal_store_zeroes_fully_poisoned_rows():
    body = section(FG, "SUBROUTINE HEAL_STORE", "RETURN\n\tEND")
    check(
        "HEAL_STORE no longer skips fully poisoned rows",
        "IF(NV .EQ. 0)CYCLE" not in body,
    )
    check(
        "HEAL_STORE zeroes the row when no sane slot exists",
        re.search(r"NV \.EQ\. 0\)THEN[\s\S]*?STORE\(KC,LS,LS2\)=0\.0_LDP", body)
        is not None,
    )


def test_sweep_store_write_gated_on_new_freq():
    block = section(FG, "Store intensities at depth grid points", "DO K=1,ND")
    check(
        "in-place I_NEXT write block is gated on NEW_FREQ",
        re.search(
            r"IF\(HAS_COMPRESSING \.AND\. NEW_FREQ \.AND\.", block
        )
        is not None,
    )


def test_rollback_restore_validates_field_values():
    body = section(
        FG,
        "SUBROUTINE RESTORE_FG_INEXT_LATEST",
        "END SUBROUTINE RESTORE_FG_INEXT_LATEST"
        if "END SUBROUTINE RESTORE_FG_INEXT_LATEST" in FG
        else "Observation-frame restore",
    )
    accept = body.index("OK=.TRUE.")
    audit = body.index("FG_GS_CAP_FILE .LT. 0")
    check(
        "rollback restore audits values before accepting",
        audit < accept
        and "I_P_DEPTH_NEXT .LT. 0.0_LDP" in body[:accept]
        and "I_M_DEPTH_NEXT .LT. 0.0_LDP" in body[:accept],
    )


def test_fail_closed_commit_defaults():
    check(
        "FG_COMMIT_CONV defaults to TRUE",
        "FG_COMMIT_CONV_V=.TRUE." in RD_CONTROL
        and "FG_COMMIT_CONV=.TRUE." in FG,
    )
    check(
        "FG_COMMIT_NC defaults to 0.05",
        "FG_COMMIT_NC_V=0.05_LDP" in RD_CONTROL
        and "FG_COMMIT_NC=0.05_LDP" in FG,
    )
    check(
        "FG_BULK_CONTR defaults to TRUE",
        "FG_BULK_CONTR_V=.TRUE." in RD_CONTROL
        and "FG_BULK_CONTR=.TRUE." in FG,
    )


def test_fg_trial_excluded_from_population_hold():
    check(
        "GS_HOLD_POPS condition excludes FG_TRIAL",
        re.search(
            r"GS_HOLD_POPS \.AND\. GS_FREEZE \.AND\. \.NOT\. MOM_TRIAL"
            r" \.AND\.\s*\n\s*\d?\s*\.NOT\. FG_TRIAL",
            POPS,
        )
        is not None,
    )


def test_rel_damp_window_default_and_raw_release_counter():
    check(
        "GS_REL_N defaults to 2",
        re.search(r"GS_REL_N=2\s+CALL RD_STORE_INT", RD_CONTROL) is not None,
    )
    check(
        "raw-cure releases counted into GS_STATE",
        "GS_RAW_RELEASES_TOT=GS_RAW_RELEASES_TOT+1" in POPS
        and "'n_raw_releases='" in POPS,
    )


def test_march_sign_flip_accumulator_is_signed():
    march = section(FG, "SUBROUTINE MARCH", "SUBROUTINE WRITE_STORE")
    check(
        "MARCH accumulates the signed update for the detector",
        "SGN_ACC=SGN_ACC+(VNEW-VOLD)/" in march,
    )
    check(
        "MARCH no longer feeds the ABS delta to the detector",
        "SGN_ACC=SGN_ACC+DCELL" not in march,
    )
    check(
        "RES/NBAD bookkeeping keeps the ABS delta",
        "DCELL=ABS(VNEW-VOLD)" in march,
    )


def test_ordered_boundary_rows_skip_out_of_window_frequencies():
    ord_body = section(
        FG, "SUBROUTINE SWEEP_FLD_ORD", "SUBROUTINE BUILD_FREQ"
    )
    row1 = section(
        ord_body, "Boundary row I=1", "CALL KN_ORD_SPW_FILL(FLD,1)"
    )
    row_ni = section(
        ord_body,
        "Inner boundary row I=NI",
        "CALL KN_ORD_SPW_FILL(FLD,NI)",
    )
    guard = "IF(HOTMODE .AND. (L .LT. W_LO .OR. L .GT. W_HI))CYCLE"
    check("row-1 boundary block has the hot-window guard", guard in row1)
    check(
        "row-NI boundary block has the hot-window guard",
        guard in row_ni,
    )


def test_march_only_called_from_sweep_fld():
    sweep = section(
        FG, "SUBROUTINE SWEEP_FLD(FLD", "SUBROUTINE KN_ORD_SPV"
    )
    check(
        "MARCH is called exactly once, inside SWEEP_FLD",
        FG.count("CALL MARCH") == 1 and "CALL MARCH" in sweep,
    )


def test_mom_outer_aa_history_rescaled_on_omega_change():
    body = section(
        MOM, "Adaptive JNEXT_RELAX", "MOM_RF_RESID_PRV=MOM_GS_MAXCHG"
    )
    check(
        "old JNEXT_RELAX captured before the smoothed update",
        "OMG_PRV=JNEXT_RELAX" in body,
    )
    check(
        "FHIST_OUT rescaled by the omega ratio after the update",
        re.search(
            r"FHIST_OUT\(1:ND_NM,1:KN_M_NCF,1:NHIST_OUT\)=\s*\n"
            r"\s*1\s+FHIST_OUT\(1:ND_NM,1:KN_M_NCF,1:NHIST_OUT\)\*\s*\n"
            r"\s*2\s+\(JNEXT_RELAX/OMG_PRV\)",
            body,
        )
        is not None,
    )
    check(
        "rescale is gated on a changed, positive old omega",
        "ABS(JNEXT_RELAX-OMG_PRV) .GT. 1.0E-6_LDP" in body
        and "OMG_PRV .GT. 0.0_LDP" in body,
    )


def test_mom_same_map_omega_estimator():
    check(
        "per-pass max residual appended to RES_PASS",
        "RES_PASS(IPASS)=PASS_CHG" in MOM,
    )
    body = section(
        MOM, "Adaptive JNEXT_RELAX", "MOM_RF_RESID_PRV=MOM_GS_MAXCHG"
    )
    check(
        "adjacent same-map ratios collected",
        "RHO_H(L)=RES_PASS(I)/RES_PASS(I-1)" in body,
    )
    check(
        "median-of-ratios selection used",
        "DTR=RHO_H((L+1)/2)" in body,
    )
    check(
        "cross-call ratio retained as the <3-ratio fallback",
        ".NOT. L_V2 .AND. MOM_RF_RESID_PRV" in body,
    )
    check(
        "V2 estimator gated on the MOM_ADAPT_OMEGA_V2 flag",
        "MOM_ADAPT_OMEGA_V2 .AND. ALLOCATED(RES_PASS)" in body,
    )


def test_solveba_global_temperature_cap_uses_max_dt_cor():
    body = section(
        SOLVEBA,
        "Only the SCALE_OPT='GLOBAL' branch",
        "Update the population levels",
    )
    check(
        "GLOBAL branch no longer hardcodes the 0.2 cap",
        "0.2_LDP/T3" not in body and "T3=0.2_LDP" not in body,
    )
    check(
        "GLOBAL branch caps on MAX_dT_COR",
        "T3=MAX_dT_COR" in body and "MAX_dT_COR/T3" in body,
    )


def test_mom_adapt_omega_v2_vadat_wiring():
    check(
        "MOM_ADAPT_OMEGA_V2 defaults to TRUE",
        "MOM_ADAPT_OMEGA_V2_V=.TRUE." in RD_CONTROL,
    )
    check(
        "MOM_ADAPT_OMEGA_V2 VADAT key is read",
        "'MOM_ADAPT_OMEGA_V2'" in RD_CONTROL,
    )
    check(
        "SET_MOM_OMEGA_V2 is wired to the refresh",
        "CALL SET_MOM_OMEGA_V2(MOM_ADAPT_OMEGA_V2_V)" in RD_CONTROL
        and "SUBROUTINE SET_MOM_OMEGA_V2" in MOM,
    )


if __name__ == "__main__":
    test_heal_store_zeroes_fully_poisoned_rows()
    test_sweep_store_write_gated_on_new_freq()
    test_rollback_restore_validates_field_values()
    test_fail_closed_commit_defaults()
    test_fg_trial_excluded_from_population_hold()
    test_rel_damp_window_default_and_raw_release_counter()
    test_march_sign_flip_accumulator_is_signed()
    test_ordered_boundary_rows_skip_out_of_window_frequencies()
    test_march_only_called_from_sweep_fld()
    test_mom_outer_aa_history_rescaled_on_omega_change()
    test_mom_same_map_omega_estimator()
    test_solveba_global_temperature_cap_uses_max_dt_cor()
    test_mom_adapt_omega_v2_vadat_wiring()
    print("all KNOP lagged-field checks passed")
