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
CTRL_MOD = (ROOT / "new_main/mod_subs/control_variable_mod.f").read_text()
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


def test_pop_precondition_gate():
    check(
        "GS_PRECOND key read with logical default",
        re.search(
            r"GS_PRECOND=L_FALSE\s+CALL RD_STORE_LOG\(GS_PRECOND,"
            r"'GS_PRECOND'",
            RD_CONTROL,
        )
        is not None,
    )
    check(
        "GS_PRECOND_MAX key read with 1.0E+03 default",
        re.search(
            r"GS_PRECOND_MAX=1\.0E\+03_LDP\s+CALL RD_STORE_DBLE"
            r"\(GS_PRECOND_MAX,'GS_PRECOND_MAX'",
            RD_CONTROL,
        )
        is not None,
    )
    check(
        "GS_PRECOND_LIM key read with 0 default",
        re.search(
            r"GS_PRECOND_LIM=0\s+CALL RD_STORE_INT\(GS_PRECOND_LIM,"
            r"'GS_PRECOND_LIM'",
            RD_CONTROL,
        )
        is not None,
    )
    check(
        "GS_PRECOND defaults T under DO_NONMONO_KNOP",
        "'GS_PRECOND'" in RD_CONTROL
        and "KNOP auto-set: GS_PRECOND=T" in RD_CONTROL,
    )
    check(
        "freeze override clears GS_FREEZE on high demand",
        re.search(
            r"GS_FREEZE \.AND\. GS_PRECOND \.AND\.[\s\S]*?"
            r"MAXCH \.GT\. GS_PRECOND_MAX[\s\S]*?"
            r"GS_FREEZE=\.FALSE\.",
            POPS,
        )
        is not None,
    )
    check(
        "override appends the pop_precond reason",
        "'pop_precond'" in POPS,
    )
    check(
        "preconditioned call is forced to LAMBDA",
        re.search(
            r"'pop_precond'[\s\S]{0,1200}?"
            r"LAMBDA_ITERATION=\.TRUE\.[\s\S]{0,200}?"
            r"T1=MAX_LAM_COR",
            POPS,
        )
        is not None,
    )
    check(
        "release-damp block excludes preconditioned calls",
        re.search(
            r"LAMBDA_ITERATION \.AND\. \.NOT\. MOM_TRIAL \.AND\.\s*\n"
            r"\s*\d?\s*\.NOT\. GS_PRECOND_ACT\)THEN",
            POPS,
        )
        is not None,
    )
    check(
        "n_precond written and read back via GS_STATE",
        "'n_precond='" in POPS
        and "READ(RSV_CH(11:),*,IOSTAT=IOS)GS_PC_RD" in POPS
        and "N_PRECOND=MAX(GS_PC_RD,0)" in POPS,
    )
    check(
        "preconditioned calls are not post-release",
        "IF(GS_FREEZE .OR. GS_PRECOND_ACT)THEN" in POPS,
    )


def test_bulk_exit_residual_floor():
    check(
        "KN_RAY_BULK_TOL DATA fallback is off (KNOP-scoped)",
        "DATA KN_RAY_BULK_TOL /0.0_LDP/" in FG,
    )
    check(
        "KN_RAY_BULK_TOL defaults 6.0E-2 under KNOP, 0 elsewhere",
        "KN_RAY_BULK_TOL_V=0.0_LDP" in RD_CONTROL
        and "IF(DO_NONMONO_KNOP)KN_RAY_BULK_TOL_V=6.0E-2_LDP"
        in RD_CONTROL
        and re.search(
            r"CALL SET_FG_KN_OPTS\(KN_RAY_TOL_V,KN_RAY_MAXPASS_V,\s*\n"
            r"\s*\d?\s*INEXT_RELAX_V,KN_RAY_BULK_TOL_V\)",
            RD_CONTROL,
        )
        is not None,
    )
    check(
        "both bulk exits carry the residual floor",
        FG.count("KN_RAY_BULK_TOL .LE. 0.0_LDP .OR.") >= 2
        and FG.count("RES .LE. KN_RAY_BULK_TOL))THEN") == 2,
    )
    check(
        "commit gate admits floored or unfloored bulk exits",
        re.search(
            r"IXIT_M \.EQ\. 2 \.AND\. \(KN_RAY_BULK_TOL \.LE\. "
            r"0\.0_LDP \.OR\.[\s\S]*?RES_M \.LE\. KN_RAY_BULK_TOL\)",
            FG,
        )
        is not None
        and re.search(
            r"IXIT_P \.EQ\. 2 \.AND\. \(KN_RAY_BULK_TOL \.LE\. "
            r"0\.0_LDP \.OR\.[\s\S]*?RES_P \.LE\. KN_RAY_BULK_TOL\)",
            FG,
        )
        is not None,
    )
    check(
        "exit census prints the bulk floor",
        "'  floor=',KN_RAY_BULK_TOL" in FG,
    )


def test_knop_observability_and_runaway_safety():
    check(
        "BA13-TRACE emits at lvl 1 under KNOP every 10th iteration",
        re.search(
            r"CONV_DIAG_LVL \.GE\. 2 \.OR\.\s*\n"
            r"\s*\d?\s*\(?DO_NONMONO_KNOP \.AND\. CONV_DIAG_LVL "
            r"\.GE\. 1 \.AND\.\s*\n"
            r"\s*\d?\s*MOD\(MAIN_COUNTER,10\) \.EQ\. 0",
            SOLVEBA,
        )
        is not None,
    )
    check(
        "release_no_progress reason is appended",
        "'release_no_progress'" in POPS,
    )
    check(
        "no-progress counter increments and prints in the GS dump",
        "N_REL_NOPROG=N_REL_NOPROG+1" in POPS
        and "'    no-progress releases" in POPS,
    )
    check(
        "release audit keys on demand fall and raw drift",
        "MAXCH .GE. 0.9_LDP*GS_REL_DEM" in POPS
        and "GS_SW_MAXCHG .GT. GS_UNSAFE_TOL" in POPS,
    )
    check(
        "release transition arms the audit",
        "GS_JUST_REL=.TRUE." in POPS and "GS_REL_DEM=MAXCH" in POPS,
    )
    check(
        "GS_FREEZE_LIM absent-key auto-set to 300 under KNOP",
        re.search(
            r"CALL GET_KEY_STRING\(TMP_STRING,'GS_FREEZE_LIM',"
            r"L_FALSE,[\s\S]*?' KNOP auto-set: GS_FREEZE_LIM=300'",
            RD_CONTROL,
        )
        is not None,
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


def test_no_progress_stops():
    check(
        "no-progress totals declared",
        "INTEGER*8, SAVE :: N_FROZEN_TOT=0" in POPS
        and "INTEGER*8, SAVE :: N_PRECOND_TOT=0" in POPS,
    )
    check(
        "frozen and preconditioned calls accumulate totals",
        "N_FROZEN_TOT=N_FROZEN_TOT+1" in POPS
        and "N_PRECOND_TOT=N_PRECOND_TOT+1" in POPS,
    )
    check(
        "release progress resets the totals",
        "N_FROZEN_TOT=0_8" in POPS and "N_PRECOND_TOT=0_8" in POPS,
    )
    check(
        "GS_FREEZE_LIM reads the no-progress total",
        "N_FROZEN_TOT .GE. GS_FREEZE_LIM" in POPS,
    )
    check(
        "GS_PRECOND_LIM reads the no-progress total",
        "N_PRECOND_TOT .GE. GS_PRECOND_LIM" in POPS,
    )
    check(
        "GS_REL_NOPROG_LIM stop exists",
        "N_REL_NOPROG .GE. GS_REL_NOPROG_LIM" in POPS
        and "INTEGER GS_REL_NOPROG_LIM" in CTRL_MOD
        and "RD_STORE_INT(GS_REL_NOPROG_LIM,'GS_REL_NOPROG_LIM'"
        in RD_CONTROL,
    )
    check(
        "GS_STATE persists the totals",
        "'n_frozen_tot='" in POPS
        and "'n_precond_tot='" in POPS
        and "'n_rel_noprog='" in POPS
        and "N_FROZEN_TOT=MAX(GS_NFT_RD" in POPS
        and "N_PRECOND_TOT=MAX(GS_PCT_RD" in POPS
        and "N_REL_NOPROG=MAX(GS_RNP_RD" in POPS,
    )
    check(
        "KNOP auto-sets the no-progress limits",
        "KNOP auto-set: GS_PRECOND_LIM=50" in RD_CONTROL
        and "KNOP auto-set: GS_REL_NOPROG_LIM=10" in RD_CONTROL,
    )


def test_fg_trial_churn_arm():
    check(
        "FG_TRIAL_CHURN declared in control_variable_mod",
        re.search(r"REAL\(KIND=LDP\) FG_TRIAL_CHURN", CTRL_MOD) is not None,
    )
    check(
        "FG_TRIAL_CHURN key read with 0.0 default",
        re.search(
            r"FG_TRIAL_CHURN=0\.0_LDP\s+CALL RD_STORE_DBLE"
            r"\(FG_TRIAL_CHURN,'FG_TRIAL_CHURN'",
            RD_CONTROL,
        )
        is not None,
    )
    check(
        "pre-churn gate snapshot taken",
        "GS_FG_PRECHURN=GS_FG_IMMATURE" in POPS,
    )
    check(
        "churn arm requires churn-only immaturity inside the factor",
        re.search(
            r"FG_CHURN_MARG=GS_FG_IMMATURE \.AND\. "
            r"\(\.NOT\. GS_FG_PRECHURN\)",
            POPS,
        )
        is not None
        and "FG_TRIAL_CHURN*MAX(50_8,GS_NTOT/1000)" in POPS,
    )
    check(
        "FG_TRIAL admits the churn arm",
        ".NOT. GS_FG_IMMATURE .OR. FG_CHURN_MARG" in POPS,
    )
    check(
        "churn arm tags its own reason",
        "'fg_trial_churn'" in POPS,
    )
    check(
        "FG_TRIAL_CHURN_MIN declared with 3 default",
        "INTEGER FG_TRIAL_CHURN_MIN" in CTRL_MOD
        and re.search(
            r"FG_TRIAL_CHURN_MIN=3\s+CALL RD_STORE_INT"
            r"\(FG_TRIAL_CHURN_MIN,'FG_TRIAL_CHURN_MIN'",
            RD_CONTROL,
        )
        is not None,
    )
    check(
        "churn-marginal streak counted and reset",
        re.search(
            r"IF\(FG_CHURN_MARG \.OR\. N_SAMECELL \.GE\. 2\)THEN\s+"
            r"N_CHURN_STREAK=N_CHURN_STREAK\+1\s+ELSE\s+"
            r"N_CHURN_STREAK=0",
            POPS,
        )
        is not None,
    )
    check(
        "same-cell LAMBDA streak updates from the dominant argmax",
        "BA_ARG_J .EQ. LST_CHG_J .AND. BA_ARG_I .EQ. LST_CHG_I" in POPS
        and "N_SAMECELL=N_SAMECELL+1" in POPS,
    )
    check(
        "same-cell streak waives the GS_NTOT>0 arm requirement",
        "GS_NTOT .GT. 0_8 .OR. N_SAMECELL .GE. 2" in POPS,
    )
    check(
        "churn arm gated on the streak not N_FROZEN",
        "N_CHURN_STREAK .GE. FG_TRIAL_CHURN_MIN" in POPS
        and re.search(
            r"\.NOT\. GS_FG_IMMATURE \.AND\. FG_TRIAL_MIN \.GT\. 0"
            r"[\s\S]*?N_FROZEN \.GE\. FG_TRIAL_MIN",
            POPS,
        )
        is not None,
    )


if __name__ == "__main__":
    test_heal_store_zeroes_fully_poisoned_rows()
    test_sweep_store_write_gated_on_new_freq()
    test_rollback_restore_validates_field_values()
    test_fail_closed_commit_defaults()
    test_fg_trial_excluded_from_population_hold()
    test_rel_damp_window_default_and_raw_release_counter()
    test_pop_precondition_gate()
    test_bulk_exit_residual_floor()
    test_knop_observability_and_runaway_safety()
    test_march_sign_flip_accumulator_is_signed()
    test_ordered_boundary_rows_skip_out_of_window_frequencies()
    test_march_only_called_from_sweep_fld()
    test_mom_outer_aa_history_rescaled_on_omega_change()
    test_mom_same_map_omega_estimator()
    test_solveba_global_temperature_cap_uses_max_dt_cor()
    test_mom_adapt_omega_v2_vadat_wiring()
    test_no_progress_stops()
    test_fg_trial_churn_arm()
    print("all KNOP lagged-field checks passed")
