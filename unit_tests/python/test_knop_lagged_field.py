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


def test_sweep_store_write_snapshot_relax():
    block = section(
        FG,
        "Store intensities at depth grid points",
        "Code-audit diagnostic",
    )
    check(
        "I_NEXT write block is not gated on NEW_FREQ (every call writes)",
        re.search(r"IF\(HAS_COMPRESSING \.AND\. NEW_FREQ", block) is None,
    )
    check(
        "stored column is snapshotted on NEW_FREQ or slot shift",
        "NEW_FREQ .OR. SNP_FREQ .NE. FREQ_IDX_ITER" in block,
    )
    check(
        "write relaxes against the snapshot, not the live store",
        "T2=SNP_P(KC,LS)" in block and "T2=SNP_M(KC,LS)" in block,
    )
    check(
        "churn metrics stay gated on NEW_FREQ",
        "IF(NEW_FREQ)FG_SW_CAPTOT" in block
        and block.count("IF(NEW_FREQ)THEN") >= 2,
    )
    check(
        "cap census stays gated on NEW_FREQ",
        "IF(NEW_FREQ .AND." in block,
    )


def test_commit_jump_veto_seed_waiver():
    check(
        "FG_SEED_COMMIT defaults to TRUE",
        "FG_SEED_COMMIT_V=.TRUE." in RD_CONTROL
        and "FG_SEED_COMMIT=.TRUE." in FG,
    )
    check(
        "per-ray seeded flag declared",
        "FG_RAY_SEEDED(:)" in FG,
    )
    check(
        "jump veto skipped while the ray is unseeded",
        re.search(
            r"L_UNSEED=\(FG_RAY_SEEDED\(LS\) \.EQ\. 0\)[\s\S]{0,200}?"
            r"IF\(\.NOT\. L_UNSEED\)JUMP=",
            FG,
        )
        is not None,
    )
    check(
        "WRITE_STORE commit marks the ray seeded",
        re.search(
            r"CALL WRITE_STORE\(FLD_P,I_P_DEPTH_NEXT,DCH_P\)[\s\S]{0,300}?"
            r"FG_RAY_SEEDED\(LS\)=1",
            FG,
        )
        is not None,
    )
    check(
        "validated restore marks all rays seeded",
        re.search(
            r"I_NEXT_GS_OK=\.NOT\. SEED_FOREIGN[\s\S]{0,300}?"
            r"FG_RAY_SEEDED=1",
            FG,
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


def test_mom_sweep_census_counts_first_call_only():
    """The comp_j_blank INACCURATE loop re-calls MOM_J_CMF_NONMONO_V1
    at the same frequency with NEW_FREQ=F, so ungated counters counted
    each cell ~14x per sweep (drift census exceeded NTOT, and
    MOM_QUAR_AGE crossed MOM_QUAR_PERSIST inside one sweep). The
    sweep census belongs to the first (fresh-solve) call only.
    """
    block = section(
        MOM,
        "Census scope: the comp_j_blank INACCURATE loop",
        "Fold the cell's drift into its caller pass",
    )
    check(
        "drift census gated on NEW_FREQ or refresh calls",
        block.count("IF(NEW_FREQ .OR. MOM_KN_REFRESHING)") >= 6,
    )
    check(
        "quarantine age ticks once per sweep, not per call",
        "IF(NEW_FREQ .AND. .NOT. MOM_KN_REFRESHING)THEN" in block,
    )
    check(
        "capped/boundary counters share the census gate",
        "MOM_INEXT_CAPPED=MOM_INEXT_CAPPED+N_QBND"
        in block.split("END IF")[0],
    )
    check(
        "significant-negative-J counter gated the same way",
        re.search(
            r"IF\(NEW_FREQ \.OR\. MOM_KN_REFRESHING\)\s*\n"
            r"\s*1\s+MOM_NEGJ_LRG=MOM_NEGJ_LRG\+1",
            MOM,
        )
        is not None,
    )
    check(
        "pass-1/pass-2 fold kept outside the census gate",
        "MOM_GS_CHG_P2=MAX(MOM_GS_CHG_P2,T5)" in MOM
        and ".NOT. MOM_KN_REFRESHING" in section(
            MOM,
            "Fold the cell's drift into its caller pass",
            "END DO",
        ),
    )


def test_mom_sweep_report_reads_sweep_drift():
    """MOM_KNOP_REFRESH_V1 overwrites MOM_GS_MAXCHG with the refresh
    residual (BEST_CHG) for the adaptive relaxation; the sweep table
    must print MOM_SW_MAXCHG alongside the restored sweep counts or
    'drift' and ITER_SUM describe different quantities.
    """
    body = section(MOM, "SUBROUTINE DUMP_MOM_JNEXT_GS", "Per-pass split")
    check(
        "sweep drift falls back to MOM_SW_MAXCHG after a refresh",
        "IF(MOM_RF_RAN)SW_DRIFT=MOM_SW_MAXCHG" in body,
    )
    check(
        "drift line prints SW_DRIFT with the sweep counts",
        "'lagged J_NEXT drift ',SW_DRIFT" in body,
    )
    check(
        "argmax history and status verdicts use SW_DRIFT",
        "MOM_ARG_HV(4)=SW_DRIFT" in body
        and body.count("SW_DRIFT .GT.") >= 2,
    )
    check(
        "refresh restores the sweep argmax with the sweep counts",
        "ARGI_SV=MOM_GS_ARGI" in MOM
        and "MOM_GS_ARGI=ARGI_SV" in MOM
        and "MOM_GS_ARGL=ARGL_SV" in MOM,
    )


def test_mom_stall_text_matches_direction_schedule():
    check(
        "compressing-stall text conditional on the schedule",
        re.search(
            r"IF\(MOM_DOM_CMP\)THEN[\s\S]{0,400}?"
            r"' stalled - dropping it after pass',IPASS,[\s\S]{0,200}?"
            r"no other direction is scheduled",
            MOM,
        )
        is not None,
    )
    check(
        "expanding-stall text conditional on the schedule",
        re.search(
            r"IF\(MOM_DOM_EXP\)THEN[\s\S]{0,400}?"
            r"no other direction is scheduled",
            MOM,
        )
        is not None,
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


def test_knop_blue_boundary_compressing_coupling():
    """L=1 (bluest frequency) must keep the compressing upwind coupling.

    Regression pin for dc41573a: at the first frequency expanding cells
    have no bluer upwind (Q=0 is the correct inflow boundary), but
    compressing cells couple to the redder lagged bin, so the sweep
    INIT branch and the refresh BUILD_FREQ L=1 branch must both apply
    the chi_hat/Q treatment to GAM<0 cells.  Forcing Q=0 there dropped
    the advective outflow and left the refresh and sweep disagreeing
    at the band edge.
    """
    build_freq = section(FG, "SUBROUTINE BUILD_FREQ", "SUBROUTINE MATCH_FREQ")
    l1 = section(build_freq, "IF(L .EQ. 1)THEN", "ELSE")
    check(
        "refresh L=1 calls MATCH_FREQ for the matched span",
        "CALL MATCH_FREQ(L)" in l1,
    )
    check(
        "refresh L=1 adds |GAM|/DLN to chi_hat for GAM<0",
        "CHI_L(I)=CHI_L(I)+GAM_MAG(I,LS)/DLN" in l1,
    )
    check(
        "refresh L=1 sets signed Q=GAM/(DLN*chi_hat) for GAM<0",
        "Q_L(I)=GAM(I,LS)/(DLN*CHI_L(I))" in l1,
    )
    check(
        "refresh L=1 keeps Q_L=0 default for expanding cells",
        re.search(r"Q_L\(I\)=0\.0_LDP[\s\S]*?IF\(GAM\(I,LS\) \.LT\.", l1)
        is not None,
    )
    src_prime = section(FG, "SUBROUTINE SRC_PRIME", "SUBROUTINE SLOPES")
    check(
        "compressing source uses ABS(Q)*redder upwind",
        "SP_L(I)=SR_L(I)+ABS(Q_L(I))*INEXT_V" in src_prime,
    )
    check(
        "compressing upwind reads the matched redder slot L+NFM",
        "FLD(JPR1(I),L+NFM_L(I))" in src_prime,
    )
    sweep_init = section(
        FG,
        "Compressing cells at the bluest frequency",
        "SOURCE_RAY(1:NI)=ETA_RAY(1:NI)/CHI_RAY(1:NI)",
    )
    check(
        "sweep INIT adds |GAM|/DLN to chi_hat for GAM<0",
        "CHI_RAY(I)=CHI_RAY(I)+T1" in sweep_init
        and "T1=GAM_MAG(I,LS)/DLN_C" in sweep_init,
    )
    check(
        "sweep INIT sets signed Q=GAM/(DLN*chi_hat) for GAM<0",
        "Q(I)=GAM(I,LS)/(DLN_C*CHI_RAY(I))" in sweep_init,
    )


def test_ord_sweep_dispatch_fallback():
    """Both field directions must dispatch SWEEP_FLD when ordered is off.

    Pillar-1 pin: ORD_SWEEP requires FG_ORD_SWEEP, a live coefficient
    pack, NI>=3, and a successful SPW allocation.  In every disabling
    configuration the fallback for BOTH the inward (I_M) and outward
    (I_P) solves is the plain alternating SWEEP_FLD -- never a blind
    SWEEP_FLD_ORD call that would read an unallocated pack or SPW.
    """
    ray = section(FG, "SUBROUTINE KNOP_RAY_SOLVE_V1", "CONTAINS")
    check(
        "ORD_SWEEP needs FG_ORD_SWEEP + pack + NI>=3",
        "ORD_SWEEP=FG_ORD_SWEEP .AND. HAVE_PACK .AND. NI .GE. 3" in ray,
    )
    check(
        "SPW allocation failure drops ordered sweep",
        re.search(
            r"ALLOCATE\(SPW\(NCF,4\),STAT=IOS\)[\s\S]*?"
            r"IF\(IOS \.NE\. 0\)ORD_SWEEP=\.FALSE\.",
            ray,
        )
        is not None,
    )
    disp = [
        d
        for d in re.findall(
            r"IF\(ORD_SWEEP\)THEN[\s\S]{0,200}?ELSE[\s\S]{0,200}?END IF",
            ray,
        )
        if "CALL SWEEP_FLD" in d
    ]
    # Four sweep-dispatch sites: production I_M and I_P solves plus the
    # two KN_Q_RESID diagnostic verify passes (M and P channels), which
    # must mirror the same effective branch so the measured residual is
    # the production map. The per-ray census IF/ELSE on ORD_SWEEP is not
    # a sweep dispatch and is filtered out above.
    check("four ordered/fallback dispatch sites", len(disp) == 4)
    for d in disp:
        check(
            "ordered branch calls SWEEP_FLD_ORD",
            "CALL SWEEP_FLD_ORD" in d,
        )
        check(
            "fallback branch calls SWEEP_FLD (not SWEEP_FLD_ORD)",
            re.search(r"ELSE\s+CALL SWEEP_FLD\(", d) is not None,
        )
    sweep = section(FG, "SUBROUTINE SWEEP_FLD(FLD", "SUBROUTINE KN_ORD_SPV")
    pack_refs = [
        ln.strip()
        for ln in sweep.splitlines()
        if re.search(r"\bPK_", ln)
    ]
    guarded = section(sweep, "IF(HAVE_PACK)THEN", "END IF")
    check(
        "SWEEP_FLD reads the pack only inside the HAVE_PACK guard",
        len(pack_refs) > 0 and all(p in guarded for p in pack_refs),
    )
    check(
        "fallback path never touches SPW",
        "SPW" not in sweep,
    )
    march = section(FG, "SUBROUTINE MARCH(FLD", "SUBROUTINE WRITE_STORE")
    check(
        "MARCH never touches the pack or SPW",
        "PK_" not in march and "SPW" not in march,
    )


def test_missing_donor_consistent_drop():
    """A missing lagged donor drops BOTH sides of the coupling.

    Pillar-2 pin: the previous code kept the +|GAM|/DLN outflow in
    chi_hat while the source dropped the upwind inflow, which wrote a
    sink the discretisation never contained.  The sweep gate
    KN_COUPLING_OK and the refresh gate KN_RF_COUPLING must now hold
    the chi_hat increment, its DCHR derivative, and the signed-Q term
    out of the operator for donor-less compressing cells, and the
    no-slot source path must fall back to SOURCE_RAY alone.
    """
    check(
        "sweep coupling gate exists",
        "LOGICAL FUNCTION KN_COUPLING_OK(I2,LS,NM2,NCF2,ND2)" in FG,
    )
    chi_hat = section(
        FG, "Generalized opacity: chi_hat", "Guard: ensure chi_hat"
    )
    check(
        "sweep chi_hat term gated on the donor",
        "KN_COUPLING_OK(I,LS,NM_RAY(I),NCF_IN,ND)" in chi_hat,
    )
    dchr = section(
        FG,
        "Update dCHIdR_RAY to include the derivative",
        "Dimensionless signed coupling",
    )
    check(
        "sweep dCHIdR term gated on the donor",
        "KN_COUPLING_OK(I,LS,NM_RAY(I),NCF_IN,ND)" in dchr,
    )
    signed_q = section(
        FG,
        "Dimensionless signed coupling: Q =",
        "SOURCE_RAY(1:NI)=ETA_RAY(1:NI)/CHI_RAY(1:NI)",
    )
    check(
        "sweep signed-Q gated and zeroed for dropped cells",
        "KN_COUPLING_OK(I,LS,NM_RAY(I),NCF_IN,ND)" in signed_q
        and "Q(I)=0.0_LDP" in signed_q,
    )
    build_freq = section(FG, "SUBROUTINE BUILD_FREQ", "SUBROUTINE MATCH_FREQ")
    check(
        "refresh coupling gate exists",
        "LOGICAL FUNCTION KN_RF_COUPLING(I2,L2)" in FG,
    )
    check(
        "refresh chi_hat/DCHR/Q gated on the donor",
        build_freq.count("KN_RF_COUPLING(I,L)") >= 5,
    )
    check(
        "refresh gate keeps boundary outflow (foot off-grid)",
        re.search(
            r"IF\(L2\+NFM_L\(I2\) \.GT\. NCF\)RETURN", FG
        )
        is not None,
    )
    check(
        "sweep no-slot source path still drops to SOURCE_RAY",
        FG.count("SOURCE_PRIME(I)=SOURCE_RAY(I)") >= 2,
    )


def test_build_freq_fail_closed():
    """Non-finite coefficient builds must fail the ray, not be repaired.

    Pillar-3 pin: BUILD_FREQ's L>1 branch previously used NaN-only
    tests so +Inf slipped into DTAU and the march coefficients and was
    committed as a plausible field.  Every value guard now uses the
    finite test, and any non-finite input flags DIVERGED so the ray
    takes the diverged-ray exit (RESID=1, no store commit) rather than
    polishing a repaired operator into a converged refresh.
    """
    build_freq = section(FG, "SUBROUTINE BUILD_FREQ", "SUBROUTINE MATCH_FREQ")
    check(
        "no NaN-only guards remain in BUILD_FREQ",
        ".NE. CHI_L(I)" not in build_freq
        and ".NE. ETA_L(I)" not in build_freq
        and ".NE. Q_L(I)" not in build_freq
        and ".NE. SR_L(I)" not in build_freq
        and ".NE. DTAU_L(I)" not in build_freq
        and ".NE. DCHR_L(I)" not in build_freq,
    )
    check(
        "non-finite values flag the ray as diverged",
        build_freq.count("DIVERGED=.TRUE.") >= 8,
    )
    check(
        "matched-span floor flags the failed build",
        build_freq.count("DIVERGED=.TRUE.\n\t      END IF") >= 3
        or build_freq.count("DIVERGED=.TRUE.") >= 8,
    )
    check(
        "march coefficients get a final finite check",
        re.search(
            r"ABS\(A0L\(I\)\) \.LE\. HUGE\(A0L\(I\)\)[\s\S]*?DIVERGED",
            build_freq,
        )
        is not None,
    )
    ray = section(FG, "SUBROUTINE KNOP_RAY_SOLVE_V1", "CONTAINS")
    check(
        "pack build aborts on a flagged ray",
        re.search(
            r"CALL BUILD_FREQ\(L\)[\s\S]{0,400}?IF\(DIVERGED\)EXIT",
            ray,
        )
        is not None,
    )
    check(
        "flagged pack build takes the diverged-ray exit",
        re.search(
            r"IF\(DIVERGED\)THEN[\s\S]{0,100}?NDIV=1[\s\S]{0,100}?GOTO 800",
            ray,
        )
        is not None,
    )
    sweep = section(FG, "SUBROUTINE SWEEP_FLD(FLD", "SUBROUTINE KN_ORD_SPV")
    check(
        "no-pack sweep aborts a flagged build",
        re.search(
            r"CALL BUILD_FREQ\(L\)[\s\S]{0,400}?IF\(DIVERGED\)EXIT",
            sweep,
        )
        is not None,
    )
    check(
        "flagged build exits the pass loop as diverged",
        re.search(
            r"IF\(DIVERGED\)THEN[\s\S]{0,60}?IXIT=5[\s\S]{0,60}?EXIT",
            sweep,
        )
        is not None,
    )


def test_mom_gam_fac_floor_discontinuity_pin():
    """Pin the current MOM_GAM_FAC_MIN floor location (audit Phase 2).

    Documents where the sign-preserving 1+sigma*K floor is applied so a
    future stabilisation separation updates this test.  The floor is
    load-bearing for fixed-point convergence (PSI collapse at the
    zero-GAM resonance); it is NOT yet split into physical vs
    regularised coefficients.
    """
    check(
        "MOM_GAM_FAC_MIN parameter present",
        "MOM_GAM_FAC_MIN=1.0E-2_LDP" in MOM,
    )
    check(
        "depth-point factor floored with sign preservation",
        MOM.count("SIGN(MAX(ABS(T1),MOM_GAM_FAC_MIN),T1)") >= 2,
    )
    check(
        "exact-zero factor reset to zero",
        MOM.count("IF(T1 .EQ. 0.0_LDP)") >= 2,
    )


def test_qcell_histogram_gated_by_verify_pass():
    """KN_QC cell statistics must sample verify passes only.

    SWEEP_FLD_ORD computes its per-cell residuals inline and never
    calls MARCH, so each residual site carries its own gated
    KN_Q_ACCUM call. Accumulation is gated by the local L_QP
    verify-pass argument, never by shared module state, so a
    production pass can never enter the histogram.
    """
    ord_body = section(FG, "SUBROUTINE SWEEP_FLD_ORD", "SUBROUTINE BUILD_FREQ")
    sites = [
        m.start()
        for m in re.finditer(r"IF\(L_QP\)CALL KN_Q_ACCUM\(T1,INWARD\)", ord_body)
    ]
    check("SWEEP_FLD_ORD has four L_QP-gated KN_Q_ACCUM sites", len(sites) == 4)
    for pos in sites:
        prefix = ord_body[max(0, pos - 400) : pos]
        check(
            "each site follows a RES_FREQ update on the same T1",
            prefix.rstrip().endswith("RES_FREQ(L)=MAX(RES_FREQ(L),T1)"),
        )
    aa_block = section(ord_body, "AA_RH(J,L,2)", "END IF")
    check(
        "Anderson census block does not feed the cell histogram",
        "KN_Q_ACCUM" not in aa_block,
    )
    march_body = section(FG, "SUBROUTINE MARCH", "SUBROUTINE KN_Q_ACCUM")
    check(
        "MARCH keeps its four L_QP-gated KN_Q_ACCUM sites",
        len(re.findall(r"IF\(L_QP\)CALL KN_Q_ACCUM\(DCELL,INWARD\)", march_body))
        == 4,
    )
    check(
        "MARCH signature takes the verify-pass flag",
        "SUBROUTINE MARCH(FLD,L,INWARD,RES,L_QP)" in FG,
    )
    check(
        "SWEEP_FLD forwards L_QP into MARCH",
        "CALL MARCH(FLD,L,INWARD,RES,L_QP)" in FG,
    )


def test_qcell_histogram_race_free_layout():
    """Phase-4 race fix: per-ray/per-channel accumulation.

    The shared gate KN_Q_ACC/KN_Q_ACC_CH let one OpenMP ray arm
    accumulation while another was still in its production sweep,
    and the shared KN_QC_* counters were updated concurrently
    with no ordering, so the old histogram, percentiles and
    channel maxima are invalid evidence. Accumulation must now
    be per-ray: each worker writes only its own LS slice, the
    channel comes from the local INWARD flag, and the serial
    totals are reduced after the ray loop joins.
    """
    accum = section(FG, "SUBROUTINE KN_Q_ACCUM", "RETURN\n\tEND")
    check(
        "accumulator takes an explicit channel argument",
        "SUBROUTINE KN_Q_ACCUM(DX,L_IN)" in FG,
    )
    check(
        "channel resolved from the local inward flag",
        "IF(.NOT.L_IN)CH=2" in accum,
    )
    for pat in (
        "KN_QC_N_RAY(LS,CH)",
        "KN_QC_MAX_RAY(LS,CH)",
        "KN_QC_NT_RAY(LS,CH,",
        "KN_QC_HIST_RAY(LS,CH,",
    ):
        check(f"accumulator writes per-ray array {pat}", pat in accum)
    for pat in ("KN_QC_N(", "KN_QC_MAX(", "KN_QC_NT(", "KN_QC_HIST("):
        check(f"accumulator never touches shared {pat}", pat not in accum)
    check(
        "no shared accumulation gate assignments remain",
        "KN_Q_ACC=" not in FG and "KN_Q_ACC_CH=" not in FG,
    )
    check(
        "verify block carries no arming comments or calls",
        "KN_Q_ACC" not in section(
            FG, "KN_QC accumulation needs no arming", "Same argmax capture"
        ),
    )
    check(
        "per-ray storage allocated for rays x channels x bins",
        "KN_QC_HIST_RAY(KN_NP,2,403)" in FG
        and "KN_QC_NT_RAY(KN_NP,2,3)" in FG
        and "KN_QC_MAX_RAY(KN_NP,2)" in FG,
    )
    check(
        "per-ray storage deallocated with the other diagnostics",
        "DEALLOCATE(KN_QC_HIST_RAY)" in FG
        and "DEALLOCATE(KN_QC_N_RAY)" in FG
        and "DEALLOCATE(KN_QC_NT_RAY)" in FG
        and "DEALLOCATE(KN_QC_MAX_RAY)" in FG,
    )
    red = section(FG, "Reduce the per-ray Q cell statistics", "T1_W=OMP_GET_WTIME")
    check(
        "serial reduction sums counts and histograms over rays",
        "KN_QC_N=SUM(KN_QC_N_RAY,DIM=1)" in red
        and "KN_QC_NT=SUM(KN_QC_NT_RAY,DIM=1)" in red
        and "KN_QC_HIST=SUM(KN_QC_HIST_RAY,DIM=1)" in red,
    )
    check(
        "serial reduction takes per-channel maxima over rays",
        "KN_QC_MAX(J)=MAXVAL(KN_QC_MAX_RAY(:,J))" in red,
    )
    check(
        "reduction runs after the parallel ray loop",
        FG.index("!$OMP END PARALLEL DO")
        < FG.index("Reduce the per-ray Q cell statistics"),
    )


def test_qargmax_store_mapping_diagnostics():
    """Phase-4 Q argmax write-eligibility audit diagnostics.

    The report must invert J_PNT for the Q argmax ray point,
    list every matching model depth with its compact store row,
    classify writability, read back the retained store value and
    print VMAX_RAW next to the floor-protected VMAXF.
    """
    check(
        "KN_QA_VM stores the raw column max",
        "KN_QA_VM(LS,1)=MAXVAL(FLD(1:NI,WRC_L))" in FG
        and "KN_QA_VM(LS,2)=MAXVAL(FLD_P(1:NI,WRC_L))" in FG,
    )
    check(
        "report prints raw and floored peaks",
        "vold,vnew,vmaxraw,vmaxf,denom=" in FG,
    )
    check(
        "inverse J_PNT scan lists every matching model depth",
        "J_PNT(IQ_KM,IQ_R) .EQ. KN_QA_J(IQ_R,IQ_C)" in FG
        and "qmap2 K/jpn/cdi/ok=" in FG,
    )
    check(
        "writability classes A/B/C are reported",
        "qmap2 class=A" in FG and "qmap2 class=B" in FG
        and "qmap2 class=C" in FG,
    )
    check(
        "retained store value read at slot L-1",
        "IQ_SLOT=KN_QA_L(IQ_R,IQ_C)-1" in FG
        and "I_P_DEPTH_NEXT(IQ_KC2,IQ_R,IQ_SLOT)" in FG
        and "I_M_DEPTH_NEXT(IQ_KC2,IQ_R,IQ_SLOT)" in FG,
    )
    check(
        "pre-refresh store snapshot armed from last argmax",
        "KN_QARG_VB=I_P_DEPTH_NEXT(KN_QARG_KC,KN_QARG_LS," in FG
        and "KN_QARG_L-1" in FG,
    )


def test_qcell_histogram_alloc_failure_path():
    """KN_QC_READY must gate accumulation on allocation success.

    The per-ray KN_QC_*_RAY arrays are allocated with STAT=IOS2.
    A failed allocation used to leave KN_Q_ACCUM writing to
    unallocated storage on the next L_QP verify pass, and the
    zero-initialised globals could read as a measured empty
    histogram. KN_QC_READY now tracks allocation end to end:
    FALSE before the attempt, TRUE only when the STAT is clean
    and every array the accumulator and reduction touch is
    allocated, FALSE again at cleanup. L_QP keeps its verify-pass
    meaning and is never reused as the readiness flag.
    """
    check(
        "KN_QC_READY declared and defaults FALSE",
        "LOGICAL, SAVE :: KN_QC_READY=.FALSE." in FG,
    )
    alloc = section(
        FG,
        "KN_QC_READY=.FALSE.\n\t  ALLOCATE(KN_QS_STAT",
        "RES_LS=0.0_LDP",
    )
    check(
        "readiness reset before the cell-stat allocation",
        "KN_QC_READY=.FALSE." in alloc,
    )
    check(
        "readiness requires clean STAT and every KN_QC array",
        "IOS2 .EQ. 0 .AND. ALLOCATED(KN_QC_HIST)" in alloc
        and "ALLOCATED(KN_QC_HIST_RAY)" in alloc
        and "ALLOCATED(KN_QC_N_RAY)" in alloc
        and "ALLOCATED(KN_QC_NT_RAY)" in alloc
        and "ALLOCATED(KN_QC_MAX_RAY)" in alloc
        and "KN_QC_READY=.TRUE." in alloc,
    )
    check(
        "failure branch zeroes scalars and reports unavailable",
        "KN_QC cell histogram unavailable" in alloc
        and "Q residual unaffected" in alloc
        and "KN_QC_N=0_8; KN_QC_NT=0_8; KN_QC_MAX=0.0_LDP" in alloc,
    )
    accum = section(FG, "SUBROUTINE KN_Q_ACCUM", "RETURN\n\tEND")
    check(
        "accumulator returns before touching storage when not ready",
        accum.index("IF(.NOT.KN_QC_READY)RETURN")
        < accum.index("KN_QC_N_RAY(LS,CH)"),
    )
    red = section(FG, "Reduce the per-ray Q cell statistics", "T1_W=OMP_GET_WTIME")
    check(
        "serial reduction gated on KN_QC_READY",
        "IF(KN_QC_READY)THEN" in red
        and "ALLOCATED(KN_QC_HIST_RAY)" not in red,
    )
    clean = section(FG, "DEALLOCATE(KN_QS_STAT)", "CONTAINS")
    check(
        "cleanup deallocates each KN_QC array separately",
        "IF(ALLOCATED(KN_QC_HIST))DEALLOCATE(KN_QC_HIST)" in clean
        and "IF(ALLOCATED(KN_QC_HIST_RAY))DEALLOCATE(KN_QC_HIST_RAY)"
        in clean
        and "IF(ALLOCATED(KN_QC_MAX_RAY))DEALLOCATE(KN_QC_MAX_RAY)"
        in clean,
    )
    check(
        "cleanup resets readiness",
        "KN_QC_READY=.FALSE." in clean,
    )
    check(
        "summary distinguishes unavailable from an empty histogram",
        "q cells unavailable: histogram" in FG
        and "IF(KN_QC_READY)THEN\n\t      DO J=1,2" in FG,
    )
    check(
        "L_QP is not reused as the readiness flag",
        "KN_QC_READY" not in section(FG, "SUBROUTINE MARCH", "END"),
    )


def test_jump_diag_veto_attribution():
    fn = section(FG, "FUNCTION FG_JUMP_MAX", "SUBROUTINE HEAL_STORE")
    # A. ordinary jump: argmax cell tracked with strict-> tie rule
    check(
        "ordinary jump records argmax cell",
        "VR=ABS(V1-V2)/" in fn
        and "MAX(V2,VMAX_S(LS2)/FG_COMMIT_JUMP)" in fn
        and "VR .GT. JT_BEST" in fn
        and "%BR=3" in fn
        and "%KS=JT_KS" in fn
        and "%L2=JT_L2" in fn,
    )
    # B. sibling-reference fallback branch 2 with ratio and threshold
    check(
        "sibling-reference branch records ratio and threshold",
        "%BR=2" in fn
        and "%RATIO=V1/REF(LS2)" in fn
        and "%THR=FG_COMMIT_JUMP*REF(LS2)" in fn,
    )
    # C. non-finite candidate branch 1 before the HUGE return
    check(
        "non-finite candidate branch 1 recorded",
        "%BR=1" in fn and "%NF=.TRUE." in fn,
    )
    br1 = section(fn, "%BR=1", "FG_JUMP_MAX=HUGE")
    check(
        "branch 1 record precedes the HUGE return",
        "%JVAL=HUGE(FG_JUMP_MAX)" in br1,
    )
    # D. no trigger reported below thresholds: report filter
    rep = section(
        FG,
        "IF(KN_JUMP_DIAG .AND. KN_JT_READY)THEN",
        "END IF\n\t  ELSE IF",
    )
    check(
        "report fires only for exceeding channels",
        "KN_JT(LS,J)%BR .LE. 2 .OR." in rep
        and "KN_JT(LS,J)%JVAL .GT. FG_COMMIT_JUMP" in rep,
    )
    # E. deterministic ties: strict > keeps first cell in scan order
    check(
        "tie-break is strict greater-than (first cell wins)",
        ".GT. JT_BEST" in fn and ".GE. JT_BEST" not in fn,
    )
    # F. diagnostic disabled by default; gated writes
    check(
        "diagnostic gated by KN_JUMP_DIAG and KN_JT_READY",
        "L_DIAG=KN_JUMP_DIAG .AND. KN_JT_READY" in fn
        and "DATA KN_JUMP_DIAG /.FALSE./" in FG,
    )
    check(
        "every record write is inside an L_DIAG guard",
        fn.count("IF(L_DIAG") >= 4,
    )
    # G. OpenMP safety: per-ray/per-channel slices, critical print
    check(
        "records indexed by ray and channel only",
        "KN_JT(LS,ICH)%" in fn
        and "KN_JT(KN_NP,2)" in FG,
    )
    check(
        "report serialized under OMP critical",
        "!$OMP CRITICAL (KN_JT_UPD)" in FG
        and "!$OMP END CRITICAL (KN_JT_UPD)" in FG,
    )
    # production behaviour preserved
    check(
        "production return values unchanged",
        "FG_JUMP_MAX=MAX(FG_JUMP_MAX,VR)" in fn
        and "FG_JUMP_MAX=HUGE(FG_JUMP_MAX)" in fn,
    )
    check(
        "both call sites pass channel indices",
        "FG_JUMP_MAX(FLD,I_M_DEPTH_NEXT,JMP_SIB_M,1)" in FG
        and "FG_JUMP_MAX(FLD_P,I_P_DEPTH_NEXT,JMP_SIB_P,2)" in FG,
    )
    check(
        "VADAT plumbing wired",
        "LOGICAL KN_JUMP_DIAG" in CTRL_MOD
        and "'KN_JUMP_DIAG'" in RD_CONTROL
        and "CALL SET_FG_JUMP_DIAG(KN_JUMP_DIAG)" in RD_CONTROL
        and "SUBROUTINE SET_FG_JUMP_DIAG" in FG,
    )
    check(
        "cleanup deallocates KN_JT and resets readiness",
        "IF(ALLOCATED(KN_JT))DEALLOCATE(KN_JT)" in FG
        and "KN_JT_READY=.FALSE." in FG,
    )


def test_kn_sib_depth_snapshot():
    """Same-depth sibling snapshot for the ray-83 veto cell.

    The snapshot must run serially before the OpenMP ray loop (the
    store is about to be rewritten by committing rays), must resolve
    J_PNT(K,RS) per ray with the NI_RAY(RS) validity bound, must
    exclude non-finite and non-positive retained values, must cover
    both channels, and must document the REF distinction.
    """
    snap = section(
        FG,
        "SUBROUTINE KN_SIB_DEPTH_SNAP_V1",
        "END SUBROUTINE KN_SIB_DEPTH_SNAP_V1",
    )
    check("snapshot subroutine exists", len(snap) > 0)
    check(
        "per-ray J_PNT resolution with NI_RAY bound",
        "JP=J_PNT(KQ,RS)" in snap
        and "JP .LT. 1 .OR. JP .GT. NI_RAY(RS)" in snap,
    )
    check(
        "compact row from COMP_DEPTH_IDX and slot offset",
        "KC=COMP_DEPTH_IDX(KQ)" in snap
        and "STORE(KC,RS,NLO+IS-1)" in snap,
    )
    check(
        "non-finite and non-positive values excluded",
        ".NOT.(ABS(V) .LE. HUGE(V))" in snap and "V .LE. 0.0_LDP" in snap,
    )
    check(
        "documented percentile convention present",
        "CEILING(0.10_LDP*NOK)" in snap
        and "CEILING(0.90_LDP*NOK)" in snap
        and "(NOK+1)/2" in snap,
    )
    check(
        "deterministic insertion sort",
        "WRK(J2+1)=VT" in snap and "DO WHILE(J2 .GE. 1)" in snap,
    )
    check(
        "ray rank counts strictly-lower siblings",
        "IF(WRK(I2) .LT. V)J2=J2+1" in snap
        and "ST(IS)%PCT=100.0_LDP*J2/NOK" in snap,
    )
    check(
        "allocation failure marks N=-1 not a fake zero",
        "ST(IS)%N=-1" in snap,
    )

    # Call site: inside FG_KNOP_REFRESH_V1, after the sibling
    # reference build and before the OMP parallel loop, gated on
    # KN_JUMP_DIAG with allocation guards.
    call = re.search(
        r"IF\(KN_JUMP_DIAG \.AND\. ALLOCATED\(J_PNT\).*?"
        r"CALL KN_SIB_DEPTH_SNAP_V1\(2,1434,1438,83,"
        r".*?I_M_DEPTH_NEXT,KN_SIB_M\)"
        r".*?CALL KN_SIB_DEPTH_SNAP_V1\(2,1434,1438,83,"
        r".*?I_P_DEPTH_NEXT,KN_SIB_P\)"
        r".*?END IF",
        FG,
        re.S,
    )
    check("gated serial snapshot call for both channels", call is not None)
    if call is not None:
        pos = call.start()
        par = FG.find("!$OMP PARALLEL DO", pos)
        check(
            "snapshot precedes the OpenMP ray loop",
            par > pos,
        )
        sibref = FG.rfind("KN_JMP_SIBREF_V1(I_P_DEPTH_NEXT", 0, pos)
        check(
            "snapshot runs after sibling references are built",
            sibref > 0 and sibref < pos,
        )

    # No cross-ray store reads inside FG_JUMP_MAX or the ray loop:
    # the snapshot is the only new consumer of sibling columns.
    jmp = section(FG, "FUNCTION FG_JUMP_MAX", "END")
    check(
        "FG_JUMP_MAX does not read sibling rays",
        "KN_SIB" not in jmp,
    )
    check(
        "report block gated on KN_JUMP_DIAG",
        "IF(KN_JUMP_DIAG)THEN" in FG
        and "KN-SIB K=" in FG
        and "ref-colmax-med=" in FG,
    )
    check(
        "REF labelled as column-max median not same-depth",
        "ref-colmax-med=" in FG
        and "column maxima" in FG,
    )
    check(
        "trigger V1 ratios reported against same-depth median",
        "v1,medP,v1/medP,v1/v2=" in FG,
    )


def test_kn_census_and_geo_diagnostics():
    """Full-store census and geometry-conditioned sibling stats.

    The census hooks must live inside FG_JUMP_MAX's existing cell
    loop (the only place V1/V2/VR are all live), stay indexed by
    (LS,ICH) for OpenMP safety, stay bounded, and stay gated on
    KN_JUMP_DIAG via KN_CEN_READY. The geo snapshot must recover
    impact parameter from R_RAY^2-Z^2 (the P dummy is not
    retained), use MU=Z/R_RAY, and report three subsets.
    """
    jmp = section(FG, "FUNCTION FG_JUMP_MAX", "END\n!")
    check("census cell counter inside the cell loop", "KN_CEN_NCELL(LS,ICH)" in jmp)
    check(
        "exceed counter and histogram gated on VR>FG_COMMIT_JUMP",
        "IF(VR .GT. FG_COMMIT_JUMP)THEN" in jmp
        and "KN_CEN_NEXC(LS,ICH)" in jmp
        and "KN_CEN_HIST(LS,ICH,KS,KN_CEN_IB)" in jmp,
    )
    check(
        "top-N records bounded and deterministic",
        "KN_CEN_TOP(LS,ICH,KN_CEN_NTOP)%VR" in jmp
        and "KN_CEN_IT=KN_CEN_NTOP" in jmp
        and "KN_CEN_TOP(LS,ICH,KN_CEN_IT)=\n\t1             KN_CEN_TOP(LS,ICH,KN_CEN_IT-1)" in jmp.replace(" ", " "),
    )
    check(
        "per-ray JMAX recorded at all three exits",
        jmp.count("KN_CEN_JMAX(LS,ICH)=") >= 3,
    )
    check(
        "V2-invalid and floor-column cells counted separately",
        "KN_CEN_NSKV(LS,ICH)" in jmp and "KN_CEN_NFLR(LS,ICH)" in jmp,
    )
    check(
        "no full N_COMP x NP x NCF census array",
        "KN_CEN_HIST(KN_NP,2,KN_ND,KN_CEN_NB)" in FG
        and "KN_CEN_TOP(KN_NP,2,KN_CEN_NTOP)" in FG,
    )
    check(
        "census gated on KN_CEN_READY not just KN_JUMP_DIAG",
        "L_CEN=KN_JUMP_DIAG .AND. KN_CEN_READY" in jmp,
    )
    check(
        "readiness requires all arrays allocated",
        "KN_CEN_READY=ALLOCATED(KN_CEN_NCELL)" in FG
        and "ALLOCATED(KN_CEN_TOP)" in FG,
    )
    check(
        "cleanup deallocates census and geo storage",
        "IF(ALLOCATED(KN_CEN_TOP))DEALLOCATE(KN_CEN_TOP)" in FG
        and "IF(ALLOCATED(KN_GEO_REPS))DEALLOCATE(KN_GEO_REPS)" in FG
        and "KN_CEN_READY=.FALSE." in FG,
    )

    geo = section(FG, "SUBROUTINE KN_SIB_GEO_V1", "END SUBROUTINE KN_SIB_GEO_V1")
    check("geo subroutine exists", len(geo) > 0)
    check(
        "impact parameter recovered from ray grid",
        "PG(N1)=SQRT(MAX(R_RAY(JP,RS)**2-Z(JP,RS)**2,0.0_LDP))" in geo,
    )
    check(
        "direction cosine and signed drift captured",
        "PM(N1)=Z(JP,RS)/R_RAY(JP,RS)" in geo
        and "PA(N1)=GAM(JP,RS)" in geo,
    )
    check(
        "J_PNT validity bound applied per ray",
        "JP .LT. 1 .OR. JP .GT. NI_RAY(RS)" in geo,
    )
    check(
        "invalid retained values excluded",
        "V .LE. 0.0_LDP" in geo and ".NOT.(ABS(V) .LE. HUGE(V))" in geo,
    )
    check(
        "three subsets reported: all, near-P, near-MU",
        "REP%SUB(1)" in geo and "REP%SUB(2)" in geo and "REP%SUB(3)" in geo
        and "ABS(PG-REP%P_T)" in geo and "ABS(PM-REP%MU_T)" in geo,
    )
    check(
        "subset size bounded by KN_GEO_NB",
        "N2=MIN(KN_GEO_NB,N1)" in geo,
    )
    check(
        "same nearest-rank convention as sibling snap",
        "CEILING(0.10_LDP*NT)" in section(FG, "SUBROUTINE KN_GEO_STATS_V1", "END SUBROUTINE KN_GEO_STATS_V1"),
    )

    # Report blocks gated and bounded.
    check(
        "geo report gated on KN_JUMP_DIAG",
        "IF(KN_JUMP_DIAG .AND. ALLOCATED(KN_GEO_REPS))THEN" in FG
        and "KN-GEO ch" in FG,
    )
    check(
        "census report gated on KN_CEN_READY",
        "IF(KN_CEN_READY)THEN" in FG and "KN-CEN refresh=" in FG,
    )
    check(
        "top list capped at 16",
        "MIN(16,KN_CEN_NZ)" in FG,
    )
    check(
        "region map is depth x log-band not full array",
        "KN_CEN_H2(KN_ND,KN_CEN_NB)" in FG and "KN-CEN map K=" in FG,
    )


if __name__ == "__main__":
    test_kn_census_and_geo_diagnostics()
    test_kn_sib_depth_snapshot()
    test_jump_diag_veto_attribution()
    test_qcell_histogram_gated_by_verify_pass()
    test_qcell_histogram_race_free_layout()
    test_qargmax_store_mapping_diagnostics()
    test_qcell_histogram_alloc_failure_path()
    test_ord_sweep_dispatch_fallback()
    test_missing_donor_consistent_drop()
    test_build_freq_fail_closed()
    test_knop_blue_boundary_compressing_coupling()
    test_mom_gam_fac_floor_discontinuity_pin()
    test_heal_store_zeroes_fully_poisoned_rows()
    test_sweep_store_write_snapshot_relax()
    test_rollback_restore_validates_field_values()
    test_fail_closed_commit_defaults()
    test_commit_jump_veto_seed_waiver()
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
    test_mom_sweep_census_counts_first_call_only()
    test_mom_sweep_report_reads_sweep_drift()
    test_mom_stall_text_matches_direction_schedule()
    test_fg_trial_churn_arm()
    print("all KNOP lagged-field checks passed")
