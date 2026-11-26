#!/usr/bin/env python3
"""Check JFNK control-flow guards that need CMFGEN state to run directly."""

from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "cmfgen"
INTERFACE = (SRC / "new_main/mod_subs/jfnk_cmfgen_interface_v1.f").read_text()
SOLVER = (SRC / "new_main/mod_subs/jfnk_solver_v1.f90").read_text()
POPS = (SRC / "new_main/mod_subs/solve_for_pops.f").read_text()
SOLVE_CORR = (SRC / "new_main/mod_subs/jfnk_solve_corr_v1.f").read_text()
INTMOD = (SRC / "new_main/mod_subs/jfnk_cmfgen_int_mod.f").read_text()
CVARMOD = (SRC / "new_main/mod_subs/control_variable_mod.f").read_text()
BLKBAND = (SRC / "new_main/subs/cmf_blkband_v3.f").read_text()
RD_CONTROL = (SRC / "new_main/mod_subs/rd_control_variables.f").read_text()
CMFGEN_SUB = (SRC / "new_main/cmfgen_sub.f").read_text()


def section(source, start, end):
    first = source.index(start)
    last = source.index(end, first)
    return source[first:last]


def check(name, passed):
    if not passed:
        raise AssertionError(name)
    print(f"[PASS] {name}")


factor = section(
    INTERFACE,
    "SUBROUTINE JFNK_FACTOR_PRECOND_V1",
    "END SUBROUTINE JFNK_FACTOR_PRECOND_V1",
)
check(
    "cached BA diagonal is restored on reuse",
    "JFNK_BA_DIAG(J,K) = JFNK_BA_A(J,J,K)" in factor,
)
check(
    "BA diagonal reset is inside rebuild branch",
    re.search(
        r"IF\(REBUILD_BA\)THEN\s+JFNK_BA_DIAG\s*=\s*0\.0_LDP",
        factor,
        re.IGNORECASE,
    )
    is not None,
)
check("banded assembly does not retain local REPLACE", "REPLACE_SAVE" not in factor)
check(
    "banded assembly refreshes REPLACE per depth",
    re.search(
        r"GENERATE_FULL_MATRIX_V3\([\s\S]*?JFNK_DIAG_INDX,\s*"
        r"JFNK_DIAG_INDX,\s*K,[\s\S]*?DO IB = 1, JFNK_NUM_BNDS",
        factor,
        re.IGNORECASE,
    )
    is not None,
)
check(
    "PTC preconditioner fold tracks dtau and log coordinates",
    "MAX(jfnk_ptc_dtau, 1.0D-8)" in factor
    and "S_PTC = S_PTC * EXP(-X_PTC)" in factor,
)
check(
    "gate BA build is reused by the real factor",
    "(DO_FACTOR_L .AND. jfnk_ptc_active)" not in factor,
)

set_scale = section(
    INTERFACE,
    "SUBROUTINE JFNK_SET_PTC_SCALE_V1",
    "END SUBROUTINE JFNK_SET_PTC_SCALE_V1",
)
check(
    "PTC scale re-arm invalidates the folded BA cache",
    re.search(r"JFNK_BA_BUILT\s*=\s*\.FALSE\.", set_scale) is not None,
)
check(
    "BA invalidation precedes the allocation-failure returns",
    set_scale.index("JFNK_BA_BUILT = .FALSE.") < set_scale.index("RETURN"),
)

refresh = section(
    INTERFACE,
    "IF(jfnk_need_precond_refresh .NE. 0)THEN",
    "jfnk_need_precond_refresh = 0",
)
check(
    "precond refresh rebuilds the BA with live fold state",
    "JFNK_BA_BUILT = .FALSE." in refresh
    and refresh.index("JFNK_BA_BUILT = .FALSE.")
    < refresh.index("CALL JFNK_FACTOR_PRECOND_V1(.TRUE.)"),
)
post_guard = section(
    factor,
    "END IF   ! JFNK_BA_BUILT reuse guard",
    "IF(.NOT. DO_FACTOR_L)THEN",
)
check(
    "real factor marks the BA cache built",
    re.search(r"JFNK_BA_BUILT\s*=\s*\.TRUE\.", post_guard) is not None,
)

init = section(
    INTERFACE,
    "SUBROUTINE JFNK_INIT_V1",
    "END SUBROUTINE JFNK_INIT_V1",
)
check(
    "JFNK_INIT_V1 resets jfnk_ptc_dtau per activation",
    re.search(r"jfnk_ptc_dtau\s*=\s*1\.0D0", init, re.IGNORECASE)
    is not None,
)
check(
    "JFNK_INIT_V1 clears adaptive-JVP state per activation",
    "JFNK_HADAPT_MODE=-1" in init
    and "JFNK_HADAPT_REARM=.FALSE." in init
    and "JFNK_HADAPT_MULT=1.0_LDP" in init
    and "JFNK_LIN_STAGE=0" in init,
)

rebuild = section(
    INTERFACE,
    "SUBROUTINE JFNK_REBUILD_BASELINE_V1",
    "END SUBROUTINE JFNK_REBUILD_BASELINE_V1",
)
pre_factor = section(
    rebuild,
    "CALL JFNK_BUILD_BA_RAD_V1()",
    "CALL JFNK_FACTOR_PRECOND_V1(.TRUE.)",
)
check(
    "baseline rebuild re-arms PTC state before the factor call",
    re.search(r"jfnk_x\s*=\s*0\.0D0", pre_factor) is not None
    and "CALL JFNK_SET_PTC_SCALE_V1" in pre_factor
    and re.search(r"jfnk_ptc_dtau\s*=\s*1\.0D0", pre_factor) is not None,
)
check(
    "baseline rebuild masks fixed unknowns before PTC re-arm",
    pre_factor.index("CALL JFNK_MASK_FIXED_UNKNOWNS")
    < pre_factor.index("jfnk_x = 0.0D0")
    and pre_factor.index("CALL JFNK_MASK_FIXED_UNKNOWNS")
    < pre_factor.index("CALL JFNK_SET_PTC_SCALE_V1"),
)
retry = section(
    SOLVE_CORR,
    "REBUILDING PRECONDITIONER (attempt",
    "JFNK_EQ_STRUCT_READY=.FALSE.",
)
check(
    "attempt-retry rebuild clears the BA cache before factoring",
    "JFNK_BA_BUILT=.FALSE." in retry
    and retry.index("JFNK_BA_BUILT=.FALSE.")
    < retry.index("CALL JFNK_FACTOR_PRECOND_V1(.TRUE.)"),
)

gate_idx = SOLVE_CORR.index("CALL JFNK_FACTOR_PRECOND_V1(.FALSE.)")
check(
    "admission gate build follows F(0) adoption and PTC arming",
    gate_idx
    > SOLVE_CORR.index("CALL JFNK_SET_PTC_SCALE_V1(JFNK_N,jfnk_fx)")
    > SOLVE_CORR.rindex("CALL JFNK_MASK_FIXED_UNKNOWNS(JFNK_N,jfnk_fx)",
                      0, gate_idx)
    > SOLVE_CORR.rindex("jfnk_fx=F0_Y", 0, gate_idx),
)
check(
    "coupled rebuild forces a refactor before returning",
    "JFNK_FACTOR_DONE=.FALSE." in rebuild
    and "JFNK_LU_READY=.FALSE." in rebuild
    and "jfnk_need_precond_refresh=1" in rebuild,
)
check(
    "baseline rebuild clears adaptive-JVP state",
    "JFNK_HADAPT_MODE=-1" in rebuild
    and "JFNK_HADAPT_REARM=.FALSE." in rebuild
    and "JFNK_HADAPT_MULT=1.0_LDP" in rebuild
    and "JFNK_LIN_STAGE=0" in rebuild,
)

jvp = section(INTERFACE, "SUBROUTINE JFNK_JVP_V1", "END SUBROUTINE JFNK_JVP_V1")
bad_residual = section(
    jvp,
    "CALL JFNK_RESIDUAL_V1(N, JFNK_W_X_PERT, JFNK_W_F_PERT)",
    "JV(I) = (JFNK_W_F_PERT(I) - jfnk_fx(I)) / H_STEP",
)
check(
    "JVP rejects repaired perturbed residuals before differencing",
    "jfnk_resid_nbad .GT. 0 .OR. JFNK_EVAL_BROKEN" in bad_residual
    and "jfnk_jvp_failed = .TRUE." in bad_residual
    and "JV(I) = (JFNK_W_F_PERT(I) - jfnk_fx(I)) / H_STEP"
    not in bad_residual,
)
check(
    "HADAPT re-arms before applying a saved multiplier",
    jvp.index("JFNK_HADAPT_MODE .NE. MODE_HA")
    < jvp.index("IF(JFNK_HADAPT_MULT .NE. 1.0_LDP)")
    and "JFNK_JVP_DIAG_DONE=.FALSE." in jvp,
)
check(
    "HADAPT sweep bypasses the cheap-report direction cap",
    "JVP_DIR_CT .LT. JFNK_JVP_NDIR .OR. JFNK_HADAPT_REARM" in jvp
    and "JFNK_JVP_DIAG_DONE" in jvp
    and "JVP_DIAG_DONE" not in jvp.replace("JFNK_JVP_DIAG_DONE", ""),
)
check(
    "HADAPT state lives in the JFNK interface module",
    "JFNK_HADAPT_MODE=-1" in INTMOD
    and "JFNK_HADAPT_REARM=.FALSE." in INTMOD
    and "JFNK_JVP_DIAG_DONE=.FALSE." in INTMOD,
)

residual = section(
    INTERFACE,
    "SUBROUTINE JFNK_RESIDUAL_V1",
    "END SUBROUTINE JFNK_RESIDUAL_V1",
)
check(
    "temperature floor crossing is counted from the unclipped value",
    "T_BASE_UNCL .LT. 1.0_LDP" in residual
    and "POPS_BASE(JFNK_NT,K)*EXP(X_BASE) .GE. 1.0_LDP" in residual,
)
check(
    "JVP step cap includes the 1 K temperature floor",
    "X_FLOOR = LOG(1.0D0/POPS_BASE(JFNK_NT,K_J))" in jvp
    and "X_MARGIN = jfnk_x(I) - X_FLOOR" in jvp,
)

gmres = section(SOLVER, "subroutine gmres_m_v1", "end subroutine gmres_m_v1")
check(
    "GMRES restart recompute retains earlier cycles and bounds total step",
    "x = jfnk_gmres_x_cycle_start + z" in gmres
    and "if(x_cycle_norm > del) x = x * (del / x_cycle_norm)" in gmres,
)
newton = section(SOLVER, "subroutine newton_hook_v1", "end subroutine newton_hook_v1")
check(
    "coupled convergence does not use maxch alone",
    "(maxch_converged .and. .not. JFNK_COUPLED_LOC)" in newton,
)
check(
    "trust-region non-monotone trials are bounded and counted",
    "if(jfnk_del > 0d0) phi_ref = min(phi_ref, 1.21d0 * phi_0)" in newton
    and "else if(ls_accepted .and. del == 0d0)" in newton
    and "jfnk_n_tr_rejected = jfnk_n_tr_rejected + 1" in newton,
)
check(
    "GS_FREEZE is asserted outside the KNOP gate",
    "IF(.NOT. DO_NONMONO_KNOP .AND. GS_FREEZE)THEN" in POPS,
)
check(
    "GS_FREEZE is only assigned inside the KNOP gate",
    re.search(r"IF\(DO_NONMONO_KNOP\)THEN\s+GS_FREEZE\s*=", POPS)
    is not None,
)

for stub in (
    "tests/f90/jfnk_int_mod_stub.f90",
    "tests/jfnk_core/jfnk_int_stub.f90",
    "unit_tests/jfnk_int_mod_stub.f90",
):
    check(
        f"int-module stub defines JFNK_BA_BUILT ({stub})",
        "JFNK_BA_BUILT" in (ROOT / stub).read_text(),
    )

tri_fail = section(BLKBAND, "9999\tCONTINUE", "\tFLAG=.FALSE.")
check(
    "TRI failure at 9999 retries via the DIAG solve",
    "TRI_RECALL" in tri_fail
    and "GO TO 8888" in tri_fail
    and "CALL CMF_DIAG_SOLVE()" in tri_fail
    and "TRI_ACTIVE" in tri_fail,
)
check(
    "TRI fallback clears the partial factor cache",
    "'BCSCRATCH'" in tri_fail and "STATUS='DELETE'" in tri_fail,
)
check(
    "BCD recall gates loop over all depths",
    re.search(
        r"DO DEPTH_INDX=1,ND\s+WRITE\(BCD_F1.{0,60}CSCRATCH",
        BLKBAND,
        re.S,
    )
    is not None
    and re.search(
        r"DO DEPTH_INDX=1,ND\s+WRITE\(BCD_F1.{0,60}BCSCRATCH",
        BLKBAND,
        re.S,
    )
    is not None
    and "DSCRATCH" in section(BLKBAND, "BCD_TRI_OK=.FALSE.", "END IF"),
)
check(
    "NUM_BNDS<3 forces METH_SOL='DIAG'",
    "IF(NUM_BNDS .LT. 3 .AND. METH_SOL .NE. 'DIAG')THEN" in RD_CONTROL,
)

# JFNK_DEFER_PROBE (21-Oct-2026): the deferral-path operator probe must
# run only on FD-gate deferrals, inside the allocated-vector scope, and
# must never revive the admission or leak the h multiplier.

defer = section(
    SOLVE_CORR,
    "Deferral-path operator probe ([JFNK_DEFER_PROBE]",
    "DEFER-PROBE end.'",
)
check(
    "defer probe gated on defer + flag + cap",
    "IF(.NOT. JFNK_REACH_OK .AND. L_EVR_FDOP" in defer
    and "JFNK_DEFER_PROBE .GT. 0" in defer
    and "JFNK_DEFER_PROBE_N .LT. JFNK_DEFER_PROBE" in defer,
)
check(
    "defer probe sits between FD gates and isolation probes",
    SOLVE_CORR.index("'  DEFER-PROBE '")
    < SOLVE_CORR.index("IF(JFNK_REACH_OK)THEN\n! Isolation probes"),
)
check(
    "defer probe never revives the admission",
    "JFNK_REACH_OK=.TRUE." not in defer.replace(" ", ""),
)
check(
    "defer probe resets HSCALE after every JVP call",
    defer.count("JFNK_PROBE_HSCALE=1.0_LDP")
    >= defer.count("CALL JFNK_JVP_V1"),
)
check(
    "defer probe reuses the gate JVP at scale 1",
    "JV_FD2=JV_FD" in defer and "L_EVP(1)=.NOT. jfnk_jvp_failed" in defer,
)
check(
    "raw dF comes from the perturbed residual, not the derivative",
    "JFNK_W_F_PERT(TPI(I))-jfnk_fx(TPI(I))" in defer,
)
check(
    "JFNK_DEFER_PROBE parsed with 0=off default",
    "JFNK_DEFER_PROBE=0" in RD_CONTROL
    and "'JFNK_DEFER_PROBE'" in RD_CONTROL
    and "INTEGER JFNK_DEFER_PROBE" in CVARMOD,
)
check(
    "probe state lives in the interface module",
    "JFNK_PROBE_H_LAST" in INTMOD
    and "JFNK_BA_SRC_IT" in INTMOD
    and "JFNK_DEFER_PROBE_N" in INTMOD,
)
check(
    "actual FD step exported after retry loop",
    "JFNK_PROBE_H_LAST=H_STEP" in INTERFACE
    and "JFNK_PROBE_H_LAST=0.0_LDP" in INTERFACE,
)
check(
    "BA source iteration stamped at build",
    "JFNK_BA_SRC_IT=MAIN_COUNTER" in POPS,
)
check(
    "unknown BA src age guarded (src=-1 on disk load)",
    "JFNK_BA_SRC_IT .GE. 0" in defer
    and "unknown (loaded from disk this run)" in defer,
)
check(
    "raw-direction write keeps A/ES pairing",
    "raw direction:',5," in defer
    and "'|BAz|/|F|','|Az|/|F|','cos(F,Az)'" in defer
    and "'cos(BAz,Az)','|hv|inf'" in defer,
)
check(
    "no implicit-interface call to OPTIONAL-arg FACTOR_PRECOND",
    "JFNK_FACTOR_PRECOND_V1()" not in SOLVE_CORR
    and "JFNK_FACTOR_PRECOND_V1()" not in INTERFACE
    and "JFNK_FACTOR_PRECOND_V1(.TRUE.)" in SOLVE_CORR
    and "JFNK_FACTOR_PRECOND_V1(.TRUE.)" in INTERFACE,
)

jvp_diag = section(
    INTERFACE,
    "Extended audit (30-Oct-2026",
    "END SUBROUTINE JFNK_JVP_V1",
)
check(
    "JVP diag h-scan uses the same masked direction",
    "HMULT*H_STEP*" in jvp_diag
    and "DBLE(JFNK_W_V_USE(I))" in jvp_diag,
)
check(
    "JVP diag h-scan forms raw residual difference at its own h",
    "JFNK_W_F_PERT(I)-jfnk_fx(I))/(HMULT*H_STEP)" in jvp_diag,
)
check(
    "JVP diag FD at each h is masked like the host JV",
    "JFNK_MASK_FIXED_UNKNOWNS(N,JV_H)" in jvp_diag,
)
check(
    "JVP diag nbad guard is an integer compare",
    "jfnk_resid_nbad .GT. 0" in jvp_diag,
)
check(
    "JVP diag reports worst |FD-BA*v| components",
    "worst |FD-BA*v| components" in jvp_diag
    and "max|BA_row|" in jvp_diag,
)
check(
    "JVP diag extended block stays one-shot unless HADAPT re-arms",
    "JFNK_JVP_DIAG_DONE = .TRUE." in INTERFACE
    and ".NOT. JFNK_JVP_DIAG_DONE .AND. JFNK_LUER .GT. 0 .AND." in INTERFACE,
)
check(
    "JVP dir norm line is VADAT-capped unless HADAPT re-arms",
    "JVP_DIR_CT .LT. JFNK_JVP_NDIR .OR. JFNK_HADAPT_REARM" in INTERFACE
    and "RD_STORE_INT(JFNK_JVP_NDIR,'JFNK_JVP_NDIR'" in RD_CONTROL
    and "INTEGER JFNK_JVP_NDIR" in CVARMOD,
)

check(
    "JVP diag re-evaluates residual at the baseline state",
    "JFNK_W_X_PERT(I)=jfnk_x(I)" in jvp_diag
    and "baseline re-eval: max|F(x0)-fx|" in jvp_diag,
)
check(
    "JVP diag probes the negative direction for switch asymmetry",
    "jfnk_x(I)-H_STEP*DBLE(JFNK_W_V_USE(I))" in jvp_diag
    and "FD(-h)=" in jvp_diag,
)
check(
    "JVP diag repeat-eval compares against reconstructed first eval",
    "jfnk_fx(I)+H_STEP*JV(I)" in jvp_diag
    and "repeat eval: max|F_rep-F_first|" in jvp_diag,
)
check(
    "JVP diag sweep reports raw dF and FD at each h",
    "DIFF_ARR(I)=JFNK_W_F_PERT(I)-jfnk_fx(I)" in jvp_diag
    and "max|DF|" in jvp_diag
    and "h-sweep on same direction" in jvp_diag,
)
check(
    "JFNK_EPS floor applies in both modes (VADAT-opt-in)",
    "EPS0 = MAX(EPS0, JFNK_EPS_JVP_LOC)" in INTERFACE
    and "IF(JFNK_COUPLED_LOC)EPS0" not in INTERFACE,
)
check(
    "JFNK_CLIP_RETRY is a VADAT key with default 2",
    "INTEGER JFNK_CLIP_RETRY" in CVARMOD
    and "RD_STORE_INT(JFNK_CLIP_RETRY,'JFNK_CLIP_RETRY'" in RD_CONTROL,
)
check(
    "JVP step is a deterministic function of (x, v-hat)",
    "DO NRETRY" not in jvp
    and "IF(JFNK_N_NEW_CLIP .GT. 0)THEN" in jvp
    and "unprevented constraint crossings -- keeping FD as-is" in jvp,
)
check(
    "JVP caps every counted crossing before the eval",
    "ELSE IF(IDX_J .EQ. JFNK_NT)THEN" in jvp
    and "IDX_J .LE. JFNK_NT-1" in jvp
    and "1.0D-100/DBLE(POPS_BASE(IDX_J,K_J))" in jvp,
)
check(
    "PC audit splits |Mz-v| by reg-substituted vs non-reg rows",
    "on reg rows=" in INTERFACE
    and "on non-reg pop rows=" in INTERFACE
    and "weighted |F|^2 share on non-reg rows" in INTERFACE,
)
check(
    "PC audit checks apply against the factored matrix",
    "factored M" in INTERFACE,
)
check(
    "JVP dir report includes M^-1 J v round-trip on reachable rows",
    "rt FD: |M^-1Jv-v|/|v|" in INTERFACE
    and "jfnk_pc_clip_event = CLIP_SAV" in INTERFACE,
)
check(
    "JFNK_GMRES_STAG is a VADAT key wired to jfnk_gmres_stag_tol",
    "REAL(KIND=LDP) JFNK_GMRES_STAG" in CVARMOD
    and "RD_STORE_DBLE(JFNK_GMRES_STAG,'JFNK_GMRES_STAG'" in RD_CONTROL
    and "jfnk_gmres_stag_tol=JFNK_GMRES_STAG" in SOLVE_CORR,
)
check(
    "JFNK_BA_JVP is a VADAT key gated on BA build",
    "LOGICAL JFNK_BA_JVP" in CVARMOD
    and "RD_STORE_LOG(JFNK_BA_JVP,'JFNK_BA_JVP'" in RD_CONTROL
    and "IF(JFNK_BA_JVP .AND. JFNK_BA_BUILT" in INTERFACE,
)
check(
    "phase-1 JVP diagnostics are opt-in VADAT keys",
    "LOGICAL JFNK_LINTEST" in CVARMOD
    and "LOGICAL JFNK_HSCAN" in CVARMOD
    and "LOGICAL JFNK_HADAPT" in CVARMOD
    and "REAL(KIND=LDP) JFNK_TRES_FAC" in CVARMOD
    and "RD_STORE_LOG(JFNK_LINTEST,'JFNK_LINTEST'" in RD_CONTROL
    and "RD_STORE_LOG(JFNK_HSCAN,'JFNK_HSCAN'" in RD_CONTROL
    and "RD_STORE_LOG(JFNK_HADAPT,'JFNK_HADAPT'" in RD_CONTROL
    and "RD_STORE_DBLE(JFNK_TRES_FAC,'JFNK_TRES_FAC'" in RD_CONTROL,
)
check(
    "BUG-001: true-residual verification cannot be switched off",
    "jfnk_gmres_tres_fac = 1d0" in SOLVER
    and "acc_fac = jfnk_gmres_tres_fac" in SOLVER
    and "if(acc_fac <= 0d0) acc_fac = 1d0" in SOLVER
    and "if(jfnk_gmres_tres_fac > 0d0) then" not in SOLVER
    and "JFNK_TRES_FAC=1.0_LDP" in RD_CONTROL
    and "JFNK_TRES_FAC .LE. 0.0_LDP" in RD_CONTROL,
)
check(
    "banded LU factorization requires the LU preconditioner",
    "JFNK_USE_BANDED_LOC .AND. JFNK_NUM_BNDS .GT. 1 .AND." in INTERFACE
    and "JFNK_PRECOND_TYPE(1:2) .EQ. 'LU'" in INTERFACE,
)
check(
    "JFNK_GMRES_RST is a VADAT key wired to jfnk_gmres_max_restart",
    "INTEGER JFNK_GMRES_RST" in CVARMOD
    and "RD_STORE_INT(JFNK_GMRES_RST,'JFNK_GMRES_RST'" in RD_CONTROL
    and "jfnk_gmres_max_restart=JFNK_GMRES_RST" in SOLVE_CORR
    and "jfnk_gmres_max_restart = 1" in SOLVER
    and "n_restart > jfnk_gmres_max_restart" in SOLVER,
)


sweep = section(
    CMFGEN_SUB,
    "SUBROUTINE JFNK_RAD_SWEEP_V1",
    "END SUBROUTINE JFNK_RAD_SWEEP_V1",
)
check(
    "BUG-006: every guarded entry save precedes the first mutation",
    sweep.index("CALL JFNK_EVAL_LSTATE_V1('S'")
    < sweep.index("ACCESS_F=JFNK_EDD_START_REC")
    and sweep.index("ACCESS_F=JFNK_EDD_START_REC")
    < sweep.index("JFNK_IN_EVAL=.TRUE.")
    and sweep.rindex("GOTO 9999", 0, sweep.index("ACCESS_F=JFNK_EDD_START_REC"))
    > -1
    and "GOTO 9999" not in section(sweep, "ACCESS_F=JFNK_EDD_START_REC", "FLUSH(LUER)"),
)


# Issue-E regression: the PTC fold writes the effective (regularized)
# diagonal into JFNK_BA_DIAG; physical metrics must price against the
# raw pre-fold diagonal or dtau shrinkage falsely reports convergence.
check(
    "Issue E: raw BA diagonal declared in int module",
    "JFNK_BA_DIAG_RAW(:,:)" in INTMOD
    and "JFNK_BA_DIAG_RAW_READY" in INTMOD,
)
check(
    "Issue E: factor fills RAW before the fold",
    "JFNK_BA_DIAG_RAW(J, K) = JFNK_C_MAT(J, J)" in factor
    and "JFNK_BA_DIAG_RAW_READY = .TRUE." in factor
    and factor.index("JFNK_BA_DIAG_RAW(J, K) = JFNK_C_MAT(J, J)")
    < factor.index("JFNK_BA_A(J,J,K) = JFNK_BA_A(J,J,K) + S_PTC"),
)
mce = section(
    INTERFACE,
    "jfnk_maxch_est = -1.0_LDP",
    "Pseudo-time term",
)
check(
    "Issue E: MAXCH prices the raw diagonal",
    "JFNK_BA_DIAG_RAW(J,K)" in mce
    and "JFNK_BA_DIAG_RAW_READY" in mce,
)
for stub in (
    "tests/f90/jfnk_int_mod_stub.f90",
    "tests/jfnk_core/jfnk_int_stub.f90",
    "unit_tests/jfnk_int_mod_stub.f90",
):
    check(
        f"int-module stub defines JFNK_BA_DIAG_RAW ({stub})",
        "JFNK_BA_DIAG_RAW" in (ROOT / stub).read_text(),
    )

# Schur S3 regression: a DPOTRS solve through a failed DPOTRF factor
# returns info=0 with a nonfinite result; the shared block solve and
# the schur form/apply must flag it rather than hand the Newton loop
# a garbage direction.
check(
    "S3: blk lsq flags nonfinite DPOTRS result",
    "ieee_is_finite(u)" in SOLVER and "info_out = np + 1" in SOLVER,
)
check(
    "S3: schur form rejects nonfinite S before DGETRF",
    "ieee_is_finite(s)" in SOLVER,
)
check(
    "S3: schur apply flags nonfinite output",
    "ieee_is_finite(v)" in SOLVER,
)
check(
    "jfnk_core links the kv-logging companion object",
    "iter_log_fmt_mod.o" in (ROOT / "tests/jfnk_core/Makefile").read_text(),
)

# Phase 2: the folded PTC diagonal must carry the u-coordinate map
# factor exp(-x) for log unknowns -- the same factor the
# JFNK_PRECOND_V1 apply rescale uses. test_ptc_coords covers the
# numerical identity at nonzero iterates; this pins the fold site.
check(
    "Phase 2: PTC fold uses the u-coordinate exp(-x) map",
    "S_PTC = S_PTC * EXP(-X_PTC)" in factor
    and "jfnk_ptc_scale(IDX_PTC)) /" in factor,
)
check(
    "Phase 2: strip routines anchor at x0, not zero",
    "(x(i) - jfnk_ptc_x0(i)) * jfnk_ptc_scale(i) / &" in SOLVER,
)
check(
    "ptc-coords regression test wired into jfnk_core",
    "test_ptc_coords" in (ROOT / "tests/jfnk_core/Makefile").read_text(),
)

# Phase 5: per-depth backward-error certificate + LSQ fallback. The
# certificate solves M_s z = M_s v on a fixed probe and rejects when the
# componentwise relative residual exceeds the tol -- it is the only
# acceptance test on the LSQ fallback attempt. JFNK_DEPTH_LSQ records
# which depths took the normal-equation factor so the apply dispatches
# DPOTRS there instead of DGETRS.
check(
    "Phase 5: cert/lsqfb flags declared as VADAT keys",
    "JFNK_PC_CERT" in CVARMOD
    and "JFNK_PC_CERT_TOL" in CVARMOD
    and "JFNK_PC_LSQFB" in CVARMOD,
)
check(
    "Phase 5: VADAT parses and routes the new keys",
    "RD_STORE_LOG(JFNK_PC_CERT" in RD_CONTROL
    and "RD_STORE_DBLE(JFNK_PC_CERT_TOL" in RD_CONTROL
    and "RD_STORE_LOG(JFNK_PC_LSQFB" in RD_CONTROL
    and "JFNK_PC_CERT, JFNK_PC_CERT_TOL, JFNK_PC_LSQFB" in RD_CONTROL,
)
check(
    "Phase 5: solver module carries the flags",
    "jfnk_pc_cert" in SOLVER
    and "jfnk_pc_cert_tol" in SOLVER
    and "jfnk_pc_lsqfb" in SOLVER,
)
check(
    "Phase 5: certificate runs on demand and on the fallback attempt",
    "jfnk_pc_cert .or. itry .gt. itry_max" in SOLVER,
)
check(
    "Phase 5: LSQ fallback is the extra ladder slot",
    "lsq_try = lsq_mode .or. (itry .gt. itry_max)" in SOLVER
    and "if(lsq_fb_l) itry_top = itry_max + 1" in SOLVER,
)
check(
    "Phase 5: cert measures the factored block, not A",
    "s_tmp + mtst(i,j) * x_tst(j)" in SOLVER,
)
check(
    "Phase 5: per-depth LSQ dispatch stored in int module",
    "JFNK_DEPTH_LSQ(:)" in INTMOD,
)
check(
    "Phase 5: apply dispatches DPOTRS on LSQ-fallback depths",
    "JFNK_DEPTH_LSQ(K)" in INTERFACE
    and "LSQ_MODE_L = LSQ_MODE_L .OR. JFNK_DEPTH_LSQ(K)" in INTERFACE,
)
check(
    "Phase 5: tests exercise cert and fallback",
    "Q8: lsqfb" in (ROOT / "tests/jfnk_core/test_pc_quality.f90").read_text()
    and "lsq_u" in (ROOT / "tests/jfnk_core/test_pc_quality.f90").read_text(),
)
