#!/usr/bin/env python3
"""Check JFNK control-flow guards that need CMFGEN state to run directly."""

from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[2]
INTERFACE = (ROOT / "new_main/mod_subs/jfnk_cmfgen_interface_v1.f").read_text()
SOLVER = (ROOT / "new_main/mod_subs/jfnk_solver_v1.f90").read_text()
POPS = (ROOT / "new_main/mod_subs/solve_for_pops.f").read_text()
SOLVE_CORR = (ROOT / "new_main/mod_subs/jfnk_solve_corr_v1.f").read_text()
INTMOD = (ROOT / "new_main/mod_subs/jfnk_cmfgen_int_mod.f").read_text()
CVARMOD = (ROOT / "new_main/mod_subs/control_variable_mod.f").read_text()
BLKBAND = (ROOT / "new_main/subs/cmf_blkband_v3.f").read_text()
RD_CONTROL = (ROOT / "new_main/mod_subs/rd_control_variables.f").read_text()


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

jvp = section(INTERFACE, "SUBROUTINE JFNK_JVP_V1", "END SUBROUTINE JFNK_JVP_V1")
bad_residual = section(
    jvp,
    "CALL JFNK_RESIDUAL_V1(N, JFNK_W_X_PERT, JFNK_W_F_PERT)",
    "IF(JFNK_N_NEW_CLIP .EQ. 0)EXIT",
)
check(
    "JVP rejects repaired perturbed residuals before differencing",
    "jfnk_resid_nbad .GT. 0 .OR. JFNK_EVAL_BROKEN" in bad_residual
    and "jfnk_jvp_failed = .TRUE." in bad_residual
    and "JV(I) = (JFNK_W_F_PERT(I) - jfnk_fx(I)) / H_STEP"
    not in bad_residual,
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
    and "BA src it=unknown" in defer,
)
check(
    "raw-direction write keeps A/ES pairing",
    "'(A,ES12.4,4(A,ES12.4))'" in defer
    and "'(A,5(A,ES12.4))'" not in defer,
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
    "JVP diag extended block stays one-shot behind JVP_DIAG_DONE",
    "JVP_DIAG_DONE = .TRUE." in INTERFACE
    and ".NOT. JVP_DIAG_DONE .AND. JFNK_LUER .GT. 0 .AND." in INTERFACE,
)
check(
    "JVP dir norm line is VADAT-capped by JFNK_JVP_NDIR",
    "JVP_DIR_CT .LT. JFNK_JVP_NDIR .AND. JFNK_LU_READY" in INTERFACE
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
    "clip retry keeps last attempt at its own h",
    "DO NRETRY = 0, JFNK_CLIP_RETRY" in INTERFACE
    and "NRETRY .EQ. JFNK_CLIP_RETRY" in INTERFACE
    and "keeping FD as-is" in INTERFACE,
)
check(
    "PC audit splits |Mz-v| by reg-substituted vs non-reg rows",
    "on reg rows=" in INTERFACE
    and "on non-reg pop rows=" in INTERFACE
    and "weighted |F|^2 share on non-reg rows" in INTERFACE,
)
