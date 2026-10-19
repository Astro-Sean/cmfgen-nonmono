#!/usr/bin/env python3
"""Check JFNK control-flow guards that need CMFGEN state to run directly."""

from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[2]
INTERFACE = (ROOT / "new_main/mod_subs/jfnk_cmfgen_interface_v1.f").read_text()
SOLVER = (ROOT / "new_main/mod_subs/jfnk_solver_v1.f90").read_text()
POPS = (ROOT / "new_main/mod_subs/solve_for_pops.f").read_text()
SOLVE_CORR = (ROOT / "new_main/mod_subs/jfnk_solve_corr_v1.f").read_text()


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
    < refresh.index("CALL JFNK_FACTOR_PRECOND_V1()"),
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
    "CALL JFNK_FACTOR_PRECOND_V1()",
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
    < retry.index("CALL JFNK_FACTOR_PRECOND_V1()"),
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

for stub in (
    "tests/f90/jfnk_int_mod_stub.f90",
    "tests/jfnk_core/jfnk_int_stub.f90",
    "unit_tests/jfnk_int_mod_stub.f90",
):
    check(
        f"int-module stub defines JFNK_BA_BUILT ({stub})",
        "JFNK_BA_BUILT" in (ROOT / stub).read_text(),
    )
