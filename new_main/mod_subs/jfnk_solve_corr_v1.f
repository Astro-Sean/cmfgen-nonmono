!
! JFNK_SOLVE_CORR_V1: Jacobian-Free Newton-Krylov population correction
! solve, replacing the direct SOLVEBA_V13 call when USE_JFNK=T.
!
! Ported from V6 SOLVE_FOR_POPS_JFNK_V1 (bak_v6revert) and slimmed to a
! drop-in replacement for the single SOLVEBA call site in
! SOLVE_FOR_POPS.  The outer mode controller (LAMBDA/full switching,
! GS field-maturity gating, hold_pops, resacc) is unchanged; this
! routine only replaces the correction solve itself.
!
! Method (Knoll & Keyes 2004; Arramy, de la Cruz Rodriguez &
! Leenaarts 2024, A&A 690, A12):
!   F(x) = SE residual at POPS_BASE*(1+x)
!   J*v  ~= [F(x+eps*v) - F(x)]/eps     (Jacobian-free matvec)
!   Newton step globalized by the GMRES hookstep trust region,
!   Eisenstat-Walker forcing, and an Armijo line search.
!
! JFNK_COUPLED=T (VADAT): the residual re-runs the formal radiation
! solve at the perturbed populations (JFNK_SWEEP_CB, registered by
! cmfgen_sub), so the finite-difference JVP probes dJ/dn -- the fully
! coupled SE+RTE system.  This is the Arramy-style evaluation and is
! the piece V6 never had (V6's residual froze the radiation field and
! therefore solved only the fixed-J sub-problem).
!
! JFNK_COUPLED=F: V6 semantics -- STEQ_coll(n) + stored radiative term
! + CMAT_RAD linearization.  Retained as a debugging/fallback mode.
!
! On non-catastrophic stagnation this routine returns to the caller's
! SOLVEBA_V13 mode.  Only an unrecoverable eval failure forces the
! conservative Lambda safety path.
!
! editted by Devin (V7-JFNK port, 02-Oct-2026)
! Altered 15-Oct-2026 : BA_DIRCHK_MAX caps the dirchk defer count so
!                          real Newton admission resumes (edited by Devin).
! Altered 16-Oct-2026 : Applied-change metrics use the BA_DEM trace
!                          mask and export BA_TR_ABS for Stage B;
!                          commits also export BA_DEM_EXM so the
!                          Stage-B exit gate sees exempted demand on
!                          JFNK iterations (edited by Devin).
! Altered 05-Oct-2026 : Demoted JFNK banners to iteration subheaders.
!
	SUBROUTINE JFNK_SOLVE_CORR_V1(SOL,POPS,POP_ATOM_A,DIAG_INDX,NT,
	1    NION,NUM_BNDS,ND,MAXCH,MAXCH_APP,TEMP_CHAR,SUCCESS,
	2    SCALE_OPT_A,LAM_SC_OPT_A,COR_LIM_IN,MDT_A,TMIN_A,
	3    COMPUTE_BA_A,WR_BA_INV_A,WR_PRT_INV_A,LAMBDA_A,
	4    MAIN_COUNTER,SET_D2_A,LST_ITERATION,JFNK_COUPLED_IN)
	USE SET_KIND_MODULE
	USE MOD_CMFGEN
	USE ANG_QW_MOD
	USE CONTROL_VARIABLE_MOD
	USE STEQ_DATA_MOD
	USE JFNK_CMFGEN_INT_MOD
	USE JFNK_SOLVER_V1_MOD
	USE FG_J_CMF_NONMONO_MOD_V1, ONLY: FREEZE_I_NEXT
	IMPLICIT NONE
!
	INTEGER NT,NION,NUM_BNDS,ND,DIAG_INDX
	INTEGER MAIN_COUNTER
	REAL(KIND=LDP) POPS(NT,ND)
	REAL(KIND=LDP) SOL(NT,ND)
	REAL(KIND=LDP) POP_ATOM_A(ND)
	REAL(KIND=LDP) MAXCH,MAXCH_APP
	REAL(KIND=LDP) COR_LIM_IN,MDT_A,TMIN_A
	CHARACTER*20 TEMP_CHAR
	CHARACTER*(*) SCALE_OPT_A,LAM_SC_OPT_A
	LOGICAL COMPUTE_BA_A,WR_BA_INV_A,WR_PRT_INV_A,LAMBDA_A
	LOGICAL SET_D2_A,SUCCESS,LST_ITERATION
! JFNK_COUPLED_IN: residual mode for THIS activation, selected by the
! FRZ_GATE/COR_GATE tier logic in SOLVE_FOR_POPS.  F = frozen-J
! residual (no sweep per eval); T = coupled eval (fresh formal solve
! per eval).  Stored in JFNK_COUPLED_LOC so residual/JVP/rebuild all
! read the same per-activation mode.
	LOGICAL JFNK_COUPLED_IN
!
! Local data
!
	REAL(KIND=LDP) POPS_PREV(NT,ND)
	REAL(KIND=LDP) Z_POP(NT)
	REAL(KIND=LDP) T1,T2
	INTEGER I,J,K,ID
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
	CHARACTER*11 JFNK_TAG
!
! Phase timers (edited by Devin, 15-Oct-2026): per-admission wall
! seconds for the expensive JFNK phases, printed as one line with the
! final diagnostics so the iteration cost can be attributed (the
! it191 admission spent ~500 s in dirchk evals that nothing time-
! stamped).  Bracketed with OMP_GET_WTIME, matching the iteration
! timer in CMFGEN_SUB.
!
	REAL(KIND=LDP) OMP_GET_WTIME
	EXTERNAL OMP_GET_WTIME
	REAL(KIND=LDP) T_PH_TOT,T_PH_CMAT,T_PH_PC,T_PH_F0
	REAL(KIND=LDP) T_PH_PREC,T_PH_NWT,TW_PH
!
! Solver settings (validated copies of the control variables)
!
	INTEGER JFNK_M_GMRES_LOC,JFNK_MAX_NEWTON_LOC
	REAL(KIND=LDP) JFNK_GTOL_LOC,JFNK_TOL_REL_LOC
	REAL(KIND=LDP) JFNK_DEL_INIT_LOC,JFNK_MNDL_LOC,JFNK_MXDL_LOC
	REAL(KIND=LDP) JFNK_EPS_LOC
!
! Newton solve state
!
	INTEGER JFNK_N,JFNK_INFO,JFNK_ATTEMPT,II_F
	LOGICAL APPLY_JFNK_UPDATE,JFNK_CATASTROPHIC
	LOGICAL BASELINE_REBUILT,F0_OK,DO_JFNK_NOW
	LOGICAL COMMITTED_BASELINE,VERIF_REJECT
	LOGICAL LAMBDA_SAVE,FIXED_T_SAVE,COMPUTE_BA_SAVE
	LOGICAL WR_BA_INV_SAVE,WR_PRT_INV_SAVE
	LOGICAL JFNK_LAG_FACTORS_OK
	CHARACTER*20 TEMP_CHAR_SAVE
	DOUBLE PRECISION F_INF,F_INF_INIT,F_W,F_W_INIT
	DOUBLE PRECISION JFNK_TOL_ABS,F_APP
	DOUBLE PRECISION, ALLOCATABLE :: BEST_X(:)
	DOUBLE PRECISION BEST_F_W,BEST_F_INF
	LOGICAL JFNK_REACH_OK,JFNK_REACH_SKIP
	DOUBLE PRECISION REACH_VAL,FF2_DP,FJZ_DP,JJZ_DP
	DOUBLE PRECISION COS_BA,JZF_BA,CLIP_FR,COS_FD,JZF_FD
	DOUBLE PRECISION X_MAX_R,X_MIN_R,X_T_MAX_R,X_T_MIN_R
	DOUBLE PRECISION ZMAX_RAW
	INTEGER NCLIP_R
	DOUBLE PRECISION JVP_N
	DOUBLE PRECISION FJZ_FD,JJZ_FD,BAFD_DP
	DOUBLE PRECISION, ALLOCATABLE :: Z_CHK(:),JZ_VEC(:),JV_FD(:),Z_CLIP(:)
	DOUBLE PRECISION, ALLOCATABLE :: JV_FD2(:),Z_BLK(:),JV_SAV(:)
	INTEGER K_TOPF,I_TOPF,I_ZMAX,II_PROBE
	DOUBLE PRECISION FMAX_LOC
!
! COLPROBE decomposition (26-Oct-2026): for each worst element print
! the row residual F_i (self-proportional equations have dF/dln n_i
! = -F_i, i.e. FD~O(F) while BA shows the raw outflow), the CMAT_RAD
! share of the BA element, the ion/level label, and the frozen
! REPLACE flag.  A separate h-scan on the topZ self element separates
! a small real derivative (FD stable vs h) from cancellation noise.
!
	INTEGER, ALLOCATABLE :: J2ID_LOC(:)
	INTEGER, ALLOCATABLE :: NZ_SP(:)
	INTEGER JROW,KROW,LEV_LOC,IHS
	DOUBLE PRECISION POP_LOC,NLTE_LOC,BADG_LOC
	CHARACTER*14 LBL_LOC
	LOGICAL REP_LOC
	DOUBLE PRECISION CMR_LOC,COLL_LOC
	DOUBLE PRECISION HS_ARR(3)
	DATA HS_ARR/1.0D2,1.0D4,1.0D6/
!
! Correction clipping state
!
	INTEGER NCLIP_POP,NCLIP_T,NCLIP_ED
	INTEGER IINC,IDEC
	INTEGER, PARAMETER :: NVJ=10
	REAL(KIND=LDP) JFNK_INC_VEC(NVJ),JFNK_DEC_VEC(NVJ)
	REAL(KIND=LDP) JFNK_MAXINC,JFNK_MAXDEC
	REAL(KIND=LDP) X_MAX,X_MIN,X_T_MAX,X_T_MIN
	REAL(KIND=LDP) MINSCALE_J,T1_RAW
!
	LOGICAL, SAVE :: JFNK_BANNER_PRINTED=.FALSE.
	INTEGER, SAVE :: JFNK_STAGNATION_COUNT=0
	INTEGER, PARAMETER :: JFNK_MAX_STAGNATE=3
	LOGICAL, SAVE :: JFNK_DISABLED=.FALSE.
!
! Cooldown re-arm: JFNK is never permanently disabled.  After
! JFNK_MAX_STAGNATE consecutive failures it pauses for
! JFNK_REARM_WAIT eligible calls, which use the outer SOLVEBA mode,
! then retries.  The wait doubles each failure cycle up to x8 so a
! persistently failing Newton is retried rarely; it resets on any
! committed JFNK update.
!
	INTEGER, SAVE :: JFNK_REARM_COUNT=0
	INTEGER, SAVE :: JFNK_REARM_CYCLE=0
	INTEGER JFNK_REARM_WAIT
!
! BA_DIRCHK_MAX arms the dirchk diagnostics for at most that many
! defer events (default 1; -1 = unlimited, the pre-cap behaviour).
! Each armed admission costs a fresh BA factorization plus ~12
! residual/JVP evals; past the cap JFNK admission proceeds normally.
! JFNK_DIRCHK_N is a CONTROL_VARIABLE_MOD counter (initialized at
! VADAT read) so the admission bypass in SOLVE_FOR_POPS drops out on
! the same event count.
!
	LOGICAL DO_DIRCHK_NOW
!
! FIX_T_AUTO scratch: the domain T-row residual norm at admission and
! the committed log-T step measured at the outcome.
!
	DOUBLE PRECISION TAU_RE_NRM,TAU_EPS_T
!
	EXTERNAL JFNK_INIT_V1
	EXTERNAL JFNK_RESIDUAL_V1,JFNK_JVP_V1,JFNK_PRECOND_V1
	EXTERNAL JFNK_DIAGNOSTIC_V1
	EXTERNAL JFNK_ASSEMBLE_STEQ_V1,JFNK_COMPUTE_STEQ_V1
	EXTERNAL JFNK_SYNC_ATM_FROM_POPS
	EXTERNAL JFNK_CAPTURE_LAG_FACTORS_V1
	EXTERNAL JFNK_MASK_FIXED_UNKNOWNS,JFNK_COMPUTE_RSCALE_V1
	EXTERNAL JFNK_BUILD_BA_RAD_V1,JFNK_FACTOR_PRECOND_V1
	EXTERNAL JFNK_REBUILD_BASELINE_V1
	EXTERNAL JFNK_SET_PTC_SCALE_V1
	EXTERNAL JFNK_SAVE_MOD_STEQ_V1,JFNK_RESTORE_MOD_STEQ_V1
	EXTERNAL JFNK_ROW_ANATOMY_V1
	DOUBLE PRECISION JFNK_DOT_V1
	EXTERNAL JFNK_DOT_V1
!
	LUER=ERROR_LU()
	JFNK_LAG_FACTORS_ACTIVE=.FALSE.
	JFNK_LAG_FACTORS_READY=.FALSE.
	VERIF_REJECT=.FALSE.
	JFNK_REACH_SKIP=.FALSE.
	LAMBDA_SAVE=LAMBDA_A
	FIXED_T_SAVE=FIXED_T
	COMPUTE_BA_SAVE=COMPUTE_BA_A
	WR_BA_INV_SAVE=WR_BA_INV_A
	WR_PRT_INV_SAVE=WR_PRT_INV_A
	TEMP_CHAR_SAVE=TEMP_CHAR
	T1=COR_LIM_IN
	DO_DIRCHK_NOW=DO_BA_DIR_CHK .AND.
	1    (BA_DIRCHK_MAX .LT. 0 .OR. JFNK_DIRCHK_N .LT. BA_DIRCHK_MAX)
!
	T_PH_TOT=OMP_GET_WTIME()
	T_PH_CMAT=0.0_LDP
	T_PH_PC=0.0_LDP
	T_PH_F0=0.0_LDP
	T_PH_PREC=0.0_LDP
	T_PH_NWT=0.0_LDP
!
! -----------------------------------------------------------------------
! Entry gating: after repeated stagnation JFNK pauses for
! JFNK_REARM_WAIT admitted calls, then re-arms.  Paused and failed
! non-catastrophic calls retain the caller's SOLVEBA mode; only a
! catastrophic eval uses the conservative Lambda safety solve.
! -----------------------------------------------------------------------
!
! JFNK_ANY_ITER (VADAT [JFNK_ANYITER], default F): bring-up flag --
! the call site fires JFNK on a LAMBDA-designated call too, so treat
! the call as a full solve for dispatch purposes.  The committed
! update still passes through solve_for_pops' hold/reject logic.
!
	DO_JFNK_NOW=(.NOT. LAMBDA_A) .OR. JFNK_ANY_ITER
!
	IF(JFNK_DISABLED .AND. DO_JFNK_NOW)THEN
	  JFNK_REARM_COUNT=JFNK_REARM_COUNT+1
	END IF
	JFNK_REARM_WAIT=JFNK_REARM_EVERY*
	1   (2**MIN(JFNK_REARM_CYCLE,3))
	IF(JFNK_DISABLED .AND. DO_JFNK_NOW .AND.
	1   JFNK_REARM_COUNT .GE. JFNK_REARM_WAIT)THEN
	  JFNK_DISABLED=.FALSE.
	  JFNK_STAGNATION_COUNT=0
	  JFNK_REARM_COUNT=0
	  WRITE(LUER,'(A,I3,A,I3,A)')
	1   'JFNK re-armed after ',JFNK_REARM_WAIT,
	2   ' skipped calls (cycle ',JFNK_REARM_CYCLE,
	3   '); retrying Newton'
	END IF
!
	IF(JFNK_DISABLED .AND. DO_JFNK_NOW)THEN
	  WRITE(LUER,'(A,I3,A,I3,A)')
	1   'JFNK paused (cooldown ',JFNK_REARM_COUNT,' of ',
	2   JFNK_REARM_WAIT,'); using caller SOLVEBA mode'
! The paused path has not changed POPS or the stored BA.  Skip Newton
! while preserving the outer solver's mode, matrix and correction cap.
	  DO_JFNK_NOW=.FALSE.
	END IF
!
	IF(DO_JFNK_NOW)THEN
!
! =======================================================================
! JFNK solve
! =======================================================================
!
! Per-activation residual mode from the tier gate.  The frozen path
! (JFNK_COUPLED_LOC=F) builds the stored radiative term + CMAT_RAD
! matvec below; the coupled path re-solves radiation inside each
! residual eval.  JFNK_COUPLED_IN=F also covers the case where the
! coupled tier was armed ([JFNK_COUPLED]=T) but only the frozen gate
! was satisfied this call.
!
	  JFNK_COUPLED_LOC=JFNK_COUPLED_IN .AND. JFNK_COUPLED
!
! Tag used by the SOLVEBA-style max-change prints so an OUTGEN
! reader can tell which residual tier produced the committed
! correction (frozen-J = cheap stored-term evals, coupled = each
! eval a full formal solve).
!
	  IF(JFNK_COUPLED_LOC)THEN
	    JFNK_TAG=' (JFNK-CPL)'
	  ELSE
	    JFNK_TAG=' (JFNK-FRZ)'
	  END IF
!
! Fresh-residual flag for the Newton entry eval (see
! jfnk_solver_v1_mod).  Cleared here so a stale flag cannot survive
! between activations; set after the F(0) test and by the baseline
! rebuild, both of which leave jfnk_fx = F(jfnk_x).
!
	  jfnk_fx_fresh=.FALSE.
!
	  WRITE(LUER,'(/,A)')' -------------------------------- JFNK SOLVER --------------------------------'
	  WRITE(LUER,'(/,A)')' --------------------------------- JFNK CONTEXT --------------------------------'
	  WRITE(LUER,'(A,I5)')
	1   '    great iteration                 :',MAIN_COUNTER
	  WRITE(LUER,'(A,A)')
	1   '    admitted tier                   : ',
	1     TRIM(ADJUSTL(JFNK_TAG))
! FIXED_T is the IN_ITS schedule flag; a FULL token leaves it false
! even when RD_FIX_T pins every T row via MOD_FIXED_T.  Report the
! effective state (schedule .OR. VADAT pin) so the banner cannot claim
! free T on a globally pinned solve.
!
	  WRITE(LUER,'(A,L2)')
	1   '    T free to vary                  :',.NOT. FIXED_T
	2   .AND. .NOT. RD_FIX_T
	  WRITE(LUER,'(A,L2)')
	1   '    Ne free to vary                 :',.NOT. FIXED_NE
	  WRITE(LUER,'(A,L2)')
	1   '    FIX_T set in VADAT              :',RD_FIX_T
	  WRITE(LUER,'(A,L2)')
	1   '    lag dissolution factors         :',JFNK_LAG_FACTORS
	  WRITE(LUER,'(A,ES10.3)')
	1   '    population change clip          :',MAX_JFNK_COR
	  WRITE(LUER,'(A,ES10.3)')
	1   '    temperature change clip         :',MDT_A
	  IF(.NOT. FIXED_T .AND. .NOT. RD_FIX_T)THEN
	    WRITE(LUER,'(A)')
	1   '    NOTE: free-T JFNK; the frozen-J residual'
	    WRITE(LUER,'(A)')
	1   '    linearizes the radiation about T0 only.'
	  END IF
	  IF(MOD_FIX_T_D_ST .GT. 0)THEN
	    WRITE(LUER,'(A,I4,A,I4)')
	1     '    T held fixed over depths ',
	2     MOD_FIX_T_D_ST,' to ',MOD_FIX_T_D_END
	  END IF
	  FLUSH(LUER)
!
!
! Local solver settings with sanity clamps (V4 audit, 09-Sep-2026).
!
	  JFNK_M_GMRES_LOC=JFNK_M_GMRES
	  JFNK_MAX_NEWTON_LOC=JFNK_MAX_NEWTON
	  JFNK_GTOL_LOC=JFNK_GTOL
	  JFNK_TOL_REL_LOC=JFNK_TOL_REL
	  JFNK_DEL_INIT_LOC=JFNK_DEL_INIT
	  JFNK_MNDL_LOC=JFNK_MNDL
	  JFNK_MXDL_LOC=JFNK_MXDL
	  JFNK_EPS_LOC=JFNK_EPS_JVP
!
	  IF(JFNK_M_GMRES_LOC .LT. 1)JFNK_M_GMRES_LOC=10
	  IF(JFNK_MAX_NEWTON_LOC .LT. 1)JFNK_MAX_NEWTON_LOC=20
	  IF(JFNK_GTOL_LOC .LE. 0.0_LDP)JFNK_GTOL_LOC=1.0E-4_LDP
	  IF(JFNK_TOL_REL_LOC .LE. 0.0_LDP)JFNK_TOL_REL_LOC=1.0E-6_LDP
	  IF(JFNK_EPS_LOC .LE. 0.0_LDP)
	1     JFNK_EPS_LOC=SQRT(EPSILON(1.0_LDP))
	  IF(JFNK_MNDL_LOC .LE. 0.0_LDP)JFNK_MNDL_LOC=1.0E-8_LDP
	  IF(JFNK_MXDL_LOC .LT. JFNK_MNDL_LOC)JFNK_MXDL_LOC=1.0E+2_LDP
	  IF(MAX_JFNK_COR .LE. 1.0_LDP)MAX_JFNK_COR=10.0_LDP
	  IF(JFNK_DEL_INIT_LOC .EQ. 0.0_LDP)JFNK_DEL_INIT_LOC=-1.0_LDP
!
	  jfnk_max_gmres_its=JFNK_MAX_GMRES
	  jfnk_use_left_pc=JFNK_LEFT_PRECOND
	  jfnk_pc_ref_frac=DBLE(JFNK_PC_REFRESH_FRAC)
	  jfnk_need_precond_refresh=0
	  jfnk_do_line_search=JFNK_LINE_SEARCH
	  jfnk_ls_max_steps=JFNK_LS_MAX
	  jfnk_ls_beta_fac=DBLE(JFNK_LS_BETA)
	  jfnk_ls_armijo_c=DBLE(JFNK_LS_C)
	  jfnk_ls_nonmono=JFNK_LS_NM
	  jfnk_ls_nmwin_s=MIN(MAX(JFNK_LS_NMWIN,1),8)
	  jfnk_do_cauchy=JFNK_CAUCHY
	  jfnk_diagnostic_flag=JFNK_DIAGNOSTIC
!
	  JFNK_FIX_T_D_ST=MOD_FIX_T_D_ST
	  JFNK_FIX_T_D_END=MOD_FIX_T_D_END
!
! Ionic charge vector for SUP_TO_FULL_V4 calls inside the residual
! (Z_POP is a cmfgen_sub local, so rebuild it here).
!
	  Z_POP=0.0_LDP
	  DO ID=1,NUM_IONS-1
	    CALL SET_Z_POP(Z_POP,ATM(ID)%ZXzV,ATM(ID)%EQXzV,
	1       ATM(ID)%NXzV,NT,ATM(ID)%XzV_PRES)
	  END DO
!
	  CALL JFNK_INIT_V1(NT,ND,NION,NUM_BNDS,DIAG_INDX,LUER,
	1      POPS,ED,T,Z_POP,DO_LEV_DISSOLUTION,
	2      FIXED_T,TEMP_CHAR,COMPUTE_BA_A,JFNK_EPS_LOC)
	  FLUSH(LUER)
!
! FIX_T_AUTO: adaptive damped-T domain.  When JFNK_TAUTO and the
! caller left T globally free, the VARFIXT depth range (or all depths
! when no range is set) stops being a hard mask: T rows stay live in
! F/JVP and every accepted Newton step has its T components scaled by
! JFNK_TAU_ALPHA and capped at JFNK_TAU_DMAX.  This is the STELLA
! shock guard -- a STELLA T profile is not in NLTE equilibrium, so a
! free first Newton step can drive T to unphysical values.
!
! RD_FIX_T is a VADAT-level global pin the IN_ITS schedule never sees:
! a FULL token leaves FIXED_T false even though every T row is fixed
! (dirchk run, it473).  TAUTO must not release that pin -- with it
! active the depth-range mask in JFNK_MASK_FIXED_UNKNOWNS and the
! T-commit skip below are bypassed, so a Newton step would update T
! under FIX_T=T.
!
	  JFNK_TAUTO_ACTIVE=JFNK_TAUTO .AND. (.NOT. JFNK_FIXED_T_LOC)
	1   .AND. (.NOT. RD_FIX_T)
	  IF(JFNK_TAUTO_ACTIVE)THEN
	    IF(JFNK_FIX_T_D_ST .GT. 0)THEN
	      JFNK_TAU_DOM_ST=JFNK_FIX_T_D_ST
	      JFNK_TAU_DOM_END=JFNK_FIX_T_D_END
	    ELSE
	      JFNK_TAU_DOM_ST=1
	      JFNK_TAU_DOM_END=ND
	    END IF
	    jfnk_tauto_on=.TRUE.
	    jfnk_tauto_alpha=DBLE(JFNK_TAU_ALPHA)
	    jfnk_tauto_dmax=DBLE(JFNK_TAU_DMAX)
	    jfnk_tauto_nt=JFNK_NT
	    jfnk_tauto_jt=JFNK_NT
	    jfnk_tauto_dst=JFNK_TAU_DOM_ST
	    jfnk_tauto_dend=JFNK_TAU_DOM_END
	  ELSE
	    jfnk_tauto_on=.FALSE.
	    JFNK_TAU_DOM_ST=0
	    JFNK_TAU_DOM_END=0
	  END IF
!
! Freeze the lagged I_NEXT store during JFNK evaluations.  In coupled
! mode the formal solver still runs, but treats the converged store as
! a fixed parameter -- the store is outer-iteration state, refreshed
! by the next real sweep after the correction is committed.
!
	  FREEZE_I_NEXT=.TRUE.
!
! edited by Devin, 05-Oct-2026: disarm any pending BA residual-
! acceptance arm. This is reachable only with [DO_JFNK]=T (default F),
! which requires [DO_NONMONO_KNOP] or [DO_NONMONO_SOB]=T; the default
! path is unchanged.
! A full SOLVEBA arms BAACC_PENDING expecting the next SOLVEBA call
! to judge its committed step. If JFNK runs instead,
! the post-JFNK residual cannot judge that BA trial: Newton may change
! the populations before the next SOLVEBA. Keep the raw-demand gate
! closed until a later BA trial is actually judged; do not label the
! displaced trial accepted. POPS_BA_SAVE stays allocated but unread.
!
	  IF(DO_BA_RES_ACC .AND. BAACC_PENDING)THEN
	    BAACC_PENDING=.FALSE.
	    IF(BA_RAW_TRIAL_PENDING)THEN
	      BA_RAW_TRIAL_PENDING=.FALSE.
	      BA_RAW_REJECTED=.TRUE.
	    END IF
	    WRITE(LUER,'(A)')
	1     ' JFNK entry: pending BA resacc arm disarmed; raw gate held'
	  END IF
!
! Assemble the baseline STEQ (accumulated by the just-completed sweep
! at the current populations) into jfnk_fx.
!
	  JFNK_N=NT*ND
	  IF(ALLOCATED(jfnk_x))THEN
	    IF(SIZE(jfnk_x) .NE. JFNK_N)DEALLOCATE(jfnk_x)
	  END IF
	  IF(.NOT. ALLOCATED(jfnk_x))ALLOCATE(jfnk_x(JFNK_N))
	  IF(ALLOCATED(jfnk_fx))THEN
	    IF(SIZE(jfnk_fx) .NE. JFNK_N)DEALLOCATE(jfnk_fx)
	  END IF
	  IF(.NOT. ALLOCATED(jfnk_fx))ALLOCATE(jfnk_fx(JFNK_N))
	  jfnk_x=0.0D0
!
! New admission baseline: the first ASSEMBLE below re-harvests the
! REPLACE/ZERO_STEQ equation structure at this state; every perturbed
! residual eval in this Newton solve replays it frozen.
!
	  JFNK_EQ_STRUCT_READY=.FALSE.
!
	  CALL JFNK_ASSEMBLE_STEQ_V1(SOL,NT,ND,NION,NUM_BNDS,
	1     DIAG_INDX,POPS)
	  DO K=1,ND
	    DO J=1,NT
	      jfnk_fx((K-1)*NT+J)=DBLE(SOL(J,K))
	    END DO
	  END DO
!
! Frozen-J mode (JFNK_COUPLED=F): split off the radiative part of the
! assembled STEQ so the residual can supply it from storage.  In
! coupled mode the field evaluator regenerates radiative terms inside
! SE%STEQ directly, so neither STEQ_RAD nor CMAT_RAD is needed.
!
	  IF(.NOT. JFNK_COUPLED_LOC)THEN
	    TW_PH=OMP_GET_WTIME()
	    CALL JFNK_SYNC_ATM_FROM_POPS(POPS)
	    IF(JFNK_LAG_FACTORS)THEN
	      CALL JFNK_CAPTURE_LAG_FACTORS_V1(JFNK_LAG_FACTORS_OK)
	      IF(JFNK_LAG_FACTORS_OK)THEN
	        JFNK_LAG_FACTORS_ACTIVE=.TRUE.
	        WRITE(LUER,'(A)')
	1         '  JFNK: POPION/dissolution factors held at activation baseline'
	      ELSE
	        WRITE(LUER,'(A)')
	1         '  JFNK: lag-factor snapshot unavailable; using refreshed factors'
	      END IF
	    END IF
!
! Order-of-operations audit (03-Oct-2026): compute_steq rewrites the
! module STEQ arrays with collisional-only values, and nothing
! restored them -- the catastrophic/stagnation SOLVEBA fallback then
! assembled its right-hand side (and measured the resacc residual)
! from a STEQ missing every radiative term.  Bracket the split so the
! module arrays keep the outer iteration's full STEQ for the rest of
! the activation; JFNK_STEQ_RAD_FLAT already holds what JFNK needs.
!
	    CALL JFNK_SAVE_MOD_STEQ_V1()
	    CALL JFNK_COMPUTE_STEQ_V1(NT,ND,NION,NUM_BNDS,DIAG_INDX,
	1       POPS,.FALSE.,FIXED_T,LST_ITERATION)
	    CALL JFNK_ASSEMBLE_STEQ_V1(SOL,NT,ND,NION,NUM_BNDS,
	1       DIAG_INDX,POPS)
	    DO K=1,ND
	      DO J=1,NT
	        JFNK_STEQ_RAD_FLAT((K-1)*NT+J)=
	1         jfnk_fx((K-1)*NT+J)-DBLE(SOL(J,K))
	      END DO
	    END DO
	    CALL JFNK_RESTORE_MOD_STEQ_V1()
	    JFNK_RAD_READY=.TRUE.
	    CALL JFNK_BUILD_BA_RAD_V1()
	    T_PH_CMAT=OMP_GET_WTIME()-TW_PH
	  END IF
	  FLUSH(LUER)
!
! Factor the BA matrix as the GMRES preconditioner (cached per-depth
! equilibrated LU, or banded LU when JFNK_USE_BANDED=T).
!
	  TW_PH=OMP_GET_WTIME()
	  CALL JFNK_FACTOR_PRECOND_V1()
	  T_PH_PC=OMP_GET_WTIME()-TW_PH
	  FLUSH(LUER)
!
	  CALL JFNK_MASK_FIXED_UNKNOWNS(JFNK_N,jfnk_fx)
	  CALL JFNK_COMPUTE_RSCALE_V1(JFNK_N,jfnk_fx)
!
! F(0) identity test: in frozen mode the JFNK residual at x=0 must
! reproduce the outer-iteration STEQ to within rounding; otherwise the
! residual is solving a different problem than the model state.  In
! coupled mode the eval re-solves J with the I_NEXT store frozen, so
! F_eval(0) legitimately differs from outer STEQ (which accumulated
! radiative terms with the store evolving mid-sweep).  The eval
! residual is adopted as the Newton baseline in both modes on the
! F0_OK path (see below); the mismatch diagnostics stay informative
! either way.
!
	  F0_OK=.TRUE.
	  TW_PH=OMP_GET_WTIME()
	  BLOCK
	    DOUBLE PRECISION F0_X(JFNK_N),F0_Y(JFNK_N)
	    DOUBLE PRECISION F0_DIFF,F0_MAX,F0_REL
	    INTEGER JJ_TEST,F0_NAN,F0_NFAR
	    F0_X=0.0D0
	    F0_Y=0.0D0
	    CALL JFNK_RESIDUAL_V1(JFNK_N,F0_X,F0_Y)
	    CALL JFNK_MASK_FIXED_UNKNOWNS(JFNK_N,F0_Y)
	    F0_DIFF=0.0D0
	    F0_MAX=0.0D0
	    F0_NAN=0
	    DO JJ_TEST=1,JFNK_N
	      F0_DIFF=MAX(F0_DIFF,ABS(F0_Y(JJ_TEST)-jfnk_fx(JJ_TEST)))
	      F0_MAX=MAX(F0_MAX,ABS(jfnk_fx(JJ_TEST)))
	      IF(F0_Y(JJ_TEST) .NE. F0_Y(JJ_TEST) .OR.
	1        jfnk_fx(JJ_TEST) .NE. jfnk_fx(JJ_TEST))F0_NAN=F0_NAN+1
	    END DO
	    F0_REL=F0_DIFF/MAX(F0_MAX,1.0D-30)
	    WRITE(LUER,'(/,A)')' ---------------------------------- F(0) TEST ----------------------------------'
	    WRITE(LUER,'(A,ES12.4)')
	1   '    max |F_eval - F_outer|          :',F0_DIFF
	    WRITE(LUER,'(A,ES12.4)')
	1   '    max |F_outer|                   :',F0_MAX
	    WRITE(LUER,'(A,ES12.4)')
	1   '    eval vs outer STEQ, rel diff    :',F0_REL
	    IF(F0_NAN .GT. 0)THEN
	      WRITE(LUER,'(A,I6,A,I6)')
	1       '    NaN residual components         :',F0_NAN,' /',
	2       JFNK_N
	    END IF
!
! In coupled mode F_eval(x) is a different operator than the outer
! STEQ: the eval re-solves the radiation with the I_NEXT store frozen
! at its post-sweep state, while the outer STEQ accumulated radiative
! terms with the store evolving mid-sweep.  The 1.0D-6 gate is only
! meaningful in frozen mode.  For [JFNK_COUPLED] Newton must zero the
! eval operator, so adopt F0_Y as the baseline residual and report
! the outer-STEQ discrepancy (and its worst equation) as a diagnostic;
! the outer iteration closes the lagged-field difference.  Fail only
! on NaN (F0_NAN) or a non-finite eval (jfnk_resid_nbad, set inside
! JFNK_RESIDUAL_V1 before the cap).
!
! A broken eval (NaN/Inf components, capped at +/-1D20 inside
! JFNK_RESIDUAL_V1) is flagged by jfnk_resid_nbad and always aborts.
! A large-but-finite mismatch is a different situation: in coupled
! mode it is dominated by the lagged I_NEXT store having evolved
! during the production sweep, which is honest physics, not a bug.
! Newton then solves the SE problem under the store the next real
! sweep will actually use.  The eval itself is the operator to
! trust; adopting it is always correct here.  The previous
! rel>0.5 abort fired on every live attempt (rel=44, 0.82) and no
! Newton step ever ran.
!
	    IF(jfnk_resid_nbad .GT. 0)THEN
	      F0_OK=.FALSE.
	      WRITE(LUER,'(A,I8,A)')
	1   '  ERROR: eval residual has',jfnk_resid_nbad,
	2   ' non-finite components -- eval is broken, not lagged.'
	      WRITE(LUER,'(A)')
	1   '  Skipping Newton solve -- falling back to LAMBDA/BA.'
	    ELSE IF(JFNK_COUPLED_LOC .AND. F0_REL .GT. 1.0D-6)THEN
	      JJ_TEST=MAXLOC(ABS(F0_Y-jfnk_fx),1)
	      WRITE(LUER,'(A)')
	1   '    eval vs outer STEQ mismatch  : expected'
	      WRITE(LUER,'(A)')
	1   '      frozen I_NEXT store lag; eval residual adopted'
	2   //' as the Newton baseline'
	      WRITE(LUER,'(A,I10)')
	1   '    worst mismatch eqn index        :',JJ_TEST
	      WRITE(LUER,'(A,I8)')
	1   '      at depth                      :',(JJ_TEST-1)/NT+1
	      WRITE(LUER,'(A,I8)')
	1   '      variable slot                 :',MOD(JJ_TEST-1,NT)+1
	      WRITE(LUER,'(A,ES12.4)')
	1   '      F_eval at that eqn            :',F0_Y(JJ_TEST)
	      WRITE(LUER,'(A,ES12.4)')
	1   '      F_outer at that eqn           :',jfnk_fx(JJ_TEST)
!
! Spread diagnostic: how many equations materially disagree.  A few
! far rows are the store-lag signature at compressing cells; if most
! of the vector moves the eval and the outer sweep are different
! operators -- reported as a warning only, the trust region and the
! committed-step residual checks bound any damage.
!
	      F0_NFAR=0
	      DO JJ_TEST=1,JFNK_N
	        IF(ABS(F0_Y(JJ_TEST)-jfnk_fx(JJ_TEST)) .GT.
	1         0.5D0*ABS(jfnk_fx(JJ_TEST))+1.0D-6*F0_MAX)
	2         F0_NFAR=F0_NFAR+1
	      END DO
	      WRITE(LUER,'(A,I8,A,I8)')
	1   '    eqns differing by more than 50% :',F0_NFAR,' /',JFNK_N
	      IF(F0_NFAR .GT. JFNK_N/2)THEN
	        WRITE(LUER,'(A)')
	1   '    WARNING: global mismatch -- eval may be'
	2   //' inconsistent, not just lagged.'
	      END IF
	      jfnk_fx=F0_Y
	      CALL JFNK_MASK_FIXED_UNKNOWNS(JFNK_N,jfnk_fx)
	      CALL JFNK_COMPUTE_RSCALE_V1(JFNK_N,jfnk_fx)
	    ELSE IF(.NOT. JFNK_COUPLED_LOC .AND. F0_REL .GT. 1.0D-6)THEN
	      F0_OK=.FALSE.
	      WRITE(LUER,'(A)')
	1   '  ERROR: F(0) from JFNK_RESIDUAL_V1 differs from outer STEQ.'
	      WRITE(LUER,'(A,ES10.2)')
	1   '  JFNK residual inconsistent with CMFGEN STEQ (rel=',F0_REL,')'
	      WRITE(LUER,'(A)')
	1   '  Skipping Newton solve -- falling back to LAMBDA/BA.'
	    END IF
	    IF(F0_NAN .GT. 0 .OR. F0_DIFF .NE. F0_DIFF .OR.
	1      F0_REL .NE. F0_REL)THEN
	      F0_OK=.FALSE.
	      WRITE(LUER,'(A)')
	1   '  ERROR: F(0) test produced NaN -- residual is undefined.'
	      WRITE(LUER,'(A)')
	1   '  Skipping Newton solve -- falling back to LAMBDA/BA.'
	    END IF
!
! Adopt the eval residual as the Newton baseline whenever it is
! finite.  The FD JVP subtracts jfnk_fx from F_eval(x+h*v), so a
! baseline of F_outer leaves (F_eval(0)-F_outer)/h inside every
! matvec: a constant vector amplified ~1/h (~7e7) that appeared as
! ~8e3 cross-depth entries in the COLPROBE single-column probes at
! it354.  The evals are bitwise deterministic (repeat-vs-pristine
! |diff|=0), so adopting F_eval(0) removes the offset exactly while
! changing the Newton RHS by at most the 1e-4 eval/outer discrepancy.
!
	    IF(F0_OK)THEN
	      jfnk_fx=F0_Y
	      CALL JFNK_MASK_FIXED_UNKNOWNS(JFNK_N,jfnk_fx)
	      CALL JFNK_COMPUTE_RSCALE_V1(JFNK_N,jfnk_fx)
	    END IF
	    FLUSH(LUER)
	  END BLOCK
	  T_PH_F0=OMP_GET_WTIME()-TW_PH
!
! jfnk_fx now holds F_eval at jfnk_x=0 (adopted in both modes on the
! F0_OK path so the FD JVP baseline is bitwise-consistent with the
! operator being solved).  NEWTON_HOOK_V1 can skip its entry re-eval.
! Only mark it fresh on the F0_OK path -- Newton is skipped otherwise
! and a stale TRUE could leak into the next activation.
!
	  IF(F0_OK)jfnk_fx_fresh=.TRUE.
!
	  JFNK_TOL_ABS=DBLE(JFNK_TOL_REL)*SQRT(JFNK_DOT_V1(JFNK_N,
	1     jfnk_fx,jfnk_fx))
	  IF(JFNK_TOL_ABS .LT. 1.0D-30)JFNK_TOL_ABS=1.0D-30
!
	  F_W_INIT=SQRT(JFNK_DOT_V1(JFNK_N,jfnk_fx,jfnk_fx))
	  F_INF_INIT=0.0D0
	  DO II_F=1,JFNK_N
	    F_INF_INIT=MAX(F_INF_INIT,ABS(jfnk_fx(II_F)))
	  END DO
	  WRITE(LUER,'(/,A)')' ---------------------------------- JFNK SETUP ---------------------------------'
	  WRITE(LUER,'(A,ES12.4)')
	1   '    weighted residual norm at entry   :',REAL(F_W_INIT)
	  WRITE(LUER,'(A,ES12.4)')
	1   '    largest residual component        :',REAL(F_INF_INIT)
	  WRITE(LUER,'(A,ES12.4)')
	1   '    convergence target, weighted norm :',REAL(JFNK_TOL_ABS)
	  WRITE(LUER,'(A,I4)')
	1   '    GMRES restart dimension           :',JFNK_M_GMRES_LOC
	  WRITE(LUER,'(A,I4)')
	1   '    max Newton iterations             :',JFNK_MAX_NEWTON_LOC
	  IF(JFNK_COUPLED_LOC)THEN
	    WRITE(LUER,'(A)')
	1   '    residual evaluation  : coupled, J re-solved per'
	2   //' eval with I_NEXT store held fixed'
	  ELSE
	    WRITE(LUER,'(A)')
	1   '    residual evaluation  : frozen-J, stored rad terms'
	2   //' + CMAT_RAD*x'
	  END IF
	  jfnk_maxch_tol=JFNK_MAXCH_CNV
	  IF(jfnk_log_unknown)THEN
	    WRITE(LUER,'(A)')'    unknown map          : ln(n/n0)'
	    jfnk_step_inf_cap=LOG(MAX(MAX_JFNK_COR,1.01_LDP))
	  ELSE
	    WRITE(LUER,'(A)')'    unknown map          : fractional (1+x)'
	    jfnk_step_inf_cap=MAX(MAX_JFNK_COR-1.0_LDP,0.01_LDP)
	  END IF
	  IF(jfnk_ptc_active)THEN
	    WRITE(LUER,'(A,ES10.3)')
	1   '    pseudo-transient continuation    : on, dtau0=',
	1     REAL(jfnk_ptc_dtau)
	  ELSE
	    WRITE(LUER,'(A)')
	1   '    pseudo-transient continuation    : off'
	  END IF
	  WRITE(LUER,'(A,ES10.3)')
	1   '    Eisenstat-Walker eta floor        :',REAL(jfnk_eta_min)
	  WRITE(LUER,'(A)')' '
!
! FIX_T_AUTO: the domain T-row residual norm RE measures how far the
! current T profile is from NLTE energy balance.  The first nonzero
! reading anchors RE0; rungs release while RE stays below
! JFNK_TAU_REL*RE0 and committed log-T steps stay under JFNK_TAU_EPS.
!
	  IF(JFNK_TAUTO_ACTIVE)THEN
	    TAU_RE_NRM=0.0D0
	    DO K=JFNK_TAU_DOM_ST,JFNK_TAU_DOM_END
	      TAU_RE_NRM=TAU_RE_NRM+
	1       jfnk_fx((K-1)*JFNK_NT+JFNK_NT)**2
	    END DO
	    TAU_RE_NRM=DSQRT(TAU_RE_NRM)
	    IF(JFNK_TAU_RE0 .LE. 0.0_LDP)
	1     JFNK_TAU_RE0=REAL(TAU_RE_NRM,KIND=LDP)
	    JFNK_TAU_RE_LAST=REAL(TAU_RE_NRM,KIND=LDP)
	    WRITE(LUER,'(A,ES10.3,A,ES10.3,A,I2,A,ES10.3,A,I4)')
	1     '    FIX_T_AUTO: |RE|=',REAL(TAU_RE_NRM),
	2     '  |RE0|=',REAL(JFNK_TAU_RE0),
	3     '  rung=',JFNK_TAU_RUNG,
	4     '  alpha=',REAL(JFNK_TAU_ALPHA),
	5     '  commits=',JFNK_TAU_NCOMMIT
	  END IF
	  FLUSH(LUER)
!
	  IF(.NOT. JFNK_BANNER_PRINTED)THEN
	    WRITE(LUER,'(/,A)')' ------------------------------ JFNK CONFIG (once) -----------------------------'
	    IF(JFNK_COUPLED_LOC)THEN
	      WRITE(LUER,'(A)')
	1     '    residual definition : SE at perturbed pops,'
	      WRITE(LUER,'(A)')
	1     '      J re-solved per eval (coupled SE+RTE)'
	    ELSE
	      WRITE(LUER,'(A)')
	1     '    residual definition : STEQ_coll + STEQ_rad +'
	2     //' CMAT_RAD*x (J frozen inside JFNK)'
	    END IF
	    WRITE(LUER,'(A)')
	1   '    preconditioner      : equilibrated depth-block LU of BA'
	    WRITE(LUER,'(A,A)')
	1   '    row residual weighting : ',TRIM(JFNK_RSCALE_MODE)
	    WRITE(LUER,'(A,I3)')
	1   '    max PC refactorizations per activation :',
	1     JFNK_MAX_RESTART
	    WRITE(LUER,'(A)')
	1   '    globalization       : Hookstep trust region +'
	2   //' line search'
	    WRITE(LUER,'(A)')
	1   '    termination         : nonlinear |F|, not GMRES'
	2   //' residual'
	    JFNK_BANNER_PRINTED=.TRUE.
	  END IF
	  FLUSH(LUER)
!
! Newton solve with preconditioner/baseline retry.
!
	  APPLY_JFNK_UPDATE=.FALSE.
	  JFNK_CATASTROPHIC=.FALSE.
	  BASELINE_REBUILT=.FALSE.
	  COMMITTED_BASELINE=.FALSE.
	  IF(F0_OK)THEN
	    IF(jfnk_ptc_active)THEN
	      CALL JFNK_SET_PTC_SCALE_V1(JFNK_N,jfnk_fx)
	      jfnk_ptc_dtau=1.0D0
	    END IF
	    ALLOCATE(BEST_X(JFNK_N))
	    BEST_X=0.0D0
	    BEST_F_W=F_W_INIT
	    BEST_F_INF=F_INF_INIT
	    JFNK_ATTEMPT=0
!
! Reachability precheck (JFNK_REACH_GATE>0, inverse PCs only): the
! preconditioned direction z=M^-1 F is projected back through the BA
! block and REACH=(2<F,Jz>-<Jz,Jz>)/<F,F> measures the fraction of
! |F|^2 a linear step can remove.  On the SOB frozen-J residual the
! radiation-pinned rows leave REACH~0 -- Newton can then only burn
! GMRES/JVP evals producing rejected ~0 steps.  Detecting that here
! costs one apply plus one block matvec instead of tens of SE evals.
! JACOBI is skipped: its diagonal z underestimates reachability and
! would gate out solvable admissions.
! The same block also hosts the direction gates (JFNK_COS_MIN,
! JFNK_JZ_F_MAX, JFNK_CLIP_FRAC_MAX; defaults 0 = off), which reuse
! the same z/Jz metrics plus one FD matvec: a preconditioned step
! that is non-descent, outsized vs |F|, or mostly clip bounds cannot
! produce a useful Newton correction, so the admission defers to the
! caller's SOLVEBA mode before GMRES/JVP evals are spent.  The block
! runs when any of the gates is armed; JFNK_FD_COS_MIN (default
! 0.05) and JFNK_FD_AUTODEF (default T) arm it by themselves so the
! FD descent test is on unless the user opts out.
!
	    JFNK_REACH_OK=.TRUE.
	    TW_PH=OMP_GET_WTIME()
	    IF((JFNK_REACH_GATE .GT. 0.0_LDP .OR.
	1       JFNK_COS_MIN .GT. 0.0_LDP .OR.
	2       JFNK_JZ_F_MAX .GT. 0.0_LDP .OR.
	3       JFNK_CLIP_FRAC_MAX .GT. 0.0_LDP .OR.
	4       JFNK_FD_COS_MIN .GT. 0.0_LDP .OR.
	4       JFNK_FD_AUTODEF .OR.
	5       DO_DIRCHK_NOW) .AND.
	5       JFNK_LU_READY .AND.
	6       (JFNK_PRECOND_TYPE(1:3) .EQ. 'LSQ' .OR.
	7        JFNK_PRECOND_TYPE(1:2) .EQ. 'LU'))THEN
	      ALLOCATE(Z_CHK(JFNK_N), JZ_VEC(JFNK_N), Z_CLIP(JFNK_N))
	      Z_CHK=jfnk_fx
	      CALL JFNK_PRECOND_V1(JFNK_N, Z_CHK)
!
! The commit path clips each component of the accepted step to the
! log-space bound [X_MIN,X_MAX] (population clip = MAX_JFNK_COR); the
! unclipped direction can carry |z|~hundreds on near-null rows that the
! commit would flatten.  Price the step the model would actually take:
! Z_CLIP applies the same bounds here so the stored-BA and FD reach
! numbers describe that step, not the raw M^-1 F.
!
	      X_MAX_R=MAX_JFNK_COR-1.0_LDP
	      X_MIN_R=1.0_LDP/MAX_JFNK_COR-1.0_LDP
	      IF(jfnk_log_unknown)THEN
	        X_MAX_R=LOG(MAX_JFNK_COR)
	        X_MIN_R=-X_MAX_R
	      END IF
	      X_T_MAX_R=LOG(MAX_JFNK_COR)
	      X_T_MIN_R=-X_T_MAX_R
	      IF(MDT_A .GT. 0.0_LDP .AND. MDT_A .LT. 1.0_LDP)THEN
	        X_T_MAX_R=MIN(X_T_MAX_R,LOG(1.0_LDP+MDT_A))
	        X_T_MIN_R=MAX(X_T_MIN_R,LOG(1.0_LDP-MDT_A))
	      END IF
	      Z_CLIP=Z_CHK
	      NCLIP_R=0
	      ZMAX_RAW=0.0_LDP
	      DO I = 1, JFNK_N
	        ZMAX_RAW=MAX(ZMAX_RAW,ABS(Z_CHK(I)))
	        JROW=MOD(I-1,JFNK_NT)+1
	        KROW=(I-1)/JFNK_NT+1
	        IF(JROW .EQ. JFNK_NT)THEN
	          IF(MOD_FIX_T_D_ST .GT. 0 .AND.
	1          KROW .GE. MOD_FIX_T_D_ST .AND.
	2          KROW .LE. MOD_FIX_T_D_END .AND.
	3          .NOT. JFNK_TAUTO_ACTIVE)THEN
	          ELSE
	            IF(Z_CLIP(I) .GT. X_T_MAX_R)THEN
	              Z_CLIP(I)=X_T_MAX_R
	              NCLIP_R=NCLIP_R+1
	            ELSE IF(Z_CLIP(I) .LT. X_T_MIN_R)THEN
	              Z_CLIP(I)=X_T_MIN_R
	              NCLIP_R=NCLIP_R+1
	            END IF
	          END IF
	        ELSE
	          IF(Z_CLIP(I) .GT. X_MAX_R)THEN
	            Z_CLIP(I)=X_MAX_R
	            NCLIP_R=NCLIP_R+1
	          ELSE IF(Z_CLIP(I) .LT. X_MIN_R)THEN
	            Z_CLIP(I)=X_MIN_R
	            NCLIP_R=NCLIP_R+1
	          END IF
	        END IF
	      END DO
	      WRITE(LUER,'(A,ES12.4,A,ES12.4,A,I8)')
	1     '  REACH clip: max|z| ',REAL(ZMAX_RAW),' -> ',
	2     REAL(MAXVAL(ABS(Z_CLIP))),', nclip=',NCLIP_R
	      JZ_VEC=0.0_LDP
	      DO K=1,JFNK_ND
	        DO I=1,JFNK_NT
	          DO J=1,JFNK_NT
	            JZ_VEC((K-1)*JFNK_NT+I)=JZ_VEC((K-1)*JFNK_NT+I)+
	1             JFNK_BA_A(I,J,K)*Z_CLIP((K-1)*JFNK_NT+J)
	          END DO
	        END DO
	      END DO
	      FF2_DP=JFNK_DOT_V1(JFNK_N,jfnk_fx,jfnk_fx)
	      FJZ_DP=JFNK_DOT_V1(JFNK_N,jfnk_fx,JZ_VEC)
	      JJZ_DP=JFNK_DOT_V1(JFNK_N,JZ_VEC,JZ_VEC)
	      REACH_VAL=(2.0_LDP*FJZ_DP-JJZ_DP)/
	1        MAX(FF2_DP,1.0E-300_LDP)
	      WRITE(LUER,'(A,ES12.4,A,ES12.4,A,ES12.4)')
	1        '  JFNK REACH: frac=',REAL(REACH_VAL),
	2        ' |Jz|/|F|=',REAL(SQRT(JJZ_DP)/SQRT(MAX(FF2_DP,1.0E-300_LDP))),
	3        ' cos(F,Jz)=',REAL(FJZ_DP/SQRT(MAX(FF2_DP*JJZ_DP,1.0E-300_LDP)))
	      CALL JFNK_ROW_ANATOMY_V1(LUER,JZ_VEC,JFNK_N)
!
! Direction gates on the stored-BA metrics (defaults 0 = off): these
! cost nothing extra -- z, the clip count, and Jz=BA*z are already
! computed above.  cos(F,Jz)<COS_MIN means even the stored operator
! reports a non-descent direction; |Jz|/|F|>JZ_F_MAX means the linear
! response dwarfs the residual; a clip fraction above CLIP_FRAC_MAX
! means the step the commit would take is mostly the clip bounds, not
! the linearization.  Deferring here avoids the COLPROBE/isolation
! JVP evals below as well as the whole Newton attempt; the admission
! is charged like a reach-gate skip (fallback, not stagnation).
!
	      COS_BA=FJZ_DP/SQRT(MAX(FF2_DP*JJZ_DP,1.0E-300_LDP))
	      JZF_BA=SQRT(JJZ_DP)/SQRT(MAX(FF2_DP,1.0E-300_LDP))
	      CLIP_FR=DBLE(NCLIP_R)/DBLE(JFNK_N)
	      IF(JFNK_COS_MIN .GT. 0.0_LDP .AND.
	1        COS_BA .LT. JFNK_COS_MIN)THEN
	        JFNK_REACH_OK=.FALSE.
	        jfnk_fx_fresh=.FALSE.
	        WRITE(LUER,'(A,ES10.2,A,ES10.2)')
	1        '  JFNK GUARD DEFER: stored-BA non-descent cos(F,Jz)=',
	2        REAL(COS_BA),' < cos_min=',REAL(JFNK_COS_MIN)
	      ELSE IF(JFNK_JZ_F_MAX .GT. 0.0_LDP .AND.
	1        JZF_BA .GT. JFNK_JZ_F_MAX)THEN
	        JFNK_REACH_OK=.FALSE.
	        jfnk_fx_fresh=.FALSE.
	        WRITE(LUER,'(A,ES10.2,A,ES10.2)')
	1        '  JFNK GUARD DEFER: stored-BA |Jz|/|F|=',
	2        REAL(JZF_BA),' > jz_f_max=',REAL(JFNK_JZ_F_MAX)
	      ELSE IF(JFNK_CLIP_FRAC_MAX .GT. 0.0_LDP .AND.
	1        CLIP_FR .GT. JFNK_CLIP_FRAC_MAX)THEN
	        JFNK_REACH_OK=.FALSE.
	        jfnk_fx_fresh=.FALSE.
	        WRITE(LUER,'(A,ES10.2,A,ES10.2,A,I8)')
	1        '  JFNK GUARD DEFER: clip fraction=',
	2        REAL(CLIP_FR),' > clip_frac_max=',REAL(JFNK_CLIP_FRAC_MAX),
	3        ' nclip=',NCLIP_R
	      END IF
!
! Pristine-vs-repeat column probe (11-Oct-2026): the residual is
! formally block-separable (no routine reads pops at other depths),
! yet a single-column FD probe produced nonzero entries at distant
! depths when measured AFTER a full-vector JVP had perturbed every
! depth.  If module state leaks across evals, later JVPs inherit
! residue from earlier perturbations and the FD operator looks
! dense.  Measure the top-|z| column FIRST (before any full-vector
! JVP in this block), then repeat it at the end: equal outputs mean
! pure evals and real global coupling; different outputs prove
! eval-order contamination.
!
	    IF(JFNK_REACH_OK)THEN
	      ALLOCATE(JV_FD(JFNK_N), JV_FD2(JFNK_N), Z_BLK(JFNK_N))
	      ALLOCATE(JV_SAV(JFNK_N))
	      IF(ALLOCATED(J2ID_LOC))DEALLOCATE(J2ID_LOC)
	      ALLOCATE(J2ID_LOC(JFNK_NT))
	      J2ID_LOC=0
	      DO ID=1,JFNK_NION
	        IF(ATM(ID)%XzV_PRES)THEN
	          DO J=ATM(ID)%EQXzV,
	1            MIN(ATM(ID)%EQXzV+ATM(ID)%NXzV-1,JFNK_NT-2)
	            IF(J .GE. 1)J2ID_LOC(J)=ID
	          END DO
	        ELSE IF(ATM(ID)%EQXzV .GE. 1 .AND.
	1            ATM(ID)%EQXzV .LE. JFNK_NT-2)THEN
	          J2ID_LOC(ATM(ID)%EQXzV)=ID
	        END IF
	      END DO
	      FMAX_LOC = 0.0_LDP
	      I_ZMAX = 1
	      DO I = 1, JFNK_N
	        IF(ABS(Z_CHK(I)) .GT. FMAX_LOC)THEN
	          FMAX_LOC = ABS(Z_CHK(I))
	          I_ZMAX = I
	        END IF
	      END DO
	      Z_BLK = 0.0_LDP
	      Z_BLK(I_ZMAX) = 1.0_LDP
	      CALL JFNK_JVP_V1(JFNK_N, Z_BLK, JV_SAV)
	      WRITE(LUER,'(A,I8,A,ES12.4)')
	1        '  COLPROBE pristine topZ idx=',I_ZMAX,
	2        ' |FDcol|=',
	3        REAL(SQRT(JFNK_DOT_V1(JFNK_N,JV_SAV,JV_SAV)))
!
! FD-vs-stored operator check (11-Oct-2026): the stored-BA REACH can
! report frac~1 while the finite-difference operator maps the same
! z orthogonally (LU admission at it297: REACH=0.99998, GMRES flat).
! One extra residual eval prices A_FD*z directly; the anatomy names
! the rows where stored BA and true dF/dx diverge.
!
	      CALL JFNK_JVP_V1(JFNK_N, Z_CLIP, JV_FD)
	      FJZ_FD=JFNK_DOT_V1(JFNK_N,jfnk_fx,JV_FD)
	      JJZ_FD=JFNK_DOT_V1(JFNK_N,JV_FD,JV_FD)
	      BAFD_DP=JFNK_DOT_V1(JFNK_N,JZ_VEC,JV_FD)
	      WRITE(LUER,'(A,ES12.4,A,ES12.4,A,ES12.4)')
	1        '  JFNK REACH-FD (clipped): |A_FD z|/|F|=',
	2        REAL(SQRT(JJZ_FD)/SQRT(MAX(FF2_DP,1.0E-300_LDP))),
	3        ' cos(F,A_FDz)=',
	4        REAL(FJZ_FD/SQRT(MAX(FF2_DP*JJZ_FD,1.0E-300_LDP))),
	5        ' cos(BAz,A_FDz)=',
	6        REAL(BAFD_DP/SQRT(MAX(JJZ_DP*JJZ_FD,1.0E-300_LDP)))
	      CALL JFNK_ROW_ANATOMY_V1(LUER,JV_FD,JFNK_N)
!
! Direction gates on the FD image (defaults 0 = off): the stored-BA
! stage-A tests can pass on states where the true Jacobian maps z
! almost orthogonally to F (the stored operator and dF/dx disagree).
! cos(F,A_FDz) and |A_FDz|/|F| are the same two metrics on the FD
! image of the clipped step -- the direction the nonlinear residual
! actually sees.  Deferring here still costs only the evals already
! spent above; the remaining isolation probes are skipped.
!
	      COS_FD=FJZ_FD/SQRT(MAX(FF2_DP*JJZ_FD,1.0E-300_LDP))
	      JZF_FD=SQRT(JJZ_FD)/SQRT(MAX(FF2_DP,1.0E-300_LDP))
	      IF(JFNK_COS_MIN .GT. 0.0_LDP .AND.
	1        COS_FD .LT. JFNK_COS_MIN)THEN
	        JFNK_REACH_OK=.FALSE.
	        jfnk_fx_fresh=.FALSE.
	        WRITE(LUER,'(A,ES10.2,A,ES10.2)')
	1        '  JFNK GUARD DEFER: FD non-descent cos(F,Az)=',
	2        REAL(COS_FD),' < cos_min=',REAL(JFNK_COS_MIN)
	      ELSE IF(JFNK_JZ_F_MAX .GT. 0.0_LDP .AND.
	1        JZF_FD .GT. JFNK_JZ_F_MAX)THEN
	        JFNK_REACH_OK=.FALSE.
	        jfnk_fx_fresh=.FALSE.
	        WRITE(LUER,'(A,ES10.2,A,ES10.2)')
	1        '  JFNK GUARD DEFER: FD |Az|/|F|=',
	2        REAL(JZF_FD),' > jz_f_max=',REAL(JFNK_JZ_F_MAX)
	      END IF
!
! FD descent floor (17-Oct-2026, [JFNK_FD_COS_MIN] default 0.05 +
! [JFNK_FD_AUTODEF] default T): the FD image is the descent test
! against the TRUE operator -- the stored-BA REACH passed a step
! whose real Jacobian maps it anti-parallel to F (KNOP_jfnkrel
! it35, cos=-0.23) and GMRES then spent 30 evals stagnating.
! Deferring here costs only the evals already spent; GMRES is
! skipped.  Independent of JFNK_COS_MIN so the stored-BA policy is
! untouched; set JFNK_FD_COS_MIN=0 and JFNK_FD_AUTODEF=F to recover
! the pre-guard behaviour.
!
	      IF(JFNK_REACH_OK .AND.
	1      ((JFNK_FD_COS_MIN .GT. 0.0_LDP .AND.
	2        COS_FD .LT. JFNK_FD_COS_MIN) .OR.
	3       (JFNK_FD_AUTODEF .AND. COS_FD .LT. 0.0_LDP)))THEN
	        JFNK_REACH_OK=.FALSE.
	        jfnk_fx_fresh=.FALSE.
	        WRITE(LUER,'(A,ES10.2,A,ES10.2)')
	1        '  JFNK GUARD DEFER: FD non-descent cos(F,Az)=',
	2        REAL(COS_FD),' < fd_cos_min=',REAL(JFNK_FD_COS_MIN)
	      END IF
!
	    IF(JFNK_REACH_OK)THEN
! Isolation probes (11-Oct-2026): A_FD*z returned ~60x|F| with
! cos(BAz,A_FDz)~0.2 -- the stored BA misses large true-J action
! excited by z's ~50-magnitude components.  JVP step h scales with
! 1/||v||_inf so scaling v does NOT shrink the perturbation; the
! informative restrictions are structural:
!   (a) z|_(K_TOPF-1): response on the top-F row (K_TOPF) prices the
!       cross-depth derivative the DIAG BA cannot store;
!   (b) z with |z_i|>1 zeroed: does the mystery response vanish when
!       the soft huge components are removed (they sit in stored-BA's
!       near-null space but may be non-null under true J);
!   (c) print top-|z| rows so the offending components are named.
!
	      FMAX_LOC = 0.0_LDP
	      I_TOPF = 1
	      DO I = 1, JFNK_N
	        IF(ABS(jfnk_fx(I)) .GT. FMAX_LOC)THEN
	          FMAX_LOC = ABS(jfnk_fx(I))
	          I_TOPF = I
	        END IF
	      END DO
	      K_TOPF = (I_TOPF-1)/JFNK_NT + 1
!
! Top-|z| components: name the rows carrying the large corrections.
!
	      WRITE(LUER,'(A)')'  REACH-FD top |z| rows:'
	      WRITE(LUER,'(A)')'    rnk    idx     K     J  z'//
	1      '            ion        lev   n'//
	2      '          n/nLTE      F_i        rad_diag    BA_diag'
	      Z_BLK = 0.0_LDP
	      DO K = 1, 50
	        FMAX_LOC = 0.0_LDP
	        I_ZMAX = 0
	        DO I = 1, JFNK_N
	          IF(ABS(Z_CHK(I)) .GT. FMAX_LOC .AND.
	1             Z_BLK(I) .EQ. 0.0_LDP)THEN
	            FMAX_LOC = ABS(Z_CHK(I))
	            I_ZMAX = I
	          END IF
	        END DO
	        IF(I_ZMAX .EQ. 0)EXIT
	        Z_BLK(I_ZMAX) = 1.0_LDP
	        JROW=MOD(I_ZMAX-1,JFNK_NT)+1
	        KROW=(I_ZMAX-1)/JFNK_NT+1
	        ID=0
	        IF(JROW .LE. JFNK_NT-2)ID=J2ID_LOC(JROW)
	        LBL_LOC='unlabeled'
	        LEV_LOC=0
	        POP_LOC=0.0_LDP
	        NLTE_LOC=0.0_LDP
	        IF(JROW .EQ. JFNK_NT)THEN
	          LBL_LOC='T (energy)'
	          POP_LOC=POPS(JROW,KROW)
	        ELSE IF(JROW .EQ. JFNK_NT-1)THEN
	          LBL_LOC='Ne (charge)'
	          POP_LOC=POPS(JROW,KROW)
	        ELSE IF(ID .GT. 0)THEN
	          LBL_LOC=ION_ID(ID)
	          LEV_LOC=JROW-ATM(ID)%EQXzV+1
	          POP_LOC=POPS(JROW,KROW)
	          IF(ALLOCATED(ATM(ID)%XzVLTE) .AND.
	1            ABS(ATM(ID)%XzVLTE(LEV_LOC,KROW)) .GT. 0.0_LDP)
	2            NLTE_LOC=POP_LOC/ATM(ID)%XzVLTE(LEV_LOC,KROW)
	        END IF
	        CMR_LOC=0.0_LDP
	        IF(ALLOCATED(JFNK_CMAT_RAD))
	1         CMR_LOC=DBLE(JFNK_CMAT_RAD(JROW,JROW,KROW))
	        BADG_LOC=JFNK_BA_A(JROW,JROW,KROW)
	        WRITE(LUER,'(I4,I9,I6,I6,ES12.4,2X,A12,I5,4ES12.4)')
	1          K,I_ZMAX,KROW,JROW,Z_CHK(I_ZMAX),LBL_LOC,LEV_LOC,
	2          POP_LOC,NLTE_LOC,jfnk_fx(I_ZMAX),CMR_LOC,BADG_LOC
	      END DO
	      IF(ALLOCATED(NZ_SP))DEALLOCATE(NZ_SP)
	      ALLOCATE(NZ_SP(JFNK_NION))
	      NZ_SP=0
	      DO I = 1, JFNK_N
	        JROW=MOD(I-1,JFNK_NT)+1
	        IF(JROW .LE. JFNK_NT-2 .AND.
	1          ABS(Z_CHK(I)) .GT. 1.0_LDP)THEN
	          ID=J2ID_LOC(JROW)
	          IF(ID .GT. 0)NZ_SP(ID)=NZ_SP(ID)+1
	        END IF
	      END DO
	      WRITE(LUER,'(A)')'    |z|>1 census by species (all depths):'
	      DO ID=1,JFNK_NION
	        IF(NZ_SP(ID) .GT. 0)
	1        WRITE(LUER,'(A,A12,I6)')'      ',ION_ID(ID),NZ_SP(ID)
	      END DO
!
! (a) column probe: JVP(e_i) returns column i of the true J
! (block-separable residual -- K69-only restriction gave zero response
! at the K70 top row, so all action is within-block).  Compare
! elementwise against the stored BA column for the top-F row and the
! top-|z| row: names the entries where stored BA and true dF/dx
! differ.  The BA column lives in block K_TOPF of JFNK_BA_A.
!
	      Z_BLK = 0.0_LDP
	      Z_BLK(I_TOPF) = 1.0_LDP
	      CALL JFNK_JVP_V1(JFNK_N, Z_BLK, JV_FD2)
	      JV_FD = 0.0_LDP
	      DO I = (K_TOPF-1)*JFNK_NT+1, K_TOPF*JFNK_NT
	        JV_FD(I) = JFNK_BA_A(MOD(I-1,JFNK_NT)+1,
	1             MOD(I_TOPF-1,JFNK_NT)+1, K_TOPF)
	      END DO
	      CALL JFNK_MASK_FIXED_UNKNOWNS(JFNK_N, JV_FD)
	      FJZ_FD = JFNK_DOT_V1(JFNK_N,JV_FD2,JV_FD2)
	      JJZ_FD = JFNK_DOT_V1(JFNK_N,JV_FD,JV_FD)
	      BAFD_DP = JFNK_DOT_V1(JFNK_N,JV_FD,JV_FD2)
	      WRITE(LUER,'(A,I8,A,ES12.4,A,ES12.4,A,ES12.4)')
	1        '  COLPROBE topF idx=',I_TOPF,' |FDcol|=',
	2        REAL(SQRT(FJZ_FD)),
	3        ' |BAcol|=',REAL(SQRT(JJZ_FD)),
	4        ' |diff|/|BA|=',
	5        REAL(SQRT(MAX(FJZ_FD+JJZ_FD-2.0_LDP*BAFD_DP,
	1             0.0_LDP))/MAX(SQRT(JJZ_FD),1.0E-300_LDP))
!
! J->ion map for the per-element labels below; same layout as
! JFNK_ROW_ANATOMY_V1 (levels EQXzV..EQXzV+NXzV-1 belong to ion ID;
! J=NT-1 is the charge row, J=NT the energy row).
!
	      IF(ALLOCATED(J2ID_LOC))DEALLOCATE(J2ID_LOC)
	      ALLOCATE(J2ID_LOC(JFNK_NT))
	      J2ID_LOC=0
	      DO ID=1,JFNK_NION
	        IF(ATM(ID)%XzV_PRES)THEN
	          DO J=ATM(ID)%EQXzV,
	1            MIN(ATM(ID)%EQXzV+ATM(ID)%NXzV-1,JFNK_NT-2)
	            IF(J .GE. 1)J2ID_LOC(J)=ID
	          END DO
	        ELSE IF(ATM(ID)%EQXzV .GE. 1 .AND.
	1            ATM(ID)%EQXzV .LE. JFNK_NT-2)THEN
	          J2ID_LOC(ATM(ID)%EQXzV)=ID
	        END IF
	      END DO
	      DO K = 1, 8
	        FMAX_LOC = 0.0_LDP
	        I_ZMAX = 0
	        DO I = 1, JFNK_N
	          IF(ABS(JV_FD2(I)-JV_FD(I)) .GT. FMAX_LOC)THEN
	            FMAX_LOC = ABS(JV_FD2(I)-JV_FD(I))
	            I_ZMAX = I
	          END IF
	        END DO
	        IF(I_ZMAX .EQ. 0)EXIT
	        WRITE(LUER,'(A,I8,A,I4,A,I6,A,ES12.4,A,ES12.4)')
	1          '    col: idx=',I_ZMAX,' K=',(I_ZMAX-1)/JFNK_NT+1,
	2          ' J=',MOD(I_ZMAX-1,JFNK_NT)+1,
	3          ' FD=',JV_FD2(I_ZMAX),' BA=',JV_FD(I_ZMAX)
!
! Row decomposition: F_i (a self-proportional equation gives
! FD ~ -F_i rather than the raw BA outflow), the CMAT_RAD share of
! the BA element (valid only for rows in the probed depth block),
! ion/level label, and the frozen REPLACE flag of the owning ion.
!
	        JROW=MOD(I_ZMAX-1,JFNK_NT)+1
	        KROW=(I_ZMAX-1)/JFNK_NT+1
	        ID=0
	        IF(JROW .LE. JFNK_NT-2)ID=J2ID_LOC(JROW)
	        CMR_LOC=0.0_LDP
	        IF(ALLOCATED(JFNK_CMAT_RAD) .AND. KROW .EQ. K_TOPF)
	1           CMR_LOC=DBLE(JFNK_CMAT_RAD(JROW,
	2             MOD(I_TOPF-1,JFNK_NT)+1,KROW))
	        COLL_LOC=JV_FD(I_ZMAX)-CMR_LOC
	        LBL_LOC='unlabeled'
	        LEV_LOC=0
	        REP_LOC=.FALSE.
	        IF(JROW .EQ. JFNK_NT)THEN
	          LBL_LOC='T (energy)'
	        ELSE IF(JROW .EQ. JFNK_NT-1)THEN
	          LBL_LOC='Ne (charge)'
	        ELSE IF(ID .GT. 0)THEN
	          LBL_LOC=ION_ID(ID)
	          LEV_LOC=JROW-ATM(ID)%EQXzV+1
	          IF(ALLOCATED(JFNK_REPLACE_SAV))
	1             REP_LOC=JFNK_REPLACE_SAV(ID,KROW)
	        END IF
	        WRITE(LUER,'(3(A,ES12.4),2X,A12,I5,A,L1)')
	1          '      F_i=',jfnk_fx(I_ZMAX),' rad=',CMR_LOC,
	2          ' coll=',COLL_LOC,LBL_LOC,LEV_LOC,' rep=',REP_LOC
	        JV_FD(I_ZMAX) = JV_FD2(I_ZMAX)
	      END DO
	!
	! h-scan on the topF self element (Issue 6): the topF column shows
	! sign-flipped off-diagonal FD-vs-BA entries (HI J=18-21 at K=70),
	! which can be genuine w-dissolution physics the BA froze out or FD
	! cancellation noise.  Stability across decades of h separates them.
	!
	      WRITE(LUER,'(A)')'  COLPROBE hscan (topF self):'
	      DO IHS=1,3
	        JFNK_PROBE_HSCALE=HS_ARR(IHS)
	        Z_BLK=0.0_LDP
	        Z_BLK(I_TOPF)=1.0_LDP
	        CALL JFNK_JVP_V1(JFNK_N,Z_BLK,JV_FD2)
	        WRITE(LUER,'(A,ES10.2,A,ES14.6)')
	1          '    h x',HS_ARR(IHS),' FD=',JV_FD2(I_TOPF)
	      END DO
	      JFNK_PROBE_HSCALE=1.0_LDP
!
! (b) same column probe for the top-|z| row (the ~147 component).
!
	      FMAX_LOC = 0.0_LDP
	      I_ZMAX = 1
	      DO I = 1, JFNK_N
	        IF(ABS(Z_CHK(I)) .GT. FMAX_LOC)THEN
	          FMAX_LOC = ABS(Z_CHK(I))
	          I_ZMAX = I
	        END IF
	      END DO
	      J = (I_ZMAX-1)/JFNK_NT + 1
	      Z_BLK = 0.0_LDP
	      Z_BLK(I_ZMAX) = 1.0_LDP
	      CALL JFNK_JVP_V1(JFNK_N, Z_BLK, JV_FD2)
	      JV_FD = 0.0_LDP
	      DO I = (J-1)*JFNK_NT+1, J*JFNK_NT
	        JV_FD(I) = JFNK_BA_A(MOD(I-1,JFNK_NT)+1,
	1             MOD(I_ZMAX-1,JFNK_NT)+1, J)
	      END DO
	      CALL JFNK_MASK_FIXED_UNKNOWNS(JFNK_N, JV_FD)
	      FJZ_FD = JFNK_DOT_V1(JFNK_N,JV_FD2,JV_FD2)
	      JJZ_FD = JFNK_DOT_V1(JFNK_N,JV_FD,JV_FD)
	      BAFD_DP = JFNK_DOT_V1(JFNK_N,JV_FD,JV_FD2)
	      WRITE(LUER,'(A,I8,A,ES12.4,A,ES12.4,A,ES12.4)')
	1        '  COLPROBE topZ idx=',I_ZMAX,' |FDcol|=',
	2        REAL(SQRT(FJZ_FD)),
	3        ' |BAcol|=',REAL(SQRT(JJZ_FD)),
	4        ' |diff|/|BA|=',
	5        REAL(SQRT(MAX(FJZ_FD+JJZ_FD-2.0_LDP*BAFD_DP,
	1             0.0_LDP))/MAX(SQRT(JJZ_FD),1.0E-300_LDP))
	      DO K = 1, 8
	        FMAX_LOC = 0.0_LDP
	        I_TOPF = 0
	        DO II_PROBE = 1, JFNK_N
	          IF(ABS(JV_FD2(II_PROBE)-JV_FD(II_PROBE)) .GT. FMAX_LOC)THEN
	            FMAX_LOC = ABS(JV_FD2(II_PROBE)-JV_FD(II_PROBE))
	            I_TOPF = II_PROBE
	          END IF
	        END DO
	        IF(I_TOPF .EQ. 0)EXIT
	        WRITE(LUER,'(A,I8,A,I4,A,I6,A,ES12.4,A,ES12.4)')
	1          '    col: idx=',I_TOPF,' K=',(I_TOPF-1)/JFNK_NT+1,
	2          ' J=',MOD(I_TOPF-1,JFNK_NT)+1,
	3          ' FD=',JV_FD2(I_TOPF),' BA=',JV_FD(I_TOPF)
!
! Same row decomposition as the topF block.  Here J holds the
! probed depth and MOD(I_ZMAX-1,NT)+1 the probed column (I_ZMAX is
! the top-|z| index selected above the loop).
!
	        JROW=MOD(I_TOPF-1,JFNK_NT)+1
	        KROW=(I_TOPF-1)/JFNK_NT+1
	        ID=0
	        IF(JROW .LE. JFNK_NT-2)ID=J2ID_LOC(JROW)
	        CMR_LOC=0.0_LDP
	        IF(ALLOCATED(JFNK_CMAT_RAD) .AND. KROW .EQ. J)
	1           CMR_LOC=DBLE(JFNK_CMAT_RAD(JROW,
	2             MOD(I_ZMAX-1,JFNK_NT)+1,KROW))
	        COLL_LOC=JV_FD(I_TOPF)-CMR_LOC
	        LBL_LOC='unlabeled'
	        LEV_LOC=0
	        REP_LOC=.FALSE.
	        IF(JROW .EQ. JFNK_NT)THEN
	          LBL_LOC='T (energy)'
	        ELSE IF(JROW .EQ. JFNK_NT-1)THEN
	          LBL_LOC='Ne (charge)'
	        ELSE IF(ID .GT. 0)THEN
	          LBL_LOC=ION_ID(ID)
	          LEV_LOC=JROW-ATM(ID)%EQXzV+1
	          IF(ALLOCATED(JFNK_REPLACE_SAV))
	1             REP_LOC=JFNK_REPLACE_SAV(ID,KROW)
	        END IF
	        WRITE(LUER,'(3(A,ES12.4),2X,A12,I5,A,L1)')
	1          '      F_i=',jfnk_fx(I_TOPF),' rad=',CMR_LOC,
	2          ' coll=',COLL_LOC,LBL_LOC,LEV_LOC,' rep=',REP_LOC
	        JV_FD(I_TOPF) = JV_FD2(I_TOPF)
	      END DO
!
! Pristine-vs-repeat: JV_SAV holds the topZ column measured before any
! full-vector JVP in this block; JV_FD2 holds the same column measured
! after several full-vector evals.  |repeat-pristine|/|pristine| near
! zero means the evals are pure and the cross-depth entries are real
! Jacobian structure; a large value proves eval-order contamination.
!
	      JV_FD = JV_FD2 - JV_SAV
	      FJZ_FD = SQRT(JFNK_DOT_V1(JFNK_N,JV_FD,JV_FD))
	      JVP_N = SQRT(JFNK_DOT_V1(JFNK_N,JV_SAV,JV_SAV))
	      WRITE(LUER,'(A,ES12.4,A,ES12.4)')
	1        '  COLPROBE repeat-vs-pristine: |diff|=',
	2        REAL(FJZ_FD),
	3        ' rel=',REAL(FJZ_FD/MAX(JVP_N,1.0E-300_LDP))
!
! h-scan on the topZ self element.  A derivative stable across
! decades of step size is real; one that drifts or flips is
! cancellation noise (row partial terms reach 1D15-1D19 against a
! nominal 1.5D-8 perturbation).
!
	      WRITE(LUER,'(A)')'  COLPROBE hscan (topZ self):'
	      DO IHS=1,3
	        JFNK_PROBE_HSCALE=HS_ARR(IHS)
	        Z_BLK=0.0_LDP
	        Z_BLK(I_ZMAX)=1.0_LDP
	        CALL JFNK_JVP_V1(JFNK_N,Z_BLK,JV_FD2)
	        WRITE(LUER,'(A,ES10.2,A,ES14.6)')
	1          '    h x',HS_ARR(IHS),' FD=',JV_FD2(I_ZMAX)
	      END DO
	      JFNK_PROBE_HSCALE=1.0_LDP
!
! (b) soft-component removal: zero every |z_i|>1 and re-price.
!
	      Z_BLK = Z_CHK
	      DO I = 1, JFNK_N
	        IF(ABS(Z_BLK(I)) .GT. 1.0_LDP) Z_BLK(I) = 0.0_LDP
	      END DO
	      CALL JFNK_JVP_V1(JFNK_N, Z_BLK, JV_FD2)
	      FJZ_FD=JFNK_DOT_V1(JFNK_N,jfnk_fx,JV_FD2)
	      JJZ_FD=JFNK_DOT_V1(JFNK_N,JV_FD2,JV_FD2)
	      BAFD_DP=JFNK_DOT_V1(JFNK_N,JZ_VEC,JV_FD2)
	      WRITE(LUER,'(A,ES12.4,A,ES12.4,A,ES12.4)')
	1        '  JFNK REACH-FD (|z|>1 zeroed): |A_FD z|/|F|=',
	2        REAL(SQRT(JJZ_FD)/SQRT(MAX(FF2_DP,1.0E-300_LDP))),
	3        ' cos(F,A_FDz)=',
	4        REAL(FJZ_FD/SQRT(MAX(FF2_DP*JJZ_FD,1.0E-300_LDP))),
	5        ' cos(BAz,A_FDz)=',
	6        REAL(BAFD_DP/SQRT(MAX(JJZ_DP*JJZ_FD,1.0E-300_LDP)))
	      CALL JFNK_ROW_ANATOMY_V1(LUER,JV_FD2,JFNK_N)
	    END IF
	      DEALLOCATE(JV_FD, JV_FD2, Z_BLK, JV_SAV)
	    END IF
! The REACH metric is priced against JFNK_BA_A, the outer
! (admission-time) BA linearization -- valid for the frozen tier,
! whose residual is that same fixed-J system.  The coupled tier
! re-solves J per eval, so its Newton differentiates a different
! operator; a stored-BA REACH would be the wrong test there and can
! veto the only tier able to reduce a J-pinned residual.  Gate the
! frozen tier only; the GMRES/Hookstep budget judges coupled
! reachability.  The COS/JZ_F/CLIP direction guards above share the
! same stored-BA basis and have the same limitation -- they are
! default-off diagnostics and are left unchanged.
!
	      IF(JFNK_REACH_GATE .GT. 0.0_LDP .AND.
	1        .NOT. JFNK_COUPLED_LOC .AND.
	2        REACH_VAL .LT. JFNK_REACH_GATE)THEN
	        JFNK_REACH_OK=.FALSE.
	        JFNK_REACH_SKIP=.TRUE.
	        WRITE(LUER,'(A,ES10.2,A)')
	1        '  JFNK: residual unreachable by population step',
	2        REAL(JFNK_REACH_GATE),' gate -- skipping Newton'
	        WRITE(LUER,'(A)')
	1        '  (frozen-J rows are radiation-pinned; LAMBDA/BA step'
	2        //' advances the outer lag instead)'
	        jfnk_fx_fresh=.FALSE.
	      END IF
	      DEALLOCATE(Z_CHK, JZ_VEC, Z_CLIP)
	    ELSE
	      CALL JFNK_ROW_ANATOMY_V1(LUER,jfnk_fx,0)
	    END IF
!
! DO_BA_DIR_CHK (VADAT [DO_BA_DIR_CHK], default F): diagnostics only.
! The precheck block above has priced the stored-BA direction the
! caller-mode SOLVEBA step is about to take (cos(F,Jz), |Jz|/|F|,
! REACH, plus the FD-vs-stored comparison when JFNK_REACH_OK stayed
! true).  Force the reach-fail path so Newton never runs: info=-2
! exits through the admission-gate branch, which is not counted as
! stagnation, and the same call's SOLVEBA_V13 fallback applies the
! real correction.  jfnk_fx_fresh must clear -- the adopted F_eval
! baseline belongs to the skipped Newton, not to a later activation.
!
	    IF(DO_DIRCHK_NOW)THEN
	      JFNK_DIRCHK_N=JFNK_DIRCHK_N+1
	      JFNK_REACH_OK=.FALSE.
	      jfnk_fx_fresh=.FALSE.
	      WRITE(LUER,'(A,I4,A)')
	1     '  JFNK DIRCHK: diagnostics complete; deferring Newton'//
	2     ' to caller-mode SOLVEBA (event ',JFNK_DIRCHK_N,')'
	    END IF
	    IF(.NOT. JFNK_REACH_OK)JFNK_INFO=-2
	    T_PH_PREC=OMP_GET_WTIME()-TW_PH
	    TW_PH=OMP_GET_WTIME()
	    FLUSH(LUER)
	    DO WHILE(JFNK_REACH_OK .AND. JFNK_ATTEMPT .LT. JFNK_MAX_RESTART)
	      JFNK_ATTEMPT=JFNK_ATTEMPT+1
	      IF(JFNK_ATTEMPT .GT. 1)THEN
	        IF(BEST_F_W .LT. F_W_INIT)THEN
	          WRITE(LUER,'(A,I3,A,I3,A,ES12.4)')
	1           '  JFNK: REBUILDING BASELINE (attempt ',
	2           JFNK_ATTEMPT,'/',JFNK_MAX_RESTART,
	3           ') best|F|=',REAL(BEST_F_W)
	          CALL JFNK_REBUILD_BASELINE_V1(BEST_X,LUER)
	          BASELINE_REBUILT=.TRUE.
	          F_W_INIT=SQRT(JFNK_DOT_V1(JFNK_N,jfnk_fx,jfnk_fx))
	          F_INF_INIT=0.0D0
	          DO II_F=1,JFNK_N
	            F_INF_INIT=MAX(F_INF_INIT,ABS(jfnk_fx(II_F)))
	          END DO
	          BEST_F_W=F_W_INIT
	          BEST_F_INF=F_INF_INIT
	          BEST_X=0.0D0
	          FLUSH(LUER)
	        ELSE
	          WRITE(LUER,'(A,I3,A,I3)')
	1           '  JFNK: REBUILDING PRECONDITIONER (attempt ',
	2           JFNK_ATTEMPT,'/',JFNK_MAX_RESTART,')'
	          CALL JFNK_FACTOR_PRECOND_V1()
!
! The rebuilt BA was assembled at the current state; let the next
! eval re-harvest its equation structure so the residual and the
! preconditioner stay on the same REPLACE/ZERO_STEQ set.
!
	          JFNK_EQ_STRUCT_READY=.FALSE.
	          FLUSH(LUER)
	        END IF
	      END IF
	      JFNK_INFO=1
	      FLUSH(LUER)
	      CALL NEWTON_HOOK_V1(
	1          JFNK_RESIDUAL_V1,
	2          JFNK_JVP_V1,
	3          JFNK_PRECOND_V1,
	4          JFNK_DIAGNOSTIC_V1,
	5          JFNK_DOT_V1,
	6          JFNK_M_GMRES_LOC,
	7          JFNK_N,
	8          DBLE(JFNK_GTOL_LOC),
	9          JFNK_TOL_ABS,
	1          DBLE(JFNK_DEL_INIT_LOC),
	1          DBLE(JFNK_MNDL_LOC),
	2          DBLE(JFNK_MXDL_LOC),
	3          JFNK_MAX_NEWTON_LOC,
	4          JFNK_INFO,
	5          LUER)
	      FLUSH(LUER)
	      F_W=SQRT(JFNK_DOT_V1(JFNK_N,jfnk_fx,jfnk_fx))
	      F_INF=0.0D0
	      DO II_F=1,JFNK_N
	        F_INF=MAX(F_INF,ABS(jfnk_fx(II_F)))
	        IF(jfnk_fx(II_F) .NE. jfnk_fx(II_F))JFNK_CATASTROPHIC=.TRUE.
	        IF(jfnk_x(II_F) .NE. jfnk_x(II_F))JFNK_CATASTROPHIC=.TRUE.
	      END DO
	      WRITE(LUER,'(/,A)')' --------------------------------- JFNK RESULT ---------------------------------'
	      WRITE(LUER,'(A,ES12.4)')
	1       '    weighted residual norm after solve :',REAL(F_W)
	      WRITE(LUER,'(A,ES12.4)')
	1       '    est max relative pop correction    :',
	2       REAL(jfnk_maxch_est)
	      WRITE(LUER,'(A,ES12.4)')
	1       '    largest residual component         :',REAL(F_INF)
	      WRITE(LUER,'(A,ES12.4)')
	1       '    convergence target, weighted norm  :',
	1       REAL(JFNK_TOL_ABS)
	      WRITE(LUER,'(A,I3)')
	1       '    solver info code                   :',JFNK_INFO
	      WRITE(LUER,'(2A)')
	1       '    solver exit reason                 : ',
	2       TRIM(jfnk_exit_reason)
	      IF(.NOT. JFNK_CATASTROPHIC .AND. F_W .LT. BEST_F_W)THEN
	        BEST_F_W=F_W
	        BEST_F_INF=F_INF
	        BEST_X=jfnk_x
	      END IF
	      IF(JFNK_CATASTROPHIC .OR. JFNK_INFO .EQ. 4 .OR.
	1         JFNK_INFO .EQ. 3)THEN
! editted by Devin, 04-Oct-2026: info=3 (Krylov/JVP/preconditioner
! abort, or trust radius exhausted) is deterministic for this
! baseline -- the retry would re-factor the same BA at the same
! point and replay the identical failure, so leave the restart
! loop immediately instead of burning JFNK_MAX_RESTART attempts.
	        WRITE(LUER,'(2A,I3)')
	1         '    solver status                    : ',
	2         'unrecoverable, info=',JFNK_INFO
	        EXIT
	      END IF
	      IF(JFNK_INFO .EQ. 0)THEN
	        WRITE(LUER,'(A)')
	1         '    solver status                    : CONVERGED'
	        APPLY_JFNK_UPDATE=.TRUE.
	        EXIT
	      END IF
	      IF(JFNK_INFO .EQ. 2)THEN
	        WRITE(LUER,'(A)')
	1         '    solver status                    : max'
	2         //' Newton iterations reached'
	      END IF
	    END DO
	    T_PH_NWT=OMP_GET_WTIME()-TW_PH
	    IF(.NOT. APPLY_JFNK_UPDATE .AND. .NOT. JFNK_CATASTROPHIC .AND.
	1     BEST_F_W .LT. F_W_INIT)THEN
	      jfnk_x=BEST_X
	      APPLY_JFNK_UPDATE=.TRUE.
	      WRITE(LUER,'(A,ES12.4)')
	1       '    best-iterate weighted residual     :',
	1       REAL(BEST_F_W)
	      WRITE(LUER,'(A,ES12.4)')
	1       '    weighted residual at entry         :',
	1       REAL(F_W_INIT)
	      WRITE(LUER,'(A)')
	1       '    solver status                    : apply'
	2       //' best attempt'
	    END IF
	    DEALLOCATE(BEST_X)
	  ELSE
	    JFNK_INFO=-1
	    JFNK_CATASTROPHIC=.TRUE.
	  END IF
!
	  FREEZE_I_NEXT=.FALSE.
!
	  WRITE(LUER,'(/,A)')' --------------------------------- JFNK DONE ---------------------------------'
	  POPS_PREV(1:NT,1:ND)=POPS(1:NT,1:ND)
	  IF(JFNK_CATASTROPHIC)APPLY_JFNK_UPDATE=.FALSE.
!
! If a baseline rebuild ran but produced no further correction, adopt
! the rebuilt POPS_BASE.  That point already reduced |F|.  Do not then
! take the Lambda SOLVEBA step as well: the two updates stack, and the
! second one is the clipped full-BA step that walks the model off the
! cliff.  A catastrophic attempt must NOT adopt POPS_BASE; put ED/T/ATM
! back on the outer populations before the safety step.
!
	  IF(.NOT. APPLY_JFNK_UPDATE .AND. BASELINE_REBUILT .AND.
	1     .NOT. JFNK_CATASTROPHIC)THEN
	    POPS(1:NT,1:ND)=POPS_BASE(1:NT,1:ND)
	    ED(1:ND)=POPS(NT-1,1:ND)
	    T(1:ND)=POPS(NT,1:ND)
	    CALL JFNK_SYNC_ATM_FROM_POPS(POPS)
!
! Same MAXCH trace exemption as SOLVEBA_V13: JFNK reports the major
! block (plus T/ED) to CMFGEN_SUB, while BA_TR_ABS carries the absolute
! committed trace change used to arm Stage B.  Without this a Newton
! commit could leave BA_TR_ABS=0 and skip the polish entirely.
! BA_DEM_EXM takes the committed exempted change as the best-effort
! proxy for exempted demand (the JFNK residual has no STEQ), so the
! Stage-B trigger and exit test see trace state on JFNK iterations.
!
	    MAXCH=0.0_LDP
	    BA_TR_ABS=0.0_LDP
	    BA_DEM_EXM=0.0_LDP
	    DO K=1,ND
	      T2=0.0_LDP
	      IF(BA_DEM_RELFRAC .GT. 0.0_LDP)
	1     T2=BA_DEM_RELFRAC*MAXVAL(POPS_PREV(1:NT-2,K))
	      DO J=1,NT
	        IF(POPS_PREV(J,K) .GT. 1.0E-100_LDP)THEN
	          SOL(J,K)=POPS(J,K)/POPS_PREV(J,K)-1.0_LDP
	        ELSE
	          SOL(J,K)=0.0_LDP
	        END IF
	        T1=100.0_LDP*ABS(SOL(J,K))
!
! Same exempted-row predicate as SOLVEBA_V13's demand/committed
! scans, on the same operands: POPS_PREV is the pre-step state, so it
! plays the POPS_OLD role and is also the RELFRAC baseline.  The
! PFLOOR>0 guard keeps an explicit 0 = off consistent with the BA
! side; without it the two solvers exempt different rows.
!
	        IF(J .LE. NT-2 .AND.
	1           ((BA_DEM_PFLOOR .GT. 0.0_LDP .AND.
	2             POPS_PREV(J,K) .LT. BA_DEM_PFLOOR) .OR.
	3            POPS_PREV(J,K) .LT. T2))THEN
	          BA_TR_ABS=MAX(BA_TR_ABS,ABS(POPS(J,K)-POPS_PREV(J,K))/
	1           MAX(POP_ATOM_A(K),1.0E-30_LDP))
	          BA_DEM_EXM=MAX(BA_DEM_EXM,T1)
	          CYCLE
	        END IF
	        IF(T1 .GT. MAXCH)MAXCH=T1
	      END DO
	    END DO
	    MAXCH_APP=MAXCH
	    SUCCESS=.TRUE.
	    COMMITTED_BASELINE=.TRUE.
! A lower residual is progress, not proof that the convergence target passed.
	    JFNK_UNCONVERGED_COMMIT=.TRUE.
	    WRITE(LUER,'(A)')
	1     '  JFNK: adopting rebuilt baseline (verified lower |F|)'
	    WRITE(LUER,'(A,ES12.4)')
	1     '  JFNK: committed rebuild, MAXCH=',MAXCH
	  ELSE IF(BASELINE_REBUILT)THEN
	    ED(1:ND)=POPS(NT-1,1:ND)
	    T(1:ND)=POPS(NT,1:ND)
	    CALL JFNK_SYNC_ATM_FROM_POPS(POPS)
	  END IF
!
	  IF(APPLY_JFNK_UPDATE)THEN
	    WRITE(LUER,'(/,A)')' --------------------------------- JFNK COMMIT ---------------------------------'
	    IF(JFNK_INFO .EQ. 0)THEN
	      WRITE(LUER,'(A)')
	1       '    applying    : nonlinear solution'
	    ELSE
	      WRITE(LUER,'(A)')
	1       '    applying    : best Newton correction,'
	2       //' not fully converged'
	    END IF
!
! [JFNK_LOGX] (default T): n = n0*exp(x), clipped to
! [1/MAX_JFNK, MAX_JFNK].  [JFNK_LOGX]=F keeps n = n0*(1+x).
! Free-T steps are always logarithmic and also honour MAX_dT.
!
	    X_MAX=MAX_JFNK_COR-1.0_LDP
	    X_MIN=1.0_LDP/MAX_JFNK_COR-1.0_LDP
	    IF(jfnk_log_unknown)THEN
	      X_MAX=LOG(MAX_JFNK_COR)
	      X_MIN=-X_MAX
	    END IF
	    X_T_MAX=LOG(MAX_JFNK_COR)
	    X_T_MIN=-X_T_MAX
	    IF(MDT_A .GT. 0.0_LDP .AND. MDT_A .LT. 1.0_LDP)THEN
	      X_T_MAX=MIN(X_T_MAX,LOG(1.0_LDP+MDT_A))
	      X_T_MIN=MAX(X_T_MIN,LOG(1.0_LDP-MDT_A))
	    END IF
	    NCLIP_POP=0
	    NCLIP_T=0
	    NCLIP_ED=0
	    MINSCALE_J=1.0_LDP
	    WRITE(LUER,'(A,L2)')
	1     '    T free to vary              :',.NOT. FIXED_T
	2     .AND. .NOT. RD_FIX_T
	    WRITE(LUER,'(A,L2)')
	1     '    Ne free to vary             :',.NOT. FIXED_NE
	    WRITE(LUER,'(A,2ES10.2)')
	1     '    pop/Ne unknown clip bounds  :',X_MIN,X_MAX
	    IF(.NOT. FIXED_T .AND. .NOT. RD_FIX_T)THEN
	      WRITE(LUER,'(A,2ES10.2)')
	1       '    T unknown clip bounds       :',X_T_MIN,X_T_MAX
	    END IF
!
	    DO K=1,ND
	      DO J=1,NT-2
	        T1=REAL(jfnk_x((K-1)*NT+J),KIND=LDP)
	        T1_RAW=T1
	        IF(T1 .GT. X_MAX)THEN
	          T1=X_MAX
	          NCLIP_POP=NCLIP_POP+1
	        ELSE IF(T1 .LT. X_MIN)THEN
	          T1=X_MIN
	          NCLIP_POP=NCLIP_POP+1
	        END IF
	        IF(jfnk_log_unknown)THEN
	          IF(T1 .NE. T1_RAW)
	1           MINSCALE_J=MIN(MINSCALE_J,EXP(T1-T1_RAW))
	          POPS(J,K)=POPS_BASE(J,K)*EXP(T1)
	        ELSE
	          IF(T1 .NE. T1_RAW .AND. 1.0_LDP+T1_RAW .NE. 0.0_LDP)
	1           MINSCALE_J=MIN(MINSCALE_J,
	2           ABS((1.0_LDP+T1)/(1.0_LDP+T1_RAW)))
	          POPS(J,K)=POPS_BASE(J,K)*(1.0_LDP+T1)
	          IF(POPS(J,K) .LT. 1.0E-100_LDP)POPS(J,K)=1.0E-100_LDP
	        END IF
	      END DO
	      IF(.NOT. FIXED_T)THEN
	        IF(MOD_FIX_T_D_ST .GT. 0 .AND.
	1          K .GE. MOD_FIX_T_D_ST .AND.
	2          K .LE. MOD_FIX_T_D_END .AND.
	3          .NOT. JFNK_TAUTO_ACTIVE)THEN
	        ELSE
	          T1=REAL(jfnk_x((K-1)*NT+NT),KIND=LDP)
	          T1_RAW=T1
	          IF(T1 .GT. X_T_MAX)THEN
	            T1=X_T_MAX
	            NCLIP_T=NCLIP_T+1
	          ELSE IF(T1 .LT. X_T_MIN)THEN
	            T1=X_T_MIN
	            NCLIP_T=NCLIP_T+1
	          END IF
	          IF(T1 .NE. T1_RAW)
	1           MINSCALE_J=MIN(MINSCALE_J,EXP(T1-T1_RAW))
	          POPS(NT,K)=POPS_BASE(NT,K)*EXP(T1)
	          IF(TMIN_A .GT. 0.0_LDP .AND.
	1             POPS(NT,K) .LT. TMIN_A)POPS(NT,K)=TMIN_A
	          T(K)=POPS(NT,K)
	        END IF
	      END IF
	      IF(.NOT. FIXED_NE)THEN
	        T1=REAL(jfnk_x((K-1)*NT+NT-1),KIND=LDP)
	        T1_RAW=T1
	        IF(T1 .GT. X_MAX)THEN
	          T1=X_MAX
	          NCLIP_ED=NCLIP_ED+1
	        ELSE IF(T1 .LT. X_MIN)THEN
	          T1=X_MIN
	          NCLIP_ED=NCLIP_ED+1
	        END IF
	        IF(jfnk_log_unknown)THEN
	          IF(T1 .NE. T1_RAW)
	1           MINSCALE_J=MIN(MINSCALE_J,EXP(T1-T1_RAW))
	          POPS(NT-1,K)=POPS_BASE(NT-1,K)*EXP(T1)
	        ELSE
	          IF(T1 .NE. T1_RAW .AND. 1.0_LDP+T1_RAW .NE. 0.0_LDP)
	1           MINSCALE_J=MIN(MINSCALE_J,
	2           ABS((1.0_LDP+T1)/(1.0_LDP+T1_RAW)))
	          POPS(NT-1,K)=POPS_BASE(NT-1,K)*(1.0_LDP+T1)
	          IF(POPS(NT-1,K) .LT. 1.0E-100_LDP)
	1             POPS(NT-1,K)=1.0E-100_LDP
	        END IF
	        ED(K)=POPS(NT-1,K)
	      ELSE
	        ED(K)=POPS(NT-1,K)
	      END IF
	    END DO
	    IF(NCLIP_POP+NCLIP_T+NCLIP_ED .GT. 0)THEN
	      WRITE(LUER,'(A,I8)')
	1       '    pop cells at clip bound     :',NCLIP_POP
	      WRITE(LUER,'(A,I8)')
	1       '    Ne cells at clip bound      :',NCLIP_ED
	      WRITE(LUER,'(A,I8)')
	1       '    T cells at clip bound       :',NCLIP_T
	      IF(NCLIP_POP+NCLIP_T+NCLIP_ED .GT. 10)THEN
	        WRITE(LUER,'(A)')
	1   '  WARNING: >10 components clipped - check preconditioner/RSCALE.'
	      END IF
	    END IF
!
! Verify the applied state against the frozen residual before the
! commit is allowed to stand.  The Newton accept tests all ran on
! unclipped iterates, so a heavily clipped commit would otherwise be
! applied blind: x_app is reconstructed from the committed POPS (it
! therefore includes every clip bound, floor, and skipped fixed row
! exactly) and F(x_app) is compared with the entry residual.  A
! committed state that does not beat |F|_w at entry is rolled back
! atomically to POPS_PREV and the caller-mode BA fallback runs on the
! restored state.  The coupled tier is skipped: one eval costs a
! full formal sweep, and coupled commits already carry
! JFNK_UNCONVERGED_COMMIT.
!
	    IF(.NOT. JFNK_COUPLED_LOC)THEN
	      BLOCK
	        DOUBLE PRECISION X_APP(JFNK_N),Y_APP(JFNK_N)
	        INTEGER JA
	        LOGICAL X_APP_OK
	        X_APP_OK=.TRUE.
	        DO K=1,ND
	          DO J=1,NT
	            JA=(K-1)*NT+J
	            IF(jfnk_log_unknown)THEN
	              IF(POPS(J,K) .GT. 0.0_LDP .AND.
     1                 POPS_BASE(J,K) .GT. 0.0_LDP)THEN
	                X_APP(JA)=LOG(DBLE(POPS(J,K))/DBLE(POPS_BASE(J,K)))
	              ELSE
	                X_APP(JA)=0.0D0
	                X_APP_OK=.FALSE.
	              END IF
	            ELSE
	              IF(POPS_BASE(J,K) .NE. 0.0_LDP)THEN
	                X_APP(JA)=DBLE(POPS(J,K))/DBLE(POPS_BASE(J,K))-1.0D0
	              ELSE
	                X_APP(JA)=0.0D0
	                X_APP_OK=.FALSE.
	              END IF
	            END IF
	          END DO
	        END DO
	        IF(X_APP_OK)THEN
	          CALL JFNK_RESIDUAL_V1(JFNK_N,X_APP,Y_APP)
	          CALL JFNK_MASK_FIXED_UNKNOWNS(JFNK_N,Y_APP)
	          CALL jfnk_ptc_strip_v1(JFNK_N,X_APP,Y_APP)
	          F_APP=SQRT(JFNK_DOT_V1(JFNK_N,Y_APP,Y_APP))
	          WRITE(LUER,'(A,ES12.4,A,ES12.4)')
     1          '    |F|_w at applied state    :',F_APP,
     2          '   (entry',F_W_INIT,')'
	        ELSE
	          F_APP=F_W_INIT
	        END IF
	        IF(.NOT. X_APP_OK .OR. jfnk_resid_nbad .GT. 0 .OR.
     1           .NOT. (F_APP .LT. F_W_INIT))THEN
	          WRITE(LUER,'(A)')
     1          '  JFNK: applied step did not improve the frozen'//
     2          ' residual -- rolling back, caller-mode BA follows'
	          POPS=POPS_PREV
	          DO K=1,ND
	            ED(K)=POPS(NT-1,K)
	            T(K)=POPS(NT,K)
	          END DO
	          APPLY_JFNK_UPDATE=.FALSE.
	          VERIF_REJECT=.TRUE.
	        END IF
	      END BLOCK
	    END IF
	  END IF
	  IF(APPLY_JFNK_UPDATE)THEN
!
! SOL = applied fractional change (same convention as SOLVEBA);
! MAXCH = 100*max|SOL|.
!
	    DO K=1,ND
	      DO J=1,NT
	        IF(POPS_PREV(J,K) .GT. 1.0E-100_LDP)THEN
	          SOL(J,K)=POPS(J,K)/POPS_PREV(J,K)-1.0_LDP
	        ELSE
	          SOL(J,K)=0.0_LDP
	        END IF
	      END DO
	    END DO
	    MAXCH=0.0_LDP
	    BA_TR_ABS=0.0_LDP
	    BA_DEM_EXM=0.0_LDP
	    DO K=1,ND
	      T2=0.0_LDP
	      IF(BA_DEM_RELFRAC .GT. 0.0_LDP)
	1     T2=BA_DEM_RELFRAC*MAXVAL(POPS_PREV(1:NT-2,K))
	      DO J=1,NT
	        T1=100.0_LDP*ABS(SOL(J,K))
	        IF(J .LE. NT-2 .AND.
	1           ((BA_DEM_PFLOOR .GT. 0.0_LDP .AND.
	2             POPS_PREV(J,K) .LT. BA_DEM_PFLOOR) .OR.
	3            POPS_PREV(J,K) .LT. T2))THEN
	          BA_TR_ABS=MAX(BA_TR_ABS,ABS(POPS(J,K)-POPS_PREV(J,K))/
	1           MAX(POP_ATOM_A(K),1.0E-30_LDP))
	          BA_DEM_EXM=MAX(BA_DEM_EXM,T1)
	          CYCLE
	        END IF
	        IF(T1 .GT. MAXCH)MAXCH=T1
	      END DO
	    END DO
	    MAXCH_APP=MAXCH
! Keep best-iterate and clipped commits out of CMFGEN_SUB's convergence gates.
	    JFNK_UNCONVERGED_COMMIT=(JFNK_INFO .NE. 0 .OR.
	1       NCLIP_POP+NCLIP_T+NCLIP_ED .GT. 0)
!
	    JFNK_MAXINC=-1.0_LDP
	    JFNK_MAXDEC=1.0_LDP
	    IINC=1
	    IDEC=1
	    JFNK_INC_VEC=0.0_LDP
	    JFNK_DEC_VEC=0.0_LDP
	    DO K=1,ND
	      DO J=1,NT
	        IF(SOL(J,K) .GT. JFNK_MAXINC)THEN
	          JFNK_MAXINC=SOL(J,K)
	          IINC=K
	        END IF
	        IF(SOL(J,K) .LT. JFNK_MAXDEC)THEN
	          JFNK_MAXDEC=SOL(J,K)
	          IDEC=K
	        END IF
	        T2=-SOL(J,K)
	        IF(T2 .LE. 0.0_LDP)THEN
	          IF(T2 .LT. JFNK_INC_VEC(NVJ))THEN
	            DO I=NVJ-1,1,-1
	              IF(T2 .GT. JFNK_INC_VEC(I))THEN
	                JFNK_INC_VEC(I+2:NVJ)=JFNK_INC_VEC(I+1:NVJ-1)
	                JFNK_INC_VEC(I+1)=T2
	                EXIT
	              END IF
	              IF(I .EQ. 1)THEN
	                JFNK_INC_VEC(2:NVJ)=JFNK_INC_VEC(1:NVJ-1)
	                JFNK_INC_VEC(1)=T2
	              END IF
	            END DO
	          END IF
	        ELSE
	          IF(T2 .GT. JFNK_DEC_VEC(NVJ))THEN
	            DO I=NVJ-1,1,-1
	              IF(T2 .LT. JFNK_DEC_VEC(I))THEN
	                JFNK_DEC_VEC(I+2:NVJ)=JFNK_DEC_VEC(I+1:NVJ-1)
	                JFNK_DEC_VEC(I+1)=T2
	                EXIT
	              END IF
	              IF(I .EQ. 1)THEN
	                JFNK_DEC_VEC(2:NVJ)=JFNK_DEC_VEC(1:NVJ-1)
	                JFNK_DEC_VEC(1)=T2
	              END IF
	            END DO
	          END IF
	        END IF
	      END DO
	    END DO
	    WRITE(LUER,'(A,I4,A,ES10.2,2A,I4)')
	1     ' Maximum % increase at depth ',IINC,' is',
	2     100.0_LDP*JFNK_MAXINC,JFNK_TAG,
	3     '  --- iteration ',MAIN_COUNTER
	    WRITE(LUER,'(A,I4,A,ES10.2,2A,I4)')
	1     ' Maximum % decrease at depth ',IDEC,' is',
	2     -100.0_LDP*JFNK_MAXDEC,JFNK_TAG,
	3     '  --- iteration ',MAIN_COUNTER
	    WRITE(LUER,'(3A,ES24.14)')
	1     'Maximm changes as returned by JFNK_SOLVE_CORR_V1',
	2     JFNK_TAG(1:LEN_TRIM(JFNK_TAG)),' is',MAXCH
	    WRITE(LUER,'(A,10ES10.2)')' DEC_VEC: ',JFNK_DEC_VEC(1:NVJ)
	    WRITE(LUER,'(A,10ES10.2)')' INC_VEC: ',JFNK_INC_VEC(1:NVJ)
	    WRITE(LUER,'(A,1PE12.4)')
	1     ' The minimum value of scale for the JFNK update is:',
	2     MINSCALE_J
	    SUCCESS=.TRUE.
	    JFNK_STAGNATION_COUNT=0
	    JFNK_REARM_CYCLE=0
	  ELSE IF(COMMITTED_BASELINE)THEN
	    JFNK_STAGNATION_COUNT=0
	    JFNK_REARM_CYCLE=0
	  ELSE IF(JFNK_CATASTROPHIC)THEN
	    WRITE(LUER,'(A,I3)')
	1     '  JFNK: CATASTROPHIC eval failure (NaN/Inf or factor '//
	2     'breakdown), info=',JFNK_INFO
! editted by Devin, 04-Oct-2026: a catastrophic failure used to reset
! the stagnation count, so a deterministic catastrophic mode retried
! forever (each admission costs a full coupled formal sweep).  Count
! it toward JFNK_MAX_STAGNATE like any other failed solve; only a
! converged Newton or a committed best-iterate resets the budget.
! W7 (17-Oct-2026): name the failure class -- a catastrophic eval
! failure (NaN/Inf, factor breakdown) is not GMRES stagnation.
	    JFNK_STAGNATION_COUNT=JFNK_STAGNATION_COUNT+1
	    IF(JFNK_STAGNATION_COUNT .GE. JFNK_MAX_STAGNATE)THEN
	      WRITE(LUER,'(A,I3,A,I3)')
	1       '  JFNK: CATASTROPHIC failure count ',JFNK_STAGNATION_COUNT,
	2       ' consecutive solves (max=',JFNK_MAX_STAGNATE,')'
	      JFNK_DISABLED=.TRUE.
	      JFNK_STAGNATION_COUNT=0
	      JFNK_REARM_CYCLE=JFNK_REARM_CYCLE+1
	    ELSE
	      WRITE(LUER,'(A)')
	1       '  JFNK remains enabled; Lambda safety step this iteration'
	    END IF
	    LAMBDA_ITERATION=.TRUE.
	    FIXED_T=.TRUE.
	    COMPUTE_BA=.TRUE.
	    TEMP_CHAR='DIAG'
! Failed Newton must not apply MAX_LAM.  The seq run's safety step
! followed |F| ~ 1e20.  Use the full-BA limit -- clamped by the
! caller's cap: a GS-damped or BA_DAMP_CAP-bounded COR_LIM_IN must
! survive the fallback (KNOP_jfnkrel it35: MAX_LAM_COR here let a
! +91% free-electron step through on the GS release call and the
! lagged field ran away, se_res -> 2E31).
	    T1=JFNK_FALLBACK_CAP(.TRUE.,COR_LIM_IN,.TRUE.)
	    CALL SOLVEBA_V13(SOL,POPS,POP_ATOM_A,DIAG_INDX,NT,NION,
	1       NUM_BNDS,ND,MAXCH,MAXCH_APP,TEMP_CHAR,SUCCESS,SCALE_OPT_A,
	2       LAM_SC_OPT_A,T1,MDT_A,TMIN_A,
	3       COMPUTE_BA,WR_BA_INV_A,WR_PRT_INV_A,LAMBDA_ITERATION,
	4       MAIN_COUNTER,SET_D2_A)
	  ELSE
! info=-2 is the reachability precheck skip: Newton never ran, so the
! admission cost ~one matvec.  Do not count it toward stagnation; the
! residual may become reachable at a later outer iterate.
	    IF(VERIF_REJECT)THEN
	      WRITE(LUER,'(A)')
     1       '  JFNK: applied correction rejected by the'//
     2       ' applied-state residual check; state restored to entry'
	    END IF
	    IF(JFNK_INFO .EQ. -2)THEN
	      WRITE(LUER,'(A)')
	1       '  JFNK: admission gate skipped Newton;'//
	2       ' retaining caller mode (not counted as stagnation)'
	    ELSE
	      JFNK_STAGNATION_COUNT=JFNK_STAGNATION_COUNT+1
	      IF(JFNK_STAGNATION_COUNT .GE. JFNK_MAX_STAGNATE)THEN
	        WRITE(LUER,'(A,I3,A,I3)')
	1         '  JFNK: stagnated ',JFNK_STAGNATION_COUNT,
	2         ' consecutive solves (max=',JFNK_MAX_STAGNATE,')'
	        JFNK_DISABLED=.TRUE.
	        JFNK_STAGNATION_COUNT=0
	        JFNK_REARM_CYCLE=JFNK_REARM_CYCLE+1
	      ELSE
	        WRITE(LUER,'(A,I3,A,I3)')
	1         '  JFNK: finite GMRES stagnation; retaining caller '//
	2         'mode, count ',
	3         JFNK_STAGNATION_COUNT,' of ',JFNK_MAX_STAGNATE
	      END IF
	    END IF
! Residual callbacks restore POPS, ED, T and STEQ.  The Newton
! preconditioner is local, so a clean FULL fallback can reuse the outer
! BA.  Restore the entry flags because called routines may update them.
	    LAMBDA_A=LAMBDA_SAVE
	    FIXED_T=FIXED_T_SAVE
	    COMPUTE_BA_A=COMPUTE_BA_SAVE
	    WR_BA_INV_A=WR_BA_INV_SAVE
	    WR_PRT_INV_A=WR_PRT_INV_SAVE
	    TEMP_CHAR=TEMP_CHAR_SAVE
! The LAMBDA fallback must not raise the caller's cap either: during
! the GS release window COR_LIM_IN carries GS_LAM_CAP and the
! field-poisoning step above applies here unchanged.
	    IF(LAMBDA_A)THEN
	      T1=JFNK_FALLBACK_CAP(.TRUE.,COR_LIM_IN,.FALSE.)
	    ELSE
	      T1=JFNK_FALLBACK_CAP(.FALSE.,COR_LIM_IN,.FALSE.)
	    END IF
	    CALL SOLVEBA_V13(SOL,POPS,POP_ATOM_A,DIAG_INDX,NT,NION,
	1       NUM_BNDS,ND,MAXCH,MAXCH_APP,TEMP_CHAR,SUCCESS,SCALE_OPT_A,
	2       LAM_SC_OPT_A,T1,MDT_A,TMIN_A,
	3       COMPUTE_BA_A,WR_BA_INV_A,WR_PRT_INV_A,LAMBDA_A,
	4       MAIN_COUNTER,SET_D2_A)
	  END IF
!
! FIX_T_AUTO ladder: measure the committed log-T step this activation
! (jfnk_x holds the accepted correction, or the best-iterate on a
! partial solve).  A step above JFNK_TAU_EPS reverts a rung and
! clears the commit streak; a rung releases only after JFNK_TAU_NLOCK
! committed activations with RE still under JFNK_TAU_REL*RE0.
!
	  IF(JFNK_TAUTO_ACTIVE)THEN
	    TAU_EPS_T=0.0D0
	    DO K=JFNK_TAU_DOM_ST,JFNK_TAU_DOM_END
	      TAU_EPS_T=MAX(TAU_EPS_T,
	1       ABS(jfnk_x((K-1)*JFNK_NT+JFNK_NT)))
	    END DO
	    IF(SUCCESS .OR. COMMITTED_BASELINE)
	1     JFNK_TAU_NCOMMIT=JFNK_TAU_NCOMMIT+1
	    IF(TAU_EPS_T .GT. DBLE(JFNK_TAU_EPS))THEN
	      IF(JFNK_TAU_RUNG .GT. 1)THEN
	        JFNK_TAU_RUNG=JFNK_TAU_RUNG-1
	        JFNK_TAU_NCOMMIT=0
	        WRITE(LUER,'(A,ES10.3,A,I2)')
	1         '  FIX_T_AUTO: committed log-T step ',
	2         REAL(TAU_EPS_T),' exceeds cap -- reverted to rung ',
	3         JFNK_TAU_RUNG
	      END IF
	    ELSE IF(JFNK_TAU_NCOMMIT .GE. JFNK_TAU_NLOCK .AND.
	1       JFNK_TAU_RE_LAST .LE. JFNK_TAU_REL*JFNK_TAU_RE0 .AND.
	2       JFNK_TAU_RUNG .LT. 4)THEN
	      JFNK_TAU_RUNG=JFNK_TAU_RUNG+1
	      JFNK_TAU_NCOMMIT=0
	      WRITE(LUER,'(A,ES10.3,A,I2)')
	1       '  FIX_T_AUTO: energy residual at ',
	2       REAL(JFNK_TAU_RE_LAST/MAX(JFNK_TAU_RE0,1.0E-30_LDP)),
	3       ' of baseline -- released to rung ',JFNK_TAU_RUNG
	    END IF
	    SELECT CASE(JFNK_TAU_RUNG)
	    CASE(1)
	      JFNK_TAU_ALPHA=JFNK_TAU_A0
	    CASE(2)
	      JFNK_TAU_ALPHA=0.2_LDP
	    CASE(3)
	      JFNK_TAU_ALPHA=0.5_LDP
	    CASE DEFAULT
	      JFNK_TAU_ALPHA=1.0_LDP
	    END SELECT
! The damp routine reads the solver-module copy; keep it on the
! current rung so a release actually frees the step.
	    jfnk_tauto_alpha=DBLE(JFNK_TAU_ALPHA)
	  END IF
!
	  WRITE(LUER,'(A)')'  JFNK final diagnostics:'
	  WRITE(LUER,'(A,I6)')'    Newton-level F evaluations:',
	1   jfnk_total_f_evals
	  WRITE(LUER,'(A,I6)')'    Residual evals incl JVPs: ',
	1   JFNK_N_RESIDUAL_EVALS
	  WRITE(LUER,'(A,I6)')'    Total JVP evaluations:',
	1   jfnk_total_jvp_evals
	  WRITE(LUER,'(A,I6)')'    Preconditioner calls: ',
	1   JFNK_N_PRECOND_APPLIES
	  WRITE(LUER,'(A,I6)')'    Preconditioner fails: ',
	1   jfnk_precond_fails
	  WRITE(LUER,'(A,I6)')'    NaN rejected steps:   ',
	1   jfnk_nan_rejects
	  WRITE(LUER,'(A,I6)')'    Residual NaN evals:   ',
	1   JFNK_N_RESIDUAL_NAN
!
! Phase-time summary: cmat = frozen-mode STEQ split + BA_rad build;
! pc = preconditioner factorization; f0 = F(0) residual test; prec =
! reach/direction-gate precheck and dirchk diagnostics; newton =
! NEWTON_HOOK_V1 attempts (all GMRES/JVP/line-search cost).  'other'
! is the residual -- entry gating, commit/rollback, and untimed
! bookkeeping.  Written on every JFNK admission so the phase timers
! can be diffed against the iteration wall clock.
!
	  WRITE(LUER,'(A,5(F8.1,A),F8.1,A)')
	1   '    JFNK phase time [s]: cmat=',T_PH_CMAT,
	2   ' pc=',T_PH_PC,' f0=',T_PH_F0,' prec=',T_PH_PREC,
	3   ' newton=',T_PH_NWT,' total=',OMP_GET_WTIME()-T_PH_TOT
	  WRITE(LUER,'(A)')' '
!
! APPLY_JFNK_UPDATE is the final commit verdict: verify-rejects,
! catastrophic paths, and every defer/fallback leave it false, so the
! NG window bookkeeping can tell a real Newton iterate from an
! admitted call that merely ran caller-mode BA.
!
	  JFNK_COMMITTED=APPLY_JFNK_UPDATE
!
! Reach-skip streak for the SOB relief de-rate (17-Oct-2026): only
! the JFNK_REACH_GATE skip counts; dirchk defers, guard defers, real
! Newton outcomes and stagnations are all non-reach outcomes and
! reset the streak, as does a committed step.
!
	  IF(JFNK_REACH_SKIP)THEN
	    JFNK_REACH_SKIP_STREAK=JFNK_REACH_SKIP_STREAK+1
	  ELSE
	    JFNK_REACH_SKIP_STREAK=0
	  END IF
!
	ELSE
!
! LAMBDA mode: standard direct solve (same as the default path).
!
! Same cap rule on the paused/LAMBDA direct-solve path: the caller's
! damped cap (GS release window, BA_DAMP_CAP) wins over MAX_LAM_COR.
	  IF(LAMBDA_A)T1=JFNK_FALLBACK_CAP(.TRUE.,COR_LIM_IN,.FALSE.)
	  CALL SOLVEBA_V13(SOL,POPS,POP_ATOM_A,DIAG_INDX,NT,NION,
	1     NUM_BNDS,ND,MAXCH,MAXCH_APP,TEMP_CHAR,SUCCESS,SCALE_OPT_A,
	2     LAM_SC_OPT_A,T1,MDT_A,TMIN_A,
	3     COMPUTE_BA,WR_BA_INV_A,WR_PRT_INV_A,LAMBDA_ITERATION,
	4     MAIN_COUNTER,SET_D2_A)
!
	END IF
	IF(JFNK_LAG_FACTORS_ACTIVE)THEN
	  JFNK_LAG_FACTORS_ACTIVE=.FALSE.
	  CALL JFNK_SYNC_ATM_FROM_POPS(POPS)
	END IF
	JFNK_LAG_FACTORS_READY=.FALSE.
!
	RETURN
!
	CONTAINS
!
! edited by Devin, 17-Oct-2026: single source for the fallback
! correction-cap mapping.  The invariant shared with SOLVE_FOR_POPS:
! LAMBDA mode caps at MAX_LAM_COR, full/BA mode at MAX_LIN_COR,
! GS_LAM_CAP is a tighter caller bound for KNOP-forced LAMBDA calls,
! and a catastrophic JFNK failure deliberately uses the full-BA limit
! (FORCE_FULL_LIM).  Every fallback answer is additionally capped by
! COR_LIM_IN, the cap the caller requested this iteration -- a
! GS-damped (GS_LAM_CAP) or BA_DAMP_CAP-bounded value must survive
! the fallback (KNOP_jfnkrel it35: an uncapped MAX_LAM_COR fallback
! passed a +91% free-electron step on the GS release call and the
! lagged field ran away, se_res -> 2E31).
!
	FUNCTION JFNK_FALLBACK_CAP(IS_LAMBDA,COR_LIM,FORCE_FULL_LIM)
	REAL(KIND=LDP) JFNK_FALLBACK_CAP
	LOGICAL, INTENT(IN) :: IS_LAMBDA,FORCE_FULL_LIM
	REAL(KIND=LDP), INTENT(IN) :: COR_LIM
	REAL(KIND=LDP) BASE
	BASE=MAX_LIN_COR
	IF(IS_LAMBDA)BASE=MAX_LAM_COR
	IF(FORCE_FULL_LIM)BASE=MIN(BASE,MAX_LIN_COR)
	JFNK_FALLBACK_CAP=MIN(BASE,COR_LIM)
	END FUNCTION JFNK_FALLBACK_CAP
!
	END
