!
! Routine designed to solve the linearized statistical and radiative equilibrium
! equations for the new population estimates. Various tricks (e.g., NG Acceleration)
! are done to help ensure convergence. Routine was cut from CMFGEN_SUB in order
! to facilitate latter changes.
!
	SUBROUTINE SOLVE_FOR_POPS(POPS,NT,NION,ND,NC,NP,NUM_BNDS,DIAG_INDX,
	1            MAXCH,MAIN_COUNTER,IREC,LU_SE,LUSCR,LST_ITERATION)
	USE SET_KIND_MODULE
	USE MOD_CMFGEN
        USE ANG_QW_MOD
	USE CONTROL_VARIABLE_MOD
	USE JFNK_CMFGEN_INT_MOD, ONLY: JFNK_REACH_SKIP_STREAK
	IMPLICIT NONE
!
! Altered:   16-Oct-2026 : Stage-B trace polish forces the BA dispatch so
!                          SOLVEBA_V13 can hold majors/T/ED (edited by Devin).
! Altered:   15-Oct-2026 : BA_DIRCHK_MAX drops the dirchk admission bypass
!                          after the diagnostic budget is spent (edited by Devin).
! Altered:   05-Oct-2026 : Added population-solve subheader to OUTGEN (edited by Devin).
! Altered:   26-Apr-2019 : Added "CALL SMOOTH_POPS_AS_WE_ITERATE" (done earlier on OSIRIS).
! Altered:   18-Nov-2016 : Before doing AV, check NG acceleration was not done recently.
! Altered:   01-Apr-2015 : Changed ESTAU call to ESTAU_V2 to take clumping to account.
! Altered:   21-Jan-2014 : Changed to CALL SOLVEBA_V12 - now pass MAIN_COUNTER.
! Aleterd:   31-Dec-2013 : No longer set COMPUTE_BA=T when T_MIN_BA_EXTRAP is true. We may need to
!                            change this.
! Altered:   26-Mar-2012 : Changed to CALL SOLVEBA_V11. Has POP_ATOM.
! Altered:   12-Mar-2012 : Changed to CALL SOLVEBA_V10. Allows populations at depth 2 to be
!                            equated to those at depth 1.
! Altered:   18-May-2010 : Change to allow BA to be held fixed after a LAMBDA iteration.
! Altered:   23-Feb-2007 : Call to SOLVEBA_V8 changed to SOLVEBA_V9; LAM_SCALE_OPT inserted
!                            into SOLVEBA_V9 call.
! Finalized: 15-Feb-2006
!
	INTEGER NT
	INTEGER NION
	INTEGER ND
	INTEGER NC
	INTEGER NP
	INTEGER NUM_BNDS
	INTEGER DIAG_INDX
!
	REAL(KIND=LDP) POPS(NT,ND)
	REAL(KIND=LDP) MAXCH
	REAL(KIND=LDP) MAXCH_APP
!
	INTEGER MAIN_COUNTER
	INTEGER IREC
	INTEGER LU_SE
	INTEGER LUSCR
	LOGICAL LST_ITERATION
!
! Local data
!
	REAL(KIND=LDP) SOL(NT,ND)
	REAL(KIND=LDP) R_OLD(ND)
	REAL(KIND=LDP) TA(ND)
	REAL(KIND=LDP) TB(ND)
	REAL(KIND=LDP) ESEC(ND)
	REAL(KIND=LDP) T1,T2
	REAL(KIND=LDP), SAVE :: MAXCH_SUM=0.0_LDP
	INTEGER, SAVE :: COUNT_FULL_LAM_BA=1
	INTEGER WLOC(2)
	REAL(KIND=LDP) SIGT
!
! Non-monotonic (DO_NONMONO_KNOP) iteration-control state.
! ITS_SINCE_START counts calls within this process: unlike
! MAIN_COUNTER (checkpointed), it reinitializes on each requeue,
! which is exactly when the process-local I_NEXT arrays in
! FG_J_CMF_NONMONO_V1 are reallocated and zeroed.
!
	INTEGER, SAVE :: ITS_SINCE_START=0
	INTEGER, SAVE :: ITER_TYPE_LAST_IT=-1
	CHARACTER(LEN=128) :: ITER_MSG
	LOGICAL FORCE_LAMBDA
	LOGICAL LAM_AT_ENTRY
	LOGICAL LAM_WAS		!LAMBDA_ITERATION on entry to mode controller
	LOGICAL GS_RWARM_SKIP
	LOGICAL GS_FG_IMMATURE
	LOGICAL GS_OK
	REAL(KIND=LDP) GS_CHG_FILE,GS_MAXCHG
	INTEGER GS_CAP_FILE,GS_CAPPED
	INTEGER*8 GS_NBAD,GS_NTOT,GS_NCONV
!
! W1-W9 (17-Oct-2026, KNOP maturity-chain repair):
!
! GS_NDIV / GS_NRAY: diverged-ray count and solved-ray count, gated
! on share not on the RESID=1.0 pin (W2).  GS_NBAD_BLK: best-pass
! bad census; the bulk rule reads MIN(GS_NBAD,GS_NBAD_BLK) (W3).
! GS_SW_MAXCHG / GS_SW_NCONV / GS_RF_DCHG: sweep-semantics drift
! and retained-field delta, kept separate from the refresh
! fixed-point metrics (W4).  GS_UNSAFE is the catastrophic
! backstop tier (>10% capped cells, >1% diverged rays, >10% raw
! bad cells); GS_PIN is the strict-freeze tier (>0.1% capped).
! FG_TRIAL mirrors MOM_TRIAL for the USE_LAM_ES=T path (W8).
! TRIAL_MIN_USED carries the fired trial's refund floor so the
! accepted-trial refund below applies to either arm (W8).
! GS_RAW_NUDGE_DONE: one-shot for the GS_RAW_MATURE=F nudge.
!
	INTEGER GS_NDIV,GS_NRAY
	INTEGER*8 GS_NBAD_BLK,GS_SW_NCONV
	REAL(KIND=LDP) GS_SW_MAXCHG,GS_RF_DCHG
	REAL(KIND=LDP) GS_MAX_ALL
	INTEGER*8 GS_NBAD_EFF,GS_NCONV_ALL
	INTEGER*8 GS_NBAD_STAL,GS_NBAD_PLAT
	LOGICAL GS_UNSAFE,GS_PIN,FG_TRIAL
	INTEGER TRIAL_MIN_USED
	LOGICAL, SAVE :: GS_RAW_NUDGE_DONE=.FALSE.
	LOGICAL MOM_RAN,MOM_OK
	REAL(KIND=LDP) MOM_CHG_FILE,MOM_MAXCHG
	INTEGER*8 MOM_CAP_FILE,MOM_CAPPED,MOM_NEGJ
	INTEGER*8 MOM_NBAD5,MOM_NBAD50,MOM_NTOT,MOM_NQUAR,MOM_NPQUAR
	REAL(KIND=LDP) MOM_SW_MAXCHG
	INTEGER*8 MOM_SW_NBAD5,MOM_SW_NBAD50,MOM_SW_NTOT
	INTEGER*8 MOM_SW_NEGJ,MOM_SW_NQUAR,MOM_SW_CAPPED
	INTEGER*8 MOM_SW_NHUGE
	INTEGER*8 MOM_SW_BNDFL
	LOGICAL MOM_RF_RAN,MOM_RF_CONV
	REAL(KIND=LDP) MOM_RF_DCHG
	INTEGER*8 MOM_RF_N5,MOM_RF_N50
	LOGICAL MOM_UNSAFE_B,MOM_FREEZE_B
	LOGICAL MOM_TRIAL,MOM_BROAD
	REAL(KIND=LDP), PARAMETER :: GS_MATURE_TOL=0.05_LDP
!
! GS_UNSAFE_TOL: in-run discriminator.
! If the last sweep's lagged coupling field is still moving by more
! than this, or cap/sign-flip events were recorded, a BA step would
! commit population corrections to a radiation field that is still
! moving -- force LAMBDA for that iteration instead.
!
	REAL(KIND=LDP), PARAMETER :: GS_UNSAFE_TOL=0.5_LDP
!
! GS_LAM_CAP / GS_LAM_TRANSIENT: maturity-scaled LAMBDA correction
! limit on frozen (immature-field) iterations.  MAX_LAM_COR=100 is
! calibrated for the cold-start LTE->NLTE relaxation, where factor-100
! corrections are legitimate.  Once the lagged field is nearly mature
! the SE solve is linearizing about a field that is still drifting,
! so the same allowance lets a single compression-depth update
! re-poison the store for ~10 iterations (MOM run 12007252, it65:
! +688% at depth 50 -> sweep nbad50 44k -> 598k).  When the bulk
! churn fraction falls below GS_LAM_TRANSIENT the field is in the
! near-mature regime and the correction is capped at factor
! GS_LAM_CAP; above it the run is in gross transient and keeps
! MAX_LAM_COR so cold start is unchanged.
!
	REAL(KIND=LDP), PARAMETER :: GS_LAM_CAP=3.0_LDP
	REAL(KIND=LDP), PARAMETER :: GS_LAM_TRANSIENT=0.02_LDP
	REAL(KIND=LDP) GS_CHURN
!
! Convergence diagnostics (dev7, 30-Sep-2026 session).  FIXT_WAS and
! BA_REJ_THIS record the state the just-run SOLVEBA actually used, so
! the ITER_SUM report survives the mode controller overwriting the
! next-iteration flags.  The HIST_* ring buffers feed the
! CONV_PLATEAU detector.
!
	LOGICAL FIXT_WAS
	LOGICAL BA_REJ_THIS
	LOGICAL BA_PART_WAS
	LOGICAL PLATEAU_NOW
	LOGICAL HIST_EXIST
	REAL(KIND=LDP) GS_DRIFT_FRAC
	CHARACTER*12 MODE_NOW,MODE_NEXT
	CHARACTER*64 ITER_REASON
	CHARACTER*7 BA_CH
	EXTERNAL APPEND_ITER_REASON
	INTEGER, PARAMETER :: PLAT_WIN=6
	REAL(KIND=LDP), SAVE :: HIST_APP(PLAT_WIN)=0.0_LDP
	REAL(KIND=LDP), SAVE :: HIST_RES(PLAT_WIN)=0.0_LDP
	REAL(KIND=LDP), SAVE :: HIST_GS(PLAT_WIN)=0.0_LDP
	INTEGER, SAVE :: HIST_POS=0,N_HIST=0,PLAT_CNT=0
!
! GS_FREEZE/N_FROZEN (Knop Paper-II nested iteration): while the
! lagged coupling field is immature, the SOLVEBA update is skipped
! entirely so the great iteration becomes a transfer-only sweep that
! converges the wavelength boundary-value problem for FIXED
! populations.  N_FROZEN drives escalation, not the update choice:
! the freeze itself is uncapped by default (the field must mature on
! its own), but GS_STATE is published for the wrapper watchdog and
! GS_FREEZE_LIM can stop the model outright.  A forced SE update on
! a field still moving by ~1e6
! kicked the lagged coupling to ~1e16 and blew up the next formal
! solve (KNOP_MOM_TFIT it16-19, 20-Sep-2026), so the SE update now
! waits for real maturity rather than a counter.
!
	LOGICAL GS_FREEZE,GS_RAW_IMMATURE
	INTEGER, SAVE :: N_FROZEN=0
	INTEGER, SAVE :: N_TRIAL_ACC=0
! Release-cycle accounting (17-Oct-2026): accepted trials refund
! N_FROZEN, so the freeze->release->refreeze loop could cycle
! forever under GS_FREEZE_LIM.  N_TRIALS_REJ counts rejected
! resacc-gated releases for the GS_RELEASE_LIM stop; N_RELEASES_TOT
! counts accepted ones for telemetry.  Both ride GS_STATE.
	INTEGER, SAVE :: N_TRIALS_REJ=0
	INTEGER, SAVE :: N_RELEASES_TOT=0
	LOGICAL TRIAL_FIRED
!
! GS_WAS_FROZEN records the previous call's freeze state so the first
! live update after a frozen stretch can be damped.  While pops are
! held the deferred SE correction keeps growing against the maturing
! field (attempted max stays ~1e9+); releasing it under the cold-start
! MAX_LAM_COR cap commits a factor-100 change at depth and re-poisons
! the field the freeze just converged (LAM_dbg it7: maxapp=1e7 ->
! se2 16 -> 5e7 -> refreeze, repeated identically twice).  The
! release step uses GS_LAM_CAP instead -- a field that just cleared
! the maturity gate is in the near-mature regime that cap was
! calibrated for.  N_REL_SINCE counts consecutive live-update calls
! after a frozen stretch so the GS_REL_N damp window covers the
! whole re-lag transient, not just the first step.
!
	LOGICAL, SAVE :: GS_WAS_FROZEN=.FALSE.
	INTEGER, SAVE :: N_REL_SINCE=0
!
! Luminosity is diagnostic only in shocked models: local deposition
! and mechanical-to-radiative conversion mean the inner and outer
! luminosities need not agree or settle together.  Keep the ratio and
! call-to-call drift in OUTGEN, but never use either in mode selection,
! correction damping, population holds, or Newton admission.
!
	LOGICAL LUM_BLOCKED,LUM_POPS_HELD
	LOGICAL, SAVE :: LUM_INIT=.FALSE.
	REAL(KIND=LDP), SAVE :: LUM_PREV=0.0_LDP
	REAL(KIND=LDP) LUM_DCHG_REL
! Runaway-detector state: the ratio (LUM_RATIO_OI) cannot see a
! divergence whose inner luminosity pegs at zero, so the guard
! tracks the absolute outer value.
	REAL(KIND=LDP), SAVE :: LUM_PREV_OUT=0.0_LDP
	INTEGER, SAVE :: LUM_DIV_NCNT=0
!
! Consecutive calls whose SCRTEMP record is a clean committed
! iterate (no reject/masked/soft-reject/pending-arm/Newton/held-pops
! record).  NG requires 4 clean records before firing -- see the
! bookkeeping block after BA_PART_WAS below.
!
	INTEGER, SAVE :: NG_CLEAN=0
!
! NG_DEFERRED marks a NEXT_NG=MAIN_COUNTER+1 target armed by the
! NG_CLEAN deferral itself.  Auxiliary SCRTEMP writes below (R-grid
! revision, AV averaging) bump MAIN_COUNTER between the deferral
! and the due call and can step over the +1 target, leaving NEXT_NG
! permanently below MAIN_COUNTER so the upstream .EQ. fire test
! never fires.  The flag widens the due test to .LE. for deferred
! targets only; sentinel values (1000/1500/2000, far future) and
! every other scheduler assignment keep the upstream .EQ.
! semantics, and any reschedule clears it.
!
	LOGICAL, SAVE :: NG_DEFERRED=.FALSE.
!
! NG_DEM_PRE (edited by Devin, 15-Oct-2026): the demanded MAXCH on
! the call that fired the last NG acceleration.  ITER_SUM prints it
! against the next call's demand so an acceleration that pushed the
! model backwards (m15: MAXCH 5.3 -> 13.0 across it181) is visible
! as a rebound instead of two unrelated numbers.
!
	REAL(KIND=LDP), SAVE :: NG_DEM_PRE=-1.0_LDP
!
! BA_MATRIX_STORED records that a call was entered with
! LAMBDA_ITERATION=F and COMPUTE_BA=T -- i.e. an iteration whose SE
! phase actually ran STORE_BA_DATA_V3 and wrote BAMAT.  The fix-BA
! logic below can switch COMPUTE_BA off so later phases recall the
! stored matrix; flipping it off in the same call that first left
! LAMBDA mode (MOM_TRIAL flips the flag mid-call, after the SE phase
! already ran under LAMBDA) recalls a BAMAT that was never written
! and stops the run (MOM_dbg it31: "cant read BA File").
!
	LOGICAL, SAVE :: BA_MATRIX_STORED=.FALSE.
	LOGICAL, SAVE :: PREV_COMPUTE_ENT=.FALSE.
	LOGICAL LAMBDA_ENT,COMPUTE_ENT
!
! BA_REUSE_STALL counts consecutive judged recall steps (armed
! steps committed under a recalled matrix) whose RES2 ratio
! exceeded BA_REUSE_TOL -- the stored BA matrix stopped producing
! meaningful contraction.  A flat judgment on a fresh-matrix step
! does not count: rebuilding an already-current matrix cannot fix
! a hard nonlinear step.  At BA_REUSE_NSTALL the FIX_BA override
! below forces one rebuild.  SUM_FIRED records that the MAXCH_SUM
! drift trigger fired this call: accumulated population movement
! large enough to invalidate the linearization is not overridden
! by residual evidence.
!
	INTEGER, SAVE :: BA_REUSE_STALL=0
	LOGICAL SUM_FIRED
!
! ITER_LEG_DONE prints the ITER_SUM field legend once per process.
! The *_ST strings carry the per-row [OK]/[WARN]/[BAD]/[HELD]
! status tag written into the summary table.
!
	LOGICAL, SAVE :: ITER_LEG_DONE=.FALSE.
	CHARACTER(LEN=7) POP_ST,SE_ST,BA_ST,LUM_ST,FG_ST,MOM_ST,JFNK_ST
!
! N_MATURE counts consecutive calls on which the raw maturity gate
! passes; GS_REL_STABLE consecutive mature sweeps are required before
! the freeze releases, so a single marginal sweep of a churning field
! cannot trigger a premature live update.
!
	INTEGER, SAVE :: N_MATURE=0
	EXTERNAL GET_FG_GS_MATURITY,GET_MOM_GS_MATURITY
!
! KN_GATE_DIVERGENCE prototype state: ring of the last GD_MMAX
! refresh sweeps' max residuals.  SAVE so the history persists
! across calls.  On a warm restart two paths exist (GS_RWARM_SKIP):
! replay N_LAM_RWARM sweeps before the ring can arm the growth
! test, or -- when the GS state file restores a mature field --
! skip the replay and let the ring refill from live sweeps.
!
	INTEGER, PARAMETER :: GD_MMAX=8
	REAL(KIND=LDP), SAVE :: GS_DIV_HIST(GD_MMAX)=0.0_LDP
	INTEGER, SAVE :: GS_DIV_NH=0
	LOGICAL GDIV_UNSAFE
!
	INTEGER I,J
	INTEGER LU_GS
	INTEGER GS_IT_RD,GS_NF_RD,GS_MT_RD,GS_RL_RD,GS_TA_RD
	INTEGER GS_TR_RD,GS_NR_RD
	LOGICAL GS_FRZ_RD
	CHARACTER*80 RSV_CH
!
	LOGICAL FILE_OPEN
	LOGICAL BA_SUCCESSFULLY_OUTPUT
	LOGICAL SUCCESS
	LOGICAL JFNK_THIS
	LOGICAL RDLAM_PROBE_THIS
!
! JFNK_COUPLED_THIS: residual mode selected for THIS call by the
! two-tier admission below (coupled eval only when the tighter
! COR_GATE tier is satisfied).  JFNK_GATE_IN is the effective
! frozen-tier gate (FRZ_GATE when set, else COR_GATE).
! JFNK_FRZ_HITS / JFNK_CPL_HITS count consecutive sub-gate solves for
! the JFNK_GATE_NREQ hysteresis (audit cmfgen_obsidian/audits/JFNK_SOLVER_AUDIT finding C1).
!
	LOGICAL JFNK_COUPLED_THIS,JFNK_DEM_CTR
	REAL(KIND=LDP) JFNK_GATE_IN
	REAL(KIND=LDP) JFNK_DEM_IN
	INTEGER JFNK_NREQ_EFF
	INTEGER, SAVE :: JFNK_FRZ_HITS=0
	INTEGER, SAVE :: JFNK_CPL_HITS=0
!
! JFNK_LAST_DEM: demanded correction at the previous call, used by
! the demand-contraction term in the gate hits.  Newton is only a
! useful admission while the outer map is demonstrably contracting;
! a sub-gate demand that is rebounding (field re-poisoned by the
! last committed step, a post-Lambda transition transient) is not a
! basin signal.  JFNK_DEM_RISE bounds the tolerated call-to-call
! demand growth: the ~5x/iter contraction on this model family is
! far above noise, and a factor 1.5 rebound is the signature of a
! rejected/re-poisoned step rather than endgame jitter.
!
	REAL(KIND=LDP), SAVE :: JFNK_LAST_DEM=0.0_LDP
	REAL(KIND=LDP), PARAMETER :: JFNK_DEM_RISE=1.5_LDP
	CHARACTER*20 TEMP_CHAR
!
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
	EXTERNAL JFNK_RETREAT_RESET_V1
!
	LUER=ERROR_LU()	
!
! Reset the process-local counter before incrementing so a second
! model in the same process restarts the warm-up count.
!
	FORCE_LAMBDA=.FALSE.
	GS_FREEZE=.FALSE.
	TRIAL_FIRED=.FALSE.
	JFNK_UNCONVERGED_COMMIT=.FALSE.
	JFNK_COMMITTED=.FALSE.
!
! Consume the rd_lam probe arm into a call-local flag.  CMFGEN_SUB
! already read BA_RDLAM_PROBE to leave LAMBDA_ITERATION false for this
! call; clearing the module flag here lets its end-of-iteration
! re-impose restore LAMBDA for the next call, and a demote before the
! solve cannot leak the arm into a later call.
!
	RDLAM_PROBE_THIS=BA_RDLAM_PROBE
	BA_RDLAM_PROBE=.FALSE.
	IF(NEWMOD)THEN
	  ITS_SINCE_START=0
	  ITER_TYPE_LAST_IT=-1
	  JFNK_FRZ_HITS=0
	  JFNK_CPL_HITS=0
	  JFNK_REACH_SKIP_STREAK=0
	  JFNK_LAST_DEM=0.0_LDP
	  MAXCH_SUM=0.0_LDP
	  COUNT_FULL_LAM_BA=1
	  N_FROZEN=0
	  N_TRIAL_ACC=0
	  N_TRIALS_REJ=0
	  N_RELEASES_TOT=0
	  TRIAL_FIRED=.FALSE.
	  GS_WAS_FROZEN=.FALSE.
	  N_REL_SINCE=0
	  N_MATURE=0
	  LUM_INIT=.FALSE.
	  LUM_PREV=0.0_LDP
	  LUM_PREV_OUT=0.0_LDP
	  LUM_DIV_NCNT=0
	  N_HIST=0
	  HIST_POS=0
	  PLAT_CNT=0
	  BA_MATRIX_STORED=.FALSE.
!
! NG_CLEAN and NG_DEFERRED are likewise process-local: a second
! model must not inherit a >=4 clean-record count (NG would fire
! before its own window fills) or a deferred due-target armed by a
! previous model's counter sequence.
!
	  NG_CLEAN=0
	  NG_DEFERRED=.FALSE.
	  BA_REUSE_STALL=0
	  PREV_COMPUTE_ENT=.FALSE.
	  BAACC_RDLAM_PROBE_ARM=.FALSE.
	  BAACC_RDLAM_PROBE=.FALSE.
!
! Retreat/ladder state is likewise process-local: a second model
! must not inherit a raised JFNK budget rung, a mid-release tau
! rung, or the KNOP divergence ring (all zeroed by a requeue).
!
	  CALL JFNK_RETREAT_RESET_V1()
	  GS_DIV_HIST=0.0_LDP
	  GS_DIV_NH=0
	END IF
	ITS_SINCE_START=ITS_SINCE_START+1
!
! Entry flags: at call time these still hold the values the previous
! call's controller assigned, i.e. the mode the just-finished
! iteration ran under -- including whether its SE phase stored
! BAMAT.  Capture before the KNOP/mode blocks rewrite them for the
! next iteration.
!
	LAMBDA_ENT=LAMBDA_ITERATION
	COMPUTE_ENT=COMPUTE_BA
!
! edited by Devin, 30-Sep-2026: N_FROZEN is process-local, so a
! requeue inside a frozen stretch silently reset both GS_FREEZE_LIM
! and the wrapper's n_frozen stall watchdog.  On the first call of a
! restarted segment, adopt the streak recorded in GS_STATE by the
! previous segment -- but only when that record describes the
! immediately preceding iteration (a stale file from an abandoned
! earlier run would carry an unrelated count).
!
	IF(.NOT. NEWMOD .AND. ITS_SINCE_START .EQ. 1)THEN
	  OPEN(NEWUNIT=LU_GS,FILE='GS_STATE',STATUS='OLD',
	1    ACTION='READ',IOSTAT=IOS)
	  IF(IOS .EQ. 0)THEN
	    GS_NF_RD=-1
	    GS_IT_RD=-1
	    GS_MT_RD=0
	    GS_RL_RD=0
	    GS_TA_RD=0
	    GS_TR_RD=0
	    GS_NR_RD=0
	    GS_FRZ_RD=.FALSE.
	    DO
	      RSV_CH=' '
	      READ(LU_GS,'(A)',IOSTAT=IOS)RSV_CH
	      IF(IOS .NE. 0)EXIT
	      IF(INDEX(RSV_CH,'n_frozen=') .EQ. 1)
	1      READ(RSV_CH(10:),*,IOSTAT=IOS)GS_NF_RD
	      IF(INDEX(RSV_CH,'iter=') .EQ. 1)
	1      READ(RSV_CH(6:),*,IOSTAT=IOS)GS_IT_RD
	      IF(INDEX(RSV_CH,'frozen= T') .GT. 0)GS_FRZ_RD=.TRUE.
	      IF(INDEX(RSV_CH,'n_mature=') .EQ. 1)
	1      READ(RSV_CH(10:),*,IOSTAT=IOS)GS_MT_RD
	      IF(INDEX(RSV_CH,'n_rel=') .EQ. 1)
	1      READ(RSV_CH(7:),*,IOSTAT=IOS)GS_RL_RD
	      IF(INDEX(RSV_CH,'n_trial_acc=') .EQ. 1)
	1      READ(RSV_CH(13:),*,IOSTAT=IOS)GS_TA_RD
	      IF(INDEX(RSV_CH,'n_trials_rej=') .EQ. 1)
	1      READ(RSV_CH(14:),*,IOSTAT=IOS)GS_TR_RD
	      IF(INDEX(RSV_CH,'n_releases=') .EQ. 1)
	1      READ(RSV_CH(12:),*,IOSTAT=IOS)GS_NR_RD
	    END DO
	    CLOSE(LU_GS)
	    IF(GS_NF_RD .GE. 0 .AND.
	1      (GS_IT_RD .EQ. MAIN_COUNTER-1 .OR.
	2       GS_IT_RD .EQ. MAIN_COUNTER))THEN
	      N_FROZEN=MAX(GS_NF_RD,0)
	      N_MATURE=MAX(GS_MT_RD,0)
	      N_REL_SINCE=MAX(GS_RL_RD,0)
	      N_TRIAL_ACC=MAX(GS_TA_RD,0)
	      N_TRIALS_REJ=MAX(GS_TR_RD,0)
	      N_RELEASES_TOT=MAX(GS_NR_RD,0)
	      GS_WAS_FROZEN=GS_FRZ_RD
	      IF(CONV_DIAG_LVL .GE. 1)WRITE(LUER,'(A,I5,A,I5,A,L2)')
	1     ' GS_STATE restart: adopted n_frozen=',N_FROZEN,
	2     ' from iter',GS_IT_RD,'  frozen=',GS_FRZ_RD
	    END IF
	  END IF
	END IF
!
! edited by Devin, 30-Sep-2026: init the GS/MOM maturity slots so the
! end-of-call diagnostics (ITER_SUM/ITER_HISTORY/GS_STATE) print
! defined values on the SOB/resacc paths; the GET_*_GS_MATURITY
! calls below only fill them under DO_NONMONO_KNOP.
!
	GS_OK=.FALSE.
	GS_CHG_FILE=-1.0_LDP
	GS_CAP_FILE=0
	GS_MAXCHG=0.0_LDP
	GS_CAPPED=0
	GS_NBAD=0
	GS_NTOT=0
	GS_NCONV=0
	GS_NDIV=0
	GS_NRAY=0
	GS_NBAD_BLK=0_8
	GS_SW_NCONV=0_8
	GS_SW_MAXCHG=0.0_LDP
	GS_RF_DCHG=0.0_LDP
	GS_MAX_ALL=0.0_LDP
	GS_NBAD_EFF=0_8
	GS_NCONV_ALL=0_8
	GS_UNSAFE=.FALSE.
	GS_PIN=.FALSE.
	FG_TRIAL=.FALSE.
	TRIAL_MIN_USED=0
	GS_CHURN=0.0_LDP
	MOM_RAN=.FALSE.
	MOM_OK=.FALSE.
	MOM_CHG_FILE=-1.0_LDP
	MOM_CAP_FILE=0
	MOM_MAXCHG=0.0_LDP
	MOM_CAPPED=0
	MOM_NEGJ=0
	MOM_NBAD5=0
	MOM_NBAD50=0
	MOM_NTOT=0
	MOM_NQUAR=0
	MOM_NPQUAR=0
	MOM_SW_MAXCHG=0.0_LDP
	MOM_SW_NBAD5=0
	MOM_SW_NBAD50=0
	MOM_SW_NTOT=0
	MOM_SW_NEGJ=0
	MOM_SW_NQUAR=0
	MOM_SW_CAPPED=0
	MOM_SW_NHUGE=0
	MOM_SW_BNDFL=0
	MOM_RF_RAN=.FALSE.
	MOM_RF_CONV=.FALSE.
	MOM_RF_DCHG=0.0_LDP
	MOM_RF_N5=0
	MOM_RF_N50=0
	MOM_TRIAL=.FALSE.
!
! Keep the former BA_LUM_GATE measurement for diagnostics and restart
! comparisons, but leave LUM_BLOCKED false.  Shocks invalidate inner-
! versus-outer luminosity agreement as a convergence condition.
!
	LUM_BLOCKED=.FALSE.
	LUM_DCHG_REL=0.0_LDP
	IF(LUM_INIT)LUM_DCHG_REL=ABS(LUM_RATIO_OI-LUM_PREV)/
	1    MAX(ABS(LUM_RATIO_OI),ABS(LUM_PREV),0.1_LDP)
	LUM_INIT=.TRUE.
	LUM_PREV=LUM_RATIO_OI
!
! edited by Devin, 17-Oct-2026: luminosity runaway guard, KNOP only
! ([LUM_DIV_FAC]=1D2, [LUM_DIV_N]=3, [LUM_DIV_CAP]=1D15; 0 disables).
! Inner/outer agreement stays diagnostic -- shocked models
! legitimately disagree -- but sustained superlinear growth of the
! absolute outer luminosity is unambiguous divergence
! (KNOP_jfnkrel: 5.7E9 -> 6.2E28 L_sun in ~5 calls, then ~60 dead
! iterations until GS_FREEZE_LIM).  LUM_BLOCKED revives the dormant
! hold/trial-veto plumbing below; past LUM_DIV_CAP the state is not
! recoverable, so stop rather than burn the frozen-cycle budget.
!
	IF(DO_NONMONO_KNOP .AND. LUM_DIV_N .GT. 0)THEN
!
! A NaN/Inf luminosity fails every .GT. test below: the streak
! counter would reset and LUM_PREV_OUT would latch NaN, leaving
! the runaway invisible while it keeps feeding POPS.  Fail closed.
!
	  IF(.NOT. (ABS(LUM_OUT_ABS) .LE. HUGE(1.0_LDP)))THEN
	    WRITE(LUER,'(A)')
	1   '  SOLVE_FOR_POPS: luminosity non-finite -- stopping'
	    STOP 'SOLVE_FOR_POPS: luminosity non-finite'
	  END IF
	END IF
	IF(DO_NONMONO_KNOP .AND. LUM_DIV_N .GT. 0 .AND.
	1   LUM_OUT_ABS .GT. 0.0_LDP)THEN
	  IF(LUM_PREV_OUT .GT. 0.0_LDP .AND.
	1    LUM_OUT_ABS .GT. LUM_DIV_FAC*LUM_PREV_OUT)THEN
	    LUM_DIV_NCNT=LUM_DIV_NCNT+1
	  ELSE
	    LUM_DIV_NCNT=0
	  END IF
	  LUM_PREV_OUT=LUM_OUT_ABS
! The cap stops on the absolute level, streak or not -- the observed
! runaway plateaued at ~6E28 with growth zero, which a streak-only
! test would leave running until GS_FREEZE_LIM.
	  IF(LUM_OUT_ABS .GT. LUM_DIV_CAP)THEN
	    WRITE(LUER,'(A,ES12.4)')
	1   '  SOLVE_FOR_POPS: luminosity diverging, L_out=',
	2   LUM_OUT_ABS
	    WRITE(LUER,'(A)')'  Non-convergent state; stopping.'
	    STOP 'SOLVE_FOR_POPS: luminosity diverging'
	  END IF
	  IF(LUM_DIV_NCNT .GE. LUM_DIV_N)THEN
	    LUM_BLOCKED=.TRUE.
	    WRITE(LUER,'(A,ES12.4,A,I2,A)')
	1   '  SOLVE_FOR_POPS: luminosity runaway, L_out=',LUM_OUT_ABS,
	2   ' streak=',LUM_DIV_NCNT,' -- LUM_BLOCKED set'
	  END IF
	END IF
!
! We can adjust the BA/STEQ equations so that certain species are
! held fixed. For all species, this now done using FIXPOP_IN_BA_V2
! called by GENERATE_FULL_MATRIX.
!
	MOD_FIXED_NE=FIXED_NE
	MOD_FIX_IMPURITY=FIX_IMPURITY
	MOD_FIXED_T=.FALSE.
!
! Determine those depths where the temperature is to be held fixed.
! If either FIXED_T of RD_FIX_T is true, the temperature is fixed independent
! of VARFIXT. We zero all elements of the R.E. Eq. except the local variation
! with respect to T. For a LAMBDA iteration this element is zero, and hence
! must be set.
!
	MOD_TAU_SCL_T=TAU_SCL_T
	MOD_FIX_T_D_ST=0
	MOD_FIX_T_D_END=0
	IF(FIXED_T .OR. RD_FIX_T)THEN
	  MOD_FIXED_T=.TRUE.
	  MOD_FIX_T_D_ST=1
	  MOD_FIX_T_D_END=ND
	  CALL ESTAU_V2(TA,R,ED,CLUMP_FAC,TB,ND)
	ELSE IF(VARFIXT .AND. CON_SCL_T .NE. 0.0_LDP)THEN
	  CALL ESTAU_V2(TA,R,ED,CLUMP_FAC,TB,ND)
	  MOD_FIXED_T=.TRUE.
	  MOD_FIX_T_D_ST=1
	  DO I=1,ND
	    IF(TA(I) .LT. TAU_SCL_T)THEN
	      MOD_FIX_T_D_END=I
	    END IF
	  END DO
	  IF(MOD_FIX_T_D_END .LT. 5)MOD_FIX_T_D_END=5  	!avoid inconsistency at outer boundary.
	  WRITE(LUER,*)' '
	  WRITE(LUER,'(A,ES14.4)')' Electron scattering optical depth at outer boundary'//
	1                        ' in SOLVE_FOR_POPS is',TA(1)
	  WRITE(LUER,'(A,I3,A,I3)')' T will be held fixed over the following depths',
	1                           MOD_FIX_T_D_ST,' to ',MOD_FIX_T_D_END
	  WRITE(LUER,*)' '
	END IF
! 
!
	WRITE(LUER,'(/,A)')
	1  ' ---------------------------- POPULATION SOLVE -------------------------------'
	IF(LAMBDA_ITERATION)THEN
	  WRITE(LUER,*)'LAMBDA iteration used'
	END IF
!
	IF(COMPUTE_BA)THEN
	  WRITE(LUER,*)'BA matrix computed'
	ELSE
	  WRITE(LUER,*)'BA matrix **NOT** computed'
	END IF
	FLUSH(LUER)
!
	NLBEGIN=0	!Initialize for next iteration
	CONTINUE
	CALL WR_ASCI_STEQ(NION,ND,'STEQ ARRAY',LU_SE)
!
! If we are currently doing a LAMBDA iteration we allow for bigger
! changes - up to a factor of 100. Otherwise changes are limited to
! a factor of 10.0D0. (Installed 28-Dec-1989). Specified by parameter
! T1 in SOLVEBA call. TEMP_CHAR is used to utilize the BLK_DIAGONAL
! nature of the BA matrix when performing a LAMBDA iteration.
!
! edited by Devin, 17-Oct-2026: the correction-cap mapping is a single
! invariant -- LAMBDA -> MAX_LAM_COR, full/BA -> MAX_LIN_COR.
! GS_LAM_CAP is an additional upper bound applied to COR_LIM_IN for
! KNOP-forced LAMBDA calls while the lagged field is immature; the
! catastrophic JFNK fallback intentionally caps at MAX_LIN_COR.
! JFNK_FALLBACK_CAP (jfnk_solve_corr_v1.f) encodes the same mapping on
! the JFNK side.  SOLVEBA_V13 only consumes the cap (CHANGE_LIM), so
! correctness lives entirely at these call sites.
!
	IF(LAMBDA_ITERATION)THEN
	  T1=MAX_LAM_COR
	  TEMP_CHAR='DIAG'
	ELSE
	  T1=MAX_LIN_COR
	  TEMP_CHAR=METH_SOL
	END IF
!
! ITER_MODE transition bookkeeping: reasons are appended by the
! forcing/decision sites below and printed once next to MODE_NEXT.
!
	ITER_REASON=' '
!
! Knop Gauss-Seidel warm-up (DO_NONMONO_KNOP only).  The lagged
! I_NEXT coupling arrays in FG_J_CMF_NONMONO_V1 are process-local and
! are reallocated+zeroed whenever the job is requeued, so the first
! N_LAM_RWARM iterations of each process must not take a BA step
! against a stale radiation field.  This gate also covers the cold
! start of a new model (ITS_SINCE_START counts from 1).
!
! The warm-up is skipped when a validated, MATURE GS generation was
! restored from the FGINEXT_GS checkpoint: I_NEXT_GS_OK plus a last-
! sweep relative change below GS_MATURE_TOL and no capped writes.
!
	LAM_AT_ENTRY=LAMBDA_ITERATION
	GS_FG_IMMATURE=.FALSE.
	IF(DO_NONMONO_KNOP)THEN
	  CALL GET_FG_GS_MATURITY(GS_OK,GS_CHG_FILE,GS_CAP_FILE,
	1      GS_MAXCHG,GS_CAPPED,GS_NBAD,GS_NTOT,GS_NCONV,
	2      GS_NDIV,GS_NRAY,GS_NBAD_BLK,
	3      GS_SW_MAXCHG,GS_SW_NCONV,GS_RF_DCHG,
	4      GS_NBAD_STAL,GS_NBAD_PLAT)
	  CALL GET_MOM_GS_MATURITY(MOM_RAN,MOM_OK,MOM_CHG_FILE,
	1      MOM_CAP_FILE,MOM_MAXCHG,MOM_CAPPED,MOM_NEGJ,
	2      MOM_NBAD5,MOM_NBAD50,MOM_NTOT,MOM_NQUAR,
	3      MOM_SW_MAXCHG,MOM_SW_NBAD5,MOM_SW_NBAD50,
	4      MOM_SW_NTOT,MOM_SW_NEGJ,MOM_SW_NQUAR,MOM_SW_CAPPED,
	5      MOM_SW_BNDFL,
	6      MOM_RF_RAN,MOM_RF_CONV,MOM_RF_DCHG,MOM_RF_N5,MOM_RF_N50,
	7      MOM_NPQUAR,MOM_SW_NHUGE)
!
! The mom-side file counts as satisfied when MOM_J_CMF_NONMONO_V1
! never ran in this process (USE_LAM_ES=T): the flag then reflects
! "not needed", not "restored".
!
	  GS_RWARM_SKIP=GS_OK .AND.
	1             GS_CHG_FILE .GE. 0.0_LDP .AND.
	1             GS_CHG_FILE .LT. GS_MATURE_TOL .AND.
	1             GS_CAP_FILE .EQ. 0 .AND.
	1             (.NOT. MOM_RAN .OR. (MOM_OK .AND.
	1               MOM_CHG_FILE .GE. 0.0_LDP .AND.
	1               MOM_CHG_FILE .LT. GS_MATURE_TOL .AND.
	1               MOM_CAP_FILE .EQ. 0))
! The warm-up applies only while DO_LAM_IT is set.  With DO_LAM_IT
! false the stock cycle owns the first iteration: a full BA step,
! large correction included.  Forcing Lambda here would skip it.
!
	  IF(RD_LAMBDA .AND. ITS_SINCE_START .LE. N_LAM_RWARM .AND.
	1     .NOT. GS_RWARM_SKIP)THEN
	    IF(CONV_DIAG_LVL .GE. 1)WRITE(LUER,'(A,I3,A,I3,A,I5)')
	1     '  KNOP GS warm-up: iter',ITS_SINCE_START,' of',
	2     N_LAM_RWARM,'; forcing LAMBDA (I_NEXT reinit), iter',
	3     MAIN_COUNTER
	    LAMBDA_ITERATION=.TRUE.
	    FIXED_T=.TRUE.
	    COMPUTE_BA=.TRUE.
	    TEMP_CHAR='DIAG'
	    T1=MAX_LAM_COR
	    FORCE_LAMBDA=.TRUE.
	    CALL APPEND_ITER_REASON(ITER_REASON,'gs_warmup')
	  END IF
!
! In-run GS immaturity gate: when a BA
! iteration is pending but the lagged I_NEXT/J_NEXT coupling field is
! still moving by more than GS_UNSAFE_TOL, or cap/negative-J events
! were recorded on the last sweep, the radiation field feeding the BA
! matrix is not the converged formal solution -- force LAMBDA for this
! iteration.  FORCE_LAMBDA also stops the mode controller below from
! re-requesting BA on the next call while the field stays immature.
!
! Bulk-aware maturity: the refresh's max-norm is pinned by a handful
! of marginally-contractive cells (|Q|->1 at strong compressions)
! even when the bulk field is converged.  When the refresh ran this
! iteration (GS_NTOT>0), the field counts as mature when fewer than
! ~0.1% of significant cells remain above KN_RAY_TOL and the max
! residual is not in runaway territory.  Without refresh data the
! original max-norm criterion applies unchanged.
!
! Store drift (GS_NCONV): the refresh can converge to its tolerance
! while still moving the lagged I_NEXT store by >50% at millions of
! cells -- the fixed point exists but sits far from the field the
! in-sweep solve just used, so the next sweep reads a substantially
! different coupling field.  That is immaturity by the same >0.1%
! bulk rule, and it is the condition the residual-only gate missed
! when a BA/stale-BA SE update committed corrections against the
! churning field (LAM run it126-131 blow-up, Sep-2026).
!
!
! KN_GATE_DIVERGENCE (default F): prototype unsafe-tier variant
! that tests residual growth instead of the absolute GS_UNSAFE_TOL
! bound.  A bounded, stationary max residual -- e.g. the period-2
! limit cycle KNDUMP showed holding ~0.53 across it125-136 on
! m15_2p82d -- reads as an asymptote, not runaway.  Unsafe
! requires a non-finite value or growth by more than
! KN_GATE_DIV_FAC against the residual from KN_GATE_DIV_M refresh
! sweeps back.  The ring fills only when the refresh ran, and the
! variant has no absolute backstop: a large but stable residual
! passes.  Dormant by design; enabling is a gate policy change.
!
	  GDIV_UNSAFE=.FALSE.
	  IF(KN_GATE_DIVERGENCE .AND. GS_NTOT .GT. 0)THEN
	    IF(GS_MAXCHG .NE. GS_MAXCHG .OR.
	1      GS_MAXCHG .GT. HUGE(GS_MAXCHG))GDIV_UNSAFE=.TRUE.
	    IF(GS_DIV_NH .GE. KN_GATE_DIV_M)THEN
	      IF(GS_MAXCHG .GT. KN_GATE_DIV_FAC*
	1        GS_DIV_HIST(GD_MMAX-KN_GATE_DIV_M+1))
	2        GDIV_UNSAFE=.TRUE.
	    END IF
	    DO I=1,GD_MMAX-1
	      GS_DIV_HIST(I)=GS_DIV_HIST(I+1)
	    END DO
	    GS_DIV_HIST(GD_MMAX)=GS_MAXCHG
	    GS_DIV_NH=GS_DIV_NH+1
	  END IF
!
!
! W1-W9 (17-Oct-2026): the gate below reads three metric families
! now kept separate by the refresh rewrite:
!   max-norm   GS_MAX_ALL = MAX(refresh residual, sweep drift,
!              retained-field delta) -- the backstop.
!   bad cells  GS_NBAD_EFF = MIN(last-pass census, best-pass census)
!              -- stall/cap exits carried their stale last-pass
!              counts forever (W3).
!   churn      GS_NCONV_ALL = MAX(refresh churn, sweep churn) --
!              the refresh used to clobber the sweep metric (W4).
!
! The share tiers W1 adds: UNSAFE is the catastrophic backstop
! (>10% capped events, >1% diverged rays, >10% raw bad cells);
! PIN adds the strict >0.1% capped tier to the marginal freeze.
! One clamp event in one cell can no longer veto the run.
!
	  GS_MAX_ALL=GS_MAXCHG
	  IF(GS_NTOT .GT. 0)THEN
	    GS_MAX_ALL=MAX(GS_MAX_ALL,GS_SW_MAXCHG,GS_RF_DCHG)
	    IF(GS_NBAD_BLK .GE. 0_8)THEN
	      GS_NBAD_EFF=MIN(GS_NBAD,GS_NBAD_BLK)
	    ELSE
	      GS_NBAD_EFF=GS_NBAD
	    END IF
	    GS_NCONV_ALL=MAX(GS_NCONV,GS_SW_NCONV)
	  ELSE
	    GS_NBAD_EFF=GS_NBAD
	    GS_NCONV_ALL=GS_NCONV
	  END IF
	  GS_UNSAFE=(GS_CAPPED .GT. MAX(50_8,GS_NTOT/10)) .OR.
	1     (GS_NDIV .GT. MAX(1,GS_NRAY/100)) .OR.
	2     (GS_NBAD .GT. MAX(50_8,GS_NTOT/10))
	  GS_PIN=GS_UNSAFE .OR.
	1     (GS_CAPPED .GT. MAX(50_8,GS_NTOT/1000))
!
	  GS_FG_IMMATURE=GS_MAX_ALL .GT. GS_MATURE_TOL
	  IF(GS_NTOT .GT. 0 .AND. GS_FG_IMMATURE)THEN
	    IF(KN_GATE_DIVERGENCE)THEN
	      GS_FG_IMMATURE=(GS_NBAD_EFF .GT.
	1        MAX(50_8,GS_NTOT/1000)) .OR. GDIV_UNSAFE
	    ELSE
	      GS_FG_IMMATURE=(GS_NBAD_EFF .GT.
	1        MAX(50_8,GS_NTOT/1000)) .OR.
	2        (GS_MAX_ALL .GT. GS_UNSAFE_TOL)
	    END IF
	  END IF
	  IF(GS_NTOT .GT. 0)THEN
	    GS_FG_IMMATURE=GS_FG_IMMATURE .OR.
	1      (GS_NCONV_ALL .GT. MAX(50_8,GS_NTOT/1000))
	  END IF
!
! Plateau exits (IXIT=6) are a different failure mode from
! stall/cap: the map stopped improving entirely, so their bad
! cells are a fixed-point floor, not a slow transient.  A plateau
! share above 1% marks the field unsafe; MIN(NBAD,NBAD_BLK)
! still prices the stall/cap class the same way as before.
!
	  IF(GS_NTOT .GT. 0)THEN
	    GS_UNSAFE=GS_UNSAFE .OR.
	1      (GS_NBAD_PLAT .GT. MAX(50_8,GS_NTOT/100))
	  END IF
!
! MOM bulk maturity (MOM_GATE_BULK, VADAT): the lagged J_NEXT field is
! judged by the fraction of cells still drifting between sweeps, not
! the refresh residual -- the refresh measures lagged-map
! self-consistency, which asymptotes on marginally contractive
! compressing cells even while the stored field settles.  The sweep
! drift (fresh solve vs stored field) is what the populations
! actually experience.  negJ and quarantined cells are counted as
! fractions: under MOM_RF_KEEPNEG repaired values never reach the
! stored field, so a bounded rate is a closure diagnostic, not field
! corruption.  Runaway backstops: capped writes, >0.1% of cells
! drifting >50%, any cell drifting >10x stored, or negJ or
! quarantined fractions >10%.  edited by Devin, 28-Sep-2026: the
! former MOM_MAXCHG>1E2 term tested a clobbered slot -- after an
! unconverged refresh MOM_GS_MAXCHG holds the retained-pass
! residual (asymptote ~1E3 on marginal cells), not the sweep drift
! the test intended, so it froze the gate unconditionally.
! Exploded refresh residuals are already caught inside the refresh
! by the 1E3*MIN_RES abort and the >1E40 retention scan.  The freeze
! additionally waits until fewer than 10% of cells still drift by
! >5% between sweeps.
!
	  IF(MOM_GATE_BULK .AND. MOM_SW_NTOT .GT. 0)THEN
	    IF(MOM_RF_CONV)THEN
!
! The refresh converged this iteration: the retained J_NEXT_FRQ is
! the fixed point of the lagged coupling map, so maturity is judged
! by how far that fixed point sits from the field the populations
! actually saw (the RF delta), not by the pre-refresh sweep drift.
! NBAD/NEGJ/NQUAR/CAPPED already hold the final converged pass values.
!
	      MOM_UNSAFE_B=(MOM_CAPPED .GT. MAX(50_8,MOM_SW_NTOT/10))
	1        .OR.
	2        (MOM_RF_N50 .GT. MAX(50_8,MOM_SW_NTOT/1000)) .OR.
	3        (MOM_RF_DCHG .GT. 1.0E+1_LDP) .OR.
	4        (MOM_NEGJ .GT. MOM_SW_NTOT/10) .OR.
	5        (MOM_NQUAR .GT. MOM_SW_NTOT/10)
	      MOM_BROAD=MOM_RF_N5 .GT. MOM_SW_NTOT/10
	      MOM_FREEZE_B=MOM_UNSAFE_B .OR. MOM_BROAD .OR.
	1        (MOM_SW_BNDFL .GT. 0)
	    ELSE
!
! edited by Devin, 30-Sep-2026: split the non-converged gate into a
! hard UNSAFE tier and a marginal tier.  UNSAFE is the catastrophic
! backstop only -- >10% of cells with capped writes, >1% of cells
! drifting >50%, a cell drifting >10x in one sweep, or >10%
! negJ/quarantined.  The strict maturity terms (>0.1% of cells
! drifting >50%, any cell >10%, >10% drifting >5%, >0.1% capped
! writes) stay in FREEZE_B: they still block the transition, but
! through the marginal-release check below rather than as an
! unconditional veto.  The old single-tier test kept the FreeNe M15
! run frozen >150 iterations on ~0.12% drifting cells that asymptote
! but never converge under the one-bin upwind.
!
! edited by Devin, 30-Sep-2026 (later): the capped-write term used
! ">0" in UNSAFE, which is unreachable on this model -- the
! quarantine-repair path writes ~2000-4000 bounded cells per sweep
! (0.02-0.35% of NTOT) by design, so UNSAFE never cleared and
! MOM_TRIAL could never fire: a guaranteed stall under the marginal
! freeze.  Capped writes now follow the same scheme as negJ/quar:
! >10% of cells is the catastrophic veto, >0.1% keeps the marginal
! tier frozen.
!
! Same treatment for MOM_SW_MAXCHG: it is a single-cell max, not a
! bulk statistic.  The marginal MOM field shows a wandering argmax
! (a different optically-thin tail cell each sweep, spikes to ~5e3)
! while 99.9% of cells are mature -- under the old >1e3 veto the
! trial could never fire.  edited by Devin, 30-Sep-2026 (later):
! the >1e6 single-cell backstop had the same defect one tier up --
! a handful of bounded-but-bad resonance cells (1+sigma*K ~ 0 at the
! inner boundary, m15_2p82d depths ~209-212) sit at 1e6-1e7 while the
! field around them contracts, so mx<1e6 stayed F indefinitely.
! Collapse means many exploded cells, so UNSAFE now counts them:
! MOM_SW_NHUGE (drift >1e6) must exceed ~0.01% of cells.  A single
! nonfinite cell still vetoes via the 1e30 sentinel.
!
	      MOM_UNSAFE_B=(MOM_SW_CAPPED .GT. MAX(50_8,MOM_SW_NTOT/10))
	1        .OR.
	2        (MOM_SW_NBAD50 .GT. MAX(50_8,MOM_SW_NTOT/100)) .OR.
	3        (MOM_SW_NHUGE .GT. MAX(50_8,MOM_SW_NTOT/10000)) .OR.
	4        (MOM_SW_MAXCHG .GT. 1.0E+25_LDP) .OR.
	5        (MOM_SW_NEGJ .GT. MOM_SW_NTOT/10) .OR.
	6        (MOM_SW_NQUAR .GT. MOM_SW_NTOT/10)
	      MOM_BROAD=MOM_SW_NBAD5 .GT. MOM_SW_NTOT/10
	      MOM_FREEZE_B=MOM_UNSAFE_B .OR. MOM_BROAD .OR.
	1        (MOM_SW_NBAD50 .GT. MAX(50_8,MOM_SW_NTOT/1000)) .OR.
	2        (MOM_SW_MAXCHG .GT. 1.0E+1_LDP) .OR.
	3        (MOM_SW_CAPPED .GT. MAX(50_8,MOM_SW_NTOT/1000)) .OR.
	4        (MOM_SW_BNDFL .GT. 0)
	    END IF
	  ELSE
	    MOM_UNSAFE_B=(MOM_CAPPED .GT. 0) .OR. (MOM_NEGJ .GT. 0) .OR.
	1      (MOM_MAXCHG .GT. GS_UNSAFE_TOL)
	    MOM_BROAD=.FALSE.
	    MOM_FREEZE_B=(MOM_CAPPED .GT. 0) .OR. (MOM_NEGJ .GT. 0) .OR.
	1      (MOM_MAXCHG .GT. GS_MATURE_TOL)
	  END IF
!
! Same limit as the warm-up.  Under DO_LAM_IT=F an immature field
! does not rewrite a scheduled BA step into Lambda; the mode
! controller inserts the Lambda iterations after the step, once
! MAXCH is known.
!
	  IF(RD_LAMBDA .AND. .NOT. LAMBDA_ITERATION)THEN
	    IF(GS_UNSAFE .OR. GS_FG_IMMATURE .OR.
	1      (MOM_RAN .AND. MOM_UNSAFE_B))THEN
	      IF(CONV_DIAG_LVL .GE. 1)THEN
	      WRITE(LUER,'(A)')
	1     'SOLVE_FOR_POPS: KNOP GS state immature'
	      WRITE(LUER,'(A,ES12.4,A,I8,A,I14,A,ES12.4,A,I8,A,I8)')
	1     '    FG: maxchg=',GS_MAXCHG,' capped=',GS_CAPPED,
	2     ' nconv50=',GS_NCONV,
	3     '  MOM: maxchg=',MOM_MAXCHG,' capped=',MOM_CAPPED,
	4     ' negJ=',MOM_NEGJ
	      WRITE(LUER,'(A)')'    forcing LAMBDA for this iteration'
	      END IF
	      LAMBDA_ITERATION=.TRUE.
	      FIXED_T=.TRUE.
	      COMPUTE_BA=.TRUE.
	      TEMP_CHAR='DIAG'
	      T1=MAX_LAM_COR
	      FORCE_LAMBDA=.TRUE.
	      CALL APPEND_ITER_REASON(ITER_REASON,'gs_immature_demote')
	    ELSE IF(DO_BA_RES_ACC .AND. BA_SNAPSHOT_MISSING)THEN
!
! Fail-closed resacc: DO_BA_RES_ACC promises rollback on residual
! growth, but the POPS_BA_SAVE arm failed (alloc failure under a
! memory cap).  A BA step taken now could not be rolled back -- the
! exact unchecked-step failure this gate exists for (LAM_dbg
! it39-48: silent alloc failure left every FULL step unguarded; one
! committed correction wrecked pops and field past recovery).  Take
! a bounded LAMBDA update instead, at the near-mature cap rather
! than MAX_LAM_COR -- under persistent demotion a factor-100 step
! would recreate the wreck the freeze/release machinery prevents.
! SOLVEBA retries the snapshot on demoted calls, so the BA path
! resumes automatically if the alloc succeeds.
!
	      IF(CONV_DIAG_LVL .GE. 1)WRITE(LUER,'(A)')
	1     'SOLVE_FOR_POPS: resacc snapshot unavailable'//
	2     ' -- demoting unchecked BA step to damped LAMBDA'
	      LAMBDA_ITERATION=.TRUE.
	      FIXED_T=.TRUE.
	      COMPUTE_BA=.TRUE.
	      TEMP_CHAR='DIAG'
	      T1=GS_LAM_CAP
	      FORCE_LAMBDA=.TRUE.
	      CALL APPEND_ITER_REASON(ITER_REASON,'ba_no_snapshot_demote')
	    ELSE IF(LUM_BLOCKED)THEN
!
! Luminosity runaway (LUM_DIV_*): the last steps drove the outer
! luminosity nonphysical; a BA or full-MAX_LAM step against this
! state can only amplify it.  Same recipe as the snapshot demote --
! damped LAMBDA at the near-mature cap so the field can settle while
! the detector below holds populations.
!
	      IF(CONV_DIAG_LVL .GE. 1)WRITE(LUER,'(A)')
	1     'SOLVE_FOR_POPS: luminosity runaway'//
	2     ' -- demoting BA step to damped LAMBDA'
	      LAMBDA_ITERATION=.TRUE.
	      FIXED_T=.TRUE.
	      COMPUTE_BA=.TRUE.
	      TEMP_CHAR='DIAG'
	      T1=GS_LAM_CAP
	      FORCE_LAMBDA=.TRUE.
	      CALL APPEND_ITER_REASON(ITER_REASON,'lum_div_demote')
	    END IF
	  END IF
!
! Knop Paper-II nested iteration: the wavelength boundary-value
! problem should be converged before the populations take a BA step
! (Paper II, Sect. 3.2).  When the just-completed sweep still shows an
! immature lagged coupling field, this call runs a forced-LAMBDA
! iteration instead of a BA step -- SOLVEBA is never skipped outright.
! Populations still
! get a conservative Lambda update and the standard SOLVEBA
! diagnostics (max % changes, DEC_VEC/INC_VEC, MINSCALE) print every
! iteration; the BA update resumes automatically on the first sweep
! that returns a mature field.
!
!
! W1 (17-Oct-2026): the capped-write veto is now share-gated
! (GS_PIN: >10% events is UNSAFE, >0.1% pins the marginal tier) --
! a single clamp/NaN event in one cell used to raw-immature every
! later iteration and run the freeze clock to GS_FREEZE_LIM.
! W2 adds the diverged-ray share (>1% of solved rays) to UNSAFE.
!
	  GS_RAW_IMMATURE=(GS_PIN .OR.
	1      GS_FG_IMMATURE .OR.
	2      (MOM_RAN .AND. MOM_FREEZE_B))
! GS_RAW_MATURE (VADAT, default F): a bulk-pass sweep with the raw
! max change still above GS_MATURE_TOL does not count toward release.
! The bulk exception tolerates a small fraction of drifting cells,
! but a max-norm of ~5E-2 still released on KNOP_jfnkrel it34 and the
! first unheld step detonated the lagged field.  Reset rather than
! merely hold the count: the GS_REL_STABLE test is consecutive.
	  IF(GS_RAW_IMMATURE .OR.
	1    (GS_RAW_MATURE .AND. GS_MAXCHG .GT. GS_MATURE_TOL))THEN
	    N_MATURE=0
	  ELSE
	    N_MATURE=N_MATURE+1
	  END IF
!
! Stability hold (GS_REL_STABLE, default 1): a single mature-looking
! sweep is not proof the lagged field is a fixed point -- a churning
! field can pass the bulk gate transiently, and the live update it
! releases then re-poisons the field (LAM_dbg it7 release cycle,
! twice).  Require the gate to pass on GS_REL_STABLE consecutive
! calls; the extra held sweeps still run the forced-LAMBDA path so
! the field keeps converging at fixed pops.
!
	  GS_FREEZE=GS_RAW_IMMATURE .OR. (N_MATURE .LT. GS_REL_STABLE)
	  IF(GS_FREEZE .AND. .NOT. GS_RAW_IMMATURE .AND.
	1     N_MATURE .LT. GS_REL_STABLE)
	2   CALL APPEND_ITER_REASON(ITER_REASON,'rel_stable_hold')
!
! GS_RAW_MATURE nudge (17-Oct-2026): the KNOP auto-set leaves
! GS_RAW_MATURE=F, so a bulk-pass sweep with raw max still in the
! 0.05-0.5 marginal band can release while the field is still
! moving (it34-class release-then-blow-up).  One-time warning;
! the default is unchanged -- the marginal-release machinery
! (W3 best-pass census, W8 FG_TRIAL) is the remediation, not a
! default change.
!
	  IF(.NOT. GS_RAW_NUDGE_DONE .AND. .NOT. GS_RAW_MATURE .AND.
	1     DO_NONMONO_KNOP .AND. .NOT. GS_RAW_IMMATURE .AND.
	2     N_MATURE .GE. GS_REL_STABLE .AND.
	3     GS_MAX_ALL .GT. 0.05_LDP .AND.
	4     GS_MAX_ALL .LT. 0.5_LDP)THEN
	    WRITE(LUER,'(A,ES10.3)')
	1   'SOLVE_FOR_POPS: NUDGE -- lagged field releasing with raw'//
	2   ' max still',GS_MAX_ALL
	    WRITE(LUER,'(A)')
	1   '    in the 0.05-0.5 marginal band under'//
	2   ' GS_RAW_MATURE=F (KNOP auto-set default); consider'//
	3   ' setting GS_RAW_MATURE=T if a released step detonates.'
	    GS_RAW_NUDGE_DONE=.TRUE.
	  END IF
	  IF(GS_FREEZE)THEN
	    N_FROZEN=N_FROZEN+1
	  ELSE
	    N_FROZEN=0
	    N_TRIAL_ACC=0
	  END IF
	  IF(GS_FREEZE)THEN
	    IF(CONV_DIAG_LVL .GE. 1)THEN
	    WRITE(LUER,'(/,A)')' -------------------------------- KNOP GS SWEEP --------------------------------'
	    WRITE(LUER,'(A)')'  KNOP inner-GS sweep'
	    WRITE(LUER,'(A,I8)')'    immature-field iteration : ',N_FROZEN
	    WRITE(LUER,'(A,2ES12.4)')
	1   '    FG / MOM maxchg          : ',GS_MAXCHG,MOM_MAXCHG
!
! edited by Devin, 17-Oct-2026: the sweep drift can reach
! 1E263 on a detonating sweep -- ES12.4 drops the 'E' for a
! 3-digit exponent (printed "1.7030+263"), so this line uses
! ES13.4.  The nbad line below passes the ' / ' literal: the
! previous (A,I14,A,I14) form had no argument for the second
! A descriptor and printed GS_NBAD_BLK's raw bytes.
!
	    WRITE(LUER,'(A,2ES13.4)')
	1   '    FG sweep drift / RF delta: ',GS_SW_MAXCHG,GS_RF_DCHG
	    WRITE(LUER,'(A,I6,A,I6)')
	1   '    FG div rays              : ',GS_NDIV,' of',GS_NRAY
	    WRITE(LUER,'(A,I14,A,I14)')
	1   '    FG nbad last / best      : ',GS_NBAD,' / ',GS_NBAD_BLK
	    WRITE(LUER,'(A)')'    MOM bulk'
	    WRITE(LUER,'(A,2I12)')
	1   '      nbad5 / nbad50         : ',MOM_NBAD5,MOM_NBAD50
	    WRITE(LUER,'(A,2I12)')
	1   '      negJ / quar            : ',MOM_NEGJ,MOM_NQUAR
	    WRITE(LUER,'(A,2I12)')
	1   '      ntot / pquar           : ',MOM_NTOT,MOM_NPQUAR
	    WRITE(LUER,'(A)')'    MOM sweep'
	    WRITE(LUER,'(A,ES12.4)')
	1   '      maxchg                 : ',MOM_SW_MAXCHG
	    WRITE(LUER,'(A,2I12)')
	1   '      nbad5 / nbad50         : ',MOM_SW_NBAD5,MOM_SW_NBAD50
	    WRITE(LUER,'(A,2I12)')
	1   '      negJ / quar            : ',MOM_SW_NEGJ,MOM_SW_NQUAR
	    IF(MOM_RF_RAN)THEN
	      WRITE(LUER,'(A)')'    MOM refresh'
	      WRITE(LUER,'(A,L2)')
	1     '      converged              : ',MOM_RF_CONV
	      WRITE(LUER,'(A,ES12.4)')
	1     '      retained delta         : ',MOM_RF_DCHG
	      WRITE(LUER,'(A,2I12)')
	1     '      n5 / n50               : ',MOM_RF_N5,MOM_RF_N50
	    END IF
	    END IF
!
! edited by Devin, 30-Sep-2026: marginal-release BA trial.  When the
! freeze is driven only by the marginal MOM drift terms the field
! has settled except for a small asymptoting subset that Lambda
! alone cannot converge.  Releasing the iteration lets the mode
! controller take a BA step with DO_BA_RES_ACC armed: residual
! growth rolls the step back and demotes to LAMBDA next call, so
! the trial is self-throttling.  The release requires the resacc
! net -- an unchecked BA step on a marginal field is the failure
! this gate exists to prevent.  N_FROZEN still counts the iteration
! as frozen, so KNOP_STALL_MAX / GS_FREEZE_LIM still bound a run
! whose field never matures.
!
! The trial is resacc-gated: an armed rollback snapshot is part of
! the promise, so a missing snapshot vetoes the trial the same way a
! missing resacc flag does.
!
	    MOM_TRIAL=(MOM_RAN .AND. (.NOT. MOM_UNSAFE_B) .AND.
	1      (.NOT. MOM_BROAD) .AND.
	2      DO_BA_RES_ACC .AND. (.NOT. BA_SNAPSHOT_MISSING) .AND.
	3      (MOM_TRIAL_MIN .GT. 0) .AND.
	4      (N_FROZEN .GE. MOM_TRIAL_MIN) .AND.
	5      (.NOT. GS_UNSAFE) .AND. (.NOT. GS_FG_IMMATURE) .AND.
	6      (.NOT. LUM_BLOCKED) .AND.
	7      (GS_TRIAL_LIM .LE. 0 .OR. N_TRIAL_ACC .LT. GS_TRIAL_LIM))
	    IF(CONV_DIAG_LVL .GE. 2)
	1   WRITE(LUER,'(A,9L2,A,I6,A,I10,A,ES12.4,5(A,I10),A,I4)')
	2   '      MOM_TRIAL terms: ran,capok,nhuge,nb50ok,'//
	3   'njok,nqok,nbroad,snap,trlim=',
	4   MOM_RAN,MOM_SW_CAPPED .LE. MAX(50_8,MOM_SW_NTOT/10),
	5   MOM_SW_NHUGE .LE. MAX(50_8,MOM_SW_NTOT/10000),
	6   MOM_SW_NBAD50 .LE. MAX(50_8,MOM_SW_NTOT/100),
	7   MOM_SW_NEGJ .LE. MOM_SW_NTOT/10,
	8   MOM_SW_NQUAR .LE. MOM_SW_NTOT/10,
	9   .NOT. MOM_BROAD,
	1   DO_BA_RES_ACC .AND. (.NOT. BA_SNAPSHOT_MISSING),
	1   (GS_TRIAL_LIM .LE. 0 .OR. N_TRIAL_ACC .LT. GS_TRIAL_LIM),
	1   '  nfroz=',N_FROZEN,
	2   '  cap=',MOM_SW_CAPPED,
	3   '  mx=',MOM_SW_MAXCHG,
	4   '  n50=',MOM_SW_NBAD50,
	5   '  ntot=',MOM_SW_NTOT,
	6   '  negj=',MOM_SW_NEGJ,
	7   '  nquar=',MOM_SW_NQUAR,
	8   '  nhuge=',MOM_SW_NHUGE,
	9   '  ntr=',N_TRIAL_ACC
!
! W8 (17-Oct-2026): FG-only marginal-release arm ([FG_TRIAL_MIN],
! default 30, 0 disables).  Mirrors MOM_TRIAL for the
! USE_LAM_ES=T path where MOM never ran: arms only when the cured
! gate sees no immaturity and no UNSAFE census (a PIN-tier cap
! share or the GS_REL_STABLE hold is the marginal condition the
! trial exists to escape), the resacc snapshot is armed, the
! diverged-ray share is inside the W2 bound (folded into
! GS_UNSAFE), luminosity is not blocked, and the accepted-trial
! pulse budget remains.
!
	    FG_TRIAL=(.NOT. MOM_RAN) .AND. (.NOT. GS_FG_IMMATURE)
	1      .AND. (.NOT. GS_UNSAFE)
	2      .AND. DO_BA_RES_ACC .AND. (.NOT. BA_SNAPSHOT_MISSING)
	3      .AND. (FG_TRIAL_MIN .GT. 0)
	4      .AND. (N_FROZEN .GE. FG_TRIAL_MIN)
	5      .AND. (.NOT. LUM_BLOCKED)
	6      .AND. (GS_TRIAL_LIM .LE. 0 .OR.
	7      N_TRIAL_ACC .LT. GS_TRIAL_LIM)
	    IF(MOM_TRIAL)THEN
	      TRIAL_FIRED=.TRUE.
	      TRIAL_MIN_USED=MOM_TRIAL_MIN
	      IF(MOM_TRIAL_LAM .GT. 0.0_LDP .AND.
	1      MAXCH .GT. MOM_TRIAL_LAM)THEN
		IF(CONV_DIAG_LVL .GE. 1)
	1       WRITE(LUER,'(A,ES10.3)')
	2       '      MOM field marginal: demand beyond BA'//
	3       ' validity -- releasing damped LAMBDA trial,'//
	4       ' maxch=',MAXCH
!
! At 1E7-scale demand the resacc-gated BA step below is rejected
! every time (the linearization extrapolates orders of magnitude
! past validity; m15 KNOP warm looped release->blow-up->rollback).
! The factor-3 LAMBDA probe is the same contractive step the
! schedule takes on LAMBDA-mode releases; on this model it cut
! the SE residual ~6x in one call.  No resacc arm rides on a
! LAMBDA call, but the step cannot run away: its correction cap
! is the near-mature limit the freeze path already uses.
!
		LAMBDA_ITERATION=.TRUE.
		FIXED_T=.TRUE.
		COMPUTE_BA=.TRUE.
		TEMP_CHAR='DIAG'
		T1=GS_LAM_CAP
		CALL APPEND_ITER_REASON(ITER_REASON,'mom_trial_lam')
	      ELSE
		IF(CONV_DIAG_LVL .GE. 1)
	1       WRITE(LUER,'(A)')'      MOM field marginal: releasing'//
	2       ' resacc-gated BA trial this iteration'
!
! The SOLVEBA call below uses the ENTRY LAMBDA_ITERATION/T1 -- the
! mode controller further down only sets the mode for the NEXT
! iteration.  Without these assignments a "released" trial would
! still take a LAMBDA solve (LAMBDA_IT=T keeps resacc disarmed) but
! with the entry T1=MAX_LAM_COR cap instead of GS_LAM_CAP.
!
		LAMBDA_ITERATION=.FALSE.
		FIXED_T=.TRUE.
		COMPUTE_BA=.TRUE.
		TEMP_CHAR=METH_SOL
		T1=MAX_LIN_COR
		CALL APPEND_ITER_REASON(ITER_REASON,'mom_trial_release')
	      END IF
	    ELSE IF(FG_TRIAL)THEN
!
! W8 (17-Oct-2026, [FG_TRIAL_MIN], default 30, 0 disables): the
! FG-only marginal-release arm.  Under USE_LAM_ES=T MOM_RAN is
! false, so MOM_TRIAL can never fire and a field frozen only by
! the marginal tier (PIN-tier cap share, or the stability hold)
! used to ride GS_FREEZE_LIM to a non-convergent stop.  The arm
! requires the W1-W3-cured gate to see no immaturity and no
! UNSAFE-tier census, plus the resacc net.  The trial is the
! damped LAMBDA probe (T1=GS_LAM_CAP), resacc-judged by the same
! accepted/rejected refund machinery as MOM_TRIAL.
!
	      TRIAL_FIRED=.TRUE.
	      TRIAL_MIN_USED=FG_TRIAL_MIN
	      IF(CONV_DIAG_LVL .GE. 1)
	1     WRITE(LUER,'(A,I6)')
	2     '      FG field marginal: releasing damped LAMBDA'//
	3     ' trial, n_frozen=',N_FROZEN
	      LAMBDA_ITERATION=.TRUE.
	      FIXED_T=.TRUE.
	      COMPUTE_BA=.TRUE.
	      TEMP_CHAR='DIAG'
	      T1=GS_LAM_CAP
	      CALL APPEND_ITER_REASON(ITER_REASON,'fg_trial_lam')
	    ELSE
	    LAMBDA_ITERATION=.TRUE.
	    FIXED_T=.TRUE.
	    COMPUTE_BA=.TRUE.
	    TEMP_CHAR='DIAG'
	    GS_CHURN=0.0_LDP
	    IF(GS_NTOT .GT. 0_8)THEN
	      GS_CHURN=DBLE(GS_NCONV)/DBLE(GS_NTOT)
	    END IF
	    IF(MOM_SW_NTOT .GT. 0_8)THEN
	      GS_CHURN=MAX(GS_CHURN,DBLE(MOM_SW_NBAD50)/DBLE(MOM_SW_NTOT))
	    END IF
	    IF(GS_CHURN .LT. GS_LAM_TRANSIENT)THEN
	      T1=GS_LAM_CAP
	    ELSE
	      T1=MAX_LAM_COR
	    END IF
!
! The churn pick above can hand back MAX_LAM_COR while a
! luminosity runaway is armed -- silently reversing the
! LUM_BLOCKED demote (lum_div_demote takes GS_LAM_CAP).
! Respect the demote regardless of the churn reading.
!
	    IF(LUM_BLOCKED)T1=MIN(T1,GS_LAM_CAP)
	    IF(CONV_DIAG_LVL .GE. 1)THEN
	      WRITE(LUER,'(A)')'    LAMBDA correction cap'
	      WRITE(LUER,'(A,ES12.4)')
	1     '      churn fraction         : ',GS_CHURN
	      WRITE(LUER,'(A,ES12.4)')
	1     '      max correction         : ',T1
	    END IF
	    FORCE_LAMBDA=.TRUE.
	    CALL APPEND_ITER_REASON(ITER_REASON,'gs_freeze_forced')
	    END IF
	  END IF
!
! Stock cycle (DO_LAM_IT=F).  A scheduled BA step is kept, including
! a large first correction: the freeze above would have rewritten it
! as Lambda.  A scheduled Lambda step is not rewritten, but while the
! lagged I_NEXT field is immature its correction stays at GS_LAM_CAP.
! Clearing that cap (m15 JFNK seq, it2) let a factor-1e5 update
! through and the next sweep's luminosity rose from 2e10 to 2e16.
! FORCE_LAMBDA is cleared so NUM_LAM can still return to BA.
! DO_LAM_IT=T leaves the freeze as written above.
!
	  IF(.NOT. RD_LAMBDA .AND. GS_FREEZE)THEN
	    IF(.NOT. LAM_AT_ENTRY)THEN
	      LAMBDA_ITERATION=.FALSE.
	      FIXED_T=RD_FIX_T
	      COMPUTE_BA=.TRUE.
	      TEMP_CHAR=METH_SOL
	      T1=MAX_LIN_COR
	      FORCE_LAMBDA=.FALSE.
	      IF(CONV_DIAG_LVL .GE. 1)WRITE(LUER,'(A)')
	1       '  stock cycle: BA step kept; immature field does not replace it'
	    ELSE
	      T1=MIN(T1,GS_LAM_CAP)
	      FORCE_LAMBDA=.FALSE.
	      IF(CONV_DIAG_LVL .GE. 1)WRITE(LUER,'(A,ES12.4)')
	1       '  stock cycle: Lambda step capped while I_NEXT is immature, factor=',
	2       T1
	    END IF
	  END IF
!
!
! Release damping (GS_REL_N, default 1): the first live LAMBDA
! update after a frozen stretch carries the deferred correction
! accumulated while pops were held -- under MAX_LAM_COR that step
! re-poisons the field the freeze just converged (LAM_dbg it7,
! twice) -- and the next few steps still see a field re-lagging to
! the moved populations (the undamped it8 step committed another
! factor-100 and re-froze the run).  Keep GS_LAM_CAP, the
! near-mature limit the freeze path already uses once churn is
! below GS_LAM_TRANSIENT, for the whole GS_REL_N-iteration window.
! MOM_TRIAL releases keep MAX_LIN_COR and non-LAMBDA releases are
! covered by the resacc net; neither is touched here.
!
	  IF((GS_WAS_FROZEN .OR.
	1    (N_REL_SINCE .GE. 1 .AND. N_REL_SINCE .LT. GS_REL_N))
	2    .AND. .NOT. GS_FREEZE .AND.
	3    LAMBDA_ITERATION .AND. .NOT. MOM_TRIAL)THEN
	    T1=MIN(T1,GS_LAM_CAP)
	    IF(CONV_DIAG_LVL .GE. 1)
	1   WRITE(LUER,'(A,ES12.4)')
	2   '      GS release: damping LAMBDA step, cap=',T1
	    CALL APPEND_ITER_REASON(ITER_REASON,'release_damped')
	  END IF
!
! edited by Devin, 29-Sep-2026 (moved post-solve 30-Sep-2026): the
! GS_STATE write lives below the SOLVEBA call so the extra keys
! report the just-completed solve rather than last call's values.
! submit_cmfgen.sh greps '^n_frozen=' with -m1, so n_frozen must
! stay the first key in the file.
!
! Optional in-model cap (VADAT [GS_FREEZE_LIM], default 0 = off):
! the field staying immature this long is a structural failure, not
! a transient -- stop loudly rather than consume the full requeue
! budget doing forced-LAMBDA iterations that cannot converge.
!
	  IF(GS_FREEZE_LIM .GT. 0 .AND. N_FROZEN .GE. GS_FREEZE_LIM)THEN
	    WRITE(LUER,'(A,I4,A)')'SOLVE_FOR_POPS: ',N_FROZEN,
	1    ' consecutive GS-frozen iterations reach GS_FREEZE_LIM'
	    WRITE(LUER,'(A)')'  the lagged GS field never matured;'//
	1    ' stopping as non-convergent'
	    STOP
	  END IF
!
! Carry the freeze state and the post-release window counter into
! the next call so released steps can be damped above.
!
	  IF(GS_FREEZE)THEN
	    N_REL_SINCE=0
	  ELSE IF(GS_WAS_FROZEN .OR. N_REL_SINCE .GT. 0)THEN
	    N_REL_SINCE=N_REL_SINCE+1
	  END IF
	  GS_WAS_FROZEN=GS_FREEZE
	END IF
!
! FIX_IMPURITY is a BA-mode option (upstream convention: only used in
! non-LAMBDA mode, cleared whenever a LAMBDA iteration is imposed).
! The KNOP gates above can force LAMBDA after a BA phase has already
! set FIX_IMPURITY=T; without this reset the forced-LAMBDA sweep would
! freeze trace ions instead of relaxing them.  MOD_FIX_IMPURITY was
! snapshotted before the gates and must be cleared directly.
!
! Same snapshot-ordering issue for MOD_FIXED_T: the gates set
! FIXED_T=.TRUE. but MOD_FIXED_T/MOD_FIX_T_D_ST/MOD_FIX_T_D_END were
! computed above from the incoming FIXED_T/RD_FIX_T.  With
! MOD_FIXED_T=F the R.E. row is not held fixed -- the LAMBDA solve
! then hits a degenerate T equation (singular BA -> SOLVEBA STOP) or
! a garbage T correction.  Force the full-depth fix the gates intend.
!
	IF(FORCE_LAMBDA)THEN
	  FIX_IMPURITY=.FALSE.
	  MOD_FIX_IMPURITY=.FALSE.
	  MOD_FIXED_T=.TRUE.
	  MOD_FIX_T_D_ST=1
	  MOD_FIX_T_D_END=ND
	END IF
!
! edited by Devin, 30-Sep-2026: GS_HOLD_POPS (VADAT, default F)
! turns the GS-freeze iteration into a true transfer-only sweep --
! the lagged coupling field is converged at fixed populations
! (Knop Paper II nested iteration) rather than against a target
! that moves inside the same iteration.  The forced-LAMBDA SOLVEBA
! still runs so the STEQ diagnostics show the attempted update;
! the committed correction is rolled back below.  The MOM_TRIAL
! release path is excluded: a released trial exists to take a real
! resacc-gated BA step.
!
! edited by Devin, 17-Oct-2026: LUM_POPS_HELD is no longer dead --
! a detected luminosity runaway (LUM_BLOCKED) holds populations so
! the sweeps can only move the lagged field, matching the
! GS_HOLD_POPS recipe of converging at fixed pops.
!
	LUM_POPS_HELD=LUM_BLOCKED
	GS_POPS_HELD=.FALSE.
	IF((GS_HOLD_POPS .AND. GS_FREEZE .AND. .NOT. MOM_TRIAL) .OR.
	1   LUM_POPS_HELD)THEN
	  IF(.NOT. ALLOCATED(POPS_GS_SAVE))THEN
	    ALLOCATE(POPS_GS_SAVE(NT,ND),STAT=IOS)
	  ELSE IF(SIZE(POPS_GS_SAVE,1) .NE. NT .OR.
	1          SIZE(POPS_GS_SAVE,2) .NE. ND)THEN
	    DEALLOCATE(POPS_GS_SAVE)
	    ALLOCATE(POPS_GS_SAVE(NT,ND),STAT=IOS)
	  ELSE
!
! No alloc ran: IOS still holds the last I/O status from the previous
! call (ITER_HISTORY append, GS_STATE read at restart), and a stale
! nonzero would silently skip the hold while reporting an alloc
! failure that never happened.  Same fix as the POPS_BA_SAVE block
! in SOLVEBA_V13.
!
	    IOS=0
	  END IF
	  IF(IOS .EQ. 0)THEN
	    POPS_GS_SAVE=POPS
	    GS_POPS_HELD=.TRUE.
	  ELSE
	    WRITE(LUER,*)'Population hold snapshot allocation failed; stopping before update'
	    STOP
	  END IF
	END IF
!
!
! USE_JFNK_SOLVER (VADAT [DO_JFNK], default F; [USE_JFNK] deprecated
! alias): route the correction solve through the Jacobian-Free
! Newton-Krylov solver instead of the direct BA factorization.  Only
! fires on a full (non-LAMBDA) call after JFNK_WARMUP calls, so the
! LAMBDA/frozen-field warmup and all GS/resacc gating above are
! untouched.  JFNK_SOLVE_CORR_V1 falls back internally to a
! LAMBDA-mode SOLVEBA_V13 on failure.  Admitted under DO_NONMONO_KNOP
! or DO_NONMONO_SOB; the coupled tier is KNOP-only (JFNK_COUPLED is
! rejected under SOB at VADAT read), so SOB always runs the frozen-J
! residual here.
!
! JFNK_ANY_ITER (VADAT [JFNK_ANYITER], default F): bring-up/debug
! only -- fire the JFNK solve even on a LAMBDA-designated call.  The
! commit/hold logic below still applies, so a hold_pops iteration
! discards the result harmlessly.
!
! Newton is not started while the lagged I_NEXT field is immature.
! The seq run called JFNK at it9 with |F| ~ 1e20 and F(0) rel=44,
! then the failure path took another MAX_LAM step.  SOLVEBA still
! runs; JFNK waits until GS_FG_IMMATURE clears.
!
! edited by Devin, 03-Oct-2026: the field test now keys on
! GS_FREEZE rather than the raw GS_FG_IMMATURE flag.  GS_FREEZE
! subsumes it and additionally holds Newton (a) through the
! N_MATURE/GS_REL_STABLE consecutive-mature-sweep stability window
! and (b) on MOM_TRIAL marginal-field releases, where the machinery
! deliberately permits only a resacc-bounded BA probe -- an
! unbounded Newton commit there defeats the trial's purpose.
!
! edited by Devin, 05-Oct-2026: admission also requires USE_LAM_ES.
! The eval-time lagged-field freeze covers only the Knop I_NEXT
! store (FREEZE_I_NEXT, armed in JFNK_SOLVE_CORR_V1, honoured in
! FG_J_CMF_NONMONO_V1).  Under USE_LAM_ES=F the J update runs
! MOM_J_CMF_NONMONO_V1, whose J_NEXT_FRQ/H_NEXT_FRQ lagged stores
! are written in place on every sweep with no freeze -- each JFNK
! residual/JVP eval would overwrite the lagged field the committed
! iterate relies on, corrupting the in-sweep coupling.
!
! edited by Devin, 17-Oct-2026: under DO_NONMONO_SOB only the frozen
! tier exists (JFNK_COUPLED is a hard STOP in the reader), and the
! frozen-J residual is the same fixed-J system SOLVEBA solves
! directly -- JFNK can therefore never accelerate a radiation-pinned
! residual (m15 SOB_l14_jfnkrel: 186 admits, 0 GMRES, every one a
! reach-gate skip).  For that model class consider disabling the
! relief/DO_JFNK, or cap the wasted admissions via
! JFNK_REL_REACH_MAX.
!
	JFNK_THIS=USE_JFNK_SOLVER .AND. ITS_SINCE_START .GT. JFNK_WARMUP
	1          .AND. (.NOT. LAMBDA_ITERATION .OR. JFNK_ANY_ITER)
	2          .AND. (.NOT. DO_NONMONO_KNOP .OR. .NOT. GS_FREEZE)
	3          .AND. USE_LAM_ES
	4          .AND. .NOT. RDLAM_PROBE_THIS
!
! JFNK_COR_GATE (VADAT [JFNK_COR_GATE], default 0=off): do not start
! Newton while the demanded correction is still large.  At entry
! MAXCH holds the previous call's raw demanded maximum correction --
! a saturated demand (1e7 sentinel or ~1e10) means the iterate is
! far outside the basin where the finite-difference JVP is usable
! (m15 JFNK_opt it11: demand ~1e10%, |F|_w=2.7e13, GMRES never got
! past its first matvec, and the unprotected MAX_LIN_COR safety
! step then pushed se_res to 2.5e23).  Deferred calls take the
! normal SOLVEBA path, which keeps resacc rollback protection on
! full steps.  JFNK stays armed and self-admits once the staged
! LAMBDA/full map brings the demand under the gate.
!
! Two-tier admission (audit cmfgen_obsidian/audits/JFNK_SOLVER_AUDIT, 03-Oct-2026): a
! coupled eval costs a full formal sweep (~250-400 s on m15), so the
! coupled residual is held for the tighter COR_GATE tier while the
! cheaper frozen-J Newton (no re-solve per eval; BA is the near-exact
! Jacobian of the frozen residual) is admitted at the looser
! JFNK_FRZ_GATE.  With FRZ_GATE=0 the frozen tier uses COR_GATE and
! behavior reduces to the original single gate.  JFNK_GATE_NREQ adds
! hysteresis: the demand must stay sub-gate for NREQ consecutive
! solves, so one lucky low-MAXCH step cannot fire Newton on a field
! that is still drifting.  Counter semantics: MAXCH at entry is the
! previous call's demanded correction, so the hits count how long
! the state has remained inside the basin.
!
	JFNK_GATE_IN=JFNK_FRZ_GATE
	IF(JFNK_GATE_IN .LE. 0.0_LDP)JFNK_GATE_IN=JFNK_COR_GATE
!
! edited by Devin, 03-Oct-2026: hysteresis counters advance only on
! calls that are otherwise JFNK-eligible.  During warmup,
! LAMBDA-designated, GS-frozen, or USE_LAM_ES=F calls the demanded
! correction is measured against a field that is not the operator
! the Newton solve would see (and under freeze the deferred SE
! correction keeps growing, and under USE_LAM_ES=F JFNK can never
! be admitted at all), so churn-era demand must not pre-arm the
! gate.
!
	JFNK_DEM_IN=MAXCH
	IF(USE_JFNK_SOLVER .AND.
	1     ITS_SINCE_START .GT. JFNK_WARMUP .AND.
	2     (.NOT. LAMBDA_ITERATION .OR. JFNK_ANY_ITER) .AND.
	3     (.NOT. DO_NONMONO_KNOP .OR. .NOT. GS_FREEZE) .AND.
	4     USE_LAM_ES .AND. .NOT. RDLAM_PROBE_THIS)THEN
!
! Demand-contraction term (edited by Devin, 03-Oct-2026): a sub-
! gate demand is only a basin signal while the demanded correction
! is not rebounding.  MAXCH here is the previous call's demand;
! JFNK_LAST_DEM is the one before it, so the ratio measures the
! outer map's actual contraction.  A rebounding demand -- field
! re-poisoned by the last committed step, or the post-Lambda-exit
! transient where the first full solve re-demands a large
! correction -- resets the count even below the gate, so Newton
! fires only on a sustained, contracting, sub-gate trajectory.
! First eligible call after a NEWMOD/ineligible stretch has no
! predecessor and tests the level only.
!
! edited by Devin, 04-Oct-2026: gate on the last non-LAMBDA demand
! (BA_LAST_RAW), not MAXCH.  MAXCH here is the previous call's
! demand; in the alternating regime that call is always a LAMBDA
! solve whose small demand reflects the LAMBDA attractor, not the
! operator Newton would see (opt run admitted JFNK at it73 on a
! 95% LAMBDA demand while the last full demand was 4.8E4% and
! |F|=2.2E12 vs target 2.2E8).  BA_LAST_RAW refreshes only at the
! end of non-LAMBDA calls, so the level and contraction tests now
! measure the same operator.  Before the first full call
! BA_LAST_RAW is 0; fall back to MAXCH there.
!
	  JFNK_DEM_IN=BA_LAST_RAW
	  IF(JFNK_DEM_IN .LE. 0.0_LDP)JFNK_DEM_IN=MAXCH
	  JFNK_DEM_CTR=JFNK_LAST_DEM .LE. 0.0_LDP .OR.
	1      JFNK_DEM_IN .LE. JFNK_DEM_RISE*JFNK_LAST_DEM
	  IF((DO_NONMONO_SOB .OR. DO_NONMONO_KNOP) .AND.
	1   (BA_RAW_REJECTED .OR. BA_RAW_TRIAL_PENDING))
	1      JFNK_DEM_CTR=.FALSE.
	  IF(JFNK_GATE_IN .GT. 0.0_LDP .AND.
	1       JFNK_DEM_IN .LT. JFNK_GATE_IN .AND. JFNK_DEM_CTR)THEN
	    JFNK_FRZ_HITS=JFNK_FRZ_HITS+1
	  ELSE
	    JFNK_FRZ_HITS=0
	  END IF
	  IF(JFNK_COR_GATE .GT. 0.0_LDP .AND.
	1       JFNK_DEM_IN .LT. JFNK_COR_GATE .AND. JFNK_DEM_CTR)THEN
	    JFNK_CPL_HITS=JFNK_CPL_HITS+1
	  ELSE
	    JFNK_CPL_HITS=0
	  END IF
	  JFNK_LAST_DEM=JFNK_DEM_IN
	ELSE
	  JFNK_FRZ_HITS=0
	  JFNK_CPL_HITS=0
	  JFNK_LAST_DEM=0.0_LDP
	END IF
!
	JFNK_COUPLED_THIS=.FALSE.
	IF(JFNK_THIS)THEN
!
! JFNK_CPL_RELWAIT (VADAT, default F; edited by Devin, 17-Oct-2026):
! hold the coupled tier through the post-release damp window.  The
! first live sweeps after a GS freeze carry the deferred SE
! correction and a field re-lagging to moved populations -- admitting
! JFNK-CPL there couples a marginal field to a perturbative solve
! (KNOP_jfnkrel it35: admitted on the release call, stagnated, and
! the LAMBDA fallback then re-poisoned the field, se_res -> 2E31).
! N_REL_SINCE was already incremented above, so 1 marks the first
! live call; the damped window spans N_REL_SINCE <= GS_REL_N.
! Frozen-tier admission below is unaffected.  NOTE: a config with a
! very large GS_REL_N defers the coupled tier indefinitely.
!
	  IF(DO_NONMONO_KNOP .AND. JFNK_CPL_RELWAIT .AND.
	1     JFNK_COUPLED .AND.
	2     (JFNK_COR_GATE .LE. 0.0_LDP .OR.
	3      JFNK_CPL_HITS .GE. JFNK_GATE_NREQ) .AND.
	4     N_REL_SINCE .GE. 1 .AND.
	5     N_REL_SINCE .LE. GS_REL_N)THEN
	    CALL APPEND_ITER_REASON(ITER_REASON,'jfnk_rel_defer')
	    WRITE(LUER,'(/,A)')' ---------------------------------- JFNK DEFER ---------------------------------'
	    WRITE(LUER,'(A)')
	1   '    reason : coupled tier held during GS release window'
	    WRITE(LUER,'(A,I6,A,I6)')
	1   '    live calls since release / window :',N_REL_SINCE,
	2   ' /',GS_REL_N
	    WRITE(LUER,'(A)')
	1   '    disposition : frozen tier / SOLVEBA path this call'
	  ELSE IF(JFNK_COUPLED .AND.
	1     (JFNK_COR_GATE .LE. 0.0_LDP .OR.
	2      JFNK_CPL_HITS .GE. JFNK_GATE_NREQ))THEN
	    JFNK_COUPLED_THIS=.TRUE.
	  END IF
	  IF(.NOT. JFNK_COUPLED_THIS)THEN
! Frozen tier: admit only when the effective gate is satisfied and
! the hysteresis count is met.  GATE_IN=0 means no gating at all
! (legacy behavior with COR_GATE=0).
!
! edited by Devin, 16-Oct-2026: hysteresis relief under
! DO_NONMONO_SOB.  In the LAMBDA-dominated regime every LAMBDA-
! designated call is JFNK-ineligible and resets JFNK_FRZ_HITS above,
! so a JFNK_GATE_NREQ>1 count is unreachable however quiet the
! demand -- the stable-exit release then runs plain SOLVEBA and JFNK
! never fires (m15 SOB_l14_jfnk OUTGEN: LAMBDA-designated defer on
! every eligible-era call).  When the LAMBDA-mode controller has
! already proven stability (MAXCH_CTRL below VAL_EXIT_LAM for
! N_STABLE_LAM consecutive calls -- the same condition that releases
! this call to full) one sub-gate, contracting full-call demand is
! enough.  The demand level and contraction tests are unchanged;
! only the multi-call hit count is relaxed.
!
	    JFNK_NREQ_EFF=JFNK_GATE_NREQ
!
! The relief cannot key on the live count: CNT_STABLE_LAM is zeroed
! by the stable-exit release that produced this eligible call, so
! CNT_STABLE_LAM >= N_STABLE_LAM was never true inside this IF.
! LAM_STABLE_PROVEN latches the earned count at the exit and stays
! set until the next instability reset (edited by Devin, 16-Oct-2026).
!
! The relief is further de-rated by reachability (17-Oct-2026): under
! SOB only the frozen tier exists, and a frozen-J Newton cannot reduce
! a radiation-pinned residual (it is the same fixed-J system SOLVEBA
! solves directly).  Once JFNK_REACH_SKIP_STREAK consecutive
! admissions have been skipped by the reach gate, the relief stands
! down and leaves the outer BA/LAMBDA loop to advance the lag.
! JFNK_REL_REACH_MAX=0 restores the original unconditional relief.
!
	    IF(DO_NONMONO_SOB .AND. LAM_STABLE_PROVEN)THEN
	      IF(JFNK_REL_REACH_MAX .GT. 0 .AND.
	1       JFNK_REACH_SKIP_STREAK .GE. JFNK_REL_REACH_MAX)THEN
	        WRITE(LUER,'(A,I4,A)')'    JFNK: stable-LAMBDA relief'//
	1       ' suppressed (reach-skip streak >=',JFNK_REL_REACH_MAX,
	2       ')'
	      ELSE
	        JFNK_NREQ_EFF=1
	        WRITE(LUER,'(A)')'    JFNK: stable-LAMBDA relief active'//
	1       ' (frozen-tier NREQ -> 1)'
	      END IF
	    END IF
	    IF(JFNK_GATE_IN .GT. 0.0_LDP .AND.
	1       (JFNK_DEM_IN .GE. JFNK_GATE_IN .OR.
	2        JFNK_FRZ_HITS .LT. JFNK_NREQ_EFF))THEN
	      JFNK_THIS=.FALSE.
	      CALL APPEND_ITER_REASON(ITER_REASON,'jfnk_defer')
	      WRITE(LUER,'(/,A)')' ---------------------------------- JFNK DEFER ---------------------------------'
	      WRITE(LUER,'(A,ES12.4)')
	1     '    last non-LAMBDA demand BA_LAST_RAW:',JFNK_DEM_IN
	      WRITE(LUER,'(A,ES12.4)')
	1     '    demanded population change MAXCH  :',MAXCH
	      WRITE(LUER,'(A,ES12.4)')
	1     '    frozen-tier admission gate        :',JFNK_GATE_IN
	      WRITE(LUER,'(A,L2)')
	1     '    demand still contracting          :',JFNK_DEM_CTR
	      WRITE(LUER,'(A,I3,A,I3)')
	1     '    consecutive sub-gate calls, need  :',
	1     JFNK_FRZ_HITS,' /',JFNK_NREQ_EFF
	      WRITE(LUER,'(A)')
	1     '    disposition : falls back to SOLVEBA path'
	    END IF
	  END IF
	END IF
!
! DO_BA_DIR_CHK (VADAT [DO_BA_DIR_CHK], default F; requires DO_JFNK):
! diagnostic-only JFNK admission on direct-BA iterations.  The demand
! gates above are bypassed -- the check exists to price the stored-BA
! direction exactly where the real BA step linearizes, which is
! precisely the high-demand first full step the gates defer.  The
! structural gates (warmup, non-LAMBDA, GS-freeze under KNOP,
! USE_LAM_ES, not a rd_lam probe) still apply: outside them the
! stored operator is not the one SOLVEBA is about to use.  Frozen-J
! tier only; JFNK_SOLVE_CORR_V1 exits through the info=-2 path after
! the direction metrics, so caller-mode SOLVEBA_V13 takes the real
! step on the same call.  The deferred call then counts toward the
! NG window like the BA iterate it is.
!
! edited by Devin, 15-Oct-2026: the [BA_DIRCHK_MAX] cap disengages
! the bypass once the diagnostic budget is spent -- after
! JFNK_DIRCHK_N defer events the admission reverts to the ordinary
! demand gates so a real Newton step can run.  -1 (default) keeps
! the uncapped behaviour.
!
	IF(DO_BA_DIR_CHK .AND. .NOT. JFNK_THIS .AND.
	1   (BA_DIRCHK_MAX .LT. 0 .OR.
	2    JFNK_DIRCHK_N .LT. BA_DIRCHK_MAX) .AND.
	3   USE_JFNK_SOLVER .AND.
	4   ITS_SINCE_START .GT. JFNK_WARMUP .AND.
	5   (.NOT. LAMBDA_ITERATION .OR. JFNK_ANY_ITER) .AND.
	6   (.NOT. DO_NONMONO_KNOP .OR. .NOT. GS_FREEZE) .AND.
	7   USE_LAM_ES .AND. .NOT. RDLAM_PROBE_THIS)THEN
	  JFNK_THIS=.TRUE.
	  JFNK_COUPLED_THIS=.FALSE.
	  CALL APPEND_ITER_REASON(ITER_REASON,'dirchk')
	  WRITE(LUER,'(/,A)')' ---------------------------------- JFNK DIRCHK --------------------------------'
	  WRITE(LUER,'(A)')
	1   '    DO_BA_DIR_CHK: direction diagnostics on the stored BA;'
	  WRITE(LUER,'(A)')
	1   '    Newton skipped -- caller-mode SOLVEBA takes the step'
	END IF
!
! Stage B is a BA-space triangular solve: SOLVEBA_V13 owns the
! POPS_OLD row restore that leaves only exempted trace rows free.
! A Newton or LAMBDA dispatch would update the full vector and undo
! the major-state freeze, so BA_TRACE_PHASE always takes the direct
! BA path.  FIX_IMPURITY is cleared for the same reason: it is the
! inverse mask and would freeze the trace block Stage B is solving.
!
	IF(BA_TRACE_PHASE)THEN
	  JFNK_THIS=.FALSE.
	  JFNK_COUPLED_THIS=.FALSE.
	  LAMBDA_ITERATION=.FALSE.
	  FIX_IMPURITY=.FALSE.
	  MOD_FIX_IMPURITY=.FALSE.
	END IF
!
! edited by Devin, 03-Oct-2026: the non-demand defer paths were
! silent -- an armed JFNK that sat out on field maturity or a
! LAMBDA designation printed nothing, which made gate verification
! from OUTGEN alone impossible.  The demand-gate deferral above
! already prints its own line; a call that reached that check had
! LAMBDA_ITERATION=F and GS_FREEZE=F by construction, so the
! branches here cannot double-report it.
!
	IF(USE_JFNK_SOLVER .AND. .NOT. JFNK_THIS)THEN
	  IF(.NOT. USE_LAM_ES)THEN
!
! Permanent hold-off: the MOM lagged store has no eval freeze (see
! the admission gate above), so this reason outranks the transient
! warmup/freeze causes below.
!
	    WRITE(LUER,'(/,A)')' ---------------------------------- JFNK DEFER ---------------------------------'
	    WRITE(LUER,'(A)')
	1   '    reason : USE_LAM_ES=F; MOM lagged J_NEXT/H_NEXT store'//
	2   ' has no JFNK eval freeze'
	  ELSE IF(ITS_SINCE_START .LE. JFNK_WARMUP)THEN
	    IF(CONV_DIAG_LVL .GE. 1)THEN
	      WRITE(LUER,'(/,A)')' ---------------------------------- JFNK DEFER ---------------------------------'
	      WRITE(LUER,'(A)')
	1     '    reason : warmup window still open'
	      WRITE(LUER,'(A,I4,A,I4)')
	1     '    eligible calls so far, needed     :',ITS_SINCE_START,
	2     ' /',JFNK_WARMUP
	    END IF
	  ELSE IF(DO_NONMONO_KNOP .AND. GS_FREEZE)THEN
	    WRITE(LUER,'(/,A)')' ---------------------------------- JFNK DEFER ---------------------------------'
	    WRITE(LUER,'(A)')
	1   '    reason : GS radiation field not ready'
	    WRITE(LUER,'(A,I3,A,I3)')
	1   '    consecutive mature sweeps, needed :',N_MATURE,
	2   ' /',GS_REL_STABLE
	    WRITE(LUER,'(A,ES12.4)')
	1   '    FG sweep max change               :',GS_MAXCHG
	    WRITE(LUER,'(A,I14,A,I14)')
	1   '    cells changed >50pct, threshold   :',MOM_SW_NBAD50,
	2   ' /',MOM_SW_NTOT
	  ELSE IF(LAMBDA_ITERATION .AND. .NOT. JFNK_ANY_ITER)THEN
	    WRITE(LUER,'(/,A)')' ---------------------------------- JFNK DEFER ---------------------------------'
	    WRITE(LUER,'(A)')
	1   '    reason : LAMBDA-designated iteration'
	  END IF
	END IF
!
! edited by Devin, 03-Oct-2026: [BA_DAMP_CAP] damped-full cap.  When
! the previous non-LAMBDA step's raw demand saturated (BA_LAST_RAW
! >= BA_SAT_RAW), every full step takes the tighter BA_DAMP_CAP
! factor bound instead of MAX_LIN_COR: the model keeps advancing
! along the BA direction but in bounded increments, so a saturated
! direction costs damped steps instead of a demote-plus-LAMBDA
! window per probe.  Enforced here -- just before the solve -- so
! the MOM_TRIAL and GS-keep-BA branches above cannot leak an
! uncapped cap through.  Telemetry-driven: BA_LAST_RAW refreshes on
! every non-LAMBDA call, so the cap releases itself once the demand
! drops below BA_SAT_RAW.  Inert unless [BA_DAMP_CAP] > 0.
!
	IF(.NOT. LAMBDA_ITERATION .AND. BA_DAMP_CAP .GT. 0.0_LDP .AND.
	1   (BA_LAST_RAW .GE. BA_SAT_RAW .OR. RDLAM_PROBE_THIS))
	2   T1=MIN(T1,BA_DAMP_CAP)
!
! Banner record of which SE solve this great iteration dispatches
! to, so the mode mix is greppable from OUTGEN without cross-
! referencing the mode-controller transitions.  Printed once per
! great iteration: additional SE calls within the same iteration
! (e.g. the T-update solve) reuse the first dispatch type.
!
	IF(MAIN_COUNTER .NE. ITER_TYPE_LAST_IT)THEN
	  ITER_TYPE_LAST_IT=MAIN_COUNTER
	  IF(JFNK_THIS)THEN
	    IF(JFNK_COUPLED_THIS)THEN
	      WRITE(ITER_MSG,'(A,I5,A)')' ITER_TYPE it=',MAIN_COUNTER,
	1       '  Running JFNK-CPL ITERATION (Newton-Krylov, J re-solved per eval)'
	    ELSE
	      WRITE(ITER_MSG,'(A,I5,A)')' ITER_TYPE it=',MAIN_COUNTER,
	1       '  Running JFNK-FRZ ITERATION (Newton-Krylov, frozen J)'
	    END IF
	  ELSE IF(LAMBDA_ITERATION)THEN
	    WRITE(ITER_MSG,'(A,I5,A)')' ITER_TYPE it=',MAIN_COUNTER,
	1     '  Running LAMBDA ITERATION (fixed-point step)'
	  ELSE
	    WRITE(ITER_MSG,'(A,I5,A)')' ITER_TYPE it=',MAIN_COUNTER,
	1     '  Running FULL ITERATION (BA linearization step)'
	  END IF
	  WRITE(LUER,'(/,A)')REPEAT('~',LEN_TRIM(ITER_MSG))
	  WRITE(LUER,'(A)')TRIM(ITER_MSG)
	  WRITE(LUER,'(A,/)')REPEAT('~',LEN_TRIM(ITER_MSG))
	END IF
!
! BA_TR_ABS/BA_DEM_EXM are exported by SOLVEBA_V13 (demanded exempted
! correction from STEQ) or, on a JFNK commit, by JFNK_SOLVE_CORR_V1
! (committed exempted change as best-effort proxy); reset before
! either solve path so a solve that exports nothing leaves no stale
! trace-delta for the Stage-B trigger in CMFGEN_SUB.
	BA_TR_ABS=0.0_LDP
	BA_DEM_EXM=0.0_LDP
!
	IF(JFNK_THIS)THEN
	  WRITE(LUER,'(/,A)')' ---------------------------------- JFNK ADMIT ---------------------------------'
	  WRITE(LUER,'(A,ES12.4)')
	1   '    demanded population change MAXCH  :',MAXCH
	  WRITE(LUER,'(A,ES12.4)')
	1   '    frozen-tier admission gate        :',JFNK_GATE_IN
	  WRITE(LUER,'(A,ES12.4)')
	1   '    coupled-tier admission gate       :',JFNK_COR_GATE
	  WRITE(LUER,'(A,I3,A,I3)')
	1   '    consecutive sub-gate calls frz/cpl:',JFNK_FRZ_HITS,
	2   ' /',JFNK_CPL_HITS
	  IF(JFNK_COUPLED_THIS)THEN
	    WRITE(LUER,'(A)')
	1   '    admitted tier : JFNK-CPL, J re-solved per eval'
	  ELSE
	    WRITE(LUER,'(A)')
	1   '    admitted tier : JFNK-FRZ, J frozen during solve'
	  END IF
	  IF(JFNK_COUPLED_THIS)THEN
	    CALL APPEND_ITER_REASON(ITER_REASON,'jfnk_cpl')
	  ELSE
	    CALL APPEND_ITER_REASON(ITER_REASON,'jfnk_frz')
	  END IF
	  CALL JFNK_SOLVE_CORR_V1(SOL,POPS,POP_ATOM,DIAG_INDX,NT,NION,
	1       NUM_BNDS,ND,MAXCH,MAXCH_APP,TEMP_CHAR,SUCCESS,SCALE_OPT,
	2       LAM_SCALE_OPT,T1,MAX_dT_COR,T_MIN,
	3       COMPUTE_BA,WR_BA_INV,WR_PART_OF_INV,LAMBDA_ITERATION,
	4       MAIN_COUNTER,SET_POPS_D2_EQ_D1,LST_ITERATION,
	5       JFNK_COUPLED_THIS)
	ELSE
	  BAACC_RDLAM_PROBE_ARM=RDLAM_PROBE_THIS
	  CALL SOLVEBA_V13(SOL,POPS,POP_ATOM,DIAG_INDX,NT,NION,NUM_BNDS,ND,
	1       MAXCH,MAXCH_APP,TEMP_CHAR,SUCCESS,SCALE_OPT,
	1       LAM_SCALE_OPT,T1,MAX_dT_COR,T_MIN,
	1       COMPUTE_BA,WR_BA_INV,WR_PART_OF_INV,LAMBDA_ITERATION,
	1       MAIN_COUNTER,SET_POPS_D2_EQ_D1)
	  BAACC_RDLAM_PROBE_ARM=.FALSE.
	END IF
!
! Population-hold rollback: the committed correction is zero -- the
! solve ran for diagnostics only. MAXCH_APP=0 is the truth (no
! population change was committed) and is safe because the current
! LAMBDA update is held while either the GS field or ratio gate is
! unstable. SOL still holds the attempted corrections, so the
! STEQ/DEC_VEC/INC_VEC dumps show what the solve proposed.  Raw
! MAXCH likewise stays at the attempted value for the reports.
!
	IF(GS_POPS_HELD)THEN
	  GS_ATTEMPT_APP=MAXCH_APP
! A DO_BA_RES_ACC rejection inside SOLVEBA_V13 has already restored
! POPS=POPS_BA_SAVE (the pre-step state of the rejected step).
! POPS_GS_SAVE was snapshotted above BEFORE the solve, i.e. on the
! post-step populations the rejection just rolled back -- restoring
! it here would re-commit the wrecked state and freeze it in.  The
! hold must not overwrite a BA rollback (m15_2p82d_knop_LAM_apu
! exploded at it64 exactly this way).
	  IF(.NOT. BA_STEP_REJECTED)POPS=POPS_GS_SAVE
	  MAXCH_APP=0.0_LDP
	  IF(LUM_POPS_HELD)THEN
	    CALL APPEND_ITER_REASON(ITER_REASON,'lum_gate_hold')
	    IF(CONV_DIAG_LVL .GE. 1)WRITE(LUER,'(A,ES12.4)')
	1     ' BA_LUM_GATE: populations held; attempted max change=',
	2     GS_ATTEMPT_APP
	  ELSE
	    CALL APPEND_ITER_REASON(ITER_REASON,'hold_pops')
	    IF(CONV_DIAG_LVL .GE. 1)WRITE(LUER,'(A,ES12.4)')
	1     ' GS_HOLD_POPS: populations held; attempted max change=',
	2     GS_ATTEMPT_APP
	  END IF
	END IF
!
! edited by Devin, 30-Sep-2026: record the temperature-fix state the
! solve actually ran under (MOD_FIXED_T folds in RD_FIX_T and the
! gate-forced fixes) before the mode controller updates FIXED_T for
! the next call.  Feeds MODE_NOW in the GS_STATE/ITER_SUM reports.
!
	FIXT_WAS=MOD_FIXED_T
!
! Temporary diagnostic: localise the origin of a bad temperature
! solution (T lives in POPS(NT,:); SOL(NT,:) is its correction).
! KNOP-path only: the stock iteration sequence is unchanged.
!
	IF(DO_NONMONO_KNOP .AND. CONV_DIAG_LVL .GE. 2)THEN
	  DO I=1,ND
	    IF(.NOT. (ABS(SOL(NT,I)) .LE. HUGE(1.0_LDP)) .OR.
	1       .NOT. (POPS(NT,I) .GT. 0.0_LDP))THEN
	      WRITE(LUER,'(A,I5,A,2ES16.6,A,I4)')' DBG-T: bad T-sol at depth',
	1         I,' SOL(NT),POPS(NT)=',SOL(NT,I),POPS(NT,I),
	2         ' MAIN_COUNTER=',MAIN_COUNTER
	    END IF
	  END DO
	END IF
!
! edited by Devin, 30-Sep-2026: report the physical context (radius,
! velocity, clamped and unclamped velocity gradient) at the depth
! where SOLVEBA proposes its largest correction.  If the worst
! corrections sit where SIGMA_TRUE is most negative (the
! compression/velocity-jump region) that localises the instability
! to the non-monotonic structure rather than a single level.
!
! edited by Devin, 06-Oct-2026: reformatted to the labelled-section
! style (----- NAME it ----- header + aligned rows); fields unchanged.
!
	IF((DO_NONMONO_SOB .OR. DO_NONMONO_KNOP .OR. DO_BA_RES_ACC)
	1    .AND. CONV_DIAG_LVL .GE. 2
	2    .AND. .NOT. BA_STEP_REJECTED)THEN
	  WLOC=MAXLOC(ABS(SOL(1:NT,1:ND)))
	  SIGT=SIGMA(WLOC(2))
	  IF(ALLOCATED(SIGMA_TRUE))SIGT=SIGMA_TRUE(WLOC(2))
	  WRITE(LUER,'(/,A,I5,A)')
	1    ' ---------------------------- BA-CTX it',
	2    MAIN_COUNTER,' -----------------------------------'
	  WRITE(LUER,'(A,L2)')
	1    '    lambda call                   :',LAMBDA_ITERATION
	  WRITE(LUER,'(A,I6,I4)')
	1    '    worst correction @(J,I)       :',WLOC(1),WLOC(2)
	  WRITE(LUER,'(A,4ES12.4)')
	1    '    R/R1, V, sigma, sigma_true    :',
	2    R(WLOC(2))/R(1),V(WLOC(2)),SIGMA(WLOC(2)),SIGT
	END IF
!
! DO_BA_RES_ACC rejection path: SOLVEBA_V13 rolled POPS back to the
! pre-step snapshot because the SE residual grew by more than
! BA_RES_GROW.  The stored BA matrix produced the bad step, so it
! must be rebuilt rather than reused next iteration.
!
! BA_REJ_THIS keeps the rejection visible to the mode controller
! and the diagnostics after the flag is cleared here.
!
	BA_REJ_THIS=BA_STEP_REJECTED
	IF(BA_STEP_REJECTED)THEN
	  COMPUTE_BA=.TRUE.
	  BA_STEP_REJECTED=.FALSE.
	  WRITE(LUER,'(A,I5)')
	1    ' BA_RES_ACC: forcing BA matrix recompute, iter ',MAIN_COUNTER
	  CALL APPEND_ITER_REASON(ITER_REASON,'ba_rejected')
	END IF
!
! edited by Devin, 03-Oct-2026: masked (partial) rollback.  SOLVEBA
! restored the failing depth columns and zeroed their fresh
! corrections; the rest of the step stands.  Force a BA recompute
! -- the stored matrix produced the rejected columns -- but do NOT
! demote: the committed columns made progress and the LAMBDA demote
! is what recreates the alternation the mask exists to break.  The
! partial step also does not count toward the N_FULL_TREL phased
! temperature release below (BA_PART_WAS is consumed there).
!
	BA_PART_WAS=BA_PART_THIS
	IF(BA_PART_THIS)THEN
	  COMPUTE_BA=.TRUE.
	  BA_PART_THIS=.FALSE.
	  CALL APPEND_ITER_REASON(ITER_REASON,'ba_partial_rej')
	END IF
!
! edited by Devin: on [DO_NONMONO_SOB] or [DO_NONMONO_KNOP], a low
! raw demand from a masked or soft-rejected BA step is not evidence
! that the BA direction contracts. With [DO_BA_RES_ACC], keep the
! LAMBDA and JFNK gates closed until a later probe passes residual
! acceptance; both flags off leave the stock path unchanged.
!
	IF((DO_NONMONO_SOB .OR. DO_NONMONO_KNOP) .AND.
	1   BA_RAW_TRIAL_PENDING)THEN
	  BA_RAW_REJECTED=BA_REJ_THIS .OR. BA_PART_WAS .OR.
	1      BA_SOFT_THIS .OR. BA_SNAPSHOT_MISSING
	  BA_RAW_TRIAL_PENDING=.FALSE.
!
! During RD_LAMBDA the only non-LAMBDA call is a [BA_PROBE_LAM] probe,
! so an unjudged-then-clean verdict here is the contraction evidence
! the stable-exit gate exists to find: arm BA_RDLAM_ACCEPT so the
! counting block below pins CNT_STABLE_LAM at N_STABLE_LAM and the
! ordinary DO_LAM_AUTO handoff fires.  A rejected/masked probe leaves
! the gate closed; the next probe is BA_PROBE_LAM calls away.
!
	  IF(RD_LAMBDA .AND. .NOT. BA_RAW_REJECTED)
	1    BA_RDLAM_ACCEPT=.TRUE.
	END IF
	IF((DO_NONMONO_SOB .OR. DO_NONMONO_KNOP) .AND.
	1   BA_SNAPSHOT_MISSING)BA_RAW_REJECTED=.TRUE.
!
! The pending marker must mirror the fresh resacc arm exactly.
! BAACC_PENDING is only armed inside SOLVEBA_V13 on non-LAMBDA calls,
! and any arm from an earlier call was already adjudicated at the top
! of this call's solve, so a set arm here means this call committed an
! unjudged full step.  The previous form also required .NOT. LAM_WAS,
! but LAM_WAS is only assigned at the mode controller below, so it
! tested the PREVIOUS call's mode: every first full call after a
! LAMBDA window (probe, stable exit, cnt_lam escape) skipped the arm
! and its rejection never reached BA_RAW_REJECTED -- dead code in the
! alternating regime the flag exists for.
!
	IF((DO_NONMONO_SOB .OR. DO_NONMONO_KNOP) .AND.
	1   DO_BA_RES_ACC .AND. BAACC_PENDING)
	1   BA_RAW_TRIAL_PENDING=.TRUE.
!
! edited by Devin, 05-Oct-2026: NG window bookkeeping.  The banded
! secant fit in NGACCEL assumes the last 4 SCRTEMP records are
! committed iterates of a (locally) linearly convergent map.  Resacc
! breaks that contract: a rejected step's record is already on disk
! before the residual judgment vetoes it, and the rejecting call
! appends a restored duplicate, so the window holds sign-reversing
! difference pairs.  A population-held call (GS_POPS_HELD) has the
! same defect: POPS is restored to the pre-solve snapshot above, so
! the SCRTEMP write below appends a duplicate of the previous
! record and the window gains a zero-difference pair the fit was
! never meant to see.  NG_CLEAN counts calls since the last
! non-contracting record -- full reject, masked rollback, LAMBDA
! follow-up soft reject, an armed-but-unjudged step (BAACC_PENDING,
! whose commit is only decided next call), a committed JFNK Newton
! iterate (not a Lambda-map iterate; extrapolating across it
! extrapolates the wrong map), or a held-pops record.  The NG fire
! test below defers until 4 clean records exist.  Upstream has no
! resacc/hold/JFNK machinery so all six conditions are always false
! and the upstream path is unchanged.
!
! edited by Devin, 05-Nov-2026: the JFNK term previously tested
! JFNK_THIS, the admission flag -- true on every eligible call while
! DO_JFNK is armed, including defers, cooldown calls, and stagnation
! fallbacks whose SCRTEMP record is an ordinary caller-mode BA
! iterate.  With an admitting gate that reset the window on every
! call and NG could never fire.  JFNK_COMMITTED marks only calls
! that ended in a committed Newton correction, so deferred JFNK
! calls now count toward the window like the BA steps they were.
!
	IF(BA_REJ_THIS .OR. BA_PART_WAS .OR. BA_SOFT_THIS .OR.
	1   BAACC_PENDING .OR. JFNK_COMMITTED .OR. GS_POPS_HELD)THEN
	  NG_CLEAN=0
	ELSE
	  NG_CLEAN=NG_CLEAN+1
	END IF
	BA_SOFT_THIS=.FALSE.
!
! edited by Devin, 30-Sep-2026: accepted-trial refund.  A resacc-netted
! marginal release that the net did NOT roll back is real progress --
! the field reached the marginal tier and the committed step was
! safe, so the GS_FREEZE_LIM/KNOP_STALL clocks should not keep
! counting it as "frozen with no progress".  Refund N_FROZEN to
! MOM_TRIAL_MIN (not zero): the post-release field is re-immature by
! construction, so the next trial still requires the field to
! re-demonstrate the marginal tier before it fires, but the run is
! not re-charged the full first-maturation cost it already paid.
! Rejected trials earn nothing -- N_FROZEN keeps counting, and
! N_TRIAL_ACC bounds consecutive accepted pulses between raw-mature
! sweeps via GS_TRIAL_LIM in the gate above.
!
! The release counters (N_RELEASES_TOT/N_TRIALS_REJ) count events,
! not the GS_TRIAL_LIM pulse budget -- they also bound the cycle
! when GS_TRIAL_LIM is 0 (unbounded pulses) but releases still fire
! through the raw-mature path.  GS_RELEASE_LIM stops the model on
! repeated rejection: the field never matures under a release.
!
	IF(TRIAL_FIRED)THEN
	  IF(.NOT. BA_REJ_THIS)THEN
	    N_RELEASES_TOT=N_RELEASES_TOT+1
	    IF(GS_TRIAL_LIM .GT. 0)THEN
!
! W8 (17-Oct-2026): TRIAL_MIN_USED carries the fired arm's refund
! floor -- MOM_TRIAL_MIN for the MOM arm, FG_TRIAL_MIN for the
! FG-only arm -- so either accepted trial refunds to its own
! re-arm threshold.
!
	      N_FROZEN=TRIAL_MIN_USED
	      N_TRIAL_ACC=N_TRIAL_ACC+1
	      IF(CONV_DIAG_LVL .GE. 1)WRITE(LUER,'(A,I4,A,I4)')
	1      ' BA_RES_ACC: trial release accepted -- n_frozen'//
	2      ' refunded to',N_FROZEN,
	3      '  consecutive accepted trials=',N_TRIAL_ACC
	    END IF
	  ELSE
	    N_TRIALS_REJ=N_TRIALS_REJ+1
	    IF(CONV_DIAG_LVL .GE. 1)WRITE(LUER,'(A)')
	1    ' BA_RES_ACC: trial release rejected --'//
	2    ' freeze clock keeps counting'
	    IF(DO_NONMONO_KNOP .AND. GS_RELEASE_LIM .GT. 0 .AND.
	1      N_TRIALS_REJ .GE. GS_RELEASE_LIM)THEN
	      WRITE(LUER,'(A,I4,A)')'SOLVE_FOR_POPS: ',N_TRIALS_REJ,
	1      ' release attempts rejected; GS field will not mature'
	      WRITE(LUER,'(A)')'  the lagged GS field never matured;'//
	1      ' stopping as non-convergent'
	      STOP
	    END IF
	  END IF
	END IF
!
! edited by Devin, 30-Sep-2026: mode this call actually ran.  The
! KNOP gates above can override the entry mode before SOLVEBA, so
! compose it from the live flags here -- before the mode controller
! below updates them for the next call.
!
	IF(LAMBDA_ITERATION)THEN
	  IF(GS_FREEZE)THEN
	    MODE_NOW='LAMBDA_GS'
	  ELSE IF(FIXT_WAS)THEN
	    MODE_NOW='LAMBDA_TFIX'
	  ELSE
	    MODE_NOW='LAMBDA'
	  END IF
	ELSE IF(FIXT_WAS)THEN
	  MODE_NOW='FULL_TFIX'
	ELSE
	  MODE_NOW='FULL'
	END IF
!
! GS_STATE (moved here 30-Sep-2026): publish the freeze state and
! the post-solve convergence metrics for the job wrapper.
! submit_cmfgen.sh reads GS_STATE after batch.sh and refuses to
! requeue -- or to report a finished run as converged -- when the
! model ends a segment still deep inside the forced-LAMBDA freeze
! (KNOP_STALL_MAX, default 60).  Written every call so the file
! always reflects the last completed iteration; n_frozen=0 is the
! unfrozen state.  n_frozen stays the first key because the wrapper
! greps '^n_frozen=' with -m1.
!
	IF(DO_NONMONO_KNOP)THEN
	  OPEN(NEWUNIT=LU_GS,FILE='GS_STATE',STATUS='REPLACE',
	1    IOSTAT=IOS)
	  IF(IOS .EQ. 0)THEN
	    WRITE(LU_GS,'(A,I0)')'n_frozen=',N_FROZEN
	    WRITE(LU_GS,'(A,I0,A,L2)')'iter=',MAIN_COUNTER,
	1    ' frozen=',GS_FREEZE
	    WRITE(LU_GS,'(A,A)')'mode=',TRIM(MODE_NOW)
	    WRITE(LU_GS,'(A,L2)')'pops_held=',GS_POPS_HELD
	    WRITE(LU_GS,'(A,ES12.4)')'maxch_app=',MAXCH_APP
	    WRITE(LU_GS,'(A,ES12.4)')'maxch_attempt=',GS_ATTEMPT_APP
	    WRITE(LU_GS,'(A,ES12.4)')'se_res=',BA_SE_RES
	    WRITE(LU_GS,'(A,ES12.4)')'se_res2=',BA_SE_RES2
	    WRITE(LU_GS,'(A,ES12.4)')'fg_bad_frac=',
	1    DBLE(GS_NBAD)/DBLE(MAX(GS_NTOT,1_8))
	    WRITE(LU_GS,'(A,ES12.4)')'fg_churn_frac=',
	1    DBLE(GS_NCONV)/DBLE(MAX(GS_NTOT,1_8))
	    WRITE(LU_GS,'(A,ES12.4)')'mom_bad5_frac=',
	1    DBLE(MOM_SW_NBAD5)/DBLE(MAX(MOM_SW_NTOT,1_8))
	    WRITE(LU_GS,'(A,ES12.4)')'mom_bad50_frac=',
	1    DBLE(MOM_SW_NBAD50)/DBLE(MAX(MOM_SW_NTOT,1_8))
	    WRITE(LU_GS,'(A,ES12.4)')'mom_negj_frac=',
	1    DBLE(MOM_SW_NEGJ)/DBLE(MAX(MOM_SW_NTOT,1_8))
	    WRITE(LU_GS,'(A,ES12.4)')'mom_quar_frac=',
	1    DBLE(MOM_SW_NQUAR)/DBLE(MAX(MOM_SW_NTOT,1_8))
	    WRITE(LU_GS,'(A,I0)')'mom_nhuge=',MOM_SW_NHUGE
	    WRITE(LU_GS,'(A,ES12.4)')'mom_maxchg=',MOM_SW_MAXCHG
	    WRITE(LU_GS,'(A,L2)')'mom_refresh_conv=',MOM_RF_CONV
	    WRITE(LU_GS,'(A,I0)')'ba_rej_consec=',BA_NREJ_CONSEC
	    WRITE(LU_GS,'(A,I0)')'ba_sat_cnt=',BA_SAT_CNT
	    WRITE(LU_GS,'(A,ES12.4)')'lum_ratio=',LUM_RATIO_OI
	    WRITE(LU_GS,'(A,I0)')'n_mature=',N_MATURE
	    WRITE(LU_GS,'(A,I0)')'n_rel=',N_REL_SINCE
	    WRITE(LU_GS,'(A,I0)')'n_trial_acc=',N_TRIAL_ACC
	    WRITE(LU_GS,'(A,I0)')'n_trials_rej=',N_TRIALS_REJ
	    WRITE(LU_GS,'(A,I0)')'n_releases=',N_RELEASES_TOT
	    CLOSE(LU_GS)
	  END IF
	END IF
!
! edited by Devin, 29-Sep-2026: under DO_NONMONO_KNOP the control
! metric is the post-limitation committed correction (MAXCH_APP), so
! that a raw 1.0E+7 sentinel from a safely capped population decrease
! does not force an unnecessary Lambda step.  Computed here -- before
! the FIX_BA block -- so the BA-recompute tests below share the
! committed-change metric; a raw sentinel would otherwise pin
! COMPUTE_BA=T and grow MAXCH_SUM forever.  NG scheduling further
! down deliberately keeps the raw MAXCH: a capped-but-huge step is a
! valid reason to defer extrapolation.
!
	IF(DO_NONMONO_KNOP)THEN
	  MAXCH_CTRL=MAXCH_APP
	ELSE
	  MAXCH_CTRL=MAXCH
	END IF
!
! Complicated algorithim to decide when to switch off BA computation.
! We only switch off BA computation in WRBAMAT_RDIN is TRUE.
!
! 1. We can switch off the BA computation, once we are no longer doing
!       LAMBDA iteartions. This is only done for N_ITS_TO_FIX_BA iterations.
!
! 2. We switch of the BA computation semi-permanently once the corrections
!       are less than VAL_FIX_BA.
!
! 3. We recompute BA if the accumulated sum of the maximum corrections
!       exceeds 3*VAL_FIX_BA.
!
! edited by Devin, 30-Sep-2026: a call entered with LAMBDA_ITERATION=F
! and COMPUTE_BA=T had its SE phase store a fresh BAMAT this run --
! only then is switching COMPUTE_BA off safe, since the recall path
! in CMFGEN_SUB trusts BAMATPNT's (unconditionally stamped) success
! flag.  Without this guard the first call to leave LAMBDA mode via
! a mid-call release (MOM_TRIAL) flips COMPUTE_BA off with no BAMAT
! ever written, and the next iteration stops on "cant read BA File".
!
	IF(.NOT. LAMBDA_ENT .AND. COMPUTE_ENT .AND. WRBAMAT)
	1   BA_MATRIX_STORED=.TRUE.
!
	SUM_FIRED=.FALSE.
	IF(N_ITS_TO_FIX_BA .GT. 0 .AND. .NOT. LAMBDA_ITERATION)THEN
	  IF(COMPUTE_BA)THEN
	    CNT_FIX_BA=0
	    MAXCH_SUM=0.0_LDP
	  END IF
	  MAXCH_SUM=MAXCH_SUM+MAXCH_CTRL
	  CNT_FIX_BA=CNT_FIX_BA+1
	  IF(CNT_FIX_BA .GT. N_ITS_TO_FIX_BA .AND. MAXCH_CTRL .GT. VAL_FIX_BA)THEN
	    COMPUTE_BA=L_TRUE
	    WRBAMAT=WRBAMAT_RDIN
	  ELSE
	    IF(WRBAMAT_RDIN .AND. BA_MATRIX_STORED)COMPUTE_BA=L_FALSE
	  END IF
	  IF(CNT_FIX_BA .GT. N_ITS_TO_FIX_BA .AND. MAXCH_SUM .GT. 3.0_LDP*VAL_FIX_BA)THEN
	    MAXCH_SUM=0.0_LDP
	    COMPUTE_BA=L_TRUE
	    WRBAMAT=WRBAMAT_RDIN
	    SUM_FIRED=.TRUE.
	  END IF
!
! This will force BA to be recomuted again, after N_ITS_TO_FIX_BA, because T was variable
! at some depths.
!
	  IF(CON_SCL_T .NE. 0.0_LDP)THEN
	    MAXCH_SUM=1000.0_LDP*VAL_FIX_BA
	  END IF
	ELSE
!
! Switch off BA computation of MAXCH if less than VAL_FIX_BA.
! We only do this provided the last iteration was not a LAMBDA
! iteration, and if T was not partially held fixed at some depths.
!
! We now only right out the BA matrix if we are nearing the correct solution,
! and assuming that it is a full solution matrix (i.e. not from a
! LAMBDA iteration and CON_SCL_T was not set).
!
	  COMPUTE_BA=L_TRUE
	  IF( MAXCH_CTRL .LT. VAL_FIX_BA .AND. WRBAMAT
	1          .AND. .NOT. LAMBDA_ITERATION
	1          .AND. CON_SCL_T .EQ. 0.0_LDP
	1          .AND. BA_MATRIX_STORED)THEN
	    COMPUTE_BA=COMPUTE_BARDIN
	  END IF
	  IF(WRBAMAT_RDIN .AND. MAXCH_CTRL .LT. 2.0_LDP*VAL_FIX_BA
	1         .AND. .NOT. LAMBDA_ITERATION
	1          .AND. CON_SCL_T .EQ. 0.0_LDP)THEN
	    WRBAMAT=.TRUE.
	  END IF
	END IF
!
! edited by Devin, 15-Oct-2026: residual-driven BA recall
! ([DO_BA_RES_REUSE], default F; requires DO_BA_RES_ACC).  The
! legacy scheme reads staleness from the demanded correction
! (MAXCH_CTRL <= VAL_FIX_BA): demand mixes refused (masked or
! trust-capped) demand with real demand, so it cannot tell a stale
! Jacobian from a stiff iterate.  The resacc judgment is the direct
! measure: an armed step whose RES2 ratio BA_JUDGE_RATIO stayed
! <= BA_REUSE_TOL contracted against the stored linearization and
! supports one more recall, even above the VAL_FIX_BA demand
! threshold.  BA_REUSE_NSTALL consecutive stalled judgments force
! one rebuild.  Rejection and the MAXCH_SUM drift trigger (both
! above) still force COMPUTE_BA=T and are honored here via
! BA_REJ_THIS/BA_PART_WAS/SUM_FIRED -- the residual test can only
! extend recall, never veto a rebuild.  Recall still requires a
! stored matrix (BA_MATRIX_STORED, cleared on grid revision),
! WRBAMAT_RDIN, a non-LAMBDA next call, and CON_SCL_T=0.  The
! recall assignment is a hard L_FALSE like the FIX_BA recall above:
! COMPUTE_BARDIN is the COMP_BA VADAT value (T in production), so
! assigning it here would leave COMPUTE_BA=T and silently defeat
! the flag.
!
	IF(DO_BA_RES_REUSE .AND. DO_BA_RES_ACC)THEN
	  IF(BA_JUDGED .AND. .NOT. JFNK_THIS .AND.
	1     .NOT. PREV_COMPUTE_ENT .AND.
	2     BA_JUDGE_RATIO .GT. BA_REUSE_TOL)THEN
	    BA_REUSE_STALL=BA_REUSE_STALL+1
	  ELSE IF(BA_JUDGED .AND. .NOT. JFNK_THIS .AND.
	1     BA_JUDGE_RATIO .LE. BA_REUSE_TOL)THEN
	    BA_REUSE_STALL=0
	  END IF
	  IF(BA_REUSE_STALL .GE. BA_REUSE_NSTALL)THEN
	    COMPUTE_BA=L_TRUE
	    WRBAMAT=WRBAMAT_RDIN
	    BA_REUSE_STALL=0
	    CALL APPEND_ITER_REASON(ITER_REASON,'ba_reuse_stall')
	    IF(CONV_DIAG_LVL .GE. 1)
	1    WRITE(LUER,'(A,I4)')' BA_RES_REUSE: stall limit reached;'//
	2    ' rebuilding BA matrix at it',MAIN_COUNTER
	  ELSE IF(.NOT. LAMBDA_ITERATION .AND. .NOT. SUM_FIRED
	1     .AND. .NOT. JFNK_THIS .AND. BA_MATRIX_STORED
	2     .AND. WRBAMAT_RDIN .AND. WRBAMAT
	3     .AND. .NOT. BA_REJ_THIS .AND. .NOT. BA_PART_WAS
	4     .AND. CON_SCL_T .EQ. 0.0_LDP
	5     .AND. ((BA_JUDGED .AND. BA_JUDGE_RATIO .LE. BA_REUSE_TOL)
	6     .OR. MAXCH_CTRL .LT. VAL_FIX_BA))THEN
	    COMPUTE_BA=L_FALSE
	    IF(CONV_DIAG_LVL .GE. 2 .AND. BA_JUDGED .AND.
	1     MAXCH_CTRL .GE. VAL_FIX_BA)
	2    WRITE(LUER,'(A,ES10.2,A,I5)')
	3    ' BA_RES_REUSE: recall on judged contraction, ratio=',
	4    BA_JUDGE_RATIO,' at it',MAIN_COUNTER
	  END IF
	  PREV_COMPUTE_ENT=COMPUTE_ENT
	END IF
!
!	IF(T_MIN_BA_EXTRAP)COMPUTE_BA=.TRUE.
!
! The STEQ array contains the percentage changes in the populations.
! Shall now determine whether the population changes are too large.
! If so we fix T. Two optical depth ranges are considered. We allow
! for bigger population changes at depth (Tau(es) > 1) so that the
! model converges to the right luminosity rapidly.
!
!**************************************************************************
!             More sophistication may be required.
!**************************************************************************
!
	IF(.NOT. LAMBDA_ITERATION)THEN
	  CON_SCL_T=0.0_LDP
	  CALL ESTAU_V2(TA,R,ED,CLUMP_FAC,TB,ND)
	  DO I=1,ND
	    DO J=1,NT
	      IF( (SOL(J,I) .GT. 0.8_LDP .OR.
	1          SOL(J,I) .LT. -5.0_LDP) .AND. TA(I) .LT. 1.0_LDP)THEN
	            TAU_SCL_T=TA(I)
	            CON_SCL_T=1000.0_LDP
	      END IF
	      IF( (SOL(J,I) .GT. 10.0_LDP .OR.
	1          SOL(J,I) .LT. -20.0_LDP) .AND. TA(I) .GE. 1.0_LDP)THEN
	            TAU_SCL_T=TA(I)
	            CON_SCL_T=1000.0_LDP
	      END IF
	    END DO
	  END DO
	END IF
!	IF(CON_SCL_T .NE. 0.0D0)THEN
!	   WRBAMAT=L_FALSE
!	   COMPUTE_BA=L_TRUE
!	END IF
!
        CALL WR2D_V2(SOL,NT,ND,'STEQ SOLUTION ARRAY','#',L_TRUE,LU_SE)
!
! GS_FREEZE iterations ran a forced LAMBDA update, not a BA step; the
! smooth/undo/average machinery assumes a real BA SOL and is skipped
! while the lagged coupling field is immature.
!
	IF(AUTO_SMOOTH_POPS .AND. .NOT. GS_FREEZE)THEN
	  CALL SMOOTH_POPS_AS_WE_ITERATE(POPS,SOL,ND,NT)
	END IF
!
! NB: The call to SUM_STEQ_SOL corrupts SOL. I is used for output --
!      the unit is closed on exit.
!
	I=7; CALL SUM_STEQ_SOL(SOL,NT,ND,I)
!
! Determine whether convergence is sufficient to consider using
! NG acceleration. The first NG acceleration is done 4 iterations after
! the iteration on which MAXCH < VAL_DO_NG. IST_PER_NG must be greater
! than, or equal to, 4. LAST_NG now always indicates the last iteration
! on which an NG acceleration was performed. NEXT_NG is used to indicate
! the next iteration on which an NG acceleration is to occur.
!
! NEXT_NG=1000 indicates a NEW model
! NEXT_NG=1500 indicates that MAXCH is again above VAL_DO_NG
! NEXT_NG=2000 indicates a CONTINUING model.
!
! A NG acceleration can be forced after 1 iteration if LAST_NG is set to
! some vale .LE. MAIN_COUNTER-ITE_PER_NG in the POINT1 file, and provide the
! change on that iteration is less than VAL_DO_NG.
!
! edited by Devin, 05-Oct-2026: a rewrite of NEXT_NG here supersedes
! any NG_CLEAN deferral target armed last call, so the deferred-due
! flag below is dropped whenever the scheduler changes the value.
!
	I=NEXT_NG
	IF(MAIN_COUNTER .LT. IT_TO_BEG_NG-3)THEN
	  NEXT_NG=1000
	ELSE IF(MAIN_COUNTER .GE. IT_TO_BEG_NG-3 .AND. NEXT_NG .EQ. 1000)THEN
	  IF( MAXCH .LE. VAL_DO_NG )NEXT_NG=MAX(LAST_NG+ITS_PER_NG,MAIN_COUNTER+4)
	ELSE IF(MAXCH .LT. VAL_DO_NG .AND. NEXT_NG .EQ. 1500)THEN
	  NEXT_NG=MAX(LAST_NG+ITS_PER_NG,MAIN_COUNTER+4)
	ELSE IF(MAXCH .LT. VAL_DO_NG .AND. NEXT_NG .EQ. 2000)THEN
	  NEXT_NG=MAX(LAST_NG+ITS_PER_NG,MAIN_COUNTER+1)
	  IF(LAST_NG .EQ. -1000)NEXT_NG=MAIN_COUNTER+4
	ELSE IF( MAXCH .GE. VAL_DO_NG )THEN
	  NEXT_NG=1500
	END IF
	IF(NEXT_NG .NE. I)NG_DEFERRED=.FALSE.
!
! We switch between LINEARIZATION and LAMBDA iterations until the
! maximum percentage change is less  than VAL_DO_LAM. RD_CNT_LAM iterations
! are performed per full linearization. RD_LAMBDA overrides this section if
! TRUE. We only do a LAMBDA iteration with T fixed. NB --- FIX_IMPURITY
! is a soft option --- is removed when convergence nearly obtained.
!
	IF(LAMBDA_ITERATION)LAST_LAMBDA=MAIN_COUNTER
	LAM_WAS=LAMBDA_ITERATION
!
! Accepted-full counter for the phased T release: increments only
! for a non-LAMBDA step not rejected, masked, or committed as an
! unconverged JFNK iterate. It resets on (re)entry to LAMBDA below.
!
! A partially rejected step committed progress at the kept columns
! but is not a clean acceptance: it does not advance the phased
! temperature-release counter.
!
	IF(.NOT. LAM_WAS .AND. .NOT. BA_REJ_THIS .AND.
	1   .NOT. BA_PART_WAS .AND. .NOT. JFNK_UNCONVERGED_COMMIT)
	2   N_ACC_FULL=N_ACC_FULL+1
!
! Saturated-demand demotion: a non-LAMBDA step whose raw correction
! demand stays at the sentinel means the applied update rides the
! limitation cap -- marching rather than contracting.  The count is
! keyed to the entry mode so the step this call committed decides.
! Residual acceptance cannot anticipate the basin cliff (SE_RES2
! stayed flat through LAM_dbg it36-47 while demand stayed pinned),
! so the demand trend is the only pre-cliff signal available.
!
	IF((DO_NONMONO_SOB .OR. DO_NONMONO_KNOP) .AND. .NOT. LAM_WAS
	1   .AND. MAXCH .GE. BA_SAT_RAW)THEN
	  BA_SAT_CNT=BA_SAT_CNT+1
	ELSE
	  BA_SAT_CNT=0
	END IF
!
! edited by Devin, 03-Oct-2026: record the raw demand of every
! non-LAMBDA step for the BA_ENTER_RAW merit-exit gate below.  The
! raw value -- not MAXCH_APP -- is what diagnoses a non-contractive
! BA direction: a saturated decrease stays huge even when the
! committed update was capped to a few percent.  Recording is
! unconditional telemetry; the gate itself is inert unless the
! VADAT [BA_ENTER_RAW] key is set and a nonmono flag
! (DO_NONMONO_SOB or DO_NONMONO_KNOP) is on.
!
! edited by Devin, 05-Oct-2026: a successful JFNK call returns the
! applied correction in MAXCH, not the raw demand -- writing it here
! would open the demand gate and JFNK admission on an applied-step
! value.  Keep the last BA-measured demand until a real SOLVEBA
! re-measures (conservative direction: gates stay held).
!
	IF(.NOT. LAM_WAS .AND. .NOT. JFNK_THIS)BA_LAST_RAW=MAXCH
! FORCE_LAMBDA (KNOP warm-up above) pins LAMBDA_ITERATION for the
! whole N_LAM_RWARM window; without this guard the generic MAXCH test
! below could flip the model into a BA step against a stale I_NEXT.
!
! edited by Devin, 29-Sep-2026: the KNOP mode decision uses the
! post-limitation applied correction (MAXCH_APP) so that a raw
! 1.0E+7 sentinel from a safely capped population decrease does not
! force an unnecessary Lambda step.  MAXCH_CTRL is set above the
! FIX_BA block; raw MAXCH is retained for NG scheduling and the
! printed diagnostic.
!
! edited by Devin, 30-Sep-2026: hysteresis replaces the single
! VAL_DO_LAM switch.  With defaults (LAM_ENTER=LAM_EXIT=LAM_VAL,
! LAM_N_STABLE=1) the behavior is identical to the legacy block: one
! below-threshold iteration exits LAMBDA and one at/above-threshold
! full step re-enters it.  When the thresholds differ, a full step
! in the dead band VAL_EXIT <= MAXCH_CTRL < VAL_ENTER stays full,
! and a LAMBDA sequence must string N_STABLE_LAM quiet iterations
! together before releasing -- a single lucky iteration can no
! longer flip the mode.  CNT_STABLE_LAM is reset on entry and by a
! forced-LAMBDA window so freeze iterations never count toward the
! stability requirement.
!
! edited by Devin: DO_LAM_AUTO's initial non-monotonic handoff uses
! the same applied-change stability count as later LAMBDA exits.  A
! forced/frozen or luminosity-drifting call cannot establish readiness.
	IF(FORCE_LAMBDA)THEN
	  CNT_STABLE_LAM=0
	  LAM_STABLE_PROVEN=.FALSE.
	END IF
	IF(RD_LAMBDA .AND. DO_LAMBDA_AUTO .AND. LAMBDA_ITERATION .AND.
	1   (DO_NONMONO_SOB .OR. DO_NONMONO_KNOP))THEN
	  IF(FORCE_LAMBDA .OR. LUM_BLOCKED .OR.
	1      MAXCH_CTRL .GE. VAL_EXIT_LAM)THEN
	    CNT_STABLE_LAM=0
	    LAM_STABLE_PROVEN=.FALSE.
	  ELSE
	    CNT_STABLE_LAM=CNT_STABLE_LAM+1
	    IF(CNT_STABLE_LAM .GE. N_STABLE_LAM)LAM_STABLE_PROVEN=.TRUE.
	  END IF
!
! edited by Devin, 11-Oct-2026: [BA_PROBE_LAM]-cadence bounded probe
! inside the RD_LAMBDA phase.  The stable-exit gate above has no
! escape when the LAMBDA attractor floors above VAL_EXIT_LAM -- the
! m15 SOB run parked at ~45% demand with stable=0/3 while every
! probe/escape path sits behind .NOT. RD_LAMBDA.  Reuse CNT_LAM
! (unused before handoff) as the cadence counter and fire one full
! call when it is armed: residual acceptance must judge the trial
! (DO_BA_RES_ACC) and a factor bound must hold it (BA_DAMP_CAP), the
! demand must sit in the band where a bounded step is informative
! (VAL_EXIT_LAM <= MAXCH_CTRL < VAL_ENTER_LAM -- below the exit a
! normal stable count may still hand off; above the enter the model
! is still churning), and no force/hold may be active.  CMFGEN_SUB
! sees BA_RDLAM_PROBE and leaves the call full; resacc adjudicates it
! next call, so a bad direction costs one capped step, not a phase.
!
	  CNT_LAM=CNT_LAM+1
	  IF(DO_BA_RES_ACC .AND. BA_PROBE_LAM .GT. 0 .AND.
	1     BA_DAMP_CAP .GT. 0.0_LDP .AND.
	2     CNT_LAM .GE. BA_PROBE_LAM .AND.
	3     MAXCH_CTRL .GE. VAL_EXIT_LAM .AND.
	4     MAXCH_CTRL .LT. VAL_ENTER_LAM .AND.
	5     .NOT. FORCE_LAMBDA .AND. .NOT. LUM_BLOCKED)THEN
	    BA_RDLAM_PROBE=.TRUE.
	    CNT_LAM=0
	    CALL APPEND_ITER_REASON(ITER_REASON,'rdlam_probe')
	  END IF
!
! A resacc-accepted probe is sustained evidence that a full step
! contracts here, so keep the exit count pinned; the demand-level
! reset above would otherwise disarm it on the same call.  cmfgen_sub
! still requires MAXCH < 50 for the switch, so the arm only matters
! while the LAMBDA demand stays sane.  Dead once RD_LAMBDA clears.
!
	  IF(BA_RDLAM_ACCEPT)THEN
	    CNT_STABLE_LAM=N_STABLE_LAM
	    LAM_STABLE_PROVEN=.TRUE.
	  END IF
	END IF
	IF(.NOT. RD_LAMBDA .AND. .NOT. FORCE_LAMBDA)THEN
	  IF(LAMBDA_ITERATION)THEN
	     IF(MAXCH_CTRL .LT. VAL_EXIT_LAM)THEN
	        CNT_STABLE_LAM=CNT_STABLE_LAM+1
	     ELSE
	        CNT_STABLE_LAM=0
	        LAM_STABLE_PROVEN=.FALSE.
	     END IF
	     IF(CNT_STABLE_LAM .GE. N_STABLE_LAM)THEN
	       LAM_STABLE_PROVEN=.TRUE.
	       IF(LUM_BLOCKED)THEN
!
! BA_LUM_GATE: the pops are quiet enough to leave LAMBDA but the
! inner-boundary luminosity is still drifting -- a BA step here is
! armed against a transient inner state (LAM_alt it13-35: lum
! excursion -28.6 -> -0.5 over ~20 its while maxraw rode the cap).
! Hold the exit; CNT_STABLE_LAM stays armed so the first
! lum-settled call exits without re-earning the stability count.
!
	        CALL APPEND_ITER_REASON(ITER_REASON,'lum_drift_hold')
	       ELSE IF((DO_NONMONO_SOB .OR. DO_NONMONO_KNOP) .AND.
	1          (BA_LAST_RAW .GE. BA_ENTER_RAW .OR.
	2           (BA_ENTER_RAW .LT. 1.0E30_LDP .AND.
	3            (BA_RAW_REJECTED .OR. BA_RAW_TRIAL_PENDING))))THEN
	        IF(CNT_LAM .GE. MAX(BA_PROBE_LAM, RD_CNT_LAM))THEN
!
! edited by Devin, 03-Oct-2026: held-gate probe (VADAT
! [BA_PROBE_LAM], default 0 = off).  BA_LAST_RAW only refreshes on a
! non-LAMBDA call, so a held gate can never learn that the BA
! direction has started to contract -- without this the only
! re-open channel is the NUM_LAM escape (SOB v2: the 1e7 probe at
! it515 held every merit exit, next probe ~it715).  After
! BA_PROBE_LAM held LAMBDA calls, release one ordinary full step
! with the same actions as lam_stable_exit -- resacc and BA_TRUST
! bound it identically.  If it demotes, CNT_LAM resets and the next
! probe is BA_PROBE_LAM calls away, so probing cannot spam; if it
! is accepted the gate re-opens on the refreshed BA_LAST_RAW.
!
	          FIX_IMPURITY=RD_FIX_IMP
	          LAMBDA_ITERATION=.FALSE.
	          FIXED_T=RD_FIX_T
	          CNT_LAM=0
	          CNT_STABLE_LAM=0
	          CALL APPEND_ITER_REASON(ITER_REASON,'ba_demand_probe')
	          IF(FRESH_BA_EXIT)THEN
	            COMPUTE_BA=.TRUE.
	            WRBAMAT=.TRUE.
	            CNT_FIX_BA=0
	            MAXCH_SUM=0.0_LDP
	          END IF
	          IF(N_FULL_TREL .GT. 0)THEN
	            T_HELD_PHASE=.TRUE.
	            N_ACC_FULL=0
	            FIXED_T=.TRUE.
	          END IF
	        ELSE
!
! edited by Devin, 03-Oct-2026: BA_ENTER_RAW demand gate (nonmono
! paths; extended SOB -> KNOP 3-Oct-2026 after the opt Knop run
! showed the identical cycle).  The LAMBDA iterate is quiet but the
! last full step's raw demand is still saturated, so quiet is the
! LAMBDA attractor, not evidence the Newton map contracts --
! releasing now just re-arms the exit/demote limit cycle.
! CNT_STABLE_LAM is reset so CNT_LAM resumes counting and the
! NUM_LAM escape still fires: each escape commits one damped step
! that re-measures BA_LAST_RAW, which is the only channel that can
! re-open the gate.
!
! edited by Devin, 04-Oct-2026: count the held call in CNT_LAM.
! CNT_LAM previously only advanced on calls whose applied change
! exceeded VAL_EXIT_LAM; under release damping every held call is
! "stable", CNT_STABLE_LAM re-armed to 1 before the else-branch
! could run, and CNT_LAM stayed 0 -- so both the NUM_LAM escape
! and the BA_PROBE_LAM probe were dead code and an armed
! BA_ENTER_RAW locked LAMBDA mode permanently.  Now each held
! call counts, and BA_PROBE_LAM (or NUM_LAM) releases a probe
! full step that re-measures BA_LAST_RAW.
!
	          CALL APPEND_ITER_REASON(ITER_REASON,'ba_demand_hold')
	          CNT_LAM=CNT_LAM+1
	          CNT_STABLE_LAM=0
	        END IF
	       ELSE
	        FIX_IMPURITY=RD_FIX_IMP
	        LAMBDA_ITERATION=.FALSE.
	        FIXED_T=RD_FIX_T
	        CNT_LAM=0
	        CNT_STABLE_LAM=0
	        CALL APPEND_ITER_REASON(ITER_REASON,'lam_stable_exit')
	        IF(FRESH_BA_EXIT)THEN
	          COMPUTE_BA=.TRUE.
	          WRBAMAT=.TRUE.
	          CNT_FIX_BA=0
	          MAXCH_SUM=0.0_LDP
	        END IF
	        IF(N_FULL_TREL .GT. 0)THEN
	          T_HELD_PHASE=.TRUE.
	          N_ACC_FULL=0
	          FIXED_T=.TRUE.
	        END IF
	       END IF
	     ELSE
	        CNT_LAM=CNT_LAM+1
	        IF(CNT_LAM .GE. RD_CNT_LAM)THEN
	          LAMBDA_ITERATION=.FALSE.
	          FIX_IMPURITY=RD_FIX_IMP
	          FIXED_T=RD_FIX_T
	          CALL APPEND_ITER_REASON(ITER_REASON,'cnt_lam_escape')
	          IF(FRESH_BA_EXIT)THEN
	            COMPUTE_BA=.TRUE.
	            WRBAMAT=.TRUE.
	            CNT_FIX_BA=0
	            MAXCH_SUM=0.0_LDP
	          END IF
	          IF(N_FULL_TREL .GT. 0)THEN
	            T_HELD_PHASE=.TRUE.
	            N_ACC_FULL=0
	            FIXED_T=.TRUE.
	          END IF
	        END IF
	        IF(MAXCH_CTRL .GT. 1.0E+05_LDP)THEN
	          LAMBDA_ITERATION=.TRUE.
	          COMPUTE_BA=.TRUE.
	          FIX_IMPURITY=.FALSE.
	          FIXED_T=.TRUE.
	        END IF
	     END IF
	  ELSE
!
! BA_REJ_THIS demotes explicitly -- the resacc rollback means the
! stored BA produced the bad step, so recover in LAMBDA rather than
! relying on the rejection sentinel clearing VAL_ENTER_LAM (which
! fails if the user sets LAM_ENTER above the ~1e2 sentinel).
!
	     IF(BA_REJ_THIS)THEN
	        LAMBDA_ITERATION=.TRUE.
	        COMPUTE_BA=.TRUE.
	        CNT_LAM=0
	        CNT_STABLE_LAM=0
	        LAM_STABLE_PROVEN=.FALSE.
	        N_ACC_FULL=0
	        T_HELD_PHASE=.FALSE.
	        FIXED_T=.TRUE.
	        FIX_IMPURITY=.FALSE.
!
! Saturated BA demand demotes before the stay-full test: under KNOP
! the applied correction rides the cap (MAXCH_CTRL < VAL_ENTER_LAM)
! even while the march walks off the basin, so a converging-looking
! committed change cannot be allowed to veto the demand check.
! Extended to DO_NONMONO_SOB 03-Oct-2026: on SOB the count provides
! the bounded leash for [BA_DAMP_CAP] damped-full runs.
!
	     ELSE IF((DO_NONMONO_SOB .OR. DO_NONMONO_KNOP) .AND.
	1        BA_SAT_LIM .GT. 0 .AND.
	2        BA_SAT_CNT .GE. BA_SAT_LIM)THEN
	        LAMBDA_ITERATION=.TRUE.
	        COMPUTE_BA=.TRUE.
	        CNT_LAM=0
	        CNT_STABLE_LAM=0
	        LAM_STABLE_PROVEN=.FALSE.
	        N_ACC_FULL=0
	        T_HELD_PHASE=.FALSE.
	        FIXED_T=.TRUE.
	        FIX_IMPURITY=.FALSE.
	        BA_SAT_CNT=0
	        BA_TRUST_CUR=MAX(1.0_LDP,0.5_LDP*BA_TRUST_CUR)
	        CALL APPEND_ITER_REASON(ITER_REASON,'ba_sat_demote')
!
! edited by Devin, 03-Oct-2026: [BA_DAMP_CAP] damped-full mode.  A
! step whose raw demand is saturated stays in full mode instead of
! demoting: the next call re-measures the BA direction under the
! BA_DAMP_CAP factor bound, so the model progresses in small
! increments rather than paying a demote-plus-LAMBDA window per
! probe.  The test is on the raw MAXCH (not MAXCH_CTRL) so a small
! damped applied correction cannot slip into the stay-full branch
! and release T.  The run length is bounded by ba_sat_demote above
! (BA_SAT_LIM, default 8 when the cap is armed) and each step still
! passes residual acceptance.  Inert unless [BA_DAMP_CAP] > 0.
!
	     ELSE IF(BA_DAMP_CAP .GT. 0.0_LDP .AND.
	1        MAXCH .GE. BA_SAT_RAW)THEN
	        FIX_IMPURITY=.FALSE.
	        CNT_LAM=0
	        CNT_STABLE_LAM=0
	        LAM_STABLE_PROVEN=.FALSE.
	        FIXED_T=.TRUE.
	        CALL APPEND_ITER_REASON(ITER_REASON,'ba_damped_full')
	     ELSE IF(MAXCH_CTRL .LT. VAL_ENTER_LAM)THEN
	        FIX_IMPURITY=.FALSE.
	        CNT_LAM=0
	        CNT_STABLE_LAM=0
	        IF(T_HELD_PHASE)THEN
	          IF(N_ACC_FULL .GE. N_FULL_TREL .AND.
	1             MAXCH_CTRL .LT. T_REL_CHG .AND.
	2             .NOT. JFNK_UNCONVERGED_COMMIT)THEN
	            FIXED_T=RD_FIX_T
	            T_HELD_PHASE=.FALSE.
	            CALL APPEND_ITER_REASON(ITER_REASON,
	1             'temperature_released')
	          ELSE
	            FIXED_T=.TRUE.
	          END IF
	        ELSE
	          FIXED_T=RD_FIX_T
	        END IF
	     ELSE
	        LAMBDA_ITERATION=.TRUE.
	        COMPUTE_BA=.TRUE.
	        CNT_LAM=0
	        CNT_STABLE_LAM=0
	        LAM_STABLE_PROVEN=.FALSE.
	        N_ACC_FULL=0
	        T_HELD_PHASE=.FALSE.
	        FIXED_T=.TRUE.
	        FIX_IMPURITY=.FALSE.
	        CALL APPEND_ITER_REASON(ITER_REASON,'maxch_above_enter')
	     END IF
	  END IF
	END IF
!
! edited by Devin, 29-Sep-2026: trace the LAMBDA/full-iteration mode
! controller each iteration.  The Dessart SOB run ping-pongs BA<->LAMBDA
! because the rescaled decrease term (100*x/(1-x)) lifts MAXCH above
! VAL_DO_LAM on every BA step; this block records the decision inputs
! so the cycle can be diagnosed from OUTGEN.
!
! edited by Devin, 06-Oct-2026: reformatted to the labelled-section
! style (----- NAME it ----- header + aligned rows); fields unchanged.
!
	IF(CONV_DIAG_LVL .GE. 2 .AND. (DO_NONMONO_SOB .OR.
	1    DO_NONMONO_KNOP .OR. DO_BA_RES_ACC))THEN
	  WRITE(LUER,'(/,A,I5,A)')
	1    ' ---------------------------- MODCTL it',
	2    MAIN_COUNTER,' -----------------------------------'
	  WRITE(LUER,'(A,L2)')
	1    '    was lambda call               :',LAM_WAS
	  WRITE(LUER,'(A,ES12.4,2X,ES12.4)')
	1    '    MAXCH demand / applied        :',MAXCH,MAXCH_APP
	  WRITE(LUER,'(A,ES12.4,2X,ES12.4)')
	1    '    lambda enter / exit           :',VAL_ENTER_LAM,VAL_EXIT_LAM
	  WRITE(LUER,'(A,I3,A,I3,4X,I3,A,I3)')
	1    '    lambda / stable counters      :',CNT_LAM,'/',RD_CNT_LAM,
	2    CNT_STABLE_LAM,'/',N_STABLE_LAM
	  WRITE(LUER,'(A,ES12.4)')
	1    '    last raw MAXCH                :',BA_LAST_RAW
	  WRITE(LUER,'(A,ES12.4)')
	1    '    exempted demand (BA_DEM floors):',BA_DEM_EXM
	  WRITE(LUER,'(A,ES12.4,2X,A,I3,A,I3,2X,A,L2)')
	1    '    BA_TR_ABS (|dn|/gas)          :',BA_TR_ABS,
	2    'sweep',BA_TR_NSWP,'/',BA_TR_MAXSWP,'phase',BA_TRACE_PHASE
	  WRITE(LUER,'(A,L2,3X,L2,3X,L2)')
	1    '    next lam / comp_ba / fix_t    :',LAMBDA_ITERATION,
	2    COMPUTE_BA,FIXED_T
	  WRITE(LUER,'(A,L2,3X,L2,3X,L2)')
	1    '    rd_lam / force / gsfrz        :',RD_LAMBDA,FORCE_LAMBDA,
	2    GS_FREEZE
	END IF
!
! Automatically adjust R grid, so that grid is uniformally spaced on the
! FLUX optical depth scale. Used for SN models with very sharp ioinization
! fronts. By doing it before the output to SCRTEMP, we ensure that
! a continuuing model starts with the revised R grid.
!
! Write pointer file and output data necessary to begin a new
! iteration.
!
	CALL SCR_RITE_V2(R,V,SIGMA,POPS,IREC,MAIN_COUNTER,RITE_N_TIMES,
	1                LAST_NG,WRITE_RVSIG,NT,ND,LUSCR,NEWMOD)
!
! Program to create a NEW_R_GRID which is equally spaced in LOG(Tau) where
! TAU is based on the FLUX mean opacity.
!
! In VADAT REVISE_R_GRID should be set to TRUE.
!
! Initialize R_OLD unconditionally: DO_NG_BAND_ACCEL_V4 reads it even
! when REVISE_R_GRID is F.
!
	R_OLD(1:ND)=R(1:ND)
	IF(REVISE_R_GRID .AND. .NOT. GS_FREEZE)THEN
	  CALL ESOPAC(ESEC,ED,ND)               !Electron scattering emission factor.
          CALL ADJUST_R_GRID_V4(POPS,ESEC,MAIN_COUNTER,R_GRID_REVISED,ND,NT)
	  IF(R_GRID_REVISED)THEN
! ADJUST_R_GRID_V4 interpolates R, V and SIGMA onto the new grid, so the
! unclamped gradient used by DO_NONMONO_SOB/DO_NONMONO_KNOP must be
! refreshed as well.
	     IF(DO_NONMONO_SOB .OR. DO_NONMONO_KNOP)
	1        CALL RECOMPUTE_SIGMA_TRUE(SIGMA_TRUE,SIGMA,R,V,ND)
!
! edited by Devin, 13-Oct-2026: a stored BAMAT was assembled on the old
! grid -- its depth-column indices no longer map to the same radii, so
! a FIX_BA/N_FIX_BA recall would apply corrections to the wrong cells.
! Clear the stored flag; the next full iteration rebuilds on the new
! grid.  Not flag-gated: a stale matrix is wrong on every path that
! revises the grid (upstream gap -- the post-regrid call is forced
! LAMBDA in cmfgen_sub.f, which does not store, so the flag only ever
! protects WRBAMAT recalls on later calls).
!
	     BA_MATRIX_STORED=.FALSE.
	     MAIN_COUNTER=MAIN_COUNTER+1
	     CALL SCR_RITE_V2(R,V,SIGMA,POPS,IREC,MAIN_COUNTER,RITE_N_TIMES,
	1                LAST_NG,WRITE_RVSIG,NT,ND,LUSCR,NEWMOD)
	  END IF
	ELSE
	  R_GRID_REVISED=.FALSE.
	END IF
! 
!
! Perform the acceleration. T1 and T2 return the percentage changes
! in the populations. I is used as an INTEGER error flag. We never do
! an NG acceleration on the last iteration.
!
! We perform a LAMBDA iteration after the NG acceleration in the
! case the where the NG acceleration has caused population changes
! > 20%. In this case, the BA matrix will also be recaluated on
! the next full iteration.
!
!	WRITE(6,*)MAIN_COUNTER,NEXT_AV,LAST_LAMBDA+NUM_OSC_AV+2
!	WRITE(6,*)NEXT_NG-4,LAST_AV+ITS_PER_AV
!	WRITE(6,*)AVERAGE_DO,LST_ITERATION
!
	AVERAGE_DONE=.FALSE.
	IF(MAIN_COUNTER .NE. LAST_LAMBDA .AND. UNDO_LAST_IT
	1          .AND. .NOT. GS_FREEZE)THEN
	  I=1; J=5
	  CALL  UNDO_IT(SOL,POPS,NT,ND,MAIN_COUNTER,I,J,AVERAGE_DONE)
	  IF(AVERAGE_DONE)THEN
	    WRITE(LUER,'(A,I4,A,I4)')'Last corrections undone for depths',I,' to',J
	  ELSE
	    WRITE(LUER,*)'Error undoing last corrections.'
	  END IF
	END IF
!
	IF(AVERAGE_DO .AND. MAIN_COUNTER .GE. NEXT_AV
	1             .AND. MAIN_COUNTER .GT. LAST_NG+NUM_OSC_AV
	1             .AND. MAIN_COUNTER .GT. LAST_LAMBDA+NUM_OSC_AV+2
	1             .AND. (MAIN_COUNTER .LT. NEXT_NG-4 .OR. .NOT. NG_DO)
	1             .AND. MAIN_COUNTER .GE. LAST_AV+ITS_PER_AV
	1             .AND. .NOT. GS_FREEZE
	1             .AND. .NOT. LST_ITERATION)THEN
	  CALL AVE_FLIPS(SOL,POPS,NT,ND,MAIN_COUNTER,NUM_OSC_AV,L_TRUE,AVERAGE_DONE)
	  IF(AVERAGE_DONE)THEN
	    WRITE(LUER,*)'Averaging of oscilating populations',
	1                   ' performed on iteration',MAIN_COUNTER
	    LAST_AV=MAIN_COUNTER
	    NEXT_AV=MAIN_COUNTER+ITS_PER_AV
	  ELSE
	    NEXT_AV=MAIN_COUNTER+8
	    WRITE(LUER,*)'Error performing AV acceleration.'
	    WRITE(LUER,*)'Will try again in 8 iterations.'
	  END IF
	END IF
!
! Save the corrected poplatons to SCREMP.
!
	IF(AVERAGE_DONE)THEN
	  MAIN_COUNTER=MAIN_COUNTER+1
	  POPS=SOL
	  CALL SCR_RITE_V2(R,V,SIGMA,POPS,IREC,MAIN_COUNTER,
	1             RITE_N_TIMES,LAST_NG,WRITE_RVSIG,NT,ND,LUSCR,NEWMOD)
	END IF
!
! FORCE_LAMBDA guard added: NEXT_NG is scheduled independently of
! LAST_LAMBDA, so without it a queued NG extrapolation could fire on
! a KNOP warm-up iteration (FORCE_LAMBDA=T, GS_FREEZE=F if the field
! was restored mature) and commit a population kick across an I_NEXT
! reinit.  The AV/UNDO paths are already protected implicitly via
! MAIN_COUNTER .NE. LAST_LAMBDA.
!
! edited by Devin, 05-Oct-2026: the deferral below arms
! NEXT_NG=MAIN_COUNTER+1, but an auxiliary SCRTEMP write (R-grid
! revision or AV averaging, which bump MAIN_COUNTER) can step over
! the +1 target before the due call, leaving NEXT_NG permanently
! below MAIN_COUNTER -- the upstream .EQ. test then never matches
! and NG silently never fires.  NG_DEFERRED marks targets this
! branch armed itself; only for those is the due test widened to
! .LE. .  Sentinel values (1000/1500/2000) and every other
! scheduling site keep the upstream .EQ. semantics, so the stock
! path is unchanged.
!
	IF(NG_DO .AND. (.NOT. LST_ITERATION) .AND.
	1       .NOT. GS_FREEZE .AND. .NOT. FORCE_LAMBDA .AND.
	1       (NEXT_NG .EQ. MAIN_COUNTER .OR.
	2        (NG_DEFERRED .AND. NEXT_NG .LE. MAIN_COUNTER)) )THEN
	  IF(NG_CLEAN .LT. 4)THEN
! The 4-record SCRTEMP window still holds a rejected, rolled-back,
! or unjudged (BAACC_PENDING) record; extrapolating it can only
! produce a vetoed acceleration.  Deferring by one call (not 8)
! lets NG fire as soon as the window ages clean -- with resacc on,
! reject-dense phases would otherwise burn every attempt.  Upstream
! path is unchanged: NG_CLEAN never stays <4 without resacc/JFNK.
! NG_DEFERRED records that this branch armed the +1 so a stepped-
! over target is still honoured (see the due test above).
	    NEXT_NG=MAIN_COUNTER+1
	    NG_DEFERRED=.TRUE.
	  ELSE
	  CALL DO_NG_BAND_ACCEL_V4(POPS,R,V,SIGMA,R_OLD,
	1         NT,ND,NG_BAND_WIDTH,NG_DONE,T1,T2,
	1         DO_NG_VALIDITY_CHECK,NG_ABSTAIN,LUSCR,LUER)
	  IF(NG_DONE)THEN
! DO_NG_BAND_ACCEL_V4 extrapolates R, V and SIGMA together with the
! populations, so the unclamped gradient used by the non-monotonic
! paths must be refreshed from the accelerated V.
	    IF(DO_NONMONO_SOB .OR. DO_NONMONO_KNOP)
	1      CALL RECOMPUTE_SIGMA_TRUE(SIGMA_TRUE,SIGMA,R,V,ND)
	    MAIN_COUNTER=MAIN_COUNTER+1
	    LAST_NG=MAIN_COUNTER
	    NEXT_NG=MAIN_COUNTER+ITS_PER_NG
	    NG_DEFERRED=.FALSE.
	    NG_DEM_PRE=MAXCH
	    CALL SCR_RITE_V2(R,V,SIGMA,POPS,IREC,MAIN_COUNTER,
	1             RITE_N_TIMES,LAST_NG,WRITE_RVSIG,NT,ND,LUSCR,NEWMOD)
	    IF(T1 .GT. 50.0_LDP .AND. T2 .GT. 50.0_LDP)THEN
	      LAMBDA_ITERATION=.TRUE.
	      CNT_LAM=0
	      FIXED_T=.TRUE.
	      FIX_IMPURITY=.FALSE.
	      COMPUTE_BA=.TRUE.
! edited by Devin, 17-Oct-2026: ng_kick is a demote branch like the
! resacc/saturation demotes below -- it must clear the same
! hysteresis state.  A surviving LAM_STABLE_PROVEN would short-
! circuit the stable-LAMBDA relief (and the stable-exit gate) on the
! next FULL call with a proof earned before the extrapolation.
	      CNT_STABLE_LAM=0
	      LAM_STABLE_PROVEN=.FALSE.
	      N_ACC_FULL=0
	      T_HELD_PHASE=.FALSE.
	      CALL APPEND_ITER_REASON(ITER_REASON,'ng_kick')
	    END IF
	  ELSE
	    NEXT_NG=MAIN_COUNTER+8
	    NG_DEFERRED=.FALSE.
	    WRITE(LUER,*)'Error performing NG acceleration.'
	    WRITE(LUER,*)'Will try again in 8 iterations.'
	    WRITE(LUER,*)'Error flag=',I
	  END IF
	  END IF
	END IF
	FLUSH(LUER)
!
! Here we check whether we keep BA fixed. Previously we always computed the BA
! matrix after a LAMBDA iteration. This will allow us to keep it fixed for
! at least 1 iteration. Since the computation of the BA matrix takes considerable time,
! this may reduce computational effort.
!
	IF(R_GRID_REVISED)THEN
	   LAMBDA_ITERATION=.TRUE.
	   COMPUTE_BA=.TRUE.
	   COUNT_FULL_LAM_BA=1000
	   CALL APPEND_ITER_REASON(ITER_REASON,'rgrid_revised')
	ELSE IF(LAST_LAMBDA .EQ. MAIN_COUNTER .AND. .NOT. LAMBDA_ITERATION)THEN
	   IF(COUNT_FULL_LAM_BA .NE. 0)THEN
	      COUNT_FULL_LAM_BA=0
	   ELSE
	     OPEN(UNIT=7, FILE='BAMATPNT',STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	     IF(IOS .EQ. 0)READ(7,*,IOSTAT=IOS)BA_SUCCESSFULLY_OUTPUT
	     INQUIRE(UNIT=7,OPENED=FILE_OPEN)
	     IF(FILE_OPEN)CLOSE(UNIT=7)
	     IF(IOS .EQ. 0 .AND. BA_SUCCESSFULLY_OUTPUT)THEN
!	       COMPUTE_BA=.FALSE.
	       COUNT_FULL_LAM_BA=COUNT_FULL_LAM_BA+1
	     END IF
	  END IF
	 END IF
!
! edited by Devin, 29-Sep-2026: final next-iteration mode after the
! post-controller overrides (NG kick, R-grid revision, averaging).
! MODCTL above records the VAL_DO_LAM decision; this block records what
! the iteration actually becomes after NG/AV/REVISE_R may have flipped
! it.
!
! edited by Devin, 06-Oct-2026: reformatted to the labelled-section
! style (----- NAME it ----- header + aligned rows); fields unchanged.
!
	IF(CONV_DIAG_LVL .GE. 2 .AND. (DO_NONMONO_SOB .OR.
	1    DO_NONMONO_KNOP .OR. DO_BA_RES_ACC))THEN
	  WRITE(LUER,'(/,A,I5,A)')
	1    ' --------------------------- MODE-OUT it',
	2    MAIN_COUNTER,' ----------------------------------'
	  WRITE(LUER,'(A,L2,3X,L2,3X,L2)')
	1    '    next lam / comp_ba / fix_t    :',LAMBDA_ITERATION,
	2    COMPUTE_BA,FIXED_T
	  WRITE(LUER,'(A,2I6)')
	1    '    next NG / AV                  :',NEXT_NG,NEXT_AV
	  WRITE(LUER,'(A,2I6)')
	1    '    last lambda / NG              :',LAST_LAMBDA,LAST_NG
	  WRITE(LUER,'(A,L2,3X,L2,3X,L2)')
	1    '    NG done / avg done / rgrid    :',NG_DONE,AVERAGE_DONE,
	2    R_GRID_REVISED
	END IF
!
! edited by Devin, 30-Sep-2026: one-line convergence record per
! great iteration plus a matching ITER_HISTORY file and a plateau
! detector.  Emitted on every call that runs a dev path, including
! CONV_DIAG_LVL=0 -- the point of the level control is to strip the
! per-call trace noise, not the summary.
!
	IF(DO_NONMONO_SOB .OR. DO_NONMONO_KNOP .OR. DO_BA_RES_ACC)THEN
!
! Mode selected for the next call, after the controller and the
! NG/R-grid/AV overrides above.
!
	  IF(LAMBDA_ITERATION)THEN
	    IF(FIXED_T)THEN
	      MODE_NEXT='LAMBDA_TFIX'
	    ELSE
	      MODE_NEXT='LAMBDA'
	    END IF
	  ELSE IF(FIXED_T)THEN
	    MODE_NEXT='FULL_TFIX'
	  ELSE
	    MODE_NEXT='FULL'
	  END IF
!
! ITER_MODE transition record: every forcing/decision site above
! appended its reason; a bare mode flip with no recorded reason gets
! 'untracked' so gaps in the coverage are visible rather than silent.
!
	  IF(CONV_DIAG_LVL .GE. 1)THEN
	    IF(ITER_REASON .EQ. ' ' .AND.
	1       MODE_NOW .NE. MODE_NEXT)ITER_REASON='untracked'
	    IF(ITER_REASON .NE. ' ')WRITE(LUER,'(A,I5,6A)')
	1     ' ITER_MODE: it=',MAIN_COUNTER,' ',
	2     TRIM(MODE_NOW),' -> ',TRIM(MODE_NEXT),
	3     ' reason=',TRIM(ITER_REASON)
	  END IF
	  IF(BA_REJ_THIS)THEN
	    BA_CH='REJECT'
	  ELSE IF(BA_PART_WAS)THEN
	    BA_CH='MASKED'
	  ELSE IF(LAM_WAS)THEN
	    BA_CH='LAMBDA'
	  ELSE IF(BAACC_PENDING)THEN
	    BA_CH='ARMED'
	  ELSE
	    BA_CH='HELD'
	  END IF
	  GS_DRIFT_FRAC=0.0_LDP
	  IF(MOM_SW_NTOT .GT. 0_8)THEN
	    GS_DRIFT_FRAC=DBLE(MOM_SW_NBAD50)/DBLE(MOM_SW_NTOT)
	  ELSE IF(GS_NTOT .GT. 0_8)THEN
	    GS_DRIFT_FRAC=DBLE(GS_NCONV)/DBLE(GS_NTOT)
	  END IF
!
! Field legend, once per process: the table rows below use short
! names so the per-iteration block stays compact.
!
	  WRITE(LUER,'(/,A)')' -------------------------------- ITER SUMMARY ---------------------------------'
	  IF(.NOT. ITER_LEG_DONE)THEN
	    ITER_LEG_DONE=.TRUE.
	    WRITE(LUER,'(A)')'  ITER_SUM legend:'
	    WRITE(LUER,'(A)')'   demand/applied = largest population'//
	1    ' correction wanted / committed this call'
	    WRITE(LUER,'(A)')'   se / se2       = max SE residual;'//
	1    ' se2 excludes the worst depth columns'
	    WRITE(LUER,'(A)')'   trust / nrej   = BA correction cap'//
	1    ' (halves on reject) / consecutive rejections'
	    WRITE(LUER,'(A)')'   lum            = inner/outer luminosity'//
	1    ' ratio and relative drift (diagnostic only for shocks)'
	    WRITE(LUER,'(A)')'   fg / mom       = lagged-field sweep'//
	1    ' metrics: d = max rel. change; bad = cells >50%;'//
	2    ' churn = cells >5%; n5,n50 = mom-cell counts,'//
	3    ' n50 printed as count/threshold; neg ='//
	4    ' significant negative-J repairs; quar/pquar ='//
	5    ' quarantined cells (pquar persistent); cap ='//
	6    ' dropped non-finite writes'
	    WRITE(LUER,'(A)')'   status         = [OK] inside the'//
	1    ' maturity gates, [WARN] marginal tier,'//
	2    ' [BAD] unsafe/failed, [HELD] pops frozen,'//
	3    ' [INFO] diagnostic only, [REJ] BA rejected'//
	4    ' (ba=MASKED = partial rollback)'
	  END IF
!
	  WRITE(LUER,'(A,I5,3(A,A),2(A,L2))')
	1    ' ITER_SUM: it=',MAIN_COUNTER,'  mode=',TRIM(MODE_NOW),
	2    ' -> ',TRIM(MODE_NEXT),'  ba=',TRIM(BA_CH),
	3    '  frozen=',GS_FREEZE,'  held=',GS_POPS_HELD
!
! Status tags mirror the gate thresholds above: [BAD] is the
! catastrophic tier and [WARN] the marginal tier.  Luminosity carries
! [INFO] because shocks make inner/outer agreement non-convergent.
!
	  POP_ST='[OK]   '
	  IF(GS_POPS_HELD)POP_ST='[HELD] '
	  SE_ST='[OK]   '
	  IF(.NOT. (BA_SE_RES .LT. 1.0E+30_LDP))SE_ST='[BAD]  '
	  BA_ST='[OK]   '
	  IF(BA_NREJ_CONSEC .GT. 0)BA_ST='[REJ]  '
	  LUM_ST='[INFO] '
!
	  WRITE(LUER,'(4X,A14,2(2X,ES10.3),2X,A,1X,A)')
	1    'pops',MAXCH,MAXCH_APP,POP_ST,
	2    'demand / applied correction'
	  IF(NG_DEM_PRE .GE. 0.0_LDP)THEN
	    WRITE(LUER,'(4X,A14,2(2X,ES10.3),2X,A)')
	1    'ng rebound',NG_DEM_PRE,MAXCH,
	2    'demand before last NG / this call'
	    NG_DEM_PRE=-1.0_LDP
	  END IF
	  WRITE(LUER,'(4X,A14,2(2X,ES10.3),2X,A,1X,A)')
	1    'se residual',BA_SE_RES,BA_SE_RES2,SE_ST,
	2    'all cells / excl. worst cols'
!
! Per-depth residual profile (edited by Devin, 15-Oct-2026): under
! CONV_DIAG_LVL list the five depths carrying the largest column-max
! SE residual and each one's share of the column sum.  BA_COL_RES
! is exported by SOLVEBA_V13 on every call, so this also covers
! LAMBDA calls and resacc-off runs where the residual used to be
! invisible (m15 dirchk_lc printed se residual 0 while the JFNK
! |F|_max was ~1E12, concentrated at depths 66-70 and He II).
!
	  IF(CONV_DIAG_LVL .GE. 1 .AND. ALLOCATED(BA_COL_RES)
	1     .AND. SIZE(BA_COL_RES) .EQ. ND)THEN
	    BLOCK
	      REAL(KIND=LDP) RES_TOP(ND),RES_SUM
	      INTEGER ITOP,ID_TOP
	      RES_TOP(1:ND)=BA_COL_RES(1:ND)
	      RES_SUM=SUM(RES_TOP)
	      DO ITOP=1,MIN(5,ND)
	        ID_TOP=MAXLOC(RES_TOP,1)
	        IF(RES_TOP(ID_TOP) .LE. 0.0_LDP)EXIT
	        WRITE(LUER,'(4X,A14,1X,I5,2X,ES10.3,2X,F5.1,A)')
	1        'res depth',ID_TOP,RES_TOP(ID_TOP),
	2        100.0_LDP*RES_TOP(ID_TOP)/MAX(RES_SUM,1.0E-30_LDP),
	3        '% of column sum'
	        RES_TOP(ID_TOP)=0.0_LDP
	      END DO
	    END BLOCK
	  END IF
	  WRITE(LUER,'(4X,A14,2X,ES10.3,2X,I4,7X,A,1X,A)')
	1    'ba',BA_TRUST_CUR,BA_NREJ_CONSEC,BA_ST,
	2    'trust / consec rejects'
	  WRITE(LUER,'(4X,A14,2(2X,ES10.3),2X,A,1X,A)')
	1    'lum diag',LUM_RATIO_OI,LUM_DCHG_REL,LUM_ST,
	2    'inner/outer ratio / relative drift'
	  IF(DO_NONMONO_KNOP)THEN
	    FG_ST='[OK]   '
	    IF(GS_CAPPED .GT. 0)THEN
	      FG_ST='[BAD]  '
	    ELSE IF(GS_FG_IMMATURE)THEN
	      FG_ST='[WARN] '
	    END IF
	    WRITE(LUER,'(4X,A14,2X,ES10.3,1X,A,I10,1X,A,I10,1X,A,I6,'//
	1    '2X,A,1X,A)')
	2    'fg',GS_MAXCHG,'bad=',GS_NBAD,'churn=',GS_NCONV,
	3    'cap=',GS_CAPPED,FG_ST,'lagged field drift'
	    IF(MOM_RAN)THEN
	      MOM_ST='[OK]   '
	      IF(MOM_UNSAFE_B)THEN
	        MOM_ST='[BAD]  '
	      ELSE IF(MOM_FREEZE_B)THEN
	        MOM_ST='[WARN] '
	      END IF
	      WRITE(LUER,'(4X,A14,2X,ES10.3,1X,A,I10,1X,A,I10,1X,A,I10,'//
	1      '1X,A,I10,1X,A,I10,1X,A,I10,1X,A,I10,2X,A)')
	2      'mom',MOM_SW_MAXCHG,'n5=',MOM_SW_NBAD5,
	3      'n50=',MOM_SW_NBAD50,'/',MOM_SW_NTOT,
	4      'neg=',MOM_SW_NEGJ,'quar=',MOM_SW_NQUAR,
	5      'pquar=',MOM_NPQUAR,'cap=',MOM_SW_CAPPED,MOM_ST
	    END IF
	  END IF
!
! JFNK gate-state row: demand vs both tier gates plus the
! consecutive-hit counters, so an OUTGEN grep of ITER_SUM blocks
! alone reconstructs why Newton did or did not fire.  [ARMED] =
! admitted this call (coupled flag disambiguates the tier),
! [DEFER] = armed but gated out, [OFF] = still inside JFNK_WARMUP.
!
	  IF(USE_JFNK_SOLVER)THEN
	    JFNK_ST='[DEFER]'
	    IF(JFNK_THIS)JFNK_ST='[ARMED]'
	    IF(ITS_SINCE_START .LE. JFNK_WARMUP)JFNK_ST='[OFF]  '
	    WRITE(LUER,'(4X,A14,2(2X,ES10.3),2X,A,1X,A,I3,A,I3,1X,A)')
	1    'jfnk',MAXCH,MAX(JFNK_GATE_IN,0.0_LDP),JFNK_ST,
	2    'hits frz=',JFNK_FRZ_HITS,'/',JFNK_GATE_NREQ,
	3    'demand vs frz gate'
	    WRITE(LUER,'(4X,A14,2(2X,ES10.3),7X,A,I3,A,I3,1X,A,1X,L2)')
	1    'jfnk_cpl',MAXCH,MAX(JFNK_COR_GATE,0.0_LDP),
	2    'hits cpl=',JFNK_CPL_HITS,'/',JFNK_GATE_NREQ,
	3    'vs cpl gate; coupled=',JFNK_COUPLED_THIS
	  END IF
	  WRITE(LUER,'(A)')' -------------------------------------------------------------------------------'
!
! Append one row per iteration; header only when the file is new.
! The file survives requeues in the run directory, which is what the
! plateau plot needs.
!
	  INQUIRE(FILE='ITER_HISTORY',EXIST=HIST_EXIST)
	  OPEN(NEWUNIT=LU_GS,FILE='ITER_HISTORY',STATUS='UNKNOWN',
	1    POSITION='APPEND',IOSTAT=IOS)
	  IF(IOS .EQ. 0)THEN
	    IF(.NOT. HIST_EXIST)WRITE(LU_GS,'(A)')'# iter mode next ba'//
	1    ' froz held maxraw maxapp maxatt se se2 trust nrej lum'//
	2    ' fg_d fg_badf fg_churnf fg_cap'//
	3    ' mom_d mom_n5f mom_n50f mom_negjf mom_quarf'//
	4    ' mom_cap mom_pquar'
	    WRITE(LU_GS,'(I6,1X,A11,1X,A11,1X,A6,2(1X,L1),'//
	1    '6(1X,ES11.4),1X,I4,12(1X,ES11.4))')
	2    MAIN_COUNTER,MODE_NOW,MODE_NEXT,BA_CH,GS_FREEZE,
	3    GS_POPS_HELD,
	4    MAXCH,MAXCH_APP,GS_ATTEMPT_APP,BA_SE_RES,BA_SE_RES2,
	5    BA_TRUST_CUR,BA_NREJ_CONSEC,LUM_RATIO_OI,
	6    GS_MAXCHG,DBLE(GS_NBAD)/DBLE(MAX(GS_NTOT,1_8)),
	7    DBLE(GS_NCONV)/DBLE(MAX(GS_NTOT,1_8)),DBLE(GS_CAPPED),
	8    MOM_SW_MAXCHG,
	9    DBLE(MOM_SW_NBAD5)/DBLE(MAX(MOM_SW_NTOT,1_8)),
	1    DBLE(MOM_SW_NBAD50)/DBLE(MAX(MOM_SW_NTOT,1_8)),
	2    DBLE(MOM_SW_NEGJ)/DBLE(MAX(MOM_SW_NTOT,1_8)),
	3    DBLE(MOM_SW_NQUAR)/DBLE(MAX(MOM_SW_NTOT,1_8)),
	4    DBLE(MOM_SW_CAPPED),DBLE(MOM_NPQUAR)
	    CLOSE(LU_GS)
	  END IF
!
! Plateau detector: flag a stall only when the committed correction,
! the second-channel SE residual and the GS drift fraction all sit
! inside a 1.5x window over the last PLAT_WIN calls.  A held-pops
! sweep (MAXCH_APP=0) and a path without resacc (BA_SE_RES2=0) drop
! the corresponding test instead of manufacturing a flat zero
! series.
!
	  HIST_POS=MOD(HIST_POS,PLAT_WIN)+1
	  HIST_APP(HIST_POS)=MAXCH_APP
	  HIST_RES(HIST_POS)=BA_SE_RES2
	  HIST_GS(HIST_POS)=GS_DRIFT_FRAC
	  N_HIST=MIN(N_HIST+1,PLAT_WIN)
	  PLATEAU_NOW=.FALSE.
	  IF(N_HIST .GE. PLAT_WIN .AND. MINVAL(HIST_APP) .GT.
	1     0.0_LDP)THEN
	    PLATEAU_NOW=MAXVAL(HIST_APP)/MINVAL(HIST_APP) .LT.
	1       1.5_LDP
	    IF(PLATEAU_NOW .AND. MAXVAL(HIST_RES) .GT. 0.0_LDP)
	1      PLATEAU_NOW=MAXVAL(HIST_RES)/
	2      MAX(MINVAL(HIST_RES),1.0E-300_LDP) .LT. 1.5_LDP
	    IF(PLATEAU_NOW .AND. MAXVAL(HIST_GS) .GT. 0.0_LDP)
	1      PLATEAU_NOW=MAXVAL(HIST_GS)/
	2      MAX(MINVAL(HIST_GS),1.0E-30_LDP) .LT. 1.5_LDP
	  END IF
	  IF(PLATEAU_NOW)THEN
	    PLAT_CNT=PLAT_CNT+1
	    IF(CONV_DIAG_LVL .GE. 1 .AND.
	1       (PLAT_CNT .EQ. 1 .OR. MOD(PLAT_CNT,PLAT_WIN) .EQ. 0))THEN
	      WRITE(LUER,'(A,I5,A,2ES11.3,A,2ES11.3,A,2ES11.3)')
	1      ' CONV_PLATEAU: it=',MAIN_COUNTER,
	2      ' maxapp range=',MINVAL(HIST_APP),MAXVAL(HIST_APP),
	3      ' se_res2 range=',MINVAL(HIST_RES),MAXVAL(HIST_RES),
	4      ' gs_drift range=',MINVAL(HIST_GS),MAXVAL(HIST_GS)
	    END IF
	  ELSE
	    PLAT_CNT=0
	  END IF
	END IF
!
	RETURN
	END
!
! Accumulate ITER_MODE transition reasons; several forcing events can
! fire on one call (e.g. gs_freeze_forced + hold_pops).
!
	SUBROUTINE APPEND_ITER_REASON(RCD,ADD)
	CHARACTER(LEN=*) RCD,ADD
	IF(RCD .EQ. ' ')THEN
	  RCD=ADD
	ELSE
	  RCD=TRIM(RCD)//'+'//ADD
	END IF
	RETURN
	END
