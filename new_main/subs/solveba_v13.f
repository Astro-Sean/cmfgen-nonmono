	SUBROUTINE SOLVEBA_V13(STEQ,POPS,POP_ATOM,DIAG_INDX,NT,NION,NUM_BNDS,ND,
	1              MAXCH,MAXCH_APP,METH_SOL,SUCCESS,SCALE_OPT,
	1              LAM_SCALE_OPT,CHANGE_LIM,MAX_dT_COR,T_MIN,
	1              BA_COMPUTED,WR_BA_INV,WR_PRT_INV,LAMBDA_IT,
	1              MAIN_COUNTER,SET_POPS_D2_EQ_D1)
	USE SET_KIND_MODULE
	USE CONTROL_VARIABLE_MOD, ONLY: LTE_MODEL,
	1     DO_BA_RES_ACC,BA_RES_GROW,BA_RES_HARD,BA_TRUST,
	2     BA_TRUST_CUR,BA_SE_RES,
	3     BAACC_PENDING,BA_STEP_REJECTED,BAACC_RES_PREV,
	4     BAACC_RES_MIN,BA_NREJ,BA_NREJ_CONSEC,POPS_BA_SAVE,
	5     DO_NONMONO_SOB,DO_NONMONO_KNOP,BA_SE_RES2,
	6     BAACC_RES2_PREV,BAACC_RES2_MIN,BA_NREJ_LIM,
	7     CONV_DIAG_LVL,BA_SNAPSHOT_MISSING,BA_RES_SAT,
	8     BAACC_COL_PREV,BAACC_COL_MIN,BA_RBK_MASK,
	9     BA_PART_THIS,BA_REJ_FRAC,BA_NREJ_SOFT,BA_SOFT_THIS,
	1     BAACC_RDLAM_PROBE_ARM,BAACC_RDLAM_PROBE,
	1     USE_ELEC_HEAT_BAL,GS_POPS_HELD,DO_BA_QUAR,BA_QUAR_MAX,BA_QUAR_MASK,
	1     BA_QUAR_COUNT,BA_DEM_PFLOOR,
	2     BA_JUDGED,BA_JUDGE_RATIO,BA_MAXCH_FAIL,BA_COL_RES,
	3     BA_DEM_EXM,BA_TRACE_PHASE,BA_TR_ABS,POPS_GS_SAVE,
	4     BA_HOLD_RESTORE_REVIVAL,BA_EX_RELFRAC,BA_EX_MAXCH,
	5     BA_DEM_EXR
	USE MOD_CMFGEN, ONLY: NUM_SPECIES,SPECIES_BEG_ID,
	1     SPECIES_END_ID,EQ_SPECIES,ATM,MOD_FIXED_NE,
	2     MOD_FIXED_T
	USE STEQ_DATA_MOD, ONLY: SE,STEQ_ED,STEQ_T,STEQ_T_EHB
	IMPLICIT NONE
!
! Altered 16-Oct-2026 : Largest exempted (floor-hidden) demand is
!                          exported as BA_DEM_EXM and named in OUTGEN;
!                          the MAXCH metric itself is unchanged
!                          (edited by Devin).
! Altered 15-Oct-2026 : SE residual measured on every call (not only
!                          under DO_BA_RES_ACC); BA_COL_RES exported
!                          for the ITER_SUM depth profile (edited by Devin).
! Altered 15-Oct-2026 : BA_MAXCH_FAIL same-call rejection of a committed
!                          catastrophic full step (edited by Devin).
! Altered 07-Oct-2026 : Per-equation quarantine (DO_BA_QUAR) under the masked residual
!                          acceptance (edited by Devin).
! Altered 07-Oct-2026 : Held-call and residual-progress accounting for BA_NREJ_CONSEC (edited by Devin).
! Altered 05-Oct-2026 : Subheader for the per-iteration population-update report (edited by Devin).
! Altered 03-Oct-2026 : Per-depth masked acceptance (BA_REJ_FRAC) and honest accounting of armed
!                          steps on LAMBDA follow-ups (edited by Devin).
! Altered 01-Oct-2026 : Residual acceptance metric reads assembled SE residual, not the STEQ dummy (edited by Devin).
! Altered 01-Oct-2026 : Saturated-residual exemption for the acceptance test (edited by Devin).
! Altered 29-Sep-2026 : Hard stop at BA_NREJ_LIM consecutive rejections (edited by Devin).
! Altered 08-Dec-2013 : No longer scale the temperature when D2 EQ D1 -- Not T(2)=T(1)
! Altered 17-Feb-2016 : Added LTE_MODEL option.
! Created 14-Feb-2014 : Changed to V13 -- added MAX_dT_COR to call.
! Created 26-Mar-2013 : Changed to V11 and added POP_ATOM to call.
! Altered 12-MAr-2013 : We can now set the populations at depth 2 equal to those at depth 1 for
!                          a FULL iteration (i.e., non=LAMBDA). This is designed to overcome a convergence
!                          issue with SN models. This should only be done when depth 2 is very close to
!                          depth 1.
! Altered 08-Nov-2011 : Updated maximum correction output with info on whether BA was computed.
! Altered 15-Nov-2010 : FIDDLE_POP_CORRECTIONS added to assist convergence. This allows
!                         -ver relaxation to be used at some depths.
! Altered 13-Oct-2026 : BA_DEM_PFLOOR / BA_DEM_RELFRAC now also gate the
!                          MAX_INC_VEC / MAX_DEC_VEC accumulation, so the
!                          DO_LEVEL_CHK NV-th-largest MAXCH respects the
!                          dead-level floors (previously bypassed them).
! Altered 17-Oct-2026 : BA_DEM_RELFRAC removed entirely -- relative
!                          exemption hid a converged-exit on an
!                          unconverged model (SOB l14 it242).
! Altered 23-Apr-2010 : DO_LEVEL_CHECK, MAX_INC_VEC, & MAX_DEC_VEC variables installed.
!                         Return 10th maximum value, rather than biggest correction.
!                         Done to help with a few non-converging levels.
! Altered 23-Nov-2007 : Changed to V9
!                       LAM_SCALE_OPT inserted in call.
! Altered 23-May-2006 : Extra column output to CORRECTION_SUM.
! Altered 19-May-2005 : T_MIN inserted into call.
! Altered 11-Feb-2004 : Bug fixed with calculation of SCALE in 'GLOBAL' scale option.
! Altered 30-Mar-2003 : Output to CORRECTION_SUM file included.
! Altered 28-Jan-2002 : Changed to V& as BA_COMPUTED and WR_BA_INV now passed in the CALL.
!                       CMF_BLKBAND_V3 now used.
!
	INTEGER NT
	INTEGER NION
	INTEGER NUM_BNDS
	INTEGER ND
	INTEGER DIAG_INDX
	INTEGER MAIN_COUNTER
!
	REAL(KIND=LDP) STEQ(NT,ND)
	REAL(KIND=LDP) POPS(NT,ND)
	REAL(KIND=LDP) POP_ATOM(ND)
	REAL(KIND=LDP) MAXCH
	REAL(KIND=LDP) MAXCH_APP
	REAL(KIND=LDP) CHANGE_LIM
	REAL(KIND=LDP) MAX_dT_COR
	REAL(KIND=LDP) T_MIN
	CHARACTER*(*) METH_SOL
	CHARACTER*(*) SCALE_OPT
	CHARACTER*(*) LAM_SCALE_OPT
	LOGICAL SUCCESS
	LOGICAL BA_COMPUTED
	LOGICAL WR_BA_INV
	LOGICAL WR_PRT_INV
	LOGICAL LAMBDA_IT
	LOGICAL SET_POPS_D2_EQ_D1
!
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
	LOGICAL FG_RB_OK
	EXTERNAL RESTORE_FG_INEXT_LATEST
	LOGICAL MOM_RB_OK
	EXTERNAL RESTORE_MOM_JNEXT_LATEST
!
	INTEGER, PARAMETER :: IZERO=0
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: ITWO=2
!
!
! Local variables.
!
	INTEGER, PARAMETER :: NV=10
	REAL(KIND=LDP) MAX_INC_VEC(NV)
	REAL(KIND=LDP) MAX_DEC_VEC(NV)
!
	REAL(KIND=LDP) SCALE,MINSCALE,T1,T2,T3
	REAL(KIND=LDP) INCREASE,DECREASE
	REAL(KIND=LDP) STEQ_MASKED_MAX
	REAL(KIND=LDP) STEQ_EXM
	INTEGER IEXM,JEXM
	INTEGER ITRC,JTRC
	REAL(KIND=LDP) EX_DEC,EX_INC,EXR_VAL,T2EX
	INTEGER IEXR,JEXR
	REAL(KIND=LDP) BIG_LIM,LIT_LIM
! Allocatable, not automatic: NT*ND is too large to trust to the stack
! on nodes with a limited stack size.
	REAL(KIND=LDP), ALLOCATABLE :: POPS_OLD(:,:)
	INTEGER I,J,K,IINC,IDEC,IOS,LU_SUM
	INTEGER JINC,JDEC
	INTEGER RES_LOC(2)
	INTEGER RES2_COL
	INTEGER ISPEC,ID,JJ
	REAL(KIND=LDP) COL_RES(ND)
	REAL(KIND=LDP) COL_RES_S(ND)
	INTEGER NFAIL
	INTEGER JK,NQUAR_NEW
	LOGICAL LCOL_OK,LMASK_OK,LQUAR_OK,L_PROG,L_ANY_MASK,
	1    LADV_OK
	REAL(KIND=LDP) STEQ_QMAX
	LOGICAL REJ_BULK,REJ_MAX,REJ_STATE,REJ_PART,WLD_REJ
	LOGICAL, SAVE :: BAACC_R2GATE_PRINTED=.FALSE.
	INTEGER, SAVE :: PV_JI=0,PV_II=0,PV_JD=0,PV_ID=0
	LOGICAL OSC_FL
	INTEGER COUNT(7)
	LOGICAL LOC_WR_BA_INV
	LOGICAL DO_LEVEL_CHK
!
	LUER=ERROR_LU()
!
! Per-equation quarantine ([DO_BA_QUAR], default F): the mask is
! allocated (NT,ND) here on first use and persists across calls --
! a parked equation stays parked for the rest of the run.
! LQUAR_OK is the per-call guard -- the mask must exist and match
! the current NT/ND.  Set unconditionally so the BA13-TRACE report
! below can read it even when DO_BA_RES_ACC is off.
!
	LQUAR_OK=.FALSE.
	IF(DO_BA_RES_ACC .AND. DO_BA_QUAR)THEN
	  IF(.NOT. ALLOCATED(BA_QUAR_MASK))THEN
	    ALLOCATE(BA_QUAR_MASK(NT,ND),STAT=IOS)
	    IF(IOS .EQ. 0)BA_QUAR_MASK(1:NT,1:ND)=.FALSE.
	  ELSE IF(SIZE(BA_QUAR_MASK,1) .NE. NT .OR.
	1     SIZE(BA_QUAR_MASK,2) .NE. ND)THEN
	    DEALLOCATE(BA_QUAR_MASK)
	    ALLOCATE(BA_QUAR_MASK(NT,ND),STAT=IOS)
	    IF(IOS .EQ. 0)BA_QUAR_MASK(1:NT,1:ND)=.FALSE.
	  ELSE
!
! No alloc ran: a stale IOS from an earlier call must not decide
! whether the mask is usable -- same hazard as POPS_BA_SAVE below.
!
	    IOS=0
	  END IF
	  LQUAR_OK=(IOS .EQ. 0)
	END IF
	NQUAR_NEW=0
	STEQ_QMAX=0.0_LDP
!
! Residual-accepted BA update (VADAT DO_BA_RES_ACC, default F).
! STEQ enters as the statistical-equilibrium residual; its max
! modulus is the acceptance metric.  If the previous BA step left
! the residual grown by more than BA_RES_GROW, that step was based
! on a stale or invalid Jacobian: restore the pre-step POPS
! snapshot, discard this iteration's correction entirely, and flag
! the rejection so SOLVE_FOR_POPS forces a fresh matrix next
! iteration.  A hard bound BA_RES_HARD*BAACC_RES_MIN (smallest
! residual seen so far) closes the slow-growth loophole in which a
! run of accepted sub-GROW steps diverges geometrically.
! Rejection halves the trust radius BA_TRUST_CUR; acceptance
! recovers it toward BA_TRUST. Ordinary LAMBDA follow-ups are not
! rolled back, to avoid repeating the same BA/LAMBDA cycle. A rejected
! BA_PROBE_LAM trial is the exception: its follow-up restores the saved
! populations and radiation state, then skips the stale Lambda solve.
!
! edited by Devin, 03-Oct-2026: per-depth masked acceptance
! ([BA_REJ_FRAC], default 0.25; <=0 restores the legacy all-or-
! nothing test).  The observed failures are local: a few stiff
! cells at compression depths drive the residual rebound while
! the rest of the step contracts (m15: H2:lev1@d70, CV:1s_2p_1Po
! @d13), and damping cannot help -- SOB v2 it509 capped the applied
! step at +50% and the residual still rebounded 120x; the rebound
! is set by the direction, not the step size.  A step that fails
! at a minority of columns is now partially committed: failing
! columns are restored from POPS_BA_SAVE and their fresh STEQ is
! zeroed below, passing columns keep their update, the call is not
! wasted, and SOLVE_FOR_POPS is not forced into the LAMBDA demote
! that recreates the FULL<->LAMBDA alternation.  The RES2 rank
! channel stays the bulk veto and rolls the whole state back.
!
! A LAMBDA follow-up still cannot roll the armed step back (the
! it33 cycle note below stands), but it now MEASURES it: a degraded
! residual halves the trust radius and increments BA_NREJ_SOFT
! instead of silently clearing the flag as an acceptance.
!
! edited by Devin, 01-Oct-2026: when the armed step's pre-step
! residual is already at/above BA_RES_SAT (default 1E30), the max
! channel is measuring a degenerate equation, not step quality --
! a bounded correction cannot move it, and any overflow above the
! floor trips the growth test spuriously, burning BA_NREJ_LIM on
! meaningless rejections (SOB run: res pinned at 1.5E262 from it 1,
! killed at 20 consec while res2 stayed ~90).  In that regime only
! the RES2 rank channel judges the step; BA_RES_SAT <= 0 restores
! the legacy all-channel test.
!
! edited by Devin, 15-Oct-2026: the residual MEASUREMENT below runs
! on every call, not only under DO_BA_RES_ACC.  The ITER_SUM 'se
! residual' channel reads BA_SE_RES/BA_SE_RES2; with resacc off it
! printed 0.000E+00 on every iteration (m15 dirchk_lc its 170+),
! leaving the run blind while the JFNK residual was ~1E12.  The
! scan is O(species*ND*N_SE) -- cheap next to the BA build.  Only
! the acceptance bookkeeping and adjudication stay gated; the
! BAACC_*_MIN/BA_RBK_MASK/RES_LOC housekeeping here is write-only
! when resacc is off and changes no solver behaviour.
!
!
! edited by Devin, 01-Oct-2026: STEQ is output-only scratch at this
! point -- the caller (SOLVE_FOR_POPS) passes uninitialized automatic
! storage (local SOL array) and the residual is only assembled into
! it inside CMF_BLKBAND_V3 below (GENERATE_FULL_MATRIX_V3 copies
! SE%STEQ, STEQ_ED and STEQ_T into each depth column).  The previous
! version measured MAXVAL(ABS(STEQ)) here, i.e. stack garbage: on
! m15_2p82d_knop_LAM_apu res pinned bit-identically at 1.5465E262
! for ~80 iterations and res2 jumped only on fresh-matrix calls, so
! the acceptance test tracked the stack layout, not the solve.
! Measure the same residual content where it lives at entry:
! per-ion SE%STEQ rows (row N_SE is the species conservation
! residual summed into STEQ_NC at EQ_SPECIES; the ionization
! residual at EQXzV folds in STEQ_ADV), plus STEQ_ED and
! STEQ_T/STEQ_T_EHB.
!
! Two deliberate differences from the assembled STEQ_VEC: rows
! FIXPOP_IN_BA_V3 would zero for held populations stay in the
! metric so LAMBDA and full calls are compared on the same
! physical quantities, and STEQ_ED/STEQ_T are skipped while NE/T
! are held fixed -- a residual the step was not asked to solve
! must not drive the acceptance test.
!
! The rollback mask is per-call state: reset before the pending
! test so a mask left over from an earlier masked rollback can
! never zero a column of an unrelated solve.
!
	  IF(ALLOCATED(BA_RBK_MASK))THEN
	    IF(SIZE(BA_RBK_MASK) .EQ. ND)BA_RBK_MASK(1:ND)=.FALSE.
	  END IF
	  BA_PART_THIS=.FALSE.
	  BA_SE_RES=0.0_LDP
!
! Per-call adjudication telemetry for the residual-driven recall
! policy [DO_BA_RES_REUSE] in SOLVE_FOR_POPS: reset each call so a
! deferred judgment (GS_POPS_HELD) or a call with no armed step
! leaves the flags false rather than replaying an old verdict.
!
	  BA_JUDGED=.FALSE.
	  BA_JUDGE_RATIO=1.0_LDP
	  RES_LOC(1)=1
	  RES_LOC(2)=1
	  COL_RES(1:ND)=0.0_LDP
	  DO ISPEC=1,NUM_SPECIES
	    DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	      DO I=1,ND
	        J=MAXLOC(ABS(SE(ID)%STEQ(1:SE(ID)%N_SE,I)),1)
	        T1=ABS(SE(ID)%STEQ(J,I))
	        IF(J .EQ. SE(ID)%N_SE)THEN
	          JJ=EQ_SPECIES(ISPEC)
	        ELSE
	          JJ=SE(ID)%EQ_IN_BA(J)
	        END IF
	        IF(ASSOCIATED(SE(ID)%STEQ_ADV))THEN
	          IF(ABS(SE(ID)%STEQ_ADV(I)) .GT. T1)THEN
	            T1=ABS(SE(ID)%STEQ_ADV(I))
	            JJ=ATM(ID)%EQXzV
	          END IF
	        END IF
!
! Per-equation quarantine [DO_BA_QUAR]: when the winning row of this
! ion column is parked, rescan the rows skipping quarantined
! equations -- a parked equation must not drive COL_RES and veto the
! acceptance test.  The fast MAXLOC path is kept whenever the winner
! is free: a quarantined non-maximum row cannot change the column
! result.  ED/T rows (NT-1,NT) are never parked and need no check.
! The mask index is tested in a nested IF: Fortran does not promise
! short-circuit .AND. evaluation, so BA_QUAR_MASK must not be
! indexed while JJ is unvalidated.
!
	        IF(LQUAR_OK .AND. JJ .GE. 1 .AND. JJ .LE. NT-2)THEN
	        IF(BA_QUAR_MASK(JJ,I))THEN
	          T1=0.0_LDP
	          DO K=1,SE(ID)%N_SE
	            IF(K .EQ. SE(ID)%N_SE)THEN
	              JK=EQ_SPECIES(ISPEC)
	            ELSE
	              JK=SE(ID)%EQ_IN_BA(K)
	            END IF
	            IF(JK .GE. 1 .AND. JK .LE. NT-2)THEN
	              IF(BA_QUAR_MASK(JK,I))CYCLE
	            END IF
	            IF(ABS(SE(ID)%STEQ(K,I)) .GT. T1)THEN
	              T1=ABS(SE(ID)%STEQ(K,I))
	              JJ=JK
	            END IF
	          END DO
	          IF(ASSOCIATED(SE(ID)%STEQ_ADV))THEN
	            JK=ATM(ID)%EQXzV
	            IF(JK .LT. 1 .OR. JK .GT. NT-2)JK=0
! Fortran does not promise short-circuit .OR. evaluation, so
! the BA_QUAR_MASK index must be staged behind the JK=0 test:
! row 0 does not exist.
!
	            LADV_OK=(JK .EQ. 0)
	            IF(.NOT. LADV_OK)
	1          LADV_OK=.NOT. BA_QUAR_MASK(JK,I)
	            IF(LADV_OK)THEN
	              IF(ABS(SE(ID)%STEQ_ADV(I)) .GT. T1)THEN
	                T1=ABS(SE(ID)%STEQ_ADV(I))
	                JJ=ATM(ID)%EQXzV
	              END IF
	            END IF
	          END IF
	        END IF
	        END IF
!
	        IF(.NOT. (T1 .LE. HUGE(1.0_LDP)))T1=HUGE(1.0_LDP)
	        IF(T1 .GT. COL_RES(I))COL_RES(I)=T1
	        IF(T1 .GT. BA_SE_RES)THEN
	          BA_SE_RES=T1
	          RES_LOC(1)=JJ
	          RES_LOC(2)=I
	        END IF
	      END DO
	    END DO
	  END DO
	  IF(.NOT. MOD_FIXED_NE)THEN
	    DO I=1,ND
	      T1=ABS(STEQ_ED(I))
	      IF(.NOT. (T1 .LE. HUGE(1.0_LDP)))T1=HUGE(1.0_LDP)
	      IF(T1 .GT. COL_RES(I))COL_RES(I)=T1
	      IF(T1 .GT. BA_SE_RES)THEN
	        BA_SE_RES=T1
	        RES_LOC(1)=NT-1
	        RES_LOC(2)=I
	      END IF
	    END DO
	  END IF
	  IF(.NOT. MOD_FIXED_T)THEN
	    DO I=1,ND
	      IF(USE_ELEC_HEAT_BAL .AND. ALLOCATED(STEQ_T_EHB))THEN
	        T1=ABS(STEQ_T_EHB(I))
	      ELSE
	        T1=ABS(STEQ_T(I))
	      END IF
	      IF(.NOT. (T1 .LE. HUGE(1.0_LDP)))T1=HUGE(1.0_LDP)
	      IF(T1 .GT. COL_RES(I))COL_RES(I)=T1
	      IF(T1 .GT. BA_SE_RES)THEN
	        BA_SE_RES=T1
	        RES_LOC(1)=NT
	        RES_LOC(2)=I
	      END IF
	    END DO
	  END IF
! Scope: the running minima feed the BA_RES_HARD bounds of the
! armed acceptance test, so they only move on calls that actually
! adjudicate one (BAACC_PENDING still armed here -- it is cleared
! below).  LAMBDA/held calls do not tighten bounds later FULL
! calls inherit.
!
	  IF(BAACC_PENDING)THEN
	    IF(BA_SE_RES .LT. BAACC_RES_MIN)BAACC_RES_MIN=BA_SE_RES
	    IF(ALLOCATED(BAACC_COL_MIN))THEN
	      IF(SIZE(BAACC_COL_MIN) .EQ. ND)THEN
	        DO I=1,ND
	          IF(COL_RES(I) .LT. BAACC_COL_MIN(I))
	1            BAACC_COL_MIN(I)=COL_RES(I)
	        END DO
	      END IF
	    END IF
	  END IF
!
! edited by Devin, 30-Sep-2026: second channel blind to a persistent
! residual floor.  BA_SE_RES saturates at the degenerate inner-
! boundary equations (~1E262 on the Dessart run), so growth anywhere
! below the floor was invisible -- an accepted BA step collapsed OI/
! HeI populations at depths 1-16 while the max never moved (FreeNe
! blow-up).  RES2 is the K-th largest per-depth column max: it
! ignores the K-1 worst depth columns wherever they sit and tracks
! the largest residual of the remaining grid.  K=ND/10 (min 4) keeps
! the channel live on small grids without dropping a large fraction
! of the model.
!
! edited by Devin, 03-Oct-2026: rank a copy, not COL_RES itself --
! the masked acceptance test below reads COL_RES per column and the
! MAXLOC decimation used to erase the top entries in place.
!
	  K=MIN(ND,MAX(4,ND/10))
	  COL_RES_S(1:ND)=COL_RES(1:ND)
	  DO J=1,K-1
	    I=MAXLOC(COL_RES_S,1)
	    COL_RES_S(I)=-HUGE(1.0_LDP)
	  END DO
	  BA_SE_RES2=MAXVAL(COL_RES_S)
	  RES2_COL=MAXLOC(COL_RES_S,1)
!
! RES2 minima obey the same pending gate as the RES/COL minima
! above: a ratchet on an unjudged call (LAMBDA, demoted, held)
! tightens the BA_RES_HARD bound at :534 on a value no acceptance
! test produced, so the next legitimate FULL call trips a phantom
! rejection and the run settles into forced RA_LAMBDA alternation.
!
	  IF(BAACC_PENDING)THEN
	    IF(.NOT. BAACC_R2GATE_PRINTED .AND.
	1      CONV_DIAG_LVL .GE. 2)THEN
	      WRITE(LUER,'(A)')
	1   '    BA_RES_ACC: RES2-min ratchet is now acceptance-gated'
	      BAACC_R2GATE_PRINTED=.TRUE.
	    END IF
	    IF(BA_SE_RES2 .LT. BAACC_RES2_MIN)
	1      BAACC_RES2_MIN=BA_SE_RES2
	  END IF
!
! Per-depth export for the ITER_SUM residual profile (15-Oct-2026).
!
	  IF(.NOT. ALLOCATED(BA_COL_RES))THEN
	    ALLOCATE(BA_COL_RES(ND),STAT=IOS)
	  ELSE IF(SIZE(BA_COL_RES) .NE. ND)THEN
	    DEALLOCATE(BA_COL_RES)
	    ALLOCATE(BA_COL_RES(ND),STAT=IOS)
	  END IF
	  IF(IOS .EQ. 0)BA_COL_RES(1:ND)=COL_RES(1:ND)
!
	  IF(DO_BA_RES_ACC)THEN
!
! edited by Devin, 29-Sep-2026: residual trajectory every call, not
! just on rejection -- needed to see whether accepted steps still
! grow the residual below the BA_RES_GROW/BA_RES_HARD thresholds.
! RES_LOC names the (variable,depth) with the largest residual --
! if one equation dominates, the residual acceptance criterion is
! measuring that equation, not the model as a whole.
!
	  IF(CONV_DIAG_LVL .GE. 1)THEN
	    WRITE(LUER,'(/,A,I5,A)')' ---------------------- BA_RES_ACC it',
	1    MAIN_COUNTER,' -------------------------------------'
	    WRITE(LUER,'(A,ES12.4)')
	1    '    max SE residual, current call   :',BA_SE_RES
	    WRITE(LUER,'(A,ES12.4)')
	1    '    residual before last full step  :',BAACC_RES_PREV
	    WRITE(LUER,'(A,ES12.4)')
	1    '    best residual seen              :',BAACC_RES_MIN
	    WRITE(LUER,'(A,L2)')
	1    '    acceptance test armed this call :',BAACC_PENDING
	    WRITE(LUER,'(A,L2)')
	1    '    lambda call                     :',LAMBDA_IT
	    WRITE(LUER,'(A,ES12.4)')
	1    '    correction cap, halves on reject:',BA_TRUST_CUR
	    WRITE(LUER,'(A,ES12.4)')
	1    '    residual excl worst columns     :',BA_SE_RES2
	    WRITE(LUER,'(A,I6)')
	1    '    depth of that residual          :',RES2_COL
	    WRITE(LUER,'(A,2I6,2A)')
	1    '    worst eqn (var,depth)           :',RES_LOC(1),
	1    RES_LOC(2),'  ',TRIM(EQ_DESC(RES_LOC(1)))
	    IF(BA_HOLD_RESTORE_REVIVAL .GT. 0_8)WRITE(LUER,'(A,I6)')
	1    '    held-snapshot restore revivals  :',
	2    BA_HOLD_RESTORE_REVIVAL
	  END IF
! edited by Devin, 29-Sep-2026: a non-finite residual (NaN or Inf in
! STEQ) fails every comparison below and so was silently accepted --
! the pending flag cleared and the trust radius grew after a poisoned
! step.  .NOT.(X .LE. HUGE) catches both Inf and NaN and counts them
! as rejections.  The same guard is applied to the RES2 channel
! (30-Sep-2026): the COL_RES clamp above already maps non-finite
! columns to HUGE, so this is belt-and-braces against any future
! change that lets an unclamped value through.
!
! An ordinary BA step followed by Lambda is not rolled back: the
! same correction can recur and form a FULL/LAMBDA cycle. A rejected
! BA_PROBE_LAM trial is different because it is a bounded measurement
! made inside RD_LAMBDA; if it degrades, the full snapshot is restored
! and this follow-up solve is skipped so its stale residual is not used.
! Other LAMBDA follow-ups are still evaluated but do not count as an
! acceptance when they show degradation.
!
! edited by Devin, 03-Oct-2026: judge the armed step per depth
! column.  Column I fails when its entry residual exceeds
! BA_RES_GROW times its own pre-step value or BA_RES_HARD times its
! own best-seen value; the BA_RES_SAT exemption applies per column
! (a column whose equation was already degenerate cannot fail on
! growth), and a column at the clamp value -- nonfinite entries map
! to HUGE above -- always fails the 0.5*HUGE floor even when its
! stored previous value was itself clamped.
!
! edited by Devin, 07-Oct-2026: skip adjudication while the caller
! discards this call's correction.  Under a population hold
! (GS_POPS_HELD, set by SOLVE_FOR_POPS before the call) POPS is
! restored to the pre-solve snapshot after SOLVEBA returns, so a
! rejection judged here acts on a state the step cannot change:
! the armed residual re-measures at frozen populations, the same
! columns mask, and BA_NREJ_CONSEC burns BA_NREJ_LIM on calls that
! could never pass (m15_2p82d_jfnk_cpl it71-75: five guaranteed
! rejections at a pinned 1.4338E11 residual killed the run).
! Defer judgement -- pending stays armed, so the armed step is
! measured at the first call that can actually move the
! populations.  The residual diagnostics above still print.
!
	  IF(.NOT. GS_POPS_HELD)THEN
	  NFAIL=0
	  LCOL_OK=ALLOCATED(BAACC_COL_PREV) .AND.
	1    ALLOCATED(BAACC_COL_MIN) .AND. ALLOCATED(BA_RBK_MASK)
	  IF(LCOL_OK)LCOL_OK=(SIZE(BAACC_COL_PREV) .EQ. ND .AND.
	1    SIZE(BAACC_COL_MIN) .EQ. ND .AND.
	2    SIZE(BA_RBK_MASK) .EQ. ND)
	  IF(BAACC_PENDING .AND. LCOL_OK)THEN
	    DO I=1,ND
	      IF(BA_RES_SAT .GT. 0.0_LDP .AND.
	1        BAACC_COL_PREV(I) .GE. BA_RES_SAT)CYCLE
	      IF(COL_RES(I) .GT. 0.5_LDP*HUGE(1.0_LDP) .OR.
	1        COL_RES(I) .GT. BA_RES_GROW*BAACC_COL_PREV(I) .OR.
	2        COL_RES(I) .GT. BA_RES_HARD*BAACC_COL_MIN(I))THEN
	        BA_RBK_MASK(I)=.TRUE.
	        NFAIL=NFAIL+1
	      END IF
	    END DO
	  END IF
!
! REJ_BULK is the grid-wide failure signal: the RES2 channel is the
! K-th largest column residual (K=ND/10), so it can only fail when
! more than K columns degraded.  REJ_MAX is the legacy single-cell
! channel, kept as part of the fallback decision when the column
! arrays are unavailable or BA_REJ_FRAC <= 0 disables the mask.
! With the mask live, a lone-cell rebound is resolved by the mask
! itself -- a whole-state rollback there would throw away a mostly
! good step.
!
	  REJ_BULK=.NOT. (BA_SE_RES2 .LE. HUGE(1.0_LDP)) .OR.
	1    BA_SE_RES2 .GT. HUGE(1.0_LDP)*0.5_LDP .OR.
	2    BA_SE_RES2 .GT. BA_RES_GROW*BAACC_RES2_PREV .OR.
	3    BA_SE_RES2 .GT. BA_RES_HARD*BAACC_RES2_MIN
	  REJ_MAX=(BA_RES_SAT .LE. 0.0_LDP .OR. BAACC_RES_PREV .LT.
	1    BA_RES_SAT) .AND.
	2    (.NOT. (BA_SE_RES .LE. HUGE(1.0_LDP)) .OR.
	3     BA_SE_RES .GT. BA_RES_GROW*BAACC_RES_PREV .OR.
	4     BA_SE_RES .GT. BA_RES_HARD*BAACC_RES_MIN)
	  REJ_STATE=.FALSE.
	  REJ_PART=.FALSE.
	  WLD_REJ=.FALSE.
	  IF(LCOL_OK .AND. BA_REJ_FRAC .GT. 0.0_LDP)THEN
	    REJ_STATE=REJ_BULK .OR.
	1      NFAIL .GT. MAX(1,NINT(BA_REJ_FRAC*DBLE(ND)))
	    REJ_PART=(.NOT. REJ_STATE) .AND. NFAIL .GT. 0
	  ELSE
	    REJ_STATE=REJ_MAX .OR. REJ_BULK
	  END IF
	  WLD_REJ=REJ_STATE .OR. REJ_PART
!
! Export the verdict of the armed step adjudicated this call for
! the residual-driven recall policy [DO_BA_RES_REUSE].  The RES2
! (rank) channel is the contraction measure: it is the bulk veto
! channel already, so a persistent degenerate floor cannot fake
! progress or fake a stall.  The column scan clamps nonfinite
! entries to HUGE, so a post-step residual at the clamp reads as a
! stall, as does growth from a zero pre-step base.
!
	  IF(BAACC_PENDING)THEN
	    BA_JUDGED=.TRUE.
	    IF(BA_SE_RES2 .GE. 0.5_LDP*HUGE(1.0_LDP))THEN
	      BA_JUDGE_RATIO=HUGE(1.0_LDP)
	    ELSE IF(BAACC_RES2_PREV .GT. 0.0_LDP)THEN
	      BA_JUDGE_RATIO=BA_SE_RES2/BAACC_RES2_PREV
	    ELSE IF(BA_SE_RES2 .GT. 0.0_LDP)THEN
	      BA_JUDGE_RATIO=HUGE(1.0_LDP)
	    ELSE
	      BA_JUDGE_RATIO=0.0_LDP
	    END IF
	  END IF
!
! edited by Devin, 07-Oct-2026: progress exemption for the
! consecutive-rejection streak.  A masked or soft-judged step
! that still lowered either residual channel is the partially-
! committed case the mask exists for (m15_2p82d_jfnk_cpl: every
! armed step it39-70 improved the max residual yet each burned
! the counter on the same 3-4 stiff inner columns, and the
! halved trust cap then zeroed most corrections).  L_PROG marks
! such a step: the streak resets instead of growing toward
! BA_NREJ_LIM, and the trust cap is left unchanged.  A masked
! step with no channel improvement still counts, so the kill
! switch survives a genuine stall.  The saturated-max exemption
! mirrors REJ_MAX.
!
	  L_PROG=BAACC_PENDING .AND.
	1    ((BA_SE_RES2 .LT. BAACC_RES2_PREV) .OR.
	2    ((BA_RES_SAT .LE. 0.0_LDP .OR.
	3      BAACC_RES_PREV .LT. BA_RES_SAT) .AND.
	4     BA_SE_RES .LT. BAACC_RES_PREV))
!
! A bounded RD_LAMBDA probe has a dedicated rollback contract only
! when the failure is grid-wide (REJ_STATE: the RES2 channel degraded
! or more than BA_REJ_FRAC of columns failed).  Restore the complete
! pre-probe state and skip this call's population solve: the residual
! and radiation stores were evaluated at the probe populations and
! would mix states.  When only a few columns fail (REJ_PART) the
! ordinary masked path below is now allowed to keep the good part of
! the probe -- the m15 SOB probes repeatedly degraded a single
! innermost column ~20x while improving the other ~69 columns ~45%,
! and the blanket full rollback discarded that progress every
! BA_PROBE_LAM calls.  At a masked column the pre-step populations
! are restored and that column's fresh STEQ correction is zeroed at
! application below, so a boundary column (ND or 1) is simply held at
! its pre-probe state for one call -- the per-depth SE system has no
! stencil to violate.  Ordinary BA/Lambda follow-ups retain the
! established soft-reject behaviour to avoid the full/lam cycle.
! The flag is armed only by BA_PROBE_LAM on non-monotonic paths.
!
! edited by Devin, 13-Oct-2026: gate on REJ_STATE, not WLD_REJ.
!
	  IF(BAACC_PENDING .AND. LAMBDA_IT .AND.
	1     BAACC_RDLAM_PROBE .AND. REJ_STATE)THEN
	    REJ_STATE=.TRUE.
	    REJ_PART=.FALSE.
	  END IF
!
	  IF(BAACC_PENDING .AND. LAMBDA_IT .AND.
	1     .NOT. (BAACC_RDLAM_PROBE .AND. WLD_REJ))THEN
!
! Honest accounting on a Lambda follow-up: the armed step is judged
! but not undone (see above).  A measured degradation halves the
! trust radius and counts as a soft rejection; previously the flag
! was cleared here untested AND accepted -- the trust radius grew
! and BA_NREJ_CONSEC reset after a step that may have wrecked the
! residual, so repeated bad full steps were invisible.
!
	    IF(WLD_REJ .AND. .NOT. L_PROG)THEN
	      BA_NREJ_SOFT=BA_NREJ_SOFT+1
	      BA_SOFT_THIS=.TRUE.
	      BA_TRUST_CUR=MAX(1.0_LDP,0.5_LDP*BA_TRUST_CUR)
! edited by Devin, 04-Oct-2026: count the soft rejection toward
! BA_NREJ_CONSEC.  In the alternating FULL/LAMBDA regime every
! rejected armed step was judged on a LAMBDA follow-up, so it was
! committed unconditionally while nrej stayed 0 and the run
! reported [OK] (jfnk_opt: 20 soft rejects by it70).  The step is
! still not undone -- the trust cap collapses toward its floor --
! but the same BA_NREJ_LIM bound that stops a same-call rejection
! run now also bounds a run of committed bad steps.
	      BA_NREJ_CONSEC=BA_NREJ_CONSEC+1
	      IF(CONV_DIAG_LVL .GE. 1)WRITE(LUER,'(A,I5,A,I4)')
	1        ' BA_RES_ACC it',MAIN_COUNTER,
	2        ': armed step degraded residual on LAMBDA'//
	3        ' follow-up (no rollback); soft rejects=',BA_NREJ_SOFT
	      IF(BA_NREJ_LIM .GT. 0 .AND.
	1         BA_NREJ_CONSEC .GE. BA_NREJ_LIM)THEN
	        WRITE(LUER,'(A,I4,A)')'BA_RES_ACC: ',BA_NREJ_CONSEC,
	1        ' consecutive rejections reach BA_NREJ_LIM -- stopping'
	        STOP
	      END IF
	    ELSE
!
! A judged-but-progressing step (WLD_REJ with L_PROG) resets the
! streak like a clean accept, but the trust radius stays put --
! the step was not clean.
!
	      BA_NREJ_CONSEC=0
	      IF(.NOT. WLD_REJ)
	1      BA_TRUST_CUR=MIN(BA_TRUST,1.25_LDP*BA_TRUST_CUR)
	    END IF
	    BAACC_PENDING=.FALSE.
	    BAACC_RDLAM_PROBE=.FALSE.
	  ELSE IF(BAACC_PENDING .AND. .NOT. REJ_STATE .AND.
	1    .NOT. REJ_PART)THEN
	    BAACC_PENDING=.FALSE.
	    BAACC_RDLAM_PROBE=.FALSE.
	    BA_NREJ_CONSEC=0
	    BA_TRUST_CUR=MIN(BA_TRUST,1.25_LDP*BA_TRUST_CUR)
	  ELSE IF(BAACC_PENDING .AND. REJ_PART)THEN
!
! Masked rollback: restore only the failed depth columns.  The
! rest of the rejected step is kept -- per-depth SE consistency
! makes the column, not the individual level, the safe rollback
! unit (rolling back one level would tear the column's ionization
! balance).  STEQ for the masked columns is zeroed after the solve
! below: it was linearized at the rejected populations.  Pending
! clears -- the step was adjudicated -- and the call continues to
! a normal solve instead of paying a wasted iteration.
!
	    DO I=1,ND
	      IF(BA_RBK_MASK(I) .AND. ALLOCATED(POPS_BA_SAVE))THEN
	        IF(SIZE(POPS_BA_SAVE,1) .EQ. NT .AND.
	1          SIZE(POPS_BA_SAVE,2) .EQ. ND)
	2          POPS(1:NT,I)=POPS_BA_SAVE(1:NT,I)
	      END IF
	    END DO
	    IF(DO_NONMONO_KNOP)THEN
!
! The lagged I_NEXT store was solved against the rejected pops at
! the masked columns; restoring the whole store is the consistent
! choice -- it is lagged by construction, so one step of staleness
! at the kept columns is the normal operating condition, while a
! store seeded at rejected populations is the documented wreck
! path (LAM_dbg it49-56).
!
	      CALL RESTORE_FG_INEXT_LATEST(FG_RB_OK)
	      IF(FG_RB_OK)WRITE(LUER,'(A)')
	1    'BA_RES_ACC: lagged FG field rolled back to pre-step state'
	    END IF
	    IF(DO_NONMONO_KNOP .OR. DO_NONMONO_SOB)THEN
	      CALL RESTORE_MOM_JNEXT_LATEST(MOM_RB_OK)
	      IF(MOM_RB_OK)WRITE(LUER,'(A)')
	1    'BA_RES_ACC: lagged MOM field rolled back to pre-step state'
	    END IF
	    BA_PART_THIS=.TRUE.
	    BA_NREJ=BA_NREJ+1
	    IF(L_PROG)THEN
	      BA_NREJ_CONSEC=0
	    ELSE
	      BA_NREJ_CONSEC=BA_NREJ_CONSEC+1
	      BA_TRUST_CUR=MAX(1.0_LDP,0.5_LDP*BA_TRUST_CUR)
	    END IF
	    WRITE(LUER,'(/,A,I5,A)')
	1      ' ---------------- BA_RES_ACC PARTIAL REJ it',
	2      MAIN_COUNTER,' ----------------------------'
	    WRITE(LUER,'(A,I6,A,I6)')
	1      '    masked depth columns         :',NFAIL,' of ',ND
	    WRITE(LUER,'(A,ES14.4)')
	1      '    SE residual before armed step:',BAACC_RES_PREV
	    WRITE(LUER,'(A,ES14.4)')
	1      '    SE residual measured this call:',BA_SE_RES
	    WRITE(LUER,'(A,ES12.4)')
	1      '    correction cap now           :',BA_TRUST_CUR
	    IF(L_PROG)WRITE(LUER,'(A)')
	1      '    residual improved; masked event not counted'
	    IF(BA_NREJ_CONSEC .GE. 3)WRITE(LUER,'(A,I3)')
	1      'BA_RES_ACC: consecutive rejections =',BA_NREJ_CONSEC
	    IF(BA_NREJ_LIM .GT. 0 .AND. BA_NREJ_CONSEC .GE. BA_NREJ_LIM)
	1    THEN
	      WRITE(LUER,'(A,I4,A)')'BA_RES_ACC: ',BA_NREJ_CONSEC,
	1      ' consecutive rejections reach BA_NREJ_LIM -- stopping'
	      STOP
	    END IF
	    BAACC_PENDING=.FALSE.
	    BAACC_RDLAM_PROBE=.FALSE.
!
! edited by Devin, 07-Oct-2026: quarantine each masked column's worst
! residual equation ([DO_BA_QUAR]).  A column that failed residual
! acceptance has at least one equation the linearized step made
! worse; parking it lets the column's other variables take their
! correction next armed solve instead of re-failing on the same
! equation every probe cycle.  SE%STEQ still holds the entry
! residual here (the solve has not run yet), so the scan identifies
! the equation that actually failed.  ED/T rows are excluded inside
! BA_QUAR_WORST_EQ.
!
	    IF(LQUAR_OK)THEN
	      DO I=1,ND
	        IF(.NOT. BA_RBK_MASK(I))CYCLE
	        IF(BA_QUAR_COUNT .GE. BA_QUAR_MAX)EXIT
	        JK=BA_QUAR_WORST_EQ(I)
	        IF(JK .GT. 0)THEN
	          BA_QUAR_MASK(JK,I)=.TRUE.
	          BA_QUAR_COUNT=BA_QUAR_COUNT+1
	          NQUAR_NEW=NQUAR_NEW+1
	          IF(CONV_DIAG_LVL .GE. 2)WRITE(LUER,'(A,I6,A,I4,2A)')
	1            'BA_RES_ACC: masked col quarantines eqn (var,depth)',
	2          JK,' @',I,'  ',TRIM(EQ_DESC(JK))
	        END IF
	      END DO
	    END IF
!
! Re-arming below copies COL_RES into BAACC_COL_PREV; at a masked
! column the measured value belongs to the rejected pops, so keep
! the pre-step residual as the best estimate of the restored state.
!
	    DO I=1,ND
	      IF(BA_RBK_MASK(I))COL_RES(I)=BAACC_COL_PREV(I)
	    END DO
	  ELSE IF(BAACC_PENDING)THEN
	    BA_STEP_REJECTED=.TRUE.
	    IF(BAACC_RDLAM_PROBE)WRITE(LUER,'(A)')
	1    'BA_RES_ACC: rejected rd_lam probe; restoring pre-probe state'
	    BA_NREJ=BA_NREJ+1
	    BA_NREJ_CONSEC=BA_NREJ_CONSEC+1
	    BA_TRUST_CUR=MAX(1.0_LDP,0.5_LDP*BA_TRUST_CUR)
!
! Snapshot-source contract: POPS_BA_SAVE is only refreshed on
! non-LAMBDA calls whose POPS equal the last committed state
! (solve_for_pops restores POPS_GS_SAVE after every held call,
! :1401-1426), so inside a held sequence the two snapshots are
! identical and either would do.  The post-:518 deferral means a
! restore can only run on an unheld call, but if the hold flag is
! ever armed while adjudication runs (contract violation), the
! caller's mandated restore target POPS_GS_SAVE is the fresher,
! committed-side state and wins.  The swap is counted so the
! event is loud rather than silently absorbed.
!
	    IF(GS_POPS_HELD .AND. ALLOCATED(POPS_GS_SAVE) .AND.
	1      SIZE(POPS_GS_SAVE,1) .EQ. NT .AND.
	2      SIZE(POPS_GS_SAVE,2) .EQ. ND)THEN
	      POPS=POPS_GS_SAVE
	      BA_HOLD_RESTORE_REVIVAL=BA_HOLD_RESTORE_REVIVAL+1_8
	      WRITE(LUER,'(A)')
	1   'BA_RES_ACC: held-sequence restore revived from'//
	2   ' POPS_GS_SAVE (snapshot overlap)'
	    ELSE IF(ALLOCATED(POPS_BA_SAVE))THEN
	      IF(SIZE(POPS_BA_SAVE,1) .EQ. NT .AND.
	1         SIZE(POPS_BA_SAVE,2) .EQ. ND)THEN
	        POPS=POPS_BA_SAVE
	      END IF
	    END IF
	    IF(DO_NONMONO_KNOP)THEN
!
! Roll the lagged radiation field back with the populations: the
! in-memory I_NEXT was solved this call on the post-step (rejected)
! populations, and re-solving from that seed can leave the GS basin
! entirely (LAM_dbg it49-56: pops rollback alone left se_res2
! growing 439 -> 2e7 at held pops).  FGINEXT_GS on disk still holds
! the generation written at the end of the rejected step's own
! iteration -- the field solved on the pre-step populations this
! block just restored -- so the pair stays consistent.  This call's
! own SCR_RITE has not run yet and cannot have overwritten it.
!
	      CALL RESTORE_FG_INEXT_LATEST(FG_RB_OK)
	      IF(FG_RB_OK)WRITE(LUER,'(A)')
	1    'BA_RES_ACC: lagged FG field rolled back to pre-step state'
	    END IF
	    IF(DO_NONMONO_KNOP .OR. DO_NONMONO_SOB)THEN
	      CALL RESTORE_MOM_JNEXT_LATEST(MOM_RB_OK)
	      IF(MOM_RB_OK)WRITE(LUER,'(A)')
	1    'BA_RES_ACC: lagged MOM field rolled back to pre-step state'
	    END IF
	    WRITE(LUER,'(/,A,I5,A)')
	1      ' -------------------- BA_RES_ACC REJ it',MAIN_COUNTER,
	2      ' -----------------------------------'
	    WRITE(LUER,'(A,ES14.4)')
	1      '    SE residual before rejected step :',BAACC_RES_PREV
	    WRITE(LUER,'(A,ES14.4)')
	1      '    SE residual measured this call   :',BA_SE_RES
	    WRITE(LUER,'(A,ES12.4)')
	1      '    correction cap after halving     :',BA_TRUST_CUR
	    IF(BA_NREJ_CONSEC .GE. 3)WRITE(LUER,'(A,I3)')
	1      'BA_RES_ACC: consecutive rejections =',BA_NREJ_CONSEC
!
! edited by Devin, 29-Sep-2026: unbounded consecutive rejections
! mean the linearized model cannot produce an acceptable step at
! any trust radius -- BA_TRUST_CUR floors at 1.0, so the trust-
! region mechanism on its own alternates LAMBDA and rejected BA
! steps forever and the job finishes un-converged while reporting
! as a completed model.  Stop loudly instead; BA_NREJ_LIM <= 0
! disables the cap.
!
	    IF(BA_NREJ_LIM .GT. 0 .AND. BA_NREJ_CONSEC .GE. BA_NREJ_LIM)
	1    THEN
	      WRITE(LUER,'(A,I4,A)')'BA_RES_ACC: ',BA_NREJ_CONSEC,
	1      ' consecutive rejections reach BA_NREJ_LIM -- stopping'
	      STOP
	    END IF
	    BAACC_PENDING=.FALSE.
	    BAACC_RDLAM_PROBE=.FALSE.
	    STEQ=0.0_LDP
	    SUCCESS=.TRUE.
! edited by Devin, 29-Sep-2026: cap the reported MAXCH -- an exploded
! residual (seen at 1E262 on the Dessart run) otherwise exceeds the
! MAX_CHNG guard and terminates the model with a misleading "bad
! initial populations" STOP instead of letting the trust radius
! contract and retry.  1.5E5 (not 1.0E5) keeps a worst-case rejection
! strictly above the 1.0E5 stay-in-Lambda guard in SOLVE_FOR_POPS;
! a cap at exactly the guard value would slip past its .GT. test.
! A non-finite residual skips the product so MIN/MAX never see NaN.
!
	    IF(BA_SE_RES .LE. HUGE(1.0_LDP))THEN
	      MAXCH=MAX(1.0E+02_LDP,MIN(1.0E+02_LDP*BA_SE_RES,1.5E+05_LDP))
	    ELSE
	      MAXCH=1.5E+05_LDP
	    END IF
! The KNOP controller reads MAXCH_APP; without this assignment a
! rejected step would report the previous call's committed change
! and could keep BA mode when it should demote.
	    MAXCH_APP=MAXCH
	    RETURN
	  END IF
	  IF(.NOT. LAMBDA_IT .OR. BA_SNAPSHOT_MISSING)THEN
	    IF(.NOT. ALLOCATED(POPS_BA_SAVE))THEN
	      ALLOCATE(POPS_BA_SAVE(NT,ND),STAT=IOS)
	    ELSE IF(SIZE(POPS_BA_SAVE,1) .NE. NT .OR.
	1       SIZE(POPS_BA_SAVE,2) .NE. ND)THEN
	      DEALLOCATE(POPS_BA_SAVE)
	      ALLOCATE(POPS_BA_SAVE(NT,ND),STAT=IOS)
	    ELSE
!
! No alloc ran: the stored IOS is stale and must not decide whether
! the snapshot is armed -- a leftover nonzero here reported "alloc
! failed" on every later call without ever attempting one, so the
! acceptance net stayed permanently disarmed (LAM_dbg it39-48).
!
	      IOS=0
	    END IF
!
! edited by Devin, 03-Oct-2026: the masked acceptance test needs a
! per-depth residual history and mask alongside the snapshot.  Each
! array gets the same alloc/realloc treatment; a failure anywhere in
! the chain leaves IOS nonzero so BA_SNAPSHOT_MISSING still reports
! the net as disarmed -- the masked test then falls back to the
! legacy scalar decision (LCOL_OK false).
!
	    IF(IOS .EQ. 0)THEN
	      IF(.NOT. ALLOCATED(BAACC_COL_PREV))THEN
	        ALLOCATE(BAACC_COL_PREV(ND),STAT=IOS)
	      ELSE IF(SIZE(BAACC_COL_PREV) .NE. ND)THEN
	        DEALLOCATE(BAACC_COL_PREV)
	        ALLOCATE(BAACC_COL_PREV(ND),STAT=IOS)
	      END IF
	    END IF
	    IF(IOS .EQ. 0)THEN
	      IF(.NOT. ALLOCATED(BAACC_COL_MIN))THEN
	        ALLOCATE(BAACC_COL_MIN(ND),STAT=IOS)
	        IF(IOS .EQ. 0)BAACC_COL_MIN(1:ND)=COL_RES(1:ND)
	      ELSE IF(SIZE(BAACC_COL_MIN) .NE. ND)THEN
	        DEALLOCATE(BAACC_COL_MIN)
	        ALLOCATE(BAACC_COL_MIN(ND),STAT=IOS)
	        IF(IOS .EQ. 0)BAACC_COL_MIN(1:ND)=COL_RES(1:ND)
	      END IF
	    END IF
	    IF(IOS .EQ. 0)THEN
	      IF(.NOT. ALLOCATED(BA_RBK_MASK))THEN
	        ALLOCATE(BA_RBK_MASK(ND),STAT=IOS)
	      ELSE IF(SIZE(BA_RBK_MASK) .NE. ND)THEN
	        DEALLOCATE(BA_RBK_MASK)
	        ALLOCATE(BA_RBK_MASK(ND),STAT=IOS)
	      END IF
	      IF(IOS .EQ. 0)BA_RBK_MASK(1:ND)=.FALSE.
	    END IF
!
! edited by Devin, 07-Oct-2026: [DO_BA_QUAR] storage.  The mask is
! persistent state -- quarantine must survive between armed calls,
! unlike the per-call BA_RBK_MASK -- so only a fresh or resized
! allocation is initialized.
!
	    IF(IOS .EQ. 0 .AND. DO_BA_QUAR)THEN
	      IF(.NOT. ALLOCATED(BA_QUAR_MASK))THEN
	        ALLOCATE(BA_QUAR_MASK(NT,ND),STAT=IOS)
	        IF(IOS .EQ. 0)BA_QUAR_MASK(1:NT,1:ND)=.FALSE.
	      ELSE IF(SIZE(BA_QUAR_MASK,1) .NE. NT .OR.
	1        SIZE(BA_QUAR_MASK,2) .NE. ND)THEN
	        DEALLOCATE(BA_QUAR_MASK)
	        ALLOCATE(BA_QUAR_MASK(NT,ND),STAT=IOS)
	        IF(IOS .EQ. 0)BA_QUAR_MASK(1:NT,1:ND)=.FALSE.
	      END IF
	    END IF
	    IF(IOS .EQ. 0)THEN
	      POPS_BA_SAVE=POPS
	      BA_SNAPSHOT_MISSING=.FALSE.
! Only non-LAMBDA steps arm the acceptance check.  A large first
! BA correction is normal (Dessart M15 epoch 1 applies thousands of
! percent and the luminosity stays on the scale of LSTAR).  Arming
! the rollback on a Lambda call would undo that class of step.
! The re-arm attempt still runs on forced-LAMBDA calls so a
! transient alloc failure can heal.
	      IF(.NOT. LAMBDA_IT)THEN
	        BAACC_PENDING=.TRUE.
	        BAACC_RDLAM_PROBE=BAACC_RDLAM_PROBE_ARM
	        BAACC_RES_PREV=BA_SE_RES
	        BAACC_RES2_PREV=BA_SE_RES2
	        BAACC_COL_PREV(1:ND)=COL_RES(1:ND)
	      END IF
	    ELSE
	      BA_SNAPSHOT_MISSING=.TRUE.
	      WRITE(LUER,*)
	1      'BA_RES_ACC: snapshot alloc failed; acceptance disabled'
	    END IF
	  END IF
	  END IF
	END IF
!
	CALL SET_CASE_UP(SCALE_OPT,IONE,IZERO)
	IF(SCALE_OPT(1:5) .NE. 'LOCAL' .AND. SCALE_OPT(1:4) .NE. 'NONE'
	1   .AND. SCALE_OPT(1:6) .NE. 'GLOBAL'
	1   .AND. SCALE_OPT(1:5) .NE. 'MAJOR')THEN
	  SCALE_OPT='MAJOR'
	  WRITE(LUER,*)'Warning - Invalid scale option in SOLVEBA',
	1          ' MAJOR scaling assumed'
	END IF
!
	WRITE(LUER,'(/,A)')
	1  ' ------------------------ SOLVEBA POPULATION UPDATE ------------------------'
!
! Solve for the perturbations. ND: The scaling of BA is now done within
! CMF_BLKBAND.
!
	IF(METH_SOL(1:4) .EQ. 'DIAG' .OR. METH_SOL(1:3) .EQ. 'TRI')THEN
	  CALL TUNE(IONE,'BLKBAND')
!
! Perform the solution for depth L. The IONE refers to DIAG_INDX, NUM_BNDS,
! and ND respectively.
!
	    LOC_WR_BA_INV=WR_BA_INV
	    IF(LAMBDA_IT)LOC_WR_BA_INV=.FALSE.
	    CALL CMF_BLKBAND_V3(STEQ,POPS,METH_SOL,SUCCESS,
	1              DIAG_INDX,NT,NION,NUM_BNDS,ND,
	1              BA_COMPUTED,LOC_WR_BA_INV,WR_PRT_INV)
!
! edited by Devin, 30-Sep-2026: record the linear-solve path --
! newmat=F with wrinv=T means the correction came from a stored
! inverse applied to a new residual, i.e. a stale Jacobian.  That
! distinction decides whether the large corrections are a matrix
! problem or a genuine Newton overshoot.
!
	    IF(CONV_DIAG_LVL .GE. 2 .AND. (DO_NONMONO_SOB .OR.
	1       DO_NONMONO_KNOP .OR. DO_BA_RES_ACC))THEN
	      WRITE(LUER,'(A,I5,A,L2,A,A4,A,L2,A,L2,A,L2,A,L2)')
	1      ' BA-SOLVE: it=',MAIN_COUNTER,' lam=',LAMBDA_IT,
	2      ' meth=',METH_SOL,' newmat=',BA_COMPUTED,
	3      ' wrinv=',LOC_WR_BA_INV,' wrprt=',WR_PRT_INV,
	4      ' ok=',SUCCESS
	    END IF
	    IF(.NOT. SUCCESS)THEN
	      WRITE(LUER,*)'Error in CMF_BLKBAND_V3 - shutting code down'
	      STOP
	    END IF
!
	  CALL TUNE(ITWO,'BLKBAND')
!
	ELSE
	  WRITE(LUER,*)'Error - invalid solution method in SOLVEBA'
	  WRITE(LUER,*)'Solution method is ',METH_SOL
	  STOP
	END IF
!
! This is done to prevent corrections at depth 2 setting the nature of future iterations.
!
	IF(SET_POPS_D2_EQ_D1 .AND. .NOT. LAMBDA_IT)THEN
	  STEQ(:,2)=0.0_LDP
	END IF
!
! BA trust-region bound (VADAT BA_TRUST, only with DO_BA_RES_ACC).
! STEQ now holds the correction; a depth whose requested change
! exceeds the current radius BA_TRUST_CUR has a linearization that
! is invalid there (observed pathology: ~1e12 entries from a stale
! stored inverse).  Zeroing the whole column keeps that depth at
! its previous populations instead of applying a scaled version
! of a garbage direction; valid depths still take their
! corrections.  BA_TRUST_CUR halves on each rejected step and
! recovers 1.25x per accepted step, up to BA_TRUST.
!
	STEQ_MASKED_MAX=0.0_LDP
	L_ANY_MASK=.FALSE.
	IF(DO_BA_RES_ACC .AND. .NOT. LAMBDA_IT)THEN
	  COUNT(2)=0
	  COUNT(3)=0
	  LMASK_OK=.FALSE.
	  IF(ALLOCATED(BA_RBK_MASK))THEN
	    IF(SIZE(BA_RBK_MASK) .EQ. ND)LMASK_OK=.TRUE.
	  END IF
	  DO I=1,ND
!
! edited by Devin, 07-Oct-2026: per-equation quarantine
! ([DO_BA_QUAR]).  Parked equations are suppressed FIRST, before
! the trust test, so a parked variable's saturated demand can no
! longer zero the whole column.  Suppression applies on non-LAMBDA
! solves only: LAMBDA calls still update the level, so a parked
! variable keeps contracting under fixed-point iteration.  The
! suppressed demand is accumulated in STEQ_QMAX for the BA13-TRACE
! report but is deliberately kept out of MAXCH -- parking declares
! the equation out of play for the demand metric, and its
! Lambda-call correction still counts there.
!
	    IF(LQUAR_OK)THEN
	      DO J=1,NT-2
	        IF(BA_QUAR_MASK(J,I))THEN
	          STEQ_QMAX=MAX(STEQ_QMAX,ABS(STEQ(J,I)))
	          STEQ(J,I)=0.0_LDP
	        END IF
	      END DO
	    END IF
!
! Masked (rolled-back) columns are held at their restored
! populations this call: their fresh STEQ was linearized at the
! rejected populations and must not be applied.  The pre-mask
! demand is recorded so it can stay visible in MAXCH below.
!
	    IF(LMASK_OK .AND. BA_RBK_MASK(I))THEN
	      STEQ_MASKED_MAX=MAX(STEQ_MASKED_MAX,
	1        MAXVAL(ABS(STEQ(1:NT,I))))
	      STEQ(:,I)=0.0_LDP
	      L_ANY_MASK=.TRUE.
	      COUNT(3)=COUNT(3)+1
	    ELSE IF(MAXVAL(ABS(STEQ(1:NT,I))) .GT. BA_TRUST_CUR)THEN
	      STEQ_MASKED_MAX=MAX(STEQ_MASKED_MAX,
	1        MAXVAL(ABS(STEQ(1:NT,I))))
!
! edited by Devin, 07-Oct-2026: an over-trust column is still
! zeroed -- this call's direction was solved jointly and cannot be
! partially trusted -- but under [DO_BA_QUAR] the population
! variables carrying the over-trust demand are parked so the next
! armed solve can move the rest of the column.  On the m15 2.82d
! SOB run ~60 of 76 columns were being zeroed each armed solve
! because one equation per column demanded annihilation.
!
	      IF(LQUAR_OK .AND. BA_QUAR_COUNT .LT. BA_QUAR_MAX)THEN
	        DO J=1,NT-2
	          IF(ABS(STEQ(J,I)) .GT. BA_TRUST_CUR .AND.
	1            .NOT. BA_QUAR_MASK(J,I))THEN
	            BA_QUAR_MASK(J,I)=.TRUE.
	            BA_QUAR_COUNT=BA_QUAR_COUNT+1
	            NQUAR_NEW=NQUAR_NEW+1
	            IF(CONV_DIAG_LVL .GE. 2)
	1              WRITE(LUER,'(A,I6,A,I4,2A)')
	2              'BA_RES_ACC: over-trust demand quarantines eqn'//
	3              ' (var,depth)',J,' @',I,'  ',TRIM(EQ_DESC(J))
	            IF(BA_QUAR_COUNT .GE. BA_QUAR_MAX)EXIT
	          END IF
	        END DO
	      END IF
	      STEQ(:,I)=0.0_LDP
	      L_ANY_MASK=.TRUE.
	      COUNT(2)=COUNT(2)+1
	    END IF
	  END DO
	  IF(LQUAR_OK .AND. (NQUAR_NEW .NE. 0 .OR.
	1     STEQ_QMAX .GT. 0.0_LDP))THEN
	    WRITE(LUER,'(A,I5,A,I4,A,ES12.4)')
	1    'BA_RES_ACC: quarantined eqns (total/new):',BA_QUAR_COUNT,
	2    '/',NQUAR_NEW,'   max parked |STEQ|:',STEQ_QMAX
	  END IF
	  IF(COUNT(2) .NE. 0)THEN
	    WRITE(LUER,'(A,I6,A,ES12.4)')
	1      'BA_RES_ACC: zeroed ',COUNT(2),
	2      ' depth corrections with |STEQ|>',BA_TRUST_CUR
	  END IF
	  IF(COUNT(3) .NE. 0)THEN
	    WRITE(LUER,'(A,I6,A)')
	1      'BA_RES_ACC: masked ',COUNT(3),
	2      ' rolled-back depth corrections'
	  END IF
	END IF
!
! edited by Devin, 13-Oct-2026: a partially rejected rdlam probe is
! adjudicated on a LAMBDA call, so the .NOT. LAMBDA_IT block above is
! skipped and its mask-zeroing never runs.  The restored columns must
! still be held here: the Lambda correction at a masked column was
! solved against the rejected populations and would move the restored
! state along a stale direction.  Zero it for the same reason the
! non-LAMBDA block does; the demand is not recorded in
! STEQ_MASKED_MAX because it was measured at the rejected state.
! Gated on BA_PART_THIS, not just the mask: the NFAIL loop also marks
! columns on a soft-rejected ordinary armed step where no rollback
! happened, and holding those (unrestored) columns would discard a
! valid Lambda correction.
!
	IF(DO_BA_RES_ACC .AND. LAMBDA_IT .AND. BA_PART_THIS)THEN
	  LMASK_OK=.FALSE.
	  IF(ALLOCATED(BA_RBK_MASK))THEN
	    IF(SIZE(BA_RBK_MASK) .EQ. ND)LMASK_OK=.TRUE.
	  END IF
	  IF(LMASK_OK)THEN
	    COUNT(3)=0
	    DO I=1,ND
	      IF(BA_RBK_MASK(I))THEN
	        STEQ(:,I)=0.0_LDP
	        L_ANY_MASK=.TRUE.
	        COUNT(3)=COUNT(3)+1
	      END IF
	    END DO
	    IF(COUNT(3) .NE. 0)THEN
	      WRITE(LUER,'(A,I6,A)')
	1        'BA_RES_ACC: masked ',COUNT(3),
	2        ' rolled-back depth corrections'
	    END IF
	  END IF
	END IF
!
! 
!
! Lambda iterations (for Ne fixed) should generally have corrections < unity.
! Due to instabilities, large -ve corrections can sometimes arrise. This
! limits these corrections, and potentilly wil help facilitate convergence.
!
	IF(LAMBDA_IT .AND. LAM_SCALE_OPT(1:5) .EQ. 'LIMIT')THEN
	  COUNT(1)=0
	  DO I=1,ND
	    DO J=1,NT-1
	      IF(STEQ(J,I) .GT. 1.1_LDP)THEN
	        STEQ(J,I)=0.999_LDP
	        COUNT(1)=COUNT(1)+1
	      END IF
	    END DO
	  END DO
	  IF(COUNT(1) .NE. 0)
	1   WRITE(LUER,*)'Warning -- using  LIMIT option for LAMBDA iteration in SOLVEBA_V13.f'
	END IF
!
! Determine maximum corrections to the 'population parameters', and output
! summary file:
!
	LU_SUM=7
	CALL GEN_ASCI_OPEN(LU_SUM,'CORRECTION_SUM','UNKNOWN',' ',' ',IZERO,IOS)
	WRITE(LU_SUM,'(A)')' '
	WRITE(LU_SUM,'(A)')' Summary of changes at each depth '
	WRITE(LU_SUM,'(A,I6)')' NT=',NT
	WRITE(LU_SUM,'(A)')' '
! Cumulative demanded-correction counts over ALL NT rows per depth.
! The histogram is deliberately unmasked: gating it on the
! BA_DEM_PFLOOR floor hid that ~99% of rows demanded
! large corrections on the m15 formal KNOP start (it1-2, Oct-2026).
! The floor still gates the DECREASE/INCREASE scan below, which feeds
! the MAXCH control metric where unsatisfiable floor-level demands
! are noise.
!
	WRITE(LU_SUM,'(3X,A,7(1X,A))')'Depth',' 100.0%','  10.0%','   1.0%',
	1        '   0.1%','  0.01%',' 0.001%','0.0001%'
!
!   STEQ is +ve when the populations are to decrease.
!   STEQ is -ve when the populations are to increase.
!
	DECREASE=0.0_LDP
	IDEC=0
	INCREASE=0.0_LDP
	IINC=0
	JDEC=0
	JINC=0
	STEQ_EXM=0.0_LDP
	IEXM=0
	JEXM=0
	BA_DEM_EXM=0.0_LDP
	BA_TR_ABS=0.0_LDP
	EX_DEC=0.0_LDP
	EX_INC=0.0_LDP
	DO I=1,ND
	  COUNT(:)=0
!
! Removed 17-Oct-2026: the [BA_DEM_RELFRAC] relative demand floor and
! its POPS(J,I).LT.T2 skip below.  Relative exemption let hundreds of
! populated-depth trace rows escape the MAXCH metric (SOB l14 it242:
! ~700/depth >0.1% while MAXCH printed 0.08%).
!
! edited by Devin, 18-Oct-2026: [BA_EX_RELFRAC] exit-only relative
! floor.  T2EX uses the old RELFRAC predicate but feeds only
! EX_DEC/EX_INC, which become BA_EX_MAXCH -- the demand the
! converged-exit gate tests.  DECREASE/INCREASE (MAXCH) are
! untouched, so the hidden rows stay in the reported metric and
! their largest demand is exported separately as BA_DEM_EXR.
!
	  T2EX=0.0_LDP
	  IF(BA_EX_RELFRAC .GT. 0.0_LDP)
	1   T2EX=BA_EX_RELFRAC*MAXVAL(POPS(1:NT-2,I))
	  DO J=1,NT
	    T1=1.0E+06_LDP*ABS(STEQ(J,I))+1.0E-31_LDP		!to ensure non-zero.
	    K=LOG10(T1)+1; K=MIN(7,K)
	    IF(K > 0)COUNT(1:K)=COUNT(1:K)+1
!
! edited by Devin, 11-Oct-2026: [BA_DEM_PFLOOR] skip floor-level rows
! in the demanded-correction metric.  A population already below the
! floor demands relative corrections that carry no physical
! information (its rates/opacities are all negligible), yet an
! unsatisfiable re-fed demand pins MAXCH at the 1E7 sentinel every
! call -- m15 fresh run: CV:1s_4p_1Po at depth 16, n~1E-55, demanded
! decrease 1E7% on every non-LAMBDA call because cascade coupling
! refills it between BA solves.  ED/T rows (J=NT-1,NT) are exempt;
! corrections on skipped rows still commit normally.
!
	    IF(BA_DEM_PFLOOR .GT. 0.0_LDP .AND. J .LE. NT-2 .AND.
	1     POPS(J,I) .LT. BA_DEM_PFLOOR)CYCLE
	    IF(STEQ(J,I) .GT. DECREASE)THEN
	      DECREASE=STEQ(J,I)
	      IDEC=I
	      JDEC=J
	    ELSE IF(STEQ(J,I) .LT. INCREASE)THEN
	      INCREASE=STEQ(J,I)
	      IINC=I
	      JINC=J
	    END IF
	    IF(J .LE. NT-2 .AND. POPS(J,I) .LT. T2EX)CYCLE
	    IF(STEQ(J,I) .GT. EX_DEC)EX_DEC=STEQ(J,I)
	    IF(STEQ(J,I) .LT. EX_INC)EX_INC=STEQ(J,I)
	  END DO
	  WRITE(LU_SUM,'(8I8)')I,COUNT(7:1:-1)
	END DO
	CLOSE(LU_SUM)
!
! edited by Devin, 16-Oct-2026: the exempted rows' largest demand is
! exported as BA_DEM_EXM so a floor-hidden correction (e.g. a trace
! species below the relative floor) is still visible in OUTGEN and
! can gate the convergence exit in CMFGEN_SUB.  STEQ is fractional;
! the export uses the same percent scale as DECREASE/INCREASE.  A
! second scan is cheaper than folding into the masked loop above:
! the exempted set is exactly the rows the two floor tests skipped.
!
	DO I=1,ND
	  DO J=1,NT-2
	    IF(BA_DEM_PFLOOR .GT. 0.0_LDP .AND.
	1       POPS(J,I) .LT. BA_DEM_PFLOOR)THEN
	      IF(ABS(STEQ(J,I)) .GT. STEQ_EXM)THEN
	        STEQ_EXM=ABS(STEQ(J,I))
	        IEXM=I
	        JEXM=J
	      END IF
	    END IF
	  END DO
	END DO
	BA_DEM_EXM=100.0_LDP*STEQ_EXM
	IF(JEXM .GT. 0)THEN
	  WRITE(LUER,9215)IEXM,BA_DEM_EXM,TRIM(EQ_DESC(JEXM)),
	1   '          --- iteration ',MAIN_COUNTER
	END IF
!
! edited by Devin, 18-Oct-2026: [BA_EX_RELFRAC] BA_DEM_EXR is the
! largest demand over the rows the exit-only relative floor skips
! (the PFLOOR rows are already covered by BA_DEM_EXM).  Exported on
! the same percent scale and printed every call so the demand the
! converged exit hides is never silent.
!
	BA_DEM_EXR=0.0_LDP
	EXR_VAL=0.0_LDP
	IEXR=0
	JEXR=0
	IF(BA_EX_RELFRAC .GT. 0.0_LDP)THEN
	  DO I=1,ND
	    T2EX=BA_EX_RELFRAC*MAXVAL(POPS(1:NT-2,I))
	    DO J=1,NT-2
	      IF(POPS(J,I) .LT. T2EX .AND.
	1        (BA_DEM_PFLOOR .LE. 0.0_LDP .OR.
	2         POPS(J,I) .GE. BA_DEM_PFLOOR))THEN
	        IF(ABS(STEQ(J,I)) .GT. EXR_VAL)THEN
	          EXR_VAL=ABS(STEQ(J,I))
	          IEXR=I
	          JEXR=J
	        END IF
	      END IF
	    END DO
	  END DO
	  BA_DEM_EXR=100.0_LDP*EXR_VAL
	  IF(JEXR .GT. 0)THEN
	    WRITE(LUER,9217)IEXR,BA_DEM_EXR,TRIM(EQ_DESC(JEXR)),
	1     '      --- iteration ',MAIN_COUNTER
	  END IF
	END IF
9215	FORMAT(' Largest exempted demand at depth ',I4,' is',ES10.2,'%  (',
	1        A,')',A,I4)
9217	FORMAT(' Largest exit-exempted demand at depth ',I4,' is',ES10.2,
	1       '%  (',A,')',A,I4)
9216	FORMAT(' Stage-B trace |dn|/gas at depth ',I4,' is',ES10.2,'  (',
	1        A,')',A,I4)
!
	DECREASE=100.0_LDP*DECREASE
	INCREASE=100.0_LDP*INCREASE
	IF(LAMBDA_IT)THEN
 	  WRITE(LUER,9000)IINC,ABS(INCREASE),'  (LAMBDA)',
	1          '          --- iteration ',MAIN_COUNTER
	  WRITE(LUER,9200)IDEC,DECREASE,'  (LAMBDA)',
	1          '          --- iteration ',MAIN_COUNTER
	ELSE IF(BA_COMPUTED)THEN
 	  WRITE(LUER,9000)IINC,ABS(INCREASE),'  (BA computed)',
	1          '     --- iteration ',MAIN_COUNTER
	  WRITE(LUER,9200)IDEC,DECREASE,'  (BA computed)',
	1          '     --- iteration ',MAIN_COUNTER
	ELSE
 	  WRITE(LUER,9000)IINC,ABS(INCREASE),'  (BA not computed)',
	1          ' --- iteration ',MAIN_COUNTER
	  WRITE(LUER,9200)IDEC,DECREASE,'  (BA not computed)',
	1          ' --- iteration ',MAIN_COUNTER
	END IF
9000	FORMAT(' Maximum % increase at depth ',I4,' is',ES10.2,A,A,I4)
9200	FORMAT(' Maximum % decrease at depth ',I4,' is',ES10.2,A,A,I4)
!
! Convert decrease into more useful form. If DECREASE=100, the
! new population value is zero and change is infinite. i.e we
! convert decrease to form 100(old-new)/new. Currently in the form
! 100(old-new)/old [same as INCREASE]. We allow for big decreases -
! limit MAXCH to 1.0D+07 - This value does not halt program execution.
!
	MAXCH=-INCREASE
	IF(DECREASE .LT. 99.999_LDP)THEN
	  DECREASE=100.0_LDP*DECREASE/(100.0_LDP-DECREASE)
	ELSE
	  DECREASE=1.0E+07_LDP
	END IF
	MAXCH=MAX(MAXCH,DECREASE)
!
! edited by Devin, 04-Oct-2026: a column zeroed by the trust mask
! (or held by the roll-back mask) is a refused correction, not a
! satisfied one -- its demand must stay visible.  Without this a
! fully masked step reported MAXCH=0 and tripped the MAXCH<EPS
! converged exit in CMFGEN_SUB (jfnk_hyb stopped at it8 with
! 69/70 columns masked, se_res=1.1E14).  MAXCH_APP stays
! post-mask: it correctly reports what was actually applied.
!
	IF(STEQ_MASKED_MAX .GT. 0.0_LDP)THEN
	  MAXCH=MAX(MAXCH,MIN(1.0D+07,1.0D+02*STEQ_MASKED_MAX))
	END IF
!
!Determine the largest set of changes:
!We determine the largest NV increases, and the largest NV decreases.
!
	MAX_INC_VEC=0.0_LDP
	MAX_DEC_VEC=0.0_LDP
	DO J=1,ND
!
! edited by Devin, 13-Oct-2026: [BA_DEM_PFLOOR] must gate
! this scan too.  When DO_LEVEL_CHK is true the NV-th largest STEQ value
! replaces MAXCH below, so floor-level rows skipped by the demand scan
! would still pin the metric whenever they fill all NV slots (the X_SOB
! 40d run: HeI singlet dead levels at depth 73 hold MAXCH ~9% while the
! spectrum is converged).  ED/T rows (I>NT-2) stay exempt; the floor
! defaults to 0, which leaves this scan unchanged.
! (17-Oct-2026: the [BA_DEM_RELFRAC] second floor was removed here as
! everywhere -- see the module comment in CONTROL_VARIABLE_MOD.)
!
	  DO I=1,NT
	    IF(BA_DEM_PFLOOR .GT. 0.0_LDP .AND. I .LE. NT-2 .AND.
	1     POPS(I,J) .LT. BA_DEM_PFLOOR)CYCLE
	    T1=STEQ(I,J)
	    IF(T1 .LE. 0.0_LDP)THEN
	      IF(T1 .LT. MAX_INC_VEC(NV))THEN
	        DO K=NV-1,1,-1
	          IF(T1 .GT. MAX_INC_VEC(K))THEN
	            MAX_INC_VEC(K+2:NV)=MAX_INC_VEC(K+1:NV-1)
	            MAX_INC_VEC(K+1)=T1
	            EXIT
	          END IF
	          IF(K .EQ. 1)THEN
	            MAX_INC_VEC(2:NV)=MAX_INC_VEC(1:NV-1)
	            MAX_INC_VEC(1)=T1
	          END IF
	        END DO
	      END IF
	    ELSE IF(STEQ(I,J) .GT. 0.0_LDP)THEN
	      IF(T1 .GT. MAX_DEC_VEC(NV))THEN
	        DO K=NV-1,1,-1
	          IF(T1 .LT. MAX_DEC_VEC(K))THEN
	            MAX_DEC_VEC(K+2:NV)=MAX_DEC_VEC(K+1:NV-1)
	            MAX_DEC_VEC(K+1)=T1
	            EXIT
	          END IF
	          IF(K .EQ. 1)THEN
	            MAX_DEC_VEC(2:NV)=MAX_DEC_VEC(1:NV-1)
	            MAX_DEC_VEC(1)=T1
	          END IF
	        END DO
	      END IF
	    END IF
	  END DO
	END DO
!
! DO_LEVEL_CHECL allows us to check for oscillating variables, or a few non-coverging levels.
! If the fast majority of the corrections (in this case all bar 2*NV) are less than 10%,
! we set DO_LEVEL_CHK as TRUE. When DO_LEVEL_CHECK is true we will do the following:
!     (a) Ignore the levels with corrections > 10% when computing scale (MAJOR only)
!     (b) Scale the corrections > 10% down by a factor of 3 (-ve relaxation).
!
	DO_LEVEL_CHK=.FALSE.
!
! edited by Devin, 12-Oct-2026: a masked column contributes STEQ=0
! here, so a fully masked step satisfies the "<10% nearly everywhere"
! test and line ~1234 overwrites MAXCH with 0 -- discarding the
! STEQ_MASKED_MAX boost above.  KNOP_mom_cold exited "converged" at
! it6 this way with refresh delta=3.1E10.  Refused demand is not a
! small correction: skip the level check when anything was masked,
! on either the armed or the Lambda mask path.
!
	IF(.NOT. L_ANY_MASK .AND.
	1   MAX_INC_VEC(NV) .GT. -1.0E-01_LDP .AND.
	2   MAX_DEC_VEC(NV) .LT. 1.0E-01_LDP)DO_LEVEL_CHK=.TRUE.
!
! Rather than return the largest change, return the max/min of the NV correction.
!
	IF(DO_LEVEL_CHK)THEN
	  MAXCH=MAX( ABS(MAX_INC_VEC(NV)), MAX_DEC_VEC(NV)/(1.0_LDP-MIN(MAX_DEC_VEC(NV),0.9999_LDP)))
	  MAXCH=100.0_LDP*MAXCH
	END IF
!
! edited by Devin, 18-Oct-2026: [BA_EX_RELFRAC] the exit demand over
! non-exempted rows, on the same new-population scale as MAXCH.  It
! is the raw max (no DO_LEVEL_CHK NV-th-largest reduction), so the
! converged exit is never looser than the reported metric.  The
! masked-demand boost is kept: a refused correction is unsatisfied
! demand on either metric.  Flag off -> BA_EX_MAXCH=MAXCH and the
! exit is identical to stock.
!
	BA_EX_MAXCH=MAXCH
	IF(BA_EX_RELFRAC .GT. 0.0_LDP)THEN
	  T1=-100.0_LDP*EX_INC
	  IF(EX_DEC .LT. 0.99999_LDP)THEN
	    T2EX=100.0_LDP*EX_DEC/(1.0_LDP-EX_DEC)
	  ELSE
	    T2EX=1.0E+07_LDP
	  END IF
	  BA_EX_MAXCH=MIN(1.0E+07_LDP,MAX(T1,T2EX))
	  IF(STEQ_MASKED_MAX .GT. 0.0_LDP)THEN
	    BA_EX_MAXCH=MAX(BA_EX_MAXCH,
	1      MIN(1.0D+07,1.0D+02*STEQ_MASKED_MAX))
	  END IF
	END IF
!
	WRITE(6,*)'Maximm changes as returned by SOLVEBA_V13 is',MAXCH
	WRITE(6,'(A,10ES10.2)')' DEC_VEC: ',MAX_DEC_VEC(1:10)
	WRITE(6,'(A,10ES10.2)')' INC_VEC: ',MAX_INC_VEC(1:10)
!
! edited by Devin, 29-Sep-2026: locate the variable index driving the
! largest corrections.  The depth-only print hides whether the same
! population variable oscillates between BA and LAMBDA iterations;
! J=NT is the temperature row.  Variable-to-species mapping uses the
! per-ion SE block offsets in MODEL_SPEC/STEQ ordering.
!
! edited by Devin, 30-Sep-2026: osc=T flags the case where this BA
! iteration's worst increase sits on the same (variable,depth) as the
! previous BA iteration's worst decrease (or vice versa) -- the
! signature of a Newton overshoot bouncing a single level.
!
! edited by Devin, 06-Oct-2026: the trace now prints as a labelled
! section (----- NAME it ----- header + aligned rows), matching the
! BA_RES_ACC / ITER SUMMARY block style; fields unchanged.
!
	IF(CONV_DIAG_LVL .GE. 2 .AND. (DO_NONMONO_SOB .OR.
	1    DO_NONMONO_KNOP .OR. DO_BA_RES_ACC))THEN
	  OSC_FL=.FALSE.
	  IF(.NOT. LAMBDA_IT)THEN
	    IF(JINC .EQ. PV_JD .AND. IINC .EQ. PV_ID)OSC_FL=.TRUE.
	    IF(JDEC .EQ. PV_JI .AND. IDEC .EQ. PV_II)OSC_FL=.TRUE.
	    PV_JI=JINC
	    PV_II=IINC
	    PV_JD=JDEC
	    PV_ID=IDEC
	  END IF
	  WRITE(LUER,'(/,A,I5,A)')
	1    ' ------------------------- BA13-TRACE it',
	2    MAIN_COUNTER,' ----------------------------------'
	  WRITE(LUER,'(A,L2,4X,L2)')
	1    '    lambda call / BA computed     :',LAMBDA_IT,BA_COMPUTED
	  WRITE(LUER,'(A,ES10.2,2X,I6,I4,2X,A)')
	1    '    largest increase  % @(J,I)    :',ABS(INCREASE),
	2    JINC,IINC,TRIM(EQ_DESC(JINC))
	  WRITE(LUER,'(A,ES10.2,2X,I6,I4,2X,A)')
	1    '    largest decrease  % @(J,I)    :',DECREASE,
	2    JDEC,IDEC,TRIM(EQ_DESC(JDEC))
	  WRITE(LUER,'(A,L2)')
	1    '    level-check                   :',DO_LEVEL_CHK
	  WRITE(LUER,'(A,I3,3X,I3)')
	1    '    rejects (hard / soft)         :',BA_NREJ,BA_NREJ_SOFT
	  WRITE(LUER,'(A,L2)')
	1    '    same-cell oscillation         :',OSC_FL
	  IF(LQUAR_OK)THEN
	    WRITE(LUER,'(A,I5,2X,ES10.2)')
	1    '    quarantined eqns / |STEQ|     :',BA_QUAR_COUNT,
	2    STEQ_QMAX
	  END IF
	END IF
!
!
!**********************************************************************************
!**********************************************************************************
!
! Adjust scale parameter so that the biggest decrease in any variable
! (except T) is a factor of 20, the biggest increase is a factor of CHANGE_LIM.
! T is limited to a maximum change of 20% . Three options are used.
! In case I (LOCAL), the changes are scaled local according to the largest
! single change. In case II (NONE) no scaling is done - the changes for EACH
! variable are limited to the values given above.  Case III (MAJOR) is a
! combination of LOCAL and NONE scaling. We use NONE scaling for species
! whose population is significantly less then the Electron density. In case IV
! (GLOBAL), the biggest change at any depth is used to scale all changes.
! This option is probaly obsolete.
!
!**********************************************************************************
!**********************************************************************************

	IF(CHANGE_LIM .LE. 1.0_LDP)THEN
          WRITE(LUER,'(A,1PE12.4)')' Error in SOLVEBA_V13'
          WRITE(LUER,'(A,1PE12.4)')' Maximum change for normal iteration must be > 1.'
	  STOP
	END IF
	BIG_LIM=(CHANGE_LIM-1.0_LDP)/CHANGE_LIM
!	LIT_LIM=-CHANGE_LIM
	LIT_LIM=1.0_LDP-CHANGE_LIM
	MINSCALE=1.0_LDP
	ALLOCATE(POPS_OLD(NT,ND),STAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Warning - POPS_OLD alloc failed in SOLVEBA_V13;'
	1    ,' MAXCH_APP will fall back to raw MAXCH'
	END IF
	IF(ALLOCATED(POPS_OLD))POPS_OLD(1:NT,1:ND)=POPS(1:NT,1:ND)
	IF(SCALE_OPT(1:5) .EQ. 'LOCAL')THEN
	  DO I=1,ND
	    T1=BIG_LIM			!Prevents division by zero and insures
	    T2=LIT_LIM 			!SCALE=1 if small changes.
	    DO J=1,NT-1
	      T1=MAX(T1,STEQ(J,I))   	!Note + means decrease
	      T2=MIN(T2,STEQ(J,I))      !Note - means increase
	    END DO
	    SCALE=MIN( BIG_LIM/T1, LIT_LIM/T2 )
!
! Limit the change in T to a maximum of MAX_dT_COR, and ensure T > T_MIN.
!
	    T3=MAX( MAX_dT_COR,ABS(STEQ(NT,I)) )
	    SCALE=MIN( MAX_dT_COR/T3,SCALE )
	    IF(STEQ(NT,I) .NE. 0 .AND. POPS(NT,I) .GT. T_MIN .AND.
	1                      POPS(NT,I)*(1.0_LDP-STEQ(NT,I)*SCALE) .LT. T_MIN)THEN
	      SCALE=(1.0_LDP-T_MIN/POPS(NT,I))/STEQ(NT,I)
	    END IF
!
	    DO J=1,NT
	      POPS(J,I)=POPS(J,I)*(1.0_LDP-STEQ(J,I)*SCALE)
	    END DO
	    MINSCALE=MIN(SCALE,MINSCALE)
	  END DO
	  WRITE(LUER,'(A,1PE12.4)')
	1  ' The local minimum value of scale is:',MINSCALE
!
!
	ELSE IF(SCALE_OPT(1:4) .EQ. 'NONE')THEN
	  DO I=1,ND
	    DO J=1,NT-1
	      IF(STEQ(J,I) .GT. BIG_LIM)THEN
	        POPS(J,I)=POPS(J,I)*(1.0_LDP-BIG_LIM)
	        MINSCALE=MIN( BIG_LIM/STEQ(J,I),MINSCALE )
	      ELSE IF(STEQ(J,I) .LT. LIT_LIM)THEN
	        POPS(J,I)=POPS(J,I)*(1.0_LDP-LIT_LIM)
	        MINSCALE=MIN( LIT_LIM/STEQ(J,I),MINSCALE )
	      ELSE
	        POPS(J,I)=POPS(J,I)*(1.0_LDP-STEQ(J,I))
	      END IF
	    END DO
	    IF(POPS(NT,I) .LT. T_MIN)POPS(NT,I)=T_MIN
!
! Limit T to a 20% change if it is a variable.
!
	    SCALE=0.2_LDP/MAX( 0.2_LDP,ABS(STEQ(NT,I)) )
	    POPS(NT,I)=POPS(NT,I)*(1.0_LDP-STEQ(NT,I)*SCALE)
	    MINSCALE=MIN(SCALE,MINSCALE)
	  END DO
!
	  WRITE(LUER,'(A,1PE12.4)')
	1  ' The minimum value of scale for all species is:',MINSCALE
!
!
	ELSE IF(SCALE_OPT(1:5) .EQ. 'MAJOR')THEN
	  CALL FIDDLE_POP_CORRECTIONS_V2(POPS,STEQ,T_MIN,CHANGE_LIM,MAX_dT_COR,
	1        SCALE_OPT,LAMBDA_IT,LU_SUM,NT,ND)
!
!
	ELSE			!Global Scaling !
	  T1=BIG_LIM		!Prevents division by zero and insures
	  T2=LIT_LIM   		!SCALE=1 if small changes.
	  DO I=1,ND
	    DO J=1,NT-1
	      T1=MAX(STEQ(J,I),T1)			!Note + means decrease
	      T2=MIN(STEQ(J,I),T2) 			!Note - means increase
	    END DO
	  END DO
	  SCALE=MIN( BIG_LIM/T1, LIT_LIM/T2)
	  DO I=1,ND
	    IF(STEQ(NT,I) .NE. 0 .AND. POPS(NT,I) .GT. T_MIN .AND.
	1                      POPS(NT,I)*(1.0_LDP-STEQ(NT,I)*SCALE) .LT. T_MIN)THEN
	        SCALE=(1.0_LDP-T_MIN/POPS(NT,I))/STEQ(NT,I)
	    END IF
	  END DO
!
! Limit the change in T to a maximum of 20%.
!
! Upstream bug fix (not flag-gated): the loop used to overwrite T3 at
! each depth, so the 20% bound was enforced only against the last
! depth STEQ(NT,ND).  T3 must accumulate the maximum T correction
! over the whole atmosphere.  Only the SCALE_OPT='GLOBAL' branch is
! affected -- LOCAL, NONE and MAJOR are unchanged.
!
	  T3=0.2_LDP
	  DO I=1,ND
	    T3=MAX( T3,ABS(STEQ(NT,I)) )
	  END DO
	  SCALE=MIN( 0.2_LDP/T3,SCALE )
	  WRITE(LUER,'(A,1PE12.4)')' The value of scale is:',SCALE
!
! Update the population levels (and the temperature) .
!
	  DO I=1,ND
	    DO J=1,NT
	      POPS(J,I)=POPS(J,I)*(1.0_LDP-STEQ(J,I)*SCALE)
	    END DO
	  END DO
	END IF
!
	IF(SET_POPS_D2_EQ_D1 .AND. .NOT. LAMBDA_IT)THEN
	  POPS(:,2)=POPS(:,1)*POP_ATOM(2)/POP_ATOM(1)
	END IF
!
! This adjust populations to their LTE values. This is only done if LTE_MODEL
! is set to TRUE, which is ascertained in SET_POPS_TO_LTE.
!
	IF(LTE_MODEL)CALL SET_POPS_TO_PURE_LTE(POPS,NT,ND)
!
! edited by Devin, 16-Oct-2026: Stage-B trace polish.  When
! BA_TRACE_PHASE is armed by CMFGEN_SUB (majors converged, trace
! deltas still above BA_TR_ABSMAX) the committed update is restricted
! to the BA_DEM-floor-exempted rows; majors, ED and T are restored
! from POPS_OLD.  Trace species are one-way coupled -- they absorb
! and emit but do not feed back on J/ED/T at their populations -- so
! solving them at a fixed major/field state is the lower triangle of
! the block-triangular system and cannot oscillate against the
! majors.  Armed only when a floor is set (otherwise the mask is
! empty and every row is a major); LAMBDA calls are exempt.
! (17-Oct-2026: BA_DEM_RELFRAC removed -- the mask is now PFLOOR-only.)
!
	IF(BA_TRACE_PHASE .AND. .NOT. LAMBDA_IT .AND.
	1   ALLOCATED(POPS_OLD) .AND.
	2   BA_DEM_PFLOOR .GT. 0.0_LDP)THEN
	  DO I=1,ND
	    DO J=1,NT-2
	      IF(POPS_OLD(J,I) .GE. BA_DEM_PFLOOR)POPS(J,I)=POPS_OLD(J,I)
	    END DO
	    POPS(NT-1,I)=POPS_OLD(NT-1,I)
	    POPS(NT,I)=POPS_OLD(NT,I)
	  END DO
	  WRITE(LUER,'(A,I4)')
	1  ' Stage-B trace polish: majors/T/ED held --- iteration ',
	2  MAIN_COUNTER
	END IF
!
! edited by Devin, 29-Sep-2026: report the maximum population change
! actually committed after FIDDLE_POP_CORRECTIONS_V2 / scale limits.
! The raw STEQ-based MAXCH can be 1.0E+7 for a pre-limit decrease that
! leaves the population at a finite floor; using the committed change
! for the KNOP mode decision prevents a raw sentinel from forcing a
! Lambda step after the update was already safely capped.  Computed for
! all paths so the argument is always defined; if the POPS_OLD
! allocation failed, MAXCH_APP falls back to the raw MAXCH.
!
! A population driven non-positive is the worst committed change
! possible and must register as the 1.0E+7 sentinel, not be skipped.
! An increase from a non-positive old value is unbounded -- likewise
! sentinel.  Increases are capped at 1.0E+7 for consistency with the
! raw metric: the committed ratio can exceed it when a post-scale
! overwrite (LTE_MODEL, D2_EQ_D1) or a tiny POPS_OLD inflates the
! ratio beyond what the scale limits imply.
!
! Unlike raw MAXCH this is a strict max over all variables -- the
! DO_LEVEL_CHK NV-th-largest reduction is not applied, so near
! convergence MAXCH_APP can read above MAXCH.  That direction is
! conservative (later Lambda exit, more BA recomputes), which suits
! the KNOP path it was built for.
!
	MAXCH_APP=MAXCH
	IF(ALLOCATED(POPS_OLD))THEN
	  INCREASE=0.0_LDP
	  DECREASE=0.0_LDP
	  DO I=1,ND
	    DO J=1,NT
!
! edited by Devin, 11-Oct-2026: same [BA_DEM_PFLOOR] exemption as the
! demand scan above -- a committed change on a floor-level row is
! equally uninformative for the mode decision.  (17-Oct-2026: the
! [BA_DEM_RELFRAC] POPS_OLD.LT.T2 skip was removed with the key.)
!
	      IF(BA_DEM_PFLOOR .GT. 0.0_LDP .AND. J .LE. NT-2 .AND.
	1     POPS_OLD(J,I) .LT. BA_DEM_PFLOOR)CYCLE
	      IF(POPS_OLD(J,I) .LE. 0.0_LDP)THEN
	        IF(POPS(J,I) .GT. 0.0_LDP)INCREASE=1.0E+07_LDP
	      ELSE IF(POPS(J,I) .LE. 0.0_LDP)THEN
	        DECREASE=1.0E+02_LDP
	      ELSE
	        T1=1.0_LDP-POPS(J,I)/POPS_OLD(J,I)
	        IF(T1 .GT. 0.0_LDP)THEN
	          DECREASE=MAX(DECREASE,100.0_LDP*T1)
	        ELSE IF(T1 .LT. 0.0_LDP)THEN
	          INCREASE=MAX(INCREASE,MIN(1.0E+07_LDP,-100.0_LDP*T1))
	        END IF
	      END IF
	    END DO
	  END DO
	  MAXCH_APP=INCREASE
	  IF(DECREASE .LT. 99.999_LDP)THEN
	    DECREASE=100.0_LDP*DECREASE/(100.0_LDP-DECREASE)
	  ELSE
	    DECREASE=1.0E+07_LDP
	  END IF
	  MAXCH_APP=MIN(1.0E+07_LDP,MAX(MAXCH_APP,DECREASE))
!
! edited by Devin, 15-Oct-2026: [BA_MAXCH_FAIL] (default 1.0D30 =
! off) same-call rejection of a committed non-LAMBDA step.  The
! scale limits bound each correction but never refuse the update,
! and residual acceptance (DO_BA_RES_ACC) can only judge it next
! call -- a run with resacc off has no path that vetoes a committed
! update before the wrecked populations feed the next FG/MOM work.
! When the committed change exceeds BA_MAXCH_FAIL, or any population
! came out non-finite, restore the pre-step POPS_OLD and raise
! BA_STEP_REJECTED: SOLVE_FOR_POPS forces a BA recompute and demotes
! through the existing BA_REJ_THIS path, so no new machinery is
! needed downstream.  No lagged-field rollback is required because
! nothing between the commit above and this point has consumed the
! rejected populations.  MAXCH keeps the 1.5E5 sentinel used by the
! resacc reject so the step stays above the 1.0E5 stay-in-LAMBDA
! guard; .NOT.(X .LE. Y) also catches a NaN MAXCH_APP, which a
! direct .GT. test would miss.  Inert at the default; LAMBDA calls
! are exempt (the Lambda map owns its own convergence behaviour).
!
	  IF(BA_MAXCH_FAIL .LT. 1.0E+29_LDP .AND. .NOT. LAMBDA_IT)THEN
	    T1=0.0_LDP
	    DO I=1,ND
	      DO J=1,NT
	        IF(POPS(J,I) .NE. POPS(J,I) .OR.
	1           ABS(POPS(J,I)) .GT. HUGE(1.0_LDP))T1=1.0_LDP
	      END DO
	    END DO
	    IF(T1 .GT. 0.0_LDP .OR.
	1       .NOT. (MAXCH_APP .LE. BA_MAXCH_FAIL))THEN
	      WRITE(LUER,'(A,I5,A,ES12.4,A,ES12.4)')
	1      ' BA_MAXCH_FAIL: rejecting committed step at it',
	2      MAIN_COUNTER,' -- committed max change',MAXCH_APP,
	3      ' vs limit',BA_MAXCH_FAIL
	      POPS(1:NT,1:ND)=POPS_OLD(1:NT,1:ND)
	      BA_STEP_REJECTED=.TRUE.
	      MAXCH=1.5E+05_LDP
	      MAXCH_APP=MAXCH
	    END IF
	  END IF
!
! edited by Devin, 16-Oct-2026: BA_TR_ABS -- the absolute
! committed change on the floor-exempted (trace) rows, normalized by
! POP_ATOM.  Relative change is meaningless across ~50 decades of
! population: 100% of 1E-55 carries no information, while a few
! percent on 1E-5 matters.  Reported every call so Stage-B polish
! progress (and hidden trace motion before it) is visible; it is
! diagnostic only and does not gate MAXCH.  Runs on the final
! committed state, so a rejected step reads 0.
!
	  ITRC=0
	  JTRC=0
	  DO I=1,ND
	    DO J=1,NT-2
	      IF(POPS_OLD(J,I) .GE. BA_DEM_PFLOOR)CYCLE
	      T1=ABS(POPS(J,I)-POPS_OLD(J,I))/
	1         MAX(POP_ATOM(I),1.0E-30_LDP)
	      IF(T1 .GT. BA_TR_ABS)THEN
	        BA_TR_ABS=T1
	        ITRC=I
	        JTRC=J
	      END IF
	    END DO
	  END DO
	  IF(JTRC .GT. 0 .AND. BA_TRACE_PHASE)THEN
	    WRITE(LUER,9216)ITRC,BA_TR_ABS,TRIM(EQ_DESC(JTRC)),
	1   '          --- iteration ',MAIN_COUNTER
	  END IF
	  DEALLOCATE(POPS_OLD)
	END IF
!
	RETURN
!
	CONTAINS
!
! edited by Devin, 30-Sep-2026: map a global equation index J to a
! readable species:level label.  Each ion's equations start at
! ATM(ID)%EQXzV (set_ba_storage.f); the first NXzV of that block are
! level populations, the rest are ionization/conservation equations
! of that ion/species.  J=NT is T and J=NT-1 is ED for the DIAG
! band layout.  Approximate for J inside a species conservation
! block -- adequate for diagnostics.
!
	FUNCTION EQ_DESC(J) RESULT(NAME)
	  USE MOD_CMFGEN, ONLY: ATM, ION_ID
	  INTEGER J
	  CHARACTER(LEN=60) NAME
	  INTEGER ID,BEST,L,NB
	  IF(J .EQ. NT)THEN
	    NAME='T'
	    RETURN
	  END IF
	  IF(J .EQ. NT-1)THEN
	    NAME='ED'
	    RETURN
	  END IF
	  BEST=0
	  L=NT
	  DO ID=1,NION
	    IF(ATM(ID)%EQXzV .GT. 0 .AND. ATM(ID)%EQXzV .LE. J .AND.
	1       J-ATM(ID)%EQXzV .LT. L)THEN
	      BEST=ID
	      L=J-ATM(ID)%EQXzV
	    END IF
	  END DO
	  IF(BEST .EQ. 0)THEN
	    NAME='?'
	    RETURN
	  END IF
	  L=L+1
	  IF(L .LE. ATM(BEST)%NXzV)THEN
	    IF(ALLOCATED(ATM(BEST)%XzVLEVNAME_F))THEN
	      NAME=TRIM(ION_ID(BEST))//':'//
	1        TRIM(ATM(BEST)%XzVLEVNAME_F(L))
!
! edited by Devin, 02-Oct-2026: XzVLEVNAME_F is allocated for every
! ion but only filled for species whose atomic file supplies level
! names; an unfilled entry is NUL bytes, which TRIM leaves in place
! so the empty-name fallback below never fired and OUTGEN received
! binary garbage ('H2:' followed by NULs).  Map any non-printable
! character to a space so the fallback can see the name as empty.
!
	      DO NB=1,LEN(NAME)
	        IF(ICHAR(NAME(NB:NB)) .LT. 32 .OR.
	1          ICHAR(NAME(NB:NB)) .GT. 126)NAME(NB:NB)=' '
	      END DO
	      IF(LEN_TRIM(NAME) .LE. LEN_TRIM(ION_ID(BEST))+1)
	1        WRITE(NAME,'(A,I6)')TRIM(ION_ID(BEST))//':lev',L
	    ELSE
	      WRITE(NAME,'(A,I6)')TRIM(ION_ID(BEST))//':lev',L
	    END IF
	  ELSE
	    NAME=TRIM(ION_ID(BEST))//':ion-eq'
	  END IF
	  RETURN
	END FUNCTION EQ_DESC
!
! edited by Devin, 07-Oct-2026: worst unquarantined residual
! equation at depth I under [DO_BA_QUAR].  Mirrors the entry column
! scan above: row N_SE is the species conservation residual mapped
! to EQ_SPECIES, and the ionization row at EQXzV folds in STEQ_ADV.
! Rows already parked are skipped so the call parks a different
! equation when a column fails twice.  ED/T rows (NT-1,NT) are
! never candidates.  Returns the global equation index, or 0 when
! no population equation remains free in the column.
!
	FUNCTION BA_QUAR_WORST_EQ(I) RESULT(JRES)
	  INTEGER I
	  INTEGER JRES
	  INTEGER ISPEC,ID,K,JLOC
	  REAL(KIND=LDP) T1
	  JRES=0
	  T1=0.0_LDP
	  DO ISPEC=1,NUM_SPECIES
	    DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	      DO K=1,SE(ID)%N_SE
	        IF(K .EQ. SE(ID)%N_SE)THEN
	          JLOC=EQ_SPECIES(ISPEC)
	        ELSE
	          JLOC=SE(ID)%EQ_IN_BA(K)
	        END IF
	        IF(JLOC .LT. 1 .OR. JLOC .GT. NT-2)CYCLE
	        IF(BA_QUAR_MASK(JLOC,I))CYCLE
	        IF(ABS(SE(ID)%STEQ(K,I)) .GT. T1)THEN
	          T1=ABS(SE(ID)%STEQ(K,I))
	          JRES=JLOC
	        END IF
	      END DO
	      IF(ASSOCIATED(SE(ID)%STEQ_ADV))THEN
	        JLOC=ATM(ID)%EQXzV
	        IF(JLOC .GE. 1 .AND. JLOC .LE. NT-2)THEN
	          IF(.NOT. BA_QUAR_MASK(JLOC,I) .AND.
	1            ABS(SE(ID)%STEQ_ADV(I)) .GT. T1)THEN
	            T1=ABS(SE(ID)%STEQ_ADV(I))
	            JRES=JLOC
	          END IF
	        END IF
	      END IF
	    END DO
	  END DO
	  RETURN
	END FUNCTION BA_QUAR_WORST_EQ
!
	END
