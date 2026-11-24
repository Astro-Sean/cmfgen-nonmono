! lam_handoff_mod.f90 -- Lambda-to-full handoff predicate and
! stability-counter step for the non-monotonic iteration cycle.
!
! Two pieces of iteration-control logic lived inline in
! CMFGEN_SUB (the one-time DO_LAM_AUTO release) and SOLVE_FOR_POPS
! (the RD_LAMBDA stability counter that feeds it). They are
! extracted here unchanged in meaning so both production sites and
! the unit driver exercise the same code (unit_tests/f90/
! drv_lam_handoff.f90).
!
! The one behavioural change versus the inline version is the
! field-readiness leg: LAM_FIELD_GATE combines the lagged-field
! maturity state (GS_RAW_IMMATURE plus the N_MATURE/GS_REL_STABLE
! window) with the hold/freeze flags into a single
! "field not ready" verdict. On the KNOP route it reduces to
! GS_POPS_HELD .OR. GS_FREEZE -- the terms the old predicate used --
! because GS_FREEZE there already carries the raw-immature and
! not-yet-mature state. On the SOB lagged-field routes
! (DO_NONMONO_SOB + USE_LAM_ES=F, DO_SOB_FG_KNOP) it adds the
! coverage those routes lacked: a quiet population demand against
! a churning lagged field can no longer release the model.
!
! Edited by Devin (handoff hardening, Nov-2026)
!
      MODULE LAM_HANDOFF_MOD
      USE SET_KIND_MODULE
      IMPLICIT NONE
!
      CONTAINS
!
! Field-readiness verdict for one SOLVE_FOR_POPS call. Without a
! lagged field (stock, plain SOB+LAM_ES=T formal-in-sweep) the
! lagged-field terms cannot be judged, so the verdict is just the
! held/frozen flags -- both already forced .FALSE. outside
! DO_NONMONO_KNOP, so stock returns .FALSE. always.
!
      LOGICAL FUNCTION LAM_FIELD_GATE(HAS_LAGGED,HELD,FREEZE,&
     &     RAW_IMMATURE,N_MATURE,REL_STABLE)
      LOGICAL, INTENT(IN) :: HAS_LAGGED,HELD,FREEZE,RAW_IMMATURE
      INTEGER, INTENT(IN) :: N_MATURE,REL_STABLE
      LAM_FIELD_GATE=HELD .OR. FREEZE .OR.&
     &     (HAS_LAGGED .AND.&
     &      (RAW_IMMATURE .OR. N_MATURE .LT. REL_STABLE))
      RETURN
      END FUNCTION LAM_FIELD_GATE
!
! One step of the RD_LAMBDA stability counter, called once per
! RD_LAMBDA LAMBDA call in SOLVE_FOR_POPS. Semantics preserved
! from the inline block:
!   * reset: forced LAMBDA, luminosity-blocked, or demand at/above
!     VAL_EXIT -- plus the new FIELD_NRDY leg, which under KNOP is
!     already covered by the force a freeze raises, and on the SOB
!     lagged-field routes is the missing protection.
!   * earn: only a call that ran free (not frozen, field ready)
!     with demand below VAL_EXIT counts toward N_STABLE.
!   * pin: a resacc-accepted probe (RDLAM_ACC) pins the count at
!     N_STABLE, but only while the field reports ready -- an
!     accepted probe is contraction evidence for the free map,
!     not for the lagged field.
!
      SUBROUTINE LAM_STABLE_STEP(FORCE_LAM,LUM_BLK,FROZEN,FIELD_NRDY,&
     &     MAXCH_CTRL,VAL_EXIT,N_STABLE,RDLAM_ACC,&
     &     CNT_STABLE,STABLE_PROVEN)
      LOGICAL, INTENT(IN) :: FORCE_LAM,LUM_BLK,FROZEN,FIELD_NRDY
      LOGICAL, INTENT(IN) :: RDLAM_ACC
      REAL(KIND=LDP), INTENT(IN) :: MAXCH_CTRL,VAL_EXIT
      INTEGER, INTENT(IN) :: N_STABLE
      INTEGER, INTENT(INOUT) :: CNT_STABLE
      LOGICAL, INTENT(INOUT) :: STABLE_PROVEN
      IF(FORCE_LAM .OR. LUM_BLK .OR. FIELD_NRDY .OR.&
     &   MAXCH_CTRL .GE. VAL_EXIT)THEN
        CNT_STABLE=0
        STABLE_PROVEN=.FALSE.
      ELSE IF(.NOT. FROZEN)THEN
        CNT_STABLE=CNT_STABLE+1
        IF(CNT_STABLE .GE. N_STABLE)STABLE_PROVEN=.TRUE.
      END IF
      IF(RDLAM_ACC .AND. .NOT. FIELD_NRDY .AND. .NOT. FORCE_LAM&
     &   .AND. .NOT. LUM_BLK)THEN
        CNT_STABLE=N_STABLE
        STABLE_PROVEN=.TRUE.
      END IF
      RETURN
      END SUBROUTINE LAM_STABLE_STEP
!
! The one-time DO_LAMBDA_AUTO handoff predicate, called from
! CMFGEN_SUB. NONMONO_ACTIVE is (DO_NONMONO_SOB .OR.
! DO_NONMONO_KNOP) at the call site; FIELD_NRDY is the
! LAM_FIELD_GATE verdict for the just-completed call and is only
! consulted on non-monotonic paths.
!
      LOGICAL FUNCTION LAM_HANDOFF_OK(DO_LAM_AUTO,RD_LAM,MAXCH_CTRL,&
     &     JFNK_UNCONV,NONMONO_ACTIVE,STABLE_OK,FIELD_NRDY)
      LOGICAL, INTENT(IN) :: DO_LAM_AUTO,RD_LAM,JFNK_UNCONV
      LOGICAL, INTENT(IN) :: NONMONO_ACTIVE,STABLE_OK,FIELD_NRDY
      REAL(KIND=LDP), INTENT(IN) :: MAXCH_CTRL
      LAM_HANDOFF_OK=DO_LAM_AUTO .AND. RD_LAM .AND.&
     &     MAXCH_CTRL .LT. 50.0_LDP .AND. .NOT. JFNK_UNCONV
      IF(LAM_HANDOFF_OK .AND. NONMONO_ACTIVE)THEN
        LAM_HANDOFF_OK=STABLE_OK .AND. .NOT. FIELD_NRDY
      END IF
      RETURN
      END FUNCTION LAM_HANDOFF_OK
!
! Nine-gate converged-exit replica for the SOLVE_FOR_POPS ITER_SUM
! report. CMFGEN_SUB evaluates the real CONV_EXIT_MET inline at both
! 9999 exit splits; this evaluator takes the same resolved state,
! recomputes the lambda-state leg itself, and returns the verdict,
! the pass count and the failure labels from one source. Routing the
! replica through one evaluator keeps it from drifting a conjunct or
! a label away from the real gate -- the missing ninth conjunct
! (.NOT. LAM_FIELD_NOT_READY, the audit F-32 finding) is exactly the
! drift class this removes. Label order matches the historical
! replica; the ninth label is appended last.
!
! edited 24-Nov-2026 [KN_LAM_PARTIAL_HANDOFF, default F]: the ninth
! conjunct reads the strict FIELD_NCONV (LAM_FIELD_NOT_CONVERGED --
! gate verdict plus the rf_not_allconv fold-in), matching the real
! exit. FIELD_NRDY stays an input so the FAILS label can still
! separate a weak-gate hold ('lam-field-not-ready') from a
! fold-only hold ('lam-field-not-converged'); FIELD_NCONV implies
! FIELD_NRDY under the flag-off wiring, so the default label text
! is unchanged.
!
      SUBROUTINE CONV_EXIT_EVAL(RD_LAM,DO_LAM_AUTO,LAST_LAM,&
     &     MAIN_COUNTER,JFNK_UNCONV,EX_MAXCH,EPS,DEM_OK,&
     &     HYDRO_REV,POPS_HELD,RESACC_PEND,LAST_NG,FIELD_NRDY,&
     &     FIELD_NCONV,LAM_OK,ALL_MET,NPASS,FAILS)
      LOGICAL, INTENT(IN) :: RD_LAM,DO_LAM_AUTO,JFNK_UNCONV,DEM_OK
      LOGICAL, INTENT(IN) :: HYDRO_REV,POPS_HELD,RESACC_PEND,FIELD_NRDY
      LOGICAL, INTENT(IN) :: FIELD_NCONV
      INTEGER, INTENT(IN) :: LAST_LAM,MAIN_COUNTER,LAST_NG
      REAL(KIND=LDP), INTENT(IN) :: EX_MAXCH,EPS
      LOGICAL, INTENT(OUT) :: LAM_OK,ALL_MET
      INTEGER, INTENT(OUT) :: NPASS
      CHARACTER(LEN=*), INTENT(OUT) :: FAILS
      LAM_OK=(RD_LAM .AND. .NOT. DO_LAM_AUTO) .OR.&
     &     (LAST_LAM .NE. MAIN_COUNTER)
      ALL_MET=LAM_OK .AND. .NOT. JFNK_UNCONV .AND.&
     &     (EX_MAXCH .LT. EPS) .AND. DEM_OK .AND.&
     &     .NOT. HYDRO_REV .AND. .NOT. POPS_HELD .AND.&
     &     .NOT. RESACC_PEND .AND. (MAIN_COUNTER .NE. LAST_NG) .AND.&
     &     .NOT. FIELD_NCONV
      NPASS=0
      IF(LAM_OK)NPASS=NPASS+1
      IF(.NOT. JFNK_UNCONV)NPASS=NPASS+1
      IF(EX_MAXCH .LT. EPS)NPASS=NPASS+1
      IF(DEM_OK)NPASS=NPASS+1
      IF(.NOT. HYDRO_REV)NPASS=NPASS+1
      IF(.NOT. POPS_HELD)NPASS=NPASS+1
      IF(.NOT. RESACC_PEND)NPASS=NPASS+1
      IF(MAIN_COUNTER .NE. LAST_NG)NPASS=NPASS+1
      IF(.NOT. FIELD_NCONV)NPASS=NPASS+1
      FAILS=' all pass'
      IF(.NOT. ALL_MET)THEN
        FAILS=' '
        IF(.NOT. LAM_OK)FAILS=TRIM(FAILS)//' lambda-state'
        IF(.NOT. (EX_MAXCH .LT. EPS))&
     &       FAILS=TRIM(FAILS)//' ex_maxch>=eps'
        IF(JFNK_UNCONV)FAILS=TRIM(FAILS)//' jfnk-unconverged'
        IF(.NOT. DEM_OK)FAILS=TRIM(FAILS)//' exempted-demand'
        IF(HYDRO_REV)FAILS=TRIM(FAILS)//' hydro-revision'
        IF(POPS_HELD)FAILS=TRIM(FAILS)//' pops-held'
        IF(RESACC_PEND)FAILS=TRIM(FAILS)//' resacc-pending'
        IF(MAIN_COUNTER .EQ. LAST_NG)&
     &       FAILS=TRIM(FAILS)//' on-ng-iteration'
        IF(FIELD_NCONV)THEN
          IF(FIELD_NRDY)THEN
            FAILS=TRIM(FAILS)//' lam-field-not-ready'
          ELSE
            FAILS=TRIM(FAILS)//' lam-field-not-converged'
          END IF
        END IF
      END IF
      RETURN
      END SUBROUTINE CONV_EXIT_EVAL
!
! Worst-margin classification for the ITER_SUM 'conv exit'/'conv
! marg' rows (round-2 phase-2 patch B). The exit gate releases
! exempted demand on three legs -- small, saturated
! (>= BA_SAT_RAW) and trace-handled (Stage-B need or release) --
! but the old reporter priced every channel against its threshold
! regardless, so a saturated BA_DEM_EXM always won the "worst
! margin" slot while the actual binding gate sat elsewhere. The
! classifier derives each channel's release state from the same
! predicates the gate applies, returns the worst margin over the
! binding channels only, and reports the largest released demand
! and the Stage-B absolute demand as separate fields. Stage-B is
! its own class: it is not a release of the exempted-demand gate
! and it never joins the SAT channel.
!
      SUBROUTINE CONV_MARGIN_EVAL(EPS,EX_MAXCH,EXM,EXR,SAT_RAW,&
     &     TRC_NEED,TRC_REL,TR_ABS,TR_ABSMAX,TR_MAXSWP,&
     &     MGN_ACT,CH_ACT,MGN_REL,CH_REL,WHY_REL,MGN_TRB)
      REAL(KIND=LDP), INTENT(IN) :: EPS,EX_MAXCH,EXM,EXR,SAT_RAW
      LOGICAL, INTENT(IN) :: TRC_NEED,TRC_REL
      REAL(KIND=LDP), INTENT(IN) :: TR_ABS,TR_ABSMAX
      INTEGER, INTENT(IN) :: TR_MAXSWP
      REAL(KIND=LDP), INTENT(OUT) :: MGN_ACT,MGN_REL,MGN_TRB
      CHARACTER(LEN=*), INTENT(OUT) :: CH_ACT,CH_REL,WHY_REL
      REAL(KIND=LDP) SMALL
      LOGICAL EXM_REL,EXR_REL
      SMALL=MAX(10.0_LDP*EPS,1.0_LDP)
      EXM_REL=(EXM .LT. SMALL) .OR. (EXM .GE. SAT_RAW) .OR.&
     &     TRC_NEED .OR. TRC_REL
      EXR_REL=(EXR .LT. SMALL) .OR. (EXR .GE. SAT_RAW)
!
! Binding set: ex_maxch always holds the exit; an exempted-demand
! channel holds it only while none of its release legs is true.
!
      MGN_ACT=-99.0_LDP
      CH_ACT='none'
      IF(EPS .GT. 0.0_LDP .AND. EX_MAXCH .GT. 0.0_LDP)THEN
        MGN_ACT=LOG10(EX_MAXCH/EPS)
        CH_ACT='ex_maxch'
      END IF
      IF(.NOT. EXM_REL .AND. EXM .GT. 0.0_LDP .AND.&
     &     LOG10(EXM/SMALL) .GT. MGN_ACT)THEN
        MGN_ACT=LOG10(EXM/SMALL)
        CH_ACT='exm'
      END IF
      IF(.NOT. EXR_REL .AND. EXR .GT. 0.0_LDP .AND.&
     &     LOG10(EXR/SMALL) .GT. MGN_ACT)THEN
        MGN_ACT=LOG10(EXR/SMALL)
        CH_ACT='exr'
      END IF
!
! Released set: the largest demand the gate let go. The reason
! string preserves which release leg fired -- SAT first, since a
! demand above the raw saturation bound dominates the report.
!
      MGN_REL=-99.0_LDP
      CH_REL='none'
      WHY_REL='-'
      IF(EXM_REL .AND. EXM .GT. 0.0_LDP .AND.&
     &     LOG10(EXM/SMALL) .GT. MGN_REL)THEN
        MGN_REL=LOG10(EXM/SMALL)
        CH_REL='exm'
        IF(EXM .GE. SAT_RAW)THEN
          WHY_REL='SAT'
        ELSE IF(TRC_NEED)THEN
          WHY_REL='trc-need'
        ELSE IF(TRC_REL)THEN
          WHY_REL='trc-rel'
        ELSE
          WHY_REL='small'
        END IF
      END IF
      IF(EXR_REL .AND. EXR .GT. 0.0_LDP .AND.&
     &     LOG10(EXR/SMALL) .GT. MGN_REL)THEN
        MGN_REL=LOG10(EXR/SMALL)
        CH_REL='exr'
        IF(EXR .GE. SAT_RAW)THEN
          WHY_REL='SAT'
        ELSE
          WHY_REL='small'
        END IF
      END IF
!
! Stage-B absolute demand vs its armed threshold; -99 when the
! Stage-B sweep is switched off (BA_TR_MAXSWP=0) or the absolute
! leg is disabled (BA_TR_ABSMAX=0).
!
      MGN_TRB=-99.0_LDP
      IF(TR_MAXSWP .GT. 0 .AND. TR_ABSMAX .GT. 0.0_LDP .AND.&
     &     TR_ABS .GT. 0.0_LDP)MGN_TRB=LOG10(TR_ABS/TR_ABSMAX)
      RETURN
      END SUBROUTINE CONV_MARGIN_EVAL
!
      END MODULE LAM_HANDOFF_MOD
