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
      END MODULE LAM_HANDOFF_MOD
