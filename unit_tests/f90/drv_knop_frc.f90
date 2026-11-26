! Driver for the WP-6 forcing-term refresh tolerance
! (KN_FORCING_TOL; audit KNOP_REMEDIATION_WP6_2026-11-26.md).
!
! Exercises the production machinery exactly as CMFGEN_SUB does:
! SET_FG_FORCING_TOL installs flag/knobs, SET_FG_FRC_RES hands in
! the last population demand, and KN_FRC_UPDATE runs one
! Eisenstat & Walker 1996 Choice-2 step per refresh:
!
!   eta_k = GAM * (F_k/F_{k-1})^ALP, safeguarded below by
!   GAM*eta_{k-1}^ALP when that floor exceeds 0.1, capped at 0.9;
!   KN_RAY_TOL_EFF = clamp(eta*TMAX/0.9, TMIN, min(TMAX,CAP)).
!
! Cases:
!   seed  first armed call -- no prior residual, so the update
!         only seeds FPRV and leaves TOL_EFF at KN_RAY_TOL
!   seq   falling demand 1.0 -> 0.5 -> 0.25 -> 0.1: eta and the
!         tolerance tighten monotonically with the contraction
!   rise  flat demand 1.0 -> 1.0 -> 1.0: eta hits the 0.9 cap and
!         the tolerance lands on the I2 cap (5E-3), not TMAX
!   sflo  safeguard floor: with eta_prev=0.6 a collapse to
!         rat=0.3 raises eta to GAM*0.36 instead of dropping
!   tmin  fast contraction rat=0.05: the TMIN floor binds
!   noup  update with no pending demand -- TOL_EFF stays at base
!   getr  GET_FG_FRC export matches the last update
!
      PROGRAM DRV_KNOP_FRC
      USE SET_KIND_MODULE
      USE FG_J_CMF_NONMONO_MOD_V1
      IMPLICIT NONE
!
      INTEGER NCASE,I
      REAL(KIND=LDP) FSEQ(4),ETA_G,RAT_G,TOL_G
      INTEGER SG_G
      CHARACTER(LEN=32) CNAME
!
      CALL GETARG(1,CNAME)
      NCASE=0
      IF(CNAME .EQ. 'seed')NCASE=1
      IF(CNAME .EQ. 'seq') NCASE=2
      IF(CNAME .EQ. 'rise')NCASE=3
      IF(CNAME .EQ. 'sflo')NCASE=4
      IF(CNAME .EQ. 'tmin')NCASE=5
      IF(CNAME .EQ. 'noup')NCASE=6
      IF(CNAME .EQ. 'getr')NCASE=7
      IF(NCASE .EQ. 0)THEN
        WRITE(*,'(A)')'RESULT_FAIL unknown_case'
        STOP
      END IF
!
! Production wiring: flag + knobs through the VADAT setter,
! I2 cap = 0.1*GS_MATURE_TOL with the default 5E-2 gate.
!
      KN_RAY_TOL=2.0E-3_LDP
      CALL SET_FG_FORCING_TOL(.TRUE.,0.9_LDP,2.0_LDP, &
      &1.0E-4_LDP,1.0E-2_LDP,5.0E-3_LDP)
!
      IF(NCASE .EQ. 1)THEN
        CALL SET_FG_FRC_RES(1.0_LDP)
        CALL KN_FRC_UPDATE()
        WRITE(*,'(A,ES16.9)')'FRC_ETA=',KN_FRC_ETA
        WRITE(*,'(A,ES16.9)')'FRC_FPRV=',KN_FRC_FPRV
        WRITE(*,'(A,ES16.9)')'FRC_TOL=',KN_RAY_TOL_EFF
        WRITE(*,'(A)')'RESULT_OK'
        STOP
      END IF
!
      IF(NCASE .EQ. 6)THEN
! No pending demand: the update must not move the tolerance.
        KN_RAY_TOL_EFF=0.0_LDP
        CALL KN_FRC_UPDATE()
        WRITE(*,'(A,ES16.9)')'FRC_TOL=',KN_RAY_TOL_EFF
        WRITE(*,'(A,ES16.9)')'FRC_ETA=',KN_FRC_ETA
        WRITE(*,'(A)')'RESULT_OK'
        STOP
      END IF
!
      IF(NCASE .EQ. 4)THEN
! Seed eta_prev by hand -- the safeguard reads the stored eta.
        KN_FRC_ETA=0.6_LDP
        KN_FRC_FPRV=1.0_LDP
        CALL SET_FG_FRC_RES(0.3_LDP)
        CALL KN_FRC_UPDATE()
        WRITE(*,'(A,ES16.9)')'FRC_ETA=',KN_FRC_ETA
        WRITE(*,'(A,ES16.9)')'FRC_RAT=',KN_FRC_RAT
        WRITE(*,'(A,I3)')'FRC_SG=',KN_FRC_SG
        WRITE(*,'(A,ES16.9)')'FRC_TOL=',KN_RAY_TOL_EFF
        WRITE(*,'(A)')'RESULT_OK'
        STOP
      END IF
!
! Sequence cases: feed the demand list through the production
! setter -> update path, one step per "iteration".
!
      IF(NCASE .EQ. 2)THEN
        FSEQ=(/1.0_LDP,0.5_LDP,0.25_LDP,0.1_LDP/)
      ELSE IF(NCASE .EQ. 3)THEN
        FSEQ=(/1.0_LDP,1.0_LDP,1.0_LDP,1.0_LDP/)
      ELSE IF(NCASE .EQ. 5)THEN
        FSEQ=(/1.0_LDP,0.05_LDP,0.01_LDP,0.005_LDP/)
      ELSE
        FSEQ=(/1.0_LDP,0.4_LDP,0.16_LDP,0.064_LDP/)
      END IF
      DO I=1,4
        CALL SET_FG_FRC_RES(FSEQ(I))
        CALL KN_FRC_UPDATE()
        WRITE(*,'(A,I3,A,ES16.9,A,ES16.9,A,I3,A,ES16.9)') &
      &'FRC_STEP ',I,' ETA=',KN_FRC_ETA,' RAT=',KN_FRC_RAT, &
      &' SG=',KN_FRC_SG,' TOL=',KN_RAY_TOL_EFF
      END DO
!
      IF(NCASE .EQ. 7)THEN
        CALL GET_FG_FRC(ETA_G,RAT_G,SG_G,TOL_G)
        WRITE(*,'(A,ES16.9)')'GET_ETA=',ETA_G
        WRITE(*,'(A,ES16.9)')'GET_RAT=',RAT_G
        WRITE(*,'(A,I3)')'GET_SG=',SG_G
        WRITE(*,'(A,ES16.9)')'GET_TOL=',TOL_G
      END IF
!
      WRITE(*,'(A)')'RESULT_OK'
      STOP
!
      END
