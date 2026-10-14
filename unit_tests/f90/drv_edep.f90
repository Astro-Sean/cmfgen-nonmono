! Driver for GET_EDEP_SHOCK_POWER.
!
!   drv_edep <mode> <vloc> <dvloc> <psp> <scl> [rloc] [drloc]
!
!   mode: F = both nonmono flags F, S = DO_NONMONO_SOB,
!         K = DO_NONMONO_KNOP, B = both T
!   scl:  T/F sets SCL_PWR_BY_FCL
!   rloc/drloc: optional; rloc>0 selects radius-space deposition
!
! Reads edep_input.dat from the cwd: line 1 = ND, then ND rows of
! R(1e10cm) V(km/s) DENSITY CLUMP_FAC.  On success prints RESULT_OK
! plus the SHOCK_POWER array.
!
      PROGRAM DRV_EDEP
      USE SET_KIND_MODULE
      USE MOD_CMFGEN, ONLY : CLUMP_FAC
      USE SHOCK_POWER_MOD
      USE CONTROL_VARIABLE_MOD, ONLY : PRESCRIBED_SHOCK_POWER, &
     &    VLOC_SHOCK_POWER, DVLOC_SHOCK_POWER, SCL_PWR_BY_FCL, &
     &    RLOC_SHOCK_POWER, DRLOC_SHOCK_POWER, &
     &    DO_NONMONO_SOB, DO_NONMONO_KNOP
      IMPLICIT NONE
!
      INTEGER, PARAMETER :: MX=64
      REAL(KIND=LDP) R(MX), V(MX), DENS(MX)
      CHARACTER(LEN=8) :: MODE, SCL
      CHARACTER(LEN=132) :: SARG
      INTEGER I, ND
!
      CALL GET_COMMAND_ARGUMENT(1,MODE)
      CALL GET_COMMAND_ARGUMENT(2,SARG)
      READ(SARG,*)VLOC_SHOCK_POWER
      CALL GET_COMMAND_ARGUMENT(3,SARG)
      READ(SARG,*)DVLOC_SHOCK_POWER
      CALL GET_COMMAND_ARGUMENT(4,SARG)
      READ(SARG,*)PRESCRIBED_SHOCK_POWER
      CALL GET_COMMAND_ARGUMENT(5,SCL)
!
      RLOC_SHOCK_POWER=0.0_LDP
      DRLOC_SHOCK_POWER=0.0_LDP
      IF(COMMAND_ARGUMENT_COUNT() .GE. 6)THEN
        CALL GET_COMMAND_ARGUMENT(6,SARG)
        READ(SARG,*)RLOC_SHOCK_POWER
      END IF
      IF(COMMAND_ARGUMENT_COUNT() .GE. 7)THEN
        CALL GET_COMMAND_ARGUMENT(7,SARG)
        READ(SARG,*)DRLOC_SHOCK_POWER
      END IF
!
      DO_NONMONO_SOB=(MODE .EQ. 'S' .OR. MODE .EQ. 'B')
      DO_NONMONO_KNOP=(MODE .EQ. 'K' .OR. MODE .EQ. 'B')
      SCL_PWR_BY_FCL=(SCL(1:1) .EQ. 'T')
!
      OPEN(21,FILE='edep_input.dat',STATUS='OLD')
      READ(21,*)ND
      ALLOCATE(CLUMP_FAC(ND))
      DO I=1,ND
        READ(21,*)R(I),V(I),DENS(I),CLUMP_FAC(I)
      END DO
      CLOSE(21)
!
      CALL GET_EDEP_SHOCK_POWER(R,V,DENS,ND)
!
      WRITE(*,'(A)')'RESULT_OK'
      DO I=1,ND
        WRITE(*,'(I4,1X,2(1X,ES24.16E3))')I,V(I),SHOCK_POWER(I)
      END DO
      END
