! Reference-tree driver for GET_EDEP_SHOCK_POWER.  Identical to
! drv_edep.f90 minus the DO_NONMONO_* variables that only exist in
! the new control_variable_mod.  The mode argument is accepted for
! command-line compatibility but ignored.
!
      PROGRAM DRV_EDEP
      USE SET_KIND_MODULE
      USE MOD_CMFGEN, ONLY : CLUMP_FAC
      USE SHOCK_POWER_MOD
      USE CONTROL_VARIABLE_MOD, ONLY : PRESCRIBED_SHOCK_POWER, &
     &    VLOC_SHOCK_POWER, DVLOC_SHOCK_POWER, SCL_PWR_BY_FCL
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
