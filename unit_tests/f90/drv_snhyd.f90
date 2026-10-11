! Driver for RV_SN_MODEL_SNHYD (VEL_LAW=13 hydro grid reader).
!
!   drv_snhyd <sn_age_days> <nd> <tau_opt>
!
! Reads SN_HYDRO_DATA from the cwd, calls the routine, and prints
! RESULT_OK plus R,V,SIGMA on success.
!
      PROGRAM DRV_SNHYD
      IMPLICIT NONE
!
      INTEGER, PARAMETER :: LDP=KIND(1.0D0)
      INTEGER, PARAMETER :: MX=80
!
      REAL(KIND=LDP) R(MX), V(MX), SIGMA(MX)
      REAL(KIND=LDP) RMAX, RCORE, SN_AGE_DAYS
      CHARACTER(LEN=132) :: SARG
      INTEGER I, ND, LU, TAU_OPT
!
      CALL GET_COMMAND_ARGUMENT(1,SARG)
      READ(SARG,*)SN_AGE_DAYS
      CALL GET_COMMAND_ARGUMENT(2,SARG)
      READ(SARG,*)ND
      CALL GET_COMMAND_ARGUMENT(3,SARG)
      READ(SARG,*)TAU_OPT
!
      R=0.0_LDP
      V=0.0_LDP
      SIGMA=0.0_LDP
      RMAX=0.0_LDP
      RCORE=0.0_LDP
      LU=11
!
      CALL RV_SN_MODEL_SNHYD(R,V,SIGMA,RMAX,RCORE,SN_AGE_DAYS,ND,LU, &
                             TAU_OPT)
!
      WRITE(*,'(A)')'RESULT_OK'
      WRITE(*,'(2ES24.16E3)')RMAX,RCORE
      DO I=1,ND
        WRITE(*,'(I4,1X,3(1X,ES24.16E3))')I,R(I),V(I),SIGMA(I)
      END DO
      END
