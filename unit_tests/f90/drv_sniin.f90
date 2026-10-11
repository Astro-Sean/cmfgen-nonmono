! Driver for RV_SN_MODEL_SNIIN (VEL_LAW=12 IIn hydro grid).
!
!   drv_sniin <rdinr T|F> <rmax> <rp> <vcore> <beta1> <vfill>
!
! RDINR=T: reads the RDINR file in cwd (beta-law velocity).
! RDINR=F: reads input_hydro.dat in cwd (hydro tau_es grid path).
!          The production code calls LIN_INTERP with the output V
!          array as the target coordinate, so what V contains on entry
!          decides what happens:
!          vfill=Z: V=0 (not monotonic) -> guard aborts.
!          vfill=I: V filled strictly increasing inside V_HYDRO's
!                   range -> guard passes, SIGMA loop runs on this V.
!          vfill=D: V filled strictly decreasing -> direction mismatch
!                   aborts.
!
! Prints RESULT_OK plus R,V,SIGMA on success.
!
      PROGRAM DRV_SNIIN
      IMPLICIT NONE
!
      INTEGER, PARAMETER :: LDP=KIND(1.0D0)
      INTEGER, PARAMETER :: ND=20
!
      REAL(KIND=LDP) R(ND), V(ND), SIGMA(ND)
      REAL(KIND=LDP) RMAX, RP, VCORE, BETA1
      LOGICAL RDINR
      CHARACTER(LEN=8) :: ARG
      CHARACTER(LEN=132) :: SARG
      INTEGER I, LU
!
      CALL GET_COMMAND_ARGUMENT(1,ARG)
      RDINR=(ARG(1:1) .EQ. 'T')
      CALL GET_COMMAND_ARGUMENT(2,SARG)
      READ(SARG,*)RMAX
      CALL GET_COMMAND_ARGUMENT(3,SARG)
      READ(SARG,*)RP
      CALL GET_COMMAND_ARGUMENT(4,SARG)
      READ(SARG,*)VCORE
      CALL GET_COMMAND_ARGUMENT(5,SARG)
      READ(SARG,*)BETA1
      CALL GET_COMMAND_ARGUMENT(6,ARG)
!
      R=0.0_LDP
      V=0.0_LDP
      SIGMA=0.0_LDP
      IF(ARG(1:1) .EQ. 'I')THEN
        DO I=1,ND
          V(I)=310.0_LDP+(8900.0_LDP-310.0_LDP)*REAL(I-1,LDP)/(ND-1)
        END DO
      ELSE IF(ARG(1:1) .EQ. 'D')THEN
        DO I=1,ND
          V(I)=8900.0_LDP-(8900.0_LDP-310.0_LDP)*REAL(I-1,LDP)/(ND-1)
        END DO
      END IF
      LU=10
      OPEN(127,FILE='fort127.log',STATUS='UNKNOWN')
!
      CALL RV_SN_MODEL_SNIIN(R,V,SIGMA,RMAX,RP,VCORE,BETA1, &
                             RDINR,LU,ND)
!
      WRITE(*,'(a)')'RESULT_OK'
      DO I=1,ND
        WRITE(*,'(i6,3(1x,es24.16e3))') I,R(I),V(I),SIGMA(I)
      END DO
!
      END
