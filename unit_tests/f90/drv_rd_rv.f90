! Driver for RD_RV_FILE_V2 NONMONO_V / VINF consistency checks.
!
!   drv_rd_rv <nonmono T|F> <vinf> <rmax> <rp>
!
! Reads RVSIG_COL from the current working directory.  Prints
! RESULT_OK plus the scaled R/V/SIGMA vectors on success.  A guard
! STOP in the production code leaves no RESULT_OK.
!
      PROGRAM DRV_RD_RV
      IMPLICIT NONE
!
      INTEGER, PARAMETER :: LDP=KIND(1.0D0)
      INTEGER, PARAMETER :: ND=20
!
      REAL(KIND=LDP) R(ND), V(ND), SIGMA(ND), RMAX, RP, VINF
      CHARACTER(LEN=8) :: ARG
      CHARACTER(LEN=132) :: SARG
      CHARACTER*80 OPTIONS(4)
      LOGICAL NONMONO_V
      INTEGER I, LUIN, N_OPT
!
      CALL GET_COMMAND_ARGUMENT(1,ARG)
      NONMONO_V=(ARG(1:1) .EQ. 'T')
      CALL GET_COMMAND_ARGUMENT(2,SARG)
      READ(SARG,*)VINF
      CALL GET_COMMAND_ARGUMENT(3,SARG)
      READ(SARG,*)RMAX
      CALL GET_COMMAND_ARGUMENT(4,SARG)
      READ(SARG,*)RP
!
      OPTIONS(1)='RVSIG_COL'
      N_OPT=1
      LUIN=10
!
      CALL RD_RV_FILE_WRAP(R,V,SIGMA,RMAX,RP,VINF,LUIN,ND,OPTIONS, &
                           N_OPT,NONMONO_V)
!
      WRITE(*,'(a)')'RESULT_OK'
      DO I=1,ND
        WRITE(*,'(i6,3(1x,es24.16e3))') I,R(I),V(I),SIGMA(I)
      END DO
!
      END
