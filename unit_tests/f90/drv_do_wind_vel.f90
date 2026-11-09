! Driver for DO_WIND_VEL_V1 (new_main/subs/do_wind_vel_v1.f).
!
! Verifies the 05-Oct-2026 guard: DO_V_PARAMS is rejected under
! DO_NONMONO_SOB / DO_NONMONO_KNOP. The STOP only fires when the
! optional DO_V_PARAMS file is present, so the harness controls file
! presence by running the binary in a per-case scratch directory.
!
! Cases:
!   sob_file    DO_NONMONO_SOB=T + file present -> expect STOP text
!   knop_file   DO_NONMONO_KNOP=T + file present -> expect STOP text
!   sob_nofile  DO_NONMONO_SOB=T + no file -> clean return (no-file
!               default path must be unchanged)
!   off_file    flags off + file present with FST_IT > counter ->
!               clean early return, DONE_V_REVISION stays F
!
      PROGRAM DRV_DO_WIND_VEL
      USE SET_KIND_MODULE
      USE CONTROL_VARIABLE_MOD
      IMPLICIT NONE
      INTEGER, PARAMETER :: ND=4
      REAL(KIND=LDP) R(ND),V(ND),SIGMA(ND),CL(ND),MSTAR
      LOGICAL DONE
      CHARACTER(LEN=32) CASE
      INTEGER NARG
      EXTERNAL DO_WIND_VEL_V1
!
      NARG=COMMAND_ARGUMENT_COUNT()
      IF(NARG .NE. 1)STOP 'usage: drv_do_wind_vel <case>'
      CALL GET_COMMAND_ARGUMENT(1,CASE)
!
      DO_NONMONO_SOB=.FALSE.; DO_NONMONO_KNOP=.FALSE.
      SELECT CASE(TRIM(CASE))
      CASE('sob_file')
        DO_NONMONO_SOB=.TRUE.
      CASE('knop_file')
        DO_NONMONO_KNOP=.TRUE.
      CASE('sob_nofile')
        DO_NONMONO_SOB=.TRUE.
      CASE('off_file')
      CASE DEFAULT
        STOP 'bad case'
      END SELECT
!
      R=(/ 4.0_LDP,3.0_LDP,2.0_LDP,1.0_LDP /)
      V=(/ 2000._LDP,1500._LDP,1000._LDP,500._LDP /)
      SIGMA=0.5_LDP; CL=1.0_LDP; MSTAR=10.0_LDP
      DONE=.FALSE.
!
      CALL DO_WIND_VEL_V1(R,V,SIGMA,CL,MSTAR,DONE,0,ND)
      WRITE(*,'(A,1X,L1)')'RETURNED',DONE
!
      END
