! Driver for GET_IBOUND (plane/get_ibound.f).
!
! The 05-Oct-2026 guards validate the INCIDENT_INTENSITY file at load
! (NCF>=2, finite positive strictly decreasing NU_GRID) and bound the
! saved index walk (CUR_NU_INDX, ML clamp, non-finite FREQ -> 0).
! The routine reads the file once on the first call, so each case is
! one binary run in a scratch dir that holds the case's file.
!
!   good     valid file; FREQ walk 35,25,15,5,NaN -- interpolation
!            and out-of-range / non-finite -> IBOUND=0 paths.
!            INIT=T on every call: CUR_NU_INDX is a non-SAVE local
!            upstream, so INIT=F cursor reuse is undefined behaviour
!            (the 05-Oct-2026 bound only clamps above; a garbage
!            negative index still OOBs below). The per-call walk,
!            interpolation and range branches all still exercise.
!   ncf1     NCF=1 file -> STOP 'need at least 2'
!   nonpos   NU_GRID contains -5 -> STOP 'non-finite or non-positive'
!   notdec   NU_GRID not strictly decreasing -> STOP
!   badmu    MU outside ANG_GRID range -> STOP 'Invalid incident MU'
!   reorder  second call with increasing FREQ -> STOP 'out of order'
!
      PROGRAM DRV_GET_IBOUND
      USE, INTRINSIC :: IEEE_ARITHMETIC
      USE SET_KIND_MODULE
      IMPLICIT NONE
      INTEGER, PARAMETER :: NP=2
      REAL(KIND=LDP) MU(NP),IBOUND(NP),FREQ
      CHARACTER(LEN=32) CASE
      INTEGER NARG,I
      EXTERNAL GET_IBOUND
!
      NARG=COMMAND_ARGUMENT_COUNT()
      IF(NARG .NE. 1)STOP 'usage: drv_get_ibound <case>'
      CALL GET_COMMAND_ARGUMENT(1,CASE)
      MU=(/ 0.5_LDP,0.75_LDP /)
!
      SELECT CASE(TRIM(CASE))
      CASE('good')
        FREQ=35.0_LDP
        CALL GET_IBOUND(IBOUND,MU,NP,FREQ,.TRUE.)
        WRITE(*,'(A,2ES20.12)')'IB',IBOUND
        FREQ=25.0_LDP
        CALL GET_IBOUND(IBOUND,MU,NP,FREQ,.TRUE.)
        WRITE(*,'(A,2ES20.12)')'IB',IBOUND
        FREQ=15.0_LDP
        CALL GET_IBOUND(IBOUND,MU,NP,FREQ,.TRUE.)
        WRITE(*,'(A,2ES20.12)')'IB',IBOUND
        FREQ=5.0_LDP
        CALL GET_IBOUND(IBOUND,MU,NP,FREQ,.TRUE.)
        WRITE(*,'(A,2ES20.12)')'IB',IBOUND
        FREQ=IEEE_VALUE(0.0_LDP,IEEE_QUIET_NAN)
        CALL GET_IBOUND(IBOUND,MU,NP,FREQ,.TRUE.)
        WRITE(*,'(A,2ES20.12)')'IB',IBOUND
      CASE('badmu')
        MU=(/ 5.0_LDP,5.5_LDP /)
        CALL GET_IBOUND(IBOUND,MU,NP,25.0_LDP,.TRUE.)
        WRITE(*,'(A,2ES20.12)')'IB',IBOUND
      CASE('reorder')
        CALL GET_IBOUND(IBOUND,MU,NP,15.0_LDP,.TRUE.)
        WRITE(*,'(A,2ES20.12)')'IB',IBOUND
        CALL GET_IBOUND(IBOUND,MU,NP,25.0_LDP,.FALSE.)
        WRITE(*,'(A,2ES20.12)')'IB',IBOUND
      CASE('ncf1','nonpos','notdec')
        CALL GET_IBOUND(IBOUND,MU,NP,25.0_LDP,.TRUE.)
        WRITE(*,'(A,2ES20.12)')'IB',IBOUND
      CASE DEFAULT
        STOP 'bad case'
      END SELECT
      WRITE(*,'(A)')'DONE'
      END
