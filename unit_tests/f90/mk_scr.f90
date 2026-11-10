! Fixture writer for test_extrap_scr.py.
!
! Writes one complete SCRTEMP/POINT1/POINT2 checkpoint via the
! production SCR_RITE_V2 so extrap_scr.exe can be exercised on real
! checkpoint bytes.  The velocity pattern is chosen by argument so
! each monotonicity guard can be triggered on demand.
!
! Usage: mk_scr CASE   (inside the model dir being fabricated)
!   mono : positive, strictly decreasing V (valid input)
!   flat : V repeats one value (ordering violation)
!   neg  : V goes non-positive at the last depth
!
      PROGRAM MK_SCR
      USE SET_KIND_MODULE
      IMPLICIT NONE
      INTEGER, PARAMETER :: ND=700, NT=4, LU=23
      REAL(KIND=LDP), ALLOCATABLE :: R(:),V(:),SIGMA(:),POPS(:)
      CHARACTER(LEN=32) CASE
      INTEGER I,IREC,NITSF,NUM_TIMES,LST_NG
      LOGICAL WRITFAIL
      CALL GETARG(1,CASE)
      ALLOCATE(R(ND),V(ND),SIGMA(ND),POPS(NT*ND))
      DO I=1,ND
        R(I)=100.0_LDP-0.1_LDP*I
        V(I)=50.0_LDP-0.05_LDP*I
        SIGMA(I)=-0.5_LDP+0.01_LDP*I
      END DO
      DO I=1,NT*ND
        POPS(I)=I*1.7_LDP+0.3_LDP
      END DO
      IF(CASE .EQ. 'flat')V(3)=V(2)
      IF(CASE .EQ. 'neg')V(ND)=-1.0_LDP
      NITSF=17; NUM_TIMES=3; LST_NG=8
      IREC=0
      CALL SCR_RITE_V2(R,V,SIGMA,POPS,&
     &     IREC,NITSF,NUM_TIMES,LST_NG,.TRUE.,NT,ND,LU,WRITFAIL)
      WRITE(*,'(A,L1,1X,I3)')'WF ',WRITFAIL,IREC
      END
