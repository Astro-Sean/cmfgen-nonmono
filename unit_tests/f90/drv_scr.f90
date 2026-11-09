! Driver for SCR_RITE_V2 / SCR_READ_V2 (subs/scr_read_v2.f).
!
! Alterations under test:
!   - WRITFAIL output flag: any IOSTAT error on the SCRTEMP or
!     POINT writes now propagates to the caller instead of leaving
!     a torn POINT/SCRTEMP pair that a requeue would silently resume.
!   - the trailing DUMP_FG_INEXT_GS / DUMP_MOM_JNEXT_GS calls that
!     checkpoint the Knop GS lagged-radiation state (they return
!     immediately here because the nonmono arrays are unallocated).
!   - INQUIRE/CLOSE so failed opens do not leak the unit.
!
! ND=700, NT=4 forces NUM_RV_RECS>1 and NUMRECS>1 on this platform
! (N_PER_REC=2047), so both multi-record read/write paths run.
!
! Usage: drv CASE  (inside a scratch cwd; SCRTEMP/POINT* live there)
!   rt     : write once (WRITE_RVSIG=T), read back last record,
!            bitwise compare of every array plus flag checks
!   rt2    : two writes; IREC_RD=0 returns the second, IREC_RD=1 the
!            first -- record-addressing across iterations
!   newmod : read with no POINT/SCRTEMP -> NEWMOD=T path
!   wrfail : write into an unwritable dir -> WRITFAIL=T
!   norv   : WRITE_RVSIG=F -> RVSIG_WRITTEN=F, R/V/SIGMA untouched
!
      PROGRAM DRV_SCR
      USE SET_KIND_MODULE
      IMPLICIT NONE
      INTEGER, PARAMETER :: ND=700, NT=4, LU=23
      REAL(KIND=LDP), ALLOCATABLE :: R(:),V(:),SIGMA(:),POPS(:)
      REAL(KIND=LDP), ALLOCATABLE :: R0(:),V0(:),S0(:),P0(:)
      REAL(KIND=LDP), ALLOCATABLE :: R1(:),V1(:),S1(:),P1(:)
      CHARACTER(LEN=32) CASE
      INTEGER I,IREC,NITSF,NUM_TIMES,LST_NG,IREC_RD
      LOGICAL WRITE_RVSIG,WRITFAIL,RVSIG_WRITTEN,NEWMOD,OK
      CALL GETARG(1,CASE)
      ALLOCATE(R(ND),V(ND),SIGMA(ND),POPS(NT*ND))
      ALLOCATE(R0(ND),V0(ND),S0(ND),P0(NT*ND))
      ALLOCATE(R1(ND),V1(ND),S1(ND),P1(NT*ND))
      NITSF=17; NUM_TIMES=3; LST_NG=8
!
      IF(CASE .EQ. 'wrfail')THEN
        IREC=0
        CALL SCR_RITE_V2(R,V,SIGMA,POPS,&
     &       IREC,NITSF,NUM_TIMES,LST_NG,.TRUE.,NT,ND,LU,WRITFAIL)
        WRITE(*,'(A,L1)')'WF ',WRITFAIL
        STOP
      END IF
!
      IF(CASE .EQ. 'newmod')THEN
        IREC_RD=0
        CALL SCR_READ_V2(R,V,SIGMA,POPS,&
     &       IREC_RD,NITSF,NUM_TIMES,LST_NG,RVSIG_WRITTEN,&
     &       NT,ND,LU,NEWMOD)
        WRITE(*,'(A,L1,1X,I3)')'NM ',NEWMOD,IREC_RD
        STOP
      END IF
!
      CALL FILL(R,V,SIGMA,POPS,ND,NT,0)
      WRITE_RVSIG=(CASE .NE. 'norv')
      IREC=0
      CALL SCR_RITE_V2(R,V,SIGMA,POPS,&
     &     IREC,NITSF,NUM_TIMES,LST_NG,WRITE_RVSIG,NT,ND,LU,WRITFAIL)
      WRITE(*,'(A,L1,1X,I3)')'WF1 ',WRITFAIL,IREC
      R0=R; V0=V; S0=SIGMA; P0=POPS
!
      IF(CASE .EQ. 'rt2')THEN
        CALL FILL(R,V,SIGMA,POPS,ND,NT,7)
        R1=R; V1=V; S1=SIGMA; P1=POPS
        NITSF=18
        CALL SCR_RITE_V2(R,V,SIGMA,POPS,&
     &       IREC,NITSF,NUM_TIMES,LST_NG,.TRUE.,NT,ND,LU,WRITFAIL)
        WRITE(*,'(A,L1,1X,I3)')'WF2 ',WRITFAIL,IREC
        NITSF=17
      END IF
!
! Read back: IREC_RD=0 selects the last written iteration.
!
      R=0.0_LDP; V=0.0_LDP; SIGMA=0.0_LDP; POPS=0.0_LDP
      IREC_RD=0; NITSF=-1; NUM_TIMES=-1; LST_NG=-1
      CALL SCR_READ_V2(R,V,SIGMA,POPS,&
     &     IREC_RD,NITSF,NUM_TIMES,LST_NG,RVSIG_WRITTEN,&
     &     NT,ND,LU,NEWMOD)
      WRITE(*,'(A,L1,1X,I3,1X,L1)')'RD ',NEWMOD,IREC_RD,RVSIG_WRITTEN
      WRITE(*,'(A,3(I4,1X))')'META ',NITSF,NUM_TIMES,LST_NG
      IF(CASE .EQ. 'rt2')THEN
        OK=ALL(R .EQ. R1) .AND. ALL(V .EQ. V1) .AND.&
     &     ALL(SIGMA .EQ. S1) .AND. ALL(POPS .EQ. P1)
      ELSE
        OK=ALL(R .EQ. R0) .AND. ALL(V .EQ. V0) .AND.&
     &     ALL(SIGMA .EQ. S0) .AND. ALL(POPS .EQ. P0)
      END IF
      WRITE(*,'(A,L1)')'EQ1 ',OK
!
      IF(CASE .EQ. 'rt2')THEN
        R=0.0_LDP; V=0.0_LDP; SIGMA=0.0_LDP; POPS=0.0_LDP
        IREC_RD=1
        CALL SCR_READ_V2(R,V,SIGMA,POPS,&
     &       IREC_RD,NITSF,NUM_TIMES,LST_NG,RVSIG_WRITTEN,&
     &       NT,ND,LU,NEWMOD)
        OK=ALL(R .EQ. R0) .AND. ALL(V .EQ. V0) .AND.&
     &     ALL(SIGMA .EQ. S0) .AND. ALL(POPS .EQ. P0)
        WRITE(*,'(A,L1,1X,I3)')'EQ2 ',OK,IREC_RD
      END IF
      END
!
      SUBROUTINE FILL(R,V,SIGMA,POPS,ND,NT,K)
      USE SET_KIND_MODULE
      IMPLICIT NONE
      INTEGER ND,NT,K
      REAL(KIND=LDP) R(ND),V(ND),SIGMA(ND),POPS(NT*ND)
      INTEGER I
      DO I=1,ND
        R(I)=100.0_LDP-0.1_LDP*I+K
        V(I)=50.0_LDP-0.05_LDP*I+0.5_LDP*K
        SIGMA(I)=-0.5_LDP+0.01_LDP*I-0.1_LDP*K
      END DO
      DO I=1,NT*ND
        POPS(I)=I*1.7_LDP+0.3_LDP+K
      END DO
      END
