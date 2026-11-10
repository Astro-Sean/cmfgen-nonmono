! Driver for ADJUST_ATM_R_GRID (misc/adjust_atm_r_grid.f).
!
! Alteration under test: the connection search and log
! interpolation use V as a coordinate, so OLD_V must be positive
! and strictly decreasing.  Upstream silently located an
! arbitrary radius branch on a repeated velocity.
!
! Usage: drv CASE
!   flat : OLD_V repeats a value -> labelled STOP
!   neg  : OLD_V goes negative -> labelled STOP
!   good : strictly decreasing positive OLD_V -> guard passes
!
      PROGRAM DRV_ADJ_GRID
      USE SET_KIND_MODULE
      IMPLICIT NONE
      INTEGER, PARAMETER :: NS=6, ND=10
      REAL(KIND=LDP) R(ND),OLD_R(NS),OLD_V(NS),OLD_TAU(NS)
      REAL(KIND=LDP) V_CON,V_RAT_MAX,IB_RAT,OB_RAT,DTAU2_ON_DTAU1
      INTEGER N_IB_INS,N_OB_INS,I
      CHARACTER(LEN=16) CASE
      CALL GETARG(1,CASE)
      DO I=1,NS
        OLD_R(I)=110.0_LDP-10.0_LDP*I
        OLD_V(I)=330.0_LDP-50.0_LDP*I
        OLD_TAU(I)=10.0_LDP**(I-3)
      END DO
      IF(CASE .EQ. 'flat')OLD_V(3)=OLD_V(2)
      IF(CASE .EQ. 'neg')OLD_V(4)=-1.0_LDP
      V_CON=50.0_LDP
      V_RAT_MAX=1.4_LDP
      IB_RAT=1.3_LDP
      OB_RAT=1.3_LDP
      DTAU2_ON_DTAU1=1.0_LDP
      N_IB_INS=1
      N_OB_INS=1
      CALL ADJUST_ATM_R_GRID(R,OLD_R,OLD_V,OLD_TAU,&
     &     V_CON,V_RAT_MAX,&
     &     IB_RAT,OB_RAT,DTAU2_ON_DTAU1,N_IB_INS,N_OB_INS,ND,NS)
      WRITE(*,'(A)')'DONE_'
      DO I=1,ND
        WRITE(*,'(I3,1X,ES16.8)')I,R(I)
      END DO
      END
