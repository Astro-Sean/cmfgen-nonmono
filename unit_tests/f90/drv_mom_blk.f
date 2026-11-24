! Driver for manufactured wavelength-depth block systems.  The python
! battery assembles the dense reference operator from the written
! coefficients and verifies test GMRES implementations against it.
! No production solver code is exercised by this stage-E1 driver.
!
! Rows obey
!   TA*x(i-1,k) + TB*x(i,k) + TC*x(i+1,k)
!     - CUP*x(i,kup) = X0.
! CUP is the positive upwind RHS coefficient used by the production
! capture, so it carries a minus sign after moving the term to the LHS.
!
! Modes extend drv_mom_nu's graph cases with nonzero depth coupling:
!   exp    expanding chain
!   cmp    compressing chain
!   ring2  adjacent expanding/compressing ring
!   ring4  stride ring
!   marg   near-marginal ring
!   mix    chains feeding two rings
!   depth  strongly nonuniform depth coefficients
!   fresh  inactive frequency columns frozen on the old field
!
! Usage: drv_mom_blk <mode> <tag>
!
      PROGRAM DRV_MOM_BLK
      USE SET_KIND_MODULE
      IMPLICIT NONE
!
      INTEGER ND,NCF
      PARAMETER (ND=5,NCF=12)
!
      REAL(KIND=LDP) TA(ND,NCF),TB(ND,NCF),TC(ND,NCF)
      REAL(KIND=LDP) X0(ND,NCF),CUP(ND,NCF),XOLD(ND,NCF)
      INTEGER KUP(ND,NCF)
      LOGICAL FRESH(NCF)
      CHARACTER(LEN=32) MODE,OUTTAG
      INTEGER I,K,IOS
!
      CALL GETARG(1,MODE)
      CALL GETARG(2,OUTTAG)
      IF(LEN_TRIM(OUTTAG) .EQ. 0)OUTTAG='blk'
!
      CALL BLK_SEED_BASE(ND,NCF,TA,TB,TC,X0,CUP,KUP,XOLD,FRESH)
!
      SELECT CASE(MODE)
      CASE('exp')
        DO K=2,NCF
          DO I=1,ND-1
            CUP(I,K)=0.80_LDP
            KUP(I,K)=K-1
          END DO
        END DO
      CASE('cmp')
        DO K=1,NCF-1
          DO I=1,ND-1
            CUP(I,K)=0.80_LDP
            KUP(I,K)=K+1
          END DO
        END DO
      CASE('ring2')
        DO I=1,ND-1
          CUP(I,2)=0.90_LDP;  KUP(I,2)=3
          CUP(I,3)=0.65_LDP;  KUP(I,3)=2
        END DO
      CASE('ring4')
        DO I=1,ND-1
          CUP(I,1)=0.80_LDP;  KUP(I,1)=4
          CUP(I,4)=0.55_LDP;  KUP(I,4)=3
          CUP(I,3)=0.55_LDP;  KUP(I,3)=2
          CUP(I,2)=0.55_LDP;  KUP(I,2)=1
        END DO
      CASE('marg')
! Near-marginal ring: unlike MOM_NUSOLV, the block solve has no freeze
! band.  It must solve the nonsingular coupled system directly.
        DO I=1,ND-1
          CUP(I,2)=4.00_LDP;  KUP(I,2)=3
          CUP(I,3)=3.90_LDP;  KUP(I,3)=2
        END DO
      CASE('mix')
        DO I=1,ND-1
          CUP(I,2)=0.45_LDP;  KUP(I,2)=1
          CUP(I,3)=0.75_LDP;  KUP(I,3)=4
          CUP(I,4)=0.55_LDP;  KUP(I,4)=3
          CUP(I,5)=0.60_LDP;  KUP(I,5)=6
          CUP(I,6)=0.75_LDP;  KUP(I,6)=7
          CUP(I,7)=0.55_LDP;  KUP(I,7)=6
          CUP(I,9)=0.50_LDP;  KUP(I,9)=11
          CUP(I,11)=0.45_LDP; KUP(I,11)=10
          CUP(I,10)=0.40_LDP; KUP(I,10)=9
        END DO
      CASE('depth')
        DO K=2,NCF
          DO I=1,ND-1
            CUP(I,K)=0.65_LDP
            KUP(I,K)=K-1
          END DO
        END DO
        DO K=1,NCF
          DO I=2,ND-1
            TA(I,K)=TA(I,K)*(1.0_LDP+0.05_LDP*K)
            TC(I,K)=TC(I,K)*(1.0_LDP+0.03_LDP*I)
            IF(MOD(I+K,5) .EQ. 0)TC(I,K)=0.25_LDP*ABS(TC(I,K))
          END DO
        END DO
      CASE('fresh')
        DO K=2,NCF
          DO I=1,ND-1
            CUP(I,K)=0.75_LDP
            KUP(I,K)=K-1
          END DO
        END DO
        FRESH(7:NCF)=.FALSE.
        DO K=7,NCF
          DO I=1,ND
            XOLD(I,K)=3.0_LDP+0.01_LDP*I
          END DO
        END DO
      CASE DEFAULT
        WRITE(*,'(A)')'RESULT unknown_mode=1'
        STOP 1
      END SELECT
!
! Fully implicit upwind weighting puts the diagonal absorption into TB
! (the Baron-Hauschildt chi-hat piece).  The marginal case models the
! production rows whose effective absorption is descaled: their diagonal
! carries only a fraction of the upwind weight, which is exactly where
! the lagged frequency map loses contraction.
!
      DO K=1,NCF
        DO I=1,ND
          IF(CUP(I,K) .NE. 0.0_LDP)THEN
            IF(MODE .EQ. 'marg')THEN
              TB(I,K)=TB(I,K)+0.25_LDP*CUP(I,K)
            ELSE
              TB(I,K)=TB(I,K)+CUP(I,K)
            END IF
          END IF
        END DO
      END DO
!
      OPEN(UNIT=301,FILE=TRIM(OUTTAG)//'_sys.txt',STATUS='REPLACE',
     1     IOSTAT=IOS)
      IF(IOS .NE. 0)STOP 'DRV_MOM_BLK: open failed'
      WRITE(301,'(A,2I8)')'# ND NCF: ',ND,NCF
      WRITE(301,'(A)')'# I K TA TB TC X0 CUP KUP FRESH XOLD'
      DO K=1,NCF
        DO I=1,ND
          WRITE(301,'(2I5,5ES25.16,I8,L3,ES25.16)')I,K,TA(I,K),
     1      TB(I,K),TC(I,K),X0(I,K),CUP(I,K),KUP(I,K),FRESH(K),
     2      XOLD(I,K)
        END DO
      END DO
      CLOSE(301)
      WRITE(*,'(A)')'RESULT_OK'
      END
!
!-----------------------------------------------------------------------
!
      SUBROUTINE BLK_SEED_BASE(ND,NCF,TA,TB,TC,X0,CUP,KUP,XOLD,FRESH)
      USE SET_KIND_MODULE
      IMPLICIT NONE
      INTEGER ND,NCF
      REAL(KIND=LDP) TA(ND,NCF),TB(ND,NCF),TC(ND,NCF)
      REAL(KIND=LDP) X0(ND,NCF),CUP(ND,NCF),XOLD(ND,NCF)
      INTEGER KUP(ND,NCF)
      LOGICAL FRESH(NCF)
      INTEGER I,K
!
      DO K=1,NCF
        FRESH(K)=.TRUE.
        DO I=1,ND
          TA(I,K)=0.0_LDP
          TB(I,K)=2.40_LDP+0.03_LDP*K+0.02_LDP*I
          TC(I,K)=0.0_LDP
          X0(I,K)=0.50_LDP+0.010_LDP*K+0.004_LDP*I
          CUP(I,K)=0.0_LDP
          KUP(I,K)=0
          XOLD(I,K)=1.0_LDP+0.02_LDP*I+0.01_LDP*K
        END DO
        DO I=2,ND
          TA(I,K)=-0.18_LDP-0.010_LDP*I
        END DO
        DO I=1,ND-1
          TC(I,K)=-0.22_LDP-0.005_LDP*K
        END DO
      END DO
      RETURN
      END
