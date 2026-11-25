! Standalone stage-E2 driver for the production matrix-free block
! GMRES module. It reads a system written by drv_mom_blk and writes
! the solver vector for comparison with the python dense reference.
!
! Usage: drv_mom_blk_solve <tag>
!
      PROGRAM DRV_MOM_BLK_SOLVE
      USE SET_KIND_MODULE
      USE MOM_BLK_GMRES_MOD
      IMPLICIT NONE
!
      INTEGER ND,NCF
      PARAMETER (ND=5,NCF=12)
!
      REAL(KIND=LDP) TA(ND,NCF),TB(ND,NCF),TC(ND,NCF)
      REAL(KIND=LDP) X0(ND,NCF),CUP(ND,NCF),XOLD(ND,NCF)
      REAL(KIND=LDP) X(ND,NCF)
      INTEGER KUP(ND,NCF)
      LOGICAL FRESH(NCF)
      CHARACTER(LEN=64) TAG,LINE
      INTEGER I,K,II,KK,IOS,NITER,NREST,INFO
      REAL(KIND=LDP) RESID,RESRAW,TOL
!
      CALL GETARG(1,TAG)
      IF(LEN_TRIM(TAG) .EQ. 0)STOP 'DRV_MOM_BLK_SOLVE: missing tag'
      OPEN(UNIT=302,FILE=TRIM(TAG)//'_sys.txt',STATUS='OLD',
     1  IOSTAT=IOS)
      IF(IOS .NE. 0)STOP 'DRV_MOM_BLK_SOLVE: open failed'
      READ(302,'(A)')LINE
      READ(302,'(A)')LINE
      DO K=1,NCF
        DO I=1,ND
          READ(302,*)II,KK,TA(I,K),TB(I,K),TC(I,K),X0(I,K),
     1      CUP(I,K),KUP(I,K),FRESH(K),XOLD(I,K)
        END DO
      END DO
      CLOSE(302)
!
      TOL=1.0E-14_LDP
      CALL MOM_BLK_GMRES_SOLVE(ND,NCF,TA,TB,TC,X0,CUP,KUP,
     1  FRESH,XOLD,X,40,TOL,6,NITER,NREST,RESID,INFO,RESRAW)
!
      OPEN(UNIT=303,FILE=TRIM(TAG)//'_xgm.txt',STATUS='REPLACE',
     1  IOSTAT=IOS)
      IF(IOS .NE. 0)STOP 'DRV_MOM_BLK_SOLVE: write failed'
      DO K=1,NCF
        DO I=1,ND
          WRITE(303,'(2I5,ES25.16)')I,K,X(I,K)
        END DO
      END DO
      CLOSE(303)
      WRITE(*,'(A,I3,A,I5,A,I3,A,ES12.5,A,ES12.5)')
     1  'GMRES INFO=',INFO,' NITER=',NITER,' NREST=',NREST,
     2  ' RESID=',RESID,' RESRAW=',RESRAW
      IF(INFO .EQ. 0)WRITE(*,'(A)')'RESULT_OK'
      END
