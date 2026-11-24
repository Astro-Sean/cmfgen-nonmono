! Matrix-free coupled wavelength-depth GMRES solver.
!
! [MOM_BLKGMRES, default F once wired] The captured row equations are
!     TA*x(i-1,k) + TB*x(i,k) + TC*x(i+1,k)
!       - CUP*x(i,kup) = X0.
! CUP is stored as the positive RHS coefficient of the upwind term, so
! moving that term to the implicit left-hand side gives the minus sign.
! The depth-tridiagonal part of TB already contains the implicit
! diagonal absorption; this solver adds no second PSI contribution.
! Inactive frequency columns are pinned to their incoming iterate by
! identity rows. The solve is matrix-free: only rowwise actions and the
! per-frequency DGTSV preconditioner are formed.
!
      MODULE MOM_BLK_GMRES_MOD
      USE SET_KIND_MODULE
      IMPLICIT NONE
!
      PRIVATE
      PUBLIC :: MOM_BLK_GMRES_SOLVE
!
      CONTAINS
!
!-----------------------------------------------------------------------
!
      SUBROUTINE MOM_BLK_GMRES_SOLVE(ND,NCF,TA,TB,TC,X0,CUP,KUP,
     1  FRESH,XIN,X,RESTART,TOL,MAXRST,NITER,NREST,RESID,INFO)
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: ND,NCF,RESTART,MAXRST
      INTEGER, INTENT(IN) :: KUP(ND,NCF)
      REAL(KIND=LDP), INTENT(IN) :: TA(ND,NCF),TB(ND,NCF)
      REAL(KIND=LDP), INTENT(IN) :: TC(ND,NCF),X0(ND,NCF)
      REAL(KIND=LDP), INTENT(IN) :: CUP(ND,NCF),XIN(ND,NCF)
      LOGICAL, INTENT(IN) :: FRESH(NCF)
      REAL(KIND=LDP), INTENT(OUT) :: X(ND,NCF)
      REAL(KIND=LDP), INTENT(IN) :: TOL
      INTEGER, INTENT(OUT) :: NITER,NREST,INFO
      REAL(KIND=LDP), INTENT(OUT) :: RESID
!
      INTEGER N,M,I,J,JUSED,K,L,IOS,IPS
      REAL(KIND=LDP) BETA,BNRM,RNRM,HN,T1,T2,RR
      REAL(KIND=LDP), ALLOCATABLE :: XV(:),B(:),Z(:),W(:),Y(:)
      REAL(KIND=LDP), ALLOCATABLE :: V(:,:),H(:,:),CS(:),SN(:),G(:)
!
      INFO=2
      NITER=0
      NREST=0
      RESID=HUGE(RESID)
      X=XIN
      IF(ND .LT. 1 .OR. NCF .LT. 1 .OR. RESTART .LT. 1 .OR.
     1   MAXRST .LT. 1 .OR. TOL .LE. 0.0_LDP)RETURN
!
      N=ND*NCF
      M=RESTART
      ALLOCATE(XV(N),B(N),Z(N),W(N),Y(M),V(N,M+1),H(M+1,M),
     1  CS(M),SN(M),G(M+1),STAT=IOS)
      IF(IOS .NE. 0)THEN
        INFO=5
        RETURN
      END IF
!
! Identity rows pin inactive columns. This keeps the vector layout
! dense while preserving the sweep's active-set semantics.
!
      DO K=1,NCF
        DO I=1,ND
          L=I+(K-1)*ND
          XV(L)=XIN(I,K)
          IF(FRESH(K))THEN
            B(L)=X0(I,K)
          ELSE
            B(L)=XIN(I,K)
          END IF
        END DO
      END DO
!
      CALL BLK_PRECOND_VEC(ND,NCF,TA,TB,TC,FRESH,B,Z,IPS)
      IF(IPS .NE. 0)THEN
        INFO=3
        GOTO 900
      END IF
      BNRM=SQRT(DOT_PRODUCT(Z,Z))
      IF(BNRM .EQ. 0.0_LDP)BNRM=1.0_LDP
!
      DO NREST=0,MAXRST-1
        CALL BLK_MATVEC(ND,NCF,TA,TB,TC,CUP,KUP,FRESH,XV,W)
        Z=B-W
        CALL BLK_PRECOND_VEC(ND,NCF,TA,TB,TC,FRESH,Z,W,IPS)
        IF(IPS .NE. 0)THEN
          INFO=3
          GOTO 900
        END IF
        BETA=SQRT(DOT_PRODUCT(W,W))
        RNRM=BETA/BNRM
        RESID=RNRM
        IF(RNRM .LE. TOL)THEN
          INFO=0
          GOTO 800
        END IF
!
        V=0.0_LDP
        H=0.0_LDP
        CS=0.0_LDP
        SN=0.0_LDP
        G=0.0_LDP
        V(:,1)=W/BETA
        G(1)=BETA
        JUSED=0
!
        DO J=1,M
          CALL BLK_MATVEC(ND,NCF,TA,TB,TC,CUP,KUP,FRESH,V(:,J),Z)
          CALL BLK_PRECOND_VEC(ND,NCF,TA,TB,TC,FRESH,Z,W,IPS)
          IF(IPS .NE. 0)THEN
            INFO=3
            GOTO 900
          END IF
!
! Two classical Gram-Schmidt passes cost little here and avoid loss of
! orthogonality when the block preconditioner clusters the spectrum.
!
          DO IPS=1,2
            DO I=1,J
              T1=DOT_PRODUCT(V(:,I),W)
              H(I,J)=H(I,J)+T1
              W=W-T1*V(:,I)
            END DO
          END DO
          HN=SQRT(DOT_PRODUCT(W,W))
          H(J+1,J)=HN
          IF(HN .GT. 0.0_LDP)V(:,J+1)=W/HN
!
! Apply the earlier Givens rotations to the new Hessenberg column.
!
          DO I=1,J-1
            T1=CS(I)*H(I,J)+SN(I)*H(I+1,J)
            H(I+1,J)=-SN(I)*H(I,J)+CS(I)*H(I+1,J)
            H(I,J)=T1
          END DO
!
! Form the rotation that zeroes H(J+1,J).
!
          IF(H(J+1,J) .EQ. 0.0_LDP)THEN
            CS(J)=1.0_LDP
            SN(J)=0.0_LDP
          ELSE IF(ABS(H(J+1,J)) .GT. ABS(H(J,J)))THEN
            T1=H(J,J)/H(J+1,J)
            SN(J)=1.0_LDP/SQRT(1.0_LDP+T1*T1)
            CS(J)=T1*SN(J)
          ELSE
            T1=H(J+1,J)/H(J,J)
            CS(J)=1.0_LDP/SQRT(1.0_LDP+T1*T1)
            SN(J)=T1*CS(J)
          END IF
          T1=CS(J)*H(J,J)+SN(J)*H(J+1,J)
          H(J+1,J)=-SN(J)*H(J,J)+CS(J)*H(J+1,J)
          H(J,J)=T1
          T1=CS(J)*G(J)
          G(J+1)=-SN(J)*G(J)
          G(J)=T1
!
          NITER=NITER+1
          JUSED=J
          RNRM=ABS(G(J+1))/BNRM
          RESID=RNRM
          IF(RNRM .LE. TOL .OR. HN .EQ. 0.0_LDP)EXIT
        END DO
!
! Back-solve the small triangular least-squares state and update the
! large vector once per restart.
!
        Y=0.0_LDP
        DO I=JUSED,1,-1
          T1=G(I)
          DO L=I+1,JUSED
            T1=T1-H(I,L)*Y(L)
          END DO
          IF(H(I,I) .EQ. 0.0_LDP)THEN
            INFO=4
            GOTO 900
          END IF
          Y(I)=T1/H(I,I)
        END DO
        DO L=1,JUSED
          XV=XV+Y(L)*V(:,L)
        END DO
!
        IF(RNRM .LE. TOL)THEN
          INFO=0
          GOTO 800
        END IF
      END DO
      INFO=1
!
  800 CONTINUE
! Report the true final preconditioned residual. The Arnoldi estimate
! is exact in exact arithmetic, but production rows are ill conditioned
! enough that the census should not rely on that identity.
      IF(INFO .EQ. 0 .OR. INFO .EQ. 1)THEN
        CALL BLK_MATVEC(ND,NCF,TA,TB,TC,CUP,KUP,FRESH,XV,W)
        Z=B-W
        CALL BLK_PRECOND_VEC(ND,NCF,TA,TB,TC,FRESH,Z,W,IPS)
        IF(IPS .EQ. 0)RESID=SQRT(DOT_PRODUCT(W,W))/BNRM
      END IF
      DO K=1,NCF
        DO I=1,ND
          X(I,K)=XV(I+(K-1)*ND)
        END DO
      END DO
!
  900 CONTINUE
      DEALLOCATE(XV,B,Z,W,Y,V,H,CS,SN,G)
      RETURN
      END SUBROUTINE MOM_BLK_GMRES_SOLVE
!
!-----------------------------------------------------------------------
!
      SUBROUTINE BLK_MATVEC(ND,NCF,TA,TB,TC,CUP,KUP,FRESH,X,Y)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ND,NCF,KUP(ND,NCF)
      REAL(KIND=LDP), INTENT(IN) :: TA(ND,NCF),TB(ND,NCF)
      REAL(KIND=LDP), INTENT(IN) :: TC(ND,NCF),CUP(ND,NCF)
      REAL(KIND=LDP), INTENT(IN) :: X(ND*NCF)
      LOGICAL, INTENT(IN) :: FRESH(NCF)
      REAL(KIND=LDP), INTENT(OUT) :: Y(ND*NCF)
      INTEGER I,K,L,LU
!
      DO K=1,NCF
        DO I=1,ND
          L=I+(K-1)*ND
          IF(.NOT. FRESH(K))THEN
            Y(L)=X(L)
          ELSE
            Y(L)=TB(I,K)*X(L)
            IF(I .GT. 1)Y(L)=Y(L)+TA(I,K)*X(L-1)
            IF(I .LT. ND)Y(L)=Y(L)+TC(I,K)*X(L+1)
            LU=KUP(I,K)
            IF(LU .GE. 1 .AND. LU .LE. NCF)
     1        Y(L)=Y(L)-CUP(I,K)*X(I+(LU-1)*ND)
          END IF
        END DO
      END DO
      RETURN
      END SUBROUTINE BLK_MATVEC
!
!-----------------------------------------------------------------------
!
      SUBROUTINE BLK_PRECOND_VEC(ND,NCF,TA,TB,TC,FRESH,R,Z,INFO)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ND,NCF
      REAL(KIND=LDP), INTENT(IN) :: TA(ND,NCF),TB(ND,NCF)
      REAL(KIND=LDP), INTENT(IN) :: TC(ND,NCF),R(ND*NCF)
      LOGICAL, INTENT(IN) :: FRESH(NCF)
      REAL(KIND=LDP), INTENT(OUT) :: Z(ND*NCF)
      INTEGER, INTENT(OUT) :: INFO
!
      INTEGER I,K,L,IOS
      REAL(KIND=LDP), ALLOCATABLE :: DL(:),DD(:),DU(:),RR(:)
!
      INFO=0
      Z=R
      ALLOCATE(DL(MAX(ND-1,1)),DD(ND),DU(MAX(ND-1,1)),RR(ND),
     1  STAT=IOS)
      IF(IOS .NE. 0)THEN
        INFO=IOS
        RETURN
      END IF
      DO K=1,NCF
        IF(.NOT. FRESH(K))CYCLE
        DO I=1,ND-1
          DL(I)=TA(I+1,K)
          DU(I)=TC(I,K)
        END DO
        DO I=1,ND
          DD(I)=TB(I,K)
          RR(I)=R(I+(K-1)*ND)
        END DO
        CALL DGTSV(ND,1,DL,DD,DU,RR,ND,IOS)
        IF(IOS .NE. 0)THEN
          INFO=IOS
          EXIT
        END IF
        DO I=1,ND
          Z(I+(K-1)*ND)=RR(I)
        END DO
      END DO
      DEALLOCATE(DL,DD,DU,RR)
      RETURN
      END SUBROUTINE BLK_PRECOND_VEC
!
      END MODULE MOM_BLK_GMRES_MOD
