!
! Subroutine to compute the quadrature weights for a modified cubic rule.
! Reference - A Nordlund (Methods in Radiative Transfer, W Kallofen).
! The routine makes specific assumptions concerning the behaviour of the
! flux at the boundary, and would need to be moified to obtain  other
! quantities other than 1st moment of the intensity. The program assumes
! that the first point corresponds to mu=1.0 .
!
! Note that these weights should not be normalized in the usual fashion.
! Physically, we dont expect V (the flux) to be constant with respect to mu.
! For small mu, we expect a that V is proportional to mu.
!
! D1, and R2 are work vectors.
!
	SUBROUTINE HWEIGHT_V2(X,dX,W,N)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Created 11-Jun-2023 - Based on HWEIGHT
!                       Altered to allow cases when X is very close to 1.
!                       Prevents overflow.
!                       Changed to V2 as dX added to call.
!
	INTEGER N
	REAL(KIND=LDP) X(N)        !MU
	REAL(KIND=LDP) dX(N)       !dMU
	REAL(KIND=LDP) W(N)        !Quadrature weight
!
	REAL(KIND=LDP) H,HN,RF,RE,SUM
	REAL(KIND=LDP) T1,T2
!
	INTEGER I
	INTEGER LS
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
	LOGICAL, PARAMETER :: CHECK=.FALSE.
!
	W(:)=0.0_LDP
!
	H=dX(1)
	HN=dX(2)
	RF=dX(1)+dX(2)
	W(1)=0.5_LDP*H*(X(1)-H/6.0_LDP)
	1       -H*H/RF/12.0_LDP*((2.0_LDP+HN/H)*X(1)-X(2))
	W(2)=0.5_LDP*H*(X(2)+H/6.0_LDP)
	1       +H*H/RF/12.0_LDP*(2.0_LDP+H/HN+HN/H)*X(1)
	W(3)=-H*H/RF/12.0_LDP*(X(2)+H/HN*X(1))
!
	DO I=3,N-1
	  H=dX(I-1)
	  RF=dX(I-2)+dX(I-1)
	  RE=dX(I-1)+dX(I)
	  W(I-2)=W(I-2)-X(I-1)*H*H/12.0_LDP/RF
	  W(I-1)=W(I-1)+0.5_LDP*H*(X(I-1)-H/6.0_LDP)
	1              +X(I)*H*H/12.0_LDP/RE
	  W(I)=W(I)+0.5_LDP*H*(X(I)+H/6.0_LDP)
	1              +X(I-1)*H*H/12.0_LDP/RF
	  W(I+1)=W(I+1)-X(I)*H*H/12.0_LDP/RE
	END DO
!
! Assumes that V(mu=0)=0 and mu.dV(mu)=0 at mu=0.
! Subtracting the X to get dX is fine for high N.
!
	IF(X(N) .EQ. 0.0_LDP)THEN
	  H=dX(N-1)
	  W(N-1)=W(N-1)+H*H/3.0_LDP
	ELSE
!
! Since V(mu=0) is zero, we dont actually need a ray with mu=0 since the weight
! is automatically zero.
!
! Integral from X(N-1) to X(N)
!
          H=X(N-1)-X(N)
          RF=X(N-2)-X(N)
          RE=X(N-1)                             !As X(N+1)=0
	  W(N-2)=W(N-2)-X(N-1)*H*H/12.0_LDP/RF
	  W(N-1)=W(N-1)+0.5_LDP*H*(X(N-1)-H/6.0_LDP)
	1              +X(N)*H*H/12.0_LDP/RE
	  W(N)=W(N)+0.5_LDP*H*(X(N)+H/6.0_LDP)
	1              +X(I-1)*H*H/12.0_LDP/RF
!
! Integral from X(N) to X(N+1)
!
	  H=X(N)
	  W(N)=W(N)+H*H/3.0_LDP
	  WRITE(lUER,*)'Warning - Extrapolation to zero required in HWEIGHT_V2'
	END IF
!
! Ensure that the weights have the correct normalization (but dont
! perform the normalization).
!
	SUM=0.0_LDP
	DO I=1,N
	  SUM=SUM+W(I)*X(I)
	END DO
	SUM=1.0_LDP/SUM/3.0_LDP
	IF(ABS(SUM-1.0_LDP) .GT. 1.0E-12_LDP)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)' Warning - weights require normalization in HWEIGHT_V2'
	  WRITE(LUER,*)' Expected SUM is 1: SUM=',SUM
	END IF
!
	IF(CHECK)THEN
	  WRITE(6,*)'Check of H weights'
          WRITE(6,'(18X,A,5X,A,6XX,A,13X,A)')'Mu','dMU(diff)','dMU(Tay)','W'
	  DO LS=1,N-1
	    T1=X(LS)-X(LS+1)
	    WRITE(6,'(F20.16,3ES14.6)')X(LS),T1,dX(LS),W(LS)
	  END DO
	  T1=0.0_LDP
	  WRITE(6,'(F20.16,3ES14.6)')X(LS),T1,dX(LS),W(LS)
	  FLUSH(UNIT=6)
	END IF
!
	RETURN
	END
