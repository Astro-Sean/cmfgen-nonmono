
! Subroutine to compute the quadrature weights for the integration of
! N. A linear variation of V (the flux variable) with mu is assumed.
! If the last point does not correspond to mu=0, the routine assumes
! that the flux for mu=0 is zero. The routine should only be used to
! calculate the third momemnt of the intensity. The program assumes
! that the first point corresponds to mu=1.0 .
!
! Note that these weights are not be normalized in the usual fashion.
! Physically, we dont expect V (the flux) to be constant with respect to mu.
! For small mu, we expect a that V is proportional to mu.
!
	SUBROUTINE NTRPWGT_V2(X,dX,W,N)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Altered 25-May-1996 - Call to DP_ZERO removed.
!                       ERROR_LU inserted.
! Created 26-Apr-1989 - Based on HWEIGHT
!
	INTEGER N
	REAL(KIND=LDP) X(N),dX(N),W(N)
!
	REAL(KIND=LDP) T1,T2,SUM,XSUM
	INTEGER I
	LOGICAL, SAVE :: CHECK=.FALSE.
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
!
	W(:)=0.0_LDP
	DO I=1,N-1
	  T1=X(I)*(2*X(I+1)*X(I+1)+3*X(I)*X(I+1)+4*X(I)*X(I))   + X(I+1)**3
	  T2=X(I+1)*(2*X(I)*X(I)+3*X(I)*X(I+1)+4*X(I+1)*X(I+1)) + X(I)**3
	  W(I)=W(I)+0.05_LDP*T1*dX(I)
	  W(I+1)=W(I+1)+0.05_LDP*T2*dX(I)
	END DO
!
! Assumes that V(mu=0)=0.
!
	IF(X(N) .NE. 0.0_LDP)THEN
!
! Since V(mu=0) is zero, we dont actually need a ray with mu=0 since the weight
! is automatically zero.
!
! Integral from X(N-1) to X(N)
!
	  W(N)=W(N)+(X(N)**4)/5.0_LDP
!
	END IF
!
! Ensure that the weights have the correct normalization (but dont
! perform the normalization). Two checks are done to insure that
! the correct answer is given for a linear varaition. Because of the
! assumption that V(mu=0)=0, we have to fiddle with the last check.
!
	XSUM=0.0_LDP
	DO I=1,N
	  XSUM=XSUM+W(I)*X(I)
	END DO
	XSUM=XSUM/0.2_LDP
!
	SUM=0.0_LDP
	DO I=1,N
	  SUM=SUM+W(I)
	END DO
	IF(X(N) .EQ. 0.0_LDP)THEN
	  T1=0.25_LDP
	ELSE
	  T1=0.25_LDP*( 1.0_LDP-(X(N)**4) )+0.2_LDP*(X(N)**4)
	END IF
	SUM=SUM/T1
	IF(ABS(SUM-1.0_LDP) .GT. 1.0E-12_LDP .OR. ABS(XSUM-1.0_LDP) .GT. 1.0E-12_LDP)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)' Warning - weights require normalization in NTRPWGT_V2'
	  WRITE(LUER,*)' Expected normalized  Sum(W) is 1:  Sum(w)-1 =',SUM-1.0D0
	  WRITE(LUER,*)' Expected normalized Sum(xW) is 1: Sum(xw)-1 =',XSUM-1.0D0
	END IF
!
	IF(CHECK)THEN
	  WRITE(6,*)'Check on N weights in NTRPWGT_V2'
	  WRITE(6,'(18X,A,11X,A,14X,A,21X,A,21X,A)')'MU','dMU','dMU(acc)','W','Wsum'
	  SUM=0.0_LDP
	  DO I=1,N-1
	    T1=X(I)-X(I+1)
	    SUM=SUM+W(I)
	    WRITE(6,'(F20.16,ES14.6,3ES22.14)')X(I),T1,dX(I),W(I),SUM
	  END DO
	  I=N; SUM=SUM+W(I)
	  WRITE(6,'(F20.16,ES14.6,3ES22.14)')X(I),0.0D0,dX(I),W(I),SUM
	END IF
!
	RETURN
	END
