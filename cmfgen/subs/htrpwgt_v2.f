!
! Subroutine to compute the quadrature weights H assiming a linear
! depndence on I.  Assumes X(I) > X(I+1).
!
	SUBROUTINE HTRPWGT_V2(X,dX,W,N)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Created 12-Jun-2023 : Basesd oh HTRPWGT
!                       Computes more accrate dX when X ~ 1
!
	INTEGER N
	REAL(KIND=LDP) X(N),dX(N),W(N)
!
	REAL(KIND=LDP) T1,T2,SUM,XSUM
	INTEGER I
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
	LOGICAL, SAVE :: CHECK=.FALSE.
!
! Compute the weights.
!
	W(:)=0.0_LDP
	DO I=1,N-1
	  W(I)  =W(I)   + dX(I)*(2.0_LDP*X(I)+X(I+1))/6.0_LDP
	  W(I+1)=W(I+1) + dX(I)*(X(I)+2.0_LDP*X(I+1))/6.0_LDP
	END DO
!
! Assumes that V(mu=0)=0 and mu.V'(mu)=0 at mu=0.
!
	IF(X(N) .NE. 0.0_LDP)THEN
!
! Since V(mu=0) is zero, we dont actually need a ray with mu=0 since the weight
! is automatically zero.
!
! Integral from X(N-1) to X(N)
!
	  W(N)=W(N)+X(N)*X(N)/3.0_LDP
!
	END IF
!
! Ensure that the weights have the correct normalization (but dont
! perform the normalization). Two checks are done to insure that
! the correct answer is given for a linear varaition. Because of the
! assumption that V(mu=0)=0, we have to fiddle with the last check.
!
	SUM=0.0_LDP
	DO I=1,N
	  SUM=SUM+W(I)*X(I)
	END DO
	SUM=1.0_LDP/SUM/3.0_LDP
!
	XSUM=0.0_LDP
	DO I=1,N
	  XSUM=XSUM+W(I)
	END DO
	IF(X(N) .EQ. 0)THEN
	  T1=0.5_LDP
	ELSE
	  T1=0.5_LDP*(1.0_LDP-X(N)*X(N))+(X(N)**2)/3.0_LDP
	END IF
	XSUM=XSUM/T1
	IF(ABS(XSUM-1.0_LDP) .GT. 1.0E-12_LDP .OR. ABS(SUM-1.0_LDP) .GT. 1.0E-12_LDP)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)' Warning - weights require normalization in HTRPWGT_V2'
	  WRITE(LUER,*)' Expected normalized  Sum(W) is 1:  Sum(w)-1 =',SUM-1.0D0
          WRITE(LUER,*)' Expected normalized Sum(xW) is 1: Sum(xw)-1 =',XSUM-1.0D0
	END IF
!
	IF(CHECK)THEN
	  WRITE(6,*)'Check on H weights in HTRPWGT_V2'
	  WRITE(6,'(18X,A,11X,A,14X,A,21X,A,18X,A)')'MU','dMU','dMU(acc)','W','Wsum'
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
