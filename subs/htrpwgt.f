!
! Subroutine to compute the quadrature weights H assiming a linear
! depndence on I.  Assumes X(I) > X(I+1).
!
	SUBROUTINE HTRPWGT(X,W,N)
	USE SET_KIND_MODULE
	IMPLICIT NONE
C
C Altered 24-May-1996 - Call to DP_ZERO removed.
C                       ERROR_LU, LUER installed.
C Altered 25-Feb-1987 - Normalization implemented.
C Altered 14-Jan-87 - Program nolonger assumes that X(N)=0, since can
C                     assume X(N+1)=0 with W(N+1)=0.0 with out loss of
C                     accuracy.
C Created 25-Nov-1986 (Based on KWEIGHT)
C
	INTEGER N,I
	REAL(KIND=LDP) X(N),W(N),T1,T2,SUM
C
	INTEGER ERROR_LU,LUER
	LOGICAL, SAVE :: CHECK = .FALSE.
	EXTERNAL ERROR_LU
C
	W(:)=0.0_LDP
	DO I=1,N-1
	  T1=0.5_LDP*(X(I)+X(I+1))
	  T2=(X(I)*X(I)+X(I)*X(I+1)+X(I+1)*X(I+1))/3.0_LDP
	  W(I)=W(I)+T2-X(I+1)*T1
	  W(I+1)=W(I+1)-T2+T1*X(I)
	END DO
C
C Assumes that V(mu=0)=0 and mu.V'(mu)=0 at mu=0.
C
	IF(X(N) .NE. 0.0_LDP)THEN
C
C Since V(mu=0) is zero, we dont actually need a ray with mu=0 since the weight
C is automatically zero.
C
C Integral from X(N-1) to X(N)
C
	  W(N)=W(N)+X(N)*X(N)/3.0_LDP
C
	END IF
C
C Ensure that the weights have the correct normalization (but dont
C perform the normalization). Two checks are done to insure that
C the correct answer is given for a linear varaition. Because of the
C assumption that V(mu=0)=0, we have to fiddle with the last check.
C
	SUM=0.0_LDP
	DO I=1,N
	  SUM=SUM+W(I)*X(I)
	END DO
	SUM=1.0_LDP/SUM/3.0_LDP
	IF(ABS(SUM-1.0_LDP) .GT. 1.0E-12_LDP)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Warning - weights require normalization in HWEIGHT'
	  WRITE(LUER,*)'Scale=',SUM
	END IF
C
	SUM=0.0_LDP
	DO I=1,N
	  SUM=SUM+W(I)
	END DO
	IF(X(N) .EQ. 0)THEN
	  T1=0.5_LDP
	ELSE
	  T1=0.5_LDP*(1.0_LDP-X(N)*X(N))+(X(N)**2)/3.0_LDP
	END IF
	SUM=T1/SUM
	IF(ABS(SUM-1.0_LDP) .GT. 1.0E-12_LDP)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)' Warning - weights require normalization in HWEIGHT'
	END IF
C
	IF(CHECK)THEN
	  WRITE(6,*)'Check on H weights in HTRPWGT'
	  WRITE(6,'(A,A,A)')'MU','dMU','W'
	  DO I=1,N-1
	    WRITE(6,'(F20.16,2ES14.6)')X(I),X(I)-X(I+1),W(I)
	  END DO
	  I=N;  WRITE(6,'(F20.16,3ES14.6)')X(I),0.0D0,W(I)
	END IF
C
	RETURN
	END
