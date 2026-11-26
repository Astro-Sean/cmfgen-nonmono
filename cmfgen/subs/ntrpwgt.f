
C Subroutine to compute the quadrature weights for the integration of
C N. A linear variation of V (the flux variable) with mu is assumed.
C If the last point does not correspond to mu=0, the routine assumes
C that the flux for mu=0 is zero. The routine should only be used to
C calculate the third momemnt of the intensity. The program assumes
C that the first point corresponds to mu=1.0 .
C
C Note that these weights are not be normalized in the usual fashion.
C Physically, we dont expect V (the flux) to be constant with respect to mu.
C For small mu, we expect a that V is proportional to mu.
C
	SUBROUTINE NTRPWGT(X,W,N)
	USE SET_KIND_MODULE
	IMPLICIT NONE
C
C Altered 25-May-1996 - Call to DP_ZERO removed.
C                       ERROR_LU inserted.
C Created 26-Apr-1989 - Based on HWEIGHT
C
	INTEGER N
	REAL(KIND=LDP) X(N),W(N)
C
	REAL(KIND=LDP) T1,T2,SUM
	INTEGER I
	LOGICAL, SAVE :: CHECK=.FALSE.
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
C
	W(:)=0.0_LDP
C
	DO I=1,N-1
	  T1=0.25_LDP*(X(I)+X(I+1)) * (X(I)*X(I)+X(I+1)*X(I+1))
	  T2=( (X(I)**4) + (X(I)**3)*X(I+1) + (X(I)*X(I+1))**2
	1      + (X(I+1)**3)*X(I) + (X(I+1)**4) )/5.0_LDP
	  W(I)=W(I)+T2-X(I+1)*T1
	  W(I+1)=W(I+1)-T2+T1*X(I)
	END DO
C
C Assumes that V(mu=0)=0.
C
	IF(X(N) .NE. 0.0_LDP)THEN
C
C Since V(mu=0) is zero, we dont actually need a ray with mu=0 since the weight
C is automatically zero.
C
C Integral from X(N-1) to X(N)
C
	  W(N)=W(N)+(X(N)**4)/5.0_LDP
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
	SUM=0.2_LDP/SUM
	IF(ABS(SUM-1.0_LDP) .GT. 1.0E-12_LDP)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)' Warning - weights require normalization in HWEIGHT'
	END IF
C
	SUM=0.0_LDP
	DO I=1,N
	  SUM=SUM+W(I)
	END DO
	IF(X(N) .EQ. 0.0_LDP)THEN
	  T1=0.25_LDP
	ELSE
	  T1=0.25_LDP*( 1.0_LDP-(X(N)**4) )+0.2_LDP*(X(N)**4)
	END IF
	SUM=T1/SUM
	IF(ABS(SUM-1.0_LDP) .GT. 1.0E-12_LDP)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)' Warning - weights require normalization in HWEIGHT'
	END IF
C
	IF(CHECK)THEN
	  WRITE(6,*)'Check on N weights in NTRPWGT'
	  WRITE(6,'(18X,A,11X,A,13X,A)')'MU','dMU','W'
	  DO I=1,N-1
	    T1=X(I)-X(I+1)
	    WRITE(6,'(F20.16,3ES14.6)')X(I),T1,W(I)
	  END DO
	  I=N;  WRITE(6,'(F20.16,3ES14.6)')X(I),0.0D0,W(I)
	END IF
C
	RETURN
	END
