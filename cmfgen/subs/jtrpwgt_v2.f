!
! Subroutine to compute the quadrature weights for the J integration.
! A trapazoidal rule is used. Assumes X(I) > X(I+1).
!
	SUBROUTINE JTRPWGT_V2(X,dX,W,N)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Aletered 24-May-1996 - Call to DP_ZERO removed
!                        EROR_LU etc installed.
! Created  17-May-1989 - Based on HWEIGHT
!
	INTEGER N
	REAL(KIND=LDP) X(N),dX(N),W(N)
!
	REAL(KIND=LDP) H,T1,T2,SUM,XSUM
	INTEGER I
	INTEGER ERROR_LU,LUER
	LOGICAL, SAVE :: CHECK=.FALSE.
	EXTERNAL ERROR_LU
!
	W(:)=0.0_LDP
	DO I=1,N-1
	  W(I)=W(I)+0.5_LDP*dX(I)
	  W(I+1)=W(I+1)+0.5_LDP*dX(I)
	END DO
!
! Assumes that j'=0
!
	IF(X(N) .NE. 0)THEN
!
! Integral from X(N) to 0
!
	  T1=X(N)/( X(N-1)*X(N-1)-X(N)*X(N) )
	  W(N-1)=W(N-1)-2.0_LDP*T1*X(N)*X(N)/3.0_LDP
	  W(N)=W(N)+T1*(X(N-1)*X(N-1)-X(N)*X(N)/3.0_LDP)
!
	END IF
!
! Ensure that the weights have the correct normalization (but dont
! perform the normalization). Two checks are done to insure that
! the correct answer is given for a linear varaition. Because of the
! assumption that du/dmu(mu=0)=0, we have to fiddle with the last check.
!
	SUM=0.0_LDP
	DO I=1,N
	  SUM=SUM+W(I)
	END DO
!
	XSUM=0.0_LDP
	DO I=1,N
	  XSUM=XSUM+W(I)*X(I)
	END DO
	IF(X(N) .EQ. 0)THEN
	  T1=0.5_LDP
	ELSE
	  T1=X(N)/( X(N-1)*X(N-1)-X(N)*X(N) )
	  T1=T1*( X(N)*(X(N-1)*X(N-1)-X(N)*X(N)/3.0_LDP)
	1        -2.0_LDP*X(N)*X(N)*X(N-1)/3.0_LDP )
	  T1=0.5_LDP*(1.0_LDP-X(N)*X(N))+T1
	END IF
	XSUM=XSUM/T1
	IF(ABS(XSUM-1.0_LDP) .GT. 1.0E-12_LDP)THEN
	   LUER=ERROR_LU()
	   WRITE(LUER,*)' Warning - weights require normalization in JTRPWGT_V2'
	   WRITE(LUER,*)' Expected normalized  Sum(W) is 1:  Sum(w)-1 =',SUM-1.0D0
	   WRITE(LUER,*)' Expected bormalized Sum(xW) is 1: Sum(xw)-1 =',XSUM-1.0D0
	END IF
!
! Check just used for initial developemnt.
!
	IF(CHECK)THEN
	  WRITE(6,*)'Check of J weights in JTRPWGT_V2'
	  WRITE(6,'(18X,A,11X,A,14X,A,21X,A,18X,A)')'MU','dMU','dMU(acc)','W','Wsum'
	  SUM=0.0_LDP
	  DO I=1,N-1
	    T1=X(I)-X(I+1)
	    SUM=SUM+W(I)
	    WRITE(6,'(F20.16,ES14.6,3ES22.14)')X(I),T1,dX(I),W(I),SUM
	  END DO
	  I=N;  SUM=SUM+W(I)
	  WRITE(6,'(F20.16,ES14.6,3ES22.14)')X(I),0.0D0,dX(I),W(I),SUM
	END IF
!
	RETURN
	END
