!
! Subroutine to compute the quadrature weights for the K integration
! (the second moment of the radiation field). A trapazoidal rule is used.
!
	SUBROUTINE KTRPWGT_V2(X,dX,W,N)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
!
! Created 12-Jun-2023 : Basesd oh KTRPWGT
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
	W(:)=0.0_LDP
	DO I=1,N-1
	  T1=X(I)*X(I)+2.0_LDP*X(I)*X(I+1)+X(I+1)*X(I+1)
	  W(I)=W(I) + dX(I)*(T1+2.0_LDP*X(I)*X(I))/12.0_LDP
	  W(I+1)=W(I+1) + dX(I)*(T1+2.0_LDP*X(I+1)*X(I+1))/12.0_LDP
	END DO
!
! Assumes that du/dmu=0 at mu=0.
!
	IF(X(N) .NE. 0.0_LDP)THEN
!
! Integral from X(N) to 0
!
	  T1=( X(N)**3 )/( X(N-1)*X(N-1)-X(N)*X(N) )
	  W(N-1)=W(N-1)-2.0_LDP*T1*X(N)*X(N)/15.0_LDP
	  W(N)=W(N)+T1*(X(N-1)*X(N-1)/3.0_LDP-0.2_LDP*X(N)*X(N))
!
	END IF
!
! Ensure that the weights have the correct normalization (but dont
! perform the normalization). Two checks are done to insure that
! the correct answer is given for a linear variation. Because of the
! assumption that du/dmu=0 at mu=0 we have to fiddle with the last check.
!
	SUM=0.0_LDP
	DO I=1,N
	  SUM=SUM+W(I)
	END DO
	SUM=1.0_LDP/SUM/3.0_LDP
!
	XSUM=0.0_LDP
	DO I=1,N
	  XSUM=XSUM+W(I)*X(I)
	END DO
	IF(X(N) .EQ. 0)THEN
	  T1=0.25_LDP
	ELSE
	  T1=( X(N)**3 )/( X(N-1)*X(N-1)-X(N)*X(N) )
	  T1=T1*( X(N)*(X(N-1)*X(N-1)/3.0_LDP-0.2_LDP*X(N)*X(N))
	1          -2.0_LDP*X(N)*X(N)*X(N-1)/15.0_LDP )
	  T1=0.25_LDP*( 1.0_LDP-X(N)**4 ) + T1
	END IF
	XSUM=XSUM/T1
	IF(ABS(XSUM-1.0_LDP) .GT. 1.0E-12_LDP .OR. ABS(SUM-1.0_LDP) .GT. 1.0E-12_LDP)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)' Warning - weights require normalization in KTRPWGT_V2'
	  WRITE(LUER,*)' Expected normalized  Sum(W) is 1:  Sum(w)-1 =',SUM-1.0D0
	  WRITE(LUER,*)' Expected normalized Sum(xW) is 1: Sum(xw)-1 =',XSUM-1.0D0
	END IF
!
        IF(CHECK)THEN
          WRITE(6,*)'Check on K weights in KTRPWGT_V2'
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
