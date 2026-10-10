C
C Subroutine to compute the weights for a "TRAPAZODIAL" or "SIMPSON"
C quadrature. The data points may be unequally spaced but with
C SIMPSON quadrature the data points must have a non zero spacing.
C The simpson rules should be used with an odd number of data points.
C If the number of points is even, the use a trapazoidal rule (SIMPUNEQ)
C of quadratic fit (SIMPEQ) to the last interval.
C
	SUBROUTINE WEIGHT(U,W,N)
	USE SET_KIND_MODULE
	IMPLICIT NONE
C
C Altered 02-Jul-1998 - Missed LUER and ERROR_LU inserted.
C Altered 28-May-1996 - ERROR_LU installed.
C                       Call to DP_ZERO removed.
C Altered 11-May-1989 - Bug fixed in SIMPEQ
C                     - Cleaned (Implicit none installed, T array removed)
C                     - Value of N is checked.
C                     - All routines checked to see if they give the correct
C                       answers for a linear function (also quadratic and cubic).
C Altered 28-SEP-1982 (6.9 in simuneq replaced by 6.0)
C
	INTEGER N,I
	REAL(KIND=LDP) U(N),W(N)
	REAL(KIND=LDP) T1,T2,TT
C
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU

C
C U is array of data points.
C W is the array the weights are returned in.
C T is a trash array which may also be u.
C
C
	ENTRY TRAPEQ(U,W,N)
	  IF(N .LT. 2)THEN
	    LUER=ERROR_LU()
	    WRITE(LUER,*)'Error in TRAPEQ: N < 2'
	    STOP
	  END IF
	  T1=ABS(U(2)-U(1))
	  W(1)=0.5_LDP*T1
	  W(N)=W(1)
	  DO I=2,N-1
	    W(I)=T1
	  END DO
	RETURN
C
C
	ENTRY TRAPUNEQ(U,W,N)
	  IF(N .LT. 2)THEN
	    LUER=ERROR_LU()
	    WRITE(LUER,*)'Error in TRAPUNEQ: N < 2'
	    STOP
	  END IF
	  DO I=2,N-1
	     W(I)=ABS(U(I+1)-U(I-1))*0.5_LDP
	  END DO
	  W(1)=ABS(U(2)-U(1))*0.5_LDP
	  W(N)=ABS(U(N)-U(N-1))*0.5_LDP
	RETURN
C
C
	ENTRY SIMPEQ(U,W,N)
	  IF(N .LT. 3)THEN
	    LUER=ERROR_LU()
	    WRITE(LUER,*)'Error in SIMPEQ: N < 3'
	    STOP
	  END IF
	  T1=ABS(U(2)-U(1))/3.0_LDP
	  T2=4.0_LDP*T1
	  W(:)=0.0_LDP
	  DO I=1,N-2,2
	    W(I)=W(I)+T1
	    W(I+1)=T2
	    W(I+2)=W(I+2)+T1
	  END DO
	  IF( (N/2)*2 .EQ. N)THEN
	    W(N)=1.25_LDP*T1
	    W(N-1)=W(N-1)+2.0_LDP*T1
	    W(N-2)=W(N-2)-0.25_LDP*T1
	  END IF
	RETURN
C
	ENTRY SIMPUNEQ(U,W,N)
	  IF(N .LT. 3)THEN
	    LUER=ERROR_LU()
	    WRITE(LUER,*)'Error in SIMPUNEQ: N < 3'
	    STOP
	  END IF
	  W(:)=0.0_LDP
	  DO I=1,N-2,2
	    T1=ABS(U(I+1)-U(I))
	    T2=ABS(U(I+2)-U(I+1))
	    TT=T1+T2
	    W(I)=W(I)+(2.0_LDP*T1-T2)*TT/(6.0_LDP*T1)
	    W(I+1)=TT*TT*TT/(6.0_LDP*T1*T2)
	    W(I+2)=(2.0_LDP*T2-T1)*TT/(6.0_LDP*T2)
	  END DO
	  IF(N/2*2 .EQ. N)THEN
	    W(N)=ABS(U(N)-U(N-1))*0.5_LDP
	    W(N-1)=W(N-1)+W(N)
	  END IF
	RETURN
C
C Compute quadrature weights using Simpons rule when data points are
C equally spaced, else uses the trapazoidal rule. The routine works
C "backwards" down the arry.
C
	ENTRY SMPTRP(U,W,N)
	  IF(N .LT. 3)THEN
	    LUER=ERROR_LU()
	    WRITE(LUER,*)'Error in SIMTRP: N < 4'
	    STOP
	  END IF
	  W(:)=0.0_LDP
	  I=N
	  DO WHILE (I .GE. 2)
	    T1=ABS(U(I-1)-U(I))
	    IF(I .NE. 2)T2=ABS(ABS(U(I-2)-U(I-1))/T1-1.0_LDP)
	    IF(T2 .LT. 1.0E-06_LDP .AND. I .NE. 2)THEN
	      W(I)=W(I)+T1/3.0_LDP
	      W(I-1)=W(I-1)+4.0_LDP*T1/3.0_LDP
	      W(I-2)=W(I-2)+T1/3.0_LDP
	      I=I-2
	    ELSE
	      W(I)=W(I)+T1*0.5_LDP
	      W(I-1)=W(I-1)+T1*0.50_LDP
	      I=I-1
	    END IF
	  END DO
	RETURN

	END
