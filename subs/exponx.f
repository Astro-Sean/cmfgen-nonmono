C
C The function EXPONX is given by (1.0-EXP(-X))/X.
C This functions is called to allow for cancellation when X is small.
C
	FUNCTION  EXPONX(X)
	USE SET_KIND_MODULE
	IMPLICIT NONE
C
C ALtered 24-May-1996 :  DOUBLE PRECISION replaced by REAL(KIND=LDP)
C Altered 26-May-1988 : Exponential no longer computed if X is greater than
C                        40. Necessary to overcome a bug with Dec software.

	REAL(KIND=LDP) EXPONX,X
C
	IF( ABS(X) .LT. 1.0E-03_LDP )THEN
	  EXPONX=1.0_LDP-X*(  0.5_LDP-X/6.0_LDP*( 1.0_LDP-X/4.0_LDP )  )
	ELSE IF(X .LT. 40)THEN
	  EXPONX=( 1.0_LDP-EXP(-X) )/X
	ELSE
	  EXPONX=1.0_LDP/X
	END IF
C
	RETURN
	END


C
C The function d_EXPONX_dX is given by d[ (1.0-EXP(-X))/X ]/dX. This
C function is called to allow for cancellation when X is small.
C
C Altered 26-May-88 - Exponential no longer computed if X is greater than
C                     40. Necessary to overcome a bug with Dec software.
C
	FUNCTION  d_EXPONX_dX(X)
	USE SET_KIND_MODULE
	IMPLICIT NONE
	REAL(KIND=LDP) d_EXPONX_dX,X,Y
C
	IF( ABS(X) .LT. 1.0E-03_LDP )THEN
	  d_EXPONX_dX=-0.5_LDP+X*( 1.0_LDP-X*(0.375_LDP-0.1_LDP*X) )/3.0_LDP
	ELSE
	  IF(X .LT. 40)THEN
	    Y=EXP(-X)
	  ELSE
	    Y=0.0_LDP
	  END IF
	  d_EXPONX_dX=( Y-(1.0_LDP-Y)/X )/X
	END IF
C
	RETURN
	END
