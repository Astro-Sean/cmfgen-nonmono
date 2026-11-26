C
C The function EXPONX is given by (1.0-EXP(-X))/X.
C This functions is called to allow for cancellation when X is small.
C
	MODULE EXPONX_GUARD_MOD
	  IMPLICIT NONE
C
C Activation counters for the X < -700 bound in EXPONX and
C d_EXPONX_dX.  Diagnostic only (SOB_check WP3): EXP(-X) overflows
C float64 for X <~ -709, which a sigma_true near -1 can reach through
C the 1e-4 arm of the tau floor (GAM ~ GAMH/1e-4 for a weakly
C inverted GAMH).  Unit-test drivers read the counters via USE;
C production runs leave them zero when the bound never engages.
C
	  INTEGER :: N_EXPONX_GUARD = 0
	  INTEGER :: N_DEXPONX_GUARD = 0
	END MODULE EXPONX_GUARD_MOD
C
	FUNCTION  EXPONX(X)
	USE SET_KIND_MODULE
	USE EXPONX_GUARD_MOD
	IMPLICIT NONE
C
C ALtered 24-May-1996 :  DOUBLE PRECISION replaced by REAL(KIND=LDP)
C Altered 26-May-1988 : Exponential no longer computed if X is greater than
C                        40. Necessary to overcome a bug with Dec software.
C Altered 27-Nov-2026 : Bounded return for X < -700 added; the
C                        exponential's argument is pinned at 700 (margin
C                        below the ~709.78 overflow) while the outer X is
C                        kept, so the return stays finite and retains the
C                        1/X shape.  Results for |X| <= 700 are unchanged.
C
	REAL(KIND=LDP) EXPONX,X
C
	IF( ABS(X) .LT. 1.0E-03_LDP )THEN
	  EXPONX=1.0_LDP-X*(  0.5_LDP-X/6.0_LDP*( 1.0_LDP-X/4.0_LDP )  )
	ELSE IF(X .LT. -700.0_LDP)THEN
	  N_EXPONX_GUARD=N_EXPONX_GUARD+1
	  EXPONX=( 1.0_LDP-EXP(700.0_LDP) )/X
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
	USE EXPONX_GUARD_MOD
	IMPLICIT NONE
	REAL(KIND=LDP) d_EXPONX_dX,X,Y
C
	IF( ABS(X) .LT. 1.0E-03_LDP )THEN
	  d_EXPONX_dX=-0.5_LDP+X*( 1.0_LDP-X*(0.375_LDP-0.1_LDP*X) )/3.0_LDP
	ELSE
	  IF(X .LT. -700.0_LDP)THEN
	    N_DEXPONX_GUARD=N_DEXPONX_GUARD+1
	    Y=EXP(700.0_LDP)
	  ELSE IF(X .LT. 40)THEN
	    Y=EXP(-X)
	  ELSE
	    Y=0.0_LDP
	  END IF
	  d_EXPONX_dX=( Y-(1.0_LDP-Y)/X )/X
	END IF
C
	RETURN
	END
