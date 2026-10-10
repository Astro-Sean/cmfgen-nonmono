	FUNCTION LAMBERT_WM_FUN(Z,LNZ)
	USE SET_KIND_MODULE
	IMPLICIT NONE
	REAL(KIND=LDP) LAMBERT_WM_FUN
!
! Function to compute the lambet function W_{-1}(Z).
!
	REAL(KIND=LDP) Z
	REAL(KIND=LDP) LNZ
	REAL(KIND=LDP) W
!
	REAL(KIND=LDP) EN
	REAL(KIND=LDP) QN
	REAL(KIND=LDP) ZN
	REAL(KIND=LDP) T1
	REAL(KIND=LDP) WOLD
!
	INTEGER J
!
! Make an initial guess -- these are god enough for the
! iteration procdures to converge in a few iterations.
!
	T1=-EXP(-1.0_LDP)
	IF(LNZ .LT. -100)THEN
	  W=LNZ
	  W=W-LOG(-W)*(1.0_LDP-1.0_LDP/W)
	ELSE IF(Z .LT. T1)THEN
	  WRITE(6,*)'Invalid value of z in LAMBERT_WM_FUN. Z=',Z
	  WRITE(6,*)'Argument must be greater (or =) to -EXP(-1.0D0)'
	  STOP
	ELSE IF(Z .EQ. T1)THEN
	  W=-1.0_LDP
	ELSE
	  IF(Z .GT. -0.15_LDP)THEN
	    W=LOG(-Z)
	    W=W-LOG(-W)*(1.0_LDP-1.0_LDP/W)
	  ELSE
	    W=-1.0_LDP-9.1_LDP*(EXP(-1.0_LDP)+Z)
	  END IF
!
	  WOLD=1
	  DO J=1,10
	    ZN=LNZ-LOG(-W)-W
	    QN=2*(1+W)*(1+W+2*ZN/3.0_LDP)
	    EN=ZN/(1+W)*(QN-ZN)/(QN-2*ZN)
	    W=W*(1.0_LDP+EN)
	    IF(W .GT. 0.0_LDP)W=-W
	    IF( ABS( (W-WOLD)/WOLD ) .LT. 1.0E-08_LDP)EXIT
	    WOLD=W
	  END DO
	END IF
!
	LAMBERT_WM_FUN=W
!
	RETURN
	END
