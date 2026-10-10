!
! Routine designed to return the flux in a rectangular aperture which is
! located against a spherical atmosphere.
!
! When the aperture crosses the X-axis, we split the integration into
! 2 parts.
!
	SUBROUTINE INT_REC_AP(IP,P,YV,X_CENT,Y_CENT,WIDTH,LENGTH,NP,NCF,NINS)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! 05-Fev-2021: All arrays made REAL(KIND=LDP).
!
	INTEGER NCF
	INTEGER NP
	REAL(KIND=LDP) IP(NP,NCF)
	REAL(KIND=LDP) P(NP)
	REAL(KIND=LDP) YV(NCF)
!
	REAL(KIND=LDP) X_CENT
	REAL(KIND=LDP) Y_CENT
	REAL(KIND=LDP) WIDTH
	REAL(KIND=LDP) LENGTH
!
	INTEGER NINS
!
! Local variables
!
	REAL(KIND=LDP) LOC_X_CENT
	REAL(KIND=LDP) LOC_Y_CENT
	REAL(KIND=LDP) LOC_LENGTH
	REAL(KIND=LDP) T1
	LOGICAL ZERO_YV
!
! As we are dealing with a spherical star, the coordinates describing the center
! of the apperture can be taken to be positive.
!
	LOC_X_CENT=ABS(X_CENT)
	LOC_Y_CENT=ABS(Y_CENT)
!
	T1=1.0E-06_LDP*(P(2)-P(1))
	IF( LOC_X_CENT .LT. T1 .AND. LOC_Y_CENT .LT. T1)THEN
!
! Box is centerd, within rounding errors, on the star.
!
	  LOC_X_CENT=0.0_LDP
	  LOC_Y_CENT=0.0_LDP
	  ZERO_YV=.TRUE.
	  CALL DO_INT_REC_AP(IP,P,YV,LOC_X_CENT,LOC_Y_CENT,WIDTH,LENGTH,
	1          ZERO_YV,NP,NCF,NINS)
!
! Aperture crosses X axis.
!
	ELSE IF(LOC_Y_CENT .LT. 0.5_LDP*LENGTH)THEN
!
! Do box above X-axis.
!
	  LOC_LENGTH=LOC_Y_CENT+0.5_LDP*LENGTH
	  LOC_Y_CENT=0.5_LDP*LOC_LENGTH
	  ZERO_YV=.TRUE.
	  CALL DO_INT_REC_AP(IP,P,YV,LOC_X_CENT,LOC_Y_CENT,WIDTH,LOC_LENGTH,
	1          ZERO_YV,NP,NCF,NINS)
	  WRITE(6,*)'App 1',YV(1:2)
!
! Do box below X axis.
!
	  ZERO_YV=.FALSE.
	  LOC_LENGTH=0.5_LDP*LENGTH-ABS(Y_CENT)
	  LOC_Y_CENT=0.5_LDP*LOC_LENGTH
	  CALL DO_INT_REC_AP(IP,P,YV,LOC_X_CENT,LOC_Y_CENT,WIDTH,LOC_LENGTH,
	1          ZERO_YV,NP,NCF,NINS)
	  WRITE(6,*)'App 2',YV(1:2)
!
	ELSE
!
! Box must be above (or just on) Y axis.
!
	  ZERO_YV=.TRUE.
	  CALL DO_INT_REC_AP(IP,P,YV,LOC_X_CENT,LOC_Y_CENT,WIDTH,LENGTH,
	1          ZERO_YV,NP,NCF,NINS)
	END IF
!
	RETURN
	END

