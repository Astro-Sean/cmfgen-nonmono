!
! Program to evaluate the error assoiciated with a multiple pseudo-Gaussian fit
! to a set of data points passed by the module GAUSS_FIT_DATA.
!
! The function has the form
!
!   Y =P1 + P2*(X-X(1))+ P5*EXP( -((X-P3)/P4)^P6) + ...
!
! This routine must be kept compatible with GAUSS_FIT_FUNC
!
! Altered 06-MAr-2023 : Fixed to use classical Gaussian with factor of 0.5 in argument of exponent.
!                       Normalize error calculation by continuum.
! Altered 09-Aug-2022 : To get consistency inthe different routines changed to use Gauss.
! Created 05-Oct-2007
!
	SUBROUTINE GAUSS_FIT_ER(PARAMS)
	USE SET_KIND_MODULE
	USE GAUSS_FIT_DATA
	IMPLICIT NONE
	REAL(KIND=LDP) PARAMS(NG_PAR)
!
	REAL(KIND=LDP) T1
	REAL(KIND=LDP) SUM
	INTEGER I,J,K
!
! Make sure YFIT is up to date.
!
	DO J=1,NG_DATA
	  SUM=PARAMS(1)+PARAMS(2)*(X_GAUSS(J)-X_GAUSS(1))
	  YCONT_FIT(J)=SUM
	  DO K=3,NG_PAR,4
	    SUM=SUM+PARAMS(K+2)*EXP(-0.5*(ABS((X_GAUSS(J)-PARAMS(K))/PARAMS(K+1)))**PARAMS(K+3))
	  END DO
	  YFIT(J)=SUM
	END DO
!
! Conservative error estimate based on fit quality. We use abs value, with 1.0-03 of
! line center. Possible influence of multiple lines is ignored.
!
	IF(ALLOCATED(EW_ERROR))DEALLOCATE(EW_ERROR,ALT_ERROR,MIN_ERROR)
	ALLOCATE (EW_ERROR(NUM_GAUSS),ALT_ERROR(NUM_GAUSS),MIN_ERROR(NUM_GAUSS))
	EW_ERROR(:)=0.0D0; ALT_ERROR(:)=0.0D0; MIN_ERROR(:)=0.0D0
	DO J=1,NG_DATA
	  DO I=1,NUM_GAUSS
	    K=3+4*(I-1)
	    T1=EXP(-(ABS((X_GAUSS(J)-PARAMS(K))/PARAMS(K+1)))**PARAMS(K+3))
	    IF(T1 .GT. 1.0D-03)THEN
	      EW_ERROR(I)=EW_ERROR(I)+(X_GAUSS(MIN(J+1,NG_DATA))-X_GAUSS(MAX(1,J-1)))*
	1                ABS(Y_GAUSS(J)-YFIT(J))/YCONT_FIT(J)
	    END IF
	    IF(ABS((X_GAUSS(J)-PARAMS(K))/PARAMS(K+1)) .LT. 4)THEN
	      ALT_ERROR(I)=ALT_ERROR(I)+(X_GAUSS(MIN(J+1,NG_DATA))-X_GAUSS(MAX(1,J-1)))*
	1                ABS(Y_GAUSS(J)-YFIT(J))/YCONT_FIT(J)
	      MIN_ERROR(I)=MIN_ERROR(I)+(X_GAUSS(MIN(J+1,NG_DATA))-X_GAUSS(MAX(1,J-1)))*
	1                (Y_GAUSS(J)-YFIT(J))/YCONT_FIT(J)
	    END IF
	  END DO
	END DO
	EW_ERROR=0.5D0*EW_ERROR
	ALT_ERROR=0.5D0*ALT_ERROR
	MIN_ERROR=0.5D0*ABS(MIN_ERROR)
!
	RETURN
	END
