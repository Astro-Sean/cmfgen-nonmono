!
! Subroutine to return the integral (in ANSWER) of a Gaussian line function
! of the form:
!
!           -HEIGHT*EXP(- 0.5D0*[ABS([X-POS]/SCALE)]**EXPONENT )/(A+B(X-X0)
!
! TOLERANCE is the absolute tolerance desired, such that the absolute error
!
!	 < PARAMS(1)*PARAMS(2)*TOL
!
! Routine uses Romberg integratio, with up to 8 levels of refinement.
! Routine should provide accuracies of better than 1.0E-04 (and much
! more generally) for
!
! EXPONENT > 0.5
!
	SUBROUTINE GAUSS_ROMB_V2(ANSWER,A,B,X0,POS,HEIGHT,SCALE,EXPONENT,TOLERANCE)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Altered 06-MAr-2023 : Fixed to use classical Gaussian with factor of 0.5 in argument of exponent.
! Altered 09-Aug-2022 : To get consistency inthe different routines changed to use Gauss.
! Altered 05-AUg-2022 : Now use rigorous definition of EW (not based on Ic at line
!                          center only). Call changed because need add.
!                          parameters.
! Altered 02-Apr-2008 : Use convention that EWs are +ve for absorption lines
! Created 28-Sep-2007
!
	REAL(KIND=LDP) ANSWER		!Returned EW
	REAL(KIND=LDP) A		!Continuum is defoned by A+B(X-X0)
	REAL(KIND=LDP) B
	REAL(KIND=LDP) X0
	REAL(KIND=LDP) POS		!Line center	
	REAL(KIND=LDP) HEIGHT		!
	REAL(KIND=LDP) SCALE		!Scale factor descibing eponential -- similar to sigma
	REAL(KIND=LDP) EXPONENT
	REAL(KIND=LDP) TOLERANCE	!Tolerence for final answer.
!
	REAL(KIND=LDP) X 		!X coordinate (as in original data)
	REAL(KIND=LDP) W		!Scaled and offset coodrinate: w=(X-POS)/SIGMA
!
! R(I,1) stores the ith trapazoidal integration.
!
	INTEGER, PARAMETER :: M=8
	REAL(KIND=LDP) R(M,M)
!
	REAL(KIND=LDP) H		!Step size
	REAL(KIND=LDP) HMAX		!Maximumstep size adopted
	REAL(KIND=LDP) RANGE		!Range of integration
!
	INTEGER J,K,L		!Loop indices
	INTEGER N		!Number of steps
	LOGICAL TESTING
!
	EXPONENT=ABS(EXPONENT)
	IF(EXPONENT .LT. 0.5)THEN
	  WRITE(6,*)'Error in GAUSS_ROMB'
	  WRITE(6,*)'Unable to integrate the modified Gaussian if'//
	1                      ' its exponent is < 0.5'
	  ANSWER=1000
	  RETURN
	END IF
!
	TESTING=.TRUE.
	RANGE=2.0*NINT(23.0**(1.0D0/EXPONENT))+1		!Set at function=1.0E-10
	HMAX=1  					!RANGE/10.0D0				!Maximum step size
	IF(TESTING)THEN
	  WRITE(6,*)'RANGE=',RANGE
	  WRITE(6,*)'HMAX=',HMAX
	END IF
!
	DO K=1,M
	  H=HMAX/2**(K-1)
	  W=-(L-1)*H
	  R(K,1)=0.0
	  N=NINT(RANGE/H)
	  IF(TESTING)THEN
	    WRITE(6,'(A,I8)')' Refinement number:',K
	    WRITE(6,'(A,I8)')' Number of steps is:',N
	    WRITE(6,'(A,ES16.6)')' Step size is:',H
	  END IF
	  DO L=-N+1,N
	    W=(L-1)*H
	    X=ABS(SCALE)*W+POS
	    R(K,1)=R(K,1)+EXP(-0.5D0*(ABS(W))**EXPONENT )/(A+B*(X-X0))
	  END DO
	  R(K,1)=R(K,1)*H                          !*1.0D0/SQRT(ACOS(-1.0D0))
	  IF(TESTING)THEN
	    WRITE(6,'(A,ES14.8)')'Current unscaled trapazoidal answer:',R(K,1)
	  END IF
	  DO L=2,K
	    J=K-L+1
	    R(J,L)=R(J+1,L-1)+(R(J+1,L-1)-R(J,L-1))/(4**(L-1)-1)
	  END DO
	  IF(K .NE. 1)THEN
	    IF(TESTING)THEN
	      WRITE(6,'(A,ES14.6)')'Current accuracy is',ABS(R(1,K)-R(1,K-1))
	    END IF
	    IF( ABS(R(1,K)-R(1,K-1)) .LT. TOLERANCE)EXIT
	  END IF	
	END DO
!
! We now use the convention that EWs are +ve for absorption lines.
!
	K=MIN(K,M)
	ANSWER=-R(1,K)*ABS(SCALE)*HEIGHT
	IF(TESTING)THEN
	  WRITE(6,'(A,ES14.6)')'Final unscaled integeral is:',R(1,K)
	  WRITE(6,'(A,ES14.6)')'          Final integral is:',ANSWER
	END IF
!
	RETURN
	END
