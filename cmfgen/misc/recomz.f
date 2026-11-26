!
! Simple program to pvide rough estimates for recombinations
! rates for hydrogenic ions. The code returns the total recombination
! rate to levels n (input) and above. n must be > 1.
!
	PROGRAM RECOMZ
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Altered 07-Jul-2008: changed input LU from 6 to 5.
!
	INTEGER I
	REAL(KIND=LDP) Z,A,GAMMA,T,N,E,ALPHA
!
	GAMMA=0.5772_LDP
!
100	WRITE(6,*)'Input Z,T and N '
	READ(5,*)Z,T,I
	IF(Z .EQ. 0)STOP
	IF(I .LE. 1)THEN
	  WRITE(6,*)'Error, N must be greater than 1'
	  GOTO 100
	END IF
	N=I
C
	N=I-1
	A=15.7892_LDP*Z*Z/N/N/T
	E=log( 1.0_LDP+EXP(-GAMMA/(1.0_LDP+2.0_LDP*A))/A )
C
	ALPHA=1.032E-13_LDP*Z*Z/SQRT(T)*(  E+log(A)+GAMMA -
	1          A/N*( 1.0_LDP/3/N+(1.0_LDP-1.0_LDP/N-A/3/N)*E )  )
C
	WRITE(6,200)ALPHA
200	FORMAT(3X,'ALPHA=',1PE9.3)
C
	GOTO 100
C
	END
