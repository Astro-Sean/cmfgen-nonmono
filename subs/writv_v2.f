!
	SUBROUTINE WRITV_V2(F,ND,NDEC,A,LU)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Altered 18-Nov-2022 : Now TRIM A before writing out.
! Altered 28-Apr-2000 : Bug fix: Length of FORM (at 15) was one character
!                                too short.
! Altered 28-May-1996 : IMPLICIT NONE installed.
! Altered 30-APR-1985 : Now writes 10 columns instead of five across a page.)
!
	INTEGER ND
	INTEGER LU
	INTEGER NDEC		!Number of decimal digits.
	REAL(KIND=LDP) F(ND)
	CHARACTER*(*) A
!
! Local variables.
!
	INTEGER I,J,L
	INTEGER N_PER_LINE
	INTEGER NX
	CHARACTER*16 FORM
!
	NX=NDEC+8
	N_PER_LINE=131/NX
!
	WRITE(FORM,'(A,I2.2,A,I2.2,A,I2.2,A)')
	1             '(1X,1P,',N_PER_LINE,'E',NX,'.',NDEC,')'
!
	L=1
	WRITE(LU,'(//,1X,A,/)')TRIM(A)
!
	IF(ND .GE. N_PER_LINE)THEN
	  DO I=1,ND+1-N_PER_LINE,N_PER_LINE
	    WRITE(LU,FORM)(F(J),J=I,I+N_PER_LINE-1,1)
	    L=I+N_PER_LINE
	  END DO
	END IF
	IF(L .LE. ND)WRITE(LU,FORM)(F(J),J=L,ND)
!
	RETURN
	END
