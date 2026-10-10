	SUBROUTINE WRITV(F,ND,A,LU)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Altered 18-Nov-2022 : Now TRIM A before writing out.
! Altered 28-May-1996 : IMPLICIT NONE installed.
! Altered 30-APR-1985 : Now writes 10 columns instead of five across a page.)
!
	INTEGER ND,LU
	REAL(KIND=LDP) F(ND)
!
! Local variables.
!
	INTEGER I,J
	CHARACTER*(*) A
!
	WRITE(LU,'(//,1X,A,/)')TRIM(A)
	DO I=1,ND,10
	  WRITE(LU,'(1X,1P,10E12.4)')(F(J),J=I,MIN(I+9,ND),1)
	END DO
!
	RETURN
	END
