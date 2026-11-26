C
C Subroutine to compute the opacity due to electron scattering at
C ND depth points.
C
	SUBROUTINE ESOPAC(RKI,ED,ND)
	USE SET_KIND_MODULE
	IMPLICIT NONE
C
C Altered 24-May-1996 : Cleaned. (IMPLICT NONE inserted)
C
	INTEGER ND
	REAL(KIND=LDP) RKI(ND),ED(ND)
C
	RKI(1:ND)=6.65E-15_LDP*ED(1:ND)
C
	RETURN
	END
