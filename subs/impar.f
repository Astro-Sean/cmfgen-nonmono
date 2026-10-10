C
C Subroutine to compute the impact parameter values .
C
	SUBROUTINE IMPAR(P,R,RP,NC,ND,NP)
	USE SET_KIND_MODULE
	IMPLICIT NONE
C
C Altered 05-Dec-1996 - END DO used to terminate DO loops.
C Altered 24-May-1996 - IONE inserted
C Altered 17-Feb-1986 - Core rays distributed equally in mu rather than p.
C
	INTEGER NC,ND,NP,I
	REAL(KIND=LDP) P(NP),R(ND),RP,DELMU
	REAL(KIND=LDP), PARAMETER :: RZERO=0.0
	REAL(KIND=LDP), PARAMETER :: RONE=1.0
C
	DELMU=RONE/NC
	P(1)=RZERO
	DO I=2,NC
	  P(I)=RP*SQRT(RONE-(DELMU*(NC-I+1))**2)
	END DO
C
	DO I=NC+1,NP
	  P(I)=R(ND-I+NC+1)
	END DO
C
	RETURN
	END
