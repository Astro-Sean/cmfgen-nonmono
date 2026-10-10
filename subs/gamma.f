C
C Routine to compute frquency independent quantities which will
C be used to compte the frequency dependent quantities Qk,d and
C Qk,d+1/2  .
C
C Altered 27-Apr-1989 - Cleaned. Implicit none installed.
C
	SUBROUTINE GAMMA(GAM,GAMH,SIGMA,Z,R,V,ND,NI)
	USE SET_KIND_MODULE
	IMPLICIT NONE
C
	INTEGER ND,NI
	REAL(KIND=LDP) GAM(NI),GAMH(NI)
	REAL(KIND=LDP) SIGMA(ND),Z(NI),R(ND),V(ND)
C
C Local variables.
C
	INTEGER I
	REAL(KIND=LDP) MU
C
C Assume (1)	SIGMAd+1/2 = 0.5*( SIGMAd+1+SIGMAd )
C 	 (2)	Vd+1/2=0.5*( Vd + Vd+1 )
C Note that V is in km/s and SIGMA=(dlnV/dlnR-1.0)
C
	DO 10 I=1,NI-1
	  MU=Z(I)/R(I)
	  GAM(I)=3.33564E-06_LDP*V(I)/R(I)*(  1.0_LDP+SIGMA(I)*(MU**2)  )
	  MU=(Z(I)+Z(I+1))/(R(I)+R(I+1))
	  GAMH(I)=(V(I+1)+V(I))*3.33564E-06_LDP*(  1.0_LDP+0.5_LDP
	1    *(MU**2)*(SIGMA(I)+SIGMA(I+1))  )/(R(I)+R(I+1))
10	CONTINUE
	GAMH(NI)=0.0_LDP
	MU=Z(NI)/R(NI)
	GAM(NI)=3.33564E-06_LDP*V(NI)*( 1.0_LDP+SIGMA(NI)*(MU**2) )/R(NI)
C
	RETURN
	END
