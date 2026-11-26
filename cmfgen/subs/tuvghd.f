C
C Routine to compute the arrays which characterize the
C Comoving-frame transfer equation.
C
C Aletered 25-MAR-82 - Sign of U variable changed.
C Altered 10-AUG-82 - Diffusion approximation and GC deleted.
C Altered 27-Apr-1989 - Cleaned - Implicit none installed.
C
	SUBROUTINE TUVGHD(TA,TB,TC,U,VB,VC,GB,H,Q,QH,
	1                   DTAU,DIFF,LS,NC,NI)
	USE SET_KIND_MODULE
	IMPLICIT NONE
C
	INTEGER LS,NC,NI
	REAL(KIND=LDP) TA(NI),TB(NI),TC(NI),VB(NI),VC(NI),GB(NI)
	REAL(KIND=LDP) H(NI),U(NI),Q(NI),QH(NI),DTAU(NI)
	LOGICAL DIFF
C
C Local variables.
C
	INTEGER I
C
	TA(1)=0.0_LDP
	VB(1)=0.0_LDP
	VC(1)=0.0_LDP
	GB(1)=-1.0_LDP/((1.0_LDP+QH(1))*DTAU(1))
	H(1)=QH(1)/(1.0_LDP+QH(1))
	TC(1)=1.0_LDP/DTAU(1)
	TB(1)=-TC(1)-1.0_LDP-Q(1)
	U(1)=-Q(1)
C
	DO I=2,NI-1
	  H(I)=QH(I)/(1.0_LDP+QH(I))
	  GB(I)=-1.0_LDP/((1.0_LDP+QH(I))*DTAU(I))
	END DO
C
	DO I=2,NI-1
	  TA(I)=-GB(I-1)
	  TC(I)=-GB(I)
	  VB(I)=-H(I-1)
	  VC(I)=H(I)
	  U(I)=-Q(I)*0.5_LDP*(DTAU(I)+DTAU(I-1))
	  TB(I)=-TA(I)-TC(I)-(1.0_LDP+Q(I))*0.5_LDP*(DTAU(I)+DTAU(I-1))
	END DO
C
	IF(LS .GT. NC)THEN
	  TA(NI)=GB(NI-1)
	  TB(NI)=(1.0_LDP+Q(NI))*0.5_LDP*DTAU(NI-1)-GB(NI-1)
	  U(NI)=Q(NI)*0.5_LDP*DTAU(NI-1)
	  VB(NI)=H(NI-1)
	ELSE IF(DIFF)THEN
	  TB(NI)=1.0_LDP/DTAU(NI-1)
	  TA(NI)=-TB(NI)
	  U(NI)=0.0_LDP
	  VB(NI)=0.0_LDP
	ELSE
	  TA(NI)=-1.0_LDP/DTAU(NI-1)
	  TB(NI)=-TA(NI)+1.0_LDP+Q(NI)
	  U(NI)=Q(NI)
	  VB(NI)=0.0_LDP
	END IF
	TC(NI)=0.0_LDP
	VC(NI)=0.0_LDP
	GB(NI)=0.0_LDP
C
	RETURN
	END
