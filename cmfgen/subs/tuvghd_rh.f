C
C Routine to compute the arrays which characterize the
C Comoving-frame transfer equation.
C
C NB: Equations are
C
C       TA(i).Y(i-1) - [ TA(i)+TB(i)+TC(i) ].Y(i) + TC(i).Y(i) = RHS(I) + ...
C
	SUBROUTINE TUVGHD_RH(TA,TB,TC,U,VB,VC,GB,H,RHS,
	1                   Q,QH,DTAU,SOURCE,DIFF,DBC,IC,LS,NC,NI)
	USE SET_KIND_MODULE
	IMPLICIT NONE
C
C Created 14-Feb-1994 - Based on TUVGHD and XVECD
C
	INTEGER LS,NC,NI
	REAL(KIND=LDP) TA(NI),TB(NI),TC(NI),VB(NI),VC(NI),GB(NI),RHS(NI)
	REAL(KIND=LDP) H(NI),U(NI),Q(NI),QH(NI),DTAU(NI),SOURCE(NI)
	REAL(KIND=LDP) DBC,IC
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
	TB(1)=1.0_LDP+Q(1)
	U(1)=-Q(1)
	RHS(1)=0.0_LDP
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
	  U(I)=-0.5_LDP*Q(I)*(DTAU(I)+DTAU(I-1))
	  TB(I)=0.5_LDP*(1.0_LDP+Q(I))*(DTAU(I)+DTAU(I-1))
	  RHS(I)=-0.5_LDP*SOURCE(I)*(DTAU(I-1)+DTAU(I))
	END DO
C
	IF(LS .GT. NC)THEN
	  TA(NI)=GB(NI-1)
	  TB(NI)=-0.5_LDP*(1.0_LDP+Q(NI))*DTAU(NI-1)
	  U(NI)=0.5_LDP*Q(NI)*DTAU(NI-1)
	  VB(NI)=H(NI-1)
	  RHS(NI)=0.5_LDP*DTAU(NI-1)*SOURCE(NI)
	ELSE IF(DIFF)THEN
	  TB(NI)=0.0_LDP
	  TA(NI)=-1.0_LDP/DTAU(NI-1)
	  U(NI)=0.0_LDP
	  VB(NI)=0.0_LDP
	  RHS(NI)=DBC
	ELSE
	  TA(NI)=-1.0_LDP/DTAU(NI-1)
	  TB(NI)=-1.0_LDP-Q(NI)
	  U(NI)=Q(NI)
	  VB(NI)=0.0_LDP
	  RHS(NI)=IC
	END IF
	TC(NI)=0.0_LDP
	VC(NI)=0.0_LDP
	GB(NI)=0.0_LDP
C
	RETURN
	END
