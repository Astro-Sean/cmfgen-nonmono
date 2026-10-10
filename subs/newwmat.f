C
C Subroutine to compute the W matrix (see notes).
C
	SUBROUTINE NEWWMAT(W,DTAU,THETA,TOR,LS,NC,NI)
	USE SET_KIND_MODULE
	IMPLICIT NONE
C
C Altered 24-May-1996 - CALL to DP_ZERO removed.
C                       Calls to EXP replaced by generic calls.
C Altered  9-Dec-1986 - Thick boundary condition installed.
C
	INTEGER LS,NC,NI,I
	REAL(KIND=LDP) W(NI,NI),DTAU(NI),THETA(NI),TOR
C
	W(:,:)=0.0_LDP
C
C Note the general non thick case, TOR is 0, hence W(1,1)=0 . The thick
C case ocurrs when TOR is non zero.
C
	 IF(TOR .GT. 0.01_LDP)THEN
	  W(1,1)=THETA(1)*(1.0_LDP-EXP(-TOR))
	ELSE
	  W(1,1)=THETA(1)*
	1   (1.0_LDP-TOR/2.0_LDP*(1.0_LDP-TOR/3.0_LDP*(1.0_LDP-TOR/4.0_LDP)))*TOR
	END IF
C
	DO 20 I=2,NI-1
	  W(I,I)=THETA(I)*(DTAU(I-1)+DTAU(I))*0.5_LDP
20	CONTINUE
	IF(LS .GT. NC)W(NI,NI)=-0.5_LDP*DTAU(NI-1)*THETA(NI)
C
	RETURN
	END
