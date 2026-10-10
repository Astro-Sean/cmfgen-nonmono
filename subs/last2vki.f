C
C Subroutine to compute the variation in intensity with opacity
C for rays with only one or two points.
C
	SUBROUTINE LAST2VKI(VK,SOURCE,CHI,DTAU,Z,TOR,NI)
	USE SET_KIND_MODULE
	IMPLICIT NONE
C
C Altered 24-May-1996 -EXP replaced bu EXP
C Altered 13-Apr-1988 - Bug fix.
C
	INTEGER NI
	REAL(KIND=LDP) VK(NI,NI),SOURCE(NI),CHI(NI),DTAU(NI),Z(NI)
	REAL(KIND=LDP) E1,E2,E3,DE1,DE2,DE3,TOR,T1,U2
C
C The incident intensity is assumed to be SOURCE(1)*(1.0-EXP(-TOR)).
C
	IF(TOR .GT. 0.01_LDP)THEN
	  T1=(1.0_LDP-EXP(-TOR))
	ELSE
	  T1=(1._LDP-TOR/2.0_LDP*(1.0_LDP-TOR/3.0_LDP*(1.0_LDP-TOR/4.0_LDP)))*TOR
	END IF
C
	IF(NI .EQ. 1)THEN
	  VK(1,1)=SOURCE(1)*(EXP(-TOR)*TOR-T1)/CHI(1)
	ELSE IF(NI .EQ. 2)THEN
C
C The D denotes derivative with respect to CHI.
C
	  E1=EXP(-DTAU(1))
	  DE1=-E1*(Z(1)-Z(2))*0.5_LDP
	  IF(DTAU(1) .LT. 1.0E-03_LDP)THEN
	    E2=DTAU(1)*0.5_LDP+DTAU(1)*DTAU(1)/6.0_LDP
	    DE2=(0.5_LDP+DTAU(1)/3.0_LDP)*(Z(1)-Z(2))*0.5_LDP
	    E3=DTAU(1)*0.5_LDP-DTAU(1)*DTAU(1)/3.0_LDP
	    DE3=(0.5_LDP-2.0_LDP*DTAU(1)/3.0_LDP)*(Z(1)-Z(2))*0.5_LDP
	  ELSE
	    E2=1.0_LDP-(1.0_LDP-E1)/DTAU(1)
	    DE2=(1.0_LDP-E1-DTAU(1)*E1)/DTAU(1)/DTAU(1)*(Z(1)-Z(2))*0.5_LDP
	    E3=(1.0_LDP-E1)/DTAU(1)-E1
	    DE3=((E1-1.0_LDP+E1*DTAU(1))+E1*DTAU(1)*DTAU(1))*
	1                    (Z(1)-Z(2))*0.5_LDP/DTAU(1)/DTAU(1)
	  END IF
C
C Note that U(2)=S1.(E1.T1+E3) + S2.E2
C and  then U(1)=0.5*[S1.(T1+E2) + U(2).E1 + S2.E3]
C
	  VK(2,2)=SOURCE(2)*(DE2-E2/CHI(2))+SOURCE(1)*(DE3+DE1*T1)
          VK(2,1)=SOURCE(1)*(DE1*T1+DE3+
	1              (E1*TOR*EXP(-TOR)-E3-E1*T1)/CHI(1))
	1              +SOURCE(2)*DE2
C
	  U2=SOURCE(1)*(T1*E1+E3)+SOURCE(2)*E2
	  VK(1,1)=0.5_LDP*(  VK(2,1)*E1
	1                 +SOURCE(1)*( DE2-(E2+T1-TOR*EXP(-TOR))/CHI(1) )
	1                 +SOURCE(2)*DE3+U2*DE1  )
	  VK(1,2)=0.5_LDP*(  VK(2,2)*E1+SOURCE(2)*(DE3-E3/CHI(2))
	1                 +SOURCE(1)*DE2+U2*DE1  )
	END IF
C
	RETURN
	END
