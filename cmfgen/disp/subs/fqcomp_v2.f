!
! This routine solves for the mean intensity as a function of depth using the
! Feautrier Technique. A Schuster or diffusion approximation is used for the
! lower boundary condition. This routine may need to be in a loop so that
! the f values are iterated to convergence.
!
! Created 24-FEB-1987
! Altered 21-Sep-1990 --- Bug fix: IBOUND was not being zeroed when THK
!                         option was off. Only effected direct computation of
!                         HBCNEW.
! Altered 04-Oct-1990 --- INBCNEW was being evaluated incorrectly.
!                         Not used if DIFF is true.
! Altered 08-Jun-2010 --- Changed to RH algorithm so that very low optical depths
!                           can be treated. DIFF replaced by INNER_BND_METH.
!                           Changed to V2.
! Altered 15-Oct-2023 - Fixed bug -- DACOS was relaced by COS instead of ACOS.
!
	SUBROUTINE FQCOMP_V2(TA,TB,TC,XM,DTAU,R,Z,P,NEWRJ,NEWRK
	1  ,SOURCE,CHI,dCHIdr,AQW,AQW3,DBB,HBCNEW
	1  ,INBCNEW,IC,S1,THK,INNER_BND_METH,NC,ND,NP,METHOD)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
	INTEGER NC,ND,NP,I,NI,LS
	REAL(KIND=LDP) TA(ND),TB(ND),TC(ND),XM(ND),R(ND),Z(ND)
	REAL(KIND=LDP) NEWRJ(ND),NEWRK(ND),DTAU(ND),dCHIdR(ND)
	REAL(KIND=LDP) SOURCE(ND),CHI(ND),AQW(ND,NP),AQW3(ND,NP),P(NP)
!
	REAL(KIND=LDP) S1,DBB,DBC,IBOUND,TOR,HBCNEW,INBCNEW,IC,E1,E2,E3
	LOGICAL THK
	CHARACTER(LEN=*) METHOD
	CHARACTER(LEN=*) INNER_BND_METH
!
	CALL DERIVCHI(dCHIdr,CHI,R,ND,METHOD)
!
! Zero all parameters.
!
	CALL DP_ZERO(NEWRJ,ND)
	CALL DP_ZERO(NEWRK,ND)
	HBCNEW=0.0_LDP
	INBCNEW=0.0_LDP
!
! ENTER LOOP FOR EACH IMPACT PARAMETER P
!
	DO 2000 LS=1,NP
	  NI=ND-(LS-NC-1)
	  IF(LS .LE. NC)THEN
	    NI=ND
	    DBC=DBB*SQRT( (R(ND)-P(LS))*(R(ND)+P(LS)) )/R(ND)/CHI(ND)
	  END IF
!
	  IF(THK)THEN
	    IF(P(LS) .GT. 0)THEN
	      TOR=CHI(1)*R(1)*R(1)*(1.570796_LDP-ACOS(P(LS)/R(1)))/P(LS)
	    ELSE
	      TOR=CHI(1)*R(1)
	    END IF
	    IBOUND=S1*(1.0_LDP-EXP(-TOR))
	  ELSE
	    IBOUND=0.0_LDP
	  END IF
!
! Compute Z and the optical depth scale DTAU for this imapct parameter.
!
	  IF(NI .GT. 1)THEN
	    DO I=1,NI
	      Z(I)=SQRT( (R(I)-P(LS))*(R(I)+P(LS)) )
	    END DO
	    DO I=1,NI-1
	      DTAU(I)=0.5_LDP*(Z(I)-Z(I+1))*(CHI(I)+CHI(I+1)+(Z(I)-Z(I+1))
	1     *(dCHIdR(I+1)*Z(I+1)/R(I+1)-dCHIdR(I)*Z(I)/R(I))/6.0_LDP)
	    END DO
	  END IF
!
! Compute XM. Compute T ( a tridiagonal matrix) and store it as three vectors
! TA,TB and TC . This code is a combined version of XVECD and TCOMPD.
!
	  IF(NI .GT. 2)THEN
	    XM(1)=-IBOUND
	    TA(1)=0.0_LDP
	    TC(1)=1.0_LDP/DTAU(1)
	    TB(1)=1.0_LDP
	    DO I=2,NI-1
	      TA(I)=TC(I-1)
	      TC(I)=1.0_LDP/DTAU(I)
	      TB(I)=0.5_LDP*(DTAU(I-1)+DTAU(I))
	      XM(I)=-SOURCE(I)*(DTAU(I-1)+DTAU(I))*0.5
	    END DO
!
	    IF(LS .GT. NC)THEN
	      TA(NI)=-TC(NI-1)
	      TB(NI)=-0.5_LDP*DTAU(NI-1)
	      XM(NI)=0.5_LDP*DTAU(NI-1)*SOURCE(NI)
	    ELSE IF(INNER_BND_METH .EQ. 'DIFFUSION')THEN
	      TA(NI)=-TC(NI-1)
	      TB(NI)=0.0_LDP
	      XM(NI)=DBC
	    ELSE IF(INNER_BND_METH .EQ. 'SCHUSTER')THEN
              TA(NI)=-TC(NI-1)
              TB(NI)=1-TA(NI)
              XM(NI)=IC
	   ELSE
!
! HOLLOW core or ZERO_FLUX. For continuum, this is the same as for LS >
! NC. Done as separate option for clarity.
!
	      TA(NI)=-TC(NI-1)
	      TB(NI)=-DTAU(NI-1)/2.0_LDP
	      XM(NI)=0.5_LDP*DTAU(NI-1)*SOURCE(NI)
	    END IF
	    TC(NI)=0.0
!
! Solve the tridiagonal system of equations.
!
	    CALL THOMAS_RH(TA,TB,TC,XM,NI,1)
	    IF(MINVAL(XM) .LT. 0)THEN
	      WRITE(6,*)'Negative intensities for LS=',LS,DBC,IBOUND
	      CALL WRITE_VEC(XM,NI,'U in FQCOMP_V2',6)
	      CALL WRITE_VEC(SOURCE,ND,'SOURCE in FQCOMP_V2',6)
	      CALL WRITE_VEC(CHI,NI,'CHI in FQCOMP_V2',6)
	      CALL WRITE_VEC(DTAU,NI-1,'DTAU in FQCOMP_V2',6)
	      STOP
	    END IF
!
	  ELSE IF(NI .EQ. 1)THEN
	    XM(1)=IBOUND
	  ELSE IF(NI .EQ. 2)THEN
	    E1=EXP(-DTAU(1))
	    E2=1.0_LDP-(1.0_LDP-E1)/DTAU(1)
	    E3=(1.0_LDP-E1)/DTAU(1)-E1
	    IF(DTAU(1) .LT. 1.0E-03_LDP)THEN
	      E2=DTAU(1)*0.5_LDP+DTAU(1)*DTAU(1)/6.0_LDP
	      E3=DTAU(1)*0.5_LDP-DTAU(1)*DTAU(1)/3.0_LDP
	    END IF
!
	    XM(2)=IBOUND*E1+SOURCE(2)*E2+SOURCE(1)*E3
            XM(1)=0.5_LDP*(IBOUND+XM(2)*E1+SOURCE(1)*E2+SOURCE(2)*E3)
	  END IF
!
! UPDATE THE FA AND FB MATRICES . (SEE NOTES)
!
	  DO I=1,NI
	    NEWRJ(I)=NEWRJ(I)+AQW(I,LS)*XM(I)
	    NEWRK(I)=NEWRK(I)+AQW3(I,LS)*XM(I)
	  END DO
!
	  HBCNEW=HBCNEW+AQW(1,LS)*(XM(1)-IBOUND)*Z(1)/R(1)
	  IF(NI .EQ. ND)THEN
	    INBCNEW=INBCNEW + AQW(ND,LS)*
	1    (XM(ND)-(XM(ND)-XM(ND-1))/DTAU(ND-1))*Z(ND)/R(ND)
	  END IF
!
2000	CONTINUE
!
! Compute the factor for the outer boundary condition.
!
	HBCNEW=HBCNEW/NEWRJ(1)
	INBCNEW=INBCNEW/(2.0_LDP*NEWRJ(ND)-IC)
!
! Compute the new Feautrier factors.
!
	DO I=1,ND
	  NEWRK(I)=NEWRK(I)/NEWRJ(I)
	END DO
!
	DO I=1,ND
	  IF(NEWRJ(I) .LE. 0 .OR. NEWRK(I) .LE. 0)THEN
	   WRITE(6,*)'Error in FQCOMP_V2 -- J or K(f) is -ve'
	   WRITE(6,*)'INNER_BND_METH=',INNER_BND_METH
	   CALL WRITE_VEC(NEWRJ,ND,'NEWRJ in FQCOMP_V2',6)
	   CALL WRITE_VEC(NEWRK,ND,'F in FQCOMP_V2',6)
	   CALL WRITE_VEC(SOURCE,ND,'SOURCE in FQCOMP_V2',6)
	   CALL WRITE_VEC(CHI,ND,'CHI in FQCOMP_V2',6)
	   STOP
	 END IF
	END DO
	WRITE(6,*)'Checked for negative J values in FQCOMP_V2'
!
! Compute the Q factors from F. The Q values are stored in NEWRJ.
! We now compute Q in JFEAUNEW so that we nolonger need to save Q (only f).
! We still compute Q here because
!    (i) Q is passed (as NEWRJ vector) and hence is corrupted by this
!        routine. We still need Q, however, for the varaition
!        routines.
!   (ii) Q is now completely consistent with the updated f value.
!
	CALL QFROMF(NEWRK,NEWRJ,R,TA,TB,ND)	!TA work vector
!
	RETURN
	END
