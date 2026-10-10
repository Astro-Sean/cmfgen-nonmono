C
C Subroutine to compute JBAR and ZNET using the escape probability
C approximation. Six vectors are returned with new values.
C
C Input:
C
C    CHIL         = Total line opacity
C    ETAL         = Total line emissivity
C    CHIL_MAT(,J) = Line opacity for line J
C    ETAL_MAT(,J) = Line emissivity for line J
C    BBCOR(,J)    = Correction to emissivity to allow for frequncy difference.
C
C Output: Note that VB_2 and VC_2 are independent of the line J
C
C    JBAR  = Line mean intensity
C    BETA  = Sobolev Escape Probability.
C    VB_2  : dZnet(,J)/dCHIL(,K)=STOT*VB_2/SRCE(,J)
C    VC_2  : dZnet(,J)/dETAL(,K)=BB_COR(,J)*STOT*VC_2/SRCE(,J)
C
C Matrices (1 vector for each line)
C
C    ZNET(,J)    = Net rate line for line J
C    VB_SIM(,J)  = dZNET(,J)/dCHIL(,J)
C    VC_SIM(,J)  = dZNET(,J)/dETAL(,J)
C    BETAC(,J)   = Sobolev continuum term appearing in net rate.
C
	SUBROUTINE SOBJBAR_SIM(SOURCE,CHI,CHIL,ETAL,V,SIGMA,R,P,AQW,
	1                  JBAR,BETA,
	1                  ZNET,VB,VC,VB_2,VC_2,BETAC,
	1                  CHIL_MAT,ETAL_MAT,BBCOR,
	1                  FL,DIF,DBB,IC,THICK,
	1                  NLF,NC,NP,ND,NSIM,METHOD)
	USE SET_KIND_MODULE
	IMPLICIT NONE
C
C Altered 26-Oct-2014 : We now write ABS(dmuV/dz) so as to handle non-monotonic velocity laws.
C Altered 22-Aug-1997 : Adjusted to handle ND=NC+NP (i.e. NI=1 or 2 is
C                         now handled --- based on FQCOMP_IBC).
C Altered 28-Oct-1996 : Bug Fix: COS converted back to ACOS in TOR expression.
C Altered 26-May-1996 : Scratch blocked removed.
C                       Dynamic allocation for scratch vectors.
C                       IONE pased in call to THOMAS.
C                       Genrical calls to EXP, COS
C
C Altered 15-Jan-1992 VB_2, VC_2 inserted. Bug fixes.
C Created 13-Nov-1992 Allows for the simultaenous treatment of several
C                     lines. Based on SOBJBAR.
C
	INTEGER NLF,NC,NP,ND,NSIM
	REAL(KIND=LDP) SOURCE(ND),CHI(ND),CHIL(ND),ETAL(ND)
	REAL(KIND=LDP) V(ND),SIGMA(ND),R(ND),P(NP),AQW(ND,NP)
C
	REAL(KIND=LDP) JBAR(ND)
	REAL(KIND=LDP) BETA(ND)
C
	REAL(KIND=LDP) ZNET(ND,NSIM)
	REAL(KIND=LDP) VB(ND,NSIM)
	REAL(KIND=LDP) VC(ND,NSIM)
	REAL(KIND=LDP) VB_2(ND)
	REAL(KIND=LDP) VC_2(ND)
	REAL(KIND=LDP) BETAC(ND,NSIM)
C
	REAL(KIND=LDP) CHIL_MAT(ND,NSIM)
	REAL(KIND=LDP) ETAL_MAT(ND,NSIM)
	REAL(KIND=LDP) BBCOR(ND,NSIM)
C
	REAL(KIND=LDP) DBB,IC,FL
	CHARACTER*(*) METHOD
	LOGICAL DIF,THICK
C
	REAL(KIND=LDP) TA(ND),TB(ND),TC(ND),AV(ND),DTAU(ND),Z(ND)
	REAL(KIND=LDP) GAM(ND),GAMH(ND),dCHIdR(ND)
	REAL(KIND=LDP) EX_VEC(ND),dEX_VEC(ND)
	REAL(KIND=LDP) ZNET_BL(ND),BETAC_BL(ND)
	REAL(KIND=LDP) VB_BL(ND),VC_BL(ND)
	REAL(KIND=LDP) ST_ON_SL
C
	REAL(KIND=LDP) EXPONX,d_EXPONX_dX
C
C Local variables.
C
	INTEGER, PARAMETER :: IONE=1
	INTEGER I,J,LS,NI
	REAL(KIND=LDP) T1,T2,DBC,TOR
	REAL(KIND=LDP) IBOUND,E1,E2,E3  		!To handle NI=2
C
C Zero arrays which are incremented as we integrate over angle.
C Evaluate the SOBOLEV optical depth (GAMH) without angle factor.
C
	DO I=1,ND
	  JBAR(I)=0.0_LDP
	  BETA(I)=0.0_LDP
	  GAMH(I)=CHIL(I)*3.0E-10_LDP*R(I)/V(I)/FL    	!C/dex(15)/dex(5)
	  ZNET_BL(I)=0.0_LDP
	  VB_BL(I)=0.0_LDP
	  VC_BL(I)=0.0_LDP
	  BETAC_BL(I)=0.0_LDP
	END DO
C
	CALL DERIVCHI(dCHIdR,CHI,R,ND,METHOD)
C
C Enter loop to perform integration along each ray.
C
	DO LS=1,NP
	  NI=ND-(LS-NC-1)
	  IF(LS .LE. NC)NI=ND
C
C Initialize AV vector. GAM is the ANGLE DEPENDENT Sobolev optical depth.
C
	  AV(1:NI)=0.0_LDP
	  CALL ZALONGP(R,Z,P(LS),NI)
	  DO I=1,NI
C
C Dessart, Audit & Hillier (2015, arXiv:1503.05463) guard: for
C sigma < 0 the factor |1+sigma*mu^2| (mu^2=Z^2/R^2) can vanish along a
C ray, so the Sobolev optical depth is floored via
C MAX(|1+sigma*mu^2|, -0.2*(1+sigma)); 1D-4 bounds the sigma=-1 case.
C Inactive for sigma >= -1, so monotonic flows are unchanged.
C
	    T1=1.0_LDP+Z(I)*Z(I)/R(I)/R(I)*SIGMA(I)
	    IF(SIGMA(I) .LT. 0.0_LDP)THEN
	      T1=MAX(ABS(T1),-0.2_LDP*(1.0_LDP+SIGMA(I)),1.0D-4)
	    ELSE
	      T1=ABS(T1)
	    END IF
	    GAM(I)=GAMH(I)/T1
	  END DO
C
C SOURCE(1) is the boundary continuum source function.
C
	  IF(THICK)THEN
	    IF(P(LS) .GT. 0)THEN
	      TOR=CHI(1)*R(1)*R(1)*(1.570796_LDP-ACOS(P(LS)/R(1)))/P(LS)
	    ELSE
 	      TOR=CHI(1)*R(1)
	    END IF
	    IBOUND=SOURCE(1)*(1.0_LDP-EXP(-TOR))
	  END IF
C
C For NI=1 and 2 AV reprents the intensity variable U (the average
C intensity on the ray). For NI>2, AV only represnts u after the call
C to THOMAS.
C
	  IF(NI .EQ. 1)THEN
	    AV(1)=IBOUND
	  ELSE IF(NI .EQ. 2)THEN
	    CALL NORDTAU(DTAU,CHI,Z,R,dCHIdR,NI)
	    E1=EXP(-DTAU(1))
	    E2=1.0_LDP-(1.0_LDP-E1)/DTAU(1)
	    E3=(1.0_LDP-E1)/DTAU(1)-E1
	    IF(DTAU(1) .LT. 1.0E-03_LDP)THEN
	      E2=DTAU(1)*0.5_LDP+DTAU(1)*DTAU(1)/6.0_LDP
	      E3=DTAU(1)*0.5_LDP-DTAU(1)*DTAU(1)/3.0_LDP
	    END IF
	    AV(2)=IBOUND*E1+SOURCE(2)*E2+SOURCE(1)*E3
            AV(1)=0.5_LDP*(IBOUND+AV(2)*E1+SOURCE(1)*E2+SOURCE(2)*E3)
	  ELSE
	    IF(DIF .AND. LS .LE. NC)THEN
	      DBC=DBB*SQRT(R(ND)*R(ND)-P(LS)*P(LS))/R(ND)/CHI(ND)
	    END IF
	    CALL NORDTAU(DTAU,CHI,Z,R,dCHIdR,NI)
	    CALL TCOMPD(TA,TB,TC,DTAU,DIF,LS,NC,ND,NI)
	    CALL XVECD(DTAU,SOURCE,AV,DIF,DBC,IC,LS,NC,ND,NI)
	    AV(1)=-IBOUND
C
C Solve for the radiation field along ray for this frequency.
C
	    CALL THOMAS(TA,TB,TC,AV,NI,IONE)
	  END IF
C
C The function EXPONX is given by (1.0-EXP(-X))/X.
C The function d_EXPONX_dX is given by d[ (1.0-EXP(-X))/X ]/dX. These
C functions are called to allow for cancellation when X is small. Note
C GAM(I)/CHIL(I) = d_GAM(I)/d_CHIL(I) .
C
	  DO I=1,NI
            EX_VEC(I)=EXPONX(GAM(I))*AQW(I,LS)
            dEX_VEC(I)=d_EXPONX_dX(GAM(I))*AQW(I,LS)
	    BETA(I)=BETA(I)+EX_VEC(I)
	  END DO
C
	  DO I=1,NI
	    T2=AV(I)*CHIL(I)/ETAL(I)
	    T1=1.0_LDP-T2
C
	    ZNET_BL(I)=ZNET_BL(I)+T1*EX_VEC(I)
	    BETAC_BL(I)=BETAC_BL(I)+EX_VEC(I)*AV(I)
	    VB_BL(I)=VB_BL(I)+( T1*dEX_VEC(I)*GAM(I)/CHIL(I)-
	1                 AV(I)*EX_VEC(I)/ETAL(I))
	    VC_BL(I)=VC_BL(I)+T2*EX_VEC(I)/ETAL(I)
C
	  END DO		!End I loop
C
	END DO			!End LS Loop
C
C Now compute the net rates, and variations of the net rates. At this
C stage it is not clear whether this should be done as below, or done
C in loop above to improve accuracy.
C
	DO J=1,NSIM
	  DO I=1,ND
            ST_ON_SL=(ETAL(I)/ETAL_MAT(I,J))*
	1               (CHIL_MAT(I,J)/CHIL(I))/BBCOR(I,J)
	    ZNET(I,J)=ST_ON_SL*ZNET_BL(I)+(1.0_LDP-ST_ON_SL)
	    BETAC(I,J)=BETAC_BL(I)*CHIL_MAT(I,J)/ETAL_MAT(I,J)/BBCOR(I,J)
            VB(I,J)=ST_ON_SL*( VB_BL(I)+(ZNET_BL(I)-1.0_LDP)*
	1                      (1.0_LDP/CHIL_MAT(I,J)-1.0_LDP/CHIL(I)) )
            VC(I,J)=ST_ON_SL*( VC_BL(I)+(ZNET_BL(I)-1.0_LDP)*
	1                      (BBCOR(I,J)/ETAL(I)-1.0_LDP/ETAL_MAT(I,J)) )
	  END DO
	END DO
C
	DO I=1,ND
          VB_2(I)=VB_BL(I)-(ZNET_BL(I)-1.0_LDP)/CHIL(I)
          VC_2(I)=VC_BL(I)+(ZNET_BL(I)-1.0_LDP)/ETAL(I)
	END DO
C
C Compute the mean line intensity. This is not very accurate, and
C fails if CHIL(I,1) .EQ. 0
C
	DO I=1,ND
	  JBAR(I)=(1.0_LDP-ZNET_BL(I))*ETAL(I)/CHIL(I)
	END DO
C
	RETURN
	END
