!
! Subroutine to compute the temperature structure for a grey atmosphere
! with a Rosseland mean opacity of ROSSMEAN. When passed, ROSSMEAN must include
! the effects (if any) of clumping.
!
! At present, two different options are available.
!
! (a) Spherical atmosphere with out velocity terms.
! (b) If JGREY_WITH_V_TERMS =.TRUE. we use a spherical atmosphere,
!       and include the effect of the first order Doppler shift.
!
	SUBROUTINE COMP_GREY_V3(POPS,TGREY,TAU_ROSS,ROSSMEAN,PLANCKMEAN,LUER,NC,ND,NP,NT)
	USE SET_KIND_MODULE
	USE ANG_QW_MOD
	USE MOD_CMFGEN
	USE CONTROL_VARIABLE_MOD
	IMPLICIT NONE
!
! Altered 03-Apr-2009 : Inserted OUT_JH call in REL section.
! Altered 20-Jan-2008 : POPS & NT inserted into call. Changed to V2.
!                       Now call JGREY_HUB_DDT_V2
!
! Altered 30_Aug-2007 : Bug fix. Invalid lower boundary condition flux for
!                         MOM_JREL and for FCOMP_PP_V2i. Forgot factor of CHI.
!
! Created 21-Dec-2004
!
	INTEGER NC
	INTEGER ND
	INTEGER NP
	INTEGER NT
	INTEGER LUER
!
	REAL(KIND=LDP) POPS(NT,ND)
	REAL(KIND=LDP) TGREY(ND)
	REAL(KIND=LDP) TAU_ROSS(ND)
	REAL(KIND=LDP) ROSSMEAN(ND)
	REAL(KIND=LDP) PLANCKMEAN(ND)
!
	REAL(KIND=LDP) RJ(ND)
	REAL(KIND=LDP) DTAU(ND)
	REAL(KIND=LDP) Z(ND)
	REAL(KIND=LDP) CHI(ND)
	REAL(KIND=LDP) dCHIdR(ND)
!
	REAL(KIND=LDP) TA(ND)
	REAL(KIND=LDP) TB(ND)
	REAL(KIND=LDP) TC(ND)
	REAL(KIND=LDP) Q(ND)
	REAL(KIND=LDP) GAM(ND)
	REAL(KIND=LDP) GAMH(ND)
	REAL(KIND=LDP) SOB(ND)
	REAL(KIND=LDP) XM(ND)
	REAL(KIND=LDP) SOURCE(ND)
	REAL(KIND=LDP) FEDD(ND)
	REAL(KIND=LDP) H_ON_J(ND)
	REAL(KIND=LDP) RSQHNU(ND)
	REAL(KIND=LDP) dlnJdlnR(ND)
!
	REAL(KIND=LDP) T1,T2
	REAL(KIND=LDP) HBC_J
	REAL(KIND=LDP) HBC_CMF
	REAL(KIND=LDP) NBC_CMF
	REAL(KIND=LDP) INBC
	REAL(KIND=LDP) DBB
	REAL(KIND=LDP) FL
	REAL(KIND=LDP) PI
	REAL(KIND=LDP) HFLUX
!
	INTEGER I,J,K,L
!
	LOGICAL LST_DEPTH_ONLY
	LOGICAL FIRST_FREQ
	LOGICAL NEW_FREQ
!
! Constants for opacity etc.
!
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	REAL(KIND=LDP) CHIBF,CHIFF,HDKT,TWOHCSQ
!
	WRITE(6,*)'BEGIN COMP_GREY_V2',ND,NC,NP
	PI=4.0_LDP*ATAN(1.0_LDP)
	CHI(1:ND)=ROSSMEAN(1:ND)
!
! Compute the Rosseland optical depth scale. The ' ' in TORSCL indicates TYPE of atmosphere,
! and here is set to ' ' so that TORSCL assumes a power law dependence. If CHI(3) .LT. CHI(1),
! an 1/r^2 density dependence at the outer boundary is assumed.
!
        CALL TORSCL(TAU_ROSS,CHI,R,TB,TC,ND,METHOD,' ')
!
	IF(PLANE_PARALLEL_NO_V)THEN
           FEDD=0.3333_LDP
	   HBC_CMF=0.7_LDP; NBC_CMF=0.0_LDP; FL=1.0_LDP; INBC=0.1_LDP
	   NEW_FREQ=.TRUE.; FIRST_FREQ=.TRUE.
!
! Note HFLUX=LUM*Lsun/16/(PI*PI)/10**2 (10**2 for 1/R**2).
! DBB - dBdR = 3. Chi. L/16(piR)**2 and is used for the lower boundary
! diffusion approximation. Since we are dealing with a plane-parallel
! atmopshere, we divide HFLUX by R*^2.
!
	   HFLUX=3.826E+13_LDP*LUM/16.0_LDP/PI**2/R(ND)/R(ND)
           DBB=3.0_LDP*CHI(ND)*HFLUX
	   T1=1000.0_LDP
!
! Compute radial (vertical) optical depth increments.
!
	   CALL DERIVCHI(TB,CHI,R,ND,METHOD)
           CALL NORDTAU(DTAU,CHI,R,R,TB,ND)
!
!
! Compute the solution vector. Note that the units need to be
! eventually included. The following follows direcly from d2K/d2Tau=0.
!
	   DO WHILE(T1 .GT. 1.0E-08_LDP)
	     T2=FEDD(1)/HBC_CMF
	     RJ(1)=HFLUX/HBC_CMF
	     DO I=2,ND
	       T2=T2+DTAU(I-1)
	       RJ(I)=T2*HFLUX/FEDD(I)
	     END DO
!
	     SOURCE(1:ND)=RJ
	     CALL FCOMP_PP_V2(R,TC,GAMH,SOURCE,CHI,IPLUS,HBC_CMF,
	1               NBC_CMF,INBC,DBB,IC,THK_CONT,DIF,ND,NC,METHOD)
	     T1=0.0_LDP
	     DO I=1,ND
	       T1=MAX(ABS(FEDD(I)-GAMH(I)),T1)
	       FEDD(I)=GAMH(I)
	     END DO
	     NEW_FREQ=.FALSE.
	   END DO
!
	   OPEN(UNIT=7,FILE='GREY_CHK',STATUS='UNKNOWN')
	     WRITE(7,'(A)')
	     WRITE(7,'(A)')'Check for plane-parallel grey atmosphere calculation'
	     WRITE(7,'(A)')
	     WRITE(7,'(4X,A,5(11X,A))')'I','    R','    J','Ray J','    f','  Tau'
	     DO I=1,ND
	       WRITE(7,'(1X,I4,5(ES16.6))')I,R(I),RJ(I),TC(I),GAMH(I),TAU_ROSS(I)
	     END DO
	   CLOSE(UNIT=7)
!
	ELSE IF(JGREY_WITH_V_TERMS)THEN
!
! This routine will supersede the one above, and included to zero
! order the effect of the velocity field.
!
	   T2=1.0E-05_LDP		!Accuracy to converge f
	   CALL JGREY_WITH_FVT(RJ,SOB,CHI,R,V,SIGMA,
	1           P,AQW,HMIDQW,KQW,NMIDQW,
	1           LUM,METHOD,DIF,IC,T2,ND,NC,NP)
!
	ELSE IF(USE_J_REL)THEN
!
	  FEDD(1:ND)=0.3_LDP		!Initial guess
	  H_ON_J(1:ND)=0.0_LDP		!Initial guess
	  GAMH(1:ND)=0.0_LDP		!Old FEDD
	  XM(1:ND)=0.0_LDP		!As grey solution, not needed (ETA)
	  dlnJdlnR=0.0_LDP
	  NEW_FREQ=.TRUE.
	  WRITE(LUER,*)'Using MOM_JREL_GREY_V1 for grey solution'
!
! Note
!   HFLUX=LUM*Lsun/16/(PI*PI)/10**2/R**2 (10**2 for 1/R**2).
!   DBB = dBdR = 3.Chi.L/16(piR)**2
!   DBB is used for the lower boundary diffusion approximation.
!
	  HFLUX=3.826E+13_LDP*LUM/(4.0_LDP*PI*R(ND))**2
          DBB=3.0_LDP*HFLUX*CHI(ND)
	  T1=1.0_LDP
	  DO WHILE(T1 .GT. 1.0E-05_LDP)
            CALL MOM_JREL_GREY_V1(XM,CHI,CHI,V,SIGMA,R,
	1              H_ON_J,FEDD,dlnJdlnR,
	1              RJ,RSQHNU,HBC_CMF,INBC,
	1              DIF,DBB,IC,METHOD,
	1              L_TRUE,L_TRUE,NEW_FREQ,ND)
!
	    CALL FGREY_NOREL_V1(FEDD,H_ON_J,RJ,CHI,R,V,SIGMA,
	1              P,AQW,HMIDQW,KQW,LUM,IC,METHOD,
	1              HBC_CMF,INBC,DIF,ND,NC,NP)
	    T1=0.0_LDP
	    DO I=1,ND
	      T1=MAX(ABS(FEDD(I)-GAMH(I)),T1)
	      GAMH(I)=FEDD(I)
	    END DO
	    NEW_FREQ=.FALSE.
	    WRITE(LUER,*)'Current grey iteration accuracy is',T1
	  END DO
!
! Output gray fluxes and boundary Eddington factors for next model in time sequence.
!     T1  is used as a dummy argument for NU.
!      I  is used as a dummy argument for NCF. Neither value is used when the
!  'GREY' option is set. The logical variable L_TRUE is also not accessed.
!
	DO I=1,ND
	  TA(I)=RJ(I)*R(I)*R(I)
	END DO
	T2=HFLUX*R(ND)*R(ND)      !R^2.H at inner boundary.
	CALL OUT_JH(TA,RSQHNU,T2,HBC_CMF,T1,I,R,V,ND,L_TRUE,'GREY')
!
	ELSE IF(USE_DJDT_RTE)THEN
!
! TA is a temporary vector with the change in enthalpy.
!
	   WRITE(6,*)'Calling JGREY_HUB_DDT_V3'
	   T2=1.0E-06_LDP		!Accuracy to converge f
	   CALL JGREY_HUB_DDT_V3(RJ,SOB,CHI,PLANCKMEAN,R,V,SIGMA,POPS,
	1              P,AQW,HMIDQW,KQW,LUM,METHOD,DIF,IC,
	1              T2,INCL_DJDT_TERMS,TIME_SEQ_NO,ND,NC,NP,NT)
!
	ELSE
!
! Will use FEDD for F, GAM for NEWRJ, GAMH for NEWRK, and T2 for NEWHBC.
! Will use HBC_J for HBC. No need to modify JGREY, as outer boundary
! will always be optically thin.
!
	    DO I=1,ND
	      FEDD(I)=1.0_LDP/3.0_LDP
	    END DO
	    HBC_J=1.0_LDP
	    T1=1000.0_LDP
	    DO WHILE(T1 .GT. 1.0E-06_LDP)
	      CALL JGREY(TA,TB,TC,XM,DTAU,R,Z,P,RJ,
	1        GAM,GAMH,Q,FEDD,CHI,dCHIdR,
	1        AQW,KQW,LUM,HBC_J,T2,NC,ND,NP,METHOD)
	      T1=0.0_LDP
	      DO I=1,ND
	        T1=MAX(ABS(FEDD(I)-GAMH(I)),T1)
	        FEDD(I)=GAMH(I)
	      END DO
	      T1=MAX(ABS(HBC_J-T2),T1)
	      HBC_J=T2
!	      WRITE(LUER,'('' Maximum change in Grey F is '',1P,E11.4)')T1
	    END DO
	END IF
!
! Compute the temperature distribution, and the Rosseland optical depth scale.
! NB sigma=5.67E-05 and the factor of 1.0E-04 is to convert T from units of
! K to units of 10^4 K. The ' ' in TORSCL indicates TYPE of atmosphere,
! and here is set to ' ' so as TORSCL assumes a 1/r^2 density dependence
! at boundary.
!
	DO I=1,ND
	  TGREY(I)=((PI/5.67E-05_LDP*RJ(I))**0.25_LDP)*1.0E-04_LDP
	END DO
!
	RETURN
	END
