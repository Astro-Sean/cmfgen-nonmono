!
! Routine to compute the equivalent width of a line and the
! radiative force of a line using the SOBOLEV approximation.
! The radiative force is expressed as a multiple of the
! force on a free electron.
!
	SUBROUTINE SOBEW_GRAD_V3(SOURCE,CHI,ESEC,CHIL,ETAL,V,SIGMA,R,P,
	1                 FORCE_MULT,RLUM,
	1                 AQW,HAQW,LINE_FLUX,EW,EMISS_EW,CONT_FLUX,
	1                 FL,INNER_BND_METH,DBB,IC,THICK,DIE,NC,NP,ND,METHOD)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Altered 01-Mar-2019 : Changed to V3; added EMISS_EW to call
! Altered 26-Oct-2014 : We now write ABS(dmuV/dz) so as to handle non-monotonic velocity laws.
! Altered 07-Jun-2010: Changed to THOMAS_RH to handle later stages of SN.
!                        Replaced DIF variable by INNER_BND_METH.
! Created 24-Oct-2001: Based on SOBEW
!
	INTEGER NC,NP,ND
	REAL(KIND=LDP) SOURCE(ND),CHI(ND),ESEC(ND),CHIL(ND),ETAL(ND)
	REAL(KIND=LDP) FORCE_MULT(ND)
	REAL(KIND=LDP) RLUM
	REAL(KIND=LDP) V(ND),SIGMA(ND),R(ND),P(NP),AQW(ND,NP),HAQW(ND,NP)
	REAL(KIND=LDP) DBB,IC,FL,CONT_FLUX,LINE_FLUX(ND)
	REAL(KIND=LDP) EW, EMISS_EW
	CHARACTER*(*) METHOD
	CHARACTER*(*) INNER_BND_METH
	LOGICAL THICK,DIE
!
! Use dynamic allocation for required vectors.
!
	REAL(KIND=LDP) TA(ND),TB(ND),TC(ND),AV(ND),DTAU(ND),Z(ND)
	REAL(KIND=LDP) GAM(ND),GAMH(ND),dCHIdR(ND),NOES(ND)
	REAL(KIND=LDP) GLINE(ND),EMISS_LINE_FLUX(ND)
!
	REAL(KIND=LDP) EXPONX
!
! Local variables.
!
	INTEGER, PARAMETER :: IONE=1
!
	INTEGER I,LS,NI
	REAL(KIND=LDP) T1,T2,T3,T4,TOR
	REAL(KIND=LDP) E1,E2,E3
	REAL(KIND=LDP) IBOUND
	REAL(KIND=LDP) JFLUX
!
	IF( INNER_BND_METH .NE. 'DIFFUSION' .AND.
	1   INNER_BND_METH .NE. 'SCHUSTER' .AND.
	1   INNER_BND_METH .NE. 'HOLLOW' .AND.
	1   INNER_BND_METH .NE. 'ZERO_FLUX' )THEN
	  WRITE(6,*)'Error in SOBEW_GRAD_V2'
	  WRITE(6,*)'Unrecognized INNER_BND_METH: ',TRIM(INNER_BND_METH)
	  STOP
	END IF
!
! Zero arrays which are incremented as we integrate over angle.
! Evaluate the SOBOLEV optical depth without angle factor (GAMH).
! Evaluate the thermal opacity.
!
	DO I=1,ND
	  EMISS_LINE_FLUX(I)=0.0_LDP
	  LINE_FLUX(I)=0.0_LDP
	  NOES(I)=CHI(I)-ESEC(I)
	  GAMH(I)=CHIL(I)*3.0E-10_LDP*R(I)/V(I)/FL    	!C/dex(15)/dex(5)
	END DO
	EW=0.0_LDP
	CONT_FLUX=0.0_LDP
	GLINE(1:ND)=0.0_LDP
	JFLUX=0.0_LDP
!
	IF(DIE)THEN
	  DO I=1,ND
	    GAMH(I)=0.0_LDP
	  END DO
	END IF
!
	CALL DERIVCHI(dCHIdR,CHI,R,ND,METHOD)
!
! Enter loop to perform integration along each ray.
!
	DO LS=1,NP
	  NI=ND-(LS-NC-1)
	  IF(LS .LE. NC)NI=ND
!
! Zero AV vector.
!
	  AV(1:NI)=0.0
!
! SOURCE(1) is the boundary continuum source function. NB: If THICK it is
!           not obvious what you mean by the FLUX going to the observer.
!           Is it Iplus or Iplus-Ibound? Hopefully the outer boundary is
!           at sufficiently large radii that it makes no difference.
!
	  IF(THICK)THEN
	    IF(P(LS) .GT. 0)THEN
	      TOR=CHI(1)*R(1)*R(1)*(1.570796_LDP-ACOS(P(LS)/R(1)))/P(LS)
	    ELSE
	      TOR=CHI(1)*R(1)
	    END IF
	    IBOUND=SOURCE(1)*(1.0_LDP-EXP(-TOR))
	  ELSE
	    IBOUND=0
	  END IF
!
	  IF(NI .EQ. 1)THEN
	    AV(1)=IBOUND
	    TB(1)=0.0_LDP			!Flux at grid point.
	  ELSE IF(NI .EQ. 2)THEN
	    CALL ZALONGP(R,Z,P(LS),NI)
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
            TB(2)=0.0_LDP
	    TB(1)=AV(1)-IBOUND
	  ELSE
!
	    CALL ZALONGP(R,Z,P(LS),NI)
!
! NB: GAM is the angle dependent SOBOLEV optical depth.
!
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
!
! Compute optical depth increments.
!
	    CALL NORDTAU(DTAU,CHI,Z,R,dCHIdR,NI)
!
! Compute tridiagonal components.
! NB: TB(I)[tri matrix] = -[TB(I)+TA(I)+TC(I)]
! This section was modified from TCOMPD.
!
	    TA(1)=0.0_LDP
	    TB(1)=1.0_LDP
	    AV(1)=-IBOUND
	    DO I=1,NI-1
	      TC(I)=1.0_LDP/DTAU(I)
	    END DO
	    DO I=2,NI-1
	      TA(I)=TC(I-1)
	      TB(I)=0.5_LDP*(DTAU(I-1)+DTAU(I))
	      AV(I)=-SOURCE(I)*(DTAU(I-1)+DTAU(I))*0.5_LDP
	    END DO
!
	    IF(LS .GT. NC)THEN
	      TA(NI)=-TC(NI-1)
	      TB(NI)=-DTAU(NI-1)/2.0_LDP
	      AV(NI)=0.5_LDP*DTAU(NI-1)*SOURCE(NI)
	    ELSE IF(INNER_BND_METH .EQ. 'DIFFUSION')THEN
	      TA(NI)=-TC(NI-1)
	      TB(NI)=0.0_LDP
	      AV(NI)=DBB*SQRT(R(ND)*R(ND)-P(LS)*P(LS))/R(ND)/CHI(ND)
	    ELSE IF(INNER_BND_METH .EQ. 'SCHUSTER')THEN
	      TA(NI)=-TC(NI-1)
	      TB(NI)=-1.0_LDP
	      AV(NI)=IC
	    ELSE
	      TA(NI)=-TC(NI-1)
	      TB(NI)=-0.5_LDP*DTAU(NI-1)
	      AV(NI)=0.5_LDP*DTAU(NI-1)*SOURCE(NI)
	    END IF
	    TC(NI)=0.0_LDP
!
!
! Solve for the radiation field along ray for this frequency.
!
	    CALL THOMAS_RH(TA,TB,TC,AV,NI,IONE)
!
! Compute the midpoint flux in order that we can compute the
! line EW.  TB and TC are not required further.
!
	    DO I=1,NI-1
	      TC(I)=(AV(I+1)-AV(I))/DTAU(I)
	    END DO
!
! Interpolate the flux onto the radius grid. We handle the end points by
! using the bondary conditions.
!
	    DO I=2,NI-1
	      TB(I)=(DTAU(I-1)*TC(I)+DTAU(I)*TC(I-1))/(DTAU(I-1)+DTAU(I))
	    END DO
	    TB(1)=AV(1)-IBOUND
	    IF(LS .LE. NC)THEN
	      IF(INNER_BND_METH .EQ. 'DIFFUSION')THEN
	        TB(NI)=DBB*SQRT(R(ND)*R(ND)-P(LS)*P(LS))/R(ND)/CHI(ND)
	       ELSE IF(INNER_BND_METH .EQ. 'SCHUSTER')THEN
	        TB(NI)=IC-AV(NI)
	      ELSE
	        TB(NI)=0.0_LDP            !Hollow core and zero_flux options.
	      END IF
	    ELSE
	      TB(NI)=0.0_LDP		!Symmetry
	    END IF
!
	  END IF			!NI > 2
!
! Update the continuum flux.
!
	  CONT_FLUX=CONT_FLUX+TB(1)*HAQW(1,LS)
	  JFLUX=JFLUX+AV(1)*AQW(1,LS)
!
! Compute the thermal optical depth scale. Can use DTAU since
! this is not required anymore.
!
	  DTAU(1)=NOES(1)*R(1)
	  DO I=2,NI
	    DTAU(I)=DTAU(I-1)+0.5_LDP*(NOES(I-1)+NOES(I))*(Z(I-1)-Z(I))
	  END DO
!
! The function EXPONX is given by (1.0-EXP(-X))/X.
!
	  DO I=1,NI
	    T2=EXPONX(GAM(I))
!
! The radiative force due to a single line is
!
!      (4pi chil)/(rohc) Int v(mu) Beta(mu) mu dmu
!
! The following just evaluates the integral part of this expression.
!
	    GLINE(I)=GLINE(I)+TB(I)*T2*HAQW(I,LS)
!
! NB: T3 is the optical depth from the oberver to the point I on the NEAR side.
!        of the atmosphere.
!
	    IF(DTAU(I) .GT. 50)THEN
	      T3=0.0_LDP
	    ELSE
	      T3=EXP(-DTAU(I))
	    END IF
!
! NB: T4 is the optical depth from the oberver to the point I on the READ side
!        of the atmosphere.
!
	    T4=DTAU(NI)+DTAU(NI)-DTAU(I)
	    IF(LS .LE. NC .OR. T4 .GT. 50)THEN
	      T4=0.0_LDP
	    ELSE
	      T4=EXP(-T4)
	    END IF
!
	    LINE_FLUX(I)=LINE_FLUX(I)
	1         +0.5_LDP*T2*AQW(I,LS)
	1         *(  (T3+T4)-( (AV(I)+TB(I))*T3+(AV(I)-TB(I))*T4 )
	1         *CHIL(I)/ETAL(I)  )
!
	    EMISS_LINE_FLUX(I)=EMISS_LINE_FLUX(I)+0.5_LDP*T2*AQW(I,LS)*(T3+T4)
!
	  END DO
!
	END DO			!LS
!
! Evaluate the continuum intensity (Jy (d=1kpc)) and compute line
! emission function such that integral gives the equivalent width
! in Angstroms.
!
	CONT_FLUX=CONT_FLUX*R(1)*R(1)*13.19868_LDP
	T1=13.19868_LDP*(2.997924E-12_LDP/FL/FL)/CONT_FLUX
	DO I=1,ND
	  LINE_FLUX(I)=T1*LINE_FLUX(I)*ETAL(I)*( R(I)**3 )
	  EMISS_LINE_FLUX(I)=T1*EMISS_LINE_FLUX(I)*ETAL(I)*( R(I)**3 )
	END DO
!
! Compute the line flux, and then evaluate the continuum intensity
! (Jy (d=1kpc)) and the line equivalent width (Angstroms).
!
	EW=0.0_LDP
	T2=(LINE_FLUX(1)-LINE_FLUX(2))/LOG(R(1)/R(2))
	DO I=1,ND-2
	  T1=T2
	  T2=(LINE_FLUX(I)-LINE_FLUX(I+2))/LOG(R(I)/R(I+2))
	  EW=EW+LOG(R(I)/R(I+1))*( LINE_FLUX(I)+LINE_FLUX(I+1)
	1      +LOG(R(I)/R(I+1))*(T2-T1)/6.0_LDP )
	END DO
	T1=T2
	T2=(LINE_FLUX(ND-1)-LINE_FLUX(ND))/LOG(R(ND-1)/R(ND))
	EW=EW+LOG(R(ND-1)/R(ND))*( LINE_FLUX(ND-1)+LINE_FLUX(ND)
	1            +LOG(R(ND-1)/R(ND))*(T2-T1)/6.0_LDP )
	EW=EW*0.5_LDP
!
	EMISS_EW=0.0_LDP
	T2=(EMISS_LINE_FLUX(1)-EMISS_LINE_FLUX(2))/LOG(R(1)/R(2))
	DO I=1,ND-2
	  T1=T2
	  T2=(EMISS_LINE_FLUX(I)-EMISS_LINE_FLUX(I+2))/LOG(R(I)/R(I+2))
	  EMISS_EW=EMISS_EW+LOG(R(I)/R(I+1))*( EMISS_LINE_FLUX(I)+EMISS_LINE_FLUX(I+1)
	1      +LOG(R(I)/R(I+1))*(T2-T1)/6.0_LDP )
	END DO
	T1=T2
	T2=(EMISS_LINE_FLUX(ND-1)-EMISS_LINE_FLUX(ND))/LOG(R(ND-1)/R(ND))
	EMISS_EW=EMISS_EW+LOG(R(ND-1)/R(ND))*( EMISS_LINE_FLUX(ND-1)+EMISS_LINE_FLUX(ND)
	1            +LOG(R(ND-1)/R(ND))*(T2-T1)/6.0_LDP )
	EMISS_EW=EMISS_EW*0.5_LDP
!
	T1=16.0E+20_LDP*(3.14159_LDP)**2/RLUM/3.845E+33_LDP
	DO I=1,ND
	  FORCE_MULT(I)=FORCE_MULT(I)+GLINE(I)*T1*R(I)*R(I)*CHIL(I)/ESEC(I)
	END DO
!
	RETURN
	END
