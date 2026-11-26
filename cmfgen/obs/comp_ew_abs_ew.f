!
! Routine to compute the the EW of a line using formal solution of the
! transfer equation using the differencing scheme of Mihalas, Kunasz, and Hummer 1975.
!
! It also computes the FLUX_EW -- this is used to check whether a line
! has an ifluence on the spectrum. The EW is not, buy itself, always
! useful since a P CYgni profile can have an EW of 0.
!
! NB - ETA when passed to this routine SHOULD NOT include the electron scattering emissivity.
!    - HQW contains the angular quadrature weights for the H computation at the MID points of the radial mesh.
!    - JQW contains the angular quadrature weights for the J computation on the radial mesh.
!    - The intrinsic line profile is assumed to be depth independent.
!
! JCONT should contain an estimat of J at line center. It will be
! updated by the code.
!
! The THK_LINE, THK_CONT options may not work.
!
	SUBROUTINE COMP_EW_ABS_EW(ETA,CHI,ESEC,CHIL,ETAL,V,SIGMA,R,P,
	1                  JCONT,EW,FLUX_EW,EW_CUT,CONT_INT,JQW,HQW,
	1                  PF,PROF,LFQW,WERFC,FL,TRANS_NAME,
	1                  DIF,DBB,IC,AMASS,
	1                  THK_LINE,THK_CONT,NLF,NC,NP,ND,METHOD)
	USE SET_KIND_MODULE
!
	IMPLICIT NONE
!
! Altered 15-Nov-2021: Fixed call to WRITE_VEC to fix vector/scaler issue.
! Altered 12-Sep-2019: Check whether the total opacity < 0. If it is we adjust
!                         TCHI to make it positive. As the routine integrates
!                         over the source function, it does not hadle lasing.
! Altered  6-Sep-2019: Now return FLUX_EW which if Int (ABS[F-F_c])/F_c dnu.
!                        Emission xomponent to EW= (FLUX_EW+EW)/2
!                        Absorption component to EW-(FLUX_EW-EW)/2
!
! Code assumes siherent scattering in the CMF. Not necessarily a good assumption
! when Vinf << 500 km [Vth(es)].
!
	INTEGER NLF,NC,NP,ND
	REAL(KIND=LDP) ETA(ND),CHI(ND),ESEC(ND),CHIL(ND),ETAL(ND)
	REAL(KIND=LDP) V(ND),SIGMA(ND),R(ND),P(NP)
	REAL(KIND=LDP) JCONT(ND)
	REAL(KIND=LDP) JQW(ND,NP)
	REAL(KIND=LDP) HQW(ND,NP)
	REAL(KIND=LDP) PF(NLF),PROF(NLF),LFQW(NLF),WERFC(NLF)
	REAL(KIND=LDP) DBB,IC,AMASS,FL
	REAL(KIND=LDP) EW		!Returned: in Angstroms
	REAL(KIND=LDP) ABS_EW   	!Returned: in Anstroms
	REAL(KIND=LDP) FLUX_EW		!Anstroms
	REAL(KIND=LDP) EW_CUT		!EW only roughly computed (no ES iteration).
	REAL(KIND=LDP) CONT_INT		!Returned: in Jy.
!
	LOGICAL DIF,THK_CONT,THK_LINE
	CHARACTER*(*) METHOD
	CHARACTER*(*) TRANS_NAME
!
	REAL(KIND=LDP) TA(ND),TB(ND),TC(ND),AV(ND),CV(ND),DTAU(ND),Z(ND)
	REAL(KIND=LDP) TCHI(ND),XM(ND),SOURCE(ND),U(ND),VB(ND),VC(ND)
	REAL(KIND=LDP) GB(ND),H(ND),GAM(ND),GAMH(ND),Q(ND),QH(ND)
	REAL(KIND=LDP) HCONT(ND)
	REAL(KIND=LDP) CVCONT(ND)
	REAL(KIND=LDP) dCHIdR(ND)
	REAL(KIND=LDP) AVM1(ND),CVM1(ND)
	REAL(KIND=LDP) WM(ND,ND),FB(ND,ND)
!
	REAL(KIND=LDP) OLD_JCONT(ND)
	REAL(KIND=LDP) JNU(ND,NLF)
	REAL(KIND=LDP) HNU(NLF)
	REAL(KIND=LDP) OLD_JNU(ND,NLF)
	REAL(KIND=LDP) FLUX_PROF(NLF)
	REAL(KIND=LDP) OLD_ABS_EW
	INTEGER NLF_PROF
!
! Local variables.
!
	INTEGER I,LS,ML,NI,NIEXT,IT
	REAL(KIND=LDP) OLDCHI,T1
	REAL(KIND=LDP) DBC,TOR,IBOUND,WERF_EXP
	REAL(KIND=LDP), PARAMETER :: EW_ACC=0.005_LDP            !i.e., 0.5%
	LOGICAL, SAVE :: FIRST=.TRUE.
	LOGICAL EQUAL
	EXTERNAL EQUAL
!
	FIRST=.FALSE.
	IF(TRANS_NAME .EQ. 'Ca2(3p6_4p_2Po[3/2]-3p6_3d_2De[3/2])')THEN
	  CALL WRITE_VEC(ESEC,ND,'ESEC',116)
	  CALL WRITE_VEC(ETA,ND,'ETA',116)
	  CALL WRITE_VEC(CHI,ND,'CHI',116)
	  CALL WRITE_VEC(ETAL,ND,'ETAL',116)
	  CALL WRITE_VEC(CHIL,ND,'CHIL',116)
	  CALL WRITE_VEC(PF,NLF,'PF',116)
	  CALL WRITE_VEC(PROF,NLF,'PROF',116)
	  CALL WRITE_VEC(LFQW,NLF,'LFQW',116)
	  TA(1)=DBB; I=1; CALL WRITE_VEC(TA,I,'DBB',116)
	  TA(1)=FL;  I=1; CALL WRITE_VEC(TA,I,'FL',116)
	  FIRST=.FALSE.
	  CALL WRITE_VEC(JCONT,ND,'JCONT',116)
	END IF
!
!	SOURCE=ETA/CHI
!	AV=ESEC/CHI	
!       CALL NEWJSOLD(TA,TB,TC,XM,WM,FB,OLD_JCONT,DTAU,R,Z,P,
!	1               SOURCE,AV,CHI,AVM1,JQW,
!	1               THK_CONT,DIF,DBB,IC,NC,ND,NP,METHOD)
!	WRITE(118,'(A2,2X,2ES14.4,ES18.8)')'JN',JCONT(1),OLD_JCONT(1),FL
!
!
! Get consistent estimate of the continuum J.
!
	IF(METHOD .NE. 'ZERO')THEN
	  CALL DERIVCHI(dCHIdR,CHI,R,ND,METHOD)
	END IF
!
! Enter loop to perform integration along each ray.
!
	DO IT=1,20
!	  WRITE(118,'(A2,3X,2ES14.4,ES18.8)')'J',JCONT(1),JCONT(ND),FL
	  OLD_JCONT=JCONT
	  JCONT=0.0_LDP
	  HCONT=0.0_LDP
	  DO I=1,ND
	    SOURCE(I)=(ETA(I)+ESEC(I)*OLD_JCONT(I))/CHI(I)
	  END DO
!	  WRITE(118,'(A2,2X,2ES14.4,ES18.8)')'S',SOURCE(1),SOURCE(ND),FL
!
	  DO LS=1,MIN(NP,ND+NC-2)
	    NI=ND-(LS-NC-1)
	    IF(LS .LE. NC)NI=ND
	    CALL ZALONGP(R,Z,P(LS),NI)
!
! Determine boundary condition for continuum intensity.
!
	    IF(THK_CONT)THEN
	      IF(P(LS) .GT. 0.0_LDP)THEN
	        TOR=CHI(1)*R(1)*R(1)*(1.570796_LDP-ACOS(P(LS)/R(1)))/P(LS)
	      ELSE
	        TOR=CHI(1)*R(1)
	      END IF
	      IBOUND=ETA(1)*(1.0_LDP-EXP(-TOR))/CHI(1)
	    ELSE
	      TOR=0.0_LDP
	      IBOUND=0.0_LDP
	    END IF
!
	    IF(DIF .AND. LS .LE. NC)THEN
	        DBC=DBB*SQRT(R(ND)*R(ND)-P(LS)*P(LS))/R(ND)/CHI(ND)
	    END IF
	    IF(METHOD .EQ. 'ZERO')THEN
	      CALL TAU(DTAU,CHI,Z,NI)
	    ELSE
	      CALL NORDTAU(DTAU,CHI,Z,R,dCHIdR,NI)
	    END IF
	    CALL TCOMPD(TA,TB,TC,DTAU,DIF,LS,NC,ND,NI)
	    CALL XVECD(DTAU,SOURCE,XM,DIF,DBC,IC,LS,NC,ND,NI)
	    XM(1)=-IBOUND
!
	    CALL THOMAS(TA,TB,TC,XM,NI,1)
	    DO I=1,NI-1
	      CV(I)=(XM(I+1)-XM(I))/DTAU(I)
	    END DO
!
	    DO I=1,NI
	       JCONT(I)=JCONT(I)+JQW(I,LS)*XM(I)
	    END DO
	    HCONT(1)=HCONT(1)+CV(1)*HQW(1,LS)
	    CVCONT(1)=CV(1)
	  END DO  		!LS
!
	  T1=0.0_LDP
	  DO I=1,ND
	    T1=MAX(T1, ABS(1.0_LDP-JCONT(I)/OLD_JCONT(I)))
	  END DO
	  IF(T1 .LT. 1.0E-06_LDP)EXIT
!	  WRITE(118,*)ETA(1),JCONT(1),HCONT(1)
	END DO			!IT
!
! 
! *************
!
! Enter loop to perform integration along each ray.
!
	DO ML=1,NLF
	  JNU(:,ML)=JCONT(:)
	END DO
!
	DO ML=2,NLF
	  NLF_PROF=ML
	  IF(PROF(ML-1) .NE. 0.0_LDP .AND. PROF(ML) .EQ. 0.0_LDP)EXIT
	END DO
!
	OLD_ABS_EW=0.0_LDP
	DO IT=1,10
	  OLD_JNU=JNU
	  JNU=0.0_LDP
	  HNU=0.0_LDP
	  EW=0.0_LDP; ABS_EW=0.0_LDP; FLUX_PROF=0.0_LDP

	  DO LS=1,MIN(NP,ND+NC-2)
	    NI=ND-(LS-NC-1)
	    IF(LS .LE. NC)NI=ND
!
! NIEXT is used to compute one extra value of TCHI so that a more accurate
! dCHIdR can be computed.
!
	    NIEXT=NI+1
	    IF(NI .EQ. ND)NIEXT=NI
!
! Zero AV, CV, dUd... and dVd... vectors.
!
	    AV(1:NI)=0.0_LDP
	    CV(1:NI)=0.0_LDP
!
	    CALL ZALONGP(R,Z,P(LS),NI)
	    CALL GAMMA(GAM,GAMH,SIGMA,Z,R,V,ND,NI)
	    GAM=GAM*FL; GAMH=GAMH*FL
!
! Determine boundary condition for continuum intensity.
!
	    IF(THK_CONT)THEN
	      IF(P(LS) .GT. 0.0_LDP)THEN
	        TOR=CHI(1)*R(1)*R(1)*(1.570796_LDP-ACOS(P(LS)/R(1)))/P(LS)
	      ELSE
	        TOR=CHI(1)*R(1)
	      END IF
	      IBOUND=ETA(1)*(1.0_LDP-EXP(-TOR))/CHI(1)
	    ELSE
	      TOR=0.0_LDP
	      IBOUND=0.0_LDP
	    END IF
!
!*********
! 
!
! Perform integration for each frequency in turn.
! This section loops over frequencies covering the intrisic line profile.
!
	    OLDCHI=CHI(NI)
	    DO ML=1,NLF_PROF
	      T1=PROF(ML)
	      DO I=1,NIEXT
	        TCHI(I)=CHI(I)+CHIL(I)*T1
	        IF(TCHI(I) .LT. 1.0E-10_LDP*CHI(I))TCHI(I)=0.001_LDP*CHI(I)
	        SOURCE(I)=(ETA(I)+ETAL(I)*T1+ESEC(I)*OLD_JNU(I,ML))/TCHI(I)
	      END DO
	      CALL QKIM(Q,QH,GAM,GAMH,TCHI,PF,ML,NI,NLF)
	      IF(DIF .AND. LS .LE. NC)THEN
	        DBC=DBB*SQRT(R(ND)*R(ND)-P(LS)*P(LS))/R(ND)/TCHI(ND)
	1            *(1.0_LDP+Q(NI)*(1.0_LDP-TCHI(NI)/OLDCHI))
	      END IF
	      IF(METHOD .EQ. 'ZERO')THEN
	        CALL TAU(DTAU,TCHI,Z,NI)
	      ELSE
	        CALL DERIVCHI(dCHIdR,TCHI,R,NIEXT,METHOD)
	        CALL NORDTAU(DTAU,TCHI,Z,R,dCHIdR,NI)
	      END IF
	      CALL TUVGHD(TA,TB,TC,U,VB,VC,GB,H,Q,QH,DTAU,DIF,LS,NC,NI)
	      CALL XVECD(DTAU,SOURCE,XM,DIF,DBC,IC,LS,NC,ND,NI)
!
!	I(incident)=ETAL(1)/CHIL(1)*(1.0D0-WERF_EXP)+IBOUND*WERF_EXP
!
	      IF(THK_LINE)THEN
	        WERF_EXP=EXP(1.0E-15_LDP*CHIL(1)/FL/GAM(1)*WERFC(ML))
	        XM(1)=(CHI(1)/TCHI(1)*WERF_EXP-1.0_LDP)*ETAL(1)/CHIL(1)
	1              -IBOUND*CHI(1)/TCHI(1)*WERF_EXP
	      ELSE
	        XM(1)=-IBOUND
	      END IF
!
! Update AV matrix.
!
	      AV(1)=XM(1)+U(1)*AV(1)
	      DO I=2,NI-1
	        AV(I)=XM(I)+U(I)*AV(I)-VB(I)*CV(I-1)-VC(I)*CV(I)
	      END DO
	      AV(NI)=XM(NI)+U(NI)*AV(NI)-VB(NI)*CV(NI-1)
!
! Solve for the radiation field along ray for this frequency.
!
	      CALL THOMAS(TA,TB,TC,AV,NI,1)
!
! Update C matrix.
!
	      DO I=1,NI-1
	        CV(I)=GB(I)*(AV(I)-AV(I+1))+H(I)*CV(I)
	      END DO
!
	      DO I=1,NI
	        JNU(I,ML)=JNU(I,ML)+AV(I)*JQW(I,LS)
	      END DO
	      HNU(ML)=HNU(ML)+CV(1)*HQW(1,LS)
	      IF(ML .EQ. 1)CVCONT(1)=CV(1)
	      EW=EW+(CV(1)-CVCONT(1))*HQW(1,LS)*LFQW(ML)
	      ABS_EW=ABS_EW+ABS(CV(1)-CVCONT(1))*HQW(1,LS)*LFQW(ML)
	      FLUX_PROF(ML)=FLUX_PROF(ML)+(CV(1)-CVCONT(1))*HQW(1,LS)
!
	      DO I=1,NI
	        AVM1(I)=AV(I)
	        CVM1(I)=CV(I)
	      END DO
!
	      OLDCHI=TCHI(NI)
	    END DO  		!End DO ML
! 
!
! Perform integration for each frequency in turn.
!
!
	    IF(METHOD .EQ. 'ZERO')THEN
	      CALL TAU(DTAU,CHI,Z,NI)
	    ELSE
	      CALL DERIVCHI(dCHIdR,CHI,R,NIEXT,METHOD)
	      CALL NORDTAU(DTAU,CHI,Z,R,dCHIdR,NI)
	    END IF
	    ML=NLF_PROF+1
	    CALL QKIM(Q,QH,GAM,GAMH,CHI,PF,ML,NI,NLF)
	    CALL TUVGHD(TA,TB,TC,U,VB,VC,GB,H,Q,QH,DTAU,DIF,LS,NC,NI)
	    SOURCE=1.0_LDP  !Dummy values
	    CALL THOMAS(TA,TB,TC,SOURCE,NI,1)
!
	    OLDCHI=CHI(NI)
	    IF(DIF .AND. LS .LE. NC)THEN
	      DBC=DBB*SQRT(R(ND)*R(ND)-P(LS)*P(LS))/R(ND)/CHI(ND)
	    END IF
!
! We now intgerate over the rest of the lines profile.
!
	    DO ML=NLF_PROF+1,NLF
	      DO I=1,NIEXT
	        SOURCE(I)=(ETA(I)+ESEC(I)*OLD_JNU(I,ML))/CHI(I)
	      END DO
	      CALL XVECD(DTAU,SOURCE,XM,DIF,DBC,IC,LS,NC,ND,NI)
!
!	I(incident)=ETAL(1)/CHIL(1)*(1.0D0-WERF_EXP)+IBOUND*WERF_EXP
!
	      IF(THK_LINE)THEN
	        WERF_EXP=EXP(1.0E-15_LDP*CHIL(1)/FL/GAM(1)*WERFC(ML))
	        XM(1)=(CHI(1)/TCHI(1)*WERF_EXP-1.0_LDP)*ETAL(1)/CHIL(1)
	1              -IBOUND*CHI(1)/TCHI(1)*WERF_EXP
	      ELSE
	        XM(1)=-IBOUND
	      END IF
!
! Update AV matrix.
!
	      AV(1)=XM(1)+U(1)*AV(1)
	      DO I=2,NI-1
	        AV(I)=XM(I)+U(I)*AV(I)-VB(I)*CV(I-1)-VC(I)*CV(I)
	      END DO
	      AV(NI)=XM(NI)+U(NI)*AV(NI)-VB(NI)*CV(NI-1)
!
! Solve for the radiation field along ray for this frequency.
!
	      CALL SIMPTH(TA,TB,TC,AV,NI,1)
!
! Update C matrix.
!
	      DO I=1,NI-1
	        CV(I)=GB(I)*(AV(I)-AV(I+1))+H(I)*CV(I)
	      END DO
!
	      DO I=1,NI
	        JNU(I,ML)=JNU(I,ML)+AV(I)*JQW(I,LS)
	      END DO
	      HNU(ML)=HNU(ML)+CV(1)*HQW(1,LS)
!
	      DO I=1,NI
	        AVM1(I)=AV(I)
	        CVM1(I)=CV(I)
	      END DO
!
	      EW=EW+(CV(1)-CVCONT(1))*HQW(1,LS)*LFQW(ML)
	      ABS_EW=ABS_EW+ABS(CV(1)-CVCONT(1))*HQW(1,LS)*LFQW(ML)
	      FLUX_PROF(ML)=FLUX_PROF(ML)+(CV(1)-CVCONT(1))*HQW(1,LS)
	      OLDCHI=CHI(NI)
	    END DO  		!End DO ML
	  END DO		!End DO LS
!
! Evaluate the line EW. The units are Angstroms. Also evaluate
! the continuum intensity ( Jys/kpc/kpc ). Note that H is
! defined midway between R(1) and R(2).
!
          EW=2.99792458E-12_LDP*EW/HCONT(1)/FL/FL
          ABS_EW=2.99792458E-12_LDP*ABS_EW/HCONT(1)/FL/FL
          CONT_INT=13.19868_LDP*HCONT(1)*( (R(1)+R(2))**2 )/4.0_LDP
	  IF(ABS_EW .LT. EW_CUT)EXIT
	  IF( EQUAL(ABS_EW,OLD_ABS_EW,EW_ACC) )EXIT
!
	END DO			!IT
!
	FLUX_EW=0.0_LDP
	DO ML=1,NLF
	  FLUX_EW=FLUX_EW+ABS(FLUX_PROF(ML))*LFQW(ML)
	END DO
        FLUX_EW=2.99792458E-12_LDP*FLUX_EW/HCONT(1)/FL/FL
!	
        WRITE(119,'(4ES14.6,F14.4,T75,A)')CONT_INT,EW,FLUX_EW,ABS_EW,2.99792458D+03/FL,TRIM(TRANS_NAME)
        FLUSH(UNIT=119)
	IF( ABS(2.99792458E+03_LDP/FL-7774.0827_LDP) .LT. 0.0001_LDP)THEN
	  DO ML=1,NLF
	    WRITE(143,'(ES14.4,F14.4,2ES14.4)')PF(ML),2.99792458D+03/PF(ML),13.19868D0*HNU(ML)*( (R(1)+R(2))**2 )/4.0D0,PROF(ML)
	  END DO
	  FLUSH(UNIT=143)
	END IF
!
	RETURN
	END
