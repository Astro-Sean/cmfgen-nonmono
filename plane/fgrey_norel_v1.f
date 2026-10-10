!
! This routine solves for the mean intensity as a function of depth using the
! Feautrier Technique for a spherical gray atmosphere. A diffusion approximation
! is used for the lower boundary condition. Radiative equilibrium is assumed.
!
! The time dependent radiative transfer equation is solved including all terms
! to first order in v/c, and assumes a Hubble flow. Thus sigma=dlnv/dlnr-1 =0.
! A Langrangian description is used to handle the time dependence. The frequency
! integrated moments J and H (and the associated boundary conditions) at the previous
! time step need to be available, and are read in from the file JH_AT_PREV_TIME.
! If desired, the full D/Dt term can be neglected, in which case the code should give
! identical answers to JGREY_WITH_FVT.
!
! Partially based on the routine ../subs/jgrey_with_fvt.f
!
	SUBROUTINE FGREY_NOREL_V1(FEDD,H_ON_J,RJ,CHI,R,VEL,SIGMA,
	1              P,JQW,HQW,KQW,LUMINOSITY,IC,METHOD,
	1              H_OUTBC,H_INBC,DIFF_APPROX,ND,NC,NP)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Created 22-July-2006
!
	INTEGER NC
	INTEGER ND
	INTEGER NP
!
	REAL(KIND=LDP) FEDD(ND)
	REAL(KIND=LDP) H_ON_J(ND)
!
	REAL(KIND=LDP) R(ND)			!Radius grid (in units of 10^10 cm)
	REAL(KIND=LDP) RJ(ND)			!Mean opacity
	REAL(KIND=LDP) CHI(ND)			!Opacity
	REAL(KIND=LDP) VEL(ND)			!Velocity (in km/s)
	REAL(KIND=LDP) SIGMA(ND)		!dlnV/dlnr-1
!
	REAL(KIND=LDP) P(NP)			!Impact parameters
	REAL(KIND=LDP) JQW(ND,NP)		!Quadrature weight for J (on grid)
	REAL(KIND=LDP) KQW(ND,NP) 		!Quadrature weight for K (on grid)
	REAL(KIND=LDP) HQW(ND,NP)		!Quadrature weight for H (at midpoints)
!
	REAL(KIND=LDP) LUMINOSITY               !Luminosity at inner boundary in Lsun.
	REAL(KIND=LDP) IC                       !Not used
	LOGICAL DIFF_APPROX		!Use a diffusion approximation (as opposed to a Schuster core)
	CHARACTER(LEN=6) METHOD
!
	REAL(KIND=LDP) H_OUTBC		!Eddington factor for H at outer boundary.
	REAL(KIND=LDP) H_INBC
!
! Local vectors & arrays
!
	REAL(KIND=LDP) Z(ND)
	REAL(KIND=LDP) TA(ND)
	REAL(KIND=LDP) TB(ND)
	REAL(KIND=LDP) TC(ND)
	REAL(KIND=LDP) XM(ND)
	REAL(KIND=LDP) CHI_MOD(ND)
	REAL(KIND=LDP) dCHIdR(ND)
	REAL(KIND=LDP) dCHI_MODdR(ND)
	REAL(KIND=LDP) Q(ND)
	REAL(KIND=LDP) F(ND)
	REAL(KIND=LDP) DTAU(ND)
	REAL(KIND=LDP) BETA(ND)
	REAL(KIND=LDP) VU(ND)
	REAL(KIND=LDP) CV(ND)
!
! FS indicates the following quantities (J, H, & K) have been computed using the
! formal solution.
!
	REAL(KIND=LDP) FS_J(ND)
	REAL(KIND=LDP) FS_H(ND)
	REAL(KIND=LDP) FS_K(ND)
!
	REAL(KIND=LDP) IBOUND
	REAL(KIND=LDP) DBB
	REAL(KIND=LDP) DBC
	REAL(KIND=LDP) C_KMS
!
	REAL(KIND=LDP) PI
	REAL(KIND=LDP) T1,T2,T3
	REAL(KIND=LDP) E1,E2,E3
	REAL(KIND=LDP) SPEED_OF_LIGHT
	EXTERNAL SPEED_OF_LIGHT
!
	INTEGER I
	INTEGER NI
	INTEGER LS
	INTEGER, PARAMETER :: IONE=1
	INTEGER LUER,ERROR_LU,LU
	EXTERNAL ERROR_LU
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
!
! Initialize values.
!
	PI=ACOS(-1.0_LDP)
	LUER=ERROR_LU()
	C_KMS=1.0E-05_LDP*SPEED_OF_LIGHT()
	DO I=1,ND
	  BETA(I)=VEL(I)/C_KMS
	END DO
	IBOUND=0.0_LDP
!
	FS_J(:)=0.0_LDP
	FS_H(:)=0.0_LDP
	FS_K(:)=0.0_LDP
	H_INBC=0.0_LDP
	H_OUTBC=0.0_LDP
!
! DBB =3L/16(piR)**2 and is used for the lower boundary diffusion approximation.
!
	T1=3.826E+13_LDP*LUMINOSITY/16.0_LDP/PI/PI
	DBB=3.0_LDP*T1/R(ND)/R(ND)
!
	CALL DERIVCHI(dCHIdr,CHI,R,ND,METHOD)
!
! Enter loop for each impact parameter P
!
	DO 2000 LS=1,NP
	  NI=ND-(LS-NC-1)
	  IF(LS .LE. NC)THEN
	    NI=ND
	    DBC=DBB*SQRT(R(ND)*R(ND)-P(LS)*P(LS))/R(ND)
	  END IF
!
! Compute Z for this impact parameter
!
	  IF(NI .GT. 1)THEN
	    DO I=1,NI
	      Z(I)=SQRT( (R(I)-P(LS))*(R(I)+P(LS)) )
	      T1=Z(I)/R(I)
	      TA(I)=BETA(I)/R(I)
	      TA(I)=4.0_LDP*TA(I)*( (1.0_LDP-T1)**2 + T1*T1*TA(I) )
	      CHI_MOD(I)=CHI(I)  + TA(I)
	    END DO
	  END IF
!
	  IF(NI .GT. 2)THEN
	    CALL DERIVCHI(TB,TA,R,NI,METHOD)
	    DO I=1,NI
	      dCHI_MODdR(I)=dCHIdR(I)+TB(I)
	    END DO
	  END IF
!
! Compute Z for this impact parameter
!
	  IF(NI .GT. 2)THEN
	    DO I=1,NI-1
	      DTAU(I)=0.5_LDP*(Z(I)-Z(I+1))*(CHI_MOD(I)+CHI_MOD(I+1)+(Z(I)-Z(I+1))
	1       *(dCHI_MODdR(I+1)*Z(I+1)/R(I+1)-dCHI_MODdR(I)*Z(I)/R(I))/6.0_LDP)
	    END DO
	  END IF
!
! Compute TA (tridiagonal matrix) and XM vector. Since it is a gray
! atmosphere, the source function is simply RJ*CHI(I)/CHI_MOD(I).
!
	  IF(NI .GT. 2)THEN
!
	    DO I=1,NI-1
	      VU(I)=1.0_LDP/DTAU(I)
	    END DO
!
	    XM(1)=-IBOUND
	    TA(1)=0.0_LDP
	    TC(1)=1.0_LDP/DTAU(1)
	    TB(1)=-1.0_LDP-TC(1)
	    DO I=2,NI-1
	      TA(I)=VU(I-1)
	      TC(I)=VU(I)
	      TB(I)=-0.5_LDP*(DTAU(I-1)+DTAU(I))-TA(I)-TC(I)
	      XM(I)=-RJ(I)*CHI(I)*(DTAU(I-1)+DTAU(I))*0.5_LDP/CHI_MOD(I)
	    END DO
!
	    IF(LS .LE. NC .AND. DIFF_APPROX)THEN
	      TB(NI)=VU(NI-1)
	      TA(NI)=-VU(NI-1)
	      XM(NI)=DBC
	    ELSE IF(LS .GT. NC)THEN
	      TA(NI)=-TC(NI-1)
	      TB(NI)=-TA(NI)+DTAU(NI-1)/2.0_LDP
	      XM(NI)=0.5_LDP*DTAU(NI-1)*RJ(NI)*CHI(NI)/CHI_MOD(NI)
	    ELSE
	      TA(NI)=-VU(NI-1)
	      TB(NI)=1.0_LDP+VU(NI-1)
	      XM(NI)=IC
	    END IF
	    TC(NI)=0.0_LDP
!
! Solve the tridiagonal system of equations.
!
	    CALL THOMAS(TA,TB,TC,XM,NI,IONE)
!
	  ELSE IF(NI .EQ. 1)THEN
	    XM(1)=0.0_LDP
	  ELSE IF(NI .EQ. 2)THEN
	    Z(1)=SQRT(R(1)*R(1)-P(LS)*P(LS))
	    DTAU(1)=0.5_LDP*Z(1)*(CHI_MOD(1)+CHI_MOD(2))		!Z(2)=0.0
	    E1=EXP(-DTAU(1))
	    E2=1.0_LDP-(1.0_LDP-E1)/DTAU(1)
	    E3=(1.0_LDP-E1)/DTAU(1)-E1
	    IF(DTAU(1) .LT. 1.0E-03_LDP)THEN
	      E2=DTAU(1)*0.5_LDP+DTAU(1)*DTAU(1)/6.0_LDP
	      E3=DTAU(1)*0.5_LDP-DTAU(1)*DTAU(1)/3.0_LDP
	    END IF
C
	    VU(1)=1.0_LDP/DTAU(1)
	    TA(2)=CHI(2)*RJ(2)/CHI_MOD(2)
	    TA(1)=CHI(1)*RJ(1)/CHI_MOD(1)
	    XM(2)=TA(2)*E2+TA(1)*E3
            XM(1)=0.5_LDP*(XM(2)*E1+TA(1)*E2+TA(2)*E3)
	  END IF
!
! Update J, H, K, & N for this angle.
!
	  DO I=1,NI
	    FS_J(I)=FS_J(I)+JQW(I,LS)*XM(I)
	    FS_K(I)=FS_K(I)+KQW(I,LS)*XM(I)
	  END DO
!
! CV is the flux (v in usual Mihalas notation).
!
	  DO I=1,NI-1
	    CV(I)=VU(I)*(XM(I+1)-XM(I))
	    FS_H(I)=FS_H(I)+HQW(I,LS)*CV(I)
	  END DO
!
	  H_OUTBC=H_OUTBC+JQW(1,LS)*(XM(1)-IBOUND)*Z(1)/R(1)
	  IF(.NOT. DIFF_APPROX)H_INBC=H_INBC+JQW(ND,LS)*(IC-XM(ND))*Z(ND)/R(ND)
!
2000	CONTINUE
!
! Compute the new Feautrier factors.
!
	FEDD=FS_K/FS_J
	DO I=1,ND-1
	  TA(I)=0.25_LDP*FS_H(I)*(R(I)+R(I+1))**2
	END DO
	CALL REGRID_H(H_ON_J,R,TA,H_OUTBC,DBB,ND,TB)
	H_ON_J=H_ON_J/FS_J/R/R
!
! Compute the factor for the outer boundary condition.
!
	H_OUTBC=H_OUTBC/FS_J(1)
	IF(.NOT. DIFF_APPROX)H_INBC=H_INBC/FS_J(ND)
!
	CALL GET_LU(LS,'MOM_JREL_GREY_V1')
        OPEN(UNIT=LS,FILE='GreyDiagnostics',STATUS='UNKNOWN',POSITION='APPEND')
	  WRITE(LS,'(A,13X,A,12X,A,10X,A,9(1X,A))')'    I','R','V','CHI','        V/cR',
	1          '          RJ','       J(FS)','    H(FS)',
	1          '           f','      H_ON_J'
	  DO I=1,ND
	    WRITE(LS,'(I5,ES14.5,11ES13.4)')I,R(I),VEL(I),CHI(I),VEL(I)/R(I)/C_KMS,
	1                                RJ(I),FS_J(I),FS_H(I),FEDD(I),H_ON_J(I)
	  END DO
	  T1=4.1274E-12_LDP*R(1)*R(1)*( FS_H(1)+BETA(1)*FS_J(1)*
	1         (1.0_LDP+FEDD(1))/(1.0_LDP-BETA(1)*BETA(1)) )
	  WRITE(LS,'(A)')' '
	  WRITE(LS,'(A,ES14.5)')' Observed flux computed by EDD routine is ',T1
	  T1=4.1274E-12_LDP*R(1)*R(1)*FS_H(1)
	  WRITE(LS,'(A,ES14.5)')' Comoving-frame flux computed by EDD routine is ',T1
	  WRITE(LS,'(A)')' '
	CLOSE(UNIT=LS)
!
	RETURN
	END
