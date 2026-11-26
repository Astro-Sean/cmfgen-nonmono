C This program will calculate an absorption profile by computing an optical
C depth at line center t0, with viogt profile.  It needs as input a
C micro-turbulent parameter b(km/sec), a central wavelength lam0(ang),
C an ossillator (sp) strength f(dimensionless), and a line width (rather a
C damping parameter) gam(sec^-1).  The temperature is also input but has little
C effect.  The voigt profile in IDL has a problem with
C negative independent variables so we do a onesided calculation and assume
C symmetry.
C Started on 1-9-92 by SRM
C
C This program is a modification of absprof.pro.  This will calculate a
C set of optical depths as a function of wavelength for the Lyman series of
C hydrogen.  When the optical depths get too large I truncate them at a
C particular maximum finite value so that when I exponentiate and
C convolve with the point spread function of the telescope I will (hopefully)
C get a well behaved transmission function which can be divided out of the
C EZ CMa spectrum to revel the ``true'' spectrum.  This modification started
C on 4-30-92 by SRM.
C
	SUBROUTINE H2ABS(WAVE,FLUX,NLAM,V_TURB,LOG_NTOT,T_IN_K)
	USE SET_KIND_MODULE
	IMPLICIT NONE
C
C Altered 16-Apr-2008: PAUSE used if cannot find file with line list.
C Altered 21-Dec-2001: Error in VTURB (**2 insetad of *2) fixed
C Altered 21-Apr-2000: H2_L_J made INTEGER.
C
	INTEGER NLAM
	REAL(KIND=LDP) WAVE(NLAM)
	REAL(KIND=LDP) FLUX(NLAM)
	REAL(KIND=LDP) V_TURB		!km/s
	REAL(KIND=LDP) LOG_NTOT
	REAL(KIND=LDP) T_IN_K
C
C Local variables
C
	INTEGER NH2
	PARAMETER (NH2=420)
	REAL(KIND=LDP) H2_LAM(NH2)
	REAL(KIND=LDP) H2_FREQ(NH2)
	REAL(KIND=LDP) H2_OSC(NH2)
	REAL(KIND=LDP) H2_GAM(NH2)
	REAL(KIND=LDP) H2_G(NH2)
	REAL(KIND=LDP) NU_DOP(NH2)
C
	INTEGER H2_L_J(NH2)
C
	REAL(KIND=LDP) G_J(0:7)
	REAL(KIND=LDP) N_J(0:7)
C
	REAL(KIND=LDP) CHIL(NH2)
C
	REAL(KIND=LDP) C_KMS,PI
	REAL(KIND=LDP) OPLIN
	REAL(KIND=LDP) TAU
	REAL(KIND=LDP) NTOT
	REAL(KIND=LDP) FREQ
	REAL(KIND=LDP) PHI
	REAL(KIND=LDP) a
	REAL(KIND=LDP) v
	REAL(KIND=LDP) T1
	REAL(KIND=LDP) HC_ON_K
	REAL(KIND=LDP) B
	REAL(KIND=LDP) POP_SUM
C
	INTEGER I,J,IOS
C
C Functions
C
	REAL(KIND=LDP) SPEED_OF_LIGHT
	REAL(KIND=LDP) FUN_PI
	REAL(KIND=LDP) VOIGT
	EXTERNAL SPEED_OF_LIGHT,FUN_PI,VOIGT
!
	CHARACTER(LEN=1) TMP_STR
C
C First read in the atomic data, kindly provided by Chuck Bowers (CWB).
C
C This is for the first 49 lines of the Lyman series.
C
C 5th column: Wavelength (Angstroms)
C 6th column: Oscilator strength
C 7th column: G (ground state)
C 8th column: Damping parameter
C
100	CONTINUE
	OPEN(UNIT=10,FILE='H2_IS_LINE_LIST',ACTION='READ',STATUS='OLD',IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(6,*)'H2_IS_LINE_LIST file not found'
	    WRITE(6,*)'Use astxt to assign file, then hit any character and return/enter'
	    READ(6,'(A)')TMP_STR
	    GOTO 100
	  END IF
	  DO I=1,NH2
	    READ(10,*)T1,T1,T1,H2_L_J(I),H2_LAM(I),H2_OSC(I),H2_G(I),H2_GAM(I)
	  END DO
	CLOSE(UNIT=10)
C
	C_KMS=1.0E-05_LDP*SPEED_OF_LIGHT()
	PI=FUN_PI()
C
        B=59.313_LDP        !molec. cnst. B [/cm]
        HC_ON_K=1.493_LDP        !hc/k [cm K]
        g_J(0)=1	!degeneracy (J" and I) of each level
        g_J(1)=9
        g_J(2)=5
        g_J(3)=21
        g_J(4)=9
        g_J(5)=33
        g_J(6)=13
        g_J(7)=45
	POP_SUM=0._LDP
        DO I=0,7
          N_J(I)=G_J(I)*EXP(-B*(I*(I+1))*HC_ON_K/T_IN_K)
	  POP_SUM=POP_SUM+N_J(I)
        END DO
C
C Normalize to the total population.
C
	NTOT=10.0E00_LDP**LOG_NTOT
	DO I=0,7
	  N_J(I)=NTOT*N_J(I)/POP_SUM
	END DO
C
	OPLIN=2.6540081E-02_LDP		!pi*e*e/m/c
	DO I=1,NH2
	  H2_FREQ(I)=0.01_LDP*C_KMS/H2_LAM(I)   		!10^15 Hz
	  NU_DOP(I)=12.85_LDP*H2_FREQ(I)*
	1             SQRT( (T_IN_K/1.0E+04_LDP)+ (V_TURB/12.85_LDP)**2 )/C_KMS
	  CHIL(I)=1.0E-15_LDP*OPLIN*H2_OSC(I)*N_J(H2_L_J(I))/SQRT(PI)/NU_DOP(I)
	END DO
C
	DO J=1,NLAM
	  TAU=0.0_LDP
	  FREQ=0.01_LDP*C_KMS/WAVE(J)
	  IF(WAVE(J) .GT. H2_LAM(NH2) .AND. WAVE(J) .LT. 1300)THEN
	    DO I=1,NH2
	      a=1.0E-15_LDP*H2_GAM(I)/4/PI/NU_DOP(I)
	      v=(FREQ-H2_FREQ(I))/NU_DOP(I)
	      PHI=VOIGT(a,v)
	      TAU=TAU+CHIL(I)*PHI
	    END DO
	    FLUX(J)=FLUX(J)*EXP(-TAU)
	  ELSE IF(WAVE(J) .LT. H2_LAM(NH2))THEN
	    FLUX(J)=0.0_LDP
	  END IF
	END DO
C
	RETURN
	END
