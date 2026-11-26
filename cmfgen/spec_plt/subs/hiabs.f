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
	SUBROUTINE HIABS(WAVE,FLUX,NLAM,V_TURB,LOG_NTOT,T_IN_K)
	USE SET_KIND_MODULE
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
C
C Altered 16-Apr-2008: PAUSE used if cannot find file with line list.
C Altered 21-Dec-2001: Error in VTURB (**2 insetad of *2) fixed
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
	INTEGER NHYD
	PARAMETER (NHYD=49)
	REAL(KIND=LDP) HYD_LAM(NHYD)
	REAL(KIND=LDP) HYD_FREQ(NHYD)
	REAL(KIND=LDP) HYD_OSC(NHYD)
	REAL(KIND=LDP) HYD_GAM(NHYD)
	REAL(KIND=LDP) HYD_G(NHYD)
	REAL(KIND=LDP) NU_DOP(NHYD)
C
	REAL(KIND=LDP) CHIL(NHYD)
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
C
	INTEGER I,J,IOS
	LOGICAL DO_SYS
	CHARACTER(LEN=100) STRING
C
C Functions
C
	REAL(KIND=LDP) SPEED_OF_LIGHT
	REAL(KIND=LDP) FUN_PI
	REAL(KIND=LDP) VOIGT
	EXTERNAL SPEED_OF_LIGHT,FUN_PI,VOIGT
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
	OPEN(UNIT=10,FILE='HI_IS_LINE_LIST',ACTION='READ',STATUS='OLD',IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(6,*)'HI_IS_LINE_LIST file not found'
	    CALL GET_ENVIRONMENT_VARIABLE('CMFDIST',STRING)
            STRING=TRIM(STRING)//'/com/assign_txt_files_tcsh.sh'
	    WRITE(6,*)'I can assign the file with the following command'
	    WRITE(6,*)'    ',TRIM(STRING)
	    DO_SYS=.TRUE.
	    CALL GEN_IN(DO_SYS,'Do you want me to issue the CMFGEN assignment command')
	    IOS=0
	    IF(DO_SYS)CALL SYSTEM(STRING)
	    IF(.NOT. DO_SYS .OR. IOS .NE. 0)THEN
	      WRITE(6,*)'Use astxt to assign file, then hit any character and return/enter'
	      READ(5,'(A)')TMP_STR
	    END IF
	    GOTO 100
	  END IF
	  DO I=1,NHYD
	    READ(10,*)T1,T1,T1,T1,HYD_LAM(I),HYD_OSC(I),HYD_G(I),HYD_GAM(I)
	  END DO
	CLOSE(UNIT=10)
C
	C_KMS=1.0E-05_LDP*SPEED_OF_LIGHT()
	PI=FUN_PI()
C
	OPLIN=2.6540081E-02_LDP		!pi*e*e/m/c
	NTOT=10.0E00_LDP**LOG_NTOT
	DO I=1,NHYD
	  HYD_FREQ(I)=0.01_LDP*C_KMS/HYD_LAM(I)   		!10^15 Hz
	  NU_DOP(I)=12.85_LDP*HYD_FREQ(I)*
	1             SQRT( (T_IN_K/1.0E+04_LDP)+ (V_TURB/12.85_LDP)**2 )/C_KMS
	  CHIL(I)=1.0E-15_LDP*OPLIN*HYD_OSC(I)*NTOT/SQRT(PI)/NU_DOP(I)
	END DO
C
	DO J=1,NLAM
	  TAU=0.0_LDP
	  FREQ=0.01_LDP*C_KMS/WAVE(J)
	  IF(WAVE(J) .GT. HYD_LAM(NHYD) .AND. WAVE(J) .LT. 1300)THEN
	    DO I=1,NHYD
	      a=1.0E-15_LDP*HYD_GAM(I)/4/PI/NU_DOP(I)
	      v=(FREQ-HYD_FREQ(I))/NU_DOP(I)
	      PHI=VOIGT(a,v)
	      TAU=TAU+CHIL(I)*PHI
	    END DO
	    FLUX(J)=FLUX(J)*EXP(-TAU)
	  ELSE IF(WAVE(J) .LT. HYD_LAM(NHYD))THEN
	    FLUX(J)=0.0_LDP
	  END IF
	END DO
C
	RETURN
	END
