!
! Routine to perform Gauss smoothing. Utilizes routines
! originally written by Jim Herald.
!
! Routine can use either:
!        (1) A gaussian of fixed width dX
!        (2) A gaussian of fixed Resolution (i.e. fixed R).
!
! X and Y data should be in linear space.
!
	SUBROUTINE DO_GAUSS_SMOOTH(XMIN,XMAX,IP,OP,T_OUT)
	USE SET_KIND_MODULE
	USE NEW_GEN_IN_INTERFACE
	USE MOD_CURVE_DATA
	USE MOD_COLOR_PEN_DEF
	USE MOD_SMEAR_PG
	USE LINE_ID_MOD
	IMPLICIT NONE
!
! Finalized: 21-Mar-2023
!
	REAL(KIND=LDP) XMIN,XMAX
	INTEGER IP,OP
	INTEGER T_OUT
!
	REAL(KIND=LDP), ALLOCATABLE :: XV(:)
	REAL(KIND=LDP), ALLOCATABLE :: YV(:)
!
	REAL(KIND=LDP) INST_RES
	REAL(KIND=LDP) RESOLUTION
	REAL(KIND=LDP) WAVE_MIN
	REAL(KIND=LDP) WAVE_MAX
	REAL(KIND=LDP) MIN_RES_KMS
	REAL(KIND=LDP) NUM_RES
	REAL(KIND=LDP) VSINI
	REAL(KIND=LDP) EPSILON
	REAL(KIND=LDP) C_CMS
	REAL(KIND=LDP) T1,T2
	REAL(KIND=LDP) XOFFSET
!
	INTEGER I,J
	INTEGER NX
	LOGICAL FFT_CONVOLVE
!
	IOS=0
	C_CMS=2.99792458D+10
	FFT_CONVOLVE=.FALSE.
	NUM_RES=5.0
	VSINI=0.0D0            !For rotational broadening, so set to zero
        EPSILON=0.0D0
	MIN_RES_KMS=1.0D0
!
! If the RESOLUTION is non zero, it gets used instead of INST_RES.
! Resolution allows a spectrum to be smoothed to a constant velocity
! resolution at all wavelengths.
!
	INST_RES=0.0D0
	RESOLUTION=0.0D0
100     CONTINUE
	CALL NEW_GEN_IN(RESOLUTION,'Resolution [X/dX(FWHM)] (km/s if -ve)')
	IF(RESOLUTION .EQ. 0.0D0)THEN
	  CALL NEW_GEN_IN(INST_RES, 'Instrumental resolution (plot units')
	END IF
!
	IF(RESOLUTION .LT. 0.0D0)THEN
	  RESOLUTION=1.0D-05*C_CMS/ABS(RESOLUTION)
	ELSE IF(RESOLUTION .EQ. 0.0D0 .AND. INST_RES .EQ. 0.0D0)THEN
	  WRITE(T_OUT,*)'Only one INST_RES and RES can be zero'
	  GOTO 100
	ELSE IF(RESOLUTION .NE. 0.0D0 .AND. INST_RES .NE. 0.0D0)THEN
	  WRITE(T_OUT,*)'Only one INST_RES and RES can be non-zero'
	  GOTO 100
	END IF
!
! Use display range as defaults.
!
	WAVE_MIN=XMIN; WAVE_MAX=XMAX
	CALL NEW_GEN_IN(WAVE_MIN,'Minimum Wavelength')
	CALL NEW_GEN_IN(WAVE_MAX,'Minimum Wavelength')
!
	ALLOCATE(XV(NPTS(IP)),YV(NPTS(IP)))
!
! SMEAR assumes that XV is passed in freuqncy units, hence the
! conversion. In practice, any units can be passed, provided
! there are no zero values
!
! This sections add an offset so that division by zero does not occur.
! Only valid when the smoothing value, dX, is constant.
!
	T1=MINVAL(CD(IP)%XVEC)
	IF(T1 .LT. 1.0D-30)THEN
	  IF(RESOLUTION .NE. 0)THEN
	    WRITE(6,*)' Error -- this routine can only be used when X > 0 for all X'
	    WRITE(6,*)'    when RESOLUTON is set'
	    IOS=1
	    RETURN
	  ELSE
	    T1=T1*1.0001D0
	    XOFFSET=ABS(MIN(CD(IP)%XVEC(2),CD(IP)%XVEC(NPTS(IP)-1)))
            XOFFSET=MAX(ABS(T1),XOFFSET)
	    XV=1.0D-07*C_CMS/(CD(IP)%XVEC+XOFFSET)
	    WAVE_MIN=WAVE_MIN+XOFFSET
	    WAVE_MAX=WAVE_MAX+XOFFSET
	  END IF
	ELSE
	  XV=1.0D-07*C_CMS/CD(IP)%XVEC
	END IF
!
	YV=CD(IP)%DATA
	NX=NPTS(IP)
!
!        CALL SMEAR_V2(XV,YV,NX,
!	1       WAVE_MAX,WAVE_MIN,
!	1       INST_RES,RESOLUTION,VSINI,EPSILON,
!	1       MIN_RES_KMS,NUM_RES,FFT_CONVOLVE)
!
        CALL SMEAR_V3_PG(XV,YV,NX,
	1       WAVE_MAX,WAVE_MIN,
	1       INST_RES,RESOLUTION,VSINI,EPSILON,
	1       MIN_RES_KMS,NUM_RES)
!
	IF(OP .EQ. IP)THEN
	  CD(OP)%DATA=YV
	ELSE
	  IF(ALLOCATED(CD(OP)%DATA))DEALLOCATE(CD(OP)%DATA,CD(OP)%XVEC)		
	  ALLOCATE(CD(OP)%DATA(NX),CD(OP)%XVEC(NX))
	  CD(OP)%XVEC=CD(IP)%XVEC
	  CD(OP)%DATA=YV
	  NPTS(OP)=NX
	  CD(OP)%CURVE_ID=' '
	  IF(OP .GT. NPLTS)NPLTS=NPLTS+1
	END IF
!
	RETURN
	END
