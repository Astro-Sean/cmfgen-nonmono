!
! Routine to compute limits on the intrinsic line profile:
!
!RETURNED:
!        VEC_STRT_FREQ - Frequency at which intrinsic line absorption profile
!                          is assumed to start. Profile limits are assumed
!                          to be symmetric about line center.
!        VEC_VDOP_MIN  - Minimum Doppler width for line
!	
!INPUT:
!        ED_IN     - Electron density (/cm^3) (Vector, length ND)
!        TEMP_IN   - Temperature (10^4 K)     (Vector, length ND)
!        VTURB_IN  - Turbulent velocity in km/s (function of depth - length ND).
!        CHIL      - Line opacity             (Vector, length ND)
!
!        PROF_TYPE - Indicate type of intrinsic line profile; Profile types
!                      curently implemented are:
!                                               DOPPLER
!                                               VOIGT
!                                               STARK (Hydrogenic species only)
!
!        AMASS_IN  - Atomic mass for Doppler profile in amu's.
!        Z_IN      - Ion charge.
!        NL        - Lower level of transition.
!        NUP       - Upper level of transition.
!        ND        - Nuber of (Ne,T) values profile is to be computed for.
!        GAM_RAD   - Natural (radiative) broadening parameter for Voigt profile
!        CAM_COL   - Collisional broadening parameter for Voigt profile
!                          (quadratic stark effect).
!
!	LIMIT_SET_BY_OPACITY - Determines the STRT_FREQ of the intrinsic
!                                 line absorption profile based on the
!                                 ratio of line to electron scattering opacity.
!       DOP_LIMIT  - Truncate the Doppler profile when the ratio of the LINE
!                      opacity to the Electron scattering opacity is DOP_LIMIT.
!                      Should be of order 10^{-3} or less.
!                      The profile CANNOT be truncated insiede 3.5 Doppler
!                      widths from line center
!       VOIGT_LIMIT  - Truncate the VOIGT profile when the ratio of the LINE
!                      opacity to the Electron scattering opacity is DOP_LIMIT.
!                      Should be of order 10^{-3} or less.
!                      The profile CANNOT be truncated inside 3.5 Doppler
!                      widths from line center
!
	SUBROUTINE SET_PROF_LIMITS_V2(VEC_STRT_FREQ,VEC_VDOP_MIN,
	1             CHIL,ED_IN,TEMP_IN,VTURB_IN,ND,
	1             PROF_TYPE,PROF_LIST_LOCATION,
	1             NU_ZERO,NL,NUP,SPECIES_IN,AMASS_IN,Z_IN,
	1             GAM_RAD,GAM_COL,VTURB_FIX,
	1             DOP_LIMIT,VOIGT_LIMIT,
	1             LIMIT_SET_BY_OPACITY)
	USE SET_KIND_MODULE
	USE MOD_STRK_LIST
!
! Altered 01-May-2003 : Minor bug fix LO_GAM_COL was being incorrectly set.
! Altered 27-Nov-2001 : Bug fix. Incorrect frequency limits for DOP_FIX option.
! Altered 06-Apr-2000 : LINE_TO_CONT_RATIO was zero at one depth, causing a divide
!                         by zero.
!
	IMPLICIT NONE
	INTEGER ND
	INTEGER NL
	INTEGER NUP
!
	REAL(KIND=LDP) VEC_STRT_FREQ
	REAL(KIND=LDP) VEC_VDOP_MIN
	REAL(KIND=LDP) DOP_LIMIT
	REAL(KIND=LDP) VOIGT_LIMIT
!
	REAL(KIND=LDP) CHIL(ND)
	REAL(KIND=LDP) ED_IN(ND)
	REAL(KIND=LDP) TEMP_IN(ND)
	REAL(KIND=LDP) VTURB_IN(ND)
!
	CHARACTER(LEN=*) SPECIES_IN
	REAL(KIND=LDP) AMASS_IN
	REAL(KIND=LDP) Z_IN
	REAL(KIND=LDP) NU_ZERO
	REAL(KIND=LDP) GAM_RAD
	REAL(KIND=LDP) GAM_COL
	REAL(KIND=LDP) VTURB_FIX
	INTEGER PROF_LIST_LOCATION
	CHARACTER*(*) PROF_TYPE
!
	LOGICAL LIMIT_SET_BY_OPACITY
!
	REAL(KIND=LDP) SPEED_OF_LIGHT
	INTEGER GET_INDX_DP
	EXTERNAL SPEED_OF_LIGHT,GET_INDX_DP
!
! Local variables
!
	INTEGER I,J
	REAL(KIND=LDP) ESEC(ND)
	REAL(KIND=LDP) TMP_VEC(ND)
	REAL(KIND=LDP) NU_DOP(ND)
	REAL(KIND=LDP) LINE_TO_CONT_RATIO(ND)
	REAL(KIND=LDP) dNU
	REAL(KIND=LDP) C_KMS
	REAL(KIND=LDP) WAVE
	REAL(KIND=LDP) LOC_GAM_RAD
	REAL(KIND=LDP) LOC_GAM_COL
	REAL(KIND=LDP) V_PROF_LIMIT
	REAL(KIND=LDP) PROF_LINE_CENTER
	REAL(KIND=LDP) T1,T2
!
	INTEGER, PARAMETER :: NUM_DOP=6
!
	C_KMS=1.0E-05_LDP*SPEED_OF_LIGHT()		!km/s
	V_PROF_LIMIT=3000.0_LDP			!km/s
	LOC_GAM_RAD=GAM_RAD
	LOC_GAM_COL=GAM_COL
	PROF_LIST_LOCATION=0
!
! The option assumes fixed width Doppler profiles, and recovers exactly the
! same option as was installed in CMFGEN prior to the installations of variable
! Doppler widths.
!
	IF(PROF_TYPE .EQ. 'DOP_FIX')THEN
	  VEC_VDOP_MIN=VTURB_FIX                                         !kms
          NU_DOP(1)=NU_ZERO*VTURB_FIX/C_KMS				!kms
          VEC_STRT_FREQ=NU_ZERO+NUM_DOP*NU_DOP(1)
	  RETURN
	END IF
!
	IF(PROF_TYPE(1:4) .EQ. 'LIST')THEN
	  WAVE=0.01_LDP*C_KMS/NU_ZERO
	  I=GET_INDX_DP(WAVE,LST_WAVE,N_LST)
	  DO J=MAX(1,I-2),MIN(I+2,N_LST)
	    IF( TRIM(SPECIES_IN) .EQ. TRIM(LST_SPECIES(J)) .AND.
	1          ABS((WAVE-LST_WAVE(J))/WAVE) .LT. 1.0E-05_LDP)THEN
	      PROF_TYPE=LST_TYPE(J)
	      PROF_LIST_LOCATION=J
	      IF(LST_V_PROF_LIMIT(J) .NE. 0)V_PROF_LIMIT=LST_V_PROF_LIMIT(J)
	      IF(PROF_TYPE .EQ. 'VOIGT')THEN
	        LOC_GAM_RAD=LST_GAM_RAD(J)
	        LOC_GAM_COL=LST_GAM_COL(J)
	      END IF
	      EXIT		!As found profile.
	    END IF
	  END DO
	END IF
!
! We assume all lines are dominated by the Doppler profile at line center.
! We ensure LINE_TO_CONT ratio is not zero, to prevent division by zero.
!
        TMP_VEC(1:ND)=( VTURB_IN(1:ND)/12.85_LDP  )**2
	ESEC(1:ND)=6.65E-15_LDP*ED_IN(1:ND)
!
        T1=1.0E-15_LDP/1.77245385095516_LDP         !1.0D-15/SQRT(PI)
	VEC_VDOP_MIN=1.0E+50_LDP			!Very large number
        DO I=1,ND
          NU_DOP(I)=12.85_LDP*SQRT( TEMP_IN(I)/AMASS_IN + TMP_VEC(I) )	!kms
	  VEC_VDOP_MIN=MIN(VEC_VDOP_MIN,NU_DOP(I))			!kms
          NU_DOP(I)=NU_DOP(I)*NU_ZERO/C_KMS				!10^15 Hz
          PROF_LINE_CENTER=T1/NU_DOP(I)
	   LINE_TO_CONT_RATIO(I)=ABS(CHIL(I))*PROF_LINE_CENTER/ESEC(I)
	   IF(LINE_TO_CONT_RATIO(I) .EQ. 0)LINE_TO_CONT_RATIO(I)=1.0E-50_LDP
	   IF(LINE_TO_CONT_RATIO(I) .GT. 1.0E+04_LDP .AND. PROF_TYPE .EQ. 'LIST_VGT')THEN
	     LOC_GAM_RAD=GAM_RAD
	     LOC_GAM_COL=0.0_LDP		!GAM_COL
	     PROF_TYPE='VOIGT'
	   END IF
	END DO
!
! If the line is not in the list, chooses an alternative profile.
!
	IF(PROF_TYPE(1:4) .EQ. 'LIST')THEN
	  PROF_TYPE='DOPPLER'
	  IF(SPECIES_IN .EQ. 'HI' .OR. SPECIES_IN .EQ. 'He2')THEN
	    PROF_TYPE='HZ_STARK'
	  END IF
	END IF
!
	VEC_STRT_FREQ=0.0_LDP
	IF(PROF_TYPE .EQ. 'DOPPLER' .OR. PROF_TYPE .EQ. 'VOIGT')THEN
	  IF(LIMIT_SET_BY_OPACITY)THEN
	    DO I=1,ND
              T1=-LOG(DOP_LIMIT/LINE_TO_CONT_RATIO(I))
	      IF(T1 .GT. 0)T1=SQRT(T1)
	      T1=MAX(3.5_LDP,T1)
              VEC_STRT_FREQ=MAX(VEC_STRT_FREQ,NU_ZERO+T1*NU_DOP(I))
	    END DO
	  ELSE
	    DO I=1,ND
              VEC_STRT_FREQ=MAX(VEC_STRT_FREQ,NU_ZERO+NUM_DOP*NU_DOP(I))
	    END DO
	  END IF
	  IF(PROF_TYPE .EQ. 'DOPPLER')RETURN
	END IF
!
! NB: 6.7005D-09=1.0D-15*SQRT(1.0D+15 /4  / PI**1.5)
! The factor of 1.0D-15 arises since dNU has to be in units of 10^15 Hz.
! The factor of 1.0D+15 arisis since NU_DOP is in units of 10^15 Hz.
!
	dNU=0.0_LDP
	IF(PROF_TYPE .EQ. 'VOIGT')THEN
	  IF(LIMIT_SET_BY_OPACITY)THEN
	    DO I=1,ND
	      T2=LOC_GAM_RAD+LOC_GAM_COL*ED_IN(I)
              T1=6.7005E-09_LDP*SQRT(T2*NU_DOP(I)*LINE_TO_CONT_RATIO(I)/VOIGT_LIMIT)
	      dNU=MAX(T1,dNU)
	    END DO
	    dNU=MIN(dNU,NU_ZERO*5.0E+03_LDP/C_KMS)
	    VEC_STRT_FREQ=MAX(VEC_STRT_FREQ,NU_ZERO+dNU)
	  ELSE
!
! Truncates profile when profile is down to 1.0E-06 of the value at line
! center.
!
	    DO I=1,ND
	      T2=LOC_GAM_RAD+LOC_GAM_COL*ED_IN(I)
              T1=6.7005E-09_LDP*SQRT(1.0E+06_LDP*T2*NU_DOP(I))
	      dNU=MAX(T1,dNU)
	    END DO
	    dNU=MIN(dNU,NU_ZERO*5.0E+03_LDP/C_KMS)
	    VEC_STRT_FREQ=MAX(VEC_STRT_FREQ,NU_ZERO+dNU)
	  END IF
	  RETURN
	END IF
!
! For all other profile types: CRUDE.
!
	VEC_STRT_FREQ=NU_ZERO*(1.0_LDP+V_PROF_LIMIT/C_KMS)
!
!	T1=(VEC_STRT_FREQ/NU_ZERO-1.0D0)*C_KMS
!	WRITE(177,'(X,A10,2E14.5)')PROF_TYPE,NU_ZERO,T1
!
	RETURN
	END
