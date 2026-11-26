	SUBROUTINE RD_EW_IDS(XPAR,ID_FILE_NAME,LU_IN,T_OUT)
	USE SET_KIND_MODULE
        USE LINE_ID_MOD
	IMPLICIT NONE
!
! Altered 01-Sep-2019 : Changed to read new-format EW file.
!
	REAL*4 XPAR(2)
	INTEGER LU_IN
	INTEGER T_OUT
	CHARACTER(LEN=*) ID_FILE_NAME
!
	REAL(KIND=LDP) LINE_FLUX
	REAL(KIND=LDP) DP_T1
	REAL(KIND=LDP) T1,T2
	INTEGER I,J
	INTEGER IOS
	CHARACTER(LEN=200) TMP_STR
	LOGICAL VAC_WAVELENGTH
	LOGICAL LOG_VAR
	LOGICAL FILE_OPEN
!
        REAL(KIND=LDP) LAM_AIR
        EXTERNAL LAM_AIR
	LOGICAL, SAVE :: FIRST_TIME=.TRUE.
!
	IF(FIRST_TIME)THEN
	  FIRST_TIME=.FALSE.
	  ID_VEC_BEG=0.7
	  ID_VEC_END=0.5
	ELSE
	  WRITE(6,'(A,2ES14.4)')'Current values of ID_VEC_BEG and ID_VEC_END are:',ID_VEC_BEG,ID_VEC_END
	  WRITE(6,'(A,2ES14.4)')'            Good values for rectified plots are:',0.7,0.5
	END IF
!
	VAC_WAVELENGTH=.FALSE.
	OPEN(UNIT=LU_IN,FILE=TRIM(ID_FILE_NAME),STATUS='OLD',IOSTAT=IOS)
	IF(IOS .EQ. 0)THEN
	    J=0
	    TMP_STR='!'
	    DO WHILE(TMP_STR(1:1) .EQ. '!' .OR. INDEX(TMP_STR,'Lam(Ang)') .NE. 0)
	      READ(LU_IN,'(A)')TMP_STR
	      IF(INDEX(TMP_STR,'Alam(Ang)') .NE. 0)VAC_WAVELENGTH=.TRUE.
	    END DO
	    BACKSPACE(LU_IN)
!
! Old format file contained AIR wavelengths (only) for Lambda > 2000A.
! In the first read, T1 is used to store the air wavelength.
! T2 is used to store for the classic EW.
!
	    DO WHILE(J+1 .LE. N_LINE_ID_MAX)
	      IF(VAC_WAVELENGTH)THEN
	        READ(LU_IN,*,END=1500)ID_WAVE(J+1),T1,ID_CONT_FLUX(J+1),T2,ID_EW(J+1),LINE_FLUX,LOG_VAR,I,I,TMP_STR
	      ELSE
	        READ(LU_IN,*,END=1500)ID_WAVE(J+1),ID_CONT_FLUX(J+1),ID_EW(J+1),LINE_FLUX,LOG_VAR,I,I,TMP_STR
	        T1=ID_WAVE(J+1)
	        ID_WAVE(J+1)=LAM_AIR(T1)
	      END IF
	      OBSERVED_WAVE(I)=.TRUE.
	      I=INDEX(TMP_STR,'(')
	      LINE_ID(J+1)=TMP_STR(1:I-1)
	      FULL_LINE_ID=TMP_STR
	      IF( (ID_WAVE(J+1)-XPAR(1))*(XPAR(2)-ID_WAVE(J+1)) .GT. 0 .AND. ABS(ID_EW(J+1)) .GT. EW_CUT)THEN
	        J=J+1
		N_EW_IDS=J
	        ID_WAVE_OFF(J)=ID_WAVE(J)
	        ID_Y_BEG(J)=0.9D0
	        ID_Y_END(J)=0.8D0
	        IF(LINE_ID(J)(2:2) .EQ. 'k')LINE_ID(J)(2:2)='i'
	        IF(LINE_ID(J)(2:2) .EQ. '2')LINE_ID(J)(2:)='II'//LINE_ID(J)(3:)
	        IF(LINE_ID(J)(3:3) .EQ. '2')LINE_ID(J)(3:)='II'//LINE_ID(J)(4:)
	        IF(LINE_ID(J)(2:5) .EQ. 'XSIX')LINE_ID(J)(2:)='XVI'//LINE_ID(J)(6:)
	        IF(LINE_ID(J)(2:5) .EQ. 'XSEV')LINE_ID(J)(2:)='XVII'//LINE_ID(J)(6:)
	        IF(LINE_ID(J)(3:6) .EQ. 'XSIX')LINE_ID(J)(3:)='XVI'//LINE_ID(J)(7:)
	        IF(LINE_ID(J)(3:6) .EQ. 'XSEV')LINE_ID(J)(3:)='XVII'//LINE_ID(J)(7:)
	        IF(LINE_ID(J)(2:4) .EQ. 'SIX')LINE_ID(J)(2:)='VI'//LINE_ID(J)(5:)
	        IF(LINE_ID(J)(2:4) .EQ. 'SEV')LINE_ID(J)(2:)='VII'//LINE_ID(J)(5:)
	        IF(LINE_ID(J)(3:5) .EQ. 'SIX')LINE_ID(J)(3:)='VI'//LINE_ID(J)(6:)
	        IF(LINE_ID(J)(3:5) .EQ. 'SEV')LINE_ID(J)(3:)='VII'//LINE_ID(J)(6:)
	        IF(LINE_ID(J)(2:2) .GE. 'a' .AND. LINE_ID(J)(2:2) .LE.  'z')THEN
	          LINE_ID(J)(3:)=' '//LINE_ID(J)(3:)
	        ELSE
	          LINE_ID(J)(2:)=' '//LINE_ID(J)(2:)
	        END IF
	      END IF
	    END DO
	ELSE
	    WRITE(T_OUT,*)'Unable top open file'
	    N_EW_IDS=0
	END IF
1500	CONTINUE
	WRITE(6,*)'Number of EW line IDs is',N_EW_IDS
	INQUIRE(UNIT=LU_IN,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(UNIT=LU_IN)
	N_LINE_IDS=N_EW_IDS
!
	RETURN
	END
