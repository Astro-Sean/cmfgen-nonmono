!
! Simple subroutine to read in different abundance data sets. This will
! allow the user to define the reference data set used in the summary file,
! MOD_SUM.
!
! NB: These abudnaces DO not affect the CMFGEN calcualtions -- only MOD_SUM.
!
	SUBROUTINE RD_SOL_ABUND_SCAL(SOL_MASS_FRAC,SOL_ABUND_HSCL,AT_MASS,AT_NO,
	1                              SOL_ABUND_REF_SET,MUST_GET_SET,NSPEC)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Finalized: 19-Jun-2023
!
	INTEGER NSPEC
	REAL(KIND=LDP) AT_MASS(NSPEC)
	REAL(KIND=LDP) AT_NO(NSPEC)
	REAL(KIND=LDP) SOL_ABUND_HSCL(NSPEC)
	REAL(KIND=LDP) SOL_MASS_FRAC(NSPEC)
	LOGICAL MUST_GET_SET
	CHARACTER(LEN=*) SOL_ABUND_REF_SET
!
! Local variabls.
!
	INTEGER, PARAMETER :: MAX_N_ABUND=10
	REAL(KIND=LDP) ABUND(MAX_N_ABUND)
	INTEGER COL(MAX_N_ABUND)
	CHARACTER(LEN=10) REF_SET(MAX_N_ABUND)
!
	REAL(KIND=LDP) T1
	REAL(KIND=LDP) LOCAL_AT_NO
	REAL(KIND=LDP) LOC_AT_MASS
	INTEGER LUIN
	INTEGER I,J,K
	INTEGER INDX
	INTEGER ST_AB_COL
	INTEGER N_ABUND
	INTEGER IOS
!
	LOGICAL ERROR
	CHARACTER(LEN=2) SYMB
	CHARACTER(LEN=20) ELEMENT
	CHARACTER(LEN=200) STRING
!
! If MUST_GET_SET is FALSE, we assume that SOL_ABUND_HSCL must have
! already been set.
!
	IF(MUST_GET_SET)SOL_ABUND_HSCL=0.0_LDP
!
	CALL GET_LU(LUIN,'In rd_sol_abund_scale')
	OPEN(UNIT=LUIN,FILE='SOL_ABUNDANCE',STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	IF(IOS .NE. 0)THEN
	  IF(MUST_GET_SET)THEN
	    WRITE(6,*)'Error -- file SOL_ABUNDANCE not found.'
	    WRITE(6,*)'File must be present as abundance scale must get set'
	    STOP
	  ELSE
	    WRITE(6,*)'Warning -- file SOL_ABUNDANCE not found.'
	    WRITE(6,*)'Using default solar abundance scale.'
	    SOL_ABUND_REF_SET='Cox2000'
	    RETURN
	  END IF
	END IF
!
	DO WHILE(INDEX(STRING,'Number of abundance columns') .EQ. 0)
	  READ(LUIN,'(A)')STRING
	END DO
	READ(STRING,*)N_ABUND
!
	READ(LUIN,'(A)')STRING
	IF(INDEX(STRING,'Format date') .EQ. 0)THEN
	   SOL_ABUND_REF_SET='Cox2000'
	   IF(MUST_GET_SET)THEN
	     WRITE(6,*)'Warning format date not found.'
	     WRITE(6,*)'Format date must be present in file.'
	     STOP
	   ELSE
	     WRITE(6,*)'Error format date not found'
	     WRITE(6,*)'Using default solar abundance scale'
	     SOL_ABUND_REF_SET='Cox2000'
	     RETURN
	   END IF
	END IF
!
! Read in reference key for each abundance set.
!
	REF_SET=' '; COL=0
	DO K=1,N_ABUND
	  STRING=' '
	  DO WHILE(STRING .EQ. ' ')
	    READ(LUIN,'(A)')STRING
	  END DO
	  I=INDEX(STRING,'Abundance reference:')
	  IF(I .NE. 0)THEN
	    I=INDEX(STRING,':')
	    STRING=STRING(I+1:)
	    STRING=ADJUSTL(STRING)
	    J=INDEX(STRING,' ')
	    REF_SET(K)=STRING(1:J)
	    READ(STRING(J:),*)COL(K)
	  ELSE
	    WRITE(6,*)'Abundance keys in SOL_ABUNDANCE not found.'
	    WRITE(6,*)'Available abundance keys are listed below.'
	    WRITE(6,'(A)')REF_SET(1:N_ABUND)
	    STOP
	  END IF
	END DO
!
	INDX=0
	DO K=1,N_ABUND
	  IF(REF_SET(K) .EQ. SOL_ABUND_REF_SET)THEN
	    INDX=COL(K)
	    EXIT
	  END IF
	END DO
	ST_AB_COL=MINVAL(COL(1:N_ABUND))
!
	IF(INDX .EQ. 0)THEN
	  WRITE(6,*)'SOL_ABUND_REF_SET not found'
	  WRITE(6,*)'Specified set is: ',TRIM(SOL_ABUND_REF_SET)
	  WRITE(6,'(10(A,2X))')' Available sets are:',(TRIM(REF_SET(K)),K=1,N_ABUND)
	  STOP
	END IF
!
	STRING=' '
	DO WHILE(INDEX(STRING,'Hydrogen') .EQ. 0)
          READ(LUIN,'(A)')STRING
	  WRITE(6,'(A)')TRIM(STRING)
	END DO
        BACKSPACE(LUIN)
!
	ABUND=0.0_LDP
	WRITE(6,*)COL
	WRITE(6,*)INDX
	DO WHILE(1 .EQ. 1)
          READ(LUIN,*,END=100)LOCAL_AT_NO,SYMB,ELEMENT,LOC_AT_MASS,
	1                        (ABUND(J),J=ST_AB_COL,N_ABUND+ST_AB_COL-1)
	  DO I=1,NSPEC
	    IF( NINT(AT_NO(I)) .EQ. NINT(LOCAL_AT_NO) )THEN
	      SOL_ABUND_HSCL(I)=ABUND(INDX)
	      EXIT
	    END IF
	  END DO
	END DO
!
100	CONTINUE
!
	ERROR=.FALSE.
	IF(MUST_GET_SET)THEN
	  DO I=1,NSPEC
	    IF(SOL_ABUND_HSCL(I) .EQ. 0.0_LDP .AND. AT_NO(I) .LE. 92)THEN
	      WRITE(6,*)'Solar abundance has not been set for atomic no',AT_NO(I)
	      ERROR=.TRUE.
	    END IF
	  END DO
	  IF(ERROR)STOP
	END IF
!
! Convert from the abundance on a logarithmic scale with H=12.0 dex, to
! mass-fractions.
!
	T1=0.0_LDP
	DO I=1,NSPEC
	  SOL_MASS_FRAC(I)=10.0_LDP**(SOL_ABUND_HSCL(I)-12.0_LDP)
	  SOL_MASS_FRAC(I)=AT_MASS(I)*SOL_MASS_FRAC(I)
	  T1=T1+SOL_MASS_FRAC(I)
	END DO
	SOL_MASS_FRAC(:)=SOL_MASS_FRAC(:)/T1
!
	RETURN
	END
