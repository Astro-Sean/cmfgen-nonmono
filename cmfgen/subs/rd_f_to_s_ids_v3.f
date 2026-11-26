C
C Routine to read in the vector describing the matching of actual atomic levels
C to that in the model atom with super levels.
C
	SUBROUTINE RD_F_TO_S_IDS_V3(F_TO_S,INT_SEQ,LEVNAME_F,N_F,N_S,
	1                LUIN,FILENAME,SL_OPTION,F_TO_S_RD_ERROR)
	USE SET_KIND_MODULE
	IMPLICIT NONE
C
C Altered 21-Apr-2023 : Changed to V3. Added F_TO_S_RD_ERROR. Routine no longer stops
C                          on error. Better error message output to F_TO_S_RD_ERROS.
C Altered 04-Feb-2010 : Handles extra blank lines/comments mixed with levels.
C Altered 30-Sep-2005 : Removed warning about mixed parity.
C Altered 31-Jan-1997 : Only a warning is output if SUPER-LEVELS are defined
C                          with mixed parity.
C                       Improved error reporting for bad super-level ID.
C Altered 29-Dec-1996 : Dynamic allocation used for CNT, PAR
C Altered 29-May-1996 : String increased to 132 to accommodate longer names.
C Altered 26-May-1996 : ACTION installed in OPEN statement.
C Altered  2-Jan-1996 : INT_SEQ variable installed.
C Altered 13-Oct-1995 : If N_F=N_S 1 to 1 mapping of levels is assumed.
C Created 31-Mar-1995
C
	INTEGER N_S,N_F
	INTEGER LUIN
	INTEGER F_TO_S(N_F),INT_SEQ(N_F)
	CHARACTER(LEN=*) SL_OPTION
	CHARACTER(LEN=*) LEVNAME_F(N_F),FILENAME
	CHARACTER*132 STRING
C
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
	LOGICAL F_TO_S_RD_ERROR
C
C Local variables.
C
	INTEGER LUER,I,J,K,IOS,ENTRY_NUM
	INTEGER CNT(N_S),PAR(N_S)
	LOGICAL, SAVE :: FIRST_TIME=.TRUE.
C
C IF N_F .EQ. N_S there must be 1 to 1 correspondence between the levels.
C Thus we do not need to open a link file.
C
	IF(N_F .EQ. N_S)THEN
	  DO I=1,N_F
	    F_TO_S(I)=I
	    INT_SEQ(I)=0
	  END DO
	  RETURN
	END IF
C
	OPEN(UNIT=LUIN,FILE=FILENAME,STATUS='OLD',IOSTAT=IOS,
	1       ACTION='READ')
	IF(IOS .NE. 0)THEN
	  CALL OPEN_F_TO_S_ERROR_FILE()
	  WRITE(LUER,*)'Error opening ',FILENAME,' in RD_F_TO_S_IDS_V3'
	  F_TO_S_RD_ERROR=.TRUE.; FLUSH(LUER)
	  RETURN
	END IF
C
	STRING=' '
	DO WHILE( INDEX(STRING,'Number of energy levels') .EQ. 0)
	  READ(LUIN,'(A)')STRING
	END DO
	BACKSPACE(LUIN)
	READ(LUIN,*)I
	IF(I .LT. N_F)THEN
	  CALL OPEN_F_TO_S_ERROR_FILE()
	  WRITE(LUER,*)'Error in RD_F_TO_S_IDS_V3'
	  WRITE(LUER,*)'Currently treating ',FILENAME
	  WRITE(LUER,*)'Insufficient levels in file'
	  F_TO_S_RD_ERROR=.TRUE.; FLUSH(LUER)
	  CLOSE(LUIN); RETURN
	END IF
C
	STRING=' '
	DO WHILE( INDEX(STRING,'Entry number of link') .EQ. 0)
	  READ(LUIN,'(A)')STRING
	END DO
	BACKSPACE(LUIN)
	READ(LUIN,*)ENTRY_NUM
	IF(ENTRY_NUM .LT. 2 .OR. ENTRY_NUM .GT. 20)THEN
	  CALL OPEN_F_TO_S_ERROR_FILE()
	  WRITE(LUER,*)'Error in RD_F_TO_S_IDS_V3'
	  WRITE(LUER,*)'Currently treating ',FILENAME
	  WRITE(LUER,*)'Bad entry number for level link'
	  F_TO_S_RD_ERROR=.TRUE.; FLUSH(LUER)
	  CLOSE(LUIN); RETURN
	END IF
C
C NB: All entries must be separated by at LEAST 2 spaces.
C
	READ(LUIN,'(A)')STRING				!Blankline
	DO I=1,N_F
	  F_TO_S(I)=0.
	  STRING=' '
	  DO WHILE (STRING .EQ . ' ' .OR. STRING(1:1) .EQ. '!')
	    READ(LUIN,'(A)')STRING
	  END DO
	  STRING=ADJUSTL(STRING)			!Strip leading blanks.
	  J=INDEX(STRING,'  ')
	  IF(STRING(1:J) .NE. LEVNAME_F(I))THEN
	    CALL OPEN_F_TO_S_ERROR_FILE()
	    WRITE(LUER,*)'Error in RD_F_TO_S_IDS_V3'
	    WRITE(LUER,*)'Currently treating ',FILENAME
	    WRITE(LUER,*)'Invalid level id'
	    WRITE(LUER,'(A,I6,X,A)')'Level is ',I,TRIM(LEVNAME_F(I))
	    WRITE(LUER,*)'F_TO_S string is',STRING(1:J)
	    F_TO_S_RD_ERROR=.TRUE.; FLUSH(LUER)
	    CLOSE(LUIN); RETURN
	  END IF
C
C NB: ENTRY_NUM-2 as
C              -1 due to 1st entry being level ID
C              -1 as only need to point to entry.
C
	  DO K=1,ENTRY_NUM-2
	    STRING(1:)=STRING(J:)
	    DO WHILE(STRING(1:1) .EQ. ' ')
	      STRING(1:)=STRING(2:)
	    END DO
	    J=INDEX(STRING,'  ')
	  END DO
	  READ(STRING(J:),*)F_TO_S(I),INT_SEQ(I)
	END DO
	CLOSE(LUIN)
!
! The ' ' in the call prevents diagnostic information from being output.
! Nominally it should be the F_TO_S FILENME.
!
        IF(INDEX(SL_OPTION,'SPLIT_LS') .NE. 0)THEN
          CALL DO_SL_LS_SPLIT(F_TO_S,INT_SEQ,N_F,N_S,LEVNAME_F,SL_OPTION)
	ELSE
	  CALL DO_SL_ADJUSTEMENT(F_TO_S,INT_SEQ,N_F,N_S,SL_OPTION,' ')
	END IF
!
	IF(MAXVAL(F_TO_S) .NE. N_S)THEN
	  CALL OPEN_F_TO_S_ERROR_FILE()
	  WRITE(LUER,*)' '
	  WRITE(LUER,'(1X,80A)')('*',I=1,80)
	  WRITE(LUER,*)'Error in RD_F_TO_S_IDS_V3: Currently treating ',TRIM(FILENAME)
	  WRITE(LUER,'(A,I4)')' NF read in is',N_F
	  WRITE(LUER,'(A,I4)')' NS read in is',N_S
	  WRITE(LUER,'(A,I4)')' Maximum NS in F_TO_S corresponding to requested NF is',MAXVAL(F_TO_S)
	  WRITE(LUER,'(A,3X,A)')' Highest level to be treated is ',TRIM(LEVNAME_F(N_F))
	  WRITE(LUER,'(1X,80A)')('*',I=1,80)
	  WRITE(LUER,*)' '
	  F_TO_S_RD_ERROR=.TRUE.; FLUSH(LUER)
	  RETURN
	END IF
!
C Now check that:
C (1) All super-levels have at least 1 level from full atom
C (2) Correct number of super levels.
C (3) Parity of super levels matched (Warning only).
C
	DO I=1,N_F
	  IF(F_TO_S(I) .LE. 0 .OR. F_TO_S(I) .GT. N_S)THEN
	    CALL OPEN_F_TO_S_ERROR_FILE()
	    WRITE(LUER,*)'Error in RD_F_TO_S_IDS_V3'
	    WRITE(LUER,*)'Currently treating ',FILENAME
	    WRITE(LUER,*)'Incorrect super level ID for level ',LEVNAME_F(I)
	    F_TO_S_RD_ERROR=.TRUE.; FLUSH(LUER)
	    RETURN
	  END IF
	END DO
C
C Check that all super levels have at least 1 level.
C
	DO I=1,N_S
	  CNT(I)=0
	END DO
	DO J=1,N_F
	  I=F_TO_S(J)
	  CNT(I)=CNT(I)+1
	END DO
	DO I=1,N_S
	  IF(CNT(I) .EQ. 0)THEN
	    CALL OPEN_F_TO_S_ERROR_FILE()
	    WRITE(LUER,*)'Error in RD_F_TO_S_IDS_V3'
	    WRITE(LUER,*)'Currently treating ',FILENAME
	    WRITE(LUER,*)'Super level ',I,' has no components'
	    F_TO_S_RD_ERROR=.TRUE.; FLUSH(LUER)
	    RETURN
	  END IF
	END DO
C
C Check parity.
C
	DO I=1,N_S
	  PAR(I)=0
	END DO
	DO J=1,N_F
	  I=F_TO_S(J)
	  IF(INDEX(LEVNAME_F(J),'e') .NE. 0)THEN
	     PAR(I)=PAR(I)+1
	  ELSE
	     PAR(I)=PAR(I)-1
	  END IF
	END DO
!
!	J=0
!	DO I=1,N_S
!	  PAR(I)=ABS(PAR(I))/CNT(I)
!	  IF(PAR(I) .NE. 1)THEN
!	    IF(J .EQ. 0)THEN
!	      WRITE(LUER,*)'Warning in RD_F_TO_S_IDS_V3'
!	      WRITE(LUER,*)'Currently treating ',FILENAME
!	    END IF
!	    J=J+1
!	    WRITE(LUER,'(4X,A,I5,A)')'Super level ',
!	1                 I,' has a mix of odd and even levels'
!	  END IF
!	END DO
!
	RETURN
!
	CONTAINS
	SUBROUTINE OPEN_F_TO_S_ERROR_FILE()
	USE SET_KIND_MODULE
	IF(FIRST_TIME)THEN
	  CALL GET_LU(LUER,'RD_F_TO_S_IDS_V3')
	  OPEN(UNIT=LUER,FILE='F_TO_S_RD_ERRORS',STATUS='UNKNOWN',ACTION='WRITE')
	  FIRST_TIME=.FALSE.
	END IF
	RETURN
	END SUBROUTINE OPEN_F_TO_S_ERROR_FILE
!
	END SUBROUTINE RD_F_TO_S_IDS_V3
