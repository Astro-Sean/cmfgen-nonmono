!
! Read in AUTOIONIZATION strengths from a file. These values are tabulated for
! levels above the ionization limit.
!
	SUBROUTINE RD_AUTO_V1(AUTO,FEDGE,STAT_WT,LEVNAME,NLEV,FILE_NAME)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Aletered 23-Sep-2005 : Bug fix -- AUTO not correctly normalized by STAT_WT
!                                     in the case on unsplit levels.
! Created  17-Jul-2001
!
	INTEGER NLEV
	REAL(KIND=LDP) AUTO(NLEV)
	REAL(KIND=LDP) FEDGE(NLEV)
	REAL(KIND=LDP) STAT_WT(NLEV)
	CHARACTER*(*) LEVNAME(NLEV)
	CHARACTER*(*) FILE_NAME
!
! External functions
!
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
!
! Local variables and arrays.
!
	REAL(KIND=LDP) G_SUM(NLEV)
	REAL(KIND=LDP) TMP_AUTO
	REAL(KIND=LDP) DEFAULT_AUTO_RATE
	REAL(KIND=LDP) T1
!
	INTEGER I,J,K
	INTEGER NL
	INTEGER IBEG
	INTEGER ENTRY_NUM
	INTEGER NUM_OF_LEVELS_RD
	INTEGER IOS
	INTEGER, PARAMETER :: LUIN=7
!
	CHARACTER*30 TMP_NAME
	CHARACTER*30 LOCNAME(NLEV)
	CHARACTER*300 STRING
!
	LUER=ERROR_LU()
!
	LOCNAME(1:NLEV)=LEVNAME(1:NLEV)
	DO NL=1,NLEV
	  K=INDEX(LOCNAME(NL),'[')
	  IF(K .NE. 0)LOCNAME(NL)(K:)=' '
	END DO
	AUTO(1:NLEV)=0.0_LDP
	TMP_AUTO=0.0_LDP
	G_SUM(1:NLEV)=0.0_LDP
!
! Find first autoionizing level.
!
	IBEG=1
	T1=1000.0_LDP*2.998E+10_LDP/1.0E+15_LDP
	DO I=1,NLEV
	  IF(FEDGE(I) .LT. T1)THEN
	    IBEG=I
	    EXIT
	  END IF
	END DO
	IF(IBEG .EQ. 0)RETURN
!
	OPEN(UNIT=LUIN,FILE=FILE_NAME,STATUS='OLD',IOSTAT=IOS,ACTION='READ')
	IF(IOS .NE. 0)THEN
          WRITE(LUER,*)'Error opening ',FILE_NAME,' in RD_AUTO_V1'
	  CLOSE(LUIN)
          STOP
        END IF
!
	STRING=' '
	DO WHILE( INDEX(STRING,'!Format date')  .EQ. 0 )
	  READ(LUIN,'(A)',IOSTAT=IOS)STRING
	  IF(IOS .NE. 0)THEN
            WRITE(LUER,*)'Error format date in RD_AUTO_V1'
	    WRITE(LUER,*)'!Format date string not found'
	    WRITE(LUER,*)'IOS=',IOS
	    CLOSE(LUIN)
	    STOP
	  END IF
	END DO
!
	STRING=' '
	DO WHILE( INDEX(STRING,'!Number of energy levels')  .EQ. 0 )
	  READ(LUIN,'(A)',IOSTAT=IOS)STRING
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error reading ',FILE_NAME
	    WRITE(LUER,*)'!Number of energy levels string not found'
	    WRITE(LUER,*)'IOS=',IOS
	    CLOSE(LUIN)
	    STOP
	  END IF
	END DO
	READ(STRING,*)NUM_OF_LEVELS_RD
!
	STRING=' '
        DO WHILE( INDEX(STRING,'!Default autoionization rate')  .EQ. 0)
	  READ(LUIN,'(A)',IOSTAT=IOS)STRING
	  IF(IOS .NE. 0)THEN
            WRITE(LUER,*)'Error in RD_AUTO_V1'
            WRITE(LUER,*)'Currently treating ',FILE_NAME
            WRITE(LUER,*)'Can''t read default autoionization rate'
	    CLOSE(LUIN)
            STOP
          END IF
	END DO
        READ(STRING,*)DEFAULT_AUTO_RATE
!
	STRING=' '
        DO WHILE( INDEX(STRING,'!Entry number for autoionization rates') .EQ. 0)
          READ(LUIN,'(A)',IOSTAT=IOS)STRING
          IF(IOS .NE. 0)THEN
            WRITE(LUER,*)'Error in RD_AUTO_V1'
            WRITE(LUER,*)'Currently treating ',FILE_NAME
            WRITE(LUER,*)'Error reading entry number for level link'
	    CLOSE(LUIN)
            STOP
	  END IF
	END DO
	READ(STRING,*)ENTRY_NUM
        IF(ENTRY_NUM .LT. 2 .OR. ENTRY_NUM .GT. 30)THEN
          WRITE(LUER,*)'Error in RD_F_TO_S_IDS'
          WRITE(LUER,*)'Currently treating ',FILE_NAME
          WRITE(LUER,*)'Bad entry number for level link'
	  CLOSE(LUIN)
          STOP
        END IF
!
	DO NL=1,NUM_OF_LEVELS_RD
!
! Skip all blank lines and comments.
!
	  STRING=' '
	  DO WHILE(STRING .EQ. ' ' .OR. STRING(1:1) .EQ. '!')
	     READ(LUIN,'(A)')STRING
	  END DO
          STRING=ADJUSTL(STRING)
!
! Values are assumed to be separated by at LEAST 2 spaces.
!
	  K=INDEX(STRING,'  ')
	  TMP_NAME=STRING(1:K-1)
	  DO J=1,ENTRY_NUM-1
	    K=INDEX(STRING,'  ')
	    STRING=STRING(K+2:)
	    STRING=ADJUSTL(STRING)
	  END DO
	  READ(STRING,*)TMP_AUTO
!
! Now find matching levels.
!
	  K=INDEX(TMP_NAME,'[')
	  IF(K  .NE. 0)THEN
	    DO I=IBEG,NLEV
	      IF(TMP_NAME .EQ. LEVNAME(I))THEN
	        AUTO(I)=TMP_AUTO
	        EXIT
	      END IF
	      IF(TMP_NAME(1:K-1) .EQ. LEVNAME(I))THEN
	        AUTO(I)=AUTO(I)+STAT_WT(I)*TMP_AUTO
                G_SUM(I)=G_SUM(I)+STAT_WT(I)
	      END IF
	    END DO
	  ELSE
	    DO I=IBEG,NLEV
	      IF(TMP_NAME .EQ. LOCNAME(I) .OR. TMP_NAME .EQ. LEVNAME(I))THEN
	        AUTO(I)=TMP_AUTO
	      END IF
	    END DO
	  END IF
	END DO            !Number of levels read in.
	CLOSE(LUIN)
!
	DO I=IBEG,NLEV
	  IF(G_SUM(I) .NE. 0.0_LDP)AUTO(I)=AUTO(I)/G_SUM(I)
	END DO
!
	J=INDEX(FILE_NAME,'_')-1
	IF(J .LE. 1)J=LEN_TRIM(FILE_NAME)
	OPEN(LUIN,FILE='AUTO_CHK_'//FILE_NAME(1:J),STATUS='UNKNOWN')
	  WRITE(LUIN,*)'Utilized autoionization rates associated with FILE',FILE_NAME
	  DO I=IBEG,NLEV
	    IF(AUTO(I) .LE. 1.0E-20_LDP .AND. FEDGE(I) .LT. 0.0_LDP)AUTO(I)=DEFAULT_AUTO_RATE
	    WRITE(LUIN,'(1X,I4,3X,A,T35,ES12.4)')I,TRIM(LEVNAME(I)),AUTO(I)
	  END DO
	CLOSE(LUIN)
!
	RETURN
	END
