!
! Program searches the given supplied atomic data directiry to create a
! file containging soft data links for all species other than H and He.
!
! The latest data directory for each ionization stage is used.
!
! Some changes may be need to the number of super and full levels in MODEL_SPEC
!    if the new model atoms differ from the old model atom.
!    The F_TO_S link file may also need to be changed
!
	PROGRAM CREATE_BATCH_INS
	USE SET_KIND_MODULE
	USE GEN_IN_INTERFACE
	USE MOD_COLOR_PEN_DEF
	IMPLICIT NONE
!
! Created 22-Apr-2023
!
	INTEGER IS
	INTEGER ID
	INTEGER L
	INTEGER J
	INTEGER IOS
!
	INTEGER, PARAMETER :: NSPEC=23
	CHARACTER(LEN=5) SPECIES(NSPEC)
	CHARACTER(LEN=3) SPECIES_ABR(NSPEC)
!
	INTEGER, PARAMETER :: NROM=20
	CHARACTER(LEN=5) ROM(NROM)
	CHARACTER(LEN=5)   GEN_ION_ID(0:NROM)

!
	CHARACTER(LEN=1) DIR_DELIM
!
	INTEGER, PARAMETER :: NDATA=8
	CHARACTER(LEN=30) DATA_FILE(NDATA)
	CHARACTER(LEN=30) LNK_NAME(NDATA)
	CHARACTER(LEN=50) F_TO_S_FILE
!
	CHARACTER(LEN=80) LATEST_DIR
	CHARACTER(LEN=80) PATH
	CHARACTER(LEN=80) SHRT_PATH
	CHARACTER(LEN=80) ATOMIC_DIRECTORY
	CHARACTER(LEN=40) SPEC_DIR
	CHARACTER(LEN=40) LNK_ID
	CHARACTER(LEN=200) FILENAME
	CHARACTER(LEN=200) BATCH_FILE
!
	LOGICAL FILE_EXISTS
	LOGICAL LNK_WRITTEN
	INTEGER, PARAMETER :: LU_DIAG=10
	INTEGER, PARAMETER :: LU_BAT=12
!
!This is dependent on the operating system.
!
	DATA DIR_DELIM/'/'/
!
! The short species abreviatons as used by CMFGEN.
!
	DATA SPECIES_ABR/'C','N','O','Ne',
	1       'Na','Mg','Al','Sk','P','S','Cl','Ar',
	1       'K','Ca','Sc','V','Tk','Cr','Mn','Fe',
	1       'Co','Nk','Ba'/
!
! The specis abreviatons used in the main atomic data directory.
!
	DATA SPECIES/'CARB','NIT','OXY','NEON',
	1       'NA','MG','AL','SIL','PHOS','SUL','CHL','ARG',
	1       'POT','CA','SCAN','VAN','TIT','CHRO','MAN','FE',
	1       'COB','NICK','BAR'/
!
! Ionization abreviations used by CMFGEN.
!
	DATA ROM/'I','II','III','IV','V','VI','VII','VIII','IX','X',
	1    'XI','XII','XIII','XIV','XV','XVI','XVII','XVIII','XIX','ZX'/
!
! Ionization abreviations used in the atomic data directory.
!
	DATA GEN_ION_ID /'0','I','2','III','IV','V',
	1                'SIX','SEV','VIII','IX','X','XI','XII',
	1                'XIII','XIV','XV','XSIX','XSEV','X8','X9','XX'/
!
! Principal atomic data file names.
!
	DATA DATA_FILE/'col_data','osc_data','auto_data','die_data',
	1              'phot_data_A','phot_data_B','phot_data_C','phot_data_D'/
!
! Used to create soft data links for CMFGEN. Their odering must be consistent
!    with the atomic data file names given above.
!
	DATA LNK_NAME/'COL_DATA','F_OSCDAT','AUTO_DATA','DIE',
	1             'PHOT_A','PHOT_B','PHOT_C','PHOT_D'/
!
	WRITE(6,*)BLUE_PEN
	WRITE(6,*)' Program searches the given supplied atomic data directory to create a'
	WRITE(6,*)'     file containging soft data links for all species other than H and He.'
	WRITE(6,*)' '
	WRITE(6,*)' The latest data directory for each ionization stage is used.'
	WRITE(6,*)' '
	WRITE(6,*)' Some changes may be need to the number of super and full levels in MODEL_SPEC'
	WRITE(6,*)'       if the new model atoms differ from the old model atom.'
	WRITE(6,*)'       The F_TO_S link file may also need to be changed.'
	WRITE(6,*)' '
	WRITE(6,*)' See BATCH_WRITE_DIAG for program diagnostics.'
	WRITE(6,*)DEF_PEN
!
	ATOMIC_DIRECTORY=' '
	CALL GEN_IN(ATOMIC_DIRECTORY,'Atomic directory (full path)')
	INQUIRE(FILE=ATOMIC_DIRECTORY,EXIST=FILE_EXISTS)
	IF(.NOT. FILE_EXISTS)THEN
	  WRITE(6,*)'Error -- ATOMIC directory not found'
	  STOP
	END IF
	L=LEN_TRIM(ATOMIC_DIRECTORY)
	IF(ATOMIC_DIRECTORY(L:L) .NE. DIR_DELIM)THEN
	  ATOMIC_DIRECTORY=TRIM(ATOMIC_DIRECTORY)//DIR_DELIM
	END IF
!
	WRITE(6,*)RED_PEN
	WRITE(6,*)' If file is present, code tries to obtain the old F_TO_S link from file'
	WRITE(6,*)DEF_PEN
	BATCH_FILE=' '
	CALL GEN_IN(BATCH_FILE,'batch.sh file with old links (can be blank)')
!
	OPEN(LU_BAT,FILE='batch_ins.sh',STATUS='NEW',ACTION='WRITE',IOSTAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(6,*)RED_PEN
	  WRITE(6,*)' Unable to open batch_ins.sh IOS= ',IOS
	  WRITE(6,*)' If file already exists you must rename it, or delete it'
	  WRITE(6,*)DEF_PEN
	  STOP
	END IF
!
	OPEN(LU_DIAG,FILE='BATCH_WRITE_DIAG',STATUS='UNKNOWN',ACTION='WRITE',IOSTAT=IOS)
!
! Loop over all species.
!
	DO IS=1,NSPEC
	  FILENAME=TRIM(ATOMIC_DIRECTORY)//TRIM(SPECIES(IS))
	  INQUIRE(FILE=FILENAME,EXIST=FILE_EXISTS)
	  SPECIES(IS)=TRIM(SPECIES(IS))//DIR_DELIM
	  IF(.NOT. FILE_EXISTS)THEN
	    WRITE(6,*)'Error -- species directory not found'
	    WRITE(6,*)'Species is ',SPECIES(IS)
	    STOP
	  END IF
	  WRITE(LU_DIAG,*)'Doing ',TRIM(SPECIES(IS))
!
! Loop over ionization state. We use the latest directory in each ionization
! state for the atomic data. The date data directories must have the form
! 1jan23 or 12jan23 etc.
!
	  DO ID=1,NROM
	    F_TO_S_FILE=' '
	    LNK_WRITTEN=.FALSE.
	    PATH=TRIM(ATOMIC_DIRECTORY)//TRIM(SPECIES(IS))//TRIM(ROM(ID))
	    INQUIRE(FILE=PATH,EXIST=FILE_EXISTS)
	    IF(FILE_EXISTS)THEN
	      CALL GET_LAST_DATA_SET(PATH,LATEST_DIR)
	      WRITE(LU_DIAG,*)'  Latest data set for ',ROM(ID)(1:4),' is ',TRIM(LATEST_DIR)
	      PATH=TRIM(PATH)//DIR_DELIM//TRIM(LATEST_DIR)
	      SHRT_PATH=TRIM(SPECIES(IS))//TRIM(ROM(ID))//DIR_DELIM//TRIM(LATEST_DIR)//DIR_DELIM
	      DO J=1,NDATA
	        FILENAME=TRIM(PATH)//DIR_DELIM//TRIM(DATA_FILE(J))
	        INQUIRE(FILE=FILENAME,EXIST=FILE_EXISTS)
	        IF(FILE_EXISTS)THEN
	          IF(INDEX(LNK_NAME(J),'PHOT') .NE. 0)THEN
	             WRITE(LU_BAT,'(1X,A,A,T60,A)')'ln -sf  ',
	1              '$ATOMIC'//DIR_DELIM//TRIM(SHRT_PATH)//TRIM(DATA_FILE(J)),
	1              'PHOT'//TRIM(SPECIES_ABR(IS))//TRIM(GEN_ION_ID(ID))//LNK_NAME(J)(5:6)
	          ELSE IF(INDEX(LNK_NAME(J),'DIE') .NE. 0)THEN
	             WRITE(LU_BAT,'(1X,A,A,T60,A)')'ln -sf  ',
	1              '$ATOMIC'//DIR_DELIM//TRIM(SHRT_PATH)//TRIM(DATA_FILE(J)),
	1              'DIE'//TRIM(SPECIES_ABR(IS))//TRIM(GEN_ION_ID(ID))
	          ELSE
	             WRITE(LU_BAT,'(1X,A,A,T60,A)')'ln -sf  ',
	1              '$ATOMIC'//DIR_DELIM//TRIM(SHRT_PATH)//TRIM(DATA_FILE(J)),
	1              TRIM(SPECIES_ABR(IS))//TRIM(GEN_ION_ID(ID))//'_'//TRIM(LNK_NAME(J))
	          END IF
	          LNK_WRITTEN=.TRUE.
	        ELSE
	          IF(J .LT. 3 .OR. J .EQ. 5)THEN
	            WRITE(6,*)' The file was not found: '//TRIM(FILENAME)
	          END IF
	          WRITE(LU_DIAG,*)'     The file was not found: '//TRIM(FILENAME)
	        END IF
	      END DO
!
! Get the exact F_TO_S file if possible, otherwise the chosen F_TO_S file is random.
!
	      IF(LNK_WRITTEN .AND. BATCH_FILE .NE. ' ')THEN
	        LNK_ID=TRIM(SPECIES_ABR(IS))//TRIM(GEN_ION_ID(ID))//'_F_TO_S'
	        CALL GET_FS_FILE_FROM_BATCH(LNK_ID,F_TO_S_FILE,DIR_DELIM,BATCH_FILE)
	        IF(F_TO_S_FILE .NE. ' ')THEN
	          FILENAME=TRIM(PATH)//DIR_DELIM//TRIM(F_TO_S_FILE)
	          INQUIRE(FILE=FILENAME,EXIST=FILE_EXISTS)
	          IF(.NOT. FILE_EXISTS)THEN
	            WRITE(6,*)' Old F_TO_S link not found in new directory: '//TRIM(FILENAME)
	            WRITE(LU_DIAG,'(10X,A)')'Old F_TO_S link not found in new directory: '//TRIM(FILENAME)
	            F_TO_S_FILE=' '
	          END IF
	        END IF
	      END IF
	      IF(F_TO_S_FILE .EQ. ' ')THEN
	         CALL GET_F_TO_S_FILE(PATH,F_TO_S_FILE,LU_DIAG)
	      END IF
	      WRITE(LU_BAT,'(1X,A,A,T60,A,/)')'ln -sf  ',
	1            '$ATOMIC'//DIR_DELIM//TRIM(SHRT_PATH)//TRIM(F_TO_S_FILE),
	1             TRIM(SPECIES_ABR(IS))//TRIM(GEN_ION_ID(ID))//'_'//'F_TO_S'
	    END IF
	  END DO
	  IF(LNK_WRITTEN)WRITE(LU_DIAG,'(A)')' '
	END DO
!
	STOP
	END
!
	SUBROUTINE GET_LAST_DATA_SET(FILENAME,LATEST_DIR)
	USE SET_KIND_MODULE
	IMPLICIT NONE
	CHARACTER(LEN=*) FILENAME,LATEST_DIR
!
	CHARACTER(LEN=3) MONTH,MONTH_SAVE
	INTEGER YEAR,DAY
	INTEGER YEAR_SAVE,DAY_SAVE
	LOGICAL FILE_EXISTS
	CHARACTER(LEN=200) DIR_FILE
	CHARACTER(LEN=200) STRING
!
	YEAR_SAVE=0; MONTH_SAVE=' '; DAY_SAVE=0
!
	INQUIRE(FILE='DIRECTORY_LISTING',EXIST=FILE_EXISTS)
	IF(FILE_EXISTS)THEN
	  WRITE(6,*)'Error -- file DIRECTORY_LISTING exists'
	  WRITE(6,*)'Rename or delete'
	  STOP
	END IF
!
	STRING='ls '//TRIM(FILENAME)//' > DIRECTORY_LISTING'
	CALL SYSTEM(STRING)
	OPEN(UNIT=20,FILE='DIRECTORY_LISTING',ACTION='READ')
	DO WHILE(1 .EQ. 1)
	  READ(20,*,END=100)DIR_FILE
	  IF(DIR_FILE(1:1) .GE. '0' .AND. DIR_FILE(1:1) .LE. '9')THEN
	    IF(DIR_FILE(2:2) .GE. '0' .AND. DIR_FILE(2:2) .LE. '9')THEN
	      READ(DIR_FILE(1:2),*)DAY
	      MONTH=DIR_FILE(3:5)
	      READ(DIR_FILE(6:7),*)YEAR
	    ELSE
	      READ(DIR_FILE(1:1),*)DAY
	      MONTH=DIR_FILE(2:4)
	      READ(DIR_FILE(5:6),*)YEAR
	    END IF
	    IF(YEAR .LT. 80)YEAR=YEAR+2000
	    IF(YEAR .GT. YEAR_SAVE)THEN
	      YEAR_SAVE=YEAR
	      MONTH_SAVE=MONTH_SAVE
	      DAY_SAVE=DAY
	      LATEST_DIR=DIR_FILE
	    ELSE IF(YEAR .EQ. YEAR_SAVE)THEN
	      IF(MONTH .GT. MONTH_SAVE)THEN
	        MONTH=MONTH_SAVE
	        DAY=DAY_SAVE
	        LATEST_DIR=DIR_FILE
	      ELSE IF(MONTH .EQ. MONTH_SAVE)THEN
	         IF(DAY .GT. DAY_SAVE)THEN
	           DAY_SAVE=DAY
	           LATEST_DIR=DIR_FILE
	         END IF
	      END IF
	    END IF
	  END IF
	END DO
100	CONTINUE
	CALL SYSTEM('rm -f DIRECTORY_LISTING')
!
	RETURN
	END

	SUBROUTINE GET_FS_FILE_FROM_BATCH(LNK_ID,F_TO_S_FILE,DIR_DELIM,BATCH_FILE)
	USE SET_KIND_MODULE
	IMPLICIT NONE
	CHARACTER(LEN=*) LNK_ID
	CHARACTER(LEN=*) F_TO_S_FILE
	CHARACTER(LEN=1) DIR_DELIM
	CHARACTER(LEN=*) BATCH_FILE
!
	INTEGER K,L
	INTEGER LUIN
	INTEGER, SAVE :: NST
	LOGICAL, SAVE :: FIRST_TIME=.TRUE.
	CHARACTER(LEN=80), SAVE :: STORE(500)
	CHARACTER(LEN=15), SAVE :: LNK_STORE(500)
	CHARACTER(LEN=200) STRING
!
	IF(FIRST_TIME)THEN
	  CALL GET_LU(LUIN,' ')
	  NST=0
	  OPEN(UNIT=LUIN,FILE=TRIM(BATCH_FILE),STATUS='OLD',ACTION='READ')
	  DO WHILE(1 .EQ. 1)
	    READ(LUIN,'(A)',END=200)STRING
	    IF(INDEX(STRING,'F_TO_S') .NE. 0)THEN
	      NST=NST+1
!
	      K=LEN_TRIM(STRING)
	      DO WHILE(K .GT. 1)
	        K=K-1
	        IF(STRING(K:K) .EQ. ' ')THEN
	          LNK_STORE(NST)=STRING(K+1:)
	          STRING(K+1:)=' '
	          EXIT
	        END IF
	      END DO
!
	      DO WHILE(K .GT. 1)
	        K=K-1
	        IF(STRING(K:K) .EQ. DIR_DELIM)THEN
	          STORE(NST)=STRING(K+1:)
	          WRITE(50,'(A)')TRIM(STORE(NST))//'     '//TRIM(LNK_STORE(NST))
	          EXIT
	        END IF
	      END DO
	    END IF
	  END DO
200	  CONTINUE
	  CLOSE(UNIT=LUIN)
	  FIRST_TIME=.FALSE.
	END IF
!
	F_TO_S_FILE=' '
	DO K=1,NST
	  IF(TRIM(LNK_ID) .EQ. TRIM(LNK_STORE(K)))THEN
	    F_TO_S_FILE=STORE(K)
	    L=INDEX(F_TO_S_FILE,'.dat')
	    IF(L .NE. 0)F_TO_S_FILE(L:)=' '
	    L=INDEX(F_TO_S_FILE,'f_to_s')
	    IF(L .NE. 0)F_TO_S_FILE=F_TO_S_FILE(L:)
	  END IF
	END DO
!
	RETURN
	END
!
! Get F_TO_S file that exist in the latest data directory. The file
! chosen is the first one obtained when doing the data listing.
!
	SUBROUTINE GET_F_TO_S_FILE(PATH,F_TO_S_FILE,LU_DIAG)
	USE SET_KIND_MODULE
	IMPLICIT NONE
	CHARACTER(LEN=*) PATH,F_TO_S_FILE
	INTEGER LU_DIAG
	INTEGER K
	LOGICAL FILE_EXISTS
	INTEGER, SAVE :: LUIN
	LOGICAL, SAVE :: FIRST_TIME=.TRUE.
	CHARACTER(LEN=2000) STRING
!
	IF(FIRST_TIME)THEN
	  CALL GET_LU(LUIN,'In GET_F_TO_S_FILE (create_batch_ins.f)')
	  FIRST_TIME=.FALSE.
	END IF
!
	F_TO_S_FILE=' '
	INQUIRE(FILE='PATH_LISTING',EXIST=FILE_EXISTS)
	IF(FILE_EXISTS)THEN
	  WRITE(6,*)'Error -- file PATH_LISTING exists'
	  WRITE(6,*)'Rename or delete'
	  STOP
	END IF
	STRING='ls '//TRIM(PATH)//' > PATH_LISTING'
	CALL SYSTEM(STRING)
	OPEN(UNIT=LUIN,FILE='PATH_LISTING',ACTION='READ')
	DO WHILE(1 .EQ. 1)
	  READ(LUIN,*,END=100)STRING
	  IF(INDEX(STRING,'f_to_s') .NE. 0)THEN
	    F_TO_S_FILE=TRIM(STRING)
	    EXIT
	  END IF
	END DO
100	CONTINUE
	CLOSE(LUIN)
	CALL SYSTEM('rm -f PATH_LISTING')
!
	IF(F_TO_S_FILE .EQ. ' ')THEN
	  WRITE(6,*)'An f_to_s file was not found. Path =',TRIM(PATH)
	END IF
!
	RETURN
	END
