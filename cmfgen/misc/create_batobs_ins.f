!
! To avoid possible conflicts, the keys
!       [TRANS_XzV]
! and
!       [SCL_DUM_ABUND]
! should not exist in CMF_FLUX_PARAM_INIT.
!
! Spectral calsulations to be run are specified by the file BAT_PARAMS
! It should have the format
!
!  RUNID   (ID)
!	VALUE [KEYWORD]
!	VALUE [KEYWORD]
!
!  RUNID   (ID)
!    etc
!
	PROGRAM CREATE_BATOBS_INS
	USE SET_KIND_MODULE
	USE GEN_IN_INTERFACE
	USE MOD_COLOR_PEN_DEF
!
! Altered 17-May-2025 : Altered SED command to work on Mac.
! Altered 22-Jun-2023 : Fixed initialization error (affected some compilers).
! Altered 10-May-2023 : Improved documentation.
! Altered 26-Sep-2022 : Option to use lower case name for IDs, _ added automatically to identifiers.
! Created 19-Jun-2022
!
	IMPLICIT NONE
	INTEGER, PARAMETER :: NUM_KEYS_MAX=200
	INTEGER NUM_KEYS
	INTEGER N_ADDS
	INTEGER I,K,K1,K2,KEY_ID
	INTEGER IOS
!
	LOGICAL VALID_KEY
	LOGICAL DONE_SED
	LOGICAL ANS
	LOGICAL FILE_EXISTS
	LOGICAL STORE_EDDFACTOR
	LOGICAL KEEP_ALL
	LOGICAL LOWER_CASE
!
	CHARACTER(LEN=30), EXTERNAL :: LC
	CHARACTER(LEN=30) NAME_MOD
	CHARACTER(LEN=80) STRING
	CHARACTER(LEN=80) FILE_NAME
	CHARACTER(LEN=80) ADD_STORE(30)
	CHARACTER(LEN=20) CMF_KEYS(NUM_KEYS_MAX)
!
	INTEGER, PARAMETER :: LUIN=10
	INTEGER, PARAMETER :: LUOUT=30
!
	KEEP_ALL=.FALSE.
	LOWER_CASE=.TRUE.
!
	WRITE(6,*)BLUE_PEN
	WRITE(6,*)'BAT_PARAMS must exist'
	WRITE(6,*)'For simple model BAT_PARAMS should have the following three lines.'
	WRITE(6,*)'It must end with the blankline'
	WRITE(6,*)'   RUNID'
	WRITE(6,*)'   T      [TRAP_J]'
	WRITE(6,*)'  '
	WRITE(6,*)'It may contain any valid key which will be set to the value you enter.'
	WRITE(6,*)'After the blank line you can enter a new RUNID etc for a new spectrum calcualtion'
	WRITE(6,*)'    with a new set  of paramters'
	WRITE(6,*)DEF_PEN
!
	WRITE(6,'(A)')RED_PEN
	WRITE(6,'(A)')' After each model:'
	WRITE(6,'(A)')'      ES_J_CONV will be deleted'
	WRITE(6,'(A)')'      EDDFACTOR will be deleted unless STORE has been set in BAT_PARAMS'
	WRITE(6,'(A)')' '
	WRITE(6,'(A)')'  MEANOPC, J_COMP, OBSFLUX (obs_cmf) HYDRO and TIMING will be deleted'
	WRITE(6,'(A)')'      They will be saved if next option is set to true.'
	WRITE(6,'(A)')DEF_PEN
	CALL GEN_IN(KEEP_ALL,'Keep unnecessary file -- MEANOPC, J_COMP, etc?')
	CALL GEN_IN(LOWER_CASE,'Lower case name identifier?')
!
! Get keys in CMF_FLUX_PARAM_INIT. These will be compared to the
! requested keys. All requested keys must exist in CMF_FLUX_PARAM_INIT
! except those related to scaling of abundances, which may/may not
! exist in the file (better if they don't exist).
!
	K=0; NUM_KEYS=0; CMF_KEYS=' '
	OPEN(UNIT=LUIN,FILE='CMF_FLUX_PARAM_INIT',STATUS='OLD',ACTION='READ')
	DO WHILE(K .LT. NUM_KEYS_MAX)
	  READ(LUIN,'(A)',END=100)STRING
	  IF(STRING .EQ. ' ' .OR. STRING(1:1) .EQ. '!')THEN
	  ELSE
	    K1=INDEX(STRING,'[')
	    K2=INDEX(STRING,']')
	    IF(K1 .GT. 0 .AND. K2 .GT. K1)THEN
	       NUM_KEYS=NUM_KEYS+1
	       CMF_KEYS(NUM_KEYS)=STRING(K1:K2)
	       IF(INDEX(CMF_KEYS(NUM_KEYS),'[TRANS_') .NE. 0)THEN
	         WRITE(6,*)'To avoid posisble conflicts, TRANS_XzV keys'//
	1                  '  should be removed from CMF_FLUX_PARAM_INIT'
	         STOP
	       END IF
	    END IF
	  END IF
	END DO
100	CONTINUE
	CLOSE(UNIT=LUIN)
	WRITE(6,*)'Number of keys in CMF_FLUX_PARAM_INIT is',NUM_KEYS
!
! Open output file, checking if it exits, and prompting to overwrite.
!
	INQUIRE(FILE='bat_ins.sh',EXIST=FILE_EXISTS)
	IF(FILE_EXISTS)THEN
	  WRITE(6,*)' '
	  WRITE(6,'(A)',ADVANCE='NO')' Warning -- bat_ins.sh exits -- do you want to continue (T,F)? '
	  READ(5,*)ANS
	  IF(.NOT. ANS)STOP
	END IF
	OPEN(UNIT=LUOUT,FILE='bat_ins.sh',STATUS='UNKNOWN',ACTION='WRITE')
!
	WRITE(LUOUT,'(/,A)')
	1     '# Ensure CMF_FLUX_PARAM is not in the directory, since it will stop the editing.'
	WRITE(LUOUT,'(/,A,/)')'rm -f CMF_FLUX_PARAM'
!
! Read in required model runs.
!	
	OPEN(UNIT=LUIN,FILE='BAT_PARAMS',STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(6,*)'BAT_PARAMS must exist'
	  STOP
	END IF
!
	DO WHILE(1 .EQ. 1)
	  READ(LUIN,'(A)',END=1000)STRING
	  STORE_EDDFACTOR=.FALSE.
!
	  K=INDEX(STRING,'RUNID')
	  IF(K .NE. 0)THEN
	    N_ADDS=0
	    NAME_MOD='_'//ADJUSTL(STRING(K+5:))
	    IF(NAME_MOD .EQ. '_')NAME_MOD=' '
	    IF(LOWER_CASE)NAME_MOD=LC(NAME_MOD)
	    DONE_SED=.FALSE.
	    IF(INDEX(STRING,'LINK') .NE. 0)THEN
	      WRITE(LUOUT,'(/,A)')'ln -sf EDDFACTOR_STORE          EDDFACTOR'
	      WRITE(LUOUT,'(A,/)')  'ln -sf EDDFACTOR_STORE_INFO     EDDFACTOR_INFO'
	      READ(LUIN,'(A)',END=1000)STRING
	    END IF
	    IF(INDEX(STRING,'STORE') .NE. 0)THEN
	      STORE_EDDFACTOR=.TRUE.
	    END IF
	    DO WHILE(1 .EQ. 1)
	      READ(LUIN,'(A)',END=200)STRING
	      IF(STRING .EQ. ' ')EXIT
!
! Check if valid key
!
	      K1=INDEX(STRING,'[')
	      K2=INDEX(STRING,']')
	      VALID_KEY=.FALSE.
	      DO KEY_ID=1,NUM_KEYS
	         IF(STRING(K1:K2) .EQ. TRIM(CMF_KEYS(KEY_ID)))THEN
	           VALID_KEY=.TRUE.
	           EXIT
	         END IF
	      END DO
!
! Check for keys we are adding to file, and which may not exits in CMF_FLUX_PARAM_INIT.
!
	      IF(.NOT. VALID_KEY)THEN
	        IF(INDEX(STRING,'[SCL_') .NE. 0)THEN
	          N_ADDS=N_ADDS+1
	          ADD_STORE(N_ADDS)="sed    -i -e '$a"//TRIM(STRING)//"'  CMF_FLUX_PARAM"
	        ELSE IF(INDEX(STRING,'[TRANS_') .NE. 0)THEN
	          N_ADDS=N_ADDS+1
	          ADD_STORE(N_ADDS)="sed    -i -e '$a"//TRIM(STRING)//"'  CMF_FLUX_PARAM"
	        ELSE
	          WRITE(6,*)'ERROR -- KEY not recognized'
	          WRITE(6,*)'Key is: ',TRIM(STRING(K1:K2))
	          STOP
	        END IF
	      END IF
!
! Create string with key repacement. New line will replace line in
! CMF_FLUX_PARAM_INT. If VALID_KEY is false, KEY involves abundace
! scaling, and has already been treated. Replacementes are
! done at the same time.
!
	      IF(.NOT. VALID_KEY)THEN
	      ELSE IF(DONE_SED)THEN
	        STRING='    -e "/\'//STRING(K1:K2-1)//'\]/s/.*/'//TRIM(STRING)//'/" \'
	        WRITE(LUOUT,'(A)')TRIM(STRING)
	      ELSE
	        STRING='sed -e "/\'//STRING(K1:K2-1)//'\]/s/.*/'//TRIM(STRING)//'/" \'
	        DONE_SED=.TRUE.
	        WRITE(LUOUT,'(A)')TRIM(STRING)
	      END IF
	    END DO
!
200	    CONTINUE
	    WRITE(LUOUT,'(A)')'CMF_FLUX_PARAM_INIT >  CMF_FLUX_PARAM'
!
! Append abudance scaling strings (if present).
!
	    WRITE(LUOUT,'(A)')' '
	    DO I=1,N_ADDS
	     WRITE(LUOUT,'(A)')TRIM(ADD_STORE(I))
	    END DO
	    WRITE(LUOUT,'(A)')' '
!
! Clean up this model in prepartion for the next model.
!
	    WRITE(LUOUT,'(A)')'$PROG_CMF_OBS < IN_FILE  >>&  ''batobs.log'' &'
	    WRITE(LUOUT,'(A)')'echo " PID of "  $PROG_CMF_OBS "is:" $! >> ''batobs.log'' '
	    WRITE(LUOUT,'(A)')'wait'
	    WRITE(LUOUT,'(/,A,/)')'echo " Program finished on:" `date` >> ''batobs.log'' '

	    WRITE(LUOUT,'(/,A,/)')'# Clean up the directory structure, removing all uneccessary files.'

	    WRITE(LUOUT,'(A)')'mv -f OBSFRAME           obs_fin'//TRIM(NAME_MOD)
	    WRITE(LUOUT,'(A,/)')'mv -f CMF_FLUX_PARAM     CMF_FLUX_PARAM'//TRIM(NAME_MOD)
	    IF(KEEP_ALL)THEN
	      WRITE(LUOUT,'(A)')'mv -f OBSFLUX            obs_cmf'//TRIM(NAME_MOD)
	      WRITE(LUOUT,'(A)')'mv -f HYDRO              hydro_fin'//TRIM(NAME_MOD)
	      WRITE(LUOUT,'(A)')'mv -f MEANOPAC           meanopac'//TRIM(NAME_MOD)
	      WRITE(LUOUT,'(A)')'mv -f TIMING             full_timing'//TRIM(NAME_MOD)
	      WRITE(LUOUT,'(A,/)')'mv -f J_COMP             J_COMP'//TRIM(NAME_MOD)
	    ELSE
	      WRITE(LUOUT,'(A)')'rm -f OBSFLUX'
	      WRITE(LUOUT,'(A)')'rm -f HYDRO'
	      WRITE(LUOUT,'(A)')'rm -f MEANOPAC'
	      WRITE(LUOUT,'(A)')'rm -f TIMING'
	      WRITE(LUOUT,'(A)')'rm -f J_COMP'
	    END IF
	    WRITE(LUOUT,'(A)')'rm -f ES_J_CONV*'
	    WRITE(LUOUT,'(A)')'rm -f fort.*'
!
	    IF(STORE_EDDFACTOR)THEN
	      WRITE(LUOUT,'(/,A)')'mv EDDFACTOR        EDDFACTOR_STORE'
	      WRITE(LUOUT,'(A,///)')'mv EDDFACTOR_INFO   EDDFACTOR_STORE_INFO'
	    ELSE
	      WRITE(LUOUT,'(/,A)')'rm -f EDDFACTOR'
	      WRITE(LUOUT,'(A,///)')'rm -f EDDFACTOR_INFO'
	    END IF
!
	    FLUSH(UNIT=LUOUT)
	  END IF
	END DO
1000	CONTINUE
	CLOSE(LUIN)
	CLOSE(LUOUT)
!
	STOP
	END
