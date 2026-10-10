!
	SUBROUTINE GET_NUM_PHOT_ROUTES_V1(N_PHOT_TABLES,NUM_ROUTES,
	1             NTERMS,N_DP,MAX_FILES,
	1             XzSIX_PRES,EDGEXzSIX_F,GXzSIX_F,F_TO_S_XzSIX,
	1             XzSIX_LEVNAME_F,NXzSIX_F,ID,DESC,LUIN,LUOUT)
	USE SET_KIND_MODULE
	USE PHOT_DATA_MOD		!Contains all photoionization data
	IMPLICIT NONE
!
	INTEGER NUM_ROUTES
	INTEGER N_PHOT_TABLES
	INTEGER LUIN,LUOUT
	INTEGER MAX_FILES
	INTEGER NTERMS(MAX_FILES)
	INTEGER N_DP(MAX_FILES)
!
	INTEGER ID             			!Species identifier
	CHARACTER(LEN=*) DESC                   !Species identifier for PHOT file.
!
	INTEGER NXzSIX_F
	LOGICAL XzSIX_PRES
	REAL(KIND=LDP) EDGEXzSIX_F(NXzSIX_F)
	REAL(KIND=LDP) GXzSIX_F(NXzSIX_F)
	INTEGER F_TO_S_XzSIX(NXzSIX_F)
	CHARACTER*(*) XzSIX_LEVNAME_F(NXzSIX_F)
!
! External Functions.
!
	REAL(KIND=LDP) RD_FREE_VAL,SPEED_OF_LIGHT
	INTEGER ICHRLEN,ERROR_LU
	EXTERNAL RD_FREE_VAL,ERROR_LU,ICHRLEN,SPEED_OF_LIGHT
!
! Local variables.
!
	INTEGER, PARAMETER :: IZERO=0
!
	INTEGER N_ALT_STATES				!Number of alternate names in PHOT file.
	INTEGER, PARAMETER :: N_ALT_MAX=5		!Up to N_ALT_MAX alternate names
!
! Maximum number of ionization routes allowing for splitting of LS states.
!
	INTEGER, PARAMETER :: MAX_NUM_ROUTES=40

!
! Many of these will be transferred to the PHOTIOINZATION module one the actual number
! of photoionization routs are known.
!
	REAL(KIND=LDP) EXCITE_FREQ(MAX_NUM_ROUTES)              !Exictation energy of final level above ground state
	REAL(KIND=LDP) STAT_WGT(MAX_NUM_ROUTES)      		!Statstical weith of final state
	REAL(KIND=LDP) GSUM(MAX_FILES)
!
	INTEGER PHOT_GRID_ID(MAX_NUM_ROUTES)                    !Indicates which cross-section grid.
	INTEGER FIN_LEV_ID(MAX_NUM_ROUTES)
	CHARACTER(LEN=40) FINAL_STATE(N_ALT_MAX)
	CHARACTER(LEN=20) ION_NAME(MAX_NUM_ROUTES)
!
	INTEGER N_PHOT
!
! Work variables
!
	REAL(KIND=LDP) T1
	REAL(KIND=LDP) GION
	REAL(KIND=LDP) EION
!
	INTEGER PHOT_ID
	INTEGER ROUTE_CNT
	INTEGER LUER
	INTEGER I,J,L,K
	INTEGER NEXT
	INTEGER F_TO_S_SAVE
	INTEGER LN
	LOGICAL FOUND_FINAL_STATE
!
! Statistical weight of resulting ion for PHOT_ID=1. Used only if
! ion is not present.
!
	LUER=ERROR_LU()
	FINAL_STATE(:)=' '
!
! 
! Open the data files to determine the number of energy levels and the
! number of data pairs. This is done for all files belonging to the current
! species so that we can allocate the necessary memory.
!
	PHOT_ID=0
	N_PHOT=1
	ROUTE_CNT=0
!
! To get the total number of photoioization routes, we need to loop over
! all files, and determines which states are split.
!
	DO WHILE (PHOT_ID .LT. N_PHOT)
	  PHOT_ID=PHOT_ID+1
!
!
! This routie returns:
!      N_PHOT       -- The number of photoionization files to be read.
!      N_ALT_STATES -- The number of alternate level names for the final state
!      FINAL_STATE  -- The list of final states
!      EION         -- The excitation freqency (10^15 Hz) of the final state
!      GION         -- The statistical weight of the final state.
!
	  CALL GET_MAIN_PHOT_DATA

! Only one alternate should match.
!
	  GSUM(PHOT_ID)=0.0_LDP
          IF(XzSIX_PRES)THEN
            FOUND_FINAL_STATE=.FALSE.
	    F_TO_S_SAVE=0
            DO K=1,N_ALT_STATES
              DO L=1,NXzSIX_F
                LN=INDEX(XzSIX_LEVNAME_F(L),'[')-1
                IF(LN .LT. 0)LN=ICHRLEN(XzSIX_LEVNAME_F(L))
                IF(XzSIX_LEVNAME_F(L)(1:LN) .EQ. FINAL_STATE(K))THEN
                  IF(FOUND_FINAL_STATE)THEN
	            IF(F_TO_S_XzSIX(L) .NE. F_TO_S_SAVE)THEN
	              ROUTE_CNT=ROUTE_CNT+1
	              FIN_LEV_ID(ROUTE_CNT)=L
	              ION_NAME(ROUTE_CNT)= FINAL_STATE(K)
	              PHOT_GRID_ID(ROUTE_CNT)=PHOT_ID
	              STAT_WGT(ROUTE_CNT)=GXzSIX_F(L)
                      GSUM(PHOT_ID)=GSUM(PHOT_ID)+GXzSIX_F(L)
	              EXCITE_FREQ(ROUTE_CNT)=EDGEXzSIX_F(1)-EDGEXzSIX_F(L)
	            ELSE
	              STAT_WGT(ROUTE_CNT)=STAT_WGT(ROUTE_CNT)+GXzSIX_F(L)
                      GSUM(PHOT_ID)=GSUM(PHOT_ID)+GXzSIX_F(L)
	            END IF
	          ELSE
	            FOUND_FINAL_STATE=.TRUE.
	            ROUTE_CNT=ROUTE_CNT+1
	            FIN_LEV_ID(ROUTE_CNT)=L
	            F_TO_S_SAVE=F_TO_S_XzSIX(L)
	            ION_NAME(ROUTE_CNT)= FINAL_STATE(K)
	            PHOT_GRID_ID(ROUTE_CNT)=PHOT_ID
	            STAT_WGT(ROUTE_CNT)=GXzSIX_F(L)
                    GSUM(PHOT_ID)=GSUM(PHOT_ID)+GXzSIX_F(L)
	            EXCITE_FREQ(ROUTE_CNT)=EDGEXzSIX_F(1)-EDGEXzSIX_F(L)
	          END IF
                END IF
	      END DO
	      IF(FOUND_FINAL_STATE)EXIT
            END DO
	    IF(.NOT. FOUND_FINAL_STATE)THEN
	      WRITE(6,*)'Final state not found ','  ',TRIM(DESC); FLUSH(UNIT=6)
              WRITE(6,*)N_ALT_STATES
	      DO K=1,N_ALT_STATES
	        WRITE(6,*)FINAL_STATE(K)
	      END DO
	    END IF
	  ELSE
	    ROUTE_CNT=ROUTE_CNT+1
	    PHOT_GRID_ID(ROUTE_CNT)=PHOT_ID
	    FIN_LEV_ID(ROUTE_CNT)=1
	    EXCITE_FREQ(ROUTE_CNT)=EION
	    STAT_WGT(ROUTE_CNT)=GION
            GSUM(PHOT_ID)=GION
	    ION_NAME(ROUTE_CNT)=FINAL_STATE(1)
	  END IF
	  CLOSE(LUIN)
!
	END DO
!
! Allocate memory associated with the photoionization module.
!
	N_PHOT_TABLES=N_PHOT
	NUM_ROUTES=ROUTE_CNT
	L=NUM_ROUTES
	ALLOCATE( PD(ID)%EXC_FREQ(L) )
	ALLOCATE( PD(ID)%DO_PHOT(L,L) )
	ALLOCATE( PD(ID)%SCALE_FAC(L) )
	ALLOCATE( PD(ID)%PHOT_GRID(L) )
	ALLOCATE( PD(ID)%ION_LEV_ID(L) )
	ALLOCATE( PD(ID)%FINAL_ION_NAME(L) )
	ALLOCATE( PD(ID)%GION(L) )
!
! If we only have the ground state of the next ion, all photoionizations
! must go to that state.
!
	PD(ID)%DO_PHOT=.FALSE.
        IF(XzSIX_PRES)THEN
	  DO L=1,NUM_ROUTES
	    PD(ID)%DO_PHOT(L,L)=.TRUE.
	  END DO
	ELSE
	  DO L=1,NUM_ROUTES
	    PD(ID)%DO_PHOT(1,L)=.TRUE.
	  END DO
	END IF
!
! We can set the modulee variable since the total number of photoionization routes is known.
!
	DO L=1,NUM_ROUTES
	  PD(ID)%SCALE_FAC(L)=STAT_WGT(L)/GSUM(PHOT_GRID_ID(L))
	  PD(ID)%PHOT_GRID(L)=PHOT_GRID_ID(L)
	  PD(ID)%EXC_FREQ(L)=EXCITE_FREQ(L)
	  PD(ID)%ION_LEV_ID(L)=FIN_LEV_ID(L)
	  PD(ID)%FINAL_ION_NAME(L)=ION_NAME(L)
	  PD(ID)%GION(L)=STAT_WGT(L)
	END DO
!
	RETURN

	CONTAINS
	SUBROUTINE GET_MAIN_PHOT_DATA
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
	INTEGER L1,L2
	INTEGER IOS
	CHARACTER(LEN=1) APPEND_CHAR
	CHARACTER(LEN=200) FILENAME
	CHARACTER(LEN=200) STRING
!
	FILENAME='PHOT'//TRIM(DESC)//'_'
	I=ICHRLEN(FILENAME)
	APPEND_CHAR=CHAR( ICHAR('A')+(PHOT_ID-1) )
	WRITE(FILENAME(I+1:I+1),'(A)')APPEND_CHAR
	CALL GEN_ASCI_OPEN(LUIN,FILENAME,'OLD',' ','READ',IZERO,IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error in RDPHOT_GEN_V1: ',DESC
	  WRITE(LUER,*)'Unable to open file',FILENAME
	  STOP
	END IF
!
! Set rubbish values so we can check if read.
!
	NTERMS(PHOT_ID)=-1
	GION=-1
	EION=-1
	IF(PHOT_ID .EQ. 1)N_PHOT=-1
	N_ALT_STATES=-1
	N_DP(PHOT_ID)=-1
!
! Read in number of energy levels.
!
	STRING=' '
	DO WHILE(INDEX(STRING,'!Type of cross-section') .EQ. 0)
	  IF( INDEX(STRING,'!Number of energy levels') .NE. 0)THEN
	    READ(STRING,*,IOSTAT=IOS)NTERMS(PHOT_ID)
	    IF(IOS .NE. 0)THEN
	      WRITE(LUER,*)' Error reading number of NTERMS in RDPHOT_GEN_V1: ',DESC
	      WRITE(LUER,*)'PHOT_ID=',PHOT_ID
	      STOP
	    END IF
!
	  ELSE IF (INDEX(STRING,'!Number of photoionization routes') .NE.  0 .AND. PHOT_ID .EQ. 1)THEN
	    READ(STRING,*,IOSTAT=IOS)N_PHOT
	    IF(IOS .NE. 0)THEN
	      WRITE(LUER,*)'Error in RDPHOT_GEN_V1: ',DESC
	      WRITE(LUER,*)'Error reading number of photoionization routes'
	      STOP
	    END IF
!
	  ELSE IF(INDEX(STRING,'!Final state in ion') .NE. 0)THEN
!
! Read in the final state. We allow for the possibility of alternate names,
! separated by a '/'. No double spaces must be present in the names.
!
!
! We can extract the Final level. To do so we frist remove
! unwanted garbage, TABS, and excess blanks.
!
	    STRING=ADJUSTL(STRING)
	    DO WHILE(INDEX(STRING,CHAR(9)) .NE. 0)
	      K=INDEX(STRING,CHAR(9))	!Remove tabs
	      STRING(K:K)='  '
	    END DO
	    L1=INDEX(STRING,'  ')
	    STRING(L1:)=' '
!
	    L1=LEN_TRIM(STRING)
	    K=INDEX(STRING(1:L1),' ')	!Remove blanks
	    DO WHILE(K .NE. 0)
	      STRING(K:)=STRING(K+1:)
	      L1=L1-1
	      K=INDEX(STRING(1:L1),' ')
	    END DO
!
! Can now get final level name
!	
	    L2=INDEX(STRING,'  ')
	    L1=0
	    DO WHILE(L1 .LT. N_ALT_MAX)
	      L1=L1+1
	      K=INDEX(STRING(1:L2),'/')-1
	      IF(K .LE. 0)K=L2
	      FINAL_STATE(L1)=STRING(1:K)
	      IF(K .EQ. L2)EXIT
	      STRING(1:)=STRING(K+2:)
	      L2=INDEX(STRING,'  ')
	    END DO
	    N_ALT_STATES=L1
!
	  ELSE IF(INDEX(STRING,'!Excitation energy of final state') .NE. 0)THEN
!
! Now read in the excitation energy of the final state above the
! ground state. Excitation energy should be in cm^-1
!
	    L2=INDEX(STRING,'  ')
	    T1=RD_FREE_VAL(STRING,1,L2-1,NEXT,' EXC_FREQ RDPHOT_GEN_V1')
	    EION=1.0E-15_LDP*T1*SPEED_OF_LIGHT()
!
	  ELSE IF(INDEX(STRING,'!Statistical weight of ion') .NE.  0)THEN
!
! Get statistical weight of ion.
!
	    L2=INDEX(STRING,'  ')
	    GION=RD_FREE_VAL(STRING,1,L2-1,NEXT,' Statistical weight of ion')
!
	  ELSE IF( INDEX(STRING,'!Total number of data pairs') .NE. 0)THEN
!
! Read in number of data pairs.
!
	    READ(STRING,*,IOSTAT=IOS)N_DP(PHOT_ID)
	    IF(IOS .NE. 0)THEN
	      WRITE(LUER,*)' Error reading number of data  pairs in RDPHOT_GEN_V1: ',DESC
	      WRITE(LUER,*)'PHOT_ID=',PHOT_ID
	      STOP
	    END IF
	  END IF
	  READ(LUIN,'(A)')STRING
	END DO
!
	IF(NTERMS(PHOT_ID) .LT. 0)THEN
	    WRITE(LUER,*)' Number of terms (levels) not read in GET_MAIN_PHOT_DATA: ',DESC
	    WRITE(LUER,*)' Make sure ''!Type of cross-section'' is not in a comment line in the PHOT_FILE'
	    WRITE(LUER,*)' PHOT_ID=',PHOT_ID
	    STOP
	ELSE IF(PHOT_ID .LT. 0)THEN
	    WRITE(LUER,*)' Number of photoionizaton routes not read in GET_MAIN_PHOT_DATA: ',DESC
	    WRITE(LUER,*)' PHOT_ID=',PHOT_ID
	    STOP
	ELSE IF(GION .LT. 0)THEN
	    WRITE(LUER,*)' Statistical weight of final level not read in GET_MAIN_PHOT_DATA: ',DESC
	    WRITE(LUER,*)' PHOT_ID=',PHOT_ID
	    STOP
	ELSE IF(EION .LT. 0)THEN
	    WRITE(LUER,*)' Energy of final level not read in GET_MAIN_PHOT_DATA: ',DESC
	    WRITE(LUER,*)' PHOT_ID=',PHOT_ID
	    STOP
	ELSE IF(N_ALT_STATES.LT. 0)THEN
	    WRITE(LUER,*)' Final state not read in GET_MAIN_PHOT_DATA: ',DESC
	    WRITE(LUER,*)' PHOT_ID=',PHOT_ID
	    STOP
	ELSE IF(N_DP(PHOT_ID).LT. 0)THEN
	    WRITE(LUER,*)' Number of data pints not read in GET_MAIN_PHOT_DATA: ',DESC
	    WRITE(LUER,*)' PHOT_ID=',PHOT_ID
	    STOP
	END IF
!
	RETURN
!
	END SUBROUTINE GET_MAIN_PHOT_DATA
	END SUBROUTINE GET_NUM_PHOT_ROUTES_V1
