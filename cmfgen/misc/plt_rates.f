!
! Designed to examine rates populating and depopulating a particular  level.
!       TOTRATE
!       MODEL
!       MODEL_SPEC
!       RVTJ
!       XzVPRRR         (XzV is the species being investigated)
!	NON_TH_RATES    (in SN model)
!       XzV_IND_COL     (collision rates for species XzV)
!
! Program can also be used to look at the NETRATE file (in which in
!   case NON_TH_RATES and XzV_IND_COL) are not needed. Options
!   will need to be updated for NETRATE.
!

! Program will run quicker (since less to read in) if you use a secondary
! file create by:
!                grep XzV NETRATE >> SMALL_FILE
!
	PROGRAM PLT_RATES
	USE SET_KIND_MODULE
	USE GEN_IN_INTERFACE
	USE MOD_COLOR_PEN_DEF
	IMPLICIT NONE
!
! Altered: 20-Jun-2024: Some cleaning.
!                       Rates for EXD (for example) can now be output to file.
! Altered: 14-Apr-2023: Extensive modifications.
!                         Collision data now read in.
!                         LID option added.
! Altered: 20-Sep-2021: Will skip over collisional and continuum terms if present (OSIRIS: 17-Nov-2021).
! Altered:  6-Jul-2021: Small bug fixis (transferred from OSIRIS 21-Jul-20201)
! Created: 24-May-2020
!
	INTEGER, PARAMETER :: NMAX=600000
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
!
	REAL(KIND=LDP), ALLOCATABLE :: NU(:)
	REAL(KIND=LDP), ALLOCATABLE :: LAM(:)
!
	REAL(KIND=LDP), ALLOCATABLE :: XV(:)
	REAL(KIND=LDP), ALLOCATABLE :: YV(:)
	REAL(KIND=LDP), ALLOCATABLE :: ZV(:)
	REAL(KIND=LDP), ALLOCATABLE :: LONG_XV(:)
	REAL(KIND=LDP), ALLOCATABLE :: LONG_YV(:)
!
	CHARACTER(LEN=60), ALLOCATABLE :: TRANS_NAME(:)
	CHARACTER(LEN=60), ALLOCATABLE :: NT_TRANS_NAME(:)
	CHARACTER(LEN=60), ALLOCATABLE :: COL_TRANS_NAME(:)
	CHARACTER(LEN=60), ALLOCATABLE :: AUTO_LEV_NAME(:)
!
	REAL(KIND=LDP), ALLOCATABLE :: R(:)
	REAL(KIND=LDP), ALLOCATABLE :: T(:)
	REAL(KIND=LDP), ALLOCATABLE :: V(:)
	REAL(KIND=LDP), ALLOCATABLE :: SIGMA(:)
	REAL(KIND=LDP), ALLOCATABLE :: ED(:)
!
	CHARACTER(LEN=10) PR_SPECIES
	CHARACTER(LEN=10) SPECIES
	CHARACTER(LEN=30) LEVEL
!
	CHARACTER(LEN=40), ALLOCATABLE :: LEVEL_NAME(:)
	REAL(KIND=LDP), ALLOCATABLE :: STAT_WT(:)
	REAL(KIND=LDP), ALLOCATABLE :: ENERGY(:)
	REAL(KIND=LDP), ALLOCATABLE :: FEDGE(:)
	INTEGER, ALLOCATABLE :: F_TO_S(:)
	INTEGER, ALLOCATABLE :: INT_SEQ(:)
	REAL(KIND=LDP) IONIZATION_ENERGY,ZION
	CHARACTER(LEN=20) OSCDATE
	CHARACTER(LEN=20) SL_OPTION
	CHARACTER(LEN=20) FL_OPTION
	CHARACTER(LEN=20) IL_OPTION
!
	INTEGER, ALLOCATABLE ::  LINK(:)
	REAL(KIND=LDP), ALLOCATABLE :: RATES(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: NEW_RATES(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: SUM_RATES(:)
!
	INTEGER, ALLOCATABLE :: NT_LINK(:)
	REAL(KIND=LDP),  ALLOCATABLE :: NT_RATES(:,:)
	REAL(KIND=LDP),  ALLOCATABLE :: NEW_NT_RATES(:,:)
!
	INTEGER, ALLOCATABLE :: COL_LINK(:)
	REAL(KIND=LDP),  ALLOCATABLE :: COL_RATES(:,:)
	REAL(KIND=LDP),  ALLOCATABLE :: NEW_COL_RATES(:,:)
!
	REAL(KIND=LDP), ALLOCATABLE :: PHOT_RATE(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: REC_RATE(:,:)
!
	REAL(KIND=LDP), ALLOCATABLE :: AUTO_RATE(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: AUTO_REC_RATE(:,:)
!
	REAL(KIND=LDP),  ALLOCATABLE :: DPTH_VEC(:)
	REAL(KIND=LDP),  ALLOCATABLE :: NORM_VEC(:)
	REAL(KIND=LDP),  ALLOCATABLE :: WRK_VEC(:)
	INTEGER, ALLOCATABLE :: INDX(:)
!
	REAL(KIND=LDP) T1,T2
	REAL(KIND=LDP) RATE_THRESHOLD
	REAL(KIND=LDP) SUM_TO,SUM_FROM
	REAL(KIND=LDP) NT_SUM_TO,NT_SUM_FROM
	REAL(KIND=LDP) COL_SUM_TO,COL_SUM_FROM
	REAL(KIND=LDP) NET_AUTO
	REAL(KIND=LDP) MAX_RATE
	REAL(KIND=LDP) CUR_SUM
	INTEGER ND
	INTEGER CURVE_TITLE_LENGTH
	INTEGER N_LINES
	INTEGER N_TRANS
	INTEGER N_NT_LINES
	INTEGER N_NT_TRANS
	INTEGER N_COL_LINES
	INTEGER N_COL_TRANS
	INTEGER N_AUTO
	INTEGER RD_COUNT
	INTEGER I,J,K,L,ML
	INTEGER LU_TMP
	INTEGER N_TRANS_LIM
	INTEGER SL_INDX
	INTEGER CNT
	INTEGER MAX_TRANS_LENGTH
	INTEGER DPTH_INDX
	INTEGER NF            	!Number of full levels
	INTEGER NLEV		!Number of super levels
	INTEGER LU
	LOGICAL FILE_OPEN
	LOGICAL TERM_OUTPUT
	LOGICAL DO_COL_RATES
	LOGICAL DO_NT_RATES
	LOGICAL DO_AUTO_RATES
	LOGICAL DOING_NETRATE
	INTEGER IOS
!
	INTEGER, PARAMETER :: LUIN=11
	INTEGER, PARAMETER :: LUOUT=12	
	INTEGER IVEC(2)
!
	CHARACTER(LEN=80) FILENAME
	CHARACTER(LEN=80) TMP_STR
	CHARACTER(LEN=80) TMP_FMT
	CHARACTER(LEN=200) STRING
        CHARACTER(LEN=30) UC
	CHARACTER(LEN=20) PLT_OPT
	CHARACTER(LEN=20) XLABEL
	CHARACTER(LEN=20) XLABEL_SAVE
	CHARACTER(LEN=20) YLABEL
	CHARACTER(LEN=60) CURVE_LAB
        EXTERNAL UC
!
	NLEV=0
	SPECIES=' ';  LEVEL=' '; PR_SPECIES=' '
	SL_OPTION=' '
	FL_OPTION=' '
	IL_OPTION='NOCHANGE'
	RATE_THRESHOLD=0.05_LDP
	CURVE_TITLE_LENGTH=120
!
	WRITE(6,'(A)')' '
	WRITE(6,'(A)')' Program to plot/examine NETRATE or TOTRATE'
	WRITE(6,'(A)')' '
!
	DOING_NETRATE=.FALSE.
100	CONTINUE
 	FILENAME='TOTRATE'
	CALL GEN_IN(FILENAME,'File with data to be plotted/examined')
	IF(FILENAME .EQ. ' ')GOTO 100
	IF(INDEX(FILENAME,'NETRATE') .NE. 0)DOING_NETRATE=.TRUE.
!
! Get number of depth points
!
	OPEN(UNIT=LUIN,FILE='MODEL',STATUS='OLD',IOSTAT=IOS)
	  IF(IOS .EQ. 0)THEN
	    DO WHILE(1 .EQ. 1)
	      READ(LUIN,'(A)',IOSTAT=IOS)STRING
	      IF(IOS .NE. 0)EXIT
	      IF(INDEX(STRING,'!Number of depth points') .NE. 0)THEN
	         READ(STRING,*)ND
	         WRITE(6,'(A,I4)')' Number of depth points in the model is:',ND
	         EXIT
	      END IF
	    END DO
	  END IF
	INQUIRE(UNIT=LUIN,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(UNIT=LUIN)
	IF(IOS .NE. 0)THEN
	  WRITE(6,*)' Unable to open MODEL file to get # of depth points'
	  CALL GEN_IN(ND,'Number of data points (must be exact)')
	END IF
!
	N_LINES=NMAX
	CALL GEN_IN(N_LINES,'Maximum number of line transitions to be read')
!
! These are updated when the data file is read in.
!
	ALLOCATE (NORM_VEC(ND))
	ALLOCATE (DPTH_VEC(ND))
	ALLOCATE (RATES(N_LINES,ND))
	ALLOCATE (TRANS_NAME(N_LINES))
	ALLOCATE (NU(N_LINES))
	ALLOCATE (LAM(N_LINES))
!
! These are used when studing a single level. LINK is used
! to access the data in th eoriginal arrays.
!
	ALLOCATE (NEW_RATES(N_LINES,ND))
	ALLOCATE (LINK(N_LINES))
	ALLOCATE (SUM_RATES(N_LINES))
!
! Used for ploting.
!
	ALLOCATE (XV(ND),YV(ND),ZV(ND))
!
! These are updated when the data file is read in.
!
	IF(DOING_NETRATE)THEN
	  DO_NT_RATES=.FALSE.
	ELSE
	  INQUIRE(FILE='NON_TH_RATES',EXIST=DO_NT_RATES)
	END IF
	IF(DO_NT_RATES)THEN
	  ALLOCATE (NT_RATES(N_LINES,ND))
	  ALLOCATE (NT_TRANS_NAME(N_LINES))
!
! These are used when studing a single level. LINK is used
! to access the data in th eoriginal arrays.
!
	  ALLOCATE (NEW_NT_RATES(N_LINES,ND))
	  ALLOCATE (NT_LINK(N_LINES))
	ELSE
	  WRITE(6,*)'File (NON_TH_RATES) containing non-thermal rates not found'
	  WRITE(6,*)'Non-thermal data will NOT be read in'
	END IF
!
! Allocate storage fore lines
!
	ALLOCATE (COL_RATES(N_LINES,ND))
	ALLOCATE (COL_TRANS_NAME(N_LINES))
!
! These are used when studing a single level. LINK is used
! to access the data in th eoriginal arrays.
!
	ALLOCATE (NEW_COL_RATES(N_LINES,ND))
	ALLOCATE (COL_LINK(N_LINES))
!
!
	OPEN(UNIT=LUIN,FILE=FILENAME,STATUS='OLD',ACTION='READ')
!
	WRITE(6,'(A)')RED_PEN
	WRITE(6,'(A)')' Reading bound-bound data -- this may take a while'
	WRITE(6,'(A)')DEF_PEN
	RD_COUNT=0
!
! Skip over collisional and continuum terms if they are present.
!
	IF(FILENAME .EQ. 'TOTRATE')THEN
	  DO WHILE(1 .EQ. 1)
	    READ(LUIN,'(A)')STRING
	    IF(INDEX(STRING,'(') .NE. 0)EXIT
	  END DO
	  BACKSPACE(LUIN)
	END IF
!
	CALL TUNE(1,'READ')
	DO ML=1,N_LINES
	  STRING=' '
	  DO WHILE(STRING .EQ. ' ' .OR. STRING(1:2) .EQ. '--')
	    READ(LUIN,'(A)',END=5000)STRING
	    IF(INDEX(STRING,'STEQ') .NE. 0)GOTO 5000
	  END DO
	  STRING=ADJUSTL(STRING)
	  K=INDEX(STRING,'  ')
	  STRING=STRING(K+2:)
	  K=INDEX(STRING,'  ')
	  IF(K .GT. LEN(TRANS_NAME(1)))THEN
	   WRITE(6,*)'Error -- TRANS_NAME too short: needed length=',K
	   STOP
	  ELSE
	    TRANS_NAME(ML)=STRING(1:K)
	  END IF
	  READ(STRING(K+1:),*)NU(ML)
	  RD_COUNT=ML
	  READ(LUIN,*)(RATES(ML,I),I=1,ND)
	END DO
!
5000	CONTINUE
	CALL TUNE(2,'READ')
	CALL TUNE(3,'TERMINAL')

	N_LINES=RD_COUNT
	WRITE(6,'(A)')BLUE_PEN
	WRITE(6,*)'Number of lines read is',N_LINES
	WRITE(6,'(A)')DEF_PEN
!
        LAM(1:N_LINES)=2.99792458E+03_LDP/NU(1:N_LINES)
	LU=10
	FILENAME='RVTJ'
	ALLOCATE (R(ND),V(ND),ED(ND),T(ND),SIGMA(ND))
	CALL RD_SING_VEC_RVTJ(R,ND,'Radius',FILENAME,LU,IOS)
	CALL RD_SING_VEC_RVTJ(V,ND,'Velocity',FILENAME,LU,IOS)
	CALL RD_SING_VEC_RVTJ(T,ND,'Temperature',FILENAME,LU,IOS)
	CALL RD_SING_VEC_RVTJ(SIGMA,ND,'dlnV/dlnr-1',FILENAME,LU,IOS)
	CALL RD_SING_VEC_RVTJ(ED,ND,'Electron',FILENAME,LU,IOS)
!
	IF(DO_NT_RATES)THEN
	  WRITE(6,'(A)')RED_PEN
	  WRITE(6,'(/,A,/)')' Reading in non-thermal rates data -- this may take a while'
	  WRITE(6,'(A)')DEF_PEN
	  OPEN(UNIT=LUIN,FILE='NON_TH_RATES',STATUS='OLD',ACTION='READ')
!
	  RD_COUNT=0
	  CALL TUNE(1,'READ_NT')
	  DO ML=1,N_LINES             !NT_LINES
	    STRING=' '
	    DO WHILE(STRING .EQ. ' ' .OR. STRING(1:2) .EQ. '--')
	      READ(LUIN,'(A)',END=6000)STRING
	    END DO
	    STRING=ADJUSTL(STRING)
	    K=INDEX(STRING,'  ')
	    NT_TRANS_NAME(ML)=STRING(1:K)
	    RD_COUNT=ML
	    READ(LUIN,*)(NT_RATES(ML,I),I=1,ND)
	  END DO
6000	  CONTINUE
	  N_NT_LINES=RD_COUNT
	  CALL TUNE(2,'READ_NT')
	  WRITE(6,'(A)')BLUE_PEN
	  WRITE(6,*)'Number of number of nonthermal lines read is',N_NT_LINES
	  WRITE(6,'(A)')DEF_PEN
	ELSE
	  N_NT_LINES=0
	END IF
	CALL TUNE(3,'TERMINAL')
!
! Set def axis.
!
	XV(1:ND)=V(1:ND)
	XLABEL='V(km/s)'
!
! Enter main plotting/examination loop.
!
	DO WHILE(1 .EQ. 1)
	  WRITE(6,*)' '
	  PLT_OPT='P'
	  CALL GEN_IN(PLT_OPT,'Plot option')
!
	  IF(UC(PLT_OPT) .EQ. 'XLOGR')THEN
	    T1=R(ND)
	    DO I=1,ND
	      XV(I)=LOG10(R(I)/T1)
	    END DO
	    XLABEL='Log R/R(ND)'
!
	  ELSE IF(UC(PLT_OPT) .EQ. 'XV')THEN
	    XV(1:ND)=V(1:ND)
	    XLABEL='V(km/s)'
!
	  ELSE IF(UC(PLT_OPT) .EQ. 'XLOGV')THEN
	    DO I=1,ND
	      IF(V(I) .LE. 0.0_LDP)THEN
	        WRITE(6,*)'Error in PLT_RATES: XLOGV requires positive V'
	        WRITE(6,*)'Depth, V=',I,V(I)
	        STOP
	      END IF
	    END DO
	    XV(1:ND)=LOG10(V(1:ND))
	    XLABEL='Log V(km/s)'
!
	  ELSE IF(UC(PLT_OPT) .EQ. 'XN')THEN
	    DO I=1,ND
	      XV(I)=I
	    END DO
	    XLABEL='Depth index'
!
! Set the species, and reads in the SPECIES//'PRRR' file.
!
	  ELSE IF(UC(PLT_OPT(1:4)) .EQ. 'SPEC')THEN
	    CALL GEN_IN(SPECIES,'Species that will be examined: e.g. He2, C2, NIII (case sensitive')
!
! Get number of super levels. This is need when we read in the the
! recombination and photoioizations rates for each super level.
!
	    WRITE(6,'(A)')' Opening MODEL_SPEC'
	    FILENAME='MODEL_SPEC'
	    OPEN(UNIT=LUIN,FILE=FILENAME,STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	      IF(IOS .NE. 0)THEN
	        WRITE(6,*)'Unable to open ',TRIM(FILENAME)
	        GOTO 1000
	      END IF
	      STRING=' '; SL_OPTION=' '
	      DO WHILE(INDEX(STRING,'[SL_OPT]') .EQ. 0)
	        READ(LUIN,'(A)',IOSTAT=IOS)STRING
	        IF(INDEX(STRING,'_ISF]') .NE. 0)EXIT		!string not present.
	      END DO
	      K=INDEX(STRING,'[SL_OPT]')
	      IF(K .NE. 0)SL_OPTION=STRING(1:K-1)
	      DO WHILE(INDEX(STRING,TRIM(SPECIES)) .EQ. 0)
	        READ(LUIN,'(A)',IOSTAT=IOS)STRING
	      END DO
	      WRITE(6,*)'SL_OPTION is ',SL_OPTION
	      IF(IOS .EQ. 0)THEN
	        READ(STRING,*)I,NLEV,NF
	        WRITE(6,*)NF,NLEV
	      ELSE
	        WRITE(6,*)'Error - species not found: IOS=',IOS
	        CLOSE(UNIT=LUIN)
	        GOTO 1000
	      END IF
	    CLOSE(UNIT=LUIN)
	    WRITE(6,'(A)')' Read MODEL_SPEC'
	    IF(ALLOCATED(PHOT_RATE))DEALLOCATE(PHOT_RATE,REC_RATE)
	    ALLOCATE (PHOT_RATE(NLEV,ND),REC_RATE(NLEV,ND))
!
! Read in the recombination and photoioizations rates for each super level.
! Note that ND is the second index.

	    WRITE(6,'(A)')' Opening PRRR file'
	    FILENAME=TRIM(SPECIES)//'PRRR'
	    OPEN(UNIT=LUIN,FILE=FILENAME,STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	    IF(IOS .NE. 0)THEN
	      WRITE(6,*)'Unable to open ',TRIM(FILENAME)
	      WRITE(6,*)'Was you species name correct: Species=',TRIM(SPECIES)
	      GOTO 1000
	    END IF
!
	    DO K=1,ND,10
	      DO WHILE(INDEX(STRING,'Photoionization') .EQ. 0)
	        READ(LUIN,'(A)')STRING
	      END DO
	      DO I=1,NLEV
	        READ(LUIN,*)(PHOT_RATE(I,L),L=K,MIN(ND,K+9))
	      END DO
	      DO WHILE(INDEX(STRING,'Recombination') .EQ. 0)
	        READ(LUIN,'(A)')STRING
	      END DO
	      DO I=1,NLEV
	        READ(LUIN,*)(REC_RATE(I,L),L=K,MIN(ND,K+9))
	      END DO
	    END DO
	    CLOSE(UNIT=LUIN)
	    PR_SPECIES=SPECIES
	    WRITE(6,'(A)')' Read PRRR file'
!
	    IF(DOING_NETRATE)THEN
	      DO_AUTO_RATES=.FALSE.
	    ELSE
	      DO_AUTO_RATES=.FALSE.
	      CALL GEN_IN(DO_AUTO_RATES,'Read in autoionization rates?')
	    END IF
	    IF(DO_AUTO_RATES)THEN
	      FILENAME=TRIM(SPECIES)//'_AUTO_RATES'
	      CALL GEN_IN(FILENAME,'File with auto rates')
	      OPEN(UNIT=LUIN,FILE=FILENAME,STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	      STRING=' '
	      DO WHILE(INDEX(STRING,'Number of ') .EQ. 0)
	        READ(LUIN,'(A)')STRING
	      END DO
	      READ(STRING,*)N_AUTO
	      READ(LUIN,'(A)')STRING; STRING=' '
	      DO WHILE(INDEX(STRING,'!') .NE. 0 .OR. STRING .EQ. ' ')
	        READ(LUIN,'(A)')STRING
	      END DO
	      BACKSPACE(UNIT=LUIN)
	      IF(ALLOCATED(AUTO_RATE))DEALLOCATE(AUTO_RATE,AUTO_REC_RATE,AUTO_LEV_NAME)
	      ALLOCATE (AUTO_RATE(N_AUTO,ND))
	      ALLOCATE (AUTO_REC_RATE(N_AUTO,ND))
	      ALLOCATE (AUTO_LEV_NAME(N_AUTO))
	      DO I=1,N_AUTO
	        STRING=' '
	        DO WHILE(INDEX(STRING,'!') .NE. 0 .OR. STRING .EQ. ' ')
	          READ(LUIN,'(A)')STRING
	        END DO
	        STRING=ADJUSTL(STRING)
	        K=INDEX(STRING,'  ');  AUTO_LEV_NAME(I)=STRING(1:K-1)
	        READ(LUIN,*)(AUTO_REC_RATE(I,K),AUTO_RATE(I,K),K=1,ND)
	      END DO
	      CLOSE(UNIT=LUIN)
	    END IF
	    LEVEL=' '
	    N_TRANS=0
!
	    IF(DOING_NETRATE)THEN
	      DO_COL_RATES=.FALSE.
	    ELSE
	      FILENAME=TRIM(SPECIES)//'_IND_RATES'
	      INQUIRE(FILE=FILENAME,EXIST=DO_COL_RATES)
	    END IF
	    IF(DO_COL_RATES)THEN
	      OPEN(UNIT=LUIN,FILE=TRIM(SPECIES)//'_IND_RATES',STATUS='OLD',ACTION='READ')
	      WRITE(6,'(A)')RED_PEN
	      WRITE(6,'(A)')' Reading collisonal data -- this may take a while'
	      WRITE(6,'(A)')DEF_PEN
	      RD_COUNT=0
	      CALL TUNE(1,'READ_COL')
	      RD_COUNT=0
	      DO ML=1,N_LINES
	        STRING=' '
	        DO WHILE(STRING .EQ. ' ' .OR. STRING(1:2) .EQ. '--')
	          READ(LUIN,'(A)',END=6200)STRING
	        END DO
	     	RD_COUNT=RD_COUNT+1
	        STRING=ADJUSTL(STRING)
	        K=INDEX(STRING,'  ')
	        COL_TRANS_NAME(RD_COUNT)=STRING(1:K)
	        READ(LUIN,*)(COL_RATES(RD_COUNT,I),I=1,ND)
!
! Because of the way we write the transitions, some may be written
! twice.
!
	        DO J=1,RD_COUNT-1
	          IF(COL_TRANS_NAME(RD_COUNT) .EQ. COL_TRANS_NAME(J))THEN
	            RD_COUNT=RD_COUNT-1
	            EXIT
	          END IF
	        END DO
	      END DO
6200	      CONTINUE
	      N_COL_LINES=RD_COUNT
	      CALL TUNE(2,'READ_COL')
	      WRITE(6,'(A)')BLUE_PEN
	      WRITE(6,*)'Number of number of collsional transitions read is',N_COL_LINES
	      WRITE(6,'(A)')DEF_PEN
	      DO I=1,N_COL_LINES
	        WRITE(35,'(A)')TRIM(COL_TRANS_NAME(I))
	      END DO
	      FLUSH(UNIT=35)
	    ELSE
	      WRITE(6,'(A)')BLUE_PEN
	      WRITE(6,*)'Collisonal data is unavailable'
	      WRITE(6,'(A)')DEF_PEN
	      N_COL_LINES=0
	    END IF
!
	    IF(ALLOCATED(LEVEL_NAME))DEALLOCATE(LEVEL_NAME,STAT_WT,ENERGY,FEDGE)
	    ALLOCATE (LEVEL_NAME(NF),STAT_WT(NF),ENERGY(NF),FEDGE(NF))
	    FILENAME=TRIM(SPECIES)//'_F_OSCDAT'
	    WRITE(6,*)'Getting energy levels',TRIM(SPECIES),'/',NF
	    CALL RD_ENERGY_V2(LEVEL_NAME,STAT_WT,ENERGY,FEDGE,NF,NF,L_TRUE,
	1             IONIZATION_ENERGY,ZION,OSCDATE,FILENAME,LUIN,LUOUT,IOS)
	    IF(IOS .NE. 0)THEN
	      WRITE(6,*)'Unable to read in energy level names'
	      GOTO 1000
	    ELSE
	      WRITE(6,*)'Got energy levels'
	    END IF
	    CLOSE(LUIN)
!
	    WRITE(6,*)'Setting SL assignent'
	    IF(ALLOCATED(F_TO_S))DEALLOCATE(F_TO_S,INT_SEQ)
	    ALLOCATE (F_TO_S(NF),INT_SEQ(NF))
	    FILENAME=TRIM(SPECIES)//'_F_TO_S'
	    CALL FDG_F_TO_S_NS_V1(NF,NLEV,K,FL_OPTION,SL_OPTION,IL_OPTION,LUIN,FILENAME)
	    CALL RD_F_TO_S_IDS_V2(F_TO_S,INT_SEQ,LEVEL_NAME,NF,NLEV,LUIN,TRIM(FILENAME),SL_OPTION)
	    IF(SL_OPTION .EQ. ' ')THEN
	      WRITE(6,*)'Read in supper level asignments'
	    ELSE
	      WRITE(6,*)'Read in supper level asignments with SL_OPTION ',TRIM(SL_OPTION)
	    END IF
!
	  ELSE IF(UC(PLT_OPT(1:3)) .EQ. 'LID')THEN
	    DO I=1,NF,20
	      WRITE(6,'(A)')' '
	      DO J=I,MIN(NF,I+9)
	        IF(J+10 .LE. NF)THEN
	          WRITE(6,'(1X,I4,3X,A,T50,I4,3X,A)')J,TRIM(LEVEL_NAME(J)),J+10,TRIM(LEVEL_NAME(J+10))
	        ELSE
	          WRITE(6,'(1X,I4,3X,A)')J,TRIM(LEVEL_NAME(J))
	        END IF
	      END DO
	      WRITE(6,'(A)')' '
	      TMP_STR='EXIT'
	      CALL GEN_IN(TMP_STR,'Any character to continue listing')
	      IF(TMP_STR .EQ. 'EXIT')exit
	    END DO
!
! Set the level for which rates will be examined, and store the
! data in a new array.
!
	  ELSE IF(UC(PLT_OPT) .EQ. 'LEVEL')THEN
	    N_TRANS=0
	    IF(SPECIES .EQ. ' ')THEN
	       WRITE(6,*)'Need to use SPECIES option first to set the species'
	       GOTO 1000
	    END IF
	    CALL GEN_IN(LEVEL,'Level -- full name -- case sensitive')
	    MAX_TRANS_LENGTH=0
!
! Find transition involvibg specified level. We as '(' to SPECIES to
! avoid degeneracies between OI and OIII.
!
	    DO ML=1,N_LINES
	       IF(INDEX(TRANS_NAME(ML),TRIM(SPECIES)//'(') .NE. 0 .AND.
	1         ( INDEX(TRANS_NAME(ML),'('//TRIM(LEVEL)) .NE. 0 .OR.
	1           INDEX(TRANS_NAME(ML),'-'//TRIM(LEVEL)) .NE. 0) )THEN
	          N_TRANS=N_TRANS+1
	          NEW_RATES(N_TRANS,:)=RATES(ML,:)
	          LINK(N_TRANS)=ML
	          MAX_TRANS_LENGTH=MAX(MAX_TRANS_LENGTH,LEN_TRIM(TRANS_NAME(ML)))
	       END IF
	    END DO
	    IF(N_TRANS .EQ. 0)THEN
	      WRITE(6,'(/,2A)')RED_PEN,' Error no transitions identified'
	      WRITE(6,'(2A,/)')' Check species and level names are consistent',DEF_PEN
	      GOTO 1000
	    ELSE
	      WRITE(6,*)'Number of bound-bound transitions identified is',N_TRANS
	    END IF
!
	    N_NT_TRANS=0
	    IF(DO_NT_RATES)THEN
	      DO ML=1,N_NT_LINES
	        IF(INDEX(NT_TRANS_NAME(ML),TRIM(SPECIES)//'(') .NE. 0 .AND.
	1          ( INDEX(NT_TRANS_NAME(ML),'('//TRIM(LEVEL)) .NE. 0 .OR.
	1            INDEX(NT_TRANS_NAME(ML),'-'//TRIM(LEVEL)) .NE. 0) )THEN
	          N_NT_TRANS=N_NT_TRANS+1
	          NEW_NT_RATES(N_NT_TRANS,:)=NT_RATES(ML,:)
	          NT_LINK(N_NT_TRANS)=ML
	          MAX_TRANS_LENGTH=MAX(MAX_TRANS_LENGTH,LEN_TRIM(NT_TRANS_NAME(ML)))
	        END IF
	      END DO
	      WRITE(6,*)'Number of non_thermal transitions identified is',N_NT_TRANS
	    END IF
!
	    N_COL_TRANS=0
	    IF(DO_COL_RATES)THEN
	      DO ML=1,N_COL_LINES
	        IF(INDEX(COL_TRANS_NAME(ML),TRIM(SPECIES)//'(') .NE. 0 .AND.
	1          ( INDEX(COL_TRANS_NAME(ML),'('//TRIM(LEVEL)) .NE. 0 .OR.
	1            INDEX(COL_TRANS_NAME(ML),'-'//TRIM(LEVEL)) .NE. 0) )THEN
	          N_COL_TRANS=N_COL_TRANS+1
	          NEW_COL_RATES(N_COL_TRANS,:)=COL_RATES(ML,:)
	          COL_LINK(N_COL_TRANS)=ML
	          MAX_TRANS_LENGTH=MAX(MAX_TRANS_LENGTH,LEN_TRIM(COL_TRANS_NAME(ML)))
	        END IF
	      END DO
	      WRITE(6,*)'Number of collisional transitions identified is',N_COL_TRANS
	    END IF
!
! Get the SL link.
!
	  DO I=1,NF
	    IF(INDEX(LEVEL_NAME(I),TRIM(LEVEL)) .NE. 0)THEN
	      SL_INDX=F_TO_S(I)
	      EXIT
	    END IF
	  END DO
!
	  WRITE(6,*)' '
	  WRITE(6,*)' The following levels have the same super level assignment'
	  DO I=1,NF
	    IF(F_TO_S(I) .EQ. SL_INDX)THEN
	      WRITE(6,'(4X,A,T50,A)')TRIM(LEVEL_NAME(I)),TRIM(LEVEL)
	    END IF
	  END DO
!
! Plot the normalized rates at each depth
!
	  ELSE IF(UC(PLT_OPT) .EQ. 'NORM')THEN
	    IF(N_TRANS .EQ. 0)THEN
	      WRITE(6,*)'Error -- no transitions to examine'
	      WRITE(6,*)'Use SPEC and LEVELS options'
	      GOTO 1000
	    END IF
!
	    DPTH_INDX=ND/2
	    CALL GEN_IN(DPTH_INDX,'Depth for choosing important transitions')
	    N_TRANS_LIM=10
	    CALL GEN_IN(N_TRANS_LIM,'Transition limit')
	    CALL GEN_IN(RATE_THRESHOLD,'Omit rates below this fractional threshold')
	    CALL GEN_IN(CURVE_TITLE_LENGTH,'As labels can be long, set long title')
	    CALL PG_SET_CURVE_TIT_LENGTH(CURVE_TITLE_LENGTH)
!
	    WRITE(6,*)' '
	    WRITE(6,'(2X,A,I5)')         '      Depth index =',DPTH_INDX
	    WRITE(6,'(2X,A,ES16.8,2X,A)')'     R(DPTH_INDX) =',R(DPTH_INDX),'(10^{10} cm)'
	    WRITE(6,'(2X,A,ES16.8,2X,A)')'     V(DPTH_INDX) =',V(DPTH_INDX),'(km/s)'
	    WRITE(6,'(2X,A,ES16.8,2X,A)')'     T(DPTH_INDX) =',T(DPTH_INDX),'(10^4 K)'
	    WRITE(6,'(2X,A,ES16.8,2X,A)')'    ED(DPTH_INDX) =',ED(DPTH_INDX),'(cm^{-3})'
	    WRITE(6,'(2X,A,ES16.8,2X,A)')'  Scaling Factor  =',MAX_RATE,'(cm^{-3} s^{-1})'
	    WRITE(6,*)' '
!
	    DO I=1,ND
	      T1=0.0_LDP
	      DO ML=1,N_TRANS
	         T1=MAX(T1,ABS(NEW_RATES(ML,I)))
	      END DO
	      DO ML=1,N_NT_TRANS
	         T1=MAX(T1,ABS(NEW_NT_RATES(ML,I)))
	      END DO
	      DO ML=1,N_COL_TRANS
	         T1=MAX(T1,ABS(NEW_COL_RATES(ML,I)))
	      END DO
	      T1=MAX(T1,REC_RATE(SL_INDX,I))
	      T1=MAX(T1,PHOT_RATE(SL_INDX,I))
	      NORM_VEC(I)=T1
	    END DO
!
	    IF(ALLOCATED(WRK_VEC))DEALLOCATE(WRK_VEC)
	    IF(ALLOCATED(INDX))DEALLOCATE(INDX)
	    K=N_TRANS+N_NT_TRANS+N_COL_TRANS
	    ALLOCATE (WRK_VEC(K),INDX(K))
	    WRK_VEC(1:N_TRANS)=ABS(NEW_RATES(1:N_TRANS,DPTH_INDX))
	    IF(N_NT_TRANS .GT. 0)THEN
	      WRK_VEC(N_TRANS+1:N_TRANS+N_NT_TRANS-1)=ABS(NEW_NT_RATES(1:N_NT_TRANS,DPTH_INDX))
	    END IF
	    IF(N_COL_TRANS .GT. 0)THEN
	      WRK_VEC(N_TRANS+N_NT_TRANS+1:K)=ABS(NEW_COL_RATES(1:N_COL_TRANS,DPTH_INDX))
	    END IF
	    CALL INDEXX(K,WRK_VEC,INDX,L_FALSE)
!
	    YV=0.0_LDP
	    DO ML=1,N_TRANS_LIM
	      L=INDX(ML)
	      IF(L .LE. N_TRANS)THEN
	        DPTH_VEC=NEW_RATES(L,:)/NORM_VEC(:)
	        J=INDEX(TRANS_NAME(LINK(L)),'-')
                IF( INDEX(TRANS_NAME(LINK(L)),TRIM(LEVEL)) .GT.  J)DPTH_VEC=-DPTH_VEC
	        CURVE_LAB=TRANS_NAME(LINK(L))
	      ELSE IF(L .LT. N_TRANS+N_NT_TRANS)THEN
	        L=L-N_TRANS
	        DPTH_VEC=NEW_NT_RATES(L,:)/NORM_VEC(:)
	        I=NT_LINK(L)
	        J=INDEX(NT_TRANS_NAME(I),'-')
                IF( INDEX(NT_TRANS_NAME(I),TRIM(LEVEL)) .GT.  J)DPTH_VEC=-DPTH_VEC
	        CURVE_LAB=TRIM(NT_TRANS_NAME(I))//'-NT'
	      ELSE
	        L=L-N_TRANS-N_NT_TRANS
	        DPTH_VEC=NEW_COL_RATES(L,:)/NORM_VEC(:)
	        I=COL_LINK(L)
	        J=INDEX(COL_TRANS_NAME(I),'-')
                IF( INDEX(COL_TRANS_NAME(I),TRIM(LEVEL)) .GT.  J)DPTH_VEC=-DPTH_VEC
	        CURVE_LAB=TRIM(COL_TRANS_NAME(I))//'-C'
	      END IF
	      YV=YV+DPTH_VEC
	      IF(MAXVAL(ABS(DPTH_VEC(1:ND))) .GT. RATE_THRESHOLD)THEN
	        WRITE(6,'(1X,A,T40,ES14.4)')TRIM(CURVE_LAB),DPTH_VEC(DPTH_INDX)
	        J=LEN_TRIM(SPECIES); CURVE_LAB=CURVE_LAB(J+2:)
	        J=LEN_TRIM(CURVE_LAB)
	        IF(CURVE_LAB(J:J) .EQ. 'C')THEN
	          CURVE_LAB(J-2:)=CURVE_LAB(J-1:)
	        ELSE IF(CURVE_LAB(J:J) .EQ. 'T')THEN
	          CURVE_LAB(J-3:)=CURVE_LAB(J-2:)
	        ELSE
	          CURVE_LAB(J:J)=' '
	        END IF
	        CALL DP_CURVE_LAB(ND,XV,DPTH_VEC,CURVE_LAB)
	      END IF
	    END DO
!
	    DPTH_VEC=-REC_RATE(SL_INDX,1:ND)/NORM_VEC(1:ND)
	    IF(MAXVAL(ABS(DPTH_VEC(1:ND))) .GT. RATE_THRESHOLD)THEN
	      WRITE(6,'(1X,A,T40,ES14.4)')'Recombination rate',ZV(DPTH_INDX)
	      CURVE_LAB='Rec. rate'
	      CALL DP_CURVE_LAB(ND,XV,DPTH_VEC,CURVE_LAB)
	    END IF
	    YV=YV+DPTH_VEC
!
	    DPTH_VEC=PHOT_RATE(SL_INDX,1:ND)/NORM_VEC(1:ND)
	    IF(MAXVAL(ABS(DPTH_VEC(1:ND))) .GT. RATE_THRESHOLD)THEN
	      WRITE(6,'(1X,A,T40,ES14.4)')'Photoionization rate',ZV(DPTH_INDX)
	      CURVE_LAB='Phot. rate'
	      CALL DP_CURVE_LAB(ND,XV,DPTH_VEC,CURVE_LAB)
	    END IF
	    YV=YV+DPTH_VEC
!
	    CALL DP_CURVE_LAB(ND,XV,YV,'Total net rate')
	    YLABEL='Normalized rates'
!
	    DEALLOCATE(INDX,WRK_VEC)
!
! Examine the rates at a specific depth.
!
	  ELSE IF(UC(PLT_OPT) .EQ. 'RATE')THEN
	    IF(N_TRANS .EQ. 0)THEN
	      WRITE(6,*)'Error -- no transitions to examine'
	      WRITE(6,*)'Use SPEC and LEVELS options'
	      GOTO 1000
	    END IF
!
	    N_TRANS_LIM=10
	    CALL GEN_IN(N_TRANS_LIM,'Transition limit')
!
	    CALL DERIVCHI(YV,XV,R,ND,'LINEAR')
	    WRITE(6,*)'Done deriv'
!
	    DO I=1,ND
	      NEW_RATES(1:N_TRANS,I)=NEW_RATES(1:N_TRANS,I)*R(I)*R(I)/XV(I)
	      IF(DO_NT_RATES)NEW_NT_RATES(:,I)=NEW_NT_RATES(:,I)*R(I)*R(I)/XV(I)
	    END DO
	    SUM_RATES=0.0_LDP
	    DO I=1,ND
	      DO ML=1,N_TRANS
	        SUM_RATES(ML)=SUM_RATES(ML)+ABS(NEW_RATES(ML,I))
	       END DO
	    END DO
	    WRITE(6,*)'Set sum rates'
	    WRITE(6,*)'Rates have bee corrupted - use LEVEL option to reset rates'
	    WRITE(6,*)' '
!
	    DO ML=1,MIN(N_TRANS_LIM,N_TRANS)
	      J=MAXLOC(SUM_RATES(1:N_TRANS),1)
	      DPTH_VEC=NEW_RATES(J,:)
	      K=INDEX(TRANS_NAME(LINK(J)),'-')
              IF( INDEX(TRANS_NAME(LINK(J)),TRIM(LEVEL)) .GT.  K )DPTH_VEC=-DPTH_VEC
	      CALL DP_CURVE(ND,XV,DPTH_VEC)
              WRITE(6,'(1X,A)')TRIM(TRANS_NAME(LINK(J)))
	      SUM_RATES(J)=0.0_LDP
	    END DO
	    WRITE(6,*)'Done line curves'
!
	    IF(DO_AUTO_RATES)THEN
	      T1=0.0_LDP
	      YV=0.0_LDP; ZV=0.0_LDP
	      DO J=1,N_AUTO
	         IF(INDEX(AUTO_LEV_NAME(J),TRIM(LEVEL)) .NE. 0)THEN
	           YV=YV+AUTO_REC_RATE(J,:)
	           ZV=ZV+AUTO_RATE(J,:)
	         END IF
	      END DO
	      CALL DP_CURVE(ND,XV,YV)
1	      CALL DP_CURVE(ND,XV,ZV)
	      WRITE(6,'(1X,A,T30,ES12.3)')'Maximum auto/anti autoionization rate',T1
	    END IF
	    YLABEL='Normalized origin'
!
! Plots actual rates at a specified depth as a function of lambda.
! Only done for line trasitions.
!
	  ELSE IF(UC(PLT_OPT) .EQ. 'AEXD')THEN
	    IF(SPECIES .EQ. ' ')THEN
	       WRITE(6,*)'Need to use SPECIES option first to set the species'
	       GOTO 1000
	    END IF
	    DPTH_INDX=ND/2
	    CALL GEN_IN(DPTH_INDX,'Depth to be examined')
	    IF(.NOT. ALLOCATED(LONG_XV))THEN
	       ALLOCATE(LONG_XV(N_LINES),LONG_YV(N_LINES))
	    END IF
	    N_TRANS=0
	    DO ML=1,N_LINES
	       IF(INDEX(TRANS_NAME(ML),TRIM(SPECIES)) .NE. 0)THEN
	         N_TRANS=N_TRANS+1; I=N_TRANS
	         LONG_XV(I)=LAM(ML)
	         LONG_YV(I)=RATES(ML,DPTH_INDX)
	       END IF
	    END DO
	    CALL DP_CURVE(N_TRANS,LONG_XV,LONG_YV)
	    XLABEL_SAVE=XLABEL
	    XLABEL='\gl(\A)'
	    YLABEL='Rate'
!
! Examine the rates at a specific depth.
!
	  ELSE IF(UC(PLT_OPT) .EQ. 'EXD')THEN
	    IF(N_TRANS .EQ. 0)THEN
	      WRITE(6,*)'Error -- no transitions to examine.'
	      WRITE(6,*)'Use SPEC and LEVELS options first.'
	      GOTO 1000
	    END IF
!
	    DPTH_INDX=ND/2
	    CALL GEN_IN(DPTH_INDX,'Depth to be examined')
	    N_TRANS_LIM=30
	    CALL GEN_IN(N_TRANS_LIM,'Limits the maximum number of transition written')
!
	    TERM_OUTPUT=.TRUE.
	    CALL GEN_IN(TERM_OUTPUT,'Terminal output (otherwise file)?')
	    IF(TERM_OUTPUT)THEN
	      LU_TMP=6
	    ELSE
	      CALL GET_LU(LU_TMP,'EXD option')
	      IF(ND .GT. 999)THEN 
	        WRITE(FILENAME,'(I4.4)')DPTH_INDX
	      ELSE IF(ND .GT. 99)THEN
	        WRITE(FILENAME,'(I3.3)')DPTH_INDX
	      ELSE
	        WRITE(FILENAME,'(I2.2)')DPTH_INDX
	      END IF
	      FILENAME='EXD_'//TRIM(SPECIES)//'_'//TRIM(FILENAME)
	      OPEN(UNIT=LU_TMP,FILE=FILENAME,STATUS='UNKNOWN',ACCESS='APPEND')
	      WRITE(6,*)'Data output (appended if already exists) to ',TRIM(FILENAME)
	    END IF
!
	    WRITE(LU_TMP,*)' '
	    WRITE(LU_TMP,*)'Data for '//TRIM(SPECIES)//'('//TRIM(LEVEL),') at depth',DPTH_INDX
	    WRITE(LU_TMP,*)' '
!
! We first deduce the maximum rate to each level.
!
	    MAX_RATE=0.0_LDP
	    IF(DO_AUTO_RATES)THEN
	      T1=0.0_LDP
	      DO J=1,N_AUTO
	         IF(INDEX(AUTO_LEV_NAME(J),TRIM(LEVEL)) .NE. 0)THEN
	           T1=MAX(T1,AUTO_REC_RATE(J,DPTH_INDX))
	           T1=MAX(T1,AUTO_RATE(J,DPTH_INDX))
	         END IF
	      END DO
	      WRITE(LU_TMP,'(1X,A,T30,ES12.3)')'Maximum auto/anti autoionization rate',T1
	      MAX_RATE=T1
	    END IF
!
	    T1=REC_RATE(SL_INDX,DPTH_INDX)
	    T1=MAX(T1,PHOT_RATE(SL_INDX,DPTH_INDX))
	    WRITE(LU_TMP,'(/,1X,A,T30,A,ES12.4)')'Maximum phot/rec. rate ','=',T1
	    MAX_RATE=MAX(MAX_RATE,T1)
!
	    I=DPTH_INDX; T1=0D0
	    DO ML=1,N_TRANS
	       T1=MAX(T1,ABS(NEW_RATES(ML,I)))
	    END DO
	    MAX_RATE=MAX(MAX_RATE,T1)
	    WRITE(LU_TMP,'(1X,A,T30,A,ES12.4)')'Maximum line rate ','=',T1
!
	    IF(N_NT_TRANS .NE. 0)THEN
	      DO ML=1,N_NT_TRANS
	        T1=MAX(T1,ABS(NEW_NT_RATES(ML,I)))
	      END DO
	      WRITE(LU_TMP,'(1X,A,T30,A,ES12.4)')'Maximum NT rate ','=',T1
	      MAX_RATE=MAX(MAX_RATE,T1)
	    END IF
!
	    IF(N_COL_TRANS .NE. 0)THEN
	      T1=0.0_LDP
	      DO ML=1,N_COL_TRANS
	        T1=MAX(T1,ABS(NEW_COL_RATES(ML,I)))
	      END DO
	      WRITE(LU_TMP,'(1X,A,T30,A,ES12.4)')'Maximum collison rate ','=',T1
	      MAX_RATE=MAX(MAX_RATE,T1)
	    END IF
	    MAX_RATE=MAX(MAX_RATE,REC_RATE(SL_INDX,I))
	    MAX_RATE=MAX(MAX_RATE,PHOT_RATE(SL_INDX,I))
!
	    WRITE(LU_TMP,*)' '
	    WRITE(LU_TMP,'(2X,A,I5)')         '      Depth index =',DPTH_INDX
	    WRITE(LU_TMP,'(2X,A,ES16.8,2X,A)')'     R(DPTH_INDX) =',R(DPTH_INDX),'(10^{10} cm)'
	    WRITE(LU_TMP,'(2X,A,ES16.8,2X,A)')'     V(DPTH_INDX) =',V(DPTH_INDX),'(km/s)'
	    WRITE(LU_TMP,'(2X,A,ES16.8,2X,A)')'     T(DPTH_INDX) =',T(DPTH_INDX),'(10^4 K)'
	    WRITE(LU_TMP,'(2X,A,ES16.8,2X,A)')'    ED(DPTH_INDX) =',ED(DPTH_INDX),'(cm^{-3})'
	    WRITE(LU_TMP,'(2X,A,ES16.8,2X,A)')'  Scaling Factor  =',MAX_RATE,'(cm^{-3} s^{-1})'
	    WRITE(LU_TMP,*)' '
	    WRITE(LU_TMP,*)'For line and non-thermal (NT) transitions, the upper level is listed'
	    WRITE(LU_TMP,*)'    first in the transition name.'
	    WRITE(LU_TMP,*)'A negative rate imples there is a net flow to the upper state'
	    WRITE(LU_TMP,*)' '
	    WRITE(LU_TMP,*)'Transitions draining ',TRIM(LEVEL)
	    WRITE(LU_TMP,*)' '
!
	    ALLOCATE (WRK_VEC(N_TRANS),INDX(N_TRANS))
	    WRK_VEC(:)=ABS(NEW_RATES(1:N_TRANS,DPTH_INDX))
	    CALL INDEXX(N_TRANS,WRK_VEC,INDX,L_TRUE)
!
	    WRITE(TMP_STR,'(I2.2)')MAX_TRANS_LENGTH+6
	    TMP_FMT='(1X,A,T'//TMP_STR(1:2)//',36X,F9.5)'
	    TMP_STR='(1X,A,T'//TMP_STR(1:2)//',ES15.4E5,3X,ES15.4E5,3X,F9.5)'
!
	    T1=0.0_LDP; I=DPTH_INDX; CNT=0
	    DO ML=N_TRANS,1,-1
	      L=INDX(ML)
	      IF( (INDEX(TRANS_NAME(LINK(L)),'('//TRIM(LEVEL)) .NE.  0 .AND. NEW_RATES(L,I) .GT. 0) .OR.
	1         (INDEX(TRANS_NAME(LINK(L)),'('//TRIM(LEVEL)) .EQ.  0 .AND. NEW_RATES(L,I) .LT. 0) )THEN
	        T1=T1+ABS(NEW_RATES(L,I))
	        CNT=CNT+1
	        T2=NEW_RATES(L,I)/MAX_RATE
	        IF(CNT .LE. N_TRANS_LIM .AND. ABS(T2) .GT.  0.00001_LDP)THEN
                  WRITE(LU_TMP,TMP_STR)TRIM(TRANS_NAME(LINK(L))),T2,LAM(LINK(L)),T1/MAX_RATE
	        END IF
	      END IF
	    END DO
	    WRITE(LU_TMP,TMP_FMT)'Total',T1/MAX_RATE
!
	    WRITE(LU_TMP,*)' '
	    WRITE(LU_TMP,*)'Transitions populating ',TRIM(LEVEL)
	    WRITE(LU_TMP,*)' '
	    T1=0.0_LDP; CNT=0.0_LDP
	    DO ML=N_TRANS,1,-1
	      L=INDX(ML)
	      IF( (INDEX(TRANS_NAME(LINK(L)),'('//TRIM(LEVEL)) .EQ.  0 .AND. NEW_RATES(L,I) .GT. 0) .OR.
	1         (INDEX(TRANS_NAME(LINK(L)),'('//TRIM(LEVEL)) .NE.  0 .AND. NEW_RATES(L,I) .LT. 0) )THEN
	        T1=T1+ABS(NEW_RATES(L,I))
	        CNT=CNT+1
	        T2=NEW_RATES(L,I)/MAX_RATE
	        IF(CNT .LE. N_TRANS_LIM .AND. ABS(T2) .GT. 0.00001_LDP)THEN
                  WRITE(LU_TMP,TMP_STR)TRIM(TRANS_NAME(LINK(L))),T2,LAM(LINK(L)),T1/MAX_RATE
	        END IF
	      END IF
	    END DO
	    WRITE(LU_TMP,TMP_FMT)'Total',T1/MAX_RATE
!
	    IF(N_NT_TRANS .NE. 0)THEN
	      WRITE(LU_TMP,'(A)')' '
	      WRITE(LU_TMP,*)'Non thermal transitions'
	      DO ML=1,N_NT_TRANS
	        T1=NEW_NT_RATES(ML,I)/MAX_RATE
	        IF(ABS(T1) .GT. 0.00001_LDP)THEN
	          WRITE(LU_TMP,TMP_STR)TRIM(NT_TRANS_NAME(NT_LINK(ML))),T1
	        END IF
	      END DO
	    END IF
!
	    IF(N_COL_TRANS .NE. 0)THEN
	      WRITE(LU_TMP,'(A)')' '
	      WRITE(LU_TMP,*)'Collisional transitions'
	      DEALLOCATE (WRK_VEC,INDX)
	      ALLOCATE (WRK_VEC(N_COL_TRANS),INDX(N_COL_TRANS))
	      WRK_VEC(:)=ABS(NEW_COL_RATES(1:N_COL_TRANS,DPTH_INDX))
	      CALL INDEXX(N_COL_TRANS,WRK_VEC,INDX,L_TRUE)
!
	      DO ML=1,N_COL_TRANS
	        L=INDX(ML)
	        STRING=TRIM(SPECIES)//'('//TRIM(LEVEL)
	        IF( INDEX(COL_TRANS_NAME(COL_LINK(L)),TRIM(STRING)) .NE.  0)THEN
	          T1=-NEW_COL_RATES(ML,I)/MAX_RATE
	        ELSE
	          T1=NEW_COL_RATES(ML,I)/MAX_RATE
	        END IF
	        IF(ABS(T1) .GT. 0.00001_LDP)THEN
	          WRITE(LU_TMP,TMP_STR)TRIM(COL_TRANS_NAME(COL_LINK(ML))),T1
	        END IF
	      END DO
	    END IF
!
! Compute the total rates.
!    Note: SUM_FROM is the net rate where the level if interest is the UPPER level.
!    Note: SUM_TO   is the net rate where the level if interest is the LOWER level.
!
	    I=DPTH_INDX
	    SUM_TO=0.0_LDP
	    SUM_FROM=0.0_LDP
	    DO ML=1,N_TRANS
	      IF(INDEX(TRANS_NAME(LINK(ML)),'('//TRIM(LEVEL)) .NE.  0)THEN
	        SUM_FROM=SUM_FROM+NEW_RATES(ML,I)
	      ELSE
	        SUM_TO=SUM_TO+NEW_RATES(ML,I)
	      END IF
	    END DO
!
	    NT_SUM_TO=0.0_LDP
	    NT_SUM_FROM=0.0_LDP
	    DO ML=1,N_NT_TRANS
	      IF(INDEX(NT_TRANS_NAME(NT_LINK(ML)),'('//TRIM(LEVEL)) .NE.  0)THEN
	        NT_SUM_FROM=NT_SUM_FROM+NEW_NT_RATES(ML,I)
	      ELSE
	        NT_SUM_TO=NT_SUM_TO+NEW_NT_RATES(ML,I)
	      END IF
	    END DO
!
	    COL_SUM_TO=0.0_LDP
	    COL_SUM_FROM=0.0_LDP
	    DO ML=1,N_COL_TRANS
	      IF(INDEX(COL_TRANS_NAME(COL_LINK(ML)),'('//TRIM(LEVEL)) .NE.  0)THEN
	        COL_SUM_FROM=COL_SUM_FROM+NEW_COL_RATES(ML,I)
	      ELSE
	        COL_SUM_TO=COL_SUM_TO+NEW_COL_RATES(ML,I)
	      END IF
	    END DO
!
	    WRITE(LU_TMP,'(A)')
	    WRITE(LU_TMP,'(A,T42,A,ES14.4,F11.5)')' Net radiative rate from higher levels ','=',SUM_TO,SUM_TO/MAX_RATE
	    WRITE(LU_TMP,'(A,T42,A,ES14.4,F11.5)')' Net radiative rate into lower levels ','=',SUM_FROM,SUM_FROM/MAX_RATE
	    WRITE(LU_TMP,'(A,T42,A,ES14.4,F11.5)')' Net thermal rate into ',   '=',NT_SUM_TO,NT_SUM_TO/MAX_RATE
	    WRITE(LU_TMP,'(A,T42,A,ES14.4,F11.5)')' Net thermal rate out of ', '=',NT_SUM_FROM,NT_SUM_FROM/MAX_RATE
	    WRITE(LU_TMP,'(A,T42,A,ES14.4,F11.5)')' Collisional rate into ',   '=',COL_SUM_TO,COL_SUM_TO/MAX_RATE
	    WRITE(LU_TMP,'(A,T42,A,ES14.4,F11.5)')' Collsion rate out of ',    '=',COL_SUM_FROM,COL_SUM_FROM/MAX_RATE
	    IF(PR_SPECIES .EQ. SPECIES)THEN
	      WRITE(LU_TMP,'(A,T42,A,ES14.4,F11.5)')' Photoionization rate ','=',
	1             PHOT_RATE(SL_INDX,DPTH_INDX),PHOT_RATE(SL_INDX,DPTH_INDX)/MAX_RATE
	      WRITE(LU_TMP,'(A,T42,A,ES14.4,F11.5)')' Recombination rate ','=',
	1             REC_RATE(SL_INDX,DPTH_INDX),REC_RATE(SL_INDX,DPTH_INDX)/MAX_RATE
	      T2=REC_RATE(SL_INDX,DPTH_INDX)-PHOT_RATE(SL_INDX,DPTH_INDX)
	      WRITE(LU_TMP,'(A,T42,A,ES14.4,F11.5)')' Net recombination rate ','=',T2,T2/MAX_RATE
	    END IF
!
	    IF(DO_AUTO_RATES)THEN
	      NET_AUTO=0.0_LDP
	      DO J=1,N_AUTO
	        IF(INDEX(AUTO_LEV_NAME(J),TRIM(LEVEL)) .NE. 0)THEN
	          WRITE(LU_TMP,'(A)')' '
	          WRITE(LU_TMP,'(A,T42,ES14.4,F11.5,5X,A)')' Autoionizaton rate ','=',
	1                AUTO_RATE(J,DPTH_INDX), AUTO_RATE(J,DPTH_INDX)/MAX_RATE, AUTO_LEV_NAME(J)
	          WRITE(LU_TMP,'(A,T42,ES14.4,F11.5)')' Rec. autoionization rate ','=',
	1                AUTO_REC_RATE(J,DPTH_INDX), AUTO_REC_RATE(J,DPTH_INDX)/MAX_RATE
	          NET_AUTO=NET_AUTO+(AUTO_REC_RATE(J,DPTH_INDX)-AUTO_RATE(J,DPTH_INDX))
	          WRITE(LU_TMP,'(A,T42,ES14.4,F11.5)')' Net recom. autoionization rate ','=',NET_AUTO,NET_AUTO/MAX_RATE
	        END IF
	      END DO
	    END IF
	    WRITE(LU_TMP,'(A)')' '
	    T2=(SUM_TO-SUM_FROM)+(NT_SUM_TO-NT_SUM_FROM)+(COL_SUM_TO-COL_SUM_FROM)+
	1            (REC_RATE(SL_INDX,I)-PHOT_RATE(SL_INDX,I))+NET_AUTO
	    WRITE(LU_TMP,'(A,T42,A,ES14.4,F11.5)')' Net Rate (into if +ve)','=',T2,T2/MAX_RATE
!
	    WRITE(LU_TMP,'(/,A)')' The last # is the fractional rate and should be small (i.e, < 0.1)'
	    WRITE(LU_TMP,'(A)')  ' If 0.1, for example, 10% of the rate into a level is missing.'
	    WRITE(LU_TMP,'(A)')  ' It is important that all processes (including collisons) are included.'
	    WRITE(LU_TMP,'(A,/)')' Use COLL option in DISPGEN to generate collisional data'
	    IF(LU_TMP .NE. 6)CLOSE(LU_TMP)
	    DEALLOCATE(INDX,WRK_VEC)
!
	  ELSE IF(UC(PLT_OPT(1:1)) .EQ. 'H')THEN
!
	    WRITE(6,*)' '
	    WRITE(6,*)' XN:      Set X-axis to depth index'
	    WRITE(6,*)' XLOGR:   Set X-axis to Log R/R(ND)'
	    WRITE(6,*)' XLOGV:   Set X-axis to Log V'
	    WRITE(6,*)' XV:      Set X-axis to V'
	    WRITE(6,*)' '
	    WRITE(6,*)' SPEC:    Set the species and reads in SPCIES//PRRR file'
	    WRITE(6,*)' LID:     Lists level names'
	    WRITE(6,*)' LEVEL:   Sets the level to be studied'
	    WRITE(6,*)' AEXD:    Plot line rates at a single depth as a function of lambda'
	    WRITE(6,*)' EXD:     Examine rates at a single depth'
	    WRITE(6,*)' NORM:    Plot rates as a function of depth'
	    WRITE(6,*)' E(X):    Exit routine'
!
	  ELSE IF(UC(PLT_OPT) .EQ. 'P')THEN
	    CALL GRAMON_PGPLOT(XLABEL,YLABEL,' ',' ')
	    XLABEL=XLABEL_SAVE
!
	  ELSE IF(UC(PLT_OPT) .EQ. 'EX' .OR. UC(PLT_OPT) .EQ. 'E')THEN
	    STOP
	  END IF
1000	  CONTINUE
	END DO
!
200	WRITE(6,*)STRING
	STOP
!
	END
