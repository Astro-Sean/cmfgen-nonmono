!
! General program to read in the photoionization cross-sections for XzV.
! Dielectronic transitions can also be read in if necessary, and they
! are then treated with the photoionization cross-sections.
!
! NB: EDGE, XzV_LEVELNAME, GION_GS, NXzV refer to the FULL atom.
!
	SUBROUTINE RDPHOT_GEN_V2(EDGE,XzV_LEVELNAME,
	1             GION_GS,AT_NO_DUM,ZXzV,NXzV,
	1             XzV_ION_LEV_ID,N_PHOT,MAX_N_PHOT,
	1             XzSIX_PRES,EDGEXzSIX_F,GXzSIX_F,F_TO_S_XzSIX,
	1             XzSIX_LEVNAME_F,NXzSIX_F,
	1             SIG_GAU_KMS,FRAC_SIG_GAU,CUT_ACCURACY,ABOVE,
	1             X_RAYS,ID,DESC,LUIN,LUOUT)
	USE SET_KIND_MODULE
!
!  Data modules required.
!
	USE PHOT_DATA_MOD		!Contains all photoionization data
!
	IMPLICIT NONE
!
! Altered 24-May-2005 : DO_SMOOTHING was not being inityialized to TRUE when
!                         smoothing required.
! Altered 17-Apr-2004 : Photoionzation cross-sections can be smoothed
!                          as they are read in. SIG_GAU etc inserted
!                          into call. Changed to V2.
! Altered 28-Mar-2003 : Check ensuring that ID < NPD_MAX installed.
! Altered 03-Jul-2001 : Only compute NEF for levels below ionization limit.
!                         May need work.
! Altered 16-Dec-1999 : FILENAME nolonger set to upper case [for UNIX].
! Altered 01-Dec-1999 : Altered to handle alternate state names in PHOT file.
!                          Done becuase there is still some inconsistency in
!                          naming conventions. Alternate names are seperated by
!                          a "/". Spaces may appear before the slash. This
!                          was done to allow earlier versions of CMFGEN to read
!                          the same data file.
! Altered 12-Dec-1997 : Bug fix. DO_PHOT was not being properly set for
!                          B level photoionizations when B level not present.
!                          In this case require DO_PHOT(1,PHOT_ID)=.TRUE.
!                          [not PHOT_ID,1] Thus when computing level 1 cross-
!                          section, we include the cross-section for PHOT_ID.
!
! Altered  5-Dec-1996 : GEN_ASCI_OPEN used for opening ASCI files.
! Altered Sep-18-1996 : Dimension of LST_U, LST_LOC and LST_CROSS changed.
!                       Name of LST_U changed to LST_FREQ.
!
!
! Passed variables.
!
	INTEGER NXzV
	INTEGER N_PHOT
	INTEGER MAX_N_PHOT
	REAL(KIND=LDP) EDGE(NXzV)
	CHARACTER(LEN=40) XzV_LEVELNAME(NXzV)
!
! Statistical weight of resulting ion for PHOT_ID=1. Used only if
! ion is not present.
!
	REAL(KIND=LDP) GION_GS
	REAL(KIND=LDP) ZXzV
	REAL(KIND=LDP) AT_NO_DUM
	INTEGER XzV_ION_LEV_ID(MAX_N_PHOT)
	INTEGER ID
	CHARACTER(LEN=*) DESC
!
	REAL(KIND=LDP) SIG_GAU_KMS		!V[FWHM]=2.34 SIG_AU_KMS
	REAL(KIND=LDP) FRAC_SIG_GAU             !Spacing in smoothed cross section
	REAL(KIND=LDP) CUT_ACCURACY		!Accuracy to delete superflous data points
	LOGICAL ABOVE                  !Smooth use data only above ionization edge.
!
	LOGICAL XzSIX_PRES
	INTEGER NXzSIX_F
	REAL(KIND=LDP) EDGEXzSIX_F(NXzSIX_F)
	REAL(KIND=LDP) GXzSIX_F(NXzSIX_F)
	INTEGER F_TO_S_XzSIX(NXzSIX_F)
	CHARACTER*(*) XzSIX_LEVNAME_F(NXzSIX_F)
!
	LOGICAL X_RAYS
	LOGICAL, SAVE :: FIRST=.TRUE.
	INTEGER LU_SUM
	INTEGER LUIN,LUOUT
!
! Common block with opacity/emissivity constants.
!
	REAL(KIND=LDP) CHIBF,CHIFF,HDKT,TWOHCSQ,OPLIN,EMLIN
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/LINE/ OPLIN,EMLIN
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
	INTEGER, PARAMETER :: MAX_FILES=5
!
	INTEGER NTERMS(MAX_FILES)
	INTEGER N_DP(MAX_FILES)
	CHARACTER(LEN=40), ALLOCATABLE :: LOC_TERM_NAME(:)
!
	REAL(KIND=LDP), PARAMETER :: RZERO=0.0
	REAL(KIND=LDP), PARAMETER :: RONE=1.0
!
	REAL(KIND=LDP) LOC_EXC_FREQ_ION,TOTAL_WT
	REAL(KIND=LDP) T1,NU_INF
	REAL(KIND=LDP) SIG_REV_KMS
	REAL(KIND=LDP) SIG_SM
!
	INTEGER NUM_ROUTES
	INTEGER PHOT_ID,IOS
	INTEGER I,J,K,L
	INTEGER GRID_ID
	INTEGER N_HD,L1,L2,LN
	INTEGER NEXT,LUER
	INTEGER END_CROSS
	INTEGER NPNTS
	INTEGER NFILES
	CHARACTER(LEN=20) CROSS_UNIT
!
	LOGICAL DO_SMOOTHING
	LOGICAL SPLIT_J
!
	INTEGER MAX_HEAD
	PARAMETER (MAX_HEAD=15)
	CHARACTER(LEN=132) HEAD_STR(MAX_HEAD)
	CHARACTER(LEN=1) APPEND_CHAR
!
	CHARACTER(LEN=132) STRING
	CHARACTER(LEN=80) FILENAME
!
	INTEGER NPNTS_TOT
	INTEGER N_CROSS
	INTEGER, PARAMETER :: N_CROSS_MAX=1000000
	REAL(KIND=LDP) TMP_NU(N_CROSS_MAX)
	REAL(KIND=LDP) TMP_CROSS(N_CROSS_MAX)
!
	LUER=ERROR_LU()
!
	CALL GET_NUM_PHOT_ROUTES_V1(NFILES,NUM_ROUTES,
	1             NTERMS,N_DP,MAX_FILES,
	1             XzSIX_PRES,EDGEXzSIX_F,GXzSIX_F,F_TO_S_XzSIX,
	1             XzSIX_LEVNAME_F,NXzSIX_F,ID,DESC,LUIN,LUOUT)
!
! 
!
	ALLOCATE (PD(ID)%A_ID(NXzV,NFILES))
	ALLOCATE (PD(ID)%NEF(NXzV,NFILES))
!
! Now perform the dynamic allocation of memory for variables determined by the
! size of the input data.
!
	PD(ID)%MAX_CROSS=0
	PD(ID)%MAX_TERMS=0
	DO PHOT_ID=1,NFILES
	  PD(ID)%MAX_CROSS=PD(ID)%MAX_CROSS+N_DP(PHOT_ID)
	  PD(ID)%MAX_TERMS=MAX(PD(ID)%MAX_TERMS,NTERMS(PHOT_ID))
	END DO
!
	ALLOCATE (PD(ID)%ST_LOC(PD(ID)%MAX_TERMS,NFILES))
	ALLOCATE (PD(ID)%END_LOC(PD(ID)%MAX_TERMS,NFILES))
	ALLOCATE (PD(ID)%CROSS_TYPE(PD(ID)%MAX_TERMS,NFILES))
	ALLOCATE (PD(ID)%LST_LOC(NXzV,NFILES))
!
	ALLOCATE (LOC_TERM_NAME(PD(ID)%MAX_TERMS))
	ALLOCATE (PD(ID)%LST_CROSS(NXzV,NUM_ROUTES))
	ALLOCATE (PD(ID)%LST_FREQ(NXzV,NUM_ROUTES))
!
! Perform default initializations.
!
	PD(ID)%A_ID(:,:)=0
	PD(ID)%DO_KSHELL_W_GS=.FALSE.
	PD(ID)%AT_NO=AT_NO_DUM
	PD(ID)%ZION=ZXzV
!
! END_CROSS provides an indication of the last storage location used. It is
! initialized here, since all photoionization routes use the same storage
! locations.
!
	END_CROSS=0
!
! If the XzSIX ground state is not present (hence XzSIX_PRES is .FALSE.) we
! include the X-RAY opacity with the regular photoionization cross-section.
! We also include the K-shell cross-section for lithium ions, since
! only one electron is ejected. Strictly speaking, should have 1s2.nl going
! to 1s nl state (rather than 1s^2) state of the HeI-like ion.
!
	IF(X_RAYS .AND. .NOT. XzSIX_PRES)THEN
	  PD(ID)%DO_KSHELL_W_GS=.TRUE.
	ELSE IF(NINT(AT_NO_DUM-ZXzV) .EQ. 2)THEN
	  PD(ID)%DO_KSHELL_W_GS=.TRUE.
	END IF
!
! 
! This section reads and stores the cross sections, looping over each final
! state. The final states should have logical name (files)
!
!               PHOTXzV_n where n=A is the ground state;
!                          where n=B, C etc denote excited states.
!
	N_CROSS=0
	NPNTS_TOT=0
	OPEN(UNIT=80,FORM='UNFORMATTED',STATUS='SCRATCH')
	DO PHOT_ID=1,NFILES
!
	  WRITE(LUOUT,'(A)')'  '
!
! Open file containing cross sections. We first read in all header info.
!
	  FILENAME='PHOT'//TRIM(DESC)//'_'
	  I=ICHRLEN(FILENAME)
	  APPEND_CHAR=CHAR( ICHAR('A')+(PHOT_ID-1) )
	  WRITE(FILENAME(I+1:I+1),'(A)')APPEND_CHAR
	  CALL GEN_ASCI_OPEN(LUIN,FILENAME,'OLD',' ','READ',IZERO,IOS)
!
! Find the record with the date the file was written. This begins the
! important information. NB: After the date the header records can be
! in any order. There must be NO BLANK lines however, as this signifies
! the start of the data.
!
	  L1=0
	  DO WHILE(L1 .EQ. 0)
	    READ(LUIN,'(A)',IOSTAT=IOS)STRING
	    L1=INDEX(STRING,'!Date')
	    IF(IOS .NE. 0)THEN

	      WRITE(LUER,*)'Error reading date from ',FILENAME
	      WRITE(LUER,*)'PHOT_ID=',PHOT_ID,'IOSTAT=',IOS
	      STOP
	    END IF
	    WRITE(LUOUT,'(A)')STRING(1:ICHRLEN(STRING))
	  END DO
!
	  N_HD=1
	  READ(LUIN,'(A)')HEAD_STR(N_HD)
	  WRITE(LUOUT,'(A)')HEAD_STR(N_HD)
	  DO WHILE(HEAD_STR(N_HD) .NE. '  ')
	    N_HD=N_HD+1
	    IF(N_HD .GT. MAX_HEAD)THEN
	      WRITE(LUER,*)'Insufficient header space in RDPHOT_GEN_V1: ',DESC
	      STOP
	    END IF
	    READ(LUIN,'(A)')HEAD_STR(N_HD)
	    WRITE(LUOUT,'(A)')HEAD_STR(N_HD)(1:ICHRLEN(HEAD_STR(N_HD)))
	  END DO
!
! Read in number of energy levels.
!
	  J=0
	  L1=0
	  DO WHILE(L1 .EQ. 0 .AND. J .LT. N_HD)
	    J=J+1
	    L1=INDEX(HEAD_STR(J),'!Number of energy levels')
	  END DO
	  IF(L1 .EQ. 0)THEN
	    WRITE(LUER,*)'Error in RDPHOT_GEN_V1: ',DESC
	    WRITE(LUER,*)'Number of energy levels string not found'
	    STOP
	  ELSE
	    READ(HEAD_STR(J),*,IOSTAT=IOS)NTERMS(PHOT_ID)
	    IF(IOS .NE. 0)THEN
	      WRITE(LUER,*)' Error reading NTERMS in RDPHOT_GEN_V1 - 2:',DESC
	      WRITE(LUER,*)' PHOT_ID=',PHOT_ID
	      STOP
	    END IF
	  END IF
!
! Read in units for cross-sections.
!
	  J=0
	  L1=0
	  DO WHILE(L1 .EQ. 0 .AND. J .LT. N_HD)
	    J=J+1
	    L1=INDEX(HEAD_STR(J),'!Cross-section unit')
	  END DO
	  IF(L1 .EQ. 0)THEN
	    WRITE(LUER,*)'Error in RDPHOT_GEN_V1: ',DESC
	    WRITE(LUER,*)'Unit for cross-sections not found'
	    STOP
	  ELSE
	    L1=INDEX(HEAD_STR(J),'  ')
	    CROSS_UNIT=HEAD_STR(J)(1:L1-1)
	    IF(CROSS_UNIT .NE. 'Megabarns')THEN
	      WRITE(LUER,*)' Error in RDPHOT_GEN_V1 - incorrect',
	1                        ' cross-section unit:',DESC
	      WRITE(LUER,*)'CROSS _UNIT=',CROSS_UNIT
	      WRITE(LUER,*)'UNIT should be Megabarns'
	      STOP
	    END IF
	  END IF
!
! Are photoioization cross-sections split for individual J states.
! This is done to faciliated determining level correspondance.
!
	  J=0
	  L1=0
	  DO WHILE(L1 .EQ. 0 .AND. J .LT. N_HD)
	    J=J+1
	    L1=INDEX(HEAD_STR(J),'!Split J levels')
	  END DO
	  IF(L1 .EQ. 0)THEN
	    WRITE(LUER,*)'Error in RDPHOT_GEN_V1: ',DESC
	    WRITE(LUER,*)'Split J levels indication not FOUND!'
	    STOP
	  ELSE
	    IOS=0
	    L1=INDEX(HEAD_STR(J),'  ')
	    READ(HEAD_STR(J)(1:L1-1),*,IOSTAT=IOS)SPLIT_J
	    IF(IOS .NE. 0)THEN
	      WRITE(LUER,*)'Error in RDPHOT_GEN_V1: ',DESC
	      WRITE(LUER,*)'Can''t read SPLIT_J value'
	    END IF
	  END IF
!
! See if data has already been smoothed.
!
	  J=0
	  L1=0
	  DO WHILE(L1 .EQ. 0 .AND. J .LT. N_HD)
	    J=J+1
	    L1=INDEX(HEAD_STR(J),'!Sigma used for Gaussian smoothing (km/s)')
	  END DO
	  SIG_REV_KMS=SIG_GAU_KMS
	  IF(L1 .EQ. 0)THEN
	      WRITE(LUER,*)' '
	      WRITE(LUER,*)'Warning in RDPHOT_GEN_V1: ',DESC
	      WRITE(LUER,*)'Sigma of smoothing Gaussian is unavailable'
	      WRITE(LUER,*)'Assuming data does not need to be smoothed'
	      DO_SMOOTHING=.FALSE.
	  ELSE
	    READ(HEAD_STR(J),*)SIG_SM
	    IF(SIG_SM .GE. 0.9_LDP*SIG_GAU_KMS)THEN
	      DO_SMOOTHING=.FALSE.
	    ELSE
	      DO_SMOOTHING=.TRUE.
	      SIG_REV_KMS=SQRT(SIG_GAU_KMS**2-SIG_SM**2)
	    END IF
	  END IF
!
! 
!
! Read in the cross-sections
!
	  DO J=1,NTERMS(PHOT_ID)
	    STRING=' '
	    DO WHILE( INDEX(STRING,'Configuration name') .EQ. 0)
	      READ(LUIN,'(A)',IOSTAT=IOS)STRING
	      IF(IOS .NE. 0)THEN
	        WRITE(LUER,*)'Error reading in configuration name from ',
	1                       FILENAME
	        WRITE(LUER,*)'Currently reading term',J
	        STOP
	      END IF
	    END DO
	    L1=INDEX(STRING,'!Configuration')-1
	    L1=ICHRLEN(STRING(1:L1))
	    LOC_TERM_NAME(J)=STRING(1:L1)
	    DO WHILE(LOC_TERM_NAME(J)(1:1) .EQ. ' ')
	      LOC_TERM_NAME(J)=LOC_TERM_NAME(J)(2:)
	    END DO
!
	    STRING=' '
	    DO WHILE( INDEX(STRING,'Type of cross-section') .EQ. 0)
	      READ(LUIN,'(A)',IOSTAT=IOS)STRING
	      IF(IOS .NE. 0)THEN
	        WRITE(LUER,*)'Error reading Type of cross-section from ',
	1         FILENAME
	        WRITE(LUER,*)'Currently reading term',J
	        STOP
	      END IF
	    END DO
	    READ(STRING,*)PD(ID)%CROSS_TYPE(J,PHOT_ID)
!
	    READ(LUIN,*)NPNTS
	    PD(ID)%ST_LOC(J,PHOT_ID)=END_CROSS+1
	    PD(ID)%END_LOC(J,PHOT_ID)=END_CROSS+NPNTS
	    END_CROSS=PD(ID)%END_LOC(J,PHOT_ID)
	    IF(PD(ID)%END_LOC(J,PHOT_ID) .GT. PD(ID)%MAX_CROSS)THEN
	      WRITE(LUER,*)'Error in RDPHOT_GEN_V1: ',DESC
	      WRITE(LUER,*)'MAX_CROSS is too small'
	      STOP
	    END IF
	    IF(PD(ID)%CROSS_TYPE(J,PHOT_ID) .EQ. 9 .AND. MOD(NPNTS,8) .NE. 0)THEN
	      WRITE(LUER,*)'Error reading incross-section for ',TRIM(FILENAME)
	      WRITE(LUER,*)'For TYPE=9, the number of points must be multiple of 8'
	      STOP
	    END IF
!
! We can either read in a tabulated data set, or parameters for
! fitting a function. For a tabulated data set there must be the frequency
! and cross-section on each line --- for the other cases one parameter
! is specified per line. For compatibility with the tabular form, a zero
! must be written after each parameter. NB: In the second case PD(ID)%NU_NORM is
! not a frequency, but rather a parameter.
!
! We output the cross-section to fort.80 so we can determined the number
! of crorss-section points for the levels that are in the model atom.
!
	    IF(PD(ID)%CROSS_TYPE(J,PHOT_ID) .LT. 20)THEN
	      READ(LUIN,*)(TMP_CROSS(K), K=1,NPNTS)
	      TMP_NU(1:NPNTS)=RZERO
	    ELSE
	      READ(LUIN,*)(TMP_NU(K),TMP_CROSS(K),K=1,NPNTS)
	      IF(DO_SMOOTHING)THEN
	        CALL SM_PHOT_V3(TMP_NU,TMP_CROSS,NPNTS,N_CROSS_MAX,
	1                          SIG_REV_KMS,FRAC_SIG_GAU,
	1                          CUT_ACCURACY,ABOVE)
	      END IF
	    END IF
	    WRITE(80)J,PHOT_ID,NPNTS
	    DO K=1,NPNTS
	      WRITE(80)TMP_NU(K),TMP_CROSS(K)
	    END DO
	    N_CROSS=N_CROSS+1
	    NPNTS_TOT=NPNTS_TOT+NPNTS
	  END DO
!	END DO
!
! Match and identify input cross sections to levels in code.
!
! Match each of the TERMS with the corresponding levels in the full atom.
! We ensure that each level has a photoionization cross-section.
!
! To get the Rydberg constant we assume that the atomic mass in AMU is just
! twice the atomic number. Only exeption is H.
!
!	DO K=1,NFILES
	  K=PHOT_ID
	  DO I=1,NXzV
	    IF(SPLIT_J)THEN
	      L1=ICHRLEN(XzV_LEVELNAME(I))
	    ELSE
	      L1=INDEX(XzV_LEVELNAME(I),'[')-1
	      IF(L1 .LE. 0)L1=ICHRLEN(XzV_LEVELNAME(I))
	    END IF
	    DO J=1,NTERMS(K)
	      IF(XzV_LEVELNAME(I)(1:L1) .EQ. LOC_TERM_NAME(J))THEN
	        IF(PD(ID)%A_ID(I,K) .EQ. 0)THEN
	          PD(ID)%A_ID(I,K)=J
	        ELSE
	          WRITE(LUER,*)'Error in handling photoionization ',
	1                'cross-sections for ',FILENAME
	          WRITE(LUER,*)'Non-unique cross-section for ',
	1                        XzV_LEVELNAME(I)
	          STOP
	        END IF
	      END IF
	    END DO
	    IF(PD(ID)%A_ID(I,K) .EQ. 0)THEN
	      WRITE(LUER,*)'Error in handling photoionization ',
	1                   'cross-sections for ',FILENAME
	      WRITE(LUER,*)'Cross section unavailable for level ',
	1                      XzV_LEVELNAME(I)
	      STOP
	    END IF
	  END DO
	END DO
!
! This simple evaluation of NEFF will be wrong for many states (since it depnds on the final state)
! however it should be right for the cross-sections where Neff is used.
!
	NU_INF=1.0E-15_LDP*109737.31_LDP*SPEED_OF_LIGHT()
	IF(PD(ID)%AT_NO .EQ. 1)THEN
	  NU_INF=NU_INF/(RONE+5.48597E-04_LDP)
	ELSE
	  NU_INF=NU_INF/(RONE+5.48597E-04_LDP/(2*PD(ID)%AT_NO))
	END IF
	PD(ID)%NEF(1,1:NFILES)=-100.0_LDP
	DO K=1,NUM_ROUTES
	  GRID_ID=PD(ID)%PHOT_GRID(K)
	  IF(PD(ID)%NEF(1,GRID_ID) .LT. -1)THEN
	    DO I=1,NXzV
	      T1=NU_INF/(EDGE(I)+PD(ID)%EXC_FREQ(K))
	      PD(ID)%NEF(I,GRID_ID)=RZERO
	      IF(T1 .GT. RZERO)PD(ID)%NEF(I,GRID_ID)=ZXzV*SQRT(T1)
	    END DO
	  END IF
	END DO
!
	REWIND(UNIT=80)
	PD(ID)%MAX_CROSS=NPNTS_TOT
	ALLOCATE (PD(ID)%NU_NORM(PD(ID)%MAX_CROSS))
	ALLOCATE (PD(ID)%CROSS_A(PD(ID)%MAX_CROSS))
	L=0
	DO I=1,N_CROSS
	  READ(80)J,PHOT_ID,NPNTS
	  DO K=1,NPNTS
	    READ(80)PD(ID)%NU_NORM(L+K),PD(ID)%CROSS_A(L+K)
	  END DO
	  PD(ID)%ST_LOC(J,PHOT_ID)=L+1
	  PD(ID)%END_LOC(J,PHOT_ID)=L+NPNTS
	  L=L+NPNTS
	END DO
	CLOSE(UNIT=80)
!
! ALPHA_BF is the constant needed to evaluate the photoionization cross-section
! for a hydrogenic ion. It is in the correct units for R in units of
! 10^10 cm, and CHI.R dimensionless.
!
	PD(ID)%ALPHA_BF=2.815E-06_LDP*(PD(ID)%ZION)**4
!
	IF(FIRST)THEN
	  CALL GET_LU(LU_SUM,'LU_SUM in RDPHOTSUM_V2')
	  OPEN(UNIT=LU_SUM,FILE='PHOT_SUMMARY_FILE',STATUS='UNKNOWN',ACTION='WRITE')
	  WRITE(LU_SUM,*)'Summary of Photoionization routes for'
	  WRITE(LU_SUM,'(A,T10,A,T40,3X,A,3X,A,4X,A,4X,A,2X,A,7X,A)')
	1      'Desc','Level name','SL','G(ion)','Exc. Freq.','Scale Fac','Grid','DP_PHOT(I,:)'
	  FIRST=.FALSE.
	ELSE
	  CALL GET_LU(LU_SUM,'LU_SUM in RDPHOTSUM_V2')
	  OPEN(UNIT=LU_SUM,FILE='PHOT_SUMMARY_FILE',STATUS='OLD',ACTION='WRITE',POSITION='APPEND')
	END IF
!
	IF(NUM_ROUTES .GT. MAX_N_PHOT)THEN
	  WRITE(LUER,*)'Error in RDPHOT_GEN_V2 -- MAX_N_PHOT needs to be increased.'
	  WRITE(LUER,*)'MAX_N_PHOT =',MAX_N_PHOT
	  WRITE(LUER,*)'Species is ',TRIM(DESC)
	  WRITE(LUER,*)'NUM_ROUTES =',NUM_ROUTES
	  STOP
	END IF
!
	DO K=1,NUM_ROUTES
	  L=PD(ID)%ION_LEV_ID(K)
	  XzV_ION_LEV_ID(K)=PD(ID)%ION_LEV_ID(K)
	  J=0
	  IF(XzSIX_PRES)THEN
	    J=F_TO_S_XzSIX(L)
	  ELSE
	    WRITE(LU_SUM,'(A)')' '
	  END IF
	  WRITE(LU_SUM,'(1X,A,T10,A,T40,I5,3X,F6.1,3X,ES12.4,3X,ES10.3,3X,I3,3X,20L5)')
	1      TRIM(DESC),TRIM(PD(ID)%FINAL_ION_NAME(K)),J,
	1      PD(ID)%GION(K),PD(ID)%EXC_FREQ(K),PD(ID)%SCALE_FAC(K),
	1      PD(ID)%PHOT_GRID(K),(PD(ID)%DO_PHOT(K,I),I=1,NUM_ROUTES)
	END DO
	CLOSE(LU_SUM)
!
! PD(ID)%NUM_PHOT_ROUTES is used to refer to the actual number of photoionization
! routes available. If the final state is not available, they go to the
! ground state.
!
	PD(ID)%NUM_PHOT_ROUTES=NUM_ROUTES
	PD(ID)%N_FILES=NFILES
!
! Set N_PHOT to the actual number of final states that are being used.
! This is returned in the call.
!
	N_PHOT=0
	DO I=1,PD(ID)%NUM_PHOT_ROUTES
	  IF(PD(ID)%DO_PHOT(I,I))N_PHOT=N_PHOT+1
	END DO
	GION_GS=PD(ID)%GION(1)
!
	DEALLOCATE (LOC_TERM_NAME)
!
! Initialize storage locations for use by PHOT_XzV. These are used to store
! photoionization cross-sections to speed up computation.
!
	PD(ID)%LST_CROSS(:,:)=RZERO
	PD(ID)%LST_FREQ(:,:)=RZERO
	DO J=1,PD(ID)%NUM_PHOT_ROUTES
	  L=PD(ID)%PHOT_GRID(J)
	  DO I=1,NXzV			!Loop over levels
	    K=PD(ID)%A_ID(I,L)			!Get pointer to cross-section
	    PD(ID)%LST_LOC(I,L)=(PD(ID)%ST_LOC(K,L) + PD(ID)%END_LOC(K,L))/2
	  END DO
	END DO
!
	RETURN
	END
