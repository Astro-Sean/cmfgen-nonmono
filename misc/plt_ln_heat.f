	PROGRAM PLT_LN_HEAT
	USE SET_KIND_MODULE
	USE GEN_IN_INTERFACE
	USE MOD_COLOR_PEN_DEF
	IMPLICIT NONE
!
! Altered: 30-Aug-2018 - All important arrays now allocatable.
!                          Cleaned up output from option MAIN.
!                          Added option IDI, and options to set X-axis.
!                          Options now done in same way as maingen.
!                          Changes done over a week (23 ? on).
! Altered:  6-Dec-2017 - Made copatible with version from osiris.
! Altered: 14-Jun-2017 - Added option to read LINEHEAT created by SOBOLEV approximation.
!
	INTEGER, PARAMETER :: NMAX=600000
	INTEGER, PARAMETER :: NS_MAX=400
!
	REAL(KIND=LDP), ALLOCATABLE :: NU(:)
	REAL(KIND=LDP), ALLOCATABLE :: LAM(:)
	REAL(KIND=LDP), ALLOCATABLE :: SCALE_FAC(:)
	REAL(KIND=LDP), ALLOCATABLE :: XV(:)
	REAL(KIND=LDP), ALLOCATABLE :: Y(:)
	INTEGER, ALLOCATABLE:: INDX(:)
	CHARACTER(LEN=60), ALLOCATABLE :: NAME(:)
!
! These are used to determine which species contribute most to the
! difference between SE_SCL and SE_NOSCL.
!
	REAL(KIND=LDP) SUMD(NS_MAX)
	CHARACTER*10 SPEC(NS_MAX)
!
	REAL(KIND=LDP), ALLOCATABLE :: LH(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: SE_SCL(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: SE_NOSCL(:,:)
	CHARACTER(LEN=120), ALLOCATABLE ::  TRANS_INFO(:)
!
	REAL(KIND=LDP) T1,T2
	REAL(KIND=LDP) CUR_SUM
	INTEGER ND
	INTEGER N_LINES
	INTEGER COUNT
	INTEGER I,J,K,L,ML,IBEG,IDEPTH
	INTEGER NL,NUP
	INTEGER IOS
	INTEGER NSPEC
	LOGICAL FILE_OPEN
	LOGICAL DO_SORT
	LOGICAL SOB_MODEL
!
	CHARACTER(LEN=80) FILENAME
	CHARACTER(LEN=200) STRING
        CHARACTER(LEN=30) UC
	CHARACTER(LEN=20) PLT_OPT
	CHARACTER(LEN=20) XLABEL
	CHARACTER(LEN=20) YLABEL
        EXTERNAL UC
!
	WRITE(6,'(A)')' '
	WRITE(6,'(A)')' Program to plot the radiative equilibrium equation.'
	WRITE(6,'(A)')' Can be used to show how RE changes as we integrate from blue to red.'
	WRITE(6,'(A)')' Designed to see effect of scaling the heating rates.'
	WRITE(6,'(A)')' Can be used to identify SL assignments that might be changed.'
	WRITE(6,'(A)')' '
!
        SOB_MODEL=.FALSE.
        CALL GEN_IN(SOB_MODEL,'Was LINEHEAT created using the SOBOLEV approximation?')
100	CONTINUE
 	FILENAME='LINEHEAT'
	CALL GEN_IN(FILENAME,'File with data to be plotted')
	IF(FILENAME .EQ. ' ')GOTO 100
!
	OPEN(UNIT=20,FILE='MODEL',STATUS='OLD',IOSTAT=IOS)
	  IF(IOS .EQ. 0)THEN
	    DO WHILE(1 .EQ. 1)
	      READ(20,'(A)',IOSTAT=IOS)STRING
	      IF(IOS .NE. 0)EXIT
	      IF(INDEX(STRING,'!Number of depth points') .NE. 0)THEN
	         READ(STRING,*)ND
	         WRITE(6,'(A,I4)')' Number of depth points in the model is:',ND
	         EXIT
	      END IF
	    END DO
	  END IF
	INQUIRE(UNIT=20,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(UNIT=20)
!
	IF(IOS .NE. 0)THEN
	  WRITE(6,*)' Unable to open MODEL file to get # of depth points'
	  CALL GEN_IN(ND,'Number of data points (must be exact)')
	END IF
	N_LINES=NMAX
	CALL GEN_IN(N_LINES,'Maximum number of lines to be read')
!
	ALLOCATE (LH(ND,N_LINES))
	ALLOCATE (SE_SCL(ND,N_LINES))
	ALLOCATE (SE_NOSCL(ND,N_LINES))
	ALLOCATE (TRANS_INFO(N_LINES))
!
	ALLOCATE (NU(N_LINES))
	ALLOCATE (LAM(N_LINES))
	ALLOCATE (SCALE_FAC(N_LINES))
	ALLOCATE (XV(N_LINES))
	ALLOCATE (Y(N_LINES))
	ALLOCATE (INDX(N_LINES))
	ALLOCATE (NAME(N_LINES))
!
	OPEN(UNIT=11,FILE=FILENAME,STATUS='OLD',ACTION='READ')
!
	WRITE(6,'(A)')RED_PEN
	WRITE(6,'(A)')' Reading data -- this may take a while'
	WRITE(6,'(A)')DEF_PEN
	COUNT=0
	CALL TUNE(1,'READ')
	DO ML=1,N_LINES
	  STRING=' '
	  DO WHILE(STRING .EQ. ' ')
	    READ(11,'(A)',END=5000)STRING
	  END DO
	  IF(INDEX(STRING,'error in L due to') .NE. 0)GOTO 5000
	    IF(SOB_MODEL)THEN
	    READ(STRING,*,ERR=200)LAM(ML)
	    READ(11,*)LH(1:ND,ML)
	    COUNT=COUNT+1
	  ELSE
	    STRING=ADJUSTL(STRING)
	    TRANS_INFO(ML)=STRING
	    K=INDEX(STRING,' ')
	    STRING=ADJUSTL(STRING(K:))
	    K=INDEX(STRING,' ')
	    NAME(ML)=STRING(1:K-1)
	    K=INDEX(STRING,')')
	    DO WHILE(K .NE. 0)
	      STRING(1:)=STRING(K+1:)
	      K=INDEX(STRING,')')
	    END DO
	    K=INDEX(STRING,'  ')
	    STRING(1:)=STRING(K:)
	    READ(STRING,*,ERR=200)NU(ML),NL,NUP,SCALE_FAC(ML)
	    READ(11,*)LH(1:ND,ML)
	    READ(11,'(A)')STRING
	    IF(STRING .NE. ' ')THEN
	      WRITE(6,'(A)')' '
	      WRITE(6,*)'Error: invalid data format'
	      WRITE(6,*)'Check ND values'
	      WRITE(6,*)'Current line count is',COUNT
	      WRITE(6,'(A)')' '
	      STOP
	    END IF
	    READ(11,*)SE_SCL(1:ND,ML)
	    READ(11,*)SE_NOSCL(1:ND,ML)
	    COUNT=COUNT+1
	  END IF
	END DO
5000	CONTINUE
	CALL TUNE(2,'READ')
	CALL TUNE(3,'TERMINAL')

	N_LINES=COUNT
	WRITE(6,'(A)')BLUE_PEN
	WRITE(6,*)'Number of lines read is',N_LINES
	WRITE(6,'(A)')DEF_PEN
!
        IF(SOB_MODEL)THEN
          NU(1:N_LINES)=2.99792458E+03_LDP/LAM(1:N_LINES)
        ELSE
          LAM(1:N_LINES)=2.99792458E+03_LDP/NU(1:N_LINES)
	END IF
!
! Set def axis.
!
	XV(1:N_LINES)=NU(1:N_LINES)
	XLABEL='\gn(10\u15 \dHz)'
!
	DO WHILE(1 .EQ. 1)
	  WRITE(6,*)' '
	  PLT_OPT='P'
	  CALL GEN_IN(PLT_OPT,'Plot option')
!
	  IF(UC(PLT_OPT) .EQ. 'XLAM')THEN
	    DO I=1,N_LINES
	      XV(I)=LAM(I)
	    END DO
	    XLABEL='\gl(\AA)'
	  ELSE IF(UC(PLT_OPT) .EQ. 'XNU')THEN
	    DO I=1,N_LINES
	      XV(I)=NU(I)
	    END DO
	    XLABEL='\gn(10\u15 \dHz)'
	  ELSE IF(UC(PLT_OPT) .EQ. 'XN')THEN
	    DO I=1,N_LINES
	      XV(I)=I
	    END DO
	    XLABEL='I'
!
	  ELSE IF(UC(PLT_OPT) .EQ. 'LH')THEN
	    K=ND
	    CALL GEN_IN(K,'Depth for plotting')
	    Y(1:N_LINES)=LH(K,1:N_LINES)
	    CALL DP_CURVE(N_LINES,XV,Y)
	    YLABEL='LH'
!
	  ELSE IF(UC(PLT_OPT) .EQ. 'SS')THEN
	    K=ND
	    CALL GEN_IN(K,'Depth for plotting')
	    Y(1:N_LINES)=SE_SCL(K,1:N_LINES)
	    WRITE(6,*)Y(N_LINES)
	    CALL DP_CURVE(N_LINES,XV,Y)
	    YLABEL='STEQ(scaled)'
!
	  ELSE IF(UC(PLT_OPT) .EQ. 'SN')THEN
	    K=ND
	    CALL GEN_IN(K,'Depth for plotting')
	    Y(1:N_LINES)=SE_NOSCL(K,1:N_LINES)
	    CALL DP_CURVE(N_LINES,XV,Y)
	    YLABEL='STEQ(not scaled)'
!
	  ELSE IF(UC(PLT_OPT) .EQ. 'DIFF')THEN
	    K=ND
	    CALL GEN_IN(K,'Depth for plotting')
	    DO L=1,N_LINES-1
	      Y(L)=(SE_NOSCL(K,L+1)-SE_NOSCL(K,L))-(SE_SCL(K,L+1)-SE_SCL(K,L))
	    END DO
	    Y(N_LINES)=0
	    CALL DP_CURVE(N_LINES,XV,Y)
	    YLABEL='Diff'
!
	  ELSE IF(UC(PLT_OPT) .EQ. 'MAIN')THEN
!
! Get the influence due to scaling of a given line on the radiative equilibrium
! equation.
!
	    IDEPTH=ND
	    CALL GEN_IN(IDEPTH,'Depth for plotting')
	    SUMD=0.0_LDP; SPEC(1:NS_MAX)=' '
!
	    NSPEC=0
	    DO ML=1,N_LINES
	      Y(ML)=LH(IDEPTH,ML)*(1.0_LDP-SCALE_FAC(ML))
	      I=INDEX(NAME(ML),'(')
	      DO J=1,200
	        IF(SPEC(J) .EQ. ' ')THEN
	          SPEC(J)=NAME(ML)(1:I-1)
	          NSPEC=NSPEC+1
	        END IF
	        IF(NAME(ML)(1:I-1) .EQ. SPEC(J))THEN
	          SUMD(J)=SUMD(J)+Y(ML)
	          EXIT
	        END IF
	      END DO
	    END DO
	    CALL INDEXX(NSPEC,SUMD,INDX,.TRUE.)
!
	    WRITE(6,'(A)')RED_PEN
	    WRITE(6,'(A)')' Output is being written to the file: HEAT_CONTRIBUTIONS'
	    WRITE(6,'(A)')DEF_PEN
	    OPEN(UNIT=27,FILE='HEAT_CONTRIBUTIONS',STATUS='UNKNOWN',ACTION='WRITE',POSITION='APPEND')
	      WRITE(27,'(A)')' '
	      WRITE(27,'(A,2X,I4)')' Contributions to Rad. Equil. Equ at depth',IDEPTH
	      WRITE(27,'(A)')' '
	      WRITE(27,'(X,A,ES14.4,5X,A,ES14.4)')'Scaled sum=',SE_SCL(IDEPTH,N_LINES),'Unscsaled sum',SE_NOSCL(IDEPTH,N_LINES)
	      WRITE(27,'(A)')' '
	      DO I=1,NSPEC
	        J=INDX(I)
	        WRITE(27,'(A10,ES14.4)')TRIM(SPEC(J)),SUMD(J)
	      END DO
!
	      DO_SORT=.TRUE.
	      CALL GEN_IN(DO_SORT,'Create a list of the 30 lines having the largest influence')
	      IF(DO_SORT)THEN
	        T1=SUM(SUMD)
	        Y(1:N_LINES)=ABS(Y(1:N_LINES))/ABS(T1)
	        CALL INDEXX(N_LINES,Y,INDX,.TRUE.)
!
	        WRITE(27,'(A)')' '
	        WRITE(27,'(A,2X,I4)')' Lines with cooling (+ve) / heating (-ve) contributions'
	        WRITE(27,'(A)')' '
	        WRITE(27,'(1X,A,T71,A,2X,A,6X,A,4X,A)')'Transition','Lam(A)','Frac. Contr.','Cur. Sum','Scale Fac.'
	        WRITE(27,'(A)')' '
	        WRITE(27,'(A)')'Total cl tr '
	        WRITE(27,'(A)')' '
	        T1=ABS(SUM(SUMD))
	        CUR_SUM=0.0_LDP
	        DO I=N_LINES,N_LINES-29,-1
	          J=INDX(I)
	          T2=LH(IDEPTH,J)*(1.0_LDP-SCALE_FAC(J))/T1
	          CUR_SUM=CUR_SUM+T2
	          WRITE(27,'(1X,A60,ES16.6,2ES14.3,ES14.5)')NAME(J),LAM(J), LH(IDEPTH,J)*(1.0D0-SCALE_FAC(J))/T1,CUR_SUM,SCALE_FAC(J)
	        END DO
	      END IF
	    CLOSE(UNIT=27)
!
	  ELSE IF(UC(PLT_OPT) .EQ. 'SF')THEN
	    CALL GEN_IN(K,'Depth for plotting')
	    Y(1:N_LINES)=SCALE_FAC(1:N_LINES)-1.0_LDP
	    CALL DP_CURVE(N_LINES,XV,Y)
	    YLABEL='SF-1.0'
!
	  ELSE IF(UC(PLT_OPT) .EQ. 'IDI')THEN
	    CALL GEN_IN(I,'Line index on plot')
	    WRITE(6,'(A)')TRANS_INFO(I)
!	
	  ELSE IF(UC(PLT_OPT) .EQ. 'FV')THEN
	    DO I=1,ND
	      XV(I)=I
	    END DO
	    Y(1:ND)=SE_SCL(1:ND,N_LINES)
	    CALL DP_CURVE(ND,XV,Y)
	    Y(1:ND)=SE_NOSCL(1:ND,N_LINES)
	    CALL DP_CURVE(ND,XV,Y)
	    CALL GRAMON_PGPLOT('Depth Index','STEQ(sc,scl)',' ',' ')
!
	  ELSE IF(UC(PLT_OPT) .EQ. 'WR')THEN
	    K=ND
	    CALL GEN_IN(K,'Depth for plotting')
	    DO I=1,N_LINES
	      WRITE(100,'(I7,4ES14.4,3X,A)')I,NU(I),SE_NOSCL(K,I)-SE_SCL(K,I),
	1                SE_NOSCL(K,I),SE_SCL(K,I),TRIM(NAME(I))
	    END DO
!
	  ELSE IF(UC(PLT_OPT(1:1)) .EQ. 'H')THEN
!
	    WRITE(6,*)' '
	    WRITE(6,*)' XN:    Set X-axis to line index'
	    WRITE(6,*)' XNU:   Set X-axis to frequency'
	    WRITE(6,*)' XLAM:  Set X-axis to wavelength'
	    WRITE(6,*)' '
	    WRITE(6,*)' IDI:   Get transition information for line I'
	    WRITE(6,*)' MAIN:  List most important lines contributing to heatinge error at given depth'
	    WRITE(6,*)'           Also lists contributions by ioization stage'
	    WRITE(6,*)' '
	    WRITE(6,*)' LH:    Plot LH  at given depth'
	    WRITE(6,*)' SS:    Plot STEQ (scaling) at given depth'
	    WRITE(6,*)' SN:    Plot STEQ (no scaling) at a given depth'
	    WRITE(6,*)' FV:    Plot final values (SCL, NO SCL) as a function of depth'
	    WRITE(6,*)' WR:    Write SS & SN data at a single depth to file'
	    WRITE(6,*)' DIFF   Plot contribution by individual lines to the difference between SS and SE.'
	    WRITE(6,*)' E(X):  Exit routine'
!
	  ELSE IF(UC(PLT_OPT) .EQ. 'P')THEN
	    CALL GRAMON_PGPLOT(XLABEL,YLABEL,' ',' ')
!
	  ELSE IF(UC(PLT_OPT) .EQ. 'EX' .OR. UC(PLT_OPT) .EQ. 'E')THEN
	    STOP
	  END IF
	END DO
!
200	WRITE(6,*)STRING
	STOP
!
	END
