
! General purpose line plotting routiine.
!
	SUBROUTINE GRAMON_PGPLOT(XLAB,YLAB,TITL,PASSED_OPT)
	USE SET_KIND_MODULE
	USE NEW_GEN_IN_INTERFACE
	USE MOD_CURVE_DATA
	USE MOD_COLOR_PEN_DEF
	USE LINE_ID_MOD
	IMPLICIT NONE
!
! Altered:  24-Sep-2023 : Added abillity to disallow data initialiazation using PASSED_OPT (05-Sep-2023).
! Altered:  26-Jul-2023 : Added (finalized?) CBAL and FBAL options.
! Altered:  22-Mar-2022 : List of changes compareed with earlier GITHUB version.
!                           SP -- step plot option installed.
!                           Better handling of line IDs
!                           Option installed to do Gauss smoothing.
! Altered:  09-Aug-2022 : Improvements to GF option
!                           (options to draw GF not yet working).
!                           Most old options now part of CGF.
! Altered:  30-Jun-2022 : Changed to allow transition name to be added to EW measurements.
! Altered:  16-Jul-2021 : Added error to link error vector to another vector.
!                           Clarifcation notes added to DC option.
!                           (based on GANNET version - 8-Jul-2021)
! Altered:  22-Nov-2020 : Now call MONBORD_V4 (updated from osiris)
!                       :   Cleaned title label writing
!                       :   Added REG option (removed UG from XAR -- inclued with REG)
!                       :   Added option to add noise (ADDN)
!                       :   Added color blind pens (udated from osiris).
!                       :   Added options to read abcissca (RDXL) and ordnate (RDYL) values (udated from osiris).
! Altered:  15-Jan-2020 : Added NMS option to help create automatic plots.
! Altered:  17-Aug-2019 : Altered to allow dashed lines
! Altered:  10-Jul-2019 : Draw errors first so curves drawn on top.
! Altered:  28-Feb-2019 : Now scale error bars for simple YAR options.
! Altered:  01-Mar-2016 : Added XN option to XAR option. This allows Y to be plotted against
!                           the index I. By default, a new plot is created [24-Feb-2016].
! Altered:  29-Jun-2015 : Changed RID o check ABS(TAU) which for pp model can -ve.
! Altered:  22-Apr-2015 : Added FILL option to fill the space between two curves that create a polygon.
!                           ANS changed to length 4 (from 3)
! Altered:  17-Feb-2015 : Can now have multi-colored titles.
! Altered:  22-Jan-2015 : Bug fix. SC option for scrolling changed to SCR.
!                         SC is reserved for entering strings by cursor
! Altered:  14-Jan-2014 : Revised LG option
! Altered:  22-Nov-2013 : Added LG option for curve type. This plots the log of the absolute
!                           value of the data but indicates where the data is -ve.
! Altered:  04-Sep-2013 : Increased MAXPEN (=MAX_PLOTS). Minor cleaning.
! Altered:  31-Aug-2013 : Added long-plot option.
! Altered:  26-Nov-2011 : Curves cycle over pen-colors 2 to 13.
!                         Dashed curve for plots > 13
!                         Marker style ignored when MARK is off (use I for invisible curve).
! Altered:  08-Apr-2006 : In NM option, normalization now correctly handles unequeally
!                            spaced data.
! Altered:  23-Aug-2005 : Grey pen option installed.
! Altered:  26-Jan-2005 : Change default margins to give more room on borders.
! Altered:  30-May-2003 : Improvements to YAR and XAR options.
! Altered:  11-Feb-2002 : Bug fixed with EW option for reversed data.
! Altered:  04-Apr-2001 : Include option for differen PLT_ST filename.
!                         Code now ouputs name of stored plots when bad plot id entered.
! Altered:  15-Jul-2000 : MOD_CURVE_DATA installed.
!                         Dynamically allocated DAta arrays to allow arbitrary
!                         size plots.
! Altered:  31-May-2000 : Now can handel smaller X and Y ranges, because
!                          exponential format has been included in MON_NUM.
! Altered:  18-Nov-1999 : Maximum number of strings increased to 100.
!                           STR_COL now written to SAVE file.
!                           Only strings inside box are output (unless
!                           PG_LOC is negative.)
! Altered:  05-Mar-1999 : Multiple titles installed. MONBORD_V3 now utilized.
! Altered:  09-Mar-1999 : EXPCHAR_SCALE etc added so that TEXT and plot
!                            symbols have same size relative to the plot,
!                            irrespective of the plotting surface.
! Altered:  30-Jul-1997 : 'W'option installed to allow lines of different
!                            weight. SOme additional work may be needed.
! Altered:  07-Jul-1997 : Revised calls to NEW_GEN_IN_MULT_? installed.
! Altered:  06-May-1997 : Corrections made so that not all curves need
!                            markers.
! Altered:  04-Mar-1997 : Bug fixed in output file name procdure due to
!                            directory conflicts.
! Finalized 07-Mar-1997 : PGPLOT version.
!                            Based on GRAMON (Mongo)
!
        INTEGER, PARAMETER :: MAXSTR=500
	INTEGER, PARAMETER :: MAXVEC=500
	INTEGER, PARAMETER :: MAXPEN=50            !Should be the same as MAX_PLTS
!
	INTEGER NDEC,GET_INDX_SP
	EXTERNAL SPACING,GET_INDX_SP
	REAL*4 SPACING
!
	CHARACTER(LEN=2) TYPE_CURVE(MAX_PLTS)
	REAL(KIND=LDP) VB_BASE(MAX_PLTS)
!
	LOGICAL DONE_NORMALIZATION
	LOGICAL DO_ERROR
	CHARACTER*5 LOG_AXIS
!
	REAL*4 XINC,XNUMST,YNUMST
	REAL*4 YINC
	INTEGER IDX,IXTICK
	INTEGER IDY,IYTICK
	REAL*4 XPAR(2),YPAR(2),XT(2),YT(2)
	REAL*4 XMIN,XMAX,YMIN,YMAX
	REAL*4 XPAR_SAV(2),YPAR_SAV(2)
	REAL(KIND=LDP) YMIN_SAV,YMAX_SAV
!
	REAL*4 YPAR_R_AX(2)
	REAL*4 YNUMST_R_AX
	REAL*4 YINC_R_AX
	LOGICAL NORMAL_R_Y_AXIS
	LOGICAL DO_BORDER
	INTEGER IDY_R_AX,IYTICK_R_AX
	CHARACTER*1 WHICH_Y_AX(MAX_PLTS)
	CHARACTER*80 YLABEL_R_AX
	CHARACTER(LEN=80) XLAB_FILE,YLAB_FILE
!
	INTEGER IFILL_PLT1(10)
	INTEGER IFILL_PLT2(10)
	INTEGER IFILL_CLR(10)
	INTEGER NFILL
	LOGICAL FILL
!
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
        LOGICAL, PARAMETER :: L_FALSE=.FALSE.
	CHARACTER(LEN=80)  XLABEL,YLABEL
	CHARACTER(LEN=*)   XLAB,YLAB,TITL,PASSED_OPT
        CHARACTER(LEN=200) FILNAME
	CHARACTER(LEN=80)  WK_STR,OPTION
	CHARACTER(LEN=200) TMP_STR
	CHARACTER(LEN=6)   TO_TEK
	CHARACTER(LEN=2)   TO_VT
	CHARACTER(LEN=3)   ADVANCE_OPT
	CHARACTER(LEN=10)  LAM_OPTION
	CHARACTER(LEN=80)  PLT_ID,RD_PLT_ID,PLT_ID_SAV
	CHARACTER(LEN=3)   XUNIT
!
	CHARACTER(LEN=80), SAVE :: ID_FILNAME=' '
        CHARACTER(LEN=80), SAVE :: EW_FILNAME=' '
        CHARACTER(LEN=80), SAVE :: EW_LINE_ID=' '
	CHARACTER(LEN=80), SAVE :: PLT_ST_FILENAME
	CHARACTER(LEN=80), SAVE :: TITLE_FILENAME
	CHARACTER(LEN=80), SAVE :: BALMER_INPUT_FILE
!
! Vector arrays
!
	REAL*4 LINEXST(MAXVEC),LINEYST(MAXVEC)
	REAL*4 LINEXEND(MAXVEC),LINEYEND(MAXVEC)
	LOGICAL VEC,FLAGLINE(MAXVEC),INIT
	LOGICAL QUERYFLAG
	INTEGER VECPEN(MAXVEC)
	INTEGER VEC_LINE_STYLE(MAXVEC)
!
	LOGICAL AIR_WAVELENGTHS
	LOGICAL DRAW_GAUSS_HARD
	EXTERNAL LAM_AIR
	REAL(KIND=LDP) LAM_AIR
	REAL(KIND=LDP) DP_T1
!
! String arrays (not labels or titles)
!
	REAL*4 ORIENTATION(MAXSTR),XSTR(MAXSTR),YSTR(MAXSTR)
	REAL*4 STR_EXP(MAXSTR),LOC_PG(MAXSTR)
	REAL*4 XSTRPOS(MAXSTR),YSTRPOS(MAXSTR)
	CHARACTER*80 STRING(MAXSTR)
	CHARACTER(LEN=500) OUTPUT_STRING
	INTEGER ISTR
	INTEGER LOC(MAXSTR)
	INTEGER STR_COL(MAXVEC)
	LOGICAL FLAGSTR(MAXSTR),STR
!
	REAL*4 DXST,DXEND,DYST,DYEND
!
! Pen type, Device ID's
!
! LINE_STYLE   : Type of curve (solid, dsahed etc drawn through data)
! MARKER_STYLE : Type of marker drawn at each datum. If negative, no
!                   curve is drawn between the markers.
!
        INTEGER LINE_STYLE(MAX_PLTS)
        INTEGER LINE_WGT(MAX_PLTS)
	INTEGER MARKER_STYLE(MAX_PLTS)
	LOGICAL HARD,TITONRHS,FIRST,FSTOPEN,DASH,MARK
	LOGICAL FILE_IS_OPEN
	LOGICAL INITIALIZE_ON_EXIT
	LOGICAL RESET_CURVE_LAB
	LOGICAL ADD_COMMA
	LOGICAL, SAVE :: REVERSE_PLOTTING_ORDER=.FALSE.
	LOGICAL REVERSE
!
! E, cursor, and continuum parameters.
!
	REAL*4, ALLOCATABLE :: CONT(:)
	REAL*4 EW,CENTROID
	REAL*4 XCUR(2),YCUR(2),SLOPE
	INTEGER PLOT_ID,CURSERR
	INTEGER L_CHAN(2)
	INTEGER, SAVE :: IP_CONT=0
	INTEGER, SAVE :: LU_EW=30
	INTEGER, SAVE :: LU_NORM=31
	LOGICAL, SAVE :: FIRST_EW=.TRUE.
	CHARACTER*1 CURSVAL
	LOGICAL CONTINUUM_DEFINED
	CHARACTER(LEN=10) DC_CURVE_OPTION
	CHARACTER(LEN=80) DC_INPUT_OPTION
!
! Functions
!
	LOGICAL END_CURS
	INTEGER PGBEG, PGCURS
	REAL(KIND=LDP) SPEED_OF_LIGHT
	REAL(KIND=LDP) LAM_VAC
	REAL(KIND=LDP) LAM_ST,LAM_END
!
	REAL(KIND=LDP) DP_CUT_ACC
	REAL(KIND=LDP) SIG_GAU_KMS
	REAL(KIND=LDP) FRAC_SIG_GAU
	INTEGER NPTS_PER_SIGMA
	REAL(KIND=LDP), ALLOCATABLE :: WRK1(:)
	REAL(KIND=LDP), ALLOCATABLE :: WRK2(:)
!
	REAL*4,  PARAMETER :: RONE=1.0
	INTEGER, PARAMETER :: IZERO=0
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: ITWO=2
	INTEGER, PARAMETER :: IFOUR=4
	INTEGER, PARAMETER :: IFIVE=5
!
! Parameters for indicating the size of the plot.
!
	REAL*4 XCM,ASR,TEMPASR,DASR
!
! CENTRAL_LAM must be REAL(KIND=LDP) as LAM_VAC is of KIND LDP function.
!
	REAL(KIND=LDP) CENTRAL_LAM
	REAL(KIND=LDP) OLD_CENTRAL_LAM
	REAL(KIND=LDP) C_KMS
	REAL(KIND=LDP) C_VAL
	LOGICAL AIR_LAM
	CHARACTER(LEN=5), SAVE :: VEL_UNIT='km/s'
	CHARACTER(LEN=10) GF_OPTION
!
! For XAR and YAR arithmetic options.
!
	CHARACTER(LEN=4) XAR_OPERATION,YAR_OPERATION,VAR_OPERATION
	CHARACTER(LEN=5) REG_OPT
	REAL*4 YAR_VAL,XAR_VAL
	INTEGER XAR_PLT,YAR_PLT
	INTEGER VAR_PLT1,VAR_PLT2,VAR_PLT3
!
! Miscellaneous
!
	CHARACTER*4 ANS
	INTEGER LENGTH
	INTEGER LAST_DP
	REAL*4 V1,IT
	REAL*4 EXPCHAR,TICK_FAC,EXPMARK
	REAL*4 EXPCHAR_SCALE,TICK_FAC_SCALE,EXPMARK_SCALE
	REAL*4 XCHAR_SIZE,YCHAR_SIZE
	REAL*4 MARGINX(2),MARGINY(2)
	REAL*4 MEAN,SIGMA
	REAL*4 T1,T2,T3,T4
	REAL*4 SCALE_FAC(MAX_PLTS)
 	REAL*4 XVAL,YVAL
	REAL*4 XVAL_SAV,YVAL_SAV
	REAL*4, ALLOCATABLE :: TA(:)
	LOGICAL TMP_LOG
	LOGICAL KEEP_YAXIS_LIMITS
	LOGICAL, SAVE :: WRITE_COMMENT
!
	INTEGER, SAVE :: PEN_OFFSET
	INTEGER BEG
	INTEGER Q	!Used for pen color
!
	INTEGER, PARAMETER :: T_IN=5
	INTEGER, PARAMETER :: T_OUT=6
!
! Loop variables
	INTEGER I,J,K,L,CNT,IP,OP
	INTEGER IP_ST,IP_END,IP_INC
!
! Color variables
	REAL*4 RED(0:15),BLUE(0:15),GREEN(0:15)
	INTEGER PEN_COL(0:MAXPEN)
	REAL*4 RTEMP,BTEMP,GTEMP
!
	INTEGER NP_OUT
	INTEGER NX_MAX
	INTEGER IPLT(MAX_PLTS)
	INTEGER IPST(MAX_PLTS)
!
! Variables to read in plot in column format
!
	REAL*4 TEMP_VAR(20)
	INTEGER COLUMN(2)
!
! Printer variables.
!
	CHARACTER(LEN=80), SAVE :: PRINTER='pgplot_1.ps/cps'
	CHARACTER(LEN=80), SAVE :: HARD_FILE
	CHARACTER(LEN=20), SAVE :: HARD_TYPE
	CHARACTER(LEN=80) CUR_HARD_FILE
	INTEGER, SAVE :: HARD_CNT
	LOGICAL, SAVE :: FIRST_HARD=.TRUE.
!
	INTEGER PGOPEN,PLT_LINE_WGT
	INTEGER ID
	SAVE ID
	REAL*4 TOTXMM,TOTYMM,SCALEFACY,SCALEFAC
	REAL*4 PRINTX1,PRINTX2,PRINTY1,PRINTY2
!
	LOGICAL SMOOTH_PLOT(MAX_PLTS)
	LOGICAL BOX_FILTER
	LOGICAL LONG_PLOT
	REAL*4 LENGTH_OF_HC_PLOT
	REAL*4 LP_ASR
!
! Variables for options 'WP' (i.e. write plot) and 'RP' (i.e. read plot).
!
	INTEGER IST,IEND,N_REC_SIZE
!
	SAVE XCM,ASR
	SAVE EXPCHAR_SCALE,EXPMARK_SCALE,TICK_FAC_SCALE
	SAVE RED,BLUE,GREEN
	SAVE FSTOPEN,PEN_COL,DASH
	SAVE MARGINX,MARGINY
	SAVE PLT_LINE_WGT
	SAVE LINE_WGT
	DATA FSTOPEN,DASH/.TRUE.,.FALSE./
	DATA PLT_LINE_WGT/1/
	DATA HARD_CNT/1/
	N_REC_SIZE=1000
!
	IF(NPLTS .EQ. 0)THEN
	  WRITE(T_OUT,*)'Error - No calls made to curve'
	  RETURN
	END IF
	N_LINE_IDS=0
	ID_SCL=1.05D0
	ID_VEC_BEG=1.01D0
	ID_VEC_END=1.04D0
	ID_EXPCHAR=1.0D0
	TAU_CUT=0.1D0
	N_OMIT_ID=0
	N_INC_ID=0
	DO_BORDER=.TRUE.
	DRAW_GAUSS_HARD=.FALSE.
	C_KMS=1.0D-05*SPEED_OF_LIGHT()
	LONG_PLOT=.FALSE.
	LENGTH_OF_HC_PLOT=200.0D0       !cm
	VB_BASE=-1000
	EW_SCALE_FAC=0.0D0
	FILL=.FALSE.
	SIG_GAU_KMS=10.0D0
	FRAC_SIG_GAU=0.25D0
	NPTS_PER_SIGMA=4
	DP_CUT_ACC=0.01D0
	DONE_NORMALIZATION=.FALSE.
	XLAB_FILE=' '; YLAB_FILE=' '
	LINE_CUT_PARAM=0.03
	USE_DEF_OFFSET=.TRUE.
	XUNIT='Ang'
	NO_DEC_DIGITS=1
!
	IF(NPLTS .GT. MAXPEN)THEN
	  WRITE(T_OUT,*)'Error n GRAMON_PLOT -- not enough pen loctions'
	  WRITE(T_OUT,*)'MAXPEN should be set to the same value as MAX_PLTS'
	  RETURN
	END IF
!
! Define character strings for switching between VT and TEK modes.
!
	TO_TEK(1:1)=CHAR(27)
	TO_TEK(2:6)='[?38h'
	TO_VT(1:1)=CHAR(27)
	TO_VT(2:2)=CHAR(3)
!
! Arithmetic options
!
	XAR_OPERATION='*'; XAR_VAL=1.0; XAR_PLT=1
	YAR_OPERATION='*'; YAR_VAL=1.0; YAR_PLT=1
	VAR_OPERATION='*'; VAR_PLT1=1;  VAR_PLT2=1
!
! Set the Aspect ratio of the plot to the default value. No longer ask
! for a new value, since can be set with N option.
!
	IF(FSTOPEN)THEN
	  ASR=0.0			!Leave as is
	  EXPCHAR_SCALE=1.5   		!Was 1.01 to prevent crashes.
	  EXPMARK_SCALE=1.5   		!Symbol size on plots.
	  TICK_FAC_SCALE=1.5D0
	  PLT_LINE_WGT=1
          PLT_ST_FILENAME=' '
          TITLE_FILENAME='TITLE.SAV'
	  LINE_WGT(:)=1
	  PEN_OFFSET=1
	  IFILL_PLT1=0; IFILL_PLT2=0
	  ID_FILNAME='LINE_ID'
	  EW_FILNAME='EWDATA'
	  EW_LINE_ID='ewdata_fin'
	  EW_CUT=1.0D0
	  WRITE_COMMENT=.FALSE.
	  ID_LINE_PEN=1
	  BALMER_INPUT_FILE='BALMER_LINE_LIMS'
	  TYPE_CURVE='L'
	END IF
	CALL GEN_ASCI_OPEN(LU_NORM,'NORM_FACTORS','UNKNOWN','APPEND',' ',IZERO,IOS)
!
	LOG_AXIS=' '
	XLABEL=XLAB
	YLABEL=YLAB
	YLABEL_R_AX=' '
	IF(N_HEADER.EQ. 0)THEN
          IF(TITL .NE. ' ')THEN
	    TITLE(1)=TITL
	    N_HEADER=1
	  END IF
	END IF
	N_TIT_SET=N_HEADER
	TITLE(N_HEADER+1:)=' '
	STR=.FALSE.
	VEC=.FALSE.
	NORMAL_R_Y_AXIS=.TRUE.
	RESET_CURVE_LAB=.TRUE.
	ADD_COMMA=.TRUE.
	DO I=1,MAXSTR
	  FLAGSTR(I)=.FALSE.
	  FLAGLINE(I)=.FALSE.
	  STR_EXP(I)=1.0	   	!Default (changed using SE option only)
	  STR_COL(I)=1			!Default (changed using SE option only)
	  LINEXST(I)=1; LINEYST(I)=1
	  LINEXEND(I)=1; LINEYEND(I)=1
	END DO
	DO_ERROR=.TRUE.
	OPTION=PASSED_OPT
	CALL SET_CASE_UP(OPTION,1,0)
	CENTRAL_LAM=1548.20			!CIV
	OLD_CENTRAL_LAM=0.0D0
	CONTINUUM_DEFINED=.FALSE.
	AIR_WAVELENGTHS=.FALSE.
	PLOT_ID=1
	KEEP_YAXIS_LIMITS=.FALSE.
!
	YPAR_R_AX(1)=0.0D0 ; YPAR_R_AX(2)=0.0D0
!
! Determine type of graph.
!
	IF(OPTION(1:4) .EQ. 'MARK')THEN
	  MARK=.TRUE.				!Connect plots
	ELSE IF(OPTION(1:4) .EQ. 'HIST')THEN
	  TYPE_CURVE(1:NPLTS)='H'
	  MARK=.FALSE.
	ELSE IF(OPTION(1:3) .EQ. 'ADJ')THEN  !Histogram, but adjacent points
	  TYPE_CURVE(1:NPLTS)='A'
	  MARK=.FALSE.
	ELSE
	  TYPE_CURVE(1:NPLTS)='L'
	  MARK=.FALSE.
	END IF
	IF(OPTION(1:3) .EQ. 'NOI')THEN  !Don't initialize plots on exit
	  INITIALIZE_ON_EXIT=.FALSE.
	ELSE
	  INITIALIZE_ON_EXIT=.TRUE.
	END IF
	WHICH_Y_AX(1:NPLTS)='L'
!
! Define default line representations (initially not dashed)
!
	DO I=1,MAX_PLTS
	  LINE_STYLE(I)=1
	  MARKER_STYLE(I)=MOD(I,5)
	END DO
	DO I=14,MAX_PLTS
	  LINE_STYLE(I)=2
	END DO
	DASH=.FALSE.
!
! Assign a color index to each pen. Keep previus assignments if they have been
! made. We avoid the white pen as it won't show on a hardcopy.
!
	IF (FSTOPEN) THEN
	  PEN_COL(0)=0
	  DO I=1,MAXPEN
	    PEN_COL(I)=MOD(I-1,14)+1
	  END DO
	END IF
!
! Assign a color index to each vecpen
!
	IF (FSTOPEN) THEN
	  DO I=1,10
	    VECPEN(I)=I+5
	  END DO
	  DO I=11,MAXVEC
	    VECPEN(I)=14
	  END DO
	  VEC_LINE_STYLE(:)=1
	END IF
!
! Define default type for "Data markers"
!
	DO I=1,NPLTS
	  MARKER_STYLE(I)=MOD(I,12)+1           !1 is not appropriate
	END DO
!
	ANS='P'
	FIRST=.TRUE.
!
! Get absica and ordinate limits.
!
	CALL GET_GRAMON_MIN_MAX(XMIN,XMAX,YMIN,YMAX,TYPE_CURVE,T_OUT)
!
	XPAR(1)=XMIN
	XPAR(2)=XMAX
	YPAR(1)=YMIN
	YPAR(2)=YMAX
	XPAR_SAV(1)=XPAR(1)		!Indicate value limits evaluated for.
	XPAR_SAV(2)=XPAR(2)
	YMAX_SAV=YMAX; YMIN_SAV=YMIN
!
! Open user set workstation (default is set into the system).
! We only ask for work-station name if first call to GRAMON.
!
	WRITE(T_OUT,4) XMIN,XMAX
4	FORMAT(' Abisca limits  :',1P2E14.4)
	WRITE(T_OUT,11) YMIN,YMAX
11	FORMAT(' Ordinate limits:',1P2E14.4)
	WRITE(T_OUT,'(A)')' '
!
	IF (FSTOPEN) THEN
	  WRITE(T_OUT,*)'NOTE: /XWINDOW supports color'
	  WRITE(T_OUT,*)'    while /XTERM allows auto-switching of the cursor'
	  BEG = PGBEG(IZERO,'?',IONE,IONE)
	  CALL PGQID(ID)
!
! Initialize viewport.
!
	  CALL PGASK(.FALSE.)                      !TURNS OFF PAGE PROMPT
	  CALL PGVSTD
	  CALL PGENV(XPAR(1),XPAR(2),YPAR(1),YPAR(2),IZERO,IZERO)
	  CALL PGQVP(IZERO,MARGINX(1),MARGINX(2),MARGINY(1),MARGINY(2))
	  IF (MARGINX(1) .LT. (.9)) MARGINX(1)=MARGINX(1)+.05
!
! Increase default border size.
!
	  MARGINX(1)=0.15; MARGINX(2)=0.9
          MARGINY(1)=0.15; MARGINY(2)=0.9
!
! Set preferred defaults for pen colors.
!
	  CALL PGSCR(0,.7,.7,.7)     !set color representations
	  CALL PGSCR(1,0.0,0.0,0.0)
	  CALL PGSCR(2,1.0,0.0,0.0)
	  CALL PGSCR(3,0.0,0.0,1.0)
	  CALL PGSCR(4,0.0,0.6,0.3)
	  CALL PGSCR(5,.5,0.0,.7)
	  CALL PGSCR(13,0.0,1.0,1.0)
	  CALL PGSCR(15,.95,.95,.95)
          DO I=0,15                  !Get these + def color representations.
            CALL PGQCR(I,RED(I),GREEN(I),BLUE(I))
          END DO
!
	  CALL DEFINE_MORE_PENS(MAXPEN)
	  DO I=26,30
	    IFILL_CLR(I-25)=I
	  END DO
!
	END IF
	FSTOPEN=.FALSE.
	HARD=.FALSE.
1000	ANS='P'
!
!	WRITE(T_OUT,*)TO_VT
!
	WRITE(T_OUT,*)' H,P,Z(ZN),A(F),L,CC,CP,N, M,W,D,C,B,'//
	1              ' VC,VF,VE, SC,SF,SE, E'
	CALL NEW_GEN_IN(ANS,'ANS')
	L=LEN_TRIM(ANS)
	CALL SET_CASE_UP(ANS,IONE,IZERO)
        IF(ANS .EQ. 'H')THEN
	  WRITE(T_OUT,*)' '
	  WRITE(T_OUT,*)'P   - Plot - default'
          WRITE(T_OUT,*)'E   - EXIT from PLOT package'
	  WRITE(T_OUT,*)'Z   - Hardcopy (ZN=Asks for new hard device)'
          WRITE(T_OUT,*)'LP  - Allow a long hard copy postscript plot to be created (with CPS)'
	  WRITE(T_OUT,*)'CUT - Cut unecessary points from a curve to reduce final plot size'
	  WRITE(T_OUT,*)' '
          WRITE(T_OUT,*)'A    - Define Axis Parameters'
          WRITE(T_OUT,*)'2A   - Define labeling of right-hand axis'
          WRITE(T_OUT,*)'F    - Change default axis parameters'
          WRITE(T_OUT,*)'N    - Define aspect ratio, plot margins, character height etc'
	  WRITE(T_OUT,*)'L    - Modify Axis Labels and Titles'
	  WRITE(T_OUT,*)'EDCL - Edit curve IDs/labels'
	  WRITE(T_OUT,*)' '
 	  WRITE(T_OUT,*)'D   - Switch dashed lines on/off'
 	  WRITE(T_OUT,*)'DE  - Edit dashed lines one by one'
	  WRITE(T_OUT,*)'W   - Change thickness (weights) of curves'
 	  WRITE(T_OUT,*)'WE  - Edit line weights one by one'
	  WRITE(T_OUT,*)'M   - Switch marking data points on/off'
	  WRITE(T_OUT,*)'C   - Indicate how curves are to be connected (L,H,A,E,I,V,B)'
	  WRITE(T_OUT,*)'CC  - Change Color setting'
	  WRITE(T_OUT,*)'CP  - Change Pen (Color Index)'
	  WRITE(T_OUT,*)'RCP - Reset default color pens'
	  WRITE(T_OUT,*)'CBP - Set pens for color blindness'
	  WRITE(T_OUT,*)'GP  - Set default for grey pens'
	  WRITE(T_OUT,*)'SFP - Adds an offset to color pen defs. (may not work with all opts)'
	  READ(T_IN,'(A)')ANS				!can use ANS here.
	  IF(ANS(1:1) .EQ. 'E' .OR. ANS(1:1) .EQ. 'e')GOTO 1000
!
	  WRITE(T_OUT,*)'B    - Switch error bars on/off'
	  WRITE(T_OUT,*)'LNKE - Link a vector as errors to another plot'
	  WRITE(T_OUT,*)'BRD  - Switch border potting on (def) or off'
	  WRITE(T_OUT,*)'FILL - Color in region between 2 curves'
	  WRITE(T_OUT,*)'OFF  - Set offsets when plotting multiple plots'
          WRITE(T_OUT,*)'RPO  - Plot curves in reverse order (switch): does not affect color'
	  WRITE(T_OUT,*)'RID  - Read line ID''s - uses file created by DISPGEN (for O stars)'
	  WRITE(T_OUT,*)'REW  - Read line ID''s - uses EWDATA file created by CMF_FLUX (emission line stars)'
	  WRITE(T_OUT,*)'SID  - Change defaults for writing line ID''s'
	  WRITE(T_OUT,*)' '
!
	  WRITE(T_OUT,*)'LX  - Switch between LINEAR/LOG labeling of X axis'
	  WRITE(T_OUT,*)'LY  - Switch between LINEAR/LOG labeling of Y axis'
	  WRITE(T_OUT,*)'LXY - Switch between LINEAR/LOG baleling of X and Y axes'
	  WRITE(T_OUT,*)' '
	  WRITE(T_OUT,*)'VC  - Define line vectors using cursor'
	  WRITE(T_OUT,*)'VF  - Define line vectors using file input'
	  WRITE(T_OUT,*)'VE  - Online edit of vectors'
          WRITE(T_OUT,*)'SC  - Define strings using cursor'
          WRITE(T_OUT,*)'SF  - Define strings using file input'
          WRITE(T_OUT,*)'SE  - Online edit of strings'
	  READ(T_IN,'(A)')ANS				!can use ANS here.
	  IF(ANS(1:1) .EQ. 'E' .OR. ANS(1:1) .EQ. 'e')GOTO 1000
!
	  WRITE(T_OUT,*)'LAM  - List wavelengths of common lines'
	  WRITE(T_OUT,*)'VEL  - Convert X axis to km/s space'
!
	  WRITE(T_OUT,*)'XAR  - Simple X axis arithmetic'
	  WRITE(T_OUT,*)'YAR  - Simple Y axis arithmetic'
	  WRITE(T_OUT,*)'VAR  - Simple arithmetic on two plots'
	  WRITE(T_OUT,*)'NM   - Scale average to 1 or to another plot'
	  WRITE(T_OUT,*)'NMS  - Similar to NM but for automatic plotting'
	  WRITE(T_OUT,*)'SP   - Step plot accross screen'
	  WRITE(T_OUT,*)' '
	  WRITE(T_OUT,*)'REP  - Simple replacment of data using a cursor'
	  WRITE(T_OUT,*)'FREP - Same as REP except nodes read in from file'
	  WRITE(T_OUT,*)'MODC - Routine for modifying/shifting a simple curve'
	  WRITE(T_OUT,*)'REG  - Regrid plot - UG, dX, R or NINS'
	  WRITE(T_OUT,*)'SIG  - Measure signal to noise of a spectral region'
	  WRITE(T_OUT,*)'ADDN - Add Poisonian noise'
	  WRITE(T_OUT,*)'SM   - Smooth data (Han) -- ignores X-spacing of data'
	  WRITE(T_OUT,*)'BXSM - Smooth data (box filter) -- ignores X-spacing of data'
	  WRITE(T_OUT,*)'GSM  - Gaussian smoothing -- set resolution'
!
	  WRITE(T_OUT,*)' '
	  WRITE(T_OUT,*)'LOC  - Use a cursor to read of (X,Y) coordinates on a plot'
	  WRITE(T_OUT,*)'CONT - Define a continuum using cursors'
	  WRITE(T_OUT,*)'DC   - Define a continuum for EW measurmenst (line segments or monotonc cubic'
	  WRITE(T_OUT,*)'MCN  - Define a continuum with repeate call using cursors'
	  WRITE(T_OUT,*)' '
!
	  WRITE(T_OUT,*)'CEW  - Measure the EW of a single line using cursors and a local continuum'
	  WRITE(T_OUT,*)'FEW  - Measure the EW of a line using locations read in from file set by CEW'
	  WRITE(T_OUT,*)'EW   - Measure the EW of a single line'
	  WRITE(T_OUT,*)'EWG  - Measure the EW of many lines use Gaussian fiting usig data from a file'
	  WRITE(T_OUT,*)'CGF  - Fit a (modfied) gaussian to an absorption or emission line'
	  WRITE(T_OUT,*)'FGF  - Fit multiple spctral regions using params read in fro GAUSS_PARAMS'
	  WRITE(T_OUT,*)'CBAL - Measure Chi^2 of model fit to observations and line EW (mainly for broad lines)'
	  WRITE(T_OUT,*)'FBAL - Use the results of CBAL to model fit to observatuobs and line EW'
!	  WRITE(T_OUT,*)'EGF  - Edit gauss-fit arameters'
!	  WRITE(T_OUT,*)'DG   - Draw gauss-fit.'
!	  WRITE(T_OUT,*)'WGF  - Write gauss-fit parameters to a file'
	  READ(T_IN,'(A)')ANS				!can use ANS here.
	  IF(ANS(1:1) .EQ. 'E' .OR. ANS(1:1) .EQ. 'e')GOTO 1000
!
	  WRITE(T_OUT,*)' '
	  WRITE(T_OUT,*)'RXY  - Read plot from asci file'
	  WRITE(T_OUT,*)'WXY  - Write plot to asci file'
	  WRITE(T_OUT,*)'WOBS - Write plot to asci file in same format as OBSFLUX'
	  WRITE(T_OUT,*)'SXY  - Write section of data to terminal'
	  WRITE(T_OUT,*)'RP   - Read labeled plots from direct accecs file'
	  WRITE(T_OUT,*)'RPF  - Similar to RP but asks for filename'
	  WRITE(T_OUT,*)'WP   - Write labeled plots to direct access file'
	  WRITE(T_OUT,*)'WPF  - Similar to WP but asks for filename'
	  WRITE(T_OUT,*)'WTIT - Write titles to ascii file'
	  WRITE(T_OUT,*)'RTIT - Read titles from ascii file'
	  WRITE(T_OUT,*)'RDXL - Read abscica values from file'
	  WRITE(T_OUT,*)'RDYL - Read ordinate values from file'
	  READ(T_IN,'(A)')ANS				!can use ANS here.
	  IF(ANS(1:1) .EQ. 'E' .OR. ANS(1:1) .EQ. 'e')GOTO 1000
!
	  WRITE(T_OUT,*)'OLF - Open LOG file'
	  WRITE(T_OUT,*)'OIF - Open input file'
	  WRITE(T_OUT,*)'CLF - Close LOG file'
	  WRITE(T_OUT,*)'CIF - Close inpit file'
	  WRITE(T_OUT,*)'CL  - Clear Graphics Screen'
!
	  WRITE(T_OUT,*)'NOI - Leave data intact on exit (switch)'
	  WRITE(T_OUT,*)'H   - Help'
	  WRITE(T_OUT,*)' '
          GOTO 1000
!
! Exit from Ploting package, saving STRING and VECTOR information.
! If the NOI option has been issued, the plots will still be
! retained. ROutine checks that CONTINUUM defition saved is it has been set.
!
	ELSE IF(ANS .EQ. 'E')THEN
	  IF(INITIALIZE_ON_EXIT)THEN
	    N_TIT_SET=0
	    IF(IP_CONT .NE. 0)THEN
	      WRITE(6,*)'Have you saved you defined continuum values'
	      QUERYFLAG=.FALSE.
	      CALL NEW_GEN_IN(QUERYFLAG,'Return to OPTION input to use WXY option?')
	      IF(QUERYFLAG)GOTO 1000
	    ELSE
	      IP_CONT=0
	    END IF
	    DO IP=1,NPLTS
	      IF(ALLOCATED(CD(IP)%XVEC))DEALLOCATE(CD(IP)%XVEC)
	      IF(ALLOCATED(CD(IP)%DATA))DEALLOCATE(CD(IP)%DATA)
	      IF(ALLOCATED(CD(IP)%EMAX))DEALLOCATE(CD(IP)%EMAX)
	      IF(ALLOCATED(CD(IP)%EMIN))DEALLOCATE(CD(IP)%EMIN)
	      ERR(IP)=.FALSE.
	      NPTS(IP)=0
	    END DO
	    NPLTS=0
	    N_TIT_SET=0
	    TITLE=' '
	  END IF
	  CLOSE(LU_NORM)
          IF(STR .OR. VEC)THEN
	    FILNAME='PGOUT'
	    CALL NEW_GEN_IN(FILNAME,'Filname for STRING/VEC storage (no ext)')
	    IF(FILNAME .EQ. ' ')FILNAME='JNK'	!Force output to jnk file as
	    L=LEN_TRIM(FILNAME)			!precaution.
	    IF(FILNAME .NE. ' ')THEN
	      IF(STR)THEN
	        FILNAME=FILNAME(:L)//'.STR'
                L=L+4
                OPEN(UNIT=33,FILE=FILNAME(:L),STATUS='UNKNOWN')
	          DO I=1,MAXSTR
                    IF(FLAGSTR(I))THEN
                      WRITE(33,17)LOC(I),XSTR(I),YSTR(I),ORIENTATION(I),
	1                   STR_EXP(I),STR_COL(I),TRIM(STRING(I))
17	              FORMAT(1X,I1,', ',4(ES12.4,','),I3,',',1X,1H',A,1H')
	           END IF
	         END DO
                CLOSE(UNIT=33)
	      END IF
	      IF(VEC)THEN
	        IF(STR)THEN
	          FILNAME(L-3:L)='.VEC'
	        ELSE
	          FILNAME=FILNAME(:L)//'.VEC'
                  L=L+4
	        END IF
                OPEN(UNIT=33,FILE=FILNAME(:L),STATUS='UNKNOWN')
	          DO I=1,MAXVEC
                    IF(FLAGLINE(I))THEN
	              WRITE(33,18)LINEXST(I),LINEYST(I),
	1                         LINEXEND(I),LINEYEND(I),VECPEN(I),VEC_LINE_STYLE(I)
18	              FORMAT(1X,1P,4E18.8,3X,2I3)
	            END IF
	          END DO
                CLOSE(UNIT=33)
	      END IF		!VEC end
            END IF		!Filename check
          END IF		!STR or VEC present
	  RETURN
	ELSE IF(ANS .EQ. 'OLF')THEN
	  FILNAME='LOG_FILE'
	  CALL NEW_GEN_IN(FILNAME,'Log file=')
	  CALL NEW_GEN_IN_OPTS('OPEN_LOG_FILE',FILNAME,20)
          GOTO 1000
	ELSE IF(ANS .EQ. 'CLF')THEN
	  CALL NEW_GEN_IN_OPTS('CLOSE_LOG_FILE','LOG_FILE',20)
          GOTO 1000
	ELSE IF(ANS .EQ. 'OIF')THEN
	  CALL NEW_GEN_IN(FILNAME,'Command input file')
	  CALL NEW_GEN_IN_OPTS('OPEN_FILE_INPUT',FILNAME,21)
          GOTO 1000
	ELSE IF(ANS .EQ. 'CIF')THEN
	  CALL NEW_GEN_IN_OPTS('CLOSE_FILE_INPUT',FILNAME,21)
          GOTO 1000
	END IF
C
	IF(FIRST)THEN
C Define xinc
	  XINC=SPACING(XPAR(1),XPAR(2))
	  YINC=SPACING(YPAR(1),YPAR(2))
	  IXTICK=2
	  IYTICK=2
C Compute Xmin,Xmax,Ymin,Ymax if not input
	  V1=1000
	  XPAR(1)=ABS(XINC)*(AINT(XMIN/ABS(XINC)+V1)-V1)
	  XPAR(2)=ABS(XINC)*(AINT(XMAX/ABS(XINC)-V1)+V1)
	  YPAR(1)=ABS(YINC)*(AINT(YMIN/ABS(YINC)+V1)-V1)
	  YPAR(2)=ABS(YINC)*(AINT(YMAX/ABS(YINC)-V1)+V1)
C
	  XNUMST=XPAR(1)
	  YNUMST=YPAR(1)
C
C Determine number of decimals
	  IDX=NDEC(XINC)
	  IT=NDEC(XPAR(1))
	  IF(IT .GT. IDX)IDX=IT
	  IDY=NDEC(YINC)
	  IT=NDEC(YPAR(1))
	  IF(IT .GT. IDY)IDY=IT
	  FIRST=.FALSE.
	END IF
C
	IF(ANS .EQ. 'A' .OR. ANS .EQ. 'F')THEN
!
! Get absica and ordinate limits.
!
	  CALL GET_GRAMON_MIN_MAX(XMIN,XMAX,YMIN,YMAX,TYPE_CURVE,T_OUT)
!
	  WRITE(T_OUT,4) XMIN,XMAX
	  CALL NEW_GEN_IN(XPAR,I,ITWO,'XST,XEND')
C
C Look for ordinate limits in X range. This will now generate a plot with
C roughly the correct scaling even though the oridinates values may be vastly
C different outside the plot window. We only use those plots that will be
C displayed.
C
	    WRITE(6,*)'Getting Y limits'
	    YMIN=1.0E+32
	    YMAX=-1.0E+32
	    DO IP=1,NPLTS
	      IF(TYPE_CURVE(IP) .EQ. 'LG')THEN
	        T1=1.0E+32;T2=-1.0E+32
	        DO I=1,NPTS(IP)
	          IF( (CD(IP)%XVEC(I) .GE. XPAR(1) .AND. CD(IP)%XVEC(I) .LE. XPAR(2)) .OR.
	1           (CD(IP)%XVEC(I) .GE. XPAR(2) .AND. CD(IP)%XVEC(I) .LE. XPAR(1)) )THEN
	             T2=MAX(T2,ABS(CD(IP)%DATA(I)))
	             T1=MIN(T1,ABS(CD(IP)%DATA(I)))
	          END IF
	        END DO
	        YMAX=MAX(LOG10(T2),YMAX)
	        IF(T1 .EQ. 0.0D0)THEN
	          YMIN=-38
	        ELSE
	          YMIN=MIN(YMIN,LOG10(T1))
	        END IF
	      ELSE IF(TYPE_CURVE(IP) .NE. 'I')THEN
	        DO I=1,NPTS(IP)
	          IF( (CD(IP)%XVEC(I) .GE. XPAR(1) .AND. CD(IP)%XVEC(I) .LE. XPAR(2)) .OR.
	1           (CD(IP)%XVEC(I) .GE. XPAR(2) .AND. CD(IP)%XVEC(I) .LE. XPAR(1)) )THEN
	             YMAX=MAX(YMAX,CD(IP)%DATA(I))
	             YMIN=MIN(YMIN,CD(IP)%DATA(I))
	          END IF
	        END DO
	      END IF
	    END DO
	    IF(YMAX .EQ. 0 .AND. YMIN .EQ. 0)THEN
	      WRITE(T_OUT,*)'Poor YMIN & YMAX values'
	      WRITE(T_OUT,*)'YMIN=',YMIN,'  YMAX=',YMAX
	      YMIN=0.0
	      YMAX=1.0
	    END IF
!
	  IF(KEEP_YAXIS_LIMITS)THEN
	  ELSE
	    WRITE(6,*)YMIN,YMAX
	    IF(XPAR(2) .NE. XPAR_SAV(2) .OR. XPAR(1) .NE. XPAR_SAV(1) .OR.
	1      YMAX .NE. YMAX_SAV .OR. YMIN .NE. YMIN_SAV)THEN
	       YINC=SPACING(YMIN,YMAX)
	       IYTICK=2
	      V1=1000
	      YPAR(1)=ABS(YINC)*(AINT(YMIN/ABS(YINC)+V1)-V1)
	      YPAR(2)=ABS(YINC)*(AINT(YMAX/ABS(YINC)-V1)+V1)
	      XPAR_SAV(1)=XPAR(1)		!Indicate value limits evaluated for.
	      X PAR_SAV(2)=XPAR(2)
	      YMIN_SAV=YMIN; YMAX_SAV=YMAX
	    END IF
	  END IF
C
	  WRITE(T_OUT,11)YMIN,YMAX
	  IF(KEEP_YAXIS_LIMITS)THEN
	    YPAR_SAV=YPAR
	    CALL NEW_GEN_IN(YPAR,I,ITWO,'YST,YEND')
	    IF(YPAR(2) .NE. YPAR_SAV(2) .OR. YPAR(1) .NE. YPAR_SAV(1))KEEP_YAXIS_LIMITS=.FALSE.
	  ELSE
	    CALL NEW_GEN_IN(YPAR,I,ITWO,'YST,YEND')
	  END IF
C
C Check that limits inserted are not absurd.
C
	  IF(ABS(YPAR(2)-YPAR(1))/MAX(ABS(YPAR(2)),ABS(YPAR(1))) .LT. 1.0E-08)THEN
	    YPAR(1)=0
	    YPAR(2)=1.00
	    WRITE(T_OUT,*)'Invalid Y limits - setting default values'
	  END IF
	  IF(ABS(XPAR(2)-XPAR(1))/MAX(ABS(XPAR(2)),ABS(XPAR(1)))  .LT. 1.0E-08)THEN
	    XPAR(1)=0
	    XPAR(2)=1.00
	    WRITE(T_OUT,*)'Invalid X limits - setting default values'
	  END IF
C
	  XNUMST=XPAR(1)
	  YNUMST=YPAR(1)
C
C Set defaults for XINC and YINC. IXTICK-1 is the number of tick
C marks between numbered ticks.
C
	  XINC=SPACING(XPAR(1),XPAR(2))
	  IXTICK=2
	  YINC=SPACING(YPAR(1),YPAR(2))
	  IYTICK=2
	  WRITE(T_OUT,*)'Indicate spacing between numberd tickmarks'
	  CALL NEW_GEN_IN(XINC,'For X tickmarks')
	  CALL NEW_GEN_IN(YINC,'For Y tickmarks')
C
C Determine number of decimals
C
	  IDX=NDEC(XINC)
	  IT=NDEC(XPAR(1))
	  IF(IT .GT. IDX)IDX=IT
	  IDY=NDEC(YINC)
	  IT=NDEC(YPAR(1))
	  IF(IT .GT. IDY)IDY=IT
C
	  IF(ANS .EQ. 'A')GOTO 1000
C
C This section is for additional axis fiddling. Can change number
C of digits after decimal point, and numbers at which axis labeling
C begins. These parameters are restored to their original values
C everytime 'A' option is called.
C
	  CALL NEW_GEN_IN(IXTICK,'X minor divisions')
	  IXTICK=MAX(IXTICK,1)
	  CALL NEW_GEN_IN(IYTICK,'Y minor divisions')
	  IYTICK=MAX(IYTICK,1)
	  CALL NEW_GEN_IN(IDX,'# of X Dec digits')
	  CALL NEW_GEN_IN(IDY,'# of Y Dec digits')
	  CALL NEW_GEN_IN(XNUMST,'Number beginning X Axis')
	  CALL NEW_GEN_IN(YNUMST,'Number beginning Y Axis')
C
	  FIRST=.FALSE.
	  GOTO 1000
C
	ELSE IF(ANS .EQ. '2A')THEN
C
C Option to allow a left hand and a right hand axis. The default is a single
C left hand axis. As this is primary for pretty plots, the options are less
C powerful than for the left axis.
C
	  NORMAL_R_Y_AXIS=.TRUE.
	  WRITE(T_OUT,'(A)')' Set Y_AX to R to use with right axis'
	  DO IP=1,NPLTS
	    WRITE(T_OUT,'(I3,'' : '',$)')IP
	    CALL NEW_GEN_IN(WHICH_Y_AX(IP),'Y_AX=')
	    CALL SET_CASE_UP(WHICH_Y_AX(IP),1,1)
	    IF(WHICH_Y_AX(IP) .EQ. 'R')NORMAL_R_Y_AXIS=.FALSE.
	  END DO
C
	  T1=1.0D+30
	  T2=-1.0D+30
	  DO IP=1,NPLTS
	    IF(WHICH_Y_AX(IP) .EQ. 'R')THEN
	      T3=MINVAL(CD(IP)%DATA)
	      T1=MIN(T1,T3)
	      T3=MAXVAL(CD(IP)%DATA)
	      T2=MAX(T2,T3)
	    END IF
	  END DO
C
	  IF(YPAR_R_AX(1) .EQ. 0 .AND. YPAR_R_AX(2) .EQ. 0)THEN
	    YINC_R_AX=SPACING(T1,T2)
	    IYTICK_R_AX=2
	    V1=1000
	    YPAR_R_AX(1)=ABS(YINC_R_AX)*(AINT(T1/ABS(YINC_R_AX)+V1)-V1)
	    YPAR_R_AX(2)=ABS(YINC_R_AX)*(AINT(T2/ABS(YINC_R_AX)-V1)+V1)
	  END IF
	  WRITE(T_OUT,11)T1,T2
	  CALL NEW_GEN_IN(YPAR_R_AX,I,ITWO,'YST,YEND for RHS axis')
C
C Check that limits inserted are not absurd.
C
	  IF(ABS(YPAR_R_AX(2)-YPAR_R_AX(1))/
	1        MAX(ABS(YPAR_R_AX(2)),ABS(YPAR_R_AX(1)))  .LT. 1.0E-08)THEN
	    YPAR_R_AX(1)=0
	    YPAR_R_AX(2)=1.00
	    WRITE(T_OUT,*)'Invalid Y limits - setting default values'
	  END IF
C
	  YINC_R_AX=SPACING(YPAR_R_AX(1),YPAR_R_AX(2))
	  IYTICK_R_AX=2
	  WRITE(T_OUT,*)'Indicate spacing between numberd tickmarks'
	  CALL NEW_GEN_IN(YINC_R_AX,'For Y tickmarks')
	  CALL NEW_GEN_IN(IYTICK_R_AX,'Y minor divisions')
	  IYTICK=MAX(IYTICK_R_AX,1)
C
	  IDY_R_AX=NDEC(YINC_R_AX)
	  IT=NDEC(YPAR_R_AX(1))
	  IF(IT .GT. IDY_R_AX)IDY_R_AX=IT
	  CALL NEW_GEN_IN(IDY_R_AX,'# of Y Dec digits')
	  YNUMST_R_AX=YPAR_R_AX(1)
	  CALL NEW_GEN_IN(YNUMST_R_AX,'Number beginning Y Axis')
	  CALL NEW_GEN_IN(YLABEL_R_AX,'RHS AXIS LABEL')
C
	  GOTO 1000
!
! We do not edit CURVE labels here -- this is done with EDCL
!
	ELSE IF(ANS .EQ. 'L')THEN
	  WRITE(6,*)'Use EDCL to edit curve labels and titles'
	  CALL NEW_GEN_IN(XLABEL,'XLAB')
	  CALL NEW_GEN_IN(YLABEL,'YLAB')
!
! This section allows us to edit a title, and delete a title but keep remaining titles.
!
	  J=1
	  TITLE(N_HEADER+1:)=' '
	  DO WHILE(J .LE. N_TITLE)
	    CALL NEW_GEN_IN(TITLE(J),'TITLE')
	    IF(TITLE(J) .EQ. ' ' .AND. TITLE(MIN(J+1,N_TITLE)) .EQ. ' ')THEN
	      N_HEADER=J-1
	      EXIT
	    ELSE IF(J .EQ. N_TITLE)THEN
	      N_HEADER=N_TITLE
	      EXIT
	    ELSE IF(TITLE(J) .EQ. ' ')THEN
	      TITLE(J:N_TITLE-1)=TITLE(J+1:N_TITLE)
	    ELSE
	      J=J+1
	    END IF
	  END DO
	  CALL NEW_GEN_IN(TITONRHS,'TITONRHS')
	  GOTO 1000
!
	ELSE IF(ANS .EQ. 'WTIT')THEN
	  TITLE_FILENAME='TITLE.SAV'
	  TMP_LOG=.TRUE.
	  DO WHILE(TMP_LOG)
	    CALL NEW_GEN_IN(TITLE_FILENAME,'File for plot titles')
	    INQUIRE(FILE=TITLE_FILENAME,EXIST=TMP_LOG)
	    IF(TMP_LOG)THEN
	      CALL NEW_GEN_IN(TMP_LOG,'Overwrite existing file?')
	      TMP_LOG=.NOT. TMP_LOG
	    END IF
	  END DO
	  OPEN(UNIT=10,FILE=TRIM(TITLE_FILENAME),STATUS='UNKNOWN',IOSTAT=IOS)
	  IF(IOS .EQ .0)THEN
	    DO J=1,N_HEADER                  !TITLE
	      IF(TITLE(J) .EQ. ' ')EXIT
	      WRITE(10,'(A)')TRIM(TITLE(J))
	    END DO
	    CLOSE(UNIT=10)
	  ELSE
	    WRITE(6,*)'Unable to open file ',TRIM(FILNAME)
	  END IF
	  GOTO 1000
!
	ELSE IF(ANS .EQ. 'RTIT')THEN
	  CALL NEW_GEN_IN(TITLE_FILENAME,'FILE')
	  CALL GET_TITLES(TITLE_FILENAME,TITLE,N_TITLE,IOS)
	  IF(IOS .EQ. 0)CALL NEW_GEN_IN(TITONRHS,'TITONRHS')
	  GOTO 1000
!
! Edit curve labels
!
        ELSE IF(ANS .EQ. 'EDCL')THEN
	  CALL NEW_GEN_IN(RESET_CURVE_LAB,'RES_CL')
	  IF(RESET_CURVE_LAB)THEN
	    CALL NEW_GEN_IN(ADD_COMMA,'ADD_COMMA')
	    DO IP=1,NPLTS
	      WRITE(TMP_STR,'(A,I2)')'IP=',IP
	      CALL NEW_GEN_IN(CD(IP)%CURVE_ID,TMP_STR)
	      IF(CD(IP)%CURVE_ID .EQ. '""')CD(IP)%CURVE_ID=' '
	    END DO
	  ELSE
	    DO J=1,N_TITLE
	      IF(CURVE_TITLE(J) .NE. ' ')WRITE(6,'(I3,3X,A)')J,TRIM(CURVE_TITLE(J))
	    END DO
	    DO J=1,N_TITLE
	      IF(CURVE_TITLE(J) .EQ. ' ')EXIT
	      CALL NEW_GEN_IN(CURVE_TITLE(J),'TITLE')
	    END DO
	    DO J=1,N_TITLE-1
	      IF(CURVE_TITLE(J) .EQ. ' ')CURVE_TITLE(J:N_TITLE-1)=CURVE_TITLE(J+1:N_TITLE)
	    END DO
	  END IF
	  GOTO 1000
!
!
! Switch to prevent inialization of data curves on exit from routine.
! By default, the plot count is set to zero on exit, and the data is lost on
! a new entry to the plotting package.
!
	ELSE IF(ANS .EQ. 'NOI')THEN
	  INITIALIZE_ON_EXIT=.NOT. INITIALIZE_ON_EXIT
	  IF(INITIALIZE_ON_EXIT)THEN
	    WRITE(T_OUT,*)'Warning: plot data will be initialized on exit.'
	  ELSE
	    WRITE(T_OUT,*)'Plot data will be retained on exit.'
	  END IF
	  GOTO 1000
!
	ELSE IF(ANS .EQ. 'CC')THEN
	  CALL CHANGE_COLOR(RED,BLUE,GREEN)
	  GOTO 1000
!
	ELSE IF(ANS .EQ. 'CP')THEN
	  CALL CHANGE_PEN(PEN_COL,MAXPEN,NPLTS)
	  GOTO 1000
!
	ELSE IF(ANS .EQ. 'SFP')THEN
	  CALL NEW_GEN_IN(PEN_OFFSET,'First colored pen: 0 [black], 1 [red]')
!
! Reset default pen color. Useful after GP option.
!
	ELSE IF(ANS .EQ. 'RCP')THEN
	  CALL PGSCR(0,.7,.7,.7)     !set color representations
	  CALL PGSCR(1,0.0,0.0,0.0)
	  CALL PGSCR(2,1.0,0.0,0.0)
	  CALL PGSCR(3,0.0,0.0,1.0)
	  CALL PGSCR(4,0.0,0.6,0.3)
	  CALL PGSCR(5,.5,0.0,.7)
	  CALL PGSCR(13,0.0,1.0,1.0)
	  CALL PGSCR(15,.95,.95,.95)
          DO I=0,15                  !Get these + def color representations.
            CALL PGQCR(I,RED(I),GREEN(I),BLUE(I))
          END DO
!
! Set pens better suited for color blindness
!
	ELSE IF(ANS .EQ. 'CBP')THEN
	  CALL PGSCR(0,1.0,1.0,1.0)     !set color representations
	  CALL PGSCR(1,0.0,0.0,0.0)     !set color representations
	  CALL PGSCR(2,0.9,0.6,0.0)
	  CALL PGSCR(3,0.35,0.70,0.90)
	  CALL PGSCR(4,0.0,0.6,0.5)
	  CALL PGSCR(5,0.95,0.9,0.25)
	  CALL PGSCR(6,0.0,0.45,.7)
	  CALL PGSCR(7,0.8,0.40,0.0)
	  CALL PGSCR(8,0.8,0.6,0.7)
	  CALL PGSCR(13,0.0,1.0,1.0)
	  DO I=0,15                  !Get these + def color representations.
	    CALL PGQCR(I,RED(I),GREEN(I),BLUE(I))
	  END DO
!
	ELSE IF(ANS .EQ. 'GP')THEN
!
! Set preferred defaults for grey pens.
!
	  CALL PGSCR(0,1.0,1.0,1.0)     !set color representations
	  CALL PGSCR(1,0.0,0.0,0.0)
	  CALL PGSCR(2,0.0,0.0,0.0)
	  CALL PGSCR(3,0.4,0.4,0.4)
	  CALL PGSCR(4,0.8,0.8,0.8)
	  CALL PGSCR(5,0.2,0.2,0.2)
	  CALL PGSCR(6,0.6,0.6,0.6)
          DO I=0,15                  !Get these + def color representations.
            CALL PGQCR(I,RED(I),GREEN(I),BLUE(I))
          END DO
	  GOTO 1000
	ELSE IF(ANS .EQ. 'SCR')THEN
	  T1=100; T2=0.2
	  CALL PGSCRL(T1,T2)
!
!Option to adjust the shape of the box, and the tick marks tec.
!
	ELSE IF(ANS .EQ. 'N')THEN
!
	  CALL NEW_GEN_IN(EXPCHAR_SCALE,'Expand Char.')
	  CALL NEW_GEN_IN(TICK_FAC_SCALE,'Expand TICK')
!
! Define the Aspect ratio of the plot.
!
	  WRITE(T_OUT,*)'Y/X > 0, 0(Device default) , X/Y < 0'
	  CALL NEW_GEN_IN(ASR,'Aspect Ratio')
!
	  IF(LONG_PLOT)THEN
	    WRITE(T_OUT,*)RED_PEN
	    WRITE(T_OUT,*)' As a long plot, 0.05 and 0.95 may be better for left/right margin'
	    WRITE(T_OUT,*)DEF_PEN
	  END IF
700	  CALL NEW_GEN_IN(MARGINX,I,ITWO,'Left and Right Margins (0:1)')
	  IF((MARGINX(1) .GT. 1) .OR. (MARGINX(1) .LT. 0)) THEN
	    WRITE(T_OUT,*)'Left-Margin is incorrect, try again.'
	    GOTO 700
	  END IF
	  IF((MARGINX(2) .GT. 1) .OR. (MARGINX(2).LT. 0)) THEN
	    WRITE(T_OUT,*)'Right-Margin is incorrect, try again.'
	    GOTO 700
	  END IF
!
710	  CALL NEW_GEN_IN(MARGINY,I,ITWO,'Bottom and Top Margins (0:1)')
	  IF((MARGINY(1) .GT. 1.0) .OR. (MARGINY(1) .LT. 0)) THEN
	    WRITE(T_OUT,*)'Bottom-Margin is incorrect, try again.'
	    GOTO 710
	  END IF
	  IF((MARGINY(2) .GT. 1.0) .OR. (MARGINY(2) .LT. 0)) THEN
	    WRITE(T_OUT,*)'Top-Margin is incorrect, try again.'
	    GOTO 710
	  END IF
	  GOTO 1000
!
	ELSE IF( ANS .EQ. 'W')THEN
	    CALL NEW_GEN_IN(LINE_WGT,I,NPLTS,'Line Weight:: 1,2 etc')
	    IF(LINE_WGT(1) .LT. 0)LINE_WGT(1:NPLTS)=ABS(LINE_WGT(1))
	ELSE IF( ANS .EQ. 'WE')THEN
	  DO IP=1,NPLTS			!Edit individually
	    CALL NEW_GEN_IN(LINE_WGT(IP),'Line weights (1,...,5)')
	  END DO
	ELSE IF( ANS .EQ. 'M')THEN
	  IF(MARK)THEN
            WRITE(T_OUT,*)'Will not mark data points'
	    MARK=.FALSE.
	  ELSE
	    MARK=.TRUE.
            WRITE(T_OUT,*)'Wil now mark data points'
	    CALL NEW_GEN_IN(EXPMARK_SCALE,'Expand Symbol Size')
	    CALL NEW_GEN_IN(MARKER_STYLE,I,NPLTS,'Marker style {+/-}1,...,31)')
	  END IF
	  GOTO 1000
!
	ELSE IF( ANS .EQ. 'C')THEN
	  DO IP=1,NPLTS
850	    WRITE(T_OUT,'(I3,'' : '',$)')IP
	    CALL NEW_GEN_IN(TYPE_CURVE(IP),'TC=')
	    CALL SET_CASE_UP(TYPE_CURVE(IP),IONE,IZERO)
	    IF( TYPE_CURVE(IP) .NE. 'L' .AND.              !Normal line
	1       TYPE_CURVE(IP) .NE. 'E' .AND.              !non-monotonic
	1       TYPE_CURVE(IP) .NE. 'EC' .AND.              !non-monotonic
	1       TYPE_CURVE(IP) .NE. 'B' .AND.              !Broken
	1       TYPE_CURVE(IP) .NE. 'I' .AND.              !Invisible
	1       TYPE_CURVE(IP) .NE. 'V' .AND.              !Verticle lines
	1       TYPE_CURVE(IP) .NE. 'VB' .AND.             !Verticle lines
	1       TYPE_CURVE(IP) .NE. 'A' .AND.              !Hist - X vert
	1       TYPE_CURVE(IP) .NE. 'H' .AND.              !Histogram
	1       TYPE_CURVE(IP) .NE. 'LG' )THEN             !Log(ABS)
!
	        WRITE(T_OUT,*)'Invalid connection specifier. Specifiers are:'
	        WRITE(T_OUT,*)'L --- Normal line'
	        WRITE(T_OUT,*)'E --- non-monotonic'
	        WRITE(T_OUT,*)'B --- Broken'
	        WRITE(T_OUT,*)'I --- Invisible'
	        WRITE(T_OUT,*)'V --- Verticle lines'
	        WRITE(T_OUT,*)'A --- Hist - X vert'
	        WRITE(T_OUT,*)'H --- Histogram'
	        WRITE(T_OUT,*)'LG --- Marked ABS(LG)'
	      GOTO 850
	    END IF
	  END DO
	  GOTO 1000
!
	ELSE IF( ANS .EQ. 'RPO')THEN
	  REVERSE_PLOTTING_ORDER= .NOT. REVERSE_PLOTTING_ORDER
	  IF(REVERSE_PLOTTING_ORDER)THEN
	    WRITE(6,*)'Order of ploting will be reversed'
	  ELSE
	    WRITE(6,*)'Normal order of ploting will be resumed'
	  END IF	
 	ELSE IF( ANS .EQ. 'D')THEN
	  IF(DASH)THEN
            WRITE(T_OUT,*)'Now all solid line plots '
	    LINE_STYLE(1:NPLTS)=1
	    DASH=.FALSE.
	  ELSE
	    WRITE(T_OUT,*)'Input dashed pens for ',NPLTS
	    DO IP=1,NPLTS,5
	      DO J=1,5
	        IF(IP .LE. NPLTS)LINE_STYLE(IP+J-1)=J
	      END DO
	    END DO
	    CALL NEW_GEN_IN(LINE_STYLE,I,NPLTS,
	1      'Dashed pen types (1,...,5)')
	    DASH=.TRUE.
	  END IF
	  GOTO 1000
!
 	ELSE IF(ANS .EQ. 'BRD')THEN
	  IF(DO_BORDER)THEN
            WRITE(T_OUT,*)'Switching off border'
	    DO_BORDER=.FALSE.
	  ELSE
	    WRITE(T_OUT,*)'Switching on border'
	    DO_BORDER=.TRUE.
	  END IF
	  GOTO 1000
!
! Edit each dashed pen separately.
C
 	ELSE IF( ANS .EQ. 'DE')THEN
	  DO IP=1,NPLTS
	    CALL NEW_GEN_IN(LINE_STYLE(IP),'Dashed pen types (1,...,5)')
	  END DO
	  DASH=.TRUE.
	  GOTO 1000
C
	ELSE IF( ANS .EQ. 'B')THEN
	  IF(DO_ERROR)THEN
            WRITE(T_OUT,*)'Swithcing off error bar plotting'
	    DO_ERROR=.FALSE.
	  ELSE
	    DO_ERROR=.TRUE.
	  END IF
	  GOTO 1000
!
	ELSE IF(ANS .EQ. 'RDXL')THEN
	  CALL NEW_GEN_IN(XLAB_FILE,'File with X abscissa values')
	  GOTO 1000
!
	ELSE IF(ANS .EQ. 'RDYL')THEN
	  CALL NEW_GEN_IN(YLAB_FILE,'File with Y ordinate values')
	  GOTO 1000
!
	ELSE IF(ANS .EQ. 'LXY')THEN
	  IF(LOG_AXIS(1:3) .EQ. 'LOG')THEN
            WRITE(T_OUT,*)'Swithcing off ALL LOG axes'
	    LOG_AXIS=' '
	  ELSE
	    LOG_AXIS='LOGXY'
	  END IF
	  GOTO 1000
	ELSE IF( ANS .EQ. 'LX')THEN
	  IF(LOG_AXIS(1:5) .EQ. 'LOGX')THEN
            WRITE(T_OUT,*)'Swithcing off X LOG axes'
	    LOG_AXIS=' '
	  ELSE IF(LOG_AXIS(1:5) .EQ. 'LOGXY')THEN
            WRITE(T_OUT,*)'Swithcing off X LOG axes'
	    LOG_AXIS='LOGY'
	  ELSE IF(LOG_AXIS(1:4) .EQ. 'LOGY')THEN
	    LOG_AXIS='LOGXY'
	  ELSE
	    LOG_AXIS='LOGX'
	  END IF
	  GOTO 1000
	ELSE IF( ANS .EQ. 'LY')THEN
	  IF(LOG_AXIS(1:4) .EQ. 'LOGY')THEN
            WRITE(T_OUT,*)'Swithcing off Y LOG axes'
	    LOG_AXIS=' '
	  ELSE IF(LOG_AXIS(1:5) .EQ. 'LOGXY')THEN
            WRITE(T_OUT,*)'Swithcing off Y LOG axes'
	    LOG_AXIS='LOGX'
	  ELSE IF(LOG_AXIS(1:4) .EQ. 'LOGX')THEN
	    LOG_AXIS='LOGXY'
	  ELSE
	    LOG_AXIS='LOGY'
	  END IF
	  GOTO 1000
!
!
! String handling section.
!
	ELSE IF (ANS .EQ. 'SF')THEN
	  INIT=.FALSE.          !Strings automatically initialized first time.
	  IF(STR)CALL NEW_GEN_IN(INIT,'Initialize string list ?')
	  IF(INIT)THEN
	    DO I=1,MAXSTR
	      FLAGSTR(I)=.FALSE.
	      STR_EXP(I)=1.0	  !Default (changed at editing stage only)
	    END DO
            STR=.FALSE.
	  END IF
!
	  L=INDEX(FILNAME,'.STR')
	  IF(L .NE. 0)FILNAME(L:L+3)='.STR'
	  ISTR=0
923	  CALL NEW_GEN_IN(FILNAME,'Input string file to be read')
	  CALL SET_CASE_UP(FILNAME,1,0)
	  L=LEN_TRIM(FILNAME)
	  IF(L .NE. 0)THEN
	    IF(INDEX(FILNAME,'.') .EQ. 0)FILNAME=FILNAME(:L)//'.STR'
            OPEN(UNIT=33,FILE=TRIM(FILNAME),STATUS='OLD',ERR=923)
	    L=0
            DO ISTR=1,MAXSTR
	      IF( .NOT. FLAGSTR(ISTR))THEN
	        READ(33,*,END=930,ERR=925)LOC(ISTR),XSTR(ISTR),YSTR(ISTR),
	1                   ORIENTATION(ISTR),STR_EXP(ISTR),STR_COL(ISTR),
	1                   STRING(ISTR)
	        FLAGSTR(ISTR)=.TRUE.
	        L=L+1
                GOTO 928
925             WRITE(T_OUT,*)' Error reading string',ISTR
928	        CONTINUE
	      END IF
            END DO
930	    CLOSE(UNIT=33)
	    IF(L .NE. 0)THEN
	      WRITE(T_OUT,932)L,'STRINGS read in from file'
              STR=.TRUE.
	    END IF
          END IF
932	  FORMAT(1X,I3,2X,(A))
!
	  GOTO 1000
!
	ELSE IF(ANS .EQ. 'SE')THEN
!
! First output strings to user, before doing online editing.
!
	  DO I=1,MAXSTR
	    IF(FLAGSTR(I))THEN
	      WRITE(T_OUT,955)I,TRIM(STRING(I))
955	      FORMAT(1X,I2,3X,A)
	    END IF
	  END DO
!
! Allow online editing of strings. NEW_GEN_IN is used so user can
! see current value.
!
	  DO WHILE (0 .EQ. 0)
	    ISTR=0			!Default terminates input
	    CALL NEW_GEN_IN(ISTR,'String number - 0 exits')
	    IF(ISTR .LE. 0 .OR. ISTR .GT. MAXSTR)GOTO 1000
	    WRITE(T_OUT,955)ISTR,TRIM(STRING(ISTR))
	    CALL NEW_GEN_IN(FLAGSTR(ISTR),'FLAG')
	    IF(FLAGSTR(ISTR))THEN
	      CALL NEW_GEN_IN(LOC(ISTR),'LOC')
	      CALL NEW_GEN_IN(XSTR(ISTR),'XPOS')
	      CALL NEW_GEN_IN(YSTR(ISTR),'YPOS')
	      CALL NEW_GEN_IN(ORIENTATION(ISTR),'ORI')
	      CALL NEW_GEN_IN(STR_EXP(ISTR),'EXP')
	      CALL NEW_GEN_IN(STR_COL(ISTR),'COL')
	      CALL NEW_GEN_IN(STRING(ISTR),'STRING')
	      STR=.TRUE.
	    END IF
	  END DO
	  GOTO 1000
!
	ELSE IF (ANS .EQ. 'SC')THEN
	  INIT=.FALSE.          !Strings automatically initialized first time.
	  IF(STR)CALL NEW_GEN_IN(INIT,'Initialize string list ?')
	  STR=.TRUE.
	  IF(INIT)THEN
	    DO I=1,MAXSTR
	      FLAGSTR(I)=.FALSE.
	    END DO
	  END IF
	      WRITE(T_OUT,820)
 820	      FORMAT(1X,'LOC: 1-9, 1=lef bot, 2=cent bot, 3=rit bot')
	      WRITE(T_OUT,821)
 821	      FORMAT(1X,'          4=lef mid, 5=cent mid, 6=rit mid')
	      WRITE(T_OUT,822)
 822	      FORMAT(1X,'          7=lef top, 8=cent top, 9=rit top')
	      WRITE(T_OUT,840)
 840	      FORMAT(1X,'ORI= AN ANGLE FROM 0.0 TO 360')
	  ISTR=1
	  DO WHILE(ISTR .LE. MAXSTR)
	    IF(.NOT. FLAGSTR(ISTR))THEN
	      CURSERR = PGCURS(XSTR(ISTR),YSTR(ISTR),CURSVAL)
	      IF(END_CURS(CURSVAL))GOTO 1000
830	      WRITE(T_OUT,810,ADVANCE='NO')
810	      FORMAT(1X,'Desc(LOC,ORI,''STRING''): ')
	      READ(T_IN,'(A)',ERR=830,END=1000)WK_STR
	      READ(WK_STR,*,ERR=830)LOC(ISTR)
	      IF(LOC(ISTR) .EQ. 0)GOTO 1000
	      READ(WK_STR,*,ERR=830,END=1000)LOC(ISTR),ORIENTATION(ISTR),STRING(ISTR)
	      FLAGSTR(ISTR)=.TRUE.
	    END IF
	    ISTR=ISTR+1
	  END DO
	  WRITE(T_OUT,*)'Maximum number of strings is',MAXSTR
	  GOTO 1000
!
	ELSE IF(ANS .EQ. 'LOC')THEN
	  CALL PGQCS(4,T1,T3)			!T3=Character Height
	  T3=T3*1.5
	  CURSERR = PGCURS(XVAL,YVAL,CURSVAL)
	  IF(END_CURS(CURSVAL))GOTO 70
	  T1=XPAR(1)-(XPAR(2)-XPAR(1))/15.0
	  T2=YPAR(2)+(YPAR(2)-YPAR(1))/15.0
	  CALL MON_NUM(T1,T2,XVAL,1,IDX+2)
	  CALL MON_NUM(T1,T2-T3,YVAL,1,IDY+2)
!
	  DO WHILE(1 .EQ. 1)
	    XVAL_SAV=XVAL
	    YVAL_SAV=YVAL
	    CURSERR = PGCURS(XVAL,YVAL,CURSVAL)
	    CALL PGSCI(0)     ! Background
	    CALL MON_NUM(T1,T2,XVAL_SAV,1,IDX+2)
	    CALL MON_NUM(T1,T2-T3,YVAL_SAV,1,IDY+2)
	    CALL PGSCI(1)     ! Axes
	    IF(END_CURS(CURSVAL))GOTO 70
	    CALL MON_NUM(T1,T2,XVAL,1,IDX+2)
	    CALL MON_NUM(T1,T2-T3,YVAL,1,IDY+2)
	  END DO
70	  CONTINUE
!
	ELSE IF(ANS .EQ. 'SLP')THEN
	  CALL PGQCS(4,T1,T3)			!T3=Character Height
	  T3=T3*1.5
	  T1=XPAR(1)-(XPAR(2)-XPAR(1))/15.0
	  T2=YPAR(2)+(YPAR(2)-YPAR(1))/15.0
!
	  DO WHILE(1 .EQ. 1)
	    CURSERR = PGCURS(XVAL,YVAL,CURSVAL)
	    IF(END_CURS(CURSVAL))EXIT
	    CURSERR = PGCURS(XVAL_SAV,YVAL_SAV,CURSVAL)
	    IF(END_CURS(CURSVAL))EXIT
	    Q=VECPEN(1)
	    CALL PGSCI(IONE)
	    T3=(YVAL_SAV-YVAL)/(XVAL_SAV-XVAL)
	    CALL PGMOVE(XVAL,YVAL)
	    CALL PGDRAW(XVAL_SAV,YVAL_SAV)
	    WRITE(6,*)'Slope of line is:',T3
	  END DO
!
!
!
! Set Line vectors
!
	ELSE IF(ANS .EQ. 'VF')THEN
	  INIT=.FALSE.
	  IF(VEC)CALL NEW_GEN_IN(INIT,'Initialize line drawing list ?')
	  VEC=.TRUE.
	  IF(INIT)THEN
	    DO I=1,MAXVEC
	      FLAGLINE(I)=.FALSE.
	    END DO
	  END IF
!
	  L=INDEX(FILNAME,'.STR')
	  IF(L .NE. 0)FILNAME(L:L+3)='.VEC'
1100	  CALL NEW_GEN_IN(FILNAME,'Line vector file forinput')
	  CALL SET_CASE_UP(FILNAME,1,0)
	  L=LEN_TRIM(FILNAME)
	  ISTR=1
	  IF(L .NE. 0)THEN
	    IF(INDEX(FILNAME,'.') .EQ. 0)FILNAME=FILNAME(:L)//'.VEC'
            OPEN(UNIT=33,FILE=TRIM(FILNAME),STATUS='OLD',ERR=1100)
	    DO ISTR=1,MAXVEC
	      READ(33,*,END=1110,ERR=1100)LINEXST(ISTR),LINEYST(ISTR),
	1              LINEXEND(ISTR),LINEYEND(ISTR),VECPEN(ISTR),VEC_LINE_STYLE(ISTR)
	      FLAGLINE(ISTR)=.TRUE.
	    END DO
1110	    CONTINUE
	    CLOSE(UNIT=33)
	  END IF
	  GOTO 1000
!
	ELSE IF(ANS .EQ. 'VE')THEN
!
! First output vectors to user, before doing online editing.
! NEW_GEN_IN is used so user can see current value. This section
! also used to input vectors by hand.
!
	  DO I=1,MAXVEC
	    IF(FLAGLINE(I))THEN
	      WRITE(T_OUT,'(1X,I2,4(3X,E14.6),2I4)')I,
	1               LINEXST(I),LINEXEND(I),
	1               LINEYST(I),LINEYEND(I),VECPEN(I),VEC_LINE_STYLE(I)
	    END IF
	  END DO
!
	  DO WHILE (0 .EQ. 0)
	    ISTR=0			!Default terminates input
	    CALL NEW_GEN_IN(ISTR,'Vector number - 0 exits')
	    IF(ISTR .LE. 0 .OR. ISTR .GT. MAXVEC)GOTO 1150
	    CALL NEW_GEN_IN(FLAGLINE(ISTR),'FLAG')
	    IF(FLAGLINE(ISTR))THEN
	      CALL NEW_GEN_IN(LINEXST(ISTR),'XST')
	      CALL NEW_GEN_IN(LINEYST(ISTR),'YST')
	      CALL NEW_GEN_IN(LINEXEND(ISTR),'XEND')
	      CALL NEW_GEN_IN(LINEYEND(ISTR),'YEND')
	      CALL NEW_GEN_IN(VEC_LINE_STYLE(ISTR),'1,2,3 or 4')
	    END IF
	  END DO
 1150	  VEC=.TRUE.
	  QUERYFLAG=.FALSE.
	  CALL NEW_GEN_IN(QUERYFLAG,'Change vector pens?')
	  IF(QUERYFLAG)CALL VECTORPEN(VECPEN,MAXVEC,FLAGLINE)
	  GOTO 1000
!
	ELSE IF(ANS .EQ. 'VC')THEN
	  INIT=.FALSE.
	  IF(VEC)CALL NEW_GEN_IN(INIT,'Initialize line drawing list?')
	  VEC=.TRUE.
	  IF(INIT)THEN
	    DO I=1,MAXVEC
	      FLAGLINE(I)=.FALSE.
	    END DO
	  END IF
	  ISTR=1
	  DO WHILE (ISTR .LT. MAXVEC)
	    IF( .NOT. FLAGLINE(ISTR) )THEN
	      CURSERR = PGCURS(LINEXST(ISTR),LINEYST(ISTR),CURSVAL)
	      IF(END_CURS(CURSVAL))GOTO 1000
	      CURSERR = PGCURS(LINEXEND(ISTR),LINEYEND(ISTR),CURSVAL)
	      IF(END_CURS(CURSVAL))GOTO 1000
	      IF(LINEXST(ISTR) .EQ. LINEXEND(ISTR) .AND.
	1        LINEYST(ISTR) .EQ. LINEYEND(ISTR) )THEN
	         ISTR=1000
	      ELSE
	        FLAGLINE(ISTR)=.TRUE.
	      END IF
	    END IF
	    ISTR=ISTR+1
	  END DO
	  GOTO 1000
!
	ELSE IF(ANS .EQ. 'RID')THEN
	  N_LINE_IDS=0
	  WRITE(6,'(A)')' '
	  CALL NEW_GEN_IN(ID_FILNAME,'File with line IDs -case sensitive')
	  CALL NEW_GEN_IN(TAU_CUT,'Omit lines with central optical depth <')
	  CALL NEW_GEN_IN(LINE_CUT_PARAM,'Omit lines whose central intensity differ by less than this amount')
	  CALL NEW_GEN_IN(AIR_WAVELENGTHS,'Air wavelengths?')
	  CALL NEW_GEN_IN(KEEP_YAXIS_LIMITS,'Keep same Y axis limits?')
	  CALL NEW_GEN_IN(ID_EXPCHAR,'Factor to scale size of ID')
	  CALL NEW_GEN_IN(XUNIT,'Ang, um, or nm?')
	  CALL NEW_GEN_IN(NO_DEC_DIGITS,'# of digits afer . :  0 to 9')
	  IF(NO_DEC_DIGITS .LE. 0 .OR. NO_DEC_DIGITS .GT.  9)NO_DEC_DIGITS=1
	  WRITE(6,'(A)')' '
!
	  XT=XPAR
	  IF(XUNIT .EQ. 'nm')XT=10.0D0*XPAR    			!Convert range to Ang
	  IF(XUNIT .EQ. 'um')XT=1.0D+04*XPAR
	  CALL RD_LINE_IDS(ID_FILNAME,AIR_WAVELENGTHS,XT)
	  CALL ADJUST_ID_WAVES(XUNIT)
	  GOTO 1000
!
	ELSE IF(ANS .EQ. 'SID')THEN
	  CALL NEW_GEN_IN(ID_VEC_BEG,'Location to start ID line')
	  CALL NEW_GEN_IN(ID_VEC_END,'Location to end ID line')
	  CALL NEW_GEN_IN(ID_EXPCHAR,'Factor to scale size of ID')
	  CALL NEW_GEN_IN(ID_LINE_PEN,'Color of pen for line labels')
	  CALL NEW_GEN_IN(USE_DEF_OFFSET,'Use same offset for all lines')
	  N_OMIT_ID=0
!
	  WRITE(T_OUT,*)' '
	  WRITE(T_OUT,*)RED_PEN,'Species name are case sensitive and have no spaces'
	  WRITE(T_OUT,*)'Use SiIII etc',DEF_PEN
	  WRITE(T_OUT,*)' '
	  DO I=1,10
	    OMIT_ID(I)='ZZ'
	    CALL NEW_GEN_IN(OMIT_ID(I),'Species ID to OMIT from labels')
	    IF(OMIT_ID(I) .EQ. 'ZZ')EXIT
	    N_OMIT_ID=I
	  END DO
	  N_INC_ID=0
	  DO I=1,10
	    INC_ID(I)='ZZ'
	    CALL NEW_GEN_IN(INC_ID(I),'Species ID to INC in labels')
	    IF(INC_ID(I) .EQ. 'ZZ')EXIT
	    N_INC_ID=I
	  END DO
!
	ELSE IF(ANS .EQ. 'REW')THEN
	  N_EW_IDS=0
	  J=0
	  CALL NEW_GEN_IN(EW_LINE_ID,'File with line EWs')
	  CALL NEW_GEN_IN(EW_CUT,'Omit lines with ABS|EW| < ?')
	  CALL NEW_GEN_IN(XUNIT,'Ang, um, or nm?')
	  CALL NEW_GEN_IN(NO_DEC_DIGITS,' 0 to 9')
!
!	  WRITE(6,*)BLUE_PEN
!	  WRITE(6,*)'Can label individual lines, of use vertical bars to indicate contributing species'
!	  WRITE(6,*)'Lines with an EW of the scaling EW are have 25% of the plot height'
!	  WRITE(6,*)DEF_PEN
!
!	  CALL NEW_GEN_IN(AIR_WAVELENGTHS,'Air wavelengths?')
	  I=33; XT=XPAR
	  IF(XUNIT .EQ. 'nm')XT=10.0D0*XPAR  		!Convert range to Ang.
	  IF(XUNIT .EQ. 'um')XT=1.0D+04*XPAR
	  CALL RD_EW_IDS(XT,EW_LINE_ID,I,T_OUT)
	  CALL ADJUST_ID_WAVES(XUNIT)
!
! 
	ELSE IF (ANS .EQ. 'REP')THEN
	  CALL PG_REPLACE_DATA_INT(
	1             XPAR,XINC,XNUMST,IXTICK,IDX,
	1             YPAR,YINC,YNUMST,IYTICK,IDY,
	1             TICK_FAC,EXPCHAR,
	1             XLABEL,YLABEL,TITONRHS,
	1             LOG_AXIS,OPTION,NORMAL_R_Y_AXIS,
	1             XLAB_FILE,YLAB_FILE,
	1             PEN_COL,PEN_OFFSET,
	1             REVERSE_PLOTTING_ORDER)
!
	ELSE IF (ANS .EQ. 'FREP')THEN
	  FILNAME='CONT_NODES'
	  CALL NEW_GEN_IN(FILNAME,'File with cursor nodes')
	  OPEN(UNIT=10,FILE=TRIM(FILNAME),STATUS='OLD',ACTION='READ')
	    READ(10,*)K,J
	    DO L=1,K
	      READ(10,*)XCUR(1),XCUR(2)
	      L_CHAN(1)=GET_INDX_SP(XCUR(1),CD(PLOT_ID)%XVEC,NPTS(PLOT_ID))
	      L_CHAN(2)=GET_INDX_SP(XCUR(2),CD(PLOT_ID)%XVEC,NPTS(PLOT_ID))
	      SLOPE=(CD(PLOT_ID)%DATA(L_CHAN(2))-CD(PLOT_ID)%DATA(L_CHAN(1)))/(XCUR(2)-XCUR(1))
	      DO I=1,NPTS(PLOT_ID)
	        IF( (CD(PLOT_ID)%XVEC(I)-XCUR(1))*(CD(PLOT_ID)%XVEC(I)-XCUR(2)) .LT. 0.0D0)THEN
                  CD(PLOT_ID)%DATA(I)=CD(PLOT_ID)%DATA(L_CHAN(1))+SLOPE*(CD(PLOT_ID)%XVEC(I)-XCUR(1))
	        END IF
	      END DO
	    END DO
	  CLOSE(UNIT=10)
!
! Allow manipulation of a plot using cursors.
!
	ELSE IF (ANS .EQ. 'MODC')THEN
	  IF(NPLTS .EQ. 1)THEN
	    PLOT_ID=1; IP=2
	  ELSE
	    PLOT_ID=NPLTS; IP=NPLTS
	  END IF
	  CALL NEW_GEN_IN(PLOT_ID,'Plot for simple modification:')
	  CALL NEW_GEN_IN(IP,'Output plot:')
	  CALL MODIFY_CURVE(PLOT_ID,IP)
	  TYPE_CURVE(IP)=TYPE_CURVE(PLOT_ID)
!
! Define a crude straight line continuum about a line, so the EW can be
! computed.
!
	ELSE IF (ANS .EQ. 'DC')THEN
	  IF(NPLTS .EQ. 1)THEN
	    PLOT_ID=1
	  ELSE
	    CALL NEW_GEN_IN(PLOT_ID,'Plot to determine EW for:')
	  END IF
!
          WRITE(6,*)RED_PEN
          WRITE(6,*)'CURSOR -- input X and Y values'
          WRITE(6,*)'CURX -- input X only, average Y between the X values'
          WRITE(6,*)'FILENAME -- first line -- number of data pairs '
          WRITE(6,*)'FILENAME -- input X1 and X2 on each line; averages Y between the X values'
          WRITE(6,*)DEF_PEN
!
	  DC_INPUT_OPTION='CURSOR'
	  CALL NEW_GEN_IN(DC_INPUT_OPTION,'CURSOR, CURX, or name of file')
	  DC_CURVE_OPTION='MONCUB'
	  IF(DC_INPUT_OPTION .EQ. 'CURSOR')DC_CURVE_OPTION='LINEAR'
	  CALL NEW_GEN_IN(DC_CURVE_OPTION,'Curve type -- LINEAR of MONCUB')
!
	  IF(ALLOCATED(CONT))DEALLOCATE(CONT)
	  ALLOCATE(CONT(NPTS(PLOT_ID)))
	  CONT=CD(PLOT_ID)%DATA
	  CD(PLOT_ID)%CURVE_ID=' '
!
	  IP=NPLTS+1
	  CALL NEW_GEN_IN(IP,'Output plot?')
	  CALL PG_DEF_CONTINUUM(CONT,PLOT_ID,IP,DC_INPUT_OPTION,DC_CURVE_OPTION,IOS)
	  IF(IOS .EQ. 0)CONTINUUM_DEFINED=.TRUE.
          IF(IP .NE. 0)TYPE_CURVE(IP)='L'
	  GOTO 1000
!
! This option allows a curve to be defined using currors. Can be used to define
! a continuum. Can be called multiples times.
!
	ELSE IF (ANS .EQ. 'MCN')THEN
	  IF(IP_CONT .EQ. 0)THEN
	    IP_CONT=NPLTS+1
	    NPLTS=NPLTS+1
	    NPTS(IP_CONT)=0
	    ERR(IP_CONT)=.FALSE.
	    TYPE_CURVE(IP_CONT)='L'
	  END IF
	  CALL NEW_GEN_IN(IP_CONT,'Plot with continuum nodes')
	  IF(IP_CONT .GT. NPLTS)THEN
	    IP_CONT=NPLTS+1
	    NPLTS=NPLTS+1
	    NPTS(IP_CONT)=0
	    WRITE(6,*)'Using ',IP_CONT,' for continuum definition'
	  END IF
	  CALL PG_MOD_CONT_NODES(IP_CONT)
          TYPE_CURVE(IP_CONT)='L'
	  GOTO 1000
!
	ELSE IF (ANS .EQ. 'CONT')THEN
!
	  WRITE(6,'(A)')' '
	  WRITE(6,'(A)')' Define continuum using X & Y cursor locations'
	  WRITE(6,'(A)')' Striaght line fit across cursor band'
	  WRITE(6,'(A)')' Routine can be called multiple times'
	  WRITE(6,'(A)')' Use E to exit cursor selection.'
	  WRITE(6,'(A)')' '
	  IF(NPLTS .EQ. 1)THEN
	    PLOT_ID=1
	  ELSE
	    IOS=-1
	    DO WHILE(IOS .NE. 0)
	      IOS=0
	      CALL NEW_GEN_IN(PLOT_ID,'Plot to determine EW for:')
	      IF(PLOT_ID .LT. 0 .OR. PLOT_ID .GT. NPLTS)THEN
	        IOS=-1
	        WRITE(6,*)'Error - invalid plot number - try again'
	      END IF
	    END DO
	  END IF
	  IF(ALLOCATED(CONT))DEALLOCATE(CONT)
	  ALLOCATE(CONT(NPTS(PLOT_ID)))
	  TMP_LOG=.FALSE.
	  CALL NEW_GEN_IN(TMP_LOG,'Zero plot outside continuum definiton?')
	  IF(TMP_LOG)THEN
	    CONT=0.0D0
	  ELSE
	    CONT=CD(PLOT_ID)%DATA
	  END IF
!
	  CURSERR = PGCURS(XCUR(2),YCUR(2),CURSVAL)
	  DO J=1,30
	    XCUR(1)=XCUR(2); YCUR(1)=YCUR(2)
	    CURSERR = PGCURS(XCUR(2),YCUR(2),CURSVAL)
	    IF(END_CURS(CURSVAL))EXIT
	    SLOPE=(YCUR(2)-YCUR(1))/(XCUR(2)-XCUR(1))
	    DO I=1,NPTS(PLOT_ID)
	      IF( (CD(PLOT_ID)%XVEC(I)-XCUR(1))*(CD(PLOT_ID)%XVEC(I)-XCUR(2)) .LT. 0)THEN
                CONT(I)=YCUR(1)+SLOPE*(CD(PLOT_ID)%XVEC(I)-XCUR(1))
	      END IF
	    END DO
	  END DO
!
	  CALL PGLINE(NPTS(PLOT_ID),CD(PLOT_ID)%XVEC,CONT)
	  CONTINUUM_DEFINED=.TRUE.
	  IP=NPLTS+1;    IOS=-1
	  DO WHILE(IOS .NE. 0)
	    IOS=0
	    CALL NEW_GEN_IN(IP,'Output plot?')
	    IF(IP .LT. 0 .OR. IP .GT. NPLTS+1)THEN
	      IOS=-1
	      WRITE(6,*)'Error - invalid plot number - try again'
	    END IF
	  END DO
          ERR(IP)=.FALSE.
	  TYPE_CURVE(IP)='L'
	  IF(IP .EQ. PLOT_ID)THEN
            CD(IP)%DATA=CONT
	  ELSE
	    IF(ALLOCATED(CD(IP)%XVEC))THEN
	      DEALLOCATE (CD(IP)%XVEC)
	      DEALLOCATE (CD(IP)%DATA)
	    END IF
	    I=NPTS(PLOT_ID)
	    ALLOCATE (CD(IP)%XVEC(I),STAT=IOS)
	    IF(IOS .EQ. 0)ALLOCATE (CD(IP)%DATA(I),STAT=IOS)
	    IF(IOS .NE. 0)THEN
	      WRITE(T_OUT,*)'Error: unable to allocate new data vectors'
	      WRITE(T_OUT,*)'IOS=',IOS
	      STOP
	    END IF
            CD(IP)%XVEC=CD(PLOT_ID)%XVEC
	    CD(IP)%DATA=CONT
            NPTS(IP)=I
            IF(IP .GT. NPLTS)NPLTS=IP
	  END IF
	  CD(IP)%CURVE_ID=' '
	  GOTO 1000
!
	ELSE IF(ANS .EQ. 'CGF' .OR. ANS .EQ. 'FGF')THEN
!
	  WRITE(6,'(A)')BLUE_PEN
	  WRITE(6,'(A)')' Calling a routine to do Gaussian Fitting using cursors to define lines'
	  WRITE(6,'(A)')' When using the CGF fit routine to define a file for FGF you should do the following:'
	  WRITE(6,'(A)')'       (a) Fit regions containing the minimum number of lines. The more lines in the fit,'
	  WRITE(6,'(A)')'              the more likely a subsequent fit is to mess up.'
	  WRITE(6,'(A)')'       (b) Choose the low turbulence model for the first fit, and also plot the high'
	  WRITE(6,'(A)')'              turbulence model so you can see how broad the lines are. This will help you '
	  WRITE(6,'(A)')'              choose a fitting region that has necessary lambda coverage.'
	  WRITE(6,'(A)')'       (c) In the FGF routines the central line wavelengths are held fixed .'
	  WRITE(6,'(A)')DEF_PEN
!
	  DRAW_GAUSS_HARD=.TRUE.
	  IF(ANS .EQ. 'CGF')GF_OPTION='CURSOR'
	  IF(ANS .EQ. 'FGF')GF_OPTION='FILE'
	  CALL GAUSS_FIT(GF_OPTION,
	1             XPAR,XINC,XNUMST,IXTICK,IDX,
	1             YPAR,YINC,YNUMST,IYTICK,IDY,
	1             TICK_FAC,EXPCHAR,
	1             XLABEL,YLABEL,TITONRHS,
	1             LOG_AXIS,OPTION,NORMAL_R_Y_AXIS,
	1             XLAB_FILE,YLAB_FILE,
	1             PEN_COL,PEN_OFFSET,
	1             REVERSE_PLOTTING_ORDER)
	  GOTO 1000
!
	ELSE IF(ANS .EQ. 'SEW')THEN
	  T1=0.002D0
	  CALL SIMP_EW(CD(1)%DATA(1),CD(1)%XVEC(1),NPTS(1),XPAR(1),XPAR(2),T1)
	  GOTO 1000
!
	ELSE IF(ANS .EQ. 'DG')THEN
!	  CALL DRAW_GAUS(L_FALSE)
	  GOTO 1000
!
	ELSE IF(ANS .EQ. 'SP')THEN
	  CALL STEP_PLOT_V1(
	1             XPAR,XINC,XNUMST,IXTICK,IDX,
	1             YPAR,YINC,YNUMST,IYTICK,IDY,
	1             TICK_FAC,EXPCHAR,
	1             XLABEL,YLABEL,TITONRHS,
	1             LOG_AXIS,OPTION,NORMAL_R_Y_AXIS,
	1             XLAB_FILE,YLAB_FILE,
	1             LINE_STYLE,LINE_WGT,
	1             PEN_COL,PEN_OFFSET,
	1             REVERSE_PLOTTING_ORDER)
!
	ELSE IF(ANS .EQ. 'CEW')THEN
	  CALL DO_CURSOR_EW_V3(
	1             XPAR,XINC,XNUMST,IXTICK,IDX,
	1             YPAR,YINC,YNUMST,IYTICK,IDY,
	1             TICK_FAC,EXPCHAR,
	1             XLABEL,YLABEL,TITONRHS,
	1             LOG_AXIS,OPTION,NORMAL_R_Y_AXIS,
	1             XLAB_FILE,YLAB_FILE,
	1             PEN_COL,PEN_OFFSET,
	1             REVERSE_PLOTTING_ORDER)
!
! This option is designed to use cursors to measure the EW of a broad line (e.g. H-gamma)
! in a model and its  corresponding observation. The observations are normalized to the
! model using a straight line fit, and done to minimize the chi^2 difference between the
! model and observations.
!
	ELSE IF(ANS .EQ. 'CBAL')THEN
	  I=NPLTS
	  CALL DO_CURSOR_BALMER(
	1             XPAR,XINC,XNUMST,IXTICK,IDX,
	1             YPAR,YINC,YNUMST,IYTICK,IDY,
	1             TICK_FAC,EXPCHAR,
	1             XLABEL,YLABEL,TITONRHS,
	1             LOG_AXIS,OPTION,NORMAL_R_Y_AXIS,
	1             XLAB_FILE,YLAB_FILE,
	1             PEN_COL,PEN_OFFSET,
	1             REVERSE_PLOTTING_ORDER)
	  IF(NPLTS .GT. I)TYPE_CURVE(NPLTS)='L'
!
	ELSE IF(ANS .EQ. 'FBAL')THEN
	  IP=1; CALL NEW_GEN_IN(IP,'Observational plot ID?')
	  CALL NEW_GEN_IN(BALMER_INPUT_FILE,'File wth data for auto chi^2 computation')
	  CALL DO_FILE_BALMER_V1(BALMER_INPUT_FILE,IP)
!
	ELSE IF(ANS .EQ. 'FEW')THEN
	  CALL DO_FILE_EW_V1(' ')
!
	ELSE IF(ANS .EQ. 'EW')THEN
!
	  IF(FIRST_EW)THEN
	    OPEN(UNIT=LU_EW,FILE='EW_FR_SPEC_PLT',STATUS='UNKNOWN',POSITION='APPEND')
	    FIRST_EW=.FALSE.
	    WRITE(LU_EW,'(A10,4A15)')'   Plot ID','  EW(X units)','  X(centroid)','     X(start)','       X(end)'
	    FLUSH(LU_EW)
	  END IF
!
	  IF(CONTINUUM_DEFINED)THEN
!
	    DO WHILE(1 .EQ. 1)
	      CURSERR = PGCURS(XCUR(1),YCUR(1),CURSVAL)
	      IF(CURSVAL .EQ. 'e' .OR. CURSVAL .EQ. 'E')GOTO 1000
	      IF(END_CURS(CURSVAL))GOTO 1000
	      XCUR(2)=XCUR(1); YCUR(2)=YCUR(1)
	      CURSERR = PGCURS(XCUR(2),YCUR(2),CURSVAL)
	      IF(CURSVAL .EQ. 'e' .OR. CURSVAL .EQ. 'E')GOTO 1000
	      IF(END_CURS(CURSVAL))GOTO 1000
!
! Find nearest channels to curser positions.
!
	      IP=PLOT_ID
	      L_CHAN(1)=GET_INDX_SP(XCUR(1),CD(IP)%XVEC,NPTS(IP))
	      L_CHAN(2)=GET_INDX_SP(XCUR(2),CD(IP)%XVEC,NPTS(IP))
!
! Draw in line limits
!
	      CALL PGMOVE(CD(IP)%XVEC(L_CHAN(1)),YPAR(1))
	      CALL PGDRAW(CD(IP)%XVEC(L_CHAN(1)),CONT(L_CHAN(1)))
	      CALL PGMOVE(CD(IP)%XVEC(L_CHAN(2)),YPAR(1))
	      CALL PGDRAW(CD(IP)%XVEC(L_CHAN(2)),CONT(L_CHAN(2)))
!
	      EW=0.0
              CENTROID=0.0
	      IF(L_CHAN(1) .LT. L_CHAN(2))THEN
	        T1=0.5D0
	      ELSE
	        K=L_CHAN(1)
	        L_CHAN(1)=L_CHAN(2)
	        L_CHAN(2)=K
	        T1=-0.5D0
	      END IF
	      DO I=L_CHAN(1),L_CHAN(2)-1
	        EW=EW+T1*( (CD(IP)%DATA(I)-CONT(I))/CONT(I) +
	1               (CD(IP)%DATA(I+1)-CONT(I+1))/CONT(I+1) )*
	1               (CD(IP)%XVEC(I+1)-CD(IP)%XVEC(I))
	        CENTROID=CENTROID+T1*( CD(IP)%XVEC(I)*(CD(IP)%DATA(I)-CONT(I))/CONT(I)
	1               + CD(IP)%XVEC(I+1)*(CD(IP)%DATA(I+1)-CONT(I+1))/CONT(I+1) )*
	1                 (CD(IP)%XVEC(I+1)-CD(IP)%XVEC(I))
	      END DO
	      IF(EW .NE. 0.0)THEN
	        CENTROID=CENTROID/EW
	      ELSE
	        CENTROID=0.0D0
	      END IF
	      WRITE(T_OUT,'(A,I3,5X,A,1PE10.3,A,A,ES14.6,2ES14.4)')
	1             ' Plot ID=',IP,'EW=',EW,' X(units);','   Centroid=',CENTROID,
	1             CD(IP)%XVEC(L_CHAN(1)),CD(IP)%XVEC(L_CHAN(2))
	      WRITE(LU_EW,'(6X,I4,4ES15.6)')IP,EW,CENTROID,CD(IP)%XVEC(L_CHAN(1)),CD(IP)%XVEC(L_CHAN(2))
	      FLUSH(LU_EW)
	      XCUR(1)=XCUR(2); YCUR(1)=YCUR(2)
	    END DO
	    GOTO 1000
!
	  ELSE
!
	    QUERYFLAG=.FALSE.
	    CALL NEW_GEN_IN(QUERYFLAG,'Compute area?')
	    CALL NEW_GEN_IN(WRITE_COMMENT,'Write a comment after each EW measurment?')
	    IF(.NOT. QUERYFLAG)THEN
              WRITE(T_OUT,*)RED_PEN
	      WRITE(T_OUT,*)'Continuum not defined - use DC to define non-unit continnum'
	      WRITE(T_OUT,*)'Assuming Ic=1 (i.e. normalized data)'
	      WRITE(T_OUT,*)'Use e or E to exit line selection'
              WRITE(T_OUT,*)DEF_PEN
!
! Show the continuum
!
	      T3=1.0
	      CALL PGMOVE(XPAR(1),T3)
	      CALL PGDRAW(XPAR(2),T3)
	    END IF
!
! Get the location of the line.
!
	    DO WHILE(1 .EQ. 1)
	      CURSERR = PGCURS(XCUR(1),YCUR(1),CURSVAL)
	      IF(CURSVAL .EQ. 'e' .OR. CURSVAL .EQ. 'E')GOTO 1000
	      IF(END_CURS(CURSVAL))GOTO 1000
	      XCUR(2)=XCUR(1); YCUR(2)=YCUR(1)
	      CURSERR = PGCURS(XCUR(2),YCUR(2),CURSVAL)
	      IF(CURSVAL .EQ. 'e' .OR. CURSVAL .EQ. 'E')GOTO 1000
	      IF(END_CURS(CURSVAL))GOTO 1000
!
! We do all plot provided the cover the range indicated by the cursors.
! Find nearest channels to curser positions.
!
! Draw in line limits
!
	      IF(QUERYFLAG)THEN
	        T3=0.3*(YPAR(2)-YPAR(1))+YPAR(1)
	      ELSE
	        T3=1.0D0
	      END IF
	      CALL PGMOVE(XCUR(1),YPAR(1))
	      CALL PGDRAW(XCUR(1),T3)
	      CALL PGMOVE(XCUR(2),YPAR(1))
	      CALL PGDRAW(XCUR(2),T3)
!
	      DO IP=1,NPLTS
	        T1=XCUR(1)-CD(IP)%XVEC(1)
	        T2=XCUR(1)-CD(IP)%XVEC(NPTS(IP))
	        T3=XCUR(2)-CD(IP)%XVEC(1)
	        T4=XCUR(2)-CD(IP)%XVEC(NPTS(IP))
	        IF(T1*T2 .LT. 0 .AND. T3*T4 .LT. 0)THEN
	          L_CHAN(1)=GET_INDX_SP(XCUR(1),CD(IP)%XVEC,NPTS(IP))
	          L_CHAN(2)=GET_INDX_SP(XCUR(2),CD(IP)%XVEC,NPTS(IP))
!
	          IF(QUERYFLAG)THEN
!
! Compute area of marked region.
!
	            EW=0.0
                    CENTROID=0.0
	            K=1
	            IF(L_CHAN(1) .LT. L_CHAN(2))THEN
	              T1=0.5D0
	            ELSE
	              K=L_CHAN(1)
	              L_CHAN(1)=L_CHAN(2)
	              L_CHAN(2)=K
	              T1=-0.5D0
	            END IF
	            DO I=L_CHAN(1),L_CHAN(2)-1
	              EW=EW+T1*( CD(IP)%DATA(I) + CD(IP)%DATA(I+1) )*
	1                 (CD(IP)%XVEC(I+1)-CD(IP)%XVEC(I))
	              CENTROID=CENTROID+T1*( CD(IP)%XVEC(I)*CD(IP)%DATA(I)
	1                 + CD(IP)%XVEC(I+1)*CD(IP)%DATA(I+1) )*
	1                  (CD(IP)%XVEC(I+1)-CD(IP)%XVEC(I))
	            END DO
	            IF(EW .NE. 0.0)THEN
	              CENTROID=CENTROID/EW
	            ELSE
	              CENTROID=0.0
	            END IF
	            WRITE(T_OUT,'(A,I3,3X,5X,A,1PE10.3,A,A,1PE14.6)')' Plot ID=',IP,'  AREA=',EW,' X(units);','  Centroid=',CENTROID
	            WRITE(LU_EW,'(A,I3,3X,5X,A,1PE10.3,A,A,1PE14.6)')' Plot ID=',IP,'  AREA=',EW,' X(units);','  Centroid=',CENTROID
	        ELSE
!
! Compute EW
!
	          EW=0.0
                  CENTROID=0.0
	          DO I=L_CHAN(1),L_CHAN(2)-1
	            EW=EW+0.5*( (CD(IP)%DATA(I)-1.0) + (CD(IP)%DATA(I+1)-1.0) )*
	1                 (CD(IP)%XVEC(I+1)-CD(IP)%XVEC(I))
	            CENTROID=CENTROID+0.5*( CD(IP)%XVEC(I)*(CD(IP)%DATA(I)-1.0)
	1               + CD(IP)%XVEC(I)*(CD(IP)%DATA(I+1)-1.0) )*
	1                 (CD(IP)%XVEC(I+1)-CD(IP)%XVEC(I))
	          END DO
	          IF(EW .NE. 0.0)THEN
	            CENTROID=CENTROID/EW
	          ELSE
	            CENTROID=0.0
	          END IF
!
	          WRITE(T_OUT,'(A,I3,5X,A,1PE10.3,A,A,ES14.6,2ES14.4)')
	1            ' Plot ID=',IP,'EW=',EW,' X(units);','   Centroid=',CENTROID,
	1                  CD(IP)%XVEC(L_CHAN(1)),CD(IP)%XVEC(L_CHAN(2))
	          IF(IP .EQ. 1 .AND. WRITE_COMMENT)THEN
	            WRITE(LU_EW,'(A)')' '
	            TMP_STR=' '
	            CALL NEW_GEN_IN(TMP_STR,'Comment=')
	            IF(TMP_STR .NE. ' ')THEN
	              WRITE(LU_EW,'(A)')TRIM(TMP_STR)
	              WRITE(LU_EW,'(A)')' '
	              TMP_STR=' '
	              CALL NEW_GEN_IN(TMP_STR,'Comment=')
	              IF(TMP_STR .NE. ' ')WRITE(LU_EW,'(A)')TRIM(TMP_STR)
	            END IF
	          END IF
	          WRITE(LU_EW,'(6X,I4,4ES15.6)')IP,EW,CENTROID,CD(IP)%XVEC(L_CHAN(1)),CD(IP)%XVEC(L_CHAN(2))
	          FLUSH(LU_EW)
	        END IF
	      END IF
	      XCUR(1)=XCUR(2); YCUR(1)=YCUR(2)
	    END DO
	    END DO		!Multiple plots
	    GOTO 1000
	  END IF
	ELSE IF(ANS .EQ. 'EWG')THEN
	  CALL DO_MANY_EW(TYPE_CURVE,NPLTS,MAX_PLTS)
!
	ELSE IF(ANS .EQ. 'LP')THEN
	  IF(LONG_PLOT)THEN
	    WRITE(T_OUT,*)'Resuming normal hard copy mode'
	    LONG_PLOT=.FALSE.
	  ELSE
	    WRITE(T_OUT,*)'Setting hard-copy to produce long plots'
	    WRITE(T_OUT,*)'This option only work CPS as the device'
	    LONG_PLOT=.TRUE.
	    CALL NEW_GEN_IN(LENGTH_OF_HC_PLOT,'Length of plot surface in cm''s')
	    LP_ASR=2.54D0*8.0D0/LENGTH_OF_HC_PLOT
	    WRITE(T_OUT,*)'Default aspect ratio of pot surface is',LP_ASR
	  END IF
!
!
!
	ELSE IF(ANS .EQ. 'Z' .OR. ANS .EQ. 'ZN')THEN
	  IF (FIRST_HARD .OR. ANS .EQ. 'ZN') THEN
	    WRITE(T_OUT,*)RED_PEN,' '
	    WRITE(T_OUT,*)'Choose a post-script device and file for printing [file/dev] '
	    WRITE(T_OUT,*)'A sensible name format is xxx_1.ps/cps'
	    WRITE(T_OUT,*)'Use a . only for the file extension'
	    WRITE(T_OUT,*)DEF_PEN,' '
	    CALL NEW_GEN_IN(PRINTER,'File and printer: enter ? for list')
	    BEG=PGOPEN(PRINTER)
	    FIRST_HARD=.FALSE.
	    CALL DEFINE_MORE_PENS(MAXPEN)
	  ELSE
	    CALL NEW_GEN_IN(HARD_FILE,'Plot file')
	    PRINTER=TRIM(HARD_FILE)//'/'//HARD_TYPE
	    BEG=PGOPEN(PRINTER)
	    CALL DEFINE_MORE_PENS(MAXPEN)		!Pen definitions are not saved.
	  END IF
!
! Save hard device and set default file name for net plot.
!
	  CALL PGQINF('TYPE',HARD_TYPE,LENGTH)
	  CALL PGQINF('FILE',HARD_FILE,LENGTH)
	  CALL GET_PGI_FILE_NUMBER(HARD_CNT,HARD_FILE)
	  CALL UPDATE_PGI_FILE_NAME(HARD_CNT,HARD_FILE)
	  PRINTER=TRIM(HARD_FILE)//'/'//HARD_TYPE
!
	  HARD=.TRUE.
!
1200	  CALL NEW_GEN_IN(PLT_LINE_WGT,'line weight')
	  IF (PLT_LINE_WGT .GT. 201) THEN
	    WRITE(T_OUT,*)'value too large'
	    GOTO 1200
	  END IF
	  GOTO 5000
!
	ELSE IF(ANS .EQ. 'CL')THEN
	  IF(.NOT. HARD)THEN			!Terminal
	    CALL PGERAS
	    CALL PGETXT
	    CALL PGPAGE
	    READ(T_IN,'(A)')ANS
	  END IF
	  GOTO 1000
!
	ELSE IF(ANS .EQ. 'WP' .OR. ANS .EQ. 'WPF')THEN
	  IF(PLT_ST_FILENAME .EQ. ' ')THEN
	    PLT_ST_FILENAME='PLT_ST'
	    CALL NEW_GEN_IN(PLT_ST_FILENAME,'Name of file with stored plots')
	  ELSE IF(ANS .EQ. 'WPF')THEN
	    CALL NEW_GEN_IN(PLT_ST_FILENAME,'Name of file with stored plots')
	  END IF
	  OPEN(UNIT=30,FORM='UNFORMATTED',FILE=PLT_ST_FILENAME,STATUS='UNKNOWN'
	1,       IOSTAT=IOS,POSITION='APPEND')
	  IF(IOS .NE. 0)THEN
	    WRITE(T_OUT,*)'Error opening file'
	    GOTO 1000
	  END IF
	  L=-1
	  IF(NPLTS .NE. 1)CALL NEW_GEN_IN(L,'Plot to output (def=-1=ALL)')
!
	  CNT=0; PLT_ID=' '
	  DO IP=1,NPLTS
	    IF(L .EQ. -1 .OR. L .EQ. IP)THEN
	      CNT=CNT+1
	      CALL NEW_GEN_IN(PLT_ID,'PLT_ID=')
	      WRITE(30)'PLT_ID=',PLT_ID
	      WRITE(30)NPTS(IP)
	      WRITE(30)ERR(IP)
	      DO K=1,(NPTS(IP)+N_REC_SIZE-1)/N_REC_SIZE
	        IST=N_REC_SIZE*(K-1)+1
	        IEND=MIN(IST+N_REC_SIZE-1,NPTS(IP))
	        WRITE(30)(CD(IP)%XVEC(I),I=IST,IEND)
	        WRITE(30)(CD(IP)%DATA(I),I=IST,IEND)
	        IF(ERR(IP))THEN
	          WRITE(30)(CD(IP)%EMAX(I),I=IST,IEND)
	          WRITE(30)(CD(IP)%EMIN(I),I=IST,IEND)
	        END IF
	      END DO
	    END IF
	  END DO
	  WRITE(T_OUT,'(I4,A)')CNT,' plots written to data file'
	  CLOSE(UNIT=30)
	  GOTO 1000
!
	ELSE IF(ANS .EQ. 'RP' .OR. ANS .EQ. 'RPF')THEN
	  IF(PLT_ST_FILENAME .EQ. ' ')THEN
	    PLT_ST_FILENAME='PLT_ST'
	    CALL NEW_GEN_IN(PLT_ST_FILENAME,'Name of file with stored plots')
	  ELSE IF(ANS .EQ. 'RPF')THEN
	    CALL NEW_GEN_IN(PLT_ST_FILENAME,'Name of file with stored plots')
	  END IF
	  OPEN(UNIT=30,FORM='UNFORMATTED',FILE=PLT_ST_FILENAME,STATUS='OLD',IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(T_OUT,*)'Error opening file'
	    GOTO 1000
	  END IF
	  PLT_ID=' '
	  PLT_ID_SAV=' '
!
	  DO WHILE(1 .EQ. 1)
!
	    WRITE(6,*)' '
	    WRITE(6,*)'The identifier must be unique, but need not be complete.'
	    WRITE(6,*)' '
	    CALL NEW_GEN_IN(PLT_ID,'PLT_ID=')
	    IF(PLT_ID .EQ. ' ' .OR. PLT_ID .EQ. PLT_ID_SAV)EXIT
	    RD_PLT_ID=' '
!
	    PLT_ID='PLT_ID='//ADJUSTL(PLT_ID)
	    I=LEN_TRIM(PLT_ID)
	    DO WHILE(RD_PLT_ID(1:I) .NE. PLT_ID(1:I))
	      READ(30,IOSTAT=IOS)RD_PLT_ID
	      IF(IOS .EQ. -1)THEN
	        WRITE(T_OUT,*)'Unable to identify plot --- IOS',IOS
	        WRITE(T_OUT,*)'Available plots are as follows:'
	        REWIND(UNIT=30)
	        IOS=0
                DO WHILE(IOS .NE. -1)
	          RD_PLT_ID=' '
	          READ(30,IOSTAT=IOS)RD_PLT_ID
	          IF(RD_PLT_ID(1:7) .EQ. 'PLT_ID=')THEN
	             WRITE(T_OUT,'(35X,A)')TRIM(RD_PLT_ID(8:))
	          END IF
	        END DO
	        CLOSE(UNIT=30)
	        GOTO 1000
	      END IF
	      I=LEN_TRIM(PLT_ID)
	    END DO
	    IP=NPLTS+1
	    IF(IP .GT. MAX_PLTS)THEN
	      WRITE(T_OUT,*)'Error in GRAMON_PGPLOT'
	      WRITE(T_OUT,*)'Too many plots have ben stored'
	      GOTO 1000
	    END IF
	    READ(30)NPTS(IP)
	    READ(30)ERR(IP)
	    ALLOCATE (CD(IP)%XVEC(NPTS(IP)))
	    ALLOCATE (CD(IP)%DATA(NPTS(IP)))
	    IF(ERR(IP))THEN
	      ALLOCATE (CD(IP)%EMIN(NPTS(IP)))
	      ALLOCATE (CD(IP)%EMAX(NPTS(IP)))
	    END IF	
	    DO K=1,(NPTS(IP)+N_REC_SIZE-1)/N_REC_SIZE
	      IST=N_REC_SIZE*(K-1)+1
	      IEND=MIN(IST+N_REC_SIZE-1,NPTS(IP))
	      READ(30)(CD(IP)%XVEC(I),I=IST,IEND)
	      READ(30)(CD(IP)%DATA(I),I=IST,IEND)
	      IF(ERR(IP))THEN
	        READ(30)(CD(IP)%EMAX(I),I=IST,IEND)
	        READ(30)(CD(IP)%EMIN(I),I=IST,IEND)
	      END IF
	    END DO
	    NPLTS=NPLTS+1			!Successfully read.
	    TYPE_CURVE(IP)='L'
	    LINE_WGT(IP)=1
	    IF(IP .NE. 1)THEN
	      LINE_STYLE(IP)=MOD(LINE_STYLE(IP-1),5)+1
              PEN_COL(IP)=PEN_COL(IP-1)+1
            ELSE
	      LINE_STYLE(IP)=1
	      PEN_COL(IP)=2
	    END IF
	    LINE_STYLE(IP)=1
	    IF(DASH)LINE_STYLE(IP)=MOD(IP-1,5)+1
	    PLT_ID_SAV=PLT_ID
	    CD(IP)%CURVE_ID=PLT_ID(8:)
	    CALL NEW_GEN_IN(CD(IP)%CURVE_ID,'Curve label')
	    REWIND(UNIT=30)
	  END DO
!
	  CLOSE(UNIT=30)
	  GOTO 1000
!
!
! Option to read simple plots in column format from file.
!
	ELSE IF(ANS .EQ. 'RXY')THEN
	  FILNAME=' '
	  CALL NEW_GEN_IN(FILNAME,'FILE=')
	  OPEN(UNIT=30,FILE=FILNAME,STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(T_OUT,*)'Error opening file'
	    GOTO 1000
	  END IF
!
! Use FILNAME as temporary STRING
!
	  CNT=0
	  FILNAME='!'
	  DO WHILE(FILNAME(1:1) .EQ. '!' .OR. FILNAME .EQ. ' ')
	    READ(30,'(A200)',IOSTAT=IOS)FILNAME
	    WRITE(T_OUT,*)TRIM(FILNAME)
	    IF(IOS .NE. 0)THEN
	      WRITE(T_OUT,*)'Error reading file'
	      CLOSE(UNIT=30)
	      GOTO 1000
	    END IF
	    IF(INDEX(FILNAME,'Number of data points:') .NE. 0)THEN
	      I=INDEX(FILNAME,'Number of data points:')
	      READ(FILNAME(I+22:),*)CNT
	      FILNAME='!'
	    END IF
	  END DO
	  BACKSPACE(UNIT=30)
!
	  IP=NPLTS+1
	  NPTS(IP)=CNT
	  IF(CNT .EQ. 0)THEN
	    CALL NEW_GEN_IN(NPTS(IP),'Number of data points')
	    FILNAME='Skip this record: '//TRIM(FILNAME)
	    CALL NEW_GEN_IN(TMP_LOG,FILNAME)
	    IF(TMP_LOG)READ(30,'(A200)',IOSTAT=IOS)FILNAME
	  END IF
	  COLUMN(1)=1; COLUMN(2)=2
	  CALL NEW_GEN_IN(COLUMN,I,ITWO,'Data columns')
	  K=MAX(COLUMN(1),COLUMN(2))
!
	  ALLOCATE (CD(IP)%XVEC(NPTS(IP)),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (CD(IP)%DATA(NPTS(IP)),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	     WRITE(T_OUT,*)'Error allocating data arrays'
	     GOTO 1000
	  END IF
	  DO I=1,NPTS(IP)
	    READ(30,*,IOSTAT=IOS,END=500)(TEMP_VAR(L),L=1,K)
	    CD(IP)%XVEC(I)=TEMP_VAR(COLUMN(1))
	    CD(IP)%DATA(I)=TEMP_VAR(COLUMN(2))
	  END DO
	  CD(IP)%CURVE_ID=' '
	  CALL NEW_GEN_IN(CD(IP)%CURVE_ID,'Curve label')
	  IF(IOS .NE. 0)THEN
	    WRITE(T_OUT,*)'Error in read plot from XY file'
	    WRITE(T_OUT,*)'Data point is',I+1
	    CLOSE(UNIT=30)
	    GOTO 1000
	  END IF
500	  CONTINUE
!
! Data successfully read
!	
	  LAST_DP=LAST_DP+NPTS(IP)
	  NPLTS=NPLTS+1; J=NPLTS
	  TYPE_CURVE(J)='L'
	  LINE_WGT(J)=1
	  IF(J .NE. 1)THEN
	    LINE_STYLE(J)=MOD(LINE_STYLE(J-1),5)+1
            PEN_COL(J)=PEN_COL(J-1)+1
          ELSE
	    LINE_STYLE(J)=1
	    PEN_COL(J)=2
	  END IF
	  CLOSE(UNIT=30)
	  GOTO 1000
!
	ELSE IF(ANS .EQ. 'WXY')THEN
	  FILNAME=' '
	  CALL NEW_GEN_IN(FILNAME,'FILE=')
	  OPEN(UNIT=30,FILE=FILNAME,STATUS='UNKNOWN',ACTION='WRITE',IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(T_OUT,*)'Error opening file'
	    GOTO 1000
	  END IF
	  LAM_ST=0.0D0; LAM_END=0.0D0
	  CALL NEW_GEN_IN(LAM_ST,'Start wavelength for output (def=ALL)')
	  CALL NEW_GEN_IN(LAM_END,'End waveliength for output (def=ALL)')
!
	  IPLT=-1
	  CALL NEW_GEN_IN(IPLT,NP_OUT,NPLTS,'Plots to output (def=-1=ALL)')
	  IF(IPLT(1) .EQ. -1)THEN
	    DO I=1,NPLTS; IPLT(I)=I; END DO
	    NP_OUT=NPLTS
	  END IF
	  DO I=1,NP_OUT
	    IF(IPLT(I) .LT. 0 .OR. IPLT(I) .GT. NPLTS)THEN
	      WRITE(6,*)'Error invalid plot -- ',IPLT(I)
	      WRITE(6,*)'Valid plot range is 1:',NPLTS
	      GOTO 1000
	    END IF
	  END DO
!
	  WRITE(30,*)NP_OUT
	  NX_MAX=0.0D0
	  DO I=1,NP_OUT
	    NX_MAX=MAX(NX_MAX,NPTS(IPLT(I)))
	  END DO
!
! Output data for all plots. Plots are padded with zero.
!
	  IF(LAM_ST .EQ. LAM_END)THEN
	    WRITE(30,*)(NPTS(IPLT(I)),I=1,NP_OUT)
	    DO J=1,NX_MAX
	      DO I=1,NP_OUT
	        IP=IPLT(I)
	        ADVANCE_OPT='NO'
	        IF(I .EQ. NP_OUT)ADVANCE_OPT='YES'
	        IF(J .LE. NPTS(IP))THEN
	          WRITE(30,'(2X,ES14.7,ES14.6)',ADVANCE=ADVANCE_OPT)CD(IP)%XVEC(J),CD(IP)%DATA(J)
	        ELSE
	          WRITE(30,'(2X,ES14.7,ES14.6)',ADVANCE=ADVANCE_OPT)0.0D0,0.0D0
	        END IF
	      END DO
	    END DO
	  ELSE
	    DO I=1,NP_OUT
	      IP=IPLT(I)
	      ADVANCE_OPT='NO'
	      IF(I .EQ. NP_OUT)ADVANCE_OPT='YES'
	      WRITE(30,'(15X,A15)',ADVANCE=ADVANCE_OPT)TRIM(CD(IP)%CURVE_ID)
	    END DO
!
	    OUTPUT_STRING=' '
	    IPST=0; NX_MAX=0
	    DO I=1,NP_OUT
	      CNT=0
	      IP=IPLT(I)
	      DO J=1,NPTS(I)
	        IF((LAM_ST-CD(IP)%XVEC(J))*(CD(IP)%XVEC(J)-LAM_END) .GE. 0)THEN
	          CNT=CNT+1
	          IF(IPST(I) .EQ. 0)IPST(I)=J
	        END IF
	      END DO
	      NX_MAX=MAX(NX_MAX,CNT)
	      ADVANCE_OPT='NO'
	      IF(I .EQ. NP_OUT)ADVANCE_OPT='YES'
	      K=LEN_TRIM(OUTPUT_STRING)
	      WRITE(30,'(20X,I10)',ADVANCE=ADVANCE_OPT)CNT
	    END DO
!
	    DO L=1,NX_MAX
	      OUTPUT_STRING=' ';  CNT=0
	      DO I=1,NP_OUT
	        J=IPST(I)+L-1
	        IP=IPLT(I)
	        K=LEN_TRIM(OUTPUT_STRING)
	        ADVANCE_OPT='NO'
	        IF(I .EQ. NP_OUT)ADVANCE_OPT='YES'
	        IF(J .LE. NPTS(IP))THEN
	          IF( (LAM_ST-CD(IP)%XVEC(J))*(CD(IP)%XVEC(J)-LAM_END) .GE. 0)THEN
	            WRITE(30,'(2X,ES14.7,ES14.6)',ADVANCE=ADVANCE_OPT)CD(IP)%XVEC(J),CD(IP)%DATA(J)
	            CNT=CNT+1
	          END IF
	        ELSE IF(NPLTS .NE. 1)THEN
	          WRITE(30,'(2X,ES14.7,ES14.6)',ADVANCE=ADVANCE_OPT)0.0D0,0.0D0
	        END IF
	      END DO
	    END DO
	  END IF
	  WRITE(T_OUT,*)NP_OUT,' plots written to ',TRIM(FILNAME)
	  CLOSE(UNIT=30)
!
!  Writes data to a file with same format as OBSFLUX. Option assumes
!  units are 10^15 Hz and Jy at 1 kps.
!
	ELSE IF(ANS .EQ. 'WOBS')THEN
	  FILNAME=' '
	  CALL NEW_GEN_IN(FILNAME,'FILE=')
	  OPEN(UNIT=30,FILE=FILNAME,STATUS='UNKNOWN',ACTION='WRITE',IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(T_OUT,*)'Error opening file'
	    GOTO 1000
	  END IF
	  J=30
	  IP=1; CALL NEW_GEN_IN(IP,'Which plot to output?')
	  WRITE(WK_STR,'(I10)')NPTS(IP)
	  WK_STR=ADJUSTL(WK_STR)
	  WRITE(30,'(/,A,/)')' Continuum Frequencies ( '//TRIM(WK_STR)//' )'
	  WRITE(30,'(8ES15.7)')(CD(IP)%XVEC(I),I=1,NPTS(IP))	
	  WRITE(30,'(//,A,/)')' Observed intensity (Janskys)'
	  WRITE(30,'(10ES12.4)')(CD(IP)%DATA(I),I=1,NPTS(IP))	
	  CLOSE(UNIT=30)
!
! 
!
	ELSE IF(ANS .EQ. 'SXY')THEN
	  DO IP=1,NPLTS
	    DO J=1,NPTS(IP)
	      IF(CD(IP)%XVEC(J) .GE. XPAR(1) .AND. CD(IP)%XVEC(J) .LE. XPAR(2))THEN
	        WRITE(T_OUT,*)IP,J,CD(IP)%XVEC(J),CD(IP)%DATA(J)
	      END IF
	    END DO
	  END DO
! 
!
	ELSE IF(ANS .EQ. 'LAM')THEN
	  LAM_OPTION=' '
	  CALL NEW_GEN_IN(LAM_OPTION,'Species -- none, He2, HE')
	  CALL SET_CASE_UP(LAM_OPTION,1,0)
	  WRITE(6,*)RED_PEN
	  WRITE(6,'(/,A,/)')' All wavelengths are vacuum'//BLUE_PEN
	  CALL WRITE_LINE_LAMBDAS(LAM_OPTION)
	  WRITE(6,*)DEF_PEN
	  GOTO 1000
!
	ELSE IF (ANS .EQ. 'CUT')THEN          !Recall ANS in always upper case
	  WRITE(6,'(A)')BLUE_PEN
	  WRITE(6,'(A)')' This option reduces the number of data points in each curve.'
	  WRITE(6,'(A)')' It is useful for creating smaller publication quality plots.'//RED_PEN
	  WRITE(6,'(A)')' This option cannot be undone.'//BLUE_PEN
	  WRITE(6,'(A)')' A negative cut accuracy does nothing'
	  WRITE(6,'(A)')' To compare with original data:'
	  WRITE(6,'(A)')'    Store original data using WP & read after cut with RP, OR'
	  WRITE(6,'(A)')'    do CUT, NOI, and redo plots'
	  WRITE(6,'(A)')DEF_PEN
	  CUT_ACCURACY=0.001D0
	  CALL NEW_GEN_IN(CUT_ACCURACY,'Fractional accuracy to retain in plots')
	  IF(CUT_ACCURACY .GT. 0.0)CALL SHRINK_VECTORS(CUT_ACCURACY)
	  GOTO 1000
!
	ELSE IF (ANS .EQ. 'KMS')THEN          !Recall ANS in always upper case
	  VEL_UNIT='km/s'
	  WRITE(6,*)'Using km/s for velocity unit'
	  GOTO 1000
!
	ELSE IF (ANS .EQ. 'MMS')THEN
	  VEL_UNIT='Mm/s'
	  WRITE(6,*)'Using Mm/s for velocity unit'
	  GOTO 1000
!
! Convert from Ang to velocity space. Data must have been originally
! in Ang. This option can be done many times, as old Ang scale is
! restored on each call.
!
	ELSE IF (ANS .EQ. 'VEL')THEN
	  C_KMS=1.0D-05*SPEED_OF_LIGHT()
	  C_VAL=C_KMS
	  IF(VEL_UNIT .EQ. 'Mm/s')C_VAL=1.0D-03*C_KMS
	  CALL NEW_GEN_IN(CENTRAL_LAM,'/\(Ang) [-ve: 10^15 Hz]')
	  IF(CENTRAL_LAM .LT. 0)THEN
	    CENTRAL_LAM=1.0D-02*C_VAL/ABS(CENTRAL_LAM)
	  ELSE IF(CENTRAL_LAM .GT. 2000)THEN
	    AIR_LAM=.FALSE.
            CALL NEW_GEN_IN(AIR_LAM,'Wavelength in air?')
	    IF(AIR_LAM)CENTRAL_LAM=LAM_VAC(CENTRAL_LAM)
	  END IF
!
! Undo previous conversion to V space.
!
	  IF(OLD_CENTRAL_LAM .NE. 0)THEN
	    DO IP=1,NPLTS
	      DO J=1,NPTS(IP)
	        CD(IP)%XVEC(J)=OLD_CENTRAL_LAM*
	1                         (1.0+CD(IP)%XVEC(J)/C_VAL)
	      END DO
	    END DO
	  END IF
!
! Puts X-axis in km/s
!
	  IF(CENTRAL_LAM .NE. 0)THEN
	    T1=C_VAL/CENTRAL_LAM
	    DO IP=1,NPLTS
	      DO J=1,NPTS(IP)
	        CD(IP)%XVEC(J)=T1*(CD(IP)%XVEC(J)-CENTRAL_LAM)
	      END DO
	    END DO
	    XLABEL='V(km\d \us\u-1\d)'
	    IF(VEL_UNIT .EQ. 'Mm/s')XLABEL='V(Mm\d \us\u-1\d)'
	  ELSE
	    XLABEL=XLAB
	  END IF
	  OLD_CENTRAL_LAM=CENTRAL_LAM
	  GOTO 1000
! 
!
! Peform simple X-axis arithmetic.
!
	ELSE IF (ANS .EQ. 'XAR')THEN
	  CALL NEW_GEN_IN(XAR_OPERATION,'Operation: *,+,-,/,LG,ALG[=10^x],R[=1/x],XN')
	  CALL SET_CASE_UP(XAR_OPERATION,IZERO,IZERO)
	  IF(XAR_OPERATION .NE. 'LG' .AND. XAR_OPERATION .NE. 'ALG' .AND.
	1           XAR_OPERATION .NE. 'XN' .AND. XAR_OPERATION .NE. 'R')THEN
	    CALL NEW_GEN_IN(XAR_VAL,'Value')
	  END IF
	  CALL NEW_GEN_IN(XAR_PLT,'Which plot? (0=ALL,-ve exits)')
	  IP=XAR_PLT
	  IF(IP .LT. 0 .OR. IP .GT. NPLTS)THEN
	    WRITE(T_OUT,*)'Invalid plot number'
	    GOTO 1000
	  END IF
	  IP_ST=IP; IF(IP .EQ. 0)IP_ST=1
	  IP_END=IP; IF(IP .EQ. 0)IP_END=NPLTS
	  DO IP=IP_ST,IP_END
	    IF(XAR_OPERATION .EQ. '+')THEN
	      CD(IP)%XVEC=CD(IP)%XVEC+XAR_VAL
	    ELSE IF(XAR_OPERATION .EQ. '-')THEN
	      CD(IP)%XVEC=CD(IP)%XVEC-XAR_VAL
	    ELSE IF(XAR_OPERATION .EQ. '*')THEN
	      CD(IP)%XVEC=CD(IP)%XVEC*XAR_VAL
	    ELSE IF(XAR_OPERATION .EQ. '/')THEN
	      CD(IP)%XVEC=CD(IP)%XVEC/XAR_VAL
	    ELSE IF(XAR_OPERATION .EQ. 'ALG')THEN
	      CD(IP)%XVEC=10.0D0**(CD(IP)%XVEC)
	      IF(XLABEL(1:3) .EQ. 'Log')THEN
	        XLABEL=XLABEL(4:)
	        XLABEL=ADJUSTL(XLABEL)
	      END IF
	    ELSE IF(XAR_OPERATION .EQ. 'SSQR')THEN
	      DO I=1,NPTS(IP)
	        T1=SQRT(ABS(CD(IP)%XVEC(I)))
	        CD(IP)%XVEC(I)=SIGN(T1,CD(IP)%XVEC(I))
	      END DO
	      XLABEL=TRIM(XLABEL)//'/SQRT(|'//TRIM(XLABEL)//'|)'
	    ELSE IF(XAR_OPERATION .EQ. 'LG')THEN
	      T1=TINY(CD(IP)%XVEC(1))
	      DO J=1,NPTS(IP)
	        IF(CD(IP)%XVEC(J) .GT. T1)THEN
	          CD(IP)%XVEC(J)=LOG10(CD(IP)%XVEC(J))
	        ELSE
	          CD(IP)%XVEC(J)=-1000.0
	        END IF
	      END DO
	      IF(XLABEL(1:3) .NE. 'Log')THEN
	        XLABEL='Log '//XLABEL
	      END IF
!
	    ELSE IF(XAR_OPERATION .EQ. 'R')THEN
	      T1=0.01D0*C_KMS
	      T2=0.01D0*C_KMS
	      CALL NEW_GEN_IN(T1,'Scale factor')
	      DO J=1,NPTS(IP)
	        IF(CD(IP)%XVEC(J) .GT. 0)THEN
	          CD(IP)%XVEC(J)=T1/CD(IP)%XVEC(J)
	        ELSE
	          CD(IP)%XVEC(J)=-1000
	  	END IF
	      END DO
	      IF(XLABEL .EQ. '\gn(10\u15 \dHz)' .AND. T1 .EQ. T2)THEN
	        XLABEL='\gl(\A)'
	      ELSE IF(XLABEL .EQ. '\gl(\A)' .AND. T1 .EQ. T2)THEN
	        XLABEL='\gn(10\u15 \dHz)'
	      END IF
!
! Set X axis to integer counter.
!
	    ELSE IF(XAR_OPERATION .EQ. 'XN')THEN
	      WRITE(6,*)RED_PEN
	      WRITE(6,*)'Curent plot is',IP
	      K=0; TMP_STR='Output plot: 0 adds new plot?'//DEF_PEN
	      CALL NEW_GEN_IN(K,TRIM(TMP_STR))
	      IF(K .GT. NPLTS .OR. K .LT. 0)THEN
	         WRITE(6,*)'Invlid plot number - maximum number of plots is',NPLTS
	         GOTO 1000
	      END IF
	      IF(K .EQ. IP)THEN
	        K=IP
	      ELSE IF(K .EQ. 0)THEN
	        CALL CURVE(NPTS(IP),CD(IP)%XVEC,CD(IP)%DATA)
	        K=NPLTS
	        WRITE(6,'(/,A,I4,/)')BLUE_PEN//' New plot number is '//DEF_PEN,K
	      ELSE
	        IF(ALLOCATED(CD(IP)%XVEC))DEALLOCATE(CD(IP)%XVEC,CD(IP)%DATA)
	        ALLOCATE(CD(IP)%XVEC(NPTS(IP)),CD(IP)%DATA(NPTS(IP)))
	        NPTS(K)=NPTS(IP)
	        NPLTS=NPLTS+1
	        CD(K)%DATA=CD(I)%DATA
	      END IF
	      DO I=1,NPTS(K)
	        CD(K)%XVEC(I)=I
	      END DO
	      TYPE_CURVE(K)=TYPE_CURVE(IP)
	    ELSE
	      WRITE(T_OUT,*)'Invalid operation: try again'
	      GOTO 1000
	    END IF
	  END DO
	  GOTO 1000
!
	ELSE IF (ANS .EQ. 'OFF')THEN
	  CALL NEW_GEN_IN(YAR_PLT,'Which plot? (0=ALL,-ve exits)')
	  IP=YAR_PLT
	  IP_ST=IP; IF(IP .EQ. 0)IP_ST=1
	  IP_END=IP; IF(IP .EQ. 0)IP_END=NPLTS
	  IF(IP .LT. 0 .OR. IP .GT. NPLTS)THEN
	    WRITE(T_OUT,*)'Invalid plot number'
	    GOTO 1000
	  END IF
	  CALL NEW_GEN_IN(YAR_VAL,'Offset')
	  DO IP=IP_ST,IP_END
	    CD(IP)%DATA=CD(IP)%DATA+YAR_VAL*(IP-IP_ST)
	  END DO
!
! Perform simple Y-axis arithmetic.
!
	ELSE IF (ANS .EQ. 'YAR')THEN
	  CALL NEW_GEN_IN(YAR_OPERATION,'Operation: *,+,-,/,LG,R[=1/Y],ALG[=10^y],ABS,MX,DX,+X,-X')
	  CALL SET_CASE_UP(YAR_OPERATION,IZERO,IZERO)
	  IF(YAR_OPERATION .EQ. 'MX')YAR_OPERATION='*X'
	  IF(YAR_OPERATION .EQ. 'DX')YAR_OPERATION='/X'
	  IF(YAR_OPERATION .EQ. '*X' .OR. YAR_OPERATION .EQ. '/X')THEN
	    YAR_VAL=0.5D0*(XPAR(1)+XPAR(2))
	    CALL NEW_GEN_IN(YAR_VAL,'X normalization')
	  ELSE IF( YAR_OPERATION .NE. 'LG' .AND.
	1          YAR_OPERATION .NE. 'ABS' .AND.
	1          YAR_OPERATION .NE. 'ALG' .AND.
	1          YAR_OPERATION(1:1) .NE. 'R')THEN
	    CALL NEW_GEN_IN(YAR_VAL,'Value')
	  END IF
	  CALL NEW_GEN_IN(YAR_PLT,'Which plot? (0=ALL,-ve exits)')
	  IP=YAR_PLT
	  IP_ST=IP; IF(IP .EQ. 0)IP_ST=1
	  IP_END=IP; IF(IP .EQ. 0)IP_END=NPLTS
	  IF(IP .LT. 0 .OR. IP .GT. NPLTS)THEN
	    WRITE(T_OUT,*)'Invalid plot number'
	    GOTO 1000
	  END IF
	  DO IP=IP_ST,IP_END
	    IF(YAR_OPERATION .EQ. '+')THEN
	      CD(IP)%DATA=CD(IP)%DATA+YAR_VAL
	      IF(ERR(IP))CD(IP)%EMIN=CD(IP)%EMIN+YAR_VAL
	      IF(ERR(IP))CD(IP)%EMAX=CD(IP)%EMAX+YAR_VAL
	    ELSE IF(YAR_OPERATION .EQ. '-')THEN
	      CD(IP)%DATA=CD(IP)%DATA-YAR_VAL
	      IF(ERR(IP))CD(IP)%EMIN=CD(IP)%EMIN-YAR_VAL
	      IF(ERR(IP))CD(IP)%EMAX=CD(IP)%EMAX-YAR_VAL
	    ELSE IF(YAR_OPERATION .EQ. '*')THEN
	      CD(IP)%DATA=CD(IP)%DATA*YAR_VAL
	      IF(ERR(IP))CD(IP)%EMIN=CD(IP)%EMIN*YAR_VAL
	      IF(ERR(IP))CD(IP)%EMAX=CD(IP)%EMAX*YAR_VAL
	    ELSE IF(YAR_OPERATION .EQ. '/')THEN
	      CD(IP)%DATA=CD(IP)%DATA/YAR_VAL
	      IF(ERR(IP))CD(IP)%EMIN=CD(IP)%EMIN/YAR_VAL
	      IF(ERR(IP))CD(IP)%EMAX=CD(IP)%EMAX/YAR_VAL
	    ELSE IF(YAR_OPERATION .EQ. 'ABS')THEN
	      CD(IP)%DATA=ABS(CD(IP)%DATA)
	      WRITE(6,*)'Warning -- error vector not adjusted'
	    ELSE IF(YAR_OPERATION .EQ. 'ALG')THEN
	      CD(IP)%DATA=10.0D0**(CD(IP)%DATA)
	      IF(YLABEL(1:3) .EQ. 'Log')THEN
	        YLABEL=YLABEL(4:)
	        YLABEL=ADJUSTL(YLABEL)
	      END IF
	    ELSE IF(YAR_OPERATION .EQ. 'LG')THEN
	      T1=TINY(CD(IP)%DATA(1))
	      DO J=1,NPTS(IP)
	        IF(CD(IP)%DATA(J) .GT. T1)THEN
	          CD(IP)%DATA(J)=LOG10(CD(IP)%DATA(J))
	        ELSE
	          CD(IP)%DATA(J)=-1000.0
	        END IF
	      END DO
	      IF(YLABEL(1:3) .NE. 'Log')THEN
	        YLABEL='Log '//YLABEL
	      END IF
	    ELSE IF(YAR_OPERATION .EQ. 'R')THEN
	      DO J=1,NPTS(IP)
	        IF(ABS(CD(IP)%DATA(J)) .GT. 1.0E-37)THEN
	          CD(IP)%DATA(J)=1.0D0/CD(IP)%DATA(J)
	        ELSE
	          CD(IP)%DATA(J)=1.0E+37
	        END IF
	      END DO
	    ELSE IF(YAR_OPERATION .EQ. '+X')THEN
	      DO J=1,NPTS(IP)
	       CD(IP)%DATA(J)=CD(IP)%DATA(J)+CD(IP)%XVEC(J)
	      END DO
	    ELSE IF(YAR_OPERATION .EQ. '-X')THEN
	      DO J=1,NPTS(IP)
	       CD(IP)%DATA(J)=CD(IP)%DATA(J)-CD(IP)%XVEC(J)
	      END DO
	    ELSE IF(YAR_OPERATION .EQ. '*X')THEN
	      DO J=1,NPTS(IP)
	       CD(IP)%DATA(J)=CD(IP)%DATA(J)*(CD(IP)%XVEC(J)/YAR_VAL)
	      END DO
	    ELSE IF(YAR_OPERATION .EQ. '/X')THEN
	      DO J=1,NPTS(IP)
	        IF(CD(IP)%XVEC(J) .NE. 0.0D0)THEN
	          CD(IP)%DATA(J)=CD(IP)%DATA(J)/(CD(IP)%XVEC(J)/YAR_VAL)
	        ELSE
	          CD(IP)%DATA(J)=1.0E+37
	        END IF
	      END DO
	    ELSE
	      WRITE(T_OUT,*)'Invalid operation: try again'
	      GOTO 1000
	    END IF
	  END DO
	  GOTO 1000
!
	ELSE IF (ANS .EQ. 'CAX')THEN
	  CALL NEW_GEN_IN(VAR_PLT1,'Plot that is to be altered')
	  CALL NEW_GEN_IN(VAR_OPERATION,'Operation: RX, RY, RXY, RYX, SXY')
	  CALL SET_CASE_UP(VAR_OPERATION,IZERO,IZERO)
	  VAR_PLT2=VAR_PLT1
	  IF(VAR_OPERATION .NE. 'SXY')CALL NEW_GEN_IN(VAR_PLT2,'Plot with new data')
	  IF(NPTS(VAR_PLT1) .NE. NPTS(VAR_PLT2))THEN
	  ELSE IF(VAR_OPERATION .EQ. 'RX')THEN
	     CD(VAR_PLT1)%XVEC=CD(VAR_PLT2)%XVEC
	  ELSE IF(VAR_OPERATION .EQ. 'RY')THEN
	     CD(VAR_PLT1)%DATA=CD(VAR_PLT2)%DATA
	  ELSE IF(VAR_OPERATION .EQ. 'RYWX')THEN
	     CD(VAR_PLT1)%DATA=CD(VAR_PLT2)%XVEC
	  ELSE IF(VAR_OPERATION .EQ. 'RXWY')THEN
	     CD(VAR_PLT1)%XVEC=CD(VAR_PLT2)%DATA
	  ELSE IF(VAR_OPERATION .EQ. 'SXY')THEN
	     DO I=1,NPTS(VAR_PLT1)
	       T1=CD(VAR_PLT1)%DATA(I)
	       CD(VAR_PLT1)%DATA(I)=CD(VAR_PLT1)%XVEC(I)
	       CD(VAR_PLT1)%XVEC(I)=T1
	     END DO
	  END IF
!
	ELSE IF (ANS .EQ. 'REG')THEN
!
	  WRITE(6,'(A)')RED_PEN
	  WRITE(6,'(A)')'Option for creating a new grid'
	  WRITE(6,'(A)')'Plot window is used for defining X-regridding range'
	  WRITE(6,'(A)')DEF_PEN
!
	  CALL NEW_GEN_IN(VAR_PLT1,'Input plot 1?')
	  VAR_PLT3=NPLTS+1
	  CALL NEW_GEN_IN(VAR_PLT3,'Output plot?')
	  TYPE_CURVE(VAR_PLT3)='L'
	  T1=0.0; T2=0.0
	  CALL NEW_GEN_IN(REG_OPT,'Regrid option: UG, dX, R, NINS')
	  CALL SET_CASE_UP(REG_OPT,1,0)
	  IF(REG_OPT .NE. 'UG')CALL NEW_GEN_IN(T1,'dX, R, or NINS (zero to quit)')
	  CALL DO_PG_REGRID(VAR_PLT1,VAR_PLT3,XPAR(1),XPAR(2),REG_OPT,T1)
!
	ELSE IF(ANS .EQ. 'ADDN')THEN
	  CALL NEW_GEN_IN(VAR_PLT1,'Input plot 1?')
	  VAR_PLT3=NPLTS+1
	  CALL NEW_GEN_IN(VAR_PLT3,'Output plot?')
	  T1=1000.0; CALL NEW_GEN_IN(T1,'Counts in continnum')
	  T2=0.2;    CALL NEW_GEN_IN(T2,'Random number 0 to 1)')
	  IF(VAR_PLT3 .NE. VAR_PLT1)THEN
	    TYPE_CURVE(VAR_PLT3)='L'
	    CD(VAR_PLT3)%CURVE_ID=CD(VAR_PLT1)%CURVE_ID
	  END IF
	  CALL PG_ADD_NOISE(VAR_PLT1,VAR_PLT3,XMIN,XMAX,T1,T2)
!	
	ELSE IF (ANS .EQ. 'VAR')THEN
	  CALL NEW_GEN_IN(VAR_PLT1,'Input plot 1?')
	  CALL NEW_GEN_IN(VAR_OPERATION,'Operation: *,+,-,/,c(opy),cc')
	  CALL SET_CASE_UP(VAR_OPERATION,IZERO,IZERO)
	  IF(VAR_OPERATION(1:2) .EQ. 'c ' .OR. VAR_OPERATION(1:2) .EQ. 'C ')THEN
	    VAR_OPERATION='C'
	    VAR_PLT2=VAR_PLT1
	  ELSE
	    CALL NEW_GEN_IN(VAR_PLT2,'Input plot 2?')
	  END IF
	  VAR_PLT3=NPLTS+1
	  IF(VAR_OPERATION(1:2) .EQ. 'cc' .OR. VAR_OPERATION(1:2) .EQ. 'CC')THEN
	    WRITE(6,*)RED_PEN
	    WRITE(6,*)'For cross-correlation of spectra plot 1 should be on a uniform grid'
	    WRITE(6,*)'and in log space. The second plot should also be in log space'
	    WRITE(6,*)'Use REG option to create a uniform grid (UG option)'
	    WRITE(6,*)'Enter 0 for output plot to exit this option'
	    WRITE(6,*)' '
	    WRITE(6,*)'Plot 1 should be the reference plot'
	    WRITE(6,*)'Use GF to measure location of peak'
	    WRITE(6,*)'Result is the required radial velocity'
	    WRITE(6,*)DEF_PEN
	    TMP_LOG=.TRUE.
	    CALL NEW_GEN_IN(TMP_LOG,'Scale X axis of output (only) to km/s')
	  END IF
	  CALL NEW_GEN_IN(VAR_PLT3,'Output plot?')
	  IF(VAR_PLT3 .EQ. 0)GOTO 1000
	  TYPE_CURVE(VAR_PLT3)='L'
	  CALL DO_VEC_OP(VAR_PLT1,VAR_PLT2,VAR_PLT3,.TRUE.,VAR_OPERATION)
	  IF(TMP_LOG)THEN
	   K=NPTS(VAR_PLT3)
	   T1=C_KMS*LOG(10.0D0)
	   CD(VAR_PLT3)%XVEC(1:K)=T1*CD(VAR_PLT3)%XVEC(1:K)
	  END IF
	  CD(VAR_PLT3)%CURVE_ID=' '
	  GOTO 1000
!
	ELSE IF (ANS .EQ. 'DER')THEN
	  CALL NEW_GEN_IN(VAR_PLT1,'Input plot 1?')
	  VAR_PLT3=NPLTS+1
	  CALL NEW_GEN_IN(VAR_PLT3,'Output plot?')
	  TMP_STR='LINMON'
	  CALL NEW_GEN_IN(TMP_STR,'Derivative option: LINMON or LOGMON')
	  CALL DO_PG_DERIV(VAR_PLT1,VAR_PLT3,TMP_STR)
	  TYPE_CURVE(VAR_PLT3)='L'
	  GOTO 1000
!
	ELSE IF (ANS .EQ. 'SM' .OR. ANS .EQ. 'BXSM')THEN
!
	  WRITE(6,*)BLUE_PEN
	  WRITE(6,*)'Note: This option asumes equally spaced data'
	  WRITE(6,*)'Use the GSM option for unequal data for Gaussian smoothing with fixed resolution (in km/s)'
	  WRITE(6,*)DEF_PEN
	  IF(ANS .EQ. 'BXSM')THEN
	    BOX_FILTER=.TRUE.
	    I=3; CALL NEW_GEN_IN(I,'NPTS > 2 (should be odd)')
	  ELSE
	    BOX_FILTER=.FALSE.
	    I=5; CALL NEW_GEN_IN(I,'NHAN > 2 (should be odd)')
	  END IF
	  IF(I .EQ. 2*(I/2))THEN
	    WRITE(6,*)'Increasing number of points by one to make odd'
	    I=I+1
	  END IF
	  SMOOTH_PLOT=.TRUE.
	  DO IP=1,NPLTS
	    WRITE(TMP_STR,'(I2)')IP
	    TMP_STR='Smooth plot ('//TMP_STR(1:2)//')'
	    CALL NEW_GEN_IN(SMOOTH_PLOT(IP),TRIM(TMP_STR))
	  END DO
	  CALL SMOOTH_PLT(SMOOTH_PLOT,BOX_FILTER,I)
!
        ELSE IF (ANS .EQ. 'GSM')THEN
!
	  WRITE(6,*)RED_PEN
	  WRITE(6,*)' In general the X and Y and axes should be in linear space'
	  WRITE(6,*)' When RESOLUTION is set, constant velocity res. at all X vlaues.'
	  WRITE(6,*)' Otherwise a fixed dX will be used'
	  WRITE(6,*)DEF_PEN
!
	  IP=1; CALL NEW_GEN_IN(IP,'Plot to smooth')
	  OP=NPLTS+1; CALL NEW_GEN_IN(OP,'Output plot')
	  LAM_ST=MIN(XPAR(1),XPAR(2)); LAM_END=MAX(XPAR(1),XPAR(2)) !Passed as R*8
	  CALL DO_GAUSS_SMOOTH(LAM_ST,LAM_END,IP,OP,T_OUT)
	  ERR(OP)=.FALSE.
	  IF(OP .NE. IP .AND .IOS .EQ. 0)THEN
	    TYPE_CURVE(OP)=TYPE_CURVE(IP)
            TYPE_CURVE(OP)='L'
            LINE_WGT(OP)=1
	  END IF
!
	ELSE IF (ANS .EQ. 'FILL')THEN
	  CALL NEW_GEN_IN(FILL,'Fill enclosed areas?')
	  IF(FILL)THEN
	    I=NFILL+1
	    CALL NEW_GEN_IN(I,'Fill identifier - default is next number')
	    IF(IFILL_PLT1(I) .EQ. 0)IFILL_PLT1(I)=1
	    IF(IFILL_PLT2(I) .EQ. 0)IFILL_PLT2(I)=2
	    CALL NEW_GEN_IN(IFILL_PLT1(I),'Input 1st plot')
	    CALL NEW_GEN_IN(IFILL_PLT2(I),'Input 2nd plot')
	    CALL NEW_GEN_IN(IFILL_CLR(I),'Color for fill region')
	    IF(I .EQ. NFILL+1)NFILL=NFILL+1
	  ELSE
	    NFILL=0
	  END IF
	  GOTO 1000
!
	ELSE IF (ANS .EQ. 'CUM')THEN
	  IP=1
	  CALL NEW_GEN_IN(IP,'Input plot 1?')
	  OP=NPLTS+1
	  CALL NEW_GEN_IN(OP,'Output plot?')
	  TYPE_CURVE(OP)='L'
	  REVERSE=.FALSE.
	  CALL NEW_GEN_IN(REVERSE,'Reverse integration direction?')
	  IF(OP .NE. IP .AND. ALLOCATED(CD(OP)%XVEC))THEN
            DEALLOCATE (CD(OP)%XVEC)
            DEALLOCATE (CD(OP)%DATA)
          END IF
          IF(.NOT. ALLOCATED(CD(OP)%XVEC))THEN
	    ALLOCATE (CD(OP)%XVEC(NPTS(IP)),STAT=IOS)
            IF(IOS .EQ. 0)ALLOCATE (CD(OP)%DATA(NPTS(IP)),STAT=IOS)
	  END IF
	  IF(REVERSE)THEN
	    T1=CD(IP)%DATA(NPTS(IP))
            CD(OP)%DATA(NPTS(IP))=0
            DO I=NPTS(IP)-1,1,-1
	      T2=CD(IP)%DATA(I)
	      CD(OP)%DATA(I)=CD(OP)%DATA(I+1)+0.5D0*(T1+T2)*
	1              ABS(CD(IP)%XVEC(I)-CD(IP)%XVEC(I+1))
	      T1=T2
	    END DO
	  ELSE
	    T1=CD(IP)%DATA(1)
            CD(OP)%DATA(1)=0
            DO I=2,NPTS(IP)
	      T2=CD(IP)%DATA(I)
	      CD(OP)%DATA(I)=CD(OP)%DATA(I-1)+0.5D0*(T1+T2)*
	1              ABS(CD(IP)%XVEC(I)-CD(IP)%XVEC(I-1))
	      T1=T2
	    END DO
	  END IF
	  IF(IP .NE. OP)CD(OP)%XVEC=CD(IP)%XVEC
	  NPTS(OP)=NPTS(IP)
	  ERR(OP)=.FALSE.
	  IF(OP .GT. NPLTS)NPLTS=OP
	  GOTO 1000
!
	ELSE IF (ANS .EQ. 'SUM')THEN
	  WRITE(6,'(A)')RED_PEN
	  WRITE(6,'(A)')'This option simply sums the data over X.'
	  WRITE(6,'(A)')'No allowance is made for the X values or spacing'
	  WRITE(6,'(A)')'Use CUM option to intergate over X'
	  WRITE(6,'(A)')DEF_PEN
	  CALL NEW_GEN_IN(IP,'Input plot 1?')
	  OP=NPLTS+1
	  CALL NEW_GEN_IN(OP,'Output plot?')
	  TYPE_CURVE(OP)='L'
	  IF(OP .NE. IP .AND. ALLOCATED(CD(OP)%XVEC))THEN
            DEALLOCATE (CD(OP)%XVEC)
            DEALLOCATE (CD(OP)%DATA)
          END IF
          IF(.NOT. ALLOCATED(CD(OP)%XVEC))THEN
	    ALLOCATE (CD(OP)%XVEC(NPTS(IP)),STAT=IOS)
            IF(IOS .EQ. 0)ALLOCATE (CD(OP)%DATA(NPTS(IP)),STAT=IOS)
	  END IF
          CD(OP)%DATA(1)=CD(IP)%DATA(1)
          DO I=2,NPTS(IP)
	    T2=CD(IP)%DATA(I)
	    CD(OP)%DATA(I)=CD(OP)%DATA(I-1)+CD(IP)%DATA(I)
	    T1=T2
	  END DO
	  IF(IP .NE. OP)CD(OP)%XVEC=CD(IP)%XVEC
	  NPTS(OP)=NPTS(IP)
	  ERR(OP)=.FALSE.
	  IF(OP .GT. NPLTS)NPLTS=OP
	  GOTO 1000
!
	ELSE IF(ANS .EQ. 'NM')THEN
	  VAR_PLT1=1
	  CALL NEW_GEN_IN(VAR_PLT1,'Reference plot for normalization (0 to norm to 1.0')
	  IF(VAR_PLT1 .LT. 0 .OR. VAR_PLT1 .GT. NPLTS)THEN
	    WRITE(T_OUT,*)'Bad plot number'
	    GOTO 1000
	  END IF
	  IF(VAR_PLT1 .EQ. 0)WRITE(T_OUT,*)'Normalizing plots to 1.0'
	  XT(1)=XPAR(1); CALL NEW_GEN_IN(XT(1),'Beginning of normalization range')
	  XT(2)=XPAR(2)
	  CALL NEW_GEN_IN(XT(2),'End of normalization range')
	  MEAN=0.0D0
	  T2=0.0D0
	  IP=VAR_PLT1
	  IF(VAR_PLT1 .NE. 0)THEN
	    DO J=2,NPTS(VAR_PLT1)-1
	      T1=(CD(IP)%XVEC(J)-XT(1))*(XT(2)-CD(IP)%XVEC(J))
	      IF(T1 .GT. 0)THEN
	         MEAN=MEAN+CD(IP)%DATA(J)*ABS(CD(IP)%XVEC(J-1)-CD(IP)%XVEC(J+1))
	         T2=T2+ABS(CD(IP)%XVEC(J-1)-CD(IP)%XVEC(J+1))
	      END IF
	    END DO
	    MEAN=0.5D0*MEAN; T2=0.5D0*T2
	    IF(MEAN .EQ. 0 .OR. T2 .EQ. 0)THEN
	      WRITE(T_OUT,*)'No normalization will be done'
	      WRITE(T_OUT,*)'Bad range of data'
	      GOTO 1000
	    ELSE
	      MEAN=MEAN/T2
	    END IF
	  ELSE
	      MEAN=1.0D0
	  END IF
	  DO IP=1,NPLTS
	    IF(IP .NE. VAR_PLT1 .AND. TYPE_CURVE(IP) .NE. 'I')THEN
	      T1=0.0D0; T2=0.0D0
	      DO J=2,NPTS(IP)-1
	        T3=(CD(IP)%XVEC(J)-XT(1))*(XT(2)-CD(IP)%XVEC(J))
	        IF(T3 .GT. 0)THEN
	          T1=T1+CD(IP)%DATA(J)*ABS(CD(IP)%XVEC(J-1)-CD(IP)%XVEC(J+1))
	          T2=T2+ABS(CD(IP)%XVEC(J-1)-CD(IP)%XVEC(J+1))
	        END IF
	      END DO
	      T1=0.5D0*T1; T2=0.5D0*T2
	      IF(T2 .EQ. 0 .OR. T1 .EQ. 0)THEN
	        WRITE(T_OUT,*)IP,' not normalized'
	      ELSE
	        T1=T1/T2
	        T3=MEAN/T1
	        WRITE(T_OUT,*)'Normalization parameter for plot',IP,' is',T3
	        WRITE(LU_NORM,'(A,I3,A,ES12.4,2X,I3,2ES12.4)')'Normalization parameters for plot',IP,' are',
	1                          T3,VAR_PLT1,XT(1),XT(2)
	        DO J=1,NPTS(IP)
	          CD(IP)%DATA(J)=CD(IP)%DATA(J)*T3
	          IF(ERR(IP))CD(IP)%EMIN(J)=CD(IP)%EMIN(J)*T3
	          IF(ERR(IP))CD(IP)%EMAX(J)=CD(IP)%EMAX(J)*T3
	        END DO
	      END IF
	    END IF
	  END DO
	  GOTO 1000
!
	ELSE IF(ANS .EQ. 'NMS')THEN
	  VAR_PLT1=1
	  CALL NEW_GEN_IN(VAR_PLT1,'Reference plot for normalization ')
	  IF(VAR_PLT1 .EQ. 0 .OR. VAR_PLT1 .GT. NPLTS)THEN
	    WRITE(T_OUT,*)'Bad plot number'
	    GOTO 1000
	  END IF
	  XT(1)=XPAR(1); CALL NEW_GEN_IN(XT(1),'Beginning of normalization range')
	  XT(2)=XPAR(2)
	  CALL NEW_GEN_IN(XT(2),'End of normalization range')
	  IF(DONE_NORMALIZATION)THEN
	    DO IP=1,NPLTS
	      T3=SCALE_FAC(IP)
	      DO J=1,NPTS(IP)
	        CD(IP)%DATA(J)=CD(IP)%DATA(J)/T3
	        IF(ERR(IP))CD(IP)%EMIN(J)=CD(IP)%EMIN(J)/T3
	        IF(ERR(IP))CD(IP)%EMAX(J)=CD(IP)%EMAX(J)/T3
	      END DO
	    END DO
	    TITLE=STORED_TITLE
	    DONE_NORMALIZATION=.FALSE.
	  ELSE
	    STORED_TITLE=TITLE
	  END IF
	  SCALE_FAC=1.0D0
!
	  MEAN=0.0D0
	  T2=0.0D0
	  IP=VAR_PLT1
	  DO J=2,NPTS(VAR_PLT1)-1
	    T1=(CD(IP)%XVEC(J)-XT(1))*(XT(2)-CD(IP)%XVEC(J))
	    IF(T1 .GT. 0)THEN
	       MEAN=MEAN+CD(IP)%DATA(J)*ABS(CD(IP)%XVEC(J-1)-CD(IP)%XVEC(J+1))
	       T2=T2+ABS(CD(IP)%XVEC(J-1)-CD(IP)%XVEC(J+1))
	    END IF
	  END DO
	  MEAN=0.5D0*MEAN; T2=0.5D0*T2
	  IF(MEAN .EQ. 0 .OR. T2 .EQ. 0)THEN
	    WRITE(T_OUT,*)'No normalization will be done'
	    WRITE(T_OUT,*)'Bad range of data'
	    GOTO 1000
	  ELSE
	    MEAN=MEAN/T2
	  END IF
!
	  IPLT=0
	  CALL NEW_GEN_IN(IPLT,K,NPLTS,'Plots to scale -- 0 un-normalizes')
	  DO L=1,K
	    IP=IPLT(L)
	    IF(IP .GT. 0 .AND. IP .LE. NPLTS .AND.
	1      IP .NE. VAR_PLT1 .AND. TYPE_CURVE(IP) .NE. 'I')THEN
	      T1=0.0D0; T2=0.0D0
	      DO J=2,NPTS(IP)-1
	        T3=(CD(IP)%XVEC(J)-XT(1))*(XT(2)-CD(IP)%XVEC(J))
	        IF(T3 .GT. 0)THEN
	          T1=T1+CD(IP)%DATA(J)*ABS(CD(IP)%XVEC(J-1)-CD(IP)%XVEC(J+1))
	          T2=T2+ABS(CD(IP)%XVEC(J-1)-CD(IP)%XVEC(J+1))
	        END IF
	      END DO
	      T1=0.5D0*T1; T2=0.5D0*T2
	      IF(T2 .EQ. 0 .OR. T1 .EQ. 0)THEN
	        WRITE(T_OUT,*)IP,' not normalized'
	      ELSE
	        T1=T1/T2
	        T3=MEAN/T1
	        WRITE(T_OUT,*)'Normalization parameter for plot',IP,' is',T3
	        WRITE(LU_NORM,'(A,I3,A,ES12.4,2X,I3,2ES12.4)')'Normalization parameters for plot',IP,' are',
	1                          T3,VAR_PLT1,XT(1),XT(2)
	        SCALE_FAC(IP)=T3
	        DO J=1,NPTS(IP)
	          CD(IP)%DATA(J)=CD(IP)%DATA(J)*T3
	          IF(ERR(IP))CD(IP)%EMIN(J)=CD(IP)%EMIN(J)*T3
	          IF(ERR(IP))CD(IP)%EMAX(J)=CD(IP)%EMAX(J)*T3
	        END DO
	        DONE_NORMALIZATION=.TRUE.
	      END IF
	    END IF
	    IF(ABS(SCALE_FAC(IP)-1.0D0) .GT. 0.005)THEN
	      J=MIN(INDEX(STORED_TITLE(IP),'\'),LEN_TRIM(STORED_TITLE(IP)))
	      WRITE(TMP_STR,'(A3,F5.3)')'\x  ',SCALE_FAC(IP)
	      IF(J .NE. 0)TITLE(IP)=STORED_TITLE(IP)(1:J-1)//TMP_STR(1:9)//STORED_TITLE(IP)(J:)
	    END IF
	  END DO
	  GOTO 1000
!
	ELSE IF(ANS .EQ. 'SIG')THEN
	  XT(1:2)=XPAR(1:2)
	  CALL NEW_GEN_IN(XT,I,ITWO,'XST,XEND')
	  DO IP=1,NPLTS
	    MEAN=0.0D0; SIGMA=0.0D0; CNT=0
	    DO J=1,NPTS(IP)
	      T1=(CD(IP)%XVEC(J)-XT(1))*(XT(2)-CD(IP)%XVEC(J))
	      IF(T1 .GT. 0)THEN
	        CNT=CNT+1
	        MEAN=MEAN+CD(IP)%DATA(J)
	      END IF
	    END DO
	    IF(CNT .GT. 1)THEN
	      MEAN=MEAN/CNT
	      DO J=1,NPTS(IP)
	        T1=(CD(IP)%XVEC(J)-XT(1))*(XT(2)-CD(IP)%XVEC(J))
	        IF(T1 .GT. 0)THEN
	          SIGMA=SIGMA+(CD(IP)%DATA(J)-MEAN)**2
	        END IF
	      END DO
	      SIGMA=SQRT( SIGMA/(CNT-1) )
	      WRITE(T_OUT,*)IP,MEAN,SIGMA,MEAN/SIGMA
	    END IF
	  END DO
	  GOTO 1000
	ELSE IF(ANS .EQ. 'LNKE')THEN
	  COLUMN(1)=1; COLUMN(2)=2
	  CALL NEW_GEN_IN(COLUMN,I,ITWO,'Data and error plot numbers')
	  IP=COLUMN(1)
	  IF(NPTS(IP) .NE. NPTS(COLUMN(2)))THEN
	    WRITE(6,*)'Can only link error to plot when number of data points is zero'
	    WRITE(6,*)'Number of points in  data vector:',NPTS(IP)
	    WRITE(6,*)'Number of points in error vector:',NPTS(COLUMN(2))
	    GOTO 1000
	  END IF
	  ERR(IP)=.TRUE.
	  IF(ALLOCATED(CD(IP)%EMAX))DEALLOCATE(CD(IP)%EMAX,CD(IP)%EMIN)
	  ALLOCATE(CD(IP)%EMAX(NPTS(IP)),CD(IP)%EMIN(NPTS(IP)))
	  J=COLUMN(2)
	  DO I=1,NPTS(IP)
	    CD(IP)%EMAX(I)=CD(IP)%DATA(I)+CD(J)%DATA(I)
	    CD(IP)%EMIN(I)=CD(IP)%DATA(I)-CD(J)%DATA(I)
	  END DO
	  TYPE_CURVE(J)='I'
	  WRITE(6,*)'Error curve made invisible'
!
	ELSE IF (ANS .NE. 'P' .AND. ANS .NE. 'R')THEN
	  WRITE(T_OUT,*)' Invalid option - try again'
	  GOTO 1000
	END IF
!
! Actual plotting section. All parameters except the plot size have been set.
!
5000	CONTINUE
!
! Open the appropriate plotting device
!
 6000	IF (HARD) CALL PGBBUF
	IF (HARD .AND. LONG_PLOT)THEN
	  T1=LENGTH_OF_HC_PLOT/2.54D0           !in inches
	  CALL PGPAP(T1,LP_ASR)
	END IF
!
! Initialize plotting parameters to the correct values
!
	CALL PGENV(XPAR(1),XPAR(2),YPAR(1),YPAR(2),0,-1)
	IF(HARD)THEN
	  CALL PGSCR(0,1.0,1.0,1.0)  !change background to white
          DO I=1,15                  !Make sure color are set for hard device.
            CALL PGSCR(I,RED(I),GREEN(I),BLUE(I))
          END DO
	  CALL PGSLW(PLT_LINE_WGT)
	ELSE
	  CALL PGQCR(0,RTEMP,GTEMP,BTEMP)
	  CALL PGSLW(IONE)
	END IF
	CALL PGSVP(MARGINX(1),MARGINX(2),MARGINY(1),MARGINY(2))
	CALL PGSWIN(XPAR(1),XPAR(2),YPAR(1),YPAR(2))
	CALL PGERAS
!
! Set the desired aspect ratio of the plot (useful so that
! aspect ratio is preserved on the screen and hardcopy device)
! 0 implies leave aspect ratio at default value
!
	CALL PGQVP(ITWO,DXST,DXEND,DYST,DYEND)
	DASR=(DYEND-DYST)/(DXEND-DXST)
	TEMPASR=DASR
	IF(ASR .LT. 0)TEMPASR=-1.0/ASR
	IF(ASR .GT. 0)TEMPASR=ASR
!
	TOTXMM=(DXEND-DXST)
	TOTYMM=(DYEND-DYST)
 350	IF(HARD) CALL NEW_GEN_IN(XCM,'Plot size (in cm)')
!
! Commented out to test the box squareing routine
! Gregson Vaux December 1996
!
	IF(HARD .AND. (XCM .NE. 0)) THEN
	  SCALEFAC=XCM*10/TOTXMM
	  SCALEFACY=XCM*10*TEMPASR/TOTYMM
	ELSE IF (DASR .LE. 1.0) THEN
	  SCALEFAC=1
	  SCALEFACY=TEMPASR/DASR
          IF(DASR .LT. TEMPASR)THEN
	    T1=DASR/TEMPASR
	    SCALEFAC=SCALEFAC*T1
	    SCALEFACY=SCALEFACY*T1
	  END IF
	ELSE
	  SCALEFAC=TEMPASR/DASR
	  SCALEFACY=1
          IF(DASR .GT. TEMPASR)THEN
	    T1=TEMPASR/DASR
	    SCALEFAC=SCALEFAC*T1
	    SCALEFACY=SCALEFACY*T1
	  END IF
	ENDIF
!
	PRINTX1=MARGINX(1)
	PRINTX2=MARGINX(1)+(MARGINX(2)-MARGINX(1))*SCALEFAC
	PRINTY1=MARGINY(1)
	PRINTY2=MARGINY(1)+(MARGINY(2)-MARGINY(1))*SCALEFACY
        IF(HARD) THEN
	  IF(PRINTX2 .GT. 1.0 .OR. PRINTY2 .GT. 1.0)THEN
	    WRITE(T_OUT,*)' Error - plot too large - try again'
	    GOTO 350
	  ENDIF
	ENDIF
	CALL PGSVP(PRINTX1,PRINTX2,PRINTY1,PRINTY2)
!
! Fiddle with the character size to guarentee that the character size has
! the same relative size to the plot, irrespective of the plot medium.
!
	CALL PGSCH(RONE)
	CALL PGQCS(IFOUR,XCHAR_SIZE,YCHAR_SIZE)
!	T1=ABS(YPAR(2)-YPAR(1))/YCHAR_SIZE/35.0
	T1=(YPAR(2)-YPAR(1))/YCHAR_SIZE/35.0
	EXPCHAR=EXPCHAR_SCALE*T1
	EXPMARK=EXPMARK_SCALE*T1
	TICK_FAC=TICK_FAC_SCALE*T1
	CALL PGSCH(EXPCHAR)
!
! Draw Graphs
!
	IP_ST=1; IP_END=NPLTS; IP_INC=1
	IF(REVERSE_PLOTTING_ORDER)THEN
	  IP_ST=NPLTS; IP_END=1; IP_INC=-1
	END IF
!
! We do the FILL option first so that the curves are drawn on TOP.
!
	IF(FILL)THEN
	  DO I=1,NFILL
	    CALL PGSCI(IFILL_CLR(I))
	    CALL DO_FILL(XPAR,YPAR,IFILL_PLT1(I),IFILL_PLT2(I),.TRUE.)
	  END DO
	END IF
!
! Draw error bars if required.
!
	IF(DO_ERROR)THEN
	  CALL PGSLS(LINE_STYLE(1))
!	  IP_ST=NPLTS; IP_END=1; IP_INC=-1
	  DO IP=IP_ST,IP_END,IP_INC
	    IF(ERR(IP) .AND. TYPE_CURVE(IP) .NE. 'I')THEN
	      Q=PEN_COL(IP+1)
	      CALL PGSCI(Q)
	      DO J=1,NPTS(IP)
	        YT(1)=CD(IP)%EMAX(J)
	        YT(2)=CD(IP)%EMIN(J)
	        XT(1)=CD(IP)%XVEC(J)
	        XT(2)=CD(IP)%XVEC(J)
	        CALL PGLINE(2,XT,YT)
	      END DO
	    END IF
	  END DO
	END IF
!
	DO IP=IP_ST,IP_END,IP_INC
	  CALL PGSLW(LINE_WGT(IP))
	  CALL PGSLS(LINE_STYLE(IP))
	  Q=PEN_COL(IP+PEN_OFFSET)
	  CALL PGSCI(Q)     ! START WITH COLOR INDEX 2 when PEN_OFFSET is 1 (default)
!
! Checks whether we are installing a right axis.
!
	  IF(WHICH_Y_AX(IP) .EQ. 'R')THEN
	    CALL PGSWIN(XPAR(1),XPAR(2),YPAR_R_AX(1),YPAR_R_AX(2))
	  END IF
!
	  IF(TYPE_CURVE(IP) .EQ. 'L' .AND. (MARKER_STYLE(IP) .GE. 0 .OR. .NOT. MARK))THEN
	    CALL PGLINE(NPTS(IP),CD(IP)%XVEC,CD(IP)%DATA)
	  ELSE IF( (TYPE_CURVE(IP) .EQ. 'E' .OR. TYPE_CURVE(IP) .EQ. 'EC')
	1              .AND. (MARKER_STYLE(IP) .GE. 0 .OR. .NOT. MARK))THEN
	    IST=1
	    IEND=2
	    T1=CD(IP)%XVEC(IEND)-CD(IP)%XVEC(IST)
	    Q=PEN_COL(IP+PEN_OFFSET)-1
	    DO WHILE(IEND .LT. NPTS(IP))
	      DO WHILE(IEND .LT. NPTS(IP))
	         IF( (CD(IP)%XVEC(IEND+1)-CD(IP)%XVEC(IEND))/T1
	1                                                   .GE. 0)THEN
	           IEND=IEND+1
	         ELSE
	           EXIT
	         END IF
	      END DO
	      IF(TYPE_CURVE(IP) .EQ. 'EC')THEN
	        Q=Q+1;
	        CALL PGSCI(Q)
	      END IF
	      J=IEND-IST+1
	      CALL PGLINE(J,CD(IP)%XVEC(IST),CD(IP)%DATA(IST))
	      IST=IEND+1
	      IEND=IST+1
	      T1=CD(IP)%XVEC(IEND)-CD(IP)%XVEC(IST)
	      IF(T1 .EQ. 0.0D0)THEN
	        WRITE(6,*)'Zero spacing in GRAMON-PGPLOT ', IST,IEND,IP
	        IEND=IEND+1
	        T1=CD(IP)%XVEC(NPTS(IEND))-CD(IP)%XVEC(IST)
	      END IF
	    END DO
!
	  ELSE IF(TYPE_CURVE(IP) .EQ. 'LG')THEN
	    IST=1
	    IEND=2
	    WRITE(6,*)'Best to use negative marker style'
	    L=ABS(MARKER_STYLE(IP))
	    CALL PGSCH(EXPMARK)
	    IF(ALLOCATED(TA))THEN
	      J=SIZE(TA)
	      IF(J .LT. NPTS(IP))DEALLOCATE(TA)
	    END IF
	    IF(.NOT. ALLOCATED(TA))THEN
	      J=MAXVAL(NPTS(1:NPLTS))
	      ALLOCATE (TA(J))
	    END IF
	    TA(1:NPTS(IP))=0.0D0
	    DO I=1,NPTS(IP)
	      IF(CD(IP)%DATA(I) .NE. 0)THEN
	        TA(I)=LOG10(ABS(CD(IP)%DATA(I)))
	      ELSE
	        TA(I)=-38.0D0
	      END IF
	    END DO
	    DO WHILE(IEND .LT. NPTS(IP))
	      DO WHILE(IEND .LT. NPTS(IP))
	         IF( CD(IP)%DATA(IEND)*CD(IP)%DATA(IEND-1) .GE. 0)THEN
	           IEND=IEND+1
	         ELSE
	           EXIT
	         END IF
	      END DO
	      J=IEND-IST
	      IF(J .NE. 1 .AND. MARKER_STYLE(IP) .GT. 0)CALL PGLINE(J,CD(IP)%XVEC(IST),TA(IST))
	      DO J=IST,IEND-1
	        IF(MARK)CALL PGPT(1,CD(IP)%XVEC(J),TA(J),L)
	        IF(CD(IP)%DATA(J) .LT. 0)CALL PGPT(1,CD(IP)%XVEC(J),TA(J),24)
	      END DO
	      IST=IEND
	      IEND=IST+1
	    END DO
	    CALL PGSCH(EXPCHAR)		!Reset character size
!
	  ELSE IF(TYPE_CURVE(IP) .EQ. 'H' .AND. (MARKER_STYLE(IP) .GE. 0 .OR. .NOT. MARK))THEN
	    CALL PGBIN(NPTS(IP),CD(IP)%XVEC,CD(IP)%DATA,.TRUE.)
!
! This routine does a histogram plot, but the X axis is assumed
! to represent the vertices of each bin, rather than the central
! position.
!
	  ELSE IF(TYPE_CURVE(IP) .EQ. 'A' .AND. (MARKER_STYLE(IP) .GE. 0 .OR. .NOT. MARK))THEN
	    CALL HIST_ADJ(CD(IP)%XVEC,CD(IP)%DATA,NPTS(IP))
!
! This routine does a series of verticle lines extending from YMIN to YV(I)
! at each X value.
!
	  ELSE IF(TYPE_CURVE(IP) .EQ. 'V' .AND. (MARKER_STYLE(IP) .GE. 0 .OR. .NOT. MARK))THEN
	    DO J=1,NPTS(IP)
	      XT(1)=CD(IP)%XVEC(J)
	      XT(2)=CD(IP)%XVEC(J)
	      YT(1)=YPAR(1)
	      YT(2)=CD(IP)%DATA(J)
	      CALL PGLINE(2,XT,YT)
	    END DO
!
	  ELSE IF(TYPE_CURVE(IP) .EQ. 'VB' .AND. (MARKER_STYLE(IP) .GE. 0 .OR. .NOT. MARK))THEN
	    CALL NEW_GEN_IN(VB_BASE(IP),'Base level')
	    DO J=1,NPTS(IP)
	      XT(1)=CD(IP)%XVEC(J)
	      XT(2)=CD(IP)%XVEC(J)
	      YT(1)=VB_BASE(IP)
	      YT(2)=CD(IP)%DATA(J)
	      CALL PGLINE(2,XT,YT)
	    END DO
!
! This routine does a broken plot. NPTS should be even. Lines are drawn
! from I to I+1 point for I odd only.
!
	  ELSE IF(TYPE_CURVE(IP) .EQ. 'B' .AND. (MARKER_STYLE(IP) .GE. 0 .OR. .NOT. MARK))THEN
	    CALL PGSLS(LINE_STYLE(IP))
	    Q=PEN_COL(IP+1)
	    CALL PGSCI(Q)     ! START WITH COLOR INDEX 2
	    DO J=1,NPTS(IP),2
	      XT(1)=CD(IP)%XVEC(J)
	      XT(2)=CD(IP)%XVEC(J+1)
	      YT(1)=CD(IP)%DATA(J)
	      YT(2)=CD(IP)%DATA(J+1)
	      CALL PGLINE(2,XT,YT)
	    END DO
	  END IF
!
! If this curve was plotted according to the right hand-axis, we need to
! reset the scle to the left axis.
!
	  IF(WHICH_Y_AX(IP) .EQ. 'R')THEN
	    CALL PGSWIN(XPAR(1),XPAR(2),YPAR(1),YPAR(2))
	  END IF
!
	END DO
	CALL PGSLW(PLT_LINE_WGT)
!
! Draw Gaussian fits if needed.
!
!	IF(HARD .AND. DRAW_GAUSS_HARD)CALL DRAW_GAUS(L_FALSE)
!
	IF(MARK)THEN
	  CALL PGSCH(EXPMARK)
	  DO IP=IP_ST,IP_END,IP_INC
	    L=ABS(MARKER_STYLE(IP))
!
! Don't draw if invisible curve.
!
	    IF(L .NE. 0 .AND. TYPE_CURVE(IP) .NE. 'I' .AND.
	1                     TYPE_CURVE(IP) .NE. 'LG')THEN
!	      CALL PGSLS(LINE_STYLE(IP))
	      Q=PEN_COL(IP+1)
	      CALL PGSCI(Q)     ! START WITH COLOR INDEX 2
	      DO J=1,NPTS(IP)
		CALL PGPT(1,CD(IP)%XVEC(J),CD(IP)%DATA(J),L)
	      END DO
	    END IF
	  END DO
	  CALL PGSCH(EXPCHAR)		!Restore to orig status.
	END IF
!
! Draw borders.
!
	I=1
	CALL PGSLS(I)
	Q=PEN_COL(1)
	CALL PGSCI(Q)   ! Draw box in color index 1
!
! *** PGBOX is commented out while MON_BORD is being tested
! September 1996 Gregson Vaux
!
!	CALL PGBOX('ABCNT',0.0,0,'ABCNT',0.0,0)
!
	CALL ARRANGE_PG_CURVE_IDS(PEN_COL,PEN_OFFSET,MAXPEN,ADD_COMMA,RESET_CURVE_LAB)
	IF(DO_BORDER)THEN
	  CALL MONBORD_V4(XPAR,XINC,XNUMST,IXTICK,IDX,
	1             YPAR,YINC,YNUMST,IYTICK,IDY,
	1             TICK_FAC,EXPCHAR,
	1             XLABEL,YLABEL,TITLE,N_TITLE,TITONRHS,
	1             LOG_AXIS,OPTION,NORMAL_R_Y_AXIS,
	1             XLAB_FILE,YLAB_FILE)
	  IF(.NOT. NORMAL_R_Y_AXIS)THEN
	    CALL DRAW_RIGHT_Y_AXIS(YPAR_R_AX,YINC_R_AX,YNUMST_R_AX,
	1          IYTICK_R_AX,IDY_R_AX,TICK_FAC,
	1          EXPCHAR,YLABEL_R_AX,LOG_AXIS)
	  END IF
	END IF
!
! Draw strings on graph. The parameters T1 and T2 ensure that only strings
! inside the BOX are plotted. Strings can be ploted outside the box by setting
! LOC_PG negative.
!
	IF(STR)THEN
	  WRITE(6,*)'Calling Justify'
	  CALL JUSTIFY_CONVERT_V2(XSTR,YSTR,LOC,LOC_PG,ORIENTATION,FLAGSTR,
     *    XSTRPOS,YSTRPOS,STRING,MAXSTR)
	  T1=(XSTRPOS(I)-XPAR(1))*(XPAR(2)-XSTRPOS(I))
	  T2=(YSTRPOS(I)-YPAR(1))*(YPAR(2)-YSTRPOS(I))
	  DO I=1,MAXSTR
	    IF(LOC(I) .LT. 0)THEN ; T1=1.0; T2=1.0; END IF
	    IF(FLAGSTR(I))THEN  ! .AND. T1 .GT. 0 .AND. T2 .GT. 0)THEN
	      CALL PGSCI(STR_COL(I))
	      CALL PGSCH(EXPCHAR*STR_EXP(I))
	      WRITE(6,*)'Calling PUT_TEXT'
!	        CALL PGPTXT(XSTRPOS(I),YSTRPOS(I),ORIENTATION(I),LOC_PG(I),STRING(I))
	      CALL PUT_TEXT(XSTRPOS(I),YSTRPOS(I),ORIENTATION(I),LOC_PG(I),STRING(I))
	    END IF
	  END DO
	END IF
!
	IF(N_LINE_IDS .NE. 0)THEN
	  CALL DRAW_LINE_IDS_V2(XPAR,YPAR,EXPCHAR,IONE,LINE_CUT_PARAM,T_OUT)
	  CALL PGSCH(EXPCHAR)		!Reset character size
	END IF
!	IF(N_EW_IDS .NE. 0)THEN
!	  IF(EW_SCALE_FAC .GT. 0.0)THEN
!	    CALL DRAW_EW_LINES(XPAR,YPAR,EXPCHAR,T_OUT)
!	  ELSE
!	    CALL DRAW_EW_IDS(XPAR,YPAR,EXPCHAR,T_OUT)
!	  END IF
!	  CALL PGSCH(EXPCHAR)		!Reset character size
!	END IF
!
! Draw vectors on graphs.
!
	IF(VEC)THEN
	  DO I=1,MAXVEC
	    CALL PGSLS(VEC_LINE_STYLE(I))
	    Q=VECPEN(I)
	    CALL PGSCI(Q)
	    IF(FLAGLINE(I))THEN
	      CALL PGMOVE(LINEXST(I),LINEYST(I))
	      CALL PGDRAW(LINEXEND(I),LINEYEND(I))
	    END IF
	  END DO
	END IF
!
! Centre text output on terminal
!
	CALL PGMOVE((XPAR(2)-XPAR(1))/2.0,(YPAR(2)-YPAR(1))/2.0)
	IF(HARD) THEN
	  CALL PGQINF('FILE',CUR_HARD_FILE,LENGTH)
	  CALL PGEBUF		!send plot to printer file
	  CALL PGCLOS
	  CALL PGSLCT(ID)
	  HARD=.FALSE.
	  IF(LONG_PLOT)CALL MODIFY_PGI_PS(CUR_HARD_FILE,LENGTH_OF_HC_PLOT,LP_ASR)
	END IF
	CALL PGSCR(0,RTEMP,GTEMP,BTEMP)
!
	GO TO 1000
	END
!
! Function to indicate when cursor input is finished. May be system dependent,
! hence function call. Mouse must be used with XWIN.
!
! On VMS 2 and 3rd mouse buttons return control to the screen.
! 1st butoon leaves control in plt window.
!
! Mouse buttons on VMS return A, D, and X.
!
! Thus use 2nd button for input.
!          3rd button to finish cursor input.
!
	LOGICAL FUNCTION END_CURS(CURSVAL)
	CHARACTER*1 CURSVAL
	END_CURS=.TRUE.
!
	IF(CURSVAL .EQ. 'A')THEN
          END_CURS=.FALSE.
	ELSE IF(CURSVAL .EQ. 'a')THEN
          END_CURS=.FALSE.
	ELSE IF(CURSVAL .EQ. 'D')THEN
          END_CURS=.FALSE.
	ELSE IF(CURSVAL .EQ. 'd')THEN
          END_CURS=.FALSE.
	ELSE IF(CURSVAL .EQ. ' ')THEN
          END_CURS=.FALSE.
	END IF
!
	RETURN
	END
