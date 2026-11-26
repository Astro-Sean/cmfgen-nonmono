!
! Subroutine to measure line EWs using direct numerical integration of the data.
! Routine is designed to be called in GRAMON_PGPLOT.
!
! Options:
!      (1) Cursor control
!            Using cursors the user specifies the range limits (X coordinates).
!            Two chocies for Y.
!                 Y location indicates the continuum
!             or  Continuum defined by average (over NPIX) of the data at the X coordinates.
!
!      (2) Read in range limits from file. In this case the continuum is ALWAYS defined
!             by an integration centered on the line limits.
!
	SUBROUTINE DO_CURSOR_EW_V3(
	1             XPAR,XINC,XNUMST,IXTICK,IDX,
	1             YPAR,YINC,YNUMST,IYTICK,IDY,
	1             TICK_FAC,EXPCHAR,
	1             XLABEL,YLABEL,TITONRHS,
	1             LOG_AXIS,OPTION,NORMAL_R_Y_AXIS,
	1             XLAB_FILE,YLAB_FILE,
	1             PEN_COL,PEN_OFFSET,REVERSE_PLOT_ORDER)
	USE SET_KIND_MODULE
	USE MOD_CURVE_DATA
	USE MOD_EW_VARIABLES
	USE GEN_IN_INTERFACE
	USE MOD_COLOR_PEN_DEF
	IMPLICIT NONE
!
! Altered: 22-Jul-2023 - FWHM was not correctly set before call to GET_LINE_ID_PG.
!                        FWHM computed using 50% points
!                        Improved EW computation.
! Altered: 22-Jul-2022 - Correct FWHM now passed to GET_LINE_ID_PG.
! Altered: 08-Jul-2022 - Only lower case options allowed.
! Altered: 07-Jul-2022 - Call updated, and cleaned.
!                          Plot can be shifted
! Altered: 30-Jun-2022 - Can no append transition name to EW file.
! Created: 27-FEb-2022
!
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: ITWO=2
	INTEGER, PARAMETER :: ISIX=6
!
	REAL*4 XPAR(2),XT(2)
	REAL*4 YPAR(2),YT(2)
        REAL*4 XINC,XNUMST
        REAL*4 YINC,YNUMST
	REAL*4 dY
        REAL*4 EXPCHAR,TICK_FAC
!
        INTEGER IDX,IXTICK
        INTEGER IDY,IYTICK
	INTEGER IY_DIR
!
	LOGICAL TITONRHS
	LOGICAL NORMAL_R_Y_AXIS
!
	CHARACTER(LEN=80) XLABEL,YLABEL
        CHARACTER(LEN=5)  LOG_AXIS
	CHARACTER(LEN=80) OPTION
	CHARACTER(LEN=80) XLAB_FILE,YLAB_FILE
!
	INTEGER, SAVE :: NPIX=3  	!Integraton band pass around line limits
	INTEGER, SAVE :: IP=1
	REAL*4,  SAVE :: LOW_CONT=0.998
	REAL*4,  SAVE :: HIGH_CONT=1.002
!
	INTEGER PEN_COL(0:MAX_PLTS)
	INTEGER PEN_OFFSET
	LOGICAL REVERSE_PLOT_ORDER
!
	INTEGER IST,IEND		!Line limits in pixel space
	REAL*4 XLOC,YLOC		!Used to read cursor location
	REAL*4 YST,YEND  		!Continuum flux at line limits
!
! In the following I use Y and F interchangably.
! Also X will normally be Lambda.
!
	REAL*4 dX               !X spacing (pixel centered)
	REAL*4 XVAL  		!Current X value
	REAL*4 YVAL  		!Current Y value
	REAL*4 YINT
	REAL*4 SLOPE  		!Continuum slope
	REAL*4 T1		!Work variable
!
	REAL*4, SAVE :: CONT_ACC=0.2		!percentage error in continuum
!
        INTEGER PGCURS
        INTEGER CURSERR
	INTEGER CUR_LW
	INTEGER, SAVE :: IPEN=1
	CHARACTER(LEN=1) CURSVAL
!
        LOGICAL, PARAMETER :: L_TRUE=.TRUE.
	LOGICAL, SAVE :: FIRST_WRITE=.TRUE.
	INTEGER, SAVE :: LUIN=0
	INTEGER, SAVE :: LUOUT=0
!
	INTEGER I,J,K
!
! External function.
!
	INTEGER GET_INDX_SP
	CHARACTER(LEN=30) UC
	EXTERNAL GET_INDX_SP
	EXTERNAL UC
!
	LOGICAL END_FILE
	LOGICAL FILE_PRES
	LOGICAL, SAVE ::  RESET_DEFAULTS=.TRUE.
!
	CHARACTER(LEN=80) LOC_FILE_WITH_LINE_LIMS
	CHARACTER(LEN=80) OUT_FILE
	CHARACTER(LEN=80) STRING
!
! XLOC, YLOC will define initial cursor location.
!
	XLOC=0.5D0*(XPAR(1)+XPAR(2))
	YLOC=0.5D0*(YPAR(1)+YPAR(2))
	dY=0.05*(YPAR(2)-YPAR(1))
	IY_DIR=1
!
	IF(.NOT. RESET_DEFAULTS)THEN
	  CALL GEN_IN(RESET_DEFAULTS,'Reset plot and parameters/settings?')
	END IF
	IF(RESET_DEFAULTS)THEN
	  IF(NPLTS .EQ. 1)THEN
	    IP=1
	  ELSE
	    CALL GEN_IN(IP,'Plot for fitting')
	  END IF
!
! Used if average data at X-limits is used to define the continuum level.
!
	  NPIX=1
	  CALL GEN_IN(NPIX,'Number of pixels at X location to average continuum (must be odd)')
	  IF(MOD(NPIX,2) .EQ. 0)THEN
	    NPIX=NPIX+1
	    WRITE(6,*)'NPIX increase by 1 to make odd; NPIX=',NPIX
	  END IF
	  CALL GEN_IN(USE_MILLI_ANG,'Output EWs in milli-Angstroms?')
	  CALL GEN_IN(CONT_ACC,'Measure accuracy for continuum -- percentage?')
	  LOW_CONT=1.0-CONT_ACC/100.0D0
	  HIGH_CONT=1.0+CONT_ACC/100.0D0
	  RESET_DEFAULTS=.FALSE.
	END IF
!
! Used as a work array -- needed for getting FWHM.
!
        ALLOCATE (ONE_MIN_FDFC(NPTS(IP)))
!
! Open output file. Data is appended if it alread exists.
!
	IF(LUOUT .EQ. 0)THEN
	  OUT_FILE='EW_DATA'
	  CALL GET_LU(LUOUT,'LUOUT in DO_CURSOR_EW_V2')
	  CALL GEN_IN(OUT_FILE,'File to OUTPUT EWs etc')
	  INQUIRE(FILE=OUT_FILE,EXIST=FILE_PRES)
	  IF(FILE_PRES)THEN
	    WRITE(6,*)'File already exists -- appending new data'
	    OPEN(UNIT=LUOUT,FILE=OUT_FILE,STATUS='OLD',ACTION='WRITE',POSITION='APPEND')
	  ELSE
	    OPEN(UNIT=LUOUT,FILE=OUT_FILE,STATUS='NEW',ACTION='WRITE')
	    CALL WRITE_EW_HEADER(LUOUT)
	  END IF
	END IF
!
	WRITE(6,'(A)')' '
	I=6; CALL WRITE_EW_HEADER(I)
	WRITE(6,'(A)')' '
!
	CALL PRINT_CURSOR_DESC
	DO WHILE(1 .EQ. 1)
!
! Change pen colors for defining limits of the next line.
!
	  IF(IPEN .EQ. 6)THEN
	    IPEN=1
	  ELSE
	    IPEN=6
	  END IF
	  CALL PGSCI(IPEN)
!
! Get start of line.
!
	  CURSERR = PGCURS(XLOC,YLOC,CURSVAL)
	  IST=GET_INDX_SP(XLOC,CD(IP)%XVEC,NPTS(IP))
	  XST=XLOC; YST=YLOC
	  IF(XST .LT. XPAR(1))THEN
	    WRITE(6,*)'Invalid X start location -- start again'
	    GOTO 1000
	  END IF
!
	  IF(UC(CURSVAL) .EQ. 'Q')THEN
	    EXIT
	  ELSE IF(UC(CURSVAL) .EQ. 'S')THEN
	    GOTO 1000
	  ELSE IF(UC(CURSVAL) .EQ. 'Y')THEN
	    CALL PGPT(IONE,XLOC,YLOC,ITWO)
!
! Y value set by average at cursor location
!
	  ELSE IF(UC(CURSVAL) .EQ. 'A')THEN
	    YST=0.0
	    DO J=IST-NPIX/2,IST+NPIX/2
	      YST=YST+CD(IP)%DATA(J)
	    END DO
	    YST=YST/NPIX
	    YT(1)=YST-dY; YT(2)=YST+dY
	    XT(1)=XLOC; XT(2)=XLOC
	    CALL PGLINE(2,XT,YT)
!
! Options to move/expand/contract plot in X and Y directions.
!
	  ELSE IF(UC(CURSVAL) .EQ. 'N' .OR.
	1         UC(CURSVAL) .EQ. 'P' .OR.
	1         UC(CURSVAL) .EQ. 'X' .OR.
	1         UC(CURSVAL) .EQ. 'C' .OR.
	1         UC(CURSVAL) .EQ. 'U' .OR.
	1         UC(CURSVAL) .EQ. 'L' .OR.
	1         UC(CURSVAL) .EQ. 'T')THEN
	    T1=XPAR(2)-XPAR(1)
	    IF(UC(CURSVAL) .EQ. 'N')THEN
	      XPAR(1)=XPAR(2)-XINC
	      XPAR(2)=XPAR(1)+T1
	      XLOC=XPAR(1)+XINC
	    ELSE IF(UC(CURSVAL) .EQ. 'P')THEN
	      XPAR(1)=XPAR(1)+XINC-(XPAR(2)-XPAR(1))
	      XPAR(2)=XPAR(1)+T1
	      XLOC=XPAR(1)+XINC
	    ELSE IF(UC(CURSVAL) .EQ. 'X')THEN
	      XPAR(2)=XPAR(2)+XINC
	      XLOC=XPAR(1)+XINC
	    ELSE IF(UC(CURSVAL) .EQ. 'C')THEN
	      XPAR(2)=XPAR(2)-XINC
	      XLOC=XPAR(1)+XINC
	    ELSE IF(UC(CURSVAL) .EQ. 'U')THEN
	      YPAR(2)=YPAR(2)-YINC*IY_DIR
	      IF(YLOC .GT. YPAR(2))YLOC=YPAR(2)-YINC
	    ELSE IF(UC(CURSVAL) .EQ. 'L')THEN
	      YPAR(1)=YPAR(1)+YINC*IY_DIR
	      IF(YLOC .LT. YPAR(1))YLOC=YPAR(1)+YINC
	    ELSE IF(UC(CURSVAL) .EQ. 'T')THEN
	      IY_DIR=-1*IY_DIR
	    END IF
	    CALL PGSWIN(XPAR(1),XPAR(2),YPAR(1),YPAR(2))
            CALL PGERAS
	    IF(REVERSE_PLOT_ORDER)THEN
	      DO J=NPLTS,1,-1
	        K=PEN_COL(J)+PEN_OFFSET
	        CALL PGSCI(K)
                CALL PGLINE(NPTS(J),CD(J)%XVEC,CD(J)%DATA)
	      END DO
	    ELSE
	      DO J=1,NPLTS
	        K=PEN_COL(J)+PEN_OFFSET
	        CALL PGSCI(K)
                CALL PGLINE(NPTS(J),CD(J)%XVEC,CD(J)%DATA)
	      END DO
	    END IF
!
! Add line IDs
!
	    CALL DRAW_LINE_IDS_V2(XPAR,YPAR,EXPCHAR,IONE,ISIX)
            CALL PGSCH(EXPCHAR)		!Reset character size
!
	    CALL PGSCI(IONE)
            CALL MONBORD_V4(XPAR,XINC,XNUMST,IXTICK,IDX,
	1           YPAR,YINC,YNUMST,IYTICK,IDY,
	1           TICK_FAC,EXPCHAR,
	1           XLABEL,YLABEL,TITLE,N_TITLE,TITONRHS,
	1           LOG_AXIS,OPTION,NORMAL_R_Y_AXIS,
	1           XLAB_FILE,YLAB_FILE)
	    GOTO 1000
	  ELSE
	    WRITE(6,*)RED_PEN,' Cursor value not recognized'	
	    CALL PRINT_CURSOR_DESC
	    GOTO 1000
	  END IF
!
! Get end of line.
!
	  CURSERR = PGCURS(XLOC,YLOC,CURSVAL)
	  IF(XLOC .GT. XPAR(2))THEN
	    WRITE(6,*)'Invalid X end location -- start again'
	    GOTO 1000
	  END IF
	  IEND=GET_INDX_SP(XLOC,CD(IP)%XVEC,NPTS(IP))
	  XEND=XLOC; YEND=YLOC
	  IF(UC(CURSVAL) .EQ. 'Q')THEN
	    EXIT
	  ELSE IF(UC(CURSVAL) .EQ. 'S')THEN
	    WRITE(6,*)'Restarting line selection again'
	    GOTO 1000
	  ELSE IF(UC(CURSVAL) .EQ. 'Y')THEN
	    CALL PGPT(IONE,XLOC,YLOC,ITWO)
	  ELSE IF(UC(CURSVAL) .EQ. 'A')THEN
	    YEND=0.0
	    DO J=IEND-NPIX/2,IEND+NPIX/2
	      YEND=YEND+CD(IP)%DATA(J)
	    END DO
	    YEND=YEND/NPIX
	    YT(1)=YEND-dY; YT(2)=YEND+dY
	    XT(1)=XLOC; XT(2)=XLOC
	    CALL PGLINE(2,XT,YT)
	  ELSE
	    WRITE(6,*)RED_PEN,' Cursor value not recognized'	
	    CALL PRINT_CURSOR_DESC
	    GOTO 1000
	  END IF
!
	  IF(IST .EQ. IEND)THEN
	    WRITE(6,*)RED_PEN
	    WRITE(6,*)'Error -- limits on line are identical'
	    WRITE(6,*)'Need to redefine the line'
	    WRITE(6,*)DEF_PEN
	    GOTO 1000
	  END IF
!
! We now determine the line parameters
!
	  EW=0.0; EWL=0.0; EWH=0.0; XMEAN=0.0
	  SLOPE=(YEND-YST)/(CD(IP)%XVEC(IEND)-CD(IP)%XVEC(IST))
	  DO I=IST,IEND
	    YVAL=YST+(CD(IP)%XVEC(I)-CD(IP)%XVEC(IST))*SLOPE
	    ONE_MIN_FDFC(I)=1.0D0-CD(IP)%DATA(I)/YVAL
	    dX=(CD(IP)%XVEC(MIN(I+1,IEND))-CD(IP)%XVEC(MAX(IST,I-1)))/2
	    EW = EW+dX*(1.0D0-CD(IP)%DATA(I)/YVAL)
	    EWL=EWL+dX*(1.0D0-CD(IP)%DATA(I)/(LOW_CONT*YVAL))
            EWH=EWH+dX*(1.0D0-CD(IP)%DATA(I)/(HIGH_CONT*YVAL))
	    XMEAN=XMEAN+dX*CD(IP)%XVEC(I)*ONE_MIN_FDFC(I)
	  END DO
!
	  IF(ABS(EW)/YST .LT. 1.0D-10)THEN
	    WRITE(6,*)'Possibe error EW is close to zero', EW/YST
	    WRITE(6,*)'Skipping this line. XST,END=',XST,XEND
	    GOTO 1000
	  END IF
	  XMEAN=XMEAN/EW
	  YCONT=YST+(XMEAN-CD(IP)%XVEC(IST))*SLOPE
	  YINT=EW
!
! These parameters are used to provide information on whether a line is blended.
!
	  SIGMA=0.0D0; SKEWNESS=0.0D0; KURTOSIS=0.0D0
	  DO I=IST,IEND
	    XVAL=CD(IP)%XVEC(I)
	    dX=(CD(IP)%XVEC(MIN(I+1,IEND))-CD(IP)%XVEC(MAX(IST,I-1)))/2
	    SIGMA=SIGMA+dX*ONE_MIN_FDFC(I)*(XVAL-XMEAN)**2
	    SKEWNESS=SKEWNESS+dX*ONE_MIN_FDFC(I)*(XVAL-XMEAN)**3
	    KURTOSIS=KURTOSIS+dX*ONE_MIN_FDFC(I)*(XVAL-XMEAN)**4
	  END DO
	  IF(SIGMA*YINT .LE. 0)THEN
	    SIGMA=-1.0; SKEWNESS=-1.0; KURTOSIS=-1.0
	  ELSE
	    SIGMA=SQRT(SIGMA/EW)
	    SKEWNESS=SKEWNESS/EW/SIGMA**3
	    KURTOSIS=KURTOSIS/EW/SIGMA**4
	  END IF
!
! Get FWHM. At present, this will only work well for
! isolated lines in theoretical spectra.
!
	  IF(EW .LT. 0.0D0)ONE_MIN_FDFC=-ONE_MIN_FDFC
	  I=MAXLOC(ONE_MIN_FDFC(IST:IEND),IONE)+IST-1
	  YINT=ONE_MIN_FDFC(I)
	  ONE_MIN_FDFC(IST:IEND)=(ONE_MIN_FDFC(IST:IEND)-0.5*YINT)
	  DO J=I+1,IEND
	    IF(ONE_MIN_FDFC(J) .LT. 0.0D0)THEN
	      T1=ONE_MIN_FDFC(J)/(ONE_MIN_FDFC(J)-ONE_MIN_FDFC(J-1))
	      XHIGH_FWHM=T1*CD(IP)%XVEC(J-1)+(1.0D0-T1)*CD(IP)%XVEC(J)
	      EXIT
	    END IF
	  END DO
	  DO J=I-1,IST,-1
	    IF(ONE_MIN_FDFC(J) .LT. 0.0D0)THEN
	      T1=ONE_MIN_FDFC(J)/(ONE_MIN_FDFC(J)-ONE_MIN_FDFC(J+1))
	      XLOW_FWHM=T1*CD(IP)%XVEC(J+1)+(1.0D0-T1)*CD(IP)%XVEC(J)
	      EXIT
	    END IF
	  END DO
	  FWHM=2.998D+05*(XHIGH_FWHM-XLOW_FWHM)/XMEAN			!km/s
	  CALL GET_LINE_ID_PG(TRANS_NAME,LINE_WAVE,EW,XMEAN,FWHM)
!	  CALL GET_LINE_ID_PG(TRANS_NAME,LINE_WAVE,EW,XMEAN,T1)
!
	  CALL WR_EW_VARIABLES(IP,LUOUT)
	  I=6; CALL WR_EW_VARIABLES(IP,I)
!
1000	  CONTINUE
	END DO
!
2000	CONTINUE
	DEALLOCATE(ONE_MIN_FDFC)
	RETURN
!
	CONTAINS
	SUBROUTINE PRINT_CURSOR_DESC
	USE SET_KIND_MODULE
!
! We restrict options to lower case to avoid accidental use of
! the mouse which returns A.
!
	  WRITE(6,'(A)')BLUE_PEN
	  WRITE(6,'(A)')' Program assumes X axis increases with X index '
	  WRITE(6,'(A)')' Cursor controls are case Insensitive'
	  WRITE(6,'(A)')' '
	  WRITE(6,'(A)')' Use cursor to define left and right side of line and continuum option:'
	  WRITE(6,'(A)')'   Y - uses Y cursor value for continuum location'
	  WRITE(6,'(A)')'   A - Y value is average over data at X location'
	  WRITE(6,'(A)')'   S - reset and start line selection again'
	  WRITE(6,'(A)')'   Q - quit line selection'
	  WRITE(6,'(A)')' '
	  WRITE(6,'(A)')'Plot movement'
	  WRITE(6,'(A)')'   N(ext)      -- shift window to left by X increment'
	  WRITE(6,'(A)')'   P(revious)  -- shift window to right by X increment'
	  WRITE(6,'(A)')'   X(extend)   -- expand plotted region by y increment but no shift'
	  WRITE(6,'(A)')'   C(ontract)  -- contract plotted region by by increment but no shift'
	  WRITE(6,'(A)')'   U(pper)     -- expand (shrink) top of window by increment but no shift'
	  WRITE(6,'(A)')'   L(lower))   -- expand (shrink) bottom of window by y increment but no shift'
	  WRITE(6,'(A)')'   T(oggle))   -- Toggle direction of U and L directions'
	  WRITE(6,'(A)')DEF_PEN
!
	  RETURN
	END SUBROUTINE PRINT_CURSOR_DESC
	END SUBROUTINE DO_CURSOR_EW_V3
