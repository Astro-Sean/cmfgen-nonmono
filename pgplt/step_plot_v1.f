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
	SUBROUTINE STEP_PLOT_V1(
	1             XPAR,XINC,XNUMST,IXTICK,IDX,
	1             YPAR,YINC,YNUMST,IYTICK,IDY,
	1             TICK_FAC,EXPCHAR,
	1             XLABEL,YLABEL,TITONRHS,
	1             LOG_AXIS,OPTION,NORMAL_R_Y_AXIS,
	1             XLAB_FILE,YLAB_FILE,
	1             LINE_STYLE,LINE_WGT,
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
	INTEGER LINE_WGT(MAX_PLTS)
	INTEGER LINE_STYLE(MAX_PLTS)
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
!
	DO WHILE(1 .EQ. 1)
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
!
! Options to move/expand/contract plot in X direction.
!
	  ELSE IF(UC(CURSVAL) .EQ. 'N' .OR.
	1         UC(CURSVAL) .EQ. 'P' .OR.
	1         UC(CURSVAL) .EQ. 'X' .OR.
	1         UC(CURSVAL) .EQ. 'C')THEN
	    T1=XPAR(2)-XPAR(1)
	    IF(UC(CURSVAL) .EQ. 'N')THEN
	      XPAR(1)=XPAR(2)-XINC
	      XPAR(2)=XPAR(1)+T1
	    ELSE IF(UC(CURSVAL) .EQ. 'P')THEN
	      XPAR(1)=XPAR(1)+XINC-(XPAR(2)-XPAR(1))
	      XPAR(2)=XPAR(1)+T1
	    ELSE IF(UC(CURSVAL) .EQ. 'X')THEN
	      XPAR(2)=XPAR(2)+XINC
	    ELSE IF(UC(CURSVAL) .EQ. 'C')THEN
	      XPAR(2)=XPAR(2)-XINC
	    END IF
	    XLOC=XPAR(1)+XINC
	    CALL PGSWIN(XPAR(1),XPAR(2),YPAR(1),YPAR(2))
            CALL PGERAS
	    IF(REVERSE_PLOT_ORDER)THEN
	      DO J=NPLTS,1,-1
	        K=PEN_COL(J)+PEN_OFFSET
	        CALL PGSCI(K)
	        CALL PGSLS(LINE_STYLE(J))
                CALL PGLINE(NPTS(J),CD(J)%XVEC,CD(J)%DATA)
	      END DO
	    ELSE
	      DO J=1,NPLTS
	        K=PEN_COL(J)+PEN_OFFSET
	        CALL PGSCI(K)
	        CALL PGSLS(LINE_STYLE(J))
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
1000	  CONTINUE
	END DO
!
2000	CONTINUE
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
	  WRITE(6,'(A)')' Cursor controls are case INsensitive'
	  WRITE(6,'(A)')' '
	  WRITE(6,'(A)')' Use cursor to define left and right side of line:'
	  WRITE(6,'(A)')'   Y - uses Y cursor value for continuum location'
	  WRITE(6,'(A)')'   S - reset and start line selection again'
	  WRITE(6,'(A)')'   A - Y value is average over data at X location'
	  WRITE(6,'(A)')'   Q - quit line selection'
	  WRITE(6,'(A)')' '
	  WRITE(6,'(A)')'Plot movement'
	  WRITE(6,'(A)')'   N(ext)      -- shift window to left by X increment'
	  WRITE(6,'(A)')'   P(revious)  -- shift window to right by X increment'
	  WRITE(6,'(A)')'   X(extend)   -- expand window by increment but no shift'
	  WRITE(6,'(A)')'   C(ontract)  -- contract window by increment but no shift'
	  WRITE(6,'(A)')DEF_PEN
!
	  RETURN
	END SUBROUTINE PRINT_CURSOR_DESC
	END SUBROUTINE STEP_PLOT_V1
