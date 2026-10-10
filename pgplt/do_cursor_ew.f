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
	SUBROUTINE DO_CURSOR_EW(XRANGE,YRANGE)
	USE SET_KIND_MODULE
	USE MOD_CURVE_DATA
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
! Created: 27-FEb-2022
!
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: ITWO=2
!
	REAL*4 XRANGE(2),XT(2)
	REAL*4 YRANGE(2),YT(2)
	REAL*4 dY
!
	INTEGER NPIX  			!Integraton band pass around line limits
	INTEGER IST,IEND		!Line limits in pixel space
	REAL*4 XLOC,YLOC		!Used to read cursor location
	REAL*4 XST,XEND                 !Line limits in world coordinates
	REAL*4 YST,YEND  		!Continuum flux at line limits
!
! In the following I use Y and F interchangably.
! Also X will normally be Lambda.
!
	REAL*4 EW    		!Line equivalent width
	REAL*4 EWL,EWH  	!Line equivalent width with a 1% shift in continuum
	REAL*4 dX               !X spacing (pixel centered)
	REAL*4 XMEAN            !
	REAL*4 YMEAN  		!Average value of (F-Fc)
	REAL*4 XVAL  		!Current X value
	REAL*4 YVAL  		!Current Y value
	REAL*4 YINT
	REAL*4 YCONT  		!Flux at line center
	REAL*4 SLOPE  		!Continuum slope
	REAL*4 T1		!Work variable
!
	REAL*4 SIGMA		!Standard deviation
	REAL*4 KURTOSIS
	REAL*4 SKEWNESS
!
	LOGICAL USE_CURSOR
        INTEGER PGCURS
        INTEGER CURSERR
	INTEGER CUR_LW
	INTEGER IPEN
	CHARACTER(LEN=1) CURSVAL
!
	INTEGER, SAVE :: IP=1
        LOGICAL, PARAMETER :: L_TRUE=.TRUE.
	LOGICAL, SAVE :: FIRST_WRITE=.TRUE.
	INTEGER I,J,K
	INTEGER, SAVE :: LUIN=0
	INTEGER, SAVE :: LUOUT=0
!
	INTEGER GET_INDX_SP
	EXTERNAL GET_INDX_SP
!
	LOGICAL, SAVE :: USE_MILLI_ANG=.TRUE.
	LOGICAL END_FILE
	LOGICAL FILE_PRES
	CHARACTER(LEN=80) FILE_WITH_LINE_LIMS
	CHARACTER(LEN=80) OUT_FILE
	CHARACTER(LEN=80) STRING
!
	XLOC=0.5D0*(XRANGE(1)+XRANGE(2))
	YLOC=0.5D0*(YRANGE(1)+YRANGE(2))
	dY=0.05*(YRANGE(2)-YRANGE(1))
!
	IF(NPLTS .EQ. 1)THEN
	  IP=1
	ELSE
	  CALL GEN_IN(IP,'Plot for fitting')
	END IF
	USE_CURSOR=.TRUE.
	CALL GEN_IN(USE_CURSOR,'Use cursor (T) or read locations from file (F)')
!
! Used if average data on X-limits to defined the continuum level.
!
	NPIX=1
	CALL GEN_IN(NPIX,'Number of pixels at X location to average continuum (must be odd)')
	IF(MOD(NPIX,2) .EQ. 0)THEN
	  NPIX=NPIX+1
	  WRITE(6,*)'NPIX increase by 1 to make odd; NPIX=',NPIX
	END IF
	CALL GEN_IN(USE_MILLI_ANG,'Output EWs in milli-Angstroms?')
!
	IF(USE_CURSOR)THEN
!
	  WRITE(6,'(A)')' '
	  WRITE(6,'(A)')' Program assumes X axis increases with X index '
	  WRITE(6,'(A)')' '
	  WRITE(6,'(A)')' Use cursor to define left and right side of line:'
	  WRITE(6,'(A)')'   Y - uses cursor Y value for continuum location'
	  WRITE(6,'(A)')'   S - reset and start line selection again'
	  WRITE(6,'(A)')'   A - Y value is average over data at X location'
	  WRITE(6,'(A)')'   E - exit line selection'
	  WRITE(6,'(A)')' '
!
	ELSE
	  IF(LUIN .EQ. 0)CALL GET_LU(LUIN,'Input file in DO_CURSOR_EW')
	  FILE_WITH_LINE_LIMS='FILE_WITH_LINE_LIMS'
	  CALL GEN_IN(FILE_WITH_LINE_LIMS,'File with XSTART, XEND defining lines')
	  OPEN(UNIT=LUIN,FILE=FILE_WITH_LINE_LIMS,STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(6,*)'Error the follwing file cannot be opened',TRIM(FILE_WITH_LINE_LIMS)
	    RETURN
	  END IF
	END IF
!
! Open output file. Data is appended if it alread exists.
!
	IF(LUOUT .EQ. 0)THEN
	  OUT_FILE='EW_DATA'
	  IF(LUOUT .EQ. 0)CALL GET_LU(LUOUT,'Output file in DO_CURSOR_EW')
	  CALL GEN_IN(OUT_FILE,'File to OUTPUT EWs etc')
	  INQUIRE(FILE=OUT_FILE,EXIST=FILE_PRES)
	  IF(FILE_PRES)THEN
	    WRITE(6,*)'File already exists -- appending new data'
	    OPEN(UNIT=LUOUT,FILE=OUT_FILE,STATUS='OLD',ACTION='WRITE',POSITION='APPEND')
	  ELSE
	    OPEN(UNIT=LUOUT,FILE=OUT_FILE,STATUS='NEW',ACTION='WRITE')
	  END IF
	  IF(USE_MILLI_ANG)THEN
	    WRITE(LUOUT,'(5(6X,A),20X,3A)')'!    XST','    XEND','Line Loc',' F(cont)',' EW(mA)',
	1                          '      FWHM',' Sig(km/s)','    Sig(A)'
	  ELSE
	    WRITE(LUOUT,'(5(6X,A),20X,3A)')'!    XST','    XEND','Line Loc',' F(cont)','  EW(A)',
	1                          '      FWHM',' Sig(km/s)','    Sig(A)'
	  END IF
	END IF
!
	WRITE(6,'(A)')' '
	IF(USE_MILLI_ANG)THEN
	  WRITE( 6,'(5(6X,A),20X,3A)')'     XST','    XEND','Line Loc',' F(cont)',' EW(mA)',
	1                          '      FWHM',' Sig(km/s)','    Sig(A)'
	ELSE
	  WRITE( 6,'(5(6X,A),20X,3A)')'     XST','    XEND','Line Loc',' F(cont)','  EW(A)',
	1                          '      FWHM',' Sig(km/s)','    Sig(A)'
	END IF
	WRITE(6,'(A)')' '
!
!	CALL PGQLW(CUR_LW)
!	J=8*CUR_LW; CALL PGSLW(J)
	DO WHILE(1 .EQ. 1)
	  IF(USE_CURSOR)THEN
!
	    CALL PGQCI(IPEN)
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
	    IF(XST .LT. XRANGE(1))THEN
	      WRITE(6,*)'Invalid X start location -- start again'
	      GOTO 1000
	    END IF
	    IF(CURSVAL .EQ. 'E' .OR. CURSVAL .EQ. 'e')THEN
	      EXIT
	    ELSE IF(CURSVAL .EQ. 'S' .OR. CURSVAL .EQ. 's')THEN
	      GOTO 1000
	    ELSE IF(CURSVAL .EQ. 'Y' .OR. CURSVAL .EQ. 'y')THEN
	      CALL PGPT(IONE,XLOC,YLOC,ITWO)
	    ELSE
	      YST=0.0
	      DO J=IST-NPIX/2,IST+NPIX/2
	        YST=YST+CD(IP)%DATA(J)
	      END DO
	      YST=YST/NPIX
	      YT(1)=YST-dY; YT(2)=YST+dY
	      XT(1)=XLOC; XT(2)=XLOC
	      CALL PGLINE(2,XT,YT)
	    END IF
!
! Get end of line.
!
	    CURSERR = PGCURS(XLOC,YLOC,CURSVAL)
	    IF(XLOC .GT. XRANGE(2))THEN
	      WRITE(6,*)'Invalid X end location -- start again'
	      GOTO 1000
	    END IF
	    IEND=GET_INDX_SP(XLOC,CD(IP)%XVEC,NPTS(IP))
	    XEND=XLOC; YEND=YLOC
	    IF(CURSVAL .EQ. 'E' .OR. CURSVAL .EQ. 'e')THEN
	      EXIT
	    ELSE IF(CURSVAL .EQ. 'S' .OR. CURSVAL .EQ. 's')THEN
	      GOTO 1000
	    ELSE IF(CURSVAL .EQ. 'Y' .OR. CURSVAL .EQ. 'y')THEN
	      CALL PGPT(IONE,XLOC,YLOC,ITWO)
	    ELSE
	      YEND=0.0
	      DO J=IEND-NPIX/2,IEND+NPIX/2
	        YEND=YEND+CD(IP)%DATA(J)
	      END DO
	      YEND=YEND/NPIX
	      YT(1)=YEND-dY; YT(2)=YEND+dY
	      XT(1)=XLOC; XT(2)=XLOC
	      CALL PGLINE(2,XT,YT)
	    END IF
!
	  ELSE
!
! File input
!
	    STRING(1:1)='!'
	    DO WHILE(STRING(1:1) .EQ. '!' .OR. STRING .EQ. ' ')
	      END_FILE=.TRUE.
	      READ(LUIN,'(A)',END=2000)STRING
	      END_FILE=.FALSE.
	    END DO
	    READ(STRING,*)XST,XEND
!
! Get continuum fluxes at line limits.
!
	    IST=GET_INDX_SP(XST,CD(IP)%XVEC,NPTS(IP))
	    YST=0.0
	    DO J=IST-NPIX/2,IST+NPIX/2
	      YST=YST+CD(IP)%DATA(J)
	    END DO
	    YST=YST/NPIX
!
	    IEND=GET_INDX_SP(XEND,CD(IP)%XVEC,NPTS(IP))
	    YEND=0.0
	    DO J=IEND-NPIX/2,IEND+NPIX/2
	      YEND=YEND+CD(IP)%DATA(J)
	    END DO
	    YEND=YEND/NPIX
!
	  END IF
!
	  IF(IST .EQ. IEND)THEN
	    WRITE(6,*)'Error -- limits on line are identical'
	    WRITE(6,*)'Need to redefine the line'
	    GOTO 1000
	  END IF
!
! We now determine the line parameters
!
	  EW=0.0;  XMEAN=0.0
	  SLOPE=(YEND-YST)/(CD(IP)%XVEC(IEND)-CD(IP)%XVEC(IST))
	  DO I=IST,IEND
	    YVAL=YST+(CD(IP)%XVEC(I)-CD(IP)%XVEC(IST))*SLOPE
	    dX=(CD(IP)%XVEC(MIN(I+1,IEND))-CD(IP)%XVEC(MAX(IST,I-1)))/2
	    EW=EW+dX*(YVAL-CD(IP)%DATA(I))
	    EWL=EW+dX*(0.99*YVAL-CD(IP)%DATA(I))
	    EWH=EW+dX*(1.01*YVAL-CD(IP)%DATA(I))
	    XMEAN=XMEAN+CD(IP)%XVEC(I)*(YVAL-CD(IP)%DATA(I))*dX
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
	  YMEAN=EW/(CD(IP)%XVEC(IEND)-CD(IP)%XVEC(IST))
	  EW=EW/YCONT
	  EWL=EWL/(YCONT*0.99)
	  EWH=EWH/(YCONT*1.01)
!
! These parameters are used to provide infomation on whether a line is blended.
!
	  SIGMA=0.0D0; SKEWNESS=0.0D0; KURTOSIS=0.0D0
	  DO I=IST,IEND
	    XVAL=CD(IP)%XVEC(I)
	    YVAL=YST+(CD(IP)%XVEC(I)-CD(IP)%XVEC(IST))*SLOPE - CD(IP)%DATA(I)
	    dX=(CD(IP)%XVEC(MIN(I+1,IEND))-CD(IP)%XVEC(MAX(IST,I-1)))/2
	    SIGMA=SIGMA+dX*YVAL*(XVAL-XMEAN)**2
	    SKEWNESS=SKEWNESS+dX*YVAL*(XVAL-XMEAN)**3
	    KURTOSIS=KURTOSIS+dX*YVAL*(XVAL-XMEAN)**4
	  END DO
	  SIGMA=SQRT(SIGMA/YINT)
	  SKEWNESS=SKEWNESS/YINT/SIGMA**3
	  KURTOSIS=KURTOSIS/YINT/SIGMA**4
!
! The FWHM assume the line is Gaussian in shape.
!
	  T1=2.998D+05*SIGMA/XMEAN
	  IF(USE_MILLI_ANG)THEN
	    WRITE(6,'(3F14.3,ES14.4,F14.2,7F10.2)')XST,XEND,XMEAN,YCONT,1000.0*EW,
	1                     100.0D0*(EW-EWL)/EW, 100.0D0*(EW-EWH)/EW,
	1                     2.355*T1,T1,SIGMA,SKEWNESS,KURTOSIS
	    WRITE(LUOUT,'(3F14.3,ES14.4,F14.2,7F10.2)')XST,XEND,XMEAN,YCONT,1000.0*EW,
	1                     100.0D0*(EW-EWL)/EW, 100.0D0*(EW-EWH)/EW,
	1                     2.355*T1,T1,SIGMA,SKEWNESS,KURTOSIS
	  ELSE
	    WRITE(6,'(3F14.3,2ES14.4,7F10.2)')XST,XEND,XMEAN,YCONT,EW,
	1                     100.0D0*(EW-EWL)/EW, 100.0D0*(EW-EWH)/EW,
	1                     2.355*T1,T1,SIGMA,SKEWNESS,KURTOSIS
	    WRITE(LUOUT,'(3F14.3,2ES14.4,7F10.2)')XST,XEND,XMEAN,YCONT,EW,
	1                     100.0D0*(EW-EWL)/EW, 100.0D0*(EW-EWH)/EW,
	1                     2.355*T1,T1,SIGMA,SKEWNESS,KURTOSIS
	  END IF
	  FLUSH(LUOUT)
!
1000	  CONTINUE
	END DO
2000	CONTINUE
!	CALL PGSLW(CUR_LW)
!
	RETURN
	END
