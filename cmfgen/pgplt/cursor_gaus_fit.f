!
! Subroutine to fits a straigt line, and a set of modified Gaussians,
! to a set of data. Routine is designed to be called in GRAMON_PGPLOT.
!
	SUBROUTINE CURSOR_GAUS_FIT(
	1             XPAR,XINC,XNUMST,IXTICK,IDX,
	1             YPAR,YINC,YNUMST,IYTICK,IDY,
	1             TICK_FAC,EXPCHAR,
	1             XLABEL,YLABEL,TITONRHS,
	1             LOG_AXIS,OPTION,NORMAL_R_Y_AXIS,
	1             XLAB_FILE,YLAB_FILE,
	1             PEN_COL,PEN_OFFSET,REVERSE_PLOT_ORDER)
	USE SET_KIND_MODULE
	USE MOD_CURVE_DATA
	USE GAUS_FIT_DATA
	USE GEN_IN_INTERFACE
	USE MOD_COLOR_PEN_DEF
	IMPLICIT NONE
!
! Altered 27-Feb-2022: Fixed bug in EW error estimates when using unormalized data.
! Altered 25-Feb-2022: Some cleaning
! Altered 01-Feb-2022: Extensive changes implmented to allow cursor input of
!                        parameters for gauss fitting.
! Altered   -Sep-07
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
	CHARACTER(LEN=80) TRANS_NAME
	CHARACTER(LEN=80) XLAB_FILE,YLAB_FILE
!
        INTEGER PEN_COL(0:MAX_PLTS)
        INTEGER PEN_OFFSET
        LOGICAL REVERSE_PLOT_ORDER
!
	INTEGER, PARAMETER :: IMONE=-1
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: ITWO=2
	INTEGER, PARAMETER :: ISIX=6
	INTEGER, PARAMETER :: NL_MAX=100
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
	REAL(KIND=LDP) LINE_CENTER(NL_MAX)
!
	REAL*4, SAVE :: FIND_LINE_TOLERANCE=0.002
!
	CHARACTER(LEN=30) UC
	REAL(KIND=LDP) TOLERANCE
	REAL(KIND=LDP) GAUS_FIT_FUNC
	EXTERNAL GAUS_FIT_FUNC,UC
!
	REAL(KIND=LDP) YST,YEND
!
	REAL(KIND=LDP) TOL
	REAL(KIND=LDP) T1,T2
	REAL(KIND=LDP) FWHM
	REAL(KIND=LDP) CONT_FLUX
	INTEGER NO_LINES
	INTEGER ITER
!
	REAL*4 SPACING
	EXTERNAL SPACING
!
	LOGICAL USE_CURSOR
	REAL*4 XVAL,YVAL
	REAL*4 DEF_SIGMA
	REAL*4 SIGMA_THIS_LINE
!
	REAL*4 SP_EW,SP_LOC,SP_FWHM,SP_LINE_CENTER
!
        INTEGER PGCURS
        INTEGER CURSERR
	CHARACTER(LEN=1) CURSVAL
!
	REAL(KIND=LDP), SAVE :: XST=0.0D0
	REAL(KIND=LDP), SAVE :: XEND=1.0D0
	REAL*4 SP_XST,SP_XEND
	INTEGER, SAVE :: IP=1
	INTEGER, SAVE :: NG_PAR_OLD
	INTEGER, SAVE :: LUGF,LUEW
        LOGICAL, PARAMETER :: L_TRUE=.TRUE.
	LOGICAL, SAVE :: FIRST_WRITE=.TRUE.
	LOGICAL NEW_FILE
	LOGICAL FILE_EXISTS
	LOGICAL DO_FIND
	LOGICAL EDIT_FIT
!
	INTEGER I,J,K
	INTEGER LW_SAVE
	INTEGER MARKER_LW
	LOGICAL GUESSED
	LOGICAL NEW_REGION
	LOGICAL SHOW_INIT_GAUS
!
	NG_PAR_MAX=4*NL_MAX+2
	IF(.NOT. ALLOCATED(PAR))THEN
	  ALLOCATE(PAR(NG_PAR_MAX))
	END IF
!
! Set plot and defaults.
!
	IP=1; DEF_SIGMA=0.1; SHOW_INIT_GAUS=.FALSE.
	CALL GEN_IN(IP,'Plot number for Gaussian fitting')
	CALL GEN_IN(DEF_SIGMA,'Default sigma for fitting line profiles')
	CALL GEN_IN(SHOW_INIT_GAUS,'Show initial Gaussian estimates on the plot?')
!
	CALL PRINT_CURSOR_DESC()
!
	MARKER_LW=4
	CALL PGQLW(LW_SAVE)
	SIGMA_THIS_LINE=0
	dY=0.05*(YPAR(2)-YPAR(1))	
!
	CALL PGSLW(MARKER_LW)
	NUM_GAUS=0
	J=1; CALL PGSCI(J)
	XVAL=0.5*(XPAR(1)+XPAR(2))
	YVAL=0.5*(YPAR(1)+YPAR(2))
!
	DO WHILE(1 .EQ. 1)
	  CURSERR = PGCURS(XVAL,YVAL,CURSVAL)
	  IF(UC(CURSVAL) .EQ. 'B')THEN
	    IF(NUM_GAUS .NE. 0)THEN
	      WRITE(6,*)'Final SUM_SQ=',SUM_SQ(1)
	      CALL WRITE_GAUS_FIT()
	      NUM_GAUS=0
	    END IF
	    YT(1)=YVAL-dY; YT(2)=YVAL+dY
	    XT(1)=XVAL; XT(2)=XVAL
	    CALL PGLINE(2,XT,YT)
	    XST=XVAL; YST=YVAL
	    CURSERR = PGCURS(XVAL,YVAL,CURSVAL)
	    YT(1)=YVAL-dY; YT(2)=YVAL+dY
	    XT(1)=XVAL; XT(2)=XVAL
	    CALL PGLINE(2,XT,YT)
	    XEND=XVAL; YEND=YVAL
	    PAR(1)=YST
	    PAR(2)=(YEND-YST)/(XEND-XST)
!
	  ELSE IF(UC(CURSVAL) .EQ. 'D')THEN
	    T1=XVAL
	    CALL PGPT(IONE,XVAL,YVAL,IMONE)
	    CURSERR = PGCURS(XVAL,YVAL,CURSVAL)
	    CALL PGPT(IONE,XVAL,YVAL,IMONE)
	    DEF_SIGMA=ABS(XVAL-T1)/2.2
!
	  ELSE IF(UC(CURSVAL) .EQ. 'S')THEN
	    CALL PGPT(IONE,XVAL,YVAL,IMONE)
	    T1=XVAL
	    CURSERR = PGCURS(XVAL,YVAL,CURSVAL)
	    CALL PGPT(IONE,XVAL,YVAL,IMONE)
	    T1=XVAL
	    SIGMA_THIS_LINE=ABS(XVAL-T1)/2.2

	  ELSE IF(UC(CURSVAL) .EQ. 'L')THEN
	    CALL PGPT(IONE,XVAL,YVAL,IMONE)
	    T1=XVAL
	    NUM_GAUS=NUM_GAUS+1
	    J=NUM_GAUS
	    K=2+(J-1)*4+1
	    PAR(K)=XVAL; PAR(K+2)=YVAL-(PAR(1)+PAR(2)*(XVAL-XST))
	    PAR(K+3)=2.0
	    PAR(K+1)=DEF_SIGMA
	    IF(SIGMA_THIS_LINE .NE. 0)PAR(K+1)=SIGMA_THIS_LINE
	    SIGMA_THIS_LINE=0
!
! Options to move/expand/contract plot in X direction.
!
	  ELSE IF(UC(CURSVAL) .EQ. 'N' .OR.
	1           UC(CURSVAL) .EQ. 'P' .OR.
	1           UC(CURSVAL) .EQ. 'R' .OR.
	1           UC(CURSVAL) .EQ. 'X' .OR.
	1           UC(CURSVAL) .EQ. 'Y' .OR.
	1           UC(CURSVAL) .EQ. 'C')THEN
	    CALL PGSLW(LW_SAVE)
	    T1=XPAR(2)-XPAR(1)
	    IF(UC(CURSVAL) .EQ. 'N')THEN
	      XPAR(1)=XPAR(2)-XINC
	      XPAR(2)=XPAR(1)+T1
	    ELSE IF(UC(CURSVAL) .EQ. 'P')THEN
	      XPAR(1)=XPAR(1)+XINC-(XPAR(2)-XPAR(1))
	      XPAR(2)=XPAR(1)+T1
	    ELSE IF(UC(CURSVAL) .EQ. 'R')THEN
	    ELSE IF(UC(CURSVAL) .EQ. 'X')THEN
	      XPAR(2)=XPAR(2)+XINC
	    ELSE IF(UC(CURSVAL) .EQ. 'C')THEN
	      XPAR(2)=XPAR(2)-XINC
	    ELSE IF(UC(CURSVAL) .EQ. 'Y')THEN
	      YPAR(1)=1.0D+30; YPAR(2)=-1.0D+30
	      DO I=1,NPTS(IP)
	        IF( (CD(IP)%XVEC(I)-XPAR(1))*(CD(IP)%XVEC(I)-XPAR(2)) .LE. 0)THEN
	          YPAR(1)=MIN(YPAR(1),CD(IP)%DATA(I))
	          YPAR(2)=MAX(YPAR(2),CD(IP)%DATA(I))
	        END IF
	      END DO
	      IF(YPAR(1) .LT. 0)YPAR(1)=0
	      YPAR(2)=1.1*YPAR(2)
	      YINC=SPACING(YPAR(1),YPAR(2))
	      dY=0.05*(YPAR(2)-YPAR(1))	
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
	    XVAL=0.5D0*(XPAR(1)+XPAR(2))
	    CALL PGSCH(XVAL)
	    YVAL=1.0
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
	    CALL PGSLW(MARKER_LW)
!
	  ELSE IF(UC(CURSVAL) .EQ. 'E')THEN
	    NG_PAR=2+4*NUM_GAUS
	    CALL ED_GAUS_FIT()
	    NG_PAR=2+4*NUM_GAUS
	    NG_PAR_OLD=NG_PAR
	    WRITE(6,*)'Exited ED_GAUS_FIT: Use F optio to redo fit'
!
	  ELSE IF(UC(CURSVAL) .EQ. 'F')THEN
	    IF(NUM_GAUS .EQ. 0)THEN
	      WRITE(6,*)'Error -- no lines selected'
	      WRITE(6,*)'Select lines, and if necssary, redefine the fitting band'
	    ELSE
!
! We can now determine the optimal fitting parameters.
!
	      NG_PAR=2+4*NUM_GAUS
!
! Draw initial Gaussians using a dashed black pen.
!
	      IF(SHOW_INIT_GAUS)THEN
	        I=2; CALL PGSLS(I)			!Dashed curve
	        CALL SET_GAUS_DATA(CD(IP)%XVEC,CD(IP)%DATA,XST,XEND,NPTS(IP),YST,YEND)
	        CALL DRAW_GAUS_V2(L_FALSE)
	        CALL PGEBUF				!Clear buffer
	        I=1; CALL PGSLS(I)			!Reset to solid lines
	      END IF
!
	      IF(ALLOCATED(SIM))DEALLOCATE(SIM,SUM_SQ,SCALE,EW,EW_CONT)
	      ALLOCATE (SIM(NG_PAR+1,NG_PAR))
	      ALLOCATE (SUM_SQ(NG_PAR+1))
	      ALLOCATE (SCALE(NG_PAR+1))
	      ALLOCATE (EW(NUM_GAUS))
	      ALLOCATE (EW_CONT(NUM_GAUS))
	      WRITE(6,*)'Done allocate', NUM_GAUS, NG_PAR; FLUSH(UNIT=6)
!
	      WRITE(6,*)CD(IP)%XVEC(1),CD(IP)%XVEC(NPTS(IP))
	      WRITE(6,*)CD(IP)%DATA(1),CD(IP)%DATA(NPTS(IP))
	      WRITE(6,*)'Call set gaus_data'; FLUSH(UNIT=6)
	      CALL SET_GAUS_DATA(CD(IP)%XVEC,CD(IP)%DATA,XST,XEND,NPTS(IP),YST,YEND)
	      WRITE(6,*)'Called set gaus_data'; FLUSH(UNIT=6)
!
	      SIM(1,1:NG_PAR)=PAR(1:NG_PAR)
!
! Set other vertices of simplex. To do this we frist set the
! the characteristic scaleover which the variable may very.
!
	      SCALE(1)=0.002			!Mean level
	      SCALE(2)=0.02/(XEND-XST)	!Slope
	      DO I=3,NG_PAR,4
	        SCALE(I)=0.2			!Wavelength
	        SCALE(I+1)=0.10		!Sigma
	        SCALE(I+2)=0.05		!Height
	        SCALE(I+3)=0.2		!Exponent
	      END DO
!
	      DO J=2,NG_PAR+1
	        DO I=1,NG_PAR
	          SIM(J,I)=SIM(1,I)
	        END DO
	      END DO
	      DO J=2,NG_PAR+1
	        SIM(J,J-1)=SIM(J,J-1)+SCALE(J-1)		!ws scale(j)
	      END DO
!
	      WRITE(6,*)'Simples vertices have been set'
	      DO J=1,NG_PAR
	        WRITE(6,'(20ES12.4)')(SIM(I,J),I=1,NG_PAR+1)
	      END DO
!
              SUM_SQ(:)=0.0D0
              DO I=1,NG_PAR+1
                PAR(1:NG_PAR)=SIM(I,1:NG_PAR)
                SUM_SQ(I)=GAUS_FIT_FUNC(PAR)
              END DO
	      WRITE(6,*)'Evaluated SUM_SQ: Will now do the fit'
!
              TOL=1.0D-08
	      I=NG_PAR+1
              CALL AMOEBA(SIM,SUM_SQ,I,NG_PAR,NG_PAR,TOL,GAUS_FIT_FUNC,ITER)
!
! Compute EW of profile. The order of passing the variables is HEIGHT,
! SIGMA, EXPONENT. The location is not needed.
!
	      TOLERANCE=1.0D-05
	      DO I=1,NUM_GAUS
	        K=2+(I-1)*4+1
	        T1=X_GAUS(1)
	        CALL GAUS_ROMB_V2(EW(I),SIM(1,1),SIM(1,2),T1,SIM(1,K),
	1                              SIM(1,K+2),SIM(1,K+1),SIM(1,K+3),TOLERANCE)
!	        T1=SIM(1,1)+SIM(1,2)*(SIM(I,K)-X_GAUS(1))				!Continuum value
!	        EW(I)=EW(I)/T1
	        EW_CONT(I)=T1
	      END DO
	      EW=EW*1000.0D0			!mAng
              PAR(1:NG_PAR)=SIM(1,1:NG_PAR)
	      CALL GAUS_FIT_ER(PAR)
!
	      DO I=1,NUM_GAUS
	        T1=1000.0D0/EW_CONT(I)
	        EW_ERROR(I)=EW_ERROR(I)*T1
	        ALT_ERROR(I)=ALT_ERROR(I)*T1
	        MIN_ERROR(I)=MIN_ERROR(I)*1000.0D0
	      END DO
!
	      WRITE(6,*)'Called AMOEBA'
	      WRITE(6,*)'Fit parameters are (NB SIGMA is not Stan. Dev.):'
	      WRITE(6,*)SIM(1,1),SIM(1,2)
	      WRITE(6,'(2X,(5X,A),5(3X,A),2(3X,A))')
	1          '        Lam','     Height','          a',
	1          '        EXP','Sigma(km/s)',' FWHM(km/s)',' EW(mA)','Err(mA)'
	      DO I=1,NUM_GAUS
	        K=2+(I-1)*4+1
	        T1=2.99794D+05*PAR(K+1)/PAR(K)
	        FWHM=2.0D0*T1*(LOG(2.0D0))**(1.0D0/PAR(K+3))
	        WRITE(6,'(2X,ES16.6,5ES14.4,4F10.2)')PAR(K),PAR(K+2),PAR(K+1),PAR(K+3),
	1                   T1/SQRT(2.0D0),FWHM,EW(I),EW_ERROR(I),ALT_ERROR(I),MIN_ERROR(I)
	      END DO
!
              SUM_SQ(:)=0.0D0
              SUM_SQ(1)=GAUS_FIT_FUNC(PAR)
!
! Draw Gaussians using a black pen.
!
	      CALL PGSLW(LW_SAVE)
	      I=1
	      CALL PGSCI(I)
	      CALL PGLINE(NG_DATA,XFIT,YFIT)
	      CALL PGEBUF				!Clear buffer
	      CALL PGSLW(MARKER_LW)
	    END IF
!
	  ELSE IF(UC(CURSVAL) .EQ. 'Q')THEN
	    IF(NUM_GAUS .NE. 0)CALL WRITE_GAUS_FIT()
	    EXIT
!
	  ELSE
	    WRITE(6,*)RED_PEN,' Cursor value not recognized'	
	    CALL PRINT_CURSOR_DESC
	  END IF
!
	END DO
!
	RETURN

	CONTAINS
!
! We output lines in wavelength order. After calling, INDEXX,
! INDX_VEC(1) will contain the first line to be output, etc.
! We don't actually perform the sort.
!
	SUBROUTINE WRITE_GAUS_FIT
	USE SET_KIND_MODULE
!
	IF(FIRST_WRITE)THEN
	  CALL GET_LU(LUGF,'In CURSOR_GAUS_EW')
	  FIRST_WRITE=.FALSE.
	  INQUIRE(FILE='GAUSS_FITS',EXIST=NEW_FILE); NEW_FILE=.NOT. NEW_FILE
	  IF(.NOT. NEW_FILE)THEN
	    WRITE(6,*)'If file GAUS_FITS exists: Default is to append new fits'
	    CALL GEN_IN(NEW_FILE,'Open new file (T) or append data (F)')
	  END IF
	  IF(NEW_FILE)THEN
	    OPEN(UNIT=LUGF,FILE='GAUSS_FITS',STATUS='UNKNOWN')
	    WRITE(LUGF,'(3X,A,1X,4(10X,A))')
	1          'NG',' XST','XEND','   A','   B    [as in A+B(x-xst)]'
	    WRITE(LUGF,'(6X,(3X,A),5(3X,A),2(3X,A))')
	1          '        Lam','     Height','          a',
	1          '        EXP','Sigma(km/s)',' FWHM(km/s)',' EW(mA)','Err(mA)'
	  ELSE
	    OPEN(UNIT=LUGF,FILE='GAUSS_FITS',STATUS='UNKNOWN',POSITION='APPEND')
	  END IF
!
	  CALL GET_LU(LUEW,'In CURSOR_GAUS_EW --  LUEW')
	  INQUIRE(FILE='GAUSS_FIT_EWS',EXIST=NEW_FILE); NEW_FILE=.NOT. NEW_FILE
	  IF(.NOT. NEW_FILE)THEN
	    WRITE(6,*)'If file GAUS_FIT)EW exists: Default is to append new fits'
	    CALL GEN_IN(NEW_FILE,'Open new file (T) or append data (F)')
	  END IF
	  IF(NEW_FILE)THEN
	    OPEN(UNIT=LUEW,FILE='GAUSS_FIT_EW',STATUS='UNKNOWN')
	    WRITE(LUEW,'(4X,A,5(4X,A),2X,A,4X,A)')
	1          'IP','Line Lam','Line Loc',' F(cont)','  EW(mA)',
	1          '   E(mA)','FWHM(km/s)','dLAM(km/s)','Transition'
	  ELSE
	    OPEN(UNIT=LUEW,FILE='GAUSS_FIT_EW',STATUS='UNKNOWN',POSITION='APPEND')
	  END IF
	ELSE
!	  OPEN(UNIT=LUGF,FILE='GAUSS_FITS',STATUS='UNKNOWN',POSITION='APPEND')
!	  OPEN(UNIT=LUEW,FILE='GAUSS_FIT_EW',STATUS='UNKNOWN',POSITION='APPEND')
	END IF
!
	IF(ALLOCATED(INDX_VEC))DEALLOCATE(INDX_VEC)
	ALLOCATE (INDX_VEC(NUM_GAUS))
        DO I=1,NUM_GAUS
	   K=3+(I-1)*4
	   LINE_CENTER(I)=PAR(K)
	END DO
	IF(NUM_GAUS .EQ. 1)THEN
	   INDX_VEC(1)=1
	ELSE
	   CALL INDEXX(NUM_GAUS,LINE_CENTER,INDX_VEC,L_TRUE)
	END IF
!
	WRITE(LUGF,'(A)')' '
	WRITE(LUGF,'(2X,I3,1X,4ES14.6)')NUM_GAUS,XST,XEND,PAR(1:2)
	DO J=1,NUM_GAUS
	   I=INDX_VEC(J)
	   K=2+(I-1)*4+1
	   T1=2.99794D+05*SIM(1,K+1)/SIM(1,K)
	   FWHM=2.0D0*T1*(LOG(2.0D0))**(1.0D0/SIM(I,K+3))
	   SP_EW=EW(I); SP_LOC=SIM(I,K); SP_FWHM=FWHM
	   CALL GET_LINE_ID_PG(TRANS_NAME,SP_LINE_CENTER,SP_EW,SP_LOC,SP_FWHM)
	   WRITE(LUGF,'(6X,ES14.6,5ES14.4,4F10.2)')SIM(I,K),SIM(1,K+2),SIM(1,K+1),SIM(I,K+3),
	1               T1/SQRT(2.0D0),FWHM,EW(I),EW_ERROR(I),ALT_ERROR(I),MIN_ERROR(I)
	   CONT_FLUX=SIM(1,1)+SIM(1,2)*(SIM(I,K)-X_GAUS(1))
	   IF(TRANS_NAME .EQ. ' ')THEN
	     WRITE(LUEW,'(2X,I4,2F12.3,F12.4,3F12.2,11A,A,4X,A)')IP,SP_LINE_CENTER,SP_LOC,CONT_FLUX,
	1                    EW(I),EW_ERROR(I),FWHM,'-'
	   ELSE
	     T2=2.998D+05*ABS(SP_LINE_CENTER-SP_LOC)/SP_LINE_CENTER
	     WRITE(LUEW,'(2X,I4,2F12.3,F12.4,4F12.2,4X,A)')IP,SP_LINE_CENTER,SP_LOC,CONT_FLUX,
	1                    EW(I),EW_ERROR(I),FWHM,T2,TRIM(TRANS_NAME)
	   END IF
	END DO
	FLUSH(UNIT=LUEW)
	END SUBROUTINE WRITE_GAUS_FIT
!
	SUBROUTINE PRINT_CURSOR_DESC
	USE SET_KIND_MODULE
!
! We restrict options to lower case to avoid accidental use of
! the mouse which returns A.
!
	WRITE(6,'(A)')BLUE_PEN
	WRITE(6,'(A)')' Program assumes X axis increases with X index '
	WRITE(6,'(A)')' '
	WRITE(6,'(A)')' Two inputs are required to define the input band continuum'
	WRITE(6,'(A)')' Two inputs are required to define sigma'
	WRITE(6,'(A)')' '
	WRITE(6,'(A)')' The following keys should be used'
	WRITE(6,'(A)')' '
	WRITE(6,'(A)')' B  -  set input band (require two inputs)'
	WRITE(6,'(A)')' L  -  indicate line center and height/depth'
	WRITE(6,'(A)')' D  -  use to set def sigma (use line mid-points)'
	WRITE(6,'(A)')' S  -  use to  sigma estimate for next line only (use line mid points)'
	WRITE(6,'(A)')' E  -  edit current line fit (setting NUM_GAUSS to 0 restarts process)'
	WRITE(6,'(A)')' F  -  finish line selection and do fit'
	WRITE(6,'(A)')' Q  -  quit selection process'
!
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
!
! Routine to draw the fit on top of an existing plot. At present
! plot is drawn in black only. Routine can also draw the difference
! betwene the data and the fitted curve.
!
	SUBROUTINE DRAW_GAUS_V2(DIFF_ALSO)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
	LOGICAL DIFF_ALSO
	REAL(KIND=LDP) GAUS_FIT_FUNC
	EXTERNAL GAUS_FIT_FUNC
	INTEGER I
	REAL(KIND=LDP) T1
!
! PAR must be set on entry.
! The Gaussians are drwan in black.
!
        T1=GAUS_FIT_FUNC(PAR)
	I=1; CALL PGSCI(I)
	CALL PGLINE(NG_DATA,XFIT,YFIT)
	IF(DIFF_ALSO)THEN
	  I=10; CALL PGSCI(I)
	  FIT_DIF=Y_GAUS-YFIT
	  CALL PGLINE(NG_DATA,XFIT,FIT_DIF)
	END IF
	RETURN
	END SUBROUTINE DRAW_GAUS_V2
	END SUBROUTINE CURSOR_GAUS_FIT
