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
	SUBROUTINE DO_CURSOR_BALMER(
	1             XPAR,XINC,XNUMST,IXTICK,IDX,
	1             YPAR,YINC,YNUMST,IYTICK,IDY,
	1             TICK_FAC,EXPCHAR,
	1             XLABEL,YLABEL,TITONRHS,
	1             LOG_AXIS,OPTION,NORMAL_R_Y_AXIS,
	1             XLAB_FILE,YLAB_FILE,
	1             PEN_COL,PEN_OFFSET,REVERSE_PLOT_ORDER)
	USE SET_KIND_MODULE
	USE MOD_CURVE_DATA
	USE GEN_IN_INTERFACE
	USE MOD_COLOR_PEN_DEF
	IMPLICIT NONE
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
	LOGICAL FIT_DONE
!
	CHARACTER(LEN=80) XLABEL,YLABEL
        CHARACTER(LEN=5)  LOG_AXIS
	CHARACTER(LEN=80) OPTION
	CHARACTER(LEN=80) XLAB_FILE,YLAB_FILE
!
	INTEGER, SAVE :: NPIX=3  	!Integraton band pass around line limits
	REAL*4,  SAVE :: LOW_CONT=0.998
	REAL*4,  SAVE :: HIGH_CONT=1.002
!
	REAL*4  OLAM_ST(20), OLAM_END(20)
	INTEGER OMIT_ST(20), OMIT_END(20)
!
	REAL*4  C_OLAM_ST(20), C_OLAM_END(20)
	INTEGER C_OMIT_ST(20), C_OMIT_END(20)
!
	INTEGER PEN_COL(0:MAX_PLTS)
	INTEGER PEN_OFFSET
	LOGICAL REVERSE_PLOT_ORDER
!
	INTEGER IST,IEND		!Line limits in pixel space
	REAL*4 XLOC,YLOC		!Used to read cursor location
	REAL*4 YST,YEND  		!Continuum flux at line limits
	REAL*4 XST,XEND
!
! In the following I use Y and F interchangably.
! Also X will normally be Lambda.
!
	REAL*4 dX               !X spacing (pixel centered)
	REAL*4 XVAL  		!Current X value
	REAL*4 YVAL  		!Current Y value
	REAL*4 YINT
	REAL*4 T1,T2,T3,T4	!Work variables (single precision)
	REAL*4 X1,X2		!Work variable
	REAL*4 D1,D2		!Work variable
	REAL(KIND=LDP) LAM_CENT
	REAL(KIND=LDP) LAM_WAVE
	REAL(KIND=LDP) MEAN
!
	REAL(KIND=LDP) SUM_OSQ, SUM_LOSQ, SUM_LSQ_OSQ
	REAL(KIND=LDP) SUM_MO, SUM_LMO
	REAL(KIND=LDP) A, B
	REAL(KIND=LDP) DET, DA, DB
	REAL(KIND=LDP) CHISQ,RAW_CHISQ,RED_CHISQ
	REAL(KIND=LDP) EW_OBS,EW_MOD
	REAL(KIND=LDP) EW_OBS_OMIT,EW_MOD_OMIT
	REAL(KIND=LDP) SLOPE,INTER,CONT
!
	REAL*4 LAM
	REAL*4 O_DATA
	REAL*4, ALLOCATABLE :: OBS_DATA(:)
	REAL*4, ALLOCATABLE :: MOD_DATA(:)
	REAL*4, ALLOCATABLE :: MASK(:)
	REAL*4, ALLOCATABLE :: C_MASK(:)
	REAL*4, ALLOCATABLE :: CHISQ_CONTRIB(:)
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
	INTEGER, SAVE :: LU_OUT=0
!
	INTEGER IPEN_SAVE
	INTEGER N_OMIT
	INTEGER NC_OMIT
	INTEGER IP_MOD
	INTEGER IP_OBS
	INTEGER IP_OUT
	INTEGER I,J,K
	INTEGER NP
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
	CHARACTER(LEN=80) TRANS_NAME
!
! XLOC, YLOC will define initial cursor location.
!
	XLOC=0.5D0*(XPAR(1)+XPAR(2))
	YLOC=0.5D0*(YPAR(1)+YPAR(2))
	dY=0.05*(YPAR(2)-YPAR(1))
	IY_DIR=1
	FIT_DONE=.FALSE.
!
	IF(.NOT. RESET_DEFAULTS)THEN
	  CALL GEN_IN(RESET_DEFAULTS,'Reset plot and parameters/settings?')
	END IF
	IF(RESET_DEFAULTS)THEN
	  IF(NPLTS .LT. 2)THEN
	    WRITE(6,*)'Error - this program requires at least two plots'
	    IOS=1
	    RETURN
	  END IF
	  IP_MOD=1; IP_OBS=2
	  CALL GEN_IN(IP_MOD,'Model plot')
	  CALL GEN_IN(IP_OBS,'Observational plot ')
	END IF
!
! Interpolate data onto OBSERVATIONAL grid.
!
	IF(ALLOCATED(MOD_DATA))THEN
	  DEALLOCATE(MOD_DATA,OBS_DATA,MASK,C_MASK,CHISQ_CONTRIB)
	END IF
	NP=NPTS(IP_OBS)
	WRITE(6,*)'NP=',NP
	ALLOCATE(MOD_DATA(NP),MASK(NP),C_MASK(NP),OBS_DATA(NP),CHISQ_CONTRIB(NP))
	CALL MON_INTERP_SP(MOD_DATA,NP,IONE,CD(IP_OBS)%XVEC,NP,CD(IP_MOD)%DATA,NPTS(IP_MOD),CD(IP_MOD)%XVEC,NPTS(IP_MOD))
	OBS_DATA=CD(IP_OBS)%DATA
	CHISQ_CONTRIB=0.0D0
!
! Open output file. Data is appended if it alread exists.
!
	IF(LU_OUT .EQ. 0)THEN
	  OUT_FILE='BALMER_DATA'
	  CALL GET_LU(LU_OUT,'LU_OUT in DO_CURSOR_BALMER_V1')
	  CALL GEN_IN(OUT_FILE,'File to OUTPUT chi^2 and EW etc')
	  INQUIRE(FILE=OUT_FILE,EXIST=FILE_PRES)
	  IF(FILE_PRES)THEN
	    WRITE(6,*)'File already exists -- appending new data'
	    OPEN(UNIT=LU_OUT,FILE=OUT_FILE,STATUS='OLD',ACTION='WRITE',POSITION='APPEND')
	  ELSE
	    OPEN(UNIT=LU_OUT,FILE=OUT_FILE,STATUS='NEW',ACTION='WRITE')
	    CALL WRITE_BALMER_HEADER(LU_OUT)
	  END IF
	END IF
!
	IST=0; IEND=0
	CALL PRINT_CURSOR_DESC
	DO WHILE(1 .EQ. 1)
!
	  CURSERR = PGCURS(XLOC,YLOC,CURSVAL)
	  IF(UC(CURSVAL) .EQ. 'Q')THEN
	    EXIT
	  ELSE IF(UC(CURSVAL) .EQ. 'S')THEN
	    FIT_DONE=.FALSE.
	    GOTO 1000
	  ELSE IF(UC(CURSVAL) .EQ. 'B')THEN
!
! This defines the full extent of the line.
! Note -- we already have the first cursor location.
!
	    IST=GET_INDX_SP(XLOC,CD(IP_OBS)%XVEC,NPTS(IP_OBS));    XST=XLOC
            YT(1)=YLOC-2*dY; YT(2)=YLOC+2*dY
            XT(1)=XLOC; XT(2)=XLOC
            CALL PGLINE(2,XT,YT)
!
	    CURSERR = PGCURS(XLOC,YLOC,CURSVAL)
	    IEND=GET_INDX_SP(XLOC,CD(IP_OBS)%XVEC,NPTS(IP_OBS));   XEND=XLOC
            YT(1)=YLOC-2*dY; YT(2)=YLOC+2*dY
            XT(1)=XLOC; XT(2)=XLOC
            CALL PGLINE(2,XT,YT)
	    N_OMIT=0
	    NC_OMIT=0
!
	  ELSE IF(UC(CURSVAL) .EQ. 'Z')THEN
!
	    IF(IST .EQ. 0 .OR. IEND .EQ. 0)THEN
	      WRITE(6,*)'Line region is undefined -- use B to defined line region'
	      GOTO 1000
	    END IF
!
	    IPEN_SAVE=IPEN
	    IPEN=7
            CALL PGSCI(IPEN)
!
! Note -- we already have the first cursor location.
!
	    I=NC_OMIT+1; C_OLAM_ST(I)=XLOC
	    C_OMIT_ST(I)=GET_INDX_SP(XLOC,CD(IP_OBS)%XVEC,NPTS(IP_OBS))
            YT(1)=YLOC-dY; YT(2)=YLOC+dY
            XT(1)=XLOC; XT(2)=XLOC
            CALL PGLINE(2,XT,YT)
!
	    CURSERR = PGCURS(XLOC,YLOC,CURSVAL); C_OLAM_END(I)=XLOC
	    C_OMIT_END(I)=GET_INDX_SP(XLOC,CD(IP_OBS)%XVEC,NPTS(IP_OBS))
            YT(1)=YLOC-dY; YT(2)=YLOC+dY
            XT(1)=XLOC; XT(2)=XLOC
            CALL PGLINE(2,XT,YT)
	    NC_OMIT=NC_OMIT+1
	    IPEN=IPEN_SAVE
!
	  ELSE IF(UC(CURSVAL) .EQ. 'O')THEN
!
	    IF(IST .EQ. 0 .OR. IEND .EQ. 0)THEN
	      WRITE(6,*)'Line region is undefined -- use B to defined line region'
	      GOTO 1000
	    END IF
!
            IF(IPEN .EQ. 6)THEN ! Change pen colors for defining limits of the next line.
              IPEN=1
            ELSE
              IPEN=6
            END IF
            CALL PGSCI(IPEN)
!
! Note -- we already have the first cursor location.
!
	    I=N_OMIT+1; OLAM_ST(I)=XLOC
	    OMIT_ST(I)=GET_INDX_SP(XLOC,CD(IP_OBS)%XVEC,NPTS(IP_OBS))
            YT(1)=YLOC-dY; YT(2)=YLOC+dY
            XT(1)=XLOC; XT(2)=XLOC
            CALL PGLINE(2,XT,YT)
!
	    CURSERR = PGCURS(XLOC,YLOC,CURSVAL); OLAM_END(I)=XLOC
	    OMIT_END(I)=GET_INDX_SP(XLOC,CD(IP_OBS)%XVEC,NPTS(IP_OBS))
            YT(1)=YLOC-dY; YT(2)=YLOC+dY
            XT(1)=XLOC; XT(2)=XLOC
            CALL PGLINE(2,XT,YT)
	    N_OMIT=N_OMIT+1
!
	  ELSE IF(UC(CURSVAL) .EQ. 'F')THEN
!
	    IF(IST .EQ. 0 .OR. IEND .EQ. 0)THEN
	      WRITE(6,*)'Line region is undefined -- use B to defined line region'
	      GOTO 1000
	    END IF
!
	    NPIX=0
	    MASK=1.0D0; C_MASK=1.0D0
	    DO I=1,N_OMIT
	      MASK(OMIT_ST(I):OMIT_END(I))=0.0D0
	    END DO
	    DO I=1,NC_OMIT
	      C_MASK(C_OMIT_ST(I):C_OMIT_END(I))=0.0D0
	    END DO
!
	    SUM_OSQ=0.0D0;    SUM_LOSQ=0.0D0;   SUM_LSQ_OSQ=0.0D0
	    SUM_MO=0.0D0;     SUM_LMO=0.0D0
	    DO I=IST,IEND
	      NPIX=NPIX+NINT(MASK(I)*1.0D0)
	      LAM=CD(IP_OBS)%XVEC(I)-CD(IP_OBS)%XVEC(IST)
	      O_DATA=CD(IP_OBS)%DATA(I)*MASK(I)*C_MASK(I)
!
	      SUM_OSQ=SUM_OSQ+O_DATA*O_DATA
	      SUM_LOSQ=SUM_LOSQ+LAM*O_DATA*O_DATA
	      SUM_LSQ_OSQ=SUM_LSQ_OSQ+LAM*LAM*O_DATA*O_DATA
!
	      SUM_MO=SUM_MO+O_DATA*MOD_DATA(I)
	      SUM_LMO=SUM_LMO+LAM*O_DATA*MOD_DATA(I)
	    END DO
	    DET=SUM_OSQ*SUM_LSQ_OSQ-SUM_LOSQ*SUM_LOSQ
	    DA=SUM_MO*SUM_LSQ_OSQ-SUM_LMO*SUM_LOSQ
	    DB=SUM_OSQ*SUM_LMO-SUM_LOSQ*SUM_MO
	    A=DA/DET; B=DB/DET
!
! Compute CHI^2. We always fit the original data.
!
	    RAW_CHISQ=0.0D0; CHISQ=0.0D0
	    DO I=IST,IEND
	      T1=CD(IP_OBS)%XVEC(I)-CD(IP_OBS)%XVEC(IST)
	      RAW_CHISQ=RAW_CHISQ+MASK(I)*(MOD_DATA(I)-OBS_DATA(I))**2/MOD_DATA(I)
	      OBS_DATA(I)=(A+B*T1)*CD(IP_OBS)%DATA(I)
	      CHISQ=CHISQ+MASK(I)*(MOD_DATA(I)-OBS_DATA(I))**2/MOD_DATA(I)
	      CHISQ_CONTRIB(I)=MASK(I)*(MOD_DATA(I)-OBS_DATA(I))**2/MOD_DATA(I)
	    END DO
!
! Simple trapzoidal rule integration.
! Compute model continuum level assuming it is defined close to the line
! bounds of the model.
!
	    T1=0.0D0; T2=0.0D0
	    DO I=IST,IST+2
	      T1=T1+MOD_DATA(I)
	    END DO
	    DO I=IEND-2,IEND
	      T2=T2+MOD_DATA(I)
	    END DO
	    T1=T1/3; T2=T2/3
	    SLOPE=(T2-T1)/(CD(IP_OBS)%XVEC(IEND-1)-CD(IP_OBS)%XVEC(IST+1))
	    INTER=T1
!
! Simple trapzoidal rule integration. These EWs include the omitted regions.
!
	    EW_OBS=0.0D0; EW_MOD=0.0D0
	    DO I=IST,IEND
	      CONT=INTER+SLOPE*(CD(IP_OBS)%XVEC(I)-CD(IP_OBS)%XVEC(IST+1))
	      T1=0.5D0*(CD(IP_OBS)%XVEC(MAX(IST,I-1))-CD(IP_OBS)%XVEC(MIN(I+1,IEND)))
	      EW_OBS=EW_OBS+(1.0D0-OBS_DATA(I)/CONT)*ABS(T1)
	      EW_MOD=EW_MOD+(1.0D0-MOD_DATA(I)/CONT)*ABS(T1)
	    END DO
!
! Simple trapzoidal rule integration. These EWs exclude the omitted regions.
!
	    EW_OBS_OMIT=0.0D0; EW_MOD_OMIT=0.0D0
	    DO I=IST,IEND
	      CONT=INTER+SLOPE*(CD(IP_OBS)%XVEC(I)-CD(IP_OBS)%XVEC(IST+1))
	      IF(MASK(I) .GT. 0.01)THEN
	        J=MAX(I-1,IST)
	        IF(MASK(J) .LT. 0.01)J=I
	        K=MIN(I+1,IEND)
	        IF(MASK(K) .LT. 0.01)K=I
	        T1=0.5D0*(CD(IP_OBS)%XVEC(J)-CD(IP_OBS)%XVEC(K))
	        EW_OBS_OMIT=EW_OBS_OMIT+(1.0D0-OBS_DATA(I)/CONT)*ABS(T1)
	        EW_MOD_OMIT=EW_MOD_OMIT+(1.0D0-MOD_DATA(I)/CONT)*ABS(T1)
	      END IF
	    END DO
!
! We ignore the clipped regions for computing the mean wavelength.
!
	    MEAN=0.0D0; LAM_CENT=0.0D0	
	    DO I=IST,IEND-1
	      X1=CD(IP_OBS)%XVEC(I); X2=CD(IP_OBS)%XVEC(I+1)
	      D1=MOD_DATA(I)-1.0D0
	      D2=MOD_DATA(I+1)-1.0D0
	      MEAN=MEAN+(D1+D2)*(X2-X1)
	      LAM_CENT=LAM_CENT+(X1*D1+X2*D2)*(X2-X1)
	   END DO
	   LAM_CENT=LAM_CENT/MEAN
!
	    WRITE(6,'(A)')' '
	    WRITE(6,'(A)')' Reduced Chi^2 value assumes a SN to 100'
	    WRITE(6,'(A)')' Reduced Chi^2 at another SN =  CHI^2 * (SN/100)^2'
	    WRITE(6,'(A)')' '
	    WRITE(6,'(1X,2(A,F12.3,3X),A,F12.4)')'Raw Chi^2=',RAW_CHISQ,'Red Raw Chi^2=',
	1                  1.0D+04*RAW_CHISQ/(NPIX-2),'Lam(cent)=',LAM_CENT
	    RED_CHISQ=CHISQ*1.0D+04/(NPIX-2)
	    WRITE(6,'(1X,2(A,F12.3,3X))')'    Chi^2=',CHISQ,    '    Red Chi^2=',RED_CHISQ
	    WRITE(6,'(A)')' '
	    FIT_DONE=.TRUE.
	    I=1; CALL PGSCI(I); J=IEND-IST+1
            CALL PGLINE(J,CD(IP_OBS)%XVEC(IST),OBS_DATA(IST))
!
	    T2=EW_MOD; T3=LAM_CENT; T4=50.0   !Rough FWHM in km/s to get closest line
	    CALL GET_LINE_ID_PG(TRANS_NAME,T1,T2,T3,T4)
	    WRITE(LU_OUT,'(I4,F12.4,6F12.3,6X,A)')IP_MOD,LAM_CENT,EW_MOD,EW_OBS,
	1                          EW_MOD_OMIT,EW_OBS_OMIT,RED_CHISQ,T1,TRIM(TRANS_NAME)
	    WRITE(LU_OUT,'(16X,2F12.4,2I12)')XST,XEND,N_OMIT,NC_OMIT
	    IF(N_OMIT .NE. 0)WRITE(LU_OUT,'(16X,8F12.3)')(OLAM_ST(J),OLAM_END(J),J=1,N_OMIT)
	    IF(NC_OMIT .NE. 0)WRITE(LU_OUT,'(16X,8F12.3)')(C_OLAM_ST(J),C_OLAM_END(J),J=1,NC_OMIT)
	    WRITE(LU_OUT,*)' '; FLUSH(LU_OUT)
!
	  ELSE IF(UC(CURSVAL) .EQ. 'H')THEN
	    CALL PRINT_CURSOR_DESC()
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
	      XPAR(1)=XPAR(1)+0.5*T1
	      XPAR(2)=XPAR(1)+T1
	      XLOC=XPAR(1)+XINC
	    ELSE IF(UC(CURSVAL) .EQ. 'P')THEN
	      XPAR(1)=XPAR(1)-0.5*T1
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
	    IF(FIT_DONE)THEN
	      I=1; CALL PGSCI(I); J=IEND-IST+1
              CALL PGLINE(J,CD(IP_OBS)%XVEC(IST),OBS_DATA(IST))
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
	IF(FIT_DONE)THEN
	  IP_OUT=NPLTS+1
	  CALL GEN_IN(IP_OUT,'Outut plot (0 not to save; can be same as IP_OBS ')
	  IF(IP_OUT .EQ. 0)THEN
	  ELSE IF(IP_OUT .EQ. IP_OBS)THEN
	    CD(IP_OBS)%DATA=OBS_DATA
	  ELSE
	    IF(NPTS(IP_OUT) .NE. 0)THEN
	      DEALLOCATE(CD(IP_OUT)%XVEC,CD(IP_OUT)%DATA)
	      IF(ALLOCATED(CD(IP_OUT)%EMIN))THEN
	        DEALLOCATE(CD(IP_OUT)%EMIN,CD(IP_OUT)%EMAX)
	        ERR(IP_OUT)=.FALSE.
	      END IF
	    ELSE
	      CD(IP_OUT)%CURVE_ID='Scaled'
	    END IF
	    ALLOCATE(CD(IP_OUT)%XVEC(NPTS(IP_OBS)),CD(IP_OUT)%DATA(NPTS(IP_OBS)))
	    NPTS(IP_OUT)=NPTS(IP_OBS)
	    CD(IP_OUT)%XVEC=CD(IP_OBS)%XVEC
	    CD(IP_OUT)%DATA=OBS_DATA
	    IF(IP_OUT .GT. NPLTS)NPLTS=NPLTS+1
	  END IF
!
	  IP_OUT=NPLTS+1
	  ALLOCATE(CD(IP_OUT)%XVEC(NPTS(IP_OBS)),CD(IP_OUT)%DATA(NPTS(IP_OBS)))
	  NPTS(IP_OUT)=NPTS(IP_OBS)
	  CD(IP_OUT)%XVEC=CD(IP_OBS)%XVEC
	  CD(IP_OUT)%DATA=CHISQ_CONTRIB
	  CD(IP_OUT)%CURVE_ID='Chi^2 cont.'
	  NPLTS=NPLTS+1
	  WRITE(6,*)'Contribution to chi^2 as a function of pixel in plot',IP_OUT
	END IF
!
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
	  WRITE(6,'(A)')' Use the cursor and B to define left and right side of the line.'
	  WRITE(6,'(A)')'   Then use O to omit regions not to be included in the fit.'
	  WRITE(6,'(A)')' '
	  WRITE(6,'(A)')RED_PEN//' NB: All cursor input should be lower case'//BLUE_PEN
	  WRITE(6,'(A)')' '
	  WRITE(6,'(A)')'   B - define full line region'
	  WRITE(6,'(A)')'   O - define regions to omit from the chi^2 & EW computation.'
	  WRITE(6,'(A)')'   Z - define regions to omit from the chi^2 computaton only.'
	  WRITE(6,'(A)')'   F - do the fitting'
	  WRITE(6,'(A)')'   S - reset and start line selection again'
	  WRITE(6,'(A)')'   Q - quit line selection'
	  WRITE(6,'(A)')'   H - print this message'
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
!
	SUBROUTINE WRITE_BALMER_HEADER(LU)
	USE SET_KIND_MODULE
	IMPLICIT NONE
	INTEGER LU,IP
!
 	WRITE(LU,'(A)')'! '
 	WRITE(LU,'(A)')'! Reduced Chi^2 value assumes a SN to 100.'
 	WRITE(LU,'(A)')'! Reduced Chi^2 at another SN =  Chi^2 * (SN/100)^2.'
 	WRITE(LU,'(A)')'! Lam refers to the centroid wavelength (no regions omitted).'
 	WRITE(LU,'(A)')'! The first  2 EW''s include the omitted regions.'
 	WRITE(LU,'(A)')'! The second 2 EW''s exclude the omitted regions.'
 	WRITE(LU,'(A)')'! '
 	WRITE(LU,'(A)')'! Omit window pairs follow each EW & Chi^2 measurement.'
 	WRITE(LU,'(A)')'! '
	DO IP=1,NPLTS
          WRITE(LU,'(A,I3,5X,A,2X,A)')'! Plot #:',IP,'Plot title:',TRIM(CD(IP)%CURVE_ID)
        END DO
        WRITE(LU,'(A)')'!'
 	WRITE(LU,'(A,1X,A,9X,A,7(5X,A))')'!','IP','Lam','EW(mod)','EW(obs)',
	1                 'EW(mod)','EW(obs)','Chi^2','Lam(ID)','Transition'
	WRITE(LU,'(A,8X,A,5X,A,4X,A)')'!','Line st','Line end','N(omit)'
 	WRITE(LU,'(A)')'! '
	FLUSH(LU)
!	
	END SUBROUTINE WRITE_BALMER_HEADER
	END SUBROUTINE DO_CURSOR_BALMER
