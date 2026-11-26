!
! Subroutine to fits a straigt line, and a set of modified Gaussians,
! to a set of data. Routine is designed to be called in GRAMON_PGPLOT.
!
	SUBROUTINE DO_GAUS_FIT(XST_PASSED,XEND_PASSED)
	USE SET_KIND_MODULE
	USE MOD_CURVE_DATA
	USE GAUS_FIT_DATA
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
! Altered 09-Aug-2022: Now allows automatic movement along a spectrum.
!                      Selection of initial Gauss parameters now done by cursor.
!                      This routine is being superceded by GAUSS_FIT.
! Altered 27-Feb-2022: Fixed bug in EW error estimates when using unormalized data.
! Altered 25-Feb-2022: Some cleaning
! Altered 01-Feb-2022: Extensive changes implmented to allow cursor input of
!                        parameters for gauss fitting.
! Altered   -Sep-07
!
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: NL_MAX=100
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
	REAL(KIND=LDP) LINE_HEIGHT(NL_MAX)
	REAL(KIND=LDP) LINE_CENTER(NL_MAX)
	REAL(KIND=LDP) LINE_SIGMA(NL_MAX)
!
	REAL*4 XST_PASSED,XEND_PASSED
	REAL*4, SAVE :: FIND_LINE_TOLERANCE=0.002
!
	REAL(KIND=LDP) TOLERANCE
	REAL(KIND=LDP) GAUS_FIT_FUNC
	EXTERNAL GAUS_FIT_FUNC
!
	REAL(KIND=LDP) YST,YEND
!
	REAL(KIND=LDP) TOL
	REAL(KIND=LDP) T1
	REAL(KIND=LDP) FWHM
	INTEGER NO_LINES
	INTEGER ITER
!
	LOGICAL USE_CURSOR
	REAL*4 XVAL,YVAL
	REAL*4 DEF_SIGMA
        INTEGER PGCURS
        INTEGER CURSERR
	CHARACTER(LEN=1) CURSVAL
!
	REAL(KIND=LDP), SAVE :: XST=0.0D0
	REAL(KIND=LDP), SAVE :: XEND=1.0D0
	REAL*4 SP_XST,SP_XEND
	REAL*4, SAVE :: XST_PASSED_SAVED,XEND_PASSED_SAVED
	INTEGER, SAVE :: IP=1
	INTEGER, SAVE :: NG_PAR_OLD
        LOGICAL, PARAMETER :: L_TRUE=.TRUE.
	LOGICAL, SAVE :: FIRST_WRITE=.TRUE.
	LOGICAL NEW_FILE
	LOGICAL WRITE_FIT
	LOGICAL FILE_EXISTS
	LOGICAL DO_FIND
!
	INTEGER I,J,K
	LOGICAL GUESSED
	LOGICAL NEW_REGION
	LOGICAL SHOW_INIT_GAUS
!
	NG_PAR_MAX=4*NL_MAX+2
	IF(.NOT. ALLOCATED(PAR))THEN
	  ALLOCATE(PAR(NG_PAR_MAX))
	END IF
!
! If we are fitting the same regon as previously, we can use the previous findings.
!
	WRITE(6,*)NPLTS,NPTS(1)
	IF(XST_PASSED_SAVED .EQ. XST_PASSED .AND. XEND_PASSED_SAVED .EQ. XEND_PASSED)THEN
	   NEW_REGION=.FALSE.
	ELSE
	  NUM_GAUS=0
	  XST_PASSED_SAVED=XST_PASSED
	  XEND_PASSED_SAVED=XEND_PASSED
	  NEW_REGION=.TRUE.
	END IF
!
	USE_CURSOR=.FALSE.
	CALL GEN_IN(IP,'Plot for fitting')
	IF(.NOT. NEW_REGION)CALL GEN_IN(NEW_REGION,'Reset Gauss selection')
	IF(NEW_REGION)THEN
	  USE_CURSOR=.TRUE.
	  CALL GEN_IN(USE_CURSOR,'Use cursor')
	END IF
	DEF_SIGMA=0.1
	CALL GEN_IN(DEF_SIGMA,'Default sigma for line profiles')
	SHOW_INIT_GAUS=.FALSE.
	CALL GEN_IN(SHOW_INIT_GAUS,'Show initial Gaussian estimates on the plot?')
!
	USE_CURSOR=.TRUE.
	IF(.NOT. NEW_REGION)THEN
	  NG_PAR=NG_PAR_OLD
!
	ELSE IF(USE_CURSOR)THEN
!
	  WRITE(6,'(A)')' '
	  WRITE(6,'(A)')' First 2 cursor inputs to define input band '
	  WRITE(6,'(A)')' Then use each cursor position to indicate height and location of line'
	  WRITE(6,'(A)')' '
!
	  XVAL=5.0; CALL PGSCH(XVAL)
	  J=1; CALL PGSCI(J)
	  XVAL=0.5D0*(XST_PASSED+XEND_PASSED)
	  YVAL=1.0
	  CURSERR = PGCURS(XVAL,YVAL,CURSVAL)
	  CALL PGPT(IONE,XVAL,YVAL,IONE)
	  XST=XVAL; YST=YVAL
	  CURSERR = PGCURS(XVAL,YVAL,CURSVAL)
	  CALL PGPT(IONE,XVAL,YVAL,IONE)
	  XEND=XVAL; YEND=YVAL
	  PAR(1)=YST
	  PAR(2)=(YEND-YST)/(XEND-XST)
	  WRITE(6,*)PAR(1),PAR(2)
	  DO J=1,NL_MAX
	    CURSERR = PGCURS(XVAL,YVAL,CURSVAL)
	    CALL PGPT(IONE,XVAL,YVAL,IONE)
	    IF(CURSVAL .EQ. 'e' .OR. CURSVAL .EQ. 'E')EXIT
	    K=2+(J-1)*4+1
	    PAR(K)=XVAL; PAR(K+2)=YVAL-(PAR(1)+PAR(2)*(XVAL-XST))
	    PAR(K+3)=2.0
	    PAR(K+1)=DEF_SIGMA
	    NUM_GAUS=J
	  END DO
	  NG_PAR=2+4*NUM_GAUS
!
	ELSE
!
	  XST_PASSED_SAVED=XST_PASSED
	  XEND_PASSED_SAVED=XEND_PASSED
	  CALL GEN_IN(XST,'Start wavelength for fitting')
	  CALL GEN_IN(XEND,'End wavelength for fitting')
	  CALL GEN_IN(NUM_GAUS,'Number of gaussians to fit: (0 to find automatically)')
	  GUESSED=.FALSE.
!
! Find Lines automatically.
!
	  IF(NUM_GAUS .EQ. 0)THEN
	    DO_FIND=.TRUE.
	    DO WHILE(DO_FIND)
	      CALL GEN_IN(FIND_LINE_TOLERANCE,'Departure from unity for lines')
	      FIND_LINE_TOLERANCE=ABS(FIND_LINE_TOLERANCE)
	      SP_XST=XST; SP_XEND=XEND
	      WRITE(6,*)IP,NPTS(IP),SP_XST,SP_XEND,T1
	      CALL FIND_LINES(LINE_CENTER, LINE_HEIGHT, LINE_SIGMA, NUM_GAUS, NL_MAX,
	1                 CD(IP)%DATA(1),CD(IP)%XVEC(1),NPTS(IP),
	1                 SP_XST,SP_XEND,FIND_LINE_TOLERANCE)
	      DO_FIND=.FALSE.
	      WRITE(6,*)'You may hand edit individual Gaussians or do auto find again'
	      CALL GEN_IN(DO_FIND,'Do automatic find again with different tolerance?')
	    END DO
	    GUESSED=.TRUE.
	    NG_PAR=4*NUM_GAUS+2
!
! If possible, we use previous fit parameters as default.
!
	    PAR(1)=1.0D0		!Mean value
	    PAR(2)=0.0D0		!Continuum slope
	    DO J=1,NUM_GAUS
	      K=3+(J-1)*4
	      PAR(K)=LINE_CENTER(J)
	      PAR(K+1)=LINE_SIGMA(J)
	      PAR(K+2)=LINE_HEIGHT(J)		!-ve if absorption line
	      PAR(K+3)=2.0D0			!i.e., assume Gaussian
	    END DO
!
! Add Lines by hand
!
	  ELSE IF(NEW_REGION)THEN
	    CALL GEN_IN(DEF_SIGMA,'Default sigma')
	    PAR=0.0D0; PAR(1)=1.0D0
	    CALL GEN_IN(PAR(1),'Mean value')
	    CALL GEN_IN(PAR(2),'Continuum slope')
	    DO J=1,NUM_GAUS
	      K=2+(J-1)*4+1
	      CALL GEN_IN(PAR(K),'Central wavlength of Gaussian')
              IF(PAR(K+1) .EQ. 0.0D0 .AND. J .GT. 1)PAR(K+1)=PAR(K-3)
	      CALL GEN_IN(DEF_SIGMA,'Modified sigma (scale) of Gaussian')
	      PAR(K+1)=DEF_SIGMA
	      IF(PAR(K+2) .EQ. 0.0D0)PAR(K+2)=-0.2D0
	      CALL GEN_IN(PAR(K+2),'Offset from continuum (-ve for absorption)')
	      IF(PAR(K+3) .EQ. 0.0D0)PAR(K+3)=2.0D0
	      CALL GEN_IN(PAR(K+3),'Guassian exponent')
	      WRITE(6,*)'You may now make corrections to the Gaussian data'
	    END DO
	  END IF
	END IF
!
! Draw initial Gaussians using a dashed black pen.
!
	IF(SHOW_INIT_GAUS)THEN
	  I=2; CALL PGSLS(I)			!Dashed curve
	  CALL SET_GAUS_DATA(CD(IP)%XVEC,CD(IP)%DATA,XST,XEND,NPTS(IP),YST,YEND)
	  CALL DRAW_GAUS(L_FALSE)
	  CALL PGEBUF				!Clear buffer
	  I=1; CALL PGSLS(I)			!Reset to solid lines
	END IF
!
! Allow Gaussians to be adjustd by hand.
!
	NG_PAR=2+4*NUM_GAUS
	CALL ED_GAUS_FIT()
	NG_PAR=2+4*NUM_GAUS
	NG_PAR_OLD=NG_PAR
!
	WRITE(6,*)'Exited ED_GAUS_FIT'; FLUSH(UNIT=6)
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
	  CALL GAUS_ROMB(EW(I),SIM(1,K+2),SIM(1,K+1),SIM(1,K+3),TOLERANCE)
	  T1=SIM(1,1)+SIM(1,2)*(SIM(I,K)-X_GAUS(1))				!Continuum value
	  EW(I)=EW(I)/T1
	  EW_CONT(I)=T1
	END DO
	EW=EW*1000.0D0			!mAng
        PAR(1:NG_PAR)=SIM(1,1:NG_PAR)
	CALL GAUS_FIT_ER(PAR)
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
	1               T1/SQRT(2.0D0),FWHM,EW(I),EW_ERROR(I),ALT_ERROR(I),MIN_ERROR(I)
	END DO
!
        SUM_SQ(:)=0.0D0
        SUM_SQ(1)=GAUS_FIT_FUNC(PAR)
!
! Draw Gaussians using a black pen.
!
	I=1
	CALL PGSCI(I)
	CALL PGLINE(NG_DATA,XFIT,YFIT)
	CALL PGEBUF				!Clear buffer
!
	CALL GEN_IN(WRITE_FIT,'Write fit to file')
	IF(WRITE_FIT)THEN
	  IF(FIRST_WRITE)THEN
	    FIRST_WRITE=.FALSE.
	    INQUIRE(FILE='GAUSS_FITS',EXIST=NEW_FILE); NEW_FILE=.NOT. NEW_FILE
	    IF(.NOT. NEW_FILE)THEN
	      WRITE(6,*)'If file GAUS_FITS exists: Default is to append new fits'
	      CALL GEN_IN(NEW_FILE,'Open new file (T) or append data (F)')
	    END IF
	    IF(NEW_FILE)THEN
	      OPEN(UNIT=35,FILE='GAUSS_FITS',STATUS='UNKNOWN')
	      WRITE(35,'(2X,(5X,A),5(5X,A),2(3X,A))')
	1          '        Lam','     Height','          a',
	1          '        EXP','Sigma(km/s)',' FWHM(km/s)',' EW(mA)','Err(mA)'
	    ELSE
	      OPEN(UNIT=35,FILE='GAUSS_FITS',STATUS='UNKNOWN',POSITION='APPEND')
	    END IF
	  ELSE
	    OPEN(UNIT=35,FILE='GAUSS_FITS',STATUS='UNKNOWN',POSITION='APPEND')
	  END IF
!
! We output lines in wavelength order. After calling, INDEXX,
! INDX_VEC(1) will contain the first line to be output, etc.
! We don't actually perform the sort.
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
	  WRITE(35,'(A)')' '
	  WRITE(35,'(I3,8X,A)')NUM_GAUS,'!Number of Gaussians'
	  WRITE(35,'(2X,2ES16.6)')XST,XEND
	  WRITE(35,'(2X,2ES16.6)')PAR(1:2)
	  DO J=1,NUM_GAUS
	     I=INDX_VEC(J)
	     K=2+(I-1)*4+1
	     T1=2.99794D+05*SIM(1,K+1)/SIM(1,K)
	     FWHM=2.0D0*T1*(LOG(2.0D0))**(1.0D0/SIM(I,K+3))
	     WRITE(35,'(2X,ES16.6,5ES16.4,4F10.2)')SIM(I,K),SIM(1,K+2),SIM(1,K+1),SIM(I,K+3),
	1               T1/SQRT(2.0D0),FWHM,EW(I),EW_ERROR(I),ALT_ERROR(I),MIN_ERROR(I)
	  END DO
	  CLOSE(UNIT=35)
	 END IF
	WRITE(6,*)'Final SUM_SQ=',SUM_SQ(1)
	WRITE(6,*)'Use P, then DG to redraw Gauss fit'
	WRITE(6,*)'Use GF to edit and redo the Gauss fit'
!	
	RETURN
	END
!
! Routine to draw the fit on top of an existing plot. At present
! plot is drawn in black only. Routine can also draw the difference
! betwene the data and the fitted curve.
!
	SUBROUTINE DRAW_GAUS(DIFF_ALSO)
	USE SET_KIND_MODULE
	USE MOD_CURVE_DATA
	USE GAUS_FIT_DATA
	USE GEN_IN_INTERFACE
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
!
	RETURN
	END
