	MODULE MOD_EW_VARIABLES
	REAL*4 XST,XEND                 !Line limits in world coordinates
	REAL*4 XMEAN
	REAL*4 YCONT
        REAL*4 EW  			!
	REAL*4 EWL,EWH			!Line equivalent width with a % shift in continuum
	REAL*4 SIGMA			!Standard deviation
	REAL*4 SKEWNESS
	REAL*4 KURTOSIS
	REAL*4 FWHM
	REAL*4 LINE_WAVE
	REAL*4 XLOW_FWHM  		!X value at 50% from line center
	REAL*4 XHIGH_FWHM		!X value at 50% from line center.
!
! Work array.
!
	REAL*4, ALLOCATABLE :: ONE_MIN_FDFC(:)
!
	LOGICAL :: USE_MILLI_ANG=.TRUE.
	CHARACTER(LEN=80) TRANS_NAME
	SAVE
	END MODULE MOD_EW_VARIABLES
!
	SUBROUTINE WRITE_EW_HEADER(LUOUT)
	USE SET_KIND_MODULE
	USE MOD_CURVE_DATA
	USE MOD_EW_VARIABLES
	IMPLICIT NONE
	INTEGER LUOUT
	INTEGER IP
!
	WRITE(LUOUT,'(A)')'!'
	WRITE(LUOUT,'(A)')'! Estimates of the error in the EW values may be poor.'
	WRITE(LUOUT,'(A)')'!'
	DO IP=1,NPLTS
	  WRITE(LUOUT,'(A,I3,5X,A,2X,A)')'! Plot #:',IP,'Plot title:',TRIM(CD(IP)%CURVE_ID)
	END DO
	WRITE(LUOUT,'(A)')'!'
!
	IF(USE_MILLI_ANG)THEN
	  WRITE(LUOUT,'(A,5(3X,A),6X,A,2X,3A,20X,3A)')
	1           '!','    XST',' XEND','   Line Loc','    F(cont)',
	1               '  EW(mA)','%E','FWHM(km/s)',' FWHM(km/s)',' Sig(A)',
	1               '  Plot','    Line Lam','   Line ID'
	ELSE
	  WRITE(LUOUT,'(A,5(3X,A),X,A,2X,3A,20X,3A)')
	1           '!','    XST',' XEND','   Line Loc','    F(cont)',
	1               '   EW(A)','%E','FWHM(km/s)',' FWHM(km/s)',' Sig(A)',
	1               '  Plot','    Line Lam','   Line ID'
	END IF
!
	RETURN
	END

	SUBROUTINE WR_EW_VARIABLES(IP,LUOUT)
	USE SET_KIND_MODULE
	USE MOD_EW_VARIABLES
	IMPLICIT NONE
	INTEGER IP,LUOUT
	REAL*4 T1,T2
!
! Since the percentage errors computed from EWL and EWH are the same
! (apart from the sign) we simply print the average value.
!
	T1=2.998D+05*SIGMA/XMEAN
	T2=(50.0D0*(EW-EWL)+50.0D0*(EW-EWH))/EW
	IF(USE_MILLI_ANG)THEN
	  WRITE(LUOUT,'(3F11.3,ES14.4,F11.2,F8.2,2F10.2,3F10.3,3X,I3,1X,F11.3,3X,A)')
	1                XST,XEND,XMEAN,YCONT,1000.0*EW,T2,
	1                FWHM,2.333*T1,
	1                SIGMA,SKEWNESS,KURTOSIS,IP,
	1                LINE_WAVE,TRIM(TRANS_NAME)
	ELSE
	  WRITE(LUOUT,'(3F14.3,2ES14.4,F8.2,2F10.3,3F10.3,3X,I3,1X,F11.3,3X,A)')
	1                XST,XEND,XMEAN,YCONT,EW,T2,
	1                FWHM,2.333*T1,
	1                SIGMA,SKEWNESS,KURTOSIS,IP,
	1                LINE_WAVE,TRIM(TRANS_NAME)
	END IF
	FLUSH(LUOUT)
!
	RETURN
	END
