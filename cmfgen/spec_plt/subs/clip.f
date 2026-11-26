!
! Simple subroutine to replace features in the continuum by a straiht
! line fite. The data is read in from DATA_CLIP and should contain the
! regions in vacuum wavelength.
!
!       WAVE1, WAVE2, WAVE3, WAVE4
!
! The clipping region is WAVE2 to WAVE3.
! WAVE1 to WAVE2 and WAVE3 to WAVE4 are use to define the best fine used
! to replace the cliiped region.
!
	SUBROUTINE CLIP(XV,YV,NCF)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Created 18-Apr-2020
!
	INTEGER NCF
	REAL(KIND=LDP) XV(NCF)
	REAL(KIND=LDP) YV(NCF)
!
	INTEGER LU
	INTEGER IOS
!
	REAL(KIND=LDP) SPEED_OF_LIGHT
	EXTERNAL SPEED_OF_LIGHT
!
	INTEGER, PARAMETER :: NCLIP_MAX=20
	REAL(KIND=LDP) FREQ(4,NCLIP_MAX)
!
	REAL(KIND=LDP) X0,T1
	REAL(KIND=LDP) X_SUM,Y_SUM
	REAL(KIND=LDP) XSQ_SUM,XY_SUM
	REAL(KIND=LDP) INTERCEPT,SLOPE
!
	INTEGER NCLIP
	INTEGER NC
	INTEGER I
	INTEGER I1,I2,I3,I4
	INTEGER NDATA
!
	LU=7
	OPEN(UNIT=LU,FILE='CLIP_DATA',STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(6,*)'Unable to open CLIP_DATA filie'
	  WRITE(6,*)'IOSTAT=',IOS
	  RETURN
	END IF
!
	READ(LU,*)NCLIP
	IF(NCLIP .GT. NCLIP_MAX)THEN
	  WRITE(6,*)'Error in clip.f'
	  WRITE(6,*)'Currently only able to read in ',NCLIP,' clipping regions'
	  NCLIP=NCLIP_MAX
	END IF
!
	DO I=1,NCLIP
	  READ(LU,*)FREQ(1:4,I)
	  FREQ(1:4,I)=1.0E-07_LDP*SPEED_OF_LIGHT()/FREQ(1:4,I)
	END DO
	CLOSE(LU)
	WRITE(6,*)'Successfully read CLIP data'
!
! Determine the cliping regions in pixel space.
!
	DO NC=1,NCLIP
	  DO I=1,NCF
	    IF(XV(I) .LT. FREQ(1,NC))THEN
	      I1=I
	      EXIT
	    END IF
	  END DO
	  DO I=I1+1,NCF
	    IF(XV(I) .LT. FREQ(2,NC))THEN
	      I2=I
	      EXIT
	    END IF
	  END DO
	  DO I=I2+1,NCF
	    IF(XV(I) .LT. FREQ(3,NC))THEN
	      I3=I
	      EXIT
	    END IF
	  END DO
	  DO I=I3+1,NCF
	    IF(XV(I) .LT. FREQ(4,NC))THEN
	      I4=I
	      EXIT
	    END IF
	  END DO
!
!	  WRITE(6,'(A,4F9.6)')'  Freq range',FREQ(1:4,NC)
!	  WRITE(6,'(A,4F9.6)')'  Freq range',XV(I1),XV(I2),XV(I3),XV(I4)
!	  WRITE(6,'(A,4I8)')' Pixel range',I1,I2,I3,I4
!
! We use X0 so that the line intercept is defined at X0, rather than 0
!
	  X_SUM=0.0_LDP; XSQ_SUM=0.0_LDP; Y_SUM=0.0_LDP;; XY_SUM=0.0_LDP
	  X0=XV(I1)
	  DO I=I1,I2
	    X_SUM=X_SUM+(XV(I)-X0)
	    Y_SUM=Y_SUM+YV(I)
	  END DO
	  DO I=I3,I4
	    X_SUM=X_SUM+(XV(I)-X0)
	    Y_SUM=Y_SUM+YV(I)
	  END DO
!
	  NDATA=(I2-I1)+(I4-I3)+2
	  DO I=I1,I2
	    T1=XV(I)-X0
	    XY_SUM=XY_SUM+(T1-X_SUM/NDATA)*(YV(I)-Y_SUM/NDATA)
	    XSQ_SUM=XSQ_SUM+(T1-X_SUM/NDATA)**2
	  END DO
	  DO I=I3,I4
	    T1=XV(I)-X0
	    XY_SUM=XY_SUM+(T1-X_SUM/NDATA)*(YV(I)-Y_SUM/NDATA)
	    XSQ_SUM=XSQ_SUM+(T1-X_SUM/NDATA)**2
	  END DO
!
	  SLOPE=XY_SUM/XSQ_SUM
	  INTERCEPT=(Y_SUM-SLOPE*X_SUM)/NDATA
!
! Replace the data
!
	  DO I=I2,I3
	    YV(I)=INTERCEPT+SLOPE*(XV(I)-X0)
	  END DO
	END DO
!
	RETURN
	END
