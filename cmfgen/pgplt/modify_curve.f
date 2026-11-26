!
! Simple routine to modify the shape of a curve using cursor input.
! Points can be replaced, deleted, and added. I also possible to
! move a sectioon up/down or left/right.
!
	SUBROUTINE MODIFY_CURVE(IP_IN,IP_OUT)
	USE SET_KIND_MODULE
	USE MOD_CURVE_DATA
        USE MOD_COLOR_PEN_DEF
	IMPLICIT NONE
!
! Altered 08-Jul-2022 : Renamed to MODIFY_CURVE.
! Altered 02-Jun-2022 : Bug fix EMIN was being deallocated twice, EMAX not deallocted.
!
	INTEGER IP_IN
	INTEGER IP_OUT
!
	REAL*4 XV(2*NPTS(IP_IN))
	REAL*4 YV(2*NPTS(IP_IN))
	REAL*4 TMP_XV(2*NPTS(IP_IN))
	REAL*4 WRK_YV(2*NPTS(IP_IN))
!
	INTEGER, PARAMETER :: IONE=1
!
	REAL*4 XVAL,YVAL
	REAL*4 XVAL2,YVAL2
	REAL*4 SYMB_EXP_FAC
	REAL*4 T1
	INTEGER ND
	INTEGER I,J
	INTEGER PGCURS
	INTEGER CURSERR
	CHARACTER(LEN=1) CURSVAL
	LOGICAL FLIP
!
	ND=NPTS(IP_IN)
	IF(CD(IP_IN)%XVEC(1) .GT. CD(IP_IN)%XVEC(ND))THEN
	  DO I=1,ND
	   XV(I)=CD(IP_IN)%XVEC(ND-I+1)
	   YV(I)=CD(IP_IN)%DATA(ND-I+1)
	  END DO
	  FLIP=.TRUE.
	ELSE
	  XV(1:ND)=CD(IP_IN)%XVEC
	  YV(1:ND)=CD(IP_IN)%DATA
	  FLIP=.FALSE.
	END IF
!
	CALL PRINT_CURSOR_DESC
	XVAL=5.0; CALL PGSCH(XVAL)
	DO WHILE(1 .EQ. 1)				!Multiple cursor entries
!
	  CURSERR = PGCURS(XVAL,YVAL,CURSVAL)
	  WRITE(6,*)'Cursor values are:',XVAL,YVAL,CURSVAL
	  IF(CURSVAL .EQ. 'e')EXIT
!
	  IF(CURSVAL .EQ. 'r')THEN
	    T1=XVAL
	    TMP_XV(1:ND)=ABS(XV(1:ND)-XVAL)
	    I=MINLOC(TMP_XV(1:ND),IONE)
	    YV(I)=YVAL
	    J=1; CALL PGSCI(J)
	    CALL PGPT(IONE,XVAL,YVAL,IONE)
	    WRITE(6,*)'Replaced YV at  XV=',XV(I)
!
	  ELSE IF(CURSVAL .EQ. 'd')THEN
	    T1=XVAL
	    TMP_XV(1:ND)=ABS(XV(1:ND)-T1)
	    I=MINLOC(TMP_XV(1:ND),IONE)
	    XV(I:ND-1)=XV(I+1:ND)
	    YV(I:ND-1)=YV(I+1:ND)
	    ND=ND-1
	    WRITE(6,*)'Deleted XV=',XV(I),' from grid'
!
	  ELSE IF(CURSVAL .EQ. 'a')THEN
	    DO I=1,ND-1
	      IF( (XV(I)-XVAL)*(XV(I+1)-XVAL) .LT. 0 )THEN
	        DO J=ND,I+1,-1
                  XV(J+1)=XV(J)
                  YV(J+1)=YV(J)
	        END DO
	        XV(I+1)=XVAL; YV(I+1)=YVAL
	        ND=ND+1
	        J=1; CALL PGSCI(J)
	        CALL PGPT(IONE,XVAL,YVAL,IONE)
	        WRITE(6,*)'Added XV=',XV(I),' to grid'
	        EXIT
	      END IF
	    END DO
!
	  ELSE IF(CURSVAL .EQ. 'm')THEN
	    T1=XVAL
	    TMP_XV(1:ND)=ABS(XV(1:ND)-T1)
	    I=MINLOC(TMP_XV(1:ND),IONE)
	    CURSERR = PGCURS(XVAL2,YVAL2,CURSVAL)
!
	    IF(XVAL2 .LT. XVAL)THEN
	      T1=XVAL2-XVAL
	      DO J=1,I
	        XV(J)=XV(J)+T1
	        XVAL=XV(J); YVAL=YV(J)
	        CALL PGPT(IONE,XVAL,YVAL,IONE)
	      END DO
	    ELSE
	      T1=XVAL2-XVAL
	      DO J=I,ND
	        XV(J)=XV(J)+T1
	        XVAL=XV(J); YVAL=YV(J)
	        CALL PGPT(IONE,XVAL,YVAL,IONE)
	      END DO
	    END IF
	    WRITE(6,*)'Shited XV grid by',T1
!
	  ELSE IF(CURSVAL .EQ. 'v')THEN
	     TMP_XV(1:ND)=ABS(XV(1:ND)-XVAL)
	     I=MINLOC(TMP_XV(1:ND),IONE)
	     J=1; CALL PGSCI(J)
	     T1=YVAL-YV(I)
	     CURSERR = PGCURS(XVAL2,YVAL2,CURSVAL)
	     IF(CURSVAL .EQ. 'l')THEN
	       DO J=1,I
	         YV(J)=YV(J)+T1
	         XVAL=XV(J); YVAL=YV(J)
	         CALL PGPT(IONE,XVAL,YVAL,IONE)
	       END DO
	     ELSE IF(CURSVAL .EQ. 'r')THEN
	       DO J=I,ND
	         YV(J)=YV(J)+T1
	         XVAL=XV(J); YVAL=YV(J)
	         CALL PGPT(IONE,XVAL,YVAL,IONE)
	       END DO
	     ELSE IF(CURSVAL .EQ. 'a')THEN
	       DO J=I,ND
	         YV(J)=YV(J)+T1
	         XVAL=XV(J); YVAL=YV(J)
	         CALL PGPT(IONE,XVAL,YVAL,IONE)
	       END DO
	     END IF 	
	     WRITE(6,*)'Shifted Y grid by',T1
!
	   ELSE
	      WRITE(6,*)RED_PEN
	      WRITE(6,*)'Error - bad cursor input'
	      CALL PRINT_CURSOR_DESC
	   END IF
	 END DO
!
	 IF(IP_OUT .GT. NPLTS+1)IP_OUT=NPLTS+1
	 IF(IP_OUT .GT. NPLTS)NPLTS=NPLTS+1
	 IF(ALLOCATED(CD(IP_OUT)%XVEC))THEN
	   DEALLOCATE(CD(IP_OUT)%XVEC,CD(IP_OUT)%DATA)
	 END IF
	 IF(ALLOCATED(CD(IP_OUT)%EMIN))THEN
	   DEALLOCATE(CD(IP_OUT)%EMIN,CD(IP_OUT)%EMAX)
	 END IF
	 NPTS(IP_OUT)=ND
	 ALLOCATE(CD(IP_OUT)%XVEC(ND))
	 ALLOCATE(CD(IP_OUT)%DATA(ND))
!
	 IF(FLIP)THEN
	   DO I=1,ND
	    CD(IP_OUT)%XVEC(I)=XV(ND-I+1)
	    CD(IP_OUT)%DATA(I)=YV(ND-I+1)
	   END DO
	 ELSE
	   CD(IP_OUT)%XVEC(1:ND)=XV(1:ND)
	   CD(IP_OUT)%DATA(1:ND)=YV(1:ND)
	 END IF
	 ERR(IP_OUT)=.FALSE.
	 CD(IP_OUT)%CURVE_ID=' '
!
	RETURN
!
	CONTAINS
	SUBROUTINE PRINT_CURSOR_DESC
	USE SET_KIND_MODULE
          WRITE(6,'(A)')BLUE_PEN
          WRITE(6,*)'Use ''r'' to replace a data point'
          WRITE(6,*)'Use ''a'' to add a data point'
          WRITE(6,*)'Use ''d'' to delete a data point'
          WRITE(6,*)'Use ''m'' move next, and points further out/inward by a fixed amount'
	  WRITE(6,*)'          Requires second cursor input to specify size and direction of shift'
          WRITE(6,*)'Use ''v'' move next, and points further  up/down by a fixed amount in V'
          WRITE(6,*)'           Requires second cursor input to specify add'
          WRITE(6,*)'             cursor = l(left), r(ight), a(all)'
          WRITE(6,'(A)')DEF_PEN
	  RETURN
	END SUBROUTINE PRINT_CURSOR_DESC
	END SUBROUTINE MODIFY_CURVE
