!
! Program to modify RVSIG_COL. Various options are available.
! Ideal for revising grid etc.
!
	SUBROUTINE CHANGE_XAXIS_GRIDDING(R,ND,ND_MAX)
	USE SET_KIND_MODULE
	USE MOD_CURVE_DATA
	USE GEN_IN_INTERFACE
	USE MOD_COLOR_PEN_DEF
	IMPLICIT NONE
!
	INTEGER ND
	INTEGER ND_MAX
	REAL(KIND=LDP) R(ND_MAX)
!
	INTEGER ND_OLD
	REAL(KIND=LDP), ALLOCATABLE :: OLD_R(:)
	REAL(KIND=LDP), ALLOCATABLE :: OLD_DATA(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: TMP_R(:)
	REAL(KIND=LDP), ALLOCATABLE :: TA(:)
!
	REAL(KIND=LDP) T1,T2
	REAL*4 XVAL,YVAL
!
	INTEGER, PARAMETER :: IONE=1
	INTEGER SGN
	INTEGER ID
	INTEGER NCUR
	INTEGER I,J,L,K
	INTEGER NEW_ND
	INTEGER PGCURS
	INTEGER CURSERR
!
	LOGICAL REPLOT
	CHARACTER(LEN=1) CURSVAL
	CHARACTER(LEN=10) OPTION
!
	ND=NPTS(1);   ND_OLD=ND
	ALLOCATE(OLD_R(ND_OLD))
	ALLOCATE(OLD_DATA(ND_OLD,NPLTS))
!
! To allow for the addtion of extra points, we allow a larger grid size.
!
	ALLOCATE(TMP_R(ND_MAX))
	ALLOCATE(TA(ND_MAX))
!
	OLD_R(1:ND)=CD(1)%XVEC(1:ND); R(1:ND)=OLD_R(1:ND)
	DO ID=1,NPLTS
	  OLD_DATA(1:ND_OLD,ID)=CD(ID)%DATA(1:ND_OLD)
	END DO
!
	XVAL=5.0; CALL PGSCH(XVAL)
	XVAL=CD(1)%XVEC(ND/2)
	YVAL=CD(1)%DATA(ND/2)
	DO WHILE(1 .EQ. 1)				!Multiple plotting
!
	  WRITE(6,'(A)')BLUE_PEN
	  WRITE(6,'(A)')' '
	  WRITE(6,*)'Click on top bar of plot window to activate cursor'
	  WRITE(6,'(A)')' '
	  WRITE(6,*)'Use ''a'' to add a data point'
	  WRITE(6,*)'Use ''d'' to delete a data point'
	  WRITE(6,*)'Use ''b'' replace band with n new data points'
	  WRITE(6,*)'Use ''e'' to exit'
	  WRITE(6,'(A)')DEF_PEN
!
	  DO WHILE(1 .EQ. 1)				!Multiple cursor entries
!
	    CURSERR = PGCURS(XVAL,YVAL,CURSVAL)
	    WRITE(6,*)'Cursor values are:',XVAL,YVAL
	    IF(CURSVAL .EQ. 'e' .OR. CURSVAL .EQ. 'E')EXIT
	    IF(CURSVAL .EQ. 'd' .OR. CURSVAL .EQ. 'D')THEN
	      T1=XVAL
	      TMP_R(1:ND)=ABS(R(1:ND)-T1)
	      I=MINLOC(TMP_R(1:ND),IONE)
	      R(I:ND-1)=R(I+1:ND)
	      ND=ND-1
	      WRITE(6,*)'Deleted R=',R(I),' from grid (ND,ND_OLD)',ND, ND_OLD
	    ELSE IF(CURSVAL .EQ. 'a' .OR. CURSVAL .EQ. 'A')THEN
	      DO I=1,ND-1
	        IF( (R(I)-XVAL)*(R(I+1)-XVAL) .LT. 0 )THEN
	          DO J=ND,I+1,-1
                    R(J+1)=R(J)
	          END DO
	          IF(ND+1 .GT. ND_MAX)THEN
	            WRITE(6,*)'Too many points added: ND_MAX is',ND_MAX
	            EXIT
	          END IF
	          R(I+1)=XVAL; ND=ND+1
	          J=2; CALL PGSCI(J)
	          CALL PGPT(IONE,XVAL,YVAL,IONE)
	          WRITE(6,*)'Add R=',R(I),' to grid (ND,ND_OLD)',ND, ND_OLD
	          EXIT
	        END IF
	      END DO
!
	    ELSE IF(CURSVAL .EQ. 'b' .OR. CURSVAL .EQ. 'B')THEN
	      T1=XVAL; TMP_R(1:ND)=ABS(R(1:ND)-T1)
	      I=MINLOC(TMP_R(1:ND),IONE)
	      CURSERR = PGCURS(XVAL,YVAL,CURSVAL)
	      WRITE(6,*)'Cursor values are:',XVAL,YVAL
	      T1=XVAL; TMP_R(1:ND)=ABS(R(1:ND)-T1)
	      J=MINLOC(TMP_R(1:ND),IONE)
	      K=MAX(I,J); I=MIN(I,J); J=K
	      WRITE(6,*)' Number of data points inside the grid is ',J-I-1
!
	      NCUR=0
	      DO WHILE(NCUR .EQ. 0)
	        CURSERR = PGCURS(XVAL,YVAL,CURSVAL)
	        WRITE(6,'(A)')'Cursor valie is',CURSVAL
	        IF(CURSVAL .EQ. '+')THEN
	          CURSERR = PGCURS(XVAL,YVAL,CURSVAL)
	          L=ICHAR(CURSVAL)-ICHAR('0')
	          CURSERR = PGCURS(XVAL,YVAL,CURSVAL)
	          NCUR=ICHAR(CURSVAL)-ICHAR('0')
	          NCUR=10*L+NCUR
	          WRITE(6,*)'NCUR=',NCUR
	        ELSE IF(CURSVAL .GT. '0' .AND. CURSVAL .LE. '9')THEN
	          NCUR=ICHAR(CURSVAL)-ICHAR('0')
	          WRITE(6,*)'NCUR=',NCUR
	        END IF
	      END DO
	      NEW_ND=ND+NCUR-(J-I-1)
	      IF(NEW_ND .GT. ND_MAX)THEN
	         WRITE(6,*)'Too many points added: ND_MAX is',ND_MAX
	         EXIT
	      END IF
!
	      TMP_R(I+NCUR+1:NEW_ND)=R(J:ND)
	      T1=(R(J)-R(I))/(NCUR+1)
	      DO L=1,NCUR
	        TMP_R(I+L)=R(I)+T1*L
	      END DO
	      R(I+1:NEW_ND)=TMP_R(I+1:NEW_ND)
	      ND=NEW_ND
!
	    END IF
	  END DO
!
! Check if grid is monotonic.
!
	  SGN=1; IF(OLD_R(ND_OLD) .GT. OLD_R(1))SGN=-1
	  DO I=1,ND-1
	    IF( SGN*(R(I)-R(I+1)) .LE. 0)THEN
	      WRITE(6,*)'Error non-montonic grid'
	      DO J=MAX(1,I-2),MIN(ND-1,I+3)
	        WRITE(6,*)J,R(J),R(J)-R(J+1)
	      END DO
	      EXIT
	    END IF
	  END DO
!
	  REPLOT=.TRUE.
	  CALL GEN_IN(REPLOT,'Replot to see grid revision (otherwise exit)')
	  IF(REPLOT)THEN
	    DO ID=1,NPLTS
	      CALL MON_INTERP(TA,ND,IONE,R,ND,OLD_DATA(1,ID),ND_OLD,OLD_R,ND_OLD)
	      IF(ND .NE. NPTS(ID))THEN
	        DEALLOCATE(CD(ID)%XVEC,CD(ID)%DATA)
	        ALLOCATE(CD(ID)%XVEC(ND),CD(ID)%DATA(ND))
	        NPTS(ID)=ND
	      END IF
	      CD(ID)%XVEC=R(1:ND)
	      CD(ID)%DATA=TA(1:ND)
	    END DO
!
! Option NOI ensures that the data remains intact after exiting GRAMON.
!
 	    CALL GRAMON_PGPLOT(' ',' ',' ','NOI')
	  ELSE
	    EXIT
	  END IF
	END DO
	IF(ALLOCATED(TA))DEALLOCATE(OLD_DATA,OLD_R,TA,TMP_R)
!
	RETURN
	END
