!
! Routine to define a continuum using file or cursor input.
! Two basic control options:
!      (a) Continuum values, X and Y, are specified
!      (b) Band  X1 and X2, are specified, and Y is determined by
!            average of data in band X1 to X2. This method is useful
!            for automatic continuum defintions.
!
! Data is output to a file, however you will need to update the number
! of data points in this file.
!
! Data can be input via CURSOR or read from a file. Comments at the top
! of the file are okay. However the last two lines before the data  must
! contain:
!                N            !Number of data points
!
!                X1   X2      !Data type
!   or
!                X     Y      !Data type
!
	SUBROUTINE PG_DEF_CONTINUUM(CONT,IPLT,OPLT,INPUT_OPTION,CURVE_OPTION,RET_IOS)
	USE SET_KIND_MODULE
	USE MOD_CURVE_DATA
	IMPLICIT NONE
!
! Altered 3-Oct-2021       Input and output options improved
!
	INTEGER RET_IOS
	INTEGER IPLT
	INTEGER OPLT
	REAL*4 CONT(NPTS(IPLT))
	CHARACTER(LEN=*) INPUT_OPTION
	CHARACTER(LEN=*) CURVE_OPTION
	CHARACTER(LEN=120) UPPER_CASE
	CHARACTER(LEN=120) LOCAL_OPTION
	CHARACTER(LEN=120) STRING
	CHARACTER(LEN=120) OUT_FILE
        CHARACTER(LEN=10) DATA_TYPE
!
	INTEGER N_NODES
	INTEGER I,J
	INTEGER I1,I2
	INTEGER ILOW,IHIGH
	INTEGER LUIN
	INTEGER LUOUT
	REAL*4 X1,X2,Y1,Y2
	REAL*4 Y_SUM,XY_SUM
!
        REAL*4 EW,CENTROID
        REAL*4 XCUR(50),YCUR(50)
	REAL*4 SLOPE
	LOGICAL FILE_EXISTS
	INTEGER PLOT_ID,CURSERR
        CHARACTER(LEN=1) CURSVAL
!
	INTEGER GET_INDX_SP,PGCURS
	EXTERNAL GET_INDX_SP
	INTEGER, PARAMETER :: IONE=1
!
	LOCAL_OPTION=UPPER_CASE(INPUT_OPTION)
	IF(LOCAL_OPTION .EQ. 'CURSOR' .OR. LOCAL_OPTION .EQ. 'CURX')THEN
	  CALL GET_LU(LUOUT,'Output LU in PG_DEF_CONTINUUM')
	  OUT_FILE='CONT_NODES'
	  INQUIRE(FILE=TRIM(OUT_FILE),EXIST=FILE_EXISTS)
	  I=0
	  DO WHILE(FILE_EXISTS)
	    I=I+1
	    WRITE(OUT_FILE(11:13),'(A1,I2.2)')'_',I
	    INQUIRE(FILE=TRIM(OUT_FILE),EXIST=FILE_EXISTS)
	  END DO
	  OPEN(UNIT=LUOUT,FILE=TRIM(OUT_FILE),ACTION='WRITE',STATUS='NEW')
	END IF
!
! Get continuum locations as defined by cursor. CURSOR
! position is used to define continum.
!
	IF (LOCAL_OPTION .EQ. 'CURSOR')THEN
	  J=1; CALL PGSCI(J)
	  WRITE(LUOUT,'(T40,A)')'Number of data points'
	  WRITE(LUOUT,'(8X,A1,15X,A1,T40,A)')'X','Y','!Data type'
	  XCUR(1)=1; YCUR(1)=1
	  DO I=1,50
	    IF(I .NE. 1)THEN
	      XCUR(I)=XCUR(I-1)
	      YCUR(I)=YCUR(I-1)
	    END IF
	    CURSERR = PGCURS(XCUR(I),YCUR(I),CURSVAL)
	    IF(CURSVAL .EQ. 'E' .OR. CURSVAL .EQ. 'e')EXIT
	    CALL PGPT(IONE,XCUR(I),YCUR(I),IONE)
	    WRITE(6,'(2ES16.7)')XCUR(I),YCUR(I)
	    N_NODES=I
	    WRITE(6,*)CURSVAL
	  END DO
	  CLOSE(LUOUT)
	  WRITE(6,*)'Nodes written to ',TRIM(OUT_FILE)
!
! Use cursors to define bands to define the continuum. In each band,
! the continuum is defined by averaging the data. BANDS are written
! to CONT_NODES and these can be reread in by the DC option. Top of
! file will need to be updated with number of nodes.
!
	ELSE IF (LOCAL_OPTION .EQ. 'CURX')THEN
	  J=1; CALL PGSCI(J)
	  X1=1; Y1=1
	  WRITE(LUOUT,'(T40,A)')'Number of data points'
	  WRITE(LUOUT,'(8X,A1,15X,A1,T40,A)')'X1','X2','!Data type'
	  WRITE(LUOUT,'(A)')' '
	  WRITE(LUOUT,'(A)')'Use cursors to define continuum bands'
	  WRITE(LUOUT,'(A)')' '
	  DO I=1,50
	    CURSERR = PGCURS(X1,Y1,CURSVAL)
	    IF(CURSVAL .EQ. 'E' .OR. CURSVAL .EQ. 'e')EXIT
	    CALL PGPT(IONE,X1,Y1,IONE)
	    X2=X1; Y2=Y1
	    CURSERR = PGCURS(X2,Y2,CURSVAL)
	    IF(CURSVAL .EQ. 'E' .OR. CURSVAL .EQ. 'e')EXIT
	    CALL PGPT(IONE,X2,Y2,IONE)
	    N_NODES=I
	    I1=GET_INDX_SP(X1,CD(IPLT)%XVEC,NPTS(IPLT))
	    I2=GET_INDX_SP(X2,CD(IPLT)%XVEC,NPTS(IPLT))
	    Y_SUM=0.0;XY_SUM=0.0
	    DO J=MIN(I1,I2),MAX(I1,I2)
	      Y_SUM=Y_SUM+CD(IPLT)%DATA(J)
	      XY_SUM=XY_SUM+CD(IPLT)%DATA(J)*CD(IPLT)%XVEC(J)
	    END DO
	    XCUR(I)=XY_SUM/Y_SUM
	    YCUR(I)=Y_SUM/ABS(I2-I1+1)
	    WRITE(LUOUT,*)X1,X2
	    X1=X2; Y1=Y2
	  END DO
	  CLOSE(LUOUT)
	  WRITE(6,*)'Nodes written to ',TRIM(OUT_FILE)
!
	ELSE
	  CALL GET_LU(LUIN,'Input LU in PG_DEF_CONTINUUM')
	  OPEN(UNIT=LUIN,FILE=TRIM(INPUT_OPTION),STATUS='OLD',ACTION='READ',IOSTAT=RET_IOS)
	  IF(RET_IOS .NE. 0)THEN
	    WRITE(6,*)'Unable to open file in PG_DEF_CONT'
	    WRITE(6,*)'IOSTAT=',RET_IOS
	    WRITE(6,*)'FILE=',TRIM(INPUT_OPTION)
	    RETURN
	  END IF
!
	  N_NODES=0; DATA_TYPE=' '
	  DO WHILE(1 .EQ. 1)
	    READ(LUIN,'(A)',IOSTAT=IOS)STRING
	    IF(IOS .NE. 0)THEN
	      WRITE(6,*)'Error reading file ',TRIM(INPUT_OPTION)
	      WRITE(6,*)'Make sure you have a string with: !Data type'
	      WRITE(6,*)'Make sure you have a string with: !Number of data points'
	      WRITE(6,*)'IOS=',IOS
	      RETURN
	    END IF
	    IF( INDEX(STRING,'!Data type') .NE. 0 )THEN
	      IF( INDEX(STRING,' X ') .NE. 0 .AND. INDEX(STRING,' Y  ') .NE. 0)THEN
	        DATA_TYPE='NODES'
	      ELSE IF( INDEX(STRING,' X1 ') .NE. 0 .AND. INDEX(STRING,' X2 ') .NE. 0)THEN
	        DATA_TYPE='BANDS'
	      ELSE
	        WRITE(6,*)'Data type not recognized'
	        RETURN
	      END IF
	      WRITE(6,'(A,T30,A)')'Data type is',TRIM(DATA_TYPE)
	    ELSE IF(INDEX(STRING,'!Number of data points') .NE. 0)THEN
	      READ(STRING,*)N_NODES
	      WRITE(6,'(A,T28,I4)')'Number of nodes is',N_NODES
	    END IF
	    IF(N_NODES .NE. 0 .AND. DATA_TYPE .NE. ' ')EXIT
	  END DO
!
	  IF(DATA_TYPE .EQ. 'NODES')THEN
	    DO I=1,N_NODES
	      READ(LUIN,*)XCUR(I),YCUR(I)
	    END DO
	  ELSE
	    DO I=1,N_NODES
	      READ(LUIN,*)X1,X2
	      I1=GET_INDX_SP(X1,CD(IPLT)%XVEC,NPTS(IPLT))
	      I2=GET_INDX_SP(X2,CD(IPLT)%XVEC,NPTS(IPLT))
	      Y_SUM=0.0;XY_SUM=0.0
	      DO J=MIN(I1,I2),MAX(I1,I2)
	        Y_SUM=Y_SUM+CD(IPLT)%DATA(J)
	        XY_SUM=XY_SUM+CD(IPLT)%DATA(J)*CD(IPLT)%XVEC(J)
	      END DO
	      XCUR(I)=XY_SUM/Y_SUM
	      YCUR(I)=Y_SUM/ABS(I2-I1+1)
	    END DO
	    WRITE(6,*)'Defined node averages'
	  END IF
	  CLOSE(UNIT=LUIN)
	END IF
!
	IF(CURVE_OPTION .EQ. 'LINEAR')THEN
	  I2=GET_INDX_SP(XCUR(1),CD(IPLT)%XVEC,NPTS(IPLT))
	  DO I=1,N_NODES-1
	    WRITE(6,*)XCUR(I),XCUR(I+1)
	    I1=I2
	    I2=GET_INDX_SP(XCUR(I+1),CD(IPLT)%XVEC,NPTS(IPLT))
	    SLOPE=(YCUR(I+1)-YCUR(I))/(XCUR(I+1)-XCUR(I))
	    DO J=MIN(I1,I2),MAX(I1,I2)
               CONT(J)=YCUR(I)+SLOPE*(CD(IPLT)%XVEC(J)-XCUR(I))
	    END DO
	  END DO
!
! Fit a monotonic cubic to a set of data bands read from a file.
!
	ELSE
	  ILOW=NPTS(IPLT); IHIGH=-1
	  DO I=1,N_NODES
	    I1=GET_INDX_SP(XCUR(I),CD(IPLT)%XVEC,NPTS(IPLT))
	    ILOW=MIN(I1,ILOW); IHIGH=MAX(I1,IHIGH)
	  END DO
	  ILOW=ILOW+1; IHIGH=IHIGH-1; J=IHIGH-ILOW+1
	  WRITE(6,*)'ILH',ILOW,IHIGH
	  CALL MON_INTERP_SP(CONT(ILOW),J,IONE,CD(IPLT)%XVEC(ILOW),J,YCUR,N_NODES,XCUR,N_NODES)
	END IF
!
! If IP is non-zero, we store plot in a regular data vector.
!
	IF(OPLT .NE. 0)THEN
	  I=NPTS(IPLT)
	  IF(ALLOCATED(CD(OPLT)%XVEC) .AND. IPLT .NE. OPLT)THEN
	    DEALLOCATE (CD(OPLT)%XVEC)
	    DEALLOCATE (CD(OPLT)%DATA)
	  END IF
	  IF(IPLT .NE. OPLT)THEN
	    ALLOCATE (CD(OPLT)%XVEC(I),STAT=RET_IOS)
	    IF(RET_IOS .EQ. 0)ALLOCATE (CD(OPLT)%DATA(I),STAT=RET_IOS)
	    IF(RET_IOS .NE. 0)THEN
	      WRITE(6,*)'Error: unable to allocate new data vectors'
	      WRITE(6,*)'RET_IOS=',RET_IOS
	      RET_IOS=1
	      RETURN
	      END IF
	    END IF
	    CD(OPLT)%XVEC(1:I)=CD(IPLT)%XVEC(1:I)
	  END IF
	CD(OPLT)%DATA(1:I)=CONT(1:I)
	NPTS(OPLT)=I
	ERR(OPLT)=.FALSE.
	IF(OPLT .GT. NPLTS)NPLTS=OPLT
!
	RETURN
	END
