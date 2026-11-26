!
! General purpose line plotting routine to label lines is a spectrum.
!
	SUBROUTINE DRAW_LINE_IDS(XPAR,YPAR,EXPCHAR,T_OUT)
	USE SET_KIND_MODULE
	USE LINE_ID_MOD
	USE MOD_COLOR_PEN_DEF
	IMPLICIT NONE
!
! Altered 22-Apr-2020 : Labeling algorithim altered to prevent overlap. Label loction may
!                         not be optimal.
!
	REAL*4 XPAR(2)
	REAL*4 YPAR(2)
	REAL*4 XSTRPOS,YSTRPOS
	REAL*4 XCHAR_SIZE,YCHAR_SIZE
	REAL*4 EXPCHAR
	INTEGER T_OUT
	LOGICAL, PARAMETER :: TRACE=.FALSE.          !TRUE.              !Set to TRUE for debugging purposes
!
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: ITWO=2
	INTEGER, PARAMETER :: ITHREE=4
	INTEGER, PARAMETER :: IFOUR=4
!
! These local varable are used for defining the label locations.
! We split the X directon into NPOS label position, each of width LAB_SIZE.
! Initally the labels are crowded together on theleft. We the sequentially
! move the labels so that they are closer to the fetaure they are trying to identify.
!
! LAB_POS is the label positon in world coorodinates.
! LAB_ID provides a link between the label positions and the line to which it refers.
!
	INTEGER LOC_NLINES
	INTEGER NPOS
	REAL*4 LAB_SIZE
	REAL*4 LAB_START
	REAL*4, ALLOCATABLE :: LAB_POS(:)
	INTEGER, ALLOCATABLE :: LAB_ID(:)
!
	INTEGER IOS
	INTEGER I,J,K,L
	REAL*4 T1,T2
	CHARACTER(LEN=80) TMP_STR
!
! We do nothing if no lines have been read in.
!
	IF(N_LINE_IDS .NE. 0)THEN
	  CALL PGSCI(IONE)
	  ID_LOC=4
	  ID_ORIENT=90.0D0
	  ID_LOC_PG=1.0D0
!
! Check which species are being identified.
!
	  DO I=1,N_LINE_IDS
	    WR_ID(I)=.TRUE.
	    DO J=1,N_OMIT_ID
	      IF( INDEX(LINE_ID(I),TRIM(OMIT_ID(J))//' ') .NE. 0)THEN
	        WR_ID(I)=.FALSE.
	        EXIT
	      END IF
	    END DO
	    DO J=1,N_INC_ID
	      IF(J .EQ. 1)WR_ID(I)=.FALSE.
	      IF( INDEX(LINE_ID(I),TRIM(INC_ID(J))//' ') .NE. 0 )THEN
	        WR_ID(I)=.TRUE.
	        EXIT
	      END IF
	    END DO
	  END DO
!
! Set the character size.
!
	  CALL PGSCH(EXPCHAR*ID_EXPCHAR)
	  CALL PGQCS(IFOUR,XCHAR_SIZE,YCHAR_SIZE)
	  IF(TRACE)WRITE(T_OUT,*)'XCHAR_SIZE=',XCHAR_SIZE
!
! Check if there are any identified lines in the current window. If not, we exit the
! program. This can occur if we have read in lines, but then changed the limits on the
! spectral window.
!
	  LOC_NLINES=0
	  DO I=1,N_LINE_IDS-1
	    T1=(ID_WAVE(I)-XPAR(1))*(XPAR(2)-ID_WAVE(I))
	    IF(WR_ID(I) .AND. T1 .GT. 0.0D0)LOC_NLINES=LOC_NLINES+1
	  END DO
	  WRITE(6,*)'Number of lines in spectral window to be indentified is:',LOC_NLINES
	  IF(LOC_NLINES .EQ. 0)RETURN
!
! Determine maximum number of label slots.
!
	  LAB_SIZE=1.1D0*XCHAR_SIZE
	  LAB_START=XPAR(1)+1.5*LAB_SIZE
	  NPOS=ABS( (XPAR(2)-XPAR(1))/LAB_SIZE )
	  WRITE(6,*)'Number of label slots is:',NPOS
	  IF(LOC_NLINES .GT. 0.8*NPOS)THEN
	    WRITE(6,*)'Error -- too many lines to label in plot window'
	    WRITE(6,*)'Change line selectrion parameters or label size'
	    RETURN
	  END IF
!
	  IF(ALLOCATED(LAB_POS))DEALLOCATE(LAB_POS,LAB_ID)
          ALLOCATE (LAB_POS(NPOS),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (LAB_ID(NPOS),STAT=IOS)
	  IF(IOS .EQ. 0)THEN
	    IF(TRACE)WRITE(6,*)'Successfully allocated LAB_POS and LAB_ID in DRAW_LINE_IDS'
	  ELSE
	    WRITE(6,*)'Unable to allocate LAB_POS or LAB_ID in DRAW_LINE_IDS'
	    WRITE(6,*)'ERROR=',IOS
	  END IF
	  LAB_ID=0; LAB_POS=0
!
! Define label positions and initialize label links.
!
	  DO I=1,NPOS
	    LAB_POS(I)=LAB_START+I*LAB_SIZE
	    LAB_ID(I)=0
	  END DO
	  IF(TRACE)WRITE(6,*)'Set LAB_POS and LAB_ID'
!
! Store lines whose ID will be written. Initially they occupy
! all slots from location 3 up to LOC_NLINES+2.
!
	  L=2; T2=2*LAB_SIZE
	  DO I=1,N_LINE_IDS-1
	    T1=(ID_WAVE(I)-XPAR(1)-T2)*(XPAR(2)+T2-ID_WAVE(I))
	    IF(WR_ID(I) .AND. T1 .GT. 0.0D0)THEN
	      IF(L+1 .GT. NPOS)THEN
	        L=L-1
	        EXIT
	      END IF
	      L=L+1
	      LAB_ID(L)=I
!	      IF(L .LT. 10)WRITE(6,*)I,ID_WAVE(I)
	    END IF
	  END DO
	  IF(TRACE)THEN
	    WRITE(6,*)'Stored lines'
	    WRITE(6,*)NPOS,LAB_SIZE,XPAR(1)
	    WRITE(6,*)LAB_START
	  END IF
!
! Spread lines out where possible.
!
	  DO I=NPOS-2,3,-1
	    IF(LAB_ID(I) .NE. 0)THEN
	      J=(ID_WAVE(LAB_ID(I))-LAB_START)/LAB_SIZE+1
	      IF(J .LE. NPOS)THEN
	        DO WHILE(LAB_ID(J) .NE. 0)
	          J=J-1
	          IF(J .LT. 1)EXIT
	        END DO
	        IF(J .LT. 1)EXIT
	        LAB_ID(J)=LAB_ID(I)
	        IF(I .NE. J)LAB_ID(I)=0
	      ELSE
	        WRITE(6,*)'Possible error'
	        WRITE(6,*)J,NPOS
	        WRITE(6,*)ID_WAVE(LAB_ID(I)),XPAR(1),XPAR(2)
	      END IF
	    END IF
	  END DO
!
	  IF(TRACE)THEN
	    WRITE(6,*)'Done spread out lines'; FLUSH(UNIT=6)
	    J=0
	    DO I=1,NPOS
	      WRITE(6,*)LAB_ID(I); FLUSH(UNIT=6)
	      IF(LAB_ID(I) .NE. 0)THEN
	        J=J+1
	        WRITE(6,*)LAB_ID(I); FLUSH(UNIT=6)
	        WRITE(6,*)ID_WAVE(LAB_ID(I)); FLUSH(UNIT=6)
	        WRITE(6,*)TRIM(LINE_ID(LAB_ID(I))); FLUSH(UNIT=6)
	        WRITE(6,'(I7,F12.4,A)')LAB_ID(I),ID_WAVE(LAB_ID(I)),TRIM(LINE_ID(LAB_ID(I)))
	        FLUSH(UNIT=6)
	      END IF
	    END DO
	    WRITE(6,*)'Number of non zero LAB_IDs is',J
	  END IF
!
! For isolated lines, we make sure that the label location is centerd on the line.
!
	  DO I=2,NPOS-1
	    IF(LAB_ID(I-1) .EQ. 0 .AND. LAB_ID(I) .NE. 0 .AND. LAB_ID(I+1) .EQ. 0)THEN
	      LAB_POS(I)=ID_WAVE(LAB_ID(I))
	    END IF
	  END DO
!
! Now set the label positions so the strings cann be written out.
!
	  ID_WAVE_OFF=XPAR(1)-(XPAR(2)-XPAR(1))
	  DO I=1,NPOS
	    IF(LAB_ID(I) .NE. 0)ID_WAVE_OFF(LAB_ID(I))=LAB_POS(I)
	  END DO
!
	  DO I=1,N_LINE_IDS
	    T1=(ID_WAVE(I)-XPAR(1))*(XPAR(2)-ID_WAVE(I))
	    IF(WR_ID(I) .AND. T1 .GT. 0)THEN
	      T1=ID_SCL*ID_Y_OFF(I)
	      ID_VEC_BEG=0.9*(T1-1.0D0)+1.0D0
	      TMP_STR=' '
	      WRITE(TMP_STR,'(F12.2)')ID_WAVE(I)
	      TMP_STR=TRIM(LINE_ID(I))//'-'//ADJUSTL(TMP_STR)
	      IF(.NOT. OBSERVED_WAVE(I))TMP_STR='*'//TMP_STR
	      CALL JUSTIFY_CONVERT_V2(ID_WAVE_OFF(I),T1,ID_LOC,ID_LOC_PG,ID_ORIENT,.TRUE.,
	1                  XSTRPOS,YSTRPOS,TMP_STR,IONE)
	      T1=(XSTRPOS-XPAR(1))*(XPAR(2)-XSTRPOS)
	      T2=(YSTRPOS-YPAR(1))*(YPAR(2)-YSTRPOS)
	      IF(T1 .GT. 0 .AND. T2 .GT. 0)THEN
	        CALL PGSCI(ITWO)
	        CALL PGSCH(EXPCHAR*ID_EXPCHAR)
	        CALL PGPTXT(XSTRPOS,YSTRPOS,ID_ORIENT,ID_LOC_PG,TMP_STR)
	        T1=ID_SCL*ID_Y_OFF(I)
	        CALL PGMOVE(ID_WAVE_OFF(I),ID_VEC_BEG)
	        CALL PGDRAW(ID_WAVE(I),ID_VEC_END)
	      END IF
	    END IF
	  END DO
!
	  DEALLOCATE (LAB_ID,LAB_POS)
	END IF
!
	RETURN
	END
