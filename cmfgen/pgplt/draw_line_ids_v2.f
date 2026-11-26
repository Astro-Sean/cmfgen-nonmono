!
! General purpose line plotting routine to label lines is a spectrum.
!
	SUBROUTINE DRAW_LINE_IDS_V2(XPAR,YPAR,EXPCHAR,IP,T_OUT)
	USE SET_KIND_MODULE
	USE MOD_CURVE_DATA
	USE LINE_ID_MOD
	USE MOD_COLOR_PEN_DEF
	IMPLICIT NONE
!
! Altered 04-Mar-2023 : Now adjustble number of decimal digits (passed in LIND_ID_MOD).
!                         Many other changes donw to improve code.
! Altered 10-Oct-2022 : Changed NPOS limit, and avoid lINE check when outside limits.
! Altered 05-Sep-2022 : Changed labeling algorithm slightly. Works better but not efficient/perfect.
! Altered 22-Apr-2020 : Labeling algorithim altered to prevent overlap. Label loction may
!                         not be optimal.
!
	REAL*4 XPAR(2)
	REAL*4 YPAR(2)
	REAL*4 XSTRPOS,YSTRPOS
	REAL*4 XCHAR_SIZE,YCHAR_SIZE
	REAL*4 EXPCHAR
	REAL*4 LOC_ID_VEC_BEG,LOC_ID_VEC_END
!
	INTEGER T_OUT
	INTEGER IP
	LOGICAL, PARAMETER :: TRACE=.FALSE.              !Set to TRUE for debugging purposes
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
	INTEGER I,J,K,L
	INTEGER CNT
	REAL*4 T1,T2
	REAL*4 XLOW,XUP
	CHARACTER(LEN=80) TMP_STR
	CHARACTER(LEN=10) TMP_FMT
	INTEGER GET_INDX_SP
	EXTERNAL GET_INDX_SP
	LOGICAL CHANGE_MADE
!
! We do nothing if no lines have been read in.
!
	IF(N_LINE_IDS .NE. 0)THEN
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
!
! We only check lines that lie in the FULL spectral window.
! This is necessary forobserved data.
!
	    IF( (ID_WAVE(I)-CD(IP)%XVEC(1))*
	1          (CD(IP)%XVEC(NPTS(IP))-ID_WAVE(I)) .LE. 0)THEN
	      WR_ID(I)=.FALSE.
	    END IF
	  END DO
!
	  IF(LINE_CUT_PARAM .GT. 0.0)THEN
	    DO I=1,N_LINE_IDS
	      IF(WR_ID(I))THEN
	        J=GET_INDX_SP(ID_WAVE(I),CD(IP)%XVEC,NPTS(IP))
	        IF(ABS(CD(IP)%DATA(J)-1.0D0) .LT. LINE_CUT_PARAM)THEN
	          WR_ID(I)=.FALSE.
	        END IF
	      END IF
	    END DO
	  END IF
!
! Set the pen and character size.
!
	  CALL PGQCI(ID_LINE_PEN_SAVE)
	  CALL PGSCI(ID_LINE_PEN)
	  CALL PGSCH(EXPCHAR*ID_EXPCHAR)
	  CALL PGQCS(IFOUR,XCHAR_SIZE,YCHAR_SIZE)
	  IF(TRACE)WRITE(T_OUT,*)'XCHAR_SIZE=',XCHAR_SIZE
!
! Determine maximum number of label slots.
!
	  LAB_SIZE=1.1D0*XCHAR_SIZE
	  LAB_START=XPAR(1)+2*LAB_SIZE
	  NPOS=ABS((XPAR(2)-XPAR(1))/LAB_SIZE)-4
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
	  XLOW=XPAR(1)+2*LAB_SIZE
	  XUP=XPAR(2)-2*LAB_SIZE
!
! Check if there are any identified lines in the current window. If not, we exit the
! program. This can occur if we have read in lines, but then changed the limits on the
! spectral window.
!
	  LOC_NLINES=0
	  DO I=1,N_LINE_IDS
	    T1=(ID_WAVE(I)-XLOW)*(XUP-ID_WAVE(I))
	    IF(WR_ID(I) .AND. T1 .GT. 0.0D0)LOC_NLINES=LOC_NLINES+1
	    IF(T1 .LT. 0)WR_ID(I)=.FALSE.
	  END DO
!
	  WRITE(6,*)RED_PEN
	  WRITE(6,*)'                                  Number of label slots is:',NPOS
	  WRITE(6,*)'   Number of lines in spectral window to be indentified is:',LOC_NLINES
	  WRITE(6,*)DEF_PEN
!
	  IF(LOC_NLINES .EQ. 0)RETURN
	  IF(LOC_NLINES .GT. NPOS)THEN
	    WRITE(6,*)'Error -- too many lines to label in plot window'
	    WRITE(6,*)'Change line selection parameters or label size'
	    RETURN
	  END IF
!
! Store lines whose ID will be written. Initially they occupy
! all slots from location 3 up to LOC_NLINES+2.
!
	  L=2; T2=2*LAB_SIZE
	  DO I=1,N_LINE_IDS-1
	    IF(WR_ID(I))THEN
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
!
! J can be > NPOS because I don't use the last two slots because of
! possile ovelap with axes markers.
!
	      ELSE IF(J .GT. NPOS+2)THEN
	        WRITE(6,*)'Possible error'
	        WRITE(6,*)J,NPOS
	        WRITE(6,*)ID_WAVE(LAB_ID(I)),XPAR(1),XPAR(2)
	      END IF
	    END IF
	  END DO
!
! The previous lop has a bug, and can mix up line ordering. This fixes
! that.
!
	 DO I=1,NPOS-1
	   K=1
	   DO WHILE(LAB_ID(I) .NE. 0 .AND. K .NE. 0)
	     K=0
	     DO J=NPOS,I+1,-1
	       IF(LAB_ID(J) .NE. 0)THEN
	         IF(ID_WAVE(LAB_ID(I)) .GT. ID_WAVE(LAB_ID(J)))THEN
	           K=LAB_ID(J)
	           LAB_ID(J)=LAB_ID(I)
	           LAB_ID(I)=K
	           EXIT
	         END IF
	       END IF
	     END DO
	   END DO
	 END DO
!
! These two loops move lines closer (when possible) to ther correct
! locations.
!
	 CHANGE_MADE=.TRUE.; CNT=0
	 DO WHILE(CHANGE_MADE .AND. CNT .LT. 10)
	   CHANGE_MADE=.FALSE.; CNT=CNT+1
	   DO I=1,NPOS-1
	     IF(LAB_ID(I) .NE. 0 .AND. LAB_ID(I+1) .EQ. 0)THEN
	       IF(ID_WAVE(LAB_ID(I)) .GT. LAB_POS(I+1)-0.5*LAB_SIZE)THEN
	         LAB_ID(I+1)=LAB_ID(I)
	         LAB_ID(I)=0
	         CHANGE_MADE=.TRUE.
	       END IF
	     END IF
	   END DO
	 END DO
!
	 CHANGE_MADE=.TRUE.; CNT=0
	 DO WHILE(CHANGE_MADE .AND. CNT .LT. 10)
	   CHANGE_MADE=.FALSE.; CNT=CNT+1
	   DO I=NPOS,2,-1
	     IF(LAB_ID(I) .NE. 0 .AND. LAB_ID(I-1) .EQ. 0)THEN
	       IF(ID_WAVE(LAB_ID(I)) .LT. LAB_POS(I-1)+0.5*LAB_SIZE)THEN
	         LAB_ID(I-1)=LAB_ID(I)
	         LAB_ID(I)=0
	         CHANGE_MADE=.TRUE.
	       END IF
	     END IF
	   END DO
	 END DO
!
	  IF(TRACE)THEN
	    WRITE(6,*)'Done spread out lines'; FLUSH(UNIT=6)
	    J=0
	    DO I=1,NPOS
	      IF(LAB_ID(I) .NE. 0)THEN
	        J=J+1
	        WRITE(6,'(1X,2I6,2F12.2,3X,2X,A)')I,LAB_ID(I),LAB_POS(I),ID_WAVE(LAB_ID(I)),TRIM(LINE_ID(LAB_ID(I)))
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
	      IF(ID_WAVE(LAB_ID(I)) .GT. LAB_POS(I-1) .AND. ID_WAVE(LAB_ID(I)) .LT. LAB_POS(I+1))THEN
	        LAB_POS(I)=ID_WAVE(LAB_ID(I))
	      ELSE IF(ID_WAVE(LAB_ID(I)) .LT. LAB_POS(I-1))THEN
	        LAB_ID(I-1)=LAB_ID(I)
	        LAB_ID(I)=0
	      ELSE IF(ID_WAVE(LAB_ID(I)) .GT. LAB_POS(I+1))THEN
	        LAB_ID(I+1)=LAB_ID(I)
	        LAB_ID(I)=0
	      END IF
	    END IF
	  END DO
!
	  IF(TRACE)THEN
	    J=0
	    DO I=1,NPOS
	      IF(LAB_ID(I) .NE. 0)THEN
	        J=J+1
	        WRITE(6,'(1X,2I6,2F12.2,3X,2X,A)')I,LAB_ID(I),LAB_POS(I),ID_WAVE(LAB_ID(I)),TRIM(LINE_ID(LAB_ID(I)))
	      END IF
	    END DO
	    FLUSH(UNIT=6)
	  END IF
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
	    T2=(ID_WAVE_OFF(I)-XPAR(1))*(XPAR(2)-ID_WAVE_OFF(I))
	    IF(WR_ID(I) .AND. T1 .GT. 0 .AND. T2. GT. 0)THEN
	      IF(USE_DEF_OFFSET)THEN
	        LOC_ID_VEC_BEG=ID_VEC_BEG
	        LOC_ID_VEC_END=ID_VEC_END
	      ELSE
	        LOC_ID_VEC_BEG=ID_Y_BEG(I)
	        LOC_ID_VEC_END=ID_Y_END(I)
	      END IF
!
	      TMP_FMT='(F12.'
	      WRITE(TMP_FMT(6:7),'(I1,A1)')NO_DEC_DIGITS,')'
	      TMP_STR=' '
	      WRITE(TMP_STR,TMP_FMT)ID_WAVE(I)
	      TMP_STR=TRIM(LINE_ID(I))//'-'//ADJUSTL(TMP_STR)
!	      IF(.NOT. OBSERVED_WAVE(I))TMP_STR='*'//TMP_STR
!
	      IF(ID_VEC_BEG .LT. ID_VEC_END)THEN
	        T1=ID_VEC_END+0.2*(ID_VEC_END-ID_VEC_BEG)
	        CALL JUSTIFY_CONVERT_V2(ID_WAVE_OFF(I),T1,ID_LOC,ID_LOC_PG,ID_ORIENT,.TRUE.,
	1                  XSTRPOS,YSTRPOS,TMP_STR,IONE)
!
	        T1=(XSTRPOS-XPAR(1))*(XPAR(2)-XSTRPOS)
	        T2=(YSTRPOS-YPAR(1))*(YPAR(2)-YSTRPOS)
	        IF(T1 .GT. 0 .AND. T2 .GT. 0)THEN
	          CALL PGSCH(EXPCHAR*ID_EXPCHAR)
	          CALL PGPTXT(XSTRPOS,YSTRPOS,ID_ORIENT,ID_LOC_PG,TMP_STR)
	          CALL PGMOVE(ID_WAVE_OFF(I),LOC_ID_VEC_END)
	          CALL PGDRAW(ID_WAVE(I),LOC_ID_VEC_BEG)
	        END IF
	      ELSE
	        ID_LOC=6
	        T1=ID_VEC_END-0.2*(ID_VEC_BEG-ID_VEC_END)
	        CALL JUSTIFY_CONVERT_V2(ID_WAVE_OFF(I),T1,ID_LOC,ID_LOC_PG,ID_ORIENT,.TRUE.,
	1                  XSTRPOS,YSTRPOS,TMP_STR,IONE)
	        ID_LOC=4
	        T1=(XSTRPOS-XPAR(1))*(XPAR(2)-XSTRPOS)
	        T2=(YSTRPOS-YPAR(1))*(YPAR(2)-YSTRPOS)
	        IF(T1 .GT. 0 .AND. T2 .GT. 0)THEN
	          CALL PGSCH(EXPCHAR*ID_EXPCHAR)
	          CALL PGPTXT(XSTRPOS,YSTRPOS,ID_ORIENT,ID_LOC_PG,TMP_STR)
	          CALL PGMOVE(ID_WAVE_OFF(I),LOC_ID_VEC_END)
	          CALL PGDRAW(ID_WAVE(I),LOC_ID_VEC_BEG)
	        END IF
	      END IF
	    END IF
	  END DO
!
	  DEALLOCATE (LAB_ID,LAB_POS)
	END IF
	CALL PGSCI(ID_LINE_PEN_SAVE)
!
	RETURN
	END
