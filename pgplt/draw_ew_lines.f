!
! General purpose iroutine designed for emission line stars.
! Program draws lines, with a length scaled by the line EW, to identify
! contibution to emission lines.
!
	SUBROUTINE DRAW_EW_LINES(XPAR,YPAR,EXPCHAR,T_OUT)
	USE SET_KIND_MODULE
	USE LINE_ID_MOD
	USE MOD_COLOR_PEN_DEF
	IMPLICIT NONE
!
! Created 28-Aug-2018
!
	REAL*4 XPAR(2)
	REAL*4 YPAR(2)
	REAL*4 XSTRPOS,YSTRPOS
	REAL*4 XCHAR_SIZE,YCHAR_SIZE
	REAL*4 EXPCHAR
	INTEGER T_OUT
!
	REAL(KIND=LDP), ALLOCATABLE :: DP_VEC(:)
	INTEGER, ALLOCATABLE :: INDX(:)
!
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: ITWO=2
	INTEGER, PARAMETER :: ITHREE=4
	INTEGER, PARAMETER :: IFOUR=4
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
!
	INTEGER LUIN
	INTEGER IOS
	INTEGER I,J,K,L
	REAL*4 T1,T2
	REAL*4 XSTART,XEND	
	REAL*4 XPOS,YPOS,ORIENT
	REAL*4 LOC_VEC_END
	CHARACTER(LEN=80) TMP_STR
!
	INTEGER NSPEC
	INTEGER, PARAMETER :: NSPEC_MAX=20
	CHARACTER(LEN=6) SPEC(NSPEC_MAX)
	INTEGER PEN(NSPEC_MAX)
!
	CALL GET_LU(LUIN,'Called in DRAW_EW_LINES')
	OPEN(UNIT=LUIN,FILE='SPECIES_LIST',STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	IF(IOS .EQ. 0)THEN
	  J=0
	  DO WHILE(1 .EQ. 1)
	    READ(LUIN,'(A)',END=100)SPEC(J+1)
	    J=J+1
	    SPEC(J)=ADJUSTL(SPEC(J))
	    K=INDEX(SPEC(J),' ')
	    IF(K .NE. LEN_TRIM(SPEC(J))+1)THEN
	      READ(SPEC(J),*)PEN(J)
	      SPEC(J)=SPEC(J)(1:K)
	    ELSE
	     PEN(J)=J+1
	    END IF
	  END DO
100	  CONTINUE
	  CLOSE(LUIN)
	  NSPEC=J
	ELSE
	  J=0
	  J=J+1; SPEC(J)='HeII'
	  J=J+1; SPEC(J)='CV'
	  J=J+1; SPEC(J)='CIV'
	  J=J+1; SPEC(J)='CIII'
	  J=J+1; SPEC(J)='OVI'
	  J=J+1; SPEC(J)='OV'
	  J=J+1; SPEC(J)='OIV'
	  J=J+1; SPEC(J)='OIII'
	  J=J+1; SPEC(J)='NeVI'
	  J=J+1; SPEC(J)='NeV'
	  J=J+1; SPEC(J)='NeIV'
	  J=J+1; SPEC(J)='NeIII'
	  J=J+1; SPEC(J)='NeII'
	  NSPEC=J
	  DO K=1,NSPEC
	    PEN(K)=K+1
	  END DO
	END IF
!
	
	WRITE(6,*)'Entered DRAW_EW_LINES'
	IF(N_EW_IDS .NE. 0)THEN
	  CALL PGSCI(IONE)
	  ID_LOC=4
	  ID_ORIENT=90.0D0
	  ID_LOC_PG=1.0D0
	  ID_VEC_END=0.05
!
! This section is left over from DRAW_LINE_IDs. It not utilized,
! but could be easily used.
!
	  DO I=1,N_EW_IDS
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
! Order the EWs so they are numerically decreasing. Assuming they are
! positive, this guaranettes that lines with small EW's are draw on top
! of lines with larger EW's, thus makig them visible.
!
	  ALLOCATE(DP_VEC(N_EW_IDS))
	  ALLOCATE(INDX(N_EW_IDS))
	  DP_VEC=ID_EW(1:N_EW_IDS)
	  CALL INDEXX(N_EW_IDS,DP_VEC,INDX,L_FALSE)

	  I=2
	  CALL PGSLW(I)
	  DO L=1,N_EW_IDS
	    I=INDX(L)
	    WRITE(70,*)I,ID_EW(I)
	    J=0
	    DO K=1,NSPEC
	      IF(INDEX(LINE_ID(I),TRIM(SPEC(K))) .NE. 0)THEN
	        J=PEN(K)
	        EXIT
	      END IF
	    END DO
	    IF(J .EQ. 0)THEN
	      J=1
	      CALL PGSLS(ITWO)		!Dashed line
	    ELSE
	      CALL PGSLS(IONE)
	    END IF
!
	    CALL PGSCI(J)
	    LOC_VEC_END=YPAR(1)+ID_VEC_END*(YPAR(2)-YPAR(1))
	    ID_VEC_BEG=LOC_VEC_END+(ID_EW(I)/EW_SCALE_FAC)*(YPAR(2)-YPAR(1))/4.0D0
	    WRITE(6,*)LOC_VEC_END,ID_VEC_BEG
	    CALL PGMOVE(ID_WAVE(I),ID_VEC_BEG)
	    CALL PGDRAW(ID_WAVE(I),LOC_VEC_END)
	  END DO
	  DEALLOCATE(DP_VEC,INDX)
!
! Draw key to color code on right hand side of plot.
!
	  WRITE(6,*)BLUE_PEN
	  WRITE(6,*)'If you cannot see the species labels, use N option'
	  WRITE(6,*)' to change right handed extent of plot'
	  WRITE(6,*)DEF_PEN
!
	  XSTART=XPAR(2)+0.05*(XPAR(2)-XPAR(1))
	  XEND=XPAR(2)+0.1*(XPAR(2)-XPAR(1))
	  XPOS=XPAR(2)+0.12*(XPAR(2)-XPAR(1))
	  YPOS=YPAR(2)
	  ORIENT=0.0
	  J=4
	  I=0;CALL PGSCLP(I)
	  CALL PGSCH(EXPCHAR*ID_EXPCHAR)
	  DO K=1,NSPEC
	    CALL PGSCI(PEN(K))
	    YPOS=YPOS-0.07*(YPAR(2)-YPAR(1))
	    CALL PGMOVE(XSTART,YPOS)
	    CALL PGDRAW(XEND,YPOS)
	    TMP_STR=SPEC(K)
	    CALL JUSTIFY_CONVERT_V2(XPOS,YPOS,J,T1,ORIENT,L_TRUE,XSTRPOS,YSTRPOS,TMP_STR,IONE)
	    CALL PGPTXT(XSTRPOS,YSTRPOS,ORIENT,T1,TMP_STR)
	  END DO
	  I=1;CALL PGSCLP(I)
	END IF
!
	RETURN
	END
