!
! Altered  05-Oct-2026: VINF consistency check uses peak speed (edited by Devin).
! Routine to read in the Radius grid to be utilized. The velocity and
! and sigma (dlnv/dlnr-1) must also be given.
!
! At present grid in must have the same number of data points as in CMFGEN,
! and must have R(ND)/R(1) (input) the same as RP/RMAX. While it is
! easy to interpolate V and SIGMA onto a new grid, it is not easy to
! decide on how a new R grid will be created.
!
	SUBROUTINE RD_RV_FILE_V2(R,V,SIGMA,RMAX,RP,VINF,LUIN,ND,OPTIONS,N_OPT,
	1                       NONMONO_V)
	USE SET_KIND_MODULE
!
! Altered  17-Apr-2024: Added check that R grid is monotonic.
! Altered  12-Jan-2024: Altered handling of SIGMA < 1
! Altered  19-Aug-2015: Added check that SIGMA < -1.0D0 (cur_hmi,21-Jun-2015).
! Altered  31-Mar-2008: Check on VINF is now set to 10% accuracy. Done as V(1) is normally
!                         < Vinf as radius grid does not extend to infinity.
! Altered  16-Jan-2007: VINF only checked of > 0.1 km/s (not important for pp models)
! Modified 02-Feb-2005: Blank lines and comments can now appear after ND string
!                         but before main data set.
! Modified 27-Apr-2000: VINF inserted in call.
! Modified 10-Feb-2000: Option string now utilized.
!                         RVSIG_COL format inserted.
!
	IMPLICIT NONE
!
	INTEGER ND
	REAL(KIND=LDP) R(ND)			!Radial grin in 10^10 cm
	REAL(KIND=LDP) V(ND)			!V in km/s
	REAL(KIND=LDP) SIGMA(ND)		!dlnV/dlnr-1
	REAL(KIND=LDP) RMAX			!Maximum radius  (in 10^10 cm)
	REAL(KIND=LDP) RP			!Core radius  (in 10^10 cm)
	REAL(KIND=LDP) VINF                     !Terminal velocity
!
! Only OPTIONS(1) is presently utilized. This indicates the format
! of the data in the file.
!
	INTEGER N_OPT
	CHARACTER*(*) OPTIONS(N_OPT)
!
! Local vectors.
!
	REAL(KIND=LDP) HT(ND-1)	
	REAL(KIND=LDP) VEL(ND-1)
	REAL(KIND=LDP) dVdR(ND-1)
!
! Local variables.
!
	REAL(KIND=LDP) T1
	INTEGER LUIN
	INTEGER LUER
	INTEGER I
	INTEGER IOS
	INTEGER ND_LOC
	INTEGER ERROR_LU
	INTEGER, PARAMETER :: IZERO=0
	LOGICAL NONMONO_V
	LOGICAL L_ERROR
	CHARACTER*132 STRING
!
	LUER=ERROR_LU()
!
! R, V, and SIGMA in column format, with simple header.
! As output by NEWRG in DISPGEN
!
	IF(OPTIONS(1) .EQ. 'RVSIG_COL')THEN
	  CALL GEN_ASCI_OPEN(LUIN,'RVSIG_COL','OLD',' ','READ',IZERO,IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error in RD_RV_FILE_V2: IOSTAT=',IOS
	    WRITE(LUER,*)'Unable to open RVSIG_COL'
	    STOP
	  END IF
	  STRING=' '
	  DO WHILE( INDEX(STRING,'!Number of depth points').EQ. 0)
	    READ(LUIN,'(A)')STRING
	  END DO
	  READ(STRING,*)ND_LOC
	  IF(ND_LOC .NE. ND)THEN
	    WRITE(LUER,*)'Error in RD_RV_FILE_V2'
	    WRITE(LUER,*)'Routine can''t yet handle a differnet number of depth points'
	    STOP
	  END IF
!
! Skip any further blank strings or comments.
!
	  STRING=' '
	  DO WHILE(STRING .EQ. ' ' .OR. STRING(1:1) .EQ. '!')
	    READ(LUIN,'(A)')STRING
	  END DO
	  BACKSPACE(LUIN)
!
	  DO I=1,ND
	    READ(LUIN,*,IOSTAT=IOS)R(I),V(I),SIGMA(I)
	    IF(IOS .NE. 0)THEN
	      WRITE(LUER,*)'Error in RD_RV_FILE_V2: IOSTAT=',IOS
	      WRITE(LUER,*)'Error reading R, V & SIGMA -- check correct # of data values'
	      STOP
	    END IF
	  END DO
	  CLOSE(LUIN)
!
	  DO I=1,ND-1
	    IF(R(I) .LE. R(I+1))THEN
	      WRITE(6,*)'Error for RD_RV_FILE_V2'
	      WRITE(6,*)'R grid must decrease monotonically with depth index'
	      WRITE(6,*)I,R(I),R(I+1)
	      STOP
	    END IF
	  END DO
!
	  DO I=1,ND-1
	    IF(V(I) .LE. V(I+1) .AND. V(I) .LT. 1.0_LDP)THEN
	      WRITE(6,*)'Warning from RD_RV_FILE_V2'
	      WRITE(6,*)'Velocity law in hydrostatic sone is not monotonic'
	      WRITE(6,*)'Starts at depth index (I) and velocities V(I),V(I+1):',I,V(I),V(I+1)
	      EXIT
	    ELSE IF(V(I) .LE. V(I+1))THEN
	      WRITE(LUER,*)'Warning in RD_RV_FILE_V2'
	      WRITE(6,*)'CMFGEN maynot handle a non-monotonic velocity law'
	      WRITE(6,*)'Starts at depth index (I) and velocities V(I),V(I+1):',I,V(I),V(I+1)
!	      STOP
	    END IF
	  END DO
!
! Now scale to match radius. V is assumed to be fixed on a r/R(ND) scale.
! SIGMA does not change with a simple scaling in r.
!
	  T1=R(ND)
	  R(1:ND)=R(1:ND)*(RP/T1)
	  IF(ABS(R(1)-RMAX) .GT. 0.1_LDP*(R(1)-R(2)))THEN
	    WRITE(LUER,*)'Error in RD_RV_FILE_V2'
	    WRITE(LUER,*)'R(1)/R(ND) must match RMAX/RP'
	    WRITE(LUER,*)'R(1)/R(ND)=',R(1)/R(ND)
	    STOP
          END IF	
	  R(1)=RMAX
!
	ELSE IF(OPTIONS(1) .EQ. 'deKOTER')THEN
C
C Read data from Alex Dekoter''s program.
C
	  CALL GEN_ASCI_OPEN(LUIN,'deKOTER','OLD',' ','READ',IZERO,IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error in RD_RV_FILE_V2: IOSTAT=',IOS
	    WRITE(LUER,*)'Unable to open RVSIG_COL'
	    STOP
	  END IF
C
	  STRING=' '
	  DO WHILE( INDEX(STRING,'STRUCTURE, GRID SIZE') .EQ. 0)
	    READ(LUIN,'(A)')STRING
	  END DO
	  READ(LUIN,'(A)')STRING
	  READ(LUIN,'(A)')STRING
	  READ(LUIN,'(A)')STRING
	  READ(LUIN,*)ND_LOC
!
	  IF(ND_LOC .NE. ND-1)THEN
	    WRITE(LUER,*)'Error in RD_RV_FILE_V2'
	    WRITE(LUER,*)'Routine can''t yet handle a differnet number of depth points'
	    WRITE(LUER,*)'ND in file should be 1 less than in CMFGEN'
	    STOP
	  END IF
!
	  STRING=' '
	  DO WHILE( INDEX(STRING,'[RCORE]') .EQ. 0)
	    READ(LUIN,'(A)')STRING
	  END DO
	  READ(LUIN,'(A)')STRING
C
	  DO I=ND_LOC,1,-1
	    READ(LUIN,'(A)')STRING
	    READ(STRING(8:20),*)HT(I)
	    READ(STRING(58:68),*)VEL(I)
	    READ(STRING(70:81),*)dVdR(I)
	  END DO
	  CLOSE(LUIN)
C
	  R(1)=(HT(1)+1.0_LDP)*RP
	  SIGMA(1)=dVDR(1)*(HT(1)+1.0_LDP)/VEL(1)-1.0_LDP
	  V(1)=VEL(1)
	  DO I=3,ND
	    R(I)=(HT(I-1)+1.0_LDP)*RP
	    SIGMA(I)=dVDR(I-1)*(HT(I-1)+1.0_LDP)/VEL(I-1)-1.0_LDP
	    V(I)=VEL(I-1)
	  END DO
!
	  CLOSE(LUIN)
! Insert extra point near outer boundary, as use first order
! boundary condition.
!
	  R(2)=R(1)-0.1_LDP*(R(1)-R(3))
	  V(2)=V(1)+0.1_LDP*(V(1)-V(3))
	  SIGMA(2)=SIGMA(1)+0.1_LDP*(SIGMA(1)-SIGMA(3))
	ELSE
	  WRITE(LUER,*)'Error in RD_RV_FILE_V2 -- OPTION not recognized'
	  WRITE(LUER,*)OPTIONS(1)
	  STOP
	END IF
!
! NONMONO_V is set only by DO_NONMONO_SOB/DO_NONMONO_KNOP.  VINF is a
! characteristic speed, not necessarily the outer boundary value: the
! fastest gas can be an interior shell, or the profile can reverse sign.
!
	IF(NONMONO_V)THEN
	  T1=MAXVAL(ABS(V(1:ND)))
	  L_ERROR=.FALSE.
	  IF(VINF .GT. 0.1_LDP)THEN
	    IF(T1 .LE. 0.0_LDP)THEN
	      L_ERROR=.TRUE.
	    ELSE IF(ABS(T1-VINF)/T1 .GT. 0.1_LDP)THEN
	      L_ERROR=.TRUE.
	    END IF
	  END IF
	  IF(L_ERROR)THEN
	    WRITE(LUER,*)'VINF in VADAT file and velocity file are inconsistent'
	    WRITE(LUER,*)'MAXVAL(ABS(V))=',T1
	    WRITE(LUER,*)'VINF',VINF
	    STOP
	  END IF
	ELSE IF( ABS(V(1)-VINF)/V(1) .GT. 0.1_LDP .AND. VINF .GT. 0.1_LDP)THEN
	  WRITE(LUER,*)'VINF in VADAT file and velocity file are inconsistent'
	  WRITE(LUER,*)'V(1)=',V(1)
	  WRITE(LUER,*)'VINF',VINF
	  STOP
	END IF
!
	RETURN
	END
