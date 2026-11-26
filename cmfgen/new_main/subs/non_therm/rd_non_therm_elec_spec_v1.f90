!
! Subroutine to read in the NON_THERM_DEGRADATION_SPEC
!
	SUBROUTINE RD_NON_THERM_ELEC_SPEC_V1(ND,LU)
	USE SET_KIND_MODULE
	USE MOD_NON_THERM
!
! Created 30-Oct-2021
!
	INTEGER ND  		!Number of depth points
	INTEGER LU
	REAL(KIND=LDP) LOC_XKT(NKT)
	CHARACTER(LEN=132) STRING
!
! Local variables.
!
	INTEGER LOC_NKT
	INTEGER LOC_ND
	INTEGER IOS
	INTEGER I
	INTEGER DPTH_INDX
!
	OPEN(UNIT=LU,FILE='NON_THERM_DEGRADATION_SPEC',STATUS='OLD',ACTION='READ')
	  STRING=' '
	  DO WHILE(INDEX(STRING,'Number of energy grid:') .EQ. 0)
	    READ(LU,'(A)')STRING
	  END DO
	  I=INDEX(STRING,':')
	  READ(STRING(I+1:),*)LOC_NKT
	  IF(LOC_NKT .NE. NKT)THEN
	    WRITE(6,*)'Error in RD_NON_THERM'
	    WRITE(6,*)'NKT read is ',LOC_NKT
	    WRITE(6,*)'NKT passed is ',NKT
	    STOP
	  END IF
!
	  STRING=' '
	  DO WHILE(INDEX(STRING,'Number of depth:') .EQ. 0)
	    READ(LU,'(A)')STRING
	  END DO
	  I=INDEX(STRING,':')
	  READ(STRING(I+1:),*)LOC_ND
	  IF(LOC_ND .NE. ND)THEN
	    WRITE(6,*)'Error in RD_NON_THERM'
	    WRITE(6,*)'ND read is ',LOC_ND
	    WRITE(6,*)'ND passed is ',ND
	    STOP
	  END IF
!
! Now get all the data. YE must be at the end of the file, otherwise the
! ordering of the data is irrlevant.
!
	  DPTH_INDX=0
	  DO WHILE(1 .EQ. 1)
	    DO WHILE(STRING .EQ. ' ')
	      READ(LU,'(A)')STRING
	    END DO
	    IF(INDEX(STRING,'Energy grids') .NE. 0)THEN
	      READ(LU,*)(LOC_XKT(I),I=1,NKT)
	      LOC_XKT=ABS(1.0_LDP-LOC_XKT/XKT)
	      IF(MAXVAL(LOC_XKT) .GT. 1.0E-07_LDP)THEN
	        WRITE(6,*)'Error with XKT -- inconsistent energy grid'
	        WRITE(6,*)'Maximum fractional difference is',MAXVAL(LOC_XKT)
	        STOP
	      END IF
	      WRITE(6,*)'Read XKT'; FLUSH(UNIT=6)
	    ELSE IF(INDEX(STRING,'Energy quadrature weights') .NE. 0)THEN
	      READ(LU,*)(dXKT(I),I=1,NKT)
	      WRITE(6,*)'Read dXKT'; FLUSH(UNIT=6)
	    ELSE IF(INDEX(STRING,'Fraction electron heating') .NE. 0)THEN
	      READ(LU,*)(FRAC_ELEC_HEATING(I),I=1,ND)
	      WRITE(6,*)'Read fraction electron heating'; FLUSH(UNIT=6)
	    ELSE IF(INDEX(STRING,'Fraction ion heating') .NE. 0)THEN
	      READ(LU,*)(FRAC_ION_HEATING(I),I=1,ND)
	      WRITE(6,*)'Read fraction ion heating'; FLUSH(UNIT=6)
	    ELSE IF(INDEX(STRING,'Fraction excitation heating') .NE. 0)THEN
	      READ(LU,*)(FRAC_ELEC_HEATING(I),I=1,ND)
	      WRITE(6,*)'Read fraction excitation heating'; FLUSH(UNIT=6)
	    ELSE IF(INDEX(STRING,'Depth:') .NE. 0)THEN
	      DPTH_INDX=DPTH_INDX+1
	      READ(LU,*)(YE(I,DPTH_INDX),I=1,NKT)
	      IF(DPTH_INDX .EQ. ND)EXIT
	    END IF
	    READ(LU,'(A)')STRING
	  END DO
	CLOSE(UNIT=LU)
!
	RETURN
	END
