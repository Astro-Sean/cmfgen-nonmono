	SUBROUTINE CHECK_LINE_OVERLAP(LINE_ST_INDX,LINE_END_INDX,VEC_TRANS_TYPE,N_LINES,MAX_SIM,NCF)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
	INTEGER N_LINES
	INTEGER MAX_SIM
	INTEGER NCF
	INTEGER LINE_ST_INDX(N_LINES)
	INTEGER LINE_END_INDX(N_LINES)
	CHARACTER(LEN=*) VEC_TRANS_TYPE(N_LINES)
!
	INTEGER MAX_OVER_RES
	INTEGER MAX_OVER_LINES
	INTEGER RES_ZONE(MAX(N_LINES,NCF))
!
	INTEGER ML,I,J,K
	LOGICAL, PARAMETER :: DO_LINE_OVERLAP=.FALSE.
!
	WRITE(6,*)'Entering check_line_overlap.f'; FLUSH(UNIT=6)
!
! This determines the maximum number of overlapping resonance zones at any frequency.
!
	RES_ZONE=0.0_LDP
	DO ML=1,N_LINES
	  IF(VEC_TRANS_TYPE(ML) .EQ. 'BLANK')THEN
	    DO I=LINE_ST_INDX(ML),LINE_END_INDX(ML)
	      RES_ZONE(I)=RES_ZONE(I)+1
	    END DO
	  END IF
	END DO
	MAX_OVER_RES=MAXVAL(RES_ZONE)
	WRITE(6,*)'Maximum number of lines with overlapping resonance zones is',MAX_OVER_RES
!
! Now determine the maximum number of overlapping lines even if some
! frequencies are not in the resonance zone.
!
	IF(DO_LINE_OVERLAP)THEN
	  RES_ZONE=0.0_LDP
	  DO ML=1,N_LINES
	    IF(VEC_TRANS_TYPE(ML) .EQ. 'BLANK')THEN
	      DO J=ML-1,1,-1
	        IF(LINE_END_INDX(J) .GE. LINE_ST_INDX(ML)
	1           .AND. LINE_ST_INDX(J)  .LE. LINE_END_INDX(ML)
	1           .AND. VEC_TRANS_TYPE(J) .EQ. 'BLANK')THEN
	           RES_ZONE(ML)=RES_ZONE(ML)+1
	        END IF
	      END DO
!
	      DO J=ML+1,N_LINES
	        IF(LINE_END_INDX(J) .GE. LINE_ST_INDX(ML)
	1             .AND. LINE_ST_INDX(J)  .LE. LINE_END_INDX(ML)
	1             .AND. VEC_TRANS_TYPE(J) .EQ. 'BLANK')THEN
	            RES_ZONE(ML)=RES_ZONE(ML)+1
	        END IF
	        IF(LINE_ST_INDX(J) .GT. LINE_END_INDX(J))EXIT
	      END DO
	   END IF
	  END DO
	  MAX_OVER_LINES=MAXVAL(RES_ZONE)
	  WRITE(6,*)'Maximum number of lines overlaping any line is:',MAX_OVER_LINES
	END IF
!
!
	IF(MAX_OVER_RES .GT. MAX_SIM)THEN
	  WRITE(6,*)'MAX_SIM in MODEL_SPEC needs to be adjusted to at least: ',MAX_OVER_RES
	  WRITE(6,*)'Current value of MAX_SIM is: ',MAX_SIM
	  STOP
	END IF
	WRITE(6,*)'Exiting check_line_overlap.f'; FLUSH(UNIT=6)
!
	RETURN
	END
