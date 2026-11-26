!
! Returns linear interpolated value for R^2 J. R must be a monotonically increasing or
! decreasing function. Code does no check on whether range is valid.
!
	FUNCTION GET_RSQJ_FROM_J(RVAL,JNU,R,NW)
	USE SET_KIND_MODULE
!
	IMPLICIT NONE
	INTEGER NW
	REAL(KIND=LDP) JNU(NW)  	!
	REAL(KIND=LDP) R(NW)		!Independent variable
	REAL(KIND=LDP) RVAL		!Value at which function is to be interp. to.
	REAL(KIND=LDP) GET_RSQJ_FROM_J
!
	INTEGER J,ILOW,IHIGH,ILOW_SAV
	REAL(KIND=LDP) FRAC
	REAL(KIND=LDP) T1,T2
	SAVE ILOW_SAV
	DATA ILOW_SAV/1/
!
! NW may have changed bewteen this and previous call, hence ILOW_SAV may
! lie outside the valid range (1:NW).
!
	IF(ILOW_SAV .GE. NW)ILOW_SAV=NW-1	
!
! Have two cases to find RVAL for, depending on whether the array R
! increases or decreases with its index.
!
	IF(R(1) .GT. R(NW))THEN
!
! Need to compare RVAL with R array (and not RSAV) in case value is
! same, but array is different.
!
	  IF( (RVAL .GE.  R(ILOW_SAV+1)) .AND.
	1       (RVAL .LE. R(ILOW_SAV)) )THEN
            ILOW=ILOW_SAV
	    IHIGH=ILOW+1
	  ELSE IF(RVAL .GT. R(1))THEN
	    GET_RSQJ_FROM_J=R(1)*R(1)*JNU(1)
	    RETURN
	  ELSE IF(RVAL .LT. R(NW))THEN
	    GET_RSQJ_FROM_J=R(NW)*R(NW)*JNU(NW)
	    RETURN
	  ELSE
!
! LOW and HIGH refer to the size of the INDEX - not R which
! gets smaller as the index decreases.
!
	    ILOW=1
	    IHIGH=NW
	    DO WHILE( (IHIGH-ILOW) .GT. 1)
	      J=(ILOW+IHIGH)/2
	      IF(RVAL .LT. R(J))THEN
	        ILOW=J
	      ELSE
	        IHIGH=J
	      END IF
	    END DO
	  END IF
!
	ELSE
	  IF( (RVAL .LE.  R(ILOW_SAV+1)) .AND.
	1       (RVAL .GE. R(ILOW_SAV)) )THEN
            ILOW=ILOW_SAV
	    IHIGH=ILOW+1
	  ELSE IF(RVAL .LT. R(1))THEN
	    GET_RSQJ_FROM_J=R(1)*R(1)*JNU(1)
	    RETURN
	  ELSE IF(RVAL .GT. R(NW))THEN
	    GET_RSQJ_FROM_J=R(NW)*R(NW)*JNU(NW)
	    RETURN
	  ELSE
!
! Here LOW and HIGH refer to the size of the INDEX - not R which
! gets smaller as the index decreases.
!
	    ILOW=1
	    IHIGH=NW
	    DO WHILE( (IHIGH-ILOW) .GT. 1)
	      J=(ILOW+IHIGH)/2
	      IF(RVAL .LT. R(J))THEN
	        IHIGH=J
	      ELSE
	        ILOW=J
	      END IF
	     END DO
!
	  END IF
	END IF
!
	ILOW_SAV=ILOW
	FRAC=(RVAL-R(IHIGH))/(R(ILOW)-R(IHIGH))
	T1=LOG(R(ILOW)*R(ILOW)*JNU(ILOW))
	T2=LOG(R(IHIGH)*R(IHIGH)*JNU(IHIGH))
	GET_RSQJ_FROM_J= EXP( FRAC*(T1-T2)+T2 )
!
	RETURN
	END
