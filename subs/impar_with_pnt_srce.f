!
! Subroutine to compute the impact parameter values .
!
	SUBROUTINE IMPAR_WITH_PNT_SRCE(P,R,RP,NC,ND,NP,R_PNT_SRCE,NC_PNT_SRCE)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Altered 06-Jun-2019 - Renamed to IMPAR_WITH_PNT_SRCE.
!                         Added R_PNT_SRCE,NC_PNT_SRCE to call.
!                         Now handles rays for point source.
! Altered 05-Dec-1996 - END DO used to terminate DO loops.
! Altered 24-May-1996 - IONE inserted
! Altered 17-Feb-1986 - Core rays distributed equally in mu rather than p.
!
	INTEGER NC,ND,NP,I
	INTEGER NC_PNT_SRCE
	REAL(KIND=LDP) P(NP),R(ND)
	REAL(KIND=LDP) RP,R_PNT_SRCE
!
	REAL(KIND=LDP) DELMU,dP,T1
	REAL(KIND=LDP), PARAMETER :: RZERO=0.0
	REAL(KIND=LDP), PARAMETER :: RONE=1.0
!
	IF(NC_PNT_SRCE .EQ. 0)THEN
	  P(1)=RZERO
	  DELMU=RONE/NC
	  DO I=2,NC
	    P(I)=RP*SQRT(RONE-(DELMU*(NC-I+1))**2)
	  END DO
	ELSE IF(R_PNT_SRCE .GE. RP)THEN
	  WRITE(6,*)'Error in IMPAR_WITH_PNT_SRCE'
	  WRITE(6,*)'Point sorce must have its radius less than shell'
	  WRITE(6,*)'R_PNT_SRCE=',R_PNT_SRCE
	  WRITE(6,*)'RP=',RP
	  STOP
	ELSE
	  P(1)=RZERO
	  dP=R_PNT_SRCE/(NC_PNT_SRCE-1)
	  DO I=2,NC_PNT_SRCE
	    P(I)=(I-1)*dP
	  END DO
	  P(NC_PNT_SRCE+1)=P(NC_PNT_SRCE)+1.0E-03_LDP*dP
!
	  T1=SQRT(RONE-(P(NC_PNT_SRCE+1)/RP)**2)
	  DELMU=T1/(NC-2)
	  DO I=NC_PNT_SRCE+2,NC
	    P(I)=RP*SQRT(RONE-(DELMU*(NC-I+1))**2)
	  END DO
	END IF	
!
	DO I=NC+1,NP
	  P(I)=R(ND-I+NC+1)
	END DO
!
	DO I=1,NC
	  WRITE(188,'(I4,3ES18.8)')I,P(I),P(I)/RP,SQRT(RONE-(P(I)/RP)**2)
	END DO
	DO I=NC+1,NP
	  WRITE(188,'(I4,2ES18.8)')I,P(I),P(I)/RP
	END DO
	FLUSH(UNIT=188)
!
	RETURN
	END
