!
! Subroutine to compute the angular quadrature weights for the
! computation of the mean intensity. The weights are based on
! the technique illustrated by Nordulund.
!
	SUBROUTINE NORDANGQW_V2(QW,R,P,NC,ND,NP,MOMWEIGHT)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Created 14-Jun-2023 - Based on NORDANGQW
!                       Call simplified. Niw use dynmic vector
!                       allocation.
!
	INTEGER NC,ND,NP
	REAL(KIND=LDP) QW(ND,NP),R(ND),P(NP)
!
	REAL(KIND=LDP) MU(NP)
	REAL(KIND=LDP) dMU(NP)
	REAL(KIND=LDP) WEIGHT(NP)
	REAL(KIND=LDP) PSQ(NP)
!
! Local variables
!
	INTEGER I,J,NW
	REAL(KIND=LDP) T1
!
	DO I=1,NP
	  PSQ(I)=P(I)*P(I)
	END DO
!
! Quadrature weights are to be stored in QW(I,J) where I=1,ND
! signifies which radius and J=1,NW signfies which ray.
!
	DO I=1,ND
	  NW=NC+ND-I+1
	  IF(NW .GT. NP)NW=NP
	  T1=R(I)*R(I)
	  DO J=1,NW
	    MU(J)=0.0_LDP
	    IF(R(I) .NE. P(J))MU(J)=SQRT(T1-PSQ(J))/R(I)
	  END DO
	  CALL SET_ACC_dMU(MU,dMU,P,R(I),NW)
	  CALL MOMWEIGHT(MU,dMU,WEIGHT,NW)
	  DO J=1,NW
	    QW(I,J)=WEIGHT(J)
	  END DO
	END DO
!
	RETURN
	END
