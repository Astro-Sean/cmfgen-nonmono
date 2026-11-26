C
C Subroutine to compute the angular quadrature weights for the
C computation of the angular moments of the radiation field
C assuming spherical geometry.
C
	SUBROUTINE GENANGQW_V2(QW,R,P,NC,ND,NP,MOMWEIGHT,AT_HALF)
	USE SET_KIND_MODULE
	IMPLICIT NONE
C
C 23-Dec-2004 : Only print out warning message when AT_HALF is FALSE.
C                 When AT_HALF is true, we always do extrapolation.
C 04-Sep-2002 : Check before doing SQRT installed (for INTEL compiler).
C 26-May-1996 : Call to DP_ZERO removed
C               ERROR_LU installed.
C               MOMWEIGHT declared external.
C 18-Jul-1991 - Commented out IMPLICIT NONE To overcome CRAY compiler bug.
C               Alternatively it could have been fixed by declaring
C               MOMWEIGHT to be EXTERNAL, and INTEGER.
C 11-Jan-1989 - Renamed to GENANGQW - Option AT_HALF installed.
C 25-Nov-1986 - Created
C
	LOGICAL AT_HALF
	INTEGER NC,ND,NP
	REAL(KIND=LDP) QW(ND,NP),R(ND),P(NP)
C
	REAL(KIND=LDP) PSQ(NP),MU(NP),dMU(NP),WEIGHT(NP)
	INTEGER ERROR_LU,LUER
	EXTERNAL MOMWEIGHT,ERROR_LU
C
	INTEGER I,J,NW,MAXND,COUNTER
	REAL(KIND=LDP) T1,T2
C
C Initialization section
C
	QW(:,:)=0.0_LDP
	DO I=1,NP
	  PSQ(I)=P(I)*P(I)
	END DO
	COUNTER=0
C
C***********************************************************************
C If AT_HALF is TRUE the angular quadrature weight are computed at the
C midpoints of the mesh.
C
C Quadrature weights are to be stored in QW(I,J) where I=1,ND-1
C signifies the radius midpoint (=I+0.5) and J=1,NW signifies which ray.
C This section of the routine is primarily used to compute the H and N
C quadrature weights for the routines where H and N are discretized
C at the mid points of the mesh.
C***********************************************************************
C***********************************************************************
C If AT_HALF is FALSE the angular quadrature weights are computed on
C the mesh.
C
C Quadrature weights are to be stored in QW(I,J) where I=1,ND
C signifies which radius and J=1,NW signifies which ray.
C
C***********************************************************************
C
	IF(AT_HALF)THEN
	  MAXND=ND-1
	ELSE
	  MAXND=ND
	END IF
C
	DO I=1,MAXND
	  NW=NC+MAXND-I+1
	  IF(NW .GT. NP)NW=NP
	  IF(AT_HALF)THEN
	    T2=0.5_LDP*( R(I)+R(I+1) )
	    T1=T2*T2
	  ELSE
	    T2=R(I)
	    T1=T2*T2
	  END IF
	  DO J=1,NW
	    MU(J)=0.0_LDP
	    IF(R(I) .NE. P(J))MU(J)=SQRT(T1-PSQ(J))/T2
	  END DO
	  CALL SET_ACC_dMU(MU,dMU,P,T2,NW)
	  CALL MOMWEIGHT(MU,dMU,WEIGHT,NW)
	  IF(MU(NW) .NE. 0.0_LDP)COUNTER=COUNTER+1
	  DO J=1,NW
	    QW(I,J)=WEIGHT(J)
	  END DO
	END DO
C
	IF(COUNTER .NE. 0 .AND. .NOT. AT_HALF)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,500)COUNTER
500	  FORMAT(1X,'Warning -',I3,' extrapolations to zero in GENANGQW')
	END IF
C
	RETURN
	END
