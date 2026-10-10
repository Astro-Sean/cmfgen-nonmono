!
! Routine to compute the Average energy of a super-level. The temperature is
! not taken into account.
!
	SUBROUTINE AVE_LEVEL_ENERGY(AVE_ENERGY,EDGE,STAT_WT,F_TO_S,
	1                  EQHYD,N_S,N_F,NT,HYD_PRES)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Input
	LOGICAL HYD_PRES
	INTEGER EQHYD,NT
	INTEGER N_S,N_F
	REAL(KIND=LDP) EDGE(N_F)
	REAL(KIND=LDP) STAT_WT(N_F)
	INTEGER F_TO_S(N_F)
!
! Output

	REAL(KIND=LDP) AVE_ENERGY(NT)
!
! Local variables
!
	REAL(KIND=LDP) G_SUM(N_S)
	REAL(KIND=LDP) EDGE_SUM(N_S)
	INTEGER I,J
!
	IF(HYD_PRES)THEN
	  G_SUM(:)=0.0_LDP
	  EDGE_SUM(:)=0.0_LDP
	  DO J=1,N_F
	    I=F_TO_S(J)
            G_SUM(I)=G_SUM(I)+STAT_WT(J)
	    EDGE_SUM(I)=EDGE_SUM(I)+EDGE(J)*STAT_WT(J)
	  END DO
	  DO I=1,N_S
	    AVE_ENERGY(EQHYD+I-1)=EDGE_SUM(I)/G_SUM(I)
	  END DO
	END IF
!
	RETURN
	END
