! This SUBROUTINE will calculate the prefactor of the KLEIN-NISHINA scattering
! redistribution function of only the frequency dependence
!
! Created June 16, 2014
!
!-------------------------------------------------------------------------------------
!
	SUBROUTINE KLEIN_SCAT(KLEIN_ARRAY,NF_GRID_PTS,E_VEC)
	USE SET_KIND_MODULE
!
	IMPLICIT NONE
!
	INTEGER :: NF_GRID_PTS
	REAL(KIND=LDP) :: KLEIN_ARRAY(NF_GRID_PTS,NF_GRID_PTS) ! KLEIN_ARRAY(INCOMING,OUTGOING)
	REAL(KIND=LDP) :: E_VEC(NF_GRID_PTS)
	INTEGER :: IN_FI,FI  ! FI=frequency index, IN denotes incoming
        REAL(KIND=LDP), PARAMETER :: a=8.09330118E-21_LDP ! a=h/mc^2 in units of seconds
	REAL(KIND=LDP), PARAMETER :: ONE=1.0_LDP
!
	KLEIN_ARRAY=0.0_LDP
	DO FI=1,NF_GRID_PTS
	  DO IN_FI=1,FI
	    IF (E_VEC(FI) .LE. E_VEC(IN_FI)) THEN
	      KLEIN_ARRAY(IN_FI,FI)=(ONE/(E_VEC(IN_FI)*(E_VEC(IN_FI)/a)))*(E_VEC(IN_FI)/E_VEC(FI) &
			+E_VEC(FI)/E_VEC(IN_FI)+2.0_LDP*(ONE/E_VEC(IN_FI)-ONE/E_VEC(FI)) &
			+(ONE/E_VEC(IN_FI)-ONE/E_VEC(FI))*(ONE/E_VEC(IN_FI)-ONE/E_VEC(FI)))
	    END IF
	  END DO
	END DO
!
	RETURN
	END SUBROUTINE KLEIN_SCAT
