! Module to store the created mu grids into a structure for every depth point
!
! Created July 25, 2014
!-----------------------------------------------------------------------------
	MODULE GAM_MU_MOD
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
	TYPE MU_GRID
	  INTEGER :: MU_PTS
	  INTEGER :: MON_MU_PTS
	  REAL(KIND=LDP), DIMENSION(:), ALLOCATABLE :: MU_VECTOR
	  REAL(KIND=LDP), DIMENSION(:), ALLOCATABLE :: MON_MU
	  REAL(KIND=LDP) :: dMU
	END TYPE
!
	TYPE (MU_GRID), ALLOCATABLE :: R_MU(:)
	END MODULE
