!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
	MODULE MOD_SPACE_GRID_V2_GAM
	USE SET_KIND_MODULE
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Module to define the r-p-z ray or characteristic ray variables.
! This save defining them in CMFGEN and having to pass them
! to subroutines.  Variables are made allocatable so array sizes can
! be passed and array sizes created at runtime.
!
! written  4/8/97  DLM  Previously used f77 and varaibles could not be
!                       allocated.
!
! altered 5/22/97  DLM  Added b_p and b_m for relativistic transfer terms
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Variables for both types of grids
!
	TYPE RAY_DATA
!
	  REAL(KIND=LDP), ALLOCATABLE :: MU(:)
	  REAL(KIND=LDP), ALLOCATABLE :: Z(:)
	  REAL(KIND=LDP), ALLOCATABLE :: R_RAY(:)
!
! Variables for characteristics
!
	  REAL(KIND=LDP), ALLOCATABLE :: S_P(:)                    !Path Length in positive direction
	  REAL(KIND=LDP), ALLOCATABLE :: S_M(:)                    !Path Length in negative direction
	  REAL(KIND=LDP), ALLOCATABLE :: MU_P(:)                   !Angle [cos(theta)] in positive direction
	  REAL(KIND=LDP), ALLOCATABLE :: MU_M(:)                   !Angle [cos(theta)] in negative direction
!
! Variables for frequency independent parts of advection and abberation terms
!
	  REAL(KIND=LDP), ALLOCATABLE :: B_P(:)                    !positive direction
	  REAL(KIND=LDP), ALLOCATABLE :: B_M(:)                    !negative direction
!
	  REAL(KIND=LDP), ALLOCATABLE :: I_P(:)
	  REAL(KIND=LDP), ALLOCATABLE :: I_P_PREV(:)
	  REAL(KIND=LDP), ALLOCATABLE :: I_P_SAVE(:)
	  REAL(KIND=LDP), ALLOCATABLE :: I_M(:)
	  REAL(KIND=LDP), ALLOCATABLE :: I_M_PREV(:)
	  REAL(KIND=LDP), ALLOCATABLE :: I_M_SAVE(:)
!
! These next two arrays are for the gamma ray routine edits
!
	  REAL(KIND=LDP), ALLOCATABLE :: ETA_M(:)
	  REAL(KIND=LDP), ALLOCATABLE :: ETA_P(:)
!
	  INTEGER, ALLOCATABLE :: LNK(:)		!Indicates link to original grid
!
	  REAL(KIND=LDP), ALLOCATABLE :: I_IN_BND_STORE(:)
	  REAL(KIND=LDP) P_RAY				      	!P for ray
	  REAL(KIND=LDP) FREQ_CONV_FAC				!Frequency conversion factor for hollow core.

	  INTEGER  NZ	 	                    	!number of grid points along a p-ray
!
!	  INTEGER, ALLOCATABLE :: NINS(:)		!Index array to keep track of # of inserted
!							!points along a ray
!	  REAL(KIND=LDP), ALLOCATABLE :: Z_SM(:)
!	  REAL(KIND=LDP), ALLOCATABLE :: Z_LG(:)
!
	END TYPE RAY_DATA
!
	TYPE (RAY_DATA) RAY(500)
!
! Arrays for "optical depth"
!
	REAL(KIND=LDP), ALLOCATABLE :: TAU(:)
	REAL(KIND=LDP), ALLOCATABLE :: DTAU(:)
!
! Defined along a ray and on the CMFGEN grid.
!
	REAL(KIND=LDP), ALLOCATABLE :: I_P_GRID(:)
	REAL(KIND=LDP), ALLOCATABLE :: I_M_GRID(:)
!
! Quadrature weights for characteristics. These are defined on the
! CMFGEN grid.
!
	REAL(KIND=LDP), ALLOCATABLE :: JQW_P(:,:)                  !Positive mu
	REAL(KIND=LDP), ALLOCATABLE :: HQW_P(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: KQW_P(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: NQW_P(:,:)
!
	REAL(KIND=LDP), ALLOCATABLE :: JQW_M(:,:)                  !NEGATIVE MU
	REAL(KIND=LDP), ALLOCATABLE :: HQW_M(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: KQW_M(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: NQW_M(:,:)
!
	REAL(KIND=LDP), ALLOCATABLE :: R_EXT_SAV(:)
	REAL(KIND=LDP), ALLOCATABLE :: FREQ_STORE(:)
!
	INTEGER CUR_LOC
	INTEGER N_STORE
	INTEGER ND_SAV
	INTEGER ND_EXT_SAV
	INTEGER NP_SAV
	LOGICAL, SAVE :: RAY_POINTS_INSERTED=.FALSE.
!
      END MODULE MOD_SPACE_GRID_V2_GAM
