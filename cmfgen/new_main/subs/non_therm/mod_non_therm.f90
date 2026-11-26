	MODULE MOD_NON_THERM
	USE SET_KIND_MODULE
!
        INTEGER NKT
!
! XKT contains the electroin energy in eV.
! dXKT is the quadrature weight for integrating over energy.
! YE contains the distribution of non-thermal electrons as a function
!   of energy and depth.
!
        REAL(KIND=LDP), ALLOCATABLE :: XKT(:)
        REAL(KIND=LDP), ALLOCATABLE :: dXKT(:)
        REAL(KIND=LDP), ALLOCATABLE :: dXKT_ON_XKT(:)
        REAL(KIND=LDP), ALLOCATABLE :: YE(:,:)
!
! These quantities provide the fraction of the total energy going into the
! three forms of heating --- electron, ionization, and excitation at each depth.
!
        REAL(KIND=LDP), ALLOCATABLE ::  FRAC_ELEC_HEATING(:)
        REAL(KIND=LDP), ALLOCATABLE ::  FRAC_ION_HEATING(:)
        REAL(KIND=LDP), ALLOCATABLE ::  FRAC_EXCITE_HEATING(:)
!
	TYPE NON_THERM_ION_DATA
	  REAL(KIND=LDP) :: N_ION_EL
	  REAL(KIND=LDP) :: ZION
	  REAL(KIND=LDP) :: PQN
	  REAL(KIND=LDP) :: ANG
	  REAL(KIND=LDP) :: ION_POT
	  REAL(KIND=LDP) :: N_ATOM
	  REAL(KIND=LDP) :: A_COL,B_COL,C_COL, D_COL
	  REAL(KIND=LDP), POINTER :: CROSS_SEC(:)
	  REAL(KIND=LDP), ALLOCATABLE :: XTAB(:)
	  REAL(KIND=LDP), ALLOCATABLE :: YTAB(:)
!
          INTEGER :: N_STATES=0
	  INTEGER :: ATOM_STATES(10)
	  INTEGER :: N_ION_ROUTES
	  INTEGER :: ION_LEV(10)
	  INTEGER :: SUM_GION
	  INTEGER :: NTAB
!
	  INTEGER :: LNK_TO_ION
	  INTEGER :: LNK_TO_SPECIES
	  LOGICAL :: PRES
	  LOGICAL :: DO_THIS_ION_ROUTE
	END TYPE NON_THERM_ION_DATA
!
	INTEGER NUM_THD
	INTEGER, PARAMETER :: MAX_NUM_THD=500
!
! These two switches are for testing purposes, and are thus hardwired into
! the code.
!
	LOGICAL, PARAMETER :: INCLUDE_NON_THERM_EXCITATION=.TRUE.
	LOGICAL, PARAMETER :: INCLUDE_NON_THERM_IONIZATION=.TRUE.
	LOGICAL, PARAMETER :: FAST_BETHE_METHOD=.TRUE.
!
	TYPE (NON_THERM_ION_DATA), SAVE :: THD(MAX_NUM_THD)
!
	END MODULE MOD_NON_THERM
