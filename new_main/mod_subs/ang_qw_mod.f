!
! Data module defining Quadrature Weights, and related variables,
! for CMFGEN.
!
	MODULE ANG_QW_MOD
	USE SET_KIND_MODULE
!
! Created 02-May-2004
!
! Angular qudrature weights for J (AQW), H, and K. These are
! frequency indepoendent.
!
	REAL(KIND=LDP), ALLOCATABLE :: P(:)              !NP
	REAL(KIND=LDP), ALLOCATABLE :: AQW(:,:)          !ND,NP
	REAL(KIND=LDP), ALLOCATABLE :: HQW(:,:)          !ND,NP
	REAL(KIND=LDP), ALLOCATABLE :: KQW(:,:)          !ND,NP
	REAL(KIND=LDP), ALLOCATABLE :: NQW(:,:)          !ND,NP
!
! Defined at the mid points of the mesh.
!
	REAL(KIND=LDP), ALLOCATABLE :: HMIDQW(:,:)       !ND,NP
	REAL(KIND=LDP), ALLOCATABLE :: NMIDQW(:,:)       !ND,NP
!
! If required, these arrays shouLd have size NDEXT*NPEXT
!
	REAL(KIND=LDP), ALLOCATABLE :: PEXT(:)            !NPEXT
	REAL(KIND=LDP), ALLOCATABLE :: AQWEXT(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: HQWEXT(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: KQWEXT(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: NQWEXT(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: HMIDQWEXT(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: NMIDQWEXT(:,:)
!
! Parameters, vectors, and arrays for computing the observed flux.
! Was 2000 (changed 20-Apr-2009).
!
        INTEGER, PARAMETER :: NST_CMF=6000
        REAL(KIND=LDP)  NU_STORE(NST_CMF)
!
        INTEGER NP_OBS_MAX
        INTEGER NP_OBS
        REAL(KIND=LDP) V_AT_RMAX                !Used if we extend the atmosphere.
        REAL(KIND=LDP) RMAX_OBS
        REAL(KIND=LDP) HFLUX_AT_OB		!In comoving frame
        REAL(KIND=LDP) HFLUX_AT_IB
!
! We allocate memory for the following vectors as we use them for the regular
! flux computation, and when extra depth points are inserted (ACCURATE=.TRUE.)
!
        REAL(KIND=LDP), ALLOCATABLE :: IPLUS_STORE(:,:)
        REAL(KIND=LDP), ALLOCATABLE :: P_OBS(:)
        REAL(KIND=LDP), ALLOCATABLE :: IPLUS(:)
        REAL(KIND=LDP), ALLOCATABLE :: MU_AT_RMAX(:)
        REAL(KIND=LDP), ALLOCATABLE :: HQW_AT_RMAX(:)
!
	END MODULE ANG_QW_MOD
