!
! Module containing basic data for:
!                   (i) Each species
!          	   (ii) Model atom data (populations etc)
!		  (iii) Basic atmospheric structure.
!
	MODULE MOD_CMF_OBS
	USE SET_KIND_MODULE
!
! Altered: 26-Jul-2024 : Increased number of species to 28.
! Altered: 20-Aug-2019 : Added PLANCK_MEAN
! Altered: 17-Jan-2017 : Added CMFGEN_TGREY and dE_RAD_DECAY.
! Altered: 19-Aug-2015 : Added MI, Max number of ions increased to 21
! Altered: 18-May-2015 : Changed GAM2, GAM4 to C4 and C6 (quadratic and Van der Waals
!                           interacton constants).
! Altered 17-Dec-2011:  Added LOG_XzVLTE_F, LOG_XzVLTE, and XzVLTE_F_ON_S to MODEL_ATOM_DATA
!
! Number of atomic species (e.g. H, C, N is 3 species).
!
	INTEGER, PARAMETER :: NUM_SPECIES=28
!
! Maximum number of ionization stages per species. For H, need at this number
! has to be 2 or higher (as I and II). A setting of 10 implies that we can treat
! full atoms for ION_IX.
!
	INTEGER, PARAMETER :: MAX_IONS_PER_SPECIES=21
	INTEGER, PARAMETER :: MAX_NUM_IONS=NUM_SPECIES*MAX_IONS_PER_SPECIES
!
! Maximum number of photoionization routes for each species.
!
	INTEGER, PARAMETER :: NPHOT_MAX=10
!
! Actual number of ions in calculation. Stored sequentially.
!
	INTEGER NUM_IONS
!
	REAL(KIND=LDP) AT_MASS(NUM_SPECIES)		!Atomic mass in amu
	REAL(KIND=LDP) AT_NO(NUM_SPECIES)		!Atomic number of species
	REAL(KIND=LDP) AT_ABUND(NUM_SPECIES)		!Fractional species abundance
	REAL(KIND=LDP) ABUND_SCALE_FAC(NUM_SPECIES)	!To scale species abundance
	REAL(KIND=LDP) SOL_MASS_FRAC(NUM_SPECIES)	!Solar mass fraction
	REAL(KIND=LDP) SOL_ABUND_HSCL(NUM_SPECIES)	!Solar abundance with H=12.0
!
! Total population density of each species (#/cm^3).
!
	REAL(KIND=LDP), ALLOCATABLE :: POP_SPECIES(:,:)
!
! Conservation equation for each species.
!
	INTEGER EQ_SPECIES(NUM_SPECIES)
!
! Indicate at what location in the ion arrays each species starts
! an ends. This includes the highest ionization stage (e.g. H+).
!
	INTEGER SPECIES_BEG_ID(NUM_SPECIES)
	INTEGER SPECIES_END_ID(NUM_SPECIES)
!
! Principal species abbreviation (e.g. CARB for carbon)
!
	CHARACTER*10 SPECIES(NUM_SPECIES)
!
! Abbreviation for species used to identify ions (e.g. Ca for calcium).
!
	CHARACTER*2  SPECIES_ABR(NUM_SPECIES)
!
	CHARACTER*5  GEN_ION_ID(MAX_IONS_PER_SPECIES)
!
! Link from ion identification to parent species.
!
	INTEGER SPECIES_LNK(MAX_NUM_IONS)
!
! Identification of ion.
!
	CHARACTER*12 ION_ID(MAX_NUM_IONS)
!
! Indicates whether a species (e.g. Carbon) is included.
!
	LOGICAL SPECIES_PRES(NUM_SPECIES)
!
! 
! Data arrays for full atom. These are ordered according to variable type.
! This helps to enure that they fall on correct boudaries etc.
!
	TYPE MODEL_ATOM_DATA
!
	  REAL(KIND=LDP), ALLOCATABLE :: XzV_F(:,:)		!Level populations in FULL atom
	  REAL(KIND=LDP), ALLOCATABLE :: XzVLTE_F(:,:)		!LTE level populations in FULL atom
	  REAL(KIND=LDP), ALLOCATABLE :: LOG_XzVLTE_F(:,:)      !Log(LTE level populations in FULL atom)
	  REAL(KIND=LDP), ALLOCATABLE :: XzVLTE_F_ON_S(:,:)     !Log(LTE level populations in FULL atom
	  REAL(KIND=LDP), ALLOCATABLE :: W_XzV_F(:,:)		!Level dissolution factors
	  REAL(KIND=LDP), ALLOCATABLE :: DXzV_F(:)		!Ion population for full atom
	  REAL(KIND=LDP), ALLOCATABLE :: AXzV_F(:,:)		!Oscillator strength (A(I,j), i<j)
	  REAL(KIND=LDP), ALLOCATABLE :: EDGEXzV_F(:)		!Ionization energy to g.s. (10^15 Hz)
	  REAL(KIND=LDP), ALLOCATABLE :: GXzV_F(:)		!Level statistical weights in full atom
	  REAL(KIND=LDP), ALLOCATABLE :: ARAD(:)		!Inverse radiative lifetime of level
	  REAL(KIND=LDP), ALLOCATABLE :: C4(:)  		!Collisional profile parameter.
	  REAL(KIND=LDP), ALLOCATABLE :: C6(:)	 		!Collisional profile parameter.
!
	  REAL(KIND=LDP), ALLOCATABLE :: DXzV(:)		!Ion population for super level
	  REAL(KIND=LDP), ALLOCATABLE :: XzV(:,:)		!Level population in SL atom
	  REAL(KIND=LDP), ALLOCATABLE :: XzVLTE(:,:)		!LTE populations in SL atom
	  REAL(KIND=LDP), ALLOCATABLE :: LOG_XzVLTE(:,:)        !LOG(LTE populations in SL at
	  REAL(KIND=LDP), ALLOCATABLE :: dlnXzVLTE_dlnT(:,:)
!
	  REAL(KIND=LDP) ZXzV			!Charge on ion (=1 for HI)
	  REAL(KIND=LDP) GIONXzV_F		!Statistical weight of ion
!
! Identifications corresponding to each photoionization route.
!
	  INTEGER XzV_ION_LEV_ID(NPHOT_MAX)
	  INTEGER, ALLOCATABLE :: F_TO_S_XzV(:)	!Link of full levels to super levels
	  INTEGER, ALLOCATABLE :: INT_SEQ_XzV(:)
!
	  INTEGER NXzV_F		!Number of levels in full atom
	  INTEGER NXzV		!Number of levels in SL atom
	  INTEGER EQXzV		!Equation in BA matrix for g.s. of atom
	  INTEGER N_XzV_PHOT		!Number of states species can ionize to.
!
	  LOGICAL, ALLOCATABLE :: OBSERVED_LEVEL(:)	!Link of full levels to super levels
	
! Indicates whether a species is present. The final ionization state is regarded
! as not present (e.g. HII_PRES is ALWAYS false, even though we treat H+ when
! we treat HI).
!
	  LOGICAL XzV_PRES		!indicates
!
! Dielectronic variables.
!
	  LOGICAL DIE_AUTO_XzV
	  LOGICAL DIE_WI_XzV
!
	  CHARACTER(LEN=40), ALLOCATABLE :: XzVLEVNAME_F(:)	!Level name
	  CHARACTER(LEN=6)  XzV_TRANS_TYPE	!Transition type (e.g. Blank, Sobolev etc)
	  CHARACTER(LEN=10) XzV_PROF_TYPE	!Type of profile (e.g. Doppler, Stark)
	  CHARACTER(LEN=11) XzV_OSCDATE
	  CHARACTER(LEN=11) NEW_XzV_OSCDATE
!
	END TYPE MODEL_ATOM_DATA
!
!
	INTEGER EQNE		!Electron conservation equation
!
	REAL(KIND=LDP), ALLOCATABLE :: R(:)		!Radius in units of 10^10 cm
	REAL(KIND=LDP), ALLOCATABLE :: V(:)		!V in units of km/s
	REAL(KIND=LDP), ALLOCATABLE :: SIGMA(:)		!dlnV/dlnR-1 (clamped >= -0.999)
	REAL(KIND=LDP), ALLOCATABLE :: SIGMA_TRUE(:)	!Unclamped dlnV/dlnR-1 for DO_NONMONO_SOB EW estimators
	REAL(KIND=LDP), ALLOCATABLE :: T(:)		!Temperature in units of 10^4 K
	REAL(KIND=LDP), ALLOCATABLE :: ED(:)		!Electron density (#/cm^3)
	REAL(KIND=LDP), ALLOCATABLE :: LANG_COORD(:)	!Langrangian coordinate
!
	REAL(KIND=LDP), ALLOCATABLE :: ROSS_MEAN(:)	!Rosseland mean opacity.
	REAL(KIND=LDP), ALLOCATABLE :: FLUX_MEAN(:)	!Flux mean opacity
	REAL(KIND=LDP), ALLOCATABLE :: PLANCK_MEAN(:)	!Flux mean opacity
	REAL(KIND=LDP), ALLOCATABLE :: POP_ATOM(:)	!Total atom density (#/cm^3)
	REAL(KIND=LDP), ALLOCATABLE :: DENSITY(:)	!Mass density (gm/cm^3)
	REAL(KIND=LDP), ALLOCATABLE :: CLUMP_FAC(:)	!Volume filling factor for clumps
	REAL(KIND=LDP), ALLOCATABLE :: POPION(:)	!Ion density
	REAL(KIND=LDP), ALLOCATABLE :: CMFGEN_TGREY(:)	!Computed by CMFGEN (named consistently with DISPGEN)
	REAL(KIND=LDP), ALLOCATABLE :: dE_RAD_DECAY(:)	!Computed by CMFGEN
!
	REAL(KIND=LDP) STARS_MASS			!In Msun
	REAL(KIND=LDP) STARS_LUM			!In Lsun
	REAL(KIND=LDP) PRESSURE_VTURB			!In km/s
!
	TYPE (MODEL_ATOM_DATA) ATM(NUM_SPECIES*MAX_IONS_PER_SPECIES)
!
! Indicates generic ionization names.
!
	DATA GEN_ION_ID /'0','I','2','III','IV','V',
	1                'SIX','SEV','VIII','IX','X','XI','XII',
	1                'XIII','XIV','XV','XSIX','XSEV','X8','X9','XX'/
!
	END MODULE MOD_CMF_OBS
