	MODULE MOD_DISP
	USE SET_KIND_MODULE
!
! Altered 26-Jul-2024 : Increased number of species to 28.
! Altered 20-Jan-2023 : Added ABS_MEAN
! Altered 15-Dec-2021 : Added ROM_ION_ID for plotting labels (osiris)
! Altered 16-Aug-2019 : Added PLACK_MEAN
! Altered 12-Jun-2015 : NION_MAX increased to 21. Added MI (cur_hmi,19-Aug-2015)
! Altered 13-May-2015 : Inserted ARAD, GAM2, GAM4 and OBSERVED_LEVEL in MODEL_ATOM_DATA.
! Altered 28-Mar-2003 : NUM_IONS variable installed.
! Altered 27-Nov-2000 : MASS changed to AT_MASS
!                       Extra species installed. More consistent with CMFGEN.
!
	INTEGER, PARAMETER :: NSPEC=28
	INTEGER, PARAMETER :: NION_MAX=21
	INTEGER, PARAMETER :: NPHOT_MAX=20
	INTEGER, PARAMETER :: MAX_ION=NSPEC*NION_MAX
	INTEGER  NUM_IONS			!Total number of ions present
!
	REAL(KIND=LDP) AT_MASS(NSPEC)
	REAL(KIND=LDP) AT_NO(NSPEC)
        REAL(KIND=LDP) SOL_MASS_FRAC(NSPEC)       	!Solar mass fraction
        REAL(KIND=LDP) SOL_ABUND_HSCL(NSPEC)      	!Solar abundance with H=12.0
!
	CHARACTER(LEN=10)  SPECIES(NSPEC)
	CHARACTER(LEN=12)  SPECIES_ABR(NSPEC)
	CHARACTER(LEN=12)  ION_ID(MAX_ION)
	CHARACTER(LEN=12)  PLT_ION_ID(MAX_ION)
	CHARACTER(LEN=5)   GEN_ION_ID(NION_MAX)
	CHARACTER(LEN=12)  ROM_ION_ID(NION_MAX)
	INTEGER  SPECIES_LNK(MAX_ION)
	INTEGER SPECIES_BEG_ID(NSPEC)
	INTEGER SPECIES_END_ID(NSPEC)
!
	TYPE MODEL_ATOM_DATA
!
	  REAL(KIND=LDP), ALLOCATABLE :: XzV_F(:,:)		!Level populations in FULL atom (NF.ND)
	  REAL(KIND=LDP), ALLOCATABLE :: XzVLTE_F(:,:)		!LTE level populations in FULL atom (NF.ND)
	  REAL(KIND=LDP), ALLOCATABLE :: LOG_XzVLTE_F(:,:)	!Log LTE level populations in FULL atom (NF.ND)
	  REAL(KIND=LDP), ALLOCATABLE :: W_XzV_F(:,:)           !Level dissolution factors (NF.ND)
	  REAL(KIND=LDP), ALLOCATABLE :: DXzV_F(:)		!Ion population for full atom (ND)
	  REAL(KIND=LDP), ALLOCATABLE :: DXzV(:)		!Ion population SL atom (ND)
	  REAL(KIND=LDP), ALLOCATABLE :: AXzV_F(:,:)		!Oscillator strength (A(I,j), i<j) (NF,NF) Ein. A (i>j)
	  REAL(KIND=LDP), ALLOCATABLE :: EDGEXzV_F(:)           !Ionization energy to g.s. (10^15 Hz) (NF)
	  REAL(KIND=LDP), ALLOCATABLE :: GXzV_F(:) 		!Level statistical weights in full atom (NF)
	  REAL(KIND=LDP), ALLOCATABLE :: ARAD(:)                !Inverse radiative lifetime of level (NF)
	  REAL(KIND=LDP), ALLOCATABLE :: GAM2(:)                !Collisional profile parameter (NF).
	  REAL(KIND=LDP), ALLOCATABLE :: GAM4(:)                !Collisional profile parameter (NF).
	  LOGICAL, ALLOCATABLE :: OBSERVED_LEVEL(:)     !Does level have a know energy (NF).
!
	  REAL(KIND=LDP) ZXzV					!Charge on ion (e.g. 1 for H)
	  REAL(KIND=LDP) GIONXzV_F				!Statistical weight of gs. of ion
!
	  INTEGER, ALLOCATABLE :: F_TO_S_XzV(:)		!Link of full levels to super levels
!
	  INTEGER NXzV_F				!What I am call NF aboev
	  INTEGER N_XzV_PHOT				!Number of photoioization routes
	  INTEGER XzV_ION_LEV_ID(NPHOT_MAX)             !
	  LOGICAL XzV_PRES 				!Indicates whether this ionization stage is present.
!
	  CHARACTER(LEN=40), ALLOCATABLE :: XzVLEVNAME_F(:)
	  CHARACTER(LEN=11) XzV_OSCDATE
	  CHARACTER(LEN=11) NEW_XzV_OSCDATE
!
	END TYPE MODEL_ATOM_DATA
	TYPE (MODEL_ATOM_DATA) ATM(NSPEC*NION_MAX)
C
C Data arrays for full atom.
C
C
	REAL(KIND=LDP), ALLOCATABLE :: R(:)
	REAL(KIND=LDP), ALLOCATABLE :: V(:)
	REAL(KIND=LDP), ALLOCATABLE :: SIGMA(:)
	REAL(KIND=LDP), ALLOCATABLE :: T(:)
	REAL(KIND=LDP), ALLOCATABLE :: ED(:)
	REAL(KIND=LDP), ALLOCATABLE :: CMFGEN_TGREY(:)
	REAL(KIND=LDP), ALLOCATABLE :: dE_RAD_DECAY(:)
!
	REAL(KIND=LDP) ABUND(NSPEC)
!
	REAL(KIND=LDP), ALLOCATABLE :: ROSS_MEAN(:)
	REAL(KIND=LDP), ALLOCATABLE :: FLUX_MEAN(:)
	REAL(KIND=LDP), ALLOCATABLE :: PLANCK_MEAN(:)
	REAL(KIND=LDP), ALLOCATABLE :: ABS_MEAN(:)
!
	REAL(KIND=LDP), ALLOCATABLE :: J_INT(:)
	REAL(KIND=LDP), ALLOCATABLE :: H_INT(:)
	REAL(KIND=LDP), ALLOCATABLE :: K_INT(:)
!
	REAL(KIND=LDP), ALLOCATABLE :: POP_ATOM(:)
	REAL(KIND=LDP), ALLOCATABLE :: MASS_DENSITY(:)
	REAL(KIND=LDP), ALLOCATABLE :: CLUMP_FAC(:)
	REAL(KIND=LDP), ALLOCATABLE :: POPION(:)
!
	REAL(KIND=LDP), ALLOCATABLE :: POPDUM(:,:)
!
! This is the same ID used in CMFGEN, and used to read files etc.
!
	DATA GEN_ION_ID /'0','I','2','III','IV','V',
	1                'SIX','SEV','VIII','IX','X','XI','XII',
	1                'XIII','XIV','XV','XSIX','XSEV','X8','X9','XX'/
!
! This is used for writing plot title etc.
!
	DATA ROM_ION_ID /'\u-1\d','I','II','III','IV','V',
	1                'VI','VII','VIII','IX','X','XI','XII',
	1                'XIII','XIV','XV','XVI','XVII','XVIII','XIX','XX'/
!
	END MODULE MOD_DISP
