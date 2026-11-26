!
! Module containing basic data for:
!                   (i) Each species
!          	   (ii) Model atom data (populations etc)
!		  (iii) Basic atmospheric structure.
!
	MODULE MOD_CMFGEN
	USE SET_KIND_MODULE
!
! Altered 25-Jul-2024 : Increased number of species from 26 to 28.
! Altered 20-Jan-2023 : Added ABS_MEAN opacity.
! Altered 17-Aug-2019 : Put on IBIS -- added dWCRXzVd 12-Jul-2019 (cur_cmf_eb)
! Altered 19-Aug-2015 : Added HMI, MAX_IONS_PER_SPECIES increased to 21 (cur_hmi,12-Jun-2015)
! Incoprated 02-Jan-2015: VTURB_VEC, ARAD, GAM2, etc added (for depth depndent profiles).
! Altered 29-Nov-2011 : OLD_LEV_POP_AVAIL added as vector.
! Altered 25-Sep-2011 : LOG_XzVLTE_F, XzVLTE_F_ON_S, LOG_XzVLTE arrays added (26-Nov-2010/5-Apr-2011).
!                         Length of level name changed to 40.
! Altered 19-Jan-2009 : SL_OPTION inserted
! Altered 06-Aug-2008 : PLANCK_MEAN inserted.
! Altered 31-Jan-2002 : FIXEDT etc inserted.
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
! Maximum number of photoionization routes for each species. This only
! dimensions a vector.
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
	REAL(KIND=LDP) SOL_MASS_FRAC(NUM_SPECIES)	!Solar mass fraction
	REAL(KIND=LDP) SOL_ABUND_HSCL(NUM_SPECIES)	!Solar abundance with H=12.0
!
! Total population density of each species (#/cm^3).
!
	REAL(KIND=LDP), ALLOCATABLE :: POP_SPECIES(:,:)
!
! Used to refer to the number of electrons/[species atom] coming from each
! species (also used a temporary store).
!
	REAL(KIND=LDP), ALLOCATABLE :: GAM_SPECIES(:,:)
!
! Conservation equation for each species.
!
	INTEGER EQBAL
	INTEGER EQ_SPECIES(NUM_SPECIES)
	INTEGER FIX_SPECIES(NUM_SPECIES)
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
! Used to indicate whether a level population is available at a previous time step.
! Used for SN calculations (1:NT).
!
	LOGICAL, ALLOCATABLE :: OLD_LEV_POP_AVAIL(:)
!
! 
! Data arrays for full atom. The names in this array should not be changed,
! although new variables can be added.
!
	TYPE MODEL_ATOM_DATA
!
	  REAL(KIND=LDP), ALLOCATABLE :: XzV_F(:,:)		!Level populations in FULL atom
	  REAL(KIND=LDP), ALLOCATABLE :: XzVLTE_F(:,:)		!LTE level populations in FULL atom
          REAL(KIND=LDP), ALLOCATABLE :: LOG_XzVLTE_F(:,:)  	!Log(LTE level populations in FULL atom)
          REAL(KIND=LDP), ALLOCATABLE :: XzVLTE_F_ON_S(:,:) 	!Log(LTE level populations in FULL atom)
	  REAL(KIND=LDP), ALLOCATABLE :: W_XzV_F(:,:)		!Level dissolution factors
	  REAL(KIND=LDP), ALLOCATABLE :: DXzV_F(:)		!Ion population for full atom
	  REAL(KIND=LDP), ALLOCATABLE :: AXzV_F(:,:)		!Oscillator strength (A(I,j), i<j)
	  REAL(KIND=LDP), ALLOCATABLE :: EDGEXzV_F(:)		!Ionization energy to g.s. (10^15 Hz)
	  REAL(KIND=LDP), ALLOCATABLE :: GXzV_F(:)		!Level statistical weights in full atom
	  REAL(KIND=LDP), ALLOCATABLE :: ARAD(:)                !Inverse radiative lifetime of level
	  REAL(KIND=LDP), ALLOCATABLE :: GAM2(:)                !Collisional profile parameter.
	  REAL(KIND=LDP), ALLOCATABLE :: GAM4(:)                !Collisional profile parameter.
	  INTEGER, ALLOCATABLE :: F_TO_S_XzV(:)		!Link of full levels to super levels
	  INTEGER, ALLOCATABLE :: INT_SEQ_XzV(:)
	  CHARACTER(LEN=40), ALLOCATABLE :: XzVLEVNAME_F(:)	!Level name
!
	  REAL(KIND=LDP), ALLOCATABLE :: DXzV(:)		!Ion population for super level
	  REAL(KIND=LDP), ALLOCATABLE :: XzV(:,:)		!Level population in SL atom
	  REAL(KIND=LDP), ALLOCATABLE :: XzVLTE(:,:)		!LTE populations in SL atom
          REAL(KIND=LDP), ALLOCATABLE :: LOG_XzVLTE(:,:)    	!LOG(LTE populations in SL atom)
	  REAL(KIND=LDP), ALLOCATABLE :: dlnXzVLTE_dlnT(:,:)
!
	  REAL(KIND=LDP), ALLOCATABLE :: WSXzV(:,:,:)		!Weights assoc. with s.e. eval.
	  REAL(KIND=LDP), ALLOCATABLE :: WCRXzV(:,:,:)		!Weights assoc. with cooling.
	  REAL(KIND=LDP), ALLOCATABLE :: dWSXzVdT(:,:,:)	!Weights assoc. with d(S.E. Eqn)/dT
	  REAL(KIND=LDP), ALLOCATABLE :: dWCRXzVdT(:,:,:)	!Weights assoc. with d(EHB)/dT
	  REAL(KIND=LDP), ALLOCATABLE :: WSE_X_XzV(:,:)		!X-ray photoionization weights.
	  REAL(KIND=LDP), ALLOCATABLE :: WCR_X_XzV(:,:)		!
!
	  REAL(KIND=LDP), ALLOCATABLE :: APRXzV(:,:)		!Arrays for confirming the SE and RE
	  REAL(KIND=LDP), ALLOCATABLE :: ARRXzV(:,:)		!equations are satisfied.
	  REAL(KIND=LDP), ALLOCATABLE :: BFCRXzV(:,:)		
!
! Data vectors to check cooling rates.
!
	  REAL(KIND=LDP), ALLOCATABLE :: FFXzV(:)
	  REAL(KIND=LDP), ALLOCATABLE :: CRRXzV(:)
	  REAL(KIND=LDP), ALLOCATABLE :: CPRXzV(:)
	  REAL(KIND=LDP), ALLOCATABLE :: COOLXzV(:)
!
	  REAL(KIND=LDP), ALLOCATABLE :: NTCXzV(:)                   !Non-thermal cooling rate for species XzV
	  REAL(KIND=LDP), ALLOCATABLE :: NTIXzV(:)                   !Non-thermal ionization rate (1 electron)
	  REAL(KIND=LDP), ALLOCATABLE :: NTIXzV_2E(:)                !Non-thermal ionization rate (2 electron)
	  REAL(KIND=LDP), ALLOCATABLE :: NT_ION_CXzV(:)              !Non thermal cooling rate (ionization)
	  REAL(KIND=LDP), ALLOCATABLE :: NT_EXC_CXzV(:)              !Non thermal cooling rate (excitation)
	  REAL(KIND=LDP), ALLOCATABLE :: NT_OMEGA(:,:)               !Non thermal collision strength
!
! Data vectors to check charge recombination rates.
!
	  REAL(KIND=LDP), ALLOCATABLE :: CHG_RRXzV(:)
	  REAL(KIND=LDP), ALLOCATABLE :: CHG_PRXzV(:)
!
	  REAL(KIND=LDP) ZXzV			!Charge on ion (=1 for HI)
	  REAL(KIND=LDP) GIONXzV_F		!Statistical weight of ion
!
	  REAL(KIND=LDP) CROSEC_NTFAC           !Factor to scale the non-thermal exictation cross-sections
	  REAL(KIND=LDP) ION_CROSEC_NTFAC       !Factor to scale non-thermal ioinzation cross-sections
!
	  INTEGER NXzV_F		!Number of levels in full atom
	  INTEGER NXzV			!Number of levels in SL atom
	  INTEGER NXzV_IV		!Number of important levels in SL atom
	  INTEGER EQXzV			!Equation in BA matrix for g.s. of atom
	  INTEGER N_XzV_PHOT		!Number of states species can ionize to.
!
	  LOGICAL, ALLOCATABLE :: OBSERVED_LEVEL(:)     !Link of full levels to super levels
!
! Ionization balance equation for each ion. Used to access BAION matrix.
!
	  INTEGER EQXzVBAL
!
! Ion identification for DIELECTRONIC and X-rays. COuld use ID, but left
! in for historical reasons.
!
	  INTEGER INDX_XzV
!
! Indicates whether to hold fixed the first N levels of an ion. If it is
! greater than the number of levels in the super-atom, it is reduced to NS.
!
	  INTEGER FIX_NXzV
!
! Identifications corresponding to each photoionization route.
!
	  INTEGER XzV_ION_LEV_ID(NPHOT_MAX)
!
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
	  CHARACTER*6  XzV_TRANS_TYPE	!Transition type (e.g. Blank, Sobolev etc)
	  CHARACTER*10 XzV_PROF_TYPE	!Type of profile (e.g. Doppler, Stark)
	  CHARACTER*11 XzV_OSCDATE
	  CHARACTER*11 NEW_XzV_OSCDATE
!
	END TYPE MODEL_ATOM_DATA
!
!
	INTEGER EQNE		!Electron conservation equation
!
	REAL(KIND=LDP), ALLOCATABLE :: R(:)		!Radius in units of 10^10 cm
	REAL(KIND=LDP), ALLOCATABLE :: V(:)		!V in units of km/s
	REAL(KIND=LDP), ALLOCATABLE :: SIGMA(:)		!dlnV/dlnR-1
	REAL(KIND=LDP), ALLOCATABLE :: SIGMA_TRUE(:)	!Unclamped dlnV/dlnR-1 (DO_NONMONO_SOB / DO_NONMONO_KNOP only)
	REAL(KIND=LDP), ALLOCATABLE :: T(:)		!Temperature in units of 10^4 K
!
! Next frequency in the sweep, set by CMFGEN_SUB before each
! COMP_J_BLANK call.  Needed at FIRST_FREQ where the usual
! dLOG_NU=LOG(FL_OLD/FL) cannot be evaluated (FL_OLD is undefined);
! the correct one-sided spacing for the bluest edge is LOG(FL/NU_NEXT).
! Only consumed when DO_NONMONO_KNOP=T.
!
	REAL(KIND=LDP) :: NU_NEXT_FREQ=0.0_LDP
	REAL(KIND=LDP), ALLOCATABLE :: ED(:)		!Electron density (#/cm^3)
	REAL(KIND=LDP), ALLOCATABLE :: LANG_COORD(:)    !Lagrangian coordinate
!
	REAL(KIND=LDP), ALLOCATABLE :: ROSS_MEAN(:)	!Rosseland mean opacity.
	REAL(KIND=LDP), ALLOCATABLE :: PLANCK_MEAN(:)	!Planck mean opacity.
	REAL(KIND=LDP), ALLOCATABLE :: ABS_MEAN(:)	!Absorption mean opacity.
	REAL(KIND=LDP), ALLOCATABLE :: FLUX_MEAN(:)	!Flux mean opacity
	REAL(KIND=LDP), ALLOCATABLE :: POP_ATOM(:)	!Total atom density (#/cm^3)
	REAL(KIND=LDP), ALLOCATABLE :: DENSITY(:)	!Mass density (gm/cm^3)
	REAL(KIND=LDP), ALLOCATABLE :: CLUMP_FAC(:)	!Volume filling factor for clumps
	REAL(KIND=LDP), ALLOCATABLE :: POPION(:)	!Ion density
	REAL(KIND=LDP), ALLOCATABLE :: VOL_EXP_FAC(:)	!Volume expansion factor (for time dependent SN models).
	REAL(KIND=LDP), ALLOCATABLE :: VTURB_VEC(:)     !Depth dependent turbulent velocity
!
	REAL(KIND=LDP) STARS_MASS			!In Msun
!
! Variables for determining whether some populations are held fixed
! when the new populations are solved for.
!
	REAL(KIND=LDP) MOD_TAU_SCL_T
	INTEGER MOD_FIX_T_D_ST
	INTEGER MOD_FIX_T_D_END
	LOGICAL MOD_FIXED_NE
	LOGICAL MOD_FIXED_T
	LOGICAL MOD_FIX_IMPURITY
	LOGICAL FIX_IN_BOUND_T
	LOGICAL LIN_INTERP_RD_SN_DATA
	INTEGER FIX_LST_X_DPTHS
!
	LOGICAL, ALLOCATABLE :: IMP_VAR(:)
!
	CHARACTER(LEN=7), PARAMETER :: NAME_CONVENTION='K_FOR_I'
	CHARACTER(LEN=20) FL_OPTION
	CHARACTER(LEN=20) SL_OPTION
	CHARACTER(LEN=20) dE_OPTION
	CHARACTER(LEN=20) IL_OPTION
!
! For unknown reasons, the DEC OSF alpha compiler requires this declaration
! to be at the end of the data module in order not to get alignment problems.
!
	TYPE (MODEL_ATOM_DATA) ATM(NUM_SPECIES*MAX_IONS_PER_SPECIES)
!
! Indicates generic ionization names.
!
	DATA GEN_ION_ID /'0','I','2','III','IV','V',
	1                'SIX','SEV','VIII','IX','X','XI','XII',
	1                'XIII','XIV','XV','XSIX','XSEV','X8','X9','XX'/
!
	END MODULE MOD_CMFGEN
