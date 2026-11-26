!
! Module containing basic data for:
!                   (i) Each species
!          	   (ii) Model atom data (populations etc)
!		  (iii) Basic atmospheric structure.
!
      MODULE STEQ_DATA_MOD
      USE SET_KIND_MODULE
!
!Alteerd 26-Jul-2024: Increased the number of species to 28.
!Altered 12-Jul-2019: Added EHB matrices (added to IBIS 17-Aug-2019).
!Altered 19-Aug-2015: BA_MAX_IONS_PER_SPECIES increaed to 21 (cur_hmi, 12-Jun-2015)
!
! Number of atomic species (e.g. H, C, N is 3 species).
!
	INTEGER, PARAMETER :: BA_NUM_SPECIES=28
!
! Maximum number of ionization stages per species. For H, need at this number
! has to be 2 or higher (as I and II). A setting of 10 implies that we can treat
! full atoms for ION_IX.
!
	INTEGER, PARAMETER :: BA_MAX_IONS_PER_SPECIES=21
	INTEGER, PARAMETER :: BA_MAX_NUM_IONS=BA_NUM_SPECIES*BA_MAX_IONS_PER_SPECIES
!
! Maximum number of photoionization routes for each species.
!
	INTEGER, PARAMETER :: BA_NPHOT_MAX=10
!
! Storage locations for the Statistical euilibrium equations, and its
! variation.
!
	TYPE STAT_EQ_DATA
          REAL(KIND=LDP), POINTER :: STEQ(:,:)  		!Statistical equilibrium Eqns. for XzV
          REAL(KIND=LDP), POINTER :: STEQ_ADV(:)  		!Statistical equilibrium ion eqns. for advection terms
          REAL(KIND=LDP), POINTER :: QFV_P(:,:)  		!
          REAL(KIND=LDP), POINTER :: QFV_R(:,:)  		!
          REAL(KIND=LDP), POINTER :: QFV_P_EHB(:,:)  		!
          REAL(KIND=LDP), POINTER :: QFV_R_EHB(:,:)  		!
          REAL(KIND=LDP), POINTER :: BA(:,:,:,:)		!BA matrix for XzV
          REAL(KIND=LDP), POINTER :: BA_PAR(:,:,:)		!BA matrix for XzV (diagonal terms)
          REAL(KIND=LDP), POINTER :: T_EHB(:) 	 		!
	  INTEGER, POINTER :: LNK_TO_IV(:)    	!
	  INTEGER, POINTER :: LNK_TO_F(:)     	!
	  INTEGER, POINTER :: EQ_IN_BA(:)     	!
	  INTEGER, POINTER :: EQ_TO_ION_LEV_PNT(:)	!
	  INTEGER, POINTER :: ION_LEV_TO_EQ_PNT(:)	!
	  INTEGER, POINTER :: STRT_ADV_ID(:)
	  INTEGER, POINTER :: END_ADV_ID(:)
	  INTEGER N_SE                        	!Number of S.E. Eqns. for XzV
	  INTEGER N_IV                        	!Number of important variables for XzV
	  INTEGER NUMBER_BAL_EQ			!Conservation equation
	  INTEGER XRAY_EQ			!Equation for Auger ionizations
	  LOGICAL Xzv_PRES                     	!Indicates whether ion is present.
          LOGICAL IMPURITY_SPECIES
	END TYPE STAT_EQ_DATA
!
        REAL(KIND=LDP), ALLOCATABLE :: STEQ_T(:)
        REAL(KIND=LDP), ALLOCATABLE :: STEQ_ED(:)
        REAL(KIND=LDP), ALLOCATABLE :: BA_T(:,:,:)
        REAL(KIND=LDP), ALLOCATABLE :: BA_ED(:,:,:)
        REAL(KIND=LDP), ALLOCATABLE :: BA_ADV_TERM(:,:)
        REAL(KIND=LDP), ALLOCATABLE :: BA_T_PAR(:,:)
!
        REAL(KIND=LDP), ALLOCATABLE :: STEQ_T_EHB(:)
        REAL(KIND=LDP), ALLOCATABLE :: BA_T_EHB(:,:,:)
        REAL(KIND=LDP), ALLOCATABLE :: BA_T_PAR_EHB(:,:)
!
        TYPE (STAT_EQ_DATA)  SE(BA_NUM_SPECIES*BA_MAX_IONS_PER_SPECIES)
!
      END MODULE STEQ_DATA_MOD
