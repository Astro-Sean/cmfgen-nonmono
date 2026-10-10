!
! Module containing DATA for charge exchange reactions.
!
	MODULE CHG_EXCH_MOD_V3
	USE SET_KIND_MODULE
	  INTEGER N_CHG_RD
	  INTEGER N_CHG
	  INTEGER LUER
!
! Altered 10-Sep-2000 : Bug fix. THI_CHG was not beeing set is not read in.
! Altered 01-Oct-1998 : Fitting range and FORMAT date is now read in.
!                       An alternative level name can also be specified
!                       (via {}).
! Altered 20-Aug-1998
! Created 24-Jun-1998
!
! Reaction data
!
          INTEGER, ALLOCATABLE :: TYPE_CHG_RD(:)
          REAL(KIND=LDP), ALLOCATABLE :: COEF_CHG_RD(:,:)
          REAL(KIND=LDP), ALLOCATABLE :: TLO_CHG_RD(:)
          REAL(KIND=LDP), ALLOCATABLE :: THI_CHG_RD(:)
          LOGICAL, ALLOCATABLE :: CHG_INCLUDED_RD(:)
	  CHARACTER(LEN=8),  ALLOCATABLE :: SPEC_ID_CHG_RD(:,:)
	  CHARACTER(LEN=30), ALLOCATABLE :: LEV_NAME_CHG_RD(:,:)
	  CHARACTER(LEN=30), ALLOCATABLE :: ALT_LEV_NAME_CHG_RD(:,:)
!
          INTEGER, ALLOCATABLE :: TYPE_CHG(:)
          REAL(KIND=LDP), ALLOCATABLE :: COEF_CHG(:,:)
          REAL(KIND=LDP), ALLOCATABLE :: G_CHG(:,:)
          REAL(KIND=LDP), ALLOCATABLE :: TLO_CHG(:)
          REAL(KIND=LDP), ALLOCATABLE :: THI_CHG(:)
	  CHARACTER(LEN=8),  ALLOCATABLE :: SPEC_ID_CHG(:,:)
	  CHARACTER(LEN=30), ALLOCATABLE :: LEV_NAME_CHG(:,:)
	  CHARACTER(LEN=0),  ALLOCATABLE :: ALT_LEV_NAME_CHG(:,:)
!
! Arrays required to update STEQ (statistical equilibrium) and BA (variation
! of STEQ) arrays.
!
          REAL(KIND=LDP), ALLOCATABLE :: Z_CHG(:,:)
          REAL(KIND=LDP), ALLOCATABLE :: AI_AR_CHG(:,:)
          REAL(KIND=LDP), ALLOCATABLE :: dlnAI_AR_CHG_dlnT(:,:)
          REAL(KIND=LDP), ALLOCATABLE :: COOL_CHG(:,:)
!
	  INTEGER, ALLOCATABLE ::  ID_ION_CHG(:,:)
	  INTEGER, ALLOCATABLE ::  LEV_IN_POPS_CHG(:,:)
	  INTEGER, ALLOCATABLE ::  LEV_IN_ION_CHG(:,:)
	  LOGICAL, ALLOCATABLE ::  CHG_REACTION_AVAILABLE(:)
!
	  INTEGER, PARAMETER :: N_COEF_MAX=5
!
	  LOGICAL INITIALIZE_ARRAYS
	  LOGICAL DO_CHG_EXCH
!
	END MODULE CHG_EXCH_MOD_V3
