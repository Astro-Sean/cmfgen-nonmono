      MODULE DUST_MOD
      USE SET_KIND_MODULE
!
      LOGICAL MOD_INCL_DUST
      LOGICAL USE_DUST_FILE
      CHARACTER(LEN=80)   DUST_FILE_NAME
      CHARACTER(LEN=3)    DUST_TYPE
      CHARACTER(LEN=30)   DUST_DIST_LAW
      LOGICAL PURE_ABS_DUST
      REAL(KIND=LDP)      FIXED_DUST_TO_GAS_RATIO

      REAL(KIND=LDP)         VMIN_DUST
      REAL(KIND=LDP)         VMAX_DUST
      REAL(KIND=LDP)         DV_DUST
      REAL(KIND=LDP)         V_SHELL_LOC
      REAL(KIND=LDP)         dV_SHELL
      REAL(KIND=LDP)         MASS_DUST
!
      REAL(KIND=LDP), ALLOCATABLE :: LOCAL_DUST_TO_GAS_RATIO(:)
!
! To store dust opacities
! 
      REAL(KIND=LDP), ALLOCATABLE  :: LAM_DUST_IN(:)
      REAL(KIND=LDP), ALLOCATABLE  :: KAP_ABS_DUST_IN(:)
      REAL(KIND=LDP), ALLOCATABLE  :: KAP_SCAT_DUST_IN(:)
      INTEGER  N_DUST_FILe
!
      REAL(KIND=LDP) HZ_TO_MICRON
!
! Local dust temp -- currently assumed to be constant.
!
      REAL(KIND=LDP) DUST_T4 ! if !=0, we add dust emission: chi_dust*BB(dust_temp)
      
      END MODULE DUST_MOD
