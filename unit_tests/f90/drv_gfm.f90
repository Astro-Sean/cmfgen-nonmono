! Driver for GENERATE_FULL_MATRIX_V3
! (new_main/subs/generate_full_matrix_v3.f).
!
! Alteration under test: the BA_ASCI diagnostic dumps at K=5,
! K=41 and K=ND are now gated on .NOT. JFNK_QUIET_IO.  The JFNK
! assembly paths (JFNK_ASSEMBLE_STEQ_V1/JFNK_JVP_V1) re-enter this
! routine once per residual eval and Krylov matvec, so an ungated
! dump rewrites ~100 MB of scratch every eval; the production
! SOLVEBA writes are unchanged.
!
! The species ranges are set empty and USE_PASSED_REP=T is passed
! so the equation-assembly body is a no-op; the dump blocks at
! the end still run when DIAG_BAND and the K value match.
!
! Usage: drv CASE
!   loud  : JFNK_QUIET_IO=F -> BA_ASCI_N_D5 + BA_ASCI_N_DND written
!   quiet : JFNK_QUIET_IO=T -> no BA_ASCI_* files written
!
      PROGRAM DRV_GFM
      USE SET_KIND_MODULE
      USE MOD_CMFGEN
      USE STEQ_DATA_MOD
      USE CONTROL_VARIABLE_MOD
      IMPLICIT NONE
      INTEGER, PARAMETER :: NT=4,ND=6,NION=1,NB=3,DI=2
      REAL(KIND=LDP) C_MAT(NT,NT),STEQ_VEC(NT),POPS(NT,ND)
      LOGICAL ZERO_STEQ(NT),REPLACE(NION)
      CHARACTER(LEN=16) CASE
      CALL GETARG(1,CASE)
!
      ALLOCATE(BA_ED(NT,NB,ND),BA_T(NT,NB,ND))
      ALLOCATE(STEQ_ED(ND),STEQ_T(ND))
      ALLOCATE(BA_ADV_TERM(NB,ND))
      SPECIES_BEG_ID=1; SPECIES_END_ID=1
      SPECIES_PRES=.FALSE.
      LTE_MODEL=.FALSE.; USE_ELEC_HEAT_BAL=.FALSE.
      JFNK_QUIET_IO=(CASE .EQ. 'quiet')
      BA_ED=0.01_LDP; BA_T=0.02_LDP
      STEQ_ED=0.03_LDP; STEQ_T=0.04_LDP
      BA_ADV_TERM=0.0_LDP
      POPS=1.0_LDP; C_MAT=0.0_LDP; STEQ_VEC=0.0_LDP
      ZERO_STEQ=.FALSE.; REPLACE=.FALSE.
!
! K=5 dump then K=ND=6 dump; each depth's first call is the
! diagonal band (DIAG_INDX=BAND_INDX=DI).
!
      CALL GENERATE_FULL_MATRIX_V3(C_MAT,STEQ_VEC,POPS,REPLACE, &
     &     ZERO_STEQ,NT,ND,NION,NB,DI,DI,5,.TRUE.,.FALSE.,.TRUE.)
      CALL GENERATE_FULL_MATRIX_V3(C_MAT,STEQ_VEC,POPS,REPLACE, &
     &     ZERO_STEQ,NT,ND,NION,NB,DI,DI,6,.FALSE.,.TRUE.,.TRUE.)
!
      WRITE(*,'(A)')'DONE_'
      END
