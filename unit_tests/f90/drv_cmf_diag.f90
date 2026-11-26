! Driver for the CMF_BLKBAND_V3 DIAG-path solve contract
! (checklist WP-7b, 27-Nov-2026).
!
! Same two-ion species fixture as drv_gfm_nc: ion 1 has level rows
! 1 (g.s.), 2 and 3 plus NC at eq 4; pseudo ion 2 sits on
! EQ_SPECIES=4; ED/T at rows 5/6.  SOL_TYPE='DIAG', BA_COMPUTED=F
! and WR_BA_INV=F route the call straight into the compute path
! (CMF_DIAG_SOLVE): GENERATE_FULL_MATRIX_V3 builds each depth
! block, then DGEEQU/DGETRF/DGECON/DGETRS solve it.
!
! Cases:
!   ok     : both depths solve -> FLAG=T
!   badrow : ion 1's row-3 rate eq is identically zero at depth 2
!            with a nonzero residual -> DGEEQU returns IFAIL=3,
!            the depth is counted failed and FLAG=F (pre-WP-7b the
!            depth zeroed STEQ and still returned ok=T).
!
      PROGRAM DRV_CMF_DIAG
      USE SET_KIND_MODULE
      USE MOD_CMFGEN
      USE STEQ_DATA_MOD
      USE CONTROL_VARIABLE_MOD
      IMPLICIT NONE
      INTEGER, PARAMETER :: NT=6,ND=2,NION=2,NB=3,DI=2
      REAL(KIND=LDP) STEQ(NT,ND),POPS(NT,ND)
      LOGICAL FLAG
      CHARACTER(LEN=16) CASE
      INTEGER J
      CALL GETARG(1,CASE)
!
      ALLOCATE(BA_ED(NT,NB,ND),BA_T(NT,NB,ND))
      ALLOCATE(STEQ_ED(ND),STEQ_T(ND))
      ALLOCATE(BA_ADV_TERM(NB,ND))
      SPECIES_BEG_ID(1)=1; SPECIES_END_ID(1)=3
      SPECIES_PRES=.FALSE.; SPECIES_PRES(1)=.TRUE.
      SPECIES_LNK(1)=1; SPECIES_LNK(2)=1
      EQ_SPECIES(1)=4
      FIX_SPECIES(1)=0
      LTE_MODEL=.FALSE.; USE_ELEC_HEAT_BAL=.FALSE.
      JFNK_QUIET_IO=.TRUE.; POP_AUDIT=.FALSE.; BA_SYM_PIN=.FALSE.
      MOD_FIXED_NE=.FALSE.; MOD_FIXED_T=.FALSE.
      FIX_IN_BOUND_T=.FALSE.; MOD_FIX_IMPURITY=.FALSE.
      BA_ADV_TERM=0.0_LDP
      STEQ_ED=0.03_LDP; STEQ_T=0.04_LDP
!
! Independent ED/T rows so the block is nonsingular.
!
      DO J=1,ND
        BA_ED(:,DI,J)=[0.01_LDP,0.02_LDP,0.03_LDP,0.01_LDP,0.5_LDP, &
     &                 0.01_LDP]
        BA_T(:,DI,J)=[0.02_LDP,0.01_LDP,0.02_LDP,0.03_LDP,0.01_LDP, &
     &                0.7_LDP]
        POPS(1,J)=1.0_LDP
        POPS(2,J)=5.0_LDP
        POPS(3,J)=10.0_LDP
        POPS(4,J)=0.5_LDP
        POPS(5,J)=2.0_LDP
        POPS(6,J)=2.0_LDP
      END DO
!
      ATM(1)%EQXzV=1; ATM(1)%NXzV=2; ATM(1)%XzV_PRES=.TRUE.
      ATM(1)%FIX_NXzV=0
      ALLOCATE(ATM(1)%XzV(2,ND),ATM(1)%DXzV(ND))
      ATM(1)%XzV=0.0_LDP; ATM(1)%DXzV=0.0_LDP
      SE(1)%N_SE=4; SE(1)%N_IV=3; SE(1)%XRAY_EQ=0
      ALLOCATE(SE(1)%EQ_IN_BA(4),SE(1)%LNK_TO_F(3))
      ALLOCATE(SE(1)%BA(4,3,NB,ND),SE(1)%STEQ(4,ND))
      ALLOCATE(SE(1)%STEQ_ADV(ND))
      SE(1)%EQ_IN_BA=[1,2,3,4]
      SE(1)%LNK_TO_F=[1,2,3]
      SE(1)%BA=0.0_LDP; SE(1)%STEQ=0.0_LDP; SE(1)%STEQ_ADV=0.0_LDP
      DO J=1,3
        SE(1)%BA(1,J,DI,:)=11.0_LDP+J
        SE(1)%BA(2,J,DI,:)=21.0_LDP+J
        SE(1)%BA(3,J,DI,:)=31.0_LDP+J
        SE(1)%BA(4,J,DI,:)=91.0_LDP+J
      END DO
      SE(1)%STEQ(1,:)=1.1_LDP
      SE(1)%STEQ(2,:)=1.2_LDP
      SE(1)%STEQ(3,:)=1.3_LDP
      SE(1)%STEQ(4,:)=1.9_LDP
!
      ATM(2)%EQXzV=4; ATM(2)%NXzV=1; ATM(2)%XzV_PRES=.FALSE.
      ATM(2)%FIX_NXzV=0
      SE(2)%N_SE=2; SE(2)%N_IV=1; SE(2)%XRAY_EQ=0
      ALLOCATE(SE(2)%EQ_IN_BA(2),SE(2)%LNK_TO_F(1))
      ALLOCATE(SE(2)%BA(2,1,NB,ND),SE(2)%STEQ(2,ND))
      ALLOCATE(SE(2)%STEQ_ADV(ND))
      SE(2)%EQ_IN_BA=[4,4]
      SE(2)%LNK_TO_F=[4]
      SE(2)%BA=0.0_LDP; SE(2)%STEQ=0.0_LDP; SE(2)%STEQ_ADV=0.0_LDP
      SE(2)%BA(1,1,DI,:)=44.0_LDP
      SE(2)%BA(2,1,DI,:)=95.0_LDP
      SE(2)%STEQ(1,:)=2.4_LDP
      SE(2)%STEQ(2,:)=2.9_LDP
!
! badrow: ion 1's row-3 rate eq vanishes at depth 2 while its
! residual stays nonzero -> DGEEQU reports row 3 exactly zero.
!
      IF(CASE .EQ. 'badrow')SE(1)%BA(3,:,DI,2)=0.0_LDP
!
      STEQ=0.0_LDP
      CALL CMF_BLKBAND_V3(STEQ,POPS,'DIAG',FLAG,DI,NT,NION,NB,ND, &
     &                    .FALSE.,.FALSE.,.FALSE.)
!
      WRITE(*,'(A,L2)')'FLAG=',FLAG
      WRITE(*,'(A)')'STEQ'
      DO J=1,ND
        WRITE(*,'(6ES24.16)')STEQ(:,J)
      END DO
      WRITE(*,'(A)')'DONE_'
      END
