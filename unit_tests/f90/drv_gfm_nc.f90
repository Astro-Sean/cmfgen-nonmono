! Driver for the [BA_NC_MAXPOP] conservation-row placement in
! GENERATE_FULL_MATRIX_V3 (checklist WP-7c, 27-Nov-2026).
!
! Seeds one species of two ions: ion 1 with three level rows
! (EQXzV=1 ground state, excited row 2, super-level row 3) and a
! pseudo top ion whose only row is EQ_SPECIES=4. ED sits at row 5
! and T at row 6. POPS at the solve depth puts the argmax on row 3
! so the armed path must move the NC row there and relocate row 3's
! rate equation to row 4.
!
! Usage: drv CASE
!   ncdef : BA_NC_MAXPOP=F -> NC at row 4, row 3 keeps its rate eq
!   ncmax : BA_NC_MAXPOP=T -> NC at row 3, row 3 rate eq at row 4
!
      PROGRAM DRV_GFM_NC
      USE SET_KIND_MODULE
      USE MOD_CMFGEN
      USE STEQ_DATA_MOD
      USE CONTROL_VARIABLE_MOD
      IMPLICIT NONE
      INTEGER, PARAMETER :: NT=6,ND=6,NION=2,NB=3,DI=2
      INTEGER, PARAMETER :: KD=3
      REAL(KIND=LDP) C_MAT(NT,NT),STEQ_VEC(NT),POPS(NT,ND)
      LOGICAL ZERO_STEQ(NT),REPLACE(NION)
      CHARACTER(LEN=16) CASE
      INTEGER J
      CALL GETARG(1,CASE)
!
      ALLOCATE(BA_ED(NT,NB,ND),BA_T(NT,NB,ND))
      ALLOCATE(STEQ_ED(ND),STEQ_T(ND))
      ALLOCATE(BA_ADV_TERM(NB,ND))
      SPECIES_BEG_ID(1)=1; SPECIES_END_ID(1)=3
      SPECIES_PRES=.FALSE.; SPECIES_PRES(1)=.TRUE.
      EQ_SPECIES(1)=4
      LTE_MODEL=.FALSE.; USE_ELEC_HEAT_BAL=.FALSE.
      JFNK_QUIET_IO=.TRUE.
      BA_NC_MAXPOP=(CASE .EQ. 'ncmax')
      BA_ED=0.01_LDP; BA_T=0.02_LDP
      STEQ_ED=0.03_LDP; STEQ_T=0.04_LDP
      BA_ADV_TERM=0.0_LDP
      C_MAT=0.0_LDP; STEQ_VEC=0.0_LDP
      ZERO_STEQ=.FALSE.; REPLACE=.FALSE.
!
! Populations: row 3 is the most populated candidate at every depth
! (ground row 1 and pseudo row 4 are excluded by construction).
!
      DO J=1,ND
        POPS(1,J)=1.0_LDP
        POPS(2,J)=5.0_LDP
        POPS(3,J)=10.0_LDP
        POPS(4,J)=0.5_LDP
        POPS(5,J)=2.0_LDP
        POPS(6,J)=2.0_LDP
      END DO
!
! Ion 1: N_SE=4 -> rate eqs at EQ_IN_BA(1..3)=(1,2,3), NC at index 4.
!
      ATM(1)%EQXzV=1; ATM(1)%NXzV=2; ATM(1)%XzV_PRES=.TRUE.
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
! Ion 2 (pseudo top ion): N_SE=2 -> rate eq at EQ_IN_BA(1)=4
! (overlaid by EQ_SPECIES), NC contribution at index 2.
!
      ATM(2)%EQXzV=4; ATM(2)%NXzV=1; ATM(2)%XzV_PRES=.FALSE.
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
      CALL GENERATE_FULL_MATRIX_V3(C_MAT,STEQ_VEC,POPS,REPLACE, &
     &     ZERO_STEQ,NT,ND,NION,NB,DI,DI,KD,.TRUE.,.FALSE.,.TRUE.)
!
      WRITE(*,'(A)')'C_MAT'
      DO J=1,NT
        WRITE(*,'(6ES24.16)')C_MAT(J,:)
      END DO
      WRITE(*,'(A)')'STEQ_VEC'
      WRITE(*,'(6ES24.16)')STEQ_VEC(:)
      WRITE(*,'(A)')'DONE_'
      END
