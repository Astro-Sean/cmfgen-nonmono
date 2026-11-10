! Driver for AUTO_ADD_ION (new_main/subs/auto_add_ion.f).
!
! Alteration under test: the high-stage and no-stage generation calls now
! pass ATM(NID)%NXzV_F; upstream passed ATM(ID)%NXzV_F where ID held the
! stale scan index (last existing ion, or 0 when no stage existed). The
! generated XzV_IN header carries the NLEV argument verbatim, so a stale
! ID wrote the wrong level count.
!
! The driver declares one species (carbon) with three ion stages CIII/CIV/
! CV and distinct NXzV_F per stage. Which _IN files exist is controlled by
! the test harness via the case directory contents; AUTO_ADD_ION_STAGES is
! toggled by the 'off' case argument.
!
      PROGRAM DRV_AUTO_ION
      USE SET_KIND_MODULE
      USE MOD_CMFGEN
      USE CONTROL_VARIABLE_MOD
      IMPLICIT NONE
!
      INTEGER :: I
      CHARACTER(LEN=32) :: CASE
!
      CALL GETARG(1, CASE)
      AUTO_ADD_ION_STAGES = (CASE .NE. 'off')
!
      SPECIES_BEG_ID = 0
      SPECIES_END_ID = 0
      SPECIES = ' '
      ION_ID = ' '
!
! One species with three stages in slots 1..3; every other species entry
! is skipped by the BEG_ID=0 early-exit in AUTO_ADD_ION.
!
      SPECIES(1) = 'CARBON'
      SPECIES_BEG_ID(1) = 1
      SPECIES_END_ID(1) = 4
      ION_ID(1) = 'CIII'
      ION_ID(2) = 'CIV'
      ION_ID(3) = 'CV'
!
! GXZV_F(1) is passed as the ion statistical weight; any positive value
! keeps the ION_POP algebra in SUB_GUESS_DC_V2 finite.
!
      DO I = 1, 4
        IF(.NOT. ALLOCATED(ATM(I)%GXZV_F))ALLOCATE(ATM(I)%GXZV_F(4))
        ATM(I)%GXZV_F = 2.0_LDP
      END DO
      ATM(1)%NXzV_F = 7
      ATM(2)%NXzV_F = 11
      ATM(3)%NXzV_F = 13
!
      CALL AUTO_ADD_ION()
      WRITE(6,*) 'DONE_AUTO_ION'
!
      END
