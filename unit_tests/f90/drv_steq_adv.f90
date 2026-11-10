! Driver for STEQ_ADVEC_V4 (new_main/subs/steq_advec_v4.f).
!
! Alteration under test: INC_ADV=T with a non-positive V is now
! rejected under DO_NONMONO_SOB/DO_NONMONO_KNOP.  The backward
! difference uses the inner neighbor as the upstream point, which
! is only a valid advection stencil for radial outflow; a field
! that crosses zero would need the opposite stencil.
!
! Species ranges are set empty so the post-guard body is a no-op;
! V lives in MOD_CMFGEN.
!
! Usage: drv CASE
!   nonpos  : SOB + INC_ADV, V(3)<0  -> STOP with label
!   zerov   : KNOP + INC_ADV, V(2)=0 -> STOP with label
!   flagoff : flags off, V(3)<0      -> guard skipped, returns
!   good    : SOB + INC_ADV, V>0     -> passes guard, returns
!   noinc   : SOB, INC_ADV=F, V<0    -> early return before guard
!
      PROGRAM DRV_STEQ_ADV
      USE SET_KIND_MODULE
      USE MOD_CMFGEN
      USE STEQ_DATA_MOD
      USE CONTROL_VARIABLE_MOD, ONLY: &
     &    DO_NONMONO_SOB,DO_NONMONO_KNOP
      IMPLICIT NONE
      INTEGER, PARAMETER :: ND=4,NB=3
      CHARACTER(LEN=16) CASE
      CALL GETARG(1,CASE)
!
      ALLOCATE(R(ND),V(ND),BA_ADV_TERM(NB,ND))
      R=(/ 20.0_LDP,15.0_LDP,10.0_LDP,5.0_LDP /)
      V=(/ 900.0_LDP,800.0_LDP,700.0_LDP,600.0_LDP /)
      NUM_IONS=0
      SPECIES_BEG_ID=1
      SPECIES_END_ID=1
      BA_ADV_TERM=0.0_LDP
!
      DO_NONMONO_SOB=(CASE .EQ. 'nonpos' .OR. CASE .EQ. 'good' &
     &     .OR. CASE .EQ. 'noinc')
      DO_NONMONO_KNOP=(CASE .EQ. 'zerov')
!
      IF(CASE .EQ. 'nonpos' .OR. CASE .EQ. 'flagoff' &
     &     .OR. CASE .EQ. 'noinc')V(3)=-50.0_LDP
      IF(CASE .EQ. 'zerov')V(2)=0.0_LDP
!
      CALL STEQ_ADVEC_V4(1.0_LDP,.TRUE.,NB,ND, &
     &     CASE .NE. 'noinc',.FALSE.,.FALSE.)
!
      WRITE(*,'(A)')'DONE_'
      END
