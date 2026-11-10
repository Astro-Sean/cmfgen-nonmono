! Driver for ADJUST_R_GRID_V4 (new_main/mod_subs/adjust_r_grid_v4.f).
!
! Alterations under test:
!   - GRID_METH=VEL_SPACE is rejected under DO_NONMONO_SOB /
!     DO_NONMONO_KNOP because the R(V) interpolation in
!     DO_VEL_REGRID_V2 assumes a monotonic velocity field.
!   - the interpolated-SIGMA < -1 clamp/STOP is skipped under the
!     same flags: sigma < -1 is expected across a velocity
!     reversal, and RECOMPUTE_SIGMA_TRUE / the KNOP formal solver
!     handle it downstream.
!
! The driver fills MOD_CMFGEN arrays and writes ADJUST_R_DEFAULTS
! in cwd (run each case in a scratch dir).  OB_ONLY nudges R(2)
! so the regrid path always reaches the sigma block.
!
! Usage: drv CASE
!   velsob : VEL_SPACE + DO_NONMONO_SOB -> gate STOP
!   velknop: VEL_SPACE + DO_NONMONO_KNOP -> gate STOP
!   sigoff : OB_ONLY + reversal V + flags off -> 'SIGMA is not
!            monotonic' STOP (upstream check preserved)
!   sigsob : same deck + DO_NONMONO_SOB -> completes
!   sigknop: same deck + DO_NONMONO_KNOP -> completes
!
      PROGRAM DRV_ADJ_R4
      USE SET_KIND_MODULE
      USE MOD_CMFGEN
      USE CONTROL_VARIABLE_MOD, ONLY: DO_NONMONO_SOB,DO_NONMONO_KNOP
      IMPLICIT NONE
      INTEGER, PARAMETER :: ND=8, NT=1
      REAL(KIND=LDP) POPS(NT,ND),ESEC(ND)
      INTEGER MAIN_COUNTER
      LOGICAL DONE_R_REV
      CHARACTER(LEN=16) CASE
      INTEGER I,LU
      CALL GETARG(1,CASE)
      ALLOCATE(R(ND),V(ND),SIGMA(ND),POP_ATOM(ND))
      ALLOCATE(ROSS_MEAN(ND),CLUMP_FAC(ND))
      DO I=1,ND
        R(I)=45.0_LDP-5.0_LDP*I
        V(I)=850.0_LDP-100.0_LDP*I
        SIGMA(I)=-0.5_LDP
        POP_ATOM(I)=1.0_LDP
        ROSS_MEAN(I)=1.0_LDP
        CLUMP_FAC(I)=1.0_LDP
        ESEC(I)=1.0_LDP
        POPS(1,I)=1.0_LDP+0.1_LDP*I
      END DO
      V(ND-1)=600.0_LDP
      V(ND)=1400.0_LDP
!
      DO_NONMONO_SOB=(CASE .EQ. 'velsob' .OR. CASE .EQ. 'sigsob')
      DO_NONMONO_KNOP=(CASE .EQ. 'velknop' .OR. CASE .EQ. 'sigknop')
!
      LU=21
      OPEN(UNIT=LU,FILE='ADJUST_R_DEFAULTS',STATUS='REPLACE')
      WRITE(LU,'(A)')'   5 [N_ITS]'
      WRITE(LU,'(A)')'   1 [STRT_ITS]'
      WRITE(LU,'(A)')'   1 [FREQ_ITS]'
      IF(CASE(1:3) .EQ. 'vel')THEN
        WRITE(LU,'(A)')'VEL_SPACE [GRID_METH]'
      ELSE
        WRITE(LU,'(A)')'OB_ONLY [GRID_METH]'
      END IF
      WRITE(LU,'(A)')'DEFAULT [OB_OPT]'
      WRITE(LU,'(A)')'DEFAULT [IB_OPT]'
      WRITE(LU,'(A)')'1.0D0 [D2OND1]'
      CLOSE(LU)
!
      MAIN_COUNTER=1
      DONE_R_REV=.FALSE.
      CALL ADJUST_R_GRID_V4(POPS,ESEC,MAIN_COUNTER,DONE_R_REV,ND,NT)
      WRITE(*,'(A,L1)')'ADJ_ DONE_R_REV=',DONE_R_REV
      WRITE(*,'(A,ES16.8)')'ADJ_ SIGMA_NDm1=',SIGMA(ND-1)
      END
