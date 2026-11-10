! Driver for COMP_GREY_V2 / COMP_GREY_V4
! (new_main/mod_subs/comp_grey_v2.f, comp_grey_v4.f).
!
! Alteration under test: under DO_NONMONO_SOB/DO_NONMONO_KNOP a
! work copy SIGMA_GREY is clamped to >= -0.999 before being
! handed to the grey solvers -- their monotonic-expansion moment
! structures have 1+sigma*g terms that sign-flip or divide near
! a pole when sigma < -1.  SIGMA and SIGMA_TRUE stay untouched;
! stock runs are bit-identical (plain copy, no message).  The
! clamp count is reported at CONV_DIAG_LVL >= 1.
!
! The message is the observable: it proves the block ran and how
! many entries were clamped.  The grey solve that follows is not
! part of this alteration.
!
! Usage: drv CASE
!   clamp : KNOP + CONV_DIAG_LVL=1, two sigmas < -0.999 -> message
!   quiet : same + CONV_DIAG_LVL=0  -> no message
!   pos   : KNOP, all sigma >= -0.999 -> no message
!   off   : flags off, clamped-needed sigma -> no message
!
      PROGRAM DRV_GREY
      USE SET_KIND_MODULE
      USE ANG_QW_MOD
      USE MOD_CMFGEN
      USE CONTROL_VARIABLE_MOD
      IMPLICIT NONE
      INTEGER, PARAMETER :: ND=4,NC=2,NP=8,NT=4
      REAL(KIND=LDP) POPS(NT,ND),TGREY(ND)
      REAL(KIND=LDP) TAU_ROSS(ND),ROSSMEAN(ND)
      REAL(KIND=LDP) PLANCKMEAN(ND)
      LOGICAL COMPUTED
      CHARACTER(LEN=16) CASE
      CALL GETARG(1,CASE)
!
      ALLOCATE(R(ND),V(ND),SIGMA(ND),SIGMA_TRUE(ND))
      ALLOCATE(CLUMP_FAC(ND))
      ALLOCATE(P(NP),AQW(ND,NP),HQW(ND,NP),KQW(ND,NP))
      ALLOCATE(NQW(ND,NP),HMIDQW(ND,NP),NMIDQW(ND,NP))
      R=(/ 20.0_LDP,15.0_LDP,10.0_LDP,5.0_LDP /)
      V=(/ 900.0_LDP,800.0_LDP,700.0_LDP,600.0_LDP /)
      SIGMA=(/ -0.5_LDP,-1.5_LDP,-2.0_LDP,-0.9_LDP /)
      SIGMA_TRUE=SIGMA
      IF(CASE .EQ. 'pos')SIGMA=(/ -0.5_LDP,-0.8_LDP, &
     &                            -0.99_LDP,-0.3_LDP /)
      POPS=1.0_LDP; TGREY=0.0_LDP
      TAU_ROSS=0.0_LDP; ROSSMEAN=0.1_LDP
      PLANCKMEAN=0.1_LDP; CLUMP_FAC=1.0_LDP
      COMPUTED=.FALSE.
!
      DO_NONMONO_SOB=.FALSE.
      DO_NONMONO_KNOP=(CASE .NE. 'off')
      CONV_DIAG_LVL=1
      IF(CASE .EQ. 'quiet')CONV_DIAG_LVL=0
      PLANE_PARALLEL_NO_V=.FALSE.; PLANE_PARALLEL=.FALSE.
      JGREY_WITH_V_TERMS=.FALSE.
      USE_J_REL=.FALSE.; USE_DJDT_RTE=.FALSE.
      DIF=.FALSE.; THK_CONT=.FALSE.
      METHOD='ZERO  '; LUM=1.0_LDP
!
      IF(CASE .EQ. 'v4')THEN
        CALL COMP_GREY_V4(POPS,TGREY,TAU_ROSS,ROSSMEAN, &
     &       PLANCKMEAN,COMPUTED,6,NC,ND,NP,NT)
      ELSE
        CALL COMP_GREY_V2(POPS,TGREY,TAU_ROSS,ROSSMEAN, &
     &       6,NC,ND,NP,NT)
      END IF
!
      WRITE(*,'(A,4ES14.6)')'SIG ',SIGMA(1:4)
      WRITE(*,'(A)')'DONE_'
      END
