! Driver for FIDDLE_POP_CORRECTIONS_V2
! (new_main/subs/fiddle_pop_corrections_v2.f).
!
! Alterations under test:
! 1. Under DO_NONMONO_SOB/DO_NONMONO_KNOP each row's contribution
!    to the SCALE scan is capped at the depth's apply-clamp bound
!    (DPTH_BIG_LIM/DPTH_LIT_LIM).  Upstream let a huge demanded
!    correction (e.g. STEQ=-1E70 from an ill-conditioned KNOP
!    cold-start row) bound SCALE to ~1E-69, freezing every other
!    correction at the depth.  The apply clamp enforces the same
!    bound on that row anyway, so capping it in the scan is
!    redundant for it and rescues the rest.
! 2. Under DO_SOB_FG_KNOP non-positive/non-finite corrected POPS
!    are floored: species to 1.0D-100, the T row to T_MIN.  A NaN
!    STEQ row otherwise wrote a permanent NaN into POPS.
!
! Depth 1 carries a 1E70 decrease demand; depth 2 carries NaN in
! the fgknop/stocknan cases.
!
! Usage: drv CASE
!   stock    : flags off        -> SCALE~1E-70, only the big row
!                                 moves (to its own clamp)
!   nonmono  : DO_NONMONO_KNOP  -> SCALE~0.80, other rows correct
!   fgknop   : KNOP+SOB_FG_KNOP -> as nonmono + NaN depth floored
!   stocknan : flags off + NaN  -> NaN propagates into POPS
!
! W5 (24-Nov-2026): ADJUST_CORRECTIONS validation cases.  These
! write a file in cwd before the call and use small uniform STEQ
! so SCALE=1 and RELAX is the observable.  Flags stay off.
!   relax0    : RELAX=0    -> warn + applied 1.0 (full correction)
!   relaxhalf : RELAX=0.5  -> half correction at every depth
!   relaxrange: RELAX=0.5, L_ST=1,L_END=1 -> depth1 half, depth2 full
!   badorder  : L_ST=2>L_END=1 -> warn + defaults on full range
!   lbeyond   : L_END=ND+1     -> warn + defaults on full range
!   lamadj    : RELAX=0.5 + LAMBDA_IT -> override to 1.0, full corr
!
      PROGRAM DRV_FPC
      USE, INTRINSIC :: IEEE_ARITHMETIC
      USE SET_KIND_MODULE
      USE CONTROL_VARIABLE_MOD, ONLY: &
     &    DO_NONMONO_SOB,DO_NONMONO_KNOP,DO_SOB_FG_KNOP
      IMPLICIT NONE
      INTEGER, PARAMETER :: NT=4,ND=2
      REAL(KIND=LDP) POPS(NT,ND),STEQ(NT,ND)
      REAL(KIND=LDP) T_MIN,QNAN
      CHARACTER(LEN=16) CASE
      CALL GETARG(1,CASE)
      QNAN=IEEE_VALUE(0.0_LDP,IEEE_QUIET_NAN)
!
      T_MIN=0.4_LDP
      POPS=1.0_LDP
      POPS(NT,:)=3.0_LDP
      STEQ=0.005_LDP
      STEQ(3,:)=-0.005_LDP
      STEQ(4,:)=0.001_LDP
      IF(CASE .EQ. 'stock' .OR. CASE .EQ. 'nonmono' .OR. &
     &   CASE .EQ. 'fgknop' .OR. CASE .EQ. 'stocknan') &
     &  STEQ(1,1)=1.0E+70_LDP
      IF(CASE .EQ. 'fgknop' .OR. CASE .EQ. 'stocknan')THEN
        STEQ(:,2)=QNAN
      END IF
!
      DO_NONMONO_SOB=.FALSE.
      DO_NONMONO_KNOP=(CASE .EQ. 'nonmono' .OR. CASE .EQ. 'fgknop')
      DO_SOB_FG_KNOP=(CASE .EQ. 'fgknop')
!
      IF(CASE .EQ. 'relax0')THEN
        OPEN(UNIT=47,FILE='ADJUST_CORRECTIONS',STATUS='REPLACE')
        WRITE(47,'(A)')'0.0   [RELAX]'
        CLOSE(47)
      ELSE IF(CASE .EQ. 'relaxhalf' .OR. CASE .EQ. 'lamadj')THEN
        OPEN(UNIT=47,FILE='ADJUST_CORRECTIONS',STATUS='REPLACE')
        WRITE(47,'(A)')'0.5   [RELAX]'
        CLOSE(47)
      ELSE IF(CASE .EQ. 'relaxrange')THEN
        OPEN(UNIT=47,FILE='ADJUST_CORRECTIONS',STATUS='REPLACE')
        WRITE(47,'(A)')'0.5   [RELAX]'
        WRITE(47,'(A)')'1     [L_ST]'
        WRITE(47,'(A)')'1     [L_END]'
        CLOSE(47)
      ELSE IF(CASE .EQ. 'badorder')THEN
        OPEN(UNIT=47,FILE='ADJUST_CORRECTIONS',STATUS='REPLACE')
        WRITE(47,'(A)')'0.5   [RELAX]'
        WRITE(47,'(A)')'2     [L_ST]'
        WRITE(47,'(A)')'1     [L_END]'
        CLOSE(47)
      ELSE IF(CASE .EQ. 'lbeyond')THEN
        OPEN(UNIT=47,FILE='ADJUST_CORRECTIONS',STATUS='REPLACE')
        WRITE(47,'(A)')'0.5   [RELAX]'
        WRITE(47,'(A)')'3     [L_END]'
        CLOSE(47)
      END IF
!
      CALL FIDDLE_POP_CORRECTIONS_V2(POPS,STEQ,T_MIN,5.0_LDP, &
     &     0.2_LDP,'MAJOR',CASE .EQ. 'lamadj',11,NT,ND)
!
      WRITE(*,'(A,4ES14.6)')'P1 ',POPS(1:4,1)
      WRITE(*,'(A,4ES14.6)')'P2 ',POPS(1:4,2)
      WRITE(*,'(A)')'DONE_'
      END
