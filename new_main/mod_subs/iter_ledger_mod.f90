! iter_ledger_mod.f90 -- per-great-iteration work ledger for the
! non-monotonic iteration cycle (audit task T-06; findings F-22/F-24:
! nothing in the tree measured whether refresh passes, JFNK evals or
! the frequency sweep dominate wall time).
!
! LEDG_ITER_RESET is armed once per great iteration from CMFGEN_SUB
! (label 20000) under DO_NONMONO_SOB/DO_NONMONO_KNOP; the call-site
! counters and timers below then accumulate over the sweep, and
! SOLVE_FOR_POPS appends one CSV record per call to ITER_LEDGER next
! to the ITER_HISTORY write. LEDG_ACTIVE stays .FALSE. on the stock
! path, so every instrumented site short-circuits and no file is
! created -- "zero output when nonmono flags are off" (T-06
! acceptance).
!
! Counter scope:
!   LEDG_N_FG / LEDG_N_MOM - calls to the continuum J solvers from
!     COMP_J_BLANK (FG_J_CMF_V13 / FG_J_CMF_NONMONO_V1 /
!     MOM_J_CMF_V11 / MOM_J_CMF_NONMONO_V1). Includes the EDDFAC
!     extended-grid passes and JFNK_RAD_SWEEP_V1 replayed solves
!     (both run inside the counted call sites); excludes the
!     PP/DDT/JREL branches, which the nonmono consistency checks
!     reject, and refresh-internal re-drives (counted in np_rf).
!   LEDG_N_RF / LEDG_NP_RF - GS refresh calls and the total sweep
!     passes they ran (per-ray passes summed for FG, IPASS for MOM).
!   LEDG_N_DEM - LAMBDA demotions decided this call, derived at write
!     time by scanning ITER_REASON for '*demote' tokens, so new
!     demote reasons are counted automatically.
!   LEDG_N_FRZ - freeze entries (GS_FREEZE false->true edge),
!     incremented where GS_WAS_FROZEN is updated.
! Timer scope:
!   t_iter - wall s from label 20000 to the ledger write
!   t_rad  - wall s inside the COMP_J_BLANK solver calls
!   t_rf   - wall s inside FG/MOM_KNOP_REFRESH_V1 (includes
!            refresh-internal solver re-drives)
!   t_pops - wall s inside SOLVE_FOR_POPS (includes JFNK evals, which
!            also feed t_rad through their replayed sweeps)
!
! Edited by Devin (T-06, 04-Nov-2026)
!
      MODULE ITER_LEDGER_MOD
      USE SET_KIND_MODULE
      IMPLICIT NONE
      SAVE
!
      LOGICAL :: LEDG_ACTIVE=.FALSE.
!
! Iteration wall origin and stage accumulators (OMP_GET_WTIME
! seconds). LEDG_T0=-1 marks "no armed iteration" so a write before
! the first reset cannot fabricate a negative time.
!
      REAL(KIND=LDP) :: LEDG_T0=-1.0_LDP
      REAL(KIND=LDP) :: LEDG_T_RAD=0.0_LDP
      REAL(KIND=LDP) :: LEDG_T_RF=0.0_LDP
      REAL(KIND=LDP) :: LEDG_T_P0=-1.0_LDP
!
      INTEGER*8 :: LEDG_N_FG=0
      INTEGER*8 :: LEDG_N_MOM=0
      INTEGER*8 :: LEDG_N_RF=0
      INTEGER*8 :: LEDG_NP_RF=0
      INTEGER*8 :: LEDG_N_FRZ=0
!
! JFNK eval counters snapshotted at iteration start; the ledger row
! carries the delta so eval counts attribute to the iteration that
! ran them even though they fire inside SOLVE_FOR_POPS.
!
      INTEGER :: LEDG_JEV0=0
      INTEGER :: LEDG_JJV0=0
!
      CONTAINS
!
! Arm the per-iteration accumulators. Called once per great
! iteration under the nonmono flags; WT is the caller's
! OMP_GET_WTIME() so the module stays OpenMP-free, and JEV/JJV are
! the JFNK interface's cumulative eval counters.
!
      SUBROUTINE LEDG_ITER_RESET(WT,JEV,JJV)
      REAL(KIND=LDP) WT
      INTEGER JEV,JJV
      LEDG_ACTIVE=.TRUE.
      LEDG_T0=WT
      LEDG_T_RAD=0.0_LDP
      LEDG_T_RF=0.0_LDP
      LEDG_T_P0=-1.0_LDP
      LEDG_N_FG=0
      LEDG_N_MOM=0
      LEDG_N_RF=0
      LEDG_NP_RF=0
      LEDG_N_FRZ=0
      LEDG_JEV0=JEV
      LEDG_JJV0=JJV
      RETURN
      END SUBROUTINE LEDG_ITER_RESET
!
! Peak resident set in MB from /proc/self/status VmHWM (Linux).
! -1 when unreadable (non-Linux or a restricted procfs); the ledger
! column then reads -1 rather than a fabricated value.
!
      SUBROUTINE LEDG_RSS_MB(RSS)
      REAL(KIND=LDP) RSS
      INTEGER LU,IOS
      REAL(KIND=LDP) TKB
      CHARACTER(LEN=96) LINE
      RSS=-1.0_LDP
      OPEN(NEWUNIT=LU,FILE='/proc/self/status',STATUS='OLD', &
           ACTION='READ',IOSTAT=IOS)
      IF(IOS .NE. 0)RETURN
      DO
        READ(LU,'(A)',IOSTAT=IOS)LINE
        IF(IOS .NE. 0)EXIT
        IF(INDEX(LINE,'VmHWM:') .EQ. 1)THEN
          READ(LINE(7:),*,IOSTAT=IOS)TKB
          IF(IOS .EQ. 0)RSS=TKB/1024.0_LDP
          EXIT
        END IF
      END DO
      CLOSE(LU)
      RETURN
      END SUBROUTINE LEDG_RSS_MB
!
! Count occurrences of TOKEN inside REASON (ITER_REASON packs
! '+demote' tokens with '+' separators, so a substring count counts
! the events; tokens share the '_demote' suffix by convention).
!
      INTEGER FUNCTION LEDG_COUNT_TOKEN(REASON,TOKEN)
      CHARACTER(LEN=*) REASON,TOKEN
      INTEGER P,L
      LEDG_COUNT_TOKEN=0
      P=1
      L=LEN(TOKEN)
      DO
        P=INDEX(REASON(P:),TOKEN)
        IF(P .LE. 0)EXIT
        LEDG_COUNT_TOKEN=LEDG_COUNT_TOKEN+1
        P=P+L
        IF(P .GT. LEN(REASON))EXIT
      END DO
      RETURN
      END FUNCTION LEDG_COUNT_TOKEN
!
! Append one CSV record to ITER_LEDGER in the run directory. Written
! next to ITER_HISTORY in SOLVE_FOR_POPS so the row shares its
! emission point; header emitted only when the file is new (the file
! survives requeues like ITER_HISTORY). WT is the caller's
! OMP_GET_WTIME() at write time.
!
! Column contract (v1):
!  iter,mode,next,ba,frozen,held,reason,
!  n_fg,n_mom,n_rf,np_rf,ba_build,ba_reuse,
!  jfnk,jfnk_cpl,jfnk_ev,jfnk_jvp,n_rej,n_dem,n_frz,
!  n_frozen,n_rel,n_raw_rel,
!  maxch,maxch_app,maxatt,se_res,se_res2,trust,
!  t_iter,t_rad,t_rf,t_pops,rss_mb
!
      SUBROUTINE LEDG_WRITE_ITER(WT, &
           ITER,MODE_NOW,MODE_NEXT,BA_CH,REASON,FROZEN,HELD, &
           N_BA_BUILD,N_BA_REUSE,JADM,JCOUPLED,JEV_NOW,JJV_NOW, &
           N_REJ,N_FROZEN,N_REL,N_RAW_REL, &
           MAXCH,MAXCH_APP,MAXATT,SE_RES,SE_RES2,TRUST)
      REAL(KIND=LDP) WT
      INTEGER ITER,N_BA_BUILD,N_BA_REUSE
      INTEGER N_REJ,N_FROZEN,N_REL
      INTEGER*8 N_RAW_REL
      INTEGER JEV_NOW,JJV_NOW
      LOGICAL FROZEN,HELD,JADM,JCOUPLED
      REAL(KIND=LDP) MAXCH,MAXCH_APP,MAXATT,SE_RES,SE_RES2,TRUST
      CHARACTER(LEN=*) MODE_NOW,MODE_NEXT,BA_CH,REASON
!
      INTEGER LU,IOS
      LOGICAL FEXIST
      REAL(KIND=LDP) RSS,T_ITER,T_POPS
      INTEGER F_I,H_I,J_I,C_I,R_I
!
! Calls before LEDG_ITER_RESET (stock path, or a pre-iteration
! SOLVE_FOR_POPS) leave no trace: no file, no row.
!
      IF(.NOT. LEDG_ACTIVE)RETURN
!
      F_I=0
      H_I=0
      J_I=0
      C_I=0
      R_I=0
      IF(FROZEN)F_I=1
      IF(HELD)H_I=1
      IF(JADM)J_I=1
      IF(JCOUPLED)C_I=1
      IF(N_REJ .GT. 0)R_I=N_REJ
!
      T_ITER=-1.0_LDP
      IF(LEDG_T0 .GE. 0.0_LDP)T_ITER=WT-LEDG_T0
      T_POPS=-1.0_LDP
      IF(LEDG_T_P0 .GE. 0.0_LDP)T_POPS=WT-LEDG_T_P0
      CALL LEDG_RSS_MB(RSS)
!
      INQUIRE(FILE='ITER_LEDGER',EXIST=FEXIST)
      OPEN(NEWUNIT=LU,FILE='ITER_LEDGER',STATUS='UNKNOWN', &
           POSITION='APPEND',IOSTAT=IOS)
      IF(IOS .NE. 0)RETURN
      IF(.NOT. FEXIST)WRITE(LU,'(A)') &
        '# iter,mode,next,ba,frozen,held,reason,'// &
        'n_fg,n_mom,n_rf,np_rf,ba_build,ba_reuse,'// &
        'jfnk,jfnk_cpl,jfnk_ev,jfnk_jvp,n_rej,n_dem,n_frz,'// &
        'n_frozen,n_rel,n_raw_rel,'// &
        'maxch,maxch_app,maxatt,se_res,se_res2,trust,'// &
        't_iter,t_rad,t_rf,t_pops,rss_mb'
!
      WRITE(LU,'(I0,",",A,",",A,",",A,",",I0,",",I0,",",A,",",'// &
        '16(I0,","),6(ES12.4,","),4(F12.2,","),F12.1)') &
        ITER,TRIM(MODE_NOW),TRIM(MODE_NEXT),TRIM(BA_CH),F_I,H_I, &
        TRIM(REASON), &
        LEDG_N_FG,LEDG_N_MOM,LEDG_N_RF,LEDG_NP_RF, &
        N_BA_BUILD,N_BA_REUSE, &
        J_I,C_I,MAX(0,JEV_NOW-LEDG_JEV0),MAX(0,JJV_NOW-LEDG_JJV0), &
        R_I,LEDG_COUNT_TOKEN(REASON,'demote'),LEDG_N_FRZ, &
        N_FROZEN,N_REL,N_RAW_REL, &
        MAXCH,MAXCH_APP,MAXATT,SE_RES,SE_RES2,TRUST, &
        T_ITER,LEDG_T_RAD,LEDG_T_RF,T_POPS,RSS
      CLOSE(LU)
      RETURN
      END SUBROUTINE LEDG_WRITE_ITER
!
      END MODULE ITER_LEDGER_MOD
