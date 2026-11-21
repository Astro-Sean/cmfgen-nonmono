# KN_HOT_SPENT source audit (Phase 2.27F)

Status: source-reviewed. All line numbers vs working tree at
e28833fe + diagnostic diffs; `plane/fg_j_cmf_nonmono_v1.f`.

## 1. Lifecycle

| Variable | Decl | Type/scope | Init | Reset | Verified |
|---|---|---|---|---|---|
| KN_RAY_HOTPASSES | 626 | INTEGER, SAVE, module | =120 | SET_FG_HOTPASS (8350-58) | yes |
| KN_HOT_TOTAL_PASSES | 627 | INTEGER, SAVE, module | =480 | SET_FG_HOTPASS | yes |
| KN_HOT_SPENT | 628 | INTEGER*8, SAVE, module, SHARED | =0 | 8905 per refresh; 4217 refresh-skip clear; 2674 JFNK-eval pop | yes |
| W_LO/W_HI/NHOT | 10298 | INTEGER, host-scope in KNOP_RAY_SOLVE_V1 (thread-private via call stack) | 10472/10700/10800 per solve | per ray solve | yes |
| HOTMODE | 11617/12217 | LOGICAL, sweep-local | .FALSE. at sweep head | per M/P sweep | yes |
| L_HOT_OVER | 11678/12402 | LOGICAL, sweep-local | set in critical | per hot pass | yes |
| KN_HOT_CNT(NCF) | 1143 | INTEGER allocatable, module | 8999-9000 =0 | per refresh (alloc) | yes |

Reset at 8905 sits in the FG_KNOP_REFRESH_V1 prologue: one reset per
production refresh, before the parallel ray loop. Scope is therefore
per-refresh, refresh-global, shared across all threads/rays.

FG_EVAL_STATE_PUSH_V1/POP_V1 (2621, 2674) bracket KN_HOT_SPENT around
JFNK residual evals, so eval sweeps neither see nor consume the
production refresh budget.

## 2. Control flow

```
FG_KNOP_REFRESH_V1 (per refresh)
  KN_HOT_SPENT = 0                                    [8905]
  !$OMP PARALLEL DO NUM_THREADS(N_TH) SCHEDULE(DYNAMIC) [9129]
    KNOP_RAY_SOLVE_V1(LS)         ! one thread per ray at a time
      M sweep (SWEEP_FLD / SWEEP_FLD_ORD, INWARD=T)
        pass loop: on stall/cap exit -> TRY_HOT()     [11989/11998, 13191/13200]
          TRY_HOT: IF KN_HOT_SPENT >= CAP RETURN F    [15024, UNSYNC READ]
                   else W_LO/W_HI = hot run +/-2, HOTMODE=T
        hot pass top:
          NHOT=NHOT+1
          !$OMP CRITICAL (KN_HOT_SPD)
            KN_HOT_SPENT=KN_HOT_SPENT+1
            L_HOT_OVER = KN_HOT_SPENT > KN_HOT_TOTAL_PASSES
          !$OMP END CRITICAL
          IF NHOT>KN_RAY_HOTPASSES .OR. L_HOT_OVER -> IXIT=4, EXIT
          march restricted to W_LO..W_HI               [11697, 12487..12730]
      P sweep: same machinery, NHOT/W_LO/W_HI RESET at 10700/10800
      -> per-ray hot allowance is per SWEEP, not per ray
```

## 3. Semantics derived from executable statements

- The pair is a per-sweep cap (NHOT <= 120 per M sweep and per P
  sweep, since NHOT resets between M and P at 10700/10800) plus a
  refresh-global cap (total consumes <= 480, with overshoot below).
- Test-and-consume is ONE atomic operation inside CRITICAL
  (KN_HOT_SPD): the increment and the over-cap test cannot interleave.
- Overshoot accounting: the spend happens BEFORE the check, so a ray
  that pushes the counter over 480 is charged one unit and exits
  without marching. KN_HOT_SPENT can end a refresh > 480 by up to the
  number of concurrently hot sweeps. 9526 prints MIN(SPENT,2**30).
- Eligibility entry (TRY_HOT, 15024) reads KN_HOT_SPENT OUTSIDE the
  critical: advisory only. A stale read can admit a ray into HOTMODE;
  the next pass-top atomic consume then charges and (if over) exits
  it after zero restricted marches... no: the consume is at the TOP
  of the pass (12398), before the march, and on L_HOT_OVER it exits
  BEFORE marching. So a stale-admitted ray burns exactly one budget
  unit, runs zero extra passes. Bounded effect.
- IXIT=4 is the hot-exhaustion exit code (both sweep variants).
- Hot passes are excluded from KN_RES_PASS and RHO_WIN samples
  (11789, 11832, 13015, 13040): they do not enter the spectral-radius
  history or the residual-per-pass record.
- KN_HOT_CNT (10753-58, 10843-49) is a per-frequency histogram of
  rays ending above tolerance, updated under CRITICAL (KN_HOT_UPD);
  write-only, feeds only the summary report at 10158-10189. No solver
  feedback.
- v2 deck: KN_RAY_HOTPASSES and KN_HOT_TOTAL_PASSES are written in
  the deck header (11524). KN_HOT_SPENT, NHOT, W_LO/W_HI and the
  allocation order are NOT in the deck.

## 4. Consequence for replay

Per-ray executed hot-pass count is the only production quantity the
global counter influences: each sweep's hot run ends at
min(per-sweep 120 cap, global allocation remaining at its consume
times). Because consumes are atomic and the counter is monotone, a
per-ray recorded executed-hot count is SUFFICIENT to reproduce the
budget's effect; the global interleaving order need not be replayed.
```
