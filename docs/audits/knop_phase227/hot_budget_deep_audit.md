# Phase 2.27G - KN_HOT_SPENT deep audit

Scope: `plane/fg_j_cmf_nonmono_v1.f` at e28833fe + local diagnostic
diffs. All line numbers verified against the current checkout.
ORD-sweep lines cited; SWEEP_FLD equivalents noted where they differ.

## 1. Executive verdict

- CONFIRMED FROM EXISTING OUTPUT: the global budget never bound in
  the archived production run. Every refresh reports
  `hot-pass budget spent=` 11-29 of 480 (OUTGEN lines 1120..8493,
  15 refreshes). Cross-ray budget contention is ELIMINATED as a
  cause of the ray-81/82 replay discrepancy.
- CONFIRMED FROM SOURCE: the consume is atomic
  (CRITICAL KN_HOT_SPD); TRY_HOT's pre-entry read is unsynchronised
  (class B). The stale-admission sequence exists in source but
  requires SPENT to reach 480, which never occurred.
- CONFIRMED FROM SOURCE: hot passes are NOT controller-neutral.
  RATIO_OM, OMEGA_RAY, NOSC, NSTALL/NPLAT, BEST_*, RES_PREV all
  update on window-restricted residuals during hot mode.
- CONFIRMED FROM SOURCE: IXIT=4 folds three distinct exits
  (itlim-denial, per-sweep 120 cap, global-cap rejection) into
  KN_CS_ITLIM. The census cannot distinguish them.
- CONFIRMED FROM SOURCE: a rejected consume exits before march,
  before omega injection (kn_replay_apply), and before any KNOMG
  row. Replay schedules contain executed passes only - correct for
  single-ray replay, incomplete for whole-refresh reconstruction.

## 2. Control-flow map (SWEEP_FLD_ORD)

| Step | Lines | Sync | State |
|---|---|---|---|
| pass top: NPO++, OMG_USED/NOP snap | 12384-12392 | thread | local |
| HOTMODE: NHOT++ | 12398-12399 | thread | local |
| consume + over test | 12400-12403 | CRITICAL KN_HOT_SPD | shared write, local verdict |
| reject: IXIT=4 EXIT | 12404-12407 | thread | local |
| omega inject (replay) | 12409-12415 | thread | local - never reached on reject |
| march (window CYCLE guards) | 12487,12516,12641,12730 | thread | local |
| BEST_* update | 13027-13031 | thread | local, runs in hot |
| RATIO_OM=RES/RES_PREV | 13032-13034 | thread | runs in hot |
| RHO_WIN sample | 13040-13044 | thread | GATED .NOT.HOTMODE |
| NOSC/KN_SGN_PREV | 13053-13059 | thread | runs in hot |
| omega controller | 13067-13100 | thread | runs in hot |
| NSTALL/NPLAT/L_PLAT | 13101-13127 | thread | runs in hot |
| RES_PREV(1)=RES | 13130 | thread | runs in hot |
| KN_OMG_DUMP row | 13131-13134 | thread file write | every pass incl. hot |
| conv exit IXIT=1 | 13161-13165 | | |
| bulk exit IXIT=2 | 13177-13184 | | |
| stall: TRY_HOT else IXIT=3/6 | 13190-13198 | | unsync read inside TRY_HOT |
| cap: TRY_HOT else IXIT=4 | 13199-13206 | | same |
| TRY_HOT body | 15019-15038 | | 15024 unsync shared read |

SWEEP_FLD duplicates the structure: consume 11674-11684, entry
points 11988-12004 (same IXIT semantics); its RES_PREV/RATIO_OM are
same-direction indexed (RES_PREV(IDIR)) at 11810-11912, KN_OMG_DUMP
at 11912.

## 3. Race analysis

- CONFIRMED: consume+verdict is one critical unit; no thread can
  steal the last unit undetected.
- CONFIRMED: TRY_HOT reads KN_HOT_SPENT at 15024 outside any sync.
  Sequence A-D (stale admit -> consume -> L_HOT_OVER -> IXIT=4 with
  zero marches) is source-possible; it requires the counter to
  actually reach 480 within the interleave window.
- CONFIRMED: TRY_HOT can deny at >=480 WITHOUT consuming (early
  return before the window scan). A denial exits with the
  UNDERLYING code: IXIT=3/6 at the stall path, IXIT=4 at the cap
  path - never a distinct "budget-denied" code.
- CONFIRMED: overshoot accounting - the increment precedes the
  test, so rejected consumes still burn a unit; SPENT can exceed
  480 by up to the number of concurrently-consuming sweeps.
- PLAUSIBLE BUT UNTESTED (moot here): under contention, per-sweep
  allowance is schedule-dependent (SCHEDULE(DYNAMIC), 9129).

## 4. Rejected-attempt side effects (IXIT=4 consume rejection)

CONFIRMED FROM SOURCE, in order: NPO incremented; NHOT incremented
(persists in counters); KN_HOT_SPENT incremented (unit burned);
IXIT=4; sweep EXIT. NOT changed: FLD (no march), RES/RES_FREQ
(last pass's values retained), RES_PREV/RATIO_OM/NOSC/NSTALL/
OMEGA_RAY (pass-end block never runs), AA/DLT histories. OMG_USED
was snapshotted but no KNOMG row is written (dump is pass-end).
kn_replay_apply is never reached, so the rejected pass is absent
from the armed-ray schedule. W_LO/W_HI/HOTMODE persist but are
dead (sweep exited; reset next sweep anyway).

## 5. Certification mapping (KN_CERT_CH 2400-2456, KN_CERTIFY 2469-2499)

CONFIRMED: per-channel classes - IXIT 4 -> ITLIM unconditionally;
IXIT 3 or 6 -> STAG; IXIT 2 + BULK floor -> BULK else INV(5);
IXIT 5 -> DIV; IXIT<0 or >6 -> INV(5); RES nonfinite/negative ->
INV(3); INVR -> INV(reason). Ray class = worse of M,P (ISM>=ISP).
Consequence: a global-cap rejection is indistinguishable from a
plain iteration-limit exit in the census. STAG<->ITLIM flips
between production and replay are possible ONLY in a contended
run (never here). There is no hot-exhaustion class.

## 6. Controller state during HOTMODE (ORD; SWEEP_FLD analogous)

| Variable | Ordinary pass | Executed hot pass | Rejected consume | M/P boundary |
|---|---|---|---|---|
| NPO | +1 | +1 | +1 | NP_M/NP_P exported |
| NHOT | - | +1 | +1 | reset 10700/10800 |
| RES_PREV | =RES | =RES (window res) | unchanged | reset 12308-09 |
| RATIO_OM | upd | upd | unchanged | reset 12310 |
| NOSC/KN_SGN_PREV | upd | upd | unchanged | reset 12299/head |
| NSTALL/NPLAT/L_PLAT | upd | upd | unchanged | reset at head |
| OMEGA_RAY | controller | controller on window res | unchanged | reseeded KN_RAY_OMEGA |
| BEST_* | upd | upd (window res) | unchanged | reset 12318 |
| RHO_WIN/KN_RES_PASS | sampled | EXCLUDED | - | - |
| FLD | marched | marched in window | unchanged | - |

CONFIRMED: hot passes perturb the controller (omega adapts to
window residuals; RATIO chains cross the full->window transition).
A rejected consume is controller-neutral but not counter-neutral.
PLAUSIBLE: window/full RES mixing in RATIO_OM can shift later
stall timing by a pass - bounded, unmeasured.

## 7. KNAPPLY/KNOMG coverage

CONFIRMED: KN_OMG_DUMP writes one row per EXECUTED pass (armed
rays), including hot marches, with OMG_USED/OMEGA_RAY/RES/
RATIO_OM/NOSC. Missing events: TRY_HOT evaluations and denials,
window bounds, consume verdicts, rejected consumes, HOTMODE
entry/exit flags, per-sweep NHOT total. Sufficiency: adequate for
single-ray Mode B (executed schedule is the contract; a replay
asking for a pass beyond the record fails loudly). NOT sufficient
for whole-refresh mechanism replay.

## 8. Existing OUTGEN evidence

CONFIRMED FROM EXISTING OUTPUT: 15 refreshes, spent 11-29/480.
Print at 9526-9528 (FG_KNOP_REFRESH_V1 epilogue, post-region,
before next refresh's reset at 8905); JFNK evals restore the
counter (2621/2674) so the printed value is the production
refresh total. Monotone within refresh -> interpretation:
cap never bound; TRY_HOT never denied on the global test;
no rejected global consume possible. The entire shared-budget
channel is dead for this run. "hot frequencies" aggregate lines
(~40-50k cells in 700-1250 runs) show hot-mode candidates are
common; spends are small, so hot runs are short.

## 9. Minimum sufficient replay record

Goal A (single-ray prescribed schedule - current Mode B):
CONFIRMED sufficient as-is: executed-pass omega schedule +
per-sweep NHOT<=120 enforced identically + fresh budget that
cannot reject. Recommended (not required) addition: HOTMODE flag
per KNOMG row so hot vs ordinary passes are explicit rather than
inferred from position.
Goal B (whole-refresh mechanism replay) would require per event:
refresh, LS, sweep, NPO, NHOT, thread, event-type (TRY_HOT
eval/deny, consume verdict, reject, march), SPENT before/after,
W_LO/W_HI, final IXIT - the ordered consume log. NOT warranted:
the mechanism never bound.

## 10. Experiment design (proposed only)

If ever needed: identical exe/input hashes, fixed node, N repeats
at production threads (>=5) + serial baseline (>=3) + optional
2/4 threads; compare final SPENT, per-sweep hot counts, IXIT,
KN_CS, FLD/store bits. Uses diagnostic binary for telemetry; an
uninstrumented binary run alongside tests instrumentation
perturbation. PRIORITY: LOW - archived data already shows no
contention for this run; revisit only if a future refresh
approaches 480.

## 11. Revised hypothesis ranking

1. M-state reconstruction - leading (unverified runtime gap;
   stall-region sensitivity per the ~1.0 verify defect).
2. Omega/controller application fidelity - Mode B bounds it.
3. Hot budget GLOBAL contention - ELIMINATED (existing output).
   Residual: per-sweep 120 cap still shapes production pass
   counts but deterministically (thread-independent).
4. FP noise - bounded by ray-84 (3.5e-8); implausible at 1e4x.

## 12. Recommendation for Phase 2.27H

GO (narrow): the next discriminating datum is replay-IM_NI vs
deck-IM_NI plus census IXIT_M/IXIT_P/NP for rays 81/82/84 from
archived KNCENS/census rows (read-only). No telemetry build is
justified by the budget evidence. Phase 2.27 remains PARTIAL
until runtime comparisons run under an approved allocation.
