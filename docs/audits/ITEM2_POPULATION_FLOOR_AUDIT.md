# Item 2 - population floor handling audit (2026-11-09)

Status: audit in progress. Scope: independent numerical/software
correctness review of the Item-2 controls `POP_MASK_FLOOR` and
`POP_TRACE_MAXFAC` (commits `f4a0f6e`, `d8900ff`, `ce68692`, `5c59a53`),
plus the masked-demand accounting and LAMBDA-path coverage the stall
requires. Testbed: `cmfgen_validation_fast/m15_2p82d_KNOP_fast`
(H/He/Fe reduced model, KNOP, `interactive`+requeue).

This note is the running checklist required by the Item-2 task. Each
finding is marked [CONFIRMED], [SUSPECT] or [OPEN]; fixes and tests
are tracked at the bottom.

## 0. Headline result

Item-2 as shipped does **not** engage on the run it was built to fix.
The m15 fast KNOP stall executes entirely on `LAMBDA` calls
(`ba=LAMBDA` x25, `ba=MASKED` x2, zero non-LAMBDA BA calls in
`OUTGEN`), while both new controls sit inside
`IF(DO_BA_RES_ACC .AND. .NOT. LAMBDA_IT)`. They are dead code on this
run: `OUTGEN` it20-25 still shows the freeze/release limit cycle,
`exempted-demand`/`ex_maxch`/`pops-held`/`lambda-state` gates stuck at
50-63%, worst margin oscillating 3.2-4.2 dex. `obs_fin` not produced.

The floor-pinned channel is unchanged from the Item-1 census:
`FeI:3d6_4s2_a5De[2]` at depth ~63, `POPS ~ 1E-31`, `rcond ~ 1E-8`,
`colmin_u ~ 1E-43`, exit-exempted demand `BA_DEM_EXR ~ 97%`.

## 1. STEQ lifecycle (mapped)

`SOLVEBA_V13` computes `STEQ(NT,ND)` = demanded fractional correction
(`demand% = 100*|STEQ|`; `STEQ>0` decrease, `<0` increase). Order of
operations per call:

1. Band/DIAG solve fills `STEQ` (the ill-conditioned block solve;
   `rcond`/`colmin_u` measured here by the Item-1 census).
2. `STEQ_RAW = STEQ` snapshot (line ~1149) - pre-mask, for RUNAWAY.
3. `IF(DO_BA_RES_ACC .AND. .NOT. LAMBDA_IT)` block (line ~1154):
   a. `POP_MASK_FLOOR`: `STEQ(J,I)=0` where `POPS(J,I)<floor`,
      level rows `J=1..NT-2` only. [Item-2]
   b. `DO_BA_QUAR` parked-equation suppression.
   c. `POP_TRACE_MAXFAC`: `STEQ=SIGN(MAXFAC,STEQ)` where `|STEQ|>
      MAXFAC`. [Item-2]
   d. NaN-column zero, `BA_RBK_MASK` column zero, `BA_TRUST_CUR`
      all-or-nothing column zero (records `STEQ_MASKED_MAX`).
4. LAMBDA-only block (`DO_BA_RES_ACC .AND. LAMBDA_IT .AND.
   BA_PART_THIS`, line ~1314): RBK-masked columns zeroed. No floor
   mask, no trace cap.
5. LAMBDA `'LIMIT'` option and the LAMBDA-path quarantine
   (`BA_LQ_HIT` persistence -> `BA_QUAR_MASK`, line ~1375). The
   persistence scan explicitly `CYCLE`s rows with
   `POPS < BA_DEM_PFLOOR`, so floor-pinned rows are never parked on
   the LAMBDA path.
6. Demand scans (line ~1454+, EVERY call incl. LAMBDA): histogram,
   `DECREASE/INCREASE` -> `MAXCH`, `BA_DEM_EXM` (PFLOOR rows),
   `BA_DEM_EXR` (RELFRAC rows), `BA_EX_MAXCH`, `BA_TR_ABS`.
7. Caller `SOLVE_FOR_POPS` applies `STEQ` through scale limits
   (`MAX_LAM_COR` LAMBDA / `MAX_LIN_COR` BA) -> `SOL` -> `POPS=SOL`,
   subject to `GS_POPS_HELD` rollback and trust rejection.

## 2. Investigation A - POP_MASK_FLOOR defensibility

- [CONFIRMED] Post-solve, not matrix-level. The module comment
  (control_variable_mod.f ~1718) still describes the intended
  `CMF_DIAG_SOLVE` design - "equation row replaced by the identity,
  column decoupled, on the POPS-scaled matrix before DGEEQU" - but
  the shipped code zeroes `STEQ` after the solve. The comment itself
  records that the matrix pin "crashed the coupled block (job
  12868544)". **Stale documentation defect**: the described and the
  shipped algorithms differ.
- [CONFIRMED] Ill-conditioning is NOT removed. The near-null column
  (`colmin_u ~ 1E-43`, `rcond ~ 1E-8`) is still factorised; the
  solve still returns contaminated corrections for the *other*
  unknowns in the same block (they were computed against the
  floor level's coupling). Zeroing one `STEQ` entry post-solve makes
  the applied update a projection, not the Newton step. Defensible
  only because the masked correction was garbage anyway; it must be
  documented as such.
- [OPEN] Conservation/ion-balance: level rows `J=1..NT-2` exclude
  the ED/T rows (`J=NT-1,NT`) and conservation rows, so the element
  and ion sums are structurally untouched. Because `STEQ=0` freezes
  the level rather than redistributing its share, a masked level
  holds its floor population - negligible by construction. Needs a
  conservation check on a populated run (TEST 9).
- [CONFIRMED] One-way mask / no recovery. `STEQ=0` -> `SOL=POPS` ->
  `POPS` stays `< POP_MASK_FLOOR` -> re-masked next call, forever,
  on the BA path. A temporarily underpopulated level cannot recover
  through BA. (LAMBDA path is currently unmasked, which is the only
  reason any recovery is even possible - and it is also the stall.)
- [SUSPECT] Threshold semantics. `POP_MASK_FLOOR=1E-20` is absolute
  and sits 10 orders above the numerical floor (`~1E-30`) and 10
  orders above `BA_DEM_PFLOOR=1E-30`. It masks every level with
  `POPS<1E-20`, a broad trace range, not only floor-pinned ones.
  Whether `1E-20` is the right cut depends on the `POPS`
  normalisation; it is not tied to `BA_DEM_PFLOOR` or to the local
  species maximum. Risk of over-masking legitimately evolving trace
  levels -> they freeze -> convergence stalls for a different
  reason. (A-E)
- [CONFIRMED] Distinct from `BA_DEM_PFLOOR`. The two floors serve
  different purposes (demand-metric exemption vs commit masking)
  but both compare against the same `POPS`. Their 10-order gap is
  not documented as intentional.

## 3. Investigation B - LAMBDA vs BA pathway

- [CONFIRMED] The mask/cap never run on the stalled model: all 25
  recorded iterations are `ba=LAMBDA`. KNOP population
  preconditioning and the GS freeze/release nested iteration are
  LAMBDA calls, so the Item-2 controls are bypassed exactly where
  the stall lives.
- [CONFIRMED] The same floor level DOES generate runaway demand on
  LAMBDA: `BA_DEM_EXR ~ 97%` is read from unmasked `STEQ` on every
  LAMBDA call.
- [CONFIRMED] `GS_HOLD_POPS=T` freezes pops on `gs_freeze`
  iterations; combined with the unresolvable exempted demand the
  release/refreeze cycle repeats (it7-9 freeze, it10-11
  `pop_precond` release, it12-13 refreeze, ... it20-25 same).
- [CONFIRMED] The LAMBDA-path quarantine (`BA_LQ_HIT`) skips
  `POPS < BA_DEM_PFLOOR` rows by design, so the floor level is
  never parked - it churns unmasked indefinitely.
- [OPEN] A recovery strategy is required, not just LAMBDA masking:
  the task forbids naively extending the mask. See section 6.

## 4. Investigation C - convergence accounting (highest priority)

- [CONFIRMED] Masked demand vanishes silently. `STEQ=0` before the
  demand scans removes the floor level from `MAXCH`, `BA_DEM_EXM`,
  `BA_DEM_EXR`, `BA_EX_MAXCH` and the applied correction with no
  separate record. `STEQ_RAW` holds the pre-mask values but feeds
  only the `POP_AUDIT`-gated RUNAWAY census - it is not a
  convergence-facing channel.
- [CONFIRMED] The convergence exit CAN be satisfied by masking:
  if every problematic row is floor-masked, `EXM`/`EXR`/`MAXCH`
  collapse and `CONV_DEM_OK`/`BA_EX_MAXCH<EPS` can pass while the
  underlying SE equations are unresolved. This is the "suppress
  unresolved equations" failure the task forbids. No masked-demand
  count, largest-residual, or location is currently reported.
- [CONFIRMED] The hard gate is `BA_DEM_EXR` (exit-exempted,
  `BA_EX_RELFRAC` rows). `CONV_DEM_OK` gives the `EXM` leg a
  trace-polish escape (`CONV_TRC_NEED/REL`) but the `EXR` leg only
  passes when `EXR < max(10*EPS,1)` or `EXR >= BA_SAT_RAW`. The
  `a5De[2]` demand ~97% sits in the dead zone - this is the precise
  blocker.

## 5. Investigation D - POP_TRACE_MAXFAC

- [CONFIRMED] Dimension: `STEQ` is fractional, so `MAXFAC` caps
  `|STEQ|` (fraction). `POP_TRACE_MAXFAC=1E4` = 1E6 percent - it
  only touches corrections above 1E6%, far above the ~97% stall
  demand and far below the 1E7-1E13% runaway. As configured it is
  essentially inert for this stall; it would only catch the extreme
  tail.
- [CONFIRMED] Sign is preserved (`SIGN(MAXFAC,STEQ)`) - correct.
- [SUSPECT] Interaction with `BA_TRUST_CUR`: the cap runs before
  the trust test, so a column whose only over-trust row is capped
  below `BA_TRUST_CUR` escapes the all-or-nothing mask - the intent.
  But `MAXFAC=1E4` ~= `BA_TRUST=1E4`, so a row capped at the max
  can still sit at the trust boundary; the margin is not explicit.
- [OPEN] Whether repeated clipping creates a slow tail, and
  whether the cap should adapt to local population/iteration
  history, is untested. No diagnostics report capped-row count,
  largest raw vs capped correction, or subsequent column rejection.

## 6. Investigation E - GS freeze/release

- [CONFIRMED] The state machine responds to an unresolvable demand,
  not physical convergence: `gs_freeze_forced+hold_pops` engages,
  pops freeze, the exempted demand stays ~97%, release
  (`pop_precond`) re-demands, refreeze. `GS_REL_NOPROG_LIM=10` /
  `GS_FREEZE_LIM` / `N_FROZEN_TOT` bound it but the cycle is a
  limit cycle, not progress.
- [OPEN] A release must produce measurable progress; a failed
  release must name the blocking equation. Currently the blocking
  equation IS identifiable (`FeI:...a5De[2] @63` via `BA_DEM_EXR`)
  but the freeze/release diagnostics do not surface it as the
  cause.

## 7. Defects to fix (justified by findings)

D1. [CONFIRMED] Item-2 controls are LAMBDA-inert; the stall is
    all-LAMBDA. -> cover the LAMBDA path with a recovery-aware
    floor handler, not a blind mask extension.
D2. [CONFIRMED] No masked-demand accounting; masked residuals can
    silently satisfy the exit. -> add masked-demand channel
    (count, largest residual, species/level/depth) that stays
    visible and gates on negligibility.
D3. [CONFIRMED] One-way mask blocks recovery. -> explicit recovery
    eligibility for masked levels.
D4. [SUSPECT] `POP_MASK_FLOOR=1E-20` absolute threshold may
    over-mask; not tied to `BA_DEM_PFLOOR` or local scale.
    -> validate and document; consider relative/local criterion.
D5. [CONFIRMED] Stale module comment describes matrix-level
    masking; shipped code is post-solve. -> correct the comment.
D6. [OPEN] `POP_TRACE_MAXFAC` diagnostics absent; value may be
    mis-scaled for the stall. -> add capped-row diagnostics; review
    the cap semantics vs `BA_TRUST_CUR`.
D7. [OPEN] Freeze/release does not name the blocking equation on a
    failed release. -> surface the exempted/masked demand location
    in the freeze/release diagnostic.

## 8. Fix design (proposal, to be validated)

- Unified floor-level handler invoked on EVERY call (BA and
  LAMBDA), placed before the demand scans so `EXM`/`EXR`/`MAXCH`
  see consistent accounting.
- For `POPS(J,I) < POP_MASK_FLOOR`: record pre-mask `|STEQ|` into a
  masked-demand accumulator (count, max, argmax row/depth/species).
  Zero `STEQ` for the commit only when the demanded correction is
  unphysically large (ill-conditioned signature); a floor level
  with a bounded demand is left to evolve (recovery).
- Recovery: masked levels carry a per-(row,depth) streak counter;
  every `POP_MASK_PROBE` masked calls one unmasked update is
  permitted (bounded by the existing `MAX_LAM`/`MAX_LIN` clamps) to
  test whether the level can populate; success clears the mask
  state, failure re-masks.
- Masked demand is exported (`BA_MASK_N`, `BA_MASK_MAX`, location)
  and reported in OUTGEN; it does NOT fold into `EXM`/`EXR`/`MAXCH`,
  and the convergence exit requires masked populations to be
  provably negligible (they are `< POP_MASK_FLOOR` by construction)
  - reported loudly, never silently.
- `POP_TRACE_MAXFAC` gains capped-row diagnostics and an explicit
  relationship to `BA_TRUST_CUR`.

## 9. Tests / validation checklist

- [ ] TEST 1 disabled controls: `POP_MASK_FLOOR=0`,
      `POP_TRACE_MAXFAC=0` -> bit-identical stock path.
- [ ] TEST 2 below-floor negligible level handled without runaway.
- [ ] TEST 3 above-floor level never masked.
- [ ] TEST 4 below-floor level recovers when rates increase.
- [ ] TEST 5 masked residual stays visible; cannot silently satisfy
      exit.
- [ ] TEST 6 single runaway correction bounded; column mates
      progress; trust behaves.
- [ ] TEST 7 LAMBDA preconditioning does not stall/artificially
      converge.
- [ ] TEST 8 freeze/release makes measurable progress or names the
      blocker and terminates.
- [ ] TEST 9 conservation within existing tolerances.
- [ ] TEST 10 regression: JFNK/Knop/Sob/stock unaffected.
- [ ] `./build.sh` clean; debug bounds-check build.
- [ ] Fast-model A (Item-2 off) / B (current Item-2) / C (revised)
      iteration + metric comparison; physical consistency, not just
      `obs_fin`.

## 10. Implementation landed (09-Nov-2026)

Shipped in the working tree (build pending at time of writing):

- `POP_FLOOR_STEP` internal procedure in `SOLVEBA_V13`, called once
  per call after the `STEQ_RAW` snapshot and before the resacc
  trust/mask block - so the floor mask now runs on LAMBDA and BA
  calls alike (was `.NOT. LAMBDA_IT`-gated = dead on the stall).
- Masked-demand accounting: `BA_MASK_N` / `BA_MASK_MAX` /
  `BA_MASK_J` / `BA_MASK_I`, reset every call, reported in OUTGEN as
  `POP_MASK: suppressed N floor-level corrections; largest <eq> @dK
  dem%= X`. A masked equation can no longer vanish silently.
- Recovery probe: per-(row,depth) `BA_MASK_STREAK`; after
  `POP_MASK_PROBE` (VADAT, default 10) consecutive masked calls the
  row is left unmasked once so a temporarily underpopulated level
  can climb back above the floor. `POP_MASK_PROBE<=0` masks always.
- `BA_MASK_STREAK(NT,ND)` allocatable module state, sized like
  `BA_LQ_HIT`.
- Stale module comment corrected: it described the unimplemented
  matrix-level pin; now documents the shipped post-solve mask.

## 11. Fast-model validation result (job 12874555, warm from it31)

Baseline preserved as `OUTGEN_baseline_stall_oldcode` (old code,
it29-31 stall). Revised code engages correctly:

- `POP_MASK:` fires every call; `suppressed 2498 floor-level
  corrections`, largest `FeI:3d6_4s2_a5De[2] @d58 dem%= ~85`.
- `BA_DEM_EXM` (PFLOOR-exempted, POPS<1E-30) dropped to `0.0` - the
  deep-floor Fe I levels are now masked on LAMBDA calls.
- **But the stall persists.** `BA_DEM_EXR` is now driven by a
  *different*, populated-trace level: `Fe2:3d6(5D)4s_a4De[7/2] @d59`,
  demand ~35-56% and growing. Gates still 50-63%, `exempted-demand`/
  `ex_maxch`/`pops-held`/`lambda-state` failing, margin ~3.3-3.8 dex.

Two conclusions:

1. **POP_MASK_FLOOR is necessary but not sufficient.** Masking
   floor-pinned levels clears the `EXM` channel, but `EXR`
   (exit-exempted, `BA_EX_RELFRAC`-relative trace rows) is driven by
   *populated* trace levels (`POPS >= POP_MASK_FLOOR`) whose demand
   is real. Floor masking cannot and should not touch those.
2. **The remaining blocker is Stage-B trace polish.** `EXR` only
   resolves via the trace-polish sweeps (`BA_TR_MAXSWP`), which arm
   only after the *major* channels (`ex_maxch`, `lambda-state`,
   `pops-held`) converge - and those are still failing in the
   freeze/release cycle. The stall is the broader KNOP
   trace-population convergence problem, not floor-pinning alone.

## 12. Outstanding / next

- [OPEN] `POP_MASK_FLOOR=1E-20` masked **2498 levels/call** - likely
      over-broad. The threshold is absolute and not tied to
      `BA_DEM_PFLOOR` or the local species maximum; it freezes a wide
      trace range. Validate against the `POPS` distribution; a value
      nearer the numerical floor (~1E-28..1E-30) or a relative
      criterion is probably correct. Over-masking risks an
      unphysical `obs_fin` (frozen trace populations).
- [OPEN] `POP_TRACE_MAXFAC` diagnostics (capped-row count, largest
      raw vs capped, subsequent column rejection) not yet added; the
      cap remains BA-path-only and inert at 1E4 for the ~1-fraction
      stall demand.
- [OPEN] GS freeze/release does not yet name the blocking equation
      on a failed release (D7).
- [OPEN] Stage-B/`EXR` engagement is the real remaining blocker -
      likely needs the major channels to mature first, or a
      trace-polish trigger that does not wait on full major
      convergence. Separate item.
- [ ] Unit tests TEST 1-10 (solver-path coverage limited; see
      report).
- [ ] `./build.sh` clean + `run_all_tests.py` regression.
- [ ] Commit on a clearly-identified branch/commits.
