# Field maturity and the GS-freeze controller

`SOLVE_FOR_POPS` evaluates the just-completed sweep's lagged-field
maturity at the top of every call and runs the mode controller. This
is the enforcement of the central ordering rule (Knop Paper II,
Sect. 3.2):

> [!important] The ordering rule
> The populations must never take an update against a radiation field
> that is not the fixed point of the transfer map at the current
> opacities.

Analogy: weighing yourself on a scale that is still rocking. The
field-maturity metrics measure whether the scale has settled;
`GS_REL_STABLE` demands it reads the same value twice in a row before
you trust it; `GS_FREEZE` stops you recording a weight (an SE update)
while the scale is still settling. A single still reading proves
nothing - the rocking pauses mid-swing too.

The two failure modes it prevents are mirror images:

1. An SE update on a still-drifting field re-poisons the field the
   sweep just converged (it is iterated on data the populations then
   move).
2. A converged field spent on an SE step that explodes forces the
   whole outer map to re-contract.

## Maturity observables

Fed by `GET_FG_GS_MATURITY` / `GET_MOM_GS_MATURITY` each sweep:

| Metric | Meaning |
|---|---|
| `GS_MAXCHG` | max normalized FG refresh residual; in-sweep store change is retained only when refresh did not run |
| `GS_NBAD` | significant cells whose same normalized refresh residual exceeds `KN_RAY_TOL` |
| `GS_NCONV` | significant stored cells whose fresh value differs by >50% (store churn) |
| `GS_CAPPED` | writes dropped as non-finite/exploding |
| `MOM_SW_MAXCHG` | max drift, fresh J solve vs stored `J_NEXT` |
| `MOM_SW_NBAD5/50` | cells drifting >5% / >50% |
| `MOM_SW_NHUGE` | cells drifting >1e6 (collapse detector) |
| `MOM_SW_NEGJ / NQUAR / NPQUAR / BNDFL / CAPPED` | repairs, quarantine (+persistent), boundary floor, dropped writes |
| `MOM_RF_CONV / DCHG / N5 / N50` | refresh converged? retained fixed-point delta |

Thresholds (`solve_for_pops.f`):

- `GS_MATURE_TOL = 0.05` - a sweep whose store moved <5% at the max
  cell is a maturity *candidate*.
- `GS_UNSAFE_TOL = 0.5` - max change >50% is runaway regardless of
  bulk.
- The FG residual is `abs(F(I)-I)/max(I, 1e-4*column_peak)`. `GS_NBAD`
  counts cells using that exact residual threshold; `GS_NTOT` counts
  cells above the 1e-4 significance floor. This keeps the residual max,
  per-ray bulk exit and population gate on the same tolerance measure.
- Bulk rule: mature when fewer than ~0.1% of **significant** cells are
  above tolerance or have >50% store churn. The two counts are checked
  separately. Optically thin tail cells below 1e-4 of the column peak
  do not dominate either fraction.
- MOM two-tier gate: **unsafe** tier (veto on everything: >10% capped,
  >1% drifting >50%, >0.01% collapsed, one non-finite, >10%
  negJ/quar) vs **marginal** tier (blocks release, eligible for
  `MOM_TRIAL`).
- `KN_GATE_DIVERGENCE` (default F, prototype, dormant): replaces the
  absolute `GS_MAXCHG > GS_UNSAFE_TOL` unsafe tier with a divergence
  test - unsafe only when the refresh max residual is non-finite or
  grows by more than `KN_GATE_DIV_FAC` (default 1.5) against the
  value `KN_GATE_DIV_M` (default 4) refresh sweeps back (ring held
  in `SOLVE_FOR_POPS`, SAVE-persisted). A bounded stationary
  residual like the 0.53 period-2 asymptote on m15_2p82d then counts
  as mature instead of runaway. There is deliberately no absolute
  backstop in this variant: a large but stable residual passes.
  Enabling is a gate policy change and needs owner sign-off; the
  p2damp validation arm leaves it F.

## The controller

```text
GS_RAW_IMMATURE = capped>0 OR FG-immature OR MOM-freeze-b
GS_FREEZE       = GS_RAW_IMMATURE OR (N_MATURE < GS_REL_STABLE)
```

- `N_MATURE` counts **consecutive** raw-mature calls; one churning
  sweep resets it. `GS_REL_STABLE` (default 1; 2 in the hybrid config)
  is the stability hold: a single mature-looking sweep can be a
  coincidence of a churning field - the premature release it enables
  re-poisons the field (the documented freeze/release limit cycle).
- During `GS_FREEZE`: `LAMBDA_ITERATION=T`, `FIXED_T=T`,
  `COMPUTE_BA=T`, `TEMP_CHAR='DIAG'` - SOLVEBA still runs (its
  demanded-correction diagnostics keep flowing) but takes a bounded
  Lambda update.
- **Churn-scaled Lambda cap**: when bulk churn
  `max(GS_NCONV, MOM_SW_NBAD50)/NTOT < GS_LAM_TRANSIENT` (2%) the field
  is near-mature and the correction is capped at `GS_LAM_CAP` (factor
  3); in gross transient it keeps `MAX_LAM_COR` (production profile
  `1e10`, effectively unclamped). A large committed step on a
  near-mature field re-poisons the store for many iterations.
- When `BA_LUM_GATE=T` and the luminosity ratio drifts beyond
  `BA_LUM_DCHG`, a scheduled BA step is demoted before solve and the
  LAMBDA correction is rolled back. Transfer then advances at fixed
  populations until the measured ratio is stable; BA_RES_ACC does not
  assess ordinary LAMBDA corrections. Near-mature release corrections
  are capped at factor 1.1. This protection is off with the default
  `BA_LUM_GATE=F`.
- `GS_HOLD_POPS` (default F): snapshot POPS before the forced-Lambda
  solve and restore after - the field converges at *literally fixed*
  populations (strict Paper-II nested order). If the snapshot cannot
  be allocated, the solver stops before updating populations rather
  than silently proceeding on the immature field. `GS_POPS_HELD` /
  `GS_ATTEMPT_APP` record either hold.
- The strict GS field-before-population ordering requires
  `DO_LAM_IT=T`. With `DO_LAM_IT=F`, the stock-cycle branch preserves a
  scheduled BA step even while `GS_FREEZE` is true; that mode is not
  safe for this ordering.

## Warmup and restart

- `N_LAM_RWARM` (default 4): the `I_NEXT` store is process-local and
  re-zeroed on every requeue; the first N calls of every process are
  forced Lambda. A BA step against a zeroed store is the classic
  requeue wreck.
- `GS_RWARM_SKIP`: a valid `FGINEXT_GS` checkpoint with saved max
  residual <5% and no capped writes may skip the process warmup. The
  checkpoint does not store the bulk census; the current iteration's
  refresh recomputes it before `SOLVE_FOR_POPS` can take a live BA step,
  so this optimization does not bypass the in-run population gate.
- `ITS_SINCE_START` (not checkpointed `MAIN_COUNTER`) drives the
  warmup count, matching the process-local store lifetime.

## Marginal release: `MOM_TRIAL`

A field frozen only on the marginal tier (no unsafe events) can sit in
a stable asymptote the Lambda map cannot improve - the FreeNe run
froze >150 iterations on ~0.12% drifting cells. The trial fires only
when **all** hold:

```text
MOM_RAN & !MOM_UNSAFE_B & !MOM_BROAD & DO_BA_RES_ACC
  & !BA_SNAPSHOT_MISSING & N_FROZEN >= MOM_TRIAL_MIN
  & GS_CAPPED==0 & !GS_FG_IMMATURE & !LUM_BLOCKED
  & N_TRIAL_ACC < GS_TRIAL_LIM
```

- The trial takes a real BA step (`T1=MAX_LIN_COR`) under the resacc
  rollback net.
- An **accepted** trial refunds the freeze clock
  (`N_FROZEN -> MOM_TRIAL_MIN`, not zero - the post-release field is
  re-immature by construction so the next trial must re-earn marginal
  status) and increments `N_TRIAL_ACC` (persisted via `GS_STATE`).
- `GS_TRIAL_LIM` (default 15) bounds accepted trials between
  raw-mature sweeps; a genuine 2-sweep mature release resets it.
- `N_FROZEN` keeps counting: `GS_FREEZE_LIM`/`KNOP_STALL_MAX` still
  kill a run that never matures.
- **JFNK never fires on a trial** - an unbounded Newton commit would
  defeat the trial's purpose.

## Watchdogs

- `GS_FREEZE_LIM` (default 0=off): in-model hard stop after N
  consecutive frozen iterations.
- `KNOP_STALL_MAX` (default 60): the slurm-wrapper stall bound, reads
  `GS_STATE` written each call.
- Plateau detector: 6-iteration ring of applied correction, second SE
  channel, drift fraction; a <1.5x window flags `CONV_PLATEAU`.

## Related

- [[Knop formal solver]] - where the field and its metrics come from
- [[JFNK admission]] - GS_FREEZE is the field gate for Newton
- [[Failure modes and guards]]
- [[Diagnostics and OUTGEN]] - the ITER_SUM fields
