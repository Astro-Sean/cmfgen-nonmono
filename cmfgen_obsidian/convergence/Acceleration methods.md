# Acceleration methods

Techniques that speed up the outer iteration without changing the
fixed point: Ng acceleration, oscillation averaging, frozen-Jacobian
reuse, and the offline equivalents.

## Ng acceleration

`subs/ngaccel.f` / `ngaccel_arb_ord.f` implement Ng's extrapolation
(Auer, *Numerical Radiative Transfer*, p. 101): a weighted Aitken-like
acceleration valid when the iteration operator is **linearly
convergent** - which is what Lambda iteration is. It is an outer-map
accelerator, applied to the population iterates (and, on the stock
path, to R, V, SIGMA together via `DO_NG_BAND_ACCEL_V4`). Under
`DO_NONMONO_SOB`/`DO_NONMONO_KNOP` the structure rows are frozen at
their current prescribed values before the fit, so NG can accelerate
populations only.

Scheduling (`solve_for_pops.f`, `NEXT_NG` state machine):

| VADAT key | Variable | Role |
|---|---|---|
| `DO_NG` | `NG_DO` | master switch |
| `BEG_NG` | `VAL_DO_NG` | MAXCH (%) below which NG arms |
| `IBEG_NG` | `IT_TO_BEG_NG` | earliest iteration |
| `BW_NG` | `NG_BAND_WIDTH` | depth-band width for the extrapolation |
| `ITS/NG` | `ITS_PER_NG` | spacing between accelerations (>=4) |
| `CHK_NG` | `DO_NG_VALIDITY_CHECK` | sanity-check the T correction |

- `NEXT_NG` sentinels: 1000 = new model, 1500 = MAXCH back above
  `VAL_DO_NG` (disarmed), 2000 = continuing model. Once armed, the
  accel fires at `NEXT_NG`, >= `LAST_NG+ITS_PER_NG` and >= 4
  iterations out (the method consumes the last 4 iterates, stored via
  the `SCRTEMP` history).
- A successful accel is a **kick**, not a normal step:
  `MAIN_COUNTER` increments, `LAST_NG` is checkpointed, and if the
  corrections were large (T1,T2 > 50) the next iteration is forced
  LAMBDA (`ITER_REASON='ng_kick'`) - exactly the re-poison the freeze
  controller exists to prevent.
- A failed accel retries in 8 iterations.
- `LAST_NG` in `POINT1` can be set <= `MAIN_COUNTER-ITS_PER_NG` to
  force an accel; `IT_SPECIFIER`'s `DO_NG` directive sets
  `NEXT_NG=MAIN_COUNTER+1`.
- **Guards added in this tree**: NG is blocked during `GS_FREEZE` and
  `FORCE_LAMBDA` (an extrapolation across an `I_NEXT` reinit commits a
  population kick against a stale field). Under non-mono,
  `DO_NG_BAND_ACCEL_V4` overwrites the four structure-history records
  with the current prescribed `R`, `V`, and `SIGMA+1`, preventing a
  stale or extrapolated grid/velocity field from being committed.
  `RECOMPUTE_SIGMA_TRUE` is still called after NG so the derived field
  remains synchronized with any future allowed structure update.
- The convergence exit is suppressed on the accel iteration itself
  (`MAIN_COUNTER .NE. LAST_NG` required to set `NUM_ITS_TO_DO=1`), and
  at least one post-accel iteration is enforced
  (`NUM_ITS_TO_DO=MAX(2,...)` at `LAST_NG`).

Acceptance is the weak point on the non-monotonic SOB model family
(05-Oct-2026 audit of `dessart_m15_2p82d_SOB_mask`, job 12084677):
mask-1 scheduled NG 18 times with every attempt cancelled
(`NUM_BAD_NG` ~50-56 of 70 depths, mostly `LOCDEC<0.09`), while v2
managed 14 successes against 39 cancels.  Two distinct causes:

- **Polluted history.** `SCR_RITE_V2` appends a record every call
  before the residual judgment that may veto it, so a rejected step's
  record is on disk and the rejecting call appends a restored
  duplicate - sign-reversing difference pairs in the 4-record window.
  Masked rollback (`BA_REJ_FRAC`) leaves half-restored records.
  `UNDO_IT`/`AVE_FLIPS` index `SCRTEMP` by iteration number, so
  records cannot simply be skipped.  The scheduler-level fix is
  `NG_CLEAN` (module-local, not a VADAT key): counts consecutive
  calls since the last reject / masked-rollback / soft-reject /
  `BAACC_PENDING` / JFNK / held-population record; the fire test
  defers by one call while `NG_CLEAN<4`. `NG_DEFERRED` marks a +1
  target created by this deferral, so an auxiliary SCRTEMP write that
  advances `MAIN_COUNTER` cannot step over and permanently miss it.
  Only such self-deferred targets use a less-than-or-equal due check;
  other `NEXT_NG` scheduling keeps the equality test. With resacc,
  JFNK, and held-population paths inactive, the upstream schedule is
  unchanged.
- **All-or-nothing per-depth veto.** The banded fit produces one
  secant per ~40-depth band shared by ~1842 components (levels + R +
  V + SIGMA+1 + T); one bad component vetoed the whole depth, and >3
  bad depths cancelled the run.  `NG_ABSTAIN=T` (new VADAT key,
  default F) instead holds per-component: zero-newest-delta
  components (restored/static rows) and non-positive extrapolations
  keep their current value and report ratio 1.  A vacuous attempt
  (every component held) returns `NG_DONE=F` so no identical record
  is appended.

Pipeline default is `DO_NG=F` since 13-Oct-2026 (X_SOB 40d): accepted
kicks overshot the interaction-depth populations by up to ~30% every
`ITS/NG` iterations, so MAXCH oscillated near ~0.5% and never crossed
`EPS_TERM` after dead-level demand was gated by `BA_DEM_RELFRAC`.  The
unaccelerated map contracted ~0.65x/iteration.  The NG machinery stays
armed for models that opt in via `vadat_overrides DO_NG=T`: `CHK_NG=T`,
`BEG_NG=5.0D0`, `IBEG_NG=15`, `BW_NG=40`, `ITS/NG=8`, and for non-mono
runs `NG_ABSTAIN=T`.  When NG is enabled, `NG_CLEAN>=4`, `GS_FREEZE`,
`FORCE_LAMBDA`, and residual-acceptance history still decide whether a
scheduled attempt is allowed to fire.

### Adaptive-gate assessment (16-Oct-2026; deferred)

Status: **future candidate update only**.  The production scheduler
remains the fixed clean-history cadence (`BEG_NG=5`, `IBEG_NG=15`,
`ITS/NG=8`); do not add adaptive firing until a side-effect-free
candidate path has been implemented and model-tested.

Magritte's adaptive Ng criterion compares the maximum relative change
between **successive Ng predictions** with the maximum relative change
between the **successive ordinary iterates**.  It accepts the latest
prediction only when the running predictions contract faster than the
regular trajectory (Ceulemans et al. 2024, arXiv:2411.03040, Sec. 4.2).

This is viable as a dry-run predictor in this tree, but it is not a
one-line scheduler edit.  `NG_MIT_BAND_OPT_V1` already builds the full
candidate in `NEWPOP` before committing it, so the candidate can be
returned without mutating the live state.  A viable fixed-order gate
would:

1. evaluate a candidate on every eligible clean-history call,
   not just on `NEXT_NG`;
2. save the previous candidate in process-local storage (a restart
   simply loses one comparison call);
3. compute `d_ng=max_rel(pred_k,pred_{k-1})` and
   `d_reg=max_rel(x_k,x_{k-1})` over positive components;
4. commit only if `d_ng < d_reg` and the existing `NG_ABSTAIN`,
   per-depth, `CHK_NG`, and post-kick checks pass;
5. reset the saved prediction after an accepted kick or a mode change.

The hard parts are the same-map condition and residual truth.  The
4-record window must not span a LAMBDA/full or solver-mode transition
(`NG_CLEAN` should reset on `MODE_NOW` changes as well as rejected
records), and a cheap prediction metric is not a true residual
acceptance test.  Evaluating `SE(candidate)` without committing would
cost a full nested solve and would need snapshots of the GS/BA state;
that is not viable as a low-cost pre-check.  Variable-order Ng
(Magritte's growing `N`) is also a larger change: `NGACCEL` and
`NG_MIT_BAND_OPT_V1` are hard-wired to four records.

Decision (16-Oct-2026): defer implementation.  Revisit only after the
fixed `BEG_NG=5` schedule has been exercised on a clean non-mono run;
the first candidate implementation should be diagnostics-only before
allowing adaptive commits.

Offline equivalent: `misc/do_ng_v2.f` runs NG / AV / SOR acceleration
on the `SCRTEMP` iterate history without running the model.

## Oscillation averaging (`DO_AV`)

`AVE_FLIPS` averages populations that oscillate iteration-to-iteration
instead of contracting. VADAT: `DO_AV`, `NOSC_AV` (consecutive
oscillations to detect, default 4), `ITS/AV` (spacing, default 8).
Deliberately kept away from other events: requires
`MAIN_COUNTER > LAST_NG+NUM_OSC_AV`, `> LAST_LAMBDA+NUM_OSC_AV+2`,
not within `NEXT_NG-4`, and not during `GS_FREEZE`. On success it
writes a checkpoint and increments `MAIN_COUNTER`; on failure it
retries in 8.

## Frozen-Jacobian reuse (`FIX_BA` / `N_FIX_BA`)

Below `FIX_BA`% demand the BA matrix can be frozen and reused for
`N_FIX_BA` iterations instead of rebuilt - a significant cost saving,
since `GENERATE_FULL_MATRIX_V3` is a substantial fraction of a full
iteration. Corrections then go through the `ORIG_POPS` conversion,
rebounded by the `MAXCH_SUM > 3*FIX_BA` recompute trigger. The run
history includes a 175-iteration stale-Jacobian plateau at
`FIX_BA=15, N_FIX_BA=0`; recommended `FIX_BA=1.0, N_FIX_BA=10`
(periodic recompute). `IT_SPECIFIER`'s `FULL(FIXBA)` directive forces
reuse for a scheduled stretch.

## Newton as accelerator

JFNK is the endgame accelerator of the outer SE map - see
[[JFNK solver]] and [[JFNK admission]]. The distinction to keep in
mind: Ng and AV assume nothing about the operator beyond smoothness;
Newton assumes local linearity and pays for it in eval cost. That is
why Ng/AV run mid-iteration and Newton is last.

## Deliberately not used

- Anderson acceleration of the outer `I_NEXT` store-update map:
  considered and deferred - it needs a contractive inner map first
  ([[JFNK solver]], "Deliberately not implemented").
- `DO_UNDO` (`UNDO_IT`, undo corrections at the last 5 depths) exists
  but is off by default.

## Related

- [[Lambda iteration]] - the operator Ng accelerates
- [[Convergence strategy]] - the staging order
- [[Warmup and staged release]] - `IT_SPECIFIER`/`IN_ITS` scheduling
- [[Convergence criteria]] - `LAST_NG` in the exit test
