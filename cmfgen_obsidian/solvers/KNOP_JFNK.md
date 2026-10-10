# KNOP + JFNK: the combined non-monotonic solver

This note describes how the Knop non-monotonic formal solver, the
Gauss-Seidel field stabilization machinery, and the JFNK population
solver fit together in this tree, and which stability criteria decide
when each piece may run.  It is the architecture reference for the
three-phase convergence structure: Lambda/GS first, frozen-J JFNK
second, coupled JFNK only at the end.  All of it is gated on
`DO_NONMONO_KNOP` / `USE_JFNK`; the default path is unchanged.

Related documents:

- [[KNOP_CONVERGENCE_STRATEGY]] -- mapping onto the Knop
  papers and the ordered-march refresh (`SWEEP_FLD_ORD`).
- [[JFNK_SOLVER_AUDIT]] -- cost asymmetry and the JVP noise
  floor that motivate the tiered admission.
- [[JFNK_SOLVER]] -- the Newton formulation (log-x unknowns,
  pseudo-transient continuation, forcing floor).
- [[convergence_issues]] -- the run log of failure modes each guard
  was built against.

## 1. The problem and the three nested maps

A non-monotonic velocity field turns the comoving-frame transfer
equation from an initial-value problem into a boundary-value problem
in wavelength: at a compression front the comoving frequency can
decrease along a ray, so intensity is transported *up* in wavelength
and the discretization couples each cell to a not-yet-marched
neighbour.  The system is solved as three nested fixed-point maps:

```
inner   : wavelength-coupled formal BVP at frozen (chi, eta)
          -> converged J on the ray   (FG_J_CMF_NONMONO_V1 + refresh)
middle  : statistical equilibrium at the lagged field
          -> populations, T, Ne        (SOLVEBA_V13 or JFNK)
outer   : great iteration
          -> opacities, source functions, new formal solve
```

The central stability requirement, from Knop et al. (2009, Paper II,
Sect. 3.2), is ordering:

> The populations must never take an update against a radiation
> field that is not the fixed point of the transfer map at the
> current opacities.

Every guard in `SOLVE_FOR_POPS` exists to enforce, approximate, or
recover from violations of that ordering.  The two failure modes it
prevents are mirror images: an SE update on a still-drifting field
re-poisons the field the sweep just converged, and a converged field
wasted on an SE step that explodes forces the whole outer map to
re-contract.

## 2. The Gauss-Seidel field solve

### 2.1 In-sweep structure (FG_J_CMF_NONMONO_V1)

Each great iteration sweeps all frequencies.  On each ray:

- Expanding cells (`sigma >= -1`): the wavelength discretization is
  purely upwind, so the march is live Gauss-Seidel -- quasi-analytic,
  one update per cell.
- Compressing cells (`sigma < -1`): the upwind-in-wavelength donor
  sits on a ray segment already marched this sweep; its intensity is
  read from the **lagged I_NEXT depth store** written by the previous
  great iteration.  The wavelength coupling is therefore a Jacobi
  term with a one-iteration lag.

`KN_RAY_ORDERED=T` (default) additionally orders each refresh pass so
rows are marched along the characteristic and frequency inside each
row is marched in the row's upwind direction; the only remaining lag
is the weak spatial corner (the Knop C-corner).  Under the legacy
alternating global passes, half of all strong links were stale every
pass and the residual pinned on marginal `|Q|/(chi+|Q|) -> 1` links.

`DO_NONMONO_FREQ_MATCH` remaps the frequency grid along the ray so
the upwind sample lands on the characteristic foot -- the
Ceulemans et al. (2025) fix for the spurious oscillatory eigenmodes
(`-1 < lambda < 0`) that survive any iteration ordering when the
per-cell Doppler drift exceeds the bin spacing.

Terminology note: "Gauss-Seidel" in this document means the Knop
in-sweep march ordering and the lagged-store fixed-point structure.
The Trujillo Bueno & Fabiani Bendicho (1995, ApJ 455, 646) GS/SOR
scheme is a different object -- it updates the source function
*inside* the formal-solution march itself, giving upper-triangular-
ALO convergence rates at local-operator cost (~2x over Jacobi, and
~sqrt(n) with optimal SOR).  Folding their population-level GS/SOR
into the march is a possible future acceleration of the inner map;
the refresh's adaptive omega already borrows the idea at the store
level.

### 2.2 The refresh (FG_KNOP_REFRESH_V1)

After each frequency sweep the lagged store is iterated toward its
fixed point on the sweep-captured opacities: alternating-direction
(or ordered) per-ray re-solves with adaptive omega damping, a
divergence abort (`1E3*MIN_RES`), a `1E40` retention scan, and a
stall exit at `KN_RAY_MAXPASS`.  A converged store means the *next*
sweep reads a self-consistent coupling field -- which is why the
maturity gate must test store drift AND sweep residual together,
not either alone.

### 2.3 In-sweep write guards

These run inside the sweep/refresh and feed the maturity metrics:

- **Capped writes**: non-finite or exploding intensities are dropped
  rather than committed to the store; the count is `GS_CAPPED` /
  `MOM_SW_CAPPED`.  A run in write-capping territory is producing a
  field the SE update must not trust.
- **Negative-J repairs** (`MOM_NEGJ`, `MOM_SW_NEGJ`): under
  `MOM_RF_KEEPNEG` repaired values never reach the stored field, so a
  bounded rate is a closure diagnostic, not corruption.
- **Quarantine** (`MOM_NQUAR`, `MOM_NPQUAR`): cells whose values
  repeatedly fail repair are excluded from the field; persistent
  quarantine is tracked separately.
- **Boundary-floor events** (`MOM_SW_BNDFL`): hard veto on maturity.
- **Significance floor**: metrics count only cells above `1e-4` of
  the column peak, so optically thin tail cells cannot dominate the
  convergence statistics.

### 2.4 The MOM path

`MOM_J_CMF_NONMONO_V1` solves a global tridiagonal system in depth
per frequency -- a moment-equation analogue, not the paper's
operator; there is no ordered march.  With `USE_LAM_ES=T` (the
paper-conform route, and the only one OBS supports) the MOM solver
does not run and its gate terms are vacuous.  When it does run, its
maturity is judged by sweep drift (fresh solve vs stored field),
which is what the populations actually experience, plus the refresh
delta when the refresh converged.

## 3. Field-maturity observables

All gates below consume these per-iteration quantities (see
`GET_FG_GS_MATURITY` / `GET_MOM_GS_MATURITY`):

| Metric | Meaning |
|---|---|
| `GS_MAXCHG` | max relative change of the lagged I_NEXT store on the last FG refresh |
| `GS_NBAD` | cells still above `KN_RAY_TOL` after refresh |
| `GS_NCONV` | cells whose stored value moved >50% (store churn) |
| `GS_CAPPED` | writes dropped as non-finite/exploding |
| `MOM_SW_MAXCHG` | max drift, fresh J solve vs stored J_NEXT |
| `MOM_SW_NBAD5/50` | cells drifting >5% / >50% |
| `MOM_SW_NHUGE` | cells drifting >1e6 (collapse detector) |
| `MOM_SW_NEGJ/NQUAR/BNDFL/CAPPED` | repairs, quarantine, boundary floor, dropped writes |
| `MOM_RF_CONV/DCHG/N5/N50` | refresh converged? retained fixed-point delta |

Thresholds (`solve_for_pops.f`):

- `GS_MATURE_TOL = 0.05` -- a sweep whose store moved <5% at the max
  cell is a maturity *candidate*.
- `GS_UNSAFE_TOL = 0.5` -- a max change >50% is runaway regardless
  of the bulk.
- Bulk rule: mature when fewer than ~0.1% of significant cells are
  still above tolerance or churning (a handful of marginally
  contractive cells must not pin the whole field immature).
- MOM two-tier gate: the **unsafe** tier (veto on everything -- >10%
  capped, >1% drifting >50%, >0.01% collapsed cells, one non-finite
  cell, >10% negJ/quar) vs the **marginal** tier (blocks release but
  eligible for `MOM_TRIAL`).

## 4. The freeze/release controller

`SOLVE_FOR_POPS` evaluates maturity at the top of each call and runs
the mode controller.  The decision chain:

```
GS_RAW_IMMATURE = capped>0  OR  FG-immature  OR  MOM-freeze-b
GS_FREEZE       = GS_RAW_IMMATURE  OR  (N_MATURE < GS_REL_STABLE)
```

- `N_MATURE` counts *consecutive* raw-mature calls; one churning
  sweep resets it.  `GS_REL_STABLE` (default 1, 2 in the hybrid
  config) is the stability hold: a single mature-looking sweep can
  be a coincidence of the churning field, and the premature release
  it enables re-poisons the field (the documented freeze/release
  limit cycle).
- During `GS_FREEZE`: `LAMBDA_ITERATION=T`, `FIXED_T=T`,
  `COMPUTE_BA=T`, `TEMP_CHAR='DIAG'` -- SOLVEBA still runs (so its
  demanded-correction diagnostics keep flowing) but takes a bounded
  Lambda update.
- The Lambda cap is churn-scaled: when the bulk churn fraction
  (`max(GS_NCONV,MOM_SW_NBAD50)/NTOT`) is below `GS_LAM_TRANSIENT`
  (2%) the field is near-mature and the correction is capped at
  `GS_LAM_CAP` (factor 3); in gross transient it keeps
  `MAX_LAM_COR` (factor 100) so cold start is unchanged.  A
  factor-100 committed step on a near-mature field re-poisons the
  store for ~10 iterations.
- `GS_HOLD_POPS` (default F) snapshots POPS before the forced-Lambda
  solve and restores after, so the field converges at *literally fixed*
  populations -- the strict Paper-II nested order. If allocation fails,
  the solver stops before updating populations. `GS_POPS_HELD` /
  `GS_ATTEMPT_APP` record the hold.
- `BA_LUM_GATE` (default F) demotes a scheduled BA step before solve
  when the luminosity ratio drifts, rolls back ordinary LAMBDA
  population corrections during the drift, blocks LAMBDA-to-BA exit
  and MOM_TRIAL, and caps near-mature release corrections at factor
  1.1. Transfer then advances at fixed populations until the ratio is
  stable. The current LAMBDA hold does not use BA_RES_ACC, which has no
  rollback for ordinary LAMBDA corrections.

### 4.1 Warmup and restart

- `N_LAM_RWARM` (default 4): the I_NEXT store is process-local and
  re-zeroed on every requeue, so the first N calls of every process
  are forced Lambda -- a BA step against a zeroed store is the
  classic requeue wreck.
- `GS_RWARM_SKIP`: a validated mature `FGINEXT_GS` checkpoint (store
  restored, last sweep <5%, no capped writes) skips the warmup --
  the field really is the fixed point already.
- `ITS_SINCE_START` (not the checkpointed `MAIN_COUNTER`) drives the
  warmup count, matching the process-local store lifetime.

### 4.2 Marginal release (`MOM_TRIAL`)

A field frozen only on the marginal tier (no unsafe events) can sit
in a stable asymptote the Lambda map cannot improve -- the FreeNe
run froze >150 iterations on ~0.12% drifting cells.  The trial
release fires only when *all* of these hold:

```
MOM_RAN & !MOM_UNSAFE_B & !MOM_BROAD & DO_BA_RES_ACC
  & !BA_SNAPSHOT_MISSING & N_FROZEN >= MOM_TRIAL_MIN
  & GS_CAPPED==0 & !GS_FG_IMMATURE & !LUM_BLOCKED
  & N_TRIAL_ACC < GS_TRIAL_LIM
```

The trial takes a real BA step (`T1=MAX_LIN_COR`) under the resacc
rollback net; an accepted trial refunds the freeze clock
(`N_FROZEN -> MOM_TRIAL_MIN`) since it demonstrated progress,
`GS_TRIAL_LIM` bounds total trials between raw-mature sweeps, and
`N_FROZEN` keeps counting so `GS_FREEZE_LIM`/`KNOP_STALL_MAX` still
kill a run that never matures.  **JFNK never fires on a trial** --
see section 6.

### 4.3 Watchdogs

- `GS_FREEZE_LIM` (default 0=off): hard stop after N consecutive
  frozen iterations -- fail loudly rather than burn requeues.
- `KNOP_STALL_MAX`: the job-wrapper stall bound.
- `GS_STATE` file: freeze state published each call for the wrapper.
- Plateau detector: 6-iteration ring of applied correction, second
  SE channel, and drift fraction; a <1.5x window flags
  `CONV_PLATEAU`.

## 5. The population-update side

When the field is live, the SE step itself is guarded:

- **Mode controller**: `LAMBDA_ITERATION` exits only after
  `LAM_N_STABLE` consecutive calls with applied correction
  `MAXCH_CTRL < LAM_EXIT`, and re-enters on one call at `LAM_ENTER`.
  `MAXCH_CTRL` is the *applied* (post-cap) correction, so a safely
  capped population decrease doesn't force unnecessary Lambda steps.
  `FRESH_BA_EXIT` rebuilds the BA on exit; `N_FULL_TREL` phases the
  temperature release (`T_HELD_PHASE`).
- **Resacc** (`DO_BA_RES_ACC=T`, mandatory on KNOP runs): snapshot
  `POPS_BA_SAVE`, take the BA step, re-evaluate the SE residual; if
  it grew, roll the step back and demote to Lambda.  `BA_TRUST_CUR`
  halves on rejection (`BA_NREJ_CONSEC`, `BA_NREJ_LIM`).  The
  rollback snapshot is part of the promise -- `BA_SNAPSHOT_MISSING`
  demotes any unchecked step (`ba_no_snapshot_demote`) and vetoes
  trials.
- **Saturated-demand demotion** (`BA_SAT_RAW`/`BA_SAT_LIM`): a full
  step whose raw demand stays at the sentinel is marching, not
  contracting -- demote and halve the trust radius.  Resacc can't
  anticipate the cliff (SE_RES2 stayed flat through the LAM_dbg
  blow-up while demand was pinned), so demand trend is the only
  pre-cliff signal.
- **Luminosity drift gate** (`BA_LUM_GATE`): when
  `|LUM_RATIO_OI - prev| > BA_LUM_DCHG * |ratio|`, demotes scheduled
  BA before solve, holds LAMBDA population updates, and blocks the
  LAMBDA->BA exit and trials. A *stability* test, not a value test --
  interacting snapshots need not equilibrate at 1. Transfer advances
  at fixed populations until the measured ratio settles.
- **FIX_BA**: after the first live BA build, later calls may recall
  the stored inverse (`BA_MATRIX_STORED` records a real write so a
  mid-call flag flip can't recall a matrix that never existed).

## 6. JFNK admission: the three phases

JFNK replaces the *middle* map's linear solve, not the iteration
structure.  Its residual has two modes with a ~100x cost asymmetry
(per-eval): the **frozen** residual uses stored radiative terms
(seconds per eval), the **coupled** residual re-runs the full Knop
formal solve (~250-400 s on m15) because `dJ/dn` is only available
through the sweep.

### 6.1 Admission chain (in order)

```
JFNK_THIS = USE_JFNK_SOLVER
        AND ITS_SINCE_START > JFNK_WARMUP        (process warmup)
        AND (.NOT. LAMBDA_ITERATION .OR. ANYITER) (live SE call)
        AND (.NOT. DO_NONMONO_KNOP .OR. .NOT. GS_FREEZE)
```

`GS_FREEZE` (not the raw `GS_FG_IMMATURE`) is the field test: it
subsumes the instantaneous flag, the `N_MATURE` stability hold, and
the trial veto.  Newton can only fire when the SE update is live and
the lagged field has been a demonstrated fixed point for
`GS_REL_STABLE` consecutive sweeps.  On a trial call the machinery
deliberately permits only a resacc-bounded BA probe -- an unbounded
Newton commit there would defeat the trial's purpose.

Then the demand tiers:

```
frozen tier : MAXCH < JFNK_FRZ_GATE  for JFNK_GATE_NREQ consecutive
              eligible calls  ->  JFNK_THIS, coupled=F
coupled tier: MAXCH < JFNK_COR_GATE  for NREQ consecutive eligible
              calls           ->  JFNK_COUPLED_THIS=T
```

`MAXCH` is the previous call's *demanded* correction -- the residual
scale of the SE problem.  The hit counters advance **only on
otherwise-eligible calls**: during warmup, Lambda, or freeze the
demand is measured against a field that is not the operator Newton
would see (and under freeze the deferred SE correction keeps
growing), so churn-era demand must not pre-arm the gate.

This is the literature-standard operator-split -> Newton handoff
made quantitative.  Hauschildt & Baron's production practice is to
iterate operator splitting until corrections are small and reserve
the full Newton linearization for the endgame; Arramy et al.
(2024, A&A 690, A12) diagnose their own Newton stalls the same way
and prescribe "a line search, a continuation method, and a few
operator-split iterations before the first Newton step".  Their
measured JFNK advantage over the RH92 operator-split reference is
only ~2x in formal-solution count -- which is precisely why this
design treats Newton as an endgame accelerator gated behind two
sustained criteria, not a wholesale replacement: the *field* must
be a converged fixed point (GS machinery) and the *iterate* must
be in the contraction basin (demand gate + hysteresis).  The
specific gate values (1e3 / 1e2 / NREQ=2) are empirical,
calibrated to this model family's demand trajectory; the two-map
structure itself is the part the literature prescribes.

### 6.2 Inside an activation

- Per-activation mode (`JFNK_COUPLED_LOC`): residual, JVP, rebuild,
  and diagnostics all read the activation flag; `INIT` resets the
  radiation readiness flags so frozen terms can't leak into a
  coupled activation.
- Frozen residual: `F(x) = STEQ_coll(n(x)) + STEQ_RAD + CMAT_RAD*g(x)`
  with `CMAT_RAD = BA_full - BA_coll` (radiative-only Jacobian, no
  double counting) and `g(x) = exp(x)-1` under `JFNK_LOGX` for *all*
  unknowns (pops, ED, T).  Exact at the baseline by construction.
- Coupled residual: each eval freezes I_NEXT, re-solves the field,
  restores all state; the eval residual is the Newton baseline.
- `JFNK_EPS` (1e-3): lower bound on the FD step -- **coupled mode
  only**.  The standard Knoll & Keyes choices (eps ~ b||x||/||v||
  with b=1e-6, or eps_mach^(1/r)*sqrt(1+||x||)/||v||) assume a
  deterministic F.  The sweep is iterative (data-dependent inner
  counts, repairs, quarantine), so F carries ~1e-4 relative noise
  and a machine-eps difference sees noise/h >> signal.  Pernice &
  Walker (1998, NITSOL) prescribe exactly this: scale the FD step
  to sqrt of the function noise level.  The frozen residual is
  deterministic arithmetic on stored arrays and keeps machine eps.
- `JFNK_MAX_EVALS` (40): per-activation eval budget; on exhaustion a
  sentinel marks the eval bad -> GMRES abort -> rejection path ->
  best-iterate/Lambda fallback.
- Newton internals (see [[JFNK_SOLVER]]): log-x unknowns,
  pseudo-transient continuation (`dtau` adaptive), Eisenstat-Walker
  eta floored at 3e-4, per-depth equilibrated-LU preconditioner
  (`BA_full + ridge`, null-row -> identity, Jacobi fallback),
  `MAX_JFNK` commit clip in ratio space.
- Fallback chain: any failure (bad F(0), GMRES abort, trust-region
  rejection, non-finite step) -> restore best iterate or LAMBDA-mode
  SOLVEBA.  `SUCCESS=F` never commits a worse state.

## 7. Per-iteration decision tree

```
top of SOLVE_FOR_POPS:
  read FG/MOM maturity from last sweep
  if requeue-warmup window and no mature checkpoint -> forced LAMBDA
  if scheduled BA step but field unsafe          -> demote to LAMBDA
  if resacc snapshot missing                     -> demote to LAMBDA
  if BA_LUM_GATE and ratio drifting              -> demote scheduled BA
  compute GS_RAW_IMMATURE, N_MATURE, GS_FREEZE
  if frozen:
      if MOM_TRIAL conditions  -> bounded resacc-gated BA probe
      else                     -> forced LAMBDA (churn-scaled cap)
  JFNK gate:
      eligible? (warmup done, live call, not frozen)
          -> accumulate frz/cpl hits on sub-gate demand
      admit  -> JFNK_SOLVE_CORR_V1 (frozen or coupled per tier)
      else   -> SOLVEBA_V13 (Lambda or full BA per mode controller)
  if GS_POPS_HELD (GS freeze or lum hold) -> restore pre-solve pops
  mode controller for NEXT call (exit/enter LAMBDA, lum gate,
      saturation demotion, trust radius, plateau detector)
  ITER_SUM / ITER_HISTORY / GS_STATE reports
```

For strict GS ordering, configure `DO_LAM_IT=T`. The stock-cycle
`DO_LAM_IT=F` branch retains scheduled BA steps during `GS_FREEZE`.

## 8. Failure modes and the guard that catches each

| Observed failure | Guard | Evidence |
|---|---|---|
| SE update on a churning field blows up the next solve | `GS_FREEZE` ordering | KNOP_MOM_TFIT it16-19 |
| Release on a lucky marginal sweep re-poisons field | `N_MATURE`/`GS_REL_STABLE` hold | LAM_dbg it7, twice |
| Deferred correction released at cold-start cap re-poisons | `GS_LAM_CAP` on release + `GS_REL_N` | same event |
| Factor-1e5 committed step -> luminosity 2e10 -> 2e16 | Lambda cap while immature | m15 JFNK seq it2 |
| BA step against zeroed I_NEXT after requeue | `N_LAM_RWARM` | requeue wrecks |
| Unchecked BA step with failed snapshot arm | `ba_no_snapshot_demote` | LAM_dbg it39-48 |
| BA step while inner flux still swinging | pre-solve demotion via `BA_LUM_GATE` | LAM_alt it13-35 |
| LAMBDA population step during flux drift | `BA_LUM_GATE` population hold + factor-1.1 release cap | JFNK_cpl it9-14 |
| Demand-pinned marching past the basin cliff | `BA_SAT_*` demotion | LAM_dbg it36-47 |
| Rollback overwriting a resacc-rejected step | `BA_STEP_REJECTED` check | m15 it64 wreck |
| JFNK at demand ~1e10 -> GMRES stall -> MAX_LAM wreck | `JFNK_*_GATE` + `GS_FREEZE` veto | m15 opt it11 |
| ~12 wasted formal sweeps on diagnostics/F(0) dup | `JFNK_DEEPDIAG`, `jfnk_fx_fresh` | audit |
| FD JVP drowned in iterative-solver noise | `JFNK_EPS` floor (coupled only) | JVP-vs-BA probe |
| Field marginal forever -> run never converges | `MOM_TRIAL` + stall watchdogs | FreeNe 150-iter freeze |

## 9. OUTGEN observability

Per iteration, `ITER_SUM` prints a tagged block:

```
ITER_SUM: it=  mode= -> ba= frozen= held=
    pops           demand / applied correction
    se residual    all cells / excl. worst cols
    ba             trust / consec rejects
    lum ratio      inner/outer luminosity
    fg             bad= churn= cap=  lagged field drift
    mom            n5= n50=/ntot neg= quar= pquar= cap=
    jfnk           demand vs frz gate, hits frz=N/M   [OFF|DEFER|ARMED]
    jfnk_cpl       demand vs cpl gate, hits cpl=N/M, coupled=T/F
```

Gate messages, grep-able:

- `JFNK gate: deferred -- warmup, call N of M`
- `JFNK gate: deferred -- GS field not ready, nmature= FG maxchg= MOM n50=`
- `JFNK gate: deferred -- LAMBDA-designated iteration`
- `JFNK gate: last demanded MAXCH= vs gate= hits= -- deferred`
- `JFNK admission: MAXCH= hits frz= cpl= coupled=`
- `iter_reason` tags: `gs_warmup`, `gs_immature_demote`,
  `gs_freeze_forced`, `mom_trial_release`, `jfnk_frz`, `jfnk_cpl`,
  `jfnk_defer`, `hold_pops`, `lam_stable_exit`, `lum_drift_hold`.

The intended run signature is: `gs_warmup` -> `gs_freeze_forced` /
Lambda steps -> live BA steps -> `JFNK gate` defers with growing
hits -> `jfnk_frz` admissions -> `jfnk_cpl` near the root.  A run
that shows `jfnk_cpl` while `frozen=T` or `nmature < GS_REL_STABLE`
in the same call is a bug.

`ITER_HISTORY` appends one row per iteration (mode, frozen/held,
raw/applied correction, se/se2, trust, lum, drift fractions) for
post-run plotting; `GS_STATE` is the wrapper-facing freeze state.

## 10. VADAT reference

`DO_NONMONO_KNOP=T` hard-requires (reader STOPs): `GLOBAL_LINE=BLANK`
(or `NONE` with all present `TRANS_*=BLANK`), `LAM_SET=F`,
`INC_GRID=F`, and no `USE_FIXED_J` / relativistic / plane-parallel
modes; `DO_BA_RES_ACC=F` warns. The OBS reader applies the same
line-treatment gates to CMF_FLUX_PARAM. Added 03-Oct-2026.

| Key | Default | Role |
|---|---|---|
| `DO_NONMONO_KNOP` | F | master flag for this path |
| `KN_RAY_ORDERED` | T | ordered-march refresh |
| `DO_NONMONO_FREQ_MATCH` | F | Ceulemans frequency matching |
| `USE_LAM_ES` | -- | FG/ray path; F enables MOM side |
| `N_LAM_RWARM` | 4 | forced-Lambda warmup per process |
| `GS_REL_STABLE` | 1 | consecutive mature sweeps to release |
| `GS_REL_N` | 1 | post-release damped-Lambda window |
| `GS_HOLD_POPS` | F | strict fixed-pops freeze |
| `GS_FREEZE_LIM` | 0 | freeze watchdog (0=off) |
| `MOM_TRIAL_MIN` | 30 | first frozen iter eligible for trial |
| `GS_TRIAL_LIM` | 15 | max accepted trials between mature sweeps |
| `DO_BA_RES_ACC` | -- | residual-acceptance rollback (mandatory here) |
| `BA_LUM_GATE`,`BA_LUM_DCHG` | F | demote BA and hold LAMBDA populations on ratio drift |
| `BA_SAT_LIM`,`BA_SAT_RAW` | -- | saturated-demand demotion |
| `LAM_ENTER/EXIT/N_STABLE` | LAM_VAL/.../1 | mode-controller hysteresis |
| `USE_JFNK` | F | arm JFNK |
| `JFNK_WARMUP` | -- | calls before JFNK may fire |
| `JFNK_COUPLED` | F | permit the coupled tier |
| `JFNK_FRZ_GATE` | 0 | frozen-tier demand gate (%; 0=use COR) |
| `JFNK_COR_GATE` | 0 | coupled-tier demand gate (%; 0=off) |
| `JFNK_GATE_NREQ` | 1 | consecutive sub-gate calls to admit |
| `JFNK_MAX_EVALS` | 0 | eval budget per activation (0=off) |
| `JFNK_EPS` | 0 | JVP FD-step floor (coupled mode) |
| `JFNK_DEEPDIAG` | F | ~12-eval diagnostic suite |
| `JFNK_M` | 10 | GMRES restart dim |
| `JFNK_NITS` | 10 | Newton iterations per activation |
| `JFNK_MAX_RESTART` | -- | activation restarts per call |
| `CONV_DIAG_LVL` | 2 | diagnostic verbosity |

Hybrid model (`configs/dessart_m15_2p82d_JFNK_hyb.yaml`):
`FRZ_GATE=1e3`, `COR_GATE=1e2`, `NREQ=2`, `MAX_EVALS=40`,
`EPS=1e-3`, `M=30`, `WARMUP=6`, `GS_REL_STABLE=2`.

## 11. Literature alignment and honest gaps

What the literature verifies, and where this implementation still
stands on engineering judgement rather than published results:

- **Verified structure.** The freeze ordering (Knop Paper II),
  the operator-split -> Newton staging (Hauschildt & Baron;
  Arramy et al. 2024), PTC globalization (Kelley & Keyes 1998 --
  succeeds where line search/trust region stagnate on
  PDE-structured problems), EW forcing with a floor inside
  Arramy's useful 3e-4..3e-2 band, and the noise-scaled FD step
  (Pernice & Walker 1998) are all straight from the literature.
- **Ordered-march convergence rate is theory-derived.** The claim
  that `KN_RAY_ORDERED` leaves only the weak C-corner lagged comes
  from the stencil derivation in
  [[KNOP_CONVERGENCE_STRATEGY]] (marked UNVERIFIED there).
  Runtime confirmation on a compressing model is still pending.
- **No line search on the committed Newton step.** Arramy names a
  line search as a missing ingredient for out-of-basin starts.
  This implementation substitutes the admission gates themselves
  (Newton is never allowed out-of-basin), plus best-iterate
  fallback, commit checks, PTC, and the resacc net on BA steps.
  A backtracking line search on the JFNK commit is a candidate
  upgrade if activations prove to fail rather than defer.
- **Preconditioner updates are deliberately rare.** Arramy's
  Table 2 shows a PC rebuild roughly doubling per-call cost; the
  design rebuilds the BA preconditioner only on baseline changes
  and falls back to Jacobi on factor failure -- consistent with
  keeping PC refreshes expensive and infrequent.
- **Block-diagonal ALO only.** Olson & Kunasz (1987) measure 2-7x
  fewer iterations from tridiagonal/pentadiagonal approximate
  operators.  `JFNK_USE_BANDED` + `NUM_BNDS>=3` exist in the
  preconditioner path but are untested on this model; the
  cross-depth radiative coupling is currently dropped from both
  the frozen matvec and the PC (self-consistently).
- **Outer-map acceleration is complementary.** CMFGEN's NG
  acceleration (`do_ng_accel_band`) and Anderson acceleration act
  on the outer fixed-point map; they accelerate a contraction that
  already exists and are not a substitute for the inner gates.
  The plateau detector addresses the false-convergence failure
  stationary iterations are known to exhibit (Lambert et al. 2015
  style stalls).
- **Thresholds are empirical.** 5% / 50% / ~0.1% bulk, the churn
  fraction, the gate values, and the stability counts were
  calibrated against the documented blow-ups in
  [[convergence_issues]], not derived.  They are the
  implementation's expression of "sustained contraction", which is
  the literature concept; the numbers may need re-calibration on a
  different model family.

## 12. Design invariants for future edits

1. **Newton only on a live SE update.** `GS_FREEZE` must remain in
   the JFNK admission test; weakening it to `GS_FG_IMMATURE` alone
   re-opens the trial-iteration hole.
2. **Ordering before accuracy.** A converged field spent on a safe
   step beats an accurate Jacobian spent on a transient one.  Do not
   add "smarter" SE updates that bypass the maturity gates.
3. **Hysteresis is load-bearing.** Every counter (`N_MATURE`,
   `CNT_STABLE_LAM`, `JFNK_*_HITS`) exists because a single lucky
   call is not evidence.  New gates need consecutive-call semantics
   and reset-on-ineligible.
4. **Every unchecked step needs a net.** Full steps go through
   resacc; JFNK commits go through best-iterate + commit checks;
   Lambda steps are capped.  If a new path can commit populations,
   it needs its own rollback.
5. **The two JFNK residuals must agree at the baseline.** Frozen
   terms are exact at `x=0`; any change to the stored-term split
   (`STEQ_RAD = STEQ_full - STEQ_coll`, `CMAT_RAD = BA_full -
   BA_coll`) must preserve `F(0) = outer STEQ` or the F(0) identity
   test will (correctly) reject the activation.
6. **Log-space consistency.** Under `JFNK_LOGX` every consumer of x
   (residual, matvec, JVP, commit clip, rebuild propagation) uses
   the same `exp(x)-1` fractional map with the same +/-40 clip.
7. **Coupled evals are the budget.** Anything that spends formal
   sweeps (diagnostics, duplicate F(0), restarts) belongs behind a
   flag or a counter, not in the default path.
