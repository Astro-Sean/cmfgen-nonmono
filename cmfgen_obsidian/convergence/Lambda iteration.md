# Lambda iteration

The bounded fixed-point update: the safe default mode of the
population solver, the state every controller falls back to, and the
only mode allowed while the radiation field is immature.

## What a Lambda step is

With `LAMBDA_ITERATION=T`, `SOLVEBA_V13` does not solve the linearized
SE system; it applies the SE residual as a bounded fixed-point
correction (`POPS *= (1 - STEQ*SCALE)` with `TEMP_CHAR='DIAG'` under
the freeze controller). No Jacobian is trusted - the update is the
Lambda-iteration map: hold the radiation field, push populations
toward SE, re-solve the field.

#status/verified properties that make it the safe mode:

- Contractive in practice on this pipeline: the pure-Lambda F baseline
  (job 11965659) contracted corrections ~10%/iter while gate-locked.
- Bounded by construction: `MAX_LAM` (`MAX_LAM_COR`) caps the
  committed factor, and the churn-scaled cap `GS_LAM_CAP` tightens it
  to factor ~3 on a near-mature field ([[Step limiting and rollback]]).
- T is held (`FIXED_T=T`), so the step cannot couple population drift
  into the temperature equation.

## The mode controller

`SOLVE_FOR_POPS` decides per call whether the next call runs Lambda or
full (`solve_for_pops.f`, the `MODCTL` block). Hysteresis, in force
since 30-Sep-2026:

- `LAM_EXIT` / `LAM_ENTER` / `LAM_N_STABLE`: exit Lambda only after
  `LAM_N_STABLE` consecutive calls with applied correction
  `MAXCH_CTRL < LAM_EXIT`; one call at `LAM_ENTER` re-enters. Defaults
  (`LAM_VAL`, 1 stable call) reproduce the legacy single-threshold
  switch; separated thresholds kill the one-lucky-iteration release.
- On the KNOP path the test uses `MAXCH_CTRL` - the **applied**
  (post-cap) correction - so a raw saturated demand that was safely
  capped does not force a pointless Lambda step.
- `NUM_LAM` sets the alternation cadence (N Lambda calls per full
  step). On the optimized SOB profile `NUM_LAM=BA_PROBE_LAM=20`, so a
  held Lambda window re-measures the BA direction every 20 calls while
  `BA_ENTER_RAW` remains armed.
- Overrides that force Lambda regardless of the controller: `RD_LAMBDA`
  (IN_ITS `DO_LAM_IT`), `FORCE_LAMBDA` (warmup/freeze), R-grid
  revision, hydro revision, clump revision, NG kick, and the
  `IT_SPECIFIER` `LAM` directive ([[Warmup and staged release]]).

## Auto-transitions out of Lambda

- `DO_LAM_AUTO` (IN_ITS, default T): the stock path flips
  `RD_LAMBDA` to F at `MAXCH_CTRL < 50`. With `DO_NONMONO_SOB` or
  `DO_NONMONO_KNOP`, the initial flip also requires `LAM_N_STABLE`
  consecutive calls with applied `MAXCH_CTRL < LAM_EXIT`; forced
  LAMBDA calls do not count. Luminosity ratio and drift are diagnostic
  only because shock deposition invalidates boundary agreement. The
  stable-wait and handoff are logged in OUTGEN. The stock behavior is unchanged when
  both non-monotonic flags are off.
- `BA_PROBE_LAM` rd_lam probe (11-Oct-2026): the stable-exit gate has
  no escape when the LAMBDA attractor floors above `LAM_EXIT` -- the
  m15 SOB run parked at ~45% demand with `stable=0/3` because every
  probe/escape path sits behind `.NOT. RD_LAMBDA`. While
  `DO_BA_RES_ACC`, `BA_PROBE_LAM > 0`, `BA_DAMP_CAP > 0`, and the
  demand sits in `LAM_EXIT <= MAXCH_CTRL < LAM_ENTER`, every
  `BA_PROBE_LAM` Lambda calls arm `BA_RDLAM_PROBE`: `cmfgen_sub`
  leaves the next call full (with the `BA_DAMP_CAP` factor bound),
  resacc adjudicates it the call after. A clean verdict arms
  `BA_RDLAM_ACCEPT`, which pins `CNT_STABLE_LAM=N_STABLE_LAM` so the
  ordinary `DO_LAM_AUTO` switch performs the handoff. A degraded probe
  is fully rolled back on its Lambda follow-up, including the Knop
  lagged field, and that call skips the stale Lambda population solve;
  the next probe is `BA_PROBE_LAM` calls away. A clean probe remains
  committed. Probe calls cannot run JFNK and a grid revision the same
  call still forces Lambda. The shared residual-acceptance guard covers
  both SOB and Knop/JFNK; inert unless the probe controls are armed.
- `USE_FIXED_J` models switch to variable J at the same `MAXCH < 50`
  threshold, optionally preceded by a grey-T iterate (`DO_GT_AUTO`).
- `DO_T_AUTO` releases `FIX_T` at `MAXCH < 50` once a non-Lambda step
  has run.

## Observed failure modes

- **Lambda<->BA alternation** ([[convergence_issues]] Sect. 1.1,
  [[Failure modes and guards]]): a BA step on a still-moving field
  demands 1e4-1e8% corrections and is forced back to Lambda; the cycle
  repeats. The maturity gate plus churn-scaled cap is the answer; on
  the Type IIn model the kicks are absorbed but not yet contracting.
- **Forced-Lambda plateau** (Sect. 10.1): on m15 the residual negJ
  population kept the MOM gate frozen ~200 iterations, and pure Lambda
  oscillated 4->27% at depth 62-63 without a trend - the protective
  mode is not guaranteed to contract, it is guaranteed to be bounded.
- **LAMBDA-attractor misreading** (04-Oct-2026): in the alternating
  regime a Lambda call's demand reflects the Lambda fixed point, not
  the BA operator - gating on `MAXCH` alone admitted JFNK at it73 on a
  95% demand while the last full demand was 4.8e4%. Gates now read
  `BA_LAST_RAW`.
- **rd_lam plateau park** (11-Oct-2026): on the m15 SOB run the
  Lambda demand contracted 166->44% by it33, then rebounded to ~49%
  and held -- the Lambda attractor floored above `LAM_EXIT=20`, so
  `CNT_STABLE_LAM` stayed 0/3 forever and the run would have burned
  all 500 iterations without one full call. A stable threshold is not
  reachable on a trajectory whose fixed point sits above it; the
  `BA_PROBE_LAM` rd_lam probe is the bounded escape.

## Related

- [[Convergence strategy]] - where Lambda sits in the staged plan
- [[BA matrix and SOLVEBA]] - the full-mode alternative and resacc
- [[Field maturity and GS freeze]] - the freeze controller that pins
  `LAMBDA_ITERATION`
- [[Acceleration methods]] - Ng/AV acceleration applied on top of
  Lambda iterates
- [[Convergence criteria]] - `LAST_LAMBDA` in the exit test
