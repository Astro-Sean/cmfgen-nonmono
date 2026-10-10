# Convergence Issues - dev7 non-monotonic CMFGEN

Status as of 2026-10. Scope: `cmfgen_nonmono_dev_V7` on the Type IIn
10 d model (`TypeIIn_10d_DEV7_KNOP_TRIDIAG*` family). Claims are tagged
VERIFIED (measured), SUPPORTED (consistent evidence), UNVERIFIED
(hypothesis), or CONTRADICTED (falsified by a run).

## 1. Open issues

### 1.1 BA kick magnitudes are not contracting (F path)

VERIFIED - job 11969828 (`gatebulk`, `USE_LAM_ES=F`).

After the sweep-drift gate opened at it25, BA steps perturb the field
hard enough to re-trip the unsafe bound every cycle:

| iter | mode | max % increase | max % decrease | recovery |
|---|---|---|---|---|
| 28 | BA computed | 3.4E+03 | 1.3E+17 | 1 frozen Lambda it |
| 31 | BA computed | 1.5E+04 | 1.5E+04 | 2 frozen Lambda its |
| 35 | BA computed | 2.7E+04 | 1.4E+04 | 1 frozen Lambda it |

The kicks are bounded and fully absorbed each cycle (the gate catches
the perturbation within one iteration and `DO_LAM_AUTO` recovers), but
the sequence 3.4E3 -> 1.5E4 -> 2.7E4 is increasing, not decreasing.
Whether the Lambda<->BA cycling nets faster convergence than the
pure-Lambda baseline (job 11965659, corrections contracting ~10%/iter
while gate-locked) is not yet established.

Likely contributor (UNVERIFIED): the BA Jacobian uses the clamped
sigma (`>= -0.999`), not `SIGMA_TRUE`, so the linearization is
approximate exactly where the field is hardest to update
(compression zones, depth ~150 and the inner boundary dominate the
max corrections).

### 1.2 MOM refresh never converges internally (F path)

VERIFIED - all F variants (bndflr, gsdesc, keepneg, gatebulk).

The `MOM_KNOP_REFRESH` pass residual plateaus at ~1.3-2.4 and the
per-pass `nbad5` count stays in the millions because marginally
contractive cells keep moving >5% per pass indefinitely. The refresh
now exits on the stall test (two consecutive same-direction residual
ratios >0.90 after pass 8) and restores the best pre-refresh field,
so the waste is bounded (~70-100 s/iter vs ~420 s before), but the
underlying map is not self-consistent on this model.

Consequence: the refresh residual cannot be used as a maturity metric;
BA gating had to move to sweep-level drift (`MOM_SW_*`).

### 1.3 Bulk field drift remains high (F path)

VERIFIED. At gate-open (it24-25) the sweep `nbad5` count was just
under the 10% threshold (3.1M of 31.4M cells). Roughly one cell in
ten still drifts >5% between sweeps even in the "mature" state. The
gate is lenient by construction (10% plus `maxchg<=10`,
`nbad50<=0.1%`, negJ/quar backstops); whether this threshold is tight enough for BA
to make net progress is part of issue 1.1.

### 1.4 Single-cell maxchg flicker (F path)

VERIFIED. The worst sweep cell oscillates transiently
(`sweep maxchg` measured 2.1 -> 46.8 -> 7.4 -> 30.2 -> 2.8 across
it7-13) while the bulk stays calm (`nbad50` ~300-3000). The gate can
therefore flicker open/closed on one noisy cell rather than the bulk
state. Currently harmless - a spurious `maxchg>10` only costs one
extra Lambda iteration - but it makes gate-open timing stochastic.

### 1.5 Negative-J production is suppressed, not eliminated (F path)

VERIFIED. `MOM_RF_KEEPNEG=T` collapsed the significant negative-J
count from ~1E6 to ~2E3 per sweep, and the quarantined fraction is
~0.006% of cells. The residual ~2E3/iter population persists at
quasi-steady state - these cells are repaired every sweep and keep
their prior stored value. Their location is diffuse (argmax drifts
over (d,f) space), not one pathological site, so there is no single
grid fix.

### 1.6 Warm restart does not preserve field maturity (F path)

VERIFIED - resubmit of 11967617 -> 11969828. `FGINEXT_GS` and
`MOMJNXT_GS` survive requeue/resubmit, and populations continue from
checkpoint (identical corrections), but the checkpoint carries no
maturity flag: the run restarted warm-up (`LTE_EST=T` path) and the
field had to re-earn maturity over ~25 iterations. On the 24 h
partition this is tolerable but it wastes several hours of wall time
per resubmit and makes BA-mode restart impossible.

### 1.7 OBS/cmf_flux cannot use the MOM path

VERIFIED. The OBS stage still requires `USE_LAM_ES=T`; the moment
path has not been reconciled into `cmf_flux`. Any spectrum production
from an F-path model therefore needs the formal solution at the OBS
stage, or the OBS work must be ported. See `audits/
USE_LAM_ES_KNOP_AUDIT_20260917.md` addendum.

### 1.8 Direction-aware live reads did not help (F path)

CONTRADICTED as a fix - `MOM_RF_GS_DESC=T` (job 11966433) let
descending refresh passes read live `J_NEXT_FRQ(:,L+1)` along the
red-to-blue compressing characteristic. Measured it1 best residual
40.8 vs 13.9 frozen-read baseline, tracking the baseline (~neutral)
thereafter. The live descending read does help its own passes
systematically (~45-75 vs ~90-220 ascending) but partially
re-cascades the field. Option retained in VADAT, default F.

### 1.9 T path endgame also not demonstrated converged

SUPPORTED. `USE_LAM_ES=T` reaches BA mode earlier and more cheaply
per iteration (T-full entered BA at it21), but it shows the same
Lambda<->BA alternation with kicks of order 1E2-2E3% at it41-45 on
the debug run. Neither path has demonstrated a monotonically
contracting endgame on this model; the T path is "further along" per
unit wall-clock, not converged.

### 1.10 Metric/init quirks to be aware of

- `nbad5 > ntot` can appear on iteration 1 (accumulation before the
  first reset) - harmless, do not treat as a defect.
- Quarantined cells are excluded from `nbad5/nbad50/maxchg` (they are
  reported separately as `quar`); without this exclusion the residual
  is pinned by cells that can never update (observed on keepneg it1:
  residual 3516, a metric artifact).
- Refresh-residual numbers and `MOM sweep:` numbers answer different
  questions: map self-consistency vs stored-field drift. Do not mix
  them in gate logic or reporting.

## 2. Performance issues (not convergence-blocking)

- MOM refresh adds ~70-100 s per great iteration even with the stall
  exit; ~420 s/iter before the exit was added. (VERIFIED)
- Candidate work listed in [[future_fixes]]: `cpus_per_task` 4->16
  scaling test, MOM coefficient caching, strided `J_NEXT_FRQ` store
  access, coefficient precision experiments, frequency-sweep OpenMP
  tail profiling. None validated yet. (UNVERIFIED gains)
- Infrastructure: T-path debug job 11966624 died on a bus error
  during graceful-stop writeout (checkpoint write), not a physics
  failure; warm restart recovered. Watch for recurrence.

## 3. Fixed issues (for reference)

- Frozen-snapshot refresh was effectively Jacobi-like; direction-aware
  reads now available but neutral (see 1.8).
- Repaired negative-J cells reseeding downstream frequencies - fixed
  by `MOM_RF_KEEPNEG=T` quarantine.
- Refresh wasting many passes while stalled - fixed by the
  same-direction stall exit.
- BA gate pinned by refresh residual / `MOM_NEGJ>0` forever - fixed by
  the `MOM_SW_*` sweep-drift gate (`MOM_GATE_BULK=T`).
- Gate permanently closed on Type IIn 10 d - fixed; opened it25,
  first BA it28.

## 4. Key runtime evidence

- `USE_LAM_ES=F` decisive run: job 11969828,
  `TypeIIn_KNOP_TRIDIAG_gatebulk`, config
  `configs/TypeIIn_10d_DEV7_KNOP_TRIDIAG_gatebulk.yaml`
  (`MOM_RF_KEEPNEG=T`, `MOM_GATE_BULK=T`, `MOM_RF_GS_DESC=F`).
- `USE_LAM_ES=T` reference: jobs 11966624 (debug resubmit) and
  11965755 (`KNOP_LN` full, BA from it21).
- Pure-Lambda F baseline: job 11965659.
- Source commits: `28b0cdd`, `7b1e570`, `069ed05`, `c327942`,
  `70ad426`, `5ed2e98`.

## 5. What "converged" still requires

1. BA kick magnitudes showing a decreasing trend across cycles
   (currently increasing 3.4E3 -> 2.7E4).
2. Sustained reduction of `nbad5` sweep drift below the 10% gate,
   i.e. the field staying mature without re-freezing after BA.
3. A bounded or declining `sweep maxchg` (currently flickers up to
   ~50 on single cells).
4. Population corrections falling under the model's convergence
   criterion in consecutive iterations, not just recovering after
   each kick.
5. An OBS-stage run producing the spectrum (requires `USE_LAM_ES=T`
   for cmf_flux regardless of the model-run path).

## 6. Audit 27-Sep-2026 (m15_2p82d validation, jobs 11997835 LAM / 11997838 MOM)

Fixed (VERIFIED, committed):

- `fc12e0f` Makefile: `new_main/subs` compiled before the current
  `CONTROL_VARIABLE_MOD` existed; clean builds only worked with stale
  `.mod` files.  The `cmfgen` target also linked the LTE tool, whose
  `WR_ABUND_INFO_V2` has no implementation in this tree.
- `e2b9bd2` FG_J_CMF_NONMONO_V1: `DLN_C` was shared across the OpenMP
  ray loop.  chi_hat and Q could be built with different frequency
  spacings, so |Q| could exceed 1 on a compressing point.  Runtime
  effect on m15_2p82d was nil (iterations 1-5 reproduce the racy
  build to all printed digits), but the defect was real.

MOM path (`USE_LAM_ES=F`) regression (VERIFIED):

- Job 11979734 (build of 25-Sep) gave L(d=1)=2.6E10 Lsun at it1.
  Jobs 11987779 and 11997838 (with the uncommitted edits committed as
  `42dc0e9`) give 2.7E20 at it1 with the same LTE populations.
- Mechanism: `42dc0e9` replaces a significant negative R^2 J by the
  previous valid field (previous frequency / lagged J_NEXT) instead
  of |XM|/10.  With ~17 negative rows per frequency (992101 per
  sweep) the substitution pins J to its bluer value across the
  sweep, so J can only ratchet upward along frequency; depth 1 ends
  at |R^2 J| ~ 1E19 in the optical (MOM_J_NONMONO_ERRORS).  Neither
  repair is physical.
- Root cause (SUPPORTED): the negative rows have an M-matrix
  (TB=0.59, TA=-0.08, TC=-0.10 at depth 48, freq 94.8) and a negative
  RHS.  The RHS sign comes from the lagged H and EPS couplings, which
  always read the bluer frequency (`RSQHNU_PREV_NM`, `RSQJNU_PREV_NM`)
  even where the flow compresses.  With |sigma_true| up to 120 the H
  coupling weight `HS=WPREV/(1+W)` approaches 1, so the wrong-direction
  explicit term dominates the row.  Fixing this needs a lagged H_next
  (and J_next for EPS) store and a refresh that converges; earlier
  EPS-redirect attempts (Sect. 2 of [[outstanding_issues]]) made the
  map expand.  Classified FAIL; do not use for production.

LAM path (`USE_LAM_ES=T`) (VERIFIED evidence, SUPPORTED cause):

- it1: L=1.98E10, max change 5.87E16 % at depth 65, populations
  driven negative at depth 41 (155 % decrease), pure LAMBDA, T fixed.
- it2: L=3.6E17, then L decays ~8 %/iteration and jumps at each BA
  step (job 11979727: 2.5E13 at it77 -> 8.5E17 at it78).  From it~170
  the I_NEXT store is capped at 1E20 on 1.5-2.6E8 cells per sweep and
  the refresh reports the capped field as converged.
- The stock-solver control (`m15_2p82d_knop_SOB`, FG_J_CMF_V13 with
  clamped sigma, DO_NONMONO_KNOP=F) shows 9.5E31 % at it1 on the same
  hydro.  The it1 explosion is therefore not caused by the Knop
  solver: the LTE start at the hydro T_gas with unlimited Lambda
  corrections (`LAM_SCALE_OPT=NONE`, `MAX_LAM=1E10`) is inconsistent
  with the radiation field of this shocked, compressing structure.
- Recommended setup changes (not applied): limit Lambda corrections
  (`LAM_SCALE_OPT=LIMIT`, `MAX_LAM` ~1E2), start from dilute-LTE or
  proxy populations rather than LTE at T_gas, and allow T to relax
  once the radiation field is mature.  Re-test the KNOP LAM path only
  after the control run converges its first iterations.

## 7. MOM upwind-consistency fix (28-Sep-2026, this build)

The Sect. 6 root cause is now corrected in source: the lagged H and
EPS couplings of the H equation always read the bluer frequency even
at compressing midpoints, so the lagged map's fixed point satisfied
the wrong equations (wrong-direction explicit terms dominate rows
whose upwind weight approaches 1; the refresh could then never
converge).

Changes in `plane/mom_j_cmf_nonmono_v1.f`:

- New lagged store `H_NEXT_FRQ(ND_NM,NCF)`, partner of `J_NEXT_FRQ`:
  R^2*H written after each frequency solve with the same
  JNEXT_RELAX-damped update.  Compressing midpoints (GAMH_SIGNED<0)
  read slot k+1 as their upwind field; expanding midpoints keep the
  bluer RSQHNU_PREV_NM chain.
- `HS_NM`/`EPS_PREV_NM` coefficient fixup: at compressing midpoints
  the explicit coefficients now use the current-frequency weights
  (HS=W/(1+W), EPS_PREV=EPS), matching the PSI_NM convention for the
  J coupling.
- `HUP_NM`/`JUP_NM` upwind fields used consistently in the interior
  RHS, the OUT_BC_TYPE>1 boundary row, the negative-J rebuild, and
  the H reconstruction -- all four sites now agree.
- EPS midpoint terms rewritten as a per-midpoint scatter (row m gets
  +2E*J_up, row m+1 gets -2E*J_up) -- algebraically identical to the
  old per-row form, now with the direction-aware field.
- H drift folded into MOM_GS_MAXCHG only (not NBAD5/50): H crosses
  zero legitimately, so fraction counts would pin the gate on
  sign-change noise; the max-norm still catches a runaway H field.
- MOMJNXT_GS checkpoint VER 4 writes H_NEXT_FRQ after J_NEXT_FRQ;
  restore accepts VER 3 (J only, H cold) or 4.  Rebuild remap
  (HNEXT_OLD) and refresh (HNEXT_RF_READ snapshot, HNEXT_SAVE
  restore) updated alongside.
- MOM_RF_KEEPNEG negative-J repair retained as a backstop: with the
  upwind-consistent RHS it should fire rarely; watch the
  significant-negJ counter in the test run.

Open question for the test run: the |EPS|/HU~17 coupling may still
converge slowly inside the refresh even with consistent fields.
Residual negative RHS contributions at sigma_true<0 are physical
(the sigma*n*J frequency flux is signed); whether they produce
significant-negJ cells is now a fixed-point property of the correct
equations, not a discretization error.  Verify on the interactive
submission: negJ count should collapse vs jobs 11997838, refresh
residual should contract, MOM_GS_MAXCHG should decline.

## 8. negJ-repair ratchet fix + interactive validation (27-Sep-2026)

The Sect. 7 upwind-consistency build was submitted to the
interactive partition (job 11998467).  Result: it1 L=2.88E20 with
significant negJ=987093 -- statistically identical to the buggy
11997838.  VERIFIED: the direction-aware lagged terms are not the
iteration-1 driver; on a cold sweep the lagged stores are zero, so
the wrong-direction reads cannot dominate the first blow-up.

Re-examination of the 42dc0e9 diff showed the actual regression.
The old repair wrote `XM(I)=ABS(XM(I))/10` at a significant-negJ
cell -- a small, local value.  42dc0e9 replaced it with
`XM(I)=XM_FB(I)`, where XM_FB preferred the lagged J_NEXT, then
fell back to `RSQJNU_PREV` -- the bluer-frequency field.  On a cold
sweep J_NEXT_FRQ is all zeros, so every failed cell imported its
bluer neighbour's value; the repairs then amplified R^2*J upward
along the frequency sweep.  That is the 2.6E10 -> 2.7E20 it1 jump.

Fix (commit cd45cf2): XM_FB now uses J_NEXT_FRQ(k) where positive,
else the local SOURCE*R^2 estimate.  The bluer-frequency fallback
was removed; a failed cell restores its same-frequency lagged
value or the local LTE estimate -- both bounded, neither importing
a neighbour's field.  Direction-aware H/J upwinding retained.

Also capped the negJ per-cell error dump at 20k rows/sweep
(commit 325d9a1): the diagnostic wrote ~3 lines per rejected cell
and reached 37 GB inside one iteration, turning the run into a
pure-I/O job.  Counts remain exact; only the file writes stop.

Runtime verification (interactive jobs 11998572 MOM / 11998353 LAM,
cold start, 8 CPU):

MOM path (USE_LAM_ES=F, SOL_METH=TRI):
  it1: L=2.64E10, max%=4.2E14, negJ=170857, nbad50=31.2M
  it2: L=1.96E10, max%=1.5E11, negJ= 55282, nbad50= 7.6M
  it3: L=1.72E10, max%=4.6E02, negJ= 30233, nbad50= 4.2M
  refresh still unconverged (best residual 5.9E3, contracting);
  retained-field rollback holds so an unconverged refresh keeps
  the converging sweep field.  MOMJNXT_GS VER-4 checkpoint writing.

LAM path (USE_LAM_ES=T, SOL_METH=DIAG):
  it1: L=1.98E10, max%=5.9E16
  it2: L=3.62E17, max%=1.4E11
  it3: L=3.91E17, max%=8.6E02, refresh unconverged cells=0
  The it2 luminosity jump is the LTE/T_gas setup transient, same
  as the stock control -- not a solver defect; the population
  correction contracts ~7 decades in 3 iterations.

Status: both KNOP solver paths now converge from a cold start.
The remaining open item is the MOM refresh lagged-map residual
(still ~6E3, slowly declining); the sweep field converges
regardless via the rollback.  Watch subsequent iterations for the
residual trend and for BA-gate opening.

## 9. Equation-level audit vs Knop Papers I+II (29-Sep-2026)

Full report: [[KNOP_MATH_AUDIT]].  Summary:

- FG path (USE_LAM_ES=T): Paper I discretization VERIFIED correct in
  nu-form (xi=1: chi_hat = chi+|GAM|/DLN, S' = eta/chi_hat + Q*I_up
  with sign-based upwind and forward-spacing DLN_N_SAV on the
  compressing branch; |Q| keeps the source contribution positive on
  both branches).  Diagonal OAB Lambda* carries the exact
  dS'/deta = 1/chi_hat, dS'/dchi = -S'/chi_hat sensitivities
  (direction-split S'_P/S'_M), matching Paper II's Lambda*
  construction.  Paper II's wavelength BVP is converged by
  FG_KNOP_REFRESH_V1's alternating-direction GS on frozen opacities,
  whose SRC_PRIME/BUILD_FREQ replicate the in-sweep discretization --
  verified bit-for-bit consistent.  Refresh + maturity gates
  reproduce the paper's "converged BVP -> source update" ordering at
  the great-iteration level.
- MOM path (USE_LAM_ES=F): moment-analogue of the Knop scheme
  VERIFIED after 3ed8f3b/cd45cf2 -- per-equation coupling signs
  (Q_NM at depth, GAMH_SIGNED at midpoints, A_BND at the boundary
  row), |GAM| in chi_hat, explicit weight = implicit weight
  (non-expansive lagged map), direction-aware JUP/HUP in all four
  read sites, same-frequency repair floor.  All RSQ*_PREV (bluer)
  reads verified on expanding branches only.
- Deviations from the papers are documented, not bugs: the
  aberration/advection term (4aI in lambda-form) is absent because
  CMFGEN's nu-form equation lacks it (stock V13 identical); diagonal
  Lambda* replaces the paper's full-bandwidth operator (fixed point
  identical, convergence machinery differs); nu-form log-spaced p
  coefficients absorb the lambda-ratio asymmetry (bounded, same
  convention as V13).
- No remaining mathematical defect found in the discretization.
  Every verified runtime failure traced to plumbing (bluer-field
  repair fallback, DLN_C OMP race, diagnostic I/O flood) or to
  setup physics (LTE start at T_gas), not to the Knop equations.
- Open items: red-edge zero-incident wavelength BC absorbs without
  re-injection at k=NCF compressing cells (boundary-condition
  choice); MOM refresh lagged-map residual still contracts without
  converging on this model (sweep field unaffected via rollback).

## 10. Non-convergence mechanism audit (27-Sep-2026)

Areas where the current code/config can fail to converge, ranked by
demonstrated impact.  VERIFIED items are backed by the m15_2p82d
run history.

### 10.1 VERIFIED live failures on m15_2p82d

- MOM gate permanently frozen (~200 consecutive iterations).  With
  the default flags (MOM_GATE_BULK=F, MOM_RF_KEEPNEG=F -- neither is
  set in the m15 MOM VADAT), MOM_UNSAFE_B/MOM_FREEZE_B trip on
  MOM_NEGJ>0, and the residual negJ population is quasi-steady
  (~300-550/sweep), never reaching 0.  BA never engages; forced
  LAMBDA runs the whole run.  solve_for_pops.f:312-315.
- Forced-LAMBDA does not contract on this model: max % increase
  oscillates 4 -> 27 % at depth 62-63 over iterations 181-200 with
  no declining trend; luminosity asymptotes ~1.26e10.  The
  protective mode perturbs populations on an immature field
  (cycle-audit issue 1).
- MOM_KNOP_REFRESH reports "unconverged -- best residual=1.0000E+00"
  every iteration: the lagged map is marginally contractive at
  best on this model's compression strengths, so J_NEXT only gets
  one-bin-per-iteration propagation inside compressing zones.
- LAM long-horizon divergence: the Sep-26 LAM history shows
  luminosity 1.98e10 (it1) -> 1.27e33 (it~190), max % ~1e100-1e200,
  frozen in forced-LAMBDA.  The capped I_NEXT store (1e20) feeds
  |Q|*I_NEXT source terms at huge magnitude once triggered -- the
  cap is a numerical guard that corrupts the fixed point after
  activation (it~170, ~2e8 capped cells/sweep).

### 10.2 Structural non-convergence areas (code-level)

- Jacobi-lagged compressing coupling has a marginal/expansive
  subspace: at near-pure-scattering cells (1-COH -> 0) the row
  diagonal approaches HL+HU+PSI so the lagged-map gain -> 1; the
  EPS frequency-lagged term (+/-2*EPS_PREV*J_up scattered into
  rows) has no matching diagonal absorption at the same frequency,
  so the convex-combination bound that protects the PSI channel
  does not cover it.  Refresh abort + best-pass restore contain
  the damage; they cannot make the map contract.
- negJ repair values reach J_NEXT_FRQ when MOM_RF_KEEPNEG=F:
  ~300-550 artificial repair values enter the lagged store each
  sweep, reseeding the lagged read with noise -- keeps sweep drift
  above gate thresholds and wobbles the SE input.  Config-level
  contributor to 10.1.
- Quarantined cells (KEEPNEG=T) hold stale values permanently and
  are excluded from the metrics: a growing quarantine hides drift
  and feeds stale J into neighbours' explicit terms.
- Eddington f-loop silent bail (comp_j_blank.f:414-423): after 10
  inner iterations INACCURATE is cleared unconditionally; a
  non-converged f is used for that frequency's J.  Currently benign
  (max residual ~5e-4) but a systematic channel if f stalls higher.
- Non-M-matrix tridiagonal + DGTSV: centered EPS stencil gives
  non-diagonal-dominant rows (observed TB=0.004 vs |TA|+|TC|=0.055);
  LAPACK solves them exactly but the result is unphysical (negative
  J) -> repair path -> noise injection (10.2 item 2).
- Boundary row branch switching: A_BND sign (HBC floored at 0) can
  flip across sweeps, alternating the row between expanding and
  compressing upwind -> oscillating outer flux.

### 10.3 Iteration-control areas

- Gate thresholds under MOM_GATE_BULK=T are lenient: freeze opens
  while ~10 % of cells still drift >5 %; BA kicks then re-trip the
  unsafe bound -> Lambda<->BA alternation with growing kick size
  (TypeIIn: 3.4e3 -> 2.7e4, increasing).  Non-contracting cycle
  risk.
- DO_LAM_AUTO NUM_LAM=2: two Lambda iterations cannot re-mature the
  field after a large BA kick -> every BA lands on an immature
  field.
- Warm-restart maturity loss: resubmits replay the 4-iteration
  warm-up and re-earn maturity (~25 its on TypeIIn); combined with
  the 500-iteration cap, repeated resubmits burn the budget.
- MAX_CHNG_LIM=1e100 disables the divergence guard (config);
  MAXCH masking via DO_LEVEL_CHK lets up to 10 levels still move
  >10 % while MAXCH reports convergence (upstream semantics).
- INEXT_RELAX/JNEXT_RELAX=0.5 damps the update but cannot change
  the fixed point; on a marginally-contractive (gain->1) map it
  only halves the noise amplitude.

### 10.4 Config state of the current submissions

m15_2p82d_knop_MOM VADAT sets neither MOM_GATE_BULK nor
MOM_RF_KEEPNEG: it carries the default (keepneg off, negJ>0 gate)
that produced the 200-iteration freeze in 10.1.  The TypeIIn
winning config (MOM_RF_KEEPNEG=T, MOM_GATE_BULK=T) addresses the
gate-freeze + noise-reseed mechanisms; whether it converges the m15
structure is untested.  Jobs 12000120/12000122 were cancelled at
~15-24 min on 27-Sep and their `_v2` run directories deleted; the
queue is currently empty.

## 11. External review triage (30-Sep-2026)

An external review of an uploaded source snapshot (knop_solver_files/)
was checked claim-by-claim against the current tree.

### 11.1 CONTRADICTED

- "dLOG_NU undeclared in COMP_J_BLANK (compile blocker)".  The
  snapshot's comp_j_blank.f lacked the USE line; the tree version
  imports RADIATION_MOD, which declares dLOG_NU
  (new_main/mod_subs/radiation_mod.f:120).  Direct ifort compile of
  comp_j_blank.f passes; the shipped binary has always linked.

### 11.2 VERIFIED and fixed in this commit

- Unbounded nonmono insertion.  MOM_J_CMF_NONMONO_V1 and
  FG_J_CMF_NONMONO_V1 both inserted INT(|dV|/T2) sub-points per
  interval with no ceiling: a single shock interval could request
  10^4+ points and, for FG, overflow the integer counter.  New caps:
  KN_NM_MAX_INS=256 (MOM) and KN_RAY_MAX_INS=256 (FG ray), applied
  identically in the count pass and the fill pass so K stays
  consistent.  MOM logs the capped-interval count; FG prints
  NRAY_MAX.
- T2<=0 guard.  T2=VDOP_FRAC*MINVAL(ABS(VDOP_VEC)) divides by zero on
  an all-zero VDOP (MOM); FG's MINVAL(VDOP_VEC) without ABS turns a
  negative VDOP entry into silent under-resolution.  Both now stop
  with an explicit message when T2<=0.
- SIGMA_TRUE silent zero.  recompute_sigma_true.f set
  SIGMA_TRUE=0 whenever V<=1e-10 or |dR|<=1e-10, marking a
  stagnation/degenerate cell as non-compressing with no record.
  Interior cells with a usable dR now regularize the 1/V pole with
  the smaller neighbour velocity (preserves the sign and the large
  magnitude appropriate to a stagnation cell); genuinely degenerate
  cells (no dR, or all neighbour V~0) still get 0 but are counted
  (N_UNDEF) and reported once per call.
- FG refresh had only a per-thread memory cap (KN_RAY_PACK_MB
  bounded one thread's coefficient pack; 16 threads could still
  reach ~tens of GB).  New KN_RF_TOT_MB=32768 bounds the aggregate:
  N_TH is scaled down by 80 B x MAXVAL(NI_RAY) x KN_NCF per thread.
  MOM refresh is serial -- not affected.
- FG refresh call site is now gated explicitly
  (DO_NONMONO_KNOP .AND. DO_KNOP_REFRESH), matching the MOM call;
  the routine's internal guard stays as a second layer.
- Grid-quality diagnostic: MOM prints ND_SM->ND_NM, T2 and the
  capped-interval count whenever the expanded grid is (re)built.

### 11.3 VERIFIED, deferred

- OBS-frame checkpoint restore allocates the whole
  I_F(N_COMP,NP,NCF_F) record transiently.  Chunked reads need the
  writer to emit per-ray records, i.e. a checkpoint-format bump;
  deferred until the current validation jobs finish so their
  FGINEXT_GS files stay readable.
- Convergence constants (gate tolerances, floors, runaway factors)
  remain distributed across routines; a shared policy block is a
  follow-up, not done here to keep the diff reviewable.

### 11.4 Second review (same-day, second reviewer)

Checked against the tree; overlaps and new items.

VERIFIED and fixed in this commit:

- Duplicate SIGMA_TRUE in obs/cmf_flux_v5.f: the EW-estimator path
  reimplemented the unclamped gradient inline with the old silent-zero
  logic.  Now calls the shared worker SIGMA_TRUE_CORE (extracted from
  RECOMPUTE_SIGMA_TRUE so it carries the stagnation regularization and
  N_UNDEF).  It must not call RECOMPUTE_SIGMA_TRUE itself: that routine
  gates on CONTROL_VARIABLE_MOD flags, which are never set inside
  cmf_flux.exe -- only the CMF_FLUX_CNTRL_VAR_MOD copies are.
- Stagnation branch caught V<0 too.  Primary test changed to
  |V|>1e-10 at all three sites (interior + both boundaries); negative
  velocities now flow through the algebraic (R/V)dV/dR-1 form, which
  is defined for V<0 and carries the infall sign into GAM.
- FG compressing mask missed V<0 cells: COMP_DEPTH_IDX and
  HAS_COMPRESSING keyed on SIGMA_TRUE<-1 only, but GAM=(V/R)(1+sigma
  mu^2) is negative for V<0 at every mu.  Both flagging sites now
  include V<0, and the rebuild consistency check uses the same
  criterion.  (MOM was already safe: its HAS_COMPRESSING keys off the
  signed GAM, which turns negative under infall on its own.)
- Stale comment: FREQ_IDX_ITER described as 1..NCF; counter is
  0..NCF-1 (slot s = counter+1).  Comment corrected.
- GAM_NM holds the physical signed gamma for one block, then is
  overwritten with the dimensionless matrix coefficient.  Works via
  GAM_SIGNED_NM, but fragile; the overwrite now carries a comment
  marking it deliberate.

ALREADY DOCUMENTED here (no new action):

- Non-M-matrix tridiagonal rows in compressing zones + DGTSV
  detection/fallback -- section 10.2.  Reviewer is right that this is
  detection-and-contain, not a provably well-posed discretization; a
  real fix is a scheme change, not a patch.
- GS fixed-point refresh not reliably convergent (stall/abort keeps
  pre-refresh field) -- section 10.2, marginal lagged map (gain -> 1
  at near-pure-scattering cells).

PARTIALLY VALID:

- "Two lagged stores, no cross-check": FG I_NEXT and MOM J_NEXT never
  coexist in one run (USE_LAM_ES picks one solver), so they cannot
  disagree within a FINAL.  The real gap is downstream: the OBS stage
  always runs the FG solver, so a MOM-mode model's OBS spectrum
  cold-starts its I_NEXT coupling (MOMJNXT_GS cannot seed it).
  Second-order for the current validation; a J/H->I mapping is a
  development task, not a patch.

### 11.5 Third review -- OUTGEN runtime log (job 12000438, MOM)

Reviewer read 55 iterations of the MOM OUTGEN.  Claim-by-claim:

1. "Luminosity never converges": OVERSTATED.  L declines 2.64e10 ->
   1.32e10, but the tail is near-flat (1.3204 -> 1.3199e10 over
   it53-55, ~0.04%/it).  The early decay is the LTE-transient plus
   population settling under forced LAMBDA, not a moving radiation
   field -- FG maxchg is 1.8e-3 and the retained MOM field is stable.
   Still immature: gate blocked by nbad50 (~50k vs 12k needed) and a
   single flickering cell (sweep maxchg ~2.5-7k).

2. "MOM_KNOP_REFRESH never converges": VERIFIED.  All 55 iterations
   exit on the stall detector (two same-direction ratios >0.90),
   residual floor ~2.5e3-4e3 oscillating, pre-refresh field retained
   every time.  Consistent with the documented marginal lagged map
   (section 10.2) -- the iteration is not contractive on this model.
   Containment works as designed: retained-field delta = 0.

3. "Chronic weak row at expanded (21-24, f~1714) = the sigma=-120
   cell": PARTIALLY WRONG.  The weak row is real and persistent
   (min TB/(|TA|+|TC|) ~0.17-0.27, TB<=0 count stays 0), but
   MOMJNXT_GS checkpoint decode (VER=4, R_PNT_NM) maps expanded
   d=21-24 to SOURCE d=20-23 -- the outer slow-ejecta plateau
   (R~6.0e14 cm, V~55-58 km/s).  The sigma=-120 cell is source d=66,
   which occupies expanded ~93-212 (the shock ramp absorbs most of
   the 142 inserted points).  Two distinct trouble spots, not one.
   Frequency ~1714 of 56853 = blue end.  Why the weakest-conditioned
   row sits in low-V outer cells is unresolved -- candidate: at small
   V the EPS upwind term dominates TB over the spatial/absorption
   terms.

4. "Blank physics-diagnostics blocks on ~half the iterations":
   VERIFIED -- reporting artifact, but it exposed two real bugs.
   Mechanism: SCR_RITE_V2/DUMP_MOM_JNEXT_GS runs AFTER the refresh.
   The in-routine accumulator reset at FREQ_INDX=1 sat AFTER the
   solve, so (a) every call erased frequency-1's own statistics --
   systematic negJ/matrix undercount -- and (b) when the refresh's
   last pass was descending (ends on FREQ_INDX=1), the reset wiped
   the whole pass's physics accumulators before the dump.  The
   ~50/50 blank/populated split is the stall-exit parity.  The
   solves themselves ran every iteration (per-pass residuals are
   real).  FIX: reset moved ahead of the first accumulator write and
   gated on .NOT.MOM_KN_REFRESHING; the refresh driver now resets
   the physics stats per pass so the dump shows the last attempted
   pass.  Rebuilt; takes effect on next requeue.

5. "Iteration 40 excursion, never fully recovers": VERIFIED.  it40
   (643.8 s vs ~136 s) is the one iteration where the maturity gate
   entered the BA branch, printed the KNOP ALO dJ/dCHI diagnostic
   ([-4e-42, 6.5e3]), got blocked, and forced LAMBDA with a capped
   +/-7 dex population swing at depth 52 (outer shelf, V~150 km/s --
   third trouble spot).  L stepped 1.380 -> 1.343e10 at it41 and
   nbad5 jumped 2-3x for ~3 iterations before relaxing.  The
   forced-LAMBDA perturbation moved the solution to a different
   unconverged branch rather than settling it.

REVIEWER'S PRIORITY CONFIRMED: the two questions worth chasing are
(a) whether the last-pass stats (now correctly reported) show the
weak row migrating or fixed -- will read from the next requeue's
OUTGEN -- and (b) why source d~20-23 produces the weakest diagonal
dominance, which needs the TA/TB/TC coefficient decomposition.

### 11.6 Fourth review -- gate feedback-loop mechanism (this commit)

Reviewer proposed the end-to-end chain: MOM refresh never converges
-> sweep nbad50 stays over NTOT/1000 -> MOM_UNSAFE_B -> GS_FREEZE ->
forced LAMBDA every iteration -> order-of-magnitude slower
convergence.  Verification:

VERIFIED, with corrections:

- The gate chain is real and confirmed in solve_for_pops.f:285-316.
  But TWO constraints bind, not one: even after nbad50 clears
  (42.5k -> 12k needed, ~3%/it decay -> ~it100), MOM_SW_MAXCHG>10
  (flicker at ~2.5-7k, no decay trend) keeps UNSAFE set alone.
- "One pathological cell": partially wrong -- nbad50 is a ~42k-cell
  population, and the refresh argmax migrates across expanded
  d=100-168 = source 59-65 (the shock ramp), not a fixed cell.
- The (21-24, f~1714) min-TB row is a DARK cell: J~1E-20, ~1E-15 of
  the column peak (MOMJNXT_GS decode).  Real weak conditioning,
  physically irrelevant to the field -- it only pins the diagnostic.
- The oscillating cells carry REAL field: argmax cells sit at
  J/colmax = 5-86%.  The NBAD50 metric's 1E-4 column-peak floor
  already guards tiny-J inflation -- the non-self-consistency is
  genuine, not a metric artifact.
- sigma=-120 at source 66 comes straight from the 70-point hydro
  profile (V=4233->6593 km/s across d 65-66); not a MON_INT
  interpolation artifact.

FIX IMPLEMENTED -- per-cell adaptive damping in MOM_KNOP_REFRESH
(reviewer's option c, the mechanism-level fix):

  Damped fixed-point iteration x+r(M-x) has multiplier
  1-r*(1-lambda), which contracts for r < 2/(1-lambda) -- any
  oscillatory mode (lambda<1) is contractible with small enough r;
  divergent modes (lambda>1) cannot be damped but none appear in
  this job (the 1000x abort never fires).  Global r=0.5 leaves
  lambda~-3 cells marginal, matching the observed stall.

  - MOM_DAMP_NM(ND_NM,NCF): persistent per-cell relax schedule,
    init JNEXT_RELAX=0.5, floor MOM_DAMP_MIN=0.02 (covers lambda to
    ~-100).  Applied only while MOM_KN_REFRESHING -- sweep writes
    keep JNEXT_RELAX unchanged.
  - Cycle detection per pass k>=3, same parity (k vs k-2 isolates
    the composite map): cell cycling iff rel step > 0.5 (the NBAD50
    floor) AND |J_k - J_{k-2}| < 0.3*|J_k - J_{k-1}|.  Flagged cells
    halve their factor; settled cells recover x1.25 toward 0.5.
  - Same factor applied to the H_NEXT write (lagged partner field).
  - Pass report gains "damped cells / min relax" when nonzero.
  - MOM_DAMP_NM is not checkpointed; a requeue re-derives it within
    ~3 passes.  Two snapshot buffers (+~190 MB) are refresh-local.

  Unlock chain: damped refresh reaches the fixed point -> retained
  J_NEXT becomes self-consistent -> next sweep's fresh-vs-stored
  drift collapses -> nbad50/maxchg fall -> GS_FREEZE releases -> BA
  resumes.  First converged refresh may still gate once: moving the
  stored field to the true fixed point is a real jump and correctly
  holds LAMBDA for that iteration.

REVIEWER OPTION (b) -- exempting chronic cells from NBAD50 --
REJECTED: it would open the gate while the stored field is
genuinely non-self-consistent at significant-J cells; BA against a
field with ~50% coupling error is the masking pattern AGENTS.md
forbids.  The gate is measuring the right thing; the map must be
fixed, not the metric.

REVIEWER OPTION (d) -- is BA needed: moot if the freeze never lifts;
under pure LAMBDA the run converges linearly at best.  L=1.319e10
is stable but population maturity is the binding constraint.

### 11.7 Fifth review -- overall non-monotonic robustness assessment

Reviewer's question: is the solver robust for non-monotonic velocity
fields.  Most items describe the pre-f70d677 snapshot; status here:

VERIFIED (and already fixed in this tree):

- "cmf_flux_v5 duplicates SIGMA_TRUE, drops the fallback": fixed in
  f70d677 -- now calls SIGMA_TRUE_CORE.
- "Stagnation fallback over-broad for V<0": fixed in f70d677 --
  |V|>tol test at all three sites.
- "FG flags compressing cells by SIGMA_TRUE<-1 only, misses V<0":
  fixed in f70d677 -- both FG flagging sites now include V<0
  (GAM<0 for V<0 at every mu).  MOM always keyed on signed GAM.
- "MOM refresh non-contractive at strong compression": mechanism
  confirmed and now addressed by per-cell adaptive damping, 2e3e563.

VERIFIED (still open):

- No runtime notice when V<0 is present.  FIXED this commit:
  SIGMA_TRUE_CORE warns once per process -- the signed path is
  algebraically consistent (flagging + sign-propagating GAM) but
  has never been exercised by a validation run.
- Multi-root resonance counts (N_RAYS_MULTI, NRES_MIN/MAX) are
  printed once as diagnostics and unused by the grid builder --
  no adaptive refinement where resonances cluster.  Design gap,
  logged; adaptive gridding is a feature, not a patch.
- obs_frame_sub_v9.f bidirectional bracket hunt: VERIFIED present
  (line 777, "second hunt direction ... non-monotonic velocity
  fields") -- the reviewer's reading of it as a real prior fix is
  correct.
- GAM=0 turning points get no artificial floor in the ray solver --
  verified at fg_j_cmf_nonmono_v1.f:1775.

SCOPE CONFIRMED: "non-monotonic" throughout means dV/dr sign change
(SIGMA_TRUE < -1), not V<0 infall.  The former is built and tested;
the latter is now sign-consistent but explicitly warned as untested.

### 11.8 Sixth review -- stall-point audit (29/30-Sep-2026)

A static review of the DO_NONMONO_KNOP stall paths
(FG_J_CMF_NONMONO_V1 / MOM_J_CMF_NONMONO_V1 / SOLVE_FOR_POPS chain)
flagged six hazards plus two blocked items.  Verified against the
live tree:

VERIFIED and fixed (commits a9ab6b7, be1213e, 72641f2, d57d105,
5e7d5cd; shell watchdog ed48071 in the main repo):

- SIGMA_TRUE flapping -> grid rebuild every iteration.  Confirmed,
  and worse than the review stated: on NEW_R_GRID the
  I_P/I_M_DEPTH_NEXT store is reallocated ZEROED (the JNEXT_OLD
  remap only exists on the MOM side), and a mid-run rebuild can
  never match the FGINEXT_GS tag anyway.  Fix: exit-side hysteresis
  on the compression-set test -- depths enter on the unchanged
  sigma < -1 / V < 0 test but leave only on sigma > -0.98 or
  V > +1 km/s.  The stored COMP_DEPTH_IDX mask still uses plain
  thresholds, so the checkpoint fingerprint is unaffected.
- Unguarded DO WHILE bracket searches (J_PNT, both H_PNT loops;
  J_INDX_NM, H_INDX_NM).  No in-loop bound check: a depth bracketing
  no interval marched K past the array end -- a silent hang in an
  unchecked build.  All five now STOP with depth/ray/endpoint
  diagnostics.
- "Permanent" cell quarantine under MOM_RF_KEEPNEG.  Partially
  confirmed: XM_BAD_NM is rebuilt every call, so the cell IS
  re-solved -- but a structurally failing cell retries
  deterministically forever holding a stale iterate.  Fixed as a
  visibility gap: MOM_QUAR_AGE counts consecutive quarantined
  sweeps per cell and MOM_GS_NPQUAR reports cells at age >=3
  (getter, MOM bulk line, DUMP_MOM_JNEXT_GS report).
- GS_FREEZE with no escape valve.  Confirmed (explicitly documented
  as intentional).  Fix is operational, per the reviewer's own
  recommendation: SOLVE_FOR_POPS writes GS_STATE every iteration
  under DO_NONMONO_KNOP; the generated slurms refuse to requeue or
  report completion when the run ends a segment >= KNOP_STALL_MAX
  (default 60) consecutive frozen iterations, touching
  .CMFGEN_KNOP_STALLED and exiting 7.  VADAT key GS_FREEZE_LIM
  (default 0 = off) adds an in-model stop for non-pipeline runs.
- BA trust radius floor: SOLVEBA_V13 does floor BA_TRUST_CUR at 1.0
  (the blocked item resolves benign), but BA_NREJ_CONSEC was
  uncapped -- a LAMBDA/rejected-BA alternation could run forever at
  the floor.  New VADAT key BA_NREJ_LIM (default 20, <=0 disables)
  stops the model loudly.

VERIFIED, no code change needed:

- MOM_KNOP_REFRESH_V1 "burns full budget then reverts": bounded --
  KN_M_MAXPASS=60, divergence abort at 1E3x min-residual,
  same-direction stall exit, bulk-drift exit (MOM_RF_BULK_EXIT),
  conditional best-pass retention vs the pre-refresh field.
- FG_KNOP_REFRESH_V1: per-ray bounded via KN_RAY_MAXPASS in
  KNOP_RAY_SOLVE_V1.
- Checkpoint fingerprint strictness: bitwise R/V/SIGMA +
  COMP_DEPTH_IDX match is deliberate -- the mask-flap trigger is
  removed by the hysteresis fix; remaining invalidations on genuine
  grid changes are correct behavior.

Follow-up (30-Sep-2026, commit 7b07fdf): MOM_TRIAL_MIN marginal
release -- see the dated entry below.


== 2026-09-27 runtime triage (job 12002633, ratcheted damping) ==

With the ratcheted schedule live, it20-45 showed: damped-cell pool
collapsed 27k -> ~500 as the field settled; min relax deepened to
0.165 then pinned -- flagged cells rotate rather than persist.
Refresh best residual stayed flat at 2.5-6.6E3 and the argmax
wandered across expanded d~100-168 (source d 59-65, shock ramp).
~70-150k cells moved >50% per sweep while only ~500-900 met the
period-2 signature -- the dominant marginal mode is same-direction
drift (lambda -> +1), which under-relaxation cannot accelerate.

FIX (this commit): per-cell Aitken-Delta2 extrapolation in the
refresh driver.  For cells that are NOT cyclers and NOT damped, the
same-parity secant LAMC = (x_k - x_{k-2})/(x_{k-2} - x_{k-4})
estimates the composite-map multiplier; for 0 < LAMC < 1 the
geometric tail is extrapolated, x += D1*LAMC/(1-LAMC), with LAMC
capped at 0.9 (correction <= 9x the same-parity step).  Sign flips
(LAMC<0) and divergent drift (LAMC>=1) are skipped -- oscillation
belongs to the damp schedule and divergence to the abort guard.
History chain extended to JPREV1..4 (x_{k-1}..x_{k-4}); +190 MB
refresh-local.  Reported per pass as "accel cells=".

If the drift diagnosis is right, watch for: accel cells >> damped
cells once past pass 5, best residual dropping below the ~3E3
floor, then refresh convergence -> retained-field delta > 0 -> the
gate's sweep-drift metric collapsing on the following iterations.


== 2026-09-28 vector Anderson acceleration on the composite map ==

Motivation: scalar Aitken engaged only ~0-30 cells/pass in job
12003796 -- the marginal mode is collective, not per-cell.  The
refresh map is affine at fixed ETA/CHI (the moment solve is linear
in the lagged coupling), so Anderson acceleration on the composite
(2-pass, same-parity) map is equivalent to GMRES on (I-A)x=b and
does not require contraction.

Implementation in the refresh driver (mom_j_cmf_nonmono_v1.f):
- XHIST_AA/FHIST_AA hold the last MOM_AA_M+1=4 composite iterates
  X_i (field at even-pass boundary, = JPREV2 pre-shuffle) and
  residuals f_i = G(X_i)-X_i.
- At even passes >= 4 the pair is pushed; once a second pair exists
  a min-norm least squares over secants dF_j = f_i - f_{i-j} gives
  gamma via DGESV, and x_new = G_i - Sum gamma_j (G_i - G_{i-j}).
- Safeguards: correction capped at 10*||f_i|| (global L2), per-cell
  fallback to the raw map output when the accelerated value would go
  non-positive, Gram solve failure skips the step entirely, and the
  1E3x divergence abort still bounds the sequence.
- The scalar Aitken block is bypassed while AA is on
  (.NOT. MOM_RF_AA in the drift branch) -- one extrapolation owner.
  Cycle-detection damping still runs: AA needs bounded raw iterates.
- Module switch MOM_RF_AA (default T) disables the whole layer.
- +770 MB refresh-local (8 field buffers of ND_NM x KN_M_NCF).

Watch for: "AA m=" lines at even passes, ||f|| shrinking across
boundaries, PASS_CHG dropping below the ~600 floor, then refresh
convergence and retained-field delta > 0.


== 2026-11-05 LAM blow-up at it127-131 and the store-drift gate fix ==

Observed (job 12002113, m15_2p82d_knop_LAM_FreeNe): healthy run to
~it123 (L~8e11, cap/bad writes=0), then
  it126: refresh residual 1.8e-3 < 2e-3 tol -> gate mature -> BA step
         (pop changes only ~50-70%) BUT store cells>50% = 2.48M
  it127: "BA matrix **NOT** computed" -- stale-matrix linearized SE
         update -> SOLVEBA returned max change 1.1e12 at depth 8;
         store cells>50% = 3.09M, refresh still converged (1.8e-3)
  it128: poisoned populations -> transfer runaway; clamp events 3221,
         L = 1.8e27
  it129-131: refresh residual 28.6 -> 47.7, unconverged cells ~5k,
         clamp events -> 47M, store cells>50% = 57M-66M,
         L -> 9.8e30 -> 1.6e32.  Gate slammed back to LAMBDA but the
         state was already diverged; per-ray NDIV rollback protects
         only rays whose refresh diverged, not a poisoned opacity
         field.  Job cancelled (checkpoint retained for post-mortem).

Root cause: GS_FG_IMMATURE in SOLVE_FOR_POPS measured only the
refresh residual (GS_MAXCHG) and refresh unconverged-cell count
(GS_NBAD).  FG_GS_NCONV -- the count of store cells the refresh had
to move >50% to reach its fixed point (or in-sweep >50% writes when
no refresh ran) -- was computed and printed ("store cells>50%=")
but never fed to the gate.  For ~10 iterations before the blowup
NCONV sat at ~2.5M (0.33% of 758M significant cells) while the
residual stayed below tolerance: the lagged field's fixed point
existed but kept moving, i.e. the field was still drifting between
great iterations.  A BA step against that field -- and worse, the
next iteration's *stale* BA matrix on a moved field -- injected
unbounded SE corrections at depths 8-60 (the compression ramp).

Fix (this commit):
- GET_FG_GS_MATURITY now also returns NCONV.
- GS_FG_IMMATURE gains the term: refresh ran (GS_NTOT>0) AND
  GS_NCONV > max(50, GS_NTOT/1000).  Same 0.1% bulk rule MOM uses
  via MOM_SW_NBAD50 / MOM_RF_N50.  Applied as an OR after the
  residual/bulk branch so a converged residual cannot mask a
  churning store.
- Immature-state diagnostic now prints nconv50.

Not fixed here (noted for follow-up): SOLVEBA_V13's "BA not
computed" path returned a 1.1e12 correction despite MAX_LIN_COR;
the LAMBDA/DIAG path carries the LIMIT option, the linearized
stale-matrix path apparently does not apply an equivalent bound.
The gate fix prevents reaching that state; an independent bound on
the returned correction is a reasonable second layer.

Reference consulted: Ceulemans, De Ceuster & Decin 2025 (JQSRT,
arXiv:2505.08309) -- perturbation analysis shows the Baron-
Hauschildt/Knop CMF discretization develops oscillatory eigenvalues
(-1 <~ lambda < 0) when the per-cell Doppler shift is large vs the
frequency spacing; their "frequency matching" stencil is a
discretization-level fix.  This is consistent with our reading that
the sigma=-120 compression zone produces marginal modes no lagged
iteration can contract -- iterative acceleration (AA) treats the
symptom, cell refinement or stencil rematching treats the cause.


== 2026-09-30 marginal-release BA trial + freeze watchdog (7b07fdf, 5e7d5cd, ed48071) ==

Two layers added for the "frozen forever" failure mode identified in
the sixth external review and observed live on the FreeNe M15 run
(>150 consecutive forced-LAMBDA iterations on ~0.12% drifting cells
that asymptote but never converge under the one-bin upwind):

- MOM gate split into tiers (solve_for_pops.f): MOM_UNSAFE_B is now
  the catastrophic backstop only (any capped write, >1% of cells
  drifting >50%, a cell drifting >10x, >10% negJ/quarantined); the
  strict maturity terms stay in MOM_FREEZE_B as a marginal tier that
  blocks BA via the release check rather than an unconditional veto.
- MOM_TRIAL_MIN (VADAT, default 30, <=0 off): after that many
  consecutive frozen iterations, with DO_BA_RES_ACC armed and the FG
  side clean, the iteration is released to a residual-accepted BA
  step -- residual growth rolls the step back and demotes to LAMBDA
  next call, so the trial is self-throttling.  N_FROZEN keeps
  counting through the trial, so the stall bounds still apply.
- GS_STATE is written every SOLVE_FOR_POPS call under
  DO_NONMONO_KNOP; the generated slurms (FINAL + SINGLE_FINAL) read
  n_frozen after batch.sh and fail the job (exit 7,
  .CMFGEN_KNOP_STALLED marker) at KNOP_STALL_MAX=60 consecutive
  frozen iterations, whether the segment ended by graceful stop or
  completed.  GS_STATE/marker are wiped on full resubmission.
- GS_FREEZE_LIM (VADAT, default 0 = off) provides the same stop
  inside the model for runs not driven by the generated slurm.

Not yet runtime-validated: the marginal release has not taken a BA
step on a real model; watch for the 'MOM field marginal' line on
the first released iteration and confirm the resacc rollback fires
if the step grows the SE residual.

## 12. Trial machinery runtime-validated + freeze-clock fix (30-Sep-2026, evening)

Section 11's "not yet runtime-validated" items now have runtime
evidence, and the validation exposed a design bug that made the
marginal tier self-defeating.

### 12.1 VERIFIED: trials release and get accepted

Job 12040589 (Dessart m15_2p82d dif, USE_LAM_ES=F, old binary): the
trial fired at it30 and it32; resacc accepted both steps (nrej=0).
Post-release the field churned (n50 9.8e-4 -> 3.2e-2) but re-matured
below the trial gate in 1-2 iterations.  lum_ratio moved
0.16 -> 1.03 -- the held-pops stall was broken by real BA progress,
exactly the nested-iteration behavior the trial tier was built for.

### 12.2 VERIFIED + FIXED: un-refunded freeze count kills progressing runs

The same run hit GS_FREEZE_LIM=60 at it~60 and stopped as
non-convergent while pulsing correctly.  Cause: N_FROZEN counted
every frozen iteration including accepted-release ones, and the
raw-mature gate (MOM_SW_MAXCHG <= 10; observed >= ~1.3e3, asymptotes
~1e3 on marginal cells) can never produce the 2 consecutive mature
sweeps that clear it.  Job 12040591 (dbg baseline) failed the same
way.  CONTRADICTED: the section-11 design assumption that
"N_FROZEN keeps counting, so the stall bounds still apply" is
sufficient -- it bounds stalled runs but also bounds healthy ones.

Fix (solve_for_pops.f, exe rebuilt 15:56): an accepted trial
(BA_REJ_THIS=F) refunds N_FROZEN to MOM_TRIAL_MIN -- the marginal
tier must be re-earned, not gifted -- and increments N_TRIAL_ACC.
Rejected trials keep counting.  New VADAT key GS_TRIAL_LIM (default
15, <=0 = old count-forever behavior) vetoes further trials between
raw-mature sweeps so pulsing-forever stays bounded.  n_trial_acc is
written to GS_STATE and adopted on warm restart; a genuine mature
release resets it.  MOM_TRIAL terms prints trlim + ntr=.

### 12.3 Paper-fidelity fixes folded into the same build

mom_j_cmf_nonmono_v1.f: removed a dead midpoint CHI_HAT loop that
read uninitialized/stale GAM_MAG_NM; expanding-row PSIPREV changed
to |GAM*GAM_FAC_PREV| at all three sites so a sign flip of
1+sigma*K_on_J cannot inject a negative upwind source (Knop I
Eq. 9 requires that coefficient >= 0 by construction).

### 12.4 Validation chain + current submission

itest (dessart_m15_2p82d_KNOP_MOM_dif_itest.yaml, interactive,
NUM_ITS=3): job 12042808 COMPLETED 0:0 in ~16 min.  Confirmed the
bordered diagnostic tables, 9-term MOM_TRIAL line, `frozen= T`
(caught and fixed a 4(A,A) -> 3(A,A) descriptor bug that printed
raw bytes), GS_STATE n_frozen-first + n_trial_acc, ITER_HISTORY,
clean RVTJ, no fatal/NaN.

Full model: FINAL 12043642, OBS 12043643 (small, 16 CPUs, 64G,
NUM_ITS=500), fresh cold start via run_model_grid.py on the rebuilt
exe.  (First run under the new CMFGEN_PROD_CPUS=16 default.)  Watch
n_frozen/n_trial_acc in GS_STATE: healthy marginal-tier behavior is
trial fires -> n_frozen refunded to 30 -> re-earns marginal status;
the run is genuinely converging when 2 raw-mature sweeps clear
n_frozen to 0 or the SE/T corrections finish under it.

## 13. JFNK F(0)-gate audit + SOB limit-cycle read (02-Oct-2026, evening)

### 13.1 Why no Newton step ever ran

Every JFNK attempt across cpl/fresh/seq/smoke aborted at the F(0)
identity test: the coupled eval residual at x=0 was compared to the
outer STEQ and aborted when rel>0.5 (observed rel=0.27, 0.82, 44, 1.0).
Audit of the machinery found the eval is *not* broken in the way the
test assumed:

- JFNK_RAD_SWEEP_V1 (cmfgen_sub.f) is a verbatim copy of the
  production continuum loop; diff is only POPS -> POPS_EVAL plus the
  intended eval guards (JFNK_IN_EVAL read-only EDDFACTOR access,
  FREEZE_I_NEXT store-write suppression, COMPUTE_BA held false so BA
  terms do not stack per eval).
- The mismatch is dominated by honest lagged-field physics: the outer
  STEQ accumulated radiative terms while the I_NEXT store evolved
  mid-sweep; the eval re-solves J reading the *updated* store
  (frozen).  rel>0.5 simply measures how far the store moved, and the
  eval operator is the one the next real sweep will actually use.
  ~39% of store cells are compressing on this model, so a moving
  store legitimately moves many residual rows.
- The smoke run's rel=1.0 was a different failure: 6 eval components
  went Inf on a diverged state (L=5.4e29, se_res=3.3e32) and the NaN
  guard overwrote the *entire* residual with 1D20, so the test read
  a masked uniform vector.

Fix (commit e764368): the residual now records jfnk_resid_nbad and
caps only the offending components at +/-1D20; the F(0) test aborts
only on nbad>0 (broken eval) or NaN, adopts the eval residual as the
Newton baseline as before, and reports a spread diagnostic (count of
equations differing >50%) with a warning if the mismatch is global.
Plus the WIP already in the tree: log-x unknown [JFNK_LOGX],
pseudo-time continuation [JFNK_PTC] (Knoll & Keyes), eta floor 3e-4
(Arramy et al. 2024), coupled baseline re-eval, stock-cycle Lambda
cap GS_LAM_CAP while the field is immature, JFNK gated on
GS_FG_IMMATURE, and BA_RES_ACC no longer rolling back on a Lambda
follow-up call.

### 13.2 Live runs (as of ~20:30)

- seq2 (job 12066258, cmfgen_validation_JFNK/m15_2p82d_jfnk_seq2):
  the intended test of the new iteration control.  it1 kept the full
  BA step (stock cycle), it2+ Lambda capped at factor 3 while churn
  was 2.8e-2; FG refresh residual fell below tol by it2.  JFNK will
  not fire until the field matures and JFNK_WARMUP=6 passes.
  Pass criterion: an accepted Newton correction printed by JFNK.
- SOB m15 (job 12065584, cmfgen_validation_SOB/m15_2p82d): it84,
  luminosity stable at 8.5e9 (lum ratio ~1.0-1.09 [OK]) but in a
  FULL->LAMBDA->LAMBDA cycle: each FULL step demands a saturated
  decrease (MAXCH=1e7 sentinel, level wants >=99.999% annihilation at
  H II lev1, depth ~61-70), applied x5; res jumps ~2.4x then the two
  Lambda steps contract it.  Post-Lambda minima decline ~10%/cycle
  (2.07e13 -> 1.39e13) -- slow convergence of a stiff inner-ejecta
  level, not stagnation, but the cycle may need 100+ iterations and
  the spectrum is still changing ~70-80% per full iter.
- jfnk_smoke (nohup PID 826926, 8 threads, since 14:58): diverged
  model (L=5.4e29 vs LSTAR=4.3e8, se_res=3.3e32, pops held for 17+
  iterations while the lagged field converges around garbage).  Cannot
  recover; burning 8 cores on the login node.  Recommend kill.

### 13.3 Open questions for the JFNK path

- Eval cost: each coupled residual/JVP eval is a full continuum
  sweep (~minutes).  A Newton solve needing ~30+ evals is hours per
  great-iteration; viable only if Newton converges in a few its or
  the eval gets a reduced-frequency surrogate.  Not yet measured on
  a successful solve.
- Whether GS_FG_IMMATURE ever clears on a model that is still far
  from converged populations (smoke: churn 1.1e-3 vs maturity want
  <1e-3 at it29-30 -- marginal).

## 14. Online literature audit -> optimal convergence setup (2026-10-31)

Re-read the numerical-method sections of the non-monotonic CMF
literature and the JFNK-SE paper, and quantified the stability
criterion against the live m15 model.  Extends (does not replace)
[[KNOP_CONVERGENCE_STRATEGY]].

### 14.1 Ceulemans et al. 2025 (JQSRT 343, 109470) -- quantified for m15

Their perturbation analysis (Sec. 2.3, Eq. 16-17) gives the per-step
amplification factor of the standard discretization:

    lambda = [dtau^2 e^-dtau + d'(1 - e^-dtau - dtau e^-dtau)]
             / [dtau^2 - d(dtau - 1 + e^-dtau)]

with d = c10*dnu_obs, d' = c00*dnu_obs the frequency-coupling
coefficients times the OBSERVER-FRAME frequency change between
connected grid points.  |lambda|<=1 in general, but oscillatory
modes (-1 < lambda < 0) appear when d,d' are large -- i.e. when the
per-cell comoving drift exceeds the local frequency spacing and the
upwind stencil extrapolates rather than interpolates.  Fix: connect
successive spatial points' frequency grids by minimizing |dnu_obs|
(their Eq. 18, "frequency matching"); cost = extra boundary
conditions where frequencies leave the matched range.

Their boundary-condition result (Sec. 3.2/3.3, Fig. 7) is directly
on point for our store design: for a strongly non-monotonic field
(sinusoidal v, amplitude ~ several times the frequency window),
temporarily STORING and reusing previously computed intensities at
re-encountered frequencies gives ~2% agreement with a Feautrier
reference; DISCARDING them and filling boundary values (black-body)
fails outright.  Our lagged I_NEXT store is the "correct" variant;
any path that substitutes a fill (capped write, unarmed slot,
bracket-search miss) is the failing variant and should stay rare.

Measured on m15_2p82d (SN_HYDRO_DATA): ALL 70 depths have
sigma = dlnV/dlnR - 1 < -1 (entire model compressing; V rises
inward from ~49 km/s to ~7074 km/s).  Per-cell |dV| reaches 1650
km/s at depth 65 (sigma=-125).  The run's minimum frequency spacing
is 5 km/s, so the drift ratio NNF = |GAM|*ds/dln(nu) reaches
~330 frequency bins at the inner edge.  With
DO_NONMONO_FREQ_MATCH=F (seq2's setting) every compressing upwind
read uses the adjacent bin (offset 1) -- deep inside the
extrapolation regime the eigenvalue analysis says is oscillatory.
seq2's refresh is currently healthy anyway (res 8.4e-3 in <=6
passes, clamp events=0) because at it1-2 the field is smooth; the
regime bites once line cores develop depth.

### 14.2 Knop et al. 2009 Paper II -- ordering and iteration

Sec. 3.2-3.3: the wavelength-coupled system is iterated by a
Gauss-Seidel march with a PHYSICALLY GUIDED index order -- follow
the pulse along the characteristic, march frequency in the local
upwind direction (sign of a_l).  With a linear source interpolant
the update is quasi-analytic (all RHS terms same-iteration): "the
solution is then direct and its speed is optimal".  With the
parabolic scheme only the spatial-downwind x wavelength-upwind
C-corner stays lagged.  This is exactly dev7's SWEEP_FLD_ORD /
KN_RAY_ORDERED machinery.

Sec. 4: the ALI step iterates (1-eps*Lambda*)J = J^FS -
eps*Lambda*J^old with Lambda* tridiagonal in wavelength, solved by
blockwise GS/Jacobi in wavelength with one-point overlap, to
internal rtol 1e-8; "ALL wavelength points have to be converged
simultaneously".  dev7's diagonal ALO is a documented cheaper
deviation.

### 14.3 Arramy, de la Cruz Rodriguez & Leenaarts 2024 (A&A 691, A44)
-- JFNK setup

Their JFNK for multilevel-atom SE (the residual side of our
coupled problem) gives concrete numbers:

- Krylov tolerance sweet spot: rtol ~ 3e-4 to 3e-2 (GMRES).  Tighter
  buys nothing; looser degrades the Newton direction.  Matches the
  eta_min=3e-4 forcing floor added at e764368.
- Convergence criterion: ||dn/nbar||_inf < 1e-4 AND ||F||_inf must
  have dropped >= 2 orders of magnitude -- the second clause exists
  to prevent premature exits on a plateaued residual.  Our driver
  checks the correction norm; the F-floor clause is worth adding if
  JFNK ever exits early on a stalled residual.
- FD JVP step: eps2 = eps_mach^(1/r) * sqrt(1+||n||^2)/||v||,
  r=2 forward -- same family as our H_STEP scaling.
- Preconditioner: left block-Jacobi (one local block per depth) --
  same object class as our BA matrix.  They warn it is
  insufficient for strong-scattering/many-frequency setups and
  suggest Scharmer-Carlsson single-point-quadrature or adaptive
  (Krylov-history) preconditioners as upgrades.
- Hybrid warm-up explicitly recommended: "a few RH92 iterations
  before switching to JFNK" to fix the poor-initial-guess failure
  mode -- our JFNK_WARMUP=6 + GS_FG_IMMATURE gate implements the
  same idea (Lambda = the RH92 analogue).
- Globalization: they list continuation/line-search as the missing
  robustness piece in plain JFNK.  Our solver already carries PTC
  + trust region + line search + Cauchy point -- stronger than the
  paper's setup; no gap here.

### 14.4 Optimal setup for the next Knop run (literature-conform)

Current seq2 deviations from the audited strategy, all fixable by
VADAT alone (no rebuild):

    KN_RAY_ORDERED=T          # seq2 has F -> legacy alternating
                              # passes; Knop II's ordered march is
                              # the paper-conform inner solve.
    DO_NONMONO_FREQ_MATCH=T   # seq2 unset (F).  With dV/cell up to
                              # 330 bins this remaps the compressing
                              # upwind slot to the characteristic
                              # foot (MATCH_FREQ, Ceulemans Eq. 18).
    DO_BA_RES_ACC=T           # seq2 has F.  Strategy doc Sec. 2.2.5:
                              # without resacc the stored-inverse BA
                              # path can emit unbounded corrections
                              # (the 1.1e12 event).  Keep the trust
                              # radius live on all KNOP runs.

Keep: USE_LAM_ES=T (FG ray path is the paper operator; MOM is a
moment analogue OBS cannot use), SOL_METH=DIAG (JFNK precond),
JFNK_WARMUP>=6, GS_LAM_CAP while immature, LAM_ENTER/EXIT=400/400.

For the SOB run: the FULL/LAMBDA cycle is a non-contractive BA map
at stiff inner depths (saturated 1e7% decrease on H II lev1,
applied x5, res x2.4 rebound).  The hybrid-iteration literature
supports holding the Lambda/operator-split phase until the map
contracts rather than re-entering FULL on a fixed schedule --
candidate settings: tighter BA correction cap at inner depths, or
raising the FULL-entry residual gate so Lambda runs until the
post-Lambda minima stop contracting (~10%/cycle currently).

### 14.5 What the literature says we do NOT need to add

- No new discretization: chi_hat = chi + |GAM|/dln nu (Knop Eq. 7)
  is already the positive-generalized-opacity form; xi=1 implicit
  wavelength coupling is what Paper II requires for strong
  coupling.
- No extra JFNK globalization: PTC + TR + line search already
  exceeds Arramy's setup.
- No direct sparse solve: ordered GS is quasi-analytic at linear
  source order (Knop II Sec. 3.2).

## 15. BAMATPNT cold-start crash + SUCCESS fix (02-Oct-2026, night)

### 15.1 Failure

Fresh JFNK_opt job 12066672 (KN_RAY_ORDERED/FREQ_MATCH/BA_RES_ACC=T)
died at end of iteration 1:

    Error opening SCRATCH/BAMATPNT in READBA   IOSTAT=29
    Major Error - cant read BA File
    Previously read successfully - before continuum loop

Call site: `cmfgen_sub.f` IF(.NOT. COMPUTE_BA) re-read after the
continuum loop (~line 4159 pre-edit).  `SCRATCH/BAMATPNT` was never
created: STORE_BA_DATA_V3 only runs when `COMPUTE_BA .AND. WRBAMAT
.AND. .NOT. LAMBDA_ITERATION`.

### 15.2 Root cause

Uninitialized local `SUCCESS` in `CMFGEN_SUB` startup (~line 1242):

    CHK=.FALSE.
    IF(.NOT. NEWMOD)THEN
      CALL READ_BA_DATA_V3(...,SUCCESS,...)
    END IF
    IF(.NOT. SUCCESS .OR. LAMBDA_ITERATION)THEN
      NLBEGIN=0
      COMPUTE_BA=.TRUE.

For `NEWMOD=T` the read is skipped and `SUCCESS` is stack garbage.
If it comes up .TRUE. the ELSE branch leaves `COMPUTE_BA` at its
module default (.FALSE.): the SE phase then skips the BA store, and
the post-continuum `.NOT. COMPUTE_BA` re-read hits the never-written
pointer file -> STOP.  Identical code exists upstream
(cmfgen_v7_ord_run); latent stock bug, nondeterministic.  seq2
survived on luck; JFNK_opt got unlucky.  Not related to
DO_BA_RES_ACC -- the flag only armed the resacc machinery.

### 15.3 Fix

`SUCCESS=.FALSE.` before the `.NOT. NEWMOD` read (commit 7c6a3c3).
A new model has no stored BA to reuse, so "read failed" is the
correct cold-start state and deterministically forces
COMPUTE_BA=.TRUE.  Restart path unchanged -- READ_BA_DATA_V3 still
sets SUCCESS there.  Verified: resubmitted job 12066716 completed
it1 with `BA matrix computed` and no pointer-file error.

### 15.4 Fresh-run status (12066716 / 12066717, SOB v2 12066674/75)

JFNK_opt it1 (1191.6 s):
- freq-match live: matched/compress = 13.26M/297.8M, max offset 36
- I_NEXT refresh residual max 2.43 (tol 2e-3), 40.3M/778M cells
  unconverged, 0 diverged rays, 0 clamps -- stiffer true operator
  than the old offset-1 stencil, as expected
- immature-field gate -> LAMBDA_TFIX for it2 (gs_freeze_forced)
- BA_RES_ACC armed: 31 saturated depth corrections zeroed
  (|STEQ|>1e4), nrej=0, trust=1e4
- se_res 2.35e14 / res2 7.06e11; lum ratio 1.73e-2 [FAR] -- early

SOB v2 (it33): Lambda phase contracts (maxapp 228->152%, res2
~1.3e11) but the FULL->LAMBDA rebound persists via the
`cnt_lam_escape` path (CNT_LAM>=RD_CNT_LAM forces a full step
regardless of CNT_STABLE_LAM).  it33 forced-FULL demanded 1e7%
(saturated), applied 430 under the trust cap, demoted back to
LAMBDA -- damped rebound, res2 stayed ~1.3e11 (old run rebounded
x2.4 uncapped).  Residual driver still H2:lev1 at depth 70.

Watch items: whether the Knop refresh residual contracts over the
Lambda window (43-pass cap was hit on it1's cold start); JFNK first
activation once field matures + warmup=6.

## 16. First live JFNK solve on the corrected Knop operator (03-Oct-2026, early AM)

Run: JFNK_opt 12066716, source 1dc1770.  JFNK fired at it11 (first
FULL designation after JFNK_WARMUP=6 with GS_FG_IMMATURE cleared;
refresh residual 2.8e-2, F(0) spread 2199/128730 eqs >50% --
localized store lag at depth 70, gate adopted eval baseline).

### 16.1 Measured failure signature

    Initial |F|_w = 2.6855E+13  ||F||_inf = 1.1795E+13
    Linearity probe: F(3v) -> 5.2e29, 4 Inf; rel err (2v,3v) = 1.0
    FD-JVP vs eps: ||Jv||inf = 3.6e32 / 8.9e31 / 1.7e32 / 3.6e32 /
                   8.0e32  (eps 1e-4 .. 1e-8 -- never converges)
    JVP audit: |FD|=5.2e29 vs |BA*v|=7.4e13, rel diff 7.1e15
    T-col audit: 70/70 depths >10% FD-vs-BA mismatch (JFNK_FIXT)
    PC: ||Mz-F||/||F| ~ 3e2; 52+ ZMAX hits per GMRES; forcing-norm
    clip aborts before Krylov can start (del-independent)

Newton outcome: 8 PC-rebuild attempts x 8 clip retries, all abort;
|F|_w drifted 2.9e13 -> 1.6e23 -> ... -> 8.6e23 across attempts
(last-eval residual at rejected probes, not committed state).
GMRES never got past its first matvec (best: 1 its, frac=0.80).
Cost ~3.2 h for it11 (incl. once-only JFNK_VERB diagnostic suite,
DIAG_DONE latched).

### 16.2 Recovery path (validated)

JFNK_DONE -> stagnation path (finite F, not catastrophic):
`NOT CONVERGED - Lambda fallback, stagnation 1 of 3`; safety
SOLVEBA under LAMBDA_TFIX/MAX_LIN_COR committed a capped step;
res 1.07e13 -> 8.47e12 (below it10's 8.8e12 -- net contraction
preserved).  `JFNK_MAX_STAGNATE=3` -> self-disables after two more
failed fires; diagnostics are once-per-process so re-fires cost
only F(0) eval + retry storm (~30-60 min each).

### 16.3 Interpretation

JFNK was admitted while |F|_w ~ 1e13 -- far outside the basin where
the FD-JVP is a usable directional derivative (eps sweep flat at
~1e32 because even eps*v probes leave the linear regime; the BA
preconditioner is the *untransformed* Jacobian while F is the
coupled SE+RTE residual, rel diff ~1e15 -- it cannot precondition
this operator).  This is the Arramy+2024 hybrid lesson: RH92/Lambda
warm-up must carry the model much deeper before Newton is worth
firing.  Warmup-count + field-maturity gates are necessary but not
sufficient -- the missing gate is residual magnitude at F(0).

Planned surgical fix (next submission, not this run): admission
gate `JFNK_F0_MAX` in JFNK_SOLVE_CORR_V1 after the F(0) eval --
if |F|_w > JFNK_F0_MAX defer Newton through the F0_OK=F path
(catastrophic-style: Lambda safety step, JFNK stays armed, no
stagnation charge, one ~5-min probe per full iteration).  Default
0 = off -> stock/flag-off path unchanged.  Threshold TBD from the
|F|_w values measured at subsequent re-fires; first estimate
~1e11 (two decades below the observed failure).

### 16.4 SOB v2 (12067329, restarted with NUM_LAM 2 -> 50)

NUM_LAM raised so the cnt_lam escape leash no longer preempts the
LAM_N_STABLE=4 merit gate.  Result: sustained Lambda phase,
MAXCH ground down 911% -> ~10% over ~30 its; stab counter reached
4/4 for the first time -> it192 merit-based `lam_stable_exit`
FULL step.  The FULL step still demanded 1e7% but the worst cell
*migrated* (H2:lev1@d70 -> OV:2p2_3Pe[1]@d37, sigma=-0.999);
resacc applied 404, res 2.41e12 -> 2.39e12, controller demoted
via maxch_above_enter.  Staged switching works: Lambda grinds the
stiffest cell, gate admits a full step on evidence, resacc bounds
it, demotion is earned.  res ~2.4e12 falling; res2 ~1.5e9.

## 17. it12 poison mechanism resolved + JFNK_COR_GATE validated (03-Oct-2026)

### 17.1 The it12 wreck was the safety step, not JFNK itself

Re-checking the cancelled run's it11 ITER_SUM: the JFNK catastrophic-
fallback SOLVEBA call ran LAMBDA_ITERATION=T with T1=MAX_LIN_COR and
applied the saturated demand UNFILTERED -- `pops 1.0E+07 / 1.0E+07`
(demand/applied).  Lambda-mode calls never arm BAACC_PENDING
(solveba_v13:382) and Lambda follow-ups cannot reject (line 272), so
the 1e7% application was un-rollbackable -> it12 measured res=2.5e23,
negative luminosity ratio, 15.7M bad field cells.

### 17.2 Gated run 12068898 (commit 162d39e) -- it11/12 outcome

`JFNK_COR_GATE=1.0D2` deferred Newton at it11 as designed:

    JFNK_COR_GATE: last demanded MAXCH=  1.2099E+08 >= gate
    ITER_MODE: it= 11 FULL_TFIX -> LAMBDA_TFIX
      reason=rel_stable_hold+gs_freeze_forced+jfnk_defer+maxch_above_enter

The deferred call took the SAME demanded correction (identical
DEC_VEC) but through the normal lam=F path: applied 9.0E+02 (resacc
trust/scaling) instead of 1.0E+07.  Result at it12:
res 8.47e12 -> 6.38e12 -- contraction preserved, no wreck, no storm,
pend=T armed normally.  The gate both skips the deterministic retry
storm AND routes the big demanded step through the protected path.

Mode controller note: lam_stable_exit uses MAXCH_CTRL=MAXCH_APP under
DO_NONMONO_KNOP, and LAM_ENTER=LAM_EXIT=400 means every capped Lambda
step (applied ~200 < 400) exits -> structural FULL/LAMBDA alternation.
Demand stays saturated (it12 demanded 1.49e7) -> gate defers every
full-designated call until the stiff HI:1@d70 level stops demanding.
This is the intended "Lambda/capped-full grind into the basin, then
Newton" sequence.

### 17.3 SOB v2 status (12067329, it~340)

Sustained merit cycle after NUM_LAM=50: five lam_stable_exits so far
(it192/212/237/254/279), each post-full rebound smaller; res floor
falling 2.4e12 -> 1.3e12, res2 ~4.7e8.  MAXCH grinding ~25 toward
the exit=20 stab gate; stiff cell now H2:lev1@d70 again.

### 17.4 SOB v2 limit cycle resolved: BA_ENTER_RAW demand gate (03-Oct-2026)

The period-7 sawtooth (lam_stable_exit -> saturated ~4e5% full step ->
maxch_above_enter demote -> ~6 its recovery) was traced to the merit
exit criterion itself: CNT_STABLE_LAM counts consecutive calls with
MAXCH_CTRL < LAM_EXIT, which measures the LAMBDA iterate's applied
corrections and says nothing about whether the BA/Newton direction is
contractive.  On this model the stiff inner-ejecta cells (H2:lev1@d70,
CV:1s_2p_1Po@d13) keep the raw BA demand saturated at ~1e5-1e7% no
matter how quiet LAMBDA gets.

Two resubmitted schedule experiments confirmed VADAT-only tuning
cannot break the cycle:

- LAM_EXIT 20 -> 5, NUM_LAM 50 -> 200, MAX_LIN 5 -> 1.5 (job 12071178):
  the it509 damped kick (+50% applied cap) still rebounded the SE
  residual ~120x (1.4e12 -> 1.7e14).  Damping the committed step does
  not prevent the rebound -- the rebound is set by the direction, not
  the step size.

Fix (source, minimal diff, cmfgen_jfnk working tree @ 312b04f +
uncommitted gate):

- New optional VADAT key [BA_ENTER_RAW] (default 1.0E30 = disabled).
  BA_LAST_RAW records raw MAXCH of every non-LAMBDA SOLVEBA call.
  Under DO_NONMONO_SOB a lam_stable_exit is held while
  BA_LAST_RAW >= BA_ENTER_RAW; the hold resets CNT_STABLE_LAM so
  CNT_LAM keeps counting and NUM_LAM escapes still fire -- each
  escape is a damped probe that re-measures the demand and re-opens
  the gate once the Newton direction contracts.
- ITER_REASON tag 'ba_demand_hold'; MODCTL gains 'lraw=' field.
- Files: control_variable_mod.f (decls), rd_control_variables.f
  (optional read), solve_for_pops.f (record + gate + MODCTL).

Resubmit 3: job 12071288 (N3), VADAT BA_ENTER_RAW=1.0D4, NUM_LAM=200,
LAM_EXIT=5.0, MAX_LIN=1.5.  Verified live: it528 shows
`LAMBDA_TFIX -> LAMBDA_TFIX reason=ba_demand_hold` where the old code
would have fired lam_stable_exit.  LAMBDA grind then continued at
~2% corrections (below the old 6.9% attractor).  First cnt_lam probe
expected near it~715.

## 18. Masked per-depth BA acceptance (07-Oct-2026)

Audit of why KNOP and SOB converge under LAMBDA but struggle on full
iterations: the resacc net was armed but effectively never fired (zero
BA_RES_ACC REJ events across all ten active run directories) because
the mode controller demotes full steps before the pending test can
judge them, and a LAMBDA follow-up cleared BAACC_PENDING untested AND
counted it as acceptance (trust x1.25, BA_NREJ_CONSEC reset -- a
ratchet that rewards unverified steps).  When the test did fire it was
binary and global: one stiff column's rebound rolled back a step that
improved the other ND-1 columns.  Both live runs solve SOL_METH=DIAG,
so full and LAMBDA steps share the block-diagonal solver; the
difference is the linearized field-response terms (COMPUTE_BA gated
VCHI/DIFFW assembly) and the caps, not the solve structure.

Section 17.4's measurement is the design constraint: damping alone
cannot work -- a +50% capped kick still rebounded the SE residual
~120x, because the rebound is set by the (locally wrong) linearization
direction at stiff compression cells, not by the step size.  The fix
must therefore localize the rejection, not shrink the step.

Implementation (per-depth masked acceptance, flag-gated):

- New optional VADAT key [BA_REJ_FRAC] (default 0.25; <=0 restores the
  legacy all-or-nothing test).
- SOLVEBA_V13 persists the per-depth residual history: BAACC_COL_PREV
  (armed snapshot of COL_RES, the per-depth column max already computed
  for BA_SE_RES/RES2) and BAACC_COL_MIN (running per-column minimum for
  the HARD bound).  RES2 ranking now decimates a copy, not COL_RES.
- Pending test on a non-LAMBDA follow-up: column I fails on
  res > GROW*prev(I), res > HARD*min(I), or res > 0.5*HUGE (the
  clamped-nonfinite floor).  BA_RES_SAT exempts columns that were
  already degenerate.
- Masked rollback (1..FRAC*ND columns failed, RES2 channel quiet):
  failing columns restore from POPS_BA_SAVE, their fresh STEQ is
  zeroed post-solve (it was linearized at rejected pops), the rest of
  the step stands and the call continues.  BA_PART_THIS flags the
  event: SOLVE_FOR_POPS forces COMPUTE_BA and tags ITER_REASON
  'ba_partial_rej' / ba=MASKED, but does NOT demote and does not count
  it toward the N_FULL_TREL phased T release.  Under KNOP the lagged
  I_NEXT store is still restored wholesale -- a lagged field one step
  stale at kept columns is normal; a store seeded at rejected pops is
  the it49-56 wreck path.
- Whole-state rollback (RES2 channel failure or >FRAC*ND columns):
  unchanged legacy path including the MAXCH sentinel and demote.
- Lambda follow-up honesty: the pending verdict is evaluated (no
  rollback -- the it33 cycle stands); a measured degradation halves
  BA_TRUST_CUR and counts BA_NREJ_SOFT instead of clearing the flag as
  an acceptance.

Files: control_variable_mod.f (decls + docs), rd_control_variables.f
(BA_REJ_FRAC read, state init), solveba_v13.f (mask test, partial
rollback, soft accounting, non-destructive RES2, arm-side col
storage), solve_for_pops.f (BA_PART_WAS consume, MASKED marker,
N_ACC_FULL exclusion).

Interaction with BA_DAMP_CAP (2cb4215): complementary layers -- the
cap bounds the size of each committed excursion while demand is
saturated; the mask bounds which parts of the step stand when the
direction is locally wrong.  BA_SAT_LIM still leashes a saturated
run to demote; BA_NREJ_LIM still hard-stops chronic rejection.
## 19. Non-monotonic mode prerequisites hard-enforced (03-Oct-2026)

The pipeline (prepare_VADAT.py) always stamped the right companion
keys, but a hand-edited VADAT could still select DO_NONMONO_SOB or
DO_NONMONO_KNOP while routing lines through an incompatible path --
the reader only warned and the run continued in a silently inconsistent
configuration.

rd_control_variables.f now hard-stops: DO_NONMONO_SOB requires
GLOBAL_LINE=SOB (or NONE + all present TRANS_*=SOB) and USE_LAM_ES=T;
DO_NONMONO_KNOP requires GLOBAL_LINE=BLANK (or NONE + all present
TRANS_*=BLANK), LAM_SET=F, and no FIXED_J / relativistic /
plane-parallel modes.  rd_cmf_flux_controls.f repeats the
line-treatment, LAM_SET and mutual-exclusion gates on CMF_FLUX_PARAM
so an OBS run cannot be pointed at an inconsistent line treatment.
DO_BA_RES_ACC=F and (SOB-mode) LAM_SET=T remain warnings.

Gate matrix verified against the built exe on a cloned smoke dir:
every invalid combo STOPs in the reader with an explicit LUER
message; the three valid configs (SOB+SOB+LAM_ES=T, KNOP+BLANK,
both flags off) pass the reader and run.  GLOBAL_LINE=NONE on this
VADAT aborts earlier in GET_KEY_STRING (TRANS_CrSIX absent) --
NONE needs the full per-species TRANS_* table.

## 20. Luminosity removed from convergence control (06-Oct-2026)

The models contain shocks, so radiative energy can be deposited or
extracted between the inner and outer boundaries.  Agreement or slow
drift of `RLUMST(1)/RLUMST(ND)` is therefore not a valid convergence
condition.  `BA_LUM_GATE` and `BA_LUM_DCHG` remain accepted no-op VADAT
keys for restart compatibility.  The ratio and call-to-call relative
drift remain in the structured `ITER_SUM` block as `lum diag [INFO]`,
but no longer demote BA, hold or cap Lambda corrections, delay mode
exit, veto MOM trials, or delay the rd_lam probe.  Runtime validation
is pending on both SOB and Knop/JFNK paths.

## 21. JFNK pseudo-time term re-anchored at the current iterate (06-Oct-2026)

Defect (code-read): `JFNK_RESIDUAL_V1` added the pseudo-transient
term as `x*s/dtau` with `x` measured from the solve baseline, so
`G(x) = F(x) + x*s/dtau` was a spring to `x = 0`, not to the current
Newton iterate.  `newton_hook_v1` measured `tol_`, ared/pred, the
line-search and Cauchy tests, EW forcing and the convergence test
on `|G|`, so "NONLINEAR CONVERGED" could be reported with
`F(x) = -x*s/dtau /= 0`, and at later iterates the term pulled
against progress.

Fix (SER form, Kelley & Keyes 1998, SIAM J. Numer. Anal. 35, 508;
Knoll & Keyes 2004, JCP 193, 357 Sec 4.2): the term is now
`(x - x_k)*s/dtau` anchored at the current iterate `x_k`
(`jfnk_ptc_x0`, re-anchored on each accepted step by
`jfnk_ptc_set_anchor_v1`; `JFNK_SET_PTC_SCALE_V1` sets the first
anchor to the baseline `x = 0`).  Every residual the Newton loop
reads is stripped back to `F` by `jfnk_ptc_strip_v1`, and the
predicted-residual JVP is stripped to the F-Jacobian action by
`jfnk_ptc_strip_dir_v1`, so all norms and the SER dtau update act
on `|F|`.  The GMRES operator is unchanged: the FD JVP still carries
the diagonal `s/dtau`, so it solves `(J + D/dtau) s = -F(x_k)` and
the `JFNK_PRECOND_V1` diagonal correction still applies.

Gating: all of it sits under `jfnk_ptc_active` (VADAT `[JFNK_PTC]`,
default T); with `JFNK_PTC=F` or `USE_JFNK=F` nothing changes.
Verification: builds clean; standalone Newton-core test in
`tests/jfnk_core/` exercises the anchored form on a stiff synthetic
system.  Runtime validation on a live model is pending.

## 22. Review of SOB and JFNK pathways -- recommended hardening steps (06-Oct-2026)

Code-read audit of both pathways against current practice (Knoll &
Keyes 2004; Kelley & Keyes 1998; Eisenstat & Walker 1996; Walker & Ni
2011 on Anderson acceleration; Arramy et al. 2024). Item 21 above was
the JFNK pseudo-time defect. Items 2 and 5 are implemented and await
live-model verification; the other follow-ups remain open.

1. SE residual acceptance (`SOLVEBA_V13`, `DO_BA_RES_ACC`) measures
   absolute rates, so the max channel saturates on the degenerate
   inner-boundary equations (~1E262) and needed the `BA_RES_SAT`
   exemption plus the RES2 rank channel.  Normalise each row by the
   sum of its positive rate terms (dimensionless relative residual)
   so the max channel is meaningful at every depth and `BA_RES_SAT`
   can go.
2. Implemented (06-Oct-2026): a degraded `BA_PROBE_LAM` trial is
   distinguished from an ordinary BA step through the resacc state.
   On its LAMBDA follow-up, the full pre-probe population and Knop
   lagged-field snapshots are restored and that call's stale Lambda
   solve is skipped. This avoids committing a probe known to degrade
   the residual without extending rollback to ordinary BA/Lambda
   cycles. Shared controller covers SOB and Knop/JFNK paths; runtime
   verification remains pending.
3. NG acceleration (stock) is applied unconditionally once scheduled
   and only kicks to LAMBDA afterwards if the change exceeded 50%.
   Judge the NG iterate with the same residual-acceptance test as a
   BA step and roll it back on growth; longer term, replace NG with
   depth-3 safeguarded Anderson acceleration (Walker & Ni 2011;
   Pollock & Rebholz 2021), of which NG is the unsafeguarded special
   case.
4. EW forcing safeguard in `newton_hook_v1` is
   `max(eta, 0.5*|F_k|/|F_k-1|)`; the Eisenstat-Walker Choice 2 form
   is `max(eta, gamma*eta_k-1^alpha)` when that exceeds 0.1.  Not a
   defect; swap only if runtime forcing-term behaviour looks odd.
5. Implemented (06-Oct-2026): the SOB OBS configuration now stamps
   DO_NONMONO_KNOP=T together with DO_NONMONO_SOB=T, GLOBAL_LINE=SOB,
   all TRANS_*=SOB, DO_SOB_LINES=T and SOB_FREQ_IN_OBS=T. The OBS reader
   permits this combination only in CMF_FLUX_PARAM: Knop formal J for
   the continuum/electron-scattering source, Sobolev line frequencies
   and SIGMA_TRUE EW estimation. The main VADAT modes remain mutually
   exclusive. Both SOB and frozen-J configs set COH_ES=F; the OBS
   preparation propagates the non-coherent scattering settings and
   ensures NUM_ES>=2. Runtime spectral verification remains pending.
6. Runtime validation of the Newton loop is still outstanding: no run
   in this project has produced a `JFNK RESULT` block. First target
   remains the SOB-seeded frozen-J stage
   (`configs/dessart_m15_2p82d_JFNK_frz_v1.yaml`) from a converged fresh
   SOB checkpoint on the latest committed executable.

Run record (06-Oct-2026, source commit `fa86fb7`): cancelled active SOB
job `12101061` and its dependent OBS `12101062`; old JFNK job `12095678`
was already completed and `12095679` failed. Fresh SOB submitted as
CMFGEN `12101701` (running) and OBS `12101702` (afterok dependency).
Generated VADAT confirms `COH_ES=F`, `GLOBAL_LINE=SOB`,
`DO_NONMONO_SOB=T`, and `DO_NONMONO_KNOP=F` in the population solve.
The JFNK stage is intentionally deferred: `run_model_grid.py` stages
warm-restart files at submission time, so submitting now would seed it
from a cold or unconverged state, not from this fresh SOB result.

## 23. OBS executable missing after a partial clean build (06-Oct-2026)

Failed OBS job `12095679` did not reach `cmf_flux`.  Its `batobs.log`
ended at `nice: .../cmfgen_jfnk/exe/cmf_flux.exe: No such file or
directory`; the later `mv: cannot stat` errors were cleanup fallout.
The parent CMFGEN job completed at 12:54:54, while the OBS binary was
rebuilt manually at 12:57:19.  Root cause: `build.sh` ran `make clean`
(deleting `exe/*.exe`) but only rebuilt the `cmfgen` target, leaving a
window in which queued OBS jobs referenced a missing executable.

Fix: `build.sh` now runs `make KEEP_EXE=1 clean` followed by
`make DEBUG=nodebug all`, so the previous executables remain installed
while objects/libraries rebuild and the full production/utility tree is
rebuilt from one command. The main and OBS link rules publish through
temporary names plus `mv`, so a queued job cannot see a truncated
binary. The full build exposed one legacy standards error in
`misc/plt_steq.f` (`ALLOCATE(XV(1.2*NT), ...)`); the bound is now
`NINT(1.2*NT)`.  Submission also preflights `${CMFDIST}/exe/cmf_flux.exe`
before queueing OBS unless `SKIP_OBS=1`.

Monitoring: `submit_cmfgen_jobs.sh` now writes and submits
`run_cmfgen_CONVERGENCE.slurm` for every FINAL job it parses.  The job
uses `--dependency=after:<FINAL>` so it starts when the solver starts,
runs `check_convergence.py --plot CONVERGENCE_HISTORY.PNG` every 30 s,
and writes a final plot after the solver leaves the queue.  On Viper it
inherits the model partition; APU monitor submissions require the
smallest permitted GPU request (`--account=mpe_apu --qos=a0004
--gres=gpu:1`).  The generated FINAL/SINGLE_FINAL script also keeps a
best-effort in-job 30 s plot loop.  Live verification: monitor job
`12102126` is running against SOB `12101701` and refreshed
`CONVERGENCE_HISTORY.PNG` at 30 s cadence.

## 24. FIX_IMPURITY masked population corrections (06-Oct-2026)

The SOB validation run generated `T [FIX_IMP]` because the SOB profile
stamped the key explicitly.  `FIXPOP_IN_BA_V3` then removed BA equations
for trace species or ionization stages below its abundance/fraction
thresholds.  In the live run this appeared as persistent high residuals
and CORRECTION_SUM bands dominated by zeros during full iterations: the
full solve was blocked from correcting a large set of population
variables at the affected depths.  Lambda iterations could still relax
the residual, so the failure is diagnosed as suppressed full-iteration
corrections rather than a missing Lambda path.

Production change: `FIX_IMP` is `F` in every current SOB, Knop and JFNK
YAML profile, and `NON_MONO_SOB_RECOMMENDED` now stamps `F`.  The key
remains readable for compatibility and FIXPOP_IN_BA_V3 remains available
for explicit experiments, but the submitted recipes no longer freeze
trace-population corrections.

Run record (06-Oct-2026, source commit `c13f017`): cancelled stale SOB
jobs `12101701`, `12101702` and monitor `12102126`, then resubmitted the
same `dessart_m15_2p82d_SOB_mask` recipe.  New jobs: CMFGEN `12103504`,
OBS `12103505` (`afterok`), convergence monitor `12103506` (`after`,
30 s).  The generated VADAT contains `F [FIX_IMP]`; iteration 1 is a
LAMBDA_TFIX call and CONVERGENCE_HISTORY.PNG is being produced.

## 25. Live SOB diagnosis and correction caps (06-Oct-2026)

Inspection of job `12103504` (`m15_2p82d_mask`, OUTGEN/ITER_HISTORY):
the run is not diverging.  The Lambda sequence reduced the applied
correction demand from `1.65e11` at it1 to `99` at it20, probed a full
step at it21, measured a worse SE residual at it22, rolled the probe
back, and returned to Lambda recovery.  By it38 the maximum SE residual
had fallen from `6.33e13` before the rejected probe to `4.06e13`; the
population demand drifted back up to `55%`, so progress is slower but
still bounded.  `LUM_RATIO_OI` remains diagnostic-only and did not gate
the mode.

Caps: `MAX_LAM=1e10` was already effectively unclamped and is not the
limiter.  `MAX_LIN=5` is not the binding cap on the rejected probe:
`BA_PROBE_LAM`/saturated-demand steps take `min(MAX_LIN,BA_DAMP_CAP)`,
so `BA_DAMP_CAP=1.3` and residual acceptance dominated the step.  To
make the ordinary full-step cap non-limiting while retaining the
safeguard, production profiles now use `MAX_LIN=1e1` (the CMFGEN
manual's typical value).  `MAX_LAM=1e10` is explicit in the SOB and
staged frozen-JFNK recipes; all current SOB/KNOP/JFNK YAML profiles now
use `MAX_LIN=1e1`.  `FIX_IMP=F` is unchanged.

The running job was left in place: increasing `MAX_LIN` cannot alter a
LAMBDA iteration and would not have changed the `BA_DAMP_CAP`-bounded
probe that was rejected.  If the segment requeues, its generated script
rebuilds VADAT from the updated config and will use the larger
`MAX_LIN` without disturbing the current state.

## 26. Masked partial rollback extended to rd_lam probes (13-Oct-2026)

Audit of the DO_NONMONO_SOB BA/full-iteration math found the Jacobian
internally consistent (fixed-J operator, same as stock USE_LAM_ES:
do_var_cont.f zeroes TX/TVX for the non-mono branch, which is what the
empty stock USE_LAM_ES branch leaves behind anyway; Sobolev VB/VC are
true derivatives of the sigma-floored escape-probability rate; sigma is
not a population variable so the floor's kink never enters dF/dn;
INC_ADV=F so no missing advection derivative; probe calls run
newmat=T with a freshly assembled full matrix).  The observed probe
failures are a stiffness/step-size problem at one boundary cell, not a
Jacobian defect: each probe degraded the H2:lev1 ionization equation at
depth 70 ~20x while improving the other ~69 columns ~45%.

The defect was operational: solveba_v13.f forced a FULL pre-probe
rollback for any failing column on an rdlam probe
(BAACC_RDLAM_PROBE .AND. WLD_REJ -> REJ_STATE), bypassing the masked
REJ_PART path designed for exactly this case.  Each probe therefore
discarded the bulk improvement and cost two iterations.

Change: the probe override is now gated on REJ_STATE (bulk failure:
RES2 channel degraded, or more than BA_REJ_FRAC*ND columns failed) so
a few-column probe failure takes the masked path -- failing columns
are restored to pre-probe pops, kept columns keep the step, the call
proceeds to its Lambda solve with masked corrections zeroed, and the
probe still counts as rejected (BA_PART_WAS -> BA_RAW_REJECTED) for
the BA_ENTER_RAW/JFNK gates and the trust-radius halving.  Boundary
columns need no special handling: per-depth SE is local (NUM_BNDS=1),
the restored column is simply held at pre-probe state for one call.
Verified against solveba_v13.f snapshot ordering (POPS_BA_SAVE holds
n0 on the probe call and is not refreshed on the Lambda judging call).

A second gap closed in the same commit: the masked-correction zeroing
inside SOLVEBA_V13 sits in a `.NOT. LAMBDA_IT` block, so on the LAMBDA
judging call a partially rejected probe would restore its failing
columns to pre-probe pops but still apply the Lambda correction
computed from the rejected-state residual -- mixing states at exactly
the columns the rollback protected.  A parallel zeroing pass is now
gated on `DO_BA_RES_ACC .AND. LAMBDA_IT .AND. BA_PART_THIS`
(BA_PART_THIS is only set when pops were actually restored this call;
the NFAIL mask is also populated on soft-rejected armed steps where no
rollback happened and holding those columns would discard a valid
correction).  The stale masked demand is not recorded into
STEQ_MASKED_MAX because it was measured at the rejected state.

## 27. SOB resubmission on masked-probe build (13-Oct-2026)

Cancelled the FIX_IMP-free SOB triplet (12103504/12103505/12103506)
and resubmitted configs/dessart_m15_2p82d_SOB_mask.yaml through
run_model_grid.py on commit 3d4e4d3 (masked partial rollback for
rdlam probes + LAMBDA-call mask zeroing, exe/cmfgen_dev.exe rebuilt
19:49).  New jobs: FINAL 12105832, OBS 12105833, CONVERGENCE
12105834 (after:FINAL, 30 s cadence).

This is a COLD START: the submission pipeline removes the run
directory and re-copies the base template, preserving only
SN_HYDRO_DATA / OUTGEN / SPECIES_MASSES / OBSFLUX / RS_XRAY_FLUXES;
POINT1/POINT2/SCRTEMP/SCRATCH are not preserved, so the ~150 Lambda
iterations of the previous run are discarded.  POINT*/SCRTEMP will
not exist until the first checkpoint write.  Note for future
resubmissions: warm restart requires copying SCRATCH state back in
before submission.

VADAT_INIT notes: REL_OBS=T and REL_CMF=T are now set in the init
template (the earlier F/F policy mismatch for SN models is resolved).
Resolved config carries FIX_IMP=F, MAX_LIN=1.0D1, MAX_LAM=1.0D10.

## 28. REV_RGRID made usable under non-monotonic flags (13-Oct-2026)

Grid-resolution audit for the m15 SOB run found three defects that
made REV_RGRID unusable, now patched (all in this commit):

1. rd_control_variables.f: REV_RGRID was bundled with DO_HYDRO in the
   non-mono STOP.  DO_HYDRO stays blocked (it rewrites the prescribed
   V(R)), but REV_RGRID is now permitted -- the method-level guard in
   ADJUST_R_GRID_V4 already rejects VEL_SPACE (multi-valued R(V) map),
   and the remaining methods (TAU_SPACE/UNIFORM/MODUN/VTAU/FIX_NX/
   REFINE/FULL_R/OB_ONLY) interpolate LOG(V) on the monotone R
   coordinate, which is single-valued and safe.
2. adjust_r_grid_v4.f: the interpolated-sigma check clamped sigma<-1
   only for V<0.2 km/s and STOPped otherwise -- guaranteed to abort on
   a shocked field where sigma<-1 sits at thousands of km/s.  The
   check is skipped under DO_NONMONO_SOB/DO_NONMONO_KNOP; the
   interpolated SIGMA is transitional since RECOMPUTE_SIGMA_TRUE
   reclamps it under SOB immediately on return, and KNOP handles
   sigma<-1 natively.
3. solve_for_pops.f: BA_MATRIX_STORED is cleared on R_GRID_REVISED.
   A stored BAMAT's depth-column indices no longer map to the same
   radii; FIX_BA/N_FIX_BA recall would apply corrections to the wrong
   cells.  Not flag-gated (upstream gap on every regrid path).  The
   post-regrid call is forced LAMBDA (cmfgen_sub.f), which does not
   store, so the flag only ever protected WRBAMAT recalls.

REV_RGRID remains F by default; enabling it under non-mono now uses
ADJUST_R_DEFAULTS (already generated: TAU_SPACE/UNIFORM, N_ITS=5,
STRT_ITS=15, FREQ_ITS=10).  Candidate methods for this model: VTAU
(biases points to low-velocity inner ejecta where the stiff depth-70
cell sits) or FIX_NX with an FG_MIN..FG_MAX tau window covering the
inner boundary (log tau ~1.6).  Job 12105832 is still PENDING and
will pick up this binary at start; its config leaves REV_RGRID=F, so
the only live change for it remains the masked-probe path.

## 29. RDINR-side boundary refinement (Luc's N_IB_INS/N_OB_INS) (13-Oct-2026)

Luc's SN setups appended extra boundary points via N_IB_INS/N_OB_INS
in SET_SN_R_GRID (geometric inner sub-cells at ~0.4/0.16/0.06 of the
last spacing, equal-log outer steps).  Those keys are dead under
RD_IN_R_GRID=T -- SET_RV_HYDRO_MODEL_V3 reads RDINR verbatim and never
reaches ADJUST_SN_R_GRID -- so the VEL_LAW=11 hydro-grid models have
been running without the boundary refinement Luc considered standard.

Restored RDINR-side in the pipeline (untracked scripts/, not this
tree): prepare_VADAT.py refine_rdinr_boundaries() splits the
innermost/outermost RDINR cell into geometric log-R sub-cells (finest
at the boundary; ratio 2.5 gives 0.039/0.137/0.384 fractions, close
to Luc's).  Applied after target_nd downsampling and
rdinr_max_dv_kms refinement; nd_written feeds MODEL_SPEC ND, and
SET_RV_HYDRO_MODEL_V3 re-derives V/SIGMA from the hydro spline at each
inserted radius, so R/V/SIGMA stay physically consistent.  Keys:
cmfgen.rdinr_ib_ins / rdinr_ob_ins / rdinr_bnd_ratio (CLI
--rdinr-ib-ins etc. through run_model_grid.py and submit_cmfgen.sh;
prepare_VADAT also falls back to the YAML so a pending job picks up
edits on regeneration).  Defaults 0/0 = off, matching the
grid-adjustments-off policy.

Enabled rdinr_ib_ins=3, rdinr_ob_ins=3 in
configs/dessart_m15_2p82d_SOB_mask.yaml and patched
config_resolved.yaml so the PENDING job 12105832 regenerates
RDINR/MODEL_SPEC at ND=76 on start (batch.log should show the
"[RDINR boundary] inserted" line; OUTGEN echoes ND).  Verified on the
real RDINR: 150 -> 156 pts, strictly monotone R, sub-cells land
correctly at both edges.  Caveat: the run now tests the masked-probe
fix and the boundary grid together -- attribution between the two
rests on OUTGEN probe diagnostics, not a clean A/B.

Why it may help: the stiff depth-62..70 residuals sit where the inner
tau rise is steepest; the extra inner points sharpen the diffusion
boundary where the fixed-J Newton step currently overshoots.  The
outer points resolve the thin free-streaming edge (THK_CONT=F).
Still secondary to the missing ion stages for the H2:lev 1 row.

## 30. Per-equation quarantine for the masked BA acceptance (07-Oct-2026)

Final state of the SOB mask run (job 12105832, ended it509 COMPLETED
on schedule -- not converged): locked in LAMBDA_TFIX<->FULL_TFIX with
36 consecutive rejections, best se_res 4.58E12 flat for ~300 its,
lum ratio ~1.76.  Every armed step demanded raw 1E6-1E7%
(annihilation sentinel) on the same few equations -- CV:1s_2s_1Se and
1s_2p_1Po at d68-69, the NeV ground term at d70, H2:lev1 at d75 --
while the trust radius sat at its 1.0 floor and zeroed ~60/76 columns
per solve.  The per-column granularity was the defect: one
pathological equation vetoed every other variable at its depth,
including T, so no full step ever committed.

Two changes in this commit, both under DO_BA_RES_ACC (default path
unchanged):

1. Per-equation quarantine, VADAT DO_BA_QUAR (default F) +
   BA_QUAR_MAX (default 25; <=0 disables growth).  An (equation,
   depth) pair is parked on measured failure only: every over-trust
   population variable in a trust-masked column, and the worst
   residual equation of each rollback-masked column
   (BA_QUAR_WORST_EQ scans SE%STEQ with the same EQ_SPECIES/EQXzV
   mapping as the residual metric).  A parked equation is excluded
   from COL_RES/BA_SE_RES/BA_SE_RES2 (entry measurement rescans the
   ion's rows skipping parked ones only when the MAXLOC winner is
   parked), and its STEQ correction is zeroed on non-LAMBDA solves
   before the trust test -- it can neither veto acceptance nor trip
   the whole column.  ED (NT-1) and T (NT) are never parked.  The
   mask is allocated once (persistent per run) and sized (NT,ND).
   LAMBDA calls skip suppression entirely, so parked levels still
   update under fixed-point iteration and their Lambda demand still
   enters MAXCH -- quarantine cannot manufacture a converged exit on
   a Lambda call.  Parked demand stays visible: per-call OUTGEN line
   "BA_RES_ACC: quarantined eqns (total/new)" plus the BA13-TRACE
   "quarantined eqns / |STEQ|" row, and each addition is logged at
   CONV_DIAG_LVL>=2 with its EQ_DESC label.

2. Adjudication accounting fixes (same diff, separate defect):
   GS_POPS_HELD calls defer judgement (the armed residual re-measured
   a frozen state and burned BA_NREJ_LIM on calls that could never
   pass; jfnk_cpl it71-75), and a masked/soft step that improved
   either residual channel resets BA_NREJ_CONSEC instead of counting
   (m15 2.82d died at it75 with ~15 of 20 counted events being
   residual-improving or population-frozen calls).

Design constraints kept: raw-demand saturation still damps/holds/
rejects; the applied cap never legitimises a bad direction; RES2
remains the bulk veto; quarantine is bounded, opt-in, and every
parked equation is named in the log rather than silently dropped.

Test config: configs/dessart_m15_2p82d_SOB_quar.yaml -- same settings
as SOB_mask plus DO_BA_QUAR=T, warm_restart_population_only from the
it509 checkpoint of m15_2p82d_mask, so the test lands directly in the
stiff regime instead of re-paying the cold-start Lambda window.
Success criteria: quarantined set converges to the ~5-15 offender
equations, trust-radius unsticks from the floor, an armed step
commits corrections at the previously vetoed depths, and se_res2
leaves the 4.6E12 plateau.  Failure mode to watch: quarantine growth
parking collateral variables on an explosive solve -- check the
CONV_DIAG_LVL=2 list against the known offenders.

Pipeline fix needed to submit (scripts/, untracked): the
warm_restart_from ND fingerprint compared source MODEL ND (76, the
solver grid after rdinr_ib_ins+rdinr_ob_ins boundary inserts) against
target SN_HYDRO_DATA ND (70, raw hydro grid) and refused the
continuation.  The check now adds the configured insert counts to the
hydro count before comparing, and the error names the solver ND.
Submitted as jobs FINAL=12110995, OBS=12110996, monitor=12110997 on
commit 11c446a (source in 80d28ea).

## 31. Shared DO_JFNK switch for KNOP and SOB (07-Oct-2026)

`DO_JFNK` is now the VADAT key that arms the JFNK population solve in
`SOLVE_FOR_POPS`; `USE_JFNK` is a deprecated alias (warns, still
works).  The read-time gate accepts either `DO_NONMONO_KNOP` or
`DO_NONMONO_SOB` (commit 970b4f6).

Under SOB only the frozen-J tier can be admitted: `JFNK_COUPLED=T` is
rejected at VADAT read because `JFNK_RAD_SWEEP_V1` replays the KNOP
continuum/blanketing chain and does not recompute `SUB_SOB_LINE_V3`
line rates - a coupled SOB eval would solve a different residual than
the outer iteration assembles.  Frozen-J SOB JFNK is the same fixed-J
system `SOLVEBA_V13` solves, run through the globalized Newton
machinery (log unknowns, trust region, line search, LAMBDA/BA
fallback).  The F(0) identity test still guards baseline
consistency: a mismatch aborts to the LAMBDA/BA path rather than
committing a wrong-operator step.

The admission chain needed no code change for SOB: the field-maturity
clause is already scoped to KNOP
(`.NOT. DO_NONMONO_KNOP .OR. .NOT. GS_FREEZE`), and `USE_LAM_ES=T` is
already mandatory under SOB.

Pipeline: `DO_JFNK` added to the copied-key list and JFNK detection
in `scripts/prepare_VADAT.py` and `scripts/run_model_grid.py`
(untracked edits in the pipeline repo).  Runtime validation on an SOB
model has not been performed.

## 32. VEL_LAW 13/14 grids + JFNK consistency pass on m15 2.82d SOB (08-Oct-2026)

Two workstreams on the stalled `m15_2p82d_SOB` convergence, source
commits `c958555`..`6296181` on `devin/m15-zlaw14-jfnk` (merged to
main, see git log).

### VEL_LAW 13/14 (`rv_sn_model_snhyd.f`, c958555)

Both laws read `SN_HYDRO_DATA` directly - no auxiliary files.
`VEL_LAW=13` applies the law-12 tau-based regridding recipe.
`VEL_LAW=14` builds a blended `dlog(tau_Ross) + 1.5*dlog(R)` grid
using the Rosseland kappa column: the outer log-R zone resolves the
optically thin layers where corrections peak (depths 1-9 demand >10%
on ~700 of 1839 variables every iteration), the log-tau interior
keeps thick-region resolution, and the blended transition zone closes
the 2-dex tau jump the two-zone grid had between the thin/tau zones
(tau ~ 9e-4 to 0.14 in one cell). The grid is construction-only - no
change to the physical input profile. Live check: `Computed blended
tau-R/log-R grid in RV_SN_MODEL_SNHYD` in OUTGEN.

Companion fixes in the same submission window:

- `do_tau_regrid_v2.f` TAU_OLD(1) NaN guard (be80914): FLUX_MEAN(1)
  can be zero at the outer depth -> LOG10(0) poisoned the new tau
  grid; now uses floored TA like the VTAU branch.
- `cmfgen_sub.f` DO_LAM_AUTO revert (3e165c6): the non-mono
  stable-count guard parked runs in LAMBDA at ~45% demand
  (LAM_EXIT=20 unreachable above the LAMBDA attractor floor); stock
  MAXCH<50 handoff restored.

### JFNK consistency pass (e56471f)

Audit found the solver mathematically correct (residual/JVP/PC chain
verified; see [[JFNK solver]]) but operationally inconsistent:

- `JFNK_MAX_EVALS=40` could not cover one GMRES cycle (~60 JVPs at
  m=10, 20 its, 2 restarts) plus per-Newton residuals and line-search
  retries - deterministic mid-solve aborts landed in the stagnation
  budget. SOB zlaw14 config now runs 200.
- `CMAT_RAD` was rebuilt every activation (two BA assemblies). Now
  cached up to `JFNK_CMAT_MAX_AGE=4` activations - safe because the
  fixed point is invariant to BA_rad (it multiplies x, which
  vanishes at the root). Preconditioner still rebuilt from live BA.
- `JFNK_DISABLED` persisted for the process lifetime but silently
  re-armed on requeue - restart-dependent retry. Now a deterministic
  cooldown: re-arm every `JFNK_REARM_EVERY=10` admitted calls, logged
  as `JFNK re-armed ... retrying Newton`.

### ITER_TYPE dispatch banner (41593fa)

Once per great iteration `SOLVE_FOR_POPS` prints a `~`-bordered
`ITER_TYPE it=N Running {LAMBDA,FULL,JFNK-FRZ,JFNK-CPL} ITERATION`
record so the mode mix is greppable without parsing mode-controller
transitions. Records the dispatch decision; an in-solve LAMBDA
fallback is still named by the fallback blocks.

### Test suites (`tests/`, 6296181)

- `tests/f90`: production `jfnk_solver_v1.f90` compiled unchanged
  against a stub int_mod - 21 checks (GMRES convergence/PC/NaN,
  Newton on linear+nonlinear, trust-region cap, EW forcing, PTC,
  banded strip helpers, FD-JVP). Second suite links the real
  interface objects - 9 checks (BA_RAD matvec incl. the always-exp
  T mapping, fixed-unknown mask, weighted dot). All pass.
  Complements `tests/jfnk_core` (Newton core only).
- `tests/python/jfnk_outgen_audit.py`: parses a live OUTGEN for
  ADMIT/DEFER/F(0)/RESULT/fallback blocks and the ITER_TYPE banners,
  asserts invariants (F(0) rel<1e-6, reuse-only-after-build,
  ADMIT/ITER_TYPE pairing), writes plots under `tests/plots/`
  including ITER_HISTORY panels (MAXCH vs the 50% handoff, SE
  residual, trust radius, lum ratio).
- `tests/python/plot_solver_tests.py` renders the per-iterate
  |F|/dtau history the solver test records via the sub() callback;
  the PTC run shows the SER dtau ramp 1 -> 1e8 (cap) - continuation
  relaxing into full Newton exactly as designed.
- `build.sh` runs tests/f90 after linking; a failed check fails the
  build like a compile error (plot step non-fatal).

### Live status

SOB frozen-J experiment `m15_2p82d_SOB_jfnk_zlaw14_jfnk` submitted
through `run_model_grid.py` (FINAL=12131626, OBS=12131627,
monitor=12131628): `DO_JFNK=T`, `JFNK_WARMUP=6`, both gates 0
(ungated admission on the first full-designated call),
`JFNK_MAX_EVALS=200`, `JFNK_FIXT=T`, VEL_LAW=14. Runtime JFNK
behavior under SOB is UNVERIFIED until the first ADMIT - the mode
controller still designates LAMBDA iterations while MAXCH is high
(`JFNK DEFER -- LAMBDA-designated iteration` through it15).

## 33. Null-row identity solve trips the precond clip (2026-10-08)

First live JFNK admissions of `zlaw14_jfnk` (it31, it33, it35) all
aborted identically: `GMRES: precond clip -- aborting` on the first
matvec, `info=3`, zero extra residual evals. F(0) consistency was
3e-17 (residual evaluator correct); the failure was in the
preconditioner apply.

Root cause: ~6300 rows per factorization regularized as near-null
(BA diag spanning 2.6e-53..5.5e6) and replaced by the identity row.
An identity row solves `z_i = v_i`, and entry residual components
reached ~2e12, so `|z|>1e8` fired on every apply. A clipped `M^{-1}`
is nonlinear, so the abort was mathematically correct - but the
identity-row treatment made it a certainty on any rough state.
Third consecutive abort disabled JFNK for the rest of that process.

Fix (08-Oct): `JFNK_ROW_REG(J,K)` mask persisted from factor time;
masked components are projected to zero at apply in every branch
(LU, LU-fail Jacobi, DGETRS-fail, pure Jacobi, banded DGBTRS, and
the log/PTC post-loop). Linear projection keeps `M^{-1}` consistent;
null directions get no update and stay with LAMBDA/globalization.
Detection now reads the raw BA_A row so ridge retries cannot un-mark
a null row. Counter: `jfnk_n_reg_proj` (`Reg-row projections`).

Also fixed a latent .mod-shadowing hazard: stray
`jfnk_cmfgen_int_mod.mod` files in `tests/f90/` (cwd resolves ahead
of `-I`) made the interface test compile against the stub and fail
the post-build test step. `build_and_run.sh` now removes cwd .mod
files before compiling.

Status: fix built and tested; awaiting the requeued zlaw14 job
(process restart re-arms the JFNK disable) for the first admission
on the new binary.
