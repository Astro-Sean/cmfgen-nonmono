# ES_ALI / ES_ALI_BA / BA_RES_ACC AUDIT

Audit date: 11-Nov-2026. Source revision audited: working tree on
`main` post `d6b48b24`. Scope: the three VADAT acceleration keys
`DO_ES_ALI`, `DO_ES_ALI_BA`, `DO_BA_RES_ACC` (plus the dependent
`DO_FG13_ALO`, `DO_FG13_ALO_TRI`, `DO_BA_RES_REUSE`) on the
production recipe `DO_NONMONO_SOB=T`, `USE_LAM_ES=T`,
`GLOBAL_LINE=SOB`, `COH_ES=T`, `DO_JFNK=F`, `DO_NG=T`, `DO_AV=T`.

Runtime evidence: `m15_2p82d_FRESH_sob` (converged it722,
maxapp 4.9e-2), `m15_2p82d_FRESH_sob_dieoff` (in progress, it212+,
tail contraction ~0.9973/it), `m15_2p82d_SOB_l14_lte_adp2_tauto`
(warm restart, `DO_ES_ALI=T`, converged it1445). All three VADATs
carry `DO_BA_RES_ACC=F` as an explicit override; the SOB absent-key
auto-set would have made it T.

This audit is read-only. No solver source was modified. The E1-E4
experiment matrix in section 7 is a proposal; nothing here authorises
activation in production.

## 1. Mechanism map

| Key | Consumer | What it changes |
|---|---|---|
| `DO_ES_ALI` | `comp_j_blank.f:829-1235` | J-loop update `RJ=TC` becomes an OAB diagonal step `RJNEW = RJ + (TC-RJ)/max(1-DJDETA*chi_scat, 0.01)`; loop exits on `ES_DJ <= ES_ALI_TOL` (1e-4) or `ES_ALI_MAX` (20) |
| `DO_ES_ALI_BA` | `do_var_cont.f:428-438` | Continuum BA Jacobian diagonal `dJ/dCHI`, `dJ/dETA` amplified by `DENB = 1/max(1-TB*chi_scat, ES_ALI_BA_MIN)` (default floor 0.1) |
| `DO_FG13_ALO_TRI` | `do_var_cont.f:451-463` | Adds stored K+-1 depth response bands to TX when `NUM_BNDS>=3`; inert unless V13 solver ran |
| `DO_BA_RES_ACC` | `solveba_v13.f:244-639+`, `solve_for_pops.f` (~150 sites) | Residual-judged acceptance of each committed population step: masked per-column or whole-state rollback to `POPS_BA_SAVE`, trust-radius halving, forced BA rebuild, LAMBDA demotion |
| `DO_BA_RES_REUSE` | `solve_for_pops.f:3242-3295` | Extends BAMAT recall past `VAL_FIX_BA` while `BA_JUDGE_RATIO` (RES2 channel) contracts <= `BA_REUSE_TOL` (0.9); `BA_REUSE_NSTALL` (2) consecutive stalls force rebuild; auto-set T under SOB when ACC is T |

All flags default F. `DO_FG13_ALO` auto-sets T under
`DO_NONMONO_SOB + USE_LAM_ES=T + .NOT. DO_SOB_FG_KNOP`
(`rd_control_variables.f:3757-3765`) and only controls BA-side
consumption; the V13 ALO store itself is written unconditionally on
every `FG_J_CMF_V13` call.

## 2. Pillar 1 - DO_ES_ALI

### 2.1 Update equation

Inside `COMP_J_BLANK`, `CONT_VEL` non-ACCURATE continuum branch, per
frequency, inside `DO WHILE(INACCURATE)`:

```
TA = ETA_CLUMP + CHI_SCAT_CLUMP * RJ          (coherent ES source)
TC = formal_solver(TA)                        (FG_J_CMF_V13 on this recipe)
DJDETA = stored dJ/dETA diagonal of that solver, same freq/depth
DEN    = max(1 - DJDETA*CHI_SCAT_CLUMP, 0.01)
RJNEW  = RJ + (TC - RJ)/DEN
if RJNEW <= 0: RJNEW = TC                     (positivity fallback)
RJ     = RJNEW
```

This is the Olson-Auer-Buchler diagonal approximate-Lambda step on
the coherent-ES fixed point `F(J) = J - Lambda[eta + chi_scat*J]`.
`DJDETA*chi_scat` is the diagonal of the ES self-coupling map, so the
division by `DEN` is the Jacobi preconditioner for the scattering
iteration. Algebra checks out: at the fixed point `TC -> RJ`,
`RJNEW -> RJ`, `ES_DJ -> 0`.

### 2.2 Timing and state

- Applied per frequency, per pass, inside the great-iteration
  frequency sweep, **before** `SOLVE_FOR_POPS` runs. Ordering is
  correct: populations always see the ES-relaxed J of the current
  sweep.
- `RJ` is the live ES iterate (seeded from EDDFACTOR, rebuilt into
  `TA` each pass). `ETA_CLUMP`/`CHI_SCAT_CLUMP` are current-iteration
  opacity/emissivity. `DJDETA` comes from the ALO store the solver
  just wrote this pass - same frequency, same iteration, never
  lagged across iterations.
- Inter-frequency coupling stays lagged by construction: pass *k*
  re-solves only frequency *k*; neighbours' fields are frozen. Same
  as the stock `RJ=TC` iteration - ES_ALI changes the contraction
  rate of the per-frequency ES subproblem, not its coupling scope.
- The moment-path twin (`ES_ALI_MOM_ON`, `USE_LAM_ES=F`) applies the
  identical formula to the `MOM_J_CMF_NONMONO_V1` result via
  `GET_MOM_ALO_DIAG`. Unreachable under KNOP (`USE_LAM_ES` forced T).

### 2.3 Convergence measurement

`ES_DJ = max_I |RJNEW - RJ_old| / |RJNEW|` per pass. Because
`RJNEW - RJ_old = (TC - RJ_old)/DEN`, the metric equals the raw
Lambda residual `|TC-RJ|/|RJ|` **divided by DEN** - i.e. the
preconditioned residual, amplified by up to 100x where
`DJDETA*chi_scat -> 1`. At the fixed point both vanish together, so
it is a valid stopping quantity for the accelerated iteration, but
it is a proxy: it is not the unaccelerated residual
`|J - Lambda[J]|/|J|`, and nothing in the loop evaluates that. The
stock Eddington-factor test (`|FOLD-FEDD|max > ACC_EDD_FAC`) remains
as the second exit criterion.

Inner-loop convergence does **not** imply outer convergence: ES_ALI
only contracts the per-frequency scattering subproblem. The
population/radiation outer iteration is untouched. Verified against
`adp2_tauto`: every sweep prints `freqs unconverged=0,
max final dJ/J~2.7e-5, ALO fallbacks=0` while the model still needed
hundreds of great iterations - the J subproblem was solved; the slow
mode lives downstream in the population step.

### 2.4 Fallbacks and edge cases

- `ALO_OK=F` (solver never stored a diagonal, e.g. `FG_OPT=DIFF` on
  V13): pass degrades to `RJ=TC`, `ES_ALI_FB` latches for the
  frequency. Exact stock behaviour.
- `RJNEW <= 0`: falls back to `TC` for that depth. Positivity
  preserved.
- `ES_ALI_MAX` clamped <= 20 at read time (`rd_control_variables.f:
  2618`); the stock 20-pass hard cap bounds the loop regardless.
  `ES_ALI_MAX <= 0` disables extra passes but the OAB update still
  applies each pass - undocumented edge, user-error class.
- `ES_ALI_TOL <= 0` forces max passes everywhere - same class.
- No lower bound on `DEN` beyond the 0.01 floor: where
  `DJDETA*chi_scat -> 1` the step is up to 100x the Lambda step.
  This is the designed strong-coupling regime; the positivity
  fallback and `ES_DJ` monitoring bound it. Under KNOP,
  `KN_NO_GAIN=T` (forced) keeps `chi >= chi_scat` so single-scatter
  albedo <= 1 backs the floor assumption. On the SOB+V13 route there
  is **no** equivalent albedo guarantee (`NEG_OPAC_OPT=ESEC_CHK`
  floors `CHI` at 0.1*ESEC, not at `CHI_SCAT`); `CHI < CHI_SCAT`
  cells could in principle push `DEN` negative - floored at 0.01,
  then positivity-fallback. Same exposure class as the original
  OAB scheme; flagged, not a demonstrated defect.
- `DO_SOB_FG_KNOP` is excluded from `ES_ALI_ON` (`.NOT.
  DO_SOB_FG_KNOP` in the gate). Inferred reason: that route lacks
  the `KN_NO_GAIN` albedo floor, so the 0.01 DEN floor's safety
  argument is weaker. Exclusion is conservative and correct.
- ES_ALI also runs inside JFNK residual/JVP eval sweeps
  (`JFNK_IN_EVAL` does not gate it): deterministic and consistent
  (lagged store frozen by `FREEZE_I_NEXT`), but each eval pays up to
  `ES_ALI_MAX` formal solves per continuum frequency - a real cost
  multiplier for `KNOP_opt`-class recipes.

### 2.5 Verdict

Mathematically correct diagonal-ALI update on the right operator,
right state, right ordering, right fallback. Runtime-verified on
this model family (`adp2_tauto`: 0 fallbacks, 0 unconverged, dJ/J
~3e-5 every sweep). Concerns are diagnostic-semantics and edge-case
hygiene, not correctness.

## 3. Pillar 2 - DO_ES_ALI_BA

### 3.1 Expression

`do_var_cont.f:428-438`, inside the `DO_NONMONO_KNOP .OR.
DO_NONMONO_SOB` continuum branch, once per frequency per BA build:

```
DENB(K) = 1 / max(1 - TB(K)*CHI_SCAT_CLUMP(K), ES_ALI_BA_MIN)
TA(K)  *= DENB(K)        ! stored dJ/dCHI diagonal
TB(K)  *= DENB(K)        ! stored dJ/dETA diagonal
```

`TA`/`TB` are the solver's stored bare-Lambda response diagonals
(dJ/dCHI < 0, dJ/dETA > 0 by construction sign convention);
`CHI_SCAT_CLUMP >= 0` is current-iteration opacity. For physical
stores `TB*chi_scat >= 0` so `DENB >= 1`: the Jacobian diagonal is
always strengthened, matching the ALI expectation that the bare
Lambda response under-corrects by the divergent factor
`1/(1-lambda_eff)` in scattering-dominated depths.

### 3.2 Checks requested

- **Signs/dimensions**: correct. `TB` is dJ/dETA per depth;
  `CHI_SCAT_CLUMP` is per depth, same grid; product dimensionless.
- **chi_scat -> 0**: `DENB -> 1`, amplification vanishes - recovers
  the unaccelerated Jacobian exactly.
- **Scattering dominated** (`TB*chi_scat -> 1`): `DENB` grows,
  capped at `1/ES_ALI_BA_MIN` = 10x default.
- **Singularity**: floored; `TB < 0` (unphysical store) would give
  `DENB < 1` - silent de-amplification, no guard. Bounded concern:
  the store comes from the same solver pass, signs are asserted by
  construction.
- **Frequency of application**: once per frequency per BA matrix
  build (`DO_VAR_CONT` under `SECTION='CONTINUUM'`,
  `COMPUTE_BA .AND. .NOT. LAMBDA_ITERATION`). The amplified Jacobian
  is then frozen inside the stored BA like every other entry -
  standard BA reuse semantics, not a per-iteration re-application.
- **Converged vs lagged field**: the ALO store reflects the **last
  (converged) J pass of that same frequency, same great iteration**
  - the solver zeroes and re-accumulates it on every call. Current,
  not lagged. `CHI_SCAT_CLUMP` likewise current.

### 3.3 Defects and concerns found

1. **No dependency on `DO_ES_ALI`** (confirmed gap, not a math
   error). The factor assumes the converged-ES response, but with
   `DO_ES_ALI=F` the J loop does `NUM_LAM`-style passes and the ES
   field may be less converged within the frequency. The amplified
   Jacobian is then a heuristic, not the exact derivative of the
   map actually iterated. E3 pairs the flags correctly; nothing
   enforces or warns about E2-without-E3 or E3-without-E2
   combinations. Recommend a VADAT-time informational warning, not a
   hard gate - the heuristic may still help.
2. **Floor asymmetry**: J loop floors `DEN` at 0.01 (up to 100x
   step); BA floors `DENB`'s denominator at `ES_ALI_BA_MIN=0.1`
   (up to 10x). If intentional (BA step more conservative) it is
   undocumented in `control_variable_mod.f:1937-1954`; if not, it is
   a 10x discrepancy between the two halves of the same approximation.
3. **Scope limit - Sobolev line Jacobian not covered** (confirmed).
   The flags modify only continuum-frequency `VJ` feeding
   `BA_UPDATE_V7`. The Sobolev line BA rows take their continuum
   coupling from the `SECTION='LINE'` `DO_VAR_CONT` call, which runs
   the stock `CONT_VEL=F` branch - `DO_ES_ALI_BA`/`TRI` never reach
   it. On an all-SOB model the line-rate Jacobian is where much of
   the population stiffness lives; the acceleration covers only the
   continuum channel. This bounds the expected E3 gain.
4. **`DO_FG13_ALO_TRI` silently inert on KNOP routes and
   `DO_SOB_FG_KNOP`**: `GET_FG13_ALO_TRI` reads only the V13 module,
   which never runs there; `TRI_OK=F` skips the fill with no
   warning, and the VADAT warning at `rd_control_variables.f:3770`
   exempts `DO_SOB_FG_KNOP`. Dead flag on those routes.
5. **Mixed-operator risk under `DO_NONMONO_SOB + USE_LAM_ES=F`**:
   `do_var_cont.f:411` takes the MOM diagonal while
   `GET_FG13_ALO_TRI` would still return V13 off-diagonals (V13 does
   run each frequency for ray moments). `DO_FG13_ALO_TRI=T` there
   yields MOM-diagonal + V13-offdiagonal with no warning. The
   warning block sits inside the `USE_LAM_ES=T` branch only. That
   route is currently non-viable anyway (moment path FAIL on m15),
   but the inconsistency is real.
6. **Pre-existing caveat (not flag-related)**: at `SECTION='LINE'`
   frequencies the continuum J and its Jacobian use stock monotonic
   solvers (`NEWJSOLD`/`JFEAU_IBC`) even under `DO_NONMONO_SOB` -
   the nonmono formal solver runs only at continuum frequencies.
   Recorded here because it bounds what any ALO flag can accelerate.

### 3.4 Verdict

Formula correct, signs correct, state correct, timing correct,
fallback-free by construction (inert when store invalid). The
concerns are coverage (lines untouched), silent-inert edges, and the
undocumented J-loop/BA floor asymmetry - all fixable with warnings
and docs, none blocking E3.

## 4. Pillar 3 - DO_BA_RES_ACC (+ DO_BA_RES_REUSE)

### 4.1 What is accelerated

Nothing is extrapolated. `DO_BA_RES_ACC` is a **guard**, not an
accelerator in the extrapolation sense: it judges each committed
population step and rolls back the ones that degrade the
statistical-equilibrium residual. Its companion `DO_BA_RES_REUSE`
recovers efficiency on the other side - while the judged residual
contracts it lets the stored BA matrix be recalled past the
`VAL_FIX_BA` demand window instead of forcing rebuilds. Together
they let the solver take full BA steps (rebuilding only on evidence
of staleness) instead of conservatively damping every step or
rebuilding blindly on a schedule.

### 4.2 The residual

Measured inside `SOLVEBA_V13` on every call (`solveba_v13.f:290-512`),
assembled from `SE(ID)%STEQ` rows (per-ion rate-balance residuals),
`STEQ_ADV` (ionization residual incl. advective terms), `STEQ_ED`
and `STEQ_T`/`STEQ_T_EHB` when NE/T are not held fixed:

- `BA_SE_RES` = global max |residual| over species x depth x eqn.
- `COL_RES(I)` = per-depth max; `BAACC_COL_PREV/MIN` track
  pre-step and best-seen per column.
- `BA_SE_RES2` = K-th largest column residual, `K=MIN(ND,MAX(4,ND/10))`
  - a rank-based bulk channel that ignores the K-1 worst depth columns
  (added 30-Sep-2026 after a degenerate ~1E262 inner-boundary floor
  pinned the max channel blind).

This is the SE residual **under the current fixed radiation field**
- `J` does not move during `SOLVEBA`. It is exactly the right metric
for judging a population step: the step is a Newton-like update on
`SE(pops; J_fixed) = 0`, so residual growth measures step quality.

### 4.3 Acceptance and rollback

Arming (`solveba_v13.f:971-1060`): on a non-LAMBDA call under
`DO_BA_RES_ACC`, `POPS_BA_SAVE=POPS` is snapshotted and
`BAACC_PENDING` set with the pre-step `BA_SE_RES`/`BA_SE_RES2`/
`COL_RES` recorded. LAMBDA calls still allocate (transient-failure
healing) but do not arm.

Judgement (`:593-668`): skipped while `GS_POPS_HELD` (deferred,
pending stays armed). Per-column failures populate `BA_RBK_MASK`;
then `REJ_BULK` (RES2 channel) vs `REJ_MAX` (scalar max, exempted
when the pre-step residual is already `>= BA_RES_SAT=1E30`) decide
whole-state vs masked rollback via `BA_REJ_FRAC` (default 0.25).
`BA_JUDGE_RATIO = BA_SE_RES2/BAACC_RES2_PREV` - the operative ratio
is the **rank channel**, not the scalar max.

- Per-column masked test: column I fails when
  `COL_RES(I)/BA_RES_GROW(=1.5) > COL_PREV(I)` or
  `/BA_RES_HARD(=1E2) > COL_MIN(I)`; `BA_RES_SAT` exempts columns
  already at degenerate residual; non-finite always fails.
- **Masked/partial reject**: only `BA_RBK_MASK` columns restored
  from `POPS_BA_SAVE`, their fresh STEQ zeroed, passing columns keep
  the update (`BA_PART_THIS`); the lagged radiation fields are also
  rolled back (`RESTORE_FG_INEXT_LATEST` under KNOP,
  `RESTORE_MOM_JNEXT_LATEST` under KNOP-or-SOB moment paths);
  worst equation per masked column is quarantined. Does not demote.
- **Full reject**: `POPS=POPS_BA_SAVE` (or `POPS_GS_SAVE` revival
  if held), lagged fields rolled back, `STEQ=0`, `MAXCH` capped to
  `[100, 1.5E5]`, `BA_STEP_REJECTED` returned.
- **LAMBDA follow-ups are soft-judged, not rolled back**: a degraded
  residual halves `BA_TRUST_CUR` and increments `BA_NREJ_CONSEC`
  toward `BA_NREJ_LIM` (=20 consecutive -> STOP); a progressing step
  resets the streak but holds trust. Only bounded `BA_PROBE_LAM`
  trials inside `RD_LAMBDA` roll back (restoring pops + radiation
  state and skipping the stale Lambda solve).
- Trust region: `BA_TRUST_CUR` (init `BA_TRUST=1E4`) is a per-column
  demand cap inside SOLVEBA - a column whose `MAXVAL(ABS(STEQ(:,I)))`
  exceeds it is zeroed before the solve; halves on reject, x1.25 on
  accept, floor 1.0. This is the mechanism that actually limits
  correction size; `MAXCH`/`MAXCH_APP` are the reported demands.
- Deferred judgement while `GS_POPS_HELD`; fail-closed demotion to
  damped LAMBDA when `BA_SNAPSHOT_MISSING` (alloc failure).
- `CONV_EXIT_MET` requires `.NOT. BAACC_PENDING` - an unjudged armed
  step blocks a converged exit, so resacc can never be bypassed at
  the finish line.
- JFNK disarms a pending arm on entry (`jfnk_solve_corr_v1.f:547-555`)
  because the post-Newton residual cannot judge the displaced BA
  trial; `JFNK_SAVE/RESTORE_MOD_STEQ_V1` brackets perturbed residual
  evals so resacc sees outer-iteration residuals, not leaked trial
  values.
- `NG_CLEAN` resets on any `BA_REJ_THIS`/`BA_PART_WAS`/`BA_SOFT_THIS`/
  pending event - rejected steps break the NG clean window.
- Populations are restored from an exact pre-step snapshot, so
  particle/charge conservation of the committed state is preserved
  by construction; masked columns revert per-column.

### 4.3a DO_BA_RES_REUSE (the recall half)

`solve_for_pops.f:3242-3295`, gated `DO_BA_RES_REUSE .AND.
DO_BA_RES_ACC` (SOB auto-sets both for absent keys; KNOP auto-sets
ACC only - deliberate asymmetry):

- While the judged residual still contracts (`BA_JUDGE_RATIO <=
  BA_REUSE_TOL=0.9`, RES2 channel), BAMAT recall is allowed past the
  `VAL_FIX_BA` demand window - the stored linearization keeps being
  used while residual evidence says it is still good.
- `BA_REUSE_STALL` counts consecutive judged recalls with ratio >
  TOL; `>= BA_REUSE_NSTALL=2` forces `COMPUTE_BA=T` rebuild
  (`ba_reuse_stall`). Only recall-steps count - a flat judgment on a
  fresh matrix cannot trigger it.
- Recall is a hard `COMPUTE_BA=L_FALSE`; rejection, partial reject,
  `SUM_FIRED`, missing/unstored matrix, `WRBAMAT_RDIN=F`,
  `CON_SCL_T!=0`, LAMBDA next-call, and JFNK all veto or preempt it.
- Intended semantic: residual evidence can override the legacy
  `MAXCH_CTRL > VAL_FIX_BA` rebuild demand within the stall budget -
  demand-based rebuild pressure is subordinate to residual evidence,
  with the stall counter as backstop. Not a veto of the rejection or
  `SUM_FIRED` rebuilds.

### 4.4 Compatibility with ES_ALI / ES_ALI_BA

Orthogonal layers targeting the same fixed point:

- ES_ALI changes the J the populations see (better ES contraction).
- ES_ALI_BA changes the continuum Jacobian diagonal (bigger, better
  population corrections).
- BA_RES_ACC judges whatever step results and rolls back failures.

No mathematical conflict: all three are consistent with
`J = Lambda[eta + chi_scat J]` and `SE(pops; J) = 0`. The interaction
risk is dynamical, not structural: ES_ALI_BA's amplified Jacobian
produces larger corrections that may trip `BA_RES_GROW` more often -
which is precisely what the rollback exists for, at the cost of
rejected-step overhead (a wasted SOLVEBA call plus forced rebuild).
`DO_BA_RES_REUSE` (auto-set with RES_ACC under SOB) mitigates the
rebuild cost via residual-driven matrix recall.

### 4.5 Concerns

- The residual is a **max-norm** on the worst equation, and the
  operative reuse/judgement ratio rides the **RES2 rank channel**,
  which deliberately ignores the worst `K-1` (`K=ND/10`) depth
  columns. On models with a degenerate floor (the ~1E262
  inner-boundary rows seen on Dessart runs) that is the right call -
  the code learned it twice (01-Oct and 30-Sep-2026 fixes) - but it
  is a tunable blind spot: a rebound confined to fewer than K
  columns is invisible to RES2 (the masked per-column test still
  catches it for rollback, but `BA_JUDGE_RATIO`/reuse recall do not
  see it).
- `BA_RES_SAT` exempts the scalar-max channel at degenerate
  residuals; RES2 remains the bulk veto. Any new degenerate
  equation class could blind both channels again.
- Judgement is one call late (arm on call N, judge on N+1); held
  calls defer, JFNK calls disarm. A pending verdict also blocks
  `CONV_EXIT_MET`, so the mechanism cannot be dodged at convergence
  - correct, but it means the last armed step always costs one extra
  call.
- Rejection rolls back **both** populations and the lagged radiation
  field (`RESTORE_FG_INEXT_LATEST`/`RESTORE_MOM_JNEXT_LATEST`) on
  the nonmono paths - stronger consistency than population-only
  rollback, at the cost of discarding a radiation sweep's work on
  every rejected step.
- `BA_NREJ_LIM=20` consecutive rejections is a hard `STOP` - a
  rejection storm does not just demote, it can kill the run. Bounded
  and intended, but E4 telemetry must watch `BA_NREJ_CONSEC`.
- During rapid radiation-field change (early iterations, or after an
  ES_ALI-driven J jump), legitimate population steps can grow the
  SE residual transiently and be rejected - self-throttling by
  design, but it can convert a would-be-accepted damped step into a
  LAMBDA demotion cycle. `BA_TRUST_CUR` halves per rejection and
  only recovers x1.25 per acceptance; a run of rejections can park
  the solver in LAMBDA mode for many iterations.
- `DO_BA_RES_ACC=F` is an **explicit override** in all three audited
  VADATs - the SOB auto-set intended T. Whatever motivated the
  override (likely the rejection-storm behaviour above on earlier
  recipes) should be re-tested rather than assumed: the masked
  rollback and RES2 channel postdate those experiments.

### 4.6 Verdict

Correct, heavily engineered, fail-closed rollback machinery -
arguably over-built relative to its job, which is the right side to
err on for a guard that can STOP the run. It does not accelerate the
map; it makes aggressive steps safe, and `DO_BA_RES_REUSE` recovers
efficiency by letting judged-contracted matrices be recalled past
the demand window. Rejection restores populations **and** the lagged
radiation field, so a rejected step leaves no half-committed state.
Compatible with the ES_ALI stack. The main risks are rejection-storm
demotion cycles and the `BA_NREJ_LIM` hard stop - both bounded,
both visible in ITER_LEDGER/OUTGEN telemetry.

## 5. Pillar 4 - Non-monotonic velocity consistency

- **Operator consistency**: every ALO diagonal used by ES_ALI /
  ES_ALI_BA is the stored response of the **actual formal solver
  that ran that frequency** - `GET_FG13_ALO_DIAG` (V13, clamped
  SIGMA) on this recipe, `GET_FG_ALO_DIAG` (Knop solver) under
  KNOP/`DO_SOB_FG_KNOP`, `GET_MOM_ALO_DIAG` on the moment path.
  The accessor dispatch mirrors the solver dispatch predicate
  exactly (the 02-Nov-2026 fix; before it, KNOP silently fell back
  to plain Lambda 56377/56377 times). No monotonic-flow assumption
  enters anywhere: the diagonals are accumulated from the solver's
  own ray-integral coefficients (`dAV/dS'/chi_hat` per ray node).
- **Sobolev lines**: `SUB_SOB_LINE_V3`/`SOBEW` use `SIGMA_TRUE`
  (unclamped) with the documented angular-factor floor - untouched
  by all three flags.
- **What the ALO does NOT capture** (inherent limitation, same as
  baseline): the operators are depth-diagonal (TRI adds K+-1 depth
  bands only) and per-frequency diagonal - the true Lambda has
  off-diagonal frequency blocks from velocity redistribution and
  multiple resonance points along a ray in non-monotonic flow.
  ES_ALI's per-frequency Jacobi therefore still leaves
  inter-frequency coupling lagged. This is not a regression: the
  stock `RJ=TC` iteration has the identical coupling scope; ES_ALI
  strictly improves the contraction rate of the part it does cover.
- **Clamped vs true SIGMA** (pre-existing, documented): V13 runs on
  clamped `SIGMA >= -0.999` while Sobolev line depths use
  `SIGMA_TRUE`; the V13 ray moments feeding `K_ON_J` inherit the
  clamp (recorded as second-order in the 13-Oct audit). The ALO
  diagonal inherits whatever the solver used - self-consistent with
  the operator actually iterated, which is the correct choice for a
  preconditioner.
- **No new monotonic approximation introduced**: verified - all
  amplification factors derive from stored response coefficients of
  the nonmono solvers; nothing assumes `dV/dr > 0` or single
  resonance per ray.

## 6. Pillar 5 - Warm restart consistency

Audited `run_model_grid.py:2864-3205` (`prepare_cmfgen_directory`)
for the prospective `FRESH_sob -> FRESH_sob_dieoff` handoff.

**`warm_restart_population_only: true` is SAFE for this pair:**

- Stages `POINT1 POINT2 SCRTEMP MODEL` byte-verbatim; population-only
  additionally rewrites the checkpoint T slot to the target hydro T
  (log-log interpolation on the checkpoint radius grid, with
  hydro-ND and radius-range guards).
- Fingerprints pass: ND=70, NT=1839 identical; population-only
  requires byte-identical `SN_HYDRO_DATA` (hard error otherwise) -
  satisfied here.
- `VADAT`, `IN_ITS`, `MODEL_SPEC` are **regenerated from the target
  YAML** - verified the dieoff VADAT carries all `DIE_*=F,F` plus
  the target's solver keys. Changed dielectronic rates are used when
  BA matrices rebuild (`BAMAT*`/LU caches live in `SCRATCH/`, wiped
  by `rmtree`; `COMPUTE_BA=.TRUE.` on first call - no contaminated
  Jacobian).
- Radiation field cold-starts: `EDDFACTOR*` is never staged
  cross-dir, `COMP_F=T` rebuilds the Eddington table, `RJ` zeroes on
  iteration 1 - expect the documented first-iteration J transient.
  `FGINEXT_GS`/`MOMJNXT_GS` are full-mode-only and absent in
  FRESH_sob anyway (SOB+LAM_ES=T formal solver).
- Acceleration histories cannot contaminate: NG/JFNK/AA state is
  runtime-only, never checkpointed. `NITSF` resumes at 722 and
  `LST_NG=710`, so NG scheduling resumes (`NEXT_NG ~ it726`) -
  appropriate for near-converged populations; mode selection stays
  demand-driven, not counter-driven.
- No convergence flags carry: `ITER_*`, `OUTGEN`, markers all wiped;
  `NEWMOD=.FALSE.` only because the checkpoint read succeeds -
  populations come from `SCRTEMP`, not LTE estimates.

**Guard gaps (edge cases, not blockers for this pair):**

- NT check compares only the integer count and is **skipped entirely
  on first submission** (no target MODEL yet); two different atomic
  models with equal NT would misalign `POPS(NT,ND)` silently. No
  species-name fingerprint exists.
- ND check degrades to warn-only when target hydro is absent at
  prepare time; full `warm_restart_from` checks ND count only - a
  same-ND different hydro passes silently and would freeze
  checkpoint T under `FIX_T`.
- No requirement that the source run converged; copying `SCRTEMP`
  from a *running* job risks a torn record (mitigated by
  `SCR_RITE_V2`'s `WRITFAIL` guard; quiescent source safest).
- Full mode carries `RVTJ`/`CUR_MODEL_DATA`/`FGINEXT_GS`/
  `MOMJNXT_GS` if present - for a *changed-physics* handoff those
  GS dumps would be stale radiation state; `population_only` avoids
  this entirely. **Recommendation: use `population_only: true` for
  all E1-E4 arms.**

## 7. Pillar 6 - Controlled experiment matrix

All four arms identical except the three keys; all warm-started
from the same `FRESH_sob` checkpoint via
`warm_restart_from` + `warm_restart_population_only: true`. The
in-flight cold `dieoff` run is the cold-start E1-equivalent baseline
(DIE flags off, all accelerators off) - let it finish; it costs
nothing extra and anchors the cold-start comparison.

| Arm | DO_ES_ALI | DO_ES_ALI_BA | DO_BA_RES_ACC | Isolates |
|---|---|---|---|---|
| E1w | F | F | F | warm baseline (control) |
| E2 | T (`ES_ALI_TOL=1e-4`, `ES_ALI_MAX=20`) | F | F | J-loop ES contraction alone |
| E3 | T | T (`ES_ALI_BA_MIN=0.1`) | F | + continuum Jacobian amplification |
| E4 | T | T | T (+auto `DO_BA_RES_REUSE`) | + residual-judged rollback safety |

Do **not** run E3 without E2 semantics confusion: E3 keeps
`DO_ES_ALI=T`, so the flags are cumulative by design. If E2 shows
ES_ALI inert on this recipe (fallbacks or `ES_DJ` never binding),
E3's amplification is unvalidated heuristic - check the OUTGEN
`ES_ALI:` summary line first.

### Record per arm

- Outer iteration count to `EPS_TERM` (same physical tolerance -
  do not relax).
- Wall-clock time and per-iteration mean (TIMING).
- Formal-solution evaluation count: proxy = iterations x continuum
  freqs x mean J-loop passes; ES_ALI adds passes - the OUTGEN
  `ES_ALI:` line gives unconverged/fallback counts per sweep.
- `maxraw`, `maxapp`, `maxatt`, `se`, `se2` from ITER_LEDGER.
- Radiation-field residual: `fg_rfd`/`fg_churnf`/`fg_cap` ledger
  columns + J_COMP drift diagnostics.
- Population residual: `BA_SE_RES`/`BA_SE_RES2` ITER_SUM channels
  (measured every call regardless of `DO_BA_RES_ACC` since
  15-Oct-2026).
- ES source-function residual: `ES_ALI_NUNC/NFALL/DJMAX` OUTGEN
  counters (E2+).
- Energy/particle/charge conservation: existing OUTGEN diagnostics.
- Rejected acceleration steps: `ba_rejected`/`ba_partial_rej`/
  `ba_reuse_stall` ITER_REASON counts, `BA_NREJ_CONSEC`/`BA_NREJ_SOFT`,
  `BA_TRUST_CUR`, `ba_build`/`ba_reuse` ledger columns (E4);
  `mom_*`/`ng` ledger columns for the armed-but-inert NG.
- Warm-start transient: first-iteration J spike magnitude and
  iterations to return below the cold run's it200 level.

### Predicted outcome (to be tested, not assumed)

ES_ALI should cut the per-frequency ES subproblem to ~2-3 passes
(proven on `adp2_tauto`), reducing per-iteration wall time if the
J-loop was pass-limited, but the **iteration-count** bottleneck is
the population/radiation outer map - ES_ALI alone (E2) may not move
MAXAPP's tail contraction much. E3's amplified continuum Jacobian is
the arm most likely to attack the 0.9973/it tail, bounded by its
10x cap and the fact that Sobolev line rows are untouched. E4 adds
safety for the larger corrections E3 produces; expect some rejection
traffic early. If E3/E4 show no tail improvement, the residual is
not in the continuum-ES channel - next suspect is the Sobolev line
Jacobian's own lagged coupling (outside all three flags' scope) or
genuine operator conditioning, which points to JFNK-CPL or grid
methods rather than more ALI.

## 8. Pillar 7 - Regression tests

Existing coverage:

- `unit_tests/python/test_dvc.py` (`esali`, `tri` cases) +
  `drv_dvc.f90`: verifies `DENB = 1/max(1-DJDETA*CHI_SCAT,
  ES_ALI_BA_MIN)` and K+-1 band fills against hand-computed values.
- `unit_tests/python/test_fg13.py` + `drv_fg13.f90` +
  `stub_fg13_alo.f90`: V13 ALO store vs finite-difference response.
- `test_iter_ledger.py`, `drv_iter_ledger.f90`: ledger schema.
- Production evidence: `adp2_tauto` (ES_ALI on, 0 fallbacks),
  `KNOP_esA/esB` A/B configs exist but runtime verification pending.

Gaps to close (proposed, not implemented):

1. **Zero-scattering limit**: `chi_scat=0` -> `DEN=DENB=1` ->
   `RJNEW=TC`, Jacobian unamplified. Unit-testable in `drv_dvc`
   style for the BA side; the J-loop side needs a `comp_j_blank`
   harness or a production A/B.
2. **Scattering-dominated limit**: `DJDETA*chi_scat -> 1` -> floor
   behaviour (0.01 J-loop, `ES_ALI_BA_MIN` BA) and positivity
   fallback on `RJNEW<=0`.
3. **Weak/strong scattering continuum**: sweep `chi_scat/chi` ratio;
   verify monotone amplification `DENB>=1` for physical stores.
4. **Rapidly changing field**: ES_ALI inner loop must still exit on
   `ES_ALI_MAX` when `ES_DJ` stays above tol - bounded passes, no
   hang. Covered implicitly by the 20-cap; worth an explicit count.
5. **Non-monotonic multi-resonance**: the ALO is depth-diagonal by
   design; a regression test should confirm no assumption of single
   resonance enters - i.e. the update uses only same-depth stored
   response (true by construction; test guards future edits).
6. **Warm restart with modified rates**: pipeline-level test that
   `population_only` staging + regenerated VADAT produces a run whose
   first BA build uses new rates (assert `COMPUTE_BA=T` on it1 and
   no `BAMAT` reuse across the handoff).
7. **Compiled coverage of the `comp_j_blank` OAB block**: currently
   only production-run evidence exists; a `drv_*` harness stubbing
   `FG_J_CMF_V13`/`GET_FG13_ALO_DIAG` would pin the update, fallback,
   and exit logic.

## 9. Confirmed defects vs concerns

### Confirmed defects (documentation/guard-level; no solver-math error)

- D1. `control_variable_mod.f:1919-1931` comment describes
  `DO_ES_ALI` as SOB+V13-only; omits the KNOP extension and the
  `ES_ALI_MOM_ON` moment branch. Stale.
- D2. `check_param_consistency.f:129` comment claims
  `DO_SOB_FG_KNOP` forces `USE_LAM_ES=T`; no such forcing exists.
- D3. `DO_FG13_ALO_TRI` silently inert on KNOP and `DO_SOB_FG_KNOP`
  routes; the VADAT warning exempts the latter - a user gets no
  warning and no effect.
- D4. `DO_FG13_ALO_TRI=T` under `DO_NONMONO_SOB + USE_LAM_ES=F`
  mixes MOM diagonal with V13 off-diagonals unwarned.
- D5. `DO_ES_ALI_BA` has no VADAT-time dependency warning: silently
  inert unless an ALO route returns valid, and silently heuristic
  when `DO_ES_ALI=F`.
- D6. J-loop DEN floor (0.01) vs BA `ES_ALI_BA_MIN` (0.1) asymmetry
  undocumented.

### Theoretical concerns (bounded, by design or inherent)

- C1. `ES_DJ` is the preconditioned residual proxy, not the raw
  `|J-Lambda[J]|` residual; exit semantics should be documented as
  such.
- C2. No albedo<=1 guarantee on the SOB+V13 route (KN_NO_GAIN is
  KNOP-only); mitigated by floor + positivity fallback.
- C3. Inter-frequency Lambda coupling stays lagged under all three
  flags - inherent to the per-frequency Jacobi scope, same as
  baseline.
- C4. ES_ALI multiplies JFNK eval cost by up to `ES_ALI_MAX` passes
  per frequency.
- C5. BA_RES_ACC rejection storms can park the solver in LAMBDA
  demotion and, at `BA_NREJ_LIM=20` consecutive, hard-STOP the run;
  bounded by `BA_TRUST_CUR`/trust recovery but worth watching in E4
  telemetry (`ba_rejected`, `BA_NREJ_CONSEC`, `trust` columns).
- C6. Warm-restart guard gaps (NT skipped on first submission, no
  species fingerprint, full-mode hydro check is count-only, no
  source-convergence requirement).
- C7. `DO_BA_RES_REUSE` recall rides the RES2 rank channel, which
  ignores the worst `K-1` depth columns by design; rebinds confined
  below that rank are invisible to `BA_JUDGE_RATIO` (masked rollback
  still catches them for the step itself).

### Explicitly verified clean

- OAB update algebra, accessor/solver dispatch parity, positivity
  fallback, stock-path bit-identity when flags are F.
- BA_RES_ACC residual definition (SE balance under fixed J, measured
  every call), masked + bulk rollback correctness including
  lagged-field restore, conservation preservation via exact
  snapshot, fail-closed behaviour on missing snapshot, JFNK
  arm-disarm isolation, `CONV_EXIT_MET` pending-verdict gate,
  `DO_BA_RES_REUSE` stall/rebuild dominance ordering.
- Warm-restart `population_only` for this pair: correct files,
  correct fingerprints, target physics authoritative, no stale
  radiation/rate state.

## 10. Recommended sequence

1. Submit E1w-E4 as warm-started clones (population_only) - four
   jobs, identical checkpoint, one-variable-at-a-time.
2. Read E2's `ES_ALI:` OUTGEN lines after ~50 iterations: if
   `ALO fallbacks` is nonzero or `freqs unconverged` stays high,
   the V13 diagonal is not serving this recipe and E3/E4 are
   testing a heuristic - record and proceed anyway (the flags are
   independent).
3. Compare tail contraction (`maxapp` slope it200+) and total
   iterations to EPS_TERM, not early-iteration behaviour - the warm
   transient is common to all arms.
4. Only after E-matrix results: fix D1-D6 (comment/guard-level,
   no solver change), then consider the Sobolev-line Jacobian
   coverage gap (section 3.3.3) as the next real acceleration
   target - it is outside all three audited flags.
