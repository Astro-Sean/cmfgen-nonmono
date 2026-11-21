# Convergence Theory Audit - Non-Monotonic CMFGEN

Source audited: `cmfgen_jfnk/` at branch `main`, commit
`cfae58528648ed138749f8748e195226ab03b0f3` (re-confirmed for this
revision; tree clean except the other agent's untracked
`docs/audits/FG_FIELD_DRIFT_AUDIT.md`). No code changed; no
compilation run.

REVISED after the critical review: all line references re-verified
at this commit; the Type-I/Type-II confusion corrected; a new
section "CRITICAL REVIEW OF LITERATURE-TO-CMFGEN ASSUMPTIONS" and a
full source-claim verification table added.

## 1. The equations the code solves

### 1.1 Inner FG fixed point (Level 1)

Under `DO_NONMONO_KNOP` (and `DO_SOB_FG_KNOP`) the formal path
stores a lagged intensity field `I_NEXT`. One great iteration:

1. In-sweep: `FG_J_CMF_NONMONO_V1` solves each ray once, then
   blends the fresh intensity into the store
   `I_NEXT <- (1 - INEXT_RELAX)*I_NEXT + INEXT_RELAX*I_fresh`
   (`WRITE_STORE`, default `INEXT_RELAX = 0.5`, adaptive when
   `KN_ADAPT_RELAX=T`). The relaxed store is what opacities,
   emissivities and SE rates see.
2. Refresh: `FG_KNOP_REFRESH_V1` (gated by `DO_KNOP_REFRESH`,
   `SOLUTION_METHOD='INTEGRAL'`) rebuilds monotone-cubic chi/eta
   coefficient packs `C` from the captured sweep fields, then per
   ray solves the frozen-coefficient fixed point

   ```
   I_ray* = G_ray(I_ray*; C)          (eq. FP-FG)
   ```

   in `KNOP_RAY_SOLVE_V1` via ordered passes `SWEEP_FLD_ORD`
   (two calls per ray solve: `I_M` at line 7126, `I_P` at 7181 -
   one direction field each). `G_ray` is the short-characteristics
   march with lagged C-corner coupling, monotone-cubic chi/eta,
   positivity guards (`VNEW<0` -> `1e-30`, lines 8159-8162),
   `NONMONO_INEXT_CAP` ceiling (8163-8166), and per-cell
   `sign(GAM)` frequency ordering (8135-8139).

The refresh is a successive-substitution iteration on the
frozen-coefficient formal-solution map, marched as a physically
ordered Gauss-Seidel-type sweep. It is NOT a Lambda iteration
(populations never enter) and NOT ALI (no approximate-operator
split of a coupled SE+RT map). See Sect. 5.

### 1.2 Moment-side analogue

`MOM_J_CMF_NONMONO_V1` keeps `J_NEXT`/`H_NEXT` under
`USE_LAM_ES=F`, with `MOM_KNOP_REFRESH_V1`, `JNEXT_RELAX = 0.5`
write-back, inner AA (`MOM_AA_M = 3`) and the outer-map AA
`MOM_AA_OUT` (default T, `MOM_AA_OUT_M = 8`).

### 1.3 The relaxed-map subtlety

The write-back makes the effective inner map the relaxed one
`x_{k+1} = x_k + omega*(G(x_k) - x_k)`, `omega = OMEGA_RAY`
inside the solve, `INEXT_RELAX`/`JNEXT_RELAX` at the store.
The AA history in `SWEEP_FLD_ORD` stores relaxed deltas
(`AA_RH = omega*(G(x)-x)`, line 8432) while `RES`/`AA_REF`/
`RATIO_OM` measure the raw residual `|VNEW - VOLD|`
(lines 8168-8173). The two conventions are consistent because
the omega-rescale (8558-8562) keeps `AA_RH` in current-omega
units and `RES` is omega-independent by construction.

## 2. Residual definitions actually in the code

Five different numbers; mixing them was a defect source:

| Name | Definition | Where | What it measures |
|---|---|---|---|
| `RES` / `KN_RF_RESID` | `max abs(VNEW - VOLD)/max(VOLD, 1e-4*col_peak)` over significant cells; `VNEW` is the RAW map output | `SWEEP_FLD_ORD` 8168-8173 (inward), symmetric outward block; per-pass history `KN_RES_PASS` 7558-7563 (per-pass MAX over rays, OMP CRITICAL) | True UNDAMPED fixed-point residual of FP-FG (max-norm) |
| `AA_RH` entries | `FLD_relaxed - x_k = omega*(G(x_k) - x_k)` | `SWEEP_FLD_ORD` 8431-8432 | RELAXED-map residual (omega units), used only inside AA fit |
| Sweep drift `FG_SW_MAXCHG` / `GS` churn | `|I_fresh - I_stored|/max(|I_stored|,floor)` | `WRITE_STORE` write-back; snapshotted before `FG_GS_MAXCHG` takes the refresh value (6589-6600) | Store churn between sweeps, DAMPED by `INEXT_RELAX` |
| `MOM_GS_MAXCHG` | fresh J solve vs stored `J_NEXT`, `max(|old|,1e-4*col_peak)` | `MOM` sweep | Same class as FG sweep drift |
| Population `MAXCH` | max demanded relative population correction | `SOLVEBA`/JFNK commit | Outer-map residual proxy (demand) |

Notes:

- The `1e-4` column-peak floor is a significance filter, not an
  acceptance floor; `NTOT_CUR`/`GS_NTOT` count cells above it
  (8175-8179). Weak-intensity errors below the floor are not
  tested - acceptable for gates, must stay out of physics
  acceptance criteria.
- `KN_RES_PASS(NPO) = max over rays at pass NPO` (7561) - a
  per-pass global-max ENVELOPE, not a per-ray history. Ratios of
  its entries are valid as a contraction estimate only if the
  argmax ray is stable across passes (the `FG_ARG_*`/`cycle`
  classifier exists to check that).

## 3. The three convergence levels as coded

- Level 1 (FG inner): per-ray `RES < KN_RAY_TOL` (2e-3),
  bulk exit `KN_RAY_BULK_TOL`, `KN_RAY_MAXPASS = 60`, commit/
  reject via `FG_COMMIT_CONV`, `FG_COMMIT_JUMP = 50`,
  `FG_COMMIT_NC = 0.05`; verdicts and `IXIT` codes in
  `KNOP_RAY_SOLVE_V1`/`FG_KNOP_REFRESH_V1`.
- Level 2 (coupled RT+SE): `GS_FREEZE` maturity controller in
  `solve_for_pops.f` - populations Lambda-held (`GS_HOLD_POPS`
  snapshot/restore) until refresh residual AND store churn pass
  (`GS_MATURE_TOL`, `GS_UNSAFE_TOL`, `GS_REL_STABLE` consecutive
  count). Knop Paper II ordering rule implemented.
- Level 3 (outer): `MAXCH < EPS` + non-Lambda/non-NG guards +
  Stage-B trace polish (`BA_TRACE_EVAL`) + `BA_EX_MAXCH` floor.

How they are mixed today:

- The refresh runs AFTER the in-sweep field already fed SE
  residuals and the BA Jacobian for that iteration - Level-2
  data is evaluated on the in-sweep relaxed field, Level-1
  certified a posteriori; `GS_FREEZE` then blocks the Level-3
  update on immature Level-1 evidence. Ordering correct, but
  one iteration can hold three field states: in-sweep relaxed
  `I_NEXT`, refresh-converged `I_ray`, committed retained store.
- `GS_NBAD` (Level-1 census) and `GS_NCONV` (Level-2 census) are
  checked separately in the gate - confirmed correct.
- JFNK writes into Level 3 through `JFNK_UNCONVERGED_COMMIT`
  and the field-maturity admission tier; frozen-J vs coupled
  modes (`JFNK_COUPLED_THIS`) differ in which field the
  residual `F` is evaluated on.

## 4. Relaxation controller analysis

### 4.1 Per-ray `OMEGA_RAY` inside `SWEEP_FLD_ORD`

- `RATIO_OM = RES/RES_PREV` (8516 ordered; `RES_PREV(IDIR)`
  per-direction in the legacy alternating `SWEEP_FLD`, 7592).
  Same ray, same pass-sequence metric - a consistent ray-local
  contraction diagnostic.
- Rules (8539-8553): `RATIO_OM > 1.05` -> `OMEGA_RAY *= 0.7`;
  `RATIO_OM < 0.5` -> `OMEGA_RAY *= 1.15`; `NOSC >= 2` sign flips
  -> `OMEGA_RAY *= 0.7`. Bounds `[0.45, 1.0]`.
- So the in-solve adaptive relaxation IS based on measured
  contraction (same-map residual ratio), supplemented by a
  sign-flip oscillation trigger. Answer to audit question H:
  measured contraction AND oscillation diagnostics.
- The `lambda ~ -1` argument in the code comment (a -1 mode
  under `omega<1` gives `|1-2*omega|<1`) is mathematically
  correct for the relaxed affine map.

### 4.2 Cross-call `INEXT_RELAX` (`FG_ADAPT_OMEGA_V2`, 6538-6588)

Two estimators:

- V2 (preferred): median of the last ~4 consecutive ratios
  `KN_RES_PASS(i)/KN_RES_PASS(i-1)`, then
  `INEXT_RELAX <- 0.5*INEXT_RELAX + 0.5*min(1,max(0.05,1/(1+rho)))`.
  `KN_RES_PASS` is the per-pass GLOBAL-MAX envelope - the ratio
  is a contraction only when the argmax ray is stable.
  Otherwise it mixes different rays' maps (spec violation:
  "no contraction from unrelated ray maxima"). The
  argmax-stability classifiers exist but do not gate this.
- Fallback: `T2 = KN_RF_RESID/KN_RF_RESID_PRV` (6582-6587) -
  global-max refresh-residual ratio across calls; same
  argmax caveat PLUS cross-coefficient-generation caveat
  (different `G` on the two sides).

### 4.3 MOM `JNEXT_RELAX` (CONFIRMED DEFECT C1)

`plane/mom_j_cmf_nonmono_v1.f` ~6931:
`DTR = MOM_GS_MAXCHG / MOM_RF_RESID_PRV`,
`JNEXT_RELAX = 0.5*MIN(1,MAX(0.05,1/(1+DTR)))`.
`MOM_GS_MAXCHG` is the CURRENT sweep drift (fresh solve vs
stored `J_NEXT`) while `MOM_RF_RESID_PRV` is the PREVIOUS
refresh max residual - different metrics, different
normalization, different pipeline points. Not a contraction
ratio. CONFIRMED DEFICIENCY.

### 4.4 Update-vs-residual confusion check

The Level-1 exit uses `RES` (raw `|VNEW - VOLD|`) directly -
undamped. Sweep drift metrics are churn gates, not convergence
tests. No instance found of a damped update alone releasing a
gate. PASS.

## 5. Classification of the refresh (corrected)

Per Sect. 10.1 of the literature review:

- NOT classical Lambda iteration (populations frozen;
  `S = eta/chi` fixed).
- NOT ALI (no `Lambda*` splitting; the lag is inside the
  discretization's corner-source terms).
- NOT operator splitting in the ALI sense.
- IS a successive-substitution fixed-point iteration on the
  frozen-coefficient formal-solution map, with
  Gauss-Seidel-type causal ordering in space (march along the
  characteristic) and frequency (march along `sign(GAM)`, the
  local direction of the `gamma*dI/dnu` advection term - Knop
  Paper II, A&A 501, 813, Sect. 3).
- The effective iterated map is `G_omega` (relaxed); its affine
  factor is `I - omega*(I - M)` where `M` is the raw-map
  derivative.

## 6. Acceleration audit vs the literature

### 6.1 `SWEEP_FLD_ORD` inner AA (FG, opt-in) - TYPE-II confirmed

Equations implemented (verified at 7949-8823; see Literature
Review Sect. 3 for the full derivation): weighted least squares
`min_gamma ||W*(r_1 - D*gamma)||` on residual differences,
update `x_mix = sum alpha_j*G_omega(x_{k-j})`,
`alpha_0 = 1 - sum gamma`. Type-II, m=2, safeguards as listed.

Gaps (corrected against Zhang/Toth-Kelley):

1. No independence test on the residual-difference columns
   `d_j = r_1 - r_{j+1}` entering the Gram (the Type-II
   analogue of Zhang's `||s_hat||` restart - NOT the same
   criterion; Zhang tests iterate-secants for a `B_k` update
   this code does not have).
2. No explicit `||alpha||_1` bound - Toth-Kelley's bounded-
   coefficient assumption is enforced only indirectly via the
   `1e3*||f||` veto.
3. Norm asymmetry: Gram fit in weighted-l2, acceptance in
   weighted-max. A mix can reduce the l2 objective while
   raising the max residual and be rejected. Stricter than
   needed but consistent and documented here; not a defect,
   a design property to keep visible.
4. History scope: per (ray, direction-field, refresh call),
   allocated/freed inside `SWEEP_FLD_ORD`. Automatic
   invalidation on coefficient generation (new C-pack => new
   arrays). CONFIRMED CORRECT SCOPE.
5. `AA_OFF` latch on rejection is for the rest of the solve;
   `NH_AA = 0` clears history; census rebuilt from slot-2 pair.
   Confirmed complete field-state rollback to `G_omega(x_k)`.

### 6.2 `MOM_AA_OUT` outer-map AA (CONFIRMED DEFECT C2)

`plane/mom_j_cmf_nonmono_v1.f` 6548-6750+:

- `XHIST_OUT/FHIST_OUT` are SAVE-persisted ACROSS refresh calls
  (depth `MOM_AA_OUT_M = 8`, default ON via `MOM_AA_OUT=T`).
- `FHIST_OUT` holds committed deltas `J_new - J_prev` between
  calls; per-cell weighting, `MOM_AA_GCAP` affine subset,
  contractive retry (`||corr|| > ||f||` => drop oldest secant,
  re-solve at m-1), step cap `10*max(|f|,floor)`, positivity
  guard, `MOM_AA_OUT_APP/REJ/SNG` counters.
- Operator-generation defect: secant pairs span calls in which
  populations - hence `chi, eta`, hence `G` - changed. The code
  rescales history on `JNEXT_RELAX` change (comment at ~6942
  correctly identifies an omega change as a map change) but
  does NOT invalidate on population commits - an operator
  change of the same kind and larger magnitude. Violates the
  same-map requirement (Walker-Ni equivalence; Zhang's
  secants-of-one-`f`; spec restriction). CONFIRMED DEFICIENCY.
- Secondary: stored `f_j` is the relaxed delta
  `omega*(G(x_j)-x_j)`; the omega-rescale fixes units across
  omega changes but the relaxed-vs-raw mixing assumption
  inside one history is unverified.
- Asymmetry note: the obsidian notes record the *deliberate
  deferral* of outer-map AA for the FG `I_NEXT` store
  ("needs a contractive inner map first") while the MOM
  analogue ships default-on. Inconsistent policy worth a
  default review.

### 6.3 Aitken / P2 damping

- `KN_ORD_AITKEN` (default F): guarded scalar Aitken inside the
  ordered march using `FLD_PRV`/`DLT_PRV`; deltas rescaled on
  omega change (8561-8562).
- `KN_ORD_P2DAMP` (default F): per-cell period-2 collapse -
  `P2_SGN` sign of last delta, `P2_FLIP` consecutive flips,
  `P2_DMP` damped state released after `FG_ORD_P2REL` flip-free
  passes. Targets `lambda ~ -1` modes where AA secants
  alternate sign and the Gram degenerates - the correct tool
  for that regime. Interaction with AA on the same cells is
  undefined (Audit S2).

## 7. JFNK side (verified faithful - see Review Sect. 6)

- EW96 Choice 2 implementation confirmed: formula, conditional
  `eta_sav > 0.1` safeguard, `eta_min = 3e-4`,
  `eta_max = 0.9`, ratio on accepted `||F||` (`jfnk_tol_prev`
  = `tol_` at 2038), `gres` in `||M^{-1}F||` units under left
  preconditioning.
- EW94 structure confirmed: true-residual acceptance at the
  trial point (`call f(n, jfnk_x, jfnk_fx)`, `ared`, `rho`),
  hookstep trust region, descent retry, NaN guard, full
  `jfnk_x/jfnk_fx/jfnk_tol` restore on rejection.
- Phase-6 items only: parameter audit; residual-mode
  consistency within a solve; `eta_prev` persistence across
  solves (harmless heuristic).

## 8. CRITICAL REVIEW OF LITERATURE-TO-CMFGEN ASSUMPTIONS

For every proposed numerical method: the theorem used, its
assumptions, whether CMFGEN satisfies them, whether the result
transfers, and what evidence is required.

| Method | Theorem/result | Assumptions | CMFGEN satisfies? | Transfers? | Status |
|---|---|---|---|---|---|
| Ordered GS march | Knop II Sect.3 ordering; standard GS-on-triangular-structure | Frequency coupling is one-directional per cell | YES by construction (`sign(GAM)` per row) | Direct - it IS the reference method | THEOREM/METHOD APPLIES |
| Under-relaxed map `G_omega` | KM iteration contraction `|1-omega(1-lambda)|<1` | Affine local map, `omega in (0,1)` | Locally, where affine | Local theorem + implemented mechanism | SAFEGUARD + LOCAL THEORY |
| Inner AA-II | Toth-Kelley Thm 2.6 (r-linear, factor c) | LOCAL contraction + bounded alpha + fixed map | Contraction NOT proven; bounded alpha only via veto; map fixed within solve | Conditional - applies where local affine regime holds | HEURISTIC WITH PARTIAL THEORY |
| AA-II residual acceptance | EW94-style sufficient decrease | True residual measure | YES - `RES` is raw undamped | Direct as safeguard; not a convergence theorem | SAFEGUARD |
| Zhang stabilization | Zhang et al. global convergence | Non-expansive `f` in l2 + Type-I `B_k` + all three mechanisms together | NO - Type-II code, no `B_k`, non-expansiveness unestablished | DOES NOT TRANSFER as theorem | ARCHITECTURE INSPIRED ONLY |
| Evans rate factor `theta_k*kappa` | Evans Thm 4.1 | GLOBAL `kappa<1`, `g in C^2`, alpha bounded AND bounded away from 0 | NO global kappa; C^2 fails at guards/caps | Conditional plausibility in local affine regime | CONDITIONAL HEURISTIC |
| `rho_k` residual ratio | empirical diagnostic | same ray, same map, stable argmax | YES when those hold | Level 1-2 evidence only, never a bound | DIAGNOSTIC |
| `KN_RES_PASS` envelope ratio | empirical | argmax-stable rays | Partially - envelope can switch rays | Level-1/2 at best | DIAGNOSTIC, argmax-gated needed |
| `MOM_AA_OUT` | none - violates same-map | one operator per history | NO - spans population updates | Violates Walker-Ni/Toth-Kelley premise | DEFECT C2 |
| `JNEXT_RELAX` estimator | contraction ratio | same metric numerator/denominator | NO - mixes sweep drift / refresh residual | Invalid as contraction | DEFECT C1 |
| EW96 Choice 2 | EW96 Thms | accepted-residual ratio, `eta_max<1` | YES (verified in source) | Direct | THEOREM APPLIES |
| EW94 globalization | EW94 | true-residual acceptance, clean rollback | YES (verified) | Direct | THEOREM APPLIES |
| Dogleg | Pawlowski | TR + `J^T F` | Hookstep TR already covers role | Redundant so far | DEFERRED |
| Multigrid | B&L17, F&B97 | coupled SE+RT map, structured mesh hierarchy | Level-1 scope wrong; no comoving-frequency variant | Does not transfer to Level-1; conceivable at Levels 2-3 | DEFERRED, RECLASSIFIED |

Hierarchy of claims (Issue G), assigned to code observables:

- Level 1 observed behaviour: single `RES` ratio; `KN_RES_PASS`
  envelope ratio; `RATIO_OM` between consecutive passes.
- Level 2 local linear regime: multi-pass stable ratios on same
  ray + stable argmax + both metrics; `cycle`/`margnl` classes.
- Level 3 local-contraction support: Level 2 + accepted AA mixes
  improving true residual + `||alpha||` bounded.
- Level 4 contraction bound: requires operator-level Lipschitz
  analysis of the march - NOT obtainable from runtime data.

## 9. Source-claim verification table (Issue F)

Every claim about the implementation in the review, verified at
commit `cfae5852`. "Confirmed" = traced in code; "inferred" =
consistent with code but not fully traced.

| # | Claim | File/Routine | Lines | Verdict |
|---|---|---|---|---|
| A | AA history is per-ray | `fg_j_cmf_nonmono_v1.f`/`SWEEP_FLD_ORD` | 7945-7975 decl, 8073-8084 alloc/push, 8829-8830 free | CONFIRMED - local to one `SWEEP_FLD_ORD` call = per (ray, direction-field, solve) |
| B | History does not persist across refreshes | same | 8073-8074, 8829-8830 | CONFIRMED - allocated/deallocated inside routine |
| C | History cleared on coefficient change | same | moot - never leaves the solve; `NH_AA=0` on reject (8446) | CONFIRMED - automatic by scope; explicit on rejection |
| D | Pending-mix acceptance uses true undamped residual | `SWEEP_FLD_ORD` | `RES` def 8168-8173 (`|VNEW-VOLD|`, `VNEW` raw at 8155-8158); test 8443-8445 | CONFIRMED - `RES`/`AA_REF` in raw-residual units; omega-independent |
| E | Fallback restores complete ray state | same | 8444-8446 restore `FLD = XH(2)+RH(2)`; census rebuilt 8455-8463 | CONFIRMED - restores `G_omega(x_k)` field + residual + census; deliberately clears AA auxiliaries |
| F | Same norm in Gram and acceptance? | same | Gram weighted-l2 8722-8731; acceptance weighted-max 8443 | CONFIRMED DIFFERENT - weighted-l2 fit vs weighted-max acceptance; documented as stricter-than-fit acceptance |
| G | Contraction diagnostic same ray/direction/coeffs? | `SWEEP_FLD_ORD` | `RATIO_OM = RES/RES_PREV` 8516; alternating `RES_PREV(IDIR)` 7592 | CONFIRMED for in-solve `RATIO_OM`; NOT for cross-call estimators (4.2, C3) |
| H | Adaptive relaxation = measured contraction or sign-flip? | same | 8539-8553 | CONFIRMED BOTH - measured `RATIO_OM` + `NOSC` flip trigger; bounds [0.45,1.0] |
| I | JFNK EW Choice 2 on accepted residual | `jfnk_solver_v1.f90` | 1019-1027 formula+safeguard; 897 `jfnk_tol`; 2038 `jfnk_tol_prev = tol_`; restores 1062 etc. | CONFIRMED faithful |
| J | JFNK acceptance = actual nonlinear residual + rollback | same | 1497 `call f`, `ared` 1502, `rho` 1515; NaN guard ~1615; restore 1062/1143/1212 | CONFIRMED |

Additional verifications requested by the audit:

| # | Claim | Lines | Verdict |
|---|---|---|---|
| K | `KN_RES_PASS` is per-pass max over rays, not per-ray history | 7558-7563 (`MAX` under OMP CRITICAL), 6755 comment | CONFIRMED - envelope, argmax-stability caveat applies |
| L | `FG_SW_*`/`FG_GS_*` metric split | 6589-6600 | CONFIRMED - churn snapshotted to `FG_SW_*`, `FG_GS_*` takes refresh residual; both exported |
| M | `MOM_AA_OUT` history SAVE-persisted across calls | `mom_j_cmf_nonmono_v1.f` 443-468, 6548-6560 alloc, 6598-6605 push | CONFIRMED - basis for defect C2 |
| N | `OMEGA_RAY` change rescales `AA_RH`/`DLT_PRV` | 8558-8562 | CONFIRMED - relaxed-units bookkeeping correct |
| O | `INEXT_RELAX` V2 = median of last ~4 envelope ratios | 6538-6578 | CONFIRMED - same argmax-envelope caveat as fallback |
| P | `JFNK_COUPLED` residual mode fixed per solve | `jfnk_solve_corr_v1.f` mode select; `JFNK_FRZ_HITS/CPL_HITS` hysteresis in `solve_for_pops.f` ~510-575 | INFERRED - admission machinery is per-call; mid-solve mode switch appears structurally impossible but full trace not re-done |
| Q | `sign(GAM)` per-row frequency ordering | 8135-8139 | CONFIRMED |
| R | OBS refresh shares `FG_KNOP_REFRESH_V1` | `obs/cmf_flux_sub_v5.f` ~2204-2216 | CONFIRMED - OBS-grid `SET_FG_GS_OBS_GRID -> RESTORE_FG_INEXT_OBS` then refresh per chunk; AA off by default so no cross-grid history issue today, but all metrics share the same module state - flagged S3 |

Unresolved uncertainties carried forward:

- U1. Whether the MOM inner-AA (even-pass boundaries, `MOM_AA_M=3`)
  has the same per-scope cleanliness as FG's (assumed by analogy;
  full trace not redone).
- U2. `DTR` symbol reuse between Gram ridge and contraction
  estimate in `MOM_KNOP_REFRESH_V1` - separate sites appear to
  use separate storage but was not fully traced.
- U3. Whether `KN_RES_PASS` envelope argmax stability is checked
  anywhere before the V2 omega update - found none (deficiency
  C3 stands).

## 10. Confirmed vs suspected inconsistencies

CONFIRMED:

- C1. `JNEXT_RELAX` estimator mixes sweep drift with refresh
  residual (Sect. 4.3; `mom_j_cmf_nonmono_v1.f` ~6931).
- C2. `MOM_AA_OUT` spans operator changes with only an omega
  rescale (Sect. 6.2); default-on contradicts the FG deferral
  policy for the same method class.
- C3. Cross-call omega estimators use global-max envelopes
  without argmax-stability gating (Sect. 4.2).

RECLASSIFIED after critical review:

- "Missing Zhang restart criterion" -> corrected to "missing
  Type-II column-independence test" (5.1 in TODO); the
  Zhang `||s_hat||` rule is Type-I-only and cannot be cited
  for this Gram system.
- "Powell regularization alternative" -> REMOVED: Powell's
  formula is specific to Type-I rank-one `B_k` updates and has
  no Type-II counterpart here.
- "Pending-mix is the Zhang safeguard" -> softened to
  "EW94-style residual-acceptance safeguard, architecture
  inspired by Zhang"; the global theorem does not transfer.
- "Multigrid not recommended (hyperbolic)" -> reclassified:
  established coupled-problem technique, wrong level for
  Level-1, deferred accurately.

SUSPECTED / TO VERIFY: S1 (`DTR` reuse), S2 (AA x P2
interaction), S3 (OBS-grid shared module state), S4
(`GS_HOLD_POPS` completeness on all pop-mutating paths).

## Section A - Nov-2026 audit-cycle findings

### A.1 Inner versus outer relaxation operators

The pooled residual ratios `RHO_M_LS`/`RHO_P_LS` measure the
per-ray, same-direction contraction of the INNER fixed-point map
`I* = G(I; chi, eta)` inside one refresh, sampled post-transient.
They do not estimate the spectral radius of the OUTER lagged-store
update (sweep write plus refresh commit), which additionally moves
when populations and coefficients change between great iterations.
Using the inner ratio to damp the outer update is an explicitly
empirical heuristic: the estimator now only runs when it has
enough same-map samples (`KN_RHO_MIN_SAM`) and is documented as
such rather than as an optimal damping law.

### A.2 Residual ratios versus spectral-radius estimates

For a local linearised fixed point with eigenvalue `lambda`, the
relaxed iterate has eigenvalue `mu = 1 - omega + omega*lambda`.
An unsigned residual ratio gives `|lambda|` at best and says
nothing about its sign: for `lambda` in `(0,1)`, under-relaxation
is slower than `omega=1`; for `lambda<0` it suppresses oscillation.
The formula `omega = 1/(1+rho)` is therefore retained strictly as
a bounded heuristic in `[0.05, 1.0]`, never advertised as the
optimal relaxation. A signed store-drift history would be needed
for anything stronger; that remains a documented limitation.

### A.3 Retained-store maturity after rejected commits

Inner-solve convergence and store retention are distinct states.
A ray can converge internally and still be vetoed (jump guard) or
rejected (non-converged class), leaving the retained field at its
previous generation. The maturity signal must therefore include
the pending proposed-vs-retained delta of rejected rays
(`FG_RF_PEND` folded into `RF_DCHG`), not only the committed
delta; otherwise the population solver releases against a field
the solve never committed.

### A.4 DTAU clamp parity

The refresh `BUILD_FREQ` operator mirrors the
`SOLUTION_METHOD='INTEGRAL'` sweep block, which applies the
min/max DTAU clamp unconditionally: `CHI_RAY`/`CHI_L` carry the
`|GAM|/DLN` frequency-derivative coupling term, which monotonic
CHI interpolation cannot regularise. The earlier refresh-side
`INSERT`/`METHOD(4:6)='MON'` skip reproduced the DIFFERENCE-path
rule on the wrong operator and is removed; under `METHOD='ZERO'`
the code path is unchanged.

## Section B - fail-safe semantics (Nov-2026, e6e8f47a)

Status-carrying routines must never encode failure in data values
the caller also uses semantically. REF=-1 meant "no sane sibling"
(cold start) and "allocation failed" simultaneously; a failed
reference then un-disabled the floor-column veto or, in the
refresh-return variant, let the maturity gate read small sweep
metrics against an unvalidated field. The corrected contract:
ISTAT carries the outcome, REF carries only physics, and an
aborted refresh sets FG_RF_FAILED which the maturity export maps
onto an immature RF_DCHG. The same pattern closes the B1 operator
question: KN_DTAU_INT_V1 makes the sweep/refresh DTAU parity a
construction fact instead of two audited-but-separate formulas.
