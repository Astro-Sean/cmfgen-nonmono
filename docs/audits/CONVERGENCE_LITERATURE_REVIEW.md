# Convergence Literature Review - Non-Monotonic CMFGEN

Scope: which fixed-point and inexact-Newton methods are
mathematically appropriate for the non-monotonic CMFGEN problem,
and what assumptions each method needs. REVISED after a critical
audit of the first version: Type-I/Type-II Anderson distinction
corrected, the Evans et al. convergence-rate formula corrected,
the Zhang-Boyd theorem's assumptions stated exactly, the multigrid
assessment qualified, and every source claim re-verified at commit
`cfae58528648ed138749f8748e195226ab03b0f3` (branch `main`).

The inner problem under study is

```
I* = G(I*; chi, eta)
```

where `I` is the stored comoving-frame intensity field (`I_NEXT`),
`G` is the formal transport map (short-characteristics solve with
monotone-cubic chi/eta coefficients built from the captured sweep
fields), and `chi`, `eta` are frozen for one call of
`FG_KNOP_REFRESH_V1`. Population updates change `chi`, `eta` and
therefore change `G` between refreshes.

Terminology discipline used throughout:

- "Theorem applies" means the stated assumptions are verified or
  reasonably establishable for the CMFGEN operator.
- "Safeguard inspired by" means the code borrows a mechanism whose
  published proof does not cover our operator.
- "Heuristic" means engineering judgment with no theorem behind it.

## 1. Toth and Kelley (2015) - convergence analysis for AA

Ref: Toth, A. and Kelley, C. T., "Convergence Analysis for Anderson
Acceleration", SIAM J. Numer. Anal. 53(2), 805-819, DOI
10.1137/130919398.

Sections consulted: 1-4 (algorithm statement, linear case,
nonlinear/local theory, conditioning discussion).

Results that matter here:

- Anderson(m) stores m+1 iterates; Anderson(0) is fixed-point
  iteration. The new iterate is a linear combination
  `sum_j alpha_j G(x_j)` with `sum alpha_j = 1`.
- Theorem 2.6 (main result): if `G` is a contraction with
  constant `c` and the Anderson coefficients stay bounded, then
  Anderson(m) is locally r-linearly convergent with factor `c`.
  This is a "do no worse than the fixed-point map" bound; it is
  not by itself a proof of acceleration.
- For `G` linear, `G(u) = Mu + b`, `||M|| <= c < 1`: the residual
  satisfies `||F(u_{k+1})|| <= c ||F(u_k)||` in `l2`; the proof
  uses only `sum alpha_j = 1` and the least-squares optimality of
  the residual combination. No coefficient bound is needed in the
  linear case.
- Anderson(1) converges q-linearly without an a-priori
  coefficient bound (the scalar optimization cannot blow up near
  the fixed point). For `m >= 2` nonlinear, bounded `alpha` is a
  genuine assumption: unbounded `alpha` is an unstable
  extrapolation on the linearized residual manifold.
- The least-squares conditioning of the residual-difference
  matrix `F_k` controls the size of `alpha`.

Application to CMFGEN:

- The ray-local map `G_ray` is not known to be contractive. It is
  piecewise defined (positivity guards at line 8159-8162, per-cell
  `NONMONO_INEXT_CAP` at 8163-8166, `sign(GAM)` direction switches
  at 8135-8139 of `plane/fg_j_cmf_nonmono_v1.f`). The theory can
  apply at most locally: fixed ray, fixed direction field, fixed
  coefficient pack, in a neighborhood of `I*` where the map is
  affine. Whether that neighborhood contains the iterate path is
  a runtime question, not a theorem.
- Under-relaxation can create local contraction where the raw map
  has a `lambda ~ -1` eigenmode: for `G_omega(x) = x +
  omega*(G(x)-x)` the affine factor is `I - omega(I - M)`, whose
  eigenvalues are `1 - omega(1 - lambda_i)`. A `lambda ~ -1` mode
  gives `|1 - 2*omega| < 1` for `omega in (0,1)` - this is the
  mechanism the code's `OMEGA_RAY` floor of 0.45 exploits
  (verified in source comment near line 8545).
- Requirements reduce to (a) a locally contractive regime and
  (b) bounded `alpha`. Both are runtime observables: (a) through
  consecutive same-ray residual ratios at fixed coefficients and
  fixed omega, (b) through `||alpha||` from the Gram solve.
- History depth: the theory gives no benefit for large `m`;
  boundedness of `alpha` is easier to keep at small `m`. The
  default `KN_ORD_AA_M = 2` is consistent with the theory.

## 2. Zhang, O'Donoghue and Boyd (2020) - globally convergent AA-I

Ref: Zhang, J., O'Donoghue, B. and Boyd, S., "Globally Convergent
Type-I Anderson Acceleration for Non-Smooth Fixed-Point
Iterations", arXiv:1808.03971, SIAM J. Optim. 30(4), 3170-3197,
DOI 10.1137/18M1232772.

Sections consulted: 1 (assumptions), 2 (Type-I vs Type-II), 3
(stabilization), 4 (convergence proof), 6 (extensions).

What the paper actually proves:

- Assumptions: `f` non-expansive in `l2`
  (`||f(x)-f(y)||_2 <= ||x-y||_2` for all x,y) and a nonempty
  fixed-point set. Non-expansive is weaker than contractive but
  is still a global Lipschitz bound.
- The algorithm has three coupled mechanisms, not independent
  tricks:
  1. Type-I multi-secant update: `B_k` approximates the Jacobian
     via rank-one secant updates
     `B_{k+1} = B_k + (y - B_k s) s_hat^T/(s_hat^T s_hat)` where
     `s_hat` is the component of the newest iterate secant
     orthogonal to the stored span.
  2. Powell-type regularization: replace `y` with
     `y~ = theta*y + (1-theta)*B*s`, `theta = phi_theta-bar(eta)`,
     giving `|det B_k| >= theta-bar^{m_k} > 0` (their Lemma 2).
     This regularization exists only because Type-I updates `B_k`
     through that rank-one formula.
  3. Restart checking (their rule 14): reset `m_k = 0` when
     `||s_hat||_2 < tau*||s||_2`, enforcing strong linear
     independence of the *iterate-secant* vectors `s` and bounding
     `||B_k||` (their Lemma 3). This is a condition on the secant
     `s = x_{k+1} - x_k`, not on the residual differences.
  4. Safeguarding steps: interleave AA steps with plain
     Krasnoselskii-Mann steps so residual decrease is guarded by
     the base iteration.
- Global convergence (`lim ||g_k|| -> 0`) is proved for the
  complete package under the non-expansiveness assumption.

Why this does not transfer as a theorem to CMFGEN:

- `G_ray` is not shown to be non-expansive in any norm. The
  short-characteristics march with lagged C-corner coupling,
  positivity clamps, per-cell caps, and direction switching can
  amplify perturbations transiently; non-expansiveness has not
  been established and probably fails at cap/guard boundaries.
- The CMFGEN implementation is Type-II (see Sect. A below), not
  Type-I. Zhang's Lemma 2/3 machinery is defined for the Type-I
  `B_k` update and has no Type-II counterpart in that paper.
- Therefore NOTHING in the current code "inherits" the Zhang
  global-convergence theorem. What is legitimately borrowed is
  the *architecture*: interleaved safeguarded iteration with
  regularization, independence control, and residual-acceptance.
  Those components in our Type-II code are safeguards inspired
  by the paper, not instances of the proved algorithm.

## 3. Issue A: Type-I vs Type-II in the CMFGEN source

Exact implementation at commit cfae5852,
`plane/fg_j_cmf_nonmono_v1.f`, routine `SWEEP_FLD_ORD` (contained
in `KNOP_RAY_SOLVE_V1`):

Stored history (lines 8072-8084, 8431-8432):
- `AA_XH(:,:,1..M+1)`: iterates `x_i` (the pre-march field of
  each pass), newest-first.
- `AA_RH(:,:,1..M+1)`: relaxed-map residuals
  `r_i = FLD_after_march - AA_XH = G_omega(x_i) - x_i`
  `= omega*(G(x_i) - x_i)`, newest-first. `G_omega` is the
  under-relaxed map (`FLD = VOLD + OMEGA_RAY*(VNEW-VOLD)`,
  line ~8422); `RES`/`AA_REF`/`RATIO_OM` use the *raw* residual
  `|VNEW - VOLD|/max(VOLD, floor)` (lines 8168-8173), so
  acceptance metrics are undamped.
- Allocation/deallocation inside `SWEEP_FLD_ORD` (8073-8074,
  8829-8830): history is local to one (ray, direction-field,
  refresh call). Never crosses rays, directions, refreshes, or
  coefficient generations. `NH_AA` cleared on rejection.

Least-squares problem (lines 8703-8752), `MP_AA = NH_AA - 1`:
- Per-cell weight `WV = 1/MAX(|x_1|, 1e-4*col_peak)` - the same
  normalized family as the march residual.
- `GRM(KK,JJ) = <W(r_1 - r_{KK+1}), W(r_1 - r_{JJ+1})>`,
  `BRH(JJ) = <W(r_1 - r_{JJ+1}), W*r_1>` - a Gram system for
  `min_gamma ||W*(r_1 - D*gamma)||^2` where `D` is the matrix of
  weighted residual *differences* `W*(r_1 - r_{j+1})`.
- Excluded cells: `max(|r_1|, max_j |r_1 - r_{j+1}|)*WV > AA_GCAP
  ( = 1.0)` - cells outside the affine regime (capped/clipped)
  are left out of the fit and keep the raw map output.
- Tikhonov ridge `+1e-8*max_diag(GRM)` on the diagonal
  (8737-8744); solve by `DGESV` (8751).

Update (lines 8758-8823):
- `x_mix = G_omega(x_k) - sum_j gamma_j*(G_omega(x_k) -
  G_omega(x_{k-j}))` (since `FLD - AA_XH(j+1) - AA_RH(j+1) =
  y_k - y_{k-j}`, the code's `XV` form is the map-output
  combination). With `alpha_0 = 1 - sum gamma_j` this is
  `x_mix = sum_j alpha_j*G_omega(x_{k-j})`.
- Vetoes/caps: `||corr||_W > 1e3*||r||_W` veto (8790-8793),
  per-cell step cap `10*|r_1|` (8807-8810), positivity fallback
  to `T4` (8811), `NONMONO_INEXT_CAP` ceiling (8812).

Acceptance (lines 8443-8463):
- `AA_REF = RES` stored at mix time. Next pass's `RES` (raw,
  undamped, weighted-max residual) must satisfy `RES <= AA_REF`,
  else restore `FLD = AA_XH(2) + AA_RH(2) = G_omega(x_k)` (the
  safe unaccelerated map output), set `RES = AA_REF`, clear
  `NH_AA`, disable AA for the rest of the solve (`AA_OFF`),
  `KN_AA_REJ++`, and rebuild the bad-cell census from the
  restored pair.

Classification: this is **Type-II Anderson acceleration** (the
Walker-Ni / Fang-Saad unconstrained reformulation: least squares
on residual differences, coefficients applied to map outputs).
It is emphatically not Zhang's Type-I: there is no `B_k`, no
rank-one update, no `s_hat` division anywhere in the code.

Consequences for the previous review's recommendations:

A. Zhang's restart criterion `||s_hat|| < tau||s||` cannot be
   applied directly: it tests independence of *iterate* secants
   `s`, which is what degenerates the Type-I `B_k` update
   (division by `s_hat^T s_hat`). Our Type-II Gram is built on
   *residual-difference* columns `d_j = r_1 - r_{j+1}`; its
   conditioning is controlled by the linear independence of
   those columns. The mathematically correct analogue is:
   orthogonalize the newest residual-difference column against
   the stored `D` span and require `||d_hat|| >= tau*||d||`
   (equivalently monitor the Gram's smallest pivot/rcond before
   solving). Same *role*, different matrix, not the published
   criterion - do not cite it as Zhang's rule 14.
B. For Type-II, the published-consistent independence controls
   are (i) Gram conditioning (smallest pivot of `GRM` after
   ridge, or LAPACK `rcond`), (ii) the `d_hat` test above, and
   (iii) an explicit `||alpha||_1` bound (Toth-Kelley's actual
   assumption). The current code has only the Tikhonov ridge
   and the `1e3*||f||` veto.
C. Sufficiency of current safeguards: the `1e-8` ridge prevents
   exact singularity but allows `||gamma||` up to ~1/ridge
   scales in the worst direction; the veto catches gross
   failures and the 10|r| cap limits per-cell damage. These are
   incomplete against *moderately* ill-conditioned histories:
   a near-dependent Gram can produce a large-but-sub-veto mix
   that fails the next-pass residual test, costing one wasted
   pass and an AA_OFF latch. A column-independence test closes
   this hole at negligible cost (MP_AA <= 8).
D. Switch to Type-I: no genuine advantage identified. Type-I
   buys Powell-regularized conditioning by construction and a
   cheaper update (no Gram solve), but its published theory
   additionally needs non-expansiveness (not established), the
   `H_k` history would still need operator-change invalidation
   on every `OMEGA_RAY` adaptation (the code already rescales
   `AA_RH` on omega change, lines 8558-8562 - a Type-I `H` has
   no such simple fix), and Type-II is the form the
   Toth-Kelley/Walker-Ni theory actually covers. Keep Type-II;
   add the correct independence control.

## 4. Evans, Pollock, Rebholz and Xiao (2020) - rate improvement

Ref: Evans, C., Pollock, S., Rebholz, L. and Xiao, M., "A Proof
That Anderson Acceleration Improves the Convergence Rate in
Linearly Converging Fixed-Point Methods (But Not in Those
Converging Quadratically)", SIAM J. Numer. Anal. 58(1), 788-810,
DOI 10.1137/19M1245384, arXiv:1810.08455.

CORRECTED EQUATIONS (previous version's
`sqrt(rho^2 - orthogonality gain)` was wrong and is removed):

- Eq 3.27: the inner optimization is
  `min_gamma ||w_{k+1} - F^k gamma^k||^2` where `F^k` has columns
  `w_{n+1} - w_n` for `n = k-m+1..k` (residual differences).
- Eq 3.28: writing the optimum as
  `||w_{k+1} - F^k delta^k|| = ||w^alpha_{k+1}|| =
  theta_k*||w_{k+1}||`, the optimization gain is
  `theta_k = sqrt(1 - ||w_R||^2/||w_{k+1}||^2)` where `w_R` is the
  projection of `w_{k+1}` onto `span{w_{n+1}-w_n}`. `theta_k` is
  the direction-sine between the current residual and the
  residual-difference subspace: how much of the residual lies
  OUTSIDE what the history can represent. `theta_k = 0` would be
  perfect acceleration but also means `alpha` unbounded.
- Assumption 3.1: `g in C^2(X)`, fixed point exists,
  `||g'(y;u)|| <= kappa||u||` for all y,u and
  `||g''(y;u,v)|| <= khat||u||||v||`.
- Assumption 3.2: `||g(y)-g(x)|| <= kappa||x-y||` for every x,y
  with `kappa < 1` - a GLOBAL contraction bound.
- Theorem 4.1 (m=1): under Assumptions 3.1-3.2, and provided the
  `alpha` coefficients are bounded AND bounded away from zero,
  `||w_{k+1}|| <= theta_k*((1-beta_{k-1}) + kappa*beta_{k-1})
  ||w_k|| + O(higher-order terms)` where `beta` is damping
  (`beta=1` undamped). For the undamped iteration the leading
  rate is `theta_k*kappa`.
- Prop. 3.3/Theorem for m=2: same structure with an extended
  history; higher-order terms grow with depth.
- The `O(...)` terms are quadratic-and-higher in the residuals
  (from the `g''` Taylor remainder). They matter near
  convergence and are the stated reason AA does not improve
  quadratically convergent iterations: there the linear term is
  already faster than `theta*kappa`, and the higher-order terms
  can dominate the bound.

Meaning of theta_k, carefully:

1. theta_k measures the fraction of `||w_{k+1}||` orthogonal to
   the span of residual differences - a per-step geometric
   quantity, not a property of `g` alone.
2. It depends on the history subspace: more secant directions
   can reduce theta_k, but near-dependent columns make `alpha`
   diverge, violating the theorem's own boundedness assumption.
3. The leading-order bound multiplies the unaccelerated rate by
   theta_k in [0,1].
4. Higher-order terms are what break the argument for
   quadratic/superlinear iterations and what can dominate when
   residuals are not yet small.
5. It does NOT establish universal acceleration: it requires
   kappa < 1 globally (Assumption 3.2), `g in C^2`, and bounded
   `alpha` bounded away from zero - the last being exactly what
   fails near rank-deficient history.

Application to CMFGEN:

- Assumption 3.2 (global `kappa < 1`) is not established for
  `G_ray`; the correct claim is conditional: IF the iterate
  sequence is in a local affine regime with dominant eigenvalue
  `|lambda| < 1`, THEN AA plausibly improves the observed rate
  by the direction-sine factor on those modes.
- A measured ratio `rho_k = ||r_k||/||r_{k-1}||` in [0.9,0.99]
  over several consecutive passes, on the same ray, same
  coefficient pack, same omega, same argmax-cell set, and both
  max and bulk metrics, is *evidence of an observed slow
  linearly-converging regime*. It is NOT a proof that
  `||G(u)-G(v)|| <= rho*||u-v||` for all u,v - transients, mode
  switching and the non-affine guards can all produce such a
  ratio on a non-contractive map. Single ratios and
  different-argmax ratios are weaker still. See Sect. 8 for the
  evidence hierarchy.
- Practical consequence unchanged: enable AA only where the
  observed regime is linear-contractive; that regime must be
  established from residual histories (Phase-1/2 diagnostics),
  not assumed.

## 5. Walker and Ni (2011) - AA/GMRES equivalence

Ref: Walker, H. F. and Ni, P., "Anderson Acceleration for
Fixed-Point Iterations", SIAM J. Numer. Anal. 49(4), 1715-1735,
DOI 10.1137/10078356X.

Results that matter here:

- For a linear fixed-point map `f(u) = Mu + b`, AA with
  unlimited depth is mathematically equivalent to GMRES applied
  to `(I-M)x = b`: the accelerated iterates minimize the
  residual over the same Krylov space. Bounded `m` corresponds
  to restarted GMRES.
- AA is a multi-secant method; Type-II updates the inverse
  Jacobian estimate `H_k` directly (satisfying `H_k y_i = s_i`
  on stored pairs); Type-I updates `B_k` and inverts.
- The equivalence is exact only for a single global linear map
  and a consistent inner product over the full unknown vector.

Application to CMFGEN:

- With frozen coefficients, `G_ray` is affine in `I` except at
  guard/cap boundaries (the short-characteristics coefficients
  are built once per refresh from captured fields). Locally on
  one ray the map is therefore affine to leading order, and the
  per-ray Type-II AA is, to first order, a low-memory GMRES on
  that ray's `(NI x NCF)` field.
- The ordered sweep is NOT a global linear iteration across
  rays or directions: `sign(GAM)` switching means the effective
  operator differs by cell regime. A global AA history over the
  whole `I_NEXT` store would mix secants of different local
  operators - the Walker-Ni equivalence fails there. The
  existing per-(ray, direction-field, call) scope is the
  mathematically defensible one; do not widen it.
- No per-direction split is needed within a sweep call: each
  `SWEEP_FLD_ORD` call already owns exactly one direction
  field (`I_M` or `I_P`), so "per direction" is already the
  implemented granularity (lines 7126, 7181).

## 6. Eisenstat and Walker - inexact Newton (papers verified vs source)

Ref A: Eisenstat, S. C. and Walker, H. F., "Choosing the Forcing
Terms in an Inexact Newton Method", SIAM J. Sci. Comput. 17(1),
16-32, DOI 10.1137/0917003.
Ref B: Eisenstat, S. C. and Walker, H. F., "Globally Convergent
Inexact Newton Methods", SIAM J. Optim. 4(2), 393-422, DOI
10.1137/0804022.

EW96 results used:

- Inexact Newton solves `J s = -F` to
  `||F + J s|| <= eta_k ||F||`. Choice 2:
  `eta_k = gamma*(||F_k||/||F_{k-1}||)^alpha`, `alpha in (1,2]`;
  safeguard `eta_k >= gamma*eta_{k-1}^alpha` WHENEVER
  `gamma*eta_{k-1}^alpha > 0.1` (conditional in the paper);
  bounds `eta in [eta_min, eta_max]`, `eta_max < 1` needed by
  EW94.

EW94 results used:

- Global convergence needs `eta_k <= eta_max < 1`, a
  sufficient-decrease acceptance test on the true nonlinear
  residual, and a fallback (line search/trust region) when the
  direction is unproductive. Rejected steps must leave the
  nonlinear state untouched.

Source verification (`new_main/mod_subs/jfnk_solver_v1.f90`):

- CONFIRMED faithful Choice 2: line 1019-1022 computes
  `jfnk_eta = gamma*(jfnk_tol/jfnk_tol_prev)^alpha` with
  `gamma = 0.9`, `alpha = 1.618034` (in (1,2]), and applies the
  EW conditional safeguard `eta_sav = gamma*jfnk_eta_prev^alpha;
  if(eta_sav > 0.1) jfnk_eta = max(jfnk_eta, eta_sav)` - exactly
  the published condition, not an unconditional max.
- CONFIRMED the ratio is on the accepted outer residual:
  `jfnk_tol = ||F(x)||` physical norm (line 897,
  `dsqrt(dp(n,jfnk_fx,jfnk_fx))`); `jfnk_tol_prev = tol_` at the
  end of each Newton step (line 2038); rejected steps restore
  `jfnk_x = x_; jfnk_fx = fx_; jfnk_tol = tol_` (e.g. 1062),
  so the ratio runs over accepted iterates only.
- CONFIRMED bounds: `eta_min = 3.0e-4`, `eta_max = 0.9`
  (lines 123-124), first step uses `gtol` clamped to the same
  interval (1007).
- Units note: the *forcing ratio* uses `||F||` (paper form); the
  *GMRES tolerance* `gres = eta*ref_norm` uses
  `ref_norm = ||M^{-1}F||` under left preconditioning
  (1040-1069) so that the inexact-Newton condition is measured
  in the same units GMRES actually tests. This is consistent
  with EW's definition - their `eta` multiplies `||F||` in the
  norm in which convergence is measured; the code documents
  this as a deliberate unit fix.
- CONFIRMED acceptance on the actual nonlinear residual:
  `call f(n, jfnk_x, jfnk_fx)` at the trial point (line ~1497),
  `ared = tol_ - jfnk_tol`, `rho = ared/pred`, NaN guard, and
  full state restore on rejection. Trust-region `jfnk_del`
  maintained in the constrained norm; descent-direction retry
  tightens `eta` (`jfnk_eta*0.1`, line ~1449) as a rescue, with
  the measured true residual short-circuiting hopeless retries.
- `JFNK_FADAPT = F` restores the fixed `eta` (line 1007 branch).
- Caveat recorded for the audit: `jfnk_eta_prev` persists
  across separate Newton solves; the EW96 safeguard is defined
  between consecutive steps of one solve - reusing it across
  solves is a harmless heuristic given the `[eta_min,eta_max]`
  clamp, not a theorem violation.
- Residual-mode caveat for Phase 6: the EW ratio is meaningful
  only if successive Newton iterates evaluate the SAME `F`.
  `JFNK_COUPLED_THIS` fixes the residual mode per solve
  (frozen-J vs coupled) so the ratio is consistent within a
  solve; confirm the admission machinery can never switch modes
  mid-solve (hysteresis counters `JFNK_FRZ_HITS/CPL_HITS`).

## 7. Pawlowski et al. (2006) - inexact Newton dogleg

Ref: Pawlowski, R. P., Simonis, J. P., Walker, H. F. and
Shadid, J. N., "Inexact Newton Dogleg Methods", SIAM J. Numer.
Anal. 46(4), 2112-2132, DOI 10.1137/050632166.

- Dogleg replaces a failed/out-of-region Krylov step by a
  Cauchy-point/dogleg combination inside the trust region - a
  cheaper rescue than a second GMRES solve, at the price of a
  `J^T F` product not currently in the residual interface.
- The existing hookstep trust region plus descent retry covers
  the same failure mode. Recommendation stands: defer unless
  the Phase-6 acceptance audit finds a specific gap.

## 8. Issue G: residual ratios vs contraction - evidence hierarchy

For a fixed map `r_k = G(I_k) - I_k`, the ratio
`rho_k = ||r_k||/||r_{k-1}||` is an empirical convergence
diagnostic. It is not a contraction bound. Legitimate claim
levels:

- Level 1 - observed behaviour: single ratio in a plausible
  band. Says only "this pass shrank the residual".
- Level 2 - local linear regime: several consecutive ratios
  near-constant in (0,1), same ray, same coefficients, same
  omega, stable argmax cell, max and bulk metrics agreeing.
  Evidence of observed linear convergence - the Evans
  precondition.
- Level 3 - local-contraction support: Level 2 plus
  direction-consistent behaviour (same sign of update), plus
  AA acceptance-tested gains consistent with `theta_k < 1`.
  Supports the hypothesis of a locally contractive affine
  regime; still not a bound.
- Level 4 - contraction proof: `||G(u)-G(v)|| <= rho*||u-v||`
  for all u,v in the domain (or a neighbourhood). Requires an
  operator-level Lipschitz analysis of the short-characteristics
  map (feasible in principle: the march is explicit in the
  upwind/source terms, so the amplification factors are
  inspectable); not obtainable from runtime residual data.

What each instrument legitimately measures:

- `rho_k` (consecutive residuals): Levels 1-2.
- `KN_RES_PASS` per-pass max envelope: a Level-1/2 global
  envelope; a contraction estimate only if the argmax cell is
  stable across the compared passes (the `FG_ARG_*`/`cycle`
  classification exists for this).
- `RATIO_OM = RES/RES_PREV` inside `SWEEP_FLD_ORD`: same ray,
  same pass-sequence metric - a Level-2 ray-local diagnostic
  (ordered sweep uses consecutive passes; the legacy
  alternating `SWEEP_FLD` keeps per-direction `RES_PREV(IDIR)`).
- `theta_k` (AA gain): Level 3 evidence when accepted mixes
  keep improving the measured residual.
- `||alpha||_1`, Gram conditioning: boundedness diagnostics for
  the Toth-Kelley assumption; large values are evidence AGAINST
  the local affine hypothesis.
- Cross-call ratios (`KN_RF_RESID/KN_RF_RESID_PRV`): Level-1
  only - different coefficient generations make them
  `||G_{k+1}(x)-x||/||G_k(x)-x||` over different maps.

## 9. Saad (2025) - review

Ref: Saad, Y., "Acceleration methods for fixed-point
iterations", Acta Numerica 34, DOI 10.1017/S0962492924000096.

- Places AA in the multisecant/low-rank extrapolation family
  with Ng/Aitken, MPE/RRE and reduced-rank methods; the shared
  requirement is a locally linear contractive regime.
- The consistent warning: extrapolation across an operator
  change mixes secants of different maps - the formal statement
  behind our "no history across coefficient generations" rule
  and behind the documented `NG_CLEAN`/polluted-history
  failures of the population-level Ng accelerator.

## 10. Issues D and E: radiative-transfer context

Bergemann and Hoppe (2026), "3D NLTE radiation transfer",
Living Rev. Comput. Astrophys. 12, art. 7, DOI
10.1007/s41115-026-00029-3, plus the primary multigrid refs:
Bjørgen & Leenaarts 2017, A&A 599, A118, DOI
10.1051/0004-6361/201630237; Fabiani Bendicho, Trujillo Bueno
& Auer 1997, A&A 324, 161-176.

### 10.1 What the refresh actually is (Issue E)

Definitions (Bergemann-Hoppe Sects. 6.2-6.3, standard
Olson-Auer-Buchler framing):

- Classical Lambda iteration: iterate the COUPLED map
  `J -> SE -> S -> Lambda -> J`: populations and field both
  move.
- ALI: split `Lambda = Lambda* + (Lambda - Lambda*)` and
  precondition the coupled update with `Lambda*` - still a
  coupled SE+RT iteration.
- `FG_KNOP_REFRESH_V1`: populations are FROZEN, hence
  `chi, eta`, the source function `S = eta/chi`, and all
  short-characteristics coefficients are fixed for the whole
  refresh. The iteration is on the INTENSITY FIELD ALONE:

  ```
  I_{k+1} = G_c(I_k)   (frozen chi, eta)
  ```

  where the only self-coupling is through the lagged source
  terms at cell corners (the C-corner / source-interpolation
  part of the discretization reads stored field values).
- Therefore the refresh is NONE of: a Lambda iteration (no
  populations in the loop), ALI (no Lambda* splitting of the
  coupled operator), or operator splitting in the ALI sense.
- The precise label: a SUCCESSIVE-SUBSTITUTION fixed-point
  iteration on the frozen-coefficient formal-solution map,
  marched as a physically ordered Gauss-Seidel-type sweep
  (cells updated in causal order within each pass, using the
  newest available values). The "lag" is in the discretization
  itself, not an operator split.
- The earlier description "ALI-like lagged-field iteration"
  is retained only in the loose sense that a stored field
  plays the role of the previous iterate; it must not be read
  as an ALI method.

### 10.2 The sign(GAM) ordering (Issue E, primary source)

- Code (verified): `plane/fg_j_cmf_nonmono_v1.f` lines
  8135-8139: per spatial row `J`, `GAM(J,LS) >= 0` marches
  frequency `L = 1..NCF` ascending; `GAM < 0` marches
  descending. Direction is chosen per depth cell by the LOCAL
  sign of the frequency-coupling coefficient.
- Physics (Knop et al. 2009b, Paper II, A&A 501, 813, Sect. 3;
  also the strategy audit `cmfgen_obsidian/audits/
  KNOP_CONVERGENCE_STRATEGY.md`): the CMF transfer equation
  contains a `gamma*dI/dnu` advection term with `gamma`
  proportional to the velocity gradient. In expansion
  (`gamma > 0` in this sign convention) comoving photons drift
  toward lower frequencies; in compression (`gamma < 0`) the
  drift reverses. Marching in the advection direction makes
  each frequency's update read only already-updated upwind
  neighbours - a causally ordered Gauss-Seidel update of the
  frequency-advection operator. Note: the local PDF copy was
  not re-extracted in this audit (no text tool); section-level
  citation only, no equation numbers claimed.
- Non-monotonic velocity changes WHERE `GAM` flips sign along
  the ray; the per-cell `LB/LE/LDI` logic handles that - this
  is exactly the machinery that makes the ordered march
  non-monotonic-capable and also the reason a single global
  linear-operator picture fails.

### 10.3 Multigrid (Issue D) - qualified assessment

What the cited methods actually do:

- Fabiani Bendicho, Trujillo Bueno & Auer (1997): NONLINEAR
  (FAS) multigrid on the COUPLED multidimensional RT +
  multilevel-atom statistical-equilibrium system. It
  accelerates the outer SE+RT coupled problem, with the local
  RT iteration as smoother.
- Bjørgen & Leenaarts (2017): nonlinear FAS multigrid in
  Multi3D where the MALI iteration is the smoother; speed-ups
  3.3-6x over MALI on 3D Cartesian atmospheres. Again the
  COUPLED rate-equation problem `W(n) n = f`.
- Steiner (1991, cited in B&L): linear multigrid on two-level
  coherent-scattering - again the coupled map.

So: multigrid in RT is established, but its published targets
are the COUPLED SE+RT system (Levels 2-3 of this audit's
level scheme), NOT the frozen-coefficient formal-solve fixed
point (Level 1). The mechanism it removes is the smooth
low-spatial-frequency error that Jacobi/GS-type smoothers
cannot - the correct error-mode identification matters.

Applicability to CMFGEN:

- For the Level-1 refresh: the dominant slow modes are not only
  spatial low-frequency but (a) frequency-direction modes set
  by `sign(GAM)` ordering and (b) `|lambda| -> 1` and
  `lambda ~ -1` modes of the lagged-coupling map. A spatial
  grid hierarchy addresses only (a)-type spatial modes; the
  published methods coarsen SPACE on Cartesian grids and do
  not treat a comoving frequency-direction march. A per-ray
  space-frequency multilevel would be a new method, not an
  import.
- For Levels 2-3: a B&L-type FAS multigrid is conceivable on
  the outer CMFGEN map but would need coarse-grid versions of
  the full opacity/population/SE stack - a major architectural
  change, and the non-monotonic frequency coupling is outside
  every published variant.
- Verdict: NOT categorically inapplicable - it is a proved
  technique for the coupled NLTE problem - but irrelevant to
  the Level-1 defect being addressed here and out of scope for
  the planned phases. Status stays: deferred research item,
  recorded accurately this time.

## 11. Direct comparison table (corrected)

| Method | Paper | Assumptions | Needs | Current code | Deficiency | Proposal | Risks | Tests |
|---|---|---|---|---|---|---|---|---|
| Ray-local AA-II | Toth-Kelley, Walker-Ni | LOCAL contraction, bounded alpha, same map | Per-ray residual history, Gram solve | `SWEEP_FLD_ORD` 7949-8823, `KN_ORD_AA` opt-in, m=2 | No column-independence test on D; no explicit `||alpha||` bound | Add `||d_hat||>=tau||d||` (Type-II analogue) + `||alpha||_1` bound | False accel on near-parallel history | Slow-converging + oscillatory rays; AA on/off compare |
| Zhang-type stabilization | Zhang et al. | Non-expansive `f` in l2 (NOT established) | Type-I `B_k` update | N/A - code is Type-II | Theorem does not apply; machinery is Type-I-only | Borrow architecture only: independence test + safeguard acceptance | Mistaken theorem claims | Gram conditioning test |
| Residual acceptance | EW94-style | True residual, unchanged state on reject | Raw `RES` metric | Pending-mix + next-pass `RES<=AA_REF` + restore 8443-8463 | Acceptance uses weighted-max norm while Gram optimizes weighted-l2 - stricter than the fit, acceptable but must be documented | Keep; document norm mismatch | - | Rejected-step counts |
| Evans rate improvement | Evans et al. | GLOBAL kappa<1, g C^2, bounded alpha != 0 | Linear regime evidence | `KN_ORD_AA` opt-in default off | Enable criterion must come from Level-2/3 evidence, not single ratios | Phase-1/2 diagnostics gate | Enabling in non-contractive regime | Slow + oscillatory A/B |
| Adaptive omega | KM/Zhang Sect.1 spirit | Same-map ratio | Same-ray residual ratio | `RATIO_OM` same-direction/consecutive + `OMEGA_RAY` bounds [0.45,1] 8539-8553; cross-call `INEXT_RELAX` median envelope 6547-6588 | Global-max envelope mixes rays when argmax moves | Gate update on argmax stability | Mis-tuned damping | Omega history audit |
| MOM outer AA | (none) - violates same-map | One operator per history | - | `MOM_AA_OUT` default T, SAVE history 6548-6750 | Secants span population updates = different maps; only omega-rescale | Invalidate on operator generation; consider default-off | Extrapolation across maps | Generation-bump test |
| JFNK forcing | EW96 Choice 2 | Accepted-residual ratio, eta_max<1 | Outer `||F||` history | `jfnk_solver_v1.f90` 1019-1027, safeguard + bounds verified | eta_prev persists across solves (harmless heuristic) | Audit only | - | JFNK unit tests |
| JFNK globalization | EW94 | True-residual accept, clean rollback | `f()` eval at trial | ared/rho + TR + retry + restore 1497-1720 | Residual mode must be fixed within a solve | Phase-6 verify only | - | Forced-rejection test |
| Dogleg | Pawlowski | Trust region | `J^T F` product | Hookstep TR covers | None found yet | Defer | - | - |
| Multigrid | B&L17, F&B97 | Coupled SE+RT map, structured grids | Coarse-grid operator stack | - | Targets Levels 2-3, not Level-1; no comoving-frequency variant | Defer, recorded accurately | - | - |

## 12. Explicit recommendations (revised)

1. Keep ordered GS (`SWEEP_FLD_ORD`, `KN_RAY_ORDERED=T`) as the
   baseline inner iteration - Knop Paper II physically ordered
   march and the contraction the accelerators rely on.
2. Keep AA scope at per-(ray, direction-field, call) - the
   implemented granularity is already finer than needed and is
   the correct scope per Walker-Ni/Toth-Kelley. Never widen.
3. Add the Type-II column-independence test
   (`||d_hat|| >= tau*||d||` on the residual-difference column
   entering the Gram, or equivalently a minimum-pivot/rcond
   test) plus an explicit `||alpha||_1` bound. These are the
   correct analogues of Zhang's restart machinery for the
   implemented Type-II system - describe them as such, not as
   the Zhang-Boyd algorithm.
4. Keep the pending-mix + next-pass raw-residual acceptance
   (verified to measure the undamped residual). Document that
   acceptance is in the weighted-max norm while the Gram fit
   is weighted-l2.
5. `KN_ORD_AA_M = 2` initial depth is literature-consistent;
   raise only with conditioning data.
6. `MOM_AA_OUT` history must be invalidated on operator
   generation change (population commits); its default-on
   status contradicts the FG-side deferral - review default.
7. Keep JFNK Choice 2 (verified faithful); Phase-6 audit of
   parameters and residual-mode consistency only.
8. Dogleg and multigrid remain deferred; multigrid now
   correctly classified as a coupled-problem technique.
9. Establish the observed contraction regime (Evidence
   hierarchy Level 2-3) from multi-pass residual diagnostics
   before enabling AA anywhere; single-pass ratios and
   global-max envelopes are Level-1 evidence only.

## Nov-2026 audit-cycle addendum

The relaxation-controller correction this cycle is consistent with
the fixed-point literature already reviewed:

- Toth & Kelley (2015), DOI 10.1137/130919398: convergence of
  Anderson-accelerated fixed-point maps is governed by the map's
  own contraction properties. The per-ray same-map residual ratio
  is a legitimate inner-map diagnostic; it cannot be re-used as an
  outer-map spectral-radius estimate because the outer operator
  includes the sweep write, the commit filter and population
  updates - a different map entirely.
- Walker & Ni (2011), DOI 10.1137/10078356X: residual-based
  progress measures are only meaningful when the map is held
  fixed. Pooling ratios across diverged rays or across fewer
  than `KN_RHO_MIN_SAM` samples violates that premise; the new
  pool filter and minimum-sample retain implement it directly.
- Evans et al. (2020), DOI 10.1137/19M1245384: acceleration
  improves linearly converging maps, not quadratically converging
  ones. The KNOP ray solve is linear-class; the heuristic stays
  justified only as long as its assumptions (enough valid samples,
  same map, no divergence) are met.
- Eisenstat & Walker (1996), DOI 10.1137/0917003: inexact-Newton
  forcing requires deterministic residual evaluations. The JFNK
  eval state bracket (`FG_EVAL_STATE_PUSH_V1`/`POP_V1`) now
  guarantees eval sweeps cannot perturb production telemetry, so
  residual determinism is no longer polluted by census resets.

No published theorem is claimed to prove convergence of the
CMFGEN implementation; these references bound what the empirical
controller may legitimately claim.
