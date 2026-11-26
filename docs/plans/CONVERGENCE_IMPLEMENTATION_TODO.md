# Convergence Implementation TODO - Non-Monotonic CMFGEN

Status file for the phased convergence work derived from
CONVERGENCE_LITERATURE_REVIEW.md and CONVERGENCE_THEORY_AUDIT.md.
Every item lists problem, source, justification, reference, proposal,
tests, acceptance criteria, status.

Statuses: NOT STARTED / IN PROGRESS / COMPLETE / BLOCKED.

Audit revision note (commit cfae5852): items whose justification
changed after the critical review are marked REVISED. Two items
were reclassified (see notes in 5.1 and the removed Powell
proposal) and no item is based on the Zhang-Boyd global theorem,
which does not apply to the CMFGEN Type-II implementation.

Baseline to preserve: `KN_RAY_ORDERED=T` ordered GS march,
`KN_ORD_AA=F`, `KN_ORD_AITKEN=F`, `KN_ORD_P2DAMP=F`,
`KN_ADAPT_RELAX=T`, `JFNK_FADAPT=T`. All new behavior behind VADAT
switches; disabled state must reproduce today's bit-path.

## PHASE 1 - Convergence diagnostics audit

### 1.1 Distinguish the five residuals in diagnostics
- Problem: refresh residual, sweep drift, store churn, retained-
  field delta and population demand are printed through overlapping
  channels; a damped update can masquerade as convergence in logs.
- Source: `cmfgen/plane/fg_j_cmf_nonmono_v1.f` (refresh prints),
  `cmfgen/plane/mom_j_cmf_nonmono_v1.f` (MOM sweep/refresh prints),
  `cmfgen/new_main/mod_subs/solve_for_pops.f` (GS_DRIFT_DIAG).
- Justification: spec 3.1/3.4; `||x_{k+1}-x_k|| = omega||r_k||` must
  never substitute for `||r_k||`.
- Reference: Toth-Kelley (residual-based convergence), spec.
- Proposal: verify every print names its metric unambiguously
  (`res=`, `sw=`, `rf=` already do in GS_DRIFT); add the same
  labeling to refresh summary lines; no tolerance change.
- Regression tests: string-level check of OUTGEN/LUER records on a
  short model run; unit test on `KN_RES_PASS` normalization.
- Acceptance: each convergence channel in diagnostics maps
  one-to-one to a defined residual in THEORY_AUDIT Sect. 2.
- Status: NOT STARTED

### 1.2 Verify 1e-4 significance floor usage
- Problem: the col-peak floor excludes weak cells from max/counts;
  must not hide weak-intensity errors.
- Source: `fg_j_cmf_nonmono_v1.f` (RES loops), MOM analogues.
- Justification: spec 3.4; floor is a significance filter, not an
  acceptance floor.
- Reference: spec; Bergemann-Hoppe Sect. on weak-field regions.
- Proposal: document floor semantics; add a diagnostic count of
  cells in [tol, floor] band so excluded-but-large cells are visible.
- Tests: synthetic field with error only in sub-floor cells.
- Acceptance: floor cells reported, not silently dropped from all
  observables.
- Status: NOT STARTED

## PHASE 2 - Ray-local residual and contraction tracking

### 2.1 Expose per-ray contraction factor
- Problem: contraction currently estimated from global max ratios
  (`KN_RF_RESID/KN_RF_RESID_PRV`); unrelated rays can dominate
  successive calls.
- Source: `KN_RES_PASS`, `KN_RAY_RES`, `KN_RF_RESID`,
  `KN_RF_RESID_PRV` in `fg_j_cmf_nonmono_v1.f`; MOM `RES_PASS`,
  `MOM_RF_RESID_PRV`.
- Justification: spec 3.2; Toth-Kelley contraction is local.
- Reference: Toth-Kelley Thm 2.6; spec "no contraction from
  unrelated ray maxima".
- Proposal: compute `rho_ray = RES_pass(k+1)/RES_pass(k)` inside
  `SWEEP_FLD_ORD` (same ray, same C-pack, same direction set);
  export a per-call distribution (median/max rho over converged-
  vs-uncoverged rays) for diagnostics; keep global ratio only as a
  divergence detector.
- Tests: manufactured two-rate rays; verify rho_ray tracks the slow
  ray, not the max.
- Acceptance: contraction numbers are ray-consistent; argmax moves
  are flagged, not folded into the ratio.
- Status: NOT STARTED

### 2.2 Fix MOM omega estimator metric mixing
- Problem: `DTR = MOM_GS_MAXCHG / MOM_RF_RESID_PRV` divides sweep
  drift by refresh residual - two different metrics (Audit C1).
- Source: `cmfgen/plane/mom_j_cmf_nonmono_v1.f` ~6931.
- Justification: contraction requires same-map, same-metric ratios.
- Reference: Toth-Kelley; Audit Sect. 4.2.
- Proposal: use `MOM_RF_RESID / MOM_RF_RESID_PRV` (refresh/refresh),
  matching the FG `KN_RF_RESID` ratio; keep sweep drift for churn
  gates only.
- Tests: unit test feeding known residual sequences; confirm omega
  responds to refresh-residual contraction.
- Acceptance: `JNEXT_RELAX` responds only to same-metric ratios.
- Status: NOT STARTED

### 2.3 Gate omega adaptation on argmax stability
- Problem: global ratio adapts omega even when the dominant cell
  moved between calls (Audit C3).
- Source: `FG_ARG_HK/HR/HL/HV`, `GET_FG_GS_CLS`; MOM argmax slots.
- Justification: ratio is a contraction only if argmax stationary.
- Reference: Audit Sect. 4.1; spec 3.2.
- Proposal: skip omega update (or hold last value) when the argmax
  tuple changed cell between calls; reuse `cycle`/`margnl`/`bndry`
  classification.
- Tests: alternating-argmax residual sequence; omega must not move.
- Acceptance: omega changes only on stable-argmax contractions.
- Status: NOT STARTED

## PHASE 3 - Relaxation control audit

### 3.1 Audit `INEXT_RELAX`/`JNEXT_RELAX` bounds and schedules
- Problem: adaptive omega can drift without a published bound;
  omega in (0,1) keeps the relaxed map non-expansive (KM iteration).
- Source: setters + adaptive blocks (FG ~6582, MOM ~6931);
  `KN_ADAPT_RELAX` default T.
- Justification: Zhang-Boyd KM safeguard uses alpha in (0,1);
  under-relaxation is the convergence backstop.
- Reference: Zhang-Boyd Sect. 1 (KM), Audit 4.
- Proposal: enforce `omega_min <= omega <= 1` explicitly (bounds
  exist in MOM `MAX(0.05,MIN(1,...))`, confirm FG parity); log every
  omega change with its measured ratio.
- Tests: omega bounds unit test; history log check.
- Acceptance: omega bounded, every change traceable to a same-map
  ratio.
- Status: NOT STARTED

## PHASE 4 - FG convergence state machine

### 4.1 Codify per-ray solve states
- Problem: ray outcomes (converged, capped, diverged, stalled,
  pass-limit, plateau) exist but as overlapping flags; a formal
  state machine prevents contradictory exits.
- Source: `KNOP_RAY_SOLVE_V1`, `IXIT` codes, commit verdicts in
  `FG_KNOP_REFRESH_V1`.
- Justification: spec Phase 4; commit logic must see one verdict.
- Reference: existing `KNOP_AUDIT_TODO.md` structure; spec.
- Proposal: single per-ray exit state {CONV, DIVERGED, STALLED,
  CAPPED, MAXPASS, UNSOLVED} mapped from current flags; commit and
  heal read only that state.
- Tests: enumerate exit paths in unit harness; each maps to one
  state.
- Acceptance: no ray committed without state CONV; verdict counts
  agree with state histogram.
- Status: NOT STARTED

## PHASE 5 - Safeguarded Anderson acceleration

### 5.1 Add column-independence control to `SWEEP_FLD_ORD` AA (REVISED)
- Problem: near-parallel residual-difference columns
  `d_j = r_1 - r_{j+1}` enter the Gram and are only caught
  post-mix by the residual rollback (Audit Sect. 6.1 gap 1).
- Source: `SWEEP_FLD_ORD` Gram build, `cmfgen/plane/fg_j_cmf_nonmono_v1.f`
  8703-8752; history push 8072-8084.
- Justification: for the implemented Type-II system the
  conditioning of `GRM = D^T W D` is what controls `||gamma||`;
  the correct control is a linear-independence test on the
  residual-difference columns, NOT Zhang's `||s_hat||` rule
  (which tests iterate-secants for the Type-I `B_k` update this
  code does not have). Toth-Kelley Thm 2.6 needs bounded
  coefficients; near-dependent `D` columns is exactly how that
  fails.
- Reference: Toth-Kelley (conditioning of `F_k`); Zhang et al.
  2020 Sect. 3.2 (role analogue only - different matrix, not the
  published criterion).
- Proposal: at history push or Gram build, compute the
  orthogonal component `d_hat` of the newest `d` column against
  the stored `D` span (MGS on <=8 columns) and require
  `||d_hat|| >= KN_ORD_AA_TAU*||d||` (default ~1e-2); on failure
  skip the push or drop the oldest column. Alternative with the
  same role: minimum-pivot/`rcond` test on `GRM` after the
  ridge. VADAT `KN_ORD_AA_TAU` + reader + consistency check +
  vadat_keys.
- Tests: rank-deficient history case; verify the test fires
  before `DGESV`.
- Acceptance: Gram solve never sees an uncontrolled
  near-dependence; `KN_AA_REJ` not spent on preventable mixes.
- Status: NOT STARTED

### 5.2 Bound Anderson coefficients explicitly
- Problem: `||alpha||` unbounded => unstable extrapolation
  (Toth-Kelley assumption); currently enforced only via the
  `1e3*||f||` correction veto.
- Source: Gram solve + correction veto, `SWEEP_FLD_ORD`
  8751-8793.
- Justification: Toth-Kelley Thm 2.6 requires bounded alpha.
- Reference: Toth-Kelley Thm 2.6 only. NOTE REVISION: the
  previous "Powell-style mixing parameter" alternative was
  removed - Powell regularization is specific to the Type-I
  rank-one `B_k` update and has no Type-II counterpart here.
- Proposal: after Gram solve, veto the mix if
  `sum_j |gamma_j| > KN_ORD_AA_ABND` (initial ~10; the
  alpha_0 = 1 - sum gamma weight is implied) and count it in
  the veto channel.
- Tests: adversarial history producing large alpha; verify veto.
- Acceptance: no committed mix with `||alpha||_1` above bound.
- Status: NOT STARTED

### 5.3 Invalidate `MOM_AA_OUT` history on operator change
- Problem: outer AA secants span population updates = different
  maps (Audit C2). Literature requires one operator per history.
- Source: `XHIST_OUT/FHIST_OUT`, `NHIST_OUT` in
  `mom_j_cmf_nonmono_v1.f` ~6548-6960; population-commit signal in
  `solve_for_pops.f`.
- Justification: Walker-Ni/GMRES equivalence needs a fixed linear
  map; Zhang-Boyd secants are defined w.r.t. one `f`.
- Reference: Walker-Ni; Zhang-Boyd; Audit Sect. 6.2.
- Proposal: (a) simplest - set `MOM_AA_OUT` default F pending (b);
  (b) proper - export an operator-generation counter bumped on any
  committed population change; clear `NHIST_OUT` when it advances.
  Decide default after Phase-1 diagnostics show whether the stored
  deltas are in fact contractive.
- Tests: history-clear unit test on generation bump; A/B run with
  `MOM_AA_OUT` on/off on the slow-convergence model.
- Acceptance: no secant pair spans a population commit; AA never
  fires on a stale-generation history.
- Status: NOT STARTED

### 5.4 Codify AA vs P2-damping precedence
- Problem: AA mix can re-excite a P2-damped cell (Audit S2).
- Source: `SWEEP_FLD_ORD` AA block + P2 damp arrays.
- Justification: period-2 = eigenvalue near -1; AA on such cells
  produces singular Grams (Evans caveat).
- Reference: Evans et al. 2020; Audit 6.3.
- Proposal: exclude P2-latched cells from the affine subset (same
  mechanism as AA_GCAP exclusion).
- Tests: oscillatory-convergence case with both on.
- Acceptance: no AA mix touches a P2-damped cell.
- Status: NOT STARTED

### 5.5 Enable-criteria study (diagnostics only) (REVISED)
- Problem: spec requires AA opt-in only after the observed regime
  is established. The Evans result is conditional on a globally
  contractive, C^2 operator - not established for `G_ray`; the
  correct enable criterion is the evidence hierarchy in
  THEORY_AUDIT Sect. 8 (multi-pass stable ratios + stable argmax +
  bounded `||alpha||`), NOT a single residual ratio.
- Source: `KN_RES_PASS` histories (per-pass max envelope, audit
  K), `KN_RES_PASS`-style per-ray pass residuals `RES_LS`,
  `GET_FG_GS_CLS` classes, `RATIO_OM` history.
- Justification: Evans Thm 4.1 improves the leading rate only
  where the linear-contractive regime holds; observed behaviour
  must reach Level-2/3 evidence first.
- Reference: Evans et al. 2020 (corrected form: rate factor
  `theta_k`, `theta_k = sqrt(1-||w_R||^2/||w_{k+1}||^2)`).
- Proposal: report per-ray rho class + argmax stability +
  `||alpha||` statistics + AA on/off counterfactual on the
  slow-convergence model; produce a documented enable criterion.
- Tests: A/B on slow + oscillatory models.
- Acceptance: documented criterion at Level-2 evidence or above;
  AA stays default-off.
- Status: NOT STARTED

## PHASE 6 - JFNK forcing and globalization audit

### 6.1 Audit EW forcing inputs and parameter values
- Problem: confirm `||F_k||/||F_{k-1}||` is measured on the same
  residual definition each step (frozen-J vs coupled modes);
  confirm gamma/eta_min/eta_max vs measured Krylov counts.
- Source: `jfnk_solver_v1.f90` ~981-1027; `jfnk_solve_corr_v1.f`
  mode selection; `JFNK_GATE_IN` hysteresis.
- Justification: EW96 Choice 2 defined on nonlinear residual
  ratios; mode mixing breaks the sequence.
- Reference: Eisenstat-Walker 1996 Sect. 2.
- Proposal: log forcing inputs; unit test fixed + adaptive paths;
  if mode can switch mid-Newton, pin residual mode per solve.
- Tests: JFNK unit tests (existing harness); Krylov-count report.
- Acceptance: eta sequence computed on one residual definition;
  documented param rationale.
- Status: NOT STARTED

### 6.2 Verify residual/state contamination on rejected steps
- Problem: rejected Krylov directions must leave nonlinear state
  and history clean (EW94 requirement).
- Source: retry/restore blocks ~1400-1469 `jfnk_solver_v1.f90`;
  resacc snapshot/restore in `solve_for_pops.f`.
- Justification: EW94 accept/reject contract.
- Reference: Eisenstat-Walker 1994.
- Proposal: audit all early-exit paths for partial state writes;
  test rejected-step leaves `F` evaluation unchanged.
- Tests: forced-rejection unit test.
- Acceptance: zero state mutation on rejected steps.
- Status: NOT STARTED

### 6.3 Deferred candidates
- Dogleg (Pawlowski): only if 6.2 finds a real acceptance gap.
  Status: NOT STARTED (deferred)
- Multigrid (B&L 2017, F&B 1997): established technique for the
  COUPLED SE+RT problem (Levels 2-3), not the Level-1 refresh;
  no published variant treats a comoving sign(GAM) frequency
  march or the CMFGEN ray-by-ray structure. Reclassified from
  "inapplicable" to "wrong level for this work".
  Status: NOT STARTED (deferred research item)

## PHASE 7 - Level integration

### 7.1 Gate population update on Level-1 convergence quality
- Problem: SE updates on immature fields re-poison the map
  (ordering rule); current gates exist but their thresholds are
  empirical.
- Source: `GS_FREEZE` controller, `GS_REL_STABLE`,
  `GS_MATURE_TOL`/`GS_UNSAFE_TOL`, `GS_HOLD_POPS`, JFNK admission
  tier.
- Justification: Knop Paper II ordering; spec 3.3.
- Reference: Knop et al. 2009b Sect. 3.2.
- Proposal: wire Level-1 per-ray convergence census (Phase 4
  states) into the maturity gate so the gate reads solved-state,
  not only residual size; keep existing thresholds.
- Tests: deliberately slow-FG model; confirm populations held until
  ray census clean.
- Acceptance: no SE update commits on a field with unconverged
  rays above census bound.
- Status: NOT STARTED

### 7.2 Final frozen-coefficient verification solve
- Problem: outer exit (`MAXCH < EPS`) can fire on a field whose
  last refresh pre-dates the last population move.
- Source: refresh call site vs `COMP_J_CONV` ordering in
  `cmfgen_sub.f`.
- Justification: fixed-point identity must hold at the accepted
  solution, spec item 10.
- Reference: spec; EW94 final-residual check analogue.
- Proposal: on converged exit, run one refresh with frozen final
  coefficients and report residual; diagnostics-first (no exit
  veto until measured).
- Tests: converged model; report residual distribution.
- Acceptance: exit-time refresh residual reported and bounded.
- Status: NOT STARTED

## PHASE 8 - Regression and stress matrix

### 8.1 Model matrix
- Cases: monotonic velocity; single reversal; multiple reversals;
  strong interaction branches; weak-intensity regions; thick/thin
  cells; deliberately slow FG; oscillatory FG; AA-rejection case;
  ill-conditioned AA history; operator-change-through-population;
  OBS-grid != main grid.
- Source: `tests/`, `unit_tests/`, slurm run templates.
- Justification: spec Sect. 8.
- Reference: spec.
- Proposal: assemble per-phase test list above into the config
  matrix; compare vs `KN_ORD_AA=F`/`KN_ADAPT_RELAX` baselines.
- Metrics: final fields, pops, spectra, conservation checks,
  convergence histories, rejected-update counts, formal-solve
  counts, runtime/memory.
- Acceptance: no regression vs baseline with features off;
  documented deltas with features on.
- Status: NOT STARTED

## Dependency notes

- 2.2 is the only pre-Phase-5 confirmed code defect; fix first
  (small, isolated).
- 5.3 depends on 2.x diagnostics (need measured contractivity
  before deciding default).
- 7.x depends on 4.1 (state machine exists to be read).
- All AA work stays behind `KN_ORD_AA`/`MOM_AA_OUT` switches;
  Phase-5 changes must not alter the disabled path.

## Change log

- 2026 audit documents created: CONVERGENCE_LITERATURE_REVIEW.md,
  CONVERGENCE_THEORY_AUDIT.md, this file. No source changes.
- Critical audit revision (commit cfae5852 re-verified): 5.1
  rewritten as the Type-II column-independence test (the
  Zhang `||s_hat||` criterion does not apply to this Gram);
  Powell-regularization alternative removed from 5.2 (Type-I
  only); 5.5 enable criterion rebuilt on the audit's evidence
  hierarchy; 6.3 multigrid reclassified as a coupled-problem
  technique deferred for scope reasons, not dismissed as
  inapplicable. No source changes.

## PHASE A - Critical source audit corrections (Nov-2026 cycle)

Independent-source audit verified against `main` at `7220b616`;
each issue implemented in its own phase and committed separately.
Unit coverage in `unit_tests/python/test_knop_ray.py` (T19-T32),
driver `unit_tests/f90/drv_knop_ray.f90`.

### A.1 OpenMP race in FG_JUMP_MAX sibling scan
- Problem: `FG_JUMP_MAX` read sibling `I_*_DEPTH_NEXT` columns from
  inside the parallel ray loop while other threads committed into
  them; the `V2 .EQ. V2` finite test also admitted `+Inf`.
- Root cause: shared store read/write race; non-finite guard that
  only excluded NaN.
- Correction: `KN_JMP_SIBREF_V1` snapshots an immutable per-
  frequency sibling reference (median of finite positive non-floor
  per-ray column peaks) before the ray loop into `JMP_SIB_M/P`;
  `FG_JUMP_MAX` reads only the snapshot and rejects non-finite
  field/store values explicitly.
- Tests: T19-T21 (1/2/4-thread determinism, poisoned siblings,
  floor columns); `sibref`, `jmpx` driver cases.
- Acceptance: verdicts, deltas and retained stores bit-identical
  across thread counts; Inf/NaN siblings cannot disable the veto.
- Status: COMPLETE (`b49b32d3`, `f41a69c0`, `8c4983f9`)

### A.2 INEXT_RELAX estimator + small-sample P90
- Problem: pooled residual-ratio P90 collapses to the pool maximum
  for L<10 samples (index `(9*L+9)/10 = L`), so one outlier set the
  global relaxation; inner-map ratios were treated as an outer-map
  contraction factor.
- Root cause: degenerate quantile on small samples.
- Correction: `KN_OMEGA_POOL_V1` pools only finite positive ratios
  from solved non-diverged rays; below `KN_RHO_MIN_SAM=10` samples
  the previous `INEXT_RELAX` is retained and coverage reported.
  `omega = 1/(1+rho_p90)` remains an explicitly empirical heuristic,
  not a claimed optimal damping factor (see THEORY_AUDIT Sect. A).
- Tests: T18 (12-sample genuine update), T22-T25 (contracting,
  expansive, L=1/8 retain, outlier rejection).
- Status: COMPLETE (`4d3b8062`, `148002b2`, `5b439f6a`)

### A.3 Jump-rejected rays reported false maturity
- Problem: `FG_RAY_VC=1/2` rays left the retained store unchanged
  and contributed nothing to `FG_RF_DCHG`; a converged inner solve
  could be vetoed yet the gate saw a mature field.
- Correction: `FG_RF_PEND` carries the max pending proposed-vs-
  retained delta over verdicts 1/2; `GET_FG_GS_MATURITY` exports
  `RF_DCHG = MAX(FG_RF_DCHG, FG_RF_PEND)` into `GS_MAX_ALL`.
- Tests: T26, T27 (`jmpx`, `jmp2`; repeated rejection stays
  immature).
- Status: COMPLETE (`b3ba07b3`, `00c4ffc8`, `61fe700a`)

### A.4 Unsolved early-return rays read as converged
- Problem: `NI<2 .OR. NCF<2` returned `RESID=0`, `NXIT=-1` -
  construction failures and legitimately inactive single-point
  geometries were indistinguishable.
- Correction: split classes. `NI<1 .OR. NCF<2` is a failed build:
  `NXIT=-2`, `RESID=1.0` fails closed into the residual max.
  `NI=1` stays `NXIT=-1`, `RESID=0` (no marchable domain). Exit
  census arrays extended to `(-2:6)`; failed builds printed
  separately.
- Tests: T28 (`ni0` fails closed), T29 (`ni1` excluded).
- Status: COMPLETE (`478280a1`, `038e0c72`, `91fd8204`)

### A.5 HEAL_STORE counted already-zero slots
- Problem: median-branch heal incremented `NHL` whenever a slot
  failed the invalid test, including already-zero slots.
- Correction: heal only fires on nonzero invalid slots; algorithm
  unchanged, count now equals actual modifications.
- Tests: T30 (`heal`; mixed poisoned store, count 15 not 21).
- Status: COMPLETE (`a58f0ba3`, `e4274a49`, `d32ffe84`)

### A.6 JFNK eval sweeps reset production FG telemetry
- Problem: a JFNK eval runs a full sweep whose INIT pass clears
  every maturity census field; entry saved only five relaxation
  scalars, so the gate and FGINEXT_GS checkpoint could read
  eval-side state.
- Correction: `FG_EVAL_STATE_PUSH_V1`/`POP_V1` in
  `FG_J_CMF_NONMONO_MOD_V1` snapshot/restore the full telemetry
  set (census, deltas, tallies, argmax history, error-record
  COMMON); called at the entry saves and epilogue restores of
  `JFNK_RAD_SWEEP_V1`.
- Tests: T31 (`esav` sentinel round-trip).
- Status: COMPLETE (`e8fc9c1f`, `aa1e7284`, `0f67b509`, `f8c04ca3`)

### A.7 DTAU clamp parity (B1 close-out)
- Problem: refresh `BUILD_FREQ` skipped the min/max DTAU clamp
  under `INSERT`/`METHOD(4:6)='MON'`, citing the DIFFERENCE-path
  rule; its actual sibling is the `SOLUTION_METHOD='INTEGRAL'`
  sweep block, which clamps unconditionally because chi_hat
  carries the `|GAM|/DLN` coupling term monotonic CHI interpolation
  cannot fix. Latent: inactive under `METHOD='ZERO'`.
- Correction: unconditional clamp in the non-ZERO branch of
  `BUILD_FREQ`, matching the INTEGRAL-path operator.
- Tests: T32 (`pmon`/`pnom` identical output).
- Status: COMPLETE (`e8fc9c1f`, `0f67b509`, `f8c04ca3`)

## PHASE B - Follow-up audit corrections (Nov-2026, e6e8f47a)

### B.1 Fail-safe sibling reference
- Problem: `KN_JMP_SIBREF_V1` returned REF=-1 on VPK allocation
  failure, indistinguishable from a genuine cold start; the refresh
  caller also returned silently on JMP_SIB allocation failure, so
  small sweep metrics could release populations against an
  unvalidated field.
- Correction: `ISTAT` out-arg (0 ok / 1 alloc fail / 2 malformed);
  REF=-1 now means cold start only. `FG_RF_FAILED` marks an aborted
  refresh and `GET_FG_GS_MATURITY` folds it into RF_DCHG=1.0 --
  fail closed. `FG_JMP_FAIL_INJ` is the deterministic injection
  hook for the unit tests.
- Tests: T33 (sfail ISTAT=1), T34 (scold ISTAT=0 all -1),
  T35 (rfail: RF_FAIL=T, GS_OK=F, GS_RFD=1.0).
- Status: COMPLETE

### B.2 DTAU parity strengthened to shared operator
- Problem: pmon/pnom compared two refresh calls, not the refresh
  against the in-sweep INTEGRAL operator.
- Correction: `KN_DTAU_INT_V1` extracted -- the Euler-Maclaurin
  cubic plus unconditional clamp -- and called by BOTH the
  `SOLUTION_METHOD='INTEGRAL'` sweep block and refresh
  `BUILD_FREQ`. Parity is now a construction fact.
- Tests: T37 (`dtau` probe, bounds enforced on monotone and
  non-monotone chi); T32 retained.
- Status: COMPLETE

### B.3 Eval bracket exit-path hardening
- Problem: `FG_EVAL_STATE_PUSH_V1` sat before the last `GOTO 9999`
  entry exit; a late entry alloc failure could leave a pushed
  state unrestored (harmless but fragile).
- Correction: push moved after the LSTATE block -- every early
  exit now precedes the armed push; the pop stays in the shared
  epilogue covering normal and broken-eval returns.
- Tests: T36 source pin (push after last early-exit, no GOTO-9999
  between push and pop).
- Status: COMPLETE
