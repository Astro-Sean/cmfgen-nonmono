# Numerical methods development plan

## Purpose and status

This plan turns the proposed convergence review into a sequence of
investigations with implementation limits and acceptance tests. It is a
research plan, not evidence that any solver change improves a physical model.
Complete one investigation at a time. Keep the current production path as the
comparison baseline.

Code snapshot reviewed: `cmfgen_jfnk` commit `e4c8f02b` on `main`. The
worktree had an in-flight change in `new_main/subs/cmf_blkband_v3.f`. With
explicit user approval, its row and column substitution was narrowed to
exactly uncoupled zero equations with zero RHS. The production helper has a
focused test; the full build and the 20-file unit suite passed on
03-Nov-2026. Recheck source and tests at the start of each phase.

## Current implementation and evidence

Several ideas from the initial review already have code or diagnostic support.
The next work is to establish their limits and validate their effect:

- `new_main/mod_subs/jfnk_solver_v1.f90` uses forward finite-difference JVPs,
  Eisenstat-Walker forcing, restarted GMRES, globalization, and true-residual
  checks. The source records cases where a Hessenberg estimate understated the
  measured linear residual.
- `new_main/mod_subs/jfnk_cmfgen_interface_v1.f` has an experimental BA
  matrix-vector path and banded BA data. The BA matrix also supplies a
  preconditioner; its accuracy as the derivative of each residual mode must be
  assessed separately.
- The preconditioner code already has block LU, banded LU, and an `n,T`
  Schur option. The current Schur option is not a coupled radiation/population
  Schur complement.
- `unit_tests/python/test_jfnk_solver.py` marks several LU solve-residual
  checks as expected failures. Treat these as open diagnostic gaps, not as
  passed tests.
- `plane/mom_j_cmf_nonmono_v1.f` contains inner and outer Anderson acceleration
  in the Knop refresh code. The formal path also has an ordered-pass
  acceleration option. Measure these implementations before proposing another
  accelerator.
- The moment solver uses `MOM_DLN_MIN` in frequency-advection coefficients.
  Existing unit tests cover a manufactured Knop ray and moment-solver pieces,
  but do not establish frequency-integrated conservation for a nearly
  degenerate frequency grid.
- The formal and moment paths have field-maturity diagnostics and gates in the
  population-solve workflow. Test whether their metrics track the consistency
  required by the subsequent population update.
- `unit_tests/README.md` lists production-source tests for JFNK, Knop rays,
  moment solves, and VADAT parsing. These tests do not establish convergence
  of a complete interacting-supernova model.

Progress: the opt-in `JFNK_HADAPT` step selector now keeps its selected
finite-difference multiplier scoped to the residual mode and baseline that
produced it. `JFNK_INIT_V1`, `JFNK_REBUILD_BASELINE_V1`, and a
frozen/coupled mode change re-arm the h sweep instead of reusing a stale
step. The Phase 1 diagnostic keys now have VADAT parser and source
regression coverage. This is an operator-state fix, not evidence of model
convergence.

The BA null-equation helper is tested by
`unit_tests/python/test_cmf_blkband_pin.py`. It pins only a fully disconnected
row and column with zero RHS. The focused test passes and the full CMFGEN
build is clean.

Useful starting points are `cmfgen_obsidian/solvers/JFNK solver.md`,
`cmfgen_obsidian/solvers/Knop formal solver.md`,
`cmfgen_obsidian/solvers/Sobolev non-mono solver.md`, and the current audit
notes under `docs/audits/`.

## Rules for all phases

1. Preserve the stock path when new flags are off. Compare against
   `../cmfgen_reference/cur_cmf` where a reference comparison is possible.
2. Define the physical residual, unknown ordering, held variables, masks,
   norms, and units before changing a solver.
3. Report the original physical residual after each accepted nonlinear update.
   A preconditioned or Hessenberg estimate is not a substitute.
4. Do not add a clamp, floor, limiter, or fallback to obtain convergence. Audit
   existing bounds at their point of use and state the equation they alter.
5. Keep experiments in unit-test drivers or developer diagnostics until they
   pass the acceptance tests. Do not add VADAT keys for internal tuning
   without a separate parameter proposal.
6. Record failed runs and safeguard activations. Do not describe a capped,
   deferred, or stalled run as converged.
7. Change one numerical component at a time. Keep the comparison model,
   compiler, inputs, and solver flags fixed between A/B runs.

## Phase 0 - Establish a reproducible baseline

### Questions

- Which failure modes reproduce in the current checkout?
- Which model and solver settings expose them?
- Are current logs sufficient to distinguish physical residual failure,
  field-map stagnation, JFNK linear-solve error, and a safeguard exit?

### Work

- Pin a small set of model inputs and capture the source revision, compiler,
  build flags, VADAT, velocity profile, frequency grid, and relevant output
  histories.
- Run the existing focused production-source tests and retain their full
  output, including expected failures.
- For the model runs, save outer residuals, field changes, formal/moment
  refresh counts, Newton and GMRES histories, true linear residuals,
  preconditioner status, rejected steps, and safeguard counts.
- Define the measured convergence quantities before examining A/B results.

### Acceptance

- A second run with the same setup reproduces the recorded convergence
  category and produces comparable histories.
- Every stop is classified as converged, stalled, iteration-limited,
  invalid-state, or failed solve. A limit exit never counts as convergence.
- The baseline includes at least one monotonic case and one non-monotonic case.

## Phase 1 - JFNK operator and preconditioner correctness

### Questions

- How much nonlinearity or noise does the production finite-difference JVP
  show across directions in one GMRES cycle?
- Does the selected perturbation scale resolve the derivative above residual
  evaluation noise and below finite-difference truncation error?
- For each frozen-J mode, which residual derivatives are represented by the
  BA matrix, and which terms are omitted?
- Which part of the known LU residual gap comes from assembly, scaling,
  factorization, or the test's surrogate system?

### Work

- Extend manufactured production-source tests to measure repeat-evaluation
  noise, directional additivity, positive/negative perturbation agreement,
  and the JVP error over a perturbation sweep. Exercise masks, near-boundary
  states, and branch changes. (Progress 03-Nov-2026: surrogate additivity
  and homogeneity now measured by T22; production-interface residual
  coverage still open.)
- Compare cycle-fixed perturbation policies with the current policy using the
  same set of Krylov directions. If a direction crosses a nonsmooth constraint,
  record the event and test a clean cycle restart. Do not quietly mix
  incompatible derivatives in an existing Arnoldi basis.
- Compare BA action against finite differences of the exact frozen-J residual
  on basis, random, and physically relevant directions. Document the residual
  mode and state dependencies for every comparison.
- Diagnose and resolve the LU expected failures. Check the assembled matrix,
  row/column scaling, pivot status, solve residual, and backward error before
  changing regularization.
- Compare diagonal-block and banded preconditioners using true preconditioned
  residuals, GMRES iterations, memory, and factor/apply time. Record the
  fraction and location of off-depth-band matrix entries.
- Keep the existing `n,T` Schur option distinct from any proposed
  radiation/population Schur method.

### Acceptance

- Tests establish expected JVP behavior over a documented range of step sizes
  and show where evaluation noise or branch changes dominate.
- The returned GMRES step passes a true-residual check using the production
  operator. Restarts and trust-region projections trigger a new check for the
  vector actually returned.
- The BA matrix is called an analytic JVP only for residual modes where the
  directional derivative test supports that claim. Otherwise it remains a
  preconditioner or explicitly approximate operator.
- ~~LU solve-residual tests no longer rely on XFAIL for cases classified as
  correctness checks. Any retained diagnostic limitation has a separate,
  explained expectation.~~ (done 03-Nov-2026: T4 machinery gates pass on
  all specs; T21 keeps the near-singular residual census as a named XFAIL.)
- The chosen preconditioner shows a measured benefit on a named case without
  worsening the physical nonlinear residual or changing the solved equations.
- With JFNK disabled, the stock population path remains unchanged.

### Defer

Do not add FGMRES unless the preconditioner actually changes during a single
Krylov solve. Do not replace GMRES or adopt BA as the production Jacobian based
only on agreement for a few directions.

## Phase 2 - Frequency operator and moment-solver conditioning

### Questions

- What limiting discrete equation does the moment solver approach as an actual
  frequency interval becomes very small?
- Does `MOM_DLN_MIN` alter frequency-integrated coupling, and which cells and
  frequency pairs activate it?
- What causes poor conditioning near `1 + sigma K ~= 0`: the closure, the
  discretization, the boundary treatment, or the iteration?

### Work

- Trace each `MOM_DLN_MIN` use from frequency-grid construction through
  `MOM_DLN_J`, `MOM_DLN_H`, and the resulting coefficients. Record units,
  interval distributions, activation counts, and affected terms.
- Build production-source tests with analytic and manufactured frequency
  grids: ordinary spacing, very small positive spacing, edge doublets, and
  repeated or invalid coordinates as rejected input.
- Measure frequency-integrated flux or the applicable discrete balance on
  both sides of the current threshold. Derive any replacement treatment from
  the limiting equation before implementation.
- Trace the near-resonant coefficient through the moment closure, angular
  terms, boundary conditions, and iteration. Compare the unmodified equation
  with any proposed numerical reformulation.
- Profile the current alternating and ordered solves, including existing
  Anderson acceleration. Add a block-frequency prototype only if the
  measured stencil and coupling show that a block solve addresses the
  diagnosed failure.

### Acceptance

- Tests verify the actual production discretization against a manufactured
  solution and check the relevant integrated balance across grid-spacing
  cases.
- Any replacement for `MOM_DLN_MIN` has a derived limiting form, explicit
  error measure, and comparisons below, at, and above the affected interval
  range. The existing guard is not removed before this test passes.
- Resonance handling is described mathematically and tested on either side of
  the condition. A larger regularization value is not accepted as a fix.
- A frequency-block prototype, if pursued, agrees with the existing solve on
  cases where both converge and improves a measured residual or cost on the
  target case without changing the discretized equations.

### Defer

Do not add generic minimum-spacing changes, positivity repairs, or new
regularization based on stability alone. Do not adopt a block solver before
the frequency stencil and boundary coupling are documented.

## Phase 3 - Radiation and population iteration

### Questions

- Do the formal/moment field maturity metrics predict a self-consistent field
  for the current population update?
- Which populations or thermodynamic variables cause coupled field/population
  oscillation?
- Does existing Anderson acceleration reduce the fixed-point residual, or
  merely change its trajectory?

### Work

- Trace the selected formal or moment solver, refresh, maturity decision,
  population update, and subsequent radiation recomputation for each supported
  mode. Record which quantities are frozen or lagged.
- Test maturity decisions against an independently recomputed field residual
  on representative monotonic and non-monotonic models. Separate field drift
  from statistical-equilibrium residual and from population correction size.
- Measure current inner and outer Anderson acceleration activation,
  acceptance, rejection, and fallback counts. Compare its residual and cost
  with the same map using the existing unaccelerated iteration.
- Test the sequence of radiation warm-up, frozen-J population solve, field
  refresh, and optional coupled correction. Base transitions on measured
  residual and consistency criteria, not a new fixed iteration-count rule.
- Prototype a coupled radiation/population Schur approximation only after
  Phase 1 establishes the population block and the relevant radiation response.

### Acceptance

- Population admission requires the documented field-maturity criteria and a
  direct check that those criteria correspond to a suitably stable radiation
  field on the tested models.
- The true statistical-equilibrium residual and radiation-field consistency
  are both reported at termination.
- Anderson acceleration preserves the fixed point on manufactured problems,
  rejects non-finite or worse candidates, and falls back to the existing map
  when its candidate fails the stated residual test.
- Any iteration-order change improves a named benchmark under an identical
  setup, or is rejected. No residual or physical acceptance tolerance is
  relaxed to make the comparison pass.

### Defer

Do not add another acceleration layer until the current Anderson paths and
maturity gates have been measured. Keep a radiation/population Schur solver as
a research prototype until its operator blocks and true residual are tested.

## Phase 4 - Physical benchmark suite and candidate architectures

### Benchmark matrix

Start with small reproducible models, then expand to full CMFGEN models:

| Case | Main question |
|---|---|
| Monotonic velocity | Does the disabled/new path preserve stock behavior? |
| Mild non-monotonic profile | Do basic sign-aware paths work? |
| One strong compression zone | Do the upwind and moment paths remain stable? |
| Multiple velocity reversals | Are direction changes and boundaries handled? |
| Interacting or shocked profile | Does field/population coupling converge? |
| Difficult ionization structure | Which population blocks limit convergence? |

### Work

- Keep inputs and initial states fixed for baseline and candidate comparisons.
- Record source revision, compiler and libraries, VADAT and model inputs,
  velocity and frequency grids, solver flags, wall time, memory, and exit
  reason.
- Record formal sweeps, refresh passes, Newton steps, GMRES iterations and
  restarts, true residuals, factor/apply cost, rejected steps, and safeguard
  counts.
- Compare radiation fields, populations, temperature, conservation measures,
  and spectra. A spectrum alone is not an acceptance test.
- Prototype physical continuation only through a documented model parameter
  or a mathematically defined homotopy. Verify that the final stage solves the
  requested physical model.
- Test deflation/recycling against identified slow modes and equal JVP or
  residual-evaluation budgets. Consider nonlinear multigrid only after a
  repeatable benchmark demonstrates a scale or cost limitation.

### Acceptance

- Each benchmark has a stored baseline and a written pass/fail rule before
  candidate runs begin.
- All claimed convergences pass the production physical residual and field
  consistency checks. Iteration caps and safeguard-dominated runs are reported
  as unresolved.
- The monotonic case preserves stock behavior with new flags disabled.
- Any claimed improvement repeats on the target case and does not conceal
  changes in conservation, observables, or safeguard activation.
- Report improvements and regressions, including wall time and memory. Do not
  promote an option to a production default until the full benchmark matrix
  supports it.

## Decision order

1. Reproduce and classify the failure in Phase 0.
2. Resolve JFNK correctness and the known preconditioner test gap in Phase 1.
3. Audit the frequency interval and resonance treatments in Phase 2 before
   changing their numerical guards.
4. Validate population admission and existing acceleration in Phase 3.
5. Compare candidate solver architectures on the Phase 4 benchmark matrix.

Stop a branch of work when its acceptance test fails. Record the failure and
return to the source of the residual, state, discretization, or conditioning
problem before adding another safeguard.

## Literature

- Knoll & Keyes 2004, "Jacobian-Free Newton-Krylov Methods: A Survey of
  Approaches and Applications," *Journal of Computational Physics* 193, 357.
- Eisenstat & Walker 1996, "Choosing the Forcing Terms in an Inexact Newton
  Method," *SIAM Journal on Scientific Computing* 17, 16.
- Knop, Hauschildt & Baron 2009a, *Astronomy & Astrophysics* 496, 295.
- Knop, Hauschildt & Baron 2009b, *Astronomy & Astrophysics* 501, 813.
- Dessart, Audit & Hillier 2015, *Monthly Notices of the Royal Astronomical
  Society* 449, 4304.

These references identify the methods in scope. Check the cited paper and its
assumptions before using it to justify a discretization or solver change.
