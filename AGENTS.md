# AGENTS.md — Devin development guide for CMFGEN non-monotonic fork

## 1. Project mission and scientific context

This repository is a development fork of CMFGEN for radiative-transfer and NLTE
supernova-ejecta models with **non-monotonic radial velocity fields**, especially
interacting Type IIn supernovae. The motivating structure is fast ejecta interacting
with slow, dense circumstellar material; reverse shocks and compression zones can
produce velocity reversals and regions where
\`sigma = d ln(V) / d ln(r)\` is below -1.

The stock CMFGEN treatment assumes a monotonic velocity law in several transfer,
Sobolev, interpolation, and related paths. The purpose of this fork is to extend
those paths while preserving the upstream solution whenever new options are disabled.

The project includes several distinct, opt-in pieces of work:

- A sign-aware Knop comoving-frame formal solver and moment-equation pathway.
- Non-monotonic Sobolev line treatment based on the Dessart, Audit & Hillier method.
- A Jacobian-Free Newton–Krylov (JFNK) population solver using GMRES, preconditioning,
  and globalization safeguards.
- Population/trace-species convergence handling and supporting numerical safeguards.
- Supernova hydrodynamic velocity-profile input, shock-power deposition, and
  observed-frame diagnostics.

Treat the code as **scientific software under active development**, not as a
production-validated solver. A successful compile, a passing unit test, or a
plausible-looking spectrum is not by itself evidence of physical correctness.

## 2. First actions for every task

Before editing:

1. Read this file completely.
2. Read \`README.md\` and the relevant build/test documentation, especially
   \`unit_tests/README.md\`.
3. Inspect the current branch, working tree, recent commits, and relevant call
   graph. Do not assume that prior audit notes still describe the current code.
4. Search for all callers, declarations, VADAT readers, module variables, and
   tests associated with the code being changed. Do not make a one-file fix when
   the feature crosses multiple layers.
5. Identify the stock/reference path and determine what must remain unchanged when
   the new feature is disabled.
6. Record a concise plan and the tests needed before making a broad or risky change.

Do not infer behavior from names or comments alone. Verify it in executable code
and tests. Comments can be stale.

## 3. Repository map

Use the following as an entry map; confirm paths against the current tree.

| Path | Role |
|---|---|
| \`new_main/\` | Main CMFGEN program and VADAT/control-variable integration |
| \`new_main/mod_subs/\` | Control-variable modules, readers, consistency checks, JFNK solver modules and related numerical infrastructure |
| \`new_main/subs/\` | New-solver support, population-solve integration, non-monotonic velocity and supernova model routines |
| \`plane/\` | Comoving-frame formal/moment solver implementations, including non-monotonic Knop paths |
| \`subs/\`, \`newsubs/\` | Core upstream and newer shared numerical/physics routines |
| \`obs/\` | Observed-frame tools, including flux and line-profile diagnostics |
| \`unit_tests/\` | Python-driven tests that compile and execute production Fortran source |
| \`tests/\` | Lower-level Fortran numerical/JFNK harnesses and compile/run checks |
| \`docs/\` | Scientific, implementation, and audit notes |
| \`build.sh\`, \`build_debug.sh\` | Supported production-style and debug build entry points |
| \`Makefile\`, \`Makefile_definitions\` | Build graph, compiler flags, libraries, and platform configuration |

Important current entry points observed in the repository:

- \`new_main/mod_subs/jfnk_solver_v1.f90\`: JFNK/GMRES/Hookstep numerical core.
- \`new_main/mod_subs/control_variable_mod.f\`: central control-variable declarations.
- \`new_main/mod_subs/rd_control_variables.f\`: VADAT input parsing and key handling.
- \`new_main/mod_subs/check_param_consistency.f\`: cross-parameter validation.
- \`plane/fg_j_cmf_nonmono_v1.f\`: non-monotonic Knop formal solver path.
- \`plane/mom_j_cmf_nonmono_v1.f\`: non-monotonic moment-equation J pathway (verify exact current path/callers before editing).
- \`unit_tests/README.md\`: inventory of currently documented production-source regression tests.

The README identifies the upstream baseline as \`cur_cmf_18jun25\` from
\`pitt-phyast/cmfsource\`. Preserve attribution and distinguish original CMFGEN
code from new fork-specific work.

## 4. Core engineering principles

### 4.1 Preserve upstream behavior

- New behavior must remain opt-in unless the task explicitly changes that policy.
- With all new flags disabled, preserve the legacy/upstream numerical path.
- For changes to existing routines, add or retain regression comparisons against
  the pristine reference where practical.
- Do not refactor large upstream routines merely for style. Prefer narrow, reviewable
  changes with explicit scientific or numerical justification.
- Do not silently alter input semantics, defaults, output formats, or physical
  assumptions.

### 4.2 Scientific correctness outranks convergence at any cost

Never make a model appear to converge by weakening or bypassing the physics or
acceptance tests. In particular:

- Do not replace non-finite values with arbitrary finite values to conceal failures.
- Do not accept a Newton/GMRES step merely because an iteration or restart cap was
  reached.
- Do not declare convergence solely from a small surrogate/Hessenberg residual.
- Do not relax population, radiative-equilibrium, or transfer tolerances without
  documenting the physical and numerical consequences.
- Do not clamp the physical velocity gradient globally to force a monotonic result.
- Distinguish numerical regularization from a change to the physical equations.
- Keep counters and diagnostics for safeguards that materially modify a step.

A robust solver must report why it stopped, whether the true residual passed, and
which safeguards fired.

### 4.3 Fortran safety and legacy conventions

- Respect the source form of each file: much of the code is fixed-form Fortran;
  selected modules and solver files use free-form Fortran.
- Check line-length and continuation rules for the relevant compiler flags.
- Preserve established kind conventions (e.g. \`LDP\`) and array layout.
- Be explicit about array bounds, leading dimensions, and workspace sizes.
- Avoid uninitialized variables, implicit interfaces, unsafe aliasing, and hidden
  dependence on global state.
- Check every allocation, deallocation, and dimension-change path.
- Ensure all error exits restore a valid state and release or preserve resources
  consistently.
- Avoid changing shared module state in helper routines unless that is part of the
  established interface and is documented.
- Do not introduce compiler-specific functionality without checking the supported
  Intel and GNU configurations.

### 4.4 Make changes maintainable

- Keep patches focused and avoid unrelated formatting churn.
- Explain non-obvious numerical choices in comments, including equations, units,
  sign conventions, and references.
- Keep VADAT key names, variable names, defaults, consistency checks, and tests in
  sync.
- Add tests alongside bug fixes. Prefer tests that fail on the old defect and pass
  with the fix.
- Do not duplicate a helper or numerical algorithm when an existing tested routine
  can be used safely.

## 5. Non-monotonic velocity physics and transfer

This is not merely a generic extension to accept a velocity array. The implementation
must preserve the physical meaning of the velocity gradient and its use in the
transfer equations.

When changing velocity-dependent code:

1. Trace the velocity profile from the input reader/model constructor through
   interpolation, grid remapping, gradient calculation, formal solution, moment
   equations, Sobolev treatment, and observed-frame processing.
2. Identify whether each routine expects radius and velocity to be ordered,
   positive, monotonic, or single-valued. State those assumptions explicitly.
3. Distinguish the raw/physical gradient from any legacy-clamped or numerically
   regularized gradient. Do not substitute one for the other without a derivation.
4. Verify the definition and sign convention of
   \`SIGMA = d ln(V)/d ln(r)\` in the exact routine being modified. Check
   compression zones, plateaus, zero/negative velocities, boundaries, and
   degenerate radial spacing.
5. Verify frequency-coupling direction and indexing for both signs of the
   relevant velocity-gradient terms. Do not assume that a monotonic-flow
   upwind stencil remains valid after a sign change.
6. Check that interpolation uses a valid monotonic independent coordinate and that
   the routine's dimensions/strides are correct.
7. Keep the legacy monotonic path bitwise or numerically equivalent within a
   justified tolerance when the new mode is off.
8. Test both monotonic and non-monotonic profiles, including multiple reversals,
   sharp compression zones, weak gradients near zero, boundaries, and pathological
   but explicitly rejected input.

For a physical input profile, do not silently reorder samples unless the associated
radius/velocity/thermodynamic quantities are reordered together and the scientific
meaning of the operation is demonstrated. Radius should remain a valid increasing
coordinate wherever the equations require it.

### 5.1 Knop formal and moment-equation paths

Treat the formal solver and moment-equation pathway as separate but coupled
implementations. Audit both whenever shared frequency coupling, SIGMA, source
functions, boundary conditions, or control flags change.

For \`USE_LAM_ES = F\`, verify that the non-monotonic moment-equation pathway is
actually selected and that the resulting J update is consistent with the
formal-solver/Knop treatment. Check parameter ordering: defaults and derived
flags must not depend on whether a VADAT key appears before or after another key.

For \`USE_LAM_ES = T\`, verify the intended lambda/formal-solution pathway and
ensure that the moment-equation path is not accidentally activated.

When changing \`DO_NONMONO_KNOP\`, \`DO_NONMONO_SOB\`, or related gates:

- Trace parsing, default assignment, alias handling, mutual exclusion, prerequisites,
  consistency checks, and final dispatch.
- Check all allowed combinations, including both flags off and any documented
  combined mode.
- Ensure unsupported combinations fail early with actionable diagnostics.
- Confirm the stock pathway remains available and unchanged.

### 5.2 Sobolev and line treatment

- Use the physically intended, unclamped gradient where the non-monotonic method
  requires it; do not conflate it with a legacy positive-gradient guard.
- Verify the optical-depth denominator and its floor against the cited method.
  A floor must prevent singular numerical behavior without masking a physically
  distinct regime.
- Check sign changes, \`sigma < -1\`, values near singular denominators, weak lines,
  masing/negative opacity, and boundaries.
- Audit every caller of shared Sobolev routines. A fix for one path must not
  silently change unrelated monotonic or non-SN use cases.

### 5.3 Hydrodynamic input and shock deposition

For SN hydrodynamic profiles:

- Validate file format, units, finite values, increasing radius, positive radius,
  velocity behavior, and all required columns before using the profile.
- Check endpoint handling and extrapolation explicitly.
- Test homologous expansion, beta-law/monotone profiles, and velocity reversals.
- Ensure derived SIGMA follows the analytic expectation for manufactured profiles.
- For shock deposition, verify the integrated deposited power against the requested
  power and test the spatial profile, normalization, clumping behavior, and invalid
  coordinates. Do not assume velocity uniquely locates a shell when V is
  non-monotonic; use the intended radius-based path where appropriate.

## 6. JFNK, GMRES, and preconditioner requirements

The JFNK path is a coupled scientific solver, not an isolated numerical library.
A change must be traced from VADAT controls through the CMFGEN interface, residual
construction, Jacobian-vector product, preconditioner assembly/application,
GMRES, globalization, state update, and convergence decision.

### 6.1 Residual and frozen-J semantics

- Write down the exact residual vector and variable ordering before modifying it.
- Verify row/column indexing, masks, scaling, units, and weighting.
- Confirm exactly which variables are solved simultaneously and which are frozen
  or lagged (including radiation field J, temperature, electron density, and
  populations).
- A frozen-J linearization is not the full derivative of a fully coupled
  radiation/population system. Document precisely which dependencies are retained
  and omitted.
- Verify that the baseline residual and every finite-difference perturbed residual
  use consistent state, masks, and scaling.
- Ensure trial residual evaluation cannot leave CMFGEN globals, populations,
  radiation fields, or caches mutated. Restore the accepted state after rejected
  trials.
- The JVP must be tested against analytic or finite-difference reference problems
  over a range of perturbation sizes. Watch for cancellation, noise, clipping,
  and branch changes.

### 6.2 GMRES correctness

- Verify Arnoldi orthogonalization, Givens rotations, Hessenberg indexing,
  restart handling, breakdown, happy breakdown, and zero RHS.
- Test the true linear residual \`b - A x\`; do not trust only the small least-squares
  residual inferred from the Hessenberg system.
- For finite-difference JVPs, recognize that the computed action may not behave
  like a perfectly fixed linear operator. Explicitly measure the true residual
  with the same production operator where possible.
- Ensure convergence tolerances use the correct norm under left or right
  preconditioning. The meaning of the forcing term and trust-region radius must
  remain consistent with the selected preconditioning side.
- Bound total iterations and restarts; expose whether limits were reached.
- Treat loss of orthogonality, stagnation, invariant Krylov subspaces, and
  breakdowns explicitly rather than returning a nominally converged step.

### 6.3 Preconditioner correctness

The current JFNK documentation describes a fixed-J BA preconditioner with
depth-block factorization/equilibration and local adaptive ridge. Do not assume
that description remains accurate after code changes: inspect the live factorization
and application routines.

Audit:

- Exact matrix blocks and their correspondence to the residual ordering.
- Row/column equilibration and the inverse transformations on solution vectors.
- LU factorization, pivoting, singularity detection, conditioning estimates, and
  solve status.
- Ridge/regularization sign and units, including rows with negative diagonal terms.
- Treatment of zero/tiny rows, masked variables, regularized rows, and their columns.
- Refresh policy: determine when the matrix is stale, what triggers rebuild, and
  whether failed refreshes are handled safely.
- Consistency between factorization and apply routines.
- Preconditioner failure must be visible and must not silently become the identity
  unless the algorithm explicitly specifies and diagnoses that fallback.

For every regularization or safeguard, report how often it fires and test its
effect on both the linear residual and the outer nonlinear convergence.

### 6.4 Newton globalization and acceptance

- Verify the sign of the Newton correction and the actual update convention.
- Check line search, non-monotone reference window, Cauchy/residual fallback,
  Hookstep/trust-region logic, step scaling, and rejection paths as one state
  machine.
- Every accepted step must satisfy its documented acceptance rule using a
  trustworthy residual evaluation.
- Rejected steps must not contaminate the accepted state or solver history.
- Enforce population positivity and physical bounds without silently changing
  the nonlinear problem.
- Distinguish convergence, stagnation, iteration cap, GMRES failure, non-finite
  residual, invalid preconditioner, and user abort in exit status and diagnostics.
- Do not use a relaxed or masked residual to claim convergence unless the final
  full/physical residual is also checked according to the intended convergence
  contract.

### 6.5 Coupled and legacy modes

Audit the complete mode matrix, not just the default configuration:

- JFNK disabled: upstream population solution unchanged.
- JFNK enabled with each supported formal/moment solver combination.
- Coupled versus lagged/frozen-J operation, as implemented and documented.
- Frozen-J and radiation-field refresh timing, including restart/checkpoint behavior.
- Temperature/electron-density fixed versus free paths where supported.
- Masked/trace species and Stage-B polish paths.
- Preconditioner variants and left/right GMRES modes if exposed.
- Unsupported flag combinations must be rejected before solver state is allocated.

Do not state that populations and J are solved simultaneously unless the actual
production call graph demonstrates this for the selected mode.

## 7. VADAT controls and dispatch

A VADAT key is an API. Changing it requires a complete audit.

For each new or modified control:

1. Add or verify its declaration and safe default.
2. Add or verify parsing, case normalization, duplicate-key behavior, and diagnostics.
3. Verify assignment order does not cause derived defaults to depend on input-key
   order.
4. Add prerequisites, mutual exclusions, aliases/deprecations, and consistency
   checks as needed.
5. Verify the final runtime dispatch uses the intended value.
6. Add parser and consistency tests for valid, invalid, omitted, reordered, and
   conflicting combinations.
7. Document the key, default, and interaction with other controls.

Never change a default silently. Preserve backward compatibility where possible.
Avoid reading an uninitialized logical or numeric control before VADAT parsing.

## 8. Testing and build workflow

### 8.1 Supported commands

From the repository root:

\`\`\`bash
# Full production-style build plus compile-time and Python regression checks
./build.sh

# Debug CMFGEN build
./build_debug.sh

# All production-source regression tests
python3 unit_tests/run_all_tests.py

# Single Python regression test
python3 unit_tests/python/test_<name>.py

# Lower-level JFNK/solver compile and run harness
tests/f90/build_and_run.sh
\`\`\`

Use the available compiler and libraries on the host. The documented build setup
uses Intel \`ifort\`/MKL on the named HPC platforms, with GNU fallbacks in selected
configurations. Check \`Makefile_definitions\` before changing flags or assuming
compiler availability.

The full build creates the output directories and checks required executables.
It also runs the lower-level JFNK harness and all Python production-source tests.
Do not replace this with a compile-only check.

### 8.2 Unit-test design

The \`unit_tests/\` suite compiles the **actual production Fortran routines** with
small Fortran CLI drivers. Keep that property: do not reimplement the numerical
routine in Python and test only the reimplementation.

When practical:

- Compare flags-off behavior against the pristine reference tree.
- Use manufactured analytic profiles/systems with known solutions.
- Include both expected-success and expected-failure cases.
- Test finite values and physical invariants, not just process exit status.
- Test boundary and degenerate dimensions.
- Assert diagnostic text and exit code for rejected inputs.
- Include regression cases reproducing the original defect.
- Keep generated plots diagnostic only; numeric assertions are authoritative.

Consult \`unit_tests/README.md\` for the live inventory before adding or renaming tests.
Do not claim a test exists or passes unless you have run it on the current working
tree and observed its result.

### 8.3 Numerical validation ladder

Choose tests proportionate to the change:

1. Syntax/compile check for touched source.
2. Focused production-source unit test.
3. Relevant neighboring regression tests.
4. Full \`unit_tests/run_all_tests.py\`.
5. \`tests/f90/build_and_run.sh\` for JFNK/interface changes.
6. Full \`./build.sh\` where the environment supports it.
7. Controlled runtime model tests spanning monotonic and non-monotonic profiles.
8. Compare convergence histories, physical residuals, radiation fields, populations,
   and spectra against known baseline/reference cases.

Do not run a large production model before focused tests pass. Do not present an
untested build or runtime result as successful.

### 8.4 Runtime evidence

For any model-run validation, retain enough information to reproduce and diagnose it:

- Commit SHA and build/compiler configuration.
- Input VADAT and model/hydrodynamic input identifiers.
- Active solver flags and defaults.
- OUTGEN and relevant diagnostic logs.
- Iteration, residual, population-change, radiation-field, and temperature histories.
- GMRES iteration/restart counts, true residual, preconditioner refresh/failure counts,
  regularization/safeguard counters, and exit reason.
- Whether the model converged under the full intended criteria.
- Any differences from the reference or previous commit.

A spectrum alone is insufficient validation.

## 9. Working with tests, data, and generated files

- Do not commit compiled objects, libraries, executables, module files, temporary
  build products, large runtime dumps, or generated diagnostics unless explicitly
  requested.
- Do not overwrite user-provided model inputs, reference outputs, or production data.
- Keep test fixtures small and deterministic.
- Respect \`.gitignore\` and inspect \`git status\` before and after testing.
- Do not delete failing tests or weaken assertions to obtain a green suite.
- If a test exposes a genuine mismatch between the expected physics and a legacy
  reference, explain the distinction and add a scientifically justified expected
  result rather than masking it.

## 10. Git and change-management rules

- Work on the branch requested by the user; otherwise create/use a dedicated feature
  branch. Never force-push or rewrite shared history.
- Inspect \`git status\` before editing. Preserve unrelated user changes.
- Keep commits focused and use descriptive messages.
- Do not merge or push changes unless the user explicitly asks.
- Do not modify upstream/reference directories to make comparisons pass.
- For broad solver changes, split work into logically reviewable commits where
  practical.
- Before finishing, inspect the full diff for accidental changes, whitespace damage,
  generated files, and modifications outside task scope.

## 11. Devin task execution protocol

For every task:

### A. Understand
- Restate the requested behavior and acceptance criteria.
- Identify affected physics, routines, modules, controls, and call paths.
- Find related tests and historical implementation notes.
- List assumptions and unknowns that need verification.

### B. Audit before patching
- Reproduce the issue or establish a failing regression test.
- Trace the path end-to-end, including callers and state mutation.
- Check the legacy path and all relevant flag combinations.
- Identify whether the defect is physical-model, discretization, indexing,
  state-management, conditioning, convergence-control, or build/test related.

### C. Implement
- Make the smallest complete change that fixes the root cause.
- Add a regression test that would have caught the defect.
- Update VADAT parsing, defaults, consistency checks, diagnostics, and documentation
  together when a control or behavior changes.
- Explain numerical choices and safeguards in code comments.

### D. Verify
- Run focused tests first, then relevant suites and builds.
- Inspect test output, exit codes, warnings, and runtime diagnostics.
- Check that disabled/new modes retain the required behavior.
- For solver changes, independently verify the final physical residual and state.
- Do not suppress or reinterpret failures.

### E. Report
Finish with a concise but technically complete report containing:

1. **What changed** — files and behavior.
2. **Why** — root cause and scientific/numerical rationale.
3. **Tests run** — exact commands and pass/fail outcomes.
4. **Unrun tests** — what could not be executed and why.
5. **Remaining risks** — unresolved assumptions, numerical caveats, or missing runtime
   validation.
6. **Git state** — branch and commit(s), without claiming a push/merge unless done.

Separate verified facts from hypotheses. Cite file paths and routine names; use line
numbers when available.

## 12. Do not do these things

- Do not perform broad refactors during a targeted numerical fix.
- Do not assume a source comment is correct without checking its implementation.
- Do not infer runtime convergence from compilation or unit tests.
- Do not call a solver robust merely because it has clipping, damping, floors,
  restarts, or fallback paths.
- Do not let a safety cap masquerade as convergence.
- Do not change physics, units, signs, or default flags without documenting and
  testing the consequences.
- Do not silently change the reference solution or weaken comparisons.
- Do not leave a half-wired feature (e.g. declaration without parser, parser without
  dispatch, dispatch without tests).
- Do not claim that the entire non-monotonic treatment is validated based on a
  single manufactured ray or one successful model.

## 13. Scientific references

Use the repository README and the relevant source comments/docs for exact citations.
Key method references currently listed in \`README.md\` include:

- Hillier & Miller (1998), CMFGEN, ApJ 496, 407.
- Knop, Hauschildt & Baron (2009), A&A 496, 295 and A&A 501, 813.
- Dessart, Audit & Hillier (2015), MNRAS 449, 4304.
- Feldmeier & Nikutta (2006), A&A 446, 661.
- Knoll & Keyes (2004), J. Comput. Phys. 193, 357.
- Arramy, de la Cruz Rodriguez & Leenaarts (2024), A&A 690, A12.
- Eisenstat & Walker (1996), SIAM J. Sci. Comput. 17, 16.
- Kelley (2003), *Solving Nonlinear Equations with Newton's Method*.

Do not add a literature claim or alter an equation based on an unverified citation.
Check the original paper and ensure the implemented discretization matches the
assumptions and conventions being used.
