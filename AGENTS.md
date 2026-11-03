# AGENTS.md — Devin development guide for CMFGEN non-monotonic fork

## 1. Project mission and scientific context

This repository is a development fork of CMFGEN for radiative-transfer and NLTE
supernova-ejecta models with **non-monotonic radial velocity fields**, especially
interacting Type IIn supernovae. Fast ejecta interacting with slow, dense
circumstellar material can produce reverse shocks, compression zones, velocity
reversals, and regions where `sigma = d ln(V) / d ln(r)` is below -1.

The purpose of this fork is to extend CMFGEN's transfer, Sobolev, interpolation,
and related paths while preserving upstream behaviour whenever new options are
disabled.

The project includes:
- A sign-aware Knop comoving-frame formal solver and moment-equation pathway.
- Non-monotonic Sobolev line treatment based on Dessart, Audit & Hillier.
- A Jacobian-Free Newton–Krylov (JFNK) population solver using GMRES,
  preconditioning, and globalization safeguards.
- Population/trace-species convergence handling and numerical safeguards.
- Supernova hydrodynamic velocity-profile input, shock-power deposition, and
  observed-frame diagnostics.

Treat this as scientific software under active development, not as a
production-validated solver. Compilation, unit tests, or a plausible spectrum
alone do not establish physical correctness.

## 2. First actions for every task

Before editing:
1. Read this file, `README.md`, and relevant build/test documentation,
   especially `unit_tests/README.md`.
2. Inspect the current branch, working tree, recent commits, and relevant call
   graph. Do not assume prior audit notes still describe the current code.
3. Search all callers, declarations, VADAT readers, module variables, and tests
   associated with the code being changed.
4. Identify the stock/reference path and what must remain unchanged when new
   features are disabled.
5. Record a concise plan and required tests before broad or risky changes.

Do not infer behaviour from names or comments alone. Verify it in executable
code and tests; comments can be stale.

## 3. Repository map

Confirm paths against the current tree. Important areas include:

| Path | Role |
|---|---|
| `new_main/` | Main CMFGEN program and VADAT/control-variable integration |
| `new_main/mod_subs/` | Control variables, readers, consistency checks, JFNK modules |
| `new_main/subs/` | Solver integration, non-monotonic velocity and SN routines |
| `plane/` | Comoving-frame formal/moment solvers, including non-monotonic Knop paths |
| `subs/`, `newsubs/` | Core upstream and newer shared numerical/physics routines |
| `obs/` | Observed-frame tools, flux and line-profile diagnostics |
| `unit_tests/` | Python-driven tests compiling and executing production Fortran |
| `tests/` | Lower-level Fortran numerical/JFNK harnesses |
| `docs/` | Scientific, implementation, and audit notes |
| `build.sh`, `build_debug.sh` | Production-style and debug build entry points |
| `Makefile`, `Makefile_definitions` | Build graph, compiler flags, libraries, platform setup |

Important entry points observed in the repository:
- `new_main/mod_subs/jfnk_solver_v1.f90`: JFNK/GMRES/Hookstep numerical core.
- `new_main/mod_subs/control_variable_mod.f`: central control declarations.
- `new_main/mod_subs/rd_control_variables.f`: VADAT parsing.
- `new_main/mod_subs/check_param_consistency.f`: cross-parameter validation.
- `plane/fg_j_cmf_nonmono_v1.f`: non-monotonic Knop formal solver path.
- `plane/mom_j_cmf_nonmono_v1.f`: non-monotonic moment-equation J pathway;
  verify exact current path and callers before editing.
- `unit_tests/README.md`: documented regression-test inventory.

The README identifies the upstream baseline as `cur_cmf_18jun25` from
`pitt-phyast/cmfsource`. Preserve attribution and distinguish original CMFGEN
code from fork-specific work.

## 4. Core engineering principles

### 4.1 Preserve upstream behaviour
- New behaviour must remain opt-in unless explicitly requested otherwise.
- With all new flags disabled, preserve the legacy/upstream numerical path.
- Add or retain regression comparisons against the pristine reference where
  practical.
- Avoid broad refactors of upstream routines. Prefer narrow, reviewable patches.
- Do not silently alter input semantics, defaults, output formats, or physics.

### 4.2 Scientific correctness outranks convergence at any cost
Never make a model appear to converge by weakening or bypassing the physics or
acceptance tests.

- Do not replace non-finite values with arbitrary finite values to conceal
  failures.
- Do not accept a Newton/GMRES step merely because an iteration/restart cap was
  reached.
- Do not declare convergence solely from a surrogate/Hessenberg residual.
- Do not relax population, radiative-equilibrium, or transfer tolerances without
  documenting the physical and numerical consequences.
- Do not clamp the physical velocity gradient globally to force monotonicity.
- Distinguish numerical regularization from a change to the physical equations.
- Keep diagnostics for safeguards that materially modify a step.

A robust solver must report why it stopped, whether the true residual passed,
and which safeguards fired.

### 4.3 Fortran safety and legacy conventions
- Respect each file's source form: much of CMFGEN is fixed-form Fortran;
  selected modules and solver files use free-form Fortran.
- Check line-length and continuation rules for the active compiler flags.
- Preserve kind conventions (e.g. `LDP`) and array layout.
- Be explicit about bounds, leading dimensions, and workspace sizes.
- Avoid uninitialized variables, implicit interfaces, unsafe aliasing, and
  hidden dependence on global state.
- Check every allocation, deallocation, and dimension-change path.
- Ensure error exits leave state and resources consistent.
- Do not introduce compiler-specific functionality without checking supported
  Intel and GNU configurations.

### 4.4 Maintainability
- Keep patches focused; avoid unrelated formatting churn.
- Explain non-obvious numerical choices, equations, units, sign conventions,
  and references in comments.
- Keep VADAT keys, declarations, defaults, consistency checks, and tests in sync.
- Add regression tests alongside bug fixes.
- Reuse existing tested helpers where safe.

## 5. Critical policy: no artificial clamps that hide physical problems

**Default rule: do not introduce a clamp, floor, ceiling, clipping operation,
limiter, saturation, positivity repair, gradient replacement, or other
artificial bound merely to make a calculation stable, reduce warnings, improve
convergence, or produce a plausible output.** Such a change can conceal a
physical regime the code is specifically intended to resolve.

Before adding or retaining any bound or limiter, Devin must:

1. **Identify the underlying quantity and its physical meaning.** State whether
   it is a physical variable, a derived physical coefficient, a discretization
   quantity, a solver variable, or a purely numerical intermediate.
2. **Locate the source of the problematic value.** Trace it upstream through
   input data, interpolation, derivatives, discretization, boundary conditions,
   residual construction, and solver state. Do not treat the location where a
   NaN, divergence, negative value, or extreme coefficient first appears as
   necessarily the root cause.
3. **Establish whether the value is physically allowed.** Check the governing
   equations, units, sign convention, limiting behaviour, and relevant
   literature. A surprising value is not automatically unphysical, especially
   in shocks, compression regions, non-monotonic flows, or optically complex
   regimes.
4. **Explain what the proposed bound changes mathematically.** Identify the
   exact expression being altered, the range affected, and the consequences
   for conservation, differentiability, monotonicity, frequency coupling,
   optical depth, energy balance, and the solution.
5. **Demonstrate that the bound is not masking a defect.** Reproduce the issue
   without the clamp where possible; test a manufactured case with known
   behaviour; and compare the unclamped and bounded results.
6. **Prefer fixing the root cause.** Correct bad indexing, invalid interpolation,
   inconsistent units, wrong derivative/sign, defective boundary treatment,
   ill-conditioned equations, stale state, or an inappropriate algorithm
   before considering an artificial bound.
7. **If a bound is genuinely required for numerical safety, make it explicit.**
   Document its mathematical justification, physical scope, threshold, units,
   and expected impact. Make it configurable or disableable where appropriate.
   Emit diagnostics and count every activation. Do not silently apply it.
8. **Test both sides of the threshold.** Include cases below, at, and above the
   threshold, plus cases where the underlying physical quantity legitimately
   crosses zero or changes sign.
9. **Quantify the impact.** Report how often it activates and how much it
   changes the solution, residual, conservation balance, and relevant physical
   observables.
10. **Never use a clamp as a convergence criterion.** A run that converges only
    when an undocumented or physically unjustified clamp is active must be
    reported as unresolved, not successful.

### 5.1 Specific prohibitions and cautions

- Never replace the physical `SIGMA = d ln(V) / d ln(r)` with a positive,
  clipped, floored, or monotonicized surrogate in a pathway intended to solve
  non-monotonic flow, unless the derivation explicitly requires a distinct
  regularized quantity and the physical SIGMA remains available.
- Do not clip negative velocity gradients, negative opacity, negative
  corrections, source terms, rates, or energy terms simply because their sign
  causes numerical difficulty. First establish whether the sign is physically
  meaningful.
- Do not floor denominators or optical depths without analysing the limiting
  equation and the effect on the physical regime. A denominator approaching
  zero may indicate a real singular limit, a change of regime, or a discretization
  defect; a floor must not decide between these possibilities silently.
- Do not impose positivity by clipping populations, intensities, rates, or
  corrections without tracing the origin and checking the governing constraints.
  If a positivity-preserving variable transformation or constrained solve is
  appropriate, justify it and verify that it solves the intended equations.
- Do not cap Newton steps, GMRES corrections, radiation fields, temperatures,
  or population changes solely to suppress divergence. Solver globalization
  must use a stated acceptance rule and must preserve the meaning of the
  nonlinear problem.
- Do not use arbitrary epsilon substitutions, `MAX/MIN` guards, absolute values,
  fallback constants, or replacement values to turn invalid states into
  apparently valid states.
- Do not add a clamp to silence diagnostic output. Fix the source of the warning
  or improve the diagnostic so it distinguishes expected physical behaviour
  from a genuine failure.
- Do not let a regularized preconditioner silently redefine the physical
  residual or accepted solution. A preconditioner may be regularized as a
  numerical aid, but the original system and its true residual must still be
  evaluated and reported.

### 5.2 Required review of existing clamps

When auditing existing code, actively search for:
- `MIN`, `MAX`, `AMAX1`, `AMIN1`, `ABS`, `SIGN`, `MAXVAL`, `MINVAL`,
  `MERGE`, and hand-written `IF` guards around physical quantities.
- Epsilon additions, denominator floors, gradient floors, optical-depth
  cutoffs, opacity/rate clipping, positivity repairs, and saturation logic.
- Branches that substitute a fallback value after NaN/Inf, singularity,
  negative values, failed factorization, or non-convergence.
- Bounds on corrections, line-search steps, trust regions, and residuals.

Do not assume every `MIN/MAX` or floor is wrong: distinguish physical model
limits, mathematically justified domain restrictions, discretization safeguards,
and solver-only regularization. For each suspicious case, record:
- source file and routine;
- quantity and units;
- original equation and modified expression;
- activation condition and frequency;
- whether the value is physically admissible;
- impact on conservation and solution observables;
- recommended action: retain with justification, expose/configure, replace,
  or remove.

Do not mechanically delete all guards. Removing a guard without understanding
the failure can be as damaging as adding one.

## 6. Non-monotonic velocity physics and transfer

This is not merely a generic extension to accept a velocity array. Preserve the
physical meaning of the velocity gradient and its use in transfer equations.

When changing velocity-dependent code:
1. Trace the velocity profile from input/model construction through interpolation,
   grid remapping, gradient calculation, formal solution, moment equations,
   Sobolev treatment, and observed-frame processing.
2. Identify whether each routine assumes radius/velocity are ordered, positive,
   monotonic, or single-valued.
3. Distinguish raw/physical gradients from legacy-clamped or regularized gradients.
4. Verify the definition and sign convention of `SIGMA = d ln(V) / d ln(r)`.
   Check compression zones, plateaus, zero/negative velocities, boundaries,
   and degenerate radial spacing.
5. Verify frequency-coupling direction and indexing for both signs of the
   relevant velocity-gradient terms.
6. Check interpolation coordinate assumptions, dimensions, and strides.
7. Preserve the legacy monotonic path when the new mode is off.
8. Test monotonic and non-monotonic profiles, multiple reversals, sharp
   compression zones, weak gradients near zero, boundaries, and rejected input.

Do not silently reorder physical samples unless all associated quantities are
reordered consistently and the scientific meaning is demonstrated. Radius must
remain increasing wherever required by the equations.

### 6.1 Knop formal and moment-equation paths

Treat formal and moment-equation pathways as separate but coupled implementations.
Audit both when shared frequency coupling, SIGMA, source functions, boundary
conditions, or control flags change.

For `USE_LAM_ES = F`, verify that the non-monotonic moment-equation pathway is
selected and that the J update is consistent with the intended Knop treatment.
Parameter defaults and derived flags must not depend on VADAT key order.

For `USE_LAM_ES = T`, verify the intended lambda/formal-solution pathway and
ensure the moment-equation path is not accidentally activated.

When changing `DO_NONMONO_KNOP`, `DO_NONMONO_SOB`, or related gates, trace
parsing, defaults, aliases, mutual exclusions, prerequisites, consistency
checks, and final dispatch. Test all supported combinations and preserve the
stock pathway.

### 6.2 Sobolev and line treatment
- Use the physically intended gradient where the non-monotonic method requires it.
- Verify optical-depth denominators and floors against the cited method.
- Test `sigma < -1`, sign changes, near-singular denominators, weak lines,
  negative opacity/masing, and boundaries.
- Audit all callers of shared Sobolev routines to avoid unintended changes to
  monotonic or non-SN cases.

### 6.3 Hydrodynamic input and shock deposition
Validate file format, units, finite values, increasing radius, positive radius,
velocity behaviour, required columns, endpoint handling, and extrapolation.
Test homologous expansion, monotone profiles, and velocity reversals. Verify
derived SIGMA against analytic manufactured profiles. For shock deposition,
check integrated deposited power, spatial profile, normalization, clumping,
and invalid coordinates. Do not assume velocity uniquely locates a shell in
a non-monotonic flow.

## 7. JFNK, GMRES, and preconditioner requirements

Trace changes from VADAT controls through the CMFGEN interface, residual,
JVP, preconditioner assembly/application, GMRES, globalization, state update,
and convergence decision.

### 7.1 Residual and frozen-J semantics
- Write down the exact residual vector and variable ordering.
- Verify row/column indexing, masks, scaling, units, and weighting.
- Establish exactly which variables are solved simultaneously and which are
  frozen/lagged (including J, temperature, electron density, and populations).
- A frozen-J linearization is not the full derivative of a fully coupled
  radiation/population system. Document retained and omitted dependencies.
- Ensure baseline and perturbed residual evaluations use consistent state,
  masks, and scaling.
- Ensure trial evaluations cannot leave globals, populations, radiation fields,
  or caches mutated. Restore accepted state after rejected trials.
- Test JVPs against analytic/reference problems over a range of perturbation
  sizes; check cancellation, noise, clipping, and branch changes.

### 7.2 GMRES correctness
- Verify Arnoldi orthogonalization, Givens rotations, Hessenberg indexing,
  restart handling, breakdown, happy breakdown, and zero RHS.
- Test the true residual `b - A x`; do not trust only the Hessenberg estimate.
- Finite-difference JVPs may not act as a perfectly fixed linear operator;
  explicitly measure the true residual with the production operator.
- Ensure tolerances use the correct norm under left/right preconditioning.
- Bound iterations and restarts; expose limit exits.
- Handle loss of orthogonality, stagnation, invariant Krylov subspaces, and
  breakdown explicitly.

### 7.3 Preconditioner correctness
The current documentation describes a fixed-J BA preconditioner with
depth-block factorization/equilibration and local adaptive ridge. Verify the
live implementation rather than relying on comments.

Audit matrix blocks and residual ordering; equilibration and inverse transforms;
LU/pivoting/conditioning/solve status; ridge sign and units; tiny/zero/masked
rows and columns; refresh policy; and consistency between factorization and
application. Preconditioner failure must be visible and must not silently become
the identity unless explicitly specified and diagnosed.

For every regularization/safeguard, report activation counts and test its effect
on both the linear residual and outer nonlinear convergence.

### 7.4 Newton globalization and acceptance
Verify correction sign and update convention. Review line search, non-monotone
reference window, Cauchy fallback, Hookstep/trust region, step scaling, and
rejection paths as one state machine. Every accepted step must satisfy its
documented acceptance rule using a trustworthy residual. Rejected steps must
not contaminate accepted state/history. Enforce physical bounds without
silently changing the nonlinear problem. Distinguish convergence, stagnation,
iteration cap, GMRES failure, non-finite residual, invalid preconditioner,
and user abort.

### 7.5 Coupled and legacy modes
Audit all supported combinations:
- JFNK disabled: upstream population solution unchanged.
- JFNK with each supported formal/moment solver.
- Coupled versus lagged/frozen-J operation.
- Frozen-J and radiation-field refresh timing.
- Fixed/free temperature and electron-density paths where supported.
- Masked/trace species and Stage-B polish paths.
- Preconditioner variants and left/right GMRES if exposed.

Do not claim populations and J are solved simultaneously unless the production
call graph demonstrates this for the selected mode.

## 8. VADAT controls and dispatch

A VADAT key is an API. For every new/modified control:
1. Verify declaration and safe default.
2. Verify parsing, case normalization, duplicate-key behaviour, and diagnostics.
3. Ensure derived defaults do not depend on input-key order.
4. Add prerequisites, mutual exclusions, aliases/deprecations, and consistency
   checks as needed.
5. Verify runtime dispatch.
6. Add tests for valid, invalid, omitted, reordered, and conflicting combinations.
7. Document key, default, and interactions.

Never change a default silently or read an uninitialized control before parsing.

### 8.1 Keep the VADAT interface small: robust defaults, optional overrides

The new numerical code must not turn VADAT into a long list of low-level solver
tuning parameters. Every additional user-facing control increases configuration
complexity, creates more combinations to validate, and makes reproducibility and
support harder. Prefer a small, stable set of scientifically meaningful
high-level options.

**Default design principle:** optional features should work with internally
selected, tested, robust default values. Users should not need to tune numerical
internals just to obtain a valid calculation.

For every proposed VADAT parameter, Devin must:

1. Explain why the setting must be user-configurable rather than derived,
   detected, or selected internally.
2. Identify the scientific or operational decision it represents. Do not expose
   implementation details (iteration counters, tiny tolerances, damping factors,
   finite-difference perturbations, ridge values, restart lengths, clipping
   thresholds, or similar internals) as routine user controls without a clear
   use case.
3. Prefer one high-level switch for an optional capability, with robust internal
   defaults, rather than many interdependent switches controlling its machinery.
4. Derive dependent settings consistently from model dimensions, scales,
   conditioning estimates, or other relevant runtime information where this is
   mathematically justified and tested.
5. Keep advanced overrides optional. If expert tuning controls are genuinely
   needed for research or diagnostics, give them safe defaults, document their
   purpose and interactions, and make the ordinary workflow independent of them.
6. Do not silently select a value that changes the physical equations. Automatic
   defaults may control numerical algorithms, but must not hide physical
   assumptions or substitute for explicit scientific choices.
7. Validate the default configuration as a first-class supported mode. Tests
   must cover omitted keys, explicit defaults, and supported combinations.
   Omitted optional keys must never leave variables uninitialized.
8. Test that the default configuration behaves robustly across representative
   monotonic and non-monotonic models, including difficult conditioning and
   boundary cases. A default is not "robust" merely because one benchmark passes.
9. Document the selected defaults, their rationale, the conditions under which
   they may be inadequate, and how an expert can override them when appropriate.
10. Avoid adding a new VADAT key solely to work around a defect. Fix the defect
    or improve the internal algorithm instead.

Before introducing a new key, Devin should prepare a short parameter proposal:
- proposed key and default;
- the user-level decision it controls;
- why an internal default or automatic selection is insufficient;
- interactions with existing keys;
- tests supporting the default and override;
- whether it changes physics or only numerical strategy.

Reject redundant, overlapping, or rarely useful keys. Consolidate related
settings where that does not obscure their meaning. Keep the common VADAT
configuration concise and readable; put specialist diagnostics and development
experiments behind optional expert controls or dedicated test harnesses rather
than exposing every internal variable to all users.

A sensible hierarchy is:
- **Normal use:** optional feature switch; tested internal defaults; minimal
  user configuration.
- **Advanced use:** a small number of documented overrides for cases where
  model properties or research requirements justify them.
- **Development/testing:** detailed instrumentation and experimental tuning
  through test drivers, diagnostic builds, or developer-only mechanisms, not
  an expanding production VADAT interface.

Do not confuse fewer parameters with less transparency: print or record the
effective numerical settings in diagnostics/logs when they materially affect
the calculation, so runs remain reproducible even when users accept defaults.

## 9. Testing and build workflow

From repository root:

```bash
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
```

Check `Makefile_definitions` for available compilers/libraries. The documented
setup uses Intel ifort/MKL on named HPC platforms, with GNU fallbacks in
selected configurations.

The unit-test suite compiles actual production Fortran routines with small
Fortran CLI drivers. Do not test only Python reimplementations. Prefer
manufactured analytic systems, expected-success and expected-failure cases,
physical invariants, boundary/degenerate dimensions, and assertions on
diagnostics and exit codes. Numeric assertions are authoritative; plots are
diagnostic only.

Validation ladder:
1. Syntax/compile check for touched source.
2. Focused production-source unit test.
3. Relevant neighbouring regression tests.
4. Full `unit_tests/run_all_tests.py`.
5. `tests/f90/build_and_run.sh` for JFNK/interface changes.
6. Full `./build.sh` where supported.
7. Controlled runtime tests spanning monotonic and non-monotonic profiles.
8. Compare convergence histories, physical residuals, radiation fields,
   populations, and spectra against reference cases.

Do not run large production models before focused tests pass. Never report an
untested build/runtime as successful.

For runtime validation, retain commit SHA, compiler/build configuration,
VADAT/model inputs, active flags, OUTGEN/logs, iteration and physical residual
histories, population/radiation/temperature histories, GMRES statistics,
preconditioner/safeguard counters, exit reason, and reference comparisons.
A spectrum alone is insufficient evidence.

## 10. Tests, data, and generated files
- Do not commit objects, libraries, executables, module files, temporary build
  products, large runtime dumps, or generated diagnostics unless requested.
- Do not overwrite user inputs or reference outputs.
- Keep fixtures small and deterministic.
- Respect `.gitignore`; inspect `git status` before and after testing.
- Do not delete failing tests or weaken assertions to obtain a green suite.
- Explain scientifically justified differences from legacy reference results.

## 11. Git and change management
- Work on the requested branch; otherwise use a dedicated feature branch.
- Never force-push or rewrite shared history.
- Inspect `git status` and preserve unrelated user changes.
- Keep commits focused and descriptive.
- Do not merge or push unless explicitly requested.
- Do not alter upstream/reference directories to make comparisons pass.
- Review the full diff for accidental/generated changes before finishing.

## 12. Devin task execution protocol

### A. Understand
Restate requested behaviour and acceptance criteria. Identify affected physics,
routines, modules, controls, call paths, tests, assumptions, and unknowns.

### B. Audit before patching
Reproduce the issue or establish a failing regression. Trace the full path,
including callers and state mutation. Check legacy behaviour and relevant flag
combinations. Classify the issue: physics, discretization, indexing,
state-management, conditioning, convergence control, or build/test.

### C. Implement
Make the smallest complete root-cause fix. Add a regression test. Update
declarations, VADAT parsing/defaults, consistency checks, diagnostics, and
documentation together when controls or behaviour change. Explain numerical
choices and safeguards in comments.

### D. Verify
Run focused tests, then relevant suites/builds. Inspect output, exit codes,
warnings, and runtime diagnostics. Check disabled/new modes. Independently
verify the final physical residual and state for solver changes. Do not
suppress or reinterpret failures.

### E. Report
Finish with:
1. What changed (files and behaviour).
2. Why (root cause and scientific/numerical rationale).
3. Tests run (exact commands and outcomes).
4. Unrun tests and why.
5. Remaining risks and missing runtime validation.
6. Git branch/commit state; do not claim a push/merge unless done.

Separate verified facts from hypotheses. Cite paths/routines and line numbers
where available.

## 13. Do not do these things
- Do not perform broad refactors during targeted numerical fixes.
- Do not assume comments are correct without checking implementation.
- Do not infer runtime convergence from compilation/unit tests.
- Do not call a solver robust merely because it has clipping, damping, floors,
  restarts, or fallback paths.
- Do not let a safety cap masquerade as convergence.
- Do not change physics, units, signs, or defaults without documenting and
  testing consequences.
- Do not weaken reference comparisons or leave a half-wired feature.
- Do not claim the entire treatment is validated from one manufactured ray
  or one successful model.

## 14. Scientific references

Use the repository README and relevant source comments/docs for exact citations.
References listed in the README include:
- Hillier & Miller (1998), CMFGEN, ApJ 496, 407.
- Knop, Hauschildt & Baron (2009), A&A 496, 295 and A&A 501, 813.
- Dessart, Audit & Hillier (2015), MNRAS 449, 4304.
- Feldmeier & Nikutta (2006), A&A 446, 661.
- Knoll & Keyes (2004), J. Comput. Phys. 193, 357.
- Arramy, de la Cruz Rodriguez & Leenaarts (2024), A&A 690, A12.
- Eisenstat & Walker (1996), SIAM J. Sci. Comput. 17, 16.
- Kelley (2003), *Solving Nonlinear Equations with Newton's Method*.

Do not add literature claims or alter equations based on unverified citations.
Check original papers and ensure implemented discretizations match their
assumptions and conventions.