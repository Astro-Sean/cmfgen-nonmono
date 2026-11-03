# CMFGEN Non-Monotonic Numerical Architecture Audit

Date: 03-Nov-2026

## Status and scope

This audit records the current `cmfgen_jfnk` source tree and its regression
results. Production-model validation is outside its scope.

- Source HEAD: `e4c8f02b` on `main`, tracking `origin/main`.
- The source worktree was already modified before this audit. Existing edits
  include JFNK h-adaptation work and its diagnostic tests.
- With user approval, the unsafe BA row/column edit was narrowed to exactly
  disconnected zero equations with zero RHS. A production helper and focused
  test were added; the focused test passes. `./build.sh` completed with zero
  errors and the helper symbol is present in `exe/cmfgen_dev.exe`.
- The numerical-method development plan at
  `docs/NUMERICAL_METHODS_DEVELOPMENT_PLAN.md` was already present and
  untracked. It was left unchanged.
- The outer repository also has unrelated local changes. They were left
  untouched.

The requested implementation paths do not all match this tree. In this tree,
`jfnk_cmfgen_interface_v1.f` and `jfnk_solve_corr_v1.f` are in
`new_main/mod_subs/`; `cmfgen_sub.f` is in `new_main/`. The production-source
Python tests are under `unit_tests/`.

## Current solver sequence

1. `new_main/cmfgen_sub.f` loops over frequency in the great iteration and
   calls `COMP_J_BLANK` for transfer and J updates. `DO_NONMONO_KNOP` selects
   the sign-aware Knop formal solver where supported. `DO_NONMONO_SOB` changes
   Sobolev line treatment; its continuum path depends on `USE_LAM_ES` and
   `DO_SOB_FG_KNOP`.
2. With the non-monotonic moment option active and `USE_LAM_ES=F`,
   `COMP_J_BLANK` calls `MOM_J_CMF_NONMONO_V1`. With `USE_LAM_ES=T`, the J
   update instead uses formal-solver ray moments. These are separate operator
   paths and need separate validation.
3. After the frequency sweep, `new_main/cmfgen_sub.f` may call
   `FG_KNOP_REFRESH_V1` or `MOM_KNOP_REFRESH_V1` to refresh lagged fields.
   The lagged state and its maturity metrics then affect population-solver
   admission.
4. `SOLVE_FOR_POPS` selects the existing LAMBDA/BA path or
   `JFNK_SOLVE_CORR_V1`. Its gate considers population demand, the selected
   frozen/coupled tier, and field maturity.
5. Frozen-J JFNK evaluates the population residual with stored radiative terms
   and the assembled radiative matrix contribution. Coupled JFNK runs the
   registered radiation sweep at each residual evaluation before forming the
   statistical-equilibrium residual. Thus the coupled tier includes a
   radiation response that the frozen-J BA matrix does not represent in full.
6. JFNK builds one vector with `NT` entries per depth. Population, electron
   density, and temperature components share this layout; fixed components
   are masked. The residual is assembled from `STEQ`, the selected radiation
   contribution, and the pseudo-transient term when enabled. A frozen-J
   derivative and a coupled derivative are not interchangeable.

## Active linear and nonlinear methods

- JFNK uses a forward finite-difference product in `JFNK_JVP_V1`, standard
  restarted GMRES, and inexact-Newton forcing. The finite-difference step is
  scaled by state and direction norms. Boundary-distance caps and retry
  reductions can change the effective step by direction.
- The main preconditioner is an equilibrated per-depth BA block solve with LU,
  least-squares, or Jacobi variants. Optional depth bands and the
  `n,T`-oriented Schur path extend this approximation. Low-rank enrichment and
  secant updates are also optional.
- Within a GMRES solve, the factorization is normally fixed. Non-finite or
  clipped preconditioner outputs are counted and cause a solve failure instead
  of being passed onward as a usable Krylov vector. The code therefore has no
  current evidence that FGMRES is needed for a changing preconditioner.
- LGMRES and harmonic-Ritz deflation/recycling are already present as optional
  paths. Their presence alone is not evidence of benefit; compare equal
  residual-evaluation budgets and validate the true residual.
- Anderson-style acceleration exists in the Knop/moment refresh machinery.
  Do not add a second acceleration layer until its activation, rejected
  candidates, fixed-point residuals, and cost are measured.
- The non-monotonic moment solver splits frequency coupling by sign and uses
  live versus lagged upwind values. It has a manufactured coefficient and
  tridiagonal-solve test, but the current suite does not establish
  frequency-integrated conservation on near-coincident grids.

## Findings and recommendations

### 1. Krylov operator consistency and residual acceptance

Classification: critical correctness investigation.

`JFNK_JVP_V1` computes a fresh finite-difference perturbation for every
direction. In addition to the base norm scaling, clipping-distance caps and
clip retries may select different effective perturbations. This can make the
finite-difference action depart from one fixed linear operator within an
Arnoldi cycle. The GMRES implementation records a Hessenberg residual and can
explicitly evaluate `b-Ax`, but `JFNK_TRES_FAC` defaults to zero, so the true
residual is not checked for acceptance by default. The current h-adaptation edits
in the worktree select a per-mode, per-baseline scale, but the regression
changes are static checks; they do not establish cycle linearity or residual
accuracy.

The logged 150-450 ratio between the Hessenberg estimate and explicit residual
is a strong reason to finish this investigation. An explicit callback
evaluation can still use direction-dependent differencing, so it may not
verify a fixed linear residual.

Recommendation:

- Add production-interface directional derivative and additivity tests across
  frozen and coupled residuals, repeat calls, positive and negative directions,
  clip margins, masks, and a perturbation sweep.
- Compare one cycle-consistent perturbation scale against current behavior.
  Detect rejected/nonsmooth directions and restart before mixing them into an
  existing basis.
- Track estimated and explicitly evaluated residuals at every GMRES exit.
  Define the acceptance rule from the inexact-Newton forcing term, not a
  separately guessed tolerance. Keep the original nonlinear residual as the
  final acceptance check.
- Do not claim the current optional `JFNK_TRES_FAC` closes the issue until its
  default, cost, and relation to the forcing term are tested.

Literature: Knoll & Keyes (2004); Eisenstat & Walker (1996); Saad & Schultz
(1986); Pernice & Walker (1998).

### 2. Exact-null handling in the BA solve

Classification: correctness fix made; production build pending.

The initial worktree edit in `new_main/subs/cmf_blkband_v3.f` replaced each
exact zero row by an identity row and set its RHS to zero. It also treated an
exact zero column as grounds to erase the corresponding row. That changed the
linear system and could discard a nonzero SE residual.

The revised `CMF_BLKBAND_PIN_NULL_V1` only pins an index when its row and column
are exactly zero and its RHS is zero. This completes a free variable with a
zero correction without changing equations for other variables. A nonzero
RHS, a zero column with a live row, a zero row with a live column, and a
near-zero entry all keep the original matrix and RHS. The existing LU failure
path handles other singular or inconsistent cases.

The helper is not flag-gated. Ordinary matrices use the same solve; only the
exact isolated `0=0` case changes. `test_cmf_blkband_pin.py` compiles and runs
the production helper on zero-equation, nonzero-RHS, zero-column,
zero-row/live-column, and near-zero cases. All six assertions pass.

Post-change verification (same day): `./build.sh` (`make KEEP_EXE=1 clean &&
make DEBUG=nodebug all`) completed with exit 0, zero errors, and produced
`exe/cmfgen_dev.exe`, `exe/cmf_flux.exe`, and `exe/obs_frame.exe`. The full
unit suite passed 20/20 test files, including the BA-adjacent solver tests
and `test_sob_prof.py` after `lib/libdisp.a` was rebuilt.

### 3. Near-coincident frequency points and the moment resonance

Classification: high-priority discretization investigation.

`MOM_J_CMF_NONMONO_V1` uses `MOM_DLN_MIN=1e-6` in `chi_hat`, `GAMH`, and `GAM`.
The source comments identify bound-free edge doublets with intervals near
`2e-11`. When active, the minimum changes the `1/dln(nu)` advection coefficient
and therefore changes the discrete transfer operator. The floor changes the
transport term as well as the iteration behavior. Existing tests confirm
coefficient assembly but do not
verify frequency-integrated conservation around the threshold.

The same solver replaces small nonzero values of
`1 + sigma*K_ON_J` with a signed `MOM_GAM_FAC_MIN=1e-2` factor, while exact
zero is set back to zero. This changes the transport coupling near the
resonance and has a special exact-zero branch. Existing manufactured tests
include resonance behavior, but a complete below/at/above threshold and
frequency-balance study remains necessary.

Recommendation:

- Derive the limiting discrete frequency equation for a true small interval,
  including edge doublets and both upwind directions. Measure activation
  counts and the integrated balance before proposing a replacement for
  `MOM_DLN_MIN`.
- Test the moment closure below, at, and above the `MOM_GAM_FAC_MIN`
  threshold, including sign changes and exact zero. Compare the original
  equation, the bounded coefficient, and the solved field.
- Do not remove either bound until the derivation and production-source tests
  pass.

Literature: Knop, Hauschildt & Baron (2009a, Paper I); Knop, Hauschildt &
Baron (2009b, Paper II); Heller (cyclic reduction), only if a derived stencil
supports a block-frequency solve.

### 4. LU and preconditioner diagnostics

Classification: high-priority numerical investigation.

`test_jfnk_solver.py` marks T4 as an expected failure with the description
"core LU solve residual <10%". In the hard surrogate, the measured `MZDIFF`
was `1.68e11` and `ZINF` was about `1e35`. This is a failed preconditioner
application-quality check and should not be reported as an ordinary pass.
However, the metric includes rows that the test preconditioner deliberately
projects to zero after regularization, so it is not by itself a clean LU
backward-error measurement. The separate least-squares apply identity checks
were near roundoff (`RELDIFF` at or below `4.1e-16` in the hard case).

The existing `test_precond_blocks` includes a backward-error calculation, but
the displayed result is not an asserted threshold. The current evidence does
not isolate factorization accuracy from row projection, regularization,
ill-conditioning, or an inconsistent test matrix.

Recommendation:

- Split the test into an exact factor/apply backward-error test on active
  rows, a separate projected-row residual report, and a preconditioner quality
  test against the full residual.
- Assert factorization and solve status, finite output, and scaled backward
  error on well-conditioned known systems. Keep near-singular systems as
  separate conditioning tests.
- Preserve the production physical residual; a poor preconditioner is not a
  reason to regularize or replace that residual.

Literature: Knoll & Keyes (2004); Benzi, Golub & Liesen (2005).

### 5. Block coupling and cost

Classification: potential efficiency improvement.

The default depth-local BA factorization does not include all cross-depth
radiation coupling. Optional depth-band storage and low-rank enrichment
provide partial corrections. The current `SCHUR` name refers to the
population/thermodynamic `n,T` coupling, not a complete
radiation/population Schur complement. A new dense radiation Schur matrix
would add substantial memory and setup cost without current measurements of
off-depth coupling.

Recommendation:

- Measure full residual/JVP coupling by depth, the stored band contribution,
  factor and apply time, GMRES iterations, and true residual reduction.
- Use current local population blocks and stored radiation bands as candidate
  approximations only after these measurements. Prototype a radiation Schur
  term separately; do not construct a full dense Jacobian.
- Retain standard GMRES while the preconditioner is fixed in each solve.
  FGMRES handles a varying preconditioner; it does not fix an inconsistent
  Jacobian-vector product.

Literature: Benzi, Golub & Liesen (2005); Saad (1993).

### 6. Outer ordering and existing acceleration

Classification: high-priority convergence measurement.

The existing order is frequency transfer and J updates, optional lagged-field
refresh, then the population update. The maturity gates address stale or
moving lagged fields, but the current regression suite does not establish the
complete population/radiation/temperature convergence sequence on a physical
interacting-supernova model.

Recommendation:

- Report separate field-refresh, statistical-equilibrium, radiative-equilibrium,
  population-change, and temperature residuals at termination.
- Measure current Anderson acceptance and rejection on the actual refresh
  maps before changing their history length or adding acceleration.
- Keep transitions based on field residuals and population demand, not
  iteration counts alone.
- Treat frozen-J as a population solve at fixed radiation, and coupled-J as a
  different residual evaluation. Validate both modes separately.

Literature: Walker & Ni (2011); Arramy et al. (2024); Rybicki & Hummer (1991).

### 7. Methods not justified for production yet

- A global block-frequency solve or cyclic reduction is a research prototype.
  First document the actual stencil, frequency boundaries, near-coincident
  points, and conservation law.
- Nonlinear continuation is a research prototype. The final stage must recover
  the original equations and input velocity field.
- New Krylov deflation or recycling options are not recommended until current
  LGMRES/deflation paths show repeatable benefit at equal residual-evaluation
  cost.
- Complex-step differentiation is not recommended for this residual path
  unless all CMFGEN state, transfer, and rate code supports complex arithmetic.
- Do not add FGMRES to compensate for finite-difference JVP variation. First
  make the Jacobian action sufficiently consistent.

## Baseline tests and benchmark status

Preflight checked the test harness and required production source files before
compilation. No physical model or STELLA input was selected for this audit.

Command:

```text
python3 unit_tests/run_all_tests.py
```

Before the BA null-equation fix, 19/19 test files passed. The suite compiled
and ran production-source
drivers with Intel Fortran 2024.0.2. It exercised the manufactured moment and
Knop-ray paths, VADAT parsing, non-monotonic velocity-gradient routines,
Sobolev paths, and JFNK surrogate systems. The JFNK wrapper classified the
mild T3 fixture check and mild/hard/shock T4 checks as expected failures.

The focused post-fix command `python3
unit_tests/python/test_cmf_blkband_pin.py` passed all six assertions. The hard
JFNK surrogate diagnostic was also run directly. It has 148 unknowns
and an infinite reported condition estimate. Its per-depth LU preconditioner
gave `MZDIFF=1.678e11`; the least-squares algebra check remained near
roundoff. This does not establish production-model behavior and does not
isolate the LU backward error because the diagnostic projects regularized
rows.

A full production build after the BA change passed with zero errors.
No `run_model_grid.py` model, baseline spectrum,
wall-time comparison, or memory comparison was performed. There is no
before/after benchmark because no numerical behavior was changed in this
audit. The older October audit reports different revisions and is not used as
a benchmark for this checkout.

## Ordered implementation plan

1. ~~Run the full build and neighboring BA solver regressions for the
   exact-null helper~~ (done 03-Nov-2026: clean build, 20/20 suite).
   Keep the isolated-zero criteria unless production
   evidence exposes a separate failure.
2. Finish JFNK operator tests: directional derivative, additivity, repeat
   noise, current h adaptation, clip/retry behavior, and explicit residual
   checks in both frozen and coupled modes.
3. Decide the default true-residual policy from those measurements and the
   Eisenstat-Walker forcing term. Verify any change in an analytic GMRES case
   and a production residual driver.
4. Replace the ambiguous LU XFAIL measurement with separate factor-solve,
   projected-row, and preconditioner-quality tests. Keep genuine failed cases
   visible.
5. Derive and test the small-frequency-interval and zero-GAM limits, including
   frequency-integrated balance and threshold sweeps.
6. Profile coupling and solver costs before extending depth bands, Schur
   methods, Anderson acceleration, or Krylov recycling.
7. Run controlled `run_model_grid.py` comparisons only after a clean,
   committed source state and a validated benchmark setup are available.
   Cover monotonic, mild reversal, strong compression, multiple reversals,
   shocked flow, and difficult ionization. Compare physical residuals,
   populations, radiation, temperature, conservation, spectra, time, and
   memory.

Acceptance requires the unchanged monotonic path, verified physical residuals,
no safeguard-dependent false convergence, and repeatable model-level results.
