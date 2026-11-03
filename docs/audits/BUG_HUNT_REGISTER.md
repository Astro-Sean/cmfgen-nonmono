# Bug Hunt Register - CMFGEN non-monotonic fork

Independent source-level bug hunt ordered by the project owner.
One entry per confirmed or suspected defect. A bug is marked fixed
only after its regression test passes on the corrected code and
fails on the pre-fix code. Severity: P0 = wrong physics/numerics
silently accepted; P1 = wrong result in exercised configurations;
P2 = latent defect on less-used paths; P3 = cosmetic/diagnostic.

Status summary:

| Bug ID  | Severity | Title                                            | Status |
|---------|----------|--------------------------------------------------|--------|
| BUG-001 | P0       | GMRES convergence claim accepted on hess_res     | FIXED  |
| BUG-003 | P3       | MKL "Parameter 4" warning on entry to DGELSY     | REFUTED|

------------------------------------------------------------------------

## BUG-001 - GMRES `info=0` convergence could bypass true-residual
##          verification (`JFNK_TRES_FAC` defaulted off)

- Severity: P0
- Status: FIXED (regression tests pass on fixed code; T4k(biased) and
  G17 confirmed to FAIL on the pre-fix solver)
- Reported: 2026-11-03 (adversarial bug hunt, Phase 1)
- Fixed: 2026-11-03 @ working tree (pre-commit)

### Affected source files and routines

- `new_main/mod_subs/jfnk_solver_v1.f90`
  - module variable `jfnk_gmres_tres_fac` (declaration ~line 174)
  - `gmres_m_v1` done-exit block (~lines 2924-3014 post-fix)
  - `gmres_m_v1` exit measurement block (~lines 3016-3054)
  - `gmres_m_v1` restart-cap exit (~lines 3160-3212)
  - `newton_hook_v1` caller gate (~lines 1182-1201)
- `new_main/mod_subs/rd_control_variables.f` (`JFNK_TRES_FAC`
  VADAT default and clamp, ~lines 1718-1721)
- `new_main/mod_subs/control_variable_mod.f` (doc comment)

### Exact source lines (pre-fix, HEAD `6410d07`)

- `jfnk_solver_v1.f90:174` - `jfnk_gmres_tres_fac = 0d0` module
  default, documented as "0 = off preserves the historical
  accept-on-hess_res path".
- `jfnk_solver_v1.f90:2933` - `if(jfnk_gmres_tres_fac > 0d0)` armed
  the true-residual measurement; with the default the block never
  ran, so `done` exits proceeded on `res <= tol` (the Hessenberg LS
  residual) alone.
- `jfnk_solver_v1.f90:3019-3028` - the exit measurement was guarded
  by `info /= 0`, so an ordinary `gmres_converged` return carried
  `jfnk_gmres_true_res_ok = .false.` and no measured residual.
- `rd_control_variables.f:1718-1721` - VADAT default `0.0` with a
  `LT 0 -> 0` clamp; no run configuration in the repository armed
  the gate, so production never measured the true residual on a
  claimed convergence.

### Suspected defect (confirmed)

A claimed GMRES convergence (`res <= tol`, `info=0`) returned
`gmres_converged` based solely on the Hessenberg least-squares
residual. The same source comments document measured hess_res values
15x-450x below `||b - A x||_w` on exits (OUTGEN runs; documented
`RESTART DIAG` lines; convergence_issues.md records gaps to 3755x).
The caller gate in `newton_hook_v1` only fires on `ginfo == 1`
(stagnation-class exits), so it could not catch an `info=0` claim.
A secondary path - the `ginfo == 2` hookstep-recompute return -
converged-coded with `true_res_ok = .false.` as well (unverified by
design of that path; see BUG-001 notes below).

### Physical / numerical consequence

For an interacting-SN population solve, an accepted phantom
convergence means the Newton direction is accepted without
satisfying the inexact-Newton forcing condition
`||M^{-1}(F - J s)||_w <= eta * ||M^{-1}F||_w`, which the
Eisenstat-Walker convergence theory requires. Consequences observed
and plausible:

- Newton commits steps whose true linear residual exceeds the
  requested forcing tolerance, degrading or falsifying the
  convergence-rate guarantees the solver relies on.
- In the documented 15-450x-gap regime the returned step can carry
  almost no real residual reduction (verified in the test fixture:
  `true_res = ||b||` exactly) while being logged as
  `GMRES: CONVERGED`.
- The false claim also truncates the GMRES early: Krylov iterations
  that would have kept reducing the residual are abandoned.

### Reproduction method

1. `tests/jfnk_core/test_guards.f90` G17: `df_phantom` returns the
   identity action to unit-norm Arnoldi basis calls, collapsing
   `h(2,1)` at iteration 1 (lucky breakdown, `hess_res ~ 0`), and
   zero action to accumulated-step calls (`|v| >= 3`), so the
   measured residual is `||b|| = sqrt(NG) = 7.746`. With
   `JFNK_TRES_FAC=0` (the pre-fix default semantics) the pre-fix
   code returns `gmres_converged` and Newton commits the phantom
   step; observed exit `too_many_rejected_steps`,
   `jfnk_gmres_true_res = 0` (never measured), 11 JVP evals.
2. `tests/jfnk_core/test_newton_core.f90` T4k(biased): `df_lin`
   plus a fixed offset `matvec_bias`; the LS residual converges
   while `||b - A x - c0||` keeps the offset floor. With `fac=0`
   the pre-fix code returns `info=0 gmres_converged`.

Verified by linking the unchanged test programs against the pre-fix
`jfnk_solver_v1.f90` (`git show HEAD:...`): G17 and T4k(biased)
both FAIL there; all other tests behave identically.

### Resolution

Unconditional true-residual measurement at every `done`-exit was
adopted (the alternative "armed-by-default" policy still leaves an
off-switch a user can set, and no safe value of the knob disables
correctness):

- The `if(jfnk_gmres_tres_fac > 0d0)` gate in `gmres_m_v1` is
  removed; the measurement runs at every done-exit for one extra
  matvec per exit.
- `JFNK_TRES_FAC` is now the acceptance-bound multiplier on the
  requested tolerance: a `res <= tol` claim is accepted only when
  `true_res <= acc_fac * tol`, where `acc_fac` floors non-positive
  input to the default `1.0`. `> 1` relaxes the verified bound
  (weakens the inexact-Newton guarantee; documented diagnostic use
  only). No value disables verification.
- Defaults updated to `1.0`: module declaration, VADAT default in
  `rd_control_variables.f` (`LE 0 -> 1.0` clamp so an explicit `0`
  cannot silently restore the pre-fix behaviour), and the
  `control_variable_mod.f` doc.
- Non-finite `true_res` still exits `gmres_tres_nonfinite`; a JVP
  soft failure in the check exits `gmres_jvp_soft_exit`; a
  preconditioner clip exits `gmres_pc_clip_tres`; a rejected claim
  restarts on the measured residual bounded by `JFNK_GMRES_RST` and
  `imx`, else `gmres_tres_reject`. All fail closed - an unverifiable
  solve is never a convergence.
- The `info /= 0` measurement block below the TRES check is now a
  safety net only (fires if a future path reaches the exit
  unmeasured); comment updated.

Weighted-norm convention (verified against the solve convention):
`w = b - A x` with one fresh `matvec`; under left preconditioning
`psolve` is applied to `w` so the measurement is
`||M^{-1}(b - A x)||_w`, matching the left-preconditioned system
GMRES actually solves. Under right preconditioning the raw
`||b - A x||_w` is used, matching the right-preconditioned
residual. The norm is `sqrt(dotprd(n,w,w))` in the same weighted
inner product as the Arnoldi basis. This convention was already
implemented and audited; the fix makes it run unconditionally.

### Regression tests

- `tests/jfnk_core/test_newton_core.f90`:
  - T4j (new): biased claim rejected at the default bound `fac=1.0`,
    restarts bounded (`evals <= 40`; observed 23).
  - T4k (new): `fac=0` floors to the default bound - consistent
    matvec still `gmres_converged` (`max|x-xs| = 1.9e-11`), biased
    matvec exits `gmres_tres_reject` (no value disables the check).
  - T4a-i (pre-existing): continue to pass, covering accept,
    reject+restart, reject-at-cap, left-PC `M^{-1}` path, NaN and
    +Inf fail-closed exits.
  - T8 (amended): the `cap_row1` fixture has a genuine FD
    true-residual floor at `|F|`; strict `fac=1.0` legitimately
    rejects every claim and both A/B runs abort before the
    comparison runs. T8 now sets `JFNK_TRES_FAC=1.5` explicitly -
    the documented `>1` relaxation - since the test exercises
    step-disposition mechanics, not solve acceptance. Result:
    capped run reaches `|F|=1007` vs uncapped `2810` (PASS).
- `tests/jfnk_core/test_guards.f90`:
  - G17 (new): end-to-end phantom-claim regression with `fac=0`
    (pre-fix default semantics). Fails on the pre-fix solver
    (`too_many_rejected_steps`, `true_res=0` unmeasured); on the
    fixed solver the claim is rejected inside GMRES, the caller
    aborts `gmres_stagnation_abort`, `jfnk_nits=0`, `jfnk_x`
    untouched, 5 JVP evals.
  - G13-G15 `tres_fac` initialisations updated `0 -> 1.0` to state
    the verified bound explicitly (behaviour identical, since `0`
    now floors to `1.0`).

### Validation evidence

- `tests/jfnk_core` under `ifort` (Intel 2024 + MKL): T1-T8 all
  PASS, G1-G17 all PASS, precond-blocks, ptc-sign, trace-state all
  PASS (full log in the build transcript).
- Negative verification: identical test binaries linked against
  `git show HEAD:new_main/mod_subs/jfnk_solver_v1.f90` - G17 FAILS
  (`too_many_rejected_steps`, `true_res=0`, 11 evals) and
  T4k(biased) FAILS (`info=0 gmres_converged` on a biased solve).
- `gfortran` configuration could not be exercised: the site gfortran
  wrapper fails with a GLIBC mismatch before compiling; noted as a
  validation gap, not a code issue (the change uses no
  compiler-specific features).
- `./build.sh` production build: see the build log; the solver
  object compiles clean under ifort and the build's embedded
  `tests/f90` + `unit_tests` suites run at the end of the build.

### Remaining notes / follow-ups

- `JFNK_TRES_FAC > 1` remains a supported relaxation; production
  runs should record `tres_rej` counts (task P4-010's runtime leg)
  to quantify how often claims fail verification on real models.
- The `ginfo == 2` hookstep-recompute return exits `info=0` with
  `true_res_ok = .false.` by design: it is a caller-requested
  trust-radius rescale on the retained Krylov basis, not a
  `res <= tol` claim. Its step is verified downstream in the caller
  by `pred_res_physical = ||F - J s||` (a fresh `df` eval -- the
  same physical residual the TRES check measures, in the
  unpreconditioned norm) plus the descent, runaway and trial
  `ared/pred` checks, so every committed step carries a measured
  residual one way or the other. Left as-designed.
- `unit_tests`/`unit_tests/systems` production-model probes with the
  fixed default are part of the final validation phase; a build-only
  verification is recorded, not a physical one.

### Surrogate suite recalibration (test_se_solver.py)

Three manufactured-system tests were calibrated on trajectories that
implicitly assumed phantom acceptance; under verified solves they
starved and failed. All were recalibrated on real verified
trajectories, not on weakened acceptance:

- T11 (`sestrike`): the LU-failure -> TSVD-commit precondition no
  longer holds on `hard` (both legs abort under verification).
  Gated to `std`+`shock`, where the commit trajectory still fires.
  Driver comment updated; `mild`/`hard` direct-call failures are
  expected (gate excludes them).
- T12 (`shock`): the T march now needs ~100 admissions to reach the
  RE floor (~2.3x longer than the phantom-accepted calibration).
  `se_driver.f90` `seqstr` was 64 chars, truncating the sequence;
  expanded to 256 (test infra only). T12 passes on verified solves:
  RE 3.817e7 -> 1.526e5, `x_T` -> 2.998 (equilibrium 3.0),
  undamped eps0 1.875 vs damped 0.144.
- T13 (`shock`): reshaped to `D*20 + P + D*3 + P + D*15`,
  `rel=0.7`. The `P` re-shocks induce commits mid-sequence, firing
  the locked -> RELEASE(rung2) -> REVERT(eps_t=1.68) -> locked arc
  on real solves. Scalar renamed (`rel` shadowed the RELEASE-lines
  list); AUTO-range lookups guarded for non-admission `P` rows.

`test_jfnk_solver.py` (wraps `test_se_solver.py --quick`):
89 passed, 0 failed. The 4 remaining `XFAIL` diagnostics in the full
suite are pre-flagged fixture-expectation gaps, not correctness
gates.

------------------------------------------------------------------------

## BUG-003 - MKL "Parameter 4 was incorrect on entry to DGELSY" warning

- Severity: P3 (diagnostic)
- Status: REFUTED - no defect; warning is MKL's arg-check on a
  deliberately NaN-poisoned input, handled correctly downstream.
- Reported: 2026-11-03 (build transcript, Phase 1)

### Finding

`Intel oneMKL ERROR: Parameter 4 was incorrect on entry to DGELSY.`
fires exactly once in `tests/f90/test_jfnk_solver` T7 (`prob_kind=3`,
`tf` returns all-NaN F). The NaN propagates into the Hessenberg
matrix handed to `dgelsy`.

### Investigation

Both production `dgelsy` call sites in `gmres_m_v1` have provably
adequate dimensions:

- Ordinary: `dgelsy(j+1, j, 1, h_, m+1, p, m+1, piv, ...)` with
  `j <= m`, so M=j+1 <= m+1=LDA and LDB >= M.
- LGMRES-augmented: `dgelsy(j+1+lg_k, j+lg_k, 1, hd_aug,
  m+1+lg_kmax, ...)` with `lg_k <= lg_kmax`, so M <= m+1+lg_kmax.

The routine copies `h_aug` into `hd_aug` before the call (dgelsy
overwrites its matrix argument) so the basis is preserved for
residual reconstruction and hookstep. "Parameter 4" = the matrix A
argument; MKL's entry check flags the NaN contents.

### Disposition

The dgelsy error return is handled: `info=-4` maps to
`gmres_dgelsy_error` -> `info=1` (not converged), and T7 asserts
`info /= 0` - PASS. Cosmetic; pre-existing in the legacy
`tests/f90` harness (fires identically on the pre-fix solver). No
production change required. A NaN-free variant of T7 could silence
the warning but is out of scope for this fix.

------------------------------------------------------------------------

## Template for subsequent findings

- Bug ID: BUG-###
- Severity: P0 | P1 | P2 | P3
- Status: SUSPECTED | CONFIRMED | FIXED | REFUTED
- Affected source files and routines:
- Exact source lines:
- Suspected defect:
- Physical or numerical consequence:
- Reproduction method:
- Regression test:
- Resolution and validation evidence:
