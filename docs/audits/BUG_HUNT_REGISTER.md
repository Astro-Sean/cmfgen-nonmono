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
| BUG-002 | P2       | SOB_RH_GEOM cache fingerprint misses SIGMA        | FIXED  |
| BUG-003 | P3       | MKL "Parameter 4" warning on entry to DGELSY     | REFUTED|
| BUG-004 | P2       | Zero-GAM resonance: floored factor vs exact-zero | REFUTED|
| BUG-005 | P2       | JFNK admissions on m15 always defer before GMRES | REFUTED|
| BUG-006 | P3       | LSTATE save armed after first eval mutations     | FIXED  |

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

## BUG-002 - SOB_RH_GEOM cache fingerprint never covers SIGMA and is
##          not injective for R,V

- Severity: P2 (latent wrong-results defect; no live trigger found in
  the current call graph - see "Production trigger" below)
- Status: FIXED (regression fails on the pre-fix source, passes on
  the fixed source)
- Reported: 2026-11-03 (adversarial bug hunt, Phase 2)
- Fixed: 2026-11-03 @ working tree (pre-commit)

### Affected source files and routines

- `subs/sob_rh_couple.f`
  - `SOB_RH_MOD` cached-signature scalars `RH_VSUM`/`RH_RSUM`
    (pre-fix lines ~68-69)
  - `SOB_RH_GEOM` fingerprint check + rebuild block (pre-fix lines
    ~99-115)

### Exact source lines (pre-fix, HEAD `addf2fcd`)

- `sob_rh_couple.f:99-103` - `VSUM=SUM(V); RSUM=SUM(R)` then early
  `RETURN` when both sums match their saved values within
  `1.0D-25*MAX(1,|sum|)` and `RH_ND.EQ.ND`.
- `sob_rh_couple.f:152-153` - the cached tables this fingerprint
  guards include `RH_T1(I,K)=|1+MU2*SIGMA(I)|` and
  `RH_T1P(I,K)=|1+MUP2*SIGMA(J)|` - i.e. the tables depend on SIGMA,
  which the fingerprint never examined.

### Suspected defect (confirmed)

Two defects in the cache-validity decision:

1. SIGMA is not fingerprinted at all. `SOB_RH_COUPLE` passes either
   `SIGMA` (DO_SOB_CLAMPED_LINES) or `SIGMA_TRUE`; the pair tables
   cache `|1+mu^2*sigma|` at both endpoints. A call with bit-identical
   R,V but a different sigma returns early and silently reuses stale
   resonance-tau tables, so `DRC` is computed with the wrong
   `t=GAMH/|1+mu^2 sigma|` and therefore the wrong remote-resonance
   absorption rate.
2. `SUM` is not injective. Two distinct V (or R) profiles sharing a
   sum also pass the fingerprint, leaving the whole pair selection
   (mu^2 acceptance) and weights stale. Essentially unreachable by
   chance at bitwise level, but correctness depended on luck, not
   construction.

### Physical / numerical consequence

The RH78 remote-pumping correction `DRC_l` feeds `T2=ETAL*ZNET-DRC`
in `SUB_SOB_LINE_V3`. A stale table multiplies each couple by the
wrong Sobolev optical depth - the error is order-unity in `t` and
exponential in `EXPONX(t)`, so even a modest sigma drift corrupts
the correction by a large factor, silently (no diagnostic; the
"SOB_RH: N pairs" rebuild banner simply never prints).

### Production trigger assessment

Audited all writers of the sigma inputs and all `RECOMPUTE_SIGMA_TRUE`
call sites (`cmfgen_sub.f` ~1135/~4650, `solve_for_pops.f` ~3664/3781):
every sigma rewrite is paired with an R or V change (model setup,
`DO_WIND_VEL_V1`, NG acceleration, `ADJUST_R_GRID_V4`), so the stale
path is not reachable on current flags. The sigma argument itself is
fixed per run (flag-selected), so no alternating-input path exists.
Severity therefore P2 (latent): any future sigma-touching call site -
a probe path, an alternating sigma source, a mid-run re-clamp - turns
this into a silent P1 wrong-results defect. The fix removes the whole
class rather than the observed instance.

### Reproduction method

`unit_tests/f90/drv_sob_rh.f90` extended: after the initial
geometry+couple it re-reads ND-line R V SIGMA blocks and re-runs
GEOM+COUPLE in the same process, each dump tagged `PHASE n`. The
`SOB_RH:` rebuild banner (printed inside GEOM only on rebuild) gives
a direct rebuild-vs-stale witness.

`test_sob_rh_couple.py` T7: phase 1 (r,v,sig=-0.5), phase 2 (same
r,v, sig=-0.9), phase 3 (identical to 2). Pair selection depends only
on R,V so the pair set must be unchanged while T1/T1P must rebuild on
the new sigma.

### Pre-fix behaviour (verified)

On the stashed pre-fix `sob_rh_couple.f`: phase 2 produces `banners=0`
(no rebuild) and T1/T1P still hold `|1+mu^2*(-0.5)|` - both new
checks FAIL. Phase 3 passes trivially (stale tables equal stale
tables).

### Resolution

`SOB_RH_MOD` now stores element-wise copies `RH_RSAV/RH_VSAV/RH_SSAV`
of all three inputs; `SOB_RH_GEOM` returns early only when
`RH_ND.EQ.ND` AND `ALL(R.EQ.RH_RSAV).AND.ALL(V.EQ.RH_VSAV).AND.
ALL(SIGMA.EQ.RH_SSAV)` (nested IF - `ALL()` on an unallocated array
is undefined). Bit-identical reuse is exactly what the old
`1D-25*MAX(1,|sum|)` tolerance implemented in practice, so the cache
behaviour on the intended fast path is unchanged; the O(ND) compare
is negligible next to the O(ND^2) rebuild it guards. Fixes defect 1
(SIGMA is now fingerprinted) and defect 2 (equality is injective -
any changed element rebuilds).

### Regression tests

`test_sob_rh_couple.py` T7 (three checks):

- phase-2 same-selection + rebuild banner: FAIL on pre-fix
  (`banners=0`), PASS on fixed.
- phase-2 T1/T1P equal `|1+mu^2*sigma_new|` per pair: FAIL on
  pre-fix (stale -0.5 values), PASS on fixed.
- phase-3 identical inputs: cache hit (`banners=0`) and bit-identical
  dumped tables - PASS on both, confirms the intended fast path is
  preserved.

Suite: 12 passed / 0 failed on fixed source; 10 passed / 2 failed on
pre-fix (the two defect-discriminating checks).

### Remaining notes / follow-ups

- The untracked `knop_solver_files/sob_rh_couple.f` staging copy is
  byte-identical to the pre-fix source and was not touched; it is
  outside git and outside `build.sh`. If it feeds an external KNOP
  package build it should be re-synced from `subs/`.
- Fingerprint discipline for the wider tree: any other cached table
  keyed on inputs should compare all inputs, not digests - worth a
  sweep in Phase 4 (DUMP/RESTORE_*_GS tag+fingerprint items already
  flagged in the TODO around P4-012/INT).

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

## BUG-004 - Zero-GAM resonance: `SIGN(MAX(ABS(T1),1e-2),T1)` floor
##          plus the `T1==0 -> fac=0` exception in MOM_J_CMF_NONMONO_V1

- Severity: P2 suspected (a discontinuous internal coefficient could
  produce a discontinuous or wrong discrete equation at the
  `1 + sigma*K_ON_J = 0` comoving-transport resonance)
- Status: REFUTED - verified regularization, not a defect; the
  exception returns the physically correct decoupled limit and the
  discontinuity is bounded at 2*MOM_GAM_FAC_MIN in the factor.
- Reported: 2026-11-03 (audit directive, Phase 3)
- Closed: 2026-11-03

### Affected source and suspected defect

`plane/mom_j_cmf_nonmono_v1.f`, four factor sites (midpoint current
~2004-2012 and lagged ~2020-2030, depth current ~2126-2135 and
lagged ~2145-2155), all computing

  T1 = 1 + SIGMA_TRUE*K_ON_J            (depth points)
  T1 = 1 + AV_SIGMA_TRUE*NMID_ON_HMID   (midpoints)
  FAC = SIGN(MAX(ABS(T1),MOM_GAM_FAC_MIN),T1)
  IF(T1 .EQ. 0.0_LDP)FAC = 0.0_LDP

with `MOM_GAM_FAC_MIN = 1.0E-2_LDP`. The suspicion: as T1 -> 0+ or
0- the factor jumps to +0.01/-0.01, while exactly at T1=0 the
exception drops it to 0 - a notch in the coefficient at the
resonance, and a sign flip across it.

### Investigation

Traced every consumer of the factors: GAM_SIGNED -> CHI_HAT -> Q_NM,
PSI/PSIPREV (J coupling), W/WPREV -> EPS/HU/HL/HS (H coupling and
eliminated diagonal), the boundary row (line ~2873), and the census
counters.

- The sign flip at T1=0 is inherent, not an artifact: sign(T1) is
  the physical comoving transport direction, so sign(Q) and the
  upwind stencil must change there. The floor only bounds the
  jump magnitude at 2*MOM_GAM_FAC_MIN.
- The `T1==0` exception returns the physically correct limit:
  GAM_SIGNED=0 -> CHI_HAT=CHI -> Q=0 -> PSI=PSIPREV=0, i.e. the row
  decouples entirely from frequency transport where the comoving
  gradient vanishes. That is the right equation at the resonance.
- It does not reintroduce the "implicit diagonal collapse" the
  floor guards against: at exact resonance the row still carries
  the absorption term and the spatial HU/HL couplings; the
  manufactured shelf solves cleanly.
- Consistent at all four sites and under MOM_GAM_PIN
  (PREV = current); the GAM_NM signed->dimensionless overwrite at
  ~line 2213 is dual-use but every post-overwrite read wants the
  coefficient, so no aliasing bug.
- Angle factors are sanitized upstream: K_ON_J clamped to
  [1e-4,1], non-finite NMID_ON_J/NMID_ON_HMID zeroed,
  |NMID_ON_HMID| capped at 1e12, so T1 is always finite and the
  resonance location is physical.
- The sign conventions check out: sign(Q)=sign(T1) marks
  compressing rows; EPS keys on sign(sigma*NMID_ON_J)
  independently of the H direction (documented flux split); the
  compressing-midpoint fixup HS=W/(1+W) keys on GAMH_SIGNED<0.

### Quantitative evidence (t_zero manufactured case)

A prescribed sigma=-3.0 shelf (K_ON_J=1/3) lands 36 cells on
bitwise T1=0: all decouple (Q=PSI=0). The shelf-edge cell computes
T1=3.3e-16 (round-off, not bitwise zero) and takes the floor
fac=+0.01 -> Q=1.6e-4, vs Q~4e-3 and 1.7e-2 at the adjacent
non-shelf cells - i.e. the worst-case notch is ~25x below the
neighbouring coupling and ~4 orders below the row diagonal
(DTAUONQ ~ 3.7). `test_mom_solve.py`: 85/85 PASS - exact resonance
present, Q=0 there, floor bounded off-resonance, no solver
failure, no step discontinuity in J.

### Disposition

No defect at the construction level. On real grids a bitwise
T1=0 is a measure-zero event (sigma must equal -1/K_ON_J exactly);
where it occurs, exact decoupling is the correct physics. The
floor is the regularization that keeps near-resonant rows coupled
at 1% of the local coefficient, and the exception honours the
exact limit. Both are exercised and verified.

Scope note: this refutes the claim that the discontinuous factor
is a wrong-result defect. The resonance-adjacent band is real in
production: m15 sweep dumps (cmfgen_obsidian, LAM_ES=F record)
show depths 1-8 at 1+sigma*K ~ +-0.04 - inside the floor band -
with the upwind direction flipping between frequencies, and those
rows dominate the negative-J quarantine dump. Whether the floor
contributes to the documented LAM_ES=F lagged-map non-contraction
is the open question - that belongs to P3-012/P3-017, not to a
coefficient-defect claim, since the floor is doing its documented
job there.

Recorded observations (no action):

- `GAMH_FAC_NM(ND_NM)=1.0` (~line 2030) is a phantom last-midpoint
  fill; W_NM(ND_NM) is never computed (all loops run to ND_NM-1),
  so the value is never consumed.
- `PASS_NRES` counts bitwise-zero cells in the <=floor census; the
  diagnostic conflates "at floor" with "exactly resonant" by ~one
  cell per event. Cosmetic.
- If a production model ever shows an ill-conditioned row at an
  exactly-resonant cell, the open design question is whether the
  exception should take +-0.01 instead of 0 - but the correct
  fix would be the decoupled row itself, which is already what the
  code assembles.

------------------------------------------------------------------------

## BUG-005 - JFNK admissions on m15 always defer before GMRES

- Bug ID: BUG-005
- Severity: P2 (measured production behaviour, investigated as a
  suspected utility/correctness gap)
- Status: REFUTED - not a defect; the guard is doing its designed
  job. The measurement matters for the open P4 runtime legs.
- Affected source files and routines:
  - `new_main/mod_subs/jfnk_solve_corr_v1.f` (admission, GUARD
    DEFER / reach test, "skipping Newton" path)
  - `new_main/mod_subs/jfnk_solver_v1.f90` (GMRES entry that is
    never reached on m15)
- Suspected defect: that the JFNK tier either runs GMRES on m15
  (so the BUG-001 fix changes production answers) or that the
  defer is a spurious rejection of a reachable solve.
- Measured behaviour: across both surviving production OUTGENs,
  every JFNK admission deferred before Newton:
  - `m15_2p82d_SOB_l14_lte_adp/.../OUTGEN`: 36 admits, 36 defers
  - `m15_2p82d_SOB_l14_lte_adp2/.../OUTGEN`: 119 admits, 119 defers
  - All 155 defers carry `GUARD DEFER: FD non-descent
    cos(F,Az) < fd_cos_min=5.00E-02`, with measured cos(F,Az)
    spanning ~1E-3 to ~1E-2 (occasionally negative). Each
    deferred admission costs ~4 residual evals (F(0) test +
    census + reach eval), is explicitly not counted against the
    FSM rungs ("retaining caller mode (not counted, rung
    unchanged)"), and the outer iteration then advances via the
    LAMBDA/BA step.
  - The defer is physical: the JFNK CENSUS at it485 shows the
    dominant residual rows are radiation-pinned HI Lyman levels
    at the outer boundary; the reach eval reports
    |Jz|/|F| = 26 - no population step can reach the residual
    while J is frozen. The FD-probed JVP direction being
    non-descent is the correct diagnosis for this state.
- Physical or numerical consequence: on m15 (SOB + USE_LAM_ES=T +
  FIX_T), the BUG-001 true-residual path is structurally
  unreachable - GMRES cannot emit `gmres_tres_reject` exits that
  never happen. Production validation of the fix on m15 therefore
  reduces to (a) stability of the admit/defer cycle under the
  fixed build and (b) surrogate/probe evidence for the GMRES path
  itself. A model where the reach test passes (earlier-stage
  demand, or a tier with T free) remains the only route to
  production tres counts.
- Reproduction method: grep of `JFNK ADMIT` vs `skipping Newton`
  in both validation OUTGENs (155/155); per-admission diagnostics
  quoted from it485.
- Regression test: none applicable - no code change.
- Resolution and validation evidence: measured 03-Nov-2026 during
  the interactive-partition stability runs (jobs 12667827 warm,
  12667984 cold, fixed build). The pattern repeats live on the
  fixed build: the warm run's first 10 admissions all carry the
  same `GUARD DEFER: FD non-descent` line (cos ~5E-3 to 2E-2), and
  the cold run's early gate defers (demand ~1E2 vs gate 3.0) show
  the pre-Newton census in action. The warm run additionally shows
  JFNK demand decaying below the frz gate (6.6 -> 0.41 vs gate
  3.0) as the LAMBDA mode re-matures - consistent with the
  admission profile moving away from JFNK territory as the model
  converges conventionally.
- Cross-reference (where GMRES does run): the vault JFNK note
  records the `jte_lgmres_ab` run (coupled tier, KNOP path)
  reaching iteration 912 with 189 admissions, one commit, and
  measured true/Hessenberg residual ratios of 1.09, 1.13 and 1.95
  on its coupled solves. Those ratios are the BUG-001 exposure:
  under the old default (`JFNK_TRES_FAC=0`, acceptance on the
  Hessenberg estimate) claims with true residual up to ~2x the
  target would have been silently accepted. GMRES does enter on
  the coupled/KNOP tier - it is the m15 frozen-SOB tier that
  structurally defers.

------------------------------------------------------------------------

## BUG-006 - `JFNK_EVAL_LSTATE_V1('S')` armed after first mutations

- Bug ID: BUG-006
- Severity: P3 (rare trigger: allocation failure inside the LSTATE
  save; persistent silent corruption when it fires)
- Status: FIXED
- Affected source files and routines:
  - `new_main/cmfgen_sub.f`, `JFNK_RAD_SWEEP_V1` entry-save block
- Exact source lines: the coupled-SOB bracket armed at ~line 5976,
  after mutations beginning at `ACCESS_F=JFNK_EDD_START_REC` (~5861).
- Defect: `JFNK_EVAL_LSTATE_V1('S',IOS_EVAL)` was the only guarded
  save (`GOTO 9999` on alloc failure) positioned after the eval's
  first mutations (`ACCESS_F`, `JFNK_IN_EVAL`, `JFNK_EVAL_BROKEN`,
  `EDDINGTON=EDD_CONT`, `COMPUTE_BA=.FALSE.`, two-photon and SL
  accumulator zeroing). Its failure jumped to the 9999 tail, which
  poisons `SE%STEQ` and returns without restoring - the tail's own
  comment claims "no production state has been mutated at that
  point", true only for the three pre-mutation determinism saves.
  Consequence if it ever fired: `JFNK_IN_EVAL=T` leaks permanently
  (only toggled inside this routine), so every later production
  sweep would run in eval mode - cold FOLD seeds, softened STOPs,
  lagged-J records never overwritten, i.e. the model evolves on a
  frozen radiation field. `EDDINGTON`/`COMPUTE_BA` self-heal at the
  next loop head; the zeroed accumulators lose one iteration of
  contributions.
- Reproduction method: ordering verified by inspection (the guarded
  save sat ~115 lines after the first mutation). Dynamic
  reproduction would require injecting an ALLOCATE failure inside
  LSTATE - not feasible in-tree.
- Fix: hoisted `JFNK_EVAL_SOB_LINES` evaluation and the
  `LSTATE('S')` arm into the entry-save block (before
  `ACCESS_F_SAVE`), keeping `COMPUTE_EW=.FALSE.` at its original
  position inside the eval mutations (LSTATE must still save the
  production flag before that write). Every guarded save now
  precedes every mutation, so the 9999 tail's claim is true on all
  paths. The success path is bit-identical - same saves, same
  mutations, same order between them.
- Regression test: source-ordering pin in
  `unit_tests/python/test_jfnk_regressions.py` ("BUG-006: every
  guarded entry save precedes the first mutation") - asserts the
  LSTATE arm precedes `ACCESS_F=JFNK_EDD_START_REC` and no
  `GOTO 9999` survives between the first mutation and the LUER
  redirect. Verified: FAILS on the pre-fix file (assertion raised),
  PASSES on the fixed file.
- Validation evidence: ifort compile clean; pin fail->pass
  confirmed by stash/restore of the pre-fix source; full suite in
  build log.

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
