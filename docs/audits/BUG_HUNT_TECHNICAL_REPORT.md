# Bug Hunt Technical Report - CMFGEN non-monotonic fork

Owner-directed independent, adversarial source-level review of the
non-monotonic solver development. This file is the running technical
record; `BUG_HUNT_REGISTER.md` holds the per-defect entries and
`FIVE_PILLAR_IMPLEMENTATION_TODO.md` tracks task status. Scientific
target: interacting-supernova structures with non-monotonic velocity
gradients, compression zones, shocks and multiple velocity reversals,
solved correctly - not merely converged.

Scope so far: Phase 0 repository reconciliation (done) and Phase 1,
the JFNK GMRES true-residual acceptance audit (fixed, BUG-001).

## 1. Phase 0 - repository reconciliation findings

Facts established against git history (authoritative):

- `cmfgen_jfnk` HEAD = `6410d07442dab6aae670a3fa476f7bb9d015e63c`,
  identical to `origin/main`. The TODO's recorded HEAD `ccc52ffd` was
  stale by exactly the seven commits `24694856 8a2a8c0e b12a9c27
  36202011 6b431a6a b9bd047e 6410d074`.
- All seven commits are documentation-only (TODO tracker creation,
  audit documents, one date correction). No production source file
  changed; the audited code at `ccc52ffd` is byte-identical to HEAD.
- Audit-document dates vs commit dates: `JFNK_PRODUCTION_OPERATOR_
  INVESTIGATION_2026-11-07.md` carries a Nov-07 filename for a Nov-03
  document (already corrected in its header, DOC-001).
- Source comments dated later than the repository revision
  (e.g. "12-Nov-2026" edits in `jfnk_solver_v1.f90` and
  `jfnk_solve_corr_v1.f`) were introduced by commits dated
  2026-10-23 - the in-source dates are systematically post-dated
  annotations, not evidence of an unaudited later revision. Git
  history is authoritative; treat in-source edit dates as unreliable.
  (Recorded as DOC-002 in the TODO.)
- Working tree at Phase-0 start showed `build_debug.sh` deleted but
  uncommitted - restored (it exists in HEAD and is a documented
  build entry point). `unit_tests/systems/_pv.txt` remains an
  untracked test artifact.

No audit conclusions had to be re-derived; the stale HEAD did not
hide any code change.

## 2. Confirmed defects

### BUG-001 (P0) - GMRES `info=0` convergence accepted on the
###           Hessenberg model residual alone

See `BUG_HUNT_REGISTER.md` for the full entry. Summary:

- Defect: `JFNK_TRES_FAC` defaulted to `0`, and the true-residual
  measurement in `gmres_m_v1` only ran under `tres_fac > 0`. Every
  ordinary `res <= tol` convergence exit therefore returned
  `gmres_converged` on the Hessenberg LS residual alone - a model
  that has been measured 15x-450x below the true weighted residual
  under the finite-difference matvec. The Newton caller gate only
  covers `ginfo == 1` exits, so `info=0` claims bypassed every
  residual check end-to-end.
- Fix: verification is unconditional at every done-exit;
  `JFNK_TRES_FAC` is now the acceptance-bound multiplier (default
  1.0, non-positive input floored to 1.0, `> 1` a documented
  weakening). All failure modes fail closed
  (`gmres_tres_nonfinite`, `gmres_jvp_soft_exit`,
  `gmres_pc_clip_tres`, `gmres_tres_reject` bounded by
  `JFNK_GMRES_RST`/`imx`).
- Proof: new G17 and T4k(biased) regression tests fail on the
  pre-fix solver (`gmres_converged` on a claim whose measured
  residual is `||b||` in full) and pass on the fixed solver
  (`gmres_tres_reject` -> `gmres_stagnation_abort`, no Newton step
  committed, bounded JVP evals).
- Scientific impact: any JFNK run at the old default could declare
  linear-solve convergence while the returned correction failed the
  inexact-Newton condition; the solver log recorded `CONVERGED` for
  it. In the measured gap regime this silently corrupts the Newton
  globalization evidence chain for exactly the stiff, coupled
  problems JFNK exists for.
- Surrogate consequence: three `test_se_solver.py` tests
  (`se:hard/T11`, `se:shock/T12`, `se:shock/T13`) were calibrated on
  trajectories that implicitly assumed phantom acceptance and
  starved under verified solves. None indicated a solver regression
  - the `shock` fixture shows the documented hess/true gap (0.136
  vs 1.041, ~8x; 9.4e-3 vs 1.455, ~155x in another admission) and
  `info=3` aborts are the correct outcome. Recalibrated on real
  verified trajectories: T11 gated to `std`+`shock` (the LU-fail /
  TSVD-commit precondition no longer holds on `hard`), T12 extended
  to 100 admissions (driver `seqstr` grew 64 -> 256 chars), T13
  reshaped to a `D*20+P+D*3+P+D*15` sequence whose re-shocks induce
  the release/revert arc honestly. `test_jfnk_solver.py`:
  89 passed, 0 failed.

### BUG-002 (P2) - SOB_RH_GEOM cache fingerprint never covers SIGMA

See `BUG_HUNT_REGISTER.md` for the full entry. Summary:

- Defect: the RH78 coupling-table cache decided validity on
  `SUM(V)`+`SUM(R)` within 1D-25 relative, but the cached `RH_T1`/
  `RH_T1P` entries are `|1+mu^2*sigma|` - SIGMA was never
  fingerprinted, and a digest cannot distinguish profiles sharing a
  sum. A same-R,V call with a new sigma reused stale resonance taus
  and computed `DRC` on the wrong `t=GAMH/|1+mu^2 sigma|`.
- Trigger audit: every current sigma rewrite pairs with an R or V
  change (setup, wind-vel, NG accel, grid revision), so the stale
  path is unreachable on present flags - P2 latent rather than
  active. Any future sigma-touching site turns it into a silent P1.
- Fix: element-wise fingerprint on stored copies of R,V,SIGMA
  (`ALL(...EQ...)`); identical-input reuse is preserved (the old
  tolerance was effectively exact-match anyway) and the compare is
  O(ND) against the O(ND^2) rebuild.
- Proof: multi-phase driver + `test_sob_rh_couple.py` T7 - pre-fix
  fails two of three checks (`banners=0`, stale T1/T1P), fixed
  passes 12/12 including the phase-3 cache-hit check.

### BUG-006 (P3) - `JFNK_EVAL_LSTATE_V1('S')` armed after first mutations

See `BUG_HUNT_REGISTER.md` for the full entry. Summary:

- Defect: the coupled-SOB entry save (`GOTO 9999` on alloc failure)
  was the only guarded save positioned after the eval's first
  mutations in `JFNK_RAD_SWEEP_V1`. Its failure path skipped all
  restores - `JFNK_IN_EVAL=T` would leak permanently (only toggled
  inside this routine), putting every later production sweep in
  eval mode: cold FOLD seeds, softened STOPs, lagged-J records
  never overwritten. `EDDINGTON`/`COMPUTE_BA` self-heal at the
  loop head; the zeroed two-photon/SL accumulators lose one
  iteration.
- Trigger: an ALLOCATE failure inside `JFNK_EVAL_LSTATE_V1('S')` -
  rare, but the consequence is persistent silent corruption rather
  than a crash.
- Fix: hoisted the `LSTATE('S')` arm (and `JFNK_EVAL_SOB_LINES`)
  into the entry-save block before `ACCESS_F_SAVE`; `COMPUTE_EW=F`
  stays in the eval mutations since LSTATE must save the flag
  first. Every guarded save now precedes every mutation, which is
  the invariant the 9999 tail assumes. Success path bit-identical.
- Proof: source-ordering pin in `test_jfnk_regressions.py`
  ("every guarded entry save precedes the first mutation") fails
  on the pre-fix ordering, passes on the fix; ifort compile clean.

## 3. Potential defects requiring further investigation

- `ginfo == 2` hookstep-recompute returns a convergence-coded
  `info=0` with `jfnk_gmres_true_res_ok = .false.` - an unmeasured
  returned step by design (it is a trust-radius rescale of an
  already-iterated direction, checked downstream by
  `pred_res_physical` and the trial eval). Not a `res <= tol` claim,
  so outside BUG-001's fix; kept under observation.
- FD-JVP consistency floor: the T8 fixture and production comments
  both show regimes where the measured residual cannot reach the
  forcing tolerance (`true_res ~ |F|`) because the finite-difference
  action is direction-dependent. Strict `fac=1.0` then legitimately
  rejects all claims and aborts to the Lambda fallback. Correct per
  the inexact-Newton criterion, but production runs should record
  `tres_rej` counts to measure how often this bites on real models
  (P4-010 runtime leg).

## 4. Refuted suspected defects

- Phase-0 suspicion that post-Nov-03 source comments implied an
  incomplete audit history: refuted - the comments are post-dated
  annotations; the code under audit is identical to the audited
  revision.
- BUG-003 (P3): the oneMKL `Parameter 4 was incorrect on entry to
  DGELSY` warning. Refuted as a defect - it fires only in the legacy
  `tests/f90` T7 all-NaN-residual scenario, where MKL's entry check
  intercepts the NaN-poisoned Hessenberg. Both production `dgelsy`
  call sites have provably adequate LDA (j <= m on the ordinary
  path, lg_k <= lg_kmax on the augmented path), the error return is
  handled (`info=-4` -> `gmres_dgelsy_error` -> not converged), and
  the warning is identical on the pre-fix solver. Cosmetic noise in
  a deliberate invalid-input test; no production change required.
- BUG-004 (P2 suspected): the zero-GAM resonance construction in
  `mom_j_cmf_nonmono_v1.f` - `SIGN(MAX(|T1|,1e-2),T1)` with a
  `T1==0 -> fac=0` exception on all four angle-factor sites.
  Refuted as a defect. The sign flip at `1+sigma*K_ON_J=0` is the
  physical transport-direction reversal (sign(Q), upwind stencil);
  the floor bounds the coefficient jump at 2*MOM_GAM_FAC_MIN; the
  exact-zero exception returns the physically correct decoupled
  row (Q=PSI=0, chi_hat=chi) which is well-conditioned via its
  absorption and spatial terms. Manufactured-shelf evidence:
  36 bitwise-resonant cells fully decouple while a shelf-edge
  T1=3.3e-16 cell takes the floor (Q=1.6e-4, ~25x below the
  neighbouring coupling, ~4 orders below the diagonal). 85/85
  moment unit tests pass with no J step at the resonance. The
  construction is a documented regularization honouring the
  exact physical limit - keep the `MOM_GAM_FAC_MIN` activation
  census (PASS_NRES) under observation in production runs, and
  see the register for two cosmetic observations (phantom
  `GAMH_FAC_NM(ND)` fill; census counts exact-zero cells).
- BUG-005 (P2 measured): every JFNK admission on the surviving m15
  production runs defers before Newton - 155/155 across the `adp`
  (36) and `adp2` (119) OUTGENs, all on `GUARD DEFER: FD
  non-descent cos(F,Az) < fd_cos_min=5.00E-02` with measured
  cos(F,Az) ~1E-3..1E-2. The dominant residual rows are
  radiation-pinned outer-boundary HI levels; the reach eval gives
  |Jz|/|F| = 26, so no population step can reach the residual
  while J is frozen. The defer is the guard working as designed
  (not counted against FSM rungs, ~4 residual evals each), but it
  means GMRES - and therefore the BUG-001 tres path - has never
  run on m15. Production tres counts are structurally zero there;
  surrogate probes remain the only GMRES exercise.

## 5. Validation gaps

- `gfortran` toolchain unusable on this host (`~/.local/bin/gfortran`
  GLIBC mismatch); all unit verification ran under `ifort`/MKL only.
  The fix itself uses no compiler-specific features.
- No production-model run yet exercises the new default; build-level
  verification plus surrogate/regression evidence is what exists.
  `tres_rej`/`tres_nonfinite` counts from real interacting-SN models
  remain the open physical-impact measurement. As of 03-Nov-2026
  three interactive-partition m15 runs are in flight on the fixed
  build: job 12667827 (warm DIAG restart at it620+, confirms the
  BUG-005 defer cycle live - 15 admits/14 defers + 1 skip-Newton,
  zero tres exits, converging through full iterations), job
  12667984 (cold start, LAMBDA phase at it~38, hours from first
  JFNK admission), and job 12669084 (warm TRI probe cloned from
  the same it620 state, NUM_BNDS=3 + SOL_METH=TRI, queued for the
  P5-013/014 paired-method comparison). Per BUG-005 the DIAG jobs
  are not expected to emit tres exits - the deferred-admission
  pattern means GMRES never enters on this model state.
- Complete unit suite + `build.sh`: recorded in the build log for
  this change set.

## 6. Production-readiness assessment (interim)

- BUG-001 was a genuine P0: convergence evidence could be fabricated
  by the linear solve. The fix removes the only hess_res-only
  acceptance route for `info=0` claims and fails closed on every
  unverifiable path. JFNK convergence claims are now backed by a
  measured residual in the same weighted norm the solve uses.
- This does not clear Pillar 4: FD-JVP operator consistency,
  hookstep-recompute measurement, restart/checkpoint restoration,
  and the coupled-J paths remain open audit items (P4-002..P4-014).
  Sobolev (P2), Knop formal/moment (P3), and DIAG/TRI (P5) are
  untouched by this fix and remain inspection-pending.
