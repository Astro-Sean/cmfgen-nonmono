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

## 5. Validation gaps

- `gfortran` toolchain unusable on this host (`~/.local/bin/gfortran`
  GLIBC mismatch); all unit verification ran under `ifort`/MKL only.
  The fix itself uses no compiler-specific features.
- No production-model run yet exercises the new default; build-level
  verification plus surrogate/regression evidence is what exists.
  `tres_rej`/`tres_nonfinite` counts from real interacting-SN models
  remain the open physical-impact measurement.
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
