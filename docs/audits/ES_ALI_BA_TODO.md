# ES_ALI_BA TODO - acceleration-stack audit follow-up

Each item lists source files/routines and acceptance criteria.
Status: [ ] open, [~] in progress, [x] done, [-] dropped/blocked.
Audit basis: ES_ALI_BA_AUDIT.md (12-Nov-2026, working tree post
`d6b48b24`). Experiment arms E1-E4 defined in audit section 7.

Naming: E1w = warm baseline, E2 = DO_ES_ALI only, E3 = +DO_ES_ALI_BA,
E4 = +DO_BA_RES_ACC. D1-D6 = audit defect list, C1-C6 = concerns.

## Experiments (no code change; YAML + VADAT keys only)

- [~] X1 - Submit E1w-E4 warm-started clones
  Done 11-Nov-2026 as `configs/dessart_m15_2p82d_sobacc_x{1..4}.yaml`
  cloned from `dessart_m15_2p82d_FRESH_sob_dieoff.yaml` (DIE-off),
  each with `cmfgen.warm_restart_from` ->
  `m15_2p82d_FRESH_sob/cmfgen/2.82days` and
  `cmfgen.warm_restart_population_only: true`.
  VADAT arms via `vadat_overrides`:
    X1: DO_ES_ALI: F  DO_ES_ALI_BA: F  DO_BA_RES_ACC: F
    X2: DO_ES_ALI: T  DO_ES_ALI_BA: F  DO_BA_RES_ACC: F
    X3: DO_ES_ALI: T  DO_ES_ALI_BA: T  ES_ALI_BA_MIN: 0.1D0  DO_BA_RES_ACC: F
    X4: DO_ES_ALI: T  DO_ES_ALI_BA: T  ES_ALI_BA_MIN: 0.1D0  DO_BA_RES_ACC: T
  Jobs: FINAL 13036072/74/76/78, OBS 13036073/75/77/79 (PENDING at
  submission). population_only staging verified: NT=1839, ND=70
  fingerprint match; T restore max |dT|=1.879e-02; SN_HYDRO_DATA
  byte-identical (md5 ea6c5d5b) pre-staged per arm.
  Acceptance: four jobs staged; population_only hydro byte-identity
  check passes at prepare time; each run's first OUTGEN shows
  `BA matrix` rebuild (COMPUTE_BA=T) and cold-J first-iteration
  transient consistent with audit section 6.

- [~] X2 - E2 early diagnostic read (~it50)
  Done on X2 run: `freqs unconverged` fell 2496->0,
  `max final dJ/J` -> ~1e-4 (at tol), ALO fallbacks=0,
  J-floor hits run=0, BA-floor hits run=0 across the whole
  run - V13 diagonal served every frequency; neither floor
  engaged on this model. X3/X4 confirm same census plus
  ES_ALI_BA floor hits run=50 (BA side only).

- [~] X3 - E-matrix convergence comparison
  X2 CONVERGED ~it752 / 33 min / ~30 new iters.
  X3 CONVERGED ~it753 / 35 min / ~28 new iters - identical
  within noise; ES_ALI_BA floor hit 50x, no tail gain.
  X1 still running ~it774 maxapp 0.65 with spikes (slow
  tail reproduced). X4 running ~it755, nrej=5, trust
  1e4->312 (CASE C forming). Cold dieoff anchor it241+.
  Metrics recorded per audit section 7 record list in
  SOB_ACCELERATION_RESULTS.md. Fixed-point agreement
  X2-vs-X3 verified at OUTGEN level (~1% on final
  corrections, same worst levels); population-file diff
  still owed.

- [~] X4 - Decision gate after E-matrix
  PRELIMINARY (X1/X4 pending): CASE A confirmed - X2
  converged ~30 iters vs X1's projected 150-700+.
  CASE B not triggered - X3 == X2 within noise.
  CASE C forming - X4 resacc catches bad steps but trust
  collapse costs iterations. If confirmed: propose
  DO_ES_ALI=T for the SOB+LAM_ES recipe (VADAT default or
  config), keep ES_ALI_BA off by default, keep
  BA_RES_ACC opt-in. JFNK-CPL NOT needed on this model -
  proposal doc written anyway (JFNK_CPL_PROPOSAL.md) for
  line-rate-dominated regimes. Still owed: second-model
  confirmation, cold-start ES_ALI arm, population-file
  fixed-point diff.

## Defect fixes (comment/guard level only; no solver-math change)

- [x] F1 - D1: refresh `DO_ES_ALI` module comment
  Done 11-Nov-2026 in control_variable_mod.f; VADAT keys note updated
  (accessor matching, fallback behaviour, KNOP KN_NO_GAIN context).

- [x] F2 - D2: fix stale forcing comment
  Done 11-Nov-2026 in check_param_consistency.f.

- [x] F3 - D3+D4: TRI inert/mixed-operator warnings
  Done 11-Nov-2026 in rd_control_variables.f (warn-only, flag not
  gated off).

- [x] F4 - D5: ES_ALI_BA dependency notice
  Done 11-Nov-2026, stronger than specced: check_param_consistency.f
  now STOPs on DO_ES_ALI_BA=T with DO_ES_ALI=F (user error, not a
  legal-but-misleading key), plus warns when DO_ES_ALI=T but no
  ES-ALI gate can engage.

- [x] F5 - D6: document floor asymmetry
  Done 11-Nov-2026: module comments state 0.01 (J loop) vs
  ES_ALI_BA_MIN=0.1 (BA) differ by design; VADAT keys note updated.
  New run-total counters ES_ALI_NFLOOR / ES_ALI_BA_NFLOOR print in
  the ES_ALI summary line and the do_var_cont ALO banner so X2-X4
  show whether either floor engages on this model.

## Regression tests (audit section 8)

- [x] T1 - Zero-scattering limit (BA side)
  Covered by existing drv_dvc cases (chi_scat=0 -> DENB=1).

- [x] T2 - Scattering-dominated limit
  New `esalif` case in drv_dvc.f90 + test_dvc.py: DJDETA values
  straddle the 0.01 J-loop floor; asserts floor factor and
  continuity. 170/170 pass 11-Nov-2026.

- [ ] T3 - J-loop OAB block compiled coverage
  New drv_* harness stubbing FG_J_CMF_V13 + GET_FG13_ALO_DIAG:
  update formula, RJNEW<=0 fallback to TC, ES_DJ exit, ES_ALI_MAX
  cap, ALO_OK=F fallback. Currently only production-run evidence
  exists (adp2_tauto). Deferred: X2 arm provides production
  evidence first; harness only if fallback census is nonzero.

- [~] T4 - Warm-restart rate-change regression
  Pipeline-level: population_only staging verified live on X1-X4
  (NT/ND fingerprint, T restore, byte-identical hydro). Remaining:
  confirm first-call COMPUTE_BA=T and no BAMAT reuse in arm OUTGEN.

- [x] T5 - Non-monotonic multi-resonance guard
  Audit section 4 established both ALO diagonals come from the
  solver that ran the frequency (GET_FG_ALO_DIAG on KNOP,
  GET_FG13_ALO_DIAG on V13, GET_MOM_ALO_DIAG on moment route);
  same-depth stored response only, no single-resonance assumption.

## Documentation

- [x] A1 - docs/audits/ES_ALI_BA_AUDIT.md (this audit)
- [ ] A2 - Mirror audit status into cmfgen_obsidian/
  Acceleration methods.md + history/convergence_issues.md once E1-E4
  results land (per AGENTS.md 1.3, status entries mirror into vault
  notes).
- [ ] A3 - AUDIT_LOG.md index entry for this audit.
