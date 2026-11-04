# Full Iteration-Cycle TODO - CMFGEN Non-Monotonic Velocity

Companion to `FULL_ITERATION_CYCLE_AUDIT.md`. Audit baseline: `cmfgen_jfnk`
commit `30b596bd744efd7f54dbaf4b5b749796bde0dfb4`.

Each task is independently actionable. Statuses: `[ ] Not started`,
`[~] In progress`, `[x] Complete`, `[!] Blocked`. A task is complete only
when its acceptance tests pass - compiling does not count.

---

## Group 1: convergence-gate defects

### T-01 Field-maturity term in the convergence exit
- Priority: P1
- Findings: F-01, F-12
- Files: `new_main/cmfgen_sub.f` (CONV_EXIT_MET at label 9999), possibly
  `new_main/mod_subs/control_variable_mod.f` for an opt-out key
- Work: add a radiation-maturity conjunct to `CONV_EXIT_MET` for the
  nonmono paths (e.g. `.NOT. (GS_FG_IMMATURE .OR. MOM_UNSAFE_B)` or the
  same predicate the freeze uses), preserving stock behaviour when
  `DO_NONMONO_*` are off. Decide whether pure-RD_LAMBDA runs should inherit
  it or keep the current semantics behind a flag.
- Acceptance: a forced-LAMBDA KNOP run with an immature field does not exit
  converged; a mature-field run exits identically to the baseline; stock
  model (`DO_NONMONO_*=F`) output unchanged.
- Tests: TEST 14 (below) plus one mature-field control.
- Dependencies: none.
- Status: [ ] Not started

### T-02 SOB + USE_LAM_ES=F field-maturity coverage
- Priority: P1
- Findings: F-16
- Files: `new_main/mod_subs/solve_for_pops.f` (GS_FREEZE assignment,
  demote path), possibly `plane/mom_j_cmf_nonmono_v1.f`
- Work: decide whether the KNOP freeze machinery should extend to
  `DO_NONMONO_SOB && .NOT. USE_LAM_ES` (the moment store exists there too)
  or whether SOB gets a lighter immature-field demote. Implement the chosen
  coverage behind the existing SOB flag scope.
- Acceptance: under SOB+LAM_ES=F, an immature `MOM_SW_*` census demotes the
  population step the same way KNOP does; the SOB+LAM_ES=T and stock paths
  are unchanged.
- Tests: unit test driving `solve_for_pops` with a synthetic immature MOM
  census under SOB; runtime check on the m15 moment model.
- Dependencies: T-01 (shared predicate).
- Status: [ ] Not started

### T-03 JFNK on the moment pathway
- Priority: P1
- Findings: F-02
- Files: `plane/mom_j_cmf_nonmono_v1.f` (eval freeze),
  `new_main/mod_subs/solve_for_pops.f` (admission gate),
  `new_main/mod_subs/check_param_consistency.f`
- Work: either implement a `FREEZE_J_NEXT`/`FREEZE_H_NEXT`-equivalent eval
  freeze on the MOM lagged store and admit JFNK under `USE_LAM_ES=F`, or
  formally close the combination with a consistency-check rejection and a
  documented limitation. Pick after estimating the freeze cost.
- Acceptance: the chosen outcome is enforced - either JFNK-CPL evals leave
  the committed store bit-identical across the solve, or VADAT
  `DO_JFNK=T` + `USE_LAM_ES=F` + nonmono stops with a clear message instead
  of silently never firing.
- Tests: eval-freeze determinism test (repeat F(x), same residual);
  consistency-check test for the closed path.
- Dependencies: none for the decision; T-09 for eval determinism testing.
- Status: [ ] Not started

### T-04 SOB frozen-tier self-suppression
- Priority: P2
- Findings: F-03
- Files: `new_main/mod_subs/solve_for_pops.f` (JFNK admission block)
- Work: after `JFNK_REACH_SKIP_STREAK >= JFNK_REL_REACH_MAX`, suppress
  further frozen-tier admissions under SOB (not just the relief) unless
  `JFNK_COUPLED_SOB` can arm the coupled fallback. Evaluate making
  `JFNK_COUPLED_SOB=T` the documented SOB default.
- Acceptance: a SOB run with a radiation-pinned residual stops burning
  frozen admissions after the streak; the coupled fallback still fires
  exactly once per activation.
- Tests: admission-count unit test; runtime telemetry on a SOB model.
- Dependencies: none.
- Status: [ ] Not started

### T-05 Resacc-off BA-on-immature-field path
- Priority: P2
- Findings: F-05
- Files: `new_main/mod_subs/solve_for_pops.f` (stock-cycle GS_FREEZE block)
- Work: decide whether a `!RD_LAMBDA` frozen BA call with resacc disarmed
  should demote to `GS_LAM_CAP` LAMBDA rather than take an unchecked
  `MAX_LIN_COR` step, or whether the resacc-off combination should be
  consistency-rejected under KNOP.
- Acceptance: resacc-off KNOP either demotes on immature field or stops at
  startup with a clear message; no unchecked full BA commit on an immature
  field.
- Tests: unit test for the gate; regression on the warm model.
- Dependencies: none.
- Status: [ ] Not started

## Group 2: numerical verification

### T-06 Runtime instrumentation
- Priority: P1
- Findings: F-22, F-24
- Files: `new_main/cmfgen_sub.f`, `new_main/mod_subs/solve_for_pops.f`,
  `new_main/mod_subs/comp_j_blank.f`, refresh routines, JFNK interface
- Work: counters and timers for great iterations, formal/moment solves, GS
  refreshes (+passes), BA builds vs reuses, JFNK admissions/solves/evals,
  rejected updates, LAMBDA demotions, freeze/release events, wall time per
  stage, peak memory. Emit one machine-readable record per great iteration
  (extend `GS_STATE`/`ITER_HISTORY` or a new file).
- Acceptance: one model run produces a per-iteration ledger whose totals
  reconcile with OUTGEN counters; zero output when nonmono flags are off.
- Tests: ledger-parse unit test; runtime smoke.
- Dependencies: none.
- Status: [x] Complete
- Implementation: `new_main/mod_subs/iter_ledger_mod.f90` (new
  `ITER_LEDGER_MOD`: counters/timers + CSV writer + `/proc/self/status`
  VmHWM reader); instrumented `cmfgen_sub.f` (arm at label 20000,
  FG/MOM refresh call+time), `comp_j_blank.f` (all six J-solver call
  sites incl. EDDFAC ACC pair), `solve_for_pops.f` (stage timer,
  freeze-edge counter, row write next to ITER_HISTORY),
  `fg_j_cmf_nonmono_v1.f`/`mom_j_cmf_nonmono_v1.f` (refresh pass
  totals). Emits `ITER_LEDGER` (one CSV row per SOLVE_FOR_POPS call):
  iter/mode/next/ba/frozen/held/reason, n_fg, n_mom, n_rf, np_rf,
  ba_build, ba_reuse, jfnk, jfnk_cpl, jfnk_ev, jfnk_jvp, n_rej, n_dem,
  n_frz, n_frozen, n_rel, n_raw_rel, maxch, maxch_app, maxatt, se_res,
  se_res2, trust, t_iter, t_rad, t_rf, t_pops, rss_mb.
- Commits: `7571268b` (module), `0d828b26` (cmfgen_sub), `c878e3c1`
  (comp_j_blank), `3fb58aab` (solve_for_pops), `0f039ef4` (fg refresh),
  `81a55d89` (mom refresh), `fa40090d`/`e02559e2`/`d5851a26` (build),
  `cbd61634`/`04d7a3df` (test), `c4f26144`/`233626a7`/`84789e95`/
  `6905ae9f` (test build lists).
- Verification: `./build.sh` clean, all three exes relinked
  (cmfgen_dev, cmf_flux, obs_frame); suite 28/28 with
  `OMP_NUM_THREADS=4` (test_interp_guards hits a 256-thread fork limit
  on the login node without the cap - environment, unrelated);
  `test_iter_ledger` 28/28 checks incl. no-file-when-unarmed.
- Deviations: JFNK evals counted via `JFNK_N_RESIDUAL_EVALS`/JVP deltas
  snapshotted at reset rather than instrumenting the eval interior -
  the sweeps already route through the counted COMP_J_BLANK sites.
  Model-run reconciliation still pending (needs a real run; covered by
  T-08 runtime smoke).

### T-07 Transient-guard influence census
- Priority: P3
- Findings: F-19
- Files: `new_main/cmfgen_sub.f`, `plane/fg_j_cmf_nonmono_v1.f`,
  `plane/mom_j_cmf_nonmono_v1.f`, `subs/sob_line` family
- Work: verify each transient guard (PREVENT_LOW_T, CHIL resets, SUB_SOB_LINE
  zeroed terms, RLUMST floors, ES_ALI DEN floor) reports activation counts
  and that committed iterates activate none of them.
- Acceptance: a converged run's last iterations show zero guard
  activations; any guard still firing at exit is flagged in the ledger.
- Tests: census unit test.
- Dependencies: T-06.
- Status: [ ] Not started

## Group 3: integration test suite (from the audit brief)

### T-08 TEST 1 + TEST 2 + TEST 3: control and compression models
- Priority: P1
- Files: `unit_tests/` (new drivers), `scripts/` run configs
- Work: (1) monotonic control proving stock path identity with all
  nonmono flags off; (2) one weak compression zone converges; (3) sigma<-1
  over a substantial region converges.
- Acceptance: (1) bit-identical or documented-tolerance output vs the
  reference tree; (2),(3) converge within the iteration budget and pass
  TEST 14's residual check.
- Tests: these are the tests.
- Dependencies: T-06 for convergence accounting.
- Status: [ ] Not started

### T-09 TEST 9 + TEST 10 + TEST 15: restart, rollback, reproducibility
- Priority: P1
- Files: `unit_tests/`, `new_main/mod_subs/jfnk_solve_corr_v1.f`,
  `plane/fg_j_cmf_nonmono_v1.f` (RESTORE_FG_INEXT_*)
- Work: restart restores lagged fields and controller state identically;
  a failed Newton step restores populations, module state, and BA matrix
  exactly; repeated F(x) evals and repeated full iterations are
  deterministic within documented tolerance.
- Acceptance: restart from FGINEXT_GS continues the same trajectory;
  forced JFNK failure leaves state bit-identical to the pre-solve
  snapshot; F(x) repeats agree to solver tolerance.
- Tests: these are the tests.
- Dependencies: none.
- Status: [ ] Not started

### T-10 TEST 5 + TEST 6 + TEST 11 + TEST 12 + TEST 13 + TEST 14
- Priority: P2
- Files: `unit_tests/` (new drivers)
- Work: formal vs moment fields on the same model; SOB vs Knop populations
  and spectra; capped/failed refresh cannot report mature; freeze/release
  cycling bounded; fixed vs variable T consistency; small-delta /
  large-residual no-exit.
- Acceptance: per test in the audit brief; TEST 11 and 14 must fail loudly
  on current code where F-01/F-17 predict.
- Tests: these are the tests.
- Dependencies: T-01 (TEST 14 acceptance), T-02 (moment coverage), T-06.
- Status: [ ] Not started

### T-11 TEST 4 + TEST 7 + TEST 8: multi-reversal and JFNK comparisons
- Priority: P2
- Files: `unit_tests/`
- Work: multi-reversal frequency coupling; JFNK vs direct BA final
  residuals/populations/cost; frozen-J vs coupled-J convergence benefit.
- Acceptance: coupled-J demonstrably beats frozen-J on a
  radiation-pinned residual or the tier is documented as diagnosis-only.
- Tests: these are the tests.
- Dependencies: T-03, T-06.
- Status: [ ] Not started

## Group 4: secondary cleanup

### T-12 BA_DEM_EXR exit-hold bounded release
- Priority: P3
- Findings: audit F-08-adjacent (documented limitation in cmfgen_sub.f)
- Files: `new_main/cmfgen_sub.f` (BA_DEM_EXIT_OK), `new_main/mod_subs/solveba_v13.f`
- Work: bound the exit-only relative-floor hold so material unsaturated
  demand cannot hold the exit indefinitely without a release mechanism.
- Acceptance: exit hold releases through Stage-B or a bounded counter, not
  a permanent stall.
- Dependencies: T-01.
- Status: [ ] Not started

### T-13 Document MAXCH semantics contract
- Priority: P3
- Findings: F-07
- Files: `new_main/mod_subs/solve_for_pops.f`, `cmfgen_obsidian/` notes
- Work: write down the three MAXCH semantics (raw demand / applied commit /
  held attempt) and the rule for which each gate must read; add the
  matching vault note.
- Acceptance: every `MAXCH` consumer in `solve_for_pops.f` is annotated or
  verified against the contract.
- Dependencies: none.
- Status: [ ] Not started
