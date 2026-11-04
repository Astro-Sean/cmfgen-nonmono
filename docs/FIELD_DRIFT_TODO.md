# FIELD_DRIFT_TODO - radiation-field drift audit task list

Sequential task list for the field-drift / radiation-field convergence
audit. One task at a time; a task is complete only when the listed
verification ran, not merely when code compiles.

Companion document: `docs/audits/FIELD_DRIFT_AUDIT.md`.

## Status log

### 1. Map the field-drift measurement chain (DONE)

What was done: traced every GS/MOM/FG drift quantity from the solvers
through the maturity getters into `SOLVE_FOR_POPS`.

Files read:

- `new_main/mod_subs/solve_for_pops.f` (getters at ~1028-1042; freeze
  state machine; JFNK admission; SOLVEBA dispatch ~2601-2644;
  population-hold rollback ~2654-2675)
- `plane/fg_j_cmf_nonmono_v1.f` (in-sweep drift ~4168-4330; INIT reset
  ~2256-2309; refresh ~6199-7070; `GET_FG_GS_MATURITY` ~5020-5073)
- `plane/mom_j_cmf_nonmono_v1.f` (per-iteration reset ~1640-1679;
  in-sweep drift ~3405-3828; `GET_MOM_GS_MATURITY` ~4427-4525;
  refresh driver `MOM_KNOP_REFRESH_V1` ~4804-6610; adaptive omega
  ~6500-6593)
- `new_main/cmfgen_sub.f` (refresh calls 3281-3293;
  `SOLVE_FOR_POPS` call 4557; `JFNK_RAD_SWEEP_V1` 5766-7349)
- `new_main/mod_subs/comp_j_blank.f` (solver dispatch 834-1089)
- `new_main/mod_subs/jfnk_solve_corr_v1.f` (`FREEZE_I_NEXT` 531/2392;
  baseline STEQ 554-625; F(0) test 654-684)

Findings recorded in the audit document, sections 3-5.

### 2. Verify the measurement is not self-referential (DONE)

Verified: every metric compares two distinct snapshots
(fresh-solve vs pre-frequency stored; refresh iterate vs previous
pass; committed vs pre-refresh save). No array is compared against
itself; no previous-state buffer is overwritten before use. The
committed store is written after the drift is measured, so the
reported sweep drift is the unrelaxed map residual, not the damped
write.

### 3. Classify the reported "stagnating drift" quantity (DONE)

The headline number (`GS_MAXCHG` / `MOM_GS_MAXCHG`, the value
printed as "largest drift" and used by the maturity gate) is the
max-norm of the *unrelaxed inner-map residual*:

- FG: after the refresh, `FG_GS_MAXCHG = KN_RF_RESID` = max over
  non-diverged rays of last-pass `|VNEW-VOLD|/max(VOLD, 1e-4*VMAX)`.
- MOM: on an unconverged refresh, `MOM_GS_MAXCHG = BEST_CHG` (best
  pass residual).

A plateau of this number is *structurally expected*: the lagged
coupling map is only marginally contractive at a small subset of
cells (strong compressions, |1+sigma*K| resonance cells near the
inner boundary), so its worst-cell residual floors out regardless
of outer progress. This is a genuine non-contraction of a minority
of cells - not a measurement bug - but it means the headline drift
is the wrong instrument for "is the production field converging".
The correct instruments exist (`FG_SW_MAXCHG`/`MOM_SW_MAXCHG` =
fresh-vs-stored sweep drift; `FG_RF_DCHG`/`MOM_RF_DCHG` = committed
retained delta) but are secondary fields, not the headline.

### 4. Audit update sequence, relaxation, JFNK, controller (DONE)

Findings (details in audit sections 6-8):

- Update order per great iteration: freq sweep (relaxed store
  writes, unrelaxed drift census) -> `FG_KNOP_REFRESH_V1` ->
  `MOM_KNOP_REFRESH_V1` -> STEQ -> `SOLVE_FOR_POPS` (gate) ->
  `SOLVEBA_V13`/`JFNK_SOLVE_CORR_V1`. The production store advances
  every iteration; the hypothesis "the field is never updated" is
  falsified at code level.
- Two relaxation systems exist and are distinguishable: the sweep
  write `INEXT_RELAX`/`JNEXT_RELAX` (adaptive, V2 same-map pass
  ratios; V1 cross-iteration ratio) and the refresh inner step
  (`KN_RAY_OMEGA` / per-cell `MOM_DAMP_NM`). The committed refresh
  field is unrelaxed by design.
- Adaptive-omega category mismatch: the V2 estimator measures the
  *inner* (frozen-source) map contraction ratio but sets the
  *outer* (sweep) write factor. A residual floor produces
  rho->1 -> omega->0.5 held indefinitely, which can slow the outer
  iterate without breaking it.
- JFNK: `FREEZE_I_NEXT` plus scalar save/restore in
  `JFNK_RAD_SWEEP_V1` protect production stores during evals; JFNK
  admission requires `USE_LAM_ES=T`, under which the MOM solver
  does not run, so `J_NEXT_FRQ` is never touched by evals. Coupled
  evals use a *frozen* lagged store - the Newton residual is a
  lagged-J operator, not the fully-implicit coupled one (documented
  caveat).

### 5. Defects found (DONE)

1. `FG_RF_DCHG` is not reset at sweep entry (`INIT` block,
   fg file ~2256-2309 resets `FG_SW_*`, `FG_GS_*`, `FG_RF_CAPPED`
   but omits `FG_RF_DCHG`). On an iteration where the refresh does
   not run, `GET_FG_GS_MATURITY` returns the previous iteration's
   retained delta. `MOM_RF_DCHG` is reset correctly (mom file 1672).
   Severity: low-medium; fixed in task 7.
2. Headline metric semantic switching: `FG_GS_MAXCHG` is the sweep
   drift before the refresh and the refresh residual after it
   (overwritten at fg ~6353). Any user monitoring one number sees a
   series whose definition changes mid-iteration; on non-refresh
   iterations it stays the sweep drift. Mitigation: the Phase-7
   diagnostic prints all three quantities separately.
3. Normalization asymmetry: the FG in-sweep denominator floors at
   `1e-4` of the *fresh* field peak; MOM floors at `1e-4` of the
   *stored* column peak. Both are defensible; they are not
   identical metrics and should not be compared cell-for-cell.

No defect was found in the comparison buffers themselves, the
update order, or the commit path.

### 6. Add the field-drift diagnostic mode (DONE - task 7)

### 7. Implement: GS_DRIFT_DIAG VADAT + FG_RF_DCHG reset fix (DONE)

Changed:

- `new_main/mod_subs/control_variable_mod.f`: new
  `LOGICAL GS_DRIFT_DIAG` (default F) with the standard
  declaration/comment block.
- `new_main/mod_subs/rd_control_variables.f`: `RD_STORE_LOG`
  registration for `GS_DRIFT_DIAG`.
- `plane/fg_j_cmf_nonmono_v1.f`: `FG_RF_DCHG=0.0_LDP` added to the
  per-iteration INIT reset block (fix for defect 5.1); new
  `GET_FG_RELAX` accessor exposing `INEXT_RELAX`,
  `KN_RF_RESID_PRV`, and the argmax indices.
- `plane/mom_j_cmf_nonmono_v1.f`: new `GET_MOM_RELAX` accessor
  exposing `JNEXT_RELAX`, `MOM_RF_RESID_PRV`, `MOM_RF_RAN`, and the
  argmax indices.
- `new_main/mod_subs/solve_for_pops.f`: gated
  `IF(GS_DRIFT_DIAG)` diagnostic block after the maturity getters,
  printing per-iteration raw sweep drift, refresh residual,
  retained (relaxed) delta, relaxation factors, census counts,
  argmax cell, and freeze/release state for both the ray and
  moment fields.

Verification: compile via unit-test Fortran build; static
regression test `unit_tests/python/test_gs_drift_diag.py`.

### 8. Add unit tests (DONE - static tier)

`unit_tests/python/test_gs_drift_diag.py` asserts the flag is
declared with default F, parsed via `RD_STORE_LOG`, the accessors
exist and export the right module state, the diagnostic block
reports raw vs relaxed quantities and the argmax cell, and the
`FG_RF_DCHG` reset is present in the INIT block.

The ten synthetic numeric cases (exact/geometric/over-damped/
oscillatory/zero-field/shock-zone/relaxation/stale-buffer/restart/
non-finite) cannot drive the embedded solver loops without
extracting the metric into a pure routine; that refactor is
deferred (audit section 9, proposed correction P2) so the diff
stays minimal. Analytic expectations for each case are recorded
in the audit so a future driver has oracles.

### 9. Validation (PARTIAL - PASS WITH LIMITATIONS)

Run and results:

- `python python/test_gs_drift_diag.py` - 16/16 PASS (new static test).
- `python python/test_vadat_keys.py` - 73/73 PASS; compiles
  `rd_control_variables.f` + `control_variable_mod.f` for real,
  so both edits are compile-verified.
- `python python/test_knop_ray.py` - 10/10 PASS; compiles
  `fg_j_cmf_nonmono_v1.f` and drives `FG_KNOP_REFRESH_V1` and the
  FGINEXT_GS checkpoint round-trip end-to-end.
- `python python/test_mom_solve.py` - 93/93 PASS; compiles
  `mom_j_cmf_nonmono_v1.f` including checkpoint and store paths.

Not run:

- Full `python unit_tests/run_all_tests.py` aggregate (the
  exercised suites above all pass individually).
- Compile of `solve_for_pops.f` - no unit test builds it; needs
  the full `./build.sh` (or an incremental make) which is a
  session-scale job. The edit uses the same fixed-form `BLOCK`
  construct already used in `jfnk_solve_corr_v1.f`, and the two
  accessor calls are implicit-EXTERNAL like the neighbouring
  `GET_*_GS_MATURITY` calls; residual risk is syntax-only.
- Model-level validation (monotonic / reversal / strong
  compression / multi-reversal / IIn shock) requires
  interactive-partition sbatch jobs per AGENTS.md.

### 10. Documentation (DONE)

`docs/audits/FIELD_DRIFT_AUDIT.md` written.
`cmfgen_obsidian/control/Field maturity and GS freeze.md` gained a
`GS_DRIFT_DIAG` section (VADAT keys and diagnostic output are
required-note changes per AGENTS.md).

## Open items

- [ ] Compile `solve_for_pops.f` via `./build.sh` or incremental
      make to close the syntax check on the diagnostic block.
- [ ] Numeric driver tests for the ten synthetic cases (needs the
      pure-metric extraction P2 or a refresh-harness driver).
- [ ] sbatch validation runs across the five velocity profiles with
      `GS_DRIFT_DIAG=T`, logging the full channel set in audit
      section 10.
- [ ] Decide whether the headline "largest drift" print should be
      relabeled raw-vs-residual (P1, cosmetic; not done).
