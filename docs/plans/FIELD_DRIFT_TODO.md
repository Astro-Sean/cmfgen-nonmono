# FIELD_DRIFT_TODO - radiation-field drift audit task list

Sequential task list for the field-drift / radiation-field convergence
audit. One task at a time; a task is complete only when the listed
verification ran, not merely when code compiles.

> **Supersession note (10-Nov-2026):** this TODO documents the
> 03-Nov-2026 audit cycle (`docs/audits/FIELD_DRIFT_AUDIT.md`). The
> authoritative current audit is
> `docs/audits/FG_FIELD_DRIFT_AUDIT.md`, which covers the
> 10-Nov-2026 critical-source-audit and fail-safe cycles
> (`7220b616`->`e6e8f47a`, `e6e8f47a`->`280c4a44`). Findings below
> that the newer audit re-classified are marked inline; the
> historical record is preserved.

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

Counts below are the state at this audit stage (03-Nov-2026);
later stages are recorded per revision in the reconciliation
section and in FG_FIELD_DRIFT_AUDIT.md.

Run and results:

- `python python/test_gs_drift_diag.py` - 16/16 PASS (new static test).
- `python python/test_vadat_keys.py` - 73/73 PASS; compiles
  `rd_control_variables.f` + `control_variable_mod.f` for real,
  so both edits are compile-verified.
- `python python/test_knop_ray.py` - 10/10 PASS at this stage
  (now 41/41 after the 10-Nov-2026 cycles).
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

`docs/audits/FIELD_DRIFT_AUDIT.md` written (the 03-Nov-2026 audit;
the 10-Nov-2026 critical-source cycle is in
`docs/audits/FG_FIELD_DRIFT_AUDIT.md`, now the authoritative doc).
`cmfgen_obsidian/control/Field maturity and GS freeze.md` gained a
`GS_DRIFT_DIAG` section (VADAT keys and diagnostic output are
required-note changes per AGENTS.md).

## Open items

Status below reconciled against FG_FIELD_DRIFT_AUDIT.md
(10-Nov-2026 cycles, HEAD `280c4a44` at time of writing).

- [x] Compile `solve_for_pops.f` - DONE: clean object builds
      throughout the 10-Nov-2026 cycles (solve_for_pops.o archived
      after each edit).
- [~] Numeric driver tests for the ten synthetic cases - PARTIAL:
      test_knop_ray.py now covers 41 cases incl. refresh/abort/
      parity probes; the explicit zero-seed, perturbed-seed and
      multi-refresh store-convergence cases remain listed as gaps
      in FG_FIELD_DRIFT_AUDIT.md Pillar 7.
- [ ] sbatch validation runs - STILL OPEN, now the dominant owed
      item: model-level run_model_grid verification of the
      Nov-10 fixes (field residuals, sweep drift, retained/pending
      deltas, rejected/diverged counts, population changes, JFNK
      residual history, emergent spectrum), plus 1-thread vs
      multi-thread reproducibility. Required before the drift
      question can be called closed.
- [ ] Decide whether the headline "largest drift" print should be
      relabeled raw-vs-residual (P1, cosmetic; not done).

## Reconciliation with FG_FIELD_DRIFT_AUDIT.md (10-Nov-2026)

The newer audit re-derived the same machinery on `7220b616` and
labels the KN_RES_PASS findings differently:

1. Task 4's "adaptive-omega category mismatch" (inner-map ratio
   feeding the outer store write) corresponds to the newer
   audit's **A2 - modelling assumption, not a proven defect**: the
   estimator now measures a well-defined inner-map dominant-mode
   rate; whether it proxies the outer map's slow mode is assumed,
   not proven. No inner/outer map equivalence is claimed anywhere.
2. The KN_RES_PASS estimator defects in the newer audit:
   - **A1** (FIXED 10-Nov-2026): successive-slot maxima compared
     different ray/direction subsets; replaced by per-ray same-
     direction `RHO_WIN` rings -> `RHO_M/RHO_P` -> pooled P90 ->
     `INEXT_RELAX`. `KN_RES_PASS` is telemetry-only now.
   - **A3** (FIXED 10-Nov-2026): `KN_RES_PASS` updated before the
     divergence check could fold NaN/huge residuals in; the update
     now sits below the divergence check in both sweeps.
3. Task 5's defect list (FG_RF_DCHG reset, headline metric
   switching, normalization asymmetry) is from this cycle's own
   audit and is preserved as historical record; the newer audit's
   ranked findings (A1-A3, C1-C4, B1, D1) supersede it for current
   status.
4. Test counts: this TODO's figures (10/10, 73/73, 16/16, 93/93)
   are the 03-Nov-2026 stage. Current: test_knop_ray.py 41/41,
   test_iter_ledger.py 31/31, test_mom_aa_state.py 62/62,
   test_vadat_keys.py 85/85 - see FG_FIELD_DRIFT_AUDIT.md and
   AUDIT_LOG.md for per-revision detail.

## Phase 5.3 cycle - P-direction false convergence (14-Nov-2026)

Continuation of the drift audit into the ordered ray sweeps. Full
record in `docs/audits/PHASE_5_3_P_DIRECTION_AUDIT.md`; the fix is
committed on `main` (`0ef42d05` and the companion entries).

### What the cycle found

1. The `.NOT. HAS_CP` fast exits in `SWEEP_FLD` (one pass) and
   `SWEEP_FLD_ORD` (two passes) measured the pass residual and then
   overwrote `RES`/`RES_FREQ` with 0 before exiting `IXIT=0`. The
   commit gate treats that as converged, so a ray whose ordered
   march still moved was committed mid-iteration. Demonstrated
   defect, not a reporting quirk.
2. The lagged C-corner terms (Knop Paper II Eqs 31-32) keep a
   genuine cross-pass dependency even with no compressing cells,
   so no single pass can be exact on a real stencil. At `J=NI-1`
   the DSU slope also reads the ring-aliased row-1 slot, so the
   lagged corner reaches both ray ends.
3. Ray 82 P-direction (m15 frz82 arm): pass-1 `R_MAP=0.482`
   argmax `(3,17881)`, pass-2 `R_MAP=0.621` argmax `(3,27360)`,
   ~27K above-tolerance cells per pass - then forced clean under
   the old exit. The corner cell alternates `3.18e-8 <-> 1.90e-8`
   across passes and refreshes, a two-state orbit pinned to ~0.3%
   - consistent with an intrinsic period-2 mode of the ordered
   map, coefficient-insensitive (frozen 0.62 vs live 0.68).
4. Store churn (`sw` 0.84 frozen vs 0.27 live) is the coefficient-
   stream divergence, a separate quantity from the corner mode.
5. Manufactured ring-faithful stencil: the map converges to a
   unique seed-independent fixed point that satisfies the
   simultaneous stencil residual at roundoff - the machinery is
   self-consistent, the early exit was the defect. Aitken and
   under-relaxation break a genuine two-cycle and reach the
   underlying fixed point; shallow Anderson does not. Both
   accelerators sit behind `NPO>=3`/`NPO>=4`, unreachable under
   the old two-pass exit.

### What was changed (5.3R, approved and implemented)

Both forced-clean `.NOT. HAS_CP` exit blocks deleted from
`plane/fg_j_cmf_nonmono_v1.f`. `RES`/`RES_FREQ` are always
measured; `IXIT=0` is retired; every success path now requires
`RES<KN_RAY_TOL` or the floored bulk exit. Maximum-pass and
stall/plateau/cap semantics unchanged; `HAS_CP=T` rays untouched.
No tolerances, stencil coefficients, ring semantics, donor
indexing, store write-back, `INEXT_RELAX`, `OMEGA_RAY`,
accelerator mathematics or population equations were changed.
Tests: `kn_psweep_source_test.py` (25 pins), `kn_ccorner_test.py`
(manufactured stencil + exit-policy + period-2 classifier),
`test_knop_lagged_field.py` re-pinned; KNOP/moment/JFNK/OBS
regressions all pass. Debug exe `c3d8b584`; production exe
`f52d8125` untouched.

### Outstanding issues

- [~] Production validation in flight: `p53s_fix82` (13158010,
  corrected exe, frz82-clone VADAT) and `p53q_probe82` (13157367)
  are queued behind `AssocMaxJobsLimit`. Needed answers: does
  ray-82 P reach stall/plateau/maxpass (`IXIT` 3/4/6), does the
  corner orbit persist as a true two-cycle at pass 3+ (d1 finite,
  d2 small), do the verdicts route to nc-reject or capped commit,
  and does `fg_rf`/`KN_RF_RESID` honestly show the ~0.6 residual.
- [ ] Ray 83 under measured exits: floor-dominated cells should
  now return a real residual; confirm the convergence
  classification stays separate from the jump-veto path and the
  ray cannot reach the converged `WRITE_STORE`.
- [ ] `KN_HOT_CNT` should now accumulate the L~17881 column;
  confirm the hot-window machinery picks it up as intended.
- [ ] Maturity-gate timing: models with corner-oscillating rays
  will report honest (larger) `fg_rf` and may hold the gate open
  longer - expected consequence of truthful reporting, flag in
  run notes.
- [ ] Corner-mode acceleration decision deferred: guarded Aitken
  or `OMEGA_RAY<1` break the manufactured two-cycle, shallow
  Anderson does not. Selection waits on the production pass-3+
  evidence; no accelerator was enabled in 5.3R.
- [ ] Whether the ray-83 floor cells and the ray-82 corner share
  a mechanism or need separate handling remains open.
- [ ] Earlier open items still stand: sbatch model-level
  validation runs (dominant owed item), headline-drift relabel
  (cosmetic), synthetic-case numeric drivers beyond the
  manufactured stencil test added in 5.3Q.
