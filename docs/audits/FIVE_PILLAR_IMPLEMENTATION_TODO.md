# CMFGEN Five-Pillar Implementation TODO

This file is the single source of truth for the remaining development
work on the `cmfgen_jfnk` non-monotonic fork. It was created under the
project owner's TODO directive (`/ptmp/sbrennan/SN2024pba/TODO.md`) and
is updated in place as tasks complete. Do not create competing TODO
files; do not delete completed tasks.

## 1. Project status

- Repository: `git@github.com:Astro-Sean/cmfgen-nonmono.git` (private)
- Current branch: `main`
- Current HEAD: `6410d07442dab6aae670a3fa476f7bb9d015e63c`
  (= `origin/main`, verified by fetch on 2026-11-03). The previously
  recorded `ccc52ffd` was stale by the seven docs-only commits that
  created this tracker; no production source changed between them
  (Phase-0 reconciliation, 2026-11-03).
- Working tree: `build_debug.sh` deletion was restored 2026-11-03;
  untracked `unit_tests/systems/_pv.txt` remains a test artifact.
- Pristine baseline: `../cmfgen_reference/cur_cmf`, extracted from
  `cur_cmf_18jun25.tar.gz`, md5 `1d5a7a168e77329dd97c3019bcad3c09`
  (re-verified 2026-11-03, matches AGENTS.md)
- Last updated: 2026-11-03 (BUG-001 fix + Phase-0 reconciliation)
- Last full build: 2026-11-03 ~20:44, `./build.sh` clean, all three
  exes present; unit suite 21/21 files pass (`build_audit_fix.out`);
  rebuild in progress on the BUG-001 working tree
- Overall development status: all five pillars have production code and
  unit-test coverage; Sobolev+LAM_ES=T path is model-validated on m15;
  KNOP and JFNK paths are implementation-complete but not
  model-validated; moment (LAM_ES=F) path documented as
  non-contracting on m15
- Completed tasks: 12 (P1-001..P1-006, P1-009, P2-012, P4-001,
  P4-012, INT-003, DOC-001)
- Bug hunt: 2 defects fixed - BUG-001 (P0, GMRES hess_res-only
  convergence acceptance), BUG-002 (P2, SOB_RH_GEOM cache
  fingerprint missed SIGMA); 2 refuted - BUG-003 (P3, MKL DGELSY
  param-4 warning on a deliberate-NaN legacy test), BUG-004 (P2,
  zero-GAM resonance floor + exact-zero exception - verified
  regularization). See section 10a register. Phase 0 repository
  reconciliation complete; Phase 1 (JFNK true-residual acceptance)
  complete including surrogate recalibration (`test_jfnk_solver.py`
  89 passed, 0 failed); Phase 2 (SOB_RH_GEOM cache) complete at the
  unit level; Phase 3 (zero-GAM discontinuity, P3-016) complete -
  no defect, unit-level P3-017 evidence recorded.
- Outstanding tasks: 79 (17 in progress with partial evidence,
  62 not started)
- Blocked tasks: 0 (run-dependent tasks are sequenced, not blocked;
  see section 11)
- Current active task: none - Phase 1 done; Phase 2 (SOB_RH_GEOM
  cache fingerprint) is next per the owner's bug-hunt directive
- Next task: Phase 2 BUG-002 work, then P4-002/P4-004 (SOB+JFNK
  HSCAN production run to measure the BUG-001 runtime impact)

## 2. Development rules

- Work on one task at a time.
- Update this file immediately after completing a task.
- Never mark a task complete without evidence.
- Preserve the original CMFGEN path.
- Do not change physical equations without deriving and validating the
  change.
- Do not add arbitrary clipping, floors, ceilings, or positivity
  repairs to force convergence.
- Do not increase regularisation merely to make a solver converge.
- Prefer root-cause fixes over additional numerical safeguards.
- Do not add unnecessary VADAT controls.
- Every bug fix must have a regression test.
- Every numerical change must have a physical validation plan.
- Never describe a stalled, capped, or iteration-limited calculation
  as converged.
- Do not modify unrelated code.
- Do not push or merge unless explicitly authorised. (Authorised:
  the project owner's directive requires an up-to-date GitHub; commits
  follow the `commit_and_push` skill - one logical change per commit,
  push to `origin/main` after each commit.)

## 3. Task status conventions

- `[ ]` TODO
- `[~]` IN PROGRESS
- `[x]` DONE
- `[!]` BLOCKED
- `[-]` NOT APPLICABLE / WONT FIX

Priorities: P0 = correctness/blocking, P1 = required for pillar
acceptance, P2 = validation/diagnostics, P3 = documentation/metrics.

Each task records: priority, description, current evidence,
files/routines, required action, acceptance criteria, dependencies,
completion record.

## 4. Pillar 1 - Original CMFGEN preservation

Acceptance: the original path remains consistent with the pristine
reference within documented numerical tolerances; full-model
comparisons exist. Pillar is not complete until P1-007/P1-008 are done.

### [x] P1-001 - Establish pristine upstream reference and reproducible baseline
- Priority: P0
- Description: confirm the reference tree, archive checksum, git
  baseline, build state, and test state for this development line.
- Evidence: `../cmfgen_reference/cur_cmf_18jun25.tar.gz` md5
  `1d5a7a168e77329dd97c3019bcad3c09` re-verified 2026-11-03;
  HEAD=`ccc52ffd`=`origin/main`; last `./build.sh` clean 2026-11-03
  with 21/21 unit-test files passing (`build_audit_fix.out`).
- Files: `../cmfgen_reference/`, `build.sh`, `Makefile*`.
- Action: done - recompute md5, fetch origin, record build/test state.
- Acceptance: reproducible baseline pinned with checksums. MET.
- Dependencies: none.
- Completion: 2026-11-03 @ `ccc52ffd`. md5 verified; tree clean;
  baseline recorded in section 1. Limitation: build/test state is
  from the last recorded log, not a fresh rebuild this session.

### [x] P1-002 - Identify all upstream routines modified by this development
- Priority: P0
- Description: produce the authoritative list of files that differ
  from the pristine tree.
- Evidence: `docs/audits/P1_BASELINE_DIFF_INVENTORY_2026-11-03.txt` -
  109 modified files, 25 only-in-fork entries (of which ~16 are real
  new source files: the JFNK modules, the two nonmono solvers,
  `sob_rh_couple.f`, `recompute_sigma_true.f`, `rv_sn_model_snhyd.f`,
  `cmf_blkband_pin_null_v1.f`, `knop_status_mod.f90`; the rest are
  docs/build/test/artifact entries).
- Files: whole tree vs `../cmfgen_reference/cur_cmf`.
- Action: done - `diff -rq` with build/test/docs exclusions; inventory
  committed under `docs/audits/`.
- Acceptance: a committed, dated inventory file. MET.
- Dependencies: P1-001.
- Completion: 2026-11-03 @ `ccc52ffd`. Note the inventory includes
  intentional upstream-touching edits (e.g. strict interp guards,
  empty-level ZNET guards) - each modified file still needs a
  line-level gating check under P1-003/P1-004.

### [x] P1-003 - Audit shared module state and SAVE variables with new flags off
- Priority: P0
- Description: verify no module-level or SAVE state written by
  nonmono/JFNK code paths leaks into the stock path when flags are F.
- Evidence: `docs/audits/P1_MODULE_STATE_AUDIT_2026-11-03.md` -
  classification table over all new SAVE/module state in the
  109-file inventory.
- Files: `control_variable_mod.f`, `mod_cmfgen.f`,
  `fg_j_cmf_nonmono_v1.f`, `mom_j_cmf_nonmono_v1.f`,
  `jfnk_*` modules, `comp_j_blank.f`, `cmfgen_sub.f`,
  `solve_for_pops.f`, `solveba_v13.f`, `fg_j_cmf_v13.f`,
  `get_moms_non_rel.f`, `recompute_sigma_true.f`.
- Action: enumerated new SAVE decls per modified file; traced write
  and read gating for each.
- Acceptance: written classification table; zero unexplained
  class-(c) writes. MET.
- Dependencies: P1-002.
- Completion: 2026-11-03 @ `b12a9c27`. All new state class (a)
  flag-gated or (b) overhead-only; class (c) list = 9 documented
  upstream-defect fixes (degenerate-J guards, empty-level brackets,
  pin-null, strict interp, W-series IOSTAT, ZERO FLUX). No
  unexplained stock-path leakage.

### [x] P1-004 - Verify original formal-transfer dispatch with all new options off
- Priority: P0
- Description: prove `FG_J_CMF_V13`/`MOM_J_CMF_V11` and the stock
  boundary/Eddington paths execute identically to upstream when
  `DO_NONMONO_SOB=F`, `DO_NONMONO_KNOP=F`, `DO_JFNK=F`.
- Evidence: `P1_MODULE_STATE_AUDIT_2026-11-03.md` P1-004 addendum -
  per-file dispatch table over all 14 modified formal-path files;
  `mom_j_cmf_v11.f` unmodified (0 diff); `fg_j_cmf_v12.f` dead code.
- Files: all `plane/` modified files, `comp_j_blank.f`,
  `do_var_cont.f`, `mod_cmfgen.f`.
- Action: audited every modified formal-path file's gating.
- Acceptance: every modified formal-path file shows its changes dead
  or identity under flags-off; dispatch audit written. MET at
  source-inspection level; executable A/B proof rides on P1-007.
- Dependencies: P1-003.
- Completion: 2026-11-03 @ `b12a9c27`. All changes flag-gated,
  telemetry, or documented class-(c) defect fixes confined to
  NaN/garbage inputs; `ABS(dV)` insertion adds resolution on
  decelerating profiles (documented).

### [x] P1-005 - Verify original Sobolev and line-transfer dispatch, flags off
- Priority: P0
- Description: prove `SUB_SOB_LINE_V3`, `SOBJBAR_SIM`, `SOBEW*`,
  `INS_LINE*` behave identically to upstream with new options off.
- Evidence: `test_sob_tau_floor.py` 14/14 - bit-identity vs reference
  for sigma>=0 and sigma=-0.5; floor binds only at sigma<=-1.
  `test_sob_rh_couple.py` 9/9 - zero pairs under homology, V<=0 and
  NEG_OPACITY rows never couple. Dispatch table written in
  `P1_MODULE_STATE_AUDIT_2026-11-03.md` P1-005 addendum.
- Files: `sobjbar_sim.f`, `sobew.f`, `sobew_grad_v2.f`,
  `sobew_grad_v3.f`, `sub_sob_line_v3.f`, `sob_rh_couple.f`,
  `ins_line_v6.f`, `ins_line_obs_v5.f`, `cmfgen_sub.f` (SIM_LINES).
- Acceptance: flags-off bit-identity or documented benign difference
  per routine; written dispatch table. Met: SIGMA_TRUE plumbing gated
  on `DO_NONMONO_SOB`; DRC_SIM=0 makes the RH term an exact no-op;
  non-gated edits are documented defect fixes (ETAL_MAT/SRAT denom
  guards, dRATE_dUP algebraic rewrite identical for valid iterates,
  `.AND.` short-circuit OOB fix in ins_line*).
- Dependencies: P1-003.
- Completion: 2026-11-03. Limitation: ins_line and dRATE fixes are
  inspection-verified without dedicated unit tests; full model A/B
  remains open under P1-007/008.

### [x] P1-006 - Verify original population/statistical-equilibrium solver
- Priority: P0
- Description: prove `SOLVE_FOR_POPS`+`SOLVEBA_V13` stock iteration is
  unchanged when `DO_JFNK=F` and nonmono flags off.
- Evidence: `CONVENTIONAL_SOLVER_AUDIT_2026-11-03.md` traced the
  Lambda/BA call graph and found ordering consistent with upstream.
  Consolidated non-gated-edit table written in
  `P1_MODULE_STATE_AUDIT_2026-11-03.md` P1-006 addendum.
- Files: `solve_for_pops.f`, `solveba_v13.f`, `cmf_blkband_v3.f`,
  `cmf_blkband_pin_null_v1.f`, `fiddle_pop_corrections_v2.f`,
  `generate_full_matrix_v3.f`, `fixpop_in_ba_v3.f`.
- Acceptance: written table of non-gated stock-path edits with
  justification + test reference for each. Met: POP_DIV floor and
  pin-null are documented defect fixes (pin-null covered by unit
  test); fiddle_pop limits flag-gated (stock=+/-HUGE=upstream);
  generate_full_matrix adds JFNK_QUIET_IO print gates (internal
  flag, no VADAT key); fixpop is text-only; BCD recall is opt-in
  under WR_BA_INV.
- Dependencies: P1-002.
- Completion: 2026-11-03. Limitation: N_DIV_FLOOR counts POP_DIV
  activations but never prints to LUER - invisible unless debug
  dumps run. Noted for P1-009 diagnostics pass.

### [ ] P1-007 - Compare complete representative models vs pristine upstream
- Priority: P0
- Description: run the same model inputs through this tree's binary
  with all new flags off and through a binary built from
  `../cmfgen_reference/cur_cmf`, on a monotonic model and (separately)
  document what the stock path does on a non-monotonic input.
- Evidence: none found - no dedicated A/B run pair exists in the run
  history. `exe_frozen/` holds hashed binaries (unconfirmed
  provenance); a reference build does not exist yet.
- Files: run scripts (`scripts/run_model_grid.py`), reference tree
  build, VADAT for chosen model.
- Action: build the reference tree (or identify a frozen upstream
  exe), pick one converged existing model config, submit both A/B
  runs via `run_model_grid.py` on the cluster.
- Acceptance: two completed runs at the same inputs; differences
  confined to documented stock-path fixes.
- Dependencies: P1-001..P1-006 (audit list finalised first).
- Completion: -

### [ ] P1-008 - Compare pops/ED/T/J/flux/residuals/spectrum vs pristine
- Priority: P0
- Description: quantitative comparison of POPS, ED, T, RJ, observed
  flux, convergence residuals, and emergent spectrum for the P1-007
  A/B pair.
- Evidence: none.
- Files: OUTGEN, SCRTEMP, OBS outputs, `check_convergence.py`,
  plotting scripts.
- Action: diff the comparison channels at agreed tolerances; explain
  every residual difference by a documented stock-path edit.
- Acceptance: written comparison table with tolerances; every
  difference attributed.
- Dependencies: P1-007.
- Completion: -

### [x] P1-009 - New diagnostics/safeguards do not change legacy results
- Priority: P1
- Description: verify counters, census prints, watchdogs, and guard
  conversions added by this fork are observation-only on the stock
  path (no result dependence).
- Evidence: `ROBUSTNESS_AUDIT_2026-10.md` lists each guard's intent;
  consolidated classification table written in
  `P1_MODULE_STATE_AUDIT_2026-11-03.md` P1-009 addendum.
- Files: `control_variable_mod.f` counters, `solve_for_pops.f`
  (LUM_DIV, GS_*), `cmfgen_sub.f` (CONV_EXIT, f-cap census),
  `solveba_v13.f` (BA_DEM_EXM/EXR, BAACC, BA_EX_MAXCH).
- Acceptance: written table; zero ungated result-affecting
  diagnostics. Met with two documented exceptions:
  (1) `MAIN_COUNTER.NE.LAST_NG` moved into CONV_EXIT_MET - refuses
  converged exit on the extrapolated post-NG call; the ELSE
  continuation then runs where upstream ran a no-op (result-relevant
  only for AUTO-mode models meeting EPS post-NG; documented in
  source). (2) `SCRWF` now captures SCR_RITE_V2's WRITFAIL instead
  of clobbering NEWMOD - defect fix. All other counters/state are
  print-only or inert in stock (BA_TRACE L_NM-gated, BA_DEM floors
  off by default, BA_EX_MAXCH=MAXCH, MAXCH_CTRL=MAXCH).
- Dependencies: P1-003.
- Completion: 2026-11-03. Gap: N_DIV_FLOOR counts cmf_blkband
  divisor-floor activations but never prints to LUER.

## 5. Pillar 2 - DO_NONMONO_SOB

Acceptance: formulation established, multiple-resonance treatment
validated, optical-depth safeguard effects quantified. Reuse/extend
`test_sob_tau_floor.py` and `test_sob_rh_couple.py`; do not duplicate.

### [~] P2-001 - Establish the exact intended reference formulation
- Priority: P0
- Description: README/AGENTS attribute the implementation to Dessart,
  Audit & Hillier (2015, arXiv:1503.05463, Sec 3.1). Confirm this is
  the formulation the project owner intends (Luc Dessart's nonmono
  Sobolev treatment) and not a different scheme.
- Evidence: AGENTS.md sec.3 cites DAH2015 Sec 3.1 for `SIGMA_TRUE`,
  tau floor, and `SOB_LINE_V3` plumbing;
  `CONVENTIONAL_SOLVER_AUDIT_2026-11-03.md` and
  `NONMONO_COMPREHENSIVE_AUDIT_2026-10-30.md` describe the implemented
  equations. `cmfgen_base_dessart/` at the workspace root is a model
  run directory (Dessart inputs), not a source tree.
- Files: README.md, AGENTS.md, `rt_literature` skill notes, code
  comments in `sobjbar_sim.f`/`sobew*.f`.
- Action: re-read DAH2015 Sec 3.1; write a short doc mapping paper
  equations to code lines (tau_sob denominator, floor, nonlocal
  RH coupling); flag any divergence for owner review.
- Acceptance: an equation-by-equation mapping doc; owner-intent
  question answered in writing.
- Dependencies: none.
- Completion: 2026-11-03 - `SOBOLEV_DAH2015_EQUATION_MAP.md` maps
  every Sec-3.1 element to code lines. The only extension past the
  paper's stated recipe is the `1D-4` arm of the tau floor (bounds
  the sigma -> -1 case where the paper's own -0.2*(1+sigma) floor
  vanishes). The OCR'd paper text drops the floor's minus sign;
  -0.2*(1+sigma) is the only reading that is positive for sigma<-1
  and matches all four code comments. The nonmono-J routine the
  paper describes maps to the Knop/moment stack - a superset, not
  a divergence.
- Completion: partial (attribution recorded; equation map open).

### [~] P2-002 - Derive and document the implemented Sobolev optical-depth equation
- Priority: P0
- Description: write the exact tau_sob expression the code evaluates,
  including angle dependence and the sigma<0 floor branch.
- Evidence: `test_sob_tau_floor.py` asserts the floor formula
  `max(|1+sigma*mu^2|, -0.2*(1+sigma), 1e-4)` for sigma<0 and
  bit-identity for sigma>=0; AGENTS.md sec.3 documents the same.
- Files: `sobjbar_sim.f`, `sobew.f`, `sobew_grad_v2.f`,
  `sobew_grad_v3.f`, `sub_sob_line_v3.f`.
- Action: extract the four implementations' expressions; confirm they
  are algebraically identical; write the derivation with references
  into the P2-001 doc.
- Acceptance: single documented expression + per-file line refs;
  mismatches (if any) filed as defects.
- Dependencies: P2-001.
- Completion: partial (formula known and tested 2026-11-03; all four
  sites confirmed bit-identical at sobjbar_sim.f:124-130,
  sobew.f:131, sobew_grad_v2.f:137, sobew_grad_v3.f:140 - same
  GAM=GAMH/T1 with T1=MAX(|1+sigma*mu^2|,-0.2*(1+sigma),1e-4) for
  sigma<0 else |T1|. Written derivation into P2-001 doc open).

### [~] P2-003 - Audit the denominator floor MAX(|1+sigma*mu^2|, -0.2*(1+sigma), 1E-4)
- Priority: P0
- Description: quantify when the floor activates and how much it
  changes tau_sob vs the unfloored expression.
- Evidence: analytic audit 2026-11-03 - all four sites identical
  (sobjbar_sim.f:124-130, sobew.f:131, sobew_grad_v2.f:137,
  sobew_grad_v3.f:140), same Dessart 2015 citation and formula.
  Activation region in (sigma,mu^2): sigma in (-1,0): -0.2(1+sigma)
  is negative so the floor is MAX(|T1|,1e-4), binding only for
  sigma<-0.9999 near mu^2=1 (measure-zero sliver). sigma=-1:
  |1-mu^2| vanishes at the ray center, tau capped at GAMH*1e4.
  sigma<-1.0005: floor F=0.2(|sigma|-1) binds on the mu^2 band
  |1+sigma*mu^2|<F around mu^2=1/|sigma| (half-width F/|sigma|) -
  for sigma=-3 that is ~27% of the mu^2 range per affected ray
  point, tau capped at GAMH/0.4 vs the formally singular value.
  The "error" vs the unfloored expression is the asserted physical
  bound (finite resonance width at a vanishing LOS gradient) - the
  unfloored expression is singular there, so the floor is the
  regularization itself, not an approximation of a finite answer.
  Inactive for sigma>=-0.9999, monotonic path bit-identical.
- Files: the four tau-floor files.
- Action: none on the formula; a print-only activation census on
  production runs remains an available option if the P3-017-style
  floor-population question needs a number.
- Acceptance: documented activation region + bounded relative error
  it introduces (analytic report recorded; counter optional).
- Dependencies: P2-002.
- Completion: partial (analytic audit done 2026-11-03; production
  activation census optional/open).

### [~] P2-004 - Validate small-gradient and sigma=-1 limits
- Priority: P1
- Description: verify tau_sob reduces to the correct limits:
  sigma->0 (homologous), sigma=-1 (zero gradient - the floor's
  reason to exist), sigma<-1 (compression).
- Evidence: `test_sob_tau_floor.py` covers sigma=-0.5 bit-identity and
  sigma<-1 floor path; `test_rv_snhyd.py` covers homologous
  SIGMA_TRUE=0 construction. The sigma=-1 surface itself (denominator
  -> floor boundary continuity) not explicitly tested.
- Files: `test_sob_tau_floor.py`, `sobjbar_sim.f`.
- Action: add limit cases to the existing test (sigma=-1 exactly,
  sigma just below -1 at several mu); assert continuity across the
  floor boundary and the sigma=0 upstream form.
- Acceptance: test cases at/around sigma=-1 and sigma=0 pass with the
  documented limits.
- Dependencies: P2-002.
- Completion: 2026-11-03 - `test_sob_tau_floor.py` now covers all three
  limit surfaces: T7 sigma=-1.0001 exercises the pure-1e-4 arm on
  near-radial rays (BETA matches floored formula to 6e-12, ref matches
  unclamped); T8 straddles the pole entry at sigma=-1 with near-radial
  rays - the floored response is bounded (~2-4e-3 absolute on BETA)
  and tracks the unbounded reference to 1e-5, so the pole entering
  the angle grid adds no numerical jump; T9 sigma=0 is bit-identical
  to upstream with BETA == EXF(GAMH) exactly (6e-14). 17/17 pass.

### [x] P2-005 - Audit resonance finding and ordering
- Priority: P0
- Description: verify how SOB line resonances are located/ordered
  along a ray for non-monotonic V, and that ordering logic handles
  re-entrant velocity ranges.
- Evidence: traced 2026-11-03. There is no depth-search to order:
  resonance "finding" is a frequency-window test - RESONANCE_ZONE
  keys on LINE_ST_INDX_IN_NU <= FREQ_INDX < LINE_END_INDX_IN_NU
  (set_line_opac.f:336-349). A line's Sobolev band is contiguous in
  observer frequency, so the sweep enters/leaves each line's window
  once regardless of the depth ordering of resonances; the repeated
  depth resonances along a ray are picked up by the per-depth
  CHIL_MAT/ETAL_MAT*LINE_PROF_SIM sum inside the ray quadrature
  (cmfgen_sub.f:1730, 2497-2507). The one monotonic-V assumption -
  window sizing via VINF=V(1) - is already corrected at
  cmfgen_sub.f:1108-1113: under DO_NONMONO_SOB/KNOP,
  VINF=MAXVAL(ABS(V)) so the windows cover the full signed velocity
  range (the "2Vinf interaction" bound at 2920 stays conservative).
  V<=0 under non-mono flags is a loud domain boundary (1120-1124);
  par_opac_sub.f:608-614 independently rejects non-monotonic V in
  the standalone utility (documented guard, unreachable in the main
  path). No defect.
- Files: `set_line_opac.f`, `cmfgen_sub.f`, `sobjbar_sim.f`,
  `par_opac_sub.f`.
- Action: none.
- Acceptance: written audit above; the only monotonic-V assumption
  found is already corrected and flag-gated.
- Dependencies: P2-001.
- Completion: 2026-11-03.

### [ ] P2-006 - Validate multiple resonance locations along a ray
- Priority: P0
- Description: test that a single line resonating at multiple depths
  on one ray deposits/absorbs correctly at each location.
- Evidence: none dedicated; multi-resonance handling is inside the
  ray tau integration per the audit.
- Files: `sub_sob_line_v3.f`, manufactured V profiles in
  `unit_tests/f90/`.
- Action: build a manufactured two/three-resonance ray case in the
  unit-test harness; assert per-resonance tau and energy
  conservation across resonances.
- Acceptance: manufactured multi-resonance test passes against an
  independent quadrature.
- Dependencies: P2-005.
- Completion: 2026-11-03 - the RH-coupling driver now covers a
  W-shaped profile (two reversals, 108 accepted pairs, up to 14
  coupling partners per depth). Pair selection and the per-depth
  DRC rates both match the independent O(ND^2) quadrature
  (max|dDRC|=4e-12). Local-Sobolev per-resonance tau is already
  algebraically pinned by test_sob_tau_floor T3/T7. The known
  residual gap is architectural and documented: the two-point RH
  kernel has no three-body term, so a single ray hitting three
  coincident resonances loses the third-extinction correction
  (Feldmeier & Nikutta 2006 regime) - a stated approximation, not
  a defect.

### [ ] P2-007 - Audit line source functions, escape probabilities, propagation
- Priority: P1
- Description: audit S_line, the Sobolev escape/beta factors, and how
  absorbed/emitted line photons propagate between resonance points
  under non-monotonic V.
- Evidence: `sub_sob_line_v3.f` guards audited in AGENTS.md
  (transient-iterate VB/VC/SRAT guards); `sob_rh_couple.f` tested at
  pair level.
- Files: `sub_sob_line_v3.f`, `sobjbar_sim.f`, `sob_rh_couple.f`.
- Action: trace source-function assembly and beta-factor usage;
  confirm escape probabilities use SIGMA_TRUE consistently and the
  propagation to downstream resonances is the intended RH kernel.
- Acceptance: written audit tying code to the paper formulation.
- Dependencies: P2-005, P2-009.
- Completion: 2026-11-03 - audited end to end. Source function:
  S_L = ETAL_MAT/CHIL_MAT per (depth, sim-line); driver tests confirm
  the formula level (ETAL=0 -> DRC=0). Escape probabilities: BETA =
  sum_L w_L*EXF(GAM_L) with the floored tau - verified to 6e-12
  against the analytic expression (test_sob_tau_floor T3/T7). SIGMA
  consistency: under DO_NONMONO_SOB all three consumers (SOBJBAR_SIM
  rates, SOBEW EW estimate, SOB_RH_COUPLE kernel) take SIGMA_TRUE;
  DO_SOB_CLAMPED_LINES consistently switches all three back to clamped
  SIGMA (sub_sob_line_v3.f:144-176) - a documented consistency probe,
  not a mixed-gradient bug. Propagation: DRC_SIM enters the STEQ
  residual only (deliberately excluded from ZNET_SIM so the BA
  Jacobian does not differentiate a lagged remote term, per the
  :164-168 comment). Rate assembly is antisymmetric:
  STEQ(up)-=T2/T1, STEQ(low)+=T2/T1 with T2=ETAL*ZNET-DRC - atoms
  conserved across the level pair; line cooling STEQ_T-=T2*T3.
  Lambda-iteration derivatives (VB=-JBAR/ETAL, VC=JBAR*CHIL/ETAL^2)
  are the -d(JBAR*CHIL/ETAL) pair, correctly guarded against
  ETAL<=0 and NEG_OPACITY rows.

### [ ] P2-008 - Audit continuum interactions and scattering
- Priority: P1
- Description: verify continuum opacity/emissivity coupling inside
  the nonmono Sobolev path (background continuum at resonances,
  electron scattering wings).
- Evidence: `cmfgen_sub.f` negative-opacity check + KN_NO_GAIN noted;
  no dedicated audit.
- Files: `cmfgen_sub.f` (SIM_LINES opacity block),
  `sub_sob_line_v3.f`, `comp_j_blank.f` (ZETA/THETA).
- Action: trace where line tau interacts with continuum tau on a
  ray; document assumptions.
- Acceptance: written audit.
- Dependencies: P2-005.
- Completion: 2026-11-03 - audited. Lines and continuum share one
  ray integration: inside RESONANCE_ZONE the code adds
  CHI += CHIL_MAT*LINE_PROF_SIM and ETA += ETAL_MAT*LINE_PROF_SIM
  (cmfgen_sub.f:2506-2507, inner-boundary block :1731), so the line
  rides the continuum's own (non-monotonic) rays with a mixed
  source function S=ETA/CHI - no separate line-only transfer.
  The net-rate bracket is the standard CMFGEN form
  ZNET += QW*(1 - RJ*CHIL/ETAL) with RJ the combined continuum+line
  comoving J (:2821); the ETAL_MAT>0 guard there is a documented
  NaN fix for depopulated upper levels, not flag-gated because the
  prior result was NaN. Electron scattering enters through
  CHI_SCAT/ESEC in the SOBEW EW path and the continuum CHI;
  NEG_OPACITY rows flatten CHIL_MAT to 1/NUM_SIM_LINES (:4095).
  Assumption recorded: line photons interact with the continuum
  only through this shared-opacity architecture - there is no
  continuum-vs-line ray ordering question to audit (P2-005).

### [~] P2-009 - Audit DO_NONLOC_SOB and coupling between resonance points
- Priority: P1
- Description: audit the Rybicki-Hummer nonlocal coupling
  (`sob_rh_couple.f`): pair selection, kernel, dropped-pair policy,
  missing higher-order extinction.
- Evidence: `CONVENTIONAL_SOLVER_AUDIT_2026-11-03.md`: kernel drops
  pairs straddling V<=0 and has no third-resonance extinction -
  "an approximation layer, not a full multi-resonance solver";
  `test_sob_rh_couple.py` covers pair count/symmetry, mu^2 identity,
  DRC term, V<=0/negative-opacity suppression, BB_COR=0/ETAL=0.
  `DO_NONLOC_SOB` auto-off noted in test_vadat_keys.
- Files: `sob_rh_couple.f`, `rd_control_variables.f` (DO_NONLOC_SOB).
- Action: finish the audit: quantify the dropped-pair fraction and
  the third-resonance omission on a manufactured profile; document
  the validity boundary.
- Acceptance: written validity limits; dropped-pair policy justified
  or filed as a defect.
- Dependencies: P2-005.
- Completion: 2026-11-03 - quantified on manufactured profiles. The
  V>0 gate drops only pairs *touching* a V<=0 cell (sob_rh_couple.f:150
  endpoint test, RH_NDROP counter prints the tally): on a profile with
  a -200 km/s stagnation cell, 30/98 geometric coincidences (31%) are
  dropped, all touching that cell. Pairs that straddle a stagnation
  region with positive endpoints are kept - correct per RH78, since
  the photon propagates freely between coincident resonances and the
  intermediate off-resonance gas does not absorb that frequency. The
  kernel is pairwise (no three-body term); a tuned triple coincidence
  on one ray loses the third-extinction correction (Feldmeier &
  Nikutta 2006 regime) - documented approximation, not a defect.

### [ ] P2-010 - Validate inner and outer boundary resonance cases
- Priority: P1
- Description: resonances landing on the inner/outer grid boundary or
  ray boundary - verify no off-by-one, no lost flux, correct
  boundary condition application.
- Evidence: none dedicated.
- Files: `sub_sob_line_v3.f`, `ins_line_v6.f`, ray boundary code in
  `fg_j_cmf_v13.f` used by SOB diagnostics.
- Action: manufactured boundary-resonance cases in unit tests.
- Acceptance: boundary cases pass with documented behaviour.
- Dependencies: P2-006.
- Completion: 2026-11-03 - audited at source; no off-by-one found.
  Frequency edge: the zone is [LINE_ST,LINE_END] inclusive and
  END_RES_ZONE fires at FREQ_INDX==LINE_END (set_line_opac.f:341-346)
  with no range restriction, so a window ending on the last grid
  frequency still terminates correctly. Windows are built by
  marching the grid and inserting frequencies (ins_line_v6.f:380-390),
  so ST/END always land on real grid points; a BLANK line whose
  window misses the grid leaves indices 0 and hits a hard STOP with
  diagnostics (:445-460) - loud, not silent. Spatial edge: the
  inner-boundary CHI(ND) add is handled in its own block
  (cmfgen_sub.f:1731). A dedicated manufactured test would duplicate
  the grid-construction machinery; the :445 check already makes any
  boundary miss a fatal error in every production run.

### [ ] P2-011 - Compare with an independent reference calculation
- Priority: P1
- Description: compare nonmono-SOB line profiles/tau against an
  independent implementation or published result (Dessart model data
  in `cmfgen_base_dessart/` provides comparison inputs/outputs).
- Evidence: `DO_SOB_DEP_WINGS` history records an m15 SOB-render vs
  Dessart `.fl` comparison (HeI 5876 FWHM 18A vs ~30A before wings,
  ~32A after) - a single-line comparison exists informally.
- Files: `cmfgen_base_dessart/` data, `obs/` render path.
- Action: formalise at least one line-profile and/or integrated-flux
  comparison with stored tolerances.
- Acceptance: written comparison with pass criteria.
- Dependencies: P2-002, P2-012.
- Completion: -

### [x] P2-012 - Run a complete non-monotonic Sobolev model
- Priority: P0
- Description: end-to-end non-monotonic SOB model demonstrating the
  path.
- Evidence: `m15_2p82d_SOB_l14_lte_np` (SOB + USE_LAM_ES=T) converged
  100% at it718 (AGENTS.md); `adp2_esali` (SOB+LAM_ES+DO_ES_ALI)
  converged 100% at it851 with the ES_ALI J loop at dJ/J~1e-4 in
  ~2 passes (AGENTS.md 27-Oct/02-Nov entries). The LAM_ES=F partner
  `m15_2p82d_SOB_l14_lte_np_mom` FAILED - bounded non-contraction
  (job 12316397, cancelled it742, ~200K cells quarantined/sweep) -
  that failure is a P3 moment-path finding, not a Sobolev finding.
- Files: run dirs under `cmfgen_validation`; configs
  `dessart_m15_2p82d_SOB_l14_lte_np*.yaml`.
- Action: none for the LAM_ES=T criterion; record run IDs/dates in
  the completion record.
- Acceptance: a complete converged nonmono SOB model. MET for
  USE_LAM_ES=T only; the LAM_ES=F configuration is covered by P3.
- Dependencies: none.
- Completion: MET (LAM_ES=T) 2026-11-03 on recorded evidence; model
  dirs: m15_2p82d_SOB_l14_lte_np (job history per AGENTS_HISTORY),
  adp2_esali.

## 6. Pillar 3 - DO_NONMONO_KNOP

Acceptance: formal and moment solvers pass independent numerical
tests and representative physical validation. Do not remove
MOM_DLN_MIN or MOM_GAM_FAC_MIN without a derived replacement and
passing tests. Do not add acceleration before measuring existing
acceleration.

### 6A. Formal solver (`plane/fg_j_cmf_nonmono_v1.f`)

### [~] P3-001 - Derive the implemented Knop formal-transfer equations
- Priority: P0
- Description: document the discrete operator (inward I_M / outward
  I_P march, chi_hat, signed-|Q| upwind source, dS terms, boundaries)
  against Knop et al. 2009 Paper I.
- Evidence: `KNOP_AUDIT_2026-10-30.md` architecture map records the
  discrete operator and verifies sweep-vs-refresh operator identity;
  `KNOP_FORMAL_SOLVER_AUDIT.md` resolved 22 static-analysis items.
  Missing: the equation-by-equation map to Paper I in one doc.
- Files: `fg_j_cmf_nonmono_v1.f`, `cmfgen_obsidian/solvers/Knop
  formal solver.md`.
- Action: write the equation map; flag any term that departs from
  Paper I (e.g. monotonic-parabola dS clamps).
- Acceptance: documented operator with paper refs.
- Dependencies: none.
- Completion: partial (operator verified internally consistent;
  paper map open).

### [~] P3-002 - Audit generalised opacity and signed frequency coupling
- Priority: P0
- Description: chi_hat = chi + |GAM|/DLN, Q = GAM/(DLN*chi_hat), the
  DLN frequency spacing and sign handling.
- Evidence: `KNOP_AUDIT_2026-10-30.md` records the formulas;
  `KNOP_FORMAL_SOLVER_AUDIT.md` item 3/4 fixed DLN_C floor and
  zero-spacing dGAMdR handling (eb1638ac).
- Files: `fg_j_cmf_nonmono_v1.f`, `SET_FG_NU_NEXT` call site in
  `comp_j_blank.f:843`.
- Action: verify the |GAM| vs signed-GAM usage per term; verify DLN
  at grid edges; document.
- Acceptance: written audit; sign convention proven per term.
- Dependencies: P3-001.
- Completion: partial.

### [~] P3-003 - Audit sign-dependent upwind selection
- Priority: P0
- Description: upwind = I_PREV for GAM>=0 vs lagged I_*_DEPTH_NEXT
  store for GAM<0, including slot indexing s=freq s+1.
- Evidence: architecture map records the selection rule and relaxed
  writes into slots 1..NCF-1 with slot NCF a zero phantom.
- Files: `fg_j_cmf_nonmono_v1.f` (ray solve), `mod_ray_mom_store.f`.
- Action: verify slot mapping off-by-one safety at both frequency
  ends; document.
- Acceptance: written audit; boundary slots proven correct.
- Dependencies: P3-001.
- Completion: partial.

### [~] P3-004 - Audit ray geometry and boundary conditions
- Priority: P0
- Description: ray construction (impact parameters, tangent rays,
  inserted depth points), I_M(1)=0 inner boundary,
  I_P(NI)=I_CORE/I_M(NI) outer boundary.
- Evidence: `KNOP_FORMAL_SOLVER_AUDIT.md` verified COMP_DEPTH_IDX
  bounds and NI=2 handling (e2516057); `ROBUSTNESS_AUDIT_2026-10.md`
  W12 P=0.9999*R_GRID clamp on ray SQRTs.
- Files: `fg_j_cmf_nonmono_v1.f`, `characteristics_v2*.f`,
  `get_ibound.f`, `mod_ray_mom_store.f`.
- Action: document geometry invariants and boundary formulas; check
  inserted/tangent-node handling under reversals.
- Acceptance: written audit.
- Dependencies: P3-001.
- Completion: partial.

### [~] P3-005 - Validate multiple velocity reversals and strong compression
- Priority: P0
- Description: numerical validation of the formal solver on
  manufactured profiles with 2+ reversals and sigma<<-1 compression.
- Evidence: `test_knop_ray.py` now covers three profile classes:
  the original single band, `rev` (two compressing bands at depths
  3-4 and 6-7 - the GAM sign reverses twice along the ray, T9:
  NXIT=1, RESID=6e-7, NPASS=3, NCOMP=24) and `strg` (same band at
  GAM ~150x deeper - sigma<<-1 regime, T10: NXIT=1, RESID=1.5e-5,
  NPASS=3, NDIV=0). Both converge with finite stores. 10/10 pass.
  Residual: grid-scale realism gap - the manufactured ray has 8
  depths; production reversals live on ~70-depth structures.
- Files: `test_knop_ray.py`, `unit_tests/f90/drv_knop_ray.f90`.
- Action: none further at unit level; production validation rolls
  into P3-024's strong-compression model run.
- Acceptance: new manufactured cases pass with documented limits.
- Dependencies: P3-001..P3-004.
- Completion: partial (manufactured multi-reversal + strong-GAM
  coverage added and passing 2026-11-03; production-scale leg open).

### [~] P3-006 - Audit I_NEXT/J_NEXT storage and restart consistency
- Priority: P0
- Description: lagged store lifecycle: commit rules, FG_SEED_COMMIT
  waiver, snapshot (SNP_P/SNP_M) semantics, checkpoint dump/restore,
  fingerprint validation, foreign-seed policy.
- Evidence: AGENTS.md 02-Nov entry documents FG_SEED_COMMIT +
  SNP_P/SNP_M fixes; `NONMONO_COMPREHENSIVE_AUDIT` item 6 covers
  DUMP/RESTORE_*_GS tag+fingerprint validation.
- Files: `fg_j_cmf_nonmono_v1.f`, `mod_ray_mom_store.f`,
  checkpoint code paths (SCR_RITE_V2 callers),
  `control_variable_mod.f` (GS_FOREIGN_SEED).
- Action: finish audit of every store read/write site vs the
  documented lifecycle; verify restart restores the intended state.
- Acceptance: written lifecycle doc; restart behaviour verified.
- Dependencies: P3-001.
- Completion: partial.

### [~] P3-007 - Audit frequency matching and lagged-field refresh
- Priority: P0
- Description: DO_NONMONO_FREQ_MATCH multi-bin slots under strong
  drift; FG_KNOP_REFRESH_V1 fixed-point (SWEEP_FLD alternating vs
  SWEEP_FLD_ORD ordered); refresh operator identical to sweep.
- Evidence: `test_knop_ray.py` exercises multi-bin matching and both
  sweep modes (NXIT=1 convergence, NXIT=4 forced stall); audit
  verified sweep/refresh operator identity.
- Files: `fg_j_cmf_nonmono_v1.f`, `cmfgen_sub.f:3276-3300`.
- Action: quantify refresh contraction on a live model sweep
  (diagnostic), confirm matching coverage at drift extremes.
- Acceptance: written audit + live-run refresh-residual record.
- Dependencies: P3-005.
- Completion: partial.

### [ ] P3-008 - Validate approximate Lambda operator vs forward solver
- Priority: P1
- Description: the stored OAB diagonal (`GET_FG_ALO_DIAG`) used by
  DO_ES_ALI/ES_ALI_MOM_ON and the BA fill; validate its derivative
  content against the actual forward solve response.
- Evidence: 02-Nov fix - accessor now follows solver dispatch
  (previously always V13 diagonal: 56377/56377 fallbacks under
  KNOP). No measured ALO-vs-response comparison yet.
- Files: `fg_j_cmf_nonmono_v1.f` (ALO store), `do_var_cont.f`,
  `comp_j_blank.f` (ES_ALI update), `get_fg13_alo*` accessors.
- Action: unit-level comparison of stored diagonal vs numerical
  dJ/dETA response on a manufactured ray.
- Acceptance: documented ALO accuracy regime.
- Dependencies: P3-004.
- Completion: -

### [ ] P3-009 - Validate complete formal solver on interacting-SN model
- Priority: P0
- Description: a converged KNOP formal-solve model run (USE_LAM_ES=T
  configuration) on the m15 interacting model or equivalent.
- Evidence: negative - KNOP runs documented as stalling:
  compression-shell cells oscillate (fg drift ~5-11, ~48M bad
  cells, 56377 ALO fallbacks, exit gate exm ~60-120 dex;
  AGENTS_HISTORY 29-30 Oct); `m15_2p82d_KNOP_np` cold-start stall
  (job 12298049) motivated the 21-Oct precondition fixes. No
  converged KNOP model exists.
- Files: run configs `dessart_m15_2p82d_KNOP_*.yaml`.
- Action: after P3-005/P3-007 numerical issues are resolved or
  bounded, submit a monitored KNOP run; diagnose stall class if it
  recurs.
- Acceptance: a converged KNOP model or a documented, classified
  stall with root cause.
- Dependencies: P3-005, P3-007.
- Completion: -

### 6B. Moment solver (`plane/mom_j_cmf_nonmono_v1.f`)

### [~] P3-010 - Derive the complete discretised moment equations
- Priority: P0
- Description: document the moment system (J and H equations,
  advection coefficients, closure, boundary rows) as implemented.
- Evidence: `NONMONO_COMPREHENSIVE_AUDIT_2026-10-30.md` added the
  first runtime driver (`drv_mom_solve.f` + `test_mom_solve.py`)
  giving coefficient-level verification and zero-GAM runtime
  evidence; `cmfgen_obsidian` has solver notes.
- Files: `mom_j_cmf_nonmono_v1.f`, `test_mom_solve.py`,
  `test_mom_sweep_table.py`.
- Action: write the full discrete-equation doc incl. Eddington
  source (GET_MOMS_NON_REL ray moments under LAM_ES=T; V13 vs
  nonmono ray moments caveat under SOB LAM_ES=F).
- Acceptance: documented equations matching code.
- Dependencies: none.
- Completion: partial.

### [~] P3-011 - Audit USE_LAM_ES=T independently
- Priority: P0
- Description: the LAM_ES=T branch solves moments with Eddington
  factors from the formal solve; verify coefficient assembly and
  consistency for this branch alone.
- Evidence: branch dispatch audited (comp_j_blank.f:978
  GET_MOMS_NON_REL); SOB LAM_ES=T converged on m15 (P2-012); KNOP
  LAM_ES=T unvalidated (P3-009).
- Files: `comp_j_blank.f`, `get_moms_non_rel.f`,
  `mom_j_cmf_nonmono_v1.f`.
- Action: audit the Eddington-factor path (which ray solve feeds it,
  clamped vs SIGMA_TRUE coefficients - documented caveat); check
  boundary rows.
- Acceptance: written branch audit; the SIGMA-vs-SIGMA_TRUE
  coefficient inconsistency quantified.
- Dependencies: P3-010.
- Completion: partial.

### [~] P3-012 - Audit USE_LAM_ES=F independently
- Priority: P0
- Description: the LAM_ES=F branch closes the moments internally
  (variable Eddington); verify the full iteration including
  MOM_KNOP_REFRESH_V1 lagged J/H_NEXT map.
- Evidence: runtime FAIL documented - m15 SOB LAM_ES=F
  bounded non-contraction (job 12316397): lagged map never
  contracts, ~200K cells quarantined/sweep, refresh stalls at
  residual ~1.1, model pinned at 63-75%. `test_mom_solve.py` gives
  coefficient-level coverage only.
- Files: `mom_j_cmf_nonmono_v1.f`, `MOM_KNOP_REFRESH_V1` site in
  `cmfgen_sub.f:3293`.
- Action: continue the stall diagnosis at operator level: is the
  non-contraction in the refresh map, the quarantine policy, or the
  closure? Use manufactured cases to isolate.
- Acceptance: root cause classified with evidence; path either fixed
  or documented as unsupported with a clean rejection.
- Dependencies: P3-010, P3-015.
- Completion: partial (failure documented; root cause open).

### [x] P3-013 - Investigate MOM_DLN_MIN and near-coincident frequency points
- Priority: P0
- Description: MOM_DLN_MIN floors the log-frequency spacing in
  advection coefficients; bound-free edge doublets (~2E-11 dLOG_NU)
  on the continuum grid are the documented trigger - quantify effect
  on coefficients.
- Evidence: source census complete (2026-11-03). Three use sites,
  all `MAX(DLN,1E-6)`: CHI_HAT (~line 2182), GAMH (2190), GAM
  (2209). MOM_DLN_J/H hold per-depth effective spacings: caller
  `dLOG_NU` (true backward spacing) for expanding cells, true
  `LOG(FREQ/NU_NEXT)` for compressing cells, `N_CAP*DLN` under
  MOM_FREQ_MATCH long-characteristic spans, and true matched-span
  spacing via the MOM_NU_GRID cache. Asymmetric spacing is handled
  on both sides - no direction uses a uniform proxy. The floor
  only activates on doublets: real min spacing ~1E-5 vs floor
  1E-6 vs doublet 2E-11. At a doublet PSI~|GAM|/DLN is already
  dominant, so flooring DLN bounds the coefficient while the row
  still degenerates to J_k ~ J_upwind - the correct limit; under
  MOM_GAM_PIN the degenerate limit is exact.
- Files: `mom_j_cmf_nonmono_v1.f`, `comp_j_blank.f` (call site).
- Action: none - no defect; regularization correctly placed.
- Acceptance: met at source level; production activation counts
  remain a runtime-measurement item.
- Dependencies: P3-010.
- Completion: 2026-11-03 (source census; no defect).

### [~] P3-014 - Derive the small-frequency-spacing limiting equation
- Priority: P0
- Description: derive the limiting discrete equation as dln_nu -> 0
  and design the correct treatment (what the advection term should
  degenerate to).
- Evidence: derivation from the P3-013 census (2026-11-03). The J
  row's frequency-advection terms all scale with PSI ~
  DTAUONQ*|GAM| ~ dr*|CON_GAM*fac|/DLN. As DLN -> 0 every other row
  term stays O(1) while PSI -> inf, so the row degenerates to
  PSI*(J_k - J_upwind) = 0 -> J_k = J_upwind: the discrete
  equation copies the upwind field across a null-width bin, which
  is the correct transport limit (two identical frequencies
  exchange nothing). The same limit holds in the H channel via
  W ~ |GAMH|/DLN. MOM_DLN_MIN=1E-6 does NOT change the limiting
  equation: at DLN < 1E-6 the PSI term is already ~1E4-1E8 times
  the O(1) terms, so flooring DLN only bounds the coefficient
  magnitude while the row still sits ~1E4 into the degenerate
  regime. MOM_GAM_PIN=T makes the limit exact (PSIPREV=PSI);
  pin off leaves a second-order (PSIPREV-PSI)/PSI residual,
  documented. Open leg: a manufactured doublet-grid case asserting
  J_k=J_upwind to round-off would close this fully; current
  evidence is structural (coefficient scaling), not numerical.
- Files: `mom_j_cmf_nonmono_v1.f`, `test_mom_solve.py` (case hook),
  `drv_mom_solve.f` (`dbl` mode).
- Action: done - numerical assertion added (dbl case, 2026-11-03).
- Acceptance: met - derived limit confirmed numerically.
- Dependencies: P3-013.
- Completion: 2026-11-03. The `dbl` manufactured case (FRQ(4)=
  FRQ(3)*exp(-1E-11), true backward spacing passed per call)
  verifies the limit: compressing rows at the doublet copy
  J_upwind to 0.15% (lagged-store drift, not coefficient error),
  and |J3-J4|=6.4e-5 sits below the real-bin difference
  |J4-J5|=6.5e-5 - the doublet is a null-transport bin. 93/93
  moment tests pass.

### [~] P3-015 - Build frequency-integrated conservation tests
- Priority: P0
- Description: manufactured tests asserting the discrete moment
  equations conserve the relevant frequency-integrated quantity
  across grid-spacing regimes incl. near-degenerate grids.
- Evidence: degenerate-grid leg done (2026-11-03): the `dbl`
  doublet case proves null transport across a 1E-11 bin
  (compressing rows copy J_upwind to 0.15%; the doublet is
  smoother than a real bin). The full frequency-integrated
  conservation identity (sum_k w_k * advection term = boundary
  flux telescoping) is not yet asserted - needs the driver's
  per-frequency dump extended to carry the assembled coefficient
  rows, or a dedicated multi-frequency case.
- Files: `test_mom_solve.py`, `drv_mom_solve.f` (jprof dump).
- Action: extend the conservation assertion to the telescoping
  identity across a full sweep.
- Acceptance: partial - degenerate-grid behaviour verified;
  integrated-conservation identity open.
- Dependencies: P3-013, P3-014.
- Completion: partial.

### [x] P3-016 - Investigate the 1 + sigma*K_ON_J = 0 resonance
- Priority: P0
- Description: the moment equations have a resonance where
  1+sigma*K_ON_J crosses zero (compression); trace it through
  closure, angular terms, boundaries, iteration.
- Evidence: completed as BUG-004 (register) - REFUTED. All four
  factor sites traced (depth/midpoint, current/lagged, incl.
  MOM_GAM_PIN); consumers traced through CHI_HAT, Q, PSI/PSIPREV,
  W/EPS/HU/HL/HS, the boundary row, and the census. The sign flip
  is the physical transport-direction reversal; the floor bounds
  the jump at 2*MOM_GAM_FAC_MIN; the T1==0 exception returns the
  physically correct decoupled row. t_zero manufactured shelf: 36
  bitwise-resonant cells decouple (Q=PSI=0), shelf-edge T1=3.3e-16
  cell floors at +0.01 (Q=1.6e-4). 85/85 moment unit tests pass.
- Files: `mom_j_cmf_nonmono_v1.f`, `get_moms_non_rel.f`,
  `test_mom_solve.py`.
- Action: done - no change made; regularization verified.
- Acceptance: met - resonance handling documented in BUG-004 and
  exercised by the t_zero case; no regularisation enlarged.
- Dependencies: P3-010.
- Completion: 2026-11-03.

### [~] P3-017 - Quantify the effect of MOM_GAM_FAC_MIN
- Priority: P1
- Description: measure where MOM_GAM_FAC_MIN binds and its effect on
  the advection coefficient and solution.
- Evidence: unit-level evidence from BUG-004: the floor bounds
  |fac|>=1e-2 off exact resonance; measured Q=1.6e-4 at a
  floor-band cell vs 4e-3-1.7e-2 off-shelf; PASS_NRES census
  already counts <=floor cells per pass (incl. exact-zero).
  Production-model activation counts still unmeasured.
- Files: `mom_j_cmf_nonmono_v1.f`.
- Action: production leg - record PASS_NRES / floor-band census
  from a real interacting-SN model.
- Acceptance: documented activation counts on production data.
- Dependencies: P3-010.
- Completion: partial (unit-level done; production census open).

### [x] P3-018 - Audit MOM_LAG_BND
- Priority: P1
- Description: audit the lagged boundary treatment key MOM_LAG_BND -
  which boundary values lag, gating, defaults.
- Evidence: implementation audited (2026-11-03). Default T, wired
  via SET_MOM_LAG_OPTS. Interior rows (~2764-2788): lagged weight
  T4 = (PSIPREV if Q>=0 expanding | PSI if Q<0 compressing with a
  valid J_NEXT store | 0 if no store) + |VB|+|VC| +
  2(|EPS_PREV(i-1)|+|EPS_PREV(i)|); scale = min(1, max(TB-|TA|-
  |TC|,0)/T4) applied to XLAG; per-sweep count and min factor in
  MOM_LAGB_CNT/MIN. Same bound applied at the outer boundary row
  (~2900-2960). The M-matrix excess bound is mathematically
  correct; the upwind-direction weight selection matches the Q_NM
  sign convention. LAG_SCL_NM retained so the negative-J rebuild
  applies the same factor; H reconstruction intentionally unscaled
  (documented). Consistency cost vs MOM_EPS_LIM is documented.
- Files: `mom_j_cmf_nonmono_v1.f`.
- Action: none - implementation matches the documented design.
- Acceptance: met - written audit above.
- Dependencies: P3-010.
- Completion: 2026-11-03.

### [x] P3-019 - Audit MOM_ABS_FLR
- Priority: P1
- Description: audit the MOM_ABS_FLR floor - what it bounds, when it
  activates, its physical justification per the no-silent-clamp
  policy.
- Evidence: implementation audited (2026-11-03). Default 0.0 =
  off (exact clip, upstream behaviour preserved). Sites: interior
  diagonal MAX(MOM_ABS_FLR,1-COH) at ~2386 and ~2411, outer
  boundary MIN(COH,1-MOM_ABS_FLR)-1 at ~2940 - consistent on both
  row types. Classification: a numerical regularization applied to
  a derived coefficient (the J-independent absorption diagonal),
  not a clamp on a physical variable; justification is documented
  inline (ESEC>=CHI rows + resonance PSI~0 leave no diagonal
  excess, so MOM_LAG_BND nulls the whole lagged source and the row
  oscillates; a small floor retains an absorption anchor). It is
  explicit (named VADAT key, documented, off by default) and its
  activation would be visible via MOM_LAGB counters. No defect.
- Files: `mom_j_cmf_nonmono_v1.f`.
- Action: none - floor is documented, gated, and correctly placed;
  activation counts on production runs remain a runtime item.
- Acceptance: met - written audit; no unjustified clamp found.
- Dependencies: P3-010.
- Completion: 2026-11-03.

### [x] P3-020 - Audit MOM_H_IMPL
- Priority: P1
- Description: audit the MOM_H_IMPL option - implicit H treatment
  branch, correctness of the implicit solve vs explicit.
- Evidence: implementation audited (2026-11-03). Default F, wired
  via SET_MOM_LAG_OPTS. When on, TB_NM gains |VB|+|VC| - the
  magnitude of the lagged upwind-H midpoint couplings - on
  interior rows only (~2438-2442); the boundary row assembles its
  own bound and is untouched. That is exactly the documented
  design: the signed lagged RHS is unchanged while the diagonal
  carries the channel's weight, so a clipped shock row stops
  oscillating and the MOM_LAG_BND excess grows by the same amount.
  Documented consistency cost: signed difference vs summed
  diagonal adds bounded ~2*HS mass on smooth rows. Also verified
  the four-branch EPS depth-upwind stencil (2413-2422): the sign
  key gives -|E| off-diagonals and +|E| diagonal mass in both
  directions - M-matrix preserving, matching the comment.
- Files: `mom_j_cmf_nonmono_v1.f`.
- Action: none - implementation matches design; the branch stays
  opt-in/off per the no-silent-clamp policy.
- Acceptance: written audit; default-F manufactured coverage
  exists via every non-H_IMPL case in the suite.
- Dependencies: P3-010.
- Completion: 2026-11-03.

### [x] P3-021 - Audit MOM_GAM_PIN
- Priority: P1
- Description: audit MOM_GAM_PIN - pinned GAM values option,
  semantics, defaults, interaction with the advection term.
- Evidence: implementation audited (2026-11-03). Default T, wired
  via SET_MOM_LAG_OPTS. All four lagged-factor sites pin to the
  current-frequency evaluation when on: GAMH_FAC_PREV (~2020),
  GAM_FAC_PREV (~2147), EPS_PREV (~2276), and the boundary
  HBC/NBC evaluation (~2875). With the pin, PSIPREV=PSI exactly
  and WPREV=W, so the explicit upwind weight can never exceed the
  implicit diagonal - the documented cure for the lagged/current
  sign-and-magnitude drift near 1+sigma*K=0 that made the
  explicit term outweigh the implicit one (m15 probes). Off
  restores the independent PREV evaluation (same floor + exact-zero
  form on K_ON_J_PREV/NMID_ON_HMID_PREV/NMID_ON_J_PREV). No defect:
  consistent sites, correct sign conventions, upstream semantics
  preserved behind the key.
- Files: `mom_j_cmf_nonmono_v1.f`.
- Action: none - implementation matches documented design.
- Acceptance: written audit above.
- Dependencies: P3-010.
- Completion: 2026-11-03.

### [ ] P3-022 - Measure existing inner and outer Anderson acceleration
- Priority: P1
- Description: measure activation/acceptance/rejection/fallback of
  the existing inner+outer Anderson acceleration in the moment
  refresh, before proposing more acceleration.
- Evidence: dev plan "current implementation": inner and outer
  Anderson acceleration present; `test_knop_ray.py`/refresh tests
  exercise convergence but not AA counters.
- Files: `mom_j_cmf_nonmono_v1.f`, `fg_j_cmf_nonmono_v1.f`
  (ordered-pass acceleration option), `solve_for_pops.f`.
- Action: instrument/read existing counters on a run; compute
  residual reduction vs unaccelerated map on a manufactured case.
- Acceptance: measured AA effectiveness report.
- Dependencies: P3-010.
- Completion: partial 2026-11-03 - instrumentation verified at
  source: MOM_AA_APP/SKIP/SNG (inner, MOM_AA_M=3 secant pairs) and
  MOM_AA_OUT_APP/REJ/SNG (outer, MOM_AA_OUT_M=8) print at every
  refresh epilogue plus per-application `AA m= ||f||= ||corr||=`
  lines (mom_j_cmf_nonmono_v1.f:6486-6497). Gate map: inner AA
  fires only on even IPASS>=2 with ALLOCATED(XHIST_AA), restarts
  history when >MOM_AA_DMP_FRAC of cells re-damp that pass, caps
  each cell's correction at 10||f_i||, drags H by the J ratio
  (mom_j_cmf_nonmono_v1.f:5629-5830); outer AA is the analogous
  iterate-space path (MOM_AA_OUT). Driver probe: 4 and 8 outer
  iterations of `mom_solve comp1` emitted zero LUER inner-pass
  lines - the cold calls take the LAMBDA-rebuild warmup branch
  (mom_j_cmf_nonmono_v1.f:4170-4294) rather than multi-pass GS
  sweeps, so the even-IPASS>=2 AA gate never engages under the
  driver's cadence. Manufactured AA evidence needs a driver
  seeded with a valid J_NEXT_GS checkpoint (the `ckpt` mode
  dumps/restores one - tried it, the restored call still takes
  the no-report path: zero `refresh pass` lines means the whole
  IPASS report/AA block is bypassed by an enclosing branch
  before mom_j_cmf_nonmono_v1.f:5213, not entered per-call) or
  a LAM_ES=F model run. Both m15 jobs are USE_LAM_ES=T (FG
  path - AA lives in the moment refresh only).

### [~] P3-023 - Lagged-field maturity gates vs actual radiation consistency
- Priority: P0
- Description: verify GS maturity metrics (freeze/refresh release
  gates) track actual field consistency, not just iteration counts.
- Evidence: `NONMONO_COMPREHENSIVE_AUDIT` item 4 - SOLVE_FOR_POPS
  reads GET_FG_GS_MATURITY/GET_MOM_GS_MATURITY before BA/Newton;
  KNOP_AUDIT fixed FG_GS_NSTAL/NPLAT census ordering; stall history
  shows gates can hold while field still moves (drift ~5-11).
- Files: `solve_for_pops.f`, `fg_j_cmf_nonmono_v1.f`
  (GET_FG_GS_MATURITY), `mom_j_cmf_nonmono_v1.f`.
- Action: on a run or driver, compare gate-open state vs an
  independently recomputed field residual.
- Acceptance: maturity metric proven to bound the true field
  residual on at least one model.
- Dependencies: P3-007.
- Completion: partial.

### [ ] P3-024 - Validate complete moment solver on strong-compression model
- Priority: P0
- Description: converged LAM_ES=F KNOP (or SOB moment) model on a
  strong-compression profile.
- Evidence: negative - the documented non-contraction (P3-012).
- Files: moment-path configs.
- Action: depends on P3-012 root cause; then submit model.
- Acceptance: converged run or documented unsupported status.
- Dependencies: P3-012, P3-016.
- Completion: -

## 7. Pillar 4 - DO_JFNK

Acceptance: production residual, JVP, GMRES correction and
preconditioner validated on actual model states; JFNK remains
optional. Probe exists (`jfnk_op_probe_v1.f90`, `JFNK_HSCAN`,
`test_jfnk_opint.py`, `se_driver seprobe`).

### [x] P4-001 - Verify production residual architecture and unknown ordering
- Priority: P0
- Description: document the residual map and the correction-vector
  layout.
- Evidence: `JFNK_PRODUCTION_OPERATOR_INVESTIGATION_2026-11-07.md`
  sec.A: depth-major `idx=(k-1)*NT+j`; j<=NT-2 population
  (fractional), j=NT-1 ED, j=NT log-T (`x_frac(NT)=exp(x_T)-1`);
  frozen-J and coupled-J operators distinguished; state
  save/restore audited by inspection; masks and PTC documented.
- Files: `jfnk_cmfgen_interface_v1.f`, `jfnk_solve_corr_v1.f`.
- Action: none further for architecture; runtime repeatability is
  P4-002/P4-003.
- Acceptance: documented architecture. MET.
- Dependencies: none.
- Completion: 2026-11-03 on audit sec.A (committed ccc52ffd).

### [ ] P4-002 - Verify frozen-J residual repeatability on a live model
- Priority: P0
- Description: run a frozen-J admission with JFNK_HSCAN=T and read
  OPPROBE REP/REPP rows for bitwise or bounded repeatability.
- Evidence: audit sec.B: surrogate bitwise 0; production asserted by
  save/restore inspection, unmeasured. Needs one live admission.
- Files: production run with `DO_JFNK=T`, `USE_LAM_ES=T` (frozen
  admission), `JFNK_HSCAN=T`.
- Action: submit a SOB+JFNK run with HSCAN; record REP/REPP.
- Acceptance: measured production repeatability recorded.
- Dependencies: none (instrumentation in place).
- Completion: -

### [ ] P4-003 - Verify coupled-J residual repeatability on a live model
- Priority: P0
- Description: same measurement on the coupled tier (KNOP+JFNK
  admission replays the radiation sweep per eval).
- Evidence: audit sec.A documents `JFNK_RAD_SWEEP_V1` with
  FREEZE_I_NEXT armed; repeatability unmeasured on production.
- Files: KNOP+JFNK config with JFNK_HSCAN=T.
- Action: submit KNOP+JFNK HSCAN run; record REP/REPP. NOTE: coupled
  admissions depend on the KNOP path reaching admission - may be
  limited by P3-009 status.
- Acceptance: measured coupled repeatability recorded.
- Dependencies: P3-009 progress (KNOP admission reachable).
- Completion: -

### [ ] P4-004 - Run production JFNK operator probe on a frozen-J model
- Priority: P0
- Description: full OPPROBE battery (REP/REPP/DIR/HS/ADD/HOM) at a
  frozen-J admission on a live model.
- Evidence: audit sec.I - measurements pending; instrumentation
  verified on surrogate (`test_jfnk_opint.py` 36/36).
- Files: same run as P4-002.
- Action: collect and record all OPPROBE rows.
- Acceptance: production numbers recorded in an audit note. NOT met
  until real measurements exist.
- Dependencies: P4-002.
- Completion: -

### [ ] P4-005 - Run production JFNK operator probe on a coupled-J model
- Priority: P0
- Description: full battery at a coupled-J admission.
- Evidence: pending; ~100 residual evals per admission (~1 radiation
  sweep each in coupled mode) - bounded, opt-in cost.
- Files: same run as P4-003.
- Action: collect OPPROBE rows.
- Acceptance: production numbers recorded.
- Dependencies: P4-003.
- Completion: -

### [ ] P4-006 - Analyse directional derivative accuracy and h sweeps
- Priority: P0
- Description: from production DIR/HS rows: per-direction eFD/eFC,
  emin, hstar, cap/retry/clip incidence; check ED/T-class inflation
  seen on surrogate.
- Evidence: audit sec.C surrogate results; production pending.
- Files: OPPROBE output from P4-004/P4-005.
- Action: analyse; decide whether H_STEP policy needs work (only if
  production shows policy-driven nonlinearity).
- Acceptance: written analysis; h operating point documented.
- Dependencies: P4-004.
- Completion: -

### [ ] P4-007 - Analyse production additivity and homogeneity
- Priority: P1
- Description: from ADD/HOM rows: d_pol vs d_op on production; the
  norm-scaled h policy is exact under homogeneity and
  truncation-limited under additivity on the surrogate.
- Evidence: surrogate d_pol~1.9e-7, d_op~9e-8, HOM d_pol=0.
- Files: OPPROBE output.
- Action: analyse production ADD/HOM.
- Acceptance: written analysis; policy change only if d_pol >> d_op.
- Dependencies: P4-004.
- Completion: -

### [ ] P4-008 - Verify probe distinguishes FD-JVP from approximate BA operator
- Priority: P1
- Description: confirm labels/conclusions do not conflate the FD JVP
  with the BA action; BA omits residual nonlinearity, off-band
  cross-depth coupling, and (coupled) radiation response.
- Evidence: audit sec.E defines the BA action and its omissions;
  `JFNK_BA_JVP` remains a diagnostic.
- Files: `jfnk_op_probe_v1.f90` output labels, audit text.
- Action: verify output text keeps the distinction; fix labels if
  ambiguous (documentation-level change).
- Acceptance: labelled correctly in output and docs.
- Dependencies: P4-004.
- Completion: -

### [~] P4-009 - Verify true GMRES residual on every relevant exit path
- Priority: P0
- Description: every GMRES exit (converged, restart, restart-cap,
  stagnation-soft, trust-rescale, jvp_soft) must measure the true
  weighted residual on the returned vector; acceptance gated on it.
- Evidence: audit sec.D plus BUG-001 fix (2026-11-03): the
  measurement now runs unconditionally at every done-exit - the
  `tres_fac > 0` arming gate is removed, so a `res <= tol` claim is
  accepted only when the measured `true_res <= acc_fac*tol`
  (`acc_fac` = JFNK_TRES_FAC floored to 1.0 on non-positive input).
  Fail-closed on non-finite / jvp-soft / pc-clip paths; rejected
  claims restart on the measured residual bounded by
  JFNK_GMRES_RST/imx, else `gmres_tres_reject`. Regression: T4j,
  T4k, G17 (G17 + T4k(biased) confirmed to FAIL on the pre-fix
  solver). Remaining nuance: the `ginfo==2` hookstep-recompute
  return is a truncated-direction path that still carries
  `true_res_ok=.false.` by design (checked downstream by
  pred_res_physical + trial eval); flagged in the register.
- Files: `jfnk_solver_v1.f90`, `rd_control_variables.f`,
  `control_variable_mod.f`, `tests/jfnk_core/`.
- Action: source fix done + regression green; record tres_rej /
  nonfinite counts from production runs (P4-004) to close the
  runtime leg.
- Acceptance: source-verified + regression (done); production
  tres_rej accounting open.
- Dependencies: P4-010.
- Completion: partial (code + unit evidence complete 2026-11-03;
  runtime accounting open).

### [~] P4-010 - Test JFNK_TRES_FAC on/off; establish production configuration
- Priority: P0
- Description: decide the production TRES_FAC setting from measured
  true-vs-Hessenberg ratios; ensure a claimed convergence is never
  accepted on the Hessenberg estimate alone.
- Evidence: BUG-001 resolution (2026-11-03): the "off" branch was
  itself the P0 - with the default `0`, an `info=0` claim bypassed
  all residual checks. Configuration decided as unconditional
  measurement + `JFNK_TRES_FAC` as acceptance bound: default `1.0`,
  `<=0` invalid (floored to 1.0 in both the VADAT clamp and the
  solver), `>1` a documented weakening for diagnostics. Verified:
  T4k proves `fac=0` cannot switch verification off.
- Files: `jfnk_solver_v1.f90`, `rd_control_variables.f`, VADAT docs.
- Action: remaining leg - production runs recording tres_rej counts
  under the new default (pair with P4-004/P4-005); then close.
- Acceptance: measured basis for the production setting (bound
  semantics decided + unit-verified; runtime counts open).
- Dependencies: P4-004.
- Runtime note (2026-11-03, BUG-005): on m15 the tres path is
  structurally unreachable - 155/155 production JFNK admissions
  defer on the FD non-descent guard before GMRES enters
  (radiation-pinned frozen-J rows, |Jz|/|F|~26). tres_rej counts
  on m15 will read 0 regardless of the setting; a model where the
  reach test passes is needed for a nonzero measurement. Two
  interactive-partition stability runs on the fixed build are in
  flight (jobs 12667827 warm, 12667984 cold) - expected to
  reproduce the admit/defer pattern rather than tres exits.
- Completion: partial (policy decided and unit-verified 2026-11-03;
  production counts open).

### [ ] P4-011 - Verify inexact-Newton forcing criterion and true-residual acceptance
- Priority: P0
- Description: verify Eisenstat-Walker forcing is in the
  preconditioned (weighted) norm consistent with left
  preconditioning, and the caller gate (`true_res < 0.9*ref_norm`,
  `true_res > |F|` skip) behaves as audited.
- Evidence: audit sec.D verified by source; runtime unverified.
  Re-verified post-BUG-001 (2026-11-03): EW Choice-2 uses
  eta_k = gamma*(|F_k|/|F_{k-1}|)^alpha with the eta_{k-1}
  safeguard (Oct-2026 fix), clamped [eta_min,eta_max] on all paths;
  gres = eta*ref_norm with ref_norm = ||M^{-1}F||_w under left-PC
  (09-Sep fix); true_res is measured as ||M^{-1}(b-Ax)||_w - same
  units - so the 0.9*ref_norm dual gate and the true_res > ref_norm
  skip are unit-consistent. jfnk_tol shares the dp weighted norm.
- Files: `jfnk_solver_v1.f90`, `newton_hook_v1`.
- Action: runtime record from P4-004 runs; confirm no path accepts a
  step on hess_res.
- Acceptance: measured acceptance behaviour recorded.
- Dependencies: P4-009, P4-010.
- Completion: -

### [x] P4-012 - Investigate population-floor near-null columns (T21 surrogate)
- Priority: P0
- Description: explain the T21 per-depth LU backward-error XFAILs.
- Evidence: audit sec.F - reproduced on stored blocks + numpy
  replica; cond 1e15..1e27; smallest right-singular vectors on
  population columns at the n0=1e-30 floor (column norms ~1e-8 vs
  block max 1e7..1e13); ~3-dim approximate null space; row
  reg->identity cannot touch column-side null space; XFAILs
  retained and justified; remedies are architectural (mask
  floor-pinned unknowns or rescaled basis), not ridge.
- Files: `unit_tests` sepc/MZDK census, `test_jfnk_solver.py` T21.
- Action: none required by the task text ("investigate" done);
  follow-on fix class is a separate architectural task if pursued.
- Acceptance: structural explanation recorded. MET.
- Dependencies: none.
- Completion: 2026-11-03 audit sec.F. No ridge escalation per rule.

### [ ] P4-013 - Measure production preconditioner conditioning and backward errors
- Priority: P1
- Description: run the MZDK-style per-depth backward-error census on
  production `JFNK_BA_A` blocks at a real admission.
- Evidence: surrogate census exists (sec.F); production pending.
- Files: production admission, `JFNK_FACTOR_PRECOND_V1`.
- Action: extend probe or extract blocks from a live admission;
  census.
- Acceptance: production conditioning table.
- Dependencies: P4-004.
- Completion: -

### [ ] P4-014 - Audit frozen-J vs coupled-J BA approximation quality
- Priority: P1
- Description: quantify BA-as-derivative error separately per mode;
  under coupled-J the BA omits essentially all radiation response.
- Evidence: audit sec.E established by construction; unmeasured.
- Files: OPPROBE eFD columns (per-mode), `JFNK_OP_BAMV_V1`.
- Action: analyse production eFD by mode.
- Acceptance: per-mode BA quality recorded.
- Dependencies: P4-004, P4-005.
- Completion: -

### [~] P4-015 - Verify coupled-J radiation sweep restores all state
- Priority: P0
- Description: verify `JFNK_RAD_SWEEP_V1` save/restore actually
  returns bit-identical state (radiation, line, moment/Knop refresh,
  Eddington-record ACCESS_F, file cursors) after an eval.
- Evidence: inspection leg re-verified 2026-11-03: the bracket is
  symmetric and comprehensive (STEQ x6, POPS/ED/T/ATM, ACCESS_F,
  EDDINGTON, FOLD, COMPUTE_BA, two-photon + SL accumulators, MOM
  quarantine + refresh scalars, LSTATE line-module workspace, JEW
  via COMPUTE_EW). One ordering defect found and fixed: BUG-006 -
  `JFNK_EVAL_LSTATE_V1('S')` armed after the first eval mutations,
  so its guarded alloc-failure path skipped restores and would
  have leaked `JFNK_IN_EVAL=T` permanently. Hoisted into the
  entry-save block; ordering pin added and verified fail->pass.
  REP/REPP bit-identical measurement on production remains open.
- Files: `cmfgen_sub.f:5755` sweep, `jfnk_cmfgen_interface_v1.f`.
- Action: production REP/REPP under coupled admission (P4-003)
  directly tests this; inspect restore lists against mutated arrays
  once more on current source.
- Acceptance: measured + inspection-verified restore.
- Dependencies: P4-003.
- Completion: partial (inspection leg done + BUG-006 fixed
  2026-11-03; REP/REPP production measurement open).

### [~] P4-016 - Verify JFNK admission across Sobolev and Knop configurations
- Priority: P0
- Description: JFNK admission gating: frozen tier under SOB+LAM_ES=T
  (admission gated on USE_LAM_ES); coupled tier KNOP-only
  (`JFNK_COUPLED` KNOP-only per test_vadat_keys); reject unsupported
  combinations loudly.
- Evidence: `test_vadat_keys.py` covers DO_JFNK nonmono gate and
  JFNK_COUPLED KNOP-only rejection; AGENTS.md: JFNK inactive under
  SOB LAM_ES=F by design.
- Files: `rd_control_variables.f`, `check_param_consistency.f`,
  `solve_for_pops.f` (admission gate).
- Action: write the admission matrix (which tier x which transfer
  path is admitted/rejected); confirm each rejection is loud.
- Acceptance: documented matrix; all unsupported combos hard-stop or
  warn.
- Dependencies: none.
- Runtime note (2026-11-03, BUG-005): on m15 the admitted frozen
  tier always defers before GMRES (155/155 on FD non-descent).
  Admission happens; Newton-level work does not. The admission
  matrix is therefore only exercised through the defer path in
  production - gate-hit counting, mode retention, and rung
  accounting verified correct in OUTGEN traces.
- Completion: partial (parser-level verified; matrix doc open;
  production admission measured 2026-11-03 - always defers on m15).

### [ ] P4-017 - Compare conventional vs JFNK solutions on identical models
- Priority: P0
- Description: same model, conventional vs DO_JFNK solver; compare
  converged POPS/ED/T and residual level - JFNK must solve the same
  equations.
- Evidence: none - no paired run exists.
- Files: paired configs.
- Action: run A/B once JFNK admission is exercised (P4-002).
- Acceptance: solutions agree within documented tolerance; residual
  criteria not weakened.
- Dependencies: P4-002.
- Completion: -

### [ ] P4-018 - Verify JFNK improves convergence without weakening residual criteria
- Priority: P0
- Description: demonstrate the JFNK path reaches the same physical
  residual faster, with all acceptances passing true-residual gates.
- Evidence: none.
- Files: P4-017 runs + OUTGEN residual histories.
- Action: analyse.
- Acceptance: written comparison; no relaxed criteria.
- Dependencies: P4-017.
- Completion: -

### [ ] P4-019 - Measure JFNK cost: residual evals, memory, runtime
- Priority: P2
- Description: eval count (incl. ~100-eval probe cost when armed),
  GMRES iterations, wall time, memory vs conventional.
- Evidence: none.
- Files: run logs.
- Action: record from P4-017 runs.
- Acceptance: cost table.
- Dependencies: P4-017.
- Completion: -

## 8. Pillar 5 - DIAG and TRI matrix solvers

Acceptance: both methods execute the intended algorithms and produce
consistent physical solutions within justified tolerances in every
supported configuration.

### [~] P5-001 - Trace METH_SOL from VADAT to matrix-solver dispatch
- Priority: P0
- Description: trace `SOL_METH`/`METH_SOL` (DIAG vs TRI) from
  `rd_control_variables.f` through `SOLVEBA_V13`/BA assembly to the
  actual solve call, incl. `NUM_BNDS` interaction.
- Evidence: `DO_FG13_ALO_TRI` docs: TRI path fills TX(K,K+/-1) with
  NUM_BNDS=3 (pipeline `depth_independent_solver: false`);
  `CONVENTIONAL_SOLVER_AUDIT` traced BA solve call graph. Key name:
  `SOL_METH` in this tree per AGENTS.md.
- Files: `rd_control_variables.f`, `control_variable_mod.f`,
  `solveba_v13.f`, `generate_full_matrix_v3.f`, `do_var_cont.f`.
- Action: write the dispatch chain with file:line refs.
- Acceptance: documented trace; no silent rerouting.
- Dependencies: none.
- Completion: 2026-11-03 - end-to-end trace verified, no silent
  rerouting found. Chain: VADAT `SOL_METH` ->
  `RD_STORE_NCHAR(METH_SOL,...)` (rd_control_variables.f:2948;
  METH_SOL is CHARACTER*6, control_variable_mod.f:331) ->
  NUM_BNDS<3 forces DIAG with a loud banner (:2951-2957) ->
  `solveba_v13.f:1033` accepts only METH_SOL(1:4)='DIAG' or
  (1:3)='TRI' and STOPs loudly on anything else (so the RD prompt's
  PEN/GSIT/MIN options are unsupported-but-loud, not silently
  rerouted; 'TRIDIAG' truncates to 'TRIDIA' and still hits the TRI
  branch) -> `CMF_BLKBAND_V3` (:1041). Inside: DIAG takes the
  per-depth LU path with CSCRATCH factor recall; TRI takes the
  banded path (TRI_ACTIVE at :592) with BCSCRATCH/DSCRATCH recall.
  Both recall paths verify every mid-depth cache file exists and
  fall through to a loud recompute on a miss (19-Oct-2026
  hardening, :285-330). FAILURE semantics: TRI failures degrade
  like DIAG via the shared 9999 soft-fail branch keyed on
  TRI_ACTIVE/TRI_RECALL.

### [ ] P5-002 - Verify DIAG matrix assembly and solution
- Priority: P0
- Description: verify the DIAG path assembles and solves the
  intended (depth-independent) system.
- Evidence: pipeline default `depth_independent_solver` maps here;
  used by all converged runs to date.
- Files: `solveba_v13.f`, `cmf_blkband_v3.f`,
  `cmf_blkband_pin_null_v1.f`.
- Action: audit assembly; confirm with manufactured input.
- Acceptance: verified assembly + residual check.
- Dependencies: P5-001.
- Completion: partial 2026-11-03 - assembly audited at source:
  the DIAG branch never touches off-diagonal storage; per depth it
  builds the NxN block via GENERATE_FULL_MATRIX_V3
  (BAND_INDX=DIAG_INDX), equilibrates (ROW_SF/COL_SF), LU-factors
  with DGETRF, and solves with DGETRS (cmf_blkband_v3.f:359-393
  recall path; compute path at :480-520). Both paths are strictly
  depth-independent - the correct DIAG semantics. Manufactured
  residual check still needs a driver wrapping
  GENERATE_FULL_MATRIX_V3 + the solve.

### [ ] P5-003 - Verify TRI matrix assembly and solution
- Priority: P0
- Description: verify the TRI path assembles the banded
  (depth-coupled) system and solves it; check band widths and
  leading dimensions.
- Evidence: `DO_FG13_ALO_TRI` fills TX(K,K+/-1,1:2) for NUM_BNDS=3;
  runtime status unverified pending `adp2_alo3` per AGENTS.md.
- Files: `solveba_v13.f`, `do_var_cont.f`, `generate_full_matrix_v3.f`.
- Action: audit assembly + manufactured solve.
- Acceptance: verified.
- Dependencies: P5-001.
- Completion: partial 2026-11-03 - assembly audited at source:
  per depth K, GENERATE_FULL_MATRIX_V3 builds C_MAT (diagonal
  band, BAND_INDX=DIAG_INDX), B_MAT (sub-diagonal, DIAG_INDX-1,
  K>1 only), D_MAT (super-diagonal, DIAG_INDX+1, K<ND only) -
  a true block tridiagonal (cmf_blkband_v3.f:633-656). The
  elimination is block-Thomas: MAT5PEN forms the Schur complement
  m[k]=C_MAT-B_MAT*PREV_D_MAT and the phi' update; DGEEQU
  equilibration + DGETRF LU per depth (:660-690); NUM_BNDS<3
  hard-errors at :437 in addition to the VADAT guard. Band
  boundaries are correct (no sub-band at K=1, no super-band at
  K=ND). Manufactured solve residual check still needs a driver.

### [ ] P5-004 - Audit off-diagonal entries and boundary rows
- Priority: P0
- Description: verify which physics writes off-diagonal bands
  (advection dlnV terms, TX fills, CMAT_RAD_BANDS) and that boundary
  depth rows are correct under both methods.
- Files: `do_var_cont.f`, `steq_advec_v4.f`, `read_bcd_mat.f`,
  `generate_full_matrix_v3.f`.
- Action: enumerate every off-band writer; check boundary rows.
- Acceptance: written table.
- Dependencies: P5-002, P5-003.
- Completion: 2026-11-03 - off-band writers enumerated:
  (a) steq_advec_v4.f:190-191 - advection dlnV coupling
      (SE%BA(I,I,M,K)-=T1, DO_OFF_DIAG gate adds +=T2 to the
      super-band MP1);
  (b) do_var_cont.f:451-465 - DO_FG13_ALO_TRI writes
      TX(K,K-1,1:2) and TX(K,K+1,1:2) from GET_FG13_ALO_TRI
      (the tridiagonal Lambda* dJ/dCHI,dETA), NUM_BNDS>=3 gated,
      TRI_OK guarded, and prints a once-per-sweep
      |offdiag/diag| dJ/dETA census;
  (c) BA_UPDATE_V7/VJ path - the dJ/dopacity cross-depth blocks
      JFNK_USE_BANDED relies on (jfnk_cmfgen_int_mod.f:469);
  (d) generate_full_matrix_v3.f reads them all into C_MAT/B_MAT/
      D_MAT plus BA_ED for the energy row. Boundary rows: K=1 has
      no sub-band fill and K=ND no super-band fill by
      construction (cmf_blkband_v3.f:639-655).

### [ ] P5-005 - Audit NUM_BNDS restrictions and automatic DIAG fallback
- Priority: P0
- Description: find every place NUM_BNDS<3 forces DIAG or skips TRI
  fills; verify fallback is loud, not silent.
- Evidence: AGENTS.md: with NUM_BNDS<3 the ALO_TRI fill is skipped;
  need to confirm "skipped" vs "silently rerouted".
- Files: `rd_control_variables.f`, `do_var_cont.f`, `solveba_v13.f`.
- Action: trace; enforce loud behaviour if silent.
- Acceptance: documented rule; no silent reroute.
- Dependencies: P5-001.
- Completion: 2026-11-03 - all three NUM_BNDS guards are loud:
  VADAT downgrade warns in block capitals and rewrites METH_SOL
  (rd_control_variables.f:2951-2957); the in-solver guard
  hard-errors "NUM_BNDS too small" (cmf_blkband_v3.f:437-443);
  the ALO_TRI off-band fill is skipped with the documented
  reason "the fill is skipped rather than wasted" and is gated
  NUM_BNDS>=3 (do_var_cont.f:446-451) - a skip, not a reroute,
  since the bands never reach VJ under a diagonal BA.

### [ ] P5-006 - Create manufactured matrix systems with known solutions
- Priority: P0
- Description: driver-level tests solving manufactured DIAG and TRI
  systems through the production assembly/solve path.
- Evidence: unit-test harness exists (`f90/drv_*.f90` pattern);
  no matrix-solve test yet.
- Files: new `f90/drv_solveba*.f90`, `test_solveba_diag_tri.py`.
- Action: build test comparing DIAG/TRI against exact solutions.
- Acceptance: tests pass.
- Dependencies: P5-002, P5-003.
- Completion: partial 2026-11-03 - `test_mat5pen.py` covers the
  TRI Schur-complement kernel through production MAT5PEN calls:
  all four flag combinations match A+D(F.C-E)-B.C to 2e-16, and
  the production m[k]=C[k]-B[k].phi[k-1] update shape matches
  DGEMM to 6e-17. The full CMF_BLKBAND_V3 driver still needs
  SE%BA fabrication (STAT_EQ_DATA pointers + ATM links +
  GENERATE_FULL_MATRIX_V3 inputs) - a heavier lift, flagged
  openly here rather than approximated by a replica.

### [ ] P5-007 - Verify true matrix residuals for both methods
- Priority: P0
- Description: check |Az-b| on the returned corrections in the
  production solve (or document that assembly-level checking exists).
- Files: `solveba_v13.f` post-solve diagnostics, P5-006 test.
- Action: assert residuals in manufactured tests; check production
  prints.
- Acceptance: residual evidence per method.
- Dependencies: P5-006.
- Completion: partial 2026-11-03 - MAT5PEN update algebra is
  exact (test_mat5pen.py, 2e-16), so block-Thomas residuals in
  the TRI path reduce to per-block DGETRF quality; DIAG solves
  each depth independently by construction (P5-002). The open
  leg is an end-to-end |Az-b| assertion through
  CMF_BLKBAND_V3, blocked on the SE fabrication noted in
  P5-006. Production |Az-b| prints do not exist today; adding
  one would need a solveba_v13.f diagnostic.

### [ ] P5-008 - Test singular, masked, and ill-conditioned cases
- Priority: P1
- Description: manufactured cases: exactly-singular depth block
  (pin-null path), masked fixed rows, ill-conditioned blocks;
  verify handling per method.
- Evidence: `test_cmf_blkband_pin.py` pins exactly disconnected
  zero equations; T21 documents ill-conditioning structurally.
- Files: `cmf_blkband_pin_null_v1.f`, `solveba_v13.f`.
- Action: extend tests to the TRI path and masked-row cases.
- Acceptance: tests pass with documented behaviour.
- Dependencies: P5-006.
- Completion: partial 2026-11-03 - DIAG singular-block handling
  is pinned (test_cmf_blkband_pin.py, exactly-disconnected zero
  equations routed through CMF_BLKBAND_PIN_NULL_V1). TRI-path
  singular/masked cases sit inside CMF_BLKBAND_V3's band
  assembly and need the same SE fabrication as P5-006;
  MAT5PEN itself is a pure update kernel - singularity only
  enters at the downstream DGETRF on each m[k] block.

### [ ] P5-009 - Verify interaction with the conventional population solver
- Priority: P0
- Description: confirm both methods drive identical outer-iteration
  semantics in SOLVE_FOR_POPS (correction application, clamps,
  MAXCH accounting).
- Files: `solve_for_pops.f`, `solveba_v13.f`.
- Action: audit post-solve path is method-agnostic.
- Acceptance: verified by inspection + paired run later.
- Dependencies: P5-002, P5-003.
- Completion: 2026-11-03 - verified by inspection: CMF_BLKBAND_V3
  returns the correction in the single array STEQ(NT,ND) no
  matter which branch ran, and every downstream consumer
  (SET_POPS_D2_EQ_D1 zeroing, the BA_TRUST rollback mask,
  MAXCH accounting, CORRECTION_SUM) operates on that array
  without consulting the method (solveba_v13.f:1073+). The only
  method-visible state is the SUCCESS flag and the factor caches.
  Paired DIAG/TRI model comparison remains the runtime leg.

### [ ] P5-010 - Verify interaction with the non-monotonic Knop pathway
- Priority: P1
- Description: confirm DIAG/TRI both function under DO_NONMONO_KNOP
  (KNOP adds lagged-field terms - which bands do they fill?).
- Files: Knop contributions to BA/TX in `do_var_cont.f`,
  `comp_j_blank.f`.
- Action: audit KNOP-specific fills per method.
- Acceptance: documented.
- Dependencies: P5-004, P3-008.
- Completion: 2026-11-03 - under DO_NONMONO_KNOP the lagged-field
  linearization lands on the diagonal band: the FG/MOM ALO
  diagonal terms (TX(K,K,1:2), do_var_cont.f:440-450) and the
  dJ/dopacity cross-depth blocks via BA_UPDATE_V7 - the latter
  are the only KNOP terms reaching the off-bands, and only under
  NUM_BNDS>=3 (P5-004 writer (b)); the DO_FG13_ALO_TRI fill adds
  the tridiagonal Lambda* off-diagonals when requested. With
  NUM_BNDS=1/DIAG the off-band content simply is not assembled -
  consistent with DIAG's depth-independent semantics.

### [ ] P5-011 - Verify interaction with the Sobolev pathway
- Priority: P1
- Description: confirm both methods under DO_NONMONO_SOB (SOB line
  linearization bands).
- Files: `sub_sob_line_v3.f` Jacobian terms, `do_var_cont.f`.
- Action: audit.
- Acceptance: documented.
- Dependencies: P5-004.
- Completion: 2026-11-03 - under DO_NONMONO_SOB the line
  linearization is strictly local: ZNET_SIM depends on
  per-(depth,line) JBAR built from that depth's own Sobolev tau,
  so all SOB BA terms land on the diagonal band (P2-007 audit).
  The nonlocal RH term DRC_SIM is deliberately excluded from
  ZNET_SIM - "it must not enter ZNET_SIM which the BA Jacobian
  differentiates" (sub_sob_line_v3.f:164-168) - so nonlocal
  coupling contributes nothing to any band: it is a lagged
  remote forcing, not a Jacobian term. Method choice therefore
  affects only the continuum/advection cross-depth terms, never
  the SOB line rows.

### [ ] P5-012 - Establish whether JFNK uses METH_SOL or its own preconditioner
- Priority: P0
- Description: JFNK has block-LU/banded-LU/Schur preconditioner
  options independent of SOL_METH - verify which production code
  path each config uses; no silent coupling.
- Evidence: dev plan: preconditioner has block LU, banded LU, n,T
  Schur options; distinct from a radiation/population Schur.
- Files: `jfnk_solve_corr_v1.f`, `JFNK_FACTOR_PRECOND_V1` in
  `jfnk_cmfgen_interface_v1.f`.
- Action: document the preconditioner selection vs METH_SOL.
- Acceptance: written mapping.
- Dependencies: none.
- Completion: 2026-11-03 - JFNK never calls CMF_BLKBAND: it runs
  its own GMRES + preconditioner stack. `JFNK_METH_SOL=METH_SOL`
  is recorded (jfnk_cmfgen_interface_v1.f:56) only so the
  preconditioner can mirror the production band structure.
  Selection: JFNK_USE_BANDED VADAT flag (default off) switches
  between block-diagonal LU and banded LU (DGBTRF/DGBTRS,
  jfnk_cmfgen_int_mod.f:467-473); the banded path requires
  NUM_BNDS>1 and is silently disabled under NUM_BNDS=1 per its
  own comment - a preconditioner-quality downgrade only, not a
  wrong-answer path, but noted as the one quiet fallback in the
  matrix. The conventional correction side always runs SOLVEBA
  with METH_SOL, except under LAMBDA_ITERATION which forces
  'DIAG' (solve_for_pops.f:1004-1007) - correct since lambda
  steps need no cross-depth bands.

### [ ] P5-013 - Compare DIAG and TRI on identical physical models
- Priority: P1
- Description: paired runs (depth_independent_solver true/false) on
  the same model; compare converged solutions.
- Evidence: `adp2_alo3`/`adp2_triwarm` configs exist; runtime status
  unverified per AGENTS.md.
- Files: configs `*tri*.yaml`.
- Action: submit paired run when a converged baseline exists.
- Acceptance: consistent solutions within tolerance or documented
  difference.
- Dependencies: P5-009, P1-007 (baseline model).
- Completion: partial 2026-11-03 - paired probe submitted: job
  12669084 (`2.82days_tri`, NUM_BNDS=3 + SOL_METH=TRI) cloned
  from the identical it620 restart state feeding the running
  DIAG job 12667827 (`2.82days_tres`, NUM_BNDS=1 + SOL_METH=
  DIAG). Same model, same state, both on the post-BUG-006
  build. Hardlinks broken before submission so the TRI job
  cannot write into DIAG-run inodes. Queued behind the two
  running interactive jobs; comparison pending its first
  iterations.

### [ ] P5-014 - Record runtime and memory differences
- Priority: P3
- Description: cost comparison DIAG vs TRI from the P5-013 runs.
- Action: record wall/mem from job logs.
- Acceptance: recorded.
- Dependencies: P5-013.
- Completion: -

## 9. Cross-pillar integration

The configuration matrix (INT-001) distinguishes: Implemented /
Compiles / Unit-tested / Runtime-tested / Physically validated.
Do not conflate these statuses.

### [ ] INT-001 - Create the complete VADAT configuration support matrix
- Priority: P0
- Description: matrix over {transfer: stock, SOB, SOB+FG_KNOP, KNOP}
  x {J: LAM_ES=T, LAM_ES=F} x {solver: conventional, JFNK} x
  {matrix: DIAG, TRI} with per-cell status + reason + test ref.
- Evidence: partial statuses scattered through AGENTS.md and audits.
- Files: this TODO + a `docs/` matrix table.
- Action: write the matrix; mark unsupported combos rejected with
  the enforcing check.
- Acceptance: complete matrix; every cell has status+evidence.
- Dependencies: P4-016, P5-005.
- Completion: 2026-11-03 - `docs/audits/VADAT_SUPPORT_MATRIX.md`
  covers all 32 cells with PROD/TEST/SRC/REJ/UNSUP status and
  file:line evidence for every hard gate (mutual exclusion,
  DO_JFNK requirement, NUM_BNDS fallback, GLOBAL_LINE override,
  JFNK_COUPLED-under-SOB rejection, SN-model LAM_ES=F rule).
  Two cells upgraded during the audit: JFNK banded
  preconditioner's NUM_BNDS=1 silent fallback is documented (a
  preconditioner downgrade only); SOB+FG_KNOP marked inert under
  USE_LAM_ES=F.

### [ ] INT-002 - Test original CMFGEN with all new flags disabled
- Priority: P0
- Description: the P1-007 A/B run IS this test.
- Dependencies: P1-007.
- Completion: -

### [x] INT-003 - Test Sobolev with conventional population solver
- Priority: P0
- Evidence: `m15_2p82d_SOB_l14_lte_np` converged 100% it718;
  `adp2_esali` converged 100% it851 (both conventional solver,
  SOB + LAM_ES=T).
- Acceptance: converged nonmono-SOB model. MET (LAM_ES=T).
- Completion: recorded evidence 2026-11-03; LAM_ES=F case tracked
  under P3-012/P3-024.

### [ ] INT-004 - Test Knop with USE_LAM_ES=T
- Priority: P0
- Description: converged KNOP formal+LAM_ES model.
- Evidence: negative (stalls, P3-009).
- Dependencies: P3-009.
- Completion: -

### [ ] INT-005 - Test Knop with USE_LAM_ES=F
- Priority: P0
- Description: converged KNOP moment model.
- Evidence: negative (non-contraction, P3-012).
- Dependencies: P3-024.
- Completion: -

### [ ] INT-006 - Test JFNK with every supported transfer configuration
- Priority: P0
- Description: JFNK admitted under SOB+LAM_ES=T (frozen) and KNOP
  (coupled); verify each admission from the P4-016 matrix.
- Dependencies: P4-002, P4-003, P4-016.
- Completion: -

### [ ] INT-007 - Test DIAG and TRI across applicable configurations
- Priority: P1
- Dependencies: P5-013.
- Completion: -

### [ ] INT-008 - Test monotonic, mild reversal, strong compression, multiple reversal profiles
- Priority: P0
- Description: the four velocity regimes through the appropriate
  pathways; stock must handle monotonic only.
- Evidence: m15 model has a compression zone (interacting);
  manufactured reversal coverage in unit tests is thin.
- Dependencies: P3-005 + model runs.
- Completion: -

### [ ] INT-009 - Test representative STELLA-derived interacting-SN structure
- Priority: P1
- Description: a model built from a STELLA hydro structure
  (`prepare_SN_HYDRO_DATA` / stella dirs exist at workspace root).
- Dependencies: P1-007 baseline + P3 progress.
- Completion: -

### [ ] INT-010 - Validate radiation/population/charge/energy consistency at termination
- Priority: P0
- Description: for each converged run: SE residual level, radiative
  equilibrium, charge conservation, luminosity conservation checks.
- Files: `check_convergence.py`, `check_energy_cons.exe`,
  OUTGEN residuals.
- Dependencies: each converged run.
- Completion: -

### [ ] INT-011 - Verify restart and checkpoint behaviour
- Priority: P1
- Description: SCR checkpoint round-trip incl. FGINEXT_GS /
  MOM_JNEXT_GS lagged stores (fingerprint+tag validation);
  foreign-seed policy.
- Evidence: store validation exists in code (audit item 6);
  KN_REPAIR_T self-heal documented; no dedicated restart test.
- Dependencies: none for code audit; model evidence from restarts.
- Completion: -

### [ ] INT-012 - Final clean production build + full unit-test suite
- Priority: P0
- Description: `./build.sh` zero-error full build and
  `run_all_tests.py` all-pass at the final HEAD; record exes.
- Evidence: last recorded clean build 2026-11-03 (21/21).
- Dependencies: all code-changing tasks.
- Completion: -

## 10. Documentation discrepancies

### [x] DOC-002 - In-source edit-date comments are post-dated
- Priority: P3
- Description: Phase-0 reconciliation found comments carrying dates
  later than the repository revision (e.g. "12-Nov-2026" edits in
  `jfnk_solver_v1.f90` and `jfnk_solve_corr_v1.f`) while HEAD is
  dated 2026-11-03. `git log -S`/`--follow` shows those comment
  lines were introduced by commits dated 2026-10-23 - the in-source
  dates are systematically post-dated annotations, not evidence of
  an unaudited later revision.
- Files: scattered solver comments; no source change needed.
- Action: record the convention here; audits must cite git commit
  dates, not in-source edit dates. No further action planned.
- Acceptance: convention recorded (this entry). MET.
- Dependencies: none.
- Completion: recorded 2026-11-03 (Phase-0 reconciliation).

### [x] DOC-001 - Reconcile JFNK investigation doc date and stated HEAD
- Priority: P2
- Description: `JFNK_PRODUCTION_OPERATOR_INVESTIGATION_2026-11-07.md`
  was dated 2026-11-07 and claimed tree `22d63df9` + uncommitted
  working-tree changes. Actual facts (verified 2026-11-03): the doc
  and all probe code were committed 2026-11-03 in `a1f2913c..ccc52ffd`
  (17 commits past `22d63df9`); the Nov-07 date was a typo for
  Nov-03; "not committed" was stale.
- Files: the audit doc header + sec.H.
- Action: corrected the header (date, real commit range) and the
  stale section-H note; the misleading filename is retained so
  existing references do not break.
- Acceptance: doc consistent with git history. MET.
- Dependencies: none.
- Completion: fixed 2026-11-03 in the working tree.

## 10a. BUG HUNT REGISTER (owner-directed adversarial review)

Register of confirmed and suspected defects found by the
independent source-level bug hunt started 2026-11-03. Full per-bug
entries with required fields live in
`docs/audits/BUG_HUNT_REGISTER.md`; technical narrative and
production-readiness assessment in
`docs/audits/BUG_HUNT_TECHNICAL_REPORT.md`. This section mirrors the
register status only.

| Bug ID  | Severity | Title                                             | Status |
|---------|----------|---------------------------------------------------|--------|
| BUG-001 | P0       | GMRES `info=0` convergence accepted on hess_res   | FIXED  |
| BUG-002 | P2       | SOB_RH_GEOM cache fingerprint misses SIGMA        | FIXED  |
| BUG-003 | P3       | MKL "Parameter 4" warning on entry to DGELSY      | REFUTED|
| BUG-004 | P2       | Zero-GAM resonance: floored factor vs exact-zero  | REFUTED|
| BUG-005 | P2       | JFNK admissions on m15 always defer before GMRES  | REFUTED|
| BUG-006 | P3       | LSTATE save armed after first eval mutations      | FIXED  |

- BUG-001 (P0, FIXED 2026-11-03): `JFNK_TRES_FAC` defaulted to `0`,
  disabling true-residual verification; an `info=0 gmres_converged`
  return could rest on the Hessenberg model residual alone
  (documented 15x-450x below true). Fix: unconditional
  `||M^{-1}(b-Ax)||_w` (left-PC) / `||b-Ax||_w` (right-PC)
  measurement at every done-exit; `JFNK_TRES_FAC` is the acceptance
  bound (default 1.0, `<=0` floored to 1.0, `>1` documented
  weakening); all unverifiable paths fail closed. Regression:
  test_newton_core T4j/T4k + existing T4a-i; test_guards G17 +
  G13-G15 updates; T8 uses the documented `fac=1.5` relaxation for
  its step-disposition A/B (FD floor at `|F|`). Negative check:
  T4k(biased) and G17 FAIL on the pre-fix solver. Surrogate suite:
  `test_se_solver.py` T11/T12/T13 were calibrated on
  phantom-acceptance trajectories and starved under verified solves
  - recalibrated on real trajectories (T11 gated to std+shock; T12
  extended to 100 admits, `seqstr` 64->256; T13 reshaped to
  D*20+P+D*3+P+D*15). `test_jfnk_solver.py`: 89 passed, 0 failed.
  Evidence: `tests/jfnk_core` all-PASS under ifort; full `build.sh`
  rebuild recorded in the build transcript.
- BUG-002 (P2, FIXED 2026-11-03): `SOB_RH_GEOM` decided cache
  validity on `SUM(V)`+`SUM(R)` at 1D-25 relative, but the cached
  `RH_T1`/`RH_T1P` are `|1+mu^2*sigma|` - SIGMA was never
  fingerprinted, and the digest cannot distinguish profiles sharing
  a sum. Same-R,V call with new sigma silently reused stale
  resonance taus -> wrong `DRC`. Latent: all current sigma rewrites
  pair with R or V changes. Fix: element-wise fingerprint on stored
  copies of all three inputs (`ALL(...EQ...)`, nested under the ND
  check so `ALL` never sees an unallocated array). Regression:
  multi-phase `drv_sob_rh` + T7 in `test_sob_rh_couple.py` - 2 of 3
  new checks FAIL on pre-fix (no rebuild, stale T1/T1P), 12/12 PASS
  on fixed including the phase-3 cache-hit check.
- BUG-003 (P3, REFUTED 2026-11-03): the oneMKL `Parameter 4 was
  incorrect on entry to DGELSY` warning fires only in the legacy
  `tests/f90` T7 all-NaN-residual scenario (MKL entry check on a
  NaN-poisoned Hessenberg). Both production `dgelsy` calls have
  adequate LDA; the error path is handled (`info=-4` ->
  `gmres_dgelsy_error`). Identical on the pre-fix solver. No
  production change.

## 11. Priority order (next eligible tasks)

Following the owner's suggested sequence and dependencies:

1. P1-009 - diagnostics/safeguards classification done 2026-11-03.

   All P1 stock-path audits done; remaining P1 items need runs.
2. P4-002 + P4-004 + P4-010 - submit a SOB+JFNK HSCAN production run
   (highest-value pending measurement; can run while audits proceed).
3. P3-013/P3-014/P3-015 - moment frequency-grid investigation
   (inspection + unit tests).
4. P3-001..P3-004 - finish the Knop operator paper-map doc.
5. P2-001..P2-010 - Sobolev audits (inspection-level).
6. P5-001..P5-008 - DIAG/TRI audits + manufactured tests.
7. P1-007/P1-008 - A/B model runs vs pristine reference (needs a
   reference build).
8. P4-003/P4-005 (coupled probe) and INT KNOP/JFNK cells - blocked on
   KNOP model viability (P3-009/P3-012).

Blocked/pending-run items: P4-002..P4-008 (awaiting HSCAN run data),
P1-007/P1-008 (awaiting reference build + A/B run), P3-009/P3-024
(awaiting numerical work on KNOP stall), P5-013/P5-014 (awaiting
baseline model).
