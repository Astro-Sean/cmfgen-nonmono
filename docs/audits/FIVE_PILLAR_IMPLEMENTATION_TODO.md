# CMFGEN Five-Pillar Implementation TODO

This file is the single source of truth for the remaining development
work on the `cmfgen_jfnk` non-monotonic fork. It was created under the
project owner's TODO directive (`/ptmp/sbrennan/SN2024pba/TODO.md`) and
is updated in place as tasks complete. Do not create competing TODO
files; do not delete completed tasks.

## 1. Project status

- Repository: `git@github.com:Astro-Sean/cmfgen-nonmono.git` (private)
- Current branch: `main`
- Current HEAD: `ccc52ffd84f5fe4527b55a211c72c20b4455e21d`
  (= `origin/main`, verified by fetch on 2026-11-03; matches the
  revision quoted in the owner's TODO directive)
- Working tree: clean except untracked `unit_tests/systems/_pv.txt`
  (test artifact, not source)
- Pristine baseline: `../cmfgen_reference/cur_cmf`, extracted from
  `cur_cmf_18jun25.tar.gz`, md5 `1d5a7a168e77329dd97c3019bcad3c09`
  (re-verified 2026-11-03, matches AGENTS.md)
- Last updated: 2026-11-03
- Last full build: 2026-11-03 ~20:44, `./build.sh` clean, all three
  exes present; unit suite 21/21 files pass (`build_audit_fix.out`)
- Overall development status: all five pillars have production code and
  unit-test coverage; Sobolev+LAM_ES=T path is model-validated on m15;
  KNOP and JFNK paths are implementation-complete but not
  model-validated; moment (LAM_ES=F) path documented as
  non-contracting on m15
- Completed tasks: 7 (P1-001, P1-002, P2-012, P4-001, P4-012,
  INT-003, DOC-001)
- Outstanding tasks: 84 (22 in progress with partial evidence,
  62 not started)
- Blocked tasks: 0 (run-dependent tasks are sequenced, not blocked;
  see section 11)
- Current active task: P1-003
- Next task: P1-004

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

### [~] P1-003 - Audit shared module state and SAVE variables with new flags off
- Priority: P0
- Description: verify no module-level or SAVE state written by
  nonmono/JFNK code paths leaks into the stock path when flags are F.
- Evidence: `CONVENTIONAL_SOLVER_AUDIT_2026-11-03.md` traced the
  conventional call graph; `JFNK_PRODUCTION_OPERATOR_INVESTIGATION`
  sec.B audited residual save/restore by inspection (production REP
  measurement pending). No dedicated SAVE-variable sweep exists.
- Files: `control_variable_mod.f`, `mod_cmfgen.f`,
  `fg_j_cmf_nonmono_v1.f`, `mom_j_cmf_nonmono_v1.f`,
  `jfnk_*` modules, `comp_j_blank.f`, `cmfgen_sub.f`,
  `solve_for_pops.f`; any `SAVE` attributes in touched files.
- Action: enumerate module variables and SAVE locals in the 109
  modified files; classify each as (a) written only under new flags,
  (b) written always but read only under new flags, (c) live in the
  stock path - class (c) needs a justification or a fix.
- Acceptance: a written classification table; zero unexplained
  class-(c) writes that alter stock results.
- Dependencies: P1-002.
- Completion: -

### [~] P1-004 - Verify original formal-transfer dispatch with all new options off
- Priority: P0
- Description: prove `FG_J_CMF_V13`/`MOM_J_CMF_V11` and the stock
  boundary/Eddington paths execute identically to upstream when
  `DO_NONMONO_SOB=F`, `DO_NONMONO_KNOP=F`, `DO_JFNK=F`.
- Evidence: `CONVENTIONAL_SOLVER_AUDIT_2026-11-03.md` sec.1: dispatch
  predicate verified - nonmono branches are flag-gated. Bit-identity
  covered at routine level for interp and SOB helpers
  (`test_interp_guards.py`, `test_sob_tau_floor.py` sigma>=0).
  Missing: an executable end-to-end stock-vs-reference J comparison.
- Files: `comp_j_blank.f` (FG_J_CMF_NONMONO_V1 :844 vs FG_J_CMF_V13
  :857; MOM_J_CMF_NONMONO_V1 :1032 vs MOM_J_CMF_V11),
  `fg_j_cmf_v13.f`, `mom_j_cmf_v11.f`, `do_var_cont.f`.
- Action: audit each edited formal-solve routine's gating predicate in
  the diff inventory; add a stock-dispatch unit test if a cheap
  driver exists; note `fg_j_cmf_v13.f:1069` SIGMA clamp is upstream.
- Acceptance: every modified formal-path file shows its changes dead
  or identity under flags-off; dispatch audit written.
- Dependencies: P1-003.
- Completion: partial (predicate audit done in conventional audit;
  per-file gating checklist open).

### [~] P1-005 - Verify original Sobolev and line-transfer dispatch, flags off
- Priority: P0
- Description: prove `SUB_SOB_LINE_V3`, `SOBJBAR_SIM`, `SOBEW*`,
  `INS_LINE*` behave identically to upstream with new options off.
- Evidence: `test_sob_tau_floor.py` asserts bit-identity vs reference
  for sigma>=0 and sigma=-0.5; `test_sob_rh_couple.py` asserts zero
  pairs under homology; AGENTS.md documents SIGMA_TRUE is only passed
  under `DO_NONMONO_SOB`.
- Files: `sobjbar_sim.f`, `sobew.f`, `sobew_grad_v2.f`,
  `sobew_grad_v3.f`, `sub_sob_line_v3.f`, `sob_rh_couple.f`,
  `ins_line_v6.f`, `ins_line_obs_v5.f`, `cmfgen_sub.f` (SIM_LINES).
- Action: confirm each SOB entry point's flag plumbing; confirm the
  empty-level ZNET/BETAC guards (29-Oct, not flag-gated) are exact
  no-ops for finite emissivity; record remaining unaudited callers.
- Acceptance: flags-off bit-identity or documented benign difference
  per routine; written dispatch table.
- Dependencies: P1-003.
- Completion: partial (tau-floor and RH-couple tests give routine
  level evidence; `sub_sob_line_v3` guard analysis open).

### [~] P1-006 - Verify original population/statistical-equilibrium solver
- Priority: P0
- Description: prove `SOLVE_FOR_POPS`+`SOLVEBA_V13` stock iteration is
  unchanged when `DO_JFNK=F` and nonmono flags off.
- Evidence: `CONVENTIONAL_SOLVER_AUDIT_2026-11-03.md` traced the
  Lambda/BA call graph and found ordering consistent with upstream.
  Known non-gated changes: `cmf_blkband_pin_null_v1` row pinning
  (exact-0 equations only), LUM_DIV telemetry, BAACC_RES2_MIN gating
  fix (W1), empty-level ZNET guards - all documented as stock-path
  defect fixes in `ROBUSTNESS_AUDIT_2026-10.md` and AGENTS_HISTORY.
- Files: `solve_for_pops.f`, `solveba_v13.f`, `cmf_blkband_v3.f`,
  `cmf_blkband_pin_null_v1.f`, `fiddle_pop_corrections_v2.f`,
  `generate_full_matrix_v3.f`, `fixpop_in_ba_v3.f`.
- Action: list every non-flag-gated edit inside these files (from the
  inventory diff); confirm each is either an upstream-defect fix with
  a test, or flag-gated; write the list.
- Acceptance: written table of non-gated stock-path edits with
  justification + test reference for each.
- Dependencies: P1-002.
- Completion: partial (audits cover much of this; consolidated table
  open).

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

### [~] P1-009 - New diagnostics/safeguards do not change legacy results
- Priority: P1
- Description: verify counters, census prints, watchdogs, and guard
  conversions added by this fork are observation-only on the stock
  path (no result dependence).
- Evidence: `ROBUSTNESS_AUDIT_2026-10.md` lists each guard's intent;
  W6 fail-closed conversions documented; f-cap census counters are
  print-only by inspection. No consolidated proof.
- Files: `control_variable_mod.f` counters, `solve_for_pops.f`
  (LUM_DIV, GS_*), `cmfgen_sub.f` (DBG-T, f-cap census),
  `solveba_v13.f` (BA quarantine - gated on DO_BA_RES_ACC+DO_BA_QUAR).
- Action: classify each diagnostic as write-only-print vs
  state-affecting; for state-affecting ones cite the VADAT gate.
- Acceptance: written table; zero ungated result-affecting
  diagnostics.
- Dependencies: P1-003.
- Completion: partial (audit lists exist; consolidated classification
  open).

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
- Completion: partial (formula known and tested; written derivation
  open).

### [ ] P2-003 - Audit the denominator floor MAX(|1+sigma*mu^2|, -0.2*(1+sigma), 1E-4)
- Priority: P0
- Description: quantify when the floor activates and how much it
  changes tau_sob vs the unfloored expression.
- Evidence: unit test pins the formula; no runtime activation census
  exists.
- Files: the four tau-floor files; OUTGEN counters if added.
- Action: analytic quantification (mu,sigma plane region where floor
  binds) plus - if a run is available - an opt-in activation counter
  (diagnostic only, no result change); decide whether the counter is
  worth code or whether the analytic report suffices.
- Acceptance: documented activation region + bounded relative error
  it introduces; any counter must be print-only.
- Dependencies: P2-002.
- Completion: -

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
- Completion: partial.

### [ ] P2-005 - Audit resonance finding and ordering
- Priority: P0
- Description: verify how SOB line resonances are located/ordered
  along a ray for non-monotonic V, and that ordering logic handles
  re-entrant velocity ranges.
- Evidence: `CONVENTIONAL_SOLVER_AUDIT_2026-11-03.md` sec.1: repeated
  resonances along a ray are handled inside the ray optical-depth
  integration; no dedicated finding/ordering audit exists.
- Files: `sub_sob_line_v3.f`, `ins_line_v6.f`,
  `cmfgen_sub.f` (SIM_LINES resonance-zone opacity), `sob_rh_couple.f`.
- Action: trace the resonance search code; document ordering
  assumptions; identify where monotonic-V assumptions survive in
  indexing or delta-function placement.
- Acceptance: written audit; any monotonic-ordering assumption either
  justified or filed as a defect.
- Dependencies: P2-001.
- Completion: -

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
- Completion: -

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
- Completion: -

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
- Completion: -

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
- Completion: partial.

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
- Completion: -

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

### [ ] P3-005 - Validate multiple velocity reversals and strong compression
- Priority: P0
- Description: numerical validation of the formal solver on
  manufactured profiles with 2+ reversals and sigma<<-1 compression.
- Evidence: `test_knop_ray.py` covers one compressing band (NI=8,
  NCF=6) - sweep convergence, store agreement ~3e-8, forced
  stall/cap exits; multiple reversals not covered.
- Files: `test_knop_ray.py`, `unit_tests/f90/drv_knop_ray*.f90`.
- Action: extend the manufactured ray test to multi-reversal and
  strong-compression V profiles.
- Acceptance: new manufactured cases pass with documented limits.
- Dependencies: P3-001..P3-004.
- Completion: -

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

### [ ] P3-013 - Investigate MOM_DLN_MIN and near-coincident frequency points
- Priority: P0
- Description: MOM_DLN_MIN floors the log-frequency spacing in
  advection coefficients; bound-free edge doublets (~2E-11 dLOG_NU)
  on the continuum grid are the documented trigger - quantify effect
  on coefficients.
- Evidence: AGENTS.md: dominant advection-coefficient inflation comes
  from edge doublets; dev plan Phase 2 scopes this investigation
  (trace every MOM_DLN_MIN use, activation counts, units).
- Files: `mom_j_cmf_nonmono_v1.f` (MOM_DLN_J/MOM_DLN_H),
  frequency grid construction.
- Action: census every MOM_DLN_MIN use site; build manufactured
  grids (ordinary, small spacing, doublets, invalid); measure the
  coefficient change per activation.
- Acceptance: documented activation conditions and coefficient
  impact; no removal without derived replacement.
- Dependencies: P3-010.
- Completion: -

### [ ] P3-014 - Derive the small-frequency-spacing limiting equation
- Priority: P0
- Description: derive the limiting discrete equation as dln_nu -> 0
  and design the correct treatment (what the advection term should
  degenerate to).
- Evidence: dev plan Phase 2 acceptance requires a derived limiting
  form before any replacement for MOM_DLN_MIN.
- Files: `mom_j_cmf_nonmono_v1.f`, Knop Paper I equations.
- Action: pen-and-paper derivation recorded in docs; then evaluate
  whether MOM_DLN_MIN matches the limit.
- Acceptance: derived limit + comparison vs current guard.
- Dependencies: P3-013.
- Completion: -

### [ ] P3-015 - Build frequency-integrated conservation tests
- Priority: P0
- Description: manufactured tests asserting the discrete moment
  equations conserve the relevant frequency-integrated quantity
  across grid-spacing regimes incl. near-degenerate grids.
- Evidence: dev plan notes existing tests do not establish
  frequency-integrated conservation for nearly degenerate grids.
- Files: `test_mom_solve.py`, new driver if needed.
- Action: add conservation assertions vs the discretized equations
  (not a weaker proxy).
- Acceptance: tests pass on ordinary grids and document behaviour on
  degenerate grids.
- Dependencies: P3-013, P3-014.
- Completion: -

### [ ] P3-016 - Investigate the 1 + sigma*K_ON_J = 0 resonance
- Priority: P0
- Description: the moment equations have a resonance where
  1+sigma*K_ON_J crosses zero (compression); trace it through
  closure, angular terms, boundaries, iteration.
- Evidence: dev plan Phase 2 question; the K_ON_J factors come live
  from ray moments; no dedicated analysis yet.
- Files: `mom_j_cmf_nonmono_v1.f`, `get_moms_non_rel.f`.
- Action: locate the coefficient, derive behaviour on each side,
  manufactured case straddling the resonance.
- Acceptance: resonance handling described mathematically and
  tested; no larger regularisation accepted as fix.
- Dependencies: P3-010.
- Completion: -

### [ ] P3-017 - Quantify the effect of MOM_GAM_FAC_MIN
- Priority: P1
- Description: measure where MOM_GAM_FAC_MIN binds and its effect on
  the advection coefficient and solution.
- Evidence: named in dev plan as a guard not to remove without
  derived replacement; no census.
- Files: `mom_j_cmf_nonmono_v1.f`.
- Action: activation census + coefficient-impact analysis like
  P3-013.
- Acceptance: documented.
- Dependencies: P3-010.
- Completion: -

### [ ] P3-018 - Audit MOM_LAG_BND
- Priority: P1
- Description: audit the lagged boundary treatment key MOM_LAG_BND -
  which boundary values lag, gating, defaults.
- Evidence: key exists in rd_control_variables (test_vadat_keys
  coverage); no dedicated audit found.
- Files: `mom_j_cmf_nonmono_v1.f`, `rd_control_variables.f`.
- Action: trace uses; document.
- Acceptance: written audit.
- Dependencies: P3-010.
- Completion: -

### [ ] P3-019 - Audit MOM_ABS_FLR
- Priority: P1
- Description: audit the MOM_ABS_FLR floor - what it bounds, when it
  activates, its physical justification per the no-silent-clamp
  policy.
- Evidence: none found.
- Files: `mom_j_cmf_nonmono_v1.f`.
- Action: trace, classify per clamp policy (physical vs numerical
  intermediate), quantify activations.
- Acceptance: written audit; unjustified clamps filed as defects.
- Dependencies: P3-010.
- Completion: -

### [ ] P3-020 - Audit MOM_H_IMPL
- Priority: P1
- Description: audit the MOM_H_IMPL option - implicit H treatment
  branch, correctness of the implicit solve vs explicit.
- Evidence: none found.
- Files: `mom_j_cmf_nonmono_v1.f`, `rd_control_variables.f`.
- Action: trace branch; verify both on manufactured case.
- Acceptance: written audit + test.
- Dependencies: P3-010.
- Completion: -

### [ ] P3-021 - Audit MOM_GAM_PIN
- Priority: P1
- Description: audit MOM_GAM_PIN - pinned GAM values option,
  semantics, defaults, interaction with the advection term.
- Evidence: none found.
- Files: `mom_j_cmf_nonmono_v1.f`, `rd_control_variables.f`.
- Action: trace and document.
- Acceptance: written audit.
- Dependencies: P3-010.
- Completion: -

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
- Completion: -

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
- Evidence: audit sec.D: all exits produce explicit true-residual
  measurement; `JFNK_TRES_FAC` gates acceptance (default 0=off);
  jvp_soft records ltres (lines ~2896-2918 verified); hess_res seen
  15-450x below true_res historically - which is why the gate
  exists.
- Files: `jfnk_solver_v1.f90`.
- Action: remaining gap - confirm TRES_FAC actually armed in
  production runs (part of P4-010); then close.
- Acceptance: source-verified (done) + production config armed.
- Dependencies: P4-010.
- Completion: partial (source-verified; runtime arming open).

### [ ] P4-010 - Test JFNK_TRES_FAC on/off; establish production configuration
- Priority: P0
- Description: decide the production TRES_FAC setting from measured
  true-vs-Hessenberg ratios; ensure a claimed convergence is never
  accepted on the Hessenberg estimate alone.
- Evidence: default 0=off; hess_res vs true_res gap documented.
- Files: `jfnk_solver_v1.f90`, VADAT docs.
- Action: runs with TRES_FAC>0 (e.g. paired with P4-004/P4-005);
  record tres_rej counts; set recommended default (recommendation
  only - changing a default needs owner sign-off).
- Acceptance: measured basis for the production setting.
- Dependencies: P4-004.
- Completion: -

### [ ] P4-011 - Verify inexact-Newton forcing criterion and true-residual acceptance
- Priority: P0
- Description: verify Eisenstat-Walker forcing is in the
  preconditioned (weighted) norm consistent with left
  preconditioning, and the caller gate (`true_res < 0.9*ref_norm`,
  `true_res > |F|` skip) behaves as audited.
- Evidence: audit sec.D verified by source; runtime unverified.
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

### [ ] P4-015 - Verify coupled-J radiation sweep restores all state
- Priority: P0
- Description: verify `JFNK_RAD_SWEEP_V1` save/restore actually
  returns bit-identical state (radiation, line, moment/Knop refresh,
  Eddington-record ACCESS_F, file cursors) after an eval.
- Evidence: audit sec.A: explicit bracket lists audited by
  inspection; REP rows on production will measure it.
- Files: `cmfgen_sub.f:5755` sweep, `jfnk_cmfgen_interface_v1.f`.
- Action: production REP/REPP under coupled admission (P4-003)
  directly tests this; inspect restore lists against mutated arrays
  once more on current source.
- Acceptance: measured + inspection-verified restore.
- Dependencies: P4-003.
- Completion: -

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
- Completion: partial (parser-level verified; matrix doc open).

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
- Completion: partial (dispatched located; end-to-end trace open).

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
- Completion: -

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
- Completion: -

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
- Completion: -

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
- Completion: -

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
- Completion: -

### [ ] P5-007 - Verify true matrix residuals for both methods
- Priority: P0
- Description: check |Az-b| on the returned corrections in the
  production solve (or document that assembly-level checking exists).
- Files: `solveba_v13.f` post-solve diagnostics, P5-006 test.
- Action: assert residuals in manufactured tests; check production
  prints.
- Acceptance: residual evidence per method.
- Dependencies: P5-006.
- Completion: -

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
- Completion: -

### [ ] P5-009 - Verify interaction with the conventional population solver
- Priority: P0
- Description: confirm both methods drive identical outer-iteration
  semantics in SOLVE_FOR_POPS (correction application, clamps,
  MAXCH accounting).
- Files: `solve_for_pops.f`, `solveba_v13.f`.
- Action: audit post-solve path is method-agnostic.
- Acceptance: verified by inspection + paired run later.
- Dependencies: P5-002, P5-003.
- Completion: -

### [ ] P5-010 - Verify interaction with the non-monotonic Knop pathway
- Priority: P1
- Description: confirm DIAG/TRI both function under DO_NONMONO_KNOP
  (KNOP adds lagged-field terms - which bands do they fill?).
- Files: Knop contributions to BA/TX in `do_var_cont.f`,
  `comp_j_blank.f`.
- Action: audit KNOP-specific fills per method.
- Acceptance: documented.
- Dependencies: P5-004, P3-008.
- Completion: -

### [ ] P5-011 - Verify interaction with the Sobolev pathway
- Priority: P1
- Description: confirm both methods under DO_NONMONO_SOB (SOB line
  linearization bands).
- Files: `sub_sob_line_v3.f` Jacobian terms, `do_var_cont.f`.
- Action: audit.
- Acceptance: documented.
- Dependencies: P5-004.
- Completion: -

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
- Completion: -

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
- Completion: -

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
- Completion: -

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

## 11. Priority order (next eligible tasks)

Following the owner's suggested sequence and dependencies:

1. P1-003/P1-004/P1-005/P1-006/P1-009 - stock-path audits (in progress;
   all inspection-level, no runs needed).
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
