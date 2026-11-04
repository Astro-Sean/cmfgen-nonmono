# Species/Population Audit - sequential task log

Started: 07-Nov-2026. Trigger: reported O V overpopulation in the
m15_2p82d validation models. Rules: no species-specific fixes; audit
before changing solver behaviour; general implementation only.

## Task 1 - Map population data model - DONE

- Files read: `mod_cmfgen.f` (MODEL_ATOM_DATA, SPECIES_*, ION_ID,
  POP_SPECIES), `steq_data_mod.f` (SE type), `cmfgen.f` (EQXzV
  allocation), `poptoion.f` (row layout), `eval_lte_v5.f` (DXzV
  aliasing), `rd_f_to_s_ids.f` (map validation), `cmfgen_sub.f`.
- Established: per-ion block = NXzV SL rows; DXzV(ID) shares row
  EQXzV+NXzV=EQXzV(ID+1) = next ion's density (deliberate); terminal
  slot = 1 row = EQ_SPECIES = terminal-ion density; POPS(NT-1)=ED,
  POPS(NT)=T. F_TO_S input-validated [1,NXzV].
- Tests: none (read-only). Concerns: none.

## Task 2 - Oxygen ionization chain (POPOXY, m15_2p82d) - DONE

- O I-VI model atoms (48/123/86/152/130/31 full levels); O VII terminal
  slot. OV fraction: 0.715 (d1) / 0.934 (d35) / 0.0042 (d70) in
  VAL_sob_nonloc; sob_jfnk/diag/tri agree to 3e-4. knop_formal run is
  O III-dominated (different thermal state) - tracked as open item.
- All 9100 OV values finite, non-negative; per-depth accounting closes.
- Tests: test_popoxy_conservation (numeric, all pass). Concerns: O V
  dominance is an ionization-state question (T/rates), not indexing.

## Task 3 - Number conservation - DONE

- Files: `generate_full_matrix_v3.f` (C_NC=sum of SE BA(N_SE,:), once
  per ion; C_ION replaces g.s. row), `fixpop_in_ba_v3.f` (EQ_SPECIES,
  SPECIES_END_ID, POP_SPECIES; historical IS bug already fixed).
- No double counting; advection via BA_ADV_TERM; XRAY_EQ channel
  present. Tests: none new. Concerns: none.

## Task 4 - Full-to-super mapping - DONE (defect found + fixed)

- Files: `sup_to_full_v3b.f`, `sup_to_full_v4.f`, `set_line_opac.f`,
  `rd_f_to_s_ids.f`. Reconstruction = C2_S * LTE-share, renormalized
  to exact super totals; LTE input reproduces LTE output.
- DEFECT: unguarded C2_S/SUM -> 0/0 NaN when a super level has zero
  population (latent, also upstream). FIXED: guard + conserving
  even-split fallback (sup_to_full_v3b.f). Same defect in uncalled
  sup_to_full_v3.f left as dead code.
- Tests: test_sup_to_full driver (A-E invariants, 7/7 pass).
  Concerns: INT_SEQ interpolation path covered by code reading, not
  numeric test.

## Task 5 - Transition indexing - DONE

- Files: `set_frequency_grid_v2.f`, `line_vec_mod.f`, `set_line_opac.f`.
- VEC_NL/NUP = same-ion POPS rows; VEC_MNL_F/MNUP_F full indices; all
  VEC_* arrays sorted on shared VEC_INDX. No cross-ion reference
  possible post-validation. Tests: none new. Concerns: none.

## Task 6 - Opacity/emissivity - DONE

- Files: `set_line_opac.f`, `sub_cmf_line.f`, `sub_sob_line_v3.f`.
- CHIL_MAT/ETAL_MAT built from POPS rows used by the rate equations;
  shared by CMF and non-mono Sobolev paths. Einstein form verified.
  Tests: none new. Concerns: none.

## Task 7 - Non-monotonic transfer path - DONE (code level)

- Both transfer paths consume identical CHIL_MAT/ETAL_MAT; Knop/Sob
  differences enter through J/transport, not population association.
  Knop-run POPOXY discrepancy tracked as open item (Task 13).
  Concerns: physical effect vs unconverged state unresolved.

## Task 8 - SE solver + rate matrix - DONE (audit level)

- `generate_full_matrix_v3.f`: each eqn assembled once; adjacent-stage
  coupling via C_ION rows NXzV+1..N_SE-1; species sums via C_NC.
- Concerns: charge-exchange/dielectronic term-by-term verification
  deferred; rate-data correctness (category D) remains the prime
  suspect if O V excess is physical.

## Task 9 - JFNK population vector - DONE

- `jfnk_cmfgen_interface_v1.f` uses same EQXzV/NXzV layout via
  POPTOION; committed vector = BA vector. O V cannot grow silently
  unless the matrix encodes it (rate question). Concerns: none.

## Task 10 - General population diagnostic - DONE

- New VADAT: POP_AUDIT (+ _SPEC/_ION/_D/_LEV filters). Per-iteration
  census at end of SOLVE_FOR_POPS: species conservation residual vs
  POP_SPECIES, ion totals, DXzV aliased-row staleness check,
  negative/non-finite censuses, optional per-level fractions and
  departure coefficients. Index-based filters, no species code.
- Files: control_variable_mod.f, rd_control_variables.f,
  solve_for_pops.f (+POP_AUDIT_V1). Compile: clean (ifort).
- Tests: test_pop_audit (static, all pass). Obsidian note: updated.

## Task 11 - Automated tests - DONE (subset)

- Added: test_sup_to_full (numeric driver), test_pop_audit (static),
  test_popoxy_conservation (model-output regression, env-configurable
  file path, skips cleanly).
- Not covered by dedicated tests yet: species-isolation perturbation,
  restart round-trip, JFNK-vs-BA vector equality (existing JFNK
  regression tests cover the solver side), non-mono multi-reversal
  line formation (covered by existing knop/sob tests).

## Task 12 - Audit doc + TODO - DONE

- docs/audits/SPECIES_POPULATION_AUDIT.md written; this log.

## Task 13 - Open verification items

- [x] Full `python unit_tests/run_all_tests.py` pass: **27/27 test
      files, 0 failures** (includes test_sup_to_full, test_pop_audit,
      test_popoxy_conservation). 07-Nov-2026.
- [ ] `build.sh` clean production build + `build_debug.sh`
      bounds-checked build.
- [ ] Run a model with `POP_AUDIT=T` + `POP_AUDIT_SPEC=<O index>` and
      attach the per-iteration output to the audit.
- [ ] Determine whether the knop_formal run's O III-dominated state is
      convergence or physics (compare ITER_HISTORY/thermal structure).
- [ ] Charge-exchange/dielectronic term-by-term rate audit if O V
      remains suspect after the above.
