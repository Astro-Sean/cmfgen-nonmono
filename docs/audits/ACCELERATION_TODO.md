# ACCELERATION TODO - DO_NONMONO_SOB

Each item lists source files/routines and acceptance criteria.
Status: [ ] open, [~] in progress, [x] done, [-] dropped/blocked.
Audit basis: ACCELERATION_AUDIT.md (rev 280c4a44, 10-Nov-2026).

Naming: D1 = initial per-iteration accelerator ledger;
D1a = ledger state-accuracy corrections (the "D1b"/"D2" label used
in earlier drafts and the AUDIT_LOG before this reconciliation);
D2 = the final OUTGEN aggregate accelerator summary (still open).

Status (10-Nov-2026): ledger v3 (aa_ran, aa_en, aao_en) implemented;
GET_MOM_AA_STATE state/activity separation implemented;
stale-counter suppression via MOM_RF_RAN implemented; outer AA
independence from MOM_RF_AA implemented; final OUTGEN aggregate
summary NOT IMPLEMENTED (D2 open).

## Diagnostics (no solver change)

- [x] D1 - Per-iteration accelerator ledger columns (done 10-Nov-2026)
  Files: new_main/mod_subs/iter_ledger_mod.f90,
         new_main/mod_subs/solve_for_pops.f,
         plane/mom_j_cmf_nonmono_v1.f (GET_MOM_AA_STATE + MOM_AA_RST),
         unit_tests/f90/drv_iter_ledger.f90,
         unit_tests/python/test_iter_ledger.py
  Delivered: v2 columns appended after rss_mb; ITER_SUM 'accel' row
  prints the same census (ng:-|d|r|A, jfnk:-|a|A, aa counts or
  aa=off, AOFAIL, freqm); MOM_AA_RST counts damp-triggered inner-AA
  history restarts; test_iter_ledger 31/31 green; stock path still
  writes nothing (LEDG_ACTIVE guard).
  Append to the ITER_LEDGER row (column contract v2, at end):
  ng_try,ng_ok,ng_def,ng_cln,jfnk_cmt,
  aa_app,aa_skip,aa_sng,aa_rst,
  aao_app,aao_rej,aao_sng,aao_fail,freqm
  Distinguishes: inner/outer AA applied vs skipped, AA history
  alloc failure, AA history restart count, NG attempted/deferred/
  fired, NG_CLEAN window state, JFNK admitted vs committed,
  FREQ_MATCH armed. t_rad/t_rf/t_pops already expose the dominant
  bottleneck by magnitude.
  Acceptance: header v2 written only on new files; restarted runs
  appending to a v1 file keep old column names at old positions;
  drv_iter_ledger + test_iter_ledger updated and green; zero output
  on the stock path (LEDG_ACTIVE guard unchanged).

- [x] D1a - Ledger state accuracy corrections (done 10-Nov-2026)
  Files: plane/mom_j_cmf_nonmono_v1.f (GET_MOM_AA_STATE + RF_RAN
         gate, outer gate decoupled to MOM_AA_OUT alone),
         new_main/mod_subs/iter_ledger_mod.f90 (v3: aa_ran,aa_en,
         aao_en appended), new_main/mod_subs/solve_for_pops.f
         (accel row states vs activity), rd_control_variables.f
         (MOM_AA_OUT help), unit_tests/f90/drv_mom_aa_state.f90,
         unit_tests/python/test_mom_aa_state.py (new, 62/62),
         drv_iter_ledger.f90, test_iter_ledger.py (v3)
  Delivered: stale-counter fix (get counters gated on MOM_RF_RAN);
  enabled vs applied separated; alloc-fail latch distinguishable;
  outer AA now independent of MOM_RF_AA (defaults unchanged);
  test_iter_ledger 31/31, test_knop_ray 41/41.

- [ ] D2 - OUTGEN accelerator summary at run end
  File: new_main/cmfgen_sub.f (final accounting block)
  Print one compact block: NG attempted/fired/deferred counts,
  JFNK admitted/committed/skipped counts, AA app/skip/sng totals,
  AA alloc-fail flag, FREQ_MATCH state, wall split
  (formal solve vs population vs refresh).
  Acceptance: block prints under DO_NONMONO_SOB/KNOP only; values
  match ITER_LEDGER sums; no output on stock path.

## Experiments (Slurm, one change per arm, same model + tolerances)

- [ ] E1 - Warm baseline: np recipe, warm from adp2_tauto it1445.
  Acceptance: converges; record its-to-converge and wall.
- [ ] E2 - E1 + DO_ES_ALI=T (ES_ALI_TOL=1e-4, ES_ALI_MAX=20).
  Acceptance: iteration count vs E1; RJ loop reaches ES_ALI_TOL in
  ~2 passes (per 27-Oct verification); solution unchanged.
- [ ] E3 - E2 + DO_ES_ALI_BA=T (ES_ALI_BA_MIN=0.1).
  Acceptance: breaks the ~0.97/it fixed-ratio decay; population
  demand decays measurably faster; no extra rejections.
- [ ] E4 - E3 + DO_BA_RES_ACC=T (the absent-key auto-set value).
  Acceptance: trial releases fire only on marginal-tier fields;
  NG_CLEAN bookkeeping stays consistent.
- [ ] E5 - E1 + JFNK_COUPLED=T + JFNK_COUPLED_SOB=T +
  JFNK_CPL_SOB_MAX_EVALS=20.
  Acceptance: A3 consistency gate |F_eval-F_prod|_w <= 1e-5 on
  first coupled admission (else the tier self-detaches and the arm
  degenerates to E1); count evals/JVPs/committed steps vs E1;
  compare WALL time not just iterations.
- [ ] E6 - LAM_ES=F characterization arm (not a speed proposal):
  FREQ_MATCH=T, AA census columns from D1 read live. Records the
  limit-cycle structure for the mom-path revival decision.
- [ ] E7 - DO_JFNK=F on np recipe (only if E5 fails to commit):
  removes ~322 inert admissions (~3-4.5 h overhead).
  Acceptance: identical iteration count to E1 minus JFNK overhead.

## Phase-8 validation (applies to the winning arm)

- [ ] V1 - Final populations match baseline within EPS_TERM envelope
       (row-by-row SCRTEMP comparison).
- [ ] V2 - Species/ion conservation unchanged (POP_AUDIT or equivalent).
- [ ] V3 - Emergent spectrum unchanged within numerical tolerance
       (obs_frame comparison vs baseline .fl).
- [ ] V4 - No new non-finite/capped/quarantined field values
       (ITER_HISTORY mom_* and fg_* columns flat at zero).
- [ ] V5 - 1 vs 8/16 OMP thread reproducibility on a short segment.

## Deferred / not recommended now

- [-] X1 - A third AA scheme. Rejected by audit: the existing
  two-level AA is correctly implemented and correctly refusing a
  non-contractive map; the limit cycle is a discretization defect.
- [-] X2 - Global MOM_RF_GS_DESC enable. Frozen-read propagation
  limit is a mom-path question only; production path unaffected.
  Revisit inside E6 if the limit cycle is ever tamed.
- [-] X3 - Raising relaxation factors / weakening maturity or
  convergence gates. Explicitly out of bounds per audit rules.
