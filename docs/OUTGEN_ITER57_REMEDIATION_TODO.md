# OUTGEN iter-57 remediation task list

Ledger format: [ ] open / [x] done (with evidence line) / [!] blocked.
A task is done only when its stated evidence exists, never on compile alone.

## Investigation (no code changes)

- [x] T1 Recon: git state, failing OUTGEN, config.
      Evidence: HEAD=3de42aae=origin/main; run dir and config recorded in
      docs/audits/OUTGEN_ITER57_FAILURE_AUDIT.md section 0.
- [x] T2 Audit MOM_KNOP_REFRESH_V1 + per-frequency solver pipeline.
      Evidence: source map in audit doc section 3; residual,
      best-pass, AA, adaptive-relax, quarantine all traced;
      XM_SAVE->J fallback verified correct (RHS restore, not J);
      runaway channel characterised (small exploding subset +
      AA extrapolation, ||f|| stays bounded).
- [x] T3 First negative-J occurrence identified and reproduced.
      Evidence: MOM_J_NONMONO_ERRORS shows XM_SAVE<0 before solve at
      compressing absorption-dominated cells with ETA~1e-133 (F1);
      driver mode 'seed' injects the same state and reproduces
      negative XM (-1.26E-9) with 41 negative solved J at band rows;
      TA/TB/TC/sigT/CHI_HAT in the per-cell dump match the
      production regime (absorp-dominated, |TA|+|TC|<<TB).
- [x] T4 solve_for_pops GS_PRECOND chain. Evidence: unsafe field
      from it1; demand trace CaI/FeI/Sk2 exempted trace levels at
      depths 1-5 and 64-65, fluctuating 3E15<->5E7 non-monotonically;
      n_precond counted each live-LAMBDA call (50 -> STOP). BA reads
      the same lagged field the formal solve stored -- the unstable
      input is the moment field, not a BA consistency defect.
- [ ] T5 J(mom)/J(ray) compare: determine if mismatch is boundary-only or
      widespread; produce per-frequency/depth worst-cell diagnostic data.
- [x] T6 Geometry/grid. Evidence: RDINR 1:1 mapping; RSTAR/RMAX
      consistent; mom grid T2=1.0*minVDOP vs formal T2=0.5*minVDOP
      (formal grid finer ~2x; both same interpolation mechanism);
      DIFFUSION at tau_inner=39.3 flagged by pipeline <50 warning,
      optically thick, not implicated.
- [x] T7 Restart/checkpoint integrity. Evidence: MOMJNXT_GS carried
      J/H_NEXT_FRQ + fingerprints + metrics but not MOM_QUAR_AGE;
      FIXED as VER5 (record added, restore accepts VER 3/4/5).
      MOM_NU_GRID/MOM_DAMP_NM re-derive per call, no checkpoint
      needed.

## Corrections (each: minimal diff, flag-gated unless pure bug)

- [x] C1 PASS_NOFF counter: now counts compressing rows with
      L+MOM_NM_J > KN_M_NCF (red-edge drop). Evidence: source edit;
      bound case reported per refresh pass.
- [x] C2 MOM_SW_FCMP stale zero: share now computed unconditionally
      when the census exists. Evidence: source edit at ~line 1670.
- [x] C3 Negative-RHS bound: VADAT key MOM_XM_POS (default T) caps the
      net lagged debit at the row's local source balance; counter
      MOM_XMNEG reported on sweep lagbnd line and per refresh pass.
      Evidence: driver mode 'seed' reproduced the production mechanism
      (XM_SAVE -1.26E-9 -> 0; NEGJ_T 41 -> 0); test_mom_solve 98/98,
      sweep_table 8/8, knop_ray 10/10, vadat_keys 77/77.
- [x] C4 J(mom)/J(ray) diagnostic widened: all-depth worst-cell
      census (JCMP_W_REL/ABS/NU/I, JCMP_N50) written to J_COMP and
      LUER at FREQ_INDX=NCF; ">1000%" label corrected to ">100%" to
      match the gate and the stock branch. Thresholds unchanged.
- [x] C5 Recovery policy assessed: existing tiers already split the
      A-D cases (UNSAFE/FREEZE/TRIAL/PRECOND_LIM); the field-side
      repair is supplied by C3+C6+VER5. Residual runaway traced to
      wandering GAM-floor resonance cells -- a conditioning
      limitation, not a gate defect. Documented in audit doc; no
      code change warranted without C7 evidence.
- [x] C6 persistent quarantine now gates: MOM_NPQUAR >1% arms
      MOM_UNSAFE_B, >0.1% keeps MOM_FREEZE_B armed, in both gate
      branches (solve_for_pops.f ~1360, ~1416).
- [x] C7 job 12719316 COMPLETE (warm restart at it57, pre-C8 exe):
      negJ 0 on every sweep (was 137K-797K); negative-RHS bound
      fires 16.7K-24K rows/sweep and falling; quar 452K->94K;
      pquar 198K->80K declining; refresh stalls ~1.3E4 (was 1.61E6
      detonation); ran it57->it63, 75% converged, graceful
      budget-exhausted end, no GS_PRECOND_LIM. The bound removed
      the first failure channel; the field stayed [BAD] via the
      separate non-contraction limitation. Run dir restart state
      was cleaned after the job; the pipeline relaunched the model
      cold-start as job 12725201 on the VER5+C8 exe (pending
      QOSGrpCpuLimit) -- the end-to-end C8 exercise.
- [~] T7 follow-on: VER5 checkpoint carries MOM_QUAR_AGE; fixed
      RESTORE_MOM_JNEXT_LATEST to accept VER5 + read the age
      record (was rejecting all VER5 rollback restores), and both
      readers now reset ages on the reject path. Driver ckpt mode
      asserts the age round-trip and the reject reset.
- [x] C8 residual gate (user choice: bulk+max hybrid):
      MOM_RF_DIV_BULK (VADAT, default T) requires the abort pass
      to regress on both metrics -- max-norm beyond 1E3*MIN_RES
      and N50 above the same direction's best. A lone
      resonance-floor argmax is logged as a tolerated excursion
      (MOM_RF_NDIVSKP); with no per-direction best yet the
      max-norm test stands alone; capped/non-finite writes still
      trip a true detonation through N50. At it58 pass 4
      N50=8198 was already under the bulk-exit bound, so a gated
      refresh can now reach the bulk exit instead of aborting.
      Interactive job 12723631: 100/100 + FP-trap; exe relinked.
      The running C7 job predates this build -- a requeue would
      exercise it on the real model.
- [~] Census growth chased: the mom-vs-stored >50% count peaked
      at it59 (1.1M) then retreated at it60 (999K); quar 143K,
      pquar 125K approaching the 119.5K release arm. Demands are
      ordinary 1E2-1E3% trace churn at varying depths. Read as
      hold-release catch-up, not a new instability -- watch only.

## Numerical tests (unit_tests/, extend drv_mom_solve/drv_mom_sweep)

- [x] A monotonic control: `v11` mode drives stock MOM_J_CMF_V11 on
      the homo atmosphere alongside the non-mono solve; the two agree
      to V11_MAXREL=0.059 (grid-expansion regrid residual). Job
      12729502, mom_solve 130/130.
- [x] B single compression region: `comp1` (sigma<-1 band, coefficient
      and convergence checks) + `zedg`/`zsrc`/`seed` variants.
- [x] C multiple reversals: `comp2` (two disjoint sigma<-1 bands)
      now exercised -- the test asserts two separate dV/dR sign
      reversals on the built profile plus coefficient/admissibility
      checks (reversal runs=2, 120 compressing rows). Job 12729448.
- [x] D ascending vs descending: `rfd` mode calls the production
      MOM_KNOP_REFRESH_V1 (argument-free, drives from KN_M_* capture)
      under forced DESC and ASC orders from the same lagged store;
      committed fields agree to RF_ASCDESC_MAXREL=0.035. Job
      12729477.
- [x] E grid refinement invariance: ND/NCF are -D defines on
      drv_mom_solve (default 24/8, -fpp/-cpp added to common.py);
      the seed band is defined in fractional depth. ND=48 vs ND=24
      gives the same mechanism and bound behaviour (SEED_XMNEG
      23 -> 31, XMSAVE_MIN=0, NEGJ_T=0, field positive) --
      job 12729428, mom_solve 106/106.
- [~] F quarantine persistence: seed mode injects invalid cells and
      asserts repair counters; NPQUAR gate coverage is at
      source-assert level in test_knop_lagged_field.
- [~] G moment/formal consistency localisation: JCMP census added
      (C4); unit-level injection test not yet written.
- [x] H population coupling rejects unsafe field: gate logic pinned
      by test_iter_conv_gate + test_knop_lagged_field; NPQUAR arms
      MOM_UNSAFE_B/MOM_FREEZE_B (C6).
- [x] I boundaries: comp1 runs under both inbm settings; boundary
      coefficient checks in the discrete-operator block.
- [~] J restart: ckpt mode covers fingerprint/tag rejection and
      J/H_NEXT_FRQ restore; VER5 adds MOM_QUAR_AGE to the record
      (tested by the same ckpt path).

## Validation

- [x] V1 Build production exe clean: cmfgen_dev.exe linked, zero
      errors (ifort, all changed objects archived into libsubs.a).
- [x] V3 Unit suite passes post-change: mom_solve 98/98 (incl. new
      seed regression), sweep_table 8/8, knop_ray 10/10,
      vadat_keys 77/77, gs_drift_diag, knop_lagged_field,
      knop_status_text, iter_conv_gate 13/13 all pass.
- [x] V4 it57 signature reproduced in controlled driver (seed mode:
      negative XM before fix, eliminated after) -- see audit doc
      section 1b. NOT yet demonstrated on the real model (C7).
- [x] V5 Full run_all_tests.py sweep: 30/30 test files pass
      (post C4/C6, incl. all gate and VADAT source-pin tests).
      Re-verified on the interactive partition (job 12720091):
      mom_solve 98/98, sweep_table 8/8, vadat_keys 77/77,
      iter_conv_gate 13/13, knop_ray 10/10.
- [~] V2 FP traps on synthetic drivers done: drv_mom_solve rebuilt
      with -check all -fpe0 on the interactive partition; seed case
      runs clean, bound fires (SEED_XMNEG=41, NEGJ_T=0), exposed and
      fixed a format item-count bug in the _sm.txt dump. Full
      build_debug.sh exe run still open.

## Run-severity findings already verified (see audit doc)

- F1 negative assembled RHS at source-free compressing cells.
- F2 diagnostic counters/labels misreport (PASS_NOFF vacuous, FCMP stale,
      J_COMP label wrong, last-freq-only compare).
- F3 runaway cap permits 1E3/pass growth -> matches 620->1.6E6.
- F4 persistent quarantine excluded from drift metric while consumed as
      upwind data.
