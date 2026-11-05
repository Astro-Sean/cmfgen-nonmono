# OUTGEN iter-57 remediation task list

Ledger format: [ ] open / [x] done (with evidence line) / [!] blocked.
A task is done only when its stated evidence exists, never on compile alone.

## Investigation (no code changes)

- [x] T1 Recon: git state, failing OUTGEN, config.
      Evidence: HEAD=3de42aae=origin/main; run dir and config recorded in
      docs/audits/OUTGEN_ITER57_FAILURE_AUDIT.md section 0.
- [~] T2 Audit MOM_KNOP_REFRESH_V1 + per-frequency solver pipeline.
      Evidence so far: source map in audit doc section 3; residual,
      best-pass, AA, adaptive-relax, quarantine all traced.
      Remaining: XM_SAVE->J fallback at ~3219, H-rec vs assembled EPS signs.
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
- [ ] T6 Geometry/grid: confirm R/V arrays consistent between hydro reader,
      mom grid expansion (212), formal grid (VFRAC_FG=0.5); check
      IB_METH=DIFFUSION appropriateness at tau_inner=39.3.
- [ ] T7 Restart/checkpoint integrity of J_NEXT_FRQ, MOM_NU_GRID,
      MOM_QUAR_AGE, MOM_DAMP_NM (VER3 vs VER4 file formats).

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
- [ ] C5 Recovery policy: on unsafe-field declaration, define bounded
      recovery (hold + reduced relax vs skip) per the task's A-D cases.
- [x] C6 persistent quarantine now gates: MOM_NPQUAR >1% arms
      MOM_UNSAFE_B, >0.1% keeps MOM_FREEZE_B armed, in both gate
      branches (solve_for_pops.f ~1360, ~1416).
- [~] C7 submitted job 12719316 (warm restart at it57, fixed exe,
      MOM_XM_POS=T default): compare negJ/quar census, refresh
      residuals, MOM_XMNEG counts, and whether the residual runaway
      persists; original evidence kept as *.it57_orig.

## Numerical tests (unit_tests/, extend drv_mom_solve/drv_mom_sweep)

- [ ] A monotonic control: KNOP solver vs monotonic on same input.
- [ ] B single compression region sigma<-1.
- [ ] C multiple reversals.
- [ ] D ascending vs descending refresh same fixed point.
- [ ] E grid refinement invariance.
- [ ] F quarantine persistence with injected invalid cells.
- [ ] G moment/formal consistency localisation.
- [ ] H population coupling rejects unsafe field.
- [ ] I diffusion vs zero-flux boundaries.
- [ ] J restart/restart lagged-field reproducibility.

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
- [ ] V2 Debug exe + FP traps on synthetic drivers.

## Run-severity findings already verified (see audit doc)

- F1 negative assembled RHS at source-free compressing cells.
- F2 diagnostic counters/labels misreport (PASS_NOFF vacuous, FCMP stale,
      J_COMP label wrong, last-freq-only compare).
- F3 runaway cap permits 1E3/pass growth -> matches 620->1.6E6.
- F4 persistent quarantine excluded from drift metric while consumed as
      upwind data.
