# CHECKLIST - cmfgen_jfnk development ledger

- [~] BLKBAND zero-row audit - PARTIAL verdict (24-Nov-2026):
      Ni10 day-40 saved BA_ASCI_N_D5 is 193x193 with no exact zero
      rows/columns or zero populations, but 12 non-finite entries:
      rows 173-176,179-180 in columns 180 (FeIII) and 193 (T).
      RHS is finite; POPS(180)=3.9648488283438e-163. D41 and DND
      captures are finite with no exact zero rows/columns. D5 SHA256:
      a1245865ccaabb1ad3f297fad09e87cb0e81f8e44dd9384bd6596980b36ab684.
      Capture is after GENERATE_FULL_MATRIX_V3 population-column scaling,
      before null pinning, DGEEQU, or DGETRF. Invalid derivatives therefore
      predate factorization; small finite populations are not proof of
      null equations. Source: DGEEQU INFO ignored; failed DIAG depths
      zero STEQ and still return FLAG=T. Consequently the output finite
      check and BA-SOLVE ok=T do not certify a valid whole-grid solve.
      Snapshot OUTGEN has 562 failures, iterations 23-106, depths
      1-6,11-12, and no null-pinning notices. The run was live; counts
      and captures are not synchronized to a single failing call.
      Frozen executable dac3e784948a2a4e has the null-pin symbol;
      manifest records git 4f69b20b with 19 dirty files, so current
      source is not proven identical to that binary.
      Production-object tests: test_cmf_blkband_pin.py 6 checks PASS;
      test_sub_sob_line.py 22 checks PASS. They do not cover tiny-positive
      emissivity overflow or this captured model. Diagnostic SciPy LAPACK
      replay of D5 returned INFO=0 but 193 non-finite solution entries;
      this is not production-MKL validation.
      FIX landed (this commit): reproduced through the production
      sub_sob_line_v3.o - the lambda-iteration derivative
      VC_SIM=JBAR*CHIL_MAT/ETAL_MAT**2 (and SRAT=ETAL/ETAL_MAT) overflows
      to Inf when ETAL_MAT is positive but pinned near the population
      floor; tiny POPS cannot rescue an already-infinite product, so
      Inf/NaN enters the BA Jacobian (the captured D5 column-180/193
      pattern). The existing ETAL_MAT sign guard now also zeroes
      non-finite VB_SIM/VC_SIM, and SRAT is zeroed when its quotient is
      non-finite - same defect class, same zero-the-unusable-derivative
      remedy already in the routine. New drv case tiny_etal
      (ETAL_MAT=1E-160): before the fix VB_SIM printed -9.3E159 and
      NNAN>0; after, VB/VC zeroed and all BA finite. test_sub_sob_line
      25/25 PASS incl. sig_clmp=base bitwise (default path unchanged).
      Remaining: full-model verification (a rerun of the Ni10 40d arm on
      this binary) and the separate diagnostics gap - DGEEQU IFAIL is
      unchecked and a failed DIAG depth returns ok=T with zeroed STEQ -
      left open, needs a user decision on the solve contract.

- [~] Declared bounded investigation K1 phase 2 - period-2
      damping arm plus NI=1/handoff gate fixes (24-Nov-2026,
      work-order on branch knop/lam-handoff-ni1). Audit findings
      that motivate it: (i) the 9-14 STAG rays per refresh are the
      wholly compressing core rays - the worst ray is a ray 1-12
      core ray in 109/122 iterations, and the depth-node estimate
      of GAM is negative at every point on rays 1-11; phase 2.28
      showed ray 81's stall is an exact period-2 border-collision
      limit cycle and phase 2.29A hardened the KN_ORD_P2DAMP state
      machine for it, while phase 2.29 rejects block-GS for this
      class. (ii) The NI=1 ray (LS=85) returns before the
      classifier, stays KN_CS_INV, and made FG_RF_ALL_CONV
      unreachable: the model sat in LAMBDA all 122 iterations
      (rf_not_allconv on every ITER_MODE line, demand already
      below LAM_EXIT, stable counter 0/5). (iii) The rf_not_allconv
      fold-in leaked an exit-grade criterion into the handoff
      gate, deadlocking LAMBDA.
      Changes (branch knop/lam-handoff-ni1, flag-gated, code at
      7abcc264..9394f6d6): KN_NI1_FULL certifies the NI=1 ray
      vacuous FULL
      (FG_NI1_FULL -> KN_CS_FULL/reason 0 in the NI=1 early
      return); KN_LAM_PARTIAL_HANDOFF splits the field verdicts -
      LAM_FIELD_NOT_CONVERGED (strict, exit + CONV_EXIT_EVAL ninth
      conjunct) always carries the rf_not_allconv fold, while
      LAM_FIELD_NOT_READY (weak, LAM_STABLE_STEP / LAM_HANDOFF_OK /
      stable counter) carries it only when the flag is off.
      Arm T1 (phase 2.28 C1, zero code):
      configs/dessart_m15_2p82d_FRESH_knop_v14_p2d.yaml - warm
      restart from the failed it122 checkpoint (same seed as the
      nbulk arm for direct attribution), baseline-pinned exe
      (same sha as nbulk), KN_ORD_P2DAMP=T + P2N=2/P2REL=3 +
      KN_RAY_BULK_TOL=0, KN_RDUMP_LS=11/12, KN_Q_RESID=T; CMFGEN
      13639702, OBS 13639703. VERDICT NEGATIVE (24-Nov-2026,
      canceled after 2 iterations): p2damped fired pervasively
      (1.0e9 / 7.8e8 writes per refresh; 23M/17M on worst ray 11)
      but the floor held at 4.31/3.95e-2 (baseline band
      2.5-4.4e-2) and cert stag went UP 9-14 -> 28/22, not down
      to <=2 - the sign-flip tracker catches noise flips on
      genuinely contracting cells and midpoint-commits land off
      the fixed point. The orbit is not a pure multiplicative
      lambda=-1 mode (fallback arm KN_ORD_AITKEN=T remains open,
      gated on this report).
      Arm T3-validation: same restart + KN_NI1_FULL=T,
      KN_LAM_PARTIAL_HANDOFF=T (no P2DAMP - T1 negative) -
      configs/dessart_m15_2p82d_FRESH_knop_v14_hand.yaml; CMFGEN
      13640311, OBS 13640312, branch pin
      bfffac62006c60d3_0e7b904d4a648954 at a75e3c7c. VERDICT:
      handoff mechanism works, first full step detonates
      (canceled after it130). (a) inv=0 from it122 - the NI=1
      ray certifies FULL under KN_NI1_FULL. (b) CNT_STABLE_LAM
      advanced 0->5 over it125-129 on the rf_not_allconv-only
      calls (flag-off ordf control: pinned 0/5 through 6+ free
      iterations); the conv-exit FAILS label correctly splits
      'lam-field-not-ready' (weak+strict hold) from
      'lam-field-not-converged' (strict-only hold). (c)
      DO_LAM_AUTO switched to full iterations at it129 - the
      LAMBDA deadlock broke, field_ready=T, demand 15.9%.
      (d) it130 demanded 1.0e7 raw, applied a bounded 400%,
      reverted FULL->LAMBDA: the SE residual sits at depths
      66-70 (41.6% of the column at depth 70, worst eq
      NeV:2s2_2p2_3Pe), exempted FeI demand 385% at depth 16.
      The strict gate's all-FULL leg exists precisely because
      an unconverged (4e-2-floor ordered march) lagged field
      poisons the first Newton step; the weak gate alone is
      not sufficient evidence on this model. Production keeps
      the strict behavior (flags default F).
      Arm T4-B (T1 failed -> triggered): ordered A/B, arm B is
      KN_RAY_ORDERED=F - configs/dessart_m15_2p82d_FRESH_knop_v14
      _ordf.yaml; CMFGEN 13640315, OBS 13640316, baseline pin
      7f4c690ab4cac684_6a33ed6eb6ecf72b (same sha as arm A =
      nbulk run 13636956, not re-run). RDUMP rays 11 (worst,
      ALLNEG), 12 (ALLNEG), 22 (mixed-sign non-committer; GAM
      census from SN_HYDRO_DATA: rays 1-15 ALLNEG, 16-36 MIXED,
      37-85 ALLPOS). VERDICT: unconditional win over 13
      refreshes (canceled) - every production ray certifies
      FULL (census full=84, stag=0, unconverged cells=0),
      residual 1.8e-3 -> 5.1e-4 (< tol 2e-3), independent q
      4.2e-2 -> 6.9e-4; only the NI=1 ray still holds
      rf_all_conv off. Deck GAM profiles confirm the census:
      ray 11 = 140/140 negative, ray 12 = 127/129 (two
      marginal blips), ray 22 = 102/128 with three sign
      crossings - wholly compressing and mixed rays alike
      converge under the frequency-major march, so no
      mixed-only restriction is motivated; the ordered march
      itself is the stall source. Next-step candidate:
      KN_RAY_ORDERED=F + KN_NI1_FULL=T (strict gate opens on
      honest all-FULL evidence; KN_LAM_PARTIAL_HANDOFF not
      needed). Aitken fallback unnecessary - orbit already
      gone. Read GAM profile from the .bin deck before
      any mixed-only
      dispatch proposal (GAM is not SIGMA: T2 flips at
      MU**2=-1/SIGMA for SIGMA<-1).
- [x] KNOP V14 cold-start A/B arms + bulk-exit arm: VERDICT
      (24-Nov-2026, serves K1). Both cold-start arms FAILED at
      it122 (13606657 DO_JFNK=F, 13606660 DO_JFNK=T, exit 6):
      SOLVE_FOR_POPS stopped on GS_REL_NOPROG_LIM (10 releases
      with no demand progress; final demand 23.3 vs release
      demand 18.2). JFNK never ran: 542 gate evaluations, all
      JFNK DEFER (demand 16.8-23.3 vs freeze gate 3.0), so the
      two arms followed the same LAMBDA path with identical
      FG-KNOP statistics. V14 does not cure the KNOP cold start.
      Refresh residual floored at 3.7-5.2e-2 (tol 2e-3) with
      fg_sw_one sawtoothing on the release/refreeze cadence and
      ending at 0.26.
      Follow-on arm (hypothesis: the floor is KN_RAY_BULK_TOL
      truncation at NPO>=6, which also bypasses TRY_HOT):
      configs/dessart_m15_2p82d_FRESH_knop_v14_nbulk.yaml, warm
      restart from the failed DO_JFNK=F run (it122 checkpoint),
      identical except KN_RAY_BULK_TOL=0.0; CMFGEN 13636956 +
      OBS 13636957, exe_frozen/7f4c690ab4cac684 (git=e08c16b8,
      dirty x8 test/audit paths). REJECTED over it122-128:
      bulk=0 every refresh, residual 2.5-4.4e-2 with no downward
      trend, tol rays 75-79 (was 70-74) oscillating, the same
      9-14 rays exit stall/plat and commit capped, unconverged
      cells ~3.7e3, rho_p90=0.89 rho_max=0.98, fg_sw_one re-pinned
      at 0.25, hot-pass spend 10-27 of 480. First release (it125,
      raw_drift 0.26) again grew the demand (16.7 -> 19.3,
      release_no_progress count=1). Both jobs cancelled once the
      census was unambiguous. The floor is bounded non-contraction
      of the lagged map on the compressing hot-frequency runs
      (852-1654 frequencies in 202-465 contiguous runs, max
      donor offset 371 bins), not an exit-policy artifact.
      Manufactured check: unit_tests/f90/drv_knop_block.f90
      (untracked, no consuming test yet) at chi~1e-5 and
      |GAM|*dZ/DLN ~ 10 bins converges in 4-5 passes and 2-4
      refreshes on both sweep paths, so a frozen-coefficient
      single ray does not reproduce the stall.
      Side finding: the frequency-form 3*GAM*I term is absent
      from the solved ray equation (inherited from stock V13);
      audit estimate ~0.025% of the implicit coupling, not a
      drift cause, no source change made.
      K1 stays open. Remaining candidates are BLOCK_GS_DESIGN
      Stage 2+3 (per-block solves and verdicts on the hot runs)
      and Stage 4 (per-block Anderson); the release controller
      overshoot is a separate defect. A further phase needs a
      user check-in (section 3.3).

- [x] KNOP V14 cold-start A/B arms SUBMITTED (23-Nov-2026): two
      fresh LTE-start (LTE_EST=T, no warm restart) KNOP formal
      (USE_LAM_ES=T, DIAG) runs on VEL_LAW=14, identical except
      DO_JFNK. Configs generated from the FRESH_knop recipe with
      vel_law 11 -> 14 in both cmfgen.vel_law and VADAT VEL_LAW:
      configs/dessart_m15_2p82d_FRESH_knop_v14.yaml (DO_JFNK=F,
      FINAL 13606657 + OBS 13606658) and
      dessart_m15_2p82d_FRESH_knop_jfnk_v14.yaml (DO_JFNK=T,
      production JFNK gates, FINAL 13606660 + OBS 13606661).
      Run dirs: cmfgen_validation/m15_2p82d_FRESH_knop{,_jfnk}_v14.
      Question (serves K1): does the V14 tau-Rosseland grid cure the
      V13 CDS-cell-jump dead-row stall that held every prior KNOP
      cold start, and does JFNK admission change the outcome? Exit
      criterion: converged write (Model Convergence) or a named
      failure mode; JFNK admission or sustained non-admission in
      OUTGEN. SN_HYDRO_DATA validated via prepare_VADAT.parse_sn_hydro
      (ND=70, R increasing, V/T/rho positive, finite, non-monotonic
      reversal present). Binary: exe_frozen/7f4c690ab4cac684_6a33ed6eb6ecf72b
      (git=dd136b82, dirty x6 doc/fixture paths).

- [~] JFNK integration audit, 23-Nov-2026: PASS WITH LIMITATIONS.
      Scope: SOB + USE_LAM_ES=T + DIAG; input, admission, residual,
      update/rollback, TRI selection, and unsupported-route guards.
      Source trace: rd_control_variables parses DO_JFNK into
      USE_JFNK_SOLVER; solve_for_pops:2827-2872 dispatches either
      JFNK_SOLVE_CORR_V1 or conventional SOLVEBA_V13. Admission requires
      warmup, an eligible non-Lambda call, contracting last full-call
      demand below the gate, and hysteresis (SOB stable-Lambda relief
      can reduce NREQ to 1). DO_JFNK does not imply Newton executed.
      Frozen residual includes stored radiative terms and CMAT_RAD
      linearization; it is not a fully coupled transfer Newton solve.
      F(0) mismatch above 1e-6 rejects frozen Newton. Applied frozen
      steps are re-evaluated against the entry residual; non-improving
      steps restore populations and fall back to SOLVEBA. Failed ran
      attempts force a BA rebuild. Fixed-T masking is executable-tested.
      Coupled sweep registration in cmfgen_sub:1463-1466 requires
      JFNK_COUPLED. SOB additionally needs JFNK_COUPLED_SOB and uses the
      reach-fail fallback, not demand admission. Neither flag is set in
      the current arm: it tests frozen-J JFNK, not coupled acceleration.
      TRI: prepare_VADAT forces DIAG unless JFNK_USE_BANDED=T, when it
      selects TRI; this audit does not validate a live banded model.
      Tests on interactive jobs 13602310 and 13602407: solver suite
      89 reported pass / 0 unexpected fail, operator battery 36/0,
      fixed-T predicate 1/0 (four modes), VADAT parser 85/0. Source-only
      test_jfnk_regressions.py also passed. The 89 includes four known
      XFAIL diagnostics and skipped inapplicable cases; three XFAILs
      are poor per-depth LU residual accuracy on ill-conditioned
      surrogate blocks (mild 7.26e8, hard 8.54e14, shock 8.29e14).
      They are numerical-quality limitations, not new regressions.
      Live job 13601645: through it16, LAMBDA_TFIX; OUTGEN explicitly
      defers JFNK on Lambda-designated iterations, with no admission.
      Effectiveness and speedup remain UNVERIFIED on this cold arm.
      Its NG=T, NT=F, and newer binary differ from the conventional
      production arm, so this pair is not a controlled timing A/B.
      Status-report correction: ITER_HISTORY se/se2 columns are SE
      diagnostics, not outer/inner luminosities; lum at it16 is 7.902.
      No solver source, run settings, or live jobs changed. Verdict:
      wiring supported by inspected code and targeted tests; full-model
      residual repeatability, accepted Newton progress, coupled/banded
      production correctness, and wall-clock benefit are not certified.
      This audit closes with DEFERRED effectiveness pending actual
      admissions in the existing run; no follow-on campaign submitted.

## Per-iteration OUTGEN layout refactor (23-Nov-2026)

- [x] Formatting-only refactor of the per-iteration console report.
      New free-form module `new_main/mod_subs/iter_log_fmt_mod.f90`
      (80-col rules/section banners, aligned `key : value` rows,
      end-of-iteration `IT N | mode -> mode | BA x | Conv y% | t s`
      banner). Consumers: `cmfgen_sub.f`, `comp_j_blank.f`,
      `solve_for_pops.f`, `solveba_v13.f`; chatter gates added in
      `se_ba_non_therm_v3.f90` (species dump -> level 2),
      `electron_non_therm_spec.f90`, `read_ba_data_v3.f`,
      `fixpop_in_ba_v3.f`, `fiddle_pop_corrections_v2.f`. RUNAWAY
      census is now one top-N table + suppressed-count line;
      DEC/INC vectors carry an index header; JFNK DEFER/ADMIT
      verdicts and BA progress lines gate at >=1. No new VADAT
      key: semantics ride on existing `CONV_DIAG_LVL` (0-2), see
      NONMONO_KEYS section 6. Parser tokens kept verbatim
      (test_iter_output.py pins).
- [x] ./build.sh clean; all touched objects recompiled;
      cmfgen_dev.exe/cmf_flux.exe/obs_frame.exe produced.
      test_iter_output.py 51/51, test_alteration_coverage.py 3/3.
      Pre-existing suite failures unchanged and unrelated:
      test_auto_ion (missing fixture dir), test_es_conv
      (drv_knop_guard misses fg_j_cmf_nonmono_v1.o link dep,
      `set_fg_q_resid_` - predates this change), knop_lagged_map
      (labeled KNOWN DEFECT), rv_snhyd 4 pinned defects,
      par_opac + phase227_io environment issues.
- [x] Live A/B: warm-restart copies of converged
      m15_2p82d_VAL_sob_diag_ngoff (it409 state, 8 OMP threads,
      interactive partition, jobs 13599622/13599623). Old-binary
      OUTGEN vs new-binary OUTGEN across it410-412: identical
      Model Convergence (89%), luminosity, ES_ALI/JCOMP counts,
      Max % inc/dec, ITER_SUM pops/residual rows; `Maximm
      changes` (MAXCH) differs at ~1e-10 relative - OpenMP
      reduction ordering noise, present between any two runs;
      zero arithmetic statements changed in the diff. NOT
      bit-identical in the strict sense; numerically equivalent
      to displayed precision. New layout verified live: section
      headers, kv rows, DEC/INC index header, IT banner
      (`IT 410 | FULL_TFIX -> FULL_TFIX | BA HELD | Conv 89% |
      376.5 s`). RUNAWAY table path not exercised live (model
      below POP_RUNAWAY_DEM); levels 0/1 covered by source pins
      only.

## Standing directive: vel_law=13 (22-Nov-2026)

User directive: use vel_law=13 for all future runs (SOB and
KNOP). Recorded in AGENTS.md section 1.2. New configs carry
`cmfgen.vel_law: 13` + explicit `vadat_overrides VEL_LAW='13'`
(the forced non-mono table defaults VEL_LAW=11 when the key is
absent). Rationale from the ND=70 production-grid A/B/C
(below): the tau-ES grid puts ~80% of depths in the CSM where
narrow cores form. VEL_LAW=14 remains available as an explicit
per-arm override (e.g. the thinT continuation 13570980 still
runs V14 - it is a V14-state warm restart, not a fresh
default-path run).

## VEL_LAW=14 robustness audit (22-Nov-2026, user-directed)

- [x] VERIFIED-with-pinned-defects: construction matches the
      Python replica node-for-node on synthetic (5e-12) and the
      real m15 file (2.3e-14); production RVTJ grid = driver
      output to 2e-10. The blend+insertion fixes the V12 CSM
      starvation (36 vs 61 CSM depths, 19 vs 7 ejecta on the
      real file); TAU_LO cap is load-bearing (literal-decades
      reading would drift back to 49/13). Feature insertion
      lands on the CDS sigma-min to the ulp; kappa steering
      verified above the cap; ltau-shift invariance exact
      (6.6e-14). LOG_DEC_CAP comment corrected (natural-log
      units, not decades). Four defects pinned failing in
      test_rv_snhyd.py (42 pass / 4 pinned): D1 non-finite
      Kappa silently absorbed by the floor MAX -> shifted grid;
      D2 non-finite Density/V/R -> LIN_INTERP bracket-walk OOB
      segfault rc=174; D3 negative interior Density silent;
      D4 missing Sigma silently kills feature detection +
      assumes homology. Fix candidate (input finiteness/
      positivity screen + missing-Sigma warning) deferred for
      sign-off. Full audit:
      docs/audits/VEL_LAW14_SNHYD_GRID_AUDIT.md.

## Cold-start CONVENTIONAL-solver campaign (22-Nov-2026, refocused)

User directive (22-Nov-2026): both DO_NONMONO_SOB and
DO_NONMONO_KNOP must produce a converged model from a cold start
(LTE_EST=T) under the CONVENTIONAL population solver - DO_JFNK=F.
JFNK is unproven and does not count toward this goal. DIAG
satisfies the matrix-option requirement for now (TRI is a
follow-up cell). Both routes must then emit a line-blanketed
observer-frame spectrum for quantitative comparison against
Luc's m15 .fl (tooling: cmfgen_validation/compare_finished_vs_luc.py).

- [x] Route audit: both routes already cold-start (LTE_EST=T,
      warm_restart=false in the resolved configs). Neither
      converges under the conventional solver: sob_diag
      (DO_JFNK=F, DIAG) exhausted NUM_ITS at it387, 88%;
      knop_formal_rd stalled at 44% (killed at it250 by the
      stall detector, lagged field never matured).
- [x] Attribution: sob_diag and the it718-converged
      m15_2p82d_SOB_l14_lte_np are the SAME model with identical
      convergence aids (BA_DEM_PFLOOR=1e-30, BA_TR_MAXSWP=3,
      BA_EX_RELFRAC=0, ESEC_CHK, EPS_TERM=5e-2). The single
      substantive difference is DO_JFNK=T vs F - JFNK converged
      where DIAG oscillates. Under the refocused goal the np
      convergence is a REFERENCE fixed point, not the deliverable.
- [~] np verification arm RUNNING (13564952 + OBS 13564953,
      submitted 18:22 22-Nov): np recipe verbatim
      (DO_JFNK=T + USE_LAM_ES=T + obs_line_treatment=blank) on
      the post-RESC_JMAX build. Reclassified role: checks that
      the recipe still reproduces it718 convergence on the
      current build and supplies the converged reference state
      + BLANK spectrum for the conventional-solver work and the
      Luc overlay. Does NOT satisfy the campaign goal.
- [~] JFNK production-level arm RUNNING (13573913 + OBS
      13573914, submitted 22-Nov, exe hash 5c27abfe):
      FRESH_sob recipe + DO_JFNK=T under the 22-Nov directive
      VEL_LAW=13, cold LTE start, pinned T, DIAG, ND=70
      (dessart_m15_2p82d_FRESH_sob_jfnk.yaml ->
      cmfgen_validation/m15_2p82d_FRESH_sob_jfnk). Question:
      does the admission-gated JFNK pathway reproduce
      it718-class convergence on the current build under the
      mandated velocity law? Exit criterion: JFNK-FRZ/JFNK
      admission visible in OUTGEN plus Model Convergence, or a
      named failure mode. Distinct from the np arm (l14) by
      velocity law only.
- [~] KNOP JFNK production-level arm RUNNING (13574052 + OBS
      13574053, submitted 22-Nov, exe hash 5c27abfe):
      FRESH_knop recipe + DO_JFNK=T under the 22-Nov directive
      VEL_LAW=13, cold LTE start, pinned T, DIAG, ND=70, 64G
      (dessart_m15_2p82d_FRESH_knop_jfnk.yaml ->
      cmfgen_validation/m15_2p82d_FRESH_knop_jfnk). Production
      JFNK gates (FRZ 3.0, COR 5.0, NREQ 2, WARMUP 6) and
      production GS maturity machinery (N_LAM_RWARM=4,
      GS_FREEZE_LIM=60, GS_REL_STABLE=2) - NOT the wide-open
      jfnk_probe settings. Question (serves K1): does the KNOP
      lagged field mature enough to admit JFNK (GS_FREEZE=F
      gate), and does JFNK change the known KNOP non-
      contraction outcome? Exit criterion: JFNK admission or
      sustained non-admission visible in OUTGEN, plus a named
      convergence/stall verdict.
- [~] Declared bounded investigation S1 - SOB-DIAG oscillation.
      Question: what makes the conventional population map
      non-contracting on the m15 structure? Evidence so far:
      nine-gate evaluator shows only ex_maxch>=eps failing
      (raw max correction ~6% vs EPS=5%) for ~250+ iterations.
      CORRECTION_SUM final state localises the instability:
      >0.1% demands exist ONLY at depths 1-13 (~28 rows/depth
      >10%); depths 14-70 are all below 0.1%. Depths 1-13 are
      the outer near-constant-V wind (sigma_true ~ -1 to -1.4;
      the near-singular Sobolev tau-floor regime).
      S1 finding 1 (22-Nov): the exempted churn DEVELOPED late -
      FeI dead rows ~0.5-3% through it240, CV Rydberg demands
      grew 5.75% (it286) -> 86% (it330) -> ~100% pinned (it374+,
      depths 1-11, cycling across 1s_nl levels). Not a static
      floor artifact; a growing regional instability.
      S1 finding 2 (22-Nov, user's exempted-rows question): the
      counted blocker set is a MIX - dead/marginal counted rows
      (SIII 1Se@13 7e-28, OIII 3Pe[2]@1 3e-14, NeIII ground@1
      1e-12, SIV@19 ~1e-8) AND genuinely populated rows
      (AlV lvl65@16 n=58 demanding +183%, MgV@14, OIV@9, NIV@8,
      He2 n=2@13). The instability crosses the 1e-30 floor
      boundary both ways - the exemption splits the metric, not
      the physics. Exempted CV Rydbergs (~1e-74 vs CV ground
      ~50 at depth 1) are individually negligible but are the
      sensor for a real defect, per the documented
      cascade-refill class. Counter-evidence note: AlV is the
      DOMINANT stage at depth 1-2 (~4.6e5); an unstable excited
      level of a dominant stage is closer to material.
      S1 finding 3: ES_ALI converges internally each call
      (0 unconverged freqs, dJ/J ~1e-4) so the field is not
      churning within a call - the non-contraction lives in the
      population update given a (nearly) fixed field. The
      persistent +6% INCREASE demand on SIII@13 for ~15
      consecutive iterations is same-sign marching, not a
      clean 2-cycle; the argmax then migrates across ~8 rows.
      NG fires with abstentions and the 75% dips land on NG
      iterations - plausible periodic re-excitation.
      Open discriminator: is the map discontinuous at the
      tau-floor branch switch (border-collision, phase-7
      classification for these same rows) or smoothly
      non-contracting? Different fixes. Pending: np-arm A/B
      on the same rows when 13564952 lands.
      Declared arm (22-Nov, user-approved): sob_diag_ngoff -
      warm restart of the stalled sob_diag it387 checkpoint
      (warm_restart_from m15_2p82d_VAL_sob_diag/cmfgen/2.82days)
      with DO_NG=F, NUM_ITS=60 additional iterations, otherwise
      identical VADAT (config:
      configs/dessart_m15_2p82d_VAL_sob_diag_ngoff.yaml,
      model_dir m15_2p82d_VAL_sob_diag_ngoff). Question: is the
      pinned ~6% counted-correction oscillation at depths 1-13
      driven by periodic NG re-excitation (the 75% conv dips
      land on NG iterations) or is the population map itself
      non-contracting? Exit criterion: corrections contract
      below EPS -> NG implicated, next arm is cold-start
      NG-off; still pinned -> NG exonerated as driver and the
      defect is map discontinuity (tau-floor border collision)
      or conditioning; partial decay then restall is also
      informative. Submitted 22-Nov: CMFGEN job 13567608 (small),
      OBS chained 13567609, tree 1fcfd3a8, pinned exe
      e031e8597875c541_9080452b4e65c32e.
      ARM VERDICT (22-Nov): CONVERGED. Warm-restart from the
      it387 stalled checkpoint with DO_NG=F: ex_maxch fell
      monotonically (21.1 -> 9.9 -> 4.3 -> 2.9 -> 2.3 -> 1.7 ->
      1.1 through it394) and the nine-gate evaluator passed at
      it409 - 'conv exit [MET] 100%, worst active margin
      -0.03 dex (ex_maxch) gates: all pass', exit reason 'clean
      convergence (MAXCH < EPS, exempted demand small)',
      exempted demand 3.9e-2% (the CV Rydberg churn vanished
      once the counted rows settled). 22 additional iterations.
      Verdict on the S1 question: NG periodic re-excitation was
      SUSTAINING the boundary-pinned oscillation at depths 1-13;
      the population map itself contracts. The stalled state was
      a limit-cycle held at the EPS boundary by the accelerator,
      not an absent fixed point. Caveat: warm restart rebuilds
      J (no FGINEXT_GS in the checkpoint) so the trajectory is
      perturbed - the deciding end-to-end test is the cold-start
      NG-off arm below. S1 closes as DECISION: production
      candidate is DO_NG=F; fix = disable NG for the SOB route
      (no code change needed - VADAT key exists).
      Declared follow-on (22-Nov): cold-start arm
      dessart_m15_2p82d_VAL_sob_ngoff_cold - identical sob_diag
      recipe (LTE_EST=T, DIAG, DO_JFNK=F) with DO_NG=F from
      scratch. THE SOB-leg acceptance test: if it converges,
      the conventional SOB pipeline cold-starts to a converged
      model. If it stalls at the same place, DO_AV=F (the other
      periodic accelerator, ITS/AV=8) is the next knob.
      Submitted 22-Nov: CMFGEN job 13569107 (small, NUM_ITS=1000),
      OBS chained 13569108, tree 8afaa8e3 +CHECKLIST, same pinned
      exe, VADAT identical to the ngoff warm arm.
      Exit criterion: named root cause + FIX commit with
      regression, or an evidenced BLOCKED verdict. No
      tolerance relaxation, no new clamps.
- [ ] Declared bounded investigation NG1 - NG sawtooth and
      adaptive acceleration (24-Nov-2026). Question: on the
      eps1em3 production deck (m15_2p82d_VAL_sob_jfnk_eps1em3,
      DO_NG=T, NG_ABSTAIN=T, ITS/NG=15, EPS_TERM=1e-3) the
      accepted NG kick every ~16 counters rebounds the
      demanded correction by ~55-60x, then the unaccelerated
      map re-damps to ~78% of the pre-kick demand per cycle -
      is the kick a wasted transient or a net-positive
      slow-mode extrapolation, and would a Magritte-style
      settle-gated or regularised fit do better on wall clock?
      Evidence so far: SCRTEMP replay of 14 historical events
      reproduces committed kicks to 2.1e-13; Gram condition
      27-190 (not singular); affine coefficient L1 ~1e3 late;
      successive 4-record shadow predictions move 2.4-3.1x
      more than the ordinary iterates, so a Magritte settle
      gate would have deferred all six late events examined.
      Evidence plan: (a) diagnostics-only shadow-window report
      (NG_SHADOW_DIAG, default F) logging per-band fit
      coefficients, conditioning and settle ratio for growing
      windows N=4..16 - the deferred 16-Oct adaptive-gate
      assessment required a diagnostics-only first step;
      (b) matched warm-restart A/B at identical checkpoint,
      binary and EPS_TERM: DO_NG=T cadence vs DO_NG=F vs one
      guarded variant if (a) motivates it; judged on
      iterations and wall time to the same nine-gate exit.
      Parameter proposal (per AGENTS 9.1): NG_SHADOW_DIAG,
      default F, diagnostic-only. User decision: opt in to
      NG fit telemetry on a live deck; no physics change,
      no behaviour change when F. No auto-set. Interactions:
      none beyond requiring DO_NG=T to have events to
      inspect. Exit criterion: settle/conditioning data
      identify the kick mechanism and either justify a
      guarded-acceptance FIX with a matched A/B run, or the
      cadence is shown net-positive and NG stays as is -
      recorded here either way.
      Status 24-Nov-2026: diagnostics landed, commits
      39bfc2c6..b81a2f54. NG_SHADOW_DIAG
      (default F) reports per band and window N=4..nwin the
      affine coefficient L1 norm, implied prediction move,
      predicted-residual ratio and Magritte settle ratio, from
      a weighted Gram/DGESV constrained fit using the
      production banding and 1/x weighting. test_ng_shadow.py
      16/16 (serving-stub SCRTEMP, numpy-verified numbers,
      bit-identical POP output F vs T); test_vadat_keys.py
      87/87; test_ng_abstain.py 17/17.
      Telemetry arm 24-Nov-2026:
      dessart_m15_2p82d_VAL_sob_jfnk_ngshadow - identical
      eps1em3 deck warm-restarted from that run's 2.82d
      checkpoint with NG_SHADOW_DIAG=T (outer-repo config
      commit 90ac2c0e; binary d47fe28b incl. flag). CMFGEN
      job 13646902 (+OBS 13646903), shares the lamesF
      restart source. A/B control same checkpoint, binary
      and EPS_TERM: dessart_m15_2p82d_VAL_sob_jfnk_ngoff
      (DO_NG=F; config commit 58cc2d65; CMFGEN job
      13646953 +OBS 13646954) - the ngshadow arm doubles
      as the accelerated side since the flag is
      POP-bit-identical. Collect per-event
      settle/conditioning series + iterations/wall time
      to exit, then decide whether a guarded-acceptance
      NG_ADAPT candidate justifies the third arm.
- [ ] Declared bounded investigation K1 - KNOP formal
      non-contraction. Question: can the corrected-exit build
      contract the lagged field from a cold start, and what
      breaks the period-2 orbit on inner rays 1-3 (~88K bad
      cells/sweep, exit=2, NI~300-340)? Evidence plan: (a)
      fix82-class honest-exit arm resubmitted on general
      partition (13158010 died on the interactive CPU cap)
      with KN_RDUMP_LS armed on the inner-ray population;
      (b) offline replay of the deferred corner-mode
      acceleration candidates (guarded Aitken / OMEGA_RAY /
      safeguarded Anderson / KN_ORD_P2DAMP) on validated decks
      per the phase-2.30 protocol before any live trial;
      (c) bounded cold-start arm. Exit criterion: converged
      cold-start run or an evidenced REJECTED/BLOCKED verdict.
      No accelerator enabled in production by default.
- [~] Lagged-store lifecycle audit (22-Nov-2026): user-requested
      source+runtime investigation of I_P/I_M_DEPTH_NEXT across
      sweep reads/writes, refresh re-solve, commit, heal, seed and
      checkpoint paths, with literature comparison (Knop Paper II
      converges the coupled BVP inside each formal solve; the dev
      lagged store is an explicit-source shortcut). New production-
      object battery test_knop_lagged_map.py (26 checks, 25 PASS):
      cold/warm/poisoned seeds on single-band, two-band-reversal and
      strong-compression profiles all reach finite seed-independent
      committed fixed points (KN_Q_RESID q max = 0.0); diverged
      poisoned seeds heal without committing. ONE CONFIRMED DEFECT:
      a converged certified solve can be permanently jump-vetoed by
      a poisoned retained generation - the store keeps NaN/Inf
      forever while cert reports FULL (HEAL_STORE runs only on
      diverged rays; pinned failing as sdab). Fix candidate
      deferred: retained-side sanity in FG_JUMP_MAX or an entry-side
      heal. docs/audits/CMFGEN_KNOP_LAGGED_STORE_POISON_AUDIT.md.
      Production-scale rho~0.97 stall not reachable on the
      manufactured grid; K1 stays open.
- [x] Discrete-operator verification (22-Nov-2026): new
      test_knop_operator.py + DBRC dump in drv_knop_ray check the
      refresh-side operator against the mathematics it discretises,
      not just sweep/refresh self-consistency. 32/32 PASS across all
      three GAM regimes - all-expanding outflow (gposo),
      all-compressing inflow (gnego, the pure lagged-store regime),
      continuous sign reversal (gcr) and a discrete
      expansion/zero/compression combination (gzb):
      A0=exp(-DTAU) to 7e-18, A0+A1+A2=1 to print precision
      (constant-source equilibrium preserved), Q*CHIH*DLNM=signed
      GAM to ~5e-9, CHIH=chi+|GAM|/DLN at coupled cells of EITHER
      sign (expanding cells couple to the live previous-frequency
      donor - verified, not assumed), CHIH=chi at uncoupled cells,
      DTAU=trapezoid(CHIH) exact under KN_METHOD_SAV='ZERO',
      SR*CHIH=analytic ETA to 5e-8, field finite/non-negative with
      vacuum outer BC and thermalisation onto the modified source
      at thick cells (1%). The ZERO-method DTAU is plain trapezoid,
      ~6.5% above the analytic int chi dz on the manufactured grid
      (discretization, not defect). gcr outer boundary cell i=1 is
      uncoupled at all frequencies (KN_COUPLING_OK drops both sides
      consistently).
- [~] BLANK spectrum + Luc comparison: renderer verified 22-Nov
      (RESC_JMAX fix; obs_fin_10 valid). First converged
      conventional-solver SOB spectrum exists 22-Nov from the
      sob_diag_ngoff arm (warm-restart convergence, OBS job
      13567609, obs_fin_10 133179 pts all finite). vs Luc's .fl
      (identical compare_finished_vs_luc.py methodology, run on
      the ngoff path): band medians FUV 0.990, UV 1.056, blue
      1.015, opt1 0.998, opt2 0.950, red 0.935;
      rms3000-1e4 = 0.042 (unconverged parent scored 0.327 with
      feature ratios 8-142x - its obs_fin_10 of Nov 22 02:57 is
      a stalled-state render, not a route product);
      feature pk/cont ratios vs Luc: C IV 0.94, C III] 1.12,
      Mg II 0.95, H-alpha 1.06, H-beta 0.93, H-gamma 0.84,
      H-delta 0.87, Ly-alpha 0.37, Na D 0.63. Overlay:
      cmfgen_validation/plots/validation/ngoff_vs_luc.png.
      compare_finished_vs_luc.py MODELS updated: ngoff entry
      added, parent relabelled UNCONVERGED. Acceptance still
      open: same test on a COLD-START converged model.
      Oxygen-residual audit (22-Nov, user question 'is the
      atomic data correct'): all O atom files (O2/OIII/OIV/OV/
      OSIX osc, coll, photo, IN) byte-identical to upstream
      cmfgen_reference - no local corruption. Checked f/A at
      the discrepant multiplets against standard values:
      O IV 1338.6 f=0.117, 1343.5 f=0.104, UV2 multiplet
      3381-3425, O III 5592 f=0.15, O V 1371.3 f=0.158,
      O VI 3811/3834 f=0.23/0.11 - all correct. The deviations
      are strength/shape, not data: O VI doublet present in
      both (Luc's slightly stronger); the ~5600 bump is a
      smooth ~8000 km/s blend, not a discrete mis-ID'd line.
      Regional medians model/Luc (ngoff): 1320-1400=0.81,
      3180-3250=1.29, 3380-3430=1.22, 5550-5700=1.23,
      3790-3850=1.01 - vs the UNCONVERGED parent (13.7, 0.76,
      6.7, 0.91, 0.63) showing the residuals are strongly
      population-state dependent and collapsed toward Luc at
      convergence. Suspects ranked: FIX_T hydro-T (12.7k K CSM
      -> 86.7k K ejecta) vs Luc's solved electron T; the
      non-mono velocity/line-overlap reprocessing itself;
      code-version microphysics. dFR_DATA empty (STAU) so no
      formation attribution this run. Pending: np-arm OBS for
      the same-data different-solver A/B.
- [ ] Declared audit A1 (22-Nov, user-directed, long-standing):
      WHY do line SHAPES differ from Luc's .fl even where
      strengths are close? Observed so far: model O IV/O III
      features broader + red-smeared; smooth ~8000 km/s ejecta
      bump at 5600 absent in Luc; UV continuum 0.81 of Luc;
      H-alpha close. Testable branches: (a) velocity mapping -
      our vel_law=14 reconstructed V(R) + ND=70 downsampled
      grid vs Luc's native structure; (b) temperature -
      FIX_T=T pins converted hydro T (12.7k->86.7k K) vs Luc's
      solved electron T - user authorized a FIX_T=F attempt;
      (c) rendering - observer-spectrum convolution/resolution
      settings; (d) populations - np JFNK A/B when OBS lands.
      Evidence plan: quantify feature widths/offsets in km/s on
      isolated lines (Ha, C IV, O VI, O IV 3400, 5600 bump);
      read .fl metadata; compare SN_HYDRO_DATA V/rho/T against
      the conversion source; then the FIX_T=F warm-restart arm
      on the converged ngoff checkpoint. Exit criterion: named
      dominant cause + fix arm verdict, or evidenced BLOCKED.
      Measurements (22-Nov, scripts in /tmp/line_shape_audit*):
      .fl col1 = wavelength A (15.7 A - 74 um), vt50 in the
      name = vturb 50 km/s matching our VTURB=50; TDOP here is
      the thermal-Doppler T (1e4 K), negligible width. Velocity
      scale AGREES: blend-peak velocities match (C IV +17.5k,
      O IV 3400 -17.2k, He II 1640 blue edge -17k km/s both);
      the ejecta V mapping is not the problem. Mismatch is in
      the narrow cores + blend weights: H-beta core FWHM 330
      (Luc) vs 960 (ours), H-alpha 438 vs 560, He II 4686 300
      vs 420, C IV 3980 vs 4560 - ours systematically broader;
      He I 10830 peaks opposite sides (+2852 vs -2992); O VI
      5290 strong in Luc, absent in ours; O V 1371/Si IV 1400
      blend maxima at very different velocities. Structure:
      CSM spans depths 1-54 (V 49-195 km/s, T 12.7-30k, rho to
      8e-14), shell 55-63, ejecta/CDS 64-70 (V to 7074, T to
      86.7k, sigma_true -157). Interpretation: cores form in
      the slow CSM; their width is set by the emission region
      (ionization) + electron-scattering wings (T_e) - both
      temperature-structure effects, pointing at FIX_T.
      Arm declared: dessart_m15_2p82d_VAL_sob_ngoff_thinT
      (config hash d0ef3098836f3f20) - warm restart of the
      converged it409 ngoff checkpoint, fix_temperature:false
      thin release (FIX_T=F + FIX_INB_T=T + FIX_X_DPTH auto at
      the tau_es=2/3 photosphere), DO_T_AUTO=T + DO_LAM_IT=T +
      DO_LAM_AUTO=T per the SOB_l14_freet recipe, DO_NG=F held,
      DO_JFNK=F (required for FIX_X_DPTH), 200 extra its.
      First submission CANCELLED (13570910/13570911): the
      launcher's NON_MONO_SOB recommended defaults inject
      FIX_T=T into vadat_overrides for unset keys, which beat
      the prepare_VADAT thin-mode stamp - the live VADAT had
      FIX_T=T + FIX_INB_T=T + FIX_X_DPTH=24, and RD_FIX_T=T pins
      T at ALL depths in solve_for_pops.f:1014 regardless of the
      inner-boundary mask. Resubmitted with explicit FIX_T=F +
      FIX_INB_T=T in vadat_overrides (13570980/13570981);
      FIX_X_DPTH left unset so prepare_VADAT auto-fills the
      tau_es=2/3 photosphere index (computed 24 on the first
      attempt = pin depths 47-70, solve T at 1-46).
      VADAT verified live: FIX_T=F, FIX_INB_T=T, FIX_X_DPTH=24,
      DO_T_AUTO=T, DO_NG=F. Early behaviour (it410-420): the T
      release IS active - FULL mode prints 'Temperature held
      fixed at 24 depths at inner boundary'. But each released
      FULL step demands an enormous change: it417 1.5e5% (BA
      rejected), it420 1.0e7% capped to 100% applied (ba=HELD),
      driven by CSIX lev1 at depth 46 (the release boundary)
      and CV Rydberg runaway at depth 18 - the SAME dead-row
      class as the parent stall. The hydro-pinned state is far
      from a free-T equilibrium; the run ping-pongs
      LAMBDA_TFIX <-> FULL with capped T steps. Watch: does
      the demand contract or cycle.
      VEL_LAW audit (22-Nov, user: 'Luc uses vel_law 12'):
      VADAT VEL_LAW selects the runtime R/V/SIGMA build; YAML
      vel_law selects the converter SIGMA formula. Paths:
      V12 = RV_SN_MODEL_SNIIN - reads input_hydro.dat, grid
      uniform in log(tau_es) + NBND_INS boundary packing,
      LIN_INTERP V in R + two-sided fd sigma. V14 =
      RV_SN_MODEL_SNHYD TAU_OPT=1 - reads SN_HYDRO_DATA,
      blended tau-Rosseland/log-R grid, monotone-cubic V +
      analytic sigma. V13 = SNHYD TAU_OPT=0 = V12 recipe on
      our data (the clean A/B). Verified at recipe level by
      offline emulation on the same 40-pt source: V12 grid
      puts 57 of 70 depths in V<100 km/s (4 in the ejecta);
      our V14 RVTJ grid has 22 in V<100 and 14 in the ejecta.
      CDS sigma spikes to -130 (V14 cubic) vs ~-53 (V12 fd).
      The two paths sample the line-forming regions
      completely differently - a live suspect for the CSM
      narrow-core broadening and the ejecta blend weights.
      Luc's own input_hydro.dat is NOT in the tree (both
      candidate files are our converter output) - only the
      .fl exists, so source-vs-recipe cannot be separated;
      V13 isolates the recipe. Arm declared:
      dessart_m15_2p82d_VAL_sob_l13_ngoff_cold (config hash
      b3c90852d4bd8291) - cold start LTE_EST=T, identical deck
      to VAL_sob_ngoff_cold but VEL_LAW=13 (tau_es grid, V12
      recipe on our data). Exit criterion: converged V13
      spectrum's narrow-core FWHMs and blend weights vs the
      V14 ngoff model's - cores sharpen toward Luc's ~330-440
      km/s -> grid recipe implicated; still ~1000 km/s ->
      grid exonerated for cores, look to T/populations.
      Submitted: CMFGEN 13572032 + OBS 13572033 (VEL_LAW=13
      preserved in VADAT; n_ob_ins forced 30->0 so the RDINR
      grid matches MODEL_SPEC ND).
      ARM VERDICT (23-Nov): FAIL - non-convergent dead-zone at
      the underresolved CDS/ejecta jump. From it1 the FULL-mode
      BA solve demands saturated corrections (MAXCH=1e7%) on the
      whole CV ion (all ~25 superlevels, drain to C3+; the model
      has no C5+ species) at cliff depths 51-52 and 60-62
      (V 70-290 km/s, sigma_true -2.0 to -4.7). LAMBDA steps run
      benignly (demand ~2%) but every cnt_lam_escape into FULL
      produces ba_damped_full -> ba=HELD with maxatt=0 (nothing
      applied); held-reject loops repeat until ba_sat_demote.
      Root cause is grid, not sigma: V14 cells carry sigma_true
      -4.8 to -5.3 across the same band and converge - the V13
      tau_es recipe crosses the CDS in one cell pair (tau_es
      17.7->160.7, V 2380->6375, ~6 cells for V=100-2400 km/s vs
      V14's ~24 graded cells), so depth-coupled BA derivatives
      saturate. Same failure class in FRESH_sob_jfnk (held
      it295-302), LUC_exact (V12, held from it68, plus 3364
      flux-accuracy warnings), FRESH_knop_jfnk (still in lambda
      warm-up). Stored STEQ residuals in the band are
      large/negative under V13/V12 vs positive/smooth under V14.
      This is the dead-row class of the parent stall - the
      AGENTS.md 6 trigger case: if a V12-recipe run is still
      needed, revisit solveba_nonmono BA damping. Consequence:
      vel_law 14 restored as the pipeline default 23-Nov-2026
      (AGENTS.md 1.2); the arm's FWHM question is answered for
      convergence before spectra - the V12 recipe cannot produce
      a converged model here.
      POST-VERDICT (23-Nov): dead conventional arms cancelled -
      LUC_exact (13573428+OBS 13573429, held it68-212, 68 held
      FULLs) and l13_ngoff_cold (13572032+OBS 13572033, held at
      it372, 56 held FULLs; fin0 13591841 finalized, no converged
      write). VADAT-knob audit for the CV dead-zone: no key
      triggers or cures it - DIE_C* = F,F in every arm,
      DC_METH=R in both a failing (LUC) and a pre-wall (knop_jfnk)
      arm, MAX_LAM/NUM_LAM/FIX_BA identical, DO_NG on AND off both
      held. Only VEL_LAW (grid coordinate) separates failing from
      converging arms. COUNTER-EVIDENCE on the solver side:
      FRESH_sob_jfnk (same V13 grid) was carried THROUGH the
      held-zone by JFNK and reached [MET] 100% all-gates at
      it413 - cancelled by user before the converged write; SCRATCH
      checkpoints at ~it413 persist (BAMAT/BAMATPNT/CSCRATCH
      09:43), so a warm restart can finish it in a few its and
      yield the only V13 converged model + spectrum for the A/B.
      FRESH_knop_jfnk cancelled at it132 (still lambda warmup, 0
      held) - restartable, low value. Remaining live: ngoff_cold
      V14 conventional (13569107, it375, 89%) = the conventional
      cold-start arm under the restored default.
      RENDER (23-Nov): all three V13 arms re-submitted as *_fin0
      warm restarts (warm_restart_from the parent's last SCRTEMP)
      and stopped via the pipeline's own graceful-stop route -
      0 [NUM_ITS] written into the live IN_ITS; CMFGEN re-reads
      IN_ITS each iteration (cmfgen_sub.f ~5759), clamps
      NUM_ITS_TO_DO to 1, takes one LST_ITERATION, writes
      POP*/RVTJ/OBSFLUX and exits clean. init_num_its: 0 in YAML
      does NOT do this - FINAL always writes NUM_ITS=10000 and
      the batch guard rewrites a startup 0 back to FINAL_NUM_ITS.
      States at write: l13 it375 78%, sob_jfnk it416 78% [HOLD]
      (the it413 [MET] all-gates pass was transient; three more
      lambda its regressed to the held dead-zone), knop_jfnk
      it134 44% (KNOP_IMMATURE, needed DO_KNOP_SPEC_GATE=F ->
      marked diagnostic render, obs 13592257). Spectra:
      outputs/dessart_m15_2p82d_{VAL_sob_l13_ngoff_cold,
      FRESH_sob_jfnk,FRESH_knop_jfnk}_fin0_CMFGEN_2.82408d.dat/png
      - all non-converged-state renders, usable for the V13 vs
      V14 shape A/B but not as converged validation products.
      FWHM A/B on the renders (23-Nov, /tmp/v13_v14_fwhm_ab.py,
      narrow cores in km/s):
        Hbeta   Luc 290 | v14_ngoff 990 | v13_l13 750 |
                v13_sobjfnk 590 | v13_knopjf 270
        Halpha  Luc 388 | v14_ngoff 550 | v13_l13 562 |
                v13_sobjfnk 400 | v13_knopjf 350
        HeII4686 Luc 270 | v14_ngoff 400 | v13_l13 440 |
                v13_sobjfnk 290 | v13_knopjf 170
      Verdict: grid recipe EXONERATED for the ~3x core broadening -
      same-solver V13 vs V14 narrows Hbeta ~24% (990->750), leaves
      Halpha unchanged (562 vs 550), HeII slightly broader. Mild
      effect at most, far short of the 3x gap; and no converged V13
      state exists to test harder. Also noted: v14_np_jfnk wrote
      UNCONVERGED state (budget exhausted, last gate ex_maxch
      0.02 dex short) yet matches Luc's cores (330/338) - the
      solver/population state moves cores more than the grid does;
      and thinT's free-T Hbeta blew out to 4190, so T-release did
      NOT sharpen cores in this realization. CIV1549 stands apart:
      Luc 2725 vs ~1600-1800 everywhere - a different-feature
      problem, not a width problem.
- [x] Production-code grid A/B/C at ND=70 on the real structure
      (unit_tests/work/grid_ab, scratch): drv_snhyd TAU_OPT=0/1 on
      SN_HYDRO_DATA + drv_sniin70 (production RV_SN_MODEL_SNIIN) on
      repackaged input_hydro.dat. Measured:
        points at V<100:        V12 55 / V13 58 / V14 22
        points 50-300 km/s:     V12 38 / V13 27 / V14 33
        points 300-2000 km/s:   V12  3 / V13  2 / V14 16
        points >2000 km/s:      V12  7 / V13  6 / V14 14
        depths tau_es 0.1-10:   V12 16 / V13 12 / V14 34
        CDS sigma min:          V12 -36 / V13 -63 / V14 -82
                                (source true sigma -157)
        max adjacent |d sigma|: V12 27 / V13 58 / V14 60
      Verdict: V13 ~ V12 in placement (confirms the V12 recipe on
      our data) but in velocity space V14 is NOT under-resolved in
      the CSM (33 pts, med dV 4.4 km/s) and is 5-8x better resolved
      through the shell/CDS/ejecta and 2-3x better in the
      line-forming tau band. The tau_es grids dump ~80% of depths
      into the optically thin outer CSM (V<100, tau<0.1). "Optimal"
      splits by feature: V12/V13 for raw CSM point count, V14 for
      CDS spike capture + ejecta/shell coverage. Spectrum decides.
      NEW DEFECT FOUND: the production V12 reader's row read is
        READ R_HYDRO(I), V_HYDRO(whole array), XD, XT, XL, XTAU,
        TAU_HYDRO(I)  -> NHYDRO+6 values per row
      so our converter's 7-column input_hydro.dat cannot feed
      VEL_LAW=12 (EOF at row 2; confirmed with both sniin_new and
      sniin70 drivers). The file also must be ordered inner->outer
      or the grid returns ascending (inner R at depth 1). V12 has
      never successfully run on our converter files - LUC_l12's
      early death is consistent with this. If V12-vs-V14 ever
      needs a real run, the converter format needs fixing first.
      RESOLVED (22-Nov-2026): two coordinated fixes. (1) Reader:
      V_HYDRO -> V_HYDRO(I) in rv_sn_model_sniin.f. The whole-array
      read could never coexist with SET_ABUND_CLUMP's 4-value row
      read of the same file (array-per-row gives it rho=V[2],
      T=V[3]), so the 7-column table is the only self-consistent
      contract; the unindexed read was a latent upstream bug that
      killed every V12 run at ~row 7 (forrtl severe(24): job
      13572087, identical to LUC_l12's death). (2) Writer:
      prepare_VADAT.py now emits rows inner->outer in both
      input_hydro.dat branches (was outer->inner); the reader's
      reversal and abund_clump's reversed-index read both require
      increasing R. generate_eddfactor.py's tau read made
      order-aware. Verified: drv on the real 64-row m15 file gives
      decreasing R grid (1.253e5 -> 1.455e4), V endpoints exact,
      sigma -1.05 outer / -31..-36 inner; test_rv_sniin.py updated
      to the 7-col contract, 14/14 pass. The array-per-row
      convention in unit_tests/work/*input_hydro.dat is a stale
      test artifact of the old quirk.
      FRESH-START DEFECT (22-Nov-2026, third): AUTO_ADD defaults
      .TRUE. (rd_control_variables.f:3100) and AUTO_ADD_ION then
      requires RVTJ/EDDFACTOR_INFO/EDDFACTOR from a previous model
      (job 13573027 died IOS=29 on EDDFACTOR_INFO). The
      m15_2p82d_VAL_sob_diag* runs never hit this because they are
      warm restarts inheriting those files; every production VADAT
      carries AUTO_ADD=F explicitly. A true cold start with
      LTE_EST=T still needs AUTO_ADD=F -> the NON-GRID init path
      generates the ion _IN files itself. Luc's VADAT lacks the
      key (his source may not have it); configs based on his file
      must set it explicitly. Fix landed in
      configs/dessart_m15_2p82d_LUC_exact.yaml; job 13573428 then
      passed init and is iterating (DIE_CIII/NIV/OV linked via
      allow_dielectronic; DIE_ScI/Sc2 forced F,F - no Sc die_data
      at any vintage).
      Blocked only on a
      converged conventional-solver model from either route;
      compare_finished_vs_luc.py supplies band ratios, log-rms
      and feature contrast vs the published .fl.

## BLANK spectrum faint-continuum root cause (22-Nov-2026)

- [x] Forensic isolation on the 13560175 render (obs_fin_10 median
      optical continuum 8.75e-9 vs 1.00e-6 in the same-populations
      Nov-11 render, ~100-270x faint, wavelength dependent, red/IR
      unaffected). Old-exe/new-populations arm (13562156, Oct-16
      frozen exe, Nov-11 params) reproduced the Nov-11 spectrum
      exactly (median 9.99e-7) - populations and OBS params
      exonerated; defect is in the Oct-16 -> Nov-22 obs code.
- [x] Stored-field comparison (EDDFACTOR, 71-double records):
      new-exe J pinned at 9.948e-4 - just under RESC_JMAX=1e-3 -
      across FL~20 -> 0.3 while the old exe legitimately reaches
      3.7e-2 at FL=10 (real peak ~0.122). J_RESCUE=T was tripping
      on physically valid fields and J_RESCUE_RESTORE latched the
      stale previous-frequency field down the warm-start chain.
- [x] A/B proof: identical exe+params+pops, J_RESCUE=F (13563
      series, obs_noresc) recovered the old spectrum to 0.004%
      (median 9.0532e-5-band, Halpha 1.22e-6), zero rescue events;
      the J field then matched the old exe (3.70e-2 at FL=10).
- [x] Param-file pitfall recorded: batobs.sh regenerates
      CMF_FLUX_PARAM from CMF_FLUX_PARAM_INIT at run start - edits
      to CMF_FLUX_PARAM alone are silently reverted. The earlier
      "J_RESCUE exonerated" arms (obs_noresc/norestore first
      submissions) still ran T [J_RESCUE] and were false negatives.
      Correct override syntax is `F [KEY]` in CMF_FLUX_PARAM_INIT,
      never `KEY=VALUE`.
- [x] Fix landed and VERIFIED: RESC_JMAX default 1e-3 -> 1e2
      (cntrl_var_mod + rd_cmf_flux_controls, comments in
      COMP_JCONT_V4 updated). Rationale: the flat rail was
      calibrated on mc170 (legit peaks ~5e-5); m15 legit peak
      ~1.2e-1. 1e2 sits ~3 orders above the m15 peak, ~23 below
      the observed runaway (1e25); the 1e4*JNU_PREV relative
      bound remains the working detector. Verification arm
      obs_rmax (job 13564117, T [J_RESCUE] + rebuilt exe):
      zero rescue events, J field unclipped (3.7048e-2 at FL=10,
      peak 1.218e-1), obs_fin_10 opt-median 9.9894e-7 -
      identical to obs_noresc and obs_oldexe at every checkpoint
      2000-50000 A (ratio 1.0000). Unit batteries
      test_jres_validity 7/7, test_jres_events pass.
      Source commits 310a9364 + f2471abf + cc5eafc2.
- [ ] Outstanding: a genuinely-divergent arm proving the rescue
      still fires under the new rail; the tauto model's
      sigma-clamp boundary artifact is no longer clipped by the
      flat rail (reverts to stock pre-rescue fields) - per-cell
      J/J_eq bound is the deferred model-independent criterion
      (threshold-attribution audit sec. 6a). CHECKLIST 13560175
      entry remains valid for undercoverage but its spectrum is
      superseded (faint, do not use).

## Ray-83 gradient and sweep-order investigation (22-Nov-2026)

- [~] User-authorized KNOP formal / USE_LAM_ES=T / DIAG investigation.
      Question: does the grazing-ray failure arise near a projected-gradient
      zero, or from lagged source coupling in the spatial-major march?
      Evidence plan: compare Dessart Sec 3.1 with the signed CMF operator;
      inspect captured ray geometry; compare the existing ordered and
      frequency-major production sweeps using the same input deck.
      Exit criterion: an evidenced correction with regression and replay
      checks, or a precise missing-state blocker and a deferred decision.
      No SIGMA/GAM clipping, tolerance relaxation or P2-on model is allowed.
      Ray 83 is not among the three captured rays; its exact trajectory
      must not be inferred from ray 81's deck.
- [~] Original-object sweep comparison 13557757 completed: frequency-major
      gives FULL / zero bad cells on 81/82/84, residuals 5.805058e-21,
      4.945920e-21, 6.681912e-23; repeated fields byte-identical.
      Both 81/82 have GAM > 0 everywhere. No SIGMA clamp applies to this
      evidence. Ordered replay still STAG on 81/82; production-trajectory
      M/IM_NI discrepancy remains open. Parser tests 29/29 pass.
- [~] One-line dispatch candidate implemented in KNOP_RAY_SOLVE_V1:
      only wholly non-negative-GAM rays use the existing frequency-major
      solver. Negative/mixed-sign selection and all physics/gates stay
      unchanged. Regression 13557822 failed 4/8 before the edit, including
      positive/zero dispatch and field-equivalence checks. Full build,
      post-edit regressions and live ray-83 verification still owed.
- [~] Earlier capture audit 13557300 validated selected inputs but failed
      before baseline replay on missing MKL runtime loader path; keep the
      original knray_replay_p230_o3 binary for a loader-corrected baseline.
      Do not rebuild that baseline against the dispatch candidate.
- [~] Original binary baseline rerun 13558443: ray 84 repeated fields
      match bytewise (sha256 111602c0...); 81/82 abort on strict runtime
      pass-count mismatch. Ray 81 requests pass 3 where the recorded
      schedule is exhausted. The strict guard is retained. Captured
      trajectory fidelity FAILS on the stalled rays; dispatch comparisons
      remain tests of the replay inputs, not proof of a live-model fix.
- [~] Full clean build 13558165: compiled and linked clean; aggregate
      suite reported 9 failures, none from the dispatch edit. Mine fixed
      (gzzo fixture, three stale pins); rest pre-existing environmental
      or infrastructure gaps (missing run fixtures, es_conv GUARD_SRCS
      link list, phase227_io harness args, par_opac missing .o, DGELSY
      under mkl/2026.2, se_ba_non_therm coverage manifest).
- [x] Candidate A/B 13559005: six focused batteries rc=0 (oneway,
      lagged_field, gs_drift, obs_dispatch, replay_p2, knop_ray).
      Base-mode deck replays, original exe vs candidate exe: rays 81/82
      STAG at 1.9e-2/4.5e-2 under the original reach 5.8e-21/4.9e-21
      FULL under the candidate; ray 84 FULL either way. Dispatch count
      confirms frequency_rays=1 per ray. Ordered-retention path verified
      only on manufactured mixed-sign cases; live ray-83 capture still
      owed. Verdict: production-object PASS, live-model open.
- [~] SOB BLANK render 13557228 FAILED (exit 3): OBS_FRAME_SUB_V9 rejected
      red-edge CMF undercoverage, required 2.93133478e-3 below grid MIN_CF
      3.0e-3. Runtime BLANK controls and NUM_ES=2 were correct. Continuum
      J_RESCUE warnings also occurred and require a separate quality check.
      Preserve failed products; investigate frequency-coverage padding and
      rerender via run_model_grid.py without weakening the coverage guard.
      Rerender submitted as 13558872 with cmf_flux_param_overrides
      MIN_CF=2.0D-3 (renderer grid only; model VADAT MIN_CF untouched).
      OBSFRAME/spectrum validation owed on completion.
- [~] 13558872 ALSO failed on the same guard: required/edge ratio was
      identical (0.9771) in both jobs, so the undercoverage is structural
      and MIN_CF cannot repair it. Re-review found the defect in
      cmf_flux_sub_v5: NU_MIN_OBS = MAX(NU(NCF-3),
      NU(NCF)*(1+2*VINF/c)) is sized by VINF, which stays V(1) on the
      stock-J BLANK path (ALLOW_NONMONO_OBS escape). Here V(1)=48.9 km/s
      (CSM) while the ejecta reaches 7073 km/s, so the OBS red anchor
      landed inside the required Doppler range by construction. The
      sob_diag render passed on the same model only because
      DO_NONMONO_SOB already fed VINF=MAXVAL(|V|). Fixed by gating the
      MAXVAL sizing on N_NONMONO_OBS>0 as well (field-driven, same
      precedent as the INC_GRID rejection; commit 86d21335, pin added
      in test_obs_nonmono_dispatch, ea64b5fc). Rerender 13560175
      submitted with the rebuilt exe; OBSFRAME validation still owed.
- [x] Census defect found and fixed in the same re-review: KN_PA_*
      channel arrays allocated only under KN_Q_RESID but kn_census_row
      reads them whenever KN_RDUMP_LS* is armed, so every KNRAY_CENSUS
      row carried uninitialised channel columns (the 4.2e-322 values in
      the refresh-249 ray-83 row). Allocation now covers either
      consumer (commit 8204b8c0).
- [x] Ray-83 production capture landed: probe 13560532 (warm restart
      from the rd it250 checkpoint, NUM_ITS=1, dump rays 81/83/84,
      exe at ae48088c). Census r1: rays 81/83/84 all exit=1 FULL with
      zero bad cells, residuals 7.7e-21 / 6.6e-18 / 1.3e-22. Deck
      parse: ray 83 has gam_min=1.096e-9 > 0, zero compressing cells
      - a wholly expanding ray, so the dispatch reroutes it like
      81/82. The NI=3 grazing geometry was incidental; the live-model
      dispatch verdict is now VERIFIED for the expanding population.
      Remaining bad cells (135579 total, fg_rf=4.05e-2) sit on inner
      rays 1-3 (nbad ~25K each, exit=2, NI~300-340) - the population
      that plausibly carries GAM<0 cells and retains the ordered
      march; a distinct mechanism, still open.
- [x] SOB BLANK render 13560175 COMPLETED (8:22, exe at 86d21335):
      obs_fin_10 written (161423 continuum frequencies, 36322 rows
      all positive), OUT_FLUX/OBSFLUX produced, zero undercoverage
      and zero J_RESCUE events in batobs.log - the VINF fix removed
      the edge-driven continuum-J divergence as well. The standing
      obligation "regenerate the SOB-population spectrum with BLANK"
      is satisfied: populations untouched (byte-identical staging),
      GLOBAL_LINE=BLANK, original SOB render preserved in
      m15_2p82d_VAL_sob_diag.

## Phase 2.30 P2 production-operator evaluation (01-Dec-2026)

Offline A/B replay evaluation of KN_ORD_P2DAMP against the
confirmed ray-81/82 period-2 limit cycle. Damper stays F; no
production change.

- [x] Part 1 source verification: KN_ORD_P2DAMP ->
      FG_ORD_P2DAMP (module DATA default F), KN_P2_TRACK/
      KN_P2_CLEAR/KN_P2_MID shared at all three march sites,
      P2_SGN/P2_FLIP/P2_DMP allocated per solve (per-thread),
      phase-2.29A stale-state fix present.
- [x] Part 2 replay-only P2 override: `p2=<v>` keyed token in
      drv_knop_replay.f90 scanned over all args >=3, malformed/
      duplicate rejected before deck open, applied as the single
      write to FG_ORD_P2DAMP after deck restore; RESULT echoes
      p2_deck/p2_eff/p2_ovr. Composes with omf. VADAT reader and
      production default untouched (pinned).
- [~] Part 6 mapd prototype (driver side only): implemented,
      but NOT validated as R(X). Sparse-store seeding need not
      reproduce an arbitrary full candidate, and KR_DBG output
      is not a raw pre-relaxation VNEW capture. Do not use this
      mode to accept the damper until these gaps are closed.
- [x] Tests: test_knop_replay_p2.py 25/25 PASS (override pins +
      arg-layer rejection); test_knop_p2_state.py 19/19;
      test_knop_ray.py 49/49 unchanged. test_knop_lagged_field
      carries a PRE-EXISTING stale pin (expects RES.GT.HUGE;
      source moved to .NOT.(ABS(RES).LE.HUGE) at 98acf207,
      ancestor of baseline) - not phase-2.30 work.
- [~] 22-Nov-2026 capture/fidelity audit opened: determine
      whether job 13488131 supplies complete v2 decks and
      authoritative M/P schedules, then run a P2-off baseline.
      Exit criterion: validated replay inputs and fidelity
      verdict, or an exact missing-state blocker. No production
      numerical changes or P2-on live run authorized by this audit.
      Direct inspection finds decks for rays 81/82/84 through
      refresh 250 and KNAPPLY files. The earlier no-decks claim
      below was false and is superseded by this audit.
- [~] Part 3 deck capture: armed COLD-start run submitted
      13488131 (config dessart_m15_2p82d_VAL_knop_formal_rd.yaml,
      KN_RDUMP_LS=81/LS2=82/LS3=84, KN_ORD_P2DAMP=F, LTE_EST=T,
      warm_restart=false, pinned exe 358814d3/8a62e7ea at
      git=92284e61). No KNRDUMP decks exist anywhere else
      (phase-2.26 archive holds replay outputs only). Run
      FAILED at it250 after 22h: the wrapper's KNOP stall
      detector counted 149 consecutive GS-frozen LAMBDA
      iterations (limit KNOP_STALL_MAX=60) - the lagged field
      never matured from LTE (fg_rf ~1.04, rf_not_allconv).
      Not a crash; the non-convergence kill worked. Old
      warm-run 13473182 (it712, refresh residual 0.107 =
      stall floor) cancelled to make room. Deck capture still
      open - no decks produced before the kill.
- [x] OBS line treatment corrected: rd config inherited
      obs_line_treatment: sob (deposited Sobolev comparison)
      -> now blank per the 14-Nov directive (blanketed CMF
      spectrum). sob_diag OBS (13474709, completed) rendered
      SOB - consistent with that comparison cell. All six
      generated VAL_* cells except sob_jfnk/sob_tri carry
      sob; generated files, flagged not hand-edited.
- [~] 22-Nov-2026: bounded P2-off baseline audit submitted as
      interactive job 13557300; selected input validation at
      refreshes 1/10/250 precedes O3 production-object replay,
      twice per selected ray at refresh 250. First attempt
      13557295 lacked numpy; replacement explicitly loads
      python-waterboa/2024.06. No P2-on jobs submitted.
      Audit: CMFGEN_PHASE230_P2_PRODUCTION_OPERATOR_EVALUATION.md.
      Infrastructure tests: test_phase230_capture_audit.py 5/5
      PASS; black/isort applied; ruff unavailable in environment.
- [~] SOB-population BLANK render: isolated OBS-only job
      13557228 submitted via run_model_grid.py, currently pending
      on QOSGrpCpuLimit. Original SOB spectrum retained.
      Config: dessart_m15_2p82d_VAL_sob_diag_blank.yaml (parent
      project). Staged VADAT/RVTJ/all POP* byte-identical to source.
      CMF parameter preflight verifies BLANK GLOBAL_LINE/TRANS_*,
      DO_SOB_FLUX_OBS=F, NUM_ES=2. Spectrum validation still owed.
- [ ] Parts 4-10: baseline fidelity remains unproven; A/B replay,
      raw-defect evaluation, mechanism and branch analysis, and
      acceptance decision await validated inputs and baseline.
- [ ] Part 11 live validation: not started; gated on offline
      replay passing first.
- Note: commits made locally per change; NOT pushed per the
      phase-2.30 task restriction ("do not push without
      approval").

Living checklist of the non-monotonic CMFGEN development. Maintained
per AGENTS.md section 11: update it in the same edit as the change it
records, keep every tick honest, and cite evidence (audit doc, commit,
or job id) rather than intent.

Legend: `[x]` done and evidenced, `[~]` in progress or done with
documented limitations, `[ ]` outstanding. Blocked/failed items say so
explicitly.

## Fresh SOB + KNOP campaign (21-Nov-2026)

Run-area cleanup removed ~336 GB of stale model dirs (all KNOP_*, old
SOB/VAL/FRESH arms, phase227b captures); LUC baselines and inputs kept.
Resubmission on the fresh build (git=5236c37a) exposed two cold-start
blockers, both diagnosed deterministically (identical death point across
two submissions on different nodes):

- [x] `arnaud_rothenflug.dat` Ti III->IV line `3d_2De` renamed to
      `3p6_3d_2De` in all three cmfgen_base* copies and staged run dirs.
      READ_ARNAUD_ION_DATA does substring matching of the arnaud name
      against atom-file level names; the TkIV data carries the `3p6_`
      prefix so the line never matched, ION_POT 27.49 eV fell outside
      the (68.9,70.8) bracket, and the routine STOPped at 22 s - the
      "missing/stale RVTJ" signature. Offline simulation of the whole
      file+atom data shows this was the only line that STOPs under this
      species set; the next candidate (FeSIX->SEV 3s, 195.8 eV) lands
      inside its (184.9,224.0) bracket and continues, confirmed in the
      run that followed. Data fix, not code: shared base input file.
- [x] `se_ba_non_therm_v3.f90` guarded copies of the
      `RADIOACTIVE_DECAY_ENERGY` pointer on ASSOCIATED (commit
      57f34a7e). The pointer is allocated only under INC_RAD_DECAYS=T;
      with TRT_NON_TE=T (SOB recommended default) and zero isotope
      blocks in SN_HYDRO_DATA the routine was still called and
      segfaulted (SIGSEGV at line 76, batch exit 174). Not flag-gated:
      zero decay energy is the correct contribution.
- [x] sob_diag model run to convergence/spectrum: 13474708 stalled
      at it387 (counted correction pinned ~6% at depths 1-13, NG
      pulses correlated); warm restart with DO_NG=F converged at
      it409 (13567608, all gates MET) and rendered obs_fin_10
      (13567609) - rms 0.042 vs Luc across 7 bands. Cold-start
      acceptance arm 13569107 running (it76, demand ~38 at last
      poll); V13 discriminator 13572032 running (it33).
- [~] knop_formal running (13473182, warm restart off the converged
      tauto checkpoint, ~575 s/iteration, it>=686 at last check);
      obs follower 13473184 chained.

## Round-2 Phase-2 diagnostic implementation (21-Nov-2026)

Three instrumentation corrections + audits; physical solver,
thresholds and production behaviour unchanged except the replica
correction itself. Head: `a2bbc85b` on main.

- [x] Patch A - nine-gate ITER_SUM replica: shared `CONV_EXIT_EVAL`
      in `lam_handoff_mod.f90` (512-combination truth table, ALL=AND,
      NPASS=popcount, `lam-field-not-ready` label, denominator 9).
      solve_for_pops replica calls the evaluator with the LIVE
      `LAM_FIELD_NOT_READY`; real `CONV_EXIT_MET` untouched.
      Commits bb71cfba..4b0c09d4; test_conv_exit 46 checks green.
- [x] Patch B - margin reporter: `CONV_MARGIN_EVAL` separates
      active-binding demand from SAT/small/trace-released demand;
      reporter prints `worst active margin` vs `largest released
      demand` with the release leg. Commits 61c11c5e..4b0c09d4.
- [x] Patch C - symmetric FG sweep drift `FG_GS/SW_MAXSYM` +
      argmax (depth,ray,freq,P/M) exported via `GET_FG_SW_SYM` into
      GS_DRIFT_DIAG and the compact fg row; diagnostic-only, never
      a gate. Bad-cell identity: per-(ray,channel) sparse packed-id
      lists (fixed-stride (ray,freq,point,ch)), serial merge diff
      with Jaccard/new/gone/untracked, `KN_BADID` VADAT key
      (default F). Commits 40c12cb7..a2bbc85b; test_fg_sw_sym 60
      checks green incl. end-to-end refresh merges == FG_GS_NBAD on
      both sweep paths; test_knop_ray 49 + test_conv_exit 46 +
      test_lam_handoff 51 green.
- [x] Assumed-shape interface bug fixed in `KN_BADID_MERGE_V1`
      (NRAY arg, explicit-shape dummies) - would have segfaulted in
      production; found by the driver test before any model run.
- [x] Phase-7 demand audit + Phase-6 block-GS assessment:
      `docs/audits/CMFGEN_BA_DEMAND_BLOCKGS_AUDIT.md` - demand
      carriers are C V/C VI quasi-dead rows (n~1e-62..1e-96) at
      depths 1-22; block-GS NOT justified (orbit is local
      border-collision, not cross-frequency lag).
- [~] Build: all exes built (cmfgen_dev 13:54, suite 91/100 in the
      build log); 9 suite fails include the pre-existing hung
      `drv_obs_frame_norm_ref` (poll on stdin, killed) - unrelated
      to this phase; individual checks for A/B/C all green.
- [ ] Model-run verification through scripts/run_model_grid.py for
      the new diagnostics (pending - needs an arm; the KN_BADID /
      GS_DRIFT outputs then land in OUTGEN).

## Phase 2.28 limit-cycle identification (21-Nov-2026)

Audit: `docs/audits/CMFGEN_PHASE228_LIMIT_CYCLE_AUDIT.md`.

- [x] Ray-81 geometry resolved: all-expanding (GAM>0 at all 5 points,
      min T2=0.85), zero compressing cells; multi-resonance and
      direction-handling hypotheses refuted for this ray.
- [x] Replay-history per-cell analysis: the L081 stall is an exact
      period-2 limit cycle (bit-identical two-state alternation,
      54K+ cells on the boundary row, omega-independent per the
      phase227 sweep). Ray 82 same class; ray 84 carries the same
      orbit at the underflow floor and certifies clean - certified
      class tracks where the orbit's amplitude meets the
      significance floor, not map quality.
- [x] Orbit localised to the lagged corner/boundary feedback
      (lagged SPW preview into KN_ORD_DS, endpoint DS(NI), relaxed
      IM_NI copy) on high-albedo floored cells; micro-channel
      inside that loop still open (PARTIAL).
- [ ] C1 evaluation pending: KN_ORD_P2DAMP on a diagnostic arm -
      the replay histories are the missing per-cell evidence the
      flag was waiting for. No solver change authorized.

## Phase 2.29A P2 state-management fix (01-Dec-2026)

Audit: `docs/audits/CMFGEN_PHASE229A_P2_STATE_FIX.md`.
KN_ORD_P2DAMP default unchanged (F).

- [x] Stale-state defect fixed at all three march paths: the P2
      tracker ran only inside the significance gate while the
      midpoint commit ran unconditionally; insignificant cells now
      clear P2_SGN/P2_FLIP/P2_DMP on the ELSE branch and take the
      ordinary relaxed update. Triplicated inline logic replaced
      by shared module procedures KN_P2_TRACK/KN_P2_CLEAR/
      KN_P2_MID (module scope for unit-test access).
- [x] Midpoint underflow fixed: SQRT(A*B) -> SQRT(A)*SQRT(B);
      the product form collapsed to zero below ~1E-154 endpoints
      while floor cells legitimately reach ~1E-283.
- [x] RES_INT/RES_BND diagnostic split added (boundary rows =
      inward I=1 zero row + outward NI/IM_NI row); KNPASS print
      under PRINT_KNOP_DEBUG. Diagnostic only; RES criterion
      unchanged.
- [x] Tests: new test_knop_p2_state.py 19/19 PASS (state
      machine, all three call sites, midpoint range); existing
      test_knop_ray.py 49/49 PASS incl. bit-identical
      repeatability and 1v4-thread determinism - the off-state
      path is bitwise unchanged.
- [~] Replay evaluation on immutable L081/L082/L084 histories:
      orbits confirmed at interior cells near l~1517-1523
      (rel ~4e-2 L081, ~8.5e-2 L082); simulated tracker would
      have armed. P2-on replay of the production operator not
      possible - the KNRDUMP input decks no longer exist;
      R(X) unevaluated. Arm KN_RDUMP_LS on the next run.

## Phase 2.27B isolated preparation (21-Nov-2026)

- [x] Rechecked both original captures: JobHeldUser. No release, cancel,
      replacement or new submission. Installed executable hashes unchanged.
- [x] Preserved dirty source, tracked diff, replay source, audit documents
      and scheduler references under external phase227b_20261121_0604.
- [~] Isolated build_v4 normal/debug candidates compile and link; same new
      solver object is linked into candidate executable, replay and IO harness.
      Other archived objects retain production flags. Not runtime-tested.
- [~] omf now uses one uninterrupted ray call with strict two-channel omega
      injection. State inventory and stage semantics are source-reviewed,
      not a demonstrated trajectory reproduction.
- [~] Separate snapshot/manifest/runner templates prepared with hash checks
      and deliberate execution block pending input/restart read-set review.
      Read-only/hash-sealed is not filesystem-enforced immutability.
- [ ] BLOCKED: actual Fortran IO tests, ray-84 control, O0/O3 field comparison
      and 81/82 identical-stage metrics require an approved isolated allocation.
      No submission is authorized in this task.
- [ ] AWAITING CAPTURE: actual omega histories, complete 84-slot census and
      per-cell transition evidence; existing row format lacks some fields.
- [ ] Review failure/heal stage coverage and full input closure before any
      separate scheduler approval. Phase remains PARTIAL; AA/P2 remain off.
      Evidence: knop_phase227/{capture_isolation_design,replay_state_inventory,
      replay_lifecycle_comparison,runtime_test_results}.md and completion_status.json.

## Phase 2.27 continuation (21-Nov-2026)

- [~] PARTIAL: omega-history replay and per-refresh ray census unavailable.
      Captures 13464217/13464311 were pending AssocMaxJobsLimit; held with
      explicit user approval after finding shared mutable run directories.
- [~] Diagnostic source only: remove census format type mismatch, increase
      omega precision and harden omf history validation. ifort syntax checks
      pass; no executable replacement or runtime regression claimed.
- [ ] Isolated, hash-pinned captures and allocated replay regressions;
      uninterrupted sweep-history injection, both-channel percentiles,
      authoritative ray timelines and cell-level attribution still required.
      See docs/audits/CMFGEN_PHASE227_ADAPTIVE_OMEGA_AND_RAY_CENSUS.md.

## 1. Roadmap (AGENTS.md 1.1)

- [x] `DO_NONMONO_SOB` on a non-monotonic model producing a spectrum
      (`m15_2p82d_SOB_l14_lte_np` converged 100% at it718;
      `m15_2p82d_FRESH_sob` running, it346+ at last check).
- [~] Non-mono SOB on an interacting (IIn-like) model - the m15
      interactings run but residual non-convergence is still under
      measurement (ES sweep studies, `CMF_FLUX_P8_2_ES_CONVERGENCE_AUDIT.md`).
- [~] Knop arbitrary-velocity solver as a plug-in - solver machinery
      complete and regression-tested; production viability blocked by
      the drift/non-contraction issue that phases 5.x are resolving
      (see section 3).

## 2. Solver pathway status

- [x] SOB production route `DO_NONMONO_SOB` + `USE_LAM_ES=T` (V13
      formal J) - validated default for new configs.
- [ ] SOB moment route `USE_LAM_ES=F` - FAIL on m15 (job 12316397:
      lagged map never contracts, ~183K persistently frozen cells).
      Kept documented, not selected. Stability audit 24-Nov-2026:
      per-solve tridiagonal operator VERIFIED sound; lagged map
      stiff by construction and stabilized by counted mechanisms
      incl. two equation-modifying bounds; defect found -
      refresh/sweep NMID_ON_HMID regularization mismatch holds a
      residual floor at sign-transition cells. Audit:
      docs/audits/MOM_J_CMF_NONMONO_STABILITY_AUDIT_2026-11-24.md.
      Retest arm (24-Nov-2026, user-directed):
      configs/dessart_m15_2p82d_VAL_sob_lamesF.yaml - warm restart
      from the converged m15_2p82d_VAL_sob_jfnk_eps1em3 checkpoint
      on the post-audit solver (commits 910c4667, 51542bbf: refresh
      regularization now matches the sweep, JUP finite guard, sweep
      off-grid census). DO_JFNK=F (admission needs USE_LAM_ES=T).
      Question: does the lagged J map contract from converged
      populations on the m15 ejecta? Evidence: MOM_GS drift census,
      refresh convergence (MOM_RF_CONV/DCHG), quarantine count,
      LAGB/XMNEG census. Exit: converges (map contractive, route
      back in contention), stalls again (same bounded class), or
      detonates (new defect). Single arm, no chain.
      VERDICT (24-Nov-2026, job 13646071, killed at it~688 after
      22 it): STALLED - same bounded non-contraction class. The
      compressing branch is the whole defect: expanding-direction
      refresh passes contract cleanly (residual ~1e4 -> 1-3, a few
      hundred cells above 5%), compressing-direction passes stall
      at residual ~1e4 with ~10K+ cells >5% every iteration.
      Retained store oscillates between best-pass candidates
      differing by 1e6-1e8 (bounded limit cycle, no fixed point).
      Containment clean: 0 negative J_NEXT, 0 quarantined, 0 solve
      failures; F4 census live (off-grid reads 0 on sweep, 184
      inside one refresh). Structural deadlock: conv-exit held on
      lam-field-not-ready (MOM_J_GS_OK can never fire while the
      refresh stalls) + ex_maxch>=eps (demand 1e5-1e6 vs eps 1e-3).
      Stall population = the audit's predicted marginal cells:
      ~10K GAM-floor (|1+sigma*K| resonance floored at 1e-2),
      ~96K eps-limit hits/pass, ~8.6K neg-RHS rows, LAG_BND bound
      scaling 170K-1M cells/sweep (smallest lag 0.0). Fixes worked
      as intended (map-mismatch floor removed, honest counting);
      the rho~1 compressing lagged coupling is not fixable by
      relaxation - needs an implicit/direct frequency-chain solve.
      OBS 13646072 cancelled on the dead dependency.
      FIX (24-Nov-2026): `MOM_NUSOLV` (VADAT key, default F)
      implements that direct solve inside the refresh driver. Per
      depth the rows' upwind-J coupling is a functional graph in
      frequency (each row owns one edge: k-1 expanding, k+MOM_NM_J
      compressing) - chains plus rings at the sign flips that pin
      the map at rho~1. `MOM_NU_CHAIN_SOLVE` resolves the captured
      equations exactly in O(NCF) per depth after each pass
      (Baron & Hauschildt 2004, A&A 427, 987: the wavelength-chain
      boundary-value structure wants a direct solve, not Jacobi
      iteration). Capture folds LAG_SCL/EPS_SCL/boundary sign into
      the coefficients so the composite map's fixed point is the
      swept system's; NU_FRESH freezes active-set-skipped columns;
      write-back reuses the sweep's damped non-negativity update;
      NU_CHG merges into PASS_CHG for the exit check. Unit
      verdict: test_mom_nu 32/32 exact vs dense solves (chains,
      2-cycles, 4-ring, mixed graph, GS-depth ordering, frozen/
      rejected paths); nurf arm lands the NUSOLV=T fixed point
      2.2% off the NUSOLV=F one on comp1; nuw admissible on compw;
      VADAT path verified both stores (test_vadat_keys +2 arms).
      Open: whether it contracts the real m15 stall - needs an
      arm run (bounded investigation, declare before submit).
- [~] `DO_NONMONO_KNOP` ray solver - retained for reference/diagnostic
      comparisons; not selected for new production configs (user
      directive 11-Nov-2026). Under active defect-resolution (sec 3).
- [x] JFNK population solve (`DO_JFNK`, admission gated on
      `USE_LAM_ES=T`) - implemented, test batteries green; coupled
      evals use a frozen lagged store (documented caveat).
- [x] `DO_ES_ALI` accelerated Lambda iteration - ~30 warm-start
      iterations to the control fixed point (X2/X3 arms);
      `DO_ES_ALI_BA` stays opt-in.

## 3. KNOP P-sweep defect cycle (phases 5.3K-5.3S, Nov-2026)

Audit: `docs/audits/PHASE_5_3_P_DIRECTION_AUDIT.md`.

- [x] Demonstrated the `HAS_CP=F` false-convergence defect: both
      sweeps measured a nonzero pass residual then forced
      `RES=0`/`IXIT=0`, committed as converged (`ray/ixit=0 npas=2
      resid=0` on every refresh, frz82 baseline).
- [x] Mechanism: the ordered stencil keeps lagged C-corner terms
      (Knop II Eqs 31-32) plus ring-aliased slope reads at both ray
      ends, so one pass is never exact; ray-82 P-dir oscillates as a
      two-state orbit (~{3.18e-8,1.90e-8}) across passes/refreshes.
- [x] Manufactured ring-faithful stencil test: unique seed-
      independent fixed point, roundoff equation residual; Aitken
      and under-relaxation break a genuine two-cycle, shallow
      Anderson does not.
- [x] 5.3R correction committed (`0ef42d05`): forced-clean exits
      deleted, `IXIT=0` retired, all success paths require measured
      `RES<KN_RAY_TOL`; debug exe `c3d8b584`, production exe
      `f52d8125` untouched.
- [~] Production validation of the fix - `p53q_probe82` (13157367)
      delivered the first pass-3+ history: ray-82 P-dir transient
      resolves by ~pass 8, then plateaus ~0.05 (stall exit ixit=3,
      would-be capped verdict). The two-cycle was the early
      transient, not a persistent orbit - audit revised. `p53s_fix82`
      (13158010, corrected exe) still in setup; will confirm the
      same exit path on all rays plus honest `fg_rf`/`KN_HOT_CNT`.
- [~] Ray-83 floor-cell behaviour under measured exits (jump veto vs
      convergence classification separation): under measured exits the
      ray neither converges nor jump-vetoes - it exits capped
      (ixit=3) on every refresh while its residual holds a period-2
      cycle (6.894E-01/6.961E-01, B-arm evidence 21-Nov-2026, 113+
      refreshes, ~88K bad cells/sweep).
- [ ] Corner-mode acceleration decision (guarded Aitken / OMEGA_RAY /
      safeguarded Anderson) - deferred until production pass history
      exists; do not enable an accelerator without it.
- [ ] Whether rays 82/83 share a mechanism or need separate handling.
- [x] Phase-2.21 A/B and repeatability arms (21-Nov-2026): pre-flag
      frozen exe 93364cac (`KNOP_B1_PREFLAG`, 13427111) vs current
      58a699d0 `FG_RF_ALL_CONV` + NI guard (`KNOP_B2_POSTFLAG`,
      13427115) produced the same pinned ray-83 residual cycle -
      the flag change does not alter this ray's stalled lagged map.
      A/A pair 13414391/13414398 ran near-identical trajectories:
      the stall is deterministic, not a seed or race artifact.  All
      four arms, their pending chain followers, and the held
      `knop_v2cap` jobs (13464311/13464217) were cancelled after the
      assessment; the arms were diagnostic-only under the 11-Nov
      KNOP scoping and their questions are answered by the pinned
      cycles themselves.
- [x] Frozen-executable experiment hygiene (the p82 drift root cause):
      `scripts/cmf_run_manifest.py` pins the resolved exes into
      `exe_frozen/<sha16>_<sha16>/` shaped as a cmfdist; submit_cmfgen.sh
      repoints CMFDIST at it, writes `manifest.json`, and the generated
      slurm re-verifies sha256 immediately before exec (mismatch aborts
      before batch.sh).  Post-run segments record exit code, termination
      class (normal / stall / signal kill) and output hashes.  Opt-out:
      `cmfgen.pin_executables: false`.  SN2024pba repo commit f0902fb9;
      tests `scripts/test_cmf_run_manifest.py` (4/4).  fix82 (13158010)
      still owes the all-rays honest-exit evidence - killed by SIGXCPU
      (interactive CPU cap) mid-iteration-1; needs a general-partition
      resubmit.
- [x] KNOP spectrum-output gate (P8.3): `obs/knop_spec_gate.f`
      classifies GS_STATE + stall marker + manifest before the
      OBSFLUX/OBSFRAME writes; non-converged KNOP models stop under
      `DO_KNOP_SPEC_GATE=T` (default) or render a marked diagnostic
      under the F opt-out.  `unit_tests/python/test_knop_spec_gate.py`
      (9 cases, production object).

## 4. Diagnostics and instrumentation (all default-off, verified)

- [x] `GS_DRIFT_DIAG` per-iteration drift reporting (`test_gs_drift_diag`).
- [x] `KN_WATCH_CELL`/ray lifecycle prints.
- [x] `KN_NOCP_PROBE`/`KN_NOP_*` pass telemetry; `KN_PSWEEP_*`
      complete-map capture incl. binary X/Y snapshots
      (`kn_psweep_source_test.py`, `kn_psweep_check.py`).
- [x] `KN_FREEZE_COEF`/`KN_MAT_FREEZE` coefficient-freeze arms with
      arm-event fingerprinting.
- [x] `ES_RESID_DIAG` OBS-side ES residual writer.
- [x] `DAMP_FJ`/`J_RESCUE` continuum-J rescue controls
      (`COMP_JCONT_V4.INC` etc.) - defaults preserve stock behaviour.
- [x] `ITER_HISTORY`/`ITER_LEDGER` machine-readable iteration rows.

## 5. Test and validation tiers

- [x] Static/unit batteries: knop_ray 49, knop_lagged_map 25/26
      (sdab pinned-fail, documented defect), lagged_field suite,
      knop_operator 32, knop_spv 93, knop_p2_state 19,
      knop_replay_p2 29, knop_spec_gate 19, knop_endpoint 6,
      mom_solve 143 (incl. nurf/nuw MOM_NUSOLV arms), mom_nu 32,
      mom_gs 34, mom_aa_state 62, mom_sweep_table 8,
      es_resid 20, jfnk opint 36 + regressions,
      obs_frame 11+15, ng_abstain 17, lam_handoff 51, gs_drift_diag,
      iter_conv_gate 13 - green at last run (22-Nov KNOP sweep).
      test_es_conv remains blocked on the known obs/ archive gap.
- [x] Manufactured ordered-stencil + exit-policy + period-2
      classifier (`kn_ccorner_test.py`).
- [x] Stale-archive gate maintained (libsubs.a member rebuilt with
      the standard ifort flags; production exe untouched).
- [~] Pre-existing gaps: `test_es_conv`/`test_es_resid` gate on obs/
      sources never archived; synthetic numeric drivers beyond the
      manufactured stencil remain listed in FG_FIELD_DRIFT_AUDIT
      Pillar 7.
- [ ] sbatch model-level validation of the 5.3R fix (dominant owed
      item): fix82 vs frz82 A/B under identical VADAT/checkpoint,
      then a KNOP model through `run_model_grid.py`.
- [ ] 1-thread vs multi-thread refresh reproducibility check on the
      corrected build.

## 6. Reporting consequences to verify on first post-fix runs

- [x] `fg_rf`/`KN_RF_RESID` now include previously suppressed
      HAS_CP=F rays (ray 82 ~0.6) - the elevated values are genuine
      non-contraction, not a reporting artifact: probe ray 83 holds a
      period-2 6.894E-01/6.961E-01 residual cycle across 113+
      refreshes on the post-flag build (21-Nov-2026, section 3).
- [ ] `KN_HOT_CNT` should accumulate the L~17881 corner column.
- [ ] Converged/rejected/capped census lines change class for
      corner-oscillating rays - record the new mix in the audit.

## 7. Hygiene

- [x] Vault notes updated for removed exit class (`Knop formal
      solver.md`, gitignored but maintained on disk).
- [x] TODO/audit documents current (`docs/plans/FIELD_DRIFT_TODO.md`,
      `docs/audits/PHASE_5_3_P_DIRECTION_AUDIT.md` 5600+ lines).
- [x] All work committed and pushed to `origin/main`.
