# CHECKLIST - cmfgen_jfnk development ledger

- [x] Source tree moved under cmfgen/ (26-Nov-2026). All upstream
      mirror dirs, the fork source dirs (plane/, pgplot_stubs/,
      obs_template/) and the upstream root build files now live in
      cmfgen/; the repo root keeps the dev tooling (unit_tests/,
      tests/, docs/), build.sh drivers and untracked outputs.
      cmfgen/Makefile_definitions redirects INSTALL_DIR one level
      up, so exe/lib/mod still land at the repo root and every
      submitted slurm script keeps its exe path. build.sh,
      unit_tests/common.py (new SRC anchor), the test harnesses,
      AGENTS.md, README.md and the submit_cmfgen.sh CMFDIST refs
      (com/restart_model.sh, obs_template) were updated to match.
      The ptmp2 worktree keeps its old layout until it pulls; a
      rebuild there writes into the new cmfgen/-relative layout.
      Follow-up same day: the generated slurms also resolve
      $cmfdist/com/out2in.sh at epilogue time, a ref missed in the
      first pass - it produced exit 127 on every CMFGEN teardown
      (Ni075_E13_m06 job 13722141, knop_jfnk 13699394). The
      cmfgen_base template, run_model_grid.py and ~148 generated
      slurms/submit copies were repointed to cmfgen/com/out2in.sh
      (project repo, uncommitted).
      Verdicts landed alongside: the Ni075_E13_m06 40 d arm NaN'd
      again at it35 with the same 193-bad STEQ signature on the
      synced deck (TRT_NON_TE=T, LAM_EXIT=5.0D1, LAM_N_STABLE=3) -
      the deck hypothesis is rejected; the failure is the open
      BLKBAND/Jacobian defect class. knop_jfnk (13699394) stopped
      non-convergent at it100 on GS_REL_NOPROG_LIM - the
      release/refreeze cycle never matured the lagged field; its
      OBS chain (13699395) was cancelled. knop_formal (13699389)
      still runs; its OBS (13699390) will need resubmission after
      teardown.
- [x] DO_JFNK + FIX_X_DPTH inner-T pin: fail closed (25-Nov-2026).
      Question: does a JFNK run silently ignore the FIX_INB_T +
      FIX_X_DPTH thin-release T pin? Evidence plan: audit which
      T-mask channels each solver reads. find showed JFNK tracks
      RD_FIX_T, JFNK_FIXT and the VARFIXT range
      (JFNK_FIX_T_D_ST/END) but nothing on the JFNK route reads
      FIX_LST_X_DPTHS - the pin lives only in fixpop_in_ba_v3,
      which zeroes those T rows of the BA system each call. The
      JFNK unknown vector would keep those rows live while the BA
      matrix holds them fixed: two solvers carrying different
      equation sets. Verdict FIX landed: check_param_consistency.f
      now rejects DO_JFNK=T + FIX_INB_T=T + FIX_X_DPTH>0 unless
      FIX_T=T (all-depth pin is honoured by both routes);
      FIX_T_AUTO stays legal (JFNK_TAUTO mirrors the VARFIXT
      range). The subroutine picks up USE MOD_CMFGEN for the FIX
      flags and the Makefile dependency is updated. prepare_VADAT
      already warned on the pairing; this makes the same contract
      fail closed at the Fortran gate. Tests: test_vadat_keys
      gains three cases - the combination aborts, the pin without
      JFNK parses legally, and the FIX_T=T exemption passes (105
      checks total). One retained probe deck
      (dessart_m15_2p82d_SOB_l14_freet_jfnk_tri.yaml) trips the
      stop by design; its header notes the abort.
- [x] Commit gate / W16 jump guard: per-column scale (25-Nov-2026).
      Question: can dim frequency columns grow by orders under a
      global-max jump scale, as the FG_JUMP_MAX port claims the
      comment intends? Evidence plan: manufactured 200-column
      12-decade J spectrum - blow-ups of 1e2..1e8x in the 50
      dimmest columns flag 0 cells under the global rule, 3000
      under a per-column one. Verdict FIX landed: both callers
      (W16 best-pass guard, MOM_RF_COMMIT_GATE) now share
      MOM_JUMP_CENSUS, which keys the 50x test on the column's own
      store peak floored at 1E-300 and prices floor columns
      (store peak <= 1E-28) against the nearest committed column's
      peak, exempt only on a true cold start - the MOM analogue of
      FG's sibling-reference rule (per-ray siblings do not exist
      on the moment grid; nearest committed frequency columns do).
      Tests: drv_mom_gate dimjmp (1e6x dim blow-up rejected),
      dimok (3x reshaping accepted), flrseed (5E-7 fill inside
      50x of a 1.8E-6 reference accepted), flrjmp (1E-2 seed
      rejected), coldstr (all-floor store exempt); the 8 prior
      cases pass unchanged, and KEEP_BEST verdicts are identical
      on the single-scale fixture because per-column and global
      scales coincide there.
- [x] Converged exit ignores f-iteration cap exits (25-Nov-2026).
      Question: can the run print "Convergence criteria met" on an
      iteration where some frequencies left the f/e-scattering
      loop at the pass cap, leaving an unconverged ES source in
      the rates? Evidence plan: FCONV_CAP_CNT was printed and
      reset each iteration but never read by CONV_EXIT_MET.
      Verdict FIX landed: the census is saved to
      FCONV_CAP_IT_CNT/T1 before the reset and joins the gate as
      a tenth conjunct in both CMFGEN_SUB copies plus the
      SOLVE_FOR_POPS ITER_SUM replica (CONV_EXIT_EVAL carries
      CAP_CNT, percent denominator now 10). A held exit names
      'f-iteration cap exits outstanding' in CONV_EXIT_WHY and
      prints a distinct "Converged exit deferred" line; on
      LST_ITERATION the budget-exhausted message reports the hold.
      Suppressed iteration-1 warnings stay inside the census (the
      gate reads the physical event, not the print policy).
      Deadlock analysis: if a frequency keeps capping the exit
      holds until the iteration budget runs out - stock CMFGEN
      behaviour is exactly this (the run ends unconverged and
      writes unconverged state); no threshold or fallback added.
      Tests: drv_conv_exit truth table extended to 1024
      combinations, named case G_cap pins cap>0 blocking with
      npass=9/10 and the 'f-cap-exits' label; test_conv_exit pins
      the save-before-reset ordering and the gate reading
      FCONV_CAP_IT_CNT rather than the suppressed subset.
- [x] Investigation: sign-blind adaptive JNEXT_RELAX (25-Nov-2026).
      Question: does rho = |r_k|/|r_{k-1}| feeding omega =
      1/(1+rho) miskey relaxation on monotone-contracting maps?
      Evidence plan: Python model of a scalar linear map and an
      alternating two-map schedule, lambda in {-0.9,-0.5,+0.5,
      +0.9,+0.99}, current rule vs undamped vs a signed estimate.
      Verdict CONFIRMED: at lambda=+0.9 the adaptive rule parks
      omega at ~0.51 and contracts at 0.949/pass vs 0.9 undamped
      (hit the 400-iter cap vs 259 undamped); at lambda=+0.5,
      79 vs 39 iters; alternating (0.5,0.9) 132 vs 69. For
      lambda<0 the rule is right (28 vs 277 iters at -0.9) - the
      formula is optimal exactly for lambda=-rho. The pass-k+1
      input IS AA-mixed (inner AA at every even pass), so rho
      measures the accelerated trajectory, not the bare map. A
      signed estimator lambda_hat = <d_k,d_prev>/<d_prev,d_prev>
      on same-direction update vectors is prototyped behind
      module flag MOM_OMEGA_SGN (default F, no VADAT key -
      [MOM_OMEGA_SGN] proposed in the PR): RES_SGN files the
      signed same-direction ratio per pass (stride 1 same
      direction, stride 2 under ALT - the LAMC pair), and
      MOM_SAME_DIR_SGN medians it into omega* = 1/(1-lambda_hat).
      Model check: signed rule matches or beats the current rule
      on every scalar case and never worsens oscillation. The FG
      sibling (FG_ADAPT_OMEGA_V2) has the same defect class but a
      different machinery (per-ray pooled RHOM/RHOP, P90
      statistics); port deferred until the MOM path is exercised
      on a real run.
- [x] Investigation: bulk census vetoes legitimate AA steps
      (25-Nov-2026). Question: does the N5F/N50F-vs-N5P/N50P churn
      test reject the large-step candidates outer AA exists to
      produce? Evidence plan: 200-cell manufactured rho=0.95 map,
      3-term Aitken every third pass, gate replica applying the
      exact counting rule; repo grep for "commit gate rejected".
      Verdict CONFIRMED, rule unchanged (proposal only per task):
      the Aitken step reaches residual ~1e-13 where the pre-AA
      field sits at ~0.6, yet the gate rejected all 36 offerings
      because the candidate churns more cells vs the retained
      store than the pre-AA field - it prices movement, not
      residual reduction. Cost in the model: 3 iters without the
      gate vs 111 with. No "commit gate rejected" lines exist in
      tests/ or docs/ logs, so the veto has not been observed on
      a real run. Manufactured case aalegit pins the behaviour.
      Proposal (not implemented): accept the candidate on its
      map residual r(F)=|G(F)-F| instead of store churn - one
      extra map evaluation per gate call (a refresh already costs
      KN_M_MAXPASS of them).
- [x] JFNK Phase 5: per-depth LU backward-error certificate +
      LSQ fallback (25-Nov-2026). Question: can the per-depth
      factorization prove fitness instead of trusting dgetrf
      info=0, and can a depth that defeats the LU ridge ladder
      be rescued instead of dropping to Jacobi? Evidence plan:
      opt-in [JFNK_PC_CERT] solves M_s z = M_s v on a fixed
      probe and rejects when the componentwise relative
      residual exceeds JFNK_PC_CERT_TOL (default 1e-8);
      [JFNK_PC_LSQFB] adds one attempt slot past the ladder top
      that factorizes B^T B + lam^2 I (SPD for lam>0) and is
      accepted only by the same certificate; JFNK_DEPTH_LSQ(K)
      records the choice so the apply dispatches DPOTRS there.
      Both default F. Verdict FIX landed, measured on
      jfnk_core: test_pc_quality Q7 cert accepts a good block
      and Q7b rejects at impossible tol; Q8 a block kept
      exactly singular by identical row+diagonal fails the one
      unridged LU attempt (ridge_floor=0, itry_cap=1), the LSQ
      slot factors it and reports lsq_used=T at itry=itry_max+1;
      Q8b same block with lsqfb off reports lu_ok=F. Harness
      note: the new OPTIONAL lsq_used dummy crashed
      implicit-interface test callers on the present() check -
      test call sites now pass a scratch logical explicitly
      (external subroutines), test_pc_quality keeps an explicit
      interface block. Batteries: jfnk_core 9/9 incl new Q7-Q8,
      opint 36/36, regressions 96 PASS incl new pins, se_solver
      unchanged (same documented XFAILs). Caveat recorded:
      dgetrf is backward stable so the certificate rarely
      rejects LU in practice - its real gate is the LSQ slot,
      where potrf failure or a bad residual is the only
      acceptance test. Remaining owed: Issue-A block-local
      derivative augmentation (folded into this phase scope per
      plan - not started); thread-count reproducibility (the
      se_driver surrogate carries no OpenMP; residual-side
      determinism is by construction - disjoint per-K writes,
      integer reduction - and only a full-model state can
      exercise it).
- [x] JFNK Phase 4: deterministic JVP step (25-Nov-2026).
      Question: can the FD step be made a pure function of
      (jfnk_x, v, POPS_BASE) so the effective operator is a
      fixed map during each Arnoldi cycle? Evidence plan:
      extend the a-priori caps to cover every crossing
      JFNK_N_NEW_CLIP counts, delete the eval-result retry,
      re-measure additivity/homogeneity/repeatability on the
      surrogate batteries. Verdict FIX landed: the deleted
      NRETRY loop shrank h on the eval RESULT - deterministic
      in (x, v) but state- and path-dependent across calls.
      All four counted crossings are now capped before the
      single trial eval: the +/-50 bound on T is now applied
      in both unknown modes (the residual clamps |x_T| at 50
      unconditionally; previously only capped under JFNK_LOGX)
      and the non-log pop/ED 1e-100 floor is capped at half
      the remaining distance. An unprevented crossing (a
      boundary the caps do not model) keeps the FD as-is and
      is counted in JFNK_N_JVP_HRETRY; JFNK_CLIP_RETRY is
      parsed but unused (deck compatibility). Measured after:
      opint 36/36 with retr=0, caps=0, clip=0, additivity
      9e-8..1.3e-7 (FD floor), homogeneity exactly 0, bitwise
      repeatability; se battery unchanged (mild T3 + T21 x3
      are the documented known-diagnostic FAILs, values
      bit-identical to pre-change); regressions 87 PASS incl
      new pins for the deterministic step and a-priori caps.
      T9 null-direction floor (rel err 3.5 at 1e-8) is
      unchanged and inherent to forward FD - the
      true-residual gate remains the mitigation, by design.
      Note: the Phase-2 PTC x~=0 test (test_ptc_coords.f90)
      landed in the prior commit batch; the entry below still
      lists it as owed. Remaining owed: Phase-5 per-depth LU
      backward-error certificate + TSVD/LSQ fallback; Issue-A
      block-local augmentation; thread-count reproducibility.
- [x] JFNK independent audit + phase-1 fixes (25-Nov-2026).
      Audit verdicts: (a) BA-vs-FD Jacobian gap CONFIRMED at
      source - JFNK_SYNC_ATM_FROM_POPS recomputes POPION,
      COMP_LEV_DIS_BLK and SUP_TO_FULL_V4 at perturbed pops
      while BA carries no derivatives of the aggregate
      dissolution factors (interface l.775 comment admits);
      all channels are depth-block-local and analytic.
      (b) JVP nonlinearity CONFIRMED measured: T9 near-null
      direction rel err 3.5 at eps=1e-8 (hard surrogate);
      additivity defect ~6x FD error; homogeneity exact.
      (c) Production LU preconditioner CONFIRMED deficient:
      GMRES residual grows 5.18e37 -> 8.91e37 under LU on
      hard surrogate while TSVD/PTC-LU/LSQ converge ~1e-12;
      per-depth LU backward err 2.4e-2 (author XFAIL);
      dgetrf info=0 does not imply a usable operator.
      (d) Equation-structure freeze, PTC coordinate fold
      (sigma*exp(-x)/dtau), residual purity (byte-identical
      repeats), globalization norms all VERIFIED correct.
      (e) MAXCH defect CONFIRMED: factor folds S_PTC into
      JFNK_BA_DIAG (in place, persistent via BA_A) so
      jfnk_maxch_est divides by the regularized diagonal.
      FIX landed this commit: JFNK_BA_DIAG_RAW +
      JFNK_BA_DIAG_RAW_READY added (int_mod), filled from
      JFNK_C_MAT before the PTC fold in BUILD_BA_RAD and
      FACTOR rebuild branch; MAXCH, RSCALE=DIAG,
      RSCALE=DIAG_RDG and the newton-diagnostic BA_diag
      column now read RAW-preferred. Schur S3 defect fixed:
      jfnk_blk_lsq_solve flags nonfinite DPOTRS results
      (info_out=np+1), schur_form rejects nonfinite S before
      DGETRF, schur_apply flags nonfinite v (S3 apply now
      ierr=5 instead of propagating garbage).
      Test-infra repairs (kv-refactor fallout, pre-existing):
      iter_log_fmt_mod.o companion added to
      unit_tests/common.py _PROD_COMPANION_SOURCES for
      jfnk_solver_v1.f90 + jfnk_op_probe_v1.f90 and to all
      tests/jfnk_core link rules incl test_trace_state;
      stale literal patterns in test_jfnk_regressions
      ('BA src it=unknown', raw-direction A/ES format) updated
      to the ILF_KVB forms; parse_opprobe absorbs ilf_kvb
      continuation lines. Tests: jfnk_core make test all 9
      binaries PASS (test_trace_state now builds); opint
      36/36; regressions all PASS incl new Issue-E/S3 checks.
      Owed: Phase-2 PTC x~=0 regression test; per-depth LU
      backward-error certificate + TSVD/LSQ fallback (Phase
      5); deterministic norm-scaled JVP step (Phase 4);
      block-local augmented preconditioner for the Issue-A
      derivative channels (Phase 5/6); thread-count
      reproducibility.
- [x] Production-route consolidation audit (25-Nov-2026, work
      order: moment solver preserved but no longer the production
      target; production = DO_NONMONO_SOB + USE_LAM_ES=T).
      Audit verdicts: J-field dispatch clean (V13 RJ=TC writes
      J/rates under LAM_ES=T); ALO accessor pairing consistent
      (GET_FG13_ALO_DIAG vs GET_MOM_ALO_DIAG); ES-ALI gating
      correct incl. DO_ES_ALI_BA hard-STOP; JFNK admission on
      genuine field maturity (jfnk_cmt ledger honest);
      BLKGMRES isolation holds (route gate, KN_M_READY
      self-guard, NU_CACHE_RDY); Nov-24 moment-audit defects
      F1-F4 already fixed in 910c4667. Five defects found and
      FIXED this commit: (1) DO_KNOP_REFRESH absent-key default
      evaluated before the absent-key USE_LAM_ES flip, resolving
      armed-T on the V13 route - misleading knop_refresh=T echo
      plus a per-iteration MOM_KNOP_REFRESH_V1 invoke (inert via
      KN_M_READY); absent key under plain SOB+LAM_ES=T now
      resolves F, matching what run_model_grid always stamps,
      and an explicit T warns it is inert. (2) Three shared-code
      MOM-store gates used the wide KNOP.OR.SOB predicate while
      the store exists only on SOB+LAM_ES=F - SET_MOM_GS_TAG
      (cmfgen_sub), RESTORE_MOM_JNEXT_LATEST x2 (solveba_v13),
      MOM_KNOP_REFRESH_V1 (cmfgen_sub) all narrowed.
      (3) Documented SOB solver-mode rejections (USE_FIXED_J,
      USE_FORMAL_REL, USE_J_REL, USE_DJDT_RTE, PP_MOD, PP_NOV)
      were never coded - STOP added mirroring the KNOP block.
      (4) MOM_BLKGMRES/MOM_NUSOLV/MOM_H_IMPL set T off the moment
      route silently dead - CHECK_PARAM_CONSISTENCY inert
      warnings added. (5) USE_LAM_ES=F banner lacked
      "unsupported for production" wording - added.
      Tests: test_vadat_keys 99->103 checks all PASS (new:
      absent-key refresh resolves F, explicit T honoured+warns,
      J_REL/PP_MOD rejections, MOM_BLKGMRES dead-key warning,
      unsupported banner); test_iter_conv_gate 18 PASS (new
      gate pins); test_obs_nonmono_dispatch PASS; test_dvc
      170 PASS; test_iter_output PASS. Pre-existing rot fixed:
      test_vadat_keys HELPERS missed iter_log_fmt_mod.f90 since
      0cea0a27 added ILF calls to mom_j_cmf_nonmono_v1 - link
      was broken before this work. Stale comment fixed in
      scripts/run_model_grid.py (still claimed USE_LAM_ES=F
      default). Remaining known gaps (deferred): DUMP_MOM_JNEXT_GS
      called unconditionally in scr_read_v2 (early-return no-op,
      upstream file left untouched); GS_STATE stale-file read at
      solve_for_pops.f:697 ungated but all consumers sit inside
      HAS_LAGGED_FIELD; JCOMP diagnostic vacuous under
      LAM_ES=T (RJ=TC by construction).

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

## Coupled wavelength-depth block solve scoping (24-Nov-2026, user-directed)

Question: viability + cost of replacing the lagged frequency
iteration with a Baron & Hauschildt coupled block solve.
Plan: docs/plans/BLOCK_SOLVE_SCOPING.md, phases A-E.
Exit criterion: per-phase verdicts; implementation only after
phase-D design sign-off. Toy prototype (phase C) is the decisive
cheap evidence.

## Coupled wavelength-depth block solve implementation (24-Nov-2026, user-directed)

Question: can the matrix-free, preconditioned wavelength-depth GMRES
solve replace the lagged frequency pass in the MOM refresh without
disturbing the default path.
Evidence plan: E1 dense-reference manufactured-graph battery; E2
standalone production-module tests to machine precision; E3
VADAT-gated integration behind MOM_BLKGMRES; E4 full unit batteries
and one bounded eps1em3 arm. Track GMRES iterations, final residual,
rejected writes and refresh residuals.
Exit criterion: converges, stalls bounded, or detonates; verdict
recorded here with job id and commit. Any negative or nonfinite write
amplification beyond existing guards stops the arm for diagnosis.

Status 24-Nov-2026: E1-E3 implemented and pushed through `66cb87c7`.
`test_mom_blk` 59/59, including dense-reference machine-precision
solves, a signed upwind edge, and the solver-error contraction probe
(lagged growth 1.847/iter on the marginal ring; block error ratio
3.27e-15). `test_mom_solve` 144/144, `test_mom_nu` 53/53, and
`test_vadat_keys` 91/91. The capture sign audit is closed analytically
and by the signed manufactured case: capture satisfies
`NU_X0 + NU_CUP*J(NU_KUP)=XM_SAVE`, so the implicit row is
`tri*x - NU_CUP*xup = NU_X0`; PSI remains only in `NU_TB`. The solver
now reports the true final preconditioned residual rather than only
the Arnoldi estimate. The production archive and `cmfgen_dev.exe`
were refreshed; frozen sha256
`3763d79869ca12c48c5bd7144bb81c876fa0c0b2ecf1099975f70a2173951cd4`.
E4 arm submitted through `run_model_grid.py`: CMFGEN `13655952`, OBS
`13655953` (dependent), same eps1em3 seed/deck with `USE_LAM_ES=F`,
`MOM_NUSOLV=F`, `MOM_BLKGMRES=T`. CMFGEN is running; verdict
pending. Early runtime evidence at it670: four refresh cycles
completed in 2.6-9.1 s; `blk-gmres` never exceeded five Krylov
iterations, used no restarts, and held final residuals below 1E-3,
though 33 of the first 55 lines needed zero iterations because the
incoming iterate already satisfied the captured-system tolerance.
The outer field response remains the open question. The full
`run_all_tests.py` pass was stopped after the solver source changed
mid-run; before that point all moment and VADAT batteries were green,
while the known unrelated pre-existing reds remained `test_auto_ion`
(missing fixture), `test_cfr4` (harness decode), `test_es_conv`
(obs archive gap), and the documented `test_knop_lagged_map` sdab
failure. `test_fg_sw_sym` and `test_fnv1a` pass after the companion-
object fix. The long fixed-source `nu-census` WRITE was split onto
continuation lines after the running arm showed only the a>0.9 prefix;
`test_mom_solve` now pins the presence of the `rmax` and `H/X0` fields.
The alteration manifest was also updated to register
`plane/mom_blk_gmres_mod.f`. A final full pass from the committed state
is still owed.

## Knop remediation WP-2: freshness contract and order modes (26-Nov-2026, task-directed)

WP-2 makes the ordered march's freshness claim measurable and
adds the per-ray order-mode selector. `KN_RAY_ORDER_MODE` (0
ordered / 1 alternating / 2 alternating seeded down the dominant
coupling sign; 3 reserved SCC, stops) resolves a -1 sentinel from
`KN_RAY_ORDERED`, so driver pokes keep their semantics while an
explicit key wins with a warning on conflict. `KN_FRESH_DIAG`
arms a per-pass census of the |Q| donor weight read from a
not-yet-updated cell, stored per (pass,dir,ray) in `KN_FRESH_FR`;
the false "never waits a pass" header comment is corrected -- the
off-row JPR share is a full pass stale.

Fixture: `ord` reads ~0.943 (=WINX) of its compressing donor
weight stale on every I_M pass and ~0 on I_P; `uno` alternates
1.0/0.0; `m2u` opens fresh and converges in 2 passes (mode 1
wastes its first pass on a fully stale read); `m2m` is bitwise
identical to `unom` -- mode 2 reduces to mode 1 on mixed rays.

Deliverables: `docs/audits/KNOP_REMEDIATION_WP2_2026-11-26.md`;
`test_knop_ordering` 50/0 (+11 WP-2 checks); `test_vadat_keys`
120/0 (+6 WP-2 checks); `test_knop_ray` 50/1 (T18 pre-existing),
`test_knop_spv` 93/0, `test_iter_conv_gate` 18/0,
`test_conv_exit` 53/0, `test_fg_gs` 36/0, `test_knop_oneway` 8/0.

## Knop remediation WP-1: defaults, tolerance-ladder invariants, domain census (26-Nov-2026, task-directed)

WP-1 of the Knop convergence remediation: config-level changes
only. Defaults flipped - `KN_RAY_ORDERED` T->F (alternating
sweep; the ordered march's inserted-point donor is a full pass
stale per the WP-0 fixture), `KN_NI1_FULL` F->T (vacuous FULL
ends the ninth-conjunct deadlock on single-point rays),
`KN_RAY_BULK_TOL` 6e-2->0 on every route (the bulk floor was
bounded non-contraction, not convergence). `GS_MATURE_TOL`
hoisted to a `CONTROL_VARIABLE_MOD` PARAMETER and exposed as a
read-only VADAT key; `KN_RAY_TOL_V`/`KN_RAY_BULK_TOL_V`/
`FG_COMMIT_NC_V`/`GS_MATURE_TOL_V` promoted for echo.

Hard invariants in `CHECK_PARAM_CONSISTENCY` (STOP) on the
FG-census routes (`DO_NONMONO_KNOP`, `DO_SOB_FG_KNOP`): I1
`KN_RAY_BULK_TOL < GS_MATURE_TOL`; I2 `KN_RAY_TOL <=
0.1*GS_MATURE_TOL`; I3 explicit `FG_COMMIT_NC /= GS_MATURE_TOL`
(absent-key legacy default exempt until WP-4; presence staged
at read time into `FG_COMMIT_NC_SET` because `CLEAN_RD_STORE`
runs before the checker); I4 domain census in the FG refresh -
`KN_NACT`/`KN_NFULL` count only `NI_RAY >= 2` rays and
`FG_RF_ALL_CONV` requires `KN_NINV=0` over the full production
mask, so a non-domain INV ray still fails closed.

Deliverables: `docs/audits/KNOP_REMEDIATION_WP1_2026-11-26.md`;
`test_vadat_keys.py` 114/0 (+9 WP-1 checks); `test_knop_ray.py`
50/1 with re-pinned census cases (T29a blocks via NINV, T29b
3-ray mixed ensemble reads TRUE off the domain rays, new T29c
`ni1s` pins the zero-domain fail-closed edge);
`test_iter_conv_gate` 18/0, `test_conv_exit` 53/0,
`test_knop_ordering` 39/0, `test_fg_gs` 36/0. Sole remaining
failures are byte-identical pre-existing ones (T18 `diag=False`,
`sdab` KNOWN DEFECT, `test_knop_lagged_field` seed_arm_probe).

Verdict: PASS. STOP per spec - Arm A warm restart from the
it122 checkpoint on `dessart_m15_2p82d_FRESH_knop_v14` with the
new defaults is the human's to submit.

## Knop remediation WP-0: ordering-stall reproduction fixture (26-Nov-2026, task-directed)

WP-0 of the Knop convergence remediation: a manufactured ray that
reproduces the SWEEP_FLD_ORD stall on the current default and
converges under KN_RAY_ORDERED=F. Deliverables:
`unit_tests/f90/drv_knop_ordering.f90`,
`unit_tests/python/test_knop_ordering.py` (39 checks, all PASS),
`docs/audits/KNOP_REMEDIATION_WP0_2026-11-26.md`. No production code
changed.

Fixture: NI=ND=140, NCF=200, every interior row inserted at
WINX=0.95 toward the deeper node so the ordered march's compressing
donor reads 0.05*own-row (fresh) + 0.95*deeper-row JPR2=J+1 (a full
pass stale). |Q|/(chi+|Q|)=0.999, matched NFM=3, store seeded 1e-6.

Measured on the current default: ord exits IXIT=6 (plateau) at
NPASS=25 with RESID=4.12e-2 and RHO_M=0.949 - inside the stalled
production-ray signature (IXIT=3/6, res 2.9-8.8e-2, np 10-37). uno
(same physics, FG_ORD_SWEEP=F) exits clean in 3 passes at
RESID=4e-16. ordx (exits off) hits the 60-pass cap at RESID=1.588e-2,
the same floor an independent Python replica produced. ordb
(drift-dominated unmatched algebra, chi_hat*dz=30, |Q|=0.999) exits
IXIT=6 at 17 with RESID=0.86.

Recorded deviations from the task spec (documented in the audit,
asserted in the test): NFM=3 + DTAU=30 is unreachable under matched
spans (|Q|*chi_hat*dz = NFM_raw/NFM <= 1), so DTAU=30 is carried by
the unmatched ordb/unob arm; the spec's "<=6 passes alternating"
mixed-ray claim fails at |Q|~1 on both signs (unom non-contracts at
~0.95/pass on the matched grid - the nbulk lagged-map finding, not
an ordering artifact); on the unmatched grid expanding rows anchor
the chain and both marches converge.

Verdict: PASS - the stall is reproduced on the default and cleared
by the unordered sweep; the fixture gates WP-1..WP-4.

## Block-solve arm verdict and outer-loop investigation (25-Nov-2026, user-directed)

E4 verdict after user-requested cancellation: bounded outer-loop stall,
not a detonation and not convergence. CMFGEN 13655952 was cancelled at
iteration 1039 after 11.0 h; OBS 13655953 was cancelled with it. The
linear solve itself was strong: median two and maximum five Krylov
iterations per pass, zero restarts, median final residual 4.18E-4, and
714 of 746 refresh writes marked converged. Recent moment quality was
orders better than the chain arm: 4.1% of cells above 5% and 0.38%
above 50%, versus 69.1% and 36.6% under MOM_NUSOLV. The production
model still did not converge: it remained in LAMBDA_TFIX, population
demand oscillated at 1E5-1E9 percent scale with cap hits, recent
accepted largest moment-store writes reached 1E10-1E14, and the
outer/inner luminosity ratio stayed flat near 0.33. The chain probe
13654121 was also cancelled at iteration 1080; all 828 refresh writes
failed and retained the pre-refresh field, a harder bounded stall.

New bounded investigation, user directive 25-Nov-2026: determine why a
linearly converged MOM_BLKGMRES refresh still leaves the outer
population/field loop non-contractive, then test the smallest justified
source or execution-policy change toward a converging eps1em3 model.
Evidence plan: audit write-back, best-pass acceptance, active-set,
adaptive relaxation, outer AA, and the lagged H/EPS channels against
the cancelled-arm output; use the existing opt-in MOM_H_IMPL arm as the
first controlled execution experiment because H is the largest known
lagged channel still outside the block solve. Exit criterion: a revised
arm converges, stalls with a sharper bounded diagnosis, or detonates;
the verdict and mechanism are recorded here before any broader campaign.

Execution experiment submitted 25-Nov-2026:
`m15_2p82d_VAL_sob_lamesF_blkgmres_himpl`, CMFGEN 13671005 and OBS
13671006, same eps1em3 seed with `MOM_BLKGMRES=T` plus `MOM_H_IMPL=T`.
Pinned `cmfgen_dev.exe` sha256
`60cdf8306f7d24f18abcb390cea66793fe78ef9306b5e148223ad6e4414c3819`;
the job is queued behind `QOSGrpCpuLimit`. The first submission attempt
was cancelled before start because the copied YAML kept the old
`model_dir`; that mistake refreshed the cancelled blkgmres run directory
and discarded its ITER_HISTORY/ITER_LEDGER files. The original OUTGEN
through iteration 1040 was preserved at
`cmfgen_validation/m15_2p82d_VAL_sob_lamesF_blkgmres/OUTGEN_cancelled_it1040`.
The corrected config uses its own `_himpl` model directory and a staged
copy of the parent SN_HYDRO_DATA.

Initial audit lead: the converged-refresh path keeps the iterated field
without applying the 50x retained-store jump guard used for best-pass
fields. That guard was written for the non-converged branch only. This
may explain how bulk-converged block refreshes can commit 1E10-1E14
largest-cell excursions while passing the bulk census. Do not change it
until the H_IMPL controlled arm separates outer H regularization from
write-back policy; otherwise the two changes cannot be attributed.

Cold-start arm authorised by the user 25-Nov-2026:
`m15_2p82d_VAL_sob_lamesF_blkgmres_cold`, same eps1em3 deck with
`MOM_BLKGMRES=T` and `MOM_H_IMPL=F`, but no warm-restart checkpoint.
CMFGEN 13671134 and OBS 13671135 were submitted 25-Nov-2026 with the
same pinned executable sha256 as the H_IMPL arm. The resolved config
has no `warm_restart_from` entry. Question: can the coupled moment
route converge as a self-consistent population/field iteration once the
USE_LAM_ES=T solver-swap shock is removed. Exit criterion: convergence,
a bounded stall reproducing the warm-start behaviour, or detonation.
This arm takes precedence over any write-back source edit because its
verdict decides whether the warm stall was mostly an initial-condition
shock or a structural outer-loop defect.

Cold evidence through iteration 13 shifted the next controlled arm to
write-back policy. The block solve remains cheap (0-5 GMRES iterations,
zero restarts), but the adaptive V2 omega estimator repeatedly read a
clamped 5.0E-2 contraction factor from alternating pass residuals and
raised `JNEXT_RELAX` toward 0.84. Iterations 12-13 then committed
`7.5E8` and `6.7E9` largest retained-field changes while the
outer/inner luminosity ratio stayed near 0.30. The warm `MOM_H_IMPL`
arm also remained lum-deprived (ratio about 0.04 at iteration 680), so
the next controlled experiment is fixed damping rather than another
implicit-H change.

Verdict (26-Nov-2026, both arms cancelled once the answer was in -
CMFGEN 13671005, 13671134 and OBS chains 13671006, 13671135):
- himpl (warm restart, MOM_H_IMPL=T): killed at it1409 after ~1400
  LAMBDA iterations with no BA handoff; population max-increases were
  still 5.2E7-6.6E7% at the last iterates (earlier peaks to 1E16).
  The implicit-H channel did not restore outer-loop contraction -
  detonation, not a bounded stall. Consistent with the augres
  finding that H is the dominant residual damage: regularising one
  lagged channel inside a J/H pair that does not solve the augmented
  system cannot fix the fixed point.
- cold (no warm restart, MOM_H_IMPL=F, adaptive relax): killed at
  it722, still LAMBDA warmup, max-increases 4.7E4-8.5E4%. Removing
  the solver-swap shock did not change the outcome - the
  non-contraction is a structural outer-loop defect on the
  USE_LAM_ES=F/MOM_BLKGMRES route, not an initial-condition shock.
  Confirms the route-table FAIL; next controlled step remains the
  augmented J+H implicit solve (FIX 4), not further regularisation.

Execution experiment submitted 25-Nov-2026:
`m15_2p82d_VAL_sob_lamesF_blkgmres_r01`, CMFGEN 13671873 and OBS
13671874. Same cold-start eps1em3 deck and pinned executable as the
cold baseline, but `KN_ADAPT_RELAX=F` and `JNEXT_RELAX=0.1D0`. Question:
does conservative fixed moment write-back remove the outer population/
field limit cycle that adaptive relaxation sustains. Exit criterion:
convergence, a lower-amplitude bounded stall, or failure to improve on
the adaptive cold arm. The submission warned of 21 pre-existing dirty
paths and reused pinned `cmfgen_dev.exe` sha256
`60cdf8306f7d24f18abcb390cea66793fe78ef9306b5e148223ad6e4414c3819`;
no solver source changed for this arm.

Literature-correctness audit declared (25-Nov-2026, user-directed):
question - does the current `USE_LAM_ES=F` implementation still solve
the radiation moment equations intended by the cited non-monotonic
literature, and are the block-GMRES linearization, boundary conditions,
Eddington closure, refresh acceptance, relaxation, and Anderson
acceleration mathematically consistent with that target. References to
check: Knop, Hauschildt and Baron 2009 Papers I/II, Baron and
Hauschildt 2004, Hillier and Miller 1998, Dessart, Audit and Hillier
2015 where applicable. Evidence plan: line-level audit of
`plane/mom_j_cmf_nonmono_v1.f`, `plane/mom_blk_gmres_mod.f`, dispatch
and capture wiring, original-paper equations, and live-arm diagnostics.
Exit criterion: per-component VERIFIED/PARTIAL/DEFECT verdicts with
prioritized fixes recorded in `docs/audits/` and here. The audit is
read-only for solver source while the three pinned-executable arms are
running.

Audit verdict recorded 25-Nov-2026 in
`docs/audits/MOMENT_LITERATURE_AUDIT_2026-11-25.md`: PARTIAL with
defects. The sign-keyed frequency upwinding, positive generalized
opacity, depth-tridiagonal solve, and GMRES kernel are usable. The
current block system is not the full literature wavelength boundary-
value problem because H and EPS coupling remain explicit in `NU_X0`,
and `H_NEXT_FRQ` is not reconstructed after the block J update. The
printed GMRES residual is measured before negative-cell rejection,
write caps, damping, and AA, so it does not certify the committed
field. LAG_BND, EPS limiting, XM_POS, the EPS depth-upwind stencil,
and MOM_H_IMPL change the equations or the fixed point unless proven
inactive. The V2 adaptive relaxation estimator mixes alternating
directional maps; fixed relaxation remains the valid control. Final
write-back acceptance is path-dependent because converged refreshes
bypass the jump guard and outer AA can write after that guard.
Decision: implementation changes DEFERRED pending user direction; the
first correctness task is a full discrete H/EPS wavelength-coupling
derivation and a post-write residual measurement. The running pinned
arms were left unchanged.

H/EPS derivation started 25-Nov-2026. Result recorded in
`docs/plans/MOMENT_H_EPS_IMPLICIT_DERIVATION_2026-11-25.md`: direct
substitution of the H recurrence is the wrong implementation because
alternating upwind directions make the eliminated frequency coupling
dense or ring-shaped. The clean target is an augmented sparse system
solving nodal `J` and midpoint `H` together. The smallest safe first
code change is a diagnostic residual evaluator for that augmented
system, before any solver rewrite. Spatial boundary rows, frequency
boundary data, EPS donor-stencil choice, row scaling, and matched-
frequency H strides remain the five derivation points that must be
closed before Phase 1.

## Augmented J+H residual evaluator phase 0 (25-Nov-2026, user-directed)

Question: does the committed J_NEXT_FRQ/H_NEXT_FRQ field solve the
full coupled moment equations, not just the bounded partial system.
Plan: implement the augmented-residual diagnostic derived in
`docs/plans/MOMENT_H_EPS_IMPLICIT_DERIVATION_2026-11-25.md`; report
RESID_IMPLEMENTED (rows as assembled, stabilisers active) and
RESID_TARGET (same rows with LAG_SCL=1, uncapped EPS, no XM_POS clip,
no MOM_H_IMPL mass, centred V11 EPS midpoint stencil) separately for
moment and H rows; extend the blk-gmres report with raw and
post-write-policy residuals. Diagnostic only - no solver rewrite.
Exit criterion: evaluator verified against an independent
reimplementation; both residuals reported on the committed field.

Status: implemented behind `MOM_AUG_RES` (default F). New routine
`MOM_AUG_RESID_V1` evaluates the augmented z(k) = [x, h interleaved]
system including both spatial boundary rows and the explicit
frequency-upwind H graph; it runs at each MOM_KNOP_REFRESH_V1 exit
after all write policies and is read-only. New capture arrays carry
the H-row and lagged-term coefficients (`NU_HU/NU_HL/NU_HS`,
`NU_VB/NU_VC`, `NU_EPS(U)`, `NU_EPSP(U)`, `NU_LAGSCL`, `NU_SRC`,
`NU_CUPU`, `NU_T1B`, `NU_KUH`, `NU_OUTBC`, `NU_HCHK`). The blk-gmres
line now prints raw `||b-Ax||/||b||` next to the preconditioned
residual plus the committed-field residuals measured after
rejection/caps. `test_mom_augres` 58/58 PASS (numpy cross-check,
machine-precision manufactured field, J/H perturbation detection,
both boundary types, EPS cap/LAG_SCL/XM_POS/H_IMPL splits, rings and
chains, off-grid reads, skipped columns, NaN counts, state untouched).
Moment regressions green: test_mom_solve 144, test_mom_nu 53,
test_mom_blk 67, test_mom_gs 34, test_mom_aa_state 62,
test_mom_sweep_table 8, test_vadat_keys 102. Default path unchanged;
flag off reaches no new code. Build caveat this session: a concurrent
unrelated kv-block conversion landed mid-build twice (jfnk interface,
mon_interp link), so full `./build.sh` completion is still owed; every
archive member the batteries consume is fresh and all 528 checks pass.

## Moment solver FIX 1 + FIX 2 (25-Nov-2026, next phase of the augres plan)

Question: can the two correctness defects the Phase-0 audit
confirmed be removed without touching the operator - the
MOM_ADAPT_OMEGA_V2 estimator ratioing residuals across
alternating directional maps (audit F), and the outer-AA
commit bypassing every pass-level guard (audit G).
Evidence plan: implement proposals FIX 1 (same-direction
residual ratios) and FIX 2 (single post-AA commit gate) from
docs/plans/MOMENT_AUG_RES_PHASE2_FIXES_2026-11-25.md; verify
on a manufactured-input driver rather than a solver sweep.
Exit criterion: estimator within 5% of known rho on
manufactured maps and no AA-to-commit path past the gate.

Verdict FIX landed, both proposals. FIX 1: RES_PASS now holds
one slot per sweep direction and MOM_SAME_DIR_RHO ratios only
consecutive same-direction entries (stride-2 under ALT), so
the estimator reads each map's own contraction instead of the
cross-map beat. Falls back to the previous cross-call ratio
only while fewer than three same-direction ratios exist.
FIX 2: JBEST/HBEST snapshot the field before the outer-AA
block and MOM_RF_COMMIT_GATE runs once after AA and before
the maturity census - non-finite or 1e40 values, J<=0 cells,
50x jumps, and bulk regression of the pre-AA n5/n50 counts
each reject and restore the snapshot atomically with a WARN
per cause. New unit pair drv_mom_gate + test_mom_gate
(45/45): rho recovered on manufactured alternating maps,
diverging map reads 1.69, the old beat false-positive is not
flagged, and the gate rejects plus restores on every
violation class. Regressions green: test_mom_solve 144,
test_mom_augres 58, test_mom_aa_state 62, test_mom_nu 53,
test_mom_blk 67, test_mom_gs 34, test_mom_sweep_table 8,
test_vadat_keys 102. Default path unchanged - V2 estimator
remains opt-in and the gate runs inside the existing
LAA_OK block. Still owed: same-direction ratio consistency
on a real capture; FIX 3-5 wait on the first real-case
MOM_AUG_RES arm for a damage ranking.

## First real-case MOM_AUG_RES arm (25-Nov-2026)

Question: on a real non-monotonic m15 iterate, how large is the
implemented-vs-target residual gap on the committed field, and
which stabiliser (LAG_BND, EPS cap, XM_POS, H_IMPL mass, EPS
donor stencil) accounts for most of it - the ordering evidence
for FIX 3-5.
Arm: configs/dessart_m15_2p82d_VAL_sob_lamesF_blkgmres_augres.yaml,
the r01 deck (USE_LAM_ES=F, MOM_BLKGMRES=T, JNEXT_RELAX=0.1 cold
start, ND=70, FIX_T pinned) plus MOM_AUG_RES=T; debug submission
on the interactive partition (job 13685499; the first attempt,
13685125, surfaced a latent internal-write overflow in
cmfgen_sub.f's F-ITER summary - TMP_STRING*80 vs an 88-char
message - fixed and rebuilt). Live arms checked first - the
three lamesF_blkgmres runs carry no MOM_AUG_RES and answer
different questions.
Evidence plan: mom-aug lines at each MOM_KNOP_REFRESH_V1 exit
give the two residual norms and the per-mechanism counters;
extract the series over the captured iterations.
Exit criterion: a residual-damage ranking with magnitudes, or a
documented statement that the diagnostic produced no usable
exits; the arm is cancelled as soon as that answer is in.

Verdict (25-Nov-2026, arm cancelled once the answer was in):
8 evaluator exits captured over the first ~50 min of a cold
start; every refresh stalled and retained the pre-refresh field
(stored-field change 0.0), so the series measures the retained
store under drifting row coefficients. Results:
- Scaled moment residual s2mom: 81 at it1 rising to ~330 and
  plateauing. IMP vs TGT gap: 18% at it1 shrinking to ~1.5-4%
  later - on a diverging field the equation modifications
  explain only a small share of the committed-field residual;
  non-convergence itself dominates.
- Scaled H residual s2H: 152 at it1, then 1e14..1e56 - the H
  rows are catastrophically inconsistent and their share grows
  without bound. This is audit finding B made quantitative:
  H is reconstructed per sweep column and never re-derived
  after the block solve, so the stored J/H pair does not solve
  the augmented system by construction. Largest measured
  residual damage by orders of magnitude; evidence for FIX 4.
- Per-mechanism site counts, stable across all 8 exits:
  lagb ~1.10-1.16M rows scaled every exit (pervasive, not
  rare), hclip ~157-172k, xmclip ~29-35k, epscap 0 every exit
  (the EPS cap never fires on this model), offgrid ~1.2-1.5k,
  bndrows ~400-590, nan 0.
- The FIX 2 commit gate fired 6 times on real data, all
  bulk-regression rejections: the outer-AA write moved
  7.8k-13.8k cells >5% on fields whose guarded path committed
  zero movement (pre-AA counts 0/0 every exit), and each write
  was vetoed and restored. The pre-AA baseline is 0/0 whenever
  the whole sweep is rejected, so AA churn on a rejected sweep
  always reads as regression - conservative and intended.
Damage ranking for FIX 3-5: (1) FIX 4 implicit H coupling -
  s2H is the dominant residual by orders of magnitude; (2)
  FIX 5 LAG_BND retirement - pervasive (~1.1M rows) but its
  measured share of s2mom is small; (3) XM_POS ~30k cells;
  (4) FIX 3 EPS work last - the cap is dormant on this model
  and the stencil share sits inside the small IMP-TGT gap.
FIX 1 note: the V2 estimator was not exercised (r01 deck pins
JNEXT_RELAX); its real-capture consistency check stays open
and needs an MOM_ADAPT_RELAX=T arm.

## FIX 4 stage-1: augmented-solve contraction prototype (25-Nov-2026)

Question: does the augmented operator (lagged H bracket and H
reconstruction as implicit rows, z = [x, h interleaved], no
diagonal mass) contract on a case where the lagged map requires
MOM_H_IMPL - and what does the mass cost the fixed point.
Evidence plan: unit_tests/python/test_mom_himpl.py builds both
operators in numpy from the drv_mom_augres captured-array dumps
(the "hcon" case = epsmix capture with NU_VB/NU_VC x30) and
measures the lagged map spectral radius with and without mass,
the massed-vs-true fixed-point shift, and GMRES iteration counts
on the augmented operator.
Exit criterion: the FIX 4 proposal's own - the augmented
prototype contracts without mass where the lagged map needs it.

Verdict: exit criterion met; FIX landed at prototype stage.
- rho(lagged, no mass) = 2.55 -> diverges; rho(massed) = 0.83.
- The massed fixed point is shifted 57% in x from the true
  solution - the uncompensated |VB|+|VC| diagonal buys
  contraction at a large model error, quantifying the cost the
  MOM_H_IMPL comment already flagged qualitatively.
- The unmassed augmented operator solves exactly (residual
  2.8e-16, manufactured field recovered to 1.8e-15); GMRES
  converges in 55 iterations bare, 31 with the per-frequency
  dense-block preconditioner.
Next: stage 2 is a production-side matrix-free matvec for the
augmented operator behind a flag, credited by RESID_TARGET on a
real capture; MOM_H_IMPL retires only after that lands.

## Moment-machinery audit (24-Nov-2026, user-directed)

Question: is the moment-equation procedure in mom_j_cmf_nonmono_v1.f
correct and robust - row assembly (TA/TB/TC/XM), H reconstruction,
EPS sign-splitting, LAG_BND scaling, XM_POS clipping, boundary rows,
upwind conventions, capture consistency.
Evidence plan: line-level static audit plus python-driven tests
closing coverage gaps beyond test_mom_nu/solve/gs/aa_state batteries.
Exit criterion: verdict per component (VERIFIED/DEFECT), fixes or
rejected hypotheses recorded in docs/audits/.
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
      A/B VERDICT (24-Nov-2026, jobs COMPLETED): ngshadow
      converged it972 in 5h49m; ngoff converged it984 in
      6h03m - NG cadence net-positive by 12 its / 14 min
      (~4%) from the same perturbed warm state. 19 events
      at cadence 16 (it672..960). Telemetry across all
      events: production N=4-7 fits NEVER settled (settle
      1.4..160; |c|1 ~ 1e3..4e4 late); committed kicks
      decayed 3.2%/2.8% -> 0.08%/0.09% tracking demand.
      Kick value asymmetric: mid-run kicks net-positive
      (rebound re-damped faster than unaccelerated tail);
      late kicks ~neutral (it960 rebounded 1.1e-3 -> 6.4e-2
      then re-damped to eps in ~10 its, about what the tail
      alone needed). Genuinely settled deep-window
      predictions exist at N~11-12 but propose ~100x
      smaller moves (1e-4 vs 2-3e-2), and it688 showed the
      degenerate failure mode (|c|1=1.00 = prediction
      copies newest iterate) - a settle gate needs a
      coefficient-norm guard to be safe. Mechanism answer:
      the cadence fires unsettled extrapolation gambles
      that pay mid-run; a Magritte settle-gate would defer
      ~all of them (equiv to ngoff, ~12 its slower on this
      model) so adaptive-acceptance alone loses the small
      win. Only variant the data justifies: a near-eps
      abstention (skip NG when demand < C*eps with C~5-10,
      or when predicted move >~ remaining demand) worth
      ~1% - below the noise floor for a third arm unless
      a new deck shows bigger kick waste. NO third arm
      submitted; keeping NG cadence as is. NG1 mechanism
      question answered; gate-design constraints recorded
      for any future NG_ADAPT proposal.
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
      ARM (24-Nov-2026): configs/dessart_m15_2p82d_VAL_sob_
      lamesF_nusolv.yaml - identical deck and identical eps1em3
      warm-restart seed to the stalled arm, MOM_NUSOLV=T the only
      change. Question: does the per-depth chain solve contract
      the compressing-direction map where Jacobi passes could
      not? Evidence: nu-solve census per pass (rings, singular,
      frozen, rejected; NU_CHG trend), refresh residual
      trajectory vs the stalled arm's plateau ~1e4, MOM_RF_CONV /
      MOM_J_GS_OK firing, lam-field-not-ready release. Exit:
      converges (route viable again), stalls identically (chain
      solve insufficient - the stiffness lives in the lagged
      H/EPS terms still explicit or in XM_POS/LAG_BND), or
      detonates (capture/consistency defect on the real grid).
      Single arm, no chain.
      VERDICT (24-Nov-2026, jobs 13651766 + 13652626, killed at
      it675 after ~9 it): STALLED with amplification - REFUTED as
      the sole fix. The chain solve runs in production exactly as
      unit-verified: 208 calls, rings resolved 5.8K -> 2.3K per
      pass, 0 singular, 0 frozen on bad coefficients, negative
      rejections declining 20K -> ~460/pass. But every refresh
      still stall-exits; residuals oscillate 37 - 2.7e4 with no
      contraction trend (the 37-250 minima are lower than the
      dead arm's ~2e4 plateau - the implicit-J map intermittently
      touches near-convergence but cannot stay). Store commits
      ESCALATE vs the dead arm's bounded 1e6-1e8 limit cycle:
      after two vetoed iterations the committed moves ran
      6.9e7 -> 4.6e11 -> 1.0e23 -> 1.8e8 with 697K-979K cells
      moved >5%. Populations dragged off the converged premise:
      outer luminosity 4.27e8 -> 8.9e7 (~5x), raw MAXCH demand
      to 1e26. JNEXT_RELAX wanders 0.48-0.80, estimated
      contraction factor 0.05-1.15 - no stable contraction.
      Mechanism: the ring solve x_j*(1-A)=S_j amplifies
      near-marginal cycles (cyclic product A -> 1 IS the rho~1
      sign-flip mode); the 1e-10 freeze threshold only catches
      exact degeneracy, so A~0.9999 rings return amplified
      exact solutions the damped write-back passes through, and
      once the >50x veto opens the commits explode. Conclusion:
      the composite map's stiffness is not the lagged-J graph
      alone - lagged H/EPS terms stay explicit, and marginal
      modes need treatment beyond an exact solve (physical
      |1-A| freeze band, per-cell nu-write cap, or implicit
      H-channel coupling). Side results: two format crashes
      fixed en route (mine at :6029; latent pre-existing
      best-pass-rejection WRITEs :6707/6713/6719, commit
      4a79500b). OBS chains 13651767 + 13652627 cancelled on
      dead dependencies.
      PHASE 2 (24-Nov-2026, user check-in granted): plan doc
      docs/plans/MOM_NUSOLV_MARGINAL_MODES_TODO.md. Question:
      do a near-marginal ring freeze band (|1-A|<1e-2 keeps old
      iterate) plus a per-cell nu-write trust region (cap ~1e2
      on the normalized ratio) stop the amplified commits while
      keeping the exact solve on non-marginal modes? Evidence:
      new NMARG/NCAP census vs NRING, commit magnitudes vs the
      1e23 excursions, residual trend, population drift. Exit:
      bounded stall without excursions (marginal-ring hypothesis
      confirmed, decide implicit-H phase next), still amplified
      (hypothesis wrong), or converges.
      LANDED (24-Nov-2026, commit 4d78adaf): freeze band +
      write cap in MOM_NU_CHAIN_SOLVE; the FSCL trust-region
      scale is the incoming column's own scale. Unit verdict:
      test_mom_nu 44/44 (marg frozen, margok solved exact
      above the band, cap clamps to the cap value to machine
      precision); test_mom_solve 143/143 fixed-point regression
      untouched (guards inert on comp1's contractive map);
      mom_gs 34/34, vadat_keys 89/89, aa_state 62/62,
      sweep_table 8/8.
      ARM (24-Nov-2026): configs/dessart_m15_2p82d_VAL_sob_
      lamesF_nusolv2.yaml - same eps1em3 seed, MOM_NUSOLV=T.
      CMFGEN 13653300, OBS 13653301. Pinned exe 3fc083cd
      carries 4d78adaf.
      VERDICT (24-Nov-2026, job 13653300, killed at it675
      after 9 it): FAIL - guards correct but the mechanism
      hypothesis is refuted. NMARG=0 on every pass of every
      iteration: no ring on the real grid is near-marginal,
      so the freeze band never engages - the amplification
      never came from near-singular 1/(1-A) division. The
      write cap binds (59-7210 cells/pass) but cannot contain
      the excursions: it limits one pass's write to 100x, and
      20-40 passes per refresh compound past it - commits hit
      7.5e23 largest change with 685K-1M cells moved, the
      same class as the unguarded arm. Populations degrade
      identically (outer L 9.1e7 by it675). Sharper mechanism:
      in the explosive regime the exact solve of stale-input
      equations amplifies through chain composition and the
      lagged source magnitude - loop gain above unity no
      matter what guards the write. The lagged splitting is
      the defect, not its ordering: Baron & Hauschildt 2004
      solve the coupled wavelength-depth system directly, no
      lagged iteration. Unit-side the guards are verified
      exact (test_mom_nu 53/53 incl. the iterated probe:
      unguarded growth measured exactly 1000.00/iter = the
      c/(1-a) amplification rate) - the failure is physical-
      regime, not implementation. OBS 13653301 cancelled on
      the dead dependency. Route stays FAIL/opt-in; the only
      contraction-capable lever left is reducing the lagged
      fraction itself (implicit H channel, or a coupled
      wavelength-depth block solve) - needs a user decision
      before a phase 3.
      PHASE 3 (24-Nov-2026, user directive "investigate and
      optimise the lagged field"): Question: which lagged
      terms dominate the loop gain in the explosive regime -
      the H-channel upwind reads, the EPS_PREV midpoint
      transfers, the other-depth J terms, or boundary reads -
      and can reducing their lag (implicit-H capture on the
      same graph, or tighter update timing) drop the composite
      loop gain below unity? Evidence plan: (1) enumerate
      every lagged term in the row assembly and its update
      timing; (2) instrument the capture to census lagged-term
      magnitudes vs the diagonal per row class; (3) implement
      the dominant term's implicitness; unit exactness tests
      as before. Exit: loop-gain census below unity on the
      stall cohort (implement, arm), lagged fraction irreducible
      within the splitting (close line, block-solve writeup),
      or user redirect.
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
      mom_solve 144 (incl. nurf/nuw MOM_NUSOLV and MOM_BLKGMRES
      refresh arms), mom_nu 53, mom_blk 59,
      mom_gs 34, mom_aa_state 62, mom_sweep_table 8,
      vadat_keys 91, es_resid 20, jfnk opint 36 + regressions,
      obs_frame 11+15, ng_abstain 17, lam_handoff 51, gs_drift_diag,
      iter_conv_gate 13. Targeted moment/VADAT batteries are green
      at 24-Nov-2026; a final full run_all_tests pass from the
      block-GMRES commit is still owed. test_es_conv remains blocked
      on the known obs/ archive gap.
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
