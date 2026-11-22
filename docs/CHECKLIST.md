# CHECKLIST - cmfgen_jfnk development ledger

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
      informative. Verdict and job id recorded here.
      Exit criterion: named root cause + FIX commit with
      regression, or an evidenced BLOCKED verdict. No
      tolerance relaxation, no new clamps.
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
- [ ] BLANK spectrum + Luc comparison: renderer verified 22-Nov
      (RESC_JMAX fix; obs_fin_10 valid). Blocked only on a
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
- [ ] sob_diag model run to convergence/spectrum - submitted
      13474708/13474709 on the fixed build, pending on `small`.
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
      Kept documented, not selected.
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

- [x] Static/unit batteries: knop_ray 43, lagged_field suite,
      mom_solve 137, mom_gs 34, jfnk opint 36 + regressions,
      obs_frame 11+15, ng_abstain 17, lam_handoff 51, gs_drift_diag,
      iter_conv_gate 13 - all green at last run.
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
