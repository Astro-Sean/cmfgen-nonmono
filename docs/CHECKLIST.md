# CHECKLIST - cmfgen_jfnk development ledger

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
