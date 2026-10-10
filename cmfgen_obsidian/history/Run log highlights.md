# Run log highlights

The experiments that shaped the current design, condensed from
[[convergence_issues]], [[outstanding_issues]] and `AGENTS.md`.
Full dated detail lives in those files - this is the orientation
layer. Job numbers are SLURM IDs on viper.

## Timeline

### 18-Sep-2026 - clean tree, first Knop port

- `cmfgen_jfnk` = `cmfgen_v7_ord_run` + V4/V6 JFNK port (`10272f6`).
  V6's solver forks deliberately not imported.
- `DO_NONMONO_KNOP` integrated; sanity on monotonic M15: chi_hat path
  confirmed, GS warmup forcing Lambda for first 4 its.
- TypeIIn 10d: sigma_min=-98, 39% of zones sigma < -1 - the real test
  case. First SOB attempt diverged at it9 (1e109 corrections);
  decision made to *diagnose, not damp*.
- Frozen-Jacobian plateau found on stock settings: `FIX_BA=15,
  N_FIX_BA=0` -> 175 iterations at MAXCH~0.357. Fix is config-level:
  `FIX_BA=1, N_FIX_BA=10`.
- Upstream `USE_J_REL=T`+`TRIDIAG` segfault root-caused to a zero
  pivot in `CMF_BLKBAND_V3` at depth 1 - upstream bug, not a dev
  regression.

### 19-Sep-2026 - the big divergence hunts

- `obs_frame_sub_v9` bidirectional frequency-hunt fix: the observed-
  frame integration assumed monotonic V; produced ~13% corrupted bins
  on TypeIIn. Boundary clamps later made OBS clean (93,906 bins, zero
  bad, Ha dominant - qualitatively Type IIn).
- **obs_fin vs OBSFLUX**: the broad blue emission in obs_fin_10 is
  REAL (approaching shocked shell); OBSFLUX/CV_BOUND was ~43-84x below
  the required bolometric flux while the GS coupling was unconverged -
  obs_fin is the trustworthy product.
- MOM-path fixes: degenerate-J guards (`K/J -> 1/3`), INEXT/JNEXT
  under-relaxation (0.5), `'ZERO FLUX'` typo. Still diverged - the
  driver was the **expansive lagged map**: mismatched explicit/implicit
  weights made J grow ~10x/great-iteration. Fixed by the PSI_NM
  equal-weight compressing coupling -> STABLE + CONTRACTING.
- RH78 non-local Sobolev coupling implemented, found physically
  correct but unable to cure the TypeIIn SOB divergence - the
  escape-probability formulation itself is anti-contractive there.

### 19-21-Sep - nested iteration recovered

- In-run GS-immaturity gate ported from dev6: force Lambda while the
  lagged field moves; frozen-population sweeps (later revised to
  forced-Lambda, matching V6 semantics).
- Per-ray/per-frequency refresh (`FG_KNOP_REFRESH_V1`,
  `MOM_KNOP_REFRESH_V1`): converge the lagged coupling inside each
  great iteration instead of waiting NCF iterations.
- Metric fixes: column-peak significance normalization, adaptive
  omega, bulk gate (~0.1% bad significant cells), `FREQ_IDX_ITER`
  off-by-one, single-precision NaN store, int32 overflows,
  MOM Jacobi-snapshot (descending live reads caused a GS cascade),
  MOM containment (1E40 cap, best-pass restore, divergence abort).
- Result: FG gate opened - refresh max residual 1.8e6 -> 1.08e-1,
  bad cells 1.3% -> 0.023%; population updates collapsed to physical
  scale. MOM stayed in forced-Lambda - its lagged map is intrinsically
  non-contractive on this model.

### 29-30-Sep - stall audit + trials

- Exit-side hysteresis on the compressing-depth set, bounded bracket
  searches, quarantine-age counters, `BA_NREJ_LIM`, `GS_STATE` +
  `KNOP_STALL_MAX` watchdog, `GS_FREEZE_LIM`.
- `MOM_TRIAL_MIN` marginal release: trial BA under resacc after 30
  frozen iterations; accepted trials refund `N_FROZEN` (the 12040589
  run died at the freeze limit while progressing - refund fixed that);
  `GS_TRIAL_LIM` bounds pulsing.
- PSIPREV positivity fix (signed -> `|GAM*GAM_FAC_PREV|` at all three
  sites); dead CHI_HAT loop removed.
- itest PASSED (job 12042808, 3 its clean); production resubmitted
  `dessart_m15_2p82d_KNOP_MOM_dbg_dif` (FINAL 12043642, 16 CPUs).

### 27-Sep / 02-03-Oct - m15 SOB + JFNK bring-up

- **SOB working setup on Dessart m15 2.82d** (job 12001805): cold LTE
  -> stable Lambda contraction; load-bearing fix was `LSTAR` =
  header luminosity (the inner flux anchor) + `DIF=T` diffusion
  boundary, matching the Dessart paper setup.
- JFNK bring-up commits: `10272f6` (coupled residual at the SOLVEBA
  site), `72a1d03` (adopt eval F(0) as baseline), `e764368` (log-x,
  PTC, eta floor, F(0) rework), `7c6a3c3` (uninit-SUCCESS cold-start
  crash), `162d39e` (`JFNK_COR_GATE` demand gate).
- **it11 first live JFNK**: admitted at |F|~1e13 - far outside the
  basin; FD-JVP flat ~1e32 at every eps; ~3.2 h wasted; Lambda
  fallback preserved contraction. The it12 wreck was traced to the
  uncapped *Lambda safety step*, not JFNK itself.
- `JFNK_COR_GATE=1e2` validated live (job 12068898): the same 1e7%
  demand routed through the resacc path applied 9e2 instead of 1e7;
  contraction preserved.

### 03-Oct - the JFNK audits and the three-phase design

- [[JFNK_SOLVER_AUDIT]]: cost asymmetry quantified (coupled
  eval ~250-400 s), JVP noise floor identified as THE structural
  obstacle for coupled JFNK, single-metric gate found insufficient.
- Result: tiered admission (frozen tier cheap + coupled tier gated),
  `JFNK_GATE_NREQ=2` hysteresis, `JFNK_MAX_EVALS=40` budget,
  `JFNK_EPS=1e-3` noise-aware FD step, `JFNK_DEEPDIAG=F`,
  `GS_REL_STABLE=2` - the `dessart_m15_2p82d_JFNK_hyb.yaml` config.
- Admission hardened further (`2df1059`, `312b04f`): gate on
  `GS_FREEZE` not `GS_FG_IMMATURE`, hits only on eligible calls,
  demand contraction required.
- Late matrix-audit fixes: consistent `exp(x)-1` in matvecs,
  `POPS_WORK/POPS_BASE` sync before PC rebuild, `jfnk_fx_fresh` only
  on F(0) success (`286150f`).
- **Order-of-operations audit** (`9956744`, `e8abf97`, `900d515`):
  `JFNK_DISABLED` was a no-op (Newton ran anyway under mutated
  LAMBDA flags); module-STEQ save/restore completed to cover
  `SE%STEQ_ADV` + `SE%T_EHB` and bracketed around every
  compute_steq/eval-sweep site (frozen mode now keeps the full
  outer STEQ for the LAMBDA fallback); coupled eval without
  `JFNK_SWEEP_CB` now fails safe; `ACCESS_F` cursor guard for the
  `EDD_CONT=F`+`ACCURATE=F` combination. See [[JFNK solver]].

### Current state (03-Oct-2026, after the order-of-ops fixes)

All running jobs were resubmitted on the fixed binary
(`cmfgen_jfnk` @ `7221de9`); the stale-binary jobs and the old
validation directories were deleted.

| Run | Dir | Job | Binary |
|-----|-----|-----|--------|
| `dessart_m15_2p82d_SOB_mask` | `cmfgen_validation_SOB/m15_2p82d_mask` | 12073333 (+OBS 12073334) | fixed |
| `dessart_m15_2p82d_KNOP_mask` | `cmfgen_validation/m15_2p82d_knop_mask` | 12073335 (+OBS 12073336) | fixed |
| `dessart_m15_2p82d_JFNK_opt` | `cmfgen_validation_JFNK/m15_2p82d_jfnk_opt` | 12073580 (+OBS 12073581) | fixed |
| `dessart_m15_2p82d_JFNK_hyb` | `cmfgen_validation_JFNK/m15_2p82d_jfnk_hyb` | 12073582 (+OBS 12073583) | fixed |
| `dessart_m15_2p82d_SOB_v2` | `cmfgen_validation_SOB/m15_2p82d_v2` | 12072559 (+dep 12072560) | pre-fix (no JFNK code exercised) |

- The JFNK_opt/hyb resubmits are the first live exercise of the
  save/restore + self-disable fixes; watch for a real LAMBDA
  fallback on stagnation instead of Newton running anyway.
- SOB run: raw correction demand (`lraw`) saturated ~1e7 while
  applied `MAXCH` stays ~2-30 - bounded corrections working, raw
  demand still extreme at the stiff inner depths.

### SOB_mask attempt 1 audit (05-Oct-2026, job 12084677, ~3432b53+)

Stopped by CMFGEN itself at it305: `BA_NREJ_LIM=20` exhausted.
Root cause per the audit: `BA_REJ_FRAC` masked retries plus
LAMBDA-follow-up soft rejects each count toward `BA_NREJ_CONSEC`
post-`f315844`, so one failed BA episode burns ~4-7 rejections
where v2 burned ~1 - v2 needed ~8 episodes to converge, the run
died after ~5 while still improving.  NG fired 18 times, vetoed
every time (`NUM_BAD_NG` ~51-56/70): the `SCRTEMP` secant history
is polluted by rollbacks so extrapolations go negative.  Demand
stayed 1e4-1e7 at depths 62-70 (CV/OV/H2/HeI) - inside the
`BA_ENTER_RAW=1e4` hold band, so ~95% of iterations were LAMBDA.
Run dir deleted 05-Oct-2026 ~09:22.

Resubmit 2 settings (`dessart_m15_2p82d_SOB_mask.yaml`):
`DO_NG=F`, `BA_NREJ_LIM=100`, `BA_SAT_RAW=1e4` (routes the
1e4-1e6 demand band into the `BA_DAMP_CAP=1.3` damped-full march
instead of instant demote).  The temporary `DO_NG=F` pipeline choice
was superseded by resubmit 3's `NG_CLEAN`/`NG_ABSTAIN` fix.

Resubmit 3 (supersedes resubmit 2, same day; source `481902b`):
the 4-agent NG audit showed v2 actually logged 14 NG successes /
39 cancels -
the earlier "NG never succeeds" was wrong for v2 - and mask-1's
0/18 traced to resacc-tainted `SCRTEMP` records plus the
all-or-nothing per-depth veto.  Source fix: `NG_CLEAN` fire gate
(defers until 4 clean committed records) and `NG_ABSTAIN=T`
(per-component hold for zero-delta/non-positive extrapolations,
vacuous attempt returns `NG_DONE=F`).  `DO_NG` back to T with
`NG_ABSTAIN=T` in the yaml and the non-mono pipeline defaults.

Optimized SOB profile (job 12090475 it90-131 trajectory):
`BA_SAT_RAW=2e3`, `NUM_LAM=BA_PROBE_LAM=20`,
`LAM_N_STABLE=3`, `FIX_BA=1.0`, `N_FIX_BA=10`, `FIX_IMP=T`,
`BA_ENTER_RAW=1e4`, `BA_DAMP_CAP=1.3`, `BA_REJ_FRAC=0.25`,
`BA_NREJ_LIM=100`, `DO_NG=T`, `NG_ABSTAIN=T`, `BEG_NG=5`,
`IBEG_NG=15`, `BW_NG=40`, `ITS/NG=8`.  The profile is now stamped by
`NON_MONO_SOB_RECOMMENDED` in `prepare_VADAT.py` and imported by
`run_model_grid.py`, so future SOB submissions get it without
repeating YAML overrides.  The same `BEG_NG=5`, `IBEG_NG=15`, and
`ITS/NG=8` schedule is now the generated-run default for the other
non-mono recipes.

## Reusable lessons

- **Diagnose, don't damp.** Every attempted clip/damp that hid a
  divergence delayed finding its cause (RH78, LAM_SCALE_OPT, the
  15-iter freeze cap).
- The lagged store makes every restart a physics question: a process
  restart re-zeroes it unless `FGINEXT_GS`/`MOMJNXT_GS` restore.
- "Converged-looking" single sweeps are coincidences; hysteresis
  counters exist because single calls lie.
- Demand (what the solve *wants*) and applied step (what it *did*)
  are different signals; watch both.
- The eval cost asymmetry (~100x) makes coupled JFNK a budget-managed
  resource, not just an algorithm.

## Related

- [[Open issues]]
- [[Failure modes and guards]]
- [[Literature]]
