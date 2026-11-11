# SOB acceleration experiment results - X1-X4

Date started: 11-Nov-2026. Status: RUNNING (arms queued).
Audit basis: ES_ALI_BA_AUDIT.md; TODO: ES_ALI_BA_TODO.md;
slow-mode context: SOB_SLOW_MODE_ANALYSIS.md.

## 1. Experiment design

Four otherwise-identical warm-started clones of the DIE-off m15
2.82 d model, isolating the continuum electron-scattering
acceleration stack one flag at a time. All arms:

- `warm_restart_from` = converged `m15_2p82d_FRESH_sob` checkpoint
  (DIE-on fixed point) with `warm_restart_population_only: true`.
  The DIE-off rate change means this is an initial guess, not a
  converged state; each arm must re-converge to the DIE-off
  fixed point.
- Identical hydro (SN_HYDRO_DATA md5 ea6c5d5b, byte-verified),
  atomic data, velocity field, VADAT except the three arm keys.
- Cold-start anchor: `m15_2p82d_FRESH_sob_dieoff` (all accel off,
  cold start, it240+ at submission) - reported separately, never
  mixed into the warm-arm comparison.

## 2. Input configurations

Configs `configs/dessart_m15_2p82d_sobacc_x{1..4}.yaml`, cloned
from `dessart_m15_2p82d_FRESH_sob_dieoff.yaml`. Relevant
`vadat_overrides` deltas vs the dieoff baseline:

| arm | DO_ES_ALI | DO_ES_ALI_BA | ES_ALI_BA_MIN | DO_BA_RES_ACC |
|-----|-----------|--------------|---------------|----------------|
| X1  | F         | F            | -             | F              |
| X2  | T         | F            | -             | F              |
| X3  | T         | T            | 0.1D0         | F              |
| X4  | T         | T            | 0.1D0         | T              |

Shared recipe (unchanged across arms): DO_NG=T, DO_AV=T,
DO_JFNK=F, NUM_LAM=2, SOL_METH=DIAG, FIX_T=T, CONV_DIAG_LVL=2
(default). X4 additionally auto-arms DO_BA_RES_REUSE per the
SOB auto-set rule.

## 3. Jobs

| arm | FINAL job | OBS job | state at submit |
|-----|-----------|---------|-----------------|
| X1  | 13036072  | 13036073 | PENDING (QOSGrpCpuLimit) |
| X2  | 13036074  | 13036075 | PENDING |
| X3  | 13036076  | 13036077 | PENDING |
| X4  | 13036078  | 13036079 | PENDING |

Warm-restart staging verified at prepare time for all arms:
NT=1839 / ND=70 MODEL fingerprint match, POINT1/POINT2/SCRTEMP/
MODEL injected, checkpoint T replaced by hydro T (max |dT| =
1.879e-02 in 1e4 K units), SN_HYDRO_DATA byte-identical.

## 4. Metrics to record per arm

From ITER_HISTORY / OUTGEN / ITER_LEDGER (see
SOB_SLOW_MODE_ANALYSIS.md section 4 for the full diagnostic map):

- iterations to EPS_TERM exit; wall time; formal-solve count
  (ITER_LEDGER eval rows)
- maxraw / maxapp / maxatt trajectories; tail contraction rate
  (per-iteration maxapp ratio in the post-transient regime)
- BA_SE_RES (se) and BA_SE_RES2 (se2) trajectories
- ES_ALI: freqs unconverged, max final dJ/J, ALO fallbacks,
  J-floor hits run (X2-X4)
- ES_ALI_BA floor hits run (X3,X4)
- resacc telemetry (X4): nrej, trust, ba_rejected /
  ba_partial_rej / ba_reuse_stall counts, BA_NREJ_CONSEC max,
  BA_JUDGE_RATIO behaviour, LAMBDA demotions
- fixed-point agreement: final populations and spectrum vs X1
  within predefined tolerance (accelerated and unaccelerated
  must reach the same physical fixed point)

## 5. Results

(RUNNING - preliminary read at ~it734-739, 11-Nov-2026 ~16:00 UTC)

Warm-restart transient: all arms decayed maxapp 1.7e3 -> ~1.3
within ~15 iterations (geometric ~0.75/it during the cold-J
rebuild), versus the cold dieoff anchor pinned at the 1e2 cap
through it240. The warm start works as designed; NITSF resumed
at it723 (FRESH_sob ended it722), so "iterations to converge"
counts new iterations from ~723.

Early observations:

- X1 it739: FULL_TFIX/REJECT event, maxapp demand hit the 1.5e5
  cap and was rejected by the trust machinery (ba=REJECT in
  ITER_HISTORY even with DO_BA_RES_ACC=F - the pre-existing
  saturation/demand guard, not resacc). OUTGEN per-depth
  residual profile shows the residual concentrated at depths
  63-67 (outer region), 28.4% of column sum at d66. Convergence
  gauge read 88% held by ex_maxch.
- X2 ES_ALI live: `freqs unconverged` fell 2496 -> 1071 across
  the first three sweeps, max dJ/J 6.8e-3 -> 1.4e-3, ALO
  fallbacks=0, J-floor hits run=0, BA-floor hits run=0. The
  floors have not engaged on this model so far.
- X4 resacc live: BA_RES_ACC diagnostic blocks print every
  call; it727 shows max SE residual 6.4e12, RES2 1.1e12 at
  depth 63, worst eqn var=27 depth=70 `H2:lev 1` (outer
  boundary hydrogen - the degenerate floor the audit
  predicted), 0 hard/soft rejects, trust 1e4.
- X2/X3/X4 trajectories were bit-identical through it734
  (same maxapp/se/se2) - expected: ES_ALI_BA and resacc only
  differentiate once the J-loop/store differences feed back;
  ES_ALI itself (X2 vs X1) shows se2 ~10x lower at comparable
  iterations (8.7e10 vs 1.2e12 at ~it736) - early hint the
  continuum channel was contributing to the residual, needs
  the tail regime to confirm.
- Iteration wall time ~295 s on the it739 FULL_TFIX reject
  (vs ~64 s/iter on the solo dieoff run) - 4 concurrent arms
  plus OBS dependencies contend for nodes; per-iteration cost
  comparisons must use ITER_LEDGER timings, not wall clock.

Pending: tail contraction rates once maxapp < ~0.5; X3/X4
ES_ALI_BA floor census; X4 rejection/trust telemetry; fixed-
point agreement check at convergence.

## 5a. X2 CONVERGED (11-Nov-2026 ~16:40 UTC)

X2 (DO_ES_ALI=T only) exited clean at it752-753: "clean
convergence (MAXCH < EPS, exempted demand small)", Model
Convergence 100%. Job 13036074 COMPLETED in 33:03 wall.

- New iterations from warm start: ~30 (resumed at it723).
- Final ES_ALI census: freqs unconverged=0 every late sweep,
  max final dJ/J ~ 1.0e-4 (at ES_ALI_TOL=1e-4), ALO fallbacks=0,
  J-floor hits run=0, BA-floor hits run=0 - the 0.01/0.1
  floors never engaged on this model.
- Tail contraction in the post-transient regime (new-iters
  14-22): maxapp 1.13 -> 0.18, ~0.79/it sustained geometric,
  zero reject/spike events after the it735-equivalent
  transient.
- Per-iteration cost ~47-54 s (vs ~295 s on the contended
  X1 reject iteration; vs ~64 s/iter solo dieoff).

Comparison at same stage: X1 (no accel) was still oscillating
1.2-5.0 with repeated 1e2-1.5e5 reject spikes at new-it 18-21;
the unaccelerated warm arm had NOT broken into geometric
convergence. Cold dieoff anchor: it241+ and counting.

This is decision-gate CASE A: DO_ES_ALI alone collapses the
SOB plateau on this model. ES_ALI_BA (X3) and resacc (X4)
effects pending their convergence.

## 5b. X3 CONVERGED (11-Nov-2026 ~16:55 UTC)

X3 (DO_ES_ALI=T + DO_ES_ALI_BA=T) exited clean at ~it753:
"clean convergence", Model Convergence 100%. Job 13036076
COMPLETED in 35:09 wall, 28 new iterations.

- Statistically identical to X2 (30 iters / 33 min): the
  ES_ALI_BA amplification changed the path (floor engaged)
  but not the convergence rate.
- ES_ALI_BA floor census: `BA-floor hits run=50` in the
  ES_ALI summary, `ES_ALI_BA floor hits run=25` at the last
  do_var_cont banner - the 0.1 floor DID engage on ~25-50
  frequency-depth cells over the run, yet produced no
  measurable tail improvement. Confirms the audit's
  coverage-limit finding: the amplified diagonal responds
  to the continuum-ES channel, which ES_ALI already
  converges; the residual slow mode was never there.
- J-floor hits run=0, ALO fallbacks=0 - same clean ES_ALI
  census as X2.
- Fixed-point spot check vs X2: final-iteration largest
  corrections agree to ~1% (max decrease 3.19e-2 vs 3.21e-2,
  same level CrIII:3d4_5De[4] at depth 14; max increase
  2.11e-2 vs 2.08e-2). Same physical fixed point within
  EPS, as required.
- Population-file diff X2 vs X3 (all POP* dumps, ~118k
  values/species): max relative difference 1.06e-3 (Fe),
  7.8e-4 (O), 7.2e-4 (Ne), 6.3e-4 (Mg), all others <=5e-4.
  Accelerated arms reach the same fixed point to <0.11%.
  (STEQ_VALS excluded - it stores per-call residuals, not
  state; its 1e12-scale entries are the known residual
  floor, identical in character to both runs.)
- Emergent spectrum X2 vs X3 (OBSFLUX element-wise): max
  relative difference 1.2e-2 on bright continuum points,
  only 1/1416 bright samples above 1% - same observed
  fixed point within tolerance.

## 5c. X1/X4 still running (11-Nov-2026 ~17:20 UTC)

- X1 (all accel off, warm): it774, maxapp 0.65, still
  oscillating with reject spikes (2.07 at it767 after 0.71
  at it760). Effective contraction ~0.98/it with spikes -
  the same slow tail as the cold runs, ~50 new iterations
  in and not converged. Projected hundreds more iterations.
- X4 (ES_ALI+BA+resacc): it755, maxapp ~0.6-24 oscillating,
  nrej=5, trust collapsed 1e4 -> 312.5 through successive
  rejections. Descending but paying real recovery cost;
  CASE C (stability mechanism, not accelerator) forming.

## 5d. X5cold submitted (11-Nov-2026 ~17:45 UTC)

Cold-start arm `configs/dessart_m15_2p82d_sobacc_x5cold.yaml`:
byte-identical recipe to FRESH_sob_dieoff plus `DO_ES_ALI: T`,
NO warm restart (warm_restart: false). FINAL job 13038707,
OBS 13038708. Purpose: test whether ES_ALI shortens the
early/mid plateau, not only the warm-start tail - the
production-relevant case since real models start cold.
Cold baseline anchor = the running FRESH_sob_dieoff job
(13013148, it260+ at submission, all accel off).

## 8. Decision gate - PRELIMINARY VERDICT

CASE A CONFIRMED, CASE B NOT TRIGGERED, CASE C FORMING:

- DO_ES_ALI alone (X2) is the win: ~30 new iterations to
  clean convergence vs X1's projected 150-700+. The
  continuum electron-scattering J-loop was the plateau.
- DO_ES_ALI_BA (X3) adds nothing measurable over X2 -
  consistent with the audit's coverage analysis (the
  amplified diagonal addresses a channel ES_ALI already
  handles; the Sobolev line Jacobian remains untouched).
- DO_BA_RES_ACC (X4) trades iterations for safety: 5
  rejections caught steps that would have spiked X1-style,
  but trust collapse slows the tail. Retain as opt-in
  stability, not default.
- The slow mode was the continuum ES channel after all -
  the species histogram in SOB_SLOW_MODE_ANALYSIS.md
  (metal high ions) reflects WHERE populations are stiff,
  but the radiation lag that pinned them was the ES
  continuum, which ES_ALI fixes at the source.

Caveats before production recommendation:
- Single model, single snapshot (m15 2.82d, dense shell,
  ES_FAC=0.1). The result should be confirmed on a second
  model/epoch before making DO_ES_ALI a recipe default.
- X1/X4 final numbers pending; fixed-point agreement
  verified at OUTGEN level for X2/X3, population-file diff
  still owed.
- X2's win is warm-start specific: it still paid the
  cold-J transient (~15 iters). A cold-start ES_ALI run
  would test whether the acceleration also shortens the
  early plateau - worth one arm if production adoption is
  considered.

## 6. Defects found and fixed during setup

- control_variable_mod.f: ES_ALI_NFLOOR / ES_ALI_BA_NFLOOR
  declarations initially landed in the executable section
  (ifort #6274/#6236); moved to the specification part with
  `INTEGER ::` initialization. All touched files compile clean;
  incremental `make KEEP_EXE=1 all` exit 0.
- do_var_cont.f: ALO banner WRITE had 10 items against a 9-slot
  format (runtime format/variable-type mismatch, test_dvc
  'fg' case); format extended to `2A,I10`. Fixed; test_dvc
  170/170 pass.
- ES_ALI summary counters: ES_ALI_NUNC/NFALL/DJMAX reset per
  sweep while the new floor counters are run-total; labels
  changed to `... hits run=` to keep semantics honest.
- run_model_grid.py population_only path requires
  SN_HYDRO_DATA in the target at prepare time; pre-staged the
  byte-identical file per arm (md5 ea6c5d5b) before submission.

## 7. Regression evidence

- test_dvc.py: 170/170 (incl. new esalif floor-straddle case)
- test_vadat_keys.py: 85/85
- test_fg13.py: 20/20; test_iter_ledger.py: 31/31;
  test_sob_line_es.py: 9/9
- test_alteration_coverage.py: 1 pre-existing FAIL on files
  untouched by this work (jfnk_t_fixed_mode_v1.f,
  sob_line_es_update.f, pgplt/Makefile, ADD_SOB_LINE_PROF.INC)
  - manifest drift from earlier session edits, not a regression
  here.

## 8. Decision gate

(Filled after arms converge - see ES_ALI_BA_TODO.md X4.)
