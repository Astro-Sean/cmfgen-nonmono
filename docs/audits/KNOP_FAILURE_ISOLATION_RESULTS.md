# KNOP Production Failure Isolation - Results

Status: IN PROGRESS - diagnostic arms launched, runtime results pending.

This document records the failure-isolation experiment defined by the
revised work order (Phases 1-7). Sections marked PENDING RUNTIME await
the Slurm arms described below; all source-level facts, freeze
semantics, and configuration are verified and final.

## 1. Experimental configuration

Four warm-restart arms from the `m15_2p82d_FRESH_knop` iteration-125
checkpoint (the run that stopped non-convergent via
`GS_REL_NOPROG_LIM` after ten no-progress population releases). All
arms run single-threaded (`cpus_per_task: 1`, so
`OMP_NUM_THREADS=1`/`KN_MAX_THREADS` effective 1) per the KN_HOT_SPENT
scheduling-race caveat. `init_num_its: 40` per Slurm invocation with
the standard wall requeue.

| Arm | VADAT keys | I_NEXT seed | Job |
|-----|-----------|-------------|-----|
| R   | `KN_FREEZE_COEF=T`, `KN_MAT_FREEZE=T`, `KN_Q_RESID=T` | it125 FGINEXT_GS | 13065894 |
| R2  | same as R | it125 FGINEXT_GS with I_P/I_M scaled by 1.001 | 13065898 |
| M   | `KN_MAT_FREEZE=T`, `KN_Q_RESID=T` | it125 FGINEXT_GS | 13065896 |
| C   | `KN_Q_RESID=T` | it125 FGINEXT_GS | 13065911 |

Configs: `configs/dessart_m15_2p82d_KNOP_iso_{R,R2,M,C}.yaml`.
Resolved configs and per-run `config_resolved.yaml` live under each
`cmfgen_validation/m15_2p82d_KNOP_iso_*/` directory.

### Reproducibility metadata

- Source commit: `11ef75c63a75818dbcfbc6a5986d017f2571c342`
  (diagnostic commits `a5b4164a`, `355f23cf`, `777b0e36`, `3b8f7959`,
  `11ef75c6` on top of `3dde64b6`).
- Executable: `exe/cmfgen_dev.exe`, md5
  `8afd859176ca0177c9de31fc6a400ffc`.
- Compiler: Intel Fortran Classic 2021.11.1 (build 20231117),
  `/mpcdf/soft/RHEL_9/packages/x86_64/intel/2024.0.2/bin/ifort`.
- Flags: `-qopenmp -mkl -nbs -diag-disable=10448,10412,10434
  -extend-source -traceback -O3 -mavx512f -mavx512dq -mavx512bw
  -mavx512vl -mavx512cd -align array64byte` (Makefile_definitions).
- Checkpoint provenance: `cmfgen_validation/m15_2p82d_FRESH_knop/
  cmfgen/2.82days` - OUTGEN records "Current great iteration count is
  125", the non-convergent stop. `MODEL` md5
  `e4fad0050aede72ccd157e29dfbaeb53`, `VADAT_INIT` md5
  `30b528c9337cc821cab41d968eb4a94c`.
- Seed-2 checkpoint: `cmfgen_validation/_ckpt_iso_seed2` - identical
  SCRTEMP/POINT1/POINT2/MODEL/CUR_MODEL_DATA; FGINEXT_GS recs 9-10
  (I_P_DEPTH_NEXT, I_M_DEPTH_NEXT; 335,443,150 values each) scaled by
  1.001, all other records bit-identical (verified: same file size,
  identical first 100 bytes, tool `perturb_fginext`).
- Bit-identical common state: arms R/R2/M/C share the same checkpoint
  files except R2's intentional I_NEXT perturbation; VADAT, MODEL,
  atomic data and solver configuration are identical apart from the
  listed diagnostic keys. Process-local SAVE variables start at their
  DATA initializers in every arm (fresh process per Slurm job), so no
  warm-restart state leaks between arms.
- Fingerprints recorded at runtime: `KN_COEF_FP` (raw KN_CHI_EXT/
  KN_ETA_EXT snapshot sums at the freeze-arm event), `KN_COEF_FPR`
  (rebuilt KN_CHI_COEF/KN_ETA_COEF sums per armed refresh),
  `KN_COEF_FP_DRIFT` (max drift vs first armed refresh), plus the
  march census and per-ray `KN_Q_LS` residuals in OUTGEN.

## 2. Checkpoint provenance

The source checkpoint is the terminal state of
`configs/dessart_m15_2p82d_FRESH_knop.yaml` (cold LTE start,
`LTE_EST=T`, no warm restart, `KN_NO_GAIN=T`, `DO_ES_ALI=T`,
`DO_ES_ALI_BA=T`, `DO_JFNK=F`, `FIX_T=T`, `NUM_LAM=2`). Its OUTGEN
ends with "the release/refreeze cycle never matured the lagged field;
stopping as non-convergent" after 125 great iterations. FGINEXT_GS is
5,367,543,560 bytes (N_COMP=70, NP=85, NCF=56119 per the run's
memory estimate; I_NEXT stores dominate).

`warm_restart_from` staged POINT1, POINT2, SCRTEMP, FGINEXT_GS,
MODEL, CUR_MODEL_DATA into each arm's run directory (verified in the
R submission log: "injecting checkpoint ... (NT=1839, ND=70)").

## 3. Actual coefficient-freeze semantics (verified from source)

`KN_FREEZE_COEF` ([KN_FREEZE_COEF], default F, setter
`SET_FG_FREEZE_COEF`):

- During the frequency sweep of `FG_J_CMF_NONMONO_V1`, each slot's
  `CHI_EXT`/`ETA_EXT` are copied into `KN_CHI_EXT`/`KN_ETA_EXT` only
  while `KN_COEF_FROZEN` is false (plane/fg_j_cmf_nonmono_v1.f,
  capture block ~line 3200). At the last frequency slot
  (`LK = KN_NCF`), `KN_STORE_READY` arms and, when the flag is set,
  `KN_COEF_FROZEN` latches true and `KN_COEF_FP` fingerprints the
  completed snapshots.
- All subsequent calls skip the snapshot copy, so every
  `FG_KNOP_REFRESH_V1` rebuilds `KN_CHI_COEF`/`KN_ETA_COEF` (the
  monotone-cubic coefficients consumed by `BUILD_FREQ` -> `CHI_L`,
  `ETA_L` -> `Q_L`, `SR_L`) from identical inputs. `SR_L`/`Q_L` are
  never independently overwritten; they remain derived from the
  frozen `CHI_L`/`ETA_L` and the unchanged velocity geometry
  (`GAM`, `GAM_MAG`, `dGAMdR`, ray table).
- The freeze covers the *refresh* coefficient path only. The main
  per-frequency formal solve still consumes the live `CHI_EXT`/
  `ETA_EXT` built fresh each call; that is the production design -
  the diagnostic targets the lagged I_NEXT store-to-store map
  `T`, not the per-frequency Lambda sweep.
- `KN_COEF_FROZEN` is a module SAVE with DATA initializer .FALSE.;
  every Slurm arm is a fresh process, so the freeze always arms on
  that arm's first complete frequency pass (the it125 matter state).

Scope boundaries (audit clarification):

- Arm R freezes exactly one thing: the `KN_CHI_EXT`/`KN_ETA_EXT`
  snapshot pair that `FG_KNOP_REFRESH_V1` uses to rebuild
  `KN_CHI_COEF`/`KN_ETA_COEF` for the lagged I_NEXT refresh
  operator.
- Arm R does **not** freeze the CMFGEN radiation iteration. The
  main per-frequency formal solve keeps its own live coefficient
  chain: `CHI_EXT`/`ETA_EXT` are recomputed on every
  `FG_J_CMF_NONMONO_V1` entry from the current matter state and
  `TA` (which still carries the live `CHI_SCAT_CLUMP*RJ[_ES]`
  term), then `CHI_COEF`/`ETA_COEF` are rebuilt from them for the
  frequency-by-frequency sweep. The freeze touches only the
  refresh-side `KN_*` arrays.
- The derived quantities `CHI_L`, `ETA_L`, `SR_L`, `Q_L` used by
  `BUILD_FREQ` are per-frequency/per-ray working values rebuilt
  inside each refresh call; under the freeze they are regenerated
  from the frozen snapshots every refresh, so they are constant
  across refreshes but not themselves stored/frozen state.
- Consequence: if arm R stalls, the defect is *not* localized to
  the fixed-coefficient refresh operator or the production ray
  geometry. The live main-solve coefficients still evolve with
  pops and RJ, so a stall is equally consistent with residual
  coupling through the live path. Only convergence of R (with C
  stalling) would support a population/matter-feedback cause; a
  stall in all arms remains INCONCLUSIVE between refresh-operator
  and live-path causes until Phase 4 geometry work separates them.

`KN_MAT_FREEZE` ([KN_MAT_FREEZE], default F): in `SOLVE_FOR_POPS`,
the hold condition gains `.OR. KN_MAT_FREEZE`, forcing the existing
`GS_POPS_HELD` path on every call. `SOLVEBA_V13` still executes (so
STEQ/demand diagnostics stay live) but the committed `POPS`
correction is rolled back. `POPS(NT,ND)` carries all statistical
unknowns including the electron-density row (NT-1) and temperature
row (NT); `POP_ATOM` is read-only inside `SOLVEBA_V13`. With pops
held, every matter-derived opacity input (`CHI`, `ETA`, `ESEC`,
clumps) stays at checkpoint values while the radiation field and
the J-dependent source iterate evolve.

`KN_Q_RESID` ([KN_Q_RESID], default F): per-ray, per-channel
independent fixed-point residual. After each channel's converged
sweep inside `KNOP_RAY_SOLVE_V1`, the code copies `FLD`/`FLD_P` to
scratch `FLD_Q`, applies exactly one effective sweep-map pass
(`SWEEP_FLD_ORD` when `ORD_SWEEP` else `SWEEP_FLD`, with a
diagnostic `.TRUE.` argument that suppresses commit/census side
effects inside the sweep), records `Q_M`/`Q_P` into `KN_Q_LS(LS)`,
and restores all shared sweep state (divergence, exit, worst-cell,
census, cap, negative-cell, p2-damping, rho exports). The extra map
application never touches production stores, commit counters, or
field diagnostics.

What `RES_Q` actually measures (source-traced, not inferred):

- `MARCH` computes, for every updated cell,
  `DCELL = ABS(VNEW-VOLD) / MAX(VOLD, 1.0E-4*VMAX_L)` where `VNEW`
  is the raw short-characteristics candidate and `VOLD` the cell's
  pre-update value; `VMAX_L` is the column peak floored at 1e-25.
  The pass residual `RES` is `MAX(RES, DCELL)` over all cells.
  `DCELL` is recorded **before** the relaxation overwrite
  `FLD = VOLD + OMEGA_RAY*(VNEW-VOLD)`, so it is the unrelaxed
  map defect, not the committed step.
- Both sweep routines set `OMEGA_RAY = KN_RAY_OMEGA` once at entry
  (`KN_RAY_OMEGA` is a module constant initialized to 1.0_LDP with
  no VADAT override). The `L_QP` exit fires at pass level
  (SWEEP_FLD ~line 8646, SWEEP_FLD_ORD ~line 9665) **before** the
  same-direction residual-ratio and sign-flip adaptation block
  (~lines 8688-8733 / ~9725-9738), so the verification pass always
  marches with `OMEGA_RAY = 1.0` for every ray and cell. With
  omega = 1 the committed update equals `VNEW`, i.e. the
  verification pass is the pure unrelaxed sweep map.
- Therefore `RES_Q` is option A of the audit question: a raw
  one-pass fixed-point defect `max_cells |T(I)-I| / norm_floor`,
  expressed in the solver's relative-change metric. It is **not**
  a relaxed residual and **not** a vector norm; the denominator is
  `max(VOLD, 1e-4*column_peak)`, a per-cell relative floor, so small
  cells are measured against the column scale rather than their own
  value.
- Because the sweeps update `FLD` in place, the pass is a
  Gauss-Seidel sweep map, not a Jacobi map: later rays see earlier
  rays' committed updates. `RES_Q` is the defect of one ordered
  sweep pass from the per-ray-converged state, sampled at pass 1.
- Zero residual still identifies exactly the same fixed point:
  `DCELL = 0` for all cells means `VNEW = VOLD`, i.e. `T(I) = I`;
  relaxation cannot move the fixed-point location for any
  `0 < omega <= 1`.
- Comparability: the production maturity residual uses the same
  `DCELL` formula, but it is sampled on the **last** pass after N
  passes of contraction and `OMEGA_RAY` adaptation (omega can fall
  to the 0.45 floor between passes). `RES_Q` samples pass 1 from the
  per-ray-converged state. A production pass-1 residual and `RES_Q`
  are directly comparable only if the starting fields match;
  `RES_Q` must not be compared magnitude-for-magnitude with the
  production terminal residual without accounting for pass number,
  omega state, and sweep direction.
- Latent hazard: `RES_Q`'s "raw defect" property depends on
  `KN_RAY_OMEGA` remaining 1.0. If that constant is ever changed,
  the verification pass silently becomes a relaxed residual
  (omega applied inside `MARCH` before subsequent rays read `FLD`)
  while the recorded `DCELL` stays pre-relaxation - the metric
  would then mix raw defects with a damped committed trajectory.
  Documented here; no code guard added per the audit-only mandate.
- Optional cheaper raw-defect supplement (proposed, not
  implemented): `MARCH` already holds `VNEW`/`VOLD`; accumulating
  `MAXVAL(ABS(VNEW-VOLD))` and an `L2` norm alongside `DCELL`
  would add an unnormalized defect scale at ~zero cost. Not added
  because the relative metric suffices for pass/fail diagnosis and
  any source change requires a rebuild the running arms must not
  pick up.

### Coefficient fingerprint strength (audit clarification)

The armed-refresh fingerprint has two tiers:

- `KN_COEF_FP`/`KN_COEF_FPR`/`KN_COEF_FPB`/`KN_COEF_FP_DRIFT`:
  global sums of the raw `KN_CHI_EXT`/`KN_ETA_EXT` snapshots (at
  arm) and of the rebuilt `KN_CHI_COEF`/`KN_ETA_COEF` (each armed
  refresh), drift = max |current sum - first-refresh sum|. Global
  sums are weak evidence: compensating local changes (one element
  up, another down) or slot permutations leave the sum unchanged,
  so `drift = 0` does not prove the coefficient arrays are
  element-wise identical across refreshes. They do prove the
  *aggregate* opacity/emissivity budget is stable, which is still
  useful for detecting gross leaks.
- Element-wise checksum, two revisions:
  - `f67034a1` (post-arm-launch, superseded): a **custom
    word-wise** checksum mislabeled as FNV-1a - it XORed each
    whole 64-bit `TRANSFER` image as one word against a truncated
    offset basis (`1469598103934665603`, not the spec value,
    because the true basis exceeds signed INTEGER*8 range). It
    detected local change but was not the published algorithm.
  - Corrected standard 64-bit FNV-1a (this commit cycle): module
    procedures `KN_FNV1A_I1`/`KN_FNV1A_COEF`/`KN_FNV1A_MIX8` in
    `plane/fg_j_cmf_nonmono_v1.f`. Offset basis
    `0xCBF29CE484222325` via BOZ literal (spec value, bit pattern
    preserved as signed int64), prime `0x100000001B3`, and every
    byte of each `REAL(LDP)` element mixed least-significant-byte
    first after `TRANSFER` to INTEGER*8 (STORAGE_SIZE=64 verified
    in the test driver). Array extents are mixed first so the
    hashed sequence is unambiguous; CHI and ETA are hashed through
    separate calls. O(1) auxiliary memory - no second
    full-size coefficient copy. Verified by
    `unit_tests/f90/drv_fnv1a.f90` + `test_fnv1a.py`: published
    FNV-1a vectors ("", "a", "foobar") pass, all 8 byte positions
    detected, compensating-pair change detected with bitwise-
    identical global sum, 12/12 driver checks.
- Integer arithmetic caveat: the recurrence needs mod-2^64
  arithmetic. Fortran leaves signed overflow processor-dependent;
  ifort implements two's-complement wrap for INTEGER*8, which is
  exactly FNV-1a's arithmetic, and the production flags carry no
  integer-overflow trap. This is compiler-specific reliance -
  stated in the source comments and verified empirically by the
  spec vectors under the production flags; a compiler change in
  overflow semantics would fail the vectors loudly.
- Collision limits: FNV-1a is a non-cryptographic checksum. A
  matching checksum is high-confidence evidence of bitwise
  equality (~2^-64 per comparison for an ideal hash; FNV-1a is
  weaker against adversarial inputs, though solver-produced
  coefficient drift is not adversarial) - it is a diagnostic, not
  a formal proof. Element-wise max-abs/max-rel diff with index
  reporting remains the fallback if a mismatch ever fires
  (offline diff of the frozen snapshots; no permanent baseline
  copy allocated).
- IMPORTANT: the running arms R/R2/M/C execute the executable
  built before either checksum existed (md5 8afd8591...). Their
  OUTGEN will contain only the sum-based fingerprints. Checksum
  lines appear in the refresh summary only for runs built after
  this commit; treat "checksum absent" in arm output as expected,
  not as a missing diagnostic. `exe/cmfgen_dev.exe` was restored
  bit-identical after the validation build for exactly this
  reason; `exe/cmf_flux.exe` was relinked by the same build
  (md5 now 526d6301) - its executed code paths are unchanged by
  this diagnostic (the flux tool never enters the KNOP refresh),
  but the pinned-md5 provenance no longer matches for OBS-stage
  invocations.

### Upstream ETA construction (Phase 2 trace)

`comp_j_blank.f` passes `TA` as the solver's ETA argument:

- `TA = ETA_CLUMP + CHI_SCAT_CLUMP*RJ` under `COHERENT_ES`, else
  `...*RJ_ES` (comp_j_blank.f ~line 841). `ETA_CLUMP = ETA*CLUMP_FAC`
  is the thermal+bound-bound emissivity assembled in `cmfgen_sub.f`
  (continuum + `ETAL_MAT` line terms + `BA_EHB_BF_UPDATE_V1`),
  population-dependent and J-independent.
- `CHI_SCAT_CLUMP*RJ[_ES]` is the coherent electron-scattering
  source: explicitly J-dependent through the Lambda/ES_ALI iterate
  `RJ`/`RJ_ES`.
- Inside `FG_J_CMF_NONMONO_V1`, `ETA_EXT` is `TA` extended to the
  outer grid; the file's own comment states "ETA does however change
  (since it depends of J) and thus ETA_EXT must be re-computed on
  each entry" (~line 3126).

Consequence for arm interpretation:

- Arm R freezes `KN_ETA_EXT` snapshots, so both the population part
  and the J-dependent scattering part of ETA are pinned at their
  it125 values inside the refresh operator: the purest
  fixed-coefficient transport test.
- Arm M holds pops (so `ETA_CLUMP`, `CHI`, `ESEC` are static) but
  leaves `CHI_SCAT_CLUMP*RJ` live: any radiation-dependent
  scattering feedback still reaches the rebuilt coefficients. If M
  stalls where R converges, the J-dependent scattering channel is
  implicated; if both converge, population-opacity feedback (arm C
  vs M/R) is implicated.
- No transfer equations were modified; the scattering term was
  already separable upstream (`ZETA + THETA*RJ` structure), so the
  experiment isolates it without touching the formal solver.

## 4. Ordered-march branch census (Phase 5 instrumentation)

`FG_ORD_SWEEP` ([KN_RAY_ORDERED], default T), `HAVE_PACK`, `NI>=3`,
and the `SPW` allocation fallback jointly determine `ORD_SWEEP` per
ray inside `KNOP_RAY_SOLVE_V1`. New `INTEGER*8` module counters
(`KN_CENS_ORD_RAY`, `KN_CENS_UNO_RAY`, `KN_CENS_SFM`,
`KN_CENS_SFOM`, `KN_CENS_SFP`, `KN_CENS_SFOP`) are zeroed once per
`FG_KNOP_REFRESH_V1` before the OpenMP ray loop and incremented
under `!$OMP ATOMIC` at each actual `SWEEP_FLD`/`SWEEP_FLD_ORD`
call, distinguishing M and P channels. The refresh summary prints:

```
    march census: ord rays=...  unordered rays=...
    solves M ord, M uno, P ord, P uno=...
```

This records executed branches, not the input flag. Verified in
`unit_tests` (`drv_knop_ray rfsh` prints the census; all 43
`test_knop_ray` checks pass, `test_knop_lagged_field` and
`test_knop_status_text` pass after updating the dispatch-site pin
for the two new diagnostic call sites).

PENDING RUNTIME: production-geometry census values from arms R/M/C.

## 5. Results for R, M and C - IN PROGRESS (refreshes 1-3, its 125-127)

Extraction tooling: `scripts/parse_knop_iso.py` (workspace, not in
this repo) parses each arm's OUTGEN into refresh-aligned records
keyed by refresh sequence with the enclosing great iteration;
ITER_SUM rows capture mode/frozen/held, pops demand/applied/status,
se residual, lum diag, fg_rf row, conv-exit gates, model
convergence and wall time.

### Interim observations (its 125-127, all arms pops-held)

- ALL four arms show `frozen=T held=T` and `pops [HELD]` with
  `pops_applied=0.0` in ITER_SUM. This is inherited checkpoint
  state, not the diagnostic flags: the parent FRESH_knop run held
  pops through it123, released at it124 (`frozen=F held=F`,
  applied=38.43), and re-froze at it125; `GS_STATE restart`
  (solve_for_pops.f ~line 741) adopts `n_frozen`/`frozen` from the
  checkpoint. Arm C is therefore NOT yet a clean coupled control -
  its populations follow the normal release cycle and will update
  only on release iterations, while arms M/R/R2 (KN_MAT_FREEZE)
  hold every call.
- Consequence for interpretation: until arm C's first pop-release
  iteration, M and C are expected to be bit-identical (confirmed:
  every refresh metric matches through refresh 3). The M-vs-C
  comparison becomes meaningful only after C releases pops.
- R-vs-M IS already informative: with pops held in both, the only
  difference is refresh-side coefficient evolution (frozen
  snapshots vs live rebuilds carrying the J-dependent
  `CHI_SCAT*RJ` term). Early trend favors R: rf_resid
  3.317e-2 -> 3.044e-2 (R) vs 3.317e-2 -> 3.109e-2 (M/C) by
  refresh 3; unconverged cells 24143 (R) vs 19760 (M/C) - mixed,
  not yet decisive.
- R-vs-R2 seed sensitivity: the 0.1% I_NEXT perturbation shows at
  refresh 1 (rf_resid 3.406e-2 vs 3.317e-2, worst_ray 5 vs 7) and
  decays to near-identity by refresh 3 (3.098e-2 vs 3.044e-2).
  One seed direction only; does not prove seed independence.
- Q residual ~3.09 in all arms at refresh 1, ~3.09 at refresh 3 -
  essentially static while production rf_resid falls ~10%. Exit
  census: exact=48, tol=18-20, bulk=14-18, zero stall/cap/div;
  1 unsolved ray throughout; worst ray passes=6-7 of max 12.
  bulk_floor=6.0e-2 vs bulk maxres=3.3e-2 - bulk exits are below
  their residual floor.
- Production refresh residual (KN_RF_RESID, the `fg_rf` row) is
  ~3.0-3.4e-2 vs tol 2.0e-3 - ~16x tolerance, [WARN] maturity
  status, conv-exit gates 50% (worst margin ~2.9 dex, ex_maxch).
- Per-refresh wall time ~480-565 s; per-iteration ~700-1000 s.
  40 requested its projects to ~10-12 h plus walltime requeues.
- Executable for all arms: cmfgen_dev.exe md5
  8afd859176ca0177c9de31fc6a400ffc (predates the corrected
  FNV-1a checksum; only sum fingerprints appear - expected).
  `coef fp ... rebuilt drift=0.0 armed refreshes=1` in arm R
  confirms the freeze armed once and sums stayed identical
  (sum-level evidence only; element-level hash unavailable in
  this exe).
- March census: ord rays=84, unordered=0 in every refresh of every
  arm - the ordered march covers the whole production ray set.

### Refreshes 4-6 (its 128-130)

- M and C remain BIT-IDENTICAL through all six refreshes; arm C is
  still `frozen=T held=T` at it130, so the M/C differentiator has
  not yet switched on. The population release cycle inherited
  from the parent run has not recurred inside the arm window.
- R-vs-R2: the 0.1% seed difference has decayed to the printed
  precision by refresh 4-6 (rf_resid 3.099e-2 both; unconverged
  19830 vs 15772 at seq4, 16464/18345 vs ~18175 at seq5-6 -
  residual-level differences of a few percent persist).
- rf_resid is NOT monotone: R dips to 2.658e-2 at seq5 then
  returns to 3.099e-2 at seq6; M/C sit at 3.109-3.110e-2 for
  seq2-6. All arms remain ~15x the 2.0e-3 tolerance.
- Q residual stays ~3.09 in every arm and refresh - static while
  rf_resid moves. Exit census drifts toward exact/tol exits
  (bulk falls 18->9-16); a single stall exit appeared in M/C at
  seq6; no cap/div exits anywhere.
- Unconverged-cell counts keep falling in all arms (R 16464,
  M/C 5258 at seq6) while the max residual stays pinned - the
  residual is controlled by a shrinking set of cells/rays.
- cell50 (cells above 50% relative change AND above the
  1e-4*peak significance floor) is 0 for R seq6 and M/C
  seq4-6. NOTE: this is the GS_NCONV churn census, not a
  floor-dominance test - the floor criterion is
  VOLD < 1e-4*VMAXF(L), a much weaker threshold. cell50=0
  only shows no LARGE relative moves persist; it does not by
  itself prove the pinned residual is floor-dominated. The
  KN_QA_* argmax output (VOLD vs 1e-4*VMAXF) is required for
  that determination and is not available in this exe.

### Refreshes 7+ (its 131+, continuing)

- C remains `frozen=T held=T` through it130 (latest complete
  ITER_SUM at time of writing); the population-release gate has
  NOT opened inside the arm window so far. M/C remain
  bit-identical; the M-vs-C comparison stays inconclusive until
  C's first `held=F`/`pops_applied>0` iteration.

### GS_MAXCHG argmax trajectories (Phase 4, from existing output)

The `argmax history (d,r,f)=v` line (CONV_DIAG_LVL>=2, already
active) records the sweep-store drift argmax per great iteration
with a 4-deep history and a cycle/margnl/bndry classifier
(`FG_GS_ARGK/ARGR/ARGL`, `FG_ARG_H*`, `FG_GS_CLS`). Parsed
trajectories:

- R and R2 converge onto ONE persistent cell: (d=54, r=30,
  f=8712) from it127/128 onward, classified `cycle` (period-2
  limit cycle at a fixed cell), value oscillating
  3.04->3.10->2.66->3.10->3.10 e-2. The same cell controls the
  maturity residual for 5+ consecutive iterations.
- M and C (identical while pops held) show a WANDERING argmax:
  (51,32,1175) -> (8,76,1372) -> (3,82,1249) -> (5,80,1249) ->
  (4,81,1249) pinned at 3.110e-2. No persistent cell; the max
  drift moves among shallow depths (d=3-5), ray ~80-82,
  frequency 1249.
- Arm differentiation: frozen refresh-side coefficients (R/R2)
  condense the drift onto a single deep-cell 2-cycle; live
  coefficient snapshots (M/C) produce a moving, non-localized
  argmax plateaued at the same ~3.1e-2 magnitude. This is the
  first clean R-vs-M behavioral split in the experiment.
- The argmax sweep-drift value (~3.1e-2) matches rf_resid
  (~3.1e-2) in magnitude, consistent with both metrics reading
  the same underlying store drift.

### Commit-path audit (Phase 4, source + output)

Commit machinery (`plane/fg_j_cmf_nonmono_v1.f`):

- `WRITE_STORE` commits converged-ray fields: per cell
  V1=fresh FLD(JP2,LS2), V2=retained STORE(KC,LS,LS2-1),
  DLOC=|V1-V2|/MAX(V2,1e-4*VMAX_F(LS2)); unrelaxed write
  (INEXT_RELAX stays on sweep sites only, W4 17-Oct-2026).
- `WRITE_STORE_CAP` commits non-converged rays under
  FG_COMMIT_CONV=T with each cell's change bounded by
  FG_COMMIT_NC x retained (seed fills capped at NC x the
  significance floor).
- `FG_COMMIT_JUMP=50` rolls back a whole ray whose committed
  store change exceeds the bound; `FG_RAY_VC` records the
  per-ray verdict (0 committed, 1 jump-rej, 2 nc-rej, 3 capped,
  4 div, -1 unsolved) and `FG_RAY_DCH` the per-ray delta.
- Per refresh the arms show 2-4 jump-rejections + ~1 capped +
  1 unsolved ray of 85; committed ~80-82. The outer residual is
  therefore NOT produced by wholesale commit failure - most
  rays commit every refresh.

MISSING QUANTITIES (existing output does not expose; no new
instrumentation added per instructions):

- `FG_RAY_VC`/`FG_RAY_DCH` for the ARGMAX ray specifically
  (r=30 in R/R2, r~80 in M/C): only aggregate verdict counts and
  the RES_LS-argmax ray's (worst_ray=5-12) values print.
- Per-cell proposed V1 vs retained V2 and relative change at
  the argmax cell (d=54,r=30,f=8712 -> store slot f-1=8711):
  computed inside WRITE_STORE/_CAP but never printed per cell.
  The KNSTOR debug line (st/sv/rel) is gated on
  PRINT_KNOP_DEBUG and a >0.25 relchg threshold that the
  ~0.031 argmax would not meet anyway.
- The argmax ray's exit class and pass count (per-ray exit
  class is not stored; only the aggregate census prints).

### R-vs-R2 seed sensitivity (Phase 5, quantified)

Refresh-aligned |R - R2| over 7 refreshes:

- rf_resid: 9e-4 down to 5e-6, rising again to ~1.3e-3 by
  seq7 (~0.02%-14% of the ~0.031 level) - the perturbation has
  NOT fully decayed; residual differences fluctuate.
- q_max: ~1.3e-3 absolute (~0.04% of 3.09), one refresh 2.8e-4.
- retained_max: 3e-4 to 1.28 absolute - seq7 shows an O(1)
  store-delta difference between seeds.
- fp_drift=0.0 both arms (frozen coefficient sums identical -
  expected; sum-level only, element-wise hash unavailable in
  this exe).
- commit_rejects identical counts every refresh.
- unconverged cells differ by 85-6400 cells refresh-to-refresh.

Conclusion: one uniform 0.1% seed perturbation produces
percent-level trajectory differences that persist and
occasionally O(1) store-delta differences - the fixed point,
if it exists, is not reached identically from these two seeds.
This tests ONLY this perturbation direction; it does not prove
or disprove general seed independence or uniqueness.

### Phase 4 instrumentation (landed, future runs only)

The running arms' exe (8afd8591) cannot emit per-ray Q detail, so
a diagnostic-only extension was added to
`plane/fg_j_cmf_nonmono_v1.f` behind the existing `KN_Q_RESID`
control (default F, zero cost when off):

- `KN_QA_Q/VOLD/VN/VM` (NM,2) and `KN_QA_J/L` (NM,2) module arrays
  store, per ray and channel (1=inward/M, 2=outward/P), the
  one-pass residual, the sweep argmax cell (depth,freq) from the
  WRC tracker, the pre-pass production field value, the post-pass
  scratch value, and the column peak VMAXF(L) reconstructed from
  the untouched production field (identical to the sweep's own
  VMAXF because FLD_Q was an exact copy on entry).
- `FG_KNOP_REFRESH_V1` prints `q argmax ray/ch/same-as-prod-
  worst/prod resid/passes`, `q argmax j,l/vold/vnew/vmaxf/
  denom/floor-dom`, and a sorted per-ray-max distribution
  (max,p50,p95,p99 + counts above 1e-2/1e-1/1). The floor-dom
  flag answers whether the Q max sits on a cell below the
  1e-4*VMAXF significance floor.
- No sweep, march, commit, tolerance or relaxation logic was
  touched; capture happens after each channel's verify pass and
  before the shared-state restore. Compiles clean under
  production flags (only pre-existing NB_Q/NT_Q kind warnings on
  the original verify call sites).

PENDING: pop-release iteration in arm C (the M/C differentiator),
longer refresh trajectories, commit/maturity evolution, and a
future instrumented run to populate the Q argmax fields.

Pre-registered interpretation (from the work order):

- R converges, M converges, C stalls -> population/matter feedback
  is the leading cause (hypothesis D/C).
- R converges, M stalls, C stalls -> radiation-dependent
  scattering-source feedback (hypothesis B).
- R stalls, M stalls, C stalls -> the production refresh operator
  or geometry itself is non-contractive (hypothesis A).

## 6. Production population conditioning (Phase 3) - NOT YET RUN

Planned: warm restart of `m15_2p82d_VAL_knop_formal` near the
maxraw=1E7 condition with `POP_AUDIT=T`,
`POP_RUNAWAY_DEM=1E3`, `POP_RUNAWAY_TOPN=10`; SVD-based rank
diagnostics on selected floor-pinned depths via
`CMF_BLKBAND_V3`/`GENERATE_FULL_MATRIX_V3`. Deliverable:
`docs/audits/KNOP_P5_T1_PRODUCTION_CONDITIONING.md`.

## 7. Production geometry dependency graph (Phase 4) - PARTIAL

Prior audits (`KNOP_PRODUCTION_DEPENDENCY_GRAPH.md`,
`KNOP_P1_T1B_LAGGED_MAP.md`) established the manufactured-model
nilpotency result and the audit-side graph tooling
(`unit_tests/python/knop_graph.py`, `knop_topo.py`). The
production-geometry Tarjan SCC census and finite-difference
validation on the it125 checkpoint remain to be executed.
Deliverable: `docs/audits/KNOP_PRODUCTION_GEOMETRY_DEPENDENCY_RESULTS.md`.

## 8. Unresolved limitations

- `KN_MAT_FREEZE` runs `SOLVEBA_V13` and discards its update rather
  than skipping the call (the work order's "do not execute" is met
  in effect - POPS is bitwise unchanged - but SOLVEBA side
  effects on diagnostic counters still occur). If a literal
  no-call variant is required it is a small follow-up.
- The q-residual is a pass-1, omega = 1.0, Gauss-Seidel-ordered
  sweep-map defect on the per-ray-converged field; it excludes the
  cross-ray store commit and the frequency-ladder shift, so it is
  a per-ray map residual, not the full refresh Jacobian norm, and
  not comparable magnitude-for-magnitude with the production
  terminal residual (different pass number and omega state). Its
  raw-defect property silently depends on `KN_RAY_OMEGA` staying
  1.0.
- `KN_FREEZE_COEF` freezes the refresh coefficient snapshots only;
  the main per-frequency formal solve still sees live
  `CHI_EXT`/`ETA_EXT`/`CHI_COEF`/`ETA_COEF` rebuilt each call. Arm
  R therefore isolates nothing about the whole radiation
  iteration - a stall in R is consistent with defects in the live
  coefficient path as well as in the refresh operator or
  geometry.
- The sum-based coefficient fingerprints can conceal local
  coefficient changes through cancellation; only the post-launch
  FNV-1a hash (absent from the running arms' executable) detects
  element-level drift.
- Seed-2 perturbation is a uniform 0.1% scale of both I_NEXT
  stores; it probes seed dependence of the fixed point but is not
  a random-direction perturbation.
- Single-threaded arms are ~4x slower per iteration than the
  production 4-thread run; 40 iterations may not reach the
  it124-equivalent maturity horizon. Watchdog/iteration limits
  were not raised; arms stop on the same production criteria.

## 9. Evidence classification

- CONFIRMED: diagnostic flags reach VADAT and the executable
  (OUTGEN echo + symbol/string check); warm restart resumes at
  it125; unit tests pass (43/43 ray, lagged-field, status-text).
- CONFIRMED: freeze semantics as implemented match the work-order
  definitions (source-verified, section 3).
- CONFIRMED: ETA's J-dependent scattering term is separable
  upstream (`CHI_SCAT_CLUMP*RJ[_ES]`), so arms R and M differ
  exactly in whether that term refreshes (source-verified).
- CONFIRMED (source audit): `KN_FREEZE_COEF` scope is limited to
  the refresh-side `KN_*` snapshots/coefficients; the main formal
  solve's live `CHI_EXT`/`ETA_EXT`/`CHI_COEF`/`ETA_COEF` chain is
  untouched (section 3).
- CONFIRMED (source audit): `RES_Q` is the raw unrelaxed
  one-pass sweep-map defect - `DCELL` uses pre-relaxation `VNEW`,
  the `L_QP` exit precedes all `OMEGA_RAY` adaptation, and the
  verification pass runs at `OMEGA_RAY = KN_RAY_OMEGA = 1.0`
  throughout (section 3). It is a relative-change metric, not a
  norm, and is not magnitude-comparable with the production
  terminal residual.
- CONFIRMED (source audit): sum-based coefficient fingerprints
  are insensitive to compensating local changes; the FNV-1a hash
  diagnostic now covers element-level drift for future builds
  (running arms predates it).
- PENDING: whether the production refresh is contractive at fixed
  coefficients (arm R/R2), whether J-dependent scattering feedback
  stalls convergence (arm M vs R), and whether population-opacity
  feedback drives the production stall (arm C vs M/R).
- NOT TESTED: production population-matrix conditioning (Phase 3),
  production-geometry SCC census (Phase 4), ordered/unordered
  production A/B (explicitly deferred by the work order).

## 10. Recommended next implementation step

Await the four arms' first refresh summaries (march census, q
residuals, coefficient fingerprints, GS_MAXCHG trajectory). If arm
R converges while C stalls, the population-opacity Jacobian is the
defect site and Phase 3 conditioning diagnostics become the
priority. If all arms stall identically, the refresh operator or
production geometry is implicated and Phase 4's dependency census
becomes the priority. No solver changes until the evidence report
is reviewed.
