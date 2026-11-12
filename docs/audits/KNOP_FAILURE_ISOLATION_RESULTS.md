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

### Refreshes 10-13 (its 134-137) - DECISIVE ARM SPLIT

The arms have separated decisively:

- **M and C (live refresh-side coefficients, pops held):
  CONVERGING.** rf_resid fell 0.0311 -> 0.0163 (it133) ->
  0.0056 (it134-137), now ~2.8x the 2.0e-3 tolerance and flat.
  Unconverged cells collapsed 5258 -> 378. Q ~3.07-3.09. Still
  `held=T pops_applied=0.0` - this convergence is achieved with
  FROZEN populations, so population release is not required for
  the radiation field to converge.
- **R and R2 (frozen refresh-side coefficient snapshots):
  STUCK/DESTABILIZING.** rf_resid pinned ~0.030-0.031
  (~15x tol) through it137; unconverged ~8-11k cells; and the
  independent Q residual JUMPED from ~3.09 to ~10.0-10.6 at
  seq11-13 - the one-pass map defect tripled while M/C's stayed
  ~3.07. R2's rf_resid (0.0297) remains slightly below R's
  (0.0310).
- M/C remain bit-identical through it136 (C's it137 ITER_SUM
  also held=T); the population gate is still closed but is now
  moot for the radiation-convergence question - M already
  demonstrates convergence without pop updates.

INTERPRETATION (restricted table): this matches classification
**D** - "R stalls while M converges: investigate whether the
frozen first-pass coefficient snapshot is itself a poor or
destabilising approximation." Under the verified semantics,
KN_FREEZE_COEF freezes ONLY the refresh-side KN_CHI_EXT/
KN_ETA_EXT snapshots (which carry the live J-dependent
CHI_SCAT*RJ term); the main per-frequency formal solve keeps
live coefficients in BOTH arms. The split therefore isolates
the REFRESH-SIDE coefficient evolution: letting the snapshots
track the evolving J-dependent scattering source is necessary
for the lagged refresh map to contract; freezing them at their
first-pass values makes the refresh operator non-contractive
(or far weaker) even though the primary solve is unchanged.

Caveats: (a) pops are still held in ALL arms, so population
feedback remains untested - C has not released; (b) M's
convergence is at the REFRESH-residual level (KN_RF_RESID/
fg_rf = 0.0056): the in-sweep store drift (sw_maxchg ~0.23,
GS_MAXCHG argmax ~1.6-3.1e-2) remains well above the 2.0e-3
field-maturity tolerance, so the production convergence gate
is NOT yet satisfied and physical validation is not
established; (c) one warm-restart window from it125 - not a
cold-start claim.

### Refreshes 7-8 (its 131-132, continuing)

- C remains `frozen=T held=T`, `pops_applied=0.0` through it131 -
  the population-release gate has NOT opened inside the arm
  window so far. M/C remain bit-identical; the M-vs-C comparison
  stays inconclusive until C's first `held=F`/`pops_applied>0`
  iteration.
- Argmax update at it131: M/C's wandering argmax has now
  LOCALIZED onto (d=4, r=81, f=1249)=3.110e-2 for two
  consecutive iterations with `cycle` class at it131 - the
  shallow-cell wander may be settling into its own persistent
  cell, distinct from R/R2's deep cell (54,30,8712).
- R/R2 retained_max jumped to ~63-64 at seq7 (was 0.4-5.9) while
  rf_resid stays ~0.030/0.031 and unconverged ~12-15k cells -
  a large committed store delta coexists with the pinned
  residual; R2's rf_resid (0.0297) now sits slightly below R's
  (0.0310), an ordering not seen earlier.
- M/C unconverged cells down to 3946 with retained_max 0.096 -
  the best cell-level convergence of any arm, still with
  rf_resid pinned at 0.0311.

### Store-slot mapping verification (Phase 7 prep)

Traced in source (`fg_j_cmf_nonmono_v1.f`):

- `FREQ_IDX_ITER` is a ZERO-BASED sweep counter: during field
  frequency k (1-based), FREQ_IDX_ITER = k-1 (lines 242-243,
  3305-3319). The printed argmax f=8712 is FREQ_IDX_ITER, so
  the field frequency is k=8713.
- `WRITE_STORE`/`WRITE_STORE_CAP` write `STORE(KC,LS,LS2-1)`
  for LS2=2..NCF - field frequency LS2 lands in store slot
  LS2-1, and frequency 1 has no slot. Therefore argmax
  FREQ_IDX_ITER=8712 -> field k=8713 -> store slot 8712. The
  printed f value coincides with the store slot index, but
  only because slot = k-1 = FREQ_IDX_ITER; the mapping is
  verified, not assumed.
- Argmax `d`=54 is the geometry index K (1..KN_ND), NOT the
  model depth: model depth = COMP_DEPTH_IDX(54), ray point =
  J_PNT(54,LS). Argmax `r`=30 is the ray index LS.
- Commits are PER-RAY inside KNOP_RAY_SOLVE_V1 (~lines
  8439-8459): WRITE_STORE for converged rays,
  WRITE_STORE_CAP for non-converged under FG_COMMIT_CONV,
  per channel (I_M_DEPTH_NEXT / I_P_DEPTH_NEXT stores).

Minimal-diagnostics feasibility (for the recommended
follow-up, NOT implemented): a gated watch-cell inside
WRITE_STORE/_CAP can capture, at (KS==K_t, LS2-1==slot_t) for
the target ray, V1 (proposed post-cap), V2 (pre-commit
retained), DLOC (relative change + denominator), the post-write
STORE value, and the call-site verdict - all read-only
captures behind a VADAT-gated logical, ~15 lines, no numerical
state touched. Validation = bit-compare OUTGEN/store dumps
instrumented vs uninstrumented on a short run. The same capture
in both WRITE_STORE and WRITE_STORE_CAP covers accepted and
capped paths; rejected rays never reach either routine, so the
verdict print at the call site covers the rejected case.

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

## 11. M-arm continuation and watch-cell telemetry - IN PROGRESS

Authorized follow-up: determine whether the live-coefficient M
pathway reaches full radiation-field maturity with populations
held, and trace the commit behaviour of the maturity-controlling
cell. No production solver, convergence, relaxation or population
changes are part of this work.

### 11.1 Watch-cell diagnostic (implemented, committed)

Source commits: `1c875cd5` (capture + call-site records),
`86fca689` (control-variable declarations), `247c9ad6` (VADAT
reads), `5a60e316` (write-format fix). All changes are
diagnostic-only and default off.

VADAT keys: `KN_WATCH_CELL` (logical, default F),
`KN_WATCH_K`, `KN_WATCH_LS`, `KN_WATCH_SLOT` (integers). When
`KN_WATCH_K>0` and the logical is T, `SET_FG_WATCH_CELL` arms the
capture; K<=0 disables.

Capture points:

- `WRITE_STORE` and `WRITE_STORE_CAP` record, at the matching
  (geometry index KS, ray LS, store slot LS2-1): the raw proposed
  field `FLD(JP2,LS2)`, the committed value `V1` after any cap,
  the retained store value `V2` read before the write,
  `VMAX_F(LS2)`, the DLOC denominator `MAX(V2,1.0E-4*VMAX_F)`,
  the resulting `DLOC`, and the post-write `STORE` value.
- `KNOP_RAY_SOLVE_V1` prints one record per solve for the armed
  ray after the whole commit decision chain, so jump-rejected and
  non-converged-rejected rays (which never reach either write
  routine) still report their `FG_RAY_VC` verdict plus exit class
  (`NXIT`), pass count (`NPASS`) and residual.
- Channel attribution uses `KN_WC_CH`, set to 1 before the M
  channel write call and 2 before the P channel call, both gated
  on `LS==KN_WATCH_LS`.

Index verification (source-traced, not assumed): the GS_MAXCHG
argmax prints `(d,r,f)=(54,30,8712)`. `d=54` is the geometry index
K into the compressing-depth arrays; the physical depth is
`COMP_DEPTH_IDX(54)` and the ray point is `J_PNT(54,30)`. The
printed `f=8712` is `FREQ_IDX_ITER` (0-based during the sweep);
the field frequency index is `FREQ_IDX_ITER+1=8713` and the store
slot written is `LS2-1 = FREQ_IDX_ITER = 8712`. The watch target
is therefore `KN_WATCH_K=54, KN_WATCH_LS=30,
KN_WATCH_SLOT=8712`.

OpenMP safety: the parallel ray loop gives each ray to exactly
one thread, and every capture write is gated on
`LS==KN_WATCH_LS`, so only the armed ray's thread touches the
`KN_WC_*` state - no cross-thread capture race. The print block
runs inside `!$OMP CRITICAL(KN_WATCH_CRIT)` so concurrent ray
solves cannot interleave OUTGEN lines. Read-only: no floating
point operation, store update or convergence decision is altered;
with `KN_WATCH_CELL=F` the added work is four predictable branches
per store write.

A format-string defect (A-descriptor/item count mismatch) was
found and fixed before any instrumented run (`5a60e316`); the
first diagnostic build `3854b948` was discarded unused.

### 11.1a Diagnostic-completeness pass (Phase 2, commit `6ce9a623`)

Review of the first telemetry revision found three gaps; all are
now closed diagnostic-only:

- Exit class: `NXIT` is the residual-setting exit, not the
  per-channel solver class. The record now prints `IXIT_M`,
  `IXIT_P` and `NXIT` separately (classes: 0 exact, 1 tol,
  2 bulk, 3 stall, 4 cap, 5 diverged, 6 plateau), plus the ray
  commit delta `FG_RAY_DCH(LS)` next to the verdict.
- Two caps separated: `KN_WC_V1I` rebuilds the value after
  `NONMONO_INEXT_CAP`/sanitize and before the
  `FG_COMMIT_NC` provisional limiter, reconstructed arithmetically
  from the raw `FLD` value inside the capture block (mirrors the
  cap logic; no extra pass over the cell loop). `KN_WC_FL` records
  bit0 = NONMONO_INEXT_CAP fired, bit1 = sanitize-to-floor path,
  bit2 = FG_COMMIT_NC fired (WRITE_STORE_CAP only; WRITE_STORE
  has no NC step so bit2 stays clear there).
- DLOC branch: `KN_WC_BR` marks 0 = `V2>0` (the printed `den` is
  the real `MAX(V2,1E-4*VMAX_F)` DLOC denominator) vs 1 =
  `V2<=0` (DLOC is the 0/1 floor test; `den` then holds the
  `1E-4*VMAX_F` significance floor, not a DLOC denominator).

### 11.2 Executable provenance

- Pinned production executable (running arms R/R2/M/C):
  `exe/cmfgen_dev.exe` md5 `8afd859176ca0177c9de31fc6a400ffc`,
  restored bit-identical after every diagnostic build and
  re-verified after the latest build.
- Diagnostic executable: md5 `b4686fa157edb8c94b0dc098636d91a3`
  (supersedes `9d64ca53`, which lacked the Phase-2 completeness
  fields; the earlier `3854b948` build was discarded for the
  format defect), staged at `/tmp/cmfgen_diag_install/exe/
  cmfgen_dev.exe` with `com/`, `lib/`, `mod/`, `obs_template/`
  symlinked from the production tree so `cmfdist` resolves. The
  batch script resolves `CMFGEN_PROG` at exec time, so the
  already-queued neutrality jobs pick up this exe automatically.
  Source state: main through `6ce9a623` plus concurrent uncommitted
  tree state at build time (see 11.4).
- Build note: a filesystem bus error during the build corrupted
  53 small members of `lib/libsubs.a`; all were recompiled from
  source and re-archived, the archive re-scanned clean (zero
  non-ELF members), and `make cmfgen` then completed with exit 0.
  Unit tests `test_knop_ray`, `test_knop_lagged_field`,
  `test_knop_status_text`: 35/35 pass.

### 11.3 Numerical-neutrality validation - RUNNING

Two identical warm-restart runs from the FRESH_knop it125
checkpoint, single thread, `KN_MAT_FREEZE=T`, `KN_Q_RESID=T`,
`init_num_its=2`, on the interactive partition:

- watch_A (control): `KN_WATCH_CELL` absent (default off), job
  `13071411`.
- watch_B (instrumented): `KN_WATCH_CELL=T`, `KN_WATCH_K=54`,
  `KN_WATCH_LS=30`, `KN_WATCH_SLOT=8712`, job `13071413`.

Both use the diagnostic executable `9d64ca53`. Pass criterion:
bitwise-identical numerical state (FGINEXT_GS store, GS_STATE,
OUTGEN numerical diagnostics, commit decisions) between A and B.
The watch_B `KN_WATCH` record itself is expected output, not a
numerical difference. If any numerical difference appears, the
continuation is blocked pending investigation.

### 11.4 Concurrent-tree caveat

A second working process is active in this repository: it
modified `subs/sobew*.f` and `unit_tests/` files during the
diagnostic build window and landed audit commits (`16adf737`
through `e6ed21a7`) plus `5a60e316` (the watch-format fix,
which swept up the pending source edit). The diagnostic
executable therefore contains whatever `sobew*` state existed at
link time; both neutrality arms share that state, so the A/B
comparison remains valid for the watch diagnostic itself. The
sobew changes are unrelated to the KNOP paths under test.

### 11.5 Arm status at time of writing

M reached it146 (rf_resid ~0.0056, unconverged cells ~200-250,
pops still `held=T`, `pops_applied=0.0`). R/R2 remain pinned near
rf_resid ~0.030 with the deep argmax (54,30,8712) in `cycle`
class. Two new measurements change the continuation picture:

- M's GS_MAXCHG argmax has localized on a shallow cell
  `(d,r)=(3,82)` since ~it139 at value ~5.6e-3, now `cycle`
  class (alternating 5.58/5.74e-3 at it145-146). The frequency
  column overflows the argmax-history I-field to `****`, so the
  printed record cannot confirm FREQ_IDX_ITER; the user-supplied
  value 1636 is used as the continuation watch slot (store slot =
  FREQ_IDX_ITER per the verified mapping). The argmax is currently
  stable on this cell; a single-cell watch suffices while it
  persists - retarget if it migrates again.
- M's independent one-pass Q residual is now growing on the same
  trajectory shape seen in R/R2: 3.09 -> 5.26 -> 10.96 -> 15.35
  -> ~20.7 and holding (it139-146), while rf_resid stays ~0.0056
  and unconverged cells stay ~200. The running executable
  predates the KN_QA_* argmax instrumentation, so the Q argmax
  ray/cell/denominator is not available for the live arms; the
  diagnostic exe used by the continuation emits it. Whether this
  Q growth is the same floor-normalized weak-column defect or a
  genuine one-pass-map divergence is exactly what the completed
  telemetry plus Q argmax output must resolve.

The M continuation will use M's final checkpoint once the arm
terminates, with `KN_FREEZE_COEF=F`, `KN_MAT_FREEZE=T`,
`KN_Q_RESID=T`, the watch cell armed at (K=3, LS=82,
SLOT=1636), single thread, and no changes to `KN_RAY_TOL`,
`KN_RAY_MAXPASS`, `KN_RAY_BULK_TOL`, `INEXT_RELAX`,
`FG_COMMIT_NC`, `FG_COMMIT_JUMP` or population solver tolerances.
It launches only after the A/B neutrality comparison passes and
is reviewed.

## 12. Phase-2 review: index mapping and the growing Q residual

### 12.1 M watch-cell index mapping - VERIFIED from source

The GS_MAXCHG argmax record `(d,r,f)` is written at
`fg_j_cmf_nonmono_v1.f:5086-5087,5109-5110,5146-5147,5163` as
`FG_GS_ARGK=K`, `FG_GS_ARGR=LS`, `FG_GS_ARGL=FREQ_IDX_ITER`,
inside the in-sweep I_NEXT write-back. Verified chain:

- `K` iterates `DO K=1,ND` with `KC=COMP_DEPTH_IDX(K)` (line
  ~5037); `KN_ND=ND` (line 2904), so `d` is the model-depth index
  in exactly the same convention `WRITE_STORE` uses for `KS`
  (`DO KS=1,KN_ND`, `KC=COMP_DEPTH_IDX(KS)`, line ~11040). The
  watch print reports `COMP_DEPTH_IDX(K)` (sparse store index)
  and `J_PNT(K,LS)` so the physical location is unambiguous.
  Correction to section 11: `d` is the model-depth index K, not
  a separate "geometry index"; the sparse store row is
  `COMP_DEPTH_IDX(d)`.
- `f` is `FREQ_IDX_ITER`, the 0-based sweep index. The lagged
  store `I_M_DEPTH_NEXT`/`I_P_DEPTH_NEXT` is allocated
  `(N_COMP,NP,NCF_IN)` (line 2871) and is written at third index
  `FREQ_IDX_ITER` directly (lines ~5168,5174).
- `WRITE_STORE`/`WRITE_STORE_CAP` are called with those arrays
  as `STORE` and write `STORE(KC,LS,LS2-1)` for `LS2=2..NCF`, so
  store slot `s = LS2-1`. Matching `s = FREQ_IDX_ITER` gives
  `LS2 = FREQ_IDX_ITER+1`, i.e. slot `s` holds field frequency
  `s+1`. For the M argmax `FREQ_IDX_ITER=1636` the watch slot is
  `KN_WATCH_SLOT=1636` (slot < NCF~56k, in range).
- Caveat: the argmax-history f column overflows its I-field to
  `****` for M/C, so 1636 rests on the user's trace, not the
  printed record; the continuation's own KN_QA/GS argmax output
  will confirm the live slot before interpretation.
- Frequency ladder: `DO_NONMONO_FREQ_MATCH` is not set in the M
  VADAT (auto-set applies only to the SOB LAM_ES=F path), so the
  grid is static; slot indices are stable across refreshes
  regardless.

### 12.2 What Q measures vs what the refresh residual measures

- `KN_Q_RESID` (Q): for each ray and channel, one scratch
  `SWEEP_FLD`/`SWEEP_FLD_ORD` pass is applied to a copy of the
  committed field at `OMEGA_RAY=1.0` (no relaxation, pre-exit
  checks bypassed for the measurement). `DCELL =
  |VNEW-VOLD|/MAX(VOLD,1E-4*VMAXF)` per cell; Q is the per-ray
  maximum over cells and frequencies. It measures the defect of
  ONE sweep-map application on the CURRENT committed field -
  the raw fixed-point residual of the inner iteration, with no
  commit, no cross-ray store write-back, and no outer
  relaxation.
- Refresh residual (`KN_RF_RESID`/`GS_MAXCHG` family): the
  store-to-store change across a complete refresh - the
  committed `I_NEXT` delta after per-ray solves, commit
  verdicts, caps and rejects. It measures the outer
  fixed-point map including cross-ray commits.
- Neither diagnostic includes population updates (held in all
  arms) nor frequency-ladder shifts (static grid). They are
  different norms on different maps: Q is per-ray max relative
  defect of one sweep pass; the refresh residual is the
  committed store drift. Absolute values are not comparable.

### 12.3 The Q growth is NOT M-specific - measured

| it | R rf_resid | R Q | M rf_resid | M Q |
|----|-----------|------|-----------|------|
| 136| 0.0308 | 10.6 | 0.00559 | 3.09 |
| 140| 0.0308 | 28.4 | 0.00558 | 10.96 |
| 142| 0.0308 | 324.4| 0.00558 | 20.66 |
| 146| 0.0308 | 325.0| 0.00558 | 20.67 |
| 147| - | - | 0.00589 | 58.9 |

Both arms now show the same phenomenon at different rates: the
one-pass map defect grows by orders of magnitude while the
committed refresh residual stays flat and unconverged cells
stay ~200 (M) or ~8k (R). `q_rays=84` is constant - every ray
reports a Q value, so the max is not a census artifact.

Working interpretation (hypothesis, not proven): the commit
machinery - `FG_COMMIT_NC` provisional caps, jump rejects,
non-converged rejects - is suppressing exactly the store updates
that would repair the cells the one-pass map flags. The
committed store therefore sits near a fixed point of the
GATED map (flat rf_resid, stable unconverged count) while the
raw sweep-map defect on that store grows (Q). If true, Q growth
is a symptom of commit gating, not of a new instability, and the
watch-cell record at the maturity cell (verdict + capped/rejected
+ V1I/V1C/V2/DLOC) is the decisive measurement. Alternatives not
yet excluded: floor-normalized weak-column dominance (needs the
KN_QA denominator), a genuine one-pass-map instability the
commit path was always masking, or transient amplification in a
non-normal operator (Trefethen & Embree 2005 - possible, not
established).

Discrimination requires the KN_QA_* output (argmax ray/channel/
depth/freq, VOLD/VNEW/VMAXF/denominator, floor-dom flag, p50/
p95/p99 distribution, threshold counts) which only the
diagnostic executable emits - the running arms' exe predates it.
The continuation will supply it.

### 12.4 M maturity status - explicit

M has NOT reached production maturity: rf_resid ~0.0056-0.0059
vs tolerance 2e-3 (~3x), GS_MAXCHG argmax ~5.6e-3 in `cycle`
class, sw_maxchg nonzero, Q now ~59 and rising, populations
held. No convergence gate is satisfied; "converging" describes
only the early refresh-residual trend, which has itself stalled.

## 13. Phase-2 work order: Q-diagnostic completeness, A/B neutrality evidence, residual definitions

Scope: diagnostic-only. No solver equation, convergence threshold,
commit rule, sweep ordering, population or opacity path was
modified. Source revision at time of writing: `88dcf220` plus the
diagnostic diff in `plane/fg_j_cmf_nonmono_v1.f` described in 13.2
(committed together with this section).

### 13.1 A/B neutrality experiment - PENDING (jobs never executed)

Scheduler state (verified 2026-11-12 ~14:30):

- watch_A job `13071411` PENDING (reason Priority); chained
  dependent `13071412` PENDING (Dependency).
- watch_B job `13071413` PENDING (reason ReqNodeNotAvail, may be
  reserved); chained dependent `13071414` PENDING (Dependency).
- Neither job has ever executed: no batch output, no OUTGEN, no
  iteration count. All runtime comparisons (1.4) are BLOCKED on
  scheduling, not on missing evidence.

Pre-execution evidence (all verified directly):

- Executable resolution: `batch.sh` sets
  `CMFGEN_PROG=$cmfdist/exe/cmfgen_dev.exe` and the submitted
  `run_cmfgen_SINGLE_FINAL.slurm` exports
  `cmfdist=/tmp/cmfgen_diag_install`. Both jobs (and their chained
  dependents, same script) will therefore exec
  `/tmp/cmfgen_diag_install/exe/cmfgen_dev.exe`, currently md5
  `b4686fa157edb8c94b0dc098636d91a3` (verified by md5sum on the
  resolved path, not inferred from the submit script). The path
  was NOT overwritten while preparing the new diagnostic build
  (13.2), so the A/B experiment keeps the intended binary.
- Slurm script diff A vs B: only `--job-name`, `--chdir`, the
  SN_HYDRO_DATA path check, the requeue job name and the
  `config_resolved.yaml` path differ - all path substitutions.
  `batch.sh` is byte-identical (same SHA-256 in both dirs).
- `config_resolved.yaml` diff: `model_name`/`model_dir` path
  strings plus exactly the four permitted keys on watch_B:
  `KN_WATCH_CELL: T`, `KN_WATCH_K: 54`, `KN_WATCH_LS: 30`,
  `KN_WATCH_SLOT: 8712`.
- Input checksums (SHA-256 over every regular file in each
  `cmfgen/2.82days/` run dir plus `config_resolved.yaml`): all
  numerical inputs are byte-identical between A and B. Key
  values: `FGINEXT_GS` =
  `7f080ec99da83f317d92a7a10e06c449f422a54f51ef1f1f89eddba46a3e23ce`,
  also identical to
  `m15_2p82d_FRESH_knop/cmfgen/2.82days/FGINEXT_GS` - same warm
  restart generation (the it125 checkpoint copy). `MODEL`,
  `VADAT_INIT`, `IN_ITS`, `SN_HYDRO_DATA`, `T_IN`, `RDINR`,
  `NUC_DECAY_DATA`, `CUR_MODEL_DATA`, `ca2_col_data_filtered`,
  `arnaud_rothenflug.dat`, `POINT1`, `POINT2`, `SCRTEMP`, `obs/`
  contents, `batch.sh`, `submit_cmfgen.sh` all match. Only
  `audit_report.*`, `MODEL_SPEC.audit.md`, the three generated
  `.slurm` files and `config_resolved.yaml` differ (job metadata,
  paths, watch keys).
- Compiler/runtime: identical batch environment (module
  intel/2024.0, mkl/2026.0, python-waterboa/2024.06, pinned
  MKL `LD_LIBRARY_PATH`), `--cpus-per-task=1`, and the same
  `CPUS_REQ`/`squeue %C` resolution logic -> single-thread
  execution expected for both arms.
- Executable provenance: `git diff 6ce9a623..88dcf220 --
  plane/fg_j_cmf_nonmono_v1.f` is EMPTY - the KNOP solver,
  KN_Q/KN_QA/KN_QC/KN_WATCH source in `b4686fa1` corresponds
  exactly to committed `6ce9a623` state for that file. The only
  `.f` changes in the range are comment-only (`sobew*.f` in
  `2e5922ff`, `sobjbar_sim.f` in `05c20d44`). Provenance of the
  complete executable (every archive member at link time,
  including the documented bus-error member rebuild) is
  **NOT PROVEN** - an MD5 cannot establish it.

Neutrality verdict per criterion:

| Criterion | Status |
|---|---|
| Identical executable MD5 at exec time | PENDING - same resolved path verified; runtime exec not yet observed |
| Verified executable provenance | PARTIAL - KNOP file PROVEN equivalent to `6ce9a623`; complete exe NOT PROVEN |
| Identical restart/input checksums | PASS (pre-execution) |
| Identical compiler/runtime/thread | PASS (by identical script; runtime confirmation pending) |
| Byte-identical persistent state | PENDING - no execution |
| Identical convergence metrics | PENDING - no execution |
| Identical commit decisions/counters | PENDING - no execution |
| No unexplained differences | PENDING - no execution |

Overall: **NOT PROVEN**. The experiment is correctly configured
and all pre-execution evidence passes, but the gate cannot be
called PASS until both jobs run and their outputs compare. If
the jobs remain starved, the gate stays open; the M-arm
continuation remains blocked regardless (13.6).

Update 2026-11-12 ~15:30: the first scheduled attempt FAILED in
7 s on both arms - batch.log shows `nice: '/tmp/cmfgen_diag_
install/exe/cmfgen_dev.exe': No such file or directory`. Cause:
/tmp is node-local; the job landed on viper07 where the install
did not exist. The identical b4686fa1 exe plus com/lib/mod/
obs_template symlinks was re-staged at the shared path
/viper/ptmp1/sbrennan/SN2024pba/cmfgen_diag_install (md5
re-verified), the cmfdist export in all four slurm scripts
(SINGLE_FINAL+OBS for A and B) was repointed, and the chain was
resubmitted: watch_A FINAL 13073518 / OBS 13073519, watch_B
FINAL 13073520 / OBS 13073521 (afterok). Original jobs
13071411/13 FAILED, dependents 13071412/14 cancelled by the
afterok dependency. Experiment design unchanged: same exe md5,
same inputs, same watch keys; only the install path moved to a
shared filesystem.

### 13.2 Q-residual reporting corrected - diagnostic-only rebuild

ChatGPT's review identified that the existing Q output could not
distinguish a localized defect from population-wide growth and
that the percentile block sorted the `-1` unreported-ray
sentinel into the distribution. `plane/fg_j_cmf_nonmono_v1.f`
now carries a diagnostic-only diff (committed together with
this section):

- Ray census `KN_QS_STAT(LS)`: 0 unreached, 1 evaluated,
  2 diverged-skip, 3 degenerate/`NI<2` early return,
  4 verify-alloc failure, 5 field not allocated,
  6 ray-local alloc failure. Printed as
  `q ray census eval,div,ni<2,allocF,nofld,rayalloc,unreached=`.
- Per-ray percentile block now excludes `KN_Q_LS < 0` before
  sorting (previously the `-1` sentinel biased p50/p95/p99
  downward whenever any ray was unreported) and is labelled
  `q dist per-ray-max` - explicitly a distribution of per-ray
  MAXIMUM residuals, not a cell distribution.
- Cell-level statistics `KN_QC_*`: every `DCELL` computed in
  `MARCH` during the two verify passes is accumulated per
  channel (1=M inward, 2=P outward) into `KN_QC_HIST` (0.1-dex
  log bins, 1e-30..1e10 plus under/overflow edges), exact
  threshold counts `KN_QC_NT` (>1e-2, >1e-1, >1), exact
  `KN_QC_MAX`, and sample count `KN_QC_N` (INTEGER*8; the full
  `NP*NI*NCF*2` domain is never materialized). Accumulation is
  gated by `KN_Q_ACC`/`KN_Q_ACC_CH`, armed only around the two
  verify `SWEEP_FLD`/`SWEEP_FLD_ORD` calls inside the
  `KN_Q_RESID` block; production passes evaluate one predictable
  branch per DCELL site and are otherwise untouched. Printed as
  `q cells ch= samples= max= p50,p95,p99= n>1e-2,1e-1,1=`; the
  quantiles are bin-centre estimates (+/-0.05 dex), thresholds
  and max are exact. NaN DCELLs count in `KN_QC_N` but are
  excluded from the histogram.
- Argmax telemetry extended: the Q argmax now prints the exact
  local normalized defect `|VNEW-VOLD|/DEN` with
  `DEN=MAX(VOLD,1E-4*VMAXF)`, the explicit floor test
  `VOLD < 1E-4*VMAXF` (recorded as `floor-dom=` with both
  operands printed - floor dominance is no longer inferred from
  VOLD being small relative to the column peak), and a
  production-coincidence classification comparing ray, ray-local
  depth and field frequency against `FG_RAY_WI/WL` (the
  production worst-cell record, same index space):
  `NO-VALID-PROD-ARGMAX` / `SAME-RAY-SAME-CELL` /
  `SAME-RAY-DIFF-CELL` / `DIFF-RAY`. Same ray alone is not
  reported as coincidence; the production residual records no
  channel so channel coincidence is n/a.
- Cycle-class comparison: the store-drift argmax
  (`FG_GS_ARGK` model-depth, `FG_GS_ARGR` ray, `FG_GS_ARGL`
  zero-based `FREQ_IDX_ITER`) is mapped into the Q cell's index
  space. `J_PNT(model_depth,ray)` returns the ray point whose
  radius equals `R(model_depth)` to 1e-12 (verified at the
  J_PNT fill, bounds-checked fill), so the GS depth maps to ray
  point `J_PNT(FG_GS_ARGK,IQ_R)`; the frequency maps as
  `KN_QA_L-1` (field index is 1-based, sweep index 0-based).
  `gs-cell-match=` is exact triple coincidence; `gs-cls=`
  prints the argmax-trajectory class. An earlier revision of
  this comparison indexed `J_PNT` with the ray-local depth
  (wrong direction, potential out-of-bounds since inserted ray
  points can exceed `ND`); corrected before build.
- Removed dead per-channel argmax arrays (`KN_QC_RAY/J/L/Q/
  VOLD/VN/VM/DEN/FLOOR`) allocated but never populated; the
  per-channel argmax is recovered from the existing `KN_QA_*`
  table.

Diff audit: every changed executable-region line is either a
gated diagnostic mark, an allocation/deallocation, an
`IF(KN_Q_ACC)CALL KN_Q_ACCUM(DCELL)` placed AFTER the existing
`RES`/`RESL` max updates at the four MARCH DCELL sites, or a
summary WRITE. The `IF(NI .EQ. 1)RETURN` early exit was
restructured to a block adding only the gated status mark -
control flow identical. No sweep equation, residual formula,
commit path, threshold, ordering or coefficient handling was
touched.

Build and validation:

- `ifort` syntax/object compile clean (only pre-existing
  #6075 dummy-arg warnings on the verify-pass SWEEP calls).
- `make cmfgen` exit 0; new diagnostic executable md5
  `34f5c1d29877f4145756dda87d25135f`, staged at
  `/tmp/cmfgen_qcells_install/exe/cmfgen_dev.exe` (com/lib/mod/
  obs_template symlinked like the existing diag install) and
  preserved as `/tmp/cmfgen_qcells_diag.exe`.
- Pinned production `exe/cmfgen_dev.exe` restored bit-identical:
  md5 `8afd859176ca0177c9de31fc6a400ffc` re-verified after the
  build. `/tmp/cmfgen_diag_install` (`b4686fa1`) untouched, so
  the pending A/B jobs keep their specified binary.
- Unit tests on the final source: `test_knop_ray` 43/43 PASS,
  `test_knop_lagged_field` all PASS (130 checks),
  `test_knop_status_text` all PASS.
- Provenance of `34f5c1d2`: built from working tree = `88dcf220`
  + the diagnostic diff described in 13.2, which is committed
  together with this documentation - the exe's source state is
  reproducible from that commit alone.

### 13.3 Residual definitions - precise statements

Production refresh residual (`RES_LS` -> `KN_RF_RESID`):

- `RES_M`/`RES_P` are the values returned by the production
  `SWEEP_FLD`/`SWEEP_FLD_ORD` calls on the committed-direction
  fields `FLD`/`FLD_P` (inward/outward). Each is the maximum
  `DCELL` of the LAST completed pass of that inner solve, where
  `DCELL = |VNEW-VOLD|/MAX(VOLD,1E-4*VMAXF)` is computed inside
  `MARCH` BEFORE the `OMEGA_RAY` relaxation is applied to the
  stored update - i.e. the raw pre-relaxation proposed change,
  normalized by the same MARCH denominator as Q.
- `RESID = MAX(RES_M,RES_P)` per ray -> `RES_LS(LS)`; diverged
  rays (`NDIV=1`) are pinned to `RESID=1.0` and excluded from
  `KN_RF_RESID` by the `NDV_LS` mask (line ~7589). The refresh
  residual is `MAXVAL(RES_LS)` - a per-ray maximum over the last
  production pass of the inner ray solve.
- It is measured on the last completed production pass of the
  ray's own field iteration. It does NOT measure store-to-store
  drift across refreshes, does NOT include cross-ray commits,
  commit rejection/capping/healing, or frequency-ladder shifts
  (the grid is static), and by itself does not establish
  convergence of the complete refresh/store-to-store map.
- Effective relaxation actually used: `KN_RAY_OMEGA` (default
  1.0, not overridden in the arm VADATs) seeds `OMEGA_RAY` at
  each sweep entry; the production sweep then adapts it
  pass-to-pass (x0.7 on growth, x1.15 on contraction, floor
  0.45, cap 1.0) - so "OMEGA_RAY=1.0" describes the verify pass
  only, not necessarily every production pass. `INEXT_RELAX`
  (default 0.5, adapted per refresh from the pooled contraction
  samples at ~line 7645) relaxes the commit/store update, a
  different knob on a different map.

Independent Q residual (`KN_Q_RESID` -> `KN_Q_LS`):

- After the production solve returns, the verify block copies
  the committed field (`FLD_Q=FLD` for the M channel,
  `FLD_Q=FLD_P` for P) and calls the SAME
  `SWEEP_FLD`/`SWEEP_FLD_ORD` routine with `L_QP=.TRUE.`.
  `L_QP` exits after exactly one completed pass (`IXIT=9`)
  BEFORE the OMEGA-adaptation block and before the
  `KN_RES_PASS`/`RHO_WIN` telemetry writes, so the verify pass
  neither pollutes production relaxation statistics nor
  iterates to convergence.
- The scratch pass therefore uses: the same sweep equations and
  ordering as production; the same module coefficient arrays
  (`KN_CHI_COEF`/`KN_ETA_COEF` state as left by the production
  solve / `KN_FREEZE_COEF` snapshots); the same
  `MATCH_FREQ`/`SRC_PRIME`/`NFM_L`/`DLNM_L` matched-frequency
  machinery (recomputed per frequency inside the sweep);
  `OMEGA_RAY = KN_RAY_OMEGA` for its single pass (no adaptation
  occurs before the `L_QP` exit); and the scratch field is
  discarded - nothing is committed to `I_M/I_P_DEPTH_NEXT`.
- Q per ray = `MAX(Q_M,Q_P)` where each is the sweep's returned
  pass-1 residual = max `DCELL` over cells/frequencies of that
  one map application on the scratch copy.
- Answers to the five dependency questions, for BOTH residuals:
  1. Cross-ray store commits: NEITHER - both measure inside one
     ray's field solve before WRITE_STORE/WRITE_STORE_CAP.
  2. Frequency-ladder shifts: NEITHER - the frequency grid is
     static within a refresh; `DO_NONMONO_FREQ_MATCH` remaps
     donor READS (`L+NFM_L`), it does not move slots.
  3. Spatial interpolation of matched-frequency donors: BOTH -
     `SRC_PRIME` reads `FLD(JPR,L+NFM_L)` with `JPR1/JPR2`
     interpolation weights inside MARCH, identically in
     production and verify passes.
  4. Full store seed -> ray solve -> commit -> next refresh
     chain: NEITHER - that is the outer fixed-point map and
     requires comparing successive committed store generations
     (e.g. `FGINEXT_GS` fingerprints or GS_STATE drift).
  5. Commit rejection/capping/healing: NEITHER directly - the
     production ray residual is computed before the commit
     verdicts; Q never reaches commit. Indirectly, rejected
     cells keep their store values, so both residuals measure a
     field that already reflects past gating.
- Neither Q nor the production ray residual is a complete outer
  fixed-point residual; do not compare their absolute values as
  the same norm.

Correction to 12.1: `DO_NONMONO_FREQ_MATCH` is not merely
"unset" - under `DO_NONMONO_KNOP` it is FORCED T
(`rd_control_variables.f` ~3453, banner-forced even over
explicit F), and M's OUTGEN prints `freq_match= T`. The
conclusion drawn there still holds (the grid is static and slot
indices are stable across refreshes) but the mechanism is
active donor remapping, present identically in production and Q
passes.

### 13.4 M watch-cell mapping - verified, one correction

Reported M argmax: model-depth `d=3`, ray `r=82`,
`f=FREQ_IDX_ITER=1636` (I-field overflow `****` in the printed
history; 1636 from the user trace).

| Quantity | Value | Basis |
|---|---|---|
| Model depth K | 3 | argmax `d` is the model-depth index (12.1) |
| Sparse store row KC | 3 | `COMP_DEPTH_IDX(K)=K`: OUTGEN prints `compressing points=70 of ND=70`, so every depth is in the compressing set and the sequential numbering is the identity map |
| Ray LS | 82 | argmax `r` |
| Ray-local point | `J_PNT(3,82)` | runtime value; not printed by the running exe - the continuation's watch print reports it (the watch record prints `J_PNT(K,LS)` per 12.1) |
| Field frequency index L | 1637 | `FREQ_IDX_ITER` is 0-based; `L = FREQ_IDX_ITER+1` |
| Matched donor slot | `L + NFM_L(i)` per cell | `DO_NONMONO_FREQ_MATCH=T` forced under KNOP; `NFM_L>=1` so the upwind donor differs from the destination slot - this is donor remapping, NOT a different store slot |
| Store slot (M and P) | 1636 | `WRITE_STORE`/`WRITE_STORE_CAP` write `STORE(KC,LS,LS2-1)`; slot `s` holds field `s+1`; `KN_WATCH_SLOT=1636` is in range (`NCF=56377`) |
| Physical frequency | nu = 7.3691905E+16 Hz (73.69E15 in table units) | `CONT_FREQ` row `ML=1637` (= field index L, same 56377-point grid ordering); the argmax history itself prints `f=1636` on the `(3,82,1636)` entries, independently corroborating the user trace |
| Wavelength | Lam = 40.68 Ang | same `CONT_FREQ` row |

`KN_WATCH_SLOT=1636` remains the correct destination-store slot
for the continuation. The matched-donor distinction matters for
interpretation: the watch record captures the DESTINATION cell;
the donor it reads may sit at a different frequency slot.

### 13.5 M-arm status at time of writing

Still RUNNING (job `13065896`, ~4.5 h wall at last check),
iteration ~149, `pops_held=T`, `KN_MAT_FREEZE=T`. Q max reached
~79.0 and is still climbing (3.09 -> 5.26 -> 10.96 -> ~20.7 ->
58.9 -> 79.0 across its139-149); `rf_resid` ~0.0056-0.0059 flat;
`q_rays=84` constant; GS argmax remains `(3,82,****)` in
`cycle`/`margnl` classes. Final metrics (5.1) remain PENDING
until the arm terminates; M is NOT converged and meets no
production maturity criterion (12.4 stands).

### 13.6 Continuation gate status - NOT AUTHORIZED

Required gates and current state:

1. A/B neutrality PASS - FAILS TO CLOSE: jobs never executed;
   verdict NOT PROVEN (13.1).
2. All four arms terminated - FAILS: R/R2/M/C all RUNNING with
   dependent requeue jobs queued.
3. Final M metrics recorded - PENDING (arm live).
4. Q argmax identified with new instrumentation - PENDING: the
   new cell-level/coincidence/census output exists only in the
   `34f5c1d2` diagnostic exe, which no running arm uses.
5. Frequency/store mapping verified - PASS: slot=1636,
   field L=1637, nu=7.369e16 Hz, Lam=40.68 Ang (13.4).
6. Diagnostic exe provenance recorded - PARTIAL: `34f5c1d2`
   provenance = `88dcf220` + the diagnostic diff in 13.2
   (committed with this section); complete-exe provenance for
   `b4686fa1` remains NOT PROVEN.
7. Diagnostic neutrality/unit checks for new code - DONE:
   compile clean, diff verified diagnostic-only, KNOP suites
   pass; a dedicated A/B neutrality run for `34f5c1d2` vs
   `b4686fa1` has NOT been performed (the Q additions are
   write-only diagnostics gated on `KN_Q_RESID`, but neutrality
   is asserted, not yet demonstrated).

The continuation is diagnostic-only by design and must not be
launched until gates 1-4 close; gate 7 ideally also wants an
explicit neutrality demonstration for the new binary.
