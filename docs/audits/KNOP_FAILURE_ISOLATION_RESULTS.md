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
field diagnostics. `q = ||T(FLD)-FLD||` style residuals are
reported per ray in the refresh summary - distinct from `GS_MAXCHG`
(the store-to-store change) and from the solver's own last-pass
residual.

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

## 5. Results for R, M and C - PENDING RUNTIME

To be filled from OUTGEN once the arms produce refresh summaries:
I_P/I_M_DEPTH_NEXT deltas, GS_MAXCHG trajectories, KN_Q_LS residuals,
per-ray pass counts, failed/diverged rays, commit rejections,
luminosity, energy balance, KN_COEF_FP/KN_COEF_FPR/KN_COEF_FP_DRIFT.

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
- The q-residual measures one sweep-map application on the
  converged per-ray field inside `KNOP_RAY_SOLVE_V1`; it does not
  include the cross-ray store commit or the frequency-ladder
  shift, so it is a per-ray map residual, not the full refresh
  Jacobian norm.
- `KN_FREEZE_COEF` freezes refresh coefficients only; the main
  per-frequency formal solve still sees live ETA (with evolving
  RJ). Arm R therefore isolates the lagged-store map T, not the
  entire radiation iteration.
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
