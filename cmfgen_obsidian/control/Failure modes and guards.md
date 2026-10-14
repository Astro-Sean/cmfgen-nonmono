# Failure modes and guards

Every guard in the tree exists because a specific blow-up was observed.
This table is the quick map; details live in [[convergence_issues]]
(numbered sections) and the dated log in `AGENTS.md`.

## Observed failure -> guard -> evidence

| Observed failure | Guard | Evidence |
|---|---|---|
| SE update on a churning field blows up the next solve | `GS_FREEZE` ordering | KNOP_MOM_TFIT it16-19 |
| Release on a lucky marginal sweep re-poisons the field | `N_MATURE`/`GS_REL_STABLE` hold | LAM_dbg it7, twice |
| Deferred correction released at cold-start cap re-poisons | `GS_LAM_CAP` on release + `GS_REL_N` | same event |
| Factor-1e5 committed step -> luminosity 2e10 -> 2e16 | Lambda cap while immature | m15 JFNK seq it2 |
| BA step against zeroed I_NEXT after requeue | `N_LAM_RWARM` warmup | requeue wrecks |
| Unchecked BA step with failed snapshot arm | `ba_no_snapshot_demote` | LAM_dbg it39-48 |
| GS_HOLD_POPS snapshot allocation fails and would let populations move | fail closed with STOP | code-path audit |
| Inner/outer luminosity mismatch in a shocked snapshot | diagnostic-only `lum diag`; never a convergence gate | shock deposition invalidates boundary agreement |
| Demand-pinned marching past the basin cliff | `BA_SAT_*` demotion | LAM_dbg it36-47 |
| Rollback overwriting a resacc-rejected step | `BA_STEP_REJECTED` check | m15 it64 wreck |
| JFNK at demand ~1e10 -> GMRES stall -> MAX_LAM wreck | `JFNK_*_GATE` + `GS_FREEZE` veto | m15 opt it11 |
| ~12 wasted formal sweeps on diagnostics/F(0) dup | `JFNK_DEEPDIAG`, `jfnk_fx_fresh` | JFNK audit B2 |
| FD JVP drowned in iterative-solver noise | `JFNK_EPS` floor on `h*||v||inf` (coupled only) | JVP-vs-BA probe |
| JFNK trial sweep mutates MOM lagged stores | Admission requires `USE_LAM_ES=T` until MOM eval freeze exists | `solve_for_pops.f` gate |
| EDD record mismatch or radiation blowup stops during trial JVP | `JFNK_EVAL_BROKEN` poisons eval residual; driver falls back | `COMP_J_BLANK` / `JFNK_RAD_SWEEP_V1` |
| Field marginal forever -> run never converges | `MOM_TRIAL` + stall watchdogs | FreeNe 150-iter freeze |
| Consecutive rejections loop Lambda/BA forever | `BA_NREJ_LIM` hard stop | stall audit |
| One stiff cell vetoes a mostly-good BA step | masked per-depth rollback (`BA_REJ_FRAC`) | rebound is directional - SOB v2 it509 +50% cap still rebounded 120x |
| Armed step untested + accepted on Lambda follow-up | soft-reject accounting (`BA_NREJ_SOFT`, trust halved) | resacc armed but never fired across all runs |
| Compressing-depth set flaps, re-zeroing the store | exit-side hysteresis (-0.98/+1km/s) | stall audit `a9ab6b7` |
| V/sigma_true drift at fixed R leaves ray grid on stale field | NEW_R_GRID fingerprints V + SIGMA_TRUE vs V_EXT/SIGMA_EXT (1E-6 rel) | mirrors MOM `V_SM_SV`; DO_V_PARAMS path |
| Ray-bracket DO WHILE walks off array (hang) | bounded searches + STOP | stall audit `be1213e` |
| JFNK catastrophic fallback committed saturated Lambda step unfiltered | `JFNK_COR_GATE` defers; fallback routes through capped path | it12 wreck, gate fix `162d39e` |
| `V(1)` used as a global speed in a non-monotonic model | non-mono `VINF=MAXVAL(ABS(V))`; RV-file consistency checks peak speed; line storage and model-summary wind-efficiency diagnostics use `VINF` | legacy audit 05-Oct-2026 |
| Velocity-defined clumping maps repeated `V` values to different radii | `DO_CL=T` STOPs under non-mono; all `CL_LAW` forms use `V` or `V(1)` | legacy audit 05-Oct-2026 |
| `DO_DDT`/`PURE_HUB` assume positive single-valued Hubble velocity | both STOP under non-mono; restart interpolation on `LOG(V)` is not valid | legacy audit 05-Oct-2026 |
| Previous-time-step `J/H` and population data are remapped on `LOG(V)` and `R(V)` | `GET_JH_AT_PREV_TIME_STEP` now requires positive, strictly decreasing old and new `V`; `GET_POPS_AT_PREV_TIME_STEP_V2` through `V5` apply the same check before their `R(V)` remap; both back up the `DO_DDT` rejection | legacy audit 05-Oct-2026 |
| Structure revision can rewrite prescribed non-monotonic `V(R)` | `DO_HYDRO` and a present `DO_V_PARAMS` file STOP under non-mono. `REV_RGRID` was unblocked 13-Oct-2026: `VEL_SPACE` is rejected inside `ADJUST_R_GRID_V4`, the interpolated-sigma `STOP` is bypassed under the flags, `SIGMA_TRUE` is refreshed post-revision, and stored `BAMAT` recall is dropped on `R_GRID_REVISED` | legacy audit 05-Oct-2026; REV_RGRID 13-Oct-2026 |
| NG history can carry or extrapolate `R`, `V`, `SIGMA` as if they were unknowns | non-mono `DO_NG_BAND_ACCEL_V4` freezes the four structure records at current values; populations only are accelerated | legacy audit 05-Oct-2026 |
| Gamma deposition/transport uses `V` as an independent coordinate or signed-`SIGMA` solver | `GAMRAY_TRANS=NONLOCAL` and `RAD_TRANS` STOP under non-mono; `LOCAL`/`ABS_TRANS` remain valid | legacy audit 05-Oct-2026 |
| `VLOC_SHOCK_POWER` places shock deposition in velocity space | `GET_EDEP_SHOCK_POWER` requires positive, strictly decreasing `V` under non-mono when depositing in V-space; repeated speeds cannot select multiple radii. `RLOC_SHOCK_POWER` (>0, internal R units) + `DRLOC_SHOCK_POWER` select a radius-space Gaussian instead, bypassing both V guards; `RLOC<=0` keeps the legacy V path | legacy audit 05-Oct-2026; RLOC 14-Oct-2026 |
| `INC_ADV` backward stencil assumes positive radial outflow | `STEQ_ADVEC_V4` STOPs under non-mono if any `V <= 0`; positive non-monotonic `V` remains valid | legacy audit 05-Oct-2026 |
| Optional physics can evaluate positive-flow formulae pointwise (`EXP(-V_SHOCK/V)`, radial advection) | model setup STOPs under non-mono if any `V <= 0`; repeated positive speeds remain supported | legacy audit 05-Oct-2026 |
| `VTURB` power law normalized signed `V` by `V(1)` | main-model non-mono path uses `ABS(V)/MAXVAL(ABS(V))`; standalone OBS_FRAME uses the same magnitude law; default path unchanged | legacy audit 05-Oct-2026 |
| SN shell-shift utility interpolates `T`, `ESEC`, `SIGMA`, `ETA`, and `CHI` on `V` | `DO_SN_SHELL_SHIFT` now requires positive, strictly decreasing `V`; utility remains homologous-flow-only as documented | legacy audit 05-Oct-2026 |
| `REVISE_RVSIG` can locate connection/refinement points through `V` | `TAU`, `NEWG`, `SCLR`, `MDOT`, and `STOP` now require positive, strictly decreasing `V`; radius-coordinate revisions remain allowed | legacy audit 05-Oct-2026 |
| `REV_RDINR` `IV` maps `R(VEL)` and accepts a replacement velocity sequence | option now requires positive, strictly decreasing input and constructed `VEL`; `OUT_RGRID` also rejects non-positive interpolated `V` before writing `sigma` | legacy audit 05-Oct-2026 |
| `GAMABS` display diagnostic inverts `R(V)` before plotting deposition | `GAM_ABS` now requires positive, strictly decreasing `V`; non-monotonic deposition plots are rejected | legacy audit 05-Oct-2026 |
| Standalone `PAR_OPAC` uses stock CMF formal/moment solvers | utility rejects non-positive or repeated `V`; use main `DO_NONMONO_*` executable for non-monotonic fields | legacy audit 05-Oct-2026 |
| Standalone `CMF_FLUX` could select stock bypass solvers or clamped-`SIGMA` EW output under a non-mono model | `RD_CMF_FLUX_CONTROLS` now mirrors the main Knop solver gates (`USE_LAM_ES=T`, no fixed-J/relativistic/DDT/PP modes, `INC_GRID=F`); `DO_CMF_EW` STOPs under either non-mono mode, and moving non-monotonic `V` requires a non-mono flag | legacy audit 05-Oct-2026 |
| Display `GRV` labels the x axis through an `X(V)` inversion | `MAINGEN` now requires strictly decreasing `V` for `GRV`; non-monotonic models are rejected for this option only | legacy audit 05-Oct-2026 |
| Display `GRED`/top-axis helper labels the x axis through an `X(ED)` inversion | `MAINGEN` and `SET_UPPER_AXIS` now require positive, strictly increasing `ED`; repeated/non-monotonic density is rejected rather than labeled on an arbitrary branch | legacy audit 05-Oct-2026 |
| `REV_HYDRO_FILE` `INT` finds a unique velocity crossing and `TST` logs a signed gradient product | utility now requires positive, strictly decreasing `V` for `INT`, bounded crossing brackets, and a positive `TST` log argument | legacy audit 05-Oct-2026 |
| `ADJUST_ATM_R_GRID` locates its connection point in `LOG(V)` | utility now requires positive, strictly decreasing `OLD_V`; this backs up the `REVISE_RVSIG` `TAU` guard | legacy audit 05-Oct-2026 |
| `TST_CLUMP_LAW` evaluates velocity-parameterized clumping laws | display `TSTF` now requires positive, strictly decreasing `V`; this matches the main `DO_CL` rejection | legacy audit 05-Oct-2026 |
| `COMPUTE_TAU_RADIAL` refines intervals and chooses a cutoff from velocity drops | spectrum plotting now requires strictly decreasing `VEL`; repeated/increasing velocities are rejected before grid insertion | legacy audit 05-Oct-2026 |
| Display `SIGMA` option takes `LOG10(SIGMA+1)`, `LOG(V)`, and a log-velocity derivative | `MAINGEN` now requires positive `V`, `SIGMA > -1`, and strictly decreasing `V` for that option | legacy audit 05-Oct-2026 |
| Display `LOGV`/`XLOGV`, `CAKT`, `XCAKT`, `VOL`, and `CHKV` take logs or fractional powers of velocity/gradient terms | `MAINGEN` now rejects non-positive `V`, non-positive `SIGMA+1`, or a non-positive expansion factor before evaluating those diagnostics; `PLT_RATES`, `PLT_NET`, `PLT_ION_MF`, `PLT_DFR`, `PLT_JH`, and `PLT_JH_CUR` apply the same positive-`V` check to `XLOGV` axes | legacy audit 05-Oct-2026 |
| Legacy callers can pass a repeated or reversed abscissa to a monotonic interpolator | `MON_INTERP`, `MON_INT_FUNS_V2`, `MON_INTERP_SP`, `SP_MON_INTERP`, `LOG_MON_INTERP`, and `LIN_INTERP` validate coordinates before interpolation; `MON_INTERP_FAST` and `_V2` now also check minimum grid lengths, strict source monotonicity, and strict output ordering before their interval walk. `LOG_MON_INTERP` checks its log domain | legacy audit 05-Oct-2026 |
| Hydro `R=R0+Vdt` can reorder the mapped radius grid or produce a non-positive density factor | `SET_RV_HYDRO_MODEL_V3` now requires positive `V`, strictly ordered input/output `R`, and a positive density expansion factor | legacy audit 05-Oct-2026 |
| `AVE_DC` maps populations between models on `LOG(V)` and expands `R` by `Vdt` | utility now requires positive, strictly decreasing `V1`/`V2` and a strictly decreasing expanded `R_NEW` | legacy audit 05-Oct-2026 |
| `MERGE_SCR`/`EXTRAP_SCR` map populations on `LOG(V)` and expand `R` by `Vdt` | utilities now require positive, strictly decreasing input `V` and a strictly decreasing expanded `R`; they also declare the pre-existing missing `L` loop index required by `IMPLICIT NONE` | legacy audit 05-Oct-2026 |
| Stock `CMF_FORM_SOL` uses the clamped signed-`SIGMA` formal treatment | `DO_NONMONO_KNOP` now bypasses that optional last-iteration path in the main and OBS continuum includes; the dedicated Knop solver remains the only formal-J path | legacy audit 05-Oct-2026 |

## The classic limit cycles

```text
freeze/release cycle:
  field matures (lucky sweep) -> BA released -> field re-poisoned
  -> freeze again -> repeat forever
  fix: N_MATURE consecutive-mature hold (GS_REL_STABLE)

Lambda<->BA alternation:
  NUM_LAM=2 alternates 2 Lambda + 1 BA; every BA demands 1e4-1e8%
  on a still-moving field -> forced back to Lambda
  fix: GS maturity gate + churn-scaled Lambda cap

FULL/LAMBDA ping-pong (KNOP opt, 04-Oct-2026 review):
  exit tests cap-saturated applied correction (<LAM_EXIT always true
  under release damping) -> every LAMBDA call releases; every full
  step saturates the cap (>LAM_ENTER) -> demotes; 36+ pairs.
  Fixes: BA_ENTER_RAW armed (holds exit while last full raw demand
  stays high), CNT_LAM now counts held calls so the probe/escape
  re-measurement channels are live, BA_PROBE_LAM sets the probe
  interval (MAX(BA_PROBE_LAM,NUM_LAM) held calls).

FULL->LAMBDA->LAMBDA (SOB):
  each full step demands saturated 1e7% decrease, applied x5,
  residual rebounds ~2.4x
  fix: NUM_LAM=BA_PROBE_LAM=20 (bounded re-measurement while the
       BA_ENTER_RAW demand hold is armed), resacc trust cap,
       BA_REJ_FRAC masked rollback (keep the good columns, hold
       the stiff ones - rebound is directional, not size-driven),
       BA_SAT_RAW=2e3 plus BA_DAMP_CAP=1.3 damped-full leash while
       demand stays saturated

rd_lam plateau park (SOB, 11-Oct-2026):
  stable-exit gate needs LAM_N_STABLE calls < LAM_EXIT; the LAMBDA
  attractor floored ~45% and rebounded, so stable=0/3 forever and no
  probe/escape path exists behind RD_LAMBDA -> permanent LAMBDA park.
  fix: BA_PROBE_LAM-cadence bounded probe inside RD_LAMBDA
  (BA_RDLAM_PROBE/BA_RDLAM_ACCEPT; resacc-judged, BA_DAMP_CAP-capped;
  accepted probe pins CNT_STABLE_LAM so DO_LAM_AUTO hands off)
```

## The bug classes that were NOT math

Every verified numerical failure in the run history traced to plumbing
or setup, not to the Knop equations ([[KNOP_MATH_AUDIT]]
verdict):

- Wrong-frequency fallback (`RSQJNU_PREV` bluer field in negJ repair
  -> it1 L=2.7e20 ratchet; removed `cd45cf2`).
- Shared OMP variable (`DLN_C` across the ray loop; `e2b9bd2`).
- Single-precision capture store underflow -> NaN refresh (KN_M_*
  promoted to LDP).
- int32 counter overflow on ~5e10-cell counts (promoted to INTEGER*8).
- Uninitialized `SUCCESS` on NEWMOD cold start -> never-written
  `BAMATPNT` re-read -> STOP (`7c6a3c3`; latent upstream bug).
- `obs_frame_sub_v9` one-directional frequency hunt assuming monotonic
  V -> 13% corrupted spectrum bins; plus CMF-grid undercoverage at the
  red edge (clamped within 2% at the boundary).
- `'ZERO FLUX'` vs `'ZERO_FLUX'` string typo in `comp_j_blank.f` -
  upstream bug, zero-flux branch never matched.
- `rv_sn_model_sniin.f` SIGMA denominator (`VEL_LAW=12` only;
  `f2f4784`).
- **False convergence via masked demand** (04-Oct-2026): `MAXCH` was
  computed after the resacc trust mask zeroed corrections, so a fully
  masked step reported demand 0 and tripped `MAXCH<EPS` - jfnk_hyb
  exited "converged" at it8 with 69/70 columns masked, se_res=1.1e14,
  L_outer=0. Fixed: masked demand stays in `MAXCH` (see
  [[BA matrix and SOLVEBA]]).
- **Deterministic identical retries** (04-Oct-2026): every
  Krylov-fail / clip / JVP-fail path set a refresh flag that was
  never consumed on reject paths and halved a trust radius that is
  only read inside a completed GMRES - each retry was bit-identical.
  Now aborts immediately (see [[JFNK solver]]).
- **Soft-reject invisibility** (04-Oct-2026): armed steps judged on
  LAMBDA follow-ups degraded but counted only `BA_NREJ_SOFT`, never
  `BA_NREJ_CONSEC` - 20 committed bad steps read `nrej=0, [OK]`.
  Now counts toward `BA_NREJ_LIM`.  Consequence seen 05-Oct-2026:
  with `BA_REJ_FRAC` masked retries plus follow-up soft rejects, one
  failed BA episode burns ~4-7 consecutive rejections, so the default
  `BA_NREJ_LIM=20` kills a masked-retry run after ~5 episodes
  (`m15_2p82d_mask` attempt 1 stopped at it305 still improving).
  Masked-retry runs need `BA_NREJ_LIM` ~100 (non-mono pipeline
  default now, override per config).
- **NG history pollution** (05-Oct-2026): `SCR_RITE_V2` appends a
  record every call, before the residual judgment that may veto it.
  A rejected step therefore leaves its record on disk and the
  rejecting call appends a restored duplicate - the 4-record
  `SCRTEMP` window consumed by the NG band fit holds sign-reversing
  difference pairs, and masked rollback leaves half-restored
  "Frankenstein" records.  NG's vetoes on such windows are correct
  rejections of poisoned secants, not NG defects.  `NG_CLEAN` (calls since last reject/partial/soft-reject/pending-arm/
  JFNK/held-population record) gates firing until the window ages
  clean; `NG_ABSTAIN` holds non-contractive components per-component.
  `NG_DEFERRED` preserves a +1 due target if an auxiliary SCRTEMP
  write steps `MAIN_COUNTER` over it. See [[Acceleration methods]].
- **Wrong-quantity JFNK gate** (04-Oct-2026): admission tested the
  previous call's demand (a LAMBDA demand in the alternating regime),
  not the last non-LAMBDA demand. Now `BA_LAST_RAW` (see
  [[JFNK admission]]).
- **Stale-mode pending arm** (audit, same class): the resacc trial
  marker `BA_RAW_TRIAL_PENDING` was armed under `.NOT. LAM_WAS`, but
  `LAM_WAS` is only assigned at the mode controller further down in
  the same call, so it tested the *previous* call's mode. Every first
  full call after a LAMBDA window (probe, stable exit, NUM_LAM
  escape) skipped the arm, and a rejected probe never reached
  `BA_RAW_REJECTED` - dead code in the alternating regime the flags
  exist for. The arm now mirrors `BAACC_PENDING` directly.
- Warmup fixups snapshot ordering: `FIX_IMPURITY`, `MOD_FIXED_T`, NG
  scheduling computed before the KNOP gates - all now force-set in the
  `FORCE_LAMBDA` block.
- **Trace-population freeze hiding demand** (06-Oct-2026): `FIX_IMP=T`
  lets `FIXPOP_IN_BA_V3` zero BA rows for low-abundance or low-fraction
  stages. On the shocked SOB run this appeared as CORRECTION_SUM bands
  dominated by zeros while the residual stayed high - the full solver
  was not allowed to correct those populations. The SOB/JFNK profiles
  now stamp `FIX_IMP=F`; the switch remains available for explicit
  experiments only.
- Partial clean build left queued OBS jobs pointing at a deleted
  `exe/cmf_flux.exe` (job `12095679`). `build.sh` now builds `make all`,
  keeps the old executables during `clean`, publishes production links
  atomically, and OBS submission preflights the executable path.

## Reading the wreck signatures

- `MAXCH` pinned at the **1e7 sentinel** = saturated demand (a level
  wants >=99.999% annihilation); marching, not contracting.
- `res` flat while `MAXCH` pinned = the resacc cliff - SE_RES2 stayed
  flat through the LAM_dbg blow-up; demand trend is the only
  pre-cliff signal.
- `LUM_RATIO_OI` and its relative drift describe radiative redistribution
  across the shocked snapshot; they are diagnostic and do not gate a step.
- JVP flat ~1e32 independent of eps = outside the linear basin, not a
  JVP bug.

## Related

- [[Field maturity and GS freeze]]
- [[JFNK admission]]
- [[BA matrix and SOLVEBA]]
- [[Open issues]]
