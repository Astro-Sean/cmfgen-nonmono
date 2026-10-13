# VADAT keys

Control-file reference for the non-monotonic and JFNK machinery.
Defaults in parentheses. All paths default OFF; unset flags give the
upstream binary.

## Master flags

| Key | Default | Role |
|---|---|---|
| `DO_NONMONO_SOB` | F | Dessart-2015 Sobolev non-mono line treatment |
| `DO_NONLOC_SOB` | F | RH78 non-local Sobolev coupling (needs SOB) |
| `DO_NONMONO_KNOP` | F | master flag: Knop formal solver + gates |
| `DO_JFNK` | F | arm JFNK (requires `DO_NONMONO_SOB` or `DO_NONMONO_KNOP`) |
| `USE_JFNK` | F | deprecated alias for `DO_JFNK`; warns |

Mode prerequisites are hard-enforced at read time (03-Oct-2026; the
line-treatment gates were warnings before):

- `DO_NONMONO_SOB=T` requires `GLOBAL_LINE=SOB` (or `NONE` with every
  present ion at `TRANS_*=SOB`) and `USE_LAM_ES=T`. `LAM_SET=T` warns.
- `DO_NONMONO_KNOP=T` requires `GLOBAL_LINE=BLANK` (or `NONE` with
  every present ion at `TRANS_*=BLANK`), `LAM_SET=F`, `INC_GRID=F`,
  and no `USE_FIXED_J`/relativistic/plane-parallel modes.
  `DO_BA_RES_ACC=F` warns.
- Both flags `T` STOP in the main VADAT population solver. CMF_FLUX
  OBS no longer supports the Knop formal-J extension; `prepare_CMF_PARAM.py`
  rejects Knop population snapshots at OBS preparation. SOB OBS uses the
  stock formal solver on clamped SIGMA with `DO_NONMONO_SOB=T` only for
  Sobolev EW estimation and profile rendering.
- `GLOBAL_LINE=NONE` needs a `TRANS_*` key for every species slot in
  VADAT; a missing key aborts in GET_KEY_STRING before these checks.
- The OBS reader `rd_cmf_flux_controls.f` validates the SOB line mode
  and the `DO_SOB_FLUX_OBS` prerequisites.

## Knop field solve

| Key | Default | Role |
|---|---|---|
| `USE_LAM_ES` | -- | T: J from FG ray formal solve; F: MOM moment solver |
| `DO_KNOP_REFRESH` | =KNOP | post-sweep lagged-store convergence pass |
| `KN_RAY_ORDERED` | T | ordered-march refresh (SWEEP_FLD_ORD) |
| `KN_ORD_AITKEN` | F | guarded Aitken delta-squared on the ordered pass iterate (failed m15 A/B; superseded by KN_ORD_AA) |
| `KN_ORD_AA` | F | safeguarded Anderson acceleration of the ordered pass iterate |
| `KN_ORD_AA_M` | 2 | AA secant depth; pairs kept = M+1 (2*NI*NCF*8 B each) |
| `KN_RAY_STALL` | 0.90 | same-direction residual ratio counting as a stall pass |
| `KN_RAY_MAXPASS` | 60 | refresh pass cap per ray |
| `KN_RAY_TOL` | 2e-3 | refresh convergence tolerance |
| `DO_NONMONO_FREQ_MATCH` | F | Ceulemans frequency matching (recommend T on compressing models) |
| `INEXT_RELAX` | 0.5 | under-relaxation on lagged writes |
| `NONMONO_INEXT_CAP` | 1e20 | overflow guard on stored I_NEXT |
| `KN_RAY_PACK_MB` | 2GB/thr | coefficient-pack memory cap |
| `MOM_RF_KEEPNEG` | F | repaired cells keep prior iterate |
| `MOM_RF_GS_DESC` | F | live descending reads (tested, neutral-to-worse) |
| `MOM_GATE_BULK` | F | gate on sweep-drift snapshot, not refresh residual |
| `MOM_GAM_PIN` | T | lagged GAM factors evaluated at current-frequency Eddingtons (PSIPREV=PSI etc.) |
| `MOM_LAG_BND` | T | lagged source weight bounded by row's diagonal excess |
| `MOM_ABS_FLR` | 0.0 | absorption-fraction floor on 1-COH in the MOM diagonal; ~1D-3 anchors pure-scattering rows |
| `MOM_TRIAL_MIN` | 30 | first frozen iter eligible for BA trial |

## Freeze / controller

| Key | Default | Role |
|---|---|---|
| `N_LAM_RWARM` | 4 | forced-Lambda warmup per process (requeue store loss) |
| `GS_RWARM_SKIP` | -- | mature checkpoint skips warmup |
| `GS_REL_STABLE` | 1 (hyb: 2) | consecutive mature sweeps to release freeze |
| `GS_REL_N` | 1 | post-release damped-Lambda window |
| `GS_LAM_CAP` / `GS_LAM_TRANSIENT` | 3 / 2% | churn-scaled Lambda cap while near-mature |
| `GS_HOLD_POPS` | F | strict fixed-pops freeze (snapshot/restore) |
| `GS_FREEZE_LIM` | 0=off | in-model frozen-iteration hard stop |
| `KNOP_STALL_MAX` | 60 | slurm-wrapper stall bound (GS_STATE) |
| `GS_TRIAL_LIM` | 15 | max accepted trials between mature sweeps |

## Population solve / resacc

| Key | Default | Role |
|---|---|---|
| `DO_BA_RES_ACC` | -- | residual-acceptance rollback; **mandatory on KNOP** |
| `BA_TRUST` | 1e4 | initial trust radius |
| `BA_NREJ_LIM` | 20 | consecutive-rejection hard stop; masked retries and LAMBDA follow-up soft rejects each count, so masked-retry runs burn ~4-7 per failed episode - pipeline default 100 for non-mono runs |
| `BA_RES_GROW` / `BA_RES_HARD` | 1.5 / 1e2 | residual-growth reject / hard cap |
| `BA_REJ_FRAC` | 0.25 | masked acceptance: <=frac*ND failing columns roll back per-column, RES2 failure or larger fraction rolls back whole; <=0 restores legacy all-or-nothing |
| `BA_SAT_RAW` / `BA_SAT_LIM` | 1e6 / 0 | saturated-demand demotion; with `BA_DAMP_CAP>0` a raw demand >=`BA_SAT_RAW` instead routes into the damped-full march (auto-arms `BA_SAT_LIM=8` leash). SOB pipeline profile: `BA_SAT_RAW=2e3`, `BA_DAMP_CAP=1.3` |
| `BA_ENTER_RAW` | 1e30=off | raw-demand hold on Lambda->full exit (SOB and KNOP; SOB pipeline profile arms `1e4`) |
| `BA_PROBE_LAM` | 0 | probe full step after N held LAMBDA calls while `BA_ENTER_RAW` holds; effective interval is `MAX(BA_PROBE_LAM,NUM_LAM)`. Also the rd_lam probe cadence during `DO_LAM_IT` (requires `DO_BA_RES_ACC` + `BA_DAMP_CAP`; see [[Lambda iteration]]). SOB pipeline profile: 20 |
| `BA_LUM_GATE` / `BA_LUM_DCHG` | F | legacy accepted no-op keys; inner/outer luminosity ratio and drift are diagnostic-only in shocked models |
| `LAM_ENTER` / `LAM_EXIT` / `LAM_N_STABLE` | LAM_VAL / -- / 1 | mode-controller hysteresis. SOB pipeline profile: 400 / 20 / 3 |
| `NUM_LAM` | -- | Lambda cadence. SOB pipeline profile: 20, paired with `BA_PROBE_LAM=20` |
| `BEG_NG` / `IBEG_NG` / `BW_NG` / `ITS/NG` | profile | Ng schedule: arm below `BEG_NG`% `MAXCH`, not before `IBEG_NG`, over `BW_NG` depths, then every `ITS/NG` iterations. Pipeline default: `5.0 / 15 / 40 / 8`; the fire test still requires four clean committed records |
| `NG_ABSTAIN` | F | per-component hold in `NG_MIT_BAND_OPT_V1`: zero-newest-delta and non-positive extrapolations keep current value instead of vetoing the depth; vacuous attempts return `NG_DONE=F`. Non-mono pipeline default T; also gated by the `NG_CLEAN>=4` fire test (module state, not a key) |
| `FIX_IMP` | T | freeze trace rows during BA via `FIXPOP_IN_BA_V3`; cleared on LAMBDA iterations. Production SOB/KNOP/JFNK profiles override it to F because frozen trace equations suppress needed full-step corrections |
| `FIX_BA` / `N_FIX_BA` | -- / -- | frozen-Jacobian reuse; use 1.0 / 10 (avoid the 15/0 stale-Jacobian plateau) |
| `MAX_LIN` / `MAX_LAM` | -- | linear / Lambda correction caps. Production profiles use `MAX_LIN=1e1` and `MAX_LAM=1e10` (effectively unclamped); saturated or probe full steps are still bounded by `BA_DAMP_CAP` and resacc |

## JFNK

| Key | Default | Role |
|---|---|---|
| `JFNK_WARMUP` | -- | calls before JFNK may fire (hyb: 6) |
| `JFNK_COUPLED` | F | permit the coupled tier (KNOP only; rejected under SOB) |
| `JFNK_FRZ_GATE` | 0=use COR | frozen-tier demand gate % (hyb: 1e3) |
| `JFNK_COR_GATE` | 0=off | coupled-tier demand gate % (hyb: 1e2) |
| `JFNK_GATE_NREQ` | 1 | consecutive sub-gate calls to admit (hyb: 2) |
| `JFNK_LOGX` | T | log-x unknowns |
| `JFNK_PTC` | T | pseudo-transient continuation |
| `JFNK_EPS` | 0 | lower bound on `h*||v||inf`, coupled only (hyb: 1e-3; cpl: 1e-2) |
| `JFNK_MAX_EVALS` | 0=off | eval budget per activation (hyb: 40; SOB zlaw14: 200 - one GMRES cycle can spend ~60 JVPs) |
| `JFNK_M` | 10 | GMRES restart dim (hyb: 30) |
| `JFNK_NITS` | 10 | Newton iterations per activation |
| `JFNK_MAX_RESTART` | -- | baseline/PC restarts per call (hyb: 2) |
| `JFNK_MAX_STAGNATE` | 3 | repeated-failure self-disable |
| `JFNK_DEEPDIAG` | F | ~12-eval diagnostic suite |
| `JFNK_FIXT` | T | hold T in evals (F is an untested free-T path) |
| `JFNK_USE_BANDED` | -- | banded preconditioner (needs NUM_BNDS>=3; untested) |
| `JFNK_ROW_TINY` | 1e-30 | near-null row regularization floor (relative to depth-block max) for the JFNK LU preconditioner; 1e-30 leaves degenerate rows that clip at |z|>1e8 - opt config uses 1e-10 (04-Oct-2026) |
| `JFNK_VERB` | -- | verbose JFNK blocks |
| `CONV_DIAG_LVL` | 2 | diagnostic verbosity |

JFNK admission additionally requires `DO_NONMONO_SOB=T` or
`DO_NONMONO_KNOP=T` plus `USE_LAM_ES=T`. `USE_LAM_ES=F` selects the
MOM path, whose lagged `J_NEXT_FRQ`/`H_NEXT_FRQ` store is not frozen
during trial evaluations. Under SOB only the frozen-J tier is
reachable: `JFNK_COUPLED=T` STOPs at read because the eval sweep does
not recompute SOB line rates.

## Model-setup keys that mattered

- `LSTAR`: ejecta runs need the SN_HYDRO_DATA luminosity (e.g.
  4.265857e8 for m15); `LSTAR=0` + `IT_ON_T=T` collapses T to the
  grey floor and wrecks the start. Keep `IT_ON_T=F`.
- `DIF` + `IB_METH`: `DIFFUSION` when inner tau ~30+ (m15 SOB);
  `ZERO_FLUX` for thin inner boundary (KNOP m15 config).
- `FIX_T=T` + `FIX_T_AUTO=F`: hold hydro T (snapshot runs).
- `FIX_NE=F`: let charge conservation set Ne.
- `COH_ES=F`: non-coherent electron scattering. The Type IIn OBS
  preparation carries this through, enforces `NUM_ES>=2`, and expands
  the electron-scattering wing range.
- `STOP_IF_BP=F`: required - the pipeline's `REL_OBS=F/USE_LAM_ES=T`
  trips the upstream SN check otherwise.
- `INC_CHG` must stay F: `CHG_EXCH_DATA` is not in the repo.

## Related

- [[Diagnostics and OUTGEN]]
- [[Model pipeline]]
- [[Field maturity and GS freeze]]
- [[JFNK admission]]
