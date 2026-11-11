# CMF_FLUX KNOP Convergence Audit

Date: 2026-11-11. Scope: observer-frame (`cmf_flux.exe`) KNOP
non-monotonic transfer, `I_NEXT` lagged-field convergence, and the
reported continuum suppression / narrow-line symptom.

Status tags: [VERIFIED] = read in source or measured in run output;
[SUSPECTED] = consistent with source but needs a runtime measurement;
[BLOCKED] = cannot be tested until a prerequisite exists.

----------------------------------------------------------------------
1. Source files and routine locations
----------------------------------------------------------------------

Population side (`new_main/`):

- `new_main/mod_subs/comp_j_blank.f` - COMP_J_BLANK dispatch. Calls
  `FG_J_CMF_NONMONO_V1` under `DO_NONMONO_KNOP` (~line 888-897),
  `MOM_J_CMF_NONMONO_V1` under `DO_NONMONO_KNOP + USE_LAM_ES=F` and
  under `DO_NONMONO_SOB + USE_LAM_ES=F` (~line 1013+), `FG_J_CMF_V13`
  under `DO_NONMONO_SOB + USE_LAM_ES=T`, stock `FG_J_CMF_V12`/`V13`
  otherwise.
- `new_main/cmfgen_sub.f` - great-iteration driver; KNOP I_NEXT
  maturity gate (`GET_FG_GS_MATURITY`), `DO_KNOP_REFRESH` /
  `MOM_KNOP_REFRESH_V1` calls, `FGINEXT_GS` checkpoint write via
  `SCR_RITE_V2`.
- `new_main/mod_subs/rd_control_variables.f` - VADAT auto-sets:
  `DO_NONMONO_FREQ_MATCH=T` under KNOP and under SOB+LAM_ES=F;
  `LAM_SET=F` override; SOB solver-mode rejections.

Observer side (`obs/`):

- `obs/cmf_flux_sub_v5.f` - main OBS driver. ES loop at ~1498
  (`DO ES_COUNTER=1,NUM_ES_ITERATIONS`); `RD_COHERENT_ES` forces
  `NUM_ES_ITERATIONS=1` at ~1473. `SET_FG_GS_OBS_GRID(NU,NCF)` at
  ~1588-1600 before the frequency loop. KNOP refresh gate at
  ~2215: `IF(DO_NONMONO_KNOP .AND. ES_COUNTER .LT.
  NUM_ES_ITERATIONS) CALL FG_KNOP_REFRESH_V1()`. Moment-path refresh
  gate at ~2228 for `DO_NONMONO_SOB .AND. .NOT. USE_LAM_ES`.
- `obs/COMP_JCONT_V4.INC` - per-frequency OBS solve. KNOP dispatch at
  ~799-822: `SET_FG_NU_NEXT(NU(MIN(ML+1,NCF)))` then
  `FG_J_CMF_NONMONO_V1(...,SIGMA_TRUE,NCF)`. SOB+LAM_ES=F moment
  dispatch at ~903-926: `SET_MOM_NU_NEXT` + `MOM_J_CMF_NONMONO_V1`.
  SOB+LAM_ES=T falls through to stock `FG_J_CMF_V12` + `RJ=TC`.
- `obs/rd_cmf_flux_controls.f` - OBS control reader + consistency.
  KNOP block at ~582-644 enforces `GLOBAL_LINE=BLANK` (or NONE with
  all `TRANS_*=BLANK`), `LAM_SET=F`, `USE_LAM_ES=T`, rejects
  `USE_FIXED_J`/relativistic/DJDT/plane-parallel modes, and
  `INC_GRID=F` (~719-725). `DO_NONMONO_FREQ_MATCH` auto-set at
  ~693-715. **No NUM_ES check anywhere.**
- `obs/check_cmf_flux_param_consis.f` - generic consistency check;
  contains no NUM_ES or KNOP-specific check.
- `obs/cmf_flux_cntrl_var_mod.f` - OBS control variables.

Solver (`plane/fg_j_cmf_nonmono_v1.f`, ~10.4k lines):

- `FG_J_CMF_NONMONO_V1` - formal solver. Lagged stores
  `I_P_DEPTH_NEXT(N_COMP,NP,NCF_IN)`, `I_M_DEPTH_NEXT(...)`;
  `FG_NU_NEXT`; `FREQ_IDX_ITER`; `NM_RAY`/`DLN_RAY` matched-offset
  arrays; `KN_CHI_EXT`/`KN_ETA_EXT`/`KN_FREQ`/`KN_DLN_P`/`KN_DLN_N`
  per-frequency capture for refresh.
- Sweep write-back at ~4707-4878: per-frequency relaxation
  `I_*_DEPTH_NEXT(KC,LS,FREQ_IDX_ITER) = T2 + INEXT_RELAX*(T1-T2)`
  with `NEW_FREQ` snapshot semantics and `FG_GS_MAXCHG`/`FG_GS_NCONV`
  churn metrics.
- Compressing upwind read at ~4096-4128 (I_M) and ~4371-4394 (I_P):
  slot `FREQ_IDX_ITER+NM_RAY(I)`, depth-interpolated between
  `KC1`/`KC2`.
- `SET_FG_GS_OBS_GRID` at ~5569-5588: registers OBS grid, sets
  `GS_OBS_MODE`.
- `RESTORE_FG_INEXT_OBS` at ~6289-6529: validates checkpoint
  (VER>=3, N_COMP/NP/ND, R/V/SIGMA/COMP_DEPTH_IDX fingerprints at
  1e-6), reads saved `KN_FREQ`, remaps `I_NEXT` to OBS grid by
  log-linear interpolation (`TGT=NU_OBS(IS+1)`, `I_FLOOR=1e-40`),
  zeroes slot `NCF_IN`, cold-starts on any failure.
- `DUMP_FG_INEXT_GS` at ~5734-5874: VER=3, writes R,V,SIGMA,
  COMP_DEPTH_IDX, metrics record, `KN_FREQ`, both intensity arrays;
  atomic tmp+rename.
- `FG_KNOP_REFRESH_V1` at ~6730+: per-ray full (space,frequency)
  re-solve over captured `KN_CHI_EXT`/`KN_ETA_EXT`; guards on
  `DO_KNOP_REFRESH`, `FREEZE_I_NEXT`, `HAS_COMPRESSING`,
  `KN_STORE_READY`, `FREQ_IDX_ITER == KN_NCF-1`,
  `SOLUTION_METHOD=='INTEGRAL'`; per-ray worker `KNOP_RAY_SOLVE_V1`
  (~7350+) with `MATCH_FREQ` machinery (`NFM_L`/`DLNM_L`/`DRF_L`,
  ~7353-7356), ordered sweeps, per-pass residual, exit census.
- `SET_FG_FREQ_MATCH` at ~6596.
- `GET_FG_GS_MATURITY` at ~5600-5666: exports `I_NEXT_GS_OK`,
  `FG_GS_CHG_FILE`, `FG_GS_CAP_FILE`, `FG_GS_MAXCHG`, `FG_INEXT_CAPPED`,
  `FG_GS_NCONV`, `FG_SW_MAXCHG`, `FG_SW_NCONV`, `FG_RF_DCHG`,
  `FG_RF_PEND`, stagnation counters - the maturity/residual metrics
  already exist but are only consumed by the population driver.

Pipeline (outside the tree):

- `scripts/prepare_CMF_PARAM.py` - writes `CMF_FLUX_PARAM` from
  `CMF_FLUX_PARAM_INIT`; companion rule forces `NUM_ES>=2`
  (~line 90, 642-644); `cmf_flux_param_overrides` YAML applied via
  `apply_overrides` (~701).

----------------------------------------------------------------------
2. Mathematical description of the current implementation
----------------------------------------------------------------------

The KNOP formal solver implements the comoving-frame transfer equation
along each ray (Knop, Hauschildt & Baron 2009, Paper I):

    dI/ds = eta - chi*I + GAM * dI/d(ln nu)

with signed frequency-drift coefficient (source lines ~3474-3486,
3541-3565):

    GAM   = (V/cR)(1 + SIGMA*mu^2)          [signed; T1*T2 ~line 2820]
    chi_hat = chi + |GAM|/dLOG_NU           [generalized opacity > 0]
    Q     = (GAM/dLOG_NU) / chi_hat         [signed, dimensionless]
    S'    = eta/chi_hat + Q * I_upwind
    dI/dtau_hat = I - S',  tau_hat = chi_hat * dz

`chi_hat` absorbs the implicit part of the frequency coupling
(the `-I_local` part of `dI/dln nu`); `Q*I_upwind` is the explicit
coupling to the neighbouring frequency bin. This split keeps the
optical depth positive for arbitrary `SIGMA` (Paper I Eq. 7
structure) and is unconditionally stable in the sweep direction.

Upwind selection (verified ~3541-3565, reads at ~4096/4371):

- `GAM >= 0` (expanding): `I_upwind = I_PREV`, the previous (bluer)
  frequency computed earlier in the same sweep. `I_P_PREV`/`I_M_PREV`
  hold it; capped at `NONMONO_INEXT_CAP`.
- `GAM < 0` (compressing, `SIGMA_TRUE < -1` reachable): `I_upwind =
  I_NEXT`, the next (redder) frequency NOT yet computed this sweep.
  Read from the lagged depth-point stores
  `I_P_DEPTH_NEXT`/`I_M_DEPTH_NEXT`, interpolated in depth between
  `KC1`/`KC2` with weight `W_INEXT`.

Slot convention (verified ~3067-3098, ~4707-4878): during frequency
`ML`, `FREQ_IDX_ITER = ML-1` (the INIT branch skips the `NEW_FREQ`
increment). The sweep writes the fresh intensity of frequency `ML`
to slot `FREQ_IDX_ITER = ML-1`; therefore **slot `s` stores the
intensity of frequency `s+1`**. A compressing cell at frequency `ML`
reads slot `ML-1+NM_RAY(I)` with `NM_RAY >= 1`, i.e. frequency
`ML+NM_RAY` - a redder bin, lagged by one full sweep (or by the last
refresh). The matched-foot remap (`DO_NONMONO_FREQ_MATCH`, forced T
under KNOP) sets `NM_RAY`/`DLN_RAY` per ray point so the read lands
on the bin containing the characteristic foot rather than the
adjacent bin (~3587-3643); the refresh worker carries the same
`MATCH_FREQ` machinery (`NFM_L`/`DLNM_L`/`DRF_L`, ~7353-7356).

Asymmetric spacing: `SET_FG_NU_NEXT(NU(ML+1))` supplies the forward
log-spacing `ln(nu_k/nu_{k+1})` used for compressing cells; the
backward `dLOG_NU` is used for expanding cells (~3578-3585,
`KN_DLN_P`/`KN_DLN_N` capture ~3087-3093).

Refresh (`FG_KNOP_REFRESH_V1`, ~6730+): after a complete sweep it
re-solves every ray's full (depth,frequency) coupled problem on the
captured coefficients `KN_CHI_EXT`/`KN_ETA_EXT` (pre-`chi_hat`
opacities/emissivities on the extended grid), alternating frequency
directions per cell sign, with ordered Gauss-Seidel sweeps, per-pass
residuals, relaxation `KN_OMEGA`, pass cap `KN_RAY_MAXPASS=60`,
tolerance `KN_RAY_TOL=2e-3`, memory budget `KN_RF_TOT_MB=32768`,
and commit vetoes. The converged field is written back into the same
`I_*_DEPTH_NEXT` stores the sweep reads. Guards: `DO_KNOP_REFRESH`,
`FREEZE_I_NEXT`, `HAS_COMPRESSING`, `KN_STORE_READY`,
`FREQ_IDX_ITER == KN_NCF-1`, `SOLUTION_METHOD=='INTEGRAL'`; each
guard failure returns silently after a one-time warning.

Checkpoint (`FGINEXT_GS`, VER=3): written by the population run via
`DUMP_FG_INEXT_GS` (R,V,SIGMA,COMP_DEPTH_IDX fingerprints, metrics,
`KN_FREQ`, both intensity arrays). OBS restore searches cwd then
`../`, requires VER>=3, validates dimensions and fingerprints
(1e-6 relative on R,V,SIGMA; exact on COMP_DEPTH_IDX), then remaps
the field from the model frequency grid to the OBS grid by
log-linear interpolation of intensity with `I_FLOOR=1e-40`
(`TGT=NU_OBS(IS+1)` for `IS=1..NCF_IN-1`, consistent with the
slot=s+1 convention). Any failure zeroes both arrays and reports a
cold start.

Observer ES loop (`cmf_flux_sub_v5.f` ~1470-2235):

    DO ES_COUNTER=1,NUM_ES_ITERATIONS
      ... frequency sweep ML=1..NCF (COMP_JCONT_V4.INC) ...
      IF(DO_NONMONO_KNOP .AND. ES_COUNTER .LT. NUM_ES_ITERATIONS)
        CALL FG_KNOP_REFRESH_V1()
      IF(DO_NONMONO_SOB .AND. .NOT. USE_LAM_ES .AND.
         ES_COUNTER .LT. NUM_ES_ITERATIONS)
        CALL MOM_KNOP_REFRESH_V1()
    END DO

`LST_ITERATION` is true only on the last pass; OBSFLUX is written
during that pass. The refresh is deliberately skipped after the last
sweep because no later sweep consumes the field.

----------------------------------------------------------------------
3. Confirmed defects
----------------------------------------------------------------------

D1. NUM_ES=1 + DO_NONMONO_KNOP=T is accepted and leaves I_NEXT cold.
[VERIFIED] The KNOP block in `rd_cmf_flux_controls.f` (~582-644)
enforces GLOBAL_LINE/LAM_SET/USE_LAM_ES/solver-mode/INC_GRID
constraints but contains no NUM_ES check;
`check_cmf_flux_param_consis.f` has no NUM_ES check at all. With
`NUM_ES_ITERATIONS=1` the gate `ES_COUNTER .LT.
NUM_ES_ITERATIONS` is never true, so `FG_KNOP_REFRESH_V1` is never
called: one sweep runs on whatever `I_NEXT` the checkpoint restore
left (zeros if absent/rejected). The only protection is the Python
pipeline (`prepare_CMF_PARAM.py` companion `NUM_ES>=2`, ~line 642);
a hand-written or template-derived `CMF_FLUX_PARAM` bypasses it.
`RD_COHERENT_ES` also forces `NUM_ES_ITERATIONS=1` (~1473), so a
coherent-ES KNOP request silently produces the same cold single
sweep. This is the cleanest confirmed route to the reported symptom.

D2. SOB populations have no KNOP checkpoint; KNOP OBS render of a
SOB model always cold-starts. [VERIFIED] `I_P_DEPTH_NEXT` is
allocated only inside `FG_J_CMF_NONMONO_V1`. A `DO_NONMONO_SOB +
USE_LAM_ES=T` population run uses `FG_J_CMF_V13` and never calls it,
so `DUMP_FG_INEXT_GS` returns early and no `FGINEXT_GS` is written.
Confirmed on disk: `cmfgen_validation/m15_2p82d_VAL_sob_jfnk/
cmfgen/2.82days/` contains no `FGINEXT_GS`. Its `obs/` run
(`DO_NONMONO_KNOP=T`, `NUM_ES=3` after a YAML override of the
companion 2) therefore ran the KNOP formal solve on a fully cold
`I_NEXT`, relying on the two refresh calls (after ES passes 1 and 2)
to build the field on the 159660-point OBS grid.

D3. KNOP+USE_LAM_ES=F populations have no matched OBS renderer.
[VERIFIED] Under `DO_NONMONO_KNOP` the OBS reader forces
`USE_LAM_ES=T` (`rd_cmf_flux_controls.f` ~627-633, STOP otherwise).
But a population run with `DO_NONMONO_KNOP=T` + `USE_LAM_ES=F`
converges its populations on `MOM_J_CMF_NONMONO_V1` (the moment
solver), while every OBS KNOP render uses `FG_J_CMF_NONMONO_V1`
(the formal solver). The two are different discretizations of the
same physics; at convergence they should agree, but nothing in the
OBS path checks or enforces that they do. This population
configuration is effectively unrenderable by its own solver.

D4. Observed continuum collapse on a real run. [VERIFIED]
`m15_2p82d_VAL_sob_jfnk/cmfgen/2.82days/obs/OBSFLUX`
(KNOP=T, NUM_ES=3, BLANK, cold I_NEXT per D2) vs
`obs_blank/OBSFLUX` (stock V12/V8, KNOP=F, SOB=F, BLANK, NUM_ES=2),
same populations, run 2026-11-10: median observed intensity
4.63e4 vs 4.41e5 (Jy-ish units of the WRITV block), maxima
2.16e7 vs 2.62e7. The KNOP render's typical flux is ~9.5x lower
than the stock render while line peaks are comparable - the
continuum is suppressed relative to the lines, matching the reported
symptom. (Flux-unit caveat: values are the raw OBSFLUX Janskys
block; the ratio is what matters.)

----------------------------------------------------------------------
4. Suspected defects requiring runtime evidence
----------------------------------------------------------------------

S1. Refresh sufficiency on the OBS grid. [SUSPECTED] Each
`FG_KNOP_REFRESH_V1` call re-solves the per-ray coupled problem on
the captured coefficients to `KN_RAY_TOL=2e-3` with up to 60 passes.
Mathematically one converged refresh should suffice for fixed
coefficients; NUM_ES=3 gives two refreshes under changing ES
scattering sources. Whether the OBS run's field actually converged
is unknown: the residual metrics (`FG_RF_DCHG`, `FG_GS_MAXCHG`,
`FG_SW_NCONV`, exit census) are computed but never printed by the
OBS path (`GET_FG_GS_MATURITY` is only called from the population
driver). The D4 collapse is consistent with either (a) refresh not
reaching tolerance on the fine OBS grid, (b) refresh skipped by a
guard (`KN_STORE_READY`, `FREQ_IDX_ITER==KN_NCF-1`,
`SOLUTION_METHOD`), or (c) an additional defect beyond I_NEXT
maturity. Needs the diagnostics exported to OBS output.

S2. Stale coefficients across ES passes. [SUSPECTED] The refresh
captures `KN_ETA_EXT` during the continuum sweep, which includes
the current electron-scattering source (`ESEC*RJ_ES` for
incoherent ES). Each ES pass updates that source; a field converged
for pass-k coefficients is slightly stale for pass-k+1. Probably
second-order, but unquantified.

S3. BLANK vs SOB line treatment interaction. [SUSPECTED] Under
`DO_NONMONO_KNOP` the OBS reader forces `GLOBAL_LINE=BLANK`
(lines are post-processed, not in the formal solve), so "BLANK looks
better" is partly tautological: the KNOP formal solve only ever
sees continuum opacity/emissivity. The perceived improvement of
BLANK-without-KNOP is a solver-path difference (stock V12/V8 vs
KNOP), not a line-treatment effect per se. The SOB+KNOP combination
is rejected by the consistency checks, so no direct test exists.

S4. Edge slot `NCF_IN`. [SUSPECTED, minor] The sweep never writes
slot `NCF_IN` (max `FREQ_IDX_ITER` is `NCF-1`), and
`RESTORE_FG_INEXT_OBS` explicitly zeroes it. Compressing cells at
`FREQ_IDX_ITER = NCF_IN-1` with `NM_RAY >= 2` (matched foot at the
red grid edge) read slot `NCF_IN` = 0 - a cold coupling for a small
set of edge cells. Bounded volume; worth a diagnostic count, not
necessarily a fix.

S5. Frequency-grid remap fidelity. [SUSPECTED] When a valid
`FGINEXT_GS` exists, OBS remaps by log-linear interpolation of
intensity between the model grid and the (much finer, differently
spaced) OBS grid. Interpolating a peaked intensity field in log
space can misplace line-core photons; for continuum-dominated
cells it should be adequate. No converged KNOP population exists to
test this end-to-end ([BLOCKED], see S6).

S6. No converged KNOP population model. [BLOCKED] All KNOP
validation runs in `cmfgen_validation/` stopped non-convergent
(`m15_2p82d_VAL_knop_formal` 75%, `m15_2p82d_FRESH_knop` 63%,
`m15_2p82d_FRESH_knop_d1a` 75%; others have no convergence record).
The intended production path - KNOP populations -> KNOP OBS render
- has never been exercised on a converged model, so the warm-start
checkpoint path (D2's remedy) is untested in production.

----------------------------------------------------------------------
5. Continuum collapse and narrow lines: mechanism
----------------------------------------------------------------------

The discretization splits the frequency-coupling term into an
implicit opacity (`|GAM|/dLOG_NU` inside `chi_hat`) and an explicit
source (`Q*I_upwind`). The pair conserves the redistributed
radiation only when `I_upwind` is the true neighbouring-frequency
intensity.

With `I_NEXT = 0` (cold start) or a stale/under-converged field:

- In compressing cells (`GAM < 0`, `SIGMA_TRUE < -1`), `chi_hat =
  chi + |GAM|/dLOG_NU` can be much larger than `chi` where the
  velocity gradient is steep (interaction regions, deceleration).
  The effective optical depth `tau_hat` is inflated, while the
  modified source collapses to `S' = eta/chi_hat` - the `Q*I_NEXT`
  source that should feed frequency-redistributed photons into the
  cell is absent.
- The cell therefore radiates a strongly diluted source function
  through a much larger opacity: emergent continuum intensity is
  suppressed. This is an artificial sink created by the lagged
  cold start, not by the Knop formulation - the same `chi_hat`
  with the correct `I_NEXT` recovers the coupling source.
- Lines survive because line `eta` in cores is orders of magnitude
  above continuum; even diluted by `chi_hat`, line peaks remain
  visible. With the frequency-coupling source missing, line photons
  are also not redistributed across the coupled bins, so features
  stay narrow. Result: suppressed continuum + strong, narrow line
  cores - exactly the reported symptom.
- Expanding cells (`GAM >= 0`) are unaffected: `I_PREV` is computed
  earlier in the same sweep and is always available. The collapse is
  therefore spatially selective, concentrated where `SIGMA_TRUE`
  is most negative - consistent with the interaction-region
  phenomenology.

BLANK line treatment makes the effect cleaner to see: the formal
solve carries continuum only, so the collapse acts purely on the
continuum while post-processed lines sit on top unchanged.

----------------------------------------------------------------------
6. Population vs OBS solver comparison
----------------------------------------------------------------------

Population dispatch (`comp_j_blank.f`) vs OBS dispatch
(`COMP_JCONT_V4.INC` + `rd_cmf_flux_controls.f`):

| Population config            | Pop solver                  | OBS flag combo                  | OBS solver               | Match? |
|------------------------------|-----------------------------|---------------------------------|--------------------------|--------|
| KNOP, USE_LAM_ES=T           | FG_J_CMF_NONMONO_V1         | DO_NONMONO_KNOP=T (LAM_ES=T req)| FG_J_CMF_NONMONO_V1      | YES    |
| KNOP, USE_LAM_ES=F           | MOM_J_CMF_NONMONO_V1        | DO_NONMONO_KNOP forces LAM_ES=T | FG_J_CMF_NONMONO_V1      | NO (D3)|
| SOB, USE_LAM_ES=T            | FG_J_CMF_V13 (clamped SIGMA)| DO_NONMONO_SOB=T, LAM_ES=T      | FG_J_CMF_V12 + RJ=TC     | YES*   |
| SOB, USE_LAM_ES=T            | FG_J_CMF_V13                | DO_NONMONO_KNOP=T               | FG_J_CMF_NONMONO_V1      | deliberate diagnostic, cold I_NEXT |
| SOB, USE_LAM_ES=F            | MOM_J_CMF_NONMONO_V1        | DO_NONMONO_SOB=T, LAM_ES=F      | MOM_J_CMF_NONMONO_V1     | YES    |
| stock                        | FG_J_CMF_V12/V13 + V11/V8   | all flags F                     | V12/V13 + V8             | YES    |

*V12 vs V13 are the standard CMFGEN formal/moment pair on clamped
SIGMA; this is the validated SOB production path.

Velocity-gradient convention: both sides pass unclamped `SIGMA_TRUE`
to the non-monotonic solvers and keep `SIGMA` clamped for
infrastructure (`rd_cmf_flux_controls.f` reads `SIGMA_TRUE` for
SOB; `COMP_JCONT_V4.INC` passes `SIGMA_TRUE` to both nonmono
solvers). Frequency-coupling direction and `SET_FG_NU_NEXT`/
`SET_MOM_NU_NEXT` next-frequency convention match. Inner boundary,
ES treatment (`COH_ERENT_ES`/`RJ_ES`), and `DO_NONMONO_FREQ_MATCH`
auto-set rules match the population auto-sets. The remaining
asymmetry is the lagged `I_NEXT` state itself: the population side
converges it over great-iterations with a maturity gate; the OBS
side has only `NUM_ES-1` refresh calls and no maturity gate.

----------------------------------------------------------------------
7. Recommended configuration per population solver
----------------------------------------------------------------------

- KNOP populations (USE_LAM_ES=T, the only OBS-renderable KNOP
  config): `DO_NONMONO_KNOP=T`, `USE_LAM_ES=T`, `GLOBAL_LINE=BLANK`,
  `LAM_SET=F`, `INC_GRID=F`, `NUM_ES>=3` pending S1 diagnostics;
  `RD_COHERENT_ES` must stay F (it forces NUM_ES=1). Requires a
  `FGINEXT_GS` from the population run for warm start; without it
  the render is a cold-start diagnostic.
- SOB populations, production spectra: `DO_NONMONO_SOB=T`,
  `USE_LAM_ES=T`, matched stock-path render (`DO_NONMONO_KNOP=F`)
  - this is the validated recipe.
- SOB populations, KNOP diagnostic: `DO_NONMONO_KNOP=T`,
  `NUM_ES>=3`, treat output as cold-start until S1/S6 resolved.
- SOB+USE_LAM_ES=F populations: `DO_NONMONO_SOB=T`, `USE_LAM_ES=F`,
  moment-matched OBS (currently FAIL on m15 per AGENTS.md 22-Oct
  entry - not viable there yet).
- KNOP+USE_LAM_ES=F populations: no correct OBS renderer exists
  (D3); either extend OBS to allow the moment path under KNOP or
  restrict population runs to USE_LAM_ES=T.

----------------------------------------------------------------------
8. A/B test results and runtime plan
----------------------------------------------------------------------

Already executed (same `m15_2p82d_VAL_sob_jfnk` populations,
2026-11-10):

- TEST A equivalent: `obs_blank/` - KNOP=F, SOB=F, BLANK, NUM_ES=2.
  OBSFLUX median flux 4.41e5.
- TEST C equivalent: `obs/` - KNOP=T, BLANK, NUM_ES=3 (YAML
  override). OBSFLUX median flux 4.63e4 (~9.5x lower), cold I_NEXT
  (no FGINEXT_GS, D2). Line peaks comparable.

Not yet run / blocked:

- TEST B (KNOP, NUM_ES=2) and TEST D (KNOP, NUM_ES=5): cheap to run
  on the same model - directly tests whether more refreshes recover
  the continuum. Requires the OBS run only (no population
  reconvergence).
- TEST E (SOB population-matched): `obs_sob/` configs exist under
  `m15_2p82d_FRESH_sob` and `m15_2p82d_VAL_sob_diag`; needs the same
  OBSFLUX-level comparison.
- True KNOP-population A/B: [BLOCKED] on S6 - no converged KNOP
  population exists. The population convergence problem (63-75%
  ceilings, GS maturity gate) is upstream of this audit.

Runtime plan (ordered):

1. Instrument OBS: call `GET_FG_GS_MATURITY` after the ES loop and
   print `I_NEXT_GS_OK`, `FG_GS_MAXCHG`, `FG_GS_NCONV`,
   `FG_RF_DCHG`, `FG_SW_NCONV`, exit census to LUER. Zero-cost
   diagnostic; distinguishes "refresh never ran" from "refresh ran
   but did not converge" from "converged but wrong".
2. Re-run `obs/` variants on `m15_2p82d_VAL_sob_jfnk` with
   NUM_ES = 2, 3, 5 (TEST B/C/D) and compare OBSFLUX continuum in
   line-free windows, integrated luminosity (RLUMST block), line
   peak/EW/FWHM, and the new I_NEXT residuals.
3. If the continuum recovers with NUM_ES, the mechanism is
   confirmed as I_NEXT maturity; if it does not, suspect S2/S4 or a
   deeper discretization issue and inspect per-ray refresh exit
   codes.
4. Physical validation ladder (Pillar 8): continuum-only transport,
   single line, overlapping lines, velocity reversal,
   SIGMA_TRUE < -1 compression, multiple resonances per ray,
   coherent/incoherent ES, frequency-grid refinement - each checked
   for luminosity conservation and absence of artificial sink,
   against the Knop Paper I/II formal solution.

----------------------------------------------------------------------
9. Proposed code corrections (small, isolated; none applied yet)
----------------------------------------------------------------------

C1. Fail closed on NUM_ES under KNOP. In
`rd_cmf_flux_controls.f` KNOP block (~after line 633): if
`NUM_ES_ITERATIONS < 2` under `DO_NONMONO_KNOP`, STOP with a clear
message (or auto-bump to 2 with a banner, matching the pipeline
companion rule). ~10 lines. Regression test: unit test asserting
the check fires; OBS run with NUM_ES=1 must stop.

C2. Export I_NEXT residual diagnostics to OBS. After the ES loop in
`cmf_flux_sub_v5.f`, call `GET_FG_GS_MATURITY` and write the
maturity fields to LUER under `DO_NONMONO_KNOP`. ~15 lines.
Regression test: production-object unit test on the getter (it
already exists for the population side); OBS run shows the block.

C3. Adaptive refresh (proposal only, needs C1+C2 evidence first).
Replace the fixed `ES_COUNTER .LT. NUM_ES_ITERATIONS` gate with a
loop that refreshes until `R_NEXT = MAX|I_NEW-I_OLD| /
MAX(|I_NEW|,|I_OLD|,I_FLOOR)` over all rays/freqs/depths falls below
a tolerance or a sweep cap is hit. `I_FLOOR` must be justified
physically (the existing `1E-4*ray_peak` significance floor used by
`FG_GS_NCONV` is the natural candidate). This is a behaviour change;
requires the Pillar-8 validation ladder first.

C4. Documentation fix: state explicitly that `DO_NONMONO_KNOP=T`
OBS on non-KNOP populations is a cold-start diagnostic, and that
KNOP+USE_LAM_ES=F populations have no matched OBS renderer (D3).
Obsidian vault update per AGENTS.md section 1.3.

C5. (Minor, optional) Count reads of slot `NCF_IN` (S4) via the
existing `FG_NOSLOT` census machinery; decide whether the red-edge
cold value needs the checkpoint's last-slot remap instead of zero.

----------------------------------------------------------------------
10. References
----------------------------------------------------------------------

1. Knop, Hauschildt & Baron (2009), "A new formal solution of the
   radiative transfer equation in arbitrary velocity fields",
   A&A 496, 295, arXiv:0902.1876, DOI 10.1051/0004-6361:200810290.
   - Paper I: the chi_hat/Q/S' split and sign-dependent upwind
     (Eq. 7 region) implemented in `FG_J_CMF_NONMONO_V1`.
2. Knop, Hauschildt & Baron (2009), "Comoving-frame radiative
   transfer in arbitrary velocity fields II", A&A 501, 813,
   arXiv:0909.4264, DOI 10.1051/0004-6361/200911793.
   - Paper II: Gauss-Seidel iteration of the lagged frequency
     coupling (Sect. 3.2, Eqs 25-28) - the algorithm
     `FG_KNOP_REFRESH_V1` implements as the post-sweep re-solve.
3. Dessart, Audit & Hillier (2015), "Numerical simulations of
   superluminous supernovae of type IIn", MNRAS 449, 4304,
   arXiv:1503.05463, DOI 10.1093/mnras/stv609.
   - Sobolev treatment of non-monotonic fields (Sect. 3.1) behind
     `DO_NONMONO_SOB` and `SIGMA_TRUE`.
4. Ceulemans, De Ceuster & Decin (2025), JQSRT 343, 109470,
   arXiv:2505.08309 - matched-frequency upwind read
     (`DO_NONMONO_FREQ_MATCH`).

The discretisation correspondence claimed above is verified
against the source lines cited in section 1; the S_tilde removal
(07-Sep-2026) is documented in the routine header (~3493-3539) as
the inherited CMFGEN nu-form truncation, not an omission.
