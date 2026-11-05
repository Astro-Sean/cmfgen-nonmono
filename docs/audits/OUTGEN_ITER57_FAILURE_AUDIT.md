# OUTGEN iteration-57 failure audit -- m15_2p82d_VAL_knop_mom

Status: IN PROGRESS (audit). This file is the live record; update in place.

## 0. Reference failure

- Run dir: `/viper/ptmp1/sbrennan/SN2024pba/cmfgen_validation/m15_2p82d_VAL_knop_mom/cmfgen/2.82days/`
- OUTGEN: 19845 lines. Terminal line: `SOLVE_FOR_POPS:       50 no-progress pop-preconditioning iterations reach GS_PRECOND_LIM` at great iteration 57.
- Repo state: `cmfgen_jfnk` at commit `3de42aae` (= `origin/main`, clean; only untracked debug dumps: `MOM_J_NONMONO_ERRORS`, `mom_{jnext,hnext,nm,sm}.txt`, `solver_hist_*.csv`, `unit_tests/MOMJNXT_GS`).

### Model configuration (from VADAT + OUTGEN)

- `DO_NONMONO_KNOP=T`, `DO_NONMONO_SOB=F`, `USE_LAM_ES=F`, `SOL_METH=TRI`, `FIX_T=T`, `DO_JFNK=F`, `DO_LAM_IT=T`.
- `DIF=T`, `IB_METH=DIFFUSION`, `N_TYPE=N_ON_J`, `COH_ES=T`, `CHK_H=MAX_VAL`, `CHK_XM=SET_POS`, `CHK_J=ABS_VAL`, `CHK_NG=T`.
- `VEL_LAW=11` (hydro `SN_HYDRO_DATA` via RDINR), `RD_IN_R_GRID=T`.
- `VFRAC_MOM=1.0`, `VFRAC_FG=0.5` (task text said 20000 -- this run's actual VADAT has 0.5).
- `MOM_RF_DIR=DESC`, `MOM_RF_KEEPNEG=T`, `MOM_H_IMPL=T`.
- KNOP auto-sets: `GS_PRECOND=T`, `GS_PRECOND_LIM=50`.
- Grid: hydro ND=70; RSTAR hydro/photosphere 4.134260e4 -> overridden 1.454688e4 (= innermost hydro radius); RMAX 8.61436267 (= outer/inner ratio; outer R = 1.253121e5, matches OUTGEN). tau_inner ~ 39.3.
- Converter: max velocity jump 1650.28 km/s vs VDOP_MIN ~20 km/s, 18 violations on original grid; refined grid ND_NM=212.
- Velocity profile: V ~ 48.9 km/s at outer edge rising to 7073 km/s at inner edge; sigma from -0.96 to +130 on the 70-pt hydro grid. V > 0 everywhere (zero sign changes) but sigma < -1 pockets exist on the refined grid (worst cells at it57: sigT = -1.04, -1.26, -1.46, -1.77).

### Failure timeline (OUTGEN)

- it 13 (handoff_smoke/knop_mom, same family): mom sweep `d=1.426E+23`, negJ 142706, quar 270590, pquar 157551, capped 127884; `fg d=1.0, bad=39M cells` -- fields unhealthy from early on.
- Every great iteration: `Error --- J(mom) and J(ray) differ by more than 1000% for last frequency` (first at OUTGEN line 1565).
- Every iteration: `Moment refresh stalled` or `diverged`; `[PRECOND] Moment field unsafe - catastrophic drift; hold disabled.`
- it54: `Moment refresh diverged at pass 16 (residual 1.6104E+06 versus the best 6.2062E+02)`; best-pass field rejected (residual runaway); pre-refresh field retained.
- it56: `stalled at pass 21 (residual 1.3740E+04)`.
- it57: `stalled at pass 11 (residual 9.6417E+03)`; then `SOLVE_FOR_POPS: pop preconditioning: demand=2.756E+06 > limit=1.000E+03; hold postponed, live LAMBDA n_precond=56` -> `50 no-progress pop-preconditioning iterations reach GS_PRECOND_LIM` -> STOP "the population demand never fell below GS_PRECOND_MAX; stopping as non-convergent".

### Sweep-report numbers (it56/57, representative)

- `drift lagged J_NEXT drift 1.374E+04 (2490230 cells >5%, 470849 >50% of 11951924 total) [BAD]`
- `repairs 4561 significant negative J_NEXT cells, 3321 round-off negatives, 452145 quarantined, 426076 physically quarantined`
- `solver 0 failures, 0 non-finite, 2505 runaway cells, 0 non-positive TB [OK]`
- `structure 0 non-monotonic rows; 10972861 forward-epsilon midpoints; largest |EPS|/HU 6.883; epsilon-cap hits 39556688`
- `lagbnd 6775206 boundary-lag cells, smallest lag 0.000`
- `matrix smallest DR 3.115E-01 at (22,1375)` -- matrix not badly conditioned
- `eq resid largest equation residual 2.148E-05 at (2,224)` -- assembled system is solved accurately
- `params largest |GAM| 9.971E+01; largest |sigT| 1.569E+02; 0 boundary HBC floors`
- refresh per-pass (it57): ~4400 negJ, ~6500 quarantined, ~4700 capped, 901338 freq matches of 11052318 compressing cells, 0 off-grid reads, largest bin offset 90, eps-limit ~1.5M cells, owned-skip ~1.3M rows, gam-floor 133732 cells.

## 1. Verified findings so far

### F1 -- First negative-J mechanism (from MOM_J_NONMONO_ERRORS in the run dir)

The per-cell dump shows cells whose **assembled RHS (XM_SAVE) is already negative**
before the solve, at compressing depths where ETA ~ 1e-133 (zero local source):
e.g. depth 39 @ FREQ=200: XM_SAVE(39) = -2.6e-26 while XM(39) solved -8.1e-27;
sigT=-2.08, absorp term =2.89 of TB=3.12. Mechanism: signed lagged coupling
(H-divergence `VB*HUP(I-1)+VC*HUP(I)` and/or the EPS scatter debit) subtracts
from a zero source. `MOM_EPS_LIM` caps only the EPS-scatter debit; the
H-divergence debit and the compressing upwind RHS are NOT positivity-capped, and
`MOM_LAG_BND` bounds magnitude vs diagonal excess, not RHS sign.
-> The negative J is a *negative RHS* problem at source-free compressing cells,
not (primarily) a tridiagonal-solve or pivoting failure. Class (A) equation
assembly / boundary-source treatment, per the task's taxonomy.

### F2 -- Diagnostics that misreport (verified in source)

- `PASS_NOFF` ("off-grid reads") is vacuous: `MOM_NM_J` is initialized to 1 and
  only ever set to `N_CAP >= 1`; `COUNT(MOM_NM_J==0)` is identically 0. Cells
  with `FREQ_INDX+MOM_NM_J > NCF_IN` silently drop the upwind read (JUP=0) --
  the reddest ~90 bins under DO_NONMONO_FREQ_MATCH. mom_j_cmf_nonmono_v1.f
  ~2554-2570, 5414-5417.
- `MOM_SW_FCMP` (printed "compressing-row share 0.000") is only assigned in the
  AUTO fall-through branch of the dominance block (line ~5213). With the
  VADAT-forced `MOM_RF_DIR='DESC'` it stays at its initial 0 -> misleading
  print only; dominance itself is correct.
- The `J(mom) vs J(ray)` check (comp_j_blank.f ~1329-1344) compares ONLY the
  outer boundary cell at the LAST frequency; the metric saturates at 200%
  (relative). The ">1000%" message fires at T1>100 -- the same threshold as the
  ">100%" branch at line 592-597, so the threshold text is wrong AND the
  diagnostic cannot distinguish one boundary cell from widespread mismatch.
  It does NOT compare interior depths or other frequencies.

### F3 -- Divergence bookkeeping

- Runaway rejection caps a cell at 1E3 x previous-valid per pass; an unstable
  cell can therefore still grow 1E3/pass -> pass-16 residual 1.6E6 after best
  6.2E2 is consistent with ~2 explosive cells, not global divergence. The abort
  at `PASS_CHG > 1E3*MIN_RES` (line ~5961) then discards the pass sequence.
- `MOM_LAG_BND`/`MOM_EPS_LIM` modify the solved equations on millions of cells
  (39.5M eps-cap hits/sweep; 6.8M lagbnd cells at lag factor 0.000). The map
  being iterated is therefore NOT the Knop equations on a large cell fraction;
  the reported "contraction" of the bounded map is not the physical one.

### F4 -- Quarantine

- `MOM_RF_KEEPNEG` freezes XM_BAD cells to previous J_NEXT_FRQ (line 3459+).
  Age `MOM_QUAR_AGE` increments on sweep calls only (not refresh); persistent =
  >=3 sweeps. Sweep report: 452145 quarantined / 426076 persistent at it56.
  Frozen cells are excluded from drift metrics (T5 path) while still being
  consumed as upwind fields by neighbours -- i.e., the solver reports a field
  that contains ~3.8% permanently frozen cells as if it were the converged
  field. This matches the task's "stale J_NEXT persists while apparent progress
  is reported" hypothesis -- verified.

## 1b. Controlled reproduction and verified correction (moment-audit)

### First numerical failure (verified)

The first invalid operation in the production radiation-field path is the
**assembly of a negative right-hand side** at interior rows satisfying all of:

- compressing (`GAM_SIGNED < 0`, `sigT < -1` pockets on the expanded grid),
- ~zero local source (ETA as small as 1E-133 in the run dump),
- absorption-dominated (`|TA|+|TC| << TB`, e.g. 0.23 vs 3.12 at the run's
  depth-39 record; `COH ~ 0.008`),

where the net lagged debit (`HS(I)*H_up(I) - HS(I-1)*H_up(I-1)` H-divergence,
EPS scatter, upwind-J term) exceeds the local source. `MOM_EPS_LIM` bounds only
the EPS-scatter debit; `MOM_LAG_BND` bounds lagged magnitude against diagonal
excess, not RHS sign. `MOM_H_IMPL=T` was already on and only damps the
diagonal. DGTSV then returns negative J at the weakly-coupled rows ->
`XM_BAD_NM` repair -> `MOM_RF_KEEPNEG` freeze -> quarantine (`452145`
quarantined, `426076` persistent at it56). Frozen cells are excluded from the
drift metric yet still feed neighbours as upwind data, so the field read as
`>50% drift on 3.9% of cells` -> `MOM_UNSAFE_B` -> `[PRECOND] Moment field
unsafe - catastrophic drift; hold disabled.` -> the `GS_PRECOND` gate
(demand 2.76E6 > 1E3) postpones the population hold -> 50 no-progress calls
-> `GS_PRECOND_LIM` STOP.

Evidence chain: run-dir `MOM_J_NONMONO_ERRORS` records XM_SAVE < 0 *before*
the solve; controlled driver mode `seed` in `unit_tests/f90/drv_mom_solve.f`
injects the same state (zero in-band source, thick compressing band, negative-
divergence upwind H) and produced `SEED_XMSAVE_MIN = -1.26E-09` with
`NEGJ_T = 41` negative solved J at the band rows.

### Correction C3: MOM_XM_POS (implemented, tested)

`MOM_XM_POS` (VADAT, default TRUE): after `MOM_LAG_BND`, cap the net lagged
debit at each interior row's local source balance so the assembled RHS stays
non-negative. Same donor-balance policy as `MOM_EPS_LIM` and the catastrophic
rebuild path (the only two places a bound already existed); F restores the old
path for A/B tests. Counted by new counter `MOM_XMNEG`, printed on the sweep
`lagbnd` report line and per refresh pass.

Files: `plane/mom_j_cmf_nonmono_v1.f` (module flag/counter, bound loop after
the LAG_BND block, per-pass reset, report lines, `SET_MOM_LAG_OPTS` arg),
`new_main/mod_subs/rd_control_variables.f` (VADAT read),
`new_main/mod_subs/control_variable_mod.f` (declaration).

Numerical evidence (seed mode): `SEED_XMSAVE_MIN = -1.26E-09 -> 0.0`,
`NEGJ_T = 41 -> 0`, solved J positive (~1E-134) at the dark absorbing band
rows. Existing suite `test_mom_solve.py`: 98/98 pass (bound is a no-op on
cases with XM >= 0); `test_mom_sweep_table.py` 8/8; `test_knop_ray.py` 10/10;
`test_vadat_keys.py` 77/77.

### Also fixed (diagnostic-only, verified defects)

- C1: `PASS_NOFF` was vacuous (`MOM_NM_J==0` can never fire). Now counts
  compressing rows whose lagged upwind slot `L+MOM_NM_J > KN_M_NCF` falls off
  the red edge (silent JUP=0 reads).
- C2: `MOM_SW_FCMP` printed stale 0.000 under a VADAT-forced `MOM_RF_DIR`;
  the compressing-share census is now computed unconditionally.

## 2. Still open (investigation continues)

- [ ] Whether MOM_XM_POS removes the it57 signature on the real model (C7
  in the TODO): negJ/quar census and refresh residuals should drop;
  if the 620->1.6E6 residual runaway persists it is a separate
  amplification channel (runaway cap 1E3/pass permits 1E3 growth per
  pass on an unstable row).
- [x] The `XM_SAVE -> boundary rows` restore at ~3219 is NOT a defect:
  it puts back the original assembled RHS at rows the rebuild loop
  skipped (boundary rows), before the second DGTSV. Correct.
- [x] solve_for_pops chain traced: `GS_FREEZE` holds pops on unsafe
  field, but `GS_PRECOND` postpones the hold whenever demand
  MAXCH > GS_PRECOND_MAX (1E3); demand 2.76E6 -> live LAMBDA on the
  unsafe field; N_PRECOND_TOT counts calls without >10% demand
  improvement; 50 -> GS_PRECOND_LIM STOP. Design is defensible
  (holding at 1E6 demand is worse), but it means an unsafe field is
  consumed live for up to 50 iterations -- the real fix must make the
  field safe, which is what MOM_XM_POS addresses upstream.
- [x] J(mom)/J(ray) compare: single boundary cell at the LAST
  frequency only; ">1000%" text fires at the same T1>100 threshold
  the stock branch labels ">100%" (mislabelled copy); cannot
  distinguish one cell vs widespread. FIXED (C4): all-depth census
  (JCMP_W_REL/ABS/NU/I, JCMP_N50) accumulated per LST_ITERATION
  sweep, worst cell + >50% count appended to J_COMP and printed to
  LUER at FREQ_INDX=NCF; label corrected to 100%. Thresholds
  unchanged.
- [x] C6 persistent-quarantine gating: MOM_NPQUAR (cells quarantined
  MOM_QUAR_PERSIST=3+ consecutive sweeps) was fetched and printed but
  never gated. Added to MOM_UNSAFE_B at >1% of cells and to
  MOM_FREEZE_B at >0.1% in both the RF-converged and non-converged
  gate branches. A persistently quarantined cell is consumed as
  upwind data while excluded from the drift metric; >1% permanent
  holes is a structural failure, not a marginal field.
- [x] T4 population chain evidence: [PRECOND] unsafe fires from
  iteration 1 (field never safe); exempted demand fluctuates
  non-monotonically 3E15 -> 5E3 -> 6E8 -> ... (CaI:4s21Se[0]
  depths 64-65 inner, then depths 1-5 outer, FeI/Sk2 trace levels);
  n_precond counts each live-LAMBDA call, 50 -> GS_PRECOND_LIM STOP.
  BA consumes the same lagged field the formal solve wrote (RJ store)
  -- consistency mechanism present, the field itself was unstable.
- [x] V5 full unit sweep: 30/30 test files pass post C4/C6.
- [x] C7 submitted: job 12719316, warm restart from SCRTEMP it57 on
  the rebuilt cmfgen_dev.exe (MOM_XM_POS=T default). OUTGEN and
  MOM_J_NONMONO_ERRORS preserved as *.it57_orig for comparison.
  Awaiting result.
- [ ] solve_for_pops chain: how demand 2.76E6 persists for 50 preconditioning
  iterations; which levels/depths dominate; whether BA consumes the unsafe field.
- [ ] Restart integrity of J_NEXT_FRQ / MOM_NU_GRID / MOM_QUAR_AGE.
- [ ] Confirm whether J(mom)/J(ray) mismatch is widespread or one boundary cell.
- [ ] Numerical tests A-J.

## 3. Source map (verified line numbers, commit 3de42aae)

- `plane/mom_j_cmf_nonmono_v1.f`
  - module state / counters: 90-940
  - solver entry `MOM_J_CMF_NONMONO_V1`: 941; grid rebuild: 1049-1400
  - freq-match tables (MOM_NM_J/H, MOM_DLN_J/H, MOM_NU_GRID): 2034-2114
  - coefficient assembly (GAM_FAC/CHI_HAT/W/EPS/DTAU/PSI/HU/HL/HS): 2124-2355
  - interior rows TA/TB/TC/XM + EPS upwind + MOM_H_IMPL: 2381-2442
  - compressing upwind read JUP_NM/HUP_NM (+refresh snapshot logic): 2516-2610
  - lagged RHS + EPS_LIM limiter + implicit rescale: 2612-2746
  - MOM_LAG_BND row bound: 2748-2788
  - outer BC (Eddington, A_BND direction) vs negated-row BC: 2790-3026
  - DGTSV solve + fallback re-solve: 3037-3236
  - negative-J classify / repair / runaway reject: 3238-3391
  - J_NEXT_FRQ write (relax, KEEPNEG, OWNW, suppression): 3393-3630
  - H reconstruction + H_NEXT_FRQ write: 3691-3828
  - per-freq KN_M_* capture: 3851-3878; GS checkpoint dump/report: 3906-4099
  - `MOM_KNOP_REFRESH_V1`: 4828-6639 (guards 4921-4958; pass loop 5238-6041;
    best-pass 5505-5567; cycle/damp/AA 5569-5868; exits 5917-6040;
    commit/rollback 6059-6186; outer AA 6188-6421; retained-delta metric
    6434-6478; adaptive JNEXT_RELAX 6529-6622)
- `new_main/mod_subs/comp_j_blank.f`: J(mom)/J(ray) boundary compare 540-608,
  1294-1345.
