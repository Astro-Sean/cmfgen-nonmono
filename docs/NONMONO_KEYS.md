# NONMONO_KEYS.md - VADAT and OBS key reference (non-monotonic fork)

Per-key reference for the control keys added or modified by this
fork. Sources of truth: `cmfgen/new_main/mod_subs/control_variable_mod.f`
(declarations), `cmfgen/new_main/mod_subs/rd_control_variables.f`
(parsing + auto-sets), `cmfgen/new_main/mod_subs/check_param_consistency.f`
(rejections). `docs/CHECKLIST.md` is authoritative for current
status; dated notes below are historical records of when/why a
behaviour landed, not a live status board.

Auto-set semantics: absent-key auto-sets apply only when the key is
omitted from VADAT; explicit values win except the structural
overrides (`LAM_SET`, the KNOP triple force) which apply even over
explicit input and print a banner.

## 1. Feature flags and gating

- `DO_NONMONO_SOB` (default F) - Dessart, Audit & Hillier (2015,
  arXiv:1503.05463, Sec 3.1) Sobolev treatment for non-monotonic
  velocity fields. Requires `GLOBAL_LINE='SOB'` (or `NONE` with
  every present ion `TRANS_*=SOB`). Rejects solver modes that
  bypass the COMP_J_BLANK non-mono branch: `USE_FIXED_J`,
  `USE_FORMAL_REL`, `USE_J_REL`, `USE_DJDT_RTE`, `PP_MOD`,
  `PP_NOV` (STOP coded 25-Nov-2026; documented as rejected since
  22-Oct-2026 but the check was never implemented until the
  production-consolidation audit). `LAM_SET` is normalized to F
  (hard-STOP backstop) since 22-Oct-2026.
- `SIGMA_TRUE` (`MOD_CMFGEN`): unclamped `dlnV/dlnR - 1` from
  centred differences on R,V. `SIGMA` keeps its infrastructure
  clamp so the moment equations and continuum formal solver are
  untouched. `RECOMPUTE_SIGMA_TRUE`
  (`cmfgen/new_main/subs/recompute_sigma_true.f`) fills it at model setup,
  after `DO_WIND_VEL_V1`, after NG acceleration, and after
  `ADJUST_R_GRID_V4` grid revision. Under `DO_NONMONO_SOB` it also
  re-clamps `SIGMA` to `>= -0.999` (the model readers leave `SIGMA`
  unclamped, which let `sigma < -1` reach stock routines that
  assume `sigma >= -1`, e.g. `MOM_J_CMF_V11`'s `W`/`EPS` terms).
  Not applied under `DO_NONMONO_KNOP` (26-Sep-2026, commit
  5c85085).
- Tau floor in `sobjbar_sim.f`, `sobew.f`, `sobew_grad_v2.f`,
  `sobew_grad_v3.f`: for sigma < 0, `|1+sigma*mu^2|` is replaced
  by `max(|1+sigma*mu^2|, -0.2*(1+sigma), 1e-4)`. Bit-identical
  to upstream for sigma >= -1 + 1e-4; in the window
  -1 < sigma < -1 + 1e-4 the 1e-4 arm can bind along near-radial
  rays (where |1+sigma*mu^2| < 1e-4) even though the flow there is
  still monotonic.
- `SUB_SOB_LINE_V3` and the dielectronic `SOBEW` call in
  `cmfgen_sub.f` pass `SIGMA_TRUE` when the flag is set.
  `SUB_SOB_LINE_V3` carries transient-iterate guards on the
  BA/Jacobian terms (26-Sep-2026, commit b107b4b): `VB_SIM`/`VC_SIM`
  and `SRAT` are set to zero when `ETAL_MAT`, `BB_COR` or `CHIL`
  go non-positive, and the self-line derivatives avoid dividing by
  population-proportional denominators (`POPS(NUP)`,
  `U_STAR_RATIO`) that can vanish on a bad Lambda iterate.
  Identical results for valid iterates.
- `DO_NONMONO_KNOP` - Knop comoving-frame arbitrary-velocity
  pathway (`cmfgen/plane/fg_j_cmf_nonmono_v1.f`,
  `cmfgen/plane/mom_j_cmf_nonmono_v1.f`). Reference/diagnostic only for
  new configs (user directive 11-Nov-2026).
- `DO_JFNK` - optional Newton-Krylov population solver; admission
  gated on `USE_LAM_ES=T`.

### USE_LAM_ES under DO_NONMONO_SOB

- `USE_LAM_ES=T` (absent-key default since 10-Nov-2026, superseding
  the 21-Oct-2026 F default; the m15 recipe, job 12001805 /
  VAL_sob_jfnk): continuum J comes from the `FG_J_CMF_V13` formal
  solve on clamped `SIGMA` - the clamped Sobolev approximation of
  Dessart+2015 Sec 3.1. `SIGMA_TRUE` still drives the Sobolev line
  optical depths and their angular-factor floor. This is the
  intended SOB product: less accurate than `DO_NONMONO_KNOP` but a
  converged initial guess for it (`warm_restart_from`).
- `USE_LAM_ES=F` (opt-in experiment since 10-Nov-2026; was default
  21-Oct -> 10-Nov-2026, rewired 12-Oct-2026; banner declares it
  experimental and unsupported for production validation since
  25-Nov-2026): continuum J comes
  from `MOM_J_CMF_NONMONO_V1` (the Knop sign-split moment solver
  native to non-monotonic fields), paired with fixed-J +
  `GET_MOM_ALO_DIAG` + `GET_MOM_CELL_QUAR` in `do_var_cont.f` - the
  same recipe as KNOP+LAM_ES=F. `FG_J_CMF_V13` still runs each
  frequency, so `N_ON_J` factors come live from its ray moments via
  `GET_MOMS_NON_REL`; `SIGMA_TRUE` drives the Sobolev lines and
  the solver GAM term. JFNK stays inactive in this mode.
  `DO_KNOP_REFRESH` defaults on here as well as under KNOP
  (`MOM_KNOP_REFRESH_V1` converges the lagged J_NEXT/H_NEXT map
  within each sweep).
  Known caveat: the V13 ray moments feeding `K_ON_J` are built on
  clamped `SIGMA` while the solver's GAM term uses `SIGMA_TRUE`
  (real but second-order per the 13-Oct audit); the dominant
  advection-coefficient inflation comes from bound-free edge
  doublets on the continuum frequency grid (~2E-11 dLOG_NU).
  Runtime status FAIL - bounded non-contraction (job 12316397,
  `m15_2p82d_SOB_l14_lte_np_mom`, cancelled at it742,
  22-Oct-2026): containment holds (no luminosity runaway) but the
  lagged map never contracts; the `USE_LAM_ES=T` partner converged
  100% at it718. The earlier V11+V12 variant is not restorable by
  design - V11 assumes expanding flow and V12 linearizes V11. The
  GS freeze/maturity/watchdog stack stays KNOP-scoped
  (`GS_FREEZE` forced F outside `DO_NONMONO_KNOP`), so the moment
  path under SOB runs bare - `MAX_CHNG` and `DO_BA_RES_ACC` are
  the only rollbacks.

### Absent-key auto-sets under DO_NONMONO_SOB

`DO_NONMONO_FREQ_MATCH=T` under `USE_LAM_ES=F`; `DO_FG13_ALO=T`
under `USE_LAM_ES=T`; `DO_KNOP_REFRESH=F` under `USE_LAM_ES=T`
when not `DO_SOB_FG_KNOP` (25-Nov-2026: the read-site default
evaluated before the absent-key `USE_LAM_ES` flip, so an absent
key previously resolved armed-T on the V13 route; an explicit T
is still honoured but warns it is inert there);
`DO_BA_RES_ACC=T`; `DO_BA_RES_REUSE=T` (when
`DO_BA_RES_ACC` ends up T); `NG_ABSTAIN=T`; `DO_ES_ALI=T`
(11-Nov-2026, when `COHERENT_ES` and not `DO_SOB_FG_KNOP` - both
`USE_LAM_ES` routes carry a stored diagonal). `LAM_SET=F` is an
override auto-set (applied even to explicit T, like KNOP, because
the key is structurally incompatible, not a knob). Overrides of
explicit user values (e.g. `GLOBAL_LINE`) print a banner warning
to OUTGEN; auto-sets print ` SOB auto-set:` lines.

## 1.5 Velocity-law dispatch (VEL_LAW)

`VEL_LAW` is a stock key, but this fork modifies two of its
readers and adds two more, so the dispatch table lives here.
Dispatch in `cmfgen/new_main/cmfgen_sub.f` (VELTYPE):

- `VEL_LAW=11` - `SET_RV_HYDRO_MODEL_V3`: reads R/V/SIGMA verbatim
  from `RDINR` when `RD_IN_R_GRID=T`. `MODEL_SPEC ND` must equal
  the RDINR count exactly, so the forced non-mono table also
  forces `N_OB_INS=0`, `N_IB_INS=0`.
- `VEL_LAW=12` - `RV_SN_MODEL_SNIIN`: reads `input_hydro.dat`
  (7-column SNIIN table, rows ordered inner->outer, tau_es column
  last). Grid uniform in the file's own log(tau_es) plus
  `NBND_INS=3` packed at each boundary; linear-interp V, two-sided
  centred-difference SIGMA. Repaired 22-Nov-2026 (commit
  `10e2b47a`): the row read used `V_HYDRO` unsubscripted, so it
  demanded NHYDRO+6 values per row and EOF'd every real file at
  ~row 7; and `prepare_VADAT.py` now writes rows inner->outer.
- `VEL_LAW=13` - `RV_SN_MODEL_SNHYD(TAU_OPT=0)`: same V12 grid
  recipe (uniform log tau_es + boundary packing) but sourced from
  `SN_HYDRO_DATA` with monotone-cubic V and analytic-derivative
  SIGMA. Standing default for new runs (user directive
  22-Nov-2026, AGENTS.md 1.2).
- `VEL_LAW=14` - `RV_SN_MODEL_SNHYD(TAU_OPT=1)`: blended
  tau-Rosseland + log-R coordinate with structure-aware insertion
  at `sigma < -1` entries/minima and density spikes
  (`U_GAP_FAC` cell-gap floor). Same V/sigma evaluation as 13.
  The exact coordinate (rv_sn_model_snhyd.f:372-379):
  `U_i = sum_j [ max(ltau_j,TAU_LO) - max(ltau_{j-1},TAU_LO)
  + ALPHA_R*|dlogR_j| ]` with `TAU_LO = max(ltau_2,
  ltau_max - 4.0)` in natural-log units, `ALPHA_R=1.5`; tau is
  trapezoid of file-Kappa*rho (floored at 1e-6*kappa_max). SIGMA =
  analytic derivative of the Steffen cubic used for V
  (`dlnV/dlnR - 1`, :517); no |dV/dR| floor is applied - the
  Sobolev-side Dessart tau floor handles sigma ~ -1 downstream.
  Correctness evidence: docs/audits/VEL_LAW14_REVIEW_2026-11-27.md.

Grid-side comparison on the real m15 structure (ND=70, production
drivers, `unit_tests/work/grid_ab` + 27-Nov audit JSON): the
tau_es grids put ~80% of depths below V=100 km/s; VEL_LAW=14
trades most of that for a measured ~4.3x node concentration in
the CDS+ejecta zone (43 vs 10 nodes below R=36460; 13 vs 3 below
R=16805 - supersedes the earlier "5-8x" estimate), and captures
the CDS sigma spike deeper (-82 vs -63/-36, source -157).

`VEL_LAW` in `vadat_overrides` always wins: the forced non-mono
table stamps 11 only when the key is absent. `RD_IN_R_GRID=F`
with 12/13/14 (the readers build the grid themselves); `T` with
11. YAML `cmfgen.vel_law` feeds `--vel-law` in prepare_VADAT,
which stamps the same value in the initial VADAT.

## 2. Acceleration and ALO keys

- `DO_FG13_ALO` (default F; absent-key auto-set T under SOB +
  `USE_LAM_ES=T` since 21-Oct-2026 - under the LAM_ES=F moment
  path the matching diagonal comes from `GET_MOM_ALO_DIAG`
  instead): `FG_J_CMF_V13` stores the same diagonal OAB
  approximate-Lambda operator the Knop solver uses
  (`GET_FG13_ALO_DIAG`); `do_var_cont.f` fills the continuum
  dJ/dCHI, dJ/dETA diagonals under `DO_NONMONO_SOB`. J stays the
  `USE_LAM_ES=T` formal solution.
- `DO_ES_ALI` (default F; `ES_ALI_TOL` 1e-4, `ES_ALI_MAX`
  20-clamped): under `DO_NONMONO_SOB` + `USE_LAM_ES=T` +
  `COHERENT_ES` on the V13 path (not `DO_SOB_FG_KNOP`), the
  `RJ=TC` Lambda update in `comp_j_blank.f` becomes an
  Olson-Auer-Buchler Jacobi step on the stored V13 diagonal:
  `RJNEW = RJ + (TC-RJ)/max(1 - DJDETA*chi_scat, 0.01)`, and the
  J loop keeps passing until the last relative update is below
  `ES_ALI_TOL` or `ES_ALI_MAX`. Targets the
  ES-scattering-dominated m15 convergence mode (~2 Lambda
  steps/great iteration under the Eddington-only exit). Added
  27-Oct-2026; verified in `adp2_esali` (J loop converges to 1e-4
  in ~2 passes, 0 unconverged freqs). Warm-arm evidence
  11-Nov-2026 (`docs/audits/SOB_ACCELERATION_RESULTS.md`): X2
  (`DO_ES_ALI=T` alone, warm restart from the converged FRESH_sob
  checkpoint) converged in ~30 new iterations where the
  unaccelerated X1 twin was still pinned at maxapp ~0.5 after
  ~80; X3 adding `DO_ES_ALI_BA` matched X2, so the BA factor stays
  opt-in. The `USE_LAM_ES=T` + `DO_NONMONO_KNOP` branch was inert
  until 02-Nov-2026: the update always read the V13 diagonal via
  `GET_FG13_ALO_DIAG`, but the KNOP formal solver is
  `FG_J_CMF_NONMONO_V1`, which stores its diagonal in
  `FG_J_CMF_NONMONO_MOD_V1`, so every frequency fell back to plain
  Lambda (`KNOP_fresh`: 56377/56377 `ES_ALI` fallbacks each
  iteration). The accessor is now selected on the same predicate
  as the solver dispatch (`GET_FG_ALO_DIAG` under KNOP /
  `DO_SOB_FG_KNOP`); `KN_NO_GAIN` is forced T under KNOP since
  08-Nov-2026, so the DEN floor (0.01) albedo <= 1 assumption
  holds there automatically. Extended 29-Oct-2026
  (`ES_ALI_MOM_ON`): same OAB Jacobi update on the `USE_LAM_ES=F`
  branch under `DO_NONMONO_SOB`, where `MOM_J_CMF_NONMONO_V1`
  solves J with the scattering source lagged
  (TA = ETA + chi_scat*RJ_prev). Uses `GET_MOM_ALO_DIAG`. Gated on
  `DO_ES_ALI` + `COHERENT_ES` + `.NOT.USE_LAM_ES` +
  `.NOT.DOB_FG_KNOP`; shares ES_DJ/ES_ALI_TOL/ES_ALI_MAX and the
  fallback/print accounting. Runtime status unverified.
- `DO_ES_ALI_BA` (default F; `ES_ALI_BA_MIN` 0.1): same Jacobi
  factor as DO_ES_ALI applied to the stored ALO when
  `do_var_cont.f` fills TX - the BA linearization then sees the
  ES-converged J response `dJ ~ ALO/(1-DJDETA*chi_scat)` instead
  of the bare-Lambda one. Without it the BA under-corrects a
  scattering-coupled mode every iteration (fixed-ratio ~0.97/it
  decay on the m15 2.82d model). Applied per depth-row to every
  filled band; composes with `DO_FG13_ALO_TRI`. Added
  28-Oct-2026; requires `DO_ES_ALI=T` since 11-Nov-2026
  (CHECK_PARAM_CONSISTENCY hard-stops the pairing; the X3 warm
  arm showed no gain over X2, so it stays opt-in). The J-loop DEN
  floor (0.01) and `ES_ALI_BA_MIN` (0.1) differ by design;
  activations are counted in `ES_ALI_NFLOOR`/`ES_ALI_BA_NFLOOR`.
- `DO_FG13_ALO_TRI` (default F): `FG_J_CMF_V13` stores the
  off-diagonal Lambda* entries `dJ(K)/dS'(K+/-1)` from the direct
  A1/A2 (+ propagated A0*A2) ray coefficients
  (`GET_FG13_ALO_TRI`); `do_var_cont.f` fills
  `TX(K,K+/-1,1:2)` so the banded BA matrix couples adjacent
  depths (Olson & Kunasz 1987). Effective only with
  `SOL_METH=TRI`/`NUM_BNDS=3` (pipeline:
  `depth_independent_solver: false`); with NUM_BNDS<3 the fill is
  skipped. Ray-node neighbors map to depth neighbors - exact on
  the non-inserted grid, approximate on inserted/tangent nodes;
  source-gradient (dS) couplings ignored, matching the diagonal.
  Added 28-Oct-2026; runtime status unverified.
- `DO_BA_RES_ACC` / `DO_BA_RES_REUSE` - BA residual-acceleration
  pair; auto-set under SOB.
- `NG_ABSTAIN` - per-component NG-extrapolation hold on
  zero-diff and non-positive extrapolation; auto-set under SOB
  (`test_ng_abstain.py`).
- `NG_SHADOW_DIAG` (default F, diagnostic-only) - at each NG
  event, `DO_NG_BAND_ACCEL_V4` reads up to 16 SCRTEMP records
  and reports per band and window size N=4..nwin the
  affine-fit coefficient L1 norm, implied prediction move,
  predicted-residual ratio, and the successive-prediction
  settle ratio that gates adaptive Ng (Magritte, Ceulemans
  et al. 2024, Sec 4.2; bounded investigation NG1). Lines are
  prefixed `ng shadow`. No physics or state change; needs
  `DO_NG=T` to have events (`test_ng_shadow.py`).

## 3. BA-solver guards and quarantine

- `BA_LQ_DEM` (default 1E12) / `BA_LQ_PAT` (default 8) -
  LAMBDA-phase quarantine in `solveba_v13.f`: a row demanding
  beyond BA_LQ_DEM% for BA_LQ_PAT consecutive LAMBDA calls parks
  into `BA_QUAR_MASK`; its correction is suppressed on every call
  so it leaves MAXCH and stops churning. Armed only with
  `DO_BA_RES_ACC`+`DO_BA_QUAR`; bounded by `BA_QUAR_MAX`; each
  parking prints its EQ_DESC and the suppressed max demand prints
  per call. Floor-level rows (`BA_DEM_PFLOOR`) and held calls are
  not counted. Reverses the deliberate "parked rows still update
  on LAMBDA calls" choice for LAMBDA-parked rows only: under the
  precond gate every call is LAMBDA, and these rows demonstrably
  do not contract. (21-Oct-2026)
- `BA_DEM_PFLOOR` - demand floor used by the quarantine census.
- `BA_EX_RELFRAC` - exit-only relative floor replacing the removed
  `BA_DEM_RELFRAC`.
- `BA_DEM_RELFRAC` - REMOVED 17-Oct-2026 (commit b398571a): was a
  relative demand floor in `SOLVEBA_V13` excluding rows below
  `RELFRAC*max(POPS)` at a depth from MAXCH/MAXCH_APP. It hid
  demand at populated depths. Do not re-add.
- `fiddle_pop_corrections_v2.f` MAJOR scale (21-Oct-2026): under
  DO_NONMONO_SOB or DO_NONMONO_KNOP each row's contribution to the
  per-depth SCALE is capped at the bound the per-element apply
  clamp already enforces (DPTH_BIG_LIM/DPTH_LIT_LIM). One
  unsatisfiable row used to veto the whole depth's corrections
  (SCALE ~1E-69); the bad row still takes its identical clamped
  step, the rest of the depth moves. Stock path unchanged.
- `solve_for_pops.f`: `N_PRECOND_TOT` counts only calls that fail
  to improve on the stretch's best demand by 10% (`GS_PC_BEST`
  tracker), matching the documented "no-progress" semantics - a
  still-contracting cold start was previously stopped at
  `GS_PRECOND_LIM` regardless of progress.
- `solveba_v13.f` (05-Nov-2026): after `CMF_BLKBAND_V3` returns,
  `STEQ` is scanned for non-finite entries; the first bad (level,
  depth, value, count) is printed and the run stops. Previously
  NaN corrections sailed through every `.GT.`/`.LT.` scale guard
  into `POPS` and `SCRTEMP`. Unconditional - a non-finite solve
  output is always a defect.
- Empty-level net-bracket guards (29-Oct-2026, not flag-gated):
  `ZNET_SIM`, `ZNET`, `BETAC` all divided by an emissivity that is
  exactly 0 where the upper level is empty (thin outer depths on
  checkpoint restarts). `RJ*CHIL/0=-Inf` then `ETAL*ZNET=0*Inf=NaN`
  poisoned `STEQ_T` -> `DEP_RAD_EQ` -> `RLUMST(ND)`, tripping the
  SOLVE_FOR_POPS non-finite stop and killing every warmS/l14_freet
  requeue at iteration 864 (trace: SIM=16, K=2, ML=1204). Now the
  `ETAL=0` resonance contributes no bracket term (ZNET_SIM/ZNET
  consumers all multiply by ~0 anyway); `BETAC=0` likewise.
  Verified: warmS it864/865 complete (~315-470 s, STEQ_T finite,
  max 8.5E12 at K=70).

## 4. KNOP lagged-field and stability keys

- KNOP forced keys (08-Nov-2026): under `DO_NONMONO_KNOP=T`,
  `RD_CONTROL_VARIABLES` forces `KN_NO_GAIN=T`, `USE_LAM_ES=T`,
  `DO_NONMONO_FREQ_MATCH=T` even over explicit F (` KNOP
  override:` banner; absent keys print ` KNOP auto-set:`).
  - `KN_NO_GAIN=T`: with F the stock `NEG_OPAC` floors can leave
    `CHI < CHI_SCAT`, so the source coefficient sum exceeds 1 and
    the cell amplifies photons; the lagged `I_P/I_M_DEPTH_NEXT`
    store compounds that gain each great iteration (~8.6x/iter
    measured; `KNOP_fresh` 02-Nov-2026 diverged to `L_out` ~1e29).
    The floor makes the lagged map non-expanding in sup norm.
  - `USE_LAM_ES=T`: the `MOM_J_CMF_NONMONO_V1` moment path never
    converged under KNOP - its lagged `J_NEXT` map is expansive
    there (m15_2p82d job 12747893, 06-Nov-2026: adaptive
    `JNEXT_RELAX` contraction factor 1.24, every refresh
    non-converged, NaN in SOLVEBA at it1478; TypeIIn same
    verdict). The MOM code stays in the tree and stays reachable
    under `DO_NONMONO_SOB`; only the KNOP route is closed.
  - `DO_NONMONO_FREQ_MATCH=T`: with F every compressing upwind
    read uses the adjacent bin; on m15 the per-cell Doppler drift
    reaches ~330 bins, so the read is at the wrong frequency in
    the oscillatory-eigenvalue regime (Ceulemans et al. 2025,
    JQSRT, arXiv:2505.08309; Knop et al. 2009 II, A&A 501, 813).
    Correctness requirement backed by theory - no failed run is
    yet attributed to F.
- WP-1 default flips and tolerance-ladder invariants
  (26-Nov-2026, `docs/audits/KNOP_REMEDIATION_WP1_2026-11-26.md`):
  - `KN_RAY_ORDERED` (default F since 26-Nov-2026, was T):
    ordered characteristic march in the FG refresh (Knop Paper
    II ordering). The WP-0 fixture measured the claimed
    freshness broken on inserted ray points: at WINX>0 the
    compressing donor share is a full pass stale and the march
    plateaus at the pass cap. F selects the alternating
    `SWEEP_FLD` passes, which update the whole donor column at
    once. T re-arms the ordered march for A/B tests.
  - `KN_RAY_BULK_TOL` (default 0 = off on every route; KNOP
    default lowered from 6.0E-2): residual floor on the bulk
    exit (IXIT=2). The measured floor is bounded
    non-contraction of the lagged map, not convergence (CHECKLIST
    nbulk arm). 0 disables the term so every production ray must
    reach `KN_RAY_TOL` to certify.
  - `GS_MATURE_TOL` (default 5.0E-2, read-only): VADAT exposure
    of the population maturity gate `GS_MATURE_TOL` (a
    `CONTROL_VARIABLE_MOD` PARAMETER since 26-Nov-2026; before
    that a `SOLVE_FOR_POPS` local). Read-only:
    `CHECK_PARAM_CONSISTENCY` stops any non-default value.
  - Hard invariants in `CHECK_PARAM_CONSISTENCY` (STOP, not
    warning), active on the FG-census routes
    (`DO_NONMONO_KNOP`, `DO_SOB_FG_KNOP`):
    - I1: `KN_RAY_BULK_TOL < GS_MATURE_TOL` -- a ray may not
      certify looser than the gate.
    - I2: `KN_RAY_TOL <= 0.1 * GS_MATURE_TOL` -- the certified
      field must be >=10x tighter than the gate.
    - I3: an explicitly-set `FG_COMMIT_NC` must not equal
      `GS_MATURE_TOL` -- a capped ray pins `GS_MAX_ALL` exactly
      at the gate (freeze/release limit cycle). The absent-key
      legacy default (0.05 == 0.05) is exempt until WP-4
      replaces the fixed cap.
    - I4 (census): `KN_NACT`/`KN_NFULL` count only rays with a
      transfer domain (`NI_RAY >= 2`); `FG_RF_ALL_CONV`
      additionally requires `KN_NINV = 0` so an early return
      left `KN_CS_INV` still fails the all-ray verdict closed.
- WP-2 freshness contract and order modes (26-Nov-2026):
  - `KN_RAY_ORDER_MODE` (default -1 = follow `KN_RAY_ORDERED`):
    march-order selector for the refresh solve. 0 = ordered
    `SWEEP_FLD_ORD`, 1 = frequency-major alternating
    `SWEEP_FLD`, 2 = frequency-major with the first pass marched
    down the ray's dominant coupling sign -- a uniform-
    compressing ray (some GAM<0, none >0) starts descending so
    its L+NFM donor column is fresh on pass 1; every other ray
    keeps the mode-1 ascending start, which reduces mode 2 to
    mode 1 exactly on expanding-only rays. 3 is reserved for
    the SCC/topological ordering and stops until mode 2 is
    measured on mixed rays. An explicit value wins over
    `KN_RAY_ORDERED` with a warning; outside 0..3 stops. The
    WP-0 fixture records the contract: `ord` reads ~WINX=0.943
    of its compressing donor weight stale on I_M and ~0 on
    I_P; `uno` alternates 1.0/0.0; `m2u` opens fresh and
    converges in 2 passes; `m2m` is bitwise identical to
    `unom` (test_knop_ordering, 50/0).
  - `KN_FRESH_DIAG` (default F, diagnostic-only): per-pass
    census of the |Q| donor weight read from a not-yet-updated
    cell -- the stale share the ordered march's header comment
    claimed was zero. Stored per (pass, dir, ray) in
    `KN_FRESH_FR` (dir 1=I_M, 2=I_P, -1 = pass never ran) and
    printed on LUER; the solve is unchanged.
- WP-3 committed-state gate (26-Nov-2026):
  - `KN_GATE_COMMIT_RES` (default T): after each ray's commit
    decision writes a store, run one L_QP verify pass on the
    committed point -- the production field at unstored cells,
    the committed store at stored cells -- and fold its max-norm
    residual (`FG_RF_COMMIT_RES`) plus the committed bad-cell
    census (per-ray `FG_NBAD_CM_RAY`/`FG_NTOT_CM_RAY`, summed to
    `FG_CM_NBAD`/`FG_CM_NTOT`) into the maturity exports. The
    `GS_MAX_ALL`/`GS_NBAD_EFF` feed in `solve_for_pops` reads
    the committed metrics when the census exists, falling back
    to the pre-commit residual/best-pass census only when the
    verify never ran (flag off, or no store written that
    refresh). Motivation: the pre-commit metrics measure the
    solve's own iterate, not what `WRITE_STORE`/`WRITE_STORE_CAP`
    actually kept -- a capped provisional commit can park the
    store far off the fixed point while the best-pass census
    reads clean. The m15-arm fixture records the contract:
    `nccm` holds at CRES=1.4e5 with 55552 committed bad cells
    on an NC-bounded write where the legacy census reported 0;
    `capok` certifies at ~4e-16; `nclg` disarms
    (test_knop_commit_gate, 10/0).
- WP-4 store line search (26-Nov-2026):
  - `KN_STORE_LS` (default T): on a non-certified commit the
    store step backtracks on the committed-state residual
    (WP-3 verify) instead of taking the fixed `FG_COMMIT_NC`
    bound. The search evaluates the committed residual at the
    retained store, then at alpha = 1, 0.5, 0.25, 0.125 blends
    (stored <- stored + alpha*(fresh-stored), both march
    directions), and accepts the first rung that strictly
    descends; `FG_COMMIT_NC` survives only as the rung floor.
    If no rung descends the previous store is retained and a
    rejection is counted (`FG_LS_REJ`, per-ray `FG_RAY_LSA=-2`);
    accepted alphas histogram into `FG_LS_HIST(5)` (rung 5 =
    the NC floor when 0 < `FG_COMMIT_NC` < 0.125) with
    `FG_LS_R0`/`FG_LS_R1` summing residual before/after --
    all printed in the refresh epilogue census. Certified
    commits, the `FG_COMMIT_JUMP` veto, and the non-finite /
    cap guards are unchanged; F restores the fixed NC cap
    exactly (legacy `WRITE_STORE_CAP` provisional write, or
    retain under NC=0). Motivation: the fixed NC bound
    ratchets O(100) refreshes toward the fixed point where a
    residual-tested full step commits in one (Eisenstat &
    Walker inexact-Newton backtracking; Arramy et al. 2024,
    A&A 690, A12, Sec. 2.2.2). The m15-arm fixture records
    both directions: `lsok` descends CRES 1.4e5 -> 3.95e-2 in
    one refresh (alpha=1); `lsrj` rejects an over-relaxed
    expansive step and keeps the prior store
    (test_knop_commit_gate, 18/0).
- WP-5 signed-Aitken store relaxation (26-Nov-2026):
  - `KN_OMEGA_SGN` (default F): replaces the FG_ADAPT_OMEGA_V2
    pool estimate inside the KN_ADAPT_RELAX adaptive branch
    with a signed vector-Aitken omega priced on the OUTER
    store-map residual sequence -- the per-slot
    (fresh - snapshot) defect each sweep write commits
    against, `x <- snap + omega*r`. omega_{k+1} =
    -<dx,dr>/<dr,dr> (Ramiere & Helfer 2015, JCP 283, eq. 36;
    scalar map gives the fixed-point weight 1/(1-lambda) in
    one step). Guards: denominator collapse holds omega, a
    non-decreasing residual halves once then falls back to
    the unaccelerated estimator, the raw price clamps to
    [KN_OSGN_MIN, KN_OSGN_MAX], and a CHI/ETA_EXT fingerprint
    move (populations changed) clears the slot history.
    T2=0 slots are skipped -- the unrelaxed T1 write there
    never saw omega. The V2 estimator keys the outer step to
    INNER per-ray pass ratios, a category error that cannot
    express an over-relaxed step; F keeps it as the default
    while the armed arm is measured. `test_knop_osgn` 14/0.
  - `KN_OSGN_MIN` (default 5.0D-2), `KN_OSGN_MAX` (default
    2.0D0): clamp range for the signed-Aitken omega.
  - `KN_FORCING_TOL` (default F): Eisenstat & Walker 1996
    Choice-2 forcing term on the refresh's inner stopping
    tolerance (WP-6). The refresh is the inner solve; the outer
    residual is the population demand MAXCH handed in before
    each refresh call. eta_k = GAM*(F_k/F_{k-1})^ALP, floored at
    GAM*eta_{k-1}^ALP when that exceeds 0.1, capped at 0.9; the
    effective tolerance KN_RAY_TOL_EFF =
    clamp(eta*TMAX/0.9, TMIN, min(TMAX, 0.1*GS_MATURE_TOL)) --
    the I2 bound means the field tolerance can never loosen past
    a tenth of the population gate. KN_RAY_TOL_EFF replaces
    KN_RAY_TOL inside the refresh machinery (certification,
    sweep exits, MARCH census, hot-run thresholds); F keeps the
    fixed KN_RAY_TOL everywhere. eta, the raw ratio, the
    safeguard flag and the effective tolerance land in the
    ITER_LEDGER frc_* columns. `test_knop_frc` 35/0.
  - `KN_FRC_GAM` (default 0.9), `KN_FRC_ALP` (default 2.0):
    the E-W forcing coefficient and exponent.
  - `KN_FRC_TMIN` (default 1.0D-4), `KN_FRC_TMAX` (default
    1.0D-2): clamp band on the priced tolerance; the I2 cap
    min(TMAX, 0.1*GS_MATURE_TOL) still applies.
  - `BA_SYM_PIN` (default F): WP-7c symmetric elimination of
    floor-pinned level unknowns in the CMF_BLKBAND_V3 DIAG block
    solve. Each level row J with POPS(J,K) <= MAX(1E-30,
    BA_DEM_PFLOOR) is eliminated before DGEEQU: the column's
    contribution is carried into the RHS, row and column are
    zeroed, the diagonal is set to 1 and the pinned correction
    solves as exactly zero -- the same result the post-solve
    POP_MASK_FLOOR mask forces, but without leaving a near-null
    direction inside the factorisation (the earlier row-only pin
    crashed the coupled block, job 12868544). ED/T and the
    conservation rows are never pinned. DIAG path only; F is the
    default path unchanged. `test_cmf_blkband_pin` manufactured
    near-null block 13/0.
  - `BA_NC_MAXPOP` (default F): WP-7c Arramy et al. 2024 (Sec.
    2.1.1) conservation-row placement in GENERATE_FULL_MATRIX_V3.
    Each species' number-conservation equation moves from the
    fixed EQ_SPECIES slot to the row of its most populated level
    at that depth; the displaced rate equation relocates to
    EQ_SPECIES so the equation count is unchanged. Ground-state
    (EQXzV) rows are excluded from candidates because the
    ionization-conservation replacement may overwrite them. F is
    the default path unchanged. `test_gfm` ncmax/ncdef 14/0.
  - `BA_REL_TRUST` (default F): WP-9 release trust seed in
    SOLVE_FOR_POPS. On the first BA call after a freeze release
    the resacc trust radius `BA_TRUST_CUR` is re-seeded from the
    committed field residual (`GS_COMMIT_RES`, the WP-3
    `KN_GATE_COMMIT_RES` measurement) as
    `MAX(1.0, MIN(BA_TRUST, GS_COMMIT_RES))` instead of riding
    its pre-freeze ratchet value, so the release-controller
    overshoot -- a ~1E7 demand committing at the MAX_LIN_COR cap
    and re-poisoning the released field -- is quarantined at a
    field-consistent radius. The step still passes the
    DO_BA_RES_ACC residual-acceptance test and the all-ray-FULL
    release gate is untouched. F is the default path unchanged.
- `GS_FREEZE` and the freeze/maturity/watchdog stack are
  KNOP-scoped: forced F outside `DO_NONMONO_KNOP`.
- `DO_KNOP_REFRESH` - converges the lagged J_NEXT/H_NEXT map
  within each sweep via `MOM_KNOP_REFRESH_V1`; defaults on under
  KNOP and under SOB + `USE_LAM_ES=F`. Under plain SOB +
  `USE_LAM_ES=T` an absent key resolves F (25-Nov-2026) and an
  explicit T is armed-but-inert: the production call site is
  gated on `.NOT. USE_LAM_ES` and `CHECK_PARAM_CONSISTENCY`
  prints an inert warning.
- `JNEXT_RELAX` - adaptive lagged-J relaxation factor.
- `GS_FOREIGN_SEED` (default F): on restart, accept a
  fingerprint-valid but tag-mismatched `FGINEXT_GS` (e.g. from a
  converged SOB run of the same model) as the lagged-field seed
  instead of a cold `I_NEXT=0` rebuild. Seed stays immature
  (`I_NEXT_GS_OK=F`) so `N_LAM_RWARM` re-converges the map first.
- `FG_SEED_COMMIT` (default T, 02-Nov-2026) - per-ray seed-fill
  waiver: a ray that has never committed via WRITE_STORE skips
  the `FG_COMMIT_JUMP` veto once, so a cold/foreign store can take
  its converged first fill instead of being rejected ray-by-ray
  (esB: 67/85 rays jump-rejected at fill gaps 1E2-6E2). Validated
  FGINEXT_GS restores count all rays seeded.
- `SNP_P`/`SNP_M` per-frequency snapshots (02-Nov-2026): the
  in-sweep store write relaxes against the frozen pre-frequency
  column on EVERY solver call, so under DO_ES_ALI (~12
  comp_j_blank re-entries/freq) the store lands on the last
  (ES-converged) pass's field rather than the first; churn
  metrics keep their NEW_FREQ one-sample census. Both changes
  touch `fg_j_cmf_nonmono_v1.f` only; default path unchanged
  (KNOP-scoped).
- `LUM_DIV_FAC3`/`LUM_DIV_N3` (KNOP only; default N3=0 = off):
  third, drift tier of the `LUM_DIV_*` runaway guard in
  SOLVE_FOR_POPS. Trips `LUM_BLOCKED` on N3 consecutive calls each
  with `L_out` growth above factor FAC3 (default 1.005). The
  FAC=1E2/N=3 and FAC2=5/N2 tiers both require per-call growth
  >= 5x, so a sustained ~2%/call drift is invisible to them - the
  m15_2p82d_KNOP_l14_freet restart climbed 1.6E9 -> 1.2E10 L_sun
  over ~300 iterations under `hold_pops` and only stopped on
  non-finite `L_out`. Opt-in: a legitimate field re-equilibration
  after a cold rebuild produces the same slow climb signature, so
  enabling is a per-run policy. Set `LUM_DIV_CAP` to a
  model-physical level (m15 ~3E9) alongside. Added 29-Oct-2026.
- `MOM_RF_DIV_BULK` - tolerated RF-divergence excursions bounded:
  max streak 3 consecutive and excursion magnitude capped at
  1E6 x MIN_RES (the 1E3 abort bound stands); committed-field
  finite/>1E40 and 50x-jump guards already existed.
- `MOM_XM_POS` / `MOM_XMWORST` - negative-RHS clipping guard plus
  a report of the largest clipped value each sweep so the
  safeguard's cost is measured, not assumed.
- `MOM_NUSOLV` (default F, 25-Nov-2026): after each MOM refresh
  pass, re-solve the captured row equations with the upwind-J
  term implicit along each depth's frequency chain
  (`MOM_NU_CHAIN_SOLVE`). At fixed depth every row owns one
  upwind edge -- bluer L-1 expanding, matched redder
  L+MOM_NM_J compressing -- so the coupling is a functional
  graph: chains plus closed rings at the sign flips (the
  expanding/compressing 2-cycles that pin the lagged map at
  rho~1; matched strides >1 make longer rings). Solved exactly
  in O(NCF) per depth -- the direct-solve treatment the CMF
  literature uses for this boundary-value structure in
  wavelength (Baron & Hauschildt 2004, A&A 427, 987). Captured
  coefficients are the swept row's own post-LAG_SCL/EPS_SCL, so
  the composite map's fixed point is the swept system's. The
  lagged H and EPS scatter terms stay explicit. Census: rings,
  singular-frozen rings, frozen rows, rejected negatives print
  per pass on the `nu-solve` line; companion `nu-census` lines report
  coupling ratios, source inflation, the H-share of `NU_X0`, and
  nonpositive-excess rows. `NU_FRESH` freezes
  active-set-skipped columns. Targets the bounded
  non-contraction measured on the m15 compressing flow
  (lamesF arm, job 13646071).
- `MOM_BLKGMRES` (default F, 24-Nov-2026): after each MOM refresh
  pass, solve the captured wavelength-depth rows in one matrix-free,
  restarted GMRES solve (`MOM_BLK_SOLVE` and `MOM_BLK_GMRES_MOD`).
  Depth-tridiagonal terms and the sign-keyed upwind wavelength term
  are implicit together. `NU_X0` already excludes the upwind term
  read by the sweep, so that coupling becomes implicit exactly once;
  the diagonal PSI absorption remains only in `NU_TB`. The left
  preconditioner is the existing per-frequency depth-tridiagonal
  inverse applied by DGTSV. The write-back keeps the chain solve's
  damped, non-negative, trust-region policy. Inactive frequency
  columns are pinned to their incoming iterate. When both
  `MOM_BLKGMRES` and `MOM_NUSOLV` are true, the block solve takes
  precedence. Census prints iterations, restarts, the final true
  preconditioned residual, rejected writes and write caps on the
  `blk-gmres` line; companion `nu-census` lines decompose the captured
  rows. Lagged H,
  EPS scatter, source functions and boundary closures remain outer
  iteration data. Diagnostic solver key.
- `MOM_AUG_RES` (default F, 25-Nov-2026): augmented-residual
  diagnostic. Captures the H-row and lagged-term coefficients
  (`NU_HU/NU_HL/NU_HS`, `NU_VB/NU_VC`, `NU_EPS`, `NU_EPSU`,
  `NU_EPSP`, `NU_EPSPU`, `NU_LAGSCL`, `NU_SRC`, `NU_CUPU`,
  `NU_T1B`, `NU_KUH`, `NU_OUTBC`, `NU_HCHK`) alongside the `NU_*`
  row capture, and at each `MOM_KNOP_REFRESH_V1` exit calls
  `MOM_AUG_RESID_V1` to evaluate the full coupled J+H equations on
  the committed `J_NEXT_FRQ`/`H_NEXT_FRQ` fields - after block
  GMRES, negative-cell rejection, write caps, damping, best-pass
  selection and outer Anderson acceleration. Two evaluations are
  reported on the `mom-aug` lines: RESID_IMPLEMENTED (the rows as
  assembled, stabilisers active) and RESID_TARGET (the same rows
  with the equation-changing stabilisers removed: centred EPS
  midpoint stencil, `LAG_SCL=1`, uncapped EPS coefficients, no
  `XM_POS` clip, no `MOM_H_IMPL` diagonal mass). Moment-row and
  H-row norms are reported separately, raw and row-scale-normalized,
  with worst cells, stabiliser-activation counts, off-grid reads and
  non-finite rows. Read-only: no solver state changes. Also forces
  the `NU_*` capture on when neither `MOM_NUSOLV` nor
  `MOM_BLKGMRES` is selected. The `blk-gmres` report additionally
  gains the raw `||b-Ax||/||b||` residual and the committed-field
  (post-write-policy) raw and preconditioned residuals.
  Diagnostic key.
- Dead-key warning (25-Nov-2026): `MOM_BLKGMRES`, `MOM_NUSOLV`,
  `MOM_H_IMPL` and `MOM_AUG_RES` only do work inside
  `MOM_J_CMF_NONMONO_V1` (SOB +
  `USE_LAM_ES=F`). Set T on any other route - production
  `USE_LAM_ES=T`, KNOP, `DO_SOB_FG_KNOP`, or stock - and
  `CHECK_PARAM_CONSISTENCY` prints an inert warning; the run is
  legal but the keys are dead in OUTGEN.
- `KN_NI1_FULL` (default T since 26-Nov-2026, was F): certify a
  production NI=1 ray KN_CS_FULL (reason 0) in
  `KNOP_RAY_SOLVE_V1` instead of leaving it KN_CS_INV. A
  single-point ray has no marchable transfer domain, so vacuous
  FULL is the honest class; with F the ray holds
  `FG_RF_ALL_CONV` FALSE forever and, through the
  rf_not_allconv fold-in, the LAMBDA handoff shut -- the
  permanent ninth-conjunct deadlock of defect (d) in
  `docs/audits/KNOP_REMEDIATION_WP1_2026-11-26.md`. F remains
  for A/B regression checks.
- `KN_LAM_PARTIAL_HANDOFF` (default F, 24-Nov-2026): split the
  field verdicts. `LAM_FIELD_NOT_CONVERGED` (strict, read by the
  converged exit in `cmfgen_sub` and the `CONV_EXIT_EVAL`
  replica) always carries the gate result plus the
  rf_not_allconv fold-in; `LAM_FIELD_NOT_READY` (weak, read by
  `LAM_STABLE_STEP`, `LAM_HANDOFF_OK` and the stable counter)
  carries the fold only when this key is F. With T a field
  mature by the share gates but not yet all-ray FULL can hand
  off while the exit keeps waiting for every ray.

## 5. OBS/render keys

- `DO_SOB_DEP_WINGS` (OBS CMF_FLUX_PARAM key, default F;
  `SOB_DEP_WING_EXT` 2000 km/s): in `ADD_SOB_LINE_PROF`, continue
  the deposited Gaussian core as 1/x^2 damping wings out to
  SOB_DEP_WING_EXT instead of truncating at 3 Doppler widths. The
  full CMF solve carries Voigt/damping wings, so the pure-Gaussian
  deposit renders optical line wings systematically too narrow
  (m15 2.82d SOB-render comparison vs the Dessart .fl: HeI 5876
  FWHM 18A vs blank-render ~30A; with wings on, 32A). Opt-in via
  `cmfgen.cmf_flux_param_overrides`. Wing shape uses the true
  per-line Voigt (`VOIGT` with
  A_VOIGT=(ARAD+1.55E4*|C4|^0.667*Ne)/(4*pi*dnu_D), same as
  SET_PROF_V6) since 28-Oct-2026.
- `DO_SOB_DEP_FLOOR` (OBS CMF_FLUX_PARAM key, default F;
  `SOB_DEP_CHI_FLOOR` 1.0): masing/negative-opacity cells in
  ADD_SOB_LINE_PROF keep their emissivity with opacity floored
  (optically-thin emitters) instead of being dropped. The dropped
  cells sit in the dense inner shell and carry the P-Cygni wing
  flux; with the floor on, injected SOB emissivity is 3.35E17 vs
  ~2/3 of that dropped before (m15 2.82d esali OBS). Added
  28-Oct-2026.
- `DO_KNOP_SPEC_GATE` (default T): `cmfgen/obs/knop_spec_gate.f`
  classifies GS_STATE + stall marker + manifest before the
  OBSFLUX/OBSFRAME writes; non-converged KNOP models stop under T
  or render a marked diagnostic under the F opt-out
  (`test_knop_spec_gate.py`, 9 cases).
- `OBS_FG_KNOP` (OBS CMF_FLUX_PARAM key, default F): diagnostic
  formal-solver selector for the observer frame. When T,
  `COMP_JCONT_V4.INC:925` dispatches the per-frequency solve to
  `FG_J_CMF_NONMONO_V1` on the OBS-built `SIGMA_TRUE`, and the
  three solver-machinery hooks arm (OBS grid registration
  `SET_FG_GS_OBS_GRID`, refresh, maturity census). Every other
  gate keeps its production `DO_NONMONO_SOB` behaviour, so an
  OBS_FG_KNOP render compares formal solvers on an otherwise
  identical state. `INC_GRID=T` + OBS_FG_KNOP aborts loudly
  (same guard as DO_NONMONO_KNOP). Opt-in via
  `cmfgen.cmf_flux_param_overrides`; not production-ready -
  audit `CMF_FLUX_V12_VS_KNOP_FORMAL_SOLVER.md` rates its
  GS-coupled field an immature iterate (direction of the
  boundary-value correction robust, magnitude pending a
  convergence program).
- `DAMP_FJ`/`J_RESCUE`/`RESC_JMAX` - continuum-J rescue controls
  (`COMP_JCONT_V4.INC` etc.); `J_RESCUE=T` restores the previous
  frequency's field when stored J exceeds
  `MIN(RESC_JMAX, 1e4*JNU_PREVmax)`. `RESC_JMAX` default `1e2`
  (22-Nov-2026; was `1e-3`) is a runaway sanity rail - the
  relative bound is the working detector. The old `1e-3` rail
  clipped legitimate m15 fields (~1e-1) and collapsed the
  spectrum ~100x. `DAMP_FJ` stays `1.0` (off).
- `ES_RESID_DIAG` - OBS-side ES residual writer.

## 6. Diagnostic-only keys

All default F; instrumentation only - they change output, not the
solve.

- `GS_DRIFT_DIAG` - per-iteration drift reporting.
- `KN_BADID` - cross-refresh bad-cell identity census (round-2
  phase-2): each ray's last-pass above-`KN_RAY_TOL` cells are
  recorded as packed (ray,freq,point,channel) ids on sparse
  per-ray lists; the serial epilogue merges and diffs consecutive
  censuses (intersection, newly bad, recovered, Jaccard, per-
  channel split). Diagnostic only; criterion and count unchanged.
- `KN_LUM_TRACE` - first-detection scan of `STEQ_T` after each
  energy-equation accumulation site (site, freq, depth, operands).
- `KN_WATCH_CELL`/ray lifecycle prints.
- `KN_NOCP_PROBE`, `KN_NOP_*` - pass telemetry.
- `KN_PSWEEP_*` - complete-map capture incl. binary X/Y
  snapshots (`kn_psweep_source_test.py`, `kn_psweep_check.py`).
- `KN_FREEZE_COEF`, `KN_MAT_FREEZE` - coefficient-freeze arms
  with arm-event fingerprinting.
- `SOB_RH_TEST_S1` (05-Nov-2026): accumulates the sob_rh_couple
  partner sum with unit source function; remote `dJ(S=1)` above
  ~1 means the two-point kernel overcounts. The `DRATE`
  finiteness stop now reports DJ, CHIL, BB_COR and partner count;
  the dominance census uses `RH_MAXDJ` normalisation (cells where
  remote `dJ` exceeds the local mean intensity) instead of the
  degenerate `|DRC|>|ETAL*ZNET|` test.
- `JFNK_HSCAN` - fences the production OPPROBE call site.
- `ITER_HISTORY`/`ITER_LEDGER` - machine-readable iteration rows
  (`test_iter_ledger.py`).
- `OUT_LEVEL` (integer, default 1, clamped 0..3) - output verbosity;
  `CONV_DIAG_LVL` is the legacy alias used when `OUT_LEVEL` is absent.
  Level 0 keeps only iteration framing/parser lines/ISSUES, 1 the
  summary tables, 2 the former full volume, and 3 full per-pass detail.
- `CONV_DIAG_LVL` (integer, legacy alias of `OUT_LEVEL`) -
  per-iteration OUTGEN verbosity. Since 23-Nov-2026 it also gates
  the shared layout in `iter_log_fmt_mod.f90` (`test_iter_output.py`):
  0 = ITER_SUM block + `IT` status banner + `Model Convergence` +
  timing; 1 = adds section banners and per-call status lines
  (solve type, max corrections, luminosity, RUNAWAY top-N table,
  ES_ALI/JCOMP census, JFNK verdicts, scale lines); 2 = everything,
  incl. trace noise (non-thermal species dump, DEC/INC vectors,
  setup chatter, held-fixed notices, BA_READ_V3 progress).
  Exception/error records print at every level. Parser tokens
  ('Current iteration count is', 'POPULATION SOLVE',
  'Luminosity of star', 'Time taken for iteration', 'Maximum
  changes...', 'Maximum %', 'Average number of calls to FG_J_CMF
  is', 'ITER_TYPE', 'Model Convergence:') keep their exact wording
  wherever they print. Round 2 (29-Nov-2026) fixed two printed
  typos: 'Maximm' -> 'Maximum' and 'Current great iteration count'
  -> 'Current iteration count'; every consumer script accepts both
  spellings (plot_cmfgen_diagnostics.py, cmfgen_convergence.py,
  check_convergence.py, check_cmfgen.py, parse_knop_iso.py,
  parse_knop_seed.py, run_model_grid.py). Round 2 also adds the
  moment-refresh pass table + one-line summary, the shared
  run-length compressor for repeated hot-loop messages, the
  per-iteration ISSUES census, and the fixed-width ES11.4E2/ES11.3E3
  numeric helper.

## 7. JFNK keys

`DO_JFNK` (admission gated on `USE_LAM_ES=T`), `JFNK_COUPLED`
(KNOP-only), `USE_JFNK` (deprecated alias), `JFNK_LAG_FACTORS`,
`JFNK_FRZ_GATE`, `JFNK_GATE_NREQ`, `JFNK_FIXT` (fixed-T lock
honoured even on non-Lambda iterations). Implementation:
`cmfgen/new_main/mod_subs/jfnk_solver_v1.f90` and the `jfnk_*` modules;
semantics per `../AGENTS.md` section 8.

WP-1 keys (JFNK acceptance audit, 27-Nov-2026):

- `JFNK_NM_ACCEPT` (logical, default F) - when set, a line-search
  step that passes the GLL reference test with |F| still rising is
  committed like any accepted iterate (the trust region neither
  grows nor restores the base on it). Default F keeps the old
  behaviour: the step is discarded and del halved. Off by default
  pending a production A/B (checklist `jfnk.md` WP-1).
- `JFNK_NM_CREDIT` (integer, default 2, clamped 0..16) - cap on
  consecutive GLL increase accepts; one credit is spent per
  `nm_armijo_increase` commit and refilled on any ared>0 commit.
  Exhaustion falls back to the rejection path.
- `JFNK_STAG_POLICY` (string, default 'ABORT') - Krylov-stagnation
  response when GMRES exits ginfo=1. 'ABORT' keeps the stock gate:
  the inexact direction is used only when both pred_res and the
  measured true_res beat 0.9*ref_norm, else Newton aborts.
  'GRADED' keeps the direction whenever measured true_res still
  beats ref_norm (descent gate decides downstream), and on a
  stale-preconditioner stall (true_res >= ref_norm or
  unverifiable) consumes the refresh flag on the rejection path,
  factors the preconditioner, and re-solves once at 10x tighter
  eta; a second failure aborts. Off by default pending a
  production A/B (checklist `jfnk.md` WP-2).
- `JFNK_GMRES_RST` (integer, default derived) - max GMRES restart
  cycles per solve. The absent-key default is now
  ceil(JFNK_MAXGMRES/JFNK_M)-1 (=4 with the stock 50/10) so the
  effective cap min(JFNK_MAXGMRES,(RST+1)*JFNK_M) equals the
  advertised JFNK_MAXGMRES budget; the old flat default of 1
  silently capped every solve at 20 matvecs. Explicit values win
  unchanged (checklist `jfnk.md` WP-2.2).
- `JFNK_LS_MAX_EVALS_FRAC` (real, default 0.25) - cap on the
  line-search eval spend per Newton pass as a fraction of the
  remaining activation eval budget (JFNK_EVAL_CAP minus
  JFNK_N_RESIDUAL_EVALS). On reaching the cap the LS keeps
  whatever trial already passed Armijo. 1.0 restores the
  unbounded ladder (checklist `jfnk.md` WP-3.4).
- `JFNK_EVAL_BUDGET_HIT` is not a VADAT key -- a solver-visible
  flag set in JFNK_RESIDUAL_V1 when the eval or wall budget
  sentinel fires. The solver maps it to info=2 'budget_stop' so
  the driver commits the best iterate instead of aborting
  (checklist `jfnk.md` WP-3.1/3.2).
- `JFNK_PC_TAU_FAC` (real, default 4, 0=off) - PTC preconditioner
  lag trigger. JFNK_FACTOR_PRECOND_V1 snapshots tau and x at
  factorization; per Newton pass, when
  max(tau/tau_f, tau_f/tau) or the relative inf-norm x drift
  exceeds the factor, jfnk_need_precond_refresh is set so the
  commit path refolds the live tau (checklist `jfnk.md` WP-5.1).
- `JFNK_DESCENT_BASIS` (string, default 'F') - which inner
  product the descent gate tests under PTC: 'F' is the stock
  stripped <F,J_F s>, 'G' the folded <F,J_G s> the step actually
  solves, 'BOTH' accepts on either positive. The Armijo slope
  uses the matching merit; the convergence test stays on |F|
  (checklist `jfnk.md` WP-6.2). Telemetry is the ledger's
  dir_deriv_F/dir_deriv_G columns.
- `JFNK_H_FREEZE` (logical, default F) - freeze the FD JVP step
  size at the value capped from the first Arnoldi direction of
  each GMRES solve, so the matvec is a fixed linear map within
  the solve (the per-direction caps otherwise make it mildly
  nonlinear: hess_res ran 150-450x below true_res). Later
  matvecs whose fresh cap would have shrunk h are counted and
  logged (checklist `jfnk.md` WP-4.2).
- `JFNK_H_FORM` (string, default 'INFNORM') - base FD step
  formula: 'INFNORM' is the stock
  h = eps*max(|x|_inf,1)/|v|_inf; 'TWO_NORM' is the Knoll &
  Keyes 2004 (Sec 3.2) h = eps*(1+|x|_2)/|v|_2;
  'COMPONENTWISE' takes the scalar min over active components
  of eps*max(|x_i|,1)/|v_i| (checklist `jfnk.md` WP-4.3).

## 8. Historical status notes

Dated records of why the machinery above exists. For current
status see `docs/CHECKLIST.md`.

- `adp2_esali` (m15 2.82d SOB+LAM_ES+DO_ES_ALI): clean convergence
  at it851, Model Convergence 100% (MAXCH<EPS=1e-3), ES_ALI J loop
  converging all freqs to dJ/J~1E-4 in ~2 passes.
- warmS run status (post 29-Oct-2026 guards): STABLE but not yet
  converged - outer/inner L ratio ~10 (4.2E9 vs 4.3E8), KNOP
  refresh residual ~4E-4 frac unconverged.
- 21-Oct-2026 KNOP precondition-phase corrections (cold-start
  stall, job 12298049): the cold KNOP start pinned the
  `GS_PRECOND` gate forever - a handful of ill-conditioned rows
  demanded +1E15..+1E72% every LAMBDA call. Fixed by the
  fiddle_pop_corrections scale cap, the BA_LQ quarantine, and the
  N_PRECOND_TOT no-progress semantics (all in sections 3-4
  above).
- 05-Nov-2026 SOB nonloc it64 NaN localisation (job 12689657):
  `GFF_VEC` indexed the 12x11 Gaunt table with `INT(NaN)` after
  `SCRTEMP` carried an all-NaN `POPS(:,11)` column; `STEQ(:,11)`
  was entirely NaN on it64's final solve while every residual
  dump before it was finite - defect inside BA matrix assembly or
  the `CMF_BLKBAND_V3` block solve, not the `DRC` injection into
  the residual (audited and cleared). Fixed by the STEQ
  non-finite scan, the sob_rh_couple DRATE/RH_MAXDJ reporting,
  `SOB_RH_TEST_S1`, and the `MOM_RF_DIV_BULK`/`MOM_XMWORST`
  bounds.
