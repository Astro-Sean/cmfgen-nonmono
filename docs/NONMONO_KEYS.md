# NONMONO_KEYS.md - VADAT and OBS key reference (non-monotonic fork)

Per-key reference for the control keys added or modified by this
fork. Sources of truth: `new_main/mod_subs/control_variable_mod.f`
(declarations), `new_main/mod_subs/rd_control_variables.f`
(parsing + auto-sets), `new_main/mod_subs/check_param_consistency.f`
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
  (`new_main/subs/recompute_sigma_true.f`) fills it at model setup,
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
  to upstream for sigma >= -1.
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
  pathway (`plane/fg_j_cmf_nonmono_v1.f`,
  `plane/mom_j_cmf_nonmono_v1.f`). Reference/diagnostic only for
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
Dispatch in `new_main/cmfgen_sub.f` (VELTYPE):

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

Grid-side comparison on the real m15 structure (ND=70, production
drivers, `unit_tests/work/grid_ab`): the tau_es grids put ~80% of
depths below V=100 km/s; VEL_LAW=14 trades most of that for 5-8x
the resolution in the 300-2000 km/s shell and the ejecta, and
captures the CDS sigma spike deeper (-82 vs -63/-36, source -157).

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
- Dead-key warning (25-Nov-2026): `MOM_BLKGMRES`, `MOM_NUSOLV` and
  `MOM_H_IMPL` only do work inside `MOM_J_CMF_NONMONO_V1` (SOB +
  `USE_LAM_ES=F`). Set T on any other route - production
  `USE_LAM_ES=T`, KNOP, `DO_SOB_FG_KNOP`, or stock - and
  `CHECK_PARAM_CONSISTENCY` prints an inert warning; the run is
  legal but the keys are dead in OUTGEN.
- `KN_NI1_FULL` (default F, 24-Nov-2026): certify a production
  NI=1 ray KN_CS_FULL (reason 0) in `KNOP_RAY_SOLVE_V1` instead
  of leaving it KN_CS_INV. A single-point ray has no marchable
  transfer domain, so vacuous FULL is the honest class; with F
  the ray holds `FG_RF_ALL_CONV` FALSE forever and, through the
  rf_not_allconv fold-in, the LAMBDA handoff shut.
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
- `DO_KNOP_SPEC_GATE` (default T): `obs/knop_spec_gate.f`
  classifies GS_STATE + stall marker + manifest before the
  OBSFLUX/OBSFRAME writes; non-converged KNOP models stop under T
  or render a marked diagnostic under the F opt-out
  (`test_knop_spec_gate.py`, 9 cases).
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
`new_main/mod_subs/jfnk_solver_v1.f90` and the `jfnk_*` modules;
semantics per `../AGENTS.md` section 8.

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
