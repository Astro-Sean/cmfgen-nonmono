# FG field drift audit: GS_MAXCHG not contracting under DO_NONMONO_KNOP

Audit of `plane/fg_j_cmf_nonmono_v1.f` (FG), `new_main/mod_subs/solve_for_pops.f`
(SFP), `new_main/cmfgen_sub.f` (SUB), commit `cfae5852`.  Scope restricted to the
FG field-update and convergence pathway.  Line numbers verified against the
named commit.

Second-pass verification (this revision): every claim below was re-read against
the live source; B1 and C3 are reclassified, runtime evidence from the FRESH
knop run (it 18-25) is incorporated.

## Pillar 1 - variable flow: what GS_MAXCHG actually is

Metric lifecycle in `fg_j_cmf_nonmono_v1.f`:

- INIT (ML=1, `FIRST_FREQ=T`, `SOLUTION_METHOD='INTEGRAL'` only, 2301-2360):
  zeros `FG_GS_MAXCHG`, `FG_GS_NCONV`, `FG_INEXT_CAPPED`, `FG_SW_*`,
  `FG_RF_DCHG`, `FG_GS_NBAD/NTOT`, `FG_GS_NDIV`, `FG_GS_NBAD_BLK/NTOT_BLK`,
  `FG_GS_NSTAL/NPLAT`, `FG_DIV_HEAL`, `FG_COMMIT_REJ`, `FG_NC_REJ/CAP`,
  `FG_SW_CAPPED/CAPTOT`, `FG_RF_CAPPED`.
- In-sweep (4331-4497, `NEW_FREQ` census only, gated `.NOT. FREEZE_I_NEXT`):
  each I_NEXT store write is relaxed, `V2 + INEXT_RELAX*(V1-V2)`, against the
  `SNP_*` snapshot (4439-4440, 4493-4494); `FG_GS_MAXCHG` accumulates the max
  UNRELAXED relative change `|V1-V2|/max(V2,floor)` (4468-4469), `FG_GS_NCONV`
  counts writes >50%, `FG_SW_CAPPED/CAPTOT` count cap/NaN events vs writes.
- Refresh (per great iteration, SUB:3313-3321, gated
  `(DO_NONMONO_KNOP .OR. DO_SOB_FG_KNOP) .AND. DO_KNOP_REFRESH`): OMP
  `PARALLEL DO` over `KN_NP` rays calls `KNOP_RAY_SOLVE_V1` on opacities
  captured during the sweep (`KN_CHI_EXT/KN_ETA_EXT`, captured 2698-2715 on
  every entry including ES re-entries).  `RES_LS(LS)=MAX(RES_M,RES_P)` (7335);
  diverged rays are excluded from the residual max by the `NDV_LS` mask
  (6517-6521) but counted in `FG_GS_NDIV` (6603).
- Post-refresh snapshot (6597-6607): `FG_SW_MAXCHG`/`FG_SW_NCONV` save the
  in-sweep drift; `FG_GS_MAXCHG = KN_RF_RESID` (max over non-diverged rays of
  the last-pass residual); `FG_GS_NCONV = SUM(NST_LS)` (committed >50%
  writes); `FG_GS_NBAD/NTOT` = last-pass bad-cell census, `NBAD_BLK` =
  best-pass.
- Commit (7217-7300): per-ray jump veto `FG_COMMIT_JUMP=50` (7224, waived for
  never-committed rays via `FG_RAY_SEEDED`, 7219-7220), then conv gate
  `FG_COMMIT_CONV=T` requiring both directions IXIT in {0,1,2} with residual
  <= KN_RAY_TOL or a floored bulk exit (IXIT=2 + `RES<=KN_RAY_BULK_TOL`,
  7230-7238); fail -> `WRITE_STORE_CAP` at `FG_COMMIT_NC=0.05` (7267-7269) or
  reject; pass -> `WRITE_STORE` (unrelaxed, 9588-9595).
- `FG_RF_DCHG` (7274, 7295): committed store vs pre-refresh store delta.
  Jump-rejected rays contribute their JUMP only to `FG_RAY_DCH` (print-only,
  7229) - NOT to `FG_RF_DCHG`.
- SOLVE_FOR_POPS reads all of it via `GET_FG_GS_MATURITY` (1068-1073).

### Variable-flow diagram

```
 formal sweep (per frequency)            refresh (per great iteration)
 --------------------------              ---------------------------
 fresh I_P/I_M  ----+                    captured KN_CHI/KN_ETA
        |          |                            |
        |     |V1 - V2store|/V2store             v
        |     -> FG_GS_MAXCHG (drift)     KNOP_RAY_SOLVE_V1 x KN_NP
        v          |                            |  RES_M, RES_P (dir sweeps)
 relaxed write     |                            v
 T2+w*(T1-T2)  ----+                       RES_LS = MAX(RES_M,RES_P)
 -> I_*_DEPTH_NEXT |                            |  (NDV mask)
                   |                            v
                   |                    KN_RF_RESID = max RES_LS
                   |                            |
                   |     commit: WRITE_STORE / WRITE_STORE_CAP / reject
                   |                            |
                   |                     FG_RF_DCHG (committed delta)
                   v                            v
         FG_SW_MAXCHG <- drift        FG_GS_MAXCHG <- KN_RF_RESID
         FG_SW_NCONV  <- churn        FG_GS_NCONV  <- SUM(NST_LS)
                                    FG_GS_NBAD/NTOT/NDIV/NSTAL/NPLAT
                                    FG_RF_CAPPED  <- last-pass + commit
                                                     caps (audit C1 fix)
                            |
                            v   GET_FG_GS_MATURITY (SFP:1068-1073)
              GS_MAX_ALL = MAX(GS_MAXCHG, GS_SW_MAXCHG, GS_RF_DCHG)
              NBAD_EFF   = MIN(GS_NBAD, GS_NBAD_BLK)
              GS_NCONV_ALL = MAX(GS_NCONV, GS_SW_NCONV)
```

Gate in SFP (1219-1300):
- `GS_MAXCHG` alone must be <= `GS_MATURE_TOL=0.05` for raw maturity; the
  composite `GS_MAX_ALL` must be <= `GS_UNSAFE_TOL=0.5` for the immature flag.
- `GS_FG_IMMATURE` fires on `GS_MAX_ALL>GS_MATURE_TOL` (then re-gated to
  bad-cell census or unsafe tolerance inside the ELSE at 1270-1279), churn
  `GS_NCONV_ALL > MAX(50,NTOT/1000)` (1296-1298), or settle-exceed.
- The gate therefore reads three different quantities max-normed: inner
  residual (a), in-sweep drift (b) which folds population-induced
  displacement (c) into itself, and committed retained-field delta.

## Pillar 2 - adaptive relaxation: the KN_RES_PASS estimator

`KN_RES_PASS(J)` (allocated `2*(KN_RAY_MAXPASS+KN_RAY_HOTPASSES)`, 6469) is
the **max over {rays still running at pass J} x {I_M pass J, I_P pass J}**
of the per-pass cell residual: `KN_RES_PASS(NPO)=MAX(KN_RES_PASS(NPO),RES)`
at 7543 (SWEEP_FLD) and 8480 (SWEEP_FLD_ORD).  Both direction sweeps write
the same slot index, so a slot mixes two fields as well as all live rays.

- CONFIRMED DEFECT A1 - FIXED 10-Nov-2026: the successive-slot ratio
  `KN_RES_PASS(I)/KN_RES_PASS(I-1)` (was 6555-6559) compared maxima over
  *different* ray/direction subsets.  Exit, IXIT, and hot-mode dropout
  (`.NOT. HOTMODE` guard) all changed the sample; a dominant ray
  leaving at pass m made slot m+1 measure a different ray.  It was not a
  contraction estimate for any ray or for a fixed operator.
  Fix: each field solve now keeps a private `RHO_WIN(6)` ring of its
  own post-transient, same-direction pass ratios (alternating: baseline
  NPO-2, admitted at NPO>=5; ordered single sequence: baseline NPO-1,
  admitted at NPO>=4; hot-mode and non-finite samples excluded; a
  diverged pass never reaches the push).  `KNOP_RAY_SOLVE_V1` exports
  the per-solve median as `RHO_M`/`RHO_P`; the refresh pools those over
  non-diverged rays (`NDV_LS=0`, `RHO>0`) and drives omega with the
  pool's 90th percentile - the dominant-mode analogue of the old
  envelope without the subset mixing or single-ray fragility.
  `KN_RES_PASS` is still updated for telemetry but no longer read by
  the estimator.
- CONFIRMED DEFECT A2 (severe, partly mitigated by the A1 fix): the
  estimator prices the *inner* per-ray fixed-point map and feeds
  `INEXT_RELAX` (the `1/(1+rho)` law), which damps the *outer*
  lagged-GS store write (4439-4440, 4493-4494).  The corrected
  estimator now measures a well-defined quantity - the inner map's
  dominant-mode rate - which is the input the original damping law was
  designed for.  Whether the inner rate proxies the outer map's slow
  mode remains a modelling assumption, now explicit: if a model still
  drifts at measured rho << 1, the outer map itself is expansive and
  the remedy is structural, not more damping.
- CONFIRMED DEFECT A3 - FIXED 10-Nov-2026: `KN_RES_PASS` was updated before
  the divergence check, so a diverging pass's RES (possibly NaN or
  >1E3*BEST_RES) was folded in; only the denominator was guarded at 6556.
  `MAX` with a NaN operand is implementation-defined.  Fix: the update now
  sits below the divergence check in both sweeps (post-edit ~7562 and
  ~8497).  Verified: `test_knop_ray.py` 18/18 pass (divergence classes
  T12/T13 unaffected, operator parity T14 unchanged).
- `OMEGA_RAY` (ray-local, adapted per same-direction pair at 7586-7600 in
  SWEEP_FLD with history NPO vs NPO-2, and 8506-8547 in SWEEP_FLD_ORD with
  single-sequence NPO vs NPO-1 plus AA/DLT history rescaling at 8557-8563)
  is correct and independent - it uses per-ray ratios on the same map, is
  floored at 0.45, and damps on same-direction growth >1.05 or on two
  consecutive signed-sum sign flips (period-2 detector, 7609-7619 /
  8516-8544).  The defect is the global `INEXT_RELAX` adaptation only.

### Implemented estimator (A1, 10-Nov-2026)

Data flow: `SWEEP_FLD`/`SWEEP_FLD_ORD` push each post-transient
same-map pass ratio into a private `RHO_WIN(6)` ring (windows are per
field solve; samples pass the same floors the damping branches use so
the seed transient never enters).  At solve end the ring median goes to
the new `KNOP_RAY_SOLVE_V1` outputs `RHO_M` (inward) / `RHO_P`
(outward); a solve that exits inside the transient exports 0 and
contributes nothing.  `FG_KNOP_REFRESH_V1` collects `RHO_M_LS`/
`RHO_P_LS` per ray, pools all samples from non-diverged rays, and takes
the 90th percentile of the sorted pool as rho_est for the existing
`1/(1+rho)` law and `0.5` smoothing.  The `CONV_DIAG_LVL>=2` line now
prints `rho_p90`, `rho_med`, `rho_max`, and the sample count `n_ray`,
so the worst-ray rate is visible next to the applied omega.
Properties: every sample measures the same map on the same ray (a true
contraction sample); per-ray medians make one noisy pass unable to
dominate; rays that exit or diverge leave the pool cleanly; direction
information is preserved in the pool composition (P90 weights whichever
direction produces the slow samples).  P90 rather than plain median:
the lagged-store damping must counter the dominant (slowest/
oscillating) mode the gate reads through `MAX(RES_LS)`, not the
typical ray.  `FG_ADAPT_OMEGA_V2=F` keeps the cross-iteration
`KN_RF_RESID/KN_RF_RESID_PRV` estimate unchanged.
Verified: `unit_tests/python/test_knop_ray.py` 21/21 - new T16/T17
(grinding alt/ord solves export finite per-solve medians; transient
exits export exactly 0) and T18 (end-to-end refresh runs the pooled
estimate, prints `omega est: rho_p90=`, INEXT_RELAX stays in
[0.05,1]); the 18 prior tests unchanged.  `test_knop_lagged_field.py`
all parser checks pass.  Driver `drv_knop_ray.f90` gained cases
`slw`/`slwo`/`rfa` (NCASE 34-36) and prints RHO_M/RHO_P.

## Pillar 3 - fixed-point iteration

Verified correct:

- `SWEEP_FLD` (alternating direction, 7422) and `SWEEP_FLD_ORD` (ordered by
  `sign(GAM)`, 7919) iterate the same map `x <- x + OMEGA_RAY*(F(x)-x)`;
  T3 in `unit_tests/python/test_knop_ray.py` asserts their committed stores
  agree.
- `SEED_FLD` (7399-7413) inverts `WRITE_STORE` exactly (slot s = frequency
  s+1, both iterate LS2=2..NCF); `MARCH` residual `|VNEW-VOLD|/max(VOLD,
  1e-4*VMAX_L)` (9391-9424) matches the WRITE_STORE NSTORE census
  normalization.
- `BUILD_FREQ` (8848) interpolates captured `KN_CHI_COEF/KN_ETA_COEF`
  cubics with exact-node fill from `KN_CHI_EXT/KN_ETA_EXT` (8874-8879);
  frequency-1 edge keeps compressing chi_hat/Q coupling to the lagged
  redder bin matching the in-sweep INIT (8882-8933).
- Discrete transfer equation internally consistent: `chi_hat =
  chi+|GAM|/DLN` (sweep 3372-3375, refresh 8930), `Q = GAM/(DLN*chi_hat)`
  (3468-3472, refresh BUILD_FREQ/SRC_PRIME), `S' = SR + |Q|*I_upwind` for
  GAM<0.  Operator A/B coverage in test T14 (ab/abx/ab2 dumps) asserts the
  two paths build identical CHI_HAT, Q, DCHR, DTAU, WINX, DLNM, A0-A4, SR.
- IXIT codes {0..6} implemented identically in both sweeps.

Defects:

- CONFIRMED DEFECT B1 (latent; INACTIVE on m15 METHOD='ZERO'): the refresh
  clamps DTAU only under `.NOT. KN_INSERT_SAV .AND.
  KN_METHOD_SAV(4:6).NE.'MON'` (9085-9094), while the INTEGRAL-path sweep
  clamps unconditionally for any `METHOD.NE.'ZERO'` (3511-3517).  Under
  the KNOP-required `FG_OPT=INT/INS` (forced at rd_control_variables:3399)
  `INSERT=T`, so for any non-ZERO METHOD (e.g. LINMON/LOGMON) the refresh
  solves an unclamped-DTAU operator the sweep did not.  The 30-Oct comment
  at 9075-9083 mis-cites the sweep's skip: `INSERT .OR. METHOD(4:6)='MON'`
  exists only in the DIFFERENCE path (2818), which the refresh can never
  mirror (refresh is guarded `SOLUTION_METHOD='INTEGRAL'`, 6404-6408).
  With METHOD='ZERO' both paths take the trapezoid branch (3500-3504 vs
  9063-9067) and agree, so this does not explain the current drift.
- SUSPECTED B2 (bounded): `KN_RF_COUPLING` requires `IKC1|IKC2` (embeds
  `J_PNT` range check, 6998-7013); `KN_COUPLING_OK` only requires
  `COMP_DEPTH_IDX>=1` (3436, 3468).  A compressing point with out-of-range
  `J_PNT` couples in the sweep but not the refresh (operator asymmetry on
  edge cells only; `J_PNT` is range-validated on the `LS<=NP-I+1`
  triangle, 1999-2015).
- CONFIRMED minor: `FG_JUMP_MAX` reads all rays' store for `VMAX_X`
  (9780-9788) while sibling OMP threads may be writing their own rays' -
  a benign but real OMP race on the store.
- CONFIRMED minor: `PK_DTAU(NI)` copies an uninitialized `DTAU_L(NI)` -
  never read anywhere, benign.
- CONFIRMED minor: `KN_CHI_EXT/KN_ETA_EXT` are `REAL` (single) while the
  captured `CHI_EXT/ETA_EXT` are `LDP` - a deliberate but real precision
  downgrade of the refresh's opacities (re-widened via `DBLE` at 8877).

## Pillar 4 - store commit pathway

- CONFIRMED DEFECT C1 (gate-relevant) - FIXED 10-Nov-2026: `FG_RF_CAPPED`
  was `SUM(NCP_LS) + SUM(NNG_LS) + SUM(NDV_LS)` (was 6617) - it mixed
  three units: `NCP_LS`/`NNG_LS` are all-pass accumulated per-cell-write
  event counts (NCAP/NNEG incremented per write across every pass) and
  `NDV_LS` is a 0/1 diverged-ray flag.  SFP divides by
  `GS_CAP_DEN=GS_NTOT` (a single-pass cell census) at 1237-1248 - an
  all-pass event numerator against a one-pass denominator.  Transient
  clamps on early passes could hold `GS_UNSAFE` (>10%) / `GS_PIN`
  (>0.1%) latched when the committed field is clean.
  Fix: `KNOP_RAY_SOLVE_V1` now exports `NCAP_LST`/`NNG_LST` = last-pass
  cap/floor events of each direction's solve (per-pass counters
  `NCP_P`/`NNG_P` reset at each pass top, captured per direction at the
  sweep tail) plus the commit-phase events (`NCP_C`, bumped in
  `WRITE_STORE`/`WRITE_STORE_CAP`).  The refresh aggregates them via
  `NCL_LS`/`NNL_LS`; `FG_RF_CAPPED = SUM(NCL_LS)+SUM(NNL_LS)` is the
  same-census numerator, `FG_INEXT_CAPPED` keeps the cumulative all-pass
  sum for checkpoint comparability, and `FG_GS_NDIV` tiers diverged
  rays separately as before.  Telemetry line `same-census cap/floor=`
  added.  Verified: 22/22 knop_ray tests including T19 ('capx' under
  NONMONO_INEXT_CAP=1e-4: NCAP=168 all-pass events vs NCAP_LST=84
  last-pass events <= NTOT=90 over NPASS=2); all lagged-field parser
  checks pass.
- CONFIRMED DEFECT C2: `HEAL_STORE` counts already-zero slots as healed in
  partially sane rows (9867-9871: `STORE>0` fails on a zero slot so the
  branch fires and `NHL++` on a no-op).  Overstates `FG_DIV_HEAL`.
- RECLASSIFIED C3: `FG_RAY_SEEDED` is set only on `WRITE_STORE` (7293) and
  on the population-path `RESTORE_FG_INEXT_GS` (5656-5657).  A
  capped-commit ray keeps the seed-fill waiver, and the OBS remap
  (`RESTORE_FG_INEXT_OBS`, 6100-6110) leaves all rays unseeded - but the
  conv gate still applies on every commit (JUMP only bypasses the
  magnitude veto, 7214-7215), and a remapped/unvalidated store is exactly
  what the waiver exists for.  Bounded exposure, mostly by-design; keep as
  minor review note, not a defect.
- SUSPECTED C4: a ray whose converged fixed point is >50x its store is
  jump-rejected every refresh; `RES_LS` reports converged but the store
  stays at the sweep-relaxed state - a blind spot invisible to every gated
  metric (`FG_RAY_DCH` is print-only, 7229).
- CONFIRMED minor: `WRITE_STORE_CAP` commits bounded by `FG_COMMIT_NC=0.05`
  per cell per refresh (9654-9655); a systematic bias ratchets
  `(1+NC)^n` toward the fixed point but cannot overshoot it in one step.
- CONFIRMED minor: `NI<2`/`NCF<2` early return (6967) reports `RESID=0`
  (converged) but `NXT=-1` (unsolved) - inconsistent census.  Live runs
  show `unsolved=1` every refresh; almost certainly this ray.

## Pillar 5 - population-to-field coupling

- The refresh re-solves on opacities captured during the sweep
  (`KN_CHI_EXT/KN_ETA_EXT`, FG:2698-2715, written on every entry including
  ES re-entries).  Internally consistent per call; populations cannot move
  between capture and measurement, so `GS_MAXCHG` post-refresh is a pure
  inner residual.
- CONFIRMED: `GS_MAX_ALL` max-norms inner residual, sweep drift, and
  committed delta.  Under the nested-iteration design this is intended -
  but the single reported `GS_MAXCHG` cannot separate (a) genuine inner
  floor from (b) population-induced re-lagging.  The `sw/res/rf` columns in
  `GS_DRIFT_DIAG` (SFP:2668-2688) are the intended discriminant.
- CONFIRMED DEFECT D1 (moderate): JFNK eval sweeps
  (`JFNK_RAD_SWEEP_V1`, SUB:5873-7464) reset the entire `FG_GS_*`/`FG_SW_*`/
  `FG_RF_*` metric set at FG INIT (2301-2360) without bracketing - only 5
  scalars are saved/restored (MGSMAX/MRFPRV/JNRLX/INRLX/KRFPRV at
  5961-5965, restored 7390-7394).  `FREEZE_I_NEXT` blocks the eval's
  store writes AND its metric accumulation (4331-4334), so post-eval the
  census sits at zero.  `DUMP_FG_INEXT_GS` then writes `FG_GS_MAXCHG`,
  `FG_INEXT_CAPPED`, `FG_GS_NCONV`, `FG_GS_NDIV`, `FG_COMMIT_REJ`,
  `FG_DIV_HEAL`, `FG_SW_MAXCHG`, `FG_SW_NCONV`, `FG_RF_DCHG` (5458-5466) -
  post-eval zeros - into the checkpoint; a restart reads
  `GS_CHG_FILE=0 < GS_MATURE_TOL` and `GS_CAP_FILE=0` and sets
  `GS_RWARM_SKIP=T` (SFP:1090-1097), skipping the `N_LAM_RWARM` warm-up
  the cold metric was meant to force.  Live gate unaffected (all gate
  inputs re-derive each iteration); restart-path only.
- CONFIRMED latent D2 (INACTIVE on m15): dielectronic `COMP_J_BLANK` calls
  (SUB:2073, before the main sweep loop at 2400) reach
  `FG_J_CMF_NONMONO_V1` (dispatch is section-blind, comp_j_blank:875-897)
  and would tick `FREQ_IDX_ITER` (2595-2596) and write store slots at
  off-grid frequencies - the store slots the continuum sweep then reads
  as lagged fields.  Inactive here: `NDIETOT=0` on m15 (all `DIE_*` VADAT
  flags F, `DIE_AS_LINE=F`).  Also a `FIRST_FREQ` hazard: `FIRST_FREQ` is
  unset before the continuum loop (only assigned at 2406/2408 inside the
  DO 10000 loop), so the dielectronic calls inherit a stale value.
  Latent defect for models that register dielectronic lines.
- CONFIRMED minor: `GS_NBAD_STAL` exported via `GET_FG_GS_MATURITY`
  (SFP:1072) but never consumed in SFP (telemetry only; `GS_NBAD_PLAT` IS
  consumed at 1345).
- CONFIRMED minor: `JFNK DEFER` block prints `GS_MAXCHG` under the label
  "FG sweep drift" (SFP:2608) - post-refresh it is the refresh residual.

## Runtime evidence (FRESH knop run, m15 2.82d, it 18-25)

Refresh residual history over 25 refreshes: 0.015-0.084 band, no monotone
contraction; the last two refreshes record 5.9E-2 and 8.4E-2 - above
`GS_MATURE_TOL=0.05` and far above `KN_RAY_TOL=2E-3`.

Per-refresh census (stable pattern):
- committed 76-81, jump-rej 2-4, nc-rej 0, capped 0-5, div 0, unsolved 1.
- Exit classes: exact ~48, tol ~9-10, bulk ~24-26, stall 0-1, cap 0-1,
  plat 0-2.
- Worst ray migrates: 12 -> 25 -> 24 -> 1 (the center ray, NI=338) for the
  last two refreshes, with `ixit=6` plateau exits and `v=3` capped
  commits; `worst-ray delta: commit=5.0E-2` (the FG_COMMIT_NC bound) while
  `max` printed 4.4E+02 on the jump-rejected channel.
- `INEXT_RELAX` drifted 0.513 -> 0.575 across the same span while the
  residual rose - the adaptive estimate is not tracking contraction.
- `retained-field delta` oscillates 0.65-18.3.
- Population side independently forced LAMBDA (preconditioning demand
  1.0E3-2.2E3 > limit); GS freeze engaged at it 25 (frozen=T, held=T).

Interpretation: not a lone diverged ray.  A small rotating set of
marginal rays (bulk/plateau exits on center and mid-impact rays) carries
the residual floor, while jump-rejected and capped commits keep the
retained store moving.  The drift is genuine inner non-convergence on
marginally contractive cells plus store motion the metrics partially
hide.

## Causal path to persistent drift

`FG_GS_MAXCHG` pins at `KN_RF_RESID` = max over non-diverged rays of
`MAX(RES_M,RES_P)` (6599).  Non-exclusive mechanisms:

1. Stall/cap/plateau rays leave residual ~0.05-0.5 -> the max-norm stays
   above `GS_MATURE_TOL` while the bulk is converged.  Genuine inner
   non-convergence on marginally contractive cells (the live pattern:
   worst ray rotates between a handful of rays, plateau exits on the
   center ray).
2. The broken `KN_RES_PASS` estimator feeds a wrong global `INEXT_RELAX`
   (cross-subset ratios -> omega drifting 0.51-0.58 while residual rises)
   -> the store the next sweep seeds from drifts or stalls -> the same
   slow rays never converge -> the metric never falls.  The defective
   estimator sits on the causal path.
3. Rejected/capped channels move the store in ways the gated metrics do
   not fully see (C4 blind spot; `FG_RF_DCHG` misses jump-rejected rays).

## Diagnostics gap analysis (Pillar 6)

Already present (refresh block, 6608-6790 + GS_DRIFT_DIAG SFP:2668):
residual max; per-IXIT exit census + unsolved; committed/jump-rej/nc-rej/
capped/div verdict counts; diverged rays; clamp events; bad/best-pass
cell censuses; retained-field delta; INEXT_RELAX; worst-ray identity,
residual, passes; worst cell i/l + same-cell-passes (within a refresh);
per-pass max history (PRINT_KNOP_DEBUG only); commit census with
committed-max delta.

Missing vs the required list:
- Median refresh residual; 90th/99th percentile (only max is kept).
- Number of rays with `RES_LS > KN_RAY_TOL` (exit census approximates but
  does not isolate it).
- ~~Median same-ray contraction ratio~~ - provided by the A1 fix:
  `CONV_DIAG_LVL>=2` prints `rho_p90/rho_med/rho_max/n_ray` each
  refresh.
- Inward vs outward split: `RES_LS=MAX(RES_M,RES_P)` collapses both; no
  per-direction median/max.
- Whether the worst cell persists ACROSS refreshes (same-cell passes only
  tracks within a refresh).
- Sweep drift maximum in the refresh block (FG_SW_MAXCHG printed only in
  GS_DRIFT_DIAG, which is off by default).

Design: extend the refresh summary under `CONV_DIAG_LVL>=1` (headline
counts already there) with a compact `res stats: med/p90/p99/max`,
`over-tol rays=N`, `dir M/P max=..`, `med rho=..`, `worst-cell streak`
fields - all computable from `RES_LS`, `RES_M/P` per ray (need to
return both directions, not just the max), and a persistent
worst-cell id.

## Pillar 7 - test coverage vs required

`unit_tests/python/test_knop_ray.py` (drives production
`KNOP_RAY_SOLVE_V1`/`FG_KNOP_REFRESH_V1` through `drv_knop_ray.f90`):

| required | coverage |
|---|---|
| 1 manufactured fixed-opacity converges | T1, T9, T10, T14 (conv/rev/strg/ab) |
| 2 converges from zero and perturbed seeds | partial: cold seed + T6 near-fixed-point; add explicit zero-seed and +x% perturbed-seed runs |
| 3 same-ray contraction correct | T16/T17: per-solve median ratios exported; transient exits export 0 |
| 4 max-envelope differs from ray-local | partial: T18 exercises the pooled estimator end-to-end; a two-ray driver swapping dominance across passes would pin the contrast directly |
| 5 period-2 oscillation detected | missing (SGN_ACC path exists, no test) |
| 6 rejected ray does not corrupt store | partial: T12/T13 diverged fail-closed; add a jump-rejected (v=1) case asserting committed slots unchanged |
| 7 ord/alt same fixed point | T3 |
| 8 sweep/refresh operators agree | T14 ab/abx/ab2 |
| 9 store converges at fixed opacities | partial: T8 single refresh; add multi-refresh loop asserting store delta -> 0 |
| 10 outer iteration distinguishes field vs population drift | missing - SFP gate-level test |

Additional pins worth adding: `FG_RAY_SEEDED` on capped commits (C3),
DTAU clamp parity under METHOD<>'ZERO' (B1 fix verification).

## Proposed corrections (not yet implemented, ranked)

| # | Defect | Fix sketch |
|---|---|---|
| A1-A3 | KN_RES_PASS estimator | DONE 10-Nov-2026: per-solve `RHO_WIN` rings -> `RHO_M/RHO_P` -> pooled P90 drives omega; KN_RES_PASS update moved after the divergence check (A3); KN_RES_PASS is now telemetry-only. |
| C1 | Cap-share denominator | DONE 10-Nov-2026: `FG_RF_CAPPED` now sums last-pass + commit-phase cap/floor events (`NCL_LS`/`NNL_LS` per ray); `FG_INEXT_CAPPED` keeps the all-pass cumulative; `NDV` tiered via `GS_NDIV`. T19 pins the census bound. |
| B1 | DTAU clamp | Mirror the INTEGRAL sweep rule: clamp whenever `KN_METHOD_SAV.NE.'ZERO'`; fix the mis-citing comment. Inactive on m15 (METHOD='ZERO') - low priority for this run. |
| D1 | JFNK eval metrics | Bracket the full FG_GS_* set around eval sweeps (or skip the checkpoint metric write after an eval touch). |
| C4 | Rejected-ray blind spot | Fold rejected-ray deltas into a gated metric or a dedicated census. |
| C2 | Heal census | Skip already-zero slots in NHL. |
| 7 | GS_NBAD_STAL | Consume it in a tier or drop the getter arg. |
| 8 | JFNK DEFER label | Print GS_SW_MAXCHG under "FG sweep drift". |

## Findings ranked by severity

1. A1-A3 (severe, FIXED 10-Nov-2026): the adaptive-omega estimator was
   not a contraction estimate and could feed NaN into INEXT_RELAX - on
   the drift causal path.  Replaced by the per-ray pooled-P90 estimator
   described above; 22/22 knop_ray tests, all lagged-field parser
   checks pass.
2. C1 (severe, FIXED 10-Nov-2026): cap-share numerator/denominator unit
   mismatch could hold GS_UNSAFE/PIN latched on transient events.  Live
   runs show clamp events=0, so this was not the active cause of the
   m15 drift; the fix removes a latched-false-positive risk and makes
   the gate read true last-pass shares (T19: NCAP=168 vs NCAP_LST=84).
3. D1 (moderate): JFNK eval sweeps zero checkpoint metrics -> restart
   skips the LAMBDA warm-up on a stale 0 file metric.
4. C4 (moderate): jump-rejected converged rays are invisible to gates.
5. B1 (latent, inactive on m15 METHOD='ZERO'): refresh/sweep DTAU
   operator mismatch under INT/INS + non-ZERO METHOD; mis-cited comment.
6. C3 (minor, mostly by-design): capped commits and the OBS remap leave
   rays unseeded -> jump veto waived on next commit.
7. B2 (low): J_PNT gating divergence between sweep and refresh on edge
   cells.
8. C2, D2 (latent, NDIETOT=0 on m15), NI<2 census, GS_NBAD_STAL, label
   defect, REAL capture, PK_DTAU(NI), OMP race (minor).

## Nov-2026 critical-source-audit cycle (main 7220b616 -> f8c04ca3)

An independent source audit was verified issue-by-issue against
`main` at `7220b616`; every confirmed item was corrected in its own
phase with regression coverage in test_knop_ray.py (T19-T32) and
drv_knop_ray.f90 (sibref/jmpx/jmp2/om*/ni0/ni1/heal/esav/pmon/pnom).

- FG_JUMP_MAX OpenMP race (CONFIRMED): sibling store columns were
  read inside the parallel ray loop, and the `.EQ.` NaN test let
  +Inf through. Immutable per-frequency sibling references
  (`JMP_SIB_M/P` via `KN_JMP_SIBREF_V1`) are now built before the
  ray loop; verdicts/stores bit-identical at 1/2/4 threads.
- INEXT_RELAX controller (valid criticism): the pooled same-map
  ratio is an inner-map diagnostic used empirically for outer
  damping; `KN_OMEGA_POOL_V1` now retains the previous omega below
  KN_RHO_MIN_SAM=10 samples (the P90 index degenerates to the pool
  max below that) and reports coverage.
- Jump-rejected false maturity (CONFIRMED): verdict-1/2 rays
  contributed nothing to FG_RF_DCHG; new FG_RF_PEND folds the max
  pending rejected delta into RF_DCHG at GET_FG_GS_MATURITY.
- Unsolved early returns (CONFIRMED): NI<1/NCF<2 now fail closed
  (NXIT=-2, RESID=1.0); NI=1 remains a plain unsolved census entry.
- HEAL_STORE census (CONFIRMED): zero slots no longer counted.
- JFNK eval telemetry (CONFIRMED): eval sweeps' INIT pass cleared
  the FG census; FG_EVAL_STATE_PUSH_V1/POP_V1 now bracket the whole
  telemetry set in JFNK_RAD_SWEEP_V1.
- B1 DTAU (latent defect, was miscited): BUILD_FREQ skipped the
  clamp under INSERT/MON per the DIFFERENCE-path rule; the
  INTEGRAL-path sibling clamps unconditionally, so the refresh now
  does too. Inactive under METHOD='ZERO'; pmon/pnom now identical.

Unit status: test_knop_ray.py 36/36 plus all pin suites green.
Model-level verification through run_model_grid.py is still owed
before the drift question can be called closed on production runs.

## Follow-up audit corrections (e6e8f47a -> next)

- KN_JMP_SIBREF_V1 alloc failure no longer masquerades as cold
  start: ISTAT out-arg (0/1/2), REF=-1 means cold start only.
- FG_RF_FAILED flags an aborted refresh (sibref or work-array
  allocation failure); GET_FG_GS_MATURITY exports RF_DCHG>=1.0 so
  the population gate fails closed instead of releasing on small
  sweep metrics.
- DTAU parity hardened: KN_DTAU_INT_V1 is now the single shared
  cubic+clamp operator behind both the INTEGRAL sweep block and
  BUILD_FREQ (T37 pins bounds on monotone and non-monotone chi).
- Eval bracket ordering: push moved after the last entry early
  exit; no exit path can bypass pop (T36 pin).
- test_knop_ray.py 41/41; model-level run still owed.
