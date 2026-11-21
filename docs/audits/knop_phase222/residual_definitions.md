# The three residual metrics (source definitions, fg_j_cmf_nonmono_v1.f)

All three use the same per-cell normalization but measure different
differences at different stages of the refresh.

## A. Production solve residual RES -> RES_LS -> KN_RF_RESID

- In MARCH (l.13651+), each cell computes the march update VNEW = T(FLD)
  from the Gauss-Seidel recurrence and measures
      DCELL = |VNEW - VOLD| / MAX(VOLD, 1e-4*VMAX_L)
  where VOLD is the previous iterate at that cell and
  VMAX_L = MAX(MAXVAL(FLD(:,L)), 1e-25) the current field's column peak
  at that frequency (l.13663).
- RES = running max over all cells/frequencies of the pass (l.13675,
  13701, 13758, 13788). RES_FREQ(L) keeps the per-frequency max
  (block-GS diagnostic, l.10116).
- The update applied is relaxed: FLD = VOLD + OMEGA_RAY*(VNEW-VOLD),
  OMEGA_RAY <= 1 per ray (damped to >=0.45 on same-direction growth,
  l.11533).
- So RES per pass is ||T_GS(x_n) - x_n||_max: the Gauss-Seidel map
  defect at the CURRENT iterate, measured cell-by-cell during the
  in-place march. It is the fixed-point defect of the transfer map,
  not a store comparison.
- RES_LS(LS) = the ray's last-pass residual. The printed
  `residual max=` is KN_RF_RESID = KN_PROD_RESID_V1(RES_LS) = the
  production-ray max (probe excluded, l.9185,9291). Value ~6.5e-2
  vs KN_RAY_TOL=2e-3 throughout A1.

## B. Independent q residual KN_Q_LS ("independent q residual max")

- After the production solve exits (l.10622-10690), exactly ONE extra
  pass is applied to a scratch copy FLD_Q of the FINAL production
  field, separately per channel (M then P). The L_QP hook exits after
  one completed map application (l.11461-11469, IXIT=9).
- Same DCELL formula, different evaluation point: the production RES
  samples the defect at each mid-pass iterate; q samples
  ||T_GS(x_final) - x_final|| on the produced field AFTER all
  relaxation/exits.
- Purpose (l.10596-10601): with OMEGA_RAY<1 the committed field is not
  the raw map output, so the production RES understates the true
  fixed-point gap. Probe measurement: last-pass RES=0.696 vs
  q=9.41 (~13x) - the produced field sits where the map defect is far
  worse than the last iterate suggested (anti-contractive / period-2
  map).
- Per-ray argmax kept in KN_QA_*; probe argmax: (raypt 2, freq 1282),
  vold=8.74e-30 vnew=1.03e-28, denom=1e-29 (1e-4*VMAX_L floor,
  vmaxraw=9.25e-30 - floor-dominated column).
- Q distribution over the 84 evaluated rays (probe included in the
  printout): exactly 4 rays >1e-2 - the probe (9.41, the printed
  max) plus exactly 3 production rays, matching the {22,45,50}
  non-committing set; 2 rays >1e-1 total; production median ~5.7e-4,
  p95 ~6.8e-3. So rays 22/45/50 are the only production rays whose
  produced fields still carry a >1% fixed-point defect.

## C. Store-update drift FG_GS_MAXCHG / retained delta FG_RF_DCHG

- Once per frequency per sweep (NEW_FREQ), over every writable
  (model-depth K, ray LS) cell (l.6049-6189):
      T4 = |T1 - T2| / MAX(T2, 1e-4*VMAX_I{P,M}(LS))
  T1 = the freshly solved (committed-candidate) intensity, T2 = the
  retained lagged store (SNP). VMAX_I = the store's ray-column peak.
- Argmax records (K, LS, FREQ_IDX_ITER) into FG_GS_ARGK/R/L; the
  4-deep ring classifies the drift as `cycle` (same cell, alternating
  values = period-2 limit cycle of the composite map), `margnl`
  (same cell, monotone decay), or `bndry` (same argmax depth pinned at
  K<=2 or >=ND-1) - l.7096-7124.
- The store then relaxes toward fresh: store += INEXT_RELAX*(T1-T2),
  omega ~0.513 adaptive (l.6137).
- FG_RF_DCHG = max over rays of the WRITE_STORE/WRITE_STORE_CAP
  committed deltas (l.10860) - bounded by FG_COMMIT_NC=0.05 on capped
  commits. It measures how much the RETAINED store moved this refresh,
  i.e. the lagged-map increment, not the solve's fixed-point defect.

## Can a small committed-store change coexist with a large transfer residual?

Yes, by construction:
1. A stalled ray carries a large solve residual (metric A ~1e-1); its
   commit is capped so the store moves <=5% of the denominator
   (metric C saturates at 0.05).
2. An unsolved ray writes nothing at all: C=0 for that ray while A
   stays large (probe qstore delta=0, vc=7).
3. A converged ray can also show C>0 without any solve problem: its
   committed field drifts because the lagged coupling/source function
   still moves refresh-over-refresh (ray 82: converges each refresh -
   q<1e-2, commits - yet its boundary cells carry the C argmax at
   ~6e-2 every refresh, class cycle/bndry).

The three metrics therefore measure: A = the solve's fixed-point
defect during iteration; B = the produced field's fixed-point defect;
C = the retained store's lagged-update increment. A and B probe the
same operator at different iterates; C is a different difference
committed-vs-retained - bounded by design and blind to whether the
committed field itself is converged.
