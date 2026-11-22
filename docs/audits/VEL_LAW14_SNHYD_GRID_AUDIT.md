# VEL_LAW=14 (RV_SN_MODEL_SNHYD TAU_OPT=1) robustness audit

Date: 22-Nov-2026. Auditor: Devin. Trigger: user request - VEL_LAW=14
is deployed as the upgrade to VEL_LAW=12 on the m15 SN_HYDRO_DATA runs;
verify it is still robust.

Scope: `new_main/subs/rv_sn_model_snhyd.f` (TAU_OPT=1 branch),
its VADAT dispatch (`cmfgen_sub.f:1094`), the flag bookkeeping
(`rd_control_variables.f:251`), and the unit battery
`unit_tests/python/test_rv_snhyd.py` driving the production objects.
The earlier Audit-v2 fixes (U_GAP_FAC floor, ordering assert,
commit 58a496dc) are included in the reviewed code.

## What VEL_LAW=14 changes vs VEL_LAW=12

- Data source: `SN_HYDRO_DATA` (as V13) instead of `input_hydro.dat`;
  `SN_MODEL`/`SN_HYDRO_MODEL` stay set so density/abundances come via
  `RD_SN_DATA` and VADAT RSTAR/RMAX are rewritten to the hydro edges -
  the V12 bookkeeping is preserved.
- Grid coordinate: not uniform in log(tau_es) but uniform in the
  blend U = integral of `dlog_tau + ALPHA_R*|dlogR|` with log_tau
  clamped below at `TAU_LO = max(ltau_2, ltau_max - LOG_DEC_CAP)`
  (`rv_sn_model_snhyd.f:284-291`). The tau is built from the file's
  Kappa (Rosseland) column times density, not constant kappa_es.
- Structure-aware insertion: the nearest uniform-U target is moved
  onto each velocity-reversal entry, each local sigma minimum with
  sigma < -1, and each density spike (>1.5x both neighbours), subject
  to the U_GAP_FAC=5% cell floor (`:305-334`).
- V(R) on the new grid comes from the Steffen monotone cubic
  (`MON_INT_FUNS_V2`) with SIGMA = cubic'(R)*R/V - 1 analytically,
  preserving sign-changing gradients (`:406-419`); V12 used LIN_INTERP
  plus a two-sided difference.
- Same NBND_INS=3 boundary packing fractions as V12 (1/20, 1/8, 1/3).

## Verified (evidence)

1. Construction matches the independent Python replica
   (`_predict_ugrid`/`_eval_v_sigma` in test_rv_snhyd.py)
   node-for-node: 4.9e-12 on the synthetic double-shock profile,
   2.3e-14 R / 9e-10 V / 9e-11 SIGMA on the real m15 SN_HYDRO_DATA.
2. All 11 guard cases abort correctly (ND bound, missing file,
   header counts, non-decreasing input and expanded R, non-positive
   V, non-positive expansion factor, all-zero Kappa under
   TAU_OPT=1, missing mass-fraction section, zero leading density).
3. The expansion step is mass-conserving: R -> R + V*dt and
   rho -> rho/(1 + dt*dV/dR)/(1 + dt*V/R)^2 with the file Sigma
   supplying dV/dR; endpoints pin RMAX/RCORE to the expanded edges.
4. Feature insertion lands exactly on hydro knots: the CDS
   sigma-minimum (file J=65, sigma=-156.9) places output depth 62 at
   R = 1.4624650000000001e4 - the file value to the last ulp.
5. The blended coordinate is exactly ltau-shift-invariant (verified
   to 6.6e-14 under kappa x7 and x0.001 rescales) because U is built
   from differences of ltau. Global kappa rescaling cannot move the
   grid.
6. Kappa steers the grid where it should: x20 inner opacity pulls
   inner-region depths 15 -> 18 on the real file (and 14 -> 17 on
   the synthetic beta profile); x0.05 suppression reverses it. The
   response exists only where cumulative tau exceeds TAU_LO - below
   the cap the blend is pure logR by construction (a 30x CSM kappa
   boost changed nothing, as designed).
7. Production: the converged m15_2p82d_VAL_sob_diag_ngoff RVTJ grid
   equals the driver output to 2.0e-10 (R) and 4.6e-8 (V) - the
   production model ran exactly this construction.
8. Edge sizes: NX=3, ND=10, and ND>NX (80 depths on 70 hydro
   points) all complete with ordered grids and pinned endpoints.
9. The upgrade claim holds on the real file: TAU_OPT=0 (the V12
   recipe on our data) puts 61 of 70 depths in the V<200 km/s CSM
   and 7 in the V>1000 ejecta; TAU_OPT=1 puts 36 and 19 - the CSM
   no longer starves the shell/ejecta.

## Confirmed defects (pinned failing in test_rv_snhyd.py)

- D1: a non-finite interior Kappa produces a valid-looking but
  shifted grid, silently. The `MAX(KAPPA_HYDRO(I),T2)` floor absorbs
  the NaN (ifort's MAX returns the second argument when the first
  comparison fails), so the corrupted cell quietly takes the floor
  opacity. Probe: NaN at the midpoint, rc=0, RESULT_OK, grid
  shifted ~3e-4-7e-3 relative.
- D2: a non-finite interior Density, Velocity, or Radius crashes
  hard (rc=174, no diagnostic). NaN propagates through kappa*rho
  into TAU_HYDRO, TAU_LO, and U_HYDRO; LIN_INTERP's bracket walk
  then reads out of bounds. Every ordering/sign guard passes on NaN
  because Fortran comparisons against NaN are all .FALSE. - no
  finiteness screening exists on any input column.
- D3: a negative interior Density is accepted silently: only
  DENSITY(1) is screened (needed to end the section read).
- D4: a file with no Sigma section runs TAU_OPT=1 with
  SIGMA_HYDRO = 0 - feature detection finds nothing and the
  expansion factor silently assumes homology. No warning.

Fix candidate (awaiting sign-off, one-change-at-a-time): an
explicit finiteness/positivity screen over the read columns at
parse time (~15 lines, aborts with a named column/index), plus a
warning when TAU_OPT=1 sees an all-zero Sigma column.

## Limitations / notes (not defects)

- LOG_DEC_CAP=4.0 acts in natural-log units: TAU_LO sits e^4 (~55x)
  under tau_max, not 4 decades. The comment said "decades"; it is
  corrected in this edit. Sensitivity measured via the replica:
  a literal 4-decade cap (9.21 ln units) would shift the census to
  49 CSM / 13 ejecta depths vs the code's 36 / 19 - i.e. back
  toward the V12 failure mode. The code value is load-bearing and
  the better choice. On the m15 file the cap lands at tau ~0.69 -
  near the photosphere by coincidence of tau_max ~ 37.
- The feature detector is blind to (a) monotone density steps -
  the real CDS is a 6-cell 100x rise that never exceeds 1.5x both
  neighbours - and (b) reversal exits (V local maxima, sigma
  crossing -1 upward). Only entries and sigma minima insert. On
  m15 the CDS is still resolved via its sigma-minimum and the
  file's own boundary packing; a shell with a monotone rise and a
  broad shallow dip could under-resolve.
- Every sigma on the real m15 file is < -1 (V peaks at the inner
  edge ~7074 km/s and declines outward to 49 km/s): "reversal
  entry" can never fire here; only local sigma minima insert.
  The field is globally decelerating, which is the regime the
  non-monotonic solvers must carry.
- Tightest output cell is 3.0e-5 (inner edge: the file's own
  boundary packing amplified by the /20 insertion); TAU_OPT=0 is
  worse (3.2e-7 at the outer edge). The dR/R < 1e-4 report fires
  as designed; fractions are the inherited V12 recipe.
- On the m15 file the Kappa column is nearly constant
  (0.290-0.327), so the realised gain vs V12 comes from the blend
  and insertion, not Rosseland-awareness; the kappa machinery is
  verified live but idle on this model.
- Output sigma at the CDS knot reproduces -82 vs the file's -157:
  the Steffen limiter bounds the knot derivative by the shallower
  adjacent secant - inherent regrid dilution, not a bug
  (SIGMA_TRUE is recomputed on the new grid downstream anyway).

## Verdict

PARTIAL->VERIFIED with pinned defects. The TAU_OPT=1 construction is
mathematically correct, replicated independently node-for-node, and
is what production actually ran; the blend plus feature insertion
measurably fixes the V12 under-resolution of the ejecta, and the
TAU_LO cap is load-bearing in the right direction. Robustness gaps
are confined to input validation: non-finite or sign-flipped
interior data either crashes opaquely (D2) or returns a quietly
wrong grid (D1/D3), and a missing Sigma column silently disables
the feature machinery (D4). None are reachable from the current
converter output, but a corrupt or hand-edited SN_HYDRO_DATA fails
without diagnostics. The D1-D4 fixes are deferred pending user
sign-off; the tests pin them until then.
