# KNOP remediation WP-6 -- forcing-term field tolerance (KN_FORCING_TOL)

Date: 26-Nov-2026. Branch: knop/wp6-forcing-tol. Status: LANDED
(flag default F; armed arm is the human A/B).

## What the package does

The Knop refresh is the inner solve of an inexact-Newton loop: the
outer map is populations -> opacities -> field -> populations, and
its residual is the population demand MAXCH. WP-6 prices the
refresh's inner stopping tolerance from the outer-residual ratio
with the Eisenstat & Walker 1996 (SIAM J. Sci. Comput. 17, 16-32)
Choice-2 forcing term:

    eta_k = GAM * (||F_k|| / ||F_{k-1}||)^ALP
    safeguard: eta_k = max(eta_k, GAM*eta_{k-1}^ALP) whenever that
      floor exceeds 0.1
    eta capped at 0.9
    KN_RAY_TOL_EFF = clamp(eta_k * TMAX/0.9, TMIN, min(TMAX, CAP))

with the recommended GAM=0.9, ALP=2. The scale TMAX/0.9 maps the
forcing term's [0,eta_max] band onto the tolerance band. A loose
inner solve when the outer residual is large is exactly what E-W
prescribes (Arramy et al. 2024 run inner tolerances near 1e-2 the
same way); a fixed 2e-3 for every iteration is not.

## Where it lives

- KN_RAY_TOL_EFF (module var, init 0 = unset) is the operative
  tolerance everywhere the refresh machinery priced KN_RAY_TOL:
  certification (KN_CERT_CH), the SWEEP_FLD / SWEEP_FLD_ORD pass
  exits, the MARCH bad-cell census and BID queue, NOP_CELL and
  TRY_HOT thresholds, the refresh report line, the KN_RAY_DUMP and
  topo-dump records (so replays reproduce the armed state).
- FG_KNOP_REFRESH_V1 sets KN_RAY_TOL_EFF=KN_RAY_TOL each call,
  then KN_FRC_UPDATE() reprices it when armed. A driver reaching
  KNOP_RAY_SOLVE_V1 directly hits the unset sentinel and falls
  back to KN_RAY_TOL -- flag-off and driver paths are identical
  to pre-WP-6. KN_CERT_CH is PURE and cannot arm the sentinel,
  so it resolves the unset state to a local LTOL=KN_RAY_TOL
  instead; a direct KN_CERTIFY caller that never passed the
  solve sees the configured tolerance (caught by the csut
  certification table, T39, before this fallback existed).
- SET_FG_FRC_RES(MAXCH) is called by CMFGEN_SUB ahead of each
  refresh. The MAXCH in scope there is the previous call's
  population demand, so the ratio is one outer step lagged --
  the standard practical E-W evaluation point. If a refresh
  no-ops, the stored value is overwritten next call and the
  ratio spans two iterations; bounded and rare.
- KN_FRC_CAP = 0.1*GS_MATURE_TOL (invariant I2 bound) arrives
  through SET_FG_FORCING_TOL because the FG module cannot see
  CONTROL_VARIABLE_MOD; unset -> TMAX alone caps (drivers).
- ITER_LEDGER v4 appends frc_eta, frc_rat, frc_sg, frc_tol.

## Flag and keys

- KN_FORCING_TOL (LOGICAL, default F): arm the forcing term.
- KN_FRC_GAM / KN_FRC_ALP (default 0.9 / 2.0): E-W coefficient
  and exponent.
- KN_FRC_TMIN / KN_FRC_TMAX (default 1E-4 / 1E-2): tolerance band.
- The I2 cap is derived, not a key: 0.1*GS_MATURE_TOL.

## Safeguards

- First armed call has no ratio: seeds FPRV only, tolerance stays
  at KN_RAY_TOL (the WP-1 default).
- ||dr||-collapse analogue: a non-positive or absent new residual
  leaves the update untouched (no division by a dead demand).
- Safeguard floor (GAM*eta_{k-1}^ALP > 0.1) prevents the tolerance
  reopening after a contraction stall -- E-W's own protection.
- eta_max 0.9, TMIN/TMAX clamps, and the I2 cap all bound the
  price; the tolerance can never loosen past a tenth of the
  population maturity gate.

## Verification (all actual outputs)

- test_knop_frc.py: 35/0. seed holds base 2E-3; seq
  1.0->0.5->0.25->0.1 gives eta 0.225,0.225,0.144 and tol
  2.5E-3,2.5E-3,1.6E-3 (monotone tightening, no inversion); rise
  1.0->1.0 caps eta at 0.9 and the tolerance lands on the I2 cap
  5E-3, not TMAX 1E-2; sflo fires the safeguard (eta_prev=0.6,
  floor 0.324 beats the raw 0.081, SG=1); tmin binds at 1E-4;
  noup leaves the base; GET_FG_FRC exports eta/rat/sg/tol
  verbatim.
- test_iter_ledger.py: 33/0 with the v4 columns (row parse,
  armed row 2.5E-3/0.5/1/2.5E-3, disarmed row -1/0/0/2E-3).
- test_vadat_keys.py: 134/0 (+5): defaults, arming, knobs, the
  staged I2 cap.
- Flag-off regressions: test_knop_ordering 50/0,
  test_knop_commit_gate 18/0, test_iter_conv_gate 18/0,
  test_conv_exit 53/0, test_knop_osgn 14/0 -- all identical, as
  required (KN_RAY_TOL_EFF=KN_RAY_TOL when disarmed).
- test_knop_ray 50/1: the csut certification table (T39) caught
  the unarmed-sentinel hole on the direct KN_CERTIFY path and
  now passes with the KN_CERT_CH local resolve; T38b's 2.2e-19
  OMP-race store noise from WP-5 returned to exact zero under
  this rebuild (codegen scheduling, as documented). T18 remains
  the documented pre-existing red.
- Full build.sh suite: stale-archive gate tripped mid-run when
  a parallel session touched jfnk_solver_v1.f90; all ~20
  freshness-gated tests reran green after the member timestamp
  was refreshed (no recompile of the parallel WIP). Remaining
  suite reds are the documented pre-existing set plus the
  parallel session's uncommitted SOB files.

## Flag-off identity

Every tolerance read inside the refresh chain switches from
KN_RAY_TOL to KN_RAY_TOL_EFF, and the refresh head resets
KN_RAY_TOL_EFF=KN_RAY_TOL unconditionally. With the flag off the
arithmetic is identical; the only other differences are module
state writes (FNEW recorded, EFF assigned) that feed no output.

## Owed measurement

Armed arm on the m15 stall model (KN_FORCING_TOL=T): watch
frc_eta/frc_tol in ITER_LEDGER track the demand fall -- tolerance
should start loose (~5E-3 cap while the demand is high), tighten
toward 1E-4 as MAXCH drops, and never exceed 5E-3 (I2). Compare
iteration cost against the fixed 2E-3 baseline.
