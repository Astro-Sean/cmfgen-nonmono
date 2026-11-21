# Replay-vs-production discrepancy reconciliation (v2 decks, block 0)

All comparisons below are same-stage: deck PHASE-2 committed end
field vs replay KR_DBG_FLDP committed end field, and deck PHASE-1
FLD_M vs replay KR_DBG_FLDM. Deck = production live state captured
by the v2 writer. Replay = unit_tests/bin/knray_replay (-O0 objects
vs production -O3).

## A. Maximum absolute discrepancy, fields

| ray | array | channel | j | l | lam_A | prod | replay | abs diff | rel diff |
|-----|-------|---------|---|---|-------|------|--------|----------|----------|
| 082 | FLD_P end | P | 1 | 18496 | 998.8 | 6.229e-03 | 5.985e-03 | 2.442e-04 | 4.1% |
| 082 | FLD_M post | M | 1 | 18497 | 998.7 | 8.795e-03 | 8.578e-03 | 2.172e-04 | 2.5% |
| 081 | FLD_P end | P | - | - | - | - | - | 2.656e-04 | - |
| 084 | FLD_P end | P | - | - | - | - | - | 3.488e-08 | ~0 |

Normalization: plain element-wise difference on the committed
intensity arrays (KNOP comoving-frame intensity, solver units).
Source: deck records 30/31 (P end), 25 (M), vs replay
KN_REPLAY_END.bin records 2/3.

## B. The ~3.5% median M-channel discrepancy - resolved

Metric: median RELATIVE difference over bright cells
(|prod| > 1e-9) of FLD_M. Values: 3.67% (ray 082, bright
median), ~3.2% on the replay hst trajectory at 30 passes.
All-cell median is 1.3e-7 (dark cells agree to FP noise).

Cause: trajectory divergence of a NON-CONTRACTING fixed-point
map, not missing deck state and not a stage mix-up. Evidence:

- v1 and v2 replays return identical maxdiff (2.656368e-4 /
  2.442209e-4 to 6 digits) - EXT reconstruction and sibling
  zero-fill contribute nothing measurable.
- Matched pass count (maxp = prod np) gives the same maxdiff -
  not an iteration-stage offset.
- Production's own bright-cell drift is tiny: refresh-to-refresh
  M medrel 2e-4, P end->seed medrel 1e-4. The replay drift
  (~3.7%) is ~200x production's own drift - it is a different
  trajectory, not lagged data.
- Control ray 084 (contracting map, ixit=1): end-state agreement
  3.5e-8. Input completeness is proven by the control.
- Mechanism: production uses adaptive omega (FBITS bit 7 set);
  OUTGEN shows applied omega evolved to 0.514 during the capture
  while the deck stores only the base KN_RAY_OMEGA=0.45. The
  replay hst loop itself documents "minus the in-sweep omega
  drift" (drv_knop_replay.f90 ~line 365). In a non-contracting
  map each pass's omega difference accumulates as a coherent
  bright-cell drift; in a contracting map it is damped out.

## C. Starting state difference

Zero by construction: replay seeds FLD_P/FLDM from deck records
25/30, which are the production committed fields of that refresh.

## D. End state difference

Equals class A: P end maxabs 2.4e-4 (082) / 2.7e-4 (081) /
3.5e-8 (084). For non-contracting rays this is the expected
spread of two different FP builds on an iteration that does not
contract; it bounds trajectory fidelity, not data fidelity.

## E. v1 reconstruction effect

Nil at maxdiff precision: identical discrepancies on v1
(reconstructed EXT, zero siblings, no DBB) and v2 decks.
CHI_EXT/ETA_EXT are r4 in v2 (default REAL decl, source line 904),
so even v2 EXT is single-precision; node values from the v1
EXP(coef) reconstruction differ by < r4 epsilon.

## Residual comparison (not a field comparison)

| ray | prod res | prod np | prod ixit | replay res | replay np | replay ixit |
|-----|----------|---------|-----------|------------|-----------|-------------|
| 081 | 3.143e-02 | 12 | 3 | 3.279e-02 | 13 | 4 (maxp12) |
| 082 | 5.581e-02 | 24 | 3 | 6.219e-02 | 24 | 4 (maxp24) |
| 084 | 2.583e-04 | 3 | 1 | 1.665e-03 | 3 | 1 |

The stall/verdict behavior is faithfully reproduced; the residual
magnitude reproduces to ~5-10% for the non-contracting rays.
