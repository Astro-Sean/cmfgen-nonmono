# CMFGEN KNOP structure-aware refresh controller (Phase 2.12)

Status: mechanism resolved; diagnostic strategies implemented
externally (per-pass omega policy over the production one-pass
map - zero production changes). Dense-graph nilpotency proof on
the all-compressing deck; mixed-sign deck tested.

Provenance: HEAD 4b3f90b2 + uncommitted driver modes; binary
`knray_p212`, production flags; decks `p1t1s_in.txt` (all-
compressing) and `p1t1s_mixed.txt` (depths 30-49 SIGT flipped to
+1, synthetic expanding region; NW drops 70 -> 52).

## 1. Controller audit (SWEEP_FLD_ORD pass end, ~lines 12200-12330)

- `RATIO_OM` = RES(pass k)/RES(pass k-1), per ray, unsigned
  residual norms - assumes monotonic decay.
- Omega update, unconditional at pass end:
  `RATIO_OM>1.05 & NPO>=4` -> omega*0.7 (floor 0.45);
  `RATIO_OM<0.5` -> omega*1.15 (cap 1.0);
  `NOSC>=2 & NPO>=4 & omega>0.45` -> omega*0.7 (period-2 branch;
  omega change rescales the AA/Aitken history).
- Exits: `RES<KN_RAY_TOL` -> IXIT=1; **bulk IXIT=2 when
  `NPO>=6 AND NB_O <= max(50, NT_O/1000)` AND
  `(KN_RAY_BULK_TOL<=0 OR RES<=KN_RAY_BULK_TOL)`**; stall IXIT=3
  (`RATIO_OM>KN_RAY_STALL` twice, NPO>=8); plateau W6; cap; div.
- `KN_RAY_BULK_TOL=0` is the production default (module DATA
  line 1131), so the residual floor on the bulk exit is
  unconstrained: a ray exits "converged" at any RES once bad
  cells drop under 50.

## 2. What actually stalls production on this deck

Cold-start production run (`p1t1d`, KN_Q_RESID armed, fresh
process): npass=14, RESID=8.3e-3, **exit class bulk=1**,
KN_Q_RESID=2.8e-3 (all-compressing AND mixed deck identical).

The same operator iterated at fixed omega=1 through the identical
one-pass production code path terminates at residual 0.0 at
pass 13 (nilpotent deck) / ~1e-14 at 15 (mixed). So the stall is
a control+exit artifact, not an operator property:

1. The bulk exit fires while the true fixed-point gap is ~e-3 -
   the refresh commits a field it then certifies as converged.
2. Within the first passes the omega=1 trajectory shows rising
   per-pass ratios and sign-flipping components (nilpotent
   transient); the NOSC>=2 and RATIO_OM>1.05 dampers can move
   OMEGA_RAY off the nilpotent manifold (floor 0.45 -> spectral
   cluster at 1-omega with a long non-normal transient tail).

## 3. Structure diagnostic (dense J, om100)

1680-dim Jacobian -> directed graph: **1680 singleton SCCs, pure
DAG, longest chain = 12** = free frequency slots. K^12=0 verified
(||J^12||=0) and the iterate terminates at exactly 0 at pass ~13.
Nilpotency proven from the graph, not from SIGT signs.

Mixed deck (18 expanding depths, NW=52, NX=1248): omega=1
collapses at pass 13 to 1.4e-11 (not exact 0) - classified
**stable but non-nilpotent** (small genuine eigenvalues or
branch-boundary effects). Still far faster than production's
stall at e-3.

## 4. Strategy comparison (driver-side pass loop over G)

All runs cold start, identical frozen coefficients; each pass =
one production `KNOP_RAY_SOLVE_V1(MAXPASS=1)` call.

| strategy | all-compressing | mixed |
|----------|-----------------|-------|
| A production adaptive | stall 8.3e-3 @14 (bulk exit) | stall 8.3e-3 @14 |
| B fixed om=1 | exact 0 @13 | ~1e-14 @15 |
| C window W=13 + damp | = B (never damped) | = B |
| D median-5 ratio ctrl | = B (no sustained growth) | = B |
| E structure-aware | = B | = B |

On both decks no strategy damps: the residual ratio never
sustains >1.05 beyond the transient window and the median
controller sees only contraction -> all equal the nilpotent
reference. Strategy B/E is optimal AND simplest.

## 5. Floor/branch handling

12/1680 outputs pinned at 1e-30 floor, count stable across all
passes and strategies (column NFLOOR in ctrl logs) - the floor
population does not churn, so it does not masquerade as smooth
residual evolution. Finite termination at omega=1 occurs despite
the pinned elements (they are fixed points of their own branch).

## 6. Answers

1. Is the controller misclassifying transient growth? YES - on
   this deck the bulk exit plus (secondary) omega dampers stop a
   nilpotent iteration at e-3 while certifying convergence; the
   fixed-omega reference reaches 0 exactly in the same pass count.
2. All-compressing nilpotent from the graph? YES - DAG, 1680
   singleton SCCs, chain length 12, ||J^12||=0, exact
   termination. Proven from dependencies, not SIGT signs.
3. Mixed-sign with genuine eigenvalues: the tested mixed deck is
   stable but non-nilpotent (no exact termination; ~1e-14 by 15).
   Within-row sign transitions and forward spatial reads are the
   remaining potential cycle sources - flagged, not exhaustively
   validated.
4. Transient-aware controller preserves convergence: yes - a
   window of ~chain-length (12-13) passes before ratio checks
   keeps omega=1 on the nilpotent manifold.
5. Structure-aware strategy on production deck: eliminates the
   e-3 stall; identical fixed point reached in ~13 passes.
6. Safe for cyclic configs: strategy E falls back to the median
   controller when nilpotency is not established; no unjustified
   omega=1 is applied.
7. Direct vs iterative fixed point: identical (Y* = committed
   store; both reached by fixed omega and by Newton).
8. Outer iteration: unchanged scope - inner fix does not claim
   outer convergence.
9. Production unchanged: all strategies are driver-external;
   no production source edits.
10. Recommended policy (diagnostic recommendation, not yet
    production): (a) run a dependency-graph/nilpotency check per
    ray before enabling adaptive damping - for DAG-structured
    unrelaxed operators hold omega=1 through a window of the
    measured chain length; (b) make the bulk exit conditional on
    a residual floor (KN_RAY_BULK_TOL>0 or the true verify
    residual), since the current bulk=1 exit can certify an
    unconverged field; (c) gate the NOSC/period-2 damper on
    sustained sign-flips rather than transient ones. Evidence:
    nilpotent proof + controller-exit trace + strategy table.

## 7. Caveats

- Single-ray/single-refresh deck; multi-ray sharing of statistics
  untested.
- Mixed deck is synthetic (SIGT flip); a real velocity-inversion
  model with genuine expanding/compressing coupling should be
  tested before production integration.
- The exact split between 'omega damper' vs 'bulk exit' blame
  within the production 14-pass run is inferred from exit class
  and structure; a per-pass omega trace inside production would
  nail it (not done - would need instrumentation).
- The store-commit relaxation (INEXT_RELAX=0.5) applies in both
  the external loop and production; it does not change the fixed
  point but enters the effective map between refreshes.
