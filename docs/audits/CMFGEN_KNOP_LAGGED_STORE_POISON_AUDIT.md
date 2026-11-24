# CMFGEN KNOP lagged-store poisoned-generation audit (22-Nov-2026)

Question (user request): is the lagged-field source in
`plane/fg_j_cmf_nonmono_v1.f` incorrect, particularly for fully
non-monotonic velocity fields and poorly converged initial stores?
Evidence plan: trace every store read/write/heal path, drive
`FG_KNOP_REFRESH_V1` on manufactured profiles through
`unit_tests/f90/drv_knop_ray.f90` on the production objects, and
compare with the literature map (subagent report, cited below).

## What the new battery measures

`unit_tests/python/test_knop_lagged_map.py` (new, 26 checks):
successive `FG_KNOP_REFRESH_V1` calls on a fixed captured problem
step the outer lagged map directly, so the committed-store
trajectory measures the map's contraction and the committed field's
seed dependence.

| case | store seed | GAM profile | refreshes |
|---|---|---|---|
| lagc | zero (cold) | one compressing band | 3 |
| lagw | warm smooth | same band | 3 |
| lagr | zero (cold) | two-band reversal | 3 |
| lags | zero (cold) | strong compression | 3 |
| sdab | warm, then garbage+NaN/Inf | one band | 2 |
| sdp | NaN/Inf/negative/huge | one band | 1 |
| qres | warm | one band, KN_Q_RESID=T | 1 |

## Results on the manufactured problem

- Cold and warm seeds converge to the identical committed field
  (final-store max rel diff 0.0); the per-ray solve reaches
  residual ~1e-15 and the store stops changing after refresh 1.
  The committed field is a fixed point of the refresh map
  (KN_Q_RESID independent q max = 0.0).
- Two-band reversal and strong-compression profiles also converge
  and stay finite; outer-map store drift is 0.00e+00 on all three
  profiles. The manufactured grid is too small to reproduce the
  production rho~0.97 stall (NCF=6 vs ~4e4, mild |Q|); this test
  pins the contract, not the production rate.
- Poisoned seed (sdp): the solve diverges honestly
  (CLS=KN_CS_DIV, VC=4), 24 poisoned slots are healed, nothing is
  committed, RF_ALL_CONV=F. Containment works for a diverged ray.

## Confirmed defect: veto-protected poisoned retained store

sdab refresh 2 seeds the retained store with garbage (finite
~1E3x oversize sign-sawtooth plus one NaN and one Inf slot). The
ray solve converges (residual 1.8E-16, certified KN_CS_FULL), but
FG_JUMP_MAX compares the candidate against the poisoned retained
generation and the commit is jump-rejected (VC=1, delta 1.8E308).
The retained store keeps the NaN and Inf slots and every wrong
finite value. On the next refresh the same thing happens:
SEED_FLD re-seeds the poison, the solve converges again, the veto
rejects again. The store has no recovery path:

- HEAL_STORE runs only on the diverged-ray branch
  (fg_j_cmf_nonmono_v1.f:11750).
- The seed path copies the retained values unchecked
  (SEED_FLD, fg_j_cmf_nonmono_v1.f:12064).
- The FG_RAY_SEEDED one-unguarded-commit waiver covers only rays
  that never committed; a SEEDED ray's veto has no retained-side
  sanity check.
- The all-ray flag still reads T (cert census full=1,
  rf_all_conv=T) while the committed store carries NaN.

The same failure class is documented in-source for the cold case
("the veto then rejects exactly the converged fill it was meant to
let through", module comment ~line 605); FG_RAY_SEEDED fixed the
cold door only.

Reachability in production: all three store write sites sanitize
(non-finite/non-positive -> 1E-30 floor, above cap ->
NONMONO_INEXT_CAP: sweep write ~line 6511, WRITE_STORE ~line
15435, WRITE_STORE_CAP), and RESTORE_FG_INEXT_GS rejects
non-finite/negative checkpoint data (~line 7965). So a poisoned
retained generation cannot form through the normal write or
documented restore paths. The residual exposures are a corrupted
or pre-validation checkpoint, a foreign-seed transplant carrying
garbage the fingerprint cannot see, or any future write path that
bypasses the same screens. Lower reach, higher trap severity: once
present, the poison is permanent and the solve side reports clean.

Secondary exposure worth noting: the veto has no retained-side
sanity check, so a legitimately-stale (finite but far-off)
committed generation can in principle be held by the same door if
the required correction exceeds FG_COMMIT_JUMP -- a converged
candidate 50x above the retained value is vetoed the same way the
garbage was.

## Literature context (subagent report)

Paper II (Knop, Hauschildt & Baron 2009b, A&A 501, 813) converges
the coupled (space x wavelength) BVP inside each formal solve to
1E-8 at frozen source function; it never uses a
one-bin-per-outer-iteration lag. The dev lagged store is an
implementation shortcut whose propagation scaling (~NCF great
iterations) the tree itself identified as the maturity bottleneck.
Lagging the donor is the explicit-source treatment Hauschildt &
Baron 2004 showed unstable at strong coupling and Ceulemans et al.
2025 (JQSRT 343, 109470) showed admits -1 < lambda < 0 oscillatory
modes when the per-cell drift exceeds one bin -- the regime
DO_NONMONO_FREQ_MATCH exists to fix. Ceulemans' store-and-reuse
boundary-value cache (blackbody cold seed, stored field
thereafter) is the published precedent for dev's seed machinery --
inside a solve that still converges.

## Verdict (22-Nov-2026): PARTIAL + one confirmed defect

- VERIFIED on production objects: refresh committed fields are
  finite, seed-independent fixed points on manufactured
  single-band, reversal, and strong-compression profiles; cold,
  warm and perturbed seeds land on the same store; poisoned seeds
  diverge and heal without committing garbage; KN_Q_RESID reports
  an exact fixed-point residual.
- CONFIRMED DEFECT -> FIXED (24-Nov-2026): a converged certified
  solve could be permanently jump-vetoed by a poisoned retained
  generation; the store kept non-finite values indefinitely while
  cert reported FULL. Mechanism as diagnosed: SEED_FLD copies the
  retained slots verbatim, so a NaN/Inf slot survived the sweep at
  un-swept cells and the veto read the certified field's own
  non-finite candidate as a branch-1 HUGE return -- a
  veto-time-only heal could not rescue a candidate that had
  already inherited the poison.
  Fix landed on branch fix/convergence-liveness-and-guards:
  row-scoped retained-side sanity (KN_STORE_ROW_SANE +
  KN_ROW_MEDIAN/KN_HEAL_ROW sharing HEAL_STORE's median and
  zero-slot convention) now runs at the seed read, at the jump
  veto (heal, then one re-evaluation), and in the retained
  FG_COMMIT_NC=0 branch, counted under FG_POISON_HEAL. A sane
  store is untouched, so the default path is bit-identical.
  sdab now passes (VC=0, committed field identical to the
  refresh-1 seed-independent fixed point); sdac pins the genuine
  factor->50 veto with a sane store (VC=1, JUMP=1E6,
  FG_POISON_HEAL=0); sdad pins slotwise healing (one NaN per row,
  FG_POISON_HEAL=2*N_COMP, siblings preserved); OMP 1 vs 4 give
  identical committed fields and heal totals.
- NOT ADDRESSED: production-scale contraction (rho~0.97, period-2
  inner rays) needs the real model's ~4E4-bin grid; the
  manufactured battery cannot see it. No solver change made.
