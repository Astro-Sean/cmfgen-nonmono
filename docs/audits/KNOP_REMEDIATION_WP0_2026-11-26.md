# KNOP_REMEDIATION_WP0 -- ordering-stall reproduction fixture, 26-Nov-2026

WP-0 of the Knop remediation campaign: a manufactured ray that
reproduces the SWEEP_FLD_ORD stall on the current default
(KN_RAY_ORDERED=T) and converges under KN_RAY_ORDERED=F. No
production code changed; this is the regression gate for WP-1..WP-4.

## Fixture

`unit_tests/f90/drv_knop_ordering.f90` + `unit_tests/python/
test_knop_ordering.py`. The driver pokes `FG_J_CMF_NONMONO_MOD_V1`
globals and calls `KNOP_RAY_SOLVE_V1`, same pattern as
`drv_knop_block.f90`.

### Geometry -- how the stale donor is manufactured

NI=ND=140, NCF=200, J_PNT(K)=K, RAY_PNT(J)=J, KN_R(K)=141-K.
Every interior row sits at
`R_RAY(J) = KN_R(J) + 0.95*(KN_R(J+1)-KN_R(J))`, i.e. WINX=0.95
uniformly. The donor-geometry block in KNOP_RAY_SOLVE_V1 then gives,
for every compressing interior row J,

    donor = (1-0.95)*FLD(J,  L+NFM)  +  0.95*FLD(J+1,L+NFM)
             JPR1=J (own row)         JPR2=J+1 (deeper row)

In the ordered march (DO J=2,NI ascending, per-row L descending for
GAM<0) the own-row share is fresh (L+NFM done earlier in the row's
L loop) while the deeper-row share is a full pass stale -- the exact
defect the diagnosis assigns to inserted ray points at reversals.
The last row self-donates (K+1 off-grid) and anchors the chain;
depths 2..140 carry the sparse compressing slots (N_COMP=139) so
SEED_FLD seeds rows 2..NI at 1e-6.

This is a manufactured geometry -- the J_PNT-vs-R consistency that
the main routine asserts is bypassed deliberately, because the
audit model needs every row's donor to be the adjacent unmarched
row. Production's node-row anchors would otherwise collapse the
stall (verified in the Python replica).

### Coefficient regime

Uniform compressing sign (ALLNEG, ray-11-like), DRF/DLN_N=2.5 so
the matched span gives NFM=3 and DLNM=3*DLN_N=0.3. chi=g/999 gives
|Q|/(chi+|Q|)=0.999 and chi_hat*dz~0.83. ETA=chi_hat keeps SR=1.
Store seeded at 1e-6 (field is O(1)); KN_RAY_TOL=2e-3,
KN_RAY_MAXPASS=60.

## Measured outcomes (ifort -O0, this commit)

| case | config | NPASS | NXIT | RESID | RHO_M |
|------|--------|-------|------|-------|-------|
| ord  | ordered, matched, exits armed | 25 | 6 (plateau) | 4.121e-2 | 0.949 |
| uno  | same physics, FG_ORD_SWEEP=F  |  3 | 1 (clean)   | 4.054e-16 | - |
| ordx | ordered, all exits off        | 60 | 4 (cap)     | 1.588e-2 | 0.979 |
| unox | unordered, exits off          |  3 | 1           | 4.054e-16 | - |
| ordb | ordered, unmatched (drift 30, chi_hat*dz=30) | 17 | 6 | 8.565e-1 | 0.891 |
| unob | unordered, unmatched          |  3 | 1           | 1.991e-16 | - |
| ordm | ordered, mixed f=0.2/0.5/0.8  | 25/25/19 | 6/6/6 | 4.2/4.1/4.3e-2 | ~0.95 |
| unom | unordered, matched mixed      | 61/55/48 | 4/6/3 | 8.6e-3/5.7e-2/0.33 | ~0.95 |
| ordmb| ordered, unmatched mixed      | 13/13/23 | 1/1/2 | 1.7/1.4/2.1e-3 | ~0.23 |
| unomb| unordered, unmatched mixed    | 13/13/23 | 1/1/2 | 1.7/1.4/2.1e-3 | ~0.23 |

The ord arm lands inside the production stall signature (rays 81/82
on the real model: IXIT=3/6, res 2.9-8.8e-2, np 10-37). ordx hits
the pass cap at 1.588e-2 -- the same residual floor the independent
Python replica of this operator produced (1.59e-2). uno is the
block-triangular exact solve: 3 passes, residual at round-off.

## Deviations from the spec's reduced model, recorded

1. **DTAU=30 + NFM=3 is unreachable under matched spans.** MATCH_FREQ
   widens DLNM to NFM*DLN_N, so |Q|*chi_hat*dz = NFM_raw/NFM <= 1.
   Near-unity |Q| therefore forces thin cells (chi_hat*dz ~ 1) on the
   production path. The spec's drift-dominated algebra is carried by
   the ordb/unob arms instead: FG_FREQ_MATCH=F, drift=30 DLN bins,
   chi_hat*dz=30.03, |Q|=0.999, offset-1 donor -- same mechanism,
   honest inputs, stall reproduced (NXIT=6, res 0.86).
2. **The spec's "<=6 passes alternating" claim on mixed rays does not
   hold at |Q|~1 on both signs.** On the matched grid both marches
   non-contract (~0.95/pass): unom stalls at 48-61 passes, res
   8.6e-3..0.33. This is the lagged-map non-contraction the nbulk
   audit measured on the real model -- not an ordering artifact --
   and is asserted as non-clean in the test (T6) so it is recorded.
   On the unmatched grid (A0~0, spatially decoupled) expanding rows
   anchor the stale chain and BOTH marches converge (~13 passes);
   a pure-descending arm is not runnable because SWEEP_FLD's
   alternation is hardwired (ASC=MOD(NPO,2)==1).
3. NI sweep: the stall needs an unbroken compressing chain of
   length O(NI); the fixture fixes NI=140 (ray-11-like).

## Assertions in test_knop_ordering.py

- uno/unox/unob: NXIT=1, RESID<tol, NPASS<=8 (fast-side contract).
- ord/ordx/ordb/ordm: no clean tolerance exit (stall reproduced).
- ord: RESID>5*tol and RHO_M>0.9 (floor + marginal contraction).
- unom: non-clean on all fractions (recorded contradiction).
- ordmb/unomb: NXIT in {1,2}, RESID<3e-3, NPASS<=30 (broken chain
  converges).
- NMAX census: 3 matched, 1 unmatched.

Verdict: PASS. The fixture reproduces the stall on the current
default and converges under KN_RAY_ORDERED=F; it is the WP-1..WP-4
regression gate.
