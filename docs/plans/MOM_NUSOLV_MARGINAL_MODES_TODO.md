# MOM_NUSOLV marginal-mode fixes - implementation checklist

Status: phase 2 of the `MOM_NUSOLV` investigation (phase 1 verdict
recorded in `docs/CHECKLIST.md`, lamesF_nusolv arm, 24-Nov-2026).
Baseline: `7cd8e3c0`. F1-F3 landed and unit-verified 24-Nov-2026
(test_mom_nu 44/44 incl. marg/margok/cap; test_mom_solve 143/143
fixed-point regression untouched; mom_gs 34/34, vadat_keys 89/89,
aa_state 62/62, sweep_table 8/8). One design refinement during
implementation: the cap scale FSCL is the incoming column's own
scale, not the solved values - a solved explosion must not inflate
the trust region it is about to be capped by.

## 1. What the arm showed

`MOM_NUSOLV=T` on the lamesF deck (jobs 13651766 + 13652626, killed
at it675 after ~9 iterations):

- The functional-graph solve is mechanically correct in production:
  ~5.8K -> 2.3K rings/pass, 0 singular, 0 frozen on bad
  coefficients, negative rejections 20K -> ~460/pass.
- Every refresh still stall-exits; residuals oscillate 37 - 2.7e4
  with no contraction trend. The 37-250 minima are below the dead
  arm's ~2e4 plateau - the implicit-J map touches near-convergence
  intermittently but cannot stay.
- Committed stores escalated: after two vetoed iterations the
  commits ran 6.9e7 -> 4.6e11 -> 1.0e23 -> 1.8e8 largest change,
  up to 979K cells moved >5% - worse than the dead arm's bounded
  1e6-1e8 limit cycle. Populations dragged off the converged
  premise (outer luminosity 4.27e8 -> 8.9e7).

## 2. Root-cause hypothesis

The ring branch solves `x_j*(1-A) = S_j` exactly, where `A` is the
cyclic product of upwind coupling ratios around the ring. The
sign-flip pairs that pin the lagged map at rho~1 produce `A ~ 1`:
amplification `1/(1-A)` multiplies the *lagged* source terms
(lagged H, EPS scatter, other-depth J already in NU_X0), so the
exact solution of the stale-input system overshoots the true
composite fixed point by orders of magnitude. The 1e-10 freeze
threshold only catches exact degeneracy; `A ~ 0.9999` rings return
amplified solutions that pass the non-negativity guard, and once
the >50x-jump cell-count veto opens, a few cells can move by 1e23.

Two complementary bounds address this at the source; a third item
attacks the residual stiffness the J-graph solve cannot reach.

## 3. Fix checklist

### [x] F1 - near-marginal ring freeze band

Raise the ring-degeneracy freeze from `|1-A| < 1e-10` to a physical
band `|1-A| < MOM_NU_RING_BAND` (module constant, ~1e-2). A ring
amplifying its lagged sources by more than ~100x is solving noise,
not signal; freezing its members at the old iterate returns those
modes to the sweep's damped Jacobi update, while non-marginal rings
and chains still get the exact solve.

- Location: `plane/mom_j_cmf_nonmono_v1.f`, `MOM_NU_CHAIN_SOLVE`
  ring branch (~:7428). The `ST=2` resolved-with-old-value path is
  unchanged - downwind chain members still back-substitute from the
  frozen ring value.
- Census: count band-frozen rings separately from exactly-singular
  ones (`MOM_NU_NSING` stays the degenerate count; add
  `MOM_NU_NMARG` + `_T` accumulator, printed on the nu-solve line).
- Trade-off to document: a frozen marginal mode converges only as
  fast as Jacobi on it - if the entire stall cohort is marginal,
  F1 alone restores bounded behaviour without fixing contraction.
- No new VADAT key: compile-time constant pending arm evidence
  (AGENTS 9.1 - written proposal first if tuning is needed).

### [x] F2 - per-cell nu-write trust region

Bound the damped write `J_NEXT_FRQ(I,K0) = XOLD + RLX*(XNEW-XOLD)`
by a per-cell ratio cap: `XNEW <= MOM_NU_WRITECAP * max(|XOLD|,
1e-4*FSCL)` with `MOM_NU_WRITECAP ~ 1e2` (same order as the
retained-store 50x jump veto). Clamp (direction preserved, magnitude
bounded), not reject-to-old: a capped cell still moves toward the
solved value, and the next call continues the approach. Shrink
moves are unbounded - they cannot explode.

- Location: write loop ~:7470-7488, applied to `XNEW` before the
  change metric so N5/N50/CHG_NU report the actual write.
- Census: `MOM_NU_NCAP` + `_T`, on the nu-solve line.
- Catches chain members inheriting ring amplification via back-
  substitution, which F1 alone does not bound.

### [x] F3 - unit coverage for the new guards (same change)

- `drv_mom_nu.f` cases: ring with A = 0.9999 (frozen under band,
  members keep old iterate, NMARG counted); chain member fed by a
  large solved value (write capped at WRITEAP bound, NCAP counted).
- `test_mom_nu.py`: dense-solve reference asserts frozen members
  equal the seeded store; capped members equal the cap value.
- Re-run `test_mom_solve` nurf/nuw arms (fixed-point regression:
  guards inactive on comp1's contractive map must leave the
  RF_NUSOLV_* tolerances untouched).

### [x] F3b - iterated-dynamics probe (added during implementation)

The static cases cannot see the m15 failure - it was dynamical.
`drv_mom_nu iter` emulates the composite sweep+nu-solve loop on a
marginal 2-cycle in the explosive regime (lagged feedback
c/(1-a)=1000 per iteration) and measures the realised growth
rates under three guard configurations. Verdict (test_mom_nu
53/53): unguarded growth measured exactly 1000.00/iter - the
amplification rate c/(1-a), confirming the marginal-ring
mechanism quantitatively; cap-only growth 102.5/iter (the trust
region binding as designed, NCAP=6/iter); freeze-band drift
~1.04/iter (the damped Jacobi rate, NMARG=3/iter). Ordering
freeze < cap < unguarded matches the m15 phenomenology: stall,
bounded-rate explosion, explosion.

### [x] F4 - rebuild, relink, re-arm

ARM VERDICT (24-Nov-2026, job 13653300, killed at it675):
FAIL - the guards are correct but the marginal-ring
hypothesis is refuted. NMARG=0 throughout: no real-grid ring
is near-marginal, so the freeze band never engages. The cap
binds per-pass but 20-40 passes/refresh compound past it -
commits reached 7.5e23 with up to 1M cells moved, same class
as the unguarded arm, populations degrading identically.
The amplification is the exact solve's loop gain through
chain composition and lagged source magnitude - a property
of the lagged splitting itself in the explosive regime, not
of ring near-degeneracy or write magnitude. Phase-3 decision
(needs user check-in): implicit H channel to shrink the
lagged fraction, or close the line - the route stays
FAIL/opt-in otherwise.

- `compile_and_link_prod` + `make -C new_main` (the relink step
  that cost job 13651386 - do not skip).
- New arm `lamesF_nusolv2`: same eps1em3 seed, `MOM_NUSOLV=T`.
  Exit criteria as before: converges / stalls bounded (now
  without the amplification excursions - the guards' census
  numbers say which fired) / detonates.
- Watch: `NMARG`/`NCAP` vs `NRING` cohort sizes, residual trend,
  commit magnitudes vs the 1e23 excursions, population drift.

## 4. Phase-2 candidates (deferred pending F1-F4 evidence)

- **Implicit H channel.** `H_NEXT_FRQ` rows have the same
  functional-graph structure (`HS*H_up` lagged reads in the
  assembled RHS). Capturing the H coefficient + edge and solving
  H on the same graph removes the largest lagged term still
  explicit inside NU_X0. Larger surgery: second capture block +
  a second solve pass (or a coupled J/H block). If F1+F2 show the
  J-graph solve now contracts *except* on H-dominated cells, this
  becomes the next lever.
- **EPS scatter implicitness.** Same argument for the EPS_PREV
  midpoint transfers folded into NU_X0.
- **VADAT keys for WRITEAP/BAND** if the arm shows the mechanism
  works but needs tuning - written proposal per AGENTS 9.1.

## 5. Out of scope

- `DO_NONMONO_KNOP` remains reference/diagnostic (directive
  11-Nov-2026); this work serves the SOB + `USE_LAM_ES=F`
  (`MOM_J_CMF_NONMONO_V1`) route, currently FAIL/opt-in.
- No default-path change: all edits stay inside the
  `MOM_NUSOLV`-gated blocks and routine.

## 6. Phase 3: lagged-field investigation (24-Nov-2026, user directive)

F4 verdict (arm 13653300): NMARG=0 everywhere, cap binding but
unable to contain - amplification is not ring inversion. Phase 3
investigates the lagged field itself.

### Mechanism analysis (pre-census)

Matrix amplification is impossible on the captured rows:
MOM_LAG_BND bounds lagged weight by diagonal excess, so every
row coupling a < 1 strictly; chain composition then decays, and
NMARG=0 bounds ring gain 1/(1-A) <= 100. Yet the write cap
(100x) fires hundreds-thousands per pass. Resolution must be in
the captured source X0, not the matrix:

1. X0 carries the lagged H-divergence + EPS scatter + upwind-J
   content. H is reconstructed algebraically from same-frequency
   J each sweep, so the H store scales with whatever field was
   committed - including inflated commits. Feedback loop closes
   through the store, not the matrix: commit inflates J -> H
   store inflates in proportion -> next X0 inflates -> solve
   output linear in X0 -> commit re-inflates. No a < 1 bound
   constrains this because the gain is in R1 = X0/(TB*|J|), the
   source-to-field ratio, which the lagbnd census cannot see.
2. xm_pos evidence (arm 13653300, 9 refreshes): thousands of
   rows per refresh with negative RHS, worst clipped 1.9e5-4.5e5
   in scaled units. Per the code comments these are compressing,
   absorption-dominated rows with ~zero local source - the rows
   where the operator split has no physical positive source at
   all. The static source cannot balance the lagged debit there
   no matter how the write is guarded.
3. lagbnd massively active (864K-1031K boundary-lag cells per
   pass): huge row populations sit at or past the M-matrix edge
   (excess ~ 0), i.e. the compressing grid is exactly where the
   diffusion-style closure is marginal. Physical reading: in
   compressing regions the CMF transport is advection-dominated,
   H is legitimately negative, and the source-positive moment
   hierarchy does not admit the solution within a lagged frame.
   This is the structural reason lambda-style iteration cannot
   contract there - consistent with the lamesF stall itself.

### Census instrumentation (committed fcde29ad)

Per-pass nu-census print decomposes the capture: per-row a-value
max, source inflation R1 = X0/(TB|J|) max, H-share of X0, EPS
share, negative-RHS row count. Discrimination logic:

- rmax small (<10) with caps firing -> ring gain at the band
  edge (A ~ 0.99) amplifies moderate sources; tighten band.
- rmax large (>50) -> source inflation confirmed; H/EPS share
  picks the lever: H-share high -> implicit-H capture (feasible:
  VB/VC row coefficients already assembled; but converges toward
  the coupled wavelength-depth solve); EPS-dominated -> EPS
  implicitness or tighter EPS bound.
- Either way the compressing-row negative-source population is
  the structural floor: rows whose only positive source is the
  clipped zero cannot carry a physical fixed point of the split.

### Probe arm

lamesF_nusolv3, CMFGEN 13654121 (+OBS 13654122), same eps1em3
seed, instrumented exe. Queued behind QOSGrpCpuLimit at writing.

### Optimization decision tree (after census)

1. If H dominates inflated X0: implicit-H capture is the only
   in-frame lever, but note it converges structurally to the
   Baron & Hauschildt coupled wavelength-depth block solve -
   at that point prefer the real block solve over a half-port.
2. If EPS dominates: EPS-implicit treatment is smaller surgery.
3. If neither (static-source/clip dominated): the lagged frame
   is wrong in compressing regions, full stop - close the route
   and concentrate on the coupled solve or USE_LAM_ES=T.
