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

### [ ] F4 - rebuild, relink, re-arm

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
