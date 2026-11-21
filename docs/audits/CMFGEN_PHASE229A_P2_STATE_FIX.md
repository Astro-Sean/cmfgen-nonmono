# Phase 2.29A - KN_ORD_P2DAMP stale-state fix and validation

Date: 2026-12-01. Tree: `cmfgen_jfnk`, branch `main`, base `0348b1d6`.
Scope: the optional ordered-sweep period-2 damper in
`plane/fg_j_cmf_nonmono_v1.f`. `KN_ORD_P2DAMP` remains `F` by
default; no convergence criterion, tolerance, formal-solution
equation, or population-solver change.

## 1. Confirmed stale-state defect (CONFIRMED BUG, FIXED)

The pre-fix tracker ran only inside the significant-cell gate

```
IF(MAX(VOLD,ABS(VNEW)) .GT. 1.0E-4_LDP*VMAXF(L))THEN
```

while the midpoint commit ran unconditionally:

```
IF(DO_P2)THEN
  IF(P2_DMP(J,L) .GT. 0_1)THEN ... FLD=SQRT(VOLD*VNEW)
```

A cell that dropped below significance kept `P2_DMP` live: the
midpoint kept being applied to an unmonitored cell, and the cell
carried its old oscillation classification back across a
significance transition.

## 2. Fix shape

Triplicated inline tracker replaced by three shared module
procedures (now module scope for direct unit testing):

- `KN_P2_TRACK(SGN,FLIP,DMP,VOLD,VNEW,NFLIP,NREL)` - the full
  state transition, called inside the significance branch.
- `KN_P2_CLEAR(SGN,FLIP,DMP)` - total reset, called on the
  significance `ELSE` branch. Insignificant cell => ordinary
  relaxation, no midpoint, fresh history on re-entry.
- `KN_P2_MID(VOLD,VNEW)` - the damped commit midpoint.

Applied identically at all three sites: inward march cells
(`P2_*(J,L)`), the `I=NI` boundary row (`P2_*(NI,L)`), and the
outward march cells (`P2_*(J,L)`).

## 3. State machine (verified against the task's 10 checks)

- First signed update records `SGN` without counting a flip.
- `FLIP` counts consecutive reversals; `DMP=1` arms at `FLIP>=NFLIP`.
- Any non-reversal (same sign or zero delta) resets `FLIP` and
  advances `DMP`; `DMP>NREL` releases to 0. A reversal while damped
  re-arms `DMP=1`, so a persistent alternation cannot escape the
  bracket.
- `int*1` cannot overflow: `FLIP` capped at 100; `DMP` resets at
  `NREL+1<=101` since `FG_ORD_P2REL` is parse-clamped to <=100.
- Alloc failure sets `DO_P2=F` -> the ordinary solver path; state
  is allocated/zero-initialised per solve and deallocated on exit,
  so nothing survives between rays or threads.

## 4. Midpoint safety (CONFIRMED NUMERICAL LIMITATION, FIXED)

`SQRT(VOLD*VNEW)` underflows to 0 below ~1E-154 per endpoint;
stored field values reach ~1E-283 on floor cells (the ray-84
orbit of phase 2.28). Replaced by `SQRT(VOLD)*SQRT(VNEW)` -
mathematically identical for positive endpoints, finite over the
whole supported range. Overflow unreachable: `NONMONO_INEXT_CAP`
bounds inputs to 1E20. Not claimed as a general fixed-point
method: it is exact only for a pure multiplicative lambda=-1
mode.

## 5. Boundary-row diagnostic

`RES_INT` and `RES_BND` now accumulate alongside `RES` (inward
I=1 row and outward `NI` row are boundary; marched cells are
interior). `PRINT_KNOP_DEBUG` emits `KNPASS res_i/res_b/res` per
pass. Diagnostic only; `RES` remains the sole criterion and the
boundary commit semantics are unchanged.

## 6. Tests

New `unit_tests/f90/drv_knop_p2.f90` +
`unit_tests/python/test_knop_p2_state.py`: 19 checks, all PASS.
Covers activation (R1), release (R2), clear-on-insignificant
(R3), re-entry with clean history (R5), all three call sites on
the significance ELSE branch (R6/R7), re-arm under a persistent
alternation (R8/R9), alloc-failure path (R10), per-solve
deallocation (R13), and midpoint finiteness at 1E-283, 1E+200,
1E-200-pair underflow, and zero/negative endpoints (R14).
R11 off-state equivalence: the full `test_knop_ray.py` suite
passes 49/49 including the bit-identical repeatability tests
(T38) and the 1-vs-4-thread determinism tests (T8/T20, covering
R12). R4 is structural: `KN_P2_MID` is reachable only under
`P2_DMP>0`, which the fix now ties to the significance gate.

## 7. Replay evaluation (immutable `KN_REPLAY_HST_L08*.bin`)

The production KNRDUMP input decks no longer exist, so the
production operator could not be re-run with damping on. What
the recorded committed-field histories establish:

- All three rays end in a persistent pass-to-pass sign
  alternation; L081/L082 pin at interior cells near l~1517-1523
  (I<=2), not the far-IR boundary-row locus reported for a
  different pass window in phase 2.28 - the orbit migrates
  between low-amplitude cells.
- Orbit endpoints vs midpoints (L081 `I2L1520`):
  7.9330e-142 / 8.2573e-142 -> geom 8.0935e-142, arith
  8.0952e-142 (differ by 0.02%). L082 `I2L1517`: geom
  2.8941e-139 vs arith 2.8965e-139 (0.08%).
- Simulated tracker on the recorded delta sequences
  (P2N=2, P2REL=3): L081 cells activate at pass 9, release 12,
  re-activate 28; L082 activates at 25 and stays armed; L084
  activates at 3 on underflow-level cells. So the corrected
  damper would have engaged on the observed orbits.
- The unrelaxed map residual `R(X)` cannot be evaluated without
  the recorded operator; midpoint correctness cannot be inferred
  from the bracket alone. UNPROVEN as acceleration: pending
  P2-on replay or a fresh `KN_RDUMP_LS` capture.
- Normalization caveat: recomputing `rel` with
  `VMAXF(L)=max_I|FLD|` reproduces L081's ~4e-2 pinning and
  L082's ~8.5e-2, but L084 shows rel~1-6 at 1e-289 fields in
  this normalization although production certified it - the
  replay history predates certification semantics and needs the
  replay driver to reconcile.

## 8. Classification

- CONFIRMED BUG FIXED: stale `P2_DMP`/`P2_SGN`/`P2_FLIP` across
  significance transitions, all three march paths.
- CONFIRMED NUMERICAL LIMITATION FIXED: `SQRT(A*B)` underflow.
- UNPROVEN HYPOTHESIS: geometric midpoint as fixed-point
  acceleration on the production operator (no replay input
  decks; `R(X)` unevaluated).
- NON-BUG: off-state production path - bitwise unchanged
  (49/49 suite incl. bit-identity tests).

## 9. Recommended next experiment

Arm `KN_RDUMP_LS=81,82` (and the P2 flags on a diagnostic arm) on
the next refresh run to capture the decks this phase lacked;
then replay std vs `KN_ORD_P2DAMP=T` and evaluate `R(X)` on the
proposed midpoints through `KN_Q_RESID`.

Build: `make -C plane DEBUG=nodebug` clean, no warnings.
"""
