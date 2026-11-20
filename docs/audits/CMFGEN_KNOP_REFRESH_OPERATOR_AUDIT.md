# CMFGEN KNOP refresh inner-operator audit (Phase 2.10)

Status: partial - operator identified, measured, and directly solved;
dependency ablations and physical/runtime comparison remain open.
All experiments are diagnostic-only; no production defaults changed.

Provenance: HEAD ~03e1a148 + uncommitted driver additions; binary
`unit_tests/bin/knray_p210` sha256
`eb6a300566faec8f68b1c42628a8042ecbccbc4faa7100b9d7bbd3a45e37ce0c`.
Input deck `docs/audits/knop_phase210/p1t1s_in.txt` = Phase 2.9
production extraction (JDMP freq window 4238-4250, ND=70, NCF=13,
all SIGT < -1, min SIGT -155.6, OUT_FLUX ESEC nearest-table).

## 1. Exact operator and state vector

Call chain: `FG_KNOP_REFRESH_V1` -> `KNOP_RAY_SOLVE_V1` ->
`SEED_FLD` + `SWEEP_FLD_ORD` (+ optional `KN_Q_RESID` verify pass).

State Y: the per-ray field `FLD(NI,NCF)` at the seeded cells -
for each compressing depth KC in the store, `SEED_FLD` places
`STORE(KC,LS,s)` into `FLD(J_PNT(KC),s+1)` for s=1..NCF-1. Free
state per channel: NW x (NCF-1). This deck: NW=70, NX=1680
(2 channels x 70 x 12). NI=248. Boundary/non-seeded cells are
parameters (frozen each pass), not state.

One pass (G): for each ray, for each direction - march J in ray
order; frequency loop direction set by `sign(GAM)` (compressing:
descending L). Compressing rows read `FLD(JPR, L+NNF2)` which is
ALREADY UPDATED this pass (true Gauss-Seidel in frequency);
expanding rows read `FLD(J,L-1)` raw (old iterate). Then
`FLD(J,L) = VOLD + OMEGA_RAY*(VNEW - VOLD)` in place.

So G is a nested composition: compressing outputs depend on
outputs computed earlier in the same pass; expanding outputs
depend on raw inputs. The non-nilpotent part of J_G lives exactly
in the expanding-row reads of raw inputs (the "-(D+L)^{-1}U" GS
structure). With OMEGA_RAY=1 the self-term (1-omega)*VOLD drops
out and the seed-point sub-block becomes compositional only.

G is invoked exactly once per pass via `KN_RAY_MAXPASS=1`; fields
exported through `KR_DBG_FLDM/FLDP` - zero production changes.

## 2. Classification: piecewise, near-affine

The row update is affine given a fixed branch set. The implemented
map is piecewise because of: VNEW floor (1e-30) and
`NONMONO_INEXT_CAP`, `KN_COUPLING_OK` gates, hot-window skips,
frequency-match donor selection, and pass-end OMEGA_RAY adaptation
(omega = f(residual history) - outside the one-pass map but inside
the production iterate).

Measured affinity defect at omega=1 (warm state):
`||G(Y1+Y2)-G(0) - [G(Y1)-G(0)] - [G(Y2)-G(0)]|| / scale`
= 4.93e-2 (`p1t1r_affine.txt`). Not globally affine; a single
global (I-M) solve is a Krylov approximation, not an exact inverse.

## 3. Jacobian measurements (matrix-free FD, eps scale-aware)

| omega | power-iter rho | Arnoldi top |lam| | top eigenvalue |
|-------|----------------|-------------|----------------|
| 1.000 | 0.0 (annihilated at k=12) | 0.585 | -0.585 |
| 0.500 | 0.787 | 0.685 | +0.685 |
| 0.507 | 0.778 | 0.681 | +0.681 |

Notes:
- omega=1: power iterate annihilated exactly at k=12 = number of
  free slots - the compositional/nilpotent chain in the frequency
  direction. Arnoldi still reports a -0.585 Ritz value, so J is
  not exactly nilpotent in FD arithmetic; the residual nonzero
  eigen-directions carry a sign-flip mode. Flagged unresolved
  rather than smoothed over.
- The 0.97 production pass-ratio is NOT the frozen single-pass
  Jacobian spectral radius (~0.68-0.79 here). It is an emergent
  property of the adaptive-omega multi-pass iterate plus outer-map
  dependencies - report rho_p90 and rho(J_G) as distinct numbers.
- Dominant eigenvalue at omega=1 is NEGATIVE (-0.585): the
  sign-flip mode the period-2 damping / NOSC machinery was built
  to catch. At omega<1 the dominant modes turn positive.

## 4. Direct solve

JFNK-Newton on Y=G(Y), GMRES(80) inner on (I-J), no shift/damping,
3 seeds (warm, 0.5x warm + 1e-3 perturbation, cold):

- omega=1.0: all seeds converge; fixed-point residual <= 4.8e-6
  absolute (state norm 1.54e-3 -> rel ~3e-3 scale-normalized);
  seed spread 6-9e-6 = Krylov tolerance, seed-independent.
- omega=0.5/0.507: same, residuals 1.6-2.5e-6.

## 5. TEST D caveat (seeded original refresh)

The warm store is already an exact G fixed point (one-pass
residual 2.2e-17). Re-running `KNOP_RAY_SOLVE_V1` seeded vs
unseeded gave 5.76e-3/npass=13 vs 2.9e-12/npass=2 - the first
call's carry-over state (RES_PREV, adapted OMEGA_RAY, committed
store) contaminates the second call, so this comparison is
confounded in the current harness. KN_Q_RESID collapse therefore
NOT established; needs per-call state reset or a fresh driver run
per seed. Honest FAIL on the comparison, not on the operator.

## 6. Answers (numbered per the brief)

1. State vector: seeded `FLD` at (depth-point, slot) per direction;
   equivalently the store slots feeding `SEED_FLD`. NX=1680 here.
2. Piecewise-affine, not globally affine (4.9% affinity defect).
3. Cycles: yes - expanding rows read raw inputs, compressing rows
   read fresh outputs; the composition creates genuine feedback
   (nonzero eigenvalues), unlike the packed lagged map.
4. Spectral radius near 0.97? NO - ~0.68-0.79 frozen single-pass.
   0.97 is a multi-pass adaptive-omega/outer property.
5. Rowwise coefficient ratio does not predict rho; never did.
6. Slow-convergence drivers: dominant negative eigenmode at
   omega=1 (sign-flip), near-translation stiffness in compressing
   chains; adaptivity converts this into the ~0.97 pass ratio.
7. Direct/Krylov exists and converges from all 3 seeds.
8. KN_Q_RESID collapse: NOT VERIFIED (harness carry-over, sect 5).
9. Outer problem: coefficient refresh + population coupling
   untouched - the inner solve converging does not imply outer
   convergence.
10. Phase 2.8 A,b/D3: not run here; remains open, independent.

## 7. Remaining gaps

- Dependency ablations (spatial-only, wavelength-only, I_PREV /
  I_NEXT isolation, clamp-on/off) not yet implemented in the
  driver - the flags exist in production but each needs a separate
  frozen-coefficient variant.
- TEST D clean re-run needs module-state reset between calls.
- Runtime/physical comparison (J, J_EQ, Eddington, rescue counts)
  not yet run on this deck.
- eps-sensitivity sweep for the Arnoldi spectrum not yet done.
