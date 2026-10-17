# KNOP Formal Solver Audit - Static Analysis Items (Oct-2026)

Scope: `plane/fg_j_cmf_nonmono_v1.f` and the KNOP control/test
machinery. 22 items from an external static analysis. Each item is
FIXED, VERIFIED (guard already present, confirmed correct),
REFUTED (the reported issue does not exist in this tree), or
DEFERRED (documented, not changed).

## Array bounds

1. `COMP_DEPTH_IDX(K)`/`K+1` with `K = RAY_PNT(I,LS)-ND_ADD` -
   VERIFIED. Both bracket indices pass explicit `K .GE. 1 .AND.
   K+1 .LE. ND` tests before use (two sites, I_M and I_P).
   `KC1`/`KC2` are clamped copies, never the raw K.
2. `I_M_DEPTH_NEXT(KC1,LS,FREQ_IDX_ITER+NM_RAY(I))` - VERIFIED.
   The frequency slot is tested `>= 1 .AND. <= NCF_IN` at both use
   sites, and `NM_RAY` is capped at `NCF_IN-FREQ_IDX_ITER-1` by the
   matched-span logic. Covered by unit-level reasoning only; no
   harness exists to drive the solver standalone.

## Division by zero

3. `Q = GAM_MAG/DLN_C` - FIXED (eb1638ac). `DLN_C` is floored
   positive; a degenerate value yields a locally opaque cell, the
   I -> I_NEXT advection limit that a zero log-spacing physically
   means.
4. `dGAMdR` central difference - FIXED (eb1638ac). Zero-spacing
   rays (J=K after clamping) skip the difference term; the
   sigma-scaled and curvature terms remain.
   Related hardening in the same commit: the `R(K+1)-R(K)` I_NEXT
   interpolation weight now uses a relative-tolerance test instead
   of exact `.NE.`.

## Uninitialized variables

5. `W_INEXT`/`INEXT_V` - REFUTED. Both are assigned on every branch
   of the guarded block before any read; the OMP clause already
   lists them PRIVATE.
6. `T1_DTAU` / `dS(NI)` - FIXED (e2516057). NI=2 rays indexed
   `DTAU(0,LS)`; now falls back to constant-slope extrapolation
   for NI<3.

## Memory management

7. JNU_STORE reallocation - FIXED (2a1b3a3f): STAT= + loud STOP.
8. ALO block reallocation - FIXED (2a1b3a3f): STAT= + loud STOP.
9. Grid-rebuild deallocations - REFUTED as silent-failure risk. A
   DEALLOCATE without STAT= aborts loudly at runtime on failure;
   the failure mode the analysis fears (silent undefined state then
   cryptic allocate failure) cannot occur. Adding 50 STAT= blocks
   would enlarge the diff without changing behavior.

## Optimization

10-12. GAM/GAM_MAG/dGAMdR/V_RAY/SIGMA_RAY caching - REFUTED.
    Already cached: all are computed inside the
    `FIRST_TIME .OR. NEW_R_GRID` block and persist across
    frequencies; the per-frequency path takes NEW_FREQ only.
    Invalidation keys off the R-grid, V, and SIGMA_TRUE
    fingerprints.
13-14. `AV_STORE=AV`, `CV_STORE=CV`, `I_P_PREV=I_P_STORE` copies -
    REFUTED. These are the lagged two-buffer scheme (prev vs
    current frequency); the copies are the mechanism, not waste.
15. `ALO_CHIHAT_RAY=CHI_RAY` per-ray copy - DEFERRED. NI elements
    per ray; measured cost negligible vs the formal integration.
16. Freq-match loop short-circuit - VERIFIED/REFUTED. The
    DO I=1,NI loop must run to set `NM_RAY`/`DLN_RAY` defaults;
    the `GAM < 0` test inside already skips all real work for
    expanding cells.
17. Fingerprint loop combination - REFUTED. V/SIGMA share one
    pass already; the R-grid loop exits on first mismatch and its
    distinct diagnostic message is worth keeping.

## OpenMP

18. PRIVATE audit - FIXED (80067173): `J` (loop index written
    inside the region) added to PRIVATE. `N_NEG_INT` updates use
    `!$OMP ATOMIC`; first-failure writes sit inside `CRITICAL`.
    No other shared writes found in the parallel sweep region.

## Tests

19. KNOP core behaviors - PARTIAL. The solver is monolithic
    (single routine needing full model state); standalone unit
    drivers are not feasible without restructuring, which violates
    the minimal-diff directive. Coverage achieved:
    - VADAT/parse path: `unit_tests/test_vadat_keys.py` exercises
      the KNOP flag contract, auto-sets, and new JFNK keys against
      the production parser (42 checks).
    - Behavioral items (chi_hat positivity, sign-dependent upwind,
      GAM sign in compression, read-before-write, refresh
      convergence, INEXT_CAP, hysteresis) are exercised by the
      running model `m15_2p82d_KNOP_w19` (jobs 12255111/12257454),
      whose per-iteration diagnostics report each mechanism's
      counters.
20. End-to-end integration test - satisfied by the model run
    itself (DO_NONMONO_KNOP=T through `run_model_grid.py`); see
    AGENTS.md verification checklist.

## Documentation

21. `FREQ_IDX_ITER` convention - VERIFIED. Comments at the KN_*
    capture site and the matched-span block already state the
    0-based iterator vs 1-based store-slot (slot s = frequency
    s+1) convention.
22. sigma=-0.999 beta-law floor - VERIFIED. Comment block at the
    THK extension already explains: BETA < 0 under shock sigma(1)
    would diverge V_EXT; flooring the exponent source at -0.999
    limits the extension to nearly constant V and warns.

## Commits (this audit)

- 2a1b3a3f dealloc STAT= checks
- eb1638ac denominator guards
- 80067173 OMP private J
- e2516057 NI<3 end-segment slope
- 79a00d5f VADAT key tests
- (earlier) 15797e32 sweep-metric denominator floor,
  478155d1 sweep report formats, 4d2259e7 JFNK_FADAPT optional key
