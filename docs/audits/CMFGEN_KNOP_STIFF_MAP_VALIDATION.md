# CMFGEN KNOP stiff production-geometry packed-map validation (Phase 2.9)

Date: 2026-11-20. Case: `knray_p29 p1t1s`. Starting commits: d90d824b,
5b685c1d, ef6bef91; HEAD during this work: 24e227d6.

## 1. Provenance

- Git HEAD: `24e227d6` (user's KN_NOCP_PROBE MODE-1 series committed:
  84fa5cba, 336013fa, 6653ce34, 24e227d6).
- Tracked modifications: `plane/fg_j_cmf_nonmono_v1.f` one uncommitted
  line `FG_GS_NCONV=SUM(NST_LS,MASK=KN_PROD_MASK)` - a probe counter
  mask only; it does not touch the packed store map or its inputs.
  `unit_tests/f90/drv_knop_ray.f90` modified by this work (DRV_P1T1S).
- Untracked: `TIMING`, `fort.4`, `unit_tests/systems/_pv.txt` (scratch).
- KN_NOCP_PROBE machinery: present but inert; `p1t1s` never calls
  `SET_FG_NOCP_PROBE`, so the probe isolation path cannot run.
- Source hashes (sha256):
  - `unit_tests/f90/drv_knop_ray.f90` before this session:
    `90c9b9ff6f74deb14eec9fa26fabac84c59d9eff4f3ce4f1531e9d70a02dab6c`
  - `plane/fg_j_cmf_nonmono_v1.f`:
    `73ad55b8597b8cba62dee0d97f30afa0fe05a4eb1544ff242209f7b12463aecf`
- Build: `ifort -O0 -g -extend-source -qopenmp -mkl` for the driver;
  production objects via `unit_tests/common.py: prod_obj` for the
  helper sources (`test_knop_ray.py` HELPERS list).
- Executable: `unit_tests/bin/knray_p29`
  sha256 `f2583656ee4da2ec21ebbeabae988c177feb8ac5eebd4e47a377ddb91e85010e`.
  It is a separate binary; prior `knray`/`p1t1*` artifacts are not
  overwritten. The binary embeds the current (probe-era) solver
  source including the uncommitted mask line - results are therefore
  not a clean-commit artifact, and are labeled accordingly.

## 2. Production case and coefficients

Input deck: `docs/audits/knop_phase29/p1t1s_in.txt`, extracted from
`caseK2_knop_nores_es4` (the INEXT_RELAX~0.507 / pass-ratio~0.97 run).

- Longest contiguous `J_RESCUE_DUMP` trigger window: JDMP indices
  4238-4250, 13 consecutive frequencies, FL 10.8711042 -> 10.8497536
  (x1e15 Hz), dlog spacing ~1.668e-4.
- Directly observed per depth and frequency: R, V, SIG, SIGT, CHI,
  ETA, ZETA (ZETA=ETA/CHI verified to 5e-5).
- Scattering: ESEC and CHI_SCAT taken from the nearest OUT_FLUX
  per-frequency table (FL=14.3276). ESEC is essentially
  frequency-independent; the offset is documented.
- Geometry: ND=70, all 70 depths SIGT < -1; SIGT(62)=-155.61.
  CHI spans 2.24e-9 .. 0.517; ESEC 2.24e-9 .. 0.513.
- Assumed, not observed: VDOP = 1% of |V|; ZERO_FLUX inner boundary
  (production uses DIF=T, DBB~1.4e-4); single central ray P=0.
- Solver-side census after the cold sweep: 3224 compressing ray
  points, 282 frequency-matched offsets > 1, max offset 12 bins,
  NI=248. The packed state is NX=1680 (70 depths x 12 active slots
  x {P,M}).

## 3. Packed-map structure: the stiffness question is misplaced

The decisive finding is structural, not parametric.

Theorem (code-verified). In the packed lagged-store map F(X)=X, every
edge goes strictly upward in store slot:

- write: `I_P_DEPTH_NEXT(KC,LS,FREQ_IDX_ITER)` - slot = current freq.
- read (compressing): `I_P_DEPTH_NEXT(KC1/KC2,LS,FREQ_IDX_ITER+K2)`
  with `K2=NM_RAY>=1` (`fg_j_cmf_nonmono_v1.f:4691-4693`, capped at the
  phantom slot); expanding cells read `I_PREV` (not packed).
- hence output slot L <- input slot L+K2 >= L+1 always: the packed
  Jacobian is strictly upper-triangular in slot index, therefore
  nilpotent, in ANY geometry - production included.

Empirical confirmation on the stiff deck:

- swp_w10 (unrelaxed): iterating F from X0V gives residuals
  2.7e-11 -> 2.2e-15 -> exactly 0.0 at k=13=NSLOT. A nilpotent map of
  degree <= NSLOT confirmed on real coefficients.
- Graph check: all packed edges slot-increasing -> acyclic by
  construction; all strongly connected components are singletons.
  (The DBSI donor arrays are freq-shifted and unfilled on this path,
  so the edge count relies on the code-level proof above plus the
  exact-zero nilpotence certificate.)
- Power iteration on J (60 map matvecs): swp_w10 and fullw05 collapse
  to 0 immediately; swp_w05 converges to rho ~ 0.61 (tail range
  0.609-0.649), i.e. slightly above 1-omega=0.5: in-sweep relaxed
  writes mix old and new values at the same slot, so the relaxed map
  is not literally (1-w)I+wF; residual non-nilpotent structure is
  mild, nowhere near 0.97.

## 4. Stiffness hypotheses (Part 4)

- A (large |GAM|/DLN vs chi): exercised - SIGT=-155.6, offsets to 12
  bins, 3224 compressing cells. Still nilpotent. Ruled out for the
  packed map.
- B (cyclic cross-depth donors): impossible in the packed map - all
  edges are slot-increasing. Ruled out.
- C (expanding/compressing block coupling): expanding reads consume
  I_PREV, which is not packed; no packed edge crosses the blocks.
  Ruled out.
- D (scattering feedback): real ESEC/CHI_SCAT loaded; map unchanged
  in structure - scattering enters JNU/outer source updates, not the
  lagged store topology.
- Conclusion: the measured production rho~0.97 pass-ratio lives in
  the refresh's inner in-place iteration (the operator behind
  KN_Q_RESID), not in the packed lagged map. Phase 2.9's premise that
  the packed map might carry the stiffness is falsified; the stiff
  object is the refresh operator.

## 5. Direct-map Newton on the stiff deck

JFNK-Newton: outer Newton step, inner unrestarted GMRES(<=80, x3
restarts) on (I-J) via map matvecs (1 eval/matvec). No damping, no
shift; linear residual recorded.

- fullw05: 1-2 outer steps to <=2e-19 from all three seeds
  (warm / 0.5*X0+perturbation / zero). Seed agreement <= 2.0e-19.
- swp_w10: converges from all seeds; seed2 5 outer steps to 5.0e-15,
  seed3 3 steps to 2.3e-14. Seed agreement <= 3.1e-14.
- swp_w05: converges for warm seed in 1 step; cold seeds take
  3-6 outer steps, residual down to 8.7e-10 (seed2) / 2.0e-6 (seed3,
  Krylov-capped). |Xfix2-Xfix1|=2.1e-9; |Xfix3-Xfix1|=2.2e-6 - a real
  but small seed spread at the relaxed fixed point.
- GMRES iteration counts 10-80 consistent with a
  near-nilpotent+identity spectrum; no singular-matrix events.

## 6. vs original iteration / KN_Q_RESID

The direct solve collapses the packed-map residual from all seeds
(lagged fixed point reached in O(1) Newton steps, not O(30) sweeps).
This is a statement about the packed map only. KN_Q_RESID is the
refresh inner residual - a different operator; the map-level solve
cannot be claimed to collapse it without a refresh-side experiment.
Flagged as an open item rather than asserted.

## 7. Previous norm comparison, resolved

Recomputed on identical vectors (p1t1b_X0.txt vs Xstar per variant,
both absolute vectors of the same 30-element packed state):

- abs L2 = 2.141e-6, rel L2 = 1.42e-8
- abs Linf = 9.907e-7, rel Linf = 1.42e-8
- per-element relative max = 9.3e-8
- refresh truncation residual (same units) = 8.6e-9 relative

So the lagged fixed point differs from the committed field by ~1.4e-8
relative - 1.65x the refresh's own truncation floor. The old
"9.9e-7 vs 8.6e-9" compared an absolute Linf distance to a relative
residual (phantom 115x). Within-norm: consistent with refresh
tolerance, concentrated in P-channel slots 1-2; M channel <= 2e-8.

## 8. Deliverables and separation from Phase 2.8

- `docs/audits/knop_phase29/p1t1s_in.txt` - extraction deck
- `docs/audits/knop_phase29/p1t1s_extract.npz` - machine-readable arrays
- `docs/audits/knop_phase29/p1t1s_*.txt` - idx, donor dump, power
  iterations, Newton logs, stdout
- `docs/audits/knop_phase29/p1t1s_results.json` - machine-readable summary
- Production physics/defaults untouched; separate binary `knray_p29`.
- Phase 2.8 per-ray Schur validation remains a separate object; this
  experiment does not substitute for A,b construction, r=b-Ax, or the
  D3 same-equations test.

## Limitations

- Refresh inner operator (the actual rho~0.97 carrier) not yet
  probed; next step is a per-pass refresh map, not the packed map.
- ESEC from nearest cadence table; boundary/VDOP documented
  approximations; single central ray; 13-frequency window.
