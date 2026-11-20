# CMFGEN KNOP refresh spectral and branch audit (Phase 2.11)

Status: spectral question resolved; clean TEST-D run in fresh
processes; adaptive-controller mechanism identified. Outer/physical
comparison and dependency ablations remain open.

Provenance: HEAD 92f04e2e; binaries `knray_p211`/`knray_p211b`
built from `drv_knop_ray.f90` + production flags; deck
`p1t1s_in.txt` (Phase 2.9 production extraction, ND=70, NCF=13,
all depths compressing, NI=248, NX=1680). All experiments
diagnostic-only; each TEST-D case ran in a separate process so no
module state carried over.

## 1. Exact spectrum of the frozen refresh Jacobian

Dense finite-difference Jacobian (scale-aware eps per column,
elementwise eps = 1e-5*max(|Y_j|,1e-3)) on the full NX=1680 state
at the warm fixed point. All eigenvalues real.

| omega | rho(J) | spectrum | sigma_max | num. abscissa | nnz |
|-------|--------|----------|-----------|---------------|-----|
| 1.000 | **0.0** (nilpotent) | all zeros | 0.982 | 0.824 | 151444 |
| 0.500 | **0.5 exactly** | {0} + {0.5} | 1.168 | 0.979 | 132405 |

Structure: J_omega = (1-omega)*I + omega*K where K is nilpotent of
degree 12 (the number of free frequency slots). At omega=1 the
(1-omega)*VOLD self-term vanishes and the map is a pure compositional
chain -> J nilpotent -> ||J^12|| = 0 measured, and the cold-start
iterate terminates at exactly zero residual at pass ~12
(`p1t1r_hist_om100.txt`).

## 2. Reconciliation of the Phase 2.10 estimates

Both prior estimates were non-normal artifacts:

- Power iteration reported 0.78-0.79: for a Jordan-type matrix,
  |J^k v| decays like rho^k * polynomial(k); during the
  polynomial-transient phase the measured ratio exceeds rho. The
  run never reached the asymptotic regime.
- Arnoldi (40 steps) reported a 0.68 positive cluster at omega<1
  and a spurious -0.585 at omega=1.0: Ritz values of a strongly
  non-normal matrix are unconverged at k=40 << NX=1680 and sit
  inside the field of values, not on the spectrum. At omega=1 the
  true J is nilpotent; the -0.585 Ritz value measured FD/branch-
  switching noise, not an eigenvalue.
- Confirmed by dense eigendecomposition: rho(om100)=0,
  rho(om050)=0.5 exactly.

## 3. Non-normality

J is strongly non-normal at all omegas: sigma_max (0.98-1.17)
exceeds rho by a factor >2 at omega=0.5; numerical abscissa
lambda_max((J+J^T)/2) = 0.82-0.98. ||J^k|| decays through a long
transient plateau (om050: 1.17, 1.02, 0.92, 0.90, 0.87, 0.83,
0.80, 0.76, ... still 0.46 at k=15) instead of rho^k.

The fixed-omega iterate histories (`p1t1r_hist_om*.txt`, cold
start) show the measured step-residual ratio RISING through early
passes (0.33 -> 0.44 -> 0.48 at om075; 0.41 -> 0.73 by k~8 at
om507) before decaying: the observed pass ratio is governed by
the numerical abscissa and nilpotent transient, not rho(J).
Asymptotic ratios at k=50 were still above rho (~0.84 vs 0.5).

## 4. Branch census and the affinity defect

Branches in one pass (source): VNEW floor 1e-30, cap
NONMONO_INEXT_CAP, KN_COUPLING_OK gates, donor/frequency-match
selection, hot-window skips. Census at the warm state:
**12/1680 outputs sit exactly on the 1e-30 floor, 0 on the cap**.
The 4.9% affinity defect is concentrated on floor-branch elements
and their downstream coupling, i.e. piecewise-linear branch
switching, not smooth nonlinearity.

## 5. Clean TEST-D (separate process per case, KN_Q_RESID armed)

| case | seed | G(Y)-Y rel | npass | RESID | KN_Q_RESID |
|------|------|-----------|-------|-------|-----------|
| A normal | committed store | 1.4e-14 | 2 | 2.9e-12 | **0.0** |
| B fixed point | Newton Y* | ( = A: store already at Y* ) | - | - | - |
| C perturbed | Y* * 1.05 + 1e-6 sin | 9.0e-2 | 13 | 1.4e-2 | 1.1e-2 |
| cold | 0 | 1.0 | 14 | 8.3e-3 | 2.8e-3 |

The production deck's committed store IS a fixed point of G
(one-pass residual 2.2e-17) and KN_Q_RESID=0 there - warm
checkpoint already converged. From cold or perturbed starts the
production refresh stalls at RESID~1e-2 within ~13 passes with
qres~e-3-e-2, while the SAME operator iterated at fixed omega=1
converges to exactly zero in 12 passes.

## 6. The adaptive controller is the non-convergence mechanism

Fixed-omega histories vs the production run:

- omega=1.0: exact termination at pass ~12 (nilpotent chain).
- omega=0.75: fastest relaxed contraction, residual 3e-14 by
  pass 50 (asymptotic ratio -> ~0.48 still decaying).
- omega=0.507/0.5/0.25: geometric tails, ratios ~0.83-0.85 at
  k=50 (transient, trending toward 1-omega).

Production `SWEEP_FLD_ORD` adapts OMEGA_RAY at pass end:
RATIO_OM>1.05 (and NPO>=4) -> omega*0.7; RATIO_OM<0.5 -> omega*1.15;
sign-flip pairs -> omega*0.7 (lines 12215-12232). The omega=1
transient produces RISING residual ratios (0.11 -> 0.63), so the
controller damps omega off the nilpotent manifold within the first
passes; the relaxed map then carries the (1-omega) spectral cluster
plus nilpotent transient, producing exactly the ~e-3 stall /
~0.9+ observed pass ratio in production. The 0.97 pass-ratio is
therefore a controller+non-normality artifact, not rho(J).

## 7. Answers

1. Verified spectral radius: rho(J_omega)=1-omega (eigenvalues
   {0,1-omega}); rho(J_1)=0, nilpotent degree 12.
2. Arnoldi/power disagreement = non-normal non-convergence at
   40-60 Krylov steps on a 1680-dim Jordan-type operator; resolved
   by dense eigendecomposition.
3. Strongly non-normal: sigma_max/rho >= 2, num abscissa ~0.98,
   transient amplification visible in ||J^k|| and rising iterate
   ratios.
4. 4.9% affinity defect: 12/1680 elements on the 1e-30 floor
   branch plus coupling through those cells.
5. Direct fixed point valid: Y* equals the committed warm store;
   the floor-branch elements are active and handled by the
   wrapper's identical branch logic.
6. Converges from 3 independent seeds (Phase 2.10, unchanged).
7. Clean TEST-D: KN_Q_RESID=0 at the fixed point in a fresh
   process (case A/B identical since the store was already the
   fixed point). From cold/perturbed starts the original refresh
   stalls at qres~e-3 - which the fixed-omega=1 map does not.
8. The adaptive controller DOES explain the ~0.97 production
   ratio: it reads the nilpotent-transient ratio rise as
   divergence and damps omega, leaving the nilpotent manifold
   before the exact-termination pass count is reached.
9. Remaining outer limitation: coefficient refresh + population
   feedback untouched; inner solve cannot proxy them.
10. Phase 2.8 A,b/D3: not run; remains open and independent.

## 8. Caveats / open items

- Single-ray channel deck; the omega=1 nilpotence relies on every
  depth being compressing (SIGT<-1 everywhere). A mixed-sign deck
  may retain genuine eigenvalues at omega=1.
- ||J^k|| and the abscissa are dense-Jacobian quantities on this
  1680-state; production rays with more seed cells would need
  matrix-free equivalents.
- Ablation suite (spatial-only, wavelength-only, I_PREV/I_NEXT,
  scattering isolation) not yet implemented.
- eps-sensitivity of the dense Jacobian not swept; elementwise
  FD at 1e-5 relative.
