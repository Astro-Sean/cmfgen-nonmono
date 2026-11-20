# KNOP lagged-map non-contraction audit (Phase 2.6)

Status: source-level analysis complete; controlled runtime
experiments N1-N5 defined but not executed. Baseline commits
7ea19241, c9229c35. Companion data: knop_maturity_k1k4.json
(Phase 2.5), docs/audits/knop_laggedmap_sources.json (this
audit).

## 1. What the iteration actually is (source trace)

Module: `plane/fg_j_cmf_nonmono_v1.f` (14,541 lines). Two
nested fixed-point processes run inside each sweep:

A. In-sweep write (per frequency, inside the production ray
   solve): at frequency slot i (0-based FREQ_IDX_ITER) the ray
   solve reads the compressing-direction donor I_NEXT[i+1]
   from the PREVIOUS generation and writes the freshly solved
   field into I_NEXT[i], under-relaxed against a per-frequency
   snapshot (SNP_P/SNP_M, taken once per NEW_FREQ):

     I_NEXT(KC,LS,i) <- T2 + INEXT_RELAX*(T1 - T2)
       (fg_j_cmf_nonmono_v1.f:5837, :5891)

   where T1 is the new solve, T2 the stored value (or the
   snapshot on DO_ES_ALI re-entries, :5718, header :262).
   Compressing cells (GAM<0, only possible where
   SIGMA_TRUE<-1, :237) use this lagged donor; expanding cells
   use the usual I_PREV previous-frequency field.

B. Refresh (FG_KNOP_REFRESH, runs once at end of a complete
   sweep, FREQ_IDX_ITER==KN_NCF-1, :8125-8135; gated on
   DO_KNOP_REFRESH, FREEZE_I_NEXT=F, HAS_COMPRESSING,
   SOLUTION_METHOD='INTEGRAL'): rebuilds monotone-cubic
   CHI/ETA coefficients on the extended grid (frozen across
   refresh calls, KN_FREEZE_COEF, :8149-8193), then solves
   each ray's full wavelength chain
   (KNOP_RAY_SOLVE_V1, one ray per thread, :8563-8571) as a
   MULTI-PASS fixed-point iteration (passes<=61 observed),
   with exit classes exact/tol/bulk/stall/cap/div/plat and
   per-cell commit policy.

So the implemented map is a doubly-lagged iteration:

   outer generation g:   I^(g) -> [frequency sweep,
      relaxed writes, lagged donors] -> I_sweep^(g)
      -> [per-ray inner fixed-point iteration, truncated at
      maxpass or a stop rule] -> I^(g+1)

Neither stage is an exact solve of the coupled system: the
sweep write is damped (INEXT_RELAX in [0.05,1]) and the inner
ray iteration is truncated before tolerance for the ray
classes that exit via stall/cap/plateau.

## 2. Residual definitions (from source, not names)

- `refresh residual` (KN_RF_RESID -> FG_GS_MAXCHG):
  MAXVAL(RES_LS) over rays, RES_LS = max unrelaxed
  |vnew-vold| / max(vold, 1e-4*column peak) of the inner
  field solve's last pass (:8811-8816, header :369-373).
  It is the LAST-PASS iterate change of the inner map, not a
  residual of the coupled system.
- `independent q residual` (KN_Q_RESID, diagnostic):
  "max over rays of the one-pass map application on each
  channel's converged field" (:9011-9015) - i.e.
  ||G(I_committed) - I_committed|| normalized per cell
  (vold denom observed: |8.68e-6 - 3.53e-6|/3.53e-6 = 1.46).
  It is a true fixed-point residual of the refresh map.
- Difference explained: refresh resid measures the inner
  solver's own stopping increment; q measures what one more
  map application would change. q ~ O(1) with refresh
  resid ~0.8 means the committed store is a TRUNCATED
  iterate, far from a fixed point - confirmed by the exit
  census (35 bulk + 1 stall + 15 plat exits of 85 rays on
  K2 refresh 1; "unconverged cells=49.1M of 728.4M, 6.7%").

## 3. The adaptive omega controller measures the spectral radius

KN_ADAPT_RELAX (default T) sets INEXT_RELAX each refresh
(:8824-8870). Under FG_ADAPT_OMEGA_V2 (default T), the input
is KN_OMEGA_POOL_V1 over RHO_M_LS/RHO_P_LS - "each field
solve's own median post-transient same-direction pass ratio
(the inner map's asymptotic rate)", pooled at the 90th
percentile (comment :8837-8852). Then

   omega <- 1/(1 + rho_p90)        (:8852-8866)

so the reported `adaptive INEXT_RELAX: omega= 5.072E-01`
implies rho_p90 ~ 0.97 for the inner map. The code's own
telemetry therefore reports the dominant-mode contraction
factor of the wavelength-coupled ray iteration as ~0.97 per
inner pass - a near-identity map. The outer lagged map
inherits this stiffness: observed outer decay ~0.96/refresh
is exactly a damped rho~0.97 iteration
(0.5 + 0.5*0.97 ~ 0.98, consistent with measured 0.977->0.757
drift over six refreshes).

## 4. Root-cause assessment (Task 4 matrix)

| hypothesis | verdict | evidence |
|---|---|---|
| A. GS ordering wrong | excluded | compressing donor is next freq by construction; slot read/write order verified (:242-246) |
| B. I_PREV/I_NEXT exchange wrong | excluded | in-place slot scheme correct for a one-generation lag; SNP snapshot prevents re-entry compounding (:254-268) |
| C. stale coefficients | excluded as primary | KN_FREEZE_COEF pins chi/eta between refreshes; fingerprint/hash guards report drift (:8165-8193) |
| D. freq-chain init wrong | not shown | cold vs warm start changes the transient (fired 2,337 vs 4,680) but both plateau at the same map class |
| E. GAM sign wrong | excluded | sign(GAM) selects donor direction; verified against Paper I form (:75-85); index-62 result confirms the sign convention produces the physically right suppression |
| F. boundary refresh inconsistent | not excluded | boundary contributions enter the ray solve each pass; no separate lag detected, but not fully traced |
| G. residual evaluated on wrong field | excluded | q resid is computed on the committed store, independent path (:9011+) - the O(1) value is real, not a metric artifact |
| H. spectral radius >= 1 / stiff map | CONFIRMED as dominant mechanism | code-measured rho_p90 ~0.97; chi_hat = chi + |GAM|/DLN with |GAM|/DLN >> chi makes the wavelength chain a near-pure translation (spectral radius -> 1 as coupling dominates); exits stall/plat confirm inner non-convergence |
| I. wrong diagonal/off-diagonal split | consistent with H | the split itself is mathematically correct (verified Phase 2.5); the STIFFNESS it creates is the problem - donor-dominated cells carry almost no local information |
| J. source-function feedback inconsistency | contributing | SOURCE=zeta+omega*J with omega~0.999 at depths 66-70 adds a second near-unity feedback to the composite map |

Root cause (dominant): the lagged map is non-contracting
because the generalized-opacity wavelength chain is stiff.
With chi_hat ~ |GAM|/DLN, each cell's discrete equation is
I_i ~ I_{donor} + O(chi*DLN/|GAM|) - an almost-exact
translation along frequency. Any lagged solver (GS, Jacobi,
damped) then has spectral radius ~ 1 - O(chi*DLN/|GAM|),
which the code itself measures as rho_p90 ~ 0.97. Two
compounding factors: (i) inner solves exit via
stall/cap/plateau above tolerance and their truncated
iterates are committed, so q resid stays O(1); (ii) the
adaptive omega controller is a feedback loop on the residual
sequence - it can sustain a bounded limit cycle and was
already flagged in-source (W5 comment :8830-8836) for
misreading a converging field as rho~1.

This is a numerical-stiffness diagnosis, not an equation
defect: the Phase 2.5 audit's equation stands. No transfer
equation change is implied or recommended.

## 5. Correspondence with observed K1-K4 behaviour

- Outer residual ~0.96/refresh: predicted by rho~0.97
  damped at omega~0.5.
- Independent q oscillating 1.46,3.35,1.04,5.95,5.63,0.99:
  truncated committed iterates re-evaluated by the verify
  map; the oscillation is the adaptive-omega feedback plus
  wandering argmax cell (d62 -> d14 -> d1 -> d63 in
  ES_RESID).
- Rescue counts growing in the warm chain: the field is
  still migrating at ~50% per sweep (ES_RESID 0.49-0.66);
  more cells cross the absolute cap as the boundary region
  redistributes, consistent with a slow near-identity map
  rather than an instability.
- OBSFLUX divergence between runs: same cause - snapshots
  of a slowly translating field, not separate solutions.

## 6. What this means for convergence strategy

Increasing NUM_ES cannot fix a rho~0.97 map (it only moves
the truncation point). The correct lever class, in order:

1. Converge the inner per-ray system to a real tolerance
   before commit (the chain is a linear bi-diagonal-ish
   system in frequency per ray point - a direct or
   Krylov-accelerated solve removes the truncation entirely
   instead of iterating at rho~0.97).
2. Replace the ad-hoc adaptive omega with a fixed or
   asymptotically-ramped relaxation tied to the measured
   rho (the feedback controller is itself a suspect for the
   bounded oscillation).
3. Treat the scattering feedback (omega~0.999, J couples
   back into SOURCE) inside the same iteration level -
   currently the J->SOURCE->I loop straddles the in-sweep
   and refresh maps.

## 7. Defined experiments (not executed this session)

N1-N5 as specified in the phase brief, plus:
- P1: freeze KN_ADAPT_RELAX (VADAT INEXT_RELAX fixed at 1.0
  and 0.5) on the K2 setup to separate map stiffness from
  controller feedback - diagnostic only, never production.
- P2: report rho_p90 (already computed) alongside q resid
  per refresh to verify the controller's input.
- P3: per-depth/ per-frequency q-resid histograms (the
  KN_QC_* machinery already accumulates them; harvest from
  existing runs before writing new code).

All require the interactive/queue pipeline; none were run.

## 8. Machine-readable data

docs/audits/knop_laggedmap_sources.json - code locations,
measured rho estimate, residual definitions, exit census,
hypothesis matrix verdicts.
