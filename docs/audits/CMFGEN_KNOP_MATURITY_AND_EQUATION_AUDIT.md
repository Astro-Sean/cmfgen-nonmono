# CMFGEN KNOP maturity and equation audit (Phase 2.5)

Baseline: Phase 2.4 report commit `72bac47b`.
Audit commits: `dd307f7d` (OBS FGINEXT_GS checkpoint dump),
`40842972` (obs-side KN_Q_RESID read),
`0a065e06` (JSCN cadence cap 150 -> 600 for 4-sweep dumps).

Frozen diagnostic executables:

| commit | exe sha256 (first 16) | frozen path |
|---|---|---|
| dd307f7d | dca8da4bd0102805 | exe_frozen/dd307f7d433f2130_dca8da4bd0102805/ |
| 40842972 | 6161a365ea2ca173 | exe_frozen/40842972_6161a365ea2ca173/ |
| 0a065e06 | b409293a14eaf3e0 | exe_frozen/0a065e06_b409293a14eaf3e0/ (K-matrix exe) |

Phase 2.4 ensemble: `obs_ensemble/` under the validation tree.
Phase 2.5 maturity runs (all on exe b409293a14eaf3e0, all with
KN_Q_RESID=T and ES_RESID_DIAG=T):

- `caseK1_knop_nores_es2`: NUM_ES=2 cold (current config, control)
- `caseK2_knop_nores_es4`: NUM_ES=4 cold (job 13344788)
- `caseK3_knop_nores_restart`: NUM_ES=2, warm-started from K2's
  FGINEXT_GS (job 13344889, afterok on K2)
- `caseK4_knop_nores_restart4`: NUM_ES=4, same warm start
  (job 13344890, afterok on K2)

All runs use identical frozen populations, RVTJ, MODEL and MODEL_SPEC
inputs. No production defaults, population numbers, RESC_JMAX or
J_RESCUE behavior were changed.

## 1. Executive summary

Phase 2.4 showed KNOP reducing the inner-boundary median field by
~100x relative to V12 but with an immature lagged I_NEXT field
(refresh residual 0.98 vs tol 2.0e-3). This audit establishes:

1. The published `4*a_lambda*I` term is genuinely absent from the
   implemented equation and has no mathematically equivalent
   contribution elsewhere in the ray solve. It is a documented,
   inherited CMFGEN truncation shared by stock V12/V13 - not a
   KNOP-specific defect. Expected effect: 3*dLOG_NU ~ 0.025% of
   the velocity coupling magnitude on the OBS continuum grid
   (dLOG_NU ~ 8e-5), sign-aware (absorption-side in expansion,
   source-side in compression). It does not explain the 40-100x
   V12/KNOP field difference.
2. The reported "I_NEXT residual=0.98" is a normalized lagged-map
   fixed-point iterate change (max over rays/cells of
   |VNEW-VOLD|/max(VOLD, 1e-4*column_peak)), not a transfer-equation
   residual and not a residual of the coupled wavelength system.
3. KNOP immaturity in Phase 2.4 is explained by: cold-started
   I_NEXT, NUM_ES=2 (one refresh), and stall/cap/plateau ray exits
   whose stores are not committed under FG_COMMIT_CONV.
4. The Phase 2.3 count of 75,738 V12 rescue events is the two-sweep
   event total (44,924 sweep-1 + 30,814 sweep-2), reconciled exactly
   against the Phase 2.4 mature-sweep count of 30,814.
5. The 1/(1-omega) ~ 435 scattering scale at depth 70 is a local
   equilibrium intensity for a source of form S = zeta + omega*J,
   not a hard bound; V12's tail sits at/above it, KNOP's ~10x below.
6. The K1-K4 chain (eight sweeps) is complete: iterate residual
   decays only ~0.96 per refresh while the independent
   fixed-point residual oscillates at O(1); rescue counts grow
   in the warm chain. The lagged map is bounded, not
   contracting - no converged KNOP spectrum exists.
7. Verified Phase 2.3 result carried in: the depth-62
   super-equilibrium spike (894 events, max 11x J_eq) is
   V12-specific. On the same 2,040 frequencies KNOP never puts
   its argmax at depth 62 and its global-max field is ~16x
   below V12's J(62) - the compression-cell spike is an artifact
   of the sigma clamp, not a physical field (section 9a).

## 2. Exact continuous transfer equation implemented

Published (Knop, Hauschildt & Baron 2009, Paper I, Eq 1):

    dI_lambda/ds = eta_lambda - chi_lambda*I_lambda
                   - 4*a_lambda*I_lambda
                   - a_lambda * d(lambda*I_lambda)/d(lambda)

with a_lambda the comoving wavelength-stretch rate (O(v/c):
a_lambda ~ mu^2 dbeta/dr + beta(1-mu^2)/r).

Implemented (source comments, fg_j_cmf_nonmono_v1.f ~line 4500 and
the V12/V13 construction it reduces to):

    dI_nu/ds = eta_nu - chi_nu*I_nu + GAM * dI_nu/d(ln nu)

with GAM = (V/(cR))(1 + SIGMA*mu^2) = -d(ln nu)/ds = +a_lambda at
O(beta) (verified: gamma^3*mu*(mu+beta) dbeta/dr + gamma*beta
(1-mu^2)/r -> mu^2 dbeta/dr + beta(1-mu^2)/r = GAM).

## 3. Exact discretized equation implemented

Decreasing-nu grid, dLOG_NU > 0 the local log spacing to the upwind
bin. Per cell i on ray LS:

    chi_hat = chi_i + |GAM_i|/dLOG_NU          (strictly positive)
    Q       = (GAM_i/dLOG_NU)/chi_hat          (signed)
    S'_i    = eta_i/chi_hat + Q * I_upwind_i

    I_M outward march (GAM >= 0): I_upwind = I_M_PREV  (bluer bin,
        already computed this sweep)
    I_M inward march (GAM < 0):  I_upwind = I_M_DEPTH_NEXT (redder
        bin, lagged I_NEXT store)

Formal step: short-characteristics integral with the monotonic
parabolic source fit (identical A0..A4/E0..E3 coefficient machinery
as V12). The equivalent difference form:

    dI/ds = eta - chi_hat*I + (GAM/dLOG_NU)*I_upwind
          = eta - chi*I - (|GAM|/dLOG_NU)*I + (GAM/dLOG_NU)*I_upwind

GAM > 0: -GAM*(I - I_prev)/dLOG_NU  = +GAM*dI/dlnnu upwind-diff
GAM < 0: +GAM*(I - I_next)/dLOG_NU  = +GAM*dI/dlnnu upwind-diff

so the discrete equation is the upwind-difference form of

    dI/ds = eta - chi*I + GAM*dI/d(ln nu)

i.e. section 2's continuous equation, differenced consistently.

## 4. The 4*a*I term analysis

Transformation (independent audit). The Lorentz-invariant specific
intensity is i = I_nu/nu^3. Along a comoving ray the invariant
equation is clean: di/ds = j - alpha*i, with the comoving frequency
drift d(nu)/ds = -a*nu (a = a_lambda):

    i = I_nu/nu^3
    di/ds = d(i)/ds|_nu + (d nu/ds) di/d nu
          = (1/nu^3) dI/ds|_nu - a*nu * d(I/nu^3)/d nu
          = (1/nu^3) [ dI/ds + 3*a*I - a*nu*dI/d nu ]

Equating di/ds = (eta - chi*I)/nu^3:

    dI_nu/ds = eta_nu - chi_nu*I_nu - 3*a*I_nu + a*dI_nu/d(ln nu)

Applying the identical construction to i = I_lambda*lambda^5/c^4
with d(lambda)/ds = +a*lambda:

    dI_lambda/ds = eta_lambda - chi_lambda*I_lambda
                   - 4*a_lambda*I_lambda
                   - a_lambda*d(lambda*I_lambda)/d(lambda)

which reproduces Paper I Eq 1 term-for-term. The two forms are the
same physics under i = I_nu/nu^3 = I_lambda*lambda^5/c^4.

Answers to the audit questions:

1. Explicitly implemented? No. No term proportional to I at
   coefficient 3*a or 4*a appears in chi_hat, Q, S', or the march.
2. Absorbed into a transformed variable? No - the field variable is
   the specific intensity itself; no invariant-intensity transform
   is applied (I_nu/nu^3 would carry no proportional term).
3. Absorbed into generalized opacity? Partially but not
   equivalently: |GAM|/dLOG_NU in chi_hat is the implicit half of
   the upwind *difference*, not the proportional term.
4. In discretized wavelength coefficients? No - the code is a
   frequency formulation; no p-/p|/p+ structure exists.
5. Canceled algebraically? Only in part: the code's own comment
   (fg_j_cmf_nonmono_v1.f ~line 4516) correctly identifies that the
   nu-form coefficient is 3 (not 4) and that an earlier S_tilde
   transplant (4*Q*dLOG_NU*I_upwind, removed 07-Sep-2026) was a
   double-count. But the -3*a*I_nu term has no counterpart in the
   solved equation.
6. Does omission change the physical equation? Yes. The complete
   equation is dI/ds = eta - chi*I - 3*a*I + a*dI/d(ln nu); the
   solved equation drops -3*a*I.

Expected numerical effect, measured on the actual OBS continuum
grid (74,704 frequencies, dLOG_NU median 8.3e-5, p95 1.7e-4):
relative to the implicit coupling opacity |GAM|/dLOG_NU, the
omitted proportional term is 3*dLOG_NU ~ 2.5e-4 (0.025%), and
relative to the local total opacity it is far smaller still
(3|GAM| ~ 1.5e-9 vs chi ~ 0.5-40 at depth 70). The omission is
sign-aware (absorption for GAM>0, source for GAM<0). On much
coarser frequency grids (dLOG_NU ~ 1e-2) it would reach ~3% of
the coupling. This is a documented inherited CMFGEN truncation
(the V12/V13 convention), orders of magnitude too small to
explain the 40-100x V12/KNOP field difference, which is set by
the upwind direction itself. Classified: known physics
deviation, not the artifact cause; a quantification run on a
coarse grid may be considered separately but is not warranted
on this grid.

## 5. Generalized-opacity derivation

Independent derivation from Paper I Eq 1 in frequency form. Write
the coupling C = -3*a*I + a*dI/d(ln nu). For the nu-grid decreasing
and information transported nu_high -> nu_low when a > 0
(expansion), the upwind difference is (I_k - I_{k-1})/dLOG_NU:

    C ~ -3*a*I_k + a*(I_k - I_{k-1})/dLOG_NU
      = -(3*a + a/dLOG_NU)*I_k + (a/dLOG_NU)*I_{k-1}

    chi_hat_correct = chi + 3*a + a/dLOG_NU
    chi_hat_code    = chi + a/dLOG_NU

a < 0 (compression): transport reverses (nu_low -> nu_high), donor
is I_{k+1} (lagged I_NEXT):

    C ~ -3*a*I_k + a*(I_{k+1} - I_k)/dLOG_NU
      = -3*a*I_k - a*I_k/dLOG_NU + a*I_{k+1}/dLOG_NU
      = -(|a|/dLOG_NU - 3|a|)*I_k - (|a|/dLOG_NU)*... 

    chi_hat_correct = chi + |a|/dLOG_NU - 3|a|   (the -3*a*I with
                        a<0 is +3|a|I, a source moved to the RHS)
    chi_hat_code    = chi + |a|/dLOG_NU

Positivity: chi_hat_code > 0 always (chi > 0, |GAM|/dLOG_NU >= 0).
The sign-aware source/absorption split of the missing term matches
the solver's own convention: for GAM>0 it belongs inside chi_hat
(sink), for GAM<0 it is a +3|GAM|I explicit source. Under the code
equation alone, chi_hat positivity does not alter the equation - the
|GAM|/dLOG_NU inside chi_hat is exactly the implicit upwind half of
GAM*(I_upwind - I)/dLOG_NU; moving it between sides is algebraically
exact for the solved (truncated) equation.

## 6. Wavelength-upwind audit

Verified against the solve loop (fg_j_cmf_nonmono_v1.f):

- OBS frequency grid strictly decreasing in nu.
- GAM(I,LS) = (V/cR)(1 + SIGMA_RAY*mu^2), signed.
- GAM == 0: no donor read, pure sink-free march.
- GAM > 0: donor I_M_PREV/I_P_PREV (previous, bluer bin).
- GAM < 0: donor I_M_DEPTH_NEXT/I_P_DEPTH_NEXT (redder bin,
  lagged store).
- DO_NONMONO_FREQ_MATCH=T (auto-set): the lagged donor foot is
  interpolated across the OBS/model grid difference (Ceulemans
  et al. 2025).
- Asymmetric-spacing correction: FG_NU_NEXT replaces dLOG_NU for
  compressing cells whose next-bin spacing differs (line edges).
- chi_hat = chi + |GAM|/dLN > 0 always; signed Q in source.
- Optical-depth increments: V12's DTAU construction on chi_hat
  (Euler-Maclaurin-corrected, clamped between endpoint chi*dZ).
- Source interpolation: monotonic parabolic dS fit, same as V12.
- Angular quadrature: unchanged ray set; mu enters via (1+sigma*mu^2).
- Boundaries: DIFFUSION inner condition on chi_hat rows;
  zero-flux/dbb outer handling as V12.
- Zero/sign-changing coupling: GAM==0 skips the donor; sign flip
  switches donor direction cell-by-cell along the ray.

## 7. GS iteration and residual definitions

Three nested iterations, reported independently:

1. Ray-local field solve (SWEEP_FLD / SWEEP_FLD_ORD): alternates
   frequency-march direction per pass; per-pass residual
   RES = max_cells |VNEW-VOLD|/max(VOLD, 1e-4*VMAX_L)
   where VMAX_L is the frequency-column peak on that ray. Cells
   below 1e-4 of the column peak are excluded from the census.
   Exit classes (IXIT): 0 exact single pass; 1 RES < KN_RAY_TOL
   (2e-3); 2 bulk (NPO>=6, bad-fraction <= max(50, NT/1000),
   contraction check, RES <= KN_RAY_BULK_TOL); 3 stall (two
   consecutive same-direction ratios > KN_RAY_STALL=0.98);
   4 pass cap (NPO >= KN_RAY_MAXPASS=60 or hot-pass budget);
   5 diverged (non-finite or RES > 1e3*best); 6 plateau (three
   consecutive same-direction ratios in [0.95,1.05]); 9 verify.
   Adaptive OMEGA_RAY: same-direction ratio > 1.05 -> *0.7 floored
   at 0.45; < 0.5 -> *1.15 toward 1; two signed-update sign flips
   damp likewise. Commit gate FG_COMMIT_CONV=T: IXIT 3/4/6 rays do
   NOT commit - their lagged store keeps the pre-refresh
   generation; FG_COMMIT_NC=0.05 bounds per-cell change for
   non-converged committed rays; FG_COMMIT_JUMP=50 bounds whole-ray
   rewrites.
2. Refresh (FG_KNOP_REFRESH_V1): once per ES sweep after the
   frequency loop (skipped on the last sweep), re-solves every
   ray's lagged field against the committed store. Reported
   KN_RF_RESID = MAXVAL over rays of the last-pass RES - i.e. the
   max normalized *iterate change*, a lagged fixed-point map
   metric, NOT a transfer-equation residual and NOT a residual of
   the coupled wavelength boundary problem. The 0.98 from Phase
   2.4 is this quantity. The independent fixed-point residual
   (KN_Q_RESID diagnostic) re-applies the map once to the
   committed field - the closest available "equation residual";
   enabled for the K runs.
3. Outer ES loop (NUM_ES_ITERATIONS): the refresh runs between
   sweeps, so NUM_ES=2 yields exactly one refresh and NUM_ES=4
   yields three. No convergence check terminates the ES loop -
   it is a fixed count. The OBS render has no outer population
   iteration beyond this. A separate full-field metric exists:
   ES_RESID_DIAG (enabled in the K runs) snapshots each sweep's
   ES_J_CONV generation and reports per-depth and global
   max/median relative change between consecutive sweeps - the
   direct measure of sweep-to-sweep field drift required by
   Task 5, independent of the refresh-internal iterate residual.

Termination in Phase 2.4: NUM_ES limit (2 sweeps, 1 refresh) - not
a convergence criterion. Blockers: cold-started I_NEXT=0; refresh
rays exiting stall/cap/plateau without commit; adaptive damping
flooring OMEGA_RAY on hot frequencies; fine OBS-grid |GAM|/dLN
modes.

## 8. KNOP maturity experiments

Run matrix (exe b409293a14eaf3e0, OBS_FG_KNOP=T, J_RESCUE_RESTORE=F,
KN_Q_RESID=T, ES_RESID_DIAG=T, identical frozen inputs):

- K1 `caseK1_knop_nores_es2`: NUM_ES=2, cold start - current-config
  control on the Phase-2.5 executable.
- K2 `caseK2_knop_nores_es4`: NUM_ES=4, cold start - jobs
  13346271/13346469, completed 16:38 (rc=0).
- K3 `caseK3_knop_nores_restart`: NUM_ES=2 warm from K2's
  FGINEXT_GS - job 13348330.
- K4 `caseK4_knop_nores_restart4`: NUM_ES=4 warm from K2's
  FGINEXT_GS - job 13348331.

Restart mechanics: the OBS restore accepts 'FGINEXT_GS' in the
working directory or '../FGINEXT_GS' (run root). The run root is
the only staging point that survives render setup - batobs.sh
runs `find . -maxdepth 1 -type f -delete` over obs/ at startup
(keeping only batobs.sh, CMF_FLUX_PARAM_INIT, cmf_OBS logs and
dated spectra), which removes any checkpoint staged inside obs/.
Both warm runs report 'I_NEXT checkpoint present: YES
../FGINEXT_GS', 'I_NEXT valid: YES' (R/V/SIGMA/COMP_DEPTH_IDX
fingerprint matched), 'RESTORED COMPLETE I_NEXT GENERATION
(tag=0)', OBS remap 74704 -> 74704.

K2 sweep history (OUT_FLUX refresh blocks + ES_RESID_DIAG):

| sweep | refresh resid (iterate) | KN_Q_RESID (indep.) | ES_RESID max rel chg | argmax depth |
|---|---|---|---|---|
| 1 | 0.9767 | 1.460 | - (first gen) | - |
| 2 | 0.9769 | 3.354 | 0.864 | 62 |
| 3 | 0.8067 | 1.043 | 0.717 | 62 |
| 4 | - (last sweep, no refresh) | - | 0.595 | 62 |

Iterate residual stays ~0.98-0.81 (tol 2e-3): not converged
after 3 refreshes. The independent fixed-point residual is
non-monotone (1.46 -> 3.35 -> 1.04) and O(1) - the committed
field is nowhere near a fixed point. The sweep-to-sweep field
change at depth 62 decays slowly (0.86 -> 0.72 -> 0.60, ~0.82
per sweep geometric factor): even if the decay held, reaching
1e-3 relative change needs roughly 20 more sweeps. K2 rescue
events: 9,340 across 4 sweeps (2,337/2,192/2,313/2,498 - not
decaying).

K3 warm continuation (sweeps 5-6 of the chain, completed rc=0):

- Refresh after chain-sweep 5: iterate residual 0.76597
  (improved from 0.98) but independent q residual 5.9497 -
  WORSE than every cold-run value. Identical numbers appear in
  K4's first refresh (same checkpoint, same sweep) - a good
  determinism cross-check.
- ES_RESID chain-sweep 6: GLOBAL_MAX_REL=0.659 at depth 14 -
  still ~66% per-sweep change; the argmax moved off the
  inversion column (K2's depth 62).
- K3 rescue events: 8,704 across 2 sweeps.

K4 warm continuation (sweeps 5-8 of the chain, completed rc=0):

- Refresh residuals: 0.7660 -> 0.7642 -> 0.7575 (tol 2e-3).
- Independent q residual: 5.950 -> 5.632 -> 0.986.
- ES_RESID GLOBAL_MAX_REL: 0.659 @ d14 (gen 2), 0.566 @ d1
  (gen 3), 0.490 @ d63 (gen 4) - slowly decaying, and the
  argmax wanders (14 -> 1 -> 63) rather than pinning on the
  inversion column.
- K4 rescue events: 17,519 across 4 sweeps
  (4,680/4,024/4,279/4,536) - rising on the last two sweeps.

K1 cold control (NUM_ES=2, completed rc=0): single refresh
residual 0.9767, independent q 1.46018 - bit-identical values
to K2's first refresh (same cold start, same generations 1-2),
which is the expected determinism check. Rescue events 4,529
(2,337/2,192).

All four runs terminate GLOBAL_CMF_GRID_UNDERCOVERAGE (rc=0)
- the expected production termination, recorded as such, not a
completed model.

Complete chain history (iterate / independent):

| chain sweep | refresh resid | KN_Q_RESID | ES_RESID max rel | run |
|---|---|---|---|---|
| 2 | 0.977 | 1.46 | 0.864 @ d62 | K2 |
| 3 | 0.977 | 3.35 | 0.717 @ d62 | K2 |
| 4 | 0.807 | 1.04 | 0.595 @ d62 | K2 |
| 5 | 0.766 | 5.95 | 0.659 @ d14 | K3/K4 |
| 6 | 0.764 | 5.63 | 0.566 @ d1 | K4 |
| 7 | 0.757 | 0.99 | 0.490 @ d63 | K4 |

The lagged map is not contracting: the iterate residual decays
by a factor ~0.96 per refresh (0.977 -> 0.757 over six
refreshes); at that rate ~150 further refresh generations would
be needed for the 2e-3 tolerance. The independent fixed-point
residual oscillates at O(1) (1.46, 3.35, 1.04, 5.95, 5.63,
0.99) with no trend - the committed field is nowhere near a
fixed point and the warm restart did not improve it (5.95 at
chain sweep 5 is the worst value in the chain). Rescue-event
counts do not decay with maturity; in the warm chain they grow
(4,024 -> 4,279 -> 4,536) and the fired argmax population
spreads inward to depths 66-68. This is the same bounded
non-contraction previously seen on the USE_LAM_ES=F moment path
- the I_NEXT lagged map on this model does not converge within
the available damping/pass budget.

Per-run summary:

| run | config | sweeps | rescue events | OBSFLUX sha256 prefix |
|---|---|---|---|---|
| K1 | es2 cold | 2 | 4,529 | 6adf9511b39b0771 |
| K2 | es4 cold | 4 | 9,340 | fbda89069b433947 |
| K3 | es2 warm (K2 ckpt) | 2 | 8,704 | 7ff4adb3502a7493 |
| K4 | es4 warm (K2 ckpt) | 4 | 17,519 | 1eb4acae1b8cf042 |

Rescue-event argmax depth histograms (all sweeps): K1
{68:1186, 70:481, 69:440, 67:65}; K2 {67:1109, 68:688,
70:514, 69:156}; K4 {67:1921, 70:1673, 68:484, 66:317,
69:97}. The KNOP fired population concentrates at depths
66-70 and shows no peak at depth 62 - see section 9a.

## 9. Remaining rescue-event analysis (Phase 2.4 data)

caseC mature sweep: 2,192 fired triggers, all trig=3 (absolute
RESC_JMAX cap). Argmax depth histogram: 68 (1186), 70 (481),
69 (440), 67 (65) - 2,172/2,192 = 99% in the inversion column
67-70; 20 events at depths 22, 55, 56, 59.

Field at argmax vs local quantities:

- J/zeta: median 2.58, p5 1.25, p95 4.88, max 1.9e5. Median
  event is a field ~2.6x the local thermal source crossing the
  absolute cap - a legitimate peak pattern, not the V12 ~300x
  runaway. Tail: 16 freqs J/zeta>435 (0.7%), 24 freqs >100 (1.1%)
  - the candidate genuine-spike class, to be re-audited at maturity.
- J/prevf: median 0.0037 - rescued J is ~270x below the
  previous-frequency field: the trigger is an absolute-cap event,
  not a chain-propagation signature.
- J/J(ND): median 1.03 - profiles are boundary-flat.
- sigma_true at argmax: median -31.3 (range -31.8..-3.4) - all in
  strongly compressing cells.

Interpretation: remaining KNOP triggers are dominated by modest
J/zeta peaks in the inversion column, not solver runaway; the
high-J/zeta tail (~1%) may be immature-map spikes and is the
population the K-runs' maturity trend will test. Against the
local scattering equilibrium J_eq = zeta/(1-omega) the fired
median is J/J_eq = 0.006 (p95 0.011, 0.7% above) - the absolute
cap is triggering on fields far inside the diffusion-limit
value, i.e. the trigger population is a ceiling artifact of the
cap's absolute scale at these zeta values, not field failure.
Whether the high-J/zeta tail shrinks with maturity is measured
in section 8.

## 9a. Index-62 investigation (verified Phase 2.3 + Phase 2.5 data)

Phase 2.3 verification established the event-depth statistics
of the V12 rescue population (audit doc
CMF_FLUX_J_RESCUE_THRESHOLD_ATTRIBUTION.md, commit 03440d56;
independently reproduced from rmax_1000/J_RESCUE_EVENTS):

- 75,738 events over 13 renders; all trigger the absolute cap
  (CAP_PREV never binds); all at sigma_true < -1 (median -31,
  worst -155.6 at depth 62).
- Per-event local equilibrium J_eq = etab/(chib - chi_scat),
  with chi_scat pure Thomson verified against the LUER dump:
  98.5% of events sit BELOW their local J_eq; 1,022 events
  exceed it; 894 of the 1,022 sit at depth 62.

Index-62 cell properties (this model, index 62 = strongest
compression cell):

- sigma_true ~= -155.6 (median -31.0 over all events);
  clamped V12 sigma = -0.999.
- J_eq event-depth distribution over the 2,040 bad=62 events:
  median J/J_eq = 0.78, p90 = 4.31, max = 11.2; 894 events
  exceed J_eq (exactly reproducing the Phase 2.3 count).

Direct V12-vs-KNOP comparison on the SAME 2,040 frequencies
(mature-sweep JFRQ rows; KNOP data from caseC, caseK2):

- KNOP argmax depth on these frequencies is NEVER 62 in caseC
  (0 of 2,040) and below the top-8 in caseK2 (<8 of 2,040).
  KNOP places its argmax at depth 70 for ~95% of them (caseC
  1948, caseK2 1684), with the remainder at 66-69 or
  uncompressing depths.
- Field magnitudes: median V12 J(62) = 7.31e-3; median KNOP
  GLOBAL maximum over all depths = 3.84e-4 on the same
  frequencies - V12's compression-cell value is ~16x KNOP's
  largest value anywhere. KNOP J(62) is therefore at least
  ~16x below V12 J(62), and below the 1.0e-3 cap.

Interpretation: the depth-62 super-equilibrium spike is a
V12-specific feature. It has no counterpart in the solver that
carries the true (unclamped) velocity gradient and the correct
GAM<0 upwind direction: in KNOP the strongest compression cell
produces no local field maximum. This supports classifying the
894-event population as a formal-solver artifact of the
sigma clamp at sigma_true << -1, not a physical radiation
field. It does NOT establish that KNOP's J(62) is correct -
KNOP remains unconverged (section 8) - only that the spike is
absent under the mathematically consistent treatment.

The converse population also differs: the 98.5% sub-equilibrium
V12 events at depths 63-70 (median J/J_eq ~ 0.04 at depth 70)
sit far below the local thermalization scale in a region with
tau_abs ~ 2e5 across the boundary cell. Whether that
sub-equilibrium field is a second V12 defect (under-populated
diffusion solution) or a legitimate non-local transfer result
is not yet decided - the KNOP field at those depths is itself
unconverged, so no two-sided conclusion is drawn.

## 10. Scattering-source analysis

Implemented decomposition (formal solve): S = ETA/CHI with
SOURCE = ZETA + THETA*RJ, THETA = CHI_SCAT/CHI (single-scattering
albedo omega), ZETA = thermal part. Depth 70 (Phase 2.4 data):
Chi(es)=0.513, Chi(tot)=0.514 -> omega ~= 0.9977.

For a source S = zeta + omega*J, the local equilibrium in a
pure-scattering optically thick medium is J_eq = zeta/(1-omega).
This is the intensity the local solution converges to, not an
upper bound: J can exceed zeta/(1-omega) where non-local
contributions (boundary radiation, advected donors) feed it, and
it is below where the medium is thin. At depth 70, 1/(1-omega) ~
435 is therefore the correct *scale* for local diffuse
equilibrium, used here as a plausibility discriminator.

J/J_eq statistics over the fired (cap-triggering) frequencies:

- A (V12+restore): median J/J_eq = 0.90, p95 2.5 - 47% of the
  rescued fields exceed the local scattering equilibrium; the
  pinned tail sits essentially AT the diffusion-limit value.
- B (V12 no restore): median 9.1, p95 112 - the unrescued field
  climbs far above local equilibrium (90% of fired freqs above).
- C (KNOP no restore): median 0.006, p95 0.011 - rescued events
  are fields orders of magnitude below local equilibrium (the
  absolute cap catches modest peaks, not runaways); 0.7% exceed.
- D (KNOP+restore): median 0.0008 - the pins sit far inside the
  equilibrium scale.

The scattering equilibrium is frequency-dependent: at depths
63-70 the continuum opacity splits between electron-scattering
windows (chi ~ 0.513, omega ~ 0.998, J_eq ~ 435*zeta) and
bound-free edges (chi ~ 40, omega ~ 0.01, J_eq ~ zeta). The
argmax events concentrate in the scattering windows. tau from
depth 68 to the boundary is ~45 - the region is optically thick,
so local equilibrium arguments are applicable.

The ratio is not a hard bound and is not used as one; the decisive
evidence remains the zeta-profile tracking (KNOP follows the
source; V12 decouples to the cap).

Reconciliation with the Phase 2.3 event-depth measurement: the
caseA figure above (median J/J_eq = 0.90 at argmax depth)
samples JFRQ rows at each frequency's argmax depth, while the
Phase 2.3 verification (98.5% below, 1.5% above) samples JEV
rows at the rescue-event depth over all sweeps. Both windows
are consistent: argmax depth 70 carries the largest J and the
largest J_eq, while the event depth selects cells where the
absolute cap actually fired. The fixed point of the
implemented source form is exact: SOURCE = ZETA +
CHI_SCAT*RJ/CHI (COMP_JCONT_V4.INC:562, :1349) gives
J -> ZETA/(1 - CHI_SCAT/CHI) = ETA/(CHI - CHI_SCAT) = J_eq
wherever the local iterate dominates non-local input.

## 11. Final converged spectral comparison

No converged KNOP spectrum exists to compare: section 8 shows
the lagged map still O(1) from a fixed point after eight
chained sweeps. What can be reported is the mutual divergence
of the K-run OBSFLUX products (sha256, all distinct):

- K1 6adf9511b39b0771, K2 fbda89069b433947,
  K3 7ff4adb3502a7493, K4 1eb4acae1b8cf042.

File-level comparison of the OBSFLUX products shows order-unity
content differences between successive generations (K2 vs K3
median relative difference ~0.29, K1 vs K4 ~1.38 on
commonly-gridded records) - consistent with the sweep-to-sweep
ES_RESID changes of ~50-66%. A V12-vs-KNOP final-spectrum
comparison is deferred until the lagged map reaches a fixed
point; drawing one now would compare a snapshot of a
non-contracting iteration, not a solution.

## 12. Phase 2.3 / Phase 2.4 count reconciliation

Reconciled exactly:

- J_RESCUE_EVENTS (JEV) records one line per rescue event per
  sweep; J_RESCUE_FREQ (JFRQ) records one summary line per
  frequency per sweep.
- caseA JEV total = 75,738 = JFRQ sweep-1 fired (44,924) +
  sweep-2 fired (30,814). Verified by per-sweep JFRQ recount.
- Phase 2.3's "75,738" was the two-sweep event total; Phase 2.4's
  "30,814" is the mature-sweep unique-frequency count. Same grid,
  same exe family, same events - different counting windows.
- Restoration shrinks the trigger set between sweeps
  (44,924 -> 30,814) because capped donors feed fewer downstream
  rescues. Without restore (caseB) the same 46,159 freqs fire on
  both sweeps (92,318 total) - the unrescued field re-triggers
  identically.
- caseB sweep-1 (46,159) > caseA sweep-1 (44,924): even within
  sweep 1 the capped-donor chain suppresses some triggers.

No definition or sampling inconsistency remains; counts may be
compared across phases only after fixing the window (mature sweep
vs all-sweep total).

## 13. Recommendation for further development

Findings:

1. Equation audit: the implemented equation is the nu-form
   drift equation dI/ds = eta - chi_hat*I + Q*I_upwind with
   chi_hat = chi + |GAM|/DLN. The explicit a*I-proportional
   piece of Knop's lambda-form 4*a*I is omitted by the same
   convention as stock V12/V13; its conservative same-bin part
   is inside the |GAM|/DLN upwind diagonal (verified against
   the in-source audit note, fg_j_cmf_nonmono_v1.f:89-101).
   Residual size ~3*dLOG_NU ~ 0.025% of the velocity-coupling
   magnitude - cannot produce the observed 40-100x V12/KNOP
   field differences. Not the defect to chase first.
2. The depth-62 super-equilibrium spike is V12-specific
   (section 9a): absent from KNOP on the same frequencies,
   consistent with sigma-clamp breakdown at sigma_true = -156.
   The V12 artifact class is now verified, not conjectured.
3. KNOP is not converged (section 8): iterate residual ~0.96
   per refresh, independent q residual oscillating O(1) after
   eight chained sweeps; rescue counts grow in the warm chain.
   KNOP J smaller than V12 J is therefore a signature of a
   different (possibly correct) sign/upwind treatment, not
   proof of correctness.
4. The remaining KNOP fired population (sections 8-9) is
   dominated by modest sub-equilibrium peaks at depths 66-70 -
   plausibly a ceiling artifact of the absolute cap's scale,
   but that verdict requires a converged field.

Recommended order of work:

- Solve the I_NEXT lagged-map non-contraction first (contract
  factor, damping schedule, or a Jacobian-coupled update);
  every downstream question - the fired population, the
  spectrum, any per-cell rescue bound - is ill-posed while the
  map does not converge.
- Retain RESC_JMAX = 1.0e-3 (user decision, Phase 2.3). Keep
  the J_eq measurements as the baseline scale for judging
  KNOP's mature field, not as a threshold.
  SUPERSEDED 22-Nov-2026: the m15 interacting-model A/B showed
  the 1e-3 rail clipping verified-legitimate fields (~1.2e-1)
  and collapsing the production spectrum ~100x; the default is
  now 1e2 (commits 310a9364/f2471abf/cc5eafc2). The Phase 2.3
  decision predates the old-exe truth reference; explicit
  RESC_JMAX=1e-3 in CMF_FLUX_PARAM still restores it per run.
  The J_eq-as-baseline-scale recommendation stands unchanged.
- Do not implement a per-cell kappa*J_eq bound now: it would
  hide the same class of solver defect the absolute cap hides.
- Once KNOP converges, redo this audit's comparisons at the
  fixed point: J(62) vs J_eq, the fired population, and the
  OBSFLUX product against the V12 baseline.
- The three-seed probe experiment
  (KNOP_THREE_SEED_EXPERIMENT.md, live) addresses
  initialisation sensitivity of the probe-ray solve, a
  different question; its results complement but do not
  replace the global maturity criterion used here.
