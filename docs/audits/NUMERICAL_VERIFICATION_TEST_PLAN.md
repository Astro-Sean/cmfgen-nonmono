# NUMERICAL_VERIFICATION_TEST_PLAN.md

Numerical verification and stress-test programme for the CMFGEN
non-monotonic velocity implementation (`cmfgen_jfnk`).  The goal is a
rigorous, evidence-based demonstration that the non-monotonic solvers
are numerically stable, physically consistent, and robust across the
configurations needed for interacting supernovae (Type IIn, shocked
ejecta, CSM interaction, strong compression, reverse shocks, multiple
velocity reversals).

Ground rules (owner directive):

- A successful build does not prove numerical correctness.
- A converged model does not prove physical correctness.
- A unit test that only checks for finite output is insufficient.
- Tests compare numerical results against known analytical solutions,
  independently assembled reference solutions, conservation
  identities, or independently evaluated residuals.
- The equations and the complete numerical pathways are tested, not
  the presence of source-code statements.
- One test at a time.  The programme is executed in the TASK order in
  section 10; do not implement the whole plan in one session.

## Status labels

- `[ ]` TODO
- `[~]` IN PROGRESS (partial evidence exists; test not complete)
- `[x]` PASS
- `[!]` FAIL
- `[-]` NOT APPLICABLE

## Evidence ladder (per-pathway maturity)

Every solver pathway is classified by the highest stage it has
demonstrated.  A pathway must not be called robust below stage 4.

| Stage | Meaning |
|---|---|
| 1 | Source inspected |
| 2 | Compiles |
| 3 | Unit tested |
| 4 | Manufactured-system tested (analytic/reference solutions) |
| 5 | Integration tested (full solver pathways) |
| 6 | Production-model tested (complete CMFGEN run) |
| 7 | Physically validated (independent physics checks on a model) |

## General test requirements

Every numerical test must include, where applicable: a positive case,
a deliberately failing (negative) case, a limiting case, a degenerate
case, an independent reference, a quantified tolerance with stated
justification, and a clear failure message.  Floating-point checks
use justified absolute/relative tolerances; randomized tests use
deterministic seeds; skips are reported, never silent.

---

## 1. Manufactured velocity-profile library - V001..V012 `[x]`

Common library: `unit_tests/velocity_profiles.py`.  Every
non-monotonic solver test draws its input velocity field from this
registry so the same profiles exercise every pathway.

Conventions: R(1) outer > R(ND) inner (CMFGEN outside-in indexing);
V in km/s; R normalized; sigma = d ln V / d ln R - 1 (CMFGEN
SIGMA_TRUE).  All profiles are smooth analytic functions of R, so R,
V, dV/dR and sigma are known in closed form.  The library computes,
independently of the production code: R, V, analytic dV/dR, sigma,
every gradient sign change (interpolated location), sigma min/max,
and every compression region (sigma < -1, i.e. dV/dR < 0) with its
width in R and ln R.

Dip family: `V = a R prod_k (1 - A_k exp(-0.5 (ln(R/R_k)/w_k)^2))`.
Jump family: `V = a R (1 - d_j/2 (1 + tanh(ln(R/R_s)/w)))`.

| ID | Definition (default parameters) | nsign | ncomp | sigma range |
|----|------------------------------|-------|-------|-------------|
| V001 | V = 8e3 R (homologous) | 0 | 0 | 0 |
| V002 | V = 1.2e4 (1 - 0.02/R)^1.5 (beta law, r_in=0.1) | 0 | 0 | [-0.97,-0.62] |
| V003 | V = 4e3 R^-0.15 (power-law deceleration) | 0 | 1 | -1.15 |
| V004 | V = 8e3 R (1 - 0.55 G(0.35,0.10)) single reversal | 2 | 1 | [-5.23,5.08] |
| V005 | V = 8e3 R (1 - 0.90 G(0.35,0.06)) strong compression | 2 | 1 | [-28.8,21.8] |
| V006 | V = 8e3 R prod(1-0.5 G_k), R_k=0.60/0.35/0.15, w=0.08 | 6 | 3 | [-5.3,5.7] |
| V007 | V = 8e3 R (1 - 0.8 S(0.35,0.03)) narrow shock | 2 | 1 | [-25.4,0] |
| V008 | V = 8e3 R (1 - 0.8 G(0.35,0.40)) broad shell | 2 | 1 | [-3.0,3.0] |
| V009 | V = 5e3 (1 + 3e-3 tanh(ln(R/0.35)/0.1)) near-flat | 0 | 0 | [-1.00,-0.97] |
| V010 | V = 8e3 R (1 - 0.985 G(0.35,0.30)) stagnant region (V~44) | 2 | 1 | [-18.3,18.7] |
| V011 | V = 8e3 R prod(1-0.6 G_k), R_k=0.42/0.30, w=0.08 | 4 | 2 | [-7.3,7.6] |
| V012 | V = 8e3 R (1 - 0.9 S(0.35,0.06)) strong smooth jump | 2 | 1 | [-15.3,0] |

Parameter sweeps supported by the generator: ND (radial resolution),
grid kind (log/linear), r_shock (shock position), width (shock
width), amp/dj/eps (compression strength), seeded smooth
perturbation, reversed array, exactly-zero-V cell, collapsed
spacing, negative-V (infall) amplitude.

### Validation results (2026-11-04, `test_velocity_profiles.py`, 103 checks)

- Library analytic dV/dR vs complex-step derivative: max rel err
  ~1e-14 on all 12 profiles (independent code path).
- Production `SIGMA_TRUE_CORE` (drv_sigma_true `core` mode) vs
  python mirror `sigma_true_core_mirror`: bitwise identical on all
  profiles and on the linear grid.
- Production vs analytic sigma on well cells (|sigma|<1E3), ND=80
  log grid, and refinement ratio (ND 160->320):

| ID | err@80 | rate | | ID | err@80 | rate |
|----|--------|------|-|----|--------|------|
| V001 | 5.7e-14 | exact | | V007 | 1.0e+01 | 3.5 |
| V002 | 1.9e-04 | 4.0 | | V008 | 6.1e-02 | 4.0 |
| V003 | 9.1e-14 | exact | | V009 | 2.2e-03 | 3.9 |
| V004 | 7.6e-01 | 3.9 | | V010 | 5.0e+00 | 3.4 |
| V005 | 1.7e+01 | 3.5 | | V011 | 1.8e+00 | 3.9 |
| V006 | 1.2e+00 | 3.9 | | V012 | 1.7e+00 | 3.8 |

  Errors are stencil truncation (quadratic in ln R): h^2 rate
  confirmed; steep-profile constants reflect sigma range, not a
  defect.  Sub-resolved shock widths are mirror-exact with the
  expected smoothing (documented in test D2).
- Census checks: sign-change count/locations correct on all
  profiles and on every sweep (position, width, strength, close
  reversals, near-merged reversals); N_UNDEF=0 for valid inputs;
  N_UNDEF=ND for reversed grid; N_UNDEF=1 for a V=0 cell; N_UNDEF=2
  for a collapsed spacing; aborts fail-closed under a non-mono
  flag; infall pocket (V<0) finite, |sigma|<=1E6, warns.
- Deterministic seeded perturbation confirmed bit-identical.

Failure modes covered: wrong stencil weights (mirror diff), wrong
sign at boundaries (analytic diff), silent monotonization (census
mismatch), uninitialized values (N_UNDEF sentinel + finiteness),
grid degeneracies (reversed/collapsed), sign mishandling at V<=0.

Test: `unit_tests/python/test_velocity_profiles.py`
(`VELOCITY-PROFILE LIBRARY + SIGMA_TRUE_CORE`, 103/103 PASS).

---

## 2. KNOP formal-transfer tests - KN001..KN014

Primary source: `plane/fg_j_cmf_nonmono_v1.f` (`FG_J_CMF_NONMONO_V1`,
`KNOP_RAY_SOLVE_V1`, `FG_KNOP_REFRESH_V1`, store/refresh modules).
Harness: `unit_tests/f90/drv_knop_ray.f90` (manufactured core ray,
NI=8, NCF=6 - extend ND/NCF and opacity laws as needed).
Required diagnostics on every run: true transfer residual, radiation
flux, frequency-integrated flux, I min/max, number of frequency
upwind reversals, number of frequency-matching events, number of
failed/capped sweeps, lagged-field residual, conservation error.
Velocity inputs come from the V-profile library (V001-V012); the
shared task for this series is passing each V-profile through the
actual solver, not only the gradient helper.

**KN001** homologous limit `[~]`
- property: non-monotonic formulation reduces to the monotonic result
- routines: FG_J_CMF_NONMONO_V1 vs FG_J_CMF_V13 (stock)
- input: V001/V002 on manufactured ray, identical CHI/ETA
- expected: I/J/H identical to stock formal solve within roundoff
- reference: stock CMFGEN formal solution on the same ray (V13 path)
- tolerance: |J_new-J_stock| <= 1e-10 rel; identical discrete
  equations in the monotonic limit by construction
- detects: sign-flipped GAM handling firing on expanding flow,
  extra upwind branches active at GAM>=0, boundary differences
- type: unit (manufactured driver, paired new/ref build)
- status note: `test_knop_ray.py` proves convergence on a
  compressing band; the V13 A/B comparison is not built yet.

**KN002** constant-source-function transfer `[ ]`
- property: formal solution reproduces I(z) for prescribed S, I_bnd
- input: V001 + constant S(z), known I(0)/I(NI); slab geometry
- expected: I(z) = I_bnd e^-tau + S(1-e^-tau) along each ray
- reference: analytic exponential attenuation evaluated in python
- tolerance: |I_num-I_ana| <= ~1e-6 rel where tau resolved; the
  discrete march is exact for constant S up to the tau quadrature
- detects: source-term indexing errors, wrong tau accumulation,
  boundary intensity misapplication
- type: unit

**KN003** pure absorption `[ ]`
- property: exponential attenuation with no emissivity
- input: V001, prescribed CHI(z) profile, ETA=0
- expected: I = I_bnd exp(-tau(z)) with the production tau convention
- reference: independent quadrature of tau from the same CHI grid
- tolerance: rel err <= 1e-6; only the opacity integral differs
  between reference and production march
- detects: opacity interpolation errors, wrong march direction,
  dropped segments
- type: unit

**KN004** conservative scattering `[ ]`
- property: flux/energy preservation under pure scattering
- input: ETA -> epsilon*CHI with epsilon -> 0, V001, closed outer BC
- expected: frequency-integrated flux conserved to the discrete
  identity; J->equilibrium value
- reference: conservation identity on the discrete system
- tolerance: conservation residual <= 1e-8 rel
- detects: unphysical gain/loss terms, wrong scattering split
- type: unit

**KN005** single compression region `[~]`
- property: sign of GAM, frequency upwinding, intensity positivity
  under sigma<-1
- input: V005 mapped onto the manufactured ray (compressing band)
- expected: compressing rows take the signed-|Q| downwind store;
  I stays positive and finite; upwind-direction census = 1 band
- reference: independent recomputation of GAM sign per depth from
  sigma_true; census count
- tolerance: census equality; I>0; residual per KN diagnostics
- detects: upwind direction wrong inside compression, sign of
  chi_hat flipped, negative intensities
- type: unit
- status note: `test_knop_ray.py` runs a compressing band to
  NXIT=1; explicit upwind-census and positivity assertions open.

**KN006** multiple velocity reversals `[~]`
- property: correct upwind direction at every reversal
- input: V006 (3 bands, 6 gradient sign changes) on the ray
- expected: per-depth upwind census alternates correctly; each
  characteristic's source upwinds from the physical side
- reference: independent sign map from analytic dV/dR
- tolerance: census equality per depth
- detects: off-by-one at reversal nodes, upwind persistence across
  a sign change
- type: unit
- status note: `rev` case (two compressing bands) passes in
  `test_knop_ray.py` T9; three-band census open.

**KN007** near-zero GAM `[ ]`
- property: behaviour at/beside the GAM=0 surface
- input: V009 (sigma ~ -1, GAM ~ 0) plus variants bracketing zero
- expected: smooth coefficient/solution variation through zero; the
  exact-zero case hits the documented decoupled row
- reference: coefficient-level comparison vs GAM->0 limit
- tolerance: continuity within discretization rounding; exact-zero
  row equals the decoupled form
- detects: pole at GAM=0, sign noise, discontinuity
- type: unit

**KN008** frequency-grid convergence `[ ]`
- property: solution converges under frequency refinement
- input: V005, NCF = 6/12/24/48 on the same atmosphere
- expected: monotone approach to a grid-independent J/H/flux
- reference: Richardson estimate of the continuum limit
- tolerance: successive-difference reduction consistent with the
  discretization order
- detects: aliasing of the frequency coupling, hidden dependence on
  grid spacing sign conventions
- type: unit

**KN009** frequency doublets `[ ]`
- property: small-dLOG_NU limit at nearly coincident frequencies
- input: manufactured grid with a 1E-11-spaced doublet (mirrors the
  bound-free edge doublets on the production continuum grid)
- expected: the doublet row degenerates to a null-transport copy
  (J_k ~ J_upwind); no coefficient blow-up
- reference: limiting equation derived for M series (P3-014) applied
  to the formal march
- tolerance: |J_k - J_upwind| <= 1e-3 of a real-bin difference
- detects: division by near-zero spacing, wrong floor placement
- type: unit

**KN010** boundary-condition consistency `[ ]`
- property: inner/outer boundaries under expanding, compressing and
  reversing V
- input: V001, V005, V006 with I_M(1)=0 inner, I_P(NI)=I_CORE outer
- expected: boundary rows apply the documented formulas at all GAM
  signs; no lost flux at a resonance landing on the boundary
- reference: boundary equations from KNOP_FORMAL_SOLVER_AUDIT
- tolerance: boundary residual <= sweep tolerance
- detects: boundary sign errors, lost inflow at compressing edges
- type: unit

**KN011** sweep-order invariance `[~]`
- property: alternating and ordered sweeps reach the same fixed point
- input: V005 compressing band on the manufactured ray
- expected: committed stores agree; both NXIT classes converge
- reference: cross-comparison of the two converged stores
- tolerance: |store_alt - store_ord| <= ~3e-8 (observed level)
- detects: order-dependent operator asymmetry, store aliasing
- type: unit
- status note: covered by `test_knop_ray.py` (both sweep modes
  converge, stores agree to ~3e-8); regression pin complete.

**KN012** lagged-field fixed-point residual `[ ]`
- property: retained J_NEXT/H_NEXT satisfies the production operator
- input: converged refresh state from KN005/KN006 runs
- expected: independently evaluated residual r = Op(J*) - J* is
  below the acceptance bound, not inferred from iteration deltas
- reference: python re-assembly of the documented discrete operator
  applied to the retained stores
- tolerance: |r| <= tol_op with tol_op stated from coefficient norms
- detects: convergence declared on a drift-only metric, stale-store
  acceptance
- type: unit
- status note: depends on dumping the assembled operator row; plan
  a driver extension (M011-style residual dump for the FG path).

**KN013** restart consistency `[ ]`
- property: checkpoint/restart reproduces an uninterrupted run
- input: run to intermediate sweep, DUMP_*_GS stores, scramble,
  RESTORE_*_GS, continue
- expected: restarted and uninterrupted committed fields identical
  (bitwise or within store precision)
- reference: uninterrupted calculation
- tolerance: bitwise after restore; drift <= store write precision
  after continuation
- detects: missing state in checkpoints, fingerprint/tag bypass,
  foreign-seed acceptance
- type: unit
- status note: store-level dump/restore + fingerprint/tag rejection
  exist for the MOM path (`ckpt` mode); the FG ray-store leg is
  open.

**KN014** grid refinement `[ ]`
- property: I, flux and J converge under radial refinement
- input: V005/V007 on NI = 8/16/32/64 manufactured rays
- expected: refinement ratios consistent with operator order
- reference: analytic or fine-grid reference on the same problem
- tolerance: monotone error decrease; rate >= operator order
- detects: non-convergent discretization, grid-scale instability
- type: unit

---

## 3. Moment-solver tests - M001..M017

Primary source: `plane/mom_j_cmf_nonmono_v1.f` via
`unit_tests/f90/drv_mom_solve.f` and
`unit_tests/python/test_mom_solve.py` (extend, do not fork).

**M001** homologous limit `[~]`
- property: non-monotonic moment solve == monotonic solution
- input: `homo` mode (sigma=0); also V001 via the library grid
- expected: J/H match the monotonic moment solution
- reference: production run on `homo` + the analytic coefficient
  re-computation already in test_mom_solve
- tolerance: coefficient-level rtol <= 1e-9 (existing checks);
  solution-level vs a monotonic reference open
- detects: compression branches firing at sigma=0, sign flips
- type: unit
- status note: coefficient verification exists; explicit monotonic
  reference comparison open.

**M002** exact zero-GAM resonance `[x]`
- property: 1+sigma*K_ON_J = 0 decouples correctly
- input: `zero` mode (bitwise-resonant shelf)
- expected: resonant cells take the decoupled row (Q=PSI=0)
- reference: measured t_zero case; BUG-004 analysis
- tolerance: exact (Q=PSI=0 bitwise on resonant cells)
- detects: pole at the resonance, regularization leakage
- type: unit
- evidence: test_mom_solve t_zero (36 resonant cells decouple;
  shelf-edge floor at +0.01); P3-016/BUG-004 completed 2026-11-03.

**M003** GAM crossing `[ ]`
- property: continuous sweep through both signs of GAM
- input: parameterized sigma_true ramp crossing zero (V003 scaled)
- expected: all coefficients continuous through the crossing; final
  solution smooth
- reference: recorded coefficient series vs sign(GAM)
- tolerance: coefficient jumps <= discretization rounding
- detects: hysteresis, sign-dependent branch errors
- type: unit

**M004** near-resonant coefficients `[~]`
- property: values immediately above/below the zero-GAM resonance
- input: t_zero shelf-edge cells (T1 ~ 1e-16 off zero)
- expected: floored factor bounds |fac|>=MOM_GAM_FAC_MIN; Q bounded
- reference: gam_fac mirror (already in test_mom_solve)
- tolerance: mirror-exact on coefficient formulas
- detects: floor applied on wrong side of resonance, sign loss
- type: unit
- status note: shelf-edge cell measured (Q=1.6e-4 at T1=3.3e-16);
  dedicated both-sides sweep open (M003).

**M005** multiple compression zones `[~]`
- property: several compression regions separated by expansion
- input: `comp2` mode + V006 through the library
- expected: correct upwind/coefficient signs per zone; converged
  J_NEXT fixed point
- reference: per-zone coefficient recomputation
- tolerance: coefficient mirror-exact; residual per M011
- detects: sign bleed between adjacent zones, zone-index offsets
- type: unit
- status note: comp2 (two bands at -4/-6 sigma) verified at
  coefficient level; three-zone V006 run open.

**M006** strong compression `[~]`
- property: solution under progressively stronger compression
- input: `strong` mode (sigma=-12) + V005 amplitude sweep
- expected: converged finite J at every strength; no runaway
- reference: run record vs strength; M011 residual at each level
- tolerance: M011 residual bound; monotonic approach to a limit
- detects: divergence, quarantine runaway, coefficient overflow
- type: unit
- status note: `strong` converges in test_mom_solve; the
  progressive sweep with per-level residuals is open.

**M007** frequency doublets `[x]`
- property: dLOG_NU spanning the smallest production spacings
- input: `dbl` mode, FRQ(4)=FRQ(3)exp(-1e-11)
- expected: compressing rows at the doublet copy J_upwind
  (null-transport bin); |J3-J4| below a real-bin difference
- reference: P3-014 limiting equation
- tolerance: J copy within 0.15% (measured); documented bound
- detects: coefficient blow-up on doublets, wrong floor arm
- type: unit
- evidence: test_mom_solve `dbl` case; P3-014 closed 2026-11-03.

**M008** radial refinement `[ ]`
- property: J/H convergence with increasing ND
- input: same physical atmosphere at ND = 24/48/96
- expected: refinement ratios consistent with the discrete order
- reference: finest grid as reference; Richardson estimate
- tolerance: monotone decrease, rate documented
- detects: grid-scale instability, inconsistent boundary scaling
- type: unit

**M009** frequency refinement `[ ]`
- property: J/H convergence with increasing NCF
- input: same atmosphere at NCF = 8/16/32
- expected: approach to a grid-independent solution
- reference: finest NCF; telescoping identity (M011/P3-015)
- tolerance: documented rate
- detects: frequency-coupling aliasing
- type: unit

**M010** boundary conditions `[~]`
- property: every supported inner/outer boundary configuration
- input: driver outbc 1/2 x inbm diff/zfx/pnt(/schu) on V001+V005
- expected: each supported combination solves with the correct
  boundary rows; unsupported combos reject loudly
- reference: boundary-row audit vs code branches
- tolerance: row formulas mirror-exact
- detects: wrong boundary assembly under compression, silent BC
  substitution
- type: unit
- status note: driver exposes the options; a systematic
  boundary x boundary x profile sweep is open.

**M011** matrix residual `[ ]`
- property: independent r = A x - b after every solve
- input: any manufactured case; dump assembled TA/TB/TC rows
- expected: independently evaluated |r| below the solver's own
  convergence estimate - never the solver's estimate itself
- reference: python re-assembly of the tridiagonal system
- tolerance: |r|_rel <= 1e-8 for converged solves
- detects: false convergence, row-assembly errors the iterate
  difference hides
- type: unit

**M012** M-matrix and positivity `[ ]`
- property: assembled-matrix sign pattern where required;
  no artificial positivity repair on permitted negative quantities
- input: V005/V007 rows
- expected: off-diagonal signs match the documented stencil;
  physical negatives (lagged H couplings) preserved
- reference: sign-pattern audit per row type
- tolerance: sign census equality
- detects: sign flips in assembly, hidden clipping
- type: unit

**M013** lagged-field fixed-point consistency `[ ]`
- property: retained J_NEXT satisfies the map F(J*) = J*
- input: converged J_NEXT stores
- expected: independently evaluated map residual below acceptance
- reference: python evaluation of the lagged map on the store dump
- tolerance: <= refresh tolerance; stated bound
- detects: maturity declared on iteration count alone
- type: unit

**M014** Anderson acceleration `[ ]`
- property: inner/outer AA actually applies, helps, and is safe
- input: driver seeded with a valid J_NEXT_GS checkpoint (current
  cadence never reaches the AA gate - see P3-022)
- expected: AA counters > 0 on the entered path; accelerated
  residual decreases; rejected acceleration falls back; history
  resets after an invalidating event; converged fixed point
  unchanged vs unaccelerated reference
- reference: unaccelerated run on the same problem
- tolerance: accelerated fixed point == unaccelerated within the
  refresh tolerance; residual reduction measured
- detects: dead AA code, wrong fallback, history corruption,
  fixed-point drift
- type: unit
- status note: P3-022 gate map documented; a valid seeded
  checkpoint or a LAM_ES=F run is needed to enter the AA paths.

**M015** quarantine behaviour `[ ]`
- property: negative-J and quarantine paths are visible, not silent
- input: cells driven to negative-J / stale-field regimes
- expected: quarantined values appear in diagnostics; global field
  cannot appear mature while quarantined cells persist
- reference: census counters (MOM_NEGJ_LRG, MOM_N_NEGJ_T) vs an
  independent bad-cell scan
- tolerance: census equality
- detects: silent repair, maturity on top of quarantine
- type: unit

**M016** maturity-gate correctness `[ ]`
- property: GS maturity gates track true field consistency
- input: false-maturity and delayed-maturity scenarios
- expected: gate decision == independently computed field residual
  threshold decision
- reference: independent field residual (M013 evaluation)
- tolerance: zero false-maturity events on the test battery
- detects: premature release of the freeze, iteration-count gating
- type: unit

**M017** restart consistency `[~]`
- property: uninterrupted vs restarted moment calculations agree
- input: `ckpt` mode (2-sweep populate, dump, scramble, restore)
- expected: restored stores bitwise; continuation identical;
  fingerprint/tag mismatches reject loudly and zero the store
- reference: uninterrupted run
- tolerance: bitwise on restore (measured CKPT_MAXDIFF=0)
- detects: incomplete checkpoint state, fingerprint bypass
- type: unit
- status note: store-level round trip + fingerprint/tag rejection
  verified in `ckpt`; full-calculation restart comparison open.

---

## 4. Sobolev tests - S001..S020

Sources: `subs/sobjbar_sim.f`, `subs/sobew.f`,
`subs/sobew_grad_v2.f`, `subs/sobew_grad_v3.f`,
`subs/sob_rh_couple.f`, `new_main/mod_subs/sub_sob_line_v3.f`.
Harness: `unit_tests/f90/drv_sob_floor.f90`, `drv_sob_rh.f90`
(existing); extend for multi-resonance ray cases.

**S001** homologous expansion `[~]`
- property: sigma=0 path reduces to the upstream result
- input: V001
- expected: bit-identical tau/BETA to reference for sigma>=0
- reference: pristine reference routine + analytic BETA
- tolerance: bitwise (documented upstream behaviour preserved)
- detects: tau-floor branch leaking into monotonic flow
- type: unit
- status note: bit-identity for sigma>=0/-0.5 is in
  test_sob_tau_floor; a V001-profile-level check is implicit.

**S002** monotonic beta-law `[ ]`
- property: beta-law sigma profile handled correctly
- input: V002
- expected: tau_sob matches the analytic expression pointwise
- reference: analytic tau + independent quadrature
- tolerance: <= 1e-6 rel where resolved
- detects: sigma mis-indexing on a smoothly varying gradient
- type: unit

**S003** sigma approaching -1 `[x]`
- property: tau near the pole entry
- input: sigma -> -1 from above and below at several mu
- expected: floored response bounded, continuous
- reference: analytic floored formula (Dessart+2015)
- tolerance: BETA vs formula <= 6e-12 (measured)
- detects: discontinuity entering the pole, wrong floor arm
- type: unit
- evidence: test_sob_tau_floor T8 straddle (2026-11-03).

**S004** sigma below -1 `[x]`
- property: tau in the compression regime
- input: sigma < -1 cells
- expected: floor MAX(|1+sigma*mu^2|, -0.2*(1+sigma), 1e-4) applied
- reference: analytic floored formula
- tolerance: as S003
- detects: floor arm selection, missing absolute value
- type: unit
- evidence: test_sob_tau_floor T3/T7.

**S005** exact optical-depth floor activation `[x]`
- property: pure-1e-4 arm at the sigma=-1 pole
- input: sigma = -1 exactly, near-radial rays
- expected: tau capped at GAMH*1e4; response bounded
- reference: analytic floor response
- tolerance: <= 2-4e-3 absolute on BETA (measured)
- detects: unbounded tau at the pole, NaN at 1+sigma*mu^2=0
- type: unit
- evidence: test_sob_tau_floor T7.

**S006** continuity of the implemented prescription `[x]`
- property: no numerical jump crossing the floor boundary
- input: sigma straddling -1, full mu grid
- expected: bounded continuous BETA response
- reference: unclamped reference off the pole + floored on it
- tolerance: tracks unbounded reference to 1e-5 off the sliver
- detects: jumps at the floor boundary, sign asymmetry
- type: unit
- evidence: test_sob_tau_floor T8.

**S007** angular integration accuracy `[ ]`
- property: the mu quadrature sums BETA/weights correctly
- input: sigma<-1 profile where tau(mu) is strongly anisotropic
- expected: angle-integrated BETA matches high-order quadrature
- reference: independent fine-mu quadrature in python
- tolerance: <= quadrature order error
- detects: wrong angular weights, index errors in the mu loop
- type: unit

**S008** escape-probability reference `[~]`
- property: BETA vs the EXF escape factor
- input: homologous + compressing cells
- expected: BETA = sum_L w_L EXF(GAM_L) with floored tau
- reference: analytic expression
- tolerance: 6e-12 (measured T3/T7)
- detects: wrong escape-probability functional form
- type: unit
- status note: verified in test_sob_tau_floor; the V-profile leg
  (sigma variation along a ray) is open.

**S009** source-function comparison `[ ]`
- property: S_line = ETAL_MAT/CHIL_MAT used correctly
- input: manufactured levels with known source function
- expected: rates/EW consistent with prescribed S
- reference: analytic S-driven solution
- tolerance: <= discretization error
- detects: swapped opacity/emissivity, wrong S index
- type: unit

**S010** multiple resonance points along a ray `[~]`
- property: one line resonating at several depths deposits correctly
- input: V004/V006 (re-entrant velocity ranges)
- expected: per-resonance tau applied at each location; energy
  balance across resonances
- reference: independent multi-resonance quadrature
- tolerance: <= 1e-6 rel per resonance
- detects: missed/double-counted resonances, ordering errors
- type: unit
- status note: multi-resonance pair geometry covered via
  test_sob_rh_couple W-profile; per-resonance tau deposition on a
  ray is open (P2-006 records the audit + the architectural
  two-point limit).

**S011** two resonances with overlapping optical-depth regions `[ ]`
- property: interacting resonance taus on one ray
- input: two resonances whose tau wings overlap
- expected: combined attenuation handled per the documented kernel
- reference: independent coupled quadrature
- tolerance: documented kernel accuracy; deviations diagnosed, not
  silent
- detects: unphysical independent-resonance assumption
- type: unit

**S012** three or more resonance points `[ ]`
- property: 3+ coincident resonances on a ray
- input: tuned triple coincidence (V006 variant)
- expected: pairwise kernel applies; the known third-extinction
  omission is diagnosed or bounded (Feldmeier & Nikutta regime)
- reference: three-resonance independent quadrature
- tolerance: documented approximation bound, explicitly reported
- detects: silent omission of higher-order extinction
- type: unit

**S013** closely spaced resonance points `[ ]`
- property: resonance spacing near the discretization limit
- input: V011 mapped to a line pair
- expected: both resonances counted or a loud resolution diagnostic
- reference: independent quadrature
- tolerance: census equality where resolved
- detects: merged-resonance miscounting
- type: unit

**S014** negative opacity `[~]`
- property: NEG_OPACITY handling
- input: rows with CHIL<0 (maser regime proxy)
- expected: the documented flatten/suppression path; no NaN rates
- reference: documented policy per code path
- tolerance: policy-exact
- detects: silent negative-tau, NaN propagation
- type: unit
- status note: suppression covered in test_sob_rh_couple; the
  line-level NEG_OPACITY path is open.

**S015** zero opacity `[~]`
- property: zero-opacity rows
- input: CHIL=0 cells
- expected: zero rate contribution; no division by CHIL
- reference: documented path
- tolerance: exact zero contribution
- detects: 0/0, spurious coupling
- type: unit
- status note: BB_COR=0/ETAL=0 safety covered in test_sob_rh_couple.

**S016** near-zero velocity gradient `[ ]`
- property: sigma ~ -1 cells inside the Sobolev evaluation
- input: V009
- expected: tau uses the documented floor; no amplification of a
  nearly zero gradient into noise
- reference: analytic floored tau
- tolerance: formula-exact
- detects: sign noise near the pole, tau jitter
- type: unit

**S017** non-local RH coupling geometry `[~]`
- property: pair selection, mu^2 identity, dropped-pair policy
- input: W-shaped V (two reversals); stagnation-cell variants
- expected: pair count/symmetry vs O(ND^2) scan; V<=0 pairs dropped
- reference: independent geometry enumeration
- tolerance: DRC vs recomputed <= 4e-12 (measured)
- detects: geometry fingerprint errors, wrong endpoint test
- type: unit
- status note: covered by test_sob_rh_couple (108 pairs, 14 max
  partners, 31% stagnation drops measured); the mathematical domain
  of validity of RH78 under arbitrary non-monotonic V is a standing
  audit item - unsupported geometries must reject or diagnose.

**S018** RH coupling cache invalidation `[x]`
- property: geometry cache rebuilds when R/V/SIGMA change
- input: same-R,V call with new sigma; phase driver
- expected: fingerprint detects the change; stale tau not reused
- reference: recomputed kernel
- tolerance: exact rebuild on any input change
- detects: stale cache, wrong DRC after grid revision
- type: unit
- evidence: BUG-002 fix + regression (2/3 checks fail on pre-fix,
  12/12 pass on fixed), 2026-11-03.

**S019** non-local rate contribution to SE residual `[~]`
- property: DRC term in the statistical-equilibrium residual
- input: manufactured coupled level pair
- expected: DRC enters STEQ residual only (not ZNET_SIM), atoms
  conserved across the pair
- reference: recomputed DRC + antisymmetry identity
- tolerance: antisymmetric assembly exact
- detects: sign error in net-rate assembly, Jacobian contamination
- type: unit
- status note: audited + pair-level tested; a direct
  residual-evaluation test is open.

**S020** restart and repeated-call consistency `[ ]`
- property: repeated calls and restarts reproduce identical rates
- input: identical inputs called twice; checkpoint round trip
- expected: bitwise-identical results; caches keyed correctly
- reference: repeated evaluation
- tolerance: bitwise
- detects: cache/state leaks between calls
- type: unit

---

## 5. JFNK tests - J001..J028

Sources: `new_main/mod_subs/jfnk_solver_v1.f90`,
`jfnk_cmfgen_interface_v1.f`, `jfnk_solve_corr_v1.f`,
`new_main/cmfgen_sub.f`, `jfnk_op_probe_v1.f90`.
Harnesses: `tests/jfnk_core/`, `unit_tests/se_driver`,
`test_jfnk_solver.py`, `test_jfnk_opint.py`,
`test_jfnk_regressions.py`.  Existing coverage is extensive; the
production operator needs live-state validation (J011/J012/J028).

**J001** exact linear residual `[x]`
- property: F(x)=Ax-b evaluated exactly
- reference: direct multiplication on manufactured systems
- tolerance: machine precision
- detects: residual indexing/sign errors
- type: unit - covered by se_driver surrogate systems.

**J002** nonlinear residual, known solution `[x]`
- property: F(x*) = 0 on a manufactured nonlinear system
- reference: analytic solution of the test problem
- tolerance: |F(x*)| <= 1e-12
- detects: residual assembly errors
- type: unit - analytic stiff problem in tests/jfnk_core.

**J003** FD JVP vs exact Jacobian `[x]`
- property: finite-difference JVP approximates J v
- reference: exact Jacobian of the surrogate problem (DIR battery)
- tolerance: error curve vs h with documented optimum; eFD at the
  truncation floor
- detects: wrong FD stencil, wrong state perturbation
- type: unit - OPPROBE DIR/HS on surrogate.

**J004** JVP step-size sweep `[x]`
- property: error(h) has the expected U-shape
- reference: per-direction h sweeps
- tolerance: monotone then roundoff-limited curve; hstar recorded
- detects: systematic h mis-scaling
- type: unit - OPPROBE HS battery.

**J005** JVP additivity `[x]`
- property: J(u+v) ~ J u + J v
- reference: three independent JVPs
- tolerance: <= FD error level (d_op ~9e-8, d_pol ~1.9e-7 measured)
- detects: nonlinear contamination, state corruption between evals
- type: unit - T22/ADD rows.

**J006** JVP homogeneity `[x]`
- property: J(a v) = a J v
- tolerance: <= roundoff (HOM d_pol = 0 measured)
- type: unit - HOM rows.

**J007** JVP direction classes `[x]`
- property: populations, ED, T and mixed directions all probed
- tolerance: per-direction eFD documented (13 classes)
- type: unit - DIR rows on surrogate.

**J008** JVP constraint crossing + retry `[~]`
- property: constraint hit triggers bounded retry, not corruption
- tolerance: retry count bounded; JVP still within FD error
- type: unit - guard paths in tests/jfnk_core; production policy
  open.

**J009** residual repeatability `[x]`
- property: same state -> identical F
- tolerance: bitwise (REP/REPP = 0 on surrogate)
- type: unit - surrogate verified; production open (J011/J012).

**J010** state restoration `[~]`
- property: eval leaves every state component bit-identical
- reference: saved-state comparison over all mutated arrays
- tolerance: bitwise
- detects: state leaks (BUG-006 class defect)
- type: unit - BUG-006 ordering pin exists; a full-state compare
  driver is open.

**J011** frozen-J residual consistency `[~]`
- property: production frozen-J residual repeatable + consistent
- tolerance: measured REP/REPP on a live admission
- type: production-model probe (JFNK_HSCAN run) - P4-002 open.

**J012** coupled-J residual consistency `[ ]`
- property: same for the coupled tier (radiation sweep per eval)
- type: production-model probe - P4-003 open; depends on KNOP
  admission.

**J013** true GMRES residual `[x]`
- property: every done-exit measures ||M^-1(b-Ax)||
- reference: BUG-001 fix; unconditional true-residual measurement
- tolerance: acceptance only when true_res <= acc_fac*tol
- detects: hess_res-only acceptance
- type: unit - T4j/T4k/G17 regressions; production accounting open
  (P4-010).

**J014** left-preconditioned GMRES `[x]` - surrogate battery
  (weighted-norm consistency verified).
**J015** right-preconditioned GMRES `[x]` - surrogate battery.
**J016** preconditioner backward error `[x]` - T21 per-depth LU
  census; documented XFAILs structurally justified (P4-012).

**J017** singular/nearly singular blocks `[~]`
- property: solver behaviour on singular population blocks
- reference: T21 analysis (cond 1e15-1e27, ~3-dim null space)
- tolerance: loud failure or documented handling
- type: unit - structural evidence; dedicated solve open.

**J018** floor-pinned unknowns `[~]`
- property: n0=1e-30 floor columns do not corrupt corrections
- type: unit - T21 surrogate evidence; masked-row policy open.

**J019** restart and recycling `[~]`
- property: GMRES restart bounds; recycled subspace correctness
- type: unit - restart logic covered; recycling leg open.

**J020** hookstep trust-region `[x]`
- property: step stays on/inside the trust boundary
- type: unit - T1/T3 far-start cases.

**J021** pseudo-transient continuation `[~]`
- property: PTC sign and scaling
- type: unit - test_ptc_sign exists; production tuning open.

**J022** non-descent rejection `[x]`
- property: non-descent directions rejected before acceptance
- evidence: BUG-005 (155/155 m15 admissions defer correctly);
  guard tests.
- type: unit + production observation.

**J023** failed residual evaluation `[x]`
- property: eval failure propagates fail-closed
- type: unit - guards G13-G15.

**J024** NaN/infinity handling `[x]`
- property: non-finite residual/JVP fails closed
- type: unit - G17 + T4k biased cases fail on pre-fix code.

**J025** evaluation budget `[~]`
- property: eval cap enforced and visible
- type: unit - budget logic exists; dedicated cap test open.

**J026** checkpoint and restart `[ ]`
- property: JFNK state round trip (GMRES history, trust radius)
- type: unit - open.

**J027** Newton convergence rate `[x]`
- property: expected order on a manufactured nonlinear system
- reference: quadratic/quasi-quadratic rate on Rosenbrock/stiff
  problems (tests/plots t2_newton_rosen)
- type: unit.

**J028** production operator probe `[~]`
- property: OPPROBE battery on a live CMFGEN state
- expected: REP/REPP/DIR/HS/ADD/HOM rows measured at an admission
- type: production-model probe - instrumentation exists
  (`jfnk_op_probe_v1.f90`, `se_driver seprobe`,
  `test_jfnk_opint.py` 36/36 on surrogate); live run pending
  (P4-004/P4-005).

---

## 6. DIAG/TRI matrix-solve tests - T001..T019

Sources: `new_main/subs/solveba_v13.f`, `cmf_blkband_v3.f`,
`cmf_blkband_pin_null_v1.f`, `generate_full_matrix_v3.f`,
`subs/mat5pen.f`, `new_main/mod_subs/do_var_cont.f`.
Existing: `test_mat5pen.py` validates only the TRI Schur-complement
kernel; `test_cmf_blkband_pin.py` validates pin-null.  A
production-path manufactured test through CMF_BLKBAND_V3 needs the
SE-fabrication driver noted in P5-006 (STAT_EQ_DATA pointers + ATM
links + GENERATE_FULL_MATRIX_V3 inputs).  For every manufactured
system report: |r|=|Ax-b| rel/abs, backward error, estimated
conditioning, factorisation failures, pin-null activations.

**T001** DIAG block system, known solution `[ ]`
- input: manufactured per-depth blocks, x* known
- reference: independent assembly + numpy solve
- tolerance: |Ax-b| rel <= 1e-10 (well-conditioned blocks)
- detects: depth indexing, RHS misalignment
- type: unit (needs the SE-fabrication driver)

**T002** TRI block system, known solution `[ ]`
- input: manufactured block-tridiagonal system
- reference: dense reference solve of the assembled matrix
- tolerance: <= conditioning bound; residual rel <= 1e-8
- detects: band fill offsets, block-Thomas ordering errors
- type: unit

**T003** TRI reduces to DIAG when off-diagonals are zero `[ ]`
- input: TRI structure with B=D=0 blocks
- expected: TRI solution == DIAG solution bitwise-comparable
- tolerance: ~1e-14 rel
- detects: phantom coupling, uninitialized band reads
- type: unit

**T004/T005** one non-zero lower/upper off-diagonal block `[ ]`
- input: single coupling block at a chosen depth
- reference: dense solve
- tolerance: <= 1e-8 rel
- detects: which band slot carries which neighbour
- type: unit

**T006** both off-diagonal blocks non-zero `[ ]`
- input: full tridiagonal coupling
- reference: dense solve
- tolerance: <= 1e-8 rel
- detects: Schur-complement assembly errors
- type: unit

**T007/T008** boundary blocks at first/last depth `[ ]`
- input: coupling only at K=1 / K=ND edges
- expected: no sub-band at K=1, no super-band at K=ND
- reference: dense solve + band census
- tolerance: <= 1e-8 rel; boundary band census exact
- detects: OOB band access, boundary row errors
- type: unit

**T009** increasing ND `[ ]`
- input: ND = 8/16/32/64 systems
- expected: cost ~ linear, residual bounded
- tolerance: residual rel <= 1e-8 throughout
- detects: scaling/alloc issues
- type: unit

**T010** increasing block size `[ ]`
- input: block sizes spanning the production range
- tolerance: <= conditioning bound
- detects: leading-dimension errors
- type: unit

**T011** exactly singular block `[~]`
- input: a depth block with an exactly disconnected zero equation
- expected: pin-null path applies per documented policy
- reference: policy + dense pseudo-inverse residual bound
- tolerance: residual within documented pin semantics
- detects: silent wrong answer on a singular block
- type: unit
- status note: pin-null unit coverage exists
  (test_cmf_blkband_pin.py); the TRI-path singular case is open.

**T012** nearly singular block `[ ]`
- input: cond ~ 1e12-1e15 blocks
- expected: documented conditioning response; residual reported
- reference: condition estimate vs numpy
- tolerance: residual <= cond*eps*|b| bound
- detects: unreported ill-conditioning
- type: unit

**T013** pinned unknowns `[~]`
- input: masked/fixed rows in manufactured blocks
- expected: pinned rows hold their values
- type: unit - DIAG side covered by pin test; TRI open.

**T014** zero columns and zero rows `[ ]`
- input: blocks with zero rows/columns and live neighbours
- expected: pin policy preserves live equations and residuals
- type: unit - partially covered by pin test.

**T015** poorly scaled blocks `[ ]`
- input: blocks with entries spanning many decades
- expected: equilibration path behaves per DGEEQU design
- tolerance: residual <= scaled bound
- detects: scaling errors, overflow
- type: unit

**T016** non-finite input rejection `[ ]`
- input: NaN/Inf in blocks or RHS
- expected: loud rejection, not a silent pin or wrong solve
- type: unit

**T017** NUM_BNDS validation `[ ]`
- input: NUM_BNDS < 3 with SOL_METH=TRI
- expected: loud rejection/downgrade per P5-005 (VADAT warning +
  in-solver hard error + documented fill skip)
- tolerance: exact policy behaviour, never a silent reroute
- type: unit - partially covered by test_vadat_keys; in-solver
  guard coverage open.

**T018** DIAG/TRI agreement on depth-independent systems `[ ]`
- input: a physically depth-independent manufactured system
- expected: both methods reproduce the known solution (each solves
  its own declared system - do NOT validate TRI through DIAG)
- reference: dense solve of each declared system
- tolerance: each method <= its own residual bound
- detects: one method silently routed through the other
- type: unit

**T019** independent residual evaluation `[ ]`
- property: |A x - b| evaluated against an independently assembled
  full matrix for every manufactured case above
- expected: residual reported on every case (rel + abs + backward
  error + condition estimate + factorisation/pin census)
- tolerance: per-case bounds stated above
- detects: assembly-vs-solve mismatches the iterate hides
- type: unit - the umbrella requirement for this series.

---

## 7. Configuration test matrix `[~]`

Every supported combination of {transfer, J-path, solver, matrix}
must either dispatch correctly or reject loudly.  Status per cell:
S=supported+dispatched correctly (test), R=rejected loudly (test),
P=pending, D=declared-unsupported (documented, gate verified).

Baseline reference: `docs/audits/VADAT_SUPPORT_MATRIX.md` (32 cells
with file:line gate evidence, 2026-11-03) and `test_vadat_keys.py`
(parser-level admission/rejection coverage).  The runtime
verification (correct solver actually dispatched per config) is the
open leg for most cells: verify dispatch via OUTGEN solver-path
markers, never by assumed routing.

| # | transfer | J path | solver | matrix | status |
|---|----------|--------|--------|--------|--------|
| C01 | stock | LAM_ES=T | conventional | DIAG | P |
| C02 | stock | LAM_ES=T | conventional | TRI | P |
| C03 | stock | LAM_ES=F | conventional | DIAG | P |
| C04 | stock | LAM_ES=F | conventional | TRI | P |
| C05 | SOB | LAM_ES=T (FG_V13 J) | conventional | DIAG | S |
| C06 | SOB | LAM_ES=T | conventional | TRI | P |
| C07 | SOB | LAM_ES=F (mom J) | conventional | DIAG | S |
| C08 | SOB | LAM_ES=F | conventional | TRI | P |
| C09 | KNOP | LAM_ES=T (formal J) | conventional | DIAG | P |
| C10 | KNOP | LAM_ES=T | conventional | TRI | P |
| C11 | KNOP | LAM_ES=F (moment J) | conventional | DIAG | P |
| C12 | KNOP | LAM_ES=F | conventional | TRI | P |
| C13 | SOB | LAM_ES=T | JFNK (frozen tier) | DIAG | P |
| C14 | SOB | LAM_ES=T | JFNK | TRI | P |
| C15 | SOB | LAM_ES=F | JFNK | any | R (gated on USE_LAM_ES - test) |
| C16 | KNOP | LAM_ES=T | JFNK (coupled) | DIAG | P |
| C17 | KNOP | LAM_ES=F | JFNK (coupled) | DIAG | P |
| C18 | KNOP | any | JFNK | TRI | P |
| C19 | stock | any | JFNK | any | R (DO_JFNK needs a nonmono mode - test) |
| C20 | SOB+KNOP both | - | - | - | R (mutual exclusion - test) |
| C21 | any | any | any | TRI & NUM_BNDS<3 | R (loud downgrade - test) |
| C22 | SOB | LAM_SET=T | - | - | R (auto-normalized to F - test) |

Rejection coverage exists at parser level in `test_vadat_keys.py`
(mutual exclusion, DO_JFNK non-mono gate, JFNK_COUPLED KNOP-only,
GLOBAL_LINE/LAM_SET prerequisites, NUM_BNDS downgrade); runtime
dispatch verification is open for the supported cells.

---

## 8. Integration model suite - I001..I010 `[ ]`

Same atmospheric structure per applicable solver comparison; model
runs via `scripts/run_model_grid.py` on the cluster.  For each
completed model, compute independently (not via OUTGEN claims):
statistical-equilibrium residual, charge conservation,
radiative-equilibrium residual, luminosity conservation,
radiation-field consistency, population positivity/normalization,
electron-density consistency, temperature consistency, emergent
spectrum, convergence history.  A model is not validated because
OUTGEN reports convergence.  Compare under radial and frequency
refinement.

| ID | Case | Status |
|----|------|--------|
| I001 | monotonic reference atmosphere (flags off; P1-007 A/B) | [ ] |
| I002 | mild velocity reversal (V003-class structure) | [ ] |
| I003 | strong compression (V005-class; m15 SOB LAM_ES=T converged it718 - independent residual checks open) | [~] |
| I004 | multiple velocity reversals (V006-class) | [ ] |
| I005 | narrow shock (V007-class) | [ ] |
| I006 | broad shocked shell (V008/V010-class) | [ ] |
| I007 | STELLA-derived interacting-SN structure (INT-009) | [ ] |
| I008 | difficult ionisation structure | [ ] |
| I009 | restarted calculation (checkpoint round trip incl. FGINEXT_GS/MOM_JNEXT_GS) | [ ] |
| I010 | grid-refined calculation (radial + frequency) | [ ] |

Independent-check harness: extend `scripts/check_convergence.py`
plus a dedicated residual evaluator on the converged SCR state; all
ten identities reported per model with tolerances and failure
messages.

---

## 9. Pathway maturity summary (evidence ladder)

| pathway | stage | evidence |
|---------|-------|----------|
| velocity-gradient kernel (SIGMA_TRUE_CORE) | 4 | manufactured V001-V012 vs analytic + mirror, refinement rate, sweep battery (this file sec.1) |
| stock CMFGEN transfer (flags off) | 3 | bit-identity vs reference in sob/interp/edep suites; model A/B open (P1-007) |
| SOB line treatment (incl. tau floor) | 4 | floor formula bit-exact vs reference for sigma>=0, floored regime vs analytic (S003-S006); RH kernel pair-level (S017/S018) |
| KNOP formal solver | 3-4 | manufactured ray sweeps converge, multi-band + strong-GAM cases; analytic transfer references (KN002-004) open |
| KNOP moment solver | 3-4 | coefficient-level verification, resonance/doublet limits; production non-contraction documented (P3-012) |
| JFNK core (GMRES/Newton/PC) | 4 | surrogate manufactured battery incl. true-residual acceptance; live-state probe open (J011/J012/J028) |
| DIAG/TRI matrix solves | 3 | kernel-level (mat5pen, pin-null); production-path manufactured solve open (T-series) |
| configuration dispatch | 3-4 | parser-level admission/rejection matrix (VADAT_SUPPORT_MATRIX); runtime dispatch verification open |

## 10. Execution order and results log

| task | scope | status |
|------|-------|--------|
| TASK 1 | velocity-profile library + derivative/stencil validation | [x] 2026-11-04 - 103/103 checks; suite 23/23 files |
| TASK 2 | KNOP formal-transfer manufactured tests (KN series) | [ ] |
| TASK 3 | moment-solver tests incl. M011/M014 (M series) | [ ] |
| TASK 4 | Sobolev tests incl. multi-resonance (S series) | [ ] |
| TASK 5 | JFNK production-operator tests (J011/J012/J028) | [ ] |
| TASK 6 | DIAG/TRI manufactured matrix-solve test (T series) | [ ] |
| TASK 7 | cross-pillar integration suite (config matrix + I) | [ ] |
| TASK 8 | STELLA-model validation suite (I007+) | [ ] |

After each task: run new tests + relevant existing tests + full
suite; record exact results here and in
FIVE_PILLAR_IMPLEMENTATION_TODO.md; commit the completed task only.

### TASK 1 record (2026-11-04)

- New files: `unit_tests/velocity_profiles.py`,
  `unit_tests/python/test_velocity_profiles.py`.
- Driver: existing `unit_tests/f90/drv_sigma_true.f90`
  (`core`/`recomp` modes) - no new Fortran needed.
- Result: 103/103 checks PASS; full suite 23/23 files PASS
  (`python3 unit_tests/run_all_tests.py`, ifort, ~90 s).
- Evidence class: stage 4 (manufactured-system) for the
  velocity-gradient kernel; the library is the input set for
  TASKS 2-8.
