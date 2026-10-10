# Mathematical audit: dev7 KNOP solvers vs Knop et al. (2009) Papers I+II

Date: 2026-10 (post commits `3ed8f3b`, `cd45cf2`, `325d9a1`, `e2b9bd2`, `fc12e0f`).
Scope: `plane/fg_j_cmf_nonmono_v1.f`, `plane/mom_j_cmf_nonmono_v1.f`,
`new_main/mod_subs/do_var_cont.f`, `comp_j_blank.f`,
`solve_for_pops.f`, vs `docs/Knop_et_al_2009_formal_solution.pdf`
(Paper I), `docs/aa11793-09.pdf` (Paper II), and pristine
`cmfgen_reference/cur_cmf`.
Claims tagged VERIFIED / SUPPORTED / UNVERIFIED / CONTRADICTED /
BLOCKED per repo policy.

## 1. The equation dev7 actually solves

Knop Eq. 1 (lambda-form): `dI/ds = eta - chi*I - 4*a*I - a*d(lambda*I)/d(lambda)`.
CMFGEN solves the nu-form O(v/c) comoving equation

  `dI/ds = eta - chi*I + a_nu * dI/dln(nu)`

with `a_nu = GAM = (V/R)(1 + sigma*mu^2)/c` per ray point.  VERIFIED:
stock `fg_j_cmf_v13.f:1408-1412` builds `chi_hat = chi + GAM/dLOG_NU`,
`Q = GAM/(dLOG_NU*chi_hat)`, `SOURCE_PRIME = eta/chi_hat + Q*I_prev`
-- the same three-line structure dev7 generalises.  The `4aI`
aberration term (which becomes `+3*a_nu*I` in nu-form) is absent in
stock CMFGEN and absent in dev7; the dev7 comment at
fg_j_cmf_nonmono_v1.f:2440-2447 documents this as deliberate.
VERIFIED: consistent with the equation CMFGEN solves; deviation from
Knop's literal Eq. 1 is inherited from upstream, not introduced by
dev7.  Adding `3*a_nu*I` would make dev7 differ from stock V13, not
improve Knop fidelity.

## 2. Paper I discretization vs dev7 (VERIFIED correct)

dev7 uses xi=1 throughout.  With xi=1 the `[1-xi]` terms of Knop
Eq. 4/10 vanish, leaving only the aberration part of S_tilde (already
absent in CMFGEN's equation).

| Paper I | dev7 | Status |
|---|---|---|
| Eq. 7: `chi_hat = chi + xi*a*p|` (always > 0 since `a*p| > 0`) | `chi_hat = chi + |GAM|/DLN_C`, floored 1e-30 | VERIFIED identical (nu-form).  `a*p|` is positive for both signs of `a` because `p|` changes sign with `a` -- dev7's `|GAM|` realises exactly that. |
| Eq. 9: `S_hat = (chi/chi_hat)(S - (a/chi)(p- I_{l-1} + p+ I_{l+1}))` | `SOURCE_PRIME = eta/chi_hat + Q*I_up`, `Q = GAM/(DLN*chi_hat)` signed expanding, `ABS(Q)` compressing | VERIFIED.  For `a>=0`, `p- < 0` gives source `+|a*p|*I_{l-1}`; for `a<0`, `p+ > 0` gives `+|a*p|*I_{l+1}`.  Both branches produce a positive source contribution at the correct upwind. |
| Eq. 3: sign of `a` selects upwind; `p` depends on position and characteristic | keyed on `GAM(I,LS)` (per ray point, mu-dependent) | VERIFIED.  The sign test is on the ray-local coupling, matching the paper's requirement that `p` varies along the characteristic. |
| `p` denominators use the upwind-side spacing | `DLN_C = dLOG_NU` expanding; `DLN_N_SAV = ln(nu_k/nu_{k+1})` compressing (Oct-2026 fix; `SET_FG_NU_NEXT` per call) | VERIFIED correct and necessary at ~10x line-edge spacing jumps. |
| `p|` vanishes at grid boundaries | blue edge: freq 1 runs INIT (Q=0).  red edge: slot `I_NEXT(:,:,NCF)` stays 0 -> `SOURCE_PRIME = S` while `chi_hat` keeps `|GAM|/dln` | Zero-incident wavelength BC at the red edge: photons advected below the grid leave it.  Defensible, documented here; it is the only boundary where `chi_hat`'s absorption lacks a compensating source term. |

Additional verified consistency:

- `d(chi_hat)/dR = dchi/dR + SIGN(1,GAM)*dGAM/dR/DLN` -- correct
  `d|GAM|/dR` for the monotonic-interpolation DTAU (line 2569).
- No negative generalized opacity can occur: `chi_hat >= chi` with
  the velocity term non-negative by construction.
- The formal solution then proceeds identically to V13 (DTAU,
  E-functions, A0-A4, monotonic slopes, same boundary options
  DIFFUSION/PNT_SRCE/ZERO_FLUX) -- spatial discretization untouched.
- I_NEXT read uses linear interpolation between the two bracketing
  compressing depths (KC1/KC2, weight W_INEXT) -- avoids a spatial
  discontinuity at strong velocity gradients.

## 3. Paper II iteration structure vs dev7

Paper II Sect. 3.2 converges the wavelength boundary-value problem
for fixed source before each ALI/source update: Eqs. 25-28 are a
Gauss-Seidel march ordered along the physical information flow
(quasi-analytic for linear S).  dev7's structure:

- Expanding cells: in-sweep reads `I_M_PREV`/`I_P_PREV` (this sweep's
  k-1) -- live GS, identical to Knop's `(n+1)` reads.  VERIFIED.
- Compressing cells: reads the lagged store `I_P/I_M_DEPTH_NEXT`
  (previous great iteration) -- a Jacobi outer iteration, then
  converged by `FG_KNOP_REFRESH_V1`: alternating-direction re-solves
  on captured (frozen) opacities.  `BUILD_FREQ`/`SRC_PRIME`
  (fg_j_cmf_nonmono_v1.f:5190-5359) replicate the in-sweep `chi_hat`,
  `Q`, spacing rule, `|Q|` compressing branch, DTAU/A-coefficients and
  boundary guards bit-for-bit -- VERIFIED the refresh solves the same
  discrete system, so its fixed point is identical.  The alternating
  passes recover the paper's information-flow ordering (its Fig. 2):
  descending passes make `I_{l+1}` reads fresh, ascending passes make
  `I_{l-1}` reads fresh.
- Under-relaxed Jacobi store write (`INEXT_RELAX`, cap
  `NONMONO_INEXT_CAP`, non-finite floor) -- pure convergence
  machinery; fixed point unchanged.
- Maturity gate in `solve_for_pops.f` (warm-up + GS_FREEZE/
  FORCE_LAMBDA on sweep-drift metrics) enforces the paper's ordering
  "converged BVP -> SE update" at the great-iteration level.

DEVIATION (documented, not a bug): the paper's ALI uses the full
spatial-bandwidth Lambda* tridiagonal-in-wavelength operator (Eqs.
29-31); dev7 uses a diagonal OAB Lambda* and merges the wavelength GS
into the CMFGEN great-iteration loop.  Same fixed-point equations,
different iteration machinery.  Domain decomposition of Paper II
Sect. 3.1 is replaced by OpenMP ray parallelism (characteristics are
independent, matching the paper) plus the sparse
`COMP_DEPTH_IDX` compressing-depth store.

## 4. Diagonal Lambda* (ALO) -- VERIFIED correct

`ALO_DJDETA`/`ALO_DJDCHI` accumulation (fg_j_cmf_nonmono_v1.f:3527-3558):

- `dAV/dS' = 0.5*(A2(K-1) + A2(K))` with correct boundary cases
  (core-ray DIFFUSION/PNT_SRCE inner BC drops `dI_P`; ZERO_FLUX and
  non-core rays reflect) accumulated with the same `JQW` weights and
  `J_PNT` mapping as `JNU_STORE` -- it is the true diagonal
  sensitivity of the computed J.
- `dJ/dETA = dAV/dS' * 1/chi_hat` and
  `dJ/dCHI = -dAV/dS' * S'/chi_hat` (split into `S'_P`/`S'_M` since
  the two directions have different upwind terms) -- exactly
  `dS'/deta = 1/chi_hat`, `dS'/dchi = -S'/chi_hat` of
  `S' = (eta + a*p*I_up)/chi_hat`, including the coupling term's
  `1/chi_hat` dependence.  This is the nu-form of Paper I's
  `(chi/chi_hat)` weighting propagated into the Lambda operator, and
  matches Paper II's "construction proceeds exactly as in Baron &
  Hauschildt (2004)".
- In `do_var_cont.f:362-419` the KNOP branch zeroes `TX/TVX`
  off-diagonals and fills `TX(K,K,1)=dJ/dCHI`, `TX(K,K,2)=dJ/dETA` --
  index convention VERIFIED against `var_mom_j_cmf_v12.f` header
  ("first 2 represent dJ/dCHI and the second dJ/dETA").  For
  `USE_LAM_ES=F` the MOM diagonal is used instead, and quarantined
  cells get zero response (their J is frozen, so zero sensitivity is
  correct).

## 5. MOM solver vs papers and V11 (VERIFIED, with caveats)

`MOM_J_CMF_NONMONO_V1` is not Knop's characteristic solver; it is
V11's variable-Eddington moment solve generalised to sign-changing
frequency coupling.  The Knop technique is applied per equation:

- Effective coupling per equation: J-equation `a_J = CON_GAM*(1 +
  sigma_true*K_ON_J)` at depth points; H-elimination
  `a_H = CON_GAMH*(1 + sigma_true*NMID_ON_HMID)` at midpoints.
  Sign tests `Q_NM(I)<0` / `GAMH_SIGNED_NM(I)<0` are per-row -- the
  correct moment-level analogue of Knop's per-point `a`.
- `chi_hat(midpoint) = chi_avg + |a_H|/DLN`,
  `chi_hat(depth) = chi + |a_J|/DLN` -- Knop Eq. 7 per equation;
  DTAU then uses `chi_hat`, so no negative optical depth.  VERIFIED.
- `PSI = DTAUONQ*|Q|` (positive implicit diagonal), `PSIPREV` signed
  explicit upwind; compressing midpoints fix `HS = W/(1+W)` and
  `EPS_PREV = EPS` so the explicit weight equals the implicit weight
  (non-expansive lagged map -- the same `|Q|` convention `PSI_NM`
  uses; the Sep-19 expansive-map bug was exactly a weight mismatch).
  UPDATE (30-Sep-2026): the expanding-branch `PSIPREV` used the
  signed lagged `GAM_FAC_PREV` factor, so a `1+sigma*K_on_J` sign
  flip across adjacent frequencies could inject a negative upwind
  source.  Fixed to `|GAM*GAM_FAC_PREV|` at all three assignment
  sites; Knop I Eq. 9 requires that coefficient >= 0 by
  construction, so "signed" above was wrong -- the upwind weight is
  positive, the *sign choice* picks which neighbour it couples.
- Upwind fields `JUP_NM`/`HUP_NM`: expanding read the k-1 chain
  (live); compressing read lagged `J_NEXT_FRQ`/`H_NEXT_FRQ`.  All
  four read sites (interior RHS, `OUT_BC_TYPE>1` row, negative-J
  rebuild, H reconstruction) use the same direction-aware fields --
  VERIFIED consistent.
- EPS scattering: per-midpoint `+2E*J_up` to row m, `-2E*J_up` to row
  m+1 -- algebraically identical to the implicit stencil, now with
  the direction-aware `J_up`.  VERIFIED.
- Frequency boundaries: same zero-incident semantics as the FG path
  (reddest slot stays 0); inner/outer spatial BCs are V11's, unchanged.
- Negative-J repair (`XM_FB`): same-frequency lagged `J_NEXT` where
  positive, else local `SOURCE*R^2`.  The `RSQJNU_PREV` (bluer-field)
  fallback was removed (`cd45cf2`) -- it had injected the wrong
  wavelength into failed cells and produced the 2.7e20 it1 ratchet.
  No remaining cross-frequency fallback was found (grep-verified all
  `RSQJNU_PREV`/`RSQHNU_PREV` uses: they only appear on expanding
  branches and the lagged-field chain, where bluer is the correct
  upwind).

## 6. Defects found and fixed during the audit arc

| Commit | Defect | Status |
|---|---|---|
| `fc12e0f` | Makefile order made clean builds depend on stale `.mod` files; `cmfgen` target linked an unrelated LTE tool | fixed, verified clean build |
| `e2b9bd2` | `DLN_C` shared across the OMP ray loop (chi_hat/Q spacing race) | fixed, now PRIVATE |
| `3ed8f3b` | lagged H/EPS terms always read bluer frequency at compressing midpoints (wrong upwind); `H_NEXT_FRQ` store added, `HS`/`EPS_PREV` fixup, checkpoint VER 4 | fixed |
| `cd45cf2` | `XM_FB` negative-J repair fell back to `RSQJNU_PREV` (bluer field) on cold start -> cross-frequency ratchet, it1 L=2.7e20 | fixed; it1 L=2.64e10 after fix |
| `325d9a1` | per-cell negJ diagnostic unbounded (37 GB in one iteration) | capped 20k rows/sweep; counters unchanged |

## 7. Known deviations and limitations (not bugs)

1. Aberration/advection `3*a_nu*I` term absent -- inherited from
   CMFGEN's nu-form equation (stock V13 identical).  Deliberate.
2. nu-form log-spaced `p` coefficients absorb the
   `lambda_l/lambda_{l+1}` asymmetry into `dln(nu)` -- same convention
   as V13; bounded by the local spacing ratio (~5-15 % typical, more
   at line-edge jumps where the DLN_N_SAV fix applies).
3. Diagonal Lambda* only; ALI/pops iteration is diagonal quasi-Newton
   with zeroed `TX/TVX` off-diagonals -- documented structural choice;
   convergence relies on the outer J<->n map.
4. Red-edge zero-incident wavelength BC: compressing cells at k=NCF
   (or with no I_NEXT slot) keep `|GAM|/dln` in `chi_hat` without a
   compensating source.  Photons advected below the grid edge are
   absorbed, not re-injected -- a boundary-condition choice, standard
   in spirit but the only place absorption lacks its source partner.
5. `NONMONO_INEXT_CAP` = 1e20 on stored I_NEXT -- numerical overflow
   guard, not physics; it fired only on an earlier divergent run.
6. MOM refresh lagged-map residual still contracts without reaching
   the convergence test on this model (~5.9e3 and declining); the
   retained-field rollback keeps the converging sweep field, so the
   production path is unaffected -- intrinsic moment-map property,
   consistent with prior findings in [[convergence_issues]] Sect. 1.2.

## 8. Runtime evidence at audit time

- MOM `12000120` (500-it cap, small): it2 L=1.96e10, capped/NaN
  writes=0, refresh pass counts contracting (cells>5% 5.87M -> 3.57M
  over 9 passes), retained-field delta 0.  Consistent with the fixed
  binary; earlier `2.88e20` it1 signature gone.
- LAM `12000122`: in iteration 1 (~400 s/iter); prior interactive run
  showed contraction 5.9e16 -> 864 over 3 its, refresh converged
  (unconverged cells=0).

## 9. Verdict

- FG/LAM path (`USE_LAM_ES=T`): Knop Paper I discretization VERIFIED
  correct in nu-form; Paper II iteration structure VERIFIED at
  great-iteration level with documented differences (diagonal ALO,
  merged GS).  Integration into do_var_cont/solve_for_pops/comp_j_blank
  is minimal-diff and flag-gated.  PASS WITH LIMITATIONS (Sect. 7
  items 1-5).
- MOM path (`USE_LAM_ES=F`): moment-analogue of the Knop upwind
  scheme VERIFIED consistent after the direction-aware and XM_FB
  fixes; per-equation sign keying and non-expansive explicit/implicit
  weight equality are in place.  Convergence of the internal refresh
  map on strongly compressing models remains the open numerical
  question (Sect. 7 item 6); sweep-level convergence and rollback
  containment are verified.  PASS WITH LIMITATIONS.
- No remaining mathematical defect was found in the discretization
  itself: every verified numerical failure in the run history traced
  to a plumbing bug (wrong-frequency fallback, shared OMP variable,
  diagnostic I/O) or to model/setup physics (LTE start at T_gas,
  intrinsic moment-map stiffness) -- not to the Knop equations.
  CORRECTION (30-Sep-2026): a re-review found two defects in the MOM
  discretization -- a dead midpoint CHI_HAT loop reading
  uninitialized/stale GAM_MAG_NM, and the signed-PSIPREV positivity
  violation noted above.  Both fixed; the verdict stands for the
  remaining terms.
