# CMF_FLUX P6 audit - observer-frame flux normalisation and angular quadrature

Scope: `obs/obs_frame_sub_v9.f`, caller `obs/cmf_flux_sub_v5.f`, quadrature
helpers `subs/hweight_v2.f`, `subs/set_acc_dmu.f`, `obs/revise_obs_p.f`.
Reference tree `cmfgen_reference/cur_cmf` carries the identical
normalisation block; the audit below is an independent verification, not a
diff.

## 1. The complete flux equation

Emergent ray intensity along impact parameter `p` at observer frequency
`nu_obs`, per ray `LS` and output bin `ML`:

    I(p,nu_obs) = I_IB * exp(-tau_tot)
                + int_ray  j_obs(nu_obs) * exp(-tau(z)) dz          (STAU)
    or the equivalent source-function formal integral                 (ETAZ)

with

    j_obs(nu_obs)  = j_CMF(nu_CMF(I)) * NU_ON_NU0(I)^2     (j_nu/nu^2 invariant)
    chi_obs(nu_obs)= chi_CMF(nu_CMF(I)) / NU_ON_NU0(I)     (nu*alpha_nu invariant)
    I_IB,obs       = I_IB,CMF * NU_ON_NU0^3                (I_nu/nu^3 invariant)
    nu_CMF(I)      = nu_obs * (1 - VMU(I)) * GAM(I)
    NU_ON_NU0(I)   = nu_obs/nu_CMF = 1/GAM/(1-VMU)

`PAR_FLUX` at the accumulation point (obs_frame_sub_v9.f :1608) is the
emergent intensity `I(p_LS, nu_obs)` in cgs units
(erg s-1 cm-2 Hz-1 sr-1).  It is not frequency-integrated and carries no
geometric factors: both boundary seed and integrand are intensities.

The relativistic boost is applied exactly once, factorised as
j -> nu^2, alpha -> nu^-1 (folded into ETA_VEC/CHI_VEC per ray point,
:1089-1090, :1292-1293), so S = ETA/CHI carries nu^3 and the boundary
term multiplies NU_ON_NU0^3 directly (:1513, :1574).  Under
DO_REL_CORRECTIONS=F all NU_ON_NU0 are 1.

Then

    OBS_FLUX(ML) = sum_LS HQW_AT_RMAX(LS) * I(p_LS,nu_obs)
                 ~= int_0^1 I(mu) mu dmu
    OBS_FLUX    *= 6.59934 * R(1)^2          ! Jy at 1 kpc

## 2. HQW_AT_RMAX and the quadrature

Comment at :1738: "HQW_AT_RMAX is MU dMU and PdP is R^2 MU dMU" - the
weights integrate V(mu) against mu dmu, i.e. sum W*V ~= int_0^1 V mu dmu.
Construction, per geometry branch in cmf_flux_sub_v5.f (:2611-2646, set
before both V9 call sites at :2688 and :3108):

- PLANE_PARALLEL: GAULEG on [0,1] over NC angles,
  W = w_gauss * mu * (R(ND)/R(1))^2.  The (RND/R1)^2 factor maps the
  unit-disk mu integral onto the slab of radius R(ND); for the standard
  PP convention R(1)=R(ND) it is 1.
- REVISE_P_GRID (non-mono default): REVISE_OBS_P builds P_OBS on
  [0, R(1)] monotone increasing, P_OBS(1)=0 (mu=1 first), the core
  boundary P=R(ND) counted once at index NC_OBS+1, NC_OBS=NC*(NC_INS+1)
  core rays uniform in mu_core, and NP_INS extra rays per shell.
  MU = sqrt(1-(P/R1)^2), dMU from SET_ACC_dMU, W from HWEIGHT_V2.
- Standard: same HWEIGHT_V2 on the caller's P grid.

HWEIGHT_V2 is the Nordlund modified-cubic first-moment rule (subs/
hweight_v2.f).  It returns W such that sum W(I)*V(I) ~= int_0^1 V mu dmu;
its internal check requires sum W*X = 1/3 (V=mu exact to 1e-12).  The
X(N)=0 limb case is handled analytically (V(0)=0 assumed); a grid ending
short of the limb extrapolates with a printed warning.

SET_ACC_dMU computes dmu = mu(i)-mu(i+1) directly, switching to the
binomial expansion of 1-sqrt(1-x^2) when x<0.999 fails - exact where
mu->1 (P->0), the numerically dangerous corner.

Ray counting: each LS accumulates exactly one term at :1608; hollow-core
back-side sweeps fold into the same PAR_FLUX before accumulation
(BACK_SIDE -> GOTO 35000, then a single +=).  No ray is counted twice or
skipped.  HQW depends only on P/R geometry - the velocity law cannot
alter the quadrature.

## 3. The Jy conversion, derived

At the outer boundary:

    F_nu(R1)   = 2 pi int_0^1 I mu dmu = (2 pi/R1^2) int_0^R1 I p dp
    F_nu(D)    = F_nu(R1) (R1/D)^2
    F_Jy(D)    = 1e23 * 2 pi (R1/D)^2 int I mu dmu

with R1 in units of 1e10 cm and D = 1 kpc = 3.0856e21 cm:

    F_Jy = [1e23 * 2 pi * 1e20 / (3.0856e21)^2] * R1^2 int I mu dmu
         = 6.5993417128 * R1^2 int I mu dmu

The code's 6.59934 agrees to 4 s.f. (rel. error 2.6e-7) - a deliberate
upstream precision choice, not a defect.

Unit chain: R,P,dZ in 1e10 cm; CHI and ETA share the same per-1e10-cm
scale factor, so tau = chi dz is dimensionless, S = ETA/CHI is cgs
intensity, and int ETA e^-tau dz is cgs intensity.  Independent
confirmation: the boundary seed PAR_FLUX = J + 3 mu H (cgs moments from
the formal solve) is added to ETA_VEC with no conversion (:1511, :1572).

Cross-check on a uniformly emitting optically thin sphere (RND core,
chord = 2 sqrt(R1^2-p^2) for p>RND, shell chord for p<=RND):

    int I p dp = (2 j/3)(R1^3 - RND^3)
    F_Jy = 6.59934 * (2j/3) * dR^3 = 4.3996 * j * dR^3

vs the luminosity route L_nu = 4 pi j (4 pi/3) dR^3:

    F_Jy = 1e23 * L/(4 pi D^2) = j dR^3 * (4 pi/3) 1e53 / 9.5209e42
         = 4.3996e10 * j_cgs * dR^3      (j stored as j_cgs*1e10)

identical.  The angular factor, geometric dilution, and Jy scale are
mutually consistent.

## 4. Distance convention

F_nu is output at D = 1 kpc - a fixed convention in every
OBS_FRAME_SUB_V* (v2-v9 all carry 6.59934*R(1)^2).  No distance control
exists upstream or in the controls module; rescaling to other distances
is the user's (1kpc/D)^2.  OBS_FREQ is in 1e15 Hz; OBS_FLUX is a
per-Hertz density, so the frequency labelling does not affect the
conversion.

## 5. Findings

1. Quadrature, geometry, boost factorisation, and the Jy constant all
   verify analytically and against the luminosity route.  No defect in
   the normalisation chain.
2. `dI_R` (WRITE_dFR) scales by the same 6.59934*R1^2 and integrates to
   OBS_FLUX minus the undistributed boundary seed - consistent within
   its diagnostic scope; RPHOT/VPHOT derive from it.
3. IP_DATA writes PAR_FLUX raw (cgs intensity, no Jy factor) - a
   diagnostic convention to document.
4. Minor, documented-not-changed: HWEIGHT_V2's extrapolation branch
   triggers only when the P grid ends short of R(1); the revised grid
   always ends at R(1) so this is dormant for non-mono runs.
5. PP-branch (RND/R1)^2 factor assumes R(1)=R(ND); harmless under the
   standard PP convention.
6. HWEIGHT_V2 is not an exact rule for V=const: sum W = 0.49960 at
   NP=33 and 0.49997 at NP=129 (O(h^2) toward 1/2), while sum W*mu =
   1/3 is enforced exactly.  This is the documented Nordlund design
   (V ~ mu at the limb), not a defect - but constant-intensity
   comparisons must use resolution-scaled tolerances, not machine
   precision.
7. Inner-boundary conventions differ by geometry and method, all
   physically consistent:
   - solid core (HOLLOW_CORE=F, V(ND)<=100): STAU seeds PAR_FLUX with
     the local source function at the deepest point - the thick-wall
     boundary condition.  A solid wall therefore emits at S and can
     dominate the small core-disk solid angle even when the envelope
     is thin.
   - hollow core (V(ND)>100): back-side sweep computes the incident
     intensity from the far hemisphere; no wall term.
   - ETAZ seeds PAR_FLUX=0 for the solid wall - the wall does not
     emit under the direct emissivity method.  This is an upstream
     convention difference between the two methods, not a bug.
8. ETAZ surface-cell resolution limit: the emissivity integral seeds
   TAU(1)=CHI(1)*HALF_DZ(2) (:1531), i.e. the surface point is treated
   as half a cell deep.  When the outermost ray cell is optically
   thick the whole integrand is suppressed by exp(-chi*dz/2) and the
   emergent flux is under-integrated (measured: ~2e4x low for a model
   with chi*dz~50).  STAU is unaffected because it works in source
   functions.  Production grids keep the outermost cells thin, so
   this is a documented resolution requirement, not a code defect;
   the same limit exists upstream (identical TAU(1) seed).
9. Relativistic boost: j -> j*D^2, chi -> chi*X, boundary I -> I*D^3
   applied pointwise, once, matching upstream obs_frame_sub_v5/v9
   byte-for-byte.  The mirror pairing (+/-VMU) cancels the O(beta)
   term; the residual scales as beta^2 (measured ratio 4.01 between
   v0=40000 and 20000 km/s on the CHI=0 pure-boost probe).  The
   beta^2 coefficient is (3 mu^2-1), which nearly cancels for a thin
   shell geometry; the ~1e-6*beta^2 residual vs a continuous-model
   prediction is consistent with discretisation of the D^2 variation
   on the subdivided ray grid, not a convention error.

## 6. Tests (drv_obs_frame_norm.f / test_obs_frame_norm.py)

- constant-I: thick sphere (chi=100, S=1) with solid wall at V(ND)=0;
  STAU flux matches 6.59934*R1^2*int I mu dmu to quadrature error;
  ETAZ run pins the section-5.8 resolution limit.
- limb law: thin hollow-core shell (V(ND)=300>100), mirrored chords
  I(p)=2j(sqrt(R1^2-p^2)-sqrt(RND^2-p^2)) inside the core radius ->
  OBS_FLUX = 6.59934*(2j/3)(R1^3-RND^3); NP series 9..129 for the
  angular convergence study; both STAU and ETAZ at NP=129.
- quadrature-only: production SET_ACC_dMU + HWEIGHT_V2 weights dumped
  at NP=33/129 and checked against int V mu dmu for V=1, mu,
  1-a+a mu (a=0.3,0.8), mu^2 with Nordlund-scaled tolerances, plus
  strict improvement with NP.
- V-independence: flat ETA/CHI with V(ND)>100 fields (uniform 200,
  monotonic, reversing, negative) agree to round-off; the solid-wall
  pair (uniform 50 vs reversing with V(ND)=50) also agrees - two
  physically distinct geometry classes, equality required only within
  each.
- rel boost: thin hollow shell, uniform V=20000/40000 km/s.  The
  CHI=0 ETAZ probe isolates the pure emissivity boost and scales as
  beta^2; the full rel flux matches an independent
  int j [D(z)^2+D(-z)^2] dz prediction to quadrature level.
- stock comparison: the same driver objects are linked once with the
  tree obs_frame_sub_v9 (obs_v9_adapter.f) and once with the frozen
  upstream file (obs_v9_adapter_ref.f, which drops the
  inner-boundary arguments); every flux_*.txt output is compared.
  Measured: bitwise identical on all 24 manufactured outputs - the
  non-mono extensions, inner-boundary arguments, and P5 coverage
  guards do not alter the arithmetic on these cases.
  Harness note: the upstream-side binary crashes inside the Intel
  OpenMP runtime teardown after STOP (all outputs are written and
  closed by then).  The test runs it with check=False and verifies
  the output set instead.
