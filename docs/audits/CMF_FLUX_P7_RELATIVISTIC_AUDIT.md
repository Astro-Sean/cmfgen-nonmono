# CMF_FLUX P7 - Relativistic Transfer Transformations Audit

Audit of the relativistic transformation chain in
`obs/obs_frame_sub_v9.f`, the observer-frame formal-integral routine
used for the CMFGEN synthetic spectrum.  The chain covers the Doppler
factor convention, the Lorentz factor construction, the three Lorentz
invariant transformations (emissivity, opacity, intensity), angular
aberration at the diffusion inner boundary, and rejection of
physically invalid velocities.

Scope: `obs/obs_frame_sub_v9.f` only, plus its callers for provenance.
No KNOP, population-solver, Sobolev-deposition, electron-scattering or
resolution code was touched.

## 1. Confirmed defects

Three defects were confirmed and fixed in this phase.

### 1.1 Superluminal velocities were floored into a finite spectrum

Old code (`obs_frame_sub_v9.f`, gamma loop, was at ~:664):

```
T2=1.0_LDP-V_RAY(I)*V_RAY(I)/C_KMS/C_KMS
IF(T2 .LE. 0.0_LDP)THEN
  WRITE(LUER,...) 'Warning OBS_FRAME_SUB_V9: V_RAY >= c ...'
  T2=1.0E-12_LDP
END IF
GAM_RAY(I)=1.0_LDP/SQRT(T2)
```

A superluminal interpolated velocity produced `GAM=1e6` and the formal
solution continued.  The warning goes to LUER amid thousands of ray
progress lines; nothing downstream re-checks the value, so a
physically impossible velocity field could yield a finite, bogus
spectrum.  Replaced with fail-closed validation (section 5).

### 1.2 Inner-boundary seed used the observer-frame direction cosine

Old code (STAU branch ~:1510, ETAZ branch ~:1571):

```
T1=Z_RAY(SM_NRAY)/R_RAY(SM_NRAY)
PAR_FLUX=T3+3.0_LDP*T1*T4          ! T3=IB_J_CMF, T4=IB_H_CMF
```

`IB_J_CMF` and `IB_H_CMF` are comoving moments (`IB_J_CMF_ST(ML)=RJ(ND)`,
`IB_H_CMF_ST(ML)=HFLUX_AT_IB` in `obs/cmf_flux_sub_v5.f`).  The
Eddington intensity law `I(mu) = J + 3*mu*H` is a comoving-frame
relation, so its angular argument must be the comoving direction
cosine, not `mu_obs = Z/R`.  The difference is O(beta): for a boundary
flow at beta the old code evaluated the seed at the wrong angle.
Corrected at both sites (section 4).

### 1.3 Plane-parallel boundary cosine was not a direction cosine

In `PLANE_PARALLEL` geometry `Z(I)=R(I)/MU` is the slant depth, so the
seed's `Z_RAY(SM_NRAY)/R_RAY(SM_NRAY)` evaluates to `1/MU` at the
inner boundary - a value in `(1,inf)`, not a direction cosine.  That
bogus value enters `J + 3*mu*H` directly on the non-relativistic path
(e.g. `PLANE_PARALLEL_NO_V` runs) and additionally passes through the
aberration under `DO_REL_CORRECTIONS`, where `1 - beta*mu` can hit
zero or produce `mu_cmf` well outside `[-1,1]`.  Production reaches
the branch via `(DO_NONMONO_SOB .or. DO_NONMONO_KNOP) .and. DIF` on a
plane-parallel model (`PLANE_PARALLEL` or `PLANE_PARALLEL_NO_V` - a
geometry property, e.g. `rev_hydro_turb.f`, not one of the solver
modes rejected under the nonmono flags).  Fixed at both sites:
`MU_AT_RMAX(LS)` is used as the boundary cosine in PP mode
(section 4.5).

## 2. Frame conventions (Part 1)

### 2.1 Quantity table

| Quantity | Frame | Units | Sign | Construction | Meaning |
|---|---|---|---|---|---|
| `V_RAY(I)` | comoving->lab | km/s | signed, + outward | interpolated from `V(R)` by `MON_INTERP_FAST` at `R_RAY(I)` (~:617) | signed radial flow speed at the ray point |
| `R_RAY(I)` | lab | 1e10 cm | positive | `sqrt(P(LS)^2 + Z_RAY(I)^2)` spherical (~:615) | spherical radius of the ray point |
| `Z_RAY(I)` | lab | 1e10 cm | signed, + toward observer | ray grid construction; negated on the mirror segment (~:637-656) | position along the ray |
| `VMU_RAY(I)` | lab | dimensionless | signed | `V_RAY*Z_RAY/R_RAY/C_KMS` spherical; `V_RAY*MU_AT_RMAX(LS)/C_KMS` plane-parallel (~:633-635) | `beta*mu_obs`, projected speed in units of c |
| `GAM_RAY(I)` | lab | dimensionless | >= 1 | `1/sqrt(1 - V_RAY^2/C_KMS^2)` (~:745) | Lorentz factor of the TOTAL radial speed |
| `NU_ON_NU0_RAY(I)` | lab | dimensionless | positive | `1/GAM_RAY/(1 - VMU_RAY)` (~:748) | `D = nu_obs/nu_cmf` |
| `mu_obs` | lab | dimensionless | signed | `Z_RAY(I)/R_RAY(I)` | cosine between photon direction and local radial direction |
| `beta` | lab | dimensionless | signed | `V_RAY/C_KMS` | total radial speed in units of c |

`C_KMS = SPEED_OF_LIGHT()*1.0E-05` at ~:251 converts the exact cgs
value `2.99792458e10 cm/s` to `299792.458 km/s`.

Mirroring (~:637-656): hollow-core and back-side ray segments negate
`Z_RAY` and `VMU_RAY` while keeping the signed radial `V_RAY`, which
is exactly what the projected-velocity definition requires: on the
far side `mu_obs<0` so the same positive `V_RAY` gives a receding
(redshifted) projected speed.

Plane-parallel: `VMU_RAY = V_RAY*MU_AT_RMAX(LS)/C_KMS` uses the fixed
direction cosine of the slab normalisation, equivalent to
`beta*mu_obs` with a constant `mu_obs`.  Same convention, different
geometry.

### 2.2 Doppler derivation from the photon four-momentum

Photon four-momentum `p^a = (h*nu/c)(1, n)` with `n` the propagation
unit vector.  The gas moves radially with velocity `v = beta*c`; the
comoving frame is boosted by `-v` along r-hat.  The frequency measured
in the comoving frame is

```
nu_cmf = nu_obs * gamma * (1 - beta * (n . r-hat))
       = nu_obs * gamma * (1 - beta * mu_obs)
```

so

```
D = nu_obs / nu_cmf = 1 / (gamma * (1 - beta*mu_obs)) .
```

For a photon heading toward the observer on the near side,
`mu_obs = Z/R > 0`: outward flow (`beta>0`) gives `1-beta*mu<1`, hence
`nu_cmf < nu_obs` - the photon is blueshifted relative to the emitter,
the standard result.  On the far side `mu_obs<0` and the same outward
flow redshifts.  Sign flip comes entirely through `mu_obs`; `beta`
keeps the sign of the radial velocity.  This is the convention in
Mihalas & Mihalas (1984, eq. 14.1-14.2 family, photon direction
cosine transformation) and in the Jerkstrand notes (Doppler factor for
a moving emitter).

`GAM_RAY` is built from `V_RAY*V_RAY` - the total radial speed - not
`VMU_RAY`.  Verified directly at ~:745.  Using the projected speed in
`gamma` would break the transverse-Doppler term (the `gamma` factor at
`mu_obs=0`); the code does not make that error.

`NU_ON_NU0_RAY = 1/GAM/(1-VMU)` at ~:748 matches the derived `D` term
for term.

The P6 pure-boost experiment independently confirmed the sign:
with `chi=0` and uniform `V`, the rel-vs-norel flux difference scales
as `beta^2` (ratio 4.01 between 40000 and 20000 km/s), which is only
possible if the `+z` and `-z` segments carry `D^2` and `D'^2` with
opposite `beta*mu` signs so the O(beta) terms cancel pairwise.

## 3. Invariant transformations (Part 2)

Invariants (Lindquist 1966, sec. 5; Mihalas & Mihalas 1984, ch. 14;
Jerkstrand notes sec. on relativistic transformations):

```
I_nu / nu^3   invariant
eta_nu / nu^2 invariant
nu * chi_nu   invariant
```

With `D = nu_obs/nu_cmf`:

```
I_OBS   = D^3 * I_CMF
ETA_OBS = D^2 * ETA_CMF
CHI_OBS = CHI_CMF / D
```

Traced every application:

1. **First ray point** (~:1089-1090): `ETA_VEC(I) = ETA_VEC(I) *
   NU_ON_NU0_RAY(I)**2`, `CHI_VEC(I) = CHI_VEC(I)/NU_ON_NU0_RAY(I)` -
   applied after interpolation at the CMF frequency, evaluated at the
   local ray index.
2. **Interior points** (~:1292-1293): identical two lines inside the
   point loop.  Applied once per point, local `NU_ON_NU0_RAY(I)`.
3. **Near-side and far-side segments**: `NU_ON_NU0_RAY` is indexed by
   the same ray-point array used to fetch `ETA/CHI`, so the Doppler
   factor always matches the point's own projected velocity.  Mirror
   segments negate `VMU_RAY` before the transform loop, so far-side
   points get `D` with their own sign.
4. **STAU and ETAZ integrations**: the transfer solvers receive
   `ETA_VEC`, `CHI_VEC` already transformed; no second factor is
   applied inside the integration.
5. **Inner-boundary seed**: `PAR_FLUX = (comoving J+3*mu_cmf*H) *
   NU_ON_NU0_RAY(SM_NRAY)**3` - the D^3 intensity boost applied once,
   evaluated at the boundary point's own Doppler factor.
6. **Hollow-core back-side handoff** (~:1488-1492): `PAR_FLUX` passes
   from the far-side sweep into the near-side sweep unchanged - it is
   already observer-frame intensity, correctly NOT re-boosted.
7. **Output accumulation** (~:1720): `OBS_FLUX += HQW(LS)*PAR_FLUX`
   once per ray; `PAR_FLUX` is emergent intensity, no additional
   factor.

No factor is omitted, doubled, applied to already-transformed data,
or evaluated at the wrong ray point.  The frequency interpolation of
`ETA`/`CHI` happens at `nu_cmf = nu_obs/D` before the invariant powers
are applied - the order is required (transform of the argument first,
then the amplitude) and is what the code does.

`dFR_DATA` scaling (`6.59934*R(1)^2` Jy-at-1-kpc) was verified in P6
and is untouched.

## 4. Angular aberration at the inner boundary (Part 3)

### 4.1 Derivation

The photon direction cosine transforms under the same radial boost as
the frequency.  For a photon with observer-frame cosine `mu_obs`
relative to the radial direction, the comoving cosine is

```
mu_cmf = (mu_obs - beta) / (1 - beta*mu_obs)
```

with the inverse

```
mu_obs = (mu_cmf + beta) / (1 + beta*mu_cmf) .
```

Sign convention: same photon-propagation convention as section 2.2 -
a photon emitted outward (`mu_cmf=+1`) in a `beta>0` flow is seen at
`mu_obs=+1`; a photon emitted at `mu_cmf=0` is aberrated forward to
`mu_obs=+beta`.  Checked against the Jerkstrand one-dimensional
aberration relation and Mihalas & Mihalas ch. 14.

### 4.2 The boundary law

Production path: `IB_J_CMF_ST(ML)=RJ(ND)`, `IB_H_CMF_ST(ML)=
HFLUX_AT_IB` in `obs/cmf_flux_sub_v5.f` - both are comoving-frame
moments of the inner-boundary radiation field.  Provenance:
`RJ(ND)` is the converged comoving mean intensity of the continuum
moment solve; `HFLUX_AT_IB` is set inside `comp_j_blank.f` to
`DBB/CHI(ND)/3` under the `DIF` inner-boundary branch - the
diffusion-limit Eddington flux `(1/3*chi)*dB/dR` applied by the same
boundary condition the moment solver itself uses (non-`DIF` runs
take the incident-intensity branch instead, and `DO_IB` callers only
set the flag when `DIF` is on).  The stored pair therefore describes
precisely the two-sided diffusion field the derivation below
requires.  The emergent intensity is the Eddington relation in the
comoving frame:

```
I_cmf(mu_cmf) = J_cmf + 3*mu_cmf*H_cmf .
```

Using `mu_obs` inside this law mixes frames: the angular argument
belongs to the frame in which `J` and `H` are measured.  Fix (both
STAU ~:1585 and ETAZ ~:1657 sites):

```
IF(PLANE_PARALLEL)THEN
  T1=MU_AT_RMAX(LS)                          ! mu_obs (defect 1.3)
ELSE
  T1=Z_RAY(SM_NRAY)/R_RAY(SM_NRAY)           ! mu_obs
END IF
IF(DO_REL_CORRECTIONS)
   T1=(T1-V_RAY(SM_NRAY)/C_KMS)/
      (1.0_LDP - V_RAY(SM_NRAY)*T1/C_KMS)   ! mu_cmf
PAR_FLUX = T3 + 3*T1*T4                      ! J + 3*mu_cmf*H
IF(DO_REL_CORRECTIONS) PAR_FLUX = PAR_FLUX * NU_ON_NU0_RAY**3
```

Under `DO_REL_CORRECTIONS=F` the line is skipped - the non-relativistic
path is bitwise unchanged.

### 4.3 Moving-boundary characteristic derivation (Part A)

The phrases "outward ray" and "entering characteristic" are not
interchangeable at a boundary whose radiation data move with the gas,
so the four statements are derived separately.

**What the boundary is.**  `R(ND)` is a fixed Eulerian sphere in the
lab frame: during the spectrum computation `R` and `V` are static
input arrays and the gas velocity `V(r)` flows through every grid
surface, including `R(ND)`.  The boundary does not move; the radiation
data at it (`IB_J_CMF=RJ(ND)`, `IB_H_CMF=HFLUX_AT_IB`) are moments in
the *comoving* frame, i.e. in the frame of the gas streaming through
the fixed surface.  "The moving boundary" is shorthand for "the
boundary data live in a moving frame".

**(i) Photon travelling outward in the observer (lab) frame.**
`mu_obs = Z_RAY/R_RAY > 0` at the boundary point.  Every seeded ray
has this: the seed sites fire only for solid-wall core rays
(`SM_NRAY=NRAY`, `LS<=NC`) at the near-side boundary crossing, whose
tangent points toward the observer.  This is a statement about the
lab-frame trajectory, not a frame-dependent quantity: a boost
changes the *comoving* cosine (and can flip its sign, below), never
the fact that the photon crosses the fixed sphere outward.

**(ii) Photon travelling outward relative to the gas.**  `mu_cmf > 0`.
Since `1 - beta*mu_obs > 0` for `|beta|<1, mu_obs<=1`,

```
sign(mu_cmf) = sign(mu_obs - beta)
```

- `beta > 0` (outflow): `mu_obs > beta` gives `mu_cmf > 0`;
  `mu_obs = beta` gives `mu_cmf = 0` (comoving-tangential);
  `0 < mu_obs < beta` gives `mu_cmf < 0` - the photon's lab-frame
  outward radial speed `c*mu_obs` is slower than the gas outflow, so
  relative to the fluid it sinks back toward the core.
- `beta = 0`: `mu_cmf = mu_obs > 0` identically.
- `beta < 0` (inflow/reversing): `mu_obs > 0 > beta` gives
  `mu_cmf > 0` for every `mu_obs > 0`.  A comoving-inward seed can
  never occur on an outward lab ray under inflow.

**(iii) Characteristic entering the computational atmosphere.**  The
domain is `r >= R(ND)` in the lab frame, and the boundary is fixed in
the lab frame, so entering vs leaving is decided by `mu_obs` alone:
`mu_obs > 0` at `R(ND)` means the characteristic proceeds to larger
`r`, i.e. it enters.  All seed rays enter; none leave.  In particular
an entering characteristic can have `mu_cmf < 0` (case `0 < mu_obs <
beta`): the photon still crosses the fixed surface outward in the lab
frame; it is only *in the gas frame* that it heads inward.

**(iv) Characteristic leaving the atmosphere.**  `mu_obs < 0` at
`R(ND)` - a photon propagating into the interior.  None of the seeded
rays are of this kind (a solid wall has no incident radiation from
below; the far-side characteristic exists only for a hollow core,
which `DO_IB_DIFF_FLUX` disables).

**Which photons supply the `0 < mu_obs < beta` lab state.**  The
aberration maps the entering characteristics `mu_obs in (0,1]` onto
the comoving interval `mu_cmf in (-beta, 1]` (`mu_obs=0 -> -beta`,
`mu_obs=1 -> 1`).  For `beta>0` the lower slice `(-beta,0)` is the
interior field's *inward* hemisphere: photons that inside the gas are
heading back toward the core, but which the faster-moving fluid
advects across the fixed boundary into the atmosphere.  They are not
envelope radiation (an envelope photon at the inner boundary would be
a lab-frame *inward* characteristic, `mu_obs<0`) and they are not
emitted by the core surface - they are ordinary members of the
two-sided interior radiation field.

### 4.4 Eddington law on both hemispheres (Part B)

The seed evaluates

```
I_cmf(mu_cmf) = J_cmf + 3*mu_cmf*H_cmf
```

on `mu_cmf in (-beta,1]`, i.e. on the inward hemisphere too when
`beta>0`.  Whether that is physically justified reduces to two
questions.

**(a) Is the field at R(ND) a two-sided diffusion field?**  Yes, by
model premise.  The seed is gated on `DIF` (via
`DO_NONMONO_SOB .or. DO_NONMONO_KNOP .and. DIF` in the caller): the
diffusion inner boundary condition asserts that `R(ND)` sits deep
enough that the comoving intensity is nearly isotropic.  `IB_J_CMF`
and `IB_H_CMF` are the model's own converged comoving moments
(`RJ(ND)`, `R^2 H` from the moment solution), and the linear-in-mu
law is the unique angular distribution at first order in anisotropy
consistent with those two moments - the same closure the moment
solver itself applies at `R(ND)`.  Crucially, at an *interior*
diffusion point both hemispheres are populated: photons stream
inward and outward with a small net outward flux `H`.  Evaluating
`J + 3*mu_cmf*H` at `mu_cmf<0` is therefore not an extrapolation -
it is the interior distribution on its natural domain, and returns
the dimmer inward-hemisphere intensity `J - 3|mu_cmf|H`, which is
exactly what the advected photons of section 4.3(iii) carry.  This is
unlike an outer boundary, where the inward hemisphere is an external
incident field, not part of `J,H`.

**(b) Can the available moments determine the required incoming
boundary intensity?**  The two moments determine the linear closure
uniquely, and under the diffusion premise that is the appropriate
closure - so for `mu_cmf in (-beta,1]` the answer is determined, and
it is `J + 3*mu_cmf*H`, no special-casing and no clamp.  What the
moments cannot determine is the `O(mu^2)` structure of the true
boundary field (a `K`-moment or petal correction, e.g. limb-darkening
of the inner sphere).  A more accurate seed on *either* hemisphere
would need the `K` moment or the full angular field at `R(ND)`;
neither is stored by the moment pipeline (`RJ`, `RSQHNU` only).  The
seed is therefore as accurate as - and exactly consistent with - the
atmosphere's own inner-boundary description.  Positivity on both
hemispheres requires `3|H| < J`, guaranteed by the diffusion premise;
no positivity clamp is applied, so a model that misplaces `R(ND)`
outside the diffusion regime produces a visibly wrong seed rather
than a silently floored one.

**Conclusion.**  Retain the current prescription: the lab-frame seed
is `D^3 * (J + 3*mu_cmf*H)` at the aberrated comoving cosine over the
full domain `mu_cmf in (-beta,1]`.  It is the lab-frame transform of
the interior diffusion field evaluated at the photon's own comoving
state - the physically correct boundary value for every characteristic
entering the computational atmosphere.

### 4.5 Plane-parallel boundary cosine (defect 1.3)

In `PLANE_PARALLEL` mode `Z(I)=R(I)/MU_AT_RMAX(LS)` is a slant depth,
so `Z_RAY/R_RAY` evaluates to `1/MU_AT_RMAX(LS)` at the boundary, not
to a direction cosine.  Passing it through the aberration gives a
meaningless `mu_cmf` (and a singular denominator for
`beta*mu_obs=1`).  Fixed at both seed sites: the boundary cosine is
`MU_AT_RMAX(LS)` under `PLANE_PARALLEL`, `Z_RAY/R_RAY` otherwise.
Verified by the `char pp` battery - per-ray seeds at fixed `mu` match
the spherical seeds at the same boundary cosine and the independent
oracle to ~3e-7.

Frequency interpolation: `IB_J_CMF`, `IB_H_CMF` are interpolated at
the CMF band frequency before the `mu_cmf` combination - verified
unchanged.

The `D^3` boost and the aberration are each applied exactly once, in
the right order (evaluate the comoving law, then boost the resulting
comoving intensity to the lab frame).

## 5. Superluminal input validation (Part 4)

The old radicand floor (section 1.1) is replaced by fail-closed
validation, active only under `DO_REL_CORRECTIONS` (the established
upstream gating convention):

1. **Grid validation** (~:302, before the ray loops): every
   `V(I), I=1..ND` must be finite (`V.NE.V` catches NaN, which would
   slip through `ABS(V).GE.C_KMS`) and strictly subluminal
   `ABS(V)<C_KMS`.  Violations stop with
   `NON_FINITE_VELOCITY` or `SUPERLUMINAL_VELOCITY` naming the depth,
   value and `C_KMS`.
2. **Ray validation** (~:700, inside the LS ray loop): the
   interpolated `V_RAY(I)` must be finite; `T2=1-V_RAY^2/c^2<=0` is
   superluminal unless `|V_RAY| < c*(1+V_RTOL)`, in which case the
   value is a ulp-level interpolation overshoot of an already-validated
   grid and `T2` is reprojected to `1-(1-V_RTOL)^2`, just inside the
   domain.  Beyond `V_RTOL` the run stops with the ray index, depth
   index, `V_RAY`, `C_KMS`, geometry and `P(LS)`.
3. **Doppler denominator** (~:735): `1-VMU_RAY` must be finite and
   positive (`BAD_DOPPLER_DENOM`).  For validated subluminal input
   `|VMU|<1` mathematically, but a non-finite interpolated `VMU` would
   otherwise pass silently.

`V_RTOL = 64*EPSILON(1.0_LDP) ~ 1.4e-14`: derived from double
precision - a validated strictly-subluminal grid can only produce a
`T2<=0` interpolation through ulp-level arithmetic overshoot, bounded
by a few ulps of `V/c`; 64 ulps is generous enough to never trip on
legitimate round-off and far too small to hide physical input (a real
model cannot have `|V|` within 1.4e-14 of c by accident - and if it
does, the resulting `GAM~8e6` sends every Doppler reach past the P5
coverage guard anyway, so no bogus spectrum can emerge).

No finite Lorentz factor is manufactured from genuinely superluminal
input; the run stops.  Negative velocities are preserved throughout -
only `ABS(V)` is tested.

One ulp-level edge case: `SPEED_OF_LIGHT()*1e-5` evaluates to
299792.45800000002 in double precision, so the literal `299792.458`
sits ~2e-8 km/s *below* `C_KMS` and a caller passing that literal
value is technically subluminal.  It then trips the P5 coverage guard
downstream (verified), so it still cannot produce a bogus spectrum;
documented here for completeness.

## 6. Manufactured tests (Part 5)

New files:

- `unit_tests/f90/drv_obs_frame_rel.f` - driver with modes
  `run` (flux battery), `quad` (angular-grid dump for the independent
  boundary sum), `prop` (closed-form property checks), `sl V`
  (single-velocity validation probe), `nan` (NaN probe),
  `char sph|pp [mu]` (single-ray boundary-seed table), `pp`
  (full-grid plane-parallel battery).
- `unit_tests/python/test_obs_frame_rel.py` - harness with analytic
  and independently-evaluated predictions.

Coverage (P7 spec TEST 1-10):

1. `V=0` under `DO_REL=T` reproduces `DO_REL=F` **bitwise** -
   `GAM=1`, `NU_ON_NU0=1` exactly.
2. Closed-form Doppler values at `beta in {+/-0.9,+/-0.5,+/-0.1}`,
   `mu in {+/-1, 0}` against `sqrt((1+b)/(1-b))`, `1/gamma` -
   `prop` mode, zero failures.
3. Aberration round trip `mu -> mu_cmf -> mu` over a 5x7
   `(beta,mu)` grid including `mu_cmf<=0` cases - recovers `mu` to
   <1e-14.
4. Lorentz invariants evaluated directly on manufactured
   `I,eta,chi` with the production powers.
5. Uniform-`V` thin shell and `chi=0.5` slab against an independent
   `I(p) = integral j*D^2*e^{-tau} dz` Python evaluation (same
   geometry, same `chi*X` opacity contraction), `rtol=3e-4`.
6. `+V`/`-V` mirror pair with `chi=0` - the `D^2` multiset must
   agree to round-off; verified.
7. Non-monotonic field (acceleration, deceleration, a maximum, a
   reversal, negative cells) completes under `DO_REL` with finite
   positive flux near the non-rel value.
8. Inner-boundary aberration: `V=+/-99.9 km/s` solid wall, `ETA=0`,
   `IB_J=1`, `IB_H=0.2`; compared against an independent Python
   evaluation `sum_{LS<=NC} W*D^3*(J+3*mu_cmf*H)*e^{-tau}` at both
   signs, `rtol=2e-5`; `beta=0` recovers the non-rel seed bitwise;
   `DO_REL=F` flux is beta-independent; STAU/ETAZ agree.
9. Superluminal rejection: `V=+299793`, `+449689`, `-299793`,
   `-449689 km/s` all abort with `SUPERLUMINAL_VELOCITY`;
   `+/-0.99c` pass validation and complete; NaN aborts with
   `NON_FINITE_VELOCITY`.  No finite spectrum is produced by invalid
   input.
10. `OMP_NUM_THREADS=1,2,4,8` full `run` mode - all 18 flux files
    bitwise identical across thread counts.
11. `char` mode (Part C): 80 single-ray seeds -
    `beta in {0,+/-0.1,+/-0.3}` x `mu_obs in {1.0,0.8,0.3,0.1}` x
    STAU/ETAZ x spherical/PP - against the independent
    forward-boosted diffusion-field oracle (section 6.1), `rtol=2e-4`
    realised at ~3e-7.  The oracle builds `I'(mu')=J+3mu'H` on a
    400001-point comoving grid, maps every direction forward
    `mu=(mu'+beta)/(1+beta mu')` carrying `I' [gam(1+beta mu')]^3`,
    and interpolates at `mu_obs`; it never evaluates the inverse
    aberration used in the code.  The `beta=+0.3, mu_obs=0.1` case
    (`mu_cmf=-0.206`) is additionally asserted to be the
    inward-hemisphere value `D^3*(J-3|mu_cmf|H)` - measurably
    distinct from clamped-`J`, zero, and outward-hemisphere
    alternatives.  PP and spherical seeds agree at the same boundary
    cosine, verifying defect 1.3's fix.
12. `pp` mode (Part D): 8-ray plane-parallel grid at
    `mu in [0.05,1]`, `NC=NP` so every ray is seeded.  `V=0` rel =
    norel bitwise; `V=+/-20000` matches the oracle-weighted boundary
    sum (Doppler + aberration), `rtol=1e-4`; `+V` exceeds `-V`;
    OMP 1 vs 8 bitwise; PP direction cosine conventionally consistent
    with the spherical path (item 11 cross-check).

## 7. Non-relativistic compatibility (Part 6)

- `DO_REL_CORRECTIONS=F` skips the aberration line and the new
  validation entirely - bitwise-identical non-rel path (TEST 1
  confirms the `DO_REL=T`,`V=0` equivalence from the other side).
- P5 frequency-window construction unchanged - the validation added
  here runs before it, and manufactured grids satisfy the P5 bounds.
- P6 bitwise stock comparison unaffected: the new code only executes
  under `DO_REL_CORRECTIONS` (grid/ray validation) or inside the
  IB-diffusion branch that upstream never reaches with `DO_IB=F`
  callers.  Frozen-reference bitwise identity for the P6 battery is
  re-verified by the P6 suite itself.

## 8. Numerical results

`test_obs_frame_rel.py`: **37 passed, 0 failed** (OMP_NUM_THREADS=1
determinism bitwise against 2/4/8 threads, `max_abs_diff=0.0`).

Part C/D selected values (boundary-condition verification):

- `char` single-ray seeds: all 80 cases
  (`beta 0/+-0.1/+-0.3` x `mu_obs 1.0/0.8/0.3/0.1` x STAU/ETAZ x
  spherical/PP) match the forward-boosted diffusion-field oracle;
  worst relative deviation `3.1e-7` (tolerance `2e-4`).
- Comoving-inward case `beta=+0.3, mu_obs=0.1` (`mu_cmf=-0.20619`):
  seed `550.04` = the interior inward-hemisphere value
  `D^3*(J-3|mu_cmf|H)` (`550.04`); the clamped-`J` alternative
  predicts `627.69`, the outward-hemisphere alternative `705.35` -
  the production value is unambiguously the advected interior field.
- PP vs spherical at the same boundary cosine: all 40 pairs agree
  (`5e-5` bound), confirming `MU_AT_RMAX(LS)` is the correct PP
  boundary direction cosine.
- PP 8-ray grid: `V=0` rel = norel bitwise (`543.24934`);
  `V=+20000` `621.85307` and `V=-20000` `473.22840` vs oracle sums
  `621.85323`/`473.22852` (`2.6e-7`); `+V > -V`; OMP 1 vs 8 bitwise.

Carried-over selected values:

- `V=0` rel vs norel: `14142961.33685997` bitwise both bands.
- Thin shell `chi=1e-6`, `V=20000/40000 km/s`: `4363.8066/4363.8160`
  vs independent boosted-attenuation prediction `4364.3331/4364.3350`
  (the residual is the code's ray-grid discretisation, bounded at
  `rtol=3e-4`).
- Slab `chi=0.5`, `V=20000`: `719.3954` vs prediction `718.8124`
  (`8.1e-4`, inside the `2e-3` discretisation tolerance at tau~1).
- Mirror pair `chi=0`, `|V|` identical: `4363.83436344` vs
  `4363.83435970`, `rel=8.6e-10` (ulp accumulation across ~10^4 ray
  points; a missing mirror would be O(1)).
- Boundary seed `V=+/-99.9 km/s` vs independent
  `sum W*D^3*(J+3*mu_cmf*H)*exp(-tau)`: `19.50890079` vs
  `19.50890585` (+V) and `19.48614054` vs `19.48614560` (-V),
  `rtol=2e-5` both signs; the two seeds sit symmetric about the
  `beta=0` value `19.49751838` (`+0.011382`/`-0.011378`) - the `D^3`
  term dominates the aberration at this beta.  STAU/ETAZ agree to
  `1e-4`.  The un-aberrated prescription predicts `19.51041411` -
  measurably different, confirming the fix is doing something real
  at O(beta).
- Superluminal: `+/-299793` and `+/-449689 km/s` all abort with
  `SUPERLUMINAL_VELOCITY` naming the depth, value and `C_KMS`;
  NaN aborts with `NON_FINITE_VELOCITY`; `+/-0.99c` pass validation
  and complete with finite flux.
- Regressions: `test_obs_frame_ib` 10/10, `test_obs_frame_cov`
  15/15, `test_sob_flux_obs` 11/11, `test_obs_nonmono_dispatch`
  all pass.

## 9. Remaining limitations

- The inner-boundary prescription is limited by the moment set it is
  given: `J,H` determine only the linear-in-mu closure (section 4.4b).
  The `O(mu^2)` angular structure of a real boundary field (a
  `K`-moment or petal term) is not available from the moment pipeline
  and cannot be reconstructed; the seed is exactly as accurate as the
  atmosphere's own diffusion inner boundary, not more.  A higher-order
  prescription would require storing `K(ND)` or the full angular
  field - out of scope of this path.
- The linear law is positive on both hemispheres only under the
  diffusion premise `3|H| << J`; no positivity clamp is applied, so a
  model whose `R(ND)` is too shallow for diffusion would produce a
  visibly wrong (possibly negative) seed rather than a silently
  floored one.  Regime enforcement stays a model-level
  responsibility, consistent with `DIF`'s contract upstream.
- `V_RTOL` reprojection only arises through ulp-level overshoot of a
  validated grid; its `GAM` value is large enough that the P5
  coverage guard aborts the run downstream in practice.
- The 2e-8 km/s window between the `299792.458` literal and the fp
  `C_KMS` is documented in section 5; it is unreachable through
  VADAT-parsed model velocities at physically meaningful precision.
- The `char`/`pp` manufactured batteries use uniform `V` and flat
  `IB_J/IB_H`; depth-dependent boundary moments and non-uniform `V`
  interpolation at the boundary point are covered only by the
  regression suites, not by dedicated manufactured cases.

## 10. References

1. Lindquist, R. W. (1966), "Relativistic Transport Theory",
   Annals of Physics, 37, 487-518.  DOI: 10.1016/0003-4916(66)90207-7
   - covariant transfer equation and the `I_nu/nu^3` invariant.
2. Mihalas, D. & Mihalas, B. W. (1984), "Foundations of Radiation
   Hydrodynamics", Oxford University Press - comoving transfer,
   Lorentz transformations, aberration.
3. Jerkstrand, A., radiative-transfer lecture notes,
   https://star.pst.qub.ac.uk/webdav/public/ajerkstrand/Courses/Radiativetransfer/main.pdf
   - one-dimensional aberration relation and the intensity/
   emissivity/opacity transformation powers.
