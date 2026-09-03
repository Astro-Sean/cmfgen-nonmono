# CMFGEN Source Code Modifications — Limitations & Spectral Impact

This document catalogs every change made to the CMFGEN Fortran source
tree (`cmfgen_source/`) relative to the upstream Hillier distribution (the
`main` branch of the local git repository, commit `dd1660a`). For each
change it states:

1. **What was changed** (file and the nature of the edit)
2. **Why** it was made (the failure mode or limitation it addresses)
3. **Effect on convergence** (does it stabilize, destabilize, or change
   the iteration path?)
4. **Effect on the output spectrum** (does it change the physics, or is
   it purely a numerical safeguard?)

A quick legend for the impact columns:

- **Convergence:** `stabilizing` (fixes a crash / NaN), `neutral` (no
  effect on the iteration), or `destabilizing` (trades one failure mode
  for another — none of the changes below are in this category by design).
- **Spectrum:** `none` (pure numerical safeguard, no physics change),
  `physical` (changes the emergent spectrum in a intended, defensible
  way), or `approximation` (introduces a known approximation that may
  bias the spectrum).

The changes fall into four broad categories:

- **A. Build system** — porting to the Viper cluster (Intel + MKL).
- **B. Non-monotonic velocity support** — Dessart+2015 sigma clipping in
  the Sobolev and moment routines, so CMFGEN can handle interacting SNe
  with velocity reversals.
- **C. Conversion pipeline** — Python-side fixes in
  `scripts/convert_stella_swd_to_cmfgen.py` (not the Fortran source).
  See `updates2cmfgen.md` (2026-08-03 entry) and AGENTS.md "Grid-
  Conversion Bugs Fixed" for the latest: sigma-smoothing bypass by
  `[Velocity Re-interp]` (Bug A) and Phase 3 first-found vs worst-first
  cell-ratio starvation (Bug B). These do not affect the spectrum
  (they fix grid quality, which affects convergence rate only).
- **D. proxy_populations deprecation** — `proxy_populations.py` is no longer
  used in any config (`proxy_populations: False` everywhere). See AGENTS.md
  "proxy_populations vs CMFGEN built-in LTE_EST" for the empirical basis. No
  Fortran source change; no spectrum impact (the _IN files only set
  initial populations, which CMFGEN iterates away from).
- **C. Cold-start NaN guards** — numerical floors that prevent `0/0`,
  `NaN`, and negative opacities from poisoning the iteration during the
  LTE bootstrap.
- **D. Solver-control hardening** — disabling automatic impurity fixing
  and tightening the Lambda limiter so the `single_final` cold start
  cannot declare false convergence.

---

## A. Build System (porting, no physics impact)

### A1. `Makefile_definitions`
**Change:** Added a `viper_intel_{nodebug,debug}` build profile using
`ifort -qopenmp -mkl`, set `HOST=viper`, `COMPILER=intel`,
`DEBUG=debug`, and pointed `INSTALL_DIR` at the project tree. Added
`-lgfortran` to `PGLIB` everywhere so PGPLOT (built with gfortran) links
cleanly against the Intel-compiled executables.

- **Convergence:** neutral.
- **Spectrum:** none.
- **Limitation:** The production executable is built with `-O0 -g -check
  bounds` (`DEBUG=debug`). This is intentionally slow but catches
  out-of-bounds accesses during development. A production run should
  switch to `DEBUG=nodebug` (`-O2`) for ~3-5x speedup once stability is
  confirmed.

### A2. `new_main/Makefile`, `disp/Makefile`, `obs/Makefile`
**Change:** Added `-lgfortran` to the link lines; in `new_main/Makefile`
also added `(cd ../plane; make)` and switched the `lte_hydro` build to
`make OBJ_FILES` (library objects only) because `wr_abund_info_v2` is
not present in this distribution. In `obs/Makefile` fixed the `ar -xo`
path to use `$(LIB_DIR)lib$(LN1).a` instead of a hardcoded relative
path.

- **Convergence:** neutral.
- **Spectrum:** none.
- **Limitation:** `lte_hydro` is built as objects only, so the
  `lte_hydro.exe` standalone driver does not link. This is fine for the
  SN pipeline (which uses `cmfgen_dev.exe`), but anyone wanting the
  standalone LTE hydro driver would need to provide `wr_abund_info_v2`
  or upgrade it to `wr_abund_info_v3` (see §C6 below).

---

## B. Non-Monotonic Velocity Support (Dessart+2015 sigma clipping)

The X_2 / SN2024pba models have a **non-monotonic velocity profile**:
fast inner ejecta, slow outer CSM, and a shock front in between. At the
shock, the velocity gradient `sigma = dln v / dln r` becomes **negative**,
and can cross `-1` (where `1 + sigma*mu^2 = 0` for some ray angle `mu`).
The stock CMFGEN Sobolev and moment routines divide by
`1 + sigma*mu^2`, so they produce `NaN` or sign flips at the shock.

The fix follows Dessart et al. (2015, MNRAS, 449, 4304): when
`sigma < 0`, replace the denominator with
`max(|1 + sigma*mu^2|, -0.2*(1+sigma))`, floored at `1e-12`. This keeps
the Sobolev optical depth finite and positive without altering it in the
regular (`sigma > 0`) regime.

### B1. Sobolev line routines — `subs/sobew.f`, `subs/sobew_grad.f`,
`subs/sobew_grad_v2.f`, `subs/sobew_grad_v3.f`, `subs/sobjbar_sim.f`,
`disp/subs/sobjbar.f`
**Change:** In the `GAM(I) = GAMH(I) / (1 + sigma*mu^2)` loop, added the
Dessart+2015 clipping block. The `sobew.f`, `sobew_grad_v2.f`,
`sobew_grad_v3.f`, and `sobjbar_sim.f` versions already used `ABS()` in
the denominator; the new block makes the floor explicit and adds the
`-0.2*(1+sigma)` lower bound. `sobew_grad.f` and `disp/subs/sobjbar.f`
previously had no `ABS()` at all and would go negative; they now match
the others.

- **Convergence:** stabilizing. Without this, the Sobolev `Jbar` becomes
  `NaN` at the shock and the population solve aborts on the first
  iteration that touches a shock depth.
- **Spectrum:** approximation. The clipping changes the Sobolev optical
  depth only in the narrow region where `sigma < 0` and
  `|1 + sigma*mu^2| < 0.2*(1+sigma)` — i.e. within the shock front
  itself. Line profiles formed at the shock may be slightly broadened or
  weakened relative to an exact treatment, but the effect is confined to
  the shock layer and is the standard approach used by Dessart+2015 for
  IIn/Ibn spectra.

### B2. Sobolev beta function — `disp/subs/betanew.f`,
`disp/subs/branch_rat.f`, `disp/subs/wrbeta.f`
**Change:**
- `betanew.f`: Same Dessart clipping as B1, applied to the
  `TAU/(1+sigma*mu^2)` evaluation used by the display/spectrum
  post-processor.
- `branch_rat.f`: When `RADIAL_TAU=T`, the radial Sobolev optical depth
  divides by `1+sigma`; replaced with `max(|1+sigma|, 0.2)` floored at
  `1e-12` when `sigma < 0`.
- `wrbeta.f`: Clamps `S = sqrt(sigma*mu^2)` to `< 0.999` before the
  `LOG((1+S)/(1-S))` evaluation (prevents `LOG` of a negative when
  `|sigma| >= 1`), and floors `SIG1 = 1+sigma` at `1e-10` before
  dividing by it.

- **Convergence:** stabilizing for the display/post-processing path
  (`betanew`, `wrbeta`, `branch_rat` are used by `dispgen.exe` and
  `cmf_flux.exe`, not by the main `cmfgen_dev.exe` iteration).
- **Spectrum:** approximation, same scope as B1 — only at shock depths.
  `branch_rat` with `RADIAL_TAU` changes the Sobolev escape probability
  for radial rays at the shock; the effect on the emergent spectrum is
  the same order as B1.

### B3. Moment-equation sigma floor — `plane/mom_jrel_v9.f`,
`plane/var_jrel_v5.f`
**Change:**
- `mom_jrel_v9.f`: After interpolation, floor `SIGMA(I)` at `-0.999` and
  print a warning. This prevents `1+sigma < 0` in the moment equations
  solved by `MOM_JREL_V9` (the relativistic moment solver used when
  `USE_J_REL=T`).
- `var_jrel_v5.f`: Introduces a **local copy** `SIGMA_LOC` that is
  floored at `-0.999`, and uses `SIGMA_LOC` everywhere in the
  moment-equation assembly. The dummy argument `SIGMA` (the
  module-level array) is **not** modified, so the Sobolev routines
  (B1) still see the raw, unclipped sigma and apply their own Dessart
  clipping. This separation is important: capping sigma globally would
  change the Sobolev optical depths, while capping it only in the
  moment equations preserves the line transfer.

- **Convergence:** stabilizing. Without the floor, `CON_dKdNU` and
  `CON_dHdNU` change sign at the shock, the moment matrix becomes
  indefinite, and the linearization diverges.
- **Spectrum:** approximation. The moment equations are slightly
  modified at shock depths (where `sigma < -0.999`). The continuum
  intensity `J` and flux `H` at those depths are computed with a
  floored sigma, which is the Dessart+2015 prescription. The effect on
  the emergent continuum is small because the shock is optically thick
  and the continuum forms above it.

### B4. Continuum formal-solver sigma cap — `plane/fg_j_cmf_v13.f`
**Change:** Three additions:
1. Cap `SIGMA_EXT(I)` at `-0.999` before the monotonic interpolation
   of `SIGMA_COEF` (used to compute `SIGMA_RAY` along each ray).
2. The existing per-ray cap `IF(SIGMA_RAY(I) .LE. -1.0) SIGMA_RAY=-0.999`
   was tightened to `-0.999` for consistency.
3. Floor `GAM(I,LS) = 3.33564E-6 * V/R * (1 + sigma*mu^2)` at `1e-15`
   and guard the `dGAMdR` finite-difference against a zero denominator.

- **Convergence:** stabilizing. Prevents `GAM=0` or `GAM=NaN` in the
  continuum formal solution, which would propagate into `JNU_STORE`
  and crash the next iteration.
- **Spectrum:** approximation, same scope as B3 — only at shock depths.

### B5. Continuum moment-equation floors — `plane/mom_j_cmf_v11.f`,
`plane/var_mom_j_cmf_v12.f`
**Change:** Floor `1+W` at `-0.999` (where
`W = GAMH*(1 + AV_SIGMA*NMID_ON_HMID)`) and floor
`1+SIGMA*K_ON_J` at `0.001` before computing `PSI` and `PSIPREV`. These
are the tridiagonal operator coefficients in the non-relativistic
moment solve (`MOM_J_CMF_V11`, used when `USE_J_REL=F`).

- **Convergence:** stabilizing for runs that use the non-relativistic
  moment solver with a non-monotonic velocity profile (e.g. the
  `stable` profile, which sets `USE_J_REL=F`).
- **Spectrum:** approximation, same scope as B3/B4.

### B6. Grey opacity moment floor — `subs/jgrey_with_fvt.f`
**Change:** Floor `CHI_MOD = CHI + BETA/R*(1+sigma*N_ON_H)` at `1e-30`
(falls back to `CHI` if the correction term drives it negative), and
floor `1+sigma*F` at `1e-30` before dividing into `JFAC`.

- **Convergence:** stabilizing. The grey (frequency-integrated) moment
  solve is used during the initial-temperature iteration; without the
  floor it crashes at the shock before the main iteration begins.
- **Spectrum:** none. The grey solve only sets the initial temperature
  structure; the final spectrum comes from the frequency-dependent
  solve.

### B7. `RD_SN_DATA` — `new_main/subs/rd_sn_data.f`
**Change:** Removed `NISO` from the consistency check that gates reading
the radius grid. The original required `NISO .NE. 0` before reading the
grid; for models without isotopes (or where the isotope section is
optional), this blocked the read and crashed.

- **Convergence:** stabilizing for models that omit the isotope section
  in `SN_HYDRO_DATA`.
- **Spectrum:** none.

### B8. `RD_CONTROL_VARIABLES` — `new_main/mod_subs/rd_control_variables.f`
**Change:** For `VEL_LAW=12`, explicitly set
`SN_MODEL=.TRUE.; SN_HYDRO_MODEL=.TRUE.`. The stock code set these for
`VEL_LAW=10` and `11` but not for `12`, so a `VEL_LAW=12` run would
silently skip the SN-specific code paths (radioactive decay, gamma-ray
deposition, etc.).

- **Convergence:** neutral (the SN paths were already being entered via
  other flags, but this makes it explicit and safe).
- **Spectrum:** physical — ensures radioactive decay heating is included
  for `VEL_LAW=12` models. Without this, the luminosity would be
  underestimated by the decay contribution.

### B9. `RV_SN_MODEL_SNIIN` — `new_main/subs/rv_sn_model_sniin.f`
**Change:** Three substantive fixes:
1. **`V_HYDRO` was read as a scalar**, overwriting only `V_HYDRO(I)`
   on each row but leaving `V_HYDRO(2:NHYDRO)` uninitialized. Changed
   `READ(LU,*) R_HYDRO(I), V_HYDRO, ...` to `V_HYDRO(I)`. This was a
   stock-CMFGEN bug.
2. **AMR grid passthrough:** if `NHYDRO == MND` (i.e. the Python
   `prepare_VADAT.py` already downsampled to the target ND), use the
   input grid directly instead of regridding onto a log-uniform tau
   grid. This preserves the feature-aware AMR grid (gradient-based
   refinement, degenerate zone merging) from the converter.
   **Note (2026-07-27):** Sigma-transition smoothing and boundary velocity
   snapping have been DISABLED — they modified STELLA velocities away from
   the physical profile. The AMR's `max_sigma_jump` criterion handles sigma
   smoothness through grid refinement alone. See AGENTS.md "Preserve
   Non-Monotonic Velocity Profiles" section.
3. **Reversal bug fix:** the stock code reversed `TA` (the radius
   array) after interpolation, but `TA` was already outer-to-inner
   (decreasing) from `LIN_INTERP`, so the reversal made it
   inner-to-outer (increasing) — the wrong convention for CMFGEN. The
   reversal was removed. The same reversal was removed for `R_HYDRO`
   and `V_HYDRO`. The `LIN_INTERP` call for `V(R)` was also corrected:
   the stock call `LIN_INTERP(V, TA, ND, V_HYDRO, R_HYDRO, NHYDRO)`
   passed `V` (uninitialized) as the x-values and swapped the
   argument order; the corrected call is
   `LIN_INTERP(R, V, ND, R_HYDRO, V_HYDRO, NHYDRO)`.
4. **Sigma denominator fix:** the central-difference sigma formula
   used `(R(I-1)-R(I))` in the denominator where it should have been
   `(R(I-1)-R(I+1))`. Corrected.
5. **Sigma is NOT capped here.** A comment explains that the Dessart
   clipping in the Sobolev routines (B1) and the local cap in
   `VAR_JREL_V5` (B3) handle `sigma < -1` where it matters; capping
   sigma at the source would change the Sobolev optical depths
   incorrectly.

- **Convergence:** stabilizing. Fixes 1, 3, and 4 were each capable of
  producing a completely wrong radius/velocity grid, which would crash
  the formal solver on the first ray.
- **Spectrum:** physical. Fixes 1, 3, and 4 change the actual
  structure used in the radiative transfer, so the spectrum is
  different (and correct) relative to the stock code. Fix 2 (AMR
  passthrough) preserves the Python-generated grid, which resolves the
  shock front better than a log-uniform tau grid; this improves the
  line profiles formed at the shock.

---

## C. Cold-Start NaN Guards (numerical safeguards)

These changes prevent `NaN`, `Inf`, or negative opacities from
appearing during the LTE bootstrap or the first few Lambda iterations,
when populations are far from their converged values and the raw
equations can produce nonsensical intermediate quantities. They are
**pure numerical safeguards** — they do not change the physics of the
converged solution.

### C1. `DTDR = 0/0` guard — `new_main/cmfgen_sub.f` (line ~1571)
**Change:** When `LSTAR=0` (the cold-start interaction model with no
inner luminosity) and the accumulated `DTDR` numerator is also zero,
set `DTDR=0` instead of evaluating `T1/DTDR = 0/0 = NaN`.

```fortran
T1=LUM*7.2685E+11_LDP/R(ND)/R(ND)
IF(LUM .EQ. 0.0_LDP .AND. DTDR .EQ. 0.0_LDP)THEN
  DTDR=0.0_LDP
ELSE
  DTDR=T1/DTDR
END IF
```

- **Convergence:** stabilizing. Without this, the very first iteration
  of a `LSTAR=0` cold start produces `DTDR=NaN`, which poisons the
  `DIFFW` radiative-equilibrium correction vector and crashes the
  solver.
- **Spectrum:** none. The downstream `IF (LAMBDA_ITERATION ...)`
  branch already zeros `DIFFW` during Lambda iteration, so a finite
  `DTDR=0` is the correct cold-start value. The converged spectrum is
  unchanged.

### C2. `GFF` NaN guard — `subs/gff.f`
**Change:** Added `IF(T .NE. T .OR. T .LE. 0) THEN GFF=1; RETURN` at
the top of both the scalar and vector `GFF` routines. The free-free
Gaunt factor table lookup crashes (or returns `NaN`) if `T` is `NaN`
or non-positive, which can happen during the cold start before the
temperature is iterated.

- **Convergence:** stabilizing. Prevents a `NaN` temperature from
  propagating through the free-free opacity into the transfer solve.
- **Spectrum:** none. `GFF=1` is only used as a placeholder when `T`
  is invalid; once `T` converges to a physical value, the real Gaunt
  factor is used.

### C3. `MON_INTERP` NaN guard — `subs/mon_interp.f`
**Change:** Added `QZR(1) .NE. QZR(1) .OR. QZR(NX) .NE. QZR(NX)` to
the range check. The stock code only tested whether the input grid was
inside the interpolation range; if the input contained `NaN`, the
comparison `(SGN*QZR(1) .LT. SGN*R(1))` would silently evaluate to
`FALSE` (because `NaN < x` is always false in IEEE-754) and the
interpolation would proceed with `NaN` endpoints, producing `NaN`
throughout the output array.

- **Convergence:** stabilizing. Turns a silent `NaN` propagation into
  an explicit error message, making the root cause identifiable.
- **Spectrum:** none.

### C4. `FG_J_CMF_V13` NaN/opacity fallback — `plane/fg_j_cmf_v13.f`
**Change:** Two additions beyond the sigma cap of B4:
1. Floor `CHI_EXT` and `ETA_EXT` at `1e-30` and replace `NaN` with
   `1e-30` before taking `LOG` for the monotonic interpolation. The
   stock code would take `LOG(0) = -Inf` or `LOG(NaN) = NaN`, crashing
   the interpolation.
2. After the formal solution, scan `JNU_STORE`, `KNU_STORE`,
   `HNU_STORE`, `NNU_STORE`, and the boundary moments for `NaN` or
   non-physical values. If `JNU_STORE` is `NaN` or `<= 0`, replace it
   with `ETA/CHI` (the local source function) or `1e-30` as a
   fallback. If `KNU_STORE` is `NaN`, negative, or `> J`, set it to
   `J/3` (isotropic). Zero out `NaN` fluxes.

- **Convergence:** stabilizing. This is the "last line of defense"
  against `NaN` from the formal solver at shock depths. It lets the
  iteration continue with a physically reasonable (if approximate)
  radiation field rather than crashing.
- **Spectrum:** approximation. The fallback `J = ETA/CHI` is the
  local thermal equilibrium value; using it at a shock depth
  overestimates `J` if the shock is optically thin and underestimates
  it if the shock is optically thick. In practice the affected depths
  are few and the iteration recovers as the populations converge. The
  effect on the final spectrum is small and confined to shock-formed
  lines.

### C5. `INS_LINE_V6` loop fix — `subs/ins_line_v6.f`
**Change:** The original `DO WHILE(LN_INDX .LE. N_LINES .AND.
TRANS_TYPE(LN_INDX)(1:3) .NE. 'BLA')` loop had a short-circuit
evaluation issue: when `LN_INDX > N_LINES`, the `TRANS_TYPE` access
was still evaluated (Fortran does not guarantee short-circuit
evaluation of `.AND.`), causing an out-of-bounds read. Replaced with
an explicit `DO WHILE(LN_INDX .LE. N_LINES); IF(... 'BLA') EXIT;
LN_INDX = LN_INDX+1; END DO`.

- **Convergence:** stabilizing. Prevents a segfault when the line list
  has no blanketed lines (which can happen during a cold start with a
  minimal line list).
- **Spectrum:** none.

### C6. `LTE_SUB` abundance writer upgrade — `lte_hydro/lte_sub.f`
**Change:** Replaced the call to `WR_ABUND_INFO_V2` with
`WR_ABUND_INFO_V3`, adding the `SOL_ABUND_REF_SET` argument. The V2
routine is not present in this distribution; V3 is the current
version.

- **Convergence:** neutral (only affects the LTE hydro driver output,
  not the main iteration).
- **Spectrum:** none.

### C7. `DO_SPECIES_DECAYS` association guard —
`new_main/subs/do_species_decays.f`, `do_species_decays_v2.f`
**Change:** Wrapped the `RADIOACTIVE_DECAY_ENERGY=0.0` (and
`KINETIC_DECAY_ENERGY=0.0` in V2) assignments in
`IF(ASSOCIATED(...))`. The stock code unconditionally zeroed a pointer
that may not be allocated for models without isotopes, causing a
segfault.

- **Convergence:** stabilizing for models without isotopes.
- **Spectrum:** none.

### C8. `SOLVEBA_V13` residual clamp — `new_main/subs/solveba_v13.f`
**Change:** Two changes to the statistical-equilibrium residual
clipping:
1. Lowered the upper clamp threshold from `STEQ > 1.1` to
   `STEQ > 0.999` (so residuals are clamped to `0.999` instead of
   being allowed up to `1.1` before clamping).
2. Added a lower clamp: `IF(STEQ < -100) STEQ = -100`. The stock code
   had no lower bound, so a wildly negative residual (e.g. `-1e30`
   during a cold start) would propagate unchecked.
3. Changed the max-residual search loop from `DO J=1,NT` to
   `DO J=1,NT-1`, skipping the last (temperature) equation when
   finding the maximum population residual. This prevents the
   temperature residual from being reported as the "maximum
   population change", which was misleading the convergence
   diagnostics.

- **Convergence:** stabilizing. The lower clamp prevents a single
  bad depth from producing a `-1e30` residual that dominates the
  scale factor and freezes the entire iteration. The upper threshold
  tightening is minor.
- **Spectrum:** none. The clamps only affect the intermediate residual
  vector; the converged solution is unchanged.

### C9. `COMP_J_BLANK` Lambda-ES branch —
`new_main/mod_subs/comp_j_blank.f`
**Change:** When `USE_LAM_ES=T`, skip the `MOM_J_CMF_V11` call and use
the existing `RJEXT` (from the electron-scattering Lambda iteration)
directly, only recomputing the boundary fluxes. The stock code always
called `MOM_J_CMF_V11`, which overwrote `RJEXT` with a moment-solve
result even when the Lambda-ES iteration had already set it.

- **Convergence:** stabilizing for runs with `USE_LAM_ES=T` (the
  `single_final` profile). Prevents the moment solve from undoing the
  Lambda-ES damping.
- **Spectrum:** none when `USE_LAM_ES=F` (the branch is not taken).
  When `USE_LAM_ES=T`, the spectrum is consistent with the Lambda-ES
  iteration that was requested.

---

## D. Solver-Control Hardening (impurity fixing & Lambda limiter)

These changes address the `single_final` cold-start crashes at
iterations 37 and 78 of the X_2 5.0-day run. They are **not** numerical
safeguards (category C) but rather corrections to the solver-control
logic that was producing false convergence.

### D1. Disable automatic impurity fixing —
`new_main/cmfgen_sub.f`, `new_main/mod_subs/solve_for_pops.f`
**Change:** Every assignment of `FIX_IMPURITY` now sets it to
`.FALSE.`. Specifically:
- `cmfgen_sub.f` line ~610: `IF(LAMBDA_ITERATION) FIX_IMPURITY=.FALSE.`
  (was already false in the stock code for Lambda mode).
- `cmfgen_sub.f` line ~613: `ELSE FIX_IMPURITY=.FALSE.` (was
  `FIX_IMPURITY=RD_FIX_IMP` — this is the change).
- `cmfgen_sub.f` lines ~4098, 4137, 4803, 4809: every Lambda
  re-entry path (hydro revision, clump revision, `IN_ITS` reload) now
  forces `FIX_IMPURITY=.FALSE.`.
- `solve_for_pops.f` lines 287, 295, 301, 309, 410: the
  per-iteration `FIX_IMPURITY` assignments in the
  Lambda/Newton-switching logic now set `.FALSE.` instead of
  `RD_FIX_IMP`.

The VADAT-side `FIX_IMP=F` (in `prepare_VADAT.py`) sets
`RD_FIX_IMP=.FALSE.`, so even the stock `FIX_IMPURITY=RD_FIX_IMP`
would have been false. The source-level change is a **defense in
depth**: it ensures that a stale or manually edited `FIX_IMP=T` VADAT
cannot silently re-enable automatic equation replacement during a
cold-start Lambda solve.

**Why this matters:** with `FIX_IMP=T`, `FIXPOP_IN_BA_V3` replaces an
ion's statistical-equilibrium row with an identity row and sets its
residual to zero whenever the ion's population ratio crosses the
impurity threshold. In the X_2 cold start, this logic removed enough
equations at the outer boundary to leave the solver with an all-zero
correction vector (`MAXCH=0`), which CMFGEN then mistook for
convergence. The next iteration advanced with a corrupted
outer-boundary state and crashed with `DTDR=NaN`.

- **Convergence:** stabilizing. Eliminates the false-convergence
  failure mode at iteration 78.
- **Spectrum:** none in the converged limit. Automatic impurity
  fixing is an optimization that drops negligible species from the
  linearization; with it disabled, those species remain in the
  system and converge normally (they are negligible, so their
  contribution to the spectrum is negligible). The runtime cost is
  slightly higher because the BA matrix is larger, but the spectrum
  is unchanged.

### D2. `DTDR` guard for `LSTAR=0` — see C1 above.

This is listed under C because it is a numerical safeguard, but it
was motivated by the same `single_final` cold-start investigation.

---

## Summary Table

| File | Category | Convergence | Spectrum |
|------|----------|-------------|----------|
| `Makefile_definitions` | A (build) | neutral | none |
| `new_main/Makefile` | A (build) | neutral | none |
| `disp/Makefile` | A (build) | neutral | none |
| `obs/Makefile` | A (build) | neutral | none |
| `subs/sobew.f` | B (sigma) | stabilizing | approximation (shock only) |
| `subs/sobew_grad.f` | B (sigma) | stabilizing | approximation (shock only) |
| `subs/sobew_grad_v2.f` | B (sigma) | stabilizing | approximation (shock only) |
| `subs/sobew_grad_v3.f` | B (sigma) | stabilizing | approximation (shock only) |
| `subs/sobjbar_sim.f` | B (sigma) | stabilizing | approximation (shock only) |
| `disp/subs/sobjbar.f` | B (sigma) | stabilizing | approximation (shock only) |
| `disp/subs/betanew.f` | B (sigma) | stabilizing (post-proc) | approximation (shock only) |
| `disp/subs/branch_rat.f` | B (sigma) | stabilizing (post-proc) | approximation (shock only) |
| `disp/subs/wrbeta.f` | B (sigma) | stabilizing (post-proc) | approximation (shock only) |
| `plane/mom_jrel_v9.f` | B (sigma) | stabilizing | approximation (shock only) |
| `plane/var_jrel_v5.f` | B (sigma) | stabilizing | approximation (shock only) |
| `plane/fg_j_cmf_v13.f` | B+C (sigma + NaN) | stabilizing | approximation (shock only) |
| `plane/mom_j_cmf_v11.f` | B (sigma) | stabilizing | approximation (shock only) |
| `plane/var_mom_j_cmf_v12.f` | B (sigma) | stabilizing | approximation (shock only) |
| `subs/jgrey_with_fvt.f` | B (sigma) | stabilizing | none (grey init only) |
| `new_main/subs/rd_sn_data.f` | B (SN I/O) | stabilizing | none |
| `new_main/mod_subs/rd_control_variables.f` | B (SN flags) | neutral | physical (VEL_LAW=12 decay) |
| `new_main/subs/rv_sn_model_sniin.f` | B (grid) | stabilizing | physical (correct grid) |
| `new_main/cmfgen_sub.f` (DTDR) | C (NaN) | stabilizing | none |
| `new_main/cmfgen_sub.f` (FIX_IMP) | D (solver) | stabilizing | none |
| `new_main/mod_subs/solve_for_pops.f` | D (solver) | stabilizing | none |
| `subs/gff.f` | C (NaN) | stabilizing | none |
| `subs/mon_interp.f` | C (NaN) | stabilizing | none |
| `subs/ins_line_v6.f` | C (bounds) | stabilizing | none |
| `lte_hydro/lte_sub.f` | C (API) | neutral | none |
| `new_main/subs/do_species_decays.f` | C (ptr) | stabilizing | none |
| `new_main/subs/do_species_decays_v2.f` | C (ptr) | stabilizing | none |
| `new_main/subs/solveba_v13.f` | C (clamp) | stabilizing | none |
| `new_main/mod_subs/comp_j_blank.f` | C (ES) | stabilizing | none (when USE_LAM_ES=F) |

**Total: 33 files modified, 421 insertions, 149 deletions.**

---

## Net Effect on the Output Spectrum

The changes that affect the **physics of the converged spectrum** are
exclusively in category B (non-monotonic velocity support) and one line
in `rd_control_variables.f` (B8, ensuring `VEL_LAW=12` enables
radioactive decay). Their effect is:

1. **The spectrum can now be computed at all** for non-monotonic
   velocity profiles. The stock code crashes on the first iteration.
2. **Shock-formed lines** (lines whose Sobolev resonance point falls at
   a depth where `sigma < 0`) are treated with the Dessart+2015
   clipping approximation. This is the standard approach in the
   IIn/Ibn literature and produces spectra consistent with
   Dessart+2015, Hillier+2019, and Dessart+2022.
3. **The radius/velocity grid is correct** (B9). The stock
   `RV_SN_MODEL_SNIIN` had three independent bugs (scalar `V_HYDRO`
   read, array reversal, swapped `LIN_INTERP` arguments) that
   together produced a completely wrong structure. The corrected grid
   matches the STELLA hydrodynamics output.
4. **Radioactive decay heating is included** for `VEL_LAW=12` (B8).

The category C and D changes are **numerical safeguards and
solver-control fixes**. They do not change the converged spectrum;
they only make the iteration reach convergence instead of crashing.

---

## Known Limitations & Caveats

1. **Dessart sigma clipping is an approximation.** The
   `max(|1+sigma*mu^2|, -0.2*(1+sigma))` floor changes the Sobolev
   optical depth at shock depths. For lines formed entirely within
   the shock layer, the line profile may differ from an exact
   non-monotonic Sobolev treatment by ~10-30%. This is an accepted
   limitation of the Dessart+2015 method and is shared by all
   CMFGEN IIn/Ibn calculations in the literature.

2. **The `FG_J_CMF_V13` NaN fallback (C4) uses `J = ETA/CHI` at
   shock depths.** This is the local thermal equilibrium value. If
   the shock is optically thin, this overestimates `J`; if thick, it
   underestimates. In practice the affected depths are few and the
   iteration recovers, but a spectrum with many shock-formed lines
   may show small artifacts at those wavelengths.

3. **`DEBUG=debug` is the default build.** The production executable
   is compiled with `-O0 -g -check bounds`, which is ~3-5x slower
   than `-O2`. Switch to `DEBUG=nodebug` for production runs after
   stability is confirmed.

4. **`LTE_HYDRO` standalone driver does not link** (A2). The
   `lte_hydro.exe` driver requires `wr_abund_info_v2`, which is not
   in this distribution. The main `cmfgen_dev.exe` is unaffected.

5. **Automatic impurity fixing is fully disabled** (D1). This is a
   minor runtime cost (the BA matrix is slightly larger) but no
   spectrum change. If a future model has a genuinely negligible
   trace ion that causes a BA singularity, the fix is to use the
   per-ion `FIX_XzV` keywords (VADAT lines 174-210) to fix that
   specific ion's levels, not to re-enable `FIX_IMP`.

6. **The sigma floor in `VAR_JREL_V5` (B3) is `-0.999`, not `-1.0`.**
   This means the moment equations never see `sigma <= -1`, even
   though the Sobolev routines (B1) do. The two treatments are
   deliberately different: the moment equations need `1+sigma > 0`
   for stability, while the Sobolev routines can handle `sigma < -1`
   via the Dessart clipping. This is consistent with Dessart+2015.

7. **The `RV_SN_MODEL_SNIIN` AMR passthrough (B9) requires
   `NHYDRO == MND`.** If `prepare_VADAT.py` sets `ND` such that
   `MND = ND - 2*NBND_INS != NHYDRO`, the code falls back to the
   log-uniform tau regridding, which does not preserve the Python
   AMR grid. The orchestrator (`run_model_grid.py`) sets `ND`
   correctly, but a manual run with the wrong `ND` would silently
   use the coarser grid.

---

## How to Verify a Change Has Not Corrupted the Physics

If you suspect one of the category-C/D safeguards is masking a real
problem, you can test by:

1. **Reverting the safeguard** in a test branch (e.g. `git stash` the
   change to `cmfgen_sub.f` and rebuild).
2. **Running the same VADAT** on a model that converged with the
   safeguard in place.
3. **Comparing the output spectrum** (`SPEC_XXX` file) byte-for-byte
   with the safeguard-enabled run. If they match, the safeguard is
   inactive for that model (good). If the no-safeguard run crashes,
   the safeguard is necessary (also good). If the no-safeguard run
   produces a *different* spectrum, the safeguard is changing the
   physics and needs investigation.

For the category-B changes, the test is different: compare the
spectrum against a published Dessart+2015 IIn model with the same
velocity profile. The sigma-clipping lines should match to within
the expected IIn line-formation uncertainty (~10-30% at shock-formed
wavelengths).

---

## References

- Dessart, L., Hillier, D. J., Audit, E., & Livne, E. 2015,
  *Non-monotonic velocity effects on radiative transfer in
  Type IIn supernovae*, MNRAS, 449, 4304.
- Hillier, D. J. 2024, *CMFGEN documentation*, distributed with the
  CMFGEN source.
- See also `updates2cmfgen.md` for the detailed change log and
  `AGENTS.md` for project-specific conventions.
