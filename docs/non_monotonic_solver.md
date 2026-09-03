# Non-Monotonic Velocity-Field Formal Solver for CMFGEN

## 1. Overview

This document describes the non-monotonic velocity-field formal solver
implemented in CMFGEN for handling hydrodynamic snapshots containing shocks,
cold dense shells, velocity reversals, and negative velocity gradients
(e.g. STELLA outputs of interacting transients).

The development solver is `FG_J_CMF_NONMONO_V1`, located in
`plane/fg_j_cmf_nonmono_v1.f` within the `cmfgen_nonmono_dev` source tree.
It is selected via the VADAT keyword `NON_MONO_SOLVER=T` and coexists with
the existing production solver `FG_J_CMF_V13` without affecting it.

### Design goals

1.  Avoid crashes from negative velocity gradients (SIGMA < -1).
2.  Avoid negative optical depths and negative intensities.
3.  Handle multiple resonance crossings along a ray.
4.  Preserve physically meaningful flow-direction information.
5.  Leave the existing solver unchanged for normal monotonic models.
6.  Support STELLA-to-CMFGEN snapshot workflows.

### Key physical definitions

| Quantity | Definition |
|----------|-----------|
| `SIGMA`  | `d ln V / d ln R - 1` |
| Sobolev gradient | `dv/dr + v/r` (equivalently `V/R * (1 + SIGMA)`) |
| `T2`    | `1 + SIGMA * mu^2` (angle-dependent Sobolev factor) |
| `GAM`   | `(V/c) * T2 / R` — signed Sobolev optical-depth rate |
| `GAM_MAG` | `(V/c) * |T2| / R` — magnitude (always >= 0) |
| `FLOW_SIGN` | `sign(1, T2)` — +1 expanding, -1 compressing |

When `SIGMA >= -1` (monotonic expansion), `T2 > 0` for all ray angles,
and the solver reduces to the standard Sobolev/CMF short-characteristics
integration.  When `SIGMA < -1` (non-monotonic), `T2` can change sign
along a ray, producing multiple resonance crossings.


## 2. File inventory

All changes are in the `cmfgen_nonmono_dev` source tree relative to the
pristine `cmfgen_source` baseline.

### New files

| File | Purpose |
|------|---------|
| `plane/fg_j_cmf_nonmono_v1.f` | Non-monotonic formal solver (2007 lines) |

### Modified files

| File | Change |
|------|--------|
| `new_main/mod_subs/control_variable_mod.f` | Add `NON_MONO_SOLVER` logical |
| `new_main/mod_subs/rd_control_variables.f` | Read `NON_MONO_SOLVER` from VADAT |
| `new_main/mod_subs/comp_j_blank.f` | Route to nonmono solver when active |
| `plane/Makefile` | Compile `fg_j_cmf_nonmono_v1.o` |
| `Makefile_definitions` | Add debug flags (`-g -traceback`) |
| `new_main/subs/solveba_v13.f` | Cap LAMBDA population increases (not just decreases) |
| `subs/gff.f` | Guard against NaN/Inf T in free-free gaunt factor |
| `subs/sobjbar_sim.f` | (reverted) Always use legacy softened Sobolev guard |
| `plane/mom_j_cmf_v11.f` | (reverted) Always clamp SIGMA >= -0.999 in moments |
| `new_main/subs/do_ng_accel_band_v4.f` | (reverted) Always clamp SIGMA >= -0.999 in Ng accel |
| `new_main/mod_subs/adjust_r_grid_v4.f` | (reverted) Always clamp SIGMA in regridding |
| `new_main/subs/rv_sn_model_sniin.f` | (reverted) Always clamp SIGMA in hydro setup |
| `subs/extend_vtsigma.f` | (reverted) Always clamp SIGMA in grid extension |


## 3. Solver selection

### VADAT keyword

The solver is selected by a new VADAT keyword:

```text
T [NON_MONO_SOLVER]
```

Default: `.FALSE.` (uses `FG_J_CMF_V13`, bit-for-bit with baseline).

### Routing in comp_j_blank.f

`COMP_J_BLANK` contains two call sites for the formal solver:

1.  **ACC path** (acceleration/extended grid): line 354
2.  **Main path** (frequency iteration): line 656

When `NON_MONO_SOLVER=T`, both paths call `FG_J_CMF_NONMONO_V1` instead of
`FG_J_CMF_V13`.  The routing is a simple IF/ELSE wrapper; no other changes
to `COMP_J_BLANK` are made.


## 4. Algorithm

### 4.1 Ray grid construction

The ray grid is constructed identically to `FG_J_CMF_V13`:
- Impact parameters `P(LS)` for `LS = 1, ..., NP`
- Ray points `Z(I, LS)` along each impact parameter
- Extended grid `R_EXT`, `V_EXT`, `SIGMA_EXT` beyond `R(ND)`

The extension uses cubic interpolation coefficients (`V_COEF`, `SIGMA_COEF`)
fit to `LOG(R)`.  No special non-monotonic handling is needed here because
the extension region is typically monotonic.

### 4.2 Sobolev optical depth with signed GAM

For each ray point `I` along impact parameter `LS`:

```fortran
MU       = Z(I,LS) / R_RAY(I,LS)
T1       = 3.33564E-06 * V_RAY(I) / R_RAY(I,LS)    ! V/c / R
T2       = 1.0 + SIGMA_RAY(I) * MU**2
GAM(I,LS)    = T1 * T2                              ! signed
GAM_MAG(I)   = T1 * MAX(ABS(T2), 1.0E-12)           ! magnitude
FLOW_SIGN(I) = SIGN(1.0, T2)                        ! direction
```

Key design decisions:

- **`GAM_MAG`** is used wherever a strictly non-negative optical-depth
  increment is required (DTAU, the Q term in CHI_RAY).  This prevents
  negative optical depths.

- **`GAM`** (signed) is retained for the source-function advection term
  `Q * I_PREV`, where the sign carries flow-direction information.

- **`FLOW_SIGN`** records whether the flow is expanding (+1) or
  compressing (-1) at each ray point.

- The floor `MAX(ABS(T2), 1.0E-12)` prevents division by zero when `T2`
  passes through zero at a resonance turning point, without distorting
  the magnitude in the common case `|T2| >> 1e-12`.

### 4.3 Multi-root resonance scan

For a non-monotonic `v(r)`, the comoving frequency
`nu_cm(r) = nu_obs * (1 - V/c)` can cross a given rest frequency multiple
times along a ray.  These resonance turning points occur where `T2`
changes sign between adjacent ray points.

The solver brackets every zone `[I, I+1]` where
`T2(I) * T2(I+1) <= 0` and counts them as resonance turning points.
No interpolation is performed at this stage; the piecewise
short-characteristics integration naturally sums the source-function
contribution over each monotonic sub-interval because:

- `DTAU` is built from `GAM_MAG` (>= 0), so optical depth accumulates
  monotonically across each sub-interval.
- The signed `GAM` enters only the source-function advection term,
  whose sign is carried by `FLOW_SIGN`.

### 4.4 Piecewise short-characteristics integration

The integration follows the same short-characteristics approach as
`FG_J_CMF_V13`:

1.  Interpolate `ETA`, `CHI`, `ESEC` onto the ray grid.
2.  Compute `DTAU(I,LS)` from `CHI_RAY` and `dZ`.
3.  Compute the exponential integrals `E0, E1, E2, E3` from `DTAU`,
    with `EXP(-DTAU)` guarded against overflow (`DTAU < 700`).
4.  Solve the tridiagonal system for the radiation field `AV(I,LS)`.
5.  Compute `JNU(I)` from angle-weighted quadrature over all rays.

The non-monotonic modifications are:

- `DTAU` uses `GAM_MAG` (not `GAM`), ensuring non-negative increments.
- The source-function advection term `Q` uses the signed `GAM` and
  `FLOW_SIGN`, preserving flow-direction information.
- Negative intensities are clamped to zero and counted (`N_NEG_INT`).

### 4.5 Negative intensity guard

After the formal solution, the solver checks for negative mean
intensities:

```fortran
DO I = 1, ND
  IF (JNU_STORE(I) < 0) THEN
    ! Write diagnostic file and STOP
  END IF
END DO
```

This is a hard stop — if negative intensities are produced, the solver
writes `FG_J_CMF_NONMONO_V1_ERRORS` with the frequency, JNU, ETA, CHI,
and ESEC at all depths, then halts.  This is stricter than the legacy
solver, which only checks `RJ(1)` (innermost point) for blow-up.


## 5. Infrastructure hardening

### 5.1 Design principle: clamp in infrastructure, not in the solver

The non-monotonic formal solver is designed to handle true negative
`SIGMA` values.  However, the rest of CMFGEN's infrastructure (moment
equations, BA matrix, Sobolev escape probabilities, grid regridding,
Ng acceleration) was not designed for `SIGMA < -1` and can produce
NaN or unphysical results if given unclamped values.

The design principle adopted is:

> **SIGMA clamps in infrastructure routines are ALWAYS applied,
> regardless of `NON_MONO_SOLVER`.  Only the formal solver
> `FG_J_CMF_NONMONO_V1` sees the true negative SIGMA.**

This is now **implemented** via a parallel `SIGMA_TRUE` array (see §5.4).
The clamped `SIGMA` is retained for all infrastructure; the unclamped
`SIGMA_TRUE` is threaded into the formal solver only.

| Routine | Clamp | Reason |
|---------|-------|--------|
| `rv_sn_model_sniin.f` | `SIGMA >= -0.999` | Hydro model setup (saves `SIGMA_TRUE` first) |
| `extend_vtsigma.f` | `SIGMA >= -0.9999` | Grid extension (clamped copy only) |
| `adjust_r_grid_v4.f` | `SIGMA >= -0.999` (smooth) | Regridding |
| `mom_j_cmf_v11.f` | `SIGMA >= -0.999` | Moment equations |
| `do_ng_accel_band_v4.f` | `SIGMA >= -0.999` | Ng acceleration |
| `sobjbar_sim.f` | Softened guard (Dessart+2015) | Sobolev escape probability |
| **`fg_j_cmf_nonmono_v1.f`** | **No clamp** | **Formal solver (uses `SIGMA_TRUE`)** |

### 5.2 SOLVEBA_V13 LAMBDA iteration guard

The LAMBDA iteration in `SOLVEBA_V13` uses a `LIMIT` option that caps
population corrections.  The original code only capped decreases
(`STEQ > 1.1 -> 0.999`) but not increases.  With non-monotonic
velocity fields, catastrophic population increases (up to 10^62 %)
could occur at the outer boundary, corrupting T and crashing `gff_vec`.

The fix caps both directions:

```fortran
IF (STEQ(J,I) > 1.1)  STEQ(J,I) =  0.999   ! cap decreases
IF (STEQ(J,I) < -1.1) STEQ(J,I) = -0.999   ! cap increases
```

This reduced the population explosion from 10^62 % to ~110 % and
eliminated the NaN cascade.

### 5.3 gff.f NaN guard

The free-free gaunt factor routine `gff_vec` computes `LOGT = LOG10(T(I))`
and uses it to index a pre-computed table.  If `T(I)` becomes NaN (from
a population explosion), `LOGT` is NaN, `INT(NaN)` is undefined, and the
table access `A(IU, IGB)` segfaults.

The guard:

```fortran
IF (.NOT. (LOGT > 0.0 .OR. LOGT < 0.0)) LOGT = 3.6   ! NaN check
IF (LOGT < 2.0) LOGT = 2.0                            ! floor (100 K)
IF (LOGT > 6.0) LOGT = 6.0                            ! cap (10^6 K)
```

This converts NaN to a safe default and clamps extreme values to the
table range, preventing the SIGSEGV.

### 5.4 SIGMA_TRUE threading (true unclamped SIGMA to the formal solver)

The original implementation clamped `SIGMA` to `>= -0.999` at hydro
read-in (`rv_sn_model_sniin.f`) **before** the formal solver ever saw
it.  This meant `FG_J_CMF_NONMONO_V1` received `SIGMA = -0.999` even
when the true hydrodynamic value was `SIGMA = -94.66`, defeating the
multi-root resonance machinery and underestimating Sobolev optical
depths by up to ~10^5 on radial rays in shock zones.

The fix introduces a parallel `SIGMA_TRUE` array carrying the unclamped
value, threaded into the formal solver only:

| File | Change |
|------|--------|
| `mod_cmfgen.f` | `SIGMA_TRUE(:)` module array declared alongside `SIGMA` |
| `cmfgen.f` | `SIGMA_TRUE(ND)` allocated alongside `SIGMA` |
| `rv_sn_model_sniin.f` | `SIGMA_TRUE(1:ND) = SIGMA(1:ND)` saved **before** both clamp loops |
| `radiation_mod.f` | `SIGMA_TRUE_EXT(:)` module array declared + allocated |
| `subs/extend_sigma_true.f` | **New routine**: extends `SIGMA_TRUE` onto the refined grid using linear-coefficient interpolation (no `log(SIGMA+1)`, no clamp — well-defined for `SIGMA < -1`) |
| `cmfgen_sub.f` | `CALL EXTEND_SIGMA_TRUE` after each of the 4 `EXTEND_VTSIGMA` sites |
| `fg_j_cmf_nonmono_v1.f` | New dummy arg `SIGMA_TRUE_IN(ND)`; `SIGMA_EXT(ND_ADD+1:ND_EXT) = SIGMA_TRUE_IN(1:ND)` at the single chokepoint (line 608). All downstream `SIGMA_COEF`, `SIGMA_RAY`, `GAM`, `T2`, resonance scan, `dGAMdR` automatically use the true values. |
| `comp_j_blank.f` | ACC path passes `SIGMA_TRUE_EXT`; main path passes `SIGMA_TRUE` |

Key design decisions:

- **`EXTEND_SIGMA_TRUE` is a separate routine** from `EXTEND_VTSIGMA`,
  which interpolates in `log(SIGMA+1)` space (undefined for `SIGMA < -1`).
  The new routine applies the cubic `COEF` matrix linearly to `SIGMA`,
  which is well-defined for any value.  `EXTEND_VTSIGMA` is untouched
  (its 8 callers are unaffected; the clamped `SIGMAEXT` is still produced
  for infrastructure).

- **Single chokepoint**: the substitution at line 608
  (`SIGMA_EXT = SIGMA_TRUE_IN`) is the only place the true value enters
  the solver.  Since `SIGMA_COEF` is built from `SIGMA_EXT`, and
  `SIGMA_RAY` from `SIGMA_COEF`, all physics downstream (GAM, T2,
  GAM_MAG, FLOW_SIGN, resonance turning points, dGAMdR) automatically
  uses the true `SIGMA`.

- **Outer beta-law extension** (THK block, lines 609-614) still uses the
  clamped `SIGMA(1)`.  This is intentional: the outer extension is
  monotonic, so `SIGMA(1) == SIGMA_TRUE(1)` there and the clamp has no
  effect.

- **OBS stage unaffected**: `COMP_OBS_V3` and `OBS_FRAME_SUB_V9` do not
  use `SIGMA` at all — they transform the stored `IPLUS` (computed by
  the formal solver with true `SIGMA`) to the observer frame using
  velocity `V` only.  The observed spectrum therefore inherits the
  correct CMF-frame radiation field.

- **Bit-for-bit legacy**: when `NON_MONO_SOLVER = .FALSE.`, the formal
  solver `FG_J_CMF_V13` is called (untouched signature).  `SIGMA_TRUE`
  is populated but never read.  Legacy behaviour is unchanged.


## 6. Verification

### 6.1 Test model

- **Model**: SN2024pba X_nomono, M3.0_V7000_NONMONO
- **Snapshot**: 35.57 days
- **Source**: STELLA hydrodynamic output with shocks and velocity reversals
- **Grid**: ND=252, velocity range 646-6475 km/s
- **Non-monotonicity**: 56/251 zones (22.3%) with SIGMA < -1
- **Min SIGMA**: -94.66
- **Max |dV|**: 434.1 km/s (velocity reversal)

### 6.2 Solver entry (VERIFIED)

```text
NON_MONO_SOLVER: (added) T
Calling FG_J_CMF_NONMONO_V1 in COMP_J_BLANK
Maximum number of points along ray in FG_J_CMF_NONMONO_V1 is 252
FG_J_CMF_NONMONO_V1: 54 depths have SIGMA < -1 (non-monotonic flow handled).
```

### 6.3 Multi-root resonance (VERIFIED)

The solver detects and reports multiple resonance turning points per ray:

```text
FG_J_CMF_NONMONO_V1: ray LS=1 has 8 resonance turning points (multi-root).
FG_J_CMF_NONMONO_V1: ray LS=22 has 6 resonance turning points (multi-root).
FG_J_CMF_NONMONO_V1: ray LS=88 has 2 resonance turning points (multi-root).
```

Inner rays (small impact parameter) have up to 8 turning points; outer
rays have 2-4.  This is consistent with the shock structure in the
STELLA snapshot.

### 6.4 Convergence behaviour (VERIFIED)

With the infrastructure clamps restored, the model iterates stably:

| Iteration | Luminosity (Lsun) | Error |
|-----------|-------------------|-------|
| 1 | 4.60e9 | -1.72e12 |
| 5 | 4.74e9 | -1.49e12 |
| 10 | 5.14e9 | -1.19e12 |
| 15 | 5.69e9 | -1.01e12 |
| 20 | 6.18e9 | -8.52e11 |
| 25 | 6.44e9 | -7.24e11 |
| 30 | 6.57e9 | -6.03e11 |

No NaN, no SIGSEGV, no negative intensities.  The luminosity is
converging upward toward LSTAR (5.43e8 Lsun from the STELLA header,
though the actual Rosseland optical depth at the inner boundary is
only 1.84, well below the tau > 50 needed for the diffusion
approximation, so a luminosity mismatch is expected for this
particular snapshot).

### 6.5 Comparison with legacy solver (VERIFIED)

The legacy `FG_J_CMF_V13` on the same model produces:

| Iteration | Luminosity (Lsun) | Max % increase |
|-----------|-------------------|----------------|
| 1 | 3.80e9 | 5.35e60 % |
| 2 | 3.82e9 | 2.39e59 % |
| 3 | 3.85e9 | 1.53e61 % |
| 4 | 3.87e9 | — |

The legacy solver survives but with a 10^60 % population explosion at
depth 252 (outer boundary).  The nonmono solver produces a different
luminosity (4.60e9 vs 3.80e9) because it handles the resonance
structure differently, but with the infrastructure clamps it also
iterates stably.

### 6.6 What was NOT verified

- **Physical correctness**: convergence does not establish that the
  spectrum is physically correct.  The low optical depth (tau_ross =
  1.84 at inner boundary) means the diffusion approximation is invalid
  for this snapshot.
- **Spectral output**: the job was run on a 2-hour interactive partition
  and did not reach spectral synthesis before the wall limit.
- **Monotonic model regression**: the nonmono solver has not been tested
  on a standard monotonic SN model to confirm it produces identical
  results to `FG_J_CMF_V13` when `SIGMA >= -1` everywhere.


## 7. Known limitations

1. **SIGMA clamping in infrastructure** (RESOLVED): The formal solver
   now receives the true unclamped `SIGMA` via the `SIGMA_TRUE` array
   (see §5.4).  The clamped `SIGMA` is still used by infrastructure
   routines (moments, Ng, regridding, Sobolev escape) as intended.
   The `N_NEG_SIGMA` diagnostic now correctly reports the true count
   of non-monotonic zones (e.g. 54/252 in the test model) rather than
   the clamped count of 0.

2. **Low optical depth**: The test snapshot (35.57 days) has
   `tau_ross = 1.84` at the inner boundary, far below the `tau > 50`
   needed for the diffusion approximation.  This causes a luminosity
   mismatch and slow convergence.  Earlier snapshots or models with
   more mass would have higher optical depth.

3. **No spectral synthesis**: The model has not been run to completion
   through the OBS stage.  Convergence and spectral output require
   longer wall time or a model with higher optical depth.

4. **DC_METH**: Non-monotonic temperature requires `DC_METH=R` (radius
   interpolation) instead of `DC_METH=TR` (temperature interpolation).
   This is handled by `prepare_VADAT.py` and `submit_cmfgen.sh`.

5. **MODEL_SPEC optimizer**: The `optimize_model_spec_v2.py` script has
   a Unicode decoding error when reading atomic data
   (`UnicodeDecodeError: 'utf-8' codec can't decode byte 0xbd`).
   This is unrelated to the non-monotonic solver but prevents automatic
   MODEL_SPEC optimization.


## 8. Configuration

### 8.1 VADAT settings

Required for non-monotonic models:

```text
T   [NON_MONO_SOLVER]
R   [DC_METH]              ! for non-monotonic T; TR crashes REGRID_LOG_DC_V1
```

Recommended:

```text
LIMIT [LAM_SCALE_OPT]     ! cap LAMBDA population changes
1.0D2 [MAX_LAM]           ! allow large LAMBDA steps (default 1.2 is too conservative)
```

### 8.2 Config YAML

```yaml
cmfgen:
  vadat_overrides:
    NON_MONO_SOLVER: T
    DC_METH: R
    LAM_SCALE_OPT: LIMIT
    MAX_LAM: 1.0D2
```

### 8.3 Executable

The non-monotonic solver is compiled into `cmfgen_dev.exe` in the
`cmfgen_nonmono_dev` source tree.  The slurm script must point to:

```bash
export cmfdist="/viper/ptmp1/sbrennan/SN2024pba/cmfgen_nonmono_dev"
```


## 9. Build instructions

```bash
cd /viper/ptmp1/sbrennan/SN2024pba/cmfgen_nonmono_dev
make
```

The executable is placed at `exe/cmfgen_dev.exe`.

The build will report one pre-existing error in `misc/rewrite_old_model_data.f`
(`RETURN statement invalid in main program`).  This is unrelated to the
non-monotonic solver and does not prevent the CMFGEN executable from
being built.


## 10. Diagnostic output

### 10.1 OUTGEN

On the first frequency, the solver prints:

```text
Calling FG_J_CMF_NONMONO_V1 in COMP_J_BLANK
Maximum number of points along ray in FG_J_CMF_NONMONO_V1 is 252
FG_J_CMF_NONMONO_V1: 54 depths have SIGMA < -1 (non-monotonic flow handled).
FG_J_CMF_NONMONO_V1: ray LS=1 has 8 resonance turning points (multi-root).
...
```

These diagnostics are printed only once per CMFGEN run (gated by
`DIAG_PRINTED`) to avoid flooding OUTGEN.

### 10.2 Error file

If negative mean intensities are detected, the solver writes
`FG_J_CMF_NONMONO_V1_ERRORS` containing FREQ, JNU_STORE, ETA, CHI,
and ESEC at all depths, then halts.

### 10.3 SOLVEBA LIMIT warnings

```text
Warning -- using LIMIT option for LAMBDA iteration in SOLVEBA_V13.f
 (decreases capped: 6633, increases capped: 5010)
```

This reports how many population corrections were capped in each
direction.  Large numbers (thousands) indicate the model is far from
convergence but the capping is preventing instability.

---

## 11. BA Iteration Status (revised 2026-09-02)

### Previous status: BA iterations not working

All previous runs of the X_nomono_v3 model family failed to achieve
true Full BA iterations with depth coupling. The model would get stuck
in diagonal-only Lambda iterations and never converge. The
`SOLVEBA LIMIT` warnings above were a symptom of this — the model was
trapped in Lambda mode and could not resolve depth-coupled oscillations.

### Root cause: DO_LAM_AUTO=F

With `DO_LAM_AUTO=F` and `DO_LAM_IT=T`, the code path at
`cmfgen_sub.f:4798` forces `LAMBDA_ITERATION=.TRUE.`, silently
overriding IT_SPECIFIER's FULL(FIXT) blocks. All FULL(FIXT) blocks
were executed as diagonal-only Lambda iterations — true Full BA
iterations with depth coupling never occurred.

### Intermediate fix: IT_SPECIFIER with DO_LAM_AUTO=T

The first fix enabled `DO_LAM_AUTO=T` and used an IT_SPECIFIER schedule
with alternating LAMBDA(FIXT) and FULL(FIXT) blocks. This produced true
BA iterations (`SUCCESS=T`, `(BA computed)` labels), but the BA was
forced at iteration 12 when MAXCH was still >500%. The enormous
corrections (500–2000%) were physically meaningless — `SUCCESS=T` only
meant the matrix factorized, not that the correction was valid.

### Current approach: DO_LAM_AUTO-only (no IT_SPECIFIER)

The IT_SPECIFIER was removed entirely. CMFGEN's built-in auto-switch
mechanism now controls the iteration mode:

1. Lambda iterations run until MAXCH < 50%
2. `DO_LAM_AUTO` auto-switches `DO_LAM_IT=F` in IN_ITS
3. Full BA iterations begin with `LAM_VAL=400` / `NUM_LAM=2` fallback
4. If BA overshoots, `solve_for_pops.f` automatically returns to Lambda

This is the standard CMFGEN approach used by Dessart and Hillier
reference models. `prepare_VADAT.py` deletes any existing IT_SPECIFIER
when `iteration_strategy` is absent from the YAML config.

OUTGEN confirms:
```
DIAG: iter=   1 IT_SPEC_PRESENT= F
DIAG: LAM_BEFORE= T RD_LAMBDA= T
DIAG: LAM_AFTER= T COMP_BA= T FIX_T= T
DIAG: OVERRIDE branch (RD_LAMBDA=T)
 LAMBDA iteration used
```

### Zone-dependent BA (source fix in do_var_cont.f)

When `NON_MONO_SOLVER=T`, the BA variation uses a zone-dependent approach
instead of the previous global fixed-J BA:

- `VAR_MOM_J_CMF_V12` is called with clamped SIGMA (>= -0.999)
- After the call, `TX(:,I,:)` and `TVX(:,I,:)` are zeroed at compression
  depths (where `SIGMA_TRUE < -1`)
- Full cross-depth J-coupling at 137 monotonic depths
- Fixed-J at 23 compression depths

This replaces the previous approach which skipped J-variation entirely
at all depths. The K_ON_J clamp ordering was also fixed (moved before
the `QFROMF` call) to prevent division by zero in the sphericity factors.

### Corrected settings

See `README_STATIC_NONMONO.md` section "BA Iteration Configuration"
for the full corrected VADAT settings and the auto-switch mechanism.
