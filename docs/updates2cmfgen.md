# CMFGEN Pipeline Changes for Non-Monotonic Velocity / Interacting SNe

## Date: 2026-07-17

## Motivation

The SN2024pba interacting supernova models have a non-monotonic velocity
profile (fast inner ejecta, slow outer CSM, shock front in between).
CMFGEN's standard SN model setup assumes homologous expansion (monotonic
velocity). The changes below adapt the CMFGEN pipeline to handle
non-monotonic velocity, following the approach of
Dessart et al. (2015, MNRAS, 449, 4304):
"An essential asset of the steady-state radiative transfer calculations
with CMFGEN is the direct treatment of the nonmonotonic velocity structure."

Key issues addressed:
1. **MOM_JREL_V9 NaN at shock fronts**: The relativistic moment solver
   produces negative J where sigma < 0 (shock compression). The previous
   CHK_J=FORM_VAL tried to use JNU_STORE (formal solution from previous
   iteration), which doesn't exist during cold-start INIT, causing NaN.
2. **ND expansion from velocity jumps**: MOM_JREL_V9 inserts extra depth
   points when |V(I)-V(I+1)| > VFRAC_MOM * VDOP_VEC, causing the grid
   to grow from 150 to 218+ points and J to explode to -1e84 at inserted
   depths.
3. **BA matrix singularity at low T**: At T < 5000 K, metal excited-state
   populations are < 1e-15 of the species total, making the BA rate
   matrix singular (DGETRF IFAIL > 0 at inner depths).
4. **Outer tau = 0 in input_hydro.dat**: RV_SN_MODEL_SNIIN builds a
   geometric tau grid from TAU(outer) to TAU(inner). If the outer tau
   is 1e-30 but the second point is ~1e-4, the geometric grid puts ~80%
   of its points in the empty gap where LIN_INTERP has no hydro data,
   producing NaN R values and SIGFPE.
5. **Composition normalization**: STELLA .swd files write Y_k*A_k (number
   fraction times atomic mass), not true mass fractions. These sum to
   <A>, not 1, and must be normalized.

---

## Files Modified

### 1. `cmfgen_base/VADAT_INIT`
**Purpose:** Master CMFGEN control file (VADAT template).

**Key changes:**

| Keyword | Old | New | Reason |
|---------|-----|-----|--------|
| `TS_NO` | 1 | 0 | Static snapshot, not time-evolved SN model |
| `NMONO_J` | T (was F) | **removed** | Not needed; USE_J_REL handles non-monotonic J |
| `DO_CO_MOV_DDT` | (absent) | F | No co-moving frame time evolution (static) |
| `MAX_CF` | 15.0 | 13.5 | Avoid high-nu edge case; 30% fewer BA matrix points |
| `V_CROSS` | 1000 | 500 | 500 km/s feature resolution |
| `SIG_GAU_KMS` | 1000 | 500 | Match V_CROSS |
| `USE_JREL` | T | F | Stale keyword (not read by CMFGEN; USE_J_REL is the real one) |
| `USE_LAM_ES` | T | F | Do NOT use J=S; solve moment equations for global transport |
| `COH_ES` | T | F | Non-coherent electron scattering (incoherent Thomson for line wings) |
| `CHK_J` | FORM_VAL | ABS_VAL | **Critical fix**: patch negative J at shock fronts (\|J\|/10); fixes NaN from FORM_VAL cold-start |
| `ACC_F` | 5e-3 | 1e-2 | Relaxed Eddington convergence for DIAG stability |
| `N_IB_INS` | 10 | 15 | More boundary points for DIAG solver stability |
| `N_OB_INS` | 10 | 15 | More boundary points for DIAG solver stability |
| `VTURB` | 75 | 20 | Dessart+2015: 20 km/s turbulent velocity |
| `LAM_SET` | T | F | Keep GLOBAL_LINE=SOB active, suppress overlap warning |
| `FIX_IMP` | T | F | Impurity populations free |
| `MAX_LAM` | 100 | 1000 | MAX_LAM controls actual corrections with SCALE_OPT=MAJOR |
| `MAX_CHNG` | 1e4 | 1e100 | Abort threshold only; MAX_LAM controls actual corrections |
| `N_FIX_BA` | 5 | 50 | Many BA recomputes (need >25 for DIAG stability) |
| `BA_CHK_FAC` | 5e-4 | 1e-5 | Check BA every iteration when frozen (tighter for DIAG) |
| `FIX_BA` | 5e-2 | 5e-4 | Freeze only at < 0.05% error (tight for DIAG stability) |

**Comments added:** Detailed Dessart+2015 references for the relativistic
moment solver configuration (USE_J_REL=T, INCL_REL=T, INCL_ADV_TRANS=T,
USE_FRM_REL=T) and the CHK_J=ABS_VAL fix for negative J at shock fronts.

---

### 2. `cmfgen_base/IN_ITS`
**Purpose:** Controls the iteration loop (number of iterations, Lambda/T auto-switches).

| Keyword | Old | New | Reason |
|---------|-----|-----|--------|
| `NUM_ITS` | 3 | 25 | More init iterations (override via --init-num-its) |
| `DO_LAM_AUTO` | F | F | Pure Lambda, no auto-switch (Dessart+2015) |
| `DO_GT_AUTO` | F | T | Auto-update GT (enabled for stiff Type IIn/CSM convergence) |
| `DO_T_AUTO` | F | F | No T auto-update (FIX_T=T: T fixed to hydro) |
| `D2_EQ_D1` | F | T | Tie depth-2 populations to depth-1 for stability |

---

### 3. `cmfgen_base/MODEL_SPEC`
**Purpose:** Grid specification (depth points, rays, frequency grid).

| Keyword | Old | New | Reason |
|---------|-----|-----|--------|
| `NC` | 15 | 8 | Fewer core rays (G_ONLY/DIAG solver needs fewer) |
| `NUM_BNDS` | 3 | 1 | 1 = diagonal (DIAG) for N_TYPE=G_ONLY; 3 = tridiagonal for N_ON_J |

---

### 4. `cmfgen_base/batch.sh`
**Purpose:** C-shell batch script for CMFGEN execution.

**Change:** Updated `cmfdist` path from `${HOME}/cur_cmf` to
`/viper/ptmp1/sbrennan/SN2024pba/cmfgen_source` (project-local CMFGEN
installation).

---

### 5. `cmfgen_base/submit_cmfgen.sh`
**Purpose:** Main CMFGEN job submission script (bash).

**Key changes:**
- **20+ new CLI arguments** added for configuration via `run_model_grid.py`:
  `--target-nd`, `--model-type`, `--inc-rad-decays`, `--no-relativistic`,
  `--keep-lstar-zero`, `--init-num-its`, `--partition`, `--time`, `--mem`,
  `--n-ob-ins`, `--final-requeue-when`, `--final-requeue-max`,
  `--auto-set-rstar`, `--auto-set-lstar`, `--phot-tau`, `--bb-teff-k`,
  `--lstar-radius-mode`, `--job-prefix`, `--model-name`, `--day-label`,
  `--project-root`, `--config-file`
- **CMFDIST path** updated to project-local installation
- **Epoch label derivation**: Automatically extracts SN age from
  SN_HYDRO_DATA for centralized output filenames
- **JOB_PREFIX** from CLI arg instead of environment variable
- **PROJECT_ROOT** from CLI arg instead of deriving from directory structure

---

### 6. `scripts/prepare_VADAT.py`
**Purpose:** HPC-side script that updates VADAT and writes RDINR from
SN_HYDRO_DATA. This is the core configuration script for CMFGEN runs.

**Major additions:**

#### a. `get_rel_config()` function
New function that returns relativistic/time-dependent solver configuration
based on `model_type`:

- **`static_snapshot`** (Dessart+2015 approach for non-monotonic velocity):
  - `USE_J_REL=T` — moment equations for J (direct treatment of
    non-monotonic velocity via MOM_JREL_V9)
  - `USE_LAM_ES=F` — do NOT use J=S; solve moment equations for global
    transport
  - `CHK_J=ABS_VAL` — patch negative J at shock fronts (|J|/10)
  - `USE_DJDT_RTE=F`, `DO_DDT=F`, `DO_CO_MOV_DDT=F` — steady state
  - `INCL_REL=T`, `INCL_ADV_TRANS=T`, `USE_FRM_REL=T` — required for
    USE_J_REL=T

- **`supernova`** (default, time-evolving):
  - `USE_DJDT_RTE=T`, `DO_CO_MOV_DDT=T`, `DO_DDT=T` — time evolution
  - Relativistic solvers off

#### b. `downsample_hydro_velocity_aware()` function
New function that regrids the hydro R/V profile to `target_nd` points
with **uniform velocity spacing** (equal dV between adjacent zones).

- Places points uniformly in velocity space: `V_target = linspace(V_max, V_min, ND)`
- Interpolates R(V) from the hydro data
- Ensures |V(I)-V(I+1)| = (V_max-V_min)/(ND-1) for all adjacent zones
- Prevents MOM_JREL_V9 from inserting extra points (which happens when
  V(I)-V(I+1) > VFRAC_MOM*VDOP_VEC)
- For non-monotonic velocity, uses sorted V for interpolation (R(V) is
  not single-valued, but the uniform-V grid is still valid)
- Enforces minimum R spacing (MIN_DR_RDINR_1E10CM = 1e-10)

#### c. `compute_optimal_vfrac_mom()` function
New function that computes the optimal VFRAC_MOM to prevent CMFGEN ND
expansion.

- Analyzes the hydro velocity profile for max velocity jump and
  non-monotonicity
- Uses conservative VDOP=20 km/s as worst case
- Returns (optimal_vfrac_mom, recommend_no_relativistic)
- Bounded to [0.5, 15.0] range

#### d. `write_rdinr()` updated
- New optional `target_nd` parameter
- When set, calls `downsample_hydro_velocity_aware()` to regrid
- Returns `(nd_written, r_oi, v_oi)` for verification
- Always regrids even when target_nd == hydro ND (SN_HYDRO_DATA may have
  non-uniform velocity spacing from STELLA tau-based conversion)

#### e. `write_input_hydro_dat_sniin()` — outer tau fix
**Critical fix**: The outer tau was set to `max(tau_oi[-1], 1e-30)`,
but RV_SN_MODEL_SNIIN builds a geometric grid from TAU(outer) to
TAU(inner). If outer tau is 1e-30 but the second point is ~1e-4, the
geometric grid puts ~80% of its points in the empty gap [1e-30, 1e-4]
where LIN_INTERP has no hydro data, producing NaN R and SIGFPE.

**Fix**: Set outer tau to 1% of the second point's tau:
```python
tau_oi[-1] = max(tau_oi[-1], tau_oi[-2] * 1.0e-2)
```

#### f. Profile-specific changes (init/stable/final/single_final)

**Relativistic config**: All profiles now use `**get_rel_config(model_type)`
instead of hardcoded USE_JREL/USE_J_REL/INCL_REL/etc. values.

**Convergence parameters** (all profiles aligned):

| Parameter | Old (init) | New (init) | Old (final) | New (final) |
|-----------|-----------|-----------|------------|------------|
| `MAX_LAM` | 1e5 | 1e3 | 1e2 | 1e3 |
| `MAX_CHNG` | 1e30 | 1e100 | 1e30 | 1e100 |
| `N_FIX_BA` | 1 (init) | 50 | 5 (final) | 50 |
| `BA_CHK_FAC` | 1e-2 | 1e-5 | 5e-4 | 1e-5 |
| `FIX_BA` | 1e-2 | 5e-4 | 5e-4 | 5e-4 |
| `LAM_VAL` | 200 | 600 | 200 | 600 |
| `NUM_LAM` | 40 | 2 | 10 | 40 |

**FIX_NE strategy** (Dessart+2017):
- `init`/`stable`: FIX_NE=T (hold Ne fixed for LTE bootstrap stability)
- `single_final`/`final`: FIX_NE=F (Ne varies self-consistently)

**FIX_IMP policy (superseded 2026-07-27)**: Automatic impurity fixing is
disabled (`FIX_IMP=F`). It can replace every statistical-equilibrium equation
at a depth with fixed-population equations when the total species population
falls below its threshold, yielding a zero correction vector rather than a
recoverable Newton/Lambda correction.

**SCALE_OPT=MAJOR**: All profiles use coordinated scaling (shared scale
factor per depth) to prevent oscillatory divergence. With
SCALE_OPT=NONE, independent per-species corrections caused J to
oscillate (residuals grew 1e34 to 1e61).

**REV_RGRID**: Changed from T to F in init/stable/final profiles.
Automatic grid revision caused ADJUST_R_GRID_V4 N_ITS error when
INC_GRID=F.

**VEL_LAW=12 handling**: When VEL_LAW=12 (RV_SN_MODEL_SNIIN for
non-monotonic flows), sets:
- `METHOD=LOGMON` (log interpolation for chi.LOGLOG)
- `VFRAC_FG=VFRAC_MOM=20000` (prevents ND expansion; 20000*23.77 =
  475,400 km/s >> any velocity jump)
- `THK_CONT=T` (thick boundary for high-tau inner grid)
- `TRAP_J=T` (trapezoidal weights for J integration)

**`--no-relativistic` for static_snapshot**: Preserves USE_J_REL=T
(Dessart+2015 direct treatment of non-monotonic velocity) while
disabling time-dependent solvers (USE_DJDT_RTE=F, DO_DDT=F,
DO_CO_MOV_DDT=F).

**VFRAC_MOM optimization**: When `--target-nd` is specified with
VEL_LAW=11, forces VFRAC_MOM=1.0 (velocity-aware RDINR grid has
dV << VDOP_MIN). When VEL_LAW=12, uses the optimized VFRAC_MOM from
`compute_optimal_vfrac_mom()`.

**vadat_overrides re-application**: Config overrides are now applied
as a final step (after all profile-specific settings including
--no-relativistic) so they take precedence.

#### g. New CLI arguments
- `--target-nd N`: Target number of depth points (velocity-aware regridding)
- `--model-type TYPE`: `static_snapshot` or `supernova` (was env var
  `CMFGEN_MODEL_TYPE`)
- `--inc-rad-decays N`: `1` (default) or `0` (was env var
  `CMFGEN_INC_RAD_DECAYS`)

---

### 7. `scripts/convert_stella_swd_to_cmfgen.py`
**Purpose:** Converts STELLA .swd output to CMFGEN SN_HYDRO_DATA format.

**Key changes:**

#### a. Composition normalization
New `normalize_swd_composition()` function that normalizes STELLA .swd
composition from Y_k*A_k to true mass fractions.

The .swd file writes `Yat(ii) = YABUN(ii,Km) * AZ(ii)` (number fraction
times atomic mass). These sum to `<A> = sum(Y_k * A_k)`, not 1. The true
mass fraction is `X_k = Y_k * A_k / <A>`.

Added STELLA element names and atomic masses arrays for all 15 elements
(H, He, C, N, O, Ne, Na, Mg, Al, Si, S, Ar, Ca, Fe, Ni).

#### b. Additional hydro quantities parsed
The .swd reader now also extracts:
- `kappa` (electron scattering opacity)
- `nbar` (mean molecular weight)
- `ne` (electron number density)
- `composition` (15-element mass fractions per zone)

#### c. New CLI arguments
- `--rho-file PATH`: *.rho file from STELLA for final evolved composition
  (preferred over .abn)
- `--abn-file PATH`: *.abn file for initial composition (fallback)
- `--enable-adaptive-grid`: Dessart-style adaptive grid refinement
- `--feature-aware-downsample`: Preserve high-gradient regions (default: on)
- `--no-feature-aware-downsample`: Disable feature-aware downsampling
- `--no-smooth-sigma`: Disable Gaussian sigma-transition smoothing (default:
  disabled 2026-07-27 — it modifies STELLA velocities, creating artificial
  reversals; AMR max_sigma_jump handles smoothness via grid refinement)
- `--no-smooth-wiggles`: Disable velocity wiggle smoothing (default: disabled
  2026-07-27 — same reason as above)
- `--outer-boundary-extra-points 0`: No artificial buffer zones (default: 0
  since 2026-07-27 — CMFGEN VEL_LAW=12 inserts its own NBND_INS=3 boundary
  points; buffer zones created sigma discontinuities)

#### d. tau-min default raised
Changed from 1e-4 to 1e-3 to avoid sparse inner grid points.

---

### 8. `scripts/submit_cmfgen_jobs.sh`
**Purpose:** Wrapper script that calls submit_cmfgen.sh with parameters
from run_model_grid.py.

**Key changes:**
- **14 new positional arguments** added to match submit_cmfgen.sh:
  `KEEP_LSTAR_ZERO`, `MODEL_TYPE`, `DAY_LABEL`, `PREV_CMFGEN_DIR`,
  `TARGET_ND`, `SLURM_MEM`, `INIT_SLURM_MEM`, `N_OB_INS`,
  `FINAL_REQUEUE_WHEN`, `FINAL_REQUEUE_MAX`, `AUTO_SET_RSTAR`,
  `AUTO_SET_LSTAR`, `PHOT_TAU`, `BB_TEFF_K`, `LSTAR_RADIUS_MODE`
- **Removed environment variable exports**: No longer exports
  `CMFGEN_INIT_PARTITION`, `CMFGEN_INIT_TIME`, `CMFGEN_INIT_NUM_ITS`,
  `CMFGEN_JOB_PREFIX`, `CMFGEN_CONFIG_FILE`, `CMFGEN_INC_RAD_DECAYS`
  (all now passed as CLI arguments)
- **BASE_DIR** hardcoded to project path
- **JOB_PREFIX** derived from model name and day label

---

### 9. `scripts/cmfgen_job_template.slurm`
**Purpose:** SLURM job template for CMFGEN runs.

**Change:** Removed the call to the deprecated
`/ptmp/sbrennan/scripts/update_VADAT_from_SN_HYDRO_DATA.py` script.
VADAT updates are now handled by `prepare_VADAT.py` called from
`submit_cmfgen.sh`.

---

## Files NOT Modified (verified correct for non-monotonic velocity)

- `cmfgen_source/` — CMFGEN Fortran source code was not modified.
  The USE_J_REL=T + CHK_J=ABS_VAL combination in the VADAT provides
  the direct treatment of non-monotonic velocity without source patches.
  (The `cmfgen_non_thermal_fix.patch` file exists but is empty.)

---

## Summary of Key Fixes

1. **CHK_J=ABS_VAL**: Patches negative J at shock fronts (|J|/10) instead
   of trying to use JNU_STORE which doesn't exist during cold-start INIT.
   This was the primary cause of NaN crashes with non-monotonic velocity.

2. **Velocity-aware grid regridding**: Uniform-dV spacing in the RDINR
   grid prevents MOM_JREL_V9 from inserting extra depth points (ND
   expansion), which caused J to explode to -1e84.

3. **VFRAC_MOM optimization**: Automatic computation of VFRAC_MOM based
   on the velocity profile, or VFRAC_MOM=20000 for VEL_LAW=12 (tau-based
   grid which is not velocity-aware).

4. **Outer tau fix in input_hydro.dat**: Set outer tau to 1% of the
   second point's tau to prevent RV_SN_MODEL_SNIIN's geometric grid from
   placing 80% of points in an empty gap.

5. **Composition normalization**: STELLA .swd Y_k*A_k values normalized
   to true mass fractions (X_k = Y_k * A_k / <A>).

6. **BA matrix singularity prevention**: Population floors and `FIX_NE=T`
   for init/stable prevent singular systems; automatic impurity fixing remains
   disabled because it can zero all corrections at a depth.

7. **Convergence stabilization**: SCALE_OPT=MAJOR (coordinated scaling)
   and MAX_LAM=1e3 (with FIX_NE=T) prevent oscillatory divergence.

---

## Date: 2026-07-27 — X_2 `single_final` solver stability fix

### Motivation

The X_2 model 5.0-day CMFGEN run repeatedly crashed with NaN propagation
and a `CMF_BLKBAND_V3` singular-matrix failure around iteration 37 of the
`single_final` LTE cold-start phase. Root-cause analysis showed two
compounding solver-control problems:

1. **`--convergence-mode stable` silently re-enabled Ng acceleration for
   `single_final`.** The `single_final` profile is intended to be a pure
   damped-Lambda LTE bootstrap (`DO_NG=F`). However, the `stable`
   convergence-mode override block in `prepare_VADAT.py` did not carry a
   `single_final` entry, so it fell through to the profile default
   (`DO_NG=T`, `IBEG_NG=20`, `ITS/NG=20`). Ng acceleration then engaged at
   iteration 20 on populations that were still being clipped by the Lambda
   limiter at every step. Two rejected Ng extrapolations were followed by a
   third that was incorrectly accepted after `CMF_BLKBAND_V3` returned zero
   corrections at depth 179; that extrapolation produced the NaN state and
   the subsequent segmentation fault.

2. **`MAX_LAM=100` was too permissive for a cold start.** The Lambda
   limiter was clipping the raw residual to ±100× / −99.9% on essentially
   every iteration (the `OUTGEN` log showed `Maximum % increase ... 1.00E+04
   (LAMBDA)` and `Maximum % decrease ... 9.99E+01 (LAMBDA)` continuously
   from iteration 1 onward). With `SCALE_OPT=MAJOR` the *applied* update is
   the clipped residual times the per-depth scale factor, so a 100× cap
   still permitted order-of-magnitude population swings per iteration —
   exactly the regime in which Ng extrapolation is unsafe.

3. **`DTDR = 0/0` when `LSTAR = 0`.** With the inner luminosity set to
   zero (`LSTAR=0`, used for the cold-start interaction model), the
   `DTDR = T1/DTDR` line in `cmfgen_sub.f` evaluated `0/0` and produced
   `DTDR = NaN`, which then poisoned the `DIFFW` radiative-equilibrium
   correction vector.

### Files Modified

#### 1. `cmfgen_source/new_main/cmfgen_sub.f`
**Purpose:** Main CMFGEN iteration driver; computes `DTDR` from the
inner-boundary luminosity.

**Change:** Added a numerical guard for the `LSTAR=0` cold-start case at
the `DTDR = T1/DTDR` evaluation:

```fortran
T1=LUM*7.2685E+11_LDP/R(ND)/R(ND)
IF(LUM .EQ. 0.0_LDP .AND. DTDR .EQ. 0.0_LDP)THEN
  DTDR=0.0_LDP
ELSE
  DTDR=T1/DTDR
END IF
```

When both the inner luminosity and the accumulated `DTDR` numerator are
zero, `DTDR` is left at zero instead of evaluating `0/0`. The downstream
`IF (LAMBDA_ITERATION ...)` branch already zeros `DIFFW` during Lambda
iteration, so a finite `DTDR=0` is the correct cold-start value.

#### 2. `scripts/prepare_VADAT.py`
**Purpose:** HPC-side VADAT/MODEL_SPEC/RDINR generator; owns the
profile-specific solver controls.

**Changes:**

##### a. `single_final` profile defaults (always-on)
The `single_final` profile block now explicitly writes the Ng controls
that were previously inherited from `init`:

| Keyword | Old (inherited) | New (`single_final`) | Reason |
|---------|-----------------|----------------------|--------|
| `DO_NG` | T (inherited)   | F | LTE bootstrap must not run Ng on clipped populations |
| `IBEG_NG` | 20            | 50 | If Ng is ever re-enabled, only after 50 Lambda iterations |
| `ITS/NG` | 20             | 4 | Conservative Ng cadence for cold start |
| `MAX_LAM` | 1.0D2         | 1.0D2 | Unchanged in profile; tightened by `stable` override (see below) |

##### b. `--convergence-mode stable` override for `single_final`
The `stable` convergence-mode override block previously had no
`single_final` entry, so `single_final` runs fell through to the profile
default (`DO_NG=T`). A dedicated `single_final` override was added that
**preserves** `DO_NG=F` and tightens the Lambda limiter:

```python
"single_final": {
    "MAX_LIN": "1.1D0",     # min 1.1 to allow non-zero correction if non-Lambda mode entered
    "MAX_LAM": "1.1D0",     # 10% maximum population increase per Lambda correction during LTE bootstrap
    "LAM_VAL": "2.0D2",
    "NUM_LAM": "40",
    "N_FIX_BA": "30",       # stay in lambda until corrections are small
    "FIX_BA": "1.0D-02",   # switch to BA only when max change < 1%
    "MAX_CHNG": "1.0D100",  # Abort threshold only; MAX_LAM controls actual corrections
    "FIX_NE": "F",          # Ne free: NLTE self-consistent (global override enforces this)
    "DO_NG": "F",           # Ng must remain disabled for the single-final LTE bootstrap
    "IBEG_NG": "50",
    "ITS/NG": "4",
    "SCALE_OPT": "LOCAL",
},
```

Key points:

- **`DO_NG=F` is now explicit** in the `stable` override, so the
  convergence-mode override can no longer silently re-enable Ng for the
  LTE bootstrap. This was the direct cause of the iteration-37 NaN crash.
- **`MAX_LAM=1.1D0`** caps the *applied* per-step Lambda population
  increase at ~10%. The large `1.00E+04` and `9.99E+01` values that still
  appear in `OUTGEN` are the **raw clipped residual diagnostics** printed
  by `SOLVEBA_V13`; the actual applied update is the clipped residual
  multiplied by the per-depth scale factor, which `SOLVEBA_V13` reports as
  `The local minimum value of scale is: 1.0000E-03`. With
  `MAX_LAM=1.1D0` and `SCALE_OPT=LOCAL`, the applied multiplicative
  population update is therefore ~10% per affected depth per iteration,
  which is in the stable regime for Ng-free Lambda iteration.
- **`SCALE_OPT=LOCAL`** replaces `MAJOR` for the cold start. `MAJOR`
  shares a single scale factor per depth across all species, which is
  correct for the NLTE phase but over-damps the LTE bootstrap where
  different species need very different step sizes. `LOCAL` lets each
  species take its own (still `MAX_LAM`-capped) step.
- **`N_FIX_BA=30`, `FIX_BA=1.0D-02`** keep the solver in Lambda iteration
  until the maximum population change drops below 1%, preventing a
  premature switch to the Newton (`BA`) step while corrections are still
  O(10–100%) and therefore outside the linear regime.

##### c. `init` / `stable` / `final` `stable` overrides aligned
The `init`, `stable`, and `final` entries in the `stable` override block
were given consistent `N_FIX_BA=30`, `FIX_BA=1.0D-02`, and
`SCALE_OPT=LOCAL` settings so that all four profiles use the same
late-stage Lambda/BA switching criterion. `final` remains the only
profile that enables Ng (`DO_NG=T`, `IBEG_NG=50`, `ITS/NG=20`), and only
after the Lambda phase has driven the maximum change below
`BEG_NG=1.0D10`.

### Verification

The corrected `single_final` configuration was resubmitted for the X_2
5.0-day snapshot (SLURM job `10727915`). The run passed iteration 37 —
the exact point at which the prior (Ng-enabled) run executed the bad
extrapolation and crashed — and remained finite through iteration 37+
with:

- no Ng acceleration engaged (`DO_NG=F` held);
- no `CMF_BLKBAND_V3` singular-matrix failure;
- no `NaN` or `MON_INTERP` errors;
- finite `DTDR = 0.000000000000000E+000` (the new `LSTAR=0` guard);
- bounded applied updates (`The local minimum value of scale is:
  1.0000E-03`, i.e. ~10% per affected population per iteration).

The run is deliberately held in the damped Lambda-settling phase (it is
not yet converged) and is being monitored to confirm that the maximum
population change decreases monotonically into the `FIX_BA=1.0D-02`
switching band before any Ng acceleration is considered.

### Summary of Key Fixes

1. **Ng acceleration is now explicitly disabled (`DO_NG=F`) in the
   `single_final` `stable` convergence-mode override**, closing the
   fall-through that re-enabled Ng at iteration 20 and caused the
   iteration-37 NaN crash.

2. **`MAX_LAM` was reduced from `1.0D2` to `1.1D0` for the `single_final`
   cold start**, capping the applied per-step Lambda population update at
   ~10% instead of ~100×. This keeps the iteration in the stable regime
   for Ng-free Lambda iteration.

3. **A numerical guard was added to `cmfgen_sub.f`** so that
   `DTDR = T1/DTDR` no longer evaluates `0/0` when `LSTAR=0`, eliminating
   the `DTDR=NaN` source that poisoned the `DIFFW`
   radiative-equilibrium correction vector.

---

## Date: 2026-07-27 — Disable automatic impurity fixing after X_2 iteration-78 failure

### Failure analysis

The corrected X_2 5.0-day `single_final` run (job `10727915`) was finite
and damped through iteration 77. At iteration 76, `CMF_BLKBAND_V3` reported
a singular diagonal block at depth 178 and zeroed the correction for that
one depth. That warning alone did not cause the crash: iterations 77 and 78
still produced finite luminosities and finite `DTDR`.

The fatal transition occurred during the iteration-78 population solve:

```text
-10 levels of 12 fixed for HI at 1 depths
-10 levels of 15 fixed for HeI at 1 depths
-10 levels of  2 fixed for He2 at 1 depths
-10 levels of  7 fixed for OI at 1 depths
...
Maximum % increase at depth 0 is 0.00E+00 (LAMBDA)
Maximum % decrease at depth 0 is 0.00E+00 (LAMBDA)
```

At the next iteration, `DTDR` became `NaN`, the inner luminosity became
`NaN`, and the radius grid was subsequently corrupted (`R(1)=0`,
`R(ND)=NaN`), ending with `MON_INTERP` outside-range failure.

The `-10` level count is a **diagnostic-print defect**, not a request to
fix negative levels: `FIXPOP_IN_BA_V3` uses `FIX_N=-10` as an internal
sentinel meaning the ion is *not* fixed, but prints the final value of
`FIX_N` rather than the count from the depth at which the fix was applied.
The important signal is `at 1 depths`: automatic impurity fixing activated
at the outermost depth for nearly every active ion system (HI, HeI, HeII,
OI, SkII, CaII, FeII, and FeIII).

With `FIX_IMP=T`, `FIXPOP_IN_BA_V3` replaces an ion's statistical-equilibrium
row with an identity row and sets its residual (`STEQ`) to zero whenever the
ion's population ratio crosses the impurity threshold. In this run that
logic removed enough equations at depth 1 to leave the solver with an
all-zero correction vector (`MAXCH=0`, scale=1). CMFGEN then treated the
invalid zero result as convergence and advanced with a corrupted
outer-boundary population/radiation state. Thus this was not physical
convergence; it was a false zero-correction state caused by `FIX_IMP`.

### Changes

#### `scripts/prepare_VADAT.py`

All four VADAT profile dictionaries (`init`, `stable`, `final`, and
`single_final`) now write:

```text
FIX_IMP = F
```

This is deliberately global rather than temperature-dependent. The X_2
failure occurred in the LTE cold start with `FIX_IMP=T`; the method cannot
safely distinguish a genuinely negligible trace ion from an active
outer-boundary ion whose rate equations must remain in the coupled system.
The existing population floor and the Lambda correction limiter remain the
appropriate safeguards against zero-pivot rate matrices.

#### `cmfgen_source/new_main/cmfgen_sub.f`

Every path that enters or resumes a Lambda iteration now sets:

```fortran
FIX_IMPURITY=.FALSE.
```

This includes initial Lambda setup, hydro/grid revision, clump revision,
and reload of `IN_ITS`. It prevents a stale or manually edited `FIX_IMP=T`
VADAT from silently re-enabling automatic equation replacement during a
cold-start Lambda solve. The normal non-Lambda assignment also remains
false in the current production executable, so automatic impurity fixing
is fully disabled.

#### `configs/X2.yaml`

The `single_final` comment now records the correct policy:
`FIX_NE=F` and `FIX_IMP=F`, so every active species remains in the
statistical-equilibrium system.

### Expected behavior after the fix

- The iteration-78 outer-boundary `FIXPOP` event cannot occur.
- `SOLVEBA_V13` continues to report nonzero clipped Lambda corrections until
  they genuinely decay; it cannot mistake an all-zero `STEQ` vector induced
  by impurity fixing for convergence.
- A diagonal singularity at an individual depth is still logged and
  locally zeroed by `CMF_BLKBAND_V3`, but it will not trigger global removal
  of ion equations at the next iteration.
- `DTDR` remains finite under the existing `LSTAR=0` guard.

## 2026-07-28: CMFGEN Convergence Parameters Aligned with Manual

### Problem: Infinite Lambda Loop

The X_2 10d CMFGEN job (10747034) was stuck in an infinite Lambda loop — 549
iterations, never a single FULL linearization, with 10000% corrections that
could never converge. The model was in pure Lambda mode for ALL iterations,
with "Maximum % increase at depth 1 is 1.00E+04 (LAMBDA)" repeating every
iteration.

### Root Cause

The `stable_overrides` in `prepare_VADAT.py` (applied via `--convergence-mode
stable` to all profiles) had overly conservative parameters:

- `MAX_LAM=1.1D0` (10% cap per Lambda) — trace species needing 10000%
  correction could only change 10% per iteration → infinite loop
- `LAM_VAL=2.0D2` (200%) — triggered Lambda at 200%, always in Lambda mode
- `NUM_LAM=40` — 40 Lambda iterations, all capped at 10%, never escaped
- `DO_NG=F` — no Ng acceleration, pure Lambda converges linearly
- `N_FIX_BA=30` + `BA_CHK_FAC=1.0D-05` — BA matrix frozen, never recomputed

### Fix: Manual Recommendations

All convergence parameters updated to match CMFGEN manual (`full_descr.pdf`):

| Parameter | Old | New | Manual Quote |
|-----------|-----|-----|-------------|
| MAX_LAM | 1.1D0 | 1.0D2 | "Typically adopt 10^5" |
| LAM_VAL | 2.0D2 | 4.0D2 | "Typically > 400" |
| NUM_LAM | 40 | 2 | "Typically choose 2" |
| DO_NG | F | T | "Important tool for NUM_BNDS=1" |
| N_FIX_BA | 30 | 5 | "Recommended value is 2" |
| BA_CHK_FAC | 1.0D-05 | 5.0D-02 | "Typically 5%" |
| IBEG_NG | 50 | 15 | "At least 15 iterations" |
| ITS/NG | 4 | 8 | "We typically adopt 8" |
| BEG_NG | 5.0D0 | 1.0D1 | "5% to 10%" |
| AMASS_DOP | 4.0D0 | 1.0D+40 | "Set 1e40, use VTURB" |

### Files Modified
- `scripts/prepare_VADAT.py`: All profile dicts + stable_overrides
- `cmfgen_base/VADAT_INIT`: Template defaults
- `parameters.md`: Parameter descriptions with manual references

### Verification

Generated VADAT for X_2 10d confirmed all values match manual recommendations.
Cancelled stuck job 10747034 and resubmitted all X_2 epochs (5d-45d).


## 2026-07-28: STELLA SIGSEGV with V_A=6000 and .swd Validation

### Problem

STELLA crashes with SIGSEGV (exit 174) for V_A=6000 km/s models. The .swd file
only contained the initial condition (no post-explosion epochs), but the convert
job (afterany dependency) created invalid SN_HYDRO_DATA from this useless data.

### Fix 1: Increase AQ

`configs/X2.yaml` and `configs/X2_noscat.yaml`: `aq: 25.0` → `aq: 40.0` (higher
artificial viscosity for Mach~1600 shell-CSM collision). X2_fixed (V_A=1000,
AQ=20) survived; AQ scales with Mach number.

### Fix 2: STELLA Models Full Evolution

`homologous_evolution_days` has been removed from all configs. STELLA now models
the full evolution for all snapshot days (`target_days = max(snapshot_days) +
stella_extra_days`). No homologous evolution fallback is used.

### Fix 3: Slurm Template .swd Validation

`scripts/stella_job_template.slurm` now:
- Counts post-explosion epochs in .swd, requires ≥ 1
- Proceeds with conversion if STELLA crashed but has useful data
- Exits with error if .swd only has initial condition
- Warns (doesn't fail) if last epoch is short of target


## 2026-08-03: Grid-Conversion Bug Fixes + proxy_populations Revert (X_1)

### Bug A: `[Velocity Re-interp]` Bypasses Sigma Smoothing

`convert_stella_swd_to_cmfgen.py`'s `[Velocity Re-interp]` step (which
re-interpolates raw STELLA velocity onto the final grid to recover extrema
like the reverse-shock dip) ran *after* `smooth_sigma_transitions()` and
`enforce_grid_spacing()`, then recomputed `sigma = sigma_velocity(r, v)`
from the unsmoothed `v` — silently discarding all earlier smoothing work.

**Fix**: Re-apply `smooth_sigma_transitions()` immediately after
`[Velocity Re-interp]`, gated on the same `--smooth-sigma` flag.

### Bug B: Cell-Ratio Smoothing Fixes First-Found, Not Worst

`enforce_grid_spacing()`'s Phase 3 (cell-ratio enforcement) scanned
left-to-right and fixed the *first* violation each pass, not the worst.
Inserting a midpoint perturbs neighbouring ratios and can spawn new mild
violations near the insertion point; with first-found strategy this
created a self-sustaining churn near the outer boundary that consumed the
entire point budget without reaching a severe, isolated 22.5x cell-ratio
violation near the inner boundary (reverse-shock region).

**Fix**: Phase 3 now scans the full array each pass and fixes the
highest-ratio violation first. Verified on X_1 5.0days: worst cell ratio
22.478 → 2.000, stable across `max_points` in [400, 430, 500, 550].

### proxy_populations Reverted to False

`configs/X1.yaml` `proxy_populations: True` → `False`. Empirically (per the
`SN_12C_proxy` finding documented in AGENTS.md "proxy_populations vs CMFGEN
built-in LTE_EST"), `proxy_populations` is worse than CMFGEN's built-in `LTE_EST`
for cold snapshots because `FIX_T=T` makes CMFGEN ignore the `_IN` file's
T column and always evaluate Saha/Boltzmann at `T_local` — so any DI
computed at `T_ion` is inherently inconsistent. Confirmed on X_1 5.0days:
MAXCH stuck at ~8.4e17 after 50+ iterations with `proxy_populations=True`.

### Impact on X_1

Full 8-day chain resubmitted with both grid fixes + `proxy_populations: False`.
Verified the live `models/X_1/cmfgen/5.0days/SN_HYDRO_DATA` now has max
cell ratio 2.0 (was 22.5). See AGENTS.md "Grid-Conversion Bugs Fixed"
and "proxy_populations Reverted to False for X_1" for full details.
