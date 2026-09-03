# SN2024pba Workflow Timeline and File Creation

## Overview
This document describes the timeline of events and file creation in the SN2024pba model grid workflow.

## Workflow Timeline

### 1. Initial Setup (run_model_grid.py)
- Reads config file (e.g., `configs/X.yaml`)
- Generates input profile using `generate_profile_from_yaml.py`
- Creates `models/X/stella/input_profile`

### 2. STELLA Job Submission (run_model_grid.py)
- Submits STELLA SLURM job
- Job ID returned (e.g., 9651339)
- No waiting - job runs asynchronously

### 3. STELLA Job Execution (stella_job_template.slurm)

#### Step 1: Build STELLA (if needed)
- Checks for build stamp: `$STELLA_ROOT/.built_nzon{NZON}`
- If not built, runs `build_stella.sh`
- Creates build stamp to skip future builds

#### Step 2: Create model-specific working directory
- Creates `models/X/stella/`
- Sets up working directory

#### Step 3: Convert input_profile -> STELLA format
- Runs `convert_profile_to_stella.py`
- Creates `models/X/stella/modmake/{MODEL_NAME}.hyd`
- Creates `models/X/stella/modmake/{MODEL_NAME}.abn`
- Creates `models/X/stella/modmake/{MODEL_NAME}.dat`

#### Step 4: Copy/link necessary STELLA files
- Creates `models/X/stella/eve/`
- Copies `eve2.exe` to `models/X/stella/eve/`
- Creates `models/X/stella/vladsf/`
- Links `lineatom.dat` to `models/X/stella/vladsf/`

#### Step 5: Generate STELLA opacity tables
- Runs `generate_stella_opacity_tables.py`
- Creates `models/X/stella/vladsf/{MODEL_NAME}.1` through `.6`
- Creates `models/X/stella/vladsf/{MODEL_NAME}.ab`

#### Step 6: Convert .hyd + .abn to .mod format via eve2
- Runs `eve2.exe` in `models/X/stella/eve/`
- Creates `models/X/stella/eve/{MODEL_NAME}.mod`
- Creates `models/X/stella/eve/{MODEL_NAME}.rho`
- Creates `models/X/stella/eve/res.{MODEL_NAME}`

#### Step 7: Run STELLA
- Sets up `models/X/stella/strad/`
- Links `.mod`, `.xni`, and opacity files
- Runs `xstella.exe`
- Creates `models/X/stella/strad/{MODEL_NAME}.swd` (main output)
- Creates `models/X/stella/strad/{MODEL_NAME}.tt` (lightcurve)
- Copies to `models/X/stella/{MODEL_NAME}.swd`
- Copies to `models/X/stella/{MODEL_NAME}.tt`

#### Step 8: Convert .swd -> SN_HYDRO_DATA
- Runs `convert_stella_swd_to_cmfgen.py`
- **Creates `models/X/cmfgen/SN_HYDRO_DATA`** (CMFGEN input)
- Uses improved opacity calculation (line opacity + composition-dependent Z/A)
- **Grid construction pipeline (updated 2026-08-03):**
  - AMR refinement → outer boundary insertion → min-ND fill → sigma
    smoothing → grid-spacing enforcement (Phases 1-3) → velocity
    re-interpolation → **re-application of sigma smoothing** (bug fix:
    previously discarded by re-interp)
  - Phase 3 cell-ratio smoothing now fixes the **worst** violation each
    pass, not the first-found (bug fix: previously starved severe
    inner-boundary violations)
- Plots lightcurve to `models/X/plots/stella_lightcurve.png`
- Plots evolution to `models/X/plots/{MODEL_NAME}_evolution.png`
- Exports lightcurve CSV to `models/X/output/lightcurve.csv`

#### Step 9: Prepare SEDONA directory
- Creates `models/X/sedona/`
- Copies files from `sedona_base/`:
  - `param.lua` (template)
  - `plot_spectrum.py`
- Converts SN_HYDRO_DATA to SEDONA .mod format:
  - Runs `convert_cmfgen_to_sedona_v2.py`
  - **Creates `models/X/sedona/{MODEL_NAME}.mod`**
- Generates `param.lua` from template:
  - Parses SN_HYDRO_DATA header for:
    - Time since explosion (days)
    - Innermost radius (cm)
    - Innermost temperature (K)
    - Luminosity (erg/s)
  - Replaces placeholders in `param.lua`
  - **Creates `models/X/sedona/param.lua`** (final)

#### Step 10: Prepare CMFGEN directory
- Creates `models/X/cmfgen/` (already exists from Step 8)
- Copies files from `cmfgen_base/`:
  - Copies all files **except** existing ones (preserves SN_HYDRO_DATA)
  - **Preserves `models/X/cmfgen/SN_HYDRO_DATA`** (created in Step 8)
  - Adds CMFGEN configuration files

### 4. STELLA Job Completion
- Job finishes with exit code 0
- All input files ready for CMFGEN and SEDONA

### 5. CMFGEN/SEDONA Job Submission (run_model_grid.py)
- Submits CMFGEN job with dependency on STELLA job ID
- Submits SEDONA job with dependency on STELLA job ID
- Jobs wait for STELLA to complete before starting

### 6. CMFGEN Job Execution
- Uses `models/X/cmfgen/SN_HYDRO_DATA` as input
- Runs NLTE spectrum synthesis
- Creates output spectrum files

### 7. SEDONA Job Execution
- Uses `models/X/sedona/{MODEL_NAME}.mod` as input
- Uses `models/X/sedona/param.lua` as configuration
- Runs Monte Carlo radiative transfer
- Creates `spectrum.h5` and `spectrum_?.dat` files
- Runs `plot_spectrum.py` to create `spectrum_auto_plot.png`

## File Creation Order

1. `models/X/stella/input_profile` (initial setup)
2. `models/X/stella/modmake/{MODEL_NAME}.hyd` (Step 3)
3. `models/X/stella/modmake/{MODEL_NAME}.abn` (Step 3)
4. `models/X/stella/modmake/{MODEL_NAME}.dat` (Step 3)
5. `models/X/stella/vladsf/{MODEL_NAME}.1-.6` (Step 5)
6. `models/X/stella/eve/{MODEL_NAME}.mod` (Step 6)
7. `models/X/stella/strad/{MODEL_NAME}.swd` (Step 7)
8. `models/X/stella/strad/{MODEL_NAME}.tt` (Step 7)
9. **`models/X/cmfgen/SN_HYDRO_DATA`** (Step 8) - Critical file
10. `models/X/sedona/{MODEL_NAME}.mod` (Step 9) - Depends on SN_HYDRO_DATA
11. `models/X/sedona/param.lua` (Step 9) - Depends on SN_HYDRO_DATA
12. CMFGEN base files (Step 10) - Don't overwrite SN_HYDRO_DATA

## Critical Dependencies

- **SN_HYDRO_DATA** must be created before SEDONA directory preparation
- **SN_HYDRO_DATA** must be preserved when copying CMFGEN base files
- **SEDONA .mod** and **param.lua** depend on SN_HYDRO_DATA
- CMFGEN/SEDONA jobs depend on STELLA job completion

## Key Settings

### X.yaml (1-day SBO test)
- `stella.target_days`: 1.5 days
- `stella.snapshot_interval`: 0.5 days
- `cmfgen_conversion.snapshot_days`: 1.0 day

### Opacity Calculation
- Electron scattering (composition-dependent)
- Bound-free (Kramers with improved Z/A)
- Free-free (Kramers)
- Line opacity (T < 10000 K, Z-dependent)

## Verification Checklist

- [ ] STELLA job completes successfully
- [ ] SN_HYDRO_DATA created in models/X/cmfgen/
- [ ] SEDONA .mod created in models/X/sedona/
- [ ] SEDONA param.lua generated with correct values
- [ ] CMFGEN base files copied without overwriting SN_HYDRO_DATA
- [ ] Lightcurve plots generated
- [ ] Evolution plots generated
- [ ] CMFGEN/SEDONA jobs submitted with correct dependencies
