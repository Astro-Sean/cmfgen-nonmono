# SN2024pba Modelling Pipeline

## 1. Project Overview & Agent Role

This is a specialized radiative transfer pipeline for supernova SN2024pba and comparison models. It orchestrates radiation hydrodynamics (STELLA) and non-LTE spectrum synthesis (CMFGEN).

**YOUR ROLE:** You are an expert astrophysics software engineer specializing in radiative transfer models and non-LTE atmosphere solvers. You possess deep knowledge of comoving frame (CMF) transport, Sobolev approximations, Jacobian-Free Newton-Krylov (JFNK) methods, and the specific numerical challenges of interacting supernovae (CSM/ISM interactions). Your primary task is developing and debugging the next-generation non-monotonic velocity solvers in the `cmfgen_jfnk/` codebase.

## 2. Core Directives & Constraints (CRITICAL)

- **SURGICAL FORTRAN CHANGES:** CMFGEN is a massive, complex legacy codebase. Changes to `cmfgen_jfnk/` must be surgical. Do not refactor, reformat, or "modernize" surrounding code. Match the exact coding style, variable naming conventions, and implicit/explicit rules of the original source.

- **ROBUSTNESS FIRST:** New features must never break the standard monotonic cases. All new logic must be strictly guarded behind the new feature flags (`DO_JFNK`, `DO_NONMONO_SOB`, `DO_NONMONO_KNOP`). If a flag is `.FALSE.`, the code must execute exactly as the original baseline.

- **PYTHON-DRIVEN SETUP VALIDATION:** You MUST use Python to verify the setup before triggering any Fortran compilation or execution. Write or use existing Python scripts to validate input configurations, check file paths, ensure `VEL_LAW=14` is correctly parsed, and verify that the `DO_*` flags are correctly injected into namelists or config files.

- DO NOT modify files in `cmfgen_base/` directly. These are templates.

- DO NOT add arbitrary clamps or hard-coded limits to the solver without explicit user approval.

- ALWAYS use the wrapper scripts (`scripts/run_model_grid.py`) to execute models.

- **WATCHDOGS MUST FAIL CLOSED:** Job submission scripts must abort with a non-zero exit code if state files are missing. Never fail open.

## 3. Directory Structure

```text
SN2024pba/
|-- AGENTS.md             # THIS FILE. Agent guidelines.
|-- configs/              # YAML model configurations
|-- scripts/              # Python orchestration, conversion, and plotting
|   |-- run_model_grid.py # MASTER ORCHESTRATOR.
|-- cmfgen_jfnk/          # *** ACTIVE DEVELOPMENT FOLDER ***
|                         # The bleeding-edge CMFGEN source with JFNK & non-mono solvers.
|-- cmfgen_validation/    # *** VALIDATION & TESTING FOLDER ***
|                         # Contains Luc Dessart's benchmark setups for code verification.
|-- cmfgen_source/        # Legacy/base CMFGEN Fortran source (reference only).
|-- cmfgen_base/          # CMFGEN base template files.
|-- atomic/               # Custom atomic data (f_to_s files).
|-- models/               # Model output directories (one per model).
|   |-- <model_name>/
|       |-- stella/       # STELLA hydro output (.swd)
|       |-- cmfgen/       # CMFGEN inputs/outputs
|       |-- plots/        # Centralized plots
|       |-- output/       # Final output files
|-- observational_data/   # Observed spectra and photometry.
|-- audits/               # Deep-dive code audits (READ THESE BEFORE TOUCHING FORTRAN).
|-- OLD/                  # Archived/deprecated models.
```

## 4. Active Development: JFNK & Non-Monotonic Velocities

**Context:** Standard CMFGEN assumes monotonic velocity fields. Interacting supernovae (e.g., SN 2024pba) feature CSM/ISM interactions that create non-monotonic velocity profiles (shocks, reverse shocks, dense shells). We are implementing new solvers to handle this.

### 4.1 New Feature Flags

You will be implementing and modifying logic controlled by three primary logical flags (likely defined in a module or header file like `cmfgen_mod.f` or similar):

1. `DO_JFNK`: Enables the Jacobian-Free Newton-Krylov solver for the radiative transfer/population coupling, replacing or augmenting the standard linearization/acceleration schemes.

2. `DO_NONMONO_SOB`: Enables a non-monotonic Sobolev-based escape probability / optical depth calculation.

3. `DO_NONMONO_KNOP`: Enables the Knop-style ordered refresh solver for non-monotonic fields (handling the complex frequency-angle coupling in compressing/expanding shocks).

4. Avoid using new VADAT paramters as possible or if they are needed allow them to be optional 

5. The code should be as streamlined as possible and not have unnecessary complexity.

### 4.2 Velocity Law Configuration

- For all interacting supernova tests and development, the velocity law parameter must be set to `VEL_LAW = 14`.
- Ensure that the `VEL_LAW=14` logic correctly feeds the non-monotonic velocity gradients into the `DO_NONMONO_*` solvers.

### 4.3 Coding Standards for `cmfgen_jfnk/`

- **Guard Clauses:** Every new block of code must be wrapped in `IF (DO_NONMONO_KNOP) THEN ... END IF` (or the respective flag). Do not alter the `ELSE` path.
- **Variable Declarations:** Follow the existing F77/F90 mixed conventions. If the surrounding code uses implicit typing, do not force `IMPLICIT NONE` in that specific subroutine unless it is a brand new file. Match the case (usually lowercase for variables, uppercase for parameters).
- **Comments:** Add clear, concise comment blocks explaining the physics and numerical rationale for the non-monotonic treatment, matching the style of the original CMFGEN comments.

## 5. Validation & Testing Strategy

You cannot consider a feature complete until it passes automated setup checks and physical validation.

- **Step 1: Python Setup Validation (Mandatory):** Before running any Fortran model, you MUST execute a Python validation script. This script must:
  1. Parse the target config/YAML and verify `VEL_LAW = 14` is set.
  2. Confirm that the `DO_JFNK`, `DO_NONMONO_SOB`, and `DO_NONMONO_KNOP` flags are correctly mapped to the Fortran namelists or input files.
  3. Verify that all required STELLA hydrodynamic input files (e.g., `.swd` and converted `SN_HYDRO_DATA`) exist and have valid, non-zero dimensions.

- **Step 2: Physical Validation Environment:** All physical validation runs must be executed in the `cmfgen_validation/` directory.

- **Step 3: The Benchmark:** The primary test case is Luc Dessart's interacting supernova spectrum and setup. This is the gold standard for verifying that the non-monotonic/JFNK code is set up correctly and producing physically meaningful results.

- **Step 4: Execution & Comparison:** Compile the code, run the Luc Dessart benchmark with the appropriate `DO_*` flags enabled, and compare the output spectrum and convergence metrics against the expected baseline. Only proceed to full model grids (`models/`) once the validation benchmark passes.

## 6. Core Pipeline Workflow

1. **STELLA:** Radiation hydrodynamics -> generates `.swd` file.
2. **Conversion:** `scripts/convert_stella_swd_to_cmfgen.py` converts `.swd` -> `SN_HYDRO_DATA`.
3. **CMFGEN:** `SN_HYDRO_DATA` -> non-LTE spectrum synthesis.

## 7. Debugging & Logs

When a job fails or stalls, check the following in this order:

1. **Python Orchestration Logs:** Check `models/<model>/logs/` for setup, conversion, or configuration parsing errors.
2. **SLURM Logs:** `slurm-<jobid>.out`. Look for OOM or segfaults.
3. **CMFGEN Logs:** Inside `models/<model>/cmfgen/` or `cmfgen_validation/`.
   - `OUTGEN`: General output and convergence metrics.
   - `KNDUMP` / `KNSTOR`: KNOP solver specific diagnostics (gated by `PRINT_KNOP_DEBUG`).
4. **State Files:** Check `GS_STATE` and checkpoint files. Ensure they are from the same iteration.
5. **Physics Checks:** If the solver stalls, check if the non-monotonic flags (`DO_NONMONO_*`) are actually being triggered for the specific depth/frequency cells experiencing the shock/reverse-shock.

## 8. Plotting Standards

All plots must be saved to `models/<model_name>/plots/` or `cmfgen_validation/plots/`.

**Naming Convention:** `{code}_{epoch}_{type}.png`

- Examples: `cmfgen_30d_spectrum.png`, `cmfgen_45d_spectropolarimetry.png`.
- Always include axis labels, units (e.g., Angstrom for wavelength, 10^-17 erg s^-1 cm^-2 Angstrom^-1 for flux), and a legend.
- When comparing to Luc Dessart's validation runs, explicitly overlay the benchmark spectrum for visual diffing.

## 9. SN_B X-Grid Geometry Convention

- The SN_A ejecta appended by `x_model_from_sn` must be kept fully
  intact and metal rich. Do NOT set `ejecta_r_inner_cm` or
  `ejecta_v_max_kms`: those options discard the Fe/Si core and the
  fast outer edge, leaving only the He layer.
- Reference config: `SN_B_X_grid/X_V10000_wind045i.yaml`.
- Wind(s) and the 4 Msun shell may be tuned freely; they occupy only
  the region between the shell surface and the SN_A inner edge
  (~7e15 cm). Overlays (`slow_shell_*`) must end below that edge so
  SN_A zones are never merged or diluted.

## 10. SEDONA Development Repository

- `sedona_source/` is its own git repo (nested; ignored by the workspace
  whitelist .gitignore). All SEDONA code development is tracked there.
- Remote: `git@github.com:Astro-Sean/pubsed.git` (user's SEDONA fork).
- `main` tracks upstream pubsed development.
- `development` is the active development branch for the
  interacting-SN modifications (non-monotonic velocity support, DDMC
  diagnostics, isotope-aware converters). Only two branches exist
  (main + development); ALL commits and pushes go to development
  by default (push.default=upstream, development is checked out).
- Build: `cd sedona_source/src && module load intel/2024.0 openmpi/4.1
  gsl/2.6 hdf5-serial/1.14.1`, set `SEDONA_HOME`, `GSL_DIR`,
  `HDF5_DIR`, `LUA_DIR=$SEDONA_HOME/src/external/lua-5.2.1`, then
  `cd build && make gomc CXX=mpicxx` (Makefile.general hardcodes
  /usr/local/bin/mpicxx; override with CXX=). Binary lands at
  `src/sedona6.ex` via install.sh, or build/gomc directly.
  NOTE: `transport.h` changes require `rm -f build/*.o` first -- the
  Makefile has no header dependencies and stale objects cause ABI
  segfaults.
  CRITICAL: `build/` holds COPIES of the sources, made once by
  install.sh. Running `make` in build/ compiles the stale copies.
  After ANY source edit, first sync:
  `cp -p transport/*.cpp transport/*.h opacity/*.cpp opacity/*.h
   utils/*.cpp utils/*.h grid/*.cpp grid/*.h hydro/*.cpp hydro/*.h
   hydro/*/*.cpp hydro/*/*.h main/*.cpp main/*.h build/`
  then `rm -f build/*.o && make gomc CXX=mpicxx`.
- `.mod` species tokens must be written as "Z.A" (e.g. "1.1 2.4 26.56");
  bare "Z" is silently parsed as A=Z by grid_1D_sphere::read_ascii_file.
