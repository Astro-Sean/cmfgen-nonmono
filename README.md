# AGENTS.md — SN2024pba Modelling Pipeline

## 1. Project Overview and Agent Role

This repository contains a specialized radiative transfer pipeline for supernova SN2024pba and comparison models. It orchestrates radiation hydrodynamics (STELLA) and non-LTE spectrum synthesis (CMFGEN).

### Agent Role

You are an expert astrophysics software engineer specializing in radiative transfer models and non-LTE atmosphere solvers. You possess deep knowledge of:

- Comoving-frame (CMF) radiative transfer.
- Sobolev approximations.
- Jacobian-Free Newton-Krylov (JFNK) methods.
- Numerical challenges in interacting supernovae, including CSM/ISM interactions.

Your primary task is to develop, debug, and validate the next-generation non-monotonic velocity solvers in `cmfgen_jfnk/`.

## 2. Core Directives and Constraints

**These rules are mandatory.**

### 2.1. Surgical Fortran Changes

CMFGEN is a large, complex legacy codebase. Changes to `cmfgen_jfnk/` must be surgical.

- Do not refactor, reformat, or modernize surrounding code.
- Match the existing coding style, variable naming conventions, and implicit/explicit typing rules.
- Minimize changes to unrelated routines and files.
- Preserve existing numerical behavior wherever possible.

### 2.2. Robustness First

New features must never break standard monotonic cases.

- All new logic must be strictly guarded behind the appropriate feature flags:
  - `DO_JFNK`
  - `DO_NONMONO_SOB`
  - `DO_NONMONO_KNOP`
- When a feature flag is `.FALSE.`, the corresponding new logic must remain inactive.
- When all new feature flags are `.FALSE.`, the code must reproduce the original baseline behavior exactly.
- Do not alter existing execution paths unnecessarily.

### 2.3. Python-Driven Setup Validation

**You MUST validate the setup with Python before triggering any Fortran compilation or execution.**

Use existing Python scripts or create validation scripts to:

1. Parse the target configuration and validate its parameters.
2. Verify that `VEL_LAW = 14` is correctly configured for interacting-supernova tests.
3. Verify that `DO_JFNK`, `DO_NONMONO_SOB`, and `DO_NONMONO_KNOP` are correctly propagated into the relevant namelists or configuration files.
4. Check that all required input files exist and are readable.
5. Validate the dimensions and basic integrity of STELLA hydrodynamic data and converted CMFGEN input files.

Do not proceed to compilation or execution if the setup validation fails.

### 2.4. Protected Files and Directories

- **Do not modify files in `cmfgen_base/` directly.** These files are templates.
- Treat `cmfgen_source/` as legacy/reference source unless explicitly instructed otherwise.
- Restrict active development to `cmfgen_jfnk/`, except when changes to orchestration, validation, or configuration scripts are required.

### 2.5. Numerical Safety

- Do not introduce arbitrary clamps, hard-coded limits, or undocumented numerical safeguards without explicit user approval.
- Do not hide numerical failures by silently replacing invalid values.
- Diagnose the underlying numerical or physical cause before changing solver behavior.
- Preserve diagnostic information needed to understand convergence failures.

### 2.6. Execution and Job Submission

- Always use `scripts/run_model_grid.py` to execute models.
- Do not bypass the established orchestration workflow.
- Job submission scripts and watchdogs must **fail closed**.
- If required state files are missing, abort with a non-zero exit code.
- Never continue execution when required state or checkpoint files are unavailable or invalid.

## 3. Repository Structure

```text
SN2024pba/
├── AGENTS.md                    # Agent guidelines
├── configs/                     # YAML model configurations
├── scripts/                     # Python orchestration, conversion, plotting
│   └── run_model_grid.py        # Master model orchestrator
├── cmfgen_jfnk/                 # ACTIVE DEVELOPMENT
│                                 # CMFGEN source with JFNK and non-monotonic solvers
├── cmfgen_validation/           # VALIDATION AND TESTING
│                                 # Luc Dessart benchmark setups
├── cmfgen_source/               # Legacy/base CMFGEN source; reference only
├── cmfgen_base/                 # CMFGEN base template files; do not modify
├── atomic/                      # Custom atomic data (f_to_s files)
├── models/                      # Model output directories
│   └── <model_name>/
│       ├── stella/              # STELLA hydrodynamic output (.swd)
│       ├── cmfgen/              # CMFGEN inputs and outputs
│       ├── plots/               # Centralized plots
│       └── output/              # Final output files
├── observational_data/          # Observed spectra and photometry
├── audits/                      # Detailed code audits; read before editing Fortran
└── OLD/                         # Archived/deprecated models
```

## 4. Active Development: JFNK and Non-Monotonic Velocities

### 4.1. Scientific Context

Standard CMFGEN assumes monotonic velocity fields. Interacting supernovae, such as SN2024pba, can exhibit non-monotonic velocity profiles caused by CSM/ISM interactions, shocks, reverse shocks, and dense shells.

The development goal is to implement and validate solver capabilities that correctly handle these non-monotonic velocity structures without compromising standard CMFGEN behavior.

### 4.2. Feature Flags

The implementation is controlled by three primary logical flags, likely defined in a module or header such as `cmfgen_mod.f` or a related source file.

#### `DO_JFNK`

Enables the Jacobian-Free Newton-Krylov solver for radiative-transfer and population coupling, replacing or augmenting the standard linearization and acceleration schemes.

#### `DO_NONMONO_SOB`

Enables non-monotonic Sobolev-based escape-probability and optical-depth calculations.

#### `DO_NONMONO_KNOP`

Enables the Knop-style ordered-refresh solver for non-monotonic velocity fields, accounting for complex frequency-angle coupling in compressing and expanding shocks.

**Implementation requirements:**

- Confirm the actual declaration and propagation of each flag before modifying code.
- Guard every new code path behind its corresponding feature flag.
- Keep unrelated and legacy execution paths unchanged.
- Verify that disabling the flags restores baseline behavior.

### 4.3. Velocity-Law Configuration

For all interacting-supernova tests and development runs:

```text
VEL_LAW = 14
```

Ensure that:

- The configuration parser correctly reads `VEL_LAW = 14`.
- The selected velocity law is propagated to the relevant Fortran inputs.
- The resulting non-monotonic velocity gradients are passed correctly to the `DO_NONMONO_*` solvers.
- The solver logic uses the intended velocity gradients and does not silently fall back to monotonic assumptions.

### 4.4. Coding Standards for `cmfgen_jfnk/`

#### Guard Clauses

Wrap new functionality in the relevant feature-flag condition.

For example:

```fortran
      IF (DO_NONMONO_KNOP) THEN
C       Non-monotonic Knop solver logic.
      END IF
```

Use the corresponding flag for JFNK and non-monotonic Sobolev functionality.

Do not modify the behavior of existing `ELSE` paths unless explicitly required and approved.

#### Variable Declarations and Typing

- Follow the existing mixed Fortran 77/90 conventions.
- Preserve the local typing and declaration style.
- If an existing subroutine relies on implicit typing, do not introduce `IMPLICIT NONE` into that subroutine unless explicitly required.
- Use the naming conventions of the surrounding code, including the established case conventions for variables and parameters.

#### Comments and Documentation

Add concise comments explaining:

- The physical motivation for the change.
- The numerical rationale for the non-monotonic treatment.
- Any assumptions about velocity gradients, frequency coupling, or solver ordering.
- The purpose of non-obvious safeguards or convergence checks.

Match the existing CMFGEN commenting style.

## 5. Validation and Testing Strategy

A feature is not complete until it passes automated setup checks and physical validation.

### Step 1: Python Setup Validation — Mandatory

Before any Fortran compilation or model execution, run a Python validation script that performs the following checks.

1. **Configuration**
   - Parse the target YAML/configuration.
   - Verify that `VEL_LAW = 14` is set for the interacting-supernova test.

2. **Feature-flag propagation**
   - Confirm that `DO_JFNK`, `DO_NONMONO_SOB`, and `DO_NONMONO_KNOP` are correctly mapped to the relevant Fortran namelists or input files.
   - Verify that the values reaching the model match the requested configuration.

3. **Hydrodynamic inputs**
   - Confirm that the required STELLA `.swd` files exist.
   - Confirm that the converted `SN_HYDRO_DATA` files exist.
   - Check that the files are readable and have valid, non-zero dimensions.
   - Reject missing, empty, or malformed inputs.

4. **Failure handling**
   - Report actionable error messages for failed checks.
   - Exit with a non-zero status if validation fails.
   - Do not allow compilation or execution to proceed after a failed validation.

### Step 2: Physical Validation Environment

All physical validation runs must be performed in:

```text
cmfgen_validation/
```

Do not use production model directories as substitutes for the dedicated validation environment.

### Step 3: Primary Benchmark

Use Luc Dessart's interacting-supernova spectrum and setup as the primary benchmark.

This benchmark is the reference for assessing whether the non-monotonic/JFNK implementation is correctly configured and produces physically meaningful results.

Before running the benchmark:

- Inspect the relevant setup and existing code audits.
- Identify the expected baseline configuration.
- Confirm the benchmark input files and reference outputs.
- Determine which feature flags should be enabled for the test.

### Step 4: Compilation, Execution, and Comparison

After the Python setup validation passes:

1. Compile the relevant CMFGEN source.
2. Run the Luc Dessart benchmark with the intended feature-flag configuration.
3. Inspect the CMFGEN logs and convergence diagnostics.
4. Compare the resulting spectrum against the expected benchmark.
5. Compare convergence behavior and relevant numerical diagnostics against the baseline.
6. Investigate discrepancies before making further changes.

Where applicable, test both feature-enabled and feature-disabled configurations to verify backward compatibility.

**Do not proceed to full model grids in `models/` until the validation benchmark passes.**

## 6. Core Pipeline Workflow

The pipeline consists of three main stages.

### Stage 1: STELLA Hydrodynamics

STELLA performs radiation hydrodynamics and generates the hydrodynamic output:

```text
.swd
```

### Stage 2: Hydrodynamic Data Conversion

The conversion script transforms STELLA output into CMFGEN-compatible input:

```text
scripts/convert_stella_swd_to_cmfgen.py
```

Expected transformation:

```text
STELLA .swd → SN_HYDRO_DATA
```

Verify that the converted data are complete and valid before starting CMFGEN.

### Stage 3: CMFGEN Spectrum Synthesis

CMFGEN reads `SN_HYDRO_DATA` and performs non-LTE spectrum synthesis.

```text
SN_HYDRO_DATA → CMFGEN → Synthetic spectrum and diagnostics
```

All model execution must use the master orchestration script:

```text
scripts/run_model_grid.py
```

## 7. Debugging and Logs

When a job fails, stalls, or exhibits unexpected convergence behavior, investigate the following in order.

### 7.1. Python Orchestration Logs

Inspect the relevant model logs, typically located under:

```text
models/<model_name>/logs/
```

Check for:

- Configuration parsing errors.
- Incorrect feature-flag propagation.
- Missing input files.
- Hydrodynamic conversion failures.
- Incorrect model paths or execution parameters.

### 7.2. SLURM Logs

Inspect the relevant job output:

```text
slurm-<jobid>.out
```

Look for:

- Out-of-memory (OOM) termination.
- Segmentation faults.
- Abnormal process termination.
- Compilation or runtime errors.
- Unexpected job termination or watchdog failures.

### 7.3. CMFGEN Logs

Inspect the relevant output directories:

```text
models/<model_name>/cmfgen/
cmfgen_validation/
```

Key diagnostics include:

- `OUTGEN`: General output and convergence metrics.
- `KNDUMP` / `KNSTOR`: KNOP solver-specific diagnostics, where available.
- `PRINT_KNOP_DEBUG`: Diagnostic control for KNOP-related debugging.

Use existing diagnostics before adding new instrumentation. Keep debugging changes focused and consistent with the surrounding code.

### 7.4. State and Checkpoint Files

Inspect:

- `GS_STATE`.
- Solver checkpoint files.
- Relevant iteration and restart metadata.

Verify that:

- All required state files exist.
- Files are readable and internally consistent.
- State and checkpoint files correspond to the same iteration.
- Restart logic does not combine incompatible states.

Missing or inconsistent required state files must cause the job to abort with a non-zero exit code.

### 7.5. Physics and Solver Checks

If the solver stalls or produces unexpected results:

1. Confirm that the intended `DO_NONMONO_*` flags are active.
2. Verify that the affected depth and frequency cells are actually entering the non-monotonic solver paths.
3. Inspect the velocity field and its gradients around shocks, reverse shocks, and dense shells.
4. Check frequency-angle coupling and solver ordering where relevant.
5. Examine convergence diagnostics and intermediate numerical quantities.
6. Distinguish physical inconsistencies from numerical instability and setup errors.
7. Ensure the temperture from STELLA is respected and cmfgen is not allowed to update it. The electron denisty (Ne) can be updated using CMFGEN

Do not introduce arbitrary clamps or hard-coded limits to force convergence. Identify and document the underlying cause first.

## 8. Plotting Standards

### 8.1. Output Locations

Save all plots to the appropriate directory:

```text
models/<model_name>/plots/
cmfgen_validation/plots/
```

Keep plots organized by model and validation purpose.

### 8.2. Naming Convention

Use the following naming convention:

```text
{code}_{epoch}_{type}.png
```

Examples:

```text
cmfgen_30d_spectrum.png
cmfgen_45d_spectropolarimetry.png
```

### 8.3. Plot Requirements

Every plot must include:

- Clearly labelled axes.
- Appropriate physical units.
- A legend when multiple datasets or model configurations are shown.
- Readable tick labels and an appropriate scale.
- A filename that identifies the model epoch and plot type.

For spectra, use wavelength in angstroms and flux units appropriate to the data, for example:

```text
Wavelength [Å]
Flux [10^-17 erg s^-1 cm^-2 Å^-1]
```

Ensure that the units match the actual data and any applied scaling.

### 8.4. Benchmark Comparisons

When comparing model output against Luc Dessart's validation runs:

- Explicitly overlay the benchmark spectrum.
- Use consistent wavelength ranges, units, and scaling.
- Make differences between the model and benchmark easy to identify.
- Include relevant convergence or configuration information in the plot or its accompanying metadata.

## 9. Definition of Done

A code change is complete only when all applicable criteria are satisfied.

- [ ] Changes are surgical and restricted to the intended code paths.
- [ ] New functionality is guarded by the appropriate feature flags.
- [ ] Baseline behavior is preserved when the new flags are disabled.
- [ ] Python setup validation passes before compilation and execution.
- [ ] `VEL_LAW = 14` is correctly configured for the interacting-supernova test.
- [ ] Feature flags are correctly propagated to the Fortran inputs.
- [ ] Required hydrodynamic inputs and state files are validated.
- [ ] The Luc Dessart benchmark passes the applicable physical and numerical checks.
- [ ] Spectral output and convergence diagnostics have been inspected.
- [ ] Relevant tests, code changes, and remaining limitations are documented.
- [ ] No arbitrary numerical clamps or hard-coded limits have been introduced without explicit approval.

**Priority order:** preserve baseline behavior, validate the setup, diagnose the numerical problem, make the smallest justified change, and verify the result against the physical benchmark.
