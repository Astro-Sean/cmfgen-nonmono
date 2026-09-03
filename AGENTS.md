````markdown
# AGENTS.md — SN2024pba Two-Event Supernova Modelling Pipeline

## 1. Purpose

This repository models a two-event interacting-supernova scenario:

- **SN_A**: an earlier Type IIb explosion producing metal-rich ejecta.
- **SN_B**: a later explosion/shell interaction with the older SN_A ejecta.

The scientific objective is to determine whether heating and ionisation of the metal-rich SN_A ejecta by the SN_B interaction can reproduce the observed spectrum, particularly strong absorption associated with:

- H
- He
- O
- Si
- Ca
- Fe

The repository contains a coupled STELLA → CMFGEN → SEDONA/SuperLite modelling workflow.

This document is the authoritative guidance for AI agents working on the repository.

Historical notes, old prompts, copied code comments and previous `AGENTS.md` files are **not authoritative** if they conflict with the current source code, configuration, input files or reproducible tests.

---

# 2. Core Execution Rule — `run_model_grid.py` Is the Single Entry Point

## 2.1 Mandatory rule

**Every model must be prepared, validated and submitted through:**

```text
scripts/run_model_grid.py
````

Agents must not normally submit individual STELLA, CMFGEN or SuperLite jobs directly.

Do not bypass the orchestrator with:

```bash
sbatch ...
```

or by manually invoking:

```text
submit_cmfgen.sh
batch.sh
superlite_job_template.slurm
stella_job_template.slurm
```

unless performing a specifically documented diagnostic test.

The orchestrator is responsible for:

1. resolving the YAML configuration;
2. determining which codes should run;
3. determining which existing outputs can be reused;
4. preparing model directories;
5. converting input data;
6. generating CMFGEN inputs;
7. constructing `MODEL_SPEC`;
8. constructing `VADAT`;
9. constructing `IN_ITS` / iteration controls where required;
10. constructing `SN_HYDRO_DATA`;
11. constructing SuperLite inputs;
12. selecting executables;
13. creating job scripts;
14. establishing job dependencies;
15. submitting jobs;
16. recording provenance and configuration;
17. collecting diagnostics.

The intended execution architecture is:

```text
YAML CONFIG
    |
    v
run_model_grid.py
    |
    +---- STELLA preparation/submission
    |          |
    |          v
    |       .swd
    |
    +---- STELLA .swd
    |          |
    |          v
    |     SN_HYDRO_DATA
    |          |
    |          +------------------+
    |          |                  |
    |          v                  v
    |       CMFGEN             SuperLite
    |          |                  |
    |          v                  v
    |       spectrum           spectrum
    |
    +---- SEDONA preparation/submission
```

## 2.2 Existing STELLA models may be reused

STELLA does **not** need to be rerun if a suitable existing `.swd` model already exists.

The orchestrator must support:

```text
run STELLA = yes
run STELLA = no / reuse existing model
```

However, even when STELLA is skipped, **all downstream preparation must still occur through `run_model_grid.py`**.

For example:

```bash
python scripts/run_model_grid.py configs/X1_fix.yaml --no-stella --run-cmfgen
```

must:

1. locate the existing `.swd`;
2. validate it;
3. convert it appropriately;
4. construct the required `SN_HYDRO_DATA`;
5. optimize CMFGEN inputs;
6. generate the CMFGEN run directory;
7. submit the CMFGEN job.

An existing model must never be assumed to be compatible merely because its filename exists.

---

# 3. Blind Initial Setup

The model preparation system must be **initially blind to the detailed input model**.

The setup code must not contain hard-coded assumptions such as:

```text
this model always uses H I–II
this model always uses Fe II–IV
this model always has these ionisation stages
this model always uses these MODEL_SPEC levels
```

Instead:

```text
INPUT HYDRO MODEL
        |
        v
ANALYSIS
        |
        +-- density structure
        +-- temperature
        +-- radiation field
        +-- velocity
        +-- composition
        +-- abundance gradients
        +-- optical-depth structure
        +-- ionisation proxies
        |
        v
MODEL-SPEC OPTIMISATION
        |
        v
VADAT / IN_ITS / CMFGEN INPUT OPTIMISATION
        |
        v
SN_HYDRO_DATA VALIDATION
        |
        v
CMFGEN
```

The preparation code must adapt to the supplied model.

Do not hard-code a single atomic configuration as the universal solution.

---

# 4. Repository Structure

The existing repository structure should be preserved unless there is a strong technical reason to change it.

Expected structure:

```text
SN2024pba/

├── configs/
│   └── *.yaml

├── scripts/
│   ├── run_model_grid.py
│   ├── run_model_grid.sh
│   ├── convert_*.py
│   ├── prepare_VADAT.py
│   ├── submit_cmfgen_jobs.sh
│   ├── plot_*.py
│   └── deprecated/

├── cmfgen_base/
│   ├── MODEL_SPEC
│   ├── VADAT_INIT
│   ├── RDINR
│   ├── T_IN
│   ├── IN_ITS
│   ├── batch.sh
│   └── submit_cmfgen.sh

├── cmfgen_source/

├── atomic/

├── sedona_base/

├── superlite_base/

├── sirocco_base/

├── models/
│   └── <model>/
│       ├── stella/
│       ├── cmfgen/
│       ├── sedona/
│       ├── superlite/
│       ├── plots/
│       └── output/

├── outputs/
├── observational_data/
├── presn_models/
├── homo_model/
├── docs/
├── OLD/
├── SN_12C*/
└── AGENTS.md
```

Generated run products should remain inside model/run directories.

Do not pollute source directories with:

* generated `SN_HYDRO_DATA`;
* spectra;
* temporary files;
* SLURM output;
* restart files;
* generated `MODEL_SPEC`;
* generated `VADAT`;
* generated `IN_ITS`;
* compiled executables;
* intermediate conversion files.

---

# 5. Evidence Policy

Every technical conclusion must be classified as:

```text
VERIFIED
SUPPORTED
UNVERIFIED
CONTRADICTED
BLOCKED
```

### VERIFIED

Demonstrated using:

* current source code;
* an automated test;
* a reproducible model run;
* direct inspection of generated input/output.

### SUPPORTED

Consistent with:

* CMFGEN documentation;
* published literature;
* known physical behaviour;

but not independently reproduced.

### UNVERIFIED

Plausible but not tested.

### CONTRADICTED

Inconsistent with current source code, input files, equations or reproducible results.

### BLOCKED

Cannot be tested because required:

* source;
* executable;
* model;
* atomic data;
* compute resource;
* input file

is unavailable.

Never convert an `UNVERIFIED` result into a `VERIFIED` result simply because it appears reasonable.

---

# 6. Model Timeline

Explicitly distinguish:

```text
t_A
t_A_snapshot
t_B
t_B_SBO
t_B_snapshot
t_B_end
```

Nominal setup:

```text
t_A = SN_A explosion

t_A_snapshot = t_A + 30 d

t_B ≈ t_A + 400 d

t_B_SBO = SN_B shock breakout

t_B_snapshot = t_B_SBO + 40 d

t_B_end = t_B_SBO + 60 d
```

Do not silently equate:

* explosion;
* shell launch;
* shock arrival;
* shock breakout.

The actual definitions must be stored in YAML configuration and run metadata.

---

# 7. STELLA

STELLA provides the hydrodynamic/radiation-hydrodynamic structure.

The relevant `.swd` model is the primary source for the downstream structure.

Validate:

* radius;
* velocity;
* density;
* gas temperature;
* radiation temperature;
* luminosity;
* composition;
* mass coordinate;
* optical-depth structure;
* shock positions;
* units;
* time.

If an existing `.swd` is reused, perform the same validation as for a newly generated model.

Never assume that an old `.swd` is valid because it previously produced a spectrum.

---

# 8. `SN_HYDRO_DATA` Is a Critical Interface

`SN_HYDRO_DATA` is not merely a format-conversion file.

It is the **physical interface between the hydrodynamic model and CMFGEN**.

The conversion must preserve the structure required by CMFGEN.

Validate:

```text
radius
velocity
density
temperature
radiation field
luminosity
composition
mass coordinate
time
```

Check:

* ordering;
* units;
* signs;
* monotonicity where required;
* duplicate points;
* extremely small/large gradients;
* shock locations;
* abundance conservation;
* mass conservation.

The conversion must not blindly copy every STELLA grid point.

---

# 9. CMFGEN Grid Optimisation

## 9.1 Fundamental rule

The CMFGEN grid must be constructed to enable robust convergence while preserving the physical structure required by the transfer calculation.

The grid must be:

> **as small as possible while still satisfying all convergence and physical-resolution constraints.**

Maximum allowed grid size:

```text
1000 points
```

The maximum of 1000 is a ceiling, **not a target**.

A model using 1000 points when 300 points are sufficient is considered poorly optimised.

---

# 10. Tau-Based AMR Grid Construction

The CMFGEN conversion should use a **tau-based adaptive mesh refinement strategy**.

The grid should be refined according to physically important quantities rather than simply inheriting the STELLA grid.

Important refinement indicators include:

```text
optical depth
velocity gradient
density gradient
temperature gradient
radiation-temperature gradient
composition gradients
shock positions
abundance discontinuities
ionisation-sensitive regions
line-forming regions
photospheric region
outer optically thin region
```

The algorithm should:

1. start from the STELLA structure;
2. calculate diagnostic gradients;
3. identify sharp physical transitions;
4. identify the relevant optical-depth range;
5. preserve shocks and abundance interfaces;
6. preserve line-forming regions;
7. coarsen smooth regions;
8. test convergence;
9. iteratively refine only where necessary;
10. stop when convergence criteria are satisfied.

Do **not** use the CMFGEN source code itself to determine the grid.

CMFGEN source inspection can identify numerical requirements, but the STELLA → CMFGEN conversion must determine the grid from the physical structure and explicit convergence criteria.

---

# 11. Grid Requirements

The grid must preserve:

### Hydrodynamic structure

* velocity;
* density;
* mass;
* radius.

### Thermal structure

* gas temperature;
* radiation temperature.

### Composition

At minimum:

```text
H
He
C
N
O
Ne
Na
Mg
Al
Si
S
Ar
Ca
Fe
Ni
```

where available.

### Transfer structure

Preserve:

* photosphere;
* optical-depth transitions;
* shocks;
* composition interfaces;
* strong velocity gradients;
* line-forming regions.

A single point must never represent a physically broad discontinuity if doing so destroys the relevant gradient.

---

# 12. Grid Convergence Test

Grid optimisation must be treated as a convergence problem.

For candidate grids:

```text
N = small
N = intermediate
N = larger
```

compare:

* temperature structure;
* electron density;
* ionisation;
* luminosity;
* emergent spectrum;
* key line equivalent widths;
* line velocities;
* continuum shape.

A grid is accepted only when further refinement produces changes below configurable tolerances.

The optimisation target is:

```text
MINIMUM N
subject to
PHYSICAL ACCURACY + CMFGEN CONVERGENCE
```

not:

```text
N = 1000
```

---

# 13. CMFGEN Model Atom Optimisation

The atomic model must be selected from the physical conditions in `SN_HYDRO_DATA`.

The following elements are particularly important:

```text
H
He
O
Si
Ca
Fe
```

The model must not blindly include every available ionisation stage.

Nor should it remove ions merely because they are difficult to converge.

The atomic configuration should be determined from:

```text
temperature
density
radiation temperature
electron density
composition
radius
velocity
optical depth
ionisation proxies
```

---

# 14. Ionisation-Stage Selection

For each important element, estimate the potentially populated ionisation stages across the model.

For example:

```text
H
    H I
    H II

He
    He I
    He II
    He III

O
    O I
    O II
    O III
    O IV
    ...

Si
    Si I
    Si II
    Si III
    Si IV
    ...

Ca
    Ca I
    Ca II
    Ca III
    Ca IV
    ...

Fe
    Fe I
    Fe II
    Fe III
    Fe IV
    Fe V
    ...
```

The actual range must be determined from the model.

Use:

* local temperature;
* radiation temperature;
* density;
* composition;
* approximate ionisation balance;
* previous CMFGEN iterations;
* diagnostic populations.

Do not assume gas temperature alone determines ionisation.

In SN models, the radiation field can differ substantially from the gas temperature.

---

# 15. Important Rule for Low-Abundance Species

A species with an extremely small mass fraction may have negligible impact on the global energy balance while still being spectroscopically important.

Therefore:

```text
low abundance != automatically removable
```

Before removing an ion/species, determine:

1. whether it contributes important opacity;
2. whether it contributes an observed line;
3. whether it influences electron density;
4. whether it affects neighbouring ionisation stages;
5. whether it is required as a recombination pathway;
6. whether removing it improves convergence without altering the spectrum.

A convergence problem must not automatically be solved by deleting the problematic species.

---

# 16. H / He / O / Si / Ca / Fe Priority

For this project, the atomic setup must explicitly optimise:

```text
H
He
O
Si
Ca
Fe
```

These elements are scientifically important and must receive higher priority than minor species.

For each one, generate an automated report containing:

```text
element
included ions
maximum/minimum estimated ion fraction
maximum/minimum abundance
temperature range
radiation-temperature range
density range
optical-depth range
important transitions
convergence behaviour
```

---

# 17. Atomic-Level Selection

The atomic model must distinguish:

```text
ionisation stages
+
explicit levels
+
superlevels
+
transitions
```

Do not simply select all levels available in the atomic database.

The objective is:

```text
minimum atomic complexity
consistent with
accurate ionisation balance
+
important line formation
+
electron density
+
energy balance
+
convergence
```

For each important element:

1. identify ions likely to exist;
2. identify ions adjacent to the dominant ion;
3. retain sufficient levels to model recombination/cascade pathways;
4. retain lower levels associated with observed transitions;
5. retain upper levels required to populate those transitions;
6. include enough levels to maintain realistic ionisation balance;
7. remove demonstrably irrelevant levels only after testing.

---

# 18. CMFGEN `MODEL_SPEC`

`MODEL_SPEC` must be generated or modified by the model-preparation pipeline.

Do not manually maintain a universal `MODEL_SPEC` for every model.

The preparation code should:

```text
SN_HYDRO_DATA
    ↓
physical analysis
    ↓
candidate ions
    ↓
candidate levels
    ↓
spectroscopic priorities
    ↓
convergence constraints
    ↓
MODEL_SPEC
```

The generated `MODEL_SPEC` must record why each important ion was included.

Where possible, save:

```text
MODEL_SPEC
MODEL_SPEC.diagnostics
MODEL_SPEC.source
MODEL_SPEC.hash
```

---

# 19. Atomic Data and `MODEL_SPEC` Settings

Atomic-data configuration must be treated as a coupled system.

Do not optimise:

```text
MODEL_SPEC
```

independently of:

```text
VADAT
IN_ITS
RDINR
T_IN
transition settings
dielectronic settings
population settings
```

The atomic configuration and numerical controls must be mutually consistent.

When changing the atomic model, inspect all related settings.

---

# 20. `TRANS_*` Settings

The model currently contains settings such as:

```text
TRANS_HI
TRANS_HeI
TRANS_OI
TRANS_O2
TRANS_CaI
TRANS_Ca2
TRANS_FeI
TRANS_Fe2
TRANS_FeIII
TRANS_FeIV
```

The agent must determine these settings from the actual atomic model.

Do not enable an ionisation stage in `TRANS_*` that is not represented consistently in the atomic data.

Likewise, do not leave important transitions enabled if the required levels are absent.

The pipeline should validate:

```text
TRANS_* ion
        ↕
MODEL_SPEC ion
        ↕
atomic data
        ↕
available levels/transitions
```

and report inconsistencies before submission.

---

# 21. Dielectronic Recombination Settings

Settings such as:

```text
DIE_HI
DIE_HeI
DIE_OI
DIE_O2
DIE_CaI
DIE_Ca2
DIE_Fe2
DIE_FeIII
DIE_FeIV
```

must be interpreted in the context of the actual atomic model.

Do not enable dielectronic treatment for an ion whose atomic data does not support the corresponding process.

Do not disable an important recombination pathway merely because convergence improves.

Any change must be physically justified and tested.

---

# 22. `FIX_*` Population Controls

The `FIX_*` controls must not be used as an arbitrary convergence hack.

Normally:

```text
FIX = 0
```

unless a specific population is deliberately fixed for a documented reason.

If a population is fixed:

* document why;
* identify the physical approximation;
* demonstrate its effect;
* test whether it changes the emergent spectrum.

A fixed population must never be used to hide an unstable atomic model.

---

# 23. `FIX_T`

Temperature fixing must be treated separately from population fixing.

If:

```text
FIX_T = T
```

then verify that this is physically intended for the SN_HYDRO_DATA input.

The gas temperature from the hydrodynamic model and CMFGEN temperature solution must not be confused.

The pipeline must explicitly document:

```text
T_gas
T_rad
CMFGEN solved temperature
```

where available.

---

# 24. Initial Population Strategy

Initial populations should be estimated from `SN_HYDRO_DATA`.

The preparation code should use:

```text
T_gas
T_rad
density
electron density estimate
composition
radiation field
```

to generate physically reasonable initial populations.

A useful hierarchy is:

```text
LTE / Saha estimate
        ↓
radiation-temperature correction
        ↓
dilution correction
        ↓
previous converged model if available
        ↓
CMFGEN iterative solution
```

Do not claim that Saha-LTE populations are the final physical solution.

They are initial estimates.

---

# 25. Population Scaling

If:

```text
POP_SCALE = T
```

or equivalent population-scaling behaviour is used, verify that it is consistent with the generated initial populations.

Population scaling must not:

* create negative populations;
* violate conservation;
* generate impossible ion fractions;
* artificially force convergence.

---

# 26. CMFGEN Numerical Settings and Atomic Model

The following settings must be considered together with the atomic model:

```text
SCALE_OPT
SOL_METH
EPS_TERM
MAX_LIN
MAX_LAM
MAX_CHNG
N_FIX_BA
BA_CHK_FAC
FIX_BA
NUM_LAM
DO_NG
ITS/NG
DO_LAM_IT
DO_AV
```

The agent must not blindly alter these settings.

If convergence is poor:

1. identify which species/ion is failing;
2. identify which depths are failing;
3. inspect populations;
4. inspect abundance;
5. inspect radiation field;
6. inspect temperature;
7. inspect atomic completeness;
8. inspect transition selection;
9. only then modify numerical controls.

---

# 27. Convergence Diagnosis

When CMFGEN fails to converge, determine whether the problem is:

```text
A. physical
B. atomic
C. grid-related
D. initial-population-related
E. numerical
F. input-format-related
```

Do not immediately increase iteration counts.

The diagnostic workflow should be:

```text
FAILED CONVERGENCE
        |
        +--> check input validity
        |
        +--> check grid
        |
        +--> check composition
        |
        +--> check T / T_rad
        |
        +--> check problematic ion
        |
        +--> check atomic completeness
        |
        +--> check initial populations
        |
        +--> check transition/dielectronic settings
        |
        +--> check numerical controls
```

---

# 28. `SN_HYDRO_DATA` → `MODEL_SPEC` Optimisation

The agent must implement a repeatable procedure:

```text
1. Read SN_HYDRO_DATA.

2. Determine physical ranges.

3. Determine abundance ranges.

4. Determine approximate ionisation ranges.

5. Identify important ions.

6. Identify important spectral transitions.

7. Determine required atomic levels.

8. Determine adjacent ionisation stages required for ionisation balance.

9. Construct candidate MODEL_SPEC.

10. Check atomic-data availability.

11. Validate TRANS_* settings.

12. Validate DIE_* settings.

13. Generate initial populations.

14. Run a reduced convergence test.

15. Inspect convergence.

16. Expand atomic model only where required.

17. Re-run.

18. Reduce unnecessary levels/ions if convergence and spectrum remain unchanged.

19. Save final optimised MODEL_SPEC and diagnostics.
```

The final model must be the **smallest scientifically adequate model**, not the smallest model that happens to converge.

---

# 29. CMFGEN Grid and Atomic Model Are Coupled

Grid optimisation and atomic-model optimisation must not be treated independently.

A grid that is adequate for continuum transfer may be inadequate for:

* sharp abundance transitions;
* strong velocity gradients;
* line formation;
* ionisation fronts;
* shocks.

Likewise, an atomic model that converges on a coarse grid may fail when the physically correct grid is restored.

The workflow should therefore iterate:

```text
GRID
 ↓
ATOMIC MODEL
 ↓
CONVERGENCE
 ↓
SPECTRUM
 ↓
GRID REFINEMENT
 ↓
ATOMIC MODEL REFINEMENT
```

until stable.

---

# 30. Spectral Validation

A successful CMFGEN run is not automatically a scientifically successful model.

Validate:

* continuum;
* H;
* He;
* O;
* Si;
* Ca;
* Fe;
* line velocities;
* equivalent widths;
* line optical depths;
* formation regions.

For important absorption features, determine:

```text
species
ion
transition
formation radius
velocity
temperature
density
abundance
optical depth
```

A visually similar line is insufficient to establish line identification.

---

# 31. SN_A / SN_B Component Attribution

Where possible, compare:

```text
SN_A only
SN_B / CSM only
SN_A + SN_B interaction
```

and perform species-ablation experiments.

For example:

```text
remove Si
remove Fe
remove Ca
remove O
```

and determine how the spectrum changes.

This establishes whether observed absorption is genuinely produced by heated SN_A material.

---

# 32. SuperLite

SuperLite supports two conversion paths.

### Legacy

```text
SN_HYDRO_DATA
    ↓
convert_sn_hydro_to_superlite.py
    ↓
input.str
```

### Preferred MC-aware path

```text
STELLA .swd
    ↓
convert_stella_swd_to_superlite.py
    ↓
input.str
```

The preferred path is:

```yaml
superlite:
  input_source: "swd"
```

The direct `.swd` path should be used where appropriate because it preserves STELLA's native opacity information and can perform MC-aware optical-depth truncation.

---

# 33. SuperLite MC-Aware Grid

The direct STELLA → SuperLite converter should use:

* native Rosseland opacity;
* electron-scattering floor;
* corrected bound-free/free-free opacity;
* hybrid optical depth;
* tau-based truncation;
* tau-weighted grid reduction.

Default concept:

```text
truncate_tau ≈ 20
max_cells ≈ 150
min_cells ≈ 20
```

These are defaults, not immutable scientific constants.

Any change must be validated.

---

# 34. Atomic Data Bug Prevention

Avoid variable shadowing in opacity calculations.

For example, do not use:

```python
Z
```

simultaneously for:

```text
metal mass fraction
```

and:

```text
atomic number
```

Use explicit names such as:

```python
Z_mass_fraction
Z_atom
```

This is particularly important in approximate opacity calculations.

---

# 35. Provenance

For every external model/code/data file record:

```text
origin
URL/repository
version
commit/tag
download date
SHA256
local path
units
coordinate convention
conversion script
modifications
validation status
```

Never silently substitute a different progenitor or atomic dataset.

---

# 36. Run Manifest

Every model run must contain a manifest with:

```text
model name
epoch
resolved YAML configuration
STELLA source
SN_HYDRO_DATA checksum
MODEL_SPEC checksum
VADAT checksum
IN_ITS checksum
CMFGEN executable checksum
SuperLite executable checksum
SEDONA executable checksum
atomic-data version
command line
scheduler/job ID
resources
status
outputs
```

---

# 37. Required Dry-Run Workflow

Before submitting a production run:

```bash
python scripts/run_model_grid.py configs/<model>.yaml --validate-only
```

then:

```bash
python scripts/run_model_grid.py configs/<model>.yaml --dry-run
```

then, where supported:

```bash
python scripts/run_model_grid.py configs/<model>.yaml --interactive-test
```

Only after these succeed should production submission occur.

---

# 38. Testing

At minimum test:

```text
configuration parsing
STELLA loading
SN_HYDRO_DATA parsing
unit conversion
radial ordering
mass conservation
elemental-mass conservation
kinetic energy
abundance normalization
shock tracking
velocity gradients
velocity reversals
zero gradients
sigma calculation
zero flux
VEL_LAW=12 propagation
MODEL_SPEC generation
atomic-data consistency
TRANS_* consistency
DIE_* consistency
initial populations
grid refinement
CMFGEN convergence
spectrum reproducibility
clean-build reproducibility
restart reproducibility
line attribution
```

---

# 39. `VEL_LAW = 12`

Use:

```text
VEL_LAW = 12
```

only after verifying that the current CMFGEN source accepts the intended arbitrary velocity profile.

Demonstrate that:

* supplied velocity values are read;
* units are preserved;
* signs are preserved;
* ordering is preserved;
* CMFGEN does not silently overwrite them;
* the final transfer calculation uses the supplied profile.

Do not silently monotonise a physical velocity profile.

---

# 40. Velocity Reversals and Shocks

If the hydrodynamic structure contains velocity reversals or near-zero velocity gradients:

* identify them explicitly;
* preserve them if physically required;
* inspect all calculations involving `dv/dr`;
* avoid divisions by zero;
* do not apply a global artificial gradient.

Any numerical floor must be:

```text
local
documented
tested
physically justified
```

---

# 41. Sigma Handling

Audit:

* sigma definition;
* finite-difference method;
* units;
* interpolation;
* floors;
* clamps;
* call sites.

Do not overwrite a shared sigma array with a numerical safeguard if another routine requires the physical value.

Test:

```text
smooth velocity
velocity reversal
shock
near-zero gradient
target model
```

---

# 42. Zero Flux

Distinguish:

```text
physical zero flux
```

from:

```text
missing luminosity
invalid boundary condition
uninitialised array
NaN
Inf
underflow
```

Do not introduce a generic epsilon floor without identifying the physical/numerical cause.

---

# 43. Minimal Source Changes

CMFGEN source changes must be:

* minimal;
* reviewable;
* reversible;
* documented.

Preserve normal CMFGEN behaviour when:

```text
VEL_LAW != 12
```

Do not modify unrelated solver components merely to obtain convergence.

Prefer:

```text
input preparation
configuration
MODEL_SPEC optimisation
grid optimisation
diagnostics
```

before modifying Fortran source.

---

# 43.5. Convergence Control Files (IT_SPECIFIER, ADJUST_CORRECTIONS, ADJUST_R_DEFAULTS)

CMFGEN reads three optional runtime control files every iteration that can
facilitate convergence without modifying Fortran source code. These are
documented in `docs/full_descr.txt` sections "IT SPECIFIER" and
"ADJUST CORRECTIONS".

## IT_SPECIFIER

Forces a particular iteration sequence (LAM, FULL, FULL(FIXT), FULL(FIXBA),
DO_NG, NORM). Read every iteration by `SPECIFY_IT_CYCLE_V3`
(`$cmfdist/subs/specify_it_cycle_v3.f`). Can be modified while the model
is running.

Default location: `cmfgen_base/IT_SPECIFIER` (contains `NORM`).

YAML config: `cmfgen.iteration_strategy` — a list of strings, each line
is one option in the IT_SPECIFIER file. Example:

```yaml
cmfgen:
  iteration_strategy:
    - FULL(FIXT) 5
    - DO_NG
    - FULL(FIXT) 5
    - DO_NG
    - NORM
```

`prepare_VADAT.py` writes the file with the YAML-configured strategy
on the compute node. **If `iteration_strategy` is absent, IT_SPECIFIER is
deleted** so that `SPECIFY_IT_CYCLE_V3` returns immediately (line 61:
`IF IOS .NE. 0 RETURN`) and `DO_LAM_AUTO` controls the iteration mode.

### Revised approach (2026-09-02): DO_LAM_AUTO-only, no IT_SPECIFIER

**VERIFIED (2026-09-02):** The IT_SPECIFIER approach was abandoned because
it forced Full BA iterations too early (at iteration 12 when MAXCH was
still >500%). The 500–2000% corrections were physically meaningless —
`SUCCESS=T` only meant the matrix factorized, not that the correction
was valid.

The current approach uses CMFGEN's built-in auto-switch (`DO_LAM_AUTO=T`)
which is the standard approach used by Dessart and Hillier reference models:

1. Lambda iterations run until MAXCH < 50%
2. `DO_LAM_AUTO` auto-switches `DO_LAM_IT=F` in IN_ITS
3. Full BA iterations begin with `LAM_VAL=400` / `NUM_LAM=2` fallback
4. If BA overshoots, `solve_for_pops.f` automatically returns to Lambda

To use this approach, omit `iteration_strategy` from the YAML config.

### Zone-dependent BA for non-monotonic models (2026-09-02)

When `NON_MONO_SOLVER=T`, the BA variation uses a zone-dependent approach:
- `VAR_MOM_J_CMF_V12` is called with clamped SIGMA (>= -0.999)
- After the call, TX/TVX are zeroed at compression depths (SIGMA_TRUE < -1)
- Full cross-depth J-coupling at monotonic depths, fixed-J at compression depths
- K_ON_J clamps are applied before `QFROMF` to prevent division by zero

### Corrected VADAT settings for non-monotonic models

The following vadat_overrides are required:

```yaml
vadat_overrides:
  FIX_NE: F           # Ne must converge self-consistently
  FIX_T_AUTO: F       # Static snapshot: T fixed to hydro
  DO_LAM_AUTO: T      # Auto-switch Lambda→BA when MAXCH < 50%
  DO_LAM_IT: T        # Start in Lambda mode (auto-switched by DO_LAM_AUTO)
  DO_NG: T            # Ng acceleration essential for DIAG solver
  COMP_BA: T          # BA corrections enabled
  SOL_METH: DIAG      # Depth-independent solver
  NUM_LAM: '2'        # Manual: "Typically choose 2"
  LAM_VAL: 4.0D2      # Fallback to Lambda if MAXCH > 400
  BEG_NG: 1.0D1       # Manual: "5% to 10%" (was 1.0D15 = disabled)
  IBEG_NG: '15'       # Manual: "at least 15 iterations"
  NON_MONO_SOLVER: T  # Non-monotonic formal solver
  USE_LAM_ES: T       # J from formal solution, not moment equations
```

**Do not set `iteration_strategy`** — let `DO_LAM_AUTO` control the mode.

### Requeue configuration

Production runs must use the `general` partition with 24h wall time
and requeue enabled:

```yaml
cmfgen:
  partition: general
  time: '24:00:00'
  final_requeue_when: wall
  final_requeue_max: 500
```

The `interactive` partition only accepts jobs < 2h and is suitable only
for debug runs (hardcoded to 01:30:00). Setting `final_requeue_max: 0`
causes the job to abort immediately on the first requeue, wasting all
convergence progress.

## ADJUST_CORRECTIONS

Provides relaxation parameters (RELAX, T_LIM, MAX_CHNG, CONSIS_CNT) to
scale population corrections at specific depth ranges. Read every
iteration by `fiddle_pop_corrections_v2.f`. Only affects BA (full)
iterations, NOT Lambda iterations.

Default location: `cmfgen_base/ADJUST_CORRECTIONS` (no scaling effect).

YAML config: `cmfgen.adjust_corrections` — a dict with keys:
`l_st`, `l_end`, `relax`, `t_lim`, `max_chng`, `consis_cnt`. Example:

```yaml
cmfgen:
  adjust_corrections:
    l_st: 1
    l_end: 394
    relax: 0.5        # under-relaxation to damp oscillations
    t_lim: 0.2        # max 20% T change per iteration
    max_chng: 3.0     # max 3x population increase
    consis_cnt: 0     # no consistency check
```

`prepare_VADAT.py` overwrites the file with the YAML-configured values
on the compute node. If `adjust_corrections` is absent, writes defaults
(no scaling effect).

## ADJUST_R_DEFAULTS

Controls runtime R-grid refinement. Already handled by
`prepare_VADAT.py` (line 4254). Activated when `REV_RGRID=T` in VADAT.

## File preservation

`submit_cmfgen.sh` performs a manual cleanup that deletes files not in
`BASE_FILES_TO_KEEP`. All three convergence control files are included
in this list to ensure they survive directory refresh.

---

# 44. Repository Hygiene

Before beginning work:

1. inspect the repository;
2. identify active code;
3. identify stale code;
4. identify duplicated scripts;
5. identify generated files;
6. preserve useful historical work;
7. archive rather than destroy uncertain material.

Never delete:

* raw models;
* reference outputs;
* provenance;
* source data

without explicit justification.

Avoid names such as:

```text
final2
final_new
fixed
latest
test2
new_version
```

Use descriptive model/version names.

---

# 45. Required Deliverables

Where relevant, generate:

```text
BASELINE_REPORT.md
PROVENANCE_REPORT.md
SN_A_CONFIGURATION.md
SN_A_2011DH_CALIBRATION.md
SN_A_SNAPSHOT_30D_REPORT.md
SN_B_CONFIGURATION.md
INTERACTION_CONFIGURATION.md
TIMELINE_AND_SNAPSHOT_REPORT.md
VEL_LAW_12_DATA_PATH.md
NONMONOTONIC_SOLVER_DESIGN.md
ZERO_FLUX_AUDIT.md
SIGMA_AUDIT.md
GRID_AND_CONSERVATION_REPORT.md
MODEL_SPEC_OPTIMISATION.md
ATOMIC_MODEL_REPORT.md
SPECTRAL_LINE_IDENTIFICATION.md
VALIDATION_MATRIX.md
MODEL_GRID_REPORT.md
PATCH_SUMMARY.md
REPRODUCTION_RUNBOOK.md
REPOSITORY_HYGIENE_REPORT.md
```

---

# 46. Agent Behaviour

When asked to modify the modelling pipeline:

### First

Inspect:

```text
AGENTS.md
repository tree
run_model_grid.py
YAML configuration
conversion scripts
current MODEL_SPEC
VADAT
IN_ITS
SN_HYDRO_DATA
relevant CMFGEN source
```

### Second

Trace the actual execution path.

Do not assume that filenames imply functionality.

### Third

Identify the smallest change that solves the problem.

### Fourth

Implement it through the orchestration system.

### Fifth

Run tests.

### Sixth

Run a dry-run.

### Seventh

Run a controlled model.

### Eighth

Compare against the baseline.

---

# 47. Never Do This

Do not:

```text
manually submit production jobs
hard-code atomic models for one SN
blindly copy STELLA grid points into CMFGEN
use 1000 CMFGEN grid points by default
remove problematic species solely to force convergence
fix populations without physical justification
blindly increase iteration counts
silently monotonise velocity
silently replace missing flux with epsilon
overwrite raw models
modify CMFGEN source unnecessarily
reuse stale generated inputs without validation
claim convergence from the existence of a spectrum file
claim a spectrum match without line attribution
```

---

# 48. Final Acceptance Criteria

A model is scientifically acceptable only if:

```text
STELLA/input provenance verified
        AND
SN_HYDRO_DATA validated
        AND
CMFGEN grid convergence demonstrated
        AND
grid is near-minimal while satisfying convergence
        AND
MODEL_SPEC is physically justified
        AND
important ions are represented
        AND
atomic-data settings are internally consistent
        AND
CMFGEN converges
        AND
spectrum is numerically stable
        AND
important lines are physically identified
        AND
model provenance is complete
```

A successful executable run alone is insufficient.

---

# 49. Final Status

Every major task must end with exactly one:

```text
PASS — scientifically and numerically validated.
```

or:

```text
PASS WITH LIMITATIONS — usable only within documented assumptions.
```

or:

```text
FAIL — fixes required.
```

or:

```text
BLOCKED — required source, input data or compute environment unavailable.
```

Do not report `PASS` if:

* the required input version is unverified;
* `SN_HYDRO_DATA` has not been validated;
* the CMFGEN grid has not been tested for convergence;
* the atomic model is unjustified;
* important ions have been removed solely for convergence;
* velocity has not been demonstrated to propagate correctly;
* stale inputs/binaries may have been used;
* production was bypassed around `run_model_grid.py`;
* the repository is left ambiguous or contaminated with stale generated files.

```
```

---

# 29. CMFGEN Restart/Requeue Stability (2026-09-01)

## 29.1 Problem

CMFGEN jobs using the `--requeue` mechanism (wall-time graceful stop →
`scontrol requeue`) exhibited three restart-related instabilities:

1. **OUTGEN unbounded growth**: CMFGEN opens `OUTGEN` with
   `'UNKNOWN','APPEND'` mode (`cmfgen.f:99`). The first run deletes
   `OUTGEN` (`rm -f batch.log OUTGEN`), but requeues do NOT truncate it.
   Result: `OUTGEN` grows by ~150 GB per 24h requeue cycle. Observed:
   **433 GB** after 3 requeues in `X_nomono_v6/cmfgen/35.00days/`.

2. **No divergence detection**: The requeue logic only checks for wall
   time and fatal patterns (NaN, segfault). It does NOT check whether
   the model is converging. A diverging model (max % change increasing
   over iterations) requeues indefinitely, wasting compute resources.
   Observed: `SMOOTH_POP_CORRECTIONS` showed max corrections **increasing**
   from 0.21 → 0.55 over 10 iterations at iteration 2929 — clear
   divergence, but the job kept requeuing.

3. **No max requeue count enforcement**: `#SBATCH --requeue` allows
   unlimited requeues. The `FINAL_REQUEUE_MAX` parameter was defined
   but not enforced in the slurm script.

## 29.2 Root Causes

* `cmfgen.f:99` — `GEN_ASCI_OPEN(LUER,'OUTGEN','UNKNOWN','APPEND',...)`
  appends to `OUTGEN` on every CMFGEN invocation, including requeues.
* `submit_cmfgen.sh` — the generated slurm templates (SINGLE_FINAL and
  FINAL) did not rotate or truncate `OUTGEN` on requeue.
* `submit_cmfgen.sh` — no convergence tracking or divergence detection
  was performed between requeues.
* The BA (statistical equilibrium) solver falls back to LAMBDA iteration
  for complex atomic models (482 levels across 17 species in the
  X_nomono_v6 model), which is very slow and can enter non-convergent
  limit cycles that Ng acceleration cannot break.

## 29.3 Fixes Applied (in `cmfgen_base/submit_cmfgen.sh`)

### OUTGEN rotation

Added `rotate_outgen()` function to both SINGLE_FINAL and FINAL slurm
templates. On requeue, if `OUTGEN` exceeds 2 GB, the last 5000 lines
are saved to `OUTGEN.tail` and the file is truncated. CMFGEN then
appends to the empty file. This prevents the 433 GB growth while
preserving recent diagnostics.

### Convergence tracking

Added `track_convergence()` function. After each `batch.sh` run, the
last "Maximum % increase at depth" value and iteration count are
extracted from `OUTGEN` (using `tail` for efficiency on large files)
and appended to a small `CONVERGENCE_HISTORY` file that persists across
requeues.

### Divergence detection

Added `check_convergence_for_divergence()` function. Before each
requeue, `CONVERGENCE_HISTORY` is checked: if 5 consecutive entries
show max % change > 50%, the model is flagged as diverging and the job
exits with error code 4 instead of requeuing. The error message
includes diagnostic guidance (simplify MODEL_SPEC, check Ng parameters,
check BA solver fallback).

### Max requeue count enforcement

Added `SLURM_RESTART_COUNT >= MAX_REQUEUE` check. If the requeue count
exceeds `FINAL_REQUEUE_MAX` (default 500, configurable via YAML
`cmfgen.final_requeue_max`), the job exits with error code 5.

## 29.4 Files Modified

* `cmfgen_base/submit_cmfgen.sh` — added `track_convergence()`,
  `check_convergence_for_divergence()`, `rotate_outgen()` functions
  and corresponding calls in both SINGLE_FINAL and FINAL slurm
  templates.

## 29.5 Verification

* `bash -n cmfgen_base/submit_cmfgen.sh` — syntax check passes.
* Generated SINGLE_FINAL and FINAL slurm templates — syntax check
  passes with `bash -n`.
* `MAX_REQUEUE` properly expanded from `FINAL_REQUEUE_MAX` at
  generation time; `\${MAX_REQUEUE}` properly escaped for runtime
  expansion in the slurm script.
* `run_model_grid.py` copies `cmfgen_base/` to model directories via
  `shutil.copytree()`, so the fix propagates to new runs.

## 29.6 Remaining Issues (Not Restart-Related)

The underlying model convergence issue in `X_nomono_v6` is NOT a
restart bug — it is an atomic model complexity problem:

* `MODEL_SPEC` has 482 levels across 17 species (Fe2: 113 levels,
  FeI: 87 levels, OI: 48 levels).
* The BA solver falls back to LAMBDA iteration (very slow for this
  many levels).
* The model diverges (SMOOTH_POP_CORRECTIONS increasing) rather than
  converging.
* Recommended fix: simplify `MODEL_SPEC` by reducing levels for Fe
  species, or use `DO_LAM_AUTO` with adjusted acceleration parameters.

The restart fixes above ensure that diverging models are detected and
stopped, rather than requeuing indefinitely and generating 433 GB
output files.
