# AGENTS.md — SN2024pba Two-Event Supernova Modelling Pipeline

## 1. Purpose

This repository is a **scientific modelling project**, not a job-running project.

It tests a two-event interacting-supernova scenario for SN 2024pba:

- **SN_A**: an earlier Type IIb explosion that produced metal-rich ejecta.
- **SN_B**: a later explosion / shell interaction that collides with the older SN_A ejecta.

The scientific question is:

> Can heating and ionisation of the metal-rich SN_A ejecta by the SN_B
> interaction reproduce the observed SN 2024pba spectrum, especially
> the strong absorption associated with H, He, O, Si, Ca and Fe?

A model answers that question only if those features can be **physically
attributed**: which ion, which transition, which radius, which velocity,
and which ejecta component (SN_A metal-rich material vs SN_B / CSM).

The repository contains a coupled STELLA → CMFGEN → SEDONA/SuperLite
modelling workflow. Those codes are tools for the question above.
Getting CMFGEN to finish, producing a plot, or lowering MAXCH is not
the scientific objective.

This document is the authoritative guidance for AI agents working on
the repository.

Historical notes, old prompts, copied code comments, chat summaries and
previous `AGENTS.md` files are **not authoritative** if they conflict
with the current source code, configuration, input files or reproducible
tests.

### 1.1 What success looks like

A scientifically useful model must do all of the following:

1. Start from a hydrodynamical structure whose provenance, time, units
   and conservation properties are known (STELLA `.swd` → `SN_HYDRO_DATA`).
2. Preserve the interaction physics that the spectrum depends on:
   shocks, velocity reversals, composition interfaces, and the
   distinction between gas temperature and radiation temperature.
3. Use an atomic model in which H, He, O, Si, Ca and Fe are present
   for physical reasons, not because they were convenient to converge.
4. Produce a numerically stable spectrum whose important absorption
   and emission features are identified with ion, transition, formation
   depth and velocity.
5. Test, where possible, whether those features arise from heated SN_A
   ejecta rather than from the SN_B shell/CSM alone.
6. Record enough provenance that another person can repeat the run from
   the YAML config through `scripts/run_model_grid.py`.

A spectrum file, a SLURM `COMPLETED` state, or a falling MAXCH history
is **not** success.

### 1.2 What this project is not

Do not treat the following as the goal:

- forcing CMFGEN (or SuperLite) to converge at any physical cost;
- deleting ions, levels or species because they are numerically awkward;
- silently monotonising a shocked or reversing velocity field;
- replacing missing flux, temperature or populations with a floor so
  that the code continues;
- matching the observed spectrum by eye without line attribution;
- assuming that STELLA continuum opacity, SuperLite LTE, and CMFGEN NLTE
  are interchangeable answers to the same question;
- using a nebular / geometrically diluted ionisation approximation in
  an optically thin, shock-heated interaction where T_gas and T_rad
  differ;
- presenting a numerical warning, or the absence of one, as a scientific
  result.

Numerical methods (Lambda iteration, DIAG BA, JFNK, Ng acceleration,
grid AMR, proxy populations) exist to **solve the physical equations
faithfully**. They are not the science.

### 1.3 Roles of the codes

Keep the codes in their physical roles:

| Code | Role | It does not |
|------|------|-------------|
| STELLA | Radiation-hydrodynamics: density, velocity, T_gas, T_rad, composition, shocks vs time | Compute an NLTE line spectrum |
| CMFGEN | NLTE statistical equilibrium and radiative transfer on a snapshot; line formation | Invent hydrodynamics or erase non-monotonic flow |
| SuperLite / SEDONA | Independent spectral / Monte Carlo cross-check | Replace CMFGEN line attribution |

YAML is **intent**. The files CMFGEN actually ran (`VADAT`, `IN_ITS`,
`SN_HYDRO_DATA`, `MODEL_SPEC`, `OUTGEN`) are **fact**. Report fact.

### 1.4 Good scientific practice — always

These rules apply to every task: diagnosis, code changes, config
changes, plots, and written conclusions. They are not optional and
they are not waived because a job is running, a deadline is close, or
a change “only affects numerics”.

**Physics before numerics.** If a model is unstable, identify whether
the cause is physical, atomic, grid, initial-population, numerical or
input-format (see §27). Do not start by increasing iteration counts,
fixing populations, or removing ions.

**Evidence, not plausibility.** Classify every technical claim as
`VERIFIED`, `SUPPORTED`, `UNVERIFIED`, `CONTRADICTED` or `BLOCKED`
(§5). A result that “looks reasonable” remains `UNVERIFIED`. Never
upgrade it to `VERIFIED` without source, a test, a run, or inspection
of generated files.

**Trace the executed path.** Read the code that ran and the files that
were written. Do not infer behaviour from a filename, a YAML comment,
a previous chat, or a default you expect CMFGEN to have used. Confirm
`VEL_LAW`, `DO_LAM_IT`, `USE_JFNK`, `LTE_EST` and similar keys from
`OUTGEN` / `VADAT`, not from the YAML alone.

**Diagnose before editing.** Inspect structure, abundances, T_gas,
T_rad, optical depth, SIGMA, populations and the actual warning. A
diagnostic that says “NaN/negative” may be a 1e-30 floor on an EUV
Wien tail; a “successful” spectrum may have no physically identified
lines. Interpret the quantity, not the English.

**Do not silently change the physical model.** No global velocity
monotonisation, no generic flux epsilon, no abundance renormalisation
that hides missing mass, no temperature floor that masquerades as
shock heating. Any numerical floor must be local, documented, tested
and physically justified.

**Conservation and units are science.** Check mass, elemental mass,
kinetic energy, abundance normalisation, radius/velocity orientation
and time definitions (`t_A` vs `t_B_SBO` vs snapshot age). A unit
error is a wrong model, not a plotting issue.

**Initial populations are guesses.** LTE, nebular (Mazzali & Lucy),
dilute-LTE and proxy `*_IN` files are starting estimates. CMFGEN’s
converged populations are the solution. Do not treat the guess as the
ionisation structure of SN 2024pba.

**Line identification is required.** A trough near an observed
wavelength is not an ID. Report species, ion, transition, formation
radius, velocity, temperature, density, abundance and optical depth
for the important H, He, O, Si, Ca and Fe features (§30–§31).

**Attribute the ejecta component.** Where possible, compare SN_A only,
SN_B / CSM only, and the interaction, and perform species-ablation
tests. The hypothesis is specifically that **heated SN_A metal-rich
ejecta** produce the observed absorption.

**Reproduce through the orchestrator.** Production, debug, restart and
re-run jobs go through `scripts/run_model_grid.py` (§2). A result
obtained by a raw `sbatch` is not a pipeline product.

**Test, then dry-run, then a controlled model.** Do not jump from an
edit to a production submission (§37–§38). After a Fortran change,
rebuild the executable that the YAML `install_dir` actually points to.

**Minimal, reversible source changes.** Prefer input preparation,
YAML, `MODEL_SPEC` and diagnostics over editing CMFGEN/STELLA/SuperLite
source. Source patches must be reviewable and must not alter monotonic
`VEL_LAW≠12` behaviour without cause (§43).

**Report limitations honestly.** Every major task ends with exactly one
of `PASS`, `PASS WITH LIMITATIONS`, `FAIL` or `BLOCKED` (§49). If the
spectrum is unattributed, the velocity path unverified, or JFNK
untested at runtime, say so.

**Leave the repository scientifically usable.** Do not overwrite raw
models, mix stale binaries with new source, pollute `cmfgen_source/`
with run products, or leave two conflicting “current” configs.

### 1.5 Observational comparison

Compare models to SN 2024pba data in `observational_data/` at the
correct epoch, in consistent flux/luminosity units, and with the same
redshift/velocity convention the observations use.

**SN_A fit window.** When judging the best-fitting SN_A lightcurve
model, only photometry with **MJD < 60544** counts. Later data belong
to the second-event / SN_B epoch and must not be held against an
SN_A-only model. The scoring tools encode this cut
(`FIT_MAX_MJD` in `scripts/plot_sn_a_pba_grid_vs_obs.py`).

**Current best SN_A fit (2026-10 rescore).** Scoring every usable
`.tt`/lightcurve under `SN_A_grid/` gives `SN_A_g_T20k_Ni075` at rms
0.19 on the SDSS (`tt_sdss`) photometric system, and `SN_A_g_T20k_Ni05`
(`archive/old_sn_a_models/SN_A/SN_A_g_T20k_Ni05/`, E=9, Ni56=0.05,
M_ext=0.05, R_ext=200, T_base=20kK) at rms 0.24 on the proxy `.tt`
system. The two systems carry a ~0.3 mag systematic offset, so compare
like with like. `SN_A_X_E07Ni05_R400M20` scores 0.60 (tt_sdss) in the
current scoring; it remains the ejecta source for the `_onR400M20`
SN_B runs.

Do not:

- shift the model in wavelength until it “looks right” without reporting
  the shift;
- scale the continuum arbitrarily and then claim a line-strength match;
- compare a 40 d post-SBO snapshot to a spectrum from a different
  reference time without stating both clocks.

---

# 2. Core Execution Rule — `run_model_grid.py` Is the Single Entry Point

## 2.1 Mandatory rule

**Every model must be prepared, validated and submitted through:**

```text
scripts/run_model_grid.py
```

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

### 2.1.1 Strict enforcement

This rule is **absolute** and applies to every submission, including:

- production runs;
- debug runs;
- single-epoch runs;
- restart runs;
- interactive-partition runs;
- tests of newly rebuilt executables;
- re-runs after source-code fixes;
- re-runs after configuration changes;
- re-runs of failed models;
- re-runs with a different grid size, temperature floor, partition, or CPU count.

There is **no** scenario in which an agent may submit a CMFGEN, STELLA, SuperLite or SEDONA job with a raw `sbatch` call and treat the result as a valid pipeline product.

If `run_model_grid.py` fails to submit, the correct response is to **fix the orchestrator or the YAML config**, not to bypass it.

If a raw `sbatch` is used as a one-off diagnostic (e.g. to test whether a script is syntactically valid SLURM), the resulting job must be clearly labelled as a diagnostic, must not be treated as a model output, and must be cancelled once the diagnostic is complete.

### 2.1.2 How to submit a model

The canonical invocation is:

```bash
python scripts/run_model_grid.py <path/to/config.yaml> [--epochs <day> ...]
```

The orchestrator:

1. resolves the YAML configuration;
2. determines which codes should run;
3. determines which existing outputs can be reused;
4. prepares model directories;
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

## 2.3 Hardcoded partition time / graceful-stop rules

SLURM time limits and graceful-stop timing are **hardcoded per partition** and must **not** be overridden by YAML config, CLI flags, or environment variables.

| Partition | Time limit | Graceful stop | USR1 signal before wall |
|-----------|------------|---------------|-------------------------|
| `interactive` | `01:59:00` | `5400s` (1 h 30 m) | `1740s` (29 m) |
| any other (e.g. `general`, `small`) | `23:59:00` | `79200s` (22 h) | `7140s` (119 m) |

These values are enforced in:

- `scripts/run_model_grid.py` — `get_time_limit()` and `get_graceful_stop_sec()`
- `cmfgen_base/submit_cmfgen.sh` — graceful stop calculation block
- `scripts/submit_cmfgen_jobs.sh` — default `TIME` positional arg

The `graceful_stop_sec` key in YAML `cmfgen:` config and the `GRACEFUL_STOP_SEC` environment variable are **ignored**.

The graceful stop fires an in-job `sleep` timer that sets `NUM_ITS=0` in `IN_ITS`, causing CMFGEN to exit cleanly after the current iteration. The `--requeue` SLURM flag then requeues the job so it continues from the checkpoint.

**Requeue is conditional on a clean finish.** A FINAL/SINGLE_FINAL job is rescheduled *only* when the graceful-stop handler set the requeue state file `.CMFGEN_REQUEUE_STATE` to `READY=1` (NUM_ITS=0) *and* `check_outgen_for_fatal` found no fatal pattern in `OUTGEN`. Any other early exit — a crash, a Fortran `STOP`, the `MAX_CHNG_LIM` guard, or simply a missing/stale `RVTJ` — fails the job with a non-zero exit instead of rescheduling, so a deterministic failure cannot loop on requeue. (Policy revised 2026-09-13; previously any missing-RVTJ exit was requeued up to 10 times. State file introduced 2026-09-15; replaces the `.CMFGEN_REQUEST_REQUEUE` sentinel with an explicit `READY=0`/`READY=1` state machine.)

The requeue state machine works as follows:

1. **Before `batch.sh`:** `write_requeue_state 0` writes `READY=0` to `.CMFGEN_REQUEUE_STATE`. This indicates CMFGEN is running and no graceful stop has been requested.
2. **On graceful stop (USR1 signal or sleep timer):** `request_graceful_stop()` sets `NUM_ITS=0` in `IN_ITS` and calls `write_requeue_state 1`, updating the state file to `READY=1`.
3. **After `batch.sh`:** The post-run logic reads `.CMFGEN_REQUEUE_STATE`:
   - `READY=1` + no fatal patterns in OUTGEN + no divergence → `scontrol requeue`
   - `READY=0` + valid RVTJ (newer than `.CMFGEN_JOB_START`) → job complete, no requeue
   - `READY=0` + missing/stale RVTJ → crash/early exit, fail with exit code 6
4. **On requeue:** The state file is deleted before `scontrol requeue`, so the next requeue cycle starts fresh.

The state file also records `TIMESTAMP` and `JOB_ID` for debugging. It is cleaned up on full resubmission (not requeue) via the cleanup block in `submit_cmfgen.sh`.

**Fresh start on resubmission (2026-09-26).** A new submission is a
cold start, not a continuation. Both wipe paths now remove warm-start
state: `POINT1`/`POINT2`/`SCRTEMP`/`SCRATCH/` (the `SCR_READ_V2`
checkpoints), species `*_IN`/`*_OUT`, `POP*`, `RVTJ`, `MODEL`,
`CUR_MODEL_DATA`, `FGINEXT_GS`/`MOMJNXT_GS`, and the
`.CMFGEN_REQUEUE_COUNT` file (a stale count >0 would otherwise route a
new submission into the requeue branch and skip VADAT prep).
Continuation survives only inside a job's own requeue cycle; deliberate
warm restarts use `submit_cmfgen.sh --restart` or
`cmfgen.warm_restart: true` in the YAML (`CMFGEN_WARM_RESTART=1` also
works). `skip_final` (OBS-only) resubmits still preserve everything
because those files are OBS inputs.

### 2.3.1 `submission:` key (partition selection)

The top-level YAML `submission:` key selects the partition (and thus the
hardcoded time limit / graceful stop) for the **whole model**:

| `submission:` | Partition | Time limit | Graceful stop |
|--------------|-----------|------------|---------------|
| `normal` | `small` | `23:59:00` | `79200s` (22 h) |
| `debug` | `interactive` | `01:59:00` | `5400s` (1 h 30 m) |

`submission` takes precedence over the legacy `debug:` boolean and over
per-section `partition:` / `time:` keys. When `submission` is set, the
per-section `partition:` / `time:` / `init_partition:` / `init_time:` keys
are stripped from the in-memory config (see `resolve_submission()` in
`scripts/run_model_grid.py`). When `submission` is absent, the legacy
behaviour (`debug:` + per-section `partition:`) is preserved.

CMFGEN requeue behaviour (`final_requeue_when: wall`, `final_requeue_max`)
applies in **both** modes; debug only disables requeue when
`final_requeue_when=never` or `final_requeue_max=0`.

### 2.3.2 STELLA concurrency cap (`MAX_RUNNING_STELLA`)

When a single `run_model_grid.py` invocation submits multiple configs,
at most `MAX_RUNNING_STELLA` (= 2, defined in `scripts/run_model_grid.py`)
STELLA hydro jobs may run concurrently. The STELLA job of config *i* is
submitted with an `afterany` dependency on the STELLA job of config
*i* − `MAX_RUNNING_STELLA` (round-robin lanes).

- `afterany` (not `afterok`) is deliberate: one failed grid point must not
  block the rest of the chain.
- Only actually-submitted STELLA jobs consume a lane — configs with
  `skip_stella`/`reuse_stella` or complete existing output do not.
- An explicit `--stella-dependency` applies to the first
  `MAX_RUNNING_STELLA` lane heads; later jobs are chained automatically.

### 2.3.3 Viper partition placement (verified 2026-09-25)

`#SBATCH --partition=general` does not mean the job runs on `general`.
The Viper job-submit filter moves a partial-node request onto `small`.
Confirm the partition from `scontrol show job`, not from the script.

Compute nodes in both partitions (`vipc`) have 256 CPUs (2 sockets x 64
cores) and 480000 MB of memory.

| Partition | MaxTime | What a job actually gets |
|-----------|---------|--------------------------|
| `interactive` | 2:00:00 | Debug jobs under 2 h |
| `small` | 1-00:00:00 | One node, shared. `MaxNodes=1`, `MaxMemPerNode=240000` MB |
| `general` | 1-00:00:00 | Exclusive whole node (`OverSubscribe=EXCLUSIVE`). The allocation is 256 CPUs |

Checked with `sbatch --test-only` on 2026-09-25:

- `--partition=general` with 8 or 16 CPUs, including `--exclusive`, was placed on `small`.
- `--partition=general` with 72 CPUs and 100 GB was placed on `general` and assigned 256 processors.
- Asking for 256 CPUs directly failed (`Requested node configuration is not available`).
- Asking for two nodes failed (`node fill grade less than 50%`).

`SLURM_CPUS_ON_NODE` on `general` is 256, the whole node, not the
`--cpus-per-task` request. `batch.sh` and the SLURM script must set
`OMP_NUM_THREADS` and `MKL_NUM_THREADS` from `SLURM_CPUS_PER_TASK`.
Using `SLURM_CPUS_ON_NODE` makes CMFGEN start one thread per core on the
node.

A 16-CPU CMFGEN job belongs on `small`. It still uses the 23:59:00 wall
and the 79200 s graceful stop from the table in §2.3. `submission: normal`
already selects `small`.


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
├── docs/                Papers, PDFs, and non-markdown reference material
├── markdown/            All .md documentation (moved from docs/ and root)
├── audits/              All audit .md files (consolidated from audit/ and root)
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

This policy applies to **scientific claims** (line IDs, ionisation,
component attribution, “the model matches the data”) as well as to
code and numerical conclusions. A plot is not evidence until the
quantity it shows has been identified.

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

## 7.1 STELLA executable: NVARS consistency fix and NCN_JUMP stability patch

The STELLA installation at `/ptmp/sbrennan/stella/OpenStella` had two critical
issues that prevented the `SN_B_v2` model from evolving past initialization:

### Issue 1: NVARS mismatch across source files

The OpenStella checkout had `PARAMETER(NVARS=4)` in `stradsep5tt.f` but
`PARAMETER(NVARS=3)` in all other source files (`begradsep.f`, `eddi.f`,
`stiffbghY12m.f`, `hcdfnrad.f`, `traneq.f`, etc.).

`NVARS` controls the number of independent variables per zone in the state
vector `Y`:
- `NVARS=3`: variables are radius, velocity, temperature (3 equations per zone)
- `NVARS=4`: would add a 4th variable (convective velocity `UC`), but `DIFMAT`
  in `hcdfnrad.f` only computes 3 equations, leaving empty Jacobian rows

The mismatch caused two distinct failures depending on which files had which
`NVARS` value:

1. **With `NVARS=3` in `begradsep.f`, `NVARS=4` in `stradsep5tt.f`** (the
   original state): the `COMMON/CUTOFF/` block had different sizes in the two
   files. `begradsep` wrote `FLOOR(4), Wacc(4), FitTau, TauTol, Rvis, CONV`
   (12 elements), but `stradsep5tt` read `FLOOR(5), Wacc(5), FitTau, TauTol,
   Rvis, CONV` (14 elements). This caused `stradsep5tt` to read `Rvis` as
   `FitTau` and `CONV` (a logical) as `TauTol`, producing `TauTol=NaN`.
   The Ncnew search condition `TAU*TauTol < FitTau` became `NaN < 1.0` which
   is always false, so the Ncnd reduction was skipped. The solver then tried
   to factorize a 500-zone matrix, causing `Y12MFF IFAIL=24`.

2. **With `NVARS=4` in all files**: the solver tried to solve for 4 variables
   per zone, but `DIFMAT` only computed 3 equations, leaving empty rows in the
   Jacobian. This caused `Y12MFF IFAIL=17` (empty row in sparse matrix).

**Fix**: Set `PARAMETER(NVARS=3)` consistently in ALL source files, including
`stradsep5tt.f`. This matches the 3-equation design of `DIFMAT` and makes the
`COMMON/CUTOFF/` block consistent across all compilation units.

### Issue 2: Ncnd decrease rate

The original code allows `Ncnd` to drop from 500 to 5 in a single step when
the optical depth is very low. This can cause solver instability at
initialization.

**Fix**: Added `PARAMETER(NCN_JUMP=20)` and the line
`IF(Ncnew>0.AND.Ncnew<Ncnd-NCN_JUMP)Ncnew=Ncnd-NCN_JUMP` after the
`IF(Ncnew<NCNMIN)Ncnew=0` check in both code paths of `stradsep5tt.f`.
This limits the Ncnd decrease to 20 zones per step, making the transition
gradual (500→480→460→...→5). When `Ncnew=0` (all-thin switch), the limit
is not applied.

### Current source parameters

| Parameter | Value | Purpose |
|-----------|-------|---------|
| `NVARS` | `3` | Variables per zone: radius, velocity, temperature |
| `NCNMIN` | `5` | Minimum diffusion zones; below this, switches to all-thin |
| `NCN_JUMP` | `20` | Maximum Ncnd decrease per step (added for stability) |
| `TREZON` | `3.D18` | Zone dropping threshold (effectively disabled) |
| `tsmoth` | `4.D6` | Smoothing time (force line commented out, as in original) |
| `Ncnew < NCNMIN` | `Ncnew = 0` | Switches to all-thin mode (ray tracing) |

### Rebuilding

All `.f` files must be recompiled with `NVARS=3`. The full rebuild procedure:

```bash
cd /ptmp/sbrennan/stella/OpenStella/strad
for f in *.f; do
    ifort -c -132 -save -zero -O3 -fp-model strict "$f"
done
cp *.o ../obj/

cd ../src
for f in *.f; do
    ifort -c -132 -save -zero -O3 -fp-model strict "$f"
done
cp *.o ../obj/

cd ../eve
ifort -c -132 -save -zero -O3 -fp-model strict evesn2.f
cp evesn2.o ../obj/

# Link xstella6new.exe (exclude evesn2.o and opazr.o to avoid symbol conflicts)
cd ../obj
rm -f evesn2.o opazr.o
ifort -o ../run/strad/xstella6new.exe *.o

# Link eve2.exe (uses opazr.o, not oparon.o)
cd ../src
ifort -c -132 -save -zero -O3 -fp-model strict opazr.f
cp opazr.o ../obj/
cd ../eve
ifort -o ../run/eve/eve2.exe evesn2.o ../obj/length.o ../obj/oparon.o \
    ../obj/stradio.o ../obj/sahaandd.o ../obj/sahaz.o ../obj/azdat.o \
    ../obj/hapsepnc.o ../obj/ubv.o ../obj/opazr.o ../obj/volenpumnoint.o \
    ../obj/vtimef90.o ../obj/words.o ../obj/burnc.o ../obj/cosetbgh.o
```

The `.trf` source requires the `trefor` preprocessor, which has environment
issues on this system. The `.f` files are the directly compilable source.

### OpenStella_conv reference

The clean reference tree at `/ptmp/sbrennan/stella/OpenStella_conv` has
`NVARS=4` in all files. This is a different code version that properly
supports 4 variables per zone. The OpenStella tree only supports 3
variables per zone and must use `NVARS=3` consistently.

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
A falling MAXCH, a finite luminosity, or a plotted flux is not a line identification.

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

### 31.1 SDEC species-ablation tool (pseudo-SDEC decomposition)

`scripts/submit_sdec_ablation.py` builds per-element ablated OBS
spectra for a converged model and a TARDIS-SDEC-style stacked
decomposition plot. It is the implemented form of the ablation
experiments above:

- creates `<cmfgen_dir>/SDEC/<Elem>/` workdirs, each a cmf_flux OBS run
  with `SCL_<species>_ABUND` scaled to ~0 (default `1.0D-10`), removing
  that element's opacity while all populations stay frozen;
- `SCL_*_ABUND` is the only ablation compatible with the non-monotonic
  OBS solver (`GLOBAL_LINE=BLANK`, all `TRANS_*=BLANK`; per-ion
  `TRANS_*=NONE` is ignored in BLANK mode);
- `../` inside each `<Elem>` dir resolves `VADAT`, `RVTJ`, `batch.sh`,
  `MODEL_SPEC`, `MODEL`, `POP*` and `FGINEXT_GS` via symlinks in
  `SDEC/` (cmf_flux reads `POP<species>` through `DIR_NAME` from the
  `../RVTJ` stdin path);
- jobs chain `afterok` with at most `--max-parallel` (default 2)
  concurrent; a final job makes the plot via
  `scripts/plot_sdec_decomposition.py`.

Orchestrator integration: `run_model_grid.py --sdec` (or YAML
`sdec.enabled: true`, optional `species`, `max_parallel`, `scale`)
submits SDEC jobs chained `afterok` on FINAL/J-RELEASE, with the plot
job also waiting on the pipeline OBS job for the baseline spectrum.

Caveats: components are not additive (line overlap/blanketing) and
populations are frozen — same approximation as TARDIS last-interaction
tagging. `SCL_HYD`/`SCL_HE` keep the e.s. continuum (`ESEC` uses `ED`
from RVTJ) but remove H/He lines and edges.

**Baseline consistency (22-Sep-2026).** The emergent flux contains an
`ESEC*J` electron-scattering emissivity evaluated on the `FGINEXT_GS`
I_NEXT field, so the baseline and every ablation must restore the same
I_NEXT generation (check the `RESTORED COMPLETE I_NEXT GENERATION
(tag=N)` line in each `OUT_FLUX`). A pipeline `obs/obs_fin_10` produced
by an OBS job that read `FGINEXT_GS` while FINAL was still writing it
carries a different generation; the resulting deltas pick up a
species-independent offset that masquerades as the same feature under
several elements (seen as fake "He and Fe emission" under H-alpha on
TypeIIn_KNOP/10.00days). `submit_sdec_ablation.py` therefore always
submits an unablated `ALL` job in the same batch (opt out with
`--no-submit-all`), and `plot_sdec_decomposition.py` prefers
`SDEC/ALL/obs_fin_10` as the baseline and warns on tag mismatches.
The plot's default `exclusive` attribution assigns each wavelength bin
to the species with the largest delta so each component is owned by a
single element; `--attribution shared` reproduces the old stacked
overlap.

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

## 39.1 `VEL_LAW=12` is downgraded to 11 under `NON_MONO_SOLVER=T`

**VERIFIED (2026-09-08):** When `NON_MONO_SOLVER=T`, `prepare_VADAT.py`
rewrites `VEL_LAW` from 12 to 11 and logs:

```text
[NON_MONO] VEL_LAW forced: 12 -> 11 (required for non-monotonic static BLANK solver)
```

This is required, not a bug: the static non-monotonic BLANK solver
`FG_J_CMF_NONMONO_V1` is reached through `SET_RV_HYDRO_MODEL_V3`, which
is the `VEL_LAW=11` code path. `VEL_LAW=12` (`RV_SN_MODEL_SNIIN`) does
not route to the non-monotonic formal solver.

Consequences that must be accounted for:

* The `vel_law` key in YAML is retained at 12 because the STELLA→CMFGEN
  converter keys its `SIGMA` formula off that value and expects the
  `VEL_LAW=12` finite-difference form (`sigma_velocity`). Under
  `VEL_LAW=11` the solver instead evaluates the analytic derivative of a
  monotone-cubic (PCHIP) interpolant, so the written `Sigma` column and
  the solver's internal `SIGMA` are computed by **different formulae**.
* The section 39 checklist above must therefore be verified against the
  `VEL_LAW=11` path, not the 12 path.
* Confirm the actual value from the `[NON_MONO] VEL_LAW forced` line in
  `OUTGEN`. Never assume the YAML value is what CMFGEN ran with.

**Enforcement pipeline (VERIFIED):** `run_model_grid.py` resolves
`NON_MONO_SOLVER=T` in `enforce_cmfgen_sn_iin_settings()`, which applies
the same forced-key table as `prepare_VADAT.NON_MONO_FORCED_VADAT`
(imported, single source of truth), fills the verified non-monotonic
iteration defaults (`DO_NG=T`, `NUM_LAM=2`, `LAM_VAL=2.0D2`,
`BEG_NG=1.0D1`, `IBEG_NG=20`, `COMP_BA=T`, `FIX_T_AUTO=F`,
`FIX_NE=F`, `FIX_IMP=F`, `MAX_LIN=5.0D0`, `FIX_BA=1.0D1`,
`N_FIX_BA=2`, `BA_CHK_FAC=1.0D-05`, `EPS_TERM=5.0D-02`,
`DO_DIS=T` (level dissolution), `INC_PEN=T` (Penning ionization),
`INC_RAY=T` (Rayleigh scattering) - adopted 2026-09-25), drops any
`iteration_strategy`, and disables `monotonize_velocity` /
`enforce_homologous` (both erase the sigma<-1 structure the solver
exists to model). The resolved config is written to
`<model_dir>/config_resolved.yaml` and *that file* — not the raw YAML —
is passed to `prepare_VADAT.py --vadat-overrides` on the compute node.
`prepare_VADAT.py` re-applies the forced keys at VADAT-write time, so
the two layers cannot disagree.

**FIX_NE=F is the project default for ALL models** (adopted
2026-09-25): `enforce_cmfgen_sn_iin_settings()` stamps `FIX_NE=F`
into `vadat_overrides` on every run unless the YAML sets it
explicitly, and `prepare_VADAT.py` writes `FIX_NE=F` +
`FIX_INB_NE=F` for every profile when no override is present. The
hydro electron density is a Saha/LTE guess, not the NLTE ionization
balance; CMFGEN must solve Ne from charge conservation. `FIX_NE=T`
remains available as a deliberate override
(`vadat_overrides: {FIX_NE: T}` or `prepare_VADAT.py --fix_ne`) —
e.g. to isolate Ne as a divergence driver — but it must not silently
be the default.

**Numerics and spectral gridding** (aligned to the reference VADAT
in `cmfgen_validation/jnk_extracted/VADAT`, adopted 2026-10): solver
numerics - `METHOD=ZERO` (zero-order depth interpolation),
`MAX_CHNG=1.0D150` (abort threshold on BA corrections), `N_PAR=2000`
(frequencies solved between BA matrix rebuilds), `COMP_F=F` (skip the
initial Eddington-factor pass), `STORE_BA_INV=F`/`WR_PRT_INV=F` (do
not store BA inverses), `ALL_FREQ=F` (apply the enhanced spatial grid
only where needed), `T_MIN=0.2` (temperature floor), `T_INIT_TAU=1.5`/
`GREY_TAU=0.5` (initial-model and grey-atmosphere optical depths);
the continuum frequency grid (`MIN_CF=3.0D-3` to `MAX_CF=2.0D2` in
1e15 Hz, `FRAC_SP=1.10`/`AMP_FAC=1.05` limit the ratio of adjacent
frequencies, `MAX_BF=0.10` caps spacing near bound-free edges); line
broadening (`TDOP=1.0`/`AMASS_DOP=4.0D+10` set the Doppler width,
`MAX_DOP=6`/`FRAC_DOP=1.0` cap the profile extent,
`dV_CMF_PROF=200`/`dV_CMF_WING=300` CMF-frame velocity sampling,
`R_CMF_WING_EXT=3` wing extent in Doppler widths); line cuts
(`GF_CUT=2.0D-3` drops weak lines, `AT_CUT=10`/`MIN_TRANS=9` trim
transitions per ion); and the OBS grid (`OBS_EXT_RAT=1.1` extends the
observed-frame range, `dV_OBS_PROF/WING/BIG=100/400/2000` sampling,
`F_LAM_BEG=0`). `F_LAM_END=2.5D4` (2.5 um, user directive).
Deliberately retained: `ES_WING_EXT=12000` (ES wing extent beyond
VINF), `SIG_GAU_KMS/V_CROSS=500` (feature resolution), and
`EPS_TERM` overrides in production configs. These defaults live in
`prepare_VADAT.py` profile dicts and both `VADAT_INIT` templates;
explicit YAML `vadat_overrides` still win.

## 39.2 `SIGMA` must be recomputed after any grid modification

**VERIFIED (2026-09-08):** `SIGMA = dlnV/dlnR - 1` is a derivative of the
grid, not a per-point sample. Any step that inserts or removes radius
points invalidates it, and subsetting it with a boolean mask is wrong —
it must be **recomputed** from the final `R`, `V`.

The `max_points` hard cap in `convert_stella_swd_to_cmfgen.py` previously
subset `r`, `v`, `T`, `rho`, `n_e`, `tau`, `kappa_arr`, `mass_fractions`
and `ni56_for_isotopes` but left `sigma` at its pre-trim length. For
`SN_B_v2_M3_V5000` at 40 d this wrote a `Sigma` block of **206 values
into a file declaring 100 data points**, whose first 100 entries belonged
to different radii — they disagreed with the file's own `R`,`V` columns at
65 of 100 depths (max discrepancy 21.94).

This parsed without error because `rd_sn_data.f:165` uses
`READ(LU,*)SIGMA_HYDRO`, a list-directed read of exactly `NX` values; the
surplus lines were then skipped by the section-header loop as
unrecognised. The corrupt column fed `VOL_EXP_FAC`
(`rd_sn_data.f:246`). `SIGMA`/`SIGMA_TRUE` were unaffected because those
are derived from `input_hydro.dat`.

`write_sn_hydro_data` now validates that every profile array and mass
fraction has exactly one value per radius point and raises otherwise.
Do not remove that guard.

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

**Attribution comment required.** Every edit to CMFGEN source code must
include a comment at the site of the change stating:

```text
edited by Sean Brennan (DATE)
```

where `DATE` is the date the edit was made. This applies to every modified
line or block so that local changes can be distinguished from the original
source during review and diffs.

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
which performs the standard Lambda-to-BA handoff:

1. Lambda iterations run until MAXCH < 50%
2. `DO_LAM_AUTO` auto-switches `DO_LAM_IT=F` in IN_ITS
3. Full BA iterations begin with `LAM_VAL=200` / `NUM_LAM=2` fallback
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
  LAM_VAL: 2.0D2      # Fall back to Lambda when MAXCH exceeds this
  BEG_NG: 1.0D1       # Manual: "5% to 10%" (was 1.0D15 = disabled)
  IBEG_NG: '20'       # Minimum Lambda iterations before Ng
  BW_NG: '40'         # Depths per Ng acceleration band
  FIX_IMP: F          # Impurity populations iterated, not held fixed
  MAX_LIN: 5.0D0      # Cap on the BA correction multiplier
  FIX_BA: 1.0D1       # Keep computing BA until MAXCH < 10%
  N_FIX_BA: '2'       # Hold BA at most 2 iterations (manual recommended)
  BA_CHK_FAC: 1.0D-05 # BA accuracy-check threshold
  EPS_TERM: 5.0D-02   # Iteration convergence threshold
  DO_DIS: T           # Level dissolution (occupation probabilities)
  INC_PEN: T          # Include Penning ionization
  INC_RAY: T          # Include Rayleigh scattering
  NON_MONO_SOLVER: T  # Non-monotonic formal solver
  USE_LAM_ES: T       # J from formal solution, not moment equations
```

**Do not set `iteration_strategy`** — let `DO_LAM_AUTO` control the mode.

#### `cmfgen.use_lam_es` — verified J-update setups (KNOP mode)

Under `non_mono_mode: KNOP`, `cmfgen.use_lam_es` selects the J source
and the verified solver pairing (TypeIIn_10d runs, 22-23 Sep 2026):

```yaml
cmfgen:
  non_mono_mode: KNOP
  use_lam_es: true    # J from FG_J_CMF_NONMONO_V1 formal solution
                      # -> USE_LAM_ES=T + SOL_METH=DIAG (NUM_BNDS=1)
  # or
  use_lam_es: false   # J from MOM_J_CMF_NONMONO_V1 moment solve
                      # -> USE_LAM_ES=F + SOL_METH=TRI (NUM_BNDS=3)
```

`use_lam_es` is the setup selector: it overrides a conflicting
`depth_independent_solver` or `vadat_overrides.USE_LAM_ES` (loudly).
Under `non_mono_mode: SOB` only `use_lam_es: true` is valid (F would
fall back to the clamped-sigma `MOM_J_CMF_V11`); F requests are forced
back to T+DIAG with a warning.  Outside non-monotonic mode the key
stamps `USE_LAM_ES` only — `SOL_METH` stays owned by
`depth_independent_solver`.  OBS (`cmf_flux`) always uses
`USE_LAM_ES=T` regardless of the main run: the OBS binary has no
`MOM_J_CMF_NONMONO_V1`, so the OBS J must come from the formal
solution either way.

### Requeue configuration

Production runs use a 24 h wall and requeue. On Viper the partition for
a partial-node CMFGEN job is `small`. A `partition: general` request
with fewer than a whole node is moved to `small` (§2.3.3). A job that
stays on `general` occupies all 256 CPUs on the node.

```yaml
cmfgen:
  partition: small
  time: '23:59:00'
  final_requeue_when: wall
  final_requeue_max: 500
```

The `interactive` partition only accepts jobs under 2 h and is suitable
only for debug runs (hardcoded to 01:59:00). Setting
`final_requeue_max: 0` causes the job to abort on the first requeue,
which discards the iterations already completed.

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

## 45.1 Writing style for reports and documents

All prose written for this repository must read like it was written by
a scientist, not generated by a model. This applies to Markdown
reports, LaTeX documents, commit messages, comments intended for
publication, and any text a collaborator will read.

- Do not use em dashes. Use full stops, commas, colons or parentheses
  instead.
- Write plainly and directly. State what was done, what was found and
  what the evidence is. Do not pad sentences with hedging or
  transitional filler.
- Avoid stock AI phrasing: "delve", "crucial", "leverage",
  "moreover", "furthermore", "it is important to note", "in
  conclusion", "not only ... but also", and similar constructions.
- Prefer concrete nouns and verbs over abstract framing. "BA
  corrections diverge at depth 57" over "the iteration exhibits
  challenging convergence behaviour".
- Do not add rhetorical emphasis, bullet-point symmetry for its own
  sake, or summary paragraphs that restate what the table already
  says.
- Keep the tone flat and technical. Enthusiasm, praise and apologetic
  language do not belong in a report.
- Do not attribute authorship to the agent. Commit messages, code
  comments, reports and documentation must not mention Devin or any
  AI assistant, and must not contain "Generated with" trailers,
  `Co-Authored-By` lines, or any other tool attribution. Commits are
  authored by the user; write them as if the user wrote them.

Existing documents do not need to be rewritten, but new or edited
sections must follow this style.

## 45.2 Plot and figure style

These rules apply to every figure produced for this repository,
whether diagnostic or for presentation.

**No titles on plots.** Do not put `ax.set_title(...)` on panels and
do not add `fig.suptitle(...)` by default. The figure must be
self-describing through axis labels, units and legends. If a figure
genuinely needs identifying context (model name, epoch), put it in
the filename and, at most, one line of suptitle text at the top.
Never repeat the same information in a title, an axis label and the
filename.

**Discrete data are drawn as steps.** Quantities defined per grid
zone or per CMFGEN depth (density, temperature, velocity, sigma,
optical depth, mass fractions, population ratios) are plotted with
`ax.step(..., where='mid')`, not `ax.plot`. Connecting lines imply a
continuity the data do not have. When replacing `ax.loglog` or
`ax.semilogx` with `ax.step`, set the scales explicitly with
`ax.set_xscale('log')` / `ax.set_yscale('log')`. Bars, histograms,
scatter markers and reference lines (`axhline`, `axvline`, `axvspan`)
stay as they are.

**Depth index on the top axis.** Radius-axis plots of CMFGEN-depth
quantities carry a secondary upper x axis showing the depth index
(1 = outermost), produced by `_add_depth_axis` in
`scripts/convert_stella_swd_to_cmfgen_plots.py`. Ticks are thinned
automatically where the grid compresses so labels do not overlap.

**Units and labels are mandatory.** Every axis states its quantity
and unit. Legends appear once per distinct meaning; do not create
duplicate entries for the same thing (for example a shaded region and
a proxy patch both labelled identically).

**No annotation clutter.** Do not add text arrows, boxed labels or
insets that restate what the axes already show. An inset or arrow is
acceptable only when it reveals structure invisible at the main scale.

**Compare on a common grid.** Before differencing or plotting a
ratio between hydro input and CMFGEN output, interpolate onto the
same radius grid. Index-by-index comparison of resampled grids
manufactures artificial features. Guard interpolation against
non-monotonic grids, out-of-range points, NaN/Inf and zero
denominators.

**No redundant figures.** Before adding a new plot, check whether an
existing diagnostic already shows the quantity. When a plot is
superseded, remove both the file and the code that generates it;
do not leave dead plotting code or stale PNGs in run directories.

**Technical defaults.** Use the matplotlib `Agg` backend (no display
is available). Save to the run's `plots/` directory. Figures are
diagnostics: they must reflect the processed data honestly and are
not evidence of convergence or correctness by themselves.

---

# 46. Agent Behaviour

§1 is the goal of every task. Pipeline edits, solver flags, grid
changes and job submissions are means. If a proposed change would
make the code finish while destroying shocks, composition interfaces,
H/He/O/Si/Ca/Fe line formation, or the SN_A vs SN_B distinction,
do not make it.

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
treat a warning string as proof of NaN, failure, or success
claim VERIFIED because a result looks physically reasonable
equate falling MAXCH with a correct NLTE solution
confuse T_gas, T_rad, T_exc and ionisation temperature
use nebular geometric dilution for this optically thin interaction
present YAML intent as the settings CMFGEN actually ran
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
        AND
the SN_A / SN_B origin of those lines has been addressed
```

A successful executable run alone is insufficient. Matching the
hypothesis in §1 is the acceptance test; the list above is how that
hypothesis is made falsifiable.

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
* the important H, He, O, Si, Ca or Fe features have not been physically identified;
* SN_A vs SN_B origin of those features has not been addressed, even qualitatively;
* a numerical method (Lambda, BA, JFNK) was substituted for a physical diagnosis.

Historical notes after this section are **not** a substitute for §§1–49.

---

# Appendix A — CMFGEN Restart/Requeue Stability (2026-09-01)

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

---

# Appendix B — SuperLite Nebular Approximation: Critical Ionisation Bug (2026-09-04)

### Problem

SuperLite models of SN_B_v2 at 40 days showed **no emission or absorption
lines** — only a smooth continuum with ~5% Monte Carlo noise.

### Root Cause

The nebular approximation (`in_nebular = t`, Mazzali & Lucy 1993) was
killing all ionisation in the outer (metal-rich) cells:

1. `in_R_phot = 2.71e15 cm` (from the converter's grey tau=2/3 estimate).
2. For the metal-rich SN_A ejecta at r ~ 2.3e16 cm, the geometric dilution
   factor is W ≈ 0.003.
3. The effective radiation temperature for ionisation is
   T_eff = T_rad × W^(1/4) ≈ 4700 × 0.24 ≈ **1100 K**.
4. At 1100 K, every element is completely neutral and in its ground state.
   No excited states → no optical line opacity.

This caused SuperLite's opacity to be **15 orders of magnitude lower**
than STELLA's:

| Quantity | STELLA .swd | SuperLite (nebular) |
|----------|------------|---------------------|
| kappa_Ross (hot metal-rich zones) | 0.45 cm²/g | 3×10⁻¹⁶ cm²/g |
| Total Rosseland tau | 0.25 | ~0 |

With zero optical depth, no photosphere and no lines can form. The
spectrum is just a sum of smooth Planck functions at the local gas
temperatures.

### Why the Nebular Approximation Is Wrong Here

The nebular approximation is designed for **optically thick** ejecta
where a well-defined photosphere exists and the radiation field above
it is geometrically diluted. In this model:

- The **total Rosseland tau is only 0.25** — the ejecta is optically thin.
- The gas temperature (set by the SN_B shock interaction, T_gas ~ 8000 K)
  is **much higher** than the radiation temperature (T_rad ~ 4700 K).
- The dilution factor drives the effective ionisation temperature to
  ~1100 K, far below the gas temperature.

LTE ionisation at T_gas = 8000 K gives significant Ca II, Fe II, and
Si II fractions — exactly the ions that produce the strong lines this
model is supposed to show.

### Fix

1. **Turn off the nebular approximation**: Set `in_nebular = f` in
   `input.par`. This uses LTE Saha ionisation at the gas temperature,
   which is much higher than the diluted radiation temperature.

2. **YAML configuration**: Add `nebular: false` under
   `superlite.numerical` in the YAML config:
   ```yaml
   superlite:
     numerical:
       nebular: false
       nlte: true
       nlte_nelem: 1
       niter: 3
   ```

3. **Keep NLTE on** (`in_nlte = t`): The NLTE iteration updates the
   radiation temperature from the transport solution, which is still
   useful. But the ionisation should use LTE (gas temperature) as the
   starting point.

### Argparse Bug Fix

The converter had `type=bool` for `--sl-nebular` and `--sl-nlte`
arguments. In Python, `bool("false")` returns `True` because any
non-empty string is truthy. This meant `--sl-nebular false` would
still enable the nebular approximation.

Fixed by adding a `str2bool()` function that properly interprets
`true/false/t/f/yes/no/1/0` strings.

### Template Override Fix

The `write_input_par` function uses a reference template
(`superlite_base/input_parameters_reference.par`) that has
`in_nebular = t` hardcoded. The template path only replaced
`{MODEL_NAME}`, `{EPOCH}`, etc., but not the nebular/NLTE settings.

Fixed by adding `re.sub()` calls to override `in_nebular`, `in_nlte`,
`in_nlte_nelem`, and `in_niter` in the template output with the actual
parameter values.

### Key Files

- `scripts/convert_stella_swd_to_superlite.py`: Converter with the fix
- `SN_B_v2/config_superlite_40d.yaml`: YAML config with `nebular: false`
- `/ptmp/sbrennan/superlite/ionsmod.f`: SuperLite nebular EOS (uses
  `in_R_phot` for dilution factor)
- `/ptmp/sbrennan/superlite/GAS/eos_update.f`: EOS dispatch (nebular vs LTE)

## Preparation consistency checks (2026-09-14)

VERIFIED: VEL_LAW=11 uses SN_HYDRO_DATA depth-dependent composition. Proxy generation must use that composition and preserve the source electron-density reference (`ne_max_iter=1`); iterating a different proxy reference is rejected. Regenerate inputs through the orchestrator after preparation changes.

VERIFIED: The converter writer may change T/kappa after grid selection. Boundary metadata and the authoritative post-write audit must use serialized profiles, not pre-write arrays. Conversion stdout/stderr and previous hydro inputs are retained in per-attempt directories; a nonzero converter exit is fatal even if an output file exists.

Verification: `python3 -m pytest scripts/test_cmfgen_preparation_consistency.py scripts/test_proxy_populations_v2.py -q`; `python3 -m pytest scripts/grid_mapping/tests -q`; `bash -n scripts/submit_cmfgen_jobs.sh`; `bash -n cmfgen_base/submit_cmfgen.sh`. These tests validate interfaces, not spectral convergence or the physical accuracy of proxy estimates.

UNVERIFIED: Source/executable equivalence for the grey-initialization patch. On 2026-09-14 Raven and Viper V4 executables had identical hashes but differing rd_control_variables.f sources. Do not infer executed grey initialization from executable strings; LIN_INT=T selects GRID, which bypasses the initial grey-temperature loop.

## STELLA executable policy — parallel OpenStella mandatory (2026-09-16)

All STELLA hydro jobs must use the OpenMP-parallel OpenStella build:

```text
/ptmp/sbrennan/stella/OpenStella_parallel
```

Enforced centrally in `scripts/run_model_grid.py`:

- `DEFAULT_STELLA_ROOT` is the parallel build. A config whose
  `stella.stella_root` points at a legacy serial build (any path whose
  basename is `OpenStella`, on any mount: `/ptmp`, `/raven`, `/viper`)
  is overridden to the parallel build with a printed warning — the
  override is reported in the submission log, so the executed build is
  always visible.
- `stella.cpus_per_task` defaults to 3 when the parallel build is used
  (matching the validated `parallel_test_p3` run), 1 otherwise. An
  explicit `cpus_per_task` in the config always wins.
- Escape hatch for deliberate serial runs (e.g. the serial-vs-parallel
  baseline): set `stella.allow_serial_stella: true` in the config.
  `configs/parallel_test/parallel_test_baseline.yaml` uses this.

The parallel build was validated by `parallel_test_p1/p2/p3`
(configs in `configs/parallel_test/`); the serial baseline run is
`parallel_test_baseline`.
