# Supernova Spectrum Pipeline

Automated workflow for generating supernova spectra using **STELLA** (or HERACLES) for hydrodynamics, **CMFGEN** for detailed radiative transfer, and **SEDONA** for Monte Carlo spectrum synthesis.

## What `run_model_grid.py` Does

`run_model_grid.py` is the master workflow script that automates the entire pipeline from progenitor structure to final spectrum. It reads a YAML configuration file and orchestrates SLURM job submissions across the full simulation chain.

### Pipeline Overview

```
YAML Config → Profile Generation → Hydro Simulation → CMFGEN Conversion → Spectrum Generation
```

### Step-by-Step Workflow

#### 1. **Profile Generation** (`generate_profile`)
Creates the initial 1D density/velocity/composition structure for the supernova progenitor.
- Supports **RSG shell + beta-law wind** profiles
- Generates `input_profile` for STELLA or HERACLES
- Creates diagnostic plots (density, velocity, composition, optical depth)

#### 2. **Hydro Simulation** (STELLA or HERACLES)
Runs the hydrodynamic radiation simulation to model the explosion and ejecta expansion.

**STELLA path** (`submit_stella_job`):
- Builds STELLA executable if needed
- Runs `eve2` to generate `.mod` file from profile
- Executes STELLA with energy injection (thermal bomb)
- Produces `.swd` snapshots and `.tt` lightcurve
- Automatically converts `.swd` → `SN_HYDRO_DATA` (CMFGEN format)
- Prepares SEDONA directory with `.mod` and `param.lua`

**HERACLES path** (`submit_heracles_job`):
- Generates Fortran namelist from config
- Submits HERACLES hydro+rad job
- Converts output → `SN_HYDRO_DATA`

#### 3. **CMFGEN Directory Setup** (`prepare_cmfgen_directory`)
Copies base CMFGEN files (VADAT, batch scripts, atomic data references) into the model's `cmfgen/` directory.

#### 4. **CMFGEN Spectrum Generation** (`submit_cmfgen_jobs`)
Submits CMFGEN jobs for detailed non-LTE radiative transfer:
- **INIT**: Initial temperature structure (optional)
- **FINAL**: Full NLTE calculation with iteration
- **OBS**: Synthetic spectrum extraction
- Supports **single-final mode** (cold start, skip INIT)
- Dependencies ensure FINAL waits for hydro, OBS waits for FINAL

#### 5. **SEDONA Spectrum Generation** (`submit_sedona_job` / `prepare_sedona_directory`)
Submits SEDONA for Monte Carlo spectrum synthesis:
- Converts `SN_HYDRO_DATA` → SEDONA `.mod` format
- Generates `param.lua` from template with inner boundary conditions from hydro data
- Runs `sedona6.ex` (main simulation, not post-processing)
- Auto-plots spectrum with `plot_spectrum.py`

### Key Features

| Feature | Description |
|---------|-------------|
| **Debug mode** (`--debug`) | Uses `interactive` partition, 1:59:00 time limit, reduced simulation days |
| **Job cancellation** | Auto-cancels existing jobs for same model to prevent conflicts |
| **Dependencies** | CMFGEN/OBS wait for hydro; SEDONA waits for hydro; proper SLURM dependency chains |
| **Skip options** | `--skip-stella` to use existing output; `--no-cmfgen`/`--no-sedona` to skip codes |
| **Config inheritance** | Debug mode creates `config_debug.yaml` with modified paths and days |

### Command-Line Options

```bash
python scripts/run_model_grid.py configs/X1.yaml [options]
```

| Option | Description |
|--------|-------------|
| `--debug` | Debug mode: interactive partition, short time limit |
| `--debug-days N` | Simulation days in debug mode (default: 1) |
| `--skip-stella` | Skip STELLA, use existing `stella/` output |
| `--no-cmfgen` | Skip CMFGEN spectrum generation |
| `--no-sedona` | Skip SEDONA spectrum generation |
| `--no-redo-heracles` | Skip hydro if `SN_HYDRO_DATA` exists |
| `--all` | Process all YAML configs in `configs/` |

### Directory Structure

```
SN2024pba/
├── configs/              # Model configuration YAMLs
│   └── X1.yaml
├── scripts/              # Python workflow scripts
│   ├── run_model_grid.py          # Master workflow
│   ├── convert_stella_swd_to_cmfgen.py
│   ├── convert_cmfgen_to_sedona_v2.py
│   ├── stella_job_template.slurm
│   └── ...
├── base/                 # CMFGEN base files (VADAT, batch.sh, etc.)
├── sedona_base/          # SEDONA templates (param.lua, job template)
├── models/
│   └── X1-test-debug/
│       ├── stella/       # STELLA output (.swd, .tt)
│       ├── cmfgen/       # CMFGEN input (SN_HYDRO_DATA)
│       ├── sedona/       # SEDONA input (.mod, param.lua)
│       ├── plots/        # Diagnostic plots
│       └── output/       # Final spectra, lightcurves
└── README.md
```

### Configuration File (YAML)

The YAML config defines:
- **Ejecta parameters**: masses, radii, compositions, wind properties
- **Hydro backend**: STELLA or HERACLES with specific settings
- **CMFGEN options**: single-final mode, iteration counts, resources
- **SEDONA settings**: partition, CPUs, memory, time
- **Conversion parameters**: tau cutoffs, grid refinement, adaptive grid

### External Dependencies

| Code | Purpose | Location |
|------|---------|----------|
| **STELLA** | Hydrodynamic radiation code | `/ptmp/sbrennan/stella/OpenStella` |
| **HERACLES** | Alternative hydro code | `/ptmp/sbrennan/heracles` |
| **CMFGEN** | NLTE radiative transfer | `/u/sbrennan/cur_cmf` |
| **SEDONA** | Monte Carlo spectrum synthesis | `/ptmp/sbrennan/sedona` |

### Workflow Example

```bash
# Full pipeline: STELLA → CMFGEN → SEDONA
cd /ptmp/sbrennan/SN2024pba
python scripts/run_model_grid.py configs/X1.yaml

# Debug mode: short run for testing
python scripts/run_model_grid.py configs/X1.yaml --debug --debug-days=20

# Skip hydro, just run spectrum codes (requires existing SN_HYDRO_DATA)
python scripts/run_model_grid.py configs/X1.yaml --skip-stella

# SEDONA only
python scripts/run_model_grid.py configs/X1.yaml --skip-stella --no-cmfgen
```

### Monitoring Jobs

```bash
squeue -u $USER                    # View all jobs
watch -n 30 squeue -u $USER       # Auto-refresh every 30s
```

### Output Files

| File | Description |
|------|-------------|
| `models/*/cmfgen/SN_HYDRO_DATA` | Hydro data in CMFGEN format |
| `models/*/sedona/*.mod` | SEDONA model file |
| `models/*/sedona/spectrum.h5` | SEDONA output spectrum |
| `models/*/cmfgen/OUTGEN` | CMFGEN output spectrum |
| `models/*/plots/*.png` | Diagnostic plots |
| `models/*/output/lightcurve.csv` | Bolometric lightcurve |

---

## Repository Contents

- `scripts/` — Workflow automation, conversion tools, templates
- `configs/` — Model parameter files
- `base/` — CMFGEN base configuration files
- `sedona_base/` — SEDONA parameter templates and plotting tools
