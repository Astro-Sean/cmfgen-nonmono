# SN2024pba Modelling Pipeline

Multi-code radiative transfer pipeline for SN2024pba and comparison supernova models.

## Directory Structure

```
SN2024pba/
├── configs/              YAML model configurations (one per model variant)
├── scripts/              Python scripts and shell wrappers
│   ├── run_model_grid.py     Master workflow orchestrator
│   ├── run_model_grid.sh     Thin wrapper for run_model_grid.py
│   ├── convert_*.py          Format conversion scripts (STELLA↔CMFGEN↔SEDONA↔SuperLite)
│   ├── plot_*.py             Plotting scripts for all codes
│   ├── prepare_VADAT.py      CMFGEN VADAT generator
│   ├── submit_cmfgen_jobs.sh CMFGEN job submission helper
│   ├── stella_job_template.slurm  STELLA SLURM template
│   └── deprecated/           Old/unused scripts
├── cmfgen_base/          CMFGEN base files (copied to each model's cmfgen/ dir)
│   ├── batch.sh              CMFGEN batch script
│   ├── submit_cmfgen.sh      CMFGEN submission orchestrator
│   ├── MODEL_SPEC            Atomic model specification
│   ├── VADAT_INIT            Default VADAT template
│   ├── RDINR, T_IN           Input templates
│   ├── tablec*.dat           Opacity tables
│   └── obs/                  Observation script helpers
├── sedona_base/          SEDONA base files (copied to each model's sedona/ dir)
│   ├── param.lua             Default SEDONA parameters
│   ├── plot_spectrum.py      SEDONA spectrum plotter
│   └── sedona_job_template.slurm  SEDONA SLURM template
├── superlite_base/       SuperLite base files (copied to each model's superlite/ dir)
│   ├── input_parameters_reference.par  Default parameter reference
│   └── superlite_job_template.slurm     SuperLite SLURM template
├── sirocco_base/         SIROCCO base files
├── atomic/               Custom atomic data (f_to_s files for CMFGEN)
├── cmfgen_source/        CMFGEN Fortran source code
├── models/               Model output directories (one per model)
│   └── X_1_fix/
│       ├── stella/           STELLA hydrodynamics output
│       ├── cmfgen/30.0days/  CMFGEN spectrum synthesis
│       ├── sedona/30.0days/  SEDONA Monte Carlo spectra
│       ├── superlite/30.0days/  SuperLite Monte Carlo spectra
│       ├── plots/            All diagnostic plots (centralized)
│       └── output/           Final output files
├── outputs/              Centralized spectrum collection
├── observational_data/   Observed spectra and photometry
├── presn_models/         Pre-supernova progenitor models
├── homo_model/           Homologous expansion test models
├── docs/                 Papers, PDFs, and non-markdown reference material
├── markdown/             All .md documentation files (consolidated)
├── audits/               All audit .md files (consolidated)
├── OLD/                  Archived/deprecated models and backups
├── SN_12C*/              SN 2011dh comparison model directories
└── AGENTS.md             Agent guidelines (do not edit without testing)
```

## Usage

```bash
# Run a single model
python scripts/run_model_grid.py configs/X1_fix.yaml

# Run with specific spectrum codes
python scripts/run_model_grid.py configs/X1_fix.yaml --run-sedona --run-superlite --no-cmfgen

# Debug mode (interactive partition, 2hr limit)
python scripts/run_model_grid.py configs/X1_fix.yaml --debug --run-sedona

# Process all configs
python scripts/run_model_grid.py --all
```

## Workflow

1. **STELLA** (or reuse existing): Radiation hydrodynamics → `.swd` file
2. **CMFGEN conversion**: `.swd` → `SN_HYDRO_DATA` (via `convert_stella_swd_to_cmfgen.py`)
3. **CMFGEN**: `SN_HYDRO_DATA` → spectrum (non-LTE atmosphere)
4. **SEDONA**: `SN_HYDRO_DATA` → `.mod` → Monte Carlo spectrum
5. **SuperLite**: `SN_HYDRO_DATA` → `input.str` → Monte Carlo spectrum

Once `SN_HYDRO_DATA` exists, CMFGEN/SEDONA/SuperLite run independently with no dependencies.

## Plot Organization

All plots go to `models/<model_name>/plots/` with filenames:
```
{code}_{epoch}_{type}.png
```
Examples: `sedona_30d_spectrum.png`, `superlite_30d_sdec_spectrum.png`
# cmfgen-nonmono
