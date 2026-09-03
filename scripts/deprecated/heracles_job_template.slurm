#!/bin/bash
# [DEPRECATED] HERACLES job template
# DEPRECATED (2026-07-09): HERACLES does not fit the science case for this project
# (interacting Type IIn supernovae). The project now uses the STELLA -> CMFGEN
# pipeline exclusively. For the active job template, see stella_job_template.slurm.
#SBATCH --job-name=heracles_{MODEL_NAME}
#SBATCH --output={MODEL_DIR}/heracles_{MODEL_NAME}_%j.out
#SBATCH --error={MODEL_DIR}/heracles_{MODEL_NAME}_%j.out
#SBATCH --partition={PARTITION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={CPUS}
#SBATCH --mem={MEMORY}
#SBATCH --time={TIME}

# HERACLES single-stage simulation for model: {MODEL_NAME}
# Mode: {MODE}
# Simulation time: {SIMULATION_DAYS} days

set -euo pipefail

# Use astro conda environment with numpy
source /mpcdf/soft/RHEL_9/packages/x86_64/python-waterboa/2024.06/bin/activate /ptmp/sbrennan/conda-envs/astro
unset PYTHONPATH PYTHONHOME

MODEL_NAME="{MODEL_NAME}"
MODEL_DIR="{MODEL_DIR}"
HERACLES_DIR="{HERACLES_DIR}"
HERACLES_ROOT="{HERACLES_ROOT}"

cd "{HERACLES_DIR}"

echo "======================================================"
echo "HERACLES Simulation: {MODEL_NAME}"
echo "======================================================"
echo "JobID: ${SLURM_JOB_ID}"
echo "HERACLES directory: {HERACLES_DIR}"
echo "HERACLES root: {HERACLES_ROOT}"
echo "Mode: {MODE}"
echo "Simulation time: {SIMULATION_DAYS} days"
echo "Time: $(date)"
echo "======================================================"

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

module purge
module load gcc/13 openmpi/4.1 mkl/2025.3
module load python-waterboa/2025.06

# Explicitly set MKL library path (workaround for module loading issues on compute nodes)
export LD_LIBRARY_PATH=/viper/u1/system/soft/RHEL_9/packages/x86_64/intel_oneapi/2025.3/mkl/2025.3/lib/intel64:$LD_LIBRARY_PATH

# Set up HERACLES paths
PIPE_PY="{HERACLES_ROOT}/scripts/supernova_model"
# Select radiation-enabled or pure-hydro executable based on mode
if [[ "{MODE}" == "hydrorad" ]]; then
  HERACLES_BIN="{HERACLES_ROOT}/heracles-2.0/bin/heracles_ray"
else
  HERACLES_BIN="{HERACLES_ROOT}/heracles-2.0/bin/heracles_noray"
fi

source "{HERACLES_ROOT}/scripts/heracles_srun_wrap.sh"
export HERACLES_LOG="raw"

# Preprocess profile if needed (for hybrid mesh)
# Disabled: input_profile is already generated with the correct n_cells
if [[ "false" == "true" ]]; then
  echo "Preprocessing profile for hybrid mesh..."
  mkdir -p modified_inputs
  python3 "${PIPE_PY}/preprocess_script.py" --apply \
    --liste liste.nml \
    --lbox-cm {LBOX_CM} --nx {N_CELLS} \
    --r-transition-cm {R_TRANSITION_CM} --nx-inner-linear {NX_INNER_LINEAR} \
    --hybrid-sigmoid-k {HYBRID_SIGMOID_K} \
    --r-inner-face-cm {R_INNER_CM} \
    --force-r-inner-face-override \
    --plot-dir "../plots/preprocess" \
    --profile input_profile \
    --profile-output input_profile \
    --profile-file-in-nml "input_profile"
fi

# Clean existing snapshots (fresh start)
echo "Cleaning existing snapshots..."
shopt -s nullglob
for d in [0-9][0-9][0-9]; do
  rm -rf "${d}"
done
shopt -u nullglob

# Run HERACLES
echo "Running HERACLES {MODE} simulation..."
heracles_run logged singlestage -- srun --kill-on-bad-exit=1 "${HERACLES_BIN}"

echo "HERACLES finished at $(date)"

# Post-process plots
echo "Generating post-process plots..."
env PYTHONNOUSERSITE=1 python3 "${PIPE_PY}/plot_hydro.py" || true

# Convert HERACLES snapshot to CMFGEN format
echo "Converting HERACLES snapshot to CMFGEN format..."
export HERACLES_ROOT="{HERACLES_ROOT}"
CONVERSION_SCRIPT="{PROJECT_DIR}/scripts/convert_heracles_to_cmfgen_wrapper.py"
# Change to model directory to resolve relative paths in config correctly
cd {MODEL_DIR}
# Use the config file path passed at job-submission time (avoids hardcoded "model_" prefix)
python3 "${CONVERSION_SCRIPT}" {CONFIG_FILE}

echo ""
echo "======================================================"
echo "HERACLES simulation complete"
echo "======================================================"
echo "Model: {MODEL_NAME}"
echo "Snapshots saved in: {HERACLES_DIR}"
echo "CMFGEN data created at: {MODEL_DIR}/cmfgen/SN_HYDRO_DATA"
echo "======================================================"
