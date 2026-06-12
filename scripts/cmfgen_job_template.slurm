#!/bin/bash
#SBATCH --job-name=cmfgen_{MODEL_NAME}
#SBATCH --output={MODEL_DIR}/{MODEL_NAME}_%j.out
#SBATCH --error={MODEL_DIR}/{MODEL_NAME}_%j.out
#SBATCH --partition={PARTITION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={CPUS}
#SBATCH --mem={MEMORY}
#SBATCH --time={TIME}

# CMFGEN simulation for model: {MODEL_NAME}
# Using HERACLES snapshot at {TARGET_DAYS} days

set -euo pipefail

# Use astro conda environment with numpy
source /mpcdf/soft/RHEL_9/packages/x86_64/python-waterboa/2024.06/bin/activate /ptmp/sbrennan/conda-envs/astro
unset PYTHONPATH PYTHONHOME

MODEL_NAME="{MODEL_NAME}"
MODEL_DIR="{MODEL_DIR}"
CMFGEN_DIR="{MODEL_DIR}/cmfgen"
CMFGEN_BASE_DIR="{BASE_DIR}"
CMFGEN_MODEL_DIR="{CMFGEN_DIR}"
PROJECT_DIR="{PROJECT_DIR}"

cd "${MODEL_DIR}"

echo "======================================================"
echo "CMFGEN Simulation: ${MODEL_NAME}"
echo "======================================================"
echo "JobID: ${SLURM_JOB_ID}"
echo "CMFGEN directory: ${CMFGEN_DIR}"
echo "CMFGEN model directory: ${CMFGEN_MODEL_DIR}"
echo "Time: $(date)"
echo "======================================================"

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

module purge
module load gcc/13 openmpi/4.1 mkl/2025.3
module load python-waterboa/2025.06

# Explicitly set MKL library path (workaround for module loading issues on compute nodes)
export LD_LIBRARY_PATH=/viper/u1/system/soft/RHEL_9/packages/x86_64/intel_oneapi/2025.3/mkl/2025.3/lib/intel64:$LD_LIBRARY_PATH

# Copy base directory to cmfgen directory (fresh setup for this run)
echo "Copying base directory from ${CMFGEN_BASE_DIR} to ${CMFGEN_DIR}..."
if [[ -d "${CMFGEN_DIR}" ]]; then
  echo "Removing existing cmfgen directory..."
  rm -rf "${CMFGEN_DIR}"
fi
cp -r "${CMFGEN_BASE_DIR}" "${CMFGEN_DIR}"

# Copy SN_HYDRO_DATA to cmfgen directory
echo "Copying SN_HYDRO_DATA to cmfgen directory..."
cp SN_HYDRO_DATA "${CMFGEN_DIR}/"

cd "${CMFGEN_DIR}"

# Run CMFGEN preparation scripts
echo "Running CMFGEN preparation..."

# Prepare IN_INTS
if [[ -f "{PROJECT_DIR}/scripts/prepare_IN_INTS.py" ]]; then
  python3 "{PROJECT_DIR}/scripts/prepare_IN_INTS.py" || true
fi

# Prepare MODEL_SPEC
if [[ -f "{PROJECT_DIR}/scripts/prepare_MODEL_SPEC.py" ]]; then
  python3 "{PROJECT_DIR}/scripts/prepare_MODEL_SPEC.py" || true
fi

# Update VADAT from SN_HYDRO_DATA
if [[ -f "/ptmp/sbrennan/scripts/update_VADAT_from_SN_HYDRO_DATA.py" ]]; then
  python3 /ptmp/sbrennan/scripts/update_VADAT_from_SN_HYDRO_DATA.py || true
fi

# Submit CMFGEN jobs based on configuration
echo "Submitting CMFGEN jobs..."
SUBMIT_SCRIPT="{SUBMIT_SCRIPT}"

if [[ "{RUN_INIT}" == "true" ]]; then
  echo "Submitting CMFGEN INIT job..."
  if [[ "{RUN_STABLE}" == "true" ]]; then
    "${SUBMIT_SCRIPT}" --all
  else
    "${SUBMIT_SCRIPT}" --all
  fi
else
  echo "Submitting CMFGEN FINAL job (restart mode)..."
  "${SUBMIT_SCRIPT}" --restart
fi

echo ""
echo "======================================================"
echo "CMFGEN jobs submitted"
echo "======================================================"
echo "Model: ${MODEL_NAME}"
echo "Monitor with: squeue -u $USER"
echo "======================================================"
