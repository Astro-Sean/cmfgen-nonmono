#!/bin/bash
# [DEPRECATED] HERACLES-to-CMFGEN conversion job template
# DEPRECATED (2026-07-09): HERACLES does not fit the science case for this project
# (interacting Type IIn supernovae). The project now uses the STELLA -> CMFGEN
# pipeline exclusively. For the active conversion path, see stella_job_template.slurm
# (STELLA -> CMFGEN conversion is embedded in the STELLA workflow).
#SBATCH --job-name=convert_{MODEL_NAME}
#SBATCH --output={MODEL_DIR}/convert_{MODEL_NAME}_%j.out
#SBATCH --error={MODEL_DIR}/convert_{MODEL_NAME}_%j.out
#SBATCH --partition={PARTITION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={CPUS}
#SBATCH --mem={MEMORY}
#SBATCH --time={TIME}
#SBATCH --dependency=afterok:{HERACLES_JOB_ID}

# CMFGEN conversion for model: {MODEL_NAME}
# Dependent on HERACLES job: {HERACLES_JOB_ID}

set -euo pipefail

# Use astro conda environment with numpy
source /mpcdf/soft/RHEL_9/packages/x86_64/python-waterboa/2024.06/bin/activate /ptmp/sbrennan/conda-envs/astro
unset PYTHONPATH PYTHONHOME

MODEL_NAME="{MODEL_NAME}"
MODEL_DIR="{MODEL_DIR}"
PROJECT_DIR="{PROJECT_DIR}"
CMFGEN_DIR="{MODEL_DIR}/cmfgen"

cd {MODEL_DIR}

echo "======================================================"
echo "CMFGEN Conversion: ${MODEL_NAME}"
echo "======================================================"
echo "JobID: ${SLURM_JOB_ID}"
echo "Dependent on HERACLES job: {HERACLES_JOB_ID}"
echo "Time: $(date)"
echo "======================================================"

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

module purge
module load gcc/13 openmpi/4.1 mkl
module load python-waterboa/2025.06

# Set up conversion
export HERACLES_ROOT="{HERACLES_ROOT}"
CONVERSION_SCRIPT="{PROJECT_DIR}/scripts/convert_heracles_to_cmfgen_wrapper.py"
# NOTE: This template is currently unused (conversion is embedded in heracles_job_template.slurm).
# CONFIG_FILE is preserved here in case this template is re-enabled in future.
CONFIG_FILE="{PROJECT_DIR}/configs/{MODEL_NAME}.yaml"

# Run conversion
echo "Converting HERACLES snapshot to CMFGEN format..."
python3 "${CONVERSION_SCRIPT}" "${CONFIG_FILE}"

# Clip sigma values to prevent CMFGEN crashes
# Sigma = dlnV/dlnr - 1; values < -0.99 cause density scaling divergence in set_rv_hydro_model_v3.f
echo "Clipping sigma values in SN_HYDRO_DATA to [-0.99, +20.0]..."
python3 << 'EOF'
import re
import numpy as np

sn_hydro = "{MODEL_DIR}/cmfgen/SN_HYDRO_DATA"
with open(sn_hydro, 'r') as f:
    content = f.read()

sigma_label = 'Sigma (dlnV/dlnr-1)'
sigma_start_idx = content.find(sigma_label)
if sigma_start_idx == -1:
    print("  WARNING: Sigma section not found in SN_HYDRO_DATA - skipping clip")
else:
    sigma_data_start = sigma_start_idx + len(sigma_label) + 1
    next_section_idx = content.find('\nTemperature', sigma_data_start)
    
    sigma_block_text = content[sigma_data_start : next_section_idx]
    vals = re.findall(r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eEdD][-+]?\d+)?', sigma_block_text)
    sigma = np.array([float(v.replace('D','E').replace('d','e')) for v in vals])
    
    SIGMA_MIN = -0.99
    SIGMA_MAX = 20.0
    sigma_clipped = np.clip(sigma, SIGMA_MIN, SIGMA_MAX)
    n_clipped = int(np.sum(sigma != sigma_clipped))
    
    if n_clipped > 0:
        print(f"  Clipped {n_clipped} sigma values. Range: [{sigma.min():.4f}, {sigma.max():.4f}] -> [{sigma_clipped.min():.4f}, {sigma_clipped.max():.4f}]")
        
        def fmt_array(arr, n_per_line=8):
            lines = []
            for i in range(0, len(arr), n_per_line):
                chunk = arr[i:i+n_per_line]
                lines.append('   ' + '   '.join(f'{v:.15E}' for v in chunk))
            return '\n'.join(lines) + '\n'
        
        new_sigma_block = fmt_array(sigma_clipped)
        new_content = (content[:sigma_data_start] + new_sigma_block + content[next_section_idx:])
        
        with open(sn_hydro, 'w') as f:
            f.write(new_content)
        print(f"  Sigma clipping applied to {sn_hydro}")
    else:
        print(f"  No sigma clipping needed (all values in range)")
EOF

# Copy SN_HYDRO_DATA to cmfgen directory
echo "Copying SN_HYDRO_DATA to cmfgen directory..."
SN_HYDRO_DATA="{MODEL_DIR}/cmfgen/SN_HYDRO_DATA"
cp "${SN_HYDRO_DATA}" "{CMFGEN_DIR}/"

echo ""
echo "======================================================"
echo "CMFGEN conversion complete"
echo "======================================================"
echo "Model: ${MODEL_NAME}"
echo "CMFGEN directory: {CMFGEN_DIR}"
echo "======================================================"
