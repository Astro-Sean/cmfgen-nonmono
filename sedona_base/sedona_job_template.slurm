#!/bin/bash
#SBATCH --job-name={MODEL_NAME}_sedona
#SBATCH --output=sedona_%j.out
#SBATCH --error=sedona_%j.err
#SBATCH --time={TIME_LIMIT}
#SBATCH --partition={PARTITION}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={CPUS}
#SBATCH --mem={MEMORY}
#SBATCH --chdir={WORK_DIR}

set -euo pipefail

# Use astro conda environment with numpy
source /mpcdf/soft/RHEL_9/packages/x86_64/python-waterboa/2024.06/bin/activate /ptmp/sbrennan/conda-envs/astro
unset PYTHONPATH PYTHONHOME

# Load required modules for SEDONA
module use /mpcdf/soft/RHEL_9/sub/intel_2024_0/modules/mpi
module use /mpcdf/soft/RHEL_9/modules/compilers
module use /mpcdf/soft/RHEL_9/modules/third-party-compilers

module purge
module load intel/2024.0
module load openmpi/4.1
module load gsl/2.6
module load hdf5-serial/1.14.1

export SEDONA_HOME=/ptmp/sbrennan/sedona

# Set library paths using loaded modules
export LD_LIBRARY_PATH=$(gsl-config --prefix)/lib:$(dirname $(dirname $(which h5cc)))/lib:$(pwd):$LD_LIBRARY_PATH

# Set OpenMP threads
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

echo "Running SEDONA for {MODEL_NAME} model"
echo "SEDONA_HOME: $SEDONA_HOME"
echo "Working directory: $(pwd)"
echo "Date: $(date)"

# Convert CMFGEN SN_HYDRO_DATA to SEDONA .mod format if needed
cmfgen_data="{CMFGEN_DATA}"
mod_file="{MOD_FILE}"
if [ -f "$cmfgen_data" ] && [ ! -f "$mod_file" ]; then
    echo "Converting CMFGEN SN_HYDRO_DATA to SEDONA .mod format..."
    python3 /ptmp/sbrennan/sedona/convert_cmfgen_to_sedona_v2.py "$cmfgen_data" "$mod_file"
    if [ $? -ne 0 ]; then
        echo "ERROR: Conversion failed"
        exit 1
    fi
    echo "Conversion complete"
fi

# Run SEDONA spectrum calculation
$SEDONA_HOME/src/sedona6.ex param.lua

echo "SEDONA completed - generating plots"

# Auto-plot the spectrum
python3 plot_spectrum.py

echo "Plots generated"
echo "Date: $(date)"
