#!/bin/bash
#SBATCH --job-name={MODEL_NAME}_sirocco
#SBATCH --output=sirocco_%j.out
#SBATCH --error=sirocco_%j.err
#SBATCH --time={TIME_LIMIT}
#SBATCH --partition={PARTITION}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={CPUS}
#SBATCH --mem={MEMORY}
#SBATCH --chdir={WORK_DIR}

set -euo pipefail

export SIROCCO_HOME=/ptmp/sbrennan/heracles/explosion_hydro/sirocco/sirocco
export PATH=$SIROCCO_HOME/bin:$PATH

# Load required modules for SIROCCO
module purge
module load intel/2024.0
module load openmpi/4.1
module load gsl/2.6

export LD_LIBRARY_PATH=$(gsl-config --prefix)/lib:$LD_LIBRARY_PATH

echo "Running SIROCCO for {MODEL_NAME} model"
echo "SIROCCO_HOME: $SIROCCO_HOME"
echo "Working directory: $(pwd)"
echo "Date: $(date)"

# Copy atomic data files if not present
if [ ! -f "standard80_sn_kurucz.dat" ]; then
    echo "Copying atomic data files..."
    cp /ptmp/sbrennan/heracles/explosion_hydro/sirocco/supernova_400d/standard80.dat standard80_sn_kurucz.dat
fi
if [ ! -d "data/atomic" ]; then
    mkdir -p data/atomic
    cp -r /ptmp/sbrennan/heracles/explosion_hydro/sirocco/sirocco/xdata/atomic/* data/atomic/
fi

# Convert CMFGEN SN_HYDRO_DATA to SIROCCO import format if needed
cmfgen_data="{CMFGEN_DATA}"
import_file="{IMPORT_FILE}"
if [ -f "$cmfgen_data" ] && [ ! -f "$import_file" ]; then
    echo "Converting CMFGEN SN_HYDRO_DATA to SIROCCO import format..."
    source /mpcdf/soft/RHEL_9/packages/x86_64/python-waterboa/2024.06/bin/activate /ptmp/sbrennan/conda-envs/astro
    unset PYTHONPATH PYTHONHOME
    python3 /ptmp/sbrennan/SN2024pba/scripts/convert_cmfgen_to_sirocco.py "$cmfgen_data" "$import_file"
    if [ $? -ne 0 ]; then
        echo "ERROR: Conversion to SIROCCO format failed"
        exit 1
    fi
    echo "Conversion complete: $import_file"
fi

if [ ! -f "$import_file" ]; then
    echo "ERROR: SIROCCO import file not found: $import_file"
    exit 1
fi

# Run SIROCCO spectrum calculation (MPI-parallel)
echo "Running SIROCCO with $SLURM_NTASKS MPI tasks..."
$SIROCCO_HOME/bin/sirocco {MODEL_NAME}.pf

echo "SIROCCO completed"
echo "Date: $(date)"
