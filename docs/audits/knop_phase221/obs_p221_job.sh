#!/bin/bash
#SBATCH --job-name=obs_p221
#SBATCH --output=obs_p221_%j.out
#SBATCH --error=obs_p221_%j.err
#SBATCH --time=23:59:00
#SBATCH --partition=general
#SBATCH --qos=n0004
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
set -euo pipefail
ARM="${OBS_ARM:?OBS_ARM must be run_diag_on or run_diag_off}"
cd "/viper/ptmp1/sbrennan/SN2024pba/obs_p221/${ARM}/obs"
module purge
module load intel/2024.0
module load mkl/2026.0
export OMP_NUM_THREADS=8
export LD_LIBRARY_PATH="/mpcdf/soft/RHEL_9/packages/x86_64/intel_oneapi/2026.0/mkl/2026.0/lib:/mpcdf/soft/RHEL_9/packages/x86_64/intel_oneapi/2026.0/compiler/2026.0/lib:${LD_LIBRARY_PATH:-}"
export cmfdist="/viper/ptmp1/sbrennan/SN2024pba/cmfgen_jfnk"
export CMF_OBS_EXE="/viper/ptmp1/sbrennan/SN2024pba/cmfgen_jfnk/exe/cmf_flux.exe"
export ATOMIC="/ptmp/sbrennan/SN2024pba/atomic"
export atomic="/ptmp/sbrennan/SN2024pba/atomic"
cp -f "$cmfdist/obs_template/batobs.sh" batobs.sh
chmod +x batobs.sh
./batobs.sh
