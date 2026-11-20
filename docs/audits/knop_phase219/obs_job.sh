#!/bin/bash
#SBATCH --job-name=obs_p219
#SBATCH --output=obs_p219_%j.out
#SBATCH --error=obs_p219_%j.err
#SBATCH --time=1:59:00
#SBATCH --partition=interactive
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
set -euo pipefail
cd /ptmp/sbrennan/SN2024pba/obs_p219/run_diag_on/obs
module purge
module load intel/2024.0
module load mkl/2026.0
export OMP_NUM_THREADS=8
export cmfdist="/viper/ptmp1/sbrennan/SN2024pba/cmfgen_jfnk"
export CMF_OBS_EXE="/viper/ptmp1/sbrennan/SN2024pba/cmfgen_jfnk/unit_tests/bin/cmf_flux_p219"
export ATOMIC="/ptmp/sbrennan/SN2024pba/atomic"
export atomic="/ptmp/sbrennan/SN2024pba/atomic"
cp -f "$cmfdist/obs_template/batobs.sh" batobs.sh
cp -f "$cmfdist/obs_template/CMF_FLUX_PARAM_INIT" CMF_FLUX_PARAM_INIT
chmod +x batobs.sh
./batobs.sh
