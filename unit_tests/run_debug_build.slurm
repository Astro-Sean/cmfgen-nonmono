#!/bin/bash
#SBATCH --job-name=dbg_build
#SBATCH --output=dbg_build_%j.out
#SBATCH --error=dbg_build_%j.err
#SBATCH --time=1:55:00
#SBATCH --partition=interactive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --chdir=/viper/ptmp1/sbrennan/SN2024pba/cmfgen_jfnk

set -uo pipefail
cd /viper/ptmp1/sbrennan/SN2024pba/cmfgen_jfnk

source /etc/profile.d/modules.sh 2>/dev/null || true
module load intel/2024.0 mkl/2026.0
export LD_LIBRARY_PATH="${MKLROOT}/lib/intel64:${MKLROOT}/lib:${LD_LIBRARY_PATH:-}"

echo "=== build_debug.sh (V2: full debug executable) ==="
./build_debug.sh > dbg_build_make.log 2>&1
rc=$?
tail -15 dbg_build_make.log
echo "DBG_BUILD_RC=$rc"
if [ "$rc" -eq 0 ]; then
  ls -la exe_dbg/
fi
