#!/bin/bash
#SBATCH --job-name=ut_xmpos
#SBATCH --output=ut_xmpos_%j.out
#SBATCH --error=ut_xmpos_%j.err
#SBATCH --time=1:00:00
#SBATCH --partition=interactive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --chdir=/viper/ptmp1/sbrennan/SN2024pba/cmfgen_jfnk/unit_tests

set -uo pipefail
cd /viper/ptmp1/sbrennan/SN2024pba/cmfgen_jfnk/unit_tests

# Compute-node python needs the same modules as the login shell
# (numpy/matplotlib for the test harness, ifort + mkl for the drivers).
source /etc/profile.d/modules.sh 2>/dev/null || true
module load intel/2024.0 mkl/2026.0 python-waterboa/2024.06

# The mkl module does not propagate its lib path to the interactive
# node's runtime loader; add it explicitly or every driver binary
# fails on libmkl_intel_lp64.so.3.
export LD_LIBRARY_PATH="${MKLROOT}/lib/intel64:${MKLROOT}/lib:${LD_LIBRARY_PATH:-}"

echo "=== unit tests covering the XM_POS bound and gate changes"
for t in test_mom_solve test_mom_sweep_table test_vadat_keys \
         test_knop_lagged_field test_iter_conv_gate \
         test_gs_drift_diag test_knop_status_text test_knop_ray \
         test_sob_rh_couple; do
  echo "================ $t"
  python3 "python/$t.py"
  rc=$?
  if [ "$rc" -ne 0 ]; then echo "FAILED: $t rc=$rc"; fi
done

echo "================ FP-trap seed driver (check all, fpe0)"
python3 - <<'PY'
import sys
sys.path.insert(0, ".")
sys.path.insert(0, "python")
import common

common.FFLAGS += ["-check", "all", "-fpe0"]
from test_mom_solve import HELPERS
from common import F90DIR, TREE, WORKDIR, compile_objs, link, run

objs = compile_objs(
    "mom_solve_fptrap",
    [TREE / h for h in HELPERS]
    + [F90DIR / "stub_vadat.f", F90DIR / "drv_mom_solve.f"],
)
exe = link("mom_solve_fptrap", objs)
r = run(exe, "seed", "nonj", "3", "1", "diff", "none", "0",
        "t_seed_fptrap", cwd=WORKDIR, check=False)
print(r.stdout[-4000:])
print(r.stderr[-2000:])
ok = r.returncode == 0 and "RESULT_OK" in r.stdout \
    and "SEED_XMNEG=" in r.stdout
print("FPTRAP_SEED_OK=", ok, " rc=", r.returncode)
sys.exit(0 if ok else 1)
PY
echo "fptrap rc=$?"
