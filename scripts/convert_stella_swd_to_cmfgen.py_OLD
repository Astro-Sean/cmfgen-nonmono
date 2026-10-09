"""Convert a STELLA *.swd snapshot to CMFGEN SN_HYDRO_DATA format.

STELLA *.swd file format (ASCII, written by PRIBAL in lbalsw.trf):
  Each epoch block: one line per zone, columns:
    t_day   km   log10(M_ext_Msun)   log10(R_cm)   V_kms
    log10(T_gas_K)   log10(T_rad_K)   log10(rho_gcc)
    log10(P_cgs)   log10(Qv)   log10(eng)   Flum_1e40   WRKX

  The km=1 (outer boundary) line has t_day as first field;
  all interior lines start with 0 or blank in that column.

  Lines with km=1 or km=Nzon are "header" lines (t_day appears);
  interior lines start with a leading zero or spaces in the t field.

Usage:
  python convert_stella_swd_to_cmfgen.py SWD_FILE OUTPUT_FILE [OPTIONS]
    --target-days T    Pick the epoch closest to T days (default: last epoch)
    --tau-max X        Outer optical depth cutoff (default: 30)
    --tau-min X        Inner optical depth cutoff (default: 1e-3)
    --nzon-target N    Target number of output grid points (default: 200)
    --time-days T      Explosion time in days (default: read from swd)
    --model-name NAME  For header comment (default: stella)
    --kappa K          Electron-scattering opacity cm^2/g (default: 0.34)
    --xni-mass X       56Ni mass in Msun for luminosity estimate (default: 0)
    --rho-file PATH    *.rho file from STELLA for final evolved composition (preferred)
    --abn-file PATH    *.abn file from STELLA for initial composition (fallback)
    --enable-adaptive-grid  Enable Dessart-style adaptive grid refinement
    --feature-aware-downsample  Preserve high-gradient regions (default: enabled)
    --no-feature-aware-downsample  Disable feature-aware downsampling

The output SN_HYDRO_DATA is in CMFGEN outer->inner (decreasing r) order.
"""

# The implementation was split into the prepare_SN_HYDRO_DATA package
# (prepare_SN_HYDRO_DATA/utils/, one module per pipeline stage). This
# file is the stable entry point: it re-exports the full public API so
# all existing callers keep working unchanged:
#   from convert_stella_swd_to_cmfgen import X          (sys.path=scripts/)
#   from scripts.convert_stella_swd_to_cmfgen import X  (project-root import)
#   importlib.spec_from_file_location on this file      (run_model_grid.py)
#   python scripts/convert_stella_swd_to_cmfgen.py ...  (CLI)

import sys
from pathlib import Path

_PKG_ROOT = Path(__file__).resolve().parent.parent
if str(_PKG_ROOT) not in sys.path:
    sys.path.insert(0, str(_PKG_ROOT))

from prepare_SN_HYDRO_DATA.utils import *  # noqa: F401,F403,E402
from prepare_SN_HYDRO_DATA.utils import __all__ as _utils_all  # noqa: E402
from prepare_SN_HYDRO_DATA.utils.driver import main  # noqa: F401,E402

__all__ = list(_utils_all) + ["main"]

if __name__ == "__main__":
    main()
