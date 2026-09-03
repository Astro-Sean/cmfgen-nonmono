"""X_newgrid: hierarchical optical-depth-based CMFGEN grid construction.

Implements the Dessart/Hillier-style grid strategy for the SN2024pba
two-event interacting-supernova model:

    1. Identify physical boundaries / protected features
    2. Build a logarithmic-tau skeleton  (~N_decades * points_per_decade)
    3. Enforce velocity-spacing constraints (incl. dv/dr=0 extrema)
    4. Protect shocks / composition boundaries / velocity extrema
    5. Refine steep T / rho regions under an explicit priority hierarchy
    6. Coarsen redundant points (anchors never removed)
    7. Write SN_HYDRO_DATA

The two-pass CMFGEN ionisation-front refinement (Pass 2) is provided by the
standalone diagnostic script ``scripts/x_newgrid_ionfront_refine.py``.

This package is a *new, standalone* implementation that reuses the
``grid_mapping`` package for I/O, conservative remapping, feature detection
and SN_HYDRO_DATA writing, but does NOT modify any existing working code.
It is selected independently of the v1/v2 converter selection in
``run_model_grid.py``.

Conventions
-----------
All radius arrays are outer-to-inner (decreasing R), matching ``grid_mapping``
and CMFGEN's internal depth indexing (depth 1 = outermost).
"""

from .config import XNewGridConfig
from .skeleton import build_tau_skeleton, n_tau_decades
from .hierarchy import hierarchical_refine, PriorityLevel
from .coarsen import coarsen_grid
from .error_metric import local_interpolation_error, ErrorReport
from .builder import build_grid, build_grid_report

__version__ = "0.1.0"
