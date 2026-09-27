"""
STELLA-to-CMFGEN grid mapping with non-monotonic velocity support.

This package implements conservative shell-overlap remapping, feature-preserving
grid selection, and validation for the locally modified CMFGEN non-monotonic
velocity solver (VEL_LAW=12, Dessart+2015 SIGMA guard).

Modules:
    source_audit       - Local modified-CMFGEN source and convention audit
    stella_io          - Parse STELLA .swd snapshots with units and metadata
    shell_model        - Conservative shell representation (STELLA + CMFGEN grids)
    solver_metrics     - Exact local CMFGEN velocity-gradient diagnostics
    feature_detection  - Shocks, interfaces, velocity extrema, tau landmarks
    conservative_remap - Shell-overlap remapping for density/composition/isotopes
    opacity_tau        - Opacity validation and optical-depth reintegration
    ionization_metrics - Ion fractions, Ne, Ye, opacity, and front diagnostics
    grid_selection     - Mandatory anchors + constrained priority refinement
    validation         - Conservation, ordering, feature, and solver checks
    cmfgen_writer      - Atomic, validated SN_HYDRO_DATA output
    reporting          - JSON, text, and plot diagnostics

See audits/grid_mapping_audit/ for the source audit and code audit that motivated
this refactor.
"""

__version__ = "0.1.0"

from .shell_model import SourceShellGrid, TargetGrid
from .stella_io import read_swd, select_epoch, read_abn, read_rho, read_xni, build_source_grid
from .solver_metrics import (sigma_velocity, sigma_cmfgen_v11, sigma_solver,
                             apply_sigma_cap, count_velocity_reversals)
from .feature_detection import detect_all_features, FeatureReport
from .conservative_remap import conservative_remap, compute_conservation_residuals
from .opacity_tau import compute_tau, find_tau_landmarks, find_photosphere
from .ionization_metrics import saha_electron_density, compute_ionization_diagnostics
from .grid_selection import GridSelector, PointBudgetPolicy
from .validation import validate_target_grid, ValidationReport
from .cmfgen_writer import write_sn_hydro_data
from .sn_hydro_io import read_sn_hydro_data, validate_sn_hydro_data, SnHydroData
from .reporting import (write_grid_audit, write_grid_features,
                        write_grid_validation_json, write_grid_validation_txt,
                        write_profile_plots)
from .source_audit import get_audit_findings, verify_source_tree, check_sigma_compatibility
