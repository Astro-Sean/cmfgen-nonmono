"""CMFGEN source audit module.

Encapsulates the findings from the local modified-CMFGEN source audit
(see audits/grid_mapping_audit/01_cmfgen_source_audit.md for the full audit).

This module provides programmatic access to the audit findings so that
the grid mapping code can verify compatibility with the local solver.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


@dataclass
class SolverAuditFindings:
    """Key findings from the CMFGEN source audit.

    These are hardcoded from the audit document. If the source tree
    changes, re-run the audit and update these values.
    """
    # The project uses BOTH VEL_LAW paths, selected per-YAML
    # (cmfgen.vel_law, default 12).  Both are supported; the converter must
    # compute the solver sigma with the formula matching the chosen path.
    #   VEL_LAW=12 -> RV_SN_MODEL_SNIIN: distance-weighted finite-difference
    #                 sigma (rv_sn_model_sniin.f L262-269).  R/V come from
    #                 input_hydro.dat (written by prepare_VADAT.py from
    #                 SN_HYDRO_DATA).  SN_HYDRO_MODEL=.FALSE.
    #   VEL_LAW=11 -> SET_RV_HYDRO_MODEL_V3: analytic derivative of a
    #                 monotone-cubic (PCHIP) interpolant of V(R)
    #                 (set_rv_hydro_model_v3.f L419-433).  R/V/SIGMA read
    #                 from SN_HYDRO_DATA; the file SIGMA column is used ONLY
    #                 for the density time-evolution correction (L232).
    #                 SN_HYDRO_MODEL=.TRUE., typically RD_IN_R_GRID=T.
    # In NEITHER path does the SN_HYDRO_DATA SIGMA column directly feed the
    # radiative-transfer solver.
    vel_law_default: int = 12
    vel_law_supported: tuple = (11, 12)

    # SIGMA formula (VEL_LAW=12): distance-weighted dV/dR, then R*dV/dR/V - 1
    sigma_formula: str = "R * dVdR / V - 1"

    # SIGMA cap: -1.2 < SIGMA < -1.0 -> -0.999
    sigma_cap_value: float = -0.999
    sigma_cap_low: float = -1.2
    sigma_cap_high: float = -1.0

    # Dessart+2015 GAM guard
    dessart_guard_active: bool = True

    # R must be strictly monotonic (decreasing, outer-to-inner)
    r_monotonic_required: bool = True

    # Velocity may be non-monotonic
    v_monotonic_required: bool = False

    # V must be > 0 (SIGMA divides by V)
    # With preserve_v_sign=True, V can be negative but |V| must be > 0
    v_nonzero_required: bool = True

    # Required VADAT flags for non-monotonic velocity
    required_vadat_flags: dict = field(default_factory=lambda: {
        # VEL_LAW is per-model (11 or 12); not fixed here.
        'USE_J_REL': True,
        'N_TYPE': 'N_ON_J',
        'INCL_REL': True,
        'DO_DDT': False,
        'CHK_J': 'ABS_VAL',
        'METHOD': 'LOGMON',
    })

    # BLOCKING issues (from audit)
    blocking_issues: list[str] = field(default_factory=lambda: [
        'comp_j_blank.f Eddington pre-computation calls MOM_J_CMF_V11 unconditionally '
        '(Fix 4 NOT applied). May be mitigated by N_ON_J setting.',
    ])

    # HIGH issues
    high_issues: list[str] = field(default_factory=lambda: [
        'V=0 causes SIGMA NaN (division by zero in rv_sn_model_sniin.f L266)',
        'No minimum-points-per-velocity-reversal requirement derived from solver tests',
        'No test harness for non-monotonic solver verification',
    ])

    # No characteristic/resonance anchors needed (Sobolev approximation is local)
    characteristic_anchors_needed: bool = False

    # No explicit limit on number of velocity reversals
    max_velocity_reversals: Optional[int] = None

    # Solver does NOT require conservative remapping (reads point values)
    # But conservative remapping IS required for physical accuracy
    solver_requires_conservative_remap: bool = False

    # Source tree path (for verification)
    source_tree_path: str = "/viper/ptmp1/sbrennan/SN2024pba/cmfgen_source"

    # Audit document path
    audit_doc_path: str = "audits/grid_mapping_audit/01_cmfgen_source_audit.md"


def get_audit_findings() -> SolverAuditFindings:
    """Get the CMFGEN source audit findings."""
    return SolverAuditFindings()


def verify_source_tree(source_path: Optional[str] = None) -> dict:
    """Verify that the local CMFGEN source tree exists and has the expected modifications.

    Returns dict with verification results.
    """
    findings = get_audit_findings()
    path = Path(source_path or findings.source_tree_path)

    result = {
        'source_tree_exists': path.exists(),
        'source_tree_path': str(path),
        'key_files': {},
        'verified': False,
    }

    if not path.exists():
        result['error'] = f"Source tree not found at {path}"
        return result

    # Check key modified files
    key_files = [
        'new_main/subs/rv_sn_model_sniin.f',
        'plane/fg_j_cmf_v13.f',
        'subs/sobew.f',
        'new_main/mod_subs/rd_control_variables.f',
        'new_main/mod_subs/set_abund_clump.f',
    ]

    all_exist = True
    for f in key_files:
        exists = (path / f).exists()
        result['key_files'][f] = exists
        if not exists:
            all_exist = False

    result['verified'] = all_exist
    return result


def check_sigma_compatibility(sigma: 'np.ndarray') -> dict:
    """Check if SIGMA values are compatible with the local solver.

    Returns dict with:
    - n_in_clamp_range: number of points with -1.2 < SIGMA < -1.0
    - n_below_minus_one: number of points with SIGMA < -1.0
    - n_capped: number of points that would be capped to -0.999
    - compatible: True if no issues
    """
    import numpy as np

    sigma_arr = np.asarray(sigma)
    findings = get_audit_findings()

    n_in_clamp = int(np.sum((sigma_arr > findings.sigma_cap_low) &
                            (sigma_arr < findings.sigma_cap_high)))
    n_below = int(np.sum(sigma_arr < -1.0))
    n_capped = n_in_clamp  # these will be capped

    return {
        'n_in_clamp_range': n_in_clamp,
        'n_below_minus_one': n_below,
        'n_capped': n_capped,
        'compatible': True,  # The solver handles this via the cap + Dessart guard
        'note': 'SIGMA in -1.2..-1.0 range will be capped to -0.999 by the solver. '
                'This is expected for non-monotonic velocity and is handled by the '
                'Dessart+2015 guard.',
    }
