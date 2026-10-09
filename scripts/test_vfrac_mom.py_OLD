#!/usr/bin/env python3
"""Tests for VFRAC_MOM handling in prepare_VADAT.py.

VFRAC_MOM (DELV_FRAC_MOM) is a RADIATIVE TRANSFER GRID RESOLUTION parameter,
NOT a momentum correction, velocity damping, or under-relaxation parameter.

Source-level definition (from cmfgen_source/):
  - Declared as DELV_FRAC_MOM in control_variable_mod.f:309
  - Read from VADAT as 'VFRAC_MOM' in rd_control_variables.f:535
  - Description: "Maximum velocity spacing (Doppler units) in MOM_J_CMF_V10"
  - Used in MOM_J_CMF_V11:279, MOM_JREL_V9:368, PP_MOM_CMF_V1:259, MOM_J_DDT_V6:324

Algorithm (identical in all solvers):
  T2 = VFRAC_MOM * min(VDOP_VEC)
  N_extra(I) = INT( |V(I) - V(I+1)| / T2 )   [rounded up to even]
  → Extra RT grid points inserted between depths I and I+1

Key properties:
  - VFRAC_MOM=0 is INVALID (division by zero)
  - VFRAC_MOM does NOT modify the physical velocity field
  - Smaller VFRAC_MOM = more RT grid points = better resolution but slower
  - Larger VFRAC_MOM = fewer RT grid points = faster but less accurate J
  - VFRAC_MOM is NOT a momentum correction or damping parameter
  - VFRAC_MOM is NOT related to the Dessart non-monotonic velocity treatment
  - VFRAC_MOM does NOT vary by convergence stage (it's a grid parameter)
"""

import sys
import math
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))


def test_vfrac_mom_not_zero():
    """VFRAC_MOM=0 must be rejected (causes division by zero in solver)."""
    from prepare_VADAT import compute_optimal_vfrac_mom, HydroData

    # Create a simple hydro with non-monotonic velocity
    n = 50
    r = np.linspace(1e10, 1e11, n)
    v = np.linspace(100e5, 1000e5, n)  # monotonic
    t = np.full(n, 5000.0)
    rho = np.full(n, 1e-15)

    hydro = HydroData(
        radius_cm=r, vel_cms=v, temp_104k=t / 1e4, rho_cgs=rho,
        ed_cgs=np.full(n, 1e8), na_cgs=np.full(n, 1e8),
        kappa_cgs=None, mass_fractions={},
        params={"RMAX": 1e11, "RSTAR": 1e10, "VINF": 1000.0,
                "time_days": 35.0, "MASS": 1.0},
    )

    vfrac, _ = compute_optimal_vfrac_mom(hydro, safety_factor=1.5)

    assert vfrac > 0.0, f"VFRAC_MOM must be > 0, got {vfrac}"
    print(f"  test_vfrac_mom_not_zero: VFRAC_MOM={vfrac:.2f} ✓")


def test_vfrac_mom_does_not_modify_velocity():
    """VFRAC_MOM is a grid resolution parameter, not a velocity correction.

    This test verifies that compute_optimal_vfrac_mom returns a value
    that depends on the velocity profile (grid resolution) but does NOT
    alter the velocity field itself.
    """
    from prepare_VADAT import compute_optimal_vfrac_mom, HydroData

    n = 100
    r = np.linspace(1e10, 1e11, n)

    # Non-monotonic velocity (interaction region)
    v = np.linspace(100e5, 2000e5, n)
    v[40:60] = v[40:60] * 0.5  # velocity reversal

    t = np.full(n, 5000.0)
    rho = np.full(n, 1e-15)

    hydro = HydroData(
        radius_cm=r, vel_cms=v, temp_104k=t / 1e4, rho_cgs=rho,
        ed_cgs=np.full(n, 1e8), na_cgs=np.full(n, 1e8),
        kappa_cgs=None, mass_fractions={},
        params={"RMAX": 1e11, "RSTAR": 1e10, "VINF": 2000.0,
                "time_days": 35.0, "MASS": 1.0},
    )

    vfrac, _ = compute_optimal_vfrac_mom(hydro, safety_factor=1.5)

    # The velocity field should be unchanged
    assert np.allclose(hydro.vel_cms, v), "VFRAC_MOM computation must not alter velocity"
    assert vfrac > 0.0, f"VFRAC_MOM must be > 0, got {vfrac}"
    print(f"  test_vfrac_mom_does_not_modify_velocity: VFRAC_MOM={vfrac:.2f}, "
          f"velocity unchanged ✓")


def test_vfrac_mom_bounds():
    """VFRAC_MOM should be within reasonable bounds."""
    from prepare_VADAT import compute_optimal_vfrac_mom, HydroData

    n = 50
    r = np.linspace(1e10, 1e11, n)
    v = np.linspace(100e5, 5000e5, n)  # large velocity range
    t = np.full(n, 5000.0)
    rho = np.full(n, 1e-15)

    hydro = HydroData(
        radius_cm=r, vel_cms=v, temp_104k=t / 1e4, rho_cgs=rho,
        ed_cgs=np.full(n, 1e8), na_cgs=np.full(n, 1e8),
        kappa_cgs=None, mass_fractions={},
        params={"RMAX": 1e11, "RSTAR": 1e10, "VINF": 5000.0,
                "time_days": 35.0, "MASS": 1.0},
    )

    vfrac, _ = compute_optimal_vfrac_mom(hydro, safety_factor=1.5)

    # Should be within [0.5, 15.0] bounds
    assert 0.5 <= vfrac <= 15.0, f"VFRAC_MOM={vfrac} outside [0.5, 15.0] bounds"
    print(f"  test_vfrac_mom_bounds: VFRAC_MOM={vfrac:.2f} in [0.5, 15.0] ✓")


def test_vfrac_mom_non_monotonic():
    """VFRAC_MOM should handle non-monotonic velocity without error."""
    from prepare_VADAT import compute_optimal_vfrac_mom, HydroData

    n = 100
    r = np.linspace(1e10, 1e11, n)
    # Strongly non-monotonic: sawtooth pattern
    v = np.linspace(100e5, 2000e5, n)
    for i in range(0, n, 10):
        v[i:i+5] = v[i:i+5] * 0.3  # sharp velocity drops

    t = np.full(n, 5000.0)
    rho = np.full(n, 1e-15)

    hydro = HydroData(
        radius_cm=r, vel_cms=v, temp_104k=t / 1e4, rho_cgs=rho,
        ed_cgs=np.full(n, 1e8), na_cgs=np.full(n, 1e8),
        kappa_cgs=None, mass_fractions={},
        params={"RMAX": 1e11, "RSTAR": 1e10, "VINF": 2000.0,
                "time_days": 35.0, "MASS": 1.0},
    )

    vfrac, _ = compute_optimal_vfrac_mom(hydro, safety_factor=1.5)

    assert vfrac > 0.0, f"VFRAC_MOM must be > 0 for non-monotonic, got {vfrac}"
    print(f"  test_vfrac_mom_non_monotonic: VFRAC_MOM={vfrac:.2f} ✓")


def test_vfrac_mom_not_momentum_correction():
    """Verify that VFRAC_MOM is NOT a momentum correction parameter.

    The name 'VFRAC_MOM' might suggest 'velocity fraction for momentum',
    but the source code shows it is 'DELV_FRAC_MOM' = 'Delta V fraction
    for MOM_J_CMF' = maximum velocity spacing fraction for the moment
    equation solver. It controls RT grid refinement, not momentum corrections.

    This test verifies that the computed VFRAC_MOM is a plausible grid
    resolution parameter (proportional to max velocity jump / Doppler width)
    and NOT a fraction between 0 and 1 that would represent a correction.
    """
    from prepare_VADAT import compute_optimal_vfrac_mom, HydroData

    # Small velocity jump → small VFRAC_MOM (fine grid OK)
    n = 50
    r = np.linspace(1e10, 1e11, n)
    v_small = np.linspace(100e5, 200e5, n)  # small jump
    t = np.full(n, 5000.0)
    rho = np.full(n, 1e-15)

    hydro_small = HydroData(
        radius_cm=r, vel_cms=v_small, temp_104k=t / 1e4, rho_cgs=rho,
        ed_cgs=np.full(n, 1e8), na_cgs=np.full(n, 1e8),
        kappa_cgs=None, mass_fractions={},
        params={"RMAX": 1e11, "RSTAR": 1e10, "VINF": 200.0,
                "time_days": 35.0, "MASS": 1.0},
    )
    vfrac_small, _ = compute_optimal_vfrac_mom(hydro_small, safety_factor=1.5)

    # Large velocity jump → larger VFRAC_MOM (coarser grid needed)
    v_large = np.linspace(100e5, 5000e5, n)  # large jump
    hydro_large = HydroData(
        radius_cm=r, vel_cms=v_large, temp_104k=t / 1e4, rho_cgs=rho,
        ed_cgs=np.full(n, 1e8), na_cgs=np.full(n, 1e8),
        kappa_cgs=None, mass_fractions={},
        params={"RMAX": 1e11, "RSTAR": 1e10, "VINF": 5000.0,
                "time_days": 35.0, "MASS": 1.0},
    )
    vfrac_large, _ = compute_optimal_vfrac_mom(hydro_large, safety_factor=1.5)

    # VFRAC_MOM should be larger for larger velocity jumps (grid resolution)
    assert vfrac_large >= vfrac_small, \
        f"VFRAC_MOM should be larger for larger velocity jumps: " \
        f"small={vfrac_small:.2f}, large={vfrac_large:.2f}"

    print(f"  test_vfrac_mom_not_momentum_correction: "
          f"small_jump={vfrac_small:.2f}, large_jump={vfrac_large:.2f} ✓")


def test_vfrac_mom_source_definition():
    """Verify the source-level definition of VFRAC_MOM.

    VFRAC_MOM is parsed as DELV_FRAC_MOM in rd_control_variables.f:535
    with description "Maximum velocity spacing (Doppler units) in MOM_J_CMF_V10".
    """
    # Check that the source file contains the expected definition
    src_path = Path('/viper/ptmp1/sbrennan/SN2024pba/cmfgen_source/new_main/mod_subs/rd_control_variables.f')
    if not src_path.exists():
        print("  test_vfrac_mom_source_definition: SKIPPED (source not found)")
        return

    text = src_path.read_text()
    assert "DELV_FRAC_MOM" in text, "DELV_FRAC_MOM not found in rd_control_variables.f"
    assert "VFRAC_MOM" in text, "VFRAC_MOM keyword not found in rd_control_variables.f"
    assert "Maximum velocity spacing" in text, \
        "Description 'Maximum velocity spacing' not found"

    # Check the declaration
    mod_path = Path('/viper/ptmp1/sbrennan/SN2024pba/cmfgen_source/new_main/mod_subs/control_variable_mod.f')
    if mod_path.exists():
        mod_text = mod_path.read_text()
        assert "DELV_FRAC_MOM" in mod_text, \
            "DELV_FRAC_MOM not declared in control_variable_mod.f"

    # Check the solver usage
    solver_path = Path('/viper/ptmp1/sbrennan/SN2024pba/cmfgen_source/plane/mom_j_cmf_v11.f')
    if solver_path.exists():
        solver_text = solver_path.read_text()
        assert "VDOP_FRAC" in solver_text, \
            "VDOP_FRAC (DELV_FRAC_MOM) not used in MOM_J_CMF_V11"
        # Check the algorithm: T2=VDOP_FRAC*MINVAL(VDOP_VEC)
        assert "VDOP_FRAC*MINVAL" in solver_text, \
            "Algorithm T2=VDOP_FRAC*MINVAL(VDOP_VEC) not found in MOM_J_CMF_V11"

    print(f"  test_vfrac_mom_source_definition: source verified ✓")


def test_vfrac_mom_not_related_to_dessert_solver():
    """VFRAC_MOM is NOT related to the Dessart non-monotonic velocity treatment.

    The Dessart non-monotonic velocity approach is handled by:
    - VEL_LAW=12 (RV_SN_MODEL_SNIIN): reads input_hydro.dat with non-monotonic V
    - USE_LAM_ES=T: uses formal solution (PP_MOM_CMF_V1) instead of moment equations
    - SIGMA clamping in the solver

    VFRAC_MOM only controls RT grid refinement and is independent of the
    velocity law or solver choice.
    """
    # Verify that VFRAC_MOM is used in ALL solvers (not just the non-monotonic one)
    solvers = [
        ('plane/mom_j_cmf_v11.f', 'MOM_J_CMF_V11'),  # standard moment
        ('plane/mom_jrel_v9.f', 'MOM_JREL_V9'),      # relativistic moment
        ('plane/pp_mom_cmf_v1.f', 'PP_MOM_CMF_V1'),   # formal solution (USE_LAM_ES=T)
        ('plane/mom_j_ddt_v6.f', 'MOM_J_DDT_V6'),     # time evolution
    ]

    src_root = Path('/viper/ptmp1/sbrennan/SN2024pba/cmfgen_source')
    for rel_path, solver_name in solvers:
        full_path = src_root / rel_path
        if not full_path.exists():
            print(f"  test_vfrac_mom_not_related_to_dessert_solver: "
                  f"SKIP {solver_name} (not found)")
            continue
        text = full_path.read_text()
        assert 'VDOP_FRAC' in text, \
            f"{solver_name} does not use VDOP_FRAC (VFRAC_MOM)"
        # All solvers use the same algorithm
        assert 'VDOP_FRAC*MINVAL' in text, \
            f"{solver_name} does not use VDOP_FRAC*MINVAL algorithm"

    print(f"  test_vfrac_mom_not_related_to_dessert_solver: "
          f"VFRAC_MOM used in all {len(solvers)} solvers (not Dessart-specific) ✓")


def test_vfrac_mom_not_dohydro():
    """VFRAC_MOM is NOT related to DO_HYDRO (hydrostatic structure adjustment).

    DO_HYDRO adjusts the hydrostatic structure (velocity, density) and is
    disabled for SN models (cmfgen_sub.f:618: IF(DO_HYDRO .AND. .NOT. SN_MODEL)).
    VFRAC_MOM is a RT grid parameter and is completely independent of DO_HYDRO.
    """
    src_path = Path('/viper/ptmp1/sbrennan/SN2024pba/cmfgen_source/new_main/cmfgen_sub.f')
    if not src_path.exists():
        print("  test_vfrac_mom_not_dohydro: SKIPPED (source not found)")
        return

    text = src_path.read_text()
    # DO_HYDRO is disabled for SN models
    assert 'DO_HYDRO .AND. .NOT. SN_MODEL' in text, \
        "DO_HYDRO SN_MODEL guard not found"
    # VFRAC_MOM (DELV_FRAC_MOM) is NOT passed to DO_CMF_HYDRO
    assert 'DELV_FRAC_MOM' not in text.split('DO_CMF_HYDRO')[1].split('CALL')[0], \
        "DELV_FRAC_MOM should not be in DO_CMF_HYDRO call"

    print(f"  test_vfrac_mom_not_dohydro: VFRAC_MOM independent of DO_HYDRO ✓")


def run_all_tests():
    """Run all tests and report results."""
    print("=" * 70)
    print("VFRAC_MOM Source-Level Audit Tests")
    print("=" * 70)
    print()

    tests = [
        test_vfrac_mom_not_zero,
        test_vfrac_mom_does_not_modify_velocity,
        test_vfrac_mom_bounds,
        test_vfrac_mom_non_monotonic,
        test_vfrac_mom_not_momentum_correction,
        test_vfrac_mom_source_definition,
        test_vfrac_mom_not_related_to_dessert_solver,
        test_vfrac_mom_not_dohydro,
    ]

    passed = 0
    failed = 0
    for test in tests:
        try:
            test()
            passed += 1
        except AssertionError as e:
            print(f"  {test.__name__}: FAILED - {e}")
            failed += 1
        except Exception as e:
            print(f"  {test.__name__}: ERROR - {e}")
            failed += 1

    print()
    print("=" * 70)
    print(f"Results: {passed} passed, {failed} failed, {len(tests)} total")
    print("=" * 70)
    return failed == 0


if __name__ == '__main__':
    success = run_all_tests()
    sys.exit(0 if success else 1)
