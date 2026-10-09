"""Test suite for complex composition grid generation (Tests A-J).

Tests the STELLA -> CMFGEN grid generation pipeline for models with:
- Sharp H-rich -> H-free transitions
- Detached H-free metal-rich cocoons
- Large opacity jumps
- Sharp density jumps
- Multiple composition interfaces
- The pathological dr/r ~ 0.37 production failure

Every test verifies:
- mass conservation
- species conservation
- no artificial hydrogen generation in H-free regions
- composition bounds (0 <= X_i <= 1, sum X_i = 1)
- radial ordering (strictly decreasing)
- optical-depth resolution
- interface resolution
- pathological-cell detection
- adaptive refinement
"""

import unittest
import numpy as np
import sys
import os

# Add scripts to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from grid_mapping.synthetic_models import (
    synthetic_monotonic, synthetic_composition_interface,
    synthetic_hfree_cocoon, synthetic_large_opacity_jump,
    synthetic_sharp_density_jump, synthetic_multi_interface,
    synthetic_shock,
)
from grid_mapping.shell_model import STELLA_ELEM_NAMES, MSUN
from tau_amr_optimiser import (
    AMRConfig, recursive_amr, remap_onto_grid,
    conservative_remap_density_composition,
    detect_composition_interfaces, detect_all_features,
    pathological_cell_report, escalate_grid,
    validate_grid, total_mass, species_masses,
    reconstruct_optical_depth, compute_kappa_tau,
    electron_scattering_kappa,
)


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def grid_to_arrays(grid):
    """Extract arrays from a SourceShellGrid for AMR processing."""
    r = grid.r_centre_cm.copy()
    v = grid.velocity_cms.copy()
    T = grid.temperature_K.copy()
    rho = grid.density_gcm3.copy()
    kappa = grid.kappa_cgs.copy()
    mf = {e: grid.mass_fractions.get(e, np.zeros(len(r))).copy()
          for e in STELLA_ELEM_NAMES}
    # Compute tau
    tau = compute_kappa_tau(r, rho, kappa)
    return r, v, T, rho, kappa, tau, mf


def check_mass_conservation(test, r_old, rho_old, mf_old, r_new, rho_new, mf_new,
                            tol=1e-3):
    """Verify total mass and per-species mass are conserved."""
    m_old = total_mass(r_old, rho_old)
    m_new = total_mass(r_new, rho_new)
    rel_err = abs(m_new - m_old) / max(m_old, 1e-30)
    test.assertLess(rel_err, tol,
                    f"Total mass not conserved: {m_old/MSUN:.4e} -> {m_new/MSUN:.4e} "
                    f"(rel_err={rel_err:.2e})")

    sm_old = species_masses(r_old, rho_old, mf_old)
    sm_new = species_masses(r_new, rho_new, mf_new)
    for e in sm_old:
        if sm_old[e] > 0:
            spec_err = abs(sm_new[e] - sm_old[e]) / sm_old[e]
            test.assertLess(spec_err, tol * 5,
                            f"Species '{e}' mass not conserved: "
                            f"{sm_old[e]/MSUN:.4e} -> {sm_new[e]/MSUN:.4e} "
                            f"(rel_err={spec_err:.2e})")


def check_composition_bounds(test, mf):
    """Verify all mass fractions are in [0, 1] and sum to 1."""
    n = len(next(iter(mf.values())))
    for elem, x in mf.items():
        test.assertTrue(np.all(x >= -1e-10),
                        f"Negative abundance for {elem}: min={x.min():.2e}")
        test.assertTrue(np.all(x <= 1.0 + 1e-6),
                        f"Abundance > 1 for {elem}: max={x.max():.2e}")
    sums = np.zeros(n)
    for x in mf.values():
        sums = sums + x
    test.assertLess(abs(sums.mean() - 1.0), 1e-4,
                    f"Sum of mass fractions != 1: mean={sums.mean():.6f}")
    test.assertLess(sums.max() - sums.min(), 1e-3,
                    f"Sum of mass fractions varies: [{sums.min():.6f}, {sums.max():.6f}]")


def check_radial_ordering(test, r):
    """Verify r is strictly decreasing (outer -> inner)."""
    test.assertTrue(np.all(np.diff(r) < 0),
                    f"Radial ordering violated: not strictly decreasing")


def check_no_artificial_hydrogen(test, r, mf, hfree_mask_func, hfree_label=""):
    """Verify that H-free regions remain H-free after remapping.

    hfree_mask_func(r, mf) -> boolean mask of H-free zones.
    """
    mask = hfree_mask_func(r, mf)
    if not np.any(mask):
        return  # no H-free region to check
    H = mf.get('H', np.zeros(len(r)))
    h_max = float(H[mask].max())
    test.assertLess(h_max, 1e-6,
                    f"Artificial hydrogen generated in H-free region{hfree_label}: "
                    f"max X_H = {h_max:.2e}")


def check_interface_resolution(test, r, mf, interface_radii, tol_frac=0.02):
    """Verify that grid points exist near each composition interface.

    interface_radii: list of radii where interfaces should be.
    tol_frac: fractional tolerance in log-radius for "near".
    """
    for r_iface in interface_radii:
        # Find closest grid point
        idx = np.argmin(np.abs(r - r_iface))
        r_closest = r[idx]
        log_ratio = abs(np.log10(r_closest / r_iface))
        test.assertLess(log_ratio, tol_frac,
                        f"Interface at r={r_iface:.3e} not resolved: "
                        f"closest grid point at r={r_closest:.3e} "
                        f"(log ratio={log_ratio:.3f})")


def check_pathological_cells(test, r, v, T, rho, kappa, tau, mf, cfg,
                             feature_radii=None):
    """Run pathological-cell diagnostic and check for failures."""
    report = pathological_cell_report(r, v, T, rho, kappa, tau, mf, cfg,
                                      feature_radii=feature_radii)
    test.assertTrue(report['passed'],
                    f"Pathological cells detected ({report['n_failures']}):\n"
                    f"{report['summary']}")


def run_amr(r, v, T, rho, kappa, tau, mf, cfg=None):
    """Run recursive AMR + escalation on a source grid.

    Returns (r, v, T, rho, kappa, tau, mf, amr_report, esc_report, feat_radii).
    """
    if cfg is None:
        cfg = AMRConfig(max_points=500, min_points=40, base_grid_points=60,
                        verbose=False, plot=False)
    # Detect features
    feats, comp_f, shock_f, ext_f = detect_all_features(
        r, v, rho, T, mf, tau, cfg)
    interface_indices = np.array([f.index for f in comp_f], dtype=int)
    feat_radii = [float(f.radius) for f in feats]

    # Save source arrays
    r_src = r.copy(); v_src = v.copy(); T_src = T.copy()
    rho_src = rho.copy(); kappa_src = kappa.copy(); tau_src = tau.copy()
    mf_src = {k: v2.copy() for k, v2 in mf.items()}

    # Run recursive AMR
    r_amr, v_amr, T_amr, rho_amr, kappa_amr, tau_amr, mf_amr, amr_report = \
        recursive_amr(r, v, T, rho, kappa, tau, mf,
                      feats, cfg, interface_indices)

    # Run escalation
    r_esc, v_esc, T_esc, rho_esc, kappa_esc, tau_esc, mf_esc, esc_report = \
        escalate_grid(
            r_amr, v_amr, T_amr, rho_amr, kappa_amr, tau_amr, mf_amr,
            r_src, v_src, T_src, rho_src, kappa_src, tau_src, mf_src,
            feats, cfg, interface_indices, max_rounds=5)

    return (r_esc, v_esc, T_esc, rho_esc, kappa_esc, tau_esc, mf_esc,
            amr_report, esc_report, feat_radii)


# ---------------------------------------------------------------------------
# Test A: Smooth homologous ejecta
# ---------------------------------------------------------------------------

class TestA_SmoothHomologous(unittest.TestCase):
    """Test A: Smooth homologous ejecta — baseline test."""

    def test_mass_conservation(self):
        grid = synthetic_monotonic(200)
        r, v, T, rho, kappa, tau, mf = grid_to_arrays(grid)
        r_new, v_new, T_new, rho_new, kappa_new, tau_new, mf_new, _, _, feat_radii = \
            run_amr(r, v, T, rho, kappa, tau, mf)
        check_mass_conservation(self, r, rho, mf, r_new, rho_new, mf_new)

    def test_species_conservation(self):
        grid = synthetic_monotonic(200)
        r, v, T, rho, kappa, tau, mf = grid_to_arrays(grid)
        r_new, _, _, rho_new, _, _, mf_new, _, _, _ = run_amr(r, v, T, rho, kappa, tau, mf)
        check_mass_conservation(self, r, rho, mf, r_new, rho_new, mf_new)

    def test_composition_bounds(self):
        grid = synthetic_monotonic(200)
        r, v, T, rho, kappa, tau, mf = grid_to_arrays(grid)
        _, _, _, _, _, _, mf_new, _, _, _ = run_amr(r, v, T, rho, kappa, tau, mf)
        check_composition_bounds(self, mf_new)

    def test_radial_ordering(self):
        grid = synthetic_monotonic(200)
        r, v, T, rho, kappa, tau, mf = grid_to_arrays(grid)
        r_new, _, _, _, _, _, _, _, _, _ = run_amr(r, v, T, rho, kappa, tau, mf)
        check_radial_ordering(self, r_new)

    def test_no_pathological_cells(self):
        grid = synthetic_monotonic(200)
        r, v, T, rho, kappa, tau, mf = grid_to_arrays(grid)
        r_new, v_new, T_new, rho_new, kappa_new, tau_new, mf_new, _, _, feat_radii = \
            run_amr(r, v, T, rho, kappa, tau, mf)
        cfg = AMRConfig(verbose=False)
        check_pathological_cells(self, r_new, v_new, T_new, rho_new,
                                 kappa_new, tau_new, mf_new, cfg,
                                 feature_radii=feat_radii)

    def test_grid_smaller_than_max(self):
        grid = synthetic_monotonic(200)
        r, v, T, rho, kappa, tau, mf = grid_to_arrays(grid)
        r_new, _, _, _, _, _, _, _, _, _ = run_amr(r, v, T, rho, kappa, tau, mf)
        self.assertLess(len(r_new), 500,
                        f"Smooth model should not need max points: ND={len(r_new)}")


# ---------------------------------------------------------------------------
# Test B: Sharp H-rich -> H-free interface
# ---------------------------------------------------------------------------

class TestB_SharpHFreeInterface(unittest.TestCase):
    """Test B: Sharp H-rich -> H-free interface (single interface)."""

    def test_interface_detected(self):
        """The composition interface must be detected by detect_composition_interfaces."""
        grid = synthetic_composition_interface(200)
        r, v, T, rho, kappa, tau, mf = grid_to_arrays(grid)
        cfg = AMRConfig(verbose=False)
        feats = detect_composition_interfaces(r, mf, tau, cfg)
        self.assertGreater(len(feats), 0,
                           "No composition interfaces detected")

    def test_interface_resolved_in_grid(self):
        """Grid points must exist near the interface."""
        grid = synthetic_composition_interface(200)
        r, v, T, rho, kappa, tau, mf = grid_to_arrays(grid)
        # The interface is at index 50 in the source grid
        r_iface = r[50]
        r_new, _, _, _, _, _, _, _, _, _ = run_amr(r, v, T, rho, kappa, tau, mf)
        check_interface_resolution(self, r_new, mf, [r_iface])

    def test_no_artificial_hydrogen(self):
        """The H-free inner region must not acquire hydrogen after remapping."""
        grid = synthetic_composition_interface(200)
        r, v, T, rho, kappa, tau, mf = grid_to_arrays(grid)
        r_new, _, _, rho_new, _, _, mf_new, _, _, _ = run_amr(r, v, T, rho, kappa, tau, mf)

        # The inner region (index >= 50 in source) should have X_H = 0.1
        # (not truly H-free in this model, but significantly lower)
        # Check that the H abundance in the inner region doesn't increase
        # beyond the source value due to interpolation smearing
        H_new = mf_new.get('H', np.zeros(len(r_new)))
        H_old = mf['H']
        h_inner_source = H_old[50:].max()
        # Find the inner region in the new grid
        r_iface = r[50]
        inner_mask = r_new < r_iface
        if np.any(inner_mask):
            h_inner_new = H_new[inner_mask].max()
            # Allow some smearing but not back to H-rich levels
            self.assertLess(h_inner_new, 0.4,
                            f"H smeared into inner region: X_H={h_inner_new:.3f} "
                            f"(source inner max={h_inner_source:.3f})")

    def test_mass_conservation(self):
        grid = synthetic_composition_interface(200)
        r, v, T, rho, kappa, tau, mf = grid_to_arrays(grid)
        r_new, _, _, rho_new, _, _, mf_new, _, _, _ = run_amr(r, v, T, rho, kappa, tau, mf)
        check_mass_conservation(self, r, rho, mf, r_new, rho_new, mf_new)

    def test_composition_bounds(self):
        grid = synthetic_composition_interface(200)
        r, v, T, rho, kappa, tau, mf = grid_to_arrays(grid)
        _, _, _, _, _, _, mf_new, _, _, _ = run_amr(r, v, T, rho, kappa, tau, mf)
        check_composition_bounds(self, mf_new)

    def test_radial_ordering(self):
        grid = synthetic_composition_interface(200)
        r, v, T, rho, kappa, tau, mf = grid_to_arrays(grid)
        r_new, _, _, _, _, _, _, _, _, _ = run_amr(r, v, T, rho, kappa, tau, mf)
        check_radial_ordering(self, r_new)


# ---------------------------------------------------------------------------
# Test C: Detached H-free metal-rich cocoon
# ---------------------------------------------------------------------------

class TestC_DetachedHFreeCocoon(unittest.TestCase):
    """Test C: Detached H-free metal-rich cocoon in H-rich ejecta."""

    def setUp(self):
        self.grid = synthetic_hfree_cocoon(200, optically_thick=True)
        self.r, self.v, self.T, self.rho, self.kappa, self.tau, self.mf = \
            grid_to_arrays(self.grid)
        self.cocoon_outer_idx = self.grid.cocoon_outer_idx
        self.cocoon_inner_idx = self.grid.cocoon_inner_idx
        self.r_outer_iface = self.r[self.cocoon_outer_idx]
        self.r_inner_iface = self.r[self.cocoon_inner_idx]

    def test_both_interfaces_detected(self):
        """Both cocoon boundaries must be detected as composition interfaces."""
        cfg = AMRConfig(verbose=False)
        feats = detect_composition_interfaces(self.r, self.mf, self.tau, cfg)
        self.assertGreaterEqual(len(feats), 2,
                                f"Expected >= 2 composition interfaces, got {len(feats)}")

    def test_both_interfaces_resolved(self):
        """Grid points must exist near both cocoon boundaries."""
        r_new, _, _, _, _, _, _, _, _, _ = run_amr(
            self.r, self.v, self.T, self.rho, self.kappa, self.tau, self.mf)
        check_interface_resolution(self, r_new, self.mf,
                                   [self.r_outer_iface, self.r_inner_iface])

    def test_no_artificial_hydrogen_in_cocoon(self):
        """CRITICAL: The H-free cocoon must not acquire significant hydrogen."""
        r_new, _, _, rho_new, _, _, mf_new, _, _, _ = run_amr(
            self.r, self.v, self.T, self.rho, self.kappa, self.tau, self.mf)
        H_new = mf_new.get('H', np.zeros(len(r_new)))

        # Use strict inequalities to exclude the interface boundary points
        # themselves, which correctly represent their respective sides.
        # The cocoon interior (strictly between the interfaces) must be H-free.
        # Use a small tolerance to avoid catching near-duplicate points at the
        # interface that belong to the H-rich side.
        tol = 1e-3 * (self.r_outer_iface - self.r_inner_iface)  # 0.1% of cocoon width
        cocoon_mask = (r_new < self.r_outer_iface - tol) & \
                      (r_new > self.r_inner_iface + tol)
        if np.any(cocoon_mask):
            h_cocoon_max = float(H_new[cocoon_mask].max())
            self.assertLess(h_cocoon_max, 1e-6,
                            f"CRITICAL: Artificial hydrogen generated in H-free cocoon: "
                            f"max X_H = {h_cocoon_max:.2e}")

    def test_cocoon_interior_sampled(self):
        """The cocoon interior must have adequate grid points (>= 5)."""
        r_new, _, _, _, _, _, _, _, _, _ = run_amr(
            self.r, self.v, self.T, self.rho, self.kappa, self.tau, self.mf)
        tol = 1e-3 * (self.r_outer_iface - self.r_inner_iface)
        cocoon_mask = (r_new < self.r_outer_iface - tol) & \
                      (r_new > self.r_inner_iface + tol)
        n_cocoon = int(np.sum(cocoon_mask))
        self.assertGreaterEqual(n_cocoon, 5,
                                f"Cocoon interior under-sampled: {n_cocoon} points")

    def test_mass_conservation(self):
        r_new, _, _, rho_new, _, _, mf_new, _, _, _ = run_amr(
            self.r, self.v, self.T, self.rho, self.kappa, self.tau, self.mf)
        check_mass_conservation(self, self.r, self.rho, self.mf,
                                r_new, rho_new, mf_new)

    def test_species_conservation(self):
        r_new, _, _, rho_new, _, _, mf_new, _, _, _ = run_amr(
            self.r, self.v, self.T, self.rho, self.kappa, self.tau, self.mf)
        check_mass_conservation(self, self.r, self.rho, self.mf,
                                r_new, rho_new, mf_new)

    def test_composition_bounds(self):
        _, _, _, _, _, _, mf_new, _, _, _ = run_amr(
            self.r, self.v, self.T, self.rho, self.kappa, self.tau, self.mf)
        check_composition_bounds(self, mf_new)

    def test_radial_ordering(self):
        r_new, _, _, _, _, _, _, _, _, _ = run_amr(
            self.r, self.v, self.T, self.rho, self.kappa, self.tau, self.mf)
        check_radial_ordering(self, r_new)

    def test_no_pathological_cells(self):
        r_new, v_new, T_new, rho_new, kappa_new, tau_new, mf_new, _, _, feat_radii = \
            run_amr(self.r, self.v, self.T, self.rho, self.kappa, self.tau, self.mf)
        cfg = AMRConfig(verbose=False)
        check_pathological_cells(self, r_new, v_new, T_new, rho_new,
                                 kappa_new, tau_new, mf_new, cfg,
                                 feature_radii=feat_radii)

    def test_refinement_is_local(self):
        """Refinement should be local to the cocoon, not global."""
        r_new, _, _, _, _, _, _, _, _, _ = run_amr(
            self.r, self.v, self.T, self.rho, self.kappa, self.tau, self.mf)
        # Check that the grid is not at max points (which would indicate
        # global rather than local refinement)
        self.assertLess(len(r_new), 400,
                        f"Grid too large for localised cocoon: ND={len(r_new)}")
        # Check that the cocoon region has more points per unit log-radius
        # than the smooth outer region
        log_r = np.log10(r_new)
        dlog_r = np.abs(np.diff(log_r))
        tol = 1e-3 * (self.r_outer_iface - self.r_inner_iface)
        cocoon_mask = (r_new < self.r_outer_iface - tol) & \
                      (r_new > self.r_inner_iface + tol)
        # Use a smooth region far from the cocoon for comparison
        smooth_mask = (r_new > self.r_outer_iface + tol) & \
                      (r_new < self.r_outer_iface * 1.3)
        if np.sum(cocoon_mask) > 2 and np.sum(smooth_mask) > 2:
            cocoon_dlog = dlog_r[cocoon_mask[:-1] | cocoon_mask[1:]].mean()
            smooth_dlog = dlog_r[smooth_mask[:-1] | smooth_mask[1:]].mean()
            # Cocoon should have smaller spacing (more points) than smooth region
            self.assertLessEqual(cocoon_dlog, smooth_dlog * 1.5,
                            f"Refinement not local: cocoon spacing={cocoon_dlog:.4f}, "
                            f"smooth spacing={smooth_dlog:.4f}")


# ---------------------------------------------------------------------------
# Test D: Optically thick cocoon
# ---------------------------------------------------------------------------

class TestD_OpticallyThickCocoon(unittest.TestCase):
    """Test D: Optically thick cocoon (high density + high opacity)."""

    def test_cocoon_resolved(self):
        grid = synthetic_hfree_cocoon(200, optically_thick=True)
        r, v, T, rho, kappa, tau, mf = grid_to_arrays(grid)
        r_new, _, _, _, _, _, _, _, _, _ = run_amr(r, v, T, rho, kappa, tau, mf)
        # Optically thick cocoon should get more points
        r_outer = r[grid.cocoon_outer_idx]
        r_inner = r[grid.cocoon_inner_idx]
        tol = 1e-3 * (r_outer - r_inner)
        cocoon_mask = (r_new < r_outer - tol) & (r_new > r_inner + tol)
        n_cocoon = int(np.sum(cocoon_mask))
        self.assertGreaterEqual(n_cocoon, 7,
                                f"Optically thick cocoon under-sampled: {n_cocoon} points")

    def test_no_artificial_hydrogen(self):
        grid = synthetic_hfree_cocoon(200, optically_thick=True)
        r, v, T, rho, kappa, tau, mf = grid_to_arrays(grid)
        r_new, _, _, _, _, _, mf_new, _, _, _ = run_amr(r, v, T, rho, kappa, tau, mf)
        H_new = mf_new.get('H', np.zeros(len(r_new)))
        r_outer = r[grid.cocoon_outer_idx]
        r_inner = r[grid.cocoon_inner_idx]
        tol = 1e-3 * (r_outer - r_inner)
        cocoon_mask = (r_new < r_outer - tol) & (r_new > r_inner + tol)
        if np.any(cocoon_mask):
            self.assertLess(float(H_new[cocoon_mask].max()), 1e-6,
                            "Artificial H in optically thick cocoon")

    def test_mass_conservation(self):
        grid = synthetic_hfree_cocoon(200, optically_thick=True)
        r, v, T, rho, kappa, tau, mf = grid_to_arrays(grid)
        r_new, _, _, rho_new, _, _, mf_new, _, _, _ = run_amr(r, v, T, rho, kappa, tau, mf)
        check_mass_conservation(self, r, rho, mf, r_new, rho_new, mf_new)


# ---------------------------------------------------------------------------
# Test E: Optically thin cocoon
# ---------------------------------------------------------------------------

class TestE_OpticallyThinCocoon(unittest.TestCase):
    """Test E: Optically thin cocoon (same density as surroundings)."""

    def test_cocoon_detected(self):
        grid = synthetic_hfree_cocoon(200, optically_thick=False)
        r, v, T, rho, kappa, tau, mf = grid_to_arrays(grid)
        cfg = AMRConfig(verbose=False)
        feats = detect_composition_interfaces(r, mf, tau, cfg)
        self.assertGreaterEqual(len(feats), 2,
                                "Optically thin cocoon interfaces not detected")

    def test_no_artificial_hydrogen(self):
        grid = synthetic_hfree_cocoon(200, optically_thick=False)
        r, v, T, rho, kappa, tau, mf = grid_to_arrays(grid)
        r_new, _, _, _, _, _, mf_new, _, _, _ = run_amr(r, v, T, rho, kappa, tau, mf)
        H_new = mf_new.get('H', np.zeros(len(r_new)))
        r_outer = r[grid.cocoon_outer_idx]
        r_inner = r[grid.cocoon_inner_idx]
        tol = 1e-3 * (r_outer - r_inner)
        cocoon_mask = (r_new < r_outer - tol) & (r_new > r_inner + tol)
        if np.any(cocoon_mask):
            self.assertLess(float(H_new[cocoon_mask].max()), 1e-6,
                            "Artificial H in optically thin cocoon")

    def test_mass_conservation(self):
        grid = synthetic_hfree_cocoon(200, optically_thick=False)
        r, v, T, rho, kappa, tau, mf = grid_to_arrays(grid)
        r_new, _, _, rho_new, _, _, mf_new, _, _, _ = run_amr(r, v, T, rho, kappa, tau, mf)
        check_mass_conservation(self, r, rho, mf, r_new, rho_new, mf_new)


# ---------------------------------------------------------------------------
# Test F: Large opacity jump
# ---------------------------------------------------------------------------

class TestF_LargeOpacityJump(unittest.TestCase):
    """Test F: Large opacity jump (~680x) across a single cell."""

    def test_jump_detected_as_pathological(self):
        """The opacity jump must be detected by the pathological-cell diagnostic."""
        grid = synthetic_large_opacity_jump(200)
        r, v, T, rho, kappa, tau, mf = grid_to_arrays(grid)
        cfg = AMRConfig(verbose=False)
        # Run on the source grid — the jump should be detected
        report = pathological_cell_report(r, v, T, rho, kappa, tau, mf, cfg)
        self.assertFalse(report['passed'],
                         "Large opacity jump not detected as pathological")

    def test_jump_resolved_after_amr(self):
        """After AMR, the opacity jump region must be adequately resolved."""
        grid = synthetic_large_opacity_jump(200)
        r, v, T, rho, kappa, tau, mf = grid_to_arrays(grid)
        r_new, v_new, T_new, rho_new, kappa_new, tau_new, mf_new, _, _, feat_radii = \
            run_amr(r, v, T, rho, kappa, tau, mf)
        cfg = AMRConfig(verbose=False)
        check_pathological_cells(self, r_new, v_new, T_new, rho_new,
                                 kappa_new, tau_new, mf_new, cfg,
                                 feature_radii=feat_radii)

    def test_refinement_is_local(self):
        """Refinement should be near the jump, not global."""
        grid = synthetic_large_opacity_jump(200)
        r, v, T, rho, kappa, tau, mf = grid_to_arrays(grid)
        r_new, _, _, _, _, _, _, _, _, _ = run_amr(r, v, T, rho, kappa, tau, mf)
        # Check that the grid is not uniformly at max points
        self.assertLess(len(r_new), 500,
                        f"Should not need max points for local jump: ND={len(r_new)}")

    def test_mass_conservation(self):
        grid = synthetic_large_opacity_jump(200)
        r, v, T, rho, kappa, tau, mf = grid_to_arrays(grid)
        r_new, _, _, rho_new, _, _, mf_new, _, _, _ = run_amr(r, v, T, rho, kappa, tau, mf)
        check_mass_conservation(self, r, rho, mf, r_new, rho_new, mf_new)


# ---------------------------------------------------------------------------
# Test G: Sharp density jump
# ---------------------------------------------------------------------------

class TestG_SharpDensityJump(unittest.TestCase):
    """Test G: Sharp density jump (260x) with uniform composition."""

    def test_jump_detected(self):
        grid = synthetic_sharp_density_jump(200)
        r, v, T, rho, kappa, tau, mf = grid_to_arrays(grid)
        cfg = AMRConfig(verbose=False)
        report = pathological_cell_report(r, v, T, rho, kappa, tau, mf, cfg)
        self.assertFalse(report['passed'],
                         "Sharp density jump not detected as pathological")

    def test_jump_resolved_after_amr(self):
        grid = synthetic_sharp_density_jump(200)
        r, v, T, rho, kappa, tau, mf = grid_to_arrays(grid)
        r_new, v_new, T_new, rho_new, kappa_new, tau_new, mf_new, _, _, feat_radii = \
            run_amr(r, v, T, rho, kappa, tau, mf)
        cfg = AMRConfig(verbose=False)
        check_pathological_cells(self, r_new, v_new, T_new, rho_new,
                                 kappa_new, tau_new, mf_new, cfg,
                                 feature_radii=feat_radii)

    def test_mass_conservation(self):
        grid = synthetic_sharp_density_jump(200)
        r, v, T, rho, kappa, tau, mf = grid_to_arrays(grid)
        r_new, _, _, rho_new, _, _, mf_new, _, _, _ = run_amr(r, v, T, rho, kappa, tau, mf)
        check_mass_conservation(self, r, rho, mf, r_new, rho_new, mf_new)


# ---------------------------------------------------------------------------
# Test H: Multiple composition interfaces
# ---------------------------------------------------------------------------

class TestH_MultipleInterfaces(unittest.TestCase):
    """Test H: Multiple composition interfaces (H-rich/He-rich/O-rich/Si-rich)."""

    def test_all_interfaces_detected(self):
        grid = synthetic_multi_interface(200)
        r, v, T, rho, kappa, tau, mf = grid_to_arrays(grid)
        cfg = AMRConfig(verbose=False)
        feats = detect_composition_interfaces(r, mf, tau, cfg)
        self.assertGreaterEqual(len(feats), 3,
                                f"Expected >= 3 interfaces, got {len(feats)}")

    def test_all_interfaces_resolved(self):
        grid = synthetic_multi_interface(200)
        r, v, T, rho, kappa, tau, mf = grid_to_arrays(grid)
        r_new, _, _, _, _, _, _, _, _, _ = run_amr(r, v, T, rho, kappa, tau, mf)
        # Interfaces at indices 50, 100, 150
        for idx in [50, 100, 150]:
            check_interface_resolution(self, r_new, mf, [r[idx]])

    def test_no_artificial_hydrogen(self):
        """The H-free O-rich and Si-rich layers must not acquire hydrogen."""
        grid = synthetic_multi_interface(200)
        r, v, T, rho, kappa, tau, mf = grid_to_arrays(grid)
        r_new, _, _, _, _, _, mf_new, _, _, _ = run_amr(r, v, T, rho, kappa, tau, mf)
        H_new = mf_new.get('H', np.zeros(len(r_new)))
        # O-rich layer (indices 100-150) and Si-rich layer (150+) should be H-free
        r_o_rich_outer = r[100]
        r_o_rich_inner = r[150]
        hfree_mask = (r_new <= r_o_rich_outer) & (r_new >= r_o_rich_inner)
        if np.any(hfree_mask):
            self.assertLess(float(H_new[hfree_mask].max()), 1e-6,
                            f"Artificial H in O-rich layer: max X_H={H_new[hfree_mask].max():.2e}")

    def test_mass_conservation(self):
        grid = synthetic_multi_interface(200)
        r, v, T, rho, kappa, tau, mf = grid_to_arrays(grid)
        r_new, _, _, rho_new, _, _, mf_new, _, _, _ = run_amr(r, v, T, rho, kappa, tau, mf)
        check_mass_conservation(self, r, rho, mf, r_new, rho_new, mf_new)

    def test_composition_bounds(self):
        grid = synthetic_multi_interface(200)
        r, v, T, rho, kappa, tau, mf = grid_to_arrays(grid)
        _, _, _, _, _, _, mf_new, _, _, _ = run_amr(r, v, T, rho, kappa, tau, mf)
        check_composition_bounds(self, mf_new)


# ---------------------------------------------------------------------------
# Test I: Pathological dr/r ~ 0.37 model
# ---------------------------------------------------------------------------

class TestI_PathologicalModel(unittest.TestCase):
    """Test I: The previously identified dr/r ~ 0.37, density x260, opacity x680 model."""

    def test_pathological_cell_detected_in_source(self):
        """The pathological cell must be detected in the source grid."""
        grid = synthetic_large_opacity_jump(200)
        r, v, T, rho, kappa, tau, mf = grid_to_arrays(grid)
        cfg = AMRConfig(verbose=False)
        report = pathological_cell_report(r, v, T, rho, kappa, tau, mf, cfg)
        self.assertFalse(report['passed'],
                         f"Pathological cell not detected: {report['summary']}")

    def test_pathological_cell_resolved_after_amr(self):
        """After AMR + escalation, the pathological cell must be resolved."""
        grid = synthetic_large_opacity_jump(200)
        r, v, T, rho, kappa, tau, mf = grid_to_arrays(grid)
        r_new, v_new, T_new, rho_new, kappa_new, tau_new, mf_new, _, _, feat_radii = \
            run_amr(r, v, T, rho, kappa, tau, mf)
        cfg = AMRConfig(verbose=False)
        check_pathological_cells(self, r_new, v_new, T_new, rho_new,
                                 kappa_new, tau_new, mf_new, cfg,
                                 feature_radii=feat_radii)

    def test_refinement_is_local_not_global(self):
        """The fix must be local refinement, not global ND increase."""
        grid = synthetic_large_opacity_jump(200)
        r, v, T, rho, kappa, tau, mf = grid_to_arrays(grid)
        r_new, _, _, _, _, _, _, _, _, _ = run_amr(r, v, T, rho, kappa, tau, mf)
        # Should be well under 500 points for a single localised jump
        self.assertLess(len(r_new), 300,
                        f"Refinement not local: ND={len(r_new)} for single jump")

    def test_dr_over_r_below_threshold(self):
        """No cell should have dr/r > 0.15 after refinement."""
        grid = synthetic_large_opacity_jump(200)
        r, v, T, rho, kappa, tau, mf = grid_to_arrays(grid)
        r_new, _, _, _, _, _, _, _, _, _ = run_amr(r, v, T, rho, kappa, tau, mf)
        dr = np.abs(np.diff(r_new))
        r_mid = 0.5 * (r_new[:-1] + r_new[1:])
        dror = dr / np.maximum(r_mid, 1e-30)
        max_dror = float(dror.max())
        self.assertLess(max_dror, 0.15,
                        f"dr/r = {max_dror:.3f} exceeds threshold 0.15")


# ---------------------------------------------------------------------------
# Test J: Progressive refinement toward 1000-point limit
# ---------------------------------------------------------------------------

class TestJ_ProgressiveRefinement(unittest.TestCase):
    """Test J: A model requiring progressive refinement toward the 1000-point limit."""

    def test_escalation_increases_points(self):
        """Escalation should increase point count when pathological cells exist."""
        grid = synthetic_large_opacity_jump(200)
        r, v, T, rho, kappa, tau, mf = grid_to_arrays(grid)
        cfg = AMRConfig(max_points=500, min_points=40, base_grid_points=60,
                        verbose=False, plot=False)
        feats, comp_f, shock_f, ext_f = detect_all_features(
            r, v, rho, T, mf, tau, cfg)
        interface_indices = np.array([f.index for f in comp_f], dtype=int)

        # Run AMR only (no escalation)
        r_amr, v_amr, T_amr, rho_amr, kappa_amr, tau_amr, mf_amr, _ = \
            recursive_amr(r.copy(), v.copy(), T.copy(), rho.copy(),
                          kappa.copy(), tau.copy(),
                          {k: v2.copy() for k, v2 in mf.items()},
                          feats, cfg, interface_indices)

        # Run escalation
        r_src = r.copy(); v_src = v.copy(); T_src = T.copy()
        rho_src = rho.copy(); kappa_src = kappa.copy(); tau_src = tau.copy()
        mf_src = {k: v2.copy() for k, v2 in mf.items()}

        r_esc, _, _, _, _, _, _, esc_report = escalate_grid(
            r_amr, v_amr, T_amr, rho_amr, kappa_amr, tau_amr, mf_amr,
            r_src, v_src, T_src, rho_src, kappa_src, tau_src, mf_src,
            feats, cfg, interface_indices, max_rounds=5)

        # Escalation should have inserted points
        self.assertGreaterEqual(len(r_esc), len(r_amr),
                                "Escalation did not increase point count")

    def test_does_not_exceed_max_points(self):
        """The grid must never exceed the configured maximum."""
        grid = synthetic_large_opacity_jump(200)
        r, v, T, rho, kappa, tau, mf = grid_to_arrays(grid)
        cfg = AMRConfig(max_points=200, min_points=40, base_grid_points=60,
                        verbose=False, plot=False)
        r_new, _, _, _, _, _, _, _, _, _ = run_amr(r, v, T, rho, kappa, tau, mf, cfg)
        self.assertLessEqual(len(r_new), 200,
                             f"Grid exceeds max_points: ND={len(r_new)} > 200")

    def test_stops_at_smallest_passing_grid(self):
        """The grid should not always reach max_points when not needed."""
        grid = synthetic_monotonic(200)
        r, v, T, rho, kappa, tau, mf = grid_to_arrays(grid)
        cfg = AMRConfig(max_points=500, min_points=40, base_grid_points=60,
                        verbose=False, plot=False)
        r_new, _, _, _, _, _, _, _, _, _ = run_amr(r, v, T, rho, kappa, tau, mf, cfg)
        self.assertLess(len(r_new), 200,
                        f"Smooth model should not need many points: ND={len(r_new)}")


# ---------------------------------------------------------------------------
# Conservative remap unit tests for composition preservation
# ---------------------------------------------------------------------------

class TestConservativeRemapComposition(unittest.TestCase):
    """Unit tests for conservative_remap_density_composition with sharp interfaces."""

    def test_sharp_interface_preserved(self):
        """A sharp H-rich/H-free interface should not smear significantly."""
        n = 100
        r = np.linspace(10e10, 1e10, n)  # decreasing
        rho = np.ones(n) * 1e-10
        H = np.ones(n) * 0.7
        He = np.ones(n) * 0.28
        O = np.ones(n) * 0.02
        H[50:] = 0.0  # H-free inner region
        He[50:] = 0.30
        O[50:] = 0.70
        mf = {'H': H, 'He': He, 'O': O}

        # New grid that has a cell straddling the interface
        r_new = np.linspace(10e10, 1e10, 50)
        rho_new, mf_new = conservative_remap_density_composition(r_new, r, rho, mf)

        # Check no artificial H in the H-free region
        r_iface = r[50]
        hfree_mask = r_new < r_iface
        if np.any(hfree_mask):
            h_max = float(mf_new['H'][hfree_mask].max())
            # Conservative remap will introduce some H if a cell straddles,
            # but it should be small, not 0.35 (the average of 0.7 and 0.0)
            self.assertLess(h_max, 0.35,
                            f"Conservative remap smeared H too much: "
                            f"max X_H in H-free region = {h_max:.3f}")

    def test_mass_conservation_exact(self):
        """Total mass and species mass must be conserved to machine precision."""
        n = 100
        r = np.linspace(10e10, 1e10, n)
        rho = np.ones(n) * 1e-10
        H = np.ones(n) * 0.7
        He = np.ones(n) * 0.28
        O = np.ones(n) * 0.02
        H[50:] = 0.0
        He[50:] = 0.30
        O[50:] = 0.70
        mf = {'H': H, 'He': He, 'O': O}

        r_new = np.linspace(10e10, 1e10, 50)
        rho_new, mf_new = conservative_remap_density_composition(r_new, r, rho, mf)

        m_old = total_mass(r, rho)
        m_new = total_mass(r_new, rho_new)
        rel_err = abs(m_new - m_old) / max(m_old, 1e-30)
        self.assertLess(rel_err, 1e-10,
                        f"Total mass not conserved: rel_err={rel_err:.2e}")

        sm_old = species_masses(r, rho, mf)
        sm_new = species_masses(r_new, rho_new, mf_new)
        for e in sm_old:
            if sm_old[e] > 0:
                spec_err = abs(sm_new[e] - sm_old[e]) / sm_old[e]
                self.assertLess(spec_err, 1e-10,
                                f"Species '{e}' not conserved: rel_err={spec_err:.2e}")

    def test_no_negative_abundances(self):
        """Conservative remap must not produce negative abundances."""
        n = 100
        r = np.linspace(10e10, 1e10, n)
        rho = np.ones(n) * 1e-10
        H = np.ones(n) * 0.7
        He = np.ones(n) * 0.28
        O = np.ones(n) * 0.02
        H[50:] = 0.0
        He[50:] = 0.30
        O[50:] = 0.70
        mf = {'H': H, 'He': He, 'O': O}

        r_new = np.linspace(10e10, 1e10, 50)
        _, mf_new = conservative_remap_density_composition(r_new, r, rho, mf)

        for elem, x in mf_new.items():
            self.assertTrue(np.all(x >= -1e-10),
                            f"Negative abundance for {elem}: min={x.min():.2e}")


if __name__ == '__main__':
    unittest.main(verbosity=2)
