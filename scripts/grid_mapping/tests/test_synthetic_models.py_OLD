"""Integration tests using synthetic models.

Tests the full grid_mapping pipeline on synthetic models with known
features. These are end-to-end tests that exercise the full path:
source grid -> feature detection -> grid selection -> conservative
remap -> validation.
"""

import sys
import os
import unittest
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from grid_mapping import (
    detect_all_features, conservative_remap, GridSelector,
    PointBudgetPolicy, validate_target_grid, sigma_velocity,
    count_velocity_reversals, write_sn_hydro_data,
)
from grid_mapping.synthetic_models import (
    synthetic_monotonic, synthetic_single_reversal, synthetic_multi_reversal,
    synthetic_shock, synthetic_composition_interface, synthetic_ionization_front,
    synthetic_all_features, SYNTHETIC_MODELS,
)


class TestSyntheticMonotonic(unittest.TestCase):
    """T1: Synthetic monotonic model."""

    def test_homologous_sigma_zero(self):
        """Homologous expansion produces sigma=0 everywhere."""
        source = synthetic_monotonic(100)
        sigma = sigma_velocity(source.r_centre_cm, source.velocity_cms)
        np.testing.assert_allclose(sigma, 0.0, atol=1e-10)

    def test_no_reversals(self):
        """Monotonic model has no velocity reversals."""
        source = synthetic_monotonic(100)
        sigma = sigma_velocity(source.r_centre_cm, source.velocity_cms)
        self.assertEqual(count_velocity_reversals(sigma), 0)

    def test_full_pipeline(self):
        """Full pipeline runs on monotonic model."""
        source = synthetic_monotonic(100)
        selector = GridSelector(source, target_points=30, max_points=40)
        target, report = selector.select_grid()
        self.assertTrue(target.n_points > 0)
        self.assertLess(target.n_points, source.n_shells)

        val = validate_target_grid(target, source, selector.feature_report)
        self.assertTrue(val.passed, f"Validation failed: {val.errors}")


class TestSyntheticSingleReversal(unittest.TestCase):
    """T2: Synthetic single-reversal model."""

    def test_reversal_detected(self):
        """Single reversal is detected."""
        source = synthetic_single_reversal(100)
        features = detect_all_features(source)
        self.assertGreaterEqual(features.n_velocity_minima, 1)

    def test_reversal_preserved_after_remap(self):
        """Reversal is preserved after conservative remapping."""
        source = synthetic_single_reversal(100)
        r_target = np.linspace(10e10, 1e10, 31)  # 30 target shells
        target = conservative_remap(source, r_target, source.time_days)

        source_sigma = sigma_velocity(source.r_centre_cm, source.velocity_cms)
        target_sigma = sigma_velocity(target.r_cm, target.velocity_cms)

        source_rev = count_velocity_reversals(source_sigma)
        target_rev = count_velocity_reversals(target_sigma)
        self.assertGreaterEqual(target_rev, 1,
                               "Reversal lost after remapping")

    def test_full_pipeline_preserves_reversal(self):
        """Full pipeline preserves the velocity reversal."""
        source = synthetic_single_reversal(100)
        selector = GridSelector(source, target_points=40, max_points=50)
        target, report = selector.select_grid()

        target_sigma = sigma_velocity(target.r_cm, target.velocity_cms)
        target_rev = count_velocity_reversals(target_sigma)
        self.assertGreaterEqual(target_rev, 1,
                               "Reversal lost in full pipeline")


class TestSyntheticMultiReversal(unittest.TestCase):
    """T3: Synthetic multi-reversal model."""

    def test_multiple_reversals_detected(self):
        """Multiple reversals are detected."""
        source = synthetic_multi_reversal(100)
        features = detect_all_features(source)
        self.assertGreaterEqual(features.n_velocity_minima, 2)

    def test_full_pipeline(self):
        """Full pipeline runs on multi-reversal model."""
        source = synthetic_multi_reversal(100)
        selector = GridSelector(source, target_points=50, max_points=60)
        target, report = selector.select_grid()
        self.assertTrue(target.n_points > 0)


class TestSyntheticShock(unittest.TestCase):
    """Shock preservation test."""

    def test_shock_detected(self):
        """Shock is detected."""
        source = synthetic_shock(100)
        features = detect_all_features(source)
        self.assertGreater(features.n_shocks + features.n_contact_discontinuities, 0)

    def test_mass_conservation_with_shock(self):
        """Mass is conserved even with a shock."""
        source = synthetic_shock(100)
        r_target = np.linspace(10e10, 1e10, 31)
        target = conservative_remap(source, r_target, source.time_days)

        source_mass = source.shell_mass_g.sum()
        target_mass = target.shell_mass_g.sum()
        residual = abs(target_mass - source_mass) / source_mass
        self.assertLess(residual, 1e-10)


class TestSyntheticCompositionInterface(unittest.TestCase):
    """Composition interface preservation test."""

    def test_interface_detected(self):
        """Composition interface is detected."""
        source = synthetic_composition_interface(100)
        features = detect_all_features(source)
        self.assertGreater(features.n_composition_interfaces, 0)

    def test_elemental_mass_conserved(self):
        """Elemental mass is conserved across the interface."""
        source = synthetic_composition_interface(100)
        r_target = np.linspace(10e10, 1e10, 31)
        target = conservative_remap(source, r_target, source.time_days)

        for elem in ['H', 'He', 'O']:
            source_mass = (source.mass_fractions[elem] * source.shell_mass_g).sum()
            target_mass = (target.mass_fractions[elem] * target.shell_mass_g).sum()
            residual = abs(target_mass - source_mass) / source_mass
            self.assertLess(residual, 1e-10,
                           f"{elem} mass not conserved: {residual}")

    def test_ni56_conserved(self):
        """Ni56 mass is conserved across the interface."""
        source = synthetic_composition_interface(100)
        r_target = np.linspace(10e10, 1e10, 31)
        target = conservative_remap(source, r_target, source.time_days)

        source_ni56 = (source.ni56_initial * source.shell_mass_g).sum()
        target_ni56 = (target.ni56_initial * target.shell_mass_g).sum()
        residual = abs(target_ni56 - source_ni56) / source_ni56
        self.assertLess(residual, 1e-10)


class TestSyntheticIonizationFront(unittest.TestCase):
    """Ionization front preservation test."""

    def test_front_detected(self):
        """Ionization front is detected."""
        source = synthetic_ionization_front(100)
        features = detect_all_features(source)
        self.assertGreater(features.n_ionization_fronts, 0)


class TestSyntheticAllFeatures(unittest.TestCase):
    """Combined model with all features."""

    def test_all_features_detected(self):
        """All feature types are detected in the combined model."""
        source = synthetic_all_features(100)
        features = detect_all_features(source)
        self.assertGreater(features.n_velocity_minima, 0)
        self.assertGreater(features.n_shocks + features.n_contact_discontinuities, 0)
        self.assertGreater(features.n_composition_interfaces, 0)
        self.assertGreater(features.n_ionization_fronts, 0)

    def test_full_pipeline_all_features(self):
        """Full pipeline runs on the combined model."""
        source = synthetic_all_features(100)
        selector = GridSelector(source, target_points=50, max_points=60)
        target, report = selector.select_grid()
        self.assertTrue(target.n_points > 0)

        val = validate_target_grid(target, source, selector.feature_report)
        # May have warnings but should not have hard errors
        self.assertEqual(len(val.errors), 0,
                        f"Validation errors: {val.errors}")


class TestAllModelsRegistry(unittest.TestCase):
    """Test that all synthetic models can be created and processed."""

    def test_all_models_create(self):
        """All synthetic models can be created."""
        for name, factory in SYNTHETIC_MODELS.items():
            with self.subTest(model=name):
                grid = factory(50)
                self.assertEqual(grid.n_shells, 50)
                errors = grid.validate()
                self.assertEqual(len(errors), 0,
                               f"{name} has validation errors: {errors}")


if __name__ == '__main__':
    unittest.main(verbosity=2)
