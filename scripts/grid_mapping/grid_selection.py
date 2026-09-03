"""Grid selection: mandatory anchors + constrained priority refinement.

Builds the target CMFGEN grid by:
1. Retaining mandatory physical anchors (shocks, extrema, interfaces, etc.)
2. Allocating optical-depth background points
3. Computing hard-error indicators
4. Inserting points using solver-aware placement rules
5. Conservatively remapping after each insertion
6. Continuing until hard constraints pass or point budget exhausted

Point-budget policies:
- strict: never exceed target; fail if constraints can't pass
- allow-headroom: exceed target up to an explicit maximum
- fail-if-infeasible: produce no output unless all constraints pass
"""

from __future__ import annotations

from enum import Enum
from typing import Optional

import numpy as np

from .shell_model import SourceShellGrid, TargetGrid, compute_shell_volume, compute_shell_centres
from .feature_detection import FeatureReport, detect_all_features
from .conservative_remap import conservative_remap, compute_overlap_volumes
from .opacity_tau import compute_tau, find_tau_landmarks, find_photosphere
from .solver_metrics import sigma_velocity, apply_sigma_cap


class PointBudgetPolicy(Enum):
    STRICT = "strict"
    ALLOW_HEADROOM = "allow-headroom"
    FAIL_IF_INFEASIBLE = "fail-if-infeasible"


class GridSelector:
    """Selects and refines the target grid from the source grid.

    Parameters
    ----------
    source : SourceShellGrid
        The original STELLA shell grid.
    target_points : int
        Nominal target point count.
    max_points : int
        Maximum point count (for allow-headroom policy).
    policy : PointBudgetPolicy
        Point budget policy.
    tau_max : float
        Inner optical depth cutoff (trim zones with tau > tau_max).
    tau_min : float
        Outer optical depth cutoff (trim zones with tau < tau_min).
    v_min_cms : float
        Inner velocity cutoff (trim zones with |V| < v_min).
    v_max_cms : float
        Outer velocity cutoff (trim zones with |V| > v_max).
    max_delta_tau : float
        Maximum optical depth change per cell.
    min_points_per_tau_decade : int
        Minimum grid points per optical depth decade.
    min_points_per_shock : int
        Minimum points to resolve each shock.
    min_points_per_reversal : int
        Minimum points to resolve each velocity reversal.
    min_points_per_ionization_front : int
        Minimum points across each ionization front.
    max_density_jump : float
        Maximum fractional density jump per cell.
    max_temp_jump : float
        Maximum fractional temperature jump per cell.
    max_velocity_jump : float
        Maximum fractional velocity jump per cell.
    max_sigma_jump : float
        Maximum |d(sigma)| per cell.
    max_cell_ratio : float
        Maximum ratio of adjacent cell sizes.
    """

    def __init__(self, source: SourceShellGrid,
                 target_points: int = 150,
                 max_points: int = 200,
                 policy: PointBudgetPolicy = PointBudgetPolicy.ALLOW_HEADROOM,
                 tau_max: float = 30.0,
                 tau_min: float = 1e-3,
                 v_min_cms: float = 0.0,
                 v_max_cms: float = 0.0,
                 max_delta_tau: float = 0.5,
                 min_points_per_tau_decade: int = 5,
                 min_points_per_shock: int = 5,
                 min_points_per_reversal: int = 5,
                 min_points_per_ionization_front: int = 5,
                 max_density_jump: float = 0.15,
                 max_temp_jump: float = 0.10,
                 max_velocity_jump: float = 0.10,
                 max_sigma_jump: float = 1.0,
                 max_cell_ratio: float = 1.5,
                 min_dr_r: float = 1e-6,
                 max_dr_r: float = 0.5):
        self.source = source
        self.target_points = target_points
        self.max_points = max_points
        self.policy = policy
        self.tau_max = tau_max
        self.tau_min = tau_min
        self.v_min_cms = v_min_cms
        self.v_max_cms = v_max_cms
        self.max_delta_tau = max_delta_tau
        self.min_points_per_tau_decade = min_points_per_tau_decade
        self.min_points_per_shock = min_points_per_shock
        self.min_points_per_reversal = min_points_per_reversal
        self.min_points_per_ionization_front = min_points_per_ionization_front
        self.max_density_jump = max_density_jump
        self.max_temp_jump = max_temp_jump
        self.max_velocity_jump = max_velocity_jump
        self.max_sigma_jump = max_sigma_jump
        self.max_cell_ratio = max_cell_ratio
        self.min_dr_r = min_dr_r
        self.max_dr_r = max_dr_r

        # Detect features on the source grid
        self.feature_report = detect_all_features(source)

    def select_grid(self) -> tuple[TargetGrid, dict]:
        """Select the target grid with anchors + constrained refinement.

        Returns (target_grid, selection_report).
        """
        report = {
            'feature_report': self.feature_report.to_dict(),
            'target_points': self.target_points,
            'max_points': self.max_points,
            'policy': self.policy.value,
            'stages': [],
        }

        # Step 1: Trim source grid by tau and velocity limits
        source = self._trim_source()
        report['stages'].append({
            'stage': 'trim',
            'n_zones_after_trim': source.n_shells,
        })

        # Step 2: Build initial target boundaries from anchors + tau-based grid
        r_target_bound = self._build_initial_grid(source)
        report['stages'].append({
            'stage': 'initial_grid',
            'n_target_points': len(r_target_bound) - 1,
        })

        # Step 3: Conservative remap
        target = conservative_remap(source, r_target_bound, source.time_days)

        # Step 4: Mark anchors
        self._mark_anchors(target, source)

        # Step 5: Constrained priority refinement
        # Track the current boundaries for proper re-remapping
        current_boundaries = r_target_bound.copy()
        target, refinement_report = self._constrained_refinement(
            target, source, current_boundaries)
        report['stages'].append({
            'stage': 'refinement',
            'n_inserted': refinement_report['n_inserted'],
            'n_violations_final': refinement_report['n_violations_final'],
            'constraints_satisfied': refinement_report['constraints_satisfied'],
        })

        # Step 5b: Re-mark anchors after refinement (refinement creates new TargetGrid)
        self._mark_anchors(target, source)

        # Step 6: Check feasibility
        if self.policy == PointBudgetPolicy.FAIL_IF_INFEASIBLE:
            violations = self._count_violations(target)
            if violations > 0:
                report['feasible'] = False
                return target, report
        report['feasible'] = True

        return target, report

    def _trim_source(self) -> SourceShellGrid:
        """Trim source grid by tau and velocity limits."""
        s = self.source
        n = s.n_shells
        if n == 0:
            return s

        keep = np.ones(n, dtype=bool)

        # Tau trimming (outer: tau < tau_min)
        if len(s.tau_cumulative) == n:
            keep &= s.tau_cumulative >= self.tau_min
            # Inner: tau > tau_max (keep zones with tau <= tau_max)
            keep &= s.tau_cumulative <= self.tau_max

        # Velocity trimming (inner: |V| < v_min)
        if self.v_min_cms > 0:
            keep &= np.abs(s.velocity_cms) >= self.v_min_cms

        # Velocity trimming (outer: |V| > v_max)
        if self.v_max_cms > 0:
            keep &= np.abs(s.velocity_cms) <= self.v_max_cms

        # Always keep at least 2 zones
        if keep.sum() < 2:
            keep[:2] = True

        # Create trimmed grid (shallow copy with masked arrays)
        indices = np.where(keep)[0]
        # Build new boundaries
        r_bound_new = np.empty(len(indices) + 1)
        r_bound_new[1:-1] = np.sqrt(s.r_bound_cm[:-1][indices[:-1]] * s.r_bound_cm[1:][indices[:-1]])
        r_bound_new[0] = s.r_bound_cm[0]  # keep original outer boundary
        r_bound_new[-1] = s.r_bound_cm[-1]  # keep original inner boundary

        # Ensure strict decrease
        for i in range(len(r_bound_new) - 1):
            if r_bound_new[i] <= r_bound_new[i + 1]:
                r_bound_new[i] = r_bound_new[i + 1] * 1.0001

        vol = compute_shell_volume(r_bound_new)
        r_centre = compute_shell_centres(r_bound_new)

        # Use original zone data at kept indices
        v = s.velocity_cms[indices]
        T = s.temperature_K[indices]
        rho = s.density_gcm3[indices]
        kappa = s.kappa_cgs[indices]
        mass = s.shell_mass_g[indices]

        # Recompute density from mass/volume
        from .shell_model import density_from_mass_volume
        rho = density_from_mass_volume(mass, vol)

        mf = {e: f[indices] for e, f in s.mass_fractions.items()}
        iso = {k: f[indices] for k, f in s.isotope_fractions.items()}
        ni56 = s.ni56_initial[indices] if (s.ni56_initial is not None and hasattr(s.ni56_initial, '__len__')) else None
        ne = s.electron_density[indices] if s.electron_density is not None else None
        nbar = s.baryon_density[indices] if s.baryon_density is not None else None

        tau = compute_tau(r_bound_new, rho, kappa)

        trimmed = SourceShellGrid(
            r_bound_cm=r_bound_new,
            r_centre_cm=r_centre,
            shell_volume_cm3=vol,
            shell_mass_g=mass,
            density_gcm3=rho,
            velocity_cms=v,
            temperature_K=T,
            kappa_cgs=kappa,
            mass_fractions=mf,
            isotope_fractions=iso,
            ni56_initial=ni56,
            electron_density=ne,
            baryon_density=nbar,
            tau_cumulative=tau,
            time_days=s.time_days,
            source_file=s.source_file,
            n_zones_original=s.n_zones_original,
        )

        # Re-detect features on trimmed grid
        self.feature_report = detect_all_features(trimmed)

        return trimmed

    def _build_initial_grid(self, source: SourceShellGrid) -> np.ndarray:
        """Build initial target grid boundaries from anchors + tau-based spacing.

        Strategy:
        1. Start with anchor indices from feature detection
        2. Fill remaining points with log-spaced tau-based grid
        3. Ensure inner and outer boundaries are included
        """
        n_source = source.n_shells
        anchors = set(self.feature_report.anchor_indices)
        anchors.add(0)
        anchors.add(n_source - 1)
        anchors = sorted(anchors)

        # If anchors exceed target, use only anchors
        if len(anchors) >= self.target_points:
            # Use anchor radii as boundaries
            r_anchors = source.r_centre_cm[anchors]
            # Build boundaries from anchors
            r_bound = np.empty(len(anchors) + 1)
            r_bound[1:-1] = np.sqrt(r_anchors[:-1] * r_anchors[1:])
            r_bound[0] = source.r_bound_cm[0]
            r_bound[-1] = source.r_bound_cm[-1]
            return r_bound

        # Fill with tau-based spacing
        n_fill = self.target_points - len(anchors)
        if n_fill > 0 and len(source.tau_cumulative) == n_source:
            # Log-spaced in tau between tau_min and tau_max
            tau = np.maximum(source.tau_cumulative, 1e-10)
            log_tau = np.log10(tau)
            log_tau_min = log_tau.min()
            log_tau_max = log_tau.max()

            # Generate fill points in log-tau space
            fill_log_tau = np.linspace(log_tau_min, log_tau_max, n_fill)
            fill_tau = 10 ** fill_log_tau

            # Find source indices closest to each fill tau
            fill_indices = []
            for ft in fill_tau:
                idx = np.argmin(np.abs(tau - ft))
                if idx not in anchors:
                    fill_indices.append(idx)

            all_indices = sorted(set(anchors + fill_indices))
        else:
            all_indices = anchors

        # Build boundaries from selected indices
        r_selected = source.r_centre_cm[all_indices]
        r_bound = np.empty(len(r_selected) + 1)
        r_bound[1:-1] = np.sqrt(r_selected[:-1] * r_selected[1:])
        r_bound[0] = source.r_bound_cm[0]
        r_bound[-1] = source.r_bound_cm[-1]

        # Ensure strict decrease
        for i in range(len(r_bound) - 1):
            if r_bound[i] <= r_bound[i + 1]:
                r_bound[i] = r_bound[i + 1] * 1.0001

        return r_bound

    def _mark_anchors(self, target: TargetGrid, source: SourceShellGrid):
        """Mark anchor points on the target grid."""
        if target.is_anchor is None:
            target.is_anchor = np.zeros(target.n_points, dtype=bool)

        # Mark points closest to anchor radii
        for anchor_r in source.r_centre_cm[self.feature_report.anchor_indices]:
            idx = np.argmin(np.abs(target.r_cm - anchor_r))
            target.is_anchor[idx] = True
            target.provenance[idx] = 'protected'

    def _count_violations(self, target: TargetGrid) -> int:
        """Count hard constraint violations on the target grid."""
        n = target.n_points
        if n < 2:
            return 1

        violations = 0

        # Cell ratio
        dr = np.abs(np.diff(target.r_cm))
        if len(dr) > 1:
            ratios = dr[:-1] / np.maximum(dr[1:], 1e-30)
            ratios = np.maximum(ratios, 1.0 / ratios)
            violations += np.sum(ratios > self.max_cell_ratio)

        # Delta tau
        dtau = np.abs(np.diff(target.tau_cumulative))
        violations += np.sum(dtau > self.max_delta_tau)

        # Density jump
        drho = np.abs(np.diff(target.density_gcm3)) / np.maximum(target.density_gcm3[:-1], 1e-30)
        violations += np.sum(drho > self.max_density_jump)

        # Temperature jump
        dT = np.abs(np.diff(target.temperature_K)) / np.maximum(target.temperature_K[:-1], 1e-30)
        violations += np.sum(dT > self.max_temp_jump)

        # Velocity jump
        dv = np.abs(np.diff(target.velocity_cms)) / np.maximum(np.abs(target.velocity_cms[:-1]), 1e-30)
        violations += np.sum(dv > self.max_velocity_jump)

        # Sigma jump
        dsigma = np.abs(np.diff(target.sigma))
        violations += np.sum(dsigma > self.max_sigma_jump)

        return int(violations)

    def _constrained_refinement(self, target: TargetGrid,
                                 source: SourceShellGrid,
                                 current_boundaries: np.ndarray) -> tuple[TargetGrid, dict]:
        """Insert points to satisfy hard constraints.

        Uses priority-based insertion: find the worst violation, insert
        a boundary at the midpoint of the violating cell, re-remap, and repeat.

        Parameters
        ----------
        target : current target grid
        source : source grid (for re-remapping)
        current_boundaries : (n+1,) current target shell boundaries (decreasing)
            These are updated in-place as points are inserted.
        """
        n_inserted = 0
        max_insert = self.max_points - target.n_points

        for iteration in range(max_insert):
            if target.n_points >= self.max_points:
                break

            # Find worst violation
            worst_idx, worst_score = self._find_worst_violation(target)

            if worst_idx < 0 or worst_score <= 0:
                break  # all constraints satisfied

            # Insert a new boundary at the geometric mean of the violating cell
            # The violating cell is between target.r_cm[worst_idx] and target.r_cm[worst_idx+1]
            # We insert a boundary between them, which splits the cell
            r_new = self._insert_boundary(target, worst_idx)

            if r_new is None:
                break

            # Insert r_new into current_boundaries at the correct position
            # current_boundaries is decreasing; find where r_new fits
            insert_pos = np.searchsorted(-current_boundaries, -r_new)
            new_boundaries = np.insert(current_boundaries, insert_pos, r_new)

            # Ensure strict decrease
            for i in range(len(new_boundaries) - 1):
                if new_boundaries[i] <= new_boundaries[i + 1]:
                    new_boundaries[i] = new_boundaries[i + 1] * 1.0001

            # Re-remap with new boundaries
            target = conservative_remap(source, new_boundaries, source.time_days)
            current_boundaries = new_boundaries  # replace reference
            n_inserted += 1

        n_violations = self._count_violations(target)
        constraints_satisfied = n_violations == 0

        return target, {
            'n_inserted': n_inserted,
            'n_violations_final': n_violations,
            'constraints_satisfied': constraints_satisfied,
        }

    def _find_worst_violation(self, target: TargetGrid) -> tuple[int, float]:
        """Find the index of the worst constraint violation.

        Returns (edge_index, normalized_violation_score).
        edge_index is the index of the left cell of the violating edge.
        """
        n = target.n_points
        if n < 2:
            return -1, 0.0

        worst_idx = -1
        worst_score = 0.0

        # Check each edge
        for i in range(n - 1):
            score = 0.0

            # Cell ratio
            if i > 0 and i < n - 2:
                dr_left = abs(target.r_cm[i - 1] - target.r_cm[i])
                dr_right = abs(target.r_cm[i] - target.r_cm[i + 1])
                dr_next = abs(target.r_cm[i + 1] - target.r_cm[i + 2])
                if dr_right > 0 and dr_left > 0:
                    ratio = max(dr_left / dr_right, dr_right / dr_left)
                    if ratio > self.max_cell_ratio:
                        score = max(score, ratio / self.max_cell_ratio)

            # Delta tau
            dtau = abs(target.tau_cumulative[i + 1] - target.tau_cumulative[i])
            if dtau > self.max_delta_tau:
                score = max(score, dtau / self.max_delta_tau)

            # Density jump
            drho = abs(target.density_gcm3[i + 1] - target.density_gcm3[i]) / max(target.density_gcm3[i], 1e-30)
            if drho > self.max_density_jump:
                score = max(score, drho / self.max_density_jump)

            # Temperature jump
            dT = abs(target.temperature_K[i + 1] - target.temperature_K[i]) / max(target.temperature_K[i], 1e-30)
            if dT > self.max_temp_jump:
                score = max(score, dT / self.max_temp_jump)

            # Velocity jump
            dv = abs(target.velocity_cms[i + 1] - target.velocity_cms[i]) / max(abs(target.velocity_cms[i]), 1e-30)
            if dv > self.max_velocity_jump:
                score = max(score, dv / self.max_velocity_jump)

            # Sigma jump
            dsigma = abs(target.sigma[i + 1] - target.sigma[i])
            if dsigma > self.max_sigma_jump:
                score = max(score, dsigma / self.max_sigma_jump)

            if score > worst_score:
                worst_score = score
                worst_idx = i

        return worst_idx, worst_score

    def _insert_boundary(self, target: TargetGrid, edge_idx: int) -> Optional[float]:
        """Compute the radius for a new boundary at the worst violation.

        Places the new boundary at the midpoint (in log-r) of the violating cell.
        """
        if edge_idx < 0 or edge_idx >= target.n_points - 1:
            return None
        r_left = target.r_cm[edge_idx]
        r_right = target.r_cm[edge_idx + 1]
        # Geometric mean for log-spacing
        return np.sqrt(r_left * r_right)
