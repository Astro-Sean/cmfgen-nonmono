"""Feature detection: shocks, interfaces, velocity extrema, ionization fronts.

Detects physical features on the original STELLA grid BEFORE any smoothing
or downsampling. Feature masks are stored on the SourceShellGrid and
protected during all later operations.

Feature types detected:
- Strong shocks (density/pressure/velocity jumps)
- Contact discontinuities (density jump, no velocity jump)
- Composition interfaces (mass fraction gradients)
- Radioactive-isotope boundaries (Ni56 gradient)
- Velocity maxima and minima
- dV/dr sign changes (velocity reversals)
- Velocity zero crossings
- Ionization/recombination fronts (Ne, Ye, T transitions)
- Optical-depth landmarks (tau = 2/3, 1, 10, etc.)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from .shell_model import SourceShellGrid
from .solver_metrics import sigma_velocity, find_velocity_extrema, find_dvdr_sign_changes
from .opacity_tau import find_tau_landmarks, find_photosphere


@dataclass
class FeatureReport:
    """Summary of detected features."""
    n_shocks: int = 0
    n_contact_discontinuities: int = 0
    n_composition_interfaces: int = 0
    n_isotope_boundaries: int = 0
    n_velocity_maxima: int = 0
    n_velocity_minima: int = 0
    n_velocity_reversals: int = 0
    n_velocity_zero_crossings: int = 0
    n_ionization_fronts: int = 0
    n_tau_landmarks: int = 0

    # Detailed locations (indices into source grid)
    shock_indices: list[int] = field(default_factory=list)
    contact_indices: list[int] = field(default_factory=list)
    composition_interface_indices: list[int] = field(default_factory=list)
    isotope_boundary_indices: list[int] = field(default_factory=list)
    velocity_maxima_indices: list[int] = field(default_factory=list)
    velocity_minima_indices: list[int] = field(default_factory=list)
    velocity_reversal_indices: list[int] = field(default_factory=list)
    velocity_zero_crossing_indices: list[int] = field(default_factory=list)
    ionization_front_indices: list[int] = field(default_factory=list)
    tau_landmark_indices: list[int] = field(default_factory=list)
    photosphere_index: Optional[int] = None

    # Anchor indices (union of all mandatory anchors)
    anchor_indices: list[int] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            'n_shocks': self.n_shocks,
            'n_contact_discontinuities': self.n_contact_discontinuities,
            'n_composition_interfaces': self.n_composition_interfaces,
            'n_isotope_boundaries': self.n_isotope_boundaries,
            'n_velocity_maxima': self.n_velocity_maxima,
            'n_velocity_minima': self.n_velocity_minima,
            'n_velocity_reversals': self.n_velocity_reversals,
            'n_velocity_zero_crossings': self.n_velocity_zero_crossings,
            'n_ionization_fronts': self.n_ionization_fronts,
            'n_tau_landmarks': self.n_tau_landmarks,
            'shock_indices': self.shock_indices,
            'contact_indices': self.contact_indices,
            'composition_interface_indices': self.composition_interface_indices,
            'isotope_boundary_indices': self.isotope_boundary_indices,
            'velocity_maxima_indices': self.velocity_maxima_indices,
            'velocity_minima_indices': self.velocity_minima_indices,
            'velocity_reversal_indices': self.velocity_reversal_indices,
            'velocity_zero_crossing_indices': self.velocity_zero_crossing_indices,
            'ionization_front_indices': self.ionization_front_indices,
            'tau_landmark_indices': self.tau_landmark_indices,
            'photosphere_index': self.photosphere_index,
            'anchor_indices': sorted(set(self.anchor_indices)),
        }


def detect_shocks(grid: SourceShellGrid,
                  density_jump_threshold: float = 2.0,
                  velocity_jump_threshold: float = 0.5,
                  pressure_jump_threshold: float = 2.0) -> list[int]:
    """Detect strong shocks using Rankine-Hugoniot-like criteria.

    A shock is detected at shell *i* if the density jump between
    shells i and i+1 exceeds the threshold AND there is a corresponding
    velocity jump.  Contact discontinuities (density jump without
    velocity jump) are detected separately.

    Parameters
    ----------
    density_jump_threshold : rho[i]/rho[i+1] ratio above which a shock is detected
    velocity_jump_threshold : |dv|/v fractional velocity jump
    pressure_jump_threshold : P[i]/P[i+1] ratio (if pressure available)
    """
    n = grid.n_shells
    if n < 2:
        return []

    rho = grid.density_gcm3
    v = grid.velocity_cms

    shock_indices = []
    for i in range(n - 1):
        rho_ratio = rho[i] / max(rho[i + 1], 1e-30)
        rho_jump = rho_ratio > density_jump_threshold or rho_ratio < 1.0 / density_jump_threshold
        dv = abs(v[i] - v[i + 1]) / max(abs(v[i]), 1e-30)
        v_jump = dv > velocity_jump_threshold

        # Pressure jump if available
        p_jump = True  # default: don't require pressure
        if grid.pressure_cgs is not None:
            p_ratio = grid.pressure_cgs[i] / max(grid.pressure_cgs[i + 1], 1e-30)
            p_jump = p_ratio > pressure_jump_threshold or p_ratio < 1.0 / pressure_jump_threshold

        # Shock: density jump + velocity jump
        if rho_jump and v_jump:
            shock_indices.append(i)

    return shock_indices


def detect_contact_discontinuities(grid: SourceShellGrid,
                                   density_jump_threshold: float = 3.0,
                                   velocity_jump_threshold: float = 0.1) -> list[int]:
    """Detect contact discontinuities: density jump without velocity jump."""
    n = grid.n_shells
    if n < 2:
        return []

    rho = grid.density_gcm3
    v = grid.velocity_cms

    contact_indices = []
    for i in range(n - 1):
        rho_ratio = rho[i] / max(rho[i + 1], 1e-30)
        rho_jump = rho_ratio > density_jump_threshold or rho_ratio < 1.0 / density_jump_threshold
        dv = abs(v[i] - v[i + 1]) / max(abs(v[i]), 1e-30)
        v_small = dv < velocity_jump_threshold

        if rho_jump and v_small:
            contact_indices.append(i)

    return contact_indices


def detect_composition_interfaces(grid: SourceShellGrid,
                                  gradient_threshold: float = 0.5,
                                  min_abundance: float = 1e-4) -> list[int]:
    """Detect composition interfaces from mass fraction gradients.

    An interface is detected at shell *i* if any element's mass fraction
    changes by more than gradient_threshold (relative) between shells
    i and i+1, AND the element is above min_abundance in at least one shell.
    """
    n = grid.n_shells
    if n < 2 or not grid.mass_fractions:
        return []

    interface_indices = set()
    for elem, frac in grid.mass_fractions.items():
        if len(frac) != n:
            continue
        for i in range(n - 1):
            f0 = frac[i]
            f1 = frac[i + 1]
            if max(f0, f1) < min_abundance:
                continue
            # Relative change
            denom = max(abs(f0), abs(f1), 1e-30)
            rel_change = abs(f0 - f1) / denom
            if rel_change > gradient_threshold:
                interface_indices.add(i)

    return sorted(interface_indices)


def detect_isotope_boundaries(grid: SourceShellGrid,
                              gradient_threshold: float = 0.5,
                              min_abundance: float = 1e-6) -> list[int]:
    """Detect radioactive-isotope boundaries (Ni56, Co56)."""
    n = grid.n_shells
    if n < 2:
        return []

    boundary_indices = set()

    # Check Ni56 initial
    if grid.ni56_initial is not None and hasattr(grid.ni56_initial, '__len__') and len(grid.ni56_initial) == n:
        for i in range(n - 1):
            f0 = grid.ni56_initial[i]
            f1 = grid.ni56_initial[i + 1]
            if max(f0, f1) < min_abundance:
                continue
            denom = max(abs(f0), abs(f1), 1e-30)
            if abs(f0 - f1) / denom > gradient_threshold:
                boundary_indices.add(i)

    # Check isotope fractions
    for (species, baryon), frac in grid.isotope_fractions.items():
        if len(frac) != n:
            continue
        for i in range(n - 1):
            f0 = frac[i]
            f1 = frac[i + 1]
            if max(f0, f1) < min_abundance:
                continue
            denom = max(abs(f0), abs(f1), 1e-30)
            if abs(f0 - f1) / denom > gradient_threshold:
                boundary_indices.add(i)

    return sorted(boundary_indices)


def detect_velocity_extrema_and_reversals(grid: SourceShellGrid) -> dict:
    """Detect velocity maxima, minima, and dV/dr sign changes."""
    r = grid.r_centre_cm
    v = grid.velocity_cms

    extrema = find_velocity_extrema(r, v)
    reversals = find_dvdr_sign_changes(r, v)

    # Velocity zero crossings
    n = len(v)
    zero_crossings = []
    for i in range(1, n):
        if v[i - 1] * v[i] < 0:  # sign change
            zero_crossings.append(i)

    return {
        'maxima': extrema['maxima'],
        'minima': extrema['minima'],
        'reversals': reversals,
        'zero_crossings': zero_crossings,
    }


def detect_ionization_fronts(grid: SourceShellGrid,
                             ne_log_threshold: float = 1.0,
                             ye_log_threshold: float = 0.5,
                             t_log_threshold: float = 0.3,
                             kappa_log_threshold: float = 0.5) -> list[int]:
    """Detect ionization/recombination fronts.

    A front is detected where electron density, electron fraction,
    temperature, or opacity changes rapidly (in log space).

    Parameters
    ----------
    ne_log_threshold : |dlog10(Ne)| above which a front is detected
    ye_log_threshold : |dlog10(Ye)| threshold
    t_log_threshold : |dlog10(T)| threshold
    kappa_log_threshold : |dlog10(kappa)| threshold
    """
    n = grid.n_shells
    if n < 2:
        return []

    front_indices = set()

    # Electron density
    if grid.electron_density is not None and len(grid.electron_density) == n:
        ne = np.maximum(grid.electron_density, 1e-30)
        log_ne = np.log10(ne)
        for i in range(n - 1):
            if abs(log_ne[i] - log_ne[i + 1]) > ne_log_threshold:
                front_indices.add(i)

    # Electron fraction Ye = Ne / n_baryon
    ye = grid.electron_fraction
    if ye is not None and len(ye) == n:
        ye_safe = np.maximum(ye, 1e-10)
        log_ye = np.log10(ye_safe)
        for i in range(n - 1):
            if abs(log_ye[i] - log_ye[i + 1]) > ye_log_threshold:
                front_indices.add(i)

    # Temperature
    T = np.maximum(grid.temperature_K, 1.0)
    log_T = np.log10(T)
    for i in range(n - 1):
        if abs(log_T[i] - log_T[i + 1]) > t_log_threshold:
            front_indices.add(i)

    # Opacity
    if len(grid.kappa_cgs) == n:
        kappa = np.maximum(grid.kappa_cgs, 1e-30)
        log_kappa = np.log10(kappa)
        for i in range(n - 1):
            if abs(log_kappa[i] - log_kappa[i + 1]) > kappa_log_threshold:
                front_indices.add(i)

    return sorted(front_indices)


def detect_all_features(grid: SourceShellGrid,
                        density_jump_threshold: float = 2.0,
                        velocity_jump_threshold: float = 0.5,
                        composition_gradient_threshold: float = 0.5,
                        ne_log_threshold: float = 1.0,
                        ye_log_threshold: float = 0.5,
                        t_log_threshold: float = 0.3,
                        kappa_log_threshold: float = 0.5) -> FeatureReport:
    """Detect all features and populate masks on the grid.

    Returns a FeatureReport with all detected features and the union
    of mandatory anchor indices.
    """
    n = grid.n_shells
    report = FeatureReport()

    # Shocks
    shocks = detect_shocks(grid, density_jump_threshold, velocity_jump_threshold)
    report.shock_indices = shocks
    report.n_shocks = len(shocks)
    grid.shock_mask = np.zeros(n, dtype=bool)
    for i in shocks:
        grid.shock_mask[i] = True

    # Contact discontinuities
    contacts = detect_contact_discontinuities(grid, density_jump_threshold=3.0)
    report.contact_indices = contacts
    report.n_contact_discontinuities = len(contacts)
    grid.contact_discontinuity_mask = np.zeros(n, dtype=bool)
    for i in contacts:
        grid.contact_discontinuity_mask[i] = True

    # Composition interfaces
    comp_ifaces = detect_composition_interfaces(grid, composition_gradient_threshold)
    report.composition_interface_indices = comp_ifaces
    report.n_composition_interfaces = len(comp_ifaces)
    grid.composition_interface_mask = np.zeros(n, dtype=bool)
    for i in comp_ifaces:
        grid.composition_interface_mask[i] = True

    # Isotope boundaries
    iso_bounds = detect_isotope_boundaries(grid)
    report.isotope_boundary_indices = iso_bounds
    report.n_isotope_boundaries = len(iso_bounds)
    grid.isotope_boundary_mask = np.zeros(n, dtype=bool)
    for i in iso_bounds:
        grid.isotope_boundary_mask[i] = True

    # Velocity extrema and reversals
    vel_features = detect_velocity_extrema_and_reversals(grid)
    report.velocity_maxima_indices = vel_features['maxima']
    report.velocity_minima_indices = vel_features['minima']
    report.velocity_reversal_indices = vel_features['reversals']
    report.velocity_zero_crossing_indices = vel_features['zero_crossings']
    report.n_velocity_maxima = len(vel_features['maxima'])
    report.n_velocity_minima = len(vel_features['minima'])
    report.n_velocity_reversals = len(vel_features['reversals'])
    report.n_velocity_zero_crossings = len(vel_features['zero_crossings'])

    grid.velocity_extremum_mask = np.zeros(n, dtype=bool)
    for i in vel_features['maxima'] + vel_features['minima']:
        grid.velocity_extremum_mask[i] = True

    grid.velocity_reversal_mask = np.zeros(n, dtype=bool)
    for i in vel_features['reversals']:
        grid.velocity_reversal_mask[i] = True

    # Ionization fronts
    ion_fronts = detect_ionization_fronts(grid, ne_log_threshold, ye_log_threshold,
                                          t_log_threshold, kappa_log_threshold)
    report.ionization_front_indices = ion_fronts
    report.n_ionization_fronts = len(ion_fronts)
    grid.ionization_front_mask = np.zeros(n, dtype=bool)
    for i in ion_fronts:
        grid.ionization_front_mask[i] = True

    # Tau landmarks
    if len(grid.tau_cumulative) == n:
        landmarks = find_tau_landmarks(grid.tau_cumulative, grid.r_centre_cm)
        tau_indices = [lm['index'] for lm in landmarks.values()]
        report.tau_landmark_indices = tau_indices
        report.n_tau_landmarks = len(tau_indices)
        phot = find_photosphere(grid.tau_cumulative, grid.r_centre_cm)
        if phot:
            report.photosphere_index = phot['index']

    # Build anchor union: inner boundary, outer boundary, all features
    anchors = set()
    anchors.add(0)       # outer boundary
    anchors.add(n - 1)   # inner boundary
    anchors.update(shocks)
    anchors.update(contacts)
    anchors.update(comp_ifaces)
    anchors.update(iso_bounds)
    anchors.update(vel_features['maxima'])
    anchors.update(vel_features['minima'])
    anchors.update(vel_features['reversals'])
    anchors.update(vel_features['zero_crossings'])
    anchors.update(ion_fronts)
    anchors.update(tau_indices)
    report.anchor_indices = sorted(anchors)

    return report
