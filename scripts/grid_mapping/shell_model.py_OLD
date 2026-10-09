"""Conservative shell representation for STELLA-to-CMFGEN grid mapping.

Provides immutable source-shell and mutable target-grid dataclasses that
carry all physical quantities needed for conservative remapping, feature
detection, and CMFGEN output generation.

Conventions
-----------
* All arrays are **outer-to-inner** (decreasing radius): ``r[0] > r[1] > ... > r[-1]``.
* Shell *boundaries* are the radius points themselves.  Shell *i* spans
  ``[r[i+1], r[i]]`` (inner edge to outer edge).
* Velocity may be **non-monotonic** — only radius must be strictly decreasing.
* Velocity sign is preserved (no ``abs()``).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

# Physical constants (CGS)
C_LIGHT = 2.99792458e10   # cm/s
K_BOLTZ = 1.380649e-16    # erg/K
M_H     = 1.6735575e-24   # g
SIGMA_SB = 5.670374419e-5 # erg/cm^2/s/K^4
LSUN    = 3.828e33        # erg/s
MSUN    = 1.989e33        # g
AVOGADRO = 6.02214076e23

# STELLA Natom=15 element order
STELLA_ELEM_NAMES = ['H', 'He', 'C', 'N', 'O', 'Ne', 'Na', 'Mg', 'Al', 'Si',
                     'S', 'Ar', 'Ca', 'Fe', 'Ni']
STELLA_AZ = np.array([1.0, 4.0, 12.0, 14.0, 16.0, 20.0, 23.0, 24.0, 27.0,
                      28.0, 32.0, 40.0, 40.0, 56.0, 58.0])


@dataclass
class SourceShellGrid:
    """Immutable conservative representation of the original STELLA shells.

    Each shell *i* occupies the volume between ``r_bound[i]`` (outer edge)
    and ``r_bound[i+1]`` (inner edge).  All quantities are cell-centred or
    cell-integrated values for that shell.

    Arrays are outer-to-inner (index 0 = outermost shell).
    """

    # --- Geometry (n_shells+1 boundaries, n_shells centres) ---
    r_bound_cm: np.ndarray         # (n_shells+1,) shell boundaries, decreasing
    r_centre_cm: np.ndarray        # (n_shells,)   shell centre radii
    shell_volume_cm3: np.ndarray   # (n_shells,)   4π/3 (r_out^3 - r_in^3)
    shell_mass_g: np.ndarray       # (n_shells,)   mass per shell

    # --- Thermodynamic state (n_shells,) ---
    density_gcm3: np.ndarray       # shell_mass / shell_volume
    velocity_cms: np.ndarray       # cm/s, sign preserved
    temperature_K: np.ndarray
    pressure_cgs: Optional[np.ndarray] = None
    internal_energy: Optional[np.ndarray] = None

    # --- Composition ---
    mass_fractions: dict[str, np.ndarray] = field(default_factory=dict)
    # Isotope mass fractions: key = (species_name, baryon_number)
    isotope_fractions: dict[tuple[str, int], np.ndarray] = field(default_factory=dict)
    ni56_initial: Optional[np.ndarray] = None  # initial Ni56 mass fraction (pre-decay)

    # --- Opacity and optical depth ---
    kappa_cgs: np.ndarray = field(default_factory=lambda: np.array([]))  # cm^2/g
    tau_cumulative: np.ndarray = field(default_factory=lambda: np.array([]))  # integrated from outside

    # --- Ionization ---
    electron_density: Optional[np.ndarray] = None  # cm^-3
    baryon_density: Optional[np.ndarray] = None    # cm^-3

    # --- Feature flags (set by feature_detection) ---
    shock_mask: Optional[np.ndarray] = None
    contact_discontinuity_mask: Optional[np.ndarray] = None
    composition_interface_mask: Optional[np.ndarray] = None
    velocity_extremum_mask: Optional[np.ndarray] = None
    velocity_reversal_mask: Optional[np.ndarray] = None  # dV/dr sign change
    ionization_front_mask: Optional[np.ndarray] = None
    isotope_boundary_mask: Optional[np.ndarray] = None

    # --- Metadata ---
    time_days: float = 0.0
    source_file: str = ""
    n_zones_original: int = 0

    @property
    def n_shells(self) -> int:
        return len(self.r_centre_cm)

    @property
    def electron_fraction(self) -> Optional[np.ndarray]:
        """Ye = Ne / n_baryon."""
        if self.electron_density is not None and self.baryon_density is not None:
            return self.electron_density / np.maximum(self.baryon_density, 1e-30)
        return None

    def validate(self) -> list[str]:
        """Return list of validation error messages (empty if OK)."""
        errors = []
        n = self.n_shells
        if n == 0:
            errors.append("Empty grid (n_shells=0)")
            return errors

        # Radius ordering
        rb = self.r_bound_cm
        if len(rb) != n + 1:
            errors.append(f"r_bound has {len(rb)} elements, expected {n+1}")
        if not np.all(np.diff(rb) < 0):
            n_bad = int(np.sum(np.diff(rb) >= 0))
            errors.append(f"r_bound not strictly decreasing ({n_bad} violations)")

        # Array lengths
        for name, arr in [("density", self.density_gcm3), ("velocity", self.velocity_cms),
                          ("temperature", self.temperature_K), ("shell_mass", self.shell_mass_g),
                          ("shell_volume", self.shell_volume_cm3), ("r_centre", self.r_centre_cm)]:
            if len(arr) != n:
                errors.append(f"{name} has {len(arr)} elements, expected {n}")

        # Finite checks
        for name, arr in [("r_bound", rb), ("r_centre", self.r_centre_cm),
                          ("density", self.density_gcm3), ("velocity", self.velocity_cms),
                          ("temperature", self.temperature_K)]:
            if not np.all(np.isfinite(arr)):
                n_bad = int(np.sum(~np.isfinite(arr)))
                errors.append(f"{name} has {n_bad} non-finite values")

        # Physical sign checks
        if np.any(self.density_gcm3 <= 0):
            n_bad = int(np.sum(self.density_gcm3 <= 0))
            errors.append(f"density has {n_bad} non-positive values")
        if np.any(self.temperature_K <= 0):
            n_bad = int(np.sum(self.temperature_K <= 0))
            errors.append(f"temperature has {n_bad} non-positive values")
        if np.any(self.shell_mass_g < 0):
            n_bad = int(np.sum(self.shell_mass_g < 0))
            errors.append(f"shell_mass has {n_bad} negative values")
        if np.any(self.shell_volume_cm3 <= 0):
            n_bad = int(np.sum(self.shell_volume_cm3 <= 0))
            errors.append(f"shell_volume has {n_bad} non-positive values")

        # Mass fraction sums
        if self.mass_fractions:
            mf_sum = np.zeros(n)
            for frac in self.mass_fractions.values():
                if len(frac) != n:
                    errors.append(f"mass_fraction array length mismatch")
                    break
                mf_sum += frac
            if not errors or "mass_fraction array length mismatch" not in errors:
                if np.any(np.abs(mf_sum - 1.0) > 0.01):
                    errors.append(f"mass fractions deviate from 1.0 (max dev={np.max(np.abs(mf_sum-1.0)):.4f})")

        return errors


@dataclass
class TargetGrid:
    """Mutable CMFGEN target grid built by grid_selection and conservative_remap.

    Points are cell boundaries (outer-to-inner, decreasing R).
    Provenance tracks whether each point was original, inserted, or protected.
    """

    r_cm: np.ndarray              # (n_points,) boundary radii, decreasing
    velocity_cms: np.ndarray      # (n_points,) cm/s, sign preserved
    temperature_K: np.ndarray
    density_gcm3: np.ndarray
    kappa_cgs: np.ndarray         # cm^2/g
    sigma: np.ndarray             # dlnV/dlnR - 1 (CMFGEN definition)
    tau_cumulative: np.ndarray    # integrated from outside

    mass_fractions: dict[str, np.ndarray] = field(default_factory=dict)
    isotope_fractions: dict[tuple[str, int], np.ndarray] = field(default_factory=dict)
    ni56_initial: Optional[np.ndarray] = None

    electron_density: Optional[np.ndarray] = None
    atom_density: Optional[np.ndarray] = None

    # Provenance: 'original', 'inserted', 'interpolated', 'extrapolated', 'protected'
    provenance: np.ndarray = field(default_factory=lambda: np.array([], dtype='<U16'))

    # Shell masses (for conservation checking)
    shell_mass_g: Optional[np.ndarray] = None
    shell_volume_cm3: Optional[np.ndarray] = None

    # Feature anchor flags
    is_anchor: Optional[np.ndarray] = None

    time_days: float = 0.0

    @property
    def n_points(self) -> int:
        return len(self.r_cm)

    def validate(self) -> list[str]:
        """Return list of validation error messages (empty if OK)."""
        errors = []
        n = self.n_points
        if n == 0:
            errors.append("Empty target grid")
            return errors

        # Radius strict monotonicity (decreasing)
        if not np.all(np.diff(self.r_cm) < 0):
            n_bad = int(np.sum(np.diff(self.r_cm) >= 0))
            errors.append(f"radius not strictly decreasing ({n_bad} violations)")

        # All arrays same length
        for name, arr in [("velocity", self.velocity_cms), ("temperature", self.temperature_K),
                          ("density", self.density_gcm3), ("kappa", self.kappa_cgs),
                          ("sigma", self.sigma), ("tau", self.tau_cumulative)]:
            if len(arr) != n:
                errors.append(f"{name} has {len(arr)} elements, expected {n}")

        # Finite
        for name, arr in [("r", self.r_cm), ("velocity", self.velocity_cms),
                          ("temperature", self.temperature_K), ("density", self.density_gcm3),
                          ("sigma", self.sigma)]:
            if not np.all(np.isfinite(arr)):
                n_bad = int(np.sum(~np.isfinite(arr)))
                errors.append(f"{name} has {n_bad} non-finite values")

        # Physical signs
        if np.any(self.density_gcm3 <= 0):
            n_bad = int(np.sum(self.density_gcm3 <= 0))
            errors.append(f"density has {n_bad} non-positive values")
        if np.any(self.temperature_K <= 0):
            n_bad = int(np.sum(self.temperature_K <= 0))
            errors.append(f"temperature has {n_bad} non-positive values")

        # Provenance length
        if len(self.provenance) != n:
            errors.append(f"provenance has {len(self.provenance)} elements, expected {n}")

        return errors


def compute_shell_volume(r_bound_cm: np.ndarray) -> np.ndarray:
    """Compute shell volumes from boundary radii.

    Shell *i* spans [r_bound[i+1], r_bound[i]] (outer-to-inner).
    V_i = 4π/3 * (r_bound[i]^3 - r_bound[i+1]^3)
    """
    r = np.asarray(r_bound_cm, dtype=float)
    if len(r) < 2:
        return np.array([])
    vol = (4.0 * np.pi / 3.0) * (r[:-1]**3 - r[1:]**3)
    return np.maximum(vol, 0.0)  # guard against roundoff


def compute_shell_centres(r_bound_cm: np.ndarray) -> np.ndarray:
    """Compute volume-weighted shell centre radii."""
    r = np.asarray(r_bound_cm, dtype=float)
    if len(r) < 2:
        return np.array([])
    # Volume-weighted centre: r_centre = ((r_out^3 + r_out^2*r_in + r_out*r_in^2 + r_in^3)/4)^(1/3)
    # This is the centroid of a spherical shell in radius.
    r_out = r[:-1]
    r_in = r[1:]
    r_cubed = (r_out**3 + r_out**2 * r_in + r_out * r_in**2 + r_in**3) / 4.0
    return np.cbrt(np.maximum(r_cubed, 0.0))


def density_from_mass_volume(mass_g: np.ndarray, volume_cm3: np.ndarray) -> np.ndarray:
    """Compute density = mass / volume with safe division."""
    vol = np.maximum(np.asarray(volume_cm3, dtype=float), 1e-99)
    return np.asarray(mass_g, dtype=float) / vol
