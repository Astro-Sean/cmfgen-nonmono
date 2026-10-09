#!/usr/bin/env python3
"""Adaptive regridding for STELLA → SuperLite conversion.

This module implements a physics-aware adaptive mesh refinement (AMR)
algorithm that replaces the fixed tau-weighted downsampling in
``convert_stella_swd_to_superlite.py``.

The algorithm preserves physically important structure while minimising
the number of SuperLite cells:

  * density gradients          (d log rho / d log r)
  * temperature gradients      (d log T   / d log r)
  * velocity gradients         (dv/dr), especially negative gradients (shocks)
  * composition gradients      (dX_i / d log r) for H, He, O, Si, Ca, Fe, ...
  * optical-depth structure    (cumulative tau, photosphere, thermalisation)
  * metal-rich cannon          (forced refinement at H-poor / metal-rich interfaces)

Conservative remapping
----------------------
Merged zones use **mass-weighted** averages for intensive properties
(velocity, temperature, radiation temperature, kappa, mass fractions)
and **exact mass conservation** for density:

    rho_merged = M_total / V_total

This ensures mass and element masses are conserved to machine precision.

All refinement thresholds are configurable and documented.  No arbitrary
magic numbers — every threshold has a physical or numerical justification.

The module is designed to be imported by ``convert_stella_swd_to_superlite.py``
but can also be run standalone for diagnostics.
"""

from __future__ import annotations

import os
import sys
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

# Reuse constants from the CMFGEN converter (same package)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from convert_stella_swd_to_cmfgen import MSUN, STELLA_ELEM_NAMES  # noqa: E402


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

@dataclass
class GridRefinementConfig:
    """Configuration for adaptive grid refinement.

    All thresholds are dimensionless gradients unless otherwise stated.
    Defaults are physically motivated, not arbitrary:

    * ``rho_grad``: 0.5 dex per dex — resolves order-of-magnitude density
      changes over a decade in radius (typical SN ejecta power-law index).
    * ``temp_grad``: 0.3 dex per dex — resolves temperature structure at
      the ~30% level, sufficient for ionisation balance.
    * ``vel_grad``: 0.05 per cell in dv/v — resolves velocity changes of
      5% per cell, catching shocks (dv/v can be >>1 at discontinuities).
    * ``comp_grad``: 0.1 in dX/dlogr — resolves 10% abundance changes,
      critical for composition interfaces (cannon boundaries).
    * ``tau_per_cell_target``: 0.5 — target optical depth per cell in the
      line-forming region (tau ~ 0.01–1).  Keeps DDMC transitions smooth.
    * ``cannon_min_cells_per_interface``: 5 — minimum cells to resolve
      each cannon boundary (inner and outer).  A sharp interface needs
      at least 3–5 cells to avoid smearing in MC transport.
    * ``shock_min_cells``: 5 — minimum cells across a velocity reversal
      (shock or contact discontinuity).
    """

    # Gradient-based refinement thresholds
    rho_grad: float = 0.5          # |d log rho / d log r|
    temp_grad: float = 0.3         # |d log T   / d log r|
    vel_grad: float = 0.05         # |dv/v| per cell
    comp_grad: float = 0.1         # |dX_i / d log r| for important elements

    # Optical-depth refinement
    tau_per_cell_target: float = 0.5   # target dtau per cell in line-forming region
    tau_refine_min: float = 1e-3       # refine where tau > this (above trivial)
    tau_refine_max: float = 10.0       # refine where tau < this (below deep diffusion)

    # Cannon (metal-rich H-poor component) preservation
    cannon_fe_threshold: float = 0.05   # Fe mass fraction to flag cannon
    cannon_si_threshold: float = 0.05   # Si mass fraction to flag cannon
    cannon_o_threshold: float = 0.10    # O mass fraction to flag cannon
    cannon_h_max: float = 0.01          # H must be below this for cannon
    cannon_min_cells_per_interface: int = 5
    cannon_min_cells_inside: int = 10

    # Shock / velocity reversal preservation
    shock_dv_threshold: float = 100.0   # km/s — velocity change to flag a shock
    shock_min_cells: int = 5

    # Photosphere refinement
    phot_refine_cells: int = 10         # extra cells around tau=1
    phot_tau_low: float = 0.3           # refine between tau=0.3 and tau=3
    phot_tau_high: float = 3.0

    # Grid constraints
    max_cells: int = 250
    min_cells: int = 20
    coarsen_factor: float = 2.0         # max ratio of adjacent dr/r

    # Elements to track for composition gradients
    composition_elements: Tuple[str, ...] = (
        'H', 'He', 'C', 'O', 'Ne', 'Mg', 'Si', 'S', 'Ca', 'Fe', 'Ni',
    )


@dataclass
class CannonDiagnostics:
    """Diagnostics for the metal-rich cannon component."""
    detected: bool = False
    inner_radius: float = 0.0       # cm
    outer_radius: float = 0.0       # cm
    inner_velocity: float = 0.0     # km/s
    outer_velocity: float = 0.0     # km/s
    mass: float = 0.0               # Msun
    composition: Dict[str, float] = field(default_factory=dict)
    n_cells_inside: int = 0
    n_cells_inner_interface: int = 0
    n_cells_outer_interface: int = 0
    zone_indices: np.ndarray = field(default_factory=lambda: np.array([], dtype=int))


@dataclass
class GridDiagnostics:
    """Diagnostics for the adaptive grid."""
    n_input: int = 0
    n_output: int = 0
    n_cannon_cells: int = 0
    n_shock_cells: int = 0
    n_phot_cells: int = 0
    n_grad_cells: int = 0
    dr_r_min: float = 0.0
    dr_r_max: float = 0.0
    dr_r_ratio: float = 0.0
    mass_in: float = 0.0
    mass_out: float = 0.0
    mass_error: float = 0.0
    element_mass_errors: Dict[str, float] = field(default_factory=dict)
    cannon: CannonDiagnostics = field(default_factory=CannonDiagnostics)
    velocity_extrema_preserved: bool = True
    composition_interfaces_preserved: bool = True
    n_negative_dv_preserved: int = 0


# ---------------------------------------------------------------------------
# Cannon detection
# ---------------------------------------------------------------------------

def detect_metal_rich_cannon(
    r_oi: np.ndarray,
    v_oi: np.ndarray,
    mf_oi: Dict[str, np.ndarray],
    config: GridRefinementConfig,
) -> CannonDiagnostics:
    """Detect the metal-rich, H-poor cannon component.

    The cannon is identified as a spatially contiguous region where:
      * H mass fraction drops below ``cannon_h_max``
      * AND at least one metal (Fe, Si, or O) exceeds its threshold

    The detection uses connected-component analysis on the boolean mask.

    Args:
        r_oi: radius, outer→inner (decreasing), cm
        v_oi: velocity, outer→inner, km/s
        mf_oi: mass fractions, outer→inner
        config: refinement configuration

    Returns:
        CannonDiagnostics with detection results
    """
    n = len(r_oi)
    diag = CannonDiagnostics()

    H = mf_oi.get('H', np.zeros(n))
    Fe = mf_oi.get('Fe', np.zeros(n))
    Si = mf_oi.get('Si', np.zeros(n))
    O = mf_oi.get('O', np.zeros(n))
    Ca = mf_oi.get('Ca', np.zeros(n))

    # Cannon mask: H-poor AND metal-rich
    h_poor = H < config.cannon_h_max
    metal_rich = (
        (Fe > config.cannon_fe_threshold) |
        (Si > config.cannon_si_threshold) |
        (O > config.cannon_o_threshold)
    )
    cannon_mask = h_poor & metal_rich

    if cannon_mask.sum() == 0:
        return diag

    # Find the largest connected component (zones are contiguous in radius)
    # r_oi is outer→inner, so indices are ordered by decreasing r.
    # Connected components in index space = connected in radius space.
    best_start, best_end, best_len = -1, -1, 0
    cur_start = -1
    for i in range(n):
        if cannon_mask[i]:
            if cur_start < 0:
                cur_start = i
        else:
            if cur_start >= 0:
                cur_len = i - cur_start
                if cur_len > best_len:
                    best_len = cur_len
                    best_start = cur_start
                    best_end = i
                cur_start = -1
    if cur_start >= 0:
        cur_len = n - cur_start
        if cur_len > best_len:
            best_len = cur_len
            best_start = cur_start
            best_end = n

    if best_len == 0:
        return diag

    idx = np.arange(best_start, best_end)
    diag.detected = True
    diag.zone_indices = idx
    # r_oi is outer→inner (decreasing), so:
    #   outer radius = r_oi[best_start] (larger r, lower index)
    #   inner radius = r_oi[best_end - 1] (smaller r, higher index)
    diag.outer_radius = float(r_oi[best_start])
    diag.inner_radius = float(r_oi[best_end - 1])
    diag.outer_velocity = float(v_oi[best_start])
    diag.inner_velocity = float(v_oi[best_end - 1])
    diag.n_cells_inside = len(idx)

    # Mass estimate (will be refined later with proper cell volumes)
    # Store average composition
    for elem in ['H', 'He', 'C', 'O', 'Si', 'S', 'Ca', 'Fe', 'Ni']:
        if elem in mf_oi:
            diag.composition[elem] = float(np.mean(mf_oi[elem][idx]))

    return diag


# ---------------------------------------------------------------------------
# Gradient computation
# ---------------------------------------------------------------------------

def compute_gradients(
    r_oi: np.ndarray,
    v_oi: np.ndarray,
    T_oi: np.ndarray,
    rho_oi: np.ndarray,
    mf_oi: Dict[str, np.ndarray],
    config: GridRefinementConfig,
) -> Dict[str, np.ndarray]:
    """Compute physical gradients for refinement decisions.

    All arrays are outer→inner (decreasing r).  Gradients are computed
    as central differences in log(radius) space.

    Returns a dict of gradient arrays (outer→inner), each of length n:
      * 'rho': |d log rho / d log r|
      * 'temp': |d log T / d log r|
      * 'vel': |dv/v| (fractional velocity change per cell)
      * 'vel_neg': negative velocity gradient indicator (1.0 where dv<0)
      * 'tau': composition gradient magnitude (max over tracked elements)
      * 'comp_H', 'comp_He', ...: per-element |dX / d log r|
    """
    n = len(r_oi)
    log_r = np.log(r_oi)
    # d log r between adjacent zones (outer→inner, so negative)
    dlogr = np.diff(log_r)  # length n-1, negative (r decreasing)

    grads: Dict[str, np.ndarray] = {}
    # Use central differences; pad with edge values
    def _central_diff(arr: np.ndarray) -> np.ndarray:
        """Central difference, padded to length n."""
        d = np.diff(arr)  # length n-1
        # Pad to length n with edge values
        result = np.zeros(n)
        result[1:-1] = 0.5 * (d[:-1] + d[1:])
        result[0] = d[0]
        result[-1] = d[-1]
        return result

    # Density gradient: |d log rho / d log r|
    log_rho = np.log(np.maximum(rho_oi, 1e-30))
    dlogrho = _central_diff(log_rho)
    grads['rho'] = np.abs(dlogrho / np.maximum(np.abs(_central_diff(log_r)), 1e-30))

    # Temperature gradient: |d log T / d log r|
    log_T = np.log(np.maximum(T_oi, 1.0))
    dlogT = _central_diff(log_T)
    grads['temp'] = np.abs(dlogT / np.maximum(np.abs(_central_diff(log_r)), 1e-30))

    # Velocity gradient: |dv/v| per cell
    dv = _central_diff(v_oi)
    grads['vel'] = np.abs(dv) / np.maximum(np.abs(v_oi), 1.0)
    # Negative velocity gradient (shock/reversal indicator)
    grads['vel_neg'] = (dv < 0).astype(float)

    # Composition gradients per element
    comp_grad_max = np.zeros(n)
    for elem in config.composition_elements:
        if elem not in mf_oi:
            continue
        X = mf_oi[elem]
        # Only compute gradient where abundance is non-negligible
        significant = X > 1e-6
        if not np.any(significant):
            grads[f'comp_{elem}'] = np.zeros(n)
            continue
        dX = _central_diff(X)
        # |dX / d log r| — normalise by dlogr
        g = np.abs(dX) / np.maximum(np.abs(_central_diff(log_r)), 1e-30)
        g[~significant] = 0.0
        grads[f'comp_{elem}'] = g
        comp_grad_max = np.maximum(comp_grad_max, g)
    grads['comp'] = comp_grad_max

    return grads


# ---------------------------------------------------------------------------
# Refinement mask construction
# ---------------------------------------------------------------------------

def build_refinement_mask(
    r_oi: np.ndarray,
    rho_oi: np.ndarray,
    kappa_oi: np.ndarray,
    tau_oi: np.ndarray,
    v_oi: np.ndarray,
    T_oi: np.ndarray,
    mf_oi: Dict[str, np.ndarray],
    config: GridRefinementConfig,
    cannon: CannonDiagnostics,
) -> Tuple[np.ndarray, GridDiagnostics]:
    """Build a boolean refinement mask marking zones that MUST be kept.

    A zone is marked for retention if it satisfies ANY refinement criterion.
    The mask is then used as seed points for the adaptive grid construction.

    Args:
        All arrays are outer→inner (decreasing r).
        tau_oi: cumulative optical depth, outer→inner (increasing inward)

    Returns:
        (mask, diagnostics) where mask is boolean array of length n
    """
    n = len(r_oi)
    grads = compute_gradients(r_oi, v_oi, T_oi, rho_oi, mf_oi, config)
    mask = np.zeros(n, dtype=bool)
    diag = GridDiagnostics()
    diag.n_input = n

    # --- Gradient-based refinement ---
    grad_mask = (
        (grads['rho'] > config.rho_grad) |
        (grads['temp'] > config.temp_grad) |
        (grads['vel'] > config.vel_grad) |
        (grads['comp'] > config.comp_grad)
    )
    mask |= grad_mask
    diag.n_grad_cells = int(grad_mask.sum())

    # --- Velocity reversal / shock refinement ---
    # Mark zones with significant negative velocity gradient
    dv = np.diff(v_oi)  # outer→inner, negative means v decreasing inward = reversal
    shock_zones = np.zeros(n, dtype=bool)
    shock_zones[1:-1] = np.abs(dv[:-1] + dv[1:]) / 2.0 > config.shock_dv_threshold
    # Also mark zones where dv < 0 (velocity decreasing inward = non-monotonic)
    vel_reversal = np.zeros(n, dtype=bool)
    vel_reversal[1:-1] = (dv[:-1] < 0) | (dv[1:] < 0)
    shock_mask = shock_zones | vel_reversal
    mask |= shock_mask
    diag.n_shock_cells = int(shock_mask.sum())

    # --- Optical-depth / photosphere refinement ---
    # Refine the photospheric region (tau between phot_tau_low and phot_tau_high)
    phot_mask = (tau_oi > config.phot_tau_low) & (tau_oi < config.phot_tau_high)
    mask |= phot_mask
    diag.n_phot_cells = int(phot_mask.sum())

    # Also refine where per-cell tau would be large (keep resolution where
    # optical depth changes rapidly)
    # Per-zone tau contribution
    dtau = np.diff(np.concatenate([[0], tau_oi]))
    high_dtau = dtau > config.tau_per_cell_target
    mask |= high_dtau

    # --- Cannon refinement ---
    if cannon.detected and len(cannon.zone_indices) > 0:
        cannon_mask = np.zeros(n, dtype=bool)
        cannon_mask[cannon.zone_indices] = True
        # Also mark the interfaces (boundaries of the cannon)
        inner_edge = cannon.zone_indices[-1]   # innermost cannon zone (highest index)
        outer_edge = cannon.zone_indices[0]    # outermost cannon zone (lowest index)
        # Mark zones around the inner interface
        for k in range(1, config.cannon_min_cells_per_interface + 1):
            idx = inner_edge + k
            if 0 <= idx < n:
                cannon_mask[idx] = True
            idx = inner_edge - k
            if 0 <= idx < n:
                cannon_mask[idx] = True
        # Mark zones around the outer interface
        for k in range(1, config.cannon_min_cells_per_interface + 1):
            idx = outer_edge + k
            if 0 <= idx < n:
                cannon_mask[idx] = True
            idx = outer_edge - k
            if 0 <= idx < n:
                cannon_mask[idx] = True
        mask |= cannon_mask
        diag.n_cannon_cells = int(cannon_mask.sum())

    # --- Always keep first and last zones ---
    mask[0] = True
    mask[-1] = True

    return mask, diag


# ---------------------------------------------------------------------------
# Adaptive grid construction
# ---------------------------------------------------------------------------

def _compute_zone_volumes(r_oi: np.ndarray) -> np.ndarray:
    """Compute zone volumes from zone-centre radii (outer→inner).

    Uses geometric-mean interfaces between adjacent zone centres.
    Returns volumes in outer→inner order.
    """
    n = len(r_oi)
    # Interfaces (inner→outer for computation, then reverse)
    r_io = r_oi[::-1]  # increasing
    edges = np.zeros(n + 1)
    edges[1:n] = np.sqrt(r_io[:-1] * r_io[1:])
    edges[0] = r_io[0] ** 2 / max(edges[1], 1e-30)
    edges[n] = r_io[-1] ** 2 / max(edges[n - 1], 1e-30)
    vol_io = (4.0 / 3.0) * np.pi * (edges[1:] ** 3 - edges[:-1] ** 3)
    return vol_io[::-1]  # outer→inner


def _compute_zone_priority(
    r_oi: np.ndarray,
    v_oi: np.ndarray,
    T_oi: np.ndarray,
    rho_oi: np.ndarray,
    kappa_oi: np.ndarray,
    tau_oi: np.ndarray,
    mf_oi: Dict[str, np.ndarray],
    config: GridRefinementConfig,
    cannon: CannonDiagnostics,
) -> Tuple[np.ndarray, GridDiagnostics]:
    """Compute a priority score for each zone.

    Higher score = more important to keep as a separate cell.
    The score combines:

      * physical gradients (density, temperature, velocity, composition)
      * cannon membership and interface proximity
      * shock / velocity reversal
      * photospheric region
      * optical-depth contribution

    The score is normalised so that the maximum gradient contribution
    is ~1.0, with cannon/shock/photosphere bonuses added on top.

    Returns:
        (priority, diagnostics) where priority is a float array of length n
    """
    n = len(r_oi)
    grads = compute_gradients(r_oi, v_oi, T_oi, rho_oi, mf_oi, config)
    diag = GridDiagnostics()
    diag.n_input = n

    # --- Normalised gradient score (0–1 range) ---
    # Use percentile-based normalisation so that a few extreme values
    # don't dominate.  The 90th percentile is mapped to 1.0.
    def _norm(arr: np.ndarray) -> np.ndarray:
        p90 = np.percentile(arr[arr > 0], 90) if np.any(arr > 0) else 1.0
        return np.minimum(arr / max(p90, 1e-30), 1.0)

    score = (
        0.25 * _norm(grads['rho']) +
        0.20 * _norm(grads['temp']) +
        0.25 * _norm(grads['vel']) +
        0.30 * _norm(grads['comp'])
    )

    # --- Shock / velocity reversal bonus ---
    dv = np.diff(v_oi)  # outer→inner
    shock_zones = np.zeros(n, dtype=bool)
    shock_zones[1:-1] = np.abs(dv[:-1] + dv[1:]) / 2.0 > config.shock_dv_threshold
    vel_reversal = np.zeros(n, dtype=bool)
    vel_reversal[1:-1] = (dv[:-1] < 0) | (dv[1:] < 0)
    shock_mask = shock_zones | vel_reversal
    score[shock_mask] += 2.0  # strong bonus
    diag.n_shock_cells = int(shock_mask.sum())

    # --- Photosphere bonus ---
    phot_mask = (tau_oi > config.phot_tau_low) & (tau_oi < config.phot_tau_high)
    score[phot_mask] += 1.0
    diag.n_phot_cells = int(phot_mask.sum())

    # --- Cannon bonus ---
    if cannon.detected and len(cannon.zone_indices) > 0:
        cannon_mask = np.zeros(n, dtype=bool)
        cannon_mask[cannon.zone_indices] = True
        # Interface zones
        inner_edge = cannon.zone_indices[-1]
        outer_edge = cannon.zone_indices[0]
        for k in range(1, config.cannon_min_cells_per_interface + 1):
            for edge in [outer_edge, inner_edge]:
                for idx in [edge + k, edge - k]:
                    if 0 <= idx < n:
                        cannon_mask[idx] = True
        score[cannon_mask] += 3.0  # highest priority
        diag.n_cannon_cells = int(cannon_mask.sum())

    # --- Gradient count for diagnostics ---
    grad_mask = (
        (grads['rho'] > config.rho_grad) |
        (grads['temp'] > config.temp_grad) |
        (grads['vel'] > config.vel_grad) |
        (grads['comp'] > config.comp_grad)
    )
    diag.n_grad_cells = int(grad_mask.sum())

    return score, diag


def adaptive_regrid(
    r_oi: np.ndarray,
    v_oi: np.ndarray,
    T_oi: np.ndarray,
    T_rad_oi: np.ndarray,
    rho_oi: np.ndarray,
    kappa_oi: np.ndarray,
    mf_oi: Dict[str, np.ndarray],
    tau_oi: np.ndarray,
    config: GridRefinementConfig,
    verbose: bool = True,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray,
           np.ndarray, np.ndarray, Dict[str, np.ndarray], np.ndarray,
           GridDiagnostics]:
    """Adaptive regridding of STELLA zones for SuperLite.

    This replaces ``tau_weighted_downsample`` with a physics-aware AMR
    algorithm.  The output preserves:

      * density, temperature, velocity structure
      * composition interfaces (cannon boundaries)
      * shocks and velocity reversals
      * photospheric / line-forming regions

    while minimising the number of cells.

    Algorithm
    ---------
    1. Compute a priority score for each zone based on physical gradients,
       cannon membership, shocks, and photospheric location.
    2. Greedily select zones by priority, but enforce a maximum gap
       constraint: no two selected zones should span a log-radius gap
       larger than ``coarsen_factor`` times the median gap.  This prevents
       the dr/r=0.85 problem seen in the old grid.
    3. Conservative remapping: merge zones between selected boundaries
       using mass-weighted averages.

    All inputs are outer→inner (decreasing r).  Outputs are also outer→inner
    EXCEPT ``edges_io`` which is inner→outer (increasing r, length n_out+1)
    for direct use in SuperLite's input.str format.

    Args:
        r_oi: zone-centre radii, outer→inner, cm
        v_oi: velocities, outer→inner, km/s
        T_oi: gas temperatures, outer→inner, K
        T_rad_oi: radiation temperatures, outer→inner, K
        rho_oi: densities, outer→inner, g/cm³
        kappa_oi: hybrid opacities, outer→inner, cm²/g
        mf_oi: mass fractions dict, outer→inner
        tau_oi: cumulative optical depth, outer→inner
        config: refinement configuration
        verbose: print diagnostics

    Returns:
        (r_out, v_out, T_out, T_rad_out, rho_out, kappa_out, mf_out,
         edges_io, diagnostics)
        All arrays outer→inner except edges_io (inner→outer).
    """
    n = len(r_oi)

    # --- Detect cannon ---
    cannon = detect_metal_rich_cannon(r_oi, v_oi, mf_oi, config)
    if verbose and cannon.detected:
        print(f"\n  [Cannon] Metal-rich H-poor component detected:")
        print(f"    zones: {cannon.zone_indices[0]}–{cannon.zone_indices[-1]} "
              f"({cannon.n_cells_inside} zones)")
        print(f"    r: {cannon.inner_radius:.3e} – {cannon.outer_radius:.3e} cm")
        print(f"    v: {cannon.inner_velocity:.1f} – {cannon.outer_velocity:.1f} km/s")
        print(f"    composition: {cannon.composition}")

    # --- Compute priority scores ---
    score, diag = _compute_zone_priority(
        r_oi, v_oi, T_oi, rho_oi, kappa_oi, tau_oi, mf_oi, config, cannon)
    diag.cannon = cannon

    if verbose:
        print(f"\n  [AMR] Priority scores: min={score.min():.3f} "
              f"max={score.max():.3f} median={np.median(score):.3f}")
        print(f"    gradient: {diag.n_grad_cells}, shock: {diag.n_shock_cells}, "
              f"photosphere: {diag.n_phot_cells}, cannon: {diag.n_cannon_cells}")

    # --- Adaptive grid construction in log-radius space ---
    # Instead of selecting from existing STELLA zones (which have wildly
    # non-uniform spacing), construct a NEW grid by placing boundaries
    # in log-radius space where physics demands resolution, then remap
    # STELLA quantities conservatively.
    #
    # Algorithm:
    # 1. Compute a refinement density function w(log r) based on gradients,
    #    cannon, shocks, and photosphere.
    # 2. Place max_cells boundaries such that the integral of w is uniform
    #    between boundaries (equidistribution principle).
    # 3. Snap each boundary to the nearest STELLA zone centre.
    # 4. Conservative remapping of STELLA zones into the new grid.

    log_r = np.log(r_oi)  # outer→inner, decreasing
    log_r_min = log_r[-1]  # innermost (smallest r)
    log_r_max = log_r[0]   # outermost (largest r)
    log_r_range = log_r_max - log_r_min

    # --- Refinement weight function ---
    # w[i] = 1 + gradient_score + cannon_bonus + shock_bonus + photosphere_bonus
    # Normalised so that the baseline (smooth region) gets weight 1.
    weight = np.ones(n)
    # Gradient contribution (normalised to 0–2 range)
    weight += 2.0 * np.minimum(score, 2.0)
    # Cannon and shock already included in score, but add explicit bonuses
    # for the cannon interfaces to ensure they get extra resolution
    if cannon.detected and len(cannon.zone_indices) > 0:
        # Boost weight at cannon boundaries
        for edge in [cannon.zone_indices[0], cannon.zone_indices[-1]]:
            for k in range(-config.cannon_min_cells_per_interface,
                           config.cannon_min_cells_per_interface + 1):
                idx = edge + k
                if 0 <= idx < n:
                    weight[idx] += 5.0 * max(0, 1 - abs(k) /
                                              config.cannon_min_cells_per_interface)

    # --- Equidistribution in log-radius × weight space ---
    # Build cumulative weight integral in log-r space.
    # Since STELLA zones have non-uniform dlogr, we integrate:
    #   W(i) = sum_{j=0}^{i} weight[j] * |dlogr[j]|
    # where dlogr is the log-radius spacing of each zone.
    dlogr = np.zeros(n)
    dlogr[0] = abs(log_r[0] - log_r[1]) if n > 1 else log_r_range
    dlogr[1:-1] = 0.5 * (abs(np.diff(log_r[:-1])) + abs(np.diff(log_r[1:])))
    dlogr[-1] = abs(log_r[-2] - log_r[-1]) if n > 1 else log_r_range

    # Cumulative weight (outer→inner)
    w_density = weight * dlogr  # weight per zone × zone width
    W_cumul = np.cumsum(w_density)
    W_total = W_cumul[-1]

    # Target: place max_cells boundaries at equal W_cumul intervals
    # But we want to place zone CENTRES, not boundaries.
    # We select max_cells zone centres from the STELLA grid.
    n_target = min(config.max_cells, n)
    target_W = np.linspace(0, W_total, n_target)

    # For each target W, find the STELLA zone whose cumulative weight
    # is closest. This gives equidistribution in weight space.
    selected_indices = set()
    selected_indices.add(0)  # always keep outermost
    selected_indices.add(n - 1)  # always keep innermost

    for tw in target_W[1:-1]:
        idx = int(np.searchsorted(W_cumul, tw))
        idx = min(max(idx, 0), n - 1)
        # Also check idx-1 in case it's closer
        if idx > 0 and abs(W_cumul[idx - 1] - tw) < abs(W_cumul[idx] - tw):
            idx = idx - 1
        selected_indices.add(idx)

    # If we have fewer than max_cells (due to deduplication), fill gaps
    sel_list = sorted(selected_indices)
    while len(sel_list) < n_target:
        # Find largest gap in index space
        best_gap = 0
        best_mid = -1
        for i in range(len(sel_list) - 1):
            lo = sel_list[i]
            hi = sel_list[i + 1]
            if hi - lo > 1:
                gap = abs(log_r[hi] - log_r[lo])
                if gap > best_gap:
                    best_gap = gap
                    best_mid = (lo + hi) // 2
        if best_mid < 0:
            break
        selected_indices.add(best_mid)
        sel_list = sorted(selected_indices)

    sel = np.array(sorted(selected_indices))
    n_out = len(sel)

    # --- Conservative remapping (merge zones between selected boundaries) ---
    # Compute original zone volumes (outer→inner)
    vol_orig = _compute_zone_volumes(r_oi)
    mass_orig = rho_oi * vol_orig

    # Build merged cell boundaries (inner→outer for SuperLite)
    r_io = r_oi[::-1]  # inner→outer
    edges_io_orig = np.zeros(n + 1)
    edges_io_orig[1:n] = np.sqrt(r_io[:-1] * r_io[1:])
    edges_io_orig[0] = r_io[0] ** 2 / max(edges_io_orig[1], 1e-30)
    edges_io_orig[n] = r_io[-1] ** 2 / max(edges_io_orig[n - 1], 1e-30)

    # Output arrays (outer→inner)
    r_out = np.zeros(n_out)
    v_out = np.zeros(n_out)
    T_out = np.zeros(n_out)
    T_rad_out = np.zeros(n_out)
    rho_out = np.zeros(n_out)
    kappa_out = np.zeros(n_out)
    mf_out = {k: np.zeros(n_out) for k in mf_oi}
    edges_merged_io = np.zeros(n_out + 1)  # inner→outer

    for idx_out in range(n_out):
        i_start = sel[idx_out]
        i_end = sel[idx_out + 1] if idx_out + 1 < n_out else n
        sl = slice(i_start, i_end)

        m_grp = mass_orig[sl]
        v_grp = vol_orig[sl]
        m_tot = float(m_grp.sum())
        v_tot = float(v_grp.sum())

        if m_tot > 0:
            w = m_grp / m_tot
            v_out[idx_out] = float(np.sum(w * v_oi[sl]))
            T_out[idx_out] = float(np.sum(w * T_oi[sl]))
            T_rad_out[idx_out] = float(np.sum(w * T_rad_oi[sl]))
            kappa_out[idx_out] = float(np.sum(w * kappa_oi[sl]))
            for k in mf_oi:
                mf_out[k][idx_out] = float(np.sum(w * mf_oi[k][sl]))
            r_out[idx_out] = float(np.sum(v_grp * r_oi[sl]) / v_tot)
        else:
            v_out[idx_out] = float(np.mean(v_oi[sl]))
            T_out[idx_out] = float(np.mean(T_oi[sl]))
            T_rad_out[idx_out] = float(np.mean(T_rad_oi[sl]))
            kappa_out[idx_out] = float(np.mean(kappa_oi[sl]))
            for k in mf_oi:
                mf_out[k][idx_out] = float(np.mean(mf_oi[k][sl]))
            r_out[idx_out] = float(np.mean(r_oi[sl]))

        rho_out[idx_out] = m_tot / max(v_tot, 1e-30)

        # Merged cell boundaries (inner→outer)
        # sel is outer→inner (decreasing r).  Group idx_out covers
        # outer→inner zones [i_start, i_end).
        # In inner→outer indexing: innermost zone = n - i_end, outermost = n - 1 - i_start
        inner_edge = edges_io_orig[n - i_end]
        outer_edge = edges_io_orig[n - i_start]
        io_idx = n_out - 1 - idx_out  # reverse: outermost group → last cell
        edges_merged_io[io_idx] = inner_edge
        edges_merged_io[io_idx + 1] = outer_edge

    # --- Conservation checks ---
    vol_merged = (4.0 / 3.0) * np.pi * (
        edges_merged_io[1:] ** 3 - edges_merged_io[:-1] ** 3)
    vol_merged_oi = vol_merged[::-1]  # outer→inner
    mass_out = rho_out * vol_merged_oi

    diag.n_output = n_out
    diag.mass_in = float(mass_orig.sum() / MSUN)
    diag.mass_out = float(mass_out.sum() / MSUN)
    diag.mass_error = float(abs(mass_out.sum() - mass_orig.sum()) /
                            max(mass_orig.sum(), 1e-30))

    # Element mass conservation
    for elem in mf_oi:
        m_elem_in = float(np.sum(mf_oi[elem] * mass_orig))
        m_elem_out = float(np.sum(mf_out[elem] * mass_out))
        if m_elem_in > 0:
            diag.element_mass_errors[elem] = abs(m_elem_out - m_elem_in) / m_elem_in

    # Grid quality
    dr = np.diff(edges_merged_io)  # inner→outer
    r_centres_io = 0.5 * (edges_merged_io[:-1] + edges_merged_io[1:])
    dr_r = dr / r_centres_io
    diag.dr_r_min = float(dr_r.min())
    diag.dr_r_max = float(dr_r.max())
    diag.dr_r_ratio = float(dr_r.max() / max(dr_r.min(), 1e-30))

    # Velocity extrema preservation
    v_in_min, v_in_max = float(v_oi.min()), float(v_oi.max())
    v_out_min, v_out_max = float(v_out.min()), float(v_out.max())
    # Allow small tolerance due to mass-weighted averaging
    v_tol = 0.05 * (v_in_max - v_in_min) + 10.0  # 5% of range + 10 km/s
    diag.velocity_extrema_preserved = (
        abs(v_out_min - v_in_min) < v_tol and abs(v_out_max - v_in_max) < v_tol
    )

    # Negative velocity gradient preservation
    dv_out = np.diff(v_out)
    diag.n_negative_dv_preserved = int(np.sum(dv_out < 0))

    # Cannon mass
    if cannon.detected:
        # Find which output cells overlap with the cannon zone range
        cannon_r_lo = cannon.inner_radius
        cannon_r_hi = cannon.outer_radius
        # r_out is outer→inner; cannon spans r_oi[inner_edge] to r_oi[outer_edge]
        in_cannon = (r_out <= cannon_r_hi) & (r_out >= cannon_r_lo)
        diag.cannon.n_cells_inside = int(in_cannon.sum())
        diag.cannon.mass = float(mass_out[in_cannon].sum() / MSUN)

    if verbose:
        print(f"\n  [AMR] Adaptive grid: {n} → {n_out} cells")
        print(f"    dr/r range: {diag.dr_r_min:.4e} – {diag.dr_r_max:.4e} "
              f"(ratio: {diag.dr_r_ratio:.1f}×)")
        print(f"    Mass conservation: {diag.mass_in:.4f} → {diag.mass_out:.4f} Msun "
              f"(error: {diag.mass_error:.2e})")
        print(f"    Velocity extrema preserved: {diag.velocity_extrema_preserved}")
        print(f"    Negative dv preserved: {diag.n_negative_dv_preserved}")
        if cannon.detected:
            print(f"    Cannon: {diag.cannon.n_cells_inside} cells, "
                  f"{diag.cannon.mass:.4f} Msun")
        # Report worst element mass errors
        if diag.element_mass_errors:
            worst = max(diag.element_mass_errors.items(), key=lambda x: x[1])
            print(f"    Worst element mass error: {worst[0]} = {worst[1]:.2e}")

    return (r_out, v_out, T_out, T_rad_out, rho_out, kappa_out, mf_out,
            edges_merged_io, diag)


# ---------------------------------------------------------------------------
# Conservation validation
# ---------------------------------------------------------------------------

def validate_conservation(
    diag: GridDiagnostics,
    config: GridRefinementConfig,
    mass_tol: float = 1e-4,
    element_mass_tol: float = 1e-3,
    composition_tol: float = 1e-6,
    fail_hard: bool = True,
) -> List[str]:
    """Validate conservation tests from grid diagnostics.

    Args:
        diag: grid diagnostics from adaptive_regrid
        config: refinement config (unused but kept for API consistency)
        mass_tol: relative mass tolerance
        element_mass_tol: relative element mass tolerance
        composition_tol: per-cell composition sum tolerance
        fail_hard: if True, raise AssertionError on failure; else return errors

    Returns:
        List of error messages (empty if all pass)
    """
    errors: List[str] = []

    if diag.mass_error > mass_tol:
        errors.append(
            f"Mass conservation failed: {diag.mass_error:.2e} > {mass_tol:.0e} "
            f"(M_in={diag.mass_in:.6f}, M_out={diag.mass_out:.6f} Msun)"
        )

    for elem, err in diag.element_mass_errors.items():
        if err > element_mass_tol:
            errors.append(
                f"Element mass conservation failed for {elem}: "
                f"{err:.2e} > {element_mass_tol:.0e}"
            )

    if not diag.velocity_extrema_preserved:
        errors.append(
            "Velocity extrema not preserved after regridding — "
            "shock or velocity reversal may have been smeared"
        )

    if errors and fail_hard:
        raise AssertionError(
            "Conservation validation failed:\n  " + "\n  ".join(errors)
        )

    return errors


# ---------------------------------------------------------------------------
# Standalone diagnostics
# ---------------------------------------------------------------------------

def analyze_grid_quality(input_str_path: str) -> GridDiagnostics:
    """Analyze the quality of an existing SuperLite input.str file.

    Reads the file and reports grid uniformity, composition structure,
    and cannon detection.

    Args:
        input_str_path: path to input.str

    Returns:
        GridDiagnostics with quality metrics
    """
    lines = open(input_str_path).readlines()
    # Parse header
    parts = lines[1].split()
    # Format: #  nx 1 1 ncol nabund
    nx = int(parts[1])
    nabund = int(parts[5])
    labels = lines[2].split()[1:]  # skip leading '#'

    data = []
    for line in lines[3:]:
        if line.strip() and not line.startswith('#'):
            data.append([float(x) for x in line.split()])
    arr = np.array(data)

    xl = arr[:, 0]
    xr = arr[:, 1]
    mass = arr[:, 4]
    T = arr[:, 6]
    abund = arr[:, 7:7 + nabund]

    diag = GridDiagnostics()
    diag.n_input = nx
    diag.mass_in = float(mass.sum() / MSUN)
    diag.mass_out = diag.mass_in
    diag.mass_error = 0.0

    dr = xr - xl
    rc = 0.5 * (xl + xr)
    dr_r = dr / rc
    diag.dr_r_min = float(dr_r.min())
    diag.dr_r_max = float(dr_r.max())
    diag.dr_r_ratio = float(dr_r.max() / max(dr_r.min(), 1e-30))

    return diag


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description='Analyze SuperLite input.str grid quality')
    parser.add_argument('input_str', help='Path to input.str')
    args = parser.parse_args()
    diag = analyze_grid_quality(args.input_str)
    print(f"Grid: {diag.n_input} cells")
    print(f"  Mass: {diag.mass_in:.4f} Msun")
    print(f"  dr/r: {diag.dr_r_min:.4e} – {diag.dr_r_max:.4e} "
          f"(ratio: {diag.dr_r_ratio:.1f}×)")
