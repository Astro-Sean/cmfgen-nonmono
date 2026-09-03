"""Configuration for the X_newgrid hierarchical grid constructor.

All thresholds are *starting points for convergence testing*, not universal
CMFGEN prescriptions (see project AGENTS.md evidence policy).  The numerical
targets follow the Dessart/Hillier literature surveyed for this experiment:

  * ~5-7 points per optical-depth decade  (default 6)
  * typical total ND ~80-200             (search range, not a fixed target)
  * maximum ND = 1000                     (hard CMFGEN ceiling, NOT a target)
  * velocity extrema / composition boundaries / shocks = mandatory points
  * explicit priority hierarchy for refinement:
        tau > velocity > discontinuities > ionisation fronts
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional


@dataclass
class XNewGridConfig:
    """All tunable parameters for the X_newgrid grid constructor.

    Priority hierarchy (Gap A)
    --------------------------
    Refinement chooses the worst-violating cell under an ordered hierarchy,
    not a weighted sum.  A higher-priority violation is always addressed
    before any lower-priority one, regardless of magnitude.  Within a single
    priority level the normalised violation magnitude breaks ties.

    The levels, in descending priority:

        TAU          - optical-depth spacing (primary transfer coordinate)
        VELOCITY     - |dv|, |dv|/v, and dv/dr extrema (frequency coupling)
        DISCONTINUITY- shocks, density jumps, composition interfaces
        ION_FRONT    - ionisation/recombination fronts (Saha estimate on
                       the STELLA source; real CMFGEN fronts handled by the
                       standalone two-pass script)
        SMOOTH       - residual T / rho / geometric smoothness
    """

    # --- Point budget ---
    max_points: int = 1000          # HARD CMFGEN ceiling (not a target)
    min_points: int = 60            # adaptive stop floor
    target_points: int = 150        # nominal search centre (80-200 range)
    points_per_tau_decade: float = 6.0   # Gap B: skeleton sizing
    skeleton_min: int = 20          # minimum skeleton points before features

    # --- Optical depth (Priority 1: TAU) ---
    tau_floor: float = 1.0e-6       # log(tau+floor) guard for thin zones
    target_dlogtau: float = 0.16    # max dlog10(tau) per cell (~1/6 decade)
    tau_max: float = 100.0          # inner optical-depth trim cutoff
    tau_min: float = 1.0e-3         # outer optical-depth trim cutoff
    kappa_source: str = "hybrid"    # 'stella' | 'hybrid' | 'kramers'

    # --- Velocity (Priority 2: VELOCITY) ---
    target_dv_kms: float = 300.0    # max |dv| per cell (km/s)
    target_dv_over_v: float = 0.30  # max |dv|/v per cell
    velocity_extremum_window: int = 5
    min_points_per_extremum: int = 5
    # dv/dr sign-change (extremum) locations are MANDATORY grid points.
    extremum_significance_kms: float = 50.0  # ignore flat noise extrema

    # --- Discontinuities (Priority 3: DISCONTINUITY) ---
    composition_jump: float = 0.30          # |dX_i| > this => interface
    composition_trace_floor: float = 1.0e-4
    min_points_per_interface: int = 5
    shock_density_jump: float = 1.5         # rho compression ratio
    shock_velocity_curv: float = 0.5        # normalised 2nd difference of v
    shock_temp_jump: float = 0.5            # dlogT
    min_points_per_shock: int = 7
    target_dlogrho: float = 0.30            # max dlog10(rho) per cell
    max_density_jump: float = 0.15          # fractional density jump per cell

    # --- Ionisation fronts (Priority 4: ION_FRONT) ---
    # On the STELLA source these are Saha-estimated; the real CMFGEN fronts
    # are refined by the standalone two-pass script.
    ne_log_threshold: float = 1.0           # |dlog10(Ne)| => front
    ye_log_threshold: float = 0.5           # |dlog10(Ye)| => front
    ionfront_t_log_threshold: float = 0.3
    ionfront_kappa_log_threshold: float = 0.5
    min_points_per_ionfront: int = 5

    # --- Smoothness (Priority 5: SMOOTH) ---
    target_dlogT: float = 0.20              # max dlog10(T) per cell
    max_temp_jump: float = 0.10             # fractional T jump per cell
    target_dr_over_r: float = 0.10          # max geometric stretch
    max_cell_ratio: float = 1.5             # max adjacent cell size ratio
    min_dr_r: float = 1.0e-6                # prevent degenerate cells

    # --- Coarsening (Gap E) ---
    enable_coarsening: bool = True
    coarsen_max_passes: int = 3
    # Hard limit = hard_factor * soft target; coarsening must stay within it.
    coarsen_hard_factor: float = 3.0

    # --- Local interpolation error (Gap D) ---
    # Accept a coarsened grid only if the max-local interpolation error
    # against the original STELLA structure stays below these per variable.
    error_tol_dlogtau: float = 0.05
    error_tol_dlogrho: float = 0.10
    error_tol_dv_frac: float = 0.05
    error_tol_dlogT: float = 0.05
    error_tol_dX: float = 0.02
    # Only score / measure error in the radiatively-important region
    # (tau > tau_metric_min).  The optically-thin outer region is covered by
    # the skeleton alone and is not over-refined; zero-mass cells there are
    # numerical artifacts, not physics errors.
    tau_metric_min: float = 1.0e-2

    # --- Conservation ---
    conserve_mass: bool = True
    mass_tol: float = 1.0e-3

    # --- Solver / output ---
    vel_law: int = 12
    preserve_v_sign: bool = True
    t_floor: float = 3000.0
    t_cap: float = 0.0
    trace_threshold: float = 1.0e-4

    # --- Diagnostics ---
    verbose: bool = True
    plot: bool = False

    # --- Velocity trimming (0 = disabled) ---
    v_min_cms: float = 0.0
    v_max_cms: float = 0.0

    def priority_order(self) -> tuple[str, ...]:
        """Return the refinement priority levels in descending priority."""
        return ("TAU", "VELOCITY", "DISCONTINUITY", "ION_FRONT", "SMOOTH")
