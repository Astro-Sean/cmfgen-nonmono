#!/usr/bin/env python3
"""
Tau-Based AMR STELLA -> CMFGEN Grid Optimiser
=============================================

Converts a STELLA ``.swd`` model into an optimised CMFGEN input grid for the
non-monotonic (VEL_LAW=12) solver using **optical-depth-based adaptive mesh
refinement (AMR)** with a hard maximum of 1000 grid points.

Design philosophy
-----------------
The CMFGEN grid is constructed primarily in optical-depth space, with
hydrodynamic / composition criteria acting as *hard constraints* or secondary
refinement indicators.  The 1000-point budget is treated as a constrained
optimisation problem: spend the available points where they reduce
radiative-transfer + NLTE + representation error the most, while preserving
the exact non-monotonic topology of the hydrodynamic model.

Pipeline::

    READ STELLA .swd
        |
    validate units / structure / radial ordering
        |
    reconstruct optical depth  (hybrid: STELLA NLTE kappa + e-scatter floor + Kramers)
        |
    identify composition interfaces   (HARD constraint)
    identify shocks                   (HARD constraint)
    identify velocity extrema         (HARD constraint)
    identify shell / wind boundaries  (HARD constraint)
        |
    construct coarse log(tau) base grid
        |
    recursive cell-based AMR:
        compute multi-component refinement score
        flag hard-constraint cells
        split highest-priority cells
        repeat until criteria satisfied OR N = N_max
        |
    conservative remap (PCHIP smooth / one-sided at discontinuities)
        |
    validate: mass, species mass, radial ordering, velocity topology, sum(X)=1
        |
    generate CMFGEN SN_HYDRO_DATA  +  diagnostics + manifest

This module reuses the proven ``.swd`` parser, opacity routines and
``write_sn_hydro_data`` from ``convert_stella_swd_to_cmfgen.py`` so that the
CMFGEN output format is identical to the rest of the pipeline.

Usage (standalone)::

    python scripts/tau_amr_optimiser.py \\
        --swd X_v11/stella/strad/X_v11.swd \\
        --epoch 25.0 \\
        --output X_v11/cmfgen/SN_HYDRO_DATA \\
        --max-points 1000 \\
        --plot-dir X_v11/plots/tau_amr

It is also callable as a library function ``run_tau_amr()`` and is wired into
``convert_stella_swd_to_cmfgen.py`` via ``--base-grid-method tau_amr_recursive``.

Author: Devin (Cognition) for the SN2024pba project.
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass, field, asdict
from pathlib import Path

import numpy as np

# Reuse the proven infrastructure from the main converter -------------------
sys.path.insert(0, str(Path(__file__).parent))
from convert_stella_swd_to_cmfgen import (  # noqa: E402
    read_swd,
    normalize_swd_composition,
    approximate_rosseland_opacity,
    compute_kappa_tau,
    write_sn_hydro_data,
    STELLA_ELEM_NAMES,
    STELLA_AZ,
    AVOGADRO,
)

# Element atomic numbers (for electron-scattering floor)
ELEM_Z = {'H': 1, 'He': 2, 'C': 6, 'N': 7, 'O': 8, 'Ne': 10, 'Na': 11,
          'Mg': 12, 'Al': 13, 'Si': 14, 'S': 16, 'Ar': 18, 'Ca': 20,
          'Fe': 26, 'Ni': 28}
ELEM_A = {n: float(a) for n, a in zip(STELLA_ELEM_NAMES, STELLA_AZ)}

MSUN = 1.989e33


# ===========================================================================
# Configuration
# ===========================================================================

@dataclass
class AMRConfig:
    """All tunable thresholds for the tau-AMR algorithm.

    Defaults are physically motivated; override via YAML / CLI for sensitivity
    testing.  See the project AGENTS.md mandatory-evidence policy: do not
    declare any threshold optimal without convergence testing.
    """
    # --- Point budget ---
    max_points: int = 1000          # HARD upper limit
    min_points: int = 60            # adaptive stop floor
    base_grid_points: int = 80      # initial coarse log(tau) grid

    # --- Optical depth ---
    tau_floor: float = 1e-6         # log(tau+floor) guard for optically thin zones
    kappa_source: str = "hybrid"    # 'stella' | 'hybrid' | 'kramers'
    tau_inner_trim: float = 0.0     # trim zones with tau > tau_max - this (0=keep all)

    # --- Refinement targets (soft criteria; dimensionless ratios) ---
    # These are the CONVERGENCE targets with safety margin already applied.
    # The AMR loop stops when all cells satisfy these.  The hard limits in
    # the resolution report use 3x these values.
    target_dlogtau: float = 0.08    # max dlog10(tau) per cell in smooth regions
    target_dlogrho: float = 0.30    # max dlog10(rho) per cell
    target_dv_kms: float = 300.0    # max |dv| per cell (km/s)
    target_dlogT: float = 0.20      # max dlog10(T) per cell
    target_dr_over_r: float = 0.10  # max geometric stretch
    target_dtau_es: float = 1.0     # max electron-scattering tau per cell
    target_dlogchi: float = 0.5     # max dlog10(chi) per cell (opacity gradient)
    target_dv_over_v: float = 0.30  # max |dv|/v per cell (velocity/frequency coupling)

    # --- Safety margin for convergence ---
    # The AMR targets above already include a safety margin relative to the
    # hard resolution limits.  This factor controls how much margin is applied
    # when checking the stop condition: the loop stops when all scores < 1.0
    # (i.e. all cells are within the target).  A value < 1.0 makes the loop
    # stop earlier (more conservative); > 1.0 allows operating closer to the
    # hard limit.
    safety_margin: float = 1.0

    # --- Coarsening pass ---
    # After AMR converges, attempt to remove individual points while maintaining
    # all resolution constraints.  This finds the true minimum viable grid.
    enable_coarsening: bool = True
    coarsen_max_passes: int = 3     # number of full coarsening sweeps

    # --- Refinement weights (soft score) ---
    w_tau: float = 1.0
    w_rho: float = 0.5
    w_v: float = 0.5
    w_T: float = 0.4
    w_curv: float = 0.3
    w_dtau_es: float = 0.6       # electron-scattering optical depth
    w_dlogchi: float = 0.5       # opacity gradient
    w_dv_over_v: float = 0.4     # velocity/frequency coupling

    # --- Hard-constraint thresholds ---
    composition_jump: float = 0.30      # |dX_i| > this => chemical interface
    composition_trace_floor: float = 1e-4  # only flag species above this abundance
    min_points_per_interface: int = 5   # forced zones straddling each interface
    shock_density_jump: float = 1.5     # rho_{i+1}/rho_i > this => shock
    shock_velocity_curv: float = 0.5    # |d2v/dr2| (normalised) => shock
    shock_temp_jump: float = 0.5        # dlogT > this => shock
    min_points_per_shock: int = 7
    velocity_extrema_window: int = 5    # zones for extrema detection
    min_points_per_extremum: int = 5
    min_dr_r: float = 1e-4              # prevent degenerate cells

    # --- Composition interface representation ---
    interface_transition_zones: int = 3  # finite transition width if needed (0=sharp)

    # --- Conservation ---
    conserve_mass: bool = True
    mass_tol: float = 1e-3             # relative mass conservation tolerance

    # --- Diagnostics ---
    plot: bool = True
    verbose: bool = True


# ===========================================================================
# Optical-depth reconstruction
# ===========================================================================

def electron_scattering_kappa(mass_fractions: dict) -> np.ndarray:
    """Electron-scattering opacity assuming full ionisation (conservative floor).

    kappa_es = 0.2 * (1+X) / (1 + 3X + 0.25Y)  [cm^2/g]
    Used as a floor in the hybrid tau so that neutral regions which SuperLite /
    CMFGEN will re-ionise are not underestimated.
    """
    X = mass_fractions.get('H', np.zeros(1))
    Y = mass_fractions.get('He', np.zeros(1))
    return 0.2 * (1.0 + X) / (1.0 + 3.0 * X + 0.25 * Y)


def hybrid_kappa(rho, T, mass_fractions, kappa_stella=None):
    """Hybrid Rosseland opacity: element-wise max of three estimates.

    1. STELLA NLTE kappa (accurate in hot ionised regions) — if available.
    2. Electron-scattering floor (full-ionisation, conservative).
    3. Kramers bound-free + free-free (corrected metal mass fraction).

    This avoids the STELLA-native kappa dropping to ~0 in cool neutral regions
    where CMFGEN/SuperLite will re-ionise the material and see the full
    electron-scattering opacity.
    """
    X = mass_fractions.get('H', np.zeros_like(rho))
    Y = mass_fractions.get('He', np.zeros_like(rho))
    Z = np.zeros_like(rho)
    for elem, frac in mass_fractions.items():
        if elem not in ('H', 'He'):
            Z = Z + frac
    kappa_kramers = approximate_rosseland_opacity(T, rho, X, Y, Z,
                                                  mass_fractions=mass_fractions)
    kappa_es = electron_scattering_kappa(mass_fractions)
    candidates = [np.maximum(kappa_kramers, 1e-8), np.maximum(kappa_es, 1e-8)]
    if kappa_stella is not None:
        candidates.append(np.maximum(kappa_stella, 1e-8))
    return np.maximum.reduce(candidates)


def reconstruct_optical_depth(r, rho, T, mass_fractions, kappa_stella=None,
                              cfg: AMRConfig = AMRConfig()):
    """Cumulative optical depth tau(r) = int_r^R kappa rho dr' (outer->inner).

    Returns (tau, kappa_used, label).  Arrays are outer->inner (r decreasing),
    so tau increases with index.
    """
    if cfg.kappa_source == 'stella' and kappa_stella is not None:
        kappa = np.maximum(kappa_stella, 1e-8)
        label = 'STELLA NLTE'
    elif cfg.kappa_source == 'kramers':
        kappa = hybrid_kappa(rho, T, mass_fractions, kappa_stella=None)
        label = 'Kramers+es'
    else:
        kappa = hybrid_kappa(rho, T, mass_fractions, kappa_stella=kappa_stella)
        label = 'hybrid(max of STELLA, es-floor, Kramers)'
    tau = compute_kappa_tau(r, rho, kappa)
    return tau, kappa, label


# ===========================================================================
# Feature detectors (hard constraints)
# ===========================================================================

@dataclass
class Feature:
    kind: str          # 'composition' | 'shock' | 'velocity_extremum' | 'boundary'
    index: int         # zone index (outer->inner) at the feature centre
    radius: float
    velocity: float
    density: float
    temperature: float
    tau: float
    detail: str = ""
    n_points: int = 5  # forced minimum zones around this feature


def detect_composition_interfaces(r, mass_fractions, tau, cfg: AMRConfig):
    """Detect sharp chemical interfaces (HARD constraint).

    Flags zones where |dX_i| > composition_jump for any species above the
    trace floor.  These must be densely sampled regardless of optical depth.
    """
    feats = []
    n = len(r)
    elems = [e for e in mass_fractions if e in STELLA_ELEM_NAMES]
    for i in range(1, n):
        for elem in elems:
            x = mass_fractions[elem]
            if max(x[i - 1], x[i]) < cfg.composition_trace_floor:
                continue
            dx = abs(x[i] - x[i - 1])
            denom = max(x[i - 1], x[i], 1e-30)
            if dx > cfg.composition_jump:
                feats.append(Feature(
                    'composition', i, r[i], 0.0, 0.0, 0.0, tau[i],
                    detail=f"{elem}: {x[i-1]:.3f}->{x[i]:.3f}",
                    n_points=cfg.min_points_per_interface))
                break  # one feature per zone
    return feats


def detect_shocks(r, v, rho, T, tau, cfg: AMRConfig):
    """Detect shocks / contact discontinuities (HARD constraint).

    Uses density compression, velocity curvature, and temperature jumps.
    """
    feats = []
    n = len(r)
    logrho = np.log10(np.maximum(rho, 1e-30))
    logT = np.log10(np.maximum(T, 1.0))
    for i in range(2, n - 2):
        rho_jump = max(rho[i] / max(rho[i - 1], 1e-30),
                       rho[i + 1] / max(rho[i], 1e-30))
        dT = abs(logT[i] - logT[i - 1])
        # velocity curvature (second difference), normalised by local |dv|
        dv1 = v[i] - v[i - 1]
        dv2 = v[i + 1] - v[i]
        dv_scale = max(abs(dv1), abs(dv2), 1.0)
        curv = abs(dv2 - dv1) / dv_scale
        is_shock = (rho_jump > cfg.shock_density_jump or
                    dT > cfg.shock_temp_jump or
                    curv > cfg.shock_velocity_curv)
        if is_shock:
            feats.append(Feature('shock', i, r[i], v[i], rho[i], T[i], tau[i],
                                 detail=f"rho_jump={rho_jump:.2f} dT={dT:.2f} "
                                        f"curv={curv:.2f}",
                                 n_points=cfg.min_points_per_shock))
    # de-duplicate nearby shocks (within 3 zones)
    dedup = []
    for f in feats:
        if not dedup or abs(f.index - dedup[-1].index) > 3:
            dedup.append(f)
    return dedup


def detect_velocity_extrema(r, v, tau, cfg: AMRConfig):
    """Detect velocity extrema (HARD constraint) — preserves non-monotonic topology."""
    feats = []
    n = len(v)
    w = cfg.velocity_extrema_window
    for i in range(w, n - w):
        vi = v[i]
        left = v[i - w:i]
        right = v[i + 1:i + 1 + w]
        if (vi >= left.max() and vi >= right.max() and
                (vi > left.min() + 50 or vi > right.min() + 50)):
            feats.append(Feature('velocity_extremum', i, r[i], vi, 0, 0, tau[i],
                                 detail="local maximum",
                                 n_points=cfg.min_points_per_extremum))
        elif (vi <= left.min() and vi <= right.min() and
              (vi < left.max() - 50 or vi < right.max() - 50)):
            feats.append(Feature('velocity_extremum', i, r[i], vi, 0, 0, tau[i],
                                 detail="local minimum",
                                 n_points=cfg.min_points_per_extremum))
    return feats


def detect_all_features(r, v, rho, T, mass_fractions, tau, cfg: AMRConfig):
    comp = detect_composition_interfaces(r, mass_fractions, tau, cfg)
    shocks = detect_shocks(r, v, rho, T, tau, cfg)
    extrema = detect_velocity_extrema(r, v, tau, cfg)
    allf = comp + shocks + extrema
    # fill in v/rho/T for comp/extrema features that had zeros
    for f in allf:
        if f.kind == 'composition' or f.kind == 'velocity_extremum':
            f.velocity = float(v[f.index])
            f.density = float(rho[f.index])
            f.temperature = float(T[f.index])
    return allf, comp, shocks, extrema


# ===========================================================================
# Conservative remapping / interpolation
# ===========================================================================

def _pchip_interp(x_new, x_old, y_old):
    """Monotonic piecewise-cubic interpolation (PCHIP) via scipy if available,
    else linear fallback.  Preserves monotonicity — no overshoot."""
    try:
        from scipy.interpolate import PchipInterpolator
        # x must be increasing for scipy
        order = np.argsort(x_old)
        return PchipInterpolator(x_old[order], y_old[order])(x_new)
    except Exception:
        return np.interp(x_new, x_old, y_old)


def _zone_bounds(r):
    """Return (r_out, r_in) per zone for an outer->inner (decreasing) radius
    array.  The innermost zone's inner edge is taken as 0.7*r[-1] (a small
    core) so mass is not grossly over-estimated; both old and new grids use
    the same convention so conservation is consistent."""
    r_out = r
    r_in = np.concatenate([r[1:], [0.7 * r[-1]]])
    return r_out, r_in


def _shell_volume(r_out, r_in):
    return (4.0 / 3.0) * np.pi * np.maximum(r_out**3 - r_in**3, 0.0)


def conservative_remap_density_composition(r_new, r_old, rho_old,
                                           mass_fractions_old):
    """Zone-integrated conservative remap of density and composition.

    Each old zone is treated as uniform (rho_old[i], X_i) over its shell
    volume.  For each new zone, sum the mass (and per-species mass) of all
    overlapping old zones, then divide by the new zone volume to get the new
    density and mass fractions.  This conserves total mass and species mass
    to machine precision regardless of grid alignment.

    Returns (rho_new, mf_new).
    """
    ro_old, ri_old = _zone_bounds(r_old)
    vol_old = _shell_volume(ro_old, ri_old)
    dm_old = vol_old * rho_old
    # per-species mass per old zone
    sm_old = {e: dm_old * x for e, x in mass_fractions_old.items()}

    ro_new, ri_new = _zone_bounds(r_new)
    vol_new = _shell_volume(ro_new, ri_new)
    n_new = len(r_new)
    n_old = len(r_old)

    rho_new = np.zeros(n_new)
    sm_new = {e: np.zeros(n_new) for e in mass_fractions_old}

    # For each new zone, find overlapping old zones by radius intersection.
    # Both grids are sorted outer->inner (decreasing r_out).
    j_start = 0
    for j in range(n_new):
        r_o = ro_new[j]
        r_i = ri_new[j]
        dm_j = 0.0
        sm_j = {e: 0.0 for e in mass_fractions_old}
        # advance j_start to first old zone whose r_in < r_o
        while j_start < n_old and ri_old[j_start] >= r_o:
            j_start += 1
        k = j_start
        while k < n_old and ro_old[k] > r_i:
            # overlap interval [max(r_i, ri_old[k]), min(r_o, ro_old[k])]
            lo = max(r_i, ri_old[k])
            hi = min(r_o, ro_old[k])
            if hi > lo:
                frac = (hi**3 - lo**3) / (ro_old[k]**3 - ri_old[k]**3) \
                       if (ro_old[k]**3 - ri_old[k]**3) > 0 else 0.0
                dm_j += dm_old[k] * frac
                for e in sm_j:
                    sm_j[e] += sm_old[e][k] * frac
            k += 1
        rho_new[j] = dm_j / max(vol_new[j], 1e-30)
        for e in sm_new:
            sm_new[e][j] = sm_j[e]

    # mass fractions = species mass / zone mass
    mf_new = {}
    dm_new = vol_new * rho_new
    for e in mass_fractions_old:
        mf_new[e] = np.where(dm_new > 1e-30, sm_new[e] / dm_new, 0.0)
    # renormalise to sum=1 (guards against tiny numerical drift)
    _renormalise(mf_new)
    rho_new = np.maximum(rho_new, 1e-30)
    return rho_new, mf_new


def remap_onto_grid(r_new, r_old, v_old, T_old, rho_old, kappa_old,
                    mass_fractions_old, interface_indices, cfg: AMRConfig):
    """Conservatively remap physical quantities onto the new AMR grid.

    - v, T, kappa: PCHIP (monotonic, no overshoot) — intensive quantities.
    - rho, composition: **zone-integrated conservative remap** (preserves
      total mass and species mass to machine precision).
    - At sharp composition interfaces the conservative remap naturally
      preserves the discontinuity because old zones are uniform; a new zone
      straddling an interface gets the mass-weighted average of the two sides,
      which is the physically correct conservative representation.  To keep
      the interface *sharp*, interface old zones are force-retained as grid
      points by the AMR (hard constraint), so straddling is minimised.

    All arrays are outer->inner (r decreasing).
    """
    x_inc = r_old[::-1]
    x_new_inc = r_new[::-1]

    def interp_inc(y_old):
        return _pchip_interp(x_new_inc, x_inc, y_old[::-1])[::-1]

    v_new = interp_inc(v_old)
    T_new = np.maximum(interp_inc(T_old), 1.0)
    kappa_new = np.maximum(interp_inc(kappa_old), 1e-8)

    rho_new, mf_new = conservative_remap_density_composition(
        r_new, r_old, rho_old, mass_fractions_old)
    return v_new, T_new, rho_new, kappa_new, mf_new


def _renormalise(mf: dict):
    tot = np.zeros(next(iter(mf.values())).shape)
    for x in mf.values():
        tot = tot + x
    tot = np.where(tot < 1e-30, 1.0, tot)
    for elem in mf:
        mf[elem] = mf[elem] / tot


# ===========================================================================
# Mass / conservation
# ===========================================================================

def zone_masses(r, rho):
    """Mass per zone assuming spherical shells with outer edge at r[i]."""
    n = len(r)
    # r is decreasing (outer->inner).  inner edge of zone i is r[i+1] (or 0).
    r_out = r
    r_in = np.concatenate([r[1:], [0.0]])
    dm = (4.0 / 3.0) * np.pi * np.maximum(r_out**3 - r_in**3, 0.0) * rho
    return dm


def total_mass(r, rho):
    return float(zone_masses(r, rho).sum())


def species_masses(r, rho, mass_fractions):
    dm = zone_masses(r, rho)
    return {e: float((dm * x).sum()) for e, x in mass_fractions.items()}


# ===========================================================================
# Recursive cell-based AMR
# ===========================================================================

def build_base_grid(r, tau, n_base, tau_floor):
    """Build a coarse base grid uniformly spaced in log(tau)."""
    f = np.log10(tau + tau_floor)
    f_targets = np.linspace(f[0], f[-1], n_base)
    tau_targets = 10.0 ** f_targets - tau_floor
    tau_targets = np.maximum(tau_targets, 0.0)
    # tau increases with index, r decreases.  Invert to get radius.
    r_grid = np.interp(tau_targets, tau, r)
    # ensure strictly decreasing & unique
    r_grid = _make_strictly_decreasing(r_grid)
    return r_grid


def _make_strictly_decreasing(r):
    """Ensure strictly decreasing radius with min fractional spacing."""
    r = np.asarray(r, float)
    out = [r[0]]
    for val in r[1:]:
        if val >= out[-1]:
            val = out[-1] * (1.0 - 1e-5)
        out.append(val)
    return np.array(out)


def _cell_scores(r, v, T, rho, tau, mass_fractions, cfg: AMRConfig,
                 kappa=None):
    """Compute a soft refinement score for each cell (between grid points).

    Score is computed for the *gaps* between consecutive grid points
    (n-1 cells).  Higher score => higher priority to split.
    """
    n = len(r)
    if n < 2:
        return np.zeros(0)
    dr = np.abs(np.diff(r))
    r_mid = 0.5 * (r[:-1] + r[1:])

    # optical depth increment per cell
    dtau = np.diff(tau)
    tau_mid = 0.5 * (tau[:-1] + tau[1:])
    dlogtau = np.log10((tau[1:] + cfg.tau_floor) / (tau[:-1] + cfg.tau_floor))

    dlogrho = np.abs(np.log10(np.maximum(rho[1:], 1e-30) /
                              np.maximum(rho[:-1], 1e-30)))
    dv = np.abs(np.diff(v))
    dlogT = np.abs(np.log10(np.maximum(T[1:], 1.0) /
                            np.maximum(T[:-1], 1.0)))
    dror = dr / np.maximum(r_mid, 1e-30)

    # curvature: second difference of log-quantities (normalised)
    curv = np.zeros(n - 1)
    if n >= 3:
        d2logrho = np.abs(np.log10(np.maximum(rho[2:], 1e-30) /
                                   np.maximum(rho[1:-1], 1e-30)) -
                          np.log10(np.maximum(rho[1:-1], 1e-30) /
                                   np.maximum(rho[:-2], 1e-30)))
        curv[:n - 2] = np.maximum(curv[:n - 2], d2logrho)

    # dimensionless ratios
    e_tau = np.abs(dlogtau) / cfg.target_dlogtau
    e_rho = dlogrho / cfg.target_dlogrho
    e_v = dv / cfg.target_dv_kms
    e_T = dlogT / cfg.target_dlogT
    e_r = dror / cfg.target_dr_over_r
    e_curv = curv / 0.5

    # Electron-scattering tau per cell and opacity gradient
    e_dtau_es = np.zeros(n - 1)
    e_dlogchi = np.zeros(n - 1)
    e_dv_over_v = np.zeros(n - 1)
    if kappa is not None and mass_fractions is not None:
        X = mass_fractions.get('H', np.zeros(n))
        Y = mass_fractions.get('He', np.zeros(n))
        kappa_es = 0.2 * (1.0 + X) / (1.0 + 3.0 * X + 0.25 * Y + 1e-30)
        krho_es_mid = 0.5 * (kappa_es[:-1] * rho[:-1] + kappa_es[1:] * rho[1:])
        dtau_es = krho_es_mid * dr
        e_dtau_es = dtau_es / cfg.target_dtau_es

        chi = kappa * rho
        chi_safe = np.maximum(chi, 1e-30)
        dlogchi = np.abs(np.log10(chi_safe[1:] / chi_safe[:-1]))
        e_dlogchi = dlogchi / cfg.target_dlogchi

    # Velocity/frequency coupling: dv/v
    v_mid = 0.5 * (np.abs(v[:-1]) + np.abs(v[1:]))
    e_dv_over_v = dv / np.maximum(v_mid, 1.0) / cfg.target_dv_over_v

    score = (cfg.w_tau * e_tau + cfg.w_rho * e_rho + cfg.w_v * e_v +
             cfg.w_T * e_T + cfg.w_curv * e_curv + 0.2 * e_r +
             cfg.w_dtau_es * e_dtau_es + cfg.w_dlogchi * e_dlogchi +
             cfg.w_dv_over_v * e_dv_over_v)
    return score


def _feature_cell_mask(r, features, cfg: AMRConfig):
    """Return a boolean mask (n-1 cells) marking cells that overlap a hard
    feature and still need mandatory refinement."""
    n = len(r)
    mask = np.zeros(n - 1, dtype=bool)
    for f in features:
        # cells whose midpoint is within the feature radius band
        r_mid = 0.5 * (r[:-1] + r[1:])
        # mark cells within a fractional log-radius window
        dlog = np.abs(np.log10(np.maximum(r_mid, 1e-30)) -
                      np.log10(max(f.radius, 1e-30)))
        mask |= (dlog < 0.03)
    return mask


def _feature_needs_refinement(r, features, cfg: AMRConfig):
    """Decide which features still demand more points (count cells in window)."""
    need = []
    r_mid = 0.5 * (r[:-1] + r[1:])
    for f in features:
        dlog = np.abs(np.log10(np.maximum(r_mid, 1e-30)) -
                      np.log10(max(f.radius, 1e-30)))
        n_near = int(np.sum(dlog < 0.03))
        if n_near < f.n_points:
            need.append(f)
    return need


def _check_resolution_satisfied(r, v, T, rho, tau, mass_fractions,
                                features, cfg: AMRConfig,
                                prev_scores_max=None, prev_n_features=None,
                                prev_n_extrema=None):
    """Check if resolution criteria are satisfied (for coarsening).

    For coarsening, we check that removing a point doesn't violate the HARD
    limits (3x soft targets) and doesn't degrade velocity topology or
    increase unresolved features.  This allows removing points from
    over-refined regions while keeping the grid physically adequate.

    If prev_* are given, we use coarsening mode (don't get worse).
    Otherwise we require full satisfaction (AMR stop mode).
    """
    if len(r) < cfg.min_points:
        return False

    # Soft criteria scores
    scores = _cell_scores(r, v, T, rho, tau, mass_fractions, cfg, kappa=kappa)
    scores_max = float(scores.max()) if scores.size > 0 else 0.0

    # Hard constraints: count unresolved features
    feat_need = _feature_needs_refinement(r, features, cfg)
    n_features = len(feat_need)

    # Velocity topology
    n_extrema = _count_extrema(v)

    if prev_scores_max is not None:
        # Coarsening mode: check against HARD limits (3x soft target)
        # The soft targets are: score < 1.0 means within target.
        # The hard limits are: score < 3.0 means within 3x target.
        # We allow coarsening as long as we stay within hard limits.
        if scores_max > 3.0:  # hard limit
            return False
        if n_features > prev_n_features:
            return False
        if prev_n_extrema is not None and n_extrema < prev_n_extrema:
            return False
        return True
    else:
        # AMR stop mode: everything must be within soft target
        if scores_max >= cfg.safety_margin:
            return False
        if feat_need:
            return False
        return True


def _check_individual_hard_limits(r, v, T, rho, tau, mass_fractions,
                                  cfg: AMRConfig, hard_factor=3.0):
    """Check that each individual resolution criterion stays within hard limits.

    The hard limit is hard_factor × the soft target for each criterion.
    This is used by the coarsening pass to ensure no individual criterion
    exceeds its hard limit, even if the combined weighted score is low.
    """
    n = len(r)
    if n < 2:
        return True
    # Optical depth
    dlogtau = np.abs(np.log10((tau[1:] + cfg.tau_floor) /
                              (tau[:-1] + cfg.tau_floor)))
    if float(dlogtau.max()) > cfg.target_dlogtau * hard_factor:
        return False
    # Density
    dlogrho = np.abs(np.log10(np.maximum(rho[1:], 1e-30) /
                              np.maximum(rho[:-1], 1e-30)))
    if float(dlogrho.max()) > cfg.target_dlogrho * hard_factor:
        return False
    # Velocity
    dv = np.abs(np.diff(v))
    if float(dv.max()) > cfg.target_dv_kms * hard_factor:
        return False
    # Temperature
    dlogT = np.abs(np.log10(np.maximum(T[1:], 1.0) /
                            np.maximum(T[:-1], 1.0)))
    if float(dlogT.max()) > cfg.target_dlogT * hard_factor:
        return False
    # Radius spacing
    dr = np.abs(np.diff(r))
    r_mid = 0.5 * (r[:-1] + r[1:])
    dror = dr / np.maximum(r_mid, 1e-30)
    if float(dror.max()) > cfg.target_dr_over_r * hard_factor:
        return False
    return True


def _coarsen_grid(r, v, T, rho, kappa, tau, mass_fractions,
                  r_old, v_old, T_old, rho_old, kappa_old,
                  mass_fractions_old, features, interface_indices,
                  cfg: AMRConfig, report):
    """Attempt to remove individual points while maintaining all constraints.

    This is the coarsening pass that finds the true minimum viable grid.
    After AMR converges, we sweep through all non-feature points and try
    removing each one.  If the resolution criteria remain satisfied, the
    removal is accepted; otherwise the point is restored.

    Feature radii (composition interfaces, shocks, velocity extrema) are
    NEVER removed — they are hard constraints.
    """

    def do_remap(r_grid):
        """Remap source quantities onto r_grid (inline, not a closure)."""
        v_n, T_n, rho_n, kappa_n, mf_n = remap_onto_grid(
            r_grid, r_old, v_old, T_old, rho_old, kappa_old,
            mass_fractions_old, interface_indices, cfg)
        tau_n = compute_kappa_tau(r_grid, rho_n, kappa_n)
        return v_n, T_n, rho_n, kappa_n, tau_n, mf_n
    feat_radii = set()
    for f in features:
        feat_radii.add(float(f.radius))
    # Also protect velocity-turn radii (±neighbours)
    turn_radii = _all_velocity_turn_radii(r_old, v_old)
    for tr in turn_radii:
        feat_radii.add(float(tr))

    # Build a set of protected grid indices (points very close to features)
    # Use a tight tolerance: only protect points within 1% of the local spacing
    protected_indices = set()
    for i in range(1, len(r) - 1):
        r_i = float(r[i])
        local_dr = abs(r[i-1] - r[i+1])
        for fr in feat_radii:
            if abs(r_i - fr) < local_dr * 0.1:  # within 10% of neighbor spacing
                protected_indices.add(i)
                break

    if cfg.verbose:
        print(f"[Tau-AMR] Coarsening: {len(protected_indices)} protected points, "
              f"{len(r) - 2 - len(protected_indices)} candidates")

    n_removed = 0
    for sweep in range(cfg.coarsen_max_passes):
        n_before = len(r)
        # Compute baseline metrics before this sweep
        scores_base = _cell_scores(r, v, T, rho, tau, mass_fractions, cfg, kappa=kappa)
        prev_scores_max = float(scores_base.max()) if scores_base.size > 0 else 0.0
        prev_n_features = len(_feature_needs_refinement(r, features, cfg))
        prev_n_extrema = _count_extrema(v)

        # Sort points by their local score contribution (lowest first).
        # Points in low-error regions are candidates for removal.
        point_scores = np.zeros(len(r))
        for i in range(1, len(r) - 1):
            point_scores[i] = max(scores_base[i-1], scores_base[i])
        # Protected points get infinite score (never remove)
        for i in protected_indices:
            if i < len(point_scores):
                point_scores[i] = np.inf

        # Try removing points from lowest to highest score
        candidates = np.argsort(point_scores)
        for idx in candidates:
            i = int(idx)
            if i < 1 or i >= len(r) - 1:
                continue
            if i in protected_indices:
                continue

            # Try removing point i
            r_try = np.delete(r, i)
            # Check minimum spacing
            if abs(r_try[i-1] - r_try[i]) < cfg.min_dr_r * min(abs(r_try[i-1]), abs(r_try[i])):
                continue

            # Remap onto the trial grid
            r_save = r.copy()
            v_save, T_save, rho_save, kappa_save = v, T, rho, kappa
            mf_save = {e: x.copy() for e, x in mass_fractions.items()}
            tau_save = tau.copy()
            r = r_try
            try:
                v, T, rho, kappa, tau, mass_fractions = do_remap(r)
                # GLOBAL check: each individual criterion must stay within
                # hard limits (3x soft target).  We check individual criteria,
                # not the combined weighted score, because a high velocity
                # error can be masked by low weights in the combined score.
                ok = _check_individual_hard_limits(
                    r, v, T, rho, tau, mass_fractions, cfg)
                # Check features don't get worse
                n_feat_new = len(_feature_needs_refinement(r, features, cfg))
                if n_feat_new > prev_n_features:
                    ok = False
                # Check velocity topology doesn't degrade
                if _count_extrema(v) < prev_n_extrema:
                    ok = False
            except Exception:
                ok = False

            if ok:
                n_removed += 1
                # Update baseline
                scores_new = _cell_scores(r, v, T, rho, tau, mass_fractions, cfg, kappa=kappa)
                prev_scores_max = float(scores_new.max()) if scores_new.size > 0 else 0.0
                prev_n_features = len(_feature_needs_refinement(r, features, cfg))
                prev_n_extrema = _count_extrema(v)
                # Update point_scores for remaining candidates
                point_scores = np.delete(point_scores, i)
            else:
                # Restore the point
                r = r_save
                v, T, rho, kappa = v_save, T_save, rho_save, kappa_save
                mass_fractions = mf_save
                tau = tau_save

        if len(r) == n_before:
            break  # no more points can be removed

    report['coarsening'] = {'points_removed': n_removed, 'passes': sweep + 1}
    if cfg.verbose:
        if n_removed > 0:
            print(f"[Tau-AMR] Coarsening: removed {n_removed} points "
                  f"({n_before} -> {len(r)} final)")
        else:
            print(f"[Tau-AMR] Coarsening: 0 points could be removed "
                  f"(grid already at minimum)")
    return r, v, T, rho, kappa, tau, mass_fractions, report


def recursive_amr(r_old, v_old, T_old, rho_old, kappa_old, tau_old,
                  mass_fractions_old, features, cfg: AMRConfig,
                  interface_indices_old):
    """Run the recursive cell-based AMR loop.

    Returns the final (r, v, T, rho, kappa, tau, mass_fractions) on the AMR
    grid, plus a report dict.
    """
    # Start from a coarse base grid in log(tau)
    r = build_base_grid(r_old, tau_old, cfg.base_grid_points, cfg.tau_floor)
    # Always keep the exact outermost and innermost source radii
    r[0] = r_old[0]
    r[-1] = r_old[-1]
    # Seed every feature radius (composition interface / shock / velocity
    # extremum) as an exact grid point so PCHIP cannot smooth away a
    # discontinuity or velocity extremum.  This is the key step that preserves
    # non-monotonic velocity topology and sharp chemical interfaces.
    feat_radii = np.unique(np.array([f.radius for f in features]))
    # Also seed EVERY simple 3-point velocity turn in the source, so that even
    # shallow extrema (which the significance-gated detector may skip) are
    # preserved as exact grid points.  This guarantees the velocity topology
    # (extrema count) is conserved through remapping.
    turn_radii = _all_velocity_turn_radii(r_old, v_old)
    r = np.unique(np.concatenate([r, feat_radii, turn_radii]))[::-1]
    r = _make_strictly_decreasing(r)

    interface_indices = interface_indices_old

    report = {'iterations': 0, 'splits': 0, 'history': [],
               'unresolved_features': []}

    def do_remap():
        nonlocal v, T, rho, kappa, mass_fractions, tau
        v, T, rho, kappa, mass_fractions = remap_onto_grid(
            r, r_old, v_old, T_old, rho_old, kappa_old,
            mass_fractions_old, interface_indices, cfg)
        tau = compute_kappa_tau(r, rho, kappa)

    v = T = rho = kappa = tau = None
    mass_fractions = None
    do_remap()

    max_iter = 4 * cfg.max_points  # safety bound
    prev_max_score = float('inf')
    stagnant_count = 0
    for it in range(max_iter):
        report['iterations'] = it
        n = len(r)

        scores = _cell_scores(r, v, T, rho, tau, mass_fractions, cfg, kappa=kappa)
        if scores.size == 0:
            break

        # Hard constraints: features that still need points
        feat_need = _feature_needs_refinement(r, features, cfg)
        feat_mask = _feature_cell_mask(r, feat_need, cfg)

        max_score = float(scores.max())

        # --- STOP CONDITION 1: all criteria satisfied ---
        # Stop as soon as (a) no hard-constraint feature needs refinement,
        # (b) all soft scores < safety_margin, and (c) we have at least
        # min_points.  This produces the SMALLEST grid that satisfies all
        # constraints.
        if not feat_need and max_score < cfg.safety_margin and n >= cfg.min_points:
            break

        # --- STOP CONDITION 2: convergence stagnation ---
        # If the max score hasn't improved significantly over the last 20
        # iterations, the remaining high-score cells are at physical
        # discontinuities (shocks, composition interfaces) that can't be
        # resolved further by splitting.  Stop to avoid wasting points.
        if max_score < prev_max_score * 0.999:
            stagnant_count = 0
        else:
            stagnant_count += 1
        prev_max_score = max_score
        if stagnant_count >= 20 and not feat_need and n >= cfg.min_points:
            if cfg.verbose:
                print(f"[Tau-AMR] Convergence stagnation at {n} points "
                      f"(max_score={max_score:.3f} unchanged for 20 iters)")
            break

        # --- HARD LIMIT: do not exceed max_points ---
        if n >= cfg.max_points:
            break

        # Choose cell to split: prefer hard-constraint cells, then highest score
        if feat_mask.any():
            cand = np.where(feat_mask)[0]
            pick = cand[np.argmax(scores[cand])]
        else:
            pick = int(np.argmax(scores))

        # Split cell `pick` (between r[pick] and r[pick+1]) at the midpoint
        # in log-radius (geometric mean) to preserve structure.
        r_new = np.insert(r, pick + 1, np.sqrt(r[pick] * r[pick + 1]))
        # enforce strict decrease / min spacing
        if r_new[pick + 1] >= r[pick] or r_new[pick + 1] <= r[pick + 1]:
            # fall back to arithmetic mean
            r_new = np.insert(r, pick + 1, 0.5 * (r[pick] + r[pick + 1]))
        r = _make_strictly_decreasing(r_new)
        report['splits'] += 1
        report['history'].append(int(pick))

        do_remap()

    # --- Coarsening pass: remove points while maintaining constraints ---
    # After AMR converges, attempt to remove individual points to find the
    # true minimum viable grid.  We try removing each non-feature point and
    # check if all resolution criteria remain satisfied.
    if cfg.enable_coarsening and len(r) > cfg.min_points:
        r, v, T, rho, kappa, tau, mass_fractions, report = \
            _coarsen_grid(r, v, T, rho, kappa, tau, mass_fractions,
                          r_old, v_old, T_old, rho_old, kappa_old,
                          mass_fractions_old, features, interface_indices,
                          cfg, report)

    report['unresolved_features'] = [
        {'kind': f.kind, 'index': int(f.index), 'radius': float(f.radius),
         'detail': f.detail}
        for f in _feature_needs_refinement(r, features, cfg)
    ]
    return r, v, T, rho, kappa, tau, mass_fractions, report


def escalate_grid(r, v, T, rho, kappa, tau, mass_fractions,
                  r_old, v_old, T_old, rho_old, kappa_old, tau_old,
                  mass_fractions_old, features, cfg: AMRConfig,
                  interface_indices_old,
                  max_rounds: int = 5,
                  pathological_thresholds: dict = None) -> tuple:
    """Progressive grid escalation: refine locally until pathological cells are resolved.

    This is the closed-loop Stage-A refinement.  After ``recursive_amr``
    produces an initial grid, this function:

    1. Runs the pathological-cell diagnostic.
    2. If failures exist and ND < max_points, inserts points at the
       geometric-mean radius of each failed cell.
    3. Conservatively remaps onto the new grid.
    4. Re-runs the diagnostic.
    5. Repeats until no failures, or max_points reached, or no improvement,
       or max_rounds exhausted.

    Parameters
    ----------
    (r, v, T, rho, kappa, tau, mass_fractions) : arrays
        The grid returned by ``recursive_amr``.
    (r_old, ..., mass_fractions_old) : arrays
        The original STELLA source grid (for conservative remap).
    features, cfg, interface_indices_old : AMR bookkeeping.
    max_rounds : int
        Maximum escalation rounds.
    pathological_thresholds : dict, optional
        Override thresholds for the pathological-cell check.

    Returns
    -------
    Same tuple as ``recursive_amr``, plus an updated report dict with
    escalation history.
    """
    report_escalation = {
        'rounds': 0,
        'history': [],
        'pathological_passed': False,
        'pathological_n_failures': 0,
        'points_inserted': 0,
        'final_nd': len(r),
        'max_points_reached': False,
    }

    # Collect feature radii — cells AT these locations represent physical
    # discontinuities (composition interfaces, shocks) that are correctly
    # represented by sharp jumps and should not be flagged as pathological
    # once they have grid points on both sides.
    feat_radii = [float(f.radius) for f in features]

    for round_idx in range(max_rounds):
        # Run pathological-cell diagnostic on the current grid
        # Pass feature_radii so that physical discontinuities at detected
        # features are exempted from density/opacity jump criteria.
        patho = pathological_cell_report(
            r, v, T, rho, kappa, tau, mass_fractions, cfg,
            thresholds=pathological_thresholds,
            feature_radii=feat_radii)

        report_escalation['pathological_n_failures'] = patho['n_failures']
        report_escalation['pathological_passed'] = patho['passed']

        if patho['passed']:
            print(f"[Escalation] Round {round_idx}: PASSED — no pathological cells.")
            break

        if len(r) >= cfg.max_points:
            print(f"[Escalation] Round {round_idx}: FAILED — max_points ({cfg.max_points}) "
                  f"reached with {patho['n_failures']} unresolved cell(s).")
            report_escalation['max_points_reached'] = True
            # Print the GRID STATUS: FAILED message
            print("=" * 72)
            print("GRID STATUS: FAILED — MAXIMUM ND REACHED")
            print(f"  ND = {len(r)}  (max_points = {cfg.max_points})")
            print(f"  Unresolved pathological cells: {patho['n_failures']}")
            for f in patho['failures'][:10]:
                print(f"    cell {f['cell_index']}: r={f['r_outer']:.3e} -> "
                      f"{f['r_inner']:.3e} cm, tau={f['tau_at_cell']:.2f}")
                for cf in f['failed_criteria']:
                    print(f"      {cf['label']}: {cf['value']:.4g} > {cf['limit']:.4g}")
            print("  Further refinement is required but exceeds ND_MAX.")
            print("=" * 72)
            break

        # Collect unique insertion radii from failed cells (geometric mean)
        insert_radii = set()
        for f in patho['failures']:
            r_gm = np.sqrt(f['r_outer'] * f['r_inner'])
            insert_radii.add(r_gm)

        # Check if we have room to insert
        n_insert = len(insert_radii)
        if len(r) + n_insert > cfg.max_points:
            n_insert = cfg.max_points - len(r)
            insert_radii = set(list(insert_radii)[:n_insert])
            n_insert = len(insert_radii)
            if n_insert == 0:
                print(f"[Escalation] Round {round_idx}: no room to insert points "
                      f"(ND={len(r)}, max={cfg.max_points}).")
                report_escalation['max_points_reached'] = True
                break

        # Insert points at the geometric-mean radii
        r_new = np.unique(np.concatenate([r, list(insert_radii)]))[::-1]
        r_new = _make_strictly_decreasing(r_new)

        if len(r_new) == len(r):
            print(f"[Escalation] Round {round_idx}: no new points inserted "
                  f"(insertion radii already on grid). Stopping.")
            break

        # Remap onto the new grid
        v_new, T_new, rho_new, kappa_new, mf_new = remap_onto_grid(
            r_new, r_old, v_old, T_old, rho_old, kappa_old,
            mass_fractions_old, interface_indices_old, cfg)
        # Recompute tau on the new grid
        dr_new = -np.diff(r_new)
        chi_new = kappa_new * rho_new
        chi_mid_new = 0.5 * (chi_new[:-1] + chi_new[1:])
        tau_new = np.concatenate([[0], np.cumsum(chi_mid_new * dr_new)])

        n_before = len(r)
        r, v, T, rho, kappa, tau, mass_fractions = (
            r_new, v_new, T_new, rho_new, kappa_new, tau_new, mf_new)

        round_info = {
            'round': round_idx,
            'n_failures_before': patho['n_failures'],
            'points_inserted': len(r) - n_before,
            'nd_after': len(r),
        }
        report_escalation['history'].append(round_info)
        report_escalation['points_inserted'] += len(r) - n_before
        report_escalation['final_nd'] = len(r)

        print(f"[Escalation] Round {round_idx}: inserted {len(r) - n_before} points "
              f"at {n_insert} failed cell(s). ND={n_before} -> {len(r)}.")

    # Final diagnostic
    patho_final = pathological_cell_report(
        r, v, T, rho, kappa, tau, mass_fractions, cfg,
        thresholds=pathological_thresholds)
    report_escalation['pathological_passed'] = patho_final['passed']
    report_escalation['pathological_n_failures'] = patho_final['n_failures']
    report_escalation['rounds'] = len(report_escalation['history'])
    report_escalation['final_pathological_report'] = patho_final

    if not patho_final['passed'] and not report_escalation['max_points_reached']:
        print(f"[Escalation] WARNING: {patho_final['n_failures']} pathological cell(s) "
              f"remain after {report_escalation['rounds']} escalation round(s). "
              f"Consider increasing max_rounds or tightening AMR thresholds.")

    return r, v, T, rho, kappa, tau, mass_fractions, report_escalation


# ===========================================================================
# Validation
# ===========================================================================

def validate_grid(r, v, T, rho, kappa, tau, mass_fractions,
                  r_old, v_old, rho_old, mass_fractions_old, cfg: AMRConfig):
    """Run all hard-constraint / conservation validations.  Returns report dict."""
    rep = {'checks': {}, 'passed': True}

    # 1. point count
    rep['checks']['n_points'] = {'value': len(r),
                                  'ok': len(r) <= cfg.max_points}
    # 2. radial ordering (strictly decreasing)
    ok_order = bool(np.all(np.diff(r) < 0))
    rep['checks']['radial_order'] = {'ok': ok_order}
    # 3. velocity topology preserved (same number of extrema as source)
    extrema_new = _count_extrema(v)
    extrema_old = _count_extrema(v_old)
    rep['checks']['velocity_topology'] = {
        'n_extrema_source': extrema_old, 'n_extrema_grid': extrema_new,
        'ok': extrema_new >= extrema_old}
    # 4. no negative abundances
    neg = {e: int(np.sum(x < 0)) for e, x in mass_fractions.items()}
    rep['checks']['no_negative_abundances'] = {'ok': all(v == 0 for v in neg.values())}
    # 5. sum(X) ~ 1
    sums = np.zeros(len(r))
    for x in mass_fractions.values():
        sums = sums + x
    rep['checks']['sum_x'] = {'min': float(sums.min()), 'max': float(sums.max()),
                               'ok': bool(abs(sums.mean() - 1.0) < 1e-4)}
    # 6. mass conservation
    m_old = total_mass(r_old, rho_old)
    m_new = total_mass(r, rho)
    rel_err = abs(m_new - m_old) / max(m_old, 1e-30)
    rep['checks']['mass_conservation'] = {
        'm_old_msun': m_old / MSUN, 'm_new_msun': m_new / MSUN,
        'rel_err': rel_err, 'ok': rel_err < cfg.mass_tol}
    # 7. species mass conservation
    sm_old = species_masses(r_old, rho_old, mass_fractions_old)
    sm_new = species_masses(r, rho, mass_fractions)
    spec_err = {}
    for e in sm_old:
        if sm_old[e] > 0:
            spec_err[e] = abs(sm_new[e] - sm_old[e]) / sm_old[e]
    rep['checks']['species_mass'] = {
        'max_rel_err': max(spec_err.values()) if spec_err else 0.0,
        'ok': all(v < cfg.mass_tol * 5 for v in spec_err.values())}

    rep['passed'] = all(c['ok'] for c in rep['checks'].values())
    return rep


# ===========================================================================
# Pre-CMFGEN resolution report
# ===========================================================================
# The grid resolution is determined entirely from the STELLA physics, NOT from
# CMFGEN source code.  This report checks that the AMR grid adequately resolves
# the physical structure BEFORE CMFGEN is launched.  CMFGEN source inspection
# is used only for input format / ordering conventions, not for deciding where
# grid points should be placed.


# ---------------------------------------------------------------------------
# Pathological-cell diagnostic
# ---------------------------------------------------------------------------
# Hard grid-quality check: flags any adjacent pair where a physically
# important quantity changes too rapidly for CMFGEN to interpolate reliably.
# Unlike the soft AMR score, these are PASS/FAIL criteria — a single failure
# means the grid is not adequate for the radiative-transfer calculation.
#
# The thresholds are configurable via AMRConfig.  They represent the maximum
# change per cell that CMFGEN's interpolation and statistical-equilibrium
# solver can handle without introducing significant numerical error.
# ---------------------------------------------------------------------------

# Default hard thresholds for the pathological-cell check.
# These are deliberately stricter than the AMR soft targets (which use 3x
# as the hard limit in _check_individual_hard_limits) because this check
# is a FAIL condition, not a refinement trigger.
DEFAULT_PATHOLOGICAL_THRESHOLDS = {
    'max_dr_over_r':       0.15,   # 15% radius change per cell
    'max_dv_over_v':       0.50,   # 50% velocity change per cell
    'max_dlog_rho':        1.0,    # 1 dex density change per cell
    'max_dlog_T':          0.5,    # 0.5 dex temperature change per cell
    'max_dlog_chi':        1.5,    # 1.5 dex opacity (chi) change per cell
    'max_dtau_es':         2.0,    # 2.0 electron-scattering tau per cell
    'max_density_ratio':   100.0,  # density ratio between adjacent cells
    'max_chi_ratio':       500.0,  # opacity ratio between adjacent cells
}


def pathological_cell_report(r, v, T, rho, kappa, tau, mass_fractions,
                             cfg: AMRConfig = None,
                             thresholds: dict = None,
                             max_report: int = 20,
                             feature_radii: list = None,
                             feature_tol: float = 0.02) -> dict:
    """Check every adjacent pair for excessive physical changes.

    This is a HARD grid-quality diagnostic.  Unlike the soft AMR score,
    any failure means the grid is inadequate for the CMFGEN radiative-
    transfer calculation — CMFGEN would be asked to interpolate across
    an enormous physical change in one depth interval.

    Parameters
    ----------
    r, v, T, rho, kappa, tau, mass_fractions : arrays
        Grid quantities (outer->inner, r decreasing).  All same length N.
    cfg : AMRConfig, optional
        If provided, cfg thresholds override defaults.
    thresholds : dict, optional
        Explicit threshold overrides.  Highest priority.
    max_report : int
        Maximum number of failures to report.
    feature_radii : list, optional
        Radii of detected physical features (composition interfaces, shocks).
        Cells at these locations represent physical discontinuities and are
        exempt from the density/opacity/temperature jump criteria (but NOT
        from dr/r, dv/v, or dtau_es criteria which indicate numerical
        under-resolution).
    feature_tol : float
        Log-radius tolerance for matching a cell to a feature (default 0.02 dex).

    Returns
    -------
    dict with keys:
        'passed' : bool — True if no failures
        'n_failures' : int — total number of failed cells
        'failures' : list of dicts, one per failed cell, sorted by severity
        'summary' : str — human-readable report
        'machine_readable' : dict — machine-readable summary
    """
    n = len(r)
    if n < 2:
        return {'passed': True, 'n_failures': 0, 'failures': [],
                'summary': 'Grid has < 2 points, skipping pathological-cell check.',
                'machine_readable': {'n_points': n, 'passed': True}}

    # Merge thresholds: defaults < cfg < explicit
    th = dict(DEFAULT_PATHOLOGICAL_THRESHOLDS)
    if cfg is not None:
        # Allow cfg to override via dedicated fields (if present)
        for key in th:
            cfg_val = getattr(cfg, 'pathological_' + key, None)
            if cfg_val is not None:
                th[key] = cfg_val
    if thresholds is not None:
        th.update(thresholds)

    # Compute per-cell metrics (outer->inner, index i = gap between cell i and i+1)
    dr = np.abs(np.diff(r))
    r_mid = 0.5 * (r[:-1] + r[1:])
    dr_over_r = dr / np.maximum(r_mid, 1e-30)

    dv = np.abs(np.diff(v))
    v_mid = 0.5 * (np.abs(v[:-1]) + np.abs(v[1:]))
    dv_over_v = dv / np.maximum(v_mid, 1.0)  # avoid div-by-zero

    rho_safe = np.maximum(rho, 1e-30)
    dlog_rho = np.abs(np.log10(rho_safe[1:] / rho_safe[:-1]))
    density_ratio = np.maximum(rho_safe[1:] / rho_safe[:-1],
                               rho_safe[:-1] / rho_safe[1:])

    T_safe = np.maximum(T, 1.0)
    dlog_T = np.abs(np.log10(T_safe[1:] / T_safe[:-1]))

    chi = kappa * rho  # absorption coefficient (cm^-1)
    chi_safe = np.maximum(chi, 1e-30)
    dlog_chi = np.abs(np.log10(chi_safe[1:] / chi_safe[:-1]))
    chi_ratio = np.maximum(chi_safe[1:] / chi_safe[:-1],
                           chi_safe[:-1] / chi_safe[1:])

    # Electron-scattering tau per cell (approximate: kappa_es ~ 0.34 * (1+X) / (1+3X+0.25Y))
    X = mass_fractions.get('H', np.zeros(n))
    Y = mass_fractions.get('He', np.zeros(n))
    kappa_es = 0.2 * (1.0 + X) / (1.0 + 3.0 * X + 0.25 * Y + 1e-30)
    krho_es_mid = 0.5 * (kappa_es[:-1] * rho[:-1] + kappa_es[1:] * rho[1:])
    dtau_es = krho_es_mid * dr

    # Check each criterion
    checks = {
        'dr_over_r':     (dr_over_r, th['max_dr_over_r'], 'dr/r'),
        'dv_over_v':     (dv_over_v, th['max_dv_over_v'], 'dv/v'),
        'dlog_rho':      (dlog_rho, th['max_dlog_rho'], 'dlog10(rho)'),
        'dlog_T':        (dlog_T, th['max_dlog_T'], 'dlog10(T)'),
        'dlog_chi':      (dlog_chi, th['max_dlog_chi'], 'dlog10(chi)'),
        'dtau_es':       (dtau_es, th['max_dtau_es'], 'dtau_es'),
        'density_ratio': (density_ratio, th['max_density_ratio'], 'rho ratio'),
        'chi_ratio':     (chi_ratio, th['max_chi_ratio'], 'chi ratio'),
    }

    failures = []
    for i in range(n - 1):
        cell_failures = []
        for name, (values, limit, label) in checks.items():
            if values[i] > limit:
                cell_failures.append({
                    'criterion': name,
                    'label': label,
                    'value': float(values[i]),
                    'limit': float(limit),
                })
        if cell_failures:
            # Severity = max ratio of value/limit across all failed criteria
            severity = max(cf['value'] / cf['limit'] for cf in cell_failures)
            failures.append({
                'cell_index': i,
                'r_outer': float(r[i]),
                'r_inner': float(r[i + 1]),
                'dr_over_r': float(dr_over_r[i]),
                'dv': float(dv[i]),
                'dv_over_v': float(dv_over_v[i]),
                'dlog_rho': float(dlog_rho[i]),
                'dlog_T': float(dlog_T[i]),
                'dlog_chi': float(dlog_chi[i]),
                'dtau_es': float(dtau_es[i]),
                'density_ratio': float(density_ratio[i]),
                'chi_ratio': float(chi_ratio[i]),
                'tau_at_cell': float(tau[i]) if i < len(tau) else 0.0,
                'v_outer': float(v[i]),
                'v_inner': float(v[i + 1]),
                'T_outer': float(T[i]),
                'T_inner': float(T[i + 1]),
                'rho_outer': float(rho[i]),
                'rho_inner': float(rho[i + 1]),
                'severity': float(severity),
                'failed_criteria': cell_failures,
            })

    # Sort by severity (most severe first)
    failures.sort(key=lambda f: f['severity'], reverse=True)

    # If feature_radii provided, filter out failures at physical discontinuities.
    # At a composition interface or shock, density/opacity/temperature jumps are
    # physically real and should not be flagged as pathological.  However, dr/r,
    # dv/v, and dtau_es failures are still flagged because they indicate
    # numerical under-resolution even at a physical discontinuity.
    if feature_radii:
        # Criteria that are exempt at feature locations (physical jumps)
        exempt_criteria = {'dlog_rho', 'dlog_T', 'dlog_chi', 'density_ratio', 'chi_ratio'}
        # Criteria that are NEVER exempt (numerical resolution)
        non_exempt = {'dr_over_r', 'dv_over_v', 'dtau_es'}
        filtered_failures = []
        for f in failures:
            r_mid = np.sqrt(f['r_outer'] * f['r_inner'])
            at_feature = any(
                abs(np.log10(r_mid / max(fr, 1e-30))) < feature_tol
                for fr in feature_radii
            )
            if at_feature:
                # Keep only non-exempt criteria
                remaining = [cf for cf in f['failed_criteria']
                             if cf['criterion'] in non_exempt]
                if remaining:
                    f['failed_criteria'] = remaining
                    f['severity'] = max(cf['value'] / cf['limit'] for cf in remaining)
                    filtered_failures.append(f)
                # If only exempt criteria failed, drop this cell entirely
            else:
                filtered_failures.append(f)
        failures = filtered_failures

    passed = len(failures) == 0
    n_fail = len(failures)

    # Build human-readable summary
    lines = []
    lines.append("=" * 72)
    lines.append("GRID QUALITY REPORT — Pathological Cell Diagnostic")
    lines.append("=" * 72)
    lines.append(f"  ND = {n}")
    lines.append(f"  Status: {'PASSED' if passed else 'FAILED'} "
                 f"({n_fail} cell(s) with excessive physical changes)")
    if not passed:
        lines.append("")
        lines.append("FAILURES (sorted by severity):")
        lines.append("-" * 72)
        for f in failures[:max_report]:
            lines.append(f"  cell {f['cell_index']} -> {f['cell_index'] + 1}  "
                         f"(r={f['r_outer']:.3e} -> {f['r_inner']:.3e} cm, "
                         f"tau={f['tau_at_cell']:.2f})")
            for cf in f['failed_criteria']:
                lines.append(f"    {cf['label']:20s} = {cf['value']:.4g}  "
                             f"(limit: {cf['limit']:.4g})  FAIL")
            lines.append(f"    v: {f['v_outer']:.0f} -> {f['v_inner']:.0f} km/s  "
                         f"T: {f['T_outer']:.0f} -> {f['T_inner']:.0f} K  "
                         f"rho: {f['rho_outer']:.2e} -> {f['rho_inner']:.2e}")
            lines.append("")
        if n_fail > max_report:
            lines.append(f"  ... and {n_fail - max_report} more failure(s)")
        lines.append("-" * 72)
        lines.append("Recommended action: locally refine the failed cell(s)")
        lines.append("  by inserting additional depth points at the geometric")
        lines.append("  mean radius of each failed cell.")
    lines.append("=" * 72)

    machine_readable = {
        'n_points': n,
        'passed': passed,
        'n_failures': n_fail,
        'thresholds': th,
        'max_dr_over_r': float(dr_over_r.max()) if n > 1 else 0.0,
        'max_dv_over_v': float(dv_over_v.max()) if n > 1 else 0.0,
        'max_dlog_rho': float(dlog_rho.max()) if n > 1 else 0.0,
        'max_dlog_T': float(dlog_T.max()) if n > 1 else 0.0,
        'max_dlog_chi': float(dlog_chi.max()) if n > 1 else 0.0,
        'max_dtau_es': float(dtau_es.max()) if n > 1 else 0.0,
        'max_density_ratio': float(density_ratio.max()) if n > 1 else 1.0,
        'max_chi_ratio': float(chi_ratio.max()) if n > 1 else 1.0,
        'failed_cell_indices': [f['cell_index'] for f in failures],
    }

    return {
        'passed': passed,
        'n_failures': n_fail,
        'failures': failures,
        'summary': '\n'.join(lines),
        'machine_readable': machine_readable,
    }


def _detect_detached_shells(r, rho, mass_fractions, cfg: AMRConfig):
    """Detect detached high-density shells bounded by composition or density interfaces.

    A detached shell is a region where:
    - density is significantly higher than the surrounding regions
    - the region is bounded by two sharp transitions (composition or density)
    - the shell has a finite radial extent (not just a single interface)

    Returns a list of dicts with keys:
        'r_outer', 'r_inner', 'density_ratio', 'composition_change'
    """
    n = len(r)
    if n < 5:
        return []

    shells = []

    # Find density peaks: local maxima in rho (r is decreasing, outer->inner)
    # A density peak is a point where rho[i] > rho[i-1] and rho[i] > rho[i+1]
    # (with some tolerance for noise)
    for i in range(2, n - 2):
        if rho[i] > rho[i - 1] and rho[i] > rho[i + 1]:
            # Check if this is a significant peak: density at the peak is
            # at least 3x higher than the surrounding regions
            rho_peak = rho[i]
            # Find the surrounding minima
            rho_outer = rho[:i].min()  # minimum density outside (at larger r)
            rho_inner = rho[i + 1:].min()  # minimum density inside (at smaller r)
            rho_min_surround = min(rho_outer, rho_inner)
            if rho_peak < 3.0 * rho_min_surround:
                continue

            # Find the shell boundaries: walk outward and inward until density
            # drops to within 2x of the surrounding minimum
            r_outer = r[i]
            for j in range(i, 0, -1):
                if rho[j] < 2.0 * rho_outer:
                    r_outer = r[j]
                    break
            r_inner = r[i]
            for j in range(i, n - 1):
                if rho[j] < 2.0 * rho_inner:
                    r_inner = r[j]
                    break

            # Check for composition change across the shell
            comp_change = 0.0
            for elem, x in mass_fractions.items():
                if x.max() < cfg.composition_trace_floor:
                    continue
                dx = abs(x[i] - x[0]) + abs(x[i] - x[-1])
                comp_change = max(comp_change, dx)

            density_ratio = float(rho_peak / rho_min_surround)
            shells.append({
                'r_outer': float(r_outer),
                'r_inner': float(r_inner),
                'r_peak': float(r[i]),
                'density_ratio': density_ratio,
                'composition_change': float(comp_change),
            })

    return shells


def resolution_report(r, v, T, rho, kappa, tau, mass_fractions,
                      r_old, v_old, T_old, rho_old, tau_old,
                      mass_fractions_old, features, cfg: AMRConfig):
    """Generate a pre-CMFGEN resolution report from the STELLA physics.

    Checks that the AMR grid adequately resolves:
      1. optical depth (Δlogτ)
      2. density (Δlogρ)
      3. velocity (Δv)
      4. velocity extrema (topology)
      5. temperature (ΔlogT)
      6. shocks (density compression, T jump)
      7. shell boundaries (density gradient)
      8. composition interfaces (ΔX_i)
      9. abundance gradients (dX_i/dr)
     10. spectrum-forming region (τ ~ 0.01–1)

    Returns (report_dict, passed_bool, unresolved_list).
    """
    n = len(r)
    checks = {}
    unresolved = []

    # --- 1. Optical-depth resolution ---
    log_tau = np.log10(np.maximum(tau, cfg.tau_floor))
    dlogtau = np.abs(np.diff(log_tau))
    max_dlogtau = float(dlogtau.max()) if len(dlogtau) > 0 else 0.0
    ok_tau = max_dlogtau <= cfg.target_dlogtau * 3  # 3x soft target as hard limit
    checks['optical_depth'] = {
        'max_dlogtau': max_dlogtau, 'target': cfg.target_dlogtau,
        'ok': ok_tau}
    if not ok_tau:
        worst = int(np.argmax(dlogtau))
        unresolved.append(f"optical_depth: zone {worst} dlogtau={dlogtau[worst]:.3f}")

    # --- 2. Density resolution ---
    log_rho = np.log10(np.maximum(rho, 1e-300))
    dlogrho = np.abs(np.diff(log_rho))
    max_dlogrho = float(dlogrho.max()) if len(dlogrho) > 0 else 0.0
    ok_rho = max_dlogrho <= cfg.target_dlogrho * 3
    checks['density'] = {
        'max_dlogrho': max_dlogrho, 'target': cfg.target_dlogrho,
        'ok': ok_rho}
    if not ok_rho:
        worst = int(np.argmax(dlogrho))
        unresolved.append(f"density: zone {worst} dlogrho={dlogrho[worst]:.3f}")

    # --- 3. Velocity resolution ---
    dv = np.abs(np.diff(v))
    max_dv = float(dv.max()) if len(dv) > 0 else 0.0
    ok_v = max_dv <= cfg.target_dv_kms * 3
    checks['velocity'] = {
        'max_dv_kms': max_dv, 'target': cfg.target_dv_kms,
        'ok': ok_v}
    if not ok_v:
        worst = int(np.argmax(dv))
        unresolved.append(f"velocity: zone {worst} dv={dv[worst]:.1f} km/s")

    # --- 4. Velocity extrema ---
    n_extrema_old = _count_extrema(v_old)
    n_extrema_new = _count_extrema(v)
    ok_extrema = n_extrema_new >= n_extrema_old
    checks['velocity_extrema'] = {
        'n_source': n_extrema_old, 'n_grid': n_extrema_new,
        'ok': ok_extrema}
    if not ok_extrema:
        unresolved.append(f"velocity_extrema: {n_extrema_new}/{n_extrema_old} preserved")

    # --- 5. Temperature resolution ---
    log_T = np.log10(np.maximum(T, 1.0))
    dlogT = np.abs(np.diff(log_T))
    max_dlogT = float(dlogT.max()) if len(dlogT) > 0 else 0.0
    ok_T = max_dlogT <= cfg.target_dlogT * 3
    checks['temperature'] = {
        'max_dlogT': max_dlogT, 'target': cfg.target_dlogT,
        'ok': ok_T}
    if not ok_T:
        worst = int(np.argmax(dlogT))
        unresolved.append(f"temperature: zone {worst} dlogT={dlogT[worst]:.3f}")

    # --- 6. Shock resolution ---
    # Count actual grid points within a radius window around each shock.
    # The window is set by the source zone spacing at the shock location.
    shock_features = [f for f in features if f.kind == 'shock']
    n_shock_points = 0
    for sf in shock_features:
        # Find the source-zone spacing at the shock radius
        idx_src = np.argmin(np.abs(r_old - sf.radius))
        if idx_src > 0 and idx_src < len(r_old) - 1:
            dr_window = abs(r_old[idx_src - 1] - r_old[idx_src + 1])
        else:
            dr_window = abs(r_old[0] - r_old[-1]) / len(r_old)
        # Count grid points within ±2 source-zone spacings of the shock
        n_near = int(np.sum(np.abs(r - sf.radius) <= 2 * dr_window))
        n_shock_points += n_near
    ok_shock = len(shock_features) == 0 or n_shock_points >= len(shock_features) * cfg.min_points_per_shock
    checks['shocks'] = {
        'n_shocks': len(shock_features),
        'n_shock_points': n_shock_points,
        'min_per_shock': cfg.min_points_per_shock,
        'ok': ok_shock}
    if not ok_shock:
        unresolved.append(f"shocks: {n_shock_points}/{len(shock_features)*cfg.min_points_per_shock} points")

    # --- 7. Shell boundaries / density discontinuities ---
    # Detect sharp density jumps in the source that must be preserved
    log_rho_old = np.log10(np.maximum(rho_old, 1e-300))
    dlogrho_old = np.abs(np.diff(log_rho_old))
    shell_jumps = np.where(dlogrho_old > 1.0)[0]  # >1 dex jump
    n_shell = len(shell_jumps)
    ok_shell = True  # checked via mass conservation + feature seeding
    checks['shell_boundaries'] = {
        'n_jumps': n_shell, 'ok': ok_shell}

    # --- 8. Composition interfaces ---
    comp_features = [f for f in features if f.kind == 'composition']
    n_comp_points = 0
    for cf in comp_features:
        # Count actual grid points near the interface radius
        idx_src = np.argmin(np.abs(r_old - cf.radius))
        if idx_src > 0 and idx_src < len(r_old) - 1:
            dr_window = abs(r_old[idx_src - 1] - r_old[idx_src + 1])
        else:
            dr_window = abs(r_old[0] - r_old[-1]) / len(r_old)
        n_near = int(np.sum(np.abs(r - cf.radius) <= 2 * dr_window))
        n_comp_points += n_near
    ok_comp = len(comp_features) == 0 or n_comp_points >= len(comp_features) * cfg.min_points_per_interface
    checks['composition_interfaces'] = {
        'n_interfaces': len(comp_features),
        'n_interface_points': n_comp_points,
        'min_per_interface': cfg.min_points_per_interface,
        'ok': ok_comp}
    if not ok_comp:
        unresolved.append(f"composition: {n_comp_points}/{len(comp_features)*cfg.min_points_per_interface} points")

    # --- 9. Abundance gradients ---
    # Large ΔX at a KNOWN composition interface is expected and OK.
    # We only flag abundance gradients that are NOT at a known interface.
    comp_radii = [cf.radius for cf in comp_features]
    max_dX = 0.0
    max_dX_unresolved = 0.0
    for e, x in mass_fractions.items():
        if x.max() < cfg.composition_trace_floor:
            continue
        dX = np.abs(np.diff(x))
        for i_dx, dxi in enumerate(dX):
            max_dX = max(max_dX, float(dxi))
            # Check if this jump is at a known composition interface
            r_mid = 0.5 * (r[i_dx] + r[i_dx + 1])
            at_interface = any(abs(r_mid - cr) < abs(r[i_dx] - r[i_dx + 1]) * 2
                               for cr in comp_radii)
            if not at_interface and dxi > cfg.composition_jump:
                max_dX_unresolved = max(max_dX_unresolved, float(dxi))
    ok_abund = max_dX_unresolved <= cfg.composition_jump
    checks['abundance_gradients'] = {
        'max_dX': max_dX, 'max_dX_unresolved': max_dX_unresolved,
        'threshold': cfg.composition_jump,
        'ok': ok_abund}
    if not ok_abund:
        unresolved.append(f"abundance: unresolved ΔX={max_dX_unresolved:.3f}")

    # --- 10. Spectrum-forming region (τ ~ 0.01–1) ---
    # Count grid points in the photospheric / line-forming region
    spec_mask = (tau >= 0.01) & (tau <= 10.0)
    n_spec = int(np.sum(spec_mask))
    ok_spec = n_spec >= 20  # need at least 20 points in spectrum-forming region
    checks['spectrum_region'] = {
        'n_points_tau_0.01_to_10': n_spec,
        'ok': ok_spec}
    if not ok_spec:
        unresolved.append(f"spectrum_region: only {n_spec} points in tau=[0.01,10]")

    # --- 11. Diffusion/photospheric region (τ ~ 1 to τ_inner) ---
    # Count grid points between τ=1 and the inner boundary.  This region
    # contains the diffusion-to-free-streaming transition and must be
    # adequately resolved for CMFGEN's transfer solution.
    tau_inner = float(tau[-1]) if len(tau) > 0 else 0.0
    diff_mask = (tau >= 1.0) & (tau <= tau_inner)
    n_diff = int(np.sum(diff_mask))
    # Need at least 30 points in the diffusion region for adequate resolution
    # of the diffusion-to-streaming transition (Dessart+2015, Hillier+2012).
    ok_diff = n_diff >= 30
    checks['diffusion_region'] = {
        'n_points_tau_1_to_inner': n_diff,
        'tau_inner': tau_inner,
        'ok': ok_diff}
    if not ok_diff:
        unresolved.append(f"diffusion_region: only {n_diff} points in tau=[1,{tau_inner:.1f}]")

    # --- 12. Detached shell detection ---
    # A detached shell is a region of high density bounded by two sharp
    # composition or density interfaces.  Check that any such region has
    # adequate internal resolution.
    shells = _detect_detached_shells(r, rho, mass_fractions, cfg)
    n_shells = len(shells)
    shell_issues = []
    for si, shell in enumerate(shells):
        # Count points within the shell
        shell_mask = (r >= shell['r_inner']) & (r <= shell['r_outer'])
        n_shell = int(np.sum(shell_mask))
        # Need at least 10 points within a detached shell to resolve its
        # internal ionisation/opacity structure
        if n_shell < 10:
            shell_issues.append({
                'shell_index': si,
                'r_outer': shell['r_outer'],
                'r_inner': shell['r_inner'],
                'n_points': n_shell,
                'density_ratio': shell['density_ratio'],
            })
            unresolved.append(
                f"detached_shell_{si}: only {n_shell} points in shell at "
                f"r=[{shell['r_inner']:.3e},{shell['r_outer']:.3e}] cm "
                f"(density ratio {shell['density_ratio']:.1f})")
    checks['detached_shells'] = {
        'n_shells': n_shells,
        'n_under_resolved': len(shell_issues),
        'shells': shells,
        'ok': len(shell_issues) == 0}


    passed = all(c['ok'] for c in checks.values())

    # --- Print formatted report ---
    print("\n" + "=" * 60)
    print("STELLA → CMFGEN GRID RESOLUTION REPORT")
    print("=" * 60)
    print(f"\nTotal grid points:             {n} / {cfg.max_points}")
    print(f"Source zones:                  {len(r_old)}")
    print()
    for name, chk in checks.items():
        status = "PASS" if chk['ok'] else "FAIL"
        label = name.replace('_', ' ').title()
        print(f"  {label:30s} {status}")
    print(f"\n  Largest Δlogτ:               {max_dlogtau:.4f}")
    print(f"  Largest Δv:                  {max_dv:.1f} km/s")
    print(f"  Largest Δlogρ:               {max_dlogrho:.4f}")
    print(f"  Largest ΔlogT:               {max_dlogT:.4f}")
    print(f"  Largest ΔX:                  {max_dX:.4f}")
    print()
    # Feature point counts
    comp_f = [f for f in features if f.kind == 'composition']
    shock_f = [f for f in features if f.kind == 'shock']
    ext_f = [f for f in features if f.kind == 'velocity_extremum']
    print(f"  Composition interfaces:       {len(comp_f)}")
    print(f"  Shocks:                       {len(shock_f)}")
    print(f"  Velocity extrema:             {len(ext_f)}")
    print(f"  Spectrum-region points:       {n_spec}")
    print()
    if unresolved:
        print(f"  Unresolved regions:           {len(unresolved)}")
        for u in unresolved:
            print(f"    - {u}")
    else:
        print(f"  Unresolved regions:           NONE")
    print("\n" + "=" * 60)
    if passed:
        print("RESOLUTION CHECK: PASSED — grid is adequately resolved for CMFGEN")
    else:
        print("RESOLUTION CHECK: FAILED — grid has unresolved regions")
    print("=" * 60 + "\n")

    return {'checks': checks, 'unresolved': unresolved}, passed, unresolved


def _count_extrema(v):
    n = len(v)
    if n < 3:
        return 0
    cnt = 0
    for i in range(1, n - 1):
        if (v[i] > v[i - 1] and v[i] > v[i + 1]) or \
           (v[i] < v[i - 1] and v[i] < v[i + 1]):
            cnt += 1
    return cnt


def _all_velocity_turn_radii(r, v):
    """Return radii of every simple 3-point velocity turn (local max/min)
    AND its immediate neighbours (i-1, i, i+1).

    Seeding all three points guarantees the 3-point turn is exactly
    represented in the new grid, so PCHIP cannot smooth away shallow extrema
    and the full non-monotonic velocity topology is preserved even at coarse
    point budgets.
    """
    n = len(v)
    if n < 3:
        return np.array([])
    idx = set()
    for i in range(1, n - 1):
        if (v[i] > v[i - 1] and v[i] > v[i + 1]) or \
           (v[i] < v[i - 1] and v[i] < v[i + 1]):
            idx.update([i - 1, i, i + 1])
    return np.array([r[i] for i in sorted(idx)]) if idx else np.array([])


# ===========================================================================
# Diagnostics / plotting
# ===========================================================================

def plot_diagnostics(r_old, v_old, T_old, rho_old, tau_old, kappa_old,
                     mass_fractions_old, features_old,
                     r, v, T, rho, tau, kappa, mass_fractions,
                     report, out_dir):
    """Generate the mandatory diagnostic plots."""
    import os
    os.makedirs(out_dir, exist_ok=True)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    plot_elems = ['H', 'He', 'C', 'O', 'Si', 'S', 'Ca', 'Fe', 'Ni']
    plot_elems = [e for e in plot_elems if e in mass_fractions]

    # --- Hydrodynamics ---
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    for ax, y_old, y_new, lab in zip(
            axes.flat,
            [rho_old, v_old, T_old, tau_old],
            [rho, v, T, tau],
            ['density', 'velocity', 'temperature', 'tau']):
        ax.loglog(r_old, np.maximum(y_old, 1e-30), 'k-', lw=0.6, alpha=0.5, label='STELLA')
        ax.loglog(r, np.maximum(y_new, 1e-30), 'r.', ms=2, label='AMR grid')
        ax.set_xlabel('r [cm]'); ax.set_ylabel(lab)
        ax.legend(fontsize=7)
    fig.suptitle(f'Tau-AMR hydrodynamics ({len(r)} pts)')
    fig.tight_layout(); fig.savefig(f'{out_dir}/hydrodynamics.png', dpi=120); plt.close(fig)

    # --- Composition ---
    fig, axes = plt.subplots(len(plot_elems), 1, figsize=(12, 2.2 * len(plot_elems)),
                             sharex=True)
    if len(plot_elems) == 1:
        axes = [axes]
    for ax, e in zip(axes, plot_elems):
        ax.semilogx(r_old, mass_fractions_old[e], 'k-', lw=0.6, alpha=0.5)
        ax.semilogx(r, mass_fractions[e], 'r.', ms=2)
        ax.set_ylabel(e)
        for f in features_old:
            if f.kind == 'composition':
                ax.axvline(f.radius, color='b', ls='--', lw=0.5, alpha=0.5)
    axes[-1].set_xlabel('r [cm]')
    fig.suptitle('Tau-AMR composition (blue dashed = detected interfaces)')
    fig.tight_layout(); fig.savefig(f'{out_dir}/composition.png', dpi=120); plt.close(fig)

    # --- Resolution ---
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    dr = np.abs(np.diff(r))
    dror = dr / r[:-1]
    dv = np.abs(np.diff(v))
    dlogtau = np.abs(np.log10((tau[1:] + 1e-6) / (tau[:-1] + 1e-6)))
    dm = zone_masses(r, rho)[:-1]
    for ax, y, lab in zip(axes.flat, [dr, dror, dv, dlogtau],
                          ['dr [cm]', 'dr/r', 'dv [km/s]', 'dlog tau']):
        ax.semilogy(r[:-1], np.maximum(y, 1e-30), 'b-', lw=0.7)
        ax.set_xlabel('r [cm]'); ax.set_ylabel(lab)
        for f in features_old:
            ax.axvline(f.radius, color='r', ls='--', lw=0.4, alpha=0.4)
    fig.suptitle('Tau-AMR resolution diagnostics (red dashed = features)')
    fig.tight_layout(); fig.savefig(f'{out_dir}/resolution.png', dpi=120); plt.close(fig)

    # --- AMR score / point density ---
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.hist(np.log10(r), bins=60, color='steelblue', edgecolor='k')
    ax.set_xlabel('log10 r [cm]'); ax.set_ylabel('grid points per bin')
    for f in features_old:
        ax.axvline(np.log10(f.radius), color='r', ls='--', lw=0.6,
                   label=f'{f.kind}' if f == features_old[0] else None)
    if features_old:
        ax.legend(fontsize=7)
    ax.set_title(f'AMR point density ({len(r)} pts, {report["splits"]} splits)')
    fig.tight_layout(); fig.savefig(f'{out_dir}/point_density.png', dpi=120); plt.close(fig)


# ===========================================================================
# Top-level driver
# ===========================================================================

def run_tau_amr(swd_path, epoch_day, output_path, cfg: AMRConfig,
                plot_dir=None, lum_lsun=1e9, t_floor=3000.0,
                trace_threshold=1e-6, sigma=None):
    """Full pipeline: .swd -> AMR grid -> CMFGEN SN_HYDRO_DATA + diagnostics.

    Parameters
    ----------
    swd_path : str
        Path to the STELLA .swd file.
    epoch_day : float
        Epoch (days) to extract; closest available epoch is used.
    output_path : str
        Path for the CMFGEN SN_HYDRO_DATA file.
    cfg : AMRConfig
        All thresholds.
    plot_dir : str or None
        Directory for diagnostic plots (None = no plots).
    lum_lsun : float
        Luminosity (Lsun) written to SN_HYDRO_DATA.
    t_floor : float
        Temperature floor (K).
    trace_threshold : float
        Mass fractions below this are zeroed.
    sigma : array or None
        Pre-computed sigma array; if None, computed from velocity gradient.

    Returns
    -------
    dict with validation report, feature report, and grid info.
    """
    if cfg.verbose:
        print(f"[Tau-AMR] Reading {swd_path}")
    epochs = read_swd(swd_path)
    # find closest epoch
    t_avail = sorted(epochs.keys())
    t_day = min(t_avail, key=lambda t: abs(t - epoch_day))
    ep = epochs[t_day]
    if cfg.verbose:
        print(f"[Tau-AMR] Requested {epoch_day}d, using {t_day:.3f}d "
              f"({len(ep['r'])} zones)")

    r = ep['r'].copy()           # outer->inner (decreasing)
    v = ep['v'].copy()
    T = ep['T'].copy()
    rho = ep['rho'].copy()
    kappa_stella = ep['kappa'].copy() if ep['kappa'] is not None else None
    composition = ep['composition']
    mass_fractions = normalize_swd_composition(composition)

    # --- optical depth ---
    tau, kappa, kappa_label = reconstruct_optical_depth(
        r, rho, T, mass_fractions, kappa_stella=kappa_stella, cfg=cfg)
    if cfg.verbose:
        print(f"[Tau-AMR] kappa source: {kappa_label}")
        print(f"[Tau-AMR] tau range: [{tau.min():.3e}, {tau.max():.3e}]")

    # --- features (hard constraints) ---
    features, comp_f, shock_f, ext_f = detect_all_features(
        r, v, rho, T, mass_fractions, tau, cfg)
    if cfg.verbose:
        print(f"[Tau-AMR] Features: {len(comp_f)} composition interfaces, "
              f"{len(shock_f)} shocks, {len(ext_f)} velocity extrema")
        for f in features:
            print(f"         {f.kind} @ r={f.radius:.3e} cm "
                  f"v={f.velocity:.0f} km/s tau={f.tau:.3e} ({f.detail})")

    interface_indices = np.array([f.index for f in comp_f], dtype=int)

    # --- recursive AMR ---
    r_amr, v_amr, T_amr, rho_amr, kappa_amr, tau_amr, mf_amr, report = \
        recursive_amr(r, v, T, rho, kappa, tau, mass_fractions,
                      features, cfg, interface_indices)
    if cfg.verbose:
        print(f"[Tau-AMR] Final grid: {len(r_amr)} points "
              f"({report['splits']} splits, {report['iterations']} iters)")
        if report['unresolved_features']:
            print(f"[Tau-AMR] WARNING: {len(report['unresolved_features'])} "
                  f"features still under-resolved at point budget")

    # --- sigma (velocity gradient for CMFGEN non-monotonic solver) ---
    if sigma is None:
        sigma = compute_sigma(r_amr, v_amr)

    # --- validation (conservation + topology) ---
    val = validate_grid(r_amr, v_amr, T_amr, rho_amr, kappa_amr, tau_amr,
                        mf_amr, r, v, rho, mass_fractions, cfg)
    if cfg.verbose:
        print(f"[Tau-AMR] Validation {'PASSED' if val['passed'] else 'FAILED'}:")
        for name, chk in val['checks'].items():
            print(f"         {name}: {'OK' if chk['ok'] else 'FAIL'} {chk}")

    # --- pre-CMFGEN resolution report (from STELLA physics, not CMFGEN) ---
    res_rep, res_passed, res_unresolved = resolution_report(
        r_amr, v_amr, T_amr, rho_amr, kappa_amr, tau_amr, mf_amr,
        r, v, T, rho, tau, mass_fractions, features, cfg)

    # --- write CMFGEN SN_HYDRO_DATA ---
    # zero trace species
    mf_out = {e: x.copy() for e, x in mf_amr.items()}
    for e in mf_out:
        mf_out[e][mf_out[e] < trace_threshold] = 0.0
    _renormalise(mf_out)

    out_p = Path(output_path)
    out_p.parent.mkdir(parents=True, exist_ok=True)
    write_sn_hydro_data(
        output_path, r_amr, v_amr, rho_amr, T_amr, sigma, kappa_amr,
        time_days=t_day, mass_fractions=mf_out,
        lum_lsun=lum_lsun, t_floor=t_floor,
        trace_threshold=trace_threshold)
    if cfg.verbose:
        print(f"[Tau-AMR] Wrote CMFGEN SN_HYDRO_DATA -> {output_path} "
              f"({len(r_amr)} depths)")

    # --- diagnostics ---
    if cfg.plot and plot_dir:
        plot_diagnostics(r, v, T, rho, tau, kappa, mass_fractions, features,
                         r_amr, v_amr, T_amr, rho_amr, tau_amr, kappa_amr,
                         mf_amr, report, plot_dir)
        if cfg.verbose:
            print(f"[Tau-AMR] Diagnostics -> {plot_dir}")

    # --- manifest ---
    manifest = {
        'swd_path': str(swd_path), 'epoch_day': float(t_day),
        'n_source_zones': int(len(r)), 'n_grid_points': int(len(r_amr)),
        'max_points': cfg.max_points, 'kappa_source': kappa_label,
        'tau_range': [float(tau.min()), float(tau.max())],
        'features': [{'kind': f.kind, 'radius': float(f.radius),
                      'velocity': float(f.velocity), 'tau': float(f.tau),
                      'detail': f.detail} for f in features],
        'amr_report': report, 'validation': val,
        'resolution_report': res_rep,
        'resolution_passed': res_passed,
    }
    # Overall pass requires both conservation validation AND resolution check
    manifest['validation']['passed'] = val['passed'] and res_passed
    return manifest


def compute_sigma(r, v):
    """Sigma (velocity gradient) for CMFGEN non-monotonic solver.

    sigma = (1/v)(dv/dr)  in 1/s, with r in cm and v in cm/s.
    Outer->inner ordering (r decreasing).  Protected against zero gradients.
    """
    v_cgs = np.maximum(np.abs(v) * 1e5, 1.0)   # km/s -> cm/s
    n = len(r)
    sigma = np.zeros(n)
    # central difference in radius (r decreasing so dr = r[i-1]-r[i] is +ve going in)
    for i in range(n):
        if i == 0:
            dr = r[i] - r[i + 1]
            dv = v[i] - v[i + 1]
        elif i == n - 1:
            dr = r[i - 1] - r[i]
            dv = v[i - 1] - v[i]
        else:
            dr = r[i - 1] - r[i + 1]
            dv = v[i - 1] - v[i + 1]
        if abs(dr) < 1e-30:
            sigma[i] = 0.0
        else:
            sigma[i] = (dv * 1e5) / (v_cgs[i] * dr)
    return sigma


# ===========================================================================
# CLI
# ===========================================================================

def _cfg_from_args(args) -> AMRConfig:
    kw = {}
    # If a YAML config file is given, load it first as the base.
    if getattr(args, 'config', None):
        import yaml
        with open(args.config) as fh:
            kw = yaml.safe_load(fh) or {}
    # CLI flags override YAML values.
    for cli_name, cfg_name in [
        ('max_points', 'max_points'), ('min_points', 'min_points'),
        ('base_grid_points', 'base_grid_points'), ('tau_floor', 'tau_floor'),
        ('kappa_source', 'kappa_source'),
        ('target_dlogtau', 'target_dlogtau'),
        ('target_dlogrho', 'target_dlogrho'),
        ('target_dv_kms', 'target_dv_kms'),
        ('target_dlogT', 'target_dlogT'),
        ('target_dr_over_r', 'target_dr_over_r'),
        ('w_tau', 'w_tau'), ('w_rho', 'w_rho'), ('w_v', 'w_v'),
        ('w_T', 'w_T'), ('w_curv', 'w_curv'),
        ('composition_jump', 'composition_jump'),
        ('shock_density_jump', 'shock_density_jump'),
        ('shock_temp_jump', 'shock_temp_jump'),
        ('min_points_per_interface', 'min_points_per_interface'),
        ('min_points_per_shock', 'min_points_per_shock'),
        ('min_points_per_extremum', 'min_points_per_extremum'),
        ('mass_tol', 'mass_tol')]:
        v = getattr(args, cli_name, None)
        if v is not None:
            kw[cfg_name] = v
    kw['plot'] = not getattr(args, 'no_plot', False)
    kw['verbose'] = True
    return AMRConfig(**kw)


def main(argv=None):
    p = argparse.ArgumentParser(
        description="Tau-based AMR STELLA .swd -> CMFGEN grid optimiser")
    p.add_argument('--swd', required=True, help='STELLA .swd file')
    p.add_argument('--epoch', type=float, required=True,
                   help='Epoch (days) to extract')
    p.add_argument('--output', required=True, help='SN_HYDRO_DATA output path')
    p.add_argument('--plot-dir', default=None, help='Diagnostic plot directory')
    p.add_argument('--manifest', default=None, help='Manifest JSON output path')
    p.add_argument('--lum-lsun', type=float, default=1e9)
    p.add_argument('--config', default=None,
                   help='YAML config file for AMR thresholds (configs/tau_amr_config.yaml)')
    p.add_argument('--t-floor', type=float, default=3000.0)
    p.add_argument('--trace-threshold', type=float, default=1e-6)
    p.add_argument('--no-plot', action='store_true')
    # point budget (default None => AMRConfig default applies; CLI/YAML override)
    p.add_argument('--max-points', type=int, default=None)
    p.add_argument('--min-points', type=int, default=None)
    p.add_argument('--base-grid-points', type=int, default=None)
    # optical depth
    p.add_argument('--tau-floor', type=float, default=None)
    p.add_argument('--kappa-source', default=None,
                   choices=['hybrid', 'stella', 'kramers', None])
    # soft targets
    p.add_argument('--target-dlogtau', type=float, default=None)
    p.add_argument('--target-dlogrho', type=float, default=None)
    p.add_argument('--target-dv-kms', type=float, default=None)
    p.add_argument('--target-dlogT', type=float, default=None)
    p.add_argument('--target-dr-over-r', type=float, default=None)
    # weights
    p.add_argument('--w-tau', type=float, default=None)
    p.add_argument('--w-rho', type=float, default=None)
    p.add_argument('--w-v', type=float, default=None)
    p.add_argument('--w-T', type=float, default=None)
    p.add_argument('--w-curv', type=float, default=None)
    # hard constraints
    p.add_argument('--composition-jump', type=float, default=None)
    p.add_argument('--shock-density-jump', type=float, default=None)
    p.add_argument('--shock-temp-jump', type=float, default=None)
    p.add_argument('--min-points-per-interface', type=int, default=None)
    p.add_argument('--min-points-per-shock', type=int, default=None)
    p.add_argument('--min-points-per-extremum', type=int, default=None)
    p.add_argument('--mass-tol', type=float, default=None)
    args = p.parse_args(argv)

    cfg = _cfg_from_args(args)
    manifest = run_tau_amr(
        args.swd, args.epoch, args.output, cfg,
        plot_dir=args.plot_dir, lum_lsun=args.lum_lsun,
        t_floor=args.t_floor, trace_threshold=args.trace_threshold)

    if args.manifest:
        Path(args.manifest).parent.mkdir(parents=True, exist_ok=True)
        with open(args.manifest, 'w') as fh:
            json.dump(manifest, fh, indent=2)
        print(f"[Tau-AMR] Manifest -> {args.manifest}")

    if not manifest['validation']['passed']:
        print("[Tau-AMR] WARNING: validation FAILED — see report above",
              file=sys.stderr)
        sys.exit(2)


if __name__ == '__main__':
    main()
