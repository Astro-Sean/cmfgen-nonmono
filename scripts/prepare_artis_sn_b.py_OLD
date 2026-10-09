#!/usr/bin/env python3
"""Prepare ARTIS input files from a SN_B STELLA .swd snapshot.

This script uses the SN_B STELLA model directly (no homologous expansion)
at the desired snapshot epoch (~40d after shell launch).  The key improvement
over the earlier test_artis_blackbody_illuminated_ejecta.py approach is:

1. Uses the actual STELLA gas temperature profile (not an artificial T_BB).
2. Computes the initial thermal energy (q column) from the STELLA luminosity
   profile — energy is deposited where STELLA says it is generated.
3. Uses the SN_B interaction profile, which already contains the shock
   interaction physics at r ~ 2e15 cm (near 1e15 cm).
4. Maps the full 15-element STELLA composition to the 7 ARTIS-available
   elements (H, He, C, O, Si, Ca, Fe).
5. Uses adaptive grid downsampling that preserves the interaction region,
   composition interfaces, and temperature gradients.

The ARTIS "q" column provides initial thermal energy per gram at t_model.
With USE_MODEL_INITIAL_ENERGY=true and INITIAL_PACKETS_ON=true, ARTIS deposits
this energy as k-packets at tmin.  These k-packets thermalize and emerge as
r-packets with a Planck spectrum at the local electron temperature T_e.

Usage:
    python3 scripts/prepare_artis_sn_b.py [--swd SWD_PATH] [--target-age DAYS]
        [--work-dir DIR] [--nrank N] [--mpkts M] [--run]

    python3 scripts/prepare_artis_sn_b.py --dry-run  # generate inputs only
"""

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
SIGMA_SB = 5.670374419e-5      # Stefan-Boltzmann constant [erg cm^-2 s^-1 K^-4]
C_LIGHT = 2.99792458e10        # speed of light [cm/s]
PI = np.pi
M_SUN = 1.989e33               # solar mass [g]
DAY_S = 86400.0                # day in seconds [s]
KB = 1.380649e-16              # Boltzmann constant [erg/K]
H_PLANCK = 6.62607015e-27      # Planck constant [erg*s]

# Default SN_B model
# Best-matching SN_A model (used as source for X/SN_B grid)
DEFAULT_SN_A_SWD = "/viper/ptmp1/sbrennan/SN2024pba/SN_A/SN_A_g_T20k_Ni06/stella/SN_A_g_T20k_Ni06.swd"
ARTIS_REPO = "/viper/ptmp1/sbrennan/artis_repo"
ARTIS_EXE = f"{ARTIS_REPO}/sn3d"
ARTIS_ATOMIC = "/viper/ptmp1/sbrennan/artisatomic/artis_files"

# STELLA 15-element composition
STELLA_ELEMENTS = ['H', 'He', 'C', 'N', 'O', 'Ne', 'Na', 'Mg', 'Al',
                   'Si', 'S', 'Ar', 'Ca', 'Fe', 'Ni']

# ARTIS uses Z=1..30 for abundances.txt
ELEMENT_Z = {'H': 1, 'He': 2, 'C': 6, 'N': 7, 'O': 8, 'Ne': 10,
             'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'S': 16,
             'Ar': 18, 'Ca': 20, 'Fe': 26, 'Ni': 28}

# ARTIS atomic data supports these elements (from compositiondata.txt):
# H(Z=1), He(Z=2), C(Z=6), O(Z=8), Si(Z=14), Ca(Z=20), Fe(Z=26)
ARTIS_ELEMENTS = {'H', 'He', 'C', 'O', 'Si', 'Ca', 'Fe'}

# Map STELLA elements not in ARTIS to closest available element
# N->C, Ne->O, Na->Si, Mg->Si, Al->Si, S->Si, Ar->Ca, Ni->Fe
ELEMENT_REMAP = {
    'N': 'C',
    'Ne': 'O',
    'Na': 'Si',
    'Mg': 'Si',
    'Al': 'Si',
    'S': 'Si',
    'Ar': 'Ca',
    'Ni': 'Fe',
}

# Default parameters
DEFAULT_SN_A_AGE_D = 400.0   # homologously evolve SN_A to 400d
DEFAULT_T_MAX_D = 50.0       # ARTIS simulation duration (days)
DEFAULT_N_TIMESTEPS = 50
DEFAULT_N_LTE = 5
DEFAULT_MAX_CELLS = 400
DEFAULT_N_RANK = 8
DEFAULT_MPKTS = 500000


# ---------------------------------------------------------------------------
# Step 1: Read SN_A .swd and homologously evolve to target age
# ---------------------------------------------------------------------------
def read_and_evolve_sn_a(swd_path, target_age_d):
    """Read the last epoch from the SN_A .swd file and homologously evolve
    it to target_age_d.

    Homologous expansion:
      r(t) = r(t0) * (t / t0)
      v(t) = v(t0)  (constant velocity)
      rho(t) = rho(t0) * (t0 / t)^3
      T(t) = T(t0) * (t0 / t)  (adiabatic)

    The luminosity from Ni56 decay is kept at the STELLA value at t0.
    For ARTIS, the q column provides the initial thermal energy; the
    Ni56 decay energy is handled separately by ARTIS through the
    X_Ni56 abundance (set to trace 1e-10 in model.txt).

    Returns the evolved epoch dict and the target time.
    """
    sys.path.insert(0, str(Path(__file__).parent))
    from convert_stella_swd_to_cmfgen import read_swd, normalize_swd_composition

    data = read_swd(swd_path)
    times = sorted([float(k) for k in data.keys()])
    last_t = times[-1]
    epoch = data[last_t]

    # Normalize composition
    comp_dict = normalize_swd_composition(epoch['composition'])
    comp = np.zeros((len(epoch['r']), len(STELLA_ELEMENTS)))
    for i, elem in enumerate(STELLA_ELEMENTS):
        comp[:, i] = comp_dict[elem]

    # Compute homologous velocity at snapshot time
    t_snap_s = last_t * DAY_S
    v_snap = epoch['r'] / t_snap_s  # cm/s

    print(f"Read SN_A .swd: {swd_path}")
    print(f"  Last epoch: t={last_t:.4f} d")
    print(f"  Zones: {len(epoch['r'])}")
    print(f"  R: {epoch['r'].min():.3e} - {epoch['r'].max():.3e} cm")
    print(f"  V: {v_snap.min()/1e5:.1f} - {v_snap.max()/1e5:.1f} km/s")
    print(f"  T_gas: {epoch['T'].min():.1f} - {epoch['T'].max():.1f} K")
    print(f"  T_rad: {epoch['T_rad'].min():.1f} - {epoch['T_rad'].max():.1f} K")
    print(f"  rho: {epoch['rho'].min():.3e} - {epoch['rho'].max():.3e} g/cm^3")
    print(f"  L_max: {epoch['flux'].max():.3e} erg/s")

    # Homologous evolution to target_age_d
    scale = target_age_d / last_t
    r_evolved = epoch['r'] * scale
    rho_evolved = epoch['rho'] / scale**3
    v_evolved = v_snap  # velocity is constant in homologous expansion
    T_evolved = epoch['T'] * (last_t / target_age_d)  # adiabatic cooling
    T_rad_evolved = epoch['T_rad'] * (last_t / target_age_d)
    # Luminosity: keep STELLA value (Ni56 decay luminosity at t0)
    # ARTIS will handle the energy deposition via the q column
    flux_evolved = epoch['flux']

    print(f"\nHomologous expansion to {target_age_d:.0f} d (scale={scale:.3f}):")
    print(f"  R: {r_evolved.min():.3e} - {r_evolved.max():.3e} cm")
    print(f"  V: {v_evolved.min()/1e5:.1f} - {v_evolved.max()/1e5:.1f} km/s")
    print(f"  T_gas: {T_evolved.min():.1f} - {T_evolved.max():.1f} K (adiabatic)")
    print(f"  T_rad: {T_rad_evolved.min():.1f} - {T_rad_evolved.max():.1f} K")
    print(f"  rho: {rho_evolved.min():.3e} - {rho_evolved.max():.3e} g/cm^3")

    # Build evolved epoch dict
    evolved = {
        'r': r_evolved,
        'v_cms': v_evolved,
        'rho': rho_evolved,
        'T': T_evolved,
        'T_rad': T_rad_evolved,
        'flux': flux_evolved,
        'composition_normalized': comp,
    }

    return evolved, target_age_d


# ---------------------------------------------------------------------------
# Step 2: Map 15-element STELLA composition to ARTIS-available elements
# ---------------------------------------------------------------------------
def map_composition_to_artis(comp_stella):
    """Map 15-element STELLA composition to 7 ARTIS-available elements.

    Elements not in ARTIS are remapped to the closest available element:
      N -> C, Ne -> O, Na/Mg/Al/S -> Si, Ar -> Ca, Ni -> Fe

    Returns an array of shape (nzon, 15) with only ARTIS elements non-zero,
    renormalized to sum=1.
    """
    n = comp_stella.shape[0]
    comp_artis = np.zeros((n, len(STELLA_ELEMENTS)))

    for i, elem in enumerate(STELLA_ELEMENTS):
        if elem in ARTIS_ELEMENTS:
            comp_artis[:, i] = comp_stella[:, i]
        elif elem in ELEMENT_REMAP:
            target = ELEMENT_REMAP[elem]
            target_idx = STELLA_ELEMENTS.index(target)
            comp_artis[:, target_idx] += comp_stella[:, i]

    # Renormalize
    totals = comp_artis.sum(axis=1, keepdims=True)
    totals = np.where(totals < 1e-30, 1.0, totals)
    comp_artis = comp_artis / totals

    # Report
    print("\nComposition mapping (STELLA -> ARTIS):")
    for elem in STELLA_ELEMENTS:
        if elem in ARTIS_ELEMENTS:
            print(f"  {elem}: kept (max={comp_stella[:, STELLA_ELEMENTS.index(elem)].max():.4f})")
        elif elem in ELEMENT_REMAP:
            print(f"  {elem} -> {ELEMENT_REMAP[elem]} (max={comp_stella[:, STELLA_ELEMENTS.index(elem)].max():.4f})")

    return comp_artis


# ---------------------------------------------------------------------------
# Step 3: Adaptive grid downsampling
# ---------------------------------------------------------------------------
def adaptive_downsample(r, rho, v, comp, T_gas, T_rad, flux, max_cells):
    """Downsample the grid using adaptive sampling.

    Preserves:
    - The interaction region (high T, high L)
    - Composition interfaces
    - Temperature gradients
    - The inner and outer boundaries

    .swd is outer-to-inner; ARTIS model.txt is inner-to-outer.
    Returns arrays in inner-to-outer order (increasing r).
    """
    n_orig = len(r)

    # First reverse everything to inner-to-outer (increasing r)
    r_io = r[::-1]
    rho_io_all = rho[::-1]
    v_io_all = v[::-1]
    comp_io_all = comp[::-1]
    T_io_all = T_gas[::-1]
    Trad_io_all = T_rad[::-1]
    flux_io_all = flux[::-1]

    if n_orig <= max_cells:
        print(f"\nNo downsampling needed: {n_orig} zones")
        return (r_io, rho_io_all, v_io_all, comp_io_all,
                T_io_all, Trad_io_all, flux_io_all)

    # Compute importance weights on the inner-to-outer array
    # 1. Temperature gradient (preserve hot/cold transitions)
    log_T = np.log10(np.maximum(T_io_all, 1.0))
    dlog_T = np.abs(np.diff(log_T, prepend=log_T[0]))

    # 2. Luminosity gradient (preserve energy generation regions)
    log_L = np.log10(np.maximum(np.abs(flux_io_all), 1.0))
    dlog_L = np.abs(np.diff(log_L, prepend=log_L[0]))

    # 3. Density gradient (preserve shell boundaries)
    log_rho = np.log10(np.maximum(rho_io_all, 1e-30))
    dlog_rho = np.abs(np.diff(log_rho, prepend=log_rho[0]))

    # 4. Composition gradient (preserve abundance discontinuities)
    comp_grad = np.zeros(n_orig)
    for i in range(1, n_orig):
        comp_grad[i] = np.sum(np.abs(comp_io_all[i] - comp_io_all[i-1]))

    # 5. Velocity gradient (preserve shock fronts)
    log_v = np.log10(np.maximum(v_io_all, 1.0))
    dlog_v = np.abs(np.diff(log_v, prepend=log_v[0]))

    # Normalize each indicator to [0, 1]
    def norm(arr):
        mx = arr.max()
        return arr / mx if mx > 0 else arr

    weights = (0.30 * norm(dlog_T) +
               0.25 * norm(dlog_L) +
               0.15 * norm(dlog_rho) +
               0.15 * norm(comp_grad) +
               0.15 * norm(dlog_v))

    # Always include first and last zones
    weights[0] = 1.0
    weights[-1] = 1.0

    # Ensure zones with T > 5000 K get high weight
    weights[T_io_all > 5000] = np.maximum(weights[T_io_all > 5000], 0.8)

    # Ensure zones with positive flux get high weight
    pos_flux = flux_io_all > 1e38
    weights[pos_flux] = np.maximum(weights[pos_flux], 0.5)

    # Convert weights to sampling probability
    prob = weights / weights.sum()

    # Sample indices using cumulative distribution
    cumprob = np.cumsum(prob)
    uniform = np.linspace(0, cumprob[-1], max_cells, endpoint=False)
    indices = np.searchsorted(cumprob, uniform)
    indices = np.clip(indices, 0, n_orig - 1)
    indices = np.unique(indices)

    # Ensure we don't exceed max_cells
    if len(indices) > max_cells:
        indices = np.linspace(0, n_orig - 1, max_cells, dtype=int)

    # Select from inner-to-outer arrays
    r_ds = r_io[indices]
    rho_ds = rho_io_all[indices]
    v_ds = v_io_all[indices]
    comp_ds = comp_io_all[indices]
    T_ds = T_io_all[indices]
    Trad_ds = Trad_io_all[indices]
    flux_ds = flux_io_all[indices]

    print(f"\nAdaptive downsampling: {n_orig} -> {len(indices)} zones")
    print(f"  R range: {r_ds.min():.3e} - {r_ds.max():.3e} cm")
    print(f"  V range: {v_ds.min()/1e5:.1f} - {v_ds.max()/1e5:.1f} km/s")
    print(f"  T range: {T_ds.min():.1f} - {T_ds.max():.1f} K")

    return r_ds, rho_ds, v_ds, comp_ds, T_ds, Trad_ds, flux_ds


# ---------------------------------------------------------------------------
# Step 4: Set temperature and compute initial energy (q) from STELLA luminosity
# ---------------------------------------------------------------------------
def set_temperature_and_energy(r_io, rho_io, T_gas_io, comp_io, flux_io,
                                t_model_d, t_max_d,
                                t_bb=7000.0, t_outer_floor=2000.0,
                                target_l_bol=1.0e42):
    """Set gas temperature and initial thermal energy (q) for the SN_B
    illumination of SN_A ejecta at 400d.

    The SN_A ejecta has been homologously evolved to 400d, so the STELLA
    temperatures are very low (~300-850 K from adiabatic cooling).  The
    SN_B shock heats the inner metal-rich zones, creating a luminous
    blackbody-like source.  We model this by:

    1. Setting the inner metal-rich zones (H < 0.01) to T_BB (SN_B shock T)
    2. Setting the outer H/He zones to max(T_adiabatic, T_outer_floor)
    3. Depositing q (initial thermal energy) in the inner zones to produce
       the target SN_B luminosity

    The q column is energy per gram at t_model.  ARTIS scales it to tmin:
        q_tmin = q_tmodel * (t_model / tmin)
    Since tmin = t_model, q_tmin = q_tmodel.

    Args:
        t_bb: Temperature of the SN_B-heated inner source [K]
        t_outer_floor: Minimum temperature for outer H/He zones [K]
        target_l_bol: Target bolometric luminosity [erg/s]
    """
    n = len(r_io)

    # Determine which zones are metal-rich (H < 0.01) — these are the
    # SN_A core zones heated by the SN_B shock
    h_idx = STELLA_ELEMENTS.index('H')
    is_metal = comp_io[:, h_idx] < 0.01

    # Find the transition zone (last metal-rich zone)
    transition_idx = np.where(is_metal)[0]
    if len(transition_idx) == 0:
        t_gas = np.full(n, t_bb)
    else:
        last_metal = transition_idx[-1]
        ramp_width = min(20, n - last_metal - 1)
        t_gas = np.full(n, t_outer_floor)
        # Outer zones: use adiabatic STELLA T, but not below floor
        t_gas[:] = np.maximum(T_gas_io, t_outer_floor)
        # Inner metal-rich zones: set to T_BB
        t_gas[:last_metal + 1] = t_bb
        # Smooth ramp at the transition
        if ramp_width > 0:
            for i in range(last_metal + 1, last_metal + 1 + ramp_width):
                if i < n:
                    frac = (i - last_metal) / ramp_width
                    t_gas[i] = t_bb * (1 - frac) + max(T_gas_io[i], t_outer_floor) * frac

    # Compute mass of each zone
    r_edges = np.concatenate([[0], r_io])
    vol = 4.0/3.0 * PI * (r_edges[1:]**3 - r_edges[:-1]**3)
    mass = vol * rho_io

    # q column: set to 0 — at late times (300-400d) the ejecta is optically
    # thin, so any initial thermal energy radiates away instantly.  The
    # sustained luminosity comes from Co56 decay heating, which ARTIS
    # handles via the X_Co56 abundance in model.txt.
    sim_duration_s = t_max_d * DAY_S
    total_energy = 0.0
    q_hot = 0.0

    q = np.zeros(n)

    # Diagnostics
    n_hot = np.sum(is_metal)
    n_outer = np.sum(~is_metal)
    n_source = np.sum(q > 0)
    r_source_min = r_io[is_metal].min() if n_hot > 0 else 0
    r_source_max = r_io[is_metal].max() if n_hot > 0 else 0

    print(f"\nTemperature and energy profile:")
    print(f"  SN_B source (T={t_bb:.0f} K): {n_hot} zones (inner metal-rich)")
    print(f"    r range: {r_source_min:.3e} - {r_source_max:.3e} cm")
    print(f"    = {r_source_min/1e15:.2f} - {r_source_max/1e15:.2f} x 1e15 cm")
    print(f"  Outer ejecta (T>{t_outer_floor:.0f} K): {n_outer} zones (H/He)")
    print(f"  T_gas range: {t_gas.min():.0f} - {t_gas.max():.0f} K")
    mass_hot = mass[is_metal].sum()
    print(f"  Mass of hot zones: {mass_hot/M_SUN:.4f} Msun")
    print(f"  q: 0 (relying on Co56 decay heating)")
    print(f"  Target L_bol: {target_l_bol:.3e} erg/s (from Co56 decay)")
    print(f"  Sim duration: {t_max_d:.1f} d = {sim_duration_s:.3e} s")

    return t_gas, q, target_l_bol


# ---------------------------------------------------------------------------
# Step 5: Write ARTIS model.txt
# ---------------------------------------------------------------------------
def write_model_txt(path, r_io, v_io, rho_io, t_gas, q, t_model_d,
                    comp_io=None, ni56_mass_msun=0.06, total_mass_msun=1.72):
    """Write ARTIS model.txt file (1D spherical format).

    Format (1D spherical, inner-to-outer):
      Line 1: ncells
      Line 2: t_model (days)
      Then per cell:
      #inputcellid velocity_outer logrho X_Fegroup X_Ni56 X_Co56 X_Fe52 X_Cr48 q

    Ni56 and Co56 mass fractions are computed from radioactive decay:
      Ni56 -> Co56 -> Fe56
      tau_Ni56 = 8.77 d (half-life 6.075 d)
      tau_Co56 = 111.5 d (half-life 77.27 d)

    The Ni56 is placed in the inner metal-rich zones (Fe-rich core).
    X_Fegroup is set to the decayed Fe56 fraction.
    """
    n = len(r_io)

    # Decay constants
    tau_ni = 8.77    # days
    tau_co = 111.5   # days

    # Remaining Ni56 and Co56 fractions at t_model
    # (assuming all Ni56 was produced at t=0, i.e. at SN_A explosion)
    frac_ni = np.exp(-t_model_d / tau_ni)
    frac_co_from_ni = (1.0 - frac_ni)  # fraction of initial Ni that decayed to Co
    frac_co_remaining = np.exp(-t_model_d / tau_co)
    # Co56 mass = initial_Ni * (1 - exp(-t/tau_ni)) * exp(-t/tau_co)
    # Ni56 mass = initial_Ni * exp(-t/tau_ni)
    # Fe56 mass (from decay) = initial_Ni - Ni56 - Co56
    ni56_msun = ni56_mass_msun * frac_ni
    co56_msun = ni56_mass_msun * frac_co_from_ni * frac_co_remaining
    fe56_msun = ni56_mass_msun - ni56_msun - co56_msun

    # Mass fractions (relative to total ejecta mass)
    x_ni56_global = ni56_msun / total_mass_msun
    x_co56_global = co56_msun / total_mass_msun
    x_fegrp_global = fe56_msun / total_mass_msun

    print(f"\nRadioactive decay at t={t_model_d:.0f}d (initial Ni56={ni56_mass_msun} Msun):")
    print(f"  Ni56: {ni56_msun:.6e} Msun (X={x_ni56_global:.6e})")
    print(f"  Co56: {co56_msun:.6e} Msun (X={x_co56_global:.6e})")
    print(f"  Fe56 (decay product): {fe56_msun:.6e} Msun (X={x_fegrp_global:.6e})")

    # Compute per-cell mass to distribute Ni56/Co56 in the Fe-rich core
    r_edges = np.concatenate([[0], r_io])
    vol = 4.0/3.0 * PI * (r_edges[1:]**3 - r_edges[:-1]**3)
    mass = vol * rho_io
    total_mass_g = mass.sum()

    # Fe abundance from composition (index of Fe in STELLA_ELEMENTS)
    fe_idx = STELLA_ELEMENTS.index('Fe')
    # Also include Ni (remapped to Fe in ARTIS)
    ni_idx = STELLA_ELEMENTS.index('Ni')

    with open(path, 'w') as f:
        f.write(f"{n}\n")
        f.write(f"{t_model_d:.6f}\n")
        f.write("#inputcellid velocity_outer logrho X_Fegroup X_Ni56 X_Co56 X_Fe52 X_Cr48 q\n")

        for i in range(n):
            cellid = i + 1
            v_outer = v_io[i] / 1e5  # km/s
            log_rho = np.log10(rho_io[i]) if rho_io[i] > 0 else -99.0

            # Place Ni56/Co56 proportional to Fe abundance in each cell
            # (Ni56 was produced in the Fe-core during the explosion)
            fe_frac = comp_io[i, fe_idx] if comp_io is not None else 0.0
            if fe_frac > 0.01 and mass[i] > 0:
                # This cell has Fe-rich material -> place Ni56/Co56 here
                # Mass fractions scaled by Fe abundance (Ni56 was produced in Fe core)
                scale = fe_frac / 0.5  # normalize to peak Fe abundance ~0.5
                x_ni56 = max(x_ni56_global * scale, 1e-10)
                x_co56 = x_co56_global * scale
                x_fegrp = x_fegrp_global * scale
                # Clamp to physical range (don't exceed Fe abundance)
                x_ni56 = min(x_ni56, fe_frac * 0.5)
                x_co56 = min(x_co56, fe_frac * 0.5)
                x_fegrp = min(x_fegrp, fe_frac)
            else:
                x_ni56 = 1e-10
                x_co56 = 0.0
                x_fegrp = 0.0

            x_fe52 = 0.0
            x_cr48 = 0.0
            q_val = q[i]
            f.write(f"   {cellid} {v_outer:.6e} {log_rho:.6f} "
                    f"{x_fegrp:.6e} {x_ni56:.6e} {x_co56:.6e} "
                    f"{x_fe52:.6e} {x_cr48:.6e} {q_val:.6e}\n")

    print(f"\nWrote {path}: {n} cells, t_model={t_model_d:.1f} d")


# ---------------------------------------------------------------------------
# Step 6: Write ARTIS abundances.txt
# ---------------------------------------------------------------------------
def write_abundances_txt(path, comp_io, elems):
    """Write ARTIS abundances.txt file.

    Format: one line per cell, cellid followed by mass fractions for Z=1..30.
    Only ARTIS-supported elements are non-zero; others are remapped.
    """
    n = len(comp_io)
    max_z = 30

    with open(path, 'w') as f:
        for i in range(n):
            cellid = i + 1
            abund = np.zeros(max_z)
            for elem in elems:
                z = ELEMENT_Z.get(elem, 0)
                if z > 0 and z <= max_z and elem in ARTIS_ELEMENTS:
                    elem_idx = STELLA_ELEMENTS.index(elem)
                    abund[z - 1] = comp_io[i, elem_idx]

            total = abund.sum()
            if total > 0:
                abund = abund / total

            vals = ' '.join(f'{a:.8e}' for a in abund)
            f.write(f"{cellid} {vals}\n")

    print(f"Wrote {path}: {n} cells, Z=1..{max_z}")


# ---------------------------------------------------------------------------
# Step 7: Write ARTIS input.txt
# ---------------------------------------------------------------------------
def write_input_txt(path, t_model_d, t_max_d, n_timesteps, n_lte_ts, mpkts):
    """Write ARTIS input.txt file."""
    tmin = t_model_d
    tmax = t_model_d + t_max_d

    lines = [
        "1281360349               # pre_zseed",
        f"{n_timesteps}                       # ntimesteps",
        f"000 {n_timesteps:03d}                  # timestep_start timestep_finish",
        f"{tmin:.6f} {tmax:.6f}        # tmin_days tmax_days",
        "1.33  1.330000001        # nusyn_min_mev nusyn_max_mev (UNUSED)",
        "80                       # nsyn_time (UNUSED)",
        "3. 0.037                 # start and end times for synthesis (UNUSED)",
        "1                        # model_type (UNUSED, autodetected)",
        "4                        # compute r-light curve (UNUSED)",
        "1                        # n_out_it (UNUSED)",
        "1.0                      # CLIGHT/CLIGHT (UNUSED)",
        "-1                       # gamma_kappagrey (UNUSED)",
        "0 0 1                    # syn_dir (UNUSED)",
        "4                        # opacity_case (UNUSED)",
        "1.0e-10                  # rho_crit_para (UNUSED)",
        "-1                       # debug_packet (UNUSED)",
        "0                        # simulation_continued_from_saved",
        "1e-6                     # rfcut (UNUSED)",
        f"{n_lte_ts}                        # num_lte_timesteps",
        "8.0 999                  # optical_depth_is_thick num_grey_timesteps",
        "-1                       # max_bf_continua (UNUSED)",
        f"{mpkts}                  # nprocs_exspec",
        "1                        # do_emission_res (UNUSED)",
        "0.001 1000               # kpktdiffusion (UNUSED)",
    ]

    with open(path, 'w') as f:
        f.write('\n'.join(lines) + '\n')

    print(f"Wrote {path}: tmin={tmin:.1f}d, tmax={tmax:.1f}d, "
          f"{n_timesteps} timesteps, {n_lte_ts} LTE")


# ---------------------------------------------------------------------------
# Step 8: Write vpkt.txt
# ---------------------------------------------------------------------------
def write_vpkt_txt(path):
    """Write ARTIS vpkt.txt (virtual packets disabled)."""
    lines = ["0", "0 10 30", "0 100", "10", "0", "1", "11.5 21.5"]
    with open(path, 'w') as f:
        f.write('\n'.join(lines) + '\n')
    print(f"Wrote {path} (virtual packets disabled)")


# ---------------------------------------------------------------------------
# Step 9: Write SLURM script
# ---------------------------------------------------------------------------
def write_slurm_script(path, work_dir, n_rank, mpkts, t_model_d, t_max_d,
                       artis_exe=ARTIS_EXE):
    """Write a SLURM batch script for the ARTIS run."""
    job_name = f"artis_sn_b_{t_model_d:.0f}d"

    script = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --partition=small
#SBATCH --qos=n0001
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks={n_rank}
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output={work_dir}/slurm_%j.out
#SBATCH --error={work_dir}/slurm_%j.err

set -e

# Load GCC 15 and OpenMPI 4.1 for ARTIS
module purge
module load gcc/15
module load openmpi/4.1
export OMPI_CXX=g++
export LD_LIBRARY_PATH=/mpcdf/soft/RHEL_9/packages/x86_64/gcc/15.1.0/lib64:/mpcdf/soft/RHEL_9/packages/znver4/gsl/gcc_15-15.1.0/2.7.1/lib:${{LD_LIBRARY_PATH}}
export LD_PRELOAD=/mpcdf/soft/RHEL_9/packages/x86_64/gcc/15.1.0/lib64/libstdc++.so.6

# Run ARTIS
cd {work_dir}
echo "Starting ARTIS at $(date)"
echo "Executable: {artis_exe}"
echo "MPI ranks: {n_rank}"
echo "Packets/rank: {mpkts}"
echo "t_model: {t_model_d:.1f} d, t_max: {t_max_d:.1f} d"

mpirun --bind-to none -np {n_rank} {artis_exe}

echo "ARTIS finished at $(date) with exit code $?"

# Check for spec.out
if [ -f spec.out ]; then
    echo "spec.out exists — spectrum produced"
else
    echo "ERROR: spec.out not found — ARTIS may have failed"
fi
"""
    with open(path, 'w') as f:
        f.write(script)
    print(f"Wrote {path}")


# ---------------------------------------------------------------------------
# Step 10: Copy atomic data and run ARTIS
# ---------------------------------------------------------------------------
def copy_atomic_data(work_dir):
    """Copy atomic data files to the working directory."""
    work_dir = Path(work_dir)

    atomic_files = ['adata.txt', 'transitiondata.txt', 'compositiondata.txt',
                    'phixsdata_v2.txt']
    for af in atomic_files:
        src = Path(ARTIS_ATOMIC) / af
        dst = work_dir / af
        if src.exists():
            shutil.copy2(src, dst)
            print(f"  Copied {af} ({src.stat().st_size / 1e6:.1f} MB)")

    data_src = Path(ARTIS_REPO) / "data"
    data_dst = work_dir / "data"
    if data_src.exists() and not data_dst.exists():
        shutil.copytree(data_src, data_dst, dirs_exist_ok=True)
        print(f"  Copied data/ directory")

    for df in ['betaminusdecays.txt', 'alphadecays.txt', 'fissiondecays.txt']:
        src = Path(ARTIS_REPO) / "data" / df
        dst = work_dir / df
        if src.exists() and not dst.exists():
            shutil.copy2(src, dst)


def run_artis(work_dir, n_rank):
    """Run ARTIS with mpirun."""
    work_dir = Path(work_dir)

    env = os.environ.copy()
    env['OMPI_CXX'] = 'g++'
    env['LD_LIBRARY_PATH'] = (
        '/mpcdf/soft/RHEL_9/packages/znver4/gsl/gcc_15-15.1.0/2.7.1/lib:'
        + env.get('LD_LIBRARY_PATH', '')
    )

    cmd = ['mpirun', '--bind-to', 'none', '-np', str(n_rank), ARTIS_EXE]

    print(f"\nRunning ARTIS in {work_dir}")
    print(f"  Executable: {ARTIS_EXE}")
    print(f"  MPI ranks: {n_rank}")
    print(f"  Command: {' '.join(cmd)}")

    result = subprocess.run(cmd, cwd=str(work_dir), env=env,
                           capture_output=True, text=True, timeout=7200)

    if result.returncode != 0:
        print(f"\nARTIS FAILED (exit code {result.returncode})")
        print("STDERR (last 2000 chars):")
        print(result.stderr[-2000:])
        for log_file in sorted(work_dir.glob("output_*.txt")):
            print(f"\n--- {log_file.name} (last 1000 chars) ---")
            content = log_file.read_text(errors='replace')
            print(content[-1000:])
    else:
        print(f"\nARTIS completed successfully")

    return result.returncode


# ---------------------------------------------------------------------------
# Step 11: Plot input structure
# ---------------------------------------------------------------------------
def plot_input_structure(work_dir, r_io, v_io, rho_io, comp_io, t_gas, q,
                          flux_io, t_model_d):
    """Plot the ARTIS input structure for diagnostic purposes."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available, skipping plots")
        return

    work_dir = Path(work_dir)
    n = len(r_io)

    # Convert r to 1e15 cm for plotting
    r_plot = r_io / 1e15

    fig, axes = plt.subplots(3, 2, figsize=(14, 15))

    # Panel 1: Temperature and q
    ax1 = axes[0, 0]
    ax1_twin = ax1.twinx()
    ax1.semilogy(r_plot, t_gas, 'r-', linewidth=1.5, label='T_gas')
    ax1.set_xlabel('r [1e15 cm]')
    ax1.set_ylabel('T_gas [K]', color='r')
    ax1.tick_params(axis='y', labelcolor='r')
    q_plot = np.where(q > 0, q, 1e-99)
    ax1_twin.semilogy(r_plot, q_plot, 'b-', linewidth=1.0, alpha=0.7, label='q')
    ax1_twin.set_ylabel('q [erg/g]', color='b')
    ax1_twin.tick_params(axis='y', labelcolor='b')
    ax1.set_title(f'Temperature and energy (t_model={t_model_d:.1f}d)')
    ax1.axvline(1.0, color='gray', linestyle=':', alpha=0.5, label='r=1e15 cm')
    ax1.legend(loc='upper left', fontsize=8)
    ax1.grid(True, alpha=0.2)

    # Panel 2: Density and velocity
    ax2 = axes[0, 1]
    ax2_twin = ax2.twinx()
    ax2.semilogy(r_plot, rho_io, 'k-', linewidth=1.5, label='rho')
    ax2.set_xlabel('r [1e15 cm]')
    ax2.set_ylabel('rho [g/cm^3]', color='k')
    ax2_twin.plot(r_plot, v_io / 1e5, 'b-', linewidth=1.5, label='v')
    ax2_twin.set_ylabel('v [km/s]', color='b')
    ax2_twin.tick_params(axis='y', labelcolor='b')
    ax2.set_title('Density and velocity')
    ax2.axvline(1.0, color='gray', linestyle=':', alpha=0.5)
    ax2.grid(True, alpha=0.2)

    # Panel 3: Composition (mass fractions)
    ax3 = axes[1, 0]
    artis_elems_sorted = ['H', 'He', 'C', 'O', 'Si', 'Ca', 'Fe']
    colors = {'H': 'blue', 'He': 'green', 'C': 'orange', 'O': 'red',
              'Si': 'purple', 'Ca': 'brown', 'Fe': 'black'}
    for elem in artis_elems_sorted:
        idx = STELLA_ELEMENTS.index(elem)
        ax3.semilogy(r_plot, np.maximum(comp_io[:, idx], 1e-10),
                     linewidth=1.0, label=elem, color=colors[elem])
    ax3.set_xlabel('r [1e15 cm]')
    ax3.set_ylabel('Mass fraction')
    ax3.set_title('Composition (ARTIS elements)')
    ax3.set_ylim(1e-5, 2)
    ax3.axvline(1.0, color='gray', linestyle=':', alpha=0.5)
    ax3.legend(fontsize=8, ncol=2)
    ax3.grid(True, alpha=0.2)

    # Panel 4: STELLA luminosity
    ax4 = axes[1, 1]
    ax4.semilogy(r_plot, np.maximum(np.abs(flux_io), 1e30), 'r-', linewidth=1.5)
    ax4.set_xlabel('r [1e15 cm]')
    ax4.set_ylabel('L [erg/s]')
    ax4.set_title('STELA luminosity profile')
    ax4.axvline(1.0, color='gray', linestyle=':', alpha=0.5)
    ax4.grid(True, alpha=0.2)

    # Panel 5: q (zoomed)
    ax5 = axes[2, 0]
    q_mask = q > 0
    if q_mask.any():
        ax5.semilogy(r_plot[q_mask], q[q_mask], 'b.-', linewidth=1.0, markersize=3)
        ax5.set_xlabel('r [1e15 cm]')
        ax5.set_ylabel('q [erg/g]')
        ax5.set_title('Initial energy (q > 0 zones only)')
        ax5.axvline(1.0, color='gray', linestyle=':', alpha=0.5)
    else:
        ax5.text(0.5, 0.5, 'No q > 0 zones', transform=ax5.transAxes,
                 ha='center', va='center')
    ax5.grid(True, alpha=0.2)

    # Panel 6: Mass per zone
    ax6 = axes[2, 1]
    r_edges = np.concatenate([[0], r_io])
    vol = 4.0/3.0 * PI * (r_edges[1:]**3 - r_edges[:-1]**3)
    mass = vol * rho_io
    ax6.semilogy(r_plot, mass / M_SUN, 'k-', linewidth=1.0)
    ax6.set_xlabel('r [1e15 cm]')
    ax6.set_ylabel('Mass per zone [Msun]')
    ax6.set_title('Zone mass distribution')
    ax6.axvline(1.0, color='gray', linestyle=':', alpha=0.5)
    ax6.grid(True, alpha=0.2)

    fig.suptitle(f'ARTIS SN_B input structure (t_model={t_model_d:.1f}d, '
                 f'{n} cells)', fontsize=14, y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    plot_path = work_dir / "input_structure.png"
    fig.savefig(plot_path, dpi=150, bbox_inches='tight')
    print(f"\nSaved input structure plot to {plot_path}")
    plt.close()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Prepare ARTIS input files from SN_B STELLA .swd")
    parser.add_argument('--swd', type=str, default=DEFAULT_SN_A_SWD,
                        help=f'Path to SN_A .swd file (default: {DEFAULT_SN_A_SWD})')
    parser.add_argument('--target-age', type=float, default=DEFAULT_SN_A_AGE_D,
                        help=f'Homologous evolution target age in days (default: {DEFAULT_SN_A_AGE_D})')
    parser.add_argument('--t-max', type=float, default=DEFAULT_T_MAX_D,
                        help=f'Simulation duration in days (default: {DEFAULT_T_MAX_D})')
    parser.add_argument('--n-timesteps', type=int, default=DEFAULT_N_TIMESTEPS,
                        help=f'Number of timesteps (default: {DEFAULT_N_TIMESTEPS})')
    parser.add_argument('--n-lte', type=int, default=DEFAULT_N_LTE,
                        help=f'Number of LTE timesteps (default: {DEFAULT_N_LTE})')
    parser.add_argument('--max-cells', type=int, default=DEFAULT_MAX_CELLS,
                        help=f'Max grid cells (default: {DEFAULT_MAX_CELLS})')
    parser.add_argument('--work-dir', type=str,
                        default='/viper/ptmp1/sbrennan/SN2024pba/tests/artis_sn_a_400d',
                        help='Working directory')
    parser.add_argument('--nrank', type=int, default=DEFAULT_N_RANK,
                        help=f'MPI ranks (default: {DEFAULT_N_RANK})')
    parser.add_argument('--mpkts', type=int, default=DEFAULT_MPKTS,
                        help=f'Packets per rank (default: {DEFAULT_MPKTS})')
    parser.add_argument('--dry-run', action='store_true',
                        help='Generate input files only, do not run ARTIS')
    parser.add_argument('--run', action='store_true',
                        help='Run ARTIS after generating input files')
    args = parser.parse_args()

    work_dir = Path(args.work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print(f"ARTIS Setup: SN_A ejecta homologously evolved to {args.target_age:.0f}d")
    print("=" * 70)
    print(f"  SN_A .swd: {args.swd}")
    print(f"  Target age: {args.target_age:.0f} d (homologous evolution)")
    print(f"  Sim duration: {args.t_max:.0f} d")
    print(f"  Max cells: {args.max_cells}")
    print(f"  MPI ranks: {args.nrank}")
    print(f"  Packets/rank: {args.mpkts}")
    print(f"  Work dir: {work_dir}")
    print(f"  Dry run: {args.dry_run}")

    # Step 1: Read SN_A profile and homologously evolve to target age
    epoch, t_snap = read_and_evolve_sn_a(args.swd, args.target_age)

    # Step 2: Map composition to ARTIS elements
    comp_artis = map_composition_to_artis(epoch['composition_normalized'])

    # Step 3: Adaptive downsampling
    r_io, rho_io, v_io, comp_io, T_io, Trad_io, flux_io = adaptive_downsample(
        epoch['r'], epoch['rho'], epoch['v_cms'], comp_artis,
        epoch['T'], epoch['T_rad'], epoch['flux'], args.max_cells)

    # Step 4: Set temperature and energy (SN_B illumination of SN_A ejecta)
    t_gas, q, l_stella = set_temperature_and_energy(
        r_io, rho_io, T_io, comp_io, flux_io, t_snap, args.t_max)

    # Step 5: Write ARTIS input files
    print(f"\nWriting ARTIS input files to {work_dir}")
    write_model_txt(work_dir / "model.txt", r_io, v_io, rho_io,
                    t_gas, q, t_snap, comp_io=comp_io)
    write_abundances_txt(work_dir / "abundances.txt", comp_io, STELLA_ELEMENTS)
    write_input_txt(work_dir / "input.txt", t_snap, args.t_max,
                    args.n_timesteps, args.n_lte, args.mpkts)
    write_vpkt_txt(work_dir / "vpkt.txt")

    # Step 6: Plot input structure
    print(f"\nPlotting input structure...")
    plot_input_structure(work_dir, r_io, v_io, rho_io, comp_io, t_gas, q,
                         flux_io, t_snap)

    # Step 7: Copy atomic data
    print(f"\nCopying atomic data...")
    copy_atomic_data(work_dir)

    # Step 8: Write SLURM script
    write_slurm_script(work_dir / "run_artis.slurm", work_dir,
                       args.nrank, args.mpkts, t_snap, args.t_max)

    # Step 9: Run ARTIS (if requested)
    if args.run and not args.dry_run:
        print(f"\nRunning ARTIS...")
        retcode = run_artis(work_dir, args.nrank)
        if retcode != 0:
            print(f"\nARTIS failed with exit code {retcode}")
            return retcode
    elif args.dry_run:
        print(f"\nDry run: input files generated, ARTIS not executed")
        print(f"  To run: cd {work_dir} && sbatch run_artis.slurm")
        print(f"  Or:     cd {work_dir} && mpirun -np {args.nrank} {ARTIS_EXE}")

    # Summary
    print(f"\n{'=' * 70}")
    print(f"Setup complete:")
    print(f"  t_model: {t_snap:.4f} d")
    print(f"  Cells: {len(r_io)}")
    print(f"  R range: {r_io.min():.3e} - {r_io.max():.3e} cm")
    print(f"  V range: {v_io.min()/1e5:.1f} - {v_io.max()/1e5:.1f} km/s")
    print(f"  T range: {t_gas.min():.0f} - {t_gas.max():.0f} K")
    print(f"  STELLA L: {l_stella:.3e} erg/s")
    print(f"  Source zones: {np.sum(q > 0)} of {len(q)}")
    if q.max() > 0:
        r_src = r_io[q > 0]
        print(f"  Source r: {r_src.min():.3e} - {r_src.max():.3e} cm")
        print(f"           = {r_src.min()/1e15:.2f} - {r_src.max()/1e15:.2f} x 1e15 cm")
    print(f"  Work dir: {work_dir}")
    print(f"{'=' * 70}")

    return 0


if __name__ == '__main__':
    sys.exit(main())
