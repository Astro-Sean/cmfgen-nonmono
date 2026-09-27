#!/usr/bin/env python3
"""
Convert a Kepler pre-supernova stellar model to STELLA .hyd + .abn files.

The Kepler presn file format (VERSION 10111) has columns:
  zone, m_outer(g), r_outer(cm), v(cm/s), rho(g/cc), T(K), P(erg/cc),
  eps, S, omega, A_bar, Ye, stability, network,
  [isotopes: n, H1, He3, He4, C12, N14, O16, Ne20, Mg24, Si28, S32, Ar36,
             Ca40, Ti44, Cr48, Fe52, Fe54, Ni56, Fe56, Fe]

This script:
  1. Reads the Kepler presn profile
  2. Excises the inner compact remnant (mass < mass_cut_msun)
  3. Remaps the remaining ejecta onto nzon STELLA zones (log-spaced in radius)
  4. Zeroes the pre-collapse infall velocities (SN_A convention: the model
     is blown up by the EBurst thermal bomb in the .dat — see below)
  5. Sets an initial temperature floor (STELLA radiation floor)
  6. Maps composition to STELLA 7-element scheme: H, He, C, O, Si, Ca, FePeak
  7. Writes .hyd, .abn, .xni, .dat, and input_profile files

Explosion convention (two-event workflow):
  This script builds the SN_A setup — the KEPLER pre-SN model is blown up
  by the EBurst thermal bomb (--explosion-energy > 0).  Keep EKO=0:
  EKO>0 makes STELLA overwrite the entire .hyd velocity column with a
  triangular-in-mass profile, destroying the constructed wind velocity.
  The --initial-velocity-energy option is a bomb-free escape hatch that
  writes the same triangular v(m) directly into the .hyd (wind zones keep
  their beta-law); never combine it with EBurst.

  The SN_B setup is the opposite convention — an embedded shell driven by
  its own mass and velocity (EBurst=0, EKO=0): see prepare_SN_B_model.py.

Usage:
    python convert_kepler_presn_to_stella.py <presn_file> <model_name> [options]

Example:
    python convert_kepler_presn_to_stella.py he8.00-presn SN_A_IIb \\
        --mass-cut 1.5 --nzon 500 --outdir /ptmp/sbrennan/SN2024pba/SN_A/
"""

import argparse
import sys
import numpy as np
from pathlib import Path
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter1d


# ── Physical constants ────────────────────────────────────────────────────────
MSUN   = 1.989e33   # g
K_B    = 1.38065e-16  # erg/K
M_H    = 1.6736e-24   # g
A_RAD  = 7.566e-15    # erg cm^-3 K^-4
SIGMA  = 5.6704e-5    # erg cm^-2 s^-1 K^-4


# ── Kepler file reader ────────────────────────────────────────────────────────
def read_kepler_presn(fpath):
    """Read a Kepler VERSION 10111 presn file.

    Returns a dict of numpy arrays with keys:
        zone, m_g, r_cm, v_cms, rho, T_K, P_cgs,
        He4, C12, N14, O16, Ne20, Mg24, Si28, S32, Ar36, Ca40,
        Ti44, Cr48, Fe52, Fe54, Ni56, Fe56, Fe_tot
    """
    records = []
    with open(fpath) as f:
        for line in f:
            cols = line.split()
            if not cols or not cols[0].rstrip(':').lstrip('-').isdigit():
                continue
            if len(cols) < 14:
                continue
            try:
                zone = int(cols[0].rstrip(':'))
            except ValueError:
                continue

            def _f(s):
                try:
                    return float(s)
                except ValueError:
                    return 0.0

            abund = cols[14:]   # isotope columns start at index 14
            records.append({
                'zone': zone,
                'm_g':  float(cols[1]),
                'r_cm': float(cols[2]),
                'v_cms': float(cols[3]),
                'rho':  float(cols[4]),
                'T_K':  float(cols[5]),
                'P_cgs': float(cols[6]),
                # isotopes (0-based index in abund list)
                'He4':  _f(abund[3])  if len(abund) > 3  else 0.0,
                'C12':  _f(abund[4])  if len(abund) > 4  else 0.0,
                'N14':  _f(abund[5])  if len(abund) > 5  else 0.0,
                'O16':  _f(abund[6])  if len(abund) > 6  else 0.0,
                'Ne20': _f(abund[7])  if len(abund) > 7  else 0.0,
                'Mg24': _f(abund[8])  if len(abund) > 8  else 0.0,
                'Si28': _f(abund[9])  if len(abund) > 9  else 0.0,
                'S32':  _f(abund[10]) if len(abund) > 10 else 0.0,
                'Ar36': _f(abund[11]) if len(abund) > 11 else 0.0,
                'Ca40': _f(abund[12]) if len(abund) > 12 else 0.0,
                'Ti44': _f(abund[13]) if len(abund) > 13 else 0.0,
                'Ni56': _f(abund[17]) if len(abund) > 17 else 0.0,
                'Fe56': _f(abund[18]) if len(abund) > 18 else 0.0,
                'Fe_tot': _f(abund[19]) if len(abund) > 19 else 0.0,
            })

    # Convert to numpy arrays
    out = {k: np.array([r[k] for r in records]) for k in records[0].keys()}
    print(f"Read {len(records)} zones from {fpath}")
    print(f"  Total mass:   {out['m_g'][-1]/MSUN:.4f} Msun")
    print(f"  Outer radius: {out['r_cm'][-1]:.3e} cm  ({out['r_cm'][-1]/6.96e10:.2f} Rsun)")
    print(f"  T range:      {out['T_K'].min():.2e} – {out['T_K'].max():.2e} K")
    print(f"  rho range:    {out['rho'].min():.2e} – {out['rho'].max():.2e} g/cc")
    return out


# ── Composition remapper ──────────────────────────────────────────────────────
def build_stella_composition(kep):
    """Map Kepler isotopes to STELLA 7-element (H, He, C, O, Si, Ca, FePeak) scheme.

    STELLA .abn columns (Natom=15):
        H, He, C, N, O, Ne, Na, Mg, Al, Si, S, Ar, Ca, Fe, Ni

    Mapping:
        H  = 0 (He star has no H, unless H surface layer added later)
        He = He4
        C  = C12 + N14
        O  = O16
        Si = Si28 + Ne20 + Mg24 (Si-group)
        Ca = Ca40
        Fe = Fe56 + Fe_tot + S32 + Ar36 + Ti44 (Fe-peak, without Ca)

    Ca is written as its own column for CMFGEN spectral modelling (Ca II H&K,
    [Ca II] NIR). STELLA's opacity tables handle Ca separately via the .abn
    Ca column.

    Ni56 is NOT placed in the .abn XNifor column (set to 0). eve2 reads
    XNifor and adds it to the Fe abundance for opacity; placing 0.251 Msun
    of Ni56 there increases the inner-zone Fe opacity by ~2.5x, which also
    causes solver stalls. Ni56 is handled exclusively via the .xni file for
    radioactive decay heating.

    Normalised to sum = 1 per zone.
    """
    H_frac   = np.zeros(len(kep['rho']))                 # no H in He star
    He_frac  = kep['He4'].copy()
    C_frac   = kep['C12'] + kep['N14']
    O_frac   = kep['O16']
    Si_frac  = kep['Si28'] + kep['Ne20'] + kep['Mg24']
    Ca_frac  = kep['Ca40'].copy()
    Fe_frac  = (kep['Fe56'] + kep['Fe_tot'] +
                kep['S32'] + kep['Ar36'] + kep['Ti44'])
    Ni56_frac = kep['Ni56'].copy()

    total = H_frac + He_frac + C_frac + O_frac + Si_frac + Ca_frac + Fe_frac
    # Normalise; guard against zero-composition zones (innermost NSE)
    bad = total < 1e-10
    total = np.where(bad, 1.0, total)
    # In the innermost NSE-dominated zones all Fe; set FePeak=1 there
    Fe_frac = np.where(bad, 1.0, Fe_frac)

    He_frac  /= total; C_frac  /= total; O_frac /= total
    Si_frac  /= total; Ca_frac /= total; Fe_frac /= total

    return H_frac, He_frac, C_frac, O_frac, Si_frac, Ca_frac, Fe_frac, Ni56_frac


def add_h_surface_layer(g, h_fraction, h_mass_msun, n_transition_zones=10):
    """Add a thin surface layer of H to the outermost zones.

    Modifies the composition in-place: sets fH = h_fraction in the outermost
    zones covering h_mass_msun, and rescales the other mass fractions by
    (1 - h_fraction) to maintain sum = 1.

    The H/He composition transition is smoothed over n_transition_zones using
    a linear ramp to prevent a sharp opacity jump that can trigger
    radiation transport instabilities in STELLA (the solver produces negative
    luminosities and the photosphere collapses to Ncnd=1).

    Parameters:
    -----------
    g : dict
        Remapped grid dict with keys 'r', 'rho', 'fH', 'fHe', 'fC', 'fO', 'fSi', 'fFe'
    h_fraction : float
        H mass fraction in the surface layer (e.g. 0.1 = 10% H by mass)
    h_mass_msun : float
        Mass of the H surface layer in Msun (controls how many outer zones get H)
    n_transition_zones : int
        Number of zones over which to smooth the H/He transition (default: 10)
    """
    if h_fraction <= 0 or h_mass_msun <= 0:
        return

    r = g['r']; rho = g['rho']
    r_lo = np.concatenate([[g['r_cut']], r[:-1]])
    dm = np.maximum((4/3)*np.pi*(r**3 - r_lo**3)*rho, 0.0)
    m_cum = np.cumsum(dm) / MSUN   # cumulative ejecta mass in Msun
    m_total = m_cum[-1]

    # Find outermost zones covering h_mass_msun
    # m_cum increases from inner to outer; we want the outer zones where
    # (m_total - m_cum) < h_mass_msun
    h_mask = (m_total - m_cum) < h_mass_msun
    n_h_zones = np.sum(h_mask)
    if n_h_zones == 0:
        # Layer too thin for even one zone; apply to outermost zone
        h_mask[-1] = True
        n_h_zones = 1

    h_mass_actual = np.sum(dm[h_mask]) / MSUN
    print(f"\n  H surface layer: {n_h_zones} outermost zones, "
          f"H fraction = {h_fraction:.4f}, "
          f"layer mass = {h_mass_actual:.4f} Msun (requested {h_mass_msun:.4f})")

    # Smooth the H/He transition over n_transition_zones to prevent
    # a sharp opacity jump that triggers radiation transport instabilities.
    # Instead of a step function (fH=0 below, fH=h_fraction above), ramp
    # linearly from 0 to h_fraction over n_transition_zones at the inner
    # edge of the H layer.
    h_indices = np.where(h_mask)[0]
    if n_h_zones > n_transition_zones and n_transition_zones > 0:
        transition_start = h_indices[0] - n_transition_zones
        transition_end = h_indices[0]
        if transition_start < 0:
            transition_start = 0

        # Build a smooth H fraction profile
        h_frac_profile = np.zeros(len(r))
        h_frac_profile[h_mask] = h_fraction

        # Linear ramp over the transition zones
        trans_idx = np.arange(transition_start, transition_end)
        if len(trans_idx) > 0:
            x = (trans_idx - transition_start) / max(transition_end - transition_start - 1, 1)
            h_frac_profile[trans_idx] = h_fraction * x

        print(f"  H/He transition: smoothed over {n_transition_zones} zones "
              f"(zones {transition_start}-{transition_end-1})")
    else:
        # Not enough H zones for smoothing; use step function
        h_frac_profile = np.zeros(len(r))
        h_frac_profile[h_mask] = h_fraction

    # Apply H fraction and rescale other elements
    # For each zone, the H fraction is h_frac_profile[i], and the remaining
    # elements are scaled by (1 - h_frac_profile[i]) relative to their
    # original values (which have fH=0).
    scale_arr = 1.0 - h_frac_profile  # per-zone scaling factor
    for key in ['fHe', 'fC', 'fO', 'fSi', 'fCa', 'fFe']:
        g[key] *= scale_arr
    g['fH'] = h_frac_profile.copy()

    # Re-normalise all zones (H layer zones and non-H zones)
    comp_sum = g['fH'] + g['fHe'] + g['fC'] + g['fO'] + g['fSi'] + g['fCa'] + g['fFe']
    comp_sum = np.maximum(comp_sum, 1e-10)
    for key in ['fH', 'fHe', 'fC', 'fO', 'fSi', 'fCa', 'fFe']:
        g[key] /= comp_sum


def add_extended_h_envelope(g, h_fraction, h_mass_msun, r_outer_rsun,
                             n_envelope_zones=50, T_surface_K=6000.0,
                             T_base_K=None, n_rho=3.5, n_T=2.0, verbose=True):
    """Add an extended H-rich envelope beyond the compact He core surface.

    This implements the Bersten et al. (2012) approach for SN 2011dh: a
    compact He core (R ~ 1 Rsun) with a low-mass (~0.1 Msun) H-rich envelope
    extended to R ~ 200 Rsun, matching the yellow supergiant progenitor
    detected in pre-explosion images.

    The envelope has a power-law density profile:
        rho(r) = rho_base * (r / r_base)^(-n_rho)
    where rho_base is set to match the He core surface density, and n_rho
    controls the falloff (default 3.5, between a wind n=2 and a steeper
    radiative envelope n~8). The total envelope mass is controlled by
    h_mass_msun — the density normalization is adjusted iteratively to
    give the requested mass on the discrete zone grid.

    Temperature follows a power law from T_base at the He core surface
    to T_surface_K at the outer edge:
        T(r) = T_base * (r / r_base)^(-n_T)
    T_base defaults to 50,000 K (appropriate for the H recombination
    temperature at the base of a YSG envelope), NOT the KEPLER He core
    surface temperature (which is ~5e6 K for a compact stripped core
    and would produce an unphysically hot envelope).

    The envelope composition is set to h_fraction H + (1-h_fraction) He,
    with trace metals rescaled from the He core values.

    Parameters:
    -----------
    g : dict
        Remapped grid dict (modified in-place — envelope zones inserted
        between existing ejecta and any CSM wind)
    h_fraction : float
        H mass fraction in the envelope (e.g. 0.7)
    h_mass_msun : float
        Total mass of the H envelope in Msun
    r_outer_rsun : float
        Outer radius of the H envelope in Rsun (e.g. 200 for SN 2011dh)
    n_envelope_zones : int
        Number of STELLA zones for the envelope (default: 50)
    T_surface_K : float
        Temperature at the outer edge of the envelope (default: 6000 K, YSG)
    T_base_K : float or None
        Temperature at the base of the envelope (default: None -> 50000 K).
        A real YSG H envelope has T_base ~ 30,000-100,000 K at the He/H
        interface (H recombination temperature). Do NOT use the KEPLER He
        core surface T (~5e6 K) — that would produce an unphysically hot
        envelope and excess UV emission.
    n_rho : float
        Density power-law exponent (default: 3.5)
    n_T : float
        Temperature power-law exponent (default: 2.0)
    verbose : bool
        Print diagnostics
    """
    if h_mass_msun <= 0 or r_outer_rsun <= 0 or n_envelope_zones <= 0:
        return

    r_base = g['r'][-1]  # He core surface radius
    r_outer = r_outer_rsun * 6.96e10  # convert Rsun to cm

    # Use physically appropriate T_base for YSG envelope (NOT KEPLER He core T)
    if T_base_K is None:
        T_base_K = 50000.0  # H recombination temperature at He/H interface
    T_base = T_base_K

    if r_outer <= r_base:
        if verbose:
            print(f"\n  [Extended H envelope] r_outer ({r_outer:.3e}) <= r_base "
                  f"({r_base:.3e}), skipping")
        return

    # Generate log-spaced envelope zones from r_base to r_outer
    r_env = np.logspace(np.log10(r_base), np.log10(r_outer), n_envelope_zones + 1)
    # Use zone outer radii (drop the first which equals r_base)
    r_env = r_env[1:]  # n_envelope_zones points

    # Density profile: rho(r) = rho_norm * (r/r_base)^(-n_rho)
    # Iteratively adjust rho_norm to match the requested mass on the
    # discrete zone grid (analytical integral overestimates by ~15%
    # due to log-spacing).
    target_mass_g = h_mass_msun * MSUN

    # Initial guess from analytical integral
    if abs(n_rho - 3.0) > 0.01:
        integral = (r_outer**(3.0 - n_rho) - r_base**(3.0 - n_rho)) / (3.0 - n_rho)
        rho_norm = target_mass_g / (4.0 * np.pi * r_base**n_rho * integral)
    else:
        integral = np.log(r_outer / r_base)
        rho_norm = target_mass_g / (4.0 * np.pi * r_base**3.0 * integral)

    # Refine with discrete zone mass calculation (2-3 iterations)
    r_lo_env = np.concatenate([[r_base], r_env[:-1]])
    for _ in range(5):
        rho_env = rho_norm * (r_env / r_base)**(-n_rho)
        dm_env = (4.0/3.0) * np.pi * (r_env**3 - r_lo_env**3) * rho_env
        m_discrete = np.sum(dm_env)
        if m_discrete > 0:
            rho_norm *= target_mass_g / m_discrete

    rho_env = rho_norm * (r_env / r_base)**(-n_rho)

    # Temperature profile: T(r) = T_base * (r/r_base)^(-n_T)
    # Clamped to T_surface_K at the outer edge
    T_env = T_base * (r_env / r_base)**(-n_T)
    T_env = np.maximum(T_env, T_surface_K)
    # Ensure outer edge is exactly T_surface_K
    T_env[-1] = T_surface_K

    # Velocity: zero (pre-explosion, static envelope)
    v_env = np.zeros(n_envelope_zones)

    # Composition: H-rich envelope
    fH_env = np.full(n_envelope_zones, h_fraction)
    fHe_env = np.full(n_envelope_zones, 1.0 - h_fraction)
    # Trace metals from He core surface (rescaled)
    fC_surf = g['fC'][-1]; fO_surf = g['fO'][-1]; fSi_surf = g['fSi'][-1]
    fCa_surf = g['fCa'][-1]; fFe_surf = g['fFe'][-1]
    metals_frac = fC_surf + fO_surf + fSi_surf + fCa_surf + fFe_surf
    if metals_frac > 0:
        # Scale metals to fill remaining mass fraction
        remaining = 1.0 - h_fraction - (1.0 - h_fraction)  # = 0, so use small trace
        # Actually, set H and He, then add trace metals by rescaling He
        trace = min(metals_frac * 0.1, 0.01)  # small trace
        fHe_env = np.full(n_envelope_zones, 1.0 - h_fraction - trace)
        fC_env = np.full(n_envelope_zones, fC_surf / metals_frac * trace)
        fO_env = np.full(n_envelope_zones, fO_surf / metals_frac * trace)
        fSi_env = np.full(n_envelope_zones, fSi_surf / metals_frac * trace)
        fCa_env = np.full(n_envelope_zones, fCa_surf / metals_frac * trace)
        fFe_env = np.full(n_envelope_zones, fFe_surf / metals_frac * trace)
    else:
        fC_env = np.zeros(n_envelope_zones)
        fO_env = np.zeros(n_envelope_zones)
        fSi_env = np.zeros(n_envelope_zones)
        fCa_env = np.zeros(n_envelope_zones)
        fFe_env = np.zeros(n_envelope_zones)
    fNi56_env = np.zeros(n_envelope_zones)

    # Compute actual envelope mass
    r_lo_env = np.concatenate([[r_base], r_env[:-1]])
    dm_env = (4.0/3.0) * np.pi * (r_env**3 - r_lo_env**3) * rho_env
    m_env_actual = np.sum(dm_env) / MSUN

    if verbose:
        print(f"\n  [Extended H envelope] {n_envelope_zones} zones, "
              f"R: {r_base/6.96e10:.2f} -> {r_outer/6.96e10:.1f} Rsun")
        print(f"    H fraction = {h_fraction:.3f}, "
              f"envelope mass = {m_env_actual:.4f} Msun (requested {h_mass_msun:.4f})")
        print(f"    rho at base: {rho_env[0]:.3e} g/cc, "
              f"rho at outer: {rho_env[-1]:.3e} g/cc")
        print(f"    T at base: {T_env[0]:.3e} K, T at outer: {T_env[-1]:.3e} K")
        print(f"    n_rho = {n_rho}, n_T = {n_T}")

    # Append envelope zones to the grid
    g['r'] = np.concatenate([g['r'], r_env])
    g['rho'] = np.concatenate([g['rho'], rho_env])
    g['T'] = np.concatenate([g['T'], T_env])
    g['v'] = np.concatenate([g['v'], v_env])
    g['fH'] = np.concatenate([g['fH'], fH_env])
    g['fHe'] = np.concatenate([g['fHe'], fHe_env])
    g['fC'] = np.concatenate([g['fC'], fC_env])
    g['fO'] = np.concatenate([g['fO'], fO_env])
    g['fSi'] = np.concatenate([g['fSi'], fSi_env])
    g['fCa'] = np.concatenate([g['fCa'], fCa_env])
    g['fFe'] = np.concatenate([g['fFe'], fFe_env])
    g['fNi56'] = np.concatenate([g['fNi56'], fNi56_env])


def smooth_steep_gradients(g, max_log_gradient=5.0, sigma_zones=3.0, verbose=True):
    """Smooth steepest density/temperature gradients to help STELLA solver converge.

    STELLA's GEAR+BGH integrator can stall (BADSTE / ifail 33 / KFLAG=-2) when
    the density profile has very steep gradients (|dlog rho/dlog r| > ~6) in the
    bomb-heated region.  This function applies a mild Gaussian smoothing to
    log(rho) and log(T) in log(r) space, but ONLY in regions where the gradient
    exceeds max_log_gradient.  The overall profile shape is preserved; only the
    sharpest features are rounded.

    Parameters:
    -----------
    g : dict
        Remapped grid dict with keys 'r', 'rho', 'T'
    max_log_gradient : float
        Threshold |dlog rho/dlog r| above which smoothing is applied (default: 5.0)
    sigma_zones : float
        Gaussian smoothing width in zone units (default: 3.0 = ~3 zones FWHM)
    verbose : bool
        Print smoothing diagnostics

    Returns:
    --------
    g : dict
        Modified in-place with smoothed rho and T
    """
    r = g['r']; rho = g['rho']; T = g['T']
    lr = np.log10(r)
    lrho = np.log10(np.maximum(rho, 1e-300))
    lT = np.log10(np.maximum(T, 1e-300))

    # Compute gradients
    dlrho = np.gradient(lrho, lr)
    dlT = np.gradient(lT, lr)

    # Find zones where either density or temperature gradient is too steep
    steep_mask_rho = np.abs(dlrho) > max_log_gradient
    steep_mask_T = np.abs(dlT) > max_log_gradient
    steep_mask = steep_mask_rho | steep_mask_T
    n_steep = np.sum(steep_mask)

    if n_steep == 0:
        if verbose:
            print(f"\n  [Smooth] No zones exceed |dlog/dlog r| > {max_log_gradient:.1f}; no smoothing needed")
        return g

    if verbose:
        worst_rho = np.argmax(np.abs(dlrho))
        worst_T = np.argmax(np.abs(dlT))
        print(f"\n  [Smooth] {n_steep} zones have |dlog/dlog r| > {max_log_gradient:.1f}")
        print(f"    Worst rho gradient: {dlrho[worst_rho]:.2f} at zone {worst_rho} (r={r[worst_rho]:.3e} cm)")
        print(f"    Worst T gradient: {dlT[worst_T]:.2f} at zone {worst_T} (r={r[worst_T]:.3e} cm)")
        print(f"    Applying Gaussian smoothing (sigma={sigma_zones:.1f} zones)")

    # Apply Gaussian smoothing to the full profile in log-log space
    # The Gaussian with sigma_zones will smooth the steep region while
    # leaving the smooth regions essentially unchanged
    lrho_smooth = gaussian_filter1d(lrho, sigma=sigma_zones, mode='nearest')
    lT_smooth = gaussian_filter1d(lT, sigma=sigma_zones, mode='nearest')

    # Only apply smoothing where the gradient is steep (blend factor)
    # Use a smooth transition: full smoothing at 2x threshold, none at threshold
    # Blend based on the MAX of density and temperature gradients
    blend = np.clip((np.maximum(np.abs(dlrho), np.abs(dlT)) - max_log_gradient) / max_log_gradient, 0, 1)
    # Extend blend to neighbors (so the transition is smooth)
    blend = gaussian_filter1d(blend, sigma=sigma_zones*0.5, mode='nearest')

    lrho_new = lrho * (1 - blend) + lrho_smooth * blend
    lT_new = lT * (1 - blend) + lT_smooth * blend

    # Check new gradients
    dlrho_new = np.gradient(lrho_new, lr)
    dlT_new = np.gradient(lT_new, lr)
    worst_rho_new = np.argmax(np.abs(dlrho_new))
    worst_T_new = np.argmax(np.abs(dlT_new))

    if verbose:
        print(f"    After smoothing: worst rho gradient = {dlrho_new[worst_rho_new]:.2f} "
              f"(was {dlrho[worst_rho]:.2f})")
        print(f"    After smoothing: worst T gradient = {dlT_new[worst_T_new]:.2f} "
              f"(was {dlT[worst_T]:.2f})")
        print(f"    Max density change: {np.max(np.abs(10**lrho_new - 10**lrho))/np.max(10**lrho)*100:.2f}%")

    g['rho'] = 10.0 ** lrho_new
    g['T'] = np.maximum(10.0 ** lT_new, np.min(T))  # Don't let T go below original min

    return g


def plot_composition_smoothing(g_before, g_after, model_name, outpath):
    """Generate a 2-panel before/after composition smoothing plot.

    Panel 1 (top):    Enclosed mass (Msun) vs mass fraction
    Panel 2 (bottom): Radius (cm) vs mass fraction

    Each panel shows all species (H, He, C, O, Si, Fe) with the
    before profile as faint solid lines and the after (smoothed)
    profile as bold dashed lines.

    Parameters:
    -----------
    g_before : dict
        Grid dict BEFORE composition smoothing (with keys fH, fHe, fC, fO, fSi, fFe, r, rho)
    g_after : dict
        Grid dict AFTER composition smoothing
    model_name : str
        Model name for the plot title
    outpath : str or Path
        Output PNG path
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    MSUN = 1.989e33
    comp_keys = [('fH', 'H', 'tab:blue'), ('fHe', 'He', 'tab:red'),
                 ('fC', 'C', 'tab:green'), ('fO', 'O', 'tab:orange'),
                 ('fSi', 'Si', 'tab:purple'), ('fFe', 'Fe', 'tab:brown')]

    def enclosed_mass(g):
        r = g['r']; rho = g['rho']
        r_lo = np.concatenate([[r[0] * 0.99], r[:-1]])
        dm = (4.0 / 3.0) * np.pi * (r**3 - r_lo**3) * rho
        m_cum = np.cumsum(dm) / MSUN
        return m_cum  # mass enclosed up to each zone (from centre outward)

    m_b = enclosed_mass(g_before)
    m_a = enclosed_mass(g_after)
    r_b = g_before['r']
    r_a = g_after['r']

    fig, axes = plt.subplots(2, 1, figsize=(12, 10), sharex=False)

    # --- Panel 1: Enclosed mass vs mass fraction ---
    ax = axes[0]
    for key, label, color in comp_keys:
        ax.plot(m_b, g_before[key], color=color, alpha=0.35, linewidth=1.2, linestyle='-')
        ax.plot(m_a, g_after[key], color=color, linewidth=2.0, linestyle='--', label=label)
    ax.set_xlabel('Enclosed mass (M$_\\odot$)', fontsize=12)
    ax.set_ylabel('Mass fraction', fontsize=12)
    ax.set_title(f'{model_name} — Composition smoothing (enclosed mass)', fontsize=13)
    ax.legend(fontsize=9, ncol=6, loc='upper center', framealpha=0.8)
    ax.set_ylim(-0.05, 1.05)
    ax.grid(True, alpha=0.3)
    # Mark before/after in legend
    ax.plot([], [], 'k-', alpha=0.35, linewidth=1.2, label='before')
    ax.plot([], [], 'k--', linewidth=2.0, label='after')
    ax.legend(fontsize=9, ncol=8, loc='upper center', framealpha=0.8)

    # --- Panel 2: Radius vs mass fraction ---
    ax = axes[1]
    for key, label, color in comp_keys:
        ax.semilogx(r_b, g_before[key], color=color, alpha=0.35, linewidth=1.2, linestyle='-')
        ax.semilogx(r_a, g_after[key], color=color, linewidth=2.0, linestyle='--', label=label)
    ax.set_xlabel('Radius (cm)', fontsize=12)
    ax.set_ylabel('Mass fraction', fontsize=12)
    ax.set_title(f'{model_name} — Composition smoothing (radius)', fontsize=13)
    ax.legend(fontsize=9, ncol=6, loc='upper center', framealpha=0.8)
    ax.set_ylim(-0.05, 1.05)
    ax.grid(True, alpha=0.3)
    ax.plot([], [], 'k-', alpha=0.35, linewidth=1.2, label='before')
    ax.plot([], [], 'k--', linewidth=2.0, label='after')
    ax.legend(fontsize=9, ncol=8, loc='upper center', framealpha=0.8)

    plt.tight_layout()
    plt.savefig(str(outpath), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Composition smoothing plot: {outpath}")


def smooth_composition(g, sigma_zones=5.0, max_fraction_jump=0.05, max_dX=0.08, verbose=True):
    """Smooth composition profiles to prevent sharp opacity jumps.

    STELLA's radiation transport solver can develop instabilities when the
    composition changes abruptly across a single zone — the opacity jumps,
    causing the photosphere to oscillate and eventually collapse (negative
    luminosities, Ncnd → 1).  This function applies Gaussian smoothing to
    each mass fraction profile in zone-index space, but ONLY in regions
    where the composition is changing rapidly (|dX/dzone| > max_fraction_jump).

    The smoothing preserves the overall composition structure (shell masses
    are approximately conserved) while rounding off the sharpest interfaces
    over ~2*sigma_zones zones.

    After smoothing, each zone is re-normalised so that the mass fractions
    sum to exactly 1.0.

    Parameters:
    -----------
    g : dict
        Remapped grid dict with composition keys 'fH', 'fHe', 'fC', 'fO',
        'fSi', 'fCa', 'fFe'
    sigma_zones : float
        Gaussian smoothing width in zone units (default: 5.0)
    max_fraction_jump : float
        Threshold |dX/dzone| above which smoothing is applied (default: 0.05)
    max_dX : float
        Maximum per-zone composition change enforced by the derivative-capping
        pass after Gaussian smoothing (default: 0.08). This is what erases
        sharp shell boundaries (e.g. the O-shell) that survive the Gaussian
        pass; raise it to preserve more shell structure.
    verbose : bool
        Print smoothing diagnostics

    Returns:
    --------
    g : dict
        Modified in-place with smoothed composition
    """
    from scipy.ndimage import gaussian_filter1d

    comp_keys = ['fH', 'fHe', 'fC', 'fO', 'fSi', 'fCa', 'fFe']
    n = len(g['r'])

    # Store original composition for comparison
    orig = {k: g[k].copy() for k in comp_keys}

    # Compute per-zone jumps for each species
    max_jump = np.zeros(n)
    for k in comp_keys:
        jumps = np.abs(np.diff(g[k], prepend=g[k][0]))
        max_jump = np.maximum(max_jump, jumps)

    # Find zones where any composition jump exceeds the threshold
    steep_mask = max_jump > max_fraction_jump
    n_steep = np.sum(steep_mask)

    if n_steep == 0:
        if verbose:
            print(f"\n  [CompSmooth] No zones exceed |dX/dzone| > {max_fraction_jump:.3f}; no smoothing needed")
        return g

    if verbose:
        # Find the worst jump for each species
        print(f"\n  [CompSmooth] {n_steep} zones have |dX/dzone| > {max_fraction_jump:.3f}")
        for k in comp_keys:
            jumps = np.abs(np.diff(g[k], prepend=g[k][0]))
            worst = np.argmax(jumps)
            print(f"    {k}: max jump = {jumps[worst]:.4f} at zone {worst} "
                  f"({g[k][worst]:.4f} -> {g[k][worst+1]:.4f})" if worst+1 < n else
                  f"    {k}: max jump = {jumps[worst]:.4f} at zone {worst}")
        print(f"    Applying Gaussian smoothing (sigma={sigma_zones:.1f} zones)")

    # Apply Gaussian smoothing to each composition species
    smoothed = {}
    for k in comp_keys:
        smoothed[k] = gaussian_filter1d(g[k], sigma=sigma_zones, mode='nearest')

    # Blend: only apply smoothing where jumps are steep
    # Extend the blend mask to neighbours so the transition is smooth
    blend = np.clip((max_jump - max_fraction_jump) / max_fraction_jump, 0, 1)
    blend = gaussian_filter1d(blend, sigma=sigma_zones * 0.5, mode='nearest')

    for k in comp_keys:
        g[k] = g[k] * (1 - blend) + smoothed[k] * blend

    # Re-normalise: ensure mass fractions sum to 1 in every zone
    comp_sum = np.zeros(n)
    for k in comp_keys:
        comp_sum += g[k]
    comp_sum = np.maximum(comp_sum, 1e-30)
    for k in comp_keys:
        g[k] /= comp_sum

    # Report results
    if verbose:
        print(f"    After smoothing:")
        for k in comp_keys:
            jumps = np.abs(np.diff(g[k], prepend=g[k][0]))
            worst = np.argmax(jumps)
            worst_val = jumps[worst] if worst < n else 0
            print(f"      {k}: max jump = {worst_val:.4f} (was "
                  f"{np.abs(np.diff(orig[k], prepend=orig[k][0])).max():.4f})")

    # Derivative-capping pass: enforce max per-zone change by iteratively
    # spreading large jumps over neighbouring zones.  This catches the
    # very sharp nuclear burning shell interfaces (e.g. O→Si→C over 2
    # zones) that Gaussian smoothing alone cannot fully round off.
    # NOTE: this is the dominant source of shell-detail loss (e.g. the
    # O-shell) — raise max_dX (function argument) to preserve more structure.
    for _ in range(20):
        worst = 0.0
        for k in comp_keys:
            jumps = np.abs(np.diff(g[k], prepend=g[k][0]))
            worst = max(worst, jumps.max())
        if worst <= max_dX:
            break
        for k in comp_keys:
            arr = g[k].copy()
            for i in range(1, n):
                d = arr[i] - arr[i-1]
                if abs(d) > max_dX:
                    # Spread the excess over the next 3 zones
                    excess = d - np.sign(d) * max_dX
                    for j in range(i, min(i+4, n)):
                        arr[j] -= excess * (1.0 - (j - i) / 4.0)
            g[k] = np.clip(arr, 0, None)
        # Re-normalise
        comp_sum = np.zeros(n)
        for k in comp_keys:
            comp_sum += g[k]
        comp_sum = np.maximum(comp_sum, 1e-30)
        for k in comp_keys:
            g[k] /= comp_sum

    if verbose:
        print(f"    After derivative capping (max_dX={max_dX}):")
        for k in comp_keys:
            jumps = np.abs(np.diff(g[k], prepend=g[k][0]))
            worst = np.argmax(jumps)
            worst_val = jumps[worst] if worst < n else 0
            print(f"      {k}: max jump = {worst_val:.4f}")

    return g


def reduce_surface_density(g, n_zones, rho_target, verbose=True):
    """Gradually reduce the density of the outermost ejecta zones.

    This creates a smoother density profile at the ejecta surface, helping
    the STELLA radiation solver transition through the optically thin regime.
    Instead of a sharp density edge, the outermost n_zones zones follow a
    log-linear ramp from their original density down to rho_target.

    The mass lost is negligible (outermost zones have very little mass),
    but the optical depth profile becomes smoother, preventing the solver
    from stalling (BADSTE / KFLAG=-2) when the photosphere approaches the
    surface.

    Parameters:
    -----------
    g : dict
        Remapped grid dict (modified in-place)
    n_zones : int
        Number of outermost zones to modify
    rho_target : float
        Target density at the outermost zone (g/cc)
    verbose : bool
        Print diagnostics
    """
    n = len(g['r'])
    if n_zones >= n:
        n_zones = n // 2

    # Original density at the boundary zone (first zone to modify)
    rho_start = g['rho'][n - n_zones]
    rho_end_orig = g['rho'][-1]

    if rho_target >= rho_start:
        if verbose:
            print(f"\n  [reduce_surface_density] rho_target ({rho_target:.2e}) >= "
                  f"rho_start ({rho_start:.2e}); skipping")
        return g

    # Log-linear ramp from rho_start to rho_target over n_zones
    log_rho_start = np.log10(rho_start)
    log_rho_target = np.log10(rho_target)
    new_log_rho = np.linspace(log_rho_start, log_rho_target, n_zones)

    # Apply the ramp to the outermost n_zones zones
    old_rho = g['rho'][n - n_zones:].copy()
    g['rho'][n - n_zones:] = 10**new_log_rho

    # Compute mass change
    r = g['r']
    r_outer_zones = r[n - n_zones:]
    r_inner_zones = np.concatenate([[r[n - n_zones - 1]], r_outer_zones[:-1]])
    dm_old = (4/3)*np.pi*(r_outer_zones**3 - r_inner_zones**3) * old_rho
    dm_new = (4/3)*np.pi*(r_outer_zones**3 - r_inner_zones**3) * g['rho'][n - n_zones:]
    M_lost = (dm_old.sum() - dm_new.sum()) / MSUN

    # Compute optical depth change
    kappa = 0.34
    dr_outer = np.diff(r_outer_zones)
    tau_old = np.sum(kappa * old_rho[:-1] * dr_outer)
    tau_new = np.sum(kappa * g['rho'][n-n_zones:-1] * dr_outer)

    if verbose:
        print(f"\n  [reduce_surface_density] Modifying outermost {n_zones} zones")
        print(f"    rho ramp: {rho_start:.2e} -> {rho_target:.2e} g/cc "
              f"({(log_rho_start - log_rho_target)/n_zones:.2f} dex/zone)")
        print(f"    Mass lost: {M_lost:.6f} Msun (negligible)")
        print(f"    Optical depth of modified region: {tau_old:.1f} -> {tau_new:.1f}")
        print(f"    Surface density: {rho_end_orig:.2e} -> {g['rho'][-1]:.2e} g/cc")

    return g


def add_csm_wind(g, n_wind_zones, mdot_msun_yr, v_wind_kms, t_wind_yr,
                 T_wind_K=5000.0, h_fraction=0.7, he_fraction=0.28,
                 density_profile='matched', power_law_index=None,
                 rho_floor=1.0e-20, verbose=True):
    """Append a low-density H envelope beyond the ejecta.

    The envelope provides optically thick material at large radii, preventing
    STELLA's radiation solver from stalling (BADSTE / Ncnd=0) when the
    photosphere exits the compact ejecta at ~10 days.

    The wind density at the inner edge matches the ejecta surface density
    (smooth transition, no density jump), then falls off as a steady-state
    1/r² wind:  rho(r) = rho_ej_surface * (r_ej / r)².

    The implied mass-loss rate is determined by the matching condition:
        Mdot = 4π r_ej² v_wind rho_ej_surface
    The user-provided mdot_msun_yr is used only as a diagnostic cross-check
    (the actual wind is set by the surface density match, not by Mdot).

    The outer radius is v_wind * t_wind, and the total wind mass is:
        M_wind = Mdot * (r_outer - r_ej) / v_wind = Mdot * t_wind

    Parameters:
    -----------
    g : dict
        Remapped grid dict (modified in-place — envelope zones appended)
    n_wind_zones : int
        Number of STELLA zones for the envelope
    mdot_msun_yr : float
        Nominal mass-loss rate (diagnostic only; actual Mdot is set by the
        surface-density matching condition).
    v_wind_kms : float
        Wind velocity.  Outer radius = v_wind * t_wind.
    t_wind_yr : float
        Envelope duration in years (determines outer radius and total mass)
    T_wind_K : float
        Envelope temperature in K (default: 5000)
    h_fraction : float
        H mass fraction in envelope (default: 0.7)
    he_fraction : float
        He mass fraction in envelope (default: 0.28)
    density_profile : str
        'matched' (default): 1/r² wind starting at ejecta surface density.
        'wind': 1/r² wind with density set by Mdot (may not match surface).
        'uniform': constant density.
        'extrapolate': power-law extension of the ejecta's OWN local
            log-log density slope at the surface (matches both density
            VALUE and GRADIENT at the interface -- no kink). This is the
            literature-standard approach (e.g. Dessart et al. 2018, A&A
            611, A6) for padding a progenitor's tenuous outer layers with
            negligible-mass material: rather than switching to a
            different (e.g. r^-2 steady-wind) slope at the boundary,
            simply extrapolate the star's own steep near-surface falloff
            outward by a modest factor in radius. Because the slope is
            steep (n~5-7 here vs. n=2 for a wind), the enclosed mass in
            the extension converges very quickly and stays negligible
            (<<0.01 Msun) even though it buys several more zones/decades
            of radius before the profile goes to hard vacuum -- avoiding
            the immediate (td~1e-10d) crash seen with a fixed r^-2 wind
            whose slope doesn't match the ejecta's at the interface.
        'crossover': upper envelope of the ejecta-slope extrapolation and
            the Mdot-set 1/r² wind -- the two extrapolated profiles join
            directly at their crossing radius (continuous density, single
            slope break).  Replaces the deprecated 'sigmoid' (which smeared
            the slope change over ~0.5 dex and produced light-curve kinks).
        'accelerating': beta-law v(r) wind, rho = max(ejecta-slope
            extrapolation, Mdot/4πr²v(r)).
        'sigmoid': DEPRECATED alias for 'crossover'.
    power_law_index : float or None
        Optional override of the fitted ejecta surface slope n (used by
        'crossover', 'extrapolate').
    rho_floor : float
        Minimum wind density in g/cc (default 1e-20).  Steep extrapolated
        tails underflow below STELLA's table limits; zones below the floor
        are clamped.
    verbose : bool
        Print diagnostics
    """
    r_ej = g['r'][-1]  # outer radius of ejecta
    rho_ej_surface = g['rho'][-1]  # density at ejecta surface
    T_ej_surface = g['T'][-1]  # temperature at ejecta surface
    v_wind = v_wind_kms * 1e5  # cm/s
    sec_per_year = 3.15e7
    r_outer = v_wind * t_wind_yr * sec_per_year
    if r_outer <= r_ej:
        r_outer = r_ej * 10  # fallback: extend 10x beyond ejecta

    kappa = 0.34  # electron scattering

    # 'sigmoid' is deprecated: the sigmoid slope-blend around the crossover
    # created an artificial ~0.5 dex transition region that showed up as
    # kinks in the decay-phase light curve.  It now maps to 'crossover',
    # which joins the two extrapolated profiles directly at r_cross.
    if density_profile == 'sigmoid':
        print("    WARNING: density_profile='sigmoid' is deprecated — using "
              "'crossover' (direct join of the extrapolated profiles, no "
              "sigmoid slope blend).")
        density_profile = 'crossover'

    # Compute the wind density at the inner edge (r_ej)
    if density_profile == 'crossover':
        # Crossover wind: extrapolate the ejecta's own steep power-law
        # density profile outward until it meets the physical 1/r² wind
        # set by Mdot, then switch directly to the wind profile.
        #
        # The previous 'matched' approach set rho_wind(r_ej) = rho_ej_surface,
        # which for a compact progenitor (rho_surf ~ 0.02 g/cc, r_ej ~ 1 R_sun)
        # gave an enormous wind mass (~2354 Msun!) because the 1/r² profile
        # integrated over a huge volume.  A real wind with Mdot=1e-5 Msun/yr
        # at v=100 km/s has rho ~ 7.5e-10 g/cc at r_ej — 2.3e7x lower.
        #
        # The two power laws cross at r_trans where:
        #   rho_ej_surf * (r/r_ej)^(-n_surf) = rho_wind_inner * (r/r_ej)^(-2)
        #   r_trans = r_ej * (rho_ej_surf / rho_wind_inner)^(1/(n_surf-2))
        #
        # Below r_trans: ejecta power law (steep, n~5-7).
        # Above r_trans: wind power law (1/r²).
        # rho(r) = max(rho_ej_ext, rho_wind): continuous at both the surface
        # and the crossover (the profiles are equal there by construction),
        # with a single slope break at r_trans where the physical wind
        # takes over from the ejecta extension — a real feature of the
        # joined profile, not an artificially smeared sigmoid transition.

        # Fit the local log-log slope from the last few real ejecta zones
        n_fit = min(6, len(g['r']) - 1)
        logr_fit = np.log(g['r'][-n_fit:])
        logrho_fit = np.log(g['rho'][-n_fit:])
        slope = np.polyfit(logr_fit, logrho_fit, 1)[0]  # d ln rho / d ln r
        n_surf = np.clip(-slope, 2.5, 12.0)

        # Allow override of the fitted slope for steeper CSM density falloff.
        # A steeper index (e.g. 12-15) makes the density drop faster from the
        # ejecta surface, reducing the CSM mass and minimizing the early-time
        # interaction bump while still providing optically thick material.
        if power_law_index is not None and power_law_index > 0:
            n_surf = float(power_law_index)

        # Physical wind density set by Mdot (NOT by surface matching)
        Mdot_implied = mdot_msun_yr * MSUN / sec_per_year  # g/s
        Mdot_implied_msun_yr = mdot_msun_yr
        rho_wind_inner = Mdot_implied / (4 * np.pi * r_ej**2 * v_wind)

        # Crossover radius where ejecta power law meets wind power law
        density_ratio = rho_ej_surface / rho_wind_inner
        log10_ratio = np.log10(density_ratio)
        r_cross = r_ej * density_ratio**(1.0 / (n_surf - 2.0))
        logr_cross = np.log10(r_cross)

        if verbose:
            print(f"    [crossover] fitted surface slope n = {n_surf:.2f} "
                  f"(d ln rho/d ln r = {slope:.2f})")
            print(f"    [crossover] wind density at r_ej: {rho_wind_inner:.3e} g/cc "
                  f"(Mdot={Mdot_implied_msun_yr:.1e} Msun/yr)")
            print(f"    [crossover] density contrast: {log10_ratio:.1f} dex")
            print(f"    [crossover] crossover at r={r_cross:.3e} cm "
                  f"({r_cross/r_ej:.0f}x r_ej, log10(r_cross)={logr_cross:.3f})")

    elif density_profile == 'accelerating':
        # Accelerating beta-law wind: rho(r) = Mdot / (4 pi r^2 v(r)) with
        # v(r) = v_wind * (1 - r_ej/r).  Because v -> 0 at the surface, the
        # wind density formally diverges there, so it crosses the ejecta's
        # extrapolated surface slope close to the edge.  The density profile
        # is the upper envelope of the two:
        #   rho(r) = max( rho_surf * (r_ej/r)^n_surf , Mdot/(4 pi r^2 v(r)) )
        # which is continuous at the interface (both start at rho_surf) and
        # at the crossover — no sigmoid stitching needed.  Only the slope
        # changes at the crossover (n_surf -> ~2).
        Mdot_implied = mdot_msun_yr * MSUN / sec_per_year  # g/s
        Mdot_implied_msun_yr = mdot_msun_yr
        rho_wind_inner = rho_ej_surface
        n_fit = min(6, len(g['r']) - 1)
        slope = np.polyfit(np.log(g['r'][-n_fit:]),
                           np.log(g['rho'][-n_fit:]), 1)[0]
        n_surf = np.clip(-slope, 2.5, 12.0)
        if verbose:
            print(f"    [accelerating] fitted surface slope n = {n_surf:.2f}")
            print(f"    [accelerating] rho continuous at edge "
                  f"(rho_surf = {rho_ej_surface:.3e} g/cc)")

    elif density_profile == 'extrapolate':
        # Fit the local log-log slope from the last few real ejecta zones
        # so the extension matches both density value AND gradient at the
        # interface (no kink -> no immediate solver crash).
        n_fit = min(6, len(g['r']) - 1)
        logr_fit = np.log(g['r'][-n_fit:])
        logrho_fit = np.log(g['rho'][-n_fit:])
        slope = np.polyfit(logr_fit, logrho_fit, 1)[0]  # d ln rho / d ln r (negative)
        n_surf = -slope
        n_surf = np.clip(n_surf, 2.5, 12.0)  # sanity bounds
        # Same optional slope override as 'crossover' (steeper index ->
        # less extended mass).
        if power_law_index is not None and power_law_index > 0:
            n_surf = float(power_law_index)
        rho_wind_inner = rho_ej_surface
        Mdot_implied = 4 * np.pi * r_ej**2 * v_wind * rho_ej_surface
        Mdot_implied_msun_yr = Mdot_implied * sec_per_year / MSUN
        if verbose:
            print(f"    [extrapolate] fitted surface slope n = {n_surf:.2f} "
                  f"(d ln rho/d ln r = {slope:.2f})")
    elif density_profile == 'matched':
        # Match ejecta surface density, then 1/r² falloff.
        # Implied Mdot = 4π r_ej² v_wind rho_ej_surface
        Mdot_implied = 4 * np.pi * r_ej**2 * v_wind * rho_ej_surface  # g/s
        Mdot_implied_msun_yr = Mdot_implied * sec_per_year / MSUN
        rho_wind_inner = rho_ej_surface
    elif density_profile == 'wind':
        # 1/r² wind with density set by the given Mdot
        Mdot_implied = mdot_msun_yr * MSUN / sec_per_year  # g/s
        Mdot_implied_msun_yr = mdot_msun_yr
        rho_wind_inner = Mdot_implied / (4 * np.pi * r_ej**2 * v_wind)
    else:  # uniform
        M_wind_nominal = mdot_msun_yr * t_wind_yr * MSUN
        vol = (4.0/3.0) * np.pi * (r_outer**3 - r_ej**3)
        rho_wind_inner = M_wind_nominal / vol
        Mdot_implied = 0.0
        Mdot_implied_msun_yr = 0.0

    # Wind mass and optical depth
    if density_profile in ('matched', 'wind'):
        # 1/r² wind: M = Mdot * (r_outer - r_ej) / v_wind = Mdot * t_wind
        M_wind = Mdot_implied * (r_outer - r_ej) / v_wind
        # tau = kappa * Mdot / (4π v) * (1/r_ej - 1/r_outer)
        tau_wind = kappa * Mdot_implied / (4 * np.pi * v_wind) * (1.0/r_ej - 1.0/r_outer)
    elif density_profile == 'crossover':
        # Crossover wind: mass is dominated by the ejecta power-law extension
        # from r_ej to r_cross (where it meets the wind), plus the wind mass
        # from r_cross to r_outer.  The ejecta extension mass is:
        #   M_ej_ext = 4π rho_surf r_ej³ ∫₁^(r_cross/r_ej) x^(2-n_surf) dx
        # The wind mass from r_cross to r_outer is:
        #   M_wind = Mdot * (r_outer - r_cross) / v_wind
        x_cross = r_cross / r_ej
        if abs(3.0 - n_surf) < 1e-6:
            integral_ej = np.log(x_cross)
        else:
            integral_ej = (x_cross**(3.0 - n_surf) - 1.0) / (3.0 - n_surf)
        M_ej_ext = 4 * np.pi * rho_ej_surface * r_ej**3 * integral_ej
        M_wind_part = Mdot_implied * (r_outer - r_cross) / v_wind
        M_wind = M_ej_ext + M_wind_part
        tau_ej_ext = kappa * rho_ej_surface * r_ej / (n_surf - 1.0) * (1.0 - x_cross**(1.0 - n_surf))
        tau_wind_part = kappa * Mdot_implied / (4 * np.pi * v_wind) * (1.0/r_cross - 1.0/r_outer)
        tau_wind = tau_ej_ext + tau_wind_part
    elif density_profile == 'accelerating':
        # Steady accelerating wind: M = Mdot * t_wind (independent of v(r));
        # tau from the terminal-velocity part of the wind.
        M_wind = Mdot_implied * t_wind_yr * sec_per_year
        tau_wind = kappa * Mdot_implied / (4 * np.pi * v_wind) * \
                   (1.0/r_ej - 1.0/r_outer)
    elif density_profile == 'extrapolate':
        # Power-law rho(r) = rho_ej_surface * (r/r_ej)^-n_surf.
        # M = 4*pi*rho_surf*r_ej^3 * integral_1^(r_outer/r_ej) x^(2-n_surf) dx
        x_max = r_outer / r_ej
        if abs(3.0 - n_surf) < 1e-6:
            integral = np.log(x_max)
        else:
            integral = (x_max**(3.0 - n_surf) - 1.0) / (3.0 - n_surf)
        M_wind = 4 * np.pi * rho_ej_surface * r_ej**3 * integral
        tau_wind = kappa * rho_ej_surface * r_ej / (n_surf - 1.0) * (1.0 - x_max**(1.0 - n_surf))
    else:
        M_wind = M_wind_nominal
        tau_wind = kappa * rho_wind_inner * (r_outer - r_ej)

    M_wind_msun = M_wind / MSUN
    density_contrast = rho_ej_surface / rho_wind_inner if rho_wind_inner > 0 else float('inf')

    if verbose:
        print(f"\n  [CSM Wind] Adding {n_wind_zones} wind zones beyond r={r_ej:.3e} cm")
        print(f"    v_wind = {v_wind_kms:.0f} km/s, t_wind = {t_wind_yr:.0f} yr")
        print(f"    Wind extent: {r_ej:.3e} -> {r_outer:.3e} cm ({r_outer/6.96e10:.0f} R_sun)")
        print(f"    Density profile: {density_profile}")
        print(f"    rho at ejecta surface: {rho_ej_surface:.3e} g/cc")
        print(f"    rho at inner (wind):   {rho_wind_inner:.3e} g/cc")
        print(f"    Density contrast: {density_contrast:.1f}x")
        if density_profile in ('matched', 'wind', 'crossover'):
            print(f"    Implied Mdot: {Mdot_implied_msun_yr:.4e} Msun/yr "
                  f"(config Mdot: {mdot_msun_yr:.4e})")
        print(f"    Wind mass (estimate): {M_wind_msun:.5f} Msun")
        print(f"    Wind optical depth (estimate): {tau_wind:.1f}")

    # Create envelope grid with smooth spacing transition from ejecta to wind.
    # The ejecta grid has very fine spacing (~0.006 dex/zone); the wind needs to
    # start at this spacing and gradually widen to cover the full wind extent.
    # dlogr starts at the ejecta spacing and ramps linearly (in zone index) to
    # the average wind spacing over a transition region.
    logr_ej = np.log10(r_ej)
    logr_outer = np.log10(r_outer)
    total_wind_dex = logr_outer - logr_ej

    # Ejecta grid spacing (dex per zone) at the surface
    n_ej = len(g['r'])
    if n_ej >= 3:
        ejecta_spacing = np.log10(g['r'][-1]) - np.log10(g['r'][-2])
    else:
        ejecta_spacing = 0.01  # fallback

    # Average wind spacing if uniformly distributed
    avg_wind_spacing = total_wind_dex / n_wind_zones

    # Build wind grid with linearly-ramped spacing: start at ejecta_spacing,
    # ramp to avg_wind_spacing over the first ~30% of wind zones.
    # This ensures no grid-spacing discontinuity at the ejecta-wind interface.
    n_trans = max(5, n_wind_zones // 3)  # transition over first 1/3 of wind
    logr_wind = []
    current_logr = logr_ej
    for i in range(n_wind_zones):
        frac = min(1.0, i / max(n_trans - 1, 1))
        spacing = ejecta_spacing + (avg_wind_spacing - ejecta_spacing) * frac
        current_logr += spacing
        if current_logr > logr_outer:
            current_logr = logr_outer
        logr_wind.append(current_logr)
    r_wind = 10**np.array(logr_wind)
    if np.any(np.diff(r_wind) <= 0):
        # Ramp overshot the outer radius before the last zone (duplicate radii
        # would give zero-mass zones).  Fall back to uniform log-spacing.
        r_wind = np.logspace(logr_ej, logr_outer, n_wind_zones + 1)[1:]

    if verbose:
        print(f"    [grid] Ejecta spacing: {ejecta_spacing:.4f} dex/zone")
        print(f"    [grid] Wind avg spacing: {avg_wind_spacing:.4f} dex/zone")
        print(f"    [grid] First wind zone at r={r_wind[0]:.3e} cm "
              f"(gap: {r_wind[0]/r_ej:.3f}x)")
        print(f"    [grid] Last wind zone at r={r_wind[-1]:.3e} cm")

    # Wind velocity: beta=1 law v(r) = v_wind * (1 - r_ej/r).  This is the
    # standard stellar-wind velocity profile: v -> 0 at the ejecta surface
    # (continuous with the static ejecta, no velocity discontinuity at the
    # interface) and v -> v_wind asymptotically at large radius.
    # Computed before the density profile because the 'accelerating' wind
    # density rho = Mdot/(4 pi r^2 v(r)) depends on it.
    v_wind_arr = v_wind * (1.0 - r_ej / r_wind)
    if verbose:
        print(f"    [velocity] v at inner wind: {v_wind_arr[0]/1e5:.2f} km/s, "
              f"outer: {v_wind_arr[-1]/1e5:.2f} km/s (beta=1 law)")

    # Envelope density profile
    if density_profile in ('matched', 'wind'):
        # 1/r² wind: rho(r) = rho_inner * (r_ej / r)²
        rho_wind = rho_wind_inner * (r_ej / r_wind)**2
    elif density_profile == 'crossover':
        # Direct extrapolation join: extend the ejecta power law outward and
        # the 1/r² wind inward, and take the upper envelope.
        #
        # rho_ej(r)   = rho_surface * (r_ej / r)^n_surf   [ejecta power law]
        # rho_wind(r) = Mdot / (4π r² v_wind)              [physical wind]
        #
        # The profiles are equal at r_cross by construction, so the envelope
        # is continuous in density with a single slope break (n_surf -> 2)
        # where the physical wind takes over.  No sigmoid stitching.
        rho_ej_ext = rho_ej_surface * (r_ej / r_wind)**n_surf
        rho_wind_w2 = rho_wind_inner * (r_ej / r_wind)**2
        rho_wind = np.maximum(rho_ej_ext, rho_wind_w2)

        if verbose:
            logr_wind = np.log10(r_wind)
            icross = np.argmin(np.abs(logr_wind - logr_cross))
            print(f"    [crossover] rho at inner: {rho_wind[0]:.3e} g/cc")
            print(f"    [crossover] rho at crossover: "
                  f"{rho_wind[icross]:.3e} g/cc (zone {icross+1})")
            print(f"    [crossover] rho at outer: {rho_wind[-1]:.3e} g/cc")

    elif density_profile == 'accelerating':
        # rho(r) = max( ejecta slope extrapolation, Mdot/(4 pi r^2 v(r)) ).
        # The accelerating-wind density diverges as v -> 0 at the surface,
        # so the upper envelope continues the ejecta slope from the edge
        # (continuous density AND slope there) until the two profiles cross,
        # then follows the wind.  Only a slope kink at the crossover — no
        # sigmoid stitching.
        rho_ej_ext = rho_ej_surface * (r_ej / r_wind)**n_surf
        rho_wind_acc = Mdot_implied / (4.0 * np.pi * r_wind**2 *
                                     np.maximum(v_wind_arr, 1.0))
        rho_wind = np.maximum(rho_ej_ext, rho_wind_acc)
        if verbose:
            cross = rho_wind_acc > rho_ej_ext
            icross = np.argmax(cross) if cross.any() else n_wind_zones - 1
            print(f"    [accelerating] ejecta-slope/wind crossover at "
                  f"r={r_wind[icross]:.3e} cm "
                  f"({r_wind[icross]/r_ej:.1f}x r_ej, zone {icross+1})")
            print(f"    [accelerating] rho at inner: {rho_wind[0]:.3e}, "
                  f"at crossover: {rho_wind[icross]:.3e}, "
                  f"outer: {rho_wind[-1]:.3e} g/cc")

    elif density_profile == 'extrapolate':
        # rho(r) = rho_ej_surface * (r_ej / r)^n_surf -- matches value AND
        # slope of the ejecta's own profile at r=r_ej (no kink).
        rho_wind = rho_ej_surface * (r_ej / r_wind)**n_surf
    else:
        rho_wind = np.full(n_wind_zones, rho_wind_inner)

    # Temperature: extrapolate the ejecta's own temperature power-law slope
    # into the wind zones, flooring at the cold wind temperature (T_wind_K)
    # at large radii.  This avoids the hard temperature discontinuity
    # (T_ej ~ 3e6 K -> T_wind = 5e3 K in one zone) that causes STELLA's
    # radiation transport solver to stall (KFLAG=-2, ifail 33) during
    # shock breakout when the photosphere reaches the ejecta/wind interface.
    #
    # Following Dessart et al. (2017, A&A 609, A6): the atmospheric/wind
    # temperature should follow the star's own temperature gradient, not a
    # hard jump to a fixed cold value.  We fit the local log-log temperature
    # slope from the last few ejecta zones and extrapolate; once the
    # extrapolation falls to the wind equilibrium temperature we hold
    # T_wind_K (single slope break where the two profiles cross — no sigmoid
    # blend region).
    n_fit = min(6, len(g['r']) - 1)
    logr_fit = np.log(g['r'][-n_fit:])
    logT_fit = np.log(np.maximum(g['T'][-n_fit:], 1.0))
    T_slope = np.polyfit(logr_fit, logT_fit, 1)[0]  # d ln T / d ln r (negative)
    # Clamp the slope to avoid unphysical extrapolation (T should decrease outward)
    T_slope = np.clip(T_slope, -12.0, -1.0)

    # Extrapolate T from the ejecta surface using the fitted slope
    T_extrap = T_ej_surface * (r_wind / r_ej)**T_slope

    # Crossover to the cold wind temperature: follow the extrapolated ejecta
    # temperature until it reaches T_wind_K, then hold the wind equilibrium
    # value.  T_extrap decreases monotonically (T_slope < 0), so the crossing
    # radius is analytic:
    #   r_Tcross = r_ej * (T_ej_surface / T_wind_K)^(-1/T_slope)
    if T_ej_surface > T_wind_K and T_slope < 0:
        r_Tcross = r_ej * (T_ej_surface / T_wind_K)**(-1.0 / T_slope)
    else:
        r_Tcross = r_ej  # already at/below wind temperature
    T_wind = np.maximum(T_extrap, T_wind_K)

    if verbose:
        print(f"    [T profile] fitted surface T slope: {T_slope:.2f} "
              f"(d ln T / d ln r)")
        print(f"    [T profile] T at ejecta surface: {T_ej_surface:.3e} K")
        print(f"    [T profile] T at inner wind (extrapolated): {T_wind[0]:.3e} K")
        print(f"    [T profile] T at outer wind (cold): {T_wind[-1]:.3e} K")
        print(f"    [T profile] crossover to T_wind at {r_Tcross:.3e} cm")

    # Density floor: STELLA cannot handle arbitrarily small densities
    # (steep 'extrapolate' tails underflow to ~1e-36).  Clamp the wind to
    # a physical vacuum-level floor (configurable via --wind-rho-floor /
    # csm_wind.rho_floor in the YAML).
    n_floored = int(np.sum(rho_wind < rho_floor))
    if n_floored > 0:
        rho_wind = np.maximum(rho_wind, rho_floor)
        if verbose:
            print(f"    [floor] {n_floored} wind zones clamped to "
                  f"rho = {rho_floor:.1e} g/cc")

    # Compute actual optical depth and mass of the wind region
    r_lo = np.concatenate([[r_ej], r_wind[:-1]])
    dm_wind = (4.0/3.0)*np.pi*(r_wind**3 - r_lo**3) * rho_wind
    M_wind_actual = np.sum(dm_wind)
    tau_actual = np.sum(kappa * rho_wind[:-1] * np.diff(r_wind))

    if verbose:
        print(f"    rho at inner: {rho_wind[0]:.3e} -> rho at outer: {rho_wind[-1]:.3e} g/cc")
        print(f"    Per-zone density drop: {(np.log10(rho_wind[0]) - np.log10(rho_wind[-1]))/n_wind_zones:.2f} dex")
        print(f"    Envelope mass (actual): {M_wind_actual/MSUN:.5f} Msun")
        print(f"    Envelope optical depth (actual): {tau_actual:.1f}")

    # Composition: H/He wind (with small traces of metals for opacity)
    fH_w = np.full(n_wind_zones, h_fraction)
    fHe_w = np.full(n_wind_zones, he_fraction)
    fC_w = np.zeros(n_wind_zones)
    fO_w = np.zeros(n_wind_zones)
    fSi_w = np.zeros(n_wind_zones)
    fCa_w = np.zeros(n_wind_zones)
    fFe_w = np.full(n_wind_zones, 1.0 - h_fraction - he_fraction)  # metals
    fNi56_w = np.zeros(n_wind_zones)

    # For the extrapolate/crossover profiles: no ejecta-zone blending needed.
    # The wind density starts at the ejecta surface value (continuous) and
    # follows the ejecta power law; for 'crossover' it switches to the wind
    # power law where the two extrapolated profiles meet.  The ejecta zones
    # keep their original densities.

    # Append wind zones to grid
    g['r'] = np.concatenate([g['r'], r_wind])
    g['rho'] = np.concatenate([g['rho'], rho_wind])
    g['T'] = np.concatenate([g['T'], T_wind])
    g['v'] = np.concatenate([g['v'], v_wind_arr])
    g['fH'] = np.concatenate([g['fH'], fH_w])
    g['fHe'] = np.concatenate([g['fHe'], fHe_w])
    g['fC'] = np.concatenate([g['fC'], fC_w])
    g['fO'] = np.concatenate([g['fO'], fO_w])
    g['fSi'] = np.concatenate([g['fSi'], fSi_w])
    g['fCa'] = np.concatenate([g['fCa'], fCa_w])
    g['fFe'] = np.concatenate([g['fFe'], fFe_w])
    g['fNi56'] = np.concatenate([g['fNi56'], fNi56_w])

    if verbose:
        print(f"    Total grid: {len(g['r'])} zones "
              f"(r=[{g['r'][0]:.3e}, {g['r'][-1]:.3e}] cm)")

    return g


# ── Initial velocity profile (kinetic explosion without EKO overwrite) ───────
def set_initial_velocity_profile(g, energy_e50, n_wind_zones=0, ai2=0.9,
                                 verbose=True):
    """Imprint a triangular-in-mass initial velocity profile on g['v'].

    Replicates the shape STELLA's EKO builds internally (begtt.f): velocity
    rises linearly in enclosed mass to a peak at mass fraction ai2 of the
    bound structure, then declines linearly to zero at the ejecta surface.
    Normalised so the kinetic energy of the bound structure (all zones
    except the wind) equals energy_e50 x 1e50 erg.

    Unlike EKO, which overwrites the entire .hyd velocity column, the
    constructed wind zones keep their beta-law velocity — the CSM is not
    launched at ejecta speed.  The ramp terminates at the ejecta surface so
    v -> 0 meets the wind's low base velocity without a jump.
    """
    r = g['r']; rho = g['rho']; n = len(r)
    r_lo = np.concatenate([[r[0]], r[:-1]])
    dm = np.maximum((4.0/3.0)*np.pi*(r**3 - r_lo**3)*rho, 0.0)
    m_enc = np.cumsum(dm)
    n_bound = n - n_wind_zones if n_wind_zones > 0 else n
    M_b = m_enc[n_bound-1]
    i2 = int(np.searchsorted(m_enc[:n_bound]/M_b, ai2))
    tri = np.zeros(n)
    tri[:i2+1] = m_enc[:i2+1]/m_enc[i2]
    tri[i2+1:n_bound] = (M_b - m_enc[i2+1:n_bound])/(M_b - m_enc[i2])
    ke = 0.5*np.sum(dm[:n_bound]*tri[:n_bound]**2)
    tri *= np.sqrt(energy_e50*1e50/ke)
    v_new = g['v'].copy()
    v_new[:n_bound] = tri[:n_bound]
    g['v'] = v_new
    if verbose:
        print(f"\n  [init-vel] triangular-in-mass initial velocity profile "
              f"(EBurst=EKO=0):")
        print(f"    KE(bound structure) = {energy_e50:.3f}e50 erg, "
              f"peak v = {tri[i2]/1e5:.0f} km/s at zone {i2} "
              f"(m={m_enc[i2]/MSUN:.3f} Msun)")
        print(f"    v -> 0 at ejecta surface; wind keeps its constructed "
              f"profile (outer edge v = {g['v'][n-1]/1e5:.1f} km/s)")


# ── Remap onto uniform STELLA grid ───────────────────────────────────────────
def remap_to_stella_grid(kep, mass_cut_msun, nzon, t_floor_K=3000.0, verbose=True,
                         h_surface_fraction=0.0, h_surface_mass_msun=0.0,
                         h_ext_fraction=0.0, h_ext_mass_msun=0.0,
                         h_ext_r_outer_rsun=0.0, h_ext_n_zones=50,
                         h_ext_T_surface_K=6000.0, h_ext_T_base_K=None,
                         h_ext_n_rho=3.5, h_ext_n_T=2.0,
                         rho_min=0.0, n_wind_zones=0, wind_mdot=0.0,
                         wind_velocity=100.0, wind_duration_yr=1000.0,
                         wind_T_K=5000.0, wind_density_profile='uniform',
                         wind_power_law_index=None, wind_rho_floor=1.0e-20,
                         reduce_surf_n_zones=0,
                         reduce_surf_rho_target=0.0,
                         initial_velocity_energy=0.0):
    """Excise inner compact remnant and remap onto nzon log-spaced STELLA zones.

    Mass coordinates above mass_cut_msun are kept.
    The new grid is log-spaced in radius from r_cut to r_max.

    If h_surface_fraction > 0, a thin surface layer of H is added to the
    outermost zones covering h_surface_mass_msun of ejecta mass.

    If rho_min > 0, the outer boundary r_max is truncated at the outermost
    KEPLER zone with density >= rho_min.  This removes the extended
    low-density atmosphere (negligible mass) that causes STELLA
    radiation-transport instabilities when the photosphere reaches it.

    Returns dict of arrays of length nzon.
    """
    mass_cut_g = mass_cut_msun * MSUN

    # Find first zone whose cumulative mass exceeds the cut
    idx_cut = np.searchsorted(kep['m_g'], mass_cut_g)
    if idx_cut >= len(kep['m_g']):
        raise ValueError(f"mass_cut_msun={mass_cut_msun} exceeds total mass!")

    r_cut   = kep['r_cm'][idx_cut]
    r_max   = kep['r_cm'][-1]
    m_cut_g = kep['m_g'][idx_cut]

    # Trim outer low-density atmosphere if requested
    if rho_min > 0:
        rho_ejecta = kep['rho'][idx_cut:]
        r_ejecta   = kep['r_cm'][idx_cut:]
        keep_rho   = rho_ejecta >= rho_min
        if np.any(keep_rho):
            last_keep = np.where(keep_rho)[0][-1]
            r_max_new = r_ejecta[last_keep]
            m_trimmed = (kep['m_g'][-1] - kep['m_g'][idx_cut + last_keep]) / MSUN
            if verbose:
                print(f"\n  rho_min trimming: {rho_min:.2e} g/cm^3")
                print(f"    r_max: {r_max:.3e} -> {r_max_new:.3e} cm "
                      f"(trimmed {m_trimmed:.5f} Msun of low-density atmosphere)")
            r_max = r_max_new
        elif verbose:
            print(f"\n  rho_min={rho_min:.2e}: entire ejecta below threshold; keeping full grid")

    if verbose:
        print(f"\nMass cut: {mass_cut_msun:.3f} Msun (zone {idx_cut+1}/{len(kep['m_g'])})")
        print(f"  r_cut = {r_cut:.3e} cm")
        print(f"  r_max = {r_max:.3e} cm")
        print(f"  Ejecta mass above cut: {(kep['m_g'][-1] - m_cut_g)/MSUN:.4f} Msun")

    # Source arrays (ejecta only, zones idx_cut : end)
    r_src   = kep['r_cm'][idx_cut:]
    rho_src = kep['rho'][idx_cut:]
    T_src   = kep['T_K'][idx_cut:]
    v_src   = kep['v_cms'][idx_cut:]   # will be zeroed for bomb mode

    # Get compositions (7-element: H, He, C, O, Si, Ca, Fe + Ni56)
    fH, fHe, fC, fO, fSi, fCa, fFe, fNi56 = build_stella_composition(kep)
    fH_src  = fH[idx_cut:];  fHe_src = fHe[idx_cut:]
    fC_src  = fC[idx_cut:];  fO_src  = fO[idx_cut:]
    fSi_src = fSi[idx_cut:]; fCa_src = fCa[idx_cut:]
    fFe_src = fFe[idx_cut:]; fNi56_src = fNi56[idx_cut:]

    # New STELLA grid: log-spaced in radius
    # If wind zones are requested, reduce ejecta zones to make room
    # If extended H envelope zones are requested, also make room for those
    n_extra = n_wind_zones + (h_ext_n_zones if h_ext_mass_msun > 0 else 0)
    n_ejecta = nzon - n_extra if n_extra > 0 else nzon
    r_new = np.logspace(np.log10(r_cut), np.log10(r_max), n_ejecta)

    def interp(y_src):
        """Log-linear interpolation onto new grid, clipped to source range."""
        ip = interp1d(np.log10(r_src), y_src, kind='linear',
                      bounds_error=False, fill_value=(y_src[0], y_src[-1]))
        return ip(np.log10(r_new))

    def log_interp(y_src):
        """Log-log interpolation for positive definite quantities."""
        log_y = np.log10(np.maximum(y_src, 1e-300))
        ip = interp1d(np.log10(r_src), log_y, kind='linear',
                      bounds_error=False, fill_value=(log_y[0], log_y[-1]))
        return 10.0 ** ip(np.log10(r_new))

    rho_new = log_interp(rho_src)
    T_new   = np.maximum(log_interp(T_src), t_floor_K)
    # Zero pre-SN infall velocities — explosion velocities are imposed later
    # (EBurst thermal bomb, EKO, or --initial-velocity-energy .hyd profile)
    v_new   = np.zeros(n_ejecta)

    fH_new  = interp(fH_src);  fHe_new = interp(fHe_src)
    fC_new  = interp(fC_src);  fO_new  = interp(fO_src)
    fSi_new = interp(fSi_src); fCa_new = interp(fCa_src)
    fFe_new = interp(fFe_src); fNi56_new = interp(fNi56_src)

    # Re-normalise composition after interpolation (7 elements)
    comp_sum = fH_new + fHe_new + fC_new + fO_new + fSi_new + fCa_new + fFe_new
    comp_sum = np.maximum(comp_sum, 1e-10)
    fH_new  /= comp_sum; fHe_new /= comp_sum; fC_new  /= comp_sum
    fO_new  /= comp_sum; fSi_new /= comp_sum; fCa_new /= comp_sum
    fFe_new /= comp_sum

    if verbose:
        print(f"\nRemapped to {nzon} STELLA zones (log-spaced in r):")
        print(f"  r:   {r_new[0]:.3e} – {r_new[-1]:.3e} cm")
        print(f"  rho: {rho_new[0]:.3e} – {rho_new[-1]:.3e} g/cc")
        print(f"  T:   {T_new[0]:.3e} – {T_new[-1]:.3e} K")
        # Mass-weighted average composition
        dr = np.concatenate([[r_new[1]-r_new[0]], np.diff(r_new)])
        r_lo = np.concatenate([[r_cut*0.9], r_new[:-1]])
        dm = (4/3)*np.pi*(r_new**3 - r_lo**3) * rho_new
        dm = np.maximum(dm, 0)
        M_tot = dm.sum()
        print(f"  Total ejecta mass (remapped): {M_tot/MSUN:.4f} Msun")
        for name, f_arr in [('H', fH_new), ('He', fHe_new),('C', fC_new),('O', fO_new),
                             ('Si', fSi_new),('Ca', fCa_new),('Fe', fFe_new)]:
            print(f"  <{name}>_m = {(f_arr * dm).sum() / M_tot:.4f}")

    # Build the return dict (7-element composition: H, He, C, O, Si, Ca, Fe + Ni56)
    g_out = {
        'r': r_new, 'rho': rho_new, 'T': T_new, 'v': v_new,
        'fH': fH_new, 'fHe': fHe_new, 'fC': fC_new,
        'fO': fO_new, 'fSi': fSi_new, 'fCa': fCa_new, 'fFe': fFe_new,
        'fNi56': fNi56_new,
        'r_cut': r_cut, 'mass_cut_g': m_cut_g,
    }

    # Add H surface layer if requested
    if h_surface_fraction > 0 and h_surface_mass_msun > 0:
        print(f"\n  Adding H surface layer (fraction={h_surface_fraction}, mass={h_surface_mass_msun} Msun)...")
        add_h_surface_layer(g_out, h_surface_fraction, h_surface_mass_msun)
        if verbose:
            dm2 = (4/3)*np.pi*(r_new**3 - np.concatenate([[r_cut], r_new[:-1]])**3) * rho_new
            dm2 = np.maximum(dm2, 0)
            M_tot2 = dm2.sum()
            for name, f_arr in [('H', g_out['fH']), ('He', g_out['fHe']),
                                ('C', g_out['fC']), ('O', g_out['fO']),
                                ('Si', g_out['fSi']), ('Ca', g_out['fCa']),
                                ('Fe', g_out['fFe'])]:
                print(f"  <{name}>_m = {(f_arr * dm2).sum() / M_tot2:.4f}")

    # Add extended H envelope if requested (Bersten 2012 approach for SN 2011dh)
    # This extends the progenitor radius from ~1 Rsun to ~200 Rsun with a
    # low-mass H-rich envelope, matching the YSG progenitor.
    if h_ext_fraction > 0 and h_ext_mass_msun > 0 and h_ext_r_outer_rsun > 0:
        print(f"\n  Adding extended H envelope (fraction={h_ext_fraction}, "
              f"mass={h_ext_mass_msun} Msun, R_outer={h_ext_r_outer_rsun} Rsun)...")
        add_extended_h_envelope(g_out, h_ext_fraction, h_ext_mass_msun,
                                h_ext_r_outer_rsun, n_envelope_zones=h_ext_n_zones,
                                T_surface_K=h_ext_T_surface_K,
                                T_base_K=h_ext_T_base_K,
                                n_rho=h_ext_n_rho, n_T=h_ext_n_T,
                                verbose=verbose)

    # Reduce surface density BEFORE adding wind (smoother outer ejecta profile)
    # This must come before add_csm_wind so the wind transition starts from
    # the reduced density, not the original high ejecta surface density.
    if reduce_surf_n_zones > 0 and reduce_surf_rho_target > 0:
        reduce_surface_density(g_out, reduce_surf_n_zones,
                               reduce_surf_rho_target, verbose=verbose)

    # Add CSM wind if requested ('extrapolate' mode needs no Mdot -- it
    # derives its density purely from the ejecta's own surface value/slope)
    if n_wind_zones > 0 and (wind_mdot > 0 or wind_density_profile
                             in ('extrapolate', 'crossover', 'sigmoid')):
        wind_h_frac = h_surface_fraction if h_surface_fraction > 0 else 0.7
        wind_he_frac = 1.0 - wind_h_frac - 0.02  # small metals fraction
        add_csm_wind(g_out, n_wind_zones, wind_mdot, wind_velocity,
                     wind_duration_yr, T_wind_K=wind_T_K,
                     h_fraction=wind_h_frac, he_fraction=wind_he_frac,
                     density_profile=wind_density_profile,
                     power_law_index=wind_power_law_index,
                     rho_floor=wind_rho_floor,
                     verbose=verbose)

    # Imprint the initial velocity field (kinetic explosion) AFTER the wind
    # is added so the wind zones keep their constructed beta-law velocity.
    # Use with EBurst=0 and EKO=0 in the .dat: the .hyd carries the explosion.
    if initial_velocity_energy > 0:
        set_initial_velocity_profile(g_out, initial_velocity_energy,
                                     n_wind_zones=n_wind_zones,
                                     verbose=verbose)

    return g_out


# ── .hyd writer (matches existing STELLA format) ─────────────────────────────
def write_hyd(fname, g, mass_cut_msun):
    """Write STELLA *.hyd initial conditions file."""
    r = g['r'].copy(); rho = g['rho']; T = g['T']; v = g['v']
    nzon = len(r)

    # Fix zero-mass zones (r[i] <= r_lo[i]) that cause STELLA's Jacobian to
    # become structurally singular (NaN velocities -> SIGSEGV in eddi.f).
    # The log-spaced grid starts at r_cut, so zone 1 has r[0]==r_cut -> dm=0.
    # Perturb inward to give non-zero mass (same approach as convert_profile_to_stella.py).
    r_lo = np.concatenate([[g['r_cut']], r[:-1]])
    for i in range(nzon):
        if r[i] <= r_lo[i]:
            r[i] = r_lo[i] * (1.0 + 1e-6)
            print(f"  WARNING: perturbed zone {i+1} radius to {r[i]:.6e} "
                  f"(was <= inner boundary {r_lo[i]:.6e})")

    r_lo = np.concatenate([[g['r_cut']], r[:-1]])
    dm   = (4/3)*np.pi * (r**3 - r_lo**3) * rho
    dm   = np.maximum(dm, 0.0)
    M_r  = np.cumsum(dm)

    rcen    = g['r_cut'] * 0.5
    rho_cen = rho[0]

    with open(fname, 'w') as f:
        # timetobegin=0 (bomb mode), BM1=mass_cut_msun
        f.write(f"   0.000E+00{nzon:6d}{mass_cut_msun:13.5E}{rcen:13.5E}{rho_cen:13.5E}\n")
        for i in range(nzon):
            dMr = dm[i] / MSUN
            Mr  = M_r[i] / MSUN + mass_cut_msun
            f.write(f" {i+1:4d}{dMr:12.4E}{r[i]:15.7E}{rho[i]:12.4E}"
                    f"{T[i]:12.4E}{v[i]:12.4E}{Mr:12.4E}{0.0:12.4E}\n")

    M_total = (M_r[-1] / MSUN) + mass_cut_msun
    print(f"Written {fname}  ({nzon} zones, M_total={M_total:.4f} Msun, "
          f"r=[{r[0]:.2e},{r[-1]:.2e}] cm)")


# ── .abn writer ──────────────────────────────────────────────────────────────
def write_abn(fname, g, ni56_mass_msun=0.0, ni56_dist_msun=0.0):
    """Write STELLA *.abn composition file.

    Columns: km  0  0  0  H  He  C  N  O  Ne  Na  Mg  Al  Si  S  Ar  Ca  Fe  Ni58  Ni56

    7 stable elements written: H, He, C, O, Si, Ca, FePeak.
    Ca is written as its own column for CMFGEN spectral modelling.

    Ni56 (XNifor last column) is set to 0 — radioactive Ni56 comes from the
    .xni file only.  eve2 reads XNifor and adds it to the Fe abundance for
    opacity (evesn2.f: YABUN(Fe) = (YABUN(Fe) + XNifor) / AZ(Fe)); placing
    0.251 Msun of Ni56 there increases inner-zone Fe opacity by ~2.5x and
    causes the radiation transport solver to stall (KFLAG=-2 at td~1.2d).
    """
    nzon = len(g['r'])
    with open(fname, 'w') as f:
        for i in range(nzon):
            f.write(f"{i+1:4d}{0.0:10.3E}{0.0:10.3E}{0.0:10.3E}"
                    f"{g['fH'][i]:10.3E}{g['fHe'][i]:10.3E}{g['fC'][i]:10.3E}"
                    f"{0.0:10.3E}{g['fO'][i]:10.3E}"    # N=0
                    f"{0.0:10.3E}{0.0:10.3E}{0.0:10.3E}{0.0:10.3E}"  # Ne Na Mg Al
                    f"{g['fSi'][i]:10.3E}"               # Si
                    f"{0.0:10.3E}{0.0:10.3E}{g['fCa'][i]:10.3E}"  # S Ar Ca
                    f"{g['fFe'][i]:10.3E}"               # FePeak (without Ca)
                    f"{0.0:10.3E}{0.0:10.3E}\n")         # Ni58 Ni56 (both 0)
    print(f"Written {fname}  ({nzon} zones, Ni56=0 - radioactive Ni56 from .xni file only)")


# ── .xni writer (Ni56 distribution) ──────────────────────────────────────────
def write_xni(fname, g, ni56_mass_msun, ni56_dist_msun, nzon):
    """Write STELLA binary .xni file with Ni56 mass fractions.

    Ni56 is placed in the innermost zones covering ni56_dist_msun of ejecta.

    The .xni file is Fortran unformatted (form='unformatted'), containing
    a single record of Nzon REAL*8 (double precision) values representing
    the Ni56 mass fraction in each zone.  STELLA reads this with:
        read(28) (XNifor(i), i=1, Nzon)
    where XNifor is REAL*8 (from IMPLICIT REAL*8(A-H,O-Z)).

    Fortran unformatted format: [4-byte rec_len] [data] [4-byte rec_len]
    """
    import struct

    # Use the same radius perturbation as write_hyd to ensure dm values match.
    # write_hyd perturbs zones where r[i] <= r_lo[i] (zone 1 starts at r_cut).
    r = g['r'].copy(); rho = g['rho']
    r_lo = np.concatenate([[g['r_cut']], r[:-1]])
    for i in range(nzon):
        if r[i] <= r_lo[i]:
            r[i] = r_lo[i] * (1.0 + 1e-6)

    r_lo = np.concatenate([[g['r_cut']], r[:-1]])
    dm   = np.maximum((4/3)*np.pi*(r**3 - r_lo**3)*rho, 0.0)
    m_cum = np.cumsum(dm) / MSUN   # cumulative ejecta mass in Msun

    # Find how many inner zones cover ni56_dist_msun
    ni56_zones = np.where(m_cum <= ni56_dist_msun)[0]
    if len(ni56_zones) == 0:
        ni56_zones = np.array([0])
    # Skip zone 1 (index 0) to avoid spurious Ni56 in eve2's core zone
    ni56_zones = ni56_zones[ni56_zones >= 1]
    if len(ni56_zones) == 0:
        ni56_zones = np.array([1])

    n_ni_zones = len(ni56_zones)
    m_ni56_region = m_cum[ni56_zones[-1]]  # actual mass in those zones

    # Ni56 mass fraction within the Ni56 region
    fni56 = ni56_mass_msun / max(m_ni56_region, 1e-10)
    fni56 = min(fni56, 1.0)

    print(f"\n.xni Ni56 distribution:")
    print(f"  Total Ni56: {ni56_mass_msun:.3f} Msun")
    print(f"  Ni56 region: {n_ni_zones} zones (2..{ni56_zones[-1]+1}) covering {m_ni56_region:.3f} Msun ejecta")
    print(f"  Ni56 mass fraction in those zones: {fni56:.4f}")
    print(f"  (zone 1 skipped — eve2 merges remnant core into .mod zone 1)")

    # Build Ni56 mass fraction array (REAL*8, one value per zone)
    xnifor = np.zeros(nzon, dtype=np.float64)
    ni56_zone_set = set(ni56_zones.tolist())
    for i in range(nzon):
        if i in ni56_zone_set:
            xnifor[i] = fni56

    # Write as Fortran unformatted: [rec_len] [Nzon * REAL*8] [rec_len]
    rec_len = nzon * 8  # bytes
    with open(fname, 'wb') as f:
        f.write(struct.pack('i', rec_len))
        f.write(struct.pack(f'{nzon}d', *xnifor))
        f.write(struct.pack('i', rec_len))

    # Verify
    actual_ni56 = (dm * xnifor).sum() / MSUN
    print(f"  Verification: Ni56 mass = {actual_ni56:.6f} Msun (target: {ni56_mass_msun:.3f})")
    print(f"Written {fname}  ({4 + rec_len + 4} bytes, Fortran unformatted REAL*8)")


# ── input_profile writer ([DEPRECATED] HERACLES-format, needed by opacity table generator) ─
def write_input_profile(fname, g):
    """Write input_profile in the 13-column format expected by convert_profile_to_stella.py.

    Columns: r_cm  rho  p_cgs  u_cms  T_K  fH  fHe  fC  fO  fSi  fCa  fFe  fNi
    C and O are kept as separate elements (not merged).

    Pressure is reconstructed from ideal gas + radiation: p = rho*kT/(mu*mH) + a*T^4/3.
    """
    r = g['r']; rho = g['rho']; T = g['T']

    # Mean molecular weight for fully ionised mix
    fH = g['fH']; fHe = g['fHe']; fC = g['fC']; fO = g['fO']
    fSi = g['fSi']; fFe = g['fFe']; fCa = g.get('fCa', np.zeros_like(fH))
    fNi = g.get('fNi56', np.zeros_like(fH))
    mu = 1.0 / (fH/1.0 + fHe/4.0 + fC/12.0 + fO/16.0 + fSi/28.0 + fCa/40.0 + fFe/56.0 + fNi/58.0)
    mu = np.maximum(mu, 0.5)
    p_gas = rho * K_B * T / (mu * M_H)
    p_rad = A_RAD * T**4 / 3.0
    p = p_gas + p_rad

    header = ("# input_profile generated from Kepler presn by "
              "convert_kepler_presn_to_stella.py\n"
              "# columns: r_cm rho_cgs p_cgs u_cms T_K fx_H fx_He fx_C fx_O fx_Si fx_Ca fx_Fe fx_Ni\n"
              "# C and O are separate elements (not merged).\n")
    with open(fname, 'w') as f:
        f.write(header)
        for i in range(len(r)):
            f.write(f"{r[i]:.8e} {rho[i]:.8e} {p[i]:.8e} {g['v'][i]:.8e} "
                    f"{T[i]:.8e} "
                    f"{fH[i]:.8e} {fHe[i]:.8e} {fC[i]:.8e} {fO[i]:.8e} "
                    f"{fSi[i]:.8e} {fCa[i]:.8e} {fFe[i]:.8e} {fNi[i]:.8e}\n")
    print(f"Written {fname}  ({len(r)} zones, 13-column input_profile)")


# ── .dat writer (STELLA control file) ────────────────────────────────────────
def write_dat(fname, model_name, target_days,
              explosion_energy=10.0, explosion_mass=0.5, explosion_time=1.0,
              explosion_delay=0.0, initial_kinetic_energy=0.0,
              nstmax=2000000,
              eps=0.003, hmin=1e-14, hmax=5e5, meth=2,
              iout=10, nout=100,
              fit_tau=15.0, tau_tol=1.3, aq=2.0, bq=1.0,
              enable_convection=False, disable_scattering=False, givdtl=True,
              wacc_r=1.0, wacc_v=1.0, wacc_t=0.5, wacc_rad=1.0,
              temp_floor=0.03, floor_rad=0.01,
              snapshot_interval=10.0,
              ni56_mass=0.075, ni56_dist=None,
              custom_output_times=None):
    """Write STELLA *.dat control file for a Type IIb explosion.

    SN_A convention: explosion_energy (EBurst) > 0 is the thermal bomb that
    blows up the KEPLER model; initial_kinetic_energy (EKO) stays 0 so the
    .hyd velocity column is preserved.  initial_velocity_energy is a
    bomb-free alternative that writes the velocity into .hyd instead —
    it must never be combined with a nonzero EBurst.

    Delegates to convert_profile_to_stella.write_dat so the .dat carries the
    full resolved parameter set (single source of truth for the .dat format).
    AMNI (Ni56 contamination mass) is set to ni56_dist if given, else
    explosion_mass — matching the historical behaviour of this writer.
    """
    from convert_profile_to_stella import write_dat as _write_dat
    _write_dat(fname, model_name, target_days,
               nout=nout,
               snapshot_interval_days=snapshot_interval,
               explosion_energy=explosion_energy,
               explosion_mass=explosion_mass,
               explosion_time=explosion_time,
               explosion_delay=explosion_delay,
               initial_kinetic_energy=initial_kinetic_energy,
               nstmax=nstmax,
               disable_scattering=disable_scattering,
               eps=eps, hmin=hmin, hmax=hmax, meth=meth,
               enable_convection=enable_convection,
               fit_tau=fit_tau, tau_tol=tau_tol,
               givdtl=givdtl, custom_output_times=custom_output_times,
               wacc_r=wacc_r, wacc_v=wacc_v, wacc_t=wacc_t, wacc_rad=wacc_rad,
               temp_floor=temp_floor, floor_rad=floor_rad,
               iout=iout, aq=aq, bq=bq,
               amni=ni56_dist if ni56_dist is not None else explosion_mass,
               xmni=ni56_mass)
    print(f"  Explosion: EBurst={explosion_energy}e50 erg, AMHT={explosion_mass} Msun, "
          f"tBurst={explosion_time} s")


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description=__doc__,
                 formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('presn_file',     help='Kepler pre-SN file (he8.00-presn)')
    parser.add_argument('model_name',     help='Model name (e.g. SN_A_IIb)')
    parser.add_argument('--mass-cut',     type=float, default=1.5,
                        help='Compact remnant mass cut in Msun (default: 1.5)')
    parser.add_argument('--nzon',         type=int,   default=500,
                        help='Number of STELLA zones (default: 500)')
    parser.add_argument('--target-days',  type=float, default=400.0,
                        help='Simulation duration in days (default: 400)')
    parser.add_argument('--explosion-energy', type=float, default=10.0,
                        help='EBurst in units of 1e50 erg (default: 10 = 1e51 erg). '
                             'SN_A convention: this thermal bomb blows up the '
                             'KEPLER model. Keep >0 for SN_A; the SN_B '
                             'embedded-shell path (prepare_SN_B_model.py) uses 0.')
    parser.add_argument('--explosion-mass',   type=float, default=0.5,
                        help='AMHT heated core mass in Msun (default: 0.5)')
    parser.add_argument('--explosion-time',   type=float, default=1.0,
                        help='tBurst energy injection duration in seconds (default: 1)')
    parser.add_argument('--initial-kinetic-energy', type=float, default=0.0,
                        help='EKO in units of 1e50 erg (default: 0 = preserve .hyd velocities). '
                             'KEEP 0 for SN_A: EKO>0 makes STELLA overwrite the whole .hyd '
                             'velocity column (including the constructed wind) with a '
                             'triangular-in-mass profile. Total energy = EKO + EBurst.')
    parser.add_argument('--initial-velocity-energy', type=float, default=0.0,
                        help='Bomb-free escape hatch: writes a triangular-in-mass '
                             'initial velocity profile into the .hyd carrying this '
                             'much KE (1e50 erg); wind zones keep their constructed '
                             'velocity. Requires --explosion-energy 0 and '
                             '--initial-kinetic-energy 0 — never combine with EBurst.')
    parser.add_argument('--ni56-mass',    type=float, default=0.1,
                        help='Total Ni56 mass in Msun (default: 0.1)')
    parser.add_argument('--ni56-dist',    type=float, default=0.2,
                        help='Ejecta mass (Msun) over which Ni56 is distributed (default: 0.2)')
    # ── .dat solver/control parameters (forwarded to the full-parameter
    #    write_dat in convert_profile_to_stella.py) ──
    parser.add_argument('--explosion-delay', type=float, default=0.0,
                        help='tstart / explosion delay in days (default: 0)')
    parser.add_argument('--snapshot-interval', type=float, default=5.0,
                        help='Days between .swd snapshot epochs (default: 5)')
    parser.add_argument('--custom-output-times', type=str, default=None,
                        help='Comma-separated list of output epochs (days); overrides '
                             '--snapshot-interval epoch generation')
    parser.add_argument('--nstmax',       type=int,   default=2000000,
                        help='NSTMAX max timesteps (default: 2000000)')
    parser.add_argument('--eps',          type=float, default=0.003,
                        help='EPS stiff-solver accuracy (default: 0.003)')
    parser.add_argument('--hmin',         type=float, default=1e-14,
                        help='HMIN min timestep in s (default: 1e-14)')
    parser.add_argument('--hmax',         type=float, default=5e5,
                        help='HMAX max timestep in s (default: 5e5)')
    parser.add_argument('--meth',         type=int,   default=2,
                        help='Integrator: 1=Adams, 2=GEAR, 3=GEAR+BGH (default: 2)')
    parser.add_argument('--iout',         type=int,   default=10,
                        help='IOUT output frequency (default: 10)')
    parser.add_argument('--nout',         type=int,   default=100,
                        help='NOUT .swd write cadence in steps (default: 100; must be <200)')
    parser.add_argument('--fit-tau',      type=float, default=15.0,
                        help='FitTau parameter for radiation transport (default: 15.0)')
    parser.add_argument('--tau-tol',      type=float, default=1.3,
                        help='TAUTOL parameter (default: 1.3)')
    parser.add_argument('--aq',           type=float, default=2.0,
                        help='AQ artificial viscosity (default: 2.0; higher for shocks)')
    parser.add_argument('--bq',           type=float, default=1.0,
                        help='BQ artificial viscosity (default: 1.0)')
    parser.add_argument('--wacc-r',       type=float, default=1.0)
    parser.add_argument('--wacc-v',       type=float, default=1.0)
    parser.add_argument('--wacc-t',       type=float, default=0.5)
    parser.add_argument('--wacc-rad',     type=float, default=1.0)
    parser.add_argument('--temp-floor',   type=float, default=0.03,
                        help='FLOOR(T) in 1e5 K units (default: 0.03 = 3000 K)')
    parser.add_argument('--floor-rad',    type=float, default=0.01,
                        help='FLOOR(RADIAT) in 1e5 K units (default: 0.01 = 1000 K)')
    parser.add_argument('--enable-convection', action='store_true', default=False)
    parser.add_argument('--disable-scattering', action='store_true', default=False)
    parser.add_argument('--no-givdtl',    dest='givdtl', action='store_false', default=True,
                        help='Disable Givdtl (default: enabled)')
    parser.add_argument('--t-floor',      type=float, default=3000.0,
                        help='Initial temperature floor in K (default: 3000)')
    parser.add_argument('--h-surface-fraction', type=float, default=0.0,
                        help='H mass fraction in the thin surface layer (default: 0 = no H layer). '
                             'Other mass fractions are rescaled by (1 - h_surface_fraction).')
    parser.add_argument('--h-surface-mass-msun', type=float, default=0.0,
                        help='Mass of the H surface layer in Msun (controls how many outer zones get H). '
                             'Default: 0 = no H layer.')
    # Extended H envelope (Bersten 2012 approach for SN 2011dh)
    parser.add_argument('--h-ext-fraction', type=float, default=0.0,
                        help='H mass fraction in the extended envelope (default: 0 = no envelope). '
                             'Used with --h-ext-mass-msun and --h-ext-r-outer-rsun.')
    parser.add_argument('--h-ext-mass-msun', type=float, default=0.0,
                        help='Total mass of the extended H envelope in Msun (default: 0 = no envelope). '
                             'Typical SN 2011dh: 0.1 Msun.')
    parser.add_argument('--h-ext-r-outer-rsun', type=float, default=0.0,
                        help='Outer radius of the extended H envelope in Rsun (default: 0 = no envelope). '
                             'Typical SN 2011dh YSG: 200 Rsun.')
    parser.add_argument('--h-ext-n-zones', type=int, default=50,
                        help='Number of STELLA zones for the extended H envelope (default: 50).')
    parser.add_argument('--h-ext-T-surface-K', type=float, default=6000.0,
                        help='Temperature at the outer edge of the envelope in K (default: 6000, YSG).')
    parser.add_argument('--h-ext-T-base-K', type=float, default=None,
                        help='Temperature at the base of the envelope in K (default: 50000). '
                             'A real YSG H envelope has T_base ~ 30-100 kK at the He/H interface. '
                             'Do NOT use the KEPLER He core surface T (~5e6 K).')
    parser.add_argument('--h-ext-n-rho', type=float, default=3.5,
                        help='Density power-law exponent for the envelope (default: 3.5). '
                             'rho(r) = rho_base * (r/r_base)^(-n_rho).')
    parser.add_argument('--h-ext-n-T', type=float, default=0.5,
                        help='Temperature power-law exponent for the envelope (default: 0.5). '
                             'T(r) = T_base * (r/r_base)^(-n_T).')
    parser.add_argument('--rho-min',      type=float, default=0.0,
                        help='Trim outer KEPLER zones with density below this threshold (g/cm^3) '
                             'before remapping onto the STELLA grid. Removes the extended '
                             'low-density atmosphere that causes STELLA radiation-transport '
                             'instabilities (BADSTE / KFLAG=-2) when the photosphere reaches '
                             'the very low density outer region. Default: 0 = no trimming.')
    parser.add_argument('--wind-zones',   type=int,   default=0,
                        help='Number of STELLA zones for CSM wind (default: 0 = no wind). '
                             'Wind zones are appended beyond the ejecta with rho ~ r^-2 profile. '
                             'Prevents radiation solver stall when photosphere exits ejecta.')
    parser.add_argument('--wind-mdot',    type=float, default=0.0,
                        help='Wind mass-loss rate in Msun/yr (default: 0 = no wind). '
                             'Typical WR wind: 1e-5; dense CSM: 1e-4 to 1e-3.')
    parser.add_argument('--wind-velocity', type=float, default=100.0,
                        help='Wind velocity in km/s (default: 100). '
                             'WR wind: ~1000-3000; RSG wind: ~10-50.')
    parser.add_argument('--wind-duration', type=float, default=1000.0,
                        help='Wind duration in years (default: 1000). Determines outer radius '
                             'of wind: r_outer = v_wind * t_wind.')
    parser.add_argument('--wind-T',       type=float, default=5000.0,
                        help='Wind temperature in K (default: 5000).')
    parser.add_argument('--wind-density-profile', type=str, default='matched',
                        help="Density profile: 'accelerating' (beta-law v(r) wind, "
                             "rho=max(ejecta-slope extrapolation, Mdot/4pir^2v) — "
                             "continuous at edge), "
                             "'crossover' (upper envelope of ejecta-slope "
                             "extrapolation and the Mdot-set 1/r² wind — the "
                             "profiles join directly at their crossing radius; "
                             "replaces the deprecated 'sigmoid'), "
                             "'matched' (default, 1/r² wind "
                             "starting at ejecta surface density), 'wind' (1/r² "
                             "with density from Mdot), 'uniform', or 'extrapolate' "
                             "(power-law extension matching the ejecta's own "
                             "local density slope -- no gradient kink; does not "
                             "require --wind-mdot).")
    parser.add_argument('--wind-power-law-index', type=float, default=None,
                        help='Override the fitted ejecta surface slope for '
                             'crossover/extrapolate profiles. Steeper index '
                             '(e.g. 12-15) reduces CSM mass and minimizes '
                             'early-time interaction bump. '
                             'Default: None (use fitted slope ~6.5).')
    parser.add_argument('--wind-rho-floor', type=float, default=1.0e-20,
                        help='Minimum wind density in g/cc (default: 1e-20). '
                             'Steep extrapolated tails can underflow below '
                             'STELLA table limits; zones below the floor are '
                             'clamped. Set 0 to disable.')
    parser.add_argument('--reduce-surf-zones', type=int, default=0,
                        help='Number of outermost ejecta zones to reduce density (default: 0). '
                             'Creates a smoother surface profile to help the radiation solver.')
    parser.add_argument('--reduce-surf-rho', type=float, default=0.0,
                        help='Target density at outermost zone in g/cc (default: 0 = disabled). '
                             'Example: 1e-6 gives a gradual ramp from ejecta density to 1e-6.')
    parser.add_argument('--smooth-profile', action='store_true', default=False,
                        help='Smooth steepest density/temperature gradients to help STELLA solver '
                             'converge during the thermal bomb phase. Only smooths zones where '
                             '|dlog rho/dlog r| > max_log_gradient. Default: False.')
    parser.add_argument('--max-log-gradient', type=float, default=5.0,
                        help='Threshold |dlog rho/dlog r| above which smoothing is applied (default: 5.0). '
                             'KEPLER pre-SN models can have gradients of -8 or steeper at the He core edge, '
                             'which causes STELLA GEAR+BGH solver to stall (BADSTE / KFLAG=-2).')
    parser.add_argument('--smooth-sigma', type=float, default=3.0,
                        help='Gaussian smoothing width in zone units (default: 3.0). '
                             'Larger values smooth more broadly; smaller values preserve more detail.')
    parser.add_argument('--comp-smooth-sigma', type=float, default=10.0,
                        help='Composition-interface Gaussian smoothing width in zone units '
                             '(default: 10.0). Applied unconditionally (independent of '
                             '--smooth-profile) to prevent sharp opacity jumps that destabilise '
                             'STELLA radiation transport. Lower values (e.g. 3-5) preserve '
                             'sharper shell boundaries (e.g. the O-shell) at the cost of less '
                             'aggressive smoothing of numerically problematic interfaces.')
    parser.add_argument('--comp-smooth-max-jump', type=float, default=0.02,
                        help='Composition per-zone jump threshold |dX/dzone| above which Gaussian '
                             'smoothing is applied (default: 0.02). Higher values leave more '
                             'shell interfaces untouched.')
    parser.add_argument('--comp-smooth-max-dx', type=float, default=0.08,
                        help='Maximum per-zone composition change enforced by the derivative-'
                             'capping pass after Gaussian smoothing (default: 0.08). This hard '
                             'cap is what erases sharp shell boundaries (e.g. O->Si->C) that the '
                             'Gaussian smoothing alone does not fully round off; raise it (e.g. '
                             '0.15-0.3) to preserve more shell structure, or set very high '
                             '(e.g. 1.0) to effectively disable the capping pass.')
    parser.add_argument('--outdir',       type=str,   default='.',
                        help='Output directory for all output files (default: .)')
    parser.add_argument('--plot-dir',     type=str,   default=None,
                        help='Directory for composition smoothing plot (default: same as --outdir)')
    parser.add_argument('--write-dat',    action='store_true', default=True,
                        help='Also write .dat STELLA control file (default: True)')
    parser.add_argument('--write-input-profile', action='store_true', default=True,
                        help='Write [DEPRECATED] HERACLES-format input_profile for opacity tables (default: True). '
                             'Project uses STELLA -> CMFGEN pipeline.')
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # 1. Read Kepler file
    print(f"\n=== Reading Kepler presn file: {args.presn_file} ===")
    kep = read_kepler_presn(args.presn_file)

    # 2. Remap onto STELLA grid
    n_wind = args.wind_zones if (args.wind_mdot > 0 or args.wind_density_profile
                               in ('extrapolate', 'crossover', 'sigmoid')) else 0
    n_ejecta = args.nzon - n_wind
    print(f"\n=== Remapping onto {args.nzon} STELLA zones ({n_ejecta} ejecta + {n_wind} wind) ===")
    g = remap_to_stella_grid(kep, args.mass_cut, args.nzon,
                             t_floor_K=args.t_floor,
                             h_surface_fraction=args.h_surface_fraction,
                             h_surface_mass_msun=args.h_surface_mass_msun,
                             h_ext_fraction=args.h_ext_fraction,
                             h_ext_mass_msun=args.h_ext_mass_msun,
                             h_ext_r_outer_rsun=args.h_ext_r_outer_rsun,
                             h_ext_n_zones=args.h_ext_n_zones,
                             h_ext_T_surface_K=args.h_ext_T_surface_K,
                             h_ext_T_base_K=args.h_ext_T_base_K,
                             h_ext_n_rho=args.h_ext_n_rho,
                             h_ext_n_T=args.h_ext_n_T,
                             rho_min=args.rho_min,
                             n_wind_zones=n_wind,
                             wind_mdot=args.wind_mdot,
                             wind_velocity=args.wind_velocity,
                             wind_duration_yr=args.wind_duration,
                             wind_T_K=args.wind_T,
                             wind_density_profile=args.wind_density_profile,
                             wind_power_law_index=args.wind_power_law_index,
                             wind_rho_floor=args.wind_rho_floor,
                             reduce_surf_n_zones=args.reduce_surf_zones,
                             reduce_surf_rho_target=args.reduce_surf_rho,
                             initial_velocity_energy=args.initial_velocity_energy)

    # 2b. Smooth steepest gradients if requested
    if args.smooth_profile:
        # Save pre-smoothing profile for diagnostic comparison plots
        pre_smooth_path = outdir / f"{args.model_name}.hyd.pre_smooth"
        print(f"\n=== Saving pre-smoothing profile: {pre_smooth_path} ===")
        write_hyd(str(pre_smooth_path), g, args.mass_cut)

        print(f"\n=== Smoothing steep gradients ===")
        g = smooth_steep_gradients(g, max_log_gradient=args.max_log_gradient,
                                   sigma_zones=args.smooth_sigma, verbose=True)

    # 2b. Smooth composition interfaces to prevent sharp opacity jumps
    #     that trigger radiation transport instabilities in STELLA
    print(f"\n=== Smoothing composition interfaces ===")
    # Save pre-smoothing composition for before/after plot
    import copy as _copy
    g_pre_comp = _copy.deepcopy(g)
    g = smooth_composition(g, sigma_zones=args.comp_smooth_sigma,
                           max_fraction_jump=args.comp_smooth_max_jump,
                           max_dX=args.comp_smooth_max_dx,
                           verbose=True)

    # Generate before/after composition smoothing plot
    plot_dir = Path(args.plot_dir) if args.plot_dir else outdir
    plot_dir.mkdir(parents=True, exist_ok=True)
    plot_path = plot_dir / f"{args.model_name}_composition_smoothing.png"
    try:
        plot_composition_smoothing(g_pre_comp, g, args.model_name, plot_path)
    except Exception as e:
        print(f"  WARNING: composition smoothing plot failed: {e}")

    # 3. Write .hyd
    hyd_path = outdir / f"{args.model_name}.hyd"
    print(f"\n=== Writing {hyd_path} ===")
    write_hyd(str(hyd_path), g, args.mass_cut)

    # 4. Write .abn (Ni56 goes in the last column -- eve2 reads it from there)
    abn_path = outdir / f"{args.model_name}.abn"
    print(f"\n=== Writing {abn_path} ===")
    write_abn(str(abn_path), g,
              ni56_mass_msun=args.ni56_mass, ni56_dist_msun=args.ni56_dist)

    # 5. Write .xni (Ni56 distribution)
    xni_path = outdir / f"{args.model_name}.xni"
    print(f"\n=== Writing {xni_path} ===")
    write_xni(str(xni_path), g, args.ni56_mass, args.ni56_dist, args.nzon)

    # 6. Write .dat
    dat_path = outdir / f"{args.model_name}.dat"
    print(f"\n=== Writing {dat_path} ===")
    custom_times = ([float(x) for x in args.custom_output_times.split(',')]
                    if args.custom_output_times else None)
    write_dat(str(dat_path), args.model_name, args.target_days,
              explosion_energy=args.explosion_energy,
              explosion_mass=args.explosion_mass,
              explosion_time=args.explosion_time,
              explosion_delay=args.explosion_delay,
              initial_kinetic_energy=args.initial_kinetic_energy,
              nstmax=args.nstmax,
              eps=args.eps, hmin=args.hmin, hmax=args.hmax, meth=args.meth,
              iout=args.iout, nout=args.nout,
              fit_tau=args.fit_tau, tau_tol=args.tau_tol,
              aq=args.aq, bq=args.bq,
              enable_convection=args.enable_convection,
              disable_scattering=args.disable_scattering,
              givdtl=args.givdtl,
              wacc_r=args.wacc_r, wacc_v=args.wacc_v,
              wacc_t=args.wacc_t, wacc_rad=args.wacc_rad,
              temp_floor=args.temp_floor, floor_rad=args.floor_rad,
              snapshot_interval=args.snapshot_interval,
              ni56_mass=args.ni56_mass, ni56_dist=args.ni56_dist,
              custom_output_times=custom_times)

    # 7. Write [DEPRECATED] HERACLES-format input_profile (needed by opacity table generator)
    # Note: HERACLES is deprecated; project uses STELLA -> CMFGEN pipeline
    ip_path = outdir / "input_profile"
    print(f"\n=== Writing {ip_path} ===")
    write_input_profile(str(ip_path), g)

    print(f"\n=== Done. Files written to {outdir}/ ===")
    print(f"  {hyd_path.name}, {abn_path.name}, {xni_path.name}, {dat_path.name}, input_profile")


if __name__ == '__main__':
    main()
