#!/usr/bin/env python3
"""
Prepare the SN_B initial model: H/He shell + connecting wind + old SN ejecta.

This script creates the initial conditions for the SN_B model, which represents
a reverse-shock (RSH) H/He shell launched at 5000 km/s into a dense wind,
which connects to a 400-day-old metal-rich supernova ejecta (the "CSM").

The SN_B model consists of three components (inner to outer):
  1. Inner zones: A new H/He shell (1.0 Msun, X=0.7, Y=0.3) at 5000 km/s
  2. Middle zones: A 1/r^2 wind (Mdot=1e-3 Msun/yr, v_wind=600 km/s)
     connecting the shell to the old ejecta
  3. Outer zones: The 400-day homologously evolved SN_12C_CSM2 ejecta
     (metal-rich, the "CSM" for the interaction)

STELLA starts at t=0 (the shell launch time). The old SN ejecta is already
400 days old at this point — its structure represents the pre-existing CSM.

Explosion convention (two-event workflow):
  This script builds the SN_B / embedded-shell setup, which is
  VELOCITY-DRIVEN.  The dynamics are carried entirely by the constructed
  initial conditions — the shell's mass and launch velocity plus the
  coasting old SN_A ejecta.  The .dat is therefore written with
  EBurst=0 (no thermal bomb may re-explode the structure) and EKO=0
  (EKO>0 would overwrite the .hyd velocity column with STELLA's
  triangular profile, destroying the shell/wind/ejecta velocity field).
  --explosion-energy exists only as an escape hatch and must stay 0.

  The SN_A setup is the opposite convention: a KEPLER pre-SN model is
  blown up by the EBurst thermal bomb — see
  convert_kepler_presn_to_stella.py.

**Recommended configuration for a solar-composition shell + wind + old SN
ejecta with maximum STELLA stability:**

    python prepare_SN_B_model.py <source> <output> \\
        --shell-composition solar --density-matched-wind

Using ``--shell-composition solar`` together with ``--density-matched-wind``
is the most stable setup for this scenario:

  * ``--shell-composition solar`` gives the shell realistic metal content
    (O, Si, Ca, Fe) matching the wind and ejecta, avoiding composition
    discontinuities that can destabilise the STELLA opacity calculation.
  * ``--density-matched-wind`` derives the wind Mdot so that
    ``rho_wind(r_interface) = rho_ejecta(r_interface)``, eliminating the
    density cliff at the wind-ejecta interface.  This is much more stable
    for STELLA than a standard 1/r^2 wind, which can produce a >10 dex
    density jump at the contact discontinuity.

Outputs:
  - {model}.hyd: STELLA hydrodynamics input
  - {model}.abn: STELLA composition input
  - {model}.xni: Ni56 distribution (zero — no radioactive decay in SN_B)
  - {model}.dat: STELLA control file
  - input_profile: 12-column profile for convert_profile_to_stella.py

Usage:
    python prepare_SN_B_model.py <source_dir> <output_dir> [options]

    --shell-mass 1.0        H/He shell mass in Msun
    --shell-velocity 5000   Shell launch velocity in km/s
    --shell-xh 0.7          H mass fraction in shell
    --shell-yhe 0.3         He mass fraction in shell
    --wind-mdot 1e-3        Wind mass-loss rate in Msun/yr
    --wind-velocity 600     Wind velocity in km/s
    --wind-T 5000           Wind temperature in K
    --target-days 45        STELLA simulation duration (days)
    --snapshot-day 30       CMFGEN snapshot day after shell launch
    --nzon 500              Number of STELLA zones
    --source-age 400        Age of source ejecta at shell launch (days)
"""

import argparse
import sys
import struct
import numpy as np
from pathlib import Path

# Import shared utilities
sys.path.insert(0, str(Path(__file__).parent))
from convert_stella_swd_to_cmfgen import read_swd, normalize_swd_composition
from homologous_evolution import find_closest_epoch, homologous_evolve

MSUN = 1.989e33
K_B = 1.3807e-16
M_H = 1.6735e-24
A_RAD = 7.5657e-15
KAPPA_ES = 0.2 * (1 + 0.7) / (1 + 0.7 + 0.25 * 0.3)  # electron scattering opacity
YR_TO_S = 365.25 * 24 * 3600  # seconds per year

# Solar composition (Asplund et al. 2009, Z=0.0134).
# Only the 6 elements needed for CMFGEN spectral modelling: H, He, O, Si, Ca, Fe.
# Other solar metals (C, N, Ne, Na, Mg, Al, S, Ar, Ni) are set to zero.
# Mass fractions are renormalized to sum=1 after combining (see below).
SOLAR_COMPOSITION = {
    'fH':  0.7381,
    'fHe': 0.2485,
    'fO':  0.00573,
    'fSi': 0.00067,
    'fCa': 0.00006,
    'fFe': 0.00129,
}


def smooth_ejecta_profile(ejecta, window=15, polyorder=3, comp_window=7):
    """Smooth numerical noise from a STELLA .swd ejecta profile.

    At late epochs (hundreds of days), the STELLA radiation-hydrodynamics
    evolution leaves zone-to-zone scatter in density, temperature, and
    velocity from:
      - adaptive grid rezoning
      - radiation transport numerical noise
      - artificial viscosity dissipation patterns
      - composition interface oscillations

    A 400-day-old ejecta should be in smooth homologous expansion with
    thermally equilibrated radiation transport.  Zone-to-zone scatter is
    therefore numerical artefact, not physical structure.

    This function applies a Savitzky-Golay filter (which preserves gradients
    and peaks better than a simple moving average) to:
      - log10(rho):  removes density spikes that cause opacity jumps
      - log10(T):    removes temperature oscillations that cause ionisation
                     fronts to jump between zones
      - v:           removes velocity noise that causes spurious shocks

    The filter is applied in log-space for rho and T (where the profiles
    span many orders of magnitude) and in linear space for v.

    Parameters:
        ejecta:    dict from read_source_ejecta() with r, rho, v, T, composition
        window:    Savitzky-Golay window size (odd, default 15).
                   Larger = more smoothing.  Must be > polyorder.
        polyorder: Savitzky-Golay polynomial order (default 3).
                   Lower = more smoothing.  3 preserves power-law gradients.

    Returns:
        (smoothed_ejecta, raw_copy) tuple:
          smoothed_ejecta: dict with smoothed rho, v, T (same structure)
          raw_copy:        dict with the original unsmoothed values (for
                           before/after plotting)
    """
    from scipy.signal import savgol_filter

    # Make sure window is odd and valid
    n = len(ejecta['r'])
    window = min(window, n if n % 2 == 1 else n - 1)
    if window <= polyorder:
        window = polyorder + 2
        if window % 2 == 0:
            window += 1
    if window > n:
        window = n if n % 2 == 1 else n - 1

    # Save raw copies for before/after comparison
    raw = {
        'r': ejecta['r'].copy(),
        'rho': ejecta['rho'].copy(),
        'v': ejecta['v'].copy(),
        'T': ejecta['T'].copy(),
    }

    # Smooth in log-space for rho and T (span many orders of magnitude)
    log_rho = np.log10(np.maximum(ejecta['rho'], 1e-300))
    log_T = np.log10(np.maximum(ejecta['T'], 1e-300))

    log_rho_smooth = savgol_filter(log_rho, window, polyorder)
    log_T_smooth = savgol_filter(log_T, window, polyorder)

    ejecta['rho'] = 10.0 ** log_rho_smooth
    ejecta['T'] = 10.0 ** log_T_smooth

    # Smooth velocity in linear space (preserves homologous gradient)
    ejecta['v'] = savgol_filter(ejecta['v'], window, polyorder)

    # Diagnostics: measure noise reduction
    def noise_metric(arr):
        """Zone-to-zone scatter as fraction of local value."""
        if len(arr) < 3:
            return 0.0
        d = np.diff(np.log10(np.maximum(arr, 1e-300)))
        return np.std(d)

    raw_noise_rho = noise_metric(raw['rho'])
    raw_noise_T = noise_metric(raw['T'])
    raw_noise_v = noise_metric(np.maximum(raw['v'], 1e-300))
    smooth_noise_rho = noise_metric(ejecta['rho'])
    smooth_noise_T = noise_metric(ejecta['T'])
    smooth_noise_v = noise_metric(np.maximum(ejecta['v'], 1e-300))

    print(f"  Ejecta smoothing (Savitzky-Golay, window={window}, order={polyorder}):")
    print(f"    Density noise:  {raw_noise_rho:.4f} -> {smooth_noise_rho:.4f} "
          f"({(1 - smooth_noise_rho/raw_noise_rho)*100:.0f}% reduction)")
    print(f"    Temperature noise: {raw_noise_T:.4f} -> {smooth_noise_T:.4f} "
          f"({(1 - smooth_noise_T/raw_noise_T)*100:.0f}% reduction)")
    print(f"    Velocity noise: {raw_noise_v:.4f} -> {smooth_noise_v:.4f} "
          f"({(1 - smooth_noise_v/raw_noise_v)*100:.0f}% reduction)")

    # --- Soften sharp composition interfaces ---
    # The source ejecta has composition steps (e.g. metal-core -> He shell
    # and He -> H envelope boundaries, and the inner edge adjoining the
    # wind) that can span only a few zones.  These produce Rosseland mean
    # opacity discontinuities that can stall STELLA's implicit solver.
    # A short moving average widens each step over ~comp_window zones
    # without erasing the compositional layering.
    comp_keys = ['fH', 'fHe', 'fC', 'fO', 'fSi', 'fCa', 'fFe', 'fNi']
    present = [k for k in comp_keys
               if k in ejecta and len(ejecta[k]) == n]
    if comp_window > 1 and present:
        kernel = np.ones(comp_window) / comp_window
        he = comp_window // 2
        for k in present:
            c = np.asarray(ejecta[k], dtype=float)
            sm = np.convolve(c, kernel, mode='same')
            # Preserve edges (convolve tapers at boundaries)
            sm[:he] = c[:he]
            sm[-he:] = c[-he:]
            ejecta[k] = np.clip(sm, 0.0, None)
        # Renormalize so mass fractions sum to 1 in every zone
        tot = np.zeros(n)
        for k in present:
            tot += ejecta[k]
        tot = np.maximum(tot, 1e-10)
        for k in present:
            ejecta[k] /= tot
        print(f"  Composition smoothing: window={comp_window} zones "
              f"({len(present)} species, renormalized)")

    # --- Drop detached outer zones ---
    # STELLA can leave isolated grid points far beyond the ejecta surface:
    # zones separated from the next zone inward by a radius gap much larger
    # than the local grid spacing, with artifact density/velocity (pile-up
    # deceleration, near-floor density).  These must be removed BEFORE the
    # outer-boundary cleanup, otherwise the tail extrapolation is anchored
    # to a bad point.
    r_arr = np.asarray(ejecta['r'])
    n_r = len(r_arr)
    desc = r_arr[0] > r_arr[-1]
    r_asc = r_arr[::-1] if desc else r_arr
    n_keep = n_r
    while n_keep > 3 and r_asc[n_keep - 1] / r_asc[n_keep - 2] > 1.15:
        n_keep -= 1
    if n_keep < n_r:
        n_drop = n_r - n_keep
        for k in list(ejecta.keys()):
            a = np.asarray(ejecta[k])
            if a.ndim == 1 and a.shape[0] == n_r:
                a_asc = a[::-1] if desc else a
                ejecta[k] = (a_asc[:n_keep])[::-1] if desc else a_asc[:n_keep]
        print(f"  Dropped {n_drop} detached outer zone(s) "
              f"(r > {r_asc[n_keep-1]:.3e} cm; gap ratio > 1.15)")
        n = n_keep

    # --- Clean outer boundary ---
    # STELLA's Lagrangian grid often produces a "pile-up" at the outer
    # boundary: the outermost zone is decelerated (velocity drops), causing
    # an artificial density enhancement.  The density also becomes
    # non-monotonic in the last few zones from numerical noise.
    # This creates a spurious feature that can reflect radiation and stall
    # the solver.  n_edge=40 covers the full noisy tail (dip + bump +
    # pile-up), not just the last few zones.
    ejecta = clean_ejecta_outer_boundary(ejecta, n_edge=40)

    return ejecta, raw


def clean_ejecta_outer_boundary(ejecta, n_edge=20):
    """Remove STELLA outer-boundary pile-up artifacts.

    STELLA's Lagrangian grid produces two known artifacts at the outer
    boundary:

    1. **Velocity pile-up**: the outermost zone(s) are decelerated relative
       to the smooth homologous profile, creating an artificial density
       enhancement.  In extreme cases the velocity actually *decreases* at
       the last zone.

    2. **Density non-monotonicity**: the density in the last few zones
       oscillates (goes up and down) instead of monotonically decreasing,
       due to numerical noise in the radiation-hydrodynamics evolution.

    This function fixes both by:
      - Extrapolating the velocity from the clean inner ejecta (where the
        profile is smooth) to overwrite the noisy outer edge.
      - Enforcing monotonically decreasing density in the outer zones.
      - Tapering the last few zones' density to smoothly approach the
        minimum, removing any spikes.

    Parameters:
        ejecta:  dict with r, rho, v, T arrays
        n_edge:  number of outer zones to clean (default 20)

    Returns:
        ejecta dict with cleaned outer boundary
    """
    r_in = np.asarray(ejecta['r'])
    n = len(r_in)

    if n < n_edge * 2:
        return ejecta

    # The source ejecta is stored outer->inner (descending r).  Normalize
    # to ascending (inner->outer) so the outer boundary is the LAST n_edge
    # zones, then write back in the original order at the end.
    descending = r_in[0] > r_in[-1]
    if descending:
        r = r_in[::-1]
        rho = np.asarray(ejecta['rho'])[::-1].copy()
        v = np.asarray(ejecta['v'])[::-1].copy()
    else:
        r = r_in
        rho = ejecta['rho']
        v = ejecta['v']

    # --- 1. Fix velocity pile-up ---
    # Use linear extrapolation in linear space from the clean inner ejecta.
    # The velocity profile is nearly linear (v ~ r/t for homologous
    # expansion), so a linear fit is stable and physical.  A quadratic fit
    # in log-log space can diverge wildly at the outer edge.
    clean_slice = slice(0, n - n_edge)
    r_clean = r[clean_slice]
    v_clean = v[clean_slice]
    coeffs = np.polyfit(r_clean, v_clean, 1)  # v = a*r + b
    v_edge_fit = np.polyval(coeffs, r[n - n_edge:])

    # Blend: smoothstep transition over the edge region (0 = actual, 1 = fit)
    blend = np.linspace(0, 1, n_edge)
    # Smoothstep (3u^2-2u^3): blends actual -> fitted with zero slope at both
    # ends of the edge region, avoiding a kink where the extrapolation starts.
    blend = blend * blend * (3.0 - 2.0 * blend)
    v[n - n_edge:] = (1 - blend) * v[n - n_edge:] + blend * v_edge_fit

    # Ensure velocity is monotonically increasing outward (v increases
    # with r in homologous expansion).  Walk inward->outward so each
    # correction propagates to the next comparison.
    dv_avg = np.mean(np.diff(v_clean)) if len(v_clean) > 1 else 0.0
    dv_avg = max(dv_avg, 0.0)
    for i in range(n - n_edge, n):
        if v[i] < v[i - 1]:
            v[i] = v[i - 1] + dv_avg

    # --- 2. Fix density non-monotonicity ---
    # Enforce monotonically decreasing density in the outer ejecta.
    # Where density increases, replace with the previous zone's density
    # times the local (clean-region) ratio, capped below 1.
    n_fixed = 0
    for i in range(n - n_edge, n - 1):
        if rho[i + 1] > rho[i]:
            # Non-monotonic: cap at previous zone's density * smooth ratio
            # Use the local density gradient from the clean region
            if i > 0:
                ratio = rho[i] / max(rho[i - 1], 1e-300)
                rho[i + 1] = rho[i] * min(ratio, 0.99)
            else:
                rho[i + 1] = rho[i] * 0.9
            n_fixed += 1

    # Taper the last 5 zones to smoothly approach the minimum density
    # This removes any remaining spikes
    n_taper = min(5, n_edge)
    for i in range(n - n_taper, n):
        if i > 0:
            # Force smooth decrease using the average gradient from the
            # clean region
            clean_ratios = rho[clean_slice][1:] / rho[clean_slice][:-1]
            avg_ratio = np.median(clean_ratios[-20:])
            rho[i] = rho[i - 1] * avg_ratio

    if n_fixed > 0:
        print(f"  Outer boundary cleanup: fixed {n_fixed} non-monotonic density "
              f"zones, tapered last {n_taper} zones")
        print(f"    Velocity: extrapolated from inner ejecta (pile-up removed)")

    # Write back in the original (descending) order if needed
    if descending:
        ejecta['rho'] = rho[::-1]
        ejecta['v'] = v[::-1]

    return ejecta


def plot_ejecta_smoothing(raw, smoothed, output_path, source_age=400):
    """Generate a before/after smoothing comparison plot.

    Shows the raw (noisy) and smoothed ejecta profiles for rho, T, and v
    in a 3-panel figure with log-log axes.

    Parameters:
        raw:       dict with r, rho, v, T (before smoothing)
        smoothed:  dict with r, rho, v, T (after smoothing)
        output_path: path to save the plot
        source_age: source age in days (for title)
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(3, 1, figsize=(14, 14), sharex=True)

    r_raw = raw['r']
    r_sm = smoothed['r']

    # --- Density ---
    ax = axes[0]
    ax.plot(r_raw, raw['rho'], 'gray', linewidth=0.5, alpha=0.7, label='Raw (noisy)')
    ax.plot(r_sm, smoothed['rho'], 'red', linewidth=1.5, label='Smoothed')
    ax.set_ylabel('Density (g/cm$^3$)', fontsize=12)
    ax.set_title(f'SN_A Ejecta Smoothing (Savitzky-Golay) — {source_age}d source', fontsize=14)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(fontsize=10, loc='best')
    ax.grid(True, alpha=0.3)

    # --- Temperature ---
    ax = axes[1]
    ax.plot(r_raw, raw['T'], 'gray', linewidth=0.5, alpha=0.7, label='Raw (noisy)')
    ax.plot(r_sm, smoothed['T'], 'blue', linewidth=1.5, label='Smoothed')
    ax.set_ylabel('Temperature (K)', fontsize=12)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(fontsize=10, loc='best')
    ax.grid(True, alpha=0.3)

    # --- Velocity ---
    ax = axes[2]
    ax.plot(r_raw, raw['v'], 'gray', linewidth=0.5, alpha=0.7, label='Raw (noisy)')
    ax.plot(r_sm, smoothed['v'], 'green', linewidth=1.5, label='Smoothed')
    ax.set_ylabel('Velocity (km/s)', fontsize=12)
    ax.set_xlabel('Radius (cm)', fontsize=12)
    ax.set_xscale('log')
    ax.legend(fontsize=10, loc='best')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Ejecta smoothing plot saved: {output_path}")


def read_source_ejecta(source_dir, source_age=400, tolerance=100, t_start_days=0.0,
                       force_homologous_velocity=False,
                       smooth_ejecta=True, smooth_window=15, smooth_polyorder=3,
                       comp_window=7):
    """Read the source SN ejecta profile at the correct age for STELLA start.

    Reads the STELLA .swd file from the source model, finds the latest epoch,
    and homologously evolves it to (source_age + t_start_days) days.

    If smooth_ejecta=True (default), applies a Savitzky-Golay filter to remove
    numerical noise from the STELLA RHD evolution.  At hundreds of days, the
    ejecta should be in smooth homologous expansion; zone-to-zone scatter is
    numerical artefact.  The raw (pre-smoothing) profile is saved in the
    returned dict as 'raw_ejecta' for before/after plotting.

    The source_age is the age of the SN_A ejecta at the time of SN_B explosion
    (t=0 for SN_B). If STELLA starts at t_start_days > 0 (pre-expanded shell),
    the SN_A ejecta has continued to expand during those t_start_days, so its
    age at the STELLA initial condition is source_age + t_start_days.

    If force_homologous_velocity=True, the ejecta velocity is overwritten with
    v = r / t (pure Hubble flow) after homologous evolution.  This removes any
    residual non-homologous structure from the source STELLA model (e.g. CSM
    deceleration, shock interactions) and ensures the ejecta is in exact free
    expansion.  This is physically appropriate for an SN_A ejecta that is
    hundreds of days old and should be in the free-expansion phase.

    Returns dict with r, rho, v, T, composition (outer->inner order).
    """
    source_dir = Path(source_dir)
    model_name = source_dir.name

    # Effective age: SN_A age at SN_B explosion + t_start_days (STELLA start offset)
    effective_age = source_age + t_start_days
    if t_start_days > 0:
        print(f"  SN_A ejecta age: {source_age}d at SN_B explosion + {t_start_days}d t_start = {effective_age}d at STELLA start")

    # Find .swd file
    swd_file = source_dir / 'stella' / 'strad' / f'{model_name}.swd'
    if not swd_file.exists():
        swd_file = source_dir / 'stella' / f'{model_name}.swd'
    if not swd_file.exists():
        raise FileNotFoundError(f"Source .swd not found: {swd_file}")

    print(f"  Reading source .swd: {swd_file}")
    epochs = read_swd(str(swd_file))
    available = sorted(epochs.keys())
    print(f"  Available epochs: {available[0]:.2f} ... {available[-1]:.2f} days "
          f"({len(available)} total)")

    # If effective_age is beyond the available range, use the latest epoch
    # and homologously evolve forward
    latest_epoch = available[-1]
    if effective_age > latest_epoch + tolerance:
        print(f"  Requested age {effective_age} days is beyond available data "
              f"(latest: {latest_epoch:.2f} days). Using latest epoch and "
              f"homologously evolving.")
        t_source = latest_epoch
        epoch_data = epochs[t_source]
    else:
        # Find closest epoch within tolerance
        t_source, epoch_data = find_closest_epoch(epochs, effective_age, tolerance)

    # Homologously evolve to effective_age (source_age + t_start_days)
    if abs(t_source - effective_age) > 1.0:
        print(f"  Homologously evolving from {t_source:.2f} to {effective_age} days...")
        epoch_data = homologous_evolve(epoch_data, t_source, effective_age)
    else:
        print(f"  Using epoch at {t_source:.2f} days (matches requested age)")

    # Force pure homologous velocity: v = r / t
    # This removes any residual non-homologous structure from the source model
    # (e.g. CSM deceleration, shock interactions) and ensures exact free
    # expansion: v(r) = r / (effective_age * 86400 s) in cm/s.
    if force_homologous_velocity:
        t_sec = effective_age * 86400.0
        v_old = epoch_data['v'].copy()
        v_new = epoch_data['r'] / t_sec  # cm/s
        # Convert to km/s for reporting (v in .swd is km/s)
        v_old_kms = v_old  # already km/s
        v_new_kms = v_new / 1e5
        ratio = np.where(v_old_kms > 0, v_new_kms / v_old_kms, 1.0)
        print(f"  Forcing homologous velocity: v = r / t (t={effective_age:.1f}d)")
        print(f"    v/ratio (new/old): min={ratio.min():.3f}, max={ratio.max():.3f}, "
              f"mean={ratio.mean():.3f}")
        print(f"    v_old: [{v_old_kms.min():.1f}, {v_old_kms.max():.1f}] km/s")
        print(f"    v_new: [{v_new_kms.min():.1f}, {v_new_kms.max():.1f}] km/s")
        epoch_data['v'] = v_new_kms

    # Normalize composition
    if epoch_data.get('composition') is not None:
        comp = normalize_swd_composition(epoch_data['composition'])
        epoch_data['fH'] = comp['H']
        epoch_data['fHe'] = comp['He']
        epoch_data['fC'] = comp['C']
        epoch_data['fO'] = comp['O']
        epoch_data['fSi'] = comp['Si']
        epoch_data['fCa'] = comp['Ca']
        epoch_data['fFe'] = comp['Fe']
        epoch_data['fNi'] = comp['Ni']
    else:
        for key in ['fH', 'fHe', 'fC', 'fO', 'fSi', 'fCa', 'fFe', 'fNi']:
            epoch_data[key] = np.zeros(len(epoch_data['r']))
        epoch_data['fHe'] = np.ones(len(epoch_data['r']))

    # --- Smooth numerical noise from the ejecta profile (default: on) ---
    # The STELLA .swd at late epochs has zone-to-zone scatter from the RHD
    # evolution.  At hundreds of days, the ejecta should be in smooth
    # homologous expansion; the scatter is numerical artefact.  Smoothing
    # removes density spikes, temperature oscillations, and velocity noise
    # that would cause STELLA solver stalls or spurious shocks.
    if smooth_ejecta:
        print(f"\n  Smoothing ejecta profile (removing STELLA RHD numerical noise)...")
        epoch_data, raw_ejecta = smooth_ejecta_profile(
            epoch_data, window=smooth_window, polyorder=smooth_polyorder,
            comp_window=comp_window)
        # Save raw and smoothed profiles for before/after plotting.
        # The smoothed copy is captured NOW (before truncation and unit
        # conversion in main()) so the plot shows a true apples-to-apples
        # comparison of the full 500-zone ejecta.
        epoch_data['raw_ejecta'] = raw_ejecta
        epoch_data['smoothed_ejecta'] = {
            'r': epoch_data['r'].copy(),
            'rho': epoch_data['rho'].copy(),
            'v': epoch_data['v'].copy(),
            'T': epoch_data['T'].copy(),
        }

    return epoch_data


def _shell_mass_factor(n_rho, frac_v0):
    """Envelope mass factor for the flat-core + power-law shell.

    M = (4pi/3) rho_0 r0^3 * factor, where the envelope contributes
    int_{r0}^{r_outer} 4pi r^2 rho_0 (r0/r)^n_rho dr. frac_v0 = r0/r_outer.
    """
    if n_rho > 3:
        return 1 + (3.0 / (n_rho - 3)) * (1 - frac_v0**(n_rho - 3))
    if n_rho < 3:
        return 1 + (3.0 / (3 - n_rho)) * (frac_v0**(n_rho - 3) - 1)
    # n_rho == 3: analytic limit, M_envelope = 4pi rho_0 r0^3 ln(1/frac_v0)
    return 1 + 3.0 * np.log(1.0 / frac_v0)


def create_hhe_shell(shell_mass_msun, shell_velocity_kms, xh, yhe, nzon_shell=100,
                     r_inner_cm=1e8, T_shell_K=1e5,
                     n_rho=8.0, v0_kms=3000.0, n_T=0.4, T_interface_K=5000.0,
                     T_max_K=65000.0, t_start_days=10.0, r_outer_override=None,
                     composition='hhe'):
    """Create the H/He shell profile using Dessart+2015 (SN 2010jl) structure.

    The shell has a power-law density profile (Dessart et al. 2015, MNRAS 449, 4304):
      rho(r) = rho_0 * (V0 / V)^N_rho   for V > V0  (outer power-law region)
      rho(r) = rho_0                    for V <= V0  (inner constant-density core)

    where V = v_shell * r / r_outer is the homologous velocity, and V0 is the
    transition velocity from the flat core to the power-law envelope.

    The temperature follows:
      T(r) = T_interface * (rho(r) / rho_interface)^N_T   (capped at T_max)

    **Critical: t_start_days sets the shell radius via r_outer = v_shell * t_start.**
    Starting at t_start > 0 is essential because at t=0 the shell density is
    ~1e4 g/cc, giving a diffusion time of ~50 million days — STELLA cannot
    evolve this. At t_start=10 days, the shell has expanded to r~4e14 cm with
    rho~1e-11 g/cc, giving a diffusion time of ~430 days (manageable for a
    45-day simulation). This matches the Dessart 2010jl approach where the
    simulation starts at t=11.6 days with the shell pre-expanded.

    Parameters:
        n_rho:        Density power-law exponent (Dessart 2010jl: N_rho=8)
        v0_kms:       Core/envelope transition velocity (Dessart: V0~3000 km/s)
        n_T:          Temperature power-law exponent in density (Dessart: N_T=0.4)
        T_interface:  Temperature at shell outer edge (matches wind T)
        T_max:        Maximum temperature in inner core (Dessart: 65000 K)
        t_start_days: Time since shell launch at which STELLA starts (days).
                      Sets r_outer = v_shell * t_start. Must be > 0.
        r_outer_override: If provided, use this r_outer instead of v_shell*t_start.
                      Used for immediate_contact mode where shell outer edge
                      must match ejecta inner edge.
    """
    M_shell = shell_mass_msun * MSUN
    v_shell = shell_velocity_kms * 1e5  # cm/s
    v0 = v0_kms * 1e5  # cm/s

    if v0_kms >= shell_velocity_kms:
        raise ValueError(
            f"v0_kms ({v0_kms}) must be < shell_velocity_kms "
            f"({shell_velocity_kms}): the flat core must end below the "
            f"outer edge, otherwise the density profile inverts.")

    # Shell outer radius from launch time (NOT from density)
    t_start_s = t_start_days * 86400.0  # seconds
    if r_outer_override is not None:
        r_outer = r_outer_override
        print(f"  H/He shell: r_outer OVERRIDE = {r_outer:.3e} cm (immediate contact)")
    else:
        r_outer = v_shell * t_start_s

    # Core radius from V0 fraction
    frac_v0 = v0_kms / shell_velocity_kms  # fraction of v_shell where core ends
    r0 = frac_v0 * r_outer

    # Inner radius: set so the innermost zone is a small fraction of r0
    # (the core extends from r_inner to r0 at constant density)
    r_inner = max(r_inner_cm, r0 * 0.01)

    # Compute rho_0 from mass conservation
    factor = _shell_mass_factor(n_rho, frac_v0)
    rho_0 = M_shell / ((4/3) * np.pi * r0**3 * factor)

    # Density at outer edge
    rho_edge = rho_0 * frac_v0**n_rho

    # Diffusion time diagnostic: use the shell mean density, not rho_0.
    # rho_0 * r_outer^2 overestimates t_diff by ~ (r_outer/r0)^2 because most
    # of the shell mass sits in the core at r0 while the envelope is thin.
    kappa_est = 0.34
    c_light = 3e10
    rho_mean = M_shell / ((4.0 / 3.0) * np.pi * r_outer**3)
    t_diff = kappa_est * rho_mean * r_outer**2 / c_light

    print(f"  H/He shell (Dessart 2010jl structure, t_start={t_start_days}d):")
    print(f"    M={shell_mass_msun} Msun, v_outer={shell_velocity_kms} km/s")
    print(f"    N_rho={n_rho}, V0={v0_kms} km/s (core/envelope transition)")
    print(f"    r_inner={r_inner:.3e} cm, r0={r0:.3e} cm, r_outer={r_outer:.3e} cm")
    print(f"    rho_0={rho_0:.3e} g/cc (core), rho_edge={rho_edge:.3e} g/cc (outer)")
    print(f"    Density contrast (core/edge) = {rho_0/rho_edge:.1f}")
    print(f"    Diffusion time (mean density): {t_diff/86400:.1f} days; "
          f">> sim length means shell internal energy stays trapped and "
          f"the lightcurve is powered by interaction, not stored heat")
    print(f"    X_H={xh}, Y_He={yhe}, composition={composition}")

    # Log-spaced zones
    r = np.logspace(np.log10(r_inner), np.log10(r_outer), nzon_shell)

    # Power-law density: constant core + power-law envelope
    rho_arr = np.where(r <= r0, rho_0, rho_0 * (r0 / r)**n_rho)

    # Homologous velocity: v = v_shell * r / r_outer
    v_arr = v_shell * r / r_outer

    # Dessart-style temperature: T ∝ rho^N_T, from T_interface at edge to T_max in core
    rho_interface = rho_edge  # density at outer edge
    T_arr = T_interface_K * (rho_arr / rho_interface)**n_T
    T_arr = np.minimum(T_arr, T_max_K)
    T_arr = np.maximum(T_arr, T_interface_K)

    print(f"    T range: [{T_arr.min():.0f}, {T_arr.max():.0f}] K")

    # Composition: only H, He, O, Si, Ca, Fe (solar) or H/He-only
    if composition == 'solar':
        sc = SOLAR_COMPOSITION
        fH = np.full(nzon_shell, sc['fH'])
        fHe = np.full(nzon_shell, sc['fHe'])
        fO = np.full(nzon_shell, sc['fO'])
        fSi = np.full(nzon_shell, sc['fSi'])
        fCa = np.full(nzon_shell, sc['fCa'])
        fFe = np.full(nzon_shell, sc['fFe'])
        print(f"    Solar composition (Asplund+2009, 6 elements): Z={1.0-sc['fH']-sc['fHe']:.4f}")
    else:
        fH = np.full(nzon_shell, xh)
        fHe = np.full(nzon_shell, yhe)
        fO = np.zeros(nzon_shell)
        fSi = np.zeros(nzon_shell)
        fCa = np.zeros(nzon_shell)
        fFe = np.zeros(nzon_shell)

    return {
        'r': r, 'rho': rho_arr, 'v': v_arr, 'T': T_arr,
        'fH': fH, 'fHe': fHe, 'fO': fO, 'fSi': fSi, 'fCa': fCa, 'fFe': fFe,
        'r_inner': r_inner, 'r_outer': r_outer,
        'rho_0': rho_0, 'rho_edge': rho_edge, 'r0': r0,
        'n_rho': n_rho, 'M_shell': M_shell,
        't_start_days': t_start_days,
    }


def create_wind(mdot_msun_yr, v_wind_kms, T_wind_K, r_inner_cm, r_outer_cm,
                nzon_wind=100, composition='solar'):
    """Create a 1/r^2 wind profile connecting the shell to the outer ejecta.

    The wind has:
    - Density: rho(r) = Mdot / (4 * pi * r^2 * v_wind)
    - Velocity: v_wind (constant)
    - Temperature: T_wind_K (isothermal)
    - Composition: solar (Asplund+2009, default) or H/He-only ('hhe')

    Boundary transitions to the shell and ejecta are handled by the
    combiner's log-linear density buffer zones (avoiding solver stalls per
    AGENTS.md lessons).
    """
    mdot = mdot_msun_yr * MSUN / YR_TO_S  # g/s
    v_wind = v_wind_kms * 1e5  # cm/s

    print(f"  Wind: Mdot={mdot_msun_yr} Msun/yr, v={v_wind_kms} km/s, "
          f"composition={composition}")
    print(f"    r_inner={r_inner_cm:.3e} cm, r_outer={r_outer_cm:.3e} cm")

    # Log-spaced zones
    r = np.logspace(np.log10(r_inner_cm), np.log10(r_outer_cm), nzon_wind)

    # 1/r^2 density profile
    rho = mdot / (4 * np.pi * r**2 * v_wind)

    # Constant velocity
    v_arr = np.full(nzon_wind, v_wind)

    # Isothermal temperature
    T_arr = np.full(nzon_wind, T_wind_K)

    # Wind composition: solar (H, He, O, Si, Ca, Fe) or H/He-only
    if composition == 'solar':
        sc = SOLAR_COMPOSITION
        fH = np.full(nzon_wind, sc['fH'])
        fHe = np.full(nzon_wind, sc['fHe'])
        fO = np.full(nzon_wind, sc['fO'])
        fSi = np.full(nzon_wind, sc['fSi'])
        fCa = np.full(nzon_wind, sc['fCa'])
        fFe = np.full(nzon_wind, sc['fFe'])
        print(f"    Solar composition (Asplund+2009, 6 elements): Z={1.0-sc['fH']-sc['fHe']:.4f}")
    else:
        fH = np.full(nzon_wind, 0.7)
        fHe = np.full(nzon_wind, 0.3)
        fO = np.zeros(nzon_wind)
        fSi = np.zeros(nzon_wind)
        fCa = np.zeros(nzon_wind)
        fFe = np.zeros(nzon_wind)

    # Compute wind mass
    r_lo = np.concatenate([[r_inner_cm * 0.99], r[:-1]])
    dm = (4/3) * np.pi * (r**3 - r_lo**3) * rho
    M_wind = np.sum(np.maximum(dm, 0))
    print(f"    Wind mass: {M_wind/MSUN:.4f} Msun")
    print(f"    rho range: [{rho.min():.3e}, {rho.max():.3e}] g/cc")
    t_flow = (r_outer_cm - r_inner_cm) / v_wind / 86400.0
    print(f"    Implied wind duration: {t_flow:.0f} d "
          f"((r_outer - r_inner) / v_wind)")

    # Compute optical depth through wind
    kappa_est = 0.34
    tau_wind = kappa_est * mdot / (4 * np.pi * v_wind * r_inner_cm)
    print(f"    Estimated optical depth (at r_inner): {tau_wind:.1f}")

    return {
        'r': r, 'rho': rho, 'v': v_arr, 'T': T_arr,
        'fH': fH, 'fHe': fHe, 'fO': fO, 'fSi': fSi, 'fCa': fCa, 'fFe': fFe,
        'r_inner': r_inner_cm, 'r_outer': r_outer_cm,
        'M_wind': M_wind,
    }


def _interface_ramp_weights(n, transition_zones=4):
    """Linear-ramp blend weights for n buffer zones.

    Returns array of length n with values:
      0.0 for inner zones (shell side)
      linear ramp over `transition_zones` centered at the midpoint
      1.0 for outer zones (ejecta side)

    This produces a sharp interface over ~transition_zones zones instead
    of a gradual blend across the entire buffer — the linear equivalent of
    the old midpoint sigmoid, without the S-curve tails.
    """
    if n <= 1:
        return np.array([1.0])
    if n <= transition_zones * 2:
        # Buffer too small for a centered transition — ramp across the whole
        # buffer so the last zone reaches the outer value.
        return np.linspace(0.0, 1.0, n)
    center = n / 2.0
    half = transition_zones / 2.0
    idx = np.arange(n, dtype=float)
    w = (idx - (center - half)) / max(transition_zones, 1)
    return np.clip(w, 0.0, 1.0)


def _interface_blend(n, val_inner, val_outer, transition_zones=4):
    """Blend two scalar values across n zones with a linear interface ramp.

    Inner zones get val_inner, outer zones get val_outer, with a linear
    ramp over ~transition_zones centered at the midpoint — the sharp
    interface between two physical regions, without sigmoid stitching.
    """
    w = _interface_ramp_weights(n, transition_zones)
    return val_inner * (1.0 - w) + val_outer * w


def clean_combined_outer_boundary(combined, n_edge=20, ej_start_idx=None):
    """Clean STELLA outer-boundary pile-up artifacts from the final combined grid.

    STELLA's Lagrangian grid often produces artifacts at the outer boundary:
    1. Velocity pile-up: the outermost zone is decelerated, causing an
       artificial density enhancement.
    2. Density non-monotonicity: the density oscillates in the last few
       zones instead of monotonically decreasing.

    This function operates on the final combined grid (shell+wind+ejecta)
    and cleans only the outermost ejecta zones.  It:
      - Extrapolates velocity from the smooth inner ejecta region.
      - Enforces monotonically decreasing density.
      - Tapers the last few zones to a smooth power-law continuation.

    Parameters:
        combined: dict with r, rho, v, T arrays (the final combined grid)
        n_edge:   number of outer zones to clean (default 20)
        ej_start_idx: optional pre-computed index where the ejecta begins.
                      If provided, the ejecta is cleaned from this index
                      outward; otherwise the inner edge is heuristically
                      identified from the largest velocity jump.

    Returns:
        combined dict with cleaned outer boundary
    """
    r = combined['r']
    rho = combined['rho']
    v = combined['v']  # cm/s
    n = len(r)

    if n < n_edge * 3:
        return combined

    # Identify the ejecta region.  If the caller has supplied the exact
    # ejecta start index (e.g. from the combiner), use it.  Otherwise fall
    # back to the largest velocity jump, not the first one, to avoid
    # misidentifying the shell core/envelope transition.
    if ej_start_idx is None:
        v_kms = v / 1e5
        dv = np.diff(v_kms)
        ej_start_idx = (int(np.argmax(dv)) if np.any(dv > 0) else n // 2) + 1
        # Safety: if the largest jump is in the inner half, take the outer
        # half maximum instead; ejecta is always the outer component.
        if ej_start_idx < n // 2:
            ej_start_idx = int(np.argmax(dv[n//2:])) + n//2 + 1
    n_ej = n - ej_start_idx

    if n_ej < n_edge * 2:
        return combined  # not enough ejecta zones to clean

    # Use the inner ejecta (excluding outer n_edge) for fitting
    clean_end = n - n_edge
    if clean_end <= ej_start_idx + 10:
        return combined  # not enough clean ejecta zones

    # --- 1. Fix velocity pile-up ---
    # Use linear extrapolation in linear space from the clean inner ejecta.
    # The velocity profile is nearly linear (v ~ r/t for homologous
    # expansion), so a linear extrapolation is stable and physical.
    # A quadratic fit in log-log space can diverge wildly, so we avoid it.
    r_clean = r[ej_start_idx:clean_end]
    v_clean = v[ej_start_idx:clean_end]
    # Linear fit: v = a * r + b
    coeffs = np.polyfit(r_clean, v_clean, 1)
    v_edge_fit = np.polyval(coeffs, r[clean_end:])

    # Smoothstep blend: actual -> fitted (smooth transition)
    blend = np.linspace(0, 1, n_edge)
    # Smoothstep (3u^2-2u^3): blends actual -> fitted with zero slope at both
    # ends of the edge region, avoiding a kink where the extrapolation starts.
    blend = blend * blend * (3.0 - 2.0 * blend)
    v[clean_end:] = (1 - blend) * v[clean_end:] + blend * v_edge_fit

    # Ensure monotonically increasing velocity in the ejecta
    for i in range(ej_start_idx + 1, n):
        if v[i] <= v[i - 1]:
            # Use the average velocity gradient from the clean region
            dv_avg = np.mean(np.diff(v_clean))
            v[i] = v[i - 1] + dv_avg

    # --- 2. Fix density non-monotonicity ---
    # Enforce monotonically decreasing density in the outer ejecta
    n_fixed = 0
    for i in range(clean_end, n - 1):
        if rho[i + 1] > rho[i]:
            if i > 0:
                ratio = rho[i] / max(rho[i - 1], 1e-300)
                rho[i + 1] = rho[i] * min(ratio, 0.99)
            else:
                rho[i + 1] = rho[i] * 0.9
            n_fixed += 1

    # Taper the last 5 zones using the median density gradient from
    # the clean ejecta region
    n_taper = min(5, n_edge)
    clean_ratios = rho[ej_start_idx + 1:clean_end] / rho[ej_start_idx:clean_end - 1]
    avg_ratio = np.median(clean_ratios[-20:])
    for i in range(n - n_taper, n):
        if i > 0:
            rho[i] = rho[i - 1] * avg_ratio

    if n_fixed > 0:
        print(f"  Combined outer boundary cleanup: fixed {n_fixed} non-monotonic "
              f"density zones, tapered last {n_taper}")
        print(f"    Velocity: extrapolated from inner ejecta (pile-up removed)")

    return combined


def _build_optimized_shell_grid(r_inner, r0, r_outer, nzon_shell,
                                core_fraction=0.25, min_dlogr=0.005):
    """Build a shell radial grid with split core/envelope resolution.

    The Dessart 2010jl shell has two distinct regions:
      - Core (r < r0): constant density, smooth T gradient
      - Envelope (r > r0): steep power-law density, steep T gradient

    A uniform log-spaced grid wastes zones in the flat core (where nothing
    interesting happens) and under-resolves the envelope (where the shock
    forms, radiation escapes, and recombination occurs).

    This function splits the grid: ``core_fraction`` of the zones go to the
    core, the rest to the envelope.  Both sections are log-spaced.

    If the envelope span is small, the envelope zone count is capped to
    respect ``min_dlogr`` and the surplus goes to the core (which has a
    much larger span and can absorb extra zones without violating the
    spacing constraint).

    Parameters:
        r_inner:       Inner radius of the shell (cm)
        r0:            Core/envelope transition radius (cm)
        r_outer:       Outer radius of the shell (cm)
        nzon_shell:    Total number of shell zones
        core_fraction: Fraction of zones for the core (default: 0.25)
        min_dlogr:     Minimum allowed dlog10(r) per zone (default: 0.005)

    Returns:
        np.ndarray of zone radii (length nzon_shell), monotonically increasing
    """
    if r0 <= r_inner or r0 >= r_outer:
        # No distinct core/envelope — fall back to uniform log spacing
        return np.logspace(np.log10(r_inner), np.log10(r_outer), nzon_shell)

    nzon_core = max(10, int(nzon_shell * core_fraction))
    nzon_envelope = nzon_shell - nzon_core

    # Cap envelope zones to respect min_dlogr (envelope has small span)
    env_span_dex = np.log10(r_outer) - np.log10(r0)
    max_env = max(int(env_span_dex / min_dlogr), 10)
    if nzon_envelope > max_env:
        surplus = nzon_envelope - max_env
        nzon_envelope = max_env
        nzon_core += surplus
        print(f"  Shell grid: envelope capped at {max_env} zones "
              f"(span={env_span_dex:.3f} dex, min_dlogr={min_dlogr:.4f}), "
              f"{surplus} extra zones -> core")

    # Also cap core zones (core has large span, rarely hits this)
    core_span_dex = np.log10(r0) - np.log10(r_inner)
    max_core = max(int(core_span_dex / min_dlogr), 10)
    if nzon_core > max_core:
        nzon_core = max_core

    r_core = np.logspace(np.log10(r_inner), np.log10(r0), nzon_core + 1)[:nzon_core]
    r_envelope = np.logspace(np.log10(r0), np.log10(r_outer), nzon_envelope + 1)[1:]
    return np.concatenate([r_core, r_envelope])


def validate_stella_grid(r, min_dlogr=0.005, max_dlogr=0.15, strict=False):
    """Validate that the radial grid meets STELLA spacing constraints.

    STELLA's implicit radiation-hydrodynamics solver is highly sensitive to
    grid spacing.  If adjacent zones are too close (dlog10(r) too small),
    the radiation diffusion timestep becomes infinitesimally small, causing
    the simulation to stall or crash.  If zones are too far apart, shocks
    and contact discontinuities are under-resolved.

    Parameters:
        r:          Radial grid (monotonically increasing, cm)
        min_dlogr:  Minimum allowed dlog10(r) between adjacent zones.
                    Default 0.005 corresponds to dr/r ~ 0.011 (1.1%).
        max_dlogr:  Maximum allowed dlog10(r) between adjacent zones.
                    Default 0.15 corresponds to dr/r ~ 0.41 (41%).
        strict:     If True, raise ValueError on min_dlogr violation
                    instead of printing a warning.

    Returns:
        True if the grid passes all checks, False if warnings were printed.
    """
    log_r = np.log10(r)
    dlogr = np.diff(log_r)
    dr_over_r = np.diff(r) / r[:-1]

    min_dl = np.min(dlogr)
    max_dl = np.max(dlogr)
    mean_dl = np.mean(dlogr)
    min_dr = np.min(dr_over_r)
    max_dr = np.max(dr_over_r)
    mean_dr = np.mean(dr_over_r)

    print(f"\n  === STELLA Grid Validation ===")
    print(f"    Zones: {len(r)}")
    print(f"    dlog10(r): min={min_dl:.5f}, max={max_dl:.5f}, mean={mean_dl:.5f}")
    print(f"    dr/r:      min={min_dr:.5f}, max={max_dr:.5f}, mean={mean_dr:.5f}")

    ok = True

    # --- Minimum spacing check (grid trap) ---
    if min_dl < min_dlogr:
        idx_min = np.argmin(dlogr)
        r_at_min = r[idx_min]
        msg = (
            f"CRITICAL: Grid spacing too small at zone {idx_min} "
            f"(r={r_at_min:.3e} cm): dlog10(r)={min_dl:.5f} < {min_dlogr:.5f} "
            f"(dr/r={min_dr:.5f}).\n"
            f"    This will kill STELLA's radiation diffusion timestep. "
            f"The solver will stall or crash.\n"
            f"    Fix: reduce the number of zones in this region, or widen "
            f"the physical span."
        )
        if strict:
            raise ValueError(msg)
        else:
            print(f"  {msg}")
            ok = False

    # --- Maximum spacing check (under-resolution) ---
    if max_dl > max_dlogr:
        idx_max = np.argmax(dlogr)
        r_at_max = r[idx_max]
        print(f"  WARNING: Grid spacing too large at zone {idx_max} "
              f"(r={r_at_max:.3e} cm): dlog10(r)={max_dl:.5f} > {max_dlogr:.5f} "
              f"(dr/r={max_dr:.5f}).")
        print(f"    Shocks and contact discontinuities may be under-resolved "
              f"in this region.")
        ok = False

    if ok:
        print(f"    Grid OK: all spacing constraints satisfied.")
    return ok


def _cap_buffer_zones(nzon_buffer, span_dex, min_dlogr=0.005, max_dlogr=None,
                      label="buffer"):
    """Cap buffer zone count to respect dlog10(r) spacing constraints.

    Two constraints are enforced:

    1. **Minimum spacing (cap DOWN):** If the physical span is small but the
       requested zone count is high, the resulting dlogr = span_dex /
       nzon_buffer can be tiny, creating a "grid trap" that stalls STELLA's
       radiation diffusion timestep.  This reduces nzon_buffer so that
       dlogr >= min_dlogr.

    2. **Maximum spacing (cap UP):** If max_dlogr is provided and the
       resulting dlogr exceeds max_dlogr, the zone count is *increased* so
       that dlogr <= max_dlogr.  This prevents grid spacing discontinuities
       at buffer boundaries — when a buffer's dlogr is much coarser than
       the adjacent wind/ejecta grids, the GEAR solver's Jacobian becomes
       structurally singular (IFAIL=3) and the timestep collapses to
       denormalized floats.  The classic failure mode is a 10x dlogr jump
       at a wind-ejecta buffer boundary when the shell shock reaches it.

    Parameters:
        nzon_buffer: Requested number of buffer zones
        span_dex:    Physical span in dex of log10(r)
        min_dlogr:   Minimum allowed dlog10(r) per zone
        max_dlogr:   Maximum allowed dlog10(r) per zone (optional).
                     When set, ensures dlogr <= max_dlogr by increasing
                     the zone count.  Typically set to 2*min_dlogr so
                     the buffer dlogr is at most 2x the adjacent grid dlogr.
        label:       Buffer name for diagnostic messages

    Returns:
        Adjusted nzon_buffer (may be > or < input)
    """
    if span_dex <= 0 or nzon_buffer <= 0:
        return nzon_buffer

    # Cap DOWN: ensure dlogr >= min_dlogr (prevent grid traps)
    required_dlogr = span_dex / nzon_buffer
    if required_dlogr < min_dlogr:
        max_zones = int(span_dex / min_dlogr)
        max_zones = max(max_zones, 2)  # need at least 2 zones
        print(f"  {label} zone cap (min): {nzon_buffer} -> {max_zones} zones "
              f"(span={span_dex:.3f} dex, min_dlogr={min_dlogr:.4f}, "
              f"was dlogr={required_dlogr:.5f})")
        nzon_buffer = max_zones
        required_dlogr = span_dex / nzon_buffer

    # Cap UP: ensure dlogr <= max_dlogr (prevent grid discontinuities)
    if max_dlogr is not None and required_dlogr > max_dlogr:
        min_zones = int(np.ceil(span_dex / max_dlogr))
        min_zones = max(min_zones, 2)
        print(f"  {label} zone cap (max): {nzon_buffer} -> {min_zones} zones "
              f"(span={span_dex:.3f} dex, max_dlogr={max_dlogr:.4f}, "
              f"was dlogr={required_dlogr:.5f})")
        nzon_buffer = min_zones

    return nzon_buffer


def combine_shell_ejecta(shell, ejecta, nzon_total=500,
                          se_transition_dex=0.5, immediate_contact=False,
                          ejecta_T_floor=1000.0, min_dlogr=0.005,
                          max_dlogr=None):
    """Combine shell (inner) + ejecta (outer) with a buffer zone (no wind).

    The buffer spans the full gap from the shell outer edge to the ejecta
    inner edge, with log-linear density ramping and smoothly blended
    composition (linear in log-r from shell H/He to ejecta metal-rich).

    When immediate_contact=True, the buffer is eliminated entirely — the
    shell outer edge directly meets the ejecta inner edge with a sharp
    density discontinuity (handled by STELLA's AQ artificial viscosity).

    Parameters:
        se_transition_dex: Unused (kept for CLI compatibility). The buffer
            always extends to the ejecta inner edge.
        immediate_contact: If True, no buffer zones — shell and ejecta are
            in direct contact with a sharp density jump.
    """
    from scipy.interpolate import PchipInterpolator

    # Maximum dlogr for buffers and adjacent regions: ensure no region's
    # dlogr exceeds 2x the minimum grid dlogr, preventing grid spacing
    # discontinuities at buffer boundaries that collapse the GEAR timestep.
    # The global max_dlogr (default 0.15) is a hard upper bound only; the
    # local smoothness constraint (2*min_dlogr) is tighter and is what
    # actually eliminates the factor-of-10 dlogr jumps seen at 80 days.
    buffer_max_dlogr = (min(max_dlogr, 2 * min_dlogr)
                        if max_dlogr is not None else 2 * min_dlogr)

    # Sort ejecta inner->outer
    ej_order = np.argsort(ejecta['r'])
    ej_r_sorted = ejecta['r'][ej_order]
    ej_rho_sorted = ejecta['rho'][ej_order]
    ej_v_sorted = ejecta['v'][ej_order]
    ej_T_sorted = ejecta['T'][ej_order]

    # Region boundaries
    r_min = shell['r_inner']
    r_shell_outer = shell['r_outer']
    # When immediate_contact=True, the shell outer edge was set to
    # ej_r_sorted[0] * 0.99 by the caller (r_outer_override).  The ejecta
    # inner edge for the grid must therefore be ej_r_sorted[0] (without the
    # 0.99 factor) so the 5-zone buffer has a non-zero radial extent.
    # Without this, r_shell_outer == r_ej_inner and the buffer zones all
    # have identical radius, producing zero-mass zones that STELLA rejects
    # with "Invalid .hyd mass or radial ordering".
    if immediate_contact:
        r_ej_inner = ej_r_sorted[0]
    else:
        r_ej_inner = ej_r_sorted[0] * 0.99
    r_max = ej_r_sorted[-1]

    log_r_min = np.log10(r_min)
    log_r_shell = np.log10(r_shell_outer)
    log_r_ej = np.log10(r_ej_inner)
    log_r_max = np.log10(r_max)

    if immediate_contact:
        # Immediate contact with a small smooth transition (5 zones).
        # The shell outer edge is at ~99% of the ejecta inner edge (set by
        # r_outer_override).  The 5-zone buffer spans the remaining 1% gap,
        # providing a smooth density/velocity/composition blend to avoid a
        # hard discontinuity that could stall STELLA.
        nzon_se_buffer = 5
        print(f"  Shell-ejecta: IMMEDIATE CONTACT (5-zone smooth transition)")
        print(f"    r_shell_outer={r_shell_outer:.3e} cm, r_ej_inner={r_ej_inner:.3e} cm")
        print(f"    buffer width={r_ej_inner - r_shell_outer:.3e} cm "
              f"({(r_ej_inner - r_shell_outer)/r_ej_inner*100:.2f}% of r_ej_inner)")
        rho_se_start = shell['rho'][-1]
        rho_se_end = ej_rho_sorted[0]
        se_contrast = abs(np.log10(max(rho_se_end, 1e-300) / max(rho_se_start, 1e-300)))
        print(f"    rho: {rho_se_start:.3e} -> {rho_se_end:.3e} g/cc "
              f"(contrast={se_contrast:.1f} dex, smooth over {nzon_se_buffer} zones)")
    else:
        # --- Buffer spans the full shell->ejecta gap ---
        rho_se_start = shell['rho'][-1]  # shell density at outer edge
        rho_se_end = ej_rho_sorted[0]    # ejecta density at inner edge
        se_contrast = abs(np.log10(max(rho_se_end, 1e-300) / max(rho_se_start, 1e-300)))
        se_span_dex = log_r_ej - log_r_shell  # full gap width in dex

        # Number of buffer zones: enough for density gradient < 16 and
        # composition resolution (~0.05 dex per zone minimum)
        nzon_se_buffer = max(30, int(se_span_dex / 0.05))
        nzon_se_buffer = min(nzon_se_buffer, 150)  # cap at 150
        # Enforce min_dlogr to prevent grid traps, and max_dlogr to prevent
        # grid discontinuities at buffer boundaries.
        nzon_se_buffer = _cap_buffer_zones(nzon_se_buffer, se_span_dex,
                                           min_dlogr=min_dlogr,
                                           max_dlogr=buffer_max_dlogr,
                                           label="SE buffer")

        se_grad = se_contrast / max(se_span_dex, 0.01)
        print(f"  Shell-ejecta buffer: {nzon_se_buffer} zones, "
              f"r=[{r_shell_outer:.3e}, {r_ej_inner:.3e}] cm "
              f"(span={se_span_dex:.2f} dex)")
        print(f"    rho: {rho_se_start:.3e} -> {rho_se_end:.3e} g/cc "
              f"(contrast={se_contrast:.1f} dex, est. gradient={se_grad:.1f})")

    # --- Zone allocation (must sum to exactly nzon_total) ---
    if immediate_contact:
        # 5-zone transition, rest split between shell and ejecta
        remaining = nzon_total - nzon_se_buffer
        nzon_shell = max(80, int(remaining * 0.30))
        nzon_ejecta = remaining - nzon_shell
    else:
        # Buffer gets its zones first, then shell/ejecta split the remainder
        remaining = nzon_total - nzon_se_buffer
        nzon_shell = max(80, int(remaining * 0.30))
        nzon_ejecta = remaining - nzon_shell
    # Cap ejecta zones to respect min_dlogr (prevent grid traps) and
    # max_dlogr (prevent grid discontinuities at buffer boundaries).
    ej_span_dex = log_r_max - log_r_ej
    nzon_ejecta = _cap_buffer_zones(nzon_ejecta, ej_span_dex,
                                    min_dlogr=min_dlogr,
                                    max_dlogr=buffer_max_dlogr,
                                    label="Ejecta")
    # Redistribute freed zones to the shell (envelope-limited)
    n_freed = nzon_total - nzon_shell - nzon_se_buffer - nzon_ejecta
    if n_freed < 0:
        # Buffer zones grew (max_dlogr enforcement) and the total exceeds
        # nzon_total.  Reduce the ejecta first, then shell.
        n_excess = -n_freed
        n_from_ej = min(n_excess, max(0, nzon_ejecta - 100))
        nzon_ejecta -= n_from_ej
        n_excess -= n_from_ej
        n_from_shell = min(n_excess, max(0, nzon_shell - 80))
        nzon_shell -= n_from_shell
        n_excess -= n_from_shell
        if n_from_ej > 0 or n_from_shell > 0:
            print(f"  Buffer max_dlogr enforcement grew buffers by {-n_freed} zones; "
                  f"reduced ejecta by {n_from_ej}, shell by {n_from_shell}")
        if n_excess > 0:
            print(f"  WARNING: {n_excess} zones could not be removed without "
                  f"violating minimums — total will be {nzon_total + n_excess}")
        n_freed = 0
    if n_freed > 0:
        shell_span_dex = np.log10(r_shell_outer) - np.log10(r_min)
        max_shell = max(int(shell_span_dex / min_dlogr), nzon_shell)
        n_to_shell = min(n_freed, max_shell - nzon_shell)
        if n_to_shell > 0:
            nzon_shell += n_to_shell
            n_freed -= n_to_shell
            print(f"  Redistributed {n_to_shell} freed zones to shell")
        if n_freed > 0:
            print(f"  WARNING: {n_freed} zones could not be placed without "
                  f"violating min_dlogr — total will be {nzon_total - n_freed}")
    n_total = nzon_shell + nzon_se_buffer + nzon_ejecta
    print(f"  Zone allocation: {nzon_shell} shell + {nzon_se_buffer} SE buffer + "
          f"{nzon_ejecta} ejecta = {n_total} total")

    # --- Build grids ---
    # Shell: split core/envelope grid for better envelope resolution.
    # The Dessart power-law shell has a flat core (r < r0) and a steep
    # envelope (r > r0) where the shock forms, radiation escapes, and
    # recombination occurs.  Concentrating zones in the envelope gives
    # 5x better resolution of the recombination front and shock.
    r0_shell = shell.get('r0', None)
    if r0_shell is not None and r0_shell > shell['r_inner']:
        shell_r = _build_optimized_shell_grid(
            r_min, r0_shell, r_shell_outer, nzon_shell, core_fraction=0.25,
            min_dlogr=min_dlogr)
        nzon_core = np.sum(shell_r < r0_shell)
        print(f"  Shell grid: {nzon_core} core + {nzon_shell - nzon_core} envelope "
              f"(split at r0={r0_shell:.2e} cm)")
    else:
        shell_r = np.logspace(log_r_min, log_r_shell, nzon_shell)
    log_r_src = np.log10(shell['r'])
    log_rho_src = np.log10(np.maximum(shell['rho'], 1e-300))
    shell_rho = 10.0 ** np.interp(np.log10(shell_r), log_r_src, log_rho_src)
    shell_v = np.interp(np.log10(shell_r), log_r_src, shell['v'])
    shell_T = np.interp(np.log10(shell_r), log_r_src, shell['T'])

    # Conserve shell mass: the log-space interpolation of a steep power-law
    # (rho ∝ r^-8) onto fewer zones under-weights the dense core, and the
    # original 100-zone discrete grid also loses mass vs the analytical
    # formula. Renormalize to the analytical mass (exact by construction).
    if 'rho_0' in shell and 'r0' in shell:
        n_rho_val = shell.get('n_rho', 8.0)
        frac_v0 = shell['r0'] / shell['r_outer']
        factor = _shell_mass_factor(n_rho_val, frac_v0)
        M_shell_target = (4/3) * np.pi * shell['r0']**3 * shell['rho_0'] * factor
    else:
        # Fallback: integrate the original discrete profile
        r_lo_src = np.concatenate([[shell['r_inner']], shell['r'][:-1]])
        M_shell_target = ((4/3) * np.pi *
                          (shell['r']**3 - r_lo_src**3) * shell['rho']).sum()
    r_lo_shell = np.concatenate([[shell_r[0] * 0.5], shell_r[:-1]])
    M_shell_interp = ((4/3) * np.pi *
                      (shell_r**3 - r_lo_shell**3) * shell_rho).sum()
    if M_shell_interp > 0:
        shell_rho *= M_shell_target / M_shell_interp
        print(f"  Shell mass conservation: target={M_shell_target/MSUN:.5f} Msun, "
              f"interp={M_shell_interp/MSUN:.5f} Msun, "
              f"renorm factor={M_shell_target/M_shell_interp:.4f}")

    if nzon_se_buffer > 0:
        # Buffer: log-spaced from shell outer edge to ejecta inner edge
        se_buf_r = np.logspace(log_r_shell, log_r_ej, nzon_se_buffer + 1)[1:]
        se_log_rho = np.linspace(np.log10(max(rho_se_start, 1e-300)),
                                 np.log10(max(rho_se_end, 1e-300)),
                                 nzon_se_buffer)
        se_buf_rho = 10.0 ** se_log_rho
        # Velocity: sharp linear-ramp transition from shell v to ejecta v.
        # The transition happens over ~4 zones at the buffer midpoint,
        # so most of the buffer has shell velocity (inner) or ejecta
        # velocity (outer) with a sharp interface.
        v_shell_edge = shell['v'][-1]
        v_ej_inner = ej_v_sorted[0]
        se_buf_v = _interface_blend(nzon_se_buffer, v_shell_edge, v_ej_inner,
                                    transition_zones=4)
        # Temperature: blend from shell T to ejecta T across the transition.
        # For immediate contact (thin 5-zone transition), this provides a
        # smooth temperature gradient at the contact discontinuity.
        # For the wide buffer case, the buffer is cold (ejecta T floor).
        if immediate_contact:
            se_buf_T = _interface_blend(nzon_se_buffer, shell['T'][-1],
                                        max(ej_T_sorted[0], ejecta_T_floor),
                                        transition_zones=4)
        else:
            se_buf_T = np.full(nzon_se_buffer, max(ej_T_sorted[0], ejecta_T_floor))
    else:
        # No buffer (immediate contact): empty arrays
        se_buf_r = np.array([])
        se_buf_rho = np.array([])
        se_buf_v = np.array([])
        se_buf_T = np.array([])

    # Ejecta: log-spaced, interpolated from source
    ej_r_new = np.logspace(log_r_ej, log_r_max, nzon_ejecta + 1)[1:]
    log_r_ej_src = np.log10(ej_r_sorted)

    # PCHIP (piecewise cubic Hermite) interpolation: shape-preserving and
    # monotone where the data are monotone, so it cannot produce the Gibbs
    # wiggles that linear-in-log-space interpolation creates at steep
    # density/temperature gradients.  Density and temperature are
    # interpolated in log-space; velocity and composition in linear space.
    # Out-of-range radii are clamped to the source bounds (constant
    # extrapolation, matching the old fill_value behaviour).
    def log_interp(y, r_target):
        log_y = np.log10(np.maximum(y, 1e-300))
        ip = PchipInterpolator(log_r_ej_src, log_y, extrapolate=False)
        lr = np.clip(np.log10(r_target), log_r_ej_src[0], log_r_ej_src[-1])
        return 10.0 ** ip(lr)

    def lin_interp(y, r_target):
        ip = PchipInterpolator(log_r_ej_src, y, extrapolate=False)
        lr = np.clip(np.log10(r_target), log_r_ej_src[0], log_r_ej_src[-1])
        return ip(lr)

    ej_rho = log_interp(ej_rho_sorted, ej_r_new)
    ej_v = lin_interp(ej_v_sorted, ej_r_new)
    ej_T = np.maximum(log_interp(ej_T_sorted, ej_r_new), ejecta_T_floor)

    # Smooth ejecta T profile: the STELLA .swd T at late epochs has
    # numerical noise from the radiation-hydrodynamics evolution (zone-to-zone
    # scatter of ~20%). At 400d the ejecta is thermally equilibrated (radiation
    # transport has smoothed gradients), so a light moving-average smooth is
    # physically justified and gives a more robust initial condition.
    if len(ej_T) > 5:
        kernel = np.ones(5) / 5.0
        ej_T_smooth = np.convolve(ej_T, kernel, mode='same')
        # Preserve edges (convolve shrinks valid range)
        ej_T_smooth[:2] = ej_T[:2]
        ej_T_smooth[-2:] = ej_T[-2:]
        ej_T = np.maximum(ej_T_smooth, ejecta_T_floor)

    # Smooth ejecta interpolation artifacts
    for _ in range(3):
        log_ej_rho = np.log10(np.maximum(ej_rho, 1e-300))
        log_ej_r = np.log10(ej_r_new)
        grads = np.abs(np.diff(log_ej_rho) / np.diff(log_ej_r))
        smoothed = ej_rho.copy()
        for i in range(1, len(ej_rho) - 1):
            g_left = grads[i - 1]
            g_right = grads[i]
            if g_left > 16 or g_right > 16:
                smoothed[i] = np.sqrt(max(ej_rho[i - 1] * ej_rho[i + 1], 1e-300))
        ej_rho = smoothed

    # --- Combine all segments ---
    r_combined = np.concatenate([shell_r, se_buf_r, ej_r_new])
    rho_combined = np.concatenate([shell_rho, se_buf_rho, ej_rho])
    v_combined = np.concatenate([shell_v, se_buf_v, ej_v])
    T_combined = np.concatenate([shell_T, se_buf_T, ej_T])

    # --- Composition (smoothly blended across buffer) ---
    shell_xh = shell.get('fH', [0.7])[0]
    shell_yhe = shell.get('fHe', [0.3])[0]

    # Ejecta composition at its inner edge (for buffer blending)
    all_elem_keys = ['fH', 'fHe', 'fO', 'fSi', 'fCa', 'fFe']
    ej_comp_inner = {}
    for key in all_elem_keys:
        if key in ejecta:
            ej_comp_inner[key] = float(ejecta[key][ej_order][0])
        else:
            ej_comp_inner[key] = 0.0

    # Ejecta composition interpolated onto ej_r_new
    ej_comp = {}
    for key in all_elem_keys:
        if key in ejecta:
            ej_comp[key] = lin_interp(ejecta[key][ej_order], ej_r_new)
        else:
            ej_comp[key] = np.zeros(nzon_ejecta)

    # Shell composition values (supports solar composition with metals)
    shell_comp_vals = {}
    for key in all_elem_keys:
        shell_comp_vals[key] = float(shell.get(key, [0.0])[0])

    # Blend weights for buffer: sharp linear ramp centered at buffer midpoint.
    # Transitions over ~4 zones so the composition interface is sharp
    # (shell H/He -> ejecta metals over 3-5 zones, not the full buffer).
    se_w = _interface_ramp_weights(nzon_se_buffer, transition_zones=4)

    comp_combined = {}
    for key in all_elem_keys:
        sv = shell_comp_vals[key]
        ev = ej_comp_inner[key]
        # Buffer: linear blend from shell to ejecta inner-edge composition
        se_buf_comp = sv * (1.0 - se_w) + ev * se_w

        comp_combined[key] = np.concatenate([
            np.full(nzon_shell, sv),
            se_buf_comp,
            ej_comp[key],
        ])

    # Renormalize composition
    comp_sum = sum(comp_combined.values())
    comp_sum = np.maximum(comp_sum, 1e-10)
    for key in comp_combined:
        comp_combined[key] /= comp_sum

    # --- Compute masses ---
    r_lo = np.concatenate([[r_combined[0] * 0.5], r_combined[:-1]])
    dm = (4/3) * np.pi * (r_combined**3 - r_lo**3) * rho_combined
    dm = np.maximum(dm, 0.0)
    M_total = np.sum(dm)

    idx_shell = nzon_shell
    idx_se_buf = idx_shell + nzon_se_buffer

    print(f"  Mass: shell={np.sum(dm[:idx_shell])/MSUN:.3f} Msun, "
          f"SE buffer={np.sum(dm[idx_shell:idx_se_buf])/MSUN:.4f} Msun, "
          f"ejecta={np.sum(dm[idx_se_buf:])/MSUN:.3f} Msun, "
          f"total={M_total/MSUN:.3f} Msun")

    # --- Check density gradients ---
    log_r = np.log10(r_combined)
    log_rho_arr = np.log10(np.maximum(rho_combined, 1e-300))
    if len(log_r) > 1:
        gradients = np.abs(np.diff(log_rho_arr) / np.diff(log_r))
        max_grad = np.max(gradients)
        grad_interior = gradients[:-5]
        max_grad_interior = np.max(grad_interior) if len(grad_interior) > 0 else 0
        print(f"  Max |dlog rho/dlog r|: {max_grad:.1f} (interior: {max_grad_interior:.1f})")
        if max_grad_interior > 16:
            print(f"  WARNING: Interior density gradient {max_grad_interior:.1f} > 16 "
                  f"— may cause solver stall")
            worst_idx = np.argmax(grad_interior)
            print(f"    Worst at zone {worst_idx}: r={r_combined[worst_idx]:.3e} cm, "
                  f"rho={rho_combined[worst_idx]:.3e} -> {rho_combined[worst_idx+1]:.3e}")
        else:
            print(f"  Interior gradients OK (< 16, buffer zones effective)")

        se_buf_grads = gradients[idx_shell-1:idx_se_buf]
        if len(se_buf_grads) > 0:
            print(f"  SE buffer gradients: max={np.max(se_buf_grads):.1f}, "
                  f"mean={np.mean(se_buf_grads):.1f}")

    return {
        'r': r_combined,
        'rho': rho_combined,
        'v': v_combined,
        'T': T_combined,
        'idx_ej': idx_se_buf,
        **comp_combined,
    }


def combine_shell_wind_ejecta(shell, wind, ejecta, nzon_total=500,
                               sw_transition_dex=0.5, sw_transition_k=3.0,
                               we_transition_dex=0.5, we_transition_k=3.0,
                               ejecta_T_floor=1000.0, min_dlogr=0.005,
                               max_dlogr=None):
    """Combine shell (inner) + wind (middle) + ejecta (outer) with buffer zones.

    Inserts log-linear buffer zones at the shell-wind and wind-ejecta boundaries
    to smooth the density cliff. The shell density is preserved unchanged; the
    buffer zone is placed OUTSIDE the shell, ramping from rho_shell to rho_wind
    in log-space over enough zones to keep dlog(rho)/dlog(r) < 16.

    This avoids the mass-loss problem of sigmoid blending (which destroys shell
    mass when the density contrast is ~10 orders of magnitude) by keeping each
    component's profile unchanged outside the buffer region.

    Parameters:
        sw_transition_dex: Target width of shell-wind buffer in dex of log10(r)
        sw_transition_k:   Unused (kept for CLI compatibility)
        we_transition_dex: Target width of wind-ejecta buffer in dex of log10(r)
        we_transition_k:   Unused (kept for CLI compatibility)
        max_dlogr:         Maximum allowed dlog10(r) per zone for buffers and
                          adjacent regions.  When None, defaults to
                          2*min_dlogr, ensuring the buffer dlogr is at most
                          2x the adjacent grid dlogr.  This prevents the GEAR
                          solver from failing with IFAIL=3 (singular matrix)
                          when the shell shock crosses a coarse buffer.
    """
    from scipy.interpolate import PchipInterpolator

    # Maximum dlogr for buffers and adjacent regions: ensure no region's
    # dlogr exceeds 2x the minimum grid dlogr, preventing grid spacing
    # discontinuities at buffer boundaries that collapse the GEAR timestep.
    # The global max_dlogr (default 0.15) is a hard upper bound only; the
    # local smoothness constraint (2*min_dlogr) is tighter and is what
    # actually eliminates the factor-of-10 dlogr jumps seen at 80 days.
    buffer_max_dlogr = (min(max_dlogr, 2 * min_dlogr)
                        if max_dlogr is not None else 2 * min_dlogr)

    # Sort ejecta inner->outer
    ej_order = np.argsort(ejecta['r'])
    ej_r_sorted = ejecta['r'][ej_order]
    ej_rho_sorted = ejecta['rho'][ej_order]
    ej_v_sorted = ejecta['v'][ej_order]
    ej_T_sorted = ejecta['T'][ej_order]

    # Region boundaries
    r_min = shell['r_inner']
    r_shell_outer = shell['r_outer']
    r_ej_inner = ej_r_sorted[0] * 0.99
    r_max = ej_r_sorted[-1]

    log_r_min = np.log10(r_min)
    log_r_shell = np.log10(r_shell_outer)
    log_r_ej = np.log10(r_ej_inner)
    log_r_max = np.log10(r_max)

    # Wind properties
    v_wind = wind['v'][0]
    T_wind = wind['T'][0]
    # Wind density on the final grid is interpolated from the supplied
    # profile rather than rebuilt from a single Mdot: a finite-age wind
    # (--wind-duration-days) concatenates two Mdot segments and its
    # leading edge must survive.  For a steady 1/r^2 wind the
    # interpolation is exact in log-log space, so behaviour is unchanged.
    _w_src_lr = np.log10(wind['r'])
    _w_src_lrho = np.log10(np.maximum(wind['rho'], 1e-300))

    def wind_rho_at(rr):
        return 10.0 ** np.interp(np.log10(rr), _w_src_lr, _w_src_lrho)

    # --- Compute buffer zone sizes ---
    # Shell-wind buffer: ramps from rho_shell_edge to rho_wind(r_buffer_end)
    # The buffer extends from r_shell_outer outward by sw_transition_dex in log(r)
    # Uses the shell density at its OUTER EDGE (not core density) for the transition
    #
    # Adaptive buffer widths: if the shell-ejecta gap is small (young source age),
    # the fixed sw_transition_dex and we_transition_dex can cause the buffers to
    # overlap, making the wind grid go backwards. Cap each buffer at 1/3 of the
    # total gap to ensure the wind always has positive radial extent.
    total_gap_dex = log_r_ej - log_r_shell
    max_buffer_dex = total_gap_dex / 3.0  # each buffer gets at most 1/3 of gap
    sw_transition_dex = min(sw_transition_dex, max_buffer_dex)
    we_transition_dex = min(we_transition_dex, max_buffer_dex)

    r_sw_buffer_end = r_shell_outer * (10 ** sw_transition_dex)
    rho_sw_start = shell['rho'][-1]  # shell density at outer edge (power-law tail)
    rho_sw_end = wind_rho_at(r_sw_buffer_end)  # wind density at buffer end
    sw_contrast = np.log10(max(rho_sw_start, 1e-300) / max(rho_sw_end, 1e-300))

    # Number of buffer zones: enough for gradient < 16
    # gradient = contrast / (n_buffer * dlog_r_per_zone)
    # We want gradient < 16, so n_buffer > contrast / (16 * dlog_r_per_zone)
    # Use dlog_r_per_zone ≈ sw_transition_dex / n_buffer (self-consistent)
    # gradient = contrast / sw_transition_dex (independent of n_buffer!)
    # So we need sw_transition_dex > contrast / 16
    sw_grad = sw_contrast / sw_transition_dex
    nzon_sw_buffer = max(20, int(sw_contrast / 16 * sw_transition_dex / 0.02))
    nzon_sw_buffer = min(nzon_sw_buffer, 60)  # cap at 60
    # Enforce min_dlogr to prevent grid traps, and max_dlogr to prevent
    # grid discontinuities at buffer boundaries.
    nzon_sw_buffer = _cap_buffer_zones(nzon_sw_buffer, sw_transition_dex,
                                       min_dlogr=min_dlogr,
                                       max_dlogr=buffer_max_dlogr,
                                       label="SW buffer")

    print(f"  Shell-wind buffer: {nzon_sw_buffer} zones, "
          f"r=[{r_shell_outer:.3e}, {r_sw_buffer_end:.3e}] cm")
    print(f"    rho: {rho_sw_start:.3e} -> {rho_sw_end:.3e} g/cc "
          f"(contrast={sw_contrast:.1f} dex, est. gradient={sw_grad:.1f})")

    # Wind-ejecta buffer: ramps from rho_wind to rho_ejecta
    r_we_buffer_start = r_ej_inner / (10 ** we_transition_dex)
    rho_we_start = wind_rho_at(r_we_buffer_start)
    rho_we_end = ej_rho_sorted[0]  # ejecta density at inner edge
    we_contrast = abs(np.log10(max(rho_we_end, 1e-300) / max(rho_we_start, 1e-300)))
    we_grad = we_contrast / we_transition_dex
    nzon_we_buffer = max(10, int(we_contrast / 16 * we_transition_dex / 0.02))
    nzon_we_buffer = min(nzon_we_buffer, 40)
    # Enforce min_dlogr to prevent grid traps, and max_dlogr to prevent
    # grid discontinuities at buffer boundaries.
    nzon_we_buffer = _cap_buffer_zones(nzon_we_buffer, we_transition_dex,
                                       min_dlogr=min_dlogr,
                                       max_dlogr=buffer_max_dlogr,
                                       label="WE buffer")

    print(f"  Wind-ejecta buffer: {nzon_we_buffer} zones, "
          f"r=[{r_we_buffer_start:.3e}, {r_ej_inner:.3e}] cm")
    print(f"    rho: {rho_we_start:.3e} -> {rho_we_end:.3e} g/cc "
          f"(contrast={we_contrast:.1f} dex, est. gradient={we_grad:.1f})")

    # --- Zone allocation ---
    # Optimized: shell gets 20% (up from 15%) — the shell is where the shock,
    # recombination, and radiation escape happen.  Wind gets fewer zones
    # (it's a smooth 1/r^2 profile).  Ejecta gets the remainder.
    nzon_shell = max(80, int(nzon_total * 0.20))
    nzon_ejecta = max(100, int(nzon_total * 0.50))
    nzon_wind = nzon_total - nzon_shell - nzon_ejecta - nzon_sw_buffer - nzon_we_buffer
    nzon_wind = max(40, nzon_wind)  # ensure minimum wind zones

    # Cap wind and ejecta zones to respect min_dlogr (prevent grid traps).
    # The wind and ejecta spans are known from the region boundaries.
    # Wind span excludes the WE buffer (which has its own zone count).
    wind_span_dex = np.log10(r_we_buffer_start) - np.log10(r_sw_buffer_end)
    ej_span_dex = log_r_max - log_r_ej
    nzon_wind_capped = _cap_buffer_zones(nzon_wind, wind_span_dex,
                                         min_dlogr=min_dlogr,
                                         max_dlogr=buffer_max_dlogr,
                                         label="Wind")
    nzon_ejecta_capped = _cap_buffer_zones(nzon_ejecta, ej_span_dex,
                                           min_dlogr=min_dlogr,
                                           max_dlogr=buffer_max_dlogr,
                                           label="Ejecta")

    # Redistribute freed zones to the shell (which has plenty of physical
    # span and can always use more envelope resolution for shock/recombination).
    # The shell is capped by its own min_dlogr constraint (envelope-limited).
    # Any zones that cannot be absorbed are dropped (the total may be < nzon_total).
    nzon_wind = nzon_wind_capped
    nzon_ejecta = nzon_ejecta_capped
    n_freed = (nzon_total - nzon_shell - nzon_sw_buffer - nzon_wind
               - nzon_we_buffer - nzon_ejecta)
    if n_freed < 0:
        # Buffer zones grew (max_dlogr enforcement) and the total exceeds
        # nzon_total.  Reduce the ejecta (most zones) first, then wind,
        # then shell, to bring the total back to nzon_total.
        n_excess = -n_freed
        # Reduce ejecta first (it has the most zones and the coarsest
        # physical resolution requirement at large radii).
        n_from_ej = min(n_excess, max(0, nzon_ejecta - 100))
        nzon_ejecta -= n_from_ej
        n_excess -= n_from_ej
        # Then reduce wind
        n_from_wind = min(n_excess, max(0, nzon_wind - 40))
        nzon_wind -= n_from_wind
        n_excess -= n_from_wind
        # Then reduce shell
        n_from_shell = min(n_excess, max(0, nzon_shell - 80))
        nzon_shell -= n_from_shell
        n_excess -= n_from_shell
        if n_from_ej > 0 or n_from_wind > 0 or n_from_shell > 0:
            print(f"  Buffer max_dlogr enforcement grew buffers by {-n_freed} zones; "
                  f"reduced ejecta by {n_from_ej}, wind by {n_from_wind}, "
                  f"shell by {n_from_shell}")
        if n_excess > 0:
            print(f"  WARNING: {n_excess} zones could not be removed without "
                  f"violating minimums — total will be {nzon_total + n_excess}")
        n_freed = 0
    if n_freed > 0:
        # Cap shell zones: the shell envelope has a small span and is the
        # limiting factor.  The _build_optimized_shell_grid function will
        # redistribute envelope-overflow to the core, but the total shell
        # zones is still limited by the combined core+envelope spans.
        shell_span_dex = np.log10(r_shell_outer) - np.log10(r_min)
        max_shell = max(int(shell_span_dex / min_dlogr), nzon_shell)
        n_to_shell = min(n_freed, max_shell - nzon_shell)
        if n_to_shell > 0:
            nzon_shell += n_to_shell
            n_freed -= n_to_shell
            print(f"  Redistributed {n_to_shell} freed zones to shell")
        if n_freed > 0:
            print(f"  WARNING: {n_freed} zones could not be placed without "
                  f"violating min_dlogr — total will be {nzon_total - n_freed}")

    # Recalculate total
    n_total = nzon_shell + nzon_sw_buffer + nzon_wind + nzon_we_buffer + nzon_ejecta
    print(f"  Zone allocation: {nzon_shell} shell + {nzon_sw_buffer} SW buffer + "
          f"{nzon_wind} wind + {nzon_we_buffer} WE buffer + {nzon_ejecta} ejecta "
          f"= {n_total} total")

    # --- Build grids ---
    # Shell: split core/envelope grid for better envelope resolution.
    # The Dessart power-law shell has a flat core (r < r0) and a steep
    # envelope (r > r0) where the shock forms, radiation escapes, and
    # recombination occurs.  Concentrating zones in the envelope gives
    # 5x better resolution of the recombination front and shock.
    r0_shell = shell.get('r0', None)
    if r0_shell is not None and r0_shell > shell['r_inner']:
        shell_r = _build_optimized_shell_grid(
            r_min, r0_shell, r_shell_outer, nzon_shell, core_fraction=0.25,
            min_dlogr=min_dlogr)
        nzon_core = np.sum(shell_r < r0_shell)
        print(f"  Shell grid: {nzon_core} core + {nzon_shell - nzon_core} envelope "
              f"(split at r0={r0_shell:.2e} cm)")
    else:
        shell_r = np.logspace(log_r_min, log_r_shell, nzon_shell)
    # Interpolate shell profile (power-law density, homologous v, T ∝ rho^N_T)
    log_r_src = np.log10(shell['r'])
    log_rho_src = np.log10(np.maximum(shell['rho'], 1e-300))
    shell_rho = 10.0 ** np.interp(np.log10(shell_r), log_r_src, log_rho_src)
    shell_v = np.interp(np.log10(shell_r), log_r_src, shell['v'])
    shell_T = np.interp(np.log10(shell_r), log_r_src, shell['T'])

    # Conserve shell mass: the log-space interpolation of a steep power-law
    # (rho ∝ r^-8) onto fewer zones under-weights the dense core, and the
    # original 100-zone discrete grid also loses mass vs the analytical
    # formula. Renormalize to the analytical mass (exact by construction).
    if 'rho_0' in shell and 'r0' in shell:
        n_rho_val = shell.get('n_rho', 8.0)
        frac_v0 = shell['r0'] / shell['r_outer']
        factor = _shell_mass_factor(n_rho_val, frac_v0)
        M_shell_target = (4/3) * np.pi * shell['r0']**3 * shell['rho_0'] * factor
    else:
        # Fallback: integrate the original discrete profile
        r_lo_src = np.concatenate([[shell['r_inner']], shell['r'][:-1]])
        M_shell_target = ((4/3) * np.pi *
                          (shell['r']**3 - r_lo_src**3) * shell['rho']).sum()
    r_lo_shell = np.concatenate([[shell_r[0] * 0.5], shell_r[:-1]])
    M_shell_interp = ((4/3) * np.pi *
                      (shell_r**3 - r_lo_shell**3) * shell_rho).sum()
    if M_shell_interp > 0:
        shell_rho *= M_shell_target / M_shell_interp
        print(f"  Shell mass conservation: target={M_shell_target/MSUN:.5f} Msun, "
              f"interp={M_shell_interp/MSUN:.5f} Msun, "
              f"renorm factor={M_shell_target/M_shell_interp:.4f}")

    # Shell-wind buffer: log-spaced, log-linear density ramp
    sw_buf_r = np.logspace(log_r_shell, np.log10(r_sw_buffer_end), nzon_sw_buffer + 1)[1:]
    sw_log_rho = np.linspace(np.log10(max(rho_sw_start, 1e-300)),
                             np.log10(max(rho_sw_end, 1e-300)),
                             nzon_sw_buffer)
    sw_buf_rho = 10.0 ** sw_log_rho
    # Velocity: linear from shell v at outer edge to wind v
    sw_buf_v = np.linspace(shell['v'][-1], v_wind, nzon_sw_buffer)
    # Temperature: linear from shell T at outer edge to wind T
    sw_buf_T = np.linspace(shell['T'][-1], T_wind, nzon_sw_buffer)

    # Wind: log-spaced, density interpolated from the supplied profile
    # (steady 1/r^2 wind, or young wind + background for finite duration)
    wind_r = np.logspace(np.log10(r_sw_buffer_end), np.log10(r_we_buffer_start), nzon_wind + 1)[1:]
    wind_rho = wind_rho_at(wind_r)
    wind_v = np.full(nzon_wind, v_wind)
    wind_T = np.full(nzon_wind, T_wind)

    # Wind-ejecta buffer: log-spaced, log-linear density ramp
    we_buf_r = np.logspace(np.log10(r_we_buffer_start), log_r_ej, nzon_we_buffer + 1)[1:]
    we_log_rho = np.linspace(np.log10(max(rho_we_start, 1e-300)),
                             np.log10(max(rho_we_end, 1e-300)), nzon_we_buffer)
    we_buf_rho = 10.0 ** we_log_rho
    we_buf_v = np.linspace(v_wind, ej_v_sorted[0], nzon_we_buffer)
    we_buf_T = np.linspace(T_wind, max(ej_T_sorted[0], ejecta_T_floor), nzon_we_buffer)

    # Ejecta: log-spaced, interpolated from source
    ej_r_new = np.logspace(log_r_ej, log_r_max, nzon_ejecta + 1)[1:]
    log_r_ej_src = np.log10(ej_r_sorted)

    # PCHIP (piecewise cubic Hermite) interpolation: shape-preserving and
    # monotone where the data are monotone, so it cannot produce the Gibbs
    # wiggles that linear-in-log-space interpolation creates at steep
    # density/temperature gradients.  Density and temperature are
    # interpolated in log-space; velocity and composition in linear space.
    # Out-of-range radii are clamped to the source bounds (constant
    # extrapolation, matching the old fill_value behaviour).
    def log_interp(y, r_target):
        log_y = np.log10(np.maximum(y, 1e-300))
        ip = PchipInterpolator(log_r_ej_src, log_y, extrapolate=False)
        lr = np.clip(np.log10(r_target), log_r_ej_src[0], log_r_ej_src[-1])
        return 10.0 ** ip(lr)

    def lin_interp(y, r_target):
        ip = PchipInterpolator(log_r_ej_src, y, extrapolate=False)
        lr = np.clip(np.log10(r_target), log_r_ej_src[0], log_r_ej_src[-1])
        return ip(lr)

    ej_rho = log_interp(ej_rho_sorted, ej_r_new)
    ej_v = lin_interp(ej_v_sorted, ej_r_new)
    ej_T = np.maximum(log_interp(ej_T_sorted, ej_r_new), ejecta_T_floor)

    # Smooth ejecta T profile: the STELLA .swd T at late epochs has
    # numerical noise from the radiation-hydrodynamics evolution (zone-to-zone
    # scatter of ~20%). At 400d the ejecta is thermally equilibrated (radiation
    # transport has smoothed gradients), so a light moving-average smooth is
    # physically justified and gives a more robust initial condition.
    if len(ej_T) > 5:
        kernel = np.ones(5) / 5.0
        ej_T_smooth = np.convolve(ej_T, kernel, mode='same')
        ej_T_smooth[:2] = ej_T[:2]
        ej_T_smooth[-2:] = ej_T[-2:]
        ej_T = np.maximum(ej_T_smooth, ejecta_T_floor)

    # Smooth ejecta interpolation artifacts: any zone where the gradient > 16
    # gets replaced by the geometric mean of its neighbors. This removes sharp
    # spikes from the source SN model interpolation without affecting the
    # overall profile shape. Only applies to the low-density outer ejecta.
    for _ in range(3):  # 3 smoothing passes
        log_ej_rho = np.log10(np.maximum(ej_rho, 1e-300))
        log_ej_r = np.log10(ej_r_new)
        grads = np.abs(np.diff(log_ej_rho) / np.diff(log_ej_r))
        smoothed = ej_rho.copy()
        for i in range(1, len(ej_rho) - 1):
            g_left = grads[i - 1]
            g_right = grads[i]
            if g_left > 16 or g_right > 16:
                # Replace with geometric mean of neighbors
                smoothed[i] = np.sqrt(max(ej_rho[i - 1] * ej_rho[i + 1], 1e-300))
        ej_rho = smoothed

    # --- Combine all segments ---
    r_combined = np.concatenate([shell_r, sw_buf_r, wind_r, we_buf_r, ej_r_new])
    rho_combined = np.concatenate([shell_rho, sw_buf_rho, wind_rho, we_buf_rho, ej_rho])
    v_combined = np.concatenate([shell_v, sw_buf_v, wind_v, we_buf_v, ej_v])
    T_combined = np.concatenate([shell_T, sw_buf_T, wind_T, we_buf_T, ej_T])

    # --- Composition ---
    shell_xh = shell.get('fH', [0.7])[0]
    shell_yhe = shell.get('fHe', [0.3])[0]

    # Ejecta composition interpolated
    all_elem_keys = ['fH', 'fHe', 'fO', 'fSi', 'fCa', 'fFe']
    ej_comp = {}
    for key in all_elem_keys:
        if key in ejecta:
            ej_comp[key] = lin_interp(ejecta[key][ej_order], ej_r_new)
        else:
            ej_comp[key] = np.zeros(nzon_ejecta)

    # Ejecta composition at its inner edge (for WE buffer blending)
    ej_comp_inner = {}
    for key in all_elem_keys:
        if key in ejecta:
            ej_comp_inner[key] = float(ejecta[key][ej_order][0])
        else:
            ej_comp_inner[key] = 0.0

    # Smooth composition blends across the buffer zones.  Constant
    # composition in the buffers produces step-function metallicity jumps at
    # the buffer boundaries, causing sudden Rosseland mean opacity jumps
    # that can stall STELLA's radiation diffusion timestep.
    #
    # Smoothstep (3u^2-2u^3) weights are used so the blend has zero slope at
    # both ends.  The wind->ejecta blend is additionally extended over the
    # first n_ej_blend ejecta zones so the composition transition is
    # distributed across a wider region rather than compressed into the
    # narrow WE buffer.
    sw_u = np.linspace(0.0, 1.0, nzon_sw_buffer)
    sw_w = sw_u * sw_u * (3.0 - 2.0 * sw_u)

    n_ej_blend = min(15, nzon_ejecta // 4)
    n_we_blend = nzon_we_buffer + n_ej_blend
    we_u = np.linspace(0.0, 1.0, n_we_blend)
    we_w_full = we_u * we_u * (3.0 - 2.0 * we_u)
    we_w = we_w_full[:nzon_we_buffer]
    we_ej_w = we_w_full[nzon_we_buffer:]

    comp_combined = {}
    for key in all_elem_keys:
        # Shell composition (supports solar with metals)
        shell_val = float(shell.get(key, [0.0])[0])
        # Wind composition: read from the wind dict (solar or H/He-only,
        # set by create_wind's composition argument).  Falls back to H/He
        # only for backward compatibility if the wind dict lacks the key.
        if key in wind and len(wind[key]) > 0:
            wind_val = float(wind[key][0])
        elif key == 'fH':
            wind_val = 0.7
        elif key == 'fHe':
            wind_val = 0.3
        else:
            wind_val = 0.0

        ej_inner_val = ej_comp_inner[key]

        sw_buf_comp = shell_val * (1.0 - sw_w) + wind_val * sw_w
        we_buf_comp = wind_val * (1.0 - we_w) + ej_inner_val * we_w

        # Extend the wind->ejecta blend into the first n_ej_blend ejecta
        # zones, blending toward the local ejecta composition.
        ejc = ej_comp[key].copy()
        ejc[:n_ej_blend] = (wind_val * (1.0 - we_ej_w)
                            + ejc[:n_ej_blend] * we_ej_w)

        comp_combined[key] = np.concatenate([
            np.full(nzon_shell, shell_val),
            sw_buf_comp,   # shell -> wind blend
            np.full(nzon_wind, wind_val),
            we_buf_comp,   # wind -> ejecta blend (buffer part)
            ejc,           # ejecta (first n_ej_blend zones continue blend)
        ])

    # Renormalize composition
    comp_sum = sum(comp_combined.values())
    comp_sum = np.maximum(comp_sum, 1e-10)
    for key in comp_combined:
        comp_combined[key] /= comp_sum

    # --- Compute masses ---
    r_lo = np.concatenate([[r_combined[0] * 0.5], r_combined[:-1]])
    dm = (4/3) * np.pi * (r_combined**3 - r_lo**3) * rho_combined
    dm = np.maximum(dm, 0.0)
    M_total = np.sum(dm)

    # Component masses by segment
    idx_shell = nzon_shell
    idx_sw_buf = idx_shell + nzon_sw_buffer
    idx_wind = idx_sw_buf + nzon_wind
    idx_we_buf = idx_wind + nzon_we_buffer

    print(f"  Mass: shell={np.sum(dm[:idx_shell])/MSUN:.3f} Msun, "
          f"SW buffer={np.sum(dm[idx_shell:idx_sw_buf])/MSUN:.4f} Msun, "
          f"wind={np.sum(dm[idx_sw_buf:idx_wind])/MSUN:.4f} Msun, "
          f"WE buffer={np.sum(dm[idx_wind:idx_we_buf])/MSUN:.4f} Msun, "
          f"ejecta={np.sum(dm[idx_we_buf:])/MSUN:.3f} Msun, "
          f"total={M_total/MSUN:.3f} Msun")

    # --- Check density gradients ---
    log_r = np.log10(r_combined)
    log_rho_arr = np.log10(np.maximum(rho_combined, 1e-300))
    if len(log_r) > 1:
        gradients = np.abs(np.diff(log_rho_arr) / np.diff(log_r))
        max_grad = np.max(gradients)
        # Exclude the outermost ejecta edge (expected sharp drop)
        grad_interior = gradients[:-5]
        max_grad_interior = np.max(grad_interior) if len(grad_interior) > 0 else 0
        print(f"  Max |dlog rho/dlog r|: {max_grad:.1f} (interior: {max_grad_interior:.1f})")
        if max_grad_interior > 16:
            print(f"  WARNING: Interior density gradient {max_grad_interior:.1f} > 16 "
                  f"— may cause solver stall")
            worst_idx = np.argmax(grad_interior)
            print(f"    Worst at zone {worst_idx}: r={r_combined[worst_idx]:.3e} cm, "
                  f"rho={rho_combined[worst_idx]:.3e} -> {rho_combined[worst_idx+1]:.3e}")
        else:
            print(f"  Interior gradients OK (< 16, buffer zones effective)")

        # Print buffer zone gradients
        sw_buf_grads = gradients[idx_shell-1:idx_sw_buf]
        if len(sw_buf_grads) > 0:
            print(f"  SW buffer gradients: max={np.max(sw_buf_grads):.1f}, "
                  f"mean={np.mean(sw_buf_grads):.1f}")
        we_buf_grads = gradients[idx_wind-1:idx_we_buf]
        if len(we_buf_grads) > 0:
            print(f"  WE buffer gradients: max={np.max(we_buf_grads):.1f}, "
                  f"mean={np.mean(we_buf_grads):.1f}")

    return {
        'r': r_combined,
        'rho': rho_combined,
        'v': v_combined,
        'T': T_combined,
        'idx_ej': idx_we_buf,
        **comp_combined,
    }


def combine_shell_wind(shell, wind, nzon_total=500,
                       sw_transition_dex=0.5,
                       min_dlogr=0.005, max_dlogr=None):
    """Combine shell (inner) + wind (outer) with no old-SN ejecta.

    Type IIn geometry: an SN ejecta shell expanding into a steady 1/r^2
    CSM wind.  Identical construction to combine_shell_wind_ejecta but
    the wind simply extends to its outer radius instead of meeting a
    coasting ejecta component, so there is no wind-ejecta buffer.

    Parameters:
        sw_transition_dex: Target width of shell-wind buffer in dex of log10(r)
        max_dlogr:         Maximum allowed dlog10(r) per zone (see
                          combine_shell_wind_ejecta for the constraint logic)
    """
    buffer_max_dlogr = (min(max_dlogr, 2 * min_dlogr)
                        if max_dlogr is not None else 2 * min_dlogr)

    # Region boundaries
    r_min = shell['r_inner']
    r_shell_outer = shell['r_outer']
    r_max = wind['r_outer']

    log_r_min = np.log10(r_min)
    log_r_shell = np.log10(r_shell_outer)
    log_r_max = np.log10(r_max)

    # Wind properties
    mdot = wind['rho'][0] * 4 * np.pi * wind['r_inner']**2 * wind['v'][0]
    v_wind = wind['v'][0]
    T_wind = wind['T'][0]

    # Shell-wind buffer: same construction as combine_shell_wind_ejecta.
    # The buffer is capped at 1/3 of the shell->wind-outer span so the
    # wind always keeps positive radial extent.
    total_gap_dex = log_r_max - log_r_shell
    max_buffer_dex = total_gap_dex / 3.0
    sw_transition_dex = min(sw_transition_dex, max_buffer_dex)

    r_sw_buffer_end = r_shell_outer * (10 ** sw_transition_dex)
    rho_sw_start = shell['rho'][-1]
    rho_sw_end = mdot / (4 * np.pi * r_sw_buffer_end**2 * v_wind)
    sw_contrast = np.log10(max(rho_sw_start, 1e-300) / max(rho_sw_end, 1e-300))

    nzon_sw_buffer = max(20, int(sw_contrast / 16 * sw_transition_dex / 0.02))
    nzon_sw_buffer = min(nzon_sw_buffer, 60)
    nzon_sw_buffer = _cap_buffer_zones(nzon_sw_buffer, sw_transition_dex,
                                       min_dlogr=min_dlogr,
                                       max_dlogr=buffer_max_dlogr,
                                       label="SW buffer")

    print(f"  Shell-wind buffer: {nzon_sw_buffer} zones, "
          f"r=[{r_shell_outer:.3e}, {r_sw_buffer_end:.3e}] cm")
    print(f"    rho: {rho_sw_start:.3e} -> {rho_sw_end:.3e} g/cc "
          f"(contrast={sw_contrast:.1f} dex)")

    # --- Zone allocation ---
    # With no ejecta component the zones split between the shell (where
    # the forward/reverse shocks and recombination live) and the wind
    # (a smooth 1/r^2 profile that needs far fewer zones).  The total
    # must equal nzon_total: the STELLA binary is compiled with Mzon =
    # nzon and the job validates .hyd rows against it, so any zones the
    # wind cannot hold are redistributed to the shell.
    nzon_shell = max(100, int(nzon_total * 0.40))
    nzon_wind = nzon_total - nzon_shell - nzon_sw_buffer
    nzon_wind = max(60, nzon_wind)

    wind_span_dex = log_r_max - np.log10(r_sw_buffer_end)
    nzon_wind = _cap_buffer_zones(nzon_wind, wind_span_dex,
                                  min_dlogr=min_dlogr,
                                  max_dlogr=buffer_max_dlogr,
                                  label="Wind")

    # Redistribute freed zones to the shell.  The shell envelope has a
    # large physical span and is where the shocks live, so extra zones
    # are always useful there.  It is capped by its own min_dlogr limit;
    # anything still unplaceable is dropped with a warning (the pipeline
    # requires the total to equal nzon_total exactly).
    n_freed = nzon_total - nzon_shell - nzon_sw_buffer - nzon_wind
    if n_freed < 0:
        # Buffers grew past the allocation: take the excess back from
        # the shell (largest component) first, then the wind.
        n_excess = -n_freed
        n_from_shell = min(n_excess, max(0, nzon_shell - 100))
        nzon_shell -= n_from_shell
        n_excess -= n_from_shell
        n_from_wind = min(n_excess, max(0, nzon_wind - 60))
        nzon_wind -= n_from_wind
        n_excess -= n_from_wind
        if n_excess > 0:
            print(f"  WARNING: {n_excess} zones could not be removed without "
                  f"violating minimums — total will be {nzon_total + n_excess}")
        n_freed = 0
    if n_freed > 0:
        shell_span_dex = log_r_shell - log_r_min
        max_shell = max(int(shell_span_dex / min_dlogr), nzon_shell)
        n_to_shell = min(n_freed, max_shell - nzon_shell)
        if n_to_shell > 0:
            nzon_shell += n_to_shell
            n_freed -= n_to_shell
            print(f"  Redistributed {n_to_shell} freed zones to shell")
        if n_freed > 0:
            print(f"  WARNING: {n_freed} zones could not be placed without "
                  f"violating min_dlogr — total will be {nzon_total - n_freed}")

    n_total = nzon_shell + nzon_sw_buffer + nzon_wind
    print(f"  Zone allocation: {nzon_shell} shell + {nzon_sw_buffer} SW buffer + "
          f"{nzon_wind} wind = {n_total} total")

    # --- Build grids ---
    r0_shell = shell.get('r0', None)
    if r0_shell is not None and r0_shell > shell['r_inner']:
        shell_r = _build_optimized_shell_grid(
            r_min, r0_shell, r_shell_outer, nzon_shell, core_fraction=0.25,
            min_dlogr=min_dlogr)
    else:
        shell_r = np.logspace(log_r_min, log_r_shell, nzon_shell)
    log_r_src = np.log10(shell['r'])
    log_rho_src = np.log10(np.maximum(shell['rho'], 1e-300))
    shell_rho = 10.0 ** np.interp(np.log10(shell_r), log_r_src, log_rho_src)
    shell_v = np.interp(np.log10(shell_r), log_r_src, shell['v'])
    shell_T = np.interp(np.log10(shell_r), log_r_src, shell['T'])

    # Conserve shell mass (same renormalisation as the ejecta combiner)
    if 'rho_0' in shell and 'r0' in shell:
        n_rho_val = shell.get('n_rho', 8.0)
        frac_v0 = shell['r0'] / shell['r_outer']
        factor = _shell_mass_factor(n_rho_val, frac_v0)
        M_shell_target = (4/3) * np.pi * shell['r0']**3 * shell['rho_0'] * factor
    else:
        r_lo_src = np.concatenate([[shell['r_inner']], shell['r'][:-1]])
        M_shell_target = ((4/3) * np.pi *
                          (shell['r']**3 - r_lo_src**3) * shell['rho']).sum()
    r_lo_shell = np.concatenate([[shell_r[0] * 0.5], shell_r[:-1]])
    M_shell_interp = ((4/3) * np.pi *
                      (shell_r**3 - r_lo_shell**3) * shell_rho).sum()
    if M_shell_interp > 0:
        shell_rho *= M_shell_target / M_shell_interp
        print(f"  Shell mass conservation: target={M_shell_target/MSUN:.5f} Msun, "
              f"interp={M_shell_interp/MSUN:.5f} Msun, "
              f"renorm factor={M_shell_target/M_shell_interp:.4f}")

    # Shell-wind buffer: log-spaced, log-linear density ramp
    sw_buf_r = np.logspace(log_r_shell, np.log10(r_sw_buffer_end), nzon_sw_buffer + 1)[1:]
    sw_log_rho = np.linspace(np.log10(max(rho_sw_start, 1e-300)),
                             np.log10(max(rho_sw_end, 1e-300)),
                             nzon_sw_buffer)
    sw_buf_rho = 10.0 ** sw_log_rho
    sw_buf_v = np.linspace(shell['v'][-1], v_wind, nzon_sw_buffer)
    sw_buf_T = np.linspace(shell['T'][-1], T_wind, nzon_sw_buffer)

    # Wind: log-spaced, 1/r^2 density, constant velocity
    wind_r = np.logspace(np.log10(r_sw_buffer_end), log_r_max, nzon_wind + 1)[1:]
    wind_rho = mdot / (4 * np.pi * wind_r**2 * v_wind)
    wind_v = np.full(nzon_wind, v_wind)
    wind_T = np.full(nzon_wind, T_wind)

    # --- Combine all segments ---
    r_combined = np.concatenate([shell_r, sw_buf_r, wind_r])
    rho_combined = np.concatenate([shell_rho, sw_buf_rho, wind_rho])
    v_combined = np.concatenate([shell_v, sw_buf_v, wind_v])
    T_combined = np.concatenate([shell_T, sw_buf_T, wind_T])

    # --- Composition ---
    # Same smoothstep blend as combine_shell_wind_ejecta: a step in
    # metallicity at the shell-wind interface produces a Rosseland
    # opacity jump that can stall the STELLA diffusion timestep.
    all_elem_keys = ['fH', 'fHe', 'fO', 'fSi', 'fCa', 'fFe']
    sw_u = np.linspace(0.0, 1.0, nzon_sw_buffer)
    sw_w = sw_u * sw_u * (3.0 - 2.0 * sw_u)

    comp_combined = {}
    for key in all_elem_keys:
        shell_val = float(shell.get(key, [0.0])[0])
        if key in wind and len(wind[key]) > 0:
            wind_val = float(wind[key][0])
        elif key == 'fH':
            wind_val = 0.7
        elif key == 'fHe':
            wind_val = 0.3
        else:
            wind_val = 0.0

        sw_buf_comp = shell_val * (1.0 - sw_w) + wind_val * sw_w
        comp_combined[key] = np.concatenate([
            np.full(nzon_shell, shell_val),
            sw_buf_comp,
            np.full(nzon_wind, wind_val),
        ])

    # Renormalize composition
    comp_sum = sum(comp_combined.values())
    comp_sum = np.maximum(comp_sum, 1e-10)
    for key in comp_combined:
        comp_combined[key] /= comp_sum

    # --- Compute masses ---
    r_lo = np.concatenate([[r_combined[0] * 0.5], r_combined[:-1]])
    dm = (4/3) * np.pi * (r_combined**3 - r_lo**3) * rho_combined
    dm = np.maximum(dm, 0.0)
    M_total = np.sum(dm)

    idx_shell = nzon_shell
    idx_sw_buf = idx_shell + nzon_sw_buffer

    print(f"  Mass: shell={np.sum(dm[:idx_shell])/MSUN:.3f} Msun, "
          f"SW buffer={np.sum(dm[idx_shell:idx_sw_buf])/MSUN:.4f} Msun, "
          f"wind={np.sum(dm[idx_sw_buf:])/MSUN:.4f} Msun, "
          f"total={M_total/MSUN:.3f} Msun")

    # --- Check density gradients ---
    log_r = np.log10(r_combined)
    log_rho_arr = np.log10(np.maximum(rho_combined, 1e-300))
    if len(log_r) > 1:
        gradients = np.abs(np.diff(log_rho_arr) / np.diff(log_r))
        max_grad = np.max(gradients)
        grad_interior = gradients[:-5]
        max_grad_interior = np.max(grad_interior) if len(grad_interior) > 0 else 0
        print(f"  Max |dlog rho/dlog r|: {max_grad:.1f} (interior: {max_grad_interior:.1f})")
        if max_grad_interior > 16:
            print(f"  WARNING: Interior density gradient {max_grad_interior:.1f} > 16 "
                  f"— may cause solver stall")

    return {
        'r': r_combined,
        'rho': rho_combined,
        'v': v_combined,
        'T': T_combined,
        'idx_ej': None,
        **comp_combined,
    }


def combine_shell_wind_ejecta_density_matched(
    shell, ejecta, nzon_total=500,
    wind_velocity_kms=100.0, wind_T_K=2000.0,
    composition='solar', ejecta_T_floor=2000.0,
    sw_transition_dex=0.5, cd_transition_zones=5,
    v_blend_zones=0, min_dlogr=0.005, max_dlogr=None,
):
    """Density-matched wind-ejecta stitching with a thin linear-ramp CD.

    Derives Mdot so that rho_wind(r_interface) = rho_ejecta(r_interface),
    i.e. the wind and ejecta are in pressure/density balance at the contact
    radius.  This eliminates the density cliff at the wind-ejecta interface
    that motivates the extended log-linear buffer in the current combiner
    (combine_shell_wind_ejecta).

    Layout (inner -> outer):
      shell | shell-wind buffer | wind | CD ramp | ejecta

    - The **shell-wind buffer** is preserved (log-linear density ramp),
      because the shell edge is much denser than the wind at that radius --
      a genuine physical discontinuity that needs a buffer to keep
      dlog(rho)/dlog(r) < 16.

    - The **wind-ejecta interface** uses a thin linear ramp
      (cd_transition_zones zones) for velocity, temperature and
      composition.  The density is continuous by construction (matched at
      r_interface via the derived Mdot), so no density buffer is needed
      there.  The ramp is deliberately thin: the CD is a genuine physical
      discontinuity, so the transition spans only a few zones.

    The wind follows rho = Mdot / (4 pi r^2 v_wind) with the derived Mdot.
    The input wind_mdot is IGNORED -- Mdot is computed from the ejecta
    inner-edge density.

    Parameters:
        wind_velocity_kms:  Wind velocity (constant).  Mdot is derived
                            consistent with this v_wind.
        wind_T_K:           Wind temperature (isothermal).
        composition:        Wind composition ('solar' or 'hhe').
        sw_transition_dex:  Shell-wind buffer width in dex of log10(r).
        cd_transition_zones: Number of zones for the CD linear ramp.
    """
    from scipy.interpolate import PchipInterpolator

    # Maximum dlogr for buffers and adjacent regions: ensure no region's
    # dlogr exceeds 2x the minimum grid dlogr, preventing grid spacing
    # discontinuities at buffer boundaries that collapse the GEAR timestep.
    # The global max_dlogr (default 0.15) is a hard upper bound only; the
    # local smoothness constraint (2*min_dlogr) is tighter and is what
    # actually eliminates the factor-of-10 dlogr jumps seen at 80 days.
    buffer_max_dlogr = (min(max_dlogr, 2 * min_dlogr)
                        if max_dlogr is not None else 2 * min_dlogr)

    # Sort ejecta inner->outer
    ej_order = np.argsort(ejecta['r'])
    ej_r_sorted = ejecta['r'][ej_order]
    ej_rho_sorted = ejecta['rho'][ej_order]
    ej_v_sorted = ejecta['v'][ej_order]
    ej_T_sorted = ejecta['T'][ej_order]

    # Region boundaries
    r_min = shell['r_inner']
    r_shell_outer = shell['r_outer']
    r_ej_inner = ej_r_sorted[0] * 0.99  # small overlap, as in production
    r_max = ej_r_sorted[-1]

    log_r_min = np.log10(r_min)
    log_r_shell = np.log10(r_shell_outer)
    log_r_ej = np.log10(r_ej_inner)
    log_r_max = np.log10(r_max)

    v_wind = wind_velocity_kms * 1e5  # cm/s
    T_wind = wind_T_K

    # --- DERIVE Mdot from density matching at the contact radius ---
    log_r_ej_src = np.log10(ej_r_sorted)
    log_rho_ej_src = np.log10(np.maximum(ej_rho_sorted, 1e-300))
    rho_ej_at_interface = 10.0 ** np.interp(
        np.log10(r_ej_inner), log_r_ej_src, log_rho_ej_src)

    # Guard: if the ejecta density at the interface is effectively zero,
    # the derived Mdot will be ~0, producing a vacuum wind that causes
    # STELLA to crash (zero opacity, zero radiation coupling).
    if rho_ej_at_interface < 1e-30:
        raise ValueError(
            f"Ejecta density at the wind-ejecta interface is extremely small "
            f"({rho_ej_at_interface:.3e} g/cc < 1e-30 g/cc).\n"
            f"  A density-matched wind derived from this would be a vacuum, "
            f"causing STELLA to crash (zero opacity, zero radiation coupling).\n"
            f"  r_interface = {r_ej_inner:.3e} cm\n"
            f"  rho_ejecta(r_interface) = {rho_ej_at_interface:.3e} g/cc\n"
            f"Fixes:\n"
            f"  1. Reduce --ejecta-rho-min so denser outer zones are retained\n"
            f"  2. Reduce --ejecta-r-max-cm so the truncation does not remove "
            f"the relevant interface zones\n"
            f"  3. Use a standard wind (--density-matched-wind off) with an "
            f"explicit --wind-mdot"
        )

    # rho_wind(r) = Mdot / (4 pi r^2 v_wind)
    # => Mdot = 4 pi r^2 v_wind rho_ej_at_interface
    mdot_matched = 4.0 * np.pi * r_ej_inner**2 * v_wind * rho_ej_at_interface
    mdot_msun_yr = mdot_matched / (MSUN / YR_TO_S)

    print(f"\n  DENSITY-MATCHED wind-ejecta combiner:")
    print(f"    r_shell_outer = {r_shell_outer:.3e} cm")
    print(f"    r_ej_inner    = {r_ej_inner:.3e} cm")
    print(f"    rho_ejecta(r_interface) = {rho_ej_at_interface:.3e} g/cc")
    print(f"    => derived Mdot = {mdot_msun_yr:.4e} Msun/yr "
          f"(v_wind={wind_velocity_kms} km/s)")
    print(f"    CD ramp: {cd_transition_zones} zones")

    # --- Shell-wind buffer (same logic as combine_shell_wind_ejecta) ---
    total_gap_dex = log_r_ej - log_r_shell
    max_buffer_dex = total_gap_dex / 3.0
    sw_transition_dex = min(sw_transition_dex, max_buffer_dex)

    r_sw_buffer_end = r_shell_outer * (10 ** sw_transition_dex)
    rho_sw_start = shell['rho'][-1]
    rho_sw_end = mdot_matched / (4.0 * np.pi * r_sw_buffer_end**2 * v_wind)
    sw_contrast = np.log10(max(rho_sw_start, 1e-300) / max(rho_sw_end, 1e-300))
    sw_grad = sw_contrast / sw_transition_dex
    nzon_sw_buffer = max(20, int(sw_contrast / 16 * sw_transition_dex / 0.02))
    nzon_sw_buffer = min(nzon_sw_buffer, 60)
    # Enforce min_dlogr to prevent grid traps, and max_dlogr to prevent
    # grid discontinuities at buffer boundaries.
    nzon_sw_buffer = _cap_buffer_zones(nzon_sw_buffer, sw_transition_dex,
                                       min_dlogr=min_dlogr,
                                       max_dlogr=buffer_max_dlogr,
                                       label="SW buffer")

    print(f"    Shell-wind buffer: {nzon_sw_buffer} zones, "
          f"r=[{r_shell_outer:.3e}, {r_sw_buffer_end:.3e}] cm")
    print(f"      rho: {rho_sw_start:.3e} -> {rho_sw_end:.3e} g/cc "
          f"(contrast={sw_contrast:.1f} dex, est. gradient={sw_grad:.1f})")

    # --- Zone allocation ---
    # Optimized: shell gets 20% (up from 15%) for better shock/recombination
    # resolution.  Wind gets fewer (smooth 1/r^2).  Ejecta gets the remainder.
    nzon_shell = max(80, int(nzon_total * 0.20))
    nzon_ejecta = max(100, int(nzon_total * 0.55))
    nzon_wind = (nzon_total - nzon_shell - nzon_ejecta
                 - nzon_sw_buffer - cd_transition_zones)
    nzon_wind = max(40, nzon_wind)

    # Cap wind and ejecta zones to respect min_dlogr (prevent grid traps).
    # In the density-matched combiner, the wind grid includes the CD ramp
    # zones, so the span covers wind + CD and the zone count is their sum.
    wind_span_dex = log_r_ej - np.log10(r_sw_buffer_end)
    ej_span_dex = log_r_max - log_r_ej
    nzon_wind_capped = _cap_buffer_zones(nzon_wind + cd_transition_zones,
                                         wind_span_dex,
                                         min_dlogr=min_dlogr,
                                         max_dlogr=buffer_max_dlogr,
                                         label="Wind+CD")
    # Subtract CD zones back out (they're allocated separately)
    nzon_wind = max(nzon_wind_capped - cd_transition_zones, 10)
    nzon_ejecta = _cap_buffer_zones(nzon_ejecta, ej_span_dex,
                                    min_dlogr=min_dlogr,
                                    max_dlogr=buffer_max_dlogr,
                                    label="Ejecta")

    # Redistribute freed zones to the shell (envelope-limited, see
    # combine_shell_wind_ejecta for details).
    n_freed = (nzon_total - nzon_shell - nzon_sw_buffer - nzon_wind
               - cd_transition_zones - nzon_ejecta)
    if n_freed < 0:
        # Buffer zones grew (max_dlogr enforcement) and the total exceeds
        # nzon_total.  Reduce the ejecta first, then wind, then shell.
        n_excess = -n_freed
        n_from_ej = min(n_excess, max(0, nzon_ejecta - 100))
        nzon_ejecta -= n_from_ej
        n_excess -= n_from_ej
        n_from_wind = min(n_excess, max(0, nzon_wind - 10))
        nzon_wind -= n_from_wind
        n_excess -= n_from_wind
        n_from_shell = min(n_excess, max(0, nzon_shell - 80))
        nzon_shell -= n_from_shell
        n_excess -= n_from_shell
        if n_from_ej > 0 or n_from_wind > 0 or n_from_shell > 0:
            print(f"    Buffer max_dlogr enforcement grew buffers by {-n_freed} zones; "
                  f"reduced ejecta by {n_from_ej}, wind by {n_from_wind}, "
                  f"shell by {n_from_shell}")
        if n_excess > 0:
            print(f"    WARNING: {n_excess} zones could not be removed without "
                  f"violating minimums — total will be {nzon_total + n_excess}")
        n_freed = 0
    if n_freed > 0:
        shell_span_dex = np.log10(r_shell_outer) - np.log10(r_min)
        max_shell = max(int(shell_span_dex / min_dlogr), nzon_shell)
        n_to_shell = min(n_freed, max_shell - nzon_shell)
        if n_to_shell > 0:
            nzon_shell += n_to_shell
            n_freed -= n_to_shell
            print(f"    Redistributed {n_to_shell} freed zones to shell")
        if n_freed > 0:
            print(f"    WARNING: {n_freed} zones could not be placed without "
                  f"violating min_dlogr — total will be {nzon_total - n_freed}")

    n_total = (nzon_shell + nzon_sw_buffer + nzon_wind
               + cd_transition_zones + nzon_ejecta)
    print(f"    Zone allocation: {nzon_shell} shell + {nzon_sw_buffer} SW buffer "
          f"+ {nzon_wind} wind + {cd_transition_zones} CD-ramp + "
          f"{nzon_ejecta} ejecta = {n_total}")

    # --- Shell grid (split core/envelope for better envelope resolution) ---
    r0_shell = shell.get('r0', None)
    if r0_shell is not None and r0_shell > shell['r_inner']:
        shell_r = _build_optimized_shell_grid(
            r_min, r0_shell, r_shell_outer, nzon_shell, core_fraction=0.25,
            min_dlogr=min_dlogr)
        nzon_core = np.sum(shell_r < r0_shell)
        print(f"    Shell grid: {nzon_core} core + {nzon_shell - nzon_core} envelope "
              f"(split at r0={r0_shell:.2e} cm)")
    else:
        shell_r = np.logspace(log_r_min, log_r_shell, nzon_shell)
    log_r_src = np.log10(shell['r'])
    log_rho_src = np.log10(np.maximum(shell['rho'], 1e-300))
    shell_rho = 10.0 ** np.interp(np.log10(shell_r), log_r_src, log_rho_src)
    shell_v = np.interp(np.log10(shell_r), log_r_src, shell['v'])
    shell_T = np.interp(np.log10(shell_r), log_r_src, shell['T'])

    # Shell mass conservation (same as combine_shell_wind_ejecta)
    if 'rho_0' in shell and 'r0' in shell:
        n_rho_val = shell.get('n_rho', 8.0)
        frac_v0 = shell['r0'] / shell['r_outer']
        factor = _shell_mass_factor(n_rho_val, frac_v0)
        M_shell_target = (4/3) * np.pi * shell['r0']**3 * shell['rho_0'] * factor
    else:
        r_lo_src = np.concatenate([[shell['r_inner']], shell['r'][:-1]])
        M_shell_target = ((4/3) * np.pi *
                          (shell['r']**3 - r_lo_src**3) * shell['rho']).sum()
    r_lo_shell = np.concatenate([[shell_r[0] * 0.5], shell_r[:-1]])
    M_shell_interp = ((4/3) * np.pi *
                      (shell_r**3 - r_lo_shell**3) * shell_rho).sum()
    if M_shell_interp > 0:
        shell_rho *= M_shell_target / M_shell_interp
        print(f"    Shell mass conservation: target={M_shell_target/MSUN:.5f} Msun, "
              f"interp={M_shell_interp/MSUN:.5f} Msun, "
              f"renorm factor={M_shell_target/M_shell_interp:.4f}")

    # --- Shell-wind buffer (log-linear density ramp) ---
    sw_buf_r = np.logspace(log_r_shell, np.log10(r_sw_buffer_end),
                           nzon_sw_buffer + 1)[1:]
    sw_log_rho = np.linspace(np.log10(max(rho_sw_start, 1e-300)),
                             np.log10(max(rho_sw_end, 1e-300)),
                             nzon_sw_buffer)
    sw_buf_rho = 10.0 ** sw_log_rho
    sw_buf_v = np.linspace(shell['v'][-1], v_wind, nzon_sw_buffer)
    sw_buf_T = np.linspace(shell['T'][-1], T_wind, nzon_sw_buffer)

    # --- Wind grid (1/r^2 with derived Mdot) ---
    # Wind spans from r_sw_buffer_end to r_ej_inner, with the last
    # cd_transition_zones points being the CD transition ramp.
    wind_r = np.logspace(np.log10(r_sw_buffer_end), log_r_ej,
                         nzon_wind + cd_transition_zones + 1)[1:]
    wind_r_pure = wind_r[:nzon_wind]
    cd_r = wind_r[nzon_wind:]

    wind_rho_pure = mdot_matched / (4.0 * np.pi * wind_r_pure**2 * v_wind)
    wind_v_pure = np.full(nzon_wind, v_wind)
    wind_T_pure = np.full(nzon_wind, T_wind)

    # --- Contact-discontinuity transition (thin linear ramp) ---
    # Density: matched by construction -- wind 1/r^2 law continues through
    # the CD points, and at r_ej_inner rho_wind = rho_ejecta by Mdot choice.
    cd_rho = mdot_matched / (4.0 * np.pi * cd_r**2 * v_wind)

    # Velocity: FLAT at v_wind through the CD.  The velocity blend from
    # the wind's constant v_wind to the ejecta's homologous v(r) profile
    # is done over a transition region at the start of the ejecta (see
    # below), not as a sharp jump at the CD.
    cd_v = np.full(cd_transition_zones, v_wind)

    # Temperature: linear ramp from T_wind to ejecta T (floored) across the
    # thin CD segment — the temperature jump at the contact discontinuity
    # is spread over cd_transition_zones zones instead of a sigmoid.
    cd_T = np.linspace(T_wind, max(ej_T_sorted[0], ejecta_T_floor),
                       cd_transition_zones)

    # --- Ejecta grid (log-spaced, interpolated from source) ---
    ej_r_new = np.logspace(log_r_ej, log_r_max, nzon_ejecta + 1)[1:]

    # PCHIP interpolation (see combine_shell_ejecta for rationale): monotone,
    # no Gibbs wiggles.  Out-of-range radii are clamped to the source bounds.
    def log_interp(y, r_target):
        log_y = np.log10(np.maximum(y, 1e-300))
        ip = PchipInterpolator(log_r_ej_src, log_y, extrapolate=False)
        lr = np.clip(np.log10(r_target), log_r_ej_src[0], log_r_ej_src[-1])
        return 10.0 ** ip(lr)

    def lin_interp(y, r_target):
        ip = PchipInterpolator(log_r_ej_src, y, extrapolate=False)
        lr = np.clip(np.log10(r_target), log_r_ej_src[0], log_r_ej_src[-1])
        return ip(lr)

    ej_rho = log_interp(ej_rho_sorted, ej_r_new)
    ej_v = lin_interp(ej_v_sorted, ej_r_new)
    ej_T = np.maximum(log_interp(ej_T_sorted, ej_r_new), ejecta_T_floor)

    # --- Homologous velocity blend at the ejecta inner edge ---
    # The wind has constant v_wind.  The ejecta has a homologous v(r)
    # profile (v ~ r/t for free expansion).
    #
    # CMFGEN now supports non-monotonic velocity profiles, so the velocity
    # blend is OPTIONAL and disabled by default (v_blend_zones=0).  A sharp
    # velocity jump at the Contact Discontinuity is fine for STELLA (which
    # only needs smooth density gradients) and is correctly handled by
    # CMFGEN's Sobolev optical-depth integration even where dv/dr < 0 or
    # changes sign.  Blending is only applied if explicitly requested for
    # specific hydrodynamic smoothing.
    #
    # When enabled, the blend weight w goes from 0 (at the inner edge,
    # v = v_wind) to 1 (at the outer edge of the blend region,
    # v = v_ej_natural), using a linear ramp across the blend zones.
    if v_blend_zones > 1:
        v_blend_zones = min(v_blend_zones, nzon_ejecta // 3)
        blend_w = np.linspace(0.0, 1.0, v_blend_zones)
        # Pad to full ejecta length (weight=1 beyond the blend region)
        full_w = np.concatenate([blend_w, np.ones(nzon_ejecta - v_blend_zones)])
        ej_v = v_wind * (1.0 - full_w) + ej_v * full_w
        print(f"    Velocity blend: {v_blend_zones} zones (SMOOTHED), "
              f"v_wind={v_wind/1e5:.1f} km/s -> v_ej(r) homologous")
    else:
        print(f"    Velocity blend: SKIPPED (sharp CD preserved for CMFGEN)")
        # ej_v remains the natural homologous v(r) from the source ejecta.
        # The wind has v_wind. The jump happens exactly at the interface.

    # Light T smoothing (same as combine_shell_wind_ejecta)
    if len(ej_T) > 5:
        kernel = np.ones(5) / 5.0
        ej_T_smooth = np.convolve(ej_T, kernel, mode='same')
        ej_T_smooth[:2] = ej_T[:2]
        ej_T_smooth[-2:] = ej_T[-2:]
        ej_T = np.maximum(ej_T_smooth, ejecta_T_floor)

    # Light rho smoothing (same as combine_shell_wind_ejecta)
    for _ in range(3):
        log_ej_rho = np.log10(np.maximum(ej_rho, 1e-300))
        log_ej_r = np.log10(ej_r_new)
        grads = np.abs(np.diff(log_ej_rho) / np.diff(log_ej_r))
        smoothed = ej_rho.copy()
        for i in range(1, len(ej_rho) - 1):
            g_left = grads[i - 1]
            g_right = grads[i]
            if g_left > 16 or g_right > 16:
                smoothed[i] = np.sqrt(max(ej_rho[i - 1] * ej_rho[i + 1], 1e-300))
        ej_rho = smoothed

    # --- Force exact density continuity at the contact discontinuity ---
    # The analytical 1/r^2 wind density at the CD can differ from the
    # *interpolated* ejecta density by 1-2%.  Even a tiny density jump at
    # the CD triggers STELLA's artificial viscosity oscillations, so pin
    # the first ejecta zone to the last CD-zone density and apply light
    # geometric smoothing over the next few ejecta zones to avoid a kink
    # in the density gradient.
    if len(ej_rho) > 0 and len(cd_rho) > 0:
        ej_rho[0] = cd_rho[-1]
        for i in range(1, min(5, len(ej_rho))):
            ej_rho[i] = np.sqrt(max(ej_rho[i - 1] * ej_rho[i], 1e-300))

    # --- Combine all segments ---
    r_combined = np.concatenate([shell_r, sw_buf_r, wind_r_pure, cd_r, ej_r_new])
    rho_combined = np.concatenate([shell_rho, sw_buf_rho, wind_rho_pure, cd_rho, ej_rho])
    v_combined = np.concatenate([shell_v, sw_buf_v, wind_v_pure, cd_v, ej_v])
    T_combined = np.concatenate([shell_T, sw_buf_T, wind_T_pure, cd_T, ej_T])

    # --- Composition ---
    all_elem_keys = ['fH', 'fHe', 'fO', 'fSi', 'fCa', 'fFe']
    ej_comp = {}
    for key in all_elem_keys:
        if key in ejecta:
            ej_comp[key] = lin_interp(ejecta[key][ej_order], ej_r_new)
        else:
            ej_comp[key] = np.zeros(nzon_ejecta)

    shell_comp_vals = {k: float(shell.get(k, [0.0])[0]) for k in all_elem_keys}

    if composition == 'solar':
        sc = SOLAR_COMPOSITION
        wind_comp_vals = {k: sc[k] for k in all_elem_keys}
    else:
        wind_comp_vals = {'fH': 0.7, 'fHe': 0.3,
                          'fO': 0.0, 'fSi': 0.0, 'fCa': 0.0, 'fFe': 0.0}

    # CD composition: linear ramp wind -> ejecta inner-edge composition
    # across the thin CD segment
    ej_comp_inner = {}
    for key in all_elem_keys:
        if key in ejecta:
            ej_comp_inner[key] = float(ejecta[key][ej_order][0])
        else:
            ej_comp_inner[key] = 0.0

    cd_w = np.linspace(0.0, 1.0, cd_transition_zones)

    # Smooth smoothstep blend across the SW buffer (shell -> wind
    # composition), avoiding a step-function metallicity jump that would
    # cause a sudden Rosseland mean opacity jump in STELLA's radiation
    # diffusion.  Smoothstep gives zero slope at both ends of the buffer.
    sw_u = np.linspace(0.0, 1.0, nzon_sw_buffer)
    sw_w = sw_u * sw_u * (3.0 - 2.0 * sw_u)

    comp_combined = {}
    for key in all_elem_keys:
        sv = shell_comp_vals[key]
        wv = wind_comp_vals[key]
        ev = ej_comp_inner[key]
        sw_buf_comp = sv * (1.0 - sw_w) + wv * sw_w
        cd_comp = wv * (1.0 - cd_w) + ev * cd_w
        comp_combined[key] = np.concatenate([
            np.full(nzon_shell, sv),
            sw_buf_comp,                   # SW buffer: shell -> wind blend
            np.full(nzon_wind, wv),
            cd_comp,
            ej_comp[key],
        ])

    # Renormalize composition
    comp_sum = sum(comp_combined.values())
    comp_sum = np.maximum(comp_sum, 1e-10)
    for key in comp_combined:
        comp_combined[key] /= comp_sum

    # --- Compute masses ---
    r_lo = np.concatenate([[r_combined[0] * 0.5], r_combined[:-1]])
    dm = (4/3) * np.pi * (r_combined**3 - r_lo**3) * rho_combined
    dm = np.maximum(dm, 0.0)
    M_total = np.sum(dm)

    idx_shell = nzon_shell
    idx_sw_buf = idx_shell + nzon_sw_buffer
    idx_wind = idx_sw_buf + nzon_wind
    idx_cd = idx_wind + cd_transition_zones

    print(f"    Mass: shell={np.sum(dm[:idx_shell])/MSUN:.3f} Msun, "
          f"SW buf={np.sum(dm[idx_shell:idx_sw_buf])/MSUN:.4f} Msun, "
          f"wind={np.sum(dm[idx_sw_buf:idx_wind])/MSUN:.4f} Msun, "
          f"CD={np.sum(dm[idx_wind:idx_cd])/MSUN:.6f} Msun, "
          f"ejecta={np.sum(dm[idx_cd:])/MSUN:.3f} Msun, "
          f"total={M_total/MSUN:.3f} Msun")

    # --- Check density gradients ---
    log_r = np.log10(r_combined)
    log_rho_arr = np.log10(np.maximum(rho_combined, 1e-300))
    if len(log_r) > 1:
        gradients = np.abs(np.diff(log_rho_arr) / np.diff(log_r))
        max_grad = np.max(gradients)
        grad_interior = gradients[:-5]
        max_grad_interior = np.max(grad_interior) if len(grad_interior) > 0 else 0
        print(f"    Max |dlog rho/dlog r|: {max_grad:.1f} "
              f"(interior: {max_grad_interior:.1f})")
        if max_grad_interior > 16:
            print(f"    WARNING: Interior density gradient {max_grad_interior:.1f} > 16 "
                  f"— may cause solver stall")
            worst_idx = np.argmax(grad_interior)
            print(f"      Worst at zone {worst_idx}: r={r_combined[worst_idx]:.3e} cm, "
                  f"rho={rho_combined[worst_idx]:.3e} -> {rho_combined[worst_idx+1]:.3e}")
        else:
            print(f"      Interior gradients OK (< 16)")

        sw_buf_grads = gradients[idx_shell-1:idx_sw_buf]
        if len(sw_buf_grads) > 0:
            print(f"    SW buffer gradients: max={np.max(sw_buf_grads):.1f}, "
                  f"mean={np.mean(sw_buf_grads):.1f}")
        cd_grads = gradients[idx_wind-1:idx_cd]
        if len(cd_grads) > 0:
            print(f"    CD ramp gradients: max={np.max(cd_grads):.2f}, "
                  f"mean={np.mean(cd_grads):.2f}")

    return {
        'r': r_combined,
        'rho': rho_combined,
        'v': v_combined,
        'T': T_combined,
        'idx_ej': idx_cd,
        **comp_combined,
    }


def write_x_hyd(fname, g, mass_cut_msun=0.0, t_start_days=0.0):
    """Write STELLA .hyd file for SN_B model.

    STELLA .hyd format (8 columns per zone row):
      km  dMr(Msun)  r(cm)  rho(g/cc)  T(K)  u(cm/s)  Mr(Msun)  0
    Header: t_start  nzon  mass_cut  rcen  rho_cen
    """
    r = g['r']; rho = g['rho']; T = g['T']; v = g['v']
    nzon = len(r)

    rcen = r[0] * 0.5
    rho_cen = rho[0]

    r_lo = np.concatenate([[rcen], r[:-1]])
    dm = (4/3) * np.pi * (r**3 - r_lo**3) * rho
    dm = np.maximum(dm, 0.0)
    M_r = np.cumsum(dm)

    with open(fname, 'w') as f:
        # timetobegin is always 0: the clock starts at the shell/wind
        # interaction; t_start_days only sets the pre-expanded structure.
        f.write(f"{0.0:12.3E}{nzon:6d}{mass_cut_msun:13.5E}{rcen:13.5E}{rho_cen:13.5E}\n")
        for km in range(nzon):
            dMr = dm[km] / MSUN
            Mr = M_r[km] / MSUN + mass_cut_msun
            # km  dMr  rHyd  rhoHyd  TpHyd  uHyd  aMr  dum
            f.write(f" {km+1:4d}{dMr:12.4E}{r[km]:15.7E}{rho[km]:12.4E}{T[km]:12.4E}"
                    f"{v[km]:12.4E}{Mr:12.4E}{0.0:12.4E}\n")
    print(f"Written {fname}  ({nzon} zones, M_total={M_r[-1]/MSUN:.4f} Msun, "
          f"t_start={t_start_days}d)")


def write_x_abn(fname, g):
    """Write STELLA .abn file for SN_B model.

    Column order (20 columns):
        km  0  0  0  H  He  C  N  O  Ne  Na  Mg  Al  Si  S  Ar  Ca  FePeak  Ni58  Ni56

    Only H, He, O, Si, Ca, Fe are written (solar composition).
    All other elements (C, N, Ne, Na, Mg, Al, S, Ar, Ni) are set to 0.
    """
    n = len(g['r'])
    with open(fname, 'w') as f:
        for km in range(n):
            f.write(f"{km+1:4d}{0.0:10.3E}{0.0:10.3E}{0.0:10.3E}"
                    f"{g['fH'][km]:10.3E}{g['fHe'][km]:10.3E}"
                    f"{0.0:10.3E}"                         # C
                    f"{0.0:10.3E}"                         # N
                    f"{g['fO'][km]:10.3E}"
                    f"{0.0:10.3E}"                         # Ne
                    f"{0.0:10.3E}"                         # Na
                    f"{0.0:10.3E}"                         # Mg
                    f"{0.0:10.3E}"                         # Al
                    f"{g['fSi'][km]:10.3E}"
                    f"{0.0:10.3E}"                         # S
                    f"{0.0:10.3E}"                         # Ar
                    f"{g['fCa'][km]:10.3E}{g['fFe'][km]:10.3E}"
                    f"{0.0:10.3E}"                         # Ni58
                    f"{0.0:10.3E}\n")                      # Ni56
    print(f"Written {fname}  ({n} zones)")


def write_x_xni(fname, nzon):
    """Write a zero .xni file (no radioactive decay in SN_B model)."""
    n_floats = 2 + 2 * nzon
    floats = [0.0] * n_floats
    floats[0] = 5.605194e-42
    floats[1] = -1.417409e-34
    data = struct.pack('f' * n_floats, *floats)
    with open(fname, 'wb') as f:
        f.write(struct.pack('i', len(data)))
        f.write(data)
        f.write(struct.pack('i', len(data)))
    print(f"Written {fname}  ({len(data)} bytes, zero Ni56)")


def write_x_dat(fname, model_name, target_days, snapshot_day,
                explosion_energy=0.0, explosion_mass=0.0, explosion_time=1.0,
                nstmax=2000000, eps=0.003, hmax=3600.0, meth=2,
                aq=10.0, bq=1.0, iout=1, snapshot_interval=5.0,
                fit_tau=5.0, temp_floor=0.03, floor_rad=0.01, hmin=1e-14,
                t_start_days=0.0, r_max_cm=0.0, v_max_cm_s=0.0,
                wacc_r=1.0, wacc_v=1.0, wacc_t=0.5, wacc_rad=0.5,
                tau_tol=1.3, disable_scattering=False, tcura_days=None):
    """Write STELLA .dat control file for SN_B model.

    The SN_B model is velocity-driven: EKO=0 and EBurst=0 are enforced
    here because the shell's launch velocity and the coasting old ejecta
    are already encoded in the .hyd file.  EKO>0 would overwrite that
    velocity column; a nonzero EBurst would add a thermal bomb on top of
    the moving shell.  Uses high AQ for non-monotonic velocity
    (shell 5000 km/s -> wind 600 km/s -> ejecta).

    t_start_days: simulation starts at this PROPER_T (shell pre-expanded).
    target_days: desired OBS.TIME (observer-frame days since explosion) to reach.
    snapshot_day: desired OBS.TIME for the CMFGEN/SuperLite snapshot.
    r_max_cm: outermost radius of the .hyd profile (cm). Used to compute the
              light-travel-time offset (r_max/c) that STELLA subtracts from
              PROPER_T to get OBS.TIME. TcurB (PROPER_T stop) is set to
              target_days + r_max/c so that OBS.TIME reaches target_days.
    v_max_cm_s: outermost velocity (cm/s). Used to account for the increasing
               light-travel-time as the ejecta expands. Without this correction,
               TcurB underestimates the needed PROPER_T and the final OBS.TIME
               falls short of target_days.
    wacc_r: Wacc(R) - radius accuracy weight (default: 1.0).
    wacc_v: Wacc(V) - velocity accuracy weight (default: 1.0).
    wacc_t: Wacc(T) - temperature accuracy weight (default: 0.5 for stability).
    wacc_rad: Wacc(RADIAT) - radiation accuracy weight (default: 0.5 for stability).
    tau_tol: TauTol - Ncnd changes when tau outside [FitTau/TauTol, FitTau*TauTol]
             (default: 1.3; use 1.1 for scattering runs to avoid large Ncnd jumps).
    disable_scattering: if True, set SCAT=F (default: False, SCAT=T).
    """
    c_light = 2.998e10  # cm/s
    light_time_days = r_max_cm / c_light / 86400.0 if r_max_cm > 0 else 0.0

    # TcurA is the STELLA clock value at run start; defaults to t_start_days
    # (the shell's coasting time) but can be overridden independently to test
    # the clock convention without changing the pre-expanded profile.
    t_start = tcura_days if tcura_days is not None else t_start_days
    # TcurB: PROPER_T at which OBS.TIME = target_days
    # OBS.TIME = PROPER_T - r_max/c, so PROPER_T = target_days + r_max/c
    #
    # CRITICAL: r_max increases as the ejecta expands.  If we only use the
    # initial light-travel-time, the run stops before OBS.TIME reaches
    # target_days.  Solve self-consistently:
    #   r_max(t) = r_max_0 + v_max * (t - t_start) * 86400
    #   light_time(t) = light_time_0 + v_max * (t - t_start) / c
    #   TcurB = target_days + light_time(TcurB)
    #   TcurB * (1 - v_max/c) = target_days + light_time_0 - v_max*t_start/c
    if v_max_cm_s > 0 and v_max_cm_s < c_light and t_start > 0:
        denom = 1.0 - v_max_cm_s / c_light
        numer = target_days + light_time_days - v_max_cm_s * t_start / c_light
        t_end = numer / denom
        light_time_final = light_time_days + v_max_cm_s * (t_end - t_start) / c_light
        print(f"  Light-travel-time correction: initial={light_time_days:.2f} d, "
              f"final={light_time_final:.2f} d (v_max={v_max_cm_s/1e5:.0f} km/s)")
        print(f"  TcurB = {t_end:.2f} d (corrected) vs {target_days + light_time_days:.2f} d (uncorrected)")
    else:
        t_end = target_days + light_time_days

    # Output epochs: STELLA compares these with OBS.TIME (t_ob = PROPER_T - r/c)
    # So the epochs should be in OBS.TIME days (relative to explosion).
    KOMAX = 80
    early = [0.5, 1.0, 2.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0]
    if snapshot_day not in early:
        early.append(snapshot_day)
    # Late epochs only if target_days > 35
    if target_days > 35:
        late = np.linspace(35.0, target_days,
                           max(5, int((target_days - 35) / snapshot_interval))).tolist()
    else:
        late = []
    epochs = sorted(set(early + late))
    # Filter to only epochs <= target_days (don't request output beyond sim end)
    epochs = [e for e in epochs if e <= target_days][:KOMAX]

    if light_time_days > 0:
        print(f"  Light-travel-time offset: {light_time_days:.2f} days (r_max={r_max_cm:.3e} cm)")
        print(f"  TcurB (PROPER_T stop) = {t_end:.2f} days (OBS.TIME target={target_days}d + offset={light_time_days:.2f}d)")

    with open(fname, 'w') as f:
        f.write(f"   EPS   (STIFF ACCURACY)   HMIN (MIN.STEP(SEC)) HMAX (MAX.STEP(SEC))\n")
        f.write(f"   {eps:.4g}                    {hmin:.4g}               {hmax:.4g}\n")
        f.write(f"   METH  ADAMS(1)-GEAR(2)-GEAR+BGH(3)   JSTART  (0-START_ELSE_ORDER) MAXORD KNadap\n")
        f.write(f"    {meth}                        0                            4        -4\n")
        f.write(f"   abs(NSTA)>abs(NSTB)=>reads last   TcurA     TcurB\n")
        f.write(f"   -9000     -1                      {t_start:.4f}   {t_end:.4f}\n")
        f.write(f"   NSTMAX  NDebug  NOUT    IOUT (0 -full <0 print zones)   MBATCH\n")
        nout_val = len(epochs)
        f.write(f"   {nstmax}   100000  {nout_val}      {iout}                             900\n")
        f.write(f"   AMNI (CONTAM. MASS)     XMNI (56NI MASS)\n")
        f.write(f"   0.0d0                  0.0d0\n")
        f.write(f"   AMHT(Solar)  EBurst(1e50) tBurst(s)  tstart -- Heated Core, Energy & time\n")
        f.write(f"   {explosion_mass:.4g}d0      {explosion_energy:.4g}d0      {explosion_time:.4g}d0      0.0d0\n")
        f.write(f"   EKO (kin.energy,1e50) mass fract.tri. u profile  us (out +1)\n")
        # EKO=0 hardcoded: preserve the .hyd velocity field (SN_B is
        # velocity-driven; EKO>0 would overwrite it with a triangular profile)
        f.write(f"   0.0d0                      0.d0  9.d-1  1.d0          +1.\n")
        f.write(f"   THRMAT   CRAP     CONV     EDTM     CHNCND  Givdtl\n")
        f.write(f"   1.D-30   1.D0     F        T        T       T\n")
        f.write(f"   FLOOR(R) FLOOR(V) FLOOR(T) FLOOR(RADIAT)    FLOOR(UCONV)\n")
        f.write(f"   1.D-12   1.D-03   {temp_floor:.4g}D+00   {floor_rad:.4g}           1.D-04\n")
        f.write(f"   Wacc(R)  Wacc(V)  Wacc(T)  Wacc(RADIAT)\n")
        f.write(f"   {wacc_r:.2g}d-0    {wacc_v:.2g}d-0    {wacc_t:.2g}D-00   {wacc_rad:.2g}D-00\n")
        scat_flag = "F" if disable_scattering else "T"
        f.write(f"   FitTau   TAUTOL   Rvis   AQ    BQ       DRT     NRT     SCAT\n")
        f.write(f"   {fit_tau:.1f}D0   {tau_tol:.1f}D0     1.D0   {aq:.1f}D0  {bq:.1f}D+00   1.      1       {scat_flag}\n")
        f.write(f"    NTO then    DATA TO\n")
        f.write(f"    {len(epochs)}\n")
        for ep in epochs:
            f.write(f"     {ep:.4f}\n")
    print(f"Written {fname}  ({len(epochs)} output epochs, OBS.TIME target={target_days} days)")


def write_x_input_profile(fname, g):
    """Write 12-column input_profile for convert_profile_to_stella.py."""
    r = g['r']; rho = g['rho']; T = g['T']

    fH = g['fH']; fHe = g['fHe']
    fO = g['fO']; fSi = g['fSi']; fFe = g['fFe']
    fCa = g.get('fCa', np.zeros_like(fH))

    mu = 1.0 / (fH/1.0 + fHe/4.0 + fO/16.0 + fSi/28.0 + fCa/40.0 + fFe/56.0)
    mu = np.maximum(mu, 0.5)
    p_gas = rho * K_B * T / (mu * M_H)
    p_rad = A_RAD * T**4 / 3.0
    p = p_gas + p_rad

    fO_profile = fO  # No C to merge

    header = ("# input_profile generated by prepare_SN_B_model.py\n"
              "# columns: r_cm rho_cgs p_cgs u_cms T_K fx_H fx_He fx_O fx_Si fx_Ca fx_Fe\n"
              "# Inner: H/He shell; Middle: wind; Outer: 400-day SN ejecta\n")
    with open(fname, 'w') as f:
        f.write(header)
        for i in range(len(r)):
            f.write(f"{r[i]:.8e} {rho[i]:.8e} {p[i]:.8e} {g['v'][i]:.8e} "
                    f"{T[i]:.8e} "
                    f"{fH[i]:.8e} {fHe[i]:.8e} {fO_profile[i]:.8e} "
                    f"{fSi[i]:.8e} {fCa[i]:.8e} {fFe[i]:.8e} {0.0:.8e}\n")
    print(f"Written {fname}  ({len(r)} zones)")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                 formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('source_dir', nargs='?', default=None,
                        help='Source SN model directory (e.g. SN_12C_CSM2). '
                             'Not required with --no-ejecta.')
    parser.add_argument('output_dir', help='SN_B model output directory (e.g. SN_B_v2)')
    parser.add_argument('--model-name', default=None,
                        help='Model name for output files (default: derived from '
                             'output_dir name). Use this to match the slurm '
                             'template MODEL_NAME when debug mode appends -debug.')
    parser.add_argument('--shell-mass', type=float, default=1.0,
                        help='H/He shell mass in Msun (default: 1.0)')
    parser.add_argument('--shell-velocity', type=float, default=5000,
                        help='Shell launch velocity in km/s (default: 5000)')
    parser.add_argument('--shell-xh', type=float, default=0.7,
                        help='H mass fraction in shell (default: 0.7)')
    parser.add_argument('--shell-yhe', type=float, default=0.3,
                        help='He mass fraction in shell (default: 0.3)')
    parser.add_argument('--shell-composition', type=str, default='hhe',
                        choices=['hhe', 'solar'],
                        help='Shell composition: hhe (H/He only, uses --shell-xh/yhe) '
                             'or solar (Asplund+2009, X=0.738 Y=0.249 Z=0.013). '
                             'For a solar-composition shell + wind + old SN ejecta '
                             'scenario with maximum STELLA stability, use '
                             '--shell-composition solar together with '
                             '--density-matched-wind: this gives the shell realistic '
                             'metal content matching the wind/ejecta and eliminates '
                             'the density cliff at the wind-ejecta interface.')
    parser.add_argument('--n-rho', type=float, default=8.0,
                        help='Shell density power-law exponent N_rho (Dessart 2010jl: 8)')
    parser.add_argument('--v0', type=float, default=3000.0,
                        help='Core/envelope transition velocity V0 in km/s (Dessart: 3000)')
    parser.add_argument('--n-T', type=float, default=0.4,
                        help='Temperature power-law exponent N_T (Dessart: 0.4)')
    parser.add_argument('--T-interface', type=float, default=5000.0,
                        help='Shell temperature at outer edge / wind interface (K)')
    parser.add_argument('--T-max', type=float, default=65000.0,
                        help='Maximum shell temperature in core (Dessart: 65000 K)')
    parser.add_argument('--t-start', type=float, default=10.0,
                        help='STELLA start time in days (shell pre-expanded). '
                             'Sets r_outer = v_shell * t_start. Must be > 0. '
                             '(default: 10.0, Dessart 2010jl: 11.6)')
    parser.add_argument('--tcura-days', type=float, default=None,
                        help='TcurA: STELLA clock start time in days written to '
                             'the .dat file. Default: use --t-start.')
    parser.add_argument('--ejecta-T-floor', type=float, default=1000.0,
                        help='Ejecta temperature floor in K (400d-old Ni56-free: ~300 R)')
    parser.add_argument('--ejecta-rho-min', type=float, default=0.0,
                        help='Truncate outer ejecta zones with rho below this threshold '
                             '(g/cc). Reduces r_max and the light-travel-time offset that '
                             'makes STELLA spend walltime evolving unobservable outer '
                             'ejecta. E.g. 1e-20 truncates zones with rho < 1e-20 g/cc.')
    parser.add_argument('--ejecta-r-max-cm', type=float, default=0.0,
                        help='Truncate outer ejecta zones with r above this threshold (cm). '
                             'Useful for quick interactive tests: the outer coasting ejecta '
                             'at r > r_max does not participate in the interaction, so '
                             'truncating it cuts the light-travel-time offset. E.g. 1e16 '
                             'cuts the offset from ~22d to ~4d for a 400d-old ejecta.')
    parser.add_argument('--ejecta-v-max-kms', type=float, default=0.0,
                        help='Truncate outer ejecta zones with v above this threshold '
                             '(km/s). Caps the outermost ejecta velocity, reducing '
                             'r_max and the light-travel-time offset. E.g. 10000 '
                             'keeps only material slower than 10000 km/s.')
    parser.add_argument('--wind-mdot', type=float, default=1e-3,
                        help='Wind mass-loss rate in Msun/yr (default: 1e-3)')
    parser.add_argument('--wind-velocity', type=float, default=600,
                        help='Wind velocity in km/s (default: 600)')
    parser.add_argument('--wind-T', type=float, default=5000,
                        help='Wind temperature in K (default: 5000)')
    parser.add_argument('--target-days', type=float, default=45.0,
                        help='STELLA simulation duration in days (default: 45)')
    parser.add_argument('--snapshot-day', type=float, default=30.0,
                        help='CMFGEN snapshot day after shell launch (default: 30)')
    parser.add_argument('--nzon', type=int, default=500,
                        help='Total number of STELLA zones (default: 500)')
    parser.add_argument('--source-age', type=float, default=400.0,
                        help='Age of source ejecta at shell launch in days (default: 400)')
    parser.add_argument('--explosion-energy', type=float, default=0.0,
                        help='EBurst in 1e50 erg. KEEP 0 for SN_B: the model is '
                             'velocity-driven (shell mass+velocity in .hyd); a '
                             'nonzero value deposits a thermal bomb on top.')
    parser.add_argument('--nstmax', type=int, default=2000000)
    parser.add_argument('--aq', type=float, default=10.0,
                        help='Artificial viscosity (default: 10.0 for strong interaction)')
    parser.add_argument('--fit-tau', type=float, default=5.0,
                        help='FitTau parameter (default: 5.0 per AGENTS.md)')
    parser.add_argument('--meth', type=int, default=3,
                        help='Integration method (default: 3=GEAR+BGH, 2=GEAR, 1=Adams)')
    parser.add_argument('--temp-floor', type=float, default=0.03,
                        help='Temperature floor in keV (default: 0.03)')
    parser.add_argument('--floor-rad', type=float, default=0.01,
                        help='Radiation floor (default: 0.01)')
    parser.add_argument('--hmin', type=float, default=1e-14,
                        help='Minimum timestep in seconds (default: 1e-14)')
    parser.add_argument('--hmax', type=float, default=3600.0,
                        help='Maximum timestep in seconds (default: 3600)')
    parser.add_argument('--wacc-r', type=float, default=1.0,
                        help='Wacc(R) radius accuracy weight (default: 1.0)')
    parser.add_argument('--wacc-v', type=float, default=1.0,
                        help='Wacc(V) velocity accuracy weight (default: 1.0)')
    parser.add_argument('--wacc-t', type=float, default=0.5,
                        help='Wacc(T) temperature accuracy weight (default: 0.5 for stability)')
    parser.add_argument('--wacc-rad', type=float, default=0.5,
                        help='Wacc(RADIAT) radiation accuracy weight (default: 0.5 for stability)')
    parser.add_argument('--tau-tol', type=float, default=1.3,
                        help='TauTol: Ncnd changes when tau outside [FitTau/TauTol, FitTau*TauTol] '
                             '(default: 1.3; use 1.1 for scattering runs)')
    parser.add_argument('--iout', type=int, default=10,
                        help='Output frequency: 1=every step, 10=every 10 steps (default: 10)')
    parser.add_argument('--disable-scattering', action='store_true',
                        help='Disable scattering (SCAT=F, default: SCAT=T)')
    parser.add_argument('--min-dlogr', type=float, default=0.005,
                        help='Minimum allowed dlog10(r) between adjacent zones '
                             '(default: 0.005, dr/r ~ 0.011). If violated, a '
                             'CRITICAL warning is printed (or ValueError with '
                             '--strict-grid). Buffer zones are auto-capped to '
                             'respect this limit.')
    parser.add_argument('--max-dlogr', type=float, default=0.15,
                        help='Maximum allowed dlog10(r) between adjacent zones '
                             '(default: 0.15, dr/r ~ 0.41). If violated, a '
                             'warning is printed about possible shock '
                             'under-resolution.')
    parser.add_argument('--strict-grid', action='store_true',
                        help='If set, raise ValueError and halt when min_dlogr '
                             'is violated, instead of printing a warning.')
    parser.add_argument('--sw-transition-dex', type=float, default=0.5,
                        help='Shell-wind buffer width in dex of log10(r) (default: 0.5)')
    parser.add_argument('--sw-transition-k', type=float, default=3.0,
                        help='Unused (kept for CLI compatibility)')
    parser.add_argument('--we-transition-dex', type=float, default=0.5,
                        help='Wind-ejecta buffer width in dex of log10(r) (default: 0.5)')
    parser.add_argument('--we-transition-k', type=float, default=3.0,
                        help='Unused (kept for CLI compatibility)')
    parser.add_argument('--no-wind', action='store_true',
                        help='Skip the wind component: shell + ejecta only (no wind in between)')
    parser.add_argument('--no-ejecta', action='store_true',
                        help='Skip the old-SN ejecta component: shell + wind only '
                             '(Type IIn geometry — an ejecta shell expanding into a '
                             'steady 1/r^2 CSM wind). source_dir is not required. '
                             'The wind outer radius is set by --wind-r-outer-cm, or '
                             'defaults to v_wind * (source_age + t_start).')
    parser.add_argument('--wind-r-outer-cm', type=float, default=0.0,
                        help='Outer radius of the CSM wind in cm (--no-ejecta mode). '
                             'If 0, defaults to v_wind * (source_age + t_start), '
                             'i.e. source_age is interpreted as the wind duration '
                             'before shell launch.')
    parser.add_argument('--wind-duration-days', type=float, default=0.0,
                        help='Finite-age wind (standard shell+wind+ejecta mode). '
                             'If >0, truncate the connecting wind at '
                             'r_lead = v_wind * wind_duration_days and fill the '
                             'remaining gap to the ejecta with a tenuous '
                             'background wind (--wind-mdot-bg). Use this when '
                             'the mass-loss episode began a known number of '
                             'days before shell launch (e.g. at SN_A), so the '
                             'wind cannot have reached the ejecta inner edge. '
                             'Default 0 fills the whole shell-ejecta gap with '
                             'the same steady Mdot (fossil-CSM assumption).')
    parser.add_argument('--wind-mdot-bg', type=float, default=1e-7,
                        help='Background wind Mdot in Msun/yr used beyond the '
                             'truncated wind leading edge (--wind-duration-days). '
                             'Represents the older, tenuous CSM the young wind '
                             'runs into.')
    parser.add_argument('--immediate-contact', action='store_true',
                        help='Shell and ejecta in direct contact (no buffer zone). '
                             'Computes self-consistent t_start so shell outer edge '
                             'meets ejecta inner edge. Requires --no-wind.')
    parser.add_argument('--se-transition-dex', type=float, default=0.5,
                        help='Shell-ejecta buffer width in dex of log10(r) (no-wind mode, default: 0.5)')
    parser.add_argument('--density-matched-wind', action='store_true',
                        help='Use density-matched wind-ejecta stitching: derive Mdot so '
                             'rho_wind(r_interface) = rho_ejecta(r_interface), with a '
                             'thin linear ramp for v/T/composition at the contact '
                             'discontinuity (instead of a log-linear density buffer). '
                             'The shell-wind buffer is preserved. The input --wind-mdot '
                             'is IGNORED in this mode (Mdot is derived). '
                             'RECOMMENDED for a solar-composition shell + wind + old '
                             'SN ejecta scenario: use --shell-composition solar '
                             '--density-matched-wind together for maximum STELLA '
                             'stability (eliminates the density cliff at the '
                             'wind-ejecta interface that a standard 1/r^2 wind '
                             'produces).')
    parser.add_argument('--cd-transition-zones', type=int, default=5,
                        help='Number of zones for the CD linear ramp '
                             '(density-matched-wind mode, default: 5)')
    parser.add_argument('--v-blend-zones', type=int, default=0,
                        help='Number of ejecta zones over which to homologously blend '
                             'velocity from v_wind to the ejecta natural v(r) profile '
                             '(density-matched-wind mode, default: 0). The wind and CD '
                             'have flat v_wind; the blend happens at the start of the '
                             'ejecta region. Default 0 keeps a sharp velocity jump at '
                             'the contact discontinuity, which CMFGEN now handles '
                             'natively via non-monotonic velocity support. Set >1 only '
                             'for specific hydrodynamic smoothing tests.')
    parser.add_argument('--force-homologous-velocity', action='store_true',
                        help='Force the SN_A ejecta velocity to pure homologous expansion '
                             '(v = r/t). Removes residual non-homologous structure from '
                             'the source STELLA model (e.g. CSM deceleration, shocks). '
                             'Physically appropriate for an SN_A ejecta that is hundreds '
                             'of days old and should be in free expansion.')
    parser.add_argument('--no-smooth-ejecta', action='store_true',
                        help='Disable Savitzky-Golay smoothing of the SN_A ejecta '
                             'profile. By default, smoothing is applied to remove '
                             'numerical noise (zone-to-zone scatter) from the STELLA '
                             'RHD evolution. At hundreds of days, the ejecta should '
                             'be in smooth homologous expansion; the scatter is '
                             'numerical artefact that can cause STELLA solver stalls.')
    parser.add_argument('--smooth-window', type=int, default=15,
                        help='Savitzky-Golay filter window size for ejecta smoothing '
                             '(odd, default 15). Larger = more smoothing.')
    parser.add_argument('--smooth-polyorder', type=int, default=3,
                        help='Savitzky-Golay polynomial order for ejecta smoothing '
                             '(default 3). Lower = more smoothing. 3 preserves '
                             'power-law gradients.')
    parser.add_argument('--comp-smooth-window', type=int, default=7,
                        help='Moving-average window (zones) used to soften sharp '
                             'composition interfaces in the source ejecta '
                             '(default 7; 0 or 1 disables). Applied to each '
                             'species mass fraction, then renormalized.')
    args = parser.parse_args()

    if not args.no_ejecta and not args.source_dir:
        parser.error('source_dir is required unless --no-ejecta is set')
    if args.no_ejecta and args.no_wind:
        parser.error('--no-ejecta with --no-wind would leave a shell-only '
                     'model, which is not supported')
    if args.no_ejecta and args.density_matched_wind:
        parser.error('--density-matched-wind derives Mdot from the ejecta '
                     'interface density; it is incompatible with --no-ejecta')

    source_dir = Path(args.source_dir) if args.source_dir else None
    output_dir = Path(args.output_dir)
    model_name = args.model_name if args.model_name else output_dir.name

    # Create output directories
    stella_dir = output_dir / 'stella'
    modmake_dir = stella_dir / 'modmake'
    plots_dir = output_dir / 'plots'
    modmake_dir.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n=== Preparing X initial model ===")
    print(f"  Source directory: {source_dir}")
    print(f"  Output directory: {output_dir}")
    print(f"  Shell: M={args.shell_mass} Msun, v={args.shell_velocity} km/s, "
          f"X={args.shell_xh}, Y={args.shell_yhe}, composition={args.shell_composition}")
    if args.no_ejecta:
        print(f"  Mode: shell + wind only (Type IIn, --no-ejecta)")
        print(f"  Wind: Mdot={args.wind_mdot} Msun/yr, v={args.wind_velocity} km/s")
    elif args.no_wind:
        print(f"  Wind: SKIPPED (--no-wind mode)")
    elif args.density_matched_wind:
        print(f"  Wind: DENSITY-MATCHED (Mdot derived, v={args.wind_velocity} km/s, "
              f"CD ramp zones={args.cd_transition_zones})")
    else:
        print(f"  Wind: Mdot={args.wind_mdot} Msun/yr, v={args.wind_velocity} km/s")
    print(f"  Source age: {args.source_age} days")
    print(f"  STELLA target: {args.target_days} days, CMFGEN snapshot: {args.snapshot_day} days")

    if args.no_ejecta:
        ejecta = None
    else:
        # Step 1: Read source ejecta and homologously evolve to source_age + t_start_days
        # The SN_A ejecta is source_age days old at SN_B explosion (t=0).
        # STELLA starts at t_start_days later, so the ejecta is (source_age + t_start_days) days old.
        effective_age = args.source_age + args.t_start
        print(f"\n--- Step 1: Reading source ejecta at {effective_age} days "
              f"({args.source_age}d at SN_B + {args.t_start}d t_start) ---")
        ejecta = read_source_ejecta(source_dir, source_age=args.source_age,
                                    t_start_days=args.t_start,
                                    force_homologous_velocity=args.force_homologous_velocity,
                                    smooth_ejecta=not args.no_smooth_ejecta,
                                    smooth_window=args.smooth_window,
                                    smooth_polyorder=args.smooth_polyorder,
                                    comp_window=args.comp_smooth_window)

        print(f"  Ejecta: {len(ejecta['r'])} zones, "
              f"r=[{ejecta['r'].min():.3e}, {ejecta['r'].max():.3e}] cm")
        print(f"  Ejecta: v=[{ejecta['v'].min():.1f}, {ejecta['v'].max():.1f}] km/s")
        print(f"  Ejecta: rho=[{ejecta['rho'].min():.3e}, {ejecta['rho'].max():.3e}] g/cc")
        print(f"  Ejecta: T=[{ejecta['T'].min():.3e}, {ejecta['T'].max():.3e}] K")

        # Convert ejecta velocities from km/s to cm/s
        ejecta['v'] = ejecta['v'] * 1e5  # km/s -> cm/s

    # Truncate outer ejecta zones with density below threshold
    # This reduces r_max and the light-travel-time offset (r_max/c) that
    # makes STELLA spend walltime evolving unobservable outer ejecta.
    # The outer ejecta zones at 400d have negligible density (1e-24 g/cc)
    # and negligible mass (<1e-4 Msun), but their large radius creates a
    # huge light-travel-time offset that delays the OBS.TIME by ~50 days.
    if ejecta is not None and args.ejecta_rho_min > 0:
        mask = ejecta['rho'] >= args.ejecta_rho_min
        n_before = len(ejecta['r'])
        n_truncated = n_before - mask.sum()
        if n_truncated > 0:
            r_max_before = ejecta['r'].max()
            mass_lost = 0.0  # approximate
            for key in ejecta:
                if isinstance(ejecta[key], np.ndarray) and len(ejecta[key]) == n_before:
                    ejecta[key] = ejecta[key][mask]
            r_max_after = ejecta['r'].max()
            c_light = 2.998e10
            lt_before = r_max_before / c_light / 86400
            lt_after = r_max_after / c_light / 86400
            print(f"\n  Ejecta truncation (rho_min={args.ejecta_rho_min:.1e} g/cc):")
            print(f"    Zones: {n_before} -> {mask.sum()} (removed {n_truncated} outer zones)")
            print(f"    r_max: {r_max_before:.3e} -> {r_max_after:.3e} cm")
            print(f"    Light-travel time: {lt_before:.2f} -> {lt_after:.2f} days")

    # Radius-based truncation: keep only ejecta zones with r <= ejecta_r_max_cm.
    # This is useful for quick interactive tests where the light-travel-time
    # offset from the extended outer ejecta (r~1e17 at 400d) would force STELLA
    # to run for ~20+ days internally. The outer coasting ejecta at r > r_max
    # does not participate in the shell-wind-ejecta interaction, so truncating
    # it is physically justified for interaction-focused snapshots.
    if ejecta is not None and args.ejecta_r_max_cm > 0:
        mask = ejecta['r'] <= args.ejecta_r_max_cm
        n_before = len(ejecta['r'])
        n_truncated = n_before - mask.sum()
        if n_truncated > 0:
            r_max_before = ejecta['r'].max()
            for key in ejecta:
                if isinstance(ejecta[key], np.ndarray) and len(ejecta[key]) == n_before:
                    ejecta[key] = ejecta[key][mask]
            r_max_after = ejecta['r'].max()
            c_light = 2.998e10
            lt_before = r_max_before / c_light / 86400
            lt_after = r_max_after / c_light / 86400
            print(f"\n  Ejecta radius truncation (r_max={args.ejecta_r_max_cm:.1e} cm):")
            print(f"    Zones: {n_before} -> {mask.sum()} (removed {n_truncated} outer zones)")
            print(f"    r_max: {r_max_before:.3e} -> {r_max_after:.3e} cm")
            print(f"    Light-travel time: {lt_before:.2f} -> {lt_after:.2f} days")

    # Velocity-based truncation: keep only ejecta zones with
    # v <= ejecta_v_max_kms.  Caps the outermost ejecta velocity, which
    # reduces r_max, the light-travel-time offset, and the velocity range
    # the CMFGEN grid must cover.  The truncated outer material is fast,
    # low-density coasting ejecta that does not participate in the
    # interaction on the timescales modelled.
    if ejecta is not None and args.ejecta_v_max_kms > 0:
        v_max_cms = args.ejecta_v_max_kms * 1e5  # v is cm/s here
        mask = ejecta['v'] <= v_max_cms
        n_before = len(ejecta['r'])
        n_truncated = n_before - mask.sum()
        if n_truncated > 0:
            r_max_before = ejecta['r'].max()
            for key in ejecta:
                if isinstance(ejecta[key], np.ndarray) and len(ejecta[key]) == n_before:
                    ejecta[key] = ejecta[key][mask]
            r_max_after = ejecta['r'].max()
            v_max_after = ejecta['v'].max()
            c_light = 2.998e10
            lt_before = r_max_before / c_light / 86400
            lt_after = r_max_after / c_light / 86400
            print(f"\n  Ejecta velocity truncation (v_max={args.ejecta_v_max_kms:.0f} km/s):")
            print(f"    Zones: {n_before} -> {mask.sum()} (removed {n_truncated} outer zones)")
            print(f"    r_max: {r_max_before:.3e} -> {r_max_after:.3e} cm")
            print(f"    v_max: {v_max_cms/1e5:.0f} -> {v_max_after/1e5:.0f} km/s")
            print(f"    Light-travel time: {lt_before:.2f} -> {lt_after:.2f} days")

    # Step 2: Create H/He shell (Dessart 2010jl power-law structure)
    # For no-wind mode: place shell outer edge directly at ejecta inner edge
    # (immediate contact — the shell is the bullet, the ejecta is the target).
    # This eliminates the large gap that would otherwise exist because the
    # shell (launched at v_shell from center) hasn't physically reached the
    # ejecta inner edge at t_start days. Instead, we start the simulation at
    # the moment of contact, like Dessart+2015.
    r_ej_inner = (ejecta['r'].min() * 0.99) if ejecta is not None else None
    if args.no_wind:
        print(f"\n--- Step 2: Creating H/He shell (immediate contact with ejecta) ---")
        print(f"  Shell outer edge set to ejecta inner edge: {r_ej_inner:.3e} cm")
        print(f"  (instead of v_shell * t_start = {args.shell_velocity*1e5*args.t_start*86400:.3e} cm)")
        shell = create_hhe_shell(args.shell_mass, args.shell_velocity,
                                 args.shell_xh, args.shell_yhe,
                                 n_rho=args.n_rho, v0_kms=args.v0,
                                 n_T=args.n_T, T_interface_K=args.T_interface,
                                 T_max_K=args.T_max, t_start_days=args.t_start,
                                 r_outer_override=r_ej_inner,
                                 composition=args.shell_composition)
    else:
        print(f"\n--- Step 2: Creating H/He shell (Dessart 2010jl structure) ---")
        shell = create_hhe_shell(args.shell_mass, args.shell_velocity,
                                 args.shell_xh, args.shell_yhe,
                                 n_rho=args.n_rho, v0_kms=args.v0,
                                 n_T=args.n_T, T_interface_K=args.T_interface,
                                 T_max_K=args.T_max, t_start_days=args.t_start,
                                 composition=args.shell_composition)

    # ------------------------------------------------------------------
    # Grid inversion safety check (after Step 2, before Step 3).
    #
    # If --t-start is too large, the shell's outer radius
    # (v_shell * t_start) can physically overtake the 400-day-old ejecta's
    # inner radius before the simulation begins.  This produces
    # wind_r_inner > wind_r_outer, an inverted grid, negative zone masses,
    # and an immediate STELLA crash.
    #
    # In --no-wind / --immediate-contact mode the shell outer edge is
    # overridden to the ejecta inner edge, so the overlap is intentional
    # and this check is skipped.
    # ------------------------------------------------------------------
    if not args.no_wind and not args.no_ejecta:
        r_shell_outer = shell['r_outer']
        r_ej_inner_check = ejecta['r'].min()
        if r_shell_outer >= r_ej_inner_check:
            v_shell_cms = args.shell_velocity * 1e5
            t_start_s = args.t_start * 86400.0
            raise ValueError(
                f"Grid inversion: shell outer radius has overtaken the ejecta "
                f"inner radius.\n"
                f"  r_shell_outer = {r_shell_outer:.3e} cm "
                f"(v_shell={args.shell_velocity} km/s * t_start={args.t_start} d)\n"
                f"  r_ejecta_inner = {r_ej_inner_check:.3e} cm "
                f"(source_age={args.source_age} d ejecta)\n"
                f"  The shell would need to travel {r_ej_inner_check / v_shell_cms / 86400.0:.2f} d "
                f"to reach the ejecta, but t_start={args.t_start} d has already "
                f"expanded it past that point.\n"
                f"Fixes:\n"
                f"  1. Reduce --t-start to < {r_ej_inner_check / v_shell_cms / 86400.0:.2f} d\n"
                f"  2. Reduce --shell-velocity (currently {args.shell_velocity} km/s)\n"
                f"  3. Use --density-matched-wind (adapts wind to the available gap)\n"
                f"  4. Use --no-wind --immediate-contact (shell directly meets ejecta)"
            )
        gap_dex = np.log10(r_ej_inner_check) - np.log10(r_shell_outer)
        print(f"\n  Grid inversion check: PASSED")
        print(f"    r_shell_outer = {r_shell_outer:.3e} cm < "
              f"r_ejecta_inner = {r_ej_inner_check:.3e} cm "
              f"(gap = {gap_dex:.2f} dex)")

    # Step 3: Create connecting wind (or skip if --no-wind)
    if args.no_ejecta:
        print(f"\n--- Step 3: Creating CSM wind (shell + wind, no ejecta) ---")
        # Wind outer radius: explicit --wind-r-outer-cm, else the radius
        # reached by wind emitted since (source_age + t_start) days, i.e.
        # source_age is the wind duration before shell launch.
        if args.wind_r_outer_cm > 0:
            r_wind_outer = args.wind_r_outer_cm
        else:
            t_wind_s = (args.source_age + args.t_start) * 86400.0
            r_wind_outer = args.wind_velocity * 1e5 * t_wind_s
        wind_r_inner = shell['r_outer'] * 1.01
        if wind_r_inner >= r_wind_outer:
            raise ValueError(
                f"Wind grid inversion (--no-ejecta): wind_r_inner >= r_wind_outer.\n"
                f"  wind_r_inner = {wind_r_inner:.3e} cm "
                f"(shell r_outer * 1.01 = {shell['r_outer']:.3e} * 1.01)\n"
                f"  r_wind_outer = {r_wind_outer:.3e} cm\n"
                f"Fixes:\n"
                f"  1. Increase --wind-r-outer-cm (currently {args.wind_r_outer_cm:.3e})\n"
                f"  2. Reduce --t-start so the shell has not expanded as far\n"
                f"  3. Increase --source-age (wind duration) or --wind-velocity"
            )
        print(f"  Wind outer radius: {r_wind_outer:.3e} cm "
              f"({r_wind_outer / (args.wind_velocity * 1e5 * 86400.0):.1f} d of wind)")
        wind = create_wind(args.wind_mdot, args.wind_velocity, args.wind_T,
                           wind_r_inner, r_wind_outer,
                           composition=args.shell_composition)
        print(f"\n--- Step 4: Combining shell + wind (no ejecta) ---")
        combined = combine_shell_wind(shell, wind, nzon_total=args.nzon,
                                      sw_transition_dex=args.sw_transition_dex,
                                      min_dlogr=args.min_dlogr,
                                      max_dlogr=args.max_dlogr)
    elif args.no_wind:
        print(f"\n--- Step 3: SKIPPING wind (--no-wind mode, immediate contact) ---")
        combined = combine_shell_ejecta(shell, ejecta, nzon_total=args.nzon,
                                        se_transition_dex=args.se_transition_dex,
                                        immediate_contact=True,
                                        ejecta_T_floor=args.ejecta_T_floor,
                                        min_dlogr=args.min_dlogr,
                                        max_dlogr=args.max_dlogr)
    elif args.density_matched_wind:
        print(f"\n--- Step 3: DENSITY-MATCHED wind-ejecta stitching ---")
        print(f"    (Mdot derived from rho_ejecta(r_interface), "
              f"CD ramp: {args.cd_transition_zones} zones)")
        combined = combine_shell_wind_ejecta_density_matched(
            shell, ejecta, nzon_total=args.nzon,
            wind_velocity_kms=args.wind_velocity,
            wind_T_K=args.wind_T,
            composition=args.shell_composition,
            ejecta_T_floor=args.ejecta_T_floor,
            sw_transition_dex=args.sw_transition_dex,
            cd_transition_zones=args.cd_transition_zones,
            v_blend_zones=args.v_blend_zones,
            min_dlogr=args.min_dlogr,
            max_dlogr=args.max_dlogr)
    else:
        print(f"\n--- Step 3: Creating connecting wind ---")
        # Wind spans from shell outer edge to ejecta inner edge
        wind_r_inner = shell['r_outer'] * 1.01
        wind_r_outer = ejecta['r'].min() * 0.99

        # ------------------------------------------------------------------
        # Wind grid boundary validation.
        #
        # Even with the Step 2 grid-inversion check, a very small gap
        # (shell nearly touching ejecta) can produce wind_r_inner >=
        # wind_r_outer after the 1% margins are applied, inverting the
        # wind grid.  Catch this explicitly with a clear error.
        # ------------------------------------------------------------------
        if wind_r_inner >= wind_r_outer:
            raise ValueError(
                f"Wind grid inversion: wind_r_inner >= wind_r_outer.\n"
                f"  wind_r_inner = {wind_r_inner:.3e} cm "
                f"(shell r_outer * 1.01 = {shell['r_outer']:.3e} * 1.01)\n"
                f"  wind_r_outer = {wind_r_outer:.3e} cm "
                f"(ejecta r.min * 0.99 = {ejecta['r'].min():.3e} * 0.99)\n"
                f"  The shell-ejecta gap is too small for a standard 1/r^2 wind "
                f"(the 1% margins on each side overlap).\n"
                f"Fixes:\n"
                f"  1. Reduce --t-start so the shell has not expanded as far\n"
                f"  2. Use --density-matched-wind (handles small gaps gracefully)\n"
                f"  3. Use --no-wind --immediate-contact (no wind needed)"
            )

        r_lead = (args.wind_velocity * 1e5
                  * args.wind_duration_days * 86400.0)
        if args.wind_duration_days > 0 and r_lead < wind_r_outer:
            # Finite-age wind: the leading edge cannot have travelled
            # farther than v_wind * duration.  The gap between the wind
            # front and the ejecta inner edge is filled with an older,
            # tenuous background wind -- a single steady Mdot spanning the
            # whole gap would claim the wind has been flowing for
            # (r_outer - r_inner) / v_wind days.
            if r_lead <= wind_r_inner:
                raise ValueError(
                    f"Wind duration too short: r_lead={r_lead:.3e} cm is "
                    f"inside the shell outer edge {wind_r_inner:.3e} cm.\n"
                    f"--wind-duration-days must exceed "
                    f"{wind_r_inner / (args.wind_velocity * 1e5 * 86400.0):.1f} d "
                    f"at v_wind={args.wind_velocity:.0f} km/s.")
            span1 = np.log10(r_lead) - np.log10(wind_r_inner)
            span2 = np.log10(wind_r_outer) - np.log10(r_lead)
            nz1 = max(20, int(round(100 * span1 / (span1 + span2))))
            nz2 = max(20, 100 - nz1)
            print(f"  Finite-age wind: r_lead={r_lead:.3e} cm "
                  f"({args.wind_duration_days:.0f} d at "
                  f"{args.wind_velocity:.0f} km/s); background beyond")
            young = create_wind(args.wind_mdot, args.wind_velocity,
                                args.wind_T, wind_r_inner, r_lead,
                                nzon_wind=nz1,
                                composition=args.shell_composition)
            # The leading edge is physically a contact discontinuity, but
            # a ~4-dex density step in one zone gives |dlog rho/dlog r|
            # ~800 and stalls the solver.  Ramp log-linearly over ~0.3
            # dex so the gradient stays under the buffer limit of 16.
            lead_span_dex = 0.30
            r_tr_end = r_lead * 10.0 ** lead_span_dex
            if r_tr_end >= wind_r_outer:
                r_tr_end = r_lead * 10.0 ** (0.8 * span2)
            nz_tr = 16
            v_w = args.wind_velocity * 1e5
            rho_bg_at_tr = (args.wind_mdot_bg * MSUN / YR_TO_S
                            / (4 * np.pi * r_tr_end**2 * v_w))
            tr = {k: np.full(nz_tr, v[0]) if isinstance(v, np.ndarray)
                  else v for k, v in young.items()}
            tr['r'] = np.logspace(np.log10(young['r'][-1]),
                                  np.log10(r_tr_end), nz_tr + 1)[1:]
            tr['rho'] = np.logspace(np.log10(young['rho'][-1]),
                                    np.log10(rho_bg_at_tr), nz_tr)
            bg = create_wind(args.wind_mdot_bg, args.wind_velocity,
                             args.wind_T, r_tr_end * 1.001, wind_r_outer,
                             nzon_wind=nz2,
                             composition=args.shell_composition)
            wind = {k: np.concatenate([young[k], tr[k], bg[k]])
                    if isinstance(v, np.ndarray) else v
                    for k, v in young.items()}
            wind['r_inner'] = wind_r_inner
            wind['r_outer'] = wind_r_outer
            tr_lo = np.concatenate([[tr['r'][0] * 0.999], tr['r'][:-1]])
            tr_mass = ((4 / 3) * np.pi
                       * (tr['r']**3 - tr_lo**3) * tr['rho']).sum()
            wind['M_wind'] = young['M_wind'] + bg['M_wind'] + tr_mass
        else:
            wind = create_wind(args.wind_mdot, args.wind_velocity,
                               args.wind_T, wind_r_inner, wind_r_outer,
                               composition=args.shell_composition)

        # Step 4: Combine shell + wind + ejecta
        print(f"\n--- Step 4: Combining shell + wind + ejecta ---")
        combined = combine_shell_wind_ejecta(shell, wind, ejecta, nzon_total=args.nzon,
                                             sw_transition_dex=args.sw_transition_dex,
                                             sw_transition_k=args.sw_transition_k,
                                             we_transition_dex=args.we_transition_dex,
                                             we_transition_k=args.we_transition_k,
                                             ejecta_T_floor=args.ejecta_T_floor,
                                             min_dlogr=args.min_dlogr,
                                             max_dlogr=args.max_dlogr)

    print(f"  Combined: {len(combined['r'])} zones")
    print(f"  r: [{combined['r'].min():.3e}, {combined['r'].max():.3e}] cm")
    print(f"  v: [{combined['v'].min():.1f}, {combined['v'].max():.1f}] cm/s "
          f"([{combined['v'].min()/1e5:.1f}, {combined['v'].max()/1e5:.1f}] km/s)")
    print(f"  rho: [{combined['rho'].min():.3e}, {combined['rho'].max():.3e}] g/cc")
    print(f"  T: [{combined['T'].min():.3e}, {combined['T'].max():.3e}] K")

    # --- Clean outer boundary pile-up artifacts ---
    # STELLA's Lagrangian grid produces velocity deceleration and density
    # non-monotonicity at the outer boundary.  Clean these from the final
    # combined grid before writing STELLA input files.  Skipped in
    # --no-ejecta mode: the outer component is a steady wind (constant v,
    # 1/r^2 rho), and the cleaner's homologous-ejecta velocity fit would
    # corrupt it.
    if not args.no_ejecta:
        combined = clean_combined_outer_boundary(
            combined, ej_start_idx=combined.get('idx_ej'))

    # --- Grid validation (before writing STELLA input files) ---
    # STELLA's implicit solver is sensitive to grid spacing: too-fine zones
    # kill the radiation diffusion timestep, too-coarse zones under-resolve
    # shocks.  Validate the final combined grid against min/max dlogr limits.
    validate_stella_grid(combined['r'], min_dlogr=args.min_dlogr,
                         max_dlogr=args.max_dlogr, strict=args.strict_grid)

    # Step 5: Write STELLA input files
    print(f"\n--- Step 5: Writing STELLA input files ---")
    hyd_path = modmake_dir / f'{model_name}.hyd'
    abn_path = modmake_dir / f'{model_name}.abn'
    xni_path = modmake_dir / f'{model_name}.xni'
    dat_path = modmake_dir / f'{model_name}.dat'
    ip_path = modmake_dir / 'input_profile'

    write_x_hyd(str(hyd_path), combined, t_start_days=args.t_start)
    write_x_abn(str(abn_path), combined)
    write_x_xni(str(xni_path), args.nzon)
    if args.explosion_energy > 0:
        print(f"  WARNING: --explosion-energy={args.explosion_energy} writes a "
              f"thermal bomb (EBurst) into the SN_B .dat.  The SN_B model is "
              f"velocity-driven; EBurst should be 0.")
    write_x_dat(str(dat_path), model_name, args.target_days, args.snapshot_day,
                explosion_energy=args.explosion_energy,
                nstmax=args.nstmax, aq=args.aq, fit_tau=args.fit_tau,
                meth=args.meth, temp_floor=args.temp_floor,
                floor_rad=args.floor_rad, hmin=args.hmin, hmax=args.hmax,
                t_start_days=args.t_start,
                tcura_days=args.tcura_days,
                r_max_cm=combined['r'].max(),
                v_max_cm_s=combined['v'].max(),
                wacc_r=args.wacc_r, wacc_v=args.wacc_v,
                wacc_t=args.wacc_t, wacc_rad=args.wacc_rad,
                tau_tol=args.tau_tol, iout=args.iout,
                disable_scattering=args.disable_scattering)
    write_x_input_profile(str(ip_path), combined)

    # Also copy input_profile to stella dir
    import shutil
    shutil.copy2(str(ip_path), str(stella_dir / 'input_profile'))

    # Step 6: Generate plots
    print(f"\n--- Step 6: Generating plots ---")

    # --- Ejecta smoothing before/after plot ---
    # Generated when smoothing was applied, showing the noise reduction
    # in the raw SN_A ejecta profile.  Uses the pre-truncation smoothed
    # copy so the plot shows the full 500-zone profile (not the truncated
    # version used for STELLA input).
    if ejecta is not None and 'raw_ejecta' in ejecta and 'smoothed_ejecta' in ejecta:
        try:
            smoothing_plot = plots_dir / 'ejecta_smoothing.png'
            plot_ejecta_smoothing(ejecta['raw_ejecta'],
                                  ejecta['smoothed_ejecta'],
                                  str(smoothing_plot),
                                  source_age=args.source_age)
        except Exception as e:
            print(f"  Ejecta smoothing plot failed: {e}")

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        ax_rho_dens, ax_r = axes[0]   # upper-left: Density, upper-right: Velocity
        ax_v, ax_rho = axes[1]        # lower-left: Temperature, lower-right: Composition

        r_s = combined['r']
        v_k = combined['v'] / 1e5

        # Component boundaries (in cm) for vspan highlighting
        r_shell_lo = shell['r_inner']
        r_shell_hi = shell['r_outer']
        r_ej_inner = (ejecta['r'].min() * 0.99) if ejecta is not None else None
        r_max_cm = combined['r'][-1]

        # --- Density profile (upper-left) with shell/wind/ejecta vspans ---
        ax_rho_dens.plot(r_s, combined['rho'], 'k-', linewidth=1.5)
        ax_rho_dens.axvspan(r_shell_lo, r_shell_hi, color='red', alpha=0.15, label='shell')
        if args.no_ejecta:
            ax_rho_dens.axvspan(r_shell_hi, r_max_cm, color='blue', alpha=0.15, label='wind')
        elif args.no_wind:
            ax_rho_dens.axvspan(r_shell_hi, r_ej_inner, color='orange', alpha=0.10, label='buffer')
            ax_rho_dens.axvspan(r_ej_inner, r_max_cm, color='green', alpha=0.15, label='ejecta')
        else:
            ax_rho_dens.axvspan(r_shell_hi, r_ej_inner, color='blue', alpha=0.15, label='wind')
            ax_rho_dens.axvspan(r_ej_inner, r_max_cm, color='green', alpha=0.15, label='ejecta')
        ax_rho_dens.set_xlabel('Radius (cm)')
        ax_rho_dens.set_ylabel('Density (g/cm^3)')
        ax_rho_dens.set_title('Density Profile')
        ax_rho_dens.set_xscale('log')
        ax_rho_dens.set_yscale('log')
        ax_rho_dens.grid(True, alpha=0.3)
        ax_rho_dens.legend(fontsize=9, loc='best')

        ax_r.plot(r_s, v_k, 'b-', linewidth=1.5)
        ax_r.set_xlabel('Radius (cm)')
        ax_r.set_ylabel('Velocity (km/s)')
        ax_r.set_title('Velocity Profile')
        ax_r.set_xscale('log')
        ax_r.grid(True, alpha=0.3)
        # Mark component boundaries
        if args.no_wind:
            ax_r.axvline(shell['r_outer'], color='r', linestyle='--', alpha=0.5, label='shell-ejecta')
        else:
            ax_r.axvline(shell['r_outer'], color='r', linestyle='--', alpha=0.5, label='shell-wind')
            if ejecta is not None:
                ax_r.axvline(ejecta['r'].min() * 0.99, color='g', linestyle='--', alpha=0.5, label='wind-ejecta')
        ax_r.legend()

        ax_v.plot(r_s, combined['T'], 'r-', linewidth=1.5)
        ax_v.set_xlabel('Radius (cm)')
        ax_v.set_ylabel('Temperature (K)')
        ax_v.set_title('Temperature Profile')
        ax_v.set_xscale('log')
        ax_v.set_yscale('log')
        ax_v.grid(True, alpha=0.3)

        # Composition (linear y-axis for mass fractions)
        for key, color, label in [('fH', 'blue', 'H'), ('fHe', 'orange', 'He'),
                                   ('fO', 'green', 'O'), ('fSi', 'red', 'Si'),
                                   ('fFe', 'purple', 'Fe'), ('fCa', 'brown', 'Ca')]:
            if np.any(combined[key] > 1e-10):
                ax_rho.plot(r_s, combined[key], color=color, label=label, linewidth=1.5)
        ax_rho.set_xlabel('Radius (cm)')
        ax_rho.set_ylabel('Mass Fraction')
        ax_rho.set_title('Composition Profile')
        ax_rho.set_xscale('log')
        ax_rho.set_ylim(-0.05, 1.05)
        ax_rho.legend(fontsize=8)
        ax_rho.grid(True, alpha=0.3)

        if args.no_ejecta:
            plt.suptitle('X Model Initial Conditions (shell + CSM wind, no ejecta)', fontsize=14)
        elif args.no_wind:
            plt.suptitle(f'X Model Initial Conditions (shell + {args.source_age}d ejecta, no wind)', fontsize=14)
        else:
            plt.suptitle(f'X Model Initial Conditions (shell + wind + {args.source_age}d ejecta)', fontsize=14)
        plt.tight_layout()
        plt.savefig(str(plots_dir / 'x_model_initial.png'), dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  Saved: {plots_dir}/x_model_initial.png")
    except Exception as e:
        print(f"  Plot generation failed: {e}")

    print(f"\n=== SN_B model preparation complete ===")
    print(f"  STELLA input files in: {modmake_dir}")
    print(f"  Run STELLA with: python scripts/run_model_grid.py {output_dir}/config.yaml")


if __name__ == '__main__':
    main()
