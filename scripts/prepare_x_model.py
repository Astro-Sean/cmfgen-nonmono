#!/usr/bin/env python3
"""
Prepare the X initial model: H/He shell + connecting wind + old SN ejecta.

This script creates the initial conditions for the X model, which represents
a reverse-shock (RSH) H/He shell launched at 5000 km/s into a dense wind,
which connects to a 400-day-old metal-rich supernova ejecta (the "CSM").

The X model consists of three components (inner to outer):
  1. Inner zones: A new H/He shell (1.0 Msun, X=0.7, Y=0.3) at 5000 km/s
  2. Middle zones: A 1/r^2 wind (Mdot=1e-3 Msun/yr, v_wind=600 km/s)
     connecting the shell to the old ejecta
  3. Outer zones: The 400-day homologously evolved SN_12C_CSM2 ejecta
     (metal-rich, the "CSM" for the interaction)

STELLA starts at t=0 (the shell launch time). The old SN ejecta is already
400 days old at this point — its structure represents the pre-existing CSM.

Outputs:
  - {model}.hyd: STELLA hydrodynamics input
  - {model}.abn: STELLA composition input
  - {model}.xni: Ni56 distribution (zero — no radioactive decay in X)
  - {model}.dat: STELLA control file
  - input_profile: 12-column profile for convert_profile_to_stella.py

Usage:
    python prepare_x_model.py <source_dir> <output_dir> [options]

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


def read_source_ejecta(source_dir, source_age=400, tolerance=100, t_start_days=0.0,
                       force_homologous_velocity=False):
    """Read the source SN ejecta profile at the correct age for STELLA start.

    Reads the STELLA .swd file from the source model, finds the latest epoch,
    and homologously evolves it to (source_age + t_start_days) days.

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

    return epoch_data


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
    if n_rho > 3:
        factor = 1 + (3.0 / (n_rho - 3)) * (1 - frac_v0**(n_rho - 3))
    else:
        factor = 1 + (3.0 / (3 - n_rho)) * (frac_v0**(n_rho - 3) - 1)
    rho_0 = M_shell / ((4/3) * np.pi * r0**3 * factor)

    # Density at outer edge
    rho_edge = rho_0 * frac_v0**n_rho

    # Diffusion time diagnostic
    kappa_est = 0.34
    c_light = 3e10
    t_diff = kappa_est * rho_0 * r_outer**2 / c_light

    print(f"  H/He shell (Dessart 2010jl structure, t_start={t_start_days}d):")
    print(f"    M={shell_mass_msun} Msun, v_outer={shell_velocity_kms} km/s")
    print(f"    N_rho={n_rho}, V0={v0_kms} km/s (core/envelope transition)")
    print(f"    r_inner={r_inner:.3e} cm, r0={r0:.3e} cm, r_outer={r_outer:.3e} cm")
    print(f"    rho_0={rho_0:.3e} g/cc (core), rho_edge={rho_edge:.3e} g/cc (outer)")
    print(f"    Density contrast (core/edge) = {rho_0/rho_edge:.1f}")
    print(f"    Diffusion time: {t_diff/86400:.1f} days (must be < ~10x target)")
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

    A sigmoid transition is applied at both boundaries for smooth density
    gradients (avoiding solver stalls per AGENTS.md lessons).
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


def sigmoid_transition(r, r_transition, width_frac, k, rho_inner, rho_outer):
    """Apply a sigmoid density transition between two profiles.

    Blends rho_inner (dominant at r << r_transition) with rho_outer (dominant
    at r >> r_transition) using a logistic sigmoid.

    The transition is centered at r_transition with a width determined by
    width_frac * r_transition and steepness k.
    """
    x = (np.log10(r) - np.log10(r_transition)) / (width_frac * np.log10(r_transition))
    weight = 1.0 / (1.0 + np.exp(k * x))  # 1 at small r, 0 at large r
    return rho_inner * weight + rho_outer * (1 - weight)


def _sigmoid_blend_weights(n, transition_zones=4):
    """Sigmoid blend weights for n buffer zones.

    Returns array of length n with values:
      0.0 for inner zones (shell side)
      sigmoid transition over `transition_zones` at the midpoint
      1.0 for outer zones (ejecta side)

    This produces a sharp interface over ~transition_zones zones instead
    of a gradual linear blend across the entire buffer.
    """
    if n <= transition_zones * 2:
        # Buffer too small for sharp transition — use linear
        return np.linspace(0.0, 1.0, n)
    center = n / 2.0
    # Steepness: transition_zones/2 zones on each side of center
    # sigmoid goes from ~0 to ~1 over ~4*steepness indices
    steepness = 10.0 / transition_zones
    idx = np.arange(n) - center
    w = 1.0 / (1.0 + np.exp(-steepness * idx))
    return np.clip(w, 0.0, 1.0)


def _sigmoid_blend(n, val_inner, val_outer, transition_zones=4):
    """Blend two scalar values across n zones using a sharp sigmoid.

    Inner zones get val_inner, outer zones get val_outer,
    with a sigmoid transition over ~transition_zones at the midpoint.
    """
    w = _sigmoid_blend_weights(n, transition_zones)
    return val_inner * (1.0 - w) + val_outer * w


def combine_shell_ejecta(shell, ejecta, nzon_total=500,
                          se_transition_dex=0.5, immediate_contact=False,
                          ejecta_T_floor=1000.0):
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
    from scipy.interpolate import interp1d

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
        se_contrast = abs(np.log10(rho_se_end / max(rho_se_start, 1e-300)))
        print(f"    rho: {rho_se_start:.3e} -> {rho_se_end:.3e} g/cc "
              f"(contrast={se_contrast:.1f} dex, smooth over {nzon_se_buffer} zones)")
    else:
        # --- Buffer spans the full shell->ejecta gap ---
        rho_se_start = shell['rho'][-1]  # shell density at outer edge
        rho_se_end = ej_rho_sorted[0]    # ejecta density at inner edge
        se_contrast = abs(np.log10(rho_se_end / max(rho_se_start, 1e-300)))
        se_span_dex = log_r_ej - log_r_shell  # full gap width in dex

        # Number of buffer zones: enough for density gradient < 16 and
        # composition resolution (~0.05 dex per zone minimum)
        nzon_se_buffer = max(30, int(se_span_dex / 0.05))
        nzon_se_buffer = min(nzon_se_buffer, 150)  # cap at 150

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
        nzon_shell = max(50, int(remaining * 0.25))
        nzon_ejecta = remaining - nzon_shell
    else:
        # Buffer gets its zones first, then shell/ejecta split the remainder
        remaining = nzon_total - nzon_se_buffer
        nzon_shell = max(50, int(remaining * 0.25))
        nzon_ejecta = remaining - nzon_shell
    n_total = nzon_shell + nzon_se_buffer + nzon_ejecta
    print(f"  Zone allocation: {nzon_shell} shell + {nzon_se_buffer} SE buffer + "
          f"{nzon_ejecta} ejecta = {n_total} total")

    # --- Build grids ---
    # Shell: log-spaced, interpolated from the Dessart power-law profile
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
        n_rho_val = 8.0  # Dessart power-law exponent (hardcoded in create_hhe_shell)
        frac_v0 = shell['r0'] / shell['r_outer']
        factor = 1 + (3.0 / (n_rho_val - 3)) * (1 - frac_v0**(n_rho_val - 3))
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
        se_log_rho = np.linspace(np.log10(rho_se_start), np.log10(max(rho_se_end, 1e-300)),
                                 nzon_se_buffer)
        se_buf_rho = 10.0 ** se_log_rho
        # Velocity: sharp sigmoid transition from shell v to ejecta v.
        # The transition happens over ~4 zones at the buffer midpoint,
        # so most of the buffer has shell velocity (inner) or ejecta
        # velocity (outer) with a sharp interface.
        v_shell_edge = shell['v'][-1]
        v_ej_inner = ej_v_sorted[0]
        se_buf_v = _sigmoid_blend(nzon_se_buffer, v_shell_edge, v_ej_inner,
                                  transition_zones=4)
        # Temperature: blend from shell T to ejecta T across the transition.
        # For immediate contact (thin 5-zone transition), this provides a
        # smooth temperature gradient at the contact discontinuity.
        # For the wide buffer case, the buffer is cold (ejecta T floor).
        if immediate_contact:
            se_buf_T = _sigmoid_blend(nzon_se_buffer, shell['T'][-1],
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

    def log_interp(y, r_target):
        log_y = np.log10(np.maximum(y, 1e-300))
        ip = interp1d(log_r_ej_src, log_y, kind='linear',
                      bounds_error=False, fill_value=(log_y[0], log_y[-1]))
        return 10.0 ** ip(np.log10(r_target))

    def lin_interp(y, r_target):
        ip = interp1d(log_r_ej_src, y, kind='linear',
                      bounds_error=False, fill_value=(y[0], y[-1]))
        return ip(np.log10(r_target))

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

    # Blend weights for buffer: sharp sigmoid centered at buffer midpoint.
    # Transitions over ~4 zones so the composition interface is sharp
    # (shell H/He -> ejecta metals over 3-5 zones, not the full buffer).
    se_w = _sigmoid_blend_weights(nzon_se_buffer, transition_zones=4)

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
        **comp_combined,
    }


def combine_shell_wind_ejecta(shell, wind, ejecta, nzon_total=500,
                               sw_transition_dex=0.5, sw_transition_k=3.0,
                               we_transition_dex=0.5, we_transition_k=3.0,
                               ejecta_T_floor=1000.0):
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
    """
    from scipy.interpolate import interp1d

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
    mdot = wind['rho'][0] * 4 * np.pi * wind['r_inner']**2 * wind['v'][0]
    v_wind = wind['v'][0]
    T_wind = wind['T'][0]

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
    rho_sw_end = mdot / (4 * np.pi * r_sw_buffer_end**2 * v_wind)  # wind density at buffer end
    sw_contrast = np.log10(rho_sw_start / max(rho_sw_end, 1e-300))

    # Number of buffer zones: enough for gradient < 16
    # gradient = contrast / (n_buffer * dlog_r_per_zone)
    # We want gradient < 16, so n_buffer > contrast / (16 * dlog_r_per_zone)
    # Use dlog_r_per_zone ≈ sw_transition_dex / n_buffer (self-consistent)
    # gradient = contrast / sw_transition_dex (independent of n_buffer!)
    # So we need sw_transition_dex > contrast / 16
    sw_grad = sw_contrast / sw_transition_dex
    nzon_sw_buffer = max(20, int(sw_contrast / 16 * sw_transition_dex / 0.02))
    nzon_sw_buffer = min(nzon_sw_buffer, 60)  # cap at 60

    print(f"  Shell-wind buffer: {nzon_sw_buffer} zones, "
          f"r=[{r_shell_outer:.3e}, {r_sw_buffer_end:.3e}] cm")
    print(f"    rho: {rho_sw_start:.3e} -> {rho_sw_end:.3e} g/cc "
          f"(contrast={sw_contrast:.1f} dex, est. gradient={sw_grad:.1f})")

    # Wind-ejecta buffer: ramps from rho_wind to rho_ejecta
    r_we_buffer_start = r_ej_inner / (10 ** we_transition_dex)
    rho_we_start = mdot / (4 * np.pi * r_we_buffer_start**2 * v_wind)
    rho_we_end = ej_rho_sorted[0]  # ejecta density at inner edge
    we_contrast = abs(np.log10(rho_we_end / max(rho_we_start, 1e-300)))
    we_grad = we_contrast / we_transition_dex
    nzon_we_buffer = max(10, int(we_contrast / 16 * we_transition_dex / 0.02))
    nzon_we_buffer = min(nzon_we_buffer, 40)

    print(f"  Wind-ejecta buffer: {nzon_we_buffer} zones, "
          f"r=[{r_we_buffer_start:.3e}, {r_ej_inner:.3e}] cm")
    print(f"    rho: {rho_we_start:.3e} -> {rho_we_end:.3e} g/cc "
          f"(contrast={we_contrast:.1f} dex, est. gradient={we_grad:.1f})")

    # --- Zone allocation ---
    # Shell: 15% of total, wind: 25%, ejecta: 60%, buffers taken from wind
    nzon_shell = max(50, int(nzon_total * 0.15))
    nzon_ejecta = max(100, int(nzon_total * 0.50))
    nzon_wind = nzon_total - nzon_shell - nzon_ejecta - nzon_sw_buffer - nzon_we_buffer
    nzon_wind = max(50, nzon_wind)  # ensure minimum wind zones

    # Recalculate total (may exceed nzon_total slightly if minimums kick in)
    n_total = nzon_shell + nzon_sw_buffer + nzon_wind + nzon_we_buffer + nzon_ejecta
    print(f"  Zone allocation: {nzon_shell} shell + {nzon_sw_buffer} SW buffer + "
          f"{nzon_wind} wind + {nzon_we_buffer} WE buffer + {nzon_ejecta} ejecta "
          f"= {n_total} total")

    # --- Build grids ---
    # Shell: log-spaced, interpolated from the Dessart power-law profile
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
        n_rho_val = 8.0  # Dessart power-law exponent (hardcoded in create_hhe_shell)
        frac_v0 = shell['r0'] / shell['r_outer']
        factor = 1 + (3.0 / (n_rho_val - 3)) * (1 - frac_v0**(n_rho_val - 3))
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
    sw_log_rho = np.linspace(np.log10(rho_sw_start), np.log10(max(rho_sw_end, 1e-300)),
                             nzon_sw_buffer)
    sw_buf_rho = 10.0 ** sw_log_rho
    # Velocity: linear from shell v at outer edge to wind v
    sw_buf_v = np.linspace(shell['v'][-1], v_wind, nzon_sw_buffer)
    # Temperature: linear from shell T at outer edge to wind T
    sw_buf_T = np.linspace(shell['T'][-1], T_wind, nzon_sw_buffer)

    # Wind: log-spaced, 1/r^2 density, constant velocity
    wind_r = np.logspace(np.log10(r_sw_buffer_end), np.log10(r_we_buffer_start), nzon_wind + 1)[1:]
    wind_rho = mdot / (4 * np.pi * wind_r**2 * v_wind)
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

    def log_interp(y, r_target):
        log_y = np.log10(np.maximum(y, 1e-300))
        ip = interp1d(log_r_ej_src, log_y, kind='linear',
                      bounds_error=False, fill_value=(log_y[0], log_y[-1]))
        return 10.0 ** ip(np.log10(r_target))

    def lin_interp(y, r_target):
        ip = interp1d(log_r_ej_src, y, kind='linear',
                      bounds_error=False, fill_value=(y[0], y[-1]))
        return ip(np.log10(r_target))

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

        comp_combined[key] = np.concatenate([
            np.full(nzon_shell, shell_val),
            np.full(nzon_sw_buffer, shell_val),  # buffer = shell composition
            np.full(nzon_wind, wind_val),
            np.full(nzon_we_buffer, wind_val),   # buffer = wind composition
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
        **comp_combined,
    }


def _steep_sigmoid_weights(n, k=5.0):
    """Steep sigmoid blend weights for n zones.

    Returns array of length n with values from 0 (inner side) to 1 (outer
    side), transitioning sharply at the midpoint.

    The steepness parameter k controls how sharp the transition is:
      k=1  : gradual (transition spans ~n/2 zones on each side)
      k=5  : steep   (transition spans ~n/5 zones on each side)
      k=10 : very steep (near-step, ~n/10 zones on each side)

    This is the CD (contact discontinuity) sigmoid used by the
    density-matched wind-ejecta combiner.  It is intentionally steeper
    than the shell-ejecta _sigmoid_blend (which uses transition_zones=4)
    because the CD is a genuine physical discontinuity: density is
    matched by construction, but velocity, temperature and composition
    jump sharply.
    """
    if n <= 1:
        return np.array([1.0])
    center = (n - 1) / 2.0
    idx = np.arange(n) - center
    # Normalize so k has a consistent meaning independent of n:
    # the transition width is ~ n / (2*k) zones on each side of center.
    w = 1.0 / (1.0 + np.exp(-k * idx / max(n / 2.0, 1.0)))
    return np.clip(w, 0.0, 1.0)


def _steep_sigmoid_blend(n, val_inner, val_outer, k=5.0):
    """Blend two scalar values across n zones using a steep sigmoid."""
    w = _steep_sigmoid_weights(n, k=k)
    return val_inner * (1.0 - w) + val_outer * w


def combine_shell_wind_ejecta_density_matched(
    shell, ejecta, nzon_total=500,
    wind_velocity_kms=100.0, wind_T_K=2000.0,
    composition='solar', ejecta_T_floor=2000.0,
    sw_transition_dex=0.5, cd_transition_zones=5, cd_sigmoid_k=5.0,
    v_blend_zones=30,
):
    """Density-matched wind-ejecta stitching with a steep sigmoid CD.

    Derives Mdot so that rho_wind(r_interface) = rho_ejecta(r_interface),
    i.e. the wind and ejecta are in pressure/density balance at the contact
    radius.  This eliminates the density cliff at the wind-ejecta interface
    that motivates the extended log-linear buffer in the current combiner
    (combine_shell_wind_ejecta).

    Layout (inner -> outer):
      shell | shell-wind buffer | wind | CD sigmoid | ejecta

    - The **shell-wind buffer** is preserved (log-linear density ramp),
      because the shell edge is much denser than the wind at that radius --
      a genuine physical discontinuity that needs a buffer to keep
      dlog(rho)/dlog(r) < 16.

    - The **wind-ejecta interface** uses a thin steep sigmoid
      (cd_transition_zones zones, steepness cd_sigmoid_k) for velocity,
      temperature and composition.  The density is continuous by
      construction (matched at r_interface via the derived Mdot), so no
      density buffer is needed there.

    The wind follows rho = Mdot / (4 pi r^2 v_wind) with the derived Mdot.
    The input wind_mdot is IGNORED -- Mdot is computed from the ejecta
    inner-edge density.

    Parameters:
        wind_velocity_kms:  Wind velocity (constant).  Mdot is derived
                            consistent with this v_wind.
        wind_T_K:           Wind temperature (isothermal).
        composition:        Wind composition ('solar' or 'hhe').
        sw_transition_dex:  Shell-wind buffer width in dex of log10(r).
        cd_transition_zones: Number of zones for the CD steep sigmoid.
        cd_sigmoid_k:       Steepness of the CD sigmoid (see
                            _steep_sigmoid_weights).
    """
    from scipy.interpolate import interp1d

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
    print(f"    CD sigmoid: {cd_transition_zones} zones, k={cd_sigmoid_k}")

    # --- Shell-wind buffer (same logic as combine_shell_wind_ejecta) ---
    total_gap_dex = log_r_ej - log_r_shell
    max_buffer_dex = total_gap_dex / 3.0
    sw_transition_dex = min(sw_transition_dex, max_buffer_dex)

    r_sw_buffer_end = r_shell_outer * (10 ** sw_transition_dex)
    rho_sw_start = shell['rho'][-1]
    rho_sw_end = mdot_matched / (4.0 * np.pi * r_sw_buffer_end**2 * v_wind)
    sw_contrast = np.log10(rho_sw_start / max(rho_sw_end, 1e-300))
    sw_grad = sw_contrast / sw_transition_dex
    nzon_sw_buffer = max(20, int(sw_contrast / 16 * sw_transition_dex / 0.02))
    nzon_sw_buffer = min(nzon_sw_buffer, 60)

    print(f"    Shell-wind buffer: {nzon_sw_buffer} zones, "
          f"r=[{r_shell_outer:.3e}, {r_sw_buffer_end:.3e}] cm")
    print(f"      rho: {rho_sw_start:.3e} -> {rho_sw_end:.3e} g/cc "
          f"(contrast={sw_contrast:.1f} dex, est. gradient={sw_grad:.1f})")

    # --- Zone allocation ---
    nzon_shell = max(50, int(nzon_total * 0.15))
    nzon_ejecta = max(100, int(nzon_total * 0.55))
    nzon_wind = (nzon_total - nzon_shell - nzon_ejecta
                 - nzon_sw_buffer - cd_transition_zones)
    nzon_wind = max(40, nzon_wind)
    n_total = (nzon_shell + nzon_sw_buffer + nzon_wind
               + cd_transition_zones + nzon_ejecta)
    print(f"    Zone allocation: {nzon_shell} shell + {nzon_sw_buffer} SW buffer "
          f"+ {nzon_wind} wind + {cd_transition_zones} CD-sigmoid + "
          f"{nzon_ejecta} ejecta = {n_total}")

    # --- Shell grid (log-spaced, interpolated from Dessart profile) ---
    shell_r = np.logspace(log_r_min, log_r_shell, nzon_shell)
    log_r_src = np.log10(shell['r'])
    log_rho_src = np.log10(np.maximum(shell['rho'], 1e-300))
    shell_rho = 10.0 ** np.interp(np.log10(shell_r), log_r_src, log_rho_src)
    shell_v = np.interp(np.log10(shell_r), log_r_src, shell['v'])
    shell_T = np.interp(np.log10(shell_r), log_r_src, shell['T'])

    # Shell mass conservation (same as combine_shell_wind_ejecta)
    if 'rho_0' in shell and 'r0' in shell:
        n_rho_val = 8.0
        frac_v0 = shell['r0'] / shell['r_outer']
        factor = 1 + (3.0 / (n_rho_val - 3)) * (1 - frac_v0**(n_rho_val - 3))
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
    sw_log_rho = np.linspace(np.log10(rho_sw_start),
                             np.log10(max(rho_sw_end, 1e-300)),
                             nzon_sw_buffer)
    sw_buf_rho = 10.0 ** sw_log_rho
    sw_buf_v = np.linspace(shell['v'][-1], v_wind, nzon_sw_buffer)
    sw_buf_T = np.linspace(shell['T'][-1], T_wind, nzon_sw_buffer)

    # --- Wind grid (1/r^2 with derived Mdot) ---
    # Wind spans from r_sw_buffer_end to r_ej_inner, with the last
    # cd_transition_zones points being the CD sigmoid.
    wind_r = np.logspace(np.log10(r_sw_buffer_end), log_r_ej,
                         nzon_wind + cd_transition_zones + 1)[1:]
    wind_r_pure = wind_r[:nzon_wind]
    cd_r = wind_r[nzon_wind:]

    wind_rho_pure = mdot_matched / (4.0 * np.pi * wind_r_pure**2 * v_wind)
    wind_v_pure = np.full(nzon_wind, v_wind)
    wind_T_pure = np.full(nzon_wind, T_wind)

    # --- Contact-discontinuity steep sigmoid ---
    # Density: matched by construction -- wind 1/r^2 law continues through
    # the CD points, and at r_ej_inner rho_wind = rho_ejecta by Mdot choice.
    cd_rho = mdot_matched / (4.0 * np.pi * cd_r**2 * v_wind)

    # Velocity: FLAT at v_wind through the CD.  The velocity blend from
    # the wind's constant v_wind to the ejecta's homologous v(r) profile
    # is done over a transition region at the start of the ejecta (see
    # below), not as a sharp jump at the CD.
    cd_v = np.full(cd_transition_zones, v_wind)

    # Temperature: steep sigmoid from T_wind to ejecta T (floored)
    cd_T = _steep_sigmoid_blend(cd_transition_zones, T_wind,
                                max(ej_T_sorted[0], ejecta_T_floor),
                                k=cd_sigmoid_k)

    # --- Ejecta grid (log-spaced, interpolated from source) ---
    ej_r_new = np.logspace(log_r_ej, log_r_max, nzon_ejecta + 1)[1:]

    def log_interp(y, r_target):
        log_y = np.log10(np.maximum(y, 1e-300))
        ip = interp1d(log_r_ej_src, log_y, kind='linear',
                      bounds_error=False,
                      fill_value=(log_y[0], log_y[-1]))
        return 10.0 ** ip(np.log10(r_target))

    def lin_interp(y, r_target):
        ip = interp1d(log_r_ej_src, y, kind='linear',
                      bounds_error=False,
                      fill_value=(y[0], y[-1]))
        return ip(np.log10(r_target))

    ej_rho = log_interp(ej_rho_sorted, ej_r_new)
    ej_v = lin_interp(ej_v_sorted, ej_r_new)
    ej_T = np.maximum(log_interp(ej_T_sorted, ej_r_new), ejecta_T_floor)

    # --- Homologous velocity blend at the ejecta inner edge ---
    # The wind has constant v_wind.  The ejecta has a homologous v(r)
    # profile (v ~ r/t for free expansion).  Rather than a sharp jump at
    # the CD, blend smoothly from v_wind to the ejecta's natural v(r)
    # over the first v_blend_zones of the ejecta.
    #
    # The blend weight w goes from 0 (at the inner edge, v = v_wind) to
    # 1 (at the outer edge of the blend region, v = v_ej_natural).  We
    # use a smooth sigmoid in log(r) for the blend.
    v_blend_zones = min(v_blend_zones, nzon_ejecta // 3)
    if v_blend_zones > 1:
        blend_w = _steep_sigmoid_weights(v_blend_zones, k=3.0)
        # Pad to full ejecta length (weight=1 beyond the blend region)
        full_w = np.concatenate([blend_w, np.ones(nzon_ejecta - v_blend_zones)])
        ej_v = v_wind * (1.0 - full_w) + ej_v * full_w
        print(f"    Velocity blend: {v_blend_zones} zones, "
              f"v_wind={v_wind/1e5:.1f} km/s -> v_ej(r) homologous")

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

    # CD sigmoid composition: steep blend wind -> ejecta inner-edge composition
    ej_comp_inner = {}
    for key in all_elem_keys:
        if key in ejecta:
            ej_comp_inner[key] = float(ejecta[key][ej_order][0])
        else:
            ej_comp_inner[key] = 0.0

    cd_w = _steep_sigmoid_weights(cd_transition_zones, k=cd_sigmoid_k)

    comp_combined = {}
    for key in all_elem_keys:
        sv = shell_comp_vals[key]
        wv = wind_comp_vals[key]
        ev = ej_comp_inner[key]
        cd_comp = wv * (1.0 - cd_w) + ev * cd_w
        comp_combined[key] = np.concatenate([
            np.full(nzon_shell, sv),
            np.full(nzon_sw_buffer, sv),   # SW buffer = shell composition
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
            print(f"    CD sigmoid gradients: max={np.max(cd_grads):.2f}, "
                  f"mean={np.mean(cd_grads):.2f}")

    return {
        'r': r_combined,
        'rho': rho_combined,
        'v': v_combined,
        'T': T_combined,
        **comp_combined,
    }


def write_x_hyd(fname, g, mass_cut_msun=0.0, t_start_days=0.0):
    """Write STELLA .hyd file for X model.

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
        f.write(f"{t_start_days:12.3E}{nzon:6d}{mass_cut_msun:13.5E}{rcen:13.5E}{rho_cen:13.5E}\n")
        for km in range(nzon):
            dMr = dm[km] / MSUN
            Mr = M_r[km] / MSUN + mass_cut_msun
            # km  dMr  rHyd  rhoHyd  TpHyd  uHyd  aMr  dum
            f.write(f" {km+1:4d}{dMr:12.4E}{r[km]:15.7E}{rho[km]:12.4E}{T[km]:12.4E}"
                    f"{v[km]:12.4E}{Mr:12.4E}{0.0:12.4E}\n")
    print(f"Written {fname}  ({nzon} zones, M_total={M_r[-1]/MSUN:.4f} Msun, "
          f"t_start={t_start_days}d)")


def write_x_abn(fname, g):
    """Write STELLA .abn file for X model.

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
    """Write a zero .xni file (no radioactive decay in X model)."""
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
                t_start_days=0.0, r_max_cm=0.0):
    """Write STELLA .dat control file for X model.

    X model uses EKO=0 (preserve velocities) since the shell already has
    its launch velocity encoded in the .hyd file. Uses high AQ for
    non-monotonic velocity (shell 5000 km/s -> wind 600 km/s -> ejecta).

    t_start_days: simulation starts at this PROPER_T (shell pre-expanded).
    target_days: desired OBS.TIME (observer-frame days since explosion) to reach.
    snapshot_day: desired OBS.TIME for the CMFGEN/SuperLite snapshot.
    r_max_cm: outermost radius of the .hyd profile (cm). Used to compute the
              light-travel-time offset (r_max/c) that STELLA subtracts from
              PROPER_T to get OBS.TIME. TcurB (PROPER_T stop) is set to
              target_days + r_max/c so that OBS.TIME reaches target_days.
    """
    c_light = 2.998e10  # cm/s
    light_time_days = r_max_cm / c_light / 86400.0 if r_max_cm > 0 else 0.0

    t_start = t_start_days  # PROPER_T at simulation start
    # TcurB: PROPER_T at which OBS.TIME = target_days
    # OBS.TIME = PROPER_T - r_max/c, so PROPER_T = target_days + r_max/c
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
        # EKO=0: preserve .hyd velocities
        f.write(f"   0.0d0                      0.d0  9.d-1  1.d0          +1.\n")
        f.write(f"   THRMAT   CRAP     CONV     EDTM     CHNCND  Givdtl\n")
        f.write(f"   1.D-30   1.D0     F        T        T       T\n")
        f.write(f"   FLOOR(R) FLOOR(V) FLOOR(T) FLOOR(RADIAT)    FLOOR(UCONV)\n")
        f.write(f"   1.D-12   1.D-03   {temp_floor:.4g}D+00   {floor_rad:.4g}           1.D-04\n")
        f.write(f"   Wacc(R)  Wacc(V)  Wacc(T)  Wacc(RADIAT)\n")
        f.write(f"   1d-0    1d-0    1.0D-00   1D-00\n")
        f.write(f"   FitTau   TAUTOL   Rvis   AQ    BQ       DRT     NRT     SCAT\n")
        f.write(f"   {fit_tau:.1f}D0   1.1D0     1.D0   {aq:.1f}D0  {bq:.1f}D+00   1.      1       T\n")
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

    header = ("# input_profile generated by prepare_x_model.py\n"
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
    parser.add_argument('source_dir', help='Source SN model directory (e.g. SN_12C_CSM2)')
    parser.add_argument('output_dir', help='X model output directory (e.g. X)')
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
                             'or solar (Asplund+2009, X=0.738 Y=0.249 Z=0.013)')
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
                        help='EBurst in 1e50 erg (0=no bomb)')
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
    parser.add_argument('--sw-transition-dex', type=float, default=0.5,
                        help='Shell-wind sigmoid width in dex of log10(r) (default: 0.5)')
    parser.add_argument('--sw-transition-k', type=float, default=3.0,
                        help='Shell-wind sigmoid steepness (default: 3.0)')
    parser.add_argument('--we-transition-dex', type=float, default=0.5,
                        help='Wind-ejecta sigmoid width in dex of log10(r) (default: 0.5)')
    parser.add_argument('--we-transition-k', type=float, default=3.0,
                        help='Wind-ejecta sigmoid steepness (default: 3.0)')
    parser.add_argument('--no-wind', action='store_true',
                        help='Skip the wind component: shell + ejecta only (no wind in between)')
    parser.add_argument('--immediate-contact', action='store_true',
                        help='Shell and ejecta in direct contact (no buffer zone). '
                             'Computes self-consistent t_start so shell outer edge '
                             'meets ejecta inner edge. Requires --no-wind.')
    parser.add_argument('--se-transition-dex', type=float, default=0.5,
                        help='Shell-ejecta buffer width in dex of log10(r) (no-wind mode, default: 0.5)')
    parser.add_argument('--density-matched-wind', action='store_true',
                        help='Use density-matched wind-ejecta stitching: derive Mdot so '
                             'rho_wind(r_interface) = rho_ejecta(r_interface), with a '
                             'thin steep sigmoid for v/T/composition at the contact '
                             'discontinuity (instead of a log-linear density buffer). '
                             'The shell-wind buffer is preserved. The input --wind-mdot '
                             'is IGNORED in this mode (Mdot is derived).')
    parser.add_argument('--cd-transition-zones', type=int, default=5,
                        help='Number of zones for the CD steep sigmoid '
                             '(density-matched-wind mode, default: 5)')
    parser.add_argument('--cd-sigmoid-k', type=float, default=5.0,
                        help='Steepness of the CD sigmoid (density-matched-wind mode, '
                             'default: 5.0). Higher = sharper transition. k=1 gradual, '
                             'k=5 steep, k=10 near-step.')
    parser.add_argument('--v-blend-zones', type=int, default=30,
                        help='Number of ejecta zones over which to homologously blend '
                             'velocity from v_wind to the ejecta natural v(r) profile '
                             '(density-matched-wind mode, default: 30). The wind and CD '
                             'have flat v_wind; the blend happens at the start of the '
                             'ejecta region.')
    parser.add_argument('--force-homologous-velocity', action='store_true',
                        help='Force the SN_A ejecta velocity to pure homologous expansion '
                             '(v = r/t). Removes residual non-homologous structure from '
                             'the source STELLA model (e.g. CSM deceleration, shocks). '
                             'Physically appropriate for an SN_A ejecta that is hundreds '
                             'of days old and should be in free expansion.')
    args = parser.parse_args()

    source_dir = Path(args.source_dir)
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
    if args.no_wind:
        print(f"  Wind: SKIPPED (--no-wind mode)")
    elif args.density_matched_wind:
        print(f"  Wind: DENSITY-MATCHED (Mdot derived, v={args.wind_velocity} km/s, "
              f"CD sigmoid k={args.cd_sigmoid_k})")
    else:
        print(f"  Wind: Mdot={args.wind_mdot} Msun/yr, v={args.wind_velocity} km/s")
    print(f"  Source age: {args.source_age} days")
    print(f"  STELLA target: {args.target_days} days, CMFGEN snapshot: {args.snapshot_day} days")

    # Step 1: Read source ejecta and homologously evolve to source_age + t_start_days
    # The SN_A ejecta is source_age days old at SN_B explosion (t=0).
    # STELLA starts at t_start_days later, so the ejecta is (source_age + t_start_days) days old.
    effective_age = args.source_age + args.t_start
    print(f"\n--- Step 1: Reading source ejecta at {effective_age} days "
          f"({args.source_age}d at SN_B + {args.t_start}d t_start) ---")
    ejecta = read_source_ejecta(source_dir, source_age=args.source_age,
                                t_start_days=args.t_start,
                                force_homologous_velocity=args.force_homologous_velocity)

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
    if args.ejecta_rho_min > 0:
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
    if args.ejecta_r_max_cm > 0:
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

    # Step 2: Create H/He shell (Dessart 2010jl power-law structure)
    # For no-wind mode: place shell outer edge directly at ejecta inner edge
    # (immediate contact — the shell is the bullet, the ejecta is the target).
    # This eliminates the large gap that would otherwise exist because the
    # shell (launched at v_shell from center) hasn't physically reached the
    # ejecta inner edge at t_start days. Instead, we start the simulation at
    # the moment of contact, like Dessart+2015.
    r_ej_inner = ejecta['r'].min() * 0.99  # ejecta inner edge (with small overlap)
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

    # Step 3: Create connecting wind (or skip if --no-wind)
    if args.no_wind:
        print(f"\n--- Step 3: SKIPPING wind (--no-wind mode, immediate contact) ---")
        combined = combine_shell_ejecta(shell, ejecta, nzon_total=args.nzon,
                                        se_transition_dex=args.se_transition_dex,
                                        immediate_contact=True,
                                        ejecta_T_floor=args.ejecta_T_floor)
    elif args.density_matched_wind:
        print(f"\n--- Step 3: DENSITY-MATCHED wind-ejecta stitching ---")
        print(f"    (Mdot derived from rho_ejecta(r_interface), "
              f"CD sigmoid: {args.cd_transition_zones} zones, k={args.cd_sigmoid_k})")
        combined = combine_shell_wind_ejecta_density_matched(
            shell, ejecta, nzon_total=args.nzon,
            wind_velocity_kms=args.wind_velocity,
            wind_T_K=args.wind_T,
            composition=args.shell_composition,
            ejecta_T_floor=args.ejecta_T_floor,
            sw_transition_dex=args.sw_transition_dex,
            cd_transition_zones=args.cd_transition_zones,
            cd_sigmoid_k=args.cd_sigmoid_k,
            v_blend_zones=args.v_blend_zones)
    else:
        print(f"\n--- Step 3: Creating connecting wind ---")
        # Wind spans from shell outer edge to ejecta inner edge
        wind_r_inner = shell['r_outer'] * 1.01
        wind_r_outer = ejecta['r'].min() * 0.99
        wind = create_wind(args.wind_mdot, args.wind_velocity, args.wind_T,
                           wind_r_inner, wind_r_outer,
                           composition=args.shell_composition)

        # Step 4: Combine shell + wind + ejecta
        print(f"\n--- Step 4: Combining shell + wind + ejecta ---")
        combined = combine_shell_wind_ejecta(shell, wind, ejecta, nzon_total=args.nzon,
                                             sw_transition_dex=args.sw_transition_dex,
                                             sw_transition_k=args.sw_transition_k,
                                             we_transition_dex=args.we_transition_dex,
                                             we_transition_k=args.we_transition_k,
                                             ejecta_T_floor=args.ejecta_T_floor)

    print(f"  Combined: {len(combined['r'])} zones")
    print(f"  r: [{combined['r'].min():.3e}, {combined['r'].max():.3e}] cm")
    print(f"  v: [{combined['v'].min():.1f}, {combined['v'].max():.1f}] cm/s "
          f"([{combined['v'].min()/1e5:.1f}, {combined['v'].max()/1e5:.1f}] km/s)")
    print(f"  rho: [{combined['rho'].min():.3e}, {combined['rho'].max():.3e}] g/cc")
    print(f"  T: [{combined['T'].min():.3e}, {combined['T'].max():.3e}] K")

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
    write_x_dat(str(dat_path), model_name, args.target_days, args.snapshot_day,
                explosion_energy=args.explosion_energy,
                nstmax=args.nstmax, aq=args.aq, fit_tau=args.fit_tau,
                meth=args.meth, temp_floor=args.temp_floor,
                floor_rad=args.floor_rad, hmin=args.hmin, hmax=args.hmax,
                t_start_days=args.t_start,
                r_max_cm=combined['r'].max())
    write_x_input_profile(str(ip_path), combined)

    # Also copy input_profile to stella dir
    import shutil
    shutil.copy2(str(ip_path), str(stella_dir / 'input_profile'))

    # Step 6: Generate plots
    print(f"\n--- Step 6: Generating plots ---")
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
        r_ej_inner = ejecta['r'].min() * 0.99
        r_max_cm = combined['r'][-1]

        # --- Density profile (upper-left) with shell/wind/ejecta vspans ---
        ax_rho_dens.plot(r_s, combined['rho'], 'k-', linewidth=1.5)
        ax_rho_dens.axvspan(r_shell_lo, r_shell_hi, color='red', alpha=0.15, label='shell')
        if args.no_wind:
            ax_rho_dens.axvspan(r_shell_hi, r_ej_inner, color='orange', alpha=0.10, label='buffer')
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

        if args.no_wind:
            plt.suptitle(f'X Model Initial Conditions (shell + {args.source_age}d ejecta, no wind)', fontsize=14)
        else:
            plt.suptitle(f'X Model Initial Conditions (shell + wind + {args.source_age}d ejecta)', fontsize=14)
        plt.tight_layout()
        plt.savefig(str(plots_dir / 'x_model_initial.png'), dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  Saved: {plots_dir}/x_model_initial.png")
    except Exception as e:
        print(f"  Plot generation failed: {e}")

    print(f"\n=== X model preparation complete ===")
    print(f"  STELLA input files in: {modmake_dir}")
    print(f"  Run STELLA with: python scripts/run_model_grid.py {output_dir}/config.yaml")


if __name__ == '__main__':
    main()
