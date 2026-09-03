#!/usr/bin/env python3
"""
Generate diagnostic plots for STELLA/CMFGEN models at each major step.

Produces:
  1. Structure plot: r, rho, v, T profiles (4 panels)
  2. Composition plot: mass fractions vs radius (H, He, C, O, Si, Fe)
  3. Lightcurve plot: luminosity vs time (from .swd or .tt file)
  4. Mass fraction distribution: mass fractions vs enclosed mass
  5. Structure evolution: rho, T, v, L profiles at multiple epochs
  6. Composition evolution: H, He, O, Fe profiles at multiple epochs

Usage:
    python plot_model.py <model_dir> [--step <step_name>] [--output <plots_dir>]

Steps for SN_A:
  - presn: KEPLER pre-SN model (before explosion)
  - stella_input: STELLA .hyd input (after KEPLER conversion)
  - stella_output: STELLA .swd output (after evolution)
  - stella_evolution: multi-epoch structure & composition evolution from .swd
  - homologous: homologously evolved snapshot
  - cmfgen_input: CMFGEN SN_HYDRO_DATA (before CMFGEN run)

Steps for SN_B:
  - sn_b_input: SN_B initial model (H/He shell + SN_A CSM)
  - stella_output: STELLA .swd output
  - cmfgen_input: CMFGEN SN_HYDRO_DATA
"""

import argparse
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

# Import read_swd from the conversion script
sys.path.insert(0, str(Path(__file__).parent))
from convert_stella_swd_to_cmfgen import read_swd, normalize_swd_composition, sigma_velocity
from convert_kepler_presn_to_stella import read_kepler_presn, build_stella_composition


MSUN = 1.989e33
K_B = 1.3807e-16
M_H = 1.6735e-24
A_RAD = 7.5657e-15


def plot_structure(r, rho, v, T, title, output_file, enclosed_mass=None,
                   ref_profile=None, x_axis='radius'):
    """Plot structure: density, velocity, temperature, enclosed mass (4 panels).

    x_axis: 'radius' (x-axis in cm, log scale) or 'mass' (x-axis in Msun).
            When x_axis='mass', enclosed_mass must be provided.

    ref_profile: optional dict with keys 'r', 'rho', 'T', 'v' for a reference
                 profile (e.g. pre-smoothing) to overlay as dashed lines.
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(title, fontsize=14, fontweight='bold')

    has_ref = ref_profile is not None

    if x_axis == 'mass':
        if enclosed_mass is None:
            raise ValueError("enclosed_mass required when x_axis='mass'")
        x = enclosed_mass / MSUN
        x_label = 'Enclosed Mass (M$_\\odot$)'
        # For mass axis, use linear scale (mass is naturally linear)
        plot_x = lambda ax, xx, yy, *args, **kw: ax.plot(xx, yy, *args, **kw)
        grid_which = 'major'
    else:
        x = r
        x_label = 'Radius (cm)'
        plot_x = lambda ax, xx, yy, *args, **kw: ax.loglog(xx, yy, *args, **kw)
        grid_which = 'both'

    # Panel 1: Density
    ax = axes[0, 0]
    if has_ref and x_axis == 'radius':
        ax.loglog(ref_profile['r'], ref_profile['rho'], 'b--',
                  linewidth=0.8, alpha=0.5, label='Pre-smoothing')
        ax.loglog(r, rho, 'b-', linewidth=1.2, label='Post-smoothing')
        ax.legend(fontsize=8)
    else:
        plot_x(ax, x, rho, 'b-', linewidth=0.8)
    ax.set_xlabel(x_label)
    ax.set_ylabel('Density (g/cm³)')
    ax.set_title('Density Profile')
    ax.grid(True, alpha=0.3, which=grid_which)

    # Panel 2: Velocity
    ax = axes[0, 1]
    if has_ref and x_axis == 'radius':
        ax.semilogx(ref_profile['r'], ref_profile['v'], 'r--',
                    linewidth=0.8, alpha=0.5, label='Pre-smoothing')
        ax.semilogx(r, v, 'r-', linewidth=1.2, label='Post-smoothing')
        ax.legend(fontsize=8)
    else:
        if x_axis == 'mass':
            ax.plot(x, v, 'r-', linewidth=0.8)
        else:
            ax.semilogx(x, v, 'r-', linewidth=0.8)
    ax.set_xlabel(x_label)
    ax.set_ylabel('Velocity (km/s)')
    ax.set_title('Velocity Profile')
    ax.grid(True, alpha=0.3, which=grid_which)

    # Panel 3: Temperature
    ax = axes[1, 0]
    if has_ref and x_axis == 'radius':
        ax.loglog(ref_profile['r'], ref_profile['T'], 'g--',
                  linewidth=0.8, alpha=0.5, label='Pre-smoothing')
        ax.loglog(r, T, 'g-', linewidth=1.2, label='Post-smoothing')
        ax.legend(fontsize=8)
    else:
        plot_x(ax, x, T, 'g-', linewidth=0.8)
    ax.set_xlabel(x_label)
    ax.set_ylabel('Temperature (K)')
    ax.set_title('Temperature Profile')
    ax.grid(True, alpha=0.3, which=grid_which)

    # Panel 4: Enclosed mass vs radius, or density vs velocity
    ax = axes[1, 1]
    if x_axis == 'radius' and enclosed_mass is not None:
        ax.semilogx(r, enclosed_mass / MSUN, 'm-', linewidth=0.8)
        ax.set_xlabel('Radius (cm)')
        ax.set_ylabel('Enclosed Mass (M$_\\odot$)')
        ax.set_title('Mass Distribution')
    elif x_axis == 'mass':
        # When x-axis is already mass, plot density vs velocity instead
        ax.loglog(v, rho, 'm-', linewidth=0.8)
        ax.set_xlabel('Velocity (km/s)')
        ax.set_ylabel('Density (g/cm³)')
        ax.set_title('Density vs Velocity')
    else:
        ax.loglog(v, rho, 'm-', linewidth=0.8)
        ax.set_xlabel('Velocity (km/s)')
        ax.set_ylabel('Density (g/cm³)')
        ax.set_title('Density vs Velocity')
    ax.grid(True, alpha=0.3, which='both')

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Structure plot: {output_file}")


def plot_sigma(r, v_kms, title, output_file, enclosed_mass=None):
    """Plot the CMFGEN sigma profile: sigma = dlnV/dlnr - 1.

    CMFGEN (rv_sn_model_sniin.f:241) computes sigma from the velocity and
    radius grid using a distance-weighted linear gradient.  This plot shows
    exactly what CMFGEN will see when it reads the SN_HYDRO_DATA file.

    Physical interpretation:
      sigma = 0  : homologous expansion (v ∝ r, free expansion)
      sigma > 0  : accelerating (v increases faster than r)
      sigma < 0  : decelerating (v increases slower than r)
      sigma = -1 : constant velocity (flat, no expansion gradient)

    CMFGEN's Sobolev line formation assumes |sigma| is not too large.
    Values |sigma| > 2-3 can cause convergence problems (negative sigma
    makes the Sobolev integral ill-conditioned; very large positive sigma
    makes the continuum source function unstable).

    The plot has two panels:
      Top:    sigma vs radius (log-log), with horizontal lines at sigma=0
              (homologous) and sigma=-1 (constant velocity).
      Bottom: velocity vs radius (log-log) for context, with a homologous
              reference line v = r/t for comparison.
    """
    sigma = sigma_velocity(r, v_kms)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), sharex=True)
    fig.suptitle(title, fontsize=14, fontweight='bold')

    # --- Top panel: sigma vs radius ---
    ax1.semilogx(r, sigma, 'b-', linewidth=0.8)
    ax1.axhline(0.0, color='k', linestyle='-', linewidth=0.5, alpha=0.5,
                label='σ = 0 (homologous)')
    ax1.axhline(-1.0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5,
                label='σ = -1 (constant v)')
    # Mark regions where |sigma| > 2 (potential CMFGEN convergence issues)
    bad = np.abs(sigma) > 2.0
    if np.any(bad):
        ax1.fill_between(r, sigma, 0, where=bad, alpha=0.15, color='red',
                         label=f'|σ| > 2 ({np.sum(bad)} zones)')
    ax1.set_ylabel('σ = dlnV/dlnr − 1')
    ax1.set_title('CMFGEN Sigma Profile (velocity gradient exponent)')
    ax1.set_ylim(min(-3, sigma.min() - 0.5), max(3, sigma.max() + 0.5))
    ax1.legend(fontsize=9, loc='best')
    ax1.grid(True, alpha=0.3, which='both')

    # Annotate key regions
    # Find the shell-wind boundary (where v drops sharply)
    dv = np.abs(np.diff(v_kms))
    if len(dv) > 0:
        idx_max_dv = np.argmax(dv)
        ax1.axvline(r[idx_max_dv], color='orange', linestyle=':', alpha=0.6,
                    linewidth=1.0)
        ax1.annotate(f'shell-wind\n(v drop)', xy=(r[idx_max_dv], sigma[idx_max_dv]),
                     fontsize=7, color='orange', ha='center',
                     xytext=(r[idx_max_dv] * 1.5, sigma[idx_max_dv] + 1),
                     arrowprops=dict(arrowstyle='->', color='orange', alpha=0.6))

    # --- Bottom panel: velocity vs radius with homologous reference ---
    ax2.loglog(r, v_kms, 'r-', linewidth=0.8, label='v(r)')
    # Homologous reference: v_hom = r / t where t is inferred from the
    # outermost point (which should be close to free expansion)
    t_ref = r[-1] / (v_kms[-1] * 1e5)  # seconds
    v_hom = r / t_ref * 1e-5  # km/s
    ax2.loglog(r, v_hom, 'k--', linewidth=0.5, alpha=0.5,
               label=f'v = r/t (homologous, t={t_ref/86400:.1f}d)')
    ax2.set_xlabel('Radius (cm)')
    ax2.set_ylabel('Velocity (km/s)')
    ax2.set_title('Velocity Profile')
    ax2.legend(fontsize=9, loc='best')
    ax2.grid(True, alpha=0.3, which='both')

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(output_file, dpi=140, bbox_inches='tight')
    plt.close()
    print(f"  Sigma plot: {output_file}")

    # Print summary statistics
    print(f"    sigma range: [{sigma.min():.2f}, {sigma.max():.2f}]")
    print(f"    sigma mean:  {sigma.mean():.3f}")
    n_bad = np.sum(np.abs(sigma) > 2.0)
    if n_bad > 0:
        print(f"    WARNING: {n_bad} zones with |sigma| > 2.0 "
              f"(may cause CMFGEN convergence issues)")
    else:
        print(f"    All zones |sigma| <= 2.0 (OK for CMFGEN)")


def plot_composition(r, comp_dict, title, output_file, x_axis='radius'):
    """Plot composition: mass fractions vs radius or enclosed mass.

    comp_dict: {'H': array, 'He': array, 'C': array, 'O': array, 'Si': array, 'Fe': array}
    x_axis: 'radius' or 'mass' or 'velocity'

    When x_axis='radius', the x-axis is radius in cm, log scale.
    """
    fig, ax = plt.subplots(figsize=(12, 6))
    fig.suptitle(title, fontsize=14, fontweight='bold')

    colors = {'H': 'blue', 'He': 'green', 'C': 'orange', 'O': 'red',
              'Si': 'purple', 'Fe': 'brown', 'Ni': 'gray', 'Ca': 'pink'}

    if x_axis == 'radius':
        x = r
        ax.set_xlabel('Radius (cm)')
        ax.set_xscale('log')
    elif x_axis == 'velocity':
        ax.set_xlabel('Velocity (km/s)')
    elif x_axis == 'mass':
        ax.set_xlabel('Enclosed Mass (Msun)')

    for elem, frac in comp_dict.items():
        if np.any(frac > 1e-10):
            ax.plot(x, frac, label=elem, color=colors.get(elem, 'black'), linewidth=1.0)

    ax.set_ylabel('Mass Fraction')
    ax.set_ylim(-0.05, 1.05)
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Composition plot: {output_file}")


def plot_composition_vs_mass(r, rho, comp_dict, title, output_file,
                             shaded_mass_msun=1.5, m_enc=None):
    """Plot mass fractions vs enclosed mass.

    The inner `shaded_mass_msun` Msun of material is shaded to highlight
    the region near the compact remnant / mass cut.

    x-axis: Enclosed Mass (Msun)
    y-axis: Mass Fraction (log scale)

    If m_enc is provided (e.g. from .hyd file, including mass cut), it is
    used directly. Otherwise it is computed from r and rho (ejecta only).
    """
    if m_enc is None:
        # Compute enclosed mass from zone masses
        r_lo = np.concatenate([[r[0] * 0.9], r[:-1]])
        dm = (4/3) * np.pi * (r**3 - r_lo**3) * rho
        dm = np.maximum(dm, 0)
        m_enc = np.cumsum(dm) / MSUN

    fig, ax = plt.subplots(figsize=(12, 6))
    fig.suptitle(title, fontsize=14, fontweight='bold')

    # Shade the inner shaded_mass_msun Msun region (mass cut / remnant)
    # Only shade if shaded_mass_msun > 0 (pre-SN model still has the remnant)
    if shaded_mass_msun > 0:
        shade_hi = min(shaded_mass_msun, m_enc[-1])
        if shade_hi > m_enc[0]:
            ax.axvspan(m_enc[0], shade_hi, color='gray', alpha=0.2,
                       label=f'Inner {shaded_mass_msun} M$_\\odot$ (mass cut)')

    colors = {'H': 'blue', 'He': 'green', 'C': 'orange', 'O': 'red',
              'Si': 'purple', 'Fe': 'brown', 'Ni': 'gray', 'Ca': 'pink'}

    for elem, frac in comp_dict.items():
        if np.any(frac > 1e-10):
            ax.plot(m_enc, frac, label=elem, color=colors.get(elem, 'black'),
                    linewidth=1.0)

    ax.set_xlabel('Enclosed Mass (M$_\\odot$)')
    ax.set_ylabel('Mass Fraction')
    ax.set_ylim(-0.05, 1.05)
    ax.legend(loc='best', fontsize=9)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Mass fraction distribution: {output_file}")


def plot_lightcurve(swd_file, title, output_file):
    """Plot lightcurve: luminosity vs time from .swd file."""
    epochs = read_swd(swd_file)
    times = sorted(epochs.keys())
    lum = []
    for t in times:
        e = epochs[t]
        # Flux is in erg/s; convert to Lsun
        L = np.sum(e['flux']) / 3.828e33  # Lsun
        lum.append(L)

    fig, ax = plt.subplots(figsize=(10, 6))
    fig.suptitle(title, fontsize=14, fontweight='bold')
    ax.semilogy(times, lum, 'b.-', markersize=4, linewidth=1.0)
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Luminosity (Lsun)')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Lightcurve plot: {output_file}")


def read_tt_file(tt_path):
    """Read a STELLA .tt file and extract multiband lightcurve data.

    Returns dict with 'time' (days) and absolute magnitudes:
    Mbol, MU, MB, MV, MI, MR, Mbolavg.
    Only post-SBO (t > 0) rows are returned.
    """
    times, Mbol, MU, MB, MV, MI, MR, Mbolavg = [], [], [], [], [], [], [], []
    in_data = False
    with open(tt_path, 'r', errors='replace') as f:
        for line in f:
            if line.strip().startswith('time') and 'Mbol' in line:
                in_data = True
                continue
            if not in_data:
                continue
            parts = line.split()
            if len(parts) < 13:
                continue
            try:
                t = float(parts[0])
            except ValueError:
                continue
            if t <= 0:
                continue
            times.append(t)
            Mbol.append(float(parts[6]))
            MU.append(float(parts[7]))
            MB.append(float(parts[8]))
            MV.append(float(parts[9]))
            MI.append(float(parts[10]))
            MR.append(float(parts[11]))
            Mbolavg.append(float(parts[12]))
    return {
        'time': np.array(times),
        'Mbol': np.array(Mbol),
        'MU': np.array(MU),
        'MB': np.array(MB),
        'MV': np.array(MV),
        'MI': np.array(MI),
        'MR': np.array(MR),
        'Mbolavg': np.array(Mbolavg),
    }


def read_2011dh_photometry(obs_dir):
    """Read SN 2011dh UBVRI photometry from Ergon et al. (2015) tablec4.dat.

    Returns dict with 'phase' (days since explosion) and apparent mags
    for U, B, V, R, I bands (with NaN for missing data).
    SN 2011dh: distance modulus = 29.62, A_V = 0.14 (Ergon et al. 2015).
    Extinction: A_U=0.21, A_B=0.16, A_V=0.14, A_R=0.11, A_I=0.09.

    Table format (fixed-width, bytes per ReadMe):
      1-8:   Epoch (JD-2400000)
      11-16: Phase (days)
      19-24: U mag
      27-32: e_U
      35-40: B mag
      43-48: e_B
      51-56: V mag
      59-64: e_V
      67-72: R mag
      75-80: e_R
      83-88: I mag
      91-96: e_I
    """
    import os
    table = os.path.join(obs_dir, 'tablec4.dat')
    phases, U, B, V, R, I = [], [], [], [], [], []
    with open(table) as f:
        for line in f:
            if len(line) < 16:
                continue
            try:
                phase = float(line[10:16])
            except ValueError:
                continue
            phases.append(phase)
            def _parse(s):
                s = s.strip()
                return float(s) if s else np.nan
            U.append(_parse(line[18:24]))
            B.append(_parse(line[34:40]))
            V.append(_parse(line[50:56]))
            R.append(_parse(line[66:72]))
            I.append(_parse(line[82:88]))
    # Distance modulus and extinction for M51/SN 2011dh
    dm = 29.62  # distance modulus (7.4 Mpc)
    ext = {'U': 0.21, 'B': 0.16, 'V': 0.14, 'R': 0.11, 'I': 0.09}
    phases = np.array(phases)
    result = {'phase': phases}
    for band, arr, ext_val in [('U', U, ext['U']), ('B', B, ext['B']),
                                ('V', V, ext['V']), ('R', R, ext['R']),
                                ('I', I, ext['I'])]:
        result[band] = np.array(arr) - dm - ext_val  # absolute magnitude
    return result


def read_2011dh_bolometric(obs_dir):
    """Read SN 2011dh bolometric luminosity from Ergon et al. (2015) tablec9.dat.

    Returns dict with 'phase' (days) and 'Mbol' (absolute bolometric magnitude).
    Luminosity in table is in units of 1e41 erg/s.
    """
    import os
    table = os.path.join(obs_dir, 'tablec9.dat')
    data = np.genfromtxt(table, usecols=[1, 2],
                         names=['phase', 'Lbol'])
    # Convert L (1e41 erg/s) to absolute bolometric magnitude
    # M_bol = -2.5 * log10(L/Lsun) + M_bol_sun
    # Lsun = 3.828e33 erg/s, M_bol_sun = 4.74
    Lsun = 3.828e33
    L_erg_s = data['Lbol'] * 1e41
    Mbol = -2.5 * np.log10(L_erg_s / Lsun) + 4.74
    return {'phase': data['phase'], 'Mbol': Mbol}


def plot_multiband_lightcurve(tt_file, model_name, output_file,
                              obs_dir=None):
    """Plot multiband (U,B,V,R,I) lightcurve from STELLA .tt file,
    compared to SN 2011dh observations (if available).

    STELLA .tt has absolute magnitudes. Observations are converted to
    absolute magnitudes using dm=29.62 and Galactic+host extinction.
    """
    stella = read_tt_file(tt_file)
    t = stella['time']

    bands = ['U', 'B', 'V', 'R', 'I']
    stella_keys = {'U': 'MU', 'B': 'MB', 'V': 'MV', 'R': 'MR', 'I': 'MI'}
    band_colors = {'U': 'purple', 'B': 'blue', 'V': 'green',
                   'R': 'orange', 'I': 'red'}

    # Read observations if available
    obs = None
    obs_bol = None
    if obs_dir is not None:
        import os
        if os.path.exists(os.path.join(obs_dir, 'tablec4.dat')):
            obs = read_2011dh_photometry(obs_dir)
        if os.path.exists(os.path.join(obs_dir, 'tablec9.dat')):
            obs_bol = read_2011dh_bolometric(obs_dir)

    n_panels = len(bands) + 1  # +1 for bolometric
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    fig.suptitle(f'{model_name} - Multiband Lightcurve vs SN 2011dh',
                 fontsize=14, fontweight='bold')
    axes_flat = axes.flatten()

    # Compute x/y limits from STELLA model data only
    t_min, t_max = t.min(), t.max()
    x_pad = 0.05 * (t_max - t_min)
    xlim = (max(t_min - x_pad, 0), t_max + x_pad)

    # Collect all STELLA magnitude ranges for shared y-limit per panel
    stella_mags = {
        'bol': np.concatenate([stella['Mbol'], stella['Mbolavg']]),
        'U': stella['MU'], 'B': stella['MB'], 'V': stella['MV'],
        'R': stella['MR'], 'I': stella['MI'],
    }
    y_ranges = {}
    for k, vals in stella_mags.items():
        y_ranges[k] = (vals.min(), vals.max())

    def _ylim(key, pad_frac=0.1):
        lo, hi = y_ranges[key]
        pad = pad_frac * (hi - lo)
        # Return (bottom, top) in normal order; invert_yaxis() flips to
        # put bright (low mag) at top
        return (lo - pad, hi + pad)

    # Bolometric panel (top-left)
    ax = axes_flat[0]
    ax.plot(t, stella['Mbol'], 'k-', linewidth=1.5, label='STELLA Mbol')
    ax.plot(t, stella['Mbolavg'], 'k--', linewidth=1.0, alpha=0.5,
            label='STELLA Mbolavg')
    if obs_bol is not None:
        ax.plot(obs_bol['phase'], obs_bol['Mbol'], 'ko', markersize=4,
                alpha=0.6, label='SN 2011dh (Ergon+2015)')
    ax.set_xlim(*xlim)
    ax.set_ylim(*_ylim('bol'))
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Absolute M_bol')
    ax.set_title('Bolometric')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.invert_yaxis()

    # Band panels
    for i, band in enumerate(bands):
        ax = axes_flat[i + 1]
        key = stella_keys[band]
        ax.plot(t, stella[key], '-', color=band_colors[band],
                linewidth=1.5, label=f'STELLA {band}')
        if obs is not None:
            obs_mag = obs[band]
            valid = ~np.isnan(obs_mag)
            ax.plot(obs['phase'][valid], obs_mag[valid], 'o',
                    color=band_colors[band], markersize=4, alpha=0.6,
                    label=f'SN 2011dh {band}')
        ax.set_xlim(*xlim)
        ax.set_ylim(*_ylim(band))
        ax.set_xlabel('Time (days)')
        ax.set_ylabel(f'Absolute M_{band}')
        ax.set_title(f'{band} band')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        ax.invert_yaxis()

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Multiband lightcurve plot: {output_file}")


def plot_structure_evolution(swd_file, title, output_file, n_epochs=8):
    """Plot structure profile evolution: rho, T, v, L vs radius at multiple epochs.

    Reads post-SBO (t > 0) epochs from the .swd file, selects n_epochs
    evenly-spaced in log time, and overlays the profiles in a 2x2 panel
    with a color gradient from early (blue) to late (red). Axis limits
    are set to the data range across all selected epochs.
    """
    epochs = read_swd(swd_file)
    all_times = sorted(epochs.keys())
    # Only post-SBO epochs (t > 0)
    times = [t for t in all_times if t > 0]
    if len(times) < 2:
        print(f"  WARNING: Only {len(times)} post-SBO epochs in .swd, need >= 2 for evolution plot")
        return

    # Select n_epochs evenly-spaced in log time
    t_min, t_max = times[0], times[-1]
    t_sel = np.logspace(np.log10(t_min), np.log10(t_max), n_epochs)
    # Snap to nearest available epoch
    sel_times = [min(times, key=lambda t: abs(t - ts)) for ts in t_sel]
    # Remove duplicates while preserving order
    seen = set()
    sel_times = [t for t in sel_times if not (t in seen or seen.add(t))]

    # Color gradient: blue (early) -> red (late)
    cmap = plt.cm.coolwarm
    colors = [cmap(i / max(len(sel_times) - 1, 1)) for i in range(len(sel_times))]

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle(title, fontsize=14, fontweight='bold')

    # Track data ranges for axis limits
    r_min = r_max = None
    rho_min = rho_max = None
    T_min = T_max = None
    v_min = v_max = None
    L_min = L_max = None

    for i, t in enumerate(sel_times):
        e = epochs[t]
        # .swd is outer->inner; reverse for plotting inner->outer
        r = e['r'][::-1]
        rho = e['rho'][::-1]
        v = e['v'][::-1]
        T = e['T'][::-1]
        flux = e['flux'][::-1]  # already erg/s per zone (luminosity)
        label = f't={t:.2f}d'
        lw = 2.0 if i == len(sel_times) - 1 else 1.0
        zorder = i + 1

        axes[0, 0].loglog(r, rho, color=colors[i], linewidth=lw, label=label, zorder=zorder)
        axes[0, 1].loglog(r, T, color=colors[i], linewidth=lw, label=label, zorder=zorder)
        axes[1, 0].semilogx(r, v, color=colors[i], linewidth=lw, label=label, zorder=zorder)
        axes[1, 1].loglog(r, flux, color=colors[i], linewidth=lw, label=label, zorder=zorder)

        # Update data ranges (only positive values for log axes)
        r_min = r.min() if r_min is None else min(r_min, r.min())
        r_max = r.max() if r_max is None else max(r_max, r.max())
        rho_pos = rho[rho > 0]
        if len(rho_pos) > 0:
            rho_min = rho_pos.min() if rho_min is None else min(rho_min, rho_pos.min())
            rho_max = rho_pos.max() if rho_max is None else max(rho_max, rho_pos.max())
        T_pos = T[T > 0]
        if len(T_pos) > 0:
            T_min = T_pos.min() if T_min is None else min(T_min, T_pos.min())
            T_max = T_pos.max() if T_max is None else max(T_max, T_pos.max())
        v_min = v.min() if v_min is None else min(v_min, v.min())
        v_max = v.max() if v_max is None else max(v_max, v.max())
        flux_pos = flux[flux > 0]
        if len(flux_pos) > 0:
            L_min = flux_pos.min() if L_min is None else min(L_min, flux_pos.min())
            L_max = flux_pos.max() if L_max is None else max(L_max, flux_pos.max())

    # Set axis limits with small padding
    def _log_lim(lo, hi, pad=0.05):
        """Add padding in log space."""
        if lo is None or hi is None or lo <= 0 or hi <= 0:
            return None
        log_lo, log_hi = np.log10(lo), np.log10(hi)
        span = log_hi - log_lo
        return (10**(log_lo - pad * span), 10**(log_hi + pad * span))

    if r_min and r_max:
        xlim = _log_lim(r_min, r_max, pad=0.02)
        for ax_row in axes:
            for ax in ax_row:
                ax.set_xlim(*xlim)

    if rho_min and rho_max:
        axes[0, 0].set_ylim(*_log_lim(rho_min, rho_max))
    if T_min and T_max:
        axes[0, 1].set_ylim(*_log_lim(T_min, T_max))
    if v_min is not None and v_max is not None:
        axes[1, 0].set_ylim(v_min * 0.9, v_max * 1.05)
    if L_min and L_max:
        axes[1, 1].set_ylim(*_log_lim(L_min, L_max))

    axes[0, 0].set_xlabel('Radius (cm)')
    axes[0, 0].set_ylabel('Density (g cm$^{-3}$)')
    axes[0, 0].set_title('Density Profile')
    axes[0, 0].legend(fontsize=7, loc='best')
    axes[0, 0].grid(True, alpha=0.3, which='both')

    axes[0, 1].set_xlabel('Radius (cm)')
    axes[0, 1].set_ylabel('Temperature (K)')
    axes[0, 1].set_title('Temperature Profile')
    axes[0, 1].grid(True, alpha=0.3, which='both')

    axes[1, 0].set_xlabel('Radius (cm)')
    axes[1, 0].set_ylabel('Velocity (km/s)')
    axes[1, 0].set_title('Velocity Profile')
    axes[1, 0].grid(True, alpha=0.3)

    axes[1, 1].set_xlabel('Radius (cm)')
    axes[1, 1].set_ylabel('Luminosity (erg s$^{-1}$)')
    axes[1, 1].set_title('Bolometric Luminosity Profile')
    axes[1, 1].grid(True, alpha=0.3, which='both')

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Structure evolution plot: {output_file}")


def plot_composition_evolution(swd_file, title, output_file, n_epochs=6):
    """Plot composition evolution: mass fractions vs radius at multiple epochs.

    Shows H, He, O, and Fe (or Ni) profiles at selected post-SBO epochs to
    illustrate how the composition structure evolves with time. Axis limits
    are set to the data range across all selected epochs.
    """
    epochs = read_swd(swd_file)
    all_times = sorted(epochs.keys())
    # Only post-SBO epochs (t > 0)
    times = [t for t in all_times if t > 0]
    if len(times) < 2:
        print(f"  WARNING: Only {len(times)} post-SBO epochs in .swd, need >= 2 for evolution plot")
        return

    # Select epochs (fewer than structure plot since composition is noisier)
    t_min, t_max = times[0], times[-1]
    t_sel = np.logspace(np.log10(t_min), np.log10(t_max), n_epochs)
    sel_times = [min(times, key=lambda t: abs(t - ts)) for ts in t_sel]
    seen = set()
    sel_times = [t for t in sel_times if not (t in seen or seen.add(t))]

    cmap = plt.cm.coolwarm
    colors = [cmap(i / max(len(sel_times) - 1, 1)) for i in range(len(sel_times))]

    # Elements to plot (4 panels)
    elem_names = ['H', 'He', 'O', 'Fe']
    elem_labels = {'H': 'H', 'He': 'He', 'O': 'O', 'Fe': 'Fe (incl. Ni56 decay)'}

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle(title, fontsize=14, fontweight='bold')
    axes_flat = axes.flatten()

    # Track radius range for shared x limits
    r_min = r_max = None

    for i, t in enumerate(sel_times):
        e = epochs[t]
        comp_raw = e.get('composition')
        if comp_raw is None:
            continue
        # Normalize to mass fractions
        comp_norm = normalize_swd_composition(comp_raw)
        r = e['r'][::-1]  # inner->outer
        label = f't={t:.2f}d'
        lw = 2.0 if i == len(sel_times) - 1 else 1.0

        r_min = r.min() if r_min is None else min(r_min, r.min())
        r_max = r.max() if r_max is None else max(r_max, r.max())

        for j, elem in enumerate(elem_names):
            ax = axes_flat[j]
            frac = comp_norm[elem][::-1]  # inner->outer
            ax.semilogx(r, frac, color=colors[i], linewidth=lw, label=label)

    # Set consistent axis limits
    if r_min and r_max:
        log_lo, log_hi = np.log10(r_min), np.log10(r_max)
        span = log_hi - log_lo
        xlim = (10**(log_lo - 0.02 * span), 10**(log_hi + 0.02 * span))
    for j in range(len(elem_names)):
        ax = axes_flat[j]
        ax.set_xlabel('Radius (cm)')
        ax.set_ylabel(f'{elem_labels[elem_names[j]]} mass fraction')
        ax.set_title(f'{elem_names[j]} Profile')
        ax.set_ylim(-0.05, 1.05)
        if r_min and r_max:
            ax.set_xlim(*xlim)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=7, loc='best')

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Composition evolution plot: {output_file}")


def read_hyd_file(hyd_path):
    """Read a STELLA .hyd file and return r, rho, T, v, enclosed_mass.

    .hyd row format (8 columns):
      col 0: zone index (1-based)
      col 1: dMr (zone mass, Msun)
      col 2: r (outer radius, cm)
      col 3: rho (density, g/cc)
      col 4: T (temperature, K)
      col 5: v (velocity, cm/s)
      col 6: Mr (enclosed mass including mass cut, Msun)
      col 7: placeholder (0.0)
    """
    with open(hyd_path) as f:
        header = f.readline()
        cols = header.split()
        nzon = int(cols[1])
        data = np.loadtxt(hyd_path, skiprows=1)
    r = data[:, 2]            # radius (cm)
    rho = data[:, 3]          # density (g/cc)
    T = data[:, 4]            # temperature (K)
    v = data[:, 5] / 1e5      # velocity (km/s, converted from cm/s)
    M_r = data[:, 6] * MSUN   # enclosed mass (g, converted from Msun)
    return r, rho, T, v, M_r


def read_abn_file(abn_path):
    """Read a STELLA .abn file and return composition dict."""
    data = np.loadtxt(abn_path)
    # Columns: km, dum, dum, dum, H, He, C, N, O, Ne, Na, Mg, Al, Si, S, Ar, Ca, Fe, Ni58, Ni56
    return {
        'H': data[:, 4],
        'He': data[:, 5],
        'C': data[:, 6],
        'O': data[:, 8],
        'Si': data[:, 13],
        'Ca': data[:, 16],
        'Fe': data[:, 17],
        'Ni': data[:, 19],
    }


def read_swd_composition(swd_path, t_day):
    """Read composition from .swd file at a specific epoch."""
    epochs = read_swd(swd_path)
    available = sorted(epochs.keys())
    closest_t = min(available, key=lambda t: abs(t - t_day))
    e = epochs[closest_t]
    if e.get('composition') is not None:
        comp = e['composition']
        # STELLA SWD composition: H, He, C, N, O, Ne, Na, Mg, Al, Si, S, Ar, Ca, Fe, Ni
        return {
            'H': comp[:, 0],
            'He': comp[:, 1],
            'C': comp[:, 2],
            'O': comp[:, 4],
            'Si': comp[:, 10],
            'Ca': comp[:, 12],
            'Fe': comp[:, 13],
            'Ni': comp[:, 14],
        }, e
    return None, e


def plot_step(model_dir, step_name, plots_dir, swd_day=None):
    """Generate all plots for a given step."""
    model_dir = Path(model_dir)
    plots_dir = Path(plots_dir)
    plots_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n  Generating plots for step: {step_name}")

    if step_name == 'presn':
        # KEPLER pre-SN model — look for any *presn* file in the model directory
        presn_file = None
        for candidate in model_dir.glob('*presn*'):
            if candidate.is_file():
                presn_file = candidate
                break
        if presn_file is None:
            print(f"  WARNING: No pre-SN file found in {model_dir}")
            return
        # Read KEPLER pre-SN file using the conversion script's reader
        kep = read_kepler_presn(str(presn_file))
        r = kep['r_cm']
        rho = kep['rho']
        T = kep['T_K']
        v = kep['v_cms'] / 1e5  # km/s
        M_r = kep['m_g']  # enclosed mass in g
        plot_structure(r, rho, v, T, f'{model_dir.name} - KEPLER Pre-SN ({presn_file.name})',
                       str(plots_dir / f'{step_name}_structure.png'), enclosed_mass=M_r)
        # Composition from KEPLER isotopes
        fH, fHe, fC, fO, fSi, fCa, fFe, fNi56 = build_stella_composition(kep)
        comp = {'H': fH, 'He': fHe, 'C': fC, 'O': fO, 'Si': fSi, 'Ca': fCa, 'Fe': fFe}
        plot_composition(r, comp, f'{model_dir.name} - Composition (KEPLER Pre-SN)',
                        str(plots_dir / f'{step_name}_composition.png'))
        plot_composition_vs_mass(r, rho, comp,
                                 f'{model_dir.name} - Mass Fractions (KEPLER Pre-SN)',
                                 str(plots_dir / f'{step_name}_mass_fractions.png'),
                                 shaded_mass_msun=1.4, m_enc=M_r / MSUN)

    elif step_name == 'stella_input':
        # STELLA .hyd input
        hyd = model_dir / 'stella' / 'modmake' / f'{model_dir.name}.hyd'
        abn = model_dir / 'stella' / 'modmake' / f'{model_dir.name}.abn'
        if not hyd.exists():
            print(f"  WARNING: .hyd not found: {hyd}")
            return
        r, rho, T, v, M_r = read_hyd_file(str(hyd))

        # Check for pre-smoothing profile (saved by convert_kepler_presn_to_stella.py
        # when --smooth-profile is enabled)
        pre_smooth_hyd = model_dir / 'stella' / 'modmake' / f'{model_dir.name}.hyd.pre_smooth'
        ref_profile = None
        if pre_smooth_hyd.exists():
            print(f"  Found pre-smoothing profile: {pre_smooth_hyd.name}")
            r_ref, rho_ref, T_ref, v_ref, _ = read_hyd_file(str(pre_smooth_hyd))
            ref_profile = {'r': r_ref, 'rho': rho_ref, 'T': T_ref, 'v': v_ref}

        title = f'{model_dir.name} - STELLA Input (.hyd)'
        if ref_profile is not None:
            title += ' [pre/post-smoothing]'
        # Plot 1: x-axis = radius (cm)
        plot_structure(r, rho, v, T, title,
                       str(plots_dir / f'{step_name}_structure.png'),
                       enclosed_mass=M_r, ref_profile=ref_profile, x_axis='radius')
        # Plot 2: x-axis = enclosed mass (Msun)
        plot_structure(r, rho, v, T, title + ' [vs Mass]',
                       str(plots_dir / f'{step_name}_structure_vs_mass.png'),
                       enclosed_mass=M_r, ref_profile=ref_profile, x_axis='mass')
        # Plot 3: CMFGEN sigma profile (velocity gradient exponent)
        plot_sigma(r, v, f'{model_dir.name} - CMFGEN Sigma Profile (STELLA Input)',
                   str(plots_dir / 'stella_sigma.png'),
                   enclosed_mass=M_r)
        if abn.exists():
            comp = read_abn_file(str(abn))
            plot_composition(r, comp, f'{model_dir.name} - Composition (STELLA Input)',
                            str(plots_dir / f'{step_name}_composition.png'))
            plot_composition_vs_mass(r, rho, comp,
                                     f'{model_dir.name} - Mass Fractions (STELLA Input)',
                                     str(plots_dir / f'{step_name}_mass_fractions.png'),
                                     shaded_mass_msun=0, m_enc=M_r / MSUN)

    elif step_name == 'stella_output':
        # STELLA .swd output
        swd = model_dir / 'stella' / 'strad' / f'{model_dir.name}.swd'
        if not swd.exists():
            swd = model_dir / 'stella' / f'{model_dir.name}.swd'
        if not swd.exists():
            print(f"  WARNING: .swd not found: {swd}")
            return
        plot_lightcurve(str(swd), f'{model_dir.name} - STELLA Lightcurve',
                       str(plots_dir / f'{step_name}_lightcurve.png'))
        # Multiband lightcurve vs SN 2011dh (from .tt file)
        tt = model_dir / 'stella' / f'{model_dir.name}.tt'
        if not tt.exists():
            tt = model_dir / 'stella' / 'strad' / f'{model_dir.name}.tt'
        if tt.exists():
            obs_dir = model_dir.parent / 'observational_data' / 'SN2011dh'
            if not obs_dir.exists():
                obs_dir = Path(__file__).resolve().parent.parent / 'observational_data' / 'SN2011dh'
            plot_multiband_lightcurve(str(tt), model_dir.name,
                                      str(plots_dir / f'{step_name}_multiband.png'),
                                      obs_dir=str(obs_dir) if obs_dir.exists() else None)
        else:
            print(f"  WARNING: .tt not found for multiband plot: {tt}")
        # Plot structure at final epoch
        epochs = read_swd(str(swd))
        t_final = sorted(epochs.keys())[-1]
        e = epochs[t_final]
        # .swd is sorted outer->inner; reverse to inner->outer for plotting
        r = e['r'][::-1]
        rho = e['rho'][::-1]
        v = e['v'][::-1]
        T = e['T'][::-1]
        plot_structure(r, rho, v, T,
                       f'{model_dir.name} - STELLA Output (t={t_final:.1f}d)',
                       str(plots_dir / f'{step_name}_structure.png'))
        # Composition at final epoch
        comp, _ = read_swd_composition(str(swd), t_final)
        if comp is not None:
            # Reverse composition arrays to match inner->outer order
            comp = {k: arr[::-1] for k, arr in comp.items()}
            plot_composition(r, comp,
                           f'{model_dir.name} - Composition (STELLA Output, t={t_final:.1f}d)',
                           str(plots_dir / f'{step_name}_composition.png'))
            plot_composition_vs_mass(r, rho, comp,
                                   f'{model_dir.name} - Mass Fractions (STELLA Output, t={t_final:.1f}d)',
                                   str(plots_dir / f'{step_name}_mass_fractions.png'),
                                   shaded_mass_msun=0)

    elif step_name == 'multiband':
        # Multiband (U,B,V,R,I) lightcurve vs SN 2011dh (from .tt file only)
        tt = model_dir / 'stella' / f'{model_dir.name}.tt'
        if not tt.exists():
            tt = model_dir / 'stella' / 'strad' / f'{model_dir.name}.tt'
        if not tt.exists():
            print(f"  WARNING: .tt not found for multiband plot: {tt}")
            return
        obs_dir = model_dir.parent / 'observational_data' / 'SN2011dh'
        if not obs_dir.exists():
            obs_dir = Path(__file__).resolve().parent.parent / 'observational_data' / 'SN2011dh'
        plot_multiband_lightcurve(str(tt), model_dir.name,
                                  str(plots_dir / 'stella_output_multiband.png'),
                                  obs_dir=str(obs_dir) if obs_dir.exists() else None)

    elif step_name == 'stella_evolution':
        # Multi-epoch structure and composition evolution from .swd
        swd = model_dir / 'stella' / 'strad' / f'{model_dir.name}.swd'
        if not swd.exists():
            swd = model_dir / 'stella' / f'{model_dir.name}.swd'
        if not swd.exists():
            print(f"  WARNING: .swd not found: {swd}")
            return
        epochs = read_swd(str(swd))
        times = sorted(epochs.keys())
        plot_structure_evolution(
            str(swd),
            f'{model_dir.name} - Structure Evolution ({len(times)} epochs, '
            f't={times[0]:.2f}-{times[-1]:.2f}d)',
            str(plots_dir / f'{step_name}_structure.png'))
        plot_composition_evolution(
            str(swd),
            f'{model_dir.name} - Composition Evolution ({len(times)} epochs, '
            f't={times[0]:.2f}-{times[-1]:.2f}d)',
            str(plots_dir / f'{step_name}_composition.png'))

    elif step_name == 'homologous':
        # Homologously evolved snapshot
        # Look for homologous .swd file
        swd = model_dir / 'stella' / 'strad' / f'{model_dir.name}.swd'
        if not swd.exists():
            swd = model_dir / 'stella' / f'{model_dir.name}.swd'
        if not swd.exists():
            print(f"  WARNING: .swd not found for homologous plot: {swd}")
            return

        # Find the homologous epoch (349 days)
        epochs = read_swd(str(swd))
        available = sorted(epochs.keys())
        # Look for epoch near 349
        homo_t = min(available, key=lambda t: abs(t - 349))
        if abs(homo_t - 349) > 20:
            print(f"  WARNING: No homologous epoch near 349 days (closest: {homo_t})")
            return
        e = epochs[homo_t]
        # .swd is sorted outer->inner; reverse to inner->outer for plotting
        r = e['r'][::-1]
        rho = e['rho'][::-1]
        v = e['v'][::-1]
        T = e['T'][::-1]
        plot_structure(r, rho, v, T,
                       f'{model_dir.name} - Homologous Evolution (t={homo_t:.1f}d)',
                       str(plots_dir / f'{step_name}_structure.png'))
        comp, _ = read_swd_composition(str(swd), homo_t)
        if comp is not None:
            comp = {k: arr[::-1] for k, arr in comp.items()}
            plot_composition(r, comp,
                           f'{model_dir.name} - Composition (Homologous, t={homo_t:.1f}d)',
                           str(plots_dir / f'{step_name}_composition.png'))
            plot_composition_vs_mass(r, rho, comp,
                                   f'{model_dir.name} - Mass Fractions (Homologous, t={homo_t:.1f}d)',
                                   str(plots_dir / f'{step_name}_mass_fractions.png'),
                                   shaded_mass_msun=0)

    elif step_name == 'cmfgen_input':
        # CMFGEN SN_HYDRO_DATA
        # Look in cmfgen/<day>days/SN_HYDRO_DATA
        for day_dir in sorted((model_dir / 'cmfgen').glob('*days')):
            sn_hydro = day_dir / 'SN_HYDRO_DATA'
            if not sn_hydro.exists():
                continue
            day_label = day_dir.name
            # Read SN_HYDRO_DATA (simple format: r, v, rho, T, ...)
            try:
                data = np.loadtxt(str(sn_hydro), comments='!')
                r = data[:, 0]
                v = data[:, 1] if data.shape[1] > 1 else np.zeros_like(r)
                rho = data[:, 2] if data.shape[1] > 2 else np.zeros_like(r)
                T = data[:, 3] if data.shape[1] > 3 else np.zeros_like(r)
                plot_structure(r, rho, v, T,
                               f'{model_dir.name} - CMFGEN Input ({day_label})',
                               str(plots_dir / f'{step_name}_structure_{day_label}.png'))
            except Exception as ex:
                print(f"  WARNING: Could not read SN_HYDRO_DATA: {ex}")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                 formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('model_dir', help='Model directory (e.g. SN_A or SN_B)')
    parser.add_argument('--step', default='all',
                        help='Plot step: presn, stella_input, stella_output, '
                             'multiband, stella_evolution, homologous, '
                             'cmfgen_input, sn_b_input, or all (default: all)')
    parser.add_argument('--output', default=None,
                        help='Output plots directory (default: <model_dir>/plots)')
    args = parser.parse_args()

    model_dir = Path(args.model_dir)
    plots_dir = Path(args.output) if args.output else model_dir / 'plots'
    plots_dir.mkdir(parents=True, exist_ok=True)

    if args.step == 'all':
        for step in ['presn', 'stella_input', 'stella_output', 'multiband',
                     'stella_evolution', 'homologous', 'cmfgen_input', 'sn_b_input']:
            try:
                plot_step(model_dir, step, plots_dir)
            except Exception as ex:
                print(f"  WARNING: Step '{step}' failed: {ex}")
    else:
        plot_step(model_dir, args.step, plots_dir)

    print(f"\nPlots saved to: {plots_dir}")


if __name__ == '__main__':
    main()
