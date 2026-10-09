#!/usr/bin/env python3
"""Plot SuperLite output files (spectrum, structure, opacity, energy).

Usage:
    python plot_superlite_output.py <superlite_dir> [--output <plots_dir>]

Reads output.* files from a SuperLite run directory and generates
diagnostic plots saved as PNG files.
"""

import argparse
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from pathlib import Path

# Distance for flux conversion: SN2024pba at 83.58 Mpc
# Model luminosity is converted to flux at this distance to match observed data
PC_CM = 3.085677581491367e18
D_MPC = 83.58
D_CM = D_MPC * 1e6 * PC_CM
FOUR_PI_D2 = 4.0 * np.pi * D_CM**2  # 4*pi*d^2 in cm^2

# Rebin factor for MC spectra: with 18000 transport groups, the group grid
# resolution (~2-4 A/bin in the optical) closely matches the flux grid (2 A/bin).
# A light rebin to ~3 A smooths residual MC noise without losing real features.
SPECTRUM_REBIN_A = 3.0


def rebin_to_grid_edges(wl_edges_in, lum_per_bin, wl_edges_out):
    """Rebin luminosity-per-bin onto new (coarser) grid edges.

    This is the physically correct way to reduce MC noise: the transport
    group grid defines the simulation's actual spectral resolution.  Flux
    bins finer than the transport group are independent MC noise.  By
    summing the luminosity from all flux bins within each transport group
    and dividing by the group width, we get the maximum-S/N spectrum at
    the simulation's native resolution.

    Parameters
    ----------
    wl_edges_in : array
        Edges of the input (fine) bins, in Angstrom.
    lum_per_bin : array
        Luminosity (or flux) per input bin (NOT per Angstrom).
    wl_edges_out : array
        Edges of the target (coarse) bins, in Angstrom.

    Returns
    -------
    wl_centers_out : array
        Centers of the output bins.
    lum_per_bin_out : array
        Luminosity per output bin.
    """
    wl_edges_in = np.asarray(wl_edges_in, dtype=float)
    wl_edges_out = np.asarray(wl_edges_out, dtype=float)
    lum_per_bin = np.asarray(lum_per_bin, dtype=float)
    n_in = len(lum_per_bin)
    n_out = len(wl_edges_out) - 1
    lum_out = np.zeros(n_out)
    # For each input bin, find ALL output bins it overlaps and distribute
    # the luminosity proportionally by overlap width.  This guarantees
    # exact flux conservation even when grids are misaligned.
    for i in range(n_in):
        in_lo = wl_edges_in[i]
        in_hi = wl_edges_in[i + 1]
        in_width = in_hi - in_lo
        if in_width <= 0:
            continue
        # First and last output bin indices that overlap this input bin
        j0 = int(np.searchsorted(wl_edges_out, in_lo, side='right')) - 1
        j1 = int(np.searchsorted(wl_edges_out, in_hi, side='right')) - 1
        j0 = max(0, j0)
        j1 = min(n_out - 1, j1)
        for j in range(j0, j1 + 1):
            overlap = min(in_hi, wl_edges_out[j + 1]) - max(in_lo, wl_edges_out[j])
            if overlap > 0:
                lum_out[j] += lum_per_bin[i] * overlap / in_width
    wl_centers_out = 0.5 * (wl_edges_out[:-1] + wl_edges_out[1:])
    return wl_centers_out, lum_out


def rebin_spectrum(wl_ang, flux_density, target_bin_A=SPECTRUM_REBIN_A,
                   wl_edges=None):
    """Rebin a spectrum to approximately target_bin_A angstroms per bin.

    Groups adjacent bins into chunks whose total width is closest to
    target_bin_A, then averages the flux density within each chunk
    (weighted by bin width when wl_edges is provided, for non-uniform
    grids).  This reduces MC noise from fine flux bins that are finer
    than the transport grid resolution.
    """
    if len(wl_ang) < 2:
        return wl_ang, flux_density
    if wl_edges is not None:
        dwl = np.diff(wl_edges)
        current_bin = dwl.mean()
    else:
        current_bin = wl_ang[1] - wl_ang[0]
    if current_bin >= target_bin_A:
        return wl_ang, flux_density  # Already coarse enough
    n_merge = max(1, int(round(target_bin_A / current_bin)))
    n_out = len(wl_ang) // n_merge
    wl_trim = wl_ang[:n_out * n_merge]
    f_trim = flux_density[:n_out * n_merge]
    wl_out = wl_trim.reshape(n_out, n_merge).mean(axis=1)
    if wl_edges is not None:
        # Width-weighted average for non-uniform bins
        dwl_trim = dwl[:n_out * n_merge]
        weights = dwl_trim.reshape(n_out, n_merge)
        f_groups = f_trim.reshape(n_out, n_merge)
        f_out = (f_groups * weights).sum(axis=1) / weights.sum(axis=1)
    else:
        f_out = f_trim.reshape(n_out, n_merge).mean(axis=1)
    return wl_out, f_out


def read_flux_grid(filepath):
    """Read output.flx_grid -> (wavelength centers, bin edges) in cm.

    Returns (bin-center wavelengths (ng values), bin edges (ng+1 values)).
    """
    with open(filepath) as f:
        hdr = f.readline().split()
        ng = int(hdr[1])
        wl_edges = np.fromstring(f.readline(), sep=' ')
    # wl array has ng+1 bin edges; compute ng bin centers
    if len(wl_edges) == ng + 1:
        wl = 0.5 * (wl_edges[:-1] + wl_edges[1:])
    else:
        wl = wl_edges[:ng]
        wl_edges = None
    return wl, wl_edges


def read_flux_luminos(filepath):
    """Read output.flx_luminos -> flux array (erg/s/Hz or similar).

    SuperLite appends one row per NLTE iteration (in_niter), so the file
    may contain multiple rows of ng values each. The LAST row is the
    converged spectrum; earlier rows are intermediate iterations. When
    the file is 2D, return the last row.
    """
    flux = np.loadtxt(filepath)
    if flux.ndim == 0:
        flux = np.array([flux])
    elif flux.ndim == 2:
        # Multiple iterations appended: keep the final (converged) row.
        flux = flux[-1]
    return flux


def read_grd_grid(filepath):
    """Read output.grd_grid -> cell-center radius array (cm).

    Returns ncell cell-center radii to match temperature/opacity arrays.
    """
    with open(filepath) as f:
        f.readline()  # igeom
        dims = f.readline().split()
        ndim = int(dims[1])
        f.readline()  # rnk
        r_edges = np.fromstring(f.readline(), sep=' ')
    # r array has ndim+1 edges; compute ndim cell centers
    if len(r_edges) == ndim + 1:
        r = 0.5 * (r_edges[:-1] + r_edges[1:])
    else:
        r = r_edges[:ndim]
    return r


def read_1d_array(filepath):
    """Read a 1D output file (one value per cell, space-separated).

    SuperLite appends one row per iteration for grid outputs (e.g.
    output.grd_capgrey), so the file may be 2D (niter, ncell). The last
    row is the most converged; return it when the file is 2D.
    """
    data = np.loadtxt(filepath)
    if data.ndim == 0:
        data = np.array([data])
    elif data.ndim == 2:
        data = data[-1]
    return data


def read_input_str(filepath):
    """Read input.str -> structure arrays."""
    with open(filepath) as f:
        data = []
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            vals = np.fromstring(line, sep=' ')
            if len(vals) > 0:
                data.append(vals)
    data = np.array(data)
    r_left = data[:, 0]
    r_right = data[:, 1]
    v_left = data[:, 2]
    v_right = data[:, 3]
    mass = data[:, 4]
    t_rad = data[:, 5]
    temp = data[:, 6]
    abund = data[:, 7:]  # H, He, C, O, Si, Ca, Ni, Co, Fe
    r_mid = (r_left + r_right) / 2
    v_mid = (v_left + v_right) / 2
    return r_mid, v_mid, mass, temp, t_rad, abund, r_left, r_right


def _trim_edge_bins(flux, n_trim=1):
    """Return a mask that drops the first and last n_trim bins.

    SuperLite's first and last frequency bins are boundary cells that
    often show extreme values (orders of magnitude above neighbours)
    due to edge effects in the Monte Carlo binning.  Trimming them
    gives a cleaner spectrum for plotting and comparison.
    """
    mask = np.ones(len(flux), dtype=bool)
    mask[:n_trim] = False
    mask[-n_trim:] = False
    return mask


def prepare_spectrum(wl, flux, wl_edges, group_wl_edges=None,
                     wl_range=None):
    """Trim edge bins, rebin, and convert to flux density.

    Parameters
    ----------
    wl : array
        Flux-grid bin centers (cm).
    flux : array
        Luminosity per flux-grid bin (erg/s).
    wl_edges : array or None
        Flux-grid bin edges (cm).  If None, assume uniform spacing.
    group_wl_edges : array or None
        Transport group-grid edges (cm).  If given, rebin onto this grid
        (the simulation's native spectral resolution) for maximum S/N.
    wl_range : (lo, hi) or None
        If given, restrict to [lo, hi] Angstrom before rebinning.

    Returns
    -------
    wl_rb : array  -- rebinned wavelengths (Angstrom)
    fd_rb : array  -- rebinned flux density (erg/s/cm^2/A)
    """
    wl_ang = wl * 1e8  # cm -> Angstrom
    edge_mask = _trim_edge_bins(flux, n_trim=1)
    wl_ang = wl_ang[edge_mask]
    flux = flux[edge_mask]

    if wl_edges is not None:
        wl_edges_ang = wl_edges * 1e8
        # Edges for the trimmed bins (drop first and last edge)
        trimmed_edges = wl_edges_ang[1:-1]
    else:
        trimmed_edges = None

    # Optional wavelength range cut
    if wl_range is not None:
        lo, hi = wl_range
        rmask = (wl_ang >= lo) & (wl_ang <= hi)
        wl_ang = wl_ang[rmask]
        flux = flux[rmask]
        if trimmed_edges is not None:
            # Align edges to the selected bins: bin i has edges
            # [trimmed_edges[i], trimmed_edges[i+1]], so for selected
            # bins s0..s1 the edges are trimmed_edges[s0:s1+2].
            selected = np.where(rmask)[0]
            if len(selected) > 0:
                trimmed_edges = trimmed_edges[selected[0]:selected[-1] + 2]
            else:
                trimmed_edges = trimmed_edges[:1]

    if group_wl_edges is not None and trimmed_edges is not None:
        # Rebin onto transport group grid (native spectral resolution)
        grp_edges = group_wl_edges * 1e8  # cm -> Angstrom
        # Restrict group edges to the wavelength range of the data
        wl_min, wl_max = wl_ang.min(), wl_ang.max()
        i0 = max(0, np.searchsorted(grp_edges, wl_min, side='right') - 1)
        i1 = min(len(grp_edges), np.searchsorted(grp_edges, wl_max, side='right'))
        grp_edges = grp_edges[i0:i1 + 1]
        if len(grp_edges) >= 2:
            wl_rb, lum_rb = rebin_to_grid_edges(trimmed_edges, flux, grp_edges)
            dwl_rb = np.diff(grp_edges)
            fd_rb = lum_rb / (FOUR_PI_D2 * dwl_rb)
            return wl_rb, fd_rb

    # Fallback: fixed-target rebin on flux density
    if trimmed_edges is not None:
        dwl = np.diff(trimmed_edges)
        if len(dwl) == len(flux):
            flux_density = flux / (FOUR_PI_D2 * dwl)
        else:
            flux_density = flux / (FOUR_PI_D2 * np.full(len(flux), np.mean(dwl)))
    else:
        flux_density = flux / (FOUR_PI_D2 * np.full(len(flux), wl_ang[1] - wl_ang[0]))
    wl_rb, fd_rb = rebin_spectrum(wl_ang, flux_density,
                                  wl_edges=trimmed_edges)
    return wl_rb, fd_rb


def plot_spectrum(wl, flux, wl_edges, outpath, model_name='',
                  group_wl_edges=None):
    """Plot the SuperLite spectrum (flux at 83.58 Mpc vs wavelength).

    First and last frequency bins are dropped (boundary artefacts).
    If *group_wl_edges* (transport group grid edges, in cm) is given the
    spectrum is rebinned onto the transport group grid — the simulation's
    native spectral resolution — which eliminates MC noise from flux bins
    finer than the transport grid.
    """
    fig, ax = plt.subplots(1, 1, figsize=(12, 6))

    wl_rb, fd_rb = prepare_spectrum(wl, flux, wl_edges, group_wl_edges)

    ax.step(wl_rb, fd_rb, 'b-', linewidth=0.8, label='SuperLite', where='mid')
    ax.set_xlabel('Wavelength ($\\AA$)')
    ax.set_ylabel('Flux (erg s$^{-1}$ cm$^{-2}$ $\\AA^{-1}$)')
    ax.set_title(f'SuperLite Spectrum (d={D_MPC} Mpc) - {model_name}')
    ax.set_xlim(wl_rb.min(), wl_rb.max())
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_yscale('log')
    # Auto-scale y-lims to the visible data within xlims
    visible = fd_rb[(wl_rb >= wl_rb.min()) & (wl_rb <= wl_rb.max())]
    visible = visible[visible > 0]
    if len(visible) > 0:
        ymax = visible.max()
        ymin = visible.min()
        ax.set_ylim(ymin * 0.3, ymax * 3)

    plt.tight_layout()
    plt.savefig(outpath, dpi=150)
    plt.close()
    print(f"  Saved spectrum plot: {outpath}")


def plot_optical_spectrum(wl, flux, wl_edges, outpath, model_name='',
                          group_wl_edges=None):
    """Plot the SuperLite optical spectrum (flux at 83.58 Mpc, 3000-10000 Angstrom).

    The first and last frequency bins are dropped — they are boundary
    cells that show extreme values from MC edge effects.
    """
    fig, ax = plt.subplots(1, 1, figsize=(12, 6))

    wl_rb, fd_rb = prepare_spectrum(wl, flux, wl_edges, group_wl_edges,
                                    wl_range=(3000, 10000))
    ax.step(wl_rb, fd_rb, 'b-', linewidth=1.0, label='SuperLite', where='mid')
    ax.set_xlabel('Wavelength ($\\AA$)')
    ax.set_ylabel('Flux (erg s$^{-1}$ cm$^{-2}$ $\\AA^{-1}$)')
    ax.set_title(f'SuperLite Optical Spectrum (d={D_MPC} Mpc) - {model_name}')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_yscale('log')
    # Auto-scale y-lims to the visible (optical) data
    visible = fd_rb
    visible = visible[visible > 0]
    if len(visible) > 0:
        ymax = visible.max()
        ymin = visible.min()
        ax.set_ylim(ymin * 0.3, ymax * 3)

    plt.tight_layout()
    plt.savefig(outpath, dpi=150)
    plt.close()
    print(f"  Saved optical spectrum plot: {outpath}")


def read_observed_spectrum(obs_file):
    """Read an observed spectrum file (wavelength, flux, error).

    Returns (wavelength_A, flux_cgs, error_cgs) or None if file missing.
    """
    data = np.loadtxt(obs_file)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    wl = data[:, 0]
    flux = data[:, 1]
    err = data[:, 2] if data.shape[1] > 2 else None
    return wl, flux, err


def plot_optical_spectrum_vs_obs(wl, flux, wl_edges, obs_file, outpath,
                                  model_name='', obs_label='SN2024pba obs',
                                  group_wl_edges=None):
    """Plot SuperLite optical spectrum overlaid with SN2024pba observed data.

    Both are in flux density (erg/s/cm^2/A) at 83.58 Mpc.
    First and last SuperLite bins are dropped (boundary artefacts).
    """
    obs = read_observed_spectrum(obs_file)
    if obs is None:
        print(f"  WARNING: observed spectrum not found: {obs_file}")
        return
    obs_wl, obs_flux, obs_err = obs

    fig, ax = plt.subplots(1, 1, figsize=(14, 7))

    # --- SuperLite model (rebinned onto transport group grid) ---
    wl_rb, fd_rb = prepare_spectrum(wl, flux, wl_edges, group_wl_edges,
                                    wl_range=(3000, 10000))
    ax.step(wl_rb, fd_rb, 'b-', linewidth=1.2,
            label=f'SuperLite model ({model_name})', where='mid', alpha=0.9)

    # --- Observed spectrum ---
    obs_mask = (obs_wl > 3000) & (obs_wl < 10000)
    ax.plot(obs_wl[obs_mask], obs_flux[obs_mask], 'k-', linewidth=0.5,
            label=obs_label, alpha=0.7)
    if obs_err is not None:
        ax.fill_between(obs_wl[obs_mask],
                        obs_flux[obs_mask] - obs_err[obs_mask],
                        obs_flux[obs_mask] + obs_err[obs_mask],
                        color='gray', alpha=0.2, step='mid')

    ax.set_xlabel('Wavelength ($\\AA$)', fontsize=12)
    ax.set_ylabel('Flux (erg s$^{-1}$ cm$^{-2}$ $\\AA^{-1}$)', fontsize=12)
    ax.set_title(f'SuperLite vs SN2024pba (d={D_MPC} Mpc) - {model_name}', fontsize=14)
    ax.legend(fontsize=10, loc='upper right')
    ax.grid(True, alpha=0.3)
    ax.set_yscale('log')
    ax.set_xlim(3000, 10000)

    # Auto-scale y from both model and obs
    all_vals = np.concatenate([
        fd_rb[fd_rb > 0],
        obs_flux[obs_mask][obs_flux[obs_mask] > 0],
    ])
    if len(all_vals) > 0:
        ax.set_ylim(all_vals.min() * 0.3, all_vals.max() * 3)

    plt.tight_layout()
    plt.savefig(outpath, dpi=150)
    plt.close()
    print(f"  Saved model-vs-obs spectrum plot: {outpath}")


def plot_structure(r, v, mass, temp, abund, r_left, r_right, outpath, model_name=''):
    """Plot ejecta structure (radius, velocity, temperature, composition)."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Top-left: Velocity vs radius
    ax = axes[0, 0]
    ax.step(r / 1e15, v / 1e5, 'b-', linewidth=1.0, where='mid')
    ax.set_xlabel('Radius (1e15 cm)')
    ax.set_ylabel('Velocity (km/s)')
    ax.set_title('Velocity profile')
    ax.grid(True, alpha=0.3)

    # Top-right: Temperature vs radius
    ax = axes[0, 1]
    ax.step(r / 1e15, temp, 'r-', linewidth=1.0, where='mid')
    ax.set_xlabel('Radius (1e15 cm)')
    ax.set_ylabel('Temperature (K)')
    ax.set_title('Temperature profile')
    ax.grid(True, alpha=0.3)

    # Bottom-left: Composition vs radius
    ax = axes[1, 0]
    labels = ['H', 'He', 'C', 'O', 'Si', 'Ca', 'Ni', 'Co', 'Fe']
    colors = ['blue', 'green', 'orange', 'red', 'purple', 'brown', 'pink', 'gray', 'black']
    for i in range(min(len(labels), abund.shape[1])):
        if abund[:, i].max() > 1e-10:
            ax.step(r / 1e15, abund[:, i], label=labels[i], color=colors[i], linewidth=1.0, where='mid')
    ax.set_xlabel('Radius (1e15 cm)')
    ax.set_ylabel('Mass fraction')
    ax.set_title('Composition profile')
    ax.legend(fontsize=8, ncol=3)
    ax.grid(True, alpha=0.3)
    ax.set_yscale('log')
    ax.set_ylim(1e-4, 1)

    # Bottom-right: Density (computed from mass and volume)
    ax = axes[1, 1]
    vol = (4.0/3.0) * np.pi * (r_right**3 - r_left**3)
    rho = mass / vol
    ax.step(r / 1e15, rho, 'k-', linewidth=1.0, where='mid')
    ax.set_xlabel('Radius (1e15 cm)')
    ax.set_ylabel('Density (g/cm^3)')
    ax.set_title('Density profile')
    ax.grid(True, alpha=0.3)
    ax.set_yscale('log')

    fig.suptitle(f'SuperLite Ejecta Structure - {model_name}', fontsize=14, y=1.01)
    plt.tight_layout()
    plt.savefig(outpath, dpi=150)
    plt.close()
    print(f"  Saved structure plot: {outpath}")


def plot_opacity(r, capgrey, capross, outpath, model_name=''):
    """Plot opacity profiles."""
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    ax.step(r / 1e15, capgrey, 'b-', linewidth=1.0, label='Grey opacity', where='mid')
    if capross is not None:
        ax.step(r / 1e15, capross, 'r--', linewidth=1.0, label='Rosseland opacity', where='mid')
    ax.set_xlabel('Radius (1e15 cm)')
    ax.set_ylabel('Opacity (cm^-1)')
    ax.set_title(f'SuperLite Opacity Profile - {model_name}')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_yscale('log')

    plt.tight_layout()
    plt.savefig(outpath, dpi=150)
    plt.close()
    print(f"  Saved opacity plot: {outpath}")


def plot_energy(tot_energy, outpath, model_name=''):
    """Plot energy budget.

    output.tot_energy is appended once per NLTE iteration; use the last
    (converged) row.
    """
    with open(tot_energy) as f:
        f.readline()  # ncol
        f.readline()  # header
        rows = [np.fromstring(ln, sep=' ') for ln in f if ln.strip()]
    vals = rows[-1]

    labels = ['E_out', 'E_velo', 'S_flux', 'S_thermal']
    colors = ['blue', 'green', 'red', 'orange']

    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    bars = ax.bar(labels, vals, color=colors)
    ax.set_ylabel('Energy (erg)')
    ax.set_title(f'SuperLite Energy Budget - {model_name}')
    ax.grid(True, alpha=0.3, axis='y')
    for bar, val in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height(),
                f'{val:.2e}', ha='center', va='bottom', fontsize=9)

    plt.tight_layout()
    plt.savefig(outpath, dpi=150)
    plt.close()
    print(f"  Saved energy plot: {outpath}")


def read_grp_data(filepath):
    """Read output.grp_* file -> (ncells x ngroups) array.

    Skips 3 header lines: title, column labels, dimensions.
    """
    with open(filepath) as f:
        f.readline()  # title
        f.readline()  # column labels
        dims = f.readline().split()
        ng = int(dims[1])
        ncell = int(dims[2])
        data = np.loadtxt(f)
    return data.reshape(ncell, ng)


# Atomic number -> element symbol mapping (for SDEC plot labels)
ELEMENT_SYMBOLS = {
    1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O',
    9: 'F', 10: 'Ne', 11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P',
    16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca', 21: 'Sc', 22: 'Ti',
    23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu',
    30: 'Zn',
}

# Distinct colors for common SN elements in SDEC plot
ELEMENT_COLORS = {
    'H': '#1f77b4', 'He': '#ff7f0e', 'C': '#2ca02c', 'N': '#98df8a',
    'O': '#d62728', 'Ne': '#9467bd', 'Mg': '#8c564b', 'Si': '#e377c2',
    'S': '#bcbd22', 'Ca': '#17becf', 'Ti': '#aec7e8', 'Cr': '#ffbb78',
    'Fe': '#ff9896', 'Co': '#c5b0d5', 'Ni': '#c49c94',
}


def read_elem_list(filepath):
    """Read output.elem_list -> list of (iz, symbol, total_cap) tuples."""
    elements = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            iz = int(parts[0])
            sym = parts[1]
            total = float(parts[2]) if len(parts) > 2 else 0.0
            elements.append((iz, sym, total))
    return elements


def plot_sdec_spectrum(sl_dir, wl_grp, wl_grp_edges, outpath, model_name='', v=None):
    """TARDIS-SDEC-style spectral decomposition plot for SuperLite.

    Following the TARDIS SDEC approach (Kromer et al. 2013, Fig. 4):
    - Emission (positive, above zero): per-element emissivity contributions,
      Doppler-shifted to observed frame, globally scaled to the MC luminosity.
    - Absorption (negative, below zero): total absorption opacity profile,
      Doppler-shifted to observed frame, showing where photons are removed.
    - MC spectrum overlaid as a dashed blue line (the actual escaping flux).
    - Blackbody photosphere overlaid as a dashed red line (T_inner, R_inner).

    Unlike TARDIS which tracks per-packet last interactions, SuperLite outputs
    per-element emissivity (total thermal emission, not escaping flux). We
    therefore use a SINGLE global scale factor (not per-wavelength normalization)
    to match the total emissivity luminosity to the total MC luminosity. This
    preserves the relative wavelength structure of each element's contribution
    without forcing artificial agreement at every wavelength.

    If v (velocity array in cm/s, one per cell) is provided, each cell's
    emissivity is Doppler-shifted to the observed frame before summing.
    """
    from astropy.modeling.models import BlackBody
    import astropy.units as u

    elem_list_file = sl_dir / 'output.elem_list'
    if not elem_list_file.exists():
        print("  WARNING: No output.elem_list found for SDEC plot")
        return

    elements = read_elem_list(elem_list_file)
    if not elements:
        print("  WARNING: No elements with non-zero opacity found for SDEC plot")
        return

    # --- Read the MC output spectrum (escaping flux) ---
    flx_lum_file = sl_dir / 'output.flx_luminos'
    flx_grid_file = sl_dir / 'output.flx_grid'
    if not (flx_lum_file.exists() and flx_grid_file.exists()):
        print("  WARNING: No output.flx_luminos/flx_grid for SDEC plot")
        return

    with open(flx_grid_file) as f:
        hdr = f.readline().split()
        ng_flx = int(hdr[1])
        wl_flx_edges = np.fromstring(f.readline(), sep=' ')
    with open(flx_lum_file) as f:
        actual_lum = np.fromstring(f.readline(), sep=' ')

    if len(wl_flx_edges) == ng_flx + 1:
        wl_flx_ang = 0.5 * (wl_flx_edges[:-1] + wl_flx_edges[1:]) * 1e8
        dwl_flx_ang = np.diff(wl_flx_edges * 1e8)
    else:
        wl_flx_ang = wl_flx_edges[:ng_flx] * 1e8
        dwl_flx_ang = np.full(ng_flx, wl_flx_ang[1] - wl_flx_ang[0])

    # Trim first/last bins (boundary artefacts)
    edge_mask = _trim_edge_bins(actual_lum, n_trim=1)
    actual_lum = actual_lum[edge_mask]
    wl_flx_ang = wl_flx_ang[edge_mask]
    dwl_flx_ang = dwl_flx_ang[edge_mask]
    actual_spec = actual_lum / (FOUR_PI_D2 * dwl_flx_ang)  # erg/s/cm^2/A
    print(f"  MC output: {len(actual_lum)} bins, "
          f"{wl_flx_ang.min():.0f}-{wl_flx_ang.max():.0f} A")

    # --- Read per-element emissivity on the group grid ---
    wl_grp_ang = wl_grp * 1e8  # cm -> Angstrom
    if wl_grp_edges is not None:
        dwl_grp_ang = np.diff(wl_grp_edges * 1e8)
    else:
        dwl_grp_ang = np.full(len(wl_grp_ang), wl_grp_ang[1] - wl_grp_ang[0])

    spec_data = {}  # symbol -> (emissivity per cell, iz)
    for iz, sym, total in elements:
        emiss_file = sl_dir / f'output.grp_emiss_Z{iz:02d}'
        if not emiss_file.exists():
            continue
        emiss = read_grp_data(emiss_file)  # (ncell, ngroup)
        label = ELEMENT_SYMBOLS.get(iz, sym.strip())
        spec_data[label] = (emiss, iz)

    if not spec_data:
        print("  WARNING: No per-element emissivity files found for SDEC plot")
        return

    # --- Read total absorption opacity for absorption panel ---
    grp_cap_file = sl_dir / 'output.grp_cap'
    cap_data = None
    if grp_cap_file.exists():
        cap_data = read_grp_data(grp_cap_file)  # (ncell, ngroup)

    # --- Read cell radii from input.str to compute volumes ---
    # The grp_emiss_Z## files contain κ_elem × B_λ (emissivity per unit volume
    # per steradian). To get total luminosity, we must multiply by cell volume
    # and 4π. Without the volume factor, small inner cells get the same weight
    # as large outer cells, which distorts the SDEC plot.
    # Also read T_inner and R_inner for the blackbody photosphere curve.
    cell_vol = None
    t_inner = 0.0
    r_inner_cm = 0.0
    input_str = sl_dir / 'input.str'
    if input_str.exists():
        r_left_list = []
        r_right_list = []
        first_data = True
        with open(input_str) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                vals = np.fromstring(line, sep=' ')
                if len(vals) >= 7:
                    r_left_list.append(vals[0])   # x_left (cm)
                    r_right_list.append(vals[1])  # x_right (cm)
                    if first_data:
                        r_inner_cm = vals[0]      # x_left of first zone
                        t_inner = vals[6]          # temperature of first zone
                        first_data = False
        if r_left_list:
            r_left = np.array(r_left_list)
            r_right = np.array(r_right_list)
            # Spherical shell volume: V = 4π/3 * (r_outer^3 - r_inner^3)
            cell_vol = (4.0 * np.pi / 3.0) * (r_right**3 - r_left**3)
            print(f"  Cell volumes: {cell_vol.min():.2e} - {cell_vol.max():.2e} cm^3")
            print(f"  Photosphere: T_inner={t_inner:.1f} K, R_inner={r_inner_cm:.3e} cm")

    # --- Doppler-shift emissivity to observed frame and sum over cells ---
    c_kms = 2.9979e5
    if v is not None:
        v_kms = v / 1e5  # cm/s -> km/s
        print(f"  Doppler-shifting emissivity to observed frame "
              f"(v={v_kms.min():.0f}-{v_kms.max():.0f} km/s)")
    else:
        v_kms = None

    n_elem = len(spec_data)
    n_flx = len(wl_flx_ang)
    emiss_on_flx = np.zeros((n_elem, n_flx))
    labels = []
    colors = []
    for i, (label, (emiss, iz)) in enumerate(sorted(spec_data.items(),
                                                     key=lambda x: x[1][0].sum(),
                                                     reverse=True)):
        if v_kms is not None:
            spectrum = np.zeros(n_flx)
            for ic in range(emiss.shape[0]):
                if emiss[ic].max() <= 0:
                    continue
                # Multiply by cell volume (converts emissivity to luminosity)
                vol_factor = cell_vol[ic] if cell_vol is not None else 1.0
                doppler = 1.0 + v_kms[ic] / c_kms
                wl_obs_cell = wl_grp_ang * doppler
                cell_density = np.maximum(emiss[ic], 0) * vol_factor / dwl_grp_ang
                spectrum += np.interp(wl_flx_ang, wl_obs_cell, cell_density,
                                      left=0, right=0)
        else:
            if cell_vol is not None:
                spectrum = (emiss * cell_vol[:, np.newaxis]).sum(axis=0)
            else:
                spectrum = emiss.sum(axis=0)
            spectrum = np.maximum(spectrum, 0) / dwl_grp_ang
            spectrum = np.interp(wl_flx_ang, wl_grp_ang, spectrum,
                                 left=0, right=0)
        emiss_on_flx[i] = spectrum
        labels.append(label)
        colors.append(ELEMENT_COLORS.get(label, plt.cm.tab20(i / max(n_elem, 1))))

    # --- Global scale factor (NOT per-wavelength normalization) ---
    # Scale the total emissivity to match the total MC luminosity over the
    # optical range. This preserves the relative wavelength structure of
    # each element's contribution without forcing artificial agreement.
    emiss_total = emiss_on_flx.sum(axis=0)
    vis_mask = (wl_flx_ang >= 3000) & (wl_flx_ang <= 10000)
    emiss_int = np.trapezoid(emiss_total[vis_mask], wl_flx_ang[vis_mask])
    spec_int = np.trapezoid(actual_spec[vis_mask], wl_flx_ang[vis_mask])
    if emiss_int > 0:
        global_scale = spec_int / emiss_int
    else:
        global_scale = 1.0
    print(f"  Global emissivity scale factor: {global_scale:.4e}")
    stacked = emiss_on_flx * global_scale

    # --- Doppler-shift absorption opacity for absorption panel ---
    # The grp_cap file contains total absorption opacity (cm^-1) per cell.
    # We weight by cell volume to get total absorption cross-section, then
    # Doppler-shift to observed frame. This shows where absorption happens
    # in wavelength space (proxy for TARDIS's absorbed luminosity).
    abs_profile = np.zeros(n_flx)
    if cap_data is not None and v_kms is not None:
        for ic in range(cap_data.shape[0]):
            if cap_data[ic].max() <= 0:
                continue
            vol_factor = cell_vol[ic] if cell_vol is not None else 1.0
            doppler = 1.0 + v_kms[ic] / c_kms
            wl_obs_cell = wl_grp_ang * doppler
            cell_cap = np.maximum(cap_data[ic], 0) * vol_factor / dwl_grp_ang
            abs_profile += np.interp(wl_flx_ang, wl_obs_cell, cell_cap,
                                      left=0, right=0)
        # Scale absorption to a reasonable display amplitude
        abs_vis = abs_profile[vis_mask]
        if abs_vis.max() > 0:
            abs_scale = actual_spec[vis_mask].max() * 0.3 / abs_vis.max()
            abs_profile *= abs_scale

    # --- Compute blackbody photosphere ---
    # Uses T_inner and R_inner already read from input.str above.
    bb_lam = None
    if t_inner > 0 and r_inner_cm > 0:
        bb = BlackBody(t_inner * u.K)
        # bb() returns spectral radiance per Hz per sr; convert to erg/(AA s)
        bb_nu = bb(wl_flx_ang * u.AA)  # erg/(Hz s cm2 sr)
        # Convert to per-Angstrom: multiply by c/lambda^2 (dnu/dlam)
        c_cgs = 2.99792458e10  # cm/s
        dnu_dlam = c_cgs / (wl_flx_ang * 1e-8)**2  # Hz/Angstrom (c/lam^2 with lam in cm)
        bb_lam_vals = (bb_nu.value * dnu_dlam * 4 * np.pi * r_inner_cm**2
                       ) / FOUR_PI_D2  # erg/(AA s)
        bb_lam = bb_lam_vals

    # --- Rebin to reduce MC noise ---
    n_merge = max(1, int(round(SPECTRUM_REBIN_A / (wl_flx_ang[1] - wl_flx_ang[0])))
                  if len(wl_flx_ang) > 1 else 1)
    if n_merge > 1:
        n_out = len(wl_flx_ang) // n_merge
        wl_rb = wl_flx_ang[:n_out * n_merge].reshape(n_out, n_merge).mean(axis=1)
        actual_spec_rb = actual_spec[:n_out * n_merge].reshape(n_out, n_merge).mean(axis=1)
        stacked_rb = stacked[:, :n_out * n_merge].reshape(n_elem, n_out, n_merge).mean(axis=2)
        abs_rb = abs_profile[:n_out * n_merge].reshape(n_out, n_merge).mean(axis=1)
        if bb_lam is not None:
            bb_rb = bb_lam[:n_out * n_merge].reshape(n_out, n_merge).mean(axis=1)
        else:
            bb_rb = None
    else:
        wl_rb, actual_spec_rb, stacked_rb = wl_flx_ang, actual_spec, stacked
        abs_rb = abs_profile
        bb_rb = bb_lam
    n_rb = len(wl_rb)

    # --- Plot (TARDIS-style: single panel, emission up, absorption down) ---
    fig, ax = plt.subplots(1, 1, figsize=(14, 8))

    # Emission: stack upward from zero
    cumulative = np.cumsum(stacked_rb, axis=0)
    lower = np.vstack([np.zeros(n_rb), cumulative[:-1]])
    for i in range(n_elem):
        ax.fill_between(wl_rb, lower[i], cumulative[i], step='mid',
                        color=colors[i], alpha=0.7, label=labels[i],
                        linewidth=0)

    # Absorption: stack downward from zero (TARDIS style)
    if abs_rb.max() > 0:
        ax.fill_between(wl_rb, 0, -abs_rb, step='mid',
                        color='silver', alpha=0.5, label='Absorption',
                        linewidth=0)

    # MC spectrum as dashed blue line (TARDIS style)
    ax.plot(wl_rb, actual_spec_rb, '--b', linewidth=1.5,
            label='MC Spectrum', zorder=10)

    # Blackbody photosphere as dashed red line
    if bb_rb is not None:
        ax.plot(wl_rb, bb_rb, '--r', linewidth=1.0,
                label='Blackbody Photosphere', zorder=9)

    ax.axhline(0, color='k', linewidth=0.5)
    ax.set_xlabel('Wavelength ($\\AA$)', fontsize=12)
    ax.set_ylabel('Flux (erg s$^{-1}$ cm$^{-2}$ $\\AA^{-1}$)', fontsize=12)
    ax.set_title(f'SuperLite SDEC (d={D_MPC} Mpc) - {model_name}',
                 fontsize=14)
    ax.legend(fontsize=9, ncol=3, loc='upper right')
    ax.grid(True, alpha=0.2)
    ax.set_xlim(3000, 10000)

    # Y-limits: based on MC spectrum max and absorption depth
    vis = (wl_rb >= 3000) & (wl_rb <= 10000)
    ymax = max(actual_spec_rb[vis].max(), cumulative[:, vis].max()) * 1.3
    ymin = -abs_rb[vis].max() * 1.3 if abs_rb.max() > 0 else 0
    ax.set_ylim(ymin, ymax)

    plt.tight_layout()
    plt.savefig(outpath, dpi=150)
    plt.close()
    print(f"  Saved SDEC spectrum: {outpath}")


def read_grp_wavelengths(filepath):
    """Read output.grp_grid -> (wavelength centers, bin edges) in cm."""
    with open(filepath) as f:
        ng = int(f.readline().split()[1])
        wl_edges = np.fromstring(f.readline(), sep=' ')
    if len(wl_edges) == ng + 1:
        wl = 0.5 * (wl_edges[:-1] + wl_edges[1:])
    else:
        wl = wl_edges[:ng]
        wl_edges = None
    return wl, wl_edges


def plot_emissivity_heatmap(wl_grp, emiss, v, outpath, model_name=''):
    """2D heatmap: emissivity vs wavelength and velocity."""
    wl_ang = wl_grp * 1e8
    v_kms = v / 1e5

    fig, ax = plt.subplots(1, 1, figsize=(14, 8))

    # Mask zeros for log scale
    em = emiss.copy()
    em[em <= 0] = np.nan

    pcm = ax.pcolormesh(wl_ang, v_kms, em, shading='auto',
                        norm=LogNorm(vmin=np.nanmin(em), vmax=np.nanmax(em)),
                        cmap='hot')
    ax.set_xlabel('Wavelength (Angstrom)')
    ax.set_ylabel('Velocity (km/s)')
    ax.set_title(f'SuperLite Emissivity Distribution - {model_name}')
    cbar = fig.colorbar(pcm, ax=ax, pad=0.02)
    cbar.set_label('Emissivity (erg/s/cm$^3$)')

    plt.tight_layout()
    plt.savefig(outpath, dpi=150)
    plt.close()
    print(f"  Saved emissivity heatmap: {outpath}")


def plot_opacity_heatmap(wl_grp, cap, v, outpath, model_name=''):
    """2D heatmap: absorption opacity vs wavelength and velocity."""
    wl_ang = wl_grp * 1e8
    v_kms = v / 1e5

    fig, ax = plt.subplots(1, 1, figsize=(14, 8))

    op = cap.copy()
    op[op <= 0] = np.nan

    pcm = ax.pcolormesh(wl_ang, v_kms, op, shading='auto',
                        norm=LogNorm(vmin=np.nanmin(op), vmax=np.nanmax(op)),
                        cmap='RdYlBu_r')
    ax.set_xlabel('Wavelength (Angstrom)')
    ax.set_ylabel('Velocity (km/s)')
    ax.set_title(f'SuperLite Absorption Opacity Distribution - {model_name}')
    cbar = fig.colorbar(pcm, ax=ax, pad=0.02)
    cbar.set_label('Opacity (cm$^{-1}$)')

    plt.tight_layout()
    plt.savefig(outpath, dpi=150)
    plt.close()
    print(f"  Saved opacity heatmap: {outpath}")


def plot_opacity_observed_frame(wl_grp, cap, v, wl_flx_ang, actual_spec,
                                 outpath, model_name=''):
    """1D plot: total line opacity projected to observed-frame wavelengths,
    overlaid with the MC spectrum.

    In homologous expansion, a line at rest wavelength λ_rest in a cell moving
    at velocity v produces absorption at observed wavelength
        λ_obs = λ_rest × (1 + v/c).
    This function Doppler-shifts each cell's opacity spectrum to the observed
    frame and sums over all cells, producing a 1D "observed opacity" profile
    that can be directly compared with the spectrum.
    """
    c_kms = 2.9979e5  # km/s
    wl_rest_ang = wl_grp * 1e8  # rest-frame Angstrom
    v_kms = v / 1e5             # km/s (input.str velocities are in cm/s)

    # Observed wavelength grid (from flux grid)
    n_flx = len(wl_flx_ang)
    obs_opacity = np.zeros(n_flx)

    # For each cell, shift its opacity spectrum to observed frame and
    # interpolate onto the flux grid.
    for ic in range(len(v_kms)):
        doppler = 1.0 + v_kms[ic] / c_kms
        wl_obs_cell = wl_rest_ang * doppler
        cell_cap = cap[ic, :]
        # Only interpolate where opacity > 0 (skip zeros for speed)
        if cell_cap.max() <= 0:
            continue
        obs_opacity += np.interp(wl_flx_ang, wl_obs_cell, cell_cap,
                                 left=0, right=0)

    fig, ax1 = plt.subplots(1, 1, figsize=(14, 8))

    # Normalize opacity for display (log scale)
    vis = (wl_flx_ang >= 3000) & (wl_flx_ang <= 10000)
    obs_vis = obs_opacity[vis]
    spec_vis = actual_spec[vis]
    pos_op = obs_vis[obs_vis > 0]
    pos_sp = spec_vis[spec_vis > 0]
    if len(pos_op) == 0 or len(pos_sp) == 0:
        print("  WARNING: No positive opacity or spectrum in visible range")
        plt.close()
        return

    # Rebin spectrum and opacity to reduce MC noise
    wl_vis = wl_flx_ang[vis]
    wl_rb, spec_rb = rebin_spectrum(wl_vis, spec_vis)
    _, obs_rb = rebin_spectrum(wl_vis, obs_vis)

    # Plot spectrum on left y-axis
    color1 = 'steelblue'
    ax1.step(wl_rb, spec_rb, color=color1, linewidth=0.8,
             label='MC spectrum', where='mid', zorder=5)
    ax1.set_xlabel('Wavelength ($\\AA$)', fontsize=12)
    ax1.set_ylabel('Flux (erg s$^{-1}$ cm$^{-2}$ $\\AA^{-1}$)', fontsize=12,
                   color=color1)
    ax1.tick_params(axis='y', labelcolor=color1)
    ax1.set_yscale('log')
    ymax_sp = pos_sp.max()
    ax1.set_ylim(ymax_sp * 1e-3, ymax_sp * 3)

    # Plot observed-frame opacity on right y-axis (inverted)
    ax2 = ax1.twinx()
    color2 = 'crimson'
    ax2.step(wl_rb, obs_rb, color=color2, linewidth=0.8,
             alpha=0.7, label='Line opacity (obs. frame)', where='mid')
    ax2.set_ylabel('Total opacity (cm$^{-1}$)', fontsize=12, color=color2)
    ax2.tick_params(axis='y', labelcolor=color2)
    ax2.set_yscale('log')
    ymax_op = pos_op.max()
    ax2.set_ylim(ymax_op * 1e-3, ymax_op * 10)
    # Invert opacity axis so peaks point down (toward absorption)
    ax2.invert_yaxis()

    ax1.set_xlim(3000, 10000)
    ax1.set_title(f'Spectrum vs Observed-Frame Opacity - {model_name}',
                  fontsize=14)

    # Combined legend
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, fontsize=10, loc='lower right')

    ax1.grid(True, alpha=0.2)
    plt.tight_layout()
    plt.savefig(outpath, dpi=150)
    plt.close()
    print(f"  Saved observed-frame opacity overlay: {outpath}")


def plot_opacity_heatmap_observed(wl_grp, cap, v, outpath, model_name=''):
    """2D heatmap: absorption opacity vs OBSERVED wavelength and velocity.

    Unlike plot_opacity_heatmap (which shows rest-frame wavelengths), this
    shifts each cell's opacity to the observed frame:
        λ_obs = λ_rest × (1 + v/c)
    so that the opacity features align with the observed spectrum.
    """
    c_kms = 2.9979e5
    wl_rest_ang = wl_grp * 1e8  # Angstrom
    v_kms = v / 1e5

    # Build observed-frame wavelength grid (use the flux grid range)
    wl_obs_min, wl_obs_max = 2500.0, 11000.0
    n_obs = 2000
    wl_obs = np.linspace(wl_obs_min, wl_obs_max, n_obs)

    ncell = len(v_kms)
    obs_cap = np.zeros((ncell, n_obs))
    for ic in range(ncell):
        doppler = 1.0 + v_kms[ic] / c_kms
        wl_obs_cell = wl_rest_ang * doppler
        cell_cap = cap[ic, :]
        if cell_cap.max() <= 0:
            continue
        obs_cap[ic, :] = np.interp(wl_obs, wl_obs_cell, cell_cap,
                                    left=0, right=0)

    fig, ax = plt.subplots(1, 1, figsize=(14, 8))
    op = obs_cap.copy()
    op[op <= 0] = np.nan

    pcm = ax.pcolormesh(wl_obs, v_kms, op, shading='auto',
                        norm=LogNorm(vmin=np.nanmin(op), vmax=np.nanmax(op)),
                        cmap='RdYlBu_r')
    ax.set_xlabel('Observed Wavelength ($\\AA$)', fontsize=12)
    ax.set_ylabel('Velocity (km/s)', fontsize=12)
    ax.set_title(f'Observed-Frame Opacity Distribution - {model_name}',
                 fontsize=14)
    ax.set_xlim(3000, 10000)
    cbar = fig.colorbar(pcm, ax=ax, pad=0.02)
    cbar.set_label('Opacity (cm$^{-1}$)')

    plt.tight_layout()
    plt.savefig(outpath, dpi=150)
    plt.close()
    print(f"  Saved observed-frame opacity heatmap: {outpath}")


def plot_stacked_spectrum(wl_grp, src_lum, v, outpath, model_name='', n_bins=8):
    """Stacked spectrum showing source luminosity from each velocity bin."""
    wl_ang = wl_grp * 1e8
    v_kms = v / 1e5

    # Bin cells into velocity shells
    v_min, v_max = v_kms.min(), v_kms.max()
    v_edges = np.linspace(v_min, v_max + 1, n_bins + 1)
    v_centers = 0.5 * (v_edges[:-1] + v_edges[1:])

    # Sum source luminosity within each velocity bin
    stacked = np.zeros((n_bins, len(wl_grp)))
    for i in range(n_bins):
        mask = (v_kms >= v_edges[i]) & (v_kms < v_edges[i + 1])
        if mask.sum() > 0:
            stacked[i] = src_lum[mask].sum(axis=0)

    # Colors: inner (hot) to outer (cool)
    colors = plt.cm.RdYlBu_r(np.linspace(0.1, 0.9, n_bins))

    fig, ax = plt.subplots(1, 1, figsize=(14, 8))

    # Only plot optical range if data is wide
    opt_mask = (wl_ang > 2000) & (wl_ang < 12000)
    if opt_mask.sum() > 0:
        wl_plot = wl_ang[opt_mask]
        stacked_plot = stacked[:, opt_mask]
    else:
        wl_plot = wl_ang
        stacked_plot = stacked

    # Stepped stacked spectrum using fill_between with step='mid'
    v_cumulative = np.cumsum(stacked_plot, axis=0)
    v_lower = np.vstack([np.zeros(len(wl_plot)), v_cumulative[:-1]])
    v_labels = [f'{v_centers[i]:.0f} km/s' for i in range(n_bins)]
    for i in range(n_bins):
        ax.fill_between(wl_plot, v_lower[i], v_cumulative[i], step='mid',
                        color=colors[i], alpha=0.8, label=v_labels[i])
    ax.set_xlabel('Wavelength (Angstrom)')
    ax.set_ylabel('Source Luminosity (erg/s)')
    ax.set_title(f'SuperLite Spectral Decomposition by Velocity - {model_name}')
    ax.legend(fontsize=8, title='Velocity bin', ncol=2, loc='upper right')
    ax.grid(True, alpha=0.3)
    ax.set_yscale('log')
    ax.set_ylim(bottom=max(stacked_plot[stacked_plot > 0].min() * 0.1, 1e30))

    plt.tight_layout()
    plt.savefig(outpath, dpi=150)
def main():
    parser = argparse.ArgumentParser(description='Plot SuperLite output files.')
    parser.add_argument('superlite_dir', help='SuperLite run directory with output.* files')
    parser.add_argument('--output', '-o', default=None, help='Output plots directory (default: <superlite_dir>/plots)')
    parser.add_argument('--model-name', default='', help='Model name for plot titles')
    parser.add_argument('--epoch', default='', help='Epoch label for filenames (e.g. 30d)')
    parser.add_argument('--obs-file', default=None,
                        help='Observed spectrum file (wl, flux, error) to overlay on optical plot')
    args = parser.parse_args()

    sl_dir = Path(args.superlite_dir)
    if not sl_dir.exists():
        print(f"ERROR: Directory not found: {sl_dir}")
        sys.exit(1)

    plots_dir = Path(args.output) if args.output else sl_dir / 'plots'
    plots_dir.mkdir(parents=True, exist_ok=True)

    model_name = args.model_name or sl_dir.parent.parent.name

    # Epoch tag for filenames: derive from --epoch, or from directory name, or from model_name
    epoch = args.epoch
    if not epoch:
        # Try to extract from directory path (e.g. .../superlite/30.0days -> 30d)
        parts = sl_dir.parts
        for p in parts:
            if p.endswith('days'):
                try:
                    day_val = float(p.replace('days', ''))
                    epoch = f"{day_val:g}d"
                except ValueError:
                    pass
                break
    if not epoch:
        # Try to extract from model_name (e.g. "X_1_fix 30d" -> "30d")
        for p in model_name.split():
            if p.endswith('d') and p[:-1].replace('.', '').isdigit():
                epoch = p
                break
    if not epoch:
        epoch = 'snapshot'

    # Filename prefix: superlite_{epoch}_{type}.png
    prefix = f'superlite_{epoch}_'

    print(f"Plotting SuperLite output from: {sl_dir}")
    print(f"Saving plots to: {plots_dir}")
    print(f"Model name: {model_name}")
    print(f"Epoch tag: {epoch}")

    # Spectrum
    flx_grid = sl_dir / 'output.flx_grid'
    flx_lum = sl_dir / 'output.flx_luminos'
    # Transport group grid — defines the simulation's native spectral
    # resolution.  Used to rebin the (finer) flux grid and eliminate MC
    # noise from sub-transport-group bins.
    grp_grid = sl_dir / 'output.grp_grid'
    group_wl_edges = None
    if grp_grid.exists():
        _, group_wl_edges = read_grp_wavelengths(grp_grid)
    if flx_grid.exists() and flx_lum.exists():
        wl, wl_edges = read_flux_grid(flx_grid)
        flux = read_flux_luminos(flx_lum)
        plot_spectrum(wl, flux, wl_edges, plots_dir / f'{prefix}spectrum.png',
                      model_name, group_wl_edges=group_wl_edges)
        plot_optical_spectrum(wl, flux, wl_edges,
                              plots_dir / f'{prefix}optical_spectrum.png',
                              model_name, group_wl_edges=group_wl_edges)
        # Model-vs-observed comparison
        obs_file = args.obs_file
        if obs_file is None:
            # Default: look for SN2024pba observed spectrum near the model epoch
            default_obs = Path('/viper/ptmp1/sbrennan/SN2024pba/observational_data/'
                               'SN2024pba_60875.145_PHASE25.1_Nordic-Optical-Telescope_ALFOSC.txt')
            if default_obs.exists():
                obs_file = str(default_obs)
        if obs_file and Path(obs_file).exists():
            obs_label = Path(obs_file).stem.replace('_', ' ')
            plot_optical_spectrum_vs_obs(
                wl, flux, wl_edges, obs_file,
                plots_dir / f'{prefix}optical_spectrum_vs_obs.png',
                model_name=model_name, obs_label=obs_label,
                group_wl_edges=group_wl_edges)
        else:
            print("  INFO: No --obs-file given and default SN2024pba spectrum not found; "
                  "skipping model-vs-obs comparison plot")
    else:
        print("  WARNING: No spectrum files found (output.flx_grid/output.flx_luminos)")

    # Structure (from input.str)
    input_str = sl_dir / 'input.str'
    if input_str.exists():
        r, v, mass, temp, t_rad, abund, r_left, r_right = read_input_str(input_str)
        plot_structure(r, v, mass, temp, abund, r_left, r_right, plots_dir / f'{prefix}structure.png', model_name)
    else:
        print("  WARNING: No input.str found")

    # Opacity (from output files)
    grd_grid = sl_dir / 'output.grd_grid'
    capgrey = sl_dir / 'output.grd_capgrey'
    capross = sl_dir / 'output.grd_capross'
    if grd_grid.exists() and capgrey.exists():
        r = read_grd_grid(grd_grid)
        cap = read_1d_array(capgrey)
        cap_r = read_1d_array(capross) if capross.exists() else None
        plot_opacity(r, cap, cap_r, plots_dir / f'{prefix}opacity.png', model_name)
    else:
        print("  WARNING: No opacity files found")

    # Energy budget
    tot_energy = sl_dir / 'output.tot_energy'
    if tot_energy.exists():
        plot_energy(tot_energy, plots_dir / f'{prefix}energy.png', model_name)
    else:
        print("  WARNING: No output.tot_energy found")

    # Spectral decomposition (emissivity/opacity heatmaps + SDEC by species)
    grp_grid = sl_dir / 'output.grp_grid'
    grp_emiss = sl_dir / 'output.grp_emiss'
    grp_cap = sl_dir / 'output.grp_cap'
    if grp_grid.exists() and input_str.exists():
        wl_grp, wl_grp_edges = read_grp_wavelengths(grp_grid)
        r, v, mass, temp, t_rad, abund, r_left, r_right = read_input_str(input_str)
        if grp_emiss.exists():
            emiss = read_grp_data(grp_emiss)
            plot_emissivity_heatmap(wl_grp, emiss, v, plots_dir / f'{prefix}emissivity_heatmap.png', model_name)
        else:
            print("  WARNING: No output.grp_emiss found (need in_io_opacdump='one' or 'all')")
        if grp_cap.exists():
            cap = read_grp_data(grp_cap)
            plot_opacity_heatmap(wl_grp, cap, v, plots_dir / f'{prefix}opacity_heatmap.png', model_name)
            # Also generate observed-frame versions that align with the spectrum
            plot_opacity_heatmap_observed(wl_grp, cap, v,
                plots_dir / f'{prefix}opacity_heatmap_observed.png', model_name)
            # 1D overlay: spectrum vs observed-frame opacity
            flx_grid_file = sl_dir / 'output.flx_grid'
            flx_lum_file = sl_dir / 'output.flx_luminos'
            if flx_grid_file.exists() and flx_lum_file.exists():
                with open(flx_grid_file) as f:
                    hdr = f.readline().split()
                    ng_flx = int(hdr[1])
                    wl_flx_edges = np.fromstring(f.readline(), sep=' ')
                with open(flx_lum_file) as f:
                    # SuperLite appends one row per NLTE iteration; the
                    # last row is the converged spectrum.
                    lines = [np.fromstring(ln, sep=' ') for ln in f if ln.strip()]
                    actual_lum = lines[-1] if lines else np.array([])
                if len(wl_flx_edges) == ng_flx + 1:
                    wl_flx_ang = 0.5 * (wl_flx_edges[:-1] + wl_flx_edges[1:]) * 1e8
                    dwl_flx_ang = np.diff(wl_flx_edges * 1e8)
                else:
                    wl_flx_ang = wl_flx_edges[:ng_flx] * 1e8
                    dwl_flx_ang = np.full(ng_flx, wl_flx_ang[1] - wl_flx_ang[0])
                edge_mask = _trim_edge_bins(actual_lum, n_trim=1)
                wl_flx_ang = wl_flx_ang[edge_mask]
                actual_lum = actual_lum[edge_mask]
                dwl_flx_ang = dwl_flx_ang[edge_mask]
                actual_spec = actual_lum / (FOUR_PI_D2 * dwl_flx_ang)
                plot_opacity_observed_frame(wl_grp, cap, v, wl_flx_ang,
                    actual_spec,
                    plots_dir / f'{prefix}opacity_vs_spectrum.png',
                    model_name)
        else:
            print("  WARNING: No output.grp_cap found (need in_io_opacdump='one' or 'all')")
        # TARDIS-SDEC-style per-species decomposition (flux at 10 pc)
        plot_sdec_spectrum(sl_dir, wl_grp, wl_grp_edges, plots_dir / f'{prefix}sdec_spectrum.png', model_name, v=v)
    else:
        print("  WARNING: No output.grp_grid or input.str found for spectral decomposition")

    print(f"\nAll plots saved to: {plots_dir}")


if __name__ == '__main__':
    main()
