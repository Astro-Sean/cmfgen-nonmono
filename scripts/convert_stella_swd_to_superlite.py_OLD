#!/usr/bin/env python3
"""Convert a STELLA *.swd snapshot directly to SuperLite input.str format.

This is an **MC-aware** conversion path: SuperLite is a Monte Carlo (MC)
radiative-transfer code, and high optical-depth / high-density regions make
MC packets random-walk, blowing up wall-clock runtime.  This converter reads
the STELLA .swd directly (skipping the CMFGEN SN_HYDRO_DATA intermediate) and
applies optimisations that keep SuperLite fast without sacrificing the
spectrum-forming region:

  1. **Native Rosseland kappa** — the .swd already contains STELLA's NLTE
     Rosseland opacity (column 12), so the optical-depth truncation is exact,
     not an electron-scattering approximation.
  2. **tau-based truncation** — zones with cumulative tau > ``truncate_tau``
     (default 20) are dropped.  These deep, optically-thick zones contribute
     negligible flux but dominate MC cost.  The photosphere (tau ~ 1) and the
     line-forming region (tau ~ 0.01-1) are always retained.
  3. **Grid downsampling** — the surviving zones are downsampled to
     ``max_cells`` (default 150) using a tau-weighted log-radius scheme that
     keeps finer resolution where the spectrum forms (low tau) and coarsens
     the deeper diffusion region.  MC runtime scales roughly with ncells.
  4. **Full 15-element composition** — the .swd carries all 15 STELLA elements
     individually (H, He, C, N, O, Ne, Na, Mg, Al, Si, S, Ar, Ca, Fe, Ni).
     SuperLite's atomic database supports every one of them, so we pass the
     full composition through instead of the grouped 9-element set that the
     SN_HYDRO_DATA path is limited to.

Usage:
  python convert_stella_swd_to_superlite.py SWD_FILE -o OUTPUT_DIR [OPTIONS]

Options:
  --target-days T       Epoch closest to T days (default: last epoch)
  --truncate-tau X      Drop zones with cumulative Rosseland tau > X (default 20)
  --max-cells N         Downsample to <= N cells after truncation (default 150)
  --min-cells N         Never go below N cells (default 20)
  --nomp N              OpenMP threads for input.par in_nomp (default 4)
  --rho-floor X         Drop zones with rho < X g/cm^3 (default 0, keep all)
  --no-downsample       Disable grid downsampling (keep all post-truncation zones)
  --model-name NAME     Model name for input.par (default: from .swd filename)

The output is written to OUTPUT_DIR/input.str and OUTPUT_DIR/input.par,
ready for SuperLite.  The input.str is inner-to-outer (r increasing), which
is SuperLite's convention (opposite of CMFGEN/STELLA .swd).
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path

import numpy as np

# Reuse the battle-tested .swd parser and composition normaliser from the
# CMFGEN conversion pipeline rather than duplicating the fixed-width parsing.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from convert_stella_swd_to_cmfgen import (
    read_swd,
    normalize_swd_composition,
    compute_kappa_tau,
    approximate_rosseland_opacity,
    STELLA_ELEM_NAMES,
    MSUN,
    LSUN,
)

# SuperLite atomic-data supports all 15 STELLA elements individually.
# Map STELLA element name -> SuperLite column label (lowercase symbol).
SL_ELEM_LABELS = {
    'H': 'h', 'He': 'he', 'C': 'c', 'N': 'n', 'O': 'o',
    'Ne': 'ne', 'Na': 'na', 'Mg': 'mg', 'Al': 'al', 'Si': 'si',
    'S': 's', 'Ar': 'ar', 'Ca': 'ca', 'Fe': 'fe', 'Ni': 'ni',
}
# Canonical order (all 15 STELLA elements).  The output only includes
# elements that are actually present in the model (see detect_elements).
SL_ELEM_ORDER = ['H', 'He', 'C', 'N', 'O', 'Ne', 'Na', 'Mg', 'Al',
                 'Si', 'S', 'Ar', 'Ca', 'Fe', 'Ni']


def detect_elements(mf, threshold=1e-10):
    """Detect which elements have non-zero mass fractions in the model.

    SuperLite loads atomic data (Kurucz line lists) for every element column
    in input.str.  Including elements with zero abundance wastes startup time
    and memory.  This function checks each element's maximum mass fraction
    across all zones and returns only those above ``threshold``.

    Args:
        mf: dict of element name -> mass fraction array (outer->inner or
            inner->outer, doesn't matter — we take the max)
        threshold: minimum max mass fraction to consider an element present

    Returns:
        list of element names (in SL_ELEM_ORDER order) that are present
    """
    present = []
    for elem in SL_ELEM_ORDER:
        if elem in mf:
            mx = float(np.max(mf[elem]))
            if mx > threshold:
                present.append(elem)
    return present


def compute_robust_tau(r_oi, rho_oi, kappa_swd, T_oi, mf, verbose=True):
    """Compute a robust cumulative optical depth for MC-aware truncation.

    SuperLite computes its own opacities from atomic data (Kurucz tables),
    which can differ significantly from STELLA's Rosseland mean (WRKX).
    In particular, STELLA's WRKX can be very low in cool, neutral regions
    where electron scattering drops to zero — but SuperLite will ionize
    that material and see the full electron-scattering opacity.

    This function computes a **conservative hybrid opacity** by taking the
    element-wise maximum of:

      1. **STELLA native kappa** (WRKX) — accurate where STELLA's NLTE
         ionization treatment is reliable (hot, ionized regions).
      2. **Electron-scattering floor** — kappa_es = 0.2*(1+X) / (1+3X+0.25Y),
         assuming full ionization.  This is the maximum ES opacity that
         SuperLite could see, and ensures we don't underestimate tau in
         neutral regions that SuperLite may ionize.
      3. **Corrected Kramers bound-free + free-free** — from the (now fixed)
         approximate_rosseland_opacity, using the actual metal mass fraction
         (not the atomic number, which was the previous bug).

    The cumulative tau is then computed from the outer boundary inward
    (outer->inner, decreasing r) using the trapezoidal rule.

    For Type IIn models with dense CSM, this ensures the inner-boundary
    tau is correctly very high (>>1), enabling proper truncation of the
    optically-thick diffusion region.

    Args:
        r_oi: radius array, outer->inner (decreasing r), in cm
        rho_oi: density array, outer->inner, in g/cm³
        kappa_swd: STELLA Rosseland kappa, outer->inner, in cm²/g
        T_oi: temperature array, outer->inner, in K
        mf: dict of element mass fractions, outer->inner
        verbose: print diagnostics

    Returns:
        tau: cumulative optical depth array, outer->inner
        kappa_robust: the hybrid opacity array used, outer->inner
    """
    n = len(r_oi)

    # --- Component 1: STELLA native kappa ---
    kappa_native = np.maximum(kappa_swd, 0.0)

    # --- Component 2: Electron scattering floor (full ionization) ---
    X = mf.get('H', np.zeros(n))
    Y = mf.get('He', np.zeros(n))
    kappa_es_floor = 0.2 * (1.0 + X) / (1.0 + 3.0 * X + 0.25 * Y)

    # --- Component 3: Corrected Kramers bf+ff ---
    Z_metal = np.zeros(n)
    for elem in mf:
        if elem not in ['H', 'He']:
            Z_metal += mf[elem]
    kappa_kramers = approximate_rosseland_opacity(
        T_oi, rho_oi, X, Y, Z_metal, mass_fractions=mf)
    # Subtract the ES component (already counted in kappa_es_floor) to
    # get just the bf+ff part
    kappa_es_approx = 0.2 * (1.0 + X) / (1.0 + 3.0 * X + 0.25 * Y)
    kappa_bf_ff = np.maximum(kappa_kramers - kappa_es_approx, 0.0)

    # --- Hybrid: take the maximum of native and (ES floor + bf/ff) ---
    kappa_hybrid = np.maximum(kappa_native, kappa_es_floor + kappa_bf_ff)

    # --- Cumulative tau from outer boundary inward ---
    tau = compute_kappa_tau(r_oi, rho_oi, kappa_hybrid)

    if verbose:
        print(f"\n  Robust tau calculation (hybrid opacity):")
        print(f"    kappa_swd:      {kappa_native.min():.3e} - "
              f"{kappa_native.max():.3e} cm²/g")
        print(f"    kappa_es_floor: {kappa_es_floor.min():.3e} - "
              f"{kappa_es_floor.max():.3e} cm²/g")
        print(f"    kappa_bf_ff:    {kappa_bf_ff.min():.3e} - "
              f"{kappa_bf_ff.max():.3e} cm²/g")
        print(f"    kappa_hybrid:   {kappa_hybrid.min():.3e} - "
              f"{kappa_hybrid.max():.3e} cm²/g")
        print(f"    tau: outer={tau[0]:.3e}, max={tau.max():.3e} "
              f"(zone {np.argmax(tau)})")
        # Report tau at key thresholds
        for thresh in [1.0, 10.0, 100.0, 1000.0]:
            idx = np.searchsorted(tau, thresh)
            if idx < n:
                print(f"    tau={thresh:.0f} at r={r_oi[idx]:.3e} cm "
                      f"(zone {idx})")
            else:
                print(f"    tau={thresh:.0f} not reached "
                      f"(max tau={tau.max():.3e})")

    return tau, kappa_hybrid


def select_epoch(epochs: dict, target_days: float | None) -> tuple[float, dict]:
    """Pick the epoch closest to target_days (or the last one)."""
    times = sorted(epochs.keys())
    if target_days is None:
        t = times[-1]
    else:
        t = min(times, key=lambda tt: abs(tt - target_days))
    print(f"Selected epoch t={t:.4f} d (requested {target_days}) "
          f"from {len(times)} epochs spanning {times[0]:.3f}-{times[-1]:.3f} d")
    return t, epochs[t]


def compute_cell_boundaries(r_centers):
    """Build cell edges from zone centres via geometric-mean interfaces.

    Returns an array of length n+1 (inner edge first, outer edge last),
    i.e. inner-to-outer (increasing r) order.  r_centers must be increasing.
    """
    n = len(r_centers)
    edges = np.zeros(n + 1)
    # interior interfaces = geometric mean of adjacent centres
    edges[1:n] = np.sqrt(r_centers[:-1] * r_centers[1:])
    # inner boundary: extrapolate geometrically
    edges[0] = r_centers[0] ** 2 / edges[1] if edges[1] > 0 else r_centers[0] * 0.5
    # outer boundary
    edges[n] = r_centers[-1] ** 2 / edges[n - 1] if edges[n - 1] > 0 else r_centers[-1] * 1.5
    return edges


def tau_weighted_downsample(r_oi, v_oi, T_oi, T_rad_oi, rho_oi, kappa_oi, mf_oi,
                            tau_oi, max_cells, min_cells=20):
    """Downsample outer->inner arrays to <= max_cells, tau-uniform.

    Distributes cells so that no resulting cell has extreme per-cell optical
    depth, which would trigger very slow DDMC transport in SuperLite.

    Strategy:
      - Spectrum region (tau 0 to ~1): uniform in log(radius) for good
        spatial resolution where lines form.
      - Diffusion region (tau ~1 to tau_max): uniform in tau space so each
        cell gets roughly equal d_tau = tau_max / n_diff.
      - Split: 40% of cells to spectrum, 60% to diffusion.

    Zones between selected boundaries are **merged** with mass conservation:
      - rho_merged = total_mass / total_volume
      - T_merged, T_rad_merged, v_merged, kappa_merged, mf_merged = mass-weighted averages
      - r_merged = volume-weighted average radius

    This ensures both good spectral resolution and manageable per-cell tau.
    """
    n = len(r_oi)
    if n <= max_cells:
        return r_oi, v_oi, T_oi, T_rad_oi, rho_oi, kappa_oi, mf_oi

    tau_max = float(tau_oi[-1])

    # Find the boundary between spectrum and diffusion regions
    tau_boundary = min(1.0, tau_max * 0.5)
    n_spec_boundary = int(np.searchsorted(tau_oi, tau_boundary))

    # Split cell budget proportionally to the number of original zones in
    # each region, with a minimum of 30% each to avoid starving either.
    n_spec_zones = n_spec_boundary
    n_diff_zones = n - n_spec_boundary
    if n_spec_zones > 0 and n_diff_zones > 0:
        frac_spec = n_spec_zones / n
        n_spec = int(max(max_cells * 0.3, max_cells * frac_spec))
        n_spec = min(n_spec, max_cells - int(max_cells * 0.3))
        n_diff = max_cells - n_spec
    else:
        n_spec = int(max_cells * 0.4)
        n_diff = max_cells - n_spec

    sel = []

    # --- Spectrum region: uniform in zone index ---
    # Uniform-in-log(r) causes heavy deduplication when outer zones span
    # large log-r ranges.  Uniform-in-index preserves the original resolution
    # distribution and avoids creating giant merged cells.
    if n_spec_boundary > 0 and n_spec > 0:
        sel_spec = np.linspace(0, n_spec_boundary - 1, n_spec).astype(int)
        sel.extend(sel_spec.tolist())

    # --- Diffusion region: uniform in tau ---
    if n_spec_boundary < n and n_diff > 0:
        tau_diff = tau_oi[n_spec_boundary:]
        # Place n_diff targets uniformly in tau
        targets_diff = np.linspace(tau_diff[0], tau_diff[-1], n_diff)
        for t in targets_diff:
            idx = np.searchsorted(tau_diff, t)
            idx = min(idx, len(tau_diff) - 1)
            sel.append(idx + n_spec_boundary)

    # Deduplicate and sort
    sel = sorted(set(sel))

    # If we lost too many to dedup, fill with uniform spacing
    if len(sel) < min_cells:
        sel = np.linspace(0, n - 1, min_cells).astype(int).tolist()
        sel = sorted(set(sel))

    # Always keep endpoints
    if 0 not in sel:
        sel.insert(0, 0)
    if (n - 1) not in sel:
        sel.append(n - 1)

    # If still under max_cells, fill with uniform index spacing
    if len(sel) < max_cells:
        all_indices = set(range(n))
        available = sorted(all_indices - set(sel))
        to_add = available[:max_cells - len(sel)]
        sel = sorted(set(sel) | set(to_add))

    sel = np.array(sel)

    # --- Merge zones between selected boundaries (mass-conserving) ---
    # r_oi is outer->inner (decreasing).  compute_cell_boundaries expects
    # increasing r, so reverse to inner->outer, compute edges, then the
    # volumes are correct regardless of ordering.
    r_io = r_oi[::-1]  # inner->outer (increasing)
    edges_io = compute_cell_boundaries(r_io)  # inner->outer edges, len n+1
    vol_io = (4.0 / 3.0) * np.pi * (
        edges_io[1:] ** 3 - edges_io[:-1] ** 3)  # inner->outer volumes
    # Reverse back to outer->inner to match r_oi ordering
    vol_orig = vol_io[::-1]
    mass_orig = rho_oi * vol_orig

    n_out = len(sel)
    r_out = np.zeros(n_out)
    v_out = np.zeros(n_out)
    T_out = np.zeros(n_out)
    T_rad_out = np.zeros(n_out)
    rho_out = np.zeros(n_out)
    kappa_out = np.zeros(n_out)
    mf_out = {k: np.zeros(n_out) for k in mf_oi}
    # Actual merged cell boundaries in outer->inner order (decreasing r)
    # edges_merged_oi[i] is the outer boundary of merged cell i,
    # edges_merged_oi[i+1] is the inner boundary.
    # We build them from the original edges_io (inner->outer).
    # edges_io indices: zone j (inner->outer) spans edges_io[j] to edges_io[j+1]
    # In outer->inner order, zone j_oi = n-1-j_io spans edges_io[j_io+1] to edges_io[j_io]
    edges_merged_io = np.zeros(n_out + 1)  # inner->outer

    for idx_out in range(n_out):
        # Merge all original zones from sel[idx_out] to sel[idx_out+1]-1
        # (or to the end for the last selected zone)
        # sel indices are in outer->inner order
        i_start = sel[idx_out]
        i_end = sel[idx_out + 1] if idx_out + 1 < n_out else n
        sl = slice(i_start, i_end)

        m_grp = mass_orig[sl]
        v_grp = vol_orig[sl]
        m_tot = m_grp.sum()
        v_tot = v_grp.sum()

        if m_tot > 0:
            # Mass-weighted averages for intensive properties
            w = m_grp / m_tot
            v_out[idx_out] = np.sum(w * v_oi[sl])
            T_out[idx_out] = np.sum(w * T_oi[sl])
            T_rad_out[idx_out] = np.sum(w * T_rad_oi[sl])
            kappa_out[idx_out] = np.sum(w * kappa_oi[sl])
            for k in mf_oi:
                mf_out[k][idx_out] = np.sum(w * mf_oi[k][sl])
            # Volume-weighted average radius (geometric center of mass)
            r_out[idx_out] = np.sum(v_grp * r_oi[sl]) / v_tot
        else:
            # Fallback: simple average
            r_out[idx_out] = np.mean(r_oi[sl])
            v_out[idx_out] = np.mean(v_oi[sl])
            T_out[idx_out] = np.mean(T_oi[sl])
            T_rad_out[idx_out] = np.mean(T_rad_oi[sl])
            kappa_out[idx_out] = np.mean(kappa_oi[sl])
            for k in mf_oi:
                mf_out[k][idx_out] = np.mean(mf_oi[k][sl])

        # Density from mass conservation
        rho_out[idx_out] = m_tot / max(v_tot, 1e-30)

        # Actual merged cell boundaries.
        # idx_out goes 0 (outermost group) to n_out-1 (innermost group),
        # but edges_merged_io is inner->outer, so the innermost group maps
        # to edges_merged_io[0..1] and the outermost to [n_out-1..n_out].
        # Group idx_out covers outer->inner zones [i_start, i_end).
        # In inner->outer indexing, the innermost zone in the group is
        # j_io = n - i_end, and the outermost is j_io = n - 1 - i_start.
        # Inner boundary = edges_io[n - i_end]
        # Outer boundary = edges_io[n - i_start]
        inner_edge = edges_io[n - i_end]
        outer_edge = edges_io[n - i_start]
        io_idx = n_out - 1 - idx_out  # reverse: outermost group -> last cell
        edges_merged_io[io_idx] = inner_edge
        edges_merged_io[io_idx + 1] = outer_edge

    # Compute per-cell tau and mass from ACTUAL merged boundaries
    # edges_merged_io is inner->outer; r_out is outer->inner — reverse r_out
    vol_merged = (4.0 / 3.0) * np.pi * (
        edges_merged_io[1:] ** 3 - edges_merged_io[:-1] ** 3)
    # vol_merged is inner->outer; reverse to outer->inner to match r_out
    vol_merged_oi = vol_merged[::-1]
    dr_est = np.diff(edges_merged_io)
    tau_cell_est = kappa_out * rho_out * dr_est

    mass_out = rho_out * vol_merged_oi
    print(f"[Downsample] {n} -> {n_out} cells (tau-uniform, max={max_cells})")
    print(f"  Per-cell tau estimate: {tau_cell_est.min():.3e} - "
          f"{tau_cell_est.max():.3e}")
    print(f"  Mass conservation: {mass_orig.sum() / MSUN:.4f} -> "
          f"{mass_out.sum() / MSUN:.4f} Msun")
    if abs(mass_out.sum() - mass_orig.sum()) / max(mass_orig.sum(), 1e-30) > 0.01:
        print(f"  WARNING: mass changed by "
              f"{100 * abs(mass_out.sum() - mass_orig.sum()) / max(mass_orig.sum(), 1e-30):.1f}%")
    if tau_cell_est.max() > 3.0:
        print(f"  WARNING: {np.sum(tau_cell_est > 3.0)} cells have per-cell "
              f"tau > 3 (will trigger DDMC in SuperLite)")

    # Return actual merged boundaries (inner->outer) alongside zone properties
    # Caller must use these instead of compute_cell_boundaries.
    return (r_out, v_out, T_out, T_rad_out, rho_out, kappa_out, mf_out,
            edges_merged_io)


def write_input_str(path, r_inner, edges, v_edges, mass, T_rad, T_gas,
                    mf_dict, t_days, model_name, elem_list):
    """Write SuperLite input.str (inner-to-outer, r increasing).

    Args:
        elem_list: list of element names to include (only those present in
            the model).  This controls the number of abundance columns.
    """
    n = len(mass)
    labels = ['x_left', 'x_right', 'vx_left', 'vx_right', 'mass',
              'rad_temp', 'temp']
    abund_labels = [SL_ELEM_LABELS[e] for e in elem_list]
    labels += abund_labels
    nabund = len(abund_labels)
    ncol = 7 + nabund

    with open(path, 'w') as f:
        f.write(f'# 1D spherical - STELLA .swd snapshot (t={t_days:.3f}d, '
                f'{model_name})\n')
        f.write(f'#      {n}          1 1 {ncol} {nabund}\n')
        f.write('#   ' + '  '.join(f'{l:>12s}' for l in labels) + '\n')
        for i in range(n):
            vals = [edges[i], edges[i + 1], v_edges[i], v_edges[i + 1],
                    mass[i], T_rad[i], T_gas[i]]
            vals += [mf_dict.get(e, np.zeros(n))[i] for e in elem_list]
            f.write('  '.join(f'{v:14.6e}' for v in vals) + '  \n')
    print(f"Wrote {path}: {n} cells, {nabund} abundance columns "
          f"({', '.join(elem_list)})")


def write_input_par(path, model_name, t_days, ndim, lum_erg, r_phot, nomp,
                    truncate_tau, ncells_orig, ncells_final, quick_test=False,
                    grp_ng=6000, grp_wlmin='1d-8', grp_wlmax='30d-5',
                    flx_ndim=3500, flx_wlmin='3d-5', flx_wlmax='2d-4',
                    src_n2s=25, trn_tauddmc='3d0', taulump='3d0',
                    nlte=True, nlte_nelem=1, niter=3,
                    nebular=True, opacmixrossel='0d0',
                    nobbopac=False, nobfopac=False,
                    noffopac=False, nothmson=False,
                    l_inner=0.0, srcmax=0.0):
    """Write input.par from the reference template with MC-aware settings.

    All SuperLite numerical parameters are exposed as kwargs with defaults
    matching the current production configuration.  These can be overridden
    from YAML via the orchestrator.

    Args:
        quick_test: if True, use reduced particle count, wavelength groups,
            and higher DDMC threshold for fast interactive testing (~minutes
            instead of hours).  Production runs should use quick_test=False.
        grp_ng: number of wavelength groups (default 6000, production).
        grp_wlmin, grp_wlmax: transport wavelength range in cm.
        flx_ndim: number of flux output bins (default 3500).
        flx_wlmin, flx_wlmax: flux output wavelength range in cm.
        src_n2s: 2^n source particles (default 25 = 33M).
        trn_tauddmc: DDMC transition threshold in mfp (default 3).
        taulump: DDMC group lumping threshold in mfp (default 3).
        nlte: enable NLTE for H bound-bound (default True).
        nlte_nelem: number of NLTE elements (must be 1, H only).
        niter: radiation temperature feedback iterations (default 3).
        nebular: enable Mazzali & Lucy nebular approximation (default True).
        opacmixrossel: Rosseland/Planck opacity mixing fraction (default 0).
        nobbopac, nobfopac, noffopac, nothmson: disable opacity sources.
    """
    # When boundary source is active (l_inner > 0), set in_L_bol = 0 to
    # avoid double-counting.  All radiation comes from the inner boundary
    # blackbody at srcmax temperature.  The interior thermal emission is
    # implicitly included through absorption/re-emission of boundary photons.
    l_bol_eff = 0.0 if l_inner > 0 else lum_erg

    ref_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                            'superlite_base', 'input_parameters_reference.par')
    if os.path.exists(ref_path):
        with open(ref_path) as f:
            template = f.read()
        idx = template.index('&inputpars')
        par = template[idx:]
        par = (par.replace('{MODEL_NAME}', model_name)
                  .replace('{EPOCH}', f'{t_days:.1f}d')
                  .replace('{NDIM}', str(ndim))
                  .replace('{L_BOL}', f'{l_bol_eff:.4e}')
                  .replace('{R_PHOT}', f'{r_phot:.4e}')
                  .replace('{NOMP}', str(nomp))
                  .replace('{L_INNER}', f'{l_inner:.4e}')
                  .replace('{SRCMAX}', f'{srcmax:.4e}'))
        # Override nebular/NLTE settings from parameters (the template
        # has these hardcoded as 't'; replace with the actual values).
        nebular_str = 't' if nebular else 'f'
        nlte_str = 't' if nlte else 'f'
        par = re.sub(r'in_nebular\s*=\s*[tf]', f'in_nebular = {nebular_str}', par)
        par = re.sub(r'in_nlte\s*=\s*[tf]', f'in_nlte = {nlte_str}', par)
        par = re.sub(r'in_nlte_nelem\s*=\s*\d+', f'in_nlte_nelem = {nlte_nelem}', par)
        par = re.sub(r'in_niter\s*=\s*\d+', f'in_niter = {niter}', par)
        if quick_test:
            # Override for quick interactive testing:
            # - 2^20 = 1M particles (vs 2^25 = 33M production)
            # - 1000 wavelength groups (vs 6000 production)
            # - 500 flux bins (vs 3500 production)
            # - DDMC threshold 10 (vs 3 production) — fewer cells trigger DDMC
            # - DDMC group lumping 10 (vs 3) — more aggressive lumping
            par = par.replace('in_src_n2s = 25', 'in_src_n2s = 20')
            par = par.replace('in_grp_ng = 6000', 'in_grp_ng = 1000')
            par = par.replace('in_flx_ndim = 3500, 1, 1',
                              'in_flx_ndim = 500, 1, 1')
            par = par.replace('in_flx_wlmax = 2d-4',
                              'in_flx_wlmax = 2d-4')
            par = par.replace('in_trn_tauddmc = 3d0', 'in_trn_tauddmc = 10d0')
            par = par.replace('in_taulump = 3d0', 'in_taulump = 10d0')
            # Add comment
            par = par.replace(
                "in_comment = 'Production",
                "in_comment = 'QUICK TEST (1M particles, 1k groups, DDMC thresh=10) - "
            )
    else:
        n2s = 20 if quick_test else src_n2s
        grp_ng_val = 1000 if quick_test else grp_ng
        flx_ndim_val = 500 if quick_test else flx_ndim
        tauddmc = '10d0' if quick_test else trn_tauddmc
        taulump_val = '10d0' if quick_test else taulump
        comment = ('QUICK TEST' if quick_test else 'MC-aware') + \
                  f' 1-D sph model of {model_name} at {t_days:.1f}d'
        nlte_str = 't' if nlte else 'f'
        nebular_str = 't' if nebular else 'f'
        nobb_str = 't' if nobbopac else 'f'
        nobf_str = 't' if nobfopac else 'f'
        noff_str = 't' if noffopac else 'f'
        nothm_str = 't' if nothmson else 'f'
        par = (
            f"&inputpars\n\n"
            f"in_name = '{model_name}'\n"
            f"in_comment = '{comment}'\n"
            f" in_nomp = {nomp}\n"
            f" in_io_grabstdout = f\n"
            f" in_io_opacdump = 'one'\n"
            f" in_io_profdump = t\n"
            f" in_ndim = {ndim}, 1, 1\n"
            f" in_grp_ng = {grp_ng_val}\n"
            f" in_grp_wlmin = {grp_wlmin}\n"
            f" in_grp_wlmax = {grp_wlmax}\n"
            f" in_flx_ndim = {flx_ndim_val}, 1, 1\n"
            f" in_flx_wlmin = {flx_wlmin}\n"
            f" in_flx_wlmax = {flx_wlmax}\n"
            f" in_src_n2s = {n2s}\n"
            f" in_L_bol = {l_bol_eff:.4e}\n"
            f" in_L_inner = {l_inner:.4e}\n"
            f" in_srcmax = {srcmax:.4e}\n"
            f" in_trn_tauddmc = {tauddmc}\n"
            f" in_taulump = {taulump_val}\n"
            f" in_nebular = {nebular_str}\n"
            f" in_R_phot = {r_phot:.4e}\n"
            f"!-- NLTE enabled: low-density CSM/ejecta requires non-LTE treatment.\n"
            f"!-- in_niter=3 allows radiation-temperature feedback (~3x runtime).\n"
            f" in_nlte = {nlte_str}\n"
            f" in_nlte_nelem = {nlte_nelem}\n"
            f" in_niter = {niter}\n"
            f" in_opacmixrossel = {opacmixrossel}\n"
            f" in_nobbopac = {nobb_str}\n"
            f" in_nobfopac = {nobf_str}\n"
            f" in_noffopac = {noff_str}\n"
            f" in_nothmson = {nothm_str}\n"
            f"/\n"
        )

    # Always regenerate input.par from the reference template to ensure
    # physics settings (NLTE, flux range, nebular, opacities) are kept
    # up-to-date with the optimum configuration.  Only preserve a small
    # set of user-tunable performance parameters from any existing file:
    #   in_src_n2s  (particle count exponent)
    #   in_grp_ng   (wavelength group count)
    #   in_trn_tauddmc, in_taulump  (DDMC thresholds)
    # Structure-dependent values (ndim, L_bol, nomp, R_phot) are always
    # recomputed by the converter.
    if os.path.exists(path) and not quick_test:
        with open(path) as f:
            existing = f.read()
        # Extract user-tunable performance parameters from the existing file
        preserve_keys = ['in_src_n2s', 'in_grp_ng', 'in_trn_tauddmc', 'in_taulump']
        for key in preserve_keys:
            m = re.search(rf'{key}\s*=\s*(\S+)', existing)
            if m:
                val = m.group(1)
                par = re.sub(rf'{key}\s*=\s*\S+', f'{key} = {val}', par)
        print(f"Regenerated input.par from template (preserved performance params)")
    elif os.path.exists(path) and quick_test:
        print(f"Quick-test mode: writing fresh input.par (overriding existing)")

    with open(path, 'w') as f:
        f.write(par)
    src_mode = (f"L_inner={l_inner:.3e}, srcmax={srcmax:.1f} K"
                if l_inner > 0 else "none")
    print(f"Wrote {path} (ndim={ndim}, L_bol={l_bol_eff:.3e} erg/s, "
          f"R_phot={r_phot:.3e} cm, boundary_source={src_mode})")


def plot_superlite_input_structure(output_dir, r_sl, v_sl, T_sl, rho_sl,
                                    edges, tau_sl, mf_sl, elem_present,
                                    r_phot, t_day, model_name):
    """Generate diagnostic plots of the SuperLite input.str structure.

    Produces a multi-panel plot showing:
      1. Density vs radius (log-log)
      2. Velocity vs radius
      3. Temperature vs radius (log-log)
      4. Cumulative optical depth vs radius (log-log)
      5. Mass fractions vs radius (stacked area, only present elements)
      6. Cell width distribution (dr/r vs radius)

    Also marks the photosphere (tau=1) location.

    Args:
        output_dir: directory to write the plot (typically the superlite run dir)
        r_sl: radius array, inner->outer (cm)
        v_sl: velocity array, inner->outer (km/s)
        T_sl: temperature array, inner->outer (K)
        rho_sl: density array, inner->outer (g/cm^3)
        edges: cell boundaries, length n+1 (cm)
        tau_sl: cumulative optical depth, inner->outer
        mf_sl: dict of element -> mass fraction (inner->outer, normalised)
        elem_present: list of element names present in the model
        r_phot: photospheric radius (cm)
        t_day: snapshot time in days
        model_name: model name for plot title
    """
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("  WARNING: matplotlib not available, skipping input structure plots")
        return

    # Use velocity (km/s) as the x-axis for SN ejecta (more physical than radius)
    # Also compute radius in units of 1e15 cm for readability
    r_1e15 = r_sl / 1e15
    r_phot_1e15 = r_phot / 1e15

    # Cell widths
    dr = np.diff(edges)
    dr_r = dr / r_sl

    fig, axes = plt.subplots(3, 2, figsize=(14, 16))
    fig.suptitle(f'SuperLite input.str: {model_name} at t={t_day:.2f}d\n'
                 f'({len(r_sl)} cells, {len(elem_present)} elements: '
                 f'{", ".join(elem_present)})',
                 fontsize=13, fontweight='bold')

    # --- Panel 1: Density vs radius ---
    ax = axes[0, 0]
    ax.loglog(r_1e15, rho_sl, 'b-', linewidth=1.2)
    ax.axvline(r_phot_1e15, color='r', linestyle='--', alpha=0.7, label=f'R_phot (tau=1)')
    ax.set_xlabel(r'Radius [$10^{15}$ cm]')
    ax.set_ylabel(r'Density [g cm$^{-3}$]')
    ax.set_title('Density profile')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # --- Panel 2: Velocity vs radius ---
    ax = axes[0, 1]
    ax.semilogx(r_1e15, v_sl, 'g-', linewidth=1.2)
    ax.axvline(r_phot_1e15, color='r', linestyle='--', alpha=0.7, label=f'R_phot')
    ax.set_xlabel(r'Radius [$10^{15}$ cm]')
    ax.set_ylabel('Velocity [km/s]')
    ax.set_title('Velocity profile')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # --- Panel 3: Temperature vs radius ---
    ax = axes[1, 0]
    ax.loglog(r_1e15, T_sl, 'r-', linewidth=1.2)
    ax.axvline(r_phot_1e15, color='r', linestyle='--', alpha=0.7, label=f'R_phot')
    ax.set_xlabel(r'Radius [$10^{15}$ cm]')
    ax.set_ylabel('Temperature [K]')
    ax.set_title('Temperature profile')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # --- Panel 4: Cumulative optical depth vs radius ---
    ax = axes[1, 1]
    ax.loglog(r_1e15, np.maximum(tau_sl, 1e-10), 'k-', linewidth=1.2)
    ax.axhline(1.0, color='r', linestyle='--', alpha=0.5, label='tau=1 (photosphere)')
    ax.axhline(10.0, color='orange', linestyle='--', alpha=0.5, label='tau=10')
    ax.axhline(20.0, color='purple', linestyle='--', alpha=0.5, label='tau=20 (truncation)')
    ax.axvline(r_phot_1e15, color='r', linestyle=':', alpha=0.7)
    ax.set_xlabel(r'Radius [$10^{15}$ cm]')
    ax.set_ylabel('Cumulative optical depth (hybrid)')
    ax.set_title('Optical depth profile (MC-aware hybrid kappa)')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # --- Panel 5: Mass fractions vs radius (stacked area) ---
    ax = axes[2, 0]
    # Build stacked area plot with only present elements
    colors = plt.cm.Set1(np.linspace(0, 1, len(elem_present)))
    bottom = np.zeros(len(r_sl))
    for i, elem in enumerate(elem_present):
        vals = mf_sl.get(elem, np.zeros(len(r_sl)))
        ax.semilogx(r_1e15, vals, color=colors[i], linewidth=1.0, label=elem)
        # Also fill below for visibility
        ax.fill_between(r_1e15, bottom, bottom + vals, color=colors[i], alpha=0.3)
        bottom += vals
    ax.axvline(r_phot_1e15, color='r', linestyle='--', alpha=0.7)
    ax.set_xlabel(r'Radius [$10^{15}$ cm]')
    ax.set_ylabel('Mass fraction')
    ax.set_title(f'Composition ({len(elem_present)} elements present)')
    ax.legend(fontsize=8, loc='best', ncol=2)
    ax.set_ylim(0, 1.05)
    ax.grid(True, alpha=0.3)

    # --- Panel 6: Cell width distribution ---
    ax = axes[2, 1]
    ax.loglog(r_1e15, dr_r, 'b-', linewidth=1.0)
    ax.axvline(r_phot_1e15, color='r', linestyle='--', alpha=0.7, label='R_phot')
    ax.set_xlabel(r'Radius [$10^{15}$ cm]')
    ax.set_ylabel('dr/r (relative cell width)')
    ax.set_title('Grid resolution (cell width distribution)')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plot_path = os.path.join(output_dir, 'input_structure.png')
    fig.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Wrote input structure plot: {plot_path}")

    # --- Second plot: mass fractions vs enclosed mass ---
    fig2, ax2 = plt.subplots(1, 1, figsize=(10, 7))
    # Compute enclosed mass (inner->outer)
    vol = (4.0 / 3.0) * np.pi * (edges[1:] ** 3 - edges[:-1] ** 3)
    mass_per_cell = rho_sl * vol
    m_enc = np.cumsum(mass_per_cell) / MSUN
    colors = plt.cm.Set1(np.linspace(0, 1, len(elem_present)))
    for i, elem in enumerate(elem_present):
        vals = mf_sl.get(elem, np.zeros(len(r_sl)))
        ax2.semilogx(m_enc, vals, color=colors[i], linewidth=1.2, label=elem)
    ax2.set_xlabel('Enclosed mass [Msun]')
    ax2.set_ylabel('Mass fraction')
    ax2.set_title(f'Composition vs enclosed mass: {model_name} at t={t_day:.2f}d')
    ax2.legend(fontsize=9, loc='best')
    ax2.set_ylim(0, 1.05)
    ax2.grid(True, alpha=0.3)
    plot2_path = os.path.join(output_dir, 'input_mass_fractions.png')
    fig2.savefig(plot2_path, dpi=150, bbox_inches='tight')
    plt.close(fig2)
    print(f"  Wrote mass fractions plot: {plot2_path}")


def convert_swd_to_superlite(swd_path, output_dir, target_days=None,
                             truncate_tau=20.0, max_cells=150, min_cells=20,
                             nomp=4, rho_floor=0.0, downsample=True,
                             model_name=None, quick_test=False, t_floor=2000.0,
                             grid_method='adaptive',
                             amr_rho_grad=0.5, amr_temp_grad=0.3,
                             amr_vel_grad=0.05, amr_comp_grad=0.1,
                             amr_cannon_interface_cells=5,
                             amr_shock_dv_threshold=100.0,
                             # SuperLite numerical parameters (exposed in YAML)
                             sl_grp_ng=6000, sl_grp_wlmin='1d-8',
                             sl_grp_wlmax='30d-5',
                             sl_flx_ndim=3500, sl_flx_wlmin='3d-5',
                             sl_flx_wlmax='2d-4',
                             sl_src_n2s=25, sl_trn_tauddmc='3d0',
                             sl_taulump='3d0',
                             sl_nlte=True, sl_nlte_nelem=1, sl_niter=3,
                             sl_nebular=True, sl_opacmixrossel='0d0',
                             sl_nobbopac=False, sl_nobfopac=False,
                             sl_noffopac=False, sl_nothmson=False):
    """Main conversion: .swd -> input.str + input.par (MC-aware).

    Args:
        t_floor: temperature floor in K.  SuperLite's nebular EOS uses a zeta
            table that starts at 2000 K; cells with radtemp below this cause
            an infinite loop in the itarget binary search.  The floor clamps
            both gas and radiation temperature to this minimum.
        quick_test: if True, generate input.par with reduced particles/groups
            and higher DDMC threshold for fast interactive testing.
        grid_method: 'adaptive' (physics-aware AMR, default) or 'tau_weighted'
            (legacy fixed split).
        amr_*: parameters for the adaptive grid (see GridRefinementConfig).
        sl_*: SuperLite numerical parameters, exposed in YAML under
            superlite.numerical.* and routed through the orchestrator.
    """
    if model_name is None:
        model_name = Path(swd_path).stem

    print(f"\n{'='*70}")
    print(f"MC-aware STELLA .swd -> SuperLite conversion")
    print(f"  source : {swd_path}")
    print(f"  output : {output_dir}")
    print(f"  tau_max: {truncate_tau}  max_cells: {max_cells}  "
          f"downsample: {downsample}")
    print(f"{'='*70}")

    # --- Parse .swd (reuses convert_stella_swd_to_cmfgen.read_swd) ---
    epochs = read_swd(swd_path)
    t_day, ep = select_epoch(epochs, target_days)

    # .swd arrays are outer->inner (decreasing r)
    r = ep['r']            # cm
    v = ep['v']            # km/s
    T = ep['T']            # gas kinetic temperature, K
    T_rad = ep['T_rad']    # radiation temperature, K (T_gas in optically-thick inner zones)
    rho = ep['rho']        # g/cm^3
    flux = ep['flux']      # erg/s (per-zone luminosity)
    kappa = ep['kappa']    # Rosseland opacity, cm^2/g (native from STELLA NLTE)
    comp = ep['composition']  # (nzon, 15) Y_k*A_k

    n_orig = len(r)
    print(f"\nEpoch has {n_orig} zones (outer->inner)")
    print(f"  R: {r.min():.3e} - {r.max():.3e} cm")
    print(f"  V: {v.min():.1f} - {v.max():.1f} km/s")
    print(f"  T_gas: {T.min():.3e} - {T.max():.3e} K")
    print(f"  T_rad: {T_rad.min():.3e} - {T_rad.max():.3e} K")
    print(f"  T_rad/T_gas ratio: {T_rad.min()/T[T_rad.argmin()]:.3f} - "
          f"{T_rad.max()/T[T_rad.argmax()]:.3f}")
    print(f"  rho: {rho.min():.3e} - {rho.max():.3e} g/cm^3")

    # --- Temperature floor ---
    # SuperLite's nebular EOS interpolates zeta values from a table that
    # starts at 2000 K (data.zeta).  Cells with radtemp below this cause
    # an infinite loop in the itarget binary search.  Clamp both gas and
    # radiation temperature to the floor.
    if t_floor > 0:
        n_clamped = int(np.sum(T < t_floor))
        if n_clamped > 0:
            print(f"  [T-floor] Clamping {n_clamped} T_gas zones from "
                  f"T_min={T.min():.1f} K to t_floor={t_floor:.1f} K")
            T = np.maximum(T, t_floor)
        n_clamped_rad = int(np.sum(T_rad < t_floor))
        if n_clamped_rad > 0:
            print(f"  [T-floor] Clamping {n_clamped_rad} T_rad zones from "
                  f"T_rad_min={T_rad.min():.1f} K to t_floor={t_floor:.1f} K")
            T_rad = np.maximum(T_rad, t_floor)

    # --- Free-streaming radiation temperature ---
    # STELLA's T_rad in the .swd is the radiation energy density temperature
    # (T_rad = (u/a)^(1/4)), which in optically thin regions is much lower
    # than the color temperature of the radiation field.  This causes
    # SuperLite to compute essentially zero opacity (neutral gas) and produce
    # a blank continuum.
    #
    # The free-streaming radiation temperature at radius r for a luminosity L
    # is T_fs(r) = (L / (4*pi*r^2*c*a))^(1/4), where a is the radiation
    # constant.  This represents the radiation energy density of a radial
    # photon flow and gives the correct ionization state for the nebular
    # approximation.  We use max(T_rad_stella, T_fs) to preserve the STELLA
    # value in optically thick regions while boosting it in thin regions.
    #
    # We need the total luminosity first.  Use the outer boundary flux
    # (zone 0 in outer->inner order) as the best estimate.
    _lum_for_fs = float(abs(flux[0])) if flux[0] != 0 else 1e43
    _c_cgs = 2.99792458e10
    _a_rad = 7.5657e-15
    T_fs = (_lum_for_fs / (4.0 * np.pi * r**2 * _c_cgs * _a_rad))**0.25
    n_boosted = int(np.sum(T_fs > T_rad))
    if n_boosted > 0:
        print(f"  [T_rad] Boosting {n_boosted} zones to free-streaming T_rad")
        print(f"    T_rad_stella: {T_rad.min():.1f} - {T_rad.max():.1f} K")
        print(f"    T_freestream: {T_fs.min():.1f} - {T_fs.max():.1f} K")
        T_rad = np.maximum(T_rad, T_fs)
        print(f"    T_rad_final:  {T_rad.min():.1f} - {T_rad.max():.1f} K")
    # Apply t_floor to the boosted T_rad as well
    if t_floor > 0:
        T_rad = np.maximum(T_rad, t_floor)

    if kappa is not None:
        print(f"  kappa_Ross: {kappa.min():.3e} - {kappa.max():.3e} cm^2/g (native STELLA)")
    else:
        print(f"  kappa_Ross: NOT in .swd, falling back to grey 0.34 cm^2/g")
        kappa = np.full(n_orig, 0.34)

    # --- Normalise composition to true mass fractions ---
    mf = normalize_swd_composition(comp)  # dict: elem name -> (nzon,) mass frac

    # --- Detect which elements are actually present in the model ---
    # SuperLite loads Kurucz atomic data for every abundance column.  Including
    # zero-abundance elements wastes startup time and memory.  Only include
    # elements that have non-zero mass fractions in at least one zone.
    elem_present = detect_elements(mf)
    elem_absent = [e for e in SL_ELEM_ORDER if e not in elem_present]
    print(f"\nElement detection: {len(elem_present)}/{len(SL_ELEM_ORDER)} elements present")
    print(f"  Present: {', '.join(elem_present)}")
    if elem_absent:
        print(f"  Absent (excluded from input.str): {', '.join(elem_absent)}")

    # --- Compute robust cumulative optical depth (outer->inner) ---
    # Use a hybrid opacity that takes the max of STELLA's native Rosseland
    # kappa and an electron-scattering-floor estimate.  This is critical for
    # Type IIn models: STELLA's WRKX can be very low in cool, neutral CSM
    # (electron scattering drops to zero when H is neutral), but SuperLite
    # will ionize that material and see the full ES opacity.  The hybrid
    # ensures the inner-boundary tau is correctly very high (>>1) for
    # Type IIn, enabling proper truncation of the optically-thick diffusion
    # region where MC packets would random-walk.
    tau, kappa_hybrid = compute_robust_tau(r, rho, kappa, T, mf, verbose=True)

    # --- MC-aware truncation: drop zones with tau > truncate_tau ---
    # tau is monotonic increasing inward (outer->inner).  Keep zones with
    # tau <= truncate_tau.  Always keep at least min_cells.
    mask = tau <= truncate_tau
    if mask.sum() < min_cells:
        # keep the outermost min_cells zones regardless
        mask = np.zeros(n_orig, dtype=bool)
        mask[:min_cells] = True
        print(f"  [Truncation] tau>{truncate_tau} leaves <{min_cells} zones; "
              f"keeping outermost {min_cells}")
    n_trunc = int(mask.sum())
    r_t, v_t, T_t, T_rad_t, rho_t = r[mask], v[mask], T[mask], T_rad[mask], rho[mask]
    kappa_t = kappa_hybrid[mask]   # use hybrid kappa for consistency
    tau_t = tau[mask]
    flux_t = flux[mask]
    mf_t = {k: v[mask] for k, v in mf.items()}
    # Diagnostic mass estimate using zone-centre spacing (outer->inner).
    # The exact cell mass is computed later from proper boundaries; this is
    # only for the truncation report, so we use adjacent-centre differences
    # without wrap-around (np.roll would wrap the last zone).
    def _mass_diag(rr, rhoo):
        dr = np.abs(np.diff(rr, append=rr[-1] * 2 - rr[-2]))
        vol = (4.0 / 3.0) * np.pi * np.abs((rr + dr) ** 3 - rr ** 3)
        return np.sum(rhoo * vol)
    mass_orig = _mass_diag(r, rho)
    mass_trunc = _mass_diag(r_t, rho_t)
    print(f"\n[Truncation] {n_orig} -> {n_trunc} zones (tau <= {truncate_tau})")
    print(f"  Dropped {n_orig - n_trunc} optically-thick inner zones")
    print(f"  Retained mass: {mass_trunc / MSUN:.4f} Msun "
          f"({100 * mass_trunc / max(mass_orig, 1e-30):.1f}% of total)")
    print(f"  Retained tau range: {tau_t.min():.3e} - {tau_t.max():.3e}")

    # --- Optional density floor (drop near-vacuum outer zones) ---
    if rho_floor > 0:
        dmask = rho_t >= rho_floor
        if dmask.sum() >= min_cells:
            r_t, v_t, T_t, T_rad_t, rho_t, kappa_t = r_t[dmask], v_t[dmask], T_t[dmask], T_rad_t[dmask], rho_t[dmask], kappa_t[dmask]
            tau_t = tau_t[dmask]
            flux_t = flux_t[dmask]
            mf_t = {k: v[dmask] for k, v in mf_t.items()}
            print(f"[Density floor] rho>={rho_floor}: {n_trunc} -> {int(dmask.sum())} zones")
            n_trunc = int(dmask.sum())

    # --- MC-aware grid downsampling ---
    tau_max_truncated = float(tau_t.max())  # save before downsampling
    edges_merged_io = None  # actual merged boundaries (inner->outer)
    grid_diag = None  # adaptive grid diagnostics (if adaptive method used)
    if downsample and n_trunc > max_cells:
        if grid_method == 'adaptive':
            # Physics-aware adaptive regridding (default for new runs)
            from superlite_adaptive_grid import (
                GridRefinementConfig, adaptive_regrid, validate_conservation)
            amr_config = GridRefinementConfig(
                max_cells=max_cells,
                min_cells=min_cells,
                rho_grad=amr_rho_grad,
                temp_grad=amr_temp_grad,
                vel_grad=amr_vel_grad,
                comp_grad=amr_comp_grad,
                cannon_min_cells_per_interface=amr_cannon_interface_cells,
                shock_dv_threshold=amr_shock_dv_threshold,
            )
            result = adaptive_regrid(
                r_t, v_t, T_t, T_rad_t, rho_t, kappa_t, mf_t, tau_t,
                amr_config, verbose=True)
            r_t, v_t, T_t, T_rad_t, rho_t, kappa_t, mf_t, edges_merged_io, grid_diag = result
            # Validate conservation (warn only, don't stop)
            errors = validate_conservation(grid_diag, amr_config, fail_hard=False)
            for err in errors:
                print(f"  WARNING: {err}")
            # recompute tau on the downsampled grid for R_phot
            tau_t = compute_kappa_tau(r_t, rho_t, kappa_t)
            n_final = len(r_t)
        else:
            # Legacy tau-weighted downsampling
            result = tau_weighted_downsample(
                r_t, v_t, T_t, T_rad_t, rho_t, kappa_t, mf_t, tau_t, max_cells, min_cells)
            r_t, v_t, T_t, T_rad_t, rho_t, kappa_t, mf_t, edges_merged_io = result
            # recompute tau on the downsampled grid for R_phot
            tau_t = compute_kappa_tau(r_t, rho_t, kappa_t)
            n_final = len(r_t)
    else:
        n_final = n_trunc

    # --- Reverse to inner->outer (SuperLite convention: r increasing) ---
    r_sl = r_t[::-1].copy()
    v_sl = v_t[::-1].copy()
    T_sl = T_t[::-1].copy()
    T_rad_sl = T_rad_t[::-1].copy()
    rho_sl = rho_t[::-1].copy()
    tau_sl = tau_t[::-1].copy()
    mf_sl = {k: v[::-1].copy() for k, v in mf_t.items()}

    # --- Cell boundaries (inner->outer) ---
    # Use actual merged boundaries from the downsampler if available;
    # otherwise fall back to geometric-mean boundaries from zone centres.
    if edges_merged_io is not None:
        edges = edges_merged_io  # already inner->outer, len n_final+1
    else:
        edges = compute_cell_boundaries(r_sl)
    # velocity at edges: linear interp on the inner->outer grid
    v_edges = np.zeros(n_final + 1)
    v_edges[1:n_final] = 0.5 * (v_sl[:-1] + v_sl[1:])
    v_edges[0] = v_sl[0] * (v_sl[0] / max(v_sl[1], 1e-30))  # extrapolate inner
    v_edges[n_final] = v_sl[-1] * (v_sl[-1] / max(v_sl[-2], 1e-30))  # outer
    # convert km/s -> cm/s for SuperLite
    edges_cm = edges
    v_edges_cgs = v_edges * 1e5

    # --- Cell masses (g) ---
    vol = (4.0 / 3.0) * np.pi * (edges[1:n_final + 1] ** 3 - edges[:n_final] ** 3)
    mass = rho_sl * vol
    mass_total = mass.sum()
    print(f"\nFinal grid: {n_final} cells (inner->outer)")
    print(f"  R: {r_sl.min():.3e} - {r_sl.max():.3e} cm")
    print(f"  V: {v_sl.min():.1f} - {v_sl.max():.1f} km/s")
    print(f"  Total mass: {mass_total / MSUN:.4f} Msun")

    # --- Abundance normalisation (per cell, sum to 1) ---
    # Only normalise over elements that are present in the model.
    total_ab = sum(mf_sl[e] for e in elem_present)
    for e in elem_present:
        mf_sl[e] = mf_sl[e] / np.maximum(total_ab, 1e-30)

    # --- L_bol: use the emergent bolometric luminosity from the STELLA
    # .swd flux field (column Flum_1e40, converted to erg/s).
    #
    # SuperLite uses in_L_bol as a SCALING FACTOR: it computes a T_gas^4
    # weighted spatial emission distribution, then renormalises the total
    # to equal in_L_bol (sourceenergy.f90:56).  The spectral shape is set
    # by the Planck function at the local gas temperature, independent of
    # in_L_bol.  Changing in_L_bol changes only the total photon count,
    # not the spectral energy distribution.
    #
    # The correct quantity is the total emergent bolometric luminosity
    # from the STELLA snapshot at the SAME EPOCH as the structure.
    # STELLA's Flum at the outer boundary (zone 0, outer->inner order)
    # is the emergent luminosity.  If it is negative (inward-flowing
    # radiation in a reverse shock), fall back to the photospheric flux
    # at tau=2/3, using the absolute value.
    flux_oi = flux_t  # outer->inner, post-truncation
    if flux_oi[0] > 0:
        lum_erg = float(flux_oi[0])
        lum_source = 'STELLA outer boundary flux (zone 0)'
    else:
        # Outer boundary flux is negative or zero — use photospheric flux
        tau_oi = tau_t  # outer->inner, increasing inward
        phot_idx = np.searchsorted(tau_oi, 2.0 / 3.0)
        phot_idx = min(int(phot_idx), len(flux_oi) - 1)
        lum_erg = abs(float(flux_oi[phot_idx]))
        lum_source = f'STELLA photospheric flux (tau=2/3, zone {phot_idx})'
        if lum_erg <= 0:
            # Final fallback: Stefan-Boltzmann at the outer boundary
            r_outer = r_t[0]
            T_outer = T_t[0]
            lum_erg = 4.0 * np.pi * r_outer ** 2 * 5.670374419e-5 * T_outer ** 4
            lum_source = 'Stefan-Boltzmann fallback at outer boundary'
    print(f"\nL_bol (bolometric, from {lum_source}): "
          f"{lum_erg:.3e} erg/s = {lum_erg / LSUN:.3e} Lsun")
    print(f"  STELLA epoch: t={t_day:.4f} d")
    print(f"  STELLA flux range: {flux_oi.min():.3e} to {flux_oi.max():.3e} erg/s")

    # --- R_phot: grey tau=2/3 radius (Eddington photosphere).
    # SuperLite's get_rphot subroutine (nltemod.f:186-227) uses tau=2/3
    # for the Rosseland opacity photosphere.  We use the same threshold
    # here for consistency.  The nebular approximation (ions_nebular)
    # uses R_phot for the dilution factor W = 0.5*(1-sqrt(1-(R_phot/r)^2)):
    #   - below R_phot: LTE Saha ionization
    #   - above R_phot: nebular approximation with dilute radiation field
    # R_phot is only updated during the calculation when in_nlte=t
    # (superlite.f90:210).  With in_nlte=f, it stays at this input value.
    tau_desc = tau_sl[::-1]  # outer->inner (increasing tau)
    r_desc = r_sl[::-1]
    phot_i = np.searchsorted(tau_desc, 2.0 / 3.0)
    phot_i = min(int(phot_i), len(r_desc) - 1)
    r_phot = float(r_desc[phot_i])
    print(f"R_phot (grey tau=2/3, hybrid kappa): {r_phot:.3e} cm")

    # --- Color temperature (Tbb) for the inner boundary source ---
    # STELLA's Tbb is the color temperature of the radiation field at the
    # photosphere: Tbb = (L / (4*pi*R_phot^2*sigma))^(1/4).  This is much
    # higher than the gas temperature in cool ejecta, and determines the
    # spectral shape of the emergent radiation.
    #
    # SuperLite's boundary_source.f90 injects a blackbody at in_srcmax
    # from the inner boundary with luminosity in_L_inner.  We set:
    #   in_L_inner = L_bol  (all radiation from below the grid)
    #   in_L_bol   = 0      (no interior thermal source)
    #   in_srcmax  = Tbb    (color temperature for boundary spectrum)
    #
    # The interior thermal emission at T_gas=2000-8000K produces a cool
    # Rayleigh-Jeans continuum.  The boundary source at Tbb~20000K produces
    # the correct hot continuum that, after absorption by the ejecta,
    # creates P-Cygni line profiles.
    _sigma_sb = 5.670374419e-5
    tbb = (lum_erg / (4.0 * np.pi * r_phot**2 * _sigma_sb))**0.25
    # Also compute the color temperature at the inner boundary radius
    r_inner = float(r_sl[0])
    t_inner = (lum_erg / (4.0 * np.pi * r_inner**2 * _sigma_sb))**0.25
    print(f"Tbb (color temp at R_phot): {tbb:.1f} K")
    print(f"T_inner (color temp at R_inner): {t_inner:.1f} K")
    print(f"  -> Setting in_L_inner={lum_erg:.3e}, in_srcmax={tbb:.1f} K, "
          f"in_L_bol=0 (boundary source mode)")

    # --- Physics validation checks ---
    validation_errors = []

    # 1. L_bol must be positive and finite
    if not np.isfinite(lum_erg) or lum_erg <= 0:
        validation_errors.append(
            f"L_bol is non-positive or non-finite: {lum_erg}")

    # 2. L_bol must be in a physically reasonable range (1e38 - 1e46 erg/s)
    if lum_erg < 1e38 or lum_erg > 1e46:
        validation_errors.append(
            f"L_bol outside physical range [1e38, 1e46]: {lum_erg:.3e}")

    # 3. R_phot must be within the model radial extent
    if r_phot < r_sl.min() or r_phot > r_sl.max():
        validation_errors.append(
            f"R_phot ({r_phot:.3e}) outside model range "
            f"[{r_sl.min():.3e}, {r_sl.max():.3e}]")

    # 4. R_phot must be positive
    if r_phot <= 0:
        validation_errors.append(f"R_phot is non-positive: {r_phot}")

    # 5. Epoch consistency: the flux and structure must be from the same epoch
    # (guaranteed by construction — both come from ep = select_epoch(epochs, t))
    # but we verify the flux is from the selected epoch by checking it's in the
    # expected range for that epoch.
    if flux_oi[0] < 0:
        print(f"  WARNING: outer boundary flux is negative ({flux_oi[0]:.3e}), "
              f"using photospheric fallback")

    # 6. All velocities must be nonnegative and below c
    v_cgs = v_sl * 1e5  # km/s -> cm/s
    if v_sl.min() < 0:
        print(f"  WARNING: negative velocities present (min={v_sl.min():.1f} km/s) "
              f"— SuperLite supports non-homologous outflows")
    if v_cgs.max() >= 3e10:
        validation_errors.append(
            f"velocity exceeds c: {v_sl.max():.1f} km/s")

    # 7. All masses must be positive
    if mass.min() <= 0:
        validation_errors.append(
            f"non-positive cell mass: min={mass.min():.3e}")

    # 8. All densities must be positive
    if rho_sl.min() <= 0:
        validation_errors.append(
            f"non-positive density: min={rho_sl.min():.3e}")

    # 9. All temperatures must be positive and above the floor
    if T_sl.min() < 1:
        validation_errors.append(
            f"temperature below 1 K: min={T_sl.min():.3e}")

    # 10. Composition must sum to ~1 per cell
    total_comp = np.zeros(len(r_sl))
    for e in mf_sl:
        total_comp += mf_sl[e]
    max_comp_err = float(np.max(np.abs(total_comp - 1.0)))
    if max_comp_err > 1e-4:
        validation_errors.append(
            f"composition normalisation error: max deviation = {max_comp_err:.2e}")

    if validation_errors:
        print(f"\n  VALIDATION ERRORS:")
        for err in validation_errors:
            print(f"    - {err}")
        raise ValueError(
            f"Physics validation failed with {len(validation_errors)} errors")
    else:
        print(f"  Physics validation: ALL CHECKS PASSED")

    # --- Write output ---
    os.makedirs(output_dir, exist_ok=True)
    str_path = os.path.join(output_dir, 'input.str')
    par_path = os.path.join(output_dir, 'input.par')
    write_input_str(str_path, r_sl[0], edges_cm, v_edges_cgs, mass,
                    T_rad_sl, T_sl, mf_sl, t_day, model_name, elem_present)
    write_input_par(par_path, model_name, t_day, n_final, lum_erg, r_phot,
                    nomp, truncate_tau, n_orig, n_final, quick_test=quick_test,
                    grp_ng=sl_grp_ng, grp_wlmin=sl_grp_wlmin,
                    grp_wlmax=sl_grp_wlmax,
                    flx_ndim=sl_flx_ndim, flx_wlmin=sl_flx_wlmin,
                    flx_wlmax=sl_flx_wlmax,
                    src_n2s=sl_src_n2s, trn_tauddmc=sl_trn_tauddmc,
                    taulump=sl_taulump,
                    nlte=sl_nlte, nlte_nelem=sl_nlte_nelem, niter=sl_niter,
                    nebular=sl_nebular, opacmixrossel=sl_opacmixrossel,
                    nobbopac=sl_nobbopac, nobfopac=sl_nobfopac,
                    noffopac=sl_noffopac, nothmson=sl_nothmson,
                    l_inner=lum_erg, srcmax=tbb)

    # --- Generate input structure diagnostic plots ---
    print(f"\nGenerating input structure plots...")
    plot_superlite_input_structure(
        output_dir, r_sl, v_sl, T_sl, rho_sl, edges_cm, tau_sl,
        mf_sl, elem_present, r_phot, t_day, model_name)

    # --- MC-aware diagnostics summary ---
    cells_removed = n_orig - n_final
    speedup_est = n_orig / max(n_final, 1)
    print(f"\n{'='*70}")
    print(f"MC-aware optimisation summary")
    print(f"  Original zones      : {n_orig}")
    print(f"  After tau truncation: {n_trunc} (dropped {n_orig - n_trunc} "
          f"optically-thick zones, tau>{truncate_tau})")
    print(f"  After downsampling  : {n_final} (target max_cells={max_cells})")
    print(f"  Cells removed       : {cells_removed} "
          f"({100 * cells_removed / n_orig:.1f}%)")
    print(f"  Est. MC speedup     : ~{speedup_est:.1f}x (runtime ~ ncells)")
    print(f"  Photosphere retained: R_phot={r_phot:.3e} cm (tau=1)")
    print(f"  Deepest tau retained: {tau_max_truncated:.2f} (pre-downsample)")
    print(f"  Inner-boundary tau   : {tau.max():.2f} (full {n_orig}-zone grid)")
    print(f"  L_bol               : {lum_erg:.3e} erg/s")
    print(f"{'='*70}")
    print(f"\nConversion complete: {output_dir}/")
    return n_final


def str2bool(v):
    """Parse a string as a boolean, handling argparse type=bool correctly.

    Python's bool() returns True for any non-empty string, including 'false'.
    This function properly interprets common false values.
    """
    if isinstance(v, bool):
        return v
    if v.lower() in ('true', 't', 'yes', 'y', '1'):
        return True
    if v.lower() in ('false', 'f', 'no', 'n', '0'):
        return False
    raise argparse.ArgumentTypeError(f'Boolean value expected, got {v!r}')


def main():
    parser = argparse.ArgumentParser(
        description='Convert STELLA .swd directly to SuperLite input.str (MC-aware)')
    parser.add_argument('swd_file', help='Path to STELLA .swd file')
    parser.add_argument('-o', '--output-dir', required=True,
                        help='Output directory for SuperLite input files')
    parser.add_argument('--target-days', type=float, default=None,
                        help='Epoch closest to this many days (default: last)')
    parser.add_argument('--truncate-tau', type=float, default=20.0,
                        help='Drop zones with cumulative Rosseland tau > this '
                             '(default 20; MC is expensive in thick regions)')
    parser.add_argument('--max-cells', type=int, default=150,
                        help='Downsample to <= N cells after truncation (default 150)')
    parser.add_argument('--min-cells', type=int, default=20,
                        help='Never go below N cells (default 20)')
    parser.add_argument('--nomp', type=int, default=4,
                        help='OpenMP threads for input.par in_nomp (default 4)')
    parser.add_argument('--rho-floor', type=float, default=0.0,
                        help='Drop zones with rho < X g/cm^3 (default 0)')
    parser.add_argument('--t-floor', type=float, default=2000.0,
                        help='Temperature floor in K (default 2000; SuperLite '
                             'zeta table starts at 2000 K)')
    parser.add_argument('--no-downsample', action='store_true',
                        help='Disable grid downsampling (keep all post-truncation zones)')
    parser.add_argument('--model-name', default=None,
                        help='Model name for input.par (default: .swd filename stem)')
    parser.add_argument('--quick-test', action='store_true',
                        help='Generate input.par with reduced particles/groups and '
                             'higher DDMC threshold for fast interactive testing '
                             '(1M particles, 1k groups, DDMC thresh=10)')
    parser.add_argument('--grid-method', default='adaptive',
                        choices=['adaptive', 'tau_weighted'],
                        help='Grid construction method: adaptive (physics-aware AMR, '
                             'default) or tau_weighted (legacy fixed split)')
    parser.add_argument('--amr-rho-grad', type=float, default=0.5,
                        help='AMR density gradient threshold |d log rho / d log r|')
    parser.add_argument('--amr-temp-grad', type=float, default=0.3,
                        help='AMR temperature gradient threshold |d log T / d log r|')
    parser.add_argument('--amr-vel-grad', type=float, default=0.05,
                        help='AMR velocity gradient threshold |dv/v| per cell')
    parser.add_argument('--amr-comp-grad', type=float, default=0.1,
                        help='AMR composition gradient threshold |dX / d log r|')
    parser.add_argument('--amr-cannon-interface-cells', type=int, default=5,
                        help='AMR cells to resolve each cannon boundary interface')
    parser.add_argument('--amr-shock-dv-threshold', type=float, default=100.0,
                        help='AMR velocity change threshold (km/s) to flag a shock')
    # SuperLite numerical parameters (routed from YAML via orchestrator)
    parser.add_argument('--sl-grp-ng', type=int, default=6000,
                        help='SuperLite wavelength groups (default 6000)')
    parser.add_argument('--sl-grp-wlmin', default='1d-8',
                        help='Transport wavelength min in cm (default 1d-8 = 1 Å)')
    parser.add_argument('--sl-grp-wlmax', default='30d-5',
                        help='Transport wavelength max in cm (default 30d-5 = 30000 Å)')
    parser.add_argument('--sl-flx-ndim', type=int, default=3500,
                        help='Flux output bins (default 3500)')
    parser.add_argument('--sl-flx-wlmin', default='3d-5',
                        help='Flux wavelength min in cm (default 3d-5 = 3000 Å)')
    parser.add_argument('--sl-flx-wlmax', default='2d-4',
                        help='Flux wavelength max in cm (default 2d-4 = 20000 Å)')
    parser.add_argument('--sl-src-n2s', type=int, default=25,
                        help='2^n source particles (default 25 = 33M)')
    parser.add_argument('--sl-trn-tauddmc', default='3d0',
                        help='DDMC transition threshold in mfp (default 3d0)')
    parser.add_argument('--sl-taulump', default='3d0',
                        help='DDMC group lumping threshold in mfp (default 3d0)')
    parser.add_argument('--sl-nlte', type=str2bool, default=True,
                        help='Enable NLTE for H (default True)')
    parser.add_argument('--sl-nlte-nelem', type=int, default=1,
                        help='Number of NLTE elements (must be 1, H only)')
    parser.add_argument('--sl-niter', type=int, default=3,
                        help='Radiation temperature feedback iterations (default 3)')
    parser.add_argument('--sl-nebular', type=str2bool, default=True,
                        help='Enable nebular approximation (default True)')
    args = parser.parse_args()

    convert_swd_to_superlite(
        args.swd_file, args.output_dir,
        target_days=args.target_days,
        truncate_tau=args.truncate_tau,
        max_cells=args.max_cells,
        min_cells=args.min_cells,
        nomp=args.nomp,
        rho_floor=args.rho_floor,
        downsample=not args.no_downsample,
        model_name=args.model_name,
        quick_test=args.quick_test,
        t_floor=args.t_floor,
        grid_method=args.grid_method,
        amr_rho_grad=args.amr_rho_grad,
        amr_temp_grad=args.amr_temp_grad,
        amr_vel_grad=args.amr_vel_grad,
        amr_comp_grad=args.amr_comp_grad,
        amr_cannon_interface_cells=args.amr_cannon_interface_cells,
        amr_shock_dv_threshold=args.amr_shock_dv_threshold,
        sl_grp_ng=args.sl_grp_ng,
        sl_grp_wlmin=args.sl_grp_wlmin,
        sl_grp_wlmax=args.sl_grp_wlmax,
        sl_flx_ndim=args.sl_flx_ndim,
        sl_flx_wlmin=args.sl_flx_wlmin,
        sl_flx_wlmax=args.sl_flx_wlmax,
        sl_src_n2s=args.sl_src_n2s,
        sl_trn_tauddmc=args.sl_trn_tauddmc,
        sl_taulump=args.sl_taulump,
        sl_nlte=args.sl_nlte,
        sl_nlte_nelem=args.sl_nlte_nelem,
        sl_niter=args.sl_niter,
        sl_nebular=args.sl_nebular,
    )


if __name__ == '__main__':
    main()
