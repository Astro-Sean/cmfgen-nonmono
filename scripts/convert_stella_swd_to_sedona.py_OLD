#!/usr/bin/env python3
"""Convert a STELLA *.swd snapshot directly to SEDONA .mod format.

This is a one-step conversion (STELLA .swd -> SEDONA .mod) that replaces the
previous two-step pipeline:

    .swd -> CMFGEN SN_HYDRO_DATA  (convert_stella_swd_to_cmfgen.py)
    SN_HYDRO_DATA -> .mod         (external convert_cmfgen_to_sedona_v2.py)

with a single, self-contained converter that reads the .swd directly.

Advantages over the two-step pipeline
-------------------------------------
1. **No CMFGEN dependency** — SEDONA no longer requires a CMFGEN SN_HYDRO_DATA
   to have been generated.  This lets SEDONA run even when CMFGEN is skipped.
2. **Native STELLA structure** — the full STELLA radial grid is preserved
   (after tau-truncation), rather than the CMFGEN AMR-regridded snapshot.
3. **Full 15-element composition** — all STELLA elements are passed through
   individually (SEDONA's cmfgen_levelcap100.hdf5 supports all 15).
4. **Native Rosseland kappa** — STELLA's NLTE Rosseland opacity is used for
   the optical-depth truncation, with an electron-scattering floor for
   robustness in cool neutral regions that SEDONA may re-ionize.

SEDONA .mod format (1D_sphere standard)
---------------------------------------
    Line 1: "1D_sphere standard"
    Line 2: nzon  rmin  texp  n_elems
    Line 3: space-separated atomic numbers (Z) of the species columns
    Line 4+: r(cm) v(cm/s) rho(g/cm^3) T(K) mass_frac_1 ... mass_frac_n
    Zones are sorted by ascending radius (inner -> outer), which is SEDONA's
    required ordering (opposite of CMFGEN/STELLA .swd outer->inner).

Usage
-----
    python convert_stella_swd_to_sedona.py SWD_FILE OUTPUT.mod [OPTIONS]
      --target-days T    Pick the epoch closest to T days (default: last epoch)
      --time-days T      Time since explosion (days) written to .mod header
                         (default: the selected epoch time)
      --tau-max X        Truncate inner zones with cumulative tau > X
                         (default: 100).  Removes the optically-thick diffusion
                         region where MC photons would random-walk.
      --min-zones N      Always keep at least N zones after truncation
                         (default: 20).
      --nzon-target N    Downsample to ~N zones (0 = no downsampling, default: 0).
      --t-floor K        Temperature floor in K (default: 2000).
      --rho-floor G      Drop outer zones with rho < G g/cm^3 (default: 0 = off).
      --species Z,Z,...  Comma-separated atomic numbers to include
                         (default: all STELLA elements present in the model).
      --homologous       Write header for hydro_module="homologous" (Sobolev
                         line expansion uses 1/t_exp).  Otherwise writes a
                         static-snapshot header (hydro_module="none").
      --use-T-rad        Use radiation temperature instead of gas temperature.
      --model-name NAME  For diagnostics (default: .swd filename stem).

The .swd arrays from read_swd() are outer->inner (decreasing r).  This script
reverses them to inner->outer (ascending r) for the SEDONA .mod output.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np

# Reuse the battle-tested .swd parser and composition normaliser from the
# CMFGEN conversion pipeline rather than duplicating the fixed-width parsing.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from convert_stella_swd_to_cmfgen import (
    read_swd,
    normalize_swd_composition,
    STELLA_ELEM_NAMES,
    MSUN,
    LSUN,
)

# STELLA Natom=15 element order -> atomic number (Z).
# Matches STELLA_ELEM_NAMES from convert_stella_swd_to_cmfgen.
STELLA_ELEM_Z = {
    'H': 1, 'He': 2, 'C': 6, 'N': 7, 'O': 8, 'Ne': 10, 'Na': 11,
    'Mg': 12, 'Al': 13, 'Si': 14, 'S': 16, 'Ar': 18, 'Ca': 20,
    'Fe': 26, 'Ni': 28,
}
# Canonical order (all 15 STELLA elements), used for species column ordering.
STELLA_ELEM_ORDER = ['H', 'He', 'C', 'N', 'O', 'Ne', 'Na', 'Mg', 'Al',
                     'Si', 'S', 'Ar', 'Ca', 'Fe', 'Ni']

SIGMA_SB = 5.670374419e-5  # erg/cm^2/s/K^4


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


def detect_elements(mf: dict, threshold: float = 1e-10) -> list[str]:
    """Detect which elements have non-zero mass fractions in the model.

    SEDONA loads atomic data for every species column in the .mod file.
    Including zero-abundance elements wastes startup time and memory, and
    can cause singular NLTE rate matrices in low-density zones.  Only
    include elements with max mass fraction above ``threshold``.

    Args:
        mf: dict of element name -> mass fraction array (any orientation;
            we take the max).
        threshold: minimum max mass fraction to consider an element present.

    Returns:
        list of element names (in STELLA_ELEM_ORDER order) that are present.
    """
    present = []
    for elem in STELLA_ELEM_ORDER:
        if elem in mf:
            mx = float(np.max(mf[elem]))
            if mx > threshold:
                present.append(elem)
    return present


def compute_hybrid_tau(r_oi: np.ndarray, rho_oi: np.ndarray,
                       kappa_swd: np.ndarray, T_oi: np.ndarray,
                       mf: dict, verbose: bool = True) -> tuple[np.ndarray, np.ndarray]:
    """Compute a robust cumulative optical depth (outer -> inner).

    A hybrid opacity is used, taking the element-wise maximum of:

      1. STELLA native Rosseland kappa (WRKX) — accurate where STELLA's NLTE
         ionization is reliable (hot, ionized regions).
      2. Electron-scattering floor — kappa_es = 0.2*(1+X) cm^2/g assuming
         full ionization.  This is the maximum ES opacity SEDONA could see,
         and ensures we don't underestimate tau in cool neutral regions that
         SEDONA may re-ionize (where STELLA's WRKX drops to ~0).

    The cumulative tau is integrated from the outer boundary inward
    (outer->inner, decreasing r) using the trapezoidal rule on zone-centre
    spacings.

    Args:
        r_oi: radius array, outer->inner (decreasing r), cm.
        rho_oi: density array, outer->inner, g/cm^3.
        kappa_swd: STELLA Rosseland kappa, outer->inner, cm^2/g.
        T_oi: temperature array, outer->inner, K.
        mf: dict of element mass fractions, outer->inner.
        verbose: print diagnostics.

    Returns:
        tau: cumulative optical depth array, outer->inner.
        kappa_hybrid: the hybrid opacity array used, outer->inner.
    """
    n = len(r_oi)

    # Component 1: STELLA native kappa (clamp non-physical negatives to 0).
    kappa_native = np.maximum(kappa_swd, 0.0)

    # Component 2: electron-scattering floor (full ionization).
    # kappa_es = sigma_T * n_e / rho = 0.2 * (1 + X) cm^2/g for H/He mix.
    X = mf.get('H', np.zeros(n))
    Y = mf.get('He', np.zeros(n))
    kappa_es = 0.2 * (1.0 + X)

    kappa_hybrid = np.maximum(kappa_native, kappa_es)

    # Cumulative tau from outer boundary inward (outer->inner, decreasing r).
    # Use zone-centre spacings; append an extrapolated spacing for the last
    # (innermost) zone so every zone contributes.
    dr = np.abs(np.diff(r_oi, append=r_oi[-1] * 2 - r_oi[-2]))
    dtau = kappa_hybrid * rho_oi * dr
    tau = np.cumsum(dtau)  # tau[0] = outermost, increasing inward

    if verbose:
        print(f"  Hybrid kappa: {kappa_hybrid.min():.3e} - {kappa_hybrid.max():.3e} cm^2/g")
        print(f"  Cumulative tau (outer->inner): {tau[0]:.3e} (outer) -> {tau[-1]:.3e} (inner)")

    return tau, kappa_hybrid


def downsample_grid(r_io, v_io, T_io, rho_io, mf_io, tau_io, nzon_target):
    """Downsample to ~nzon_target zones using uniform sampling in log(tau).

    Preserves the inner and outer boundaries.  All arrays are inner->outer
    (ascending r).  Mass fractions are averaged density-weighted within each
    merged cell; r, v, T are density-weighted averages; rho is volume-weighted.

    Args:
        *_io: inner->outer (ascending r) arrays.
        tau_io: cumulative tau, inner->outer (ascending, so tau increases
            outward from the inner boundary).
        nzon_target: target number of zones.

    Returns:
        Downsampled arrays (r, v, T, rho, mf dict, tau), inner->outer.
    """
    n = len(r_io)
    if nzon_target <= 0 or n <= nzon_target:
        return r_io, v_io, T_io, rho_io, mf_io, tau_io

    # Sample uniformly in log(tau) space, always keeping the first and last.
    # tau_io is inner->outer and DESCENDING (high tau at inner boundary, low
    # at outer).  np.interp requires ascending x-coordinates, so sort first
    # and map target log-tau values back to original indices via interpolation.
    log_tau = np.log10(np.maximum(tau_io, 1e-30))
    asc_order = np.argsort(log_tau)          # indices sorting log_tau ascending
    log_tau_asc = log_tau[asc_order]         # ascending log-tau
    # Target log-tau positions spanning the full range (ascending).
    target_log_tau = np.linspace(log_tau_asc[0], log_tau_asc[-1], nzon_target)
    # Interpolate: for each target log-tau, the original index where it occurs.
    idx_float = np.interp(target_log_tau, log_tau_asc, asc_order)
    idx = np.unique(np.round(idx_float).astype(int))
    idx = np.clip(idx, 0, n - 1)
    # Ensure endpoints are included.
    idx = np.unique(np.concatenate(([0], idx, [n - 1])))

    # Build merged cells: each selected index is a representative zone.
    # Cell boundaries are midpoints between selected radii.
    r_sel = r_io[idx]
    n_sel = len(idx)
    edges = np.zeros(n_sel + 1)
    edges[1:n_sel] = np.sqrt(r_sel[:-1] * r_sel[1:])  # geometric mean
    edges[0] = r_sel[0] ** 2 / edges[1] if edges[1] > 0 else r_sel[0] * 0.5
    edges[n_sel] = r_sel[-1] ** 2 / edges[n_sel - 1] if edges[n_sel - 1] > 0 else r_sel[-1] * 1.5

    r_out = np.zeros(n_sel)
    v_out = np.zeros(n_sel)
    T_out = np.zeros(n_sel)
    rho_out = np.zeros(n_sel)
    tau_out = np.zeros(n_sel)
    mf_out = {k: np.zeros(n_sel) for k in mf_io}

    for j in range(n_sel):
        # Zones in this merged cell: from edges[j] to edges[j+1].
        mask = (r_io >= edges[j]) & (r_io <= edges[j + 1])
        if not np.any(mask):
            # Fall back to the representative index.
            mask = np.array([False] * n)
            mask[idx[j]] = True
        w = rho_io[mask]  # density weighting
        wsum = np.sum(w)
        if wsum <= 0:
            w = np.ones(np.sum(mask))
            wsum = float(len(w))
        # Volume of each constituent zone (approximate from centre spacing).
        r_cell = r_io[mask]
        n_cell = len(r_cell)
        if n_cell >= 2:
            dr = np.abs(np.diff(r_cell))
            dr = np.append(dr, dr[-1])  # extrapolate last zone spacing
        else:
            dr = np.array([r_cell[0] * 0.01])  # single-zone cell
        vol = np.abs((r_cell + dr) ** 3 - r_cell ** 3)
        vol = np.maximum(vol, 1e-30)
        vw = vol * rho_io[mask]
        vwsum = np.sum(vw)
        if vwsum <= 0:
            vw = np.ones(np.sum(mask))
            vwsum = float(len(vw))
        r_out[j] = np.sum(r_cell * vw) / vwsum
        v_out[j] = np.sum(v_io[mask] * vw) / vwsum
        T_out[j] = np.sum(T_io[mask] * vw) / vwsum
        rho_out[j] = np.sum(rho_io[mask] * vol) / np.sum(vol)
        tau_out[j] = tau_io[idx[j]]
        for k in mf_io:
            mf_out[k][j] = np.sum(mf_io[k][mask] * vw) / vwsum

    # Renormalise mass fractions to sum to 1 per zone (guard against drift).
    for j in range(n_sel):
        total = sum(mf_out[k][j] for k in mf_out)
        if total > 0:
            for k in mf_out:
                mf_out[k][j] /= total

    print(f"  Downsampled {n} -> {n_sel} zones (uniform in log-tau)")
    return r_out, v_out, T_out, rho_out, mf_out, tau_out


def convert_swd_to_sedona(swd_path: str, output_file: str,
                          target_days: float | None = None,
                          time_days: float | None = None,
                          tau_max: float = 100.0,
                          min_zones: int = 20,
                          nzon_target: int = 0,
                          t_floor: float = 2000.0,
                          rho_floor: float = 0.0,
                          species_list: list[int] | None = None,
                          homologous: bool = False,
                          use_T_rad: bool = False,
                          model_name: str | None = None):
    """Main conversion: STELLA .swd -> SEDONA .mod (one step).

    See module docstring for the .mod format and the conversion strategy.
    """
    if model_name is None:
        model_name = Path(swd_path).stem

    print(f"\n{'='*70}")
    print(f"STELLA .swd -> SEDONA .mod (direct conversion)")
    print(f"  source : {swd_path}")
    print(f"  output : {output_file}")
    print(f"  tau_max: {tau_max}  nzon_target: {nzon_target}  "
          f"t_floor: {t_floor}  homologous: {homologous}")
    print(f"{'='*70}")

    # --- Parse .swd (reuses convert_stella_swd_to_cmfgen.read_swd) ---
    epochs = read_swd(swd_path)
    t_day, ep = select_epoch(epochs, target_days)

    # .swd arrays are outer->inner (decreasing r).
    r = ep['r']            # cm
    v_kms = ep['v']        # km/s
    T_gas = ep['T']        # gas kinetic temperature, K
    T_rad = ep['T_rad']    # radiation temperature, K (= T_gas in optically-thick inner zones)
    rho = ep['rho']        # g/cm^3
    flux = ep['flux']      # erg/s (per-zone luminosity)
    kappa = ep['kappa']    # Rosseland opacity, cm^2/g (native from STELLA NLTE)
    comp = ep['composition']  # (nzon, 15) Y_k*A_k

    n_orig = len(r)
    print(f"\nEpoch has {n_orig} zones (outer->inner)")
    print(f"  R: {r.min():.3e} - {r.max():.3e} cm")
    print(f"  V: {v_kms.min():.1f} - {v_kms.max():.1f} km/s")
    print(f"  T_gas: {T_gas.min():.3e} - {T_gas.max():.3e} K")
    print(f"  T_rad: {T_rad.min():.3e} - {T_rad.max():.3e} K")
    print(f"  rho: {rho.min():.3e} - {rho.max():.3e} g/cm^3")
    if kappa is not None:
        print(f"  kappa_Ross: {kappa.min():.3e} - {kappa.max():.3e} cm^2/g (native STELLA)")
    else:
        print(f"  kappa_Ross: NOT in .swd, falling back to grey 0.34 cm^2/g")
        kappa = np.full(n_orig, 0.34)

    # --- Choose temperature field ---
    if use_T_rad:
        T = T_rad.copy()
        print(f"  Using RADIATION temperature (T_rad) for SEDONA .mod")
    else:
        T = T_gas.copy()
        print(f"  Using GAS temperature (T_gas) for SEDONA .mod")

    # --- Temperature floor ---
    # SEDONA's LTE Saha-Boltzmann populations and opacity tables are
    # undefined/unreliable below ~2000 K.  Clamp to the floor.
    if t_floor > 0:
        n_clamped = int(np.sum(T < t_floor))
        if n_clamped > 0:
            print(f"  [T-floor] Clamping {n_clamped} zones from "
                  f"T_min={T.min():.1f} K to t_floor={t_floor:.1f} K")
            T = np.maximum(T, t_floor)

    # --- Normalise composition to true mass fractions ---
    mf = normalize_swd_composition(comp)  # dict: elem name -> (nzon,) mass frac

    # --- Detect which elements are actually present in the model ---
    elem_present = detect_elements(mf)
    elem_absent = [e for e in STELLA_ELEM_ORDER if e not in elem_present]
    print(f"\nElement detection: {len(elem_present)}/{len(STELLA_ELEM_ORDER)} elements present")
    print(f"  Present: {', '.join(elem_present)}")
    if elem_absent:
        print(f"  Absent (excluded from .mod): {', '.join(elem_absent)}")

    # --- Determine species list (atomic numbers) ---
    if species_list is not None:
        # User-specified: validate against present elements, keep order.
        present_z = {STELLA_ELEM_Z[e] for e in elem_present}
        species_z = [z for z in species_list if z in present_z]
        skipped = [z for z in species_list if z not in present_z]
        if skipped:
            print(f"  WARNING: requested species {skipped} have zero abundance; excluded")
    else:
        # Default: all present STELLA elements, in canonical order.
        species_z = [STELLA_ELEM_Z[e] for e in elem_present]
    # Reverse map: Z -> element name (for mass-fraction column lookup).
    z_to_elem = {STELLA_ELEM_Z[e]: e for e in STELLA_ELEM_ORDER}
    n_species = len(species_z)
    print(f"  Species (Z): {species_z}  ({[z_to_elem[z] for z in species_z]})")

    # --- Compute robust cumulative optical depth (outer->inner) ---
    tau, kappa_hybrid = compute_hybrid_tau(r, rho, kappa, T, mf, verbose=True)

    # --- Tau-truncation: drop inner zones with tau > tau_max ---
    # tau is monotonic increasing inward (outer->inner).  Keep zones with
    # tau <= tau_max.  Always keep at least min_zones.
    mask = tau <= tau_max
    if mask.sum() < min_zones:
        mask = np.zeros(n_orig, dtype=bool)
        mask[:min_zones] = True
        print(f"  [Truncation] tau>{tau_max} leaves <{min_zones} zones; "
              f"keeping outermost {min_zones}")
    n_trunc = int(mask.sum())
    r_t = r[mask]
    v_t = v_kms[mask]
    T_t = T[mask]
    rho_t = rho[mask]
    kappa_t = kappa_hybrid[mask]
    tau_t = tau[mask]
    flux_t = flux[mask]
    mf_t = {k: v[mask] for k, v in mf.items()}

    # Diagnostic mass estimate using zone-centre spacing (outer->inner).
    def _mass_diag(rr, rhoo):
        dr = np.abs(np.diff(rr, append=rr[-1] * 2 - rr[-2]))
        vol = (4.0 / 3.0) * np.pi * np.abs((rr + dr) ** 3 - rr ** 3)
        return np.sum(rhoo * vol)
    mass_orig = _mass_diag(r, rho)
    mass_trunc = _mass_diag(r_t, rho_t)
    print(f"\n[Truncation] {n_orig} -> {n_trunc} zones (tau <= {tau_max})")
    print(f"  Dropped {n_orig - n_trunc} optically-thick inner zones")
    print(f"  Retained mass: {mass_trunc / MSUN:.4f} Msun "
          f"({100 * mass_trunc / max(mass_orig, 1e-30):.1f}% of total)")
    print(f"  Retained tau range: {tau_t.min():.3e} - {tau_t.max():.3e}")

    # --- Optional density floor (drop near-vacuum outer zones) ---
    if rho_floor > 0:
        dmask = rho_t >= rho_floor
        if dmask.sum() >= min_zones:
            r_t, v_t, T_t, rho_t, kappa_t = r_t[dmask], v_t[dmask], T_t[dmask], rho_t[dmask], kappa_t[dmask]
            tau_t = tau_t[dmask]
            flux_t = flux_t[dmask]
            mf_t = {k: v[dmask] for k, v in mf_t.items()}
            print(f"[Density floor] rho>={rho_floor}: {n_trunc} -> {int(dmask.sum())} zones")
            n_trunc = int(dmask.sum())
        else:
            print(f"[Density floor] rho>={rho_floor} would leave <{min_zones} zones; keeping all")

    # --- Reverse to inner->outer (ascending r) for SEDONA .mod ---
    r_io = r_t[::-1].copy()
    v_io = v_t[::-1].copy()
    T_io = T_t[::-1].copy()
    rho_io = rho_t[::-1].copy()
    tau_io = tau_t[::-1].copy()
    mf_io = {k: v[::-1].copy() for k, v in mf_t.items()}
    flux_io = flux_t[::-1].copy()

    # --- Optional downsampling ---
    if nzon_target > 0 and len(r_io) > nzon_target:
        r_io, v_io, T_io, rho_io, mf_io, tau_io = downsample_grid(
            r_io, v_io, T_io, rho_io, mf_io, tau_io, nzon_target)

    n_final = len(r_io)
    v_cms = v_io * 1.0e5  # km/s -> cm/s

    # --- Build species mass-fraction matrix (n_final, n_species) ---
    species_mf = np.zeros((n_final, n_species))
    for j, z in enumerate(species_z):
        elem = z_to_elem[z]
        if elem in mf_io:
            species_mf[:, j] = mf_io[elem]

    # Remove species columns that are all zeros (no atoms of this element).
    non_zero_mask = np.any(species_mf > 0, axis=0)
    if not np.all(non_zero_mask):
        removed_z = [species_z[i] for i in range(n_species) if not non_zero_mask[i]]
        print(f"  Removing zero-abundance species from .mod: {removed_z}")
        species_z = [species_z[i] for i in range(n_species) if non_zero_mask[i]]
        species_mf = species_mf[:, non_zero_mask]
        n_species = len(species_z)
        print(f"  Remaining species (Z): {species_z}")

    # Normalise mass fractions to sum to 1 per zone.
    for i in range(n_final):
        total = np.sum(species_mf[i, :])
        if total > 0:
            species_mf[i, :] /= total

    # --- Header values ---
    t_model = time_days if time_days is not None else t_day
    t_sec = t_model * 24.0 * 3600.0
    r_max = r_io[-1]

    # rmin: inner boundary radius.  SEDONA creates zone 0 from rmin to r_out[0].
    # Must NOT be 0 for a trimmed shell (creates a huge fake zone 0 with
    # ~Msun of fake mass and tau~10^3, trapping all photons).
    if homologous:
        # Homologous: rmin is not used by SEDONA (it computes from t and v),
        # but write a sensible value just inside the first zone.
        rmin = r_io[0] * 0.99
    else:
        # Static snapshot: rmin just inside r[0] (half a zone width inward).
        # Guard against rmin=0: when zones are widely spaced (e.g. after heavy
        # downsampling) r[0]-dr0/2 can be negative.  In that case fall back to
        # 99% of the inner radius so the fake zone 0 stays small.
        if n_final > 1:
            dr0 = abs(r_io[1] - r_io[0])
        else:
            dr0 = r_io[0] * 0.01
        rmin = r_io[0] - 0.5 * dr0
        if rmin <= 0 or rmin < 0.5 * r_io[0]:
            rmin = 0.99 * r_io[0]

    # --- Data verification ---
    print(f"\n  Data verification before writing:")
    print(f"    n_zones: {n_final}")
    print(f"    r_cm: {r_io.min():.3e} - {r_io.max():.3e} cm (ascending: {bool(np.all(np.diff(r_io) > 0))})")
    print(f"    v_cms: {v_cms.min():.3e} - {v_cms.max():.3e} cm/s")
    print(f"    rho: {rho_io.min():.3e} - {rho_io.max():.3e} g/cm^3")
    print(f"    T: {T_io.min():.1f} - {T_io.max():.1f} K")
    print(f"    rmin: {rmin:.6e} cm, t_exp: {t_sec:.1f} s ({t_model:.3f} d)")

    # Monotonicity check: SEDONA requires strictly ascending r.
    if n_final >= 2 and not np.all(np.diff(r_io) > 0):
        n_degen = int(np.sum(np.diff(r_io) <= 0))
        raise ValueError(
            f"SEDONA .mod requires strictly ascending radius, but {n_degen} "
            f"non-increasing step(s) found.  Check the .swd grid for duplicate "
            f"radii (possible at shock discontinuities).")

    # --- Write .mod file ---
    with open(output_file, 'w') as f:
        # Line 1: format identifier.
        f.write("1D_sphere standard\n")
        # Line 2: nzon  rmin  texp  n_elems
        #   rmin: inner boundary radius (cm).  See notes above.
        #   texp: time since explosion in seconds (used for t_now in SEDONA).
        f.write(f"{n_final}\t{rmin:.8e}\t{t_sec:.1f}\t{n_species}\n")
        if homologous:
            print(f"  Header (homologous): {n_final} zones, rmin={rmin:.3e} cm, "
                  f"t_now={t_sec:.1f}s ({t_model:.3f}d)")
        else:
            print(f"  Header (static): {n_final} zones, rmin={rmin:.3e} cm, "
                  f"t_now={t_sec:.1f}s ({t_model:.3f}d)")
        # Line 3: species atomic numbers.
        f.write(" ".join(str(z) for z in species_z) + "\n")
        # Data lines: r, v, rho, T, mass_fracs...
        # Use .8e precision for r, v, rho, T to avoid zero-thickness zones
        # when adjacent zones differ by less than 1e-4 * r.
        for i in range(n_final):
            line = f"{r_io[i]:.8e} {v_cms[i]:.8e} {rho_io[i]:.8e} {T_io[i]:.8e}"
            for j in range(n_species):
                line += f" {species_mf[i, j]:.4f}"
            line += "\n"
            f.write(line)

    print(f"\nWritten {output_file}: {n_final} zones, {n_species} species")
    print(f"  Time: {t_model:.3f} days")
    print(f"  Radius range: [{r_io[0]:.3e}, {r_io[-1]:.3e}] cm")
    print(f"  Velocity range: [{v_cms[0]:.3e}, {v_cms[-1]:.3e}] cm/s")
    print(f"  Density range: [{rho_io[0]:.3e}, {rho_io[-1]:.3e}] g/cm^3")
    print(f"  Temperature range: [{T_io[0]:.1f}, {T_io[-1]:.1f}] K")

    # --- Inner-boundary diagnostics for the orchestrator ---
    # The orchestrator (run_model_grid.py) reads rmin, T[0] and computes
    # L = 4*pi*r^2*sigma*T^4 from the first zone for the core photon source.
    r_inner = rmin
    T_inner = float(T_io[0])
    luminosity = 4.0 * np.pi * r_io[0] ** 2 * SIGMA_SB * T_io[0] ** 4
    print(f"\n  Inner boundary (for param.lua core source):")
    print(f"    r_inner (core_radius): {r_inner:.6e} cm")
    print(f"    T[0] (core_temperature): {T_inner:.1f} K")
    print(f"    L (4*pi*r^2*sigma*T^4): {luminosity:.3e} erg/s")
    # STELLA per-zone flux at the innermost retained zone (sanity check).
    print(f"    STELLA flux at innermost retained zone: {flux_io[0]:.3e} erg/s")

    return {
        'n_zones': n_final,
        'n_species': n_species,
        'species_z': species_z,
        'rmin': rmin,
        't_sec': t_sec,
        't_days': t_model,
        'r_inner': r_inner,
        'T_inner': T_inner,
        'luminosity': luminosity,
    }


def main():
    parser = argparse.ArgumentParser(
        description='Convert STELLA .swd directly to SEDONA .mod (one step).')
    parser.add_argument('swd_file', help='Path to STELLA .swd file')
    parser.add_argument('output_file', help='Path to output .mod file')
    parser.add_argument('--target-days', type=float, default=None,
                        help='Pick epoch closest to T days (default: last epoch)')
    parser.add_argument('--time-days', type=float, default=None,
                        help='Time since explosion in days written to .mod header '
                             '(default: the selected epoch time)')
    parser.add_argument('--tau-max', type=float, default=100.0,
                        help='Truncate inner zones with cumulative tau > X (default: 100)')
    parser.add_argument('--min-zones', type=int, default=20,
                        help='Always keep at least N zones after truncation (default: 20)')
    parser.add_argument('--nzon-target', type=int, default=0,
                        help='Downsample to ~N zones (0 = no downsampling, default: 0)')
    parser.add_argument('--t-floor', type=float, default=2000.0,
                        help='Temperature floor in K (default: 2000)')
    parser.add_argument('--rho-floor', type=float, default=0.0,
                        help='Drop outer zones with rho < G g/cm^3 (default: 0 = off)')
    parser.add_argument('--species', type=str, default=None,
                        help='Comma-separated atomic numbers (e.g. 1,2,6,8,14,20,26)')
    parser.add_argument('--homologous', action='store_true',
                        help='Write header for hydro_module="homologous" (Sobolev line '
                             'expansion). Otherwise writes a static-snapshot header.')
    parser.add_argument('--use-T-rad', action='store_true',
                        help='Use radiation temperature instead of gas temperature.')
    parser.add_argument('--model-name', type=str, default=None,
                        help='Model name for diagnostics (default: .swd filename stem)')
    args = parser.parse_args()

    species_list = None
    if args.species:
        species_list = [int(x.strip()) for x in args.species.split(',') if x.strip()]

    convert_swd_to_sedona(
        args.swd_file, args.output_file,
        target_days=args.target_days,
        time_days=args.time_days,
        tau_max=args.tau_max,
        min_zones=args.min_zones,
        nzon_target=args.nzon_target,
        t_floor=args.t_floor,
        rho_floor=args.rho_floor,
        species_list=species_list,
        homologous=args.homologous,
        use_T_rad=args.use_T_rad,
        model_name=args.model_name,
    )


if __name__ == '__main__':
    main()
