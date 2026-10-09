#!/usr/bin/env python3
"""
Generate STELLA vladsf opacity tables with physically-motivated absorption.

STELLA requires pre-tabulated opacity files test.1 through test.6 (and test.ab)
in run/vladsf/. For opacityHomo.inc models, these contain simple electron
scattering opacity (kappa_es ~ 0.34 cm^2/g) independent of frequency,
temperature, density, and zone.

Table format (each file = one Fortran unformatted WRITE record):
  WRITE(nu) nw, Stime, Nfreq, Msta, Nrho, NTp,
            Wavel(Nfreq), Freq(Nfreq), Freqmn(Nfreq),
            TpTab(Ntab), RhoTab(Ntab),
            hpbansc1(Nfreq, Ntab, Ntab, Mzon/skip+1)

- nw: integer, table index (1..Ns)
- Stime: float64, time in days for this table
- Nfreq: integer, number of frequency bins (=100 for Nfreq=96 + 4 guard)
- Msta: integer = Ns (number of time slices = 6)
- Nrho, NTp: integer = Ntab = 14
- Wavel: float64(Nfreq), wavelengths in Angstrom
- Freq, Freqmn: float64(Nfreq)
- TpTab: float64(Ntab), log(T) [natural log] for tabulation
- RhoTab: float64(Ntab), log(rho) [natural log] for tabulation
- hpbansc1: float32(Nfreq, Ntab, Ntab, Mzon/skip+1), log(kappa * URHO * rho / rho_stella)
            = log(kappa_es * 10^ULGR)   (constant for homogeneous e-scattering)

The .ab file contains hpbanab1 (Rosseland-mean absorption opacity), computed
from three physical sources on the T/rho table grid:

  1. Free-free (Kramers):
       kappa_ff = 3.68e22 * (1 - 0.75*fHe) * rho * T^{-3.5}  [cm^2/g]
     Dominates at T > 15000 K or very high density.

  2. H bound-free (Saha + hydrogenic cross-section):
       kappa_bf = CHIBF * (N_H0/rho) * <alpha_bf>_Ross
     where N_H0/rho is the neutral H fraction from the Saha equation and
     CHIBF = 2.815e-6 is the CMFGEN bound-free constant.  Dominates at
     T = 3000-15000 K (photospheric regime) where H is partially neutral.
     This is the dominant term that was missing from the old Kramers-only
     treatment and caused Mbol >> optical bands.

  3. H- free-free:
       kappa_Hminus = 1.38e-25 * P_e * T^{-0.5} * rho^{-1}  [cm^2/g]
     where P_e = Ne * k * T is the electron pressure.  Dominates at
     T < 8000 K in the dense cool RSG shell.

All three are summed and the result is clamped to [kappa_es * floor, kappa_es].

If --opal-tables is given (recommended for composition-gradient models), the
absorption opacity is instead read from OPAL Rosseland mean tables (Iglesias &
Rogers 1996) via the compiled opal_opacity.so extension (built from SNEC's
opal_opacity.F90).  Per-zone composition is read from --input-profile, giving
a zone-dependent hpbanab1 that correctly handles both H-rich inner zones and
H-poor metal-rich outer zones in the same model.

Usage:
    python generate_stella_opacity_tables.py <vladsf_dir> <model_name> \\
        [--nzon NZON] [--kappa KAPPA] [--nfreq NFREQ] [--ntab NTAB] [--ns NS]
        [--fh FH] [--fhe FHE] [--abs-fraction-floor F]
        [--opal-tables] [--opal-table-dir DIR] [--input-profile FILE]
"""

import numpy as np
import struct
import os
import sys
import argparse
import importlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


def write_fortran_record(f, *arrays):
    """Write a Fortran unformatted sequential record (little-endian 4-byte markers)."""
    data = b''
    for arr in arrays:
        if isinstance(arr, np.ndarray):
            data += arr.tobytes()
        elif isinstance(arr, (np.int32, np.int64)):
            data += struct.pack('<i', int(arr))
        elif isinstance(arr, (np.float32, np.float64)):
            data += arr.tobytes()
        elif isinstance(arr, int):
            data += struct.pack('<i', arr)
        elif isinstance(arr, float):
            data += struct.pack('<d', arr)
        else:
            raise TypeError(f"Unsupported type: {type(arr)}")
    length = len(data)
    marker = struct.pack('<i', length)
    f.write(marker)
    f.write(data)
    f.write(marker)


def detect_stella_skip(stella_root, nzon):
    """
    Read @skip from STELLA source (opacityHomo.inc or opacity45.inc) and
    Mzon from zone.inc, then evaluate the skip expression to return the
    correct nim = Mzon/@skip + 1 matching the compiled binary.

    The .inc files use a macro language where consecutive non-commented lines
    override previous ones, so the last active (non-comment) @skip line wins.

    Falls back to nzon//2 if the STELLA source files are not found.
    """
    import re

    # Locate the source directory
    src_dir = None
    for candidate in [
        os.path.join(stella_root, 'src'),
        os.path.join(stella_root, '..', 'src'),
    ]:
        if os.path.isdir(candidate):
            src_dir = os.path.abspath(candidate)
            break
    if src_dir is None:
        return max(1, nzon // 2)

    # Try opacityHomo.inc first (used by homogeneous runs), then opacity45.inc
    skip_val = None
    for inc_name in ('opacityHomo.inc', 'opacity45.inc'):
        inc_path = os.path.join(src_dir, inc_name)
        if not os.path.exists(inc_path):
            continue
        with open(inc_path) as f:
            lines = f.readlines()
        # Find all non-comment @skip lines (comments start with --)
        for line in lines:
            stripped = line.strip()
            if stripped.startswith('--'):
                continue
            # Match "@skip <expr>" where expr may be Mzon/N or just a number
            m = re.search(r'@skip\s+([\w/()]+)', stripped)
            if m:
                expr = m.group(1).strip()
                # Evaluate: replace Mzon with nzon
                try:
                    skip_val = int(eval(expr.replace('Mzon', str(nzon))))
                except Exception:
                    continue
        if skip_val is not None:
            break

    if skip_val is None or skip_val < 1:
        skip_val = max(1, nzon // 2)

    nim = nzon // skip_val + 1
    print(f"  STELLA @skip={skip_val} (Mzon={nzon}) -> nim={nim} zone groups  "
          f"[from {inc_name if skip_val else 'fallback'}]")
    return skip_val


def saha_neutral_H_fraction(T_grid, rho_grid, fH=0.70, fHe=0.28):
    """
    Compute the neutral hydrogen number density fraction x_H0 = N_H0 / N_H
    on the (Ntab_rho, Ntab_T) grid using the Saha equation.

    Saha:  N_H+ * Ne / N_H0 = (2 * u_H+ / u_H0) * (2*pi*m_e*k*T / h^2)^1.5
                              * exp(-chi_H / kT)
    where chi_H = 13.6 eV, u_H+/u_H0 = 1 (ground state partition functions).

    Parameters
    ----------
    T_grid   : (Ntab_T,)   temperature grid [K]
    rho_grid : (Ntab_rho,) density grid [g/cm^3]
    fH       : hydrogen mass fraction (default 0.70)
    fHe      : helium mass fraction (default 0.28)

    Returns
    -------
    x_H0 : (Ntab_rho, Ntab_T) array, fraction of H that is neutral [0..1]
    """
    T   = np.asarray(T_grid,   dtype=np.float64)    # (Ntab_T,)
    rho = np.asarray(rho_grid, dtype=np.float64)    # (Ntab_rho,)

    # Physical constants (CGS)
    k_B  = 1.3807e-16   # erg/K
    h    = 6.6261e-27   # erg s
    m_e  = 9.1094e-28   # g
    m_H  = 1.6726e-24   # g
    m_He = 6.6465e-24   # g
    eV   = 1.6022e-12   # erg
    chi_H = 13.6 * eV   # H ionisation potential

    # Number density of H nuclei per cm^3:  N_H = fH * rho / m_H
    N_H = fH * rho[:, None] / m_H          # (Ntab_rho, Ntab_T)

    # Saha RHS: Phi(T) = 2*(2*pi*m_e*k*T/h^2)^1.5 * exp(-chi_H/kT)
    Phi = 2.0 * (2.0 * np.pi * m_e * k_B * T[None, :] / h**2)**1.5 \
          * np.exp(-chi_H / (k_B * T[None, :]))   # (1, Ntab_T) -> broadcast

    # Saha: x^2/(1-x) * N_H = Phi  where x = ionisation fraction
    # Solve quadratic: x^2 + (Phi/N_H)*x - Phi/N_H = 0
    # x = [-b + sqrt(b^2 + 4b)] / 2  where b = Phi/N_H
    b = Phi / np.maximum(N_H, 1e-300)       # (Ntab_rho, Ntab_T)
    x_ion = 0.5 * (-b + np.sqrt(b**2 + 4.0 * b))   # ionisation fraction
    x_ion = np.clip(x_ion, 0.0, 1.0)
    x_H0  = 1.0 - x_ion                    # neutral fraction
    return x_H0


def physical_absorption(T_grid, rho_grid, kappa_es=0.34, fH=0.70, fHe=0.28,
                        floor=1e-4):
    """
    Compute total Rosseland-mean absorption opacity / kappa_es on the
    (Ntab_rho, Ntab_T) table grid, summing three physical sources:

      1. Free-free (Kramers):  kappa_ff = 3.68e22*(1-0.75*fHe)*rho*T^-3.5
      2. H bound-free (Saha):  kappa_bf = CHIBF * N_H0 * <alpha>_Ross / rho
      3. H- free-free:         kappa_Hm = 1.38e-25 * P_e * T^-0.5 / rho

    Constants from CMFGEN source:
      CHIBF = 2.815e-6  (bound-free Gaunt-factor-weighted cross-section constant)
      CHIFF = 3.69e-29  (free-free constant, -> Kramers via Rosseland integral)

    Returns
    -------
    frac_2d : (Ntab_rho, Ntab_T) array of (kappa_abs / kappa_es), clipped to [floor, 1]
    kappa_ff_2d, kappa_bf_2d, kappa_hm_2d : individual components [cm^2/g]
    """
    T   = np.asarray(T_grid,   dtype=np.float64)    # (Ntab_T,)
    rho = np.asarray(rho_grid, dtype=np.float64)    # (Ntab_rho,)

    # Physical constants
    k_B  = 1.3807e-16   # erg/K
    m_H  = 1.6726e-24   # g
    m_e  = 9.1094e-28   # g
    eV   = 1.6022e-12   # erg

    # --- 1. Free-free (Kramers Rosseland mean) ---
    # kappa_ff [cm^2/g]; factor (1-0.75*fHe) accounts for He dilution
    kappa_ff = (3.68e22 * (1.0 - 0.75 * fHe)
                * rho[:, None] * T[None, :] ** (-3.5))   # (Ntab_rho, Ntab_T)

    # --- 2. H bound-free (Saha + Kramers Rosseland-mean cross-section) ---
    # Rosseland-mean BF opacity (Kramers approximation, Cox & Giuli 1968 eq.17.48):
    #   kappa_bf = 4.34e25 * (x_H0 * fH) * rho * T^{-3.5}  [cm^2/g]
    # x_H0 * fH is the neutral H mass fraction from the Saha equation.
    # Same T^{-3.5} Rosseland weighting as free-free; rho scaling is identical.
    # At the photosphere (T~8000K, rho~1e-12, x_H0~0.04): kappa_bf ~ 40*kappa_ff.
    x_H0 = saha_neutral_H_fraction(T_grid, rho_grid, fH=fH, fHe=fHe)
    # x_H0 * fH is the neutral H mass fraction; rho broadcasts over T axis.
    # Coefficient 4.34e25 subsumes m_H, the Gaunt factor, and the Rosseland
    # frequency integral (Cox & Giuli 1968, eq. 17.48).
    kappa_bf = (4.34e25 * x_H0 * fH
                * rho[:, None] * T[None, :] ** (-3.5))  # (Ntab_rho, Ntab_T)

    # --- 3. H- free-free ---
    # Dominant at T < 8000 K.  Uses electron pressure P_e = Ne * k * T.
    # Ne estimated from H ionisation (Saha) + He first ionisation.
    # kappa_Hm = 1.38e-25 * P_e / (rho * T^0.5)  [cm^2/g]  (Bell & Berger 1987)
    x_ion_H = 1.0 - x_H0  # H ionisation fraction
    # Ne from H: N_e_H = x_ion * N_H
    Ne_H = x_ion_H * fH * rho[:, None] / m_H
    # He is singly ionised for T > 10000 K, neutral below: simple step approx
    Ne_He = np.where(T[None, :] > 1e4,
                     fHe * rho[:, None] / (4.0 * m_H), 0.0)
    Ne = Ne_H + Ne_He   # (Ntab_rho, Ntab_T)
    P_e = Ne * k_B * T[None, :]   # erg/cm^3
    kappa_Hm = (1.38e-25 * P_e
                / (np.maximum(rho[:, None], 1e-300) * T[None, :] ** 0.5))

    # --- Sum and normalise ---
    kappa_abs = kappa_ff + kappa_bf + kappa_Hm
    frac_2d = np.clip(kappa_abs / kappa_es, floor, 1.0)

    return frac_2d, kappa_ff, kappa_bf, kappa_Hm


def load_opal_module(scripts_dir=None):
    """
    Load the compiled opal_opacity extension module.
    Searches: scripts_dir, directory of this file, then sys.path.
    Returns the module or raises ImportError with a helpful message.
    """
    search_dirs = []
    if scripts_dir:
        search_dirs.append(scripts_dir)
    search_dirs.append(os.path.dirname(os.path.abspath(__file__)))
    for d in search_dirs:
        if d not in sys.path:
            sys.path.insert(0, d)
    try:
        return importlib.import_module('opal_opacity')
    except ImportError as e:
        raise ImportError(
            f"Cannot import opal_opacity module: {e}\n"
            f"Searched: {search_dirs}\n"
            f"Build it with: scripts/build_opal_module.sh"
        ) from e


def opal_absorption(T_grid, rho_grid, profile_path, opal_table_dir,
                    kappa_es=0.34, floor=1e-4, ntab=14, nim=1,
                    scripts_dir=None):
    """
    Compute absorption opacity using OPAL Rosseland mean tables (Iglesias &
    Rogers 1996), with per-zone composition read from the STELLA input profile.

    The OPAL tables give the total Rosseland mean opacity (scattering +
    absorption).  Absorption is extracted as kappa_abs = kappa_Ross - kappa_es,
    clamped to [kappa_es * floor, kappa_es].

    Returns hpbanab1 of shape (ntab_rho, ntab_T, nim) ready to broadcast into
    the 4D STELLA table.

    Parameters
    ----------
    T_grid        : (ntab,) temperature grid [K]
    rho_grid      : (ntab,) density grid [g/cm^3]
    profile_path  : path to STELLA input_profile
                    columns: r_cm rho p u fH fHe fC fO fSi fFe
    opal_table_dir: directory containing the tables/ subdirectory (codataa..e)
    kappa_es      : electron scattering opacity [cm^2/g]
    floor         : minimum kappa_abs/kappa_es
    ntab          : number of T/rho grid points (must match T_grid/rho_grid)
    nim           : number of zone groups (Mzon/skip + 1)
    scripts_dir   : directory containing opal_opacity.so
    """
    opal = load_opal_module(scripts_dir)

    # --- Load composition profile ---
    # Profile format (12 columns): r_cm rho_cgs p_cgs u_cgs T_K fx_H fx_He fx_O fx_Si fx_Ca fx_Fe fx_Ni
    # Old format (10 columns, no T): r_cm rho p u fH fHe fC fO fSi fFe
    # Auto-detect by checking whether column 4 looks like T (~1e3-1e6) or fH (~0-1)
    prof = np.loadtxt(profile_path, comments='#')
    n_prof   = prof.shape[0]
    r_prof   = prof[:, 0]
    rho_prof = prof[:, 1]

    # Detect column layout: if col 4 values are >> 1, it's T (12-col format)
    if np.max(np.abs(prof[:, 4])) > 10.0:
        # 12-column format: T at index 4, composition starts at index 5
        fH_prof  = prof[:, 5]
        fHe_prof = prof[:, 6] if prof.shape[1] > 6 else np.full(n_prof, 0.28)
        fC_prof  = np.zeros(n_prof)  # no C column in 12-col format
        fO_prof  = prof[:, 7] if prof.shape[1] > 7 else np.zeros(n_prof)
    else:
        # 10-column format (legacy): composition starts at index 4
        fH_prof  = prof[:, 4]
        fHe_prof = prof[:, 5] if prof.shape[1] > 5 else np.full(n_prof, 0.28)
        fC_prof  = prof[:, 6] if prof.shape[1] > 6 else np.zeros(n_prof)
        fO_prof  = prof[:, 7] if prof.shape[1] > 7 else np.zeros(n_prof)

    # --- Composition-aware zone group sampling ---
    # STELLA maps spatial zone iz -> opacity group iz//skip, so groups are
    # uniform in zone-index space, which is approximately uniform in log-r.
    # We therefore divide the profile into nim equal log-r intervals, then
    # pin additional boundaries at every major composition interface (|dfH|>0.05)
    # so that H-rich and H-poor layers always land in separate groups.
    dr      = np.abs(np.diff(r_prof, prepend=r_prof[0]))
    dm      = rho_prof * dr                         # mass proxy per cell

    log_r   = np.log10(np.maximum(r_prof, 1e-30))
    lr_min, lr_max = log_r[0], log_r[-1]

    # Base equal-log-r boundaries
    base_edges_lr = np.linspace(lr_min, lr_max, nim + 1)

    # Detect composition interfaces: cells where |dfH| > 0.05
    dfH = np.abs(np.diff(fH_prof, prepend=fH_prof[0]))
    iface_idx = np.where(dfH > 0.05)[0]
    iface_lr  = log_r[iface_idx]

    # Insert interface positions into boundary list and re-select nim+1 boundaries
    all_edges = np.unique(np.concatenate([base_edges_lr, iface_lr]))
    # Greedily select nim+1 boundaries keeping interface positions, spaced > 0
    final_edges = [lr_min]
    for e in all_edges[1:]:
        if e > final_edges[-1] + 1e-6:
            final_edges.append(e)
    final_edges = np.array(final_edges)
    # Reduce to nim+1 by thinning non-interface edges if we have too many
    while len(final_edges) > nim + 1:
        # Find smallest gap that is NOT adjacent to an interface
        gaps = np.diff(final_edges)
        iface_set = set(np.searchsorted(final_edges, iface_lr, side='left'))
        remove_candidates = [i for i in range(len(gaps))
                             if i not in iface_set and (i+1) not in iface_set]
        if not remove_candidates:
            remove_candidates = list(range(len(gaps)))
        i_remove = remove_candidates[int(np.argmin(gaps[remove_candidates]))]
        # Merge by removing final_edges[i_remove+1]
        final_edges = np.delete(final_edges, i_remove + 1)
    # Pad if too few
    while len(final_edges) < nim + 1:
        gaps = np.diff(final_edges)
        i_split = int(np.argmax(gaps))
        new_edge = 0.5 * (final_edges[i_split] + final_edges[i_split + 1])
        final_edges = np.insert(final_edges, i_split + 1, new_edge)

    # Convert log-r edges to profile-cell boundaries
    cell_edges = np.searchsorted(log_r, final_edges)
    cell_edges[0]  = 0
    cell_edges[-1] = n_prof

    zone_fH  = np.zeros(nim)
    zone_fHe = np.zeros(nim)
    zone_fC  = np.zeros(nim)
    zone_fO  = np.zeros(nim)
    for iz in range(nim):
        i0 = cell_edges[iz]
        i1 = cell_edges[iz + 1]
        if i1 <= i0:
            i1 = i0 + 1
        i1 = min(i1, n_prof)
        w    = dm[i0:i1]
        wsum = w.sum() if w.sum() > 0 else 1.0
        zone_fH[iz]  = np.sum(fH_prof[i0:i1]  * w) / wsum
        zone_fHe[iz] = np.sum(fHe_prof[i0:i1] * w) / wsum
        zone_fC[iz]  = np.sum(fC_prof[i0:i1]  * w) / wsum
        zone_fO[iz]  = np.sum(fO_prof[i0:i1]  * w) / wsum

    print(f"  OPAL zone compositions (log-r + interface-aware, {nim} groups):")
    for iz in range(nim):
        i0, i1 = int(cell_edges[iz]), int(min(cell_edges[iz+1], n_prof))
        print(f"    zone {iz}: cells [{i0},{i1})  "
              f"r=[{r_prof[i0]:.2e},{r_prof[min(i1,n_prof-1)]:.2e}]  "
              f"fH={zone_fH[iz]:.3f} fHe={zone_fHe[iz]:.3f} "
              f"fC={zone_fC[iz]:.4f} fO={zone_fO[iz]:.4f}")

    # OPAL parameters: Z=0.02 (solar metallicity tables)
    # xxC, xxO are enhancements ABOVE the solar amount already in Z
    Z_opal = np.float32(0.02)
    C_sol  = 0.00245   # solar C mass fraction in Z=0.02
    O_sol  = 0.00653   # solar O mass fraction in Z=0.02

    # Build flat (T, rho) grid for vectorised OPAL call
    # OPAL variables: t6 = T * 1e-6,  R = rho / t6^3
    # OPAL valid range: logT in [3.75, 8.7]  => t6 in [5.62e-3, 501]
    #                   logR in [-8, +1]      => R  in [1e-8, 10]
    T6_MIN = 10**3.75 * 1e-6   # 0.005623
    T6_MAX = 10**8.70 * 1e-6   # 501.19
    R_MIN  = 10**-8.0           # 1e-8
    R_MAX  = 10**+1.0           # 10.0

    T   = np.asarray(T_grid,   dtype=np.float64)
    rho = np.asarray(rho_grid, dtype=np.float64)
    T_2d, rho_2d = np.meshgrid(T, rho, indexing='ij')   # (ntab_T, ntab_rho)
    t6_flat  = np.clip((T_2d.ravel() * 1e-6), T6_MIN, T6_MAX).astype(np.float32)
    rho_flat = rho_2d.ravel().astype(np.float32)
    R_flat   = np.clip(rho_flat / np.maximum(t6_flat, 1e-30)**3,
                       R_MIN, R_MAX).astype(np.float32)
    n_pts    = len(t6_flat)

    urho = 1.0e-6   # STELLA URHO
    ur   = 1.0e14   # STELLA UR

    # Output: (ntab_rho, ntab_T, nim)
    logkabs_3d = np.zeros((ntab, ntab, nim), dtype=np.float32)

    print(f"  OPAL: {n_pts} (T,rho) points x {nim} zone groups")
    for iz in range(nim):
        xh  = np.float32(zone_fH[iz])
        xxc = np.float32(max(0.0, float(zone_fC[iz]) - float(Z_opal) * C_sol))
        xxo = np.float32(max(0.0, float(zone_fO[iz]) - float(Z_opal) * O_sol))

        log_kappa = opal.opal_grid(
            Z_opal, xh, xxc, xxo,
            t6_flat, R_flat,
            opal_table_dir
        )  # log10(kappa_Ross [cm^2/g]), shape (n_pts,)

        kappa_tot = 10.0 ** log_kappa.reshape(ntab, ntab)  # (ntab_T, ntab_rho)

        # Absorption = total Rosseland - electron scattering, clamped
        kappa_abs = np.clip(kappa_tot - kappa_es,
                            kappa_es * floor, kappa_es)

        # Transpose to (ntab_rho, ntab_T) for STELLA index convention
        logk = np.log(np.maximum(kappa_abs.T * urho * ur, 1e-30)).astype(np.float32)
        logkabs_3d[:, :, iz] = logk

        if iz == 0 or iz == nim - 1 or iz % max(1, nim // 4) == 0:
            iT_ref, ir_ref = 3, 6   # T~1e4K, rho~1e-11
            kabs_ref = kappa_abs[iT_ref, ir_ref]
            print(f"    zone {iz:3d}: xH={xh:.3f} xxC={xxc:.3f} xxO={xxo:.3f}  "
                  f"kappa_abs(T=1e4,rho=1e-11)={kabs_ref:.3e} cm2/g")

    return logkabs_3d


def make_opacity_plots(T_grid, rho_grid, hpbanab1, hpbansc1,
                       model_name, plot_dir, kappa_es=0.34,
                       urho=1e-6, ur=1e14, nim=1):
    """
    Generate diagnostic plots of the STELLA opacity tables.

    Plots produced:
    - opacity_heatmap.png  : 2D (T, rho) kappa_abs and kappa_es maps per zone
    - opacity_vs_T.png     : kappa_abs vs T at fixed rho for all zones
    - opacity_vs_rho.png   : kappa_abs vs rho at fixed T for all zones
    - opacity_ratio.png    : kappa_abs/kappa_es ratio map zone 0 vs last zone
    """
    os.makedirs(plot_dir, exist_ok=True)
    logT  = np.log10(T_grid)
    logrho = np.log10(rho_grid)
    colors = plt.cm.plasma(np.linspace(0.1, 0.9, max(nim, 2)))

    # Decode kappa_abs from hpbanab1: shape (nfreq, ntab_T, ntab_rho, nim)
    # hpbanab1 = log(kappa_abs * urho * ur)  ->  kappa_abs = exp(val) / (urho*ur)
    # Use freq=0 slice
    logkabs_Tr = hpbanab1[0, :, :, :]   # (ntab_T, ntab_rho, nim)
    kappa_abs  = np.exp(logkabs_Tr) / (urho * ur)  # cm^2/g

    # kappa_es is constant (from hpbansc1)
    kappa_sc = np.exp(hpbansc1[0, 0, 0, 0]) / (urho * ur)

    # ---- 1. 2D heatmaps: kappa_abs per zone --------------------------------
    ncols = min(nim, 4)
    nrows = (nim + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols,
                              figsize=(4.5 * ncols, 4.0 * nrows),
                              squeeze=False)
    for iz in range(nim):
        ax = axes[iz // ncols][iz % ncols]
        data = kappa_abs[:, :, iz]  # (ntab_T, ntab_rho)
        data_safe = np.where(data > 0, data, np.nan)
        pcm = ax.pcolormesh(logrho, logT, data_safe,
                            norm=LogNorm(vmin=1e-5, vmax=1.0),
                            cmap='inferno', shading='auto')
        ax.set_xlabel(r'$\log_{10}(\rho\ [\mathrm{g/cm}^3])$', fontsize=8)
        ax.set_ylabel(r'$\log_{10}(T\ [\mathrm{K}])$', fontsize=8)
        ax.set_title(f'Zone {iz}', fontsize=9)
        plt.colorbar(pcm, ax=ax, label=r'$\kappa_{\rm abs}\ [\mathrm{cm}^2/\mathrm{g}]$')
    # Hide unused axes
    for iz in range(nim, nrows * ncols):
        axes[iz // ncols][iz % ncols].set_visible(False)
    fig.suptitle(f'{model_name}: Absorption opacity $\\kappa_{{\\rm abs}}(T, \\rho)$ per zone',
                 fontsize=11)
    plt.tight_layout()
    out = os.path.join(plot_dir, 'opacity_heatmap.png')
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Diagnostic plot: {out}")

    # ---- 2. kappa_abs vs T at photospheric rho (rho ~ 1e-11 -> ir=7) -------
    ir_phot = min(6, len(rho_grid) - 1)   # rho ~ 1e-11 g/cm^3
    ir_cool = min(10, len(rho_grid) - 1)  # rho ~ 1e-7 g/cm^3
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    for iz in range(nim):
        lbl = f'Zone {iz}'
        axes[0].plot(logT, kappa_abs[:, ir_phot, iz], color=colors[iz], lw=2, label=lbl)
        axes[1].plot(logT, kappa_abs[:, ir_cool, iz], color=colors[iz], lw=2, label=lbl)
    for ax, ir, rlab in zip(axes, [ir_phot, ir_cool],
                             [f'$\\rho=10^{{{logrho[ir_phot]:.0f}}}$ g/cm³ (photosphere)',
                              f'$\\rho=10^{{{logrho[ir_cool]:.0f}}}$ g/cm³ (cool shell)']):
        ax.axhline(kappa_sc, color='gray', ls='--', lw=1.5, label=f'$\\kappa_{{es}}$={kappa_sc:.3f}')
        ax.set_xlabel(r'$\log_{10}(T\ [\mathrm{K}])$')
        ax.set_ylabel(r'$\kappa_{\rm abs}\ [\mathrm{cm}^2/\mathrm{g}]$')
        ax.set_yscale('log')
        ax.set_title(rlab)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
    fig.suptitle(f'{model_name}: $\\kappa_{{\\rm abs}}$ vs T', fontsize=11)
    plt.tight_layout()
    out = os.path.join(plot_dir, 'opacity_vs_T.png')
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Diagnostic plot: {out}")

    # ---- 3. kappa_abs vs rho at fixed T ------------------------------------
    iT_phot = min(3, len(T_grid) - 1)   # T ~ 1e4 K
    iT_hot  = min(6, len(T_grid) - 1)   # T ~ 1e5 K
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    for iz in range(nim):
        lbl = f'Zone {iz}'
        axes[0].plot(logrho, kappa_abs[iT_phot, :, iz], color=colors[iz], lw=2, label=lbl)
        axes[1].plot(logrho, kappa_abs[iT_hot,  :, iz], color=colors[iz], lw=2, label=lbl)
    for ax, iT, Tlab in zip(axes, [iT_phot, iT_hot],
                              [f'$T=10^{{{logT[iT_phot]:.1f}}}$ K (photosphere)',
                               f'$T=10^{{{logT[iT_hot]:.1f}}}$ K (hot)']):
        ax.axhline(kappa_sc, color='gray', ls='--', lw=1.5, label=f'$\\kappa_{{es}}$={kappa_sc:.3f}')
        ax.set_xlabel(r'$\log_{10}(\rho\ [\mathrm{g/cm}^3])$')
        ax.set_ylabel(r'$\kappa_{\rm abs}\ [\mathrm{cm}^2/\mathrm{g}]$')
        ax.set_yscale('log')
        ax.set_title(Tlab)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
    fig.suptitle(f'{model_name}: $\\kappa_{{\\rm abs}}$ vs $\\rho$', fontsize=11)
    plt.tight_layout()
    out = os.path.join(plot_dir, 'opacity_vs_rho.png')
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Diagnostic plot: {out}")

    # ---- 4. kappa_abs/kappa_es ratio: zone 0 vs last zone -------------------
    if nim > 1:
        fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), sharey=True)
        for ax, iz, title in zip(axes, [0, nim-1],
                                  [f'Zone 0 (innermost)', f'Zone {nim-1} (outermost)']):
            ratio = kappa_abs[:, :, iz] / kappa_sc
            ratio_safe = np.where(ratio > 0, ratio, np.nan)
            pcm = ax.pcolormesh(logrho, logT, ratio_safe,
                                norm=LogNorm(vmin=1e-5, vmax=2.0),
                                cmap='RdYlGn_r', shading='auto')
            ax.set_xlabel(r'$\log_{10}(\rho)$')
            ax.set_ylabel(r'$\log_{10}(T)$')
            ax.set_title(title)
            plt.colorbar(pcm, ax=ax, label=r'$\kappa_{\rm abs}/\kappa_{\rm es}$')
        fig.suptitle(f'{model_name}: absorption/scattering ratio', fontsize=11)
        plt.tight_layout()
        out = os.path.join(plot_dir, 'opacity_ratio.png')
        fig.savefig(out, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f"  Diagnostic plot: {out}")


def generate_opus_tables(vladsf_dir, model_name, nzon=200, kappa_es=0.34,
                         nfreq=100, ntab=14, ns=6, skip=None,
                         wlmin=1.0, wlmax=5.0e4, ulgr=14.0,
                         abs_fraction=None, abs_fraction_floor=1e-4,
                         fH=0.70, fHe=0.28,
                         opal_tables=False, opal_table_dir=None,
                         input_profile=None, plot_dir=None,
                         stella_root=None):
    """
    Generate STELLA opacity table files.

    Parameters
    ----------
    vladsf_dir        : directory to write tables into
    model_name        : base name for files (model_name.1 ... model_name.6, model_name.ab)
    nzon              : number of STELLA zones (Mzon)
    kappa_es          : electron scattering opacity in cm^2/g (default 0.34)
    nfreq             : number of frequency bins (NFREQ in zone.inc = @Nfreq+4 = 100)
    ntab              : number of T/rho table points (@Ntab = 14)
    ns                : number of time slices (@Ns = 6)
    skip              : zone skip for opacity table (@skip). If None, auto-detected
                        from opacityHomo.inc in stella_root/src; falls back to Mzon/2.
    stella_root       : path to OpenStella root (used to auto-detect @skip).
    wlmin,wlmax       : wavelength range in Angstrom
    ulgr              : log10(UR) where UR = STELLA length unit in cm (default 14 -> UR=1e14 cm)
    abs_fraction      : if set, override full physics with a constant kappa_abs/kappa_es scalar.
                        Default None = use full physical model (ff + H bf + H-).
    abs_fraction_floor: minimum kappa_abs/kappa_es (prevents log(0), default 1e-4).
    fH                : hydrogen mass fraction for Kramers/Saha model (default 0.70).
    fHe               : helium mass fraction for Kramers/Saha model (default 0.28).
    opal_tables       : if True, use OPAL Rosseland mean tables with per-zone composition.
    opal_table_dir    : path to directory containing tables/codataa..codatae.
    input_profile     : path to STELLA input_profile with per-cell composition.
    plot_dir          : if set, save diagnostic opacity plots to this directory.
    """
    if skip is None:
        if stella_root and os.path.isdir(stella_root):
            skip = detect_stella_skip(stella_root, nzon)
        else:
            skip = max(1, nzon // 2)
            print(f"  @skip: stella_root not provided, defaulting to Mzon/2={skip}")

    nim = nzon // skip + 1  # number of mass zone groups: Mzon/@skip+1

    # --- Frequency grid (matching ronfitab.trf exactly) ---
    # Speed of light in cm/s (physical, not STELLA units)
    clight = 2.998e10
    # FREQ has Nfreq+1 bin edges, geometric progression
    # freq(1) = c/wlmax, freq(nfreq+1) = c/wlmin (in cm/s, Angstrom -> cm via *1e-8)
    freq_arr = np.zeros(nfreq + 1, dtype=np.float64)
    freq_arr[0] = clight * 1e8 / wlmax   # c (cm/s) * 1e8 (AA/cm) / wlmax(AA)
    freq_arr[nfreq] = clight * 1e8 / wlmin
    basis = (wlmax / wlmin) ** (1.0 / nfreq)
    for i in range(1, nfreq):
        freq_arr[i] = freq_arr[i-1] * basis
    # FREQMN(L) = sqrt(FREQ(L)*FREQ(L+1))  -- geometric midpoints
    freqmn_arr = np.sqrt(freq_arr[:nfreq] * freq_arr[1:nfreq+1])
    # WAVEL(L) = c/FREQMN(L) in Angstrom
    wavel = clight * 1e8 / freqmn_arr

    # --- Temperature grid ---
    # TpTab(itp) = log(T) for itp=1..Ntab, T = 10^(3.4 + 0.2*(itp-1)) K
    # Natural log (from ronfitab.trf)
    TpTab = np.array([np.log(10**(3.4 + 0.2*(itp-1))) for itp in range(1, ntab+1)],
                     dtype=np.float64)

    # --- Density grid ---
    # RhoTab(ir) = log(rho) for ir=1..Ntab, rho = 10^(-18+ir) g/cm^3
    RhoTab = np.array([np.log(10**(-18.0 + ir)) for ir in range(1, ntab+1)],
                      dtype=np.float64)

    # --- Opacity value ---
    # hpbansc1 = ln( (scatop + opac) * UR / (rho/URHO) )  where UR=10^ULGR
    # For pure electron scattering: scatop = kappa_es * rho [cm^{-1}], opac = 0
    # URHO = 10^{-6} g/cm^3, UR = 10^14 cm
    # hpbansc1 = ln(kappa_es * rho * UR / (rho / URHO))
    #          = ln(kappa_es * UR * URHO)
    #          = ln(kappa_es * 1e8)  [~ 17.34 for kappa_es=0.34]
    urho = 1.0e-6  # g/cm^3
    ur = 10**ulgr   # cm
    hapval = np.float32(np.log(kappa_es * urho * ur))

    print(f"  Opacity value (ln(kappa_es * URHO * UR)): {hapval:.4f}")
    print(f"  kappa_es = {kappa_es} cm^2/g, UR = {ur:.2e} cm, URHO = {urho:.2e} g/cm^3")

    # hpbansc1 shape: (Nfreq, Ntab, Ntab, nim) - 4D, one per time slice file
    hpbansc1 = np.full((nfreq, ntab, ntab, nim), hapval, dtype=np.float32)

    # hpbanab1 (absorption opacity)
    T_grid   = np.array([10**(3.4 + 0.2*(itp-1)) for itp in range(1, ntab+1)])
    rho_grid = np.array([10**(-18.0 + ir)         for ir  in range(1, ntab+1)])

    scripts_dir = os.path.dirname(os.path.abspath(__file__))

    if opal_tables:
        # --- OPAL path: per-zone composition from input profile ---
        if not input_profile or not os.path.exists(input_profile):
            raise FileNotFoundError(
                f"--opal-tables requires --input-profile pointing to an existing file. "
                f"Got: {input_profile}")
        if not opal_table_dir or not os.path.isdir(opal_table_dir):
            raise FileNotFoundError(
                f"--opal-table-dir must be the directory containing tables/codataa..e. "
                f"Got: {opal_table_dir}")
        print(f"  Absorption opacity: OPAL Rosseland mean (Iglesias & Rogers 1996)")
        print(f"    profile: {input_profile}")
        print(f"    tables:  {opal_table_dir}/tables/codataa..e")
        print(f"    floor:   {abs_fraction_floor:.0e}")
        # Returns (ntab_rho, ntab_T, nim)
        logkabs_3d = opal_absorption(
            T_grid, rho_grid,
            profile_path=input_profile,
            opal_table_dir=opal_table_dir,
            kappa_es=kappa_es, floor=abs_fraction_floor,
            ntab=ntab, nim=nim,
            scripts_dir=scripts_dir)
        # Broadcast to (Nfreq, Ntab_T, Ntab_rho, nim)
        # logkabs_3d is (ntab_rho, ntab_T, nim) -> transpose to (ntab_T, ntab_rho, nim)
        logkabs_3d_T = logkabs_3d.transpose(1, 0, 2)  # (ntab_T, ntab_rho, nim)
        hpbanab1 = np.broadcast_to(
            logkabs_3d_T[np.newaxis, :, :, :],
            (nfreq, ntab, ntab, nim)).copy().astype(np.float32)

    elif abs_fraction is not None:
        # User override: constant fraction across all T/rho and zones
        frac_2d = np.full((ntab, ntab), float(abs_fraction))
        print(f"  Absorption opacity: constant abs_fraction={abs_fraction:.2e} (user override)")
        kappa_abs_2d = kappa_es * frac_2d
        logkabs_2d = np.log(np.maximum(kappa_abs_2d * urho * ur, 1e-30)).astype(np.float32)
        logkabs_4d = logkabs_2d.T[np.newaxis, :, :, np.newaxis]
        hpbanab1 = np.broadcast_to(logkabs_4d,
                                    (nfreq, ntab, ntab, nim)).copy().astype(np.float32)

    else:
        # Full physical model: ff + H bound-free (Saha) + H-
        frac_2d, kappa_ff, kappa_bf, kappa_Hm = physical_absorption(
            T_grid, rho_grid, kappa_es=kappa_es,
            fH=fH, fHe=fHe, floor=abs_fraction_floor)
        print(f"  Absorption opacity: physical model (ff + H bf + H-)")
        print(f"    fH={fH}, fHe={fHe}, floor={abs_fraction_floor:.0e}")
        print(f"    kappa_ff  range: {kappa_ff.min():.2e} -- {kappa_ff.max():.2e} cm2/g")
        print(f"    kappa_bf  range: {kappa_bf.min():.2e} -- {kappa_bf.max():.2e} cm2/g")
        print(f"    kappa_H-  range: {kappa_Hm.min():.2e} -- {kappa_Hm.max():.2e} cm2/g")
        print(f"    kappa_abs/kappa_es range: {frac_2d.min():.2e} -- {frac_2d.max():.2e}")
        kappa_abs_2d = kappa_es * frac_2d
        logkabs_2d = np.log(np.maximum(kappa_abs_2d * urho * ur, 1e-30)).astype(np.float32)
        logkabs_4d = logkabs_2d.T[np.newaxis, :, :, np.newaxis]
        hpbanab1 = np.broadcast_to(logkabs_4d,
                                    (nfreq, ntab, ntab, nim)).copy().astype(np.float32)
        for (ir_ref, iT_ref, label) in [
                (min(6,ntab-1),  min(3,ntab-1),  'photosphere (rho~1e-11, T~1e4K)'),
                (min(10,ntab-1), min(1,ntab-1),  'cool shell  (rho~1e-7,  T~4e3K)')]:
            frac_ref = frac_2d[ir_ref, iT_ref]
            print(f"  {label}: kappa_abs/kappa_es = {frac_ref:.3e}  "
                  f"(kappa_abs = {kappa_es*frac_ref:.3e} cm^2/g)")

    # --- Time grid ---
    # daybegl=0, dayendl=2 in ronfitab.trf
    # time(nw) = 10^(daybegl + (dayendl-daybegl)/(ns-1) * (nw-1))
    daybegl = 0.0
    dayendl = 2.0
    times = [10**(daybegl + (dayendl - daybegl) / (ns - 1) * (nw - 1))
             for nw in range(1, ns + 1)]

    os.makedirs(vladsf_dir, exist_ok=True)

    # Write files .1 through .6
    for nw, stime in enumerate(times, start=1):
        fname = os.path.join(vladsf_dir, f"{model_name}.{nw}")
        print(f"  Writing {fname}  (Stime={stime:.3f} days)")
        with open(fname, 'wb') as f:
            write_fortran_record(f,
                np.int32(nw),
                np.float64(stime),
                np.int32(nfreq),
                np.int32(ns),      # Msta
                np.int32(ntab),    # Nrho
                np.int32(ntab),    # NTp
                wavel.astype(np.float64),       # Wavel(Nfreq) = 100 elements
                freq_arr.astype(np.float64),    # Freq(Nfreq+1) = 101 elements
                freqmn_arr.astype(np.float64),  # Freqmn(Nfreq) = 100 elements
                TpTab,
                RhoTab,
                hpbansc1
            )

    # Write .ab file (absorption opacity)
    fname_ab = os.path.join(vladsf_dir, f"{model_name}.ab")
    print(f"  Writing {fname_ab}")
    with open(fname_ab, 'wb') as f:
        write_fortran_record(f, hpbanab1)

    print(f"Done. Tables written to {vladsf_dir}/")
    total_size = sum(
        os.path.getsize(os.path.join(vladsf_dir, f"{model_name}.{i}"))
        for i in range(1, ns + 1)
    )
    print(f"  Total size of .1-.{ns} files: {total_size/1024:.1f} KB")

    if plot_dir:
        print(f"  Generating diagnostic opacity plots in {plot_dir}")
        make_opacity_plots(
            T_grid, rho_grid, hpbanab1, hpbansc1,
            model_name=model_name, plot_dir=plot_dir,
            kappa_es=kappa_es, urho=urho, ur=ur, nim=nim)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('vladsf_dir', help='Path to run/vladsf/ directory')
    parser.add_argument('model_name', help='Model name (e.g. test)')
    parser.add_argument('--nzon', type=int, default=200, help='Number of STELLA zones (Mzon)')
    parser.add_argument('--kappa', type=float, default=0.34,
                        help='Electron scattering opacity in cm^2/g (default 0.34)')
    parser.add_argument('--nfreq', type=int, default=100,
                        help='Number of frequency bins (NFREQ = @Nfreq+4 = 100)')
    parser.add_argument('--ntab', type=int, default=14,
                        help='Number of T/rho table points (@Ntab = 14)')
    parser.add_argument('--ns', type=int, default=6,
                        help='Number of time slices (@Ns = 6)')
    parser.add_argument('--skip', type=int, default=None,
                        help='Zone skip (@skip, default=Mzon/2)')
    parser.add_argument('--abs-fraction', type=float, default=None,
                        help='Constant kappa_abs/kappa_es override (default: use Kramers '
                             'free-free computed on T/rho table grid). '
                             'Set to a fixed value (e.g. 0.01) to revert to old behaviour.')
    parser.add_argument('--abs-fraction-floor', type=float, default=1e-4,
                        help='Minimum kappa_abs/kappa_es when using Kramers mode (default 1e-4). '
                             'Prevents log(0) at high-T/low-rho extremes.')
    parser.add_argument('--fh', type=float, default=0.70,
                        help='H mass fraction for bound-free + H- opacity (default 0.70).')
    parser.add_argument('--fhe', type=float, default=0.28,
                        help='He mass fraction for free-free + Saha correction (default 0.28).')
    parser.add_argument('--opal-tables', action='store_true', default=False,
                        help='Use OPAL Rosseland mean tables with per-zone composition '
                             '(requires --opal-table-dir and --input-profile).')
    parser.add_argument('--opal-table-dir', type=str, default=None,
                        help='Directory containing the tables/ subdirectory with '
                             'codataa..codatae OPAL table files.')
    parser.add_argument('--input-profile', type=str, default=None,
                        help='Path to STELLA input_profile with per-cell composition '
                             '(columns: r rho p u fH fHe fC fO fSi fFe).')
    parser.add_argument('--plot-dir', type=str, default=None,
                        help='Directory to save diagnostic opacity plots. '
                             'If omitted, no plots are generated.')
    parser.add_argument('--stella-root', type=str, default=None,
                        help='Path to OpenStella root directory. Used to auto-detect '
                             '@skip from opacityHomo.inc so nim matches the compiled binary.')
    args = parser.parse_args()

    print(f"Generating STELLA opacity tables:")
    print(f"  vladsf_dir: {args.vladsf_dir}")
    print(f"  model_name: {args.model_name}")
    print(f"  Nzon={args.nzon}, Nfreq={args.nfreq}, Ntab={args.ntab}, Ns={args.ns}")
    print(f"  kappa_es={args.kappa} cm^2/g, fH={args.fh}, fHe={args.fhe}")
    if args.opal_tables:
        print(f"  absorption model: OPAL (per-zone composition), floor={args.abs_fraction_floor:.0e}")
    elif args.abs_fraction is not None:
        print(f"  abs_fraction={args.abs_fraction:.2e} (constant override)")
    else:
        print(f"  absorption model: physical (ff + H bf + H-), floor={args.abs_fraction_floor:.0e}")

    generate_opus_tables(
        args.vladsf_dir,
        args.model_name,
        nzon=args.nzon,
        kappa_es=args.kappa,
        nfreq=args.nfreq,
        ntab=args.ntab,
        ns=args.ns,
        skip=args.skip,
        abs_fraction=args.abs_fraction,
        abs_fraction_floor=args.abs_fraction_floor,
        fH=args.fh,
        fHe=args.fhe,
        opal_tables=args.opal_tables,
        opal_table_dir=args.opal_table_dir,
        input_profile=args.input_profile,
        plot_dir=args.plot_dir,
        stella_root=args.stella_root,
    )


if __name__ == '__main__':
    main()
