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
  4. Zeroes the pre-collapse infall velocities (STELLA uses thermal bomb)
  5. Sets an initial temperature floor (STELLA radiation floor)
  6. Maps composition to STELLA 6-element scheme: H, He, C, O, Si, FePeak
  7. Writes .hyd, .abn, .xni, .dat, and input_profile files

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
    """Map Kepler isotopes to STELLA 6-element (H, He, C, O, Si, FePeak) scheme.

    STELLA composition columns:
        H, He, C, O, Si, FePeak
    where FePeak absorbs Fe54 + Fe56 + Ti44 + Cr48 + Ca40 + S32 + Ar36
    and Si absorbs Si28 + Mg24 + Ne20 (Si-group)
    and O  absorbs O16
    and C  absorbs C12 + N14
    and He absorbs He4
    and H  = 0 (He star has no H)

    Normalised to sum = 1 per zone.
    """
    H_frac   = np.zeros(len(kep['rho']))                 # no H in He star
    He_frac  = kep['He4'].copy()
    C_frac   = kep['C12'] + kep['N14']
    O_frac   = kep['O16']
    Si_frac  = kep['Si28'] + kep['Ne20'] + kep['Mg24']
    Fe_frac  = (kep['Fe56'] + kep['Fe_tot'] +
                kep['S32'] + kep['Ar36'] + kep['Ca40'] +
                kep['Ti44'])

    total = H_frac + He_frac + C_frac + O_frac + Si_frac + Fe_frac
    # Normalise; guard against zero-composition zones (innermost NSE)
    bad = total < 1e-10
    total = np.where(bad, 1.0, total)
    # In the innermost NSE-dominated zones all Fe; set FePeak=1 there
    Fe_frac = np.where(bad, 1.0, Fe_frac)

    He_frac  /= total; C_frac  /= total; O_frac /= total
    Si_frac  /= total; Fe_frac /= total

    return H_frac, He_frac, C_frac, O_frac, Si_frac, Fe_frac


# ── Remap onto uniform STELLA grid ───────────────────────────────────────────
def remap_to_stella_grid(kep, mass_cut_msun, nzon, t_floor_K=3000.0, verbose=True):
    """Excise inner compact remnant and remap onto nzon log-spaced STELLA zones.

    Mass coordinates above mass_cut_msun are kept.
    The new grid is log-spaced in radius from r_cut to r_max.

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

    if verbose:
        print(f"\nMass cut: {mass_cut_msun:.3f} Msun (zone {idx_cut+1}/{len(kep['m_g'])})")
        print(f"  r_cut = {r_cut:.3e} cm")
        print(f"  Ejecta mass above cut: {(kep['m_g'][-1] - m_cut_g)/MSUN:.4f} Msun")

    # Source arrays (ejecta only, zones idx_cut : end)
    r_src   = kep['r_cm'][idx_cut:]
    rho_src = kep['rho'][idx_cut:]
    T_src   = kep['T_K'][idx_cut:]
    v_src   = kep['v_cms'][idx_cut:]   # will be zeroed for bomb mode

    # Get compositions
    fH, fHe, fC, fO, fSi, fFe = build_stella_composition(kep)
    fH_src  = fH[idx_cut:];  fHe_src = fHe[idx_cut:]
    fC_src  = fC[idx_cut:];  fO_src  = fO[idx_cut:]
    fSi_src = fSi[idx_cut:]; fFe_src = fFe[idx_cut:]

    # New STELLA grid: log-spaced in radius
    r_new = np.logspace(np.log10(r_cut), np.log10(r_max), nzon)

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
    # Zero pre-SN infall velocities — STELLA thermal bomb sets explosion velocities
    v_new   = np.zeros(nzon)

    fH_new  = interp(fH_src);  fHe_new = interp(fHe_src)
    fC_new  = interp(fC_src);  fO_new  = interp(fO_src)
    fSi_new = interp(fSi_src); fFe_new = interp(fFe_src)

    # Re-normalise composition after interpolation
    comp_sum = fH_new + fHe_new + fC_new + fO_new + fSi_new + fFe_new
    comp_sum = np.maximum(comp_sum, 1e-10)
    fH_new  /= comp_sum; fHe_new /= comp_sum; fC_new  /= comp_sum
    fO_new  /= comp_sum; fSi_new /= comp_sum; fFe_new /= comp_sum

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
        for name, f_arr in [('He', fHe_new),('C', fC_new),('O', fO_new),
                             ('Si', fSi_new),('FePeak', fFe_new)]:
            print(f"  <{name}>_m = {(f_arr * dm).sum() / M_tot:.4f}")

    return {
        'r': r_new, 'rho': rho_new, 'T': T_new, 'v': v_new,
        'fH': fH_new, 'fHe': fHe_new, 'fC': fC_new,
        'fO': fO_new, 'fSi': fSi_new, 'fFe': fFe_new,
        'r_cut': r_cut, 'mass_cut_g': m_cut_g,
    }


# ── .hyd writer (matches existing STELLA format) ─────────────────────────────
def write_hyd(fname, g, mass_cut_msun):
    """Write STELLA *.hyd initial conditions file."""
    r = g['r']; rho = g['rho']; T = g['T']; v = g['v']
    nzon = len(r)

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
def write_abn(fname, g):
    """Write STELLA *.abn composition file.

    Columns: km  0  0  0  H  He  C  N  O  Ne  Na  Mg  Al  Si  S  Ar  Ca  FePeak  Ni58  Ni56
    Ni56 is set to 0 here - radioactive Ni56 comes from the .xni file only.
    """
    nzon = len(g['r'])
    with open(fname, 'w') as f:
        for i in range(nzon):
            f.write(f"{i+1:4d}{0.0:10.3E}{0.0:10.3E}{0.0:10.3E}"
                    f"{g['fH'][i]:10.3E}{g['fHe'][i]:10.3E}{g['fC'][i]:10.3E}"
                    f"{0.0:10.3E}{g['fO'][i]:10.3E}"    # N=0
                    f"{0.0:10.3E}{0.0:10.3E}{0.0:10.3E}{0.0:10.3E}"  # Ne Na Mg Al
                    f"{g['fSi'][i]:10.3E}"               # Si
                    f"{0.0:10.3E}{0.0:10.3E}{0.0:10.3E}"  # S Ar Ca
                    f"{g['fFe'][i]:10.3E}"               # FePeak
                    f"{0.0:10.3E}{0.0:10.3E}\n")         # Ni58 Ni56 (both 0)
    print(f"Written {fname}  ({nzon} zones, Ni56=0 - radioactive Ni56 from .xni file only)")


# ── .xni writer (Ni56 distribution) ──────────────────────────────────────────
def write_xni(fname, g, ni56_mass_msun, ni56_dist_msun, nzon):
    """Write STELLA binary .xni file with Ni56 mass fractions.

    Ni56 is placed in the innermost zones covering ni56_dist_msun of ejecta.
    The xni format: 2 header floats + nzon*(fNi56, placeholder) pairs.
    """
    import struct

    r = g['r']; rho = g['rho']
    r_lo = np.concatenate([[g['r_cut']], r[:-1]])
    dm   = np.maximum((4/3)*np.pi*(r**3 - r_lo**3)*rho, 0.0)
    m_cum = np.cumsum(dm) / MSUN   # cumulative ejecta mass in Msun

    # Find how many inner zones cover ni56_dist_msun
    ni56_zones = np.where(m_cum <= ni56_dist_msun)[0]
    if len(ni56_zones) == 0:
        ni56_zones = np.array([0])

    n_ni_zones = len(ni56_zones)
    m_ni56_region = m_cum[ni56_zones[-1]]  # actual mass in those zones

    # Ni56 mass fraction within the Ni56 region
    fni56 = ni56_mass_msun / max(m_ni56_region, 1e-10)
    fni56 = min(fni56, 1.0)

    print(f"\n.xni Ni56 distribution:")
    print(f"  Total Ni56: {ni56_mass_msun:.3f} Msun")
    print(f"  Ni56 region: {n_ni_zones} inner zones covering {m_ni56_region:.3f} Msun ejecta")
    print(f"  Ni56 mass fraction in those zones: {fni56:.4f}")

    # Build float array: 2 header + nzon pairs (value, placeholder)
    PLACEHOLDER = -1.417409e-34
    floats = [5.605194e-42, PLACEHOLDER]  # header
    for i in range(nzon):
        if i < n_ni_zones:
            floats.append(float(fni56))
        else:
            floats.append(0.0)
        floats.append(PLACEHOLDER)

    data = struct.pack('f' * len(floats), *floats)
    with open(fname, 'wb') as f:
        f.write(data)
    print(f"Written {fname}  ({len(data)} bytes)")


# ── input_profile writer (HERACLES-format, needed by opacity table generator) ─
def write_input_profile(fname, g):
    """Write a HERACLES-format input_profile from remapped grid data.

    Columns: r_cm  rho  p_cgs  v_cms  fH  fHe  fC  fO  fSi  fFe
    This is the format consumed by generate_stella_opacity_tables.py
    when --opal-tables and --input-profile are requested.
    Pressure is reconstructed from ideal gas + radiation: p = rho*kT/(mu*mH) + a*T^4/3.
    """
    r = g['r']; rho = g['rho']; T = g['T']

    # Mean molecular weight for fully ionised He/CO mix (approximate)
    # mu ~ 4/(2+3*fHe + 2*fC + 2*fO + 2*fSi + 2*fFe) for full ionisation
    # For a He-dominated envelope mu ~ 4/3 ~ 1.33; use 2.0 as safe floor
    fHe = g['fHe']; fC = g['fC']; fO = g['fO']
    mu = 4.0 / np.maximum(2.0 + 3*fHe + 2*(fC + fO + g['fSi'] + g['fFe']), 0.5)
    p_gas = rho * K_B * T / (mu * M_H)
    p_rad = A_RAD * T**4 / 3.0
    p = p_gas + p_rad

    header = ("# HERACLES-format input_profile generated from Kepler presn by "
              "convert_kepler_presn_to_stella.py\n"
              "# columns: r_cm rho_cgs p_cgs u_cgs fx_H fx_He fx_C fx_O fx_Si fx_Fe\n")
    with open(fname, 'w') as f:
        f.write(header)
        for i in range(len(r)):
            f.write(f"{r[i]:.8e} {rho[i]:.8e} {p[i]:.8e} {g['v'][i]:.8e} "
                    f"{g['fH'][i]:.8e} {fHe[i]:.8e} {fC[i]:.8e} {g['fO'][i]:.8e} "
                    f"{g['fSi'][i]:.8e} {g['fFe'][i]:.8e}\n")
    print(f"Written {fname}  ({len(r)} zones, HERACLES-format input_profile)")


# ── .dat writer (STELLA control file) ────────────────────────────────────────
def write_dat(fname, model_name, target_days,
              explosion_energy=10.0, explosion_mass=0.5, explosion_time=1.0,
              explosion_delay=0.0, nstmax=2000000,
              eps=0.003, hmax=4.32e4, meth=3,
              wacc_t=0.5, temp_floor=0.002,
              snapshot_interval=10.0, iout=1):
    """Write STELLA *.dat control file for a Type IIb explosion."""
    import math

    t_start = 0.0
    t_end   = target_days

    # Output epochs (log-spaced after first few days, then linear)
    KOMAX = 80
    n_epochs = min(KOMAX, max(10, int(target_days / snapshot_interval)))
    # First few days more densely sampled (shock breakout), then uniform
    early = [0.5, 1.0, 2.0, 5.0, 10.0, 15.0, 20.0]
    late  = np.linspace(30.0, target_days, n_epochs - len(early)).tolist()
    epochs = sorted(set(early + late))[:KOMAX]

    with open(fname, 'w') as f:
        f.write(f"   EPS   (STIFF ACCURACY)   HMIN (MIN.STEP(SEC)) HMAX (MAX.STEP(SEC))\n")
        f.write(f"   {eps:.4g}                    1e-14               {hmax:.4g}\n")
        f.write(f"   METH  ADAMS(1)-GEAR(2)-GEAR+BGH(3)   JSTART  (0-START_ELSE_ORDER) MAXORD KNadap\n")
        f.write(f"    {meth}                        0                            4        -4\n")
        f.write(f"   abs(NSTA)>abs(NSTB)=>reads last   TcurA     TcurB\n")
        f.write(f"   -9000     -1                      {t_start:.4f}   {t_end:.4f}\n")
        f.write(f"   NSTMAX  NDebug  NOUT    IOUT (0 -full <0 print zones)   MBATCH\n")
        f.write(f"   {nstmax}   100000  {len(epochs)}      {iout}                             900\n")
        f.write(f"   AMNI (CONTAM. MASS)     XMNI (56NI MASS)\n")
        f.write(f"   8.D-01                  0.0d0\n")
        f.write(f"   AMHT(Solar)  EBurst(1e50) tBurst(s)  tstart -- Heated Core, Energy & time\n")
        f.write(f"   {explosion_mass:.4g}d0      {explosion_energy:.4g}d0      {explosion_time:.4g}d0      {explosion_delay:.4g}d0\n")
        f.write(f"   EKO (kin.energy,1e50) mass fract.tri. u profile  us (out +1)\n")
        f.write(f"   0d0                      0.d0  9.d-1  1.d0          +1.\n")
        f.write(f"   THRMAT   CRAP     CONV     EDTM     CHNCND  Givdtl\n")
        f.write(f"   1.D-30   1.D0     F        T        T       T\n")
        f.write(f"   FLOOR(R) FLOOR(V) FLOOR(T) FLOOR(RADIAT)    FLOOR(UCONV)\n")
        f.write(f"   1.D-12   1.D-03   {temp_floor:.3g}D+00   1.D-02           1.D-04\n")
        f.write(f"   Wacc(R)  Wacc(V)  Wacc(T)  Wacc(RADIAT)\n")
        f.write(f"   1d-0    1d-0    {wacc_t:.2g}D-00   1D-00\n")
        f.write(f"   FitTau   TAUTOL   Rvis   AQ    BQ       DRT     NRT     SCAT\n")
        f.write(f"   15.0D0   1.1D0     1.D0   2.D0  1.D+00   1.      1       T\n")
        f.write(f"    NTO then    DATA TO\n")
        f.write(f"    {len(epochs)}\n")
        for ep in epochs:
            f.write(f"     {ep:.4f}\n")

    print(f"Written {fname}  ({len(epochs)} output epochs, target={target_days} days)")
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
                        help='EBurst in units of 1e50 erg (default: 10 = 1e51 erg)')
    parser.add_argument('--explosion-mass',   type=float, default=0.5,
                        help='AMHT heated core mass in Msun (default: 0.5)')
    parser.add_argument('--explosion-time',   type=float, default=1.0,
                        help='tBurst energy injection duration in seconds (default: 1)')
    parser.add_argument('--ni56-mass',    type=float, default=0.1,
                        help='Total Ni56 mass in Msun (default: 0.1)')
    parser.add_argument('--ni56-dist',    type=float, default=0.2,
                        help='Ejecta mass (Msun) over which Ni56 is distributed (default: 0.2)')
    parser.add_argument('--t-floor',      type=float, default=3000.0,
                        help='Initial temperature floor in K (default: 3000)')
    parser.add_argument('--outdir',       type=str,   default='.',
                        help='Output directory for all output files (default: .)')
    parser.add_argument('--write-dat',    action='store_true', default=True,
                        help='Also write .dat STELLA control file (default: True)')
    parser.add_argument('--write-input-profile', action='store_true', default=True,
                        help='Write HERACLES-format input_profile for opacity tables (default: True)')
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # 1. Read Kepler file
    print(f"\n=== Reading Kepler presn file: {args.presn_file} ===")
    kep = read_kepler_presn(args.presn_file)

    # 2. Remap onto STELLA grid
    print(f"\n=== Remapping onto {args.nzon} STELLA zones ===")
    g = remap_to_stella_grid(kep, args.mass_cut, args.nzon,
                             t_floor_K=args.t_floor)

    # 3. Write .hyd
    hyd_path = outdir / f"{args.model_name}.hyd"
    print(f"\n=== Writing {hyd_path} ===")
    write_hyd(str(hyd_path), g, args.mass_cut)

    # 4. Write .abn
    abn_path = outdir / f"{args.model_name}.abn"
    print(f"\n=== Writing {abn_path} ===")
    write_abn(str(abn_path), g)

    # 5. Write .xni (Ni56 distribution)
    xni_path = outdir / f"{args.model_name}.xni"
    print(f"\n=== Writing {xni_path} ===")
    write_xni(str(xni_path), g, args.ni56_mass, args.ni56_dist, args.nzon)

    # 6. Write .dat
    dat_path = outdir / f"{args.model_name}.dat"
    print(f"\n=== Writing {dat_path} ===")
    write_dat(str(dat_path), args.model_name, args.target_days,
              explosion_energy=args.explosion_energy,
              explosion_mass=args.explosion_mass,
              explosion_time=args.explosion_time)

    # 7. Write HERACLES-format input_profile (needed by opacity table generator)
    ip_path = outdir / "input_profile"
    print(f"\n=== Writing {ip_path} ===")
    write_input_profile(str(ip_path), g)

    print(f"\n=== Done. Files written to {outdir}/ ===")
    print(f"  {hyd_path.name}, {abn_path.name}, {xni_path.name}, {dat_path.name}, input_profile")


if __name__ == '__main__':
    main()
