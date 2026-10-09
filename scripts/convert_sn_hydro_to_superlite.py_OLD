#!/usr/bin/env python3
"""Convert CMFGEN SN_HYDRO_DATA to SuperLite input.str format.

SuperLite input.str format (1D spherical):
  Line 1: comment
  Line 2: #  nx ny nz ncol nabund
  Line 3: column labels
  Data: one row per cell (outer to inner for CMFGEN, inner to outer for SuperLite)

SuperLite expects:
  x_left, x_right, vx_left, vx_right, mass, rad_temp, temp, [abundances...]

SuperLite uses cgs units for radius (cm) and velocity (cm/s).
Mass is in grams per cell.
Temperature in Kelvin.
Abundances are mass fractions.
"""

import argparse
import numpy as np
import sys
import os


def read_sn_hydro_data(filepath):
    """Read CMFGEN SN_HYDRO_DATA file."""
    with open(filepath, 'r') as f:
        lines = f.readlines()

    header = {}
    for line in lines:
        if line.startswith('Number of data points'):
            header['nd'] = int(line.split(':')[1].strip())
        elif line.startswith('Number of mass fractions'):
            header['nmf'] = int(line.split(':')[1].strip())
        elif line.startswith('Time(days)'):
            header['t_days'] = float(line.split(':')[1].strip())
        elif line.startswith('Energy Lum'):
            header['lum_lsun'] = float(line.split(':')[1].strip())
        elif line.startswith('Inner boundary radius'):
            header['r_inner'] = float(line.split(':')[1].strip())
        elif line.startswith('Inner boundary temperature'):
            header['t_inner'] = float(line.split(':')[1].strip())

    nd = header['nd']
    
    def read_section(keyword, n_per_line=6):
        """Read a data section starting with keyword line."""
        for i, line in enumerate(lines):
            if line.strip().startswith(keyword):
                data = []
                j = i + 1
                while len(data) < nd:
                    parts = lines[j].split()
                    for p in parts:
                        data.append(float(p))
                    j += 1
                return np.array(data[:nd])
        raise ValueError(f"Section '{keyword}' not found")

    # Read sections (each has nd values, 6 per line)
    r = read_section('Radius grid') * 1e10  # 10^10 cm -> cm
    v = read_section('Velocity') * 1e5      # km/s -> cm/s
    sigma = read_section('Sigma')
    T = read_section('Temperature') * 1e4   # 10^4 K -> K
    rho = read_section('Density')           # g/cm^3

    # Read kappa (Rosseland opacity) if present
    try:
        kappa = read_section('Kappa')       # cm^2/g
    except ValueError:
        kappa = None
    
    # Read mass fractions
    mf_labels = ['HYD', 'HE', 'CARB', 'OXY', 'SIL', 'CAL', 'NICK', 'COB', 'IRON']
    mass_fracs = {}
    for label in mf_labels:
        mass_fracs[label] = read_section(f'{label} mass fraction')
    
    # Also read Ni56, Co56, Fe56 if present
    for label in ['NICK 56', 'COB 56', 'IRON 56']:
        try:
            mass_fracs[label] = read_section(f'{label} mass fraction')
        except ValueError:
            pass

    return header, r, v, sigma, T, rho, mass_fracs, kappa


def convert_to_superlite(header, r, v, sigma, T, rho, mass_fracs, output_dir,
                          t_snap=None, truncate_tau=None, nomp=4, kappa=None):
    """Convert to SuperLite input.str format.
    
    CMFGEN indexes outer-to-inner (r decreasing with index).
    SuperLite expects inner-to-outer (r increasing with index).
    
    SuperLite 1D spherical geometry (igeom=11):
      - x_left, x_right: cell boundaries in cm
      - vx_left, vx_right: velocities at boundaries in cm/s
      - mass: cell mass in grams
      - rad_temp: radiation temperature (optional, use T if not available)
      - temp: electron temperature in K
      - abundances: mass fractions per element
    """
    nd = len(r)
    
    # Reverse from CMFGEN (outer-to-inner) to SuperLite (inner-to-outer)
    r_sl = r[::-1].copy()
    v_sl = v[::-1].copy()
    T_sl = T[::-1].copy()
    rho_sl = rho[::-1].copy()
    
    # Reverse mass fractions
    mf_sl = {}
    for key, val in mass_fracs.items():
        mf_sl[key] = val[::-1].copy()
    
    # Compute cell boundaries (x_left, x_right)
    # r(i) is the cell center; we need boundaries
    # Use geometric mean of adjacent cell centers for boundaries
    x_left = np.zeros(nd + 1)
    x_right = np.zeros(nd + 1)
    vx_left = np.zeros(nd + 1)
    vx_right = np.zeros(nd + 1)
    
    # Inner boundary
    x_left[0] = r_sl[0] * (r_sl[0] / r_sl[1])  # extrapolate inner edge
    vx_left[0] = v_sl[0] * (v_sl[0] / v_sl[1]) if v_sl[1] > 0 else v_sl[0]
    
    # Cell boundaries between cells
    for i in range(nd - 1):
        x_right[i] = np.sqrt(r_sl[i] * r_sl[i + 1])
        vx_right[i] = 0.5 * (v_sl[i] + v_sl[i + 1])
        x_left[i + 1] = x_right[i]
        vx_left[i + 1] = vx_right[i]
    
    # Outer boundary
    x_right[nd - 1] = r_sl[nd - 1] * (r_sl[nd - 1] / r_sl[nd - 2])
    vx_right[nd - 1] = v_sl[nd - 1] * (v_sl[nd - 1] / v_sl[nd - 2])
    
    # Compute cell masses: m = rho * 4/3 * pi * (r_right^3 - r_left^3)
    vol = (4.0 / 3.0) * np.pi * (x_right[:nd]**3 - x_left[:nd]**3)
    mass = rho_sl * vol
    
    # Truncate if requested (by Rosseland optical depth)
    if truncate_tau is not None:
        # Use actual Rosseland kappa from SN_HYDRO_DATA if available;
        # fall back to grey electron-scattering opacity (0.34 cm^2/g)
        if kappa is not None:
            kappa_sl = kappa[::-1].copy()  # reverse to inner-to-outer
            print(f"  Using Rosseland kappa from SN_HYDRO_DATA: {kappa_sl.min():.3e} - {kappa_sl.max():.3e} cm^2/g")
        else:
            kappa_sl = np.full(nd, 0.34)
            print(f"  No kappa in SN_HYDRO_DATA, using grey kappa=0.34 cm^2/g")
        tau = np.cumsum(kappa_sl[:nd] * rho_sl * (x_right[:nd] - x_left[:nd]))
        mask = tau <= truncate_tau
        if mask.sum() < 10:
            print(f"Warning: truncation at tau={truncate_tau} leaves only {mask.sum()} cells, keeping all")
            mask = np.ones(nd, dtype=bool)
        nd_use = mask.sum()
        x_left = x_left[:nd_use + 1]
        x_right = x_right[:nd_use + 1]
        vx_left = vx_left[:nd_use + 1]
        vx_right = vx_right[:nd_use + 1]
        mass = mass[:nd_use]
        T_sl = T_sl[:nd_use]
        vol = vol[:nd_use]
        for key in mf_sl:
            mf_sl[key] = mf_sl[key][:nd_use]
        print(f"Truncated to {nd_use} cells (tau < {truncate_tau})")
    else:
        nd_use = nd
    
    # Map mass fractions to SuperLite element names
    # SuperLite uses element symbols: h, he, c, n, o, ne, na, mg, si, s, ar, ca, fe, ni
    # CMFGEN SN_HYDRO_DATA has: HYD, HE, CARB, OXY, SIL, CAL, NICK, COB, IRON
    # SIL = Si-group (Si, S, Ar, Ca combined in STELLA)
    # CAL = Ca (sometimes separate)
    # We'll map: H->h, He->he, C->c, O->o, Si-group->si, Ca->ca, Ni->ni, Co->co, Fe->fe
    
    # For SuperLite, we need individual elements.
    # Our SN_HYDRO_DATA has combined groups. Let's split them:
    # SIL = Si + S + Ar (STELLA groups them)
    # CAL = Ca
    # NICK = Ni (stable)
    # COB = Co (stable)
    # IRON = Fe (stable)
    
    # For simplicity, map SIL to si, CAL to ca
    # SuperLite's Kurucz data has all these elements
    
    abund_labels = ['h', 'he', 'c', 'o', 'si', 'ca', 'ni', 'co', 'fe']
    abund_data = [
        mf_sl.get('HYD', np.zeros(nd_use)),
        mf_sl.get('HE', np.zeros(nd_use)),
        mf_sl.get('CARB', np.zeros(nd_use)),
        mf_sl.get('OXY', np.zeros(nd_use)),
        mf_sl.get('SIL', np.zeros(nd_use)),   # Si-group
        mf_sl.get('CAL', np.zeros(nd_use)),    # Ca
        mf_sl.get('NICK', np.zeros(nd_use)),   # Ni
        mf_sl.get('COB', np.zeros(nd_use)),    # Co
        mf_sl.get('IRON', np.zeros(nd_use)),   # Fe
    ]
    
    # Normalize abundances to sum to 1 (in case of rounding)
    total = sum(abund_data)
    for i in range(len(abund_data)):
        abund_data[i] = abund_data[i] / np.maximum(total, 1e-30)
    
    nabund = len(abund_labels)
    
    # Write input.str
    os.makedirs(output_dir, exist_ok=True)
    str_file = os.path.join(output_dir, 'input.str')
    
    with open(str_file, 'w') as f:
        t_label = header.get('t_days', 0.0)
        f.write(f'# 1D spherical - SN snapshot (t={t_label:.1f}d)\n')
        ncol = 7 + nabund  # x_left, x_right, vx_left, vx_right, mass, rad_temp, temp + abundances
        f.write(f'#      {nd_use}          1 1 {ncol} {nabund}\n')
        
        # Column labels
        labels = ['x_left', 'x_right', 'vx_left', 'vx_right', 'mass', 'rad_temp', 'temp']
        labels += abund_labels
        f.write('#   ' + '  '.join(f'{l:>12s}' for l in labels) + '\n')
        
        # Data: one row per cell
        for i in range(nd_use):
            vals = [x_left[i], x_right[i], vx_left[i], vx_right[i], mass[i], T_sl[i], T_sl[i]]
            vals += [abund_data[j][i] for j in range(nabund)]
            f.write('  '.join(f'{v:14.6e}' for v in vals) + '  \n')
    
    print(f"Wrote {str_file} with {nd_use} cells, {nabund} abundance columns")
    
    # Write input.par from reference template
    t_days = t_snap if t_snap else header.get('t_days', 30.0)
    lum_lsun = header.get('lum_lsun')
    if lum_lsun is None:
        # Ejecta-only files have no Energy Lum header. Estimate the emergent
        # luminosity as the blackbody at the tau=2/3 surface (kappa*rho
        # integrated inward), else the inner-boundary blackbody.
        # r/T are outer-to-inner; r already in cm, T in K.
        if kappa is not None and np.any(kappa > 0):
            tau = np.zeros(nd)
            for i in range(1, nd):
                tau[i] = tau[i - 1] + 0.5 * (
                    kappa[i] * rho[i] + kappa[i - 1] * rho[i - 1]
                ) * (r[i - 1] - r[i])
            if tau[-1] >= 2.0 / 3.0:
                i_ph = min(int(np.searchsorted(tau, 2.0 / 3.0)), nd - 1)
                lum_lsun = (
                    4.0 * np.pi * r[i_ph] ** 2 * 5.6704e-5 * T[i_ph] ** 4
                    / 3.828e33
                )
        if lum_lsun is None:
            lum_lsun = 4.0 * np.pi * r[-1] ** 2 * 5.6704e-5 * T[-1] ** 4 / 3.828e33
        print(f"No Energy Lum header - using ejecta luminosity estimate: "
              f"{lum_lsun:.3e} Lsun")
    lum_erg = lum_lsun * 3.828e33  # lsun to erg/s

    # Model name from SN_HYDRO_DATA or default
    model_name = header.get('model_name', 'SN')

    # Compute photospheric radius (grey tau=1) for nebular approximation
    # SN_HYDRO_DATA is outer-to-inner (r decreasing). Integrate from outside in.
    # Use mean kappa from the data, or fallback to 0.4 cm^2/g.
    # NOTE: r is already in cm (converted from 10^10 cm units at read time,
    # see read_sn_hydro_data line 62). Do NOT multiply by 1e10 again.
    kappa_mean = float(np.mean(kappa[kappa > 0])) if np.any(kappa > 0) else 0.4
    r_cm = r  # already in cm
    tau_cum = 0.0
    r_phot = r_cm[0]  # default: outermost (SN_HYDRO_DATA is outer-to-inner)
    for i in range(len(r_cm) - 1, 0, -1):
        dr = r_cm[i - 1] - r_cm[i]
        tau_cum += kappa_mean * rho[i] * dr
        if tau_cum >= 1.0:
            r_phot = r_cm[i]
            break
    print(f"Photospheric radius (grey tau=1, kappa={kappa_mean:.3f}): {r_phot:.4e} cm")

    # Read the reference template and substitute placeholders
    ref_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                            'superlite_base', 'input_parameters_reference.par')
    if os.path.exists(ref_path):
        with open(ref_path) as f:
            template = f.read()
        # Strip comment lines (everything before &inputpars)
        idx = template.index('&inputpars')
        par_template = template[idx:]
        par_content = (par_template
                        .replace('{MODEL_NAME}', model_name)
                        .replace('{EPOCH}', f'{t_days:.1f}d')
                        .replace('{NDIM}', str(nd_use))
                        .replace('{L_BOL}', f'{lum_erg:.4e}')
                        .replace('{R_PHOT}', f'{r_phot:.4e}')
                        .replace('{NOMP}', str(nomp)))
    else:
        # Fallback: hardcoded optical-only defaults
        par_content = f"""&inputpars

in_name = '{model_name}'
in_comment = 'Production nebular, 1-D sph model of {model_name} at {t_days:.1f}d'

!-- misc (in_nomp = OpenMP threads, must match --cpus-per-task in SLURM)
 in_nomp = {nomp}
 in_io_grabstdout = f
 in_io_opacdump = 'one'
 in_io_profdump = t

!-- grid settings
 in_ndim = {nd_use}, 1, 1

!-- group settings (WIDE: 1-30000 Angstrom for correct UV transport)
 in_grp_ng = 6000
 in_grp_wlmin = 1d-8
 in_grp_wlmax = 30d-5

!-- flux settings (OPTICAL+NIR: 3000-20000 Angstrom, ~5 A/bin for good S/N)
!-- Extended beyond 10000 A to capture the true spectral peak in cool,
!-- interacting-supernova models where the color temperature may be < 3000 K.
 in_flx_ndim = 3500, 1, 1
 in_flx_wlmin = 3d-5
 in_flx_wlmax = 2d-4

!-- particle numbers (2^26 = 67M particles for high S/N: ~19k packets/bin)
 in_src_n2s = 26
 in_L_bol = {lum_erg:.4e}

!-- ddmc mfp per cell
 in_trn_tauddmc = 3d0
 in_taulump = 3d0

!-- EOS / ionization: nebular approximation (Mazzali & Lucy 1993)
 in_nebular = t
 in_R_phot = {r_phot:.4e}    ! grey tau=1 radius from SN_HYDRO_DATA

!-- NLTE enabled: iterative radiation-temperature feedback.
!-- in_nlte=t with in_niter=3 allows the radiation field from the hot shock
!-- region to heat and ionize the outer CSM/ejecta during iteration. This is
!-- essential for interacting supernovae where the ionization front is
!-- maintained by the radiation field, not by the gas temperature alone.
 in_nlte = t
 in_nlte_nelem = 1
 in_niter = 3

!-- opacity settings (all physical opacities enabled)
 in_nobbopac = f
 in_nobfopac = f
 in_noffopac = f
 in_nothmson = f

/"""
    
    par_file = os.path.join(output_dir, 'input.par')

    # Always regenerate input.par from the reference template to ensure
    # physics settings (NLTE, flux range, nebular, opacities) are kept
    # up-to-date with the optimum configuration.  Only preserve a small
    # set of user-tunable performance parameters from any existing file.
    if os.path.exists(par_file):
        import re
        with open(par_file) as f:
            existing = f.read()
        # Extract user-tunable performance parameters from the existing file
        preserve_keys = ['in_src_n2s', 'in_grp_ng', 'in_trn_tauddmc', 'in_taulump']
        for key in preserve_keys:
            m = re.search(rf'{key}\s*=\s*(\S+)', existing)
            if m:
                val = m.group(1)
                par_content = re.sub(rf'{key}\s*=\s*\S+', f'{key} = {val}', par_content)
        print(f"Regenerated input.par from template (preserved performance params)")
    else:
        print(f"Wrote {par_file}")

    with open(par_file, 'w') as f:
        f.write(par_content)

    return nd_use


def main():
    parser = argparse.ArgumentParser(description='Convert SN_HYDRO_DATA to SuperLite input.str')
    parser.add_argument('sn_hydro_data', help='Path to SN_HYDRO_DATA file')
    parser.add_argument('-o', '--output-dir', required=True, help='Output directory for SuperLite input files')
    parser.add_argument('--truncate-tau', type=float, default=None, help='Truncate at this Rosseland optical depth (using kappa from SN_HYDRO_DATA)')
    parser.add_argument('--nomp', type=int, default=4, help='Number of OpenMP threads for SuperLite (in_nomp)')
    args = parser.parse_args()
    
    header, r, v, sigma, T, rho, mass_fracs, kappa = read_sn_hydro_data(args.sn_hydro_data)
    
    print(f"Read SN_HYDRO_DATA: {header['nd']} points, t={header['t_days']:.1f}d")
    print(f"  R range: {r.min():.3e} - {r.max():.3e} cm")
    print(f"  V range: {v.min():.3e} - {v.max():.3e} cm/s")
    print(f"  T range: {T.min():.3e} - {T.max():.3e} K")
    print(f"  rho range: {rho.min():.3e} - {rho.max():.3e} g/cm^3")
    
    nd = convert_to_superlite(header, r, v, sigma, T, rho, mass_fracs,
                               args.output_dir, truncate_tau=args.truncate_tau,
                               nomp=args.nomp, kappa=kappa)
    
    print(f"\nConversion complete. {nd} cells written to {args.output_dir}/")


if __name__ == '__main__':
    main()
