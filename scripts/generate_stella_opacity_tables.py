#!/usr/bin/env python3
"""
Generate STELLA vladsf opacity tables for homogeneous electron scattering.

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

The .ab file contains hpbanab1 (absorption opacity), also = hpbansc1 for pure e-scattering.

Usage:
    python generate_stella_opacity_tables.py <vladsf_dir> <model_name> \\
        [--nzon NZON] [--kappa KAPPA] [--nfreq NFREQ] [--ntab NTAB] [--ns NS]
"""

import numpy as np
import struct
import os
import sys
import argparse


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


def generate_opus_tables(vladsf_dir, model_name, nzon=200, kappa_es=0.34,
                         nfreq=100, ntab=14, ns=6, skip=None,
                         wlmin=1.0, wlmax=5.0e4, ulgr=14.0):
    """
    Generate STELLA opacity table files for homogeneous electron scattering.

    Parameters
    ----------
    vladsf_dir  : directory to write tables into
    model_name  : base name for files (model_name.1 ... model_name.6, model_name.ab)
    nzon        : number of STELLA zones (Mzon)
    kappa_es    : electron scattering opacity in cm^2/g (default 0.34)
    nfreq       : number of frequency bins (NFREQ in zone.inc = @Nfreq+4 = 100)
    ntab        : number of T/rho table points (@Ntab = 14)
    ns          : number of time slices (@Ns = 6)
    skip        : zone skip for opacity table (@skip, default = Mzon/2)
    wlmin,wlmax : wavelength range in Angstrom
    ulgr        : log10(UR) where UR = STELLA length unit in cm (default 14 -> UR=1e14 cm)
    """
    if skip is None:
        skip = max(1, nzon // 2)

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
    # hpbanab1 (absorption only) — for pure electron scattering set = hpbansc1
    # (assuming all extinction is scattering; absorption = very small fraction)
    hapabs_val = np.float32(np.log(max(kappa_es * 1e-4 * urho * ur, 1e-30)))
    hpbanab1 = np.full((nfreq, ntab, ntab, nim), hapabs_val, dtype=np.float32)

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
    args = parser.parse_args()

    print(f"Generating STELLA electron scattering opacity tables:")
    print(f"  vladsf_dir: {args.vladsf_dir}")
    print(f"  model_name: {args.model_name}")
    print(f"  Nzon={args.nzon}, Nfreq={args.nfreq}, Ntab={args.ntab}, Ns={args.ns}")
    print(f"  kappa_es={args.kappa} cm^2/g")

    generate_opus_tables(
        args.vladsf_dir,
        args.model_name,
        nzon=args.nzon,
        kappa_es=args.kappa,
        nfreq=args.nfreq,
        ntab=args.ntab,
        ns=args.ns,
        skip=args.skip,
    )


if __name__ == '__main__':
    main()
