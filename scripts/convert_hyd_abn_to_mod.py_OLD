#!/usr/bin/env python3
"""
Convert STELLA .hyd + .abn files to .mod format (binary unformatted Fortran).

This bypasses the need for eve2, which requires complex trefor preprocessing.
The .mod format is defined in stradio.trf:
- First record: t, NSTEP (2 reals)
- Second record: KTAIL,NZON,NQUSED,KRAD,NCND,NFRUS,NREG (7 integers)
- Third record: XA,YA,XYZA,EKO,RADBEG,Rce,ELOST,ULGCAP,ULGEPS,ULGP,ULGV,
               ULGE,UEPRI,UTIME,UP,UPI,UE,UEI,CK1,CK2,CFR,CRAOLD,
               UM,UV,URM,DM,AM,DMIN,AMINI,DMOUT,AMOUT,YABUN (28 reals)
- Fourth record: (Y(ISAVE,1),ISAVE=1,NZON*NVARS+2*KRAD) (array of reals)
"""

import numpy as np
import struct


def read_hyd(fname: str):
    """Read STELLA .hyd file."""
    with open(fname, 'r') as f:
        # First line: timetobegin  Nzone  BM1  Rcen  rhoCen
        line1 = f.readline().split()
        t_start = float(line1[0])
        n = int(line1[1])
        bm1 = float(line1[2])
        rcen = float(line1[3])
        rhocen = float(line1[4])
        
        # Data lines: km  dMr  rHyd  rhoHyd  TpHyd  uHyd  aMr  dum
        data = []
        for _ in range(n):
            line = f.readline().split()
            km = int(line[0])
            dMr = float(line[1])
            r = float(line[2])
            rho = float(line[3])
            Tp = float(line[4])
            u = float(line[5])
            aMr = float(line[6])
            dum = float(line[7])
            data.append([km, dMr, r, rho, Tp, u, aMr, dum])
    
    return t_start, n, bm1, rcen, rhocen, np.array(data)


def read_abn(fname: str):
    """Read STELLA .abn file."""
    with open(fname, 'r') as f:
        lines = f.readlines()
    
    data = []
    for line in lines:
        parts = line.split()
        km = int(parts[0])
        dum1 = float(parts[1])
        dum2 = float(parts[2])
        dum3 = float(parts[3])
        H = float(parts[4])
        He = float(parts[5])
        C = float(parts[6])
        N = float(parts[7])
        O = float(parts[8])
        Ne = float(parts[9])
        Na = float(parts[10])
        Mg = float(parts[11])
        Al = float(parts[12])
        Si = float(parts[13])
        S = float(parts[14])
        Ar = float(parts[15])
        Ca = float(parts[16])
        Fe = float(parts[17])
        Ni58 = float(parts[18])
        Ni56 = float(parts[19])
        data.append([km, dum1, dum2, dum3, H, He, C, N, O, Ne, Na, Mg, Al, Si, S, Ar, Ca, Fe, Ni58, Ni56])
    
    return np.array(data)


def write_mod(fname: str, t_start: float, n: int, hyd_data: np.ndarray, abn_data: np.ndarray):
    """Write STELLA .mod file (binary unformatted Fortran format)."""
    
    # Constants based on STELLA source
    NVARS = 7  # Number of variables per zone (from stradio.trf)
    KRAD = 1   # Number of radiation groups
    
    # First record: t, NSTEP (2 reals)
    # t: time (s), NSTEP: step number (0 for initial model)
    record1 = np.array([t_start, 0.0], dtype=np.float64)
    
    # Second record: KTAIL,NZON,NQUSED,KRAD,NCND,NFRUS,NREG (7 integers)
    # KTAIL: tail of model, NZON: number of zones, NQUSED: number of zones used
    # KRAD: number of radiation groups, NCND: number of conditions, NFRUS: number of frequencies used, NREG: number of regions
    record2 = np.array([n, n, n, KRAD, 1, 1, 1], dtype=np.int32)
    
    # Third record: XA,YA,XYZA,EKO,RADBEG,Rce,ELOST,ULGCAP,ULGEPS,ULGP,ULGV,
    #               ULGE,UEPRI,UTIME,UP,UPI,UE,UEI,CK1,CK2,CFR,CRAOLD,
    #               UM,UV,URM,DM,AM,DMIN,AMINI,DMOUT,AMOUT,YABUN (28 reals)
    # These are mostly zero for initial models
    record3 = np.zeros(28, dtype=np.float64)
    record3[0] = 0.0  # XA
    record3[1] = 0.0  # YA
    record3[2] = 0.0  # XYZA
    record3[3] = 0.0  # EKO (kinetic energy)
    record3[4] = 0.0  # RADBEG
    record3[5] = 0.0  # Rce (core radius)
    record3[6] = 0.0  # ELOST
    record3[7] = 0.0  # ULGCAP
    record3[8] = 0.0  # ULGEPS
    record3[9] = 0.0  # ULGP
    record3[10] = 0.0  # ULGV
    record3[11] = 0.0  # ULGE
    record3[12] = 0.0  # UEPRI
    record3[13] = 0.0  # UTIME
    record3[14] = 0.0  # UP
    record3[15] = 0.0  # UPI
    record3[16] = 0.0  # UE
    record3[17] = 0.0  # UEI
    record3[18] = 0.0  # CK1
    record3[19] = 0.0  # CK2
    record3[20] = 0.0  # CFR
    record3[21] = 0.0  # CRAOLD
    record3[22] = 0.0  # UM
    record3[23] = 0.0  # UV
    record3[24] = 0.0  # URM
    record3[25] = 0.0  # DM
    record3[26] = 0.0  # AM
    record3[27] = 0.0  # YABUN
    
    # Fourth record: Y array (zone data)
    # Y(ISAVE,1) for ISAVE=1 to NZON*NVARS+2*KRAD
    # This contains the physical variables for each zone
    n_y = n * NVARS + 2 * KRAD
    record4 = np.zeros(n_y, dtype=np.float64)
    
    # Fill Y array with data from .hyd and .abn files
    # Order based on STELLA source: likely r, rho, T, u, abundances...
    for i in range(n):
        base_idx = i * NVARS
        record4[base_idx + 0] = hyd_data[i, 2]  # r
        record4[base_idx + 1] = hyd_data[i, 3]  # rho
        record4[base_idx + 2] = hyd_data[i, 4]  # T
        record4[base_idx + 3] = hyd_data[i, 5]  # u
        record4[base_idx + 4] = hyd_data[i, 1]  # dMr
        record4[base_idx + 5] = hyd_data[i, 6]  # aMr
        record4[base_idx + 6] = abn_data[i, 4]  # H abundance
    
    # Write binary file in Fortran unformatted format
    # Each record is prefixed and suffixed with record length (4 bytes)
    with open(fname, 'wb') as f:
        # Record 1
        len1 = record1.nbytes
        f.write(struct.pack('>i', len1))
        f.write(record1.tobytes())
        f.write(struct.pack('>i', len1))
        
        # Record 2
        len2 = record2.nbytes
        f.write(struct.pack('>i', len2))
        f.write(record2.tobytes())
        f.write(struct.pack('>i', len2))
        
        # Record 3
        len3 = record3.nbytes
        f.write(struct.pack('>i', len3))
        f.write(record3.tobytes())
        f.write(struct.pack('>i', len3))
        
        # Record 4
        len4 = record4.nbytes
        f.write(struct.pack('>i', len4))
        f.write(record4.tobytes())
        f.write(struct.pack('>i', len4))
    
    print(f"Written {fname}  ({n} zones, {n_y} Y values)")


def main():
    import sys
    
    if len(sys.argv) != 4:
        print("Usage: convert_hyd_abn_to_mod.py <hyd_file> <abn_file> <mod_file>")
        sys.exit(1)
    
    hyd_file = sys.argv[1]
    abn_file = sys.argv[2]
    mod_file = sys.argv[3]
    
    print(f"Reading {hyd_file}...")
    t_start, n, bm1, rcen, rhocen, hyd_data = read_hyd(hyd_file)
    
    print(f"Reading {abn_file}...")
    abn_data = read_abn(abn_file)
    
    print(f"Writing {mod_file}...")
    write_mod(mod_file, t_start, n, hyd_data, abn_data)
    
    print("Conversion complete.")


if __name__ == '__main__':
    main()
