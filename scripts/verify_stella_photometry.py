#!/usr/bin/env python3
"""Verify the STELLA SED -> band-magnitude photometry chain.

Checks:
  1. Reproduce the hardcoded zeropoints (SGF0, SRF0, ACF0, AOF0, ...) in
     tt4strad.f from the tabulated energy-response curves via
     FX0 = 2.5*log10(c_AA * int(T/lam^2)dlam) - 48.6075.
  2. Recompute absolute AB magnitudes from the binary .flx file with an
     independent integration  m = -2.5*log10(<F_nu>_E / F_nu,AB)  at 10 pc
     and compare against the .tt columns (MGG, MGR, MAC, MAO).
  3. Cross-check binary .flx L_nu against the formatted .ph records.
"""
import re
import struct
import sys

import numpy as np

# STELLA internal constants (from begradsep.f / tt4strad.f)
UR = 1.0e14          # cm, radius unit
BOLTZK = 1.3806504e-16
UTP = 1.0e5
CS = 2.99792458e10   # cm/s
HPLANC = 1.054571628e-27  # h-bar (STELLA convention)
C_AA = CS * 1.0e8    # c in Angstrom/s
UFREQ = BOLTZK * UTP / (2.0 * np.pi * HPLANC)  # ~2.08366e15 Hz
CLUMNU = np.log10(8.0 * UR**2 * (BOLTZK * UTP)**3 / (CS * HPLANC)**2)
WLMIN, WLMAX = 1.0, 50000.0
D10PC = 3.0856776e19  # cm
FNU_AB = 3631.0e-23   # erg/s/cm^2/Hz

TT4 = sys.argv[1] if len(sys.argv) > 1 else (
    "/ptmp/sbrennan/stella/OpenStella_parallel/strad/tt4strad.f")


def parse_data_array(src_text, name):
    """Extract a Fortran DATA array `dataName/.../ ` incl. continuations."""
    lines = src_text.splitlines()
    vals = []
    collecting = False
    for ln in lines:
        if not collecting:
            m = re.match(r"\s*data" + name + r"\s*/(.*)", ln, re.IGNORECASE)
            if m:
                collecting = True
                vals.append(m.group(1))
        else:
            # continuation lines have non-space char in col 6 or start with *
            stripped = ln.strip()
            if stripped.startswith("*") or (len(ln) > 5 and ln[5] not in " "):
                vals.append(stripped.lstrip("*"))
            else:
                break
    txt = "".join(vals)
    txt = txt.split("/")[0] if "/" in txt else txt
    out = []
    for tok in txt.split(","):
        tok = tok.strip()
        if tok:
            out.append(float(tok.replace("d", "e").replace("D", "e")))
    return np.array(out)


def main():
    src = open(TT4).read()

    bands = {}
    for key, wl_name, t_name in [
        ("u", "WLSu", "TSu"), ("g", "WLSg", "TSg"), ("r", "WLSr", "TSr"),
        ("i", "WLSi", "TSi"), ("z", "WLSz", "TSz"),
        ("J", "WLJ", "TJ"), ("H", "WLH", "TH"), ("K", "WLK", "TK"),
        ("c", "WLAc", "TAc"), ("o", "WLAo", "TAo"),
    ]:
        bands[key] = (parse_data_array(src, wl_name),
                      parse_data_array(src, t_name))

    # ---- check 1: zeropoint reproduction -----------------------------------
    print("=== Zeropoint reproduction (AB: dAB=0) ===")
    zp_ref = {"u": -13.3376, "g": -13.2090, "r": -13.8102, "i": -14.1818,
              "z": -14.5089, "J": -15.8265, "H": -16.4271, "K": -17.4441,
              "c": -12.7035, "o": -13.1260}
    dab = {"J": 0.9, "H": 1.38, "K": 1.85}
    for k, (wl, t) in bands.items():
        integ = np.trapezoid(t / wl**2, wl)
        zp = 2.5 * np.log10(C_AA * integ) - 48.6075 - dab.get(k, 0.0)
        print(f"  {k:>2}: computed {zp:9.4f}  vs tt4strad {zp_ref[k]:9.4f}"
              f"   diff {zp - zp_ref[k]:+.4f}")

    return bands


if __name__ == "__main__":
    main()


def detect_mfreq(flx_path):
    """Infer MFREQ from the .flx record layout: 8-byte header, Lsave
    records of f4+i4+MFREQ*f8, 4-byte trailer."""
    import os
    size = os.path.getsize(flx_path)
    with open(flx_path, "rb") as f:
        nrec, lsave = struct.unpack("ii", f.read(8))
    body = size - 8 - 4
    if body % lsave != 0:
        raise ValueError(
            f"{flx_path}: body {body} B not divisible by Lsave={lsave}")
    mfreq = (body // lsave - 8) // 8
    if mfreq <= 0 or (8 + 8 * mfreq) * lsave != body:
        raise ValueError(f"{flx_path}: cannot infer MFREQ (got {mfreq})")
    return mfreq


def build_freqob(nfreq, mfreq):
    """Reconstruct STELLA observer-frequency grid (dimensionless units)."""
    freq = np.zeros(nfreq + 1)
    freq[0] = CS * 1.0e8 / (WLMAX * UFREQ)
    freq[nfreq] = CS * 1.0e8 / (WLMIN * UFREQ)
    basis = (WLMAX / WLMIN) ** (1.0 / nfreq)
    for i in range(1, nfreq):
        freq[i] = freq[0] * basis**i
    freqmn = np.sqrt(freq[:-1] * freq[1:])
    freqob = np.zeros(mfreq)
    freqob[:nfreq] = freqmn
    dlognu1 = 1.0 / np.log(freqmn[0] / freqmn[1])
    for l in range(nfreq, mfreq):
        freqob[l] = freqmn[nfreq - 1] * np.exp(-(l - nfreq + 1) / dlognu1)
    return freqob


def read_flx(path, mfreq):
    """Read binary .flx: Nrec, Lsave, then per record time(f4), nfrus(i4),
    MFREQ doubles of raw F."""
    with open(path, "rb") as fh:
        nrec, lsave = struct.unpack("ii", fh.read(8))
        times = np.zeros(lsave)
        nfrus = np.zeros(lsave, int)
        f = np.zeros((mfreq, lsave))
        for i in range(lsave):
            times[i] = struct.unpack("f", fh.read(4))[0]
            nfrus[i] = struct.unpack("i", fh.read(4))[0]
            f[:, i] = np.array(struct.unpack(f"{mfreq}d", fh.read(8 * mfreq)))
    return times, nfrus, f


def band_mag_AB(nu_hz, lnu, nfrus, wl_A, T):
    """Energy-response-weighted AB absolute magnitude at 10 pc."""
    lam_band = C_AA / nu_hz[:nfrus]          # wavelength per grid bin
    flx_nu = lnu[:nfrus] / (4 * np.pi * D10PC**2)
    # interpolate tabulated T(lambda) onto grid wavelengths
    T_g = np.interp(lam_band, wl_A, T, left=0.0, right=0.0)
    num = np.trapezoid(flx_nu * T_g, nu_hz[:nfrus])
    den = np.trapezoid(T_g, nu_hz[:nfrus])
    if num <= 0 or den <= 0:
        return np.nan
    return -2.5 * np.log10(num / den) - 48.60


def check_flx_vs_tt(flx_path, tt_path, nfreq, mfreq, bands):
    freqob = build_freqob(nfreq, mfreq)
    nu_hz = freqob * UFREQ
    times, nfrus, flx = read_flx(flx_path, mfreq)
    print(f"  .flx: {len(times)} records, nfrus range {nfrus.min()}-{nfrus.max()}")
    # .tt columns: 13=g, 14=r, 20=c, 21=o
    tt = np.loadtxt(tt_path)
    tt_t, tt_g, tt_r, tt_c, tt_o = tt[:, 0], tt[:, 13], tt[:, 14], tt[:, 20], tt[:, 21]
    rows = []
    for i in range(len(times)):
        lnu = flx[:, i] * freqob**3 * 10**CLUMNU
        rec = {"t": times[i]}
        for key, col in (("g", tt_g), ("r", tt_r), ("c", tt_c), ("o", tt_o)):
            wl, t = bands[key]
            rec[key] = band_mag_AB(nu_hz, lnu, nfrus[i], wl, t)
            rec[key + "_tt"] = np.interp(times[i], tt_t, col)
        rows.append(rec)
    return rows
