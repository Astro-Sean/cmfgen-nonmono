#!/usr/bin/env python3
"""Extract a complete multi-band lightcurve from STELLA using pysynphot.

Reads the binary *.flx file (frequency-dependent emergent flux) and computes
UBVRI (Vega) and SDSS ugriz (AB) magnitudes with pysynphot synthetic photometry.
Where the *.flx file has gaps, the script falls back to a blackbody approximation
from the *.swd snapshots.

Late-time epochs are optionally calibrated to STELLA's native *.tt file.

Output format (matching existing lightcurve.csv with ugriz appended):
    time,U,B,V,R,I,u,g,r,i,z,bol_mag,bol_ergs

Usage:
    python extract_lightcurve_pysynphot.py SWD_FILE FLX_FILE [OPTIONS]
      --tt-file PATH     STELLA *.tt for zero-point calibration
      --output CSV       Output CSV path (default: <model>.lightcurve.csv)
      --plot PNG         Optional diagnostic plot
"""

import argparse
import contextlib
import io
import os
import struct
import sys
import warnings
from pathlib import Path

import numpy as np

# ── Ensure pysynphot can find its CDBS ────────────────────────────────────
if 'PYSYN_CDBS' not in os.environ:
    _cdbs_candidates = [
        '/ptmp/sbrennan/grp/redcat/trds',
        os.path.expanduser('~/grp/redcat/trds'),
        os.path.expanduser('~/.local/share/pysynphot/trds'),
    ]
    for _p in _cdbs_candidates:
        if Path(_p).exists():
            os.environ['PYSYN_CDBS'] = _p
            break
    else:
        print("ERROR: PYSYN_CDBS not found. Install CDBS or set the env var.")
        sys.exit(1)

# Suppress warnings BEFORE importing pysynphot
warnings.filterwarnings('ignore')

import pysynphot as S
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ── STELLA constants ──────────────────────────────────────────────────────
UR        = 1.0e14
BOLTZK    = 1.3806504e-16
UTP       = 1.0e5
CS        = 2.99792458e10
CSIGM     = 5.6704004778e-5  # Stefan-Boltzmann constant (STELLA value)
HPLANC    = 1.054571628e-27
CLUMNU    = np.log10(8.0 * UR**2 * (BOLTZK * UTP)**3 / (CS * HPLANC)**2)
UFREQ     = 2.08366e15
UPC       = 3.0856776e18
NFREQ     = 100
MFREQ     = 130
D_CM      = 10.0 * UPC
LSUN      = 3.828e33
MBOL_SUN  = 4.74
SIGMA_SB  = 5.670374e-5

# ── Frequency grid ───────────────────────────────────────────────────────
freq = np.zeros(NFREQ + 1)
freq[0]     = CS * 1.0e8 / (50000.0 * UFREQ)
freq[NFREQ] = CS * 1.0e8 / (1.0 * UFREQ)
Basis       = (50000.0 / 1.0)**(1.0 / NFREQ)
for i in range(1, NFREQ):
    freq[i] = freq[0] * Basis**i

freqmn = np.sqrt(freq[:-1] * freq[1:])

freqob = np.zeros(MFREQ)
freqob[:NFREQ] = freqmn
DLOGNU1 = 1.0 / np.log(freqmn[0] / freqmn[1])
for L in range(NFREQ, MFREQ):
    freqob[L] = freqmn[NFREQ - 1] * np.exp(-(L - NFREQ + 1) / DLOGNU1)

nu_hz  = freqob * UFREQ
lam_A  = CS / nu_hz * 1e8

dfreq = np.zeros(NFREQ)
dfreq[:-1] = np.diff(freq[:NFREQ]) * UFREQ
dfreq[-1]  = (freq[NFREQ] - freq[NFREQ - 1]) * UFREQ

# DLOGNU per bin (log spacing between consecutive freqmn bins, negative)
DLOGNU = np.zeros(NFREQ)
for L in range(NFREQ - 1):
    DLOGNU[L] = 1.0 / np.log(freqmn[L] / freqmn[L + 1])

# ── STELLA native UBVRI filter curves (from obsubvri.f) ──────────────────
# These are the exact transmission curves STELLA uses internally.
# Wavelengths in Angstroms, transmission normalized to 1.0 at peak.
STELLA_FILTERS = {
    'U': {
        'wl': np.array([3000,3050,3100,3150,3200,3250,3300,3350,3400,3450,
                        3500,3550,3600,3650,3700,3750,3800,3850,3900,3950,
                        4000,4050,4100,4150,4200], dtype=float),
        'T': np.array([0.000,0.016,0.068,0.167,0.287,0.423,0.560,0.673,0.772,0.841,
                       0.905,0.943,0.981,0.993,1.000,0.989,0.916,0.804,0.625,0.423,
                       0.238,0.114,0.051,0.019,0.000]),
        'F0': -13.90,
    },
    'B': {
        'wl': np.array([3600,3700,3800,3900,4000,4100,4200,4300,4400,4500,
                        4600,4700,4800,4900,5000,5100,5200,5300,5400,5500,
                        5600], dtype=float),
        'T': np.array([0.000,0.030,0.134,0.567,0.920,0.978,1.000,0.978,0.935,0.853,
                       0.740,0.640,0.536,0.424,0.325,0.235,0.150,0.095,0.043,0.009,
                       0.000]),
        'F0': -13.00,
    },
    'V': {
        'wl': np.array([4700,4800,4900,5000,5100,5200,5300,5400,5500,5600,
                        5700,5800,5900,6000,6100,6200,6300,6400,6500,6600,
                        6700,6800,6900,7000], dtype=float),
        'T': np.array([0.000,0.030,0.163,0.458,0.780,0.967,1.000,0.973,0.898,0.792,
                       0.684,0.574,0.461,0.359,0.270,0.197,0.135,0.081,0.045,0.025,
                       0.017,0.013,0.009,0.000]),
        'F0': -13.72,
    },
    'R': {
        'wl': np.array([5500,5600,5700,5800,5900,6000,6100,6200,6300,6400,
                        6500,6600,6700,6800,6900,7000,7100,7200,7300,7400,
                        7500,8000,8500,9000], dtype=float),
        'T': np.array([0.00,0.23,0.74,0.91,0.98,1.00,0.98,0.96,0.93,0.90,
                       0.86,0.81,0.78,0.72,0.67,0.61,0.56,0.51,0.46,0.40,
                       0.35,0.14,0.03,0.00]),
        'F0': -13.66,
    },
    'I': {
        'wl': np.array([7000,7100,7200,7300,7400,7500,7600,7700,7800,7900,
                        8000,8100,8200,8300,8400,8500,8600,8700,8800,8900,
                        9000,9100,9200], dtype=float),
        'T': np.array([0.000,0.024,0.232,0.555,0.785,0.910,0.965,0.985,0.990,0.995,
                       1.000,1.000,0.990,0.980,0.950,0.910,0.860,0.750,0.560,0.330,
                       0.150,0.030,0.000]),
        'F0': -14.42,
    },
}

# STELLA SCALE: distance modulus for 10 Mpc in STELLA internal units
# SCALE = -2.5 * log10(CFLUX * (UR / (10*UPC))^2)
CFLUX = 60.0 * CSIGM * (UTP / np.pi)**4
STELLA_SCALE = -2.5 * np.log10(CFLUX * (UR / (10.0 * UPC))**2)

# Pre-compute frequency bin index for each filter wavelength
# (replicates the indfreq lookup in obsubvri.f)
for _band, _info in STELLA_FILTERS.items():
    _inds = []
    for _wl in _info['wl']:
        # Find largest L such that FREQMN(L) > CCL / wl
        # (CCL = speed of light in cm/s, wl in Angstroms)
        _nu_target = CS * 1e8 / _wl  # Hz
        _nu_mn = freqmn * UFREQ      # Hz
        _L = NFREQ - 1
        while _L > 0 and _nu_mn[_L] > _nu_target:
            _L -= 1
        _inds.append(_L)
    _info['indfreq'] = np.array(_inds)


def stella_ubv(Flsave_col, nfrus):
    """Compute UBVRI magnitudes replicating STELLA's obsubvri.f subroutine.

    Uses the same filter transmission curves, log-log interpolation,
    and zero points as STELLA's internal computation.

    Parameters
    ----------
    Flsave_col : array (MFREQ,)
        Flsave values for one epoch (column from .flx data).
    nfrus : int
        Number of active frequency bins for this epoch.

    Returns
    -------
    dict with 'U', 'B', 'V', 'R', 'I' magnitudes (or np.nan if failed).
    """
    # CCL in STELLA = CS * 1e8 / UFREQ (speed of light in Angstroms/s per UFREQ)
    # This converts wavelength in Angstroms to frequency in UFREQ units:
    #   nu[UFREQ] = CCL / wl[Å]
    CCL = CS * 1e8 / UFREQ

    mags = {}
    for band, info in STELLA_FILTERS.items():
        wl = info['wl']
        T = info['T']
        inds = info['indfreq']
        F0 = info['F0']

        flxavg = 0.0
        wlstep = (wl[1] - wl[0]) / CCL  # Å / (cm/s) — mixed units, matches Fortran
        sumfhx = 0.0

        for i in range(1, len(wl) - 1):
            Lfr = inds[i]
            FHobs = max(abs(Flsave_col[Lfr]), 1e-100)
            FHobs1 = max(abs(Flsave_col[Lfr + 1]), 1e-100)

            if Lfr < nfrus:
                # Log-log interpolation in frequency
                # nu_target = CCL / wl[i] (in mixed units: cm/s / Å)
                # nu_mn_next = FREQMN(Lfr+1) (in UFREQ units)
                # DLOGNU(Lfr) = 1 / log(FREQMN(Lfr) / FREQMN(Lfr+1))
                # The Fortran uses FREQMN directly (in UFREQ units), and
                # CCL/wl also needs to be in UFREQ units for the log ratio.
                # CCL/wl[Å] = CS / wl_angstrom
                # In UFREQ units: (CS / wl_angstrom) / UFREQ
                nu_target = CCL / wl[i] / UFREQ
                nu_mn_next = freqmn[Lfr + 1]
                FHX = np.exp(
                    np.log(FHobs1) +
                    np.log(FHobs / FHobs1) * DLOGNU[Lfr] *
                    np.log(nu_target / nu_mn_next)
                )
            else:
                FHX = abs(FHobs)

            # (CCL/wl)^5 in mixed units (cm/s / Å)^5 — matches Fortran exactly
            FHX = FHX * T[i] * (CCL / wl[i])**5

            stepn = (wl[i + 1] - wl[i]) / CCL
            if abs(wlstep - stepn) / wlstep < 1e-4 and i != len(wl) - 2:
                sumfhx += FHX
            else:
                flxavg += sumfhx * wlstep + FHX * (wlstep + stepn) / 2.0
                wlstep = stepn
                sumfhx = 0.0

        if flxavg > 0:
            mags[band] = -2.5 * np.log10(flxavg) + F0 + STELLA_SCALE
        else:
            mags[band] = np.nan

    return mags


def parse_stella_native_lightcurve(res_path):
    """Parse STELLA's native lightcurve from the .res/.out/.tt file.

    Extracts the per-epoch photosphere parameters and UBVRI magnitudes
    that STELLA computes internally (written by PRIBAL every NOUT/10 steps).

    Returns dict with arrays: time, teff, rph, vph, mbol, U, B, V, R, I
    """
    import re
    data = []
    with open(res_path) as f:
        for line in f:
            m = re.match(
                r'^\s*(-?\d+\.\d+)\s+E\s+'
                r'(-?\d+\.\d+)\s+'
                r'(-?\d+\.\d+)\s+'
                r'(-?\d+\.\d+)\s+'
                r'(-?\d+\.\d+)\s+'
                r'(\d+\.\d+E[+-]\d+)\s+'
                r'(-?\d+\.\d+)\s+'
                r'(-?\d+\.\d+)\s+'
                r'(-?\d+\.\d+)\s+'
                r'(-?\d+\.\d+)\s+'
                r'(-?\d+\.\d+)\s+'
                r'(-?\d+\.\d+)',
                line.strip()
            )
            if m:
                data.append([float(m.group(k)) for k in range(1, 13)])

    if not data:
        return None

    arr = np.array(data)
    # Deduplicate (each line appears twice in .res)
    _, unique_idx = np.unique(arr[:, 0], return_index=True)
    arr = arr[np.sort(unique_idx)]

    return {
        'time':  arr[:, 0],
        'teff':  arr[:, 1],
        'rph':   arr[:, 2],
        'vph':   arr[:, 3],
        'xph':   arr[:, 4],
        'hused': arr[:, 5],
        'mbol':  arr[:, 6],
        'U':     arr[:, 7],
        'B':     arr[:, 8],
        'V':     arr[:, 9],
        'ub':    arr[:, 10],
        'bv':    arr[:, 11],
    }


# ── Filter definitions ───────────────────────────────────────────────────
FILTER_DEFS = {
    'U': {'pysyn_name': 'johnson,u', 'system': 'vegamag'},
    'B': {'pysyn_name': 'johnson,b', 'system': 'vegamag'},
    'V': {'pysyn_name': 'johnson,v', 'system': 'vegamag'},
    'R': {'pysyn_name': 'cousins,r', 'system': 'vegamag'},
    'I': {'pysyn_name': 'cousins,i', 'system': 'vegamag'},
    'u': {'pysyn_name': 'sdss,u',   'system': 'abmag'},
    'g': {'pysyn_name': 'sdss,g',   'system': 'abmag'},
    'r': {'pysyn_name': 'sdss,r',   'system': 'abmag'},
    'i': {'pysyn_name': 'sdss,i',   'system': 'abmag'},
    'z': {'pysyn_name': 'sdss,z',   'system': 'abmag'},
}

# ── Parsers ────────────────────────────────────────────────────────────────
def read_flx(flx_path: str):
    with open(flx_path, 'rb') as f:
        Nrec  = struct.unpack('i', f.read(4))[0]
        Lsave = struct.unpack('i', f.read(4))[0]
        times = np.zeros(Lsave)
        nfrus = np.zeros(Lsave, dtype=int)
        flsave = np.zeros((MFREQ, Lsave))
        for i in range(Lsave):
            times[i] = struct.unpack('f', f.read(4))[0]
            nfrus[i] = struct.unpack('i', f.read(4))[0]
            flsave[:, i] = np.array(struct.unpack(f'{MFREQ}d', f.read(8 * MFREQ)))
        end = struct.unpack('i', f.read(4))[0]
        assert end == Nrec
    return times, nfrus, flsave


def read_swd(swd_path: str, ncols: int = 14):
    """Read STELLA .swd file.

    Parses per-zone data from PRIBAL (FORMAT 93 in lbalsw.f):
      col 0:  time (days, observer frame)
      col 1:  Km (zone number, 1=NZON=outermost)
      col 2:  AMPR (Lagrangian mass, solar units)
      col 3:  Ry (radius, in UR=1e14 cm units)
      col 4:  Uy (velocity, in 1e6 cm/s / (UTIME*CRAP) units)
      col 5:  Ty (gas temperature, in UTP=1e5 K units)
      col 6:  TpRAD (radiation temperature, in UTP units)
      col 7:  PLLOG (log10 pressure)
      col 8:  PLOG (log10 density)
      col 9:  QVLOG
      col 10: QRTLG
      col 11: WRK (work)
      col 12: ENG (energy)
      col 13: Flum (luminosity, in 1e50 erg/s * R^2 units — see lbalsw.f:261)
      col 14+: WRKX, barionn, Xe, ... (not needed for photometry)

    Note: tau is NOT in the per-zone data.  It is in a separate %B: section.
    """
    epochs = {}
    current_t = None
    buf = []
    with open(swd_path) as fh:
        for line in fh:
            # STELLA .swd files can have numbers concatenated without spaces
            # (e.g., "10.91649213663.1725" should be "10.916492 13663.1725").
            # Use regex to split on transitions between digits and signs.
            import re
            cols = re.findall(r'[-+]?\d*\.?\d+[eE][-+]?\d+|[-+]?\d*\.?\d+', line.strip())
            if len(cols) < ncols:
                continue
            try:
                t = float(cols[0])
                km = int(float(cols[1]))
            except (ValueError, IndexError):
                continue
            if km == 1:
                if current_t is not None and buf:
                    arr = np.array(buf)
                    order = np.argsort(arr[:, 3])[::-1]
                    epochs[current_t] = arr[order]
                current_t = t
                buf = []
            if current_t is None:
                continue
            try:
                row = [float(c) for c in cols[:ncols]]
            except ValueError:
                continue
            buf.append(row)
    if current_t is not None and buf:
        arr = np.array(buf)
        order = np.argsort(arr[:, 3])[::-1]
        epochs[current_t] = arr[order]
    return epochs


# ── Spectrum helpers ─────────────────────────────────────────────────────
def flx_to_spectrum(Flsave, nfrus):
    """Convert STELLA Flsave to a pysynphot ArraySpectrum."""
    L_nu = Flsave * freqob**3 * 10.0**CLUMNU
    L_nu[nfrus:] = 0.0

    F_nu = L_nu / (4.0 * np.pi * D_CM**2)
    c_A = CS * 1e8
    F_lam = F_nu * c_A / lam_A**2

    mask = (lam_A > 0) & np.isfinite(F_lam) & (F_lam > 0)
    if not np.any(mask):
        return None, 0.0

    w = lam_A[mask]
    f = F_lam[mask]
    order = np.argsort(w)
    sp = S.ArraySpectrum(wave=w[order], flux=f[order], fluxunits='flam')

    L_bol = np.sum(L_nu[:NFREQ] * dfreq)
    return sp, L_bol


def swd_to_spectrum(snap):
    """Convert .swd snapshot to a blackbody pysynphot spectrum.

    Uses the zone with the highest radiation temperature as the photosphere
    proxy.  The .swd per-zone data does NOT contain tau (tau is in a separate
    %B: section), so we cannot interpolate to tau=2/3.  The T_rad maximum is
    a good proxy for the photosphere location in expanding ejecta.

    If Flum (col 13) is available, uses STELLA's own luminosity estimate
    instead of sigma*T^4 (more accurate when Flum is present).
    """
    # snap columns (see read_swd docstring):
    #   0=t, 1=km, 2=AMPR, 3=Ry, 4=Uy, 5=Ty(UTP), 6=TpRAD(UTP),
    #   7=PLLOG, 8=PLOG, 9=QVLOG, 10=QRTLG, 11=WRK, 12=ENG, 13=Flum
    ncols = snap.shape[1]

    # Temperatures are in UTP=1e5 K units (log10 is NOT applied — these are
    # linear values in UTP, e.g. Ty=9.528 means T=9.528e5 K).
    T_rad = snap[:, 6] * UTP  # radiation temperature in K
    T_gas = snap[:, 5] * UTP  # gas temperature in K

    # Use T_rad maximum as photosphere proxy (most reliable in optically thick
    # ejecta; T_rad peaks near the photosphere where radiation and matter are
    # most strongly coupled).
    i_photo = int(np.argmax(T_rad))
    T_photo = T_rad[i_photo]
    if T_photo < 10.0:
        T_photo = T_gas[i_photo]
    if T_photo <= 0 or not np.isfinite(T_photo):
        return None, 0.0

    # Radius in cm (col 3 is in UR=1e14 cm units)
    r_photo = snap[i_photo, 3] * UR

    # Bolometric luminosity: prefer STELLA's Flum if available
    if ncols > 13:
        # Flum (col 13) is in units of 1e50 erg/s * R^2 (see lbalsw.f:261:
        # Flum = (CLUMF * 1e-50) * Flum_internal * Ry(NZON-1)^2)
        # But the .swd writes Flum per-zone, and the total luminosity is
        # the value at the surface zone (km=1, outermost).
        # Actually, Flum in the .swd is the luminosity at that zone.
        # The surface luminosity is at the outermost zone (km=1 in STELLA
        # numbering, which is snap[0] after our outer->inner sort).
        Flum_surface = snap[0, 13]  # outermost zone's Flum
        # Flum is in 1e50 erg/s units (after the CLUMF*1e-50*R^2 normalization
        # in lbalsw.f:261, the .swd value is L / 1e50)
        L_bol = Flum_surface * 1e50
    else:
        # Fallback: Stefan-Boltzmann
        L_bol = 4.0 * np.pi * r_photo**2 * SIGMA_SB * T_photo**4

    if L_bol <= 0 or not np.isfinite(L_bol):
        return None, 0.0

    # Blackbody spectrum at the photosphere
    lam = np.logspace(2, 5, 1000)
    lam_cm = lam * 1e-8
    h = 6.62607015e-27
    k = BOLTZK
    B = 2.0 * h * CS**2 / lam_cm**5 / (np.expm1(h * CS / (lam_cm * k * T_photo)))
    F = np.pi * B * (r_photo / D_CM)**2 * 1e-8

    sp = S.ArraySpectrum(wave=lam, flux=F, fluxunits='flam')
    return sp, L_bol


# ── Quiet pysynphot wrapper ────────────────────────────────────────────────
def _quiet_observation(sp, bp):
    """Create Observation while suppressing pysynphot's stdout chatter."""
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        return S.Observation(sp, bp)


# ── Main extraction ────────────────────────────────────────────────────────
def extract(swd_path: str, flx_path: str, tt_path: str | None = None,
            max_epochs: int = 500, dt_min_days: float = 0.05):
    """Extract multi-band lightcurve from STELLA .swd + .flx files.

    Parameters
    ----------
    max_epochs : int
        Maximum number of epochs to process.  If the .flx file has more
        epochs, they are downsampled uniformly.  Default 500 (~12 epochs/day
        over 40 days, enough for smooth lightcurve plots without overwhelming
        pysynphot).
    dt_min_days : float
        Minimum time spacing between consecutive epochs in the output.
        Epochs closer than this are merged (the first one is kept).
        Default 0.05 days = 1.2 hours.
    """
    swd = read_swd(swd_path)
    flx_times, flx_nfrus, flx_data = read_flx(flx_path)
    print(f"[pysynphot] .swd epochs: {len(swd)}, .flx epochs: {len(flx_times)}")

    # ── Downsample .flx epochs if too many ──────────────────────────────
    # STELLA writes .flx every timestep (IOUT=1), producing thousands of
    # epochs per day.  For plotting we only need ~10/day.
    if len(flx_times) > max_epochs:
        step = max(1, len(flx_times) // max_epochs)
        flx_times = flx_times[::step]
        flx_nfrus = flx_nfrus[::step]
        flx_data = flx_data[:, ::step]
        print(f"[pysynphot] Downsampled .flx to {len(flx_times)} epochs "
              f"(every {step}th)")

    # ── Enforce minimum time spacing ────────────────────────────────────
    # Remove epochs closer than dt_min_days to the previous kept epoch.
    kept = [0]
    for i in range(1, len(flx_times)):
        if flx_times[i] - flx_times[kept[-1]] >= dt_min_days:
            kept.append(i)
    if len(kept) < len(flx_times):
        flx_times = flx_times[kept]
        flx_nfrus = flx_nfrus[kept]
        flx_data = flx_data[:, kept]
        print(f"[pysynphot] After dt_min={dt_min_days}d spacing: "
              f"{len(flx_times)} epochs")

    flx_idx = {round(t, 6): i for i, t in enumerate(flx_times)}

    # Unified time grid: .flx times + .swd times without nearby .flx match
    flx_time_arr = np.array(flx_times)
    all_times = list(flx_times)
    for t_swd in sorted(swd.keys()):
        if len(flx_time_arr) == 0 or np.min(np.abs(flx_time_arr - t_swd)) > 0.01:
            all_times.append(t_swd)

    all_times = sorted(set([round(t, 6) for t in all_times]))
    n = len(all_times)
    print(f"[pysynphot] Unified time grid: {n} epochs")

    # Pre-load bandpass objects quietly
    bandpasses = {}
    for band, info in FILTER_DEFS.items():
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            bandpasses[band] = S.ObsBandpass(info['pysyn_name'])

    results = {k: np.full(n, np.nan) for k in
               ['time', 'U', 'B', 'V', 'R', 'I', 'u', 'g', 'r', 'i', 'z',
                'bol_mag', 'bol_ergs']}

    flx_times_arr = np.array(flx_times)
    for i, t in enumerate(all_times):
        results['time'][i] = t
        sp = None
        L_bol = np.nan

        flx_key = round(t, 6)
        if flx_key in flx_idx:
            j = flx_idx[flx_key]
            sp, L_bol = flx_to_spectrum(flx_data[:, j], flx_nfrus[j])

        # If no exact .flx match, try linear interpolation between neighbouring .flx epochs
        if sp is None and len(flx_times_arr) > 1:
            idx = np.searchsorted(flx_times_arr, t)
            if 0 < idx < len(flx_times_arr):
                t0, t1 = flx_times_arr[idx - 1], flx_times_arr[idx]
                # Only interpolate if gap is small (< 3 days) to avoid bad extrapolation
                if t1 - t0 < 3.0:
                    w = (t - t0) / (t1 - t0)
                    sp0, L0 = flx_to_spectrum(flx_data[:, idx - 1], flx_nfrus[idx - 1])
                    sp1, L1 = flx_to_spectrum(flx_data[:, idx], flx_nfrus[idx])
                    if sp0 is not None and sp1 is not None:
                        # Interpolate flux arrays directly on a common wavelength grid
                        lam_common = np.logspace(2, 5, 1000)
                        f0 = np.interp(lam_common, sp0.wave, sp0.flux, left=0, right=0)
                        f1 = np.interp(lam_common, sp1.wave, sp1.flux, left=0, right=0)
                        f_interp = f0 * (1 - w) + f1 * w
                        sp = S.ArraySpectrum(wave=lam_common, flux=f_interp, fluxunits='flam')
                        L_bol = L0 * (1 - w) + L1 * w

        if sp is None and t in swd:
            sp, L_bol = swd_to_spectrum(swd[t])

        if sp is None:
            continue

        results['bol_ergs'][i] = L_bol
        results['bol_mag'][i] = MBOL_SUN - 2.5 * np.log10(L_bol / LSUN)

        for band, bp in bandpasses.items():
            try:
                obs = _quiet_observation(sp, bp)
                mag = obs.effstim(FILTER_DEFS[band]['system'])
                if np.isfinite(mag) and mag < 50:
                    results[band][i] = mag
            except Exception:
                pass

    # ── Calibrate against .tt ──────────────────────────────────────────
    if tt_path and Path(tt_path).exists():
        tt_data = {}
        with open(tt_path) as fh:
            for line in fh:
                s = line.strip()
                if not s or s.startswith('time') or s.startswith('XNianm'):
                    continue
                try:
                    parts = s.split()
                    if len(parts) >= 14:
                        t = float(parts[0])
                        if t > 0:
                            tt_data[t] = {
                                'U': float(parts[7]), 'B': float(parts[8]),
                                'V': float(parts[9]), 'R': float(parts[11]),
                                'I': float(parts[10]),
                                'bol': float(parts[6]),
                            }
                except ValueError:
                    continue

        if tt_data:
            print(f"[pysynphot] Calibrating against {len(tt_data)} .tt epochs")
            for band in ['U', 'B', 'V', 'R', 'I', 'bol_mag']:
                offsets = []
                for t_tt, vals in tt_data.items():
                    idx = np.argmin(np.abs(results['time'] - t_tt))
                    if band == 'bol_mag':
                        our_val = results['bol_mag'][idx]
                        ref = vals['bol']
                    else:
                        our_val = results[band][idx]
                        ref = vals[band]
                    if np.isfinite(our_val) and np.isfinite(ref):
                        offsets.append(ref - our_val)
                if offsets:
                    zp = np.median(offsets)
                    results[band] = results[band] + zp
                    print(f"  {band} offset = {zp:.3f} mag")

    return results


# ── Output ─────────────────────────────────────────────────────────────────
def write_csv(results: dict, outpath: str):
    cols = ['time', 'U', 'B', 'V', 'R', 'I', 'u', 'g', 'r', 'i', 'z', 'bol_mag', 'bol_ergs']
    with open(outpath, 'w') as f:
        f.write(','.join(cols) + '\n')
        for i in range(len(results['time'])):
            row = [f"{results['time'][i]:.6f}"]
            for c in cols[1:-1]:
                v = results[c][i]
                row.append(f"{v:.4f}" if np.isfinite(v) else "")
            row.append(f"{results['bol_ergs'][i]:.6e}" if np.isfinite(results['bol_ergs'][i]) else "")
            f.write(','.join(row) + '\n')
    print(f"CSV saved: {outpath}")


# ── Plotting ───────────────────────────────────────────────────────────────
def plot(results: dict, outpath: str, title: str, sn2011dh: bool = True):
    # Load SN 2011dh reference if requested
    sn2011dh_data = None
    sn2011dh_bol = None
    if sn2011dh:
        try:
            from load_sn2011dh import load_sn2011dh_ubvri, load_sn2011dh_bolometric
            sn2011dh_data = load_sn2011dh_ubvri()
            sn2011dh_bol = load_sn2011dh_bolometric()
        except Exception as e:
            print(f"WARNING: Could not load SN 2011dh reference: {e}")

    fig, axes = plt.subplots(2, 1, figsize=(11, 10), gridspec_kw={'height_ratios': [3, 1]})
    ax_top = axes[0]
    ax_bot = axes[1]

    t = results['time']
    pre = t < 0
    post = t >= 0

    colors = {
        'U': '#7B00D4', 'B': '#0000FF', 'V': '#00AA00',
        'R': '#FF4400', 'I': '#AA0000',
        'u': '#A020F0', 'g': '#228B22', 'r': '#FF4500',
        'i': '#8B0000', 'z': '#000000',
        'bol_mag': 'grey',
    }

    for band in ['U', 'B', 'V', 'R', 'I', 'u', 'g', 'r', 'i', 'z', 'bol_mag']:
        mag = results[band]
        valid = np.isfinite(mag)
        if not np.any(valid):
            continue
        c = colors.get(band, 'black')
        lbl = band
        if np.any(valid & post):
            ax_top.plot(t[valid & post], mag[valid & post],
                       color=c, ls='-', lw=1.2, marker='o', markersize=1.5,
                       alpha=0.7, label=lbl)
        if np.any(valid & pre):
            ax_top.plot(t[valid & pre], mag[valid & pre],
                       color=c, ls='--', lw=0.8, marker='s', markersize=1.5,
                       alpha=0.4)

    # Overlay SN 2011dh reference as scatter points
    if sn2011dh_data is not None:
        ref_bands = {
            'U': ('MU', '#7B00D4'), 'B': ('MB', '#0000FF'),
            'V': ('MV', '#00AA00'), 'R': ('MR', '#FF4400'),
            'I': ('MI', '#AA0000'),
        }
        for band, (key, color) in ref_bands.items():
            if band not in results:
                continue
            ref_mag = sn2011dh_data[key]
            ref_t = sn2011dh_data['phase']
            mask = ~np.isnan(ref_mag)
            if np.any(mask):
                ax_top.scatter(ref_t[mask], ref_mag[mask], color=color, marker='o',
                              s=15, alpha=0.4, edgecolors='none', zorder=5)

    ax_top.invert_yaxis()
    ax_top.set_xlim(0, max(t[-1] * 1.02, 40))
    ax_top.set_ylabel('Absolute Magnitude', fontsize=12)
    ax_top.set_xlabel('Time (days)', fontsize=12)
    ax_top.legend(fontsize=8, ncol=4, loc='lower right')
    ax_top.grid(True, alpha=0.3)
    plot_title = title + " vs SN 2011dh" if sn2011dh_data is not None else title
    ax_top.set_title(plot_title, fontsize=13)

    lum = results['bol_ergs']
    valid = np.isfinite(lum)
    if np.any(valid & post):
        ax_bot.semilogy(t[valid & post], lum[valid & post],
                       color='black', lw=1.2, marker='o', markersize=1.5, alpha=0.7,
                       label='STELLA')
    if np.any(valid & pre):
        ax_bot.semilogy(t[valid & pre], lum[valid & pre],
                       color='black', ls='--', lw=0.8, marker='s', markersize=1.5, alpha=0.4)

    # Overlay SN 2011dh bolometric
    if sn2011dh_bol is not None:
        ref_t = sn2011dh_bol['phase']
        ref_l = sn2011dh_bol['Lbol']
        mask = ref_t > 0
        ax_bot.semilogy(ref_t[mask], ref_l[mask], 'rx', markersize=4, alpha=0.5,
                       label='SN 2011dh', zorder=5)
        ax_bot.legend(fontsize=9, loc='upper right')

    ax_bot.set_xlim(0, max(t[-1] * 1.02, 40))
    ax_bot.set_ylabel(r'$L_\mathrm{bol}$ [erg/s]', fontsize=12)
    ax_bot.set_xlabel('Time (days)', fontsize=12)
    ax_bot.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    print(f"Plot saved: {outpath}")
    plt.close(fig)


# ── Native lightcurve output ──────────────────────────────────────────────
def write_native_csv(native: dict, outpath: str):
    """Write STELLA native lightcurve to CSV."""
    cols = ['time', 'teff', 'rph', 'vph', 'mbol', 'U', 'B', 'V']
    # Add R, I only if present (STELLA .res only writes UBV by default)
    if 'R' in native:
        cols.append('R')
    if 'I' in native:
        cols.append('I')
    with open(outpath, 'w') as f:
        f.write(','.join(cols) + '\n')
        for i in range(len(native['time'])):
            row = [f"{native['time'][i]:.6f}"]
            for c in cols[1:]:
                v = native[c][i] if c in native else np.nan
                row.append(f"{v:.4f}" if np.isfinite(v) else "")
            f.write(','.join(row) + '\n')
    print(f"Native CSV saved: {outpath}")


def plot_native(native: dict, outpath: str, title: str, sn2011dh: bool = True):
    """Plot STELLA native lightcurve from .res file."""
    sn2011dh_data = None
    sn2011dh_bol = None
    if sn2011dh:
        try:
            from load_sn2011dh import load_sn2011dh_ubvri, load_sn2011dh_bolometric
            sn2011dh_data = load_sn2011dh_ubvri()
            sn2011dh_bol = load_sn2011dh_bolometric()
        except Exception as e:
            print(f"WARNING: Could not load SN 2011dh reference: {e}")

    fig, axes = plt.subplots(2, 1, figsize=(11, 10), gridspec_kw={'height_ratios': [3, 1]})
    ax_top = axes[0]
    ax_bot = axes[1]

    t = native['time']
    pre = t < 0
    post = t >= 0

    colors = {'U': '#7B00D4', 'B': '#0000FF', 'V': '#00AA00',
              'R': '#FF4400', 'I': '#AA0000', 'mbol': 'grey'}

    for band in ['U', 'B', 'V', 'R', 'I', 'mbol']:
        if band not in native:
            continue
        mag = native[band]
        valid = np.isfinite(mag)
        if not np.any(valid):
            continue
        c = colors.get(band, 'black')
        lbl = 'M_bol' if band == 'mbol' else band
        ls = '--' if band == 'mbol' else '-'
        lw = 1.5 if band == 'mbol' else 1.2
        if np.any(valid & post):
            ax_top.plot(t[valid & post], mag[valid & post],
                       color=c, ls=ls, lw=lw, marker='o', markersize=2,
                       alpha=0.7, label=lbl)
        if np.any(valid & pre):
            ax_top.plot(t[valid & pre], mag[valid & pre],
                       color=c, ls='--', lw=0.8, marker='s', markersize=2, alpha=0.4)

    if sn2011dh_data is not None:
        ref_bands = {
            'U': ('MU', '#7B00D4'), 'B': ('MB', '#0000FF'),
            'V': ('MV', '#00AA00'), 'R': ('MR', '#FF4400'),
            'I': ('MI', '#AA0000'),
        }
        for band, (key, color) in ref_bands.items():
            ref_mag = sn2011dh_data[key]
            ref_t = sn2011dh_data['phase']
            mask = ~np.isnan(ref_mag)
            if np.any(mask):
                ax_top.scatter(ref_t[mask], ref_mag[mask], color=color, marker='o',
                              s=15, alpha=0.4, edgecolors='none', zorder=5)

    ax_top.invert_yaxis()
    ax_top.set_xlim(0, max(t[-1] * 1.02, 40))
    ax_top.set_ylabel('Absolute Magnitude', fontsize=12)
    ax_top.set_xlabel('Time (days)', fontsize=12)
    ax_top.legend(fontsize=8, ncol=4, loc='lower right')
    ax_top.grid(True, alpha=0.3)
    plot_title = title + " vs SN 2011dh" if sn2011dh_data is not None else title
    ax_top.set_title(plot_title, fontsize=13)

    # Bolometric luminosity from MBOL
    mbol = native['mbol']
    L_bol = 10**((4.75 - mbol) / 2.5) * 3.86e33  # erg/s (STELLA uses 3.86e33)
    valid = np.isfinite(L_bol) & (L_bol > 0)
    if np.any(valid & post):
        ax_bot.semilogy(t[valid & post], L_bol[valid & post],
                       color='black', lw=1.2, marker='o', markersize=2, alpha=0.7,
                       label='STELLA (SB)')
    if np.any(valid & pre):
        ax_bot.semilogy(t[valid & pre], L_bol[valid & pre],
                       color='black', ls='--', lw=0.8, marker='s', markersize=2, alpha=0.4)

    if sn2011dh_bol is not None:
        ref_t = sn2011dh_bol['phase']
        ref_l = sn2011dh_bol['Lbol']
        mask = ref_t > 0
        ax_bot.semilogy(ref_t[mask], ref_l[mask], 'rx', markersize=4, alpha=0.5,
                       label='SN 2011dh', zorder=5)
        ax_bot.legend(fontsize=9, loc='upper right')

    ax_bot.set_xlim(0, max(t[-1] * 1.02, 40))
    ax_bot.set_ylabel(r'$L_\mathrm{bol}$ [erg/s]', fontsize=12)
    ax_bot.set_xlabel('Time (days)', fontsize=12)
    ax_bot.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    print(f"Native plot saved: {outpath}")
    plt.close(fig)


# ── CLI ────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('swd_file', help='STELLA *.swd ASCII file')
    parser.add_argument('flx_file', help='STELLA *.flx binary file')
    parser.add_argument('--res-file', default=None,
                        help='STELLA *.res file for native lightcurve extraction. '
                             'If provided, extracts STELLA-native UBVRI/MBOL directly.')
    parser.add_argument('--tt-file', default=None, help='STELLA *.tt for calibration')
    parser.add_argument('--output', default=None, help='Output CSV path')
    parser.add_argument('--plot', default=None, help='Optional PNG plot')
    parser.add_argument('--title', default='STELLA Lightcurve (pysynphot)', help='Plot title')
    parser.add_argument('--no-sn2011dh', dest='sn2011dh', action='store_false', default=True,
                        help='Do not overlay SN 2011dh photometry')
    parser.add_argument('--max-epochs', type=int, default=500,
                        help='Max epochs to process (downsamples .flx if exceeded). '
                             'Default 500. Use 0 for no downsampling.')
    parser.add_argument('--dt-min', type=float, default=0.05,
                        help='Minimum time spacing between output epochs (days). '
                             'Default 0.05 (1.2h). Use 0 for no spacing filter.')
    args = parser.parse_args()

    for p in [args.swd_file, args.flx_file]:
        if not Path(p).exists():
            print(f"ERROR: {p} not found")
            sys.exit(1)

    # ── Native STELLA lightcurve from .res file (primary) ──────────────
    if args.res_file and Path(args.res_file).exists():
        native = parse_stella_native_lightcurve(args.res_file)
        if native is not None:
            print(f"\n[native] STELLA native lightcurve: {len(native['time'])} epochs "
                  f"from {native['time'][0]:.3f} to {native['time'][-1]:.3f} days")
            out_native_csv = args.output
            if out_native_csv is None:
                out_native_csv = str(Path(args.res_file).with_suffix('.native_lightcurve.csv'))
            else:
                out_native_csv = out_native_csv.replace('.csv', '.native.csv')
            write_native_csv(native, out_native_csv)

            if args.plot:
                native_plot = args.plot.replace('.png', '.native.png')
                plot_native(native, native_plot, title=args.title + ' (native)',
                           sn2011dh=args.sn2011dh)
        else:
            print(f"[native] WARNING: No native lightcurve found in {args.res_file}")

    # ── pysynphot lightcurve from .flx file (secondary) ────────────────
    results = extract(args.swd_file, args.flx_file, tt_path=args.tt_file,
                      max_epochs=args.max_epochs if args.max_epochs > 0 else 999999,
                      dt_min_days=args.dt_min if args.dt_min > 0 else 0.0)
    n_valid = np.sum(np.isfinite(results['bol_mag']))
    print(f"Extracted {n_valid}/{len(results['time'])} valid epochs "
          f"from {results['time'][0]:.3f} to {results['time'][-1]:.3f} days")

    out_csv = args.output
    if out_csv is None:
        out_csv = str(Path(args.swd_file).with_suffix('.lightcurve.csv'))
    write_csv(results, out_csv)

    if args.plot:
        plot(results, args.plot, title=args.title, sn2011dh=args.sn2011dh)


if __name__ == '__main__':
    main()
