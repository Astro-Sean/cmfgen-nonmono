#!/usr/bin/env python3
"""Load SN 2011dh multiband photometry from Ergon et al. (2014, 2015).

Data source: VizieR catalogs J/A+A/562/A17 (first 100 days) and
J/A+A/580/A142 (two years), downloaded to observational_data/SN2011dh/.

The photometry is in apparent magnitudes (Johnson-Cousins UBVRI, SDSS
ugriz, 2MASS JHK). To compare with STELLA model magnitudes (which are
absolute magnitudes at 10 pc), apply the distance modulus and extinction
correction:

    M_band = m_band - (mu - A_band)

For SN 2011dh in M51:
    Distance: 8.4 Mpc (Ergon et al. 2015)
    Distance modulus: mu = 29.61
    Extinction: A_V = 0.14 mag (E(B-V) = 0.045, R_V = 3.1)
        A_U = 0.21, A_B = 0.17, A_V = 0.14, A_R = 0.11, A_I = 0.08

Usage:
    from load_sn2011dh import load_sn2011dh_ubvri, load_sn2011dh_bolometric
    data = load_sn2011dh_ubvri()  # returns dict with 'phase', 'U', 'B', ...
    bol = load_sn2011dh_bolometric()  # returns dict with 'phase', 'Lbol'
"""

import os
import numpy as np

# ── SN 2011dh parameters ───────────────────────────────────────────────────
DISTANCE_MPC = 8.4          # M51 distance (Ergon et al. 2015)
DISTANCE_MODULUS = 5 * np.log10(DISTANCE_MPC * 1e6) - 5  # 29.61
EBV = 0.045                 # Galactic + host extinction
RV = 3.1
# Cardelli et al. (1989) extinction at Johnson-Cousins bands
A_V = EBV * RV              # 0.14
A_U = 1.530 * A_V           # ~0.21
A_B = 1.210 * A_V           # ~0.17
A_R = 0.780 * A_V           # ~0.11
A_I = 0.590 * A_V           # ~0.08

# Data directory (relative to project root)
_DATA_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "observational_data", "SN2011dh",
)


def _get_data_path(filename):
    """Return absolute path to a data file, or None if not found."""
    path = os.path.join(_DATA_DIR, filename)
    if os.path.isfile(path):
        return path
    return None


def load_sn2011dh_ubvri():
    """Load Johnson-Cousins UBVRI photometry for SN 2011dh.

    Returns
    -------
    dict with keys:
        'phase' : np.ndarray, days since explosion
        'U', 'B', 'V', 'R', 'I' : np.ndarray, apparent magnitudes (nan if missing)
        'eU', 'eB', 'eV', 'eR', 'eI' : np.ndarray, magnitude errors
        'MU', 'MB', 'MV', 'MR', 'MI' : np.ndarray, absolute magnitudes (extinction-corrected)
        'tel' : list of str, telescope/instrument
    """
    path = _get_data_path("tablec4.dat")
    if path is None:
        raise FileNotFoundError(
            f"SN 2011dh UBVRI data not found in {_DATA_DIR}. "
            "Download from VizieR J/A+A/580/A142/tablec4.dat"
        )

    phase, U, eU, B, eB, V, eV, R, eR, I, eI, tel = [], [], [], [], [], [], [], [], [], [], [], []

    with open(path) as f:
        for line in f:
            if len(line) < 17:
                continue
            try:
                p = float(line[10:16])
            except ValueError:
                continue

            def _f(s):
                s = s.strip()
                return float(s) if s else np.nan

            phase.append(p)
            U.append(_f(line[18:24]));  eU.append(_f(line[26:32]))
            B.append(_f(line[34:40]));  eB.append(_f(line[42:48]))
            V.append(_f(line[50:56]));  eV.append(_f(line[58:64]))
            R.append(_f(line[66:72]));  eR.append(_f(line[74:80]))
            I.append(_f(line[82:88]));  eI.append(_f(line[90:96]))
            tel.append(line[98:114].strip() if len(line) > 98 else "")

    phase = np.array(phase)
    U = np.array(U); eU = np.array(eU)
    B = np.array(B); eB = np.array(eB)
    V = np.array(V); eV = np.array(eV)
    R = np.array(R); eR = np.array(eR)
    I = np.array(I); eI = np.array(eI)

    # Convert to absolute magnitudes: M = m - mu + A_band
    # (subtract distance modulus, add extinction correction)
    MU = U - DISTANCE_MODULUS + A_U
    MB = B - DISTANCE_MODULUS + A_B
    MV = V - DISTANCE_MODULUS + A_V
    MR = R - DISTANCE_MODULUS + A_R
    MI = I - DISTANCE_MODULUS + A_I

    return {
        'phase': phase,
        'U': U, 'eU': eU, 'MU': MU,
        'B': B, 'eB': eB, 'MB': MB,
        'V': V, 'eV': eV, 'MV': MV,
        'R': R, 'eR': eR, 'MR': MR,
        'I': I, 'eI': eI, 'MI': MI,
        'tel': tel,
        'distance_modulus': DISTANCE_MODULUS,
        'A_V': A_V,
    }


def load_sn2011dh_bolometric():
    """Load bolometric luminosity for SN 2011dh.

    Returns
    -------
    dict with keys:
        'phase' : np.ndarray, days since explosion
        'Lbol' : np.ndarray, bolometric luminosity in erg/s
        'e_Lbol' : np.ndarray, statistical error in erg/s
        'Mbol' : np.ndarray, absolute bolometric magnitude
    """
    path = _get_data_path("tablec9.dat")
    if path is None:
        raise FileNotFoundError(
            f"SN 2011dh bolometric data not found in {_DATA_DIR}. "
            "Download from VizieR J/A+A/580/A142/tablec9.dat"
        )

    phase, lbol, elbol = [], [], []
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                p = float(parts[1])
                l = float(parts[2]) * 1e41  # convert from 1e41 erg/s to erg/s
                el = float(parts[3]) * 1e41 if len(parts) > 3 else 0.0
            except (ValueError, IndexError):
                continue
            phase.append(p)
            lbol.append(l)
            elbol.append(el)

    phase = np.array(phase)
    lbol = np.array(lbol)
    elbol = np.array(elbol)

    # Absolute bolometric magnitude: Mbol = -2.5*log10(L) + 4.74 - 5*log10(d/10pc)
    # But L is already the luminosity, so Mbol = -2.5*log10(L/Lsun) + Mbol_sun
    # where Mbol_sun = 4.74, Lsun = 3.828e33 erg/s
    LSUN = 3.828e33
    MBOL_SUN = 4.74
    mbol = -2.5 * np.log10(lbol / LSUN) + MBOL_SUN

    return {
        'phase': phase,
        'Lbol': lbol,
        'e_Lbol': elbol,
        'Mbol': mbol,
    }


def load_sn2011dh_ugriz():
    """Load SDSS ugriz photometry for SN 2011dh.

    Returns dict with 'phase', 'u', 'g', 'r', 'i', 'z' (apparent mags).
    """
    path = _get_data_path("tablec5.dat")
    if path is None:
        raise FileNotFoundError(
            f"SN 2011dh ugriz data not found in {_DATA_DIR}."
        )

    phase, u, g, r, i, z = [], [], [], [], [], []
    with open(path) as f:
        for line in f:
            if len(line) < 17:
                continue
            try:
                p = float(line[10:16])
            except ValueError:
                continue

            def _f(s):
                s = s.strip()
                return float(s) if s else np.nan

            phase.append(p)
            u.append(_f(line[18:24]))
            g.append(_f(line[34:40]))
            r.append(_f(line[50:56]))
            i.append(_f(line[66:72]))
            z.append(_f(line[82:88]))

    return {
        'phase': np.array(phase),
        'u': np.array(u), 'g': np.array(g),
        'r': np.array(r), 'i': np.array(i), 'z': np.array(z),
    }


if __name__ == "__main__":
    # Quick test
    print("=== SN 2011dh UBVRI (Ergon et al. 2014, 2015) ===")
    ubvri = load_sn2011dh_ubvri()
    print(f"  {len(ubvri['phase'])} epochs, phase {ubvri['phase'].min():.1f}-{ubvri['phase'].max():.1f} d")
    v_mask = ~np.isnan(ubvri['V'])
    v_peak = np.argmin(ubvri['V'][v_mask])
    print(f"  V peak: {ubvri['V'][v_mask][v_peak]:.2f} (apparent) / "
          f"{ubvri['MV'][v_mask][v_peak]:.2f} (absolute) at "
          f"t={ubvri['phase'][v_mask][v_peak]:.1f} d")

    print("\n=== SN 2011dh Bolometric ===")
    bol = load_sn2011dh_bolometric()
    peak = np.argmax(bol['Lbol'])
    print(f"  {len(bol['phase'])} epochs, phase {bol['phase'].min():.1f}-{bol['phase'].max():.1f} d")
    print(f"  Peak: L={bol['Lbol'][peak]:.2e} erg/s, Mbol={bol['Mbol'][peak]:.2f} "
          f"at t={bol['phase'][peak]:.1f} d")
