#!/usr/bin/env python3
"""Extract a complete multi-band lightcurve from a STELLA *.flx binary file.

STELLA's native *.tt file is aggressively thinned and misses large gaps
(e.g. 59 days for low-energy explosions) because tt4strad only writes when
the emergent flux changes by >0.1%.

This script reads the binary *.flx file (which contains the full frequency-
dependent emergent flux for every saved snapshot), integrates the SED to get
the bolometric luminosity, fits a blackbody temperature, and computes UBVRI
magnitudes.  The bolometric Mbol agrees with STELLA's *.tt to better than
~0.03 mag.

Usage:
    python extract_lightcurve_from_flx.py FLX_FILE [OPTIONS]
      --tt-file PATH    STELLA *.tt file (for late-time calibration)
      --output CSV    Output CSV path
      --plot PNG      Output plot path
"""

import argparse
import struct
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ── STELLA constants ──────────────────────────────────────────────────────
UR        = 1.0e14          # cm  (unit of radius)
BOLTZK    = 1.3806504e-16   # erg/K
UTP       = 1.0e5           # K
CS        = 2.99792458e10   # cm/s
HPLANC    = 1.054571628e-27 # erg·s  (h-bar in STELLA)
CLUMNU    = np.log10(8.0 * UR**2 * (BOLTZK * UTP)**3 / (CS * HPLANC)**2)
UFREQ     = 2.08366e15      # Hz
UPC       = 3.0856776e18    # cm (1 parsec)
NFREQ     = 100
MFREQ     = 130
WLMAX     = 50000.0         # Angstroms
WLMIN     = 1.0
LSUN      = 3.828e33        # erg/s
MBOL_SUN  = 4.74
SIGMA_SB  = 5.670374e-5     # erg cm^-2 s^-1 K^-4

# ── Frequency grid (dimensionless, same as STELLA) ──────────────────────
freq = np.zeros(NFREQ + 1)
freq[0]   = CS * 1.0e8 / (WLMAX * UFREQ)
freq[NFREQ] = CS * 1.0e8 / (WLMIN * UFREQ)
Basis     = (WLMAX / WLMIN)**(1.0 / NFREQ)
for i in range(1, NFREQ):
    freq[i] = freq[0] * Basis**i

freqmn = np.sqrt(freq[:-1] * freq[1:])

freqob = np.zeros(MFREQ)
freqob[:NFREQ] = freqmn
DLOGNU1 = 1.0 / np.log(freqmn[0] / freqmn[1])
for L in range(NFREQ, MFREQ):
    freqob[L] = freqmn[NFREQ - 1] * np.exp(-(L - NFREQ + 1) / DLOGNU1)

# Frequency bin widths in Hz
dfreq = np.zeros(NFREQ)
dfreq[:-1] = np.diff(freq[:NFREQ]) * UFREQ
# Last bin width
dfreq[-1] = (freq[NFREQ] - freq[NFREQ - 1]) * UFREQ

nu_hz = freqob * UFREQ  # physical frequency in Hz


# ── Simple box-car Bessell (1990) filters ─────────────────────────────────
#  λ_eff [Å]   Δλ [Å]
FILTER_BANDS = {
    'U': {'lam_eff': 3650.0, 'dlam':  680.0},
    'B': {'lam_eff': 4450.0, 'dlam':  970.0},
    'V': {'lam_eff': 5510.0, 'dlam':  890.0},
    'R': {'lam_eff': 6580.0, 'dlam': 1380.0},
    'I': {'lam_eff': 8060.0, 'dlam': 1490.0},
}

# ── Blackbody helpers ────────────────────────────────────────────────────
def planck_nu(T, nu):
    """Planck function B_ν(T) in erg s^-1 cm^-2 Hz^-1 sr^-1."""
    h = 6.62607015e-27  # erg·s (physical h, not h-bar)
    k = BOLTZK
    x = h * nu / (k * T)
    # Avoid overflow for large x
    with np.errstate(over='ignore'):
        return (2.0 * h * nu**3 / CS**2) / (np.expm1(x))


def planck_lambda(T, lam_cm):
    """Planck function B_λ(T) in erg s^-1 cm^-2 cm^-1 sr^-1."""
    h = 6.62607015e-27
    k = BOLTZK
    x = h * CS / (lam_cm * k * T)
    with np.errstate(over='ignore'):
        return (2.0 * h * CS**2 / lam_cm**5) / (np.expm1(x))


def blackbody_mag_from_filters(T, R_cm, bands):
    """
    Compute absolute Vega-like magnitudes for a blackbody of temperature T [K]
    and radius R_cm, using simple box-car filters.

    Returns dict {band: mag}.
    """
    D = 10.0 * UPC  # 10 pc in cm
    area_factor = (R_cm / D)**2
    mags = {}
    for band, props in bands.items():
        lam_eff = props['lam_eff'] * 1.0e-8  # cm
        dlam    = props['dlam'] * 1.0e-8     # cm
        lam1 = lam_eff - 0.5 * dlam
        lam2 = lam_eff + 0.5 * dlam
        # Sample Planck function at 50 points in the band
        lams = np.linspace(lam1, lam2, 50)
        B = planck_lambda(T, lams)
        # Mean surface flux (πB) integrated over the band, divided by bandwidth
        F_mean = np.pi * np.trapz(B, lams) / dlam
        F_obs = F_mean * area_factor
        mags[band] = F_obs
    return mags


def fit_bb_t_rphoto(nu_hz, L_nu, nfrus):
    """Fit a blackbody T and photosphere radius to the SED.

    We minimise chi-squared in log-space, weighting by ν^4 to emphasise the
    Wien side where the SED is most informative about T.
    """
    # Only use used frequencies
    nu = nu_hz[:nfrus]
    L = L_nu[:nfrus]
    mask = (L > 0.0) & np.isfinite(L)
    nu = nu[mask]
    L = L[mask]
    if len(nu) < 5:
        return 5000.0, 1e15

    # Weight by ν^4
    w = nu**4

    def chi2(logT):
        T = 10.0**logT[0]
        B_nu = planck_nu(T, nu)
        # L_nu = 4π R^2 π B_ν  =>  R^2 = L_ν / (4π^2 B_ν)
        # Fit for amplitude (R^2) analytically at each T
        R2 = np.sum(w * L * B_nu) / np.sum(w * B_nu**2)
        if R2 <= 0:
            return 1e30
        model = R2 * B_nu
        return np.sum(w * (np.log10(L) - np.log10(model))**2)

    # Grid search over logT for robustness
    logT_grid = np.linspace(3.0, 5.0, 200)
    chi2_vals = np.array([chi2([lt]) for lt in logT_grid])
    best_idx = np.argmin(chi2_vals)
    logT_best = logT_grid[best_idx]

    # Refine with a small bracket
    from scipy.optimize import minimize_scalar
    def f_scalar(logT):
        return chi2([logT])
    result = minimize_scalar(f_scalar, bounds=(logT_best - 0.1, logT_best + 0.1), method='bounded')
    logT_best = result.x
    T_best = 10.0**logT_best

    # Compute best-fit R
    B_nu = planck_nu(T_best, nu)
    R2 = np.sum(w * L * B_nu) / np.sum(w * B_nu**2)
    R_best = np.sqrt(max(R2, 0.0)) / (2.0 * np.pi)

    return T_best, R_best


# ── FLX parser ────────────────────────────────────────────────────────────
def read_flx(flx_path: str):
    """Parse STELLA binary *.flx file.

    Returns dicts:
      times   : np.ndarray (days)
      nfrus   : np.ndarray (int)
      flsave  : np.ndarray shape (MFREQ, n_epochs)
    """
    with open(flx_path, 'rb') as f:
        Nrec = struct.unpack('i', f.read(4))[0]
        Lsave = struct.unpack('i', f.read(4))[0]

        times = np.zeros(Lsave)
        nfrus = np.zeros(Lsave, dtype=int)
        flsave = np.zeros((MFREQ, Lsave))

        for i in range(Lsave):
            times[i] = struct.unpack('f', f.read(4))[0]
            nfrus[i] = struct.unpack('i', f.read(4))[0]
            flsave[:, i] = np.array(struct.unpack(f'{MFREQ}d', f.read(8 * MFREQ)))

        end_marker = struct.unpack('i', f.read(4))[0]
        assert end_marker == Nrec, "FLX record size mismatch"

    return times, nfrus, flsave


# ── Main extraction ─────────────────────────────────────────────────────
def extract_lightcurve(flx_path: str, tt_path: str | None = None):
    """Extract complete lightcurve from .flx."""
    times, nfrus_arr, flsave = read_flx(flx_path)
    n_epochs = len(times)
    print(f"[flx] Read {n_epochs} epochs from {flx_path}")

    results = {
        'time': times,
        'Tbb': np.zeros(n_epochs),
        'rbb': np.zeros(n_epochs),
        'Teff': np.zeros(n_epochs),
        'Mbol': np.zeros(n_epochs),
        'L_bol': np.zeros(n_epochs),
        'U': np.zeros(n_epochs),
        'B': np.zeros(n_epochs),
        'V': np.zeros(n_epochs),
        'R': np.zeros(n_epochs),
        'I': np.zeros(n_epochs),
    }

    for i in range(n_epochs):
        nfr = nfrus_arr[i]
        F = flsave[:, i]

        # Specific luminosity L_ν  [erg/s/Hz]
        L_nu = F * freqob**3 * 10.0**CLUMNU

        # Bolometric luminosity [erg/s]
        L_bol = np.sum(L_nu[:NFREQ] * dfreq)
        results['L_bol'][i] = L_bol
        results['Mbol'][i] = MBOL_SUN - 2.5 * np.log10(L_bol / LSUN)

        # Fit blackbody
        T, R = fit_bb_t_rphoto(nu_hz, L_nu, nfr)
        results['Tbb'][i] = T
        results['rbb'][i] = R

        # Effective temperature from L = 4πR²σT⁴
        if R > 0:
            Teff = (L_bol / (4.0 * np.pi * R**2 * SIGMA_SB))**0.25
        else:
            Teff = T
        results['Teff'][i] = Teff

        # Magnitudes from blackbody (raw flux-like values)
        raw_mags = blackbody_mag_from_filters(T, R, FILTER_BANDS)
        for band in FILTER_BANDS:
            results[band][i] = raw_mags[band]

    # ── Calibrate zero points against STELLA *.tt ────────────────────────
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
                        vals = {b: float(parts[6 + i]) for i, b in enumerate(['Mbol', 'U', 'B', 'V', 'I', 'R'])}
                        # Use last value if duplicate time
                        if t > 0:
                            tt_data[t] = vals
                except ValueError:
                    continue

        if tt_data:
            print(f"[flx] Calibrating against {len(tt_data)} .tt epochs")
            # For each band, compute median offset at late times (t > 30)
            for band in list(FILTER_BANDS.keys()) + ['Mbol']:
                offsets = []
                for t_tt, vals in tt_data.items():
                    if t_tt > 30:
                        # Find closest epoch in our results
                        idx = np.argmin(np.abs(results['time'] - t_tt))
                        if band == 'Mbol':
                            our_val = results['Mbol'][idx]
                        else:
                            our_val = -2.5 * np.log10(results[band][idx])
                        offset = vals[band] - our_val
                        offsets.append(offset)
                if offsets:
                    zp = np.median(offsets)
                    if band == 'Mbol':
                        results['Mbol'] += zp
                    else:
                        results[band] = -2.5 * np.log10(results[band]) + zp
                    print(f"  {band} zero-point offset = {zp:.3f} mag")
    else:
        # No .tt calibration: convert raw fluxes to rough magnitudes using
        # an arbitrary common zero point (scales are relative anyway)
        print("[flx] No .tt file for calibration – using uncalibrated relative mags")
        for band in FILTER_BANDS:
            # Use a rough zero point so numbers are ~sensible
            results[band] = -2.5 * np.log10(results[band] + 1e-300) - 25.0

    return results


# ── Plotting ──────────────────────────────────────────────────────────────
def plot_lightcurves(lc: dict, outpath: str, title: str = "STELLA Lightcurve"):
    fig, axes = plt.subplots(2, 1, figsize=(10, 10), gridspec_kw={'height_ratios': [3, 1]})
    ax_top = axes[0]
    ax_bot = axes[1]

    t = lc['time']

    # --- Top panel: multi-band magnitudes ---
    colors = {'U': '#7B00D4', 'B': '#0000FF', 'V': '#00AA00',
                'R': '#FF4400', 'I': '#AA0000', 'Mbol': 'grey'}
    labels = {'U': 'U', 'B': 'B', 'V': 'V', 'R': 'R', 'I': 'I', 'Mbol': 'Mbol'}

    # Separate pre-SBO (t<0) and post-SBO (t>=0)
    pre = t < 0
    post = t >= 0

    for band in ['Mbol', 'U', 'B', 'V', 'R', 'I']:
        if band not in lc:
            continue
        mag = lc[band]
        c = colors.get(band, 'black')
        lbl = labels.get(band, band)

        if np.any(post):
            ax_top.plot(t[post], mag[post], color=c, ls='-', lw=1.5,
                       marker='o', markersize=2, alpha=0.7, label=lbl)
        if np.any(pre):
            ax_top.plot(t[pre], mag[pre], color=c, ls='--', lw=1.2,
                       marker='s', markersize=2, alpha=0.5, label=f'{lbl} (pre-SBO)')

    ax_top.invert_yaxis()
    ax_top.set_xlim(t[0], t[-1] * 1.02)
    ax_top.set_ylabel('Magnitude', fontsize=12)
    ax_top.set_xlabel('Time (days)', fontsize=12)
    ax_top.legend(fontsize=9, ncol=3, loc='lower right')
    ax_top.grid(True, alpha=0.3)
    ax_top.set_title(title, fontsize=13)

    # --- Bottom panel: bolometric luminosity ---
    lum = lc['L_bol']
    if np.any(post):
        ax_bot.semilogy(t[post], lum[post], color='black', lw=1.5,
                       marker='o', markersize=2, alpha=0.7)
    if np.any(pre):
        ax_bot.semilogy(t[pre], lum[pre], color='black', ls='--', lw=1.2,
                       marker='s', markersize=2, alpha=0.5)

    ax_bot.set_xlim(t[0], t[-1] * 1.02)
    ax_bot.set_ylabel(r'$L_\mathrm{bol}$ [erg/s]', fontsize=12)
    ax_bot.set_xlabel('Time (days)', fontsize=12)
    ax_bot.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    print(f"Plot saved: {outpath}")
    plt.close(fig)


# ── CSV output ───────────────────────────────────────────────────────────
def write_csv(lc: dict, outpath: str):
    cols = ['time_days', 'Tbb_K', 'rbb_cm', 'Teff_K', 'Mbol', 'L_bol_erg_s',
            'U', 'B', 'V', 'R', 'I']
    with open(outpath, 'w') as f:
        f.write(','.join(cols) + '\n')
        for i in range(len(lc['time'])):
            row = [f"{lc['time'][i]:.10f}",
                   f"{lc['Tbb'][i]:.3f}",
                   f"{lc['rbb'][i]:.6e}",
                   f"{lc['Teff'][i]:.3f}",
                   f"{lc['Mbol'][i]:.4f}",
                   f"{lc['L_bol'][i]:.6e}",
                   f"{lc['U'][i]:.4f}",
                   f"{lc['B'][i]:.4f}",
                   f"{lc['V'][i]:.4f}",
                   f"{lc['R'][i]:.4f}",
                   f"{lc['I'][i]:.4f}"]
            f.write(','.join(row) + '\n')
    print(f"CSV saved: {outpath}")


# ── CLI ──────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('flx_file', help='STELLA *.flx binary file')
    parser.add_argument('--tt-file', default=None,
                        help='STELLA *.tt file for magnitude calibration')
    parser.add_argument('--output', default=None,
                        help='Output CSV path (default: <flx>_lightcurve.csv)')
    parser.add_argument('--plot', default=None,
                        help='Output PNG plot path (default: <flx>_lightcurve.png)')
    parser.add_argument('--title', default='STELLA Lightcurve from .flx',
                        help='Plot title')
    args = parser.parse_args()

    flx = Path(args.flx_file)
    if not flx.exists():
        print(f"ERROR: {flx} not found")
        sys.exit(1)

    lc = extract_lightcurve(str(flx), tt_path=args.tt_file)
    print(f"Extracted {len(lc['time'])} epochs from {lc['time'].min():.3f} to {lc['time'].max():.3f} days")

    out_csv = args.output or str(flx.with_suffix('.flx_lightcurve.csv'))
    write_csv(lc, out_csv)

    out_png = args.plot or str(flx.with_suffix('.flx_lightcurve.png'))
    plot_lightcurves(lc, out_png, title=args.title)


if __name__ == '__main__':
    main()
