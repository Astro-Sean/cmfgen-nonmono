#!/usr/bin/env python3
"""Extract a complete multi-band lightcurve from STELLA .swd + .flx files.

STELLA's native *.tt output has a ~59-day gap for low-energy explosions
because tt4strad only writes when the emergent flux changes by >0.1%.  The
*.flx file has the same gap (182 epochs vs 216 in *.swd).

This script uses a hybrid approach:
  1.  *.flx  -> actual frequency-dependent SED for UBVRI where available
  2.  *.swd  -> photosphere T, R to fill gaps with blackbody approximation
  3.  Late-time magnitudes are calibrated to STELLA's *.tt file.

Usage:
    python extract_complete_lightcurve.py SWD_FILE FLX_FILE [OPTIONS]
      --tt-file PATH   STELLA *.tt for calibration
      --output CSV     Output CSV path
      --plot PNG       Output plot path
"""

import argparse
import struct
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar

# ── STELLA constants ──────────────────────────────────────────────────────
UR        = 1.0e14
BOLTZK    = 1.3806504e-16
UTP       = 1.0e5
CS        = 2.99792458e10
HPLANC    = 1.054571628e-27
CLUMNU    = np.log10(8.0 * UR**2 * (BOLTZK * UTP)**3 / (CS * HPLANC)**2)
UFREQ     = 2.08366e15
UPC       = 3.0856776e18
NFREQ     = 100
MFREQ     = 130
WLMAX     = 50000.0
WLMIN     = 1.0
LSUN      = 3.828e33
MBOL_SUN  = 4.74
SIGMA_SB  = 5.670374e-5

# ── Frequency grid ─────────────────────────────────────────────────────────
freq = np.zeros(NFREQ + 1)
freq[0]       = CS * 1.0e8 / (WLMAX * UFREQ)
freq[NFREQ]   = CS * 1.0e8 / (WLMIN * UFREQ)
Basis         = (WLMAX / WLMIN)**(1.0 / NFREQ)
for i in range(1, NFREQ):
    freq[i] = freq[0] * Basis**i

freqmn = np.sqrt(freq[:-1] * freq[1:])

freqob = np.zeros(MFREQ)
freqob[:NFREQ] = freqmn
DLOGNU1 = 1.0 / np.log(freqmn[0] / freqmn[1])
for L in range(NFREQ, MFREQ):
    freqob[L] = freqmn[NFREQ - 1] * np.exp(-(L - NFREQ + 1) / DLOGNU1)

dfreq = np.zeros(NFREQ)
dfreq[:-1] = np.diff(freq[:NFREQ]) * UFREQ
dfreq[-1]  = (freq[NFREQ] - freq[NFREQ - 1]) * UFREQ

nu_hz = freqob * UFREQ

# ── Filter bands ───────────────────────────────────────────────────────────
FILTER_BANDS = {
    'U': {'lam_eff': 3650.0, 'dlam':  680.0},
    'B': {'lam_eff': 4450.0, 'dlam':  970.0},
    'V': {'lam_eff': 5510.0, 'dlam':  890.0},
    'R': {'lam_eff': 6580.0, 'dlam': 1380.0},
    'I': {'lam_eff': 8060.0, 'dlam': 1490.0},
}

# ── Blackbody helpers ──────────────────────────────────────────────────────
def planck_nu(T, nu):
    h = 6.62607015e-27
    k = BOLTZK
    x = h * nu / (k * T)
    with np.errstate(over='ignore'):
        return (2.0 * h * nu**3 / CS**2) / (np.expm1(x))


def planck_lambda(T, lam_cm):
    h = 6.62607015e-27
    k = BOLTZK
    x = h * CS / (lam_cm * k * T)
    with np.errstate(over='ignore'):
        return (2.0 * h * CS**2 / lam_cm**5) / (np.expm1(x))


def bb_filter_fluxes(T, R_cm, bands):
    D = 10.0 * UPC
    area = (R_cm / D)**2
    out = {}
    for band, props in bands.items():
        lam_eff = props['lam_eff'] * 1.0e-8
        dlam    = props['dlam'] * 1.0e-8
        lams = np.linspace(lam_eff - 0.5*dlam, lam_eff + 0.5*dlam, 50)
        B = planck_lambda(T, lams)
        F_mean = np.pi * np.trapz(B, lams) / dlam
        out[band] = F_mean * area
    return out


def fit_bb_t_r(nu_hz, L_nu, nfrus):
    nu = nu_hz[:nfrus]
    L  = L_nu[:nfrus]
    mask = (L > 0.0) & np.isfinite(L)
    nu = nu[mask]
    L  = L[mask]
    if len(nu) < 5:
        return 5000.0, 1e15
    w = nu**4

    def chi2(logT):
        T = 10.0**logT[0]
        B = planck_nu(T, nu)
        R2 = np.sum(w * L * B) / np.sum(w * B**2)
        if R2 <= 0:
            return 1e30
        return np.sum(w * (np.log10(L) - np.log10(R2 * B))**2)

    logT_grid = np.linspace(3.0, 5.0, 200)
    chi2_vals = np.array([chi2([lt]) for lt in logT_grid])
    best = logT_grid[np.argmin(chi2_vals)]
    res = minimize_scalar(lambda lt: chi2([lt]), bounds=(best-0.1, best+0.1), method='bounded')
    T_best = 10.0**res.x
    B = planck_nu(T_best, nu)
    R2 = np.sum(w * L * B) / np.sum(w * B**2)
    R_best = np.sqrt(max(R2, 0.0)) / (2.0 * np.pi)
    return T_best, R_best


# ── Parsers ────────────────────────────────────────────────────────────────
def read_swd(swd_path: str):
    """Read ASCII .swd file.  Returns dict keyed by t_day -> snapshot arrays."""
    epochs = {}
    current_t = None
    buf = []
    with open(swd_path) as fh:
        for line in fh:
            cols = line.strip().split()
            if len(cols) < 13:
                continue
            try:
                t = float(cols[0])
                km = int(cols[1])
            except ValueError:
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
            row = [float(c) for c in cols[:13]]
            buf.append(row)
    if current_t is not None and buf:
        arr = np.array(buf)
        order = np.argsort(arr[:, 3])[::-1]
        epochs[current_t] = arr[order]
    return epochs


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


# ── Main extraction ────────────────────────────────────────────────────────
def extract(swd_path: str, flx_path: str, tt_path: str | None = None):
    swd = read_swd(swd_path)
    flx_times, flx_nfrus, flx_data = read_flx(flx_path)
    print(f"[complete] .swd epochs: {len(swd)}, .flx epochs: {len(flx_times)}")

    # Build lookup for .flx
    flx_idx = {round(t, 6): i for i, t in enumerate(flx_times)}

    # Unified time grid: all unique times from .swd and .flx, sorted
    all_times = sorted(set(list(swd.keys()) + [round(t, 6) for t in flx_times]))
    n = len(all_times)
    print(f"[complete] Unified time grid: {n} epochs")

    results = {k: np.zeros(n) for k in ['time','Tbb','rbb','Teff','Mbol','L_bol',
                                           'U','B','V','R','I']}

    for i, t in enumerate(all_times):
        results['time'][i] = t
        use_flx = False
        flx_key = round(t, 6)
        if flx_key in flx_idx:
            j = flx_idx[flx_key]
            F = flx_data[:, j]
            L_nu = F * freqob**3 * 10.0**CLUMNU
            L_bol = np.sum(L_nu[:NFREQ] * dfreq)
            T, R = fit_bb_t_r(nu_hz, L_nu, flx_nfrus[j])
            raw = bb_filter_fluxes(T, R, FILTER_BANDS)
            for band in FILTER_BANDS:
                results[band][i] = raw[band]
            results['Tbb'][i] = T
            results['rbb'][i] = R
            if R > 0:
                Teff = (L_bol / (4.0 * np.pi * R**2 * SIGMA_SB))**0.25
            else:
                Teff = T
            results['Teff'][i] = Teff
            results['L_bol'][i] = L_bol
            results['Mbol'][i] = MBOL_SUN - 2.5 * np.log10(L_bol / LSUN)
            use_flx = True

        if not use_flx and t in swd:
            # Fallback to .swd blackbody
            snap = swd[t]
            r_photo = (10.0**snap[:, 3])[0]
            T_photo = (10.0**snap[:, 6])[0] if (10.0**snap[:, 6])[0] > 0 else (10.0**snap[:, 5])[0]
            L_bol = 4.0 * np.pi * r_photo**2 * SIGMA_SB * T_photo**4
            results['L_bol'][i] = L_bol
            results['Mbol'][i] = MBOL_SUN - 2.5 * np.log10(L_bol / LSUN)
            results['Tbb'][i] = T_photo
            results['rbb'][i] = r_photo
            results['Teff'][i] = T_photo
            raw = bb_filter_fluxes(T_photo, r_photo, FILTER_BANDS)
            for band in FILTER_BANDS:
                results[band][i] = raw[band]

    # ── Calibrate zero points against .tt ──────────────────────────────────
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
                                'Mbol': float(parts[6]),
                                'U': float(parts[7]), 'B': float(parts[8]),
                                'V': float(parts[9]), 'R': float(parts[11]),
                                'I': float(parts[10]),
                            }
                except ValueError:
                    continue

        if tt_data:
            print(f"[complete] Calibrating against {len(tt_data)} .tt epochs")
            for band in list(FILTER_BANDS.keys()) + ['Mbol']:
                offsets = []
                for t_tt, vals in tt_data.items():
                    idx = np.argmin(np.abs(results['time'] - t_tt))
                    if band == 'Mbol':
                        our_val = results['Mbol'][idx]
                    else:
                        our_val = -2.5 * np.log10(results[band][idx])
                    offsets.append(vals[band] - our_val)
                if offsets:
                    zp = np.median(offsets)
                    if band == 'Mbol':
                        results['Mbol'] += zp
                    else:
                        results[band] = -2.5 * np.log10(results[band]) + zp
                    print(f"  {band} offset = {zp:.3f} mag")
    else:
        for band in FILTER_BANDS:
            results[band] = -2.5 * np.log10(results[band] + 1e-300) - 25.0

    return results


# ── Plotting ───────────────────────────────────────────────────────────────
def plot(lc: dict, outpath: str, title: str):
    fig, axes = plt.subplots(2, 1, figsize=(10, 10), gridspec_kw={'height_ratios': [3, 1]})
    ax_top = axes[0]
    ax_bot = axes[1]

    t = lc['time']
    pre = t < 0
    post = t >= 0

    colors = {'U': '#7B00D4', 'B': '#0000FF', 'V': '#00AA00',
              'R': '#FF4400', 'I': '#AA0000', 'Mbol': 'grey'}

    for band in ['Mbol', 'U', 'B', 'V', 'R', 'I']:
        mag = lc[band]
        c = colors.get(band, 'black')
        lbl = band
        if np.any(post):
            ax_top.plot(t[post], mag[post], color=c, ls='-', lw=1.3,
                       marker='o', markersize=1.5, alpha=0.7, label=lbl)
        if np.any(pre):
            ax_top.plot(t[pre], mag[pre], color=c, ls='--', lw=1.0,
                       marker='s', markersize=1.5, alpha=0.5)

    ax_top.invert_yaxis()
    ax_top.set_xlim(t[0], t[-1] * 1.02)
    ax_top.set_ylabel('Magnitude', fontsize=12)
    ax_top.set_xlabel('Time (days)', fontsize=12)
    ax_top.legend(fontsize=9, ncol=3, loc='lower right')
    ax_top.grid(True, alpha=0.3)
    ax_top.set_title(title, fontsize=13)

    lum = lc['L_bol']
    if np.any(post):
        ax_bot.semilogy(t[post], lum[post], color='black', lw=1.3,
                       marker='o', markersize=1.5, alpha=0.7)
    if np.any(pre):
        ax_bot.semilogy(t[pre], lum[pre], color='black', ls='--', lw=1.0,
                       marker='s', markersize=1.5, alpha=0.5)

    ax_bot.set_xlim(t[0], t[-1] * 1.02)
    ax_bot.set_ylabel(r'$L_\mathrm{bol}$ [erg/s]', fontsize=12)
    ax_bot.set_xlabel('Time (days)', fontsize=12)
    ax_bot.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    print(f"Plot saved: {outpath}")
    plt.close(fig)


# ── CSV ────────────────────────────────────────────────────────────────────
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


# ── CLI ────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('swd_file', help='STELLA *.swd ASCII file')
    parser.add_argument('flx_file', help='STELLA *.flx binary file')
    parser.add_argument('--tt-file', default=None, help='STELLA *.tt for calibration')
    parser.add_argument('--output', default=None, help='Output CSV')
    parser.add_argument('--plot', default=None, help='Output PNG')
    parser.add_argument('--title', default='STELLA Complete Lightcurve', help='Plot title')
    args = parser.parse_args()

    for p in [args.swd_file, args.flx_file]:
        if not Path(p).exists():
            print(f"ERROR: {p} not found")
            sys.exit(1)

    lc = extract(args.swd_file, args.flx_file, tt_path=args.tt_file)
    print(f"Extracted {len(lc['time'])} epochs from {lc['time'].min():.3f} to {lc['time'].max():.3f} days")

    out_csv = args.output or str(Path(args.swd_file).with_suffix('.complete_lightcurve.csv'))
    write_csv(lc, out_csv)

    out_png = args.plot or str(Path(args.swd_file).with_suffix('.complete_lightcurve.png'))
    plot(lc, out_png, title=args.title)


if __name__ == '__main__':
    main()
