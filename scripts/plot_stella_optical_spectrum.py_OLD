#!/usr/bin/env python3
"""Plot the optical spectrum (L_lambda vs wavelength) from a STELLA *.flx file.

STELLA stores the frequency-dependent emergent specific luminosity L_nu for
every saved snapshot in the binary *.flx file.  This script reads it, converts
to L_lambda, and plots the optical (3000-10000 A) spectrum for a selection of
epochs.

Usage:
    python plot_stella_optical_spectrum.py FLX_FILE [OPTIONS]
      --output PNG          Output plot path
      --epochs T1 T2 ...    Times (days) to plot (default: auto selection)
      --wmin WMIN           Minimum wavelength [A] (default 3000)
      --wmax WMAX           Maximum wavelength [A] (default 10000)
      --log                 Plot log y-axis
"""

import argparse
import struct
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ── STELLA constants (must match extract_lightcurve_from_flx.py) ──────────
UR        = 1.0e14          # cm  (unit of radius)
BOLTZK    = 1.3806504e-16   # erg/K
UTP       = 1.0e5           # K
CS        = 2.99792458e10   # cm/s
HPLANC    = 1.054571628e-27 # erg·s  (h-bar in STELLA)
CLUMNU    = np.log10(8.0 * UR**2 * (BOLTZK * UTP)**3 / (CS * HPLANC)**2)
UFREQ     = 2.08366e15      # Hz
NFREQ     = 100
MFREQ     = 130
WLMAX     = 50000.0         # Angstroms
WLMIN     = 1.0

# Frequency grid (dimensionless, same as STELLA)
freq = np.zeros(NFREQ + 1)
freq[0]    = CS * 1.0e8 / (WLMAX * UFREQ)
freq[NFREQ] = CS * 1.0e8 / (WLMIN  * UFREQ)
Basis = (WLMAX / WLMIN) ** (1.0 / NFREQ)
for i in range(1, NFREQ):
    freq[i] = freq[0] * Basis**i

freqmn = np.sqrt(freq[:-1] * freq[1:])           # bin centres (dimensionless)

freqob = np.zeros(MFREQ)
freqob[:NFREQ] = freqmn
DLOGNU1 = 1.0 / np.log(freqmn[0] / freqmn[1])
for L in range(NFREQ, MFREQ):
    freqob[L] = freqmn[NFREQ - 1] * np.exp(-(L - NFREQ + 1) / DLOGNU1)

nu_hz = freqob * UFREQ                            # physical frequency [Hz]
wl_A  = 1.0e8 * CS / nu_hz                        # wavelength [Angstrom]


def read_flx(flx_path: str):
    """Parse STELLA binary *.flx file.

    Returns:
      times  : np.ndarray (days), shape (Lsave,)
      nfrus  : np.ndarray (int),  shape (Lsave,)
      flsave : np.ndarray shape (MFREQ, Lsave)  -- raw STELLA flux variable F
    """
    with open(flx_path, 'rb') as f:
        Nrec  = struct.unpack('i', f.read(4))[0]
        Lsave = struct.unpack('i', f.read(4))[0]

        times  = np.zeros(Lsave)
        nfrus  = np.zeros(Lsave, dtype=int)
        flsave = np.zeros((MFREQ, Lsave))

        for i in range(Lsave):
            times[i]        = struct.unpack('f', f.read(4))[0]
            nfrus[i]        = struct.unpack('i', f.read(4))[0]
            flsave[:, i]    = np.array(struct.unpack(f'{MFREQ}d', f.read(8 * MFREQ)))

        end_marker = struct.unpack('i', f.read(4))[0]
        if end_marker != Nrec:
            print(f"WARNING: FLX record size mismatch ({end_marker} != {Nrec})")

    return times, nfrus, flsave


def flx_to_L_lambda(F, nfr):
    """Convert raw STELLA flux variable F to L_lambda [erg/s/A].

    L_nu  = F * freqob^3 * 10^CLUMNU   [erg/s/Hz]
    L_lam = L_nu * c / lambda^2        [erg/s/A]   (lambda in cm)
    """
    L_nu = F * freqob**3 * 10.0**CLUMNU          # erg/s/Hz
    lam_cm = wl_A * 1.0e-8                        # cm
    L_lam = L_nu * CS / lam_cm**2                 # erg/s/A
    # Mask out unused bins
    L_lam[nfr:] = np.nan
    return L_lam


def find_epoch_indices(times, requested_times):
    """Return list of (index, time) closest to each requested time."""
    out = []
    for t_req in requested_times:
        idx = int(np.argmin(np.abs(times - t_req)))
        out.append((idx, times[idx]))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('flx_file', help='STELLA *.flx binary file')
    ap.add_argument('--output', default=None,
                    help='Output PNG path (default: <flx>_optical_spectrum.png)')
    ap.add_argument('--epochs', nargs='*', type=float, default=None,
                    help='Times (days) to plot. Default: auto selection.')
    ap.add_argument('--wmin', type=float, default=3000.0, help='Min wavelength [A]')
    ap.add_argument('--wmax', type=float, default=10000.0, help='Max wavelength [A]')
    ap.add_argument('--log', action='store_true', help='Log y-axis')
    ap.add_argument('--title', default=None, help='Plot title')
    args = ap.parse_args()

    flx = Path(args.flx_file)
    if not flx.exists():
        print(f"ERROR: {flx} not found")
        sys.exit(1)

    times, nfrus, flsave = read_flx(str(flx))
    n_epochs = len(times)
    print(f"[flx] Read {n_epochs} epochs, t in [{times.min():.2f}, {times.max():.2f}] d")

    # Default epoch selection: a spread of post-SBO times
    if args.epochs is None:
        tmax = max(times.max(), 1.0)
        candidates = [1.0, 3.0, 5.0, 10.0, 20.0, 30.0, 40.0, 60.0]
        args.epochs = [t for t in candidates if t <= tmax]
        if not args.epochs:
            args.epochs = [times[times >= 0][0]] if np.any(times >= 0) else [times[-1]]

    epoch_list = find_epoch_indices(times, args.epochs)

    # Optical mask (wavelength decreasing with index, so just boolean mask)
    opt_mask = (wl_A >= args.wmin) & (wl_A <= args.wmax)

    fig, ax = plt.subplots(figsize=(11, 7))

    cmap = plt.cm.plasma
    n = len(epoch_list)
    for k, (idx, t_act) in enumerate(epoch_list):
        L_lam = flx_to_L_lambda(flsave[:, idx], nfrus[idx])
        w = wl_A[opt_mask]
        L = L_lam[opt_mask]
        # Sort by wavelength for plotting
        order = np.argsort(w)
        w = w[order]
        L = L[order]
        color = cmap(k / max(n - 1, 1))
        ax.plot(w, L, color=color, lw=1.6, alpha=0.9,
                label=f"t = {t_act:.2f} d")

    ax.set_xlim(args.wmin, args.wmax)
    if args.log:
        ax.set_yscale('log')
    ax.set_xlabel(r'Wavelength [$\mathrm{\AA}$]', fontsize=12)
    ax.set_ylabel(r'$L_\lambda$ [erg s$^{-1}$ $\mathrm{\AA}^{-1}$]', fontsize=12)
    title = args.title or f"STELLA optical spectrum — {flx.parent.parent.name}/{flx.stem}"
    ax.set_title(title, fontsize=12)
    ax.legend(fontsize=9, ncol=2, loc='upper right')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    out_png = args.output or str(flx.with_suffix('.optical_spectrum.png'))
    plt.savefig(out_png, dpi=150, bbox_inches='tight')
    print(f"Plot saved: {out_png}")
    plt.close(fig)


if __name__ == '__main__':
    main()
