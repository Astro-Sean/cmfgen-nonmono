#!/usr/bin/env python3
"""Extract a lightcurve from STELLA .lbol file (fallback when .flx is unavailable).

STELLA's .lbol file contains:
    time       L_ubvri    L_bol    XEUV<325   IR>890
All luminosities are in log10(erg/s).

This script converts to absolute magnitudes and writes a CSV compatible
with extract_lightcurve_pysynphot.py output format.

Usage:
    python extract_lightcurve_from_lbol.py LBOL_FILE [--output CSV] [--plot PNG]
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

LSUN = 3.828e33
MBOL_SUN = 4.74


def load_lbol(path):
    """Load STELLA .lbol file. Returns dict with time, L_ubvri, L_bol, L_xeuv, L_ir (all in log10 erg/s)."""
    data = np.loadtxt(path, skiprows=1)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return {
        'time': data[:, 0],
        'L_ubvri': data[:, 1],
        'L_bol': data[:, 2],
        'L_xeuv': data[:, 3] if data.shape[1] > 3 else np.full(len(data), np.nan),
        'L_ir': data[:, 4] if data.shape[1] > 4 else np.full(len(data), np.nan),
    }


def lbol_to_mbol(lbol_log):
    """Convert log10(L_erg/s) to absolute bolometric magnitude."""
    lbol = 10.0 ** lbol_log
    return -2.5 * np.log10(lbol / LSUN) + MBOL_SUN


def main():
    parser = argparse.ArgumentParser(description="Extract lightcurve from STELLA .lbol file")
    parser.add_argument("lbol_file", help="Path to STELLA .lbol file")
    parser.add_argument("--output", "-o", default="lightcurve_from_lbol.csv", help="Output CSV path")
    parser.add_argument("--plot", "-p", default=None, help="Output plot PNG path")
    parser.add_argument("--days-max", type=float, default=None, help="Max days to plot")
    args = parser.parse_args()

    if not Path(args.lbol_file).exists():
        print(f"ERROR: .lbol file not found: {args.lbol_file}", file=sys.stderr)
        sys.exit(1)

    data = load_lbol(args.lbol_file)
    t = data['time']
    mbol = lbol_to_mbol(data['L_bol'])
    lbol_ergs = 10.0 ** data['L_bol']

    # Write CSV (compatible with extract_lightcurve_pysynphot.py format)
    # We only have bolometric from .lbol; UBVRI would need .tt
    with open(args.output, 'w') as f:
        f.write("time,bol_mag,bol_ergs\n")
        for i in range(len(t)):
            f.write(f"{t[i]:.6f},{mbol[i]:.4f},{lbol_ergs[i]:.6e}\n")
    print(f"Written {len(t)} epochs to {args.output}")
    print(f"  Time range: {t.min():.3f} to {t.max():.3f} days")
    print(f"  Mbol range: {mbol.min():.2f} to {mbol.max():.2f}")
    peak = np.argmin(mbol)
    print(f"  Peak: t={t[peak]:.1f}d, Mbol={mbol[peak]:.2f}, L={lbol_ergs[peak]:.2e} erg/s")

    # Plot
    if args.plot:
        # Try to load SN 2011dh reference
        sn2011dh_bol = None
        try:
            from load_sn2011dh import load_sn2011dh_bolometric
            sn2011dh_bol = load_sn2011dh_bolometric()
        except Exception as e:
            print(f"WARNING: Could not load SN 2011dh reference: {e}")

        fig, ax = plt.subplots(1, 1, figsize=(10, 7))

        # Plot bolometric
        mask = t >= 0 if np.any(t >= 0) else np.ones_like(t, dtype=bool)
        ax.plot(t[mask], mbol[mask], 'k-', label='STELLA Mbol', linewidth=2)

        # Overlay SN 2011dh bolometric
        if sn2011dh_bol is not None:
            ref_t = sn2011dh_bol['phase']
            ref_mbol = sn2011dh_bol['Mbol']
            ref_mask = ref_t > 0
            ax.scatter(ref_t[ref_mask], ref_mbol[ref_mask], c='red', marker='x',
                      s=30, alpha=0.6, label='SN 2011dh', zorder=5)

        if args.days_max:
            ax.set_xlim(0, args.days_max)
        else:
            ax.set_xlim(0, max(t.max(), 40))
        ax.set_ylim(ax.get_ylim()[::-1])  # Magnitudes inverted
        ax.set_xlabel('Time since explosion (days)', fontsize=13)
        ax.set_ylabel('Absolute bolometric magnitude', fontsize=13)
        title = 'SN_12C STELLA Lightcurve (from .lbol)'
        if sn2011dh_bol is not None:
            title += ' vs SN 2011dh'
        ax.set_title(title, fontsize=14)
        ax.legend(loc='best', fontsize=11)
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(args.plot, dpi=150)
        print(f"Plot saved: {args.plot}")


if __name__ == "__main__":
    main()
