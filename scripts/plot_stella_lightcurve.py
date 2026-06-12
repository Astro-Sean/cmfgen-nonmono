#!/usr/bin/env python3
"""Plot STELLA lightcurve from a *.tt file.

STELLA is a 1D radiation-hydrodynamics code (explicit multi-group radiation
transport coupled to Lagrangian hydro), so its *.tt output contains full
multi-band photometric lightcurves:

  time  Tbb  rbb  Teff  Rlast_sc  R(tau2/3)  Mbol  MU  MB  MV  MI  MR  Mbolavg  gdepos

Usage:
    python plot_stella_lightcurve.py TT_FILE [OPTIONS]
      --output PATH    Save figure to this path (default: alongside tt file)
      --distance D     Distance modulus in mag (default: 0 = absolute)
      --extinction E   V-band extinction A_V in mag (default: 0)
      --days-max N     X-axis upper limit in days (default: auto)
      --bands B        Comma-separated bands to plot (default: Mbol,V,B,U,R,I)
      --no-invert      Do not invert y-axis (default: inverted, brighter up)
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


BAND_COLS = {
    'Mbol':   (6,  'grey',    '--',  'Mbol'),
    'U':      (7,  '#7B00D4', '-',   'U'),
    'B':      (8,  '#0000FF', '-',   'B'),
    'V':      (9,  '#00AA00', '-',   'V'),
    'I':      (10, '#AA0000', '-',   'I'),
    'R':      (11, '#FF4400', '-',   'R'),
    'Mbolavg':(12, 'black',   ':',   r'$\langle$Mbol$\rangle$'),
}

# Bessell (1990) standard extinction coefficients R_lambda = A_lambda / E(B-V)
# using R_V = 3.1
EXTINCTION_COEFF = {
    'U': 4.78,
    'B': 4.10,
    'V': 3.10,
    'R': 2.32,
    'I': 1.49,
    'Mbol': 3.10,   # approximate: treat as V
    'Mbolavg': 3.10,
}


def read_tt(tt_path: str) -> dict:
    """Read a STELLA *.tt file and return a dict of arrays."""
    data = []
    header_found = False
    with open(tt_path) as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith('time'):
                header_found = True
                continue
            if stripped.startswith('XNianm'):
                continue
            try:
                row = [float(x) for x in stripped.split()]
                if len(row) >= 14:
                    data.append(row)
            except ValueError:
                continue
    if not data:
        raise ValueError(f"No data rows found in {tt_path}. "
                         f"{'Header line found.' if header_found else 'No header line found either.'}")
    arr = np.array(data)
    return {
        'time':    arr[:, 0],
        'Tbb':     arr[:, 1],
        'rbb':     arr[:, 2],
        'Teff':    arr[:, 3],
        'Rlast_sc':arr[:, 4],
        'Rtau23':  arr[:, 5],
        'Mbol':    arr[:, 6],
        'U':       arr[:, 7],
        'B':       arr[:, 8],
        'V':       arr[:, 9],
        'I':       arr[:, 10],
        'R':       arr[:, 11],
        'Mbolavg': arr[:, 12],
        'gdepos':  arr[:, 13],
    }


def write_lightcurve_csv(tt: dict, outpath: str, distance_modulus: float = 0.0):
    """Write lightcurve data to a CSV with columns: time,U,B,V,R,I,bol_mag,bol_ergs."""
    import csv

    LSUN = 3.828e33  # erg/s per solar luminosity
    MBOL_SUN = 4.74

    t = tt['time']
    mbol = tt['Mbol'] + distance_modulus
    lum_erg_s = LSUN * 10.0 ** ((MBOL_SUN - mbol) / 2.5)

    csv_path = Path(outpath)
    with open(csv_path, 'w', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow(['time', 'U', 'B', 'V', 'R', 'I', 'bol_mag', 'bol_ergs'])
        for i in range(len(t)):
            writer.writerow([
                f"{t[i]:.6f}",
                f"{tt['U'][i]:.4f}",
                f"{tt['B'][i]:.4f}",
                f"{tt['V'][i]:.4f}",
                f"{tt['R'][i]:.4f}",
                f"{tt['I'][i]:.4f}",
                f"{mbol[i]:.4f}",
                f"{lum_erg_s[i]:.6e}",
            ])
    print(f"Lightcurve CSV saved: {csv_path}")


def plot_lightcurves(tt: dict, bands: list[str], distance_modulus: float,
                     av_extinction: float, days_max: float | None,
                     invert_y: bool, model_name: str, outpath: str):
    """Make a multi-panel lightcurve figure."""
    fig, axes = plt.subplots(2, 1, figsize=(9, 10),
                             gridspec_kw={'height_ratios': [3, 1]})
    ax_lc = axes[0]
    ax_bol = axes[1]

    t = tt['time']
    rv = 3.1
    ebv = av_extinction / rv

    # Separate pre-SBO (negative) and post-SBO (positive) times
    pre_sbo_mask = t < 0
    post_sbo_mask = t >= 0

    # --- Top panel: individual band lightcurves ---
    for band in bands:
        if band not in BAND_COLS:
            print(f"WARNING: Unknown band '{band}', skipping.")
            continue
        col_idx, color, ls, label = BAND_COLS[band]
        mag = tt[band].copy()
        # Apply distance modulus + extinction
        ext = EXTINCTION_COEFF.get(band, 3.1) * ebv
        mag = mag + distance_modulus + ext
        
        # Plot post-SBO (positive times) with solid lines
        if np.any(post_sbo_mask):
            ax_lc.plot(t[post_sbo_mask], mag[post_sbo_mask], color=color, ls='-', lw=1.8, label=label)
        
        # Plot pre-SBO (negative times) with dashed lines
        if np.any(pre_sbo_mask):
            ax_lc.plot(t[pre_sbo_mask], mag[pre_sbo_mask], color=color, ls='--', lw=1.8, alpha=0.7, label=f'{label} (pre-SBO)')

    if invert_y:
        ax_lc.invert_yaxis()
    if days_max is not None:
        ax_lc.set_xlim(t[0], days_max)
    else:
        ax_lc.set_xlim(t[0], t[-1] * 1.05)
    ax_lc.set_ylabel('Magnitude', fontsize=12)
    ax_lc.set_xlabel('Time (days)', fontsize=12)
    ax_lc.legend(fontsize=10, ncol=3, loc='lower right')
    ax_lc.grid(True, alpha=0.3)

    title = f"STELLA Lightcurves — {model_name}"
    if distance_modulus != 0:
        title += f"  (μ={distance_modulus:.2f})"
    if av_extinction != 0:
        title += f"  A_V={av_extinction:.2f}"
    ax_lc.set_title(title, fontsize=13)

    # --- Bottom panel: bolometric luminosity in erg/s ---
    LSUN = 3.828e33  # erg/s per solar luminosity
    mbol = tt['Mbol'].copy() + distance_modulus
    # M_bol_sun = 4.74
    lum_erg_s = LSUN * 10.0 ** ((4.74 - mbol) / 2.5)
    
    # Plot post-SBO (positive times) with solid line
    if np.any(post_sbo_mask):
        ax_bol.semilogy(t[post_sbo_mask], lum_erg_s[post_sbo_mask], color='black', lw=1.8)
    
    # Plot pre-SBO (negative times) with dashed line
    if np.any(pre_sbo_mask):
        ax_bol.semilogy(t[pre_sbo_mask], lum_erg_s[pre_sbo_mask], color='black', ls='--', lw=1.8, alpha=0.7)
    
    if days_max is not None:
        ax_bol.set_xlim(t[0], days_max)
    else:
        ax_bol.set_xlim(t[0], t[-1] * 1.05)
    ax_bol.set_ylabel(r'$L_\mathrm{bol}\ [\mathrm{erg/s}]$', fontsize=12)
    ax_bol.set_xlabel('Time (days)', fontsize=12)
    ax_bol.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    print(f"Lightcurve figure saved: {outpath}")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('tt_file', help='STELLA *.tt lightcurve file')
    parser.add_argument('--output', default=None,
                        help='Output figure path (default: <tt_file>.png)')
    parser.add_argument('--distance', type=float, default=0.0,
                        help='Distance modulus μ = 5log(d/10pc) in mag')
    parser.add_argument('--extinction', type=float, default=0.0,
                        help='V-band extinction A_V in mag')
    parser.add_argument('--days-max', type=float, default=None,
                        help='X-axis upper limit in days')
    parser.add_argument('--bands', default='Mbol,V,B,U,R,I',
                        help='Comma-separated bands to plot')
    parser.add_argument('--no-invert', action='store_true',
                        help='Do not invert y-axis (default: brighter = up)')
    parser.add_argument('--csv', default=None,
                        help='Also write lightcurve CSV to this path')
    args = parser.parse_args()

    tt_path = Path(args.tt_file)
    if not tt_path.exists():
        print(f"ERROR: tt file not found: {tt_path}", file=sys.stderr)
        sys.exit(1)

    model_name = tt_path.stem
    outpath = args.output or str(tt_path.with_suffix('.png'))
    bands = [b.strip() for b in args.bands.split(',')]

    print(f"Reading {tt_path}...")
    tt = read_tt(str(tt_path))
    print(f"  {len(tt['time'])} epochs: t=[{tt['time'][0]:.3f}, {tt['time'][-1]:.3f}] days")
    print(f"  Mbol range: [{tt['Mbol'].min():.2f}, {tt['Mbol'].max():.2f}]")

    plot_lightcurves(
        tt, bands=bands,
        distance_modulus=args.distance,
        av_extinction=args.extinction,
        days_max=args.days_max,
        invert_y=not args.no_invert,
        model_name=model_name,
        outpath=outpath,
    )

    if args.csv:
        write_lightcurve_csv(tt, args.csv, distance_modulus=args.distance)


if __name__ == '__main__':
    main()
