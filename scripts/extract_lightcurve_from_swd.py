#!/usr/bin/env python3
"""Extract a complete lightcurve from a STELLA *.swd snapshot file.

STELLA's native *.tt lightcurve file often has large gaps (tens of days)
for low-energy explosions because the internal tt4strad thinning routine
only writes when the emergent flux changes by >0.1% or when >1 day has
passed. During shock propagation the flux is trapped, so nothing is written
until shock breakout.

This script reads the *.swd snapshots (which DO contain complete time
coverage) and recomputes the photosphere and blackbody luminosity at
each epoch, producing a lightcurve CSV with no gaps.

Usage:
    python extract_lightcurve_from_swd.py SWD_FILE [--abn-file ABN] \
        [--output OUT.csv] [--kappa K]

Output columns:
    time_days, Tbb_K, rbb_cm, Teff_K, Rlast_sc_cm, Rtau23_cm,
    Mbol, L_bol_erg_s
"""

import argparse
import sys
from pathlib import Path

import numpy as np

SIGMA_SB = 5.670374e-5  # erg cm^-2 s^-1 K^-4
LSUN = 3.828e33         # erg/s
MBOL_SUN = 4.74


def read_swd(swd_path: str) -> dict[float, dict]:
    """Parse a STELLA *.swd file.

    Returns a dict keyed by t_day -> {'r': [], 'v_kms': [], 'rho': [], 'T': [], 'flux': [], 'km': []}.

    'flux' is the STELLA Flum column (swd col index 11, written as
    Flum*1e-40 in lbalsw.f) -- the local luminosity in erg/s, NOT a
    flux per unit area.  Col 12 is WRKX (Rosseland opacity), not flux.
    All values are numpy arrays sorted outer->inner (decreasing r).
    """
    epochs: dict[float, dict] = {}
    current_t = None
    buf_km = []
    buf_r = []
    buf_v = []
    buf_T = []
    buf_rho = []
    buf_flux = []

    def _flush():
        nonlocal current_t, buf_km, buf_r, buf_v, buf_T, buf_rho, buf_flux
        if current_t is not None and buf_r:
            order = np.argsort(buf_r)[::-1]  # outer->inner
            epochs[current_t] = {
                'km':   np.array(buf_km)[order],
                'r':    np.array(buf_r)[order],
                'v':    np.array(buf_v)[order],
                'T':    np.array(buf_T)[order],
                'rho':  np.array(buf_rho)[order],
                'flux': np.array(buf_flux)[order],
            }
        buf_km   = []
        buf_r    = []
        buf_v    = []
        buf_T    = []
        buf_rho  = []
        buf_flux = []

    with open(swd_path) as fh:
        for raw_line in fh:
            line = raw_line.rstrip()
            if not line.strip():
                continue
            cols = line.split()
            try:
                t_candidate = float(cols[0])
                km          = int(cols[1])
                log_r       = float(cols[3])
                v_kms       = float(cols[4]) * 1000.0  # .swd unit: 1e8 cm/s -> km/s
                log_T       = float(cols[5])
                log_rho     = float(cols[7])
                flux_1e40   = float(cols[11])  # Flum: luminosity in units of 1e40 erg/s
            except (ValueError, IndexError):
                continue

            # In STELLA .swd, the first zone line of each epoch carries the
            # time in column 0; continuation zone lines carry 0 there.
            # Some .swd variants repeat the time on every line. Treat a
            # nonzero time that differs from the current epoch as a new
            # epoch, and 0 as a continuation line.
            if t_candidate != 0.0 and (current_t is None or abs(t_candidate - current_t) > 1e-6):
                _flush()
                current_t = t_candidate
            elif current_t is None:
                current_t = t_candidate

            r_cm  = 10.0 ** log_r
            T_K   = 10.0 ** log_T
            # STELLA SWD uses internal units: rho is log10(rho / 1e-6 g/cc)
            # Convert to cgs: rho_cgs = 10^log_rho * 1e-6
            rho   = (10.0 ** log_rho) * 1e-6
            flux  = flux_1e40 * 1e40  # luminosity L(r) in erg/s

            buf_km.append(km)
            buf_r.append(r_cm)
            buf_v.append(v_kms)
            buf_T.append(T_K)
            buf_rho.append(rho)
            buf_flux.append(flux)

    _flush()
    if not epochs:
        raise ValueError(f"No epochs found in {swd_path}")
    return epochs


def compute_kappa_tau(r_cm: np.ndarray, rho: np.ndarray, kappa_cm2_per_g: float) -> np.ndarray:
    """Compute cumulative optical depth from outside in (tau=0 at surface)."""
    tau = np.zeros(len(r_cm))
    for i in range(len(r_cm) - 1, 0, -1):
        rho_avg = 0.5 * (rho[i] + rho[i - 1])
        dr = abs(r_cm[i] - r_cm[i - 1])
        dtau = kappa_cm2_per_g * rho_avg * dr
        tau[i - 1] = tau[i] + dtau
    return tau


def extract_lightcurve(swd_path: str, kappa: float = 0.34, tau_min: float = None) -> dict:
    """Extract lightcurve data from .swd snapshots.

    Returns dict with arrays: time, Tbb, rbb, Teff, Rlast_sc, Rtau23, Mbol, L_bol
    
    Args:
        swd_path: Path to .swd file
        kappa: Electron scattering opacity in cm^2/g
        tau_min: Optional minimum optical depth for zone filtering (external trimming)
                 Zones with tau < tau_min are excluded from analysis
    """
    epochs = read_swd(swd_path)
    times = sorted(epochs.keys())

    results = {
        'time': [],
        'Tbb': [],
        'rbb': [],
        'Teff': [],
        'Rlast_sc': [],
        'Rtau23': [],
        'Mbol': [],
        'L_bol': [],
    }

    for t in times:
        snap = epochs[t]
        r = snap['r']
        rho = snap['rho']
        T = snap['T']
        flux = snap['flux']

        # Photosphere at tau = 2/3
        tau = compute_kappa_tau(r, rho, kappa)
        
        # External trimming: filter zones with tau < tau_min
        if tau_min is not None:
            mask = tau >= tau_min
            r = r[mask]
            rho = rho[mask]
            T = T[mask]
            flux = flux[mask]
            tau = tau[mask]
            
            # Recompute tau after filtering
            tau = compute_kappa_tau(r, rho, kappa)
        
        # Find photosphere after any filtering
        i_phot = np.searchsorted(tau[::-1], 2.0 / 3.0)
        i_phot = len(tau) - 1 - i_phot
        i_phot = max(0, min(i_phot, len(r) - 1))

        r_phot = r[i_phot]
        T_phot = T[i_phot]

        # Blackbody luminosity at photosphere
        L_bol = 4.0 * np.pi * r_phot**2 * SIGMA_SB * T_phot**4
        mbol = MBOL_SUN - 2.5 * np.log10(L_bol / LSUN)

        # Rlast_sc: outermost radius where scattering dominates
        # Approximate: outermost zone with significant flux
        # For simplicity, use the outermost zone radius
        r_last_sc = r[0]

        # Teff: effective temperature from the emergent luminosity at
        # the outer boundary: F_surf = L_out / (4*pi*r_out^2),
        # T_eff = (F_surf/sigma)^0.25
        l_surf = flux[0]
        if l_surf > 0:
            teff = (l_surf / (4.0 * np.pi * r[0] ** 2 * SIGMA_SB)) ** 0.25
        else:
            teff = T[0]

        results['time'].append(t)
        results['Tbb'].append(T_phot)
        results['rbb'].append(r_phot)
        results['Teff'].append(teff)
        results['Rlast_sc'].append(r_last_sc)
        results['Rtau23'].append(r_phot)
        results['Mbol'].append(mbol)
        results['L_bol'].append(L_bol)

    for k in results:
        results[k] = np.array(results[k])

    return results


def write_csv(lc: dict, outpath: str):
    """Write lightcurve to CSV."""
    with open(outpath, 'w') as f:
        f.write("time_days,Tbb_K,rbb_cm,Teff_K,Rlast_sc_cm,Rtau23_cm,Mbol,L_bol_erg_s\n")
        for i in range(len(lc['time'])):
            f.write(f"{lc['time'][i]:.10f},{lc['Tbb'][i]:.3f},{lc['rbb'][i]:.6e},"
                    f"{lc['Teff'][i]:.3f},{lc['Rlast_sc'][i]:.6e},{lc['Rtau23'][i]:.6e},"
                    f"{lc['Mbol'][i]:.4f},{lc['L_bol'][i]:.6e}\n")
    print(f"Lightcurve CSV saved: {outpath}")


def plot_lightcurve(lc: dict, outpath: str):
    """Quick diagnostic plot of bolometric lightcurve."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 1, figsize=(10, 8))

    t = lc['time']
    mbol = lc['Mbol']
    lum = lc['L_bol']

    # Top: magnitude vs time
    ax = axes[0]
    ax.plot(t, mbol, 'ko-', markersize=3, lw=1)
    ax.invert_yaxis()
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Mbol')
    ax.set_title('Bolometric Lightcurve from .swd snapshots')
    ax.grid(True, alpha=0.3)

    # Bottom: luminosity vs time
    ax = axes[1]
    ax.semilogy(t, lum, 'ko-', markersize=3, lw=1)
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('L_bol (erg/s)')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    print(f"Lightcurve plot saved: {outpath}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('swd_file', help='STELLA *.swd snapshot file')
    parser.add_argument('--kappa', type=float, default=0.34,
                        help='Electron scattering opacity in cm^2/g (default: 0.34)')
    parser.add_argument('--tau-min', type=float, default=None,
                        help='Minimum optical depth for zone filtering (external trimming, default: None)')
    parser.add_argument('--output', default=None,
                        help='Output CSV path (default: <swd_file>_lightcurve.csv)')
    parser.add_argument('--plot', default=None,
                        help='Optional output PNG plot path')
    args = parser.parse_args()

    swd = Path(args.swd_file)
    if not swd.exists():
        print(f"ERROR: {swd} not found")
        sys.exit(1)

    print(f"Reading {swd} ...")
    lc = extract_lightcurve(str(swd), kappa=args.kappa, tau_min=args.tau_min)
    print(f"Extracted {len(lc['time'])} epochs from {lc['time'].min():.3f} to {lc['time'].max():.3f} days")
    if args.tau_min is not None:
        print(f"External trimming: zones with tau < {args.tau_min} excluded")

    # Write CSV
    out_csv = args.output or str(swd.with_suffix('.swd_lightcurve.csv'))
    write_csv(lc, out_csv)

    # Optional plot
    if args.plot:
        plot_lightcurve(lc, args.plot)


if __name__ == '__main__':
    main()
