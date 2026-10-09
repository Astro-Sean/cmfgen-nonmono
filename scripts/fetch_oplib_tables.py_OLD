#!/usr/bin/env python3.12
"""
Fetch OPLIB Rosseland mean opacity tables from the LANL TOPS web service
(Colgan et al. 2016, ApJ 817, 116) for the compositions present in a
STELLA input_profile and cache them as a single .npz file.

Run this ONCE on the login node before submitting a STELLA job.
The cached .npz can then be used on the compute node without any internet
access required.  (Note: generate_stella_opacity_tables.py currently uses
OPAL Rosseland mean tables via the --opal-tables path; the OPLIB cache
is a standalone output for diagnostic or future use.)

Usage:
    python3.12 fetch_oplib_tables.py <input_profile> <cache_npz> [--skip N]
                                     [--plot] [--plot-dir DIR]

Arguments:
    input_profile  STELLA input_profile (columns: r rho p u fH fHe fC fO ...)
    cache_npz      Output .npz cache file path

Options:
    --skip N       Map profile cells to N representative compositions
                   (default: same as nim = nzon//skip+1 used by STELLA,
                    but you can pass any integer; 3 is enough for X1).
    --plot         Generate diagnostic plots of the fetched opacity tables.
    --plot-dir DIR Directory for diagnostic plots (default: same dir as cache).

Output npz keys:
    zone_compositions  float32 (N_zones, 5): [xH, xHe, xC, xO, Z] per zone group
    LogT               float64 (70,): log10(T/K) grid
    LogR               float64 (19,): log10(R) grid  where R=rho/T6^3
    LogKappa           float32 (N_zones, 70, 19): log10(kappa_Ross [cm^2/g])
"""

import numpy as np
import argparse
import os
import sys
import tempfile
import shutil
import warnings

try:
    from pyTOPSScrape.api import call as tops_call
    from pyTOPSScrape.api.convert import (
        load_tops, convert_rho_2_LogR
    )
    HAS_TOPS = True
except ImportError:
    HAS_TOPS = False


# ---------------------------------------------------------------------------
# Composition helpers
# ---------------------------------------------------------------------------

def read_profile_compositions(profile_path, n_groups):
    """
    Read the STELLA input_profile and return n_groups representative
    compositions sampled uniformly along the profile.

    Returns array of shape (n_groups, 5): [xH, xHe, xC, xO, Z] per group.
    Z is computed as 1 - xH - xHe.
    """
    prof = np.loadtxt(profile_path, comments='#')
    n_cells = prof.shape[0]
    # columns: r rho p u fH fHe fC fO ...
    fH  = prof[:, 4]
    fHe = prof[:, 5] if prof.shape[1] > 5 else np.full(n_cells, 0.28)
    fC  = prof[:, 6] if prof.shape[1] > 6 else np.zeros(n_cells)
    fO  = prof[:, 7] if prof.shape[1] > 7 else np.zeros(n_cells)

    stride = max(1, n_cells // n_groups)
    comps = []
    for i in range(n_groups):
        idx = min(i * stride, n_cells - 1)
        xH  = float(fH[idx])
        xHe = float(fHe[idx])
        xC  = float(fC[idx])
        xO  = float(fO[idx])
        Z   = max(0.0, 1.0 - xH - xHe)
        comps.append([xH, xHe, xC, xO, Z])
    return np.array(comps, dtype=np.float64)


def write_abundance_file(xH, xHe, xC, xO, path):
    """
    Write a pyTOPSScrape-format abundance table file for a given composition.

    Format (from pyTOPSScrape module docs):
      Row 0: STD [Fe/H] [alpha/Fe] [C/Fe] [N/Fe] [O/Fe] [r/Fe] [s/Fe] C/O X Y,Z
      Row 1: H He Li Be B C N O F Ne   (log a(i), with H=12)
      Row 2: Na Mg Al Si P S Cl Ar K Ca
      Row 3: Sc Ti V Cr Mn Fe Co Ni Cu Zn
      (rows 4-10: heavier elements, all set to solar)

    a(i) = log10( n(i)/n(H) ) + 12
    """
    Z = max(1e-10, 1.0 - xH - xHe)
    xHe_use = max(1e-10, xHe)
    xH_use  = max(1e-10, xH)

    # Atomic masses
    am = {'H':1.008,'He':4.003,'Li':6.941,'Be':9.012,'B':10.81,
          'C':12.01,'N':14.01,'O':16.00,'F':19.00,'Ne':20.18,
          'Na':22.99,'Mg':24.31,'Al':26.98,'Si':28.09,'P':30.97,
          'S':32.07,'Cl':35.45,'Ar':39.95,'K':39.10,'Ca':40.08,
          'Sc':44.96,'Ti':47.87,'V':50.94,'Cr':52.00,'Mn':54.94,
          'Fe':55.85,'Co':58.93,'Ni':58.69,'Cu':63.55,'Zn':65.38}

    # Solar log a(i) values (Asplund 2009, with H=12)
    sol_a = {'H':12.00,'He':10.93,'Li':1.05,'Be':1.38,'B':2.70,
             'C':8.43,'N':7.83,'O':8.69,'F':4.56,'Ne':7.93,
             'Na':6.24,'Mg':7.60,'Al':6.45,'Si':7.51,'P':5.41,
             'S':7.12,'Cl':5.50,'Ar':6.40,'K':5.03,'Ca':6.34,
             'Sc':3.15,'Ti':4.95,'V':3.93,'Cr':5.64,'Mn':5.43,
             'Fe':7.50,'Co':4.99,'Ni':6.22,'Cu':4.19,'Zn':4.56}

    # Convert solar a(i) -> number fractions
    n_sol = {k: 10**(v - 12.0) for k, v in sol_a.items()}

    # Build number fractions for this composition.
    # H and He are set directly; C and O are set from input.
    # Remaining metals scale as solar Z distribution over remaining Z budget.
    n_H  = xH_use  / am['H']
    n_He = xHe_use / am['He']
    n_C  = xC      / am['C']  if xC  > 0 else 0.0
    n_O  = xO      / am['O']  if xO  > 0 else 0.0

    # Budget remaining Z after C and O
    Z_CO     = xC + xO
    Z_rest   = max(0.0, Z - Z_CO)
    sol_Z_rest = sum(n_sol[k] * am[k] for k in sol_a
                     if k not in ('H','He','C','O'))
    scale = Z_rest / sol_Z_rest if sol_Z_rest > 1e-30 else 0.0

    n = {}
    n['H']  = n_H
    n['He'] = n_He
    n['C']  = n_C
    n['O']  = n_O
    for k in sol_a:
        if k not in ('H','He','C','O'):
            n[k] = n_sol[k] * scale

    # a(i) = log10( n(i)/n(H) ) + 12
    def a(elem):
        ni = n.get(elem, 0.0)
        if ni <= 0 or n['H'] <= 0:
            return -99.00
        return np.log10(ni / n['H']) + 12.0

    # Header row: STD [Fe/H] [alpha/Fe] [C/Fe] [N/Fe] [O/Fe] [r/Fe] [s/Fe] C/O X Y,Z
    feh  = a('Fe') - sol_a['Fe']
    aFe  = 0.0
    cFe  = a('C')  - sol_a['C']
    nFe  = a('N')  - sol_a['N']
    oFe  = a('O')  - sol_a['O']
    rFe  = -99.0
    sFe  = -99.0
    CO   = (a('C') - a('O')) if (n['O'] > 0 and n['C'] > 0) else -99.0
    Y    = xHe_use
    co_ratio = xC / xO if xO > 0 else 0.0

    with open(path, 'w') as f:
        f.write(f"F {feh:.2f} {aFe:.2f} {cFe:.2f} {nFe:.2f} {oFe:.2f} "
                f"{rFe:.2f} {sFe:.2f} {CO:.2f} {xH_use:.4f} {Y:.4f},{Z:.4e}\n")
        # Row 1: H He Li Be B C N O F Ne
        row1 = ['H','He','Li','Be','B','C','N','O','F','Ne']
        f.write(' '.join(f"{a(e):.2f}" for e in row1) + '\n')
        # Row 2: Na Mg Al Si P S Cl Ar K Ca
        row2 = ['Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca']
        f.write(' '.join(f"{a(e):.2f}" for e in row2) + '\n')
        # Row 3: Sc Ti V Cr Mn Fe Co Ni Cu Zn
        row3 = ['Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn']
        f.write(' '.join(f"{a(e):.2f}" for e in row3) + '\n')
        # Rows 4-10: heavier elements (not needed by TOPS, set to -99)
        for row in [
            ['Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr'],
            ['Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn'],
            ['Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd'],
            ['Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb'],
            ['Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg'],
            ['Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th'],
            ['Pa','U'],
        ]:
            f.write(' '.join('-99.00' for _ in row) + '\n')


def write_map_file(xH, xHe, Z, path):
    """Write a single-row comma-delimited X,Y,Z abundance map."""
    with open(path, 'w') as f:
        f.write(f"{xH:.6f},{xHe:.6f},{Z:.6f}\n")


# ---------------------------------------------------------------------------
# Fetch one composition via TOPS
# ---------------------------------------------------------------------------

def fetch_one_composition(xH, xHe, xC, xO, cache_dir, idx, jobs=5):
    """
    Query LANL TOPS for a single composition (xH, xHe, xC, xO).
    Returns (LogR, LogT, LogKappa) arrays on the OPAL (19, 70) grid.
    """
    abund_path = os.path.join(cache_dir, f'comp_{idx}.abun')
    map_path   = os.path.join(cache_dir, f'map_{idx}.dat')
    raw_dir    = os.path.join(cache_dir, f'raw_{idx}')
    os.makedirs(raw_dir, exist_ok=True)

    write_abundance_file(xH, xHe, xC, xO, abund_path)

    # Single-composition map: X,Y,Z (comma-delimited)
    Z = max(0.0, 1.0 - xH - xHe)
    write_map_file(xH, xHe, Z, map_path)

    print(f"  Fetching TOPS table {idx}: xH={xH:.4f} xHe={xHe:.4f} "
          f"xC={xC:.4f} xO={xO:.4f} Z={Z:.4f}")
    tops_call(map_path, abund_path, raw_dir, jobs=jobs, opacityType='plank')

    # Find the downloaded file
    raw_files = [f for f in os.listdir(raw_dir) if f.startswith('OP:')]
    if not raw_files:
        raise RuntimeError(f"TOPS query returned no files for composition {idx}")
    raw_path = os.path.join(raw_dir, raw_files[0])

    rho, LogT, log_kappa_rho = load_tops(raw_path)
    LogR, LogT_opal, LogKappa = convert_rho_2_LogR(rho, LogT, log_kappa_rho)

    return LogR, LogT_opal, LogKappa  # (19,), (70,), (70, 19)


# ---------------------------------------------------------------------------
# Diagnostic plots
# ---------------------------------------------------------------------------

def make_diagnostic_plots(comps, LogT, LogR, LogKappa, plot_dir):
    """
    Generate diagnostic plots of the fetched OPLIB opacity tables.
    - 2D heatmap per zone group
    - kappa vs T at fixed rho for all zones
    - kappa vs rho at fixed T for all zones
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    from matplotlib.cm import ScalarMappable

    os.makedirs(plot_dir, exist_ok=True)
    n_zones = comps.shape[0]
    colors = plt.cm.plasma(np.linspace(0.1, 0.9, n_zones))

    # --- 1. 2D heatmaps per zone ---
    fig, axes = plt.subplots(1, n_zones, figsize=(5 * n_zones, 4.5),
                              sharey=True)
    if n_zones == 1:
        axes = [axes]
    vmin = np.nanpercentile(LogKappa[np.isfinite(LogKappa)], 2)
    vmax = np.nanpercentile(LogKappa[np.isfinite(LogKappa)], 98)
    for iz, ax in enumerate(axes):
        xH, xHe, xC, xO, Z = comps[iz]
        im = ax.pcolormesh(LogT, LogR,
                           LogKappa[iz].T,
                           vmin=vmin, vmax=vmax,
                           cmap='inferno', shading='auto')
        ax.set_xlabel('log T [K]')
        if iz == 0:
            ax.set_ylabel('log R')
        ax.set_title(f'Zone {iz}\n'
                     f'xH={xH:.3f} xC={xC:.3f} xO={xO:.3f}',
                     fontsize=9)
        plt.colorbar(im, ax=ax, label=r'$\log_{10}\,\kappa_{\rm Ross}$')
    fig.suptitle('OPLIB Rosseland Mean Opacity Tables (Colgan et al. 2016)',
                 fontsize=11)
    plt.tight_layout()
    out = os.path.join(plot_dir, 'oplib_heatmaps.png')
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {out}")

    # --- 2. kappa vs T at fixed rho (pick middle R index) ---
    iR_mid = len(LogR) // 2
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for iz in range(n_zones):
        xH, xHe, xC, xO, Z = comps[iz]
        kap = LogKappa[iz, :, iR_mid]
        ax.plot(LogT, kap, color=colors[iz], lw=2,
                label=f'Zone {iz}: xH={xH:.3f}, xO={xO:.3f}')
    ax.set_xlabel(r'$\log_{10}(T/\mathrm{K})$')
    ax.set_ylabel(r'$\log_{10}(\kappa_{\rm Ross}\ [\mathrm{cm}^2/\mathrm{g}])$')
    ax.set_title(f'OPLIB opacity vs T at log R = {LogR[iR_mid]:.1f}')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.axhline(np.log10(0.34), color='gray', ls='--', lw=1.5,
               label='κ_es=0.34')
    plt.tight_layout()
    out = os.path.join(plot_dir, 'oplib_kappa_vs_T.png')
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {out}")

    # --- 3. kappa vs rho at fixed T ---
    iT_phot = np.argmin(np.abs(LogT - 4.0))  # T ~ 1e4 K
    iT_hot  = np.argmin(np.abs(LogT - 5.0))  # T ~ 1e5 K
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharey=True)
    for iz in range(n_zones):
        xH, xHe, xC, xO, Z = comps[iz]
        lbl = f'Zone {iz}: xH={xH:.3f}, xO={xO:.3f}'
        axes[0].plot(LogR, LogKappa[iz, iT_phot, :],
                     color=colors[iz], lw=2, label=lbl)
        axes[1].plot(LogR, LogKappa[iz, iT_hot, :],
                     color=colors[iz], lw=2, label=lbl)
    for ax, iT, T_label in zip(axes,
                                [iT_phot, iT_hot],
                                [f'log T={LogT[iT_phot]:.2f} (photosphere)',
                                 f'log T={LogT[iT_hot]:.2f} (hot)']):
        ax.set_xlabel(r'$\log_{10}(R)$')
        ax.set_ylabel(r'$\log_{10}(\kappa_{\rm Ross})$')
        ax.set_title(T_label)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        ax.axhline(np.log10(0.34), color='gray', ls='--', lw=1.5)
    fig.suptitle('OPLIB opacity vs log R (Colgan et al. 2016)', fontsize=11)
    plt.tight_layout()
    out = os.path.join(plot_dir, 'oplib_kappa_vs_R.png')
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {out}")

    # --- 4. Zone comparison at photospheric conditions ---
    fig, ax = plt.subplots(figsize=(7, 4.5))
    kappa_phot = 10 ** LogKappa[:, iT_phot, iR_mid]
    zone_labels = [f'Zone {i}\nxH={comps[i,0]:.3f}\nxO={comps[i,3]:.3f}'
                   for i in range(n_zones)]
    bars = ax.bar(zone_labels, kappa_phot, color=colors)
    ax.axhline(0.34, color='gray', ls='--', lw=1.5, label='κ_es=0.34')
    ax.set_ylabel(r'$\kappa_{\rm Ross}\ [\mathrm{cm}^2/\mathrm{g}]$')
    ax.set_title(f'Rosseland mean opacity per zone at\n'
                 f'T=10$^{{{LogT[iT_phot]:.1f}}}$ K, log R={LogR[iR_mid]:.1f}')
    ax.legend()
    for bar, kap in zip(bars, kappa_phot):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() * 1.05,
                f'{kap:.3f}', ha='center', va='bottom', fontsize=9)
    plt.tight_layout()
    out = os.path.join(plot_dir, 'oplib_zone_comparison.png')
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {out}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('input_profile', help='STELLA input_profile path')
    parser.add_argument('cache_npz',     help='Output .npz cache file path')
    parser.add_argument('--skip', type=int, default=None,
                        help='Number of zone groups (default: auto from profile)')
    parser.add_argument('--nzon', type=int, default=500,
                        help='STELLA Mzon (used to compute default skip, default=500)')
    parser.add_argument('--jobs', type=int, default=5,
                        help='Number of parallel TOPS query threads (default=5)')
    parser.add_argument('--plot', action='store_true', default=False,
                        help='Generate diagnostic plots')
    parser.add_argument('--plot-dir', type=str, default=None,
                        help='Directory for diagnostic plots')
    args = parser.parse_args()

    if not HAS_TOPS:
        print("ERROR: pyTOPSScrape not installed. Run: pip install pyTOPSScrape")
        sys.exit(1)

    # Determine number of zone groups
    skip_val = args.skip or max(1, args.nzon // 2)
    n_groups = args.nzon // skip_val + 1
    print(f"OPLIB table fetch for {args.input_profile}")
    print(f"  Zone groups: {n_groups} (nzon={args.nzon}, skip={skip_val})")

    # Read compositions
    comps = read_profile_compositions(args.input_profile, n_groups)
    print(f"  Compositions:")
    for i, (xH, xHe, xC, xO, Z) in enumerate(comps):
        print(f"    Zone {i}: xH={xH:.4f} xHe={xHe:.4f} xC={xC:.4f} "
              f"xO={xO:.4f} Z={Z:.4f}")

    # Fetch from TOPS
    with tempfile.TemporaryDirectory(prefix='tops_fetch_') as tmpdir:
        LogR_ref = None
        LogT_ref = None
        all_LogKappa = []

        for i, (xH, xHe, xC, xO, Z) in enumerate(comps):
            LogR, LogT, LogKappa = fetch_one_composition(
                xH, xHe, xC, xO, tmpdir, i, jobs=args.jobs)
            if LogR_ref is None:
                LogR_ref = LogR
                LogT_ref = LogT
            all_LogKappa.append(LogKappa)  # (70, 19)

    LogKappa_arr = np.array(all_LogKappa, dtype=np.float32)  # (N, 70, 19)

    # Save cache
    os.makedirs(os.path.dirname(os.path.abspath(args.cache_npz)), exist_ok=True)
    np.savez(args.cache_npz,
             zone_compositions=comps.astype(np.float32),
             LogT=LogT_ref,
             LogR=LogR_ref,
             LogKappa=LogKappa_arr)
    print(f"\nSaved OPLIB cache: {args.cache_npz}")
    print(f"  LogT: {LogT_ref.min():.3f} -- {LogT_ref.max():.3f} ({len(LogT_ref)} points)")
    print(f"  LogR: {LogR_ref.min():.1f} -- {LogR_ref.max():.1f} ({len(LogR_ref)} points)")
    print(f"  LogKappa shape: {LogKappa_arr.shape}")
    kmin = np.nanmin(LogKappa_arr[np.isfinite(LogKappa_arr)])
    kmax = np.nanmax(LogKappa_arr[np.isfinite(LogKappa_arr)])
    print(f"  log10(kappa) range: {kmin:.2f} -- {kmax:.2f}")

    # Diagnostic plots
    if args.plot:
        plot_dir = args.plot_dir or os.path.dirname(os.path.abspath(args.cache_npz))
        print(f"\nGenerating diagnostic plots in {plot_dir}")
        make_diagnostic_plots(comps, LogT_ref, LogR_ref, LogKappa_arr, plot_dir)

    print("\nDone. Use --oplib-cache in generate_stella_opacity_tables.py.")


if __name__ == '__main__':
    main()
