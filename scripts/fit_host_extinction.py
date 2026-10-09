#!/usr/bin/env python3
"""Joint (model, host-extinction) fit for the Dessart-inspired grid.

The Milky Way component is fixed at E(B-V)=0.15 with the Gordon+2023
Rv=3.1 curve.  The host component is scanned over E(B-V)_host and
R_V,host independently - the two components are applied as separate
screenings:

    A_band = EBV_MW * coeff_G23(Rv=3.1)_band
           + EBV_host * coeff_G23(Rv=Rv_host)_band

For every grid point the observations are re-dereddened and each
model re-scored with the same best_phase_shift metric used by
plot_dessart_grid_vs_obs.py.  Output: per-model best (EBV_host,
Rv_host) and a score heatmap over the grid for the best model.

Usage:
    python scripts/fit_host_extinction.py [--n-models 8]
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import plot_dessart_grid_vs_obs as base

EBV_MW = 0.15  # S&F foreground at the SN position (IRSA)

# G23 band-averaged A/E(B-V) over the SVO SDSS g, SDSS r and ATLAS o
# transmission curves (lam_eff 4719/6185/6866 A), tabulated vs Rv.
RV_TABLE = np.array([2.3, 2.5, 2.7, 3.1, 4.0, 5.5])
A_G_TABLE = np.array([2.971, 3.172, 3.374, 3.776, 4.682, 6.191])
A_R_TABLE = np.array([1.941, 2.130, 2.318, 2.695, 3.543, 4.957])
A_O_TABLE = np.array([1.627, 1.808, 1.989, 2.350, 3.163, 4.518])

COEFF = {'g': A_G_TABLE, 'r': A_R_TABLE, 'o': A_O_TABLE}
COEFF_MW = {'g': 3.776, 'r': 2.695, 'o': 2.350}  # G23 Rv=3.1


def host_extinction(band, rv):
    """G23 A/E(B-V) for `band` at host Rv, interpolated on the table."""
    return float(np.interp(rv, RV_TABLE, COEFF[band]))


def deredden(obs0, ebv_host, rv_host):
    """Apply MW + host extinction to the unreddened obs magnitudes."""
    obs = {}
    for band, det in obs0.items():
        A = EBV_MW * COEFF_MW[band] + ebv_host * host_extinction(band, rv_host)
        d = det.copy()
        d['Mabs'] = det['Mabs'] - A
        obs[band] = d
    return obs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--grid-dir', default=base.GRID_DIR)
    ap.add_argument('--n-models', type=int, default=8,
                    help='top-N models (by MW-only score) to refit')
    ap.add_argument('--fit-end', type=float, default=100.0)
    ap.add_argument('--ebv-max', type=float, default=0.40)
    ap.add_argument('--ebv-step', type=float, default=0.025)
    args = ap.parse_args()

    import glob
    models = []
    for mdir in sorted(d for d in
                       glob.glob(os.path.join(args.grid_dir, '*'))
                       if os.path.isdir(d)):
        tt = base.find_tt(mdir)
        if tt is None:
            continue
        df = base.load_tt(tt)
        if df is not None and len(df):
            models.append((os.path.basename(mdir), df))
    if not models:
        raise SystemExit(f'no completed models in {args.grid_dir}')
    print(f'loaded {len(models)} models')

    # Observations with zero applied extinction (apparent mags minus
    # the distance modulus) - the MW+host screening is then applied
    # per grid point.
    base.set_ebv(0.0)
    obs0 = base.load_obs()

    # first pass at MW-only to pick the best models
    obs_mw = deredden(obs0, 0.0, 3.1)
    first = []
    for name, df in models:
        tpk = df.loc[df['R'].idxmin(), 'time']
        score, dt, cr, ct, cc, pen = base.best_phase_shift(
            df, obs_mw, tpk, args.fit_end)
        first.append((score, np.isfinite(pen), name, df, tpk))
    first.sort()
    top = [f for f in first if f[1]][:args.n_models]

    ebv_grid = np.arange(0.0, args.ebv_max + 1e-9, args.ebv_step)
    rv_grid = np.array([2.3, 2.5, 2.7, 3.1, 4.0, 5.5])

    print(f'\n{"model":40s} {"MW-only":>7s} -> {"best":>6s} '
          f'{"EBV_h":>6s} {"Rv":>5s}   components (rise/tail/col)')
    results = {}
    for score0, ok, name, df, tpk in top:
        best = (np.inf, None, None, None)
        grid = np.full((len(rv_grid), len(ebv_grid)), np.nan)
        for i, rv in enumerate(rv_grid):
            for j, eh in enumerate(ebv_grid):
                obs = deredden(obs0, eh, rv)
                s = base.best_phase_shift(df, obs, tpk, args.fit_end)
                grid[i, j] = s[0] if np.isfinite(s[5]) else np.nan
                if s[0] < best[0] and np.isfinite(s[5]):
                    best = (s[0], eh, rv, s)
        results[name] = (score0, best, grid)
        s = best[3]
        print(f'{name:40s} {score0:7.1f} -> {best[0]:6.1f} '
              f'{best[1]:6.3f} {best[2]:5.1f}   '
              f'({s[2]:.1f}/{s[3]:.1f}/{s[4]:.1f})')

    # heatmap of score over the (EBV_host, Rv) grid for the best model
    name0 = min(results, key=lambda n: results[n][1][0])
    grid = results[name0][2]
    fig, ax = plt.subplots(figsize=(8, 5))
    vmax = np.nanpercentile(grid, 95)
    im = ax.imshow(grid, origin='lower', aspect='auto', cmap='viridis_r',
                   extent=[ebv_grid[0], ebv_grid[-1], rv_grid[0],
                           rv_grid[-1]],
                   vmin=np.nanmin(grid), vmax=vmax)
    ax.set_xlabel('E(B-V) host')
    ax.set_ylabel('Rv host')
    ax.set_title(f'{name0}: score over host extinction grid')
    ij = np.unravel_index(np.nanargmin(grid), grid.shape)
    ax.plot(ebv_grid[ij[1]], rv_grid[ij[0]], 'r*', ms=18)
    ax.grid(True, alpha=0.3)
    plt.colorbar(im, ax=ax, label='score')
    out = os.path.join(args.grid_dir, 'host_extinction_grid.png')
    plt.tight_layout()
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'\nSaved: {out}')


if __name__ == '__main__':
    main()
