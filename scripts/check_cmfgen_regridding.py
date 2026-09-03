#!/usr/bin/env python3
"""Plot CMFGEN regridding diagnostics from R_REGRIDDING_LOG."""
import argparse, os, sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


def parse_r_regridding_log(filepath):
    """Parse R_REGRIDDING_LOG. Returns dict with old/new grid arrays."""
    lines = Path(filepath).read_text().splitlines()
    sec1, sec2 = [], []
    section = None
    for line in lines:
        s = line.strip()
        if s.startswith('!Index') and 'Log(TAU)' in s:
            section = 1; continue
        if s.startswith('!Index') and 'R' in s and 'Rold' in s:
            section = 2; continue
        if not s or s.startswith('!'): continue
        parts = s.split()
        try:
            vals = [float(p) for p in parts]
        except ValueError:
            continue
        if section == 1 and len(vals) >= 8:
            sec1.append(vals)
        elif section == 2 and len(vals) >= 9:
            sec2.append(vals)

    sec1 = np.array(sec1)
    sec2 = np.array(sec2)
    d = {
        'index': sec2[:, 0].astype(int),
        'log_tau': sec1[:, 1], 'dlog_tau': sec1[:, 2], 'ratio': sec1[:, 3],
        'log_tau_old': sec1[:, 4], 'dlog_tau_old': sec1[:, 5],
        'log_r_old_sec1': sec1[:, 6], 'v_old_sec1': sec1[:, 7],
        'r': sec2[:, 1], 'r_old': sec2[:, 2],
        'log_r': sec2[:, 3], 'log_r_old': sec2[:, 4],
        'v': sec2[:, 5], 'v_old': sec2[:, 6],
        'sigma': sec2[:, 7], 'sigma_old': sec2[:, 8],
    }
    return d


def parse_sn_hydro_data(filepath):
    """Parse SN_HYDRO_DATA. Returns dict with r, v, t, rho, ne or None."""
    p = Path(filepath)
    if not p.exists():
        return None
    lines = p.read_text().splitlines()
    nd = None
    for line in lines:
        if 'Number of data points' in line:
            nd = int(line.split(':')[1].strip())
            break
    if nd is None:
        return None

    def read_section(keyword, n):
        for i, line in enumerate(lines):
            if keyword in line:
                vals = []
                for j in range(i + 1, len(lines)):
                    vals.extend([float(x) for x in lines[j].split()])
                    if len(vals) >= n:
                        return np.array(vals[:n])
                return None
        return None

    r = read_section('Radius grid', nd)
    v = read_section('Velocity', nd)
    t = read_section('Temperature', nd)
    rho = read_section('Density', nd)
    ne = read_section('Electron density', nd)
    return {'r': r, 'v': v, 't': t, 'rho': rho, 'ne': ne, 'nd': nd}


def _interp_onto_old_grid(hydro, r_old):
    """Interpolate hydro quantities onto old grid radii.

    The regridding resamples the same physical structure at different radii,
    so interpolating the new-grid profiles onto the old-grid radii gives
    the old grid's sampling of the same profile.
    """
    r_new = hydro['r']
    result = {}
    for key in ('v', 't', 'rho', 'ne'):
        if hydro.get(key) is None:
            result[key] = None
            continue
        result[key] = np.interp(r_old, r_new, hydro[key])
    return result


def plot_regridding(data, hydro, outdir):
    """Create diagnostic plots showing full structure profiles."""
    idx = data['index']
    os.makedirs(outdir, exist_ok=True)

    r_new = data['r']
    r_old = data['r_old']
    r_label = r'R (10^{15} cm)'

    # Interpolate hydro quantities onto old grid radii
    hydro_old = None
    if hydro and hydro.get('rho') is not None and len(hydro['r']) == len(r_new):
        hydro_old = _interp_onto_old_grid(hydro, r_old)
    elif hydro and hydro.get('rho') is not None and len(hydro['r']) != len(r_new):
        print(f'Warning: R_REGRIDDING_LOG has {len(r_new)} points but '
              f'SN_HYDRO_DATA has {len(hydro["r"])} — skipping old-grid interpolation')

    # --- Figure 1: Full structure profiles vs radius ---
    n_panels = 6 if hydro and hydro.get('rho') is not None else 4
    ncols = 3
    nrows = (n_panels + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 5 * nrows),
                             squeeze=False)
    fig.suptitle('CMFGEN Regridding — Full Structure Profiles vs Radius',
                 fontsize=14)

    panel = 0

    # Radius vs depth (as markers, to show grid sampling)
    ax = axes[panel // ncols, panel % ncols]
    ax.plot(r_new / 1e15, idx, 'b-', lw=1.5, label='New')
    ax.plot(r_old / 1e15, idx, 'r--', lw=1.5, label='Old')
    ax.set_xlabel(r_label); ax.set_ylabel('Depth index')
    ax.set_title('Grid mapping (R vs depth)')
    ax.legend(); ax.invert_yaxis()
    panel += 1

    # Velocity vs radius
    ax = axes[panel // ncols, panel % ncols]
    ax.plot(r_new / 1e15, data['v'], 'b-', lw=1.5, label='New')
    ax.plot(r_old / 1e15, data['v_old'], 'r--', lw=1.5, label='Old')
    ax.set_xlabel(r_label); ax.set_ylabel('V (km/s)')
    ax.set_title('Velocity')
    ax.legend()
    panel += 1

    # Sigma vs radius
    ax = axes[panel // ncols, panel % ncols]
    ax.plot(r_new / 1e15, data['sigma'], 'b-', lw=1.5, label='New')
    ax.plot(r_old / 1e15, data['sigma_old'], 'r--', lw=1.5, label='Old')
    ax.axhline(-1.0, color='k', ls=':', alpha=0.5, label='Sigma=-1')
    ax.set_xlabel(r_label); ax.set_ylabel('Sigma (dlnv/dlnr)')
    ax.set_title('Velocity gradient (Sigma)')
    ax.legend()
    panel += 1

    # Log(Tau) vs radius
    ax = axes[panel // ncols, panel % ncols]
    ax.plot(r_new / 1e15, data['log_tau'], 'b-', lw=1.5, label='New')
    ax.plot(r_old / 1e15, data['log_tau_old'], 'r--', lw=1.5, label='Old')
    ax.set_xlabel(r_label); ax.set_ylabel('log10(Tau)')
    ax.set_title('Optical depth')
    ax.legend()
    panel += 1

    # Density vs radius (from SN_HYDRO_DATA)
    if hydro and hydro.get('rho') is not None:
        hr = hydro['r']  # use hydro's own radii
        ax = axes[panel // ncols, panel % ncols]
        ax.semilogy(hr / 1e15, hydro['rho'], 'b-', lw=1.5, label='New')
        if hydro_old and hydro_old.get('rho') is not None:
            ax.semilogy(r_old / 1e15, hydro_old['rho'], 'r--', lw=1.5,
                        label='Old (interp)')
        ax.set_xlabel(r_label); ax.set_ylabel('Density (g/cm^3)')
        ax.set_title('Density')
        ax.legend()
        panel += 1

        # Temperature vs radius
        ax = axes[panel // ncols, panel % ncols]
        ax.plot(hr / 1e15, hydro['t'], 'b-', lw=1.5, label='New')
        if hydro_old and hydro_old.get('t') is not None:
            ax.plot(r_old / 1e15, hydro_old['t'], 'r--', lw=1.5,
                    label='Old (interp)')
        ax.set_xlabel(r_label); ax.set_ylabel('T (10^4 K)')
        ax.set_title('Temperature')
        ax.legend()
        panel += 1

    # Hide unused panels
    while panel < nrows * ncols:
        axes[panel // ncols, panel % ncols].set_visible(False)
        panel += 1

    plt.tight_layout()
    p1 = os.path.join(outdir, 'regridding_profiles_vs_radius.png')
    fig.savefig(p1, dpi=150)
    plt.close(fig)
    print(f'Saved: {p1}')

    # --- Figure 2: Full structure profiles vs depth index ---
    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 5 * nrows),
                             squeeze=False)
    fig.suptitle('CMFGEN Regridding — Full Structure Profiles vs Depth Index',
                 fontsize=14)

    panel = 0

    # Radius vs depth
    ax = axes[panel // ncols, panel % ncols]
    ax.plot(idx, r_new / 1e15, 'b-', lw=1.5, label='New')
    ax.plot(idx, r_old / 1e15, 'r--', lw=1.5, label='Old')
    ax.set_xlabel('Depth index'); ax.set_ylabel(r_label)
    ax.set_title('Radius')
    ax.legend(); ax.invert_xaxis()
    panel += 1

    # Velocity vs depth
    ax = axes[panel // ncols, panel % ncols]
    ax.plot(idx, data['v'], 'b-', lw=1.5, label='New')
    ax.plot(idx, data['v_old'], 'r--', lw=1.5, label='Old')
    ax.set_xlabel('Depth index'); ax.set_ylabel('V (km/s)')
    ax.set_title('Velocity')
    ax.legend(); ax.invert_xaxis()
    panel += 1

    # Sigma vs depth
    ax = axes[panel // ncols, panel % ncols]
    ax.plot(idx, data['sigma'], 'b-', lw=1.5, label='New')
    ax.plot(idx, data['sigma_old'], 'r--', lw=1.5, label='Old')
    ax.axhline(-1.0, color='k', ls=':', alpha=0.5, label='Sigma=-1')
    ax.set_xlabel('Depth index'); ax.set_ylabel('Sigma (dlnv/dlnr)')
    ax.set_title('Velocity gradient (Sigma)')
    ax.legend(); ax.invert_xaxis()
    panel += 1

    # Log(Tau) vs depth
    ax = axes[panel // ncols, panel % ncols]
    ax.plot(idx, data['log_tau'], 'b-', lw=1.5, label='New')
    ax.plot(idx, data['log_tau_old'], 'r--', lw=1.5, label='Old')
    ax.set_xlabel('Depth index'); ax.set_ylabel('log10(Tau)')
    ax.set_title('Optical depth')
    ax.legend(); ax.invert_xaxis()
    panel += 1

    # Density vs depth
    if hydro and hydro.get('rho') is not None:
        hidx = np.arange(1, hydro['nd'] + 1)
        ax = axes[panel // ncols, panel % ncols]
        ax.semilogy(hidx, hydro['rho'], 'b-', lw=1.5, label='New')
        if hydro_old and hydro_old.get('rho') is not None:
            ax.semilogy(idx, hydro_old['rho'], 'r--', lw=1.5, label='Old (interp)')
        ax.set_xlabel('Depth index'); ax.set_ylabel('Density (g/cm^3)')
        ax.set_title('Density')
        ax.legend(); ax.invert_xaxis()
        panel += 1

        # Temperature vs depth
        ax = axes[panel // ncols, panel % ncols]
        ax.plot(hidx, hydro['t'], 'b-', lw=1.5, label='New')
        if hydro_old and hydro_old.get('t') is not None:
            ax.plot(idx, hydro_old['t'], 'r--', lw=1.5, label='Old (interp)')
        ax.set_xlabel('Depth index'); ax.set_ylabel('T (10^4 K)')
        ax.set_title('Temperature')
        ax.legend(); ax.invert_xaxis()
        panel += 1

    # Hide unused panels
    while panel < nrows * ncols:
        axes[panel // ncols, panel % ncols].set_visible(False)
        panel += 1

    plt.tight_layout()
    p2 = os.path.join(outdir, 'regridding_profiles_vs_depth.png')
    fig.savefig(p2, dpi=150)
    plt.close(fig)
    print(f'Saved: {p2}')


def main():
    ap = argparse.ArgumentParser(description='Plot CMFGEN regridding diagnostics')
    ap.add_argument('logfile', help='Path to R_REGRIDDING_LOG')
    ap.add_argument('--output-dir', default=None,
                    help='Output directory (default: same as logfile)')
    args = ap.parse_args()

    logpath = Path(args.logfile)
    if not logpath.exists():
        print(f'Error: {logpath} not found', file=sys.stderr)
        sys.exit(1)

    outdir = args.output_dir or str(logpath.parent / 'plots')
    data = parse_r_regridding_log(logpath)
    print(f'Parsed {len(data["index"])} grid points from {logpath}')

    hydro_path = logpath.parent / 'SN_HYDRO_DATA'
    hydro = parse_sn_hydro_data(hydro_path)
    if hydro:
        print(f'Parsed SN_HYDRO_DATA: {hydro["nd"]} points')
    else:
        print('No SN_HYDRO_DATA found (density/temperature plots skipped)')

    plot_regridding(data, hydro, outdir)


if __name__ == '__main__':
    main()
