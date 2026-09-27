#!/usr/bin/env python3
"""
Compare parallel STELLA test runs: baseline (Y12M) vs PARDISO (1/2/3 CPUs).

Produces comparison plots:
  1. Lightcurve (bolometric luminosity) vs time from .swd flux column
  2. Density profile vs radius at common epochs
  3. Velocity profile vs radius at common epochs
  4. Runtime comparison

Usage:
  python scripts/plot_parallel_stella_comparison.py
  python scripts/plot_parallel_stella_comparison.py \
      --models "Y12M=parallel_test_baseline" "PAR1=parallel_test_p1" \
      --out parallel_test_plots --baseline "Y12M"
"""

import argparse
import sys
import os
import re
import importlib.util
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent.parent
spec = importlib.util.spec_from_file_location(
    "conv", PROJECT_DIR / "scripts" / "convert_stella_swd_to_cmfgen.py")
conv = importlib.util.module_from_spec(spec)
spec.loader.exec_module(conv)


def find_swd(model_dir):
    """Find .swd file in model directory."""
    for pattern in ['stella/*.swd', 'stella/**/*.swd']:
        files = list(Path(model_dir).glob(pattern))
        if files:
            return files[0]
    return None


def find_job_id(model_dir):
    """Find SLURM job ID from output file."""
    out_files = list(Path(model_dir).glob('stella/stella_*.out'))
    for f in out_files:
        m = re.search(r'_(\d+)\.out$', f.name)
        if m:
            return m.group(1)
    return None


def get_runtime_sacct(job_id):
    """Get elapsed time from sacct."""
    import subprocess
    try:
        result = subprocess.run(
            ['sacct', '-j', str(job_id), '--format=Elapsed', '-n'],
            capture_output=True, text=True, timeout=10)
        for line in result.stdout.strip().split('\n'):
            line = line.strip()
            if line and ':' in line:
                parts = line.split(':')
                if len(parts) == 3:
                    return int(parts[0]) * 3600 + int(parts[1]) * 60 + int(parts[2])
    except Exception:
        pass
    return None


def read_swd_lightcurve(swd_path):
    """Read bolometric luminosity from .swd flux column (already L in 1e40 erg/s)."""
    epochs = conv.read_swd(str(swd_path))
    times = []
    lum = []
    for t in sorted(epochs.keys()):
        d = epochs[t]
        if d is None:
            continue
        times.append(t)
        # flux column = luminosity in 1e40 erg/s
        lum.append(d['flux'][0])  # outermost zone
    return np.array(times), np.array(lum)


def read_swd_profile(swd_path, target_days=None):
    """Read .swd at a specific epoch (or last)."""
    epochs = conv.read_swd(str(swd_path), target_days=target_days)
    valid = [t for t in sorted(epochs.keys()) if epochs[t] is not None]
    if not valid:
        return None
    t = valid[-1]
    d = epochs[t]
    return {
        'time': t,
        'r': d['r'],
        'rho': d['rho'],
        'v': d['v'],
        'T': d['T'],
        'T_rad': d['T_rad'],
    }


def read_all_epochs(swd_path):
    """Return {time: profile_dict} for all valid epochs."""
    epochs = conv.read_swd(str(swd_path))
    out = {}
    for t in sorted(epochs.keys()):
        d = epochs[t]
        if d is None:
            continue
        out[t] = {'r': d['r'], 'rho': d['rho'], 'v': d['v'],
                  'T': d['T'], 'T_rad': d['T_rad']}
    return out


def nearest_time(times, target):
    """Return the time in `times` closest to target."""
    times = np.asarray(times)
    return times[np.argmin(np.abs(times - target))]


DEFAULT_MODELS = {
    'Y12M (serial)':   {'dir': 'Stella_parrallel/spar_base',
                        'color': 'black', 'linestyle': '-'},
    'PARDISO (1 CPU)': {'dir': 'Stella_parrallel/spar_p1',
                        'color': 'blue', 'linestyle': '--'},
    'PARDISO (2 CPU)': {'dir': 'Stella_parrallel/spar_p2',
                        'color': 'red', 'linestyle': '-.'},
    'PARDISO (3 CPU)': {'dir': 'Stella_parrallel/spar_p3',
                        'color': 'green', 'linestyle': ':'},
}
DEFAULT_OUT = 'Stella_parrallel/plots'


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--models', nargs='*', default=None,
                    help='label=dir pairs; defaults to the Stella_parrallel set')
    ap.add_argument('--baseline', default=None,
                    help='label of the baseline model for matched-time plots '
                         '(default: first model)')
    ap.add_argument('--out', default=DEFAULT_OUT,
                    help='output directory for plots, relative to project dir')
    args = ap.parse_args()

    models = dict(DEFAULT_MODELS)
    if args.models:
        models = {}
        for spec_ in args.models:
            label, _, d = spec_.partition('=')
            models[label] = {'dir': d, 'color': None, 'linestyle': '-'}

    palette = ['black', 'blue', 'red', 'green', 'purple', 'orange']
    for i, info in enumerate(models.values()):
        info['dir'] = PROJECT_DIR / info['dir']
        info['color'] = info['color'] or palette[i % len(palette)]

    base_label = args.baseline or next(iter(models))

    # Find output files
    swd_paths = {}
    job_ids = {}
    for label, info in models.items():
        model_dir = info['dir']
        swd = find_swd(model_dir)
        jid = find_job_id(model_dir)
        if swd:
            swd_paths[label] = swd
            print(f"  {label}: {swd}")
        if jid:
            job_ids[label] = jid
        if not swd:
            print(f"  {label}: NO .swd found")

    if not swd_paths:
        print("ERROR: No .swd files found.")
        sys.exit(1)

    plots_dir = PROJECT_DIR / args.out
    plots_dir.mkdir(parents=True, exist_ok=True)

    # ============================================================
    # Plot 1: Bolometric lightcurve
    # ============================================================
    print("\n--- Plot 1: Bolometric lightcurve ---")
    fig, ax = plt.subplots(figsize=(10, 6))
    any_data = False
    for label, info in models.items():
        if label not in swd_paths:
            continue
        try:
            times, lum = read_swd_lightcurve(str(swd_paths[label]))
            if len(times) > 0:
                ax.plot(times, lum, label=label, color=info['color'],
                        linestyle=info['linestyle'], linewidth=1.0, alpha=0.8)
                any_data = True
                print(f"  {label}: {len(times)} epochs, t={times[0]:.1f}-{times[-1]:.1f}d, "
                      f"L_last={lum[-1]:.2e} erg/s")
        except Exception as e:
            print(f"  {label}: error reading lightcurve: {e}")
    if any_data:
        ax.set_xlabel('Time (days)')
        ax.set_ylabel('Bolometric Luminosity (10$^{40}$ erg/s)')
        ax.set_yscale('log')
        ax.set_title('STELLA Bolometric Lightcurve: Y12M vs PARDISO')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(plots_dir / 'lightcurve_comparison.png', dpi=150)
        print(f"  Saved: {plots_dir / 'lightcurve_comparison.png'}")
    plt.close(fig)

    # ============================================================
    # Plot 2: Density profile at final epoch
    # ============================================================
    print("\n--- Plot 2: Density profile ---")
    fig, ax = plt.subplots(figsize=(10, 6))
    any_data = False
    for label, info in models.items():
        if label not in swd_paths:
            continue
        prof = read_swd_profile(str(swd_paths[label]))
        if prof is None:
            continue
        ax.plot(prof['r'] / 1e15, prof['rho'],
                label=f"{label} (t={prof['time']:.2f}d)",
                color=info['color'], linestyle=info['linestyle'], linewidth=1.5)
        any_data = True
        print(f"  {label}: t={prof['time']:.2f}d, rho_range=[{prof['rho'].min():.1e},{prof['rho'].max():.1e}]")
    if any_data:
        ax.set_xlabel('Radius (10$^{15}$ cm)')
        ax.set_ylabel('Density (g/cm$^3$)')
        ax.set_yscale('log')
        ax.set_title('Density Profile at Final Epoch')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(plots_dir / 'density_comparison.png', dpi=150)
        print(f"  Saved: {plots_dir / 'density_comparison.png'}")
    plt.close(fig)

    # ============================================================
    # Plot 3: Velocity profile at final epoch
    # ============================================================
    print("\n--- Plot 3: Velocity profile ---")
    fig, ax = plt.subplots(figsize=(10, 6))
    any_data = False
    for label, info in models.items():
        if label not in swd_paths:
            continue
        prof = read_swd_profile(str(swd_paths[label]))
        if prof is None:
            continue
        ax.plot(prof['r'] / 1e15, prof['v'],
                label=f"{label} (t={prof['time']:.2f}d)",
                color=info['color'], linestyle=info['linestyle'], linewidth=1.5)
        any_data = True
        print(f"  {label}: t={prof['time']:.2f}d, v_range=[{prof['v'].min():.0f},{prof['v'].max():.0f}] km/s")
    if any_data:
        ax.set_xlabel('Radius (10$^{15}$ cm)')
        ax.set_ylabel('Velocity (km/s)')
        ax.set_title('Velocity Profile at Final Epoch')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(plots_dir / 'velocity_comparison.png', dpi=150)
        print(f"  Saved: {plots_dir / 'velocity_comparison.png'}")
    plt.close(fig)

    # ============================================================
    # Plot 4: Matched-time agreement (baseline vs each PARDISO run)
    # This is the real equivalence test: same epoch, compare rho & v.
    # ============================================================
    print("\n--- Plot 4: Matched-time agreement (baseline vs PARDISO) ---")
    if base_label in swd_paths:
        base_ep = read_all_epochs(str(swd_paths[base_label]))
        base_times = sorted(base_ep.keys())
        for label, info in models.items():
            if label == base_label or label not in swd_paths:
                continue
            test_ep = read_all_epochs(str(swd_paths[label]))
            test_times = sorted(test_ep.keys())
            if not test_times:
                continue
            # Common overlap window
            t_lo = max(min(base_times), min(test_times))
            t_hi = min(max(base_times), max(test_times))
            overlap = [t for t in base_times if t_lo <= t <= t_hi]
            if not overlap:
                print(f"  {label}: no time overlap with baseline")
                continue
            # Pick up to 4 matched epochs
            pick = np.linspace(0, len(overlap)-1, min(4, len(overlap)), dtype=int)
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
            for idx in pick:
                tb = overlap[idx]
                tt = nearest_time(test_times, tb)
                db, dt = base_ep[tb], test_ep[tt]
                a = 0.35 + 0.65*(list(pick).index(idx)/max(1, len(pick)-1))
                ax1.plot(db['r']/1e15, db['rho'], 'k-', alpha=a, lw=1.5,
                         label=f'Y12M t={tb:.2f}')
                ax1.plot(dt['r']/1e15, dt['rho'], 'r--', alpha=a, lw=1.0,
                         label=f'PAR t={tt:.2f}')
                ax2.plot(db['r']/1e15, db['v'], 'k-', alpha=a, lw=1.5)
                ax2.plot(dt['r']/1e15, dt['v'], 'r--', alpha=a, lw=1.0)
            ax1.set_yscale('log'); ax1.set_xlabel('Radius (1e15 cm)')
            ax1.set_ylabel('Density'); ax1.set_title('Density: Y12M(-) vs PARDISO(--)')
            ax1.legend(fontsize=6); ax1.grid(alpha=0.3)
            ax2.set_xlabel('Radius (1e15 cm)'); ax2.set_ylabel('Velocity (km/s)')
            ax2.set_title('Velocity: Y12M(-) vs PARDISO(--)'); ax2.grid(alpha=0.3)
            safe = re.sub(r'[^A-Za-z0-9]+', '_', label)
            fig.tight_layout(); fig.savefig(plots_dir/f'matched_{safe}.png', dpi=150)
            plt.close(fig)
            # Quantitative agreement at latest common epoch
            tb = overlap[-1]; tt = nearest_time(test_times, tb)
            db, dt = base_ep[tb], test_ep[tt]
            nb = min(len(db['rho']), len(dt['rho']))
            rho_rel = np.abs(db['rho'][:nb]-dt['rho'][:nb])/(np.abs(db['rho'][:nb])+1e-300)
            v_rel = np.abs(db['v'][:nb]-dt['v'][:nb])/(np.abs(db['v'][:nb])+1e-300)
            print(f"  {label}: {len(overlap)} common epochs, t=[{t_lo:.2f},{t_hi:.2f}]")
            print(f"     latest t={tb:.2f}: max rel diff rho={rho_rel.max():.2e}, v={v_rel.max():.2e}")

    # ============================================================
    # Plot 5: Density & velocity evolution
    # ============================================================
    print("\n--- Plot 5: Density & velocity evolution ---")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    any_data = False
    for label, info in models.items():
        if label not in swd_paths:
            continue
        try:
            epochs = conv.read_swd(str(swd_paths[label]))
        except Exception as e:
            print(f"  {label}: {e}")
            continue
        times = sorted(t for t in epochs.keys() if epochs[t] is not None)
        if len(times) < 2:
            continue
        # Pick ~4 evenly spaced epochs
        indices = np.linspace(0, len(times) - 1, min(5, len(times)), dtype=int)
        for ii, idx in enumerate(indices):
            t = times[idx]
            d = epochs[t]
            alpha = 0.3 + 0.7 * (ii / max(1, len(indices) - 1))
            ax1.plot(d['r'] / 1e15, d['rho'], color=info['color'],
                     linestyle=info['linestyle'], alpha=alpha, linewidth=0.8)
            ax2.plot(d['r'] / 1e15, d['v'], color=info['color'],
                     linestyle=info['linestyle'], alpha=alpha, linewidth=0.8)
        ax1.plot([], [], color=info['color'], linestyle=info['linestyle'],
                 label=label)
        ax2.plot([], [], color=info['color'], linestyle=info['linestyle'],
                 label=label)
        any_data = True
    if any_data:
        ax1.set_xlabel('Radius (10$^{15}$ cm)')
        ax1.set_ylabel('Density (g/cm$^3$)')
        ax1.set_yscale('log')
        ax1.set_title('Density Evolution')
        ax1.legend(fontsize=8)
        ax1.grid(True, alpha=0.3)

        ax2.set_xlabel('Radius (10$^{15}$ cm)')
        ax2.set_ylabel('Velocity (km/s)')
        ax2.set_title('Velocity Evolution')
        ax2.legend(fontsize=8)
        ax2.grid(True, alpha=0.3)

        fig.tight_layout()
        fig.savefig(plots_dir / 'evolution_comparison.png', dpi=150)
        print(f"  Saved: {plots_dir / 'evolution_comparison.png'}")
    plt.close(fig)

    # ============================================================
    # Plot 6: Runtime comparison
    # ============================================================
    print("\n--- Plot 6: Runtime ---")
    labels_rt = []
    runtimes = []
    colors_rt = []
    for label, info in models.items():
        jid = job_ids.get(label)
        if jid:
            rt = get_runtime_sacct(jid)
            if rt:
                labels_rt.append(label)
                runtimes.append(rt)
                colors_rt.append(info['color'])
                print(f"  {label}: {rt}s (job {jid})")

    if labels_rt:
        fig, ax = plt.subplots(figsize=(8, 5))
        bars = ax.bar(range(len(labels_rt)), runtimes, color=colors_rt)
        ax.set_xticks(range(len(labels_rt)))
        ax.set_xticklabels(labels_rt, fontsize=8, rotation=15, ha='right')
        ax.set_ylabel('Runtime (seconds)')
        ax.set_title('STELLA Runtime Comparison')
        ax.grid(True, alpha=0.3, axis='y')
        for bar, rt in zip(bars, runtimes):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                    f'{rt/60:.1f}m', ha='center', va='bottom')
        fig.tight_layout()
        fig.savefig(plots_dir / 'runtime_comparison.png', dpi=150)
        print(f"  Saved: {plots_dir / 'runtime_comparison.png'}")
    else:
        print("  No runtime data (jobs still running or sacct unavailable)")
    plt.close(fig)

    # ============================================================
    # Summary table
    # ============================================================
    print("\n=== Summary ===")
    print(f"{'Model':<25} {'Epochs':>7} {'t_last':>10} {'Runtime':>10} {'Job':>10}")
    print("-" * 65)
    for label, info in models.items():
        swd = swd_paths.get(label)
        jid = job_ids.get(label)
        n_ep = 0
        last_t = 'N/A'
        if swd:
            try:
                epochs = conv.read_swd(str(swd))
                valid = [t for t in epochs if epochs[t] is not None]
                n_ep = len(valid)
                if valid:
                    last_t = f"{max(valid):.2f}d"
            except:
                pass
        rt = get_runtime_sacct(jid) if jid else None
        rt_str = f"{rt/60:.1f}m" if rt else 'N/A'
        print(f"{label:<25} {n_ep:>7} {last_t:>10} {rt_str:>10} {jid or 'N/A':>10}")

    print(f"\nAll plots saved to: {plots_dir}")


if __name__ == '__main__':
    main()
