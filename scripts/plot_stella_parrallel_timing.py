#!/usr/bin/env python3
"""
Compare wall-clock time to reach a target simulation day across the
Stella_parrallel benchmark runs (serial Y12M vs PARDISO at 1/2/3 CPUs).

Timing comes from STELLA's own log: each `time=HHMMSS.SSS` stamp is a
REAL_CLOCK print that precedes an `ipout` block writing one .swd record
(counts match ~1:1 with .swd epochs). Stamp i pairs with epoch i by
order; a few records of index drift are possible and noted on the plot.

Elapsed time is measured from the FIRST stamp, so input preparation
(eve2, modmake, ronfict) is excluded — this is solver wall time.
sacct Elapsed (whole job, prep included) is shown as a second series
when available.

Usage:
  python scripts/plot_stella_parrallel_timing.py
  python scripts/plot_stella_parrallel_timing.py --target 8.0 \
      --models "Y12M=parallel_test_baseline" "PAR1=parallel_test_p1"
"""

import argparse
import importlib.util
import re
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

PROJECT_DIR = Path(__file__).resolve().parent.parent
spec = importlib.util.spec_from_file_location(
    "conv", PROJECT_DIR / "scripts" / "convert_stella_swd_to_cmfgen.py")
conv = importlib.util.module_from_spec(spec)
spec.loader.exec_module(conv)

DEFAULT_MODELS = {
    'Y12M (serial)':      {'dir': 'Stella_parrallel/spar_base', 'color': 'black',  'ls': '-'},
    'PARDISO (1 CPU)':    {'dir': 'Stella_parrallel/spar_p1',   'color': 'blue',   'ls': '--'},
    'PARDISO (2 CPU)':    {'dir': 'Stella_parrallel/spar_p2',   'color': 'red',    'ls': '-.'},
    'PARDISO (3 CPU)':    {'dir': 'Stella_parrallel/spar_p3',   'color': 'green',  'ls': ':'},
}
DEFAULT_TARGET = 10.0
DEFAULT_OUT = 'Stella_parrallel/plots/time_to_target.png'


def find_log(model_dir):
    """Latest stella_<model>_<jobid>.out in the model dir (excludes -debug)."""
    d = Path(model_dir) / 'stella'
    files = [f for f in d.glob('stella_*_*.out') if '-debug' not in f.name]
    files = [f for f in files if re.search(r'_\d+\.out$', f.name)]
    if not files:
        return None, None
    f = max(files, key=lambda p: p.stat().st_mtime)
    return f, re.search(r'_(\d+)\.out$', f.name).group(1)


def parse_time_stamps(out_path):
    """Wall-clock seconds-of-day for every `time=` stamp, midnight-unwrapped."""
    raw = subprocess.run(['strings', str(out_path)],
                         capture_output=True, text=True, check=True).stdout
    stamps = []
    for m in re.finditer(r'time=(\d{2})(\d{2})(\d{2}\.\d+)', raw):
        stamps.append(int(m.group(1)) * 3600 + int(m.group(2)) * 60
                      + float(m.group(3)))
    # jobs crossing midnight would go backwards without the unwrap
    for i in range(1, len(stamps)):
        while stamps[i] < stamps[i - 1]:
            stamps[i] += 86400.0
    return np.array(stamps)


def read_epochs(model_dir):
    swds = sorted(Path(model_dir).glob('stella/*.swd'),
                  key=lambda p: p.stat().st_mtime)
    if not swds:
        return None
    epochs = conv.read_swd(str(swds[-1]))
    return np.array(sorted(t for t in epochs if epochs[t] is not None))


def job_elapsed_sacct(job_id):
    try:
        out = subprocess.run(
            ['sacct', '-j', str(job_id), '--format=Elapsed', '-n', '-P'],
            capture_output=True, text=True, timeout=15).stdout
    except Exception:
        return None
    for line in out.splitlines():
        line = line.strip().strip('|')
        m = re.match(r'(?:(\d+)-)?(\d+):(\d+):(\d+)$', line)
        if m:
            days = int(m.group(1) or 0)
            return days * 86400 + int(m.group(2)) * 3600 \
                + int(m.group(3)) * 60 + int(m.group(4))
    return None


def collect(model_dir):
    """Return (wall_seconds, epochs, job_id, job_elapsed) or None."""
    out_file, job_id = find_log(model_dir)
    if out_file is None:
        return None
    stamps = parse_time_stamps(out_file)
    if len(stamps) == 0:
        return None
    epochs = read_epochs(model_dir)
    elapsed = job_elapsed_sacct(job_id) if job_id else None
    return stamps, epochs, job_id, elapsed


def time_to_target(stamps, epochs, target):
    """Wall seconds from first stamp to the stamp paired with epoch >= target.

    Pairs are index-aligned (stamp i <-> epoch i), truncated to the
    shorter side. Returns (seconds, epoch_reached, reached_target).
    """
    if epochs is None or len(epochs) == 0:
        return None, None, False
    n = min(len(stamps), len(epochs))
    stamps, epochs = stamps[:n], epochs[:n]
    hit = np.nonzero(epochs >= target)[0]
    if len(hit):
        j = hit[0]
        return float(stamps[j] - stamps[0]), float(epochs[j]), True
    return float(stamps[-1] - stamps[0]), float(epochs[-1]), False


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--target', type=float, default=DEFAULT_TARGET,
                    help='simulation day the runs must reach (default 10)')
    ap.add_argument('--models', nargs='*', default=None,
                    help='label=dir pairs; defaults to the Stella_parrallel set')
    ap.add_argument('--out', default=DEFAULT_OUT)
    args = ap.parse_args()

    models = dict(DEFAULT_MODELS)
    if args.models:
        models = {}
        for spec_ in args.models:
            label, _, d = spec_.partition('=')
            models[label] = {'dir': d, 'color': None, 'ls': '-'}

    palette = ['black', 'blue', 'red', 'green', 'purple', 'orange']
    data = {}
    for i, (label, info) in enumerate(models.items()):
        model_dir = PROJECT_DIR / info['dir']
        got = collect(model_dir)
        info['color'] = info['color'] or palette[i % len(palette)]
        if got is None:
            print(f"  {label}: no usable log in {model_dir}")
            continue
        stamps, epochs, job_id, job_elapsed = got
        # plot and target in days-since-sim-start so the axis is positive
        if epochs is not None and len(epochs):
            epochs = epochs - epochs[0]
        t_target, ep_reached, ok = time_to_target(stamps, epochs, args.target)
        n = min(len(stamps), len(epochs)) if epochs is not None else 0
        data[label] = {
            'stamps': stamps, 'epochs': epochs, 'n': n,
            'job_id': job_id, 'job_elapsed': job_elapsed,
            't_target': t_target, 'epoch_reached': ep_reached, 'ok': ok,
        }
        status = f"{t_target/60:.1f} min to {ep_reached:.2f}d" if ok else \
            f"INCOMPLETE: {ep_reached:.2f}d after {t_target/60:.1f} min" \
            if t_target else 'no epochs yet'
        print(f"  {label}: job={job_id} stamps={len(stamps)} "
              f"epochs={0 if epochs is None else len(epochs)} -> {status}")

    if not data:
        print("ERROR: nothing to plot (jobs may not have started)")
        sys.exit(1)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

    # Panel 1: wall-clock vs simulation day. Index-aligned pairs.
    for label, info in models.items():
        d = data.get(label)
        if d is None or d['epochs'] is None or d['n'] == 0:
            continue
        wall_min = (d['stamps'][:d['n']] - d['stamps'][0]) / 60.0
        ax1.plot(wall_min, d['epochs'][:d['n']],
                 color=info['color'], linestyle=info['ls'],
                 linewidth=1.2, label=label)
    ax1.set_xlabel('Elapsed wall time since solver start (min)')
    ax1.set_ylabel('Simulation time since start (days)')
    ax1.set_title('STELLA progress: wall clock vs simulation day')
    ax1.legend(fontsize=8)
    ax1.grid(True, alpha=0.3)

    # Panel 2: solver wall time to reach target day.
    labels = [l for l in models if l in data and data[l]['t_target'] is not None]
    xs = np.arange(len(labels))
    solver = [data[l]['t_target'] / 60.0 for l in labels]
    ax2.bar(xs, solver, 0.6,
            color=[models[l]['color'] for l in labels])
    for i, l in enumerate(labels):
        ax2.text(xs[i], solver[i], f"{solver[i]:.0f}m",
                 ha='center', va='bottom', fontsize=9)
    ax2.set_xticks(xs)
    ax2.set_xticklabels(labels, fontsize=8, rotation=15, ha='right')
    ax2.set_ylabel('Minutes')
    ax2.grid(True, alpha=0.3, axis='y')

    fig.tight_layout()
    out = PROJECT_DIR / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"\nSaved: {out}")

    print(f"\n{'Model':<18} {'Job':>9} {'Solver>10d':>11} {'JobElapsed':>11} "
          f"{'Epoch':>8} {'Reached':>8}")
    for l in labels:
        d = data[l]
        print(f"{l:<18} {str(d['job_id']):>9} "
              f"{d['t_target']/60:>9.1f}m "
              f"{(d['job_elapsed']/60) if d['job_elapsed'] else float('nan'):>9.1f}m "
              f"{d['epoch_reached']:>7.2f}d {str(d['ok']):>8}")


if __name__ == '__main__':
    main()
