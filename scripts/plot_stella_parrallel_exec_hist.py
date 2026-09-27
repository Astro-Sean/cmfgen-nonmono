#!/usr/bin/env python3
"""
Histogram of execution times for the Stella_parrallel benchmark runs.

Two panels:
  1. Distribution of per-epoch solver wall time (seconds between
     successive `time=` REAL_CLOCK stamps in the STELLA log). This shows
     how much each epoch costs under Y12M vs PARDISO.
  2. Total execution time per run: solver wall time (first to last stamp)
     and sacct job elapsed for context.

Usage:
  python scripts/plot_stella_parrallel_exec_hist.py
  python scripts/plot_stella_parrallel_exec_hist.py \
      --models "Y12M=Stella_parrallel/spar_base" \
      --out Stella_parrallel/plots/exec_time_hist.png
"""

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

PROJECT_DIR = Path(__file__).resolve().parent.parent

sys.path.insert(0, str(PROJECT_DIR / 'scripts'))
from plot_stella_parrallel_timing import (
    DEFAULT_MODELS,
    find_log,
    job_elapsed_sacct,
    parse_time_stamps,
)

DEFAULT_OUT = 'Stella_parrallel/plots/exec_time_hist.png'


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--models', nargs='*', default=None,
                    help='label=dir pairs; defaults to the Stella_parrallel set')
    ap.add_argument('--out', default=DEFAULT_OUT,
                    help='output path for the plot, relative to project dir')
    ap.add_argument('--bins', type=int, default=25,
                    help='number of histogram bins (default 25)')
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
        info['color'] = info.get('color') or palette[i % len(palette)]
        out_file, job_id = find_log(PROJECT_DIR / info['dir'])
        if out_file is None:
            print(f"  {label}: no usable log")
            continue
        stamps = parse_time_stamps(out_file)
        if len(stamps) < 2:
            print(f"  {label}: only {len(stamps)} stamps")
            continue
        deltas = np.diff(stamps)
        # deltas include midnight-unwrap padding; drop non-positive values
        deltas = deltas[deltas > 0]
        data[label] = {
            'deltas': deltas,
            'total_s': float(stamps[-1] - stamps[0]),
            'job_s': job_elapsed_sacct(job_id),
            'job_id': job_id,
        }
        print(f"  {label}: job={job_id} epochs={len(deltas)} "
              f"median={np.median(deltas):.1f}s mean={np.mean(deltas):.1f}s "
              f"total={stamps[-1] - stamps[0]:.0f}s")

    if not data:
        print("ERROR: nothing to plot")
        sys.exit(1)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))

    all_d = np.concatenate([d['deltas'] for d in data.values()])
    hi = np.percentile(all_d, 99)
    bins = np.linspace(0, hi, args.bins + 1)
    for label, info in models.items():
        d = data.get(label)
        if d is None:
            continue
        ax1.hist(d['deltas'], bins=bins, alpha=0.45, color=info['color'],
                 label=f"{label} (med {np.median(d['deltas']):.0f}s)",
                 edgecolor=info['color'])
    ax1.set_xlabel('Wall time per epoch (s)')
    ax1.set_ylabel('Count')
    ax1.set_title('Distribution of per-epoch solver time')
    ax1.legend(fontsize=8)
    ax1.grid(True, alpha=0.3, axis='y')

    labels = [l for l in models if l in data]
    xs = np.arange(len(labels))
    solver = [data[l]['total_s'] / 60.0 for l in labels]
    job = [data[l]['job_s'] / 60.0 if data[l]['job_s'] else np.nan
           for l in labels]
    w = 0.38
    ax2.bar(xs - w / 2, solver, w, label='solver wall time',
            color=[models[l]['color'] for l in labels])
    ax2.bar(xs + w / 2, job, w, label='job elapsed (sacct)',
            color='lightgray', edgecolor='gray')
    for i in range(len(labels)):
        ax2.text(xs[i] - w / 2, solver[i], f"{solver[i]:.0f}m",
                 ha='center', va='bottom', fontsize=8)
        if not np.isnan(job[i]):
            ax2.text(xs[i] + w / 2, job[i], f"{job[i]:.0f}m",
                     ha='center', va='bottom', fontsize=8, color='gray')
    ax2.set_xticks(xs)
    ax2.set_xticklabels(labels, fontsize=8, rotation=15, ha='right')
    ax2.set_ylabel('Minutes')
    ax2.set_title('Total execution time')
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.3, axis='y')

    fig.tight_layout()
    out = PROJECT_DIR / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"\nSaved: {out}")


if __name__ == '__main__':
    main()
