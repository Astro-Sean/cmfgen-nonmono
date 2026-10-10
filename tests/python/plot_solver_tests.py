#!/usr/bin/env python3
"""Plot JFNK solver test convergence histories.

Reads the solver_hist_*.csv files that test_jfnk_solver writes (one
residual + dtau per accepted Newton iterate, recorded in the sub()
sync callback) and writes PNGs under tests/plots/solver/.

Run from tests/f90 (where the CSVs land) or pass the directory:

    python3 plot_solver_tests.py [csv_dir]
"""

import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

CASES = {
    "solver_hist_near.csv": "near start, PTC off",
    "solver_hist_far_ptc.csv": "far start, PTC on",
}


def load(path):
    its, res, dtau, info = [], [], [], None
    for row in csv.reader(open(path)):
        if not row:
            continue
        if row[0].startswith("#"):
            if "solver_info" in row[0]:
                info = int(row[0].split("=")[1])
            continue
        its.append(int(row[0]))
        res.append(float(row[1]))
        dtau.append(float(row[2]))
    return its, res, dtau, info


def main():
    csv_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.cwd()
    if not csv_dir.is_dir():
        csv_dir = Path(__file__).resolve().parents[1] / "f90"
    outdir = Path(__file__).resolve().parents[1] / "plots" / "solver"
    outdir.mkdir(parents=True, exist_ok=True)

    series = []
    for fname, label in CASES.items():
        p = csv_dir / fname
        if not p.exists():
            print(f"missing {p} - run test_jfnk_solver first")
            continue
        its, res, dtau, info = load(p)
        series.append((label, its, res, dtau, info))

    if not series:
        sys.exit("no solver history CSVs found")

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(8, 6), sharex=True,
        gridspec_kw={"height_ratios": [3, 1]})
    for label, its, res, dtau, info in series:
        tag = f"{label} (info={info})"
        ax1.semilogy(its, res, "o-", ms=4, label=tag)
        ax2.semilogy(its, dtau, "s-", ms=4, label=tag)
    ax1.set_ylabel(r"weighted $|F|$")
    ax1.set_title("newton_hook_v1 residual per accepted iterate "
                  "(production solver, synthetic problems)")
    ax1.grid(alpha=0.3)
    ax1.legend()
    ax2.set_ylabel(r"$\Delta\tau$ (PTC)")
    ax2.set_xlabel("accepted Newton iterate")
    ax2.grid(alpha=0.3)
    fig.tight_layout()
    out = outdir / "newton_convergence.png"
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
