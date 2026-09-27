#!/usr/bin/env python3
"""Compare STELLA lightcurves of two models with an optional time shift.

Reads the ``<model>_lightcurve.csv`` (or ``*_lightcurve_complete.csv``) files
in each model's ``output/`` directory and plots bolometric luminosity plus
absolute magnitudes.  The second model's time axis can be shifted to align
features; by default the shift is chosen so the two bolometric peaks coincide.

Usage:
    python scripts/plot_two_model_lightcurves.py MODEL_A MODEL_B
    python scripts/plot_two_model_lightcurves.py SN_B_v2 SN_B_v2_94W --shift 28
"""

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BANDS = ["U", "B", "V", "R", "I"]
BAND_COLORS = {
    "U": "#7B00D4",
    "B": "#0000FF",
    "V": "#00AA00",
    "R": "#FF4400",
    "I": "#AA0000",
}


def find_csv(model_dir):
    """Return the lightcurve CSV in model_dir/output/ with the most bands."""
    out = Path(model_dir) / "output"
    candidates = sorted(out.glob("*_lightcurve*.csv"))
    if not candidates:
        raise FileNotFoundError(f"no lightcurve CSV under {out}")

    def n_bands(path):
        with open(path) as f:
            header = f.readline().split(",")
        return sum(b in header for b in BANDS)

    return max(candidates, key=n_bands)


def read_csv(path):
    """Read a lightcurve CSV into column arrays."""
    cols = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            for key, val in row.items():
                cols.setdefault(key, []).append(float(val) if val else np.nan)
    return {k: np.asarray(v) for k, v in cols.items()}


def peak_time(time, lum):
    """Time of the maximum bolometric luminosity."""
    return float(time[np.nanargmax(lum)])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("model_a", type=Path, help="reference model directory")
    ap.add_argument("model_b", type=Path, help="model to shift onto A")
    ap.add_argument(
        "--shift",
        type=float,
        default=None,
        help="days added to model B times (default: align bolometric peaks)",
    )
    ap.add_argument("--output", "-o", type=Path, default=None)
    args = ap.parse_args()

    csv_a = find_csv(args.model_a)
    csv_b = find_csv(args.model_b)
    a = read_csv(csv_a)
    b = read_csv(csv_b)

    t_a = a["time"]
    t_b = b["time"]
    shift = args.shift
    if shift is None:
        shift = peak_time(t_a, a["bol_ergs"]) - peak_time(t_b, b["bol_ergs"])
    t_b = t_b + shift

    label_a = args.model_a.name
    label_b = f"{args.model_b.name} (shift {shift:+.1f} d)"

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(9, 8), sharex=True, gridspec_kw={"height_ratios": [1, 1.4]}
    )

    ax1.plot(t_a, a["bol_ergs"], color="k", lw=1.6, label=label_a)
    ax1.plot(t_b, b["bol_ergs"], color="tab:red", lw=1.6, label=label_b)
    for t, l, c in (
        (peak_time(t_a, a["bol_ergs"]), a["bol_ergs"].max(), "k"),
        (peak_time(t_b, b["bol_ergs"]), b["bol_ergs"].max(), "tab:red"),
    ):
        ax1.plot(t, l, marker="v", color=c, ms=7, zorder=5)
    ax1.set_yscale("log")
    ax1.set_ylabel(r"$L_{\rm bol}$ [erg s$^{-1}$]")
    ax1.legend(loc="upper right", fontsize=9)
    ax1.grid(alpha=0.3)
    ax1.set_title(
        f"{args.model_a.name} vs {args.model_b.name} lightcurves", fontsize=11
    )

    for band in BANDS:
        if band in a:
            ax2.plot(
                t_a,
                a[band],
                color=BAND_COLORS[band],
                lw=1.4,
                label=f"{label_a} {band}",
            )
        if band in b:
            ax2.plot(
                t_b,
                b[band],
                color=BAND_COLORS[band],
                lw=1.4,
                ls="--",
                label=f"{args.model_b.name} {band}",
            )
    ax2.invert_yaxis()
    ax2.set_xlabel("observer-frame time [days]")
    ax2.set_ylabel("absolute magnitude")
    ax2.legend(ncol=2, fontsize=8, loc="upper right")
    ax2.grid(alpha=0.3)

    fig.tight_layout()
    out = args.output or (
        args.model_b / "plots" / f"{args.model_b.name}_vs_{args.model_a.name}_lc.png"
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150)
    print(f"peak {label_a}: t={peak_time(t_a, a['bol_ergs']):.2f} d")
    print(f"peak {args.model_b.name}: t={peak_time(t_b, b['bol_ergs']):.2f} d")
    print(f"applied shift to {args.model_b.name}: {shift:+.2f} d")
    print(f"saved {out}")


if __name__ == "__main__":
    main()
