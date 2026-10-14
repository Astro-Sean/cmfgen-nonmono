#!/usr/bin/env python3
"""Compare SN_B X models with and without the SN_A ejecta.

Top panel: emergent bolometric lightcurve -- the .swd Flum column
(col 11, luminosity in erg/s) evaluated at the outermost zone.
Bottom panel: local luminosity L(r) (the .swd Flum column, erg/s)
at matched epochs.
Solid = with SN_A ejecta (X); dashed = without (X_wind).

X_wind is shown over its computed range only; the run stalled at
~2.9 d (timestep collapse near shock emergence) so its curve is
necessarily partial.
"""

import argparse
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import AutoMinorLocator
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from extract_lightcurve_from_swd import read_swd

GRID_DIR = "/ptmp/sbrennan/SN2024pba/SN_B_X_grid"

# (model, linestyle, label, color)
WITH = [
    ("X", "-", r"X: with SN_A ejecta", "k"),
]
WITHOUT = [
    ("X_wind", "--", r"X_wind: no SN_A, CSM to 6e15", "tab:blue"),
    ("X_windext", "--", r"X_windext: no SN_A, wind to 3.2e16",
     "tab:red"),
]

plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["DejaVu Serif"],
        "mathtext.fontset": "dejavuserif",
        "font.size": 8,
        "axes.labelsize": 8,
        "legend.fontsize": 6.5,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "axes.linewidth": 0.7,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


def find_swd(mdir):
    name = os.path.basename(mdir)
    for c in (
        os.path.join(mdir, "stella", f"{name}.swd"),
        os.path.join(mdir, "stella", "strad", f"{name}.swd"),
    ):
        if os.path.exists(c):
            return c
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grid-dir", default=GRID_DIR)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    out = args.out or os.path.join(args.grid_dir, "x_with_without_sna.png")

    fig, (ax_lc, ax_lum) = plt.subplots(
        2, 1, figsize=(3.54, 4.6), height_ratios=[1.2, 1.0]
    )

    # --- top: emergent luminosity L_out(t) from Flum ---
    epoch_cache = {}
    for name, ls, lab, color in WITH + WITHOUT:
        swd = find_swd(os.path.join(args.grid_dir, name))
        if swd is None:
            print(f"  {name}: no .swd, skipped")
            continue
        eps = read_swd(swd)
        epoch_cache[name] = eps
        ts = np.array(sorted(eps.keys()))
        l_out = np.array([eps[t]["flux"][0] for t in ts])
        pos = l_out > 0
        ax_lc.plot(
            ts[pos],
            l_out[pos],
            color=color,
            ls=ls,
            lw=1.1,
            zorder=10,
        )
        ax_lc.plot(
            ts[pos][-1],
            l_out[pos][-1],
            "o",
            ms=3,
            color=color,
            zorder=11,
        )
        if name == "X_wind" and ts[-1] < 10:
            ax_lc.annotate(
                f"stalled at\n{ts[-1]:.1f} d",
                xy=(ts[pos][-1], l_out[pos][-1]),
                xytext=(12, 3e39),
                fontsize=5.8,
                color=color,
                arrowprops=dict(arrowstyle="-", lw=0.6, color=color),
            )
        print(f"  {name}: {len(ts)} epochs, "
              f"{ts[0]:.1f}-{ts[-1]:.1f} d")

    ax_lc.set_yscale("log")
    ax_lc.set_ylabel(r"$L_{\rm bol}$ (erg s$^{-1}$)")
    ax_lc.set_xlabel("days since shell launch")

    # --- bottom: local luminosity L(r) at matched epochs ---
    # .swd 'flux' field is STELLA's Flum: luminosity in erg/s, not a
    # flux per area -- plotted directly against r.
    # pick epochs present in the shortest-coverage run, preferring
    # runs that actually reached t>0 so matched epochs are post-launch
    post = [n for n in epoch_cache if max(epoch_cache[n]) > 0]
    ref = min(post or epoch_cache, key=lambda n: max(epoch_cache[n]))
    ref_times = np.array(sorted(epoch_cache[ref].keys()))
    if len(ref_times) < 5:
        ref_times = np.array(sorted(
            epoch_cache[min(epoch_cache, key=lambda n: len(epoch_cache[n]))].keys()))
    # sample evenly in time, not epoch index -- a run whose timestep
    # collapses packs most .swd records into its last fraction of a day
    span = ref_times[-1] - ref_times[0]
    want = np.linspace(
        ref_times[0] + 0.05 * span,
        ref_times[-1] - 0.02 * span,
        5,
    )
    epoch_ts = [ref_times[np.argmin(np.abs(ref_times - w))] for w in want]
    cmap = plt.get_cmap("viridis")
    colors = cmap(np.linspace(0.05, 0.9, len(epoch_ts)))

    for t, c in zip(epoch_ts, colors):
        for name, ls, lab, col in WITH[:1] + WITHOUT:
            if name not in epoch_cache:
                continue
            eps = epoch_cache[name]
            ts = np.array(sorted(eps.keys()))
            i = np.argmin(np.abs(ts - t))
            if abs(ts[i] - t) > 1.0:
                continue
            d = eps[ts[i]]
            ax_lum.plot(
                d["r"],
                d["flux"],
                color=c,
                ls=ls,
                lw=1.0,
                zorder=10,
            )
        print(f"  epoch t={t:.2f} d plotted")

    ax_lum.set_xscale("log")
    ax_lum.set_yscale("log")
    ax_lum.set_xlabel("radius (cm)")
    ax_lum.set_ylabel(r"$L(r)$ (erg s$^{-1}$)")

    # epoch colour legend
    ep_handles = [
        Line2D([], [], color=c, lw=1.4, label=f"{t:.2f} d")
        for t, c in zip(epoch_ts, colors)
    ]
    style_handles = [
        Line2D([], [], color="0.3", ls="-", lw=1.2, label="with SN_A"),
        Line2D([], [], color="0.3", ls="--", lw=1.2, label="without SN_A"),
    ]
    leg1 = ax_lum.legend(
        handles=ep_handles, fontsize=5.8, loc="lower left",
        title="epoch", title_fontsize=6, frameon=False,
        handlelength=1.5, labelspacing=0.3,
    )
    ax_lum.add_artist(leg1)
    ax_lum.legend(handles=style_handles, fontsize=6, loc="upper right",
                  frameon=False)

    handles = [
        Line2D([], [], color=col, ls=ls, lw=1.1, label=lab)
        for _, ls, lab, col in WITH + WITHOUT
    ]
    ax_lc.legend(handles=handles, fontsize=6.5, loc="lower left",
                 frameon=False)

    for ax in (ax_lc, ax_lum):
        ax.grid(True, alpha=0.3)
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.grid(True, which="minor", alpha=0.12)

    fig.align_ylabels((ax_lc, ax_lum))
    plt.tight_layout()
    fig.savefig(out, dpi=300, bbox_inches="tight")
    pdf = os.path.splitext(out)[0] + ".pdf"
    fig.savefig(pdf, bbox_inches="tight")
    print(f"Saved: {out}")
    print(f"Saved: {pdf}")


if __name__ == "__main__":
    main()
