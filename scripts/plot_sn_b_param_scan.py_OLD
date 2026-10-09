#!/usr/bin/env python3
"""SN_B double-peak parameter scan on model time.

Rows = photometric bands (SDSS g, r, i + ATLAS o where the .tt provides
them), columns = parameter axes (shell mass, shell velocity, wind mdot,
t_start). Curves are plotted vs days since model start so the internal
double peak - shell/wind collision SBO, dip, ejecta-heating peak - is
visible directly. No Ni56 in any model; the late peak is interaction
heating of the outer SN_A ejecta.

Old 14-col .tt files carry Johnson UBVRI only; those models are drawn
with dotted lines and labelled (proxy).
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_sn_a_pba_grid_vs_obs import BASE, DM, load_model

OUT = os.path.join(BASE, "SN_B_double_peak_grid",
                   "sn_b_param_scan_time.png")

DP = os.path.join(BASE, "SN_B_double_peak_grid")
V2 = os.path.join(BASE, "SN_B_v2_T6k_grid")

# (group title, [(model_dir, label, color, linestyle)])
AXES = [
    ("shell mass  (V=10k, W=1e-3)", [
        (f"{DP}/SN_B_T6k_M5.0_V10000", "M=5", "tab:blue", "-"),
        (f"{DP}/SN_B_T6k_M7.0_V10000", "M=7", "k", "-"),
        (f"{DP}/SN_B_T6k_M10.0_V10000", "M=10", "tab:red", "-"),
        (f"{DP}/SN_B_T6k_M12.0_V10000", "M=12", "tab:purple", "-"),
    ]),
    ("shell velocity  (M=7, W=1e-3)", [
        (f"{DP}/SN_B_T6k_M7.0_V8000", "v=8000", "tab:cyan", "--"),
        (f"{DP}/SN_B_T6k_M7.0_V10000", "v=10000", "k", "-"),
        (f"{DP}/SN_B_T6k_M7.0_V12000", "v=12000", "tab:orange", "-"),
    ]),
    ("wind mdot  (M=7, V=10k)", [
        (f"{DP}/SN_B_T6k_M7.0_V10000_W1e-4", "W=1e-4", "tab:gray", "-"),
        (f"{DP}/SN_B_T6k_M7.0_V10000_W3e-4", "W=3e-4", "tab:blue", "-"),
        (f"{DP}/SN_B_T6k_M7.0_V10000", "W=1e-3", "k", "-"),
        (f"{DP}/SN_B_T6k_M7.0_V10000_W2e-3", "W=2e-3", "tab:green", "-"),
        (f"{V2}/X_M7p0_V1e4_W3e-3", "W=3e-3", "tab:olive", "-"),
        (f"{V2}/X_M7p0_V1e4_W5e-3", "W=5e-3", "tab:brown", "-"),
        (f"{V2}/X_M7p0_V1e4_W1e-2", "W=1e-2", "tab:red", "-"),
        (f"{V2}/X_M7p0_V1e4_W3e-2", "W=3e-2", "tab:pink", "-"),
    ]),
    ("t_start  (M=7, V=10k, W=1e-3)", [
        (f"{DP}/SN_B_T6k_M7.0_V10000", "ts=5", "k", "-"),
        (f"{DP}/SN_B_T6k_M7.0_V10000_ts20", "ts=20", "tab:orange", "-."),
        (f"{DP}/SN_B_T6k_M7.0_V10000_ts60", "ts=60", "tab:brown", ":"),
    ]),
]

BANDS = [("g", "SDSS g"), ("r", "SDSS r"), ("i", "SDSS i"),
         ("o", "ATLAS o")]


def main():
    fig, axes = plt.subplots(len(BANDS), len(AXES), figsize=(19, 12),
                             sharex=True)
    for j, (axtitle, models) in enumerate(AXES):
        for i, (band, btitle) in enumerate(BANDS):
            ax = axes[i, j]
            any_curve = False
            for path, label, color, ls in models:
                if not os.path.isdir(path):
                    continue
                df, src = load_model(path)
                if df is None or band not in df or df.empty:
                    continue
                proxy = src != "tt_sdss"
                ax.plot(df["time"], df[band], color=color,
                        ls=":" if proxy else ls, lw=1.3,
                        label=label + (" (proxy)" if proxy else ""))
                any_curve = True
            ax.set_yscale("linear")
            ax.invert_yaxis()
            ax.grid(alpha=0.3)
            if i == 0:
                ax.set_title(axtitle, fontsize=10)
            if j == 0:
                ax.set_ylabel(f"M_{band} (AB)", fontsize=11)
                ax.annotate(btitle, xy=(0.02, 0.04),
                            xycoords="axes fraction", fontsize=9,
                            color="0.4")
            if i == len(BANDS) - 1:
                ax.set_xlabel("days since model start", fontsize=10)
            if any_curve:
                ax.legend(fontsize=7, loc="best")
            else:
                ax.text(0.5, 0.5, "awaiting\nruns", ha="center",
                        va="center", transform=ax.transAxes,
                        fontsize=9, color="0.6")
    fig.suptitle(
        "SN_B double-peak parameter scan - model time axis\n"
        "pk1 = shell/wind collision SBO, pk2 = heated outer SN_A ejecta "
        "(no Ni56)", fontsize=12)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(OUT, dpi=150, bbox_inches="tight")
    print("saved", OUT)


if __name__ == "__main__":
    main()
