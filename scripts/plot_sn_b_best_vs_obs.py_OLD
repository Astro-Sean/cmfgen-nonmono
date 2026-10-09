#!/usr/bin/env python3
"""Zoomed SN_B-window plot: best-fit double-peak model vs SN 2024pba obs.

Best model is the top RMS scorer among models with true double-peak
morphology (see score_sn_b_grid_rms.py): SB_M10_ts8_W3e-4.
Context curves: SB_M8_ts8_W3e-4 (runner-up double peak),
SN_B_T6k_M7.0_V10000_ts17 (earlier morphology pick, too bright),
v2 T6k_M7.0_V10000 (serial-generation reference).

Models are time-shifted only (o-peak anchored to observed ATLAS o peak);
magnitudes stay absolute.
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_sn_a_pba_grid_vs_obs import BASE, DM, load_model
from plot_sn_b_grid_vs_obs import A_O, PEAK_MJD, load_atlas
from plot_sn_b_par_vs_obs import SNB_OBS, load_obs_det

GRID = os.path.join(BASE, "SN_B_double_peak_grid")
V2 = os.path.join(BASE, "SN_B_v2_T6k_grid")

MODELS = [
    (os.path.join(GRID, "SB_M10_ts8_W3e-4"),
     "SB_M10_ts8_W3e-4 (best)", "tab:blue", 2.4, "-"),
    (os.path.join(GRID, "SB_M8_ts8_W3e-4"),
     "SB_M8_ts8_W3e-4", "tab:cyan", 1.4, "-"),
    (os.path.join(GRID, "SN_B_T6k_M7.0_V10000_ts17"),
     "M7 ts17 (too bright)", "tab:purple", 1.4, "-"),
    (os.path.join(GRID, "SB_M10_ts16_W3e-4"),
     "M10 ts16 W3e-4 (single)", "0.6", 1.2, "-"),
    (os.path.join(GRID, "SB_M10_ts8_nr12"),
     "M10 ts8 nr12 (n_rho=12)", "tab:brown", 1.6, "-"),
    (os.path.join(V2, "T6k_M7.0_V10000"),
     "v2 M7 ref", "k", 1.2, ":"),
]

OUT = os.path.join(GRID, "plots", "sn_b_best_model_vs_obs.png")


def main():
    os.makedirs(os.path.dirname(OUT), exist_ok=True)

    atlas = load_atlas(
        os.path.join(BASE, "observational_data",
                     "SN2024pba_atclean_cleaned.csv"), band="o")
    o_det = atlas["det"].copy()
    o_det["MJD"] = o_det["phase"] + PEAK_MJD
    ob_peak_mjd = o_det.loc[o_det["Mabs"].idxmin(), "MJD"]
    ob_peak_mag = o_det["Mabs"].min()
    r_det, _ = load_obs_det(SNB_OBS, "r")
    g_det, _ = load_obs_det(SNB_OBS, "g")

    fig, ax = plt.subplots(figsize=(12, 7))

    for mdir, lbl, col, lw, ls in MODELS:
        df, src = load_model(mdir)
        if df is None:
            print(f"skip {lbl}: no .tt")
            continue
        pk = df["o"].idxmin() if "o" in df else df["r"].idxmin()
        band = "o" if "o" in df else "r"
        t0 = ob_peak_mjd - df.loc[pk, "time"]
        mjd = t0 + df["time"].values
        ax.plot(mjd, df[band], color=col, lw=lw, ls=ls,
                label=f"{lbl} o (pk {df.loc[pk, band]:.2f})", zorder=5)
        if "r" in df:
            ax.plot(mjd, df["r"], color=col, lw=lw * 0.5, ls="--",
                    alpha=0.55, zorder=4)

    ax.errorbar(o_det["MJD"], o_det["Mabs"], yerr=o_det["Mabs_err"],
                fmt="^", ms=5, capsize=2, lw=0.8, color="darkorange",
                label="obs o (ATLAS)", zorder=10)
    ax.errorbar(r_det["MJD"], r_det["Mabs"], yerr=r_det["Mabs_err"],
                fmt="o", ms=4, capsize=2, lw=0.8, color="tab:red",
                mfc="none", label="obs r (SN_B)", zorder=10)
    ax.errorbar(g_det["MJD"], g_det["Mabs"], yerr=g_det["Mabs_err"],
                fmt="s", ms=4, capsize=2, lw=0.8, color="tab:green",
                mfc="none", label="obs g (SN_B)", zorder=10)

    ax.axhline(-15.6, color="tab:red", ls=":", lw=1, alpha=0.6)
    ax.annotate("first-peak target ~-15.6", xy=(60795, -15.5), fontsize=8,
                color="tab:red")
    ax.axhline(ob_peak_mag, color="darkorange", ls=":", lw=1, alpha=0.6)
    ax.annotate(f"o peak {ob_peak_mag:.2f}", xy=(60795, ob_peak_mag - 0.15),
                fontsize=8, color="darkorange")
    ax.axvline(ob_peak_mjd, color="darkorange", ls=":", lw=1, alpha=0.5)

    ax.set_ylim(-13.5, -20)
    ax.set_xlim(60790, 60960)
    ax.set_xlabel("MJD (model o-peak aligned to obs o peak)")
    ax.set_ylabel("Absolute magnitude (AB)")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, ncol=2, loc="lower left")
    ax.set_title(
        "SN_B best-fit: SB_M10_ts8_W3e-4 "
        "(M=10 Msun, v=10000 km/s, ts=8 d, mdot=3e-4)\n"
        "o-band solid, r-band dashed; model o-peaks time-aligned to obs")
    plt.tight_layout()
    plt.savefig(OUT, dpi=150, bbox_inches="tight")
    print("saved", OUT)


if __name__ == "__main__":
    main()
