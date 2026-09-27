#!/usr/bin/env python3
"""Compare all completed SN_B model lightcurves against SN 2024pba obs.

Each model's o-band curve is shifted in time only, so its o-peak lands on
the observed ATLAS o peak MJD; magnitudes stay absolute so peak-level and
first-peak differences remain visible.  r-band shown dashed for the
first-peak/plateau comparison against the ~-15.5 mag target.
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_sn_a_pba_grid_vs_obs import BASE, DM, load_model
from plot_sn_b_grid_vs_obs import A_O, A_R, PEAK_MJD, load_atlas
from plot_sn_b_par_vs_obs import SNA_OBS, SNB_OBS, load_obs_det

GRID = os.path.join(BASE, "SN_B_double_peak_grid")
V2 = os.path.join(BASE, "SN_B_v2_T6k_grid")
BEST = "SB_M10_ts8_W3e-4"

OUT = os.path.join(GRID, "plots", "completed_models_vs_obs.png")


def discover_completed():
    """All grid dirs with a usable .tt, plus the v2 reference pair.

    Returns list of (mdir, label, color, lw, ls, alpha, zorder).
    """
    import glob
    import itertools
    models = []
    cyc = itertools.cycle(plt.cm.tab20.colors)
    for mdir in sorted(glob.glob(os.path.join(GRID, "*"))):
        name = os.path.basename(mdir)
        if not os.path.isdir(mdir) or name.startswith(("plots", "recon")):
            continue
        df, src = load_model(mdir)
        if df is None or "o" not in df or df["o"].notna().sum() < 5:
            continue
        if name == BEST:
            models.append((mdir, name + " (best)", "tab:blue", 2.6, "-", 1.0, 8))
        else:
            models.append((mdir, name, next(cyc), 1.0, "-", 0.55, 3))
    for ref, col in ((os.path.join(V2, "T6k_M7.0_V10000"), "k"),
                     (os.path.join(V2, "T6k_M10.0_V8000"), "0.55")):
        df, src = load_model(ref)
        if df is not None:
            models.append((ref, "v2 " + os.path.basename(ref) + " (ref)",
                           col, 1.4, ":", 0.9, 4))
    return models


def main():
    os.makedirs(os.path.dirname(OUT), exist_ok=True)

    snb_r, _ = load_obs_det(SNB_OBS, "r")
    snb_g, _ = load_obs_det(SNB_OBS, "g")
    atlas = load_atlas(
        os.path.join(BASE, "observational_data",
                     "SN2024pba_atclean_cleaned.csv"), band="o")
    o_det = atlas["det"].copy()
    if len(o_det):
        o_det["MJD"] = o_det["phase"] + PEAK_MJD
    ob_peak_mjd = o_det.loc[o_det["Mabs"].idxmin(), "MJD"]

    fig, ax = plt.subplots(figsize=(13, 7))

    for mdir, lbl, col, lw, ls, alpha, zo in discover_completed():
        df, src = load_model(mdir)
        if df is None:
            print(f"skip {lbl}: no .tt")
            continue
        band = "o" if "o" in df else "r"
        pk = df[band].idxmin()
        t0 = ob_peak_mjd - df.loc[pk, "time"]
        mjd = t0 + df["time"].values
        ax.plot(mjd, df[band], color=col, lw=lw, ls=ls, alpha=alpha,
                zorder=zo, label=f"{lbl} o (pk {df.loc[pk, band]:.2f})")
        if "r" in df:
            ax.plot(mjd, df["r"], color=col, lw=lw * 0.55, ls="--",
                    alpha=alpha * 0.6, zorder=zo)

    if len(o_det):
        ax.errorbar(o_det["MJD"], o_det["Mabs"], yerr=o_det["Mabs_err"],
                    fmt="^", ms=5, capsize=2, lw=0.8, color="darkorange",
                    label="obs o (ATLAS)", zorder=10)
    if len(snb_r):
        ax.errorbar(snb_r["MJD"], snb_r["Mabs"], yerr=snb_r["Mabs_err"],
                    fmt="o", ms=4, capsize=2, lw=0.8, color="tab:red",
                    mfc="none", label="obs r (SN_B)", zorder=10)
    if len(snb_g):
        ax.errorbar(snb_g["MJD"], snb_g["Mabs"], yerr=snb_g["Mabs_err"],
                    fmt="s", ms=4, capsize=2, lw=0.8, color="tab:green",
                    mfc="none", label="obs g (SN_B)", zorder=10)

    ax.axhline(-15.5, color="tab:red", ls=":", lw=1, alpha=0.6)
    ax.annotate("r plateau target -15.5", xy=(60812, -15.6), fontsize=8,
                color="tab:red")
    ax.axhline(-18.54, color="darkorange", ls=":", lw=1, alpha=0.6)
    ax.annotate("o peak target -18.54", xy=(60812, -18.7), fontsize=8,
                color="darkorange")
    ax.axvline(ob_peak_mjd, color="darkorange", ls=":", lw=1, alpha=0.5)

    ax.set_ylim(-8, -20)
    ax.set_xlim(60785, 60960)
    ax.set_xlabel("MJD (model o-peak aligned to obs peak)")
    ax.set_ylabel("Absolute magnitude (AB)")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, ncol=2, loc="lower left")
    ax.set_title("Completed SN_B models vs SN 2024pba "
                 "(o solid, r dashed; time-aligned at o peak)")
    plt.tight_layout()
    plt.savefig(OUT, dpi=150, bbox_inches="tight")
    print("saved", OUT)


if __name__ == "__main__":
    main()
