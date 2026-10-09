#!/usr/bin/env python3
"""r-band lightcurves of the SN_B_cocoon models vs SN 2024pba SN_B data.

For every model dir under SN_B_cocoon with a STELLA .tt lightcurve,
the SDSS r-band curve (Mr col of the 24-col .tt; Johnson MR proxy for
old 14-col files) is plotted against the observed r-band photometry.
Each model is time-anchored with the campaign convention: its own
ATLAS o-band (Mo) peak is placed on the observed ATLAS o peak, and the
x-axis is days since that peak.

Best matches are ranked by the inverse-variance weighted RMS between
the model r curve and the observed r detections inside the model's
time coverage (same band_rms metric as score_sn_b_grid_rms.py).  The
combined o/r/g score is also printed for cross-checking against the
campaign table in SN_B_cocoon/MODEL_INDEX.md.

Left panel: ejecta-bearing models.  Right panel: paired _windctl
controls (shell+wind only) - they should sit far below the data if
the ejecta interaction supplies the flux.

Usage:
    python scripts/plot_sn_b_cocoon_rband.py [--top N] [--out FIG.png]
"""
import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_sn_a_pba_grid_vs_obs import BASE, load_model
from plot_sn_b_grid_vs_obs import PEAK_MJD, load_atlas
from plot_sn_b_par_vs_obs import SNB_OBS, load_obs_det

GRID = os.path.join(BASE, "SN_B_cocoon")
ATLAS = os.path.join(BASE, "observational_data",
                     "SN2024pba_atclean_cleaned.csv")
OUT_FILE = os.path.join(GRID, "sn_b_cocoon_rband_vs_obs.png")

SNB_EXP_MJD = 60805.0          # adopted SN_B explosion epoch (reference)
TOP_N = 3                    # highlighted best matches per panel
MIN_ERR = 0.05               # mag floor in the weighted RMS


def band_rms(model_mjd, model_mag, obs_mjd, obs_mag, obs_err):
    """Weighted RMS of model vs obs at obs epochs inside model coverage."""
    ok = np.isfinite(model_mag)
    if ok.sum() < 2:
        return np.nan, 0
    sel = (obs_mjd >= model_mjd[ok].min()) & (obs_mjd <= model_mjd[ok].max())
    if sel.sum() < 2:
        return np.nan, int(sel.sum())
    pred = np.interp(obs_mjd[sel], model_mjd[ok], model_mag[ok])
    w = 1.0 / np.maximum(obs_err[sel], MIN_ERR) ** 2
    rms = float(np.sqrt(np.average((pred - obs_mag[sel]) ** 2, weights=w)))
    return rms, int(sel.sum())


def discover_models():
    """(name, df, is_ctl) for every cocoon model dir with a usable .tt."""
    models = []
    for mdir in sorted(glob.glob(os.path.join(GRID, "*"))):
        name = os.path.basename(mdir)
        if not os.path.isdir(mdir) or name.startswith("_"):
            continue
        df, src = load_model(mdir)
        if df is None or not len(df):
            continue
        models.append((name, df, name.endswith("_windctl"), src))
    return models


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--top", type=int, default=TOP_N,
                    help="best matches highlighted per panel")
    ap.add_argument("--out", default=OUT_FILE)
    args = ap.parse_args()

    # ---- observations ------------------------------------------------
    r_det, r_nd = load_obs_det(SNB_OBS, "r")
    g_det, _ = load_obs_det(SNB_OBS, "g")
    # load_atlas gives phase vs PEAK_MJD=60850; add it back for MJD
    o_det = load_atlas(ATLAS, band="o")["det"].copy()
    o_det["MJD"] = o_det["phase"] + PEAK_MJD
    ob_peak_mag = o_det["Mabs"].min()
    ob_peak_mjd = float(o_det.loc[o_det["Mabs"].idxmin(), "MJD"])
    for df in (r_det, g_det):
        df["phase"] = df["MJD"] - ob_peak_mjd
    r_nd = r_nd.assign(phase=r_nd["MJD"] - ob_peak_mjd) if len(r_nd) else r_nd
    print(f"obs ATLAS o peak: {ob_peak_mag:.2f} at MJD {ob_peak_mjd:.1f}")
    print(f"obs r peak: {r_det['Mabs'].min():.2f} at "
          f"phase {r_det.loc[r_det['Mabs'].idxmin(), 'phase']:.1f} d")
    print(f"obs det: r={len(r_det)} g={len(g_det)} o={len(o_det)}")

    # ---- models ------------------------------------------------------
    models = discover_models()
    rows = []
    for name, df, is_ctl, src in models:
        # anchor: model o peak -> obs o peak
        if "o" in df and df["o"].notna().sum() >= 5:
            t_pk = float(df.loc[df["o"].idxmin(), "time"])
        else:
            t_pk = float(df.loc[df["r"].idxmin(), "time"])
        mjd = df["time"].values - t_pk + ob_peak_mjd

        r_rms, n_r = band_rms(mjd, df["r"].values, r_det["MJD"].values,
                              r_det["Mabs"].values, r_det["Mabs_err"].values)
        o_rms, n_o = band_rms(mjd, df["o"].values, o_det["MJD"].values,
                              o_det["Mabs"].values, o_det["Mabs_err"].values) \
            if "o" in df else (np.nan, 0)
        g_rms, n_g = band_rms(mjd, df["g"].values, g_det["MJD"].values,
                              g_det["Mabs"].values, g_det["Mabs_err"].values)

        # combined weighted RMS over the three bands (campaign metric)
        r2 = wsum = 0.0
        for det, col in ((o_det, "o"), (r_det, "r"), (g_det, "g")):
            if col not in df:
                continue
            sel = ((det["MJD"].values >= np.nanmin(mjd))
                   & (det["MJD"].values <= np.nanmax(mjd)))
            if sel.sum() < 1:
                continue
            pred = np.interp(det["MJD"].values[sel], mjd, df[col].values)
            w = 1.0 / np.maximum(det["Mabs_err"].values[sel], MIN_ERR) ** 2
            r2 += float((w * (pred - det["Mabs"].values[sel]) ** 2).sum())
            wsum += float(w.sum())
        score = np.sqrt(r2 / wsum) if wsum else np.nan

        rows.append(dict(name=name, df=df, is_ctl=is_ctl, mjd=mjd,
                         phase=mjd - ob_peak_mjd,
                         r_rms=r_rms, n_r=n_r, o_rms=o_rms, g_rms=g_rms,
                         score=score))

    prim = sorted((r for r in rows if not r["is_ctl"]),
                  key=lambda r: (np.isnan(r["r_rms"]), r["r_rms"]))
    ctl = sorted((r for r in rows if r["is_ctl"]),
                 key=lambda r: (np.isnan(r["r_rms"]), r["r_rms"]))

    hdr = (f"\n{'model':28s} {'r_RMS':>6s} {'n_r':>4s} {'o_RMS':>6s} "
           f"{'g_RMS':>6s} {'comb':>5s}")
    for tag, group in (("ejecta models", prim), ("wind controls", ctl)):
        print(f"\n== {tag} (ranked by r-band RMS) =={hdr}")
        for r in group:
            f = lambda v: f"{v:6.2f}" if np.isfinite(v) else "    --"
            print(f"{r['name']:28s} {f(r['r_rms'])} {r['n_r']:4d} "
                  f"{f(r['o_rms'])} {f(r['g_rms'])} {f(r['score'])}")

    # ---- plot --------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(15, 7), sharex=True, sharey=True)
    x_lo, x_hi = -60.0, 90.0
    palette = ["tab:blue", "tab:purple", "teal", "brown"]

    def draw(ax, group, title):
        det = r_det[(r_det["phase"] >= x_lo) & (r_det["phase"] <= x_hi)]
        ax.errorbar(det["phase"], det["Mabs"], yerr=det["Mabs_err"],
                    fmt="o", color="tab:red", ms=6, mec="k", mew=0.8,
                    capsize=2, elinewidth=0.8, alpha=0.95, zorder=30)
        gd = g_det[(g_det["phase"] >= x_lo) & (g_det["phase"] <= x_hi)]
        ax.errorbar(gd["phase"], gd["Mabs"], yerr=gd["Mabs_err"],
                    fmt="s", color="tab:green", ms=4, mec="k", mew=0.5,
                    capsize=2, elinewidth=0.6, alpha=0.45, zorder=25)
        od = o_det[(o_det["MJD"] - ob_peak_mjd >= x_lo)
                   & (o_det["MJD"] - ob_peak_mjd <= x_hi)]
        ax.errorbar(od["MJD"] - ob_peak_mjd, od["Mabs"],
                    yerr=od["Mabs_err"], fmt="^", color="darkorange",
                    ms=4, mec="k", mew=0.5, capsize=2, elinewidth=0.6,
                    alpha=0.35, zorder=24)
        if len(r_nd):
            nd = r_nd[(r_nd["phase"] >= x_lo) & (r_nd["phase"] <= x_hi)]
            ax.scatter(nd["phase"], nd["Mabs_lim"], marker="v",
                       facecolor="none", edgecolor="tab:red", s=28,
                       alpha=0.5, zorder=20)
        ax.axvline(SNB_EXP_MJD - ob_peak_mjd, color="k", ls=":", lw=1.1,
                   zorder=1)

        best = [r for r in group if np.isfinite(r["r_rms"])][:args.top]
        best_names = {r["name"] for r in best}
        for r in group:
            if r["name"] in best_names:
                continue
            ax.plot(r["phase"], r["df"]["r"].values, color="0.6", lw=0.8,
                    alpha=0.5, zorder=5)
        for i, r in enumerate(best):
            ax.plot(r["phase"], r["df"]["r"].values,
                    color=palette[i % len(palette)], lw=2.2, zorder=15,
                    label=f"{r['name']}  (r RMS {r['r_rms']:.2f})")
        ax.set_xlim(x_lo, x_hi)
        ax.set_title(title, fontsize=11)
        ax.set_xlabel("days since ATLAS o peak (MJD "
                      f"{ob_peak_mjd:.0f})")
        ax.legend(fontsize=8, loc="upper left", framealpha=0.9)

    draw(axes[0], prim, "SN_B_cocoon: shell+wind+ejecta - model r band")
    draw(axes[1], ctl, "paired wind controls (no ejecta)")
    axes[0].set_ylabel("absolute magnitude")

    y_top = min(r_det["Mabs"].min(), o_det["Mabs"].min()) - 0.4
    for r in rows:
        in_win = r["df"]["r"][(r["phase"] >= x_lo) & (r["phase"] <= x_hi)]
        if in_win.notna().any():
            y_top = min(y_top, in_win.min() - 0.3)
    axes[0].set_ylim(-10.0, max(y_top, -20.5))

    handles = [
        Line2D([], [], marker="o", color="tab:red", mec="k", ms=7, ls="",
               label="r det (PSF)"),
        Line2D([], [], marker="v", color="tab:red", mfc="none", ms=7,
               ls="", label="r 5-sigma limit"),
        Line2D([], [], marker="s", color="tab:green", mec="k", ms=6,
               ls="", alpha=0.5, label="g det (PSF)"),
        Line2D([], [], marker="^", color="darkorange", mec="k", ms=6,
               ls="", alpha=0.5, label="o det (ATLAS)"),
        Line2D([], [], color="0.6", lw=1.0, label="other models r"),
        Line2D([], [], color="k", ls=":", lw=1.1,
               label="SN_B explosion (est.)"),
    ]
    fig.legend(handles=handles, ncol=6, loc="upper center",
               frameon=False, fontsize=9, bbox_to_anchor=(0.5, 0.955))
    fig.suptitle("SN_B_cocoon r-band lightcurves vs SN 2024pba "
                 "(SDSS Mr, o-peak anchored)", fontsize=12, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(args.out, dpi=160)
    print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()
