#!/usr/bin/env python3
"""Per-model SN_B comparison plots for a grid directory.

For every model dir containing a usable STELLA .tt, writes
<model_dir>/plots/<name>_vs_obs.png:

  - model o (solid) and r (dashed) absolute magnitudes
  - observed ATLAS o / PSF r / PSF g detections with errors
  - model time axis shifted so the model o-peak lands within +/-10 d of
    the observed o-peak MJD 60850.6 -- the bounded offset is chosen by
    RMS against the obs o detections; magnitudes are never rescaled
  - annotated fit targets: r plateau -15.59 (MJD 60810-60830 window),
    o peak -18.54, rise ~19 d, decay +2.3 mag/100 d guide

Usage:
  python scripts/plot_broad_model_vs_obs.py [GRID_DIR]
  (default: SN_B_broad_pilot)
"""
import glob
import os
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_sn_a_pba_grid_vs_obs import BASE, load_model  # noqa: E402
from plot_sn_b_grid_vs_obs import PEAK_MJD, load_atlas  # noqa: E402
from plot_sn_b_par_vs_obs import SNB_OBS, load_obs_det  # noqa: E402
from plot_sn_b_all_vs_obs import fit_phase  # noqa: E402

DEFAULT_GRID = os.path.join(BASE, "SN_B_broad_pilot")
R_PLATEAU_TARGET = -15.59
O_PEAK_TARGET = -18.54
DECAY_TARGET = 2.3  # mag/100 d
PLATEAU_WIN = (60810.0, 60830.0)


def plot_one(mdir, obs):
    name = os.path.basename(mdir)
    df, _ = load_model(mdir)
    if df is None or "o" not in df or df["o"].notna().sum() < 5:
        print(f"  {name}: no usable .tt (skipped)")
        return None
    o_det, r_det, g_det, ob_peak_mjd = obs

    t = df["time"].values.astype(float)
    o = df["o"].values.astype(float)
    ipk = int(np.nanargmin(o))
    det_o = {"o": (o_det["MJD"].values, o_det["Mabs"].values)}
    d_phase = fit_phase(t, {"o": o}, det_o, o_det)
    xday = d_phase - t[ipk] + t     # days relative to the observed o peak
    x_win = (-60.0, 150.0)

    fig, ax = plt.subplots(figsize=(10.5, 6.5))
    model_handles = []
    if "r" in df:
        model_handles += ax.plot(xday, df["r"].values, color="tab:red",
                                 lw=2.2, zorder=5, label=name)

    band_handles = [ax.errorbar(
        o_det["MJD"] - ob_peak_mjd, o_det["Mabs"], yerr=o_det["Mabs_err"],
        fmt="^", ms=4.5, capsize=2, lw=0.7, mfc="orange", mec="k",
        ecolor="orange", label="obs o", zorder=10)]
    band_handles.append(ax.errorbar(
        r_det["MJD"] - ob_peak_mjd, r_det["Mabs"], yerr=r_det["Mabs_err"],
        fmt="o", ms=3.5, capsize=2, lw=0.7, mfc="tab:red", mec="k",
        ecolor="tab:red", label="obs r", zorder=10))
    band_handles.append(ax.errorbar(
        g_det["MJD"] - ob_peak_mjd, g_det["Mabs"], yerr=g_det["Mabs_err"],
        fmt="s", ms=3.5, capsize=2, lw=0.7, mfc="tab:green", mec="k",
        ecolor="tab:green", label="obs g", zorder=10))

    ax.axhspan(R_PLATEAU_TARGET - 0.1, R_PLATEAU_TARGET + 0.1,
               color="tab:red", alpha=0.12, zorder=1)
    ax.annotate("r plateau -15.6", xy=(-58, R_PLATEAU_TARGET - 0.3),
                fontsize=8, color="tab:red")
    ax.axhline(O_PEAK_TARGET, color="darkorange", ls=":", lw=1, alpha=0.7)
    ax.axvline(0.0, color="darkorange", ls=":", lw=1, alpha=0.5)
    ax.annotate(f"o peak {O_PEAK_TARGET:.2f}",
                xy=(-58, O_PEAK_TARGET + 0.08), fontsize=8,
                color="darkorange")
    xg = np.array([0.0, 60.0])
    ax.plot(xg, O_PEAK_TARGET + DECAY_TARGET * xg / 100.0,
            color="darkorange", ls="-.", lw=1, alpha=0.6)
    ax.annotate("decay +2.3 mag/100d", xy=(15, -18.1),
                fontsize=8, color="darkorange", rotation=4)

    # decay of the model itself over covered post-peak window
    # model names: single column top right; band keys: one row above axes
    leg_m = ax.legend(
        handles=model_handles,
        labels=[h.get_label() for h in model_handles],
        loc="upper right", fontsize=8, framealpha=0.9)
    ax.add_artist(leg_m)
    ax.legend(handles=band_handles,
              labels=[h.get_label() for h in band_handles],
              loc="lower center", bbox_to_anchor=(0.5, 1.0),
              ncol=len(band_handles), fontsize=8, frameon=False,
              columnspacing=1.2, handletextpad=0.5)
    # y limits framed on the observed detections inside the x window
    obs_mags = np.concatenate([
        d["Mabs"].values[(d["MJD"] - ob_peak_mjd >= x_win[0])
                         & (d["MJD"] - ob_peak_mjd <= x_win[1])]
        for d in (o_det, r_det, g_det)])
    obs_mags = obs_mags[np.isfinite(obs_mags)]
    ax.set_ylim(obs_mags.max() + 0.4, obs_mags.min() - 0.4)
    ax.set_xlim(*x_win)
    ax.set_xlabel("days since observed o peak")
    ax.set_ylabel("Absolute magnitude (AB)")
    ax.grid(alpha=0.3)
    plt.tight_layout(rect=[0, 0, 1, 0.92])
    outdir = os.path.join(mdir, "plots")
    os.makedirs(outdir, exist_ok=True)
    out = os.path.join(outdir, f"{name}_vs_obs.png")
    plt.savefig(out, dpi=140)
    plt.close(fig)
    return out


def main():
    grid = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_GRID
    if not os.path.isabs(grid):
        grid = os.path.join(BASE, grid)
    atlas = load_atlas(os.path.join(BASE, "observational_data",
                                    "SN2024pba_atclean_cleaned.csv"),
                       band="o")
    o_det = atlas["det"].copy()
    o_det["MJD"] = o_det["phase"] + PEAK_MJD
    ob_peak_mjd = float(o_det.loc[o_det["Mabs"].idxmin(), "MJD"])
    r_det, _ = load_obs_det(SNB_OBS, "r")
    g_det, _ = load_obs_det(SNB_OBS, "g")
    obs = (o_det, r_det, g_det, ob_peak_mjd)

    n = 0
    for mdir in sorted(glob.glob(os.path.join(grid, "SB_*"))):
        if not os.path.isdir(mdir):
            continue
        out = plot_one(mdir, obs)
        if out:
            print("saved", out)
            n += 1
    print(f"{n} plots written")


if __name__ == "__main__":
    main()
