#!/usr/bin/env python3
"""Two-event lightcurve of SN 2024pba: SN_A model + SN_B model vs obs.

Follows the user's fit strategy:
  1. SN_A model (STELLA .tt) is time-shifted to best-fit the initial
     r-band detections (SN_A epoch, MJD < FIT_MAX_MJD).
  2. SN_B model .tt o-band peak is anchored to the observed ATLAS o-band
     peak (peak alignment, as in plot_sn_b_grid_vs_obs.py).
  3. The post-peak decay is then compared without further adjustment.

Conventions (shared with plot_sn_a_pba_grid_vs_obs.py /
plot_sn_b_grid_vs_obs.py):
  - DM = 34.52 (DL = 80.4 Mpc), E(B-V) = 0.15
  - A_g = 0.495, A_r = 0.343, A_i = 0.255, A_o = 0.36
  - Observed second-event peak MJD = 60850 (TNS astronote 2025-224)
  - SN_A explosions anchored near MJD 60488 unless the r-band fit
    prefers a shift.
  - 24-col .tt: SDSS g/r/i at cols 13-15, ATLAS c/o at cols 20-21 (AB).
    14-col .tt: Johnson V/R/I used as g/r/i proxies (dotted).

Usage:
  python scripts/plot_sn_b_par_vs_obs.py <sn_b_model_dir> [--out PNG]
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_sn_a_pba_grid_vs_obs import (
    BASE, DM, EXP_MJD, FIT_MAX_MJD, load_model,
)
from plot_sn_b_grid_vs_obs import (
    A_O, A_G, A_R, PEAK_MJD, load_atlas,
)

A_I = 1.698 * 0.15  # ~0.255, same convention as plot_sn_a_pba_grid_vs_obs

SNA_MODEL = os.path.join(
    BASE, "archive/old_sn_a_models/SN_A/SN_A_g_T20k_Ni06")
SNA_OBS = os.path.join(BASE, "observational_data", "SN_A_lightcurve.dat")
SNB_OBS = os.path.join(BASE, "observational_data", "SN_B_lightcurve.dat")
ATLAS = os.path.join(BASE, "observational_data",
                     "SN2024pba_atclean_cleaned.csv")

# SN_B: adopted explosion epoch used only for reference marking; the
# model curve is peak-aligned to the ATLAS o peak instead.
SNB_EXP_MJD = 60805.0
# ATLAS o band passes only o-band points through this fit comparison.
DECAY_WINDOW = (30.0, 90.0)      # days after o peak used for decay check
SNA_R_FIT_WINDOW = (60488.0, FIT_MAX_MJD)


def load_obs_det(path, band):
    """Load one band from a lightcurve CSV -> (detections, limits).

    Same Mabs convention as plot_sn_a_pba_grid_vs_obs.load_observations
    but tolerant of bands absent from the file (SN_B file has r/g only).
    """
    A = {"g": A_G, "r": A_R, "i": A_I}[band]
    raw = pd.read_csv(path)
    sub = raw[raw["Filter"] == band].copy()
    det = sub[sub["Mag"].apply(
        lambda x: x != "-" and not pd.isna(x))].copy()
    if len(det):
        det["Mag"] = det["Mag"].astype(float)
        det["Error"] = det["Error"].astype(float)
        det["Mabs"] = det["Mag"] - DM - A
        det["Mabs_err"] = det["Error"]
    else:
        det = det.assign(Mabs=pd.Series(dtype=float),
                         Mabs_err=pd.Series(dtype=float))
    nd = sub[sub["Mag"].apply(
        lambda x: x == "-" or pd.isna(x))].copy()
    if len(nd):
        nd["Mabs_lim"] = nd["Limit_5p0S2N"].astype(float) - DM - A
    else:
        nd = nd.assign(Mabs_lim=pd.Series(dtype=float))
    return det, nd


def best_time_shift(model_t, model_mag, obs_mjd, obs_mag, obs_t0,
                    dmax=8.0, dstep=0.1):
    """Shift s so that model MJD = obs_t0 + s + model_t best matches obs.

    Returns (shift, rms, obs_used).  Interpolates the model onto each
    observation MJD inside the model's time coverage.
    """
    best = None
    for s in np.arange(-dmax, dmax + dstep, dstep):
        mjd = obs_t0 + s + model_t
        lo, hi = mjd.min(), mjd.max()
        sel = (obs_mjd >= lo) & (obs_mjd <= hi)
        if sel.sum() < 3:
            continue
        pred = np.interp(obs_mjd[sel], mjd, model_mag)
        rms = float(np.sqrt(np.mean((pred - obs_mag[sel]) ** 2)))
        if best is None or rms < best[0]:
            best = (rms, s, sel)
    if best is None:
        return 0.0, np.nan, None
    return best[1], best[0], best[2]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("sn_b_dir", help="SN_B model dir (contains stella/*.tt)")
    ap.add_argument("--sn-a-dir", default=SNA_MODEL)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    snb_name = os.path.basename(args.sn_b_dir.rstrip("/"))
    out = args.out or os.path.join(
        args.sn_b_dir, "plots", f"{snb_name}_vs_obs.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)

    # ---- load models -------------------------------------------------
    dfa, srca = load_model(args.sn_a_dir)
    dfb, srcb = load_model(args.sn_b_dir)
    if dfb is None:
        sys.exit(f"no .tt lightcurve under {args.sn_b_dir}")

    # ---- load observations -------------------------------------------
    sna_r, sna_r_nd = load_obs_det(SNA_OBS, "r")
    sna_g, sna_g_nd = load_obs_det(SNA_OBS, "g")
    sna_i, sna_i_nd = load_obs_det(SNA_OBS, "i")
    snb_r, snb_r_nd = load_obs_det(SNB_OBS, "r")
    snb_g, snb_g_nd = load_obs_det(SNB_OBS, "g")
    atlas = load_atlas(ATLAS, band="o")
    o_det, o_nd = atlas["det"].copy(), atlas["nd"].copy()
    # load_atlas returns phase = MJD - PEAK_MJD (no MJD col); restore it
    if len(o_det):
        o_det["MJD"] = o_det["phase"] + PEAK_MJD
    if len(o_nd):
        o_nd["MJD"] = o_nd["phase"] + PEAK_MJD

    # ---- SN_A: best-fit shift to initial r-band ----------------------
    sna_shift = sna_rms = np.nan
    sna_sel = None
    if dfa is not None:
        sel = ((sna_r["MJD"] >= SNA_R_FIT_WINDOW[0])
               & (sna_r["MJD"] <= SNA_R_FIT_WINDOW[1]))
        sna_shift, sna_rms, sna_sel = best_time_shift(
            dfa["time"].values, dfa["r"].values,
            sna_r["MJD"].values[sel], sna_r["Mabs"].values[sel],
            EXP_MJD)
        dfa_mjd = EXP_MJD + sna_shift + dfa["time"].values
        print(f"[SN_A] r-band best shift = {sna_shift:+.1f} d, "
              f"rms = {sna_rms:.3f} mag "
              f"({int(sna_sel.sum())} dets in window)")

    # ---- SN_B: o-peak alignment --------------------------------------
    ob_peak_idx = o_det["Mabs"].idxmin() if len(o_det) else None
    ob_peak_mjd = (o_det.loc[ob_peak_idx, "MJD"]
                   if ob_peak_idx is not None else np.nan)
    ob_peak_mag = (o_det.loc[ob_peak_idx, "Mabs"]
                   if ob_peak_idx is not None else np.nan)

    mod_peak_idx = dfb["o"].idxmin() if "o" in dfb else dfb["r"].idxmin()
    mod_peak_t = dfb.loc[mod_peak_idx, "time"]
    mod_peak_mag = dfb.loc[mod_peak_idx,
                          "o" if "o" in dfb else "r"]
    snb_t0 = (ob_peak_mjd - mod_peak_t) if np.isfinite(ob_peak_mjd) \
        else SNB_EXP_MJD
    dfb_mjd = snb_t0 + dfb["time"].values
    print(f"[SN_B] model o peak {mod_peak_mag:.2f} at t={mod_peak_t:.1f}d"
          f" -> MJD {snb_t0 + mod_peak_t:.1f}; "
          f"obs o peak {ob_peak_mag:.2f} at MJD {ob_peak_mjd:.1f}; "
          f"implied SN_B start MJD {snb_t0:.1f} "
          f"(adopted {SNB_EXP_MJD:.0f})")

    # ---- decay metric ------------------------------------------------
    ph_lo, ph_hi = DECAY_WINDOW
    sel_o = ((o_det["MJD"] - ob_peak_mjd >= ph_lo)
             & (o_det["MJD"] - ob_peak_mjd <= ph_hi))
    dec_msg = "n/a"
    if sel_o.sum() >= 2:
        ph = o_det["MJD"].values[sel_o] - ob_peak_mjd
        mo = o_det["Mabs"].values[sel_o]
        obs_rate = np.polyfit(ph, mo, 1)[0] * 100.0
        mph = dfb_mjd - (snb_t0 + mod_peak_t)
        sel_m = (mph >= ph_lo) & (mph <= ph_hi)
        if sel_m.sum() >= 2:
            mod_rate = np.polyfit(
                mph[sel_m], dfb["o"].values[sel_m], 1)[0] * 100.0
            dec_msg = (f"obs {obs_rate:+.2f} vs model {mod_rate:+.2f} "
                       f"mag/100d (o, +{ph_lo:.0f}-{ph_hi:.0f}d)")
        else:
            dec_msg = f"obs {obs_rate:+.2f} mag/100d; model n/a"
    print(f"[decay] {dec_msg}")

    # ---- plot --------------------------------------------------------
    fig, ax = plt.subplots(figsize=(13, 7))

    # models
    if dfa is not None:
        pxy = srca != "tt_sdss"
        for band, col, lbl in (("r", "tab:red", "SN_A r" + (" (V/R proxy)" if pxy else "")),
                               ("g", "tab:green", None),
                               ("i", "goldenrod", None)):
            if band in dfa:
                ax.plot(dfa_mjd, dfa[band], color=col,
                        ls=":" if pxy else "-", lw=1.2,
                        alpha=0.5 if lbl is None else 0.9,
                        label=lbl)
    for band, col, lbl in (("o", "darkorange", f"{snb_name} o"),
                           ("r", "tab:red", f"{snb_name} r"),
                           ("g", "tab:green", None)):
        if band in dfb:
            ax.plot(dfb_mjd, dfb[band], color=col, lw=1.6,
                    ls="-" if band == "o" else "--",
                    alpha=0.9 if lbl else 0.45, label=lbl)

    # obs detections
    def det(df, mjdcol, ycol, ecol, **kw):
        if len(df):
            ax.errorbar(df[mjdcol], df[ycol], yerr=df[ecol],
                        fmt="o", ms=4, capsize=2, lw=0.8, **kw)

    det(sna_r, "MJD", "Mabs", "Mabs_err", color="tab:red",
        label="obs r (SN_A)")
    det(sna_g, "MJD", "Mabs", "Mabs_err", color="tab:green",
        label="obs g (SN_A)")
    det(sna_i, "MJD", "Mabs", "Mabs_err", color="goldenrod",
        label="obs i (SN_A)")
    det(snb_r, "MJD", "Mabs", "Mabs_err", color="tab:red",
        mfc="none", label="obs r (SN_B)")
    det(snb_g, "MJD", "Mabs", "Mabs_err", color="tab:green",
        mfc="none", label="obs g (SN_B)")
    if len(o_det):
        ax.errorbar(o_det["MJD"], o_det["Mabs"], yerr=o_det["Mabs_err"],
                    fmt="^", ms=4.5, capsize=2, lw=0.8,
                    color="darkorange", label="obs o (ATLAS)")
    # upper limits (5-sigma) - sparse, low alpha
    for df, col, ycol in ((sna_r_nd, "tab:red", "Mabs_lim"),
                          (snb_r_nd, "tab:red", "Mabs_lim"),
                          (o_nd, "darkorange", "Mabs_lim")):
        if len(df):
            ax.scatter(df["MJD"], df[ycol], marker="v", s=12,
                       color=col, alpha=0.25, lw=0)

    # fit-target shading
    ax.axvspan(SNA_R_FIT_WINDOW[0], SNA_R_FIT_WINDOW[1],
               color="tab:red", alpha=0.05)
    ax.annotate("SN_A r fit window", xy=(60516, 0.02),
                xycoords=("data", "axes fraction"), fontsize=8,
                color="tab:red", ha="center")
    ax.axvline(ob_peak_mjd, color="darkorange", ls=":", lw=1)
    ax.annotate(f"ATLAS o peak\nMJD {ob_peak_mjd:.0f}",
                xy=(ob_peak_mjd, 0.96), xycoords=("data", "axes fraction"),
                fontsize=8, color="darkorange", ha="center", va="top")

    ax.set_ylim(-8, -20)  # inverted: bright (-20) at top
    ax.set_xlim(60480, 61210)
    ax.set_xlabel("MJD")
    ax.set_ylabel("Absolute magnitude (AB)")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, ncol=3, loc="lower left")

    txt = (f"SN_A shift {sna_shift:+.1f} d, rms {sna_rms:.2f} mag | "
           f"o peak: model {mod_peak_mag:.2f} vs obs {ob_peak_mag:.2f} | "
           f"decay: {dec_msg}")
    ax.set_title(f"{snb_name} + SN_A_g_T20k_Ni06 vs SN 2024pba\n{txt}",
                 fontsize=10)
    plt.tight_layout()
    plt.savefig(out, dpi=150, bbox_inches="tight")
    print("saved", out)


if __name__ == "__main__":
    main()
