#!/usr/bin/env python3
"""Plot ALL completed SN_B models (grey) with the best-fitting 5 in colour.

Every model with a usable STELLA .tt under the given grid dirs is scored
by region-wise RMS residuals against the observed detections: the early
r-band plateau (~ -15.8), the ATLAS o main-peak region (~ -18), the
r-band second-peak region (~ -18.2), and a weaker all-data term over
every SN_B-epoch detection in o, c, g and r.  r/g decay and rise time
survive only as tiebreakers.  A model whose light curve does not cover
or match the r plateau AT ITS FITTED PHASE is vetoed from the ranking
entirely (it cannot be shifted to match that region); vetoed models are
still drawn grey.

Model x-axis convention: the model o-peak is anchored to the observed
o-peak, then a bounded +/-10 d offset is fitted by minimising the mean
per-band RMS over all available detections (o, c, g, r).  Magnitudes are
never rescaled.

Usage:
    python scripts/plot_sn_b_all_vs_obs.py [grid_dir ...]
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import glob

import matplotlib.pyplot as plt
import numpy as np

from plot_sn_a_pba_grid_vs_obs import load_model  # noqa: E402
from plot_sn_b_grid_vs_obs import PEAK_MJD, load_atlas  # noqa: E402
from plot_sn_b_par_vs_obs import SNB_OBS, load_obs_det  # noqa: E402
from analyze_sn_b_trends import model_metrics, obs_targets  # noqa: E402

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GRID_DIRS = sys.argv[1:] or [
    os.path.join(BASE, "SN_B_broad_grid"),
]
OUT = os.path.join(BASE, "SN_B_broad_grid", "plots",
                   "sn_b_all_models_vs_obs.png")
TOP_N = 5
# Scoring is region-wise trend agreement, not single-value matching: the
# model light curve is interpolated at the OBSERVED detection MJDs inside
# each region and scored by the RMS magnitude residual.
#   Region 1 (x4): r-band plateau detections, MJD 60810-60831 (~ -15.8)
#   Region 2 (x4): ATLAS o detections around the main peak, +/-15 d
#                  (~ -18.2 .. -18.5)
#   Region 3 (x4): ZTF r second-peak detections, MJD 60860-60905
#                  (~ -18.2 .. -17.7; the r data have a 60831-60867 gap,
#                  so this is an explicit window on the 5 peak-cluster
#                  dets, not a +/-half window).  Without this a model can
#                  fit the plateau and the o peak while its r second peak
#                  sits ~1 mag too faint - the failure the plot shows.
#   All-data (x2): mean per-band RMS over every SN_B-epoch detection in
#                  o, c, g and r - the global-trend term
# r/g decay (x0.5) and rise time (x0.05/d) survive only as tiebreakers.
W_REG1 = 4.0
W_REG2 = 4.0
W_REG3 = 4.0
W_ALL = 3.0
W_RISE = 0.05
W_DECAY = 0.5
# Hard gate on the early r plateau: at the model's fitted phase it must
# cover at least VETO_MIN_COV of the 13 plateau detections AND reach
# rms_r1 <= VETO_RMS.  A model that only matches the plateau at a shift
# inconsistent with the rest of the light curve cannot be shifted to
# match this region - it is removed from the ranking, not penalised.
VETO_RMS = 1.5
VETO_MIN_COV = 6
# Region windows in MJD (observer frame, after the bounded phase fit).
REG1_WIN = (60810.0, 60831.0)   # early r detections (~ -15.8)
REG2_HALF = 15.0              # +/- days around the ATLAS o peak (~ -18)
REG3_WIN = (60860.0, 60905.0)   # r second-peak cluster (~ -18.2..-17.7)
# Penalty applied when a model's shifted light curve does not cover enough
# detections in a region to be tested at all.
REG_MISS = 2.5
# Models whose shifted light curve does not reach back into the obs plateau
# window have no prediction for the early r emission at all -- that is a
# failed criterion, not an excused one. Count it as a PLAT_MISS-mag miss.
PLAT_MISS = 2.5
# SN_B-only plotting: the ATLAS o-band file also covers the SN_A event and
# the inter-event plateau (o detections end at MJD 60717 and resume at
# 60839.6, a 122 d gap).  Everything below this cut is SN_A data.
SNB_MIN_MJD = 60800.0


def fit_phase(t_mod, curves, det_sets, o_det, max_shift=10.0):
    """Bounded phase fit on all available observed data.

    The model o-peak is anchored to the observed o-peak, then the offset
    (<= max_shift days either way) minimising the mean of the per-band RMS
    residuals - model interpolated at each detection MJD - is returned.
    Per-band averaging keeps the 71 ATLAS o points from swamping the
    13 plateau r points; a band contributes only where >=3 detections are
    covered at that shift.
    """
    anchor0 = o_det["MJD"].iloc[int(np.argmin(o_det["Mabs"]))] \
        - t_mod[int(np.nanargmin(curves["o"]))]
    best_d, best_rms = 0.0, np.inf
    for d in np.linspace(-max_shift, max_shift, 81):
        per_band = []
        for band, mag in curves.items():
            ok = np.isfinite(mag)
            rms, _ = region_rms(t_mod[ok] + anchor0 + d, mag[ok],
                                det_sets[band][0], det_sets[band][1],
                                win=None)
            if np.isfinite(rms):
                per_band.append(rms)
        if not per_band:
            continue
        rms = float(np.mean(per_band))
        if rms < best_rms:
            best_rms, best_d = rms, float(d)
    return best_d


def region_rms(mod_mjd, mod_mag, det_mjd, det_mag, win, min_pts=3):
    """RMS magnitude residual of a model vs detections inside a window.

    mod_mjd is the model's phase-shifted MJD axis; the model is interpolated
    at each detection MJD inside win (win=None uses every detection).
    Returns (rms, n_covered); rms is nan when the model covers fewer than
    min_pts of the window's detections.
    """
    dm, dg = np.asarray(det_mjd), np.asarray(det_mag)
    if win is not None:
        sel = (dm >= win[0]) & (dm <= win[1])
        dm, dg = dm[sel], dg[sel]
    if len(dm) == 0:
        return np.nan, 0
    cov = (dm >= mod_mjd.min()) & (dm <= mod_mjd.max())
    if cov.sum() < min_pts:
        return np.nan, int(cov.sum())
    interp = np.interp(dm[cov], mod_mjd, mod_mag)
    return float(np.sqrt(np.mean((interp - dg[cov]) ** 2))), int(cov.sum())


def region_bias(mod_mjd, mod_mag, det_mjd, det_mag, win):
    """Mean signed residual (model - obs) over the covered detections.

    A uniform offset is a luminosity error that RMS understates when it is
    consistent; the sign tells which way the model is off.
    """
    dm, dg = np.asarray(det_mjd), np.asarray(det_mag)
    sel = (dm >= win[0]) & (dm <= win[1])
    dm, dg = dm[sel], dg[sel]
    cov = (dm >= mod_mjd.min()) & (dm <= mod_mjd.max())
    if cov.sum() < 3:
        return np.nan
    return float(np.mean(np.interp(dm[cov], mod_mjd, mod_mag) - dg[cov]))


def score(met, tg):
    r1 = met["rms_r1"] if np.isfinite(met["rms_r1"]) else REG_MISS
    r2 = met["rms_o2"] if np.isfinite(met["rms_o2"]) else REG_MISS
    r3 = met["rms_r3"] if np.isfinite(met["rms_r3"]) else REG_MISS
    ra = met["rms_all"] if np.isfinite(met["rms_all"]) else REG_MISS
    terms = [W_REG1 * r1, W_REG2 * r2, W_REG3 * r3, W_ALL * ra,
             W_RISE * abs(met["t_rise"] - tg["t_rise"])]
    if np.isfinite(met["decay"]):
        terms.append(W_DECAY * abs(met["decay"] - tg["decay"]))
    return float(np.sqrt(np.sum(np.square(terms))))


def collect_models():
    out = {}
    for grid in GRID_DIRS:
        for mdir in sorted(glob.glob(os.path.join(grid, "*"))):
            if not os.path.isdir(mdir):
                continue
            name = os.path.basename(mdir)
            if name.startswith(("plots", "recon")) or name in out:
                continue
            df, _ = load_model(mdir)
            if df is None or "o" not in df or df["o"].notna().sum() < 10:
                continue
            out[name] = (mdir, df)
    return out


def main():
    tg = obs_targets()
    atlas = load_atlas(os.path.join(BASE, "observational_data",
                                    "SN2024pba_atclean_cleaned.csv"), band="o")
    o = atlas["det"].copy()
    o["MJD"] = o["phase"] + PEAK_MJD
    o = o[o["MJD"] >= SNB_MIN_MJD].copy()
    cdet = load_atlas(os.path.join(BASE, "observational_data",
                                   "SN2024pba_atclean_cleaned.csv"),
                      band="c")["det"].copy()
    cdet["MJD"] = cdet["phase"] + PEAK_MJD
    cdet = cdet[cdet["MJD"] >= SNB_MIN_MJD]
    r_det, _ = load_obs_det(SNB_OBS, "r")
    r_det = r_det[r_det["MJD"] >= SNB_MIN_MJD]
    g_det, _ = load_obs_det(SNB_OBS, "g")
    g_det = g_det[g_det["MJD"] >= SNB_MIN_MJD]
    # All available SN_B-epoch detections, model band -> (MJD, Mabs)
    det_sets = {"o": (o["MJD"].values, o["Mabs"].values),
                "c": (cdet["MJD"].values, cdet["Mabs"].values),
                "r": (r_det["MJD"].values, r_det["Mabs"].values),
                "g": (g_det["MJD"].values, g_det["Mabs"].values)}
    reg2_win = (tg["pk_mjd"] - REG2_HALF, tg["pk_mjd"] + REG2_HALF)
    reg3_win = REG3_WIN
    n_reg1 = int(((r_det["MJD"] >= REG1_WIN[0])
                  & (r_det["MJD"] <= REG1_WIN[1])).sum())
    n_reg2 = int(((o["MJD"] >= reg2_win[0]) & (o["MJD"] <= reg2_win[1])).sum())
    n_reg3 = int(((r_det["MJD"] >= reg3_win[0])
                  & (r_det["MJD"] <= reg3_win[1])).sum())
    n_all = sum(len(v[0]) for v in det_sets.values())
    models = collect_models()
    print(f"{len(models)} models with usable .tt")
    print(f"obs: r_pk1={tg['r_pk1']:.2f} r_pk2={tg['r_pk2']:.2f} "
          f"o_pk={tg['o_pk']:.2f} t_rise={tg['t_rise']:.0f}d "
          f"decay={tg['decay']:+.2f}")
    print(f"regions: r plateau {REG1_WIN[0]:.0f}-{REG1_WIN[1]:.0f} "
          f"({n_reg1} dets), o peak +/-{REG2_HALF:.0f}d ({n_reg2} dets), "
          f"r peak {reg3_win[0]:.0f}-{reg3_win[1]:.0f} ({n_reg3} dets), "
          f"all-data {n_all} dets in o/c/r/g")
    print(f"veto: plateau needs cov>={VETO_MIN_COV}/{n_reg1} and "
          f"rms_r1<={VETO_RMS} at the fitted phase")

    scored, vetoed = [], []
    for name, (mdir, df) in models.items():
        t = df["time"].values.astype(float)
        curves = {b: df[b].values.astype(float)
                  for b in det_sets if b in df}
        d_phase = fit_phase(t, curves, det_sets, o)
        pk_frame = tg["pk_mjd"] + d_phase      # shifted obs frame
        met = model_metrics(df, pk_frame)
        met["d_phase"] = d_phase
        # Model MJD axis under the fitted phase
        mjd_mod = pk_frame - met["t_pk"] + t
        met["rms_r1"], met["n_r1"] = region_rms(
            mjd_mod[np.isfinite(curves["r"])],
            curves["r"][np.isfinite(curves["r"])],
            *det_sets["r"], win=REG1_WIN)
        met["rms_o2"], met["n_o2"] = region_rms(
            mjd_mod[np.isfinite(curves["o"])],
            curves["o"][np.isfinite(curves["o"])],
            *det_sets["o"], win=reg2_win)
        met["rms_r3"], met["n_r3"] = region_rms(
            mjd_mod[np.isfinite(curves["r"])],
            curves["r"][np.isfinite(curves["r"])],
            *det_sets["r"], win=reg3_win)
        met["bias_r1"] = region_bias(
            mjd_mod[np.isfinite(curves["r"])],
            curves["r"][np.isfinite(curves["r"])],
            *det_sets["r"], win=REG1_WIN)
        met["bias_o2"] = region_bias(
            mjd_mod[np.isfinite(curves["o"])],
            curves["o"][np.isfinite(curves["o"])],
            *det_sets["o"], win=reg2_win)
        met["bias_r3"] = region_bias(
            mjd_mod[np.isfinite(curves["r"])],
            curves["r"][np.isfinite(curves["r"])],
            *det_sets["r"], win=reg3_win)
        per_band = []
        for band, mag in curves.items():
            ok = np.isfinite(mag)
            rms, _ = region_rms(mjd_mod[ok], mag[ok],
                                *det_sets[band], win=None)
            if np.isfinite(rms):
                per_band.append(rms)
        met["rms_all"] = float(np.mean(per_band)) if per_band else np.nan
        met["veto"] = not (met["n_r1"] >= VETO_MIN_COV
                           and np.isfinite(met["rms_r1"])
                           and met["rms_r1"] <= VETO_RMS)
        entry = (score(met, tg), name, df, met)
        (vetoed if met["veto"] else scored).append(entry)
    scored.sort(key=lambda s: s[0])
    vetoed.sort(key=lambda s: s[0])
    top = scored[:TOP_N]
    # Best r-peak matches among the vetoed (plateau failures): drawn in
    # purple so the trade-off behind the survivors' faint r second peak
    # is visible - matching the peak costs the early plateau.
    peak_picks = sorted((e for e in vetoed
                         if np.isfinite(e[3]["rms_r3"])),
                        key=lambda e: e[3]["rms_r3"])[:3]

    print(f"\n{len(vetoed)}/{len(models)} models vetoed "
          f"(cannot match the r plateau at their fitted phase)")
    print("\nTop 5 (surviving the plateau veto):")
    for s, name, _, met in top:
        print(f"  {name:34s} score={s:.2f}  "
              f"rms_r1={met['rms_r1']:.2f}({met['n_r1']}d,"
              f"bias{met['bias_r1']:+.2f}) "
              f"rms_o2={met['rms_o2']:.2f}({met['n_o2']}d,"
              f"bias{met['bias_o2']:+.2f}) "
              f"rms_r3={met['rms_r3']:.2f}({met['n_r3']}d,"
              f"bias{met['bias_r3']:+.2f}) "
              f"rms_all={met['rms_all']:.2f} "
              f"t_rise={met['t_rise']:.0f}d "
              f"decay={met['decay']:+.1f}  shift={met['d_phase']:+.1f}d")

    print("\nBest r-peak fits among vetoed models (purple in figure):")
    for s, name, _, met in peak_picks:
        r1 = met["rms_r1"] if np.isfinite(met["rms_r1"]) else -1
        o2 = met["rms_o2"] if np.isfinite(met["rms_o2"]) else -1
        print(f"  {name:34s} rms_r3={met['rms_r3']:.2f}({met['n_r3']}d,"
              f"bias{met['bias_r3']:+.2f})  rms_r1={r1:.2f}({met['n_r1']}d)"
              f"  rms_o2={o2:.2f}  shift={met['d_phase']:+.1f}d")

    fig, ax = plt.subplots(figsize=(9.5, 6.0))
    pk = tg["pk_mjd"]
    x_win = (-60.0, 150.0)      # days relative to the observed o peak

    # all models grey (vetoed included - visible but not ranked)
    for _, name, df, met in scored + vetoed:
        t = df["time"].values.astype(float)
        om = df["o"].values.astype(float)
        rm = df["r"].values.astype(float)
        mjd = (pk + met["d_phase"]
               - t[int(np.nanargmin(om))] + t)
        ax.plot(mjd - pk, rm, color="0.75", lw=0.8, alpha=0.55, zorder=1)

    # scoring regions, so what is being fit is visible
    ax.axvspan(REG1_WIN[0] - pk, REG1_WIN[1] - pk, color="tab:red",
               alpha=0.07, zorder=0)
    ax.axvspan(reg2_win[0] - pk, reg2_win[1] - pk, color="darkorange",
               alpha=0.07, zorder=0)
    ax.axvspan(reg3_win[0] - pk, reg3_win[1] - pk, color="tab:purple",
               alpha=0.07, zorder=0)
    ax.text(REG1_WIN[0] + 0.5 - pk, 0.98, "r plateau region",
            fontsize=7, color="tab:red", va="top",
            transform=ax.get_xaxis_transform())
    ax.text(reg2_win[0] + 0.5 - pk, 0.93, "o peak region",
            fontsize=7, color="darkorange", va="top",
            transform=ax.get_xaxis_transform())
    ax.text(reg3_win[0] + 0.5 - pk, 0.88, "r peak region",
            fontsize=7, color="tab:purple", va="top",
            transform=ax.get_xaxis_transform())

    # top-5 in colour: model r band, solid
    colors = plt.cm.tab10(np.arange(TOP_N))
    model_handles = []
    for (s, name, df, met), c in zip(top, colors):
        t = df["time"].values.astype(float)
        om = df["o"].values.astype(float)
        rm = df["r"].values.astype(float)
        mjd = (pk + met["d_phase"]
               - t[int(np.nanargmin(om))] + t)
        model_handles += ax.plot(mjd - pk, rm, color=c, lw=2.0,
                                 zorder=4, label=name)
    # vetoed models with the best r-peak fits, in purple
    pk_colors = plt.cm.Purples(np.linspace(0.55, 0.9, len(peak_picks)))
    for (s, name, df, met), c in zip(peak_picks, pk_colors):
        t = df["time"].values.astype(float)
        om = df["o"].values.astype(float)
        rm = df["r"].values.astype(float)
        mjd = (pk + met["d_phase"]
               - t[int(np.nanargmin(om))] + t)
        model_handles += ax.plot(mjd - pk, rm, color=c, lw=1.8,
                                 zorder=3, label=f"{name} [vetoed]")
    band_handles = [ax.errorbar(
        o["MJD"] - pk, o["Mabs"], yerr=o["Mabs_err"], fmt="o",
        ms=6, mfc="orange", mec="k", ecolor="orange", elinewidth=0.9,
        capsize=2, zorder=6, label="ATLAS o")]
    band_handles.append(ax.errorbar(
        cdet["MJD"] - pk, cdet["Mabs"], yerr=cdet["Mabs_err"],
        fmt="d", ms=5, mfc="cyan", mec="k", ecolor="cyan",
        elinewidth=0.8, capsize=2, zorder=6, label="ATLAS c"))
    for band, marker, col, lab in [("g", "s", "tab:green", "ZTF g"),
                                   ("r", "^", "tab:red", "ZTF r")]:
        det, _ = load_obs_det(SNB_OBS, band)
        det = det[det["MJD"] >= SNB_MIN_MJD]
        if det is not None and len(det):
            band_handles.append(ax.errorbar(
                det["MJD"] - pk, det["Mabs"], yerr=det["Mabs_err"],
                fmt=marker, ms=4.5, mfc=col, mec="k",
                ecolor=col, elinewidth=0.7, capsize=2,
                zorder=5, label=lab))

    ax.axhline(tg["o_pk"], ls=":", lw=0.8, color="0.4")
    ax.axhline(tg["r_plateau"], ls=":", lw=0.8, color="0.6")
    ax.text(x_win[0] + 1, tg["o_pk"] + 0.12, "obs o peak", fontsize=8,
            color="0.35")
    ax.text(x_win[0] + 1, tg["r_plateau"] - 0.18, "obs r plateau",
            fontsize=8, color="0.45")

    # y limits framed on the observed detections inside the x window
    obs_mags = np.concatenate([
        d["Mabs"].values[(d["MJD"] - pk >= x_win[0])
                         & (d["MJD"] - pk <= x_win[1])]
        for d in (o, cdet, g_det, r_det)])
    obs_mags = obs_mags[np.isfinite(obs_mags)]
    ax.set_ylim(obs_mags.max() + 0.4, obs_mags.min() - 0.4)
    ax.set_xlim(*x_win)
    ax.set_xlabel("days since observed o peak (model o-peak aligned)")
    ax.set_ylabel("Absolute magnitude")
    # model names: single column top right; band keys: one row above axes
    leg_m = ax.legend(
        handles=model_handles,
        labels=[h.get_label() for h in model_handles],
        loc="upper right", fontsize=7.5, framealpha=0.9)
    ax.add_artist(leg_m)
    ax.legend(handles=band_handles,
              labels=[h.get_label() for h in band_handles],
              loc="lower center", bbox_to_anchor=(0.5, 1.0),
              ncol=len(band_handles), fontsize=8, frameon=False,
              columnspacing=1.2, handletextpad=0.5)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    fig.savefig(OUT, dpi=200)
    print(f"\nsaved {OUT}")


if __name__ == "__main__":
    main()
