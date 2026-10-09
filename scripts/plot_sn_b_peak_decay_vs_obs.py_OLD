#!/usr/bin/env python3
"""Plot SN_B broad-grid models scored by direct o/g/r band comparison.

Companion to plot_sn_b_all_vs_obs.py.  Each model is phase-aligned by a
bounded fit on the o peak (+/-10 d) and then scored by the RMS between
the model light curve and the observed detections in each band:

  score = sqrt(rms_o^2 + rms_r^2 + rms_g^2)

where rms_X is measured at every observed detection in band X that the
model's shifted light curve covers (min 3 detections per band, else a
miss penalty).  SN_A / inter-event photometry (MJD < 60800) is excluded.

The figure shows every model's o band in grey and the top 5 in colour,
with the coloured models' r (dashed) and g (dotted) curves drawn as
well so the per-band match can be judged.

Usage:
    python scripts/plot_sn_b_peak_decay_vs_obs.py [grid_dir ...]
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
from plot_sn_b_all_vs_obs import fit_phase, region_rms  # noqa: E402

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GRID_DIRS = sys.argv[1:] or [
    os.path.join(BASE, "SN_B_broad_grid"),
]
OUT = os.path.join(BASE, "SN_B_broad_grid", "plots",
                   "sn_b_peak_decay_vs_obs.png")
TOP_N = 5
REG_MISS = 2.5    # penalty per band with too few covered detections
# SN_B-only plotting: the ATLAS o-band file also covers the SN_A event and
# the inter-event plateau (o detections end at MJD 60717 and resume at
# 60839.6, a 122 d gap).  Everything below this cut is SN_A data.
SNB_MIN_MJD = 60800.0


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
                                    "SN2024pba_atclean_cleaned.csv"),
                       band="o")
    o = atlas["det"].copy()
    o["MJD"] = o["phase"] + PEAK_MJD
    o = o[o["MJD"] >= SNB_MIN_MJD].copy()
    r_det, _ = load_obs_det(SNB_OBS, "r")
    g_det, _ = load_obs_det(SNB_OBS, "g")
    r_det = r_det[r_det["MJD"] >= SNB_MIN_MJD]
    g_det = g_det[g_det["MJD"] >= SNB_MIN_MJD]
    full_win = (SNB_MIN_MJD, np.inf)
    models = collect_models()
    print(f"{len(models)} models with usable .tt")
    print(f"obs: o_pk={tg['o_pk']:.2f} at MJD {tg['pk_mjd']:.1f}, "
          f"dets: o={len(o)} r={len(r_det)} g={len(g_det)}")

    scored = []
    for name, (mdir, df) in models.items():
        t = df["time"].values.astype(float)
        om = df["o"].values.astype(float)
        rm = df["r"].values.astype(float)
        gm = (df["g"].values.astype(float)
              if "g" in df else np.full_like(om, np.nan))
        det_o = {"o": (o["MJD"].values, o["Mabs"].values)}
        d_phase = fit_phase(t, {"o": om}, det_o, o)
        pk_frame = tg["pk_mjd"] + d_phase
        met = model_metrics(df, pk_frame)
        met["d_phase"] = d_phase
        mjd_mod = pk_frame - met["t_pk"] + t
        for band, mm, det in (("o", om, o), ("r", rm, r_det),
                              ("g", gm, g_det)):
            met[f"rms_{band}"], met[f"n_{band}"] = region_rms(
                mjd_mod[np.isfinite(mm)], mm[np.isfinite(mm)],
                det["MJD"].values, det["Mabs"].values, full_win)
        terms = [met[f"rms_{b}"] if np.isfinite(met[f"rms_{b}"])
                 else REG_MISS for b in "org"]
        s = float(np.sqrt(np.sum(np.square(terms))))
        scored.append((s, name, df, met))
    scored.sort(key=lambda s: s[0])
    top = scored[:TOP_N]

    print("\nTop 5:")
    for s, name, _, met in top:
        print(f"  {name:34s} score={s:.2f}  "
              f"rms_o={met['rms_o']:.2f}({met['n_o']}d) "
              f"rms_r={met['rms_r']:.2f}({met['n_r']}d) "
              f"rms_g={met['rms_g']:.2f}({met['n_g']}d) "
              f"shift={met['d_phase']:+.1f}d")

    fig, ax = plt.subplots(figsize=(9.5, 6.0))
    pk = tg["pk_mjd"]
    x_win = (-60.0, 150.0)      # days relative to the observed o peak

    # all models grey (r band)
    for _, name, df, met in scored:
        t = df["time"].values.astype(float)
        om = df["o"].values.astype(float)
        rm = df["r"].values.astype(float)
        mjd = (pk + met["d_phase"]
               - t[int(np.nanargmin(om))] + t)
        ax.plot(mjd - pk, rm, color="0.75", lw=0.8, alpha=0.55, zorder=1)

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

    band_handles = [ax.errorbar(
        o["MJD"] - pk, o["Mabs"], yerr=o["Mabs_err"], fmt="o",
        ms=6, mfc="orange", mec="k", ecolor="orange", elinewidth=0.9,
        capsize=2, zorder=6, label="ATLAS o")]
    for det, marker, col, lab in [(g_det, "s", "tab:green", "ZTF g"),
                                  (r_det, "^", "tab:red", "ZTF r")]:
        if det is not None and len(det):
            band_handles.append(ax.errorbar(
                det["MJD"] - pk, det["Mabs"], yerr=det["Mabs_err"],
                fmt=marker, ms=4.5, mfc=col, mec="k",
                ecolor=col, elinewidth=0.7, capsize=2,
                zorder=5, label=lab))

    ax.axhline(tg["o_pk"], ls=":", lw=0.8, color="0.4")

    # y limits framed on the observed detections inside the x window
    obs_mags = np.concatenate([
        d["Mabs"].values[(d["MJD"] - pk >= x_win[0])
                         & (d["MJD"] - pk <= x_win[1])]
        for d in (o, g_det, r_det)])
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
