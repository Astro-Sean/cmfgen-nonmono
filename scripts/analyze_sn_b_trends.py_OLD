#!/usr/bin/env python3
"""Trend analysis of the SN_B_double_peak_grid against four lightcurve criteria.

All grid models share the fixed SN_A_g_T20k_Ni06@400d ejecta; only the
shell (M, v, n_rho, t_start) and wind (Mdot, v_w) vary.

Criteria (targets measured from the observations in this script):
  1. early r plateau       - mean r over the obs plateau window (MJD 60810-60830)
                             with the model o-peak aligned to the obs o-peak
  2. rise time to 2nd peak - t_rise = t(o peak) - t(dip minimum before it);
                             also t_peak from model t=0 (explosion)
  3. 2nd peak brightness   - min o
  4. decay                 - o slope +20..+60 d after peak (mag/100 d),
                             limited by model coverage

Writes a CSV table and trend figures to SN_B_double_peak_grid/plots/.
Alignment convention: time shift only, magnitudes absolute (same as
score_sn_b_grid_rms.py).
"""
import glob
import os
import sys

import matplotlib
import numpy as np
import pandas as pd
import yaml

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_sn_a_pba_grid_vs_obs import BASE, load_model  # noqa: E402
from plot_sn_b_grid_vs_obs import PEAK_MJD, load_atlas  # noqa: E402
from plot_sn_b_par_vs_obs import SNB_OBS, load_obs_det  # noqa: E402

GRID = os.path.join(BASE, "SN_B_double_peak_grid")
OUTDIR = os.path.join(GRID, "plots")
PLATEAU_WIN = (60810.0, 60830.0)
DECAY_WIN = (20.0, 60.0)  # days after o peak


def read_params(mdir):
    c = os.path.join(mdir, "config.yaml")
    if not os.path.exists(c):
        return None
    x = (yaml.safe_load(open(c)) or {}).get("x_model_from_sn") or {}
    return dict(M=x.get("shell_mass_msun"), v=x.get("shell_velocity_kms"),
                nrho=x.get("n_rho"), ts=x.get("t_start_days"),
                mdot=x.get("wind_mdot_msun_yr"), vw=x.get("wind_velocity_kms"))


def decay_slope(t, m, t_pk):
    """Band slope after peak in mag/100 d (used for r and g).

    Uses +20..+60 d when covered; otherwise falls back to +10 d .. end of
    coverage so short double-peak runs (peak at ~57 d of a 100 d run) still
    get a measurement. Returns (slope, post-peak coverage, window used).
    """
    cov = float(t.max() - t_pk)
    lo, hi = DECAY_WIN
    if cov < hi:
        lo, hi = 10.0, cov
    sel = (t >= t_pk + lo) & (t <= t_pk + hi)
    if sel.sum() < 4 or (t[sel].max() - t[sel].min()) < 10:
        return np.nan, cov, f"{lo:.0f}-{hi:.0f}"
    p = np.polyfit(t[sel], m[sel], 1)
    return float(p[0] * 100.0), cov, f"{lo:.0f}-{hi:.0f}"


def model_metrics(df, ob_peak_mjd):
    t = df["time"].values.astype(float)
    o = df["o"].values.astype(float)
    r = df["r"].values.astype(float)
    g = df["g"].values.astype(float) if "g" in df else np.full_like(o, np.nan)
    ok = np.isfinite(o) & np.isfinite(r)
    t, o, r = t[ok], o[ok], r[ok]
    g = g[ok]
    ipk = int(np.argmin(o))
    t_pk, o_pk = t[ipk], o[ipk]
    # dip: faintest o between the first local brightness max and the peak
    pre = o[:ipk]
    if len(pre) > 5:
        i_first = int(np.argmin(pre[: max(2, ipk // 2)]))
        seg = o[i_first:ipk + 1]
        i_dip = i_first + int(np.argmax(seg))
        dip = float(o[i_dip] - o[i_first])
        t_rise = float(t_pk - t[i_dip])
    else:
        dip, t_rise = 0.0, float(t_pk)
    mjd = ob_peak_mjd - t_pk + t
    sel = (mjd >= PLATEAU_WIN[0]) & (mjd <= PLATEAU_WIN[1])
    r_plateau = float(np.mean(r[sel])) if sel.sum() >= 3 else np.nan
    # r-band first-peak magnitude: brightest r in the obs plateau window
    r_pk1 = float(np.min(r[sel])) if sel.sum() >= 3 else np.nan
    # r-band second-peak magnitude: brightest r after the o-band dip
    i_dip_end = i_dip if len(pre) > 5 else 0
    r_pk2 = float(np.min(r[i_dip_end:])) if i_dip_end < len(r) else np.nan
    # Post-peak decay measured in the ZTF bands (r and g); combined value is
    # the mean of whichever band slopes are measurable.
    slope_r, cov, win = decay_slope(t, r, t_pk)
    gok = np.isfinite(g)
    slope_g = (decay_slope(t[gok], g[gok], t_pk)[0]
               if gok.sum() > 4 else np.nan)
    slopes = [s for s in (slope_r, slope_g) if np.isfinite(s)]
    slope = float(np.mean(slopes)) if slopes else np.nan
    return dict(o_pk=float(o_pk), t_pk=float(t_pk), t_rise=t_rise, dip=dip,
                r_plateau=r_plateau, r_pk1=r_pk1, r_pk2=r_pk2,
                decay=slope, decay_r=slope_r,
                decay_g=slope_g, cov=cov, decay_win=win)


def obs_targets():
    atlas = load_atlas(os.path.join(BASE, "observational_data",
                                    "SN2024pba_atclean_cleaned.csv"), band="o")
    o = atlas["det"].copy()
    o["MJD"] = o["phase"] + PEAK_MJD
    ipk = o["Mabs"].idxmin()
    pk_mjd, pk_mag = float(o.loc[ipk, "MJD"]), float(o.loc[ipk, "Mabs"])
    r, _ = load_obs_det(SNB_OBS, "r")
    sel = (r["MJD"] >= PLATEAU_WIN[0]) & (r["MJD"] <= PLATEAU_WIN[1])
    r_plat = float(r.loc[sel, "Mabs"].mean())
    # r-band first peak: brightest detection inside the plateau window
    r_pk1 = float(r.loc[sel, "Mabs"].min())
    # r-band second peak: brightest r detection after the plateau (the main
    # peak near MJD 60867, ~ -18.2)
    r_pk2 = float(r.loc[r["MJD"] > PLATEAU_WIN[1], "Mabs"].min())
    # dip proxy: last r detection before the rise
    pre = r[(r["MJD"] < pk_mjd) & (r["MJD"] > PLATEAU_WIN[0])]
    dip_mjd = float(pre["MJD"].max())
    t_rise = pk_mjd - dip_mjd
    # Decay target from the r/g detections (same +20..+60 d window after the
    # o peak), matching the band choice in model_metrics.
    band_decays = []
    for band in ("r", "g"):
        det, _ = load_obs_det(SNB_OBS, band)
        sel = ((det["MJD"] >= pk_mjd + DECAY_WIN[0])
               & (det["MJD"] <= pk_mjd + DECAY_WIN[1]))
        if sel.sum() >= 3:
            p = np.polyfit(det.loc[sel, "MJD"], det.loc[sel, "Mabs"], 1)
            band_decays.append(float(p[0] * 100.0))
    decay = float(np.mean(band_decays)) if band_decays else np.nan
    return dict(pk_mjd=pk_mjd, o_pk=pk_mag, r_plateau=r_plat, t_rise=t_rise,
                r_pk1=r_pk1, r_pk2=r_pk2, decay=decay, dip_mjd=dip_mjd)


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    tg = obs_targets()
    print("OBS targets:")
    print(f"  r plateau (MJD {PLATEAU_WIN[0]:.0f}-{PLATEAU_WIN[1]:.0f}): "
          f"{tg['r_plateau']:.2f}")
    print(f"  o peak: {tg['o_pk']:.2f} at MJD {tg['pk_mjd']:.1f}")
    print(f"  rise from last plateau r det (MJD {tg['dip_mjd']:.1f}) to o "
          f"peak: {tg['t_rise']:.1f} d")
    print(f"  o decay +{DECAY_WIN[0]:.0f}..+{DECAY_WIN[1]:.0f} d: "
          f"{tg['decay']:+.2f} mag/100d\n")

    rows = []
    for mdir in sorted(glob.glob(os.path.join(GRID, "*"))):
        if not os.path.isdir(mdir):
            continue
        name = os.path.basename(mdir)
        if name.startswith(("plots", "recon")):
            continue
        prm = read_params(mdir)
        df, _ = load_model(mdir)
        if prm is None or df is None or "o" not in df or "r" not in df:
            continue
        if df["o"].notna().sum() < 10:
            continue
        m = model_metrics(df, tg["pk_mjd"])
        rows.append(dict(model=name, **prm, **m))
    tab = pd.DataFrame(rows)
    tab["d_plateau"] = tab["r_plateau"] - tg["r_plateau"]
    tab["d_peak"] = tab["o_pk"] - tg["o_pk"]
    tab["d_rise"] = tab["t_rise"] - tg["t_rise"]
    tab["d_decay"] = tab["decay"] - tg["decay"]
    tab = tab.sort_values(["M", "ts", "mdot"]).reset_index(drop=True)
    csv = os.path.join(OUTDIR, "sn_b_trend_metrics.csv")
    tab.to_csv(csv, index=False, float_format="%.3f")
    print(f"saved {csv}  ({len(tab)} models)\n")

    cols = ["model", "M", "v", "nrho", "ts", "mdot", "vw", "r_plateau",
            "dip", "t_pk", "t_rise", "o_pk", "decay", "decay_win", "cov"]
    show = tab[cols].copy()
    show["mdot"] = show["mdot"].map(lambda x: f"{x:.0e}")
    with pd.option_context("display.width", 220, "display.max_rows", 200,
                           "display.float_format", "{:.2f}".format):
        print(show.to_string(index=False))

    # ---- trend figure: metrics vs t_start, coloured by mdot, marker by M
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    mk = {5.0: "o", 7.0: "s", 8.0: "D", 10.0: "^", 12.0: "v"}
    cm = {0.0001: "tab:purple", 0.0003: "tab:blue", 0.001: "tab:green",
          0.002: "tab:olive", 0.003: "tab:red"}
    base = tab[(tab["v"] == 10000) & (tab["nrho"] == 8.0) & (tab["vw"] == 100)
               & (tab["ts"] <= 25)].copy()
    # decay slopes from <20 d of post-peak coverage are not trustworthy
    base.loc[base["cov"] < 20, "decay"] = np.nan
    panels = [("r_plateau", tg["r_plateau"], "early r plateau (mag)", True),
              ("t_rise", tg["t_rise"], "rise time dip->peak (d)", False),
              ("o_pk", tg["o_pk"], "o peak (mag)", True),
              ("decay", tg["decay"], "o decay +20..60d (mag/100d)", False)]
    for ax, (key, target, lbl, inv) in zip(axes.flat, panels):
        for (M, md), g in base.groupby(["M", "mdot"]):
            ax.plot(g["ts"], g[key], marker=mk.get(M, "x"), ls="-", lw=0.8,
                    color=cm.get(md, "k"), ms=6, alpha=0.85,
                    label=f"M={M:g} mdot={md:g}")
        ax.axhline(target, color="k", ls=":", lw=1.2)
        ax.set_xlabel("t_start (d)")
        ax.set_ylabel(lbl)
        if inv:
            ax.invert_yaxis()
        ax.grid(alpha=0.3)
    axes[0, 0].legend(fontsize=6, ncol=3)
    fig.suptitle("SN_B grid trends vs t_start (v=10000, n_rho=8, v_w=100; "
                 "dotted = obs target)")
    plt.tight_layout()
    f1 = os.path.join(OUTDIR, "sn_b_trends_vs_tstart.png")
    plt.savefig(f1, dpi=140)
    print("\nsaved", f1)

    # ---- 4-criteria residual figure: every model, |delta| per criterion
    fig, ax = plt.subplots(figsize=(14, 6))
    crit = ["d_plateau", "d_peak", "d_rise", "d_decay"]
    scale = [1.0, 1.0, 10.0, 2.0]  # units per "1 mag-equivalent"
    lbls = ["r plateau (mag)", "o peak (mag)", "rise time (/10 d)",
            "decay (/2 mag/100d)"]
    t2 = tab.copy()
    t2["score"] = sum((t2[c].abs() / s).fillna(3.0)
                      for c, s in zip(crit, scale))
    t2 = t2.sort_values("score")
    x = np.arange(len(t2))
    bottom = np.zeros(len(t2))
    for c, s, lb in zip(crit, scale, lbls):
        h = (t2[c].abs() / s).fillna(3.0).values
        ax.bar(x, h, bottom=bottom, label=lb, width=0.8)
        bottom += h
    ax.set_xticks(x)
    ax.set_xticklabels(t2["model"], rotation=90, fontsize=6)
    ax.set_ylabel("sum of |residual| (mag-equivalent; NaN -> 3)")
    ax.set_title("Four-criteria residual per model (lower = better; "
                 "sorted)")
    ax.legend(fontsize=8)
    ax.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    f2 = os.path.join(OUTDIR, "sn_b_four_criteria_residuals.png")
    plt.savefig(f2, dpi=140)
    print("saved", f2)
    print("\nTop 10 by combined residual:")
    with pd.option_context("display.width", 200,
                           "display.float_format", "{:.2f}".format):
        print(t2[["model", "score"] + crit + ["cov"]].head(10)
              .to_string(index=False))


if __name__ == "__main__":
    main()
