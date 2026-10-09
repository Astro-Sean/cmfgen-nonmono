#!/usr/bin/env python3
"""Rank completed SN_B_double_peak_grid models by direct lightcurve RMS.

For each model .tt:
  - time-shift only: model o-band peak is anchored to the observed ATLAS
    o-band peak MJD (same convention as plot_sn_b_completed_vs_obs.py /
    plot_sn_b_par_vs_obs.py).  Magnitudes stay absolute.
  - inverse-variance weighted RMS between the model curve and every
    observed detection that falls inside the model's time coverage,
    per band (ATLAS o, PSF r, PSF g).
  - combined score = weighted RMS over all bands/points.

Also reports the double-peak morphology metrics from
score_sn_b_double_peak.py so the ranking can be checked against the
first-peak / dip / decay targets (obs: r 1st ~-15.6, o peak -18.54,
decay +1.97 mag/100d).
"""
import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_sn_a_pba_grid_vs_obs import BASE, DM, load_model
from plot_sn_b_grid_vs_obs import A_O, PEAK_MJD, load_atlas
from plot_sn_b_par_vs_obs import SNB_OBS, load_obs_det
from score_sn_b_double_peak import score_model

GRID = os.path.join(BASE, "SN_B_double_peak_grid")


def band_rms(model_mjd, model_mag, obs_mjd, obs_mag, obs_err):
    """Weighted RMS of model vs obs at obs epochs inside model coverage."""
    sel = (obs_mjd >= model_mjd.min()) & (obs_mjd <= model_mjd.max())
    if sel.sum() < 2:
        return np.nan, 0, np.nan
    pred = np.interp(obs_mjd[sel], model_mjd, model_mag)
    w = 1.0 / np.maximum(obs_err[sel], 0.05) ** 2
    rms = float(np.sqrt(np.average((pred - obs_mag[sel]) ** 2, weights=w)))
    bias = float(np.average(pred - obs_mag[sel], weights=w))
    return rms, int(sel.sum()), bias


def main():
    # ---- observations ------------------------------------------------
    atlas = load_atlas(
        os.path.join(BASE, "observational_data",
                     "SN2024pba_atclean_cleaned.csv"), band="o")
    o_det = atlas["det"].copy()
    o_det["MJD"] = o_det["phase"] + PEAK_MJD
    ob_peak_mjd = o_det.loc[o_det["Mabs"].idxmin(), "MJD"]
    ob_peak_mag = o_det["Mabs"].min()
    print(f"obs ATLAS o peak: {ob_peak_mag:.2f} at MJD {ob_peak_mjd:.1f}")

    r_det, _ = load_obs_det(SNB_OBS, "r")
    g_det, _ = load_obs_det(SNB_OBS, "g")
    print(f"obs det: o={len(o_det)} r={len(r_det)} g={len(g_det)}\n")

    # ---- models ------------------------------------------------------
    dirs = sorted(d for d in glob.glob(os.path.join(GRID, "*"))
                  if os.path.isdir(d)
                  and not os.path.basename(d).startswith(("plots", "recon")))
    rows = []
    for mdir in dirs:
        name = os.path.basename(mdir)
        df, src = load_model(mdir)
        if df is None or "o" not in df or df["o"].notna().sum() < 5:
            continue
        pk = df["o"].idxmin()
        t0 = ob_peak_mjd - df.loc[pk, "time"]
        mjd = t0 + df["time"].values

        res = {}
        tot_r2 = tot_w = 0.0
        n_all = 0
        for band, det in (("o", o_det), ("r", r_det), ("g", g_det)):
            if band not in df or not len(det):
                res[band] = (np.nan, 0, np.nan)
                continue
            rms, n, bias = band_rms(
                mjd, df[band].values, det["MJD"].values,
                det["Mabs"].values, det["Mabs_err"].values)
            res[band] = (rms, n, bias)
            if n:
                w = 1.0 / np.maximum(det["Mabs_err"].values, 0.05) ** 2
                sel = (det["MJD"].values >= mjd.min()) & \
                      (det["MJD"].values <= mjd.max())
                pred = np.interp(det["MJD"].values[sel], mjd,
                                 df[band].values)
                tot_r2 += float((w[sel] * (pred - det["Mabs"].values[sel])
                                 ** 2).sum())
                tot_w += float(w[sel].sum())
                n_all += int(sel.sum())
        score = np.sqrt(tot_r2 / tot_w) if tot_w else np.nan

        morph = score_model(
            os.path.join(mdir, "stella", f"{name}.tt")) \
            if os.path.exists(os.path.join(mdir, "stella", f"{name}.tt")) \
            else None
        rows.append(dict(name=name, score=score, n=n_all,
                         res=res, morph=morph,
                         pk_t=float(df.loc[pk, "time"]),
                         pk_m=float(df.loc[pk, "o"])))

    rows.sort(key=lambda r: (np.isnan(r["score"]), r["score"]))
    hdr = (f"{'model':30s} {'RMS':>5s} {'n':>3s} "
           f"{'o_rms':>6s} {'r_rms':>6s} {'g_rms':>6s} "
           f"{'m1st':>7s} {'dip':>6s} {'mPk':>7s} {'dec':>5s}")
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        m = r["morph"] or {}
        mf = f"{m['m_first']:.2f}" if m.get("m_first") else "  --"
        dp = f"{m['dip_mag']:.2f}" if m.get("dip_mag") else "  --"
        mp = f"{m['m_peak']:.2f}" if m.get("m_peak") else "  --"
        dc = f"{m['decay']:.2f}" if m.get("decay") else " --"
        o_rms = f"{r['res']['o'][0]:.2f}" if r["res"]["o"][1] else "  --"
        r_rms = f"{r['res']['r'][0]:.2f}" if r["res"]["r"][1] else "  --"
        g_rms = f"{r['res']['g'][0]:.2f}" if r["res"]["g"][1] else "  --"
        print(f"{r['name']:30s} {r['score']:5.2f} {r['n']:3d} "
              f"{o_rms:>6s} {r_rms:>6s} {g_rms:>6s} "
              f"{mf:>7s} {dp:>6s} {mp:>7s} {dc:>5s}")


if __name__ == "__main__":
    main()
