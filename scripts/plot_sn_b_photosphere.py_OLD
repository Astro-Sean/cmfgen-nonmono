#!/usr/bin/env python3
"""Photospheric radius vs time for an SN_B model.

Uses STELLA's own output columns from the .tt file --
time, Tbb, rbb, Teff, Rlast_sc, R(tau2/3), Mbol -- nothing is
recomputed from the .swd structure.

Three stacked panels on a shared time axis:
  1. lightcurve (o band + Mbol)
  2. R(tau2/3) photospheric radius, with the initial structure
     boundaries (core/envelope/buffer/wind/ejecta) marked, plus rbb
     and Rlast_sc for reference
  3. Tbb and Teff

Usage:
  python scripts/plot_sn_b_photosphere.py [MODEL_DIR]
  (default: SN_B_double_peak_grid/SB_M10_ts8_W3e-4)
"""
import glob
import os
import sys

import matplotlib
import numpy as np
import yaml

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_sn_a_pba_grid_vs_obs import BASE  # noqa: E402

KM = 1e5
DAY = 86400.0

C_CORE = "#4477AA"
C_ENV = "#88CCEE"
C_BUF = "#DDCC77"
C_WIND = "#CC6677"
C_EJ = "#AA4499"


def load_tt(path):
    """STELLA .tt -> dict of named columns (time in days)."""
    rows = []
    for line in open(path):
        s = line.strip()
        if not s or not (s[0].isdigit() or s[0] in "-+"):
            continue
        try:
            row = [float(v) for v in s.split()]
        except ValueError:
            continue
        if len(row) >= 22:
            rows.append(row)
    arr = np.array(rows)
    cols = ["time", "Tbb", "rbb", "Teff", "Rlast_sc", "R_ph", "Mbol",
            "MU", "MB", "MV", "MI", "MR", "Mu", "Mg", "Mr", "Mi", "Mz",
            "MJ", "MH", "MKs", "Mc", "Mo"]
    return {c: arr[:, i] for i, c in enumerate(cols)}


def boundaries(mdir):
    """Initial structure boundaries from config + input_profile (cm)."""
    x = yaml.safe_load(open(os.path.join(mdir, "config.yaml"))
                       )["x_model_from_sn"]
    a = np.loadtxt(os.path.join(mdir, "stella", "input_profile"),
                   comments="#")
    r, v = a[:, 0], a[:, 3]
    r_out = x["shell_velocity_kms"] * KM * x["t_start_days"] * DAY
    r0 = x["v0_kms"] * KM * x["t_start_days"] * DAY
    vw = x["wind_velocity_kms"] * KM
    iw = np.where((r > r_out) & (np.abs(v - vw) < 0.10 * vw))[0]
    return dict(r_in=r[0], r0=r0, r_out=r_out,
                r_w0=r[iw[0]], r_ej=r[iw[-1]], r_end=r[-1])


def main():
    mdir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        BASE, "SN_B_double_peak_grid", "SB_M10_ts8_W3e-4")
    if not os.path.isabs(mdir):
        mdir = os.path.join(BASE, mdir)
    name = os.path.basename(mdir)

    tt = glob.glob(os.path.join(mdir, "stella", "*.tt"))
    if not tt:
        sys.exit(f"{name}: no .tt")
    d = load_tt(tt[0])
    b = boundaries(mdir)

    fig, (ax_l, ax_r, ax_t) = plt.subplots(
        3, 1, figsize=(8.5, 8.0), sharex=True,
        gridspec_kw=dict(hspace=0.09))

    # ---- 1: lightcurve (STELLA o band + Mbol)
    ax_l.plot(d["time"], d["Mo"], color="k", lw=1.4, label="o band")
    ax_l.plot(d["time"], d["Mbol"], color="tab:purple", lw=1.1,
              ls="--", alpha=0.85, label="Mbol")
    ax_l.plot(d["time"], d["Mr"], color="tab:red", lw=1.0, ls=":",
              alpha=0.8, label="r band")
    ax_l.set_ylim(-14, -20)
    ax_l.set_ylabel("abs. mag")
    ax_l.legend(fontsize=8, loc="lower left", frameon=False, ncol=3)
    ax_l.set_title(f"{name}: STELLA photospheric evolution", fontsize=10)

    # ---- 2: photospheric radius with initial component bands
    bands = [(b["r_in"], b["r0"], C_CORE, "core"),
             (b["r0"], b["r_out"], C_ENV, "envelope"),
             (b["r_out"], b["r_w0"], C_BUF, "buffer"),
             (b["r_w0"], b["r_ej"], C_WIND, "wind"),
             (b["r_ej"], b["r_end"], C_EJ, "SN$_A$ ejecta")]
    for ra, rb, col, lab in bands:
        ax_r.axhspan(ra, rb, color=col, alpha=0.25, lw=0)
        ax_r.text(0.012, np.sqrt(ra * rb), lab, fontsize=7,
                  va="center", transform=ax_r.get_yaxis_transform())
    ax_r.plot(d["time"], d["R_ph"], color="k", lw=1.4,
              label="$R(\\tau=2/3)$")
    ax_r.plot(d["time"], d["rbb"], color="tab:blue", lw=1.0, ls="--",
              alpha=0.8, label="$R_{\\rm bb}$")
    ax_r.plot(d["time"], d["Rlast_sc"], color="tab:green", lw=1.0,
              ls=":", alpha=0.8, label="$R_{\\rm last\\,sc}$")
    ax_r.set_yscale("log")
    r_all = np.concatenate([d["R_ph"], d["rbb"], d["Rlast_sc"]])
    r_all = r_all[np.isfinite(r_all) & (r_all > 0)]
    ax_r.set_ylim(r_all.min() / 1.8, r_all.max() * 1.8)
    ax_r.set_ylabel("$R$ [cm]")
    ax_r.legend(loc="lower right", fontsize=8, frameon=False)

    # ---- 3: temperatures
    ax_t.plot(d["time"], d["Tbb"], color="k", lw=1.4,
              label="$T_{\\rm bb}$")
    ax_t.plot(d["time"], d["Teff"], color="tab:orange", lw=1.1,
              ls="--", alpha=0.85, label="$T_{\\rm eff}$")
    ax_t.set_yscale("log")
    ax_t.set_ylabel("$T$ [K]")
    ax_t.set_xlabel("time since shell start [d]")
    ax_t.legend(loc="upper right", fontsize=8, frameon=False)

    outdir = os.path.join(mdir, "plots")
    os.makedirs(outdir, exist_ok=True)
    out = os.path.join(outdir, f"{name}_photosphere.png")
    fig.savefig(out, dpi=160, bbox_inches="tight")
    print("saved", out)


if __name__ == "__main__":
    main()
