#!/usr/bin/env python3
"""Plot SN_B grid trends: one figure per parameter lever.

Each figure shows the model family (one parameter varied, rest fixed at
the SB_M10_ts8_W3e-4 winner values) against the SN 2024pba observations,
with the four fit criteria annotated:

  1. early r plateau target  -15.59  (red dotted band)
  2. rise plateau-end -> o peak      (annotated arrow on obs)
  3. o peak target           -18.54  (orange dotted line + vline)
  4. decay target            +2.3 mag/100d (orange dotted guide line)

Convention (same as plot_sn_b_best_vs_obs.py): model o-peaks are
time-shifted to the observed o peak; magnitudes stay absolute.
o-band solid, r-band dashed (same colour, thinner).

Outputs to SN_B_double_peak_grid/plots/trends/.
"""
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

GRID = os.path.join(BASE, "SN_B_double_peak_grid")
OUTDIR = os.path.join(GRID, "plots", "trends")

R_PLATEAU_TARGET = -15.59
O_PEAK_TARGET = -18.54
DECAY_TARGET = 2.3  # mag/100d

# model dir name -> legend label
FAMILIES = [
    ("tstart",
     "t_start (M=10, mdot=3e-4): later launch -> brighter plateau & peak",
     [("SB_M10_ts8_W3e-4", "ts8"),
      ("SB_M10_ts12_W3e-4", "ts12"),
      ("SB_M10_ts16_W3e-4", "ts16"),
      ("SB_M10_ts20_W3e-4", "ts20")]),
    ("mdot",
     "wind Mdot (M=10, ts8): denser wind -> brighter single early peak",
     [("SB_M10_ts8_W1e-4", "mdot=1e-4"),
      ("SB_M10_ts8_W3e-4", "mdot=3e-4"),
      ("SB_M10_ts8_W1e-3", "mdot=1e-3"),
      ("SB_M10_ts8_W3e-3", "mdot=3e-3")]),
    ("mass",
     "shell mass (ts8, mdot=3e-4): M barely moves plateau or peak",
     [("SB_M5_ts8_W3e-4", "M=5"),
      ("SB_M8_ts8_W3e-4", "M=8"),
      ("SB_M10_ts8_W3e-4", "M=10")]),
    ("vwind",
     "wind velocity (M=10, ts8, mdot=3e-4): faster wind -> thinner "
     "density D=Mdot/v_w",
     [("SB_M10_ts8_W3e-4", "v_w=100"),
      ("SB_M10_ts8_vw300", "v_w=300"),
      ("SB_M10_ts8_vw850_t180", "v_w=850 (obs)")]),
    ("nrho",
     "envelope slope (M=10, ts8, mdot=3e-4): n_rho=12 merges the peaks",
     [("SB_M10_ts8_W3e-4", "n_rho=8"),
      ("SB_M10_ts8_nr12", "n_rho=12")]),
    ("corridor",
     "wind density D=mdot/v_w sequence (single-parameter family): "
     "morphology transition",
     [("SB_M10_ts8_W1e-4", "D=1e-6 (mdot=1e-4)"),
      ("SB_M10_ts8_vw300", "D=1e-6 (v_w=300)"),
      ("SB_M10_ts8_W3e-4", "D=3e-6 (parent)"),
      ("SB_M10_ts8_W1e-3", "D=1e-5 (mdot=1e-3)"),
      ("SB_M10_ts8_W3e-3", "D=3e-5 (mdot=3e-3)")]),
]

COLORS = ["tab:blue", "tab:orange", "tab:green", "tab:red",
          "tab:purple", "tab:brown"]


def plot_family(fig_name, title, models, obs):
    o_det, r_det, g_det, ob_peak_mjd = obs
    fig, ax = plt.subplots(figsize=(10.5, 6.5))
    any_model = False
    for i, (dname, lbl) in enumerate(models):
        df, _ = load_model(os.path.join(GRID, dname))
        if df is None or "o" not in df or df["o"].notna().sum() < 5:
            print(f"  {dname}: no .tt yet (skipped)")
            continue
        any_model = True
        t = df["time"].values.astype(float)
        o = df["o"].values.astype(float)
        ipk = int(np.nanargmin(o))
        t0 = ob_peak_mjd - t[ipk]
        mjd = t0 + t
        col = COLORS[i % len(COLORS)]
        ax.plot(mjd, o, color=col, lw=2.0, zorder=5,
                label=f"{lbl} (o pk {o[ipk]:.2f})")
        if "r" in df:
            ax.plot(mjd, df["r"].values, color=col, lw=1.0, ls="--",
                    alpha=0.5, zorder=4)
    if not any_model:
        plt.close(fig)
        return None

    ax.errorbar(o_det["MJD"], o_det["Mabs"], yerr=o_det["Mabs_err"],
                fmt="^", ms=4.5, capsize=2, lw=0.7, color="darkorange",
                label="obs o", zorder=10)
    ax.errorbar(r_det["MJD"], r_det["Mabs"], yerr=r_det["Mabs_err"],
                fmt="o", ms=3.5, capsize=2, lw=0.7, color="tab:red",
                mfc="none", label="obs r", zorder=10)
    ax.errorbar(g_det["MJD"], g_det["Mabs"], yerr=g_det["Mabs_err"],
                fmt="s", ms=3.5, capsize=2, lw=0.7, color="tab:green",
                mfc="none", label="obs g", zorder=10)

    # criterion 1: plateau target
    ax.axhspan(R_PLATEAU_TARGET - 0.1, R_PLATEAU_TARGET + 0.1,
               color="tab:red", alpha=0.12, zorder=1)
    ax.annotate("r plateau target -15.6", xy=(60792, R_PLATEAU_TARGET - 0.28),
                fontsize=8, color="tab:red")
    # criterion 3: o peak target
    ax.axhline(O_PEAK_TARGET, color="darkorange", ls=":", lw=1, alpha=0.7)
    ax.axvline(ob_peak_mjd, color="darkorange", ls=":", lw=1, alpha=0.5)
    ax.annotate(f"o peak {O_PEAK_TARGET:.2f}", xy=(60792, O_PEAK_TARGET + 0.08),
                fontsize=8, color="darkorange")
    # criterion 2: rise arrow (obs plateau end -> o peak)
    ax.annotate("", xy=(ob_peak_mjd, O_PEAK_TARGET),
                xytext=(60831.5, R_PLATEAU_TARGET),
                arrowprops=dict(arrowstyle="->", color="k", lw=1.2))
    ax.annotate("rise ~19 d", xy=(60834, -16.6), fontsize=8, color="k")
    # criterion 4: decay guide (+2.3 mag/100d from peak)
    xg = np.array([ob_peak_mjd, ob_peak_mjd + 60])
    ax.plot(xg, O_PEAK_TARGET + DECAY_TARGET * (xg - ob_peak_mjd) / 100.0,
            color="darkorange", ls="-.", lw=1, alpha=0.6)
    ax.annotate("decay +2.3 mag/100d", xy=(ob_peak_mjd + 18, -18.0),
                fontsize=8, color="darkorange", rotation=4)

    ax.set_ylim(-13.5, -20)
    ax.set_xlim(60790, 60930)
    ax.set_xlabel("MJD (model o-peaks aligned to obs o peak)")
    ax.set_ylabel("Absolute magnitude (AB)")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, ncol=2, loc="lower left")
    ax.set_title(title, fontsize=10)
    plt.tight_layout()
    out = os.path.join(OUTDIR, f"trend_{fig_name}.png")
    plt.savefig(out, dpi=140)
    plt.close(fig)
    return out


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    atlas = load_atlas(os.path.join(BASE, "observational_data",
                                    "SN2024pba_atclean_cleaned.csv"), band="o")
    o_det = atlas["det"].copy()
    o_det["MJD"] = o_det["phase"] + PEAK_MJD
    ob_peak_mjd = float(o_det.loc[o_det["Mabs"].idxmin(), "MJD"])
    r_det, _ = load_obs_det(SNB_OBS, "r")
    g_det, _ = load_obs_det(SNB_OBS, "g")
    obs = (o_det, r_det, g_det, ob_peak_mjd)

    for fig_name, title, models in FAMILIES:
        out = plot_family(fig_name, title, models, obs)
        if out:
            print("saved", out)


if __name__ == "__main__":
    main()
