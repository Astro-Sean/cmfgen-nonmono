#!/usr/bin/env python3
"""Per-model two-event comparison plots for SN 2024pba.

One physical timeline per model, one free parameter: dt, the time
offset in days applied to the model, bounded to +/-10 d about the
alignment of the model's luminosity peak with the observed peak of
the event it is meant to reproduce.  Magnitudes are never rescaled.

Outputs into <model_dir>/plots/:

  <name>_vs_sn_a.png -- days since the SN_A explosion (MJD 60488),
      zoomed on the SN_A window (0-80 d).  SN_A g/r/i detections and
      upper limits; the model g/r/i curve on the same absolute axis.
      For cocoon / x_model_from_sn models the model t=0 is contact
      with the SN_A ejecta, placed at SN_A age = source_age_days +
      t_start_days; other models are treated as SN_A explosions
      (t=0 at the SN_A epoch).

  <name>_vs_sn_b.png -- the SN_B rebrightening window (~300-470 d
      after SN_A, i.e. days since the observed ATLAS-o peak shown on a
      secondary convention).  Model o/r vs ATLAS o, PSF g/r detections,
      with the r-plateau / o-peak / decay guides.

  <name>_structure.png -- the .hyd structure with the three physical
      components shaded: shell (SN_B), wind/CSM, SN_A ejecta.

Usage:
  python scripts/plot_model_vs_sn_ab.py MODEL_DIR [MODEL_DIR ...]
"""
import os
import sys

import matplotlib
import numpy as np
import yaml

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_sn_a_pba_grid_vs_obs import (  # noqa: E402
    BASE,
    EXP_MJD,
    FIT_MAX_MJD,
    OBS_FILE,
    load_model,
    load_observations,
)
from plot_sn_b_grid_vs_obs import PEAK_MJD, load_atlas  # noqa: E402
from plot_sn_b_par_vs_obs import SNB_OBS, load_obs_det  # noqa: E402

SN_A_COLORS = {"g": "tab:green", "r": "tab:red", "i": "goldenrod"}
R_PLATEAU_TARGET = -15.59
O_PEAK_TARGET = -18.54
DECAY_TARGET = 2.3  # mag/100 d
DT_MAX = 10.0


def sn_a_offset_days(mdir):
    """SN_A age (days) at the model's t=0 (cocoon contact for
    x_model_from_sn models; 0 for SN_A explosions)."""
    cfg_path = os.path.join(mdir, "config.yaml")
    if not os.path.exists(cfg_path):
        return 0.0
    try:
        with open(cfg_path) as fh:
            cfg = yaml.safe_load(fh)
    except Exception:
        return 0.0
    x = (cfg or {}).get("x_model_from_sn") or {}
    if not x:
        return 0.0
    return float(x.get("source_age_days", 0.0)) + float(
        x.get("t_start_days", 0.0)
    )


def load_snb_obs():
    """SN_B obs on the absolute MJD axis: ATLAS o, PSF r, PSF g."""
    atlas = load_atlas(
        os.path.join(BASE, "observational_data",
                     "SN2024pba_atclean_cleaned.csv"),
        band="o",
    )
    o_det = atlas["det"].copy()
    o_det["MJD"] = o_det["phase"] + PEAK_MJD
    r_det, _ = load_obs_det(SNB_OBS, "r")
    g_det, _ = load_obs_det(SNB_OBS, "g")
    ob_peak_mjd = float(o_det.loc[o_det["Mabs"].idxmin(), "MJD"])
    return o_det, r_det, g_det, ob_peak_mjd


def fit_dt(df, offset, snb_obs, sna_obs):
    """The single free parameter: dt in [-10, +10] d.

    If the model is a cocoon / SN_B model (offset > 0), align the model
    o-peak to the observed ATLAS-o peak and pick dt by RMS over the o
    detections within +/-30 d of peak.  Otherwise align the model i-peak
    to the SN_A i-peak and fit dt by RMS over the SN_A i detections.
    """
    t = df["time"].values + offset
    if offset > 0:
        o_det = snb_obs[0]
        obs_day = o_det["MJD"].values - EXP_MJD
        obs_mag = o_det["Mabs"].values
        mod_mag = df["o"].values if "o" in df else df["r"].values
        obs_peak_day = obs_day[np.argmin(obs_mag)]
        dt0 = obs_peak_day - t[np.nanargmin(mod_mag)]
        band = "o" if "o" in df else "r"
        sel = np.abs(obs_day - obs_peak_day) < 30
    else:
        i_det = sna_obs["i"]["det"]
        obs_day = i_det["MJD"].values - EXP_MJD
        obs_mag = i_det["Mabs"].values
        mod_mag = df["i"].values
        obs_peak_day = obs_day[np.argmin(obs_mag)]
        dt0 = obs_peak_day - t[np.nanargmin(mod_mag)]
        band = "i"
        sel = (obs_day < FIT_MAX_MJD - EXP_MJD)

    def rms(dt):
        interp = np.interp(obs_day[sel], t + dt, df[band].values)
        return float(np.sqrt(np.mean((interp - obs_mag[sel]) ** 2)))

    dts = np.arange(-DT_MAX, DT_MAX + 0.25, 0.25)
    best = min(dts, key=lambda d: rms(dt0 + d))
    return dt0 + best


def plot_vs_sn_a(mdir, df, x_abs, outdir):
    name = os.path.basename(mdir)
    obs = load_observations(OBS_FILE)
    fig, ax = plt.subplots(figsize=(10.5, 6.5))

    for band in ("g", "r", "i"):
        det, nd = obs[band]["det"], obs[band]["nd"]
        c = SN_A_COLORS[band]
        ax.errorbar(det["MJD"] - EXP_MJD, det["Mabs"], yerr=det["Mabs_err"],
                    fmt="o", ms=4, capsize=2, lw=0.7, mfc=c, mec="k",
                    ecolor=c, label=f"obs {band}", zorder=10)
        if len(nd):
            ax.errorbar(nd["MJD"] - EXP_MJD, nd["Mabs_lim"], fmt="v",
                        ms=3.5, mfc="none", mec=c, ecolor=c, alpha=0.5,
                        zorder=9)
    for band, ls in (("g", ":"), ("r", "-"), ("i", "--")):
        if band in df:
            ax.plot(x_abs, df[band].values, ls, lw=1.8, color="tab:blue",
                    zorder=6)
    ax.plot([], [], "-", color="tab:blue", lw=1.8, label=f"{name} g/r/i")

    ax.axvline(FIT_MAX_MJD - EXP_MJD, color="k", ls=":", lw=0.8,
               alpha=0.4)
    ax.annotate("SN_A fit window ends", xy=(FIT_MAX_MJD - EXP_MJD + 2,
                -13.0), fontsize=8, alpha=0.7)
    ax.set_ylim(-12.5, -19.0)
    ax.set_xlim(-10, 90)
    ax.set_xlabel("days since SN_A explosion (MJD 60488)")
    ax.set_ylabel("Absolute magnitude (AB)")
    ax.set_title(f"{name} vs SN_A (dt fixed by model-peak alignment)")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, ncol=2, loc="lower left")
    plt.tight_layout()
    out = os.path.join(outdir, f"{name}_vs_sn_a.png")
    plt.savefig(out, dpi=140)
    plt.close(fig)
    return out


def plot_vs_sn_b(mdir, df, x_abs, snb_obs, outdir):
    name = os.path.basename(mdir)
    o_det, r_det, g_det, ob_peak_mjd = snb_obs
    fig, ax = plt.subplots(figsize=(10.5, 6.5))

    if "o" in df:
        ax.plot(x_abs, df["o"].values, "-", lw=2.0, color="tab:red",
                label=f"{name} o", zorder=6)
    if "r" in df:
        ax.plot(x_abs, df["r"].values, "--", lw=1.5, color="tab:red",
                alpha=0.7, label=f"{name} r", zorder=6)
    ax.errorbar(o_det["MJD"] - EXP_MJD, o_det["Mabs"],
                yerr=o_det["Mabs_err"], fmt="^", ms=4.5, capsize=2,
                lw=0.7, mfc="orange", mec="k", ecolor="orange",
                label="obs o", zorder=10)
    ax.errorbar(r_det["MJD"] - EXP_MJD, r_det["Mabs"],
                yerr=r_det["Mabs_err"], fmt="o", ms=3.5, capsize=2,
                lw=0.7, mfc="tab:red", mec="k", ecolor="tab:red",
                label="obs r", zorder=10)
    ax.errorbar(g_det["MJD"] - EXP_MJD, g_det["Mabs"],
                yerr=g_det["Mabs_err"], fmt="s", ms=3.5, capsize=2,
                lw=0.7, mfc="tab:green", mec="k", ecolor="tab:green",
                label="obs g", zorder=10)

    ax.axhspan(R_PLATEAU_TARGET - 0.1, R_PLATEAU_TARGET + 0.1,
               color="tab:red", alpha=0.12, zorder=1)
    ax.axhline(O_PEAK_TARGET, color="darkorange", ls=":", lw=1, alpha=0.7)
    ax.axvline(PEAK_MJD - EXP_MJD, color="darkorange", ls=":", lw=1,
               alpha=0.5)
    pk = PEAK_MJD - EXP_MJD
    xg = np.array([pk, pk + 60.0])
    ax.plot(xg, O_PEAK_TARGET + DECAY_TARGET * (xg - pk) / 100.0,
            color="darkorange", ls="-.", lw=1, alpha=0.6)
    ax.annotate("r plateau -15.6", xy=(300, R_PLATEAU_TARGET - 0.35),
                fontsize=8, color="tab:red")
    ax.annotate(f"o peak {O_PEAK_TARGET:.2f}", xy=(300, O_PEAK_TARGET + 0.08),
                fontsize=8, color="darkorange")
    ax.annotate("decay +2.3 mag/100d", xy=(pk + 12, -18.1), fontsize=8,
                color="darkorange", rotation=4)

    ax.set_ylim(-14.5, -19.0)
    ax.set_xlim(280, 500)
    ax.set_xlabel("days since SN_A explosion (MJD 60488)")
    ax.set_ylabel("Absolute magnitude (AB)")
    ax.set_title(f"{name} vs SN_B (dt fixed by model-peak alignment)")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, ncol=2, loc="lower left")
    plt.tight_layout()
    out = os.path.join(outdir, f"{name}_vs_sn_b.png")
    plt.savefig(out, dpi=140)
    plt.close(fig)
    return out


def plot_structure(mdir, outdir):
    """Radius-density-velocity structure with the physical components
    shaded: shell (SN_B), wind/CSM, SN_A ejecta."""
    name = os.path.basename(mdir)
    candidates = [
        os.path.join(mdir, "stella", "modmake", f"{name}.hyd"),
        os.path.join(mdir, "stella", "modmake", f"{name}_mod.hyd"),
    ]
    hyd = next((c for c in candidates if os.path.exists(c)), None)
    if hyd is None:
        print(f"  {name}: no .hyd for structure plot")
        return None
    lines = open(hyd).read().splitlines()
    nz = int(lines[0].split()[1])
    arr = np.array([[float(x) for x in l.split()]
                    for l in lines[1:1 + nz] if len(l.split()) >= 6])
    r, rho, v = arr[:, 2], arr[:, 3], arr[:, 5]
    with np.errstate(divide="ignore", invalid="ignore"):
        tv = r / v / 86400.0  # homologous age per zone

    fig, ax = plt.subplots(figsize=(10.5, 6.5))
    ax.loglog(r, rho, color="k", lw=1.2, label="rho")
    ax2 = ax.twinx()
    ax2.semilogx(r, v / 1e5, color="tab:purple", lw=1.0, alpha=0.7,
                 label="v")

    # Component boundaries from the velocity structure:
    #   shell:  outer part with r/v ~ t_start (< ~150 d)
    #   ejecta: zones with r/v in the SN_A-age range (~200-500 d)
    #   wind/CSM: the 850 km/s detached component at large r
    is_shell = tv < 150
    is_csm = np.abs(v / 1e5 - 850.0) < 1.0
    is_ejecta = ~is_shell & ~is_csm
    for mask, c, lab in ((is_shell, "tab:red", "SN_B shell"),
                       (is_ejecta, "tab:blue", "SN_A ejecta"),
                       (is_csm, "tab:green", "wind/CSM 850 km/s")):
        if mask.any():
            ax.axvspan(r[mask].min(), r[mask].max(), color=c, alpha=0.08)
            ax.annotate(lab, xy=(np.sqrt(r[mask].min() * r[mask].max()),
                        0.9), xycoords=("data", "axes fraction"),
                        fontsize=8, color=c, ha="center")
    ax.set_xlabel("radius (cm)")
    ax.set_ylabel("density (g/cc)")
    ax2.set_ylabel("velocity (km/s)", color="tab:purple")
    ax.set_title(f"{name} structure: shell / wind / ejecta")
    ax.grid(alpha=0.3, which="both")
    plt.tight_layout()
    out = os.path.join(outdir, f"{name}_structure.png")
    plt.savefig(out, dpi=140)
    plt.close(fig)
    return out


def process(mdir, snb_obs, sna_obs):
    df, _ = load_model(mdir)
    name = os.path.basename(mdir)
    outdir = os.path.join(mdir, "plots")
    os.makedirs(outdir, exist_ok=True)
    outs = []
    s = plot_structure(mdir, outdir)
    if s:
        outs.append(s)
    if df is None or not len(df):
        print(f"  {name}: no usable .tt (structure only)")
        return outs
    offset = sn_a_offset_days(mdir)
    dt = fit_dt(df, offset, snb_obs, sna_obs)
    x_abs = df["time"].values + offset + dt
    outs.append(plot_vs_sn_a(mdir, df, x_abs, outdir))
    outs.append(plot_vs_sn_b(mdir, df, x_abs, snb_obs, outdir))
    return outs


def main():
    mdirs = [a for a in sys.argv[1:] if os.path.isdir(a)]
    if not mdirs:
        print("usage: plot_model_vs_sn_ab.py MODEL_DIR [MODEL_DIR ...]")
        sys.exit(1)
    snb_obs = load_snb_obs()
    sna_obs = load_observations(OBS_FILE)
    for mdir in mdirs:
        for out in process(mdir, snb_obs, sna_obs):
            print("saved", out)


if __name__ == "__main__":
    main()
