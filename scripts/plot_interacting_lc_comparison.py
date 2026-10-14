#!/usr/bin/env python3
"""Compare the SN 2024pba lightcurve to interacting-SNe analogues.

Plots absolute magnitude vs rest-frame days since peak for SN 2024pba
(r band, heavy) and a set of interacting transients that show either a
post-peak break or interaction-driven structure: SNe 1994W, 1998S,
2005ip, 2006jd, 2010jl, iPTF13z and 2009ip.

Photometry is taken from the Open Supernova Catalog dumps in
observational_data/*_osc.csv (columns mjd,band,mag,e_mag).  Distances
are the OSC luminosity distances (or host values at low z); phases are
rest frame relative to each object's own OSC maxdate, except SN 2024pba
which uses the measured r-band peak MJD 60867.4.
"""

import argparse
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import AutoMinorLocator
import numpy as np
import pandas as pd

OBS_DIR = "/viper/ptmp1/sbrennan/SN2024pba/observational_data"

# (name, band, distance Mpc, redshift, peak MJD or None->min-mag,
#  marker colour, zorder)
ANALOGS = [
    ("SN1994W", "V", 26.0, 0.00403, 49576.0, "tab:purple", 5),
    ("SN1998S", "R", 15.0, 0.00302, None, "tab:orange", 5),
    ("SN2005ip", "r", 35.0, 0.007138, None, "tab:olive", 4),
    ("SN2006jd", "r", 79.0, 0.01855, None, "tab:cyan", 4),
    ("iPTF13z", "R", 148.8, 0.0328, None, "tab:green", 4),
    ("SN2010jl", "R", 49.0, 0.0107, 55477.5, "tab:blue", 4),
    ("SN2009ip", "r", 26.426, 0.005944, 56204.5, "tab:pink", 3),
]

# SN 2024pba
PBA_PEAK_MJD = 60867.4  # measured r-band maximum
PBA_DL_MPC = 80.4
PBA_Z = 0.0186
PBA_EBV = 0.03
A_R = 2.695 * PBA_EBV

BREAK_PHASE = 80.0

plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["DejaVu Serif"],
        "mathtext.fontset": "dejavuserif",
        "font.size": 8,
        "axes.labelsize": 8,
        "legend.fontsize": 6,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "axes.linewidth": 0.7,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


def load_osc(path, band):
    """Return (mjd, mag, emag) detections in the requested band."""
    d = pd.read_csv(path)
    sub = d[d["band"] == band].copy()
    sub = sub[np.isfinite(sub["mag"])]
    return (
        sub["mjd"].to_numpy(float),
        sub["mag"].to_numpy(float),
        pd.to_numeric(sub["e_mag"], errors="coerce").to_numpy(),
    )


def load_analog(name, band):
    """Per-file adapter; returns (mjd, mag, emag, flag).

    SN2010jl.csv stores days-since-max in its 'mjd' column; SN2009ip.csv
    is a wide multi-band table of absolute magnitudes.  Everything else
    is the OSC mjd,band,mag convention.  flag is True for 'already
    phase', 'abs' for 'already absolute magnitude', False otherwise.
    """
    if name == "SN2010jl":
        path = os.path.join(OBS_DIR, "SN2010jl.csv")
        d = pd.read_csv(path)
        # 'mjd' column is days relative to an unspecified epoch;
        # the caller re-zeros at maximum light
        return (
            d["mjd"].to_numpy(float),
            d["R"].to_numpy(float),
            d["R_err"].to_numpy(float),
            True,
        )
    if name == "SN2009ip":
        path = os.path.join(OBS_DIR, "SN2009ip.csv")
        d = pd.read_csv(path)
        sub = d[["mjd", band]].dropna()
        emag = d.loc[sub.index, f"{band}_err"].to_numpy(float)
        # this file already stores absolute magnitudes
        return (
            sub["mjd"].to_numpy(float),
            sub[band].to_numpy(float),
            emag,
            "abs",
        )
    mjd, mag, emag = load_osc(os.path.join(OBS_DIR, f"{name}_osc.csv"),
                            band)
    return mjd, mag, emag, False


def bin_lightcurve(phase, mag, emag=None, bin_days=1.0,
                   clip_sigma=3.0, spike_mag=0.75):
    """1-day binning with two-pass outlier rejection.

    Pass 1: within each bin, reject points beyond clip_sigma of the
    median (scaled MAD, floored at 0.05 mag).  Pass 2: reject isolated
    single-bin spikes -- bins deviating from both neighbours while the
    neighbours themselves agree.  Real slope breaks (e.g. SN 1994W)
    survive because consecutive bins trend together.
    """
    df = pd.DataFrame({"phase": phase, "mag": mag})
    df["emag"] = np.nan if emag is None else emag
    df["bin"] = np.floor(df["phase"] / bin_days)
    rows = []
    for _, g in df.groupby("bin"):
        m = g["mag"].to_numpy()
        if len(m) >= 4:
            med = np.median(m)
            sig = max(np.median(np.abs(m - med)) * 1.4826, 0.05)
            g = g[np.abs(m - med) <= max(clip_sigma * sig, 0.3)]
        e = np.where(np.isfinite(g["emag"]) & (g["emag"] > 0),
                     g["emag"], 0.2)
        w = 1.0 / e**2
        rows.append(
            {
                "phase": np.average(g["phase"], weights=w),
                "mag": np.average(g["mag"], weights=w),
                "emag": 1.0 / np.sqrt(w.sum()),
                "n": len(g),
            }
        )
    out = pd.DataFrame(rows).sort_values("phase").reset_index(drop=True)
    if len(out) >= 5:
        m = out["mag"].to_numpy()
        drop = np.zeros(len(out), dtype=bool)
        for i in range(1, len(out) - 1):
            step_in = abs(m[i] - m[i - 1])
            step_out = abs(m[i] - m[i + 1])
            back = abs(m[i + 1] - m[i - 1])
            if (step_in > spike_mag and step_out > spike_mag
                    and back < 0.3):
                drop[i] = True
        out = out[~drop]
    return (
        out["phase"].to_numpy(),
        out["mag"].to_numpy(),
        out["emag"].to_numpy(),
    )


def load_pba_rband():
    """SN 2024pba r band in absolute magnitude, days from r-peak."""
    raw = pd.read_csv(os.path.join(OBS_DIR, "SN_B_lightcurve.dat"))
    sub = raw[raw["Filter"] == "r"].copy()
    det = sub[sub["Mag"].apply(lambda x: x != "-" and not pd.isna(x))]
    det["Mag"] = det["Mag"].astype(float)
    det["Error"] = det["Error"].astype(float)
    det["phase"] = (det["MJD"] - PBA_PEAK_MJD) / (1.0 + PBA_Z)
    det["Mabs"] = det["Mag"] - 5.0 * np.log10(PBA_DL_MPC * 1e5) - A_R
    return det["phase"].to_numpy(), det["Mabs"].to_numpy(), det[
        "Error"
    ].to_numpy()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(
        os.path.dirname(OBS_DIR), "SN_B_X_grid",
        "interacting_lc_comparison.png"))
    args = ap.parse_args()

    fig, ax = plt.subplots(figsize=(3.8, 3.4))

    # --- SN 2024pba ---
    ph, m, e = load_pba_rband()
    ax.errorbar(
        ph, m, yerr=e, fmt="o", ms=3.5, mfc="k", mec="k", mew=0.4,
        ecolor="k", elinewidth=0.5, capsize=0, zorder=20,
        label="SN 2024pba ($r$)",
    )

    handles = [
        Line2D([], [], marker="o", ls="", ms=3.5, color="k",
               label="SN 2024pba ($r$)")
    ]

    for name, band, D, z, peak, color, zo in ANALOGS:
        try:
            mjd, mag, emag, flag = load_analog(name, band)
        except FileNotFoundError:
            print(f"  {name}: no file, skipped")
            continue
        if len(mjd) == 0:
            print(f"  {name}: no {band}-band detections, skipped")
            continue
        n_raw = len(mjd)
        if flag is True:
            phase = mjd - mjd[np.argmin(mag)]
            Mabs = mag - 5.0 * np.log10(D * 1e5)
        else:
            if peak is None:
                peak = mjd[np.argmin(mag)]
            phase = (mjd - peak) / (1.0 + z)
            if flag == "abs":
                Mabs = mag
            else:
                Mabs = mag - 5.0 * np.log10(D * 1e5)
        phase, Mabs, emag = bin_lightcurve(phase, Mabs, emag)
        ax.plot(
            phase, Mabs, "o-", ms=2.2, lw=0.7, color=color,
            mec="none", alpha=0.85, zorder=zo,
        )
        handles.append(
            Line2D([], [], marker="o", ls="-", ms=2.2, lw=0.7,
                   color=color, label=f"{name} (${band}$)")
        )
        print(f"  {name}: {n_raw}->{len(phase)} {band}-band points, "
              f"phase {phase.min():.0f}..{phase.max():.0f} d")

    ax.axvline(BREAK_PHASE, color="tab:red", ls=":", lw=0.9, zorder=1)
    ax.text(
        BREAK_PHASE + 3, -20.2, "+80 d",
        color="tab:red", fontsize=6, rotation=90, va="top",
    )

    ax.set_xlabel("rest-frame days since peak")
    ax.set_ylabel("absolute magnitude")
    ax.invert_yaxis()
    ax.set_xlim(-60, 200)
    ax.set_ylim(-11.8, -20.6)
    ax.legend(handles=handles, fontsize=6, ncol=2, loc="lower left",
              frameon=False, handlelength=1.6, columnspacing=0.8,
              labelspacing=0.3)

    ax.grid(True, alpha=0.3)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.grid(True, which="minor", alpha=0.12)

    plt.tight_layout()
    fig.savefig(args.out, dpi=300, bbox_inches="tight")
    pdf = os.path.splitext(args.out)[0] + ".pdf"
    fig.savefig(pdf, bbox_inches="tight")
    print(f"Saved: {args.out}")
    print(f"Saved: {pdf}")


if __name__ == "__main__":
    main()
