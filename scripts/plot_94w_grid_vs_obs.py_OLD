#!/usr/bin/env python3
"""Plot SN_B_v2_94W grid lightcurves vs SN 2024pba SN_B observations.

Reads the parallel-STELLA .tt files (extended band set: SDSS g/r and
ATLAS o are native .tt columns), converts to apparent magnitudes with
mu=34.61 (83.5 Mpc), and aligns each model's r-band peak to the observed
SN_B r-band peak (MJD 60867.446).

Usage:
    python scripts/plot_94w_grid_vs_obs.py
    python scripts/plot_94w_grid_vs_obs.py --models 94W_V4k_W1e-2 94W_V4k_Wmatch
"""
import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BASE = Path("/viper/ptmp1/sbrennan/SN2024pba")
GRID_DIR = BASE / "SN_B_v2_94W_grid"
OBS_FILE = BASE / "observational_data" / "SN_B_lightcurve.dat"
BASELINE_TT = BASE / "SN_B_v2_94W" / "stella" / "SN_B_v2_94W.tt"

DISTANCE_MPC = 83.5
DM = 5.0 * np.log10(DISTANCE_MPC * 1.0e6) - 5.0
SN_B_R_PEAK_MJD = 60867.446

# .tt column indices (0-based) for the parallel build's extended bands
TT_COLS = {"g": 13, "r": 14, "o": 21, "bol": 6}


def read_tt(tt_path):
    """Parse a parallel-build STELLA .tt file -> dict of band arrays."""
    rows = []
    with open(tt_path) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("time") or s.startswith("day"):
                continue
            parts = s.split()
            if len(parts) < 23:
                continue
            try:
                rows.append([float(x) for x in parts])
            except ValueError:
                continue
    arr = np.asarray(rows)
    out = {"time": arr[:, 0]}
    for band, idx in TT_COLS.items():
        out[band] = arr[:, idx]
    return out


def read_obs(path):
    """Read SN_B_lightcurve.dat -> (mjd, mag, err, filt), detections only."""
    mjd, mag, err, filt = [], [], [], []
    with open(path) as f:
        f.readline()
        for line in f:
            p = line.strip().split(",")
            if len(p) < 6 or p[2] == "-" or p[3] == "-":
                continue
            band = p[4].strip()
            if band not in ("g", "r"):
                continue
            mjd.append(float(p[0]))
            mag.append(float(p[2]))
            err.append(float(p[3]))
            filt.append(band)
    return np.array(mjd), np.array(mag), np.array(err), np.array(filt)


def align(model, shift=0.0, t_min=0.0):
    """Model days -> MJD, aligning model r peak to the observed r peak.

    t_min > 0 ignores the early startup flash when finding the peak.
    """
    t, r = model["time"], model["r"]
    sel = t > t_min if t_min > 0 else np.full(t.shape, True)
    t_peak = t[sel][np.argmin(r[sel])]
    return SN_B_R_PEAK_MJD - t_peak + shift


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--models",
        nargs="*",
        default=None,
        help="grid model names (default: all dirs with a .tt)",
    )
    ap.add_argument("--baseline", action="store_true", default=True)
    ap.add_argument("--no-baseline", dest="baseline", action="store_false")
    ap.add_argument(
        "--align-broad",
        action="store_true",
        help="align the post-flash broad peak (t>10d) instead of the spike",
    )
    ap.add_argument("-o", "--output", default=None)
    args = ap.parse_args()

    obs_mjd, obs_mag, obs_err, obs_filt = read_obs(OBS_FILE)

    models = {}
    names = args.models
    if names is None:
        names = sorted(
            d.name for d in GRID_DIR.iterdir() if (d / "stella").exists()
        )
    for name in names:
        tts = list((GRID_DIR / name / "stella").glob("*.tt"))
        if not tts:
            print(f"  {name}: no .tt yet, skipped")
            continue
        models[name] = read_tt(tts[0])
        print(f"  {name}: {len(models[name]['time'])} epochs")

    if args.baseline and BASELINE_TT.exists():
        models = {"SN_B_v2_94W (baseline)": read_tt(BASELINE_TT), **models}

    cmap = plt.get_cmap("tab10")
    colors = {name: cmap(i % 10) for i, name in enumerate(models)}

    fig, axes = plt.subplots(3, 1, figsize=(13, 11), sharex=True)
    band_titles = {"g": "SDSS g", "r": "SDSS r", "o": "ATLAS o (model only)"}

    for row, band in enumerate(["g", "r", "o"]):
        ax = axes[row]
        mask = obs_filt == band
        if mask.any():
            ax.errorbar(
                obs_mjd[mask],
                obs_mag[mask],
                yerr=obs_err[mask],
                fmt="s",
                color="k",
                ms=4,
                capsize=2,
                mfc="white",
                label="SN_B obs",
                zorder=8,
            )
        for name, m in models.items():
            t0 = align(m, t_min=10.0 if args.align_broad else 0.0)
            ls = "--" if "baseline" in name else "-"
            ax.plot(
                t0 + m["time"],
                m[band] + DM,
                ls=ls,
                color=colors[name],
                lw=1.6,
                label=name,
            )
        ax.invert_yaxis()
        ax.set_ylabel(f"{band} mag")
        ax.set_title(band_titles[band], fontsize=12)
        ax.grid(alpha=0.3)
        ax.axvline(SN_B_R_PEAK_MJD, color="gray", ls=":", alpha=0.4)

    axes[0].legend(fontsize=8, ncol=2, loc="lower left")
    axes[2].set_xlabel("MJD")
    axes[0].set_xlim(obs_mjd.min() - 15, obs_mjd.min() + 200)

    fig.suptitle(
        "SN2024pba SN_B event vs SN_B_v2_94W grid "
        "(r-peak aligned, mu=%.2f)" % DM,
        fontsize=13,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    out = args.output or str(GRID_DIR / "grid_vs_obs_gro.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"saved {out}")

    # Peak summary
    print("\nmodel                    r_peak(abs)  g-r@peak  decay50d(r)")
    for name, m in models.items():
        i = np.argmin(m["r"])
        t = m["time"]
        r = m["r"]
        t50 = t[i] + 50.0
        j = np.searchsorted(t, t50)
        decay = r[min(j, len(r) - 1)] - r[i] if j < len(r) else np.nan
        print(
            f"{name:24s} {r[i]:8.2f}   {m['g'][i]-r[i]:+6.2f}   "
            f"{decay:+6.2f} mag/50d"
        )


if __name__ == "__main__":
    main()
