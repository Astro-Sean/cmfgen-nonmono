#!/usr/bin/env python3
"""Plot r-band lightcurves for the SN_B TcurA clock-start test models.

The three models in SN_B_tcurA_test/ are identical except for tcura_days
(the STELLA TcurA clock value at run start): 0, 5 and 10 days.

Model magnitudes are Johnson-Cousins R absolute magnitudes read from the
STELLA .tt files (cols: time Tbb rbb Teff Rlast_sc Rtau23 Mbol MU MB MV
MI MR Mbolavg gdepos).  Model time is OBS.TIME, days since shell launch.

The observed SN 2024pba r-band (SDSS r PSF photometry) is overlaid for
context, phased to the adopted SN_B explosion epoch MJD 60805 so that
both axes read "days since launch".  Johnson R is not identical to SDSS
r; expect ~0.1-0.3 mag systematic offsets.
"""
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

GRID_DIR = "/viper/ptmp1/sbrennan/SN2024pba/SN_B_tcurA_test"
OBS_FILE = "/viper/ptmp1/sbrennan/SN2024pba/observational_data/SN_B_lightcurve.dat"
OUT_FILE = os.path.join(GRID_DIR, "tcura_rband_lightcurves.png")

MODELS = [
    ("SB_M10_V10k_n10_W1e-5_TcurA00", "TcurA = 0 d", "tab:blue"),
    ("SB_M10_V10k_n10_W1e-5_TcurA05", "TcurA = 5 d", "tab:orange"),
    ("SB_M10_V10k_n10_W1e-5_TcurA10", "TcurA = 10 d", "tab:green"),
]

# Same photometric constants as plot_sn_b_grid_vs_obs.py.
DM = 34.52
EBV = 0.15
A_R = 2.285 * EBV
SNB_EXP_MJD = 60805.0


def load_model_tt(tt_path):
    """Load (time, R) from a STELLA .tt file."""
    rows = []
    with open(tt_path) as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("time") or s.startswith("XNianm"):
                continue
            try:
                row = [float(x) for x in s.split()]
            except ValueError:
                continue
            if len(row) >= 14:
                rows.append((row[0], row[11]))
    arr = np.array(rows)
    return pd.DataFrame({"time": arr[:, 0], "R": arr[:, 1]})


def main():
    fig, ax = plt.subplots(figsize=(10, 6.5))

    for name, label, color in MODELS:
        tt = os.path.join(GRID_DIR, name, "stella", f"{name}.tt")
        df = load_model_tt(tt)
        df = df[df["time"] >= -20].copy()
        ax.plot(df["time"], df["R"], color=color, lw=1.6, label=label)
        i = df["R"].idxmin()
        print(
            f"{name}: R peak M={df['R'].min():.2f} at t={df.loc[i, 'time']:.1f} d"
        )

    raw = pd.read_csv(OBS_FILE)
    obs = raw[(raw["Filter"] == "r") & (raw["Mag"] != "-")].copy()
    obs["Mag"] = obs["Mag"].astype(float)
    obs["Mabs"] = obs["Mag"] - DM - A_R
    obs["phase"] = obs["MJD"] - SNB_EXP_MJD
    ax.errorbar(
        obs["phase"],
        obs["Mabs"],
        yerr=obs["Error"].astype(float),
        fmt="o",
        color="tab:red",
        ms=5,
        mec="k",
        mew=0.7,
        capsize=2,
        elinewidth=0.7,
        alpha=0.85,
        zorder=20,
        label="SN 2024pba r (PSF)",
    )

    ax.set_xlabel("days since shell launch (obs: MJD - 60805)")
    ax.set_ylabel("absolute magnitude (model: Johnson R, obs: SDSS r)")
    ax.set_title(
        "SN_B TcurA clock-start test: "
        "M$_{\\rm shell}$=10 M$_\\odot$, v$_{\\rm shell}$=10000 km/s, "
        "n$_{\\rho}$=10, $\\dot{M}$=1e-5 M$_\\odot$/yr"
    )
    ax.invert_yaxis()
    ax.set_xlim(-20, 120)
    ax.grid(alpha=0.3)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUT_FILE, dpi=150)
    print(f"Wrote {OUT_FILE}")


if __name__ == "__main__":
    main()
