#!/usr/bin/env python3
"""Compare two CMFGEN observed spectra (obs_fin_10 format).

Reproduces the TypeIIn KNOP vs KNOP_ES comparison figure:
full optical spectrum, H-alpha zoom, and flux ratio.

The H-alpha zoom y-limit is set from the data so the emission peak is
not clipped.
"""

import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_cmfgen_obs_spectrum import read_obs_fin_10

KNOP_FILE = "/ptmp/sbrennan/SN2024pba/TypeIIn_KNOP/cmfgen/10.00days/obs/obs_fin_10"
KNOP_ES_FILE = (
    "/ptmp/sbrennan/SN2024pba/TypeIIn_KNOP_ES/cmfgen/10.00days/obs/obs_fin_10"
)
OUT_FILE = (
    "/ptmp/sbrennan/SN2024pba/TypeIIn_KNOP_ES/cmfgen/10.00days/plots/"
    "TypeIIn_10d_KNOP_ES_vs_KNOP_obs_spectrum.png"
)

# Line identifications (rest wavelengths, Angstrom)
LINE_IDS = [
    (3933.7, "Ca II K"),
    (3968.5, "Ca II H"),
    (4101.7, "H$\\delta$"),
    (4340.5, "H$\\gamma$"),
    (4861.3, "H$\\beta$"),
    (4923.9, "Fe II"),
    (5018.4, "Fe II"),
    (5169.0, "Fe II"),
    (5875.6, "He I"),
    (6562.8, "H$\\alpha$"),
    (7774.0, "O I"),
    (8498.0, "Ca II"),
]


def main():
    wl1, fl1 = read_obs_fin_10(KNOP_FILE)
    wl2, fl2 = read_obs_fin_10(KNOP_ES_FILE)

    # Sort by wavelength (obs_fin is frequency-ordered)
    s1, s2 = np.argsort(wl1), np.argsort(wl2)
    wl1, fl1 = wl1[s1], fl1[s1]
    wl2, fl2 = wl2[s2], fl2[s2]

    fig, axes = plt.subplots(2, 1, figsize=(14, 9))

    # Panel 1: full optical spectrum.  The y-limit is set from the 99th
    # percentile of the visible flux - i.e. the continuum level - so narrow
    # emission-line cores (H-alpha, H-beta, Fe II) may clip while the
    # continuum stays readable.
    ax = axes[0]
    ax.plot(wl1, fl1, color="k", lw=0.6, label="TypeIIn_KNOP (no ES)")
    ax.plot(wl2, fl2, color="r", lw=0.6, alpha=0.7, label="TypeIIn_KNOP_ES")
    ax.set_xlim(3200, 9000)
    wl_vis = (wl1 > 3200) & (wl1 < 9000)
    ytop = np.nanpercentile(fl1[wl_vis], 99) * 1.15
    ax.set_ylim(0, ytop)
    for wl_id, name in LINE_IDS:
        ax.axvline(wl_id, color="b", ls="--", lw=0.8, alpha=0.5)
        ax.text(wl_id, ytop * 0.97, name, rotation=90, fontsize=7,
                color="b", ha="right", va="top")
    ax.set_ylabel(r"$F_\lambda$ (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)")
    ax.set_title("CMFGEN Observed Spectrum, 10 days - ES vs non-ES")
    ax.legend(loc="upper right")
    ax.grid(True, alpha=0.3)

    # Panel 2: H-alpha zoom - y-limit from data so the peak is not clipped
    ax = axes[1]
    ha_lo, ha_hi = 6200.0, 6900.0
    m1 = (wl1 > ha_lo) & (wl1 < ha_hi)
    m2 = (wl2 > ha_lo) & (wl2 < ha_hi)
    ymax_ha = max(np.nanmax(fl1[m1]), np.nanmax(fl2[m2]))
    ax.plot(wl1[m1], fl1[m1], color="k", lw=0.9, label="KNOP")
    ax.plot(wl2[m2], fl2[m2], color="r", lw=0.9, alpha=0.8, label="KNOP_ES")
    ax.axvline(6562.8, color="r", ls="--", lw=0.8, alpha=0.5)
    ax.set_xlim(ha_lo, ha_hi)
    ax.set_ylim(0, ymax_ha * 1.08)
    ax.set_xlabel(r"Wavelength ($\AA$)")
    ax.set_ylabel(r"$F_\lambda$ (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)")
    ax.set_title("H$\\alpha$ region zoom")
    ax.legend(loc="upper right")
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(OUT_FILE, dpi=150)
    print(f"Saved {OUT_FILE}")
    print(f"H-alpha peaks: KNOP={np.nanmax(fl1[m1]):.0f}, "
          f"KNOP_ES={np.nanmax(fl2[m2]):.0f} erg/s/cm2/A")


if __name__ == "__main__":
    main()
