#!/usr/bin/env python3
"""Plot observed SN_B lightcurve versus a STELLA SN_B model.

Computes SDSS g and r synthetic photometry from the STELLA .flx file and
overlays the observed g/r photometry of the SN_B event of SN 2024pba.
SN_A detections are shown faintly for context.

The STELLA model time axis is aligned to the observed SN_B r-band peak
(MJD 60867.446, brightest r-band detection) unless --shift is given.

Usage:
    python scripts/plot_sn_b_94w_vs_obs.py SN_B_v2_94W
    python scripts/plot_sn_b_94w_vs_obs.py SN_B_v2_94W --shift -5
"""
import argparse
import os
import struct

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

# STELLA constants
UR = 1.0e14
BOLTZK = 1.3806504e-16
UTP = 1.0e5
CS = 2.99792458e10
HPLANC = 1.054571628e-27
CLUMNU = np.log10(8.0 * UR**2 * (BOLTZK * UTP) ** 3 / (CS * HPLANC) ** 2)
UFREQ = 2.08366e15
NFREQ = 100
MFREQ = 130

freq = np.zeros(NFREQ + 1)
freq[0] = CS * 1.0e8 / (50000.0 * UFREQ)
freq[NFREQ] = CS * 1.0e8 / (1.0 * UFREQ)
Basis = (50000.0 / 1.0) ** (1.0 / NFREQ)
for i in range(1, NFREQ):
    freq[i] = freq[0] * Basis**i
freqmn = np.sqrt(freq[:-1] * freq[1:])
freqob = np.zeros(MFREQ)
freqob[:NFREQ] = freqmn
DLOGNU1 = 1.0 / np.log(freqmn[0] / freqmn[1])
for L in range(NFREQ, MFREQ):
    freqob[L] = freqmn[NFREQ - 1] * np.exp(-(L - NFREQ + 1) / DLOGNU1)

nu_hz = freqob * UFREQ
lam_A = CS * 1.0e8 / nu_hz

# SDSS filter response curves (Doi et al. 2010)
sdss_g_wave = np.array(
    [
        3000, 3500, 3700, 3800, 3900, 4000, 4100, 4200, 4300, 4400,
        4500, 4600, 4700, 4800, 4900, 5000, 5100, 5200, 5300, 5400,
        5500, 5600, 5700, 5800, 5900, 6000, 6100, 6200, 6300, 6400,
        6500, 6600,
    ]
)
sdss_g_resp = np.array(
    [
        0.00, 0.00, 0.00, 0.01, 0.04, 0.10, 0.20, 0.35, 0.52, 0.68,
        0.80, 0.88, 0.93, 0.96, 0.97, 0.96, 0.93, 0.88, 0.80, 0.70,
        0.58, 0.45, 0.32, 0.20, 0.10, 0.04, 0.01, 0.00, 0.00, 0.00,
        0.00, 0.00,
    ]
)

sdss_r_wave = np.array(
    [
        5000, 5200, 5400, 5600, 5800, 5900, 6000, 6100, 6200, 6300,
        6400, 6500, 6600, 6700, 6800, 6900, 7000, 7100, 7200, 7300,
        7400, 7500, 7600, 7700, 7800, 7900, 8000, 8100, 8200, 8300,
        8400, 8500, 8600, 8700, 8800, 8900, 9000, 9100, 9200, 9300,
        9400, 9500, 9600, 9700, 9800, 9900, 10000,
    ]
)
sdss_r_resp = np.array(
    [
        0.00, 0.00, 0.01, 0.05, 0.12, 0.18, 0.26, 0.35, 0.45, 0.55,
        0.65, 0.73, 0.80, 0.85, 0.89, 0.92, 0.94, 0.95, 0.96, 0.96,
        0.95, 0.94, 0.92, 0.89, 0.85, 0.80, 0.74, 0.67, 0.59, 0.50,
        0.41, 0.32, 0.23, 0.15, 0.08, 0.03, 0.01, 0.00, 0.00, 0.00,
        0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
    ]
)

g_interp = interp1d(sdss_g_wave, sdss_g_resp, bounds_error=False, fill_value=0.0)
r_interp = interp1d(sdss_r_wave, sdss_r_resp, bounds_error=False, fill_value=0.0)
g_trans = g_interp(lam_A)
r_trans = r_interp(lam_A)

# Distance modulus: SN2024pba is at 83.5 Mpc
DISTANCE_MPC = 83.5
DM = 5.0 * np.log10(DISTANCE_MPC * 1.0e6) - 5.0
PC_CM = 3.0856776e18
D_CM = DISTANCE_MPC * 1.0e6 * PC_CM

# Brightest r-band detection of the SN_B event
SN_B_R_PEAK_MJD = 60867.446


def read_flx(flx_path):
    """Parse STELLA binary *.flx file."""
    with open(flx_path, "rb") as f:
        Nrec = struct.unpack("i", f.read(4))[0]
        Lsave = struct.unpack("i", f.read(4))[0]
        times = np.zeros(Lsave)
        nfrus = np.zeros(Lsave, dtype=int)
        flsave = np.zeros((MFREQ, Lsave))
        for i in range(Lsave):
            times[i] = struct.unpack("f", f.read(4))[0]
            nfrus[i] = struct.unpack("i", f.read(4))[0]
            flsave[:, i] = np.array(struct.unpack(f"{MFREQ}d", f.read(8 * MFREQ)))
        end_marker = struct.unpack("i", f.read(4))[0]
        assert end_marker == Nrec, "FLX record size mismatch"
    return times, nfrus, flsave


def compute_sdss_mags(flx_path):
    """Compute SDSS g and r apparent magnitudes from a STELLA .flx file."""
    times, nfrus_arr, flsave = read_flx(flx_path)
    n_epochs = len(times)

    g_mag = np.full(n_epochs, np.nan)
    r_mag = np.full(n_epochs, np.nan)

    for i in range(n_epochs):
        F = flsave[:, i]
        L_nu = F * freqob**3 * 10.0**CLUMNU
        F_nu = L_nu / (4.0 * np.pi * D_CM**2)

        for band, trans in [("g", g_trans), ("r", r_trans)]:
            mask = trans > 1e-10
            if mask.sum() < 3:
                continue
            nu_band = nu_hz[mask]
            F_band = F_nu[mask]
            T_band = trans[mask]
            num = np.trapezoid(F_band * T_band, nu_band)
            den = np.trapezoid(T_band, nu_band)
            if den > 0 and num > 0:
                m = -2.5 * np.log10(num / den) - 48.6
                if band == "g":
                    g_mag[i] = m
                else:
                    r_mag[i] = m

    return times, g_mag, r_mag


def read_obs_lightcurve(filepath):
    """Read observational lightcurve CSV file (detections only)."""
    mjd, mag, err, filt = [], [], [], []
    with open(filepath) as f:
        f.readline()
        for line in f:
            parts = line.strip().split(",")
            if len(parts) < 6:
                continue
            try:
                m = float(parts[0])
            except ValueError:
                continue
            mag_val, err_val = parts[2], parts[3]
            band = parts[4].strip()
            if mag_val == "-" or err_val == "-" or band not in ("g", "r"):
                continue
            mjd.append(m)
            mag.append(float(mag_val))
            err.append(float(err_val))
            filt.append(band)
    return np.array(mjd), np.array(mag), np.array(err), np.array(filt)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("model_dir", help="model directory (e.g. SN_B_v2_94W)")
    ap.add_argument(
        "--shift",
        type=float,
        default=None,
        help="extra days added to the model MJD axis after peak alignment",
    )
    ap.add_argument("--output", "-o", default=None)
    args = ap.parse_args()

    model_dir = args.model_dir
    model_name = os.path.basename(model_dir.rstrip("/"))
    flx_path = os.path.join(model_dir, "stella", f"{model_name}.flx")
    obs_dir = "/viper/ptmp1/sbrennan/SN2024pba/observational_data"
    outpath = args.output or os.path.join(
        model_dir, "plots", f"{model_name}_vs_obs.png"
    )

    obs_b_mjd, obs_b_mag, obs_b_err, obs_b_filt = read_obs_lightcurve(
        f"{obs_dir}/SN_B_lightcurve.dat"
    )
    print(f"SN_B: {len(obs_b_mjd)} detections")
    sn_a_path = f"{obs_dir}/SN_A_lightcurve.dat"
    obs_a = None
    if os.path.exists(sn_a_path):
        obs_a = read_obs_lightcurve(sn_a_path)
        print(f"SN_A: {len(obs_a[0])} detections")

    print(f"Computing STELLA SDSS photometry from {flx_path}...")
    stella_t, stella_g, stella_r = compute_sdss_mags(flx_path)
    print(f"  {len(stella_t)} epochs, t=[{stella_t[0]:.3f}, {stella_t[-1]:.3f}] d")

    valid_r = np.isfinite(stella_r)
    if np.any(valid_r):
        r_peak_idx = np.argmin(stella_r[valid_r])
        stella_r_peak_t = stella_t[valid_r][r_peak_idx]
        stella_r_peak_mag = stella_r[valid_r][r_peak_idx]
        print(
            f"  STELLA r-band peak: t={stella_r_peak_t:.2f} d, "
            f"m_r={stella_r_peak_mag:.3f} (apparent)"
        )
        t0_mjd = SN_B_R_PEAK_MJD - stella_r_peak_t
    else:
        t0_mjd = SN_B_R_PEAK_MJD - 8.0
        print("  WARNING: no valid r-band epochs, using default offset")

    if args.shift is not None:
        t0_mjd += args.shift
        print(f"  extra shift applied: {args.shift:+.2f} d")

    stella_mjd = t0_mjd + stella_t

    fig, axes = plt.subplots(2, 1, figsize=(14, 10), sharex=True)

    for row, band in enumerate(["g", "r"]):
        ax = axes[row]
        stella_mag = stella_g if band == "g" else stella_r

        if obs_a is not None:
            mask_a = obs_a[3] == band
            ax.errorbar(
                obs_a[0][mask_a],
                obs_a[1][mask_a],
                yerr=obs_a[2][mask_a],
                fmt="o",
                color="crimson",
                ms=4,
                capsize=2,
                mfc="white",
                alpha=0.5,
                label="SN_A obs",
                zorder=4,
            )

        mask_b = obs_b_filt == band
        ax.errorbar(
            obs_b_mjd[mask_b],
            obs_b_mag[mask_b],
            yerr=obs_b_err[mask_b],
            fmt="s",
            color="royalblue",
            ms=5,
            capsize=3,
            mfc="white",
            label="SN_B obs",
            zorder=5,
        )

        valid = np.isfinite(stella_mag)
        ax.plot(
            stella_mjd[valid],
            stella_mag[valid],
            "-",
            color="darkgreen",
            lw=2.5,
            alpha=0.9,
            label=f"STELLA {model_name}",
            zorder=6,
        )

        ax.axvline(SN_B_R_PEAK_MJD, color="royalblue", ls=":", alpha=0.3, zorder=1)

        ax.set_ylabel(f"{band} magnitude", fontsize=13)
        ax.set_title(f"SDSS {band} band", fontsize=14)
        ax.invert_yaxis()
        ax.legend(loc="lower left", fontsize=10, ncol=3)
        ax.grid(True, alpha=0.3)

        obs_b = obs_b_mag[mask_b]
        all_mags = np.concatenate([obs_b, stella_mag[np.isfinite(stella_mag)]])
        all_mags = all_mags[np.isfinite(all_mags) & (all_mags < 25)]
        if len(all_mags) > 0:
            ax.set_ylim(max(all_mags.max() + 1.0, 22.0), all_mags.min() - 1.0)

    axes[1].set_xlabel("MJD", fontsize=13)
    axes[0].set_xlim(obs_b_mjd.min() - 20, obs_b_mjd.max() + 10)

    fig.suptitle(
        f"SN2024pba SN_B event: g and r band lightcurves - obs vs {model_name}",
        fontsize=15,
        y=0.98,
    )
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    plt.savefig(outpath, dpi=150, bbox_inches="tight")
    print(f"Plot saved to: {outpath}")
    plt.close(fig)


if __name__ == "__main__":
    main()
