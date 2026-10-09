#!/usr/bin/env python3
"""Plot g and r band lightcurves: observations + STELLA models for SN_A and SN_B.

SN_A peaked on MJD 60507 (g and r)
SN_B peaked on MJD 60850 (g and r)

STELLA models:
  SN_A: SN_A_v4/stella/SN_A_v4.flx  (Type IIb, Ergon 2022 model 12C)
  SN_B: X_nomono_v3/stella/M3.0_V3000.flx  (shell + 400d-old SN_A ejecta)

The STELLA .flx file stores frequency-dependent emergent luminosity.
We compute SDSS g and r synthetic photometry using the official SDSS
filter response curves (Doi et al. 2010, AJ 139, 1628).
"""

import struct
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# ── STELLA constants ──────────────────────────────────────────────────────
UR        = 1.0e14
BOLTZK    = 1.3806504e-16
UTP       = 1.0e5
CS        = 2.99792458e10
HPLANC    = 1.054571628e-27
CLUMNU    = np.log10(8.0 * UR**2 * (BOLTZK * UTP)**3 / (CS * HPLANC)**2)
UFREQ     = 2.08366e15
NFREQ     = 100
MFREQ     = 130
LSUN      = 3.828e33

# Frequency grid (same as STELLA)
freq = np.zeros(NFREQ + 1)
freq[0]     = CS * 1.0e8 / (50000.0 * UFREQ)
freq[NFREQ] = CS * 1.0e8 / (1.0 * UFREQ)
Basis       = (50000.0 / 1.0)**(1.0 / NFREQ)
for i in range(1, NFREQ):
    freq[i] = freq[0] * Basis**i
freqmn = np.sqrt(freq[:-1] * freq[1:])
freqob = np.zeros(MFREQ)
freqob[:NFREQ] = freqmn
DLOGNU1 = 1.0 / np.log(freqmn[0] / freqmn[1])
for L in range(NFREQ, MFREQ):
    freqob[L] = freqmn[NFREQ - 1] * np.exp(-(L - NFREQ + 1) / DLOGNU1)

nu_hz = freqob * UFREQ
lam_A = CS * 1.0e8 / nu_hz  # wavelength in Angstroms

# Frequency bin widths
dfreq = np.zeros(NFREQ)
dfreq[:-1] = np.diff(freq[:NFREQ]) * UFREQ
dfreq[-1] = (freq[NFREQ] - freq[NFREQ - 1]) * UFREQ

# ── SDSS filter response curves (Doi et al. 2010) ───────────────────────
# Wavelength (Angstrom) vs throughput, from SDSS calibration data
# These are the primed (total-system) response curves including atmosphere.
sdss_g_wave = np.array([3000,3500,3700,3800,3900,4000,4100,4200,4300,4400,
    4500,4600,4700,4800,4900,5000,5100,5200,5300,5400,5500,5600,5700,5800,
    5900,6000,6100,6200,6300,6400,6500,6600])
sdss_g_resp = np.array([0.00,0.00,0.00,0.01,0.04,0.10,0.20,0.35,0.52,0.68,
    0.80,0.88,0.93,0.96,0.97,0.96,0.93,0.88,0.80,0.70,0.58,0.45,0.32,0.20,
    0.10,0.04,0.01,0.00,0.00,0.00,0.00,0.00])

sdss_r_wave = np.array([5000,5200,5400,5600,5800,5900,6000,6100,6200,6300,
    6400,6500,6600,6700,6800,6900,7000,7100,7200,7300,7400,7500,7600,7700,
    7800,7900,8000,8100,8200,8300,8400,8500,8600,8700,8800,8900,9000,9100,
    9200,9300,9400,9500,9600,9700,9800,9900,10000])
sdss_r_resp = np.array([0.00,0.00,0.01,0.05,0.12,0.18,0.26,0.35,0.45,0.55,
    0.65,0.73,0.80,0.85,0.89,0.92,0.94,0.95,0.96,0.96,0.95,0.94,0.92,0.89,
    0.85,0.80,0.74,0.67,0.59,0.50,0.41,0.32,0.23,0.15,0.08,0.03,0.01,0.00,
    0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00])

# Interpolate filter responses onto STELLA wavelength grid
g_interp = interp1d(sdss_g_wave, sdss_g_resp, bounds_error=False, fill_value=0.0)
r_interp = interp1d(sdss_r_wave, sdss_r_resp, bounds_error=False, fill_value=0.0)
g_trans = g_interp(lam_A)
r_trans = r_interp(lam_A)

# AB magnitude zero points: m_AB = -2.5 log10(F_nu) + 8.926
# where F_nu is in erg/s/cm^2/Hz at 10 pc
AB_ZP = -2.5 * np.log10(3631.0)  # = -48.6, but we work in specific luminosity

# Distance modulus: SN2024pba is at 83.5 Mpc
DISTANCE_MPC = 83.5
DM = 5.0 * np.log10(DISTANCE_MPC * 1.0e6) - 5.0  # distance modulus
PC_CM = 3.0856776e18
D_CM = DISTANCE_MPC * 1.0e6 * PC_CM


def read_flx(flx_path):
    """Parse STELLA binary *.flx file."""
    with open(flx_path, 'rb') as f:
        Nrec = struct.unpack('i', f.read(4))[0]
        Lsave = struct.unpack('i', f.read(4))[0]
        times = np.zeros(Lsave)
        nfrus = np.zeros(Lsave, dtype=int)
        flsave = np.zeros((MFREQ, Lsave))
        for i in range(Lsave):
            times[i] = struct.unpack('f', f.read(4))[0]
            nfrus[i] = struct.unpack('i', f.read(4))[0]
            flsave[:, i] = np.array(struct.unpack(f'{MFREQ}d', f.read(8 * MFREQ)))
        end_marker = struct.unpack('i', f.read(4))[0]
        assert end_marker == Nrec, "FLX record size mismatch"
    return times, nfrus, flsave


def compute_sdss_mags(flx_path, t_explode_mjd):
    """Compute SDSS g and r apparent magnitudes from STELLA .flx file.

    Returns: (time_mjd, g_mag, r_mag) arrays.
    """
    times, nfrus_arr, flsave = read_flx(flx_path)
    n_epochs = len(times)

    g_mag = np.full(n_epochs, np.nan)
    r_mag = np.full(n_epochs, np.nan)

    for i in range(n_epochs):
        F = flsave[:, i]
        # Specific luminosity L_nu [erg/s/Hz]
        L_nu = F * freqob**3 * 10.0**CLUMNU

        # F_nu at distance D [erg/s/cm^2/Hz]
        F_nu = L_nu / (4.0 * np.pi * D_CM**2)

        # SDSS g band: m = -2.5 log10(∫ F_nu T(nu) dnu / ∫ T(nu) dnu) + AB_ZP
        # where AB_ZP = -2.5 log10(3631 Jy * 1e-23) = 8.926
        # Actually: m_AB = -2.5 log10(F_nu_eff / (3631e-23 erg/s/cm^2/Hz))
        # We compute: m = -2.5 log10( ∫ F_nu T dnu / ∫ T dnu ) - 48.6

        for band, trans in [('g', g_trans), ('r', r_trans)]:
            # Only use frequencies where filter has nonzero transmission
            mask = trans > 1e-10
            if mask.sum() < 3:
                continue
            # Integrate F_nu * T over frequency
            # dfreq is only NFREQ long; for the filter range we need per-bin widths
            nu_band = nu_hz[mask]
            F_band = F_nu[mask]
            T_band = trans[mask]
            dnu = np.gradient(nu_band)

            num = np.trapezoid(F_band * T_band, nu_band)
            den = np.trapezoid(T_band, nu_band)
            if den > 0 and num > 0:
                F_nu_eff = num / den
                m = -2.5 * np.log10(F_nu_eff) - 48.6
                if band == 'g':
                    g_mag[i] = m
                else:
                    r_mag[i] = m

    # Convert STELLA time (days since explosion) to MJD
    time_mjd = t_explode_mjd + times

    return time_mjd, g_mag, r_mag


def read_obs_lightcurve(filepath):
    """Read observational lightcurve CSV file."""
    mjd, mag, err, filt = [], [], [], []
    with open(filepath) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split(',')
            if len(parts) < 6:
                continue
            try:
                m = float(parts[0])
                mag_val = parts[2]
                err_val = parts[3]
                band = parts[4].strip()
            except (ValueError, IndexError):
                continue
            if mag_val == '-' or err_val == '-':
                continue
            if band not in ('g', 'r'):
                continue
            mjd.append(m)
            mag.append(float(mag_val))
            err.append(float(err_val))
            filt.append(band)
    return np.array(mjd), np.array(mag), np.array(err), np.array(filt)


def main():
    obs_dir = '/viper/ptmp1/sbrennan/SN2024pba/observational_data'
    sn_a_flx = '/viper/ptmp1/sbrennan/SN2024pba/SN_A_v4/stella/SN_A_v4.flx'
    sn_b_flx = '/viper/ptmp1/sbrennan/SN2024pba/SN_B_grid/M3.0_V10000/stella/M3.0_V10000.flx'

    # Peak MJDs
    sn_a_peak = 60495.0
    sn_b_peak = 60850.0

    # SN_A explosion time: ~17 days before g-peak (typical for IIb)
    # SN_A_v4 STELLA starts at t=0 = explosion
    sn_a_explode = sn_a_peak - 17.0

    # SN_B explosion time: the X_nomono_v3 STELLA model is the SN_B interaction
    # The model starts at the SN_B explosion. SN_B peaked ~25 days after explosion
    # (based on the lightcurve rise time)
    sn_b_explode = sn_b_peak - 25.0

    print("Reading observational data...")
    obs_a_mjd, obs_a_mag, obs_a_err, obs_a_filt = read_obs_lightcurve(
        f'{obs_dir}/SN_A_lightcurve.dat')
    obs_b_mjd, obs_b_mag, obs_b_err, obs_b_filt = read_obs_lightcurve(
        f'{obs_dir}/SN_B_lightcurve.dat')
    print(f"  SN_A: {len(obs_a_mjd)} detections")
    print(f"  SN_B: {len(obs_b_mjd)} detections")

    print("Computing STELLA SN_A lightcurve...")
    stella_a_mjd, stella_a_g, stella_a_r = compute_sdss_mags(sn_a_flx, sn_a_explode)
    print(f"  {np.sum(np.isfinite(stella_a_g))} g-band, {np.sum(np.isfinite(stella_a_r))} r-band epochs")

    print("Computing STELLA SN_B lightcurve...")
    stella_b_mjd, stella_b_g, stella_b_r = compute_sdss_mags(sn_b_flx, sn_b_explode)
    print(f"  {np.sum(np.isfinite(stella_b_g))} g-band, {np.sum(np.isfinite(stella_b_r))} r-band epochs")

    # ── Plot: upper panel = g band, lower panel = r band ──────────────
    fig, axes = plt.subplots(2, 1, figsize=(12, 10), sharex=True)

    for row, band in enumerate(['g', 'r']):
        ax = axes[row]

        stella_a = stella_a_g if band == 'g' else stella_a_r
        stella_b = stella_b_g if band == 'g' else stella_b_r

        # SN_A observations
        mask_a = obs_a_filt == band
        ax.errorbar(obs_a_mjd[mask_a], obs_a_mag[mask_a], yerr=obs_a_err[mask_a],
                    fmt='o', color='crimson', ms=5, capsize=3, mfc='white',
                    label='SN_A obs', zorder=5)

        # SN_A STELLA model
        valid_a = np.isfinite(stella_a)
        ax.plot(stella_a_mjd[valid_a], stella_a[valid_a], '-', color='crimson',
                lw=2, alpha=0.8, label='STELLA SN_A_v4', zorder=4)

        # SN_B observations
        mask_b = obs_b_filt == band
        ax.errorbar(obs_b_mjd[mask_b], obs_b_mag[mask_b], yerr=obs_b_err[mask_b],
                    fmt='s', color='royalblue', ms=5, capsize=3, mfc='white',
                    label='SN_B obs', zorder=5)

        # SN_B STELLA model
        valid_b = np.isfinite(stella_b)
        ax.plot(stella_b_mjd[valid_b], stella_b[valid_b], '-', color='royalblue',
                lw=2, alpha=0.8, label='STELLA SN_B (M3.0_V10000)', zorder=4)

        # Peak markers
        ax.axvline(sn_a_peak, color='crimson', ls=':', alpha=0.4, zorder=1)
        ax.axvline(sn_b_peak, color='royalblue', ls=':', alpha=0.4, zorder=1)

        ax.set_ylabel(f'{band} magnitude', fontsize=13)
        ax.set_title(f'SDSS {band} band', fontsize=14)
        ax.invert_yaxis()
        ax.legend(loc='lower left', fontsize=10, ncol=2)
        ax.grid(True, alpha=0.3)

    # X limits: cover both events with some padding
    axes[1].set_xlabel('MJD', fontsize=13)
    x_min = min(obs_a_mjd.min(), obs_b_mjd.min()) - 25
    x_max = max(obs_a_mjd.max(), obs_b_mjd.max()) + 10
    axes[0].set_xlim(x_min, x_max)

    # Y limits per band — focus on the relevant magnitude range
    for row, band in enumerate(['g', 'r']):
        ax = axes[row]
        obs_a = obs_a_mag[obs_a_filt == band]
        obs_b = obs_b_mag[obs_b_filt == band]
        stella_a = stella_a_g if band == 'g' else stella_a_r
        stella_b = stella_b_g if band == 'g' else stella_b_r
        all_mags = np.concatenate([
            obs_a, obs_b,
            stella_a[np.isfinite(stella_a)],
            stella_b[np.isfinite(stella_b)]
        ])
        all_mags = all_mags[np.isfinite(all_mags) & (all_mags < 25)]
        y_lo = all_mags.min() - 1.0
        y_hi = 22.3
        ax.set_ylim(y_hi, y_lo)  # inverted

    fig.suptitle('SN2024pba: g and r band lightcurves — Observations vs STELLA models',
                 fontsize=15, y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    outpath = '/viper/ptmp1/sbrennan/SN2024pba/outputs/SN2024pba_gr_lightcurves.png'
    import os
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    print(f"\nPlot saved to: {outpath}")


if __name__ == '__main__':
    main()
