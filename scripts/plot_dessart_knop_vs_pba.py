#!/usr/bin/env python3
"""Compare the dessart_m15_2p82d KNOP OBSFLUX spectrum against Luc's .fl model.

Both spectra are F_lambda in erg/s/cm2/A normalized to 1 kpc.  OBSFLUX is
read in obs_fin format (Jy at 1 kpc) and converted; Luc's .fl file is already
wavelength(A) + F_lambda at 1 kpc.
"""
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("/viper/ptmp1/sbrennan/SN2024pba")
RUN = ROOT / "cmfgen_validation/m15_2p82d_knop"
LUC = ROOT / "cmfgen_validation/jnk_extracted/m15_z1p4em2_ekin1p2_mdot0p01_rcsm6e14_nb10a_vt50.fl"
# SN 2020pni +1.5 d KeckII/DEIMOS spectrum (Final Moments II flash-spectra
# sample, Zenodo 11154246). Already rest-framed (Ha peaks at 6563 A);
# UGC 9684 redshift z = 0.01687.
PNI = ROOT / "cmfgen_validation/SN2020pni_59047.32_KeckII+DEIMOS.cat"
# Continuum band shared by the observed and model spectra, used to scale
# the observed flux onto the model continuum level.
NORM_BAND_A = (5500.0, 6500.0)


def parse_fortran_float(s):
    s = re.sub(r"([0-9])([+-][0-9]{2,})$", r"\1E\2", s.strip())
    return float(s)


def read_obsflux(path):
    """obs_fin format: 'Continuum Frequencies (N)' block (1e15 Hz), then
    'Observed intensity (Janskys)' block. Stop after N values so trailing
    log text does not break parsing."""
    lines = open(path).readlines()
    h_freq = next(i for i, l in enumerate(lines) if "Continuum Frequencies" in l)
    nfreq = int(lines[h_freq].split("(")[1].split(")")[0])
    h_flux = next(i for i, l in enumerate(lines) if "Observed intensity" in l)

    freq = []
    for l in lines[h_freq + 1:h_flux]:
        freq.extend(parse_fortran_float(v) for v in l.split())
    freq = np.array(freq[:nfreq])

    flux = []
    for l in lines[h_flux + 1:]:
        for v in l.split():
            try:
                flux.append(parse_fortran_float(v))
            except ValueError:
                break
        if len(flux) >= nfreq:
            break
    flux = np.array(flux[:nfreq])

    c_cm = 2.99792458e10
    wl_A = c_cm / (freq * 1e15) * 1e8
    idx = np.argsort(wl_A)
    return wl_A[idx], flux[idx]


def jy_to_flam(wl_A, flux_jy):
    """Jy at 1 kpc -> F_lambda [erg/s/cm2/A] at 1 kpc."""
    f_nu = flux_jy * 1e-23
    lam_cm = wl_A * 1e-8
    return f_nu * 2.99792458e10 / lam_cm**2 * 1e-8


def read_luc_fl(path):
    d = np.loadtxt(path)
    m = (d[:, 0] > 0) & (d[:, 1] > 0)
    return d[m, 0], d[m, 1]


def main():
    # obs_fin_10 (cmf_flux, accurate obs-frame routine) is the definitive
    # product; OBSFLUX is the per-iteration in-loop spectrum fallback.
    spec_src = RUN / "obs" / "obs_fin_10"
    if not spec_src.is_file():
        spec_src = RUN / "OBSFLUX"
    print(f"model spectrum: {spec_src}")
    wl_m, jy_m = read_obsflux(spec_src)
    fl_m = jy_to_flam(wl_m, jy_m)
    wl_l, fl_l = read_luc_fl(LUC)

    wl_p, fl_p = read_luc_fl(PNI)  # already rest frame
    m_p = (wl_p >= NORM_BAND_A[0]) & (wl_p <= NORM_BAND_A[1]) & (fl_p > 0)
    m_l = (wl_l >= NORM_BAND_A[0]) & (wl_l <= NORM_BAND_A[1]) & (fl_l > 0)
    fl_p = fl_p * np.median(fl_l[m_l]) / np.median(fl_p[m_p])

    fig, ax = plt.subplots(figsize=(12, 6))

    ax.plot(wl_l, fl_l, color="k", lw=0.8, label="Luc m15 (converged, powered)")
    ax.plot(wl_m, fl_m, color="crimson", lw=0.8, alpha=0.85,
            label="ours m15_2p82d KNOP obs_fin_10 (iter 78, non-converged, unpowered)")
    ax.plot(wl_p, fl_p, color="steelblue", lw=0.8, alpha=0.85,
            label="SN 2020pni +1.5 d Keck/DEIMOS (continuum-scaled)")
    ax.set_yscale("log")
    ax.set_xlim(2000, 9000)
    ax.set_ylim(5e-8, 3e-4)
    ax.set_xlabel(r"Wavelength [$\AA$]")
    ax.set_ylabel(r"$F_\lambda$ @1kpc [erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    fig.tight_layout()
    out = RUN / "plots" / "knop_m15_vs_luc.svg"
    out.parent.mkdir(exist_ok=True)
    fig.savefig(out)
    fig.savefig(out.with_suffix(".png"), dpi=150)
    print(f"wrote {out}")


if __name__ == "__main__":
    sys.exit(main())
