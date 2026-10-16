#!/usr/bin/env python3
"""Compare the m15_2p82d_SOB_l14_hydrot obs_fin_10 spectrum against Luc's .fl model.

Both spectra are F_lambda in erg/s/cm2/A normalized to 1 kpc.  obs_fin_10 is
read in obs_fin format (Jy at 1 kpc) and converted; Luc's .fl file is already
wavelength(A) + F_lambda at 1 kpc.  Adapted from
plot_dessart_knop_vs_pba.py.
"""
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("/viper/ptmp1/sbrennan/SN2024pba")
RUN = ROOT / "cmfgen_validation/m15_2p82d_SOB_l14_hydrot/cmfgen/2.82days"
LUC = ROOT / "cmfgen_validation/m15_z1p4em2_ekin1p2_mdot0p01_rcsm6e14_nb10a_vt50.fl"


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
    spec_src = RUN / "obs" / "obs_fin_10"
    if not spec_src.is_file():
        spec_src = RUN / "OBSFLUX"
    print(f"model spectrum: {spec_src}")
    wl_m, jy_m = read_obsflux(spec_src)
    fl_m = jy_to_flam(wl_m, jy_m)
    # BATOBS can emit NaN/divergent points at isolated frequencies; mask so
    # they render as gaps instead of breaking the log axis.
    fl_m = np.where(np.isfinite(fl_m) & (fl_m > 0), fl_m, np.nan)
    wl_l, fl_l = read_luc_fl(LUC)

    fig, axes = plt.subplots(2, 1, figsize=(12, 10), sharex=False)

    ax = axes[0]
    ax.plot(wl_l, fl_l, color="k", lw=0.8, label="Luc m15 (converged, powered)")
    ax.plot(wl_m, fl_m, color="crimson", lw=0.8, alpha=0.85,
            label="ours m15_2p82d SOB_l14_hydrot obs_fin_10")
    ax.set_yscale("log")
    ax.set_xlim(2000, 10000)
    ax.set_ylim(1e-8, 1e-3)
    ax.set_xlabel(r"Wavelength [$\AA$]")
    ax.set_ylabel(r"$F_\lambda$ @1kpc [erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    ax.set_title("optical")

    ax = axes[1]
    ax.plot(wl_l, fl_l, color="k", lw=0.8)
    ax.plot(wl_m, fl_m, color="crimson", lw=0.8, alpha=0.85)
    ax.set_yscale("log")
    ax.set_xlim(900, 2000)
    ax.set_ylim(1e-8, 1e-3)
    ax.set_xlabel(r"Wavelength [$\AA$]")
    ax.set_ylabel(r"$F_\lambda$ @1kpc [erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]")
    ax.grid(alpha=0.3)
    ax.set_title("UV")

    fig.tight_layout()
    out = RUN / "plots" / "cmfgen_2p82d_spectrum_vs_luc.png"
    out.parent.mkdir(exist_ok=True)
    fig.savefig(out, dpi=150)
    fig.savefig(out.with_suffix(".svg"))
    print(f"wrote {out}")


if __name__ == "__main__":
    sys.exit(main())
