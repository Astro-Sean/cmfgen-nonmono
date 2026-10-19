#!/usr/bin/env python3
"""Compare m15_2p82d CMFGEN output observed spectra against Luc's .fl.

Only post-processed obs_fin products (output/*_obs_fin.dat, Jy at 1 kpc)
are plotted, converted to F_lambda [erg/s/cm2/A].  Luc's .fl file is
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
VAL = ROOT / "cmfgen_validation"
LUC = VAL / "m15_z1p4em2_ekin1p2_mdot0p01_rcsm6e14_nb10a_vt50.fl"

# (label, color, path, linestyle, zorder)
# Only post-processed output spectra (output/*obs_fin.dat) - files
# produced by the batobs/cmf_flux chain.  In-code OBSFLUX writes are
# excluded: they are diagnostic dumps at whatever iteration wrote them,
# not the observed-spectrum product.
SPECS = [
    (
        "SOB l14 lte (conv exit, 88%)",
        "crimson",
        VAL
        / "m15_2p82d_SOB_l14_lte/cmfgen/2.82days/output"
        / "dessart_m15_2p82d_SOB_l14_lte_2.82days_obs_fin.dat",
        "-",
        6,
    ),
    (
        "SOB l14 jfnk (unconverged)",
        "tab:blue",
        VAL
        / "m15_2p82d_SOB_l14_jfnk/cmfgen/2.82days/output"
        / "dessart_m15_2p82d_SOB_l14_jfnk_2.82days_obs_fin.dat",
        "--",
        5,
    ),
]

OUTDIR = VAL / "plots"


def parse_fortran_float(s):
    s = re.sub(r"([0-9])([+-][0-9]{2,})$", r"\1E\2", s.strip())
    return float(s)


def read_obsflux(path):
    """OBSFLUX/obs_fin format: 'Continuum Frequencies (N)' block (1e15 Hz),
    then 'Observed intensity (Janskys)' block. Stop after N values."""
    lines = open(path).readlines()
    h_freq = next(i for i, l in enumerate(lines) if "Continuum Frequencies" in l)
    nfreq = int(lines[h_freq].split("(")[1].split(")")[0])
    h_flux = next(i for i, l in enumerate(lines) if "Observed intensity" in l)

    freq = []
    for l in lines[h_freq + 1 : h_flux]:
        freq.extend(parse_fortran_float(v) for v in l.split())
    freq = np.array(freq[:nfreq])

    flux = []
    for l in lines[h_flux + 1 :]:
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
    curves = []
    for label, color, path, ls, zo in SPECS:
        if not path.is_file():
            print(f"skip (missing): {path}")
            continue
        wl, jy = read_obsflux(path)
        fl = jy_to_flam(wl, jy)
        fl = np.where(np.isfinite(fl) & (fl > 0), fl, np.nan)
        curves.append((label, color, wl, fl, ls, zo))
        print(f"{label}: {len(wl)} pts  {wl.min():.0f}-{wl.max():.0f} A")

    wl_l, fl_l = read_luc_fl(LUC)

    fig, axes = plt.subplots(3, 1, figsize=(13, 14), sharex=False)

    def draw(ax, lo, hi, title):
        ax.plot(wl_l, fl_l, color="k", lw=1.0, zorder=1,
                label="Luc m15 (converged, powered)")
        for label, color, wl, fl, ls, zo in curves:
            ax.plot(wl, fl, color=color, lw=0.7, ls=ls, alpha=0.85,
                    zorder=zo, label=label)
        ax.set_yscale("log")
        ax.set_xlim(lo, hi)
        ax.set_ylim(1e-8, 1e-3)
        ax.set_xlabel(r"Wavelength [$\AA$]")
        ax.set_ylabel(r"$F_\lambda$ @1kpc [erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]")
        ax.grid(alpha=0.3)
        ax.set_title(title)
        ax.legend(fontsize=8, loc="lower right")

    draw(axes[0], 2000, 10000, "optical")
    draw(axes[1], 900, 2000, "UV")
    draw(axes[2], 4500, 7000, "optical zoom (C IV 5800s, N II, He II 6560)")

    fig.tight_layout()
    OUTDIR.mkdir(exist_ok=True)
    out = OUTDIR / "cmfgen_2p82d_models_vs_luc.png"
    fig.savefig(out, dpi=150)
    fig.savefig(out.with_suffix(".svg"))
    print(f"wrote {out}")


if __name__ == "__main__":
    sys.exit(main())
