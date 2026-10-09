#!/usr/bin/env python3
"""Plot Luc's m15 spectrum with the fresh USE_LAM_ES=T and =F runs.

A model spectrum that is missing, still being written, or too short to
parse is skipped. Luc's spectrum is plotted on its own if neither model
file is ready.
"""
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("/viper/ptmp1/sbrennan/SN2024pba")
LUC = ROOT / "cmfgen_validation/jnk_extracted/m15_z1p4em2_ekin1p2_mdot0p01_rcsm6e14_nb10a_vt50.fl"
RUNS = (
    ("LAM", ROOT / "cmfgen_validation/m15_2p82d_knop_LAM", "C0"),
    ("MOM", ROOT / "cmfgen_validation/m15_2p82d_knop_MOM", "C3"),
)
OUT = ROOT / "cmfgen_validation/plots/dessart_luc_lam_mom.png"
MIN_POINTS = 20
XLIM_A = (2000.0, 9000.0)


def parse_fortran_float(token):
    token = re.sub(r"([0-9])([+-][0-9]{2,})$", r"\1E\2", token.strip())
    return float(token)


def read_obsflux(path):
    """obs_fin / OBSFLUX: frequency block in 1e15 Hz, then flux in Jy."""
    lines = path.read_text(errors="replace").splitlines()
    h_freq = next(i for i, line in enumerate(lines) if "Continuum Frequencies" in line)
    nfreq = int(lines[h_freq].split("(")[1].split(")")[0])
    h_flux = next(i for i, line in enumerate(lines) if "Observed intensity" in line)
    if nfreq < MIN_POINTS:
        raise ValueError(f"only {nfreq} frequencies")

    freq = []
    for line in lines[h_freq + 1:h_flux]:
        freq.extend(parse_fortran_float(v) for v in line.split())
    if len(freq) < nfreq:
        raise ValueError(f"frequency block {len(freq)}/{nfreq}")
    freq = np.asarray(freq[:nfreq], dtype=float)

    flux = []
    for line in lines[h_flux + 1:]:
        for token in line.split():
            try:
                flux.append(parse_fortran_float(token))
            except ValueError:
                break
        if len(flux) >= nfreq:
            break
    if len(flux) < nfreq:
        raise ValueError(f"flux block {len(flux)}/{nfreq}")
    flux = np.asarray(flux[:nfreq], dtype=float)
    if not np.isfinite(freq).all() or not np.isfinite(flux).all():
        raise ValueError("non-finite wavelength or flux")

    c_cm = 2.99792458e10
    wl_A = c_cm / (freq * 1e15) * 1e8
    order = np.argsort(wl_A)
    return wl_A[order], flux[order]


def jy_to_flam(wl_A, flux_jy):
    """Jy at 1 kpc to F_lambda [erg/s/cm2/A] at 1 kpc."""
    f_nu = flux_jy * 1e-23
    lam_cm = wl_A * 1e-8
    return f_nu * 2.99792458e10 / lam_cm**2 * 1e-8


def read_luc_fl(path):
    data = np.loadtxt(path)
    keep = (data[:, 0] > 0) & (data[:, 1] > 0) & np.isfinite(data[:, 0]) & np.isfinite(data[:, 1])
    return data[keep, 0], data[keep, 1]


def model_spectrum(run_dir):
    """Prefer the finished obs-frame file. Fall back to the in-loop OBSFLUX."""
    candidates = (run_dir / "obs" / "obs_fin_10", run_dir / "OBSFLUX")
    for path in candidates:
        if not path.is_file() or path.stat().st_size == 0:
            continue
        try:
            wl_A, flux_jy = read_obsflux(path)
        except (OSError, ValueError, IndexError, StopIteration) as exc:
            print(f"skip {path}: {exc}")
            continue
        return path, wl_A, jy_to_flam(wl_A, flux_jy)
    print(f"skip {run_dir.name}: no readable spectrum yet")
    return None


def flux_in_window(wl_A, fl):
    keep = (
        (wl_A >= XLIM_A[0]) & (wl_A <= XLIM_A[1])
        & (fl > 0) & np.isfinite(fl)
    )
    return fl[keep]


def set_ylim_to_data(ax, samples):
    """Log limits from the plotted flux inside the wavelength window."""
    positive = np.concatenate([s for s in samples if s.size])
    if positive.size == 0:
        return
    lo = float(np.min(positive))
    hi = float(np.max(positive))
    if hi <= lo:
        hi = lo * 10.0
    span = np.log10(hi / lo)
    pad = max(0.08 * span, 0.08)
    ax.set_ylim(lo / 10.0**pad, hi * 10.0**pad)


def main():
    if not LUC.is_file():
        print(f"Luc spectrum not found: {LUC}", file=sys.stderr)
        return 1
    wl_luc, fl_luc = read_luc_fl(LUC)
    samples = [flux_in_window(wl_luc, fl_luc)]

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(wl_luc, fl_luc, color="k", lw=0.8, label="Luc m15")
    n_model = 0
    for label, run_dir, color in RUNS:
        loaded = model_spectrum(run_dir)
        if loaded is None:
            continue
        path, wl_A, fl = loaded
        ax.plot(wl_A, fl, color=color, lw=0.8, alpha=0.85, label=f"{label} ({path.name})")
        samples.append(flux_in_window(wl_A, fl))
        print(f"plotted {label}: {path}")
        n_model += 1

    ax.set_yscale("log")
    ax.set_xlim(*XLIM_A)
    set_ylim_to_data(ax, samples)
    ax.set_xlabel(r"Wavelength [$\AA$]")
    ax.set_ylabel(r"$F_\lambda$ at 1 kpc [erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]")
    ax.legend(fontsize=9, frameon=False)
    ax.grid(alpha=0.3)
    fig.tight_layout()

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=150)
    print(f"wrote {OUT} ({n_model} model spectra)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
