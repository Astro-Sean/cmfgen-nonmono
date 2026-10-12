#!/usr/bin/env python3
"""
Synthesize a resolved observed-frame spectrum for a SOB (Sobolev line)
CMFGEN model from the OBS EWDATA file.

Rationale: with TRANS_*=SOB, cmf_flux computes each line's flux through
the Sobolev escape formalism (SOBEW_GRAD_V2, SIGMA_TRUE under
DO_NONMONO_SOB) and writes it to EWDATA -- it is NOT deposited into
obs_fin_10, which then contains continuum only. The EWs use the same
line-flux treatment the populations converged under, so depositing them
over each line's observed-frame resonance span produces the spectrum
consistent with the population solve. This is also the morphology in
the Dessart .fl reference spectra: modest peaks spread over a velocity
span, not the tall narrow cores the CMF/BLANK transfer emits (CMF line
transfer over-emits on non-monotonic velocity fields by ~10x because a
photon re-encounters the same resonance at each velocity reversal).

Deposition model: a shell emitting at radial speed v presents every
projected velocity in [-v, +v], so its observed profile is a box over
lam0*(1-v/c) to lam0*(1+v/c). Lines form over a range of shell speeds;
stacking boxes from slow to fast gas gives a profile peaked at rest,
approximated here by a single symmetric triangle of half-width
lam0*vdep/c. Net-absorption lines (EW<0) absorb only gas in front of
the continuum source, so their lobe covers [lam0*(1-vdep/c), lam0].
The shape is approximate -- the EW normalization is exact; the native
cmf_flux deposition path weights each depth by its line emissivity and
needs no fixed vdep. Default vdep=3000 km/s was chosen by matching the
m15 Dessart .fl peak heights for Ha, HeII 5412 and HeI 5876; a single
global value underpredicts lines forming deep in the slow CSM (HeII
4686).

Output is OBSFRAME format (same as obs_fin_10) so the existing
plot/read tooling works on it unchanged.

Usage:
    make_sob_spectrum.py EWDATA --out spec.dat [--vdep 3000]
                             [--lam-min 1000] [--lam-max 30000]
                             [--dlam 0.3]
"""

import argparse
import sys

import numpy as np

C_KMS = 2.99792458e5
C_A_S = C_KMS * 1e13  # speed of light in Angstrom/s


def read_ewdata(path):
    """Read EWDATA: Lam(Ang), C.Flux(Jy), EW(Ang) per line."""
    lam, cflux, ew = [], [], []
    with open(path) as f:
        for line in f:
            if not line.strip() or line.lstrip().startswith("!"):
                continue
            parts = line.split()
            try:
                lam.append(float(parts[0]))
                cflux.append(float(parts[2]))
                ew.append(float(parts[3]))
            except (ValueError, IndexError):
                continue
    a = np.argsort(lam)
    return np.array(lam)[a], np.array(cflux)[a], np.array(ew)[a]


def deposit_emission(lam_grid, lam0, ew_a, cflux, vdep):
    """Symmetric triangle peaked at rest over lam0*(1+-vdep/c).

    The triangle approximates the sum of per-shell boxes from slow to
    fast gas. The integral of (F-Fc)/Fc over the span equals ew_a.
    """
    half_span = lam0 * vdep / C_KMS
    m = (lam_grid >= lam0 - half_span) & (lam_grid <= lam0 + half_span)
    if m.sum() == 0:
        return None
    shape = 1.0 - np.abs(lam_grid[m] - lam0) / half_span
    shape /= shape.sum() * np.median(np.diff(lam_grid))
    return m, ew_a * cflux * shape


def deposit_absorption(lam_grid, lam0, ew_a, cflux, vdep):
    """Flat trough over [lam0*(1-vdep/c), lam0].

    Front-hemisphere gas at projected velocities in [-v, 0] absorbs the
    continuum. A box keeps the trough broad and shallow; a peaked shape
    puts the full EW in a narrow dip unlike an observed P-Cygni trough.
    """
    blue = lam0 * (1.0 - vdep / C_KMS)
    m = (lam_grid >= blue) & (lam_grid <= lam0)
    if m.sum() == 0:
        return None
    shape = np.ones(int(m.sum()))
    shape /= shape.sum() * np.median(np.diff(lam_grid))
    return m, ew_a * cflux * shape


def write_obsframe(path, lam, flux_jy):
    """Write OBSFRAME format (freq in 1e15 Hz descending, flux Jy)."""
    freq = C_A_S / lam[::-1] / 1e15
    flux = flux_jy[::-1]
    with open(path, "w") as f:
        f.write(f"\n Continuum Frequencies ( {len(freq)} )\n\n")
        for i in range(0, len(freq), 8):
            f.write("  ".join(f"{v:14.7E}" for v in freq[i : i + 8]) + "\n")
        f.write("\n Observed intensity (Janskys)\n\n")
        for i in range(0, len(flux), 4):
            f.write("  ".join(f"{v:14.7E}" for v in flux[i : i + 4]) + "\n")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("ewdata", help="EWDATA file from the OBS run")
    p.add_argument("--out", required=True, help="Output OBSFRAME spectrum")
    p.add_argument(
        "--vdep",
        type=float,
        default=3000.0,
        help="Deposition half-width km/s (default 3000, see docstring)",
    )
    p.add_argument("--lam-min", type=float, default=1000.0)
    p.add_argument("--lam-max", type=float, default=30000.0)
    p.add_argument("--dlam", type=float, default=0.3,
                   help="Output grid spacing in Angstrom")
    p.add_argument("--ew-cut", type=float, default=0.01,
                   help="Ignore lines with |EW| below this (Angstrom)")
    args = p.parse_args()

    lam_l, cflux_l, ew_l = read_ewdata(args.ewdata)
    print(f"{len(lam_l)} lines from {args.ewdata}")

    grid = np.arange(args.lam_min, args.lam_max, args.dlam)

    # Continuum: EWDATA's C.Flux sampled at line positions, interpolated
    # onto the output grid in log-log space (smooth wind continuum).
    ok = np.isfinite(cflux_l) & (cflux_l > 0)
    cont = np.exp(
        np.interp(np.log(grid), np.log(lam_l[ok]), np.log(cflux_l[ok]))
    )

    flux = cont.copy()
    n_dep = 0
    for i in range(len(lam_l)):
        if not ok[i] or abs(ew_l[i]) < args.ew_cut:
            continue
        if ew_l[i] > 0:
            res = deposit_emission(grid, lam_l[i], ew_l[i], cflux_l[i],
                                   args.vdep)
        else:
            res = deposit_absorption(grid, lam_l[i], ew_l[i], cflux_l[i],
                                     args.vdep)
        if res is None:
            continue
        m, dF = res
        flux[m] += dF
        n_dep += 1
    # Overlapping absorption lobes can drive bins below zero where the
    # deposited EW exceeds the local continuum; a negative flux is
    # unphysical and also breaks the OBSFRAME log plots.
    flux = np.maximum(flux, 1e-6 * cont)
    print(f"deposited {n_dep} lines with |EW|>={args.ew_cut} A "
          f"over {args.vdep:.0f} km/s half-width")

    write_obsframe(args.out, grid, flux)
    print(f"wrote {args.out} ({len(grid)} pts)")


if __name__ == "__main__":
    sys.exit(main())
