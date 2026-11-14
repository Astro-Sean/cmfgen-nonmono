#!/usr/bin/env python3
"""P8.2 spectral convergence comparison tool.

Compares obs_fin_10 spectra produced at different NUM_ES_ITERATIONS
values.  Field convergence of RJ_ES is necessary but not sufficient:
this tool measures whether the emergent spectrum itself stops changing
as sweeps are added.

Usage:
    es_sweep_compare.py FILE1 FILE2 [FILE3 ...]
                        [--labels L1 L2 ...]
                        [--out OUTDIR]
                        [--floor F] [--cont A:B,A:B,...]
                        [--lines NAME:A:B,...]
                        [--no-plots]

Each FILE is an obs_fin_10 written by cmf_flux.exe (OBSFRAME output).
Consecutive pairs (1-2, 2-3, ...) are compared; metrics are written
per pair and per wavelength window to compare_metrics.csv, and
residual-vs-sweep-count summaries plus overlays/difference spectra to
PNG files.

Checks performed per spec:
  * identical wavelength grids (mismatch -> hard error)
  * non-finite rows rejected
  * global residual R_F = max|dF| / max(|F_A|,|F_B|) (0-safe)
  * relative residual on an active mask (max(|F|) > floor * scale)
  * continuum windows: rel flux diff, max dev, integrated diff
  * line windows: core, blue-wing, red-wing, EW, centroid shift
  * global integrated spectral difference + argmax location
"""
import argparse
import csv
import re
import sys
from pathlib import Path

import numpy as np

C_AA = 2.99792458e18  # speed of light, Angstrom/s

# Default windows (rest Angstroms) chosen in line-free and line-rich
# regions typical of the SN II / IIn optical range.  All windows are
# reported in the CSV and named on the plots; override with --cont and
# --lines to pin exact windows for a given experiment.
DEFAULT_CONT = "4100:4300,5300:5600,6100:6300,7000:7400"
DEFAULT_LINES = "Hb:4700:5000,HeI5876:5800:5960,Ha:6400:6720,CaII:8400:8750"


def read_obs_fin_10(path):
    """Parse obs_fin_10 -> (wavelength Angstrom ascending, flux Jy)."""
    with open(path) as f:
        lines = f.readlines()
    nfreq = None
    for i, ln in enumerate(lines):
        if "Continuum Frequencies" in ln:
            nfreq = int(ln.split("(")[1].split(")")[0])
            start = i + 1
            break
    if nfreq is None:
        raise SystemExit(f"{path}: no 'Continuum Frequencies' header")
    vals = []
    i = start
    while len(vals) < nfreq and i < len(lines):
        vals.extend(float(x) for x in lines[i].split())
        i += 1
    if len(vals) != nfreq:
        raise SystemExit(f"{path}: short frequency block")
    freq = np.array(vals[:nfreq])
    while i < len(lines) and "Observed intensity" not in lines[i]:
        i += 1
    if i == len(lines):
        raise SystemExit(f"{path}: no 'Observed intensity' header")
    i += 1
    fvals = []
    while len(fvals) < nfreq and i < len(lines):
        fvals.extend(float(x) for x in lines[i].split())
        i += 1
    if len(fvals) != nfreq:
        raise SystemExit(f"{path}: short flux block")
    flux = np.array(fvals[:nfreq])
    lam = C_AA / freq
    order = np.argsort(lam)
    return lam[order], flux[order]


def load_set(files, labels):
    spectra = []
    for f, lab in zip(files, labels):
        lam, flx = read_obs_fin_10(f)
        bad = ~np.isfinite(flx)
        if bad.any():
            print(f"{f}: {bad.sum()} non-finite flux rows removed",
                  file=sys.stderr)
        keep = np.isfinite(flx)
        lam, flx = lam[keep], flx[keep]
        spectra.append((lab, lam, flx))
    lam0 = spectra[0][1]
    for lab, lam, _ in spectra[1:]:
        if lam.shape != lam0.shape or not np.array_equal(lam, lam0):
            raise SystemExit(
                f"wavelength grid mismatch: {lab} differs from "
                f"{spectra[0][0]} - comparison requires identical grids")
    return lam0, spectra


def global_residual(a, b):
    """max|a-b| / max(|a|,|b|) with explicit 0-scale handling."""
    amax = np.maximum(np.abs(a), np.abs(b))
    dmax = np.abs(a - b).max()
    fscale = amax.max()
    if fscale > 0:
        return dmax / fscale, dmax, fscale
    return (0.0 if dmax == 0 else 1.0), dmax, fscale


def active_rel(a, b, floor):
    """max relative change on nodes above floor * global scale."""
    amax = np.maximum(np.abs(a), np.abs(b))
    fscale = amax.max()
    mask = amax > floor * fscale
    nmask = int((~mask).sum())
    if not mask.any():
        return 0.0, 0, nmask
    rel = np.abs(a[mask] - b[mask]) / amax[mask]
    return float(rel.max()), int(mask.sum()), nmask


def window_mask(lam, lo, hi):
    return (lam >= lo) & (lam <= hi)


def continuum_metrics(lam, a, b, lo, hi):
    m = window_mask(lam, lo, hi)
    if not m.any():
        return dict(rel=np.nan, maxdev=np.nan, integ=np.nan, n=0)
    fa, fb = a[m], b[m]
    sc = np.maximum(np.abs(fa), np.abs(fb)).max()
    da = np.abs(fa - fb)
    rel = da.max() / sc if sc > 0 else np.nan
    integ = float(np.trapezoid(fb - fa, lam[m]))
    return dict(rel=rel, maxdev=float(da.max()), integ=integ,
                n=int(m.sum()))


def line_metrics(lam, a, b, lo, hi):
    m = window_mask(lam, lo, hi)
    if not m.any():
        return dict(core=np.nan, blue=np.nan, red=np.nan,
                    dew=np.nan, dcent=np.nan)
    wl, fa, fb = lam[m], a[m], b[m]
    # local continuum = straight line through the window edges
    cont_a = np.interp(wl, [wl[0], wl[-1]], [fa[0], fa[-1]])
    cont_b = np.interp(wl, [wl[0], wl[-1]], [fb[0], fb[-1]])
    sc = np.maximum(np.abs(fa), np.abs(fb)).max()
    da = np.abs(fa - fb)
    n3 = max(1, m.sum() // 3)
    core = da.max() / sc if sc > 0 else np.nan
    blue = da[:n3].max() / sc if sc > 0 else np.nan
    red = da[-n3:].max() / sc if sc > 0 else np.nan
    ew_a = np.trapezoid(1.0 - fa / np.where(cont_a == 0, np.nan, cont_a), wl)
    ew_b = np.trapezoid(1.0 - fb / np.where(cont_b == 0, np.nan, cont_b), wl)
    # centroid displacement of the absorption dip
    depth_a = cont_a - fa
    depth_b = cont_b - fb
    pos_a = np.clip(depth_a, 0, None)
    pos_b = np.clip(depth_b, 0, None)
    cent_a = (wl * pos_a).sum() / pos_a.sum() if pos_a.sum() > 0 else np.nan
    cent_b = (wl * pos_b).sum() / pos_b.sum() if pos_b.sum() > 0 else np.nan
    dcent = abs(cent_b - cent_a) if np.isfinite(cent_a + cent_b) else np.nan
    return dict(core=core, blue=blue, red=red,
                dew=float(ew_b - ew_a), dcent=dcent)


def parse_windows(spec, named):
    out = []
    for tok in spec.split(","):
        parts = tok.strip().split(":")
        if named:
            name, lo, hi = parts[0], float(parts[1]), float(parts[2])
        else:
            name, lo, hi = f"{parts[0]}-{parts[1]}", float(parts[0]), float(parts[1])
        out.append((name, lo, hi))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("files", nargs="+")
    ap.add_argument("--labels", nargs="*", default=None,
                    help="labels (e.g. sweep counts); default = file stems")
    ap.add_argument("--out", default=None,
                    help="output dir (default: alongside first file)")
    ap.add_argument("--floor", type=float, default=1.0e-6,
                    help="active-mask floor fraction of global scale")
    ap.add_argument("--cont", default=DEFAULT_CONT,
                    help="continuum windows lo:hi,lo:hi (Angstrom)")
    ap.add_argument("--lines", default=DEFAULT_LINES,
                    help="line windows name:lo:hi (Angstrom)")
    ap.add_argument("--no-plots", action="store_true")
    args = ap.parse_args()

    files = [Path(f) for f in args.files]
    labels = args.labels or [re.sub(r"\W+", "_", p.parent.name or p.stem)
                             for p in files]
    outdir = Path(args.out) if args.out else files[0].parent
    outdir.mkdir(parents=True, exist_ok=True)

    lam, spectra = load_set(files, labels)
    cont_w = parse_windows(args.cont, named=False)
    line_w = parse_windows(args.lines, named=True)

    rows = []
    pair_summ = []
    for (la, lam_a, fa), (lb, lam_b, fb) in zip(spectra, spectra[1:]):
        rg, dmax, fscale = global_residual(fa, fb)
        rrel, nact, nmask = active_rel(fa, fb, args.floor)
        argm = int(np.abs(fa - fb).argmax())
        integ = float(np.trapezoid(fb - fa, lam_a))
        pair_summ.append((la, lb, rg))
        rows.append(dict(pair=f"{la}->{lb}", metric="GLOBAL_MAX_R",
                         window="all", value=f"{rg:.6e}"))
        rows.append(dict(pair=f"{la}->{lb}", metric="GLOBAL_MAX_ABS",
                         window="all", value=f"{dmax:.6e}"))
        rows.append(dict(pair=f"{la}->{lb}", metric="GLOBAL_MAX_REL",
                         window="all", value=f"{rrel:.6e}"))
        rows.append(dict(pair=f"{la}->{lb}", metric="ARGMAX_LAM",
                         window="all", value=f"{lam_a[argm]:.2f}"))
        rows.append(dict(pair=f"{la}->{lb}", metric="INTEGRATED_DF",
                         window="all", value=f"{integ:.6e}"))
        rows.append(dict(pair=f"{la}->{lb}", metric="N_ACTIVE",
                         window="all", value=str(nact)))
        rows.append(dict(pair=f"{la}->{lb}", metric="N_MASKED",
                         window="all", value=str(nmask)))
        for name, lo, hi in cont_w:
            cm = continuum_metrics(lam, fa, fb, lo, hi)
            for key, val in cm.items():
                rows.append(dict(pair=f"{la}->{lb}",
                                 metric=f"CONT_{key.upper()}",
                                 window=name, value=f"{val:.6e}"))
        for name, lo, hi in line_w:
            lm = line_metrics(lam, fa, fb, lo, hi)
            for key, val in lm.items():
                rows.append(dict(pair=f"{la}->{lb}",
                                 metric=f"LINE_{key.upper()}",
                                 window=name, value=f"{val:.6e}"))

    csv_path = outdir / "compare_metrics.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["pair", "metric", "window", "value"])
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {csv_path}")
    for la, lb, rg in pair_summ:
        print(f"  {la} -> {lb}: GLOBAL_MAX_R = {rg:.4e}")

    if not args.no_plots:
        make_plots(lam, spectra, pair_summ, outdir, cont_w, line_w)


def make_plots(lam, spectra, pair_summ, outdir, cont_w, line_w):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels = [s[0] for s in spectra]
    xs = list(range(len(spectra)))
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(xs[1:], [p[2] for p in pair_summ], "o-", lw=1.4)
    ax.set_xticks(xs, labels)
    ax.set_xlabel("NUM_ES_ITERATIONS")
    ax.set_ylabel("GLOBAL_MAX_R vs previous run")
    ax.set_yscale("log")
    ax.set_title("Emergent-spectrum residual vs ES sweep count")
    fig.tight_layout()
    fig.savefig(outdir / "residual_vs_sweep.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 5))
    for lab, _, flx in spectra:
        ax.plot(lam, flx, lw=0.7, label=lab)
    for name, lo, hi in cont_w:
        ax.axvspan(lo, hi, color="0.9", alpha=0.5, zorder=0)
    ax.set_xlabel("Wavelength (A)")
    ax.set_ylabel("Flux (Jy)")
    ax.set_xlim(lam[lam > 3000].min() if (lam > 3000).any() else lam.min(),
                min(9500.0, lam.max()))
    ax.legend(fontsize=8, ncol=2)
    ax.set_title("Overlaid spectra across ES sweep counts")
    fig.tight_layout()
    fig.savefig(outdir / "overlaid_spectra.png", dpi=160)
    plt.close(fig)

    fig, axes = plt.subplots(len(spectra) - 1, 1, sharex=True,
                             figsize=(9, 2.2 * (len(spectra) - 1)))
    if len(spectra) - 1 == 1:
        axes = [axes]
    for ax, (la, _, fa), (lb, _, fb) in zip(axes, spectra, spectra[1:]):
        ax.plot(lam, fb - fa, lw=0.7)
        ax.axhline(0, color="k", lw=0.5)
        for name, lo, hi in line_w:
            ax.axvspan(lo, hi, color="orange", alpha=0.15, zorder=0)
        ax.set_ylabel(f"{lb}-{la}\nJy", fontsize=8)
    axes[-1].set_xlabel("Wavelength (A)")
    axes[0].set_title("Difference spectra (consecutive sweep counts)")
    fig.tight_layout()
    fig.savefig(outdir / "difference_spectra.png", dpi=160)
    plt.close(fig)
    print(f"wrote {outdir}/residual_vs_sweep.png, overlaid_spectra.png,"
          f" difference_spectra.png")


if __name__ == "__main__":
    main()
