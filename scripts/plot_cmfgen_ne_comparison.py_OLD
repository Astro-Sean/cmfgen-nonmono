#!/usr/bin/env python3
"""
Plot input vs CMFGEN-final Ne and T for a converged CMFGEN model.

Reads:
  - SN_HYDRO_DATA  (ASCII: Radius grid, Temperature, Electron density sections)
  - CUR_MODEL_DATA (Fortran unformatted sequential: R, V, SIGMA, T, ED, ...)

The hydro file and the CMFGEN model can use different depth grids (the
pipeline resamples the input grid, e.g. 169 -> 150 depths).  The input
quantities are therefore interpolated (log-log in radius) onto the
CMFGEN radius grid before any ratio is computed -- comparing index by
index would pair points at different radii.

Produces a 2x2 figure:
  1. Ne vs radius: hydro input (points + interpolated) vs CMFGEN
  2. Ne ratio (CMFGEN/input) vs radius
  3. T vs radius: input vs CMFGEN
  4. T ratio (CMFGEN/input) vs radius

Left-column panels carry a secondary top axis with the CMFGEN depth
index (1 = outer boundary) for cross-referencing OUTGEN diagnostics.

Usage:
  python3 plot_cmfgen_ne_comparison.py <cmfgen_dir> [--out <path>] [--title ...]
"""

import argparse
import struct
import sys
from pathlib import Path

try:
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import ScalarFormatter, NullLocator
except ImportError as e:
    print(f"ERROR: Missing required Python package: {e}", file=sys.stderr)
    sys.exit(1)


# ---------------------------------------------------------------------------
# Readers
# ---------------------------------------------------------------------------

def read_sn_hydro(path):
    """Read radius, temperature and electron density from SN_HYDRO_DATA.

    Section lengths are taken from the 'Number of data points' header --
    never hardcode ND, the input grid length varies between models.
    """
    lines = open(path).readlines()

    nd = None
    for line in lines:
        if "Number of data points" in line:
            nd = int(line.split(":")[1])
            break

    def _parse_section(header_str):
        vals = []
        for i, line in enumerate(lines):
            if header_str in line:
                for j in range(i + 1, len(lines)):
                    try:
                        vals.extend(float(x) for x in lines[j].split())
                    except ValueError:
                        break
                break
        arr = np.array(vals)
        if nd is not None:
            arr = arr[:nd]
        return arr

    r = _parse_section("Radius grid")
    t = _parse_section("Temperature")
    ne = _parse_section("Electron density")
    return r, t, ne


def read_cur_model_data(path):
    """Read R, T and ED from CUR_MODEL_DATA (Fortran unformatted sequential).

    Record layout: date, (ND, NUM_SPECIES), SN_AGE, then R, V, SIGMA,
    T, ED, POP_ATOM, DENSITY -- each an ND-length float64 record.
    """
    with open(path, "rb") as f:
        def _read_record():
            n = struct.unpack("i", f.read(4))[0]
            data = f.read(n)
            end = struct.unpack("i", f.read(4))[0]
            if n != end:
                raise ValueError(f"Fortran record marker mismatch: {n} != {end}")
            return data

        _read_record()                              # date
        nd, _nsp = struct.unpack("ii", _read_record())
        _read_record()                              # SN_AGE
        r = np.frombuffer(_read_record(), dtype=np.float64).copy()
        _ = np.frombuffer(_read_record(), dtype=np.float64).copy()   # V
        _ = np.frombuffer(_read_record(), dtype=np.float64).copy()   # SIGMA
        t = np.frombuffer(_read_record(), dtype=np.float64).copy()
        ed = np.frombuffer(_read_record(), dtype=np.float64).copy()
        if len(r) != nd or len(t) != nd or len(ed) != nd:
            raise ValueError(f"CUR_MODEL_DATA record length != ND ({nd})")
        return r, t, ed


def interp_loglog(x_new, x_old, y_old):
    """Interpolate y_old(x_old) onto x_new in log-log space.

    Both grids are monotonically decreasing in radius -- flip to the
    increasing order np.interp expects, then flip back.
    """
    xo, yo = x_old[::-1], y_old[::-1]
    xn = x_new[::-1]
    out = np.interp(np.log(xn), np.log(xo), np.log(yo))
    return np.exp(out)[::-1]


def _set_ratio_axis(ax, ratio):
    """Linear axis near unity, log axis when the spread is large."""
    rmin, rmax = ratio.min(), ratio.max()
    ax.axhline(1.0, color="gray", linestyle="--", alpha=0.6)
    if rmax / max(rmin, 1e-30) < 4.0:
        pad = 0.05 * max(rmax - 1.0, 1.0 - rmin, 0.01)
        ax.set_ylim(min(1.0, rmin) - pad, max(1.0, rmax) + pad)
    else:
        ax.set_yscale("log")
        ax.set_ylim(max(rmin * 0.5, 1e-3), min(rmax * 2.0, 1e4))
        ax.yaxis.set_major_formatter(ScalarFormatter())
        ax.yaxis.set_minor_locator(NullLocator())


def _add_depth_axis(ax, r_c):
    """Secondary top axis: CMFGEN depth index (1 = outer boundary).

    Tick depths are kept only where adjacent labels stay separated in
    radius -- the grid compresses near shocks, and evenly spaced depth
    ticks would overprint there.
    """
    nd = len(r_c)
    sec = ax.secondary_xaxis("top")
    lo, hi = ax.get_xlim()
    candidates = np.arange(1, nd + 1, 10, dtype=int)
    if candidates[-1] != nd:
        candidates = np.append(candidates, nd)
    log_pos = np.log(r_c[np.clip(candidates - 1, 0, nd - 1)])
    span = abs(np.log(max(lo, hi)) - np.log(min(lo, hi)))
    keep = [0]
    for i in range(1, len(candidates)):
        if abs(log_pos[i] - log_pos[keep[-1]]) > 0.07 * span:
            keep.append(i)
    ticks = candidates[keep]
    pos = np.clip(r_c[ticks - 1], min(lo, hi), max(lo, hi))
    sec.set_xticks(pos)
    sec.set_xticklabels([str(t) for t in ticks])
    sec.set_xlabel("CMFGEN depth index (1 = outer)", fontsize=9)


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

def plot_comparison(r_h, t_h, ne_h, r_c, t_c, ed_c, outpath, title=""):
    # Interpolate input onto the CMFGEN grid before any ratio.
    ne_in = interp_loglog(r_c, r_h, ne_h)
    t_in = interp_loglog(r_c, r_h, t_h)

    ne_ratio = ed_c / ne_in
    t_ratio = t_c / t_in

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # --- 1. Ne vs radius -------------------------------------------------
    ax = axes[0, 0]
    ax.step(r_h, ne_h, ".", color="0.6", markersize=3, where="mid",
            label="Input hydro grid (raw)", zorder=1)
    ax.step(r_c, ne_in, "b-", linewidth=1.5, where="mid",
            label="Input interp. to CMFGEN grid", zorder=2)
    ax.step(r_c, ed_c, "r-", linewidth=1.5, where="mid",
            label="CMFGEN ED (FIX_NE=T: frozen)", zorder=3)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"Radius ($10^{10}$ cm)")
    ax.set_ylabel(r"Electron density (cm$^{-3}$)")
    ax.set_title("Electron density vs radius")
    ax.legend(fontsize=8, loc="lower left")
    ax.grid(True, alpha=0.3, which="major")
    ax.grid(True, alpha=0.1, which="minor")
    _add_depth_axis(ax, r_c)

    # --- 2. Ne ratio vs radius -------------------------------------------
    ax = axes[0, 1]
    ax.step(r_c, ne_ratio, "k-", linewidth=1.5, where="mid")
    ax.set_xscale("log")
    ax.set_xlabel(r"Radius ($10^{10}$ cm)")
    ax.set_ylabel(r"$N_e$(CMFGEN) / $N_e$(input)")
    ax.set_title("Ne ratio vs radius")
    _set_ratio_axis(ax, ne_ratio)
    ax.grid(True, alpha=0.3, which="major")

    # --- 3. T vs radius ---------------------------------------------------
    ax = axes[1, 0]
    ax.step(r_h, t_h * 1.0e4, ".", color="0.6", markersize=3, where="mid",
            label="Input hydro grid (raw)", zorder=1)
    ax.step(r_c, t_in * 1.0e4, "b-", linewidth=1.5, where="mid",
            label="Input interp. to CMFGEN grid", zorder=2)
    ax.step(r_c, t_c * 1.0e4, "r-", linewidth=1.5, where="mid",
            label="CMFGEN T (FIX_T=T: input)", zorder=3)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"Radius ($10^{10}$ cm)")
    ax.set_ylabel("Temperature (K)")
    ax.set_title("Temperature vs radius")
    ax.legend(fontsize=8, loc="lower left")
    ax.grid(True, alpha=0.3, which="major")
    ax.grid(True, alpha=0.1, which="minor")
    _add_depth_axis(ax, r_c)

    # --- 4. T ratio vs radius --------------------------------------------
    ax = axes[1, 1]
    ax.step(r_c, t_ratio, "k-", linewidth=1.5, where="mid")
    ax.set_xscale("log")
    ax.set_xlabel(r"Radius ($10^{10}$ cm)")
    ax.set_ylabel(r"$T$(CMFGEN) / $T$(input)")
    ax.set_title("T ratio vs radius")
    _set_ratio_axis(ax, t_ratio)
    ax.grid(True, alpha=0.3, which="major")

    fig.suptitle(title or "Input vs CMFGEN final: Ne and T", fontsize=14)
    fig.tight_layout()
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return outpath


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("cmfgen_dir", help="CMFGEN run directory (contains SN_HYDRO_DATA and CUR_MODEL_DATA)")
    parser.add_argument("--out", default=None, help="Output PNG path")
    parser.add_argument("--title", default="", help="Plot title")
    args = parser.parse_args()

    cmfgen_dir = Path(args.cmfgen_dir)
    hydro_file = cmfgen_dir / "SN_HYDRO_DATA"
    cur_file = cmfgen_dir / "CUR_MODEL_DATA"

    if not hydro_file.exists():
        print(f"ERROR: {hydro_file} not found", file=sys.stderr)
        sys.exit(1)
    if not cur_file.exists():
        print(f"ERROR: {cur_file} not found (CMFGEN may not have finished)", file=sys.stderr)
        sys.exit(1)

    r_h, t_h, ne_h = read_sn_hydro(hydro_file)
    r_c, t_c, ed_c = read_cur_model_data(cur_file)

    for name, arr in (("input radius", r_h), ("input Ne", ne_h), ("input T", t_h)):
        if len(arr) == 0:
            print(f"ERROR: could not parse {name} from SN_HYDRO_DATA", file=sys.stderr)
            sys.exit(1)

    if not np.all(np.diff(r_h) < 0) or not np.all(np.diff(r_c) < 0):
        print("ERROR: radius grids are not monotonically decreasing", file=sys.stderr)
        sys.exit(1)
    if r_c[0] > r_h[0] or r_c[-1] < r_h[-1]:
        print("ERROR: CMFGEN radius range outside the input grid", file=sys.stderr)
        sys.exit(1)

    print(f"  Input grid: {len(r_h)} depths;  CMFGEN grid: {len(r_c)} depths")
    if len(r_h) != len(r_c) or not np.allclose(r_h[: len(r_c)], r_c):
        print("  Grids differ: input interpolated onto the CMFGEN radius grid")

    # Default output path
    if args.out:
        outpath = args.out
    else:
        model_name = cmfgen_dir.parent.parent.name  # models/<model>/cmfgen/<epoch>
        epoch = cmfgen_dir.name
        plots_dir = cmfgen_dir / "plots"
        plots_dir.mkdir(parents=True, exist_ok=True)
        outpath = plots_dir / f"{model_name}_{epoch}_ne_comparison.png"

    title = args.title or f"{cmfgen_dir.parent.parent.name} {cmfgen_dir.name}: input vs CMFGEN"
    plot_comparison(r_h, t_h, ne_h, r_c, t_c, ed_c, str(outpath), title=title)
    print(f"Saved: {outpath}")

    # Summary
    ne_in = interp_loglog(r_c, r_h, ne_h)
    t_in = interp_loglog(r_c, r_h, t_h)
    print(f"  Input Ne:  {ne_in.min():.3e} to {ne_in.max():.3e}")
    print(f"  CMFGEN Ne: {ed_c.min():.3e} to {ed_c.max():.3e}")
    print(f"  Ne ratio:  {(ed_c / ne_in).min():.3f} to {(ed_c / ne_in).max():.3f}")
    print(f"  Input T:   {(t_in * 1e4).min():.0f} to {(t_in * 1e4).max():.0f} K")
    print(f"  CMFGEN T:  {(t_c * 1e4).min():.0f} to {(t_c * 1e4).max():.0f} K")
    print(f"  T ratio:   {(t_c / t_in).min():.3f} to {(t_c / t_in).max():.3f}")


if __name__ == "__main__":
    main()
