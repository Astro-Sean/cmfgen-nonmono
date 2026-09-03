#!/usr/bin/env python3
"""
Plot input Ne (from SN_HYDRO_DATA) vs CMFGEN final Ne (from CUR_MODEL_DATA).

Reads:
  - SN_HYDRO_DATA  (ASCII, electron density section)
  - CUR_MODEL_DATA (Fortran unformatted binary, ED field)

Produces a 4-panel figure:
  1. Ne vs depth (input vs CMFGEN)
  2. Ne ratio (CMFGEN/input) vs depth
  3. Ne vs radius
  4. Ne ratio vs radius

Usage:
  python3 plot_cmfgen_ne_comparison.py <cmfgen_dir> [--out <path>]

  <cmfgen_dir>  Directory containing SN_HYDRO_DATA and CUR_MODEL_DATA
  --out         Output PNG path (default: <cmfgen_dir>/plots/<model>_ne_comparison.png)
"""

import argparse
import os
import struct
import sys
from pathlib import Path

try:
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError as e:
    print(f"ERROR: Missing required Python package: {e}", file=sys.stderr)
    sys.exit(1)


# ---------------------------------------------------------------------------
# Readers
# ---------------------------------------------------------------------------

def read_sn_hydro_ne(path):
    """Read electron density from SN_HYDRO_DATA (6 values per line, 150 total)."""
    with open(path, "r") as f:
        lines = f.readlines()

    def _parse_section(header_str, n=150):
        vals = []
        for i, line in enumerate(lines):
            if header_str in line:
                for j in range(i + 1, len(lines)):
                    try:
                        row = [float(x) for x in lines[j].split()]
                        if len(row) == 6:
                            vals.extend(row)
                            if len(vals) >= n:
                                return np.array(vals[:n])
                    except (ValueError, IndexError):
                        break
                break
        return np.array(vals[:n]) if vals else np.array([])

    ne = _parse_section("Electron density")
    r = _parse_section("Radius grid")
    return r, ne


def read_cur_model_data(path):
    """Read R and ED from CUR_MODEL_DATA (Fortran unformatted sequential)."""
    with open(path, "rb") as f:
        def _read_record():
            n = struct.unpack("i", f.read(4))[0]
            data = f.read(n)
            end = struct.unpack("i", f.read(4))[0]
            if n != end:
                raise ValueError(f"Fortran record marker mismatch: {n} != {end}")
            return data

        # Record 1: date string
        _read_record()
        # Record 2: ND, NUM_SPECIES
        rec = _read_record()
        nd, _nsp = struct.unpack("ii", rec)
        # Record 3: SN_AGE
        _read_record()
        # Records 4-10: R, V, SIGMA, T, ED, POP_ATOM, DENSITY
        r = np.frombuffer(_read_record(), dtype=np.float64).copy()
        _ = np.frombuffer(_read_record(), dtype=np.float64).copy()  # V
        _ = np.frombuffer(_read_record(), dtype=np.float64).copy()  # SIGMA
        _ = np.frombuffer(_read_record(), dtype=np.float64).copy()  # T
        ed = np.frombuffer(_read_record(), dtype=np.float64).copy()
        return r, ed


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

def plot_ne_comparison(r_input, ne_input, r_cmfgen, ed_cmfgen, outpath, title=""):
    nd = len(ne_input)
    depth = np.arange(1, nd + 1)
    ratio = ed_cmfgen / ne_input

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # 1. Ne vs depth
    ax = axes[0, 0]
    ax.semilogy(depth, ne_input, "b-", label="Input (SN_HYDRO_DATA)", alpha=0.7, linewidth=1.5)
    ax.semilogy(depth, ed_cmfgen, "r-", label="CMFGEN final iteration", alpha=0.7, linewidth=1.5)
    ax.axvline(x=nd // 2, color="gray", linestyle="--", alpha=0.5)
    ax.set_xlabel("Depth (1=outer, {}=inner)".format(nd))
    ax.set_ylabel(r"Electron density (cm$^{-3}$)")
    ax.set_title("Electron Density: Input vs CMFGEN Final")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # 2. Ne ratio vs depth
    ax = axes[0, 1]
    ax.semilogy(depth, ratio, "k-", linewidth=1.5)
    ax.axhline(y=1.0, color="gray", linestyle="--", alpha=0.5)
    ax.set_xlabel("Depth")
    ax.set_ylabel(r"Ne(CMFGEN) / Ne(input)")
    ax.set_title("Ne Ratio (CMFGEN / Input)")
    ax.grid(True, alpha=0.3)
    lo = max(ratio.min() * 0.5, 1e-3)
    hi = min(ratio.max() * 2.0, 1e3)
    ax.set_ylim(lo, hi)

    # 3. Ne vs radius
    ax = axes[1, 0]
    ax.loglog(r_input, ne_input, "b-", label="Input", alpha=0.7, linewidth=1.5)
    ax.loglog(r_cmfgen, ed_cmfgen, "r-", label="CMFGEN final", alpha=0.7, linewidth=1.5)
    ax.set_xlabel(r"Radius ($10^{10}$ cm)")
    ax.set_ylabel(r"Electron density (cm$^{-3}$)")
    ax.set_title("Ne vs Radius")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # 4. Ne ratio vs radius
    ax = axes[1, 1]
    ax.semilogy(r_input, ratio, "k-", linewidth=1.5)
    ax.axhline(y=1.0, color="gray", linestyle="--", alpha=0.5)
    ax.set_xlabel(r"Radius ($10^{10}$ cm)")
    ax.set_ylabel(r"Ne(CMFGEN) / Ne(input)")
    ax.set_title("Ne Ratio vs Radius")
    ax.grid(True, alpha=0.3)
    ax.set_ylim(lo, hi)

    fig.suptitle(title or "Input Ne vs CMFGEN Final Ne", fontsize=14, y=1.01)
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

    r_input, ne_input = read_sn_hydro_ne(hydro_file)
    r_cmfgen, ed_cmfgen = read_cur_model_data(cur_file)

    if len(ne_input) == 0:
        print("ERROR: could not parse electron density from SN_HYDRO_DATA", file=sys.stderr)
        sys.exit(1)
    if len(ed_cmfgen) == 0:
        print("ERROR: could not parse ED from CUR_MODEL_DATA", file=sys.stderr)
        sys.exit(1)
    if len(ne_input) != len(ed_cmfgen):
        print(f"WARNING: length mismatch (input={len(ne_input)}, cmfgen={len(ed_cmfgen)}); truncating to min")
        n = min(len(ne_input), len(ed_cmfgen))
        ne_input, ed_cmfgen = ne_input[:n], ed_cmfgen[:n]
        r_input, r_cmfgen = r_input[:n], r_cmfgen[:n]

    # Default output path
    if args.out:
        outpath = args.out
    else:
        model_name = cmfgen_dir.parent.parent.name  # models/<model>/cmfgen/<epoch>
        epoch = cmfgen_dir.name
        plots_dir = cmfgen_dir / "plots"
        plots_dir.mkdir(parents=True, exist_ok=True)
        outpath = plots_dir / f"{model_name}_{epoch}_ne_comparison.png"

    title = args.title or f"{cmfgen_dir.parent.parent.name} {cmfgen_dir.name}: Input Ne vs CMFGEN Final Ne"
    plot_ne_comparison(r_input, ne_input, r_cmfgen, ed_cmfgen, str(outpath), title=title)
    print(f"Saved: {outpath}")

    # Summary
    ratio = ed_cmfgen / ne_input
    print(f"  Input Ne:  {ne_input.min():.3e} to {ne_input.max():.3e}")
    print(f"  CMFGEN Ne: {ed_cmfgen.min():.3e} to {ed_cmfgen.max():.3e}")
    print(f"  Ratio:     {ratio.min():.3f} to {ratio.max():.3f}")


if __name__ == "__main__":
    main()
