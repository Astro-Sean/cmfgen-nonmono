"""Command-line entry point for the X_newgrid one-pass grid builder.

Usage:
    python -m x_newgrid <swd_file> <epoch_days> <output_dir> [options]

Example:
    python -m x_newgrid X/stella/strad/M3.0_V7000.swd 35.57 \\
        /tmp/x_newgrid_test --lum 92.4

This is a STANDALONE diagnostic tool.  It does NOT submit CMFGEN jobs and
does NOT go through run_model_grid.py.  It produces a SN_HYDRO_DATA file
that can be fed into a CMFGEN run directory manually, or compared against
the v1/v2 converter outputs.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

# allow running as `python -m x_newgrid` or `python scripts/x_newgrid/cli.py`
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from grid_mapping.stella_io import read_swd, select_epoch, normalize_swd_composition, build_source_grid

try:
    from .config import XNewGridConfig
    from .builder import build_grid, build_grid_report
except ImportError:
    from x_newgrid.config import XNewGridConfig
    from x_newgrid.builder import build_grid, build_grid_report


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        description="X_newgrid: hierarchical optical-depth CMFGEN grid builder")
    p.add_argument("swd", type=Path, help="STELLA .swd file")
    p.add_argument("epoch_days", type=float, help="snapshot epoch (days)")
    p.add_argument("output_dir", type=Path, help="output directory")
    p.add_argument("--lum", type=float, default=0.0, help="luminosity (Lsun)")
    p.add_argument("--max-points", type=int, default=None)
    p.add_argument("--target-points", type=int, default=None)
    p.add_argument("--points-per-decade", type=float, default=None)
    p.add_argument("--tau-max", type=float, default=None)
    p.add_argument("--tau-min", type=float, default=None)
    p.add_argument("--kappa", type=float, default=0.34, help="e-scattering floor")
    p.add_argument("--quiet", action="store_true")
    args = p.parse_args(argv)

    cfg = XNewGridConfig(verbose=not args.quiet)
    if args.max_points is not None:
        cfg.max_points = args.max_points
    if args.target_points is not None:
        cfg.target_points = args.target_points
    if args.points_per_decade is not None:
        cfg.points_per_tau_decade = args.points_per_decade
    if args.tau_max is not None:
        cfg.tau_max = args.tau_max
    if args.tau_min is not None:
        cfg.tau_min = args.tau_min

    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Read STELLA snapshot
    data = read_swd(str(args.swd))
    best_t, ep = select_epoch(data, target_days=args.epoch_days)
    print(f"Epoch: {best_t:.3f}d, {len(ep['r'])} zones")

    mf = normalize_swd_composition(ep.get('composition'))
    source = build_source_grid(ep, best_t, mass_fractions=mf,
                               source_file=str(args.swd))

    outpath = args.output_dir / "SN_HYDRO_DATA"
    target, report = build_grid(source, cfg, outpath=str(outpath),
                                lum_lsun=args.lum)

    print(build_grid_report(report))
    with open(args.output_dir / "x_newgrid_report.json", "w") as f:
        json.dump(report, f, indent=2, default=str)
    print(f"\nReport: {args.output_dir / 'x_newgrid_report.json'}")
    print(f"Output: {outpath}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
