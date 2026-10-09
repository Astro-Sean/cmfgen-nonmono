#!/usr/bin/env python3
"""
Copy CMFGEN _OUT departure coefficient files to _IN files for a restart.

After a short bootstrap CMFGEN run (e.g. 10 iterations with LTE_EST=T),
the _OUT files contain partially-converged departure coefficients that
are much closer to the NLTE solution than a fresh LTE_EST initialization.

This script:
1. Copies all *_OUT files to *_IN files in the same directory.
2. Removes SCRTEMP, POINT1, POINT2 to force CMFGEN to read _IN files
   (rather than restarting from the scratch file).
3. Optionally removes EDDFACTOR to force recomputation of the radiation
   field, or keeps it to preserve the radiation field (USE_FIXED_J).

Usage:
    python cmfgen_restart_from_out.py <cmfgen_dir> [--keep-eddfactor]
                                       [--keep-scrtemp] [--dry-run]

Arguments:
    cmfgen_dir          Directory containing _OUT files (e.g. X_edd/cmfgen/35.57days)
    --keep-eddfactor    Keep EDDFACTOR file (use with USE_FIXED_J=T)
    --keep-scrtemp      Keep SCRTEMP/POINT1/POINT2 (restart from scratch file)
    --dry-run           Show what would be done without modifying files

The typical workflow is:

    Phase 1 (bootstrap):
        Run CMFGEN with LTE_EST=T, MAX_LAM=10, NUM_ITS=10
        → produces _OUT files with partially converged populations
        → produces EDDFACTOR with radiation field
        → produces RVTJ with full model state

    Phase 2 (this script):
        Copy _OUT → _IN
        Remove SCRTEMP/POINT1/POINT2 (force _IN read)
        Keep or remove EDDFACTOR

    Phase 3 (restart):
        Run CMFGEN with LTE_EST=F, LIN_INT=T, MAX_LAM=100
        If EDDFACTOR kept: USE_FIXED_J=T, RD_LAMBDA=T (fixed J, update pops)
        If EDDFACTOR removed: normal radiation field computation
"""

from __future__ import annotations

import argparse
import re
import shutil
import sys
from pathlib import Path


def find_out_files(cmfgen_dir: Path) -> list[Path]:
    """Find all *_OUT files in the directory."""
    out_files = sorted(cmfgen_dir.glob("*OUT"))
    # Exclude files like "SCRTEMP_OUT" or other non-DC outputs
    # CMFGEN DC output files are named like HIOUT, HeIOUT, Fe2OUT, etc.
    # They should have corresponding species in MODEL_SPEC
    dc_out_files = []
    for f in out_files:
        # Skip non-DC output files
        if f.name in ("SCRTEMP_OUT", "RVTJ_OUT", "EDDFACTOR_OUT"):
            continue
        dc_out_files.append(f)
    return dc_out_files


def out_to_in_name(out_file: Path) -> Path:
    """Convert _OUT filename to _IN filename.

    CMFGEN naming convention:
        HIOUT  -> HI_IN
        HeIOUT -> HeI_IN
        Fe2OUT -> Fe2_IN
        FeIIIOUT -> FeIII_IN
    """
    name = out_file.name
    if name.endswith("OUT"):
        base = name[:-3]
        return out_file.parent / f"{base}_IN"
    return out_file


def remove_restart_files(cmfgen_dir: Path, dry_run: bool = False) -> list[Path]:
    """Remove SCRTEMP, POINT1, POINT2 to force _IN file reading."""
    removed = []
    for fname in ["SCRTEMP", "POINT1", "POINT2", "MODEL_SCR"]:
        fpath = cmfgen_dir / fname
        if fpath.exists():
            removed.append(fpath)
            if not dry_run:
                fpath.unlink()
    return removed


def handle_eddfactor(
    cmfgen_dir: Path, keep: bool, dry_run: bool = False
) -> str:
    """Handle EDDFACTOR file."""
    edd = cmfgen_dir / "EDDFACTOR"
    edd_info = cmfgen_dir / "EDDFACTOR_INFO"
    if keep:
        if edd.exists() and edd_info.exists():
            return "kept"
        return "not_present"
    removed = []
    for fpath in [edd, edd_info]:
        if fpath.exists():
            removed.append(fpath)
            if not dry_run:
                fpath.unlink()
    return "removed" if removed else "not_present"


def verify_out_file_format(out_file: Path) -> bool:
    """Quick check that an _OUT file has the expected format."""
    try:
        with open(out_file) as f:
            lines = [f.readline() for _ in range(5)]
        # First non-empty line should contain a date or format date
        for line in lines:
            if "Format date" in line or "FEB" in line or "Dec" in line or "Jan" in line:
                return True
            if line.strip() and not line.strip().startswith("!"):
                # Could be the R(ND), LUM, NHYD, ND line
                parts = line.split()
                if len(parts) >= 4:
                    try:
                        float(parts[0])
                        float(parts[1])
                        int(parts[2])
                        int(parts[3])
                        return True
                    except (ValueError, IndexError):
                        pass
        return True  # Don't block on format check
    except Exception:
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Copy CMFGEN _OUT files to _IN files for restart"
    )
    parser.add_argument(
        "cmfgen_dir",
        type=Path,
        help="Directory containing _OUT files",
    )
    parser.add_argument(
        "--keep-eddfactor",
        action="store_true",
        help="Keep EDDFACTOR (use with USE_FIXED_J=T)",
    )
    parser.add_argument(
        "--keep-scrtemp",
        action="store_true",
        help="Keep SCRTEMP/POINT1/POINT2 (restart from scratch file)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be done without modifying files",
    )
    args = parser.parse_args()

    cmfgen_dir = args.cmfgen_dir.resolve()
    if not cmfgen_dir.is_dir():
        print(f"ERROR: {cmfgen_dir} is not a directory", file=sys.stderr)
        sys.exit(1)

    # Find _OUT files
    out_files = find_out_files(cmfgen_dir)
    if not out_files:
        print(f"ERROR: No *_OUT files found in {cmfgen_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"CMFGEN restart from _OUT files")
    print(f"  Directory: {cmfgen_dir}")
    print(f"  Found {len(out_files)} _OUT files:")
    for f in out_files:
        in_name = out_to_in_name(f)
        print(f"    {f.name} -> {in_name.name}")
        if not verify_out_file_format(f):
            print(f"      WARNING: format check failed for {f.name}")

    if args.dry_run:
        print("\n[DRY RUN - no files modified]")
    else:
        # Copy _OUT -> _IN
        copied = 0
        for out_file in out_files:
            in_file = out_to_in_name(out_file)
            shutil.copy2(out_file, in_file)
            copied += 1
        print(f"\nCopied {copied} _OUT -> _IN files")

    # Handle restart files
    if not args.keep_scrtemp:
        removed = remove_restart_files(cmfgen_dir, args.dry_run)
        if removed:
            action = "Would remove" if args.dry_run else "Removed"
            print(f"{action} restart files:")
            for f in removed:
                print(f"  {f.name}")
        else:
            print("No restart files (SCRTEMP/POINT1/POINT2) found")
    else:
        print("Keeping SCRTEMP/POINT1/POINT2 (restart from scratch)")

    # Handle EDDFACTOR
    edd_status = handle_eddfactor(cmfgen_dir, args.keep_eddfactor, args.dry_run)
    if args.keep_eddfactor:
        print(f"EDDFACTOR: {edd_status}")
        if edd_status == "kept":
            print("  → Use USE_FIXED_J=T, RD_LAMBDA=T in VADAT for fixed J restart")
    else:
        print(f"EDDFACTOR: {edd_status}")
        if edd_status == "removed":
            print("  → Radiation field will be recomputed on restart")

    # Print restart instructions
    print("\nRestart instructions:")
    if args.keep_eddfactor and edd_status == "kept":
        print("  1. Set in VADAT: LTE_EST=F, LIN_INT=T, USE_FIXED_J=T")
        print("  2. Set in IN_ITS: DO_LAM_IT=T (Lambda iterations)")
        print("  3. Run CMFGEN - it will read _IN files + fixed J from EDDFACTOR")
    else:
        print("  1. Set in VADAT: LTE_EST=F, LIN_INT=T")
        print("  2. Set in IN_ITS: DO_LAM_IT=T (Lambda iterations)")
        print("  3. Run CMFGEN - it will read _IN files + recompute J")

    print("\nDone.")


if __name__ == "__main__":
    main()
