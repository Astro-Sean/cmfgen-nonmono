#!/usr/bin/env python3
"""Regression check: every cos(a,b) value printed to OUTGEN must lie in
[-1, 1] (plus float slack).

Guards the 11-Oct-2026 bug where the JJZ_DP accumulator was reused as
scratch inside the colprobe block, so the second REACH-FD print divided
by |JV_SAV| instead of |BAz|^2 and reported cos(BAz,A_FDz)=1.85e11.

Usage: verify_cosine_bounds.py <OUTGEN> [<OUTGEN> ...]
Exit 0 when all printed cosines are bounded, 1 otherwise.
"""

import re
import sys

COS_RE = re.compile(r"cos\([^)]*\)\s*=\s*([+-]?\d+\.?\d*[Ee][+-]?\d+)")
TOL = 1.0 + 1.0e-4


def main(paths):
    bad = 0
    seen = 0
    for path in paths:
        try:
            lines = open(path, errors="replace").read().splitlines()
        except OSError as exc:
            print(f"cannot read {path}: {exc}")
            bad += 1
            continue
        for ln, line in enumerate(lines, 1):
            for m in COS_RE.finditer(line):
                seen += 1
                val = float(m.group(1))
                if abs(val) > TOL:
                    bad += 1
                    print(f"{path}:{ln}: |cos| = {val:.4e}  {line.strip()}")
    print(f"checked {seen} cosine prints, {bad} out of bounds")
    return 1 if bad else 0


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(2)
    sys.exit(main(sys.argv[1:]))
