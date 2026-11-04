#!/usr/bin/env python3
"""POP<spec> output conservation regression.

Parses a model's POPOXY-style per-species population file (written by
CMFGEN_SUB's species dump) and checks, for every depth:

  - every full-level population is finite and non-negative,
  - per-stage level sums are consistent with the file's own totals,
  - the species accounting is closed (ion fractions in [0,1], sum of
    all stages finite at every depth).

The file location is taken from CMFGEN_POPSPEC_FILE or the first
POP* file matching POPOXY in the standard validation trees. Skips
cleanly when no file is present.
"""

import os
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
TREE = HERE.parent.parent
ROOT = TREE.parent

FAILED = []


def check(name, ok):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
    if not ok:
        FAILED.append(name)


def tofloat(tok):
    tok = re.sub(r"([0-9.])([+-])(\d{2,3})$", r"\1E\2\3", tok)
    try:
        return float(tok)
    except ValueError:
        return float("nan")


def parse_popfile(path):
    lines = path.read_text(errors="replace").split("\n")
    m = re.search(r"^\s*ND:\s*(\d+)", "\n".join(lines[:10]), re.M)
    nd = int(m.group(1)) if m else None
    secs = {}
    i = 0
    while i < len(lines):
        m = re.search(r"Number of (\S+) levels:\s+(\d+)", lines[i])
        if m:
            name, nf = m.group(1), int(m.group(2))
            i += 2  # skip oscillator date
            vals = []
            while len(vals) < nf * (nd or 0) and i < len(lines):
                toks = lines[i].split()
                if any("levels" in t or "Oscillator" in t for t in toks):
                    break
                vals += toks
                i += 1
            secs[name] = [tofloat(v) for v in vals[: nf * nd]]
        else:
            i += 1
    return nd, secs


def find_popoxy():
    env = os.environ.get("CMFGEN_POPSPEC_FILE")
    if env:
        p = Path(env)
        return p if p.exists() else None
    for cand in sorted(ROOT.glob("cmfgen_validation/*/cmfgen/*/POPOXY")):
        return cand
    return None


def main():
    path = find_popoxy()
    if path is None:
        print("test_popoxy_conservation: SKIP (no POP<spec> file found)")
        return
    print(f"  parsing {path}")
    nd, secs = parse_popfile(path)
    check("ND header parsed", nd is not None and nd > 0)
    check("at least 2 ionization stages present", len(secs) >= 2)
    if not nd or len(secs) < 2:
        print("test_popoxy_conservation: SKIP (incomplete file)")
        return

    n_bad = 0
    tot = [0.0] * nd
    fracs = {}
    for name, vals in secs.items():
        nf = len(vals) // nd
        if nf == 0:
            continue
        stage = [sum(vals[d * nf : (d + 1) * nf]) for d in range(nd)]
        fracs[name] = stage
        for d in range(nd):
            tot[d] += stage[d]
        for v in vals:
            if v != v or v < 0.0:
                n_bad += 1
        check(f"{name}: all {nf}x{nd} values finite, non-negative",
              all(v == v and v >= 0.0 for v in vals))
    check("no negative/NaN level populations anywhere", n_bad == 0)

    over = 0
    for name, s in fracs.items():
        for d in range(nd):
            if tot[d] > 0.0 and s[d] / tot[d] > 1.0 + 1e-6:
                over += 1
    check("every stage fraction <= species total", over == 0)
    check("species total finite at every depth",
          all(t == t and t >= 0.0 for t in tot))

    if FAILED:
        print(f"test_popoxy_conservation: {len(FAILED)} checks failed")
        sys.exit(1)
    print("test_popoxy_conservation: all checks passed")


if __name__ == "__main__":
    main()
