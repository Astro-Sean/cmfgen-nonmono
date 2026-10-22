#!/usr/bin/env python3
"""Run the bake-off: methods x instances -> CSV of RESULT records.

Usage:
  python3 run_bakeoff.py --set TUNE --methods s0,and2,s2f [--out tune.csv]
  python3 run_bakeoff.py --set HELDOUT --methods-file methods.txt

Never run --set HELDOUT before DECISION_RULE.md is committed.  Raw
counters (sweeps, BA factorizations, cheap evals, JVPs) are recorded;
sweep-equivalent cost is computed downstream for c_ba in {0.1,.3,1}.
"""

import argparse
import csv
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
EXE = os.path.join(HERE, "bakeoff")
IDIR = os.path.join(HERE, "instances")
TIMEOUT = 300  # seconds per (method, instance) run

FIELDS = [
    "instance", "method", "sweeps", "ba", "cheap", "jvp",
    "outer", "fallback", "andrej", "status", "rmax", "rrms",
]


def instances(split_tag):
    out = []
    with open(os.path.join(IDIR, "split.txt")) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            name, tag = line.split()
            if split_tag == "ALL" or tag == split_tag:
                out.append(name)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", required=True,
                    choices=["TUNE", "HELDOUT", "ALL"])
    ap.add_argument("--methods", default="s0")
    ap.add_argument("--methods-file", default=None)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    if args.methods_file:
        with open(args.methods_file) as fh:
            methods = [l.strip() for l in fh if l.strip()
                       and not l.startswith("#")]
    else:
        methods = args.methods.split(",")

    inst = instances(args.set)
    out = args.out or os.path.join(
        HERE, f"results_{args.set.lower()}.csv")

    rows = []
    for name in inst:
        pin = os.path.join(IDIR, name + ".in")
        for m in methods:
            try:
                r = subprocess.run(
                    [EXE, m, pin], capture_output=True, text=True,
                    timeout=TIMEOUT, cwd=HERE)
            except subprocess.TimeoutExpired:
                rows.append(dict(zip(FIELDS,
                    [name, m, "", "", "", "", "", "", "",
                     "timeout", "", ""])))
                print(f"{name} {m} TIMEOUT", flush=True)
                continue
            res = [l for l in r.stdout.splitlines()
                   if l.startswith("RESULT,")]
            if not res:
                rows.append(dict(zip(FIELDS,
                    [name, m, "", "", "", "", "", "", "",
                     "crash", "", ""])))
                print(f"{name} {m} CRASH rc={r.returncode}",
                      flush=True)
                continue
            f = res[0].split(",")[1:]
            rows.append(dict(zip(FIELDS, [name] + f)))
            print(f"{name} {m} {f[8]} r={f[9].strip()}", flush=True)

    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {len(rows)} rows -> {out}")


if __name__ == "__main__":
    main()
