#!/usr/bin/env python3
"""Aggregate bake-off RESULT CSVs into per-method statistics.

Sweep-equivalents = sweeps + c_ba * ba_factorizations.  Cheap evals and
JVPs are reported raw; they do not add sweeps (a frozen JVP shares the
eval's transfer operator on the toy, matching the production frozen
tier where the JVP is a cheap linear response).

Usage: python3 analyze.py results_tune.csv [results_heldout.csv ...]
"""

import csv
import statistics
import sys

C_BA = (0.1, 0.3, 1.0)


def load(path):
    rows = []
    with open(path) as fh:
        for r in csv.DictReader(fh):
            for k in ("sweeps", "ba", "cheap", "jvp", "outer",
                      "fallback", "andrej"):
                r[k] = int(r[k]) if r[k] else 0
            r["rmax"] = float(r["rmax"]) if r["rmax"] else float("nan")
            rows.append(r)
    return rows


def main():
    rows = []
    for p in sys.argv[1:]:
        rows += load(p)
    methods = sorted({r["method"] for r in rows})
    insts = sorted({r["instance"] for r in rows})
    n_inst = len(insts)

    print(f"instances={n_inst} methods={len(methods)} rows={len(rows)}")
    hdr = (
        f"{'method':<14} {'ok':>5} {'rate':>5}"
        + "".join(f" {'med@'+str(c):>8}" for c in C_BA)
        + f" {'worst@0.3':>9} {'medcheap':>8} {'medjvp':>6}"
        + "  statuses"
    )
    print(hdr)
    print("-" * len(hdr))
    for m in methods:
        rs = [r for r in rows if r["method"] == m]
        ok = [r for r in rs if r["status"] == "converged"]
        med, worst = [], []
        for c in C_BA:
            costs = [r["sweeps"] + c * r["ba"] for r in ok]
            med.append(statistics.median(costs) if costs else float("nan"))
        if ok:
            worst = max(r["sweeps"] + 0.3 * r["ba"] for r in ok)
        st = {}
        for r in rs:
            st[r["status"]] = st.get(r["status"], 0) + 1
        sts = ",".join(f"{k}:{v}" for k, v in sorted(st.items()))
        print(
            f"{m:<14} {len(ok):>3}/{n_inst} {len(ok)/n_inst:>5.2f}"
            + "".join(f" {x:>8.1f}" for x in med)
            + f" {worst:>9.1f}"
            f" {statistics.median([r['cheap'] for r in rs]):>8.0f}"
            f" {statistics.median([r['jvp'] for r in rs]):>6.0f}"
            + f"  {sts}"
        )

    # per-instance matrix for the two reference methods
    print("\nper-instance sweep-equivalents @ c_ba=0.3 (blank = failed)")
    sel = [m for m in methods]
    print(f"{'instance':<14}" + "".join(f" {m[:9]:>9}" for m in sel))
    for inst in insts:
        line = f"{inst:<14}"
        for m in sel:
            r = next(
                (x for x in rows
                 if x["instance"] == inst and x["method"] == m),
                None,
            )
            if r is None or r["status"] != "converged":
                line += f" {'--':>9}"
            else:
                line += f" {r['sweeps'] + 0.3*r['ba']:>9.1f}"
        print(line)


if __name__ == "__main__":
    main()
