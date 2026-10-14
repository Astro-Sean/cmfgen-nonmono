#!/usr/bin/env python3
"""Extract per-iteration solver metrics from a CMFGEN OUTGEN file.

Phase-0 baseline instrumentation: parses the structured blocks the
non-mono development tree emits (ITER_TYPE, ITER_SUM, MODE-OUT,
BA_RES_ACC reports, NG/AV events, TUNE-independent per-iteration
timings) plus the stock "Maximum % increase/decrease" lines, and
prints one row per great iteration plus a summary of event counts.

Usage:
    iter_phase0_metrics.py OUTGEN [--csv FILE] [--summary]

The parser tolerates missing blocks: fields default to None/NaN so a
file from a tree without the instrumented blocks still parses via the
stock correction lines.
"""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass, field


RE_ITERTYPE = re.compile(
    r"ITER_TYPE it=\s*(\d+)\s+Running (\S+)"
)
RE_MAXINC = re.compile(
    r"Maximum % increase at depth\s+(\d+) is\s+(\S+)\s+"
    r"\(([^)]*)\).*--- iteration\s+(\d+)"
)
RE_MAXDEC = re.compile(
    r"Maximum % decrease at depth\s+(\d+) is\s+(\S+)\s+"
    r"\(([^)]*)\).*--- iteration\s+(\d+)"
)
RE_ITERSUM = re.compile(
    r"ITER_SUM: it=\s*(\d+)\s+mode=(\S+) -> (\S+)\s+ba=(\S+)"
    r"\s+frozen= (\S)\s+held= (\S)"
)
RE_POPS = re.compile(
    r"pops\s+(\S+)\s+(\S+)\s+\[(\w+)\]"
)
RE_SERES = re.compile(
    r"se residual\s+(\S+)\s+(\S+)\s+\[(\w+)\]"
)
RE_TRUST = re.compile(
    r"ba\s+(\S+)\s+(\d+)\s+\[(\w+)\]"
)
RE_LUM = re.compile(
    r"lum diag\s+(\S+)\s+(\S+)"
)
RE_MODEOUT = re.compile(r"MODE-OUT it\s+(\d+)")
RE_NEXT_LAM = re.compile(
    r"next lam / comp_ba / fix_t\s+:\s+(\S)\s+(\S)\s+(\S)"
)
RE_TIME = re.compile(
    r"Time taken for iteration\s+(\d+):\s+(\S+)"
)
RE_REJ = re.compile(r"BA_RES_ACC (PARTIAL )?REJ it")
RE_MASKED = re.compile(
    r"masked depth columns\s*:\s*(\d+) of\s+(\d+)"
)
RE_SOFT = re.compile(r"soft rejects=\s*(\d+)")
RE_NG_FAIL = re.compile(r"Error performing NG acceleration")
RE_AV = re.compile(r"Averaging of oscilating")
RE_LUMSTAR = re.compile(
    r"Luminosity of star \(d=1,ND\)\(iteration\s+(\d+)\) is:\s+(\S+)\s+(\S+)"
)
RE_FG_CALLS = re.compile(
    r"Average number of calls to FG_J_CMF is\s+(\S+)"
)
RE_BAMAT_RECALL = re.compile(r"BAMAT recall requested")
RE_BAMAT_RECOMPUTE = re.compile(r"recomputing LU from stored BA")
RE_NG_DONE = re.compile(
    r"NG acceleration|Max change after NG|NG_BAND", re.IGNORECASE
)


@dataclass
class IterRow:
    it: int
    iter_type: str = ""
    mode: str = ""
    mode_next: str = ""
    ba_mode: str = ""
    frozen: str = ""
    held: str = ""
    pops_demand: float = float("nan")
    pops_applied: float = float("nan")
    se_res: float = float("nan")
    se_res2: float = float("nan")
    trust: float = float("nan")
    nrej: int = 0
    lum_ratio: float = float("nan")
    next_lam: str = ""
    comp_ba_next: str = ""
    fix_t_next: str = ""
    time_s: float = float("nan")
    lum_inner: float = float("nan")
    lum_outer: float = float("nan")
    fg_calls: float = float("nan")
    max_inc_pct: float = float("nan")
    max_dec_pct: float = float("nan")
    inc_tag: str = ""
    rej_partial: bool = False
    rej_full: bool = False
    masked_cols: int = 0
    extra: list = field(default_factory=list)


def _f(text):
    try:
        return float(text.replace("D", "E").replace("d", "e"))
    except (ValueError, AttributeError):
        return float("nan")


def parse_outgen(path):
    rows = {}
    events = {
        "rej_full": 0,
        "rej_partial": 0,
        "soft_rej_max": 0,
        "ng_fail": 0,
        "ng_done": 0,
        "av_done": 0,
        "bamat_recall": 0,
        "bamat_recompute": 0,
    }
    pending_rej_partial = False
    pending_rej_full = False
    cur_it = None
    saw_modeout = False

    def row(n):
        if n not in rows:
            rows[n] = IterRow(it=n)
        return rows[n]

    with open(path, errors="replace") as fh:
        for line in fh:
            m = RE_ITERTYPE.search(line)
            if m:
                cur_it = int(m.group(1))
                row(cur_it).iter_type = m.group(2)
                continue
            m = RE_ITERSUM.search(line)
            if m:
                r = row(int(m.group(1)))
                r.mode, r.mode_next = m.group(2), m.group(3)
                r.ba_mode, r.frozen, r.held = m.group(4), m.group(5), m.group(6)
                continue
            m = RE_POPS.search(line)
            if m and cur_it is not None:
                r = row(cur_it)
                r.pops_demand, r.pops_applied = _f(m.group(1)), _f(m.group(2))
                continue
            m = RE_SERES.search(line)
            if m and cur_it is not None:
                r = row(cur_it)
                r.se_res, r.se_res2 = _f(m.group(1)), _f(m.group(2))
                continue
            m = RE_TRUST.search(line)
            if m and cur_it is not None:
                r = row(cur_it)
                r.trust = _f(m.group(1))
                r.nrej = int(m.group(2))
                continue
            m = RE_LUM.search(line)
            if m and cur_it is not None:
                row(cur_it).lum_ratio = _f(m.group(1))
                continue
            if RE_MODEOUT.search(line):
                saw_modeout = True
                continue
            m = RE_NEXT_LAM.search(line)
            if m and cur_it is not None:
                r = row(cur_it)
                r.next_lam, r.comp_ba_next, r.fix_t_next = m.groups()
                continue
            m = RE_TIME.search(line)
            if m:
                row(int(m.group(1))).time_s = _f(m.group(2))
                continue
            m = RE_MAXINC.search(line)
            if m:
                r = row(int(m.group(4)))
                r.max_inc_pct = _f(m.group(2))
                r.inc_tag = m.group(3).strip()
                continue
            m = RE_MAXDEC.search(line)
            if m:
                row(int(m.group(4))).max_dec_pct = _f(m.group(2))
                continue
            m = RE_LUMSTAR.search(line)
            if m:
                r = row(int(m.group(1)))
                r.lum_inner, r.lum_outer = _f(m.group(2)), _f(m.group(3))
                continue
            m = RE_FG_CALLS.search(line)
            if m and cur_it is not None:
                row(cur_it).fg_calls = _f(m.group(1))
                continue
            m = RE_REJ.search(line)
            if m:
                if m.group(1):
                    pending_rej_partial = True
                    events["rej_partial"] += 1
                    if cur_it is not None:
                        row(cur_it).rej_partial = True
                else:
                    pending_rej_full = True
                    events["rej_full"] += 1
                    if cur_it is not None:
                        row(cur_it).rej_full = True
                continue
            m = RE_MASKED.search(line)
            if m and pending_rej_partial and cur_it is not None:
                row(cur_it).masked_cols = int(m.group(1))
                pending_rej_partial = False
                continue
            m = RE_SOFT.search(line)
            if m:
                events["soft_rej_max"] = max(
                    events["soft_rej_max"], int(m.group(1))
                )
                continue
            if RE_NG_FAIL.search(line):
                events["ng_fail"] += 1
                continue
            if RE_AV.search(line):
                events["av_done"] += 1
                continue
            if RE_BAMAT_RECALL.search(line):
                events["bamat_recall"] += 1
                continue
            if RE_BAMAT_RECOMPUTE.search(line):
                events["bamat_recompute"] += 1
                continue
            if pending_rej_full:
                pending_rej_full = False

    # NG success is inferred from MODE-OUT "NG done ... T" flags.
    ng_re = re.compile(r"NG done / avg done / rgrid\s+:\s+(\S)")
    ng_done_iters = 0
    prev_ng_done = False
    with open(path, errors="replace") as fh:
        for line in fh:
            m = ng_re.search(line)
            if m:
                if m.group(1) == "T" and not prev_ng_done:
                    ng_done_iters += 1
                prev_ng_done = m.group(1) == "T"
    events["ng_done"] = ng_done_iters
    events["saw_modeout"] = saw_modeout
    return rows, events


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("outgen", help="Path to OUTGEN")
    ap.add_argument("--csv", help="Write per-iteration table to CSV")
    ap.add_argument("--summary", action="store_true",
                    help="Print only the event summary")
    args = ap.parse_args()

    rows, events = parse_outgen(args.outgen)
    its = sorted(rows)

    if not args.summary:
        hdr = (
            "it    type      mode            ba       dem       app"
            "     se_res    se_res2     trust  nrej  time_s"
        )
        print(hdr)
        for n in its:
            r = rows[n]
            print(
                f"{r.it:>4}  {r.iter_type or '-':<9} "
                f"{(r.mode_next or r.mode) or '-':<14} "
                f"{r.ba_mode or '-':<8} "
                f"{r.pops_demand:9.2e} {r.pops_applied:9.2e} "
                f"{r.se_res:9.2e} {r.se_res2:9.2e} "
                f"{r.trust:9.2e} {r.nrej:>4} {r.time_s:>8.1f}"
            )

    n_lam = sum(1 for n in its if "LAMBDA" in (rows[n].iter_type or ""))
    n_full = sum(1 for n in its if "FULL" in (rows[n].iter_type or ""))
    n_ba_comp = sum(
        1 for n in its if rows[n].inc_tag.startswith("BA computed")
    )
    n_ba_recall = sum(
        1 for n in its if rows[n].inc_tag.startswith("BA not computed")
    )
    total_t = sum(rows[n].time_s for n in its if rows[n].time_s == rows[n].time_s)

    print("\n--- summary ---")
    print(f"iterations parsed        : {len(its)}")
    print(f"  LAMBDA / FULL / other  : {n_lam} / {n_full} / "
          f"{len(its) - n_lam - n_full}")
    print(f"BA computed / recalled   : {n_ba_comp} / {n_ba_recall}")
    print(f"hard rejections          : {events['rej_full']}")
    print(f"partial (masked) rejects : {events['rej_partial']}")
    print(f"soft rejects (max count) : {events['soft_rej_max']}")
    print(f"NG completed / failed    : {events['ng_done']} / "
          f"{events['ng_fail']}")
    print(f"AV averaging events      : {events['av_done']}")
    print(f"BAMAT recalls / fallback recomputes: "
          f"{events['bamat_recall']} / {events['bamat_recompute']}")
    print(f"total iteration time [s] : {total_t:.1f}")

    if args.csv:
        import csv

        cols = [
            "it", "iter_type", "mode", "mode_next", "ba_mode",
            "frozen", "held", "pops_demand", "pops_applied",
            "se_res", "se_res2", "trust", "nrej", "next_lam",
            "comp_ba_next", "fix_t_next", "time_s", "lum_inner",
            "lum_outer", "lum_ratio", "fg_calls", "max_inc_pct",
            "max_dec_pct", "inc_tag", "rej_partial", "rej_full",
            "masked_cols",
        ]
        with open(args.csv, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(cols)
            for n in its:
                w.writerow([getattr(rows[n], c) for c in cols])
        print(f"wrote {args.csv}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
