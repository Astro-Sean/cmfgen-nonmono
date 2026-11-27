#!/usr/bin/env python3
"""JFNK admission harvest for the acceptance-path audit (checklists/jfnk.md).

Parses a CMFGEN run directory (or a bare OUTGEN) and prints:
  * the D5 before/after row: admissions, Newton attempts, commits,
    residual evals, JVPs, preconditioner calls, exit-reason histogram,
    wall time in the JFNK phase;
  * the WP-7.1 gate histogram from JFNK_GATE_REC lines: defer-class
    counts, demand vs gate, hysteresis counters, BA_LAST_RAW age,
    stable-LAMBDA relief share;
  * LINTEST operator-consistency defects (repeat/homogeneity/additivity)
    when the probe is armed;
  * JFNK_LEDGER rows (WP-0 instrumentation) when present: per-Newton-pass
    dispositions and per-activation exit reasons.

Usage:
  jfnk_harvest.py RUN_DIR [RUN_DIR ...]        summary table per run
  jfnk_harvest.py --gate RUN_DIR [...]         include gate histograms
  jfnk_harvest.py --ledger RUN_DIR             dump ledger step rows

RUN_DIR is either the .../cmfgen/<days>days work directory (containing
OUTGEN and optionally JFNK_LEDGER/ITER_LEDGER) or a direct path to an
OUTGEN file.
"""

import re
import sys
from pathlib import Path

GATE_RE = re.compile(
    r"JFNK_GATE_REC it=\s*(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+"
    r"(\d+)\s+(\d+)\s+(\d+)\s+([TF])\s+(\d+)\s+([TF])\s+(\S+)\s+(\S+)"
)
EXIT_RE = re.compile(r"solver exit reason\s*:\s*(\S+)")
NUM = {
    "fev": re.compile(r"Newton-level F evaluations:\s*(\d+)"),
    "ev": re.compile(r"Residual evals incl JVPs:\s*(\d+)"),
    "jvp": re.compile(r"Total JVP evaluations:\s*(\d+)"),
    "pc": re.compile(r"Preconditioner calls:\s*(\d+)"),
    "pcf": re.compile(r"Preconditioner fails:\s*(\d+)"),
    "nan": re.compile(r"NaN rejected steps:\s+(\d+)"),
    "t": re.compile(r"newton=\s*([\d.]+)\s+total=\s*([\d.]+)"),
}
LINTEST_RE = re.compile(
    r"(repeat|homogene|additive)[^:=]*[:=]\s*([-0-9.Ee+]+)"
)


def find_outgen(arg):
    p = Path(arg)
    if p.is_file():
        return p, p.parent / "JFNK_LEDGER"
    for cand in (p / "OUTGEN", p / "cmfgen" / "2.82days" / "OUTGEN"):
        if cand.is_file():
            return cand, cand.parent / "JFNK_LEDGER"
    raise SystemExit(f"no OUTGEN under {arg}")


def harvest(path):
    txt = path.read_text(errors="replace")
    gates = [m.groups() for m in GATE_RE.finditer(txt)]
    exits = [m.group(1) for m in EXIT_RE.finditer(txt)]
    fin = {"fev": 0, "ev": 0, "jvp": 0, "pc": 0, "pcf": 0, "nan": 0,
           "newton_s": 0.0, "total_s": 0.0, "blocks": 0,
           "newton0": 0}
    for blk in txt.split("JFNK final diagnostics:")[1:]:
        seg = blk[:1500]
        fin["blocks"] += 1
        for k in ("fev", "ev", "jvp", "pc", "pcf", "nan"):
            m = NUM[k].search(seg)
            if m:
                fin[k] += int(m.group(1))
        m = NUM["t"].search(seg)
        if m:
            tn, tt = float(m.group(1)), float(m.group(2))
            fin["newton_s"] += tn
            fin["total_s"] += tt
            fin["newton0"] += int(tn == 0.0)
    lin = [(k, float(v)) for k, v in LINTEST_RE.findall(txt)]
    return {
        "outgen": str(path),
        "admissions": txt.count("JFNK ADMIT"),
        "skipped": txt.count("admission gate skipped Newton"),
        "stag": txt.count("finite GMRES stagnation"),
        "commits": txt.count("JFNK COMMIT"),
        "guard_defers": txt.count("JFNK GUARD DEFER"),
        "catastrophic": txt.count("CATASTROPHIC eval failure"),
        "verif_rejects": txt.count("applied-state residual check"),
        "budget_stops": txt.count("budget_stop"),
        "exits": exits,
        "fin": fin,
        "gates": gates,
        "lintest": lin,
    }


def gate_histogram(gates):
    from collections import Counter

    why = Counter(g[11] for g in gates)
    n = len(gates) or 1
    out = [f"  gate calls={len(gates)}"]
    for k, c in why.most_common():
        out.append(f"    {k:<14}{c:>6}  ({100.0 * c / n:.1f}%)")
    dems = [float(g[1]) for g in gates]
    gated = [float(g[2]) for g in gates]
    admit = [g for g in gates if g[10] != "-"]
    out.append(f"  admits={len(admit)} "
               f"(frz={sum(1 for g in admit if g[10] == 'frz')}, "
               f"cpl={sum(1 for g in admit if g[10] == 'cpl')})")
    if admit:
        rel = sum(1 for g in admit if g[9] == "T")
        out.append(f"  admits under stable-LAMBDA relief: {rel}"
                   f" ({100.0 * rel / len(admit):.1f}%)")
        ages = sorted(int(g[8]) for g in admit)
        out.append(f"  BA_LAST_RAW age at admit: min={ages[0]}"
                   f" med={ages[len(ages) // 2]} max={ages[-1]}")
    if dems:
        s = sorted(dems)
        out.append(f"  demand: min={s[0]:.3g} med={s[len(s) // 2]:.3g}"
                   f" p90={s[int(0.9 * len(s))]:.3g} max={s[-1]:.3g}")
        below = sum(1 for d, g in zip(dems, gated) if d < g)
        out.append(f"  calls below frz gate: {below} "
                   f"({100.0 * below / len(dems):.1f}%)")
    hits = Counter(int(g[4]) for g in gates if g[11] == "hits_lt_nreq")
    if hits:
        out.append(f"  hits_lt_nreq frz_hits dist: {dict(sorted(hits.items()))}")
    return out


def ledger_rows(path):
    if not path or not path.is_file():
        return []
    rows = []
    for ln in path.read_text(errors="replace").splitlines():
        if ln.startswith("#") or ln.lstrip().startswith("nits"):
            continue
        if ln.strip():
            rows.append(ln.strip())
    return rows


def main():
    args = sys.argv[1:]
    do_gate = "--gate" in args
    do_led = "--ledger" in args
    args = [a for a in args if a not in ("--gate", "--ledger")]
    if not args:
        print(__doc__)
        return
    for arg in args:
        outgen, led = find_outgen(arg)
        h = harvest(outgen)
        f = h["fin"]
        print(f"== {outgen}")
        print(f"  admissions={h['admissions']} commits={h['commits']} "
              f"skipped(info=-2)={h['skipped']} stagnation={h['stag']} "
              f"guard_defers={h['guard_defers']} "
              f"verif_rejects={h['verif_rejects']} "
              f"catastrophic={h['catastrophic']}")
        print(f"  fin-diagnostic blocks={f['blocks']} fev={f['fev']} "
              f"evals={f['ev']} jvp={f['jvp']} pc={f['pc']} "
              f"pcfail={f['pcf']} nan_rej={f['nan']}")
        print(f"  newton_s={f['newton_s']:.0f} total_s={f['total_s']:.0f} "
              f"admissions_with_zero_newton={f['newton0']}")
        if h["exits"]:
            from collections import Counter
            print(f"  exit reasons: {dict(Counter(h['exits']))}")
        if h["lintest"]:
            from collections import defaultdict
            agg = defaultdict(list)
            for k, v in h["lintest"]:
                agg[k].append(v)
            for k, v in agg.items():
                v = sorted(v)
                print(f"  lintest {k}: n={len(v)} min={v[0]:.3g} "
                      f"med={v[len(v) // 2]:.3g} max={v[-1]:.3g}")
        if do_gate:
            print("\n".join(gate_histogram(h["gates"])))
        if do_led:
            rows = ledger_rows(led)
            print(f"  ledger rows={len(rows)}")
            for r in rows[-15:]:
                print(f"    {r}")
        print()


if __name__ == "__main__":
    main()
