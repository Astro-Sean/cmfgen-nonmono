#!/usr/bin/env python3
"""check_convergence — compact status report from a KNOP ITER_HISTORY file.

Usage:
    check_convergence.py [RUN_DIR_OR_ITER_HISTORY] [--tail N] [--watch SEC]
                         [--plot [FILE]]

Accepts either the ITER_HISTORY file itself or a CMFGEN run directory
containing one — run it from inside a run folder with no arguments.
If an OUTGEN sits alongside, the ITER_MODE reason strings are scanned for
controller events (lum_drift_hold, release_damped, gs_freeze_forced, mode
transitions).

--plot writes a convergence-history figure next to ITER_HISTORY
(default FILE: CONVERGENCE_HISTORY.PNG).  Uses matplotlib when available;
otherwise falls back to a pure-stdlib SVG (FILE.SVG).  Under --watch the
figure is regenerated every cycle so it tracks the live run.  Watch mode
stops on its own once the stop-converge criterion is met (the latest
demanded % corrections fall below EPS_TERM): it writes a final plot and
exits with 'MODEL CONVERGED'.  It also stops if the run goes quiet --
no writes to ITER_HISTORY/OUTGEN for --stale seconds (default 1800) --
exiting with 'RUN APPEARS STOPPED'.

Pure stdlib — no numpy needed.

ITER_HISTORY columns (written by SOLVE_FOR_POPS):
    iter mode next ba froz held maxraw maxapp maxatt se se2 trust nrej lum
    fg_d fg_badf fg_churnf fg_cap mom_d mom_n5f mom_n50f mom_negjf
    mom_quarf mom_cap mom_pquar
"""

import argparse
import math
import re
import sys
import time
from pathlib import Path

SATURATED_RES = 1.0e30      # BA_RES_SAT — degenerate inner-boundary rows
LAM_EXIT_DEFAULT = 100.0    # used when VADAT has no LAM_EXIT key
ZERO_CORR_TOL = 1.0e-12     # demanded correction ~0 = vetoed step, not conv.
# Demand values pinned at this level are the MODCTL record cap (HELD
# rows report maxraw=1e7 regardless of the true demand); they carry no
# rate information for the contraction estimate.
DEMAND_CAP = 9.0e6
# BA states where nothing was committed, so the row must not count as
# progress evidence.
BA_VETOED = ("MASKED", "REJECT")
BA_NONCOMMIT = ("LAMBDA", "HELD") + BA_VETOED

# Bounded OUTGEN read: the head carries the control echo (EPS_TERM,
# VADAT), the tail carries the recent iterations the watch cares about.
_OUTGEN_HEAD = 8 * 1024 * 1024
_OUTGEN_TAIL = 256 * 1024 * 1024
_OUTGEN_CACHE = {}

# ITER_HISTORY mode/ba columns -> marker style on the demand panel.
_MODE_STYLE = {
    "LAMBDA_GS":   ("o", "tab:orange", "LAMBDA (GS freeze)"),
    "LAMBDA_TFIX": ("o", "tab:blue",   "LAMBDA (T fixed)"),
    "LAMBDA":      ("o", "tab:cyan",   "LAMBDA"),
    "FULL_TFIX":   ("s", "tab:green",  "FULL (T fixed)"),
    "FULL":        ("D", "tab:red",    "FULL"),
}
# BA-state overlay on the demanded series: marker, edge color, label.
# Hollow markers keep the mode marker underneath readable.
_BA_STYLE = {
    "ARMED":  ("*", "tab:purple", 11, "BA armed (unjudged step)"),
    "HELD":   ("s", "k",          8,  "BA held"),
    "MASKED": ("x", "tab:orange", 8,  "BA masked"),
    "REJECT": ("X", "tab:red",    9,  "BA rejected"),
}
_SVG_MODE_COLOR = {
    "LAMBDA_GS": "#ff7f0e", "LAMBDA_TFIX": "#1f77b4", "LAMBDA": "#17becf",
    "FULL_TFIX": "#2ca02c", "FULL": "#d62728",
}
# tab:* / "k" are matplotlib names, not CSS; map ba states explicitly
_SVG_BA_COLOR = {
    "ARMED": "#9467bd", "HELD": "#000000",
    "MASKED": "#ff7f0e", "REJECT": "#d62728",
}


def _outgen_lines(hist_path):
    """OUTGEN as a list of lines; [] when absent.

    Whole file when small, otherwise head + tail (the VADAT/EPS echo is
    at the top, the iterations of interest at the bottom).  Cached on
    (size, mtime) so a watch cycle reads the file at most once no matter
    how many scanners run.
    """
    p = Path(hist_path).parent / "OUTGEN"
    try:
        st = p.stat()
    except OSError:
        return []
    key = str(p)
    ent = _OUTGEN_CACHE.get(key)
    if ent and ent[0] == st.st_size and ent[1] == st.st_mtime_ns:
        return ent[2]
    if st.st_size <= _OUTGEN_HEAD + _OUTGEN_TAIL:
        text = p.read_text(errors="replace")
    else:
        with open(p, "rb") as f:
            head = f.read(_OUTGEN_HEAD)
            f.seek(st.st_size - _OUTGEN_TAIL)
            tail = f.read()
        text = (head.decode(errors="replace") + "\n"
                + tail.decode(errors="replace"))
    lines = text.splitlines()
    _OUTGEN_CACHE.clear()
    _OUTGEN_CACHE[key] = (st.st_size, st.st_mtime_ns, lines)
    return lines


def fnum(s):
    """Parse a Fortran ES value; tolerate the dropped 'E' (1.5465+262)."""
    try:
        return float(s)
    except ValueError:
        m = re.match(r"^(.*?)([+-]\d+)$", s)
        if m:
            return float(m.group(1) + "E" + m.group(2))
        return float("nan")


def vadat_number(run_dir, key, default):
    path = run_dir / "VADAT"
    if not path.is_file():
        return default
    marker = f"[{key}]"
    for line in path.read_text(errors="replace").splitlines():
        if marker in line:
            try:
                return fnum(line.split()[0].replace("D", "E").replace("d", "e"))
            except (IndexError, ValueError):
                return default
    return default


def load_history(path):
    cols, rows = None, []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith("#"):
            cols = line[1:].split()
            continue
        toks = line.split()
        if not toks[0].isdigit():
            continue
        rows.append(toks)
    if cols is None:
        cols = ("iter mode next ba froz held maxraw maxapp maxatt se se2 "
                "trust nrej lum fg_d fg_badf fg_churnf fg_cap").split()
    idx = {c: i for i, c in enumerate(cols)}
    return idx, rows


def col(row, idx, name, kind="f"):
    i = idx.get(name)
    if i is None or i >= len(row):
        return float("nan") if kind == "f" else ""
    v = row[i]
    return fnum(v) if kind == "f" else v


def geomean(vals):
    vals = [v for v in vals if v and v > 0]
    return math.exp(sum(math.log(v) for v in vals) / len(vals)) if vals else float("nan")


def _seg_start(rows, idx):
    """Index where the last segment begins.

    A restarted run appends to ITER_HISTORY with the iteration counter
    reset to 1, so a row whose iter is <= the previous row's starts a new
    segment.  Trend stats and milestones only make sense within the
    latest segment (OUTGEN also only reflects it).
    """
    seg0 = 0
    prev = None
    for i, r in enumerate(rows):
        v = int(float(col(r, idx, "iter")))
        if prev is not None and v <= prev:
            seg0 = i
        prev = v
    return seg0


def report(hist_path, tail, lines=None):
    idx, rows = load_history(hist_path)
    if not rows:
        print(f"no data rows in {hist_path}")
        return 1

    n = len(rows)
    last = rows[-1]
    it = int(float(col(last, idx, "iter")))
    mode = col(last, idx, "mode", "s")
    nxt = col(last, idx, "next", "s")
    ba = col(last, idx, "ba", "s")
    froz = col(last, idx, "froz", "s") == "T"
    held = col(last, idx, "held", "s") == "T"
    raw = col(last, idx, "maxraw")
    se = col(last, idx, "se")
    se2 = col(last, idx, "se2")
    nrej = col(last, idx, "nrej")
    lum = col(last, idx, "lum")
    fgd = col(last, idx, "fg_d")
    fgbad = col(last, idx, "fg_badf")
    fgcap = col(last, idx, "fg_cap")

    lums = [col(r, idx, "lum") for r in rows]
    seg0 = _seg_start(rows, idx)
    seg_rows = rows[seg0:]
    lam_exit = vadat_number(hist_path.parent, "LAM_EXIT", LAM_EXIT_DEFAULT)

    print(f"file    : {hist_path}")
    nresets = sum(1 for i in range(1, n)
                  if int(float(col(rows[i], idx, "iter")))
                  <= int(float(col(rows[i - 1], idx, "iter"))))
    print(f"iters   : {n}   latest: it{it}  {mode} -> {nxt}  ba={ba}  "
          f"frozen={froz} held={held}"
          + (f"   [{nresets} restart(s), current segment has "
             f"{len(seg_rows)} iters]" if seg0 else ""))

    # ---- modes seen -----------------------------------------------------
    from collections import Counter
    mc = Counter(col(r, idx, "mode", "s") for r in rows)
    bac = Counter(col(r, idx, "ba", "s") for r in rows)
    nfroz = sum(1 for r in rows if col(r, idx, "froz", "s") == "T")
    nheld = sum(1 for r in rows if col(r, idx, "held", "s") == "T")
    print(f"modes   : {dict(mc)}")
    print(f"ba states: {dict(bac)}   frozen {nfroz}/{n}   held {nheld}/{n}")

    # ---- pops demand trend ----------------------------------------------
    # stay inside the last segment: crossing a restart makes the ratio junk
    win = rows[max(seg0, n - tail):]
    wba = Counter(col(r, idx, "ba", "s") for r in win)
    # Masked/rejected steps commit nothing and rows whose demand is pinned
    # at the MODCTL record cap (HELD reports maxraw=1e7 flat) carry no rate
    # information.  Excluding them keeps the contraction and ETA honest:
    # they describe the iterates that actually moved.
    wr_pairs = [(col(r, idx, "maxraw"), col(r, idx, "ba", "s")) for r in win]
    wr = [v for v, b in wr_pairs
          if math.isfinite(v) and v > 0
          and b not in BA_VETOED and v < DEMAND_CAP]
    nskip = len(wr_pairs) - len(wr)
    if len(wr) > 1 and all(math.isfinite(v) and v > 0 for v in wr):
        rate = geomean([wr[i + 1] / wr[i] for i in range(len(wr) - 1)])
        note = ""
        if nskip:
            note = (f"  [{nskip}/{len(wr_pairs)} rows excluded: "
                    "vetoed/cap-pinned]")
        if held:
            note = "  [HELD — demand from field conv. only, ETA not meaningful]"
        print(f"maxraw  : {raw:.3e}   window-{tail} contraction {rate:.2f}x/iter"
              + note)
        if 0 < rate < 1 and raw > lam_exit and not held:
            eta = math.log(lam_exit / raw) / math.log(rate)
            print(f"          ETA to maxraw<{lam_exit:.0f}: ~{eta:.0f} its "
                  f"(~{eta * 5.5 / 60:.1f} h at 330 s/iter)")
    else:
        print(f"maxraw  : {raw:.3e}"
              + (f"   [window ba: {dict(wba)}]" if wba else ""))

    # ---- luminosity ------------------------------------------------------
    dl = abs(lums[-1] - lums[-2]) if len(lums) > 1 else float("nan")
    scale = max(abs(lums[-1]), abs(lums[-2]), 0.1) if len(lums) > 1 else float("nan")
    rel_dl = dl / scale if math.isfinite(dl) else float("nan")
    lmax = max(abs(v) for v in lums if math.isfinite(v))
    exc = [i + 1 for i, v in enumerate(lums)
           if math.isfinite(v) and i > 0 and abs(v - lums[i - 1]) > 3 * max(abs(v), 1)]
    print(f"lum     : {lum:+.3e}   rel_drift={rel_dl:.3e}  diagnostic-only  "
          f"peak|ratio|={lmax:.2e}")
    if exc:
        print(f"          excursion iters: {exc[:8]}{'...' if len(exc) > 8 else ''}")

    # ---- field maturity ---------------------------------------------------
    print(f"fg      : d={fgd:.2e}  bad={fgbad * 100:.3f}%  cap={fgcap:.0f}   "
          f"{'MATURE' if fgbad < 1e-4 and fgd < 5e-3 else 'immature'}")

    # ---- residuals --------------------------------------------------------
    sat = se > SATURATED_RES
    print(f"se_res  : {se:.3e} {'(saturated channel — judged on res2)' if sat else ''}"
          f"   res2={se2:.3e}")

    # ---- rejects -----------------------------------------------------------
    trej = sum(int(col(r, idx, "nrej")) for r in rows)
    print(f"rejects : nrej now {int(nrej)}, cumulative {trej}")

    # ---- wall time ---------------------------------------------------------
    tt = scan_outgen_time(hist_path, lines)
    if tt:
        tot = sum(v for _, v in tt)
        unit, scale = _time_unit(tot)
        print(f"wtime   : total {tot / scale:.2f} {unit} over {len(tt)} iters"
              f"   last iter {tt[-1][1]:.0f} s")

    # ---- milestones (within the current segment) -------------------------
    seg_iter = lambda r: int(float(col(r, idx, "iter")))
    first_unfrozen = next((seg_iter(r) for r in seg_rows
                           if col(r, idx, "froz", "s") == "F"), None)
    first_full = next((seg_iter(r) for r in seg_rows
                       if col(r, idx, "ba", "s") not in ("LAMBDA", "")), None)
    first_mature = next((seg_iter(r) for r in seg_rows
                         if math.isfinite(col(r, idx, "fg_d"))
                         and col(r, idx, "fg_d") < 5e-3), None)
    ms = []
    if seg0:
        ms.append(f"iteration counter reset at row {seg0 + 1} "
                  f"(it{seg_iter(seg_rows[0])})")
    if first_mature:
        ms.append(f"field matured ~it{first_mature}")
    if first_unfrozen:
        ms.append(f"pops released it{first_unfrozen}")
    if first_full:
        ms.append(f"first BA step it{first_full}")
    print("milestones: " + ("; ".join(ms) if ms else "none yet"))

    # ---- last N rows --------------------------------------------------------
    print(f"\nlast {min(tail, n)} rows:")
    fmt = ("{:>4} {:<12.12} {:<12.12} {:<7.7} {:>6} {:>5}"
           " {:>9} {:>9} {:>9} {:>9}")
    hdr = fmt.format("iter", "mode", "next_mode", "ba", "frozen", "held",
                     "maxraw", "maxapp", "lum_ratio", "fg_drift")
    print(hdr)
    print("-" * len(hdr))
    prev_it = None
    for r in rows[-tail:]:
        itv = int(float(col(r, idx, "iter")))
        if prev_it is not None and itv <= prev_it:
            print("    " + "~" * 10 + " restart " + "~" * 10)
        prev_it = itv
        print(fmt.format(
            int(float(col(r, idx, "iter"))),
            col(r, idx, "mode", "s"), col(r, idx, "next", "s"),
            col(r, idx, "ba", "s"), col(r, idx, "froz", "s"),
            col(r, idx, "held", "s"),
            f"{col(r, idx, 'maxraw'):.2e}", f"{col(r, idx, 'maxapp'):.2e}",
            f"{col(r, idx, 'lum'):.2e}", f"{col(r, idx, 'fg_d'):.2e}"))

    return 0


def scan_outgen(hist_path, lines=None):
    if lines is None:
        lines = _outgen_lines(hist_path)
    if not lines:
        return
    from collections import Counter
    reasons = Counter()
    transitions = 0
    for line in lines:
        if "ITER_MODE" in line:
            transitions += 1
            m = re.search(r"reason=(\S+)", line)
            if m:
                for t in m.group(1).split("+"):
                    reasons[t] += 1
    if transitions:
        print(f"\nOUTGEN reasons ({transitions} mode lines):")
        for k, v in reasons.most_common():
            print(f"  {k:<24} {v}")


def scan_outgen_chg(hist_path, lines=None):
    """Per-iteration 'Maximum % increase/decrease' diagnostics from OUTGEN.

    Returns two dicts {iteration: (value, depth)}.  These are percent
    changes reported by the population solver, distinct units from the
    fractional maxraw/maxapp columns of ITER_HISTORY.
    """
    inc, dec = {}, {}
    if lines is None:
        lines = _outgen_lines(hist_path)
    if not lines:
        return inc, dec
    pat = re.compile(
        r"Maximum % (increase|decrease) at depth\s+(\d+)\s+is\s+(\S+)"
        r".*?iteration\s+(\d+)")
    for line in lines:
        m = pat.search(line)
        if m:
            tgt = inc if m.group(1) == "increase" else dec
            tgt[int(m.group(4))] = (fnum(m.group(3)), int(m.group(2)))
    return inc, dec


def scan_outgen_solveba(hist_path, lines=None):
    values = {}
    if lines is None:
        lines = _outgen_lines(hist_path)
    if not lines:
        return values
    iter_pat = re.compile(r"--- iteration\s+(\d+)")
    maxch_pat = re.compile(
        r"Maximm changes as returned by SOLVEBA_V13 is\s+(\S+)")
    current_iter = None
    for line in lines:
        iter_match = iter_pat.search(line)
        if iter_match:
            current_iter = int(iter_match.group(1))
        maxch_match = maxch_pat.search(line)
        if maxch_match and current_iter is not None:
            values[current_iter] = fnum(maxch_match.group(1))
    return values


def scan_outgen_eps(hist_path, lines=None):
    """Effective EPS_TERM from OUTGEN's control-echo block.

    The model terminates when MAXCH (the maximum fractional % change)
    falls below EPS -- i.e. EPS is in the same percent units as the
    correction series on panel 1.  Returns None if not found.
    """
    if lines is None:
        lines = _outgen_lines(hist_path)
    if not lines:
        return None
    num = r"([-+0-9.]+[dDeE][-+0-9]+|[-+0-9.]+)"
    override = re.compile(rf"EPS_TERM.*?override\s*=\s*{num}")
    change = re.compile(rf"EPS_TERM[^=]*?{num}\s*->\s*{num}")
    bare = re.compile(rf"EPS_TERM[^=]*?{num}")
    eps = None
    for line in lines:
        if "EPS_TERM" not in line:
            continue
        m = override.search(line) or change.search(line) or bare.search(line)
        if m:
            eps = float(m.group(m.lastindex).replace("D", "E")
                        .replace("d", "e"))
    return eps


def converged_status(hist_path, lines=None):
    """(converged, maxch, eps) for the latest ITER_HISTORY row.

    The stop criterion is MAXCH < EPS_TERM: the largest demanded percent
    correction (max % increase, and the % decrease expressed in
    demand-equivalent units) below EPS_TERM.  maxraw is included since it
    tracks the same demanded-correction maximum.  Returns converged=False
    when EPS_TERM could not be read from OUTGEN.
    """
    idx, rows = load_history(hist_path)
    if not rows:
        return False, float("nan"), None
    eps = scan_outgen_eps(hist_path, lines)
    if eps is None:
        return False, float("nan"), None
    it = int(float(col(rows[-1], idx, "iter")))
    inc, dec = scan_outgen_chg(hist_path, lines)
    cands = [col(rows[-1], idx, "maxraw")]
    if it in inc and _fin(inc[it][0]):
        cands.append(inc[it][0])
    if it in dec and _fin(dec[it][0]):
        v = dec[it][0]
        cands.append(1.0e7 if v >= 99.999
                     else (100.0 * v / (100.0 - v) if v > 0 else 0.0))
    finite = [v for v in cands if _fin(v)]
    maxch = max(finite) if finite else float("nan")
    # False-convergence veto (m15_2p82d_jfnk_hyb, job 12073582): a
    # BA_RES_ACC rejection masks/zeros every depth correction, so the
    # demanded change reports as exactly 0 — that is a vetoed step, not
    # a converged model.  Require a nonzero demanded correction, and
    # treat a latest iteration whose BA step was masked or counted
    # consecutive rejections as not converged outright.
    last_ba = col(rows[-1], idx, "ba", "s").upper()
    last_nrej = col(rows[-1], idx, "nrej")
    rejected = (last_ba == "MASKED"
                or (_fin(last_nrej) and last_nrej > 0))
    converged = (_fin(maxch) and ZERO_CORR_TOL < maxch < eps
                 and not rejected)
    return converged, maxch, eps


def _run_signature(hist_path):
    """(mtime_ns, size) of the files the model writes (ITER_HISTORY + OUTGEN).

    Deliberately not a whole-directory scan: the watch's own plot output
    would keep the directory looking alive forever.  None when
    ITER_HISTORY is absent.
    """
    try:
        h = hist_path.stat()
    except OSError:
        return None
    sig = [(h.st_mtime_ns, h.st_size)]
    try:
        o = (hist_path.parent / "OUTGEN").stat()
        sig.append((o.st_mtime_ns, o.st_size))
    except OSError:
        pass
    return tuple(sig)


def scan_outgen_time(hist_path, lines=None):
    """Per-iteration wall times from OUTGEN 'Time taken for iteration' lines.

    Returns a sorted list of (iteration, seconds).
    """
    if lines is None:
        lines = _outgen_lines(hist_path)
    if not lines:
        return []
    pat = re.compile(
        r"Time taken for iteration\s+(\d+)\s*:\s*([-+0-9.]+[dDeE][-+0-9]+|"
        r"[-+0-9.]+)")
    tt = {}
    for line in lines:
        m = pat.search(line)
        if m:
            tt[int(m.group(1))] = float(
                m.group(2).replace("D", "E").replace("d", "e"))
    return sorted(tt.items())


def _time_unit(total_s):
    """Smallest unit keeping the total under ~200 of itself."""
    for unit, scale in (("days", 86400.0), ("hours", 3600.0),
                        ("min", 60.0)):
        if total_s / scale >= 2.0:
            return unit, scale
    return "s", 1.0


# ---------------------------------------------------------------------------
# plotting
# ---------------------------------------------------------------------------

def _extract(rows, idx):
    """Pull the plottable columns out of ITER_HISTORY rows.

    x = row position (a monotonic axis immune to iteration-counter
    resets); rst = x positions where a new run segment begins.
    """
    g = lambda name, kind="f": [col(r, idx, name, kind) for r in rows]
    iters = [int(float(v)) for v in g("iter")]
    return {
        "it": iters, "x": list(range(len(iters))),
        "rst": [i for i in range(1, len(iters)) if iters[i] <= iters[i - 1]],
        "mode": g("mode", "s"), "nxt": g("next", "s"), "ba": g("ba", "s"),
        "froz": [v == "T" for v in g("froz", "s")],
        "held": [v == "T" for v in g("held", "s")],
        "raw": g("maxraw"), "app": g("maxapp"),
        "se": g("se"), "se2": g("se2"),
        "lum": g("lum"), "fgd": g("fg_d"), "fgb": g("fg_badf"),
        "nrej": g("nrej"),
    }


def _spans(flags, xs):
    """Contiguous (x0, x1) ranges where flags is True."""
    spans, start, last = [], None, None
    for x, f in zip(xs, flags):
        if f:
            if start is None:
                start = x
            last = x
        elif start is not None:
            spans.append((start, last))
            start = None
    if start is not None:
        spans.append((start, last))
    return spans


def _milestones(d):
    """(x position, label) pairs — positions are row indices, not iters."""
    seg0 = d["rst"][-1] if d["rst"] else 0
    ms = []
    for name, pred in (
            ("Pops released", lambda i: not d["froz"][i]),
            ("First BA step", lambda i: d["ba"][i] not in ("LAMBDA", "")),
            ("Field mature", lambda i: math.isfinite(d["fgd"][i]) and d["fgd"][i] < 5e-3)):
        x = next((i for i in range(seg0, len(d["it"])) if pred(i)), None)
        if x is not None:
            ms.append((x, name))
    return ms


def _fin(v):
    return math.isfinite(v)


def _lum_is_big(d):
    """True when the luminosity ratio spans >50x in magnitude.

    Near-unity values (the common converging case) read better on a
    linear axis; symlog is only needed when shocks/deposition push the
    ratio through orders of magnitude or negative excursions.
    """
    mags = [abs(v) for v in d["lum"] if _fin(v) and v != 0]
    return bool(mags) and max(mags) / min(mags) > 50


def make_plot(hist_path, out_path, lines=None):
    """Write a convergence-history figure. Returns the path written."""
    idx, rows = load_history(hist_path)
    if not rows:
        print("no data rows — nothing to plot")
        return None
    if lines is None:
        lines = _outgen_lines(hist_path)
    d = _extract(rows, idx)
    inc_map, dec_map = scan_outgen_chg(hist_path, lines)
    solveba_map = scan_outgen_solveba(hist_path, lines)
    nan = float("nan")
    # OUTGEN series are keyed by iteration number, which repeats across
    # restarted segments — attach each value to the LAST row carrying
    # that iter, since OUTGEN describes the newest segment.
    last_pos = {v: i for i, v in enumerate(d["it"])}
    keyed = lambda m: [m[v][0] if v in m and i == last_pos[v] else nan
                       for i, v in enumerate(d["it"])]
    d["pinc"] = keyed(inc_map)
    d["pdec"] = keyed(dec_map)
    d["solveba"] = [
        solveba_map[v] if v in solveba_map and i == last_pos[v] else nan
        for i, v in enumerate(d["it"])]
    # Demand-equivalent of the decrease channel, using the same transform as
    # solveba_v13.f: a demanded decrease of d% reports as 100*d/(100-d),
    # saturated at 1e7.  Overlaps maxraw exactly while the pin is active.
    d["pdec_eq"] = [
        min(1.0e7, 100.0 * v / (100.0 - v))
        if _fin(v) and 0 < v < 99.999
        else (1.0e7 if _fin(v) and v >= 99.999 else nan)
        for v in d["pdec"]]
    d["eps"] = scan_outgen_eps(hist_path, lines)
    d["lam_exit"] = vadat_number(hist_path.parent, "LAM_EXIT", LAM_EXIT_DEFAULT)
    tt_map = dict(scan_outgen_time(hist_path, lines))
    d["tt"] = [(i, tt_map[v]) for i, v in enumerate(d["it"])
               if v in tt_map and i == last_pos[v]]
    ms = _milestones(d)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        _plot_mpl(plt, d, ms, hist_path, out_path)
        return out_path
    except ImportError:
        svg = out_path.with_suffix(".SVG")
        _plot_svg(d, ms, hist_path, svg)
        return svg


def _plot_mpl(plt, d, ms, hist_path, out_path):
    xs = d["x"]            # monotonic row positions — the plot x-axis
    it = d["it"]           # actual iteration numbers, used for labels
    tt = d.get("tt") or []
    npan = 5 if tt else 4
    fig, axes = plt.subplots(npan, 1, figsize=(10, 2.2 * (npan + 1)),
                             sharex=True,
                             gridspec_kw=dict(hspace=0.14,
                                              height_ratios=[2] + [1] * (npan - 1)))
    fig.align_ylabels(axes)

    from matplotlib import ticker
    pos = lambda vs: [v if (_fin(v) and v > 0) else float("nan") for v in vs]
    htags = []   # (ax, y, label, linestyle) — drawn after the first draw
    TBBOX = dict(boxstyle="round,pad=0.12", facecolor="white",
                 edgecolor="none", alpha=0.85)

    def _legend(ax, *a, **k):
        # borderless white halo: text stays readable over data, no heavy box
        leg = ax.legend(*a, frameon=True, framealpha=0.85,
                        edgecolor="none", **k)
        leg.set_zorder(12)
        return leg

    # Frozen / held shading on every panel.  The x-limits are pinned to
    # one half-cell beyond the first/last iteration so the vspan edges
    # coincide exactly with the axis limits on every panel.
    for ax in axes:
        for x0, x1 in _spans(d["froz"], xs):
            ax.axvspan(x0 - 0.45, x1 + 0.45, color="orange", alpha=0.13, lw=0)
        for x0, x1 in _spans(d["held"], xs):
            ax.axvspan(x0 - 0.45, x1 + 0.45, color="tab:blue", alpha=0.07, lw=0)
        for r in d["rst"]:
            ax.axvline(r - 0.5, ls=":", c="tab:red", alpha=0.6, lw=1.2)
        ax.set_xlim(-0.5, len(xs) - 0.5)
    axes[0].plot([], [], color="orange", alpha=0.4, lw=6,
                 label="Frozen iterations")
    axes[0].plot([], [], color="tab:blue", alpha=0.4, lw=6,
                 label="Populations held")

    def _empty_note(ax, has_data):
        # a silent empty panel reads as "converged"; say why it is blank
        if not has_data:
            ax.text(0.5, 0.5, "no data on this solver path",
                    transform=ax.transAxes, ha="center", va="center",
                    fontsize=9, color="0.5", style="italic")

    def _mode_scatter(ax, xs, modes, vs, logpos=True):
        """Per-mode markers over a series."""
        ok = lambda v: _fin(v) and (v > 0 or not logpos)
        put = ax.semilogy if logpos else ax.plot
        for mode, (mk, c, lab) in _MODE_STYLE.items():
            pts = [(x, v) for x, m, v in zip(xs, modes, vs)
                   if m == mode and ok(v)]
            if pts:
                put(*zip(*pts), marker=mk, ms=4.5, ls="none", c=c,
                    label=lab)
        # unknown/future modes still get marked, neutrally
        pts = [(x, v) for x, m, v in zip(xs, modes, vs)
               if m not in _MODE_STYLE and ok(v)]
        if pts:
            put(*zip(*pts), "o", ms=4.5, ls="none", c="tab:blue")

    def _full_frame(ax, modes, vs, logpos=True):
        """Hollow black square at every full-iteration x on a panel.

        Frames the series value where it exists; where the panel has no
        usable point at that x the square is pinned to the bottom edge.
        Callers skip this entirely on panels with no data at all.
        """
        ok = lambda v: _fin(v) and (v > 0 or not logpos)
        pts = [(x, v) for x, m, v in zip(xs, modes, vs)
               if m.startswith("FULL")]
        hit = [(x, v) for x, v in pts if ok(v)]
        if hit:
            (ax.semilogy if logpos else ax.plot)(
                *zip(*hit), "s", ms=11, ls="none", mfc="none",
                mec="k", mew=1.2, label="Full iteration")
        miss = [x for x, v in pts if not ok(v)]
        if miss:
            ax.plot(miss, [0.05] * len(miss), "s", ms=8, ls="none",
                    mfc="none", mec="k", mew=1.2,
                    transform=ax.get_xaxis_transform())

    # 1: pops demand
    ax = axes[0]
    ax.semilogy(xs, pos(d["raw"]), "-", lw=1.2, c="tab:blue",
                label="Demanded correction (maxraw)")
    _mode_scatter(ax, xs, d["mode"], d["raw"])
    _full_frame(ax, d["mode"], d["raw"])
    ax.semilogy(xs, pos(d["app"]), "s--", ms=3, c="0.45", alpha=0.7,
                label="Applied correction (maxapp)")
    if any(_fin(v) for v in d["pinc"]):
        ax.semilogy(xs, pos(d["pinc"]), "^-", ms=4, c="tab:green",
                    alpha=0.75, label="Max % increase")
    if any(_fin(v) for v in d["pdec_eq"]):
        ax.semilogy(xs, pos(d["pdec_eq"]), "v-", ms=4, c="tab:red",
                    alpha=0.75,
                    label="Max decrease demand (as correction)")
    if any(_fin(v) for v in d["solveba"]):
        ax.semilogy(xs, pos(d["solveba"]), "x--", ms=5, c="tab:purple",
                    alpha=0.85, label="SOLVEBA_V13 returned max")
    for ba, (mk, c, msz, lab) in _BA_STYLE.items():
        pts = [(x, v) for x, b, v in zip(xs, d["ba"], d["raw"])
               if b == ba and _fin(v) and v > 0]
        if pts:
            ax.semilogy(*zip(*pts), marker=mk, ms=msz, ls="none",
                        mfc="none", mec=c, mew=1.4, label=lab)
    htags.append((ax, 1.0e7, "Reporting cap (decrease ≥ 99.999%)", ":"))
    htags.append((ax, d["lam_exit"], f"BA exit threshold {d['lam_exit']:.0f}", ":"))
    if d.get("eps") is not None:
        htags.append((ax, d["eps"],
                      f"Stop converge (EPS_TERM={d['eps']:g}%)", "--"))
    ax.set_ylabel("Max population correction")
    # Iteration-type / BA-state labels are carried by the luminosity
    # panel's legend; keep only the data series here.
    _marker_labs = ({lab for _, _, lab in _MODE_STYLE.values()}
                    | {t[3] for t in _BA_STYLE.values()}
                    | {"Full iteration"})
    _hh, _ll = ax.get_legend_handles_labels()
    _da_h = [h for h, l in zip(_hh, _ll) if l not in _marker_labs]
    _da_l = [l for l in _ll if l not in _marker_labs]
    _legend(ax, _da_h, _da_l, loc="upper center", fontsize=7, ncol=3)

    # 2: luminosity ratio — symlog absorbs the occasional negative excursion
    ax = axes[1]
    has_lum = any(_fin(v) for v in d["lum"])
    ax.plot(xs, d["lum"], "-", lw=1.2, c="tab:green",
            label="Luminosity ratio")
    if has_lum:
        _mode_scatter(ax, xs, d["mode"], d["lum"], logpos=False)
        _full_frame(ax, d["mode"], d["lum"], logpos=False)
    # BA-state glyphs on the luminosity series too: this panel carries
    # the iteration-type legend, and excursions are easier to interpret
    # when the solver state that produced them is visible
    for ba, (mk, c, msz, lab) in _BA_STYLE.items():
        pts = [(x, v) for x, b, v in zip(xs, d["ba"], d["lum"])
               if b == ba and _fin(v)]
        if pts:
            ax.plot(*zip(*pts), marker=mk, ms=msz, ls="none",
                    mfc="none", mec=c, mew=1.4, label=lab)
    _empty_note(ax, has_lum)
    if _lum_is_big(d):
        ax.set_yscale("symlog", linthresh=0.1)
    else:
        ax.set_yscale("linear")
    htags.append((ax, 1.0, "Unity", ":"))
    ax.set_ylabel("Luminosity ratio (outer / inner)")
    if has_lum:
        _legend(ax, loc="best", fontsize=7, ncol=3)

    # 3: field maturity
    ax = axes[2]
    has_fg = any(_fin(v) and v > 0 for v in d["fgd"])
    ax.semilogy(xs, pos(d["fgd"]), "o-", ms=4, c="tab:red",
                label="Field drift (max rel. change)")
    if has_fg:
        _full_frame(ax, d["mode"], d["fgd"])
    _empty_note(ax, has_fg)
    htags.append((ax, 5e-3, "Maturity threshold 5e-3", ":"))
    ax.set_ylabel("FG field drift")
    fgb_pct = [v * 100 if _fin(v) else float("nan") for v in d["fgb"]]
    if any(_fin(v) and v > 1e-4 for v in d["fgb"]):
        ax2 = ax.twinx()
        ax2.plot(xs, fgb_pct, "v--", ms=3, c="tab:purple", alpha=0.6,
                 label="Fraction of bad cells")
        ax2.set_ylabel("Fraction of bad field cells (%)", color="tab:purple")
        ax2.tick_params(axis="y", colors="tab:purple")
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        _legend(ax, h1 + h2, l1 + l2, loc="upper center", fontsize=8)
    elif has_fg:
        _legend(ax, loc="upper center", fontsize=8)

    # 4: SE residual + cumulative rejects
    ax = axes[3]
    has_se = any(_fin(v) and v > 0 for v in d["se"] + d["se2"])
    ax.semilogy(xs, pos(d["se"]), "s--", ms=3, c="0.45", alpha=0.7,
                label="SE residual (se)")
    ax.semilogy(xs, pos(d["se2"]), "o-", ms=4, c="tab:brown",
                label="SE residual excl. worst (res2)")
    if has_se:
        _full_frame(ax, d["mode"], d["se2"])
    _empty_note(ax, has_se)
    if any(_fin(v) and v > SATURATED_RES for v in d["se"]):
        htags.append((ax, SATURATED_RES, "Saturated channel (1e30)", ":"))
    ax.set_ylabel("SE residual")
    cum, s = [], 0
    for v in d["nrej"]:
        s += int(v) if _fin(v) else 0
        cum.append(s)
    if cum and max(cum) > 0:
        rjax = ax.twinx()
        rjax.step(xs, cum, where="post", c="grey", alpha=0.7,
                  label="Cumulative BA rejects")
        rjax.set_ylim(-0.5, max(cum[-1] + 1, 2))
        rjax.set_ylabel("Cumulative BA rejects", color="grey")
        rjax.tick_params(axis="y", colors="grey")
        rjax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = rjax.get_legend_handles_labels()
        _legend(ax, h1 + h2, l1 + l2, loc="upper center", fontsize=8)
    else:
        _legend(ax, loc="upper center", fontsize=8)

    # 5: wall time — cumulative on the left, per-iteration on a twin axis,
    # each axis in its own most-readable unit
    if npan == 5:
        from itertools import accumulate
        tx = [x for x, _ in tt]
        per = [v for _, v in tt]
        tcum = list(accumulate(per))
        u_cum, s_cum = _time_unit(tcum[-1])
        u_per, s_per = _time_unit(max(per))
        ax = axes[4]
        ax.plot(tx, [v / s_cum for v in tcum], "o-", ms=4, c="tab:cyan",
                label="Cumulative wall time")
        cum_by_x = dict(zip(tx, [v / s_cum for v in tcum]))
        _full_frame(ax, d["mode"],
                    [cum_by_x.get(x, float("nan")) for x in xs],
                    logpos=False)
        ax.set_ylabel(f"Cumulative ({u_cum})")
        ax2 = ax.twinx()
        ax2.plot(tx, [v / s_per for v in per], "v--", ms=3, c="grey",
                 alpha=0.55, label="Per-iteration")
        ax2.set_ylabel(f"Per iteration ({u_per})", color="grey")
        ax2.tick_params(axis="y", colors="grey")
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        _legend(ax, h1 + h2, l1 + l2, loc="upper center", fontsize=8)

    axes[-1].set_xlabel("Iteration"
                        + (" (counter reset — segments run left to right)"
                           if d["rst"] else ""))

    # tick labels show the actual iteration number at each position;
    # repeated numbers mark restarted segments
    step = max(1, len(xs) // 10)
    tk = list(range(0, len(xs), step))
    axes[-1].set_xticks(tk)
    axes[-1].set_xticklabels([str(it[i]) for i in tk], fontsize=8)

    # Resolve final limits, then draw threshold tags: label at the left edge
    # centred on the level, line running from the label's right end to the
    # axes edge.
    fig.canvas.draw()
    for ax, y, label, ls in htags:
        t = ax.text(0.01, y, f" {label} ",
                    transform=ax.get_yaxis_transform(),
                    va="center", ha="left", fontsize=8,
                    zorder=10, bbox=TBBOX)
        x1 = ax.transAxes.inverted().transform(
            t.get_window_extent().max)[0]
        ax.axhline(y, xmin=min(x1 + 0.01, 0.98), xmax=1.0,
                   ls=ls, c="k", lw=1, zorder=5)
    for r in d["rst"]:
        axes[0].text(r / float(len(xs)), 0.985, "restart",
                     transform=axes[0].transAxes, fontsize=7, va="top",
                     ha="center", color="tab:red", zorder=10, bbox=TBBOX)

    last = len(it) - 1
    run = "/".join(hist_path.resolve().parent.parts[-3:])
    axes[0].set_title(
        f"{run}\n"
        f"Iteration {it[-1]}  {d['mode'][last]} → {d['nxt'][last]}  "
        f"BA={d['ba'][last]}  Frozen={d['froz'][last]}  "
        f"Held={d['held'][last]}   ({time.strftime('%Y-%m-%d %H:%M')})",
        fontsize=10)
    fig.savefig(out_path, dpi=110, bbox_inches="tight")
    plt.close(fig)


def _plot_svg(d, ms, hist_path, out_path):
    """Dependency-free fallback: stacked-panel SVG."""
    from xml.sax.saxutils import escape as esc
    it = d["it"]           # actual iteration numbers (tick labels)
    xs = d["x"]            # monotonic row positions (the x-axis)
    chg_series = []
    if any(_fin(v) for v in d["pinc"]):
        chg_series.append(("Max % increase", d["pinc"], "#2ca02c"))
    if any(_fin(v) for v in d["pdec_eq"]):
        chg_series.append(("Max decrease demand", d["pdec_eq"], "#d62728"))
    if any(_fin(v) for v in d["solveba"]):
        chg_series.append(("SOLVEBA_V13 returned max", d["solveba"], "#9467bd"))
    panels = [
        {"ylabel": "Max population correction", "log": True,
         "series": [("Demanded (maxraw)", d["raw"], "#1f77b4"),
                    ("Applied (maxapp)", d["app"], "#7f7f7f")]
                   + chg_series,
         "hlines": ([(1.0e7, "Report cap"), (d["lam_exit"], "BA exit")]
                    + ([(d["eps"], "EPS_TERM")]
                       if d.get("eps") is not None else []))},
        {"ylabel": "Luminosity ratio (outer / inner)",
         "log": (_lum_is_big(d)
                 and all(_fin(v) and v > 0 for v in d["lum"])),
         "series": [("Luminosity ratio", d["lum"], "#2ca02c")],
         "hlines": [(1.0, "Unity")]},
        {"ylabel": "FG field drift", "log": True,
         "series": [("Field drift", d["fgd"], "#d62728")],
         "hlines": [(5e-3, "Maturity")]},
    ]
    if any(_fin(v) and v > 1e-4 for v in d["fgb"]):
        panels.append({"ylabel": "Bad field cells (%)", "log": False,
                       "series": [("Bad cells",
                                   [v * 100 if _fin(v) else float("nan")
                                    for v in d["fgb"]], "#8c2d8c")],
                       "hlines": []})
    panels.append({"ylabel": "SE residual", "log": True,
                   "series": [("se", d["se"], "#7f7f7f"),
                              ("res2 (excl. worst)", d["se2"], "#8c564b")],
                   "hlines": ([(SATURATED_RES, "Saturated")]
                              if any(_fin(v) and v > SATURATED_RES
                                     for v in d["se"]) else [])})
    if d.get("tt"):
        from itertools import accumulate
        tmap = dict(d["tt"])
        cum = list(accumulate(tmap.get(x, 0.0) for x in xs))
        unit, scale = _time_unit(cum[-1])
        panels.append({"ylabel": f"Cumulative time ({unit})", "log": False,
                       "series": [("Cumulative time",
                                   [v / scale for v in cum], "#17becf")],
                       "hlines": []})

    W, PH, GAP, TOP, BOT, ML, MR = 960, 140, 12, 46, 30, 72, 14
    n = len(panels)
    heights = [2 * PH] + [PH] * (n - 1)          # top panel double height
    H = TOP + sum(heights) + (n - 1) * GAP + BOT
    xmin, xmax = xs[0], xs[-1]
    if xmax == xmin:
        xmax += 1
    X = lambda x: ML + (x - xmin) / (xmax - xmin) * (W - ML - MR)

    run = "/".join(hist_path.resolve().parent.parts[-3:])
    out = ['<svg xmlns="http://www.w3.org/2000/svg" '
           f'width="{W}" height="{H}" font-family="monospace" font-size="11">',
           f'<rect width="{W}" height="{H}" fill="white"/>',
           f'<text x="{W / 2}" y="16" text-anchor="middle" font-size="12">'
           f'{esc(run)}</text>',
           f'<text x="{W / 2}" y="32" text-anchor="middle" font-size="11">'
           f'Iteration {it[-1]} {esc(d["mode"][-1])} → {esc(d["nxt"][-1])} '
           f'BA={esc(d["ba"][-1])}  ({time.strftime("%Y-%m-%d %H:%M")})</text>']

    for pi, p in enumerate(panels):
        ph = heights[pi]
        y0 = TOP + sum(heights[:pi]) + pi * GAP
        ys = [v for _, s, _ in p["series"] for v in s
              if _fin(v) and (v > 0 or not p["log"])]
        if not ys:
            # a missing panel reads as "converged"; say why it is blank
            out.extend([
                f'<rect x="{ML}" y="{y0}" width="{W - ML - MR}" '
                f'height="{ph}" fill="none" stroke="#999"/>',
                f'<text x="{ML + 4}" y="{y0 + 12}" font-size="10" '
                f'fill="#555">{esc(p["ylabel"])}</text>',
                f'<text x="{(ML + W - MR) / 2:.0f}" '
                f'y="{y0 + ph / 2:.0f}" text-anchor="middle" '
                'font-size="10" fill="#999">'
                'no data on this solver path</text>'])
            continue
        lo, hi = min(ys), max(ys)
        if p["log"]:
            lo, hi = math.floor(math.log10(lo)), math.ceil(math.log10(hi))
            if hi <= lo:
                hi = lo + 1
            Y = lambda v: y0 + ph - (math.log10(v) - lo) / (hi - lo) * ph
            ticks = list(range(lo, hi + 1))
            tlab = lambda e: f"1e{e}"
            tval = lambda e: Y(10.0 ** e)
        else:
            if hi <= lo:
                hi = lo + 1
            pad = 0.08 * (hi - lo)
            lo, hi = lo - pad, hi + pad
            Y = lambda v: y0 + ph - (v - lo) / (hi - lo) * ph
            ticks = [lo + (hi - lo) * i / 4 for i in range(5)]
            tlab = lambda v: f"{v:.3g}"
            tval = Y

        out.append(f'<rect x="{ML}" y="{y0}" width="{W - ML - MR}" '
                   f'height="{ph}" fill="none" stroke="#999"/>')
        # frozen shading + restart boundaries
        for x0, x1 in _spans(d["froz"], xs):
            xa = max(X(x0 - 0.45), ML)
            xb = min(X(x1 + 0.45), W - MR)
            out.append(f'<rect x="{xa:.1f}" y="{y0}" '
                       f'width="{xb - xa:.1f}" '
                       f'height="{ph}" fill="orange" opacity="0.12"/>')
        for r in d["rst"]:
            out.append(f'<line x1="{X(r - 0.5):.1f}" x2="{X(r - 0.5):.1f}" '
                       f'y1="{y0}" y2="{y0 + ph}" stroke="#d62728" '
                       f'stroke-dasharray="2,3" opacity="0.6"/>')
            if pi == 0:
                top_rst = X(r - 0.5)
                out.append(f'<text x="{top_rst:.1f}" y="{y0 + 22}" '
                           f'text-anchor="middle" font-size="9" '
                           f'fill="#d62728">restart</text>')
        for t in ticks:
            out.append(f'<line x1="{ML}" x2="{W - MR}" y1="{tval(t):.1f}" '
                       f'y2="{tval(t):.1f}" stroke="#eee"/>'
                       f'<text x="{ML - 4}" y="{tval(t) + 4:.1f}" '
                       f'text-anchor="end" font-size="9">{tlab(t)}</text>')
        # overlay text collected separately and emitted last, with a white
        # stroke halo (paint-order), so labels sit on top of the data
        HALO = ('paint-order="stroke" stroke="white" stroke-width="3" '
                'stroke-linejoin="round"')
        top = []
        for yv, lab in p["hlines"]:
            if p["log"]:
                if yv <= 0:
                    continue
                yy = Y(yv) if lo <= math.log10(yv) <= hi else None
            else:
                yy = Y(yv) if lo <= yv <= hi else None
            if yy is not None:
                # label centred on the level at the left; line runs from
                # just past the label to the right edge
                tw = len(lab) * 9 * 0.62 + 6
                out.append(f'<line x1="{ML + 4 + tw:.1f}" x2="{W - MR}" '
                           f'y1="{yy:.1f}" y2="{yy:.1f}" stroke="#000" '
                           f'stroke-dasharray="3,3"/>')
                top.append(f'<text x="{ML + 4}" y="{yy + 3:.1f}" '
                           f'font-size="9" {HALO}>{esc(lab)}</text>')
        top.append(f'<text x="{ML + 4}" y="{y0 + 12}" font-size="10" '
                   f'fill="#555" {HALO}>{esc(p["ylabel"])}</text>')

        # series
        lx = W - MR - 4
        for name, s, color in p["series"]:
            segs, cur = [], []
            for x, v in zip(xs, s):
                ok = _fin(v) and (v > 0 or not p["log"])
                if ok:
                    cur.append(f"{X(x):.1f},{Y(v):.1f}")
                elif cur:
                    segs.append(cur)
                    cur = []
            if cur:
                segs.append(cur)
            for seg in segs:
                out.append(f'<polyline points="{" ".join(seg)}" fill="none" '
                           f'stroke="{color}" stroke-width="1.4"/>')
            top.append(f'<line x1="{lx - 40}" x2="{lx - 24}" y1="{y0 + 12}" '
                       f'y2="{y0 + 12}" stroke="{color}" stroke-width="2"/>'
                       f'<text x="{lx - 20}" y="{y0 + 15}" font-size="9" '
                       f'text-anchor="end" {HALO}>{esc(name)}</text>')
            lx -= 100
        # full iterations on every panel: transparent black-edged square
        # at the first series' value, pinned to the bottom edge when that
        # panel has no usable point at the x position
        s0 = p["series"][0][1]
        for x, m, v in zip(xs, d["mode"], s0):
            if not m.startswith("FULL"):
                continue
            py = Y(v) if _fin(v) and (v > 0 or not p["log"]) \
                else y0 + ph - 7
            out.append(f'<rect x="{X(x) - 5.5:.1f}" y="{py - 5.5:.1f}" '
                       'width="11" height="11" fill="none" '
                       'stroke="#000" stroke-width="1.5"/>')
        if pi == 0:
            # great-iteration type: per-mode dots on the demanded series
            for x, m, v in zip(xs, d["mode"], d["raw"]):
                if _fin(v) and v > 0:
                    out.append(
                        f'<circle cx="{X(x):.1f}" cy="{Y(v):.1f}" r="3" '
                        f'fill="{_SVG_MODE_COLOR.get(m, "#1f77b4")}" '
                        'fill-opacity="0.85"/>')
            # BA state: hollow glyph outlines over the mode dots
            for ba, (mk, _c, msz, _lab) in _BA_STYLE.items():
                c = _SVG_BA_COLOR[ba]
                for x, b, v in zip(xs, d["ba"], d["raw"]):
                    if b != ba or not (_fin(v) and v > 0):
                        continue
                    px, py, r = X(x), Y(v), msz / 2
                    if mk == "s":
                        out.append(f'<rect x="{px - r:.1f}" y="{py - r:.1f}" '
                                   f'width="{2 * r:.1f}" height="{2 * r:.1f}" '
                                   f'fill="none" stroke="{c}" '
                                   'stroke-width="1.5"/>')
                    elif mk in ("x", "X"):
                        out.append(f'<path d="M {px - r:.1f} {py - r:.1f} '
                                   f'L {px + r:.1f} {py + r:.1f} '
                                   f'M {px + r:.1f} {py - r:.1f} '
                                   f'L {px - r:.1f} {py + r:.1f}" '
                                   f'stroke="{c}" stroke-width="1.5" '
                                   'fill="none"/>')
                    else:
                        out.append(f'<circle cx="{px:.1f}" cy="{py:.1f}" '
                                   f'r="{r:.1f}" fill="none" stroke="{c}" '
                                   'stroke-width="1.5"/>')
            # mode legend strip below the series legend
            lxm = W - MR - 4
            for m, c in _SVG_MODE_COLOR.items():
                if m in set(d["mode"]):
                    top.append(f'<circle cx="{lxm - 8}" cy="{y0 + 24}" r="4" '
                               f'fill="{c}"/>'
                               f'<text x="{lxm - 16}" y="{y0 + 27}" '
                               f'font-size="9" text-anchor="end" {HALO}>'
                               f'{esc(m)}</text>')
                    lxm -= len(m) * 6.5 + 34
        out.extend(top)

    # x ticks on bottom panel — labels show the actual iteration number
    yb = TOP + sum(heights) + (n - 1) * GAP
    step = max(1, len(xs) // 8)
    for i in range(0, len(xs), step):
        out.append(f'<text x="{X(xs[i]):.1f}" y="{yb + 14}" '
                   f'text-anchor="middle" font-size="9">{it[i]}</text>')
    out.append(f'<text x="{(ML + W - MR) / 2:.0f}" y="{yb + 27}" '
               f'text-anchor="middle" font-size="10">Iteration</text>')
    out.append("</svg>")
    out_path.write_text("\n".join(out))


def resolve(arg):
    p = Path(arg)
    if p.is_dir():
        p = p / "ITER_HISTORY"
    return p


def main():
    ap = argparse.ArgumentParser(description="KNOP ITER_HISTORY convergence report")
    ap.add_argument("path", nargs="?", default=".",
                    help="ITER_HISTORY file or run directory containing one")
    ap.add_argument("--tail", type=int, default=10, help="window size (default 10)")
    ap.add_argument("--watch", type=float, default=0,
                    help="re-print every SEC seconds; stops when the model "
                         "converges (MAXCH < EPS_TERM), when output goes "
                         "stale, or on Ctrl-C")
    ap.add_argument("--stale", type=float, default=1800.0, metavar="SEC",
                    help="in watch mode, stop if ITER_HISTORY/OUTGEN have "
                         "not been written for SEC seconds (default 1800; "
                         "0 disables)")
    ap.add_argument("--plot", nargs="?", const="CONVERGENCE_HISTORY.PNG",
                    metavar="FILE", default=None,
                    help="write convergence-history figure next to "
                         "ITER_HISTORY (default: CONVERGENCE_HISTORY.PNG)")
    a = ap.parse_args()

    last_sig, last_change = None, time.time()
    while True:
        hist = resolve(a.path)
        if not hist.exists():
            print(f"no ITER_HISTORY at {hist}", file=sys.stderr)
            if not a.watch:
                return 1
            time.sleep(a.watch)
            continue
        sig = _run_signature(hist)
        if a.watch and sig == last_sig:
            # nothing written since the last cycle: skip the report and
            # the plot, only the stale timer still matters
            idle = time.time() - last_change
            if a.stale and idle > a.stale:
                unit, scale = _time_unit(idle)
                print(f"\n*** RUN APPEARS STOPPED: no ITER_HISTORY/OUTGEN "
                      f"writes for {idle / scale:.0f} {unit} — "
                      f"stopping watch ***")
                break
            time.sleep(a.watch)
            continue
        last_sig, last_change = sig, time.time()
        lines = _outgen_lines(hist)
        print(f"\n===== check_convergence  {time.strftime('%H:%M:%S')} =====")
        report(hist, a.tail, lines)
        scan_outgen(hist, lines)
        if a.plot:
            out = Path(a.plot)
            if not out.is_absolute():
                out = hist.parent / out
            try:
                written = make_plot(hist, out, lines)
                if written:
                    print(f"\nplot    : {written}")
            except Exception as e:
                print(f"\nplot failed: {e}", file=sys.stderr)
        if not a.watch:
            break
        conv, maxch, eps = converged_status(hist, lines)
        if conv:
            print(f"\n*** MODEL CONVERGED: max demanded correction "
                  f"{maxch:.3e} < EPS_TERM={eps:g}% — stopping watch ***")
            break
        time.sleep(a.watch)
    return 0


if __name__ == "__main__":
    sys.exit(main())
