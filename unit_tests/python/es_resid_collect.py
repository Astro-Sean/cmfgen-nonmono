#!/usr/bin/env python3
"""P8.2 ES_RESID_DIAG collector: residual histories + final audit.

Reads a sweep directory tree of the form

    <root>/obs_p82_<sob|knop>_n<N>[_diagoff]/ES_RESID_DIAG
    <root>/obs_p82_<sob|knop>_n<N>[_diagoff]/batobs.log   (runtime)

and writes

    <out>/p82_residual_history.csv   one row per completed ES
                                     generation per run (includes
                                     partial histories from crashed
                                     runs)
    <out>/p82_final_audit.csv        configured/completed sweeps,
                                     refresh generation, render-input
                                     generation, render match
    <out>/p82_runtimes.csv           wall time per run (diag-on and
                                     diag-off controls)

Spectral metrics are produced separately by
unit_tests/python/es_sweep_compare.py (pairwise obs_fin_10
comparison); this collector only handles the Fortran diagnostic
files so partial runs still contribute residual data.

Usage:
    es_resid_collect.py --root <model dir> --out <outdir>
"""

import argparse
import csv
import re
from datetime import datetime
from pathlib import Path

TIME_FMT = "%a %b %d %H:%M:%S %Z %Y"


def _grab(ln, key, cast=float):
    m = re.search(key + r"=\s*([^\s]+)", ln)
    return cast(m.group(1)) if m else None


def parse_diag(path):
    """Return (global rows, extra FINAL/RENDER fields)."""
    rows, extra = [], {}
    depth_sums = {}
    if not Path(path).exists():
        return rows, extra
    for ln in Path(path).read_text().splitlines():
        if ln.startswith("DEPTH SWEEP="):
            sw = int(re.search(r"SWEEP=(\d+)", ln).group(1))
            nm = re.search(r"NMASK=(\d+)", ln)
            depth_sums[sw] = depth_sums.get(sw, 0) + (
                int(nm.group(1)) if nm else 0
            )
            continue
        if ln.startswith("GLOBAL ES_SWEEP="):
            sw = _grab(ln, "ES_SWEEP", int)
            rows.append(dict(
                sweep=sw,
                gen=_grab(ln, "REFRESH_GENERATION", int),
                r=_grab(ln, "GLOBAL_MAX_R", float),
                k=_grab(ln, "DEPTH", int),
                dmax=_grab(ln, "GLOBAL_MAX_ABS", float),
                rrel=_grab(ln, "GLOBAL_MAX_REL", float),
                nact=_grab(ln, "TOTAL_ACTIVE", int),
                nmask=depth_sums.get(sw, 0),
                nnonfin=_grab(ln, "TOTAL_INVALID", int),
                nneg=_grab(ln, "TOTAL_NEG", int),
                first=("FIRST_GENERATION" in ln),
            ))
        elif ln.startswith("FINAL "):
            for key in ("CONFIGURED_ES_SWEEPS", "COMPLETED_ES_SWEEPS",
                        "FINAL_REFRESH_GENERATION",
                        "FINAL_RENDER_INPUT_GENERATION"):
                m = re.search(key + r"=\s*([^\s]+)", ln)
                if m:
                    extra[key] = m.group(1)
            extra["RENDER_MATCH"] = "RENDER_MATCH=T" in ln
        elif ln.startswith("RENDER "):
            m = re.search(r"FIELD_MATCH=([TF])", ln)
            if m:
                extra["FIELD_MATCH"] = m.group(1) == "T"
    return rows, extra


def parse_runtime(path):
    if not Path(path).exists():
        return None
    lines = Path(path).read_text().splitlines()
    started = finished = None
    for i, ln in enumerate(lines):
        if "Program started on:" in ln:
            started = lines[i + 1].strip()
        if "Program finished on:" in ln:
            finished = lines[i + 1].strip()
    if started is None or finished is None:
        return None
    try:
        return (datetime.strptime(finished, TIME_FMT)
                - datetime.strptime(started, TIME_FMT)).total_seconds()
    except ValueError:
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True,
                    help="dir containing obs_p82_* run dirs")
    ap.add_argument("--out", required=True, help="output dir for CSV/PNG")
    args = ap.parse_args()
    root = Path(args.root)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    runs = sorted(root.glob("obs_p82_*"))
    if not runs:
        raise SystemExit(f"no obs_p82_* dirs under {root}")

    diags, timing = {}, {}
    for d in runs:
        m = re.match(r"obs_p82_(sob|knop)_n(\d+)(_diagoff)?", d.name)
        if not m:
            continue
        key = (m.group(1), int(m.group(2)), bool(m.group(3)))
        diags[key] = parse_diag(d / "ES_RESID_DIAG")
        secs = parse_runtime(d / "batobs.log")
        if secs is not None:
            timing[key] = secs

    with open(out / "p82_residual_history.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([
            "pathway", "num_es", "sweep", "refresh_generation",
            "global_max_r", "depth_at_max", "global_max_abs_d",
            "global_max_active_rel", "total_active", "total_masked",
            "total_nonfinite", "total_negative", "first_generation",
        ])
        for (pw, n, off), (rows, _x) in sorted(diags.items()):
            for r in rows:
                w.writerow([pw, n, r["sweep"], r["gen"], r["r"], r["k"],
                            r["dmax"], r["rrel"], r["nact"], r["nmask"],
                            r["nnonfin"], r["nneg"], r["first"]])

    with open(out / "p82_final_audit.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["pathway", "num_es", "diag", "configured",
                    "completed", "final_refresh_gen", "render_input_gen",
                    "render_field_match"])
        for (pw, n, off), (_rows, extra) in sorted(diags.items()):
            w.writerow([pw, n, not off,
                        extra.get("CONFIGURED_ES_SWEEPS"),
                        extra.get("COMPLETED_ES_SWEEPS"),
                        extra.get("FINAL_REFRESH_GENERATION"),
                        extra.get("FINAL_RENDER_INPUT_GENERATION"),
                        extra.get("RENDER_MATCH")])

    with open(out / "p82_runtimes.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["pathway", "num_es", "diag", "runtime_s"])
        for key, secs in sorted(timing.items()):
            w.writerow([key[0], key[1], not key[2], secs])

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib unavailable - CSV only")
        return

    fig, ax = plt.subplots(figsize=(6, 4))
    for (pw, n, off), (rows, _x) in sorted(diags.items()):
        if off or not rows:
            continue
        gs = [r["gen"] for r in rows if r["r"] is not None]
        rs = [r["r"] for r in rows if r["r"] is not None]
        if gs:
            ax.plot(gs, rs, "o-", label=f"{pw} NUM_ES={n}")
    ax.set_yscale("log")
    ax.set_xlabel("RJ_ES refresh generation")
    ax.set_ylabel("R_GLOBAL")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(out / "p82_residual_history.png", dpi=150)
    plt.close(fig)
    print("plots written")


if __name__ == "__main__":
    main()
