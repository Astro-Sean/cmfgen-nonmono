#!/usr/bin/env python3
"""Check running and recently completed CMFGEN jobs and report convergence progress.

Uses the robust multi-iteration convergence diagnostics from cmfgen_convergence.py
which distinguishes:
  CONVERGING, DAMPED_OSCILLATION, PERSISTENT_OSCILLATION,
  STAGNATING, DIVERGING, CONVERGED, CRASHED

CMFGEN convergence is NOT necessarily monotonic iteration-by-iteration.
A model with occasional increases in corrections can still be converging.
"""

import subprocess
import re
import os
import glob
import sys
from pathlib import Path

# Import the robust convergence diagnostics
sys.path.insert(0, str(Path(__file__).parent))
from cmfgen_convergence import (
    parse_outgen as robust_parse_outgen,
    assess_convergence as robust_assess_convergence,
    format_assessment,
    ConvergenceHistory,
)

SEARCH_ROOTS = [
    "/ptmp/sbrennan/SN2024pba",
    "/ptmp/sbrennan/jwst_cycle6",
]

# Regex patterns for OUTGEN parsing
ITER_RE = re.compile(rb'Current great iteration count is\s+(\d+)')
LUM_RE = re.compile(rb'Luminosity of star \(d=1,ND\)\(iteration\s+(\d+)\) is:\s+([\d.eE+\-]+)')
MAX_INC_RE = re.compile(rb'Maximum % increase at depth\s+(\d+)\s+is\s+([\d.eE+\-]+)')
MAX_DEC_RE = re.compile(rb'Maximum % decrease at depth\s+(\d+)\s+is\s+([\d.eE+\-]+)')
NT_RE = re.compile(rb'NT=\s+(\d+)')
# CORRECTION_SUM: last line has convergence counts at various thresholds
CORR_DEPTH_RE = re.compile(rb'^\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)', re.MULTILINE)
# Check for crash/error patterns
CRASH_RE = re.compile(rb'forrtl: severe| segmentation fault| SIGSEGV|NaN|NaN detected')
COMPLETE_RE = re.compile(rb'CMFGEN completed|Normal termination|END OF PROGRAM')


def get_cmfgen_jobs():
    """Get all CMFGEN jobs (running + recently completed) from sacct/squeue."""
    jobs = []

    # Running/pending from squeue
    try:
        result = subprocess.run(
            ["squeue", "-u", os.environ.get("USER", "sbrennan"),
             "--format=%.10i %.40j %.8T %.10M"],
            capture_output=True, text=True, timeout=10
        )
        for line in result.stdout.splitlines()[1:]:
            parts = line.split(None, 3)
            if len(parts) < 4:
                continue
            job_id, job_name, state, elapsed = parts
            if "cmfgen" in job_name.lower():
                jobs.append({
                    "id": job_id, "name": job_name, "state": state,
                    "elapsed": elapsed, "source": "squeue"
                })
    except Exception:
        pass

    return jobs


def parse_job_name(name):
    """Parse CMFGEN job name into model, epoch, profile.
    e.g. X_2_5.0days_cmfgen_INIT -> model=X_2, epoch=5.0days, profile=INIT
         SN_A_30.0days_cmfgen_STABLE -> model=SN_A, epoch=30.0days, profile=STABLE
    """
    # Remove _cmfgen_ and split
    m = re.match(r'(.+?)_(\d+\.?\d*days)_cmfgen_(\w+)', name)
    if m:
        return m.group(1), m.group(2), m.group(3)
    return name, "?", "?"


def find_cmfgen_dir(model, epoch):
    """Find the CMFGEN directory for a model/epoch."""
    for root in SEARCH_ROOTS:
        # Standard pattern: models/MODEL/cmfgen/EPOCH/
        pattern = f"{root}/**/cmfgen/{epoch}/"
        matches = glob.glob(pattern, recursive=True)
        # Filter to matching model
        for match in matches:
            if model in match:
                return match
        if matches:
            return matches[0]

        # Also check MODEL/cmfgen/EPOCH/ (for SN_A, SN_B)
        pattern2 = f"{root}/{model}/cmfgen/{epoch}/"
        if os.path.isdir(pattern2):
            return pattern2

    return None


def parse_outgen(cmfgen_dir):
    """Parse OUTGEN file for convergence progress.

    Uses the robust multi-iteration parser from cmfgen_convergence.py
    in addition to the legacy single-value extraction.
    """
    outgen = Path(cmfgen_dir) / "OUTGEN"
    if not outgen.exists():
        return None

    # Use the robust parser for full iteration history
    history = robust_parse_outgen(outgen)

    try:
        with open(outgen, 'rb') as f:
            f.seek(0, 2)
            size = f.tell()
            f.seek(max(0, size - 100000))
            tail = f.read()
    except Exception:
        return None

    # Find last iteration count
    iter_matches = ITER_RE.findall(tail)
    last_iter = int(iter_matches[-1]) if iter_matches else None

    # Find last luminosity
    lum_matches = LUM_RE.findall(tail)
    last_lum_iter = int(lum_matches[-1][0]) if lum_matches else None
    last_lum = float(lum_matches[-1][1]) if lum_matches else None

    # Find first luminosity (for convergence tracking)
    first_lum = float(lum_matches[0][1]) if lum_matches else None

    # Find last max % increase/decrease
    inc_matches = MAX_INC_RE.findall(tail)
    last_inc_depth = int(inc_matches[-1][0]) if inc_matches else None
    last_inc_val = float(inc_matches[-1][1]) if inc_matches else None

    dec_matches = MAX_DEC_RE.findall(tail)
    last_dec_depth = int(dec_matches[-1][0]) if dec_matches else None
    last_dec_val = float(dec_matches[-1][1]) if dec_matches else None

    # Check for crash
    crashed = bool(CRASH_RE.search(tail))
    completed = bool(COMPLETE_RE.search(tail))

    # Track luminosity convergence over last few iterations
    lum_history = []
    for m in lum_matches[-5:]:
        lum_history.append(float(m[1]))

    return {
        "iteration": last_iter,
        "lum_iter": last_lum_iter,
        "lum": last_lum,
        "first_lum": first_lum,
        "lum_history": lum_history,
        "max_inc_depth": last_inc_depth,
        "max_inc": last_inc_val,
        "max_dec_depth": last_dec_depth,
        "max_dec": last_dec_val,
        "crashed": crashed,
        "completed": completed,
        "history": history,  # Full multi-iteration history for robust assessment
    }


def parse_elapsed_to_minutes(elapsed_str):
    """Parse SLURM elapsed time to minutes."""
    s = elapsed_str.strip()
    if "-" in s:
        days, s = s.split("-", 1)
        days = int(days)
    else:
        days = 0
    parts = s.split(":")
    if len(parts) == 3:
        h, m, sec = int(parts[0]), int(parts[1]), int(parts[2])
    elif len(parts) == 2:
        h, m, sec = 0, int(parts[0]), int(parts[1])
    else:
        h, m, sec = 0, 0, int(parts[0])
    return days * 1440 + h * 60 + m + sec / 60.0


def format_lum(lum):
    """Format luminosity in scientific notation."""
    if lum is None:
        return "?"
    if lum >= 1e40:
        return f"{lum:.2e}"
    return f"{lum:.3e}"


def assess_convergence(info):
    """Assess CMFGEN convergence status using robust multi-iteration analysis.

    This replaces the old single-iteration threshold approach with a
    rolling-window trend analysis that correctly handles:
    - convergence with occasional increases
    - damped oscillation
    - persistent oscillation
    - stagnation
    - divergence

    CMFGEN convergence is NOT necessarily monotonic iteration-by-iteration.
    """
    if info is None:
        return "no OUTGEN"

    if info["crashed"]:
        return "CRASHED"

    if info["completed"]:
        return "COMPLETED"

    # Use the robust assessment if we have enough iterations
    if info.get("history") is not None and info["history"].n_iterations >= 3:
        assessment = robust_assess_convergence(info["history"], window=10)
        return f"{assessment.classification} (max%={info['max_inc']:.0e} at d{info['max_inc_depth']}, slope={assessment.slope_log10:.3f})"

    # Fallback for early iterations
    if info["max_inc"] is not None:
        if info["max_inc"] >= 1e4:
            if len(info["lum_history"]) >= 3:
                lums = info["lum_history"]
                rel_change = abs(lums[-1] - lums[-2]) / max(abs(lums[-2]), 1e-30)
                if rel_change < 0.01:
                    return f"converging (L stable, max%={info['max_inc']:.0e} clipped)"
                else:
                    return f"iterating (L changing {rel_change*100:.1f}%, max%={info['max_inc']:.0e})"
            return f"iterating (max%={info['max_inc']:.0e} at d{info['max_inc_depth']})"
        elif info["max_inc"] < 10:
            return f"near convergence (max%={info['max_inc']:.1f})"
        elif info["max_inc"] < 1000:
            return f"converging (max%={info['max_inc']:.0f})"
        else:
            return f"slow (max%={info['max_inc']:.0e})"

    return "starting"


def main():
    import argparse
    parser = argparse.ArgumentParser(
        description="Check running and recently completed CMFGEN jobs."
    )
    parser.add_argument('--detailed', action='store_true',
                        help='Show full multi-iteration convergence report for each job.')
    parser.add_argument('--window', type=int, default=10,
                        help='Rolling window size for convergence assessment (default: 10).')
    args = parser.parse_args()

    jobs = get_cmfgen_jobs()

    if not jobs:
        print("No CMFGEN jobs found (running or completed in last 24h).")
        return

    # Sort: running first, then by name
    jobs.sort(key=lambda j: (0 if j["state"] == "RUNNING" else 1, j["name"]))

    print(f"{'Model':<16} {'Epoch':<10} {'Profile':<10} {'JobID':<10} {'State':<10} "
          f"{'Runtime':<10} {'Iter':<6} {'Luminosity':<14} {'Max%inc':<10} {'Status'}")
    print("-" * 130)
    # Profile = CMFGEN phase: INIT (fixed T/Ne) -> STABLE (relax) -> FINAL (spectrum)
    # Iter = current great iteration count (Lambda iterations)
    # Luminosity = L(d=1,ND) in erg/s at current iteration
    # Max%inc = max population change this iteration (goal: <10% = converged)

    for job in jobs:
        model, epoch, profile = parse_job_name(job["name"])
        job_id = job["id"]
        state = job["state"]
        elapsed = job["elapsed"]

        cmfgen_dir = find_cmfgen_dir(model, epoch)
        info = parse_outgen(cmfgen_dir) if cmfgen_dir else None

        if info:
            it = info["iteration"] or "?"
            lum = format_lum(info["lum"])
            if info["max_inc"] is not None:
                if info["max_inc"] >= 1e4:
                    max_str = f"{info['max_inc']:.0e}"
                else:
                    max_str = f"{info['max_inc']:.1f}"
            else:
                max_str = "?"
            status = assess_convergence(info)
        else:
            it = "?"
            lum = "?"
            max_str = "?"
            if state == "PENDING":
                status = "pending"
            elif state == "COMPLETED":
                status = "completed (no OUTGEN)"
            elif state == "FAILED":
                ec = job.get("exitcode", "?")
                status = f"FAILED (exit={ec})"
            else:
                status = "no OUTGEN"

        print(f"{model:<16} {epoch:<10} {profile:<10} {job_id:<10} {state:<10} "
              f"{elapsed:<10} {it:<6} {lum:<14} {max_str:<10} {status}")

        # Detailed convergence report
        if args.detailed and info and info.get("history") is not None:
            assessment = robust_assess_convergence(info["history"], window=args.window)
            print()
            print(format_assessment(assessment))
            print()

    print()
    print("Legend: Profile=INIT/STABLE/FINAL phase, Iter=Lambda iteration #, "
          "Luminosity=erg/s, Max%inc=max pop change (goal<10%)")
    print("Use --detailed for full multi-iteration convergence report.")


if __name__ == "__main__":
    main()
