#!/usr/bin/env python3
"""Check running STELLA jobs and report their simulation time progress."""

import subprocess
import re
import os
import glob
import sys
from pathlib import Path

# Search roots for STELLA output files
SEARCH_ROOTS = [
    "/ptmp/sbrennan/SN2024pba",
    "/ptmp/sbrennan/jwst_cycle6",
]

# Regex for parsing td= values from STELLA output
TD_RE = re.compile(rb'td=\s*([\d.eE+\-]+)')
ISTIM_RE = re.compile(rb'istim=\s*(\d+)')
NZMOD_RE = re.compile(rb'Nzmod=\s*(\d+)')
GY_RE = re.compile(rb'GY:.*?H=\s*([\d.eE+\-]+)\s+Hu=\s*([\d.eE+\-]+)')
KFLAG_RE = re.compile(rb'GY:.*?kflag,\s*H=.*?kflag\s*=\s*(-?\d+)')
# Also match the inline format: "6636  1444  1 Hu= 3.379E-04 ..."
INLINE_RE = re.compile(rb'\d+\s+\d+\s+\d+\s+Hu=\s*([\d.eE+\-]+)\s+Ro=\s*([\d.eE+\-]+)\s+td=\s*([\d.eE+\-]+)')
DAT_LAST_RE = re.compile(r'last=\s*([\d.]+)\s*days')
ERROR_HMIN_RE = re.compile(rb'error test fails for H=HMIN')


def get_running_stella_jobs():
    """Get list of running STELLA jobs from squeue."""
    try:
        result = subprocess.run(
            ["squeue", "-u", os.environ.get("USER", "sbrennan"),
             "--format=%.10i %.40j %.8T %.10M"],
            capture_output=True, text=True, timeout=10
        )
    except Exception:
        return []

    jobs = []
    for line in result.stdout.splitlines()[1:]:
        parts = line.split(None, 3)
        if len(parts) < 4:
            continue
        job_id, job_name, state, elapsed = parts
        if "stella" in job_name.lower() and state == "RUNNING":
            jobs.append({"id": job_id, "name": job_name, "elapsed": elapsed})
    return jobs


def find_out_file(job_id, job_name):
    """Find the STELLA .out file for a given job."""
    # Job name like "X_1_stella" -> model dir could be models/X_1/stella/
    # or "SN_A_stella" -> SN_A/stella/
    # or "jwst_highNi_stella" -> models/jwst_highNi/stella/
    patterns = []
    for root in SEARCH_ROOTS:
        patterns.extend([
            f"{root}/**/stella/stella_*_{job_id}.out",
            f"{root}/**/stella/*{job_id}*.out",
        ])

    for pattern in patterns:
        matches = glob.glob(pattern, recursive=True)
        if matches:
            return matches[0]

    # Fallback: search by job name (without _stella suffix)
    model_name = job_name.replace("_stella", "")
    for root in SEARCH_ROOTS:
        pattern = f"{root}/**/stella/stella_{model_name}_{job_id}.out"
        matches = glob.glob(pattern, recursive=True)
        if matches:
            return matches[0]

    return None


def get_target_days(out_file):
    """Extract target days from the .out file (last=XX.XX days in the Written line)."""
    try:
        with open(out_file, 'rb') as f:
            f.seek(0, 2)
            size = f.tell()
            f.seek(max(0, size - 100000))
            tail = f.read()
        m = DAT_LAST_RE.search(tail.decode('utf-8', errors='replace'))
        if m:
            return float(m.group(1))
    except Exception:
        pass

    # Fallback: also check beginning of file
    try:
        with open(out_file, 'rb') as f:
            head = f.read(50000)
        m = DAT_LAST_RE.search(head.decode('utf-8', errors='replace'))
        if m:
            return float(m.group(1))
    except Exception:
        pass

    return None


def parse_stella_progress(out_file):
    """Parse the last td=, istim=, Nzmod=, and Hu= values from STELLA output."""
    try:
        with open(out_file, 'rb') as f:
            # Read last 50KB for efficiency
            f.seek(0, 2)
            size = f.tell()
            f.seek(max(0, size - 50000))
            tail = f.read()
    except Exception:
        return None

    # Find all td= values
    td_matches = INLINE_RE.findall(tail)
    last_td = float(td_matches[-1][2]) if td_matches else None

    # Also try standalone td= pattern
    if last_td is None:
        td_matches = TD_RE.findall(tail)
        last_td = float(td_matches[-1]) if td_matches else None

    # Find last istim=
    istim_matches = ISTIM_RE.findall(tail)
    last_istim = int(istim_matches[-1]) if istim_matches else None

    # Find last Nzmod=
    nzmod_matches = NZMOD_RE.findall(tail)
    last_nzmod = int(nzmod_matches[-1]) if nzmod_matches else None

    # Find last Hu= from GY lines
    gy_matches = GY_RE.findall(tail)
    last_hu = float(gy_matches[-1][1]) if gy_matches else None

    # Check for HMIN error (stuck)
    hmin_error = bool(ERROR_HMIN_RE.search(tail))

    # Count total td= lines for rate estimation
    total_td_lines = len(td_matches) + len(INLINE_RE.findall(tail))

    return {
        "td": last_td,
        "istim": last_istim,
        "nzmod": last_nzmod,
        "hu": last_hu,
        "hmin_error": hmin_error,
    }


def parse_elapsed_to_minutes(elapsed_str):
    """Parse SLURM elapsed time (H:MM:SS or D-H:MM:SS) to minutes."""
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


def main():
    jobs = get_running_stella_jobs()

    if not jobs:
        print("No running STELLA jobs found.")
        return

    print(f"{'Model':<25} {'JobID':<10} {'Runtime':<10} {'td (days)':<12} "
          f"{'Target':<10} {'Progress':<10} {'istim':<7} {'Nzmod':<10} {'Status'}")
    print("-" * 120)
    # td = simulation time in days (current physical time STELLA has reached)
    # istim = output epoch index (how many snapshots written; 0 = pre-SBO)
    # Nzmod = total timesteps taken so far

    for job in jobs:
        name = job["name"].replace("_stella", "")
        job_id = job["id"]
        elapsed = job["elapsed"]
        elapsed_min = parse_elapsed_to_minutes(elapsed)

        out_file = find_out_file(job_id, job["name"])

        if not out_file:
            print(f"{name:<25} {job_id:<10} {elapsed:<10} {'?':<12} "
                  f"{'?':<10} {'?':<10} {'?':<7} {'?':<10} No .out file found")
            continue

        progress = parse_stella_progress(out_file)
        target = get_target_days(out_file)

        if progress is None:
            print(f"{name:<25} {job_id:<10} {elapsed:<10} {'?':<12} "
                  f"{target or '?':<10} {'?':<10} {'?':<7} {'?':<10} Cannot parse output")
            continue

        td = progress["td"]
        istim = progress["istim"]
        nzmod = progress["nzmod"]
        hu = progress["hu"]

        # Calculate progress percentage
        if target and td is not None and target > 0:
            pct = f"{td / target * 100:.1f}%"
        else:
            pct = "?"

        # Determine status
        status_parts = []

        if progress["hmin_error"]:
            status_parts.append("STUCK (HMIN error)")

        if td is not None:
            if td < 0.01 and elapsed_min > 30:
                status_parts.append("STUCK (td~0 after >30min)")
            elif td < 0.1 and elapsed_min > 60:
                status_parts.append("SLOW (td<0.1 after >1hr)")
            elif target and td / target < 0.05 and elapsed_min > 120:
                status_parts.append("SLOW (<5% after >2hr)")

        if hu is not None and hu < 1e-12:
            status_parts.append("tiny timesteps")

        if not status_parts:
            if td is not None and target and td >= target * 0.9:
                status_parts.append("near completion")
            elif td is not None and td > 0:
                status_parts.append("OK")
            else:
                status_parts.append("starting")

        status = ", ".join(status_parts)

        td_str = f"{td:.3f}" if td is not None else "?"
        istim_str = str(istim) if istim is not None else "?"
        nzmod_str = f"{nzmod:,}" if nzmod is not None else "?"

        print(f"{name:<25} {job_id:<10} {elapsed:<10} {td_str:<12} "
              f"{target or '?':<10} {pct:<10} {istim_str:<7} {nzmod_str:<10} {status}")

    print()
    print("Legend: td=sim time (days), istim=snapshot # (0=pre-SBO), "
          "Nzmod=total timesteps, Hu=timestep size (s)")


if __name__ == "__main__":
    main()
