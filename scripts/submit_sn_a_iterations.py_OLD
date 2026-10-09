#!/usr/bin/env python3
"""Submit SN_A iteration models one at a time via SLURM dependencies.

Submits the first model, captures its job ID, then submits each subsequent
model with --stella-dependency=afterok:PREV_JOBID so they run sequentially.

Usage:
    python scripts/submit_sn_a_iterations.py
    python scripts/submit_sn_a_iterations.py --dry-run
"""
import subprocess
import sys
import os

BASE = "/viper/ptmp1/sbrennan/SN2024pba"

# Models in order — each depends on the previous completing
MODELS = [
    "baseline",
    "Menv010",
    "R300",
    "Menv010_R300",
    "T15k",
    "Menv010_R300_T15k",
    "Menv010_R300_T15k_k40",
    "Menv010_R300_T15k_k40_Ni07",
]

def main():
    dry_run = "--dry-run" in sys.argv

    print(f"Submitting {len(MODELS)} SN_A iteration models sequentially")
    print(f"(each model depends on the previous completing successfully)\n")

    prev_job_id = None
    submitted = []

    for i, name in enumerate(MODELS):
        config = f"{BASE}/SN_A/{name}/config.yaml"
        if not os.path.exists(config):
            print(f"ERROR: config not found: {config}")
            sys.exit(1)

        cmd = ["python", "scripts/run_model_grid.py", config]
        if prev_job_id:
            cmd.extend(["--stella-dependency", f"afterok:{prev_job_id}"])

        dep_str = f" (depends on job {prev_job_id})" if prev_job_id else ""
        print(f"[{i+1}/{len(MODELS)}] Submitting {name}{dep_str}...")

        if dry_run:
            print(f"  DRY RUN: would run: {' '.join(cmd)}")
            prev_job_id = "DRY_RUN_JOB_ID"
            submitted.append((name, "DRY"))
            continue

        result = subprocess.run(cmd, capture_output=True, text=True,
                                cwd=BASE)
        # Parse job ID from "STELLA job submitted: XXXXX"
        job_id = None
        for line in result.stdout.splitlines():
            if "STELLA job submitted:" in line:
                job_id = line.split(":")[-1].strip()
                break

        if job_id:
            print(f"  Submitted: job {job_id}")
            prev_job_id = job_id
            submitted.append((name, job_id))
        else:
            print(f"  ERROR: could not parse job ID")
            print(f"  stdout: {result.stdout[-500:]}")
            print(f"  stderr: {result.stderr[-500:]}")
            sys.exit(1)

    print(f"\n{'='*60}")
    print(f"All {len(submitted)} models submitted with dependencies:")
    for name, jid in submitted:
        print(f"  {name:<30} -> job {jid}")
    if not dry_run:
        print(f"\nMonitor: squeue -u $USER")
        print(f"Chain: {' -> '.join(j for _, j in submitted)}")

if __name__ == "__main__":
    main()
