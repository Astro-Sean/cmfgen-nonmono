#!/usr/bin/env bash
# watchme_status.sh — print a single SLURM status frame for `watchme`.
#
# Called by:  watch -n <interval> watchme_status.sh [name_width]
# (the watchme bash function in ~/.bashrc wraps this with `watch`).
#
# Output:
#   - free CPUs per partition (idle/total)
#   - RUNNING jobs (sorted by time)
#   - PENDING jobs (sorted by JobID)
set -euo pipefail

name_width="${1:-40}"
user="${SLURM_JOB_USER:-$(whoami)}"
fmt="%.10i %.${name_width}j %.12P %.10T %.12M %.20R"

echo "══════════════════════════════════════════════════════════════════════"
echo "  SLURM Job Monitor — user: ${user}  —  $(date '+%Y-%m-%d %H:%M:%S')"
echo "══════════════════════════════════════════════════════════════════════"
echo ""

# Running jobs
echo "── RUNNING ──────────────────────────────────────────────────────────"
running=$(squeue -u "$user" -t RUNNING -h -o "$fmt" 2>/dev/null | sort -k5 -r)
if [[ -n "$running" ]]; then
    printf "  %-10s %-${name_width}s %-12s %-10s %-12s %s\n" \
           "JOBID" "NAME" "PARTITION" "STATE" "TIME" "REASON/NODE"
    echo "$running" | while IFS= read -r line; do echo "  $line"; done
else
    echo "  (none)"
fi
echo ""

# Pending jobs
echo "── PENDING ──────────────────────────────────────────────────────────"
pending=$(squeue -u "$user" -t PENDING -h -o "$fmt" 2>/dev/null | sort -k1 -n)
if [[ -n "$pending" ]]; then
    printf "  %-10s %-${name_width}s %-12s %-10s %-12s %s\n" \
           "JOBID" "NAME" "PARTITION" "STATE" "TIME" "REASON"
    echo "$pending" | while IFS= read -r line; do echo "  $line"; done
else
    echo "  (none)"
fi
echo ""
echo "─────────────────────────────────────────────────────────────────────"
echo "  Usage: watchme [name_width] [interval]"
