#!/bin/bash
# Submit SN_B grid models through run_model_grid.py with dependency chaining.
# Only 2 STELLA models run at a time: jobs 3,4 wait for 1,2; jobs 5,6 wait for 3,4; etc.
#
# Fault-tolerant: if one model fails to submit or its STELLA job fails, the
# remaining models still run. Dependencies use afterany (not afterok) so a
# failed job does not cancel downstream jobs.
#
# Grid: M{3,5,7,10} x V{3000,5000,7000,10000} x {wind, no-wind}
#   4 masses x 4 velocities x 2 wind settings = 32 models
# Max velocity: 10,000 km/s
# STELLA evolution: 60 days
# Codes: STELLA + SuperLite
#
# Already complete (STELLA + SuperLite done, skipped):
#   M3.0_V3000, M3.0_V3000_NW, M3.0_V5000
#
# STELLA backup restored (orchestrator auto-detects complete .swd, SuperLite only):
#   M3.0_V7000
set -uo pipefail

cd /viper/ptmp1/sbrennan/SN2024pba

# Ordered list of SN_B grid configs (28 models remaining)
# Order: mass-major, then velocity, then wind before no-wind.
# Skipped: M3.0_V3000, M3.0_V3000_NW, M3.0_V5000 (STELLA + SuperLite already done)
# M3.0_V7000 has a restored backup .swd — orchestrator will auto-skip STELLA.
CONFIGS=(
  SN_B_grid/M3.0_V5000_NW/config.yaml
  SN_B_grid/M3.0_V7000/config.yaml
  SN_B_grid/M3.0_V7000_NW/config.yaml
  SN_B_grid/M3.0_V10000/config.yaml
  SN_B_grid/M3.0_V10000_NW/config.yaml
  SN_B_grid/M5.0_V3000/config.yaml
  SN_B_grid/M5.0_V3000_NW/config.yaml
  SN_B_grid/M5.0_V5000/config.yaml
  SN_B_grid/M5.0_V5000_NW/config.yaml
  SN_B_grid/M5.0_V7000/config.yaml
  SN_B_grid/M5.0_V7000_NW/config.yaml
  SN_B_grid/M5.0_V10000/config.yaml
  SN_B_grid/M5.0_V10000_NW/config.yaml
  SN_B_grid/M7.0_V3000/config.yaml
  SN_B_grid/M7.0_V3000_NW/config.yaml
  SN_B_grid/M7.0_V5000/config.yaml
  SN_B_grid/M7.0_V5000_NW/config.yaml
  SN_B_grid/M7.0_V7000/config.yaml
  SN_B_grid/M7.0_V7000_NW/config.yaml
  SN_B_grid/M7.0_V10000/config.yaml
  SN_B_grid/M7.0_V10000_NW/config.yaml
  SN_B_grid/M10.0_V3000/config.yaml
  SN_B_grid/M10.0_V3000_NW/config.yaml
  SN_B_grid/M10.0_V5000/config.yaml
  SN_B_grid/M10.0_V5000_NW/config.yaml
  SN_B_grid/M10.0_V7000/config.yaml
  SN_B_grid/M10.0_V7000_NW/config.yaml
  SN_B_grid/M10.0_V10000/config.yaml
  SN_B_grid/M10.0_V10000_NW/config.yaml
)

declare -a JOB_IDS
declare -a JOB_TYPES  # "stella" or "superlite" or "failed"
FAILED_COUNT=0
SUCCESS_COUNT=0

for i in "${!CONFIGS[@]}"; do
  cfg="${CONFIGS[$i]}"
  model=$(dirname "$cfg")
  echo "=== [$((i+1))/${#CONFIGS[@]}] Submitting $model ==="

  # Submit through orchestrator and capture output
  output=$(python scripts/run_model_grid.py "$cfg" 2>&1)
  rc=$?

  # Try to find a STELLA job ID first
  job_id=$(echo "$output" | grep -oP "STELLA job submitted: \K\d+" | head -1)
  job_type="stella"

  # If no STELLA job, check for SuperLite job (e.g. M3.0_V7000 with restored .swd)
  if [ -z "$job_id" ]; then
    job_id=$(echo "$output" | grep -oP "SuperLite job submitted: \K\d+" | head -1)
    job_type="superlite"
  fi

  if [ -z "$job_id" ]; then
    echo "  WARNING: Failed to submit $model (orchestrator rc=$rc)"
    echo "$output" | tail -20
    JOB_IDS[$i]=""
    JOB_TYPES[$i]="failed"
    FAILED_COUNT=$((FAILED_COUNT + 1))
    echo ""
    continue
  fi

  JOB_IDS[$i]=$job_id
  JOB_TYPES[$i]=$job_type
  SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
  echo "  $job_type job: $job_id"

  # For STELLA jobs after the first 2, add dependency on the pair 2 positions back.
  # This limits concurrency to 2 STELLA jobs at a time.
  # Use afterany (not afterok) so a failed job does not cancel downstream jobs.
  if [ "$job_type" = "stella" ] && [ $i -ge 2 ]; then
    dep_idx=$((i - 2))
    dep_job="${JOB_IDS[$dep_idx]}"
    dep_type="${JOB_TYPES[$dep_idx]}"

    if [ -n "$dep_job" ] && [ "$dep_type" = "stella" ]; then
      echo "  Setting dependency: afterany:$dep_job (wait for $dep_type job $dep_job)"
      scontrol update jobid="$job_id" dependency=afterany:"$dep_job" 2>/dev/null || \
        echo "  WARNING: Failed to set dependency on $dep_job"
    else
      echo "  No STELLA dependency to chain on (previous was $dep_type or empty) — runs immediately"
    fi
  fi

  echo ""
done

echo "=============================================="
echo "SN_B grid submission complete"
echo "Grid: M{3,5,7,10} x V{3000,5000,7000,10000} x {wind, NW}"
echo "  Submitted: $SUCCESS_COUNT models"
echo "  Failed:    $FAILED_COUNT models"
echo "  Skipped:   3 models (M3.0_V3000, M3.0_V3000_NW, M3.0_V5000 — already complete)"
echo "STELLA evolution: 60 days | Codes: STELLA + SuperLite"
echo "Concurrency: 2 STELLA jobs at a time (afterany dependencies)"
echo "=============================================="
echo ""
echo "Job chain:"
for i in "${!CONFIGS[@]}"; do
  model=$(dirname "${CONFIGS[$i]}")
  job_id="${JOB_IDS[$i]}"
  job_type="${JOB_TYPES[$i]}"
  if [ "$job_type" = "failed" ]; then
    echo "  $model -> FAILED TO SUBMIT"
  elif [ $i -ge 2 ] && [ "$job_type" = "stella" ]; then
    dep_idx=$((i - 2))
    dep="${JOB_IDS[$dep_idx]}"
    dep_type="${JOB_TYPES[$dep_idx]}"
    if [ -n "$dep" ] && [ "$dep_type" = "stella" ]; then
      echo "  $model -> $job_id ($job_type, afterany:$dep)"
    else
      echo "  $model -> $job_id ($job_type, runs immediately)"
    fi
  else
    echo "  $model -> $job_id ($job_type, runs immediately)"
  fi
done
