#!/bin/bash
# Submit grid models through run_model_grid.py with dependency chaining.
# Only 2 STELLA models run at a time: jobs 3,4 wait for 1,2; jobs 5,6 wait for 3,4; etc.
set -euo pipefail

cd /viper/ptmp1/sbrennan/SN2024pba

# Ordered list of grid configs
CONFIGS=(
  SN_A_g_M02/config.yaml
  SN_A_g_M05/config.yaml
  SN_A_g_R050/config.yaml
  SN_A_g_R100/config.yaml
  SN_A_g_T20k/config.yaml
  SN_A_g_T35k/config.yaml
  SN_A_g_E06/config.yaml
  SN_A_g_E12/config.yaml
  SN_A_g_N02/config.yaml
  SN_A_g_N06/config.yaml
  SN_A_g_H04/config.yaml
  SN_A_g_Ergon_lo/config.yaml
  SN_A_g_Compact/config.yaml
  SN_A_g_E12_M03/config.yaml
)

declare -a JOB_IDS

for i in "${!CONFIGS[@]}"; do
  cfg="${CONFIGS[$i]}"
  model=$(dirname "$cfg")
  echo "=== [$((i+1))/${#CONFIGS[@]}] Submitting $model ==="

  # Submit through orchestrator and capture the STELLA job ID
  output=$(python scripts/run_model_grid.py "$cfg" 2>&1)
  job_id=$(echo "$output" | grep -oP "STELLA job submitted: \K\d+" | head -1)

  if [ -z "$job_id" ]; then
    echo "ERROR: Failed to submit $model"
    echo "$output" | tail -20
    exit 1
  fi

  JOB_IDS[$i]=$job_id
  echo "  STELLA job: $job_id"

  # For models after the first 2, add dependency on the pair 2 positions back
  if [ $i -ge 2 ]; then
    dep_idx=$((i - 2))
    dep_job="${JOB_IDS[$dep_idx]}"
    echo "  Setting dependency: afterok:$dep_job (wait for job $dep_job)"
    scontrol update jobid="$job_id" dependency=afterok:"$dep_job"
  fi

  echo ""
done

echo "=============================================="
echo "Grid submitted with dependency chaining (2 at a time)"
echo "=============================================="
echo ""
echo "Job chain (each waits for the job 2 positions before it):"
for i in "${!CONFIGS[@]}"; do
  model=$(dirname "${CONFIGS[$i]}")
  job_id="${JOB_IDS[$i]}"
  if [ $i -ge 2 ]; then
    dep_idx=$((i - 2))
    dep="${JOB_IDS[$dep_idx]}"
    echo "  $model -> $job_id (afterok:$dep)"
  else
    echo "  $model -> $job_id (runs immediately)"
  fi
done
