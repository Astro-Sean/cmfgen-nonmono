#!/bin/bash
# Submit Ni56 grid models through run_model_grid.py with dependency chaining.
# Only 1 STELLA model runs at a time: job 2 waits for job 1, job 3 waits for job 2, etc.
set -euo pipefail

cd /viper/ptmp1/sbrennan/SN2024pba

# Ordered list of Ni56 grid configs (lowest Ni56 first)
CONFIGS=(
  SN_A_g_T20k_Ni04/config.yaml
  SN_A_g_T20k_Ni05/config.yaml
  SN_A_g_T20k_Ni06/config.yaml
  SN_A_g_T20k_Ni075/config.yaml
  SN_A_g_T20k_Ni09/config.yaml
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

  # For models after the first, add dependency on the previous job
  if [ $i -ge 1 ]; then
    dep_idx=$((i - 1))
    dep_job="${JOB_IDS[$dep_idx]}"
    echo "  Setting dependency: afterok:$dep_job (wait for job $dep_job)"
    scontrol update jobid="$job_id" dependency=afterok:"$dep_job"
  fi

  echo ""
done

echo "=============================================="
echo "Ni56 grid submitted with dependency chaining (1 at a time)"
echo "=============================================="
echo ""
echo "Job chain (each waits for the previous job):"
for i in "${!CONFIGS[@]}"; do
  model=$(dirname "${CONFIGS[$i]}")
  job_id="${JOB_IDS[$i]}"
  if [ $i -ge 1 ]; then
    dep_idx=$((i - 1))
    dep="${JOB_IDS[$dep_idx]}"
    echo "  $model -> $job_id (afterok:$dep)"
  else
    echo "  $model -> $job_id (runs immediately)"
  fi
done
