#!/bin/bash
# Submit broad SN_B grid models through run_model_grid.py with a hard cap on
# concurrent STELLA jobs, enforced by SLURM afterany dependencies (lanes).
#
# Lane rule: new STELLA job i is submitted with --dependency afterany:J where J
# is the job MAXN positions earlier in submission order. The lane list is seeded
# with any STELLA jobs from the grid that are ALREADY queued/running, so the cap
# holds across separate invocations of this script.
#
# Submission order is RANDOMISED (user directive 22-Sep-2026): which of the
# queued models run first is left to shuffling, not directory order.
#
# Skip rules (user directive 21-Sep-2026):
#   - model already has stella/<name>.tt            -> skip (do not resubmit)
#   - model already has a queued/running STELLA job -> skip submission, but the
#     job id still occupies a lane slot
#
# Usage:
#   scripts/submit_broad_grid.sh [GRID_DIR ...]      # default: SN_B_broad_pilot
#   MAXN=4 MAXSUB=280 scripts/submit_broad_grid.sh SN_B_broad_grid
#
# MAXSUB stops after that many new submissions (assoc MaxSubmit=300 on the
# queue - pending dependency-held jobs still count).
set -euo pipefail

cd /viper/ptmp1/sbrennan/SN2024pba

MAXN=${MAXN:-4}
MAXSUB=${MAXSUB:-0}   # 0 = unlimited
NSUB=0
GRID_DIRS=("$@")
if [ ${#GRID_DIRS[@]} -eq 0 ]; then
  GRID_DIRS=(SN_B_broad_pilot)
fi

# Currently active SLURM jobs for this user: "<jobid> <name>"
mapfile -t ACTIVE < <(squeue -u "$USER" -h -o "%i %j" | awk '{print $1, $2}')
NACTIVE=$(squeue -u "$USER" -h | wc -l)
active_job_id() {
  local name="$1"
  [ "$NACTIVE" -eq 0 ] && return 1
  for entry in "${ACTIVE[@]}"; do
    if [ "${entry#* }" = "${name}_stella" ]; then
      echo "${entry%% *}"
      return 0
    fi
  done
  return 1
}

declare -a LANES   # job ids in submission order (seeds + new submissions)
NLANES=0

# Pass 1: seed the lane list with every active grid job BEFORE submitting
# anything. If seeds were interleaved with new submissions in scan order,
# early submissions could find NLANES < MAXN and run alongside all the
# active jobs, breaching the cap.
for grid in "${GRID_DIRS[@]}"; do
  for dir in "$grid"/*/; do
    cfg="${dir}config.yaml"
    [ -f "$cfg" ] || continue
    mname=$(grep -m1 '^model_name:' "$cfg" | awk '{print $2}' | tr -d '"'"'"'')
    mname=${mname:-$(basename "$dir")}
    if jid=$(active_job_id "$mname"); then
      echo "ACTIVE $mname (job $jid occupies a lane)"
      LANES+=("$jid")
      NLANES=$((NLANES + 1))
    fi
  done
done

# Pass 2: iterate configs in random order; skip .tt-complete and active.
mapfile -t CFGS < <(
  for grid in "${GRID_DIRS[@]}"; do
    for dir in "$grid"/*/; do
      [ -f "${dir}config.yaml" ] && printf '%s\n' "${dir}config.yaml"
    done
  done | shuf
)
echo "Candidate configs: ${#CFGS[@]}"

for cfg in "${CFGS[@]}"; do
  dir=$(dirname "$cfg")/
  model=$(basename "$dir")
  mname=$(grep -m1 '^model_name:' "$cfg" | awk '{print $2}' | tr -d '"'"'"'')
  mname=${mname:-$model}

  # Rule 1: .tt already exists -> do not resubmit
  if [ -f "${dir}stella/${mname}.tt" ]; then
    continue
  fi

  # Rule 2: already queued/running -> lane already seeded in pass 1
  if active_job_id "$mname" >/dev/null; then
    continue
  fi

  dep_arg=()
  dep_job=""
  if [ "$NLANES" -ge "$MAXN" ]; then
    dep_idx=$(( NLANES - MAXN ))
    dep_job="${LANES[$dep_idx]}"
    dep_arg=(--stella-dependency "afterany:${dep_job}")
  fi

  echo "=== Submitting $model ${dep_job:+[afterany:$dep_job]} ==="
  output=$(python scripts/run_model_grid.py "$cfg" ${dep_arg[@]+"${dep_arg[@]}"} 2>&1) || {
    echo "ERROR: run_model_grid.py failed for $model"
    echo "$output" | tail -20
    continue
  }
  job_id=$(echo "$output" | grep -oP "STELLA job submitted: \K\d+" | head -1)
  if [ -z "$job_id" ]; then
    echo "ERROR: no job id for $model"
    echo "$output" | tail -15
    continue
  fi
  LANES+=("$job_id")
  NLANES=$((NLANES + 1))
  NSUB=$((NSUB + 1))
  echo "  STELLA job: $job_id ${dep_job:+afterany:$dep_job}  (submission $NSUB)"
  if [ "$MAXSUB" -gt 0 ] && [ "$NSUB" -ge "$MAXSUB" ]; then
    echo "Reached MAXSUB=$MAXSUB - stopping."
    break
  fi
done

echo "=============================================="
echo "Submission complete ($NSUB new). Lane order (cap=$MAXN):"
if [ "$NLANES" -gt 0 ]; then
  for i in "${!LANES[@]}"; do
    echo "  [$i] ${LANES[$i]}"
  done
fi
