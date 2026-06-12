#!/bin/bash
# Submit CMFGEN jobs (init, final, obs) with proper dependencies

set -euo pipefail

# Parse optional --debug flag
DEBUG_FLAG=""
if [[ "$1" == "--debug" ]]; then
    DEBUG_FLAG="--debug"
    shift
fi

MODEL_NAME="$1"
MODEL_DIR="$2"
CMFGEN_DIR="$3"
HERACLES_JOB_ID="$4"
INIT_NUM_ITS="${5:-3}"
PARTITION="${6:-general}"
TIME="${7:-24:00:00}"
SINGLE_FINAL="${8:-1}"

# Export environment variables for submit_cmfgen.sh
export CMFGEN_INIT_PARTITION="${PARTITION}"
export CMFGEN_INIT_TIME="${TIME}"

# Export model name for use in submit_cmfgen.sh
export MODEL_NAME
# Export a unique job name prefix based on model name
export CMFGEN_JOB_PREFIX="${MODEL_NAME}"

# Derive BASE_DIR from MODEL_DIR (MODEL_DIR is <project>/models/<model>)
# Go up 2 directories to get to project root
BASE_DIR="$(dirname "$(dirname "$MODEL_DIR")")"
echo "DEBUG: MODEL_DIR=$MODEL_DIR"
echo "DEBUG: BASE_DIR=$BASE_DIR"
echo "DEBUG: Expected submit_cmfgen.sh path: ${BASE_DIR}/base/submit_cmfgen.sh"

if [[ $# -lt 4 ]]; then
    echo "Usage: $0 <MODEL_NAME> <MODEL_DIR> <CMFGEN_DIR> <HERACLES_JOB_ID> [INIT_NUM_ITS] [PARTITION] [TIME] [SINGLE_FINAL]"
    exit 1
fi

echo "======================================================"
echo "Submitting CMFGEN jobs for model: ${MODEL_NAME}"
echo "======================================================"
echo "Model directory: ${MODEL_DIR}"
echo "CMFGEN directory: ${CMFGEN_DIR}"
echo "Dependent on HERACLES job: ${HERACLES_JOB_ID}"
echo "Base directory: ${BASE_DIR}"
echo "INIT_NUM_ITS: ${INIT_NUM_ITS}"
echo "SINGLE_FINAL: ${SINGLE_FINAL} (1=skip INIT, LTE cold start; 0=INIT->FINAL chain)"
echo "Temperature control: always allowed to vary (FIX_T=F, FIX_X_DPTH=0)"
echo "======================================================"

cd "${CMFGEN_DIR}"

# Update local submit_cmfgen.sh with INIT_NUM_ITS from config
if [[ -f "submit_cmfgen.sh" ]]; then
    echo "Updating local submit_cmfgen.sh with INIT_NUM_ITS=${INIT_NUM_ITS}"
    sed -i "s/^INIT_NUM_ITS=.*/INIT_NUM_ITS=${INIT_NUM_ITS}/" submit_cmfgen.sh
fi

# Call submit_cmfgen.sh to generate and submit all CMFGEN jobs
echo "Calling submit_cmfgen.sh to generate and submit CMFGEN jobs..."
echo "DEBUG: Attempting to call: ${BASE_DIR}/base/submit_cmfgen.sh"
echo "DEBUG: Current directory: $(pwd)"
echo "DEBUG: Checking for local submit_cmfgen.sh..."
# Build the single_final flag to pass to submit_cmfgen.sh
if [[ "${SINGLE_FINAL}" == "1" ]]; then
    SINGLE_FINAL_FLAG="--single-final"
else
    SINGLE_FINAL_FLAG="--with-init"
fi

# Always use the base copy to ensure latest version with --dependency support
if [[ ! -f "${BASE_DIR}/base/submit_cmfgen.sh" ]]; then
    echo "ERROR: submit_cmfgen.sh not found at ${BASE_DIR}/base/submit_cmfgen.sh"
    echo "Looking for base directory..."
    for d in "${BASE_DIR}/base" "${BASE_DIR}/../base" "${BASE_DIR}/../../base"; do
        if [[ -f "$d/submit_cmfgen.sh" ]]; then
            echo "Found at: $d/submit_cmfgen.sh"
            BASE_DIR="$d"
            break
        fi
    done
fi
"${BASE_DIR}/base/submit_cmfgen.sh" \
    ${DEBUG_FLAG} \
    --dependency="${HERACLES_JOB_ID}" \
    ${SINGLE_FINAL_FLAG} \
    "${CMFGEN_DIR}"

# Parse job IDs from submit.log
FINAL_JOB=$(grep "FINAL Job ID" submit.log | awk '{print $4}' | head -1)
OBS_JOB=$(grep "BATOBS Job ID" submit.log | awk '{print $4}' | head -1)
if [[ "${SINGLE_FINAL}" == "1" ]]; then
    INIT_JOB="(none)"
else
    INIT_JOB=$(grep "INIT Job ID" submit.log | awk '{print $NF}' | head -1)
fi

echo ""
echo "======================================================"
echo "CMFGEN jobs submitted successfully"
echo "======================================================"
echo "INIT:  ${INIT_JOB}"
echo "FINAL: ${FINAL_JOB}"
echo "OBS:   ${OBS_JOB}"
echo "======================================================"

# Output job IDs for dependency tracking
echo "${INIT_JOB}:${FINAL_JOB}:${OBS_JOB}"
