#!/bin/bash
# Submit CMFGEN jobs (init, final, obs) with proper dependencies

set -euo pipefail

# Parse optional --debug flag
DEBUG_FLAG=""
if [[ "$1" == "--debug" ]]; then
    DEBUG_FLAG="--debug"
    shift
fi

if [[ $# -lt 4 ]]; then
    echo "Usage: $0 [--debug] <MODEL_NAME> <MODEL_DIR> <CMFGEN_DIR> <HERACLES_JOB_ID> [INIT_NUM_ITS] [PARTITION] [TIME] [SINGLE_FINAL] [CONFIG_FILE] [INC_RAD_DECAYS] [VEL_LAW] [NO_RELATIVISTIC]"
    exit 1
fi

MODEL_NAME="$1"
MODEL_DIR="$2"
CMFGEN_DIR="$3"
HERACLES_JOB_ID="$4"
INIT_NUM_ITS="${5:-3}"
PARTITION="${6:-general}"
TIME="${7:-24:00:00}"
SINGLE_FINAL="${8:-1}"
CONFIG_FILE="${9:-}"
INC_RAD_DECAYS="${10:-0}"
VEL_LAW="${11:-12}"
NO_RELATIVISTIC="${12:-0}"

# Export environment variables for submit_cmfgen.sh
export CMFGEN_INIT_PARTITION="${PARTITION}"
export CMFGEN_INIT_TIME="${TIME}"
export CMFGEN_INIT_NUM_ITS="${INIT_NUM_ITS}"

# Export model name for use in submit_cmfgen.sh
export MODEL_NAME
# Export a unique job name prefix based on model name
export CMFGEN_JOB_PREFIX="${MODEL_NAME}"
# Export config file path (avoids hardcoded "model_" prefix assumption)
# Only export if CONFIG_FILE is non-empty
if [[ -n "${CONFIG_FILE}" ]]; then
  export CMFGEN_CONFIG_FILE="${CONFIG_FILE}"
  echo "Config file path: ${CMFGEN_CONFIG_FILE}"
else
  echo "Config file path not provided; will use default model_${MODEL_NAME}.yaml pattern"
fi

# Export radioactive decays flag for prepare_VADAT.py
export CMFGEN_INC_RAD_DECAYS="${INC_RAD_DECAYS}"

# Export no-relativistic flag for prepare_VADAT.py
export CMFGEN_NO_RELATIVISTIC="${NO_RELATIVISTIC}"

# Derive BASE_DIR from MODEL_DIR (MODEL_DIR is <project>/models/<model>/cmfgen)
# Use the actual project root path
BASE_DIR="/viper/ptmp1/sbrennan/SN2024pba/cmfgen_base"
echo "======================================================"
echo "Submitting CMFGEN jobs for model: ${MODEL_NAME}"
echo "======================================================"
echo "Model directory: ${MODEL_DIR}"
echo "CMFGEN directory: ${CMFGEN_DIR}"
echo "Dependent on HERACLES job: ${HERACLES_JOB_ID}"
echo "Base directory: ${BASE_DIR}"
echo "INIT_NUM_ITS: ${INIT_NUM_ITS}"
echo "SINGLE_FINAL: ${SINGLE_FINAL} (1=skip INIT, LTE cold start; 0=INIT->FINAL chain)"
echo "VEL_LAW: ${VEL_LAW} (12=IIn tau-based grid; 11=direct hydro R/V/SIGMA)"
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
# Build the single_final flag to pass to submit_cmfgen.sh
if [[ "${SINGLE_FINAL}" == "1" ]]; then
    SINGLE_FINAL_FLAG="--single-final"
else
    SINGLE_FINAL_FLAG="--with-init"
fi

# Always use the base copy to ensure latest version with --dependency support
if [[ ! -f "${BASE_DIR}/submit_cmfgen.sh" ]]; then
    echo "ERROR: submit_cmfgen.sh not found at ${BASE_DIR}/submit_cmfgen.sh"
    echo "Looking for base directory..."
    for d in "${BASE_DIR}" "${BASE_DIR}/../cmfgen_base" "${BASE_DIR}/../../cmfgen_base"; do
        if [[ -f "$d/submit_cmfgen.sh" ]]; then
            echo "Found at: $d/submit_cmfgen.sh"
            BASE_DIR="$d"
            break
        fi
    done
fi
# Ensure CMFGEN_JOB_PREFIX is passed through to submit_cmfgen.sh
env CMFGEN_JOB_PREFIX="${CMFGEN_JOB_PREFIX}" "${BASE_DIR}/submit_cmfgen.sh" \
    ${DEBUG_FLAG} \
    --dependency="${HERACLES_JOB_ID}" \
    ${SINGLE_FINAL_FLAG} \
    --vel-law="${VEL_LAW}" \
    "${CMFGEN_DIR}"

# Parse job IDs from submit.log (in CMFGEN_DIR where submit_cmfgen.sh writes it)
FINAL_JOB=$(grep "FINAL Job ID" "${CMFGEN_DIR}/submit.log" | awk '{print $4}' | head -1)
OBS_JOB=$(grep "BATOBS Job ID" "${CMFGEN_DIR}/submit.log" | awk '{print $4}' | head -1)
if [[ "${SINGLE_FINAL}" == "1" ]]; then
    INIT_JOB="(none)"
else
    INIT_JOB=$(grep "INIT Job ID" "${CMFGEN_DIR}/submit.log" | awk '{print $NF}' | head -1)
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
