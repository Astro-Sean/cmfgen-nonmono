#!/bin/bash
# Submit CMFGEN jobs (final + obs) with proper dependencies.
# INIT/STABLE stages have been removed — only single_final (LTE cold start)
# or final (warm restart from previous day) are used.
# All configuration is passed as positional arguments from run_model_grid.py (which reads YAML).

set -euo pipefail

# Parse optional --debug flag
DEBUG_FLAG=""
if [[ "$1" == "--debug" ]]; then
    DEBUG_FLAG="--debug"
    shift
fi

if [[ $# -lt 4 ]]; then
    echo "Usage: $0 [--debug] <MODEL_NAME> <MODEL_DIR> <CMFGEN_DIR> <HERACLES_JOB_ID> [INIT_NUM_ITS] [PARTITION] [TIME] [SINGLE_FINAL] [CONFIG_FILE] [INC_RAD_DECAYS] [VEL_LAW] [NO_RELATIVISTIC] [MODEL_TYPE] [DAY_LABEL] [PREV_CMFGEN_DIR] [TARGET_ND] [SLURM_MEM] [INIT_SLURM_MEM] [N_OB_INS] [FINAL_REQUEUE_WHEN] [FINAL_REQUEUE_MAX] [AUTO_SET_RSTAR] [AUTO_SET_LSTAR] [PHOT_TAU] [BB_TEFF_K] [LSTAR_RADIUS_MODE] [PROXY_POPULATIONS] [PROXY_POPULATIONS_MIN_ION_FRAC] [PROXY_POPULATIONS_NLTE_SOLVER] [PROXY_POPULATIONS_GREY_MAX_ITER] [PROXY_POPULATIONS_GREY_TOL] [OBS_TIME] [USE_HYDRO_GRID] [SMOOTH_VELOCITY_JUMPS] [SMOOTH_MAX_DV_RATIO] [SMOOTH_N_INSERT] [SMOOTH_RHO_SIGMA] [SMOOTH_TEMP_SIGMA] [BLEND_COMPOSITION] [BLEND_SIGMA] [BLEND_SOLAR_WEIGHT] [MIN_MASS_FRACTION] [MONOTONIZE_VELOCITY] [CHEMICAL_N_INSERT] [CHEMICAL_JUMP_THRESHOLD] [REVERSAL_N_INSERT] [PROXY_POPULATIONS_T_EXC_MIN]"
    exit 1
fi

MODEL_NAME="$1"
MODEL_DIR="$2"
CMFGEN_DIR="$3"
HERACLES_JOB_ID="$4"
INIT_NUM_ITS="${5:-3}"
PARTITION="${6:-small}"
TIME="${7:-23:59:00}"
SINGLE_FINAL="${8:-1}"
CONFIG_FILE="${9:-}"
INC_RAD_DECAYS="${10:-0}"
VEL_LAW="${11:-12}"
NO_RELATIVISTIC="${12:-0}"
MODEL_TYPE="${13:-supernova}"
DAY_LABEL="${14:-}"
PREV_CMFGEN_DIR="${15:-}"
TARGET_ND="${16:-}"
SLURM_MEM="${17:-16G}"
INIT_SLURM_MEM="${18:-16G}"
N_OB_INS="${19:-30}"
FINAL_REQUEUE_WHEN="${20:-wall}"
FINAL_REQUEUE_MAX="${21:-500}"
AUTO_SET_RSTAR="${22:-0}"
AUTO_SET_LSTAR="${23:-0}"
PHOT_TAU="${24:-0.6666667}"
BB_TEFF_K="${25:-9100}"
LSTAR_RADIUS_MODE="${26:-RMAX}"
PROXY_POPULATIONS="${27:-0}"
PROXY_POPULATIONS_MIN_ION_FRAC="${28:-}"
PROXY_POPULATIONS_NLTE_SOLVER="${29:-0}"
PROXY_POPULATIONS_GREY_MAX_ITER="${30:-}"
PROXY_POPULATIONS_GREY_TOL="${31:-}"
OBS_TIME="${32:-}"
USE_HYDRO_GRID="${33:-0}"
SMOOTH_VELOCITY_JUMPS="${34:-0}"
SMOOTH_MAX_DV_RATIO="${35:-2.0}"
SMOOTH_N_INSERT="${36:-8}"
SMOOTH_RHO_SIGMA="${37:-0.0}"
SMOOTH_TEMP_SIGMA="${38:-0.0}"
BLEND_COMPOSITION="${39:-0}"
BLEND_SIGMA="${40:-5.0}"
BLEND_SOLAR_WEIGHT="${41:-0.0}"
MIN_MASS_FRACTION="${42:-0.0}"
MONOTONIZE_VELOCITY="${43:-0}"
CHEMICAL_N_INSERT="${44:-0}"
CHEMICAL_JUMP_THRESHOLD="${45:-0.1}"
REVERSAL_N_INSERT="${46:-0}"
PROXY_POPULATIONS_T_EXC_MIN="${47:-}"
PROXY_POPULATIONS_V2="${48:-0}"
CPUS_PER_TASK="${49:-}"

# Named args appended after the positional args by run_model_grid.py
# (otherwise silently dropped -> node-side proxy falls back to defaults).
PROXY_IONIZATION_MODE=""
PROXY_T_EXC_SOURCE=""
PROXY_T_ION_FLOOR=""
PROXY_COMPLEX_MODEL=""
PROXY_POPULATIONS_MASS_FRAC_FLOOR=""
PROXY_DC_MAX=""
PROXY_EXCITATION_MODE=""
PROXY_DC_LOG_CAP=""
PROXY_NLTE_UPSILON=""
PROXY_NLTE_UPSILON_FORB=""
PROXY_NLTE_F_THRESHOLD=""
PROXY_NLTE_ELEMS=""
PROXY_STELLA_FILE=""
PROXY_STELLA_RADIUS_COLUMN=""
PROXY_STELLA_VELOCITY_COLUMN=""
PROXY_STELLA_TAU_COLUMN=""
PROXY_STELLA_TRAD_COLUMN=""
PROXY_STELLA_TION_COLUMN=""
PROXY_STELLA_LUMINOSITY_COLUMN=""
PROXY_STELLA_ENERGY_DENSITY_COLUMN=""
PROXY_STELLA_RADIUS_UNIT=""
PROXY_STELLA_TEMPERATURE_UNIT=""
PROXY_STELLA_INTERP_COORD=""
PROXY_V2_STRATEGY=""
PROXY_PREV_MODEL=""
PROXY_NE_MAX_ITER=""
PROXY_NE_TOL=""
PROXY_TAU_TRANS=""
PROXY_NO_NON_LOCAL_REVERSALS=""
GT_AUTO=""
for _arg in "$@"; do
  case "${_arg}" in
    --gt-auto=*) GT_AUTO="${_arg#*=}" ;;
    --proxy-ionization-mode=*) PROXY_IONIZATION_MODE="${_arg#*=}" ;;
    --proxy-t-exc-source=*) PROXY_T_EXC_SOURCE="${_arg#*=}" ;;
    --proxy-t-ion-floor=*) PROXY_T_ION_FLOOR="${_arg#*=}" ;;
    --proxy-complex-model=*) PROXY_COMPLEX_MODEL="${_arg#*=}" ;;
    --proxy-populations-mass-frac-floor=*) PROXY_POPULATIONS_MASS_FRAC_FLOOR="${_arg#*=}" ;;
    --proxy-dc-max=*) PROXY_DC_MAX="${_arg#*=}" ;;
    --proxy-excitation-mode=*) PROXY_EXCITATION_MODE="${_arg#*=}" ;;
    --proxy-dc-log-cap=*) PROXY_DC_LOG_CAP="${_arg#*=}" ;;
    --proxy-ne-max-iter=*) PROXY_NE_MAX_ITER="${_arg#*=}" ;;
    --proxy-ne-tol=*) PROXY_NE_TOL="${_arg#*=}" ;;
    --proxy-tau-trans=*) PROXY_TAU_TRANS="${_arg#*=}" ;;
    --proxy-no-non-local-reversals=*) PROXY_NO_NON_LOCAL_REVERSALS="${_arg#*=}" ;;
    --proxy-nlte-upsilon=*) PROXY_NLTE_UPSILON="${_arg#*=}" ;;
    --proxy-nlte-upsilon-forb=*) PROXY_NLTE_UPSILON_FORB="${_arg#*=}" ;;
    --proxy-nlte-f-threshold=*) PROXY_NLTE_F_THRESHOLD="${_arg#*=}" ;;
    --proxy-nlte-elems=*) PROXY_NLTE_ELEMS="${_arg#*=}" ;;
    --proxy-stella-file=*) PROXY_STELLA_FILE="${_arg#*=}" ;;
    --proxy-stella-radius-column=*) PROXY_STELLA_RADIUS_COLUMN="${_arg#*=}" ;;
    --proxy-stella-velocity-column=*) PROXY_STELLA_VELOCITY_COLUMN="${_arg#*=}" ;;
    --proxy-stella-tau-column=*) PROXY_STELLA_TAU_COLUMN="${_arg#*=}" ;;
    --proxy-stella-trad-column=*) PROXY_STELLA_TRAD_COLUMN="${_arg#*=}" ;;
    --proxy-stella-tion-column=*) PROXY_STELLA_TION_COLUMN="${_arg#*=}" ;;
    --proxy-stella-luminosity-column=*) PROXY_STELLA_LUMINOSITY_COLUMN="${_arg#*=}" ;;
    --proxy-stella-energy-density-column=*) PROXY_STELLA_ENERGY_DENSITY_COLUMN="${_arg#*=}" ;;
    --proxy-stella-radius-unit=*) PROXY_STELLA_RADIUS_UNIT="${_arg#*=}" ;;
    --proxy-stella-temperature-unit=*) PROXY_STELLA_TEMPERATURE_UNIT="${_arg#*=}" ;;
    --proxy-stella-interpolation-coordinate=*) PROXY_STELLA_INTERP_COORD="${_arg#*=}" ;;
    --proxy-v2-strategy=*) PROXY_V2_STRATEGY="${_arg#*=}" ;;
    --proxy-prev-model=*) PROXY_PREV_MODEL="${_arg#*=}" ;;
    --proxy-*) echo "ERROR: unknown proxy option: ${_arg}" >&2; exit 2 ;;
  esac
done

BLEND_FLAG=""
if [[ "${BLEND_COMPOSITION}" == "1" ]]; then
  BLEND_FLAG="--blend-composition --blend-sigma=${BLEND_SIGMA} --blend-solar-weight=${BLEND_SOLAR_WEIGHT}"
  if [[ "${MIN_MASS_FRACTION}" != "0" && "${MIN_MASS_FRACTION}" != "0.0" ]]; then
    BLEND_FLAG="${BLEND_FLAG} --min-mass-fraction=${MIN_MASS_FRACTION}"
  fi
fi
MONOTONIZE_FLAG=""
if [[ "${MONOTONIZE_VELOCITY}" == "1" ]]; then
  MONOTONIZE_FLAG="--monotonize-velocity"
fi
CHEMICAL_FLAG=""
if [[ "${CHEMICAL_N_INSERT}" != "0" ]]; then
  CHEMICAL_FLAG="--chemical-n-insert=${CHEMICAL_N_INSERT} --chemical-jump-threshold=${CHEMICAL_JUMP_THRESHOLD}"
fi
REVERSAL_FLAG=""
if [[ "${REVERSAL_N_INSERT}" != "0" ]]; then
  REVERSAL_FLAG="--reversal-n-insert=${REVERSAL_N_INSERT}"
fi
USE_PROXY_V2="${USE_PROXY_V2:-0}"

# Derive BASE_DIR from MODEL_DIR (MODEL_DIR is <project>/models/<model>/cmfgen)
# Cluster-aware: detect RAVEN vs VIPER from hostname to avoid hardcoded paths.
_host="$(hostname 2>/dev/null || echo)"
case "$_host" in
  raven*|rav*) BASE_DIR="/raven/ptmp/sbrennan/SN2024pba/cmfgen_base" ;;
  *)            BASE_DIR="/viper/ptmp1/sbrennan/SN2024pba/cmfgen_base" ;;
esac
# Allow explicit override via environment
BASE_DIR="${CMFGEN_BASE_DIR:-$BASE_DIR}"
PROJECT_ROOT="$(dirname "$BASE_DIR")"

# Derive job prefix from model name and day label
if [[ -n "${DAY_LABEL}" ]]; then
  JOB_PREFIX="${MODEL_NAME}_${DAY_LABEL}"
else
  JOB_PREFIX="${MODEL_NAME}"
fi

echo "======================================================"
echo "Submitting CMFGEN jobs for model: ${MODEL_NAME}"
echo "======================================================"
echo "Model directory: ${MODEL_DIR}"
echo "CMFGEN directory: ${CMFGEN_DIR}"
echo "Dependent on HERACLES job: ${HERACLES_JOB_ID}"
echo "Base directory: ${BASE_DIR}"
echo "Mode: ${SINGLE_FINAL:+single_final (LTE cold start)}${PREV_CMFGEN_DIR:+restart (warm from previous day)}"
echo "VEL_LAW: ${VEL_LAW} (12=IIn tau-based grid; 11=direct hydro R/V/SIGMA)"
echo "MODEL_TYPE: ${MODEL_TYPE}"
if [[ "${PROXY_POPULATIONS}" == "1" ]]; then
  echo "PROXY_POPULATIONS: 1 (use proxy dilute-LTE _IN files instead of CMFGEN built-in LTE)"
  if [[ -n "${PROXY_POPULATIONS_MIN_ION_FRAC}" ]]; then
    echo "PROXY_POPULATIONS_MIN_ION_FRAC: ${PROXY_POPULATIONS_MIN_ION_FRAC}"
  fi
fi
echo "======================================================"

cd "${CMFGEN_DIR}"

# Build the restart flag: if PREV_CMFGEN_DIR is set, use warm restart (--restart);
# otherwise use single_final (LTE cold start, now the default in submit_cmfgen.sh).
if [[ -n "${PREV_CMFGEN_DIR}" ]]; then
    RESTART_FLAG="--restart"
else
    RESTART_FLAG=""
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

# Pass SKIP_OBS to submit_cmfgen.sh via environment (read from config)
export SKIP_OBS="${SKIP_OBS:-0}"
export SKIP_FINAL="${SKIP_FINAL:-0}"

# Pass all config as CLI args to submit_cmfgen.sh (no env vars)
"${BASE_DIR}/submit_cmfgen.sh" \
    ${DEBUG_FLAG} \
    --dependency="${HERACLES_JOB_ID}" \
    ${RESTART_FLAG} \
    --vel-law="${VEL_LAW}" \
    --use-hydro-grid="${USE_HYDRO_GRID}" \
    --smooth-velocity-jumps="${SMOOTH_VELOCITY_JUMPS}" \
    --smooth-max-dv-ratio="${SMOOTH_MAX_DV_RATIO}" \
    --smooth-n-insert="${SMOOTH_N_INSERT}" \
    --smooth-rho-sigma="${SMOOTH_RHO_SIGMA}" \
    --smooth-temp-sigma="${SMOOTH_TEMP_SIGMA}" \
    ${BLEND_FLAG} \
    ${MONOTONIZE_FLAG} \
    ${CHEMICAL_FLAG} \
    ${REVERSAL_FLAG} \
    --init-num-its="${INIT_NUM_ITS}" \
    ${GT_AUTO:+--gt-auto="${GT_AUTO}"} \
    --model-type="${MODEL_TYPE}" \
    --inc-rad-decays="${INC_RAD_DECAYS}" \
    --no-relativistic="${NO_RELATIVISTIC}" \
    --partition="${PARTITION}" \
    --time="${TIME}" \
    ${OBS_TIME:+--obs-time="${OBS_TIME}"} \
    --mem="${SLURM_MEM}" \
    --init-mem="${INIT_SLURM_MEM}" \
    ${CPUS_PER_TASK:+--cpus="${CPUS_PER_TASK}"} \
    --n-ob-ins="${N_OB_INS}" \
    --final-requeue-when="${FINAL_REQUEUE_WHEN}" \
    --final-requeue-max="${FINAL_REQUEUE_MAX}" \
    --auto-set-rstar="${AUTO_SET_RSTAR}" \
    --auto-set-lstar="${AUTO_SET_LSTAR}" \
    --phot-tau="${PHOT_TAU}" \
    --bb-teff-k="${BB_TEFF_K}" \
    --lstar-radius-mode="${LSTAR_RADIUS_MODE}" \
    --job-prefix="${JOB_PREFIX}" \
    --model-name="${MODEL_NAME}" \
    ${DAY_LABEL:+--day-label="${DAY_LABEL}"} \
    --project-root="${PROJECT_ROOT}" \
    ${CONFIG_FILE:+--config-file="${CONFIG_FILE}"} \
    ${TARGET_ND:+--target-nd="${TARGET_ND}"} \
    --proxy-populations="${PROXY_POPULATIONS}" \
    ${PROXY_POPULATIONS_MIN_ION_FRAC:+--proxy-populations-min-ion-frac="${PROXY_POPULATIONS_MIN_ION_FRAC}"} \
    ${PROXY_POPULATIONS_NLTE_SOLVER:+--proxy-populations-nlte-solver="${PROXY_POPULATIONS_NLTE_SOLVER}"} \
    ${PROXY_POPULATIONS_GREY_MAX_ITER:+--proxy-populations-grey-max-iter="${PROXY_POPULATIONS_GREY_MAX_ITER}"} \
    ${PROXY_POPULATIONS_GREY_TOL:+--proxy-populations-grey-tol="${PROXY_POPULATIONS_GREY_TOL}"} \
    ${PROXY_POPULATIONS_T_EXC_MIN:+--proxy-populations-t-exc-min="${PROXY_POPULATIONS_T_EXC_MIN}"} \
    --proxy-populations-v2="${PROXY_POPULATIONS_V2}" \
    ${PROXY_IONIZATION_MODE:+--proxy-ionization-mode="${PROXY_IONIZATION_MODE}"} \
    ${PROXY_T_EXC_SOURCE:+--proxy-t-exc-source="${PROXY_T_EXC_SOURCE}"} \
    ${PROXY_T_ION_FLOOR:+--proxy-t-ion-floor="${PROXY_T_ION_FLOOR}"} \
    ${PROXY_COMPLEX_MODEL:+--proxy-complex-model="${PROXY_COMPLEX_MODEL}"} \
    ${PROXY_POPULATIONS_MASS_FRAC_FLOOR:+--proxy-populations-mass-frac-floor="${PROXY_POPULATIONS_MASS_FRAC_FLOOR}"} \
    ${PROXY_DC_MAX:+--proxy-dc-max="${PROXY_DC_MAX}"} \
    ${PROXY_EXCITATION_MODE:+--proxy-excitation-mode="${PROXY_EXCITATION_MODE}"} \
    ${PROXY_DC_LOG_CAP:+--proxy-dc-log-cap="${PROXY_DC_LOG_CAP}"} \
    ${PROXY_NE_MAX_ITER:+--proxy-ne-max-iter="${PROXY_NE_MAX_ITER}"} \
    ${PROXY_NE_TOL:+--proxy-ne-tol="${PROXY_NE_TOL}"} \
    ${PROXY_TAU_TRANS:+--proxy-tau-trans="${PROXY_TAU_TRANS}"} \
    ${PROXY_NO_NON_LOCAL_REVERSALS:+--proxy-no-non-local-reversals="${PROXY_NO_NON_LOCAL_REVERSALS}"} \
    ${PROXY_NLTE_UPSILON:+--proxy-nlte-upsilon="${PROXY_NLTE_UPSILON}"} \
    ${PROXY_NLTE_UPSILON_FORB:+--proxy-nlte-upsilon-forb="${PROXY_NLTE_UPSILON_FORB}"} \
    ${PROXY_NLTE_F_THRESHOLD:+--proxy-nlte-f-threshold="${PROXY_NLTE_F_THRESHOLD}"} \
    ${PROXY_NLTE_ELEMS:+--proxy-nlte-elems="${PROXY_NLTE_ELEMS}"} \
    ${PROXY_STELLA_FILE:+--proxy-stella-file="${PROXY_STELLA_FILE}"} \
    ${PROXY_STELLA_RADIUS_COLUMN:+--proxy-stella-radius-column="${PROXY_STELLA_RADIUS_COLUMN}"} \
    ${PROXY_STELLA_VELOCITY_COLUMN:+--proxy-stella-velocity-column="${PROXY_STELLA_VELOCITY_COLUMN}"} \
    ${PROXY_STELLA_TAU_COLUMN:+--proxy-stella-tau-column="${PROXY_STELLA_TAU_COLUMN}"} \
    ${PROXY_STELLA_TRAD_COLUMN:+--proxy-stella-trad-column="${PROXY_STELLA_TRAD_COLUMN}"} \
    ${PROXY_STELLA_TION_COLUMN:+--proxy-stella-tion-column="${PROXY_STELLA_TION_COLUMN}"} \
    ${PROXY_STELLA_LUMINOSITY_COLUMN:+--proxy-stella-luminosity-column="${PROXY_STELLA_LUMINOSITY_COLUMN}"} \
    ${PROXY_STELLA_ENERGY_DENSITY_COLUMN:+--proxy-stella-energy-density-column="${PROXY_STELLA_ENERGY_DENSITY_COLUMN}"} \
    ${PROXY_STELLA_RADIUS_UNIT:+--proxy-stella-radius-unit="${PROXY_STELLA_RADIUS_UNIT}"} \
    ${PROXY_STELLA_TEMPERATURE_UNIT:+--proxy-stella-temperature-unit="${PROXY_STELLA_TEMPERATURE_UNIT}"} \
    ${PROXY_STELLA_INTERP_COORD:+--proxy-stella-interpolation-coordinate="${PROXY_STELLA_INTERP_COORD}"} \
    ${PROXY_V2_STRATEGY:+--proxy-v2-strategy="${PROXY_V2_STRATEGY}"} \
    ${PROXY_PREV_MODEL:+--proxy-prev-model="${PROXY_PREV_MODEL}"} \
    "${CMFGEN_DIR}"

# Parse job IDs from submit.log (in CMFGEN_DIR where submit_cmfgen.sh writes it)
FINAL_JOB=$(grep "CMFGEN Job ID" "${CMFGEN_DIR}/submit.log" | awk '{print $4}' | head -1)
OBS_JOB=$(grep "OBS Job ID" "${CMFGEN_DIR}/submit.log" | awk '{print $4}' | head -1)
# If FINAL was skipped, FINAL_JOB will be "(skipped" — set to empty
if [[ "${FINAL_JOB}" == "(skipped" || "${FINAL_JOB}" == "" ]]; then
  FINAL_JOB=""
fi
# If OBS was skipped, OBS_JOB will be "(skipped" — set to empty
if [[ "${OBS_JOB}" == "(skipped" || "${OBS_JOB}" == "" ]]; then
  OBS_JOB=""
fi

echo ""
echo "======================================================"
echo "CMFGEN jobs submitted successfully"
echo "======================================================"
echo "CMFGEN: ${FINAL_JOB}"
echo "OBS:    ${OBS_JOB}"
echo "======================================================"

# Output job IDs for dependency tracking (INIT removed; format: :FINAL:OBS)
echo ":${FINAL_JOB}:${OBS_JOB}"
