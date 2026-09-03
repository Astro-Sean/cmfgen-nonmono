#!/usr/bin/env bash
# submit_cmfgen.sh
#
# Creates Slurm scripts in the model directory:
#   - run_cmfgen_SINGLE_FINAL.slurm  (LTE cold start; default for first day or single snapshot)
#   - run_cmfgen_FINAL.slurm         (warm restart from *_IN files; for time-series subsequent days)
#   - run_cmfgen_OBS.slurm           (afterok FINAL)
#
# Also creates: <MODEL_DIR>/obs/ and copies templates from:
#   $CMFDIST/obs_template/{batobs.sh,CMF_FLUX_PARAM_INIT}
#
# Options:
#   --debug         Use interactive partition preset
#   --all           Submit full chain (same as default now; INIT/STABLE removed)
#   --single-final  Use SINGLE_FINAL slurm (LTE cold start; this is the default)
#   --quick         Request minimal resources (4 CPUs, 32G mem, 12h FINAL) for faster queue turnaround
#   --no-submit     Only write files; do not call sbatch
#   --restart       Use FINAL slurm (warm restart from *_IN files; for time-series)
#   --vel-law=N     Set CMFGEN VEL_LAW value (default: 12 for IIn tau-based grid; 11 for direct hydro)
#   -h|--help       Show help
#
# Notes:
# - Slurm dependencies use afterok.
# - In --restart mode, FINAL is submitted with no dependency, OBS depends on FINAL.
# - INIT script is still written for consistency, but it is not submitted in --restart mode.
# - INIT has no graceful timer / requeue. FINAL requeue policy is set by CMFGEN_FINAL_REQUEUE_WHEN (see below).
# - All submission output is logged to <MODEL_DIR>/submit.log
# - When SN_HYDRO_DATA exists, writes overview PNGs to <MODEL_DIR>/plots/ (structure + mass fractions).
# - Full runs (not --restart): runs CMFDIST/com/restart_model.sh --all before regenerating Slurm;
#   INIT jobs also run it at start so every fresh allocation matches a clean model tree.

set -euo pipefail

usage() { sed -n '1,23p' "$0"; }


rm -f -- cmf_INIT_*.out cmf_INIT_*.err cmf_STABLE_*.out cmf_STABLE_*.err cmf_FINAL_*.out cmf_FINAL_*.err cmf_OBS_*.out cmf_OBS_*.err

debug=0
all_jobs=0
single_final=1
quick=0
no_submit=0
restart=0
VEL_LAW="12"
TARGET_ND=""
MODEL_DIR=""
DEPENDENCY_JOB=""
MODEL_TYPE="supernova"
INC_RAD_DECAYS="1"
NO_RELATIVISTIC="0"
PROXY_POPULATIONS="0"
PROXY_POPULATIONS_MIN_ION_FRAC=""
PROXY_POPULATIONS_T_EXC_MIN=""
PROXY_POPULATIONS_NLTE_SOLVER="0"
PROXY_POPULATIONS_GREY_MAX_ITER=""
PROXY_POPULATIONS_GREY_TOL=""
PROXY_POPULATIONS_V2="0"
PROXY_COMPLEX_MODEL="0"
PROXY_POPULATIONS_MASS_FRAC_FLOOR=""
INIT_PARTITION="general"
FINAL_PARTITION=""
INIT_TIME="23:59:00"
INIT_NUM_ITS="${INIT_NUM_ITS:-1000}"
SLURM_MEM="16G"
INIT_SLURM_MEM="16G"
CPUS_OVERRIDE=""
N_OB_INS="30"
FINAL_REQUEUE_WHEN="wall"
FINAL_REQUEUE_MAX="500"
AUTO_SET_RSTAR="0"
AUTO_SET_LSTAR="0"
PHOT_TAU="0.6666667"
BB_TEFF_K="9100"
LSTAR_RADIUS_MODE="RMAX"
JOB_PREFIX=""
MODEL_NAME=""
DAY_LABEL=""
PREV_CMFGEN_DIR=""
PROJECT_ROOT=""
CONFIG_FILE=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --debug) debug=1; shift ;;
    --all) all_jobs=1; shift ;;
    --single-final) single_final=1; shift ;;
    --quick) quick=1; shift ;;
    --no-submit) no_submit=1; shift ;;
    --restart) restart=1; shift ;;
    --dependency=*) DEPENDENCY_JOB="${1#*=}"; shift ;;
    --vel-law=*) VEL_LAW="${1#*=}"; shift ;;
    --use-hydro-grid=*) USE_HYDRO_GRID="${1#*=}"; shift ;;
    --smooth-velocity-jumps=*) SMOOTH_VELOCITY_JUMPS="${1#*=}"; shift ;;
    --smooth-max-dv-ratio=*) SMOOTH_MAX_DV_RATIO="${1#*=}"; shift ;;
    --smooth-n-insert=*) SMOOTH_N_INSERT="${1#*=}"; shift ;;
    --smooth-rho-sigma=*) SMOOTH_RHO_SIGMA="${1#*=}"; shift ;;
    --smooth-temp-sigma=*) SMOOTH_TEMP_SIGMA="${1#*=}"; shift ;;
    --blend-composition) BLEND_COMPOSITION=1; shift ;;
    --blend-sigma=*) BLEND_SIGMA="${1#*=}"; shift ;;
    --blend-solar-weight=*) BLEND_SOLAR_WEIGHT="${1#*=}"; shift ;;
    --min-mass-fraction=*) MIN_MASS_FRACTION="${1#*=}"; shift ;;
    --monotonize-velocity) MONOTONIZE_VELOCITY=1; shift ;;
    --chemical-n-insert=*) CHEMICAL_N_INSERT="${1#*=}"; shift ;;
    --chemical-jump-threshold=*) CHEMICAL_JUMP_THRESHOLD="${1#*=}"; shift ;;
    --reversal-n-insert=*) REVERSAL_N_INSERT="${1#*=}"; shift ;;
    --target-nd=*) TARGET_ND="${1#*=}"; shift ;;
    --init-num-its=*) INIT_NUM_ITS="${1#*=}"; shift ;;
    --model-type=*) MODEL_TYPE="${1#*=}"; shift ;;
    --inc-rad-decays=*) INC_RAD_DECAYS="${1#*=}"; shift ;;
    --no-relativistic=*) NO_RELATIVISTIC="${1#*=}"; shift ;;
    --proxy-populations=*) PROXY_POPULATIONS="${1#*=}"; shift ;;
    --proxy-populations-min-ion-frac=*) PROXY_POPULATIONS_MIN_ION_FRAC="${1#*=}"; shift ;;
    --proxy-populations-t-exc-min=*) PROXY_POPULATIONS_T_EXC_MIN="${1#*=}"; shift ;;
    --proxy-populations-nlte-solver=*) PROXY_POPULATIONS_NLTE_SOLVER="${1#*=}"; shift ;;
    --proxy-populations-grey-max-iter=*) PROXY_POPULATIONS_GREY_MAX_ITER="${1#*=}"; shift ;;
    --proxy-populations-grey-tol=*) PROXY_POPULATIONS_GREY_TOL="${1#*=}"; shift ;;
    --proxy-populations-v2=*) PROXY_POPULATIONS_V2="${1#*=}"; shift ;;
    --proxy-complex-model=*) PROXY_COMPLEX_MODEL="${1#*=}"; shift ;;
    --proxy-populations-mass-frac-floor=*) PROXY_POPULATIONS_MASS_FRAC_FLOOR="${1#*=}"; shift ;;
    --partition=*) INIT_PARTITION="${1#*=}"; FINAL_PARTITION="${1#*=}"; shift ;;
    --time=*) INIT_TIME="${1#*=}"; TIME_DEFAULT="${1#*=}"; shift ;;
    --obs-time=*) OBS_TIME_OVERRIDE="${1#*=}"; shift ;;
    --mem=*) SLURM_MEM="${1#*=}"; shift ;;
    --init-mem=*) INIT_SLURM_MEM="${1#*=}"; shift ;;
    --cpus=*) CPUS_OVERRIDE="${1#*=}"; shift ;;
    --n-ob-ins=*) N_OB_INS="${1#*=}"; shift ;;
    --final-requeue-when=*) FINAL_REQUEUE_WHEN="${1#*=}"; shift ;;
    --final-requeue-max=*) FINAL_REQUEUE_MAX="${1#*=}"; shift ;;
    --auto-set-rstar=*) AUTO_SET_RSTAR="${1#*=}"; shift ;;
    --auto-set-lstar=*) AUTO_SET_LSTAR="${1#*=}"; shift ;;
    --phot-tau=*) PHOT_TAU="${1#*=}"; shift ;;
    --bb-teff-k=*) BB_TEFF_K="${1#*=}"; shift ;;
    --lstar-radius-mode=*) LSTAR_RADIUS_MODE="${1#*=}"; shift ;;
    --job-prefix=*) JOB_PREFIX="${1#*=}"; shift ;;
    --model-name=*) MODEL_NAME="${1#*=}"; shift ;;
    --day-label=*) DAY_LABEL="${1#*=}"; shift ;;
    --prev-cmfgen-dir=*) PREV_CMFGEN_DIR="${1#*=}"; shift ;;
    --project-root=*) PROJECT_ROOT="${1#*=}"; shift ;;
    --config-file=*) CONFIG_FILE="${1#*=}"; shift ;;
    -h|--help) usage; exit 0 ;;
    --) shift; break ;;
    -*)
      echo "ERROR: unknown option: $1" >&2
      usage >&2
      exit 1
      ;;
    *)
      if [[ -z "$MODEL_DIR" ]]; then
        MODEL_DIR="$1"
      else
        echo "ERROR: unexpected extra argument: $1" >&2
        usage >&2
        exit 1
      fi
      shift
      ;;
  esac
done

MODEL_DIR="${MODEL_DIR:-$PWD}"
MODEL_DIR="$(cd "$MODEL_DIR" && pwd)"

# ---------------------------
# CMFGEN environment paths (needed early for manual cleanup)
# ---------------------------
CMFDIST="${CMFDIST:-/viper/ptmp1/sbrennan/SN2024pba/cmfgen_source}"

# ---------------------------
# Manual cleanup (if restart_model.sh not found)
# Must happen BEFORE submit.log is created
# ---------------------------
RESTART_HELPER="${CMFDIST}/com/restart_model.sh"
if [[ $restart -eq 0 && $no_submit -eq 0 && ! -f "$RESTART_HELPER" ]]; then
  echo "WARNING: restart_model.sh not found at $RESTART_HELPER"
  echo "  Performing manual cleanup: keeping base files and SN_HYDRO_DATA, deleting everything else..."
  
  # Manual cleanup: keep base files and SN_HYDRO_DATA, delete everything else
  # Base files to keep (copied from /ptmp/sbrennan/SN2024pba/cmfgen_base)
  # IT_SPECIFIER, ADJUST_CORRECTIONS, ADJUST_R_DEFAULTS are convergence
  # control files read by CMFGEN every iteration (see docs/full_descr.txt).
  BASE_FILES_TO_KEEP="batch.sh VADAT_INIT IN_ITS MODEL_SPEC NUC_DECAY_DATA arnaud_rothenflug.dat IT_SPECIFIER ADJUST_CORRECTIONS ADJUST_R_DEFAULTS RDINR T_IN VADAT submit_cmfgen.sh tablec1.dat tablec2.dat tablec3.dat tablec4.dat tablec5.dat tablec6.dat tablec7.dat tablec8.dat ca2_col_data_filtered"
  # Also keep SN_HYDRO_DATA, obs directory, and plots/ (SN_HYDRO_DATA diagnostic
  # plots from convert_stella_swd_to_cmfgen.py — structure_overview.png, mass_fractions.png)
  # Also keep EDDFACTOR + EDDFACTOR_INFO when present (for USE_FIXED_J initialization)
  FILES_TO_KEEP="$BASE_FILES_TO_KEEP SN_HYDRO_DATA obs plots EDDFACTOR EDDFACTOR_INFO"
  
  cd "${MODEL_DIR}"
  
  # Create temporary list of files to delete
  for item in *; do
    keep=0
    for keep_pattern in $FILES_TO_KEEP; do
      if [[ "$item" == "$keep_pattern" ]] || [[ "$item" == "$keep_pattern".* ]]; then
        keep=1
        break
      fi
    done
    if [[ $keep -eq 0 ]]; then
      echo "  Removing: $item"
      rm -rf "$item" 2>/dev/null || true
    fi
  done
  
  # Clean obs directory but keep the directory itself
  if [[ -d "obs" ]]; then
    for item in obs/*; do
      if [[ -e "$item" ]]; then
        echo "  Removing: $item"
        rm -rf "$item" 2>/dev/null || true
      fi
    done
  fi
  
  echo "  Manual cleanup complete"
  cd - > /dev/null
fi

# ---------------------------
# Redirect all output to submit.log (and terminal via tee)
# (Named submit.log to avoid conflict with CMFGEN's batch.log)
# ---------------------------
LOGFILE="${MODEL_DIR}/submit.log"
exec > >(tee "$LOGFILE") 2>&1

# Use model directory name for SLURM job names
FOLDER_NAME="$(basename "$MODEL_DIR")"
# Use JOB_PREFIX from CLI arg, otherwise fall back to FOLDER_NAME
JOB_PREFIX="${JOB_PREFIX:-${FOLDER_NAME}}"

# Model name and epoch label for centralized output filenames
CMFGEN_MODEL_BASE="${MODEL_NAME:-${MODEL_NAME_FROM_DIR}}"
# Derive epoch from SN_HYDRO_DATA if available (e.g. "4.90200E+00" -> "4.9d")
EPOCH_LABEL=""
if [[ -f "${MODEL_DIR}/SN_HYDRO_DATA" ]]; then
  _sn_age=$(grep "Time(days) since explosion" "${MODEL_DIR}/SN_HYDRO_DATA" 2>/dev/null | awk '{print $NF}')
  if [[ -n "$_sn_age" ]]; then
    EPOCH_LABEL=$(python3 -c "print(f'{float(\"${_sn_age}\"):g}d')" 2>/dev/null || echo "")
  fi
fi
# Fallback to day label (e.g. "5.0days" -> "5d")
if [[ -z "$EPOCH_LABEL" && -n "${DAY_LABEL:-}" ]]; then
  EPOCH_LABEL=$(echo "${DAY_LABEL}" | sed 's/days/d/' | sed 's/day/d/')
fi
if [[ -z "$EPOCH_LABEL" ]]; then
  EPOCH_LABEL="${FOLDER_NAME}"
fi


echo "======================================================"
echo "CMFGEN Submission Script"
echo "Model Directory: $MODEL_DIR"
echo "Log File: $LOGFILE"
echo "Debug Mode: $debug"
echo "All Jobs General: $all_jobs"
echo "Quick (minimal resources): $quick"
echo "No-Submit: $no_submit"
echo "Restart: $restart"
echo "Single-Final (LTE cold start, skip INIT): $single_final"
echo "======================================================"

# ---------------------------
# CMFGEN environment paths
# ---------------------------
CMFDIST="${CMFDIST:-/viper/ptmp1/sbrennan/SN2024pba/cmfgen_source}"
ATOMIC="${HOME}/atomic"

PROJECT_ROOT="${PROJECT_ROOT:-$(dirname "$(dirname "$(dirname "$MODEL_DIR")")")}"
VADAT_SCRIPT="${PROJECT_ROOT}/scripts/prepare_VADAT.py"
INITS_SCRIPT="${PROJECT_ROOT}/scripts/prepare_IN_INTS.py"
PROXY_POPULATIONS_SCRIPT="${PROJECT_ROOT}/scripts/proxy_populations.py.deprecate"
# Model name from CLI arg; fallback to parent directory name
MODEL_NAME_FROM_DIR="${MODEL_NAME:-$(basename "$(dirname "$MODEL_DIR")")}"
# Use CONFIG_FILE from CLI arg if provided; otherwise find the config by name.
# Config files are named <MODEL_NAME>.yaml (not model_<MODEL_NAME>.yaml).
if [[ -n "${CONFIG_FILE:-}" && -f "${CONFIG_FILE}" ]]; then
  CONFIG_FILE="${CONFIG_FILE}"
elif [[ -f "${PROJECT_ROOT}/configs/${MODEL_NAME_FROM_DIR}.yaml" ]]; then
  CONFIG_FILE="${PROJECT_ROOT}/configs/${MODEL_NAME_FROM_DIR}.yaml"
else
  # Strip day label suffix (e.g. "X_2_5.0days" -> "X_2") for config lookup.
  # Day labels are appended to MODEL_NAME in some submission paths but config
  # files are named after the base model name only.
  # Use % (non-greedy) to strip only the day label, not parts of the model name.
  MODEL_NAME_BASE="${MODEL_NAME_FROM_DIR%_*days}"
  # Try: base name, base name without underscores (e.g. X_2 -> X2), legacy model_ prefix
  if [[ "${MODEL_NAME_BASE}" != "${MODEL_NAME_FROM_DIR}" && -f "${PROJECT_ROOT}/configs/${MODEL_NAME_BASE}.yaml" ]]; then
    CONFIG_FILE="${PROJECT_ROOT}/configs/${MODEL_NAME_BASE}.yaml"
  elif [[ -f "${PROJECT_ROOT}/configs/${MODEL_NAME_BASE//_}.yaml" ]]; then
    # Config naming convention may omit underscores (e.g. X_2 -> X2.yaml)
    CONFIG_FILE="${PROJECT_ROOT}/configs/${MODEL_NAME_BASE//_}.yaml"
  elif [[ -f "${PROJECT_ROOT}/configs/model_${MODEL_NAME_FROM_DIR}.yaml" ]]; then
    # Legacy fallback: model_<name>.yaml
    CONFIG_FILE="${PROJECT_ROOT}/configs/model_${MODEL_NAME_FROM_DIR}.yaml"
  else
    CONFIG_FILE="${PROJECT_ROOT}/configs/model_${MODEL_NAME_BASE}.yaml"
  fi
fi
# Hard error: CONFIG_FILE must exist — vadat_overrides are critical for convergence.
if [[ ! -f "${CONFIG_FILE}" ]]; then
  echo "ERROR: CONFIG_FILE not found: ${CONFIG_FILE}" >&2
  echo "  Looked for: ${PROJECT_ROOT}/configs/${MODEL_NAME_FROM_DIR}.yaml" >&2
  echo "  And base:   ${PROJECT_ROOT}/configs/${MODEL_NAME_BASE:-${MODEL_NAME_FROM_DIR}}.yaml" >&2
  echo "  And no-us:  ${PROJECT_ROOT}/configs/${MODEL_NAME_BASE:-${MODEL_NAME_FROM_DIR}}//_.yaml" >&2
  echo "  And legacy: ${PROJECT_ROOT}/configs/model_${MODEL_NAME_FROM_DIR}.yaml" >&2
  echo "  vadat_overrides (THK_CONT, OB_RAT, etc.) are required for correct CMFGEN setup." >&2
  echo "  Pass --config-file=/path/to/config.yaml explicitly." >&2
  exit 1
fi
echo "  CONFIG_FILE: ${CONFIG_FILE}"
NOTIFICATION_EMAIL=$(python3 - "${CONFIG_FILE}" <<'PY'
import sys
import yaml
with open(sys.argv[1]) as stream:
    config = yaml.safe_load(stream) or {}
print(config.get('cmfgen', {}).get('notification_email', ''))
PY
)

OBS_TEMPLATE_DIR="${CMFDIST}/obs_template"
if [[ ! -d "${OBS_TEMPLATE_DIR}" ]]; then
  OBS_TEMPLATE_DIR="${HOME}/cur_cmf_og/obs_template"
fi

echo ""
echo "Environment:"
echo "  CMFDIST: $CMFDIST"
echo "  ATOMIC:  $ATOMIC"
echo "  VADAT:   $VADAT_SCRIPT"
echo "  INITS:   $INITS_SCRIPT"
echo "  OBS_TPL: $OBS_TEMPLATE_DIR"
if [[ "${PROXY_POPULATIONS}" == "1" || "${PROXY_POPULATIONS_V2}" == "1" ]]; then
  echo "  PROXY:   $PROXY_POPULATIONS_SCRIPT (proxy dilute-LTE _IN files)"
fi

# ---------------------------
# Initialize optional flags for prepare_VADAT.py (set -u safe for heredoc use)
# ---------------------------
NO_REL_FLAG=""
if [[ "${NO_RELATIVISTIC}" == "1" ]]; then
  NO_REL_FLAG="--no-relativistic"
fi
MODEL_TYPE_FLAG="--model-type=${MODEL_TYPE}"
INC_RAD_DECAYS_FLAG="--inc-rad-decays=${INC_RAD_DECAYS}"
TARGET_ND_FLAG=""
if [[ -n "${TARGET_ND:-}" ]]; then
  TARGET_ND_FLAG="--target-nd=${TARGET_ND}"
fi
USE_HYDRO_GRID_FLAG=""
if [[ "${USE_HYDRO_GRID:-0}" == "1" ]]; then
  USE_HYDRO_GRID_FLAG="--use-hydro-grid"
fi
SMOOTH_JUMPS_FLAG=""
if [[ "${SMOOTH_VELOCITY_JUMPS:-0}" == "1" ]]; then
  SMOOTH_JUMPS_FLAG="--smooth-velocity-jumps --smooth-max-dv-ratio=${SMOOTH_MAX_DV_RATIO:-2.0} --smooth-n-insert=${SMOOTH_N_INSERT:-8}"
  if [[ -n "${SMOOTH_RHO_SIGMA:-}" && "${SMOOTH_RHO_SIGMA}" != "0" && "${SMOOTH_RHO_SIGMA}" != "0.0" ]]; then
    SMOOTH_JUMPS_FLAG="${SMOOTH_JUMPS_FLAG} --smooth-rho-sigma=${SMOOTH_RHO_SIGMA}"
  fi
  if [[ -n "${SMOOTH_TEMP_SIGMA:-}" && "${SMOOTH_TEMP_SIGMA}" != "0" && "${SMOOTH_TEMP_SIGMA}" != "0.0" ]]; then
    SMOOTH_JUMPS_FLAG="${SMOOTH_JUMPS_FLAG} --smooth-temp-sigma=${SMOOTH_TEMP_SIGMA}"
  fi
fi
BLEND_FLAG=""
if [[ "${BLEND_COMPOSITION:-0}" == "1" ]]; then
  BLEND_FLAG="--blend-composition --blend-sigma=${BLEND_SIGMA:-5.0} --blend-solar-weight=${BLEND_SOLAR_WEIGHT:-0.0}"
fi

# ---------------------------
# Resource presets (smaller requests = faster queue scheduling on most clusters)
# Increase MEM_DEFAULT / CPUS_DEFAULT if FINAL fails with OOM or needs more threads.
# ---------------------------
PARTITION_DEFAULT="general"
TIME_DEFAULT="${TIME_DEFAULT:-24:00:00}"
CPUS_DEFAULT=8
# CMFGEN BA / matrix steps often need many GB; too-small requests can OOM during the first BA build.
# Default 32G for production; debug uses 30G on interactive partition.
MEM_DEFAULT="${SLURM_MEM}"

# INIT can often run with less than FINAL, but still needs enough for the first BA build.
INIT_MEM_DEFAULT="${INIT_SLURM_MEM}"

PARTITION_DEBUG="interactive"
TIME_DEBUG="01:59:00"
CPUS_DEBUG=4
MEM_DEBUG="30G"

# Quick preset: minimal resources for faster backfill / shorter queue wait
TIME_QUICK="24:00:00"
CPUS_QUICK=8
MEM_QUICK="${SLURM_MEM}"

# INIT: default general partition, 24h. Override with CMFGEN_INIT_* env vars.
# Keep memory aligned with main MEM_DEFAULT to avoid early BA OOM.
INIT_PARTITION="${INIT_PARTITION}"
INIT_TIME="${INIT_TIME}"

# FINAL/SINGLE_FINAL: high-iteration mode; graceful chunking/requeue handles walltime.
FINAL_NUM_ITS=100000

OBS_TIME_DEFAULT="24:00:00"
OBS_TIME_DEBUG="01:30:00"

FINAL_SLURM="run_cmfgen_FINAL.slurm"
SINGLE_FINAL_SLURM="run_cmfgen_SINGLE_FINAL.slurm"
OBS_SLURM="run_cmfgen_OBS.slurm"

if [[ $debug -eq 1 ]]; then
  part="$PARTITION_DEBUG"
  timelimit="$TIME_DEBUG"
  cpus="$CPUS_DEBUG"
  mem="$MEM_DEBUG"
  init_mem="$MEM_DEBUG"
  obs_timelimit="$OBS_TIME_DEBUG"
  tag="debug"
elif [[ $quick -eq 1 ]]; then
  part="$PARTITION_DEFAULT"
  timelimit="$TIME_QUICK"
  cpus="$CPUS_QUICK"
  mem="$MEM_QUICK"
  init_mem="$MEM_QUICK"
  obs_timelimit="$TIME_QUICK"
  tag=""
else
  part="${FINAL_PARTITION:-$PARTITION_DEFAULT}"
  timelimit="$TIME_DEFAULT"
  cpus="$CPUS_DEFAULT"
  mem="$MEM_DEFAULT"
  init_mem="${INIT_MEM_DEFAULT:-$mem}"
  obs_timelimit="$OBS_TIME_DEFAULT"
  tag=""
fi

# Allow explicit CPU override from YAML config (via run_model_grid.py)
if [[ -n "${CPUS_OVERRIDE}" ]]; then
  cpus="$CPUS_OVERRIDE"
  echo "  CPU override from config: ${cpus} CPUs/task"
fi

# Allow explicit OBS time override (e.g. for interactive partition with <2h limit)
if [[ -n "${OBS_TIME_OVERRIDE:-}" ]]; then
  obs_timelimit="$OBS_TIME_OVERRIDE"
fi

# If using interactive partition (but not via --debug), cap OBS time to <2h
# since the interactive partition requires less than 2h runtime.
if [[ "$part" == "interactive" && "$obs_timelimit" == "24:00:00" ]]; then
  obs_timelimit="$OBS_TIME_DEBUG"
  echo "  Interactive partition detected: capping OBS time to $obs_timelimit"
fi

# FINAL graceful chunking (walltime-dependent):
# - Normal partition (24h): graceful stop at 20h (72000s), signal 4h before wall.
# - Debug mode (01:30): graceful stop disabled (FINAL_WALL_GRACEFUL=0).
# The sleep fallback fires AFTER the USR1 signal would fire, acting as a safety
# net in case SLURM misses the signal.
#
# User-configurable via GRACEFUL_STOP_SEC (overrides the automatic calculation).
GRACEFUL_STOP_SEC="${GRACEFUL_STOP_SEC:-}"

# Compute wall time in seconds for the current timelimit
_wall_sec_graceful=$(echo "$timelimit" | awk -F: '{if (NF==3) print $1*3600+$2*60+$3; else if (NF==2) print $1*60+$2; else print $1}')

if [[ -n "${GRACEFUL_STOP_SEC}" ]]; then
  FINAL_GRACEFUL_SLEEP_SEC="${GRACEFUL_STOP_SEC}"
else
  # Default: 83.3% of wall time for >=12h jobs, 75% for shorter (e.g. 2h interactive)
  if [[ "$_wall_sec_graceful" -ge 43200 ]]; then
    FINAL_GRACEFUL_SLEEP_SEC=$(echo "$_wall_sec_graceful" | awk '{printf "%d", $1 * 0.833}')
  else
    FINAL_GRACEFUL_SLEEP_SEC=$(echo "$_wall_sec_graceful" | awk '{printf "%d", $1 * 0.75}')
  fi
fi
FINAL_SIGNAL_BEFORE_SEC=$(( _wall_sec_graceful - FINAL_GRACEFUL_SLEEP_SEC ))
echo "  Graceful stop: ${FINAL_GRACEFUL_SLEEP_SEC}s into ${_wall_sec_graceful}s wall (signal ${FINAL_SIGNAL_BEFORE_SEC}s before wall)"

# FINAL continuation policy (written into run_cmfgen_FINAL.slurm):
#   wall         — default: Slurm SIGUSR1 before walltime + sleep fallback; NUM_ITS=0 then scontrol requeue.
#   finish       — no wall preemption; after each clean tcsh -f ./batch.sh exit, scontrol requeue until max restarts
#                  or you create .CMFGEN_FINAL_NO_REQUEUE in the model dir. Does not auto-submit OBS (see submit block).
#   recoverable  — no wall preemption; requeue only if OUTGEN matches a few soft-fail strings. OBS not auto-submitted.
FINAL_REQUEUE_WHEN="${FINAL_REQUEUE_WHEN}"
FINAL_WALL_GRACEFUL="${FINAL_WALL_GRACEFUL:-}"
if [[ -z "${FINAL_WALL_GRACEFUL}" ]]; then
  if [[ "${FINAL_REQUEUE_WHEN}" == "finish" || "${FINAL_REQUEUE_WHEN}" == "recoverable" ]]; then
    FINAL_WALL_GRACEFUL=0
  else
    FINAL_WALL_GRACEFUL=1
  fi
fi
FINAL_REQUEUE_MAX="${FINAL_REQUEUE_MAX}"

# Safety: if wall time < FINAL_SIGNAL_BEFORE_SEC, SLURM sends USR1 at job start (immediate requeue).
# This should not happen with the proportional calculation above, but keep as guard.
if [[ "${FINAL_WALL_GRACEFUL}" != "0" ]]; then
  _wall_sec=$(echo "$timelimit" | awk -F: '{if (NF==3) print $1*3600+$2*60+$3; else if (NF==2) print $1*60+$2; else print $1}')
  if [[ -n "$_wall_sec" ]] && [[ "$_wall_sec" -lt "$FINAL_SIGNAL_BEFORE_SEC" ]]; then
    FINAL_SIGNAL_BEFORE_SEC=$(( _wall_sec * 4 / 5 ))
    FINAL_GRACEFUL_SLEEP_SEC=$(( _wall_sec / 5 ))
    echo "  Adjusted FINAL_SIGNAL_BEFORE_SEC to ${FINAL_SIGNAL_BEFORE_SEC}s (wall time ${_wall_sec}s too short)"
  fi
fi

# Debug mode: limited iterations, no graceful stop, no requeue — UNLESS the
# config explicitly requests requeue (final_requeue_when != never and
# final_requeue_max > 0).  This lets debug/interactive runs use the same
# graceful-stop + requeue mechanism as production runs, which is essential
# for models that need more than one time-chunk to converge.
# Must come AFTER the graceful timer defaults so these take final effect.
if [[ $debug -eq 1 ]]; then
  FINAL_NUM_ITS=${INIT_NUM_ITS}
  if [[ "${FINAL_REQUEUE_WHEN}" == "never" || "${FINAL_REQUEUE_MAX}" == "0" ]]; then
    FINAL_GRACEFUL_SLEEP_SEC=999999999  # Disable sleep-based fallback timer
    FINAL_WALL_GRACEFUL=0               # Disable #SBATCH --signal + USR1 handler
  else
    echo "  Debug mode: graceful stop + requeue ENABLED (final_requeue_when=${FINAL_REQUEUE_WHEN}, max=${FINAL_REQUEUE_MAX})"
  fi
fi

echo ""
echo "Resources:"
echo "  Partition:       $part"
echo "  Time:            $timelimit"
echo "  CPUs/Task:       $cpus"
echo "  Memory (FINAL/OBS): $mem"
echo "  Partition:       $INIT_PARTITION"
echo "  Time:            $INIT_TIME"
echo "  FINAL NUM_ITS:   $FINAL_NUM_ITS"
echo "  OBS Time:        $obs_timelimit"
if [[ "${FINAL_WALL_GRACEFUL}" != "0" ]]; then
  echo "  FINAL wall-graceful: USR1 ${FINAL_SIGNAL_BEFORE_SEC}s before wall; timer fallback ${FINAL_GRACEFUL_SLEEP_SEC}s"
else
  echo "  FINAL wall-graceful: OFF (no SIGUSR1 / no in-job sleep timer)"
fi
echo "  FINAL requeue when: ${FINAL_REQUEUE_WHEN} (max restarts ${FINAL_REQUEUE_MAX}; stop file .CMFGEN_FINAL_NO_REQUEUE for finish mode)"
echo "======================================================"

# ---------------------------
# Sanity checks
# ---------------------------
die() { echo "ERROR: $*" >&2; exit 1; }

[[ -d "$MODEL_DIR" ]] || die "Model directory not found: $MODEL_DIR"
# SN_HYDRO_DATA may not exist yet in the time-series workflow (created by the
# post-STELLA conversion job).  Warn at submission time; the generated slurm
# script will check at runtime and exit gracefully if still missing.
if [[ ! -f "${MODEL_DIR}/SN_HYDRO_DATA" ]]; then
  echo "WARNING: SN_HYDRO_DATA not found in ${MODEL_DIR} at submission time."
  echo "  In the time-series workflow, it will be created by the conversion job."
  echo "  The slurm script will check again at runtime and skip if still missing."
fi
[[ -f "$VADAT_SCRIPT" ]] || die "Missing: $VADAT_SCRIPT"
[[ -f "$INITS_SCRIPT" ]] || die "Missing: $INITS_SCRIPT"
[[ -x "${MODEL_DIR}/batch.sh" ]] || die "Missing or not executable: ${MODEL_DIR}/batch.sh"
[[ -x "${CMFDIST}/com/out2in.sh" ]] || die "Missing or not executable: ${CMFDIST}/com/out2in.sh"

[[ -d "$OBS_TEMPLATE_DIR" ]] || die "Missing obs template dir: $OBS_TEMPLATE_DIR"
[[ -f "${OBS_TEMPLATE_DIR}/batobs.sh" ]] || die "Missing: ${OBS_TEMPLATE_DIR}/batobs.sh"
[[ -f "${OBS_TEMPLATE_DIR}/CMF_FLUX_PARAM_INIT" ]] || die "Missing: ${OBS_TEMPLATE_DIR}/CMF_FLUX_PARAM_INIT"

if [[ $no_submit -eq 0 ]]; then
  command -v sbatch >/dev/null 2>&1 || die "sbatch not found in PATH (use --no-submit to just generate files)"
fi

# Pre-create obs/ so SLURM --chdir won't fail at job start
mkdir -p "${MODEL_DIR}/obs"

# SN_HYDRO_DATA overview: structure + mass fractions -> plots/ (non-fatal if matplotlib or file missing)
PLOT_SN_HYDRO_SCRIPT="${PROJECT_ROOT}/scripts/plot_sn_hydro_cmfgen_inputs.py"
if [[ ! -f "$PLOT_SN_HYDRO_SCRIPT" ]]; then
  # Fallback: try to find in VADAT_SCRIPT directory (legacy)
  PLOT_SN_HYDRO_SCRIPT="$(dirname "$VADAT_SCRIPT")/plot_sn_hydro_cmfgen_inputs.py"
fi
HYDRO_DATA="${MODEL_DIR}/SN_HYDRO_DATA"
if [[ -f "$HYDRO_DATA" ]]; then
  if [[ -f "$PLOT_SN_HYDRO_SCRIPT" ]]; then
    mkdir -p "${MODEL_DIR}/plots"
    export PREPARE_VADAT_PY="$VADAT_SCRIPT"
    echo "Writing SN_HYDRO_DATA overview plots -> ${MODEL_DIR}/plots"
    if command -v python3 >/dev/null 2>&1; then
      python3 "$PLOT_SN_HYDRO_SCRIPT" "$MODEL_DIR" || echo "WARNING: hydro overview plots failed (non-fatal)"
    else
      python "$PLOT_SN_HYDRO_SCRIPT" "$MODEL_DIR" || echo "WARNING: hydro overview plots failed (non-fatal)"
    fi
  else
    echo "WARNING: plot script not found (skip hydro plots): $PLOT_SN_HYDRO_SCRIPT"
  fi
else
  echo "WARNING: no SN_HYDRO_DATA in model dir (skip hydro overview plots): $HYDRO_DATA"
fi

# Full (non-restart) pipeline: hard-reset stale CMFGEN run products so INIT matches a fresh model.
# restart_model.sh is csh; keeps VADAT / SN_HYDRO_DATA / batch.sh; see CMFDIST/com/restart_model.sh.
#
# IMPORTANT: Do not hard-reset when --no-submit is used. Users (and the time-series setup tooling)
# often call --no-submit to regenerate Slurm files while preserving an existing converged model
# for later --restart runs. INIT itself does its own hard reset inside the job script.
if [[ $restart -eq 0 && $no_submit -eq 0 ]]; then
  RESTART_HELPER="${CMFDIST}/com/restart_model.sh"
  if [[ -f "$RESTART_HELPER" ]]; then
    command -v csh >/dev/null 2>&1 || die "csh not in PATH (required for restart_model.sh)"
    echo "Running restart_model.sh --all in ${MODEL_DIR}"
    ( cd "${MODEL_DIR}" && csh "$RESTART_HELPER" --all ) || die "restart_model.sh --all failed"
  fi
  # Note: Manual cleanup for missing restart_model.sh is now done BEFORE submit.log creation (see lines 91-135)

  # Full resubmission safety: wipe all local "prep" markers and requeue sentinels so
  # no stale state can skip prep steps or leave NUM_ITS=0 from a prior graceful-stop.
  rm -f "${MODEL_DIR}/.CMFGEN_INIT_PREPPED" "${MODEL_DIR}/.CMFGEN_STABLE_PREPPED" "${MODEL_DIR}/.CMFGEN_FINAL_PREPPED" 2>/dev/null || true
  rm -f "${MODEL_DIR}/.CMFGEN_SINGLE_FINAL_PREPPED" 2>/dev/null || true
  rm -f "${MODEL_DIR}/.CMFGEN_REQUEST_REQUEUE" "${MODEL_DIR}/.CMFGEN_STABLE_REQUEST_REQUEUE" 2>/dev/null || true
  rm -f "${MODEL_DIR}/.CMFGEN_FINAL_NEEDS_NUM_ITS_RESET" 2>/dev/null || true
fi

# If we're doing a full (re)run, force prep to re-run.
# In --restart mode, we preserve the existing prep markers.
if [[ $restart -eq 0 ]]; then
  rm -f "${MODEL_DIR}/.CMFGEN_INIT_PREPPED" "${MODEL_DIR}/.CMFGEN_STABLE_PREPPED" "${MODEL_DIR}/.CMFGEN_FINAL_PREPPED" 2>/dev/null || true
  rm -f "${MODEL_DIR}/.CMFGEN_SINGLE_FINAL_PREPPED" 2>/dev/null || true
fi

# Always rewrite all scripts
rm -f "${MODEL_DIR}/"{${FINAL_SLURM},${SINGLE_FINAL_SLURM},${OBS_SLURM}} 2>/dev/null || true

# ---------------------------------------------------------
# Write FINAL slurm (restore+delete)
# ---------------------------------------------------------
# Define optional flags for prepare_VADAT.py (set -u safe for heredoc use)
MONOTONIZE_FLAG=""
if [[ "${MONOTONIZE_VELOCITY:-0}" == "1" ]]; then
  MONOTONIZE_FLAG="--monotonize-velocity"
fi
CHEMICAL_FLAG=""
if [[ "${CHEMICAL_N_INSERT:-0}" != "0" ]]; then
  CHEMICAL_FLAG="--chemical-n-insert=${CHEMICAL_N_INSERT} --chemical-jump-threshold=${CHEMICAL_JUMP_THRESHOLD:-0.1}"
fi
REVERSAL_FLAG=""
if [[ "${REVERSAL_N_INSERT:-0}" != "0" ]]; then
  REVERSAL_FLAG="--reversal-n-insert=${REVERSAL_N_INSERT}"
fi

cat > "${MODEL_DIR}/${FINAL_SLURM}" <<EOF
#!/bin/bash
#SBATCH --job-name=${JOB_PREFIX}_cmfgen${tag:+_${tag}}
#SBATCH --output=cmf_%j.out
#SBATCH --error=cmf_%j.err
#SBATCH --time=${timelimit}
#SBATCH --partition=${part}
$([ "${part}" != "interactive" ] && echo "#SBATCH --qos=n0004")

$([ "${FINAL_WALL_GRACEFUL}" != "0" ] && echo "#SBATCH --signal=B:USR1@${FINAL_SIGNAL_BEFORE_SEC}")
#SBATCH --requeue
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${cpus}
#SBATCH --mem=${mem}
# #SBATCH --mail-type=END
# #SBATCH --mail-user=sbrennan@mpe.mpg.de
#SBATCH --chdir=${MODEL_DIR}

set -euo pipefail

# Runtime check: SN_HYDRO_DATA must exist (created by post-STELLA conversion job)
if [[ ! -f "${MODEL_DIR}/SN_HYDRO_DATA" ]]; then
  echo "ERROR: SN_HYDRO_DATA not found in ${MODEL_DIR} at runtime." >&2
  echo "  This snapshot may be beyond the epoch STELLA reached before crashing." >&2
  exit 1
fi

check_outgen_for_fatal() {
  local outgen="OUTGEN"
  if [[ ! -f "\$outgen" ]]; then
    echo "ERROR: missing \$outgen after batch.sh (CMFGEN did not write OUTGEN?)" >&2
    return 2
  fi
  if grep -v "recommend --no-relativistic to prevent" "\$outgen" | grep -Eq "NaN|Error in MON_INTERP|Segmentation fault|SIG(ABRT|SEGV)|floating invalid|Error in comp_j_blank\\.f|Mean intensity blowing up|bad initial population guesses|Predicted changes are too large|Error in SETPOP_WLD_V2|zero RGU"; then
    echo "ERROR: fatal pattern detected in \$outgen (NaN / MON_INTERP / crash). Aborting stage." >&2
    echo "---- OUTGEN tail (last 80 lines) ----" >&2
    tail -n 80 "\$outgen" >&2 || true
    echo "------------------------------------" >&2
    return 3
  fi
}

# Track convergence across requeues in a small file. CMFGEN appends to OUTGEN
# with APPEND mode, so OUTGEN grows unboundedly across requeues (observed
# 433 GB after 3 requeues). We rotate OUTGEN each requeue and track the max
# population change in CONVERGENCE_HISTORY to detect divergence.
CONV_HIST="CONVERGENCE_HISTORY"
MAX_REQUEUE=${FINAL_REQUEUE_MAX:-500}

# Extract the last "Maximum % increase" value from OUTGEN using tail (fast
# even on multi-hundred-GB files because tail seeks to end-of-file).
track_convergence() {
  local outgen="OUTGEN"
  [[ -f "\$outgen" ]] || return 0
  local val
  val=\$(tail -n 200 "\$outgen" 2>/dev/null | grep -m1 "Maximum % increase at depth" | sed 's/.*is *\\([0-9.E+]*\\).*/\\1/' 2>/dev/null || true)
  if [[ -n "\$val" ]]; then
    local iter
    iter=\$(tail -n 200 "\$outgen" 2>/dev/null | grep -m1 "Current great iteration count" | sed 's/.*count is *\\([0-9]*\\).*/\\1/' 2>/dev/null || echo "?")
    echo "\${iter} \${val}" >> "\${CONV_HIST}"
    echo "Convergence tracking: iteration \${iter}, max % increase = \${val}%"
  fi
}

# Check CONVERGENCE_HISTORY for divergence. If the last N entries all show
# max change > DIVERGENCE_THRESHOLD and are not decreasing, the model is
# diverging and requeuing will not help.
check_convergence_for_divergence() {
  [[ -f "\${CONV_HIST}" ]] || return 0
  local threshold=50.0  # % change above which we consider the model diverging
  local n_required=5    # need this many consecutive high values to flag divergence
  local n_high=0
  local last_vals=""
  while read -r iter val; do
    [[ -z "\$val" ]] && continue
    last_vals="\${last_vals} \${val}"
    if (( \$(echo "\$val > \$threshold" | bc -l 2>/dev/null || echo 0) )); then
      n_high=\$((n_high + 1))
    else
      n_high=0  # reset on any non-high value
    fi
  done < "\${CONV_HIST}"
  if [[ \$n_high -ge \$n_required ]]; then
    echo "ERROR: CMFGEN model is DIVERGING — last max % changes:\${last_vals}" >&2
    echo "  \${n_high} consecutive iterations with max change > \${threshold}%." >&2
    echo "  Requeuing will not help. Stopping to prevent resource waste." >&2
    echo "  Possible causes:" >&2
    echo "    - Atomic model too complex (too many levels in MODEL_SPEC)" >&2
    echo "    - Bad initial populations (check proxy_populations)" >&2
    echo "    - Non-convergent Ng acceleration parameters (DO_NG, BEG_NG, IBEG_NG)" >&2
    echo "    - BA solver falling back to LAMBDA iteration (check OUTGEN)" >&2
    echo "  Fix: simplify MODEL_SPEC (reduce levels), check VADAT convergence params." >&2
    return 4
  fi
  return 0
}

# Rotate OUTGEN to prevent unbounded growth across requeues.
# CMFGEN opens OUTGEN with APPEND mode (cmfgen.f:99), so without rotation
# the file grows by ~150 GB per 24h requeue cycle.
rotate_outgen() {
  local outgen="OUTGEN"
  [[ -f "\$outgen" ]] || return 0
  local size
  size=\$(stat -c%s "\$outgen" 2>/dev/null || echo 0)
  # Rotate if larger than 2 GB
  if [[ \$size -gt 2147483648 ]]; then
    echo "Rotating OUTGEN (size: \${size} bytes) — saving tail to OUTGEN.tail"
    tail -n 5000 "\$outgen" > OUTGEN.tail 2>/dev/null || true
    : > "\$outgen"  # Truncate (CMFGEN will append to the empty file)
  fi
}

auto_set_rstar_lstar_from_hydro() {
  python3 - <<'PY'
from __future__ import annotations

import math
import os
from pathlib import Path

hydro = Path("SN_HYDRO_DATA")
vadat = Path("VADAT")
if not hydro.is_file():
    raise SystemExit(0)

auto_rstar = "${AUTO_SET_RSTAR}".strip() not in {"0", "F", "FALSE", "false"}
auto_lstar = "${AUTO_SET_LSTAR}".strip() not in {"0", "F", "FALSE", "false"}
tau_phot = float("${PHOT_TAU}")
teff = float("${BB_TEFF_K}")
lstar_rmode = ("${LSTAR_RADIUS_MODE}" or "RMAX").strip().upper()
if tau_phot <= 0:
    raise SystemExit("auto-RSTAR/LSTAR: PHOT_TAU must be > 0")
if teff <= 0:
    raise SystemExit("auto-RSTAR/LSTAR: BB_TEFF_K must be > 0")
if lstar_rmode not in {"RMAX", "PHOT"}:
    raise SystemExit(f"auto-RSTAR/LSTAR: invalid LSTAR_RADIUS_MODE={lstar_rmode!r} (use RMAX or PHOT)")

lines = hydro.read_text().splitlines()

def _read_block(marker: str) -> list[float]:
    try:
        i0 = lines.index(marker)
    except ValueError:
        return []
    out: list[float] = []
    for ln in lines[i0 + 1 :]:
        if not ln.strip():
            break
        for tok in ln.split():
            try:
                out.append(float(tok))
            except ValueError:
                pass
    return out

R_1e10 = _read_block("Radius grid (10^10cm)")
rho = _read_block("Density (gm/cm^3)")
kappa = _read_block("Kappa (cm^2/gm)")
if not (R_1e10 and rho and kappa) or not (len(R_1e10) == len(rho) == len(kappa)):
    print("auto-RSTAR/LSTAR: missing hydro blocks; skip")
    raise SystemExit(0)
n = len(R_1e10)

tau = [0.0] * n
for i in range(1, n):
    dr_cm = (R_1e10[i - 1] - R_1e10[i]) * 1.0e10
    if dr_cm <= 0:
        continue
    chi0 = kappa[i - 1] * rho[i - 1]
    chi1 = kappa[i] * rho[i]
    tau[i] = tau[i - 1] + 0.5 * (chi0 + chi1) * dr_cm

rphot_1e10 = R_1e10[-1]
if tau[-1] > tau_phot:
    j = 0
    while j < n and tau[j] < tau_phot:
        j += 1
    j = min(max(j, 1), n - 1)
    t0, t1 = tau[j - 1], tau[j]
    r0, r1 = R_1e10[j - 1], R_1e10[j]
    if t1 > t0:
        f = (tau_phot - t0) / (t1 - t0)
        rphot_1e10 = r0 + f * (r1 - r0)

sigma = 5.670374419e-5
Lsun = 3.828e33
rmax_1e10 = max(R_1e10)
r_for_lstar_1e10 = (rmax_1e10 if lstar_rmode == "RMAX" else rphot_1e10)
R_cm = max(r_for_lstar_1e10, 1.0e-99) * 1.0e10
Lbb_lsun = max((4.0 * math.pi * R_cm * R_cm * sigma * teff**4) / Lsun, 1.0e-30)

def _patch_vadat(tag: str, value: float) -> None:
    if not vadat.is_file():
        return
    vlines = vadat.read_text().splitlines()
    for i, ln in enumerate(vlines):
        if f"[{tag}]" in ln:
            vlines[i] = f"{value:.8g}  [{tag}]"
            vadat.write_text("\\n".join(vlines) + "\\n")
            print(f"auto-{tag}: set VADAT {tag}={value}")
            return

def _patch_header_lstar() -> None:
    for i, ln in enumerate(lines[:40]):
        if ln.strip().startswith("Energy Lum [lsun]:"):
            # Always overwrite (unless CMFGEN_AUTO_SET_LSTAR=0).
            lines[i] = f"Energy Lum [lsun]:                   {Lbb_lsun:.6E}"
            hydro.write_text("\\n".join(lines) + "\\n")
            print(
                "auto-LSTAR: wrote SN_HYDRO_DATA Energy Lum [lsun]={:.6e} (Teff={:.1f}K, mode={}, R_for_LSTAR={:.6e} 1e10cm)".format(
                    Lbb_lsun, teff, lstar_rmode, r_for_lstar_1e10
                )
            )
            return

if auto_lstar:
    _patch_header_lstar()
    _patch_vadat("LSTAR", Lbb_lsun)
if auto_rstar:
    _patch_vadat("RSTAR", rphot_1e10)
    # Keep outer radius fixed in physical units when RSTAR changes.
    rmax_ratio = max(rmax_1e10 / max(rphot_1e10, 1.0e-99), 1.0)
    _patch_vadat("RMAX", rmax_ratio)

PY
}

require_restart_inputs_or_abort() {
  # If CMFGEN cannot read POINT1/POINT2 it "starts a new model" and requires species *_IN files
  # (e.g. HI_IN) for population estimates. In --restart workflows this typically means the
  # directory is not actually restartable (e.g. was cleaned by restart_model.sh or never ran).
  local have_point=0
  local have_in=0
  if [[ -f POINT1 || -f POINT2 ]]; then
    have_point=1
  fi
  # Require at least one species input file (HI_IN is the first one CMFGEN asks for).
  # CFDAT__IN/T_IN/GAMMAS_IN alone are not sufficient for a new-model start.
  if ls *_IN >/dev/null 2>&1; then
    if [[ -f HI_IN ]]; then
      have_in=1
    fi
  fi
  if [[ \$have_point -eq 0 && \$have_in -eq 0 ]]; then
    echo "ERROR: directory is not restartable: missing POINT1/POINT2 and missing HI_IN (species *_IN files)." >&2
    echo "  Fix: run a full chain first (submit_cmfgen.sh --all), or restore restart state files from a previous run." >&2
    echo "  Hint: running submit_cmfgen.sh with --no-submit should not clean state anymore (this script was updated)." >&2
    exit 2
  fi
}

echo "======================================"
echo "CMFGEN FINAL Job Started"
echo "Job ID: \$SLURM_JOB_ID"
echo "Partition: \$SLURM_JOB_PARTITION"
echo "CPUs: \$SLURM_CPUS_PER_TASK"
echo "Node: \$SLURMD_NODENAME"
echo "======================================"

# Label this FINAL (re)queue instance as N1, N2, ... for easier tracking.
# Slurm provides SLURM_RESTART_COUNT during requeues/restarts (initial run typically 0).
restart_count="\${SLURM_RESTART_COUNT:-0}"
requeue_idx=\$((restart_count + 1))
final_job_name="${JOB_PREFIX}_cmfgen_N\${requeue_idx}"
if command -v scontrol >/dev/null 2>&1; then
  scontrol update JobId="\$SLURM_JOB_ID" JobName="\$final_job_name" >/dev/null 2>&1 || true
fi
echo "CMFGEN iteration label (for job-name tracking): \${final_job_name}"

# Do not delete Slurm stdout/stderr logs; they are essential for debugging failures.
module purge
module load intel/2024.0
module load mkl/2026.0
module load python-waterboa/2024.06
# Explicitly set MKL library path (workaround for module loading issues on compute nodes)
export LD_LIBRARY_PATH=/viper/u1/system/soft/RHEL_9/packages/x86_64/intel_oneapi/2026.0/mkl/2026.0/lib/intel64:\${LD_LIBRARY_PATH:+:\$LD_LIBRARY_PATH}

export cmfdist="${CMFDIST}"
export ATOMIC="${ATOMIC}"
export CMFGEN_MODEL_TYPE="${MODEL_TYPE}"
export SCRATCH_DIR="SCRATCH/"
mkdir -p SCRATCH/

# Sentinel indicating we previously requested a clean stop (NUM_ITS=0)
# and need to restore NUM_ITS back to 500 before starting CMFGEN again.
NUM_ITS_RESET_SENTINEL=".CMFGEN_FINAL_NEEDS_NUM_ITS_RESET"

# Graceful chunking (FINAL only): ~18h into a 24h walltime, stop cleanly (NUM_ITS=0),
# then out2in.sh and scontrol requeue so the model can continue until convergence.
request_graceful_stop() {
  echo ""
  echo "======================================================"
  echo "Pre-timeout signal received: requesting graceful stop."
  echo "Setting NUM_ITS=0 in IN_ITS so CMFGEN exits cleanly."
  echo "======================================================"
  python - <<'PY'
from pathlib import Path
p = Path("IN_ITS")
lines = p.read_text().splitlines() if p.exists() else []
out = []
done = False
for ln in lines:
    if "[NUM_ITS]" in ln:
        out.append("0 [NUM_ITS] ")
        done = True
    else:
        out.append(ln)
if not done:
    out.insert(0, "0 [NUM_ITS] ")
p.write_text("\n".join(out) + "\n")
PY
  touch .CMFGEN_REQUEST_REQUEUE
  touch "\${NUM_ITS_RESET_SENTINEL}"
}
trap request_graceful_stop USR1

# Fallback for clusters that do not reliably deliver pre-timeout signals:
# request graceful stop ~18h after start (for a 24h walltime).
graceful_timer_pid=""
# stop_graceful_timer: always defined; no-op when timer is disabled.
stop_graceful_timer() {
  if [[ -n "\${graceful_timer_pid}" ]] && kill -0 "\${graceful_timer_pid}" 2>/dev/null; then
    kill "\${graceful_timer_pid}" 2>/dev/null || true
  fi
}
# Start fallback sleep timer unless disabled (debug mode sets to 999999999).
if [[ ${FINAL_GRACEFUL_SLEEP_SEC} -lt 100000000 ]]; then
  (
    sleep ${FINAL_GRACEFUL_SLEEP_SEC}
    request_graceful_stop
  ) &
  graceful_timer_pid=\$!
  echo "Graceful timer started: fires in ${FINAL_GRACEFUL_SLEEP_SEC}s (PID \${graceful_timer_pid})"
else
  echo "Graceful timer disabled (FINAL_GRACEFUL_SLEEP_SEC=${FINAL_GRACEFUL_SLEEP_SEC})"
fi

# ---------------------------
# Resource settings (from SLURM)
# ---------------------------
CPUS=\${SLURM_CPUS_PER_TASK:-8}
export OMP_NUM_THREADS=\$CPUS
export MKL_NUM_THREADS=\$CPUS
export OMP_STACKSIZE=512M
export KMP_STACKSIZE=512M

echo "CPUs: \$CPUS"

# Check for divergence before continuing on requeue. If the model is
# diverging, stop instead of wasting more compute on a non-convergent cycle.
if [[ "\${SLURM_RESTART_COUNT:-0}" -gt 0 ]]; then
  check_convergence_for_divergence
  div_rc=\$?
  if [[ \$div_rc -ne 0 ]]; then
    echo "ERROR: Divergence detected on requeue \${requeue_idx}. Aborting." >&2
    exit \$div_rc
  fi
  # Check max requeue count to prevent infinite requeueing.
  if [[ "\${SLURM_RESTART_COUNT:-0}" -ge "\${MAX_REQUEUE}" ]]; then
    echo "ERROR: Maximum requeue count (\${MAX_REQUEUE}) reached. Stopping." >&2
    echo "  The model has not converged after \${SLURM_RESTART_COUNT} requeues." >&2
    echo "  Check CONVERGENCE_HISTORY and OUTGEN.tail for diagnostics." >&2
    exit 5
  fi
  # Rotate OUTGEN to prevent unbounded growth across requeues.
  rotate_outgen
fi

echo "Running out2in.sh..."
"\$cmfdist/com/out2in.sh"

init="INIT_SN_HYDRO_DATA"
cur="SN_HYDRO_DATA"
next="SN_HYDRO_FOR_NEXT_MODEL"
tmp="\${cur}.tmp.\$\$"

if [[ -f "\$init" ]]; then
  cp -p -- "\$init" "\$tmp"
  mv -- "\$tmp" "\$cur"
  echo "Restored \$init -> \$cur"
else
  echo "No \$init found, skipping restore"
fi

rm -f -- "\$next"

echo "Running final Python scripts..."
# Print flag status for logging
if [[ "${NO_RELATIVISTIC}" == "1" ]]; then
  echo "NO_RELATIVISTIC=1 (passed --no-relativistic to prepare_VADAT.py)"
fi
python "${VADAT_SCRIPT}" --final --convergence-mode stable --vel-law "${VEL_LAW}" ${MODEL_TYPE_FLAG} ${INC_RAD_DECAYS_FLAG} ${TARGET_ND_FLAG} ${USE_HYDRO_GRID_FLAG} ${SMOOTH_JUMPS_FLAG} ${BLEND_FLAG} ${MONOTONIZE_FLAG:-} ${CHEMICAL_FLAG:-} ${REVERSAL_FLAG:-} ${NO_REL_FLAG} --vadat-overrides "${CONFIG_FILE}" > 'batch.log'
python "${INITS_SCRIPT}" --final --num-its ${FINAL_NUM_ITS} >> 'batch.log'
export CMFGEN_AUTO_SET_RSTAR="${AUTO_SET_RSTAR}"
export CMFGEN_AUTO_SET_LSTAR="${AUTO_SET_LSTAR}"
auto_set_rstar_lstar_from_hydro

# Preserve restart state unconditionally.
#
# Deleting POINT1/POINT2/SCRTEMP/BAMAT* causes CMFGEN to "Starting a new model",
# resets the great-iteration counter, and can break --restart workflows (e.g. HI_IN missing).
#
# INIT already does a hard reset via restart_model.sh --all, so keeping restart files
# here is safe and makes FINAL reliably continue from the most recent state.
echo "FINAL: preserving restart state (POINT1/POINT2, SCRTEMP, BAMAT*, etc.)"

# If we previously requested a clean stop, restore NUM_ITS to FINAL_NUM_ITS now
# so the requeued FINAL run continues instead of immediately exiting.
if [[ -f "\${NUM_ITS_RESET_SENTINEL}" ]]; then
  echo "Resetting IN_ITS NUM_ITS to ${FINAL_NUM_ITS} for requeued FINAL run."
  python - <<PY
from pathlib import Path
p = Path("IN_ITS")
lines = p.read_text().splitlines() if p.exists() else []
out = []
done = False
for ln in lines:
    if "[NUM_ITS]" in ln:
        out.append("${FINAL_NUM_ITS} [NUM_ITS] ")
        done = True
    else:
        out.append(ln)
if not done:
    out.insert(0, "${FINAL_NUM_ITS} [NUM_ITS] ")
p.write_text("\\n".join(out) + "\\n")
PY
  rm -f "\${NUM_ITS_RESET_SENTINEL}"
fi

# If IN_ITS still has NUM_ITS=0 (e.g. no sentinel), bump to FINAL_NUM_ITS so FINAL does not exit immediately.
python - <<PY
from pathlib import Path

p = Path("IN_ITS")
lines = p.read_text().splitlines() if p.is_file() else []
num = None
for ln in lines:
    if "[NUM_ITS]" in ln:
        parts = ln.split()
        if parts:
            try:
                num = int(parts[0])
            except ValueError:
                num = None
        break
if num != 0:
    raise SystemExit(0)
out = []
done = False
for ln in lines:
    if "[NUM_ITS]" in ln:
        out.append("${FINAL_NUM_ITS} [NUM_ITS] ")
        done = True
    else:
        out.append(ln)
if not done:
    out.insert(0, "${FINAL_NUM_ITS} [NUM_ITS] ")
p.write_text("\\n".join(out) + "\\n")
print("IN_ITS had NUM_ITS=0; set to ${FINAL_NUM_ITS} before FINAL batch.sh.")
PY

echo "Running batch.sh..."
# Skip restart-state check for single-final cold starts (LTE_EST=T; no species _IN needed)
if [ "${single_final:-0}" != "1" ]; then
  require_restart_inputs_or_abort
fi
# Setup SCRATCH/ symlinks for files that don't use SCRATCH_DIR
# For FINAL (restart): move existing POINT1/POINT2/SCRTEMP to SCRATCH/ to preserve restart state
mkdir -p SCRATCH
for _f in POINT1 POINT2 SCRTEMP STEQ_VALS J_COMP JEW EWDATA CFDAT_OUT CONT_FREQ LINEHEAT NETRATE TOTRATE TRANS_INFO ADIABAT_CHK; do
  [[ -e "\$_f" && ! -L "\$_f" ]] && mv "\$_f" SCRATCH/ 2>/dev/null || true
  ln -sf "SCRATCH/\$_f" "\$_f" 2>/dev/null || true
done
tcsh -f ./batch.sh
check_outgen_for_fatal
stop_graceful_timer

# Track convergence: extract last max % change from OUTGEN and append to
# CONVERGENCE_HISTORY. This small file persists across requeues and is
# used to detect divergence before the next requeue.
track_convergence

echo "Preparing restart inputs (out2in.sh)..."
"\$cmfdist/com/out2in.sh"

if [[ -f .CMFGEN_REQUEST_REQUEUE ]]; then
  # Check for divergence before requeuing. If the model is diverging,
  # stop instead of wasting more compute on a non-convergent cycle.
  check_convergence_for_divergence
  div_rc=\$?
  if [[ \$div_rc -ne 0 ]]; then
    echo "ERROR: Divergence detected after batch.sh. NOT requeuing." >&2
    exit \$div_rc
  fi
  echo "Requeue requested; requeueing job \$SLURM_JOB_ID..."
  rm -f .CMFGEN_REQUEST_REQUEUE
  scontrol requeue "\$SLURM_JOB_ID" || echo "WARNING: scontrol requeue failed; re-submit manually."
  exit 0
fi

echo "FINAL Job Complete"
echo "======================================"
EOF

# ---------------------------------------------------------
# Write OBS slurm
# ---------------------------------------------------------
cat > "${MODEL_DIR}/${OBS_SLURM}" <<EOF
#!/bin/bash
#SBATCH --job-name=${JOB_PREFIX}_cmfgen_OBS${tag:+_${tag}}
#SBATCH --output=cmf_OBS_%j.out
#SBATCH --error=cmf_OBS_%j.err
#SBATCH --time=${obs_timelimit}
#SBATCH --partition=${part}
$([ "${part}" != "interactive" ] && echo "#SBATCH --qos=n0004")

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${cpus}
#SBATCH --mem=${mem}
#SBATCH --chdir=${MODEL_DIR}/obs

set -euo pipefail

# Runtime check: SN_HYDRO_DATA must exist (created by post-STELLA conversion job)
if [[ ! -f "${MODEL_DIR}/SN_HYDRO_DATA" ]]; then
  echo "ERROR: SN_HYDRO_DATA not found in ${MODEL_DIR} at runtime." >&2
  echo "  This snapshot may be beyond the epoch STELLA reached before crashing." >&2
  exit 1
fi

echo "======================================"
echo "CMFGEN OBS Job Started"
echo "Job ID: \$SLURM_JOB_ID"
echo "Partition: \$SLURM_JOB_PARTITION"
echo "======================================"

module purge
module load intel/2024.0
module load mkl/2026.0
module load python-waterboa/2024.06
# Explicitly set MKL library path (workaround for module loading issues on compute nodes)
export LD_LIBRARY_PATH=/viper/u1/system/soft/RHEL_9/packages/x86_64/intel_oneapi/2026.0/mkl/2026.0/lib/intel64:\${LD_LIBRARY_PATH:+:\$LD_LIBRARY_PATH}

export cmfdist="${CMFDIST}"
export ATOMIC="${ATOMIC}"
export CMFGEN_MODEL_TYPE="${MODEL_TYPE}"
export SCRATCH_DIR="SCRATCH/"
mkdir -p SCRATCH/
export atomic="${ATOMIC}"

echo "Preparing batobs.sh..."
cp -f "${OBS_TEMPLATE_DIR}/batobs.sh" batobs.sh
cp -f "${OBS_TEMPLATE_DIR}/CMF_FLUX_PARAM_INIT" CMF_FLUX_PARAM_INIT

# N_TYPE: Use N_ON_J for OBS formal solution with complex non-monotonic velocity.
# G_ONLY was tested and caused "Mean intensity blowing up" at the inner boundary
# because it ignores the N/J anisotropy term that is essential for non-monotonic
# velocity profiles (Dessart+2015 interacting SNe). N_ON_J includes the full
# radiation field anisotropy, which stabilizes the formal solution.
# The template CMF_FLUX_PARAM_INIT already has N_ON_J — no patching needed.

chmod +x batobs.sh

# ---------------------------
# Resource settings (from SLURM)
# ---------------------------
CPUS=\${SLURM_CPUS_PER_TASK:-8}
export OMP_NUM_THREADS=\$CPUS
export MKL_NUM_THREADS=\$CPUS
export OMP_STACKSIZE=512M
export KMP_STACKSIZE=512M

echo "CPUs: \$CPUS"

echo "Checking required inputs for BATOBS..."
if [[ ! -f "../RVTJ" ]]; then
  echo "ERROR: missing ../RVTJ (BATOBS cannot run). This usually means FINAL did not complete or outputs were cleaned." >&2
  echo "  In obs/ workdir: $(pwd)" >&2
  echo "  Parent dir: $(cd .. && pwd)" >&2
  exit 2
fi

echo "Running batobs.sh..."
# Setup SCRATCH/ symlinks for OBS scratch files
mkdir -p SCRATCH
for _f in POINT1 POINT2 SCRTEMP STEQ_VALS J_COMP JEW EWDATA CFDAT_OUT CONT_FREQ; do
  [[ -e "\$_f" && ! -L "\$_f" ]] && mv "\$_f" SCRATCH/ 2>/dev/null || true
  ln -sf "SCRATCH/\$_f" "\$_f" 2>/dev/null || true
done
./batobs.sh

# batobs.sh should create OBSFRAME and then rename it -> obs_fin_10. Catch silent failures.
if [[ ! -s "obs_fin_10" ]]; then
  echo "ERROR: BATOBS did not produce obs_fin_10 (missing OBSFRAME output)." >&2
  if [[ -f OUT_FLUX ]]; then
    echo "---- OUT_FLUX (tail) ----" >&2
    tail -n 80 OUT_FLUX >&2 || true
    echo "-------------------------" >&2
  fi
  if [[ -f batobs.log ]]; then
    echo "---- batobs.log (tail) ----" >&2
    tail -n 80 batobs.log >&2 || true
    echo "---------------------------" >&2
  fi
  exit 3
fi

# Plot observed spectrum and copy output after BATOBS finishes (non-fatal if plotting fails).
PLOTS_DIR="${MODEL_DIR}/plots"
OUTPUT_DIR="${MODEL_DIR}/output"
CENTRAL_OUTPUTS="${PROJECT_ROOT}/outputs"
mkdir -p "\$PLOTS_DIR" "\$OUTPUT_DIR" "\$CENTRAL_OUTPUTS"

echo "Plotting observed spectrum (post-BATOBS)..."
plotter="${PROJECT_ROOT}/scripts/plot_cmfgen_obs_spectrum.py"
spectrum_plot="\$PLOTS_DIR/${JOB_PREFIX}_obs_spectrum.png"
if [[ -f "\$plotter" ]]; then
  if command -v python3 >/dev/null 2>&1; then
    python3 "\$plotter" "obs_fin_10" --out "\$spectrum_plot" || echo "WARNING: spectrum plot failed (non-fatal)"
  else
    python "\$plotter" "obs_fin_10" --out "\$spectrum_plot" || echo "WARNING: spectrum plot failed (non-fatal)"
  fi
else
  echo "WARNING: plotter not found at \$plotter (skip spectrum plot)"
fi

echo "Copying output spectrum to output directory..."
cp -f obs_fin_10 "\$OUTPUT_DIR/${JOB_PREFIX}_obs_fin.dat" 2>/dev/null || echo "WARNING: failed to copy obs_fin_10 to output dir (non-fatal)"

# Copy to central outputs/ directory with standardized filename
# CMFGEN_MODEL_BASE and EPOCH_LABEL are shell variables set earlier in submit_cmfgen.sh;
# the heredoc expands them when writing the slurm script.
if [[ -s "obs_fin_10" ]]; then
  cp -f obs_fin_10 "\$CENTRAL_OUTPUTS/${CMFGEN_MODEL_BASE}_CMFGEN_${EPOCH_LABEL}.dat" 2>/dev/null \
    && echo "Copied obs_fin_10 -> \$CENTRAL_OUTPUTS/${CMFGEN_MODEL_BASE}_CMFGEN_${EPOCH_LABEL}.dat" \
    || echo "WARNING: failed to copy obs_fin_10 to central outputs/ (non-fatal)"
else
  echo "WARNING: obs_fin_10 missing or empty; skipping central outputs copy"
fi
if [[ -f "\$spectrum_plot" ]]; then
  cp -f "\$spectrum_plot" "\$CENTRAL_OUTPUTS/${CMFGEN_MODEL_BASE}_CMFGEN_${EPOCH_LABEL}.png" 2>/dev/null \
    && echo "Copied plot -> \$CENTRAL_OUTPUTS/${CMFGEN_MODEL_BASE}_CMFGEN_${EPOCH_LABEL}.png" \
    || echo "WARNING: failed to copy plot to central outputs/ (non-fatal)"
fi

if [[ -n "${NOTIFICATION_EMAIL}" ]]; then
  if [[ -s "obs_fin_10" ]]; then
    if command -v mailx >/dev/null 2>&1; then
      printf 'CMFGEN spectrum completed successfully.\nModel: %s\nEpoch: %s\nSpectrum: %s\nJob: %s\n' \
        '${CMFGEN_MODEL_BASE}' '${EPOCH_LABEL}' "\$CENTRAL_OUTPUTS/${CMFGEN_MODEL_BASE}_CMFGEN_${EPOCH_LABEL}.dat" "\$SLURM_JOB_ID" \
        | mailx -s 'CMFGEN spectrum complete: ${CMFGEN_MODEL_BASE} ${EPOCH_LABEL}' '${NOTIFICATION_EMAIL}'
      echo "Notification email sent to ${NOTIFICATION_EMAIL}"
    else
      echo "ERROR: notification requested but mailx is unavailable" >&2
    fi
  else
    echo "WARNING: obs_fin_10 missing or empty; skipping notification email" >&2
  fi
fi

# Plot input Ne vs CMFGEN final Ne (non-fatal if CUR_MODEL_DATA missing or plot fails)
echo "Plotting Ne comparison (input vs CMFGEN final)..."
ne_plotter="${PROJECT_ROOT}/scripts/plot_cmfgen_ne_comparison.py"
cmfgen_run_dir="${MODEL_DIR}"
if [[ -f "\$ne_plotter" ]] && [[ -f "\$cmfgen_run_dir/CUR_MODEL_DATA" ]]; then
  python3 "\$ne_plotter" "\$cmfgen_run_dir" || echo "WARNING: Ne comparison plot failed (non-fatal)"
else
  echo "WARNING: skipping Ne comparison (plotter or CUR_MODEL_DATA missing)"
fi

echo "BATOBS Job Complete"
echo "======================================"
EOF

# ---------------------------------------------------------
# Write SINGLE_FINAL slurm (LTE cold start; no prior INIT)
# ---------------------------------------------------------
# Define optional flags for prepare_VADAT.py (set -u safe for heredoc use)
MONOTONIZE_FLAG=""
if [[ "${MONOTONIZE_VELOCITY:-0}" == "1" ]]; then
  MONOTONIZE_FLAG="--monotonize-velocity"
fi
CHEMICAL_FLAG=""
if [[ "${CHEMICAL_N_INSERT:-0}" != "0" ]]; then
  CHEMICAL_FLAG="--chemical-n-insert=${CHEMICAL_N_INSERT} --chemical-jump-threshold=${CHEMICAL_JUMP_THRESHOLD:-0.1}"
fi
REVERSAL_FLAG=""
if [[ "${REVERSAL_N_INSERT:-0}" != "0" ]]; then
  REVERSAL_FLAG="--reversal-n-insert=${REVERSAL_N_INSERT}"
fi

cat > "${MODEL_DIR}/${SINGLE_FINAL_SLURM}" <<EOF
#!/bin/bash
#SBATCH --job-name=${JOB_PREFIX}_cmfgen${tag:+_${tag}}
#SBATCH --output=cmf_%j.out
#SBATCH --error=cmf_%j.err
#SBATCH --time=${timelimit}
#SBATCH --partition=${part}
$([ "${part}" != "interactive" ] && echo "#SBATCH --qos=n0004")

$([ "${FINAL_WALL_GRACEFUL}" != "0" ] && echo "#SBATCH --signal=B:USR1@${FINAL_SIGNAL_BEFORE_SEC}")
#SBATCH --requeue
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${cpus}
#SBATCH --mem=${mem}
#SBATCH --chdir=${MODEL_DIR}

set -euo pipefail

# Runtime check: SN_HYDRO_DATA must exist (created by post-STELLA conversion job)
if [[ ! -f "${MODEL_DIR}/SN_HYDRO_DATA" ]]; then
  echo "ERROR: SN_HYDRO_DATA not found in ${MODEL_DIR} at runtime." >&2
  echo "  This snapshot may be beyond the epoch STELLA reached before crashing." >&2
  exit 1
fi

check_outgen_for_fatal() {
  local outgen="OUTGEN"
  if [[ ! -f "\$outgen" ]]; then
    echo "ERROR: missing \$outgen after batch.sh (CMFGEN did not write OUTGEN?)" >&2
    return 2
  fi
  if grep -v "recommend --no-relativistic to prevent" "\$outgen" | grep -Eq "NaN|Error in MON_INTERP|Segmentation fault|SIG(ABRT|SEGV)|floating invalid|Error in comp_j_blank\\.f|Mean intensity blowing up|bad initial population guesses|Predicted changes are too large|Error in SETPOP_WLD_V2|zero RGU"; then
    echo "ERROR: fatal pattern detected in \$outgen (NaN / MON_INTERP / crash). Aborting stage." >&2
    echo "---- OUTGEN tail (last 80 lines) ----" >&2
    tail -n 80 "\$outgen" >&2 || true
    echo "------------------------------------" >&2
    return 3
  fi
}

# Track convergence across requeues in a small file. CMFGEN appends to OUTGEN
# with APPEND mode, so OUTGEN grows unboundedly across requeues (observed
# 433 GB after 3 requeues). We rotate OUTGEN each requeue and track the max
# population change in CONVERGENCE_HISTORY to detect divergence.
CONV_HIST="CONVERGENCE_HISTORY"
MAX_REQUEUE=${FINAL_REQUEUE_MAX:-500}

# Extract the last "Maximum % increase" value from OUTGEN using tail (fast
# even on multi-hundred-GB files because tail seeks to end-of-file).
track_convergence() {
  local outgen="OUTGEN"
  [[ -f "\$outgen" ]] || return 0
  local val
  val=\$(tail -n 200 "\$outgen" 2>/dev/null | grep -m1 "Maximum % increase at depth" | sed 's/.*is *\\([0-9.E+]*\\).*/\\1/' 2>/dev/null || true)
  if [[ -n "\$val" ]]; then
    local iter
    iter=\$(tail -n 200 "\$outgen" 2>/dev/null | grep -m1 "Current great iteration count" | sed 's/.*count is *\\([0-9]*\\).*/\\1/' 2>/dev/null || echo "?")
    echo "\${iter} \${val}" >> "\${CONV_HIST}"
    echo "Convergence tracking: iteration \${iter}, max % increase = \${val}%"
  fi
}

# Check CONVERGENCE_HISTORY for divergence. If the last N entries all show
# max change > DIVERGENCE_THRESHOLD and are not decreasing, the model is
# diverging and requeuing will not help.
check_convergence_for_divergence() {
  [[ -f "\${CONV_HIST}" ]] || return 0
  local threshold=50.0  # % change above which we consider the model diverging
  local n_required=5    # need this many consecutive high values to flag divergence
  local n_high=0
  local last_vals=""
  while read -r iter val; do
    [[ -z "\$val" ]] && continue
    last_vals="\${last_vals} \${val}"
    if (( \$(echo "\$val > \$threshold" | bc -l 2>/dev/null || echo 0) )); then
      n_high=\$((n_high + 1))
    else
      n_high=0  # reset on any non-high value
    fi
  done < "\${CONV_HIST}"
  if [[ \$n_high -ge \$n_required ]]; then
    echo "ERROR: CMFGEN model is DIVERGING — last max % changes:\${last_vals}" >&2
    echo "  \${n_high} consecutive iterations with max change > \${threshold}%." >&2
    echo "  Requeuing will not help. Stopping to prevent resource waste." >&2
    echo "  Possible causes:" >&2
    echo "    - Atomic model too complex (too many levels in MODEL_SPEC)" >&2
    echo "    - Bad initial populations (check proxy_populations)" >&2
    echo "    - Non-convergent Ng acceleration parameters (DO_NG, BEG_NG, IBEG_NG)" >&2
    echo "    - BA solver falling back to LAMBDA iteration (check OUTGEN)" >&2
    echo "  Fix: simplify MODEL_SPEC (reduce levels), check VADAT convergence params." >&2
    return 4
  fi
  return 0
}

# Rotate OUTGEN to prevent unbounded growth across requeues.
# CMFGEN opens OUTGEN with APPEND mode (cmfgen.f:99), so without rotation
# the file grows by ~150 GB per 24h requeue cycle.
rotate_outgen() {
  local outgen="OUTGEN"
  [[ -f "\$outgen" ]] || return 0
  local size
  size=\$(stat -c%s "\$outgen" 2>/dev/null || echo 0)
  # Rotate if larger than 2 GB
  if [[ \$size -gt 2147483648 ]]; then
    echo "Rotating OUTGEN (size: \${size} bytes) — saving tail to OUTGEN.tail"
    tail -n 5000 "\$outgen" > OUTGEN.tail 2>/dev/null || true
    : > "\$outgen"  # Truncate (CMFGEN will append to the empty file)
  fi
}

echo "======================================"
echo "CMFGEN SINGLE-FINAL Job Started (LTE cold start)"
echo "Job ID: \$SLURM_JOB_ID"
echo "Partition: \$SLURM_JOB_PARTITION"
echo "CPUs: \$SLURM_CPUS_PER_TASK"
echo "Node: \$SLURMD_NODENAME"
echo "======================================"

restart_count="\${SLURM_RESTART_COUNT:-0}"
requeue_idx=\$((restart_count + 1))
final_job_name="${JOB_PREFIX}_cmfgen_N\${requeue_idx}"
if command -v scontrol >/dev/null 2>&1; then
  scontrol update JobId="\$SLURM_JOB_ID" JobName="\$final_job_name" >/dev/null 2>&1 || true
fi

module purge
module load intel/2024.0
module load mkl/2026.0
module load python-waterboa/2024.06
# Explicitly set MKL library path (workaround for module loading issues on compute nodes)
export LD_LIBRARY_PATH=/viper/u1/system/soft/RHEL_9/packages/x86_64/intel_oneapi/2026.0/mkl/2026.0/lib/intel64:\${LD_LIBRARY_PATH:+:\$LD_LIBRARY_PATH}

export cmfdist="${CMFDIST}"
export ATOMIC="${ATOMIC}"
export CMFGEN_MODEL_TYPE="${MODEL_TYPE}"
export SCRATCH_DIR="SCRATCH/"
mkdir -p SCRATCH/

NUM_ITS_RESET_SENTINEL=".CMFGEN_FINAL_NEEDS_NUM_ITS_RESET"

request_graceful_stop() {
  echo ""
  echo "======================================================"
  echo "Pre-timeout signal received: requesting graceful stop."
  echo "Setting NUM_ITS=0 in IN_ITS so CMFGEN exits cleanly."
  echo "======================================================"
  python - <<'PY'
from pathlib import Path
p = Path("IN_ITS")
lines = p.read_text().splitlines() if p.exists() else []
out = []
done = False
for ln in lines:
    if "[NUM_ITS]" in ln:
        out.append("0 [NUM_ITS] ")
        done = True
    else:
        out.append(ln)
if not done:
    out.insert(0, "0 [NUM_ITS] ")
p.write_text("\n".join(out) + "\n")
PY
  touch .CMFGEN_REQUEST_REQUEUE
  touch "\${NUM_ITS_RESET_SENTINEL}"
}
trap request_graceful_stop USR1

graceful_timer_pid=""
# stop_graceful_timer: always defined; no-op when timer is disabled.
stop_graceful_timer() {
  if [[ -n "\${graceful_timer_pid}" ]] && kill -0 "\${graceful_timer_pid}" 2>/dev/null; then
    kill "\${graceful_timer_pid}" 2>/dev/null || true
  fi
}
# Start fallback sleep timer unless disabled (debug mode sets to 999999999).
if [[ ${FINAL_GRACEFUL_SLEEP_SEC} -lt 100000000 ]]; then
  (
    sleep ${FINAL_GRACEFUL_SLEEP_SEC}
    request_graceful_stop
  ) &
  graceful_timer_pid=\$!
  echo "Graceful timer started: fires in ${FINAL_GRACEFUL_SLEEP_SEC}s (PID \${graceful_timer_pid})"
else
  echo "Graceful timer disabled (FINAL_GRACEFUL_SLEEP_SEC=${FINAL_GRACEFUL_SLEEP_SEC})"
fi

CPUS=\${SLURM_CPUS_PER_TASK:-8}
export OMP_NUM_THREADS=\$CPUS
export MKL_NUM_THREADS=\$CPUS
export OMP_STACKSIZE=512M
export KMP_STACKSIZE=512M
echo "CPUs: \$CPUS"

# On the very first run (not a requeue), do LTE cold-start prep.
# No hard reset needed for fresh models (copied from base/).
# On requeues (SLURM_RESTART_COUNT>0), skip prep so we continue from checkpoint.
if [[ "\${SLURM_RESTART_COUNT:-0}" -eq 0 && ! -f .CMFGEN_SINGLE_FINAL_PREPPED ]]; then
  echo "First run: LTE cold-start prep..."
  rm -f batch.log OUTGEN
  echo "MODEL_TYPE=${MODEL_TYPE} INC_RAD_DECAYS=${INC_RAD_DECAYS}"
  # Print flag status for logging
  if [[ "${NO_RELATIVISTIC}" == "1" ]]; then
    echo "NO_RELATIVISTIC=1 (passed --no-relativistic to prepare_VADAT.py)"
  fi
  python "${VADAT_SCRIPT}" --single-final --convergence-mode stable --vel-law "${VEL_LAW}" ${MODEL_TYPE_FLAG} ${INC_RAD_DECAYS_FLAG} ${TARGET_ND_FLAG} ${USE_HYDRO_GRID_FLAG} ${SMOOTH_JUMPS_FLAG} ${BLEND_FLAG} ${MONOTONIZE_FLAG:-} ${CHEMICAL_FLAG:-} ${REVERSAL_FLAG:-} ${NO_REL_FLAG} --vadat-overrides "${CONFIG_FILE}" > 'batch.log'
  python "${INITS_SCRIPT}" --single-final --num-its ${FINAL_NUM_ITS} >> 'batch.log'
  if [[ "${PROXY_POPULATIONS}" == "1" || "${PROXY_POPULATIONS_V2}" == "1" ]]; then
    echo "PROXY_POPULATIONS: generating proxy dilute-LTE _IN files and overriding VADAT..."
    # Create atomic data symlinks needed by proxy_populations (normally done by batch.sh)
    if [[ ! -e HI_F_OSCDAT ]]; then
      echo "  Creating atomic data symlinks for proxy_populations..."
      ln -sf "\$ATOMIC/HYD/I/5dec96/hi_osc.dat"              HI_F_OSCDAT
      ln -sf "\$ATOMIC/HE/I/5dec96/heioscdat_a7.dat"         HeI_F_OSCDAT
      ln -sf "\$ATOMIC/HE/II/5dec96/he2_osc.dat"             He2_F_OSCDAT
      ln -sf "\$ATOMIC/OXY/I/19apr23/osc_data"               OI_F_OSCDAT
      ln -sf "\$ATOMIC/OXY/II/19apr23/osc_data"              O2_F_OSCDAT
      ln -sf "\$ATOMIC/OXY/III/19apr23/osc_data"             OIII_F_OSCDAT
      ln -sf "\$ATOMIC/SIL/I/19apr23/osc_data"               SkI_F_OSCDAT
      ln -sf "\$ATOMIC/SIL/II/19apr23/osc_data"              Sk2_F_OSCDAT
      ln -sf "\$ATOMIC/SIL/III/19apr23/osc_data"             SkIII_F_OSCDAT
      ln -sf "\$ATOMIC/SIL/IV/19apr23/osc_data"              SkIV_F_OSCDAT
      ln -sf "\$ATOMIC/CA/I/19apr23/osc_data"                CaI_F_OSCDAT
      ln -sf "\$ATOMIC/CA/II/19apr23/osc_data"               Ca2_F_OSCDAT
      ln -sf "\$ATOMIC/CA/III/19apr23/osc_data"              CaIII_F_OSCDAT
      ln -sf "\$ATOMIC/FE/I/19apr23/osc_data"                FeI_F_OSCDAT
      ln -sf "\$ATOMIC/FE/II/19apr23/osc_data"               Fe2_F_OSCDAT
      ln -sf "\$ATOMIC/FE/III/7feb05/FeIII_OSC"              FeIII_F_OSCDAT
      ln -sf "\$ATOMIC/FE/IV/19apr23/feiv_osc"               FeIV_F_OSCDAT
      ln -sf "\$ATOMIC/FE/V/19apr23/osc_data"                FeV_F_OSCDAT
      ln -sf "\$ATOMIC/CARB/I/19apr23/osc_data"              CI_F_OSCDAT
      ln -sf "\$ATOMIC/CARB/II/19apr23/osc_data"             C2_F_OSCDAT
      # F_TO_S symlinks for proxy_populations (needed before batch.sh runs)
      ln -sf "\$ATOMIC/HYD/I/5dec96/hi_f_to_s_15.dat"        HI_F_TO_S
      ln -sf "\$ATOMIC/HE/I/5dec96/hei_f_to_s_a7_ext.dat"    HeI_F_TO_S
      ln -sf "\$ATOMIC/HE/II/5dec96/he2_f_to_s_10.dat"       He2_F_TO_S
      ln -sf "\$ATOMIC/OXY/I/19apr23/f_to_s_58"              OI_F_TO_S
      ln -sf "\$ATOMIC/OXY/II/19apr23/f_to_s_137"            O2_F_TO_S
      ln -sf "\$ATOMIC/OXY/III/19apr23/f_to_s_107"           OIII_F_TO_S
      ln -sf "\$ATOMIC/SIL/I/19apr23/f_to_s_ls"              SkI_F_TO_S
      ln -sf "/viper/ptmp1/sbrennan/SN2024pba/atomic/Sk2_f_to_s_no_3s3p2"  Sk2_F_TO_S
      ln -sf "\$ATOMIC/SIL/III/19apr23/f_to_s_split"         SkIII_F_TO_S
      ln -sf "\$ATOMIC/SIL/IV/19apr23/f_to_s_split"          SkIV_F_TO_S
      ln -sf "\$ATOMIC/CA/I/19apr23/f_to_s_term"             CaI_F_TO_S
      ln -sf "\$ATOMIC/CA/II/19apr23/f_to_s_sm"              Ca2_F_TO_S
      ln -sf "\$ATOMIC/CA/III/19apr23/f_to_s"                CaIII_F_TO_S
      ln -sf "\$ATOMIC/FE/I/19apr23/f_to_s_181"              FeI_F_TO_S
      ln -sf "\$ATOMIC/FE/II/19apr23/f_to_s_135"             Fe2_F_TO_S
      ln -sf "\$ATOMIC/FE/III/19apr23/f_to_s_105"            FeIII_F_TO_S
      ln -sf "\$ATOMIC/FE/IV/19apr23/f_to_s_110"             FeIV_F_TO_S
    fi
    if [[ -f ./SN_HYDRO_DATA ]]; then
      PROXY_HYDRO="./SN_HYDRO_DATA"
    else
      PROXY_HYDRO="../SN_HYDRO_DATA"
    fi
    if [[ "${PROXY_POPULATIONS_V2}" == "1" ]]; then
      _v2_strategy="LTE"
      if [[ "${PROXY_POPULATIONS_NLTE_SOLVER}" == "1" ]]; then
        _v2_strategy="GREY"
      fi
      python "${PROJECT_ROOT}/scripts/proxy_populations_v2.py" . --hydro-data "\${PROXY_HYDRO}" --output-dir . --plots --strategy "\${_v2_strategy}" --min-ion-frac "${PROXY_POPULATIONS_MIN_ION_FRAC:-1.0e-30}" --mass-frac-floor "${PROXY_POPULATIONS_MASS_FRAC_FLOOR:-1.0e-8}" ${PROXY_POPULATIONS_T_EXC_MIN:+--t-exc-min "${PROXY_POPULATIONS_T_EXC_MIN}"} --write-t-ion --smooth --t-ion-source stella ${PROXY_COMPLEX_MODEL:+--complex-model} >> 'batch.log' || true
    else
      python "${PROXY_POPULATIONS_SCRIPT}" . --hydro-data "\${PROXY_HYDRO}" --plot ${PROXY_POPULATIONS_MIN_ION_FRAC:+--min-ion-frac "${PROXY_POPULATIONS_MIN_ION_FRAC}"} ${PROXY_POPULATIONS_NLTE_SOLVER:+--nlte-solver} ${PROXY_POPULATIONS_GREY_MAX_ITER:+--grey-max-iter "${PROXY_POPULATIONS_GREY_MAX_ITER}"} ${PROXY_POPULATIONS_GREY_TOL:+--grey-tol "${PROXY_POPULATIONS_GREY_TOL}"} >> 'batch.log'
    fi
    python - <<'PYPROXY'
from pathlib import Path
import numpy as np
p = Path("VADAT")
# DC_METH selection: TR interpolates in T space, but REGRID_LOG_DC_V1
# has a known issue with non-monotonic T (regrid_log_dc_v1.f:284:
# "This code may have issues if T is very nonmonotonic"). With multiple
# T matches per depth, the TR option picks the wrong interval, producing
# pathological DCs and 10^22+ population corrections on iteration 1.
# For non-monotonic T models, use DC_METH=R (radius interpolation) instead.
dc_meth = 'TR'
hydro_path = Path("SN_HYDRO_DATA")
if hydro_path.exists():
    try:
        lines_h = hydro_path.read_text().split('\n')
        in_sec = False
        t_data = []
        for ln in lines_h:
            s = ln.strip()
            if not s: continue
            if 'Temperature' in s:
                in_sec = True
                continue
            if in_sec:
                if s.startswith('!') and 'mass' not in s.lower():
                    in_sec = False
                    continue
                try:
                    t_data.extend([float(x) for x in s.split()])
                except ValueError:
                    in_sec = False
        if t_data and len(t_data) > 1:
            t_arr = np.array(t_data)
            n_non_mono = int(np.sum(np.diff(t_arr) > 0))
            if n_non_mono > len(t_arr) * 0.1:
                dc_meth = 'R'
                print(f"  [DC_METH] T is non-monotonic ({n_non_mono}/{len(t_arr)-1} zones), "
                      f"using DC_METH=R instead of TR to avoid REGRID_LOG_DC_V1 interpolation failure")
    except Exception:
        pass
lines = p.read_text().splitlines()
# Pre-scan: check if LTE_EST=T was explicitly set by prepare_VADAT.py (via vadat_overrides)
lte_est_is_t = any(
    ln.split("[")[0].strip().upper() == "T" and "[LTE_EST]" in ln
    for ln in lines
)
out = []
for ln in lines:
    if "[LTE_EST]" in ln:
        if lte_est_is_t:
            out.append(ln)  # Keep LTE_EST=T (explicitly requested via vadat_overrides)
        else:
            out.append("F [LTE_EST]")
    elif "[LIN_INT]" in ln:
        if not lte_est_is_t:
            out.append("T [LIN_INT]")
        else:
            out.append("F [LIN_INT]")  # LIN_INT=F with LTE_EST=T (LTE cold start, no _IN interpolation)
    elif "[DO_DDT]" in ln:
        out.append("F [DO_DDT]")
    elif "[USE_DJDT_RTE]" in ln:
        out.append("F [USE_DJDT_RTE]")
    elif "[DC_METH]" in ln:
        out.append(f"{dc_meth} [DC_METH]")  # Always apply dc_meth (REGRID_LOG_DC_V1 uses DC_METH even with LTE_EST=T)
    # NOTE: USE_J_REL, USE_JREL, MAX_LAM, MAX_LIN, IBEG_NG, DO_NG are NOT overridden here.
    # prepare_VADAT.py sets them correctly per profile and get_rel_config().
    else:
        out.append(ln)
p.write_text("\n".join(out) + "\n")
if lte_est_is_t:
    print(f"  VADAT: LTE_EST=T preserved (LTE cold start; proxy _IN files generated but ignored by LTE_EST)")
else:
    print(f"  VADAT: LTE_EST=F, LIN_INT=T, DO_DDT=F, USE_DJDT_RTE=F, DC_METH={dc_meth} (proxy _IN files; MAX_LAM/MAX_LIN/DO_NG from prepare_VADAT.py)")
PYPROXY
  fi
  touch .CMFGEN_SINGLE_FINAL_PREPPED
else
  echo "Requeue \${requeue_idx}: restoring checkpoint state..."

  # Check for divergence before continuing. If the model is diverging,
  # stop instead of wasting more compute on a non-convergent cycle.
  check_convergence_for_divergence
  div_rc=\$?
  if [[ \$div_rc -ne 0 ]]; then
    echo "ERROR: Divergence detected on requeue \${requeue_idx}. Aborting." >&2
    exit \$div_rc
  fi

  # Check max requeue count to prevent infinite requeueing.
  if [[ "\${SLURM_RESTART_COUNT:-0}" -ge "\${MAX_REQUEUE}" ]]; then
    echo "ERROR: Maximum requeue count (\${MAX_REQUEUE}) reached. Stopping." >&2
    echo "  The model has not converged after \${SLURM_RESTART_COUNT} requeues." >&2
    echo "  Check CONVERGENCE_HISTORY and OUTGEN.tail for diagnostics." >&2
    exit 5
  fi

  # Rotate OUTGEN to prevent unbounded growth across requeues.
  # Without this, OUTGEN can grow to hundreds of GB (observed 433 GB).
  rotate_outgen

  "\$cmfdist/com/out2in.sh"
  # Switch from LTE cold-start to warm-start: read departure coefficients
  # from the _IN files (created by out2in.sh from the previous run's OUT files).
  # LTE_EST=T forces DC_INTERP_METHOD='LTE' (rd_control_variables.f:847),
  # which ignores _IN files and restarts from LTE — losing all convergence.
  # LIN_INT=T (GRID=T) enables direct interpolation of _IN files.
  # DC_METH=TR interpolates in T space (manual p.96: preferred when GRID=F),
  # but for non-monotonic T, use DC_METH=R (radius interpolation) instead.
  python - <<'PYREQUEUE'
from pathlib import Path
import numpy as np
p = Path("VADAT")
# Check T monotonicity from SN_HYDRO_DATA
dc_meth = 'TR'
hydro_path = Path("SN_HYDRO_DATA")
if hydro_path.exists():
    try:
        lines_h = hydro_path.read_text().split('\n')
        in_sec = False
        t_data = []
        for ln in lines_h:
            s = ln.strip()
            if not s: continue
            if 'Temperature' in s:
                in_sec = True
                continue
            if in_sec:
                if s.startswith('!') and 'mass' not in s.lower():
                    in_sec = False
                    continue
                try:
                    t_data.extend([float(x) for x in s.split()])
                except ValueError:
                    in_sec = False
        if t_data and len(t_data) > 1:
            t_arr = np.array(t_data)
            n_non_mono = int(np.sum(np.diff(t_arr) > 0))
            if n_non_mono > len(t_arr) * 0.1:
                dc_meth = 'R'
    except Exception:
        pass
lines = p.read_text().splitlines() if p.exists() else []
out = []
for ln in lines:
    if "[LTE_EST]" in ln:
        out.append("F [LTE_EST]")
    elif "[LIN_INT]" in ln:
        out.append("T [LIN_INT]")
    elif "[DC_METH]" in ln:
        out.append(f"{dc_meth} [DC_METH]")
    else:
        out.append(ln)
p.write_text("\n".join(out) + "\n")
print(f"  Requeue: LTE_EST=F, LIN_INT=T, DC_METH={dc_meth} (warm start from _IN files; MAX_LAM/MAX_LIN unchanged)")
PYREQUEUE
  # Restore NUM_ITS if a graceful stop set it to 0
  if [[ -f "\${NUM_ITS_RESET_SENTINEL}" ]]; then
    echo "Resetting IN_ITS NUM_ITS to ${FINAL_NUM_ITS} for requeued run."
    python - <<PY
from pathlib import Path
p = Path("IN_ITS")
lines = p.read_text().splitlines() if p.exists() else []
out = []
done = False
for ln in lines:
    if "[NUM_ITS]" in ln:
        out.append("${FINAL_NUM_ITS} [NUM_ITS] ")
        done = True
    else:
        out.append(ln)
if not done:
    out.insert(0, "${FINAL_NUM_ITS} [NUM_ITS] ")
p.write_text("\\n".join(out) + "\\n")
PY
    rm -f "\${NUM_ITS_RESET_SENTINEL}"
  fi
fi

# Safety: if IN_ITS still has NUM_ITS=0 bump it
python - <<PY
from pathlib import Path
p = Path("IN_ITS")
lines = p.read_text().splitlines() if p.is_file() else []
num = None
for ln in lines:
    if "[NUM_ITS]" in ln:
        parts = ln.split()
        if parts:
            try: num = int(parts[0])
            except ValueError: num = None
        break
if num != 0:
    raise SystemExit(0)
out = []
done = False
for ln in lines:
    if "[NUM_ITS]" in ln:
        out.append("${FINAL_NUM_ITS} [NUM_ITS] ")
        done = True
    else:
        out.append(ln)
if not done:
    out.insert(0, "${FINAL_NUM_ITS} [NUM_ITS] ")
p.write_text("\\n".join(out) + "\\n")
print("IN_ITS had NUM_ITS=0; set to ${FINAL_NUM_ITS} before batch.sh.")
PY

echo "Running batch.sh..."
# Setup SCRATCH/ symlinks for files that don't use SCRATCH_DIR
mkdir -p SCRATCH
for _f in POINT1 POINT2 SCRTEMP STEQ_VALS J_COMP JEW EWDATA CFDAT_OUT CONT_FREQ LINEHEAT NETRATE TOTRATE TRANS_INFO ADIABAT_CHK; do
  [[ -e "\$_f" && ! -L "\$_f" ]] && mv "\$_f" SCRATCH/ 2>/dev/null || true
  ln -sf "SCRATCH/\$_f" "\$_f" 2>/dev/null || true
done
tcsh -f ./batch.sh
check_outgen_for_fatal
stop_graceful_timer

# Track convergence: extract last max % change from OUTGEN and append to
# CONVERGENCE_HISTORY. This small file persists across requeues and is
# used to detect divergence before the next requeue.
track_convergence

echo "Preparing restart inputs (out2in.sh)..."
"\$cmfdist/com/out2in.sh"

if [[ -f .CMFGEN_REQUEST_REQUEUE ]]; then
  # Check for divergence before requeuing. If the model is diverging,
  # stop instead of wasting more compute on a non-convergent cycle.
  check_convergence_for_divergence
  div_rc=\$?
  if [[ \$div_rc -ne 0 ]]; then
    echo "ERROR: Divergence detected after batch.sh. NOT requeuing." >&2
    exit \$div_rc
  fi
  echo "Requeue requested; requeueing job \$SLURM_JOB_ID..."
  rm -f .CMFGEN_REQUEST_REQUEUE
  scontrol requeue "\$SLURM_JOB_ID" || echo "WARNING: scontrol requeue failed; re-submit manually."
  exit 0
fi

echo "SINGLE-FINAL Job Complete"
echo "======================================"
EOF

chmod +x "${MODEL_DIR}/${FINAL_SLURM}" "${MODEL_DIR}/${SINGLE_FINAL_SLURM}" "${MODEL_DIR}/${OBS_SLURM}"

echo ""
echo "Generated Slurm scripts:"
echo "  ${MODEL_DIR}/${FINAL_SLURM}"
echo "  ${MODEL_DIR}/${SINGLE_FINAL_SLURM}"
echo "  ${MODEL_DIR}/${OBS_SLURM}"
echo ""

if [[ $no_submit -eq 1 ]]; then
  echo "Use --no-submit flag: Files ready, run sbatch manually."
  exit 0
fi

# ---------------------------------------------------------
# Ensure only one set of CMFGEN jobs for this model at a time
# ---------------------------------------------------------
# Match CMFGEN and OBS jobs for this model prefix.
# CMFGEN jobs: ${JOB_PREFIX}_cmfgen or ${JOB_PREFIX}_cmfgen_N{idx} (requeued)
# OBS jobs: ${JOB_PREFIX}_cmfgen_OBS
existing=$(squeue -u "${USER:-$(whoami)}" -h -o "%j" 2>/dev/null | grep -E "^${JOB_PREFIX}_cmfgen(_N[0-9]+|_OBS|_${tag})?$" || true)
if [[ -n "$existing" ]]; then
  echo "ERROR: CMFGEN jobs already running or pending for this model (${JOB_PREFIX}):"
  echo "$existing"
  echo "Wait for them to finish or cancel with: scancel <job_id>"
  echo "Then run submit_cmfgen again."
  exit 1
fi

# ---------------------------------------------------------
# Submit with dependencies
# ---------------------------------------------------------
# Build dependency flag for first job if HERACLES dependency is provided
DEP_FLAG=""
if [[ -n "$DEPENDENCY_JOB" ]]; then
  DEP_FLAG="--dependency=afterok:${DEPENDENCY_JOB}"
  echo "First CMFGEN job depends on HERACLES job: $DEPENDENCY_JOB"
fi

if [[ $restart -eq 1 ]]; then
  # Warm restart: use FINAL_SLURM (reads *_IN populations from previous run)
  echo "Restart mode: submitting CMFGEN (warm restart) and OBS..."

  jid_final=$(sbatch --parsable $DEP_FLAG "${MODEL_DIR}/${FINAL_SLURM}")
  echo "CMFGEN Job ID: $jid_final"

  if [[ "${SKIP_OBS:-0}" != "1" ]]; then
    jid_obs=$(sbatch --parsable --dependency=afterok:"$jid_final" "${MODEL_DIR}/${OBS_SLURM}")
    echo "OBS Job ID:    $jid_obs (depends on $jid_final)"
  else
    echo "OBS Job ID:    (skipped -- SKIP_OBS=1)"
  fi
else
  # Default: LTE cold-start via SINGLE_FINAL_SLURM
  echo "Single-final mode: submitting CMFGEN (LTE cold start) and OBS..."

  jid_final=$(sbatch --parsable $DEP_FLAG "${MODEL_DIR}/${SINGLE_FINAL_SLURM}")
  echo "CMFGEN Job ID: $jid_final"

  if [[ "${SKIP_OBS:-0}" != "1" ]]; then
    jid_obs=$(sbatch --parsable --dependency=afterok:"$jid_final" "${MODEL_DIR}/${OBS_SLURM}")
    echo "OBS Job ID:    $jid_obs (depends on $jid_final)"
  else
    echo "OBS Job ID:    (skipped -- SKIP_OBS=1)"
  fi
fi