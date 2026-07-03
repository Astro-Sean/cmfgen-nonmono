#!/usr/bin/env bash
# submit_cmfgen.sh
#
# Creates three Slurm scripts in the model directory:
#   - run_cmfgen_INIT.slurm
#   - run_cmfgen_FINAL.slurm   (afterok INIT)
#   - run_cmfgen_OBS.slurm     (afterok FINAL)
#
# Also creates: <MODEL_DIR>/obs/ and copies templates from:
#   $CMFDIST/obs_template/{batobs.sh,CMF_FLUX_PARAM_INIT}
#
# Options:
#   --debug         Use interactive partition preset
#   --all           Submit full INIT->FINAL->OBS chain (STABLE is skipped by default)
#   --with-stable   Include STABLE phase (INIT->STABLE->FINAL->OBS)
#   --single-final  Skip INIT; start FINAL from LTE (default: OFF)
#   --with-init     Disable single-final mode; use INIT->FINAL chain
#   --quick         Request minimal resources (4 CPUs, 32G mem, 12h FINAL) for faster queue turnaround
#   --no-submit     Only write files; do not call sbatch
#   --restart       Only (re)submit FINAL and OBS Slurm jobs (INIT is not submitted)
#   --vary_T        Allow temperature to vary (default: T is held fixed)
#   --fix-T         Fix temperature (default: T is held fixed). Use --no-fix-T to allow T to vary
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
with_stable=0
single_final=0
quick=0
no_submit=0
restart=0
vary_t=0
fix_t=1
VEL_LAW="12"
MODEL_DIR=""
DEPENDENCY_JOB=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --debug) debug=1; shift ;;
    --all) all_jobs=1; shift ;;
    --with-stable) with_stable=1; shift ;;
    --single-final) single_final=1; shift ;;
    --with-init) single_final=0; shift ;;
    --quick) quick=1; shift ;;
    --no-submit) no_submit=1; shift ;;
    --restart) restart=1; shift ;;
    --vary_T) vary_t=1; shift ;;
    --fix-T) fix_t=1; shift ;;
    --no-fix-T) fix_t=0; shift ;;
    --dependency=*) DEPENDENCY_JOB="${1#*=}"; shift ;;
    --vel-law=*) VEL_LAW="${1#*=}"; shift ;;
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
CMFDIST="${CMFDIST:-${HOME}/cur_cmf}"

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
  BASE_FILES_TO_KEEP="batch.sh VADAT_INIT IN_ITS MODEL_SPEC NUC_DECAY_DATA arnaud_rothenflug.dat"
  # Also keep SN_HYDRO_DATA and obs directory
  FILES_TO_KEEP="$BASE_FILES_TO_KEEP SN_HYDRO_DATA obs"
  
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
# Use CMFGEN_JOB_PREFIX if set (from run_model_grid.py), otherwise fall back to FOLDER_NAME
JOB_PREFIX="${CMFGEN_JOB_PREFIX:-${FOLDER_NAME}}"


echo "======================================================"
echo "CMFGEN Submission Script"
echo "Model Directory: $MODEL_DIR"
echo "Log File: $LOGFILE"
echo "Debug Mode: $debug"
echo "All Jobs General: $all_jobs"
echo "With STABLE stage: $with_stable"
echo "Quick (minimal resources): $quick"
echo "No-Submit: $no_submit"
echo "Restart: $restart"
echo "Single-Final (LTE cold start, skip INIT): $single_final"
echo "Fix-T: $fix_t (temperature fixed unless --no-fix-T)"
echo "======================================================"

# ---------------------------
# CMFGEN environment paths
# ---------------------------
CMFDIST="${CMFDIST:-${HOME}/cur_cmf}"
ATOMIC="${HOME}/atomic"

PROJECT_ROOT="$(dirname "$(dirname "$(dirname "$MODEL_DIR")")")"
VADAT_SCRIPT="${PROJECT_ROOT}/scripts/prepare_VADAT.py"
INITS_SCRIPT="${PROJECT_ROOT}/scripts/prepare_IN_INTS.py"
# Model name = name of the models/<MODEL_NAME>/cmfgen directory's parent (i.e. models/<MODEL_NAME>)
MODEL_NAME_FROM_DIR="$(basename "$(dirname "$MODEL_DIR")")"
# Use CMFGEN_CONFIG_FILE if set by submit_cmfgen_jobs.sh; otherwise find the config by name.
# Config files are named <MODEL_NAME>.yaml (not model_<MODEL_NAME>.yaml).
if [[ -n "${CMFGEN_CONFIG_FILE:-}" && -f "${CMFGEN_CONFIG_FILE}" ]]; then
  CONFIG_FILE="${CMFGEN_CONFIG_FILE}"
elif [[ -f "${PROJECT_ROOT}/configs/${MODEL_NAME_FROM_DIR}.yaml" ]]; then
  CONFIG_FILE="${PROJECT_ROOT}/configs/${MODEL_NAME_FROM_DIR}.yaml"
else
  # Legacy fallback: model_<name>.yaml
  CONFIG_FILE="${PROJECT_ROOT}/configs/model_${MODEL_NAME_FROM_DIR}.yaml"
fi
echo "  CONFIG_FILE: ${CONFIG_FILE}"

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

# ---------------------------
# Resource presets (smaller requests = faster queue scheduling on most clusters)
# Increase MEM_DEFAULT / CPUS_DEFAULT if FINAL fails with OOM or needs more threads.
# ---------------------------
PARTITION_DEFAULT="general"
TIME_DEFAULT="24:00:00"
CPUS_DEFAULT=16
# CMFGEN BA / matrix steps often need many GB; too-small requests can OOM during the first BA build.
# Defaults tuned to be "reasonable" while still avoiding the prior 8G OOM failures.
MEM_DEFAULT="${CMFGEN_SLURM_MEM:-32G}"

# INIT can often run with less than FINAL, but still needs enough for the first BA build.
INIT_MEM_DEFAULT="${CMFGEN_INIT_SLURM_MEM:-32G}"

PARTITION_DEBUG="interactive"
TIME_DEBUG="01:59:00"
CPUS_DEBUG=4
MEM_DEBUG="8G"

# Quick preset: minimal resources for faster backfill / shorter queue wait
TIME_QUICK="24:00:00"
CPUS_QUICK=16
MEM_QUICK="${CMFGEN_SLURM_MEM_QUICK:-32G}"

# INIT: default general partition, 24h. Override with CMFGEN_INIT_* env vars.
# Keep memory aligned with main MEM_DEFAULT to avoid early BA OOM.
INIT_PARTITION="${CMFGEN_INIT_PARTITION:-general}"
INIT_TIME="${CMFGEN_INIT_TIME:-24:00:00}"

# STABLE always runs on general partition (longer runtime needed for pure Lambda settling)
STABLE_PARTITION="general"
STABLE_TIME="24:00:00"
# STABLE graceful stop: request clean stop ~1h before walltime (23h into 24h job)
# Fallback at 20h to ensure 4h margin for graceful stop completion
STABLE_SIGNAL_BEFORE_SEC=3600
STABLE_GRACEFUL_SLEEP_SEC=72000

# Keep INIT short to avoid interactive/general timeout pressure.
# Always run 3 iterations in --init mode.
INIT_NUM_ITS=3
# STABLE phase: pure Lambda iterations to settle populations before Non-Lambda
# Default to a short pre-conditioning run; override by editing here if needed.
STABLE_NUM_ITS=10
# Keep FINAL at high-iteration mode; graceful chunking/requeue handles walltime.
FINAL_NUM_ITS=100000

OBS_TIME_DEFAULT="24:00:00"
OBS_TIME_DEBUG="01:30:00"

INIT_SLURM="run_cmfgen_INIT.slurm"
STABLE_SLURM="run_cmfgen_STABLE.slurm"
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
  part="$PARTITION_DEFAULT"
  timelimit="$TIME_DEFAULT"
  cpus="$CPUS_DEFAULT"
  mem="$MEM_DEFAULT"
  init_mem="${INIT_MEM_DEFAULT:-$mem}"
  obs_timelimit="$OBS_TIME_DEFAULT"
  tag=""
fi

varyt_arg=""
# If fix_t=0 (not fixing T), allow T to vary by passing --vary_T
# If fix_t=1 (default), keep T fixed by not passing --vary_T
if [[ $fix_t -eq 0 ]]; then
  varyt_arg="--vary_T"
fi

# FINAL graceful chunking (walltime-dependent):
# - For 24h jobs: request clean stop after ~18h (6h before walltime) via USR1 + 64800s sleep fallback.
# - For --quick (12h): use proportionally earlier stop (~9h runtime) so requeue still fits.
FINAL_SIGNAL_BEFORE_SEC=21600
FINAL_GRACEFUL_SLEEP_SEC=64800
if [[ $quick -eq 1 ]]; then
  FINAL_SIGNAL_BEFORE_SEC=10800
  FINAL_GRACEFUL_SLEEP_SEC=32400
fi

# FINAL continuation policy (written into run_cmfgen_FINAL.slurm):
#   wall         — default: Slurm SIGUSR1 before walltime + sleep fallback; NUM_ITS=0 then scontrol requeue.
#   finish       — no wall preemption; after each clean tcsh -f ./batch.sh exit, scontrol requeue until max restarts
#                  or you create .CMFGEN_FINAL_NO_REQUEUE in the model dir. Does not auto-submit OBS (see submit block).
#   recoverable  — no wall preemption; requeue only if OUTGEN matches a few soft-fail strings. OBS not auto-submitted.
FINAL_REQUEUE_WHEN="${CMFGEN_FINAL_REQUEUE_WHEN:-wall}"
FINAL_WALL_GRACEFUL="${CMFGEN_FINAL_WALL_GRACEFUL:-}"
if [[ -z "${FINAL_WALL_GRACEFUL}" ]]; then
  if [[ "${FINAL_REQUEUE_WHEN}" == "finish" || "${FINAL_REQUEUE_WHEN}" == "recoverable" ]]; then
    FINAL_WALL_GRACEFUL=0
  else
    FINAL_WALL_GRACEFUL=1
  fi
fi
FINAL_REQUEUE_MAX="${CMFGEN_FINAL_REQUEUE_MAX:-500}"

# Debug mode: 2 iterations only, no graceful stop, no requeue.
# Must come AFTER the graceful timer defaults so these take final effect.
if [[ $debug -eq 1 ]]; then
  FINAL_NUM_ITS=2
  FINAL_GRACEFUL_SLEEP_SEC=999999999  # Disable sleep-based fallback timer
  FINAL_WALL_GRACEFUL=0               # Disable #SBATCH --signal + USR1 handler
fi

echo ""
echo "Resources:"
echo "  Partition:       $part"
echo "  Time:            $timelimit"
echo "  CPUs/Task:       $cpus"
echo "  Memory (FINAL/OBS): $mem"
echo "  INIT Memory:     $init_mem"
echo "  INIT Partition:  $INIT_PARTITION"
echo "  INIT Time:       $INIT_TIME"
echo "  INIT NUM_ITS:    $INIT_NUM_ITS"
echo "  STABLE Partition: $STABLE_PARTITION"
echo "  STABLE Time:     $STABLE_TIME"
echo "  STABLE NUM_ITS:  $STABLE_NUM_ITS"
echo "  STABLE graceful: USR1 ${STABLE_SIGNAL_BEFORE_SEC}s before wall; requeue enabled"
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
rm -f "${MODEL_DIR}/"{${INIT_SLURM},${STABLE_SLURM},${FINAL_SLURM},${SINGLE_FINAL_SLURM},${OBS_SLURM}} 2>/dev/null || true

# ---------------------------------------------------------
# Write INIT slurm
# ---------------------------------------------------------
cat > "${MODEL_DIR}/${INIT_SLURM}" <<EOF
#!/bin/bash
#SBATCH --job-name=${JOB_PREFIX}_cmfgen_INIT${tag:+_${tag}}
#SBATCH --output=cmf_INIT_%j.out
#SBATCH --error=cmf_INIT_%j.err
#SBATCH --time=${INIT_TIME}
#SBATCH --partition=${INIT_PARTITION}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${cpus}
#SBATCH --mem=${init_mem}
#SBATCH --chdir=${MODEL_DIR}

# INIT: no #SBATCH --signal, no in-job graceful timer, no requeue (runs until exit or walltime).

set -euo pipefail

check_outgen_for_fatal() {
  # CMFGEN can continue for a while after NaNs appear, producing corrupted restart state.
  # Catch the common hard-failure patterns early so Slurm marks this stage as failed.
  local outgen="OUTGEN"
  if [[ ! -f "\$outgen" ]]; then
    echo "ERROR: missing \$outgen after batch.sh (CMFGEN did not write OUTGEN?)" >&2
    return 2
  fi
  if grep -Eq "NaN|Error in MON_INTERP|Segmentation fault|SIG(ABRT|SEGV)|floating invalid|Error in comp_j_blank\\.f|Mean intensity blowing up|Error in SET_CONT_FREQ_V4|MIN_FREQ too big|FREQ\\(1\\)\\s*=" "\$outgen"; then
    echo "ERROR: fatal pattern detected in \$outgen (NaN / MON_INTERP / crash). Aborting stage." >&2
    echo "---- OUTGEN tail (last 80 lines) ----" >&2
    tail -n 80 "\$outgen" >&2 || true
    echo "------------------------------------" >&2
    return 3
  fi
}

auto_set_rstar_lstar_from_hydro() {
  # Optional blackbody override AFTER prepare_VADAT.py (default OFF: prepare owns LSTAR + RSTAR/RMAX).
  # Set CMFGEN_AUTO_SET_LSTAR=1 and/or CMFGEN_AUTO_SET_RSTAR=1 to restore blackbody tuning.
  #
  # Controls (env vars):
  # - CMFGEN_PHOT_TAU: target optical depth for "photosphere" (default 0.6666667)
  # - CMFGEN_BB_TEFF_K: blackbody Teff in K (default 9100)
  # - CMFGEN_LSTAR_RADIUS_MODE: RMAX or PHOT (default RMAX)
  # - CMFGEN_AUTO_SET_RSTAR: 1/0 (default 0)
  # - CMFGEN_AUTO_SET_LSTAR: 1/0 (default 0)
  python3 - <<'PY'
from __future__ import annotations

import math
import os
import re
from pathlib import Path

hydro = Path("SN_HYDRO_DATA")
vadat = Path("VADAT")
if not hydro.is_file():
    raise SystemExit(0)

auto_rstar = os.environ.get("CMFGEN_AUTO_SET_RSTAR", "0").strip() not in {"0", "F", "FALSE", "false"}
auto_lstar = os.environ.get("CMFGEN_AUTO_SET_LSTAR", "0").strip() not in {"0", "F", "FALSE", "false"}
tau_phot = float(os.environ.get("CMFGEN_PHOT_TAU", "0.6666667"))
teff = float(os.environ.get("CMFGEN_BB_TEFF_K", "9100"))
lstar_rmode = (os.environ.get("CMFGEN_LSTAR_RADIUS_MODE", "RMAX") or "RMAX").strip().upper()
if tau_phot <= 0:
    raise SystemExit("auto-RSTAR/LSTAR: CMFGEN_PHOT_TAU must be > 0")
if teff <= 0:
    raise SystemExit("auto-RSTAR/LSTAR: CMFGEN_BB_TEFF_K must be > 0")
if lstar_rmode not in {"RMAX", "PHOT"}:
    raise SystemExit(f"auto-RSTAR/LSTAR: invalid CMFGEN_LSTAR_RADIUS_MODE={lstar_rmode!r} (use RMAX or PHOT)")

lines = hydro.read_text().splitlines()

def _read_block(marker: str) -> list[float]:
    try:
        i0 = lines.index(marker)
    except ValueError as e:
        raise ValueError(f"missing marker {marker!r} in {hydro}") from e
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

# SN_HYDRO_DATA is outer->inner (radius decreasing).
R_1e10 = _read_block("Radius grid (10^10cm)")
rho = _read_block("Density (gm/cm^3)")
kappa = _read_block("Kappa (cm^2/gm)")
if not (len(R_1e10) == len(rho) == len(kappa)):
    raise SystemExit(
        f"auto-RSTAR/LSTAR: block length mismatch: R={len(R_1e10)} rho={len(rho)} kappa={len(kappa)}"
    )
n = len(R_1e10)
if n < 3:
    raise SystemExit("auto-RSTAR/LSTAR: too few zones")

# Compute cumulative tau inward with tau[0]=0 at outer boundary.
tau = [0.0] * n
for i in range(1, n):
    r0 = R_1e10[i - 1]
    r1 = R_1e10[i]
    dr_cm = (r0 - r1) * 1.0e10
    if dr_cm <= 0:
        continue
    chi0 = kappa[i - 1] * rho[i - 1]
    chi1 = kappa[i] * rho[i]
    tau[i] = tau[i - 1] + 0.5 * (chi0 + chi1) * dr_cm

# Find photosphere radius by interpolation where tau crosses tau_phot.
rphot_1e10 = R_1e10[-1]
if tau[-1] <= tau_phot:
    # Entire grid optically thin; fall back to inner radius.
    rphot_1e10 = R_1e10[-1]
else:
    j = 0
    while j < n and tau[j] < tau_phot:
        j += 1
    j = min(max(j, 1), n - 1)
    t0, t1 = tau[j - 1], tau[j]
    r0, r1 = R_1e10[j - 1], R_1e10[j]
    if t1 <= t0:
        rphot_1e10 = r1
    else:
        f = (tau_phot - t0) / (t1 - t0)
        rphot_1e10 = r0 + f * (r1 - r0)

sigma = 5.670374419e-5  # erg / cm^2 / s / K^4
Lsun = 3.828e33
rmax_1e10 = max(R_1e10)
r_for_lstar_1e10 = (rmax_1e10 if lstar_rmode == "RMAX" else rphot_1e10)
R_cm = max(r_for_lstar_1e10, 1.0e-99) * 1.0e10
Lbb_lsun = max((4.0 * math.pi * R_cm * R_cm * sigma * teff**4) / Lsun, 1.0e-30)

def _update_header_lstar() -> None:
    for i, ln in enumerate(lines[:40]):
        if ln.strip().startswith("Energy Lum [lsun]:"):
            # Always overwrite (unless CMFGEN_AUTO_SET_LSTAR=0).
            lines[i] = f"Energy Lum [lsun]:                   {Lbb_lsun:.6E}"
            hydro.write_text("\n".join(lines) + "\n")
            print(
                "auto-LSTAR: wrote SN_HYDRO_DATA Energy Lum [lsun]={:.6e} "
                "(Teff={:.1f}K, mode={}, R_for_LSTAR={:.6e} 1e10cm, Rphot={:.6e} 1e10cm, tau={})".format(
                    Lbb_lsun, teff, lstar_rmode, r_for_lstar_1e10, rphot_1e10, tau_phot
                )
            )
            return

def _patch_vadat_scalar(tag: str, value: float) -> None:
    if not vadat.is_file():
        return
    vlines = vadat.read_text().splitlines()
    pat = re.compile(rf"^\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[DEde][+-]?\d+)?)\s+\[{re.escape(tag)}\]\s*$")
    for i, ln in enumerate(vlines):
        if f"[{tag}]" in ln:
            vlines[i] = f"{value:.8g}  [{tag}]"
            vadat.write_text("\n".join(vlines) + "\n")
            print(f"auto-{tag}: set VADAT {tag}={value}")
            return

if auto_lstar:
    _update_header_lstar()
    _patch_vadat_scalar("LSTAR", Lbb_lsun)

if auto_rstar:
    _patch_vadat_scalar("RSTAR", rphot_1e10)
    # Keep outer radius fixed in physical units when RSTAR changes.
    rmax_ratio = max(rmax_1e10 / max(rphot_1e10, 1.0e-99), 1.0)
    _patch_vadat_scalar("RMAX", rmax_ratio)

PY
}

echo "======================================"
echo "CMFGEN INIT Job Started"
echo "Job ID: \$SLURM_JOB_ID"
echo "Partition: \$SLURM_JOB_PARTITION"
echo "CPUs: \$SLURM_CPUS_PER_TASK"
echo "Node: \$SLURMD_NODENAME"
echo "======================================"

module purge
module load intel/2024.0
module load mkl/2026.0
module load python-waterboa/2024.06
# Explicitly set MKL library path (workaround for module loading issues on compute nodes)
export LD_LIBRARY_PATH=/viper/u1/system/soft/RHEL_9/packages/x86_64/intel_oneapi/2026.0/mkl/2026.0/lib/intel64:\${LD_LIBRARY_PATH:+:\$LD_LIBRARY_PATH}

# ---------------------------
# Resource settings (from SLURM)
# ---------------------------
CPUS=\${SLURM_CPUS_PER_TASK:-2}
export OMP_NUM_THREADS=\$CPUS
export MKL_NUM_THREADS=\$CPUS
export OMP_STACKSIZE=512M
export KMP_STACKSIZE=512M

echo "CPUs: \$CPUS"

export cmfdist="${CMFDIST}"
export ATOMIC="${ATOMIC}"

# Match a fresh model: clear partial run trash (OUTGEN, POINT*, *_IN, etc.) before prep.
command -v csh >/dev/null 2>&1 || { echo "ERROR: csh not in PATH (restart_model.sh)" >&2; exit 1; }
echo "Hard reset for INIT (restart_model.sh --all)..."
csh "\$cmfdist/com/restart_model.sh" --all || { echo "ERROR: restart_model.sh --all failed" >&2; exit 1; }

# Clean old CMFGEN output files for fresh INIT run
echo "Cleaning old CMFGEN files..."
if [[ -f "${CMFDIST}/com/clean.sh" ]]; then
  bash "${CMFDIST}/com/clean.sh" || echo "WARNING: clean.sh failed (non-fatal)"
else
  echo "WARNING: clean.sh not found at ${CMFDIST}/com/clean.sh"
fi

# INIT has no Slurm requeue: always run Python prep for this submission so a stale
# .CMFGEN_INIT_PREPPED (e.g. left after OOM during batch.sh) cannot skip OUTGEN wipe.
echo "Running Python prep scripts (INIT chunk)..."
rm -f .CMFGEN_INIT_PREPPED
rm -f batch.log OUTGEN
python "${VADAT_SCRIPT}" --init --convergence-mode stable --fix_ne --vel-law "${VEL_LAW}" ${varyt_arg} --vadat-overrides "${CONFIG_FILE}" > 'batch.log'
python "${INITS_SCRIPT}" --init --num-its ${INIT_NUM_ITS} ${varyt_arg} >> 'batch.log'
export CMFGEN_AUTO_SET_RSTAR="\${CMFGEN_AUTO_SET_RSTAR:-0}"
export CMFGEN_AUTO_SET_LSTAR="\${CMFGEN_AUTO_SET_LSTAR:-0}"
auto_set_rstar_lstar_from_hydro

echo "Running out2in.sh and batch.sh..."
"\$cmfdist/com/out2in.sh"


tcsh -f ./batch.sh
check_outgen_for_fatal

echo "Preparing restart inputs (out2in.sh)..."
"\$cmfdist/com/out2in.sh"

touch .CMFGEN_INIT_PREPPED
echo "INIT Job Complete"
echo "======================================"
EOF

# ---------------------------------------------------------
# Write STABLE slurm (pure Lambda settling between init and final)
# ---------------------------------------------------------
cat > "${MODEL_DIR}/${STABLE_SLURM}" <<EOF
#!/bin/bash
#SBATCH --job-name=${JOB_PREFIX}_cmfgen_STABLE${tag:+_${tag}}
#SBATCH --output=cmf_STABLE_%j.out
#SBATCH --error=cmf_STABLE_%j.err
#SBATCH --time=${STABLE_TIME}
#SBATCH --partition=${STABLE_PARTITION}
#SBATCH --signal=B:USR1@${STABLE_SIGNAL_BEFORE_SEC}
#SBATCH --requeue
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${cpus}
#SBATCH --mem=${mem}
#SBATCH --chdir=${MODEL_DIR}

set -euo pipefail

check_outgen_for_fatal() {
  local outgen="OUTGEN"
  if [[ ! -f "\$outgen" ]]; then
    echo "ERROR: missing \$outgen after batch.sh (CMFGEN did not write OUTGEN?)" >&2
    return 2
  fi
  if grep -Eq "NaN|Error in MON_INTERP|Segmentation fault|SIG(ABRT|SEGV)|floating invalid|Error in comp_j_blank\\.f|Mean intensity blowing up|Error in SET_CONT_FREQ_V4|MIN_FREQ too big|FREQ\\(1\\)\\s*=" "\$outgen"; then
    echo "ERROR: fatal pattern detected in \$outgen (NaN / MON_INTERP / crash). Aborting stage." >&2
    echo "---- OUTGEN tail (last 80 lines) ----" >&2
    tail -n 80 "\$outgen" >&2 || true
    echo "------------------------------------" >&2
    return 3
  fi
}

echo "======================================"
echo "CMFGEN STABLE Job Started"
echo "Job ID: \$SLURM_JOB_ID"
echo "Partition: \$SLURM_JOB_PARTITION"
echo "CPUs: \$SLURM_CPUS_PER_TASK"
echo "Node: \$SLURMD_NODENAME"
echo "======================================"

module purge
module load intel/2024.0
module load mkl/2026.0
module load python-waterboa/2024.06
# Explicitly set MKL library path (workaround for module loading issues on compute nodes)
export LD_LIBRARY_PATH=/viper/u1/system/soft/RHEL_9/packages/x86_64/intel_oneapi/2026.0/mkl/2026.0/lib/intel64:\${LD_LIBRARY_PATH:+:\$LD_LIBRARY_PATH}

# ---------------------------
# Resource settings (from SLURM)
# ---------------------------
CPUS=\${SLURM_CPUS_PER_TASK:-2}
export OMP_NUM_THREADS=\$CPUS
export MKL_NUM_THREADS=\$CPUS
export OMP_STACKSIZE=512M
export KMP_STACKSIZE=512M

echo "CPUs: \$CPUS"

# Label this STABLE requeue instance for easier tracking.
restart_count="\${SLURM_RESTART_COUNT:-0}"
requeue_idx=\$((restart_count + 1))
stable_job_name="${JOB_PREFIX}_STABLE_N\${requeue_idx}_cmfgen_${tag}"
if command -v scontrol >/dev/null 2>&1; then
  scontrol update JobId="\$SLURM_JOB_ID" JobName="\$stable_job_name" >/dev/null 2>&1 || true
fi

export cmfdist="${CMFDIST}"
export ATOMIC="${ATOMIC}"

# Graceful stop for STABLE: before walltime, set NUM_ITS=0 and requeue.
request_graceful_stop() {
  echo ""
  echo "======================================================"
  echo "Pre-timeout signal received: requesting graceful stop."
  echo "Setting NUM_ITS=0 in IN_ITS so CMFGEN exits cleanly."
  echo "======================================================"
  python - <<'PY'
from pathlib import Path
import re
p = Path("IN_ITS")
if not p.exists():
    exit(0)
lines = p.read_text().splitlines()
out = []
for ln in lines:
    if re.match(r"\s*\d+\s+\[NUM_ITS\]", ln):
        out.insert(0, "0 [NUM_ITS] ")
    else:
        out.append(ln)
if not any("[NUM_ITS]" in ln for ln in out):
    out.insert(0, "0 [NUM_ITS] ")
p.write_text("\n".join(out) + "\n")
PY
  touch .CMFGEN_STABLE_REQUEST_REQUEUE
}
trap request_graceful_stop USR1

# Fallback timer: request graceful stop ~22h after start (for 24h walltime).
graceful_timer_pid=""
graceful_after_sec="\${CMFGEN_STABLE_GRACEFUL_AFTER_SEC:-${STABLE_GRACEFUL_SLEEP_SEC}}"
(
  sleep "\$graceful_after_sec"
  if command -v scontrol >/dev/null 2>&1; then
    scontrol signal JobId="\$SLURM_JOB_ID" USR1 2>/dev/null || true
  fi
) &
graceful_timer_pid="\$!"

# New Slurm allocation (not a requeue continuation): drop stale prep marker so Python
# prep can re-run after rsync/OOM. Requeued STABLE chunks keep SLURM_RESTART_COUNT>0.
if [[ "\${SLURM_RESTART_COUNT:-0}" -eq 0 ]]; then
  rm -f .CMFGEN_STABLE_PREPPED
fi

if [[ ! -f .CMFGEN_STABLE_PREPPED ]]; then
  echo "Running out2in.sh (get populations from INIT)..."
  "\$cmfdist/com/out2in.sh"

  echo "Running Python prep scripts for STABLE phase..."
  python "${VADAT_SCRIPT}" --stable --convergence-mode stable --vel-law "${VEL_LAW}" ${varyt_arg} --vadat-overrides "${CONFIG_FILE}" > 'batch.log'
  python "${INITS_SCRIPT}" --stable --num-its ${STABLE_NUM_ITS} ${varyt_arg} >> 'batch.log'

  # Outer-boundary refinement to reduce J-instabilities (comp_j_blank blow-ups).
  # Allow override via CMFGEN_N_OB_INS, default 30.
  CMFGEN_N_OB_INS="\${CMFGEN_N_OB_INS:-30}"
  python - <<PY
from pathlib import Path

p = Path("VADAT")
if not p.exists():
    raise SystemExit(0)

desired = int("\${CMFGEN_N_OB_INS}")
lines = p.read_text().splitlines()
out = []
done = False
for ln in lines:
    if "[N_OB_INS]" in ln:
        out.append(f"{desired} [N_OB_INS]")
        done = True
    else:
        out.append(ln)
if not done:
    inserted = False
    for i, ln in enumerate(out):
        if "[N_IB_INS]" in ln:
            out.insert(i + 1, f"{desired} [N_OB_INS]")
            inserted = True
            break
    if not inserted:
        out.insert(0, f"{desired} [N_OB_INS]")
p.write_text("\n".join(out) + "\n")
print(f"Applied outer-boundary refinement: N_OB_INS={desired}")
PY

  touch .CMFGEN_STABLE_PREPPED
else
  echo "Skipping Python prep scripts (already prepped)."
fi

echo "Running batch.sh (pure Lambda settling)..."
tcsh -f ./batch.sh
check_outgen_for_fatal

echo "Preparing restart inputs (out2in.sh)..."
"\$cmfdist/com/out2in.sh"

# Requeue STABLE if it didn't converge (hit walltime)
if [[ -f .CMFGEN_STABLE_REQUEST_REQUEUE ]]; then
  echo "Requeue requested; requeueing job \$SLURM_JOB_ID..."
  rm -f .CMFGEN_STABLE_REQUEST_REQUEUE
  scontrol requeue "\$SLURM_JOB_ID" || echo "WARNING: scontrol requeue failed; re-submit manually."
  exit 0
fi

echo "STABLE Job Complete"
echo "======================================"
EOF

# ---------------------------------------------------------
# Write FINAL slurm (restore+delete)
# ---------------------------------------------------------
cat > "${MODEL_DIR}/${FINAL_SLURM}" <<EOF
#!/bin/bash
#SBATCH --job-name=${JOB_PREFIX}_cmfgen_FINAL${tag:+_${tag}}
#SBATCH --output=cmf_FINAL_%j.out
#SBATCH --error=cmf_FINAL_%j.err
#SBATCH --time=${timelimit}
#SBATCH --partition=${part}
$([ "${FINAL_WALL_GRACEFUL}" != "0" ] && echo "#SBATCH --signal=B:USR1@${FINAL_SIGNAL_BEFORE_SEC}")
#SBATCH --requeue
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${cpus}
#SBATCH --mem=${mem}
# #SBATCH --mail-type=END
# #SBATCH --mail-user=sbrennan@mpe.mpg.de
#SBATCH --chdir=${MODEL_DIR}

set -euo pipefail

check_outgen_for_fatal() {
  local outgen="OUTGEN"
  if [[ ! -f "\$outgen" ]]; then
    echo "ERROR: missing \$outgen after batch.sh (CMFGEN did not write OUTGEN?)" >&2
    return 2
  fi
  if grep -Eq "NaN|Error in MON_INTERP|Segmentation fault|SIG(ABRT|SEGV)|floating invalid|Error in comp_j_blank\\.f|Mean intensity blowing up|bad initial population guesses|Predicted changes are too large" "\$outgen"; then
    echo "ERROR: fatal pattern detected in \$outgen (NaN / MON_INTERP / crash). Aborting stage." >&2
    echo "---- OUTGEN tail (last 80 lines) ----" >&2
    tail -n 80 "\$outgen" >&2 || true
    echo "------------------------------------" >&2
    return 3
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

auto_rstar = os.environ.get("CMFGEN_AUTO_SET_RSTAR", "0").strip() not in {"0", "F", "FALSE", "false"}
auto_lstar = os.environ.get("CMFGEN_AUTO_SET_LSTAR", "0").strip() not in {"0", "F", "FALSE", "false"}
tau_phot = float(os.environ.get("CMFGEN_PHOT_TAU", "0.6666667"))
teff = float(os.environ.get("CMFGEN_BB_TEFF_K", "9100"))
lstar_rmode = (os.environ.get("CMFGEN_LSTAR_RADIUS_MODE", "RMAX") or "RMAX").strip().upper()
if tau_phot <= 0:
    raise SystemExit("auto-RSTAR/LSTAR: CMFGEN_PHOT_TAU must be > 0")
if teff <= 0:
    raise SystemExit("auto-RSTAR/LSTAR: CMFGEN_BB_TEFF_K must be > 0")
if lstar_rmode not in {"RMAX", "PHOT"}:
    raise SystemExit(f"auto-RSTAR/LSTAR: invalid CMFGEN_LSTAR_RADIUS_MODE={lstar_rmode!r} (use RMAX or PHOT)")

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
final_job_name="${JOB_PREFIX}_FINAL_N\${requeue_idx}_cmfgen"
if command -v scontrol >/dev/null 2>&1; then
  scontrol update JobId="\$SLURM_JOB_ID" JobName="\$final_job_name" >/dev/null 2>&1 || true
fi
echo "FINAL iteration label (for job-name tracking): \${final_job_name}"

# Do not delete Slurm stdout/stderr logs; they are essential for debugging failures.
module purge
module load intel/2024.0
module load mkl/2026.0
module load python-waterboa/2024.06
# Explicitly set MKL library path (workaround for module loading issues on compute nodes)
export LD_LIBRARY_PATH=/viper/u1/system/soft/RHEL_9/packages/x86_64/intel_oneapi/2026.0/mkl/2026.0/lib/intel64:\${LD_LIBRARY_PATH:+:\$LD_LIBRARY_PATH}

export cmfdist="${CMFDIST}"
export ATOMIC="${ATOMIC}"

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
CPUS=\${SLURM_CPUS_PER_TASK:-2}
export OMP_NUM_THREADS=\$CPUS
export MKL_NUM_THREADS=\$CPUS
export OMP_STACKSIZE=512M
export KMP_STACKSIZE=512M

echo "CPUs: \$CPUS"

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
python "${VADAT_SCRIPT}" --final --convergence-mode stable --vel-law "${VEL_LAW}" ${varyt_arg} --vadat-overrides "${CONFIG_FILE}" > 'batch.log'
python "${INITS_SCRIPT}" --final --num-its ${FINAL_NUM_ITS} ${varyt_arg} >> 'batch.log'
export CMFGEN_AUTO_SET_RSTAR="\${CMFGEN_AUTO_SET_RSTAR:-0}"
export CMFGEN_AUTO_SET_LSTAR="\${CMFGEN_AUTO_SET_LSTAR:-0}"
auto_set_rstar_lstar_from_hydro

# Outer-boundary refinement to reduce J-instabilities (comp_j_blank blow-ups).
# Allow override via CMFGEN_N_OB_INS, default 30.
CMFGEN_N_OB_INS="\${CMFGEN_N_OB_INS:-30}"
python - <<PY
from pathlib import Path

p = Path("VADAT")
if not p.exists():
    raise SystemExit(0)

desired = int("\${CMFGEN_N_OB_INS}")
lines = p.read_text().splitlines()
out = []
done = False
for ln in lines:
    if "[N_OB_INS]" in ln:
        out.append(f"{desired} [N_OB_INS]")
        done = True
    else:
        out.append(ln)
if not done:
    inserted = False
    for i, ln in enumerate(out):
        if "[N_IB_INS]" in ln:
            out.insert(i + 1, f"{desired} [N_OB_INS]")
            inserted = True
            break
    if not inserted:
        out.insert(0, f"{desired} [N_OB_INS]")
p.write_text("\n".join(out) + "\n")
print(f"Applied outer-boundary refinement: N_OB_INS={desired}")
PY

# Prevent VAR_MOM_J_CMF NV<ND crash (observed ND=218, NV=200) by capping
# *insertion depth count* after FINAL prep for this model.
# IMPORTANT: Do not cap outer-boundary resolution (N_OB_INS) here; CMFGEN may need
# a finer outer boundary to avoid J-instability (comp_j_blank blow-ups).
python - <<'PY'
from pathlib import Path
p = Path("VADAT")
lines = p.read_text().splitlines() if p.exists() else []
out = []
for ln in lines:
    if "[N_INS]" in ln:
        out.append("1 [N_INS]")
    else:
        out.append(ln)
p.write_text("\n".join(out) + "\n")
print("Applied SN13 FINAL grid cap: N_INS=1 (kept N_OB_INS/N_IB_INS from prepare_VADAT)")
PY

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
require_restart_inputs_or_abort
tcsh -f ./batch.sh
check_outgen_for_fatal
stop_graceful_timer

echo "Preparing restart inputs (out2in.sh)..."
"\$cmfdist/com/out2in.sh"

if [[ -f .CMFGEN_REQUEST_REQUEUE ]]; then
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
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${cpus}
#SBATCH --mem=${mem}
#SBATCH --chdir=${MODEL_DIR}/obs

set -euo pipefail

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
export atomic="${ATOMIC}"

echo "Preparing batobs.sh..."
cp -f "${OBS_TEMPLATE_DIR}/batobs.sh" batobs.sh
cp -f "${OBS_TEMPLATE_DIR}/CMF_FLUX_PARAM_INIT" CMF_FLUX_PARAM_INIT

chmod +x batobs.sh

# ---------------------------
# Resource settings (from SLURM)
# ---------------------------
CPUS=\${SLURM_CPUS_PER_TASK:-2}
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

# Plot observed spectrum after BATOBS finishes (non-fatal if plotting fails).
echo "Plotting observed spectrum (post-BATOBS)..."
plotter="${PROJECT_ROOT}/scripts/plot_cmfgen_obs_spectrum.py"
if [[ -f "\$plotter" ]]; then
  if command -v python3 >/dev/null 2>&1; then
    python3 "\$plotter" "obs_fin_10" --out "obs_spectrum.png" || echo "WARNING: spectrum plot failed (non-fatal)"
  else
    python "\$plotter" "obs_fin_10" --out "obs_spectrum.png" || echo "WARNING: spectrum plot failed (non-fatal)"
  fi
else
  echo "WARNING: plotter not found at \$plotter (skip spectrum plot)"
fi

echo "BATOBS Job Complete"
echo "======================================"
EOF

# ---------------------------------------------------------
# Write SINGLE_FINAL slurm (LTE cold start; no prior INIT)
# ---------------------------------------------------------
cat > "${MODEL_DIR}/${SINGLE_FINAL_SLURM}" <<EOF
#!/bin/bash
#SBATCH --job-name=${JOB_PREFIX}_cmfgen_FINAL${tag:+_${tag}}
#SBATCH --output=cmf_FINAL_%j.out
#SBATCH --error=cmf_FINAL_%j.err
#SBATCH --time=${timelimit}
#SBATCH --partition=${part}
$([ "${FINAL_WALL_GRACEFUL}" != "0" ] && echo "#SBATCH --signal=B:USR1@${FINAL_SIGNAL_BEFORE_SEC}")
#SBATCH --requeue
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${cpus}
#SBATCH --mem=${mem}
#SBATCH --chdir=${MODEL_DIR}

set -euo pipefail

check_outgen_for_fatal() {
  local outgen="OUTGEN"
  if [[ ! -f "\$outgen" ]]; then
    echo "ERROR: missing \$outgen after batch.sh (CMFGEN did not write OUTGEN?)" >&2
    return 2
  fi
  if grep -Eq "NaN|Error in MON_INTERP|Segmentation fault|SIG(ABRT|SEGV)|floating invalid|Error in comp_j_blank\\.f|Mean intensity blowing up|bad initial population guesses|Predicted changes are too large" "\$outgen"; then
    echo "ERROR: fatal pattern detected in \$outgen (NaN / MON_INTERP / crash). Aborting stage." >&2
    echo "---- OUTGEN tail (last 80 lines) ----" >&2
    tail -n 80 "\$outgen" >&2 || true
    echo "------------------------------------" >&2
    return 3
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
final_job_name="${JOB_PREFIX}_FINAL_N\${requeue_idx}_cmfgen"
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

CPUS=\${SLURM_CPUS_PER_TASK:-2}
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
  # Pass CMFGEN_INC_RAD_DECAYS environment variable to prepare_VADAT.py
  if [[ -n "\${CMFGEN_INC_RAD_DECAYS:-}" ]]; then
    echo "CMFGEN_INC_RAD_DECAYS=\${CMFGEN_INC_RAD_DECAYS} (passed to prepare_VADAT.py)"
  fi
  # Build optional flags for prepare_VADAT.py
  NO_REL_FLAG=""
  if [[ "${CMFGEN_NO_RELATIVISTIC:-0}" == "1" ]]; then
    NO_REL_FLAG="--no-relativistic"
    echo "CMFGEN_NO_RELATIVISTIC=1 (passed --no-relativistic to prepare_VADAT.py)"
  fi
  python "${VADAT_SCRIPT}" --single-final --convergence-mode stable --vel-law "${VEL_LAW}" ${varyt_arg} ${NO_REL_FLAG} --vadat-overrides "${CONFIG_FILE}" > 'batch.log'
  python "${INITS_SCRIPT}" --single-final --num-its ${FINAL_NUM_ITS} ${varyt_arg} >> 'batch.log'
  touch .CMFGEN_SINGLE_FINAL_PREPPED
else
  echo "Requeue \${requeue_idx}: restoring checkpoint state..."
  "\$cmfdist/com/out2in.sh"
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
tcsh -f ./batch.sh
check_outgen_for_fatal
stop_graceful_timer

echo "Preparing restart inputs (out2in.sh)..."
"\$cmfdist/com/out2in.sh"

if [[ -f .CMFGEN_REQUEST_REQUEUE ]]; then
  echo "Requeue requested; requeueing job \$SLURM_JOB_ID..."
  rm -f .CMFGEN_REQUEST_REQUEUE
  scontrol requeue "\$SLURM_JOB_ID" || echo "WARNING: scontrol requeue failed; re-submit manually."
  exit 0
fi

echo "SINGLE-FINAL Job Complete"
echo "======================================"
EOF

chmod +x "${MODEL_DIR}/${INIT_SLURM}" "${MODEL_DIR}/${STABLE_SLURM}" "${MODEL_DIR}/${FINAL_SLURM}" "${MODEL_DIR}/${SINGLE_FINAL_SLURM}" "${MODEL_DIR}/${OBS_SLURM}"

echo ""
echo "Generated Slurm scripts:"
echo "  ${MODEL_DIR}/${INIT_SLURM}"
echo "  ${MODEL_DIR}/${STABLE_SLURM}"
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
# Match INIT/STABLE/OBS (*_cmfgen_<tag>) and FINAL (*_FINAL_N*_cmfgen); STABLE was missing
# from the old pattern, and FINAL uses _FINAL_N1_cmfgen not _FINAL_cmfgen_<tag>.
# Also match dynamically renamed FINAL jobs (e.g., _FINAL_N1_cmfgen, _FINAL_N2_cmfgen)
existing=$(squeue -u "${USER:-$(whoami)}" -h -o "%j" 2>/dev/null | grep -E "^${JOB_PREFIX}_(INIT|STABLE|FINAL|OBS)_|${JOB_PREFIX}_FINAL_N[0-9]+_cmfgen" || true)
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
  echo "Restart mode: submitting FINAL and OBS only (no INIT)..."

  jid_final=$(sbatch --parsable $DEP_FLAG "${MODEL_DIR}/${FINAL_SLURM}")
  echo "FINAL Job ID:  $jid_final"

  jid_obs=$(sbatch --parsable --dependency=afterok:"$jid_final" "${MODEL_DIR}/${OBS_SLURM}")
  echo "BATOBS Job ID: $jid_obs (depends on $jid_final)"
elif [[ $single_final -eq 1 ]]; then
  echo "Single-final mode: submitting FINAL (LTE cold start) and OBS only (no INIT)..."

  jid_final=$(sbatch --parsable $DEP_FLAG "${MODEL_DIR}/${SINGLE_FINAL_SLURM}")
  echo "FINAL Job ID:  $jid_final"

  jid_obs=$(sbatch --parsable --dependency=afterok:"$jid_final" "${MODEL_DIR}/${OBS_SLURM}")
  echo "BATOBS Job ID: $jid_obs (depends on $jid_final)"
else
  echo "Submitting jobs with dependencies..."

  jid_init=$(sbatch --parsable $DEP_FLAG "${MODEL_DIR}/${INIT_SLURM}")
  echo "INIT Job ID:   $jid_init"

  if [[ $with_stable -eq 1 ]]; then
    jid_stable=$(sbatch --parsable --dependency=afterok:"$jid_init" "${MODEL_DIR}/${STABLE_SLURM}")
    echo "STABLE Job ID: $jid_stable (depends on $jid_init)"

    jid_final=$(sbatch --parsable --dependency=afterok:"$jid_stable" "${MODEL_DIR}/${FINAL_SLURM}")
    echo "FINAL Job ID:  $jid_final (depends on $jid_stable)"
  else
    echo "Skipping STABLE stage (default). Use --with-stable to include it."
    jid_final=$(sbatch --parsable --dependency=afterok:"$jid_init" "${MODEL_DIR}/${FINAL_SLURM}")
    echo "FINAL Job ID:  $jid_final (depends on $jid_init)"
  fi

  jid_obs=$(sbatch --parsable --dependency=afterok:"$jid_final" "${MODEL_DIR}/${OBS_SLURM}")
  echo "BATOBS Job ID: $jid_obs (depends on $jid_final)"
fi

echo ""
echo "======================================================"
echo "All jobs submitted successfully!"
echo "Monitor with: squeue -u \$USER"
echo "======================================================"