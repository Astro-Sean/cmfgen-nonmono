#!/bin/bash
# =============================================================================
# submit_40d_staged.sh — Submit the X_1_fix 40-day staged CMFGEN workflow
#
# Staged warm-restart sequence for physically-motivated convergence:
#   Stage 1 (s1_init):   LTE cold start, T+Ne fixed to STELLA hydro
#   Stage 2 (s2_final):  Warm restart, T varies via radiative equilibrium
#   Stage 3 (s3_nlte):   Full NLTE, T+Ne both vary [optional, --with-s3]
#   OBS:                 Formal solution (spectrum synthesis)
#
# Usage:
#   ./submit_40d_staged.sh              # Submit s1 -> s2 -> OBS
#   ./submit_40d_staged.sh --with-s3    # Submit s1 -> s2 -> s3 -> OBS
#   ./submit_40d_staged.sh --s2-only    # Skip s1 (use existing converged s1)
#   ./submit_40d_staged.sh --obs-only   # Just run OBS from existing RVTJ
#
# Naming convention:
#   SLURM jobs: X1fix_40d_s1_init, X1fix_40d_s2_final, X1fix_40d_s3_nlte, X1fix_40d_obs
#   Directory:  models/X_1_fix/cmfgen/40d/
#   Logs:       s1_init.log, s2_final.log, s3_nlte.log, obs.log
# =============================================================================

set -euo pipefail

PROJECT_DIR="/viper/ptmp1/sbrennan/SN2024pba"
MODEL_DIR="${PROJECT_DIR}/models/X_1_fix/cmfgen/40d"
CONFIGS="${PROJECT_DIR}/configs"
SCRIPTS="${PROJECT_DIR}/scripts"
SUBMIT_SCRIPT="${MODEL_DIR}/submit_cmfgen.sh"

# Stage selection
WITH_S3=false
S2_ONLY=false
OBS_ONLY=false

for arg in "$@"; do
    case "$arg" in
        --with-s3)   WITH_S3=true ;;
        --s2-only)   S2_ONLY=true ;;
        --obs-only)  OBS_ONLY=true ;;
        *) echo "Unknown arg: $arg"; exit 1 ;;
    esac
done

# Job name prefix
PREFIX="X1fix_40d"

echo "======================================================"
echo "X_1_fix 40-day Staged CMFGEN Workflow"
echo "  Directory: ${MODEL_DIR}"
echo "  Prefix:    ${PREFIX}"
echo "  Stages:    $([ "$OBS_ONLY" = true ] && echo "OBS only" || { [ "$S2_ONLY" = true ] && echo "s2 -> OBS" || { [ "$WITH_S3" = true ] && echo "s1 -> s2 -> s3 -> OBS" || echo "s1 -> s2 -> OBS"; } })"
echo "======================================================"

# -----------------------------------------------------------------------------
# Helper: run out2in (rename *OUT -> *_IN for warm restart)
# -----------------------------------------------------------------------------
out2in() {
    local dir="$1"
    local count=0
    for f in "$dir"/*OUT "$dir"/*_OUT; do
        [ -f "$f" ] || continue
        local base="${f%OUT}"
        base="${base%_}"  # strip trailing _ from *_OUT
        local dest="${base}_IN"
        mv "$f" "$dest"
        count=$((count + 1))
    done
    echo "  out2in: renamed ${count} files (*OUT -> *_IN)"
}

# -----------------------------------------------------------------------------
# Helper: prepare VADAT for a stage
# -----------------------------------------------------------------------------
prepare_vadat() {
    local stage="$1"
    local profile="$2"
    local config="$3"

    echo "  Preparing VADAT for ${stage} (--${profile})..."
    cd "${MODEL_DIR}"
    python3 "${SCRIPTS}/prepare_VADAT.py" \
        --${profile} \
        --vadat-overrides "${config}" \
        --model-spec MODEL_SPEC \
        --inc-rad-decays T \
        --vel-law 12 \
        2>&1 | tail -5
}

# -----------------------------------------------------------------------------
# Helper: submit a SLURM job
# -----------------------------------------------------------------------------
submit_job() {
    local stage_name="$1"   # e.g. s1_init
    local stage_dir="$2"    # working directory
    local cmd="$3"          # command to run
    local dep="$4"          # dependency job ID (or empty)

    local job_name="${PREFIX}_${stage_name}"
    local slurm_script="${stage_dir}/run_${stage_name}.slurm"
    local log_file="${stage_dir}/${stage_name}.log"

    cat > "$slurm_script" << EOF
#!/bin/bash
#SBATCH --job-name=${job_name}
#SBATCH --partition=general
#SBATCH --time=23:59:00
#SBATCH --mem=16G
#SBATCH --output=${log_file}
#SBATCH --error=${log_file}

cd "${stage_dir}"
${cmd}
EOF

    local dep_flag=""
    if [ -n "$dep" ]; then
        dep_flag="--dependency=afterok:${dep}"
    fi

    local jid
    jid=$(sbatch --parsable $dep_flag "$slurm_script")
    echo "  Submitted ${job_name} (JobID: ${jid})"
    echo "$jid"
}

# -----------------------------------------------------------------------------
# OBS job (always last)
# -----------------------------------------------------------------------------
submit_obs() {
    local dep="$1"
    local obs_dir="${MODEL_DIR}/obs"
    local job_name="${PREFIX}_obs"
    local slurm_script="${obs_dir}/run_obs.slurm"
    local log_file="${MODEL_DIR}/obs.log"

    # Ensure OBS reads the correct RVTJ
    cat > "$slurm_script" << EOF
#!/bin/bash
#SBATCH --job-name=${job_name}
#SBATCH --partition=general
#SBATCH --time=23:59:00
#SBATCH --mem=16G
#SBATCH --output=${log_file}
#SBATCH --error=${log_file}

cd "${obs_dir}"
# Ensure RVTJ exists from the final converged stage
if [ ! -f "../RVTJ" ]; then
    echo "ERROR: ../RVTJ not found — previous stage must converge first"
    exit 1
fi

# Run OBS
./batobs.sh
EOF

    local dep_flag=""
    if [ -n "$dep" ]; then
        dep_flag="--dependency=afterok:${dep}"
    fi

    local jid
    jid=$(sbatch --parsable $dep_flag "$slurm_script")
    echo "  Submitted ${job_name} (JobID: ${jid})"
    echo "$jid"
}

# -----------------------------------------------------------------------------
# Main workflow
# -----------------------------------------------------------------------------

LAST_JOB_ID=""

if [ "$OBS_ONLY" = true ]; then
    # Just submit OBS
    echo ""
    echo "[OBS] Submitting OBS job only"
    submit_obs ""
    echo ""
    echo "Done. Monitor with: squeue --me"
    exit 0
fi

if [ "$S2_ONLY" = false ]; then
    # ---------------------------------------------------------------------------
    # Stage 1: s1_init (LTE cold start, T+Ne fixed)
    # ---------------------------------------------------------------------------
    echo ""
    echo "[Stage 1] s1_init: LTE cold start (FIX_T=T, FIX_NE=T)"
    prepare_vadat "s1_init" "single-final" "${CONFIGS}/X1fix_40d_s1.yaml"

    S1_CMD="export cmfdist=\$(dirname \$(readlink -f batch.sh | xargs dirname)/..)/cmfgen_source && export ATOMIC=\${cmfdist}/AtomicData && ./batch.sh"
    # Actually, use the batch.sh which already handles environment
    S1_CMD="./batch.sh"
    LAST_JOB_ID=$(submit_job "s1_init" "${MODEL_DIR}" "${S1_CMD}" "")
fi

# ---------------------------------------------------------------------------
# Stage 2: s2_final (warm restart, T varies via radiative equilibrium)
# ---------------------------------------------------------------------------
echo ""
echo "[Stage 2] s2_final: Warm restart (FIX_T=F, DO_T_AUTO=T, FIX_NE=T)"

# Create a wrapper script that does out2in + prepare VADAT + run
S2_WRAPPER="${MODEL_DIR}/run_s2_final.sh"
cat > "$S2_WRAPPER" << 'SCRIPT'
#!/bin/bash
set -e
cd "$(dirname "$0")"

echo "[s2_final] Step 1: out2in (rename *OUT -> *_IN)"
for f in *OUT *_OUT; do
    [ -f "$f" ] || continue
    base="${f%OUT}"
    base="${base%_}"
    mv "$f" "${base}_IN"
    echo "  $f -> ${base}_IN"
done

echo "[s2_final] Step 2: Prepare VADAT (final profile, FIX_T=F)"
python3 /viper/ptmp1/sbrennan/SN2024pba/scripts/prepare_VADAT.py \
    --final \
    --vadat-overrides /viper/ptmp1/sbrennan/SN2024pba/configs/X1fix_40d_s2.yaml \
    --model-spec MODEL_SPEC \
    --inc-rad-decays T \
    --vel-law 12 \
    2>&1 | tail -10

echo "[s2_final] Step 3: Run CMFGEN"
./batch.sh
SCRIPT
chmod +x "$S2_WRAPPER"

S2_CMD="./run_s2_final.sh"
LAST_JOB_ID=$(submit_job "s2_final" "${MODEL_DIR}" "${S2_CMD}" "${LAST_JOB_ID}")

# ---------------------------------------------------------------------------
# Stage 3: s3_nlte (optional - full NLTE, T+Ne vary)
# ---------------------------------------------------------------------------
if [ "$WITH_S3" = true ]; then
    echo ""
    echo "[Stage 3] s3_nlte: Full NLTE (FIX_T=F, FIX_NE=F)"

    S3_WRAPPER="${MODEL_DIR}/run_s3_nlte.sh"
    cat > "$S3_WRAPPER" << 'SCRIPT'
#!/bin/bash
set -e
cd "$(dirname "$0")"

echo "[s3_nlte] Step 1: out2in (rename *OUT -> *_IN)"
for f in *OUT *_OUT; do
    [ -f "$f" ] || continue
    base="${f%OUT}"
    base="${base%_}"
    mv "$f" "${base}_IN"
    echo "  $f -> ${base}_IN"
done

echo "[s3_nlte] Step 2: Prepare VADAT (final profile, FIX_T=F, FIX_NE=F)"
python3 /viper/ptmp1/sbrennan/SN2024pba/scripts/prepare_VADAT.py \
    --final \
    --vadat-overrides /viper/ptmp1/sbrennan/SN2024pba/configs/X1fix_40d_s3.yaml \
    --model-spec MODEL_SPEC \
    --inc-rad-decays T \
    --vel-law 12 \
    2>&1 | tail -10

echo "[s3_nlte] Step 3: Run CMFGEN"
./batch.sh
SCRIPT
    chmod +x "$S3_WRAPPER"

    S3_CMD="./run_s3_nlte.sh"
    LAST_JOB_ID=$(submit_job "s3_nlte" "${MODEL_DIR}" "${S3_CMD}" "${LAST_JOB_ID}")
fi

# ---------------------------------------------------------------------------
# OBS: Formal solution (spectrum synthesis)
# ---------------------------------------------------------------------------
echo ""
echo "[OBS] Submitting OBS job"
submit_obs "${LAST_JOB_ID}"

echo ""
echo "======================================================"
echo "Workflow submitted. Monitor with:"
echo "  squeue --me"
echo ""
echo "Job names:"
echo "  ${PREFIX}_s1_init   — Stage 1: LTE cold start"
echo "  ${PREFIX}_s2_final  — Stage 2: Warm restart (T varies)"
if [ "$WITH_S3" = true ]; then
    echo "  ${PREFIX}_s3_nlte   — Stage 3: Full NLTE (T+Ne vary)"
fi
echo "  ${PREFIX}_obs       — OBS: Spectrum synthesis"
echo "======================================================"
