#!/bin/bash
# [DEPRECATED] Master workflow script for single SN2024pba model run (HERACLES-based)
# DEPRECATED (2026-07-09): HERACLES does not fit the science case for this project
# (interacting Type IIn supernovae). The project now uses the STELLA -> CMFGEN
# pipeline exclusively. For the active workflow, see run_model_grid.py.
# Usage: ./scripts/run_single_model.sh configs/model_params.yaml

set -euo pipefail

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

# Parse command line arguments
if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <config_file.yaml>"
    echo "Example: $0 configs/model_params.yaml"
    exit 1
fi

CONFIG_FILE="$1"
CONFIG_FILE_ABS="$(cd "$(dirname "$CONFIG_FILE")" && pwd)/$(basename "$CONFIG_FILE")"

# Check if config file exists
if [[ ! -f "$CONFIG_FILE_ABS" ]]; then
    echo "ERROR: Config file not found: $CONFIG_FILE_ABS"
    exit 1
fi

# Function to read YAML value (simple implementation)
# For production, consider using a proper YAML parser
read_yaml() {
    local key="$1"
    local file="$2"
    grep "^${key}:" "$file" | awk '{print $2}' | sed 's/"//g' | sed "s/'//g"
}

# Read configuration
MODEL_NAME=$(read_yaml "model_name" "$CONFIG_FILE_ABS")
DESCRIPTION=$(read_yaml "description" "$CONFIG_FILE_ABS")

# Create output directories
OUTPUT_BASE="${PROJECT_DIR}/models"
MODEL_DIR="${OUTPUT_BASE}/${MODEL_NAME}"
HERACLES_OUTPUT="${MODEL_DIR}/heracles"
CMFGEN_OUTPUT="${MODEL_DIR}/cmfgen"
PLOTS_DIR="${MODEL_DIR}/plots"

echo "======================================================"
echo "SN2024pba Single Model Run"
echo "======================================================"
echo "Model: ${MODEL_NAME}"
echo "Description: ${DESCRIPTION}"
echo "Config: ${CONFIG_FILE_ABS}"
echo "Output directory: ${MODEL_DIR}"
echo "======================================================"

mkdir -p "${HERACLES_OUTPUT}"
mkdir -p "${CMFGEN_OUTPUT}"
mkdir -p "${PLOTS_DIR}"

# Step 1: Create HERACLES initial profile
echo ""
echo "Step 1: Creating HERACLES initial profile..."
echo "------------------------------------------------------"

M_EJECTA_A=$(read_yaml "M_ejecta_A" "$CONFIG_FILE_ABS")
M_WIND=$(read_yaml "M_wind" "$CONFIG_FILE_ABS")
M_EJECTA_B=$(read_yaml "M_ejecta_B" "$CONFIG_FILE_ABS")
FX_H_A=$(read_yaml "fx_H_A" "$CONFIG_FILE_ABS")
FX_HE_A=$(read_yaml "fx_He_A" "$CONFIG_FILE_ABS")
FX_H_WIND=$(read_yaml "fx_H_wind" "$CONFIG_FILE_ABS")
FX_HE_WIND=$(read_yaml "fx_He_wind" "$CONFIG_FILE_ABS")
FX_H_B=$(read_yaml "fx_H_B" "$CONFIG_FILE_ABS")
FX_HE_B=$(read_yaml "fx_He_B" "$CONFIG_FILE_ABS")
FX_FE_B=$(read_yaml "fx_Fe_B" "$CONFIG_FILE_ABS")
FX_SI_B=$(read_yaml "fx_Si_B" "$CONFIG_FILE_ABS")
FX_C_B=$(read_yaml "fx_C_B" "$CONFIG_FILE_ABS")
N_CELLS=$(read_yaml "n_cells" "$CONFIG_FILE_ABS")

PROFILE_OUTPUT="${PROJECT_DIR}/heracles/input_profile"

cd "${PROJECT_DIR}"
python3 create_embedded_hhe_shell.py \
    --days_ejecta 400 \
    --n_cells ${N_CELLS} \
    --output "${PROFILE_OUTPUT}" \
    --M_ejecta_A ${M_EJECTA_A} \
    --M_wind ${M_WIND} \
    --M_ejecta_B ${M_EJECTA_B} \
    --fx_H_A ${FX_H_A} \
    --fx_He_A ${FX_HE_A} \
    --fx_H_wind ${FX_H_WIND} \
    --fx_He_wind ${FX_HE_WIND} \
    --fx_H_B ${FX_H_B} \
    --fx_He_B ${FX_HE_B} \
    --fx_Fe ${FX_FE_B} \
    --fx_Si ${FX_SI_B} \
    --fx_C ${FX_C_B}

echo "Profile created: ${PROFILE_OUTPUT}"

# Step 2: Run HERACLES simulation
echo ""
echo "Step 2: Running HERACLES simulation..."
echo "------------------------------------------------------"
echo "Note: HERACLES simulation must be run manually or via separate script"
echo "HERACLES directory: ${PROJECT_DIR}/heracles"
echo "HERACLES output: ${HERACLES_OUTPUT}"
echo ""
echo "To run HERACLES:"
echo "  cd ${PROJECT_DIR}/heracles"
echo "  sbatch run.slurm"
echo ""
read -p "Press Enter to continue after HERACLES completes (or skip with Ctrl+C)..."
echo "Skipping HERACLES run for now..."

# Step 3: Convert HERACLES output to CMFGEN format
echo ""
echo "Step 3: Converting HERACLES output to CMFGEN format..."
echo "------------------------------------------------------"

HERACLES_DIR=$(read_yaml "heracles_dir" "$CONFIG_FILE_ABS")
SNAPSHOT_DAYS=$(read_yaml "snapshot_days" "$CONFIG_FILE_ABS")
CMFGEN_OUTPUT_FILE=$(read_yaml "output_file" "$CONFIG_FILE_ABS")

cd "${SCRIPT_DIR}"
python3 convert_heracles_to_cmfgen.py \
    "${HERACLES_DIR}" \
    ${SNAPSHOT_DAYS} \
    "${CMFGEN_OUTPUT_FILE}"

echo "CMFGEN input file created: ${CMFGEN_OUTPUT_FILE}"

# Step 4: Run CMFGEN
echo ""
echo "Step 4: Running CMFGEN..."
echo "------------------------------------------------------"

CMFGEN_MODEL_DIR=$(read_yaml "model_dir" "$CONFIG_FILE_ABS")
SUBMIT_SCRIPT=$(read_yaml "submit_script" "$CONFIG_FILE_ABS")

RUN_INIT=$(read_yaml "run_init" "$CONFIG_FILE_ABS")
RUN_STABLE=$(read_yaml "run_stable" "$CONFIG_FILE_ABS")
RUN_FINAL=$(read_yaml "run_final" "$CONFIG_FILE_ABS")
RUN_OBS=$(read_yaml "run_obs" "$CONFIG_FILE_ABS")
FIX_TEMPERATURE=$(read_yaml "fix_temperature" "$CONFIG_FILE_ABS")

SUBMIT_ARGS=""
if [[ "$RUN_STABLE" == "true" ]]; then
    SUBMIT_ARGS="--with-stable"
fi

# Temperature control is always allowed to vary (FIX_T=F, FIX_X_DPTH=0)

cd "${CMFGEN_MODEL_DIR}"
if [[ "$RUN_INIT" == "true" ]]; then
    echo "Submitting CMFGEN INIT job..."
    "${SUBMIT_SCRIPT}" --all ${SUBMIT_ARGS}
else
    echo "Submitting CMFGEN FINAL job (restart mode)..."
    "${SUBMIT_SCRIPT}" --restart ${SUBMIT_ARGS}
fi

echo ""
echo "======================================================"
echo "Single model run workflow complete"
echo "======================================================"
echo "Model: ${MODEL_NAME}"
echo "CMFGEN jobs submitted. Monitor with squeue."
echo "Output directory: ${MODEL_DIR}"
echo "======================================================"
