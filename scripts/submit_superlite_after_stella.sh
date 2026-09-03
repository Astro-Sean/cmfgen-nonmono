#!/bin/bash
# Submit SuperLite jobs for all 9 SN_B_grid models after their STELLA jobs complete.
# Usage: bash submit_superlite_after_stella.sh
#
# This script checks if STELLA .swd files have been updated (reaching 60d),
# then submits SuperLite jobs with --skip-stella --run-superlite.

set -euo pipefail

PROJECT_DIR="/viper/ptmp1/sbrennan/SN2024pba"
GRID_DIR="/ptmp/sbrennan/SN2024pba/SN_B_grid"
MODELS=(
    M1.0_V3000 M1.0_V5000 M1.0_V7000
    M3.0_V3000 M3.0_V5000 M3.0_V7000
    M5.0_V3000 M5.0_V5000 M5.0_V7000
)

echo "=== SuperLite submission script ==="
echo "Checking STELLA completion for ${#MODELS[@]} models..."
echo ""

all_done=true
for model in "${MODELS[@]}"; do
    swd="${GRID_DIR}/${model}/stella/${model}.swd"
    tt="${GRID_DIR}/${model}/stella/${model}.tt"

    if [ ! -f "$swd" ]; then
        echo "  ${model}: .swd not found — STELLA still running or failed"
        all_done=false
        continue
    fi

    # Check if STELLA reached 60d by looking at the .tt file
    if [ ! -f "$tt" ]; then
        echo "  ${model}: .tt not found — STELLA may still be running"
        all_done=false
        continue
    fi

    # Get the last time value from .tt (first column of last data line)
    last_time=$(grep -a "^[0-9]" "$tt" 2>/dev/null | tail -1 | awk '{print $1}')
    if [ -z "$last_time" ]; then
        echo "  ${model}: cannot parse .tt — STELLA may still be running"
        all_done=false
        continue
    fi

    # Check if reached 60d (with 5d tolerance)
    reached=$(python3 -c "print(1 if float('${last_time}') >= 55.0 else 0)")
    if [ "$reached" = "1" ]; then
        echo "  ${model}: STELLA reached ${last_time}d — ready for SuperLite"
    else
        echo "  ${model}: STELLA at ${last_time}d — not yet at 60d"
        all_done=false
    fi
done

echo ""
if [ "$all_done" = "false" ]; then
    echo "Not all STELLA runs are complete. Exiting."
    echo "Re-run this script later once STELLA jobs finish."
    exit 1
fi

echo "All STELLA runs complete. Submitting SuperLite jobs..."
echo ""

# Submit SuperLite for each model (skip STELLA, only run SuperLite)
for model in "${MODELS[@]}"; do
    cfg="${GRID_DIR}/${model}/config.yaml"
    echo "--- Submitting SuperLite for ${model} ---"
    cd "$PROJECT_DIR"
    python scripts/run_model_grid.py "$cfg" --skip-stella --run-superlite 2>&1 | tail -20
    echo ""
done

echo "=== All SuperLite jobs submitted ==="
echo "Monitor with: squeue -u sbrennan"
