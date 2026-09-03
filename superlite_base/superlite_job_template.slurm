#!/bin/bash
#SBATCH --job-name={MODEL_NAME}_superlite
#SBATCH --output=superlite_%j.out
#SBATCH --error=superlite_%j.err
#SBATCH --time={TIME_LIMIT}
#SBATCH --partition={PARTITION}
#SBATCH --qos=interactive
#SBATCH --ntasks={NTASKS}
#SBATCH --cpus-per-task={CPUS_PER_TASK}
#SBATCH --mem={MEMORY}
#SBATCH --chdir={WORK_DIR}

set -uo pipefail

# Activate conda environment for Python (conversion script needs numpy)
source /mpcdf/soft/RHEL_9/packages/x86_64/python-waterboa/2024.06/bin/activate /ptmp/sbrennan/conda-envs/astro
unset PYTHONPATH PYTHONHOME

# Load GCC compiler for SuperLite (built with gfortran + AOCL for AMD EPYC)
module use /mpcdf/soft/RHEL_9/modules/compilers
module purge
module load gcc/14

# Set OpenMP threads to match in_nomp parameter
export OMP_NUM_THREADS={CPUS_PER_TASK}

echo "Running SuperLite for {MODEL_NAME} model"
echo "Working directory: $(pwd)"
echo "Date: $(date)"

# SuperLite executable
SUPERLITE_EXE={SUPERLITE_EXE}

# Convert hydro data to SuperLite input.str if needed.
# Two paths, selected by INPUT_SOURCE:
#   swd           : direct STELLA .swd -> input.str (MC-aware, 15 elements, native kappa)
#   sn_hydro_data : CMFGEN SN_HYDRO_DATA -> input.str (legacy, 9 grouped elements)
input_str="input.str"
cmfgen_data="{CMFGEN_DATA}"
input_source="{INPUT_SOURCE}"
stella_swd="{STELLA_SWD}"
target_days="{TARGET_DAYS}"
truncate_tau="{TRUNCATE_TAU}"
max_cells="{MAX_CELLS}"
t_floor="{T_FLOOR}"

if [ ! -f "$input_str" ]; then
    echo "No input.str found - attempting conversion (input_source=$input_source)..."
    if [ "$input_source" = "swd" ] && [ -f "$stella_swd" ]; then
        echo "Converting STELLA .swd -> SuperLite input.str (MC-aware direct path)..."
        python3 /ptmp/sbrennan/SN2024pba/scripts/convert_stella_swd_to_superlite.py \
            "$stella_swd" -o "$(pwd)" --model-name "{MODEL_BASE}" \
            --target-days "$target_days" --truncate-tau "$truncate_tau" \
            --max-cells "$max_cells" --t-floor "$t_floor" --nomp {CPUS_PER_TASK}
        if [ $? -eq 0 ]; then
            echo "Conversion complete (direct .swd path)"
        else
            echo "ERROR: .swd conversion failed"
            exit 1
        fi
    elif [ -f "$cmfgen_data" ]; then
        echo "Converting CMFGEN SN_HYDRO_DATA to SuperLite input.str format..."
        python3 /ptmp/sbrennan/SN2024pba/scripts/convert_sn_hydro_to_superlite.py "$cmfgen_data" -o "$(pwd)" --nomp {CPUS_PER_TASK}
        if [ $? -eq 0 ]; then
            echo "Conversion complete"
        else
            echo "ERROR: Conversion failed"
            exit 1
        fi
    else
        echo "ERROR: No input source found (swd=$stella_swd, sn_hydro_data=$cmfgen_data)"
        exit 1
    fi
else
    echo "Using existing input.str: $input_str"
fi

# Symlink SuperLite data files if not present
SUPERLITE_DATA={SUPERLITE_DATA}
for f in Atoms data.bf_verner data.ff_sutherland data.ion data.zeta Tinit.dat; do
    if [ ! -e "$f" ] && [ -e "$SUPERLITE_DATA/$f" ]; then
        ln -sf "$SUPERLITE_DATA/$f" .
    fi
done
# Also need input.fluxwl and input.wlgrid
SUPERLITE_INPUT={SUPERLITE_INPUT}
for f in input.fluxwl input.wlgrid; do
    if [ ! -f "$f" ] && [ -f "$SUPERLITE_INPUT/$f" ]; then
        ln -sf "$SUPERLITE_INPUT/$f" .
    fi
done

# Run SuperLite
SUPERLITE_EXIT=0
$SUPERLITE_EXE || SUPERLITE_EXIT=$?
echo "SuperLite completed (exit=$SUPERLITE_EXIT)"

# Copy spectrum to central outputs/ directory
OUTPUT_DIR="{OUTPUT_DIR}"
MODEL_BASE="{MODEL_BASE}"
EPOCH_LABEL="{EPOCH_LABEL}"

if [[ -n "$OUTPUT_DIR" && -d "$OUTPUT_DIR" ]]; then
    # SuperLite outputs spectrum in output.flx_luminos (text, wavelength-binned flux)
    # and grid in output.flx_grid (wavelength bins)
    if [[ -f "output.flx_luminos" ]]; then
        OUTNAME="$OUTPUT_DIR/$MODEL_BASE"_SUPERLITE_"$EPOCH_LABEL"
        cp -f output.flx_luminos "$OUTNAME.flx_luminos"
        echo "Copied output.flx_luminos -> $OUTNAME.flx_luminos"
        if [[ -f "output.flx_grid" ]]; then
            cp -f output.flx_grid "$OUTNAME.flx_grid"
            echo "Copied output.flx_grid -> $OUTNAME.flx_grid"
        fi
    else
        echo "WARNING: No output.flx_luminos found to copy to outputs/"
    fi
else
    echo "WARNING: Output directory $OUTPUT_DIR does not exist; skipping copy to outputs/"
fi

# Generate diagnostic plots — ALL go to model-level plots/ dir
PLOTS_DIR="{PLOTS_DIR}"
mkdir -p "$PLOTS_DIR"
EPOCH_LABEL="{EPOCH_LABEL}"

# Input structure plot (always generate - shows the input.str that was used)
# For the direct .swd path, the converter already generates input_structure.png
# and input_mass_fractions.png.  For the legacy CMFGEN path, use the comparison
# plot script that overlays SN_HYDRO_DATA and input.str.
if [[ -f "input_structure.png" ]]; then
    echo "Copying converter-generated input structure plots..."
    cp input_structure.png "$PLOTS_DIR/superlite_{EPOCH_LABEL}_input_structure.png"
    cp input_mass_fractions.png "$PLOTS_DIR/superlite_{EPOCH_LABEL}_input_mass_fractions.png" 2>/dev/null || true
else
    echo "Generating input structure plots (legacy path)..."
    python3 /ptmp/sbrennan/SN2024pba/scripts/plot_input_structure.py \
        --sn-hydro "{CMFGEN_DATA}" \
        --superlite-str "$(pwd)/input.str" \
        -o "$PLOTS_DIR" --model-name "{MODEL_BASE} {EPOCH_LABEL}" 2>/dev/null || true
fi

# Rename input structure plots to include superlite_{epoch}_ prefix (avoid double-prefix)
for f in input_structure_comparison.png input_mass_fractions_comparison.png \
         sn_hydro_data_input_structure.png; do
    if [[ -f "$PLOTS_DIR/$f" ]]; then
        base="${f%.png}"
        mv "$PLOTS_DIR/$f" "$PLOTS_DIR/superlite_{EPOCH_LABEL}_${base}.png" 2>/dev/null || true
    fi
done
# superlite_input_structure.png already has the prefix
if [[ -f "$PLOTS_DIR/superlite_input_structure.png" ]]; then
    mv "$PLOTS_DIR/superlite_input_structure.png" "$PLOTS_DIR/superlite_{EPOCH_LABEL}_input_structure.png" 2>/dev/null || true
fi

if [[ -f "output.flx_luminos" ]]; then
    echo "Generating SuperLite diagnostic plots..."
    python3 /ptmp/sbrennan/SN2024pba/scripts/plot_superlite_output.py "$(pwd)" --output "$PLOTS_DIR" --model-name "{MODEL_BASE} {EPOCH_LABEL}" --epoch "{EPOCH_LABEL}"
    if [ $? -eq 0 ]; then
        echo "Plots saved to $PLOTS_DIR"
    else
        echo "WARNING: Plot generation failed"
    fi
fi

echo "Date: $(date)"
