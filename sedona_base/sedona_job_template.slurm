#!/bin/bash
#SBATCH --job-name={MODEL_NAME}_sedona
#SBATCH --output=sedona_%j.out
#SBATCH --error=sedona_%j.err
#SBATCH --time={TIME_LIMIT}
#SBATCH --partition={PARTITION}
#SBATCH --qos=n0008
#SBATCH --ntasks={NTASKS}
#SBATCH --cpus-per-task=1
#SBATCH --mem={MEMORY}
#SBATCH --chdir={WORK_DIR}

set -euo pipefail

# Use system python (numpy available via waterboa module)
unset PYTHONPATH PYTHONHOME

# Load required modules for SEDONA
module purge
module load intel/2024.0
module load openmpi/4.1
module load gsl/2.6
module load hdf5-serial/1.14.1
module load python-waterboa/2025.06

export SEDONA_HOME=/ptmp/sbrennan/sedona

# Set library paths using loaded modules
export LD_LIBRARY_PATH=$(gsl-config --prefix)/lib:$(dirname $(dirname $(which h5cc)))/lib:$(pwd):$LD_LIBRARY_PATH

# Set OpenMP threads (1 per MPI task to avoid oversubscription)
export OMP_NUM_THREADS=1

echo "Running SEDONA for {MODEL_NAME} model"
echo "SEDONA_HOME: $SEDONA_HOME"
echo "Working directory: $(pwd)"
echo "Date: $(date)"

# Convert to SEDONA .mod format if needed
# Priority: 1) Use existing .mod file,
#           2) Convert from STELLA .swd directly (convert_stella_swd_to_sedona.py),
#           3) Convert from CMFGEN SN_HYDRO_DATA (convert_cmfgen_to_sedona_v2.py).
mod_file="{MOD_FILE}"
cmfgen_data="{CMFGEN_DATA}"
stella_swd="{STELLA_SWD}"
target_days="{SEDONA_TARGET_DAYS}"
sedona_tau_max="{SEDONA_TAU_MAX}"
sedona_t_floor="{SEDONA_T_FLOOR}"

if [ ! -f "$mod_file" ]; then
    echo "No .mod file found - attempting conversion..."

    # 2) Direct STELLA .swd -> .mod (one step, no CMFGEN dependency)
    if [ -f "$stella_swd" ]; then
        echo "Converting STELLA .swd -> SEDONA .mod (direct)..."
        python3 {WORK_DIR}/../../scripts/convert_stella_swd_to_sedona.py \
            "$stella_swd" "$mod_file" \
            --target-days "$target_days" --time-days "$target_days" \
            --tau-max "$sedona_tau_max" --t-floor "$sedona_t_floor"
        if [ $? -eq 0 ] && [ -f "$mod_file" ]; then
            echo "Direct .swd conversion complete"
        else
            echo "WARNING: Direct .swd conversion failed; trying CMFGEN SN_HYDRO_DATA path"
            # 3) Fallback: CMFGEN SN_HYDRO_DATA -> .mod (v2, correct rmin/texp header)
            if [ -f "$cmfgen_data" ]; then
                echo "Converting CMFGEN SN_HYDRO_DATA to SEDONA .mod format..."
                python3 /ptmp/sbrennan/sedona/convert_cmfgen_to_sedona_v2.py "$cmfgen_data" "$mod_file"
                if [ $? -eq 0 ]; then
                    echo "CMFGEN conversion complete"
                else
                    echo "ERROR: CMFGEN conversion failed"
                    exit 1
                fi
            else
                echo "ERROR: No SN_HYDRO_DATA found at $cmfgen_data and .swd conversion failed"
                echo "  Skipping this SEDONA job."
                exit 1
            fi
        fi
    # 3) CMFGEN SN_HYDRO_DATA -> .mod (when no .swd available)
    elif [ -f "$cmfgen_data" ]; then
        echo "Converting CMFGEN SN_HYDRO_DATA to SEDONA .mod format..."
        python3 /ptmp/sbrennan/sedona/convert_cmfgen_to_sedona_v2.py "$cmfgen_data" "$mod_file"
        if [ $? -eq 0 ]; then
            echo "CMFGEN conversion complete"
        else
            echo "ERROR: CMFGEN conversion failed"
            exit 1
        fi
    else
        echo "ERROR: No .swd found at $stella_swd and no SN_HYDRO_DATA at $cmfgen_data"
        echo "  This snapshot may be beyond the epoch STELLA reached before crashing."
        echo "  Skipping this SEDONA job."
        exit 1
    fi
else
    echo "Using existing .mod file: $mod_file"
fi

# Fix param.lua (zero-abundance filtering, core radius safety)
python3 /ptmp/sbrennan/SN2024pba/scripts/fix_sedona_param.py "$(pwd)"

# Run SEDONA spectrum calculation with MPI parallelism
# --bind-to none: disable core binding to avoid "binding more processes than cpus" errors
# on nodes where the binding domain has fewer cores than MPI tasks
# Capture exit code (don't let set -e abort before plotting)
SEDONA_EXIT=0
mpirun --bind-to none -np $SLURM_NTASKS $SEDONA_HOME/src/sedona6.ex param.lua || SEDONA_EXIT=$?
echo "SEDONA completed (exit=$SEDONA_EXIT) - generating plots"

# All plots go to model-level plots/ dir
PLOTS_DIR="{PLOTS_DIR}"
mkdir -p "$PLOTS_DIR"

# Generate input structure plots (shows the .mod file that was used)
echo "Generating input structure plots..."
python3 {WORK_DIR}/../../scripts/plot_input_structure.py \
    --sedona-mod "$(pwd)/{MOD_FILE}" \
    -o "$PLOTS_DIR" --model-name "{MODEL_BASE}" 2>/dev/null || true

# Rename input structure plots to include sedona_ prefix (avoid double-prefix)
for f in input_structure_comparison.png input_mass_fractions_comparison.png \
         sn_hydro_data_input_structure.png; do
    if [[ -f "$PLOTS_DIR/$f" ]]; then
        base="${f%.png}"
        mv "$PLOTS_DIR/$f" "$PLOTS_DIR/sedona_${base}.png" 2>/dev/null || true
    fi
done
# sedona_input_structure.png already has the prefix
if [[ -f "$PLOTS_DIR/sedona_input_structure.png" ]]; then
    mv "$PLOTS_DIR/sedona_input_structure.png" "$PLOTS_DIR/sedona_input_structure.png.tmp" 2>/dev/null || true
fi

# Auto-plot the spectrum (saves to model-level plots/ via get_plots_dir())
python3 plot_spectrum.py

echo "Plots generated"

# Copy final spectrum and plot to central outputs/ directory
# Python template placeholders are substituted before bash runs:
#   OUTPUT_DIR, MODEL_BASE are set as bash variables below.
# EPOCH_LABEL is read from param.lua at runtime for the actual snapshot epoch.
# Do NOT use \${VAR} bash syntax for these — it conflicts with Python substitution.
OUTPUT_DIR="{OUTPUT_DIR}"
MODEL_BASE="{MODEL_BASE}"

# Read actual epoch from param.lua (tstep_time_start in days)
EPOCH_DAYS=$(grep "tstep_time_start" param.lua 2>/dev/null | head -1 | sed 's/.*= *//;s/ .*//;s/\*.*//')
if [[ -n "$EPOCH_DAYS" ]]; then
    EPOCH_LABEL=$(python3 -c "print(f'{float(\"$EPOCH_DAYS\"):g}d')" 2>/dev/null || echo "{EPOCH_LABEL}")
else
    EPOCH_LABEL="{EPOCH_LABEL}"
fi
echo "Epoch label for outputs: $EPOCH_LABEL"

if [[ -n "$OUTPUT_DIR" && -d "$OUTPUT_DIR" ]]; then
    # Find the highest-numbered spectrum_N.dat (final iteration)
    FINAL_SPEC=$(ls -v spectrum_*.dat 2>/dev/null | tail -1)
    if [[ -n "$FINAL_SPEC" ]]; then
        OUTNAME="$OUTPUT_DIR/$MODEL_BASE"_SEDONA_"$EPOCH_LABEL"
        cp -f "$FINAL_SPEC" "$OUTNAME.dat"
        echo "Copied $FINAL_SPEC -> $OUTNAME.dat"
    else
        echo "WARNING: No spectrum_*.dat found to copy to outputs/"
    fi
    # Copy the auto-plot
    if [[ -f "spectrum_auto_plot.png" ]]; then
        cp -f "spectrum_auto_plot.png" "$OUTNAME.png"
        echo "Copied spectrum_auto_plot.png -> $OUTNAME.png"
    else
        echo "WARNING: No spectrum_auto_plot.png found to copy to outputs/"
    fi
else
    echo "WARNING: Output directory $OUTPUT_DIR does not exist; skipping copy to outputs/"
fi

echo "Date: $(date)"
