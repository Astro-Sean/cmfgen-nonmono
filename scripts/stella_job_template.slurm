#!/bin/bash
#SBATCH --job-name={MODEL_NAME}_stella
#SBATCH --output=stella_{MODEL_NAME}_%j.out
#SBATCH --error=stella_{MODEL_NAME}_%j.err
#SBATCH --time={TIME}
#SBATCH --partition={PARTITION}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem={MEMORY}
#SBATCH --chdir={STELLA_RUN_DIR}

set -euo pipefail

# Use astro conda environment with numpy
source /mpcdf/soft/RHEL_9/packages/x86_64/python-waterboa/2024.06/bin/activate /ptmp/sbrennan/conda-envs/astro
unset PYTHONPATH PYTHONHOME

echo "======================================"
echo "STELLA Job Started: {MODEL_NAME}"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURMD_NODENAME"
echo "======================================"

STELLA_ROOT="{STELLA_ROOT}"
SCRIPTS_DIR="{SCRIPTS_DIR}"
PROJECT_DIR="{PROJECT_DIR}"
MODEL_NAME="{MODEL_NAME}"
MODEL_DIR="{MODEL_DIR}"
STELLA_RUN_DIR="{STELLA_RUN_DIR}"
TARGET_DAYS="{TARGET_DAYS}"
NZON="{NZON}"
KAPPA="{KAPPA}"
TAU_MAX="{TAU_MAX}"
TAU_MIN="{TAU_MIN}"
R_MIN="{R_MIN}"
CONFIG_FILE="{CONFIG_FILE}"

INPUT_PROFILE="{INPUT_PROFILE}"

module purge
module load intel/2024.0
module load python-waterboa/2025.06
export SYSTYPE="ifort"
export HOMEStella="$STELLA_ROOT"
export PATH="$STELLA_ROOT/bin:$PATH"

# -------------------------------------------------------
# Step 1: Build STELLA if needed (only done once per NZON)
# -------------------------------------------------------
BUILD_STAMP="$STELLA_ROOT/.built_nzon{NZON}"
if [[ ! -f "$BUILD_STAMP" ]]; then
    echo "Building OpenStella for NZON=$NZON..."
    BUILD_SCRIPT="$(dirname "$STELLA_ROOT")/build_stella.sh"
    bash "$BUILD_SCRIPT" "$NZON" 2>&1 | tee stella_build.log
    touch "$BUILD_STAMP"
    echo "Build complete."
else
    echo "STELLA already built for NZON=$NZON, skipping build."
fi

XSTELLA="$STELLA_ROOT/bin/xstella.exe"
if [[ ! -f "$XSTELLA" ]]; then
    echo "ERROR: STELLA executable not found at $XSTELLA" >&2
    exit 1
fi

# -------------------------------------------------------
# Step 2: Create model-specific working directory
# -------------------------------------------------------
echo ""
echo "--- Setting up model-specific working directory ---"
MODEL_WORK_DIR="$STELLA_RUN_DIR"
mkdir -p "$MODEL_WORK_DIR"
cd "$MODEL_WORK_DIR"

# -------------------------------------------------------
# Step 3: Convert input_profile -> STELLA *.hyd + *.abn
# -------------------------------------------------------
echo ""
echo "--- Converting input profile to STELLA format ---"
MODMAKE_DIR="$MODEL_WORK_DIR/modmake"
mkdir -p "$MODMAKE_DIR"

# Use python from loaded module (has numpy)
python "$SCRIPTS_DIR/convert_profile_to_stella.py" \
    "$INPUT_PROFILE" \
    "$MODEL_NAME" \
    "$MODMAKE_DIR" \
    --nzon="$NZON" \
    --target-days="$TARGET_DAYS" \
    --snapshot-interval="{SNAPSHOT_INTERVAL}" \
    --r-min="$R_MIN" \
    --explosion-energy="{EXPLOSION_ENERGY}" \
    --explosion-mass="{EXPLOSION_MASS}" \
    --explosion-time="{EXPLOSION_TIME}" \
    --explosion-delay="{EXPLOSION_DELAY}" \
    --initial-kinetic-energy="{INITIAL_KINETIC_ENERGY}" \
    --nstmax="{NSTMAX}" \
    --eps="{EPS}" \
    --hmax="{HMAX}" \
    {ENABLE_CONV_ARG} \
    {DISABLE_SCATTERING_ARG} \
    --fit-tau="{FIT_TAU}" \
    {GIVDTL_ARG} \
    {CUSTOM_OUTPUT_TIMES_ARG} \
    --outer-fraction="{OUTER_FRACTION}" \
    --transition-radius="{TRANSITION_RADIUS}" \
    2>&1 | tee -a stella_build.log

# Create strad directory and copy .dat file
STRAD_DIR="$MODEL_WORK_DIR/strad"
mkdir -p "$STRAD_DIR"
DAT_FILE="$STRAD_DIR/{MODEL_NAME}.dat"
cp "$MODMAKE_DIR/{MODEL_NAME}.dat" "$DAT_FILE"
echo "STELLA input files ready."

# -------------------------------------------------------
# Step 4: Copy/link necessary STELLA files to model directory
# -------------------------------------------------------
echo ""
echo "--- Copying STELLA files to model directory ---"

# Create eve directory and copy eve2.exe
mkdir -p "$MODEL_WORK_DIR/eve"
if [[ -f "$STELLA_ROOT/run/eve/eve2.exe" ]]; then
    cp "$STELLA_ROOT/run/eve/eve2.exe" "$MODEL_WORK_DIR/eve/"
else
    echo "ERROR: eve2.exe not found" >&2
    exit 1
fi

# Link lineatom.dat to model vladsf directory
mkdir -p "$MODEL_WORK_DIR/vladsf"
cd "$MODEL_WORK_DIR/vladsf"
if [[ -f "$STELLA_ROOT/vladsf/lineatom.dat" ]]; then
    ln -sf "$STELLA_ROOT/vladsf/lineatom.dat" lineatom.dat
elif [[ -f "$STELLA_ROOT/run/vladsf/lineatom.dat" ]]; then
    ln -sf "$STELLA_ROOT/run/vladsf/lineatom.dat" lineatom.dat
else
    echo "WARNING: lineatom.dat not found" >&2
fi

# Generate STELLA opacity tables (test.1 through test.6) for homogeneous electron scattering
# BINARY_NZON must match zone.inc Mzon used to compile the STELLA binary (not model zone count)
BINARY_NZON=$(grep -oP '^PARAMETER\(Mzon=\K[0-9]+' "$STELLA_ROOT/src/zone.inc" 2>/dev/null | head -1 || echo "{NZON}")
echo ""
echo "--- Generating vladsf opacity tables (Mzon=$BINARY_NZON, kappa={KAPPA}) ---"
rm -f "$MODEL_WORK_DIR/vladsf/{MODEL_NAME}."{1,2,3,4,5,6} "$MODEL_WORK_DIR/vladsf/{MODEL_NAME}.ab" 2>/dev/null || true
python "$SCRIPTS_DIR/generate_stella_opacity_tables.py" \
    "$MODEL_WORK_DIR/vladsf" "{MODEL_NAME}" \
    --nzon "$BINARY_NZON" --kappa {KAPPA} 2>&1
if [ ! -s "$MODEL_WORK_DIR/vladsf/{MODEL_NAME}.1" ]; then
    echo "ERROR: Failed to generate opacity tables"
    exit 1
fi
FILE_SIZE=$(wc -c < "$MODEL_WORK_DIR/vladsf/{MODEL_NAME}.1")
echo "Opacity tables generated (${FILE_SIZE} bytes each)."

# -------------------------------------------------------
# Step 5: Convert .hyd + .abn to .mod format via eve2
# -------------------------------------------------------
echo ""
echo "--- Converting .hyd + .abn to .mod format via eve2 ---"
cd "$MODEL_WORK_DIR/eve"

# Read total mass from .hyd header for eve parameters
HYD="$MODMAKE_DIR/{MODEL_NAME}.hyd"
BMK_SOL=$(python -c "
import numpy as np
d=np.loadtxt('$HYD', skiprows=1, usecols=(6,))
print(f'{d[-1]:.4f}')
" 2>/dev/null || echo "12.0")

# Write model-specific .eve parameter file (based on sn2.eve, Kadapt=5 = homogeneous from .hyd)
cat > {MODEL_NAME}.eve <<EVEPARM
 ANTROP             DENSCE          DENSCE-old
  1.00               1.d-17              1.
 Rbeg    CRHO    RWIND   Tpwind   pwind
1.d0     1.d14   -3.D3    2.5D3   0.d0
 Rce  -- Core Radius (Solar units)
0.d0
 TPROC -- CENTRAL TEMPERATURE
1.d-3
 VSTRE -- FITTING POINT
0.55
 ALPHA
0.3333333333333
 BMK       BM1   AMeveNi
${BMK_SOL}d0  0.d0   0.
 STEP
-0.98
 expfac   velofac   rhofac   xnifac
1.d0     1.0d0     1.d0     1.d0
 Kadapt   1-FPRho   2-Polytropes  3-Isothermal  4-foreign  5-homo  6-woo87a
5
 US
+1.
 TC  ISOTP. CENTRE   TS   SURFACE TEMPERATURE
0.3                 0.3
 AMISO  GAMT
1.D0   1.6666666667
 AS       ZS
12.      6.
 XH       XHE        Depl
0.7      0.296      0.0d0
 AVERG    IWOO
28.       2
 AMZinn   AMZout   AMHein   AMHeout   BXY
0.0      0.0      0.0      0.0       0.0001
 IPRINT
2
 ITER1   ITER2
100     100
 EPS1    EPSA    EPSD
-14.0    -14.0   -10.0
 ULGCAP  ULGEPS
0.0     12.0
 H
0.00001
 DM1         DM2         PLATO        CENTR   DMID    MIDDLE
1.D-02      1.D-02      0.4D0        0.4D0   1.D-2    4.550
EVEPARM

# Write eve.1: path table for eve2 (use relative paths to avoid Fortran path length limits)
cat > eve.1 <<EVEOF
Data             Result          Model                Opacity            Rho               ExtMod                     CompMod                    Ni56 file
'{MODEL_NAME}.eve'  'res.{MODEL_NAME}'  '{MODEL_NAME}.mod'  'DummyOpafile'  '{MODEL_NAME}.rho'  '../modmake/{MODEL_NAME}.hyd'  '../modmake/{MODEL_NAME}.abn'  '{MODEL_NAME}.xni'
EVEOF

# Remove stale outputs (but keep pre-generated .xni file if it exists)
rm -f "{MODEL_NAME}.mod" "res.{MODEL_NAME}" "{MODEL_NAME}.rho"
# Only remove .xni if it wasn't pre-generated by run_model_grid.py
if [ ! -f "{MODEL_NAME}.xni" ]; then
    # .xni doesn't exist, will be generated by eve2
    :
else
    # .xni exists (pre-generated), keep it
    echo "Using pre-generated .xni file"
fi

./eve2.exe 2>&1 | tee -a "$MODEL_WORK_DIR/stella_build.log"
EXIT_EVE=$?

if [ ! -s "$MODEL_WORK_DIR/eve/{MODEL_NAME}.mod" ]; then
    echo "ERROR: eve2 failed to generate .mod file (exit=$EXIT_EVE)"
    exit 1
fi
echo ".mod file created successfully ($(wc -c < {MODEL_NAME}.mod) bytes)"

# -------------------------------------------------------
# Step 6: Set up STELLA run directory
# -------------------------------------------------------
echo ""
echo "--- Setting up STELLA run directory ---"
STRAD_DIR="$MODEL_WORK_DIR/strad"
mkdir -p "$STRAD_DIR"
cd "$STRAD_DIR"

# Link .mod/.xni files from eve directory
ln -sf "../eve/{MODEL_NAME}.mod" .
ln -sf "../eve/{MODEL_NAME}.xni" .

# Link opacity file
ln -sf "../vladsf/{MODEL_NAME}.1" .

# -------------------------------------------------------
# Step 7: Run STELLA
# -------------------------------------------------------
echo ""
echo "--- Running STELLA (target: $TARGET_DAYS days) ---"
# STELLA checks for {runname}.prf and {runname}.crv: both must be absent for BEGRUN=.TRUE.
# STELLA also cannot overwrite .res, .bm, .balance, .swd, .tt files
# Clean up files that STELLA must overwrite from both strad/ and model work directories
# Use absolute paths to ensure files are deleted from correct locations
echo "--- Cleaning up files that STELLA must overwrite ---"
for cleanup_dir in "$STRAD_DIR" "$MODEL_WORK_DIR"; do
    rm -f "$cleanup_dir/$MODEL_NAME.res" "$cleanup_dir/$MODEL_NAME.bm" \
          "$cleanup_dir/$MODEL_NAME.balance" "$cleanup_dir/$MODEL_NAME.swd" \
          "$cleanup_dir/$MODEL_NAME.tt" "$cleanup_dir/$MODEL_NAME.prf" \
          "$cleanup_dir/$MODEL_NAME.crv" 2>/dev/null || true
done

# strad.1: header + one entry (5 columns)
# STELLA always constructs Sumprf={runname}.prf and Sumcur={runname}.crv internally
# strad.1 columns 6-7 are ignored; what matters is the files don't exist (fresh run)
STRAD1="strad.1"
cat > "$STRAD1" <<STRAD
Run             Results             Model            Nickel             Opacity
{MODEL_NAME}  {MODEL_NAME}.res  {MODEL_NAME}  {MODEL_NAME}  {MODEL_NAME}.1
STRAD

echo "  strad.1 written; dat file: $DAT_FILE"

# Run STELLA executable
"$XSTELLA" < /dev/null 2>&1 | tee stella_run.log
echo "STELLA run complete."

# STELLA moves .swd and .tt to $HOMEStella/res/ at end of run
# Copy them back from STELLA res/ directory to model work directory
cp "$STELLA_ROOT/res/$MODEL_NAME.swd" "$MODEL_WORK_DIR/" 2>/dev/null || true
cp "$STELLA_ROOT/res/$MODEL_NAME.tt" "$MODEL_WORK_DIR/" 2>/dev/null || true
# Also try strad/ in case STELLA behavior changes
cp "$STRAD_DIR/$MODEL_NAME.swd" "$MODEL_WORK_DIR/" 2>/dev/null || true
cp "$STRAD_DIR/$MODEL_NAME.tt" "$MODEL_WORK_DIR/" 2>/dev/null || true

# -------------------------------------------------------
# Step 8: Locate .swd output and convert to SN_HYDRO_DATA
# -------------------------------------------------------
echo ""
echo "--- Converting STELLA .swd -> SN_HYDRO_DATA ---"

SWD_FILE="$MODEL_WORK_DIR/$MODEL_NAME.swd"
if [[ ! -f "$SWD_FILE" ]]; then
    SWD_FILE="$STRAD_DIR/$MODEL_NAME.swd"
fi
if [[ ! -f "$SWD_FILE" ]]; then
    SWD_FILE="$STELLA_ROOT/res/$MODEL_NAME.swd"
fi
if [[ ! -f "$SWD_FILE" ]]; then
    echo "ERROR: STELLA .swd output not found for model $MODEL_NAME" >&2
    exit 1
fi
echo "Found .swd: $SWD_FILE"

ABN_FILE="$MODMAKE_DIR/{MODEL_NAME}.abn"
CMFGEN_DIR="$MODEL_DIR/cmfgen"
mkdir -p "$CMFGEN_DIR"
OUTPUT_HYDRO="$CMFGEN_DIR/SN_HYDRO_DATA"

python "$SCRIPTS_DIR/convert_stella_swd_to_cmfgen.py" \
    "$SWD_FILE" \
    "$OUTPUT_HYDRO" \
    --abn-file="$ABN_FILE" \
    --target-days="{CMFGEN_SNAPSHOT_DAYS}" \
    --kappa="$KAPPA" \
    --tau-max="$TAU_MAX" \
    --tau-min="$TAU_MIN" \
    --model-name="$MODEL_NAME" \
    --grid-thinning-threshold="{GRID_THINNING_THRESHOLD}" \
    {ENABLE_GRID_THINNING_ARG} \
    {ENABLE_ADAPTIVE_GRID_ARG} \
    --base-grid-target="{BASE_GRID_TARGET}" \
    --max-radius-jump="{MAX_RADIUS_JUMP}" \
    --max-temp-jump="{MAX_TEMP_JUMP}" \
    --max-density-jump="{MAX_DENSITY_JUMP}" \
    --max-ne-jump="{MAX_NE_JUMP}" \
    --max-logtau-jump="{MAX_LOGTAU_JUMP}" \
    2>&1 | tee -a stella_build.log

echo ""
echo "SN_HYDRO_DATA written to: $OUTPUT_HYDRO"

# Copy lightcurve file (.tt) if it exists
# STELLA writes .tt alongside .swd in the strad/ subdirectory of the run directory
TT_FILE=""
for d in "$STRAD_DIR" "$MODEL_WORK_DIR/strad" "$STELLA_ROOT/run/strad/res" "$STELLA_ROOT/res" "$STELLA_ROOT/run/strad"; do
    if [[ -f "$d/{MODEL_NAME}.tt" ]]; then
        TT_FILE="$d/{MODEL_NAME}.tt"
        break
    fi
done
if [[ -z "$TT_FILE" ]]; then
    TT_FILE=$(find "$MODEL_WORK_DIR" "$STELLA_ROOT" -name "{MODEL_NAME}.tt" 2>/dev/null | head -1)
fi

# Create output directory and copy lightcurve
OUTPUT_DIR="$MODEL_DIR/output"
mkdir -p "$OUTPUT_DIR"

# Always create plots directory and plot evolution from .swd
PLOTS_DIR="$MODEL_DIR/plots"
mkdir -p "$PLOTS_DIR"

# Plot evolution (density, temperature, pressure, etc.) from .swd snapshots
python "$SCRIPTS_DIR/plot_stella_evolution.py" \
    "$SWD_FILE" \
    "$PLOTS_DIR" \
    --n-epochs 10 \
    2>&1 | tee -a stella_build.log
echo "Evolution plots: $PLOTS_DIR/{MODEL_NAME}_evolution.png, $PLOTS_DIR/{MODEL_NAME}_radial_profiles.png"

# Lightcurve and .tt handling
if [[ -n "$TT_FILE" ]]; then
    cp "$TT_FILE" "$CMFGEN_DIR/"
    cp "$TT_FILE" "$OUTPUT_DIR/"
    echo "Lightcurve copied: $TT_FILE -> $CMFGEN_DIR/"
    echo "Lightcurve copied: $TT_FILE -> $OUTPUT_DIR/"

    # Plot lightcurve
    python "$SCRIPTS_DIR/plot_stella_lightcurve.py" \
        "$CMFGEN_DIR/{MODEL_NAME}.tt" \
        --output "$PLOTS_DIR/stella_lightcurve.png" \
        --days-max "$TARGET_DAYS" \
        --csv "$OUTPUT_DIR/lightcurve.csv" \
        2>&1 | tee -a stella_build.log
    echo "Lightcurve plot: $PLOTS_DIR/stella_lightcurve.png"
    echo "Lightcurve CSV: $OUTPUT_DIR/lightcurve.csv"
else
    echo "Note: No .tt lightcurve file found for {MODEL_NAME}"
fi

# -------------------------------------------------------
# Step 9: Prepare SEDONA directory
# -------------------------------------------------------
echo ""
echo "--- Preparing SEDONA directory ---"

SEDONA_DIR="$MODEL_DIR/sedona"
mkdir -p "$SEDONA_DIR"

# Copy sedona_base files
SEDONA_BASE="$PROJECT_DIR/sedona_base"
if [[ -d "$SEDONA_BASE" ]]; then
    echo "Copying SEDONA base files from $SEDONA_BASE..."
    for f in "$SEDONA_BASE"/*; do
        if [[ -f "$f" ]]; then
            cp "$f" "$SEDONA_DIR/"
            echo "  Copied $(basename $f)"
        fi
    done
else
    echo "WARNING: sedona_base directory not found at $SEDONA_BASE"
fi

# Convert CMFGEN SN_HYDRO_DATA to SEDONA .mod format
MOD_FILE="$SEDONA_DIR/{MODEL_NAME}.mod"
if [[ -f "$OUTPUT_HYDRO" ]]; then
    echo "Converting CMFGEN SN_HYDRO_DATA to SEDONA .mod format..."
    python /ptmp/sbrennan/sedona/convert_cmfgen_to_sedona_v2.py "$OUTPUT_HYDRO" "$MOD_FILE" 2>&1 | tee -a stella_build.log
    if [[ $? -eq 0 ]]; then
        echo "Successfully created $MOD_FILE"
    else
        echo "ERROR: Conversion failed"
        exit 1
    fi
else
    echo "ERROR: SN_HYDRO_DATA not found at $OUTPUT_HYDRO"
    exit 1
fi

# Generate param.lua from template
PARAM_TEMPLATE="$SEDONA_DIR/param.lua"
if [[ -f "$PARAM_TEMPLATE" ]]; then
    # Read values from SN_HYDRO_DATA header
    # Handle negative times (pre-SBO) - these are valid and represent time before shock breakout
    R_INNER=$(python -c "
import sys
try:
    with open('$OUTPUT_HYDRO') as f:
        lines = f.readlines()
    for line in lines:
        if 'Time(days) since explosion' in line:
            parts = line.split()
            for j, p in enumerate(parts):
                if p == 'explosion:':
                    time_val = float(parts[j+1])
                    # Use absolute value for snapshot_days in param.lua
                    # (SEDONA expects positive time since explosion)
                    print(abs(time_val))
                    break
            break
except:
    print('13.0')
" 2>/dev/null || echo "13.0")

    T_INNER=$(python -c "
import sys
try:
    with open('$OUTPUT_HYDRO') as f:
        lines = f.readlines()
    section = None
    temp_data = []
    for line in lines:
        stripped = line.strip()
        if stripped.startswith('Temperature'):
            section = 'temperature'
            continue
        elif stripped == 'Density (gm/cm^3)':
            section = 'density'
            continue
        elif stripped.startswith('Sigma'):
            section = 'sigma'
            continue
        elif stripped.startswith('Electron density'):
            section = 'electron_density'
            continue
        elif stripped.startswith('Atom density'):
            section = 'atom_density'
            continue
        elif stripped == '':
            section = None
            continue
        elif section == 'temperature' and stripped:
            values = [float(x) for x in stripped.split()]
            temp_data.extend(values)
    if temp_data:
        print(f'{temp_data[-1] * 1e4:.2e}')  # Innermost T in K
    else:
        print('5.0e3')
except:
    print('5.0e3')
" 2>/dev/null || echo "5.0e3")

    R_INNER_CM=$(python -c "
import sys
try:
    with open('$OUTPUT_HYDRO') as f:
        lines = f.readlines()
    section = None
    radius_data = []
    for line in lines:
        stripped = line.strip()
        if stripped.startswith('Radius grid'):
            section = 'radius'
            continue
        elif stripped.startswith('Velocity'):
            section = 'velocity'
            continue
        elif stripped == '':
            section = None
            continue
        elif section == 'radius' and stripped and not stripped.startswith('Number'):
            values = [float(x) for x in stripped.split()]
            radius_data.extend(values)
    if radius_data:
        print(f'{radius_data[-1] * 1e10:.3e}')  # Innermost R in cm
    else:
        print('6.0e14')
except:
    print('6.0e14')
" 2>/dev/null || echo "6.0e14")

    LUMINOSITY=$(python -c "
import sys
try:
    with open('$OUTPUT_HYDRO') as f:
        lines = f.readlines()
    for line in lines:
        if 'Energy Lum [lsun]' in line:
            parts = line.split()
            for j, p in enumerate(parts):
                if p == '[lsun]:':
                    print(f'{float(parts[j+1]) * 3.828e33:.3e}')
                    break
            break
except:
    print('8.32e42')
" 2>/dev/null || echo "8.32e42")

    # Replace placeholders in param.lua
    sed -i "s/@model_name@/{MODEL_NAME}/g" "$PARAM_TEMPLATE"
    sed -i "s/@snapshot_days@/$R_INNER/g" "$PARAM_TEMPLATE"
    sed -i "s/@r_inner@/$R_INNER_CM/g" "$PARAM_TEMPLATE"
    sed -i "s/@temperature@/$T_INNER/g" "$PARAM_TEMPLATE"
    sed -i "s/@luminosity@/$LUMINOSITY/g" "$PARAM_TEMPLATE"
    
    echo "Generated param.lua with R=$R_INNER_CM cm, T=$T_INNER K, L=$LUMINOSITY erg/s"
else
    echo "WARNING: param.lua template not found in sedona_base"
fi

echo "SEDONA directory prepared: $SEDONA_DIR"

# -------------------------------------------------------
# Step 10: Prepare CMFGEN directory
# -------------------------------------------------------
echo ""
echo "--- Preparing CMFGEN directory ---"

CMFGEN_DIR="$MODEL_DIR/cmfgen"
CMFGEN_BASE="$PROJECT_DIR/cmfgen_base"

if [[ -d "$CMFGEN_BASE" ]]; then
    echo "Copying CMFGEN base files from $CMFGEN_BASE to $CMFGEN_DIR..."
    # Create cmfgen directory if it doesn't exist (don't remove - SN_HYDRO_DATA is already there)
    mkdir -p "$CMFGEN_DIR"
    # Copy base files, but don't overwrite existing files (like SN_HYDRO_DATA)
    for f in "$CMFGEN_BASE"/*; do
        if [[ -f "$f" ]]; then
            dest="$CMFGEN_DIR/$(basename $f)"
            if [[ ! -f "$dest" ]]; then
                cp "$f" "$dest"
                echo "  Copied $(basename $f)"
            else
                echo "  Skipping $(basename $f) (already exists)"
            fi
        fi
    done
    echo "CMFGEN directory prepared: $CMFGEN_DIR"
else
    echo "WARNING: cmfgen_base directory not found at $CMFGEN_BASE"
fi

echo ""
echo "======================================"
echo "STELLA Job Complete: $MODEL_NAME"
echo "All input files ready for CMFGEN and SEDONA"
echo "======================================"
