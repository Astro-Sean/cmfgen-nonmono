#!/bin/bash
#SBATCH --job-name={MODEL_NAME}_stella
#SBATCH --output=stella_{MODEL_NAME}_%j.out
#SBATCH --error=stella_{MODEL_NAME}_%j.err
#SBATCH --time={TIME}
#SBATCH --partition={PARTITION}
#SBATCH --qos=a0008
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={CPUS_PER_TASK}
#SBATCH --mem={MEMORY}
#SBATCH --chdir={STELLA_RUN_DIR}

set -euo pipefail

# Use astro conda environment with numpy.  Paths differ between clusters:
# viper has waterboa under RHEL_9 and the env at /ptmp/sbrennan/conda-envs;
# raven has waterboa under SLE_15 and the env on the /viper mount.
if [ -f /mpcdf/soft/SLE_15/packages/x86_64/python-waterboa/2024.06/bin/activate ]; then
    WATERBOA_ROOT=/mpcdf/soft/SLE_15/packages/x86_64/python-waterboa/2024.06
else
    WATERBOA_ROOT=/mpcdf/soft/RHEL_9/packages/x86_64/python-waterboa/2024.06
fi
for _env in /viper/ptmp1/sbrennan/conda-envs/astro \
            /ptmp/sbrennan/conda-envs/astro \
            /raven/ptmp/sbrennan/conda-envs/astro; do
    if [ -d "$_env" ]; then ASTRO_ENV="$_env"; break; fi
done
if [ -n "${ASTRO_ENV:-}" ]; then
    source "$WATERBOA_ROOT/bin/activate" "$ASTRO_ENV"
else
    # No conda env on this cluster — the python-waterboa module (loaded
    # below) provides numpy/scipy/yaml directly.
    source "$WATERBOA_ROOT/bin/activate"
fi
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
RESTART_STELLA="{RESTART_STELLA}"
STELLA_RUN_DIR="{STELLA_RUN_DIR}"
TARGET_DAYS="{TARGET_DAYS}"
NZON="{NZON}"
KAPPA="{KAPPA}"
STELLA_TAU_MAX="{TAU_MAX}"
STELLA_TAU_MIN="{TAU_MIN}"
R_MIN="{R_MIN}"
CONFIG_FILE="{CONFIG_FILE}"

INPUT_PROFILE="{INPUT_PROFILE}"

module purge
module load intel/2024.0
# intel/2024.0 does not set MKLROOT; the PARDISO build needs the mkl module.
# Must be 2026.0: xstella.exe links libmkl_rt.so.3 (2024.x only has .so.2).
module load mkl/2026.0
module load python-waterboa/2024.06
export SYSTYPE="ifort"
export HOMEStella="$STELLA_ROOT"
export PATH="$STELLA_ROOT/bin:$PATH"
export OMP_NUM_THREADS={OMP_NUM_THREADS}
export MKL_NUM_THREADS={OMP_NUM_THREADS}
export MKL_THREADING_LAYER=intel
# Bitwise reproducibility of the PARDISO solves: MKL_DYNAMIC=TRUE (default)
# lets MKL shrink the thread count under node load, and threaded PARDISO then
# returns different results run to run (verified: identical inputs diverged at
# ~1e-15 within ~45 steps). CBWR selects a deterministic code path.
export MKL_DYNAMIC=FALSE
export MKL_CBWR=AUTO
# Ensure MKL runtime library path is available for PARDISO-linked builds
# (module purge may clear LD_LIBRARY_PATH before module load can restore it)
for _mklpath in \
    "${MKLROOT:-}/lib/intel64" \
    "/viper/u1/system/soft/RHEL_9/packages/x86_64/intel_oneapi/2026.0/mkl/2026.0/lib/intel64" \
    "/mpcdf/soft/RHEL_9/packages/x86_64/intel_oneapi/2024.0/mkl/latest/lib/intel64"; do
    if [[ -d "$_mklpath" ]]; then
        export LD_LIBRARY_PATH="$_mklpath:${LD_LIBRARY_PATH:-}"
        break
    fi
done

# -------------------------------------------------------
# Step 1: Build STELLA if needed (only done once per NZON)
# -------------------------------------------------------
BUILD_STAMP="$STELLA_ROOT/.built_nzon{NZON}"
if [[ ! -f "$BUILD_STAMP" ]]; then
    echo "Building OpenStella for NZON=$NZON..."
    BUILD_SCRIPT="$(dirname "$STELLA_ROOT")/build_stella.sh"
    # Use parallel build script if STELLA_ROOT is the parallel variant
    if [[ "$STELLA_ROOT" == *"_parallel" ]]; then
        BUILD_SCRIPT="$(dirname "$STELLA_ROOT")/build_stella_parallel.sh"
    fi
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
#         (or use pre-generated KEPLER .hyd/.abn if present)
# -------------------------------------------------------
echo ""
echo "--- Converting input profile to STELLA format ---"
MODMAKE_DIR="$MODEL_WORK_DIR/modmake"
mkdir -p "$MODMAKE_DIR"

# Check if .hyd and .abn already exist (KEPLER pre-SN mode pre-generates them)
# KEPLER_MODE is set by run_model_grid.py based on whether 'kepler_presn' is in the config
HYD_PREGEN="$MODMAKE_DIR/{MODEL_NAME}.hyd"
ABN_PREGEN="$MODMAKE_DIR/{MODEL_NAME}.abn"
KEPLER_MODE="{KEPLER_MODE}"

if [[ "$KEPLER_MODE" == "true" && -f "$HYD_PREGEN" && -f "$ABN_PREGEN" ]]; then
    echo "Using pre-generated .hyd and .abn files"
    DAT_PREGEN="$MODMAKE_DIR/{MODEL_NAME}.dat"
    # Always regenerate the .dat via convert_profile_to_stella.py: the KEPLER
    # write_dat only carries a reduced parameter set (hardcoded FitTau/AQ and
    # default hmax/iout), so a stale pre-generated .dat would silently override
    # the resolved YAML stella: parameters.  Regenerating here guarantees the
    # .dat matches the resolved config.  A pre-existing .dat is preserved as
    # .kepler_bak for provenance but is NOT restored.
    # Back up the KEPLER .hyd, .abn, .dat (convert_profile_to_stella.py will
    # overwrite them; only .hyd/.abn are restored afterwards).
    cp "$HYD_PREGEN" "${HYD_PREGEN}.kepler_bak"
    cp "$ABN_PREGEN" "${ABN_PREGEN}.kepler_bak"
    if [[ -f "$DAT_PREGEN" ]]; then
        cp "$DAT_PREGEN" "${DAT_PREGEN}.kepler_bak"
    fi
    # Run convert_profile_to_stella.py to generate .dat with all correct parameters
    # (the KEPLER write_dat has fewer parameters; this ensures full parameter coverage)
    python "$SCRIPTS_DIR/convert_profile_to_stella.py" \
        "$INPUT_PROFILE" \
        "$MODEL_NAME" \
        "$MODMAKE_DIR" \
        --nzon="$NZON" \
        --nzon-max="{NZON_MAX}" \
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
        --hmin="{HMIN}" \
        --meth="{METH}" \
        {ENABLE_CONV_ARG} \
        {DISABLE_SCATTERING_ARG} \
        --fit-tau="{FIT_TAU}" \
        --tau-tol="{TAU_TOL}" \
        {GIVDTL_ARG} \
        {CUSTOM_OUTPUT_TIMES_ARG} \
        {IOUT_ARG} \
        {NOUT_ARG} \
        --outer-fraction="{OUTER_FRACTION}" \
        --transition-radius="{TRANSITION_RADIUS}" \
        --tau-min="{TAU_MIN}" \
        {RHO_MIN_ARG} \
        --kappa="{KAPPA}" \
        --wacc-r="{WACC_R}" \
        --wacc-v="{WACC_V}" \
        --wacc-t="{WACC_T}" \
        --wacc-rad="{WACC_RAD}" \
        --temp-floor="{TEMP_FLOOR}" \
        --floor-rad="{FLOOR_RAD}" \
        --aq="{AQ}" \
        --bq="{BQ}" \
        --min-dr-r="{MIN_DR_R}" \
        --amni="{NI56_DIST_MASS}" \
        --xmni="{NI56_TOTAL_MASS}" \
        --tcura-days="{TCURA_DAYS}" \
        2>&1 | tee -a stella_build.log
    # If the conversion left a missing/empty .dat, restore the pre-generated
    # one so the run can still proceed (with a loud warning).
    if [[ ! -s "$DAT_PREGEN" && -f "${DAT_PREGEN}.kepler_bak" ]]; then
        echo "WARNING: .dat regeneration failed — restoring pre-generated .dat" >&2
        cp "${DAT_PREGEN}.kepler_bak" "$DAT_PREGEN"
    fi
    # Restore the pre-generated .hyd/.abn only.  The freshly generated .dat
    # (full resolved-config parameters) is kept; the pre-generated .dat is
    # deliberately NOT restored.  Backups are removed.
    if [[ -f "${HYD_PREGEN}.kepler_bak" ]]; then
        cp "${HYD_PREGEN}.kepler_bak" "$HYD_PREGEN"
        cp "${ABN_PREGEN}.kepler_bak" "$ABN_PREGEN"
        rm -f "${HYD_PREGEN}.kepler_bak" "${ABN_PREGEN}.kepler_bak" "${DAT_PREGEN}.kepler_bak"
        echo "Pre-generated .hyd/.abn restored (.dat regenerated from resolved config)"
    fi
else
    # Standard mode: convert input_profile -> .hyd + .abn + .dat
    python "$SCRIPTS_DIR/convert_profile_to_stella.py" \
        "$INPUT_PROFILE" \
        "$MODEL_NAME" \
        "$MODMAKE_DIR" \
        --nzon="$NZON" \
        --nzon-max="{NZON_MAX}" \
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
        --hmin="{HMIN}" \
        --meth="{METH}" \
        {ENABLE_CONV_ARG} \
        {DISABLE_SCATTERING_ARG} \
        --fit-tau="{FIT_TAU}" \
        --tau-tol="{TAU_TOL}" \
        {GIVDTL_ARG} \
        {CUSTOM_OUTPUT_TIMES_ARG} \
        {IOUT_ARG} \
        {NOUT_ARG} \
        --outer-fraction="{OUTER_FRACTION}" \
        --transition-radius="{TRANSITION_RADIUS}" \
        --tau-min="{TAU_MIN}" \
        {RHO_MIN_ARG} \
        --kappa="{KAPPA}" \
        --wacc-r="{WACC_R}" \
        --wacc-v="{WACC_V}" \
        --wacc-t="{WACC_T}" \
        --wacc-rad="{WACC_RAD}" \
        --temp-floor="{TEMP_FLOOR}" \
        --floor-rad="{FLOOR_RAD}" \
        --aq="{AQ}" \
        --bq="{BQ}" \
        --min-dr-r="{MIN_DR_R}" \
        --amni="{NI56_DIST_MASS}" \
        --xmni="{NI56_TOTAL_MASS}" \
        --tcura-days="{TCURA_DAYS}" \
        2>&1 | tee -a stella_build.log
fi

python - "$HYD_PREGEN" "$ABN_PREGEN" "$NZON" <<'PY'
import sys
import numpy as np
hyd = np.loadtxt(sys.argv[1], skiprows=1)
abn = np.loadtxt(sys.argv[2])
nzon = int(sys.argv[3])
if hyd.shape != (nzon, 8):
    raise SystemExit(f'Invalid .hyd shape {hyd.shape}; expected ({nzon}, 8)')
if abn.shape != (nzon, 20):
    raise SystemExit(f'Invalid .abn shape {abn.shape}; expected ({nzon}, 20)')
if np.any(hyd[:, 1] <= 0) or np.any(np.diff(hyd[:, 2]) <= 0) or np.any(np.diff(hyd[:, 6]) < 0):
    raise SystemExit('Invalid .hyd mass or radial ordering')
composition_sum = abn[:, 4:20].sum(axis=1)
if np.any(abn[:, 4:20] < -1e-12) or np.any(np.abs(composition_sum - 1.0) > 2e-3):
    raise SystemExit(f'Invalid .abn composition sums: [{composition_sum.min():.6f}, {composition_sum.max():.6f}]')
print(f'STELLA input validation passed: {nzon} zones, ejecta mass={hyd[:, 1].sum():.5f} Msun, X(Ni56) max={abn[:, 19].max():.5f}')
PY

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

# Prepare vladsf directory and link lineatom.dat (needed by both grey and ronfict paths)
# LINEATOM_FILE (stella.line_list in YAML) overrides the default VVF/Kurucz list.
mkdir -p "$MODEL_WORK_DIR/vladsf"
cd "$MODEL_WORK_DIR/vladsf"
if [[ -n "{LINEATOM_FILE}" && -f "{LINEATOM_FILE}" ]]; then
    ln -sf "{LINEATOM_FILE}" lineatom.dat
    echo "  Using custom line list: {LINEATOM_FILE}"
elif [[ -f "$STELLA_ROOT/vladsf/lineatom.dat" ]]; then
    ln -sf "$STELLA_ROOT/vladsf/lineatom.dat" lineatom.dat
elif [[ -f "$STELLA_ROOT/run/vladsf/lineatom.dat" ]]; then
    ln -sf "$STELLA_ROOT/run/vladsf/lineatom.dat" lineatom.dat
else
    echo "WARNING: lineatom.dat not found" >&2
fi

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
# Back up pre-generated .xni before eve2 (eve2 overwrites it with 1e-6-floored values)
XNI_BACKUP=""
if [ -f "{MODEL_NAME}.xni" ]; then
    cp "{MODEL_NAME}.xni" "{MODEL_NAME}.xni.pregen"
    XNI_BACKUP="{MODEL_NAME}.xni.pregen"
    echo "Backed up pre-generated .xni file before eve2"
fi

./eve2.exe 2>&1 | tee -a "$MODEL_WORK_DIR/stella_build.log"
EXIT_EVE=$?

if [ ! -s "$MODEL_WORK_DIR/eve/{MODEL_NAME}.mod" ]; then
    echo "ERROR: eve2 failed to generate .mod file (exit=$EXIT_EVE)"
    exit 1
fi
echo ".mod file created successfully ($(wc -c < {MODEL_NAME}.mod) bytes)"

# Overwrite the .xni file with the correct Ni56 distribution.
# eve2 floors Ni56 at 1e-6 for all zones (evesn2.f: XNifor=max(xnifac*XNifor,1.d-6)),
# so when .abn has Ni56=0, eve2 writes 1e-6 everywhere -> STELLA sees ~1.7e-6 Msun
# instead of the intended Ni56 mass.  Restore the pre-generated .xni (backed up
# before eve2, or from modmake/ in KEPLER mode) after eve2 overwrites it.
if [ "{ZERO_NI56}" = "true" ]; then
    echo "Zeroing .xni file (no radioactive decay — Ni56 disabled)"
    python3 -c "
import struct
nzon = {NZON}
rec_len = nzon * 8
with open('{MODEL_NAME}.xni', 'wb') as f:
    f.write(struct.pack('i', rec_len))
    f.write(struct.pack(f'{nzon}d', *([0.0] * nzon)))
    f.write(struct.pack('i', rec_len))
print(f'  Written zero .xni: {{4 + rec_len + 4}} bytes, {nzon} zones (Fortran unformatted REAL*8)')
"
elif [ -f "$MODMAKE_DIR/{MODEL_NAME}.xni" ]; then
    echo "Restoring correct .xni file (Ni56 distribution) from modmake/ after eve2"
    cp "$MODMAKE_DIR/{MODEL_NAME}.xni" "{MODEL_NAME}.xni"
    python3 -c "
import struct, numpy as np
with open('{MODEL_NAME}.xni', 'rb') as f:
    data = f.read()
rec_len = struct.unpack('i', data[:4])[0]
nzon = rec_len // 8
xnifor = np.array(struct.unpack(f'{nzon}d', data[4:4+nzon*8]))
print(f'  .xni restored: {nzon} zones, max XNifor={xnifor.max():.4f}, non-zero={(xnifor > 1e-30).sum()}')
"
elif [ -n "$XNI_BACKUP" ] && [ -f "$XNI_BACKUP" ]; then
    echo "Restoring pre-generated .xni file (Ni56 distribution) after eve2"
    cp "$XNI_BACKUP" "{MODEL_NAME}.xni"
    python3 -c "
import struct, numpy as np
with open('{MODEL_NAME}.xni', 'rb') as f:
    data = f.read()
rec_len = struct.unpack('i', data[:4])[0]
nzon = rec_len // 8
xnifor = np.array(struct.unpack(f'{nzon}d', data[4:4+nzon*8]))
print(f'  .xni restored: {nzon} zones, max XNifor={xnifor.max():.4f}, non-zero={(xnifor > 1e-30).sum()}')
"
else
    echo "WARNING: No pre-generated .xni found — using eve2 output (Ni56~1e-6)"
fi
rm -f "{MODEL_NAME}.xni.pregen" 2>/dev/null || true

# -------------------------------------------------------
# Step 5b: Generate vladsf opacity tables
#   USE_RONFICT=true  -> xronfict.exe (full Saha + per-frequency bf/ff/line physics)
#   USE_RONFICT=false -> generate_stella_opacity_tables.py (grey Rosseland mean)
# -------------------------------------------------------
echo ""
USE_RONFICT="{USE_RONFICT}"
RONFICT_EXE="{RONFICT_EXE}"
YAKOVLEV_DIR="{YAKOVLEV_DIR}"
BINARY_NZON=$(grep -oP '^PARAMETER\(Mzon=\K[0-9]+' "$STELLA_ROOT/src/zone.inc" 2>/dev/null | head -1 || echo "{NZON}")

if [[ "$USE_RONFICT" == "true" ]]; then
    echo "--- Generating multigroup vladsf opacity tables via xronfict.exe ---"
    if [[ ! -x "$RONFICT_EXE" ]]; then
        echo "ERROR: xronfict.exe not found or not executable at $RONFICT_EXE" >&2
        exit 1
    fi

    cd "$MODEL_WORK_DIR/vladsf"

    # Link yakovlev cross-section data directory
    if [[ -d "$YAKOVLEV_DIR" ]]; then
        ln -sfn "$YAKOVLEV_DIR" yakovlev
    else
        echo "WARNING: yakovlev dir not found at $YAKOVLEV_DIR" >&2
    fi

    # Prefer pre-compiled linedata.dump if available (much faster than lineatom.dat)
    # SKIP_LINEDATA_DUMP=true when a custom line list is in use: the shared
    # dump is a binary image of the DEFAULT list and would silently shadow it.
    if [[ "{SKIP_LINEDATA_DUMP}" == "true" ]]; then
        rm -f linedata.dump
        echo "  Custom line list in use: skipping shared linedata.dump (fresh dump will be built from lineatom.dat)"
    else
        LINEDATA_DUMP="$STELLA_ROOT/models/X_STELA/vladsf/linedata.dump"
        if [[ -f "$LINEDATA_DUMP" ]]; then
            ln -sf "$LINEDATA_DUMP" linedata.dump
            echo "  Using pre-compiled linedata.dump"
        fi
    fi

    # Remove stale opacity files
    rm -f "{MODEL_NAME}."{1,2,3,4,5,6} "{MODEL_NAME}.ab" 2>/dev/null || true

    # Write ronfict.1 input file (blank comment line + data line)
    # Paths are relative to the vladsf run directory
    cat > ronfict.1 <<RONFEOF
  Model                          NiDist                        results          opafile
'../eve/{MODEL_NAME}.mod'  '../eve/{MODEL_NAME}.xni'  'ronfict.res'  '{MODEL_NAME}.0'
RONFEOF

    # xronfict.exe runs for ~3 min per time slice * 6 slices = ~18-20 min total.
    # It writes one file per slice (.1 first, then .2-.6) only after completing all
    # 50 zones * 14 rho * 14 T = 9800 opacity calls per slice.
    # Each call invokes the full Saha solver + bf/ff/line physics.
    # DO NOT interrupt: all 6 files + .ab must be present and non-empty for STELLA.
    echo "  Running xronfict.exe (expected ~18-20 min for 6 time slices)..."
    echo "  Progress is logged to ronfict.log; each 'done time=' line = 1 completed slice."
    echo "  xronfict start: $(date)"

    # Run xronfict.exe synchronously (no background monitoring to avoid pipefail issues)
    "$RONFICT_EXE" 2>&1 | tee ronfict.log
    EXIT_RONF=$?

    echo "  xronfict end: $(date) (exit=$EXIT_RONF)"

    # Validate: KNadap=-4 requires ALL 6 frequency files + .ab to be non-empty
    RONF_OK=true
    for nf in 1 2 3 4 5 6; do
        if [[ ! -s "{MODEL_NAME}.$nf" ]]; then
            echo "ERROR: xronfict.exe produced empty/missing opacity file {MODEL_NAME}.$nf" >&2
            RONF_OK=false
        fi
    done
    if [[ ! -s "{MODEL_NAME}.ab" ]]; then
        echo "ERROR: xronfict.exe produced empty/missing absorption file {MODEL_NAME}.ab" >&2
        RONF_OK=false
    fi
    if [[ "$RONF_OK" != "true" ]]; then
        echo "  Check $MODEL_WORK_DIR/vladsf/ronfict.log for details" >&2
        echo "  Last 20 lines of ronfict.log:" >&2
        tail -20 ronfict.log >&2
        exit 1
    fi
    FILE_SIZE=$(wc -c < "{MODEL_NAME}.1")
    echo "Multigroup opacity tables generated via xronfict: all 7 files present (${FILE_SIZE} bytes each)."

else
    echo "--- Generating vladsf opacity tables via generate_stella_opacity_tables.py (Mzon=$BINARY_NZON, kappa={KAPPA}) ---"
    cd "$MODEL_WORK_DIR/vladsf"
    rm -f "{MODEL_NAME}."{1,2,3,4,5,6} "{MODEL_NAME}.ab" 2>/dev/null || true
    "$WATERBOA_ROOT/bin/python3.12" "$SCRIPTS_DIR/generate_stella_opacity_tables.py" \
        "$MODEL_WORK_DIR/vladsf" "{MODEL_NAME}" \
        --nzon "$BINARY_NZON" --kappa {KAPPA} \
        --fh {FH} --fhe {FHE} --abs-fraction-floor {ABS_FRACTION_FLOOR} \
        --stella-root "{STELLA_ROOT}" \
        --plot-dir "{PLOTS_DIR}" \
        {ABS_FRACTION_ARG} {OPAL_TABLES_ARG} 2>&1
    if [[ ! -s "$MODEL_WORK_DIR/vladsf/{MODEL_NAME}.1" ]]; then
        echo "ERROR: Failed to generate opacity tables"
        exit 1
    fi
    FILE_SIZE=$(wc -c < "{MODEL_NAME}.1")
    echo "Grey opacity tables generated (${FILE_SIZE} bytes each)."
fi

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
# Backup previous output before cleaning (preserves lightcurve for comparison)
echo "--- Backing up previous STELLA output ---"
BACKUP_DIR="$MODEL_WORK_DIR/backup_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$BACKUP_DIR"
for backup_dir in "$STELLA_ROOT/res" "$STRAD_DIR" "$MODEL_WORK_DIR"; do
    for ext in swd swd_rad tt flx lbol ph res; do
        if [[ -s "$backup_dir/$MODEL_NAME.$ext" ]]; then
            cp "$backup_dir/$MODEL_NAME.$ext" "$BACKUP_DIR/" 2>/dev/null || true
        fi
    done
done
if ls "$BACKUP_DIR"/* >/dev/null 2>&1; then
    echo "  Backed up previous output to $BACKUP_DIR"
else
    rmdir "$BACKUP_DIR" 2>/dev/null || true
    echo "  No previous output to back up"
fi

# Clean up files that STELLA must overwrite from both strad/ and model work directories
# Use absolute paths to ensure files are deleted from correct locations
# RESTART_STELLA: if "true", preserve .prf and .crv files so STELLA
# restarts from them (BEGRUN=.FALSE. in stradsep5tt.f). STELLA reads the
# restart state from .prf (profile) and .crv (curve) files, NOT from .res.
# The .res file is output-only (opened with status='unknown' after patch).
# Also preserve .swd and .tt so the lightcurve is continuous across restart.
echo "--- Cleaning up files that STELLA must overwrite ---"
mkdir -p "$STELLA_ROOT/res"
rm -f "$STELLA_ROOT/res/$MODEL_NAME.swd" "$STELLA_ROOT/res/$MODEL_NAME.swd_rad" "$STELLA_ROOT/res/$MODEL_NAME.tt" \
      "$STELLA_ROOT/res/$MODEL_NAME.lbol" "$STELLA_ROOT/res/$MODEL_NAME.ph" \
      "$STELLA_ROOT/res/$MODEL_NAME.flx"
for cleanup_dir in "$STRAD_DIR" "$MODEL_WORK_DIR"; do
    if [[ "$RESTART_STELLA" == "true" && -f "$cleanup_dir/$MODEL_NAME.prf" && -f "$cleanup_dir/$MODEL_NAME.crv" ]]; then
        echo "  RESTART MODE: preserving .prf, .crv, .swd, .tt for checkpoint restart"
    else
        rm -f "$cleanup_dir/$MODEL_NAME.res" "$cleanup_dir/$MODEL_NAME.bm" \
              "$cleanup_dir/$MODEL_NAME.balance" "$cleanup_dir/$MODEL_NAME.swd" \
              "$cleanup_dir/$MODEL_NAME.swd_rad" \
              "$cleanup_dir/$MODEL_NAME.tt" "$cleanup_dir/$MODEL_NAME.flx" \
              "$cleanup_dir/$MODEL_NAME.lbol" "$cleanup_dir/$MODEL_NAME.ph" \
              "$cleanup_dir/$MODEL_NAME.prf" \
              "$cleanup_dir/$MODEL_NAME.crv" 2>/dev/null || true
    fi
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
# Capture exit code without aborting (set -e would exit before .swd copy)
STELLA_EXIT=0
"$XSTELLA" < /dev/null 2>&1 | tee stella_run.log || STELLA_EXIT=$?
echo "STELLA run complete (exit=$STELLA_EXIT)."

# STELLA moves .swd and .tt to $HOMEStella/res/ at end of run
# Copy them back even if STELLA crashed — partial .swd output is still useful
# for CMFGEN/SEDONA snapshots at epochs before the crash.
cp "$STELLA_ROOT/res/$MODEL_NAME.swd" "$MODEL_WORK_DIR/" 2>/dev/null || true
cp "$STELLA_ROOT/res/$MODEL_NAME.swd_rad" "$MODEL_WORK_DIR/" 2>/dev/null || true
cp "$STELLA_ROOT/res/$MODEL_NAME.tt" "$MODEL_WORK_DIR/" 2>/dev/null || true
cp "$STELLA_ROOT/res/$MODEL_NAME.lbol" "$MODEL_WORK_DIR/" 2>/dev/null || true
cp "$STELLA_ROOT/res/$MODEL_NAME.ph" "$MODEL_WORK_DIR/" 2>/dev/null || true
cp "$STELLA_ROOT/res/$MODEL_NAME.flx" "$MODEL_WORK_DIR/" 2>/dev/null || true
# Also try strad/ in case STELLA behavior changes
cp "$STRAD_DIR/$MODEL_NAME.swd" "$MODEL_WORK_DIR/" 2>/dev/null || true
cp "$STRAD_DIR/$MODEL_NAME.swd_rad" "$MODEL_WORK_DIR/" 2>/dev/null || true
cp "$STRAD_DIR/$MODEL_NAME.tt" "$MODEL_WORK_DIR/" 2>/dev/null || true
cp "$STRAD_DIR/$MODEL_NAME.flx" "$MODEL_WORK_DIR/" 2>/dev/null || true
cp "$STRAD_DIR/$MODEL_NAME.lbol" "$MODEL_WORK_DIR/" 2>/dev/null || true

# Check what we got
if [[ ! -s "$MODEL_WORK_DIR/$MODEL_NAME.swd" ]]; then
    echo "ERROR: No usable .swd file found after STELLA run" >&2
    exit 1
fi
SWD_SIZE=$(stat -c%s "$MODEL_WORK_DIR/$MODEL_NAME.swd")
echo "  .swd file: $MODEL_WORK_DIR/$MODEL_NAME.swd ($SWD_SIZE bytes)"
if [[ -f "$MODEL_WORK_DIR/$MODEL_NAME.swd_rad" ]]; then
    SWD_RAD_SIZE=$(stat -c%s "$MODEL_WORK_DIR/$MODEL_NAME.swd_rad")
    echo "  .swd_rad file: $MODEL_WORK_DIR/$MODEL_NAME.swd_rad ($SWD_RAD_SIZE bytes)"
else
    echo "  WARNING: No .swd_rad file found after STELLA run"
fi
if [[ -f "$MODEL_WORK_DIR/$MODEL_NAME.tt" ]]; then
    echo "  .tt file: $MODEL_WORK_DIR/$MODEL_NAME.tt"
else
    echo "  WARNING: No .tt file found after STELLA run"
fi

if [[ $STELLA_EXIT -ne 0 ]]; then
    echo "WARNING: STELLA exited with code $STELLA_EXIT; checking for usable partial output." >&2
fi

# Count epochs and find last epoch in .swd
SWD_INFO=$(python - "$MODEL_WORK_DIR/$MODEL_NAME.swd" <<'PY'
import sys, re
epochs = set()
last = None
with open(sys.argv[1]) as fh:
    for line in fh:
        # Use regex to handle concatenated numbers (missing spaces in .swd)
        fields = re.findall(r'[-+]?\d*\.?\d+[eE][-+]?\d+|[-+]?\d*\.?\d+', line.strip())
        if len(fields) < 13:
            continue
        try:
            t = float(fields[0])
            km = int(float(fields[1]))
            if km == 1:
                epochs.add(round(t, 2))
            last = max(t, last) if last is not None else t
        except (ValueError, IndexError):
            pass
if last is None:
    print("0|None|0")
else:
    # Count post-explosion epochs (t > 0)
    post_exp = [t for t in epochs if t > 0]
    print(f"{len(epochs)}|{last:.4f}|{len(post_exp)}")
PY
)
SWD_EPOCH_COUNT=$(echo "$SWD_INFO" | cut -d'|' -f1)
LAST_SWD_DAY=$(echo "$SWD_INFO" | cut -d'|' -f2)
POST_EXP_EPOCHS=$(echo "$SWD_INFO" | cut -d'|' -f3)

echo "  .swd contains $SWD_EPOCH_COUNT epochs ($POST_EXP_EPOCHS post-explosion), last = $LAST_SWD_DAY d"

# Require at least 1 post-explosion epoch for useful output
if [[ "$POST_EXP_EPOCHS" -lt 1 ]]; then
    echo "ERROR: .swd has no post-explosion epochs (only initial condition); cannot produce useful output." >&2
    exit 1
fi

# Warn if last epoch is far from target
MIN_SWD_DAY=$(python - "$TARGET_DAYS" <<'PY'
import sys
end = float(sys.argv[1])
print(end - max(0.25, 0.01 * end))
PY
)
if ! python - "$LAST_SWD_DAY" "$MIN_SWD_DAY" <<'PY'
import sys
raise SystemExit(0 if float(sys.argv[1]) >= float(sys.argv[2]) else 1)
PY
then
    echo "WARNING: STELLA .swd last epoch ($LAST_SWD_DAY d) is short of target ($MIN_SWD_DAY d)." >&2
    echo "         Late snapshot days may be skipped by epoch_tolerance filter." >&2
fi

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
CMFGEN_DIR="$MODEL_DIR/cmfgen/{CMFGEN_DAY_LABEL}"
mkdir -p "$CMFGEN_DIR"
OUTPUT_HYDRO="$CMFGEN_DIR/SN_HYDRO_DATA"

python "$SCRIPTS_DIR/convert_stella_swd_to_cmfgen.py" \
    "$SWD_FILE" \
    "$OUTPUT_HYDRO" \
    --abn-file="$ABN_FILE" \
    --target-days="{CMFGEN_SNAPSHOT_DAYS}" \
    --kappa="$KAPPA" \
    --tau-max="{CMFGEN_TAU_MAX}" \
    --tau-min="{CMFGEN_TAU_MIN}" \
    --model-name="$MODEL_NAME" \
    --grid-thinning-threshold="{GRID_THINNING_THRESHOLD}" \
    {ENABLE_GRID_THINNING_ARG} \
    {ENABLE_ADAPTIVE_GRID_ARG} \
    --base-grid-target="{BASE_GRID_TARGET}" \
    --min-points="{MIN_POINTS}" \
    --max-radius-jump="{MAX_RADIUS_JUMP}" \
    --max-temp-jump="{MAX_TEMP_JUMP}" \
    --max-density-jump="{MAX_DENSITY_JUMP}" \
    --max-ne-jump="{MAX_NE_JUMP}" \
    --max-logtau-jump="{MAX_LOGTAU_JUMP}" \
    --max-points="{MAX_POINTS}" \
    --min-nd="{MIN_ND}" \
    {SMOOTH_VELOCITY_ARG} \
    --sigma-threshold="{SIGMA_THRESHOLD}" \
    --smooth-width="{SMOOTH_WIDTH}" \
    {COARSEN_NEGATIVE_SIGMA_ARG} \
    {XNI_ARG} \
    {INC_RAD_DECAYS_ARG} \
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

# Remove only the stale plots this stage regenerates; leave plots owned by
# other stages (presn/stella_input, sedona_*, superlite_*, cmfgen, etc.) intact
rm -f "$PLOTS_DIR/{MODEL_NAME}_"*evolution*.png \
      "$PLOTS_DIR/{MODEL_NAME}_"*radial_profiles*.png \
      "$PLOTS_DIR/{MODEL_NAME}_composition_evolution.png" \
      "$PLOTS_DIR/{MODEL_NAME}_composition_profile.png" \
      "$PLOTS_DIR/{MODEL_NAME}_stella_lightcurve.png" \
      "$PLOTS_DIR/{MODEL_NAME}_lightcurve_pysynphot.png" \
      "$PLOTS_DIR/stella_output_multiband.png" \
      "$PLOTS_DIR/nickel_evolution.png" \
      "$PLOTS_DIR/opacity_heatmap.png" "$PLOTS_DIR/opacity_vs_T.png" \
      "$PLOTS_DIR/opacity_vs_rho.png" "$PLOTS_DIR/opacity_ratio.png" \
      2>/dev/null || true

# Plot evolution (density, temperature, pressure, etc.) from .swd snapshots
python "$SCRIPTS_DIR/plot_stella_evolution.py" \
    "$SWD_FILE" \
    "$PLOTS_DIR" \
    --n-epochs 8 \
    2>&1 | tee -a stella_build.log || echo "WARNING: evolution plotting failed (non-fatal)"
echo "Evolution plots: $PLOTS_DIR/{MODEL_NAME}_evolution.png, $PLOTS_DIR/{MODEL_NAME}_radial_profiles.png"

# Plot radial evolution of mass fractions for each element
if [[ -f "$SWD_FILE" && -f "$SCRIPTS_DIR/plot_composition_evolution.py" ]]; then
    echo "Plotting composition evolution from .swd file (now includes composition data)..."
    python "$SCRIPTS_DIR/plot_composition_evolution.py" \
        "$ABN_FILE" \
        "$PLOTS_DIR" \
        --swd_file "$SWD_FILE" \
        --xni-file "$MODEL_WORK_DIR/eve/{MODEL_NAME}.xni" \
        --n-epochs 10 \
        2>&1 | tee -a stella_build.log || echo "WARNING: composition evolution plotting failed (non-fatal)"
    echo "Composition evolution plot: $PLOTS_DIR/{MODEL_NAME}_composition_evolution.png"
else
    echo "Skipping composition evolution plot (no .swd file or script not found)"
fi

# Lightcurve and .tt handling
if [[ -n "$TT_FILE" ]]; then
    cp "$TT_FILE" "$CMFGEN_DIR/"
    cp "$TT_FILE" "$OUTPUT_DIR/{MODEL_NAME}.tt"
    echo "Lightcurve copied: $TT_FILE -> $CMFGEN_DIR/"
    echo "Lightcurve copied: $TT_FILE -> $OUTPUT_DIR/{MODEL_NAME}.tt"

    # Plot lightcurve (with SN 2011dh reference overlay)
    PYTHONPATH="$SCRIPTS_DIR:${PYTHONPATH:-}" python "$SCRIPTS_DIR/plot_stella_lightcurve.py" \
        "$CMFGEN_DIR/{MODEL_NAME}.tt" \
        --output "$PLOTS_DIR/{MODEL_NAME}_stella_lightcurve.png" \
        --days-max "$TARGET_DAYS" \
        --csv "$OUTPUT_DIR/{MODEL_NAME}_lightcurve.csv" \
        --sn2011dh \
        2>&1 | tee -a stella_build.log || echo "WARNING: lightcurve plotting failed (non-fatal)"
    echo "Lightcurve plot: $PLOTS_DIR/{MODEL_NAME}_stella_lightcurve.png"
    echo "Lightcurve CSV: $OUTPUT_DIR/{MODEL_NAME}_lightcurve.csv"

    # Multiband (U,B,V,R,I) lightcurve vs SN 2011dh observations
    # Uses plot_model.py --step multiband (only needs .tt file)
    PYTHONPATH="$SCRIPTS_DIR:${PYTHONPATH:-}" python "$SCRIPTS_DIR/plot_model.py" \
        "$MODEL_DIR" --step multiband \
        2>&1 | tee -a stella_build.log || echo "WARNING: multiband plot failed (non-fatal)"
    echo "Multiband plot: $PLOTS_DIR/stella_output_multiband.png"

    # Two-event comparison plots: model vs SN_A obs (absolute timeline)
    # and vs SN_B obs (peak-aligned window)
    PYTHONPATH="$SCRIPTS_DIR:${PYTHONPATH:-}" python "$SCRIPTS_DIR/plot_model_vs_sn_ab.py" \
        "$MODEL_DIR" \
        2>&1 | tee -a stella_build.log || echo "WARNING: SN_A/SN_B comparison plots failed (non-fatal)"
    echo "Two-event plots: $PLOTS_DIR/{MODEL_NAME}_vs_sn_a.png, $PLOTS_DIR/{MODEL_NAME}_vs_sn_b.png"
else
    echo "Note: No .tt lightcurve file found for {MODEL_NAME}"
fi

# -------------------------------------------------------
# Step 8b: Extract complete multi-band lightcurve with pysynphot
# -------------------------------------------------------
echo ""
echo "--- Extracting complete pysynphot lightcurve ---"
FLX_FILE="$STRAD_DIR/{MODEL_NAME}.flx"
if [[ ! -f "$FLX_FILE" ]]; then
    FLX_FILE="$MODEL_WORK_DIR/strad/{MODEL_NAME}.flx"
fi
if [[ ! -f "$FLX_FILE" ]]; then
    FLX_FILE="$MODEL_WORK_DIR/{MODEL_NAME}.flx"
fi
if [[ ! -s "$FLX_FILE" ]]; then
    echo "WARNING: STELLA did not produce a usable .flx SED."
    echo "  Falling back to .lbol-based lightcurve extraction."

    # Fallback: extract lightcurve from .lbol (STELLA native bolometric + UBVRI)
    LBOL_FILE="$MODEL_WORK_DIR/{MODEL_NAME}.lbol"
    if [[ ! -f "$LBOL_FILE" ]]; then
        LBOL_FILE="$STELLA_ROOT/res/{MODEL_NAME}.lbol"
    fi
    if [[ -s "$LBOL_FILE" ]]; then
        echo "  Using .lbol file: $LBOL_FILE"
        PYTHONPATH="$SCRIPTS_DIR:${PYTHONPATH:-}" python "$SCRIPTS_DIR/extract_lightcurve_from_lbol.py" \
            "$LBOL_FILE" \
            --output "$OUTPUT_DIR/{MODEL_NAME}_lightcurve_complete.csv" \
            2>&1 | tee -a stella_build.log || echo "WARNING: .lbol extraction failed"
        echo "Complete lightcurve CSV (from .lbol): $OUTPUT_DIR/{MODEL_NAME}_lightcurve_complete.csv"
    else
        echo "ERROR: No .lbol file found either; cannot extract lightcurve." >&2
        exit 1
    fi
else
# Use waterboa python (has pysynphot in user .local); conda env may lack it
PYTHONPATH="$SCRIPTS_DIR:${PYTHONPATH:-}" \
"$WATERBOA_ROOT/bin/python3.12" \
    "$SCRIPTS_DIR/extract_lightcurve_pysynphot.py" \
    "$SWD_FILE" \
    "$FLX_FILE" \
    --tt-file "$TT_FILE" \
    --output "$OUTPUT_DIR/{MODEL_NAME}_lightcurve_complete.csv" \
    2>&1 | tee -a stella_build.log || echo "WARNING: pysynphot lightcurve extraction failed (non-fatal; NumPy 2.0 incompatibility in pysynphot)"
echo "Complete lightcurve CSV: $OUTPUT_DIR/{MODEL_NAME}_lightcurve_complete.csv"
fi
echo "Complete lightcurve CSV: $OUTPUT_DIR/{MODEL_NAME}_lightcurve_complete.csv"

# -------------------------------------------------------
# Step 9: Prepare SEDONA directory
# -------------------------------------------------------
echo ""
echo "--- Preparing SEDONA directory ---"

SEDONA_DIR="$MODEL_DIR/sedona/{CMFGEN_DAY_LABEL}"
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
    
    # Fix zero-abundance species and core radius safety
    python3 /ptmp/sbrennan/SN2024pba/scripts/fix_sedona_param.py "$SEDONA_DIR"
else
    echo "WARNING: param.lua template not found in sedona_base"
fi

echo "SEDONA directory prepared: $SEDONA_DIR"

# -------------------------------------------------------
# Step 10: Prepare CMFGEN directory
# -------------------------------------------------------
echo ""
echo "--- Preparing CMFGEN directory ---"

CMFGEN_DIR="$MODEL_DIR/cmfgen/{CMFGEN_DAY_LABEL}"
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

# -------------------------------------------------------
# Step 11: Generate nickel evolution plot (if Ni56 enabled)
# -------------------------------------------------------
echo ""
echo "--- Generating nickel evolution plot ---"

PLOTS_DIR="$MODEL_DIR/plots"
mkdir -p "$PLOTS_DIR"

XNI_FILE="$MODEL_WORK_DIR/eve/{MODEL_NAME}.xni"
ABN_DIR="$MODEL_WORK_DIR/modmake"
SWD_FILE="$MODEL_WORK_DIR/{MODEL_NAME}.swd"
NICKEL_PLOT="$PLOTS_DIR/nickel_evolution.png"

# When Ni56 is present, the decay-chain panels are merged into
# {MODEL_NAME}_composition_evolution.png by plot_composition_evolution.py
# (--xni-file), so the standalone nickel plot is skipped.  The standalone
# plot is still produced when Ni56 was not detected (uniform all-zero .xni)
# or when the composition script is unavailable.
NI56_DETECTED=0
if [[ -f "$XNI_FILE" ]]; then
    NI56_DETECTED=$(python - "$XNI_FILE" "$SCRIPTS_DIR" <<'PY' 2>/dev/null || echo 0
import sys
sys.path.insert(0, sys.argv[2])
from plot_nickel_evolution import read_xni
print(1 if read_xni(sys.argv[1]).max() > 1e-30 else 0)
PY
)
fi

if [[ "$NI56_DETECTED" == "1" ]]; then
    echo "Ni56 detected: decay-chain panels merged into composition_evolution.png; skipping standalone nickel plot"
elif [[ -f "$XNI_FILE" ]] && [[ -d "$ABN_DIR" ]]; then
    echo "Found .xni file and .abn directory, generating nickel evolution plot..."
    python "$SCRIPTS_DIR/plot_nickel_evolution.py" \
        "$MODEL_DIR" \
        --xni-file "$XNI_FILE" \
        --abn-dir "$ABN_DIR" \
        --swd-file "$SWD_FILE" \
        --n-epochs 10 \
        --output "$NICKEL_PLOT" || echo "WARNING: Nickel evolution plot generation failed"
else
    echo "Skipping nickel evolution plot (no .xni or .abn files found)"
fi

echo ""
echo "======================================"
echo "STELLA Job Complete: $MODEL_NAME"
echo "All input files ready for CMFGEN and SEDONA"
echo "======================================"
