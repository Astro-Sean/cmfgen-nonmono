#!/bin/bash
# Build the opal_opacity Python extension module (opal_opacity.so) from
# SNEC's opal_opacity.F90 using ifort + f2py.
#
# Run this once on the cluster login/compute node before using --opal-tables.
# The .so is placed in the same directory as this script (scripts/).
#
# Usage:
#   bash scripts/build_opal_module.sh [SNEC_DIR]
#
# SNEC_DIR defaults to /ptmp/sbrennan/SNEC/SNEC

set -e

SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SNEC_DIR="${1:-/ptmp/sbrennan/SNEC/SNEC}"
OPAL_SRC="$SNEC_DIR/src/opal_opacity.F90"

if [ ! -f "$OPAL_SRC" ]; then
    echo "ERROR: $OPAL_SRC not found. Set SNEC_DIR correctly."
    exit 1
fi

# Locate Python 3.12 (the system install with numpy 1.26)
PYTHON=$(which python3.12 2>/dev/null || which python3)
echo "Using Python: $PYTHON ($($PYTHON --version))"

NP_INC=$($PYTHON -c "import numpy; print(numpy.get_include())")
F2PY_SRC=$($PYTHON -c "import numpy.f2py, os; print(os.path.join(os.path.dirname(numpy.f2py.__file__), 'src'))")
PY_INC=$($PYTHON -c "import sysconfig; print(sysconfig.get_path('include'))")
PY_LIB=$($PYTHON -c "import sysconfig; print(sysconfig.get_config_var('LIBDIR'))")
PY_VER=$($PYTHON -c "import sysconfig; print(sysconfig.get_config_var('LDVERSION'))")

cd "$SCRIPTS_DIR"

echo "Generating f2py wrappers..."
$PYTHON -m numpy.f2py -m opal_opacity opal_wrapper.f90 "$OPAL_SRC" 2>/dev/null

echo "Compiling..."
ifort -O2 -fPIC -free -diag-disable=10448 -c "$OPAL_SRC" -o opal_opacity_core.o
ifort -O2 -fPIC -free -diag-disable=10448 -c opal_wrapper.f90 -o opal_wrapper.o
ifort -O2 -fPIC -diag-disable=10448      -c opal_opacity-f2pywrappers.f -o opal_f2py_wrappers.o
gcc -O2 -fPIC \
    -c opal_opacitymodule.c \
    -I"$NP_INC" -I"$F2PY_SRC" -I"$PY_INC" \
    -o opal_opacitymodule.o
gcc -O2 -fPIC \
    -c "$F2PY_SRC/fortranobject.c" \
    -I"$NP_INC" -I"$F2PY_SRC" -I"$PY_INC" \
    -o fortranobject.o

echo "Linking..."
ifort -shared -fPIC -diag-disable=10448 \
    opal_opacity_core.o opal_wrapper.o opal_f2py_wrappers.o \
    opal_opacitymodule.o fortranobject.o \
    -L"$PY_LIB" -lpython${PY_VER} \
    -Wl,-rpath,"$PY_LIB" \
    -o opal_opacity.so

# Clean up intermediate files
rm -f opal_opacity_core.o opal_wrapper.o opal_f2py_wrappers.o \
      opal_opacitymodule.o fortranobject.o \
      opal_opacitymodule.c opal_opacity-f2pywrappers.f

echo ""
echo "SUCCESS: opal_opacity.so built in $SCRIPTS_DIR"
ls -lh "$SCRIPTS_DIR/opal_opacity.so"
