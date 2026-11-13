#!/bin/bash
# Build all CMFGEN executables (optimized), including cmf_flux.exe used by OBS.
# KEEP_EXE keeps the previous executables live while objects/libraries rebuild;
# production links publish via .tmp + mv. Fails immediately on a make error.
set -euo pipefail

cd "$(dirname "$0")"

# Serialize builds in this tree: two concurrent runs race object
# files and archive updates (ar reads a .o the other run is still
# writing or has just cleaned), which fails with "file truncated".
# flock is broken in this environment, so use a mkdir lock - mkdir
# is atomic on GPFS.  A lock older than 2h is treated as stale (the
# holder was killed) and reclaimed; waiting gives up after 1h.
LOCK_DIR=".build.lock.d"
waited=0
while ! mkdir "$LOCK_DIR" 2>/dev/null; do
    if find "$LOCK_DIR" -maxdepth 0 -mmin +120 | grep -q .; then
        rmdir "$LOCK_DIR" 2>/dev/null || true
        continue
    fi
    waited=$((waited + 15))
    if [[ $waited -ge 3600 ]]; then
        echo "ERROR: build lock $LOCK_DIR held longer than 3600s" >&2
        exit 1
    fi
    echo "build.sh: waiting for build lock $LOCK_DIR" >&2
    sleep 15
done
trap 'rmdir "$LOCK_DIR" 2>/dev/null; rm -f exe/*.tmp' EXIT

# Makefiles never create the output dirs; a fresh clone fails without them.
mkdir -p mod lib exe mod_dbg lib_dbg exe_dbg

make KEEP_EXE=1 clean
make DEBUG=nodebug all

for exe in exe/cmfgen_dev.exe exe/cmf_flux.exe exe/obs_frame.exe exe/plt_steq.exe; do
    [[ -x "${exe}" ]] || { echo "ERROR: missing build product: ${exe}" >&2; exit 1; }
done

# Compile-time verification: the JFNK solver/interface suites compile
# the production solver unchanged and exit non-zero on any failed
# check, so a broken solver fails the build like a compile error.
tests/f90/build_and_run.sh

# Parser/unit regression suite: each python/test_*.py compiles the
# production sources it exercises, so a parser or guard regression
# fails the build here rather than in a model run.
python3 unit_tests/run_all_tests.py

echo "Build succeeded."
