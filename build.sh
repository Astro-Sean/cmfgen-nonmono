#!/bin/bash
# Build all CMFGEN executables (optimized), including cmf_flux.exe used by OBS.
# KEEP_EXE keeps the previous executables live while objects/libraries rebuild;
# production links publish via .tmp + mv. Fails immediately on a make error.
set -euo pipefail

cd "$(dirname "$0")"
trap 'rm -f exe/*.tmp' EXIT

make KEEP_EXE=1 clean
make DEBUG=nodebug all

for exe in exe/cmfgen_dev.exe exe/cmf_flux.exe exe/obs_frame.exe exe/plt_steq.exe; do
    [[ -x "${exe}" ]] || { echo "ERROR: missing build product: ${exe}" >&2; exit 1; }
done

# Compile-time verification: the JFNK solver/interface suites compile
# the production solver unchanged and exit non-zero on any failed
# check, so a broken solver fails the build like a compile error.
tests/f90/build_and_run.sh

echo "Build succeeded."
