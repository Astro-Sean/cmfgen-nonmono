#!/bin/bash
# Build CMFGEN dev7 executable (debug).
# Fails immediately if any make step returns non-zero.
set -euo pipefail

cd "$(dirname "$0")"

# Makefiles never create the output dirs; a fresh clone fails without them.
mkdir -p mod_dbg lib_dbg exe_dbg

make dbg_clean
make DEBUG=debug cmfgen

echo "Debug build succeeded."
