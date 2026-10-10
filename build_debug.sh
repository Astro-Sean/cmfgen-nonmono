#!/bin/bash
# Build CMFGEN dev7 executable (debug).
# Fails immediately if any make step returns non-zero.
set -euo pipefail

cd "$(dirname "$0")"

make dbg_clean
make DEBUG=debug cmfgen

echo "Debug build succeeded."
