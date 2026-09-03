#!/bin/bash
# Thin wrapper around run_model_grid.py
# Usage: ./scripts/run_model_grid.sh [run_model_grid.py args...]
# Examples:
#   ./scripts/run_model_grid.sh configs/X1.yaml
#   ./scripts/run_model_grid.sh --all
#   ./scripts/run_model_grid.sh configs/X1.yaml --no-redo-heracles --no-cmfgen

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

exec python3 "${SCRIPT_DIR}/run_model_grid.py" "$@"
