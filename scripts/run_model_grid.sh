#!/bin/bash
# Grid runner template for SN2024pba model grid calculations
# Usage: ./scripts/run_model_grid.sh configs/parameter_grid.yaml

set -euo pipefail

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

# Parse command line arguments
if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <grid_config.yaml>"
    echo "Example: $0 configs/parameter_grid.yaml"
    exit 1
fi

GRID_CONFIG_FILE="$1"

# Check if config file exists
if [[ ! -f "$GRID_CONFIG_FILE" ]]; then
    echo "ERROR: Grid config file not found: $GRID_CONFIG_FILE"
    exit 1
fi

echo "======================================================"
echo "SN2024pba Model Grid Run"
echo "======================================================"
echo "Grid config: $GRID_CONFIG_FILE"
echo ""

# This is a template for grid runs
# For actual grid runs, implement parameter grid parsing and parallel execution
# Example structure:
#   - Read grid config file
#   - Generate list of parameter combinations
#   - For each combination:
#     - Create temporary config file
#     - Run single model workflow
#     - Collect results
#   - Optionally run in parallel using GNU parallel or Slurm array jobs

echo "Grid runner template - implement parameter grid parsing here"
echo ""
echo "Example grid config structure:"
echo "  grid_parameters:"
echo "    - name: M_ejecta_A"
echo "      values: [9.0, 9.5, 10.0]"
echo "    - name: M_wind"
echo "      values: [0.05, 0.1, 0.15]"
echo ""
echo "This will create 3 x 3 = 9 model runs"
echo ""
echo "To implement:"
echo "  1. Parse grid config YAML"
echo "  2. Generate parameter combinations"
echo "  3. Create config files for each combination"
echo "  4. Run models (serial or parallel)"
echo "  5. Collect and organize results"
echo ""
echo "======================================================"
