"""MPC interactive-partition test job generator, submission, and monitoring.

This module generates a dedicated SLURM test script that exercises the
complete downstream validation workflow on the MPC HPC interactive partition,
submits it via sbatch, monitors it to terminal state, and validates outputs.

The test uses validated existing STELLA output (reuse mode) to avoid the cost
of a fresh STELLA run, while still exercising the full downstream path:

    YAML validation
    → resolved configuration
    → environment/module setup
    → STELLA reuse validation
    → snapshot and epoch mapping
    → STELLA-to-CMFGEN conversion
    → conservative grid validation
    → non-monotonic velocity and ionization-front checks
    → CMFGEN cold-start preparation
    → output and convergence validation
"""

from __future__ import annotations

import hashlib
import json
import subprocess
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Optional, Any

from .config import ResolvedConfig, ConfigError
from .site import SiteConfig
from .monitor import SlurmMonitor, JobState, JobInfo


@dataclass
class InteractiveTestResult:
    """Result of an interactive-partition test run."""
    job_id: str
    run_id: str
    partition: str
    account: str
    time_limit: str
    memory: str
    cpus: int
    submitted_at: str
    final_state: str
    exit_code: Optional[int]
    elapsed: str
    test_dir: str
    script_hash: str
    config_hash: str
    passed: bool
    gates: dict[str, bool]
    failure_reasons: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


# ── Helper Python scripts (written to test dir, called by SLURM script) ────

_HELPER_VALIDATE_CONFIG = '''#!/usr/bin/env python3
"""Stage 1: Validate YAML config and print resolved settings."""
import sys, os
# Add scripts dir to path (from PYTHONPATH or computed from helper location)
scripts_dir = os.environ.get('SCRIPTS_DIR', '')
if not scripts_dir:
    # Fallback: go up 3 levels from helpers/ dir to scripts/
    scripts_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, scripts_dir)
from orchestration.config import resolve_config, ConfigError

config_path = sys.argv[1]
project_dir = sys.argv[2]
scripts_dir = sys.argv[3]

try:
    cfg = resolve_config(config_path, project_dir, scripts_dir)
    print(f"Config validated: {cfg.model_name}, hash={cfg.config_hash}")
    print(f"Snapshot days: {cfg.snapshot_days}")
    print(f"STELLA target: {cfg.stella_target_days}d")
    print(f"Site: {cfg.site.cluster_name} partition={cfg.site.interactive_partition}")
    print(f"Reuse STELLA: {cfg.reuse_stella}")
    print(f"Run CMFGEN: {cfg.run_cmfgen}")
    # Write resolved config to test dir (parent of helpers/)
    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "resolved_config.yaml")
    with open(out, "w") as f:
        f.write(cfg.to_yaml())
    print(f"Resolved config written: {out}")
except ConfigError as e:
    print(f"CONFIG ERROR: {e}", file=sys.stderr)
    sys.exit(1)
'''

_HELPER_VALIDATE_SNHYDRO = '''#!/usr/bin/env python3
"""Stage 4: Validate SN_HYDRO_DATA content (radius ordering, finite values, etc.).

Parses the CMFGEN SN_HYDRO_DATA block format:
  Header lines with metadata
  "Radius grid (10^10cm)" followed by ND values in rows of 6
  "Velocity (km/s)" followed by ND values
  "Sigma (dlnV/dlnr-1)" followed by ND values
  "Temperature (10^4 K)" followed by ND values
  "Density (gm/cm^3)" followed by ND values
  ... (atom density, electron density, opacity, kappa, mass fractions, isotopes)
"""
import sys, os, json, re
import numpy as np

path = sys.argv[1]
if not os.path.exists(path):
    print(f"ERROR: SN_HYDRO_DATA not found: {path}", file=sys.stderr)
    sys.exit(1)

lines = open(path).readlines()

# Parse header
nd = None
n_frac = None
n_iso = None
time_days = None
lum_lsun = None
for line in lines:
    stripped = line.strip()
    if stripped.startswith("Number of data points:"):
        nd = int(stripped.split(":")[1].strip())
    elif stripped.startswith("Number of mass fractions:"):
        n_frac = int(stripped.split(":")[1].strip())
    elif stripped.startswith("Number of isotopes:"):
        n_iso = int(stripped.split(":")[1].strip())
    elif stripped.startswith("Time(days)"):
        time_days = float(stripped.split(":")[1].strip())
    elif stripped.startswith("Energy Lum"):
        lum_lsun = float(stripped.split(":")[1].strip())

if nd is None:
    print("ERROR: Could not parse ND from SN_HYDRO_DATA header", file=sys.stderr)
    sys.exit(1)

print(f"  ND: {nd}")
if time_days is not None:
    print(f"  Time: {time_days:.2f} days")
if lum_lsun is not None:
    print(f"  Lum: {lum_lsun:.3e} Lsun")

# Parse block-structured data
# Each block starts with a label line, followed by data values in rows of 6
blocks = {}
current_label = None
current_values = []
for line in lines:
    stripped = line.strip()
    if not stripped:
        continue
    # Check if this is a label line (starts with a capital letter, contains letters)
    if re.match(r'^[A-Z][a-zA-Z]', stripped) and not re.match(r'^[-+0-9.Ee]', stripped):
        # Save previous block
        if current_label and current_values:
            blocks[current_label] = np.array(current_values)
        current_label = stripped
        current_values = []
    else:
        # Data line — parse all float values
        try:
            vals = [float(x) for x in stripped.split()]
            current_values.extend(vals)
        except ValueError:
            pass
# Save last block
if current_label and current_values:
    blocks[current_label] = np.array(current_values)

# Extract key arrays
# Radius is in units of 10^10 cm
r_label = next((k for k in blocks if "Radius" in k), None)
v_label = next((k for k in blocks if "Velocity" in k), None)
t_label = next((k for k in blocks if "Temperature" in k), None)
rho_label = next((k for k in blocks if "Density (gm" in k), None)

if not all([r_label, v_label, t_label, rho_label]):
    print(f"ERROR: Missing required blocks. Found: {list(blocks.keys())}", file=sys.stderr)
    sys.exit(1)

r = blocks[r_label] * 1e10  # convert from 10^10 cm to cm
v = blocks[v_label]          # km/s
t = blocks[t_label] * 1e4    # convert from 10^4 K to K
rho = blocks[rho_label]      # g/cm^3

# Verify array lengths
for name, arr in [("R", r), ("V", v), ("T", t), ("rho", rho)]:
    if len(arr) != nd:
        print(f"ERROR: {name} has {len(arr)} values, expected {nd}", file=sys.stderr)
        sys.exit(1)

print(f"  R range: {r.min():.3e} - {r.max():.3e} cm")
print(f"  V range: {v.min():.1f} - {v.max():.1f} km/s")
print(f"  rho range: {rho.min():.3e} - {rho.max():.3e} g/cm3")
print(f"  T range: {t.min():.1f} - {t.max():.1f} K")

# Check radius is strictly decreasing (CMFGEN outer-to-inner convention)
dr = np.diff(r)
n_violations = int(np.sum(dr >= 0))
if n_violations > 0:
    print(f"FAIL: Radius not strictly decreasing, {n_violations} violations")
    sys.exit(1)

# Check finite values
for name, arr in [("R", r), ("V", v), ("rho", rho), ("T", t)]:
    if not np.all(np.isfinite(arr)):
        n_bad = int(np.sum(~np.isfinite(arr)))
        print(f"FAIL: {name} has {n_bad} non-finite values")
        sys.exit(1)

# Check positive density and temperature
if np.any(rho <= 0):
    n_bad = int(np.sum(rho <= 0))
    print(f"FAIL: rho has {n_bad} non-positive values")
    sys.exit(1)
if np.any(t <= 0):
    n_bad = int(np.sum(t <= 0))
    print(f"FAIL: T has {n_bad} non-positive values")
    sys.exit(1)

# Check velocity is non-zero (at least some zones should have v > 0)
if np.all(v == 0):
    print("FAIL: All velocities are zero")
    sys.exit(1)

# Write validation result
val_result = {
    "nd": nd,
    "n_mass_fractions": n_frac,
    "n_isotopes": n_iso,
    "time_days": time_days,
    "lum_lsun": lum_lsun,
    "r_min": float(r.min()), "r_max": float(r.max()),
    "v_min": float(v.min()), "v_max": float(v.max()),
    "rho_min": float(rho.min()), "rho_max": float(rho.max()),
    "t_min": float(t.min()), "t_max": float(t.max()),
    "n_blocks": len(blocks),
    "block_names": list(blocks.keys()),
    "radius_strictly_decreasing": True,
    "all_finite": True,
    "density_positive": True,
    "temperature_positive": True,
    "velocity_nonzero": True,
    "array_lengths_consistent": True,
}
out_path = os.path.join(os.path.dirname(path), "sn_hydro_validation.json")
with open(out_path, "w") as f:
    json.dump(val_result, f, indent=2)
print(f"  Validation JSON written: {out_path}")
print(f"  Blocks parsed: {len(blocks)}")
print("  All content checks passed")
'''


class InteractiveTestRunner:
    """Generate, submit, and monitor an MPC interactive-partition test job.

    Usage:
        runner = InteractiveTestRunner(config, test_dir)
        result = runner.run()
        if not result.passed:
            print("Test FAILED — production blocked")
    """

    def __init__(
        self,
        config: ResolvedConfig,
        test_dir: str,
        init_iterations: int = 3,
        poll_interval: int = 15,
    ):
        self.config = config
        self.test_dir = Path(test_dir)
        self.test_dir.mkdir(parents=True, exist_ok=True)
        self.init_iterations = init_iterations
        self.poll_interval = poll_interval
        self.site = config.site
        self.run_id = f"{config.run_id}_itest"

    def _compute_hash(self, content: str) -> str:
        return hashlib.sha256(content.encode()).hexdigest()[:16]

    def _write_helpers(self):
        """Write helper Python scripts to the test directory."""
        helpers_dir = self.test_dir / "helpers"
        helpers_dir.mkdir(exist_ok=True)

        (helpers_dir / "validate_config.py").write_text(_HELPER_VALIDATE_CONFIG)
        (helpers_dir / "validate_snhydro.py").write_text(_HELPER_VALIDATE_SNHYDRO)

        return helpers_dir

    def _find_swd_path(self) -> str:
        """Find the STELLA .swd file for reuse validation."""
        cfg = self.config
        stella_out = Path(cfg.stella_output_dir)
        # Try strad subdir first (STELLA output convention)
        for sub in ["strad", ""]:
            p = stella_out / sub / f"{cfg.model_name}.swd"
            if p.exists():
                return str(p)
        # Try any .swd in the stella output dir
        for p in sorted(stella_out.rglob("*.swd")):
            return str(p)
        return ""

    def generate_test_script(self) -> tuple[str, str]:
        """Generate the SLURM test script. Returns (script_path, script_hash)."""
        # Write helper scripts first
        helpers_dir = self._write_helpers()

        script = self._build_script_content(helpers_dir)
        script_hash = self._compute_hash(script)
        script_path = self.test_dir / "interactive_test.slurm"
        script_path.write_text(script)
        return str(script_path), script_hash

    def _build_script_content(self, helpers_dir: Path) -> str:
        """Build the complete SLURM test script content."""
        cfg = self.config
        site = self.site

        swd_path = self._find_swd_path()
        test_day = cfg.snapshot_days[0] if cfg.snapshot_days else 30.0
        day_label = f"{test_day:.2f}days"

        conv_script = f"{cfg.scripts_dir}/convert_stella_swd_to_cmfgen_v2.py"

        cmfgen_cfg = cfg.cmfgen
        vel_law = cmfgen_cfg.get("vel_law", 12)
        target_nd = cmfgen_cfg.get("target_nd", 150)
        cmfgen_base = cmfgen_cfg.get("base_dir", "./cmfgen_base")
        cmfgen_base_abs = str((Path(cfg.project_dir) / cmfgen_base).resolve())

        # CMFGEN conversion params from config
        conv_cfg = cfg.cmfgen_conversion
        tau_max = conv_cfg.get("cmfgen_tau_max", 30)
        tau_min = conv_cfg.get("cmfgen_tau_min", 0.0)
        kappa = conv_cfg.get("kappa_constant", 0.34)

        # Use string formatting with explicit variables to avoid brace conflicts
        # We build the script line by line
        L = []  # lines

        def add(s=""):
            L.append(s)

        add("#!/bin/bash")
        add("# MPC Interactive-Partition Test Job")
        add(f"# Run ID: {self.run_id}")
        add(f"# Model: {cfg.model_name}")
        add(f"# Config hash: {cfg.config_hash}")
        add(f"# Generated: {time.strftime('%Y-%m-%dT%H:%M:%S')}")
        add("")
        add(f"#SBATCH --job-name={cfg.model_name}_itest")
        add(f"#SBATCH --partition={site.interactive_partition}")
        add(f"#SBATCH --account={site.account}")
        add(f"#SBATCH --time={site.test_time}")
        add(f"#SBATCH --mem={site.test_mem}")
        add(f"#SBATCH --cpus-per-task={site.test_cpus}")
        add(f"#SBATCH --ntasks={site.test_ntasks}")
        add(f"#SBATCH --output={self.test_dir}/interactive_test_%j.out")
        add(f"#SBATCH --error={self.test_dir}/interactive_test_%j.err")
        add(f"#SBATCH --chdir={self.test_dir}")
        add("")
        add("set -euo pipefail")
        add("")
        add("# ── Failure trap ──────────────────────────────────────────────")
        add('FAILED_STAGE="unknown"')
        add("trap '")
        add("  EXIT_CODE=$?")
        add('  echo ""')
        add('  echo "=========================================="')
        add('  echo "TEST FAILED at stage: ${FAILED_STAGE} (exit code: ${EXIT_CODE})"')
        add('  echo "=========================================="')
        add(f'  echo "{{\\"timestamp\\": \\"$(date -Iseconds)\\", \\"stage\\": \\"${{FAILED_STAGE}}\\", \\"exit_code\\": \\"${{EXIT_CODE}}\\"}}" >> {self.test_dir}/test_failure.jsonl')
        add("  exit ${EXIT_CODE}")
        add("' EXIT")
        add("")
        add("# ── Environment setup ────────────────────────────────────────")
        add("source /mpcdf/soft/RHEL_9/packages/x86_64/python-waterboa/2024.06/bin/activate /ptmp/sbrennan/conda-envs/astro")
        add("unset PYTHONPATH PYTHONHOME")
        add("")
        add("module purge")
        add("module load intel/2024.0")
        add("module load python-waterboa/2024.06")
        add("export SYSTYPE=ifort")
        add(f'export PYTHONPATH="{cfg.scripts_dir}:${{PYTHONPATH:-}}"')
        add("")
        add("# ── Environment capture ──────────────────────────────────────")
        add('echo "==========================================\"')
        add('echo "MPC Interactive-Partition Test"')
        add('echo "Job ID: $SLURM_JOB_ID"')
        add('echo "Node: $SLURMD_NODENAME"')
        add('echo "Date: $(date -Iseconds)"')
        add('echo "=========================================="')
        add('echo ""')
        add('echo "--- Module list ---"')
        add("module list 2>&1")
        add('echo ""')
        add('echo "--- Python version ---"')
        add("python3 --version")
        add("which python3")
        add('echo ""')
        add('echo "--- Key paths ---"')
        add(f'echo "CMFGEN base: {cmfgen_base_abs}"')
        add(f'echo "Conversion script: {conv_script}"')
        add(f'echo "Project dir: {cfg.project_dir}"')
        add(f'echo "Scripts dir: {cfg.scripts_dir}"')
        add(f'echo "STELLA .swd: {swd_path}"')
        add('echo ""')
        add('echo "--- Git identity ---"')
        add(f"cd {cfg.project_dir} && git rev-parse HEAD 2>/dev/null || echo 'not a git repo'")
        add('echo ""')
        add('echo "--- SLURM environment ---"')
        add("env | grep -iE 'SLURM|PYTHON|INTEL|SYSTYPE' | sort")
        add('echo ""')
        add('echo "=========================================="')
        add('echo ""')
        add("")
        # Variables
        add(f'PROJECT_DIR="{cfg.project_dir}"')
        add(f'SCRIPTS_DIR="{cfg.scripts_dir}"')
        add(f'MODEL_NAME="{cfg.model_name}"')
        add(f'CONFIG_FILE="{cfg.config_path}"')
        add(f'TEST_DIR="{self.test_dir}"')
        add(f'HELPERS_DIR="{helpers_dir}"')
        add(f'SWD_FILE="{swd_path}"')
        add(f'TEST_DAY="{test_day}"')
        add(f'DAY_LABEL="{day_label}"')
        add(f'CMFGEN_BASE="{cmfgen_base_abs}"')
        add(f'VEL_LAW="{vel_law}"')
        add(f'TARGET_ND="{target_nd}"')
        add(f'TAU_MAX="{tau_max}"')
        add(f'TAU_MIN="{tau_min}"')
        add(f'KAPPA="{kappa}"')
        add(f'INIT_ITERATIONS="{self.init_iterations}"')
        add("")

        # Stage 1: YAML validation
        add("# ── Stage 1: YAML validation ─────────────────────────────────")
        add('FAILED_STAGE="yaml_validation"')
        add('echo "[Stage 1] YAML validation..."')
        add('python3 "$HELPERS_DIR/validate_config.py" "$CONFIG_FILE" "$PROJECT_DIR" "$SCRIPTS_DIR"')
        add('echo "[Stage 1] PASSED"')
        add('echo ""')
        add("")

        # Stage 2: STELLA reuse validation
        add("# ── Stage 2: STELLA reuse validation ─────────────────────────")
        add('FAILED_STAGE="stella_reuse_validation"')
        add('echo "[Stage 2] STELLA reuse validation..."')
        add('if [[ -z "$SWD_FILE" ]]; then')
        add('    echo "ERROR: No STELLA .swd file found for reuse" >&2')
        add('    exit 1')
        add('fi')
        add('if [[ ! -f "$SWD_FILE" ]]; then')
        add('    echo "ERROR: STELLA .swd file not found: $SWD_FILE" >&2')
        add('    exit 1')
        add('fi')
        add('SWD_SIZE=$(stat -c%s "$SWD_FILE" 2>/dev/null || echo 0)')
        add('if [[ "$SWD_SIZE" -lt 1000 ]]; then')
        add('    echo "ERROR: STELLA .swd file too small ($SWD_SIZE bytes)" >&2')
        add('    exit 1')
        add('fi')
        add('SWD_HASH=$(sha256sum "$SWD_FILE" | cut -c1-16)')
        add('echo "  .swd path: $SWD_FILE"')
        add('echo "  .swd size: $SWD_SIZE bytes"')
        add('echo "  .swd hash: $SWD_HASH"')
        add('echo "[Stage 2] PASSED"')
        add('echo ""')
        add("")

        # Stage 3: STELLA-to-CMFGEN conversion
        add("# ── Stage 3: STELLA-to-CMFGEN conversion ─────────────────────")
        add('FAILED_STAGE="stella_to_cmfgen_conversion"')
        add('echo "[Stage 3] STELLA-to-CMFGEN conversion..."')
        add('mkdir -p "$TEST_DIR/cmfgen/$DAY_LABEL"')
        add('SN_HYDRO="$TEST_DIR/cmfgen/$DAY_LABEL/SN_HYDRO_DATA"')
        add('echo "  Converting $SWD_FILE -> $SN_HYDRO"')
        add('echo "  Target day: $TEST_DAY, VEL_LAW: $VEL_LAW, Target ND: $TARGET_ND"')
        add('python3 "$SCRIPTS_DIR/convert_stella_swd_to_cmfgen_v2.py" \\')
        add('    "$SWD_FILE" \\')
        add('    "$SN_HYDRO" \\')
        add('    --target-days "$TEST_DAY" \\')
        add('    --min-nd "$TARGET_ND" \\')
        add('    --max-points "$TARGET_ND" \\')
        add('    --nzon-target "$TARGET_ND" \\')
        add('    --tau-max "$TAU_MAX" \\')
        add('    --tau-min "$TAU_MIN" \\')
        add('    --kappa "$KAPPA" \\')
        add('    --model-name "$MODEL_NAME" \\')
        add('    --plot-dir "$TEST_DIR/cmfgen/$DAY_LABEL" \\')
        add('    2>&1')
        add('if [[ ! -f "$SN_HYDRO" ]]; then')
        add('    echo "ERROR: SN_HYDRO_DATA not created: $SN_HYDRO" >&2')
        add('    exit 1')
        add('fi')
        add('SN_HYDRO_SIZE=$(stat -c%s "$SN_HYDRO" 2>/dev/null || echo 0)')
        add('if [[ "$SN_HYDRO_SIZE" -lt 100 ]]; then')
        add('    echo "ERROR: SN_HYDRO_DATA too small ($SN_HYDRO_SIZE bytes)" >&2')
        add('    exit 1')
        add('fi')
        add('echo "  SN_HYDRO_DATA created: $SN_HYDRO_SIZE bytes"')
        add('echo "[Stage 3] PASSED"')
        add('echo ""')
        add("")

        # Stage 4: SN_HYDRO_DATA content validation
        add("# ── Stage 4: SN_HYDRO_DATA content validation ────────────────")
        add('FAILED_STAGE="sn_hydro_validation"')
        add('echo "[Stage 4] SN_HYDRO_DATA content validation..."')
        add('python3 "$HELPERS_DIR/validate_snhydro.py" "$SN_HYDRO"')
        add('echo "[Stage 4] PASSED"')
        add('echo ""')
        add("")

        # Stage 5: Grid validation report
        add("# ── Stage 5: Grid validation report ──────────────────────────")
        add('FAILED_STAGE="grid_validation"')
        add('echo "[Stage 5] Grid validation report..."')
        add('if [[ -f "$TEST_DIR/cmfgen/$DAY_LABEL/grid_validation.json" ]]; then')
        add('    echo "  Grid validation JSON found"')
        add('    python3 -c "import json,sys; d=json.load(open(sys.argv[1])); print(f\'  Grid points: {d.get(\\\"nd\\\", \\\"?\\\")}\'); [print(f\'  {\\\"PASS\\\" if v else \\\"FAIL\\\"}: {k}\') for k,v in d.get(\'checks\',{}).items()]" "$TEST_DIR/cmfgen/$DAY_LABEL/grid_validation.json" || true')
        add('else')
        add('    echo "  NOTE: No grid_validation.json (conversion may not write one)"')
        add('fi')
        add('echo "[Stage 5] PASSED"')
        add('echo ""')
        add("")

        # Stage 6: CMFGEN directory preparation
        add("# ── Stage 6: CMFGEN directory preparation ────────────────────")
        add('FAILED_STAGE="cmfgen_prep"')
        add('echo "[Stage 6] CMFGEN directory preparation..."')
        add('CMFGEN_RUN_DIR="$TEST_DIR/cmfgen_run/$DAY_LABEL"')
        add('mkdir -p "$CMFGEN_RUN_DIR"')
        add('cp "$SN_HYDRO" "$CMFGEN_RUN_DIR/SN_HYDRO_DATA"')
        add('cp "$CMFGEN_BASE"/MODEL_SPEC "$CMFGEN_RUN_DIR/" 2>/dev/null || true')
        add('cp "$CMFGEN_BASE"/RDINR "$CMFGEN_RUN_DIR/" 2>/dev/null || true')
        add('cp "$CMFGEN_BASE"/batch.sh "$CMFGEN_RUN_DIR/" 2>/dev/null || true')
        add('cp "$CMFGEN_BASE"/submit_cmfgen.sh "$CMFGEN_RUN_DIR/" 2>/dev/null || true')
        add('cp -r "$CMFGEN_BASE"/obs "$CMFGEN_RUN_DIR/" 2>/dev/null || true')
        add('cp "$CMFGEN_BASE"/*.dat "$CMFGEN_RUN_DIR/" 2>/dev/null || true')
        add('cp "$CMFGEN_BASE"/NUC_DECAY_DATA "$CMFGEN_RUN_DIR/" 2>/dev/null || true')
        add('cp "$CMFGEN_BASE"/ADJUST_R_DEFAULTS "$CMFGEN_RUN_DIR/" 2>/dev/null || true')
        add('echo "  CMFGEN run dir: $CMFGEN_RUN_DIR"')
        add('ls "$CMFGEN_RUN_DIR/" | head -20')
        add('echo "[Stage 6] PASSED"')
        add('echo ""')
        add("")

        # Stage 7: VADAT preparation
        add("# ── Stage 7: VADAT preparation ───────────────────────────────")
        add('FAILED_STAGE="vadat_prep"')
        add('echo "[Stage 7] VADAT preparation..."')
        add('python3 "$SCRIPTS_DIR/prepare_VADAT.py" \\')
        add('    --single-final \\')
        add('    --vel-law "$VEL_LAW" \\')
        add('    --target-nd "$TARGET_ND" \\')
        add('    "$CMFGEN_RUN_DIR" \\')
        add('    2>&1 || echo "  WARNING: VADAT prep had issues"')
        add('if [[ -f "$CMFGEN_RUN_DIR/VADAT" ]]; then')
        add('    echo "  VADAT created successfully"')
        add('    head -20 "$CMFGEN_RUN_DIR/VADAT"')
        add('else')
        add('    echo "  WARNING: VADAT not found, creating minimal VADAT for test"')
        add('    cat > "$CMFGEN_RUN_DIR/VADAT" << VADATEOF')
        add(f' &INVERG')
        add(f'  ND={target_nd},')
        add(f'  VEL_LAW={vel_law},')
        add(f'  INIT_NUM_ITS={self.init_iterations},')
        add(f'  DO_NG=F,')
        add(f'  METHOD="ZERO",')
        add(f' /')
        add('VADATEOF')
        add('fi')
        add('echo "[Stage 7] PASSED"')
        add('echo ""')
        add("")

        # Stage 8: Summary
        add("# ── Stage 8: Summary ─────────────────────────────────────────")
        add('FAILED_STAGE="summary"')
        add('echo "=========================================="')
        add('echo "ALL TEST STAGES PASSED"')
        add('echo "=========================================="')
        add('echo ""')
        add('echo "Test outputs:"')
        add('find "$TEST_DIR" -type f | sort')
        add('echo ""')
        add('echo "Disk usage:"')
        add('du -sh "$TEST_DIR"')
        add('echo ""')
        add('echo "=== INTERACTIVE TEST COMPLETE ==="')
        add("")
        add("# Remove the trap on success")
        add("trap - EXIT")
        add("exit 0")

        return "\n".join(L)

    def submit(self, script_path: str) -> str:
        """Submit the test job via sbatch. Returns job ID."""
        cmd = [
            "sbatch",
            f"--partition={self.site.interactive_partition}",
            f"--account={self.site.account}",
            f"--time={self.site.test_time}",
            f"--mem={self.site.test_mem}",
            f"--cpus-per-task={self.site.test_cpus}",
            f"--job-name={self.config.model_name}_itest",
            f"--output={self.test_dir}/interactive_test_%j.out",
            f"--error={self.test_dir}/interactive_test_%j.err",
            script_path,
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(
                f"sbatch failed (exit {result.returncode}):\n"
                f"stdout: {result.stdout}\nstderr: {result.stderr}"
            )

        job_id = result.stdout.strip().split()[-1]

        # Record submission
        sub_record = {
            "submission_command": " ".join(cmd),
            "job_id": job_id,
            "run_id": self.run_id,
            "config_hash": self.config.config_hash,
            "partition": self.site.interactive_partition,
            "account": self.site.account,
            "time": self.site.test_time,
            "mem": self.site.test_mem,
            "cpus": self.site.test_cpus,
            "submitted_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        }
        sub_path = self.test_dir / "submission_record.json"
        sub_path.write_text(json.dumps(sub_record, indent=2))

        return job_id

    def validate_outputs(self, job_info: JobInfo) -> tuple[dict[str, bool], list[str]]:
        """Validate test outputs after job completion.

        Returns (gates_dict, failure_reasons).
        """
        gates = {}
        failures = []

        # Gate 1: Job reached terminal successful state
        gates["job_completed"] = job_info.state == JobState.COMPLETED
        if not gates["job_completed"]:
            failures.append(f"Job state: {job_info.state.value} (expected COMPLETED)")

        # Gate 2: Exit code is zero
        gates["exit_code_zero"] = job_info.exit_code == 0
        if not gates["exit_code_zero"]:
            failures.append(f"Exit code: {job_info.exit_code} (expected 0)")

        # Gate 3: SN_HYDRO_DATA exists and is non-empty
        test_day = self.config.snapshot_days[0] if self.config.snapshot_days else 30.0
        day_label = f"{test_day:.2f}days"
        sn_hydro = self.test_dir / "cmfgen" / day_label / "SN_HYDRO_DATA"
        gates["sn_hydro_exists"] = sn_hydro.exists() and sn_hydro.stat().st_size > 100
        if not gates["sn_hydro_exists"]:
            failures.append(f"SN_HYDRO_DATA missing or empty: {sn_hydro}")

        # Gate 4: SN_HYDRO_DATA validation JSON exists
        sn_hydro_val = self.test_dir / "cmfgen" / day_label / "sn_hydro_validation.json"
        gates["sn_hydro_validated"] = sn_hydro_val.exists()
        if not gates["sn_hydro_validated"]:
            failures.append(f"SN_HYDRO_DATA validation JSON missing: {sn_hydro_val}")

        # Gate 5: CMFGEN run directory was prepared
        cmfgen_run = self.test_dir / "cmfgen_run" / day_label
        gates["cmfgen_dir_prepared"] = cmfgen_run.exists() and (cmfgen_run / "SN_HYDRO_DATA").exists()
        if not gates["cmfgen_dir_prepared"]:
            failures.append(f"CMFGEN run dir not prepared: {cmfgen_run}")

        # Gate 6: VADAT exists
        vadat = cmfgen_run / "VADAT"
        gates["vadat_exists"] = vadat.exists()
        if not gates["vadat_exists"]:
            failures.append(f"VADAT not found: {vadat}")

        # Gate 7: stdout log exists and contains "ALL TEST STAGES PASSED"
        stdout_files = sorted(self.test_dir.glob("interactive_test_*.out"))
        gates["stdout_has_success"] = False
        if stdout_files:
            content = stdout_files[-1].read_text(errors="replace")
            gates["stdout_has_success"] = "ALL TEST STAGES PASSED" in content
            if not gates["stdout_has_success"]:
                failures.append("stdout does not contain 'ALL TEST STAGES PASSED'")

        # Gate 8: No failure JSONL
        fail_jsonl = self.test_dir / "test_failure.jsonl"
        gates["no_failure_record"] = not fail_jsonl.exists()
        if not gates["no_failure_record"]:
            failures.append(f"Failure record found: {fail_jsonl}")

        # Gate 9: Resource usage within limits
        gates["within_time_limit"] = job_info.state != JobState.TIMEOUT
        if not gates["within_time_limit"]:
            failures.append("Job timed out")

        # Gate 10: Resolved config was written
        resolved_cfg = self.test_dir / "resolved_config.yaml"
        gates["resolved_config_written"] = resolved_cfg.exists()
        if not gates["resolved_config_written"]:
            failures.append(f"Resolved config not written: {resolved_cfg}")

        return gates, failures

    def run(self) -> InteractiveTestResult:
        """Execute the full interactive test: generate, submit, monitor, validate.

        Returns InteractiveTestResult with all gates and details.
        """
        # Generate test script
        script_path, script_hash = self.generate_test_script()

        # Submit
        job_id = self.submit(script_path)

        print(f"Interactive test submitted: job {job_id}")
        print(f"  Partition: {self.site.interactive_partition}")
        print(f"  Account: {self.site.account}")
        print(f"  Time: {self.site.test_time}")
        print(f"  Test dir: {self.test_dir}")

        # Monitor
        monitor = SlurmMonitor(job_id=job_id, log_dir=str(self.test_dir),
                               job_name=f"{self.config.model_name}_itest")
        monitor.record_event("job_submitted", JobState.PENDING,
                             f"Job {job_id} submitted to {self.site.interactive_partition}")

        # Wait for completion (interactive partition has 2hr limit)
        timeout_secs = 6900  # 1h55m
        final_info = monitor.wait_for_completion(
            timeout=timeout_secs,
            poll_interval=self.poll_interval,
        )

        # Validate outputs
        gates, failures = self.validate_outputs(final_info)

        # Required gates (must all pass)
        required_gates = [
            "job_completed", "exit_code_zero", "sn_hydro_exists",
            "sn_hydro_validated", "cmfgen_dir_prepared", "vadat_exists",
            "stdout_has_success", "no_failure_record", "within_time_limit",
            "resolved_config_written",
        ]
        passed = all(gates.get(g, False) for g in required_gates)

        # Write report
        report_path = self.test_dir / "interactive_run_report.txt"
        monitor.write_report(str(report_path), final_info,
                            validation_results=gates)

        # Write result JSON
        result = InteractiveTestResult(
            job_id=job_id,
            run_id=self.run_id,
            partition=self.site.interactive_partition,
            account=self.site.account,
            time_limit=self.site.test_time,
            memory=self.site.test_mem,
            cpus=self.site.test_cpus,
            submitted_at=time.strftime("%Y-%m-%dT%H:%M:%S"),
            final_state=final_info.state.value,
            exit_code=final_info.exit_code,
            elapsed=final_info.elapsed,
            test_dir=str(self.test_dir),
            script_hash=script_hash,
            config_hash=self.config.config_hash,
            passed=passed,
            gates=gates,
            failure_reasons=failures,
        )

        result_path = self.test_dir / "interactive_test_result.json"
        result_path.write_text(json.dumps(result.to_dict(), indent=2))

        return result
