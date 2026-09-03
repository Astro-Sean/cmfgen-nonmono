# REPRODUCTION_RUNBOOK.md (Proxy Populations)

**Date:** 2026-08-21

---

## 1. Prerequisites

- Access to the Viper SLURM cluster at MPCDF
- CMFGEN source compiled at `cmfgen_source/exe/cmfgen_dev.exe`
- Atomic data at `/u/sbrennan/atomic/` (symlinked from model directories)
- Python 3.12+ with numpy, pytest
- A CMFGEN model directory with MODEL, VADAT, MODEL_SPEC, SN_HYDRO_DATA,
  F_OSCDAT/F_TO_S/COL_DATA symlinks, T_IN, GAMMAS_IN

---

## 2. Generate Proxy Populations

```bash
cd /viper/ptmp1/sbrennan/SN2024pba

# Generate LTE proxy _IN files (conservative default)
python scripts/proxy_populations_v2.py SN_12C/cmfgen/0.5days \
  --output-dir runs/proxy_validation/SN_12C_0.5d \
  --strategy LTE

# Verify with independent parser
python tests/test_independent_in_parser.py \
  runs/proxy_validation/SN_12C_0.5d \
  SN_12C/cmfgen/0.5days/MODEL_SPEC \
  runs/proxy_validation/SN_12C_0.5d/independent_parse_report.json

# Build manifests
python tests/build_manifests.py \
  SN_12C/cmfgen/0.5days \
  runs/proxy_validation/SN_12C_0.5d
```

---

## 3. Run Tests

```bash
cd /viper/ptmp1/sbrennan/SN2024pba
python -m pytest scripts/test_proxy_populations_v2.py -q
# Expected: 57 passed
```

---

## 4. Download Reference Files

```bash
mkdir -p data/reference/M15-He-O-Ca-evo-new-4-all/33
cd data/reference/M15-He-O-Ca-evo-new-4-all/33
BASE="https://ttt.astro.su.se/~maer0651/M15-He-O-Ca-evo-new-4-all/33"
# Download key files (see REFERENCE_COMPARABILITY_REPORT.md for full list)
for f in MODEL_SPEC MODEL VADAT SN_HYDRO_DATA GAMMAS_IN IN_ITS \
         HIOUT HeIOUT He2OUT OIOUT CaIOUT Ca2OUT \
         CORRECTION_SUM CORRECTION_LINK OUTGEN; do
  curl -sSL -o "$f" "$BASE/$f"
done
sha256sum *  # Record checksums
```

---

## 5. Compare OUT Files (Structural)

```bash
cd /viper/ptmp1/sbrennan/SN2024pba
python scripts/compare_out_files.py
# Outputs: OUT_FILE_INVENTORY.csv, OUT_FILE_COMPARISON.csv, OUT_FILE_COMPARISON.md
```

---

## 6. CMFGEN Startup Test (Documented Diagnostic)

```bash
cd /viper/ptmp1/sbrennan/SN2024pba

# Create test directory
TESTDIR=runs/proxy_validation/cmfgen_startup_test
SRC=SN_12C/cmfgen/0.5days
mkdir -p "$TESTDIR/obs"

# Copy essential files
cp "$SRC"/{MODEL,VADAT,MODEL_SPEC,IN_ITS,SN_HYDRO_DATA,T_IN,GAMMAS_IN, \
     batch.sh,ADJUST_R_DEFAULTS,NUC_DECAY_DATA,FIN_CAL_GRID,RVTJ, \
     input_hydro.dat,arnaud_rothenflug.dat} "$TESTDIR/"
cp "$SRC/obs/"* "$TESTDIR/obs/"

# Copy atomic data symlinks
cp -a "$SRC"/*F_OSCDAT "$SRC"/*F_TO_S "$SRC"/*COL_DATA "$SRC"/*AUTO_DATA "$TESTDIR/"

# Copy proxy _IN files
cp runs/proxy_validation/SN_12C_0.5d/*_IN "$TESTDIR/"
cp runs/proxy_validation/SN_12C_0.5d/GAMMAS_IN "$TESTDIR/"

# Set IN_ITS to 2 iterations
cat > "$TESTDIR/IN_ITS" << 'EOF'
2                [NUM_ITS]
T [DO_LAM_IT]
F [DO_LAM_AUTO]
F [DO_GT_AUTO]
F [DO_T_AUTO]
T [D2_EQ_D1]
EOF

# Submit startup test
cd "$TESTDIR"
sbatch run_startup_test.slurm
# Monitor: squeue -u sbrennan
# Check: cat batch.log, head OUTGEN, head CORRECTION_SUM
```

---

## 7. Full Convergence Run (Production — via orchestrator)

**NOTE:** The orchestrator integration is currently incomplete (see
PROXY_POPULATIONS_AUDIT.md §6). Until fixed, production runs require manual
staging as in §6 above, but with IN_ITS set to a large number (e.g. 100000)
and submitted to the general partition.

To fix the orchestrator:
1. Add `--proxy-populations` argument to `scripts/prepare_VADAT.py`
2. Have `prepare_VADAT.py` call `proxy_populations_v2.py` to generate _IN files
3. Stage generated _IN files in the CMFGEN run directory
4. Set LTE_EST=F in VADAT when proxy populations are used
5. Submit via: `python scripts/run_model_grid.py configs/<model>.yaml --skip-stella --run-cmfgen`

---

## 8. Verify Deliverables

All deliverables should be present in `docs/` and `runs/proxy_validation/`:

| Deliverable | Location |
|---|---|
| PROXY_POPULATIONS_AUDIT.md | docs/ |
| ATOMIC_DATA_MANIFEST.json | runs/proxy_validation/SN_12C_0.5d/ |
| MODEL_ATOM_MANIFEST.json | runs/proxy_validation/SN_12C_0.5d/ |
| INPUT_GRID_MANIFEST.json | runs/proxy_validation/SN_12C_0.5d/ |
| PROXY_FORMAT_SPECIFICATION.md | docs/ |
| PROXY_ROUNDTRIP_REPORT.md | docs/ |
| REFERENCE_COMPARABILITY_REPORT.md | docs/ |
| OUT_FILE_INVENTORY.csv | runs/proxy_validation/ |
| OUT_FILE_COMPARISON.csv | runs/proxy_validation/ |
| OUT_FILE_COMPARISON.md | runs/proxy_validation/ |
| POPULATION_NORMALIZATION_REPORT.md | docs/ |
| FRESH_MODEL_CONVERGENCE_REPORT.md | docs/ |
| REPOSITORY_HYGIENE_REPORT_PROXY.md | docs/ |
| REPRODUCTION_RUNBOOK.md | docs/ (this file) |
