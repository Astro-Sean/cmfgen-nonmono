# REPRODUCTION_RUNBOOK.md

## Prerequisites

- Access to viper cluster (MPCDF)
- Intel ifort 2024.0.2 compiler
- CMFGEN source at `/viper/ptmp1/sbrennan/SN2024pba/cmfgen_source`
- Clean baseline at `/viper/ptmp1/sbrennan/SN2024pba/cmfgen_baseline/clean_cmfgen/cur_cmf`
- X_v5 model config at `/viper/ptmp1/sbrennan/SN2024pba/X_v5/config_debug.yaml`

## Step 1: Verify Clean Baseline

```bash
cd /viper/ptmp1/sbrennan/SN2024pba/cmfgen_baseline
sha256sum cur_cmf_18jun25.tar.gz
# Expected: fc81d2dcf2df6538058ee46eab27ab231e01152e02c12666906d583c4dae4a00
```

## Step 2: Verify Working Tree

```bash
cd /viper/ptmp1/sbrennan/SN2024pba
diff -rq cmfgen_baseline/clean_cmfgen/cur_cmf/new_main cmfgen_source/new_main | grep differ
# Should show 9 modified files (after cleanup)
```

## Step 3: Build CMFGEN

```bash
cd /viper/ptmp1/sbrennan/SN2024pba/cmfgen_source/new_main
make clean && make
# Verify executable:
ls -la /viper/ptmp1/sbrennan/SN2024pba/cmfgen_source/exe/cmfgen_dev.exe
sha256sum /viper/ptmp1/sbrennan/SN2024pba/cmfgen_source/exe/cmfgen_dev.exe
```

## Step 4: Submit X_v5 Model (Stable Config)

Use the stable bootstrap configuration:
```yaml
cmfgen:
  inc_rad_decays: false
  vadat_overrides:
    FIX_NE: T
    FIX_T: T
    FIX_IMP: T
    MAX_LAM: 1.0D3
    MAX_LIN: 1.0D1
    LAM_VAL: 1.0D2
    NUM_LAM: 5
    DO_NG: F
    DO_LAM_IT: T
    DO_LAM_AUTO: F
    N_OB_INS: '0'
    N_IB_INS: '0'
    REV_RGRID: F
  vel_law: 12
```

Submit:
```bash
cd /viper/ptmp1/sbrennan/SN2024pba
python scripts/run_model_grid.py X_v5/config_debug.yaml
```

## Step 5: Monitor

```bash
# Job status
squeue -u sbrennan | grep X_v5

# Iteration count
cat X_v5/cmfgen/19.59days/SCRATCH/POINT1

# Moment errors
wc -l X_v5/cmfgen/19.59days/MOM_J_ERRORS

# Convergence
tail -5 X_v5/cmfgen/19.59days/CORRECTION_SUM

# Internal grid
head -5 X_v5/cmfgen/19.59days/fort.127
wc -l X_v5/cmfgen/19.59days/fort.127
```

## Step 6: Verify VEL_LAW=12 Active

```bash
grep "VEL_LAW" X_v5/cmfgen/19.59days/VADAT
# Should show: 12 [VEL_LAW]

grep "RD_IN_R_GRID" X_v5/cmfgen/19.59days/VADAT
# Should show: F [RD_IN_R_GRID]

grep "SN_HYDRO_MODEL" X_v5/cmfgen/19.59days/VADAT
# Should show: T [SN_HYDRO_MODEL]
```

## Step 7: Apply Patch to Clean Baseline (for verification)

```bash
cd /viper/ptmp1/sbrennan/SN2024pba
cp -a cmfgen_baseline/clean_cmfgen/cur_cmf cmfgen_baseline/patched_cmfgen
cd cmfgen_baseline/patched_cmfgen
patch -p1 < /viper/ptmp1/sbrennan/SN2024pba/docs/nonmonotonic_solver/patches/minimal_patch.diff
# Build and test
```

## Step 8: Revert Patch (for regression test)

```bash
cd cmfgen_baseline/patched_cmfgen
patch -R -p1 < /viper/ptmp1/sbrennan/SN2024pba/docs/nonmonotonic_solver/patches/minimal_patch.diff
# Should produce identical source to clean baseline
diff -rq . ../clean_cmfgen/cur_cmf
```
