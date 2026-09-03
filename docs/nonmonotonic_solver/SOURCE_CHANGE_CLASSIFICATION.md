# SOURCE_CHANGE_CLASSIFICATION.md

## Summary

| Category | Files | Action |
|----------|-------|--------|
| REQUIRED (kept) | 13 | Retained — directly needed for VEL_LAW=12 |
| RELATED (kept) | 3 | Retained — needs verification but consistent |
| UNRELATED (reverted) | 3 | Reverted to clean baseline |
| DANGEROUS (reverted) | 1 | Reverted to clean baseline |
| Build config (kept) | 4 | Retained — necessary for viper cluster |

## Reverted Changes

### misc/rewrite_old_model_data.f — DANGEROUS
- **Change**: `RETURN` → `STOP`
- **Reason**: Would terminate CMFGEN when this code path is reached, breaking model restart
- **Action**: Reverted to clean baseline

### new_main/subs/solveba_v13.f — UNRELATED
- **Change**: Added JINC/JDEC diagnostic printing for first 5 iterations
- **Reason**: Debug diagnostic only, spams stdout
- **Action**: Reverted to clean baseline

### subs/ltepop_wld_v2.f — UNRELATED
- **Change**: Added LTE population diagnostic at depth 176 for Fe2
- **Reason**: Hardcoded debug for specific depth/species
- **Action**: Reverted to clean baseline

### subs/regrid_log_dc_v1.f — UNRELATED
- **Change**: Added regrid diagnostic for Fe2_IN at depth 176
- **Reason**: Hardcoded debug for specific depth/species
- **Action**: Reverted to clean baseline

## Retained Changes (REQUIRED)

All 13 REQUIRED files implement the non-monotonic velocity solver:

1. **rv_sn_model_sniin.f** — VEL_LAW=12 grid construction, V_HYDRO index fix, ordering fix, sigma clamp
2. **cmfgen_sub.f** — Pass N_OB_INS parameter
3. **rd_control_variables.f** — Set SN_HYDRO_MODEL for VEL_LAW=12
4. **fg_j_cmf_v13.f** — Dessart+2015 GAM guard, derivative clamp
5. **fg_j_cmf_v12.f** — Dessart+2015 GAM guard
6. **mom_j_cmf_v11.f** — Sigma clamp for moment equation coefficients
7-12. **Six Sobolev routines** — Dessart+2015 GAM guard (replaces ABS/unguarded division)
13. **do_ng_accel_band_v4.f** — Post-Ng sigma clamp

## Retained Changes (RELATED)

3 RELATED files are retained pending verification:

1. **set_abund_clump.f** — Ordering fix + T/ED interpolation from hydro grid
2. **rd_sn_data.f** — Population fraction check for species replacement
3. **generate_full_matrix_v3.f** — POPS zero guard before C_MAT multiplication

These are consistent with the non-monotonic solver changes but have not been independently verified to preserve standard CMFGEN behavior.

## Build Verification

After cleanup (reverting 4 files), the source compiles successfully:
```
ifort -qopenmp -qmkl ... -o cmfgen_dev.exe
Exit code: 0
```

Executable: `/viper/ptmp1/sbrennan/SN2024pba/cmfgen_source/exe/cmfgen_dev.exe`
SHA-256: `bff8c45725a0c07b182042b7b55b13442aa5336f881dafebda78f116d53257db`
