# CLEANUP_PLAN.md

## Overview

The existing working tree contains 20 modified Fortran source files. 14 are required for VEL_LAW=12 non-monotonic support, 3 are related but need verification, 3 are unrelated debug diagnostics, and 1 is dangerous.

## Cleanup Actions

### 1. Revert DANGEROUS change (immediate)

**File:** `misc/rewrite_old_model_data.f`
**Action:** Revert `STOP` back to `RETURN`
**Reason:** The `STOP` was likely a debug artifact. It would terminate CMFGEN when this code path is reached, breaking model restart functionality.

### 2. Remove UNRELATED debug diagnostics

**Files:**
- `new_main/subs/solveba_v13.f` — remove JINC/JDEC diagnostic block
- `subs/ltepop_wld_v2.f` — remove Fe2 LTE diagnostic at depth 176
- `subs/regrid_log_dc_v1.f` — remove Fe2_IN regrid diagnostic at depth 176

**Reason:** These are hardcoded debug prints for specific depths/species that spam stdout and provide no general diagnostic value. They should not be in production code.

**Action:** Revert these three files to the clean baseline.

### 3. Verify RELATED changes

**Files:**
- `new_main/mod_subs/set_abund_clump.f` — verify T/ED interpolation is correct
- `new_main/subs/rd_sn_data.f` — verify population fraction check doesn't affect standard models
- `new_main/subs/generate_full_matrix_v3.f` — verify POPS zero guard doesn't change results

**Action:** Run the monotonic regression test with and without these changes. If results are identical within tolerance, keep them. If not, investigate.

### 4. Keep REQUIRED changes

All 14 REQUIRED files are kept as-is. They implement:
- VEL_LAW=12 grid construction (rv_sn_model_sniin.f)
- Dessart+2015 GAM guard formula (6 Sobolev routines + 2 fg_j_cmf routines)
- Sigma clamp for moment solver and derivatives (mom_j_cmf_v11.f, fg_j_cmf_v13.f)
- Ng extrapolation clamp (do_ng_accel_band_v4.f)
- SN_HYDRO_MODEL activation for VEL_LAW=12 (rd_control_variables.f)

### 5. Archive existing working tree

Before applying any cleanup, archive the current working tree:
```
cp -a cmfgen_source cmfgen_source_archive_20260821
```

### 6. Create clean implementation branch

Start from the clean baseline, apply only the REQUIRED changes, then verify:
1. Monotonic regression test passes (same as baseline within tolerance)
2. VEL_LAW=12 with monotonic velocity passes
3. VEL_LAW=12 with non-monotonic velocity runs without crash

## Rejected Changes

| File | Change | Reason for rejection |
|------|--------|---------------------|
| `misc/rewrite_old_model_data.f` | `RETURN` → `STOP` | Dangerous; breaks restart |
| `solveba_v13.f` | JINC/JDEC debug | Unrelated; log spam |
| `ltepop_wld_v2.f` | Fe2 depth-176 debug | Unrelated; hardcoded |
| `regrid_log_dc_v1.f` | Fe2_IN depth-176 debug | Unrelated; hardcoded |
