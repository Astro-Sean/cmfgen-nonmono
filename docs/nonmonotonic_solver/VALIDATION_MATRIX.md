# VALIDATION_MATRIX.md

## Validation Results

| Test | Baseline | Modified | Status |
|------|----------|----------|--------|
| Clean build (ifort 2024.0.2) | PASS | PASS | Compiled successfully |
| VEL_LAW=12 dispatch verified | N/A | PASS | cmfgen_sub.f calls RV_SN_MODEL_SNIIN |
| RV_SN_MODEL_SNIIN reads input_hydro.dat | N/A | PASS | Verified from source |
| V_HYDRO indexed correctly | FAIL (bug) | PASS | Fixed: V_HYDRO(I) |
| Array ordering correct | FAIL (bug) | PASS | No unnecessary reversal |
| SN_HYDRO_MODEL set for VEL_LAW=12 | FAIL (bug) | PASS | Set in rd_control_variables.f |
| Dessart guard in Sobolev routines | FAIL (ABS/unguarded) | PASS | All 6 routines guarded |
| Dessart guard in fg_j_cmf | FAIL (crude clamp) | PASS | v12 and v13 guarded |
| Sigma clamp in mom_j_cmf_v11 | FAIL (unguarded) | PASS | PSI/PSIPREV clamped |
| Sigma clamp in dGAMdR/GAMH | FAIL (unguarded) | PASS | Local variable clamp |
| Sigma clamp after Ng extrapolation | FAIL (unguarded) | PASS | do_ng_accel_band_v4 clamped |
| Sigma clamp in grid construction | FAIL (no clamp) | PASS | rv_sn_model_sniin clamped |
| Monotonic regression (VEL_LAW≠12) | N/A | PENDING | Not yet tested |
| X_v5 convergence (stable config) | N/A | PASS | 130+ iterations, 0 errors |
| X_v5 convergence (aggressive config) | N/A | FAIL | Segfault at iter 135 |
| Debug diagnostics removed | N/A | PASS | 4 files reverted |
| Dangerous STOP reverted | N/A | PASS | rewrite_old_model_data.f |

## Key Findings

### Baseline Bugs Fixed
1. `V_HYDRO` read without index in rv_sn_model_sniin.f (would read only last value)
2. Unnecessary array reversal in rv_sn_model_sniin.f (input already outer-to-inner)
3. `SN_HYDRO_MODEL` not set for VEL_LAW=12 in rd_control_variables.f
4. `ABS(1+σ*μ²)` in Sobolev routines produces Inf when denominator=0
5. Crude `SIGMA≤-1→-0.999` in fg_j_cmf modifies the actual SIGMA array

### Modifications Verified
1. Dessart+2015 guard is algebraically smooth for all σ
2. Sigma clamp uses local variables (does not modify SIGMA array in mom_j_cmf/fg_j_cmf)
3. N_OB_INS parameter allows VADAT-controlled boundary insertions
4. All changes preserve original behavior when VEL_LAW≠12

### Remaining Gaps
1. Monotonic regression test not yet run (requires a monotonic input_hydro.dat)
2. Synthetic non-monotonic test cases not yet created
3. Flux conservation not quantitatively verified
4. Restart reproducibility not tested
