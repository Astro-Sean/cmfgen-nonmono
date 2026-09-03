# PATCH_SUMMARY.md

## Patch Statistics

| Metric | Value |
|--------|-------|
| Modified files (REQUIRED) | 13 |
| Modified files (RELATED) | 3 |
| Modified files (reverted) | 4 |
| Total changed lines (full patch) | 500 |
| Total changed lines (minimal patch) | 428 |

## Minimal Patch (13 files, REQUIRED only)

| File | Lines changed | Change description |
|------|--------------|-------------------|
| `rv_sn_model_sniin.f` | ~80 | V_HYDRO index fix, ordering fix, N_OB_INS param, boundary insertion for NBND_INS>3, sigma clamp |
| `cmfgen_sub.f` | 3 | Pass N_OB_INS to RV_SN_MODEL_SNIIN |
| `rd_control_variables.f` | 3 | Set SN_HYDRO_MODEL=.TRUE. for VEL_LAW=12 |
| `fg_j_cmf_v13.f` | ~40 | Dessart guard for GAM, clamp for dGAMdR/GAMH |
| `fg_j_cmf_v12.f` | ~8 | Dessart guard for GAM |
| `mom_j_cmf_v11.f` | ~40 | Clamp for PSI/PSIPREV/AV_SIGMA |
| `sobew.f` | 4 | Dessart guard |
| `sobew_grad.f` | 4 | Dessart guard |
| `sobew_grad_v2.f` | 4 | Dessart guard |
| `sobew_grad_v3.f` | 4 | Dessart guard |
| `sobjbar_sim.f` | 4 | Dessart guard |
| `sobjbar.f` (disp) | 4 | Dessart guard |
| `do_ng_accel_band_v4.f` | 7 | Post-Ng sigma clamp |

## Full Patch (16 files, REQUIRED + RELATED)

Adds 3 RELATED files:
| File | Lines changed | Change description |
|------|--------------|-------------------|
| `set_abund_clump.f` | ~19 | Ordering fix, T/ED interpolation from hydro grid |
| `rd_sn_data.f` | ~11 | Population fraction check for species replacement |
| `generate_full_matrix_v3.f` | ~10 | POPS zero guard before C_MAT multiplication |

## Reverted Changes (4 files)

| File | Reason |
|------|--------|
| `misc/rewrite_old_model_data.f` | DANGEROUS: RETURN→STOP would break restart |
| `new_main/subs/solveba_v13.f` | UNRELATED: debug diagnostic printing |
| `subs/ltepop_wld_v2.f` | UNRELATED: hardcoded Fe2 depth-176 debug |
| `subs/regrid_log_dc_v1.f` | UNRELATED: hardcoded Fe2_IN depth-176 debug |

## File-by-File Justification

### rv_sn_model_sniin.f
**Why:** This is the VEL_LAW=12 entry point. The baseline had three bugs: (1) `V_HYDRO` read without index, (2) unnecessary array reversal of already-correct outer-to-inner data, (3) no sigma clamp for non-monotonic velocities. The N_OB_INS parameter addition allows VADAT to control boundary insertions instead of hardcoding 3.

### cmfgen_sub.f
**Why:** Must pass N_OB_INS to match the updated RV_SN_MODEL_SNIIN signature.

### rd_control_variables.f
**Why:** The baseline did not set SN_HYDRO_MODEL=.TRUE. for VEL_LAW=12, which meant set_abund_clump.f would not read the hydro data path. This is a baseline bug.

### fg_j_cmf_v13.f / fg_j_cmf_v12.f
**Why:** The baseline used `GAM=T1*(1+SIGMA*μ²)` with a crude `SIGMA≤-1→-0.999` clamp that modified the actual SIGMA array. The Dessart+2015 guard is smooth, always-positive, and does not modify SIGMA. The derivative clamp (dGAMdR, GAMH) uses local variables.

### mom_j_cmf_v11.f
**Why:** The moment equations use `(1+SIGMA*K_ON_J)` terms that go negative when SIGMA<-1, breaking the tridiagonal solver. The clamp uses local variables (T1) so the original SIGMA array is not modified.

### sobew.f / sobew_grad*.f / sobjbar*.f
**Why:** The baseline used `ABS(1+Z²/R²*SIGMA)` or direct division, which produces Inf when the denominator is zero. The Dessart guard is smooth and well-defined for all SIGMA.

### do_ng_accel_band_v4.f
**Why:** Ng acceleration extrapolates SIGMA across iterations; the extrapolated value can go below -1. The clamp prevents this from corrupting the moment solver.

### set_abund_clump.f (RELATED)
**Why:** The ordering fix is consistent with rv_sn_model_sniin.f. The T/ED interpolation is needed because SN_HYDRO_MODEL=.TRUE. now triggers the hydro data path.

### rd_sn_data.f (RELATED)
**Why:** Prevents replacing species when both populations are negligibly small, avoiding numerical noise.

### generate_full_matrix_v3.f (RELATED)
**Why:** Prevents NaN propagation when a population is exactly zero.

## Backward Compatibility

- All changes preserve original behavior when VEL_LAW ≠ 12
- The Dessart guard reduces to the baseline expression when SIGMA > -1 (verified algebraically)
- The sigma clamp only activates when SIGMA < -0.999 (never in monotonic models)
- The N_OB_INS parameter defaults to 3 when VADAT specifies 0 (fallback)
- No input file format changes
- No output file format changes
