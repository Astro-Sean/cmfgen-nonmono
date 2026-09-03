# EXISTING_DIFF_CLASSIFICATION.md

## Classification Criteria

- **REQUIRED**: Directly needed for VEL_LAW=12 non-monotonic velocity support
- **RELATED**: Potentially related but needs verification
- **UNRELATED**: Not related to non-monotonic support
- **DANGEROUS**: Could break standard CMFGEN behavior

## File-by-File Classification

### 1. `new_main/subs/rv_sn_model_sniin.f` — REQUIRED

**Changes:**
1. Added `N_OB_INS_IN` parameter to subroutine signature (allows VADAT to control boundary insertions)
2. Changed `NBND_INS` from hardcoded 3 to `N_OB_INS_IN` with fallback to 3
3. Fixed `V_HYDRO(I)` read (was `V_HYDRO` without index — bug in baseline)
4. Removed array reversal of `TA` and `R_HYDRO`/`V_HYDRO` (input_hydro.dat is already outer-to-inner)
5. Changed velocity interpolation from `LIN_INTERP(V,TA,ND,...)` to `LIN_INTERP(TA,V_TEMP,MND,...)` with boundary extrapolation
6. Added NBND_INS>3 branch for boundary point insertion
7. Added sigma clamp: `IF(SIGMA(I) .LT. -0.999_LDP)SIGMA(I)=-0.999_LDP`

**Assessment:** All changes are required for VEL_LAW=12. The baseline had a bug (V_HYDRO without index), wrong array ordering (reversing already-correct data), and no sigma clamp. The N_OB_INS parameter addition is needed for grid control.

### 2. `new_main/cmfgen_sub.f` — REQUIRED

**Change:** Added `N_OB_INS` argument to `RV_SN_MODEL_SNIIN` call.

**Assessment:** Required — matches the subroutine signature change in rv_sn_model_sniin.f.

### 3. `new_main/mod_subs/rd_control_variables.f` — REQUIRED

**Change:** Added `SN_MODEL=.TRUE.` and `SN_HYDRO_MODEL=.TRUE.` for VELTYPE=12.

**Assessment:** Required — the baseline did not set SN_HYDRO_MODEL for VEL_LAW=12, which meant `set_abund_clump.f` would not read the hydro data path. This is a baseline bug.

### 4. `plane/fg_j_cmf_v13.f` — REQUIRED

**Changes:**
1. Removed `IF(SIGMA_RAY(I) .LE. -1.0_LDP)SIGMA_RAY(I)=-0.999_LDP` (baseline clamp)
2. Replaced `GAM=T1*(1+SIGMA*MU^2)` with Dessart+2015 guard: `MAX(SQRT(T2*T2+0.01),-0.2*(1+SIGMA),0.1)`
3. Added SIGMA clamp to -0.999 for dGAMdR and GAMH computations (derivative guard)

**Assessment:** Required. The Dessart+2015 formula is algebraically safe for any SIGMA. The derivative guard prevents pathological dGAMdR when SIGMA<-1. The baseline had a crude `SIGMA<=-1 → -0.999` clamp that was too aggressive (it modified the actual SIGMA array, not just a local copy).

### 5. `plane/fg_j_cmf_v12.f` — REQUIRED

**Changes:** Same Dessart+2015 guard as v13. Removed baseline `SIGMA<=-1 → -0.999` clamp.

**Assessment:** Required — same rationale as v13.

### 6. `plane/mom_j_cmf_v11.f` — REQUIRED

**Changes:**
1. Added `IF(SIGMA(I) .LT. -0.999_LDP)SIGMA(I)=-0.999_LDP` at line 386
2. Added clamp for `AV_SIGMA` (cell-averaged sigma)
3. Added clamp for `PSI`/`PSIPREV` (moment equation coefficients)
4. Added clamp for boundary-condition PSI terms

**Assessment:** Required. The moment equations use `(1+SIGMA*K_ON_J)` terms that go negative when SIGMA<-1, breaking the tridiagonal solver. The clamp uses local variables (T1) so the original SIGMA array is not modified.

### 7-12. Sobolev routines (`sobew.f`, `sobew_grad.f`, `sobew_grad_v2.f`, `sobew_grad_v3.f`, `sobjbar_sim.f`, `disp/subs/sobjbar.f`) — REQUIRED

**Changes:** All six routines replaced `GAMH/ABS(1+Z^2/R^2*SIGMA)` or `GAMH/(1+Z^2/R^2*SIGMA)` with the Dessart+2015 guard: `GAMH/MAX(SQRT(T1*T1+0.01),-0.2*(1+SIGMA),0.1)`.

**Assessment:** Required. The baseline used `ABS()` or direct division, which produces Inf when `1+Z^2/R^2*SIGMA=0` (i.e., SIGMA=-R^2/Z^2). The Dessart guard is smooth and well-defined for all SIGMA values.

### 13. `new_main/subs/do_ng_accel_band_v4.f` — REQUIRED

**Change:** Added sigma clamp after Ng extrapolation: `IF(SIGMA(I) .LT. -0.999_LDP)SIGMA(I)=-0.999_LDP`

**Assessment:** Required. Ng acceleration extrapolates SIGMA across iterations; the extrapolated value can go below -1. The clamp prevents this from corrupting the moment solver.

### 14. `new_main/mod_subs/set_abund_clump.f` — RELATED

**Changes:**
1. Changed `R_HYDRO(NHYDRO-I+1)` read to `R_HYDRO(I)` (removed reversal)
2. Added T interpolation from hydro grid and ED estimation

**Assessment:** The reversal fix is RELATED (consistent with rv_sn_model_sniin.f ordering fix). The T/ED interpolation is needed because VEL_LAW=12 now sets SN_HYDRO_MODEL=.TRUE., which triggers the hydro data path in set_abund_clump. **Needs verification** that the T/ED initialization is correct.

### 15. `new_main/subs/rd_sn_data.f` — RELATED

**Change:** Added population fraction check to the species replacement condition.

**Assessment:** RELATED — prevents replacing species when both the old and new populations are negligibly small (< 1e-15). This avoids numerical noise in the species equation selection. **Needs verification** that it doesn't affect standard SN models.

### 16. `new_main/subs/generate_full_matrix_v3.f` — RELATED

**Change:** Added `IF(ABS(POPS(J,L)) .GT. 0.0_LDP)` guard before multiplying C_MAT by POPS.

**Assessment:** RELATED — prevents NaN propagation when a population is exactly zero. **Needs verification** that this doesn't change results for standard models.

### 17. `new_main/subs/solveba_v13.f` — UNRELATED (debug)

**Change:** Added JINC/JDEC diagnostic printing for first 5 iterations.

**Assessment:** UNRELATED — debug diagnostic only. Should be removed or placed behind a debug flag.

### 18. `subs/ltepop_wld_v2.f` — UNRELATED (debug)

**Change:** Added LTE population diagnostic at depth 176 for Fe2.

**Assessment:** UNRELATED — hardcoded debug diagnostic for a specific depth and species. Should be removed.

### 19. `subs/regrid_log_dc_v1.f` — UNRELATED (debug)

**Change:** Added regrid diagnostic for Fe2_IN at depth 176.

**Assessment:** UNRELATED — hardcoded debug diagnostic. Should be removed.

### 20. `misc/rewrite_old_model_data.f` — DANGEROUS

**Change:** Changed `RETURN` to `STOP`.

**Assessment:** DANGEROUS — this changes a normal return to a hard stop, which would terminate CMFGEN when this code path is reached. This appears to be a debugging artifact that was accidentally left in. **Must be reverted.**

### Build files — UNRELATED

`Makefile_definitions`, `new_main/Makefile`, `disp/Makefile`, `obs/Makefile` — build configuration for the viper cluster. Not related to non-monotonic solver.

## Summary

| Category | Count | Files |
|----------|-------|-------|
| REQUIRED | 14 | rv_sn_model_sniin.f, cmfgen_sub.f, rd_control_variables.f, fg_j_cmf_v13.f, fg_j_cmf_v12.f, mom_j_cmf_v11.f, sobew.f, sobew_grad.f, sobew_grad_v2.f, sobew_grad_v3.f, sobjbar_sim.f, sobjbar.f, do_ng_accel_band_v4.f |
| RELATED | 3 | set_abund_clump.f, rd_sn_data.f, generate_full_matrix_v3.f |
| UNRELATED | 3 | solveba_v13.f, ltepop_wld_v2.f, regrid_log_dc_v1.f |
| DANGEROUS | 1 | rewrite_old_model_data.f |
| UNRELATED (build) | 4 | Makefiles |
