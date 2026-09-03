# X_v5_RUN_REPORT.md

## Run History

### Run 1: Job 10973784 (original, INC_RAD_DECAYS=T, FIX_NE=F)
- **Config**: MAX_LAM=1.2D0, MAX_LIN=10, INC_RAD_DECAYS=T, FIX_NE=F
- **Result**: MOM_J error at iteration 126, luminosity blow-up (-4.78e6)
- **Cause**: Radioactive decay heating with free electron density during cold start
- **Status**: CANCELLED

### Run 2: Job 10974035 (bootstrap, INC_RAD_DECAYS=F, FIX_NE=T, MAX_LAM=1.05)
- **Config**: MAX_LAM=1.05D0, MAX_LIN=5, INC_RAD_DECAYS=F, FIX_NE=T
- **Result**: 55 iterations, zero MOM_J errors, stable convergence
- **Status**: CANCELLED (user requested faster convergence)

### Run 3: Job 10974089 (MAX_LAM=1.2)
- **Config**: MAX_LAM=1.2D0, MAX_LIN=10, INC_RAD_DECAYS=F, FIX_NE=T
- **Result**: Cancelled immediately (user requested MAX_LAM=1.0D3)
- **Status**: CANCELLED

### Run 4: Job 10974111 (MAX_LAM=1.0D3)
- **Config**: MAX_LAM=1.0D3, MAX_LIN=10, INC_RAD_DECAYS=F, FIX_NE=T, LAM_VAL=100, NUM_LAM=5
- **Result**: 130+ iterations, zero MOM_J errors, stable convergence (~15s/iter)
- **Status**: CANCELLED (user changed config aggressively)

### Run 5: Job 10974565 (user aggressive config changes)
- **Config**: MAX_LAM=1.2D0, MAX_LIN=1.2D0, INC_RAD_DECAYS=T, FIX_NE removed,
  USE_LAM_ES=T, DIF=T, IB_METH=DIFFUSION, T_INIT_TAU=20, DO_AV=T, ALL_FREQ=T,
  AUTO_SMOOTH=T, N_OB_INS=5, N_FIX_BA=5, MIN_CF=1.0D-2, NUM_LAM=6,
  base_grid_target=100, cmfgen_tau_max=30, max_points=150, min_nd=80
- **Result**: Segfault at iteration 135 (exit code 46)
- **Cause**: Likely from aggressive VADAT changes, not from non-monotonic solver
- **Status**: FAILED

## Analysis of Run 5 Segfault

The segfault occurred at iteration 135, well after grid construction (which happens at startup). The non-monotonic solver code was working correctly in Runs 2-4 with the same source but different VADAT parameters.

Key differences in Run 5 that could cause segfault:
1. **USE_LAM_ES=T**: AGENTS.md states the current setup uses USE_LAM_ES=F. Changing to T activates a different formal solution path that may have array dimension issues with the smaller grid.
2. **DIF=T + IB_METH=DIFFUSION**: Enables diffusion mode which may require different array sizing.
3. **T_INIT_TAU=20.0D0**: Was 1.0D-10. This changes the initial optical depth scale dramatically.
4. **Smaller grid (ND=134 vs 219)**: With N_OB_INS=5, MND=124. Some arrays may be dimensioned for larger grids.
5. **INC_RAD_DECAYS=T with FIX_NE removed**: Re-introduces the original blowup condition.

## Recommendation

Revert to the stable Run 4 configuration (MAX_LAM=1.0D3, INC_RAD_DECAYS=F, FIX_NE=T) and apply changes one at a time to isolate the segfault cause.

## Grid Diagnostics (Run 4, stable)

| Metric | Value |
|--------|-------|
| ND | 219 |
| N_OB_INS | 0 (fallback to 3) |
| R range | 1.0177e+05 to 2.3167e+06 (10^10 cm) |
| V range | 113.02 to 5839.20 km/s |
| SIGMA range | -0.999 to 12.12 |
| N(SIGMA<0) | 77 |
| N(SIGMA==-0.999) | 74 |
| Cell ratio range | 1.000 to 7.097 |
| Negative-sigma regions | 4 |

## Grid Diagnostics (Run 5, segfaulted)

| Metric | Value |
|--------|-------|
| ND | 134 |
| N_OB_INS | 5 |
| NHYDRO | 124 |
| fort.127 lines | 509 |
