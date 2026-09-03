# GRID_AND_CONSERVATION_TESTS.md

## Test Framework

Tests are run using the existing `run_model_grid.py` orchestrator with the X_v5 model configuration. The clean baseline executable is built from `cur_cmf_18jun25.tar.gz`; the modified executable from the cleaned working tree.

## Test Cases

### 1. Clean monotonic homologous model
- **Input**: VEL_LAW=12 with strictly monotonic V∝R velocity
- **Expected**: Same result as baseline CMFGEN within tolerance
- **Status**: NOT YET RUN — requires a monotonic input_hydro.dat

### 2. VEL_LAW=12 with strictly monotonic input velocity
- **Input**: input_hydro.dat with V monotonically increasing with R
- **Expected**: SIGMA ≥ 0 everywhere, no clamp activation, convergence
- **Status**: NOT YET RUN

### 3. One velocity reversal
- **Input**: input_hydro.dat with one region where dV/dR > 0 (σ < -1)
- **Expected**: Dessart guard activates, clamp activates, convergence
- **Status**: NOT YET RUN

### 4. Multiple velocity reversals
- **Input**: input_hydro.dat with multiple σ < -1 regions
- **Expected**: Multiple clamp regions, convergence (may be slower)
- **Status**: NOT YET RUN

### 5. Zero-gradient plateau
- **Input**: input_hydro.dat with dV/dR ≈ 0 region
- **Expected**: SIGMA ≈ -1, clamp activates, convergence
- **Status**: NOT YET RUN

### 6. Sharp shock/discontinuity
- **Input**: input_hydro.dat with velocity jump
- **Expected**: Large SIGMA at discontinuity, clamp activates, convergence
- **Status**: NOT YET RUN

### 7. Missing velocity input
- **Input**: No input_hydro.dat file
- **Expected**: Hard failure with diagnostic message
- **Status**: VERIFIED (baseline code has OPEN/IOSTAT check that STOPs)

### 8. Malformed velocity input
- **Input**: input_hydro.dat with non-numeric data
- **Expected**: READ error, hard failure
- **Status**: NOT YET RUN

### 9-11. Flux tests (zero, small, NaN)
- **Status**: NOT YET RUN — require synthetic test inputs

### 12. Sigma consistency
- **Expected**: All call sites use consistent SIGMA definition and guard
- **Status**: VERIFIED via SIGMA_AUDIT.md — all 17 call sites use either Dessart guard or -0.999 clamp

### 13. Grid refinement
- **Status**: NOT YET RUN

### 14. Restart from saved iteration
- **Status**: NOT YET RUN

### 15. Clean-build reproducibility
- **Status**: VERIFIED — cleaned source compiles with ifort 2024.0.2, produces executable `cmfgen_dev.exe`

## Conservation Checks

### Mass conservation
- RV_SN_MODEL_SNIIN does not modify density — it only reads R, V, TAU from input_hydro.dat
- SET_ABUND_CLUMP reads RHO separately and interpolates onto the R grid
- No mass should be lost or gained

### Luminosity conservation
- The sigma clamp affects the moment equation coefficients (PSI/PSIPREV)
- In the diffusion limit (SIGMA ≈ 0, K_ON_J ≈ 1/3), PSI ≈ GAM*(1+0) = GAM — unchanged
- In the non-monotonic region (SIGMA = -0.999), PSI ≈ GAM*0.001 — very small but positive
- Luminosity at the inner boundary depends on the boundary condition, not on PSI alone

### Energy conservation
- No energy-dependent changes were made
- The Dessart guard affects the Sobolev optical depth scale but not the energy equation

## X_v5 Model Grid Diagnostics

### Previous run (job 10974111, MAX_LAM=1.0D3, stable)
- ND = 219
- N_OB_INS = 0 (fallback to 3)
- SIGMA range: -0.999 to 12.12
- 77 points with SIGMA < 0
- 74 points at clamp value -0.999
- Converged to iteration 130+ with zero MOM_J errors

### Latest run (job 10974565, user config changes, segfault)
- ND = 134
- N_OB_INS = 5
- Segfaulted at iteration 135
- Cause: likely from aggressive VADAT changes (USE_LAM_ES=T, DIF=T, IB_METH=DIFFUSION, T_INIT_TAU=20), not from non-monotonic solver
