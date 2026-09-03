# ZERO_FLUX_AUDIT.md

## Problem Statement

The "zero flux" issue refers to situations where the radiative flux is zero or near-zero at certain depths or frequencies, causing numerical problems in the moment solver or formal solution.

## Audit of Flux-Related Code Paths

### 1. Flux initialization

In `mom_j_cmf_v11.f`, the flux-like quantities PSI and PSIPREV are computed from GAM, DTAUONQ, and SIGMA. When SIGMA is clamped to -0.999, the term `(1 + SIGMA*K_ON_J)` approaches zero when K_ON_J approaches 1 (which is the Eddington factor at optically thick depths). This produces near-zero PSI, which is physically correct (zero flux in the diffusion limit) but can cause numerical issues if division by PSI occurs elsewhere.

### 2. Flux clipping

No explicit flux clipping was found in the modified code. The baseline does not clip flux to zero.

### 3. Flux as boundary condition

In `mom_j_cmf_v11.f` lines 668-671, the inner boundary condition uses:
```fortran
PSI(1) = GAM(1)*(HBC + NBC*SIGMA(1))
```
When SIGMA is clamped to -0.999, this becomes `GAM(1)*(HBC - 0.999*NBC)`. If HBC ≈ NBC (as in the Eddington approximation), PSI(1) ≈ 0, which is the correct diffusion-limit behavior.

### 4. Luminosity calculation

Luminosity is computed from the flux at the inner boundary. If PSI(1) ≈ 0, the luminosity can be near-zero or negative, which was the original "luminosity blow-up" symptom. However, the blow-up was traced to `INC_RAD_DECAYS=T` with `FIX_NE=F` during cold start, not to a zero-flux issue.

### 5. Comoving-to-observer-frame transfer

The `disp/sobjbar.f` routine uses the Dessart guard for GAM computation. No zero-flux issue was identified in the observer-frame path.

## Classification

The zero-flux issue in this model is:

**C. A convergence problem** — not an input, boundary-condition, underflow, or indexing bug.

The near-zero PSI at depths where SIGMA ≈ -0.999 and K_ON_J ≈ 1 is physically correct (diffusion limit). The convergence problem arises when the moment solver tries to iterate on populations near these depths — the near-zero PSI makes the tridiagonal system ill-conditioned, causing large population corrections that destabilize the iteration.

## Fix

No explicit zero-flux fix is needed. The sigma clamp at -0.999 (rather than -1.0) ensures that `(1 + SIGMA*K_ON_J)` remains slightly positive even when K_ON_J = 1:
- With SIGMA = -0.999 and K_ON_J = 1: PSI factor = 1 + (-0.999)(1) = 0.001 > 0
- With SIGMA = -1.0 and K_ON_J = 1: PSI factor = 1 + (-1.0)(1) = 0.0 (singular)

The 0.001 margin is small but sufficient for the tridiagonal solver. If ill-conditioning persists, the clamp could be tightened to -0.99 (giving a 0.01 margin), but this would further distort the physics.

## Regression Tests

The following tests should verify zero-flux behavior:
1. **Exact zero flux**: Set SIGMA=-1.0 at one depth → should be clamped to -0.999, PSI ≈ 0.001*GAM
2. **Very small flux**: Set SIGMA=-0.999 at one depth → PSI = 0.001*GAM*DTAUONQ
3. **NaN/Inf flux**: Set SIGMA=-1e10 → should be clamped to -0.999, no NaN
4. **Physical zero flux**: Diffusion limit at inner boundary → PSI ≈ 0, luminosity from boundary condition
