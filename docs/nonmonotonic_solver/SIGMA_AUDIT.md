# SIGMA_AUDIT.md

## Sigma Definition

**σ = dlnV/dlnR - 1 = (R/V)·(dV/dR) - 1**

- Dimensionless
- Frame: comoving (CMF)
- Sign convention: σ=0 for homologous expansion (V∝R), σ=-1 for constant velocity, σ<0 for decelerating, σ>0 for accelerating
- For non-monotonic velocity: σ can be any real value, including σ<-1 (velocity decreasing inward faster than 1/R)

## Complete Call-Site Table

| File | Routine | Line | Sigma usage | Guard applied? |
|------|---------|------|-------------|----------------|
| `rv_sn_model_sniin.f` | RV_SN_MODEL_SNIIN | 262-269 | Computed via finite difference | Clamp to >= -0.999 (L274) |
| `fg_j_cmf_v13.f` | FG_J_CMF_V13 | 956 | GAM = T1*(1+σ*μ²) | Dessart guard: MAX(SQRT(T2²+0.01),-0.2*(1+σ),0.1) |
| `fg_j_cmf_v13.f` | FG_J_CMF_V13 | 979 | dGAMdR = GAM*σ/R | Local clamp T2=MAX(σ,-0.999) |
| `fg_j_cmf_v13.f` | FG_J_CMF_V13 | 986 | GAMH = (σ(I)+σ(I+1))/2 | Local clamp T1,T2=MAX(σ,-0.999) |
| `fg_j_cmf_v12.f` | FG_J_CMF_V12 | 944 | GAM = T1*(1+σ*μ²) | Dessart guard (same as v13) |
| `mom_j_cmf_v11.f` | MOM_J_CMF_V11 | 386 | Direct use in moment equations | Clamp to >= -0.999 |
| `mom_j_cmf_v11.f` | MOM_J_CMF_V11 | 545 | AV_SIGMA = (σ(I)+σ(I+1))/2 | Local clamp T1 |
| `mom_j_cmf_v11.f` | MOM_J_CMF_V11 | 584 | PSI = DTAUONQ*GAM*(1+σ*K_ON_J) | Local clamp T1 |
| `mom_j_cmf_v11.f` | MOM_J_CMF_V11 | 668 | PSI(1) = GAM*(HBC+NBC*σ) | Local clamp T1 |
| `mom_j_cmf_v11.f` | MOM_J_CMF_V11 | 679 | PSI(1) = GAM*(1+σ*K_ON_J) | Local clamp T1 |
| `sobew.f` | SOBEW | 118 | GAM = GAMH/ABS(1+Z²/R²*σ) | Dessart guard (replaces ABS) |
| `sobew_grad.f` | SOBEW_GRAD | 113 | GAM = GAMH/(1+Z²/R²*σ) | Dessart guard |
| `sobew_grad_v2.f` | SOBEW_GRAD_V2 | 128 | GAM = GAMH/ABS(1+Z²/R²*σ) | Dessart guard |
| `sobew_grad_v3.f` | SOBEW_GRAD_V3 | 131 | GAM = GAMH/ABS(1+Z²/R²*σ) | Dessart guard |
| `sobjbar_sim.f` | SOBJBAR_SIM | 113 | GAM = GAMH/ABS(1+Z²/R²*σ) | Dessart guard |
| `sobjbar.f` (disp) | SOBJBAR | 95 | GAM = GAMH/(1+Z²/R²*σ) | Dessart guard |
| `do_ng_accel_band_v4.f` | DO_NG_ACCEL | 178 | Post-Ng extrapolation | Clamp to >= -0.999 |
| `adjust_r_grid_v4.f` | ADJUST_R_GRID | 206 | Grid revision spline | Clamp: if σ<-1 and V<0.2, σ=-0.999; if σ<-1, ERROR |

## Two Guard Mechanisms

### 1. Dessart+2015 GAM Guard (Sobolev + formal solution)

```fortran
T1 = 1.0 + SIGMA*(MU**2)  ! or 1.0 + Z²/R²*SIGMA for Sobolev
T1 = MAX(SQRT(T1*T1 + 0.01), -0.2*(1.0+SIGMA), 0.1)
GAM = GAMH / T1  ! or GAM = T1_velocity * T1 for fg_j
```

This replaces the baseline `ABS(1+σ*μ²)` or `1+σ*μ²` with a smooth, always-positive expression:
- `SQRT(T1²+0.01)`: regularizes the denominator near T1=0
- `-0.2*(1+σ)`: ensures GAM doesn't grow unboundedly for large negative σ
- `0.1`: absolute floor

**This guard is algebraically safe for ALL σ values, including σ<-1.** It does not modify the SIGMA array itself.

### 2. Sigma Clamp (moment solver + derivatives)

```fortran
T1 = SIGMA(I)
IF(T1 .LT. -0.999_LDP) T1 = -0.999_LDP
! Use T1 instead of SIGMA(I)
```

This clamps σ to >= -0.999 for:
- `mom_j_cmf_v11.f`: PSI/PSIPREV coefficients (tridiagonal solver stability)
- `fg_j_cmf_v13.f`: dGAMdR and GAMH (derivative stability)
- `do_ng_accel_band_v4.f`: post-Ng extrapolation
- `rv_sn_model_sniin.f`: initial grid construction

**The clamp uses local variables (T1, T2, T3, T4) so the original SIGMA array is NOT modified**, except in `rv_sn_model_sniin.f` (initial construction) and `do_ng_accel_band_v4.f` (Ng extrapolation) where the SIGMA array itself is clamped.

## Clamp Validation

### Why the clamp is required

The moment equations use `(1+σ*K_ON_J)` terms. When σ<-1 and K_ON_J>0, this term goes negative, producing:
- Negative PSI coefficients → negative J from the tridiagonal solver
- Oscillating solutions → divergence
- NaN propagation through the rate equations

The Dessart+2015 guard protects GAM (the Sobolev optical depth scale), but does NOT protect the moment equation coefficients (PSI/PSIPREV). The clamp is therefore required as a complementary protection.

### AGENTS.md warning (2026-08-10)

Removing the clamp was tested and caused **catastrophic divergence** (MAXCH ~1E+18 within ~100 iterations). The blowup occurred at depths with **positive** σ (~0.9-1.0), not just the σ<-1 region. This indicates that letting σ cross -1 smoothly destabilizes the Ng acceleration's polynomial extrapolation across stored iterations.

### Classification: **B. Numerically stabilizing but physically approximate**

The clamp is:
- **Required** for convergence with the current moment solver
- **Not physically correct** — it creates a flat σ plateau at -0.999
- **Not needed** for the Dessart GAM guard (which is algebraically safe)
- **Needed** for PSI/PSIPREV/dGAMdR (which use σ directly)

## Consistency Check

All call sites now use one of two consistent guards:
1. **Dessart guard** for GAM computation (6 Sobolev routines + 2 fg_j routines)
2. **-0.999 clamp** for moment coefficients and derivatives (mom_j_cmf_v11 + fg_j_cmf_v13 dGAMdR/GAMH + do_ng_accel + rv_sn_model_sniin)

No call site uses the baseline `ABS()` or unguarded division. No call site modifies the SIGMA array except rv_sn_model_sniin.f (initial construction) and do_ng_accel_band_v4.f (Ng extrapolation).

## Stale Sigma Risk

After velocity remapping (e.g., by adjust_r_grid_v4.f), SIGMA is recomputed from the new V(R) via the cubic spline derivative. The recomputed SIGMA is then clamped in adjust_r_grid_v4.f. No stale SIGMA arrays persist because:
1. RV_SN_MODEL_SNIIN computes SIGMA from the final R/V grid
2. ADJUST_R_GRID_V4 recomputes SIGMA from the revised grid
3. DO_NG_ACCEL recomputes SIGMA from the extrapolated populations (then clamps)

## Recommended Diagnostic

For every depth point, record:
- R, V, dV/dR, dlnV/dlnR, raw σ, clamped σ, clamp activation flag

This diagnostic should be written to a separate file (e.g., `SIGMA_DIAGNOSTIC`) behind a debug flag to avoid impacting production I/O.
