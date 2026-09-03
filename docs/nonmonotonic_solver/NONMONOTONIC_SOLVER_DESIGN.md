# NONMONOTONIC_SOLVER_DESIGN.md

## Approach

The implementation follows the Dessart+Hillier approach for applying CMFGEN to non-monotonic velocity fields:

1. **Retain hydrodynamic radial structure** — VEL_LAW=12 reads R, V from input_hydro.dat
2. **Do not force homologous expansion** — SIGMA computed from actual dV/dR, not assumed V∝R
3. **Use local velocity gradient correctly** — dV/dR computed via distance-weighted finite difference
4. **Support velocity reversals** — Dessart+2015 guard handles σ<-1 algebraically
5. **Preserve population/RT separation** — changes only affect GAM/SIGMA computation, not rate equations
6. **Document limitations** — Sobolev approximation is local; multiple resonance points not handled

## Mathematical Convention

### Sigma definition

σ = dlnV/dlnR - 1 = (R/V)·(dV/dR) - 1

### GAM (Sobolev optical depth scale)

For the Sobolev approximation with non-monotonic velocity:

**Baseline (monotonic):**
```
GAM = GAMH / (1 + (Z/R)² · σ)
```
or
```
GAM = GAMH / |1 + (Z/R)² · σ|
```

**Modified (Dessart+2015 guard, non-monotonic):**
```
T1 = 1 + (Z/R)² · σ
T1 = MAX(SQRT(T1² + 0.01), -0.2·(1+σ), 0.1)
GAM = GAMH / T1
```

This is smooth, positive-definite, and well-defined for all σ including σ<-1.

### Moment equation coefficients (PSI)

**Baseline:**
```
PSI = DTAUONQ · GAM · (1 + σ · K_ON_J)
```

**Modified (clamp):**
```
T1 = MAX(σ, -0.999)
PSI = DTAUONQ · GAM · (1 + T1 · K_ON_J)
```

The clamp ensures `(1 + T1·K_ON_J) > 0` when `K_ON_J ≤ 1` (which is guaranteed by the Eddington factor closure).

## Limitations

1. **Sobolev approximation**: The implementation uses the Sobolev (large velocity gradient) approximation. This is valid when the velocity scale height is much larger than the photon mean free path. For velocity plateaus (dV/dR ≈ 0) or discontinuities, the Sobolev approximation breaks down.

2. **Single resonance**: The Sobolev approximation assumes a single resonance point along each ray. For velocity reversals, multiple resonance points can exist. The current implementation does not handle multiple resonances — it uses the local velocity gradient at each depth point.

3. **Sigma clamp**: The -0.999 clamp is a numerical stabilization, not a physical correction. It creates a flat σ plateau in regions where σ would naturally be < -1. This affects the moment equation coefficients but not the GAM computation (which uses the Dessart guard).

4. **No full non-monotonic transfer**: This is NOT a full non-monotonic radiative transfer solver. It is a local Sobolev/single-resonance approximation with algebraic guards for σ<-1. A full solver would require ray-by-ray integration with multiple resonance detection.

5. **Grid reconstruction**: RV_SN_MODEL_SNIIN reconstructs the grid from input_hydro.dat using a log-tau spacing. This means the final CMFGEN grid is NOT the same as the converter grid. The reconstruction can introduce interpolation artifacts in regions with steep velocity gradients.

## Changed Routines

| Routine | Change | Reversible? |
|---------|--------|-------------|
| RV_SN_MODEL_SNIIN | V_HYDRO index fix, ordering fix, N_OB_INS param, sigma clamp | Yes |
| CMFGEN_SUB | Pass N_OB_INS to RV_SN_MODEL_SNIIN | Yes |
| RD_CONTROL_VARIABLES | Set SN_HYDRO_MODEL for VEL_LAW=12 | Yes |
| FG_J_CMF_V13 | Dessart guard for GAM, clamp for dGAMdR/GAMH | Yes |
| FG_J_CMF_V12 | Dessart guard for GAM | Yes |
| MOM_J_CMF_V11 | Clamp for PSI/PSIPREV/AV_SIGMA | Yes |
| SOBEW | Dessart guard | Yes |
| SOBEW_GRAD | Dessart guard | Yes |
| SOBEW_GRAD_V2 | Dessart guard | Yes |
| SOBEW_GRAD_V3 | Dessart guard | Yes |
| SOBJBAR_SIM | Dessart guard | Yes |
| SOBJBAR (disp) | Dessart guard | Yes |
| DO_NG_ACCEL_V4 | Post-Ng sigma clamp | Yes |
| SET_ABUND_CLUMP | Ordering fix, T/ED interpolation | Yes |

## Unchanged Routines

All other CMFGEN routines are unchanged from the clean baseline. Specifically:
- Rate equations (SOLVE_FOR_POPS, etc.)
- Atomic data reading
- Opacity computation
- Boundary conditions (except PSI clamp in mom_j_cmf_v11)
- Time-dependent physics
- Grid adjustment (ADJUST_R_GRID_V4 — already had sigma clamp in baseline)
