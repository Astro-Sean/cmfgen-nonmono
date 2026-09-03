# VEL_LAW_12_DATA_PATH.md

## Complete Data Path

```
input_hydro.dat
  → rv_sn_model_sniin.f (RV_SN_MODEL_SNIIN)
    → reads R_HYDRO, V_HYDRO, TAU_HYDRO (7 columns per point)
    → constructs log-spaced tau grid (MND = ND - 2*NBND_INS points)
    → LIN_INTERP: interpolates R onto tau grid → TA(1:MND)
    → inserts NBND_INS boundary points at each end → R(1:ND)
    → LIN_INTERP: interpolates V from hydro grid onto TA grid → V_TEMP(1:MND)
    → linear extrapolation of V at boundary points
    → computes SIGMA via distance-weighted finite difference
    → clamps SIGMA to >= -0.999
  → cmfgen_sub.f
    → R(1:ND), V(1:ND), SIGMA(1:ND) are the model arrays
    → SET_ABUND_CLUMP: reads input_hydro.dat again for rho/T, interpolates onto R grid
    → FG_J_CMF_V13: uses R_RAY, V_RAY, SIGMA_RAY for GAM computation
    → MOM_J_CMF_V11: uses SIGMA for moment equation coefficients (PSI/PSIPREV)
    → SOBEW/SOBEW_GRAD/SOBJBAR: uses SIGMA for Sobolev GAM
  → observer frame (disp/sobjbar.f): uses SIGMA for Sobolev GAM in spectrum
```

## Input File: `input_hydro.dat`

| Field | Value |
|-------|-------|
| Format | Free-format ASCII |
| Line 1 | `ADD_ACC_ZONE` (logical: T/F) |
| Lines 2-3 (if ADD_ACC_ZONE) | RINT, HRHO, VRAT, VMIN, BETA_ACC |
| Next line | Header string (ignored) |
| Next line | `NHYDRO` (integer: number of hydro points) |
| Next line | Bogus string (from Python writer) |
| Next NHYDRO lines | `R, V, RHO, T, XL, XTAU, TAU` (7 floats per line) |

### Column descriptions

| Column | Variable | Units | Notes |
|--------|----------|-------|-------|
| 1 | R | 10^10 cm | Radius, outer-to-inner (decreasing) |
| 2 | V | km/s | Velocity, sign preserved (no abs) |
| 3 | RHO | g/cm^3 | Density |
| 4 | T | K | Temperature (Kelvin, NOT 10^4 K) |
| 5 | XL | — | Unused (dummy) |
| 6 | XTAU | — | Unused (dummy) |
| 7 | TAU | — | Electron-scattering optical depth (outer=small, inner=large) |

### Ordering

Data is **outer-to-inner** (decreasing R, increasing TAU). This matches the CMFGEN convention where index 1 = outermost (RMAX) and index ND = innermost (RP).

### Critical requirements

1. **R must be strictly decreasing** — duplicate or non-decreasing R causes LIN_INTERP failure
2. **V must be nonzero** — SIGMA computation divides by V
3. **TAU must be positive** — RV_SN_MODEL_SNIIN takes LOG(TAU) for grid construction
4. **Outer TAU must be > 0** — LOG(TAU_HYDRO(1)) is used; zero causes NaN/SIGFPE
5. **TAU must increase inward** — geometric grid spans TAU(outer) to TAU(inner)
6. **NHYDRO must equal ND** — mismatch causes dimension errors (not explicitly checked in baseline)

## Internal CMFGEN Units

| Quantity | Input units | Internal units | Conversion |
|----------|-------------|----------------|------------|
| R | 10^10 cm | 10^10 cm | None (same) |
| V | km/s | km/s | None (same) |
| RHO | g/cm^3 | g/cm^3 | None (same) |
| T | K | 10^4 K | Divided by 1e4 in SET_ABUND_CLUMP |
| TAU | dimensionless | dimensionless | None (same) |
| SIGMA | dimensionless | dimensionless | Computed internally |

## SIGMA Computation (VEL_LAW=12)

From `rv_sn_model_sniin.f` lines 262-269:

```fortran
! Interior points (distance-weighted central difference)
T1=(V(I-1)-V(I))/(R(I-1)-R(I))
T2=(V(I)-V(I+1))/(R(I)-R(I+1))
SIGMA(I)=(T1*(R(I)-R(I+1))+T2*(R(I-1)-R(I)))/(R(I-1)-R(I))
SIGMA(I)=R(I)*SIGMA(I)/V(I)-1.0_LDP

! Outer boundary (one-sided)
SIGMA(1)=R(1)*(V(1)-V(2))/(R(1)-R(2))/V(1)-1.0_LDP

! Inner boundary (one-sided)
SIGMA(ND)=R(ND)*(V(ND-1)-V(ND))/(R(ND-1)-R(ND))/V(ND)-1.0_LDP
```

Mathematical definition: **σ = (R/V)·(dV/dR) - 1 = dlnV/dlnR - 1**

## Grid Construction

The VEL_LAW=12 grid is NOT the converter grid. RV_SN_MODEL_SNIIN constructs a new grid:

1. **Tau grid**: MND = ND - 2*NBND_INS points, log-spaced from TAU(outer) to TAU(inner)
2. **R grid**: R interpolated from hydro grid onto tau grid (LIN_INTERP in log-tau space)
3. **Boundary points**: NBND_INS points inserted at each boundary with fractional spacing
4. **V grid**: V interpolated from hydro grid onto R grid, with linear extrapolation at boundaries

The converter's SN_HYDRO_DATA grid is only used indirectly — `prepare_VADAT.py` writes `input_hydro.dat` from SN_HYDRO_DATA, and RV_SN_MODEL_SNIIN reconstructs the grid from that file.

## RD_IN_R_GRID Behavior

| RD_IN_R_GRID | Behavior |
|--------------|----------|
| F (current) | Reads input_hydro.dat, constructs tau-based grid, interpolates V from hydro |
| T | Reads RDINR file for R grid, computes V=VCORE*(R/R(ND))^BETA1, SIGMA=BETA1-1 (power law, ignores hydro) |

**RD_IN_R_GRID=F is correct for VEL_LAW=12.** Setting T overrides the hydro data with a power law.

## Startup Diagnostic (recommended addition)

Add to `rv_sn_model_sniin.f` after line 280:

```fortran
WRITE(LUER,'(A)') 'VEL_LAW=12 ACTIVE: RV_SN_MODEL_SNIIN'
WRITE(LUER,'(A,I6)') '  ND = ', ND
WRITE(LUER,'(A,I6)') '  NHYDRO = ', NHYDRO
WRITE(LUER,'(A,2ES14.4)') '  R range = ', R(1), R(ND)
WRITE(LUER,'(A,2ES14.4)') '  V range = ', V(1), V(ND)
WRITE(LUER,'(A,2ES14.4)') '  SIGMA range = ', MINVAL(SIGMA), MAXVAL(SIGMA)
WRITE(LUER,'(A,I6)') '  N(SIGMA<0) = ', COUNT(SIGMA .LT. 0.0_LDP)
WRITE(LUER,'(A,I6)') '  N(SIGMA<-1) = ', COUNT(SIGMA .LT. -1.0_LDP)
WRITE(LUER,'(A,I6)') '  N(SIGMA==-0.999) = ', COUNT(ABS(SIGMA+0.999_LDP) .LT. 1.0E-10_LDP)
WRITE(LUER,'(A,I6)') '  NBND_INS = ', NBND_INS
```
