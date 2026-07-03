# CMFGEN Source Code Changes

This document lists all modifications made to the CMFGEN source code to fix runtime errors and enable successful execution with SN_HYDRO_DATA input.

## Changes Made

### 1. RD_SN_DATA - Allow Zero Isotopes

**File**: `/ptmp/sbrennan/SN2024pba/cmfgen_source/new_main/subs/rd_sn_data.f`

**Issue**: CMFGEN stopped with error "ND, NSP, NISO or model time undefined" when `NISO=0` (no isotopes present in SN_HYDRO_DATA).

**Fix**: Removed the `NISO .EQ. 0` check from the error condition to allow models with zero isotopes.

**Location**: Lines 124-126

**Before**:
```fortran
ELSE IF(INDEX(STRING,'Radius grid') .NE. 0)THEN
  IF(NSP .EQ. 0 .OR. NX .EQ. 0 .OR. OLD_SN_AGE_DAYS.EQ. 0.0_LDP .OR. NISO .EQ. 0)THEN
    WRITE(LUER,*)'Error in RD_SN_DATA'
    WRITE(LUER,*)'ND, NSP, NISO or model time undefined'
    STOP
```

**After**:
```fortran
ELSE IF(INDEX(STRING,'Radius grid') .NE. 0)THEN
  IF(NSP .EQ. 0 .OR. NX .EQ. 0 .OR. OLD_SN_AGE_DAYS.EQ. 0.0_LDP)THEN
    WRITE(LUER,*)'Error in RD_SN_DATA'
    WRITE(LUER,*)'ND, NSP or model time undefined'
    STOP
```

**Impact**: Allows CMFGEN to process SN_HYDRO_DATA files that have no isotope populations (e.g., when radioactive decay is disabled).

---

### 2. SET_EDGE_FREQ_V3 - Filter Negative Edge Frequencies

**File**: `/ptmp/sbrennan/SN2024pba/cmfgen_source/newsubs/set_edge_freq_v3.f`

**Issue**: Negative edge frequencies were being added to the frequency grid, causing "Error in SET_CONT_FREQ_V4 - MIN_FREQ too big" when `FREQ(1) < MIN_FREQ`.

**Fix**: Added check to skip negative edge frequencies by requiring `EDGE_FREQ(J) .GT. 0.0_LDP` before adding to the frequency array.

**Location**: Line 68

**Before**:
```fortran
IF( T1 .GT. 0)THEN
```

**After**:
```fortran
IF( T1 .GT. 0 .AND. EDGE_FREQ(J) .GT. 0.0_LDP)THEN
```

**Impact**: Prevents invalid (negative) frequencies from being added to the continuum frequency grid, avoiding the MIN_FREQ check failure.

---

## Configuration Changes

### 3. VADAT_INIT - Lower MIN_CF

**File**: `/ptmp/sbrennan/SN2024pba/cmfgen_base/VADAT_INIT`

**Issue**: `MIN_CF=0.003D0` was too high for the actual frequency range in SN_HYDRO_DATA, causing "Error in SET_CONT_FREQ_V4 - MIN_FREQ too big".

**Fix**: Reduced `MIN_CF` from `0.003D0` to `1.0D-4` to accommodate the actual frequency range.

**Location**: Line 53

**Before**:
```
0.003D0         [MIN_CF]       ! [OPT] was 0.00149897; raise to avoid sparse far-IR
```

**After**:
```
1.0D-4         [MIN_CF]       ! [OPT] was 0.003; lower to accommodate SN_HYDRO_DATA frequency range
```

---

### 4. prepare_VADAT.py - Lower MIN_CF in All Profiles

**File**: `/ptmp/sbrennan/SN2024pba/scripts/prepare_VADAT.py`

**Issue**: Profile dictionaries had `MIN_CF=0.01D0`, which was too high for SN_HYDRO_DATA frequency range.

**Fix**: Updated `MIN_CF` from `0.01D0` to `1.0D-4` in all 4 profile dictionaries (init, stable, final, single_final).

**Locations**: Lines 1305, 1422, 1573, 1703

**Before**:
```python
"MIN_CF": "0.01D0",
```

**After**:
```python
"MIN_CF": "1.0D-4",
```

**Impact**: Ensures all CMFGEN runs (init, stable, final, single-final) use a MIN_CF value compatible with SN_HYDRO_DATA frequency ranges.

---

### 5. CMFGEN_SUB - Reset LTE_EST after initial cold-start setup

**File**: `/ptmp/sbrennan/SN2024pba/cmfgen_source/new_main/cmfgen_sub.f`

**Issue**: When the `single_final` profile starts from LTE (`LTE_EST=T`, `LIN_INT=F`), CMFGEN uses `SET_NEW_MODEL_ESTIMATES` to set departure coefficients to LTE. On subsequent Slurm requeues/restarts, the `VADAT` file still contained `LTE_EST=T`, so the model could re-initialize populations to LTE instead of continuing from the existing `SCRTEMP` / `*_IN` files, producing large convergence jumps (e.g. iteration 57 jumping to ~1e57 % changes).

**Fix**: After the initial `SET_NEW_MODEL_ESTIMATES` call, if `LTE_EST` was used, reset the internal flag `SET_LTE_AS_INIT_ESTIMATES=.FALSE.` and update the `VADAT` file to `LTE_EST=F` and `LIN_INT=T`. This makes later restarts read from the saved model files rather than re-doing the LTE cold-start.

**Location**: Lines 1132-1142

**Before**:
```fortran
IF(NEWMOD)THEN
  NITSF=0
  IREC=0
  LAST_NG=-1000  			!Must be -1000
  NEXT_NG=1000				!Must be initialized to 1000
  CALL SET_NEW_MODEL_ESTIMATES(POPS,Z_POP,NU,NU_EVAL_CONT,FQW,
	1            LUER,LUIN,NC,ND,NP,NT,NCF,N_LINE_FREQ,MAX_SIM)
END IF
```

**After**:
```fortran
IF(NEWMOD)THEN
  NITSF=0
  IREC=0
  LAST_NG=-1000  			!Must be -1000
  NEXT_NG=1000				!Must be initialized to 1000
  CALL SET_NEW_MODEL_ESTIMATES(POPS,Z_POP,NU,NU_EVAL_CONT,FQW,
	1            LUER,LUIN,NC,ND,NP,NT,NCF,N_LINE_FREQ,MAX_SIM)
!
! LTE_EST is only needed for the initial cold-start model setup. Once
! the initial estimates have been installed, disable it in the internal
! state and in VADAT so that any subsequent restart (e.g. Slurm requeue)
! continues from the existing SCRTEMP / *_IN files rather than re-setting
! departure coefficients to LTE and causing a convergence jump.
!
  IF(SET_LTE_AS_INIT_ESTIMATES)THEN
    SET_LTE_AS_INIT_ESTIMATES=.FALSE.
    CALL UPDATE_KEYWORD(.FALSE.,'[LTE_EST]','VADAT',L_TRUE,L_TRUE,LUIN)
    CALL UPDATE_KEYWORD(.TRUE.,'[LIN_INT]','VADAT',L_TRUE,L_TRUE,LUIN)
  END IF
END IF
```

**Impact**: Prevents single-final LTE cold-start models from re-initializing to LTE on restart, allowing them to continue convergence smoothly from the previous state.

---

## Summary

These changes enable CMFGEN to successfully process SN_HYDRO_DATA files with:
- Zero isotope populations (when radioactive decay is disabled)
- Negative edge frequencies (filtered out before use)
- Lower minimum frequency threshold (MIN_CF=1.0D-4) to accommodate actual frequency ranges
- Stable restarts for single-final LTE cold-start models (LTE_EST reset to F after the initial setup)

The model X1-noscat successfully reached iteration 2 after applying these fixes.
