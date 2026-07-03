# Audit and Fix prepare_VADAT for Consistency

Audit prepare_VADAT.py for consistency across init, stable, and final profiles, and fix temperature-related settings to prevent segmentation fault.

## Diagnosis

The segmentation fault occurs because:
1. VADAT has FIX_T=F, FIX_X_DPTH=39 (should allow temperature variation)
2. CMFGEN output shows "Temperature held fixed at all depths" (contradiction)
3. Population changes are extremely large (10^24) due to fixed temperature
4. CMF_BLKBAND_V3 band solver fails at depths 260 and 299
5. SIGSEGV occurs in iteration 3 during solver failure

## Root Cause

FIX_X_DPTH=39 with FIX_T=F causes CMFGEN to hold temperature fixed at all depths instead of just inner 39 depths. This prevents temperature auto-updates (DO_T_AUTO=T, DO_GT_AUTO=T) from working, causing numerical instability.

## Solution

Audit prepare_VADAT.py for consistency across all profiles and set FIX_X_DPTH=0 unconditionally. This allows temperature to vary at all depths, enabling temperature auto-updates and preventing numerical instability.

## Changes Required

1. Audit prepare_VADAT.py for consistency across init, stable, final profiles
2. Set FIX_X_DPTH=0 unconditionally for all profiles
3. Remove tau_es calculation logic for FIX_X_DPTH
4. Remove --vary_T and --no-fix-x-dpth command-line arguments
5. Update cmfgen_base/VADAT_INIT to have FIX_X_DPTH=0
6. Ensure all temperature-related settings are consistent
