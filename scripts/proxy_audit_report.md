# Proxy Populations Audit Report

## Date: 2026-08-20
## File: scripts/proxy_populations.py (~2922 lines)

## Executive Summary

Total issues found: 67
- BLOCKING: 5
- HIGH: 18
- MEDIUM: 28
- LOW: 13 (not detailed here)

## Verified CMFGEN Conventions

### _IN File Format (verified from regrid_log_dc_v1.f)

- **Header**: `R(ND)  RLUM  NOLD  ND`
- **Per-depth row**: `R  DI  ED  T  0.0  V  CLUMP_FAC  I`
- **DC block**: NOLD values, 5 per line, E17.7 format
- **Depth ordering**: inner-to-outer (I=1 = innermost, R decreases with index)
- **R units**: 10^10 cm
- **T units**: 10^4 K
- **DI units**: cm^-3 (next ion stage density, NOT current ion)
- **ED units**: cm^-3
- **DC format**: LINEAR DC (CMFGEN takes LOG internally via TAKE_LOGS=.TRUE.)
  - Date string "07-Jul-1997" or any non-specific date → TAKE_LOGS=.TRUE.
  - CMFGEN does `DPOP = LOG(DPOP)` to convert linear DC to ln(DC)
  - Then `cnvt_fr_dc_v2.f`: `population = EXP(ln(DC) + ln(LTE_pop)) = DC * LTE_pop`
  - **Proxy writes linear DC → CORRECT**

### LTE Population Formula (verified from ltepop_wld_v2.f)

```
n_LTE(J) = W(J) * g(J) * exp(-EDGE(J)*HDKT/T) * (2.07078e-22 * ED * DI * T^-1.5) / GION
```

- **T**: local gas temperature (10^4 K), passed as argument, NOT from _IN file
- **DI**: next ion stage density (cm^-3), from _IN file
- **ED**: electron density (cm^-3), from GAMMAS_IN or _IN file
- **W**: occupation probability (level dissolution)
- **GION**: statistical weight of ion ground state
- **Constant**: 2.07078e-22 (cgs, verified)

### DC Application (verified from cnvt_fr_dc_v2.f)

```
population(J) = EXP(XzV_F(J) + LOG_XzVLTE(J))
             = EXP(ln(DC) + ln(LTE_pop))
             = DC * LTE_pop
```

- **XzV_F** = ln(DC) after regrid processing
- **LOG_XzVLTE** = ln(LTE population) from LTEPOP_WLD_V2

### SET_DC_LTE_V2 (verified from set_dc_lte_v2.f)

- T > TMIN_LOC: ln(DC) = 0 (DC = 1, full LTE)
- T <= TMIN_LOC: ln(DC) = HDKT * EDGE(J) * (1/TMIN_LOC - 1/T)
- TMIN_LOC = TMIN * (1 + max(0, EDGE(1)-3.5)/20)
- This boosts excited-level populations at low T to prevent zeros

### F_TO_S Mapping (verified from regrid_log_dc_v1.f lines 609-619)

When SUM(INT_SEQ) == 0:
- All levels in the same superlevel get the SAME DC
- DC taken from the LOWEST state in the group

### Charge Conservation (verified from steq_multi_v10.f)

```
STEQ_ED = sum over species of [(ZION-1)*N_levels + ZION*DI] - ED
```

## BLOCKING Issues

### B1: Two Saha implementations with different behavior
**Location**: compute_saha_ion_frac() L1097-1236 vs inline in iterate_grey_lambda() L1482-1527
**Issue**: The inline version renormalizes fractions to sum=1 (L1525-1526), the function version does not. This causes charge conservation drift between the grey-lambda iteration and the final DI computation.
**Fix**: Unify into one Saha solver, always renormalize.

### B2: No charge conservation check in grey-lambda convergence
**Location**: iterate_grey_lambda() L1575-1578
**Issue**: Convergence checks only relative change in ne and kappa. No absolute charge conservation error check. ED can converge to a value inconsistent with the ionization state.
**Fix**: Add charge conservation residual to convergence criteria.

### B3: Ionization energy fallback can silently use wrong atomic data
**Location**: IONIZATION_ENERGIES_CM_FALLBACK L981-995, used at L1152-1159
**Issue**: When F_OSCDAT is not found, falls back to hard-coded ionization energies without failing. These may not match the actual atomic data used by CMFGEN.
**Fix**: Fail loudly when F_OSCDAT is missing, or require explicit --allow-fallback flag.

### B4: VADAT abundance fallback
**Location**: L1720-1722
**Issue**: Falls back to depth-varying hydro abundances if VADAT not found. CMFGEN uses depth-averaged VADAT abundances in SOLVEBA. This can cause ~1e28 mismatch at composition boundaries (H/He shell).
**Fix**: Require VADAT, fail if not found. Report all abundance differences.

### B5: T_write = T_ion violates "preserve physical temperature" rule
**Location**: L1947, write_in_file() L1322
**Issue**: Default mode writes T_ion (radiation temperature) to the T column of _IN files. The audit requirement says "Do not write a synthetic T_ion into an _IN temperature column unless the exact CMFGEN reader proves that this column cannot alter the physical model."
**Verification**: LTEPOP_WLD_V2 uses T from the _IN file (via REGRID_T_ED_V3 which reads T_IN). However, the T column in species _IN files (e.g. HI_IN) is NOT used for the physical temperature structure — that comes from T_IN (read separately). The species _IN T column IS used by LTEPOP_WLD_V2 for the LTE population computation.
**Status**: NEEDS VERIFICATION — does REGRID_LOG_DC_V1 use the T from the _IN file or from T_IN? From the code, REGRID_LOG_DC_V1 reads OLD_T from the _IN file and interpolates it, but the calling code in set_new_model_estimates.f passes the model T (from T_IN) to LTEPOP_WLD_V2, not the species _IN T. The species _IN T is only used for interpolation reference (OLD_T).
**Fix**: Default to T_local (Mode A conservative). Make T_ion opt-in (Mode B).

## HIGH Priority Issues

### H1: Four DC computation methods with inconsistent physics
**Location**: compute_dc() L760-792 (orphaned), compute_dc_dilute_lte() L795-851, compute_dc_lte_cmfgen() L716-757, compute_dc_nlte_reduced() L854-969
**Issue**: Different exponent signs, different floors/caps (1e10 vs 1e30), different physical models.
**Fix**: Remove orphaned compute_dc(). Unify caps. Document each method clearly.

### H2: Statistical weight fallback (gion_default=1.0)
**Location**: get_atomic_info() L198-215, build_saha_atomic_data() L1077-1086
**Issue**: Defaults g_ion=1.0 for unknown species. Affects Saha ratio by factor of g_ion.
**Fix**: Require g_ion from F_OSCDAT. Fail if not available.

### H3: DC caps inconsistent (1e10 vs 1e30)
**Location**: compute_dc() L789 (1e10), compute_dc_dilute_lte() L848 (1e30), compute_dc_nlte_reduced() L961 (1e30)
**Fix**: Unify to a single configurable cap.

### H4: Smoothing affects DI computation
**Location**: smooth_depth_profile() L1239-1302, called at L1966
**Issue**: 3-point median filter on ion_frac. Can erase real ionization fronts.
**Fix**: Make opt-in, not default. Report all smoothed depths.

### H5: SAHA_CONST and CMFGEN_C should be verified
**Location**: SAHA_CONST=2.414e15 (L1004), CMFGEN_C=2.07078e-22 (L912)
**Issue**: Should match standard values and CMFGEN Fortran exactly.
**Verification**: CMFGEN_C=2.07078e-22 matches ltepop_wld_v2.f line 67. SAHA_CONST should be 2*(2*pi*m_e*k/h^2)^1.5 = 2.414e15 in cgs — verified correct.

### H6: Hard-coded Kramers coefficient
**Location**: L1378: `3.7e22`
**Issue**: No comment explaining origin. Should be documented or replaced with verified value.

### H7: Multiple g_ion fallbacks in build_saha_atomic_data
**Location**: L1077-1086
**Issue**: Multiple fallback paths for statistical weights.

### H8: Renormalization inconsistency between Saha implementations
**Location**: L1525-1526 (inline) vs compute_saha_ion_frac() (no renorm)
**Issue**: Can cause charge conservation drift.

### H9-H18: Various flooring and fallback issues
See detailed list in code audit. All should be made opt-in or fail loudly.

## MEDIUM Priority Issues (summary)

- SPEED_LIGHT precision mismatch (L19 vs L910)
- Two W computation methods (L488-500 vs L623-654)
- Local constants should be module-level
- Duplicated SIGMA_T
- Temperature floors
- Carbon abundance fallback
- Kappa fallback
- Orphaned functions (get_excite_temp, compute_dc)
- No adaptive damping in grey-lambda
- ED update without validation
- Malformed file parsing risks

## Architecture Recommendations

### Required Layers

1. **io_model**: Parse MODEL, MODEL_SPEC, VADAT, GAMMAS_IN, SN_HYDRO_DATA, F_OSCDAT
2. **conventions**: Single config object for depth ordering, units, DI meaning, DC format
3. **atomic_data**: Parse and validate level energies, weights, ionization edges
4. **ionization**: Stable log-space Saha, complete chain, omitted-stage reporting
5. **excitation**: LTE populations using exact CMFGEN formula, DC from target/LTE ratio
6. **radiation_field**: Grey W, T_rad, T_ion as optional diagnostics
7. **population_builder**: Separate current/next ion populations, abundance closure
8. **cmfgen_roundtrip**: Write synthetic _IN, read back, verify reconstruction
9. **validation**: Abundance closure, charge conservation, finite values, level counts
10. **reporting**: JSON diagnostics, text summaries, per-species/per-depth reports

### Population Strategies

- **Mode A (default)**: Conservative LTE at T_local, DC=1 (or SET_DC_LTE_V2 low-T correction), no grey, no smoothing
- **Mode B (opt-in)**: Grey dilute-LTE diagnostic, T_ion in _IN, DC=W, auto-fallback to Mode A
- **Mode C**: Previous-model continuation (validate species, grid, levels, abundances)

### Critical Fix Priority

1. Fix charge conservation in grey-lambda (B1, B2)
2. Make T_local the default T_write (B5)
3. Require F_OSCDAT and VADAT, fail if missing (B3, B4)
4. Remove orphaned code (H1)
5. Make smoothing opt-in (H4)
6. Unify DC caps (H3)
7. Add round-trip validation
8. Add JSON diagnostics

## Unresolved Assumptions Requiring Confirmation

1. Does REGRID_LOG_DC_V1 use the species _IN T column or T_IN for LTEPOP_WLD_V2?
   - Evidence suggests T_IN is used for LTE, species _IN T is only for interpolation reference
   - Needs confirmation by tracing the full call path in set_new_model_estimates.f

2. Does the species _IN ED column override GAMMAS_IN?
   - REGRID_LOG_DC_V1 reads ED from _IN, but set_new_model_estimates.f may use GAMMAS_IN ED
   - Needs confirmation

3. How does F_TO_S handle the DC values written by the proxy?
   - regrid_log_dc_v1.f L609-619: assigns DC of lowest state to all levels in superlevel
   - But sup_to_full_v3b.f uses LTE-weighted interpolation for full→super mapping
   - The proxy writes NS values — are these superlevel DCs or full-level DCs?

4. What happens when NOLD < N (proxy writes fewer levels than model expects)?
   - regrid_log_dc_v1.f L557-602: extrapolates using excitation temperature
   - Needs verification that this doesn't cause issues
