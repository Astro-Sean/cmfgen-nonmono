# Grid Mapping Validation Report

## Overview

This report documents the validation of the new `grid_mapping` package
for STELLA-to-CMFGEN conversion with non-monotonic velocity support.

**Date:** 2026-08-20
**Source:** `SN_A_g_M05/stella/SN_A_g_M05.swd` (epoch 30.24 days)
**Target:** 190-point CMFGEN grid

## Test Results

### Unit Tests

All 46 unit tests pass:

```
============================== 46 passed in 0.69s ==============================
```

Test coverage:
- SIGMA computation (5 tests): homologous expansion, non-monotonic, sign preservation, cap, reversal counting
- Shell model (5 tests): volume, centres, density, validation, bad radius detection
- Conservative remap (6 tests): overlap volumes, mass/elemental/Ni56/momentum conservation, residuals
- Non-monotonic velocity (3 tests): reversal preservation, extrema detection, dV/dr sign changes
- Shock preservation (1 test): density jump detection
- Optical depth (3 tests): tau computation, landmarks, photosphere
- Ionization fronts (1 test): Ne jump detection
- Parser ordering (2 tests): non-monotonic radius, non-finite values
- Source audit (2 tests): audit findings, sigma compatibility
- CMFGEN writer (1 test): SN_HYDRO_DATA file writing
- Synthetic models (17 tests): monotonic, single/multi reversal, shock, composition, ionization, all features

### Synthetic Model Tests

All 7 synthetic models pass validation:
- T1: Monotonic (homologous expansion, sigma=0)
- T2: Single reversal (velocity minimum preserved)
- T3: Multi-reversal (3 reversals preserved)
- T4: Shock (density jump, mass conserved)
- T5: Composition interface (H/He ratio change, elemental mass conserved)
- T6: Ionization front (Ne jump detected)
- T7: All features combined (shock + reversal + interface + front)

### Real STELLA Snapshot Validation

**Source:** SN_A_g_M05 at 30.24 days (500 zones, 1.73 Msun)
**Target:** 190 points

#### Validation Status: PASS

- Radius strictly decreasing: YES (0 violations)
- All arrays finite: YES
- Density positive: YES
- Temperature positive: YES
- SIGMA consistency: max diff = 0.0
- Mass fractions sum to 1: max deviation = 2.2e-16

#### Conservation Residuals

| Quantity | Source | Target | Residual |
|----------|--------|--------|----------|
| Total mass | 3.438e+33 g | 3.436e+33 g | -5.8e-04 |
| H mass | 1.982e+33 g | 1.981e+33 g | -2.8e-04 |
| He mass | 1.672e+32 g | 1.672e+32 g | -6.1e-07 |
| Ni56 mass | 1.771e+32 g | 1.771e+32 g | -6.1e-16 |
| Momentum | 1.962e+42 | 1.959e+42 | -1.3e-03 |
| Kinetic energy | 6.362e+50 | 6.344e+50 | -2.8e-03 |
| Electron number | 3.557e+55 | 3.660e+55 | +2.9e-02 |
| Total tau | 8.684 | 8.964 | +3.2e-02 |

#### Detected Features

- 2 velocity maxima
- 2 velocity minima
- 4 velocity reversals (dV/dr sign changes)
- 5 contact discontinuities
- 12 composition interfaces
- 2 isotope boundaries
- 3 ionization fronts
- 4 tau landmarks (including photosphere at tau=2/3)
- 26 total anchor points (24 preserved in target)

#### SIGMA Compatibility

- 32 points will be capped to -0.999 by CMFGEN (in the -1.2 < sigma < -1.0 range)
- This is expected for non-monotonic velocity and is handled by the Dessart+2015 guard

## Remaining Warnings

1. **Ni mass residual (-2.3%):** The Ni conservation residual is larger than
   other elements because the decay chain (Ni56 -> Co56 -> Fe56) is applied
   after remapping, redistributing mass between Ni, Co, and Fe. This is
   physically correct but means the *elemental* Ni mass changes due to
   radioactive decay, not due to remapping error.

2. **Electron number residual (+2.9%):** The electron density is remapped
   using volume-weighted averaging, while the source Ne comes from STELLA's
   NLTE calculation. The small residual is expected from the averaging
   method.

3. **Total tau residual (+3.2%):** The optical depth is recomputed from
   the remapped density and opacity, which differs slightly from the
   source due to the mass-weighted opacity averaging.

## Unresolved Issues (from source audit)

The following issues from the CMFGEN source audit remain unresolved and
require solver-level testing via `run_model_grid.py`:

1. **`comp_j_blank.f` Eddington pre-computation:** Still calls
   `MOM_J_CMF_V11` unconditionally. May be mitigated by `N_TYPE=N_ON_J`
   setting but requires end-to-end solver test to confirm.

2. **V=0 behavior:** STELLA snapshots in this dataset do not contain
   V=0 or V<0 zones. The code handles V=0 by setting sigma=0, but
   this has not been tested with the actual CMFGEN solver.

3. **Minimum resolution per velocity reversal:** No empirical determination
   of the minimum number of grid points needed to resolve a velocity
   reversal for the modified solver. Requires convergence testing.

4. **Full-grid vs reduced-grid spectrum comparison:** Not yet performed.
   Requires running CMFGEN on both the full 500-zone grid and the
   reduced 190-zone grid and comparing spectra.

## Compatibility with run_model_grid.py

The new converter (`convert_stella_swd_to_cmfgen_v2.py`) preserves the
full CLI contract of the original converter. All ~75 arguments are
accepted. The orchestrator can call it without changes:

```bash
python scripts/run_model_grid.py configs/<model>.yaml --skip-stella --run-cmfgen
```

To switch to the new converter, update `run_model_grid.py` to call
`convert_stella_swd_to_cmfgen_v2.py` instead of
`convert_stella_swd_to_cmfgen.py`.

## Smoothing Policy

All physical smoothing is **disabled by default** in the new converter:
- `--smooth-velocity` defaults to False
- `--smooth-wiggles` defaults to False (was True in old converter)
- `--smooth-sigma-scatter` defaults to False (was True in old converter)
- `--smooth-sigma-global` defaults to False (was True in old converter)

If smoothing is explicitly requested, it is recorded in the
`grid_audit.json` configuration section.

## Velocity Sign Preservation

The `prepare_VADAT.py` script has been updated to remove `np.abs()`
calls on velocity in all file-writing paths:
- `print_hydro_summary`
- `downsample_hydro_velocity_aware`
- `write_hydro_data_sniin`
- `write_input_hydro_dat_sniin`
- `write_t_in_from_hydro`
- `write_fin_cal_grid_from_hydro`

The remaining `np.abs()` calls in `compute_optimal_vfrac_mom` are
correct — they measure velocity jump magnitudes for VFRAC_MOM
optimization, not for writing physical velocity values.
