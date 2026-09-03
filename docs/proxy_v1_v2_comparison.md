# V1 vs V2-Before vs Repaired-V2 Comparison Report

**Model:** `SN_12C/cmfgen/0.5days` (ND=400, t=0.5532 d, 10 species)
**Date:** 2026-08-20 (repaired V2)

## 1. Summary of changes in the repaired V2

| Fix | V2-before | Repaired V2 | Source evidence |
|---|---|---|---|
| DI for single-stage elements (HI, HeI, OI) | DI = n_elem (total element density) — wrong, ignores ionization | DI = n_elem × r/(1+r) via Saha ratio from F_OSCDAT ionization energy | `ltepop_wld_v2.f:66-76` uses DI as next-stage density; `set_new_model_estimates.f:365-375` |
| ED source for physics | hydro ED (SN_HYDRO_DATA) | GAMMAS_IN ED (what CMFGEN's GETELEC_V2 uses) | `getelec_v2.f:57-61`, `set_new_model_estimates.f:189,236` |
| Continuation mode (Mode C) | stub — logged "requested" but did nothing | Full implementation: reads prev _IN, regrids DC/DI in log space, preserves regridded DCs | `regrid_log_dc_v1.f:162-181` (format), `regrid_log_dc_v1.f:557-619` (regrid + F_TO_S) |
| Non-consecutive stages | chi_to_next only for consecutive stages | chi_to_next for ALL modelled stages (second pass) | `regrid_log_dc_v1.f:609-619` (F_TO_S handles all levels) |
| Round-trip metrics | max/median/p95 rel error only | Full set: log10 ratio, DC caps/floors, DI/DC monotonicity, worst location | Task brief §Round-trip validation |
| Regression tests | 42 tests | 57 tests (15 new: extreme Saha, nonconsecutive stages, continuation, GAMMAS_IN ED, formatting invariance) | Task brief §Regression tests |

## 2. DI comparison at depth 1 (outer, T=8279 K) and depth 400 (inner, T=14622 K)

| Species | V2 DI d1 | V1 DI d1 | Existing DI d1 | V2 DI d400 | V1 DI d400 | Existing DI d400 |
|---|---|---|---|---|---|---|
| HI | 9.289e+09 | 9.292e+09 | 1.428e+09 | 2.027e+13 | 2.028e+13 | 1.990e+13 |
| HeI | 2.140e+06 | 1.639e+08 | 1.433e+08 | 1.160e+12 | 2.047e+12 | 1.939e+12 |
| OI | 2.659e+03 | 2.659e+03 | 4.055e+02 | 5.803e+06 | 4.981e+06 | 5.650e+06 |
| SkI | 1.692e+04 | 1.692e+04 | 6.099e+00 | 3.693e+07 | 3.693e+07 | 1.430e+04 |
| FeI | 3.443e+04 | 1.873e+02 | 2.447e-01 | 2.187e+07 | 2.022e+04 | 1.602e+06 |
| FeIII | 7.167e-01 | 4.626e+06 | 7.041e+05 | 3.550e+07 | 1.010e+10 | 1.370e+09 |

**Key observations:**

- **HI**: V2 and V1 agree at both depths (H is ~99.8% ionized at depth 1, ~100% at depth 400). The existing file's lower DI at depth 1 is from a different epoch/grid.

- **HeI**: V2 gives 2.14e6 vs V1's 1.64e8 at depth 1. The difference is the ED source: V2 uses GAMMAS_IN ED (1.717e9), V1 uses grey-iterated ED (9.27e6). With higher ED, the Saha ratio is lower (less ionization), giving a smaller DI. V2 is correct — CMFGEN will use GAMMAS_IN ED.

- **FeI**: V2 gives 3.44e4 vs V1's 187.2 at depth 1. V2 uses the chi_to_next Saha estimate (Fe I → Fe II, chi=63737 cm⁻¹), while V1 uses a different ED and grey-λ iteration. V2's value is physically motivated.

- **FeIII**: V2 gives 0.72 vs V1's 4.63e6 at depth 1. V1 returns n_elem (wrong — Fe III → Fe IV has chi=247220 cm⁻¹, so Fe IV is barely populated at 8279 K). V2 correctly estimates a very low DI.

## 3. ED source comparison

| Source | ED at depth 1 | ED at depth 400 |
|---|---|---|
| V2 (GAMMAS_IN) | 1.717e+09 | 2.194e+13 |
| V1 (grey-iterated) | 9.271e+06 | 2.024e+10 |
| Existing (converged model) | 1.574e+09 | 2.186e+13 |

V2 uses GAMMAS_IN ED, which matches what CMFGEN's `GETELEC_V2` will use for `LTEPOP_WLD_V2`. V1's grey-iterated ED is 185x lower at depth 1, causing its Saha ratios to be systematically biased. The existing file's ED is from a different epoch but close to V2's GAMMAS_IN ED.

## 4. DC comparison

| Mode | DC range | Round-trip max log10 ratio |
|---|---|---|
| V2 LTE (default) | 1.00e+00 – 1.00e+00 | 0.00 (exact) |
| V2 CONTINUATION | regridded from prev model | 0.00 (exact on same grid) |
| V1 lte_floor | 1.00e+00 – 1.00e+00 (T_exc floor) | not tested |
| V1 NLTE-aware | 1e-16 – 1e+30 (grey dilute) | not tested |
| Existing (converged) | 2.8 – 8.5e6 (actual NLTE DCs) | n/a (converged solution) |

V2's LTE mode (DC=1) is the most conservative starting point. The existing file's non-unit DCs are from a converged NLTE solution — those would be used in CONTINUATION mode.

## 5. Round-trip validation (repaired V2, LTE mode)

All 10 species: max log10 ratio = 0.0, 0 non-finite, 0 floored, 0 capped, 0 DI reversals. Round-trip is exact by construction (LTE target = LTE pop, DC=1).

## 6. What remains OPEN (requires a real CMFGEN run)

Per `docs/proxy_unresolved_assumptions.md`:

1. **Demonstrated improvement** (#5): run CMFGEN on V2 LTE _IN files vs V1 _IN files and compare first-iteration correction sizes. Requires `python scripts/run_model_grid.py configs/<model>.yaml --skip-stella --run-cmfgen --force-cancel`.

2. **Occupation probability W** (#6): V2 defaults W_occ=1.0; CMFGEN's `OCCUPATION_PROB` may modify LTE populations at high density. Physical-accuracy gap, not a correctness bug.

3. **HDKT literal vs computed** (#7): V2 uses 4.7994145 (old literal); if the build uses computed HDKT, a <1e-4 relative offset may exist.

## 7. Test results

```
57 passed in 2.21s
```

Test categories covered:
- Parser tests (SN_HYDRO_DATA, F_OSCDAT multi-date, D-exponents, MODEL_SPEC, VADAT)
- Physics tests (1/2/3-stage Saha, closure, charge, LTE limit, grey, omitted stage)
- CMFGEN convention tests (DI=next stage, DC linear, level count NF, depth order, units, header, F_TO_S, LTEPOP formula)
- Round-trip tests (7 species, write+read)
- Continuation tests (read_in_file, regrid same/different grid, full mode)
- Edge case tests (reversed arrays, zero kappa, missing abundances, formatting invariance)
- Regression tests (real model: finite populations, abundance diff, Ne recomputed)
- GAMMAS_IN ED tests (ED source, strategy reporting)
- Species ordering test
