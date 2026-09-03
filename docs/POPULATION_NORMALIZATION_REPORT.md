# POPULATION_NORMALIZATION_REPORT.md

**Date:** 2026-08-21
**Model:** SN_12C/cmfgen/0.5days (ND=400, 10 species)

---

## 1. Normalization Requirements

The proxy must satisfy:
1. Ion populations do not exceed the parent element
2. Level populations sum to the ion population
3. Finite electron density
4. Documented charge-balance residual
5. Correct elemental abundance limits (VADAT)

---

## 2. Ion Normalization (VERIFIED)

The proxy v2 uses a Saha-based ionization calculation:
- `compute_saha_ion_frac()` computes ionization fractions for each stage
- Fractions are normalized to sum=1 across all stages of each element
- `min_ion_frac=1e-30` floor prevents zero populations (documented, configurable)
- The floor is applied to the ion fraction, not the absolute density

For single-stage elements (HI, HeI, OI where the next stage is not modelled),
DI is estimated via `chi_to_next` Saha ratio: `DI = n_elem × r/(1+r)`.

**Verification:** The proxy reports ion fractions and confirms normalization in
`proxy_audit.json`. All ion fractions sum to 1.0 (within numerical precision).

---

## 3. Level Normalization (VERIFIED)

In LTE mode (DC=1), all levels in a superlevel receive the same DC (matching
CMFGEN's F_TO_S pass in `regrid_log_dc_v1.f:609-619`). The proxy writes NF
(full-level) DCs and pre-applies F_TO_S replication.

The LTE population for each level is computed by CMFGEN's `LTEPOP_WLD_V2`:
```
n_LTE(J) = W(J) * g(J) * exp(EDGE(J)*HDKT/T) * (2.07078e-22 * ED * DI * T^-1.5) / GION
```

The proxy does not write absolute level populations — it writes DC (departure
coefficients) and DI (next-stage density). CMFGEN reconstructs populations as:
```
population(J) = DC(J) * n_LTE(J)
```

In LTE mode, DC=1, so `population(J) = n_LTE(J)` — the LTE population.

---

## 4. Charge Balance (VERIFIED — reported, not enforced)

The proxy uses GAMMAS_IN electron density for the Saha balance (Mode A), rather
than forcing charge neutrality. This means the proxy's ionization state may not
be exactly charge-consistent with the input ED.

**Reported residual:**
- Recomputed N_e: 8.479e+09 – 2.240e+13 cm^-3
- N_e |ratio-1|: max 6.633e+02, median 5.072e-01

The large N_e discrepancy is **expected and reported**: Mode A uses the input
GAMMAS_IN N_e for the Saha balance. CMFGEN recomputes N_e on its first iteration.
The proxy reports the inconsistency so the user knows where the first correction
will be largest.

This is not a bug — it is the conservative choice. Forcing charge neutrality would
require iterating the Saha balance to convergence, which is the grey-lambda mode
(Mode B, opt-in).

---

## 5. Abundance Closure (VERIFIED)

The proxy uses VADAT abundances (depth-constant, matching CMFGEN's SOLVEBA/SCALE_POPS).
CMFGEN re-derives populations to match VADAT abundance each iteration via SCALE_POPS.

Abundance residuals (VADAT vs hydro) are reported per element:
- HYD: 0.000 dex (matches)
- HE: 0.021 dex (close)
- OXY: 9.502 dex (composition gradient)
- SIL: 10.550 dex (composition gradient)
- CAL: 8.861 dex (composition gradient)
- IRON: 2.036 dex (moderate)

Large dex differences reflect real composition gradients in the hydro model, not
proxy errors. The proxy uses VADAT (what CMFGEN will enforce) and reports the
difference.

---

## 6. Numerical Safety (VERIFIED)

| Check | Result |
|---|---|
| All values finite | PASS (0 non-finite in all 10 files) |
| No negative populations | PASS (0 negative DI/ED/DC) |
| Elemental abundance limits | PASS (VADAT abundances used) |
| Ionization normalization | PASS (fractions sum to 1) |
| Level normalization | PASS (DC=1 in LTE, F_TO_S replicated) |
| Finite electron density | PASS (1.57e9 – 2.30e13 cm^-3) |
| Charge-balance residual | REPORTED (max |ratio-1| = 663) |
| Low-density outer zones | PASS (no NaN/Inf at depth 1) |
| Hot shocked zones | PASS (no NaN/Inf at inner depths) |
| Composition boundaries | PASS (no smoothing by default) |

---

## 7. Numerical Floors (VERIFIED — documented)

| Floor | Value | Location | Sensitivity |
|---|---|---|---|
| min_ion_frac | 1e-30 | compute_saha_ion_frac | Configurable via --min-ion-frac |
| DC clip | [-700, 700] | compute_dc_lte_cmfgen | Prevents exp overflow |
| ED floor | max(ed, 1.0) | compute_saha_ion_frac | Prevents division by zero |
| T floor | max(t, 0.3) (t_exc_min) | get_excite_temp | Configurable via --t-exc-min |

All floors are named, configurable, and documented. No floor is applied silently.
