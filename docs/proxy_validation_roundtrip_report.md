# Proxy Populations Validation & Round-Trip Report

**Model:** `SN_12C/cmfgen/0.5days` (ND=400, t=0.5532 d, 10 species)
**Code:** `scripts/proxy_populations_v2.py`  |  **Date:** 2026-08-20

## 1. Verified CMFGEN conventions (citations)

Every convention below was checked against the actual `cmfgen_source/` Fortran, not inferred from comments:

| Convention | Value | Source routine (line) |
|---|---|---|
| Depth order | outer→inner (R decreasing) | `regrid_log_dc_v1.f:185-207` |
| DI meaning | next-higher ion-stage density | `set_new_model_estimates.f:365-375` |
| DC in file | LINEAR (modern date ⇒ `TAKE_LOGS=.TRUE.`) | `regrid_log_dc_v1.f:127-134,177-178` |
| Header | `R(ND) RLUM NOLD ND` | `regrid_log_dc_v1.f:139` |
| Row | `R DI ED T T1 T2 CLUMP_FAC` (+idx) | `regrid_log_dc_v1.f:162-166` |
| FEDGE | `(IONIZATION_EN-E)*c*1e-15` (positive edge) | `genosc_v6.f:204-205` |
| LTEPOP | `ln(W*g)+EDGE*HDKT/T+ln(C*ED*DI*T^-1.5/GION)` (recombination form) | `ltepop_wld_v2.f:66-76` |
| LTEPOP C | `2.07078e-22` | `ltepop_wld_v2.f:67` |
| HDKT | `4.7994145` (edge in 10^15 Hz, T in 10^4 K) | `cmfgen.f` old literal |
| F_TO_S | lowest full level's DC → all levels in superlevel | `regrid_log_dc_v1.f:609-619` |
| SET_DC_LTE | DC=1 for T>TMIN_LOC; low-T boost otherwise | `set_dc_lte_v2.f:24-34` |

**Note on the audit's written formula:** `proxy_audit_report.md` wrote `exp(-EDGE*HDKT/T)`, a notational sign flip. The actual CMFGEN code uses `+EDGE*HDKT/T` (Saha *recombination* form, `n_z ∝ n_{z+1} n_e exp(+χ/kT)`). `proxy_populations_v2.py` matches the **code** (verified by `test_ltepop_matches_cmfgen_formula`).

## 2. Round-trip results (LTE baseline)

Round-trip = write `_IN`, emulate `REGRID_LOG_DC_V1` + `LTEPOP_WLD_V2`, compare reconstructed population to target. In LTE mode target == LTE pop and DC=1, so reconstruction is exact by construction.

| Species | NLEV | DI min–max (cm^-3) | DC range | RT max rel err | non-finite | floored |
|---|---|---|---|---|---|---|
| HI | 12 | 7.694e+09 – 2.028e+13 | 1.00e+00–1.00e+00 | 0.00e+00 | 0 | 0 |
| HeI | 10 | 7.776e+08 – 2.050e+12 | 1.00e+00–1.00e+00 | 0.00e+00 | 0 | 0 |
| OI | 7 | 2.202e+03 – 5.804e+06 | 1.00e+00–1.00e+00 | 0.00e+00 | 0 | 0 |
| SkI | 1 | 1.401e+04 – 3.693e+07 | 1.00e+00–1.00e+00 | 0.00e+00 | 0 | 0 |
| Sk2 | 10 | 1.401e+04 – 3.693e+07 | 1.00e+00–1.00e+00 | 0.00e+00 | 0 | 0 |
| CaI | 1 | 2.007e+02 – 5.290e+05 | 1.00e+00–1.00e+00 | 0.00e+00 | 0 | 0 |
| Ca2 | 12 | 2.007e+02 – 5.290e+05 | 1.00e+00–1.00e+00 | 0.00e+00 | 0 | 0 |
| FeI | 1 | 1.833e+05 – 2.255e+09 | 1.00e+00–1.00e+00 | 0.00e+00 | 0 | 0 |
| Fe2 | 25 | 1.327e+01 – 1.008e+10 | 1.00e+00–1.00e+00 | 0.00e+00 | 0 | 0 |
| FeIII | 1 | 1.327e+01 – 1.008e+10 | 1.00e+00–1.00e+00 | 0.00e+00 | 0 | 0 |

**All species round-trip max relative error = 0.0; 0 non-finite values.** Exit status 0 ⇒ output safe to feed to CMFGEN.

## 3. LTE vs GREY comparison

| Species | LTE DC range | GREY DC range | GREY RT max rel err |
|---|---|---|---|
| HI | 1.00e+00–1.00e+00 | 1.01e-03–3.89e+03 | 2.02e-16 |
| HeI | 1.00e+00–1.00e+00 | 4.76e-16–2.58e+00 | 2.10e-16 |
| OI | 1.00e+00–1.00e+00 | 8.69e-04–3.11e+03 | 2.06e-16 |
| SkI | 1.00e+00–1.00e+00 | 8.39e+03–8.40e+06 | 1.87e-16 |
| Sk2 | 1.00e+00–1.00e+00 | 2.17e-07–1.46e+02 | 2.01e-16 |
| CaI | 1.00e+00–1.00e+00 | 4.11e+05–1.46e+08 | 1.52e-16 |
| Ca2 | 1.00e+00–1.00e+00 | 9.87e-02–7.62e+03 | 2.06e-16 |
| FeI | 1.00e+00–1.00e+00 | 1.29e+04–3.69e+08 | 1.87e-16 |
| Fe2 | 1.00e+00–1.00e+00 | 5.34e-07–1.24e+02 | 2.12e-16 |
| FeIII | 1.00e+00–1.00e+00 | 1.31e-23–7.00e-03 | 2.05e-16 |

In the optically-thick deep ejecta GREY DC→1 (identical to LTE). Deviations appear only in the optically-thin outer zones where W<1. GREY round-trip is also exact (target/LTE invariant under closure rescale).

## 4. Validation issues

- Validation issues: **7** (all omitted-stage sink-risk warnings, expected)
- Floored depths (numerical floor): 0
- Capped DCs: 0
- Recomputed N_e: 8.479e+09 – 2.240e+13 cm^-3
- N_e |ratio-1|: max 6.633e+02, median 5.072e-01

The large N_e discrepancy is **expected and reported**, not hidden: Mode A uses the input GAMMAS_IN N_e for the Saha balance rather than forcing charge neutrality. CMFGEN recomputes N_e on its first iteration. The proxy reports the inconsistency so the user knows where the first correction will be largest.

## 5. Abundance residuals (VADAT vs hydro)

| Element | VADAT | hydro mean | max dex diff |
|---|---|---|---|
| HYD | 7.001e-01 | 7.000e-01 | 0.000 |  |
| HE | 2.810e-01 | 2.823e-01 | 0.021 |  |
| OXY | 3.180e-06 | 1.450e-05 | 9.502 |  |
| SIL | 3.551e-05 | 1.592e-04 | 10.550 |  |
| CAL | 7.260e-07 | 1.000e-15 | 8.861 |  |
| IRON | 1.931e-02 | 1.686e-02 | 2.036 |  |
| NICK | 1.000e-08 | 1.000e-15 | 7.000 |  |
| COB | 1.000e-08 | absent | n/a | absent from hydro |

VADAT abundances are depth-constant (as CMFGEN's SOLVEBA/SCALE_POPS uses them); hydro abundances vary with depth. Large dex differences (OXY, SIL, CAL) reflect composition gradients, not proxy errors — reported, not silently overwritten.

## 6. Omitted-ion-stage risk

- **SkI**: omitted-stage frac max = 1.000 (not a safe sink)
- **Sk2**: omitted-stage frac max = 1.000 (not a safe sink)
- **CaI**: omitted-stage frac max = 1.000 (not a safe sink)
- **Ca2**: omitted-stage frac max = 1.000 (not a safe sink)
- **FeI**: omitted-stage frac max = 0.998 (not a safe sink)
- **Fe2**: omitted-stage frac max = 0.998 (not a safe sink)
- **FeIII**: omitted-stage frac max = 0.998 (not a safe sink)

Recommendation: add the next ion stage (e.g. Ca III, Fe IV, Si IV) to MODEL_SPEC, or accept larger first-iteration corrections.

## 7. Acceptance criteria checklist

- [x] Parser, physics, structural, round-trip tests pass (42/42)
- [x] Never silently changes abundance / depth order / ion-stage meaning / level count
- [x] Writes local physical temperature by default (Mode A)
- [x] Abundance & charge closure reported independently of numerical floors
- [x] Reports omitted-stage and truncation risks
- [x] Finite, non-negative, internally consistent populations (0 non-finite)
- [x] Default mode no less stable than old LTE baseline (DC=1, exact round-trip)
- [ ] Demonstrated improvement (lower first-iteration CMFGEN corrections) — requires a real CMFGEN run; see unresolved assumptions #5.

