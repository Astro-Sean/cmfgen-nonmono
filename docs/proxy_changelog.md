# Changelog: proxy_populations.py → proxy_populations_v2.py

**Date:** 2026-08-20
**Old:** `scripts/proxy_populations.py` (~2922 lines)
**New:** `scripts/proxy_populations_v2.py` (~1509 lines, 10 explicit layers)
**Audit:** `scripts/proxy_audit_report.md` (67 issues: 5 BLOCKING, 18 HIGH, 28 MEDIUM, 13 LOW)

Every behaviour change is listed below with the audit issue ID it resolves.
Unchanged physical constants (HDKT, LTEPOP_C, SAHA_CONST) are NOT behaviour changes.

## BLOCKING fixes (physical-model-affecting)

| ID | Old behaviour | New behaviour |
|---|---|---|
| B1 | Two Saha implementations: `compute_saha_ion_frac()` (no renormalization) and an inline copy in `iterate_grey_lambda()` (renormalizes to sum=1). Charge-conservation drift between grey iteration and final DI. | One log-space Saha solver (`solve_saha_chain`) using `scipy.special.logsumexp`; fractions always sum to 1 (pre-floor). No duplicated formula. |
| B2 | Grey-lambda convergence checked only relative Δne and Δκ; no absolute charge-conservation residual. ED could converge to a value inconsistent with the ionization state. | `validate_run` recomputes N_e from the physical (pre-floor) ion fractions and reports `ne_ratio_max` / `ne_ratio_median` vs input GAMMAS_IN ED. |
| B3 | Missing `F_OSCDAT` silently fell back to hard-coded `IONIZATION_ENERGIES_CM_FALLBACK`, which may not match the atomic data CMFGEN uses. | Missing `F_OSCDAT` skips the species with an explicit `"skipped"` strategy-log entry; no silent fallback to invented ionization energies. |
| B4 | Missing VADAT silently fell back to depth-varying hydro abundances, causing ~1e28 mismatch at composition boundaries (CMFGEN uses depth-constant VADAT abundances in SOLVEBA). | VADAT abundances are required (depth-constant); all VADAT-vs-hydro differences are reported in `abundance_diff_report` / `proxy_validation.txt`. |
| B5 | Default mode wrote `T_ion` (radiation temperature) into the `_IN` T column. | Default (Mode A) writes **T_local** (physical gas temperature). `T_ion` is opt-in via `--write-t-ion` and only meaningful in GREY mode. |

## HIGH fixes (consistency / numerical)

| ID | Old behaviour | New behaviour |
|---|---|---|
| H1 | Four DC methods with inconsistent exponent signs, floors, and caps (1e10 vs 1e30): `compute_dc` (orphaned), `compute_dc_dilute_lte`, `compute_dc_lte_cmfgen`, `compute_dc_nlte_reduced`. | One `dc_lte_cmfgen` (Mode A, DC=1 / SET_DC_LTE_V2 low-T boost) and one `dc_from_target` (Mode B/C, target/LTE with a single configurable cap). Orphaned `compute_dc` removed. |
| H2 | Statistical weight defaulted to `gion_default=1.0` for unknown species, biasing Saha by factor g_ion. | `gion` read from F_OSCDAT/MODEL_SPEC; missing g_ion falls back to `sp.gion` and is recorded, not silently set to 1. |
| H3 | DC caps inconsistent (1e10 vs 1e30). | Single `dc_max=1e300` in `dc_from_target` (effectively uncapped; caps are reported via `capped_mask`). |
| H4 | 3-point median filter on `ion_frac` applied by default; could erase real ionization fronts. | Smoothing is **opt-in** (`--smooth`); off by default; all smoothed depths would be reported. |
| H5 | Constants `SAHA_CONST`/`CMFGEN_C` unverified. | Verified against source: `LTEPOP_C=2.07078e-22` = `ltepop_wld_v2.f:67`; `SAHA_CONST=2.414e15` = standard cgs Saha constant. |
| H6 | Hard-coded Kramers coefficient `3.7e22` with no origin. | Documented; grey opacity is diagnostic-only and not on the default path. |
| H7 | Multiple g_ion fallback paths in `build_saha_atomic_data`. | Single fallback to `sp.gion`, recorded. |
| H8 | Renormalization inconsistency between the two Saha implementations. | Eliminated (single solver, see B1). |
| H9–H18 | Various silent floors/fallbacks. | Floors are opt-in (`--min-ion-frac`, default 1e-30), reported per species, and **excluded** from abundance/charge closure. |

## MEDIUM fixes (architecture / hygiene)

- **M1 (level count):** old code wrote an ambiguous level count; new code defaults to **NF** (full levels) and pre-applies F_TO_S replication so the file is self-consistent regardless of DC_METH. `--write-levels NS` is available for the legacy RTX path.
- **M2:** SPEED_LIGHT precision unified (one module constant).
- **M3:** Two W computation methods collapsed into `grey_w_eff`.
- **M4:** No CARB abundance fallback invented; caller decides and differences are reported.
- **M5:** Fragile regex transition parser removed (transitions unused by the population build).
- **M6:** GAMMAS_IN is copied atomically to the output dir; source never overwritten without `--force-overwrite-source`.
- Removed orphaned `get_excite_temp` / `compute_dc` dead code.
- All `_IN` writes are atomic (`os.replace` of a `.tmp` file).

## Architecture (new)

Refactored into 10 explicit, testable layers: `io_model`, `conventions`,
`atomic_data`, `ionization`, `excitation`, `radiation_field`,
`population_builder`, `cmfgen_roundtrip`, `validation`, `reporting`.
Each layer has a single responsibility and is covered by
`scripts/test_proxy_populations_v2.py` (42 tests).

## New diagnostics (written every run)

- `proxy_audit.json` — machine-readable: conventions, per-species DI/DC/round-trip,
  abundance diff, validation.
- `proxy_validation.txt` — human-readable: per-species summary, omitted-stage risk,
  recomputed N_e vs input, abundance residuals.
- `proxy_strategy.json` — selected strategy + abundance provenance.

## Selection logic (new, conservative)

1. Validated previous CMFGEN model present ⇒ CONTINUATION (Mode C).
2. CMFGEN conventions cannot be verified ⇒ **stop**, no production output.
3. Otherwise ⇒ LTE baseline (Mode A).
4. GREY (Mode B) only if it reduces predicted corrections without violating
   closure/charge/correction limits.
5. Omitted-stage fraction above threshold ⇒ warn + exclude from NLTE enhancement.
6. Any round-trip/structural test fails ⇒ non-zero exit, no production output.

## CLI changes

- Positional `model_dir` (was `--model-dir`).
- Added: `--strategy {LTE,GREY,CONTINUATION}`, `--write-t-ion`, `--smooth`,
  `--write-levels {NF,NS}`, `--min-ion-frac`, `--t-exc-min`, `--tau-trans`,
  `--prev-model`, `--force-overwrite-source`.
- Exit status: 0 = round-trip OK; 2 = round-trip failed (do not use output).

## Unchanged by design

- Physical constants (HDKT, LTEPOP_C, SAHA_CONST, SPEED_LIGHT) — verified equal.
- Depth ordering (outer→inner), DI = next-higher stage, DC linear in file,
  R in 10^10 cm, T in 10^4 K — all preserved and now explicitly configured.
- The proxy remains a **preconditioner**; no claim of NLTE accuracy is made.
