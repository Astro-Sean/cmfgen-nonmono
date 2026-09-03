# PROXY_POPULATIONS_AUDIT.md

**Date:** 2026-08-21
**Scope:** Full audit of the CMFGEN proxy populations workflow in the SN2024pba project.
**Evidence labels:** Each conclusion is tagged VERIFIED / SUPPORTED / UNVERIFIED / CONTRADICTED / BLOCKED.

---

## 1. Executive Summary

The proxy populations workflow uses `scripts/proxy_populations_v2.py` (a refactored,
testable version of the original `scripts/proxy_populations.py`). The v2 script
generates CMFGEN `_IN` population files using a conservative LTE baseline (Mode A,
DC=1) by default, with optional grey dilute-LTE (Mode B) and continuation (Mode C)
strategies.

**Key findings:**

| Finding | Status |
|---|---|
| _IN file format conventions (DI=next stage, linear DC, outer-to-inner, units) | VERIFIED (source-traced) |
| Round-trip serialization (write → independent parse → reconstruct) | VERIFIED (exact, 0 error) |
| Independent parser validates all 10 _IN files | VERIFIED (all pass) |
| MODEL_SPEC NF matches proxy NOLD for all species | VERIFIED |
| Atomic data manifests built and checksummed (65 files) | VERIFIED |
| Input grid identity (ND=400, outer-to-inner) | VERIFIED |
| Reference URL accessible; 34 files downloaded with checksums | VERIFIED |
| Reference model atom NOT comparable (different ND, species, epoch, version) | VERIFIED |
| Orchestrator integration incomplete (proxy not wired through run_model_grid.py) | VERIFIED |
| CMFGEN startup test (confirm proxy files read) | IN PROGRESS |
| Fresh-model convergence comparison | BLOCKED (requires full CMFGEN run) |

---

## 2. Existing Prior Work

The following prior work was found in the repository (all untracked, never committed):

- `scripts/proxy_populations.py` (v1, 133KB) — original implementation
- `scripts/proxy_populations_v2.py` (101KB) — refactored, repaired version
- `scripts/test_proxy_populations_v2.py` (45KB, 57 tests) — all passing
- `docs/proxy_populations_audit.md` — detailed v1 audit (67 issues, source-verified)
- `docs/proxy_v1_v2_comparison.md` — v1 vs v2 comparison
- `docs/proxy_validation_roundtrip_report.md` — round-trip validation
- `docs/proxy_unresolved_assumptions.md` — 7 assumptions, 3 resolved, 4 open
- `docs/proxy_changelog.md` — change log

The prior work verified all CMFGEN conventions from source inspection:
- DI = next-higher ion-stage density (`set_new_model_estimates.f:365-375`)
- DC = linear in file (modern date → `TAKE_LOGS=.TRUE.`) (`regrid_log_dc_v1.f:127-134`)
- Physical T/ED from T_IN/GAMMAS_IN, not species _IN (`set_new_model_estimates.f:158-163`)
- LTE_EST=T bypasses proxy files entirely (`set_new_model_estimates.f:261-270`)
- LTEPOP formula and constant 2.07078e-22 (`ltepop_wld_v2.f:66-76`)
- F_TO_S: lowest full level's DC → all levels in superlevel (`regrid_log_dc_v1.f:609-619`)

---

## 3. CMFGEN Convention Verification (VERIFIED)

All conventions were verified against `cmfgen_source/` Fortran by the prior audit
and re-confirmed this session. See `docs/proxy_populations_audit.md` §0 for
full citations.

| Convention | Value | Source | Status |
|---|---|---|---|
| Depth order | outer→inner (R decreasing) | `regrid_log_dc_v1.f:185-207` | VERIFIED |
| DI meaning | next-higher ion-stage density | `set_new_model_estimates.f:365-375` | VERIFIED |
| DC in file | LINEAR (modern date) | `regrid_log_dc_v1.f:127-134` | VERIFIED |
| Header | R(ND) RLUM NOLD ND | `regrid_log_dc_v1.f:139` | VERIFIED |
| Row columns | R DI ED T T1 T2 CLUMP_FAC | `regrid_log_dc_v1.f:162-166` | VERIFIED |
| LTEPOP C | 2.07078e-22 | `ltepop_wld_v2.f:67` | VERIFIED |
| HDKT (built binary) | 4.79924307 (computed) | OUTGEN from startup test | VERIFIED |
| HDKT (proxy uses) | 4.7994145 (old literal) | proxy_populations_v2.py | VERIFIED |
| HDKT discrepancy | ~3e-5 relative | — | VERIFIED (negligible) |

**Note:** The CMFGEN startup test OUTGEN shows `HDKT=4.79924307E+00` (computed from
physical constants), while the proxy uses the old literal `4.7994145`. The relative
difference is ~3e-5, which is negligible for LTE populations. This resolves
unresolved assumption #7 from `docs/proxy_unresolved_assumptions.md`.

---

## 4. Independent Parser Validation (VERIFIED)

An independent parser (`tests/test_independent_in_parser.py`) was written that does
NOT import the proxy generator. It parses all 10 generated _IN files and validates:

- All values finite (0 non-finite)
- No negative populations (0 negative)
- Radius monotonically decreasing (outer-to-inner confirmed)
- Consistent ND/R/ED/T across all species files
- NOLD matches MODEL_SPEC NF for all species
- CLUMP_FAC = 1.0 (no clumping)
- Header R(ND) matches innermost radius

**Result:** All 10 files pass all checks. See
`runs/proxy_validation/SN_12C_0.5d/independent_parse_report.json`.

---

## 5. Round-Trip Validation (VERIFIED)

The proxy v2 performs an internal round-trip: write _IN → emulate REGRID_LOG_DC_V1 +
LTEPOP_WLD_V2 → compare reconstructed population to target.

In LTE mode (DC=1), the round-trip is exact by construction:
- max log10 ratio = 0.0 for all 10 species
- 0 non-finite values
- 0 floored/capped DCs
- 0 DI reversals

See `docs/proxy_validation_roundtrip_report.md` for full per-species tables.

---

## 6. Orchestrator Integration (VERIFIED — INCOMPLETE)

**Finding:** The proxy populations workflow is NOT wired through the mandatory
execution path (`run_model_grid.py → proxy_populations.py → CMFGEN`).

Evidence:
1. `run_model_grid.py` reads `proxy_populations` from YAML config and passes it to
   `submit_cmfgen_jobs.sh` (line 4856, 4870).
2. `submit_cmfgen_jobs.sh` passes `--proxy-populations` to `prepare_VADAT.py` (line 137).
3. `prepare_VADAT.py` uses strict `p.parse_args()` (line 2584) and has NO
   `--proxy-populations` argument → would crash with argparse error.
4. No script in the repository calls `proxy_populations_v2.py` or
   `proxy_populations.py` to generate _IN files.

**Consequence:** The proxy works as a standalone tool but cannot be invoked through
the orchestrator. Production runs using proxy populations require manual staging.

**Required fix:** Wire `proxy_populations_v2.py` into the orchestrator:
- Add `--proxy-populations` argument to `prepare_VADAT.py`
- Have `prepare_VADAT.py` (or `submit_cmfgen_jobs.sh`) call
  `proxy_populations_v2.py` to generate _IN files
- Stage the generated _IN files in the CMFGEN run directory
- Set LTE_EST=F in VADAT when proxy populations are used

---

## 7. Test Results (VERIFIED)

```
57 passed in 3.08s
```

Test categories (from `scripts/test_proxy_populations_v2.py`):
- Parser tests (SN_HYDRO_DATA, F_OSCDAT, MODEL_SPEC, VADAT)
- Physics tests (1/2/3-stage Saha, closure, charge, LTE limit, grey)
- CMFGEN convention tests (DI=next stage, DC linear, level count, depth order, units)
- Round-trip tests (7 species, write+read)
- Continuation tests (read_in_file, regrid, full mode)
- Edge case tests (reversed arrays, zero kappa, missing abundances)
- Regression tests (real model: finite populations, abundance diff, Ne recomputed)
- GAMMAS_IN ED tests

---

## 8. Validation Issues (VERIFIED — expected)

The proxy reports 7 validation issues, all omitted-stage sink-risk warnings:
- SkI, Sk2, CaI, Ca2: highest modelled stage carries >50% (max 1.000)
- FeI, Fe2, FeIII: highest modelled stage carries >50% (max 0.998)

These are expected: the reduced model atom does not include the next ion stage
for these species (e.g. Ca III, Fe IV, Si IV). This means the proxy's DI estimate
for the highest modelled stage is less reliable. Recommendation: add the next ion
stage to MODEL_SPEC, or accept larger first-iteration corrections.

---

## 9. Abundance Residuals (VERIFIED)

VADAT abundances (depth-constant, as CMFGEN uses) vs hydro abundances (depth-varying):

| Element | VADAT | hydro mean | max dex diff |
|---|---|---|---|
| HYD | 7.001e-01 | 7.000e-01 | 0.000 |
| HE | 2.810e-01 | 2.823e-01 | 0.021 |
| OXY | 3.180e-06 | 1.450e-05 | 9.502 |
| SIL | 3.551e-05 | 1.592e-04 | 10.550 |
| CAL | 7.260e-07 | 1.000e-15 | 8.861 |
| IRON | 1.931e-02 | 1.686e-02 | 2.036 |

Large dex differences (OXY, SIL, CAL) reflect composition gradients, not proxy
errors. The proxy uses VADAT abundances (matching CMFGEN's SOLVEBA/SCALE_POPS) and
reports the differences.

---

## 10. Unresolved Assumptions Status

| # | Assumption | Status |
|---|---|---|
| 1 | Species _IN T column vs T_IN for LTEPOP | RESOLVED (T_IN) |
| 2 | Species _IN ED column vs GAMMAS_IN | RESOLVED (GAMMAS_IN) |
| 3 | F_TO_S: NF vs NS DCs | RESOLVED (NF safe) |
| 4 | NOLD < N extrapolation | OPEN (low risk for LTE) |
| 5 | Demonstrated improvement vs old proxy | BLOCKED (needs CMFGEN run) |
| 6 | Occupation probability W=1 omission | OPEN (physical-accuracy gap) |
| 7 | HDKT literal vs computed | RESOLVED (3e-5 relative, negligible) |

---

## 11. Final Decision

### PASS criteria evaluation

| Criterion | Status | Evidence |
|---|---|---|
| Exact atomic data used by generator and CMFGEN match | VERIFIED | ATOMIC_DATA_MANIFEST.json (65 files, same symlinks) |
| Exact model atom matches | VERIFIED | NOLD = MODEL_SPEC NF for all 10 species |
| Exact input-grid structure and ordering match | VERIFIED | ND=400, outer-to-inner, INPUT_GRID_MANIFEST.json |
| `_IN` semantics verified from current CMFGEN source | VERIFIED | DI=next stage, linear DC, LTEPOP formula (source-traced) |
| Generated proxy files pass independent parsing and round-trip tests | VERIFIED | 10/10 files pass, round-trip exact (max log10 ratio=0.0) |
| CMFGEN reads the intended proxy files | VERIFIED | Startup test: "*_IN files will be used to start model" |
| Fresh-model convergence is measured | PARTIAL | 2-iteration startup test done; full convergence BLOCKED (compute time) |
| Newly created `*_OUT` files compared with demonstrably comparable reference | STRUCTURAL ONLY | Reference not numerically comparable (different ND, atom, epoch, version) |
| Repository remains clean and well structured | VERIFIED | REPOSITORY_HYGIENE_REPORT_PROXY.md |

### Limitations

1. **Full convergence not measured:** The startup test ran only 2 iterations
   (47 seconds). Full convergence (~15,000 iterations, ~4 days) is BLOCKED by
   compute time. The proxy init is stable (no crash/NaN) and produces non-zero
   luminosity from iteration 1, but final converged populations/spectra were not
   compared between proxy and built-in initializations.

2. **Reference comparison is structural only:** The reference model
   (M15-He-O-Ca, ND=150, 6 species, 2014 CMFGEN) is NOT numerically comparable
   to the current model (SN_12C, ND=400, 10 species, 2026 CMFGEN). Only the file
   format was validated as identical. Per AGENTS.md: "If the reference grid or
   atom differs, produce a structural comparison but do not call it a validation
   of the current model."

3. **Orchestrator integration incomplete:** `proxy_populations_v2.py` works
   standalone but is NOT wired through `run_model_grid.py → prepare_VADAT.py →
   CMFGEN`. `prepare_VADAT.py` does not accept `--proxy-populations` (would
   crash). Production runs require manual staging until this is fixed.

4. **Occupation probability W=1:** The proxy omits the level-dissolution factor
   W that CMFGEN's `OCCUPATION_PROB` computes. This is a physical-accuracy gap
   (not a correctness bug) — the round-trip is exact because both sides use W=1
   in the emulation, but CMFGEN's first-iteration LTE will differ by the W factor
   in high-density layers.

### Status

**PASS WITH LIMITATIONS — proxy populations are correctly formatted, atomic-data
consistent, and validated against CMFGEN startup. Reference comparison is
structural only (not numerically comparable). Full convergence evidence is
limited (2-iteration startup test). Orchestrator integration is incomplete
(manual staging required for production).**

The proxy LTE populations are **usable** as a CMFGEN initialization within the
documented assumptions. To upgrade to PASS:
1. Wire `proxy_populations_v2.py` through the orchestrator
2. Run full convergence comparison (proxy init vs built-in init)
3. Find or generate a numerically comparable reference (same model atom + grid)
