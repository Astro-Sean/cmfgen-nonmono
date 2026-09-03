# REFERENCE_COMPARABILITY_REPORT.md

**Date:** 2026-08-21
**Reference URL:** <https://ttt.astro.su.se/~maer0651/M15-He-O-Ca-evo-new-4-all/33/>
**Retrieval timestamp:** 2026-08-21 12:04 UTC

---

## 1. Retrieval Status (VERIFIED)

The reference URL was successfully accessed via HTTP. A directory listing was
retrieved (16990 bytes HTML). 34 key files were downloaded to
`data/reference/M15-He-O-Ca-evo-new-4-all/33/` with SHA-256 checksums recorded.

Downloaded files include:
- Population files: HIOUT, HeIOUT, He2OUT, OIOUT, CaIOUT, Ca2OUT (+ PRRR files)
- Model definition: MODEL_SPEC, MODEL, VADAT, VADAT_INIT, IN_ITS
- Grid: SN_HYDRO_DATA, NEW_SN_R_GRID, OLD_SN_R_GRID
- Diagnostics: CORRECTION_SUM, CORRECTION_LINK, OUTGEN, MOD_SUM, WARNINGS
- Other: GAMMAS_IN, GAMMAS, TRANS_INFO, batch.sh, SN_DATA_INPUT_CHK

Full checksum list in `data/reference/M15-He-O-Ca-evo-new-4-all/33/` (run `sha256sum *`).

---

## 2. Reference Model Description (VERIFIED)

| Property | Value |
|---|---|
| Model name | M15-He-O-Ca-evo-new-4-all, epoch 33 |
| Date | 2016-04-23 |
| CMFGEN version | Main program last changed 27-Mar-2014 |
| ND | 150 |
| NC | 15 |
| NP | 165 |
| NUM_BNDS | 3 |
| Species | HI[14/14], HeI[19/19], He2[10/10], OI[13/7], CaI[1/1], Ca2[17/12] |
| Total NF levels | 74 |
| Composition (by number) | H=0.907, He=9.197e-2, O=8.13e-3, Ca=4.85e-5 |
| V at Tau=278 | 494.8 km/s |
| V at Tau=20 | 3702.9 km/s |
| R* | 2.301e4 (10^10 cm) |
| RMAX/R* | 20.430 |
| Mean atomic mass | 1.2858 amu |
| DC content | Converged NLTE departure coefficients |

---

## 3. Current Model Description (VERIFIED)

| Property | Value |
|---|---|
| Model name | SN_12C/cmfgen/0.5days |
| Date | 2026-08-06 |
| CMFGEN version | Current build (2026, cmfgen_source) |
| ND | 400 |
| Species | HI[12], HeI[10], OI[7], SkI[1], Sk2[10], CaI[1], Ca2[12], FeI[1], Fe2[25], FeIII[1] |
| Total NF levels | 80 |
| Time | 0.5532 days post-explosion |
| LTE_EST | F (proxy _IN files are read) |
| DC_METH | R (R-based interpolation) |
| DC content | LTE (DC=1, proxy default) |

---

## 4. Comparability Assessment (VERIFIED — NOT COMPARABLE)

| Criterion | Proxy (SN_12C) | Reference (M15) | Comparable? |
|---|---|---|---|
| ND | 400 | 150 | **NO** |
| Species count | 10 | 6 | **NO** |
| Shared species | 5 (HI, HeI, OI, CaI, Ca2) | 5 | partial |
| HI levels (NF) | 12 | 14 | **NO** |
| HeI levels (NF) | 10 | 19 | **NO** |
| OI levels (NF) | 7 | 13 | **NO** |
| CaI levels (NF) | 1 | 1 | YES |
| Ca2 levels (NF) | 12 | 17 | **NO** |
| He2 | absent | present (10) | **NO** |
| Si (Sk) | present (11) | absent | **NO** |
| Fe | present (27) | absent | **NO** |
| Epoch | 0.55 d | evolution seq. epoch 33 | **NO** |
| CMFGEN version | 2026 build | 2014 | **NO** |
| DC content | LTE (DC=1) | converged NLTE | **NO** |
| File format | identical | identical | **YES** |

**Conclusion:** The reference and current models are **NOT numerically comparable**.
They differ in grid size, model atom, epoch, physical model, CMFGEN version, and
DC content. Only the file format is identical.

---

## 5. Structural Comparison (VERIFIED)

Both the proxy _IN files and the reference _OUT files were parsed by the same
independent parser (`tests/test_independent_in_parser.py`). Both parse
successfully with 0 errors, confirming identical file format:

| File | ND | NOLD | Format date | DC range | Parse errors |
|---|---|---|---|---|---|
| Proxy HI_IN | 400 | 12 | 21-Aug-2026 | [1.0, 1.0] | 0 |
| Ref HIOUT | 150 | 14 | 24-FEB-2004 | [0.496, 35.65] | 0 |
| Proxy CaI_IN | 400 | 1 | 21-Aug-2026 | [1.0, 1.0] | 0 |
| Ref CaIOUT | 150 | 1 | 24-FEB-2004 | [0.139, 2.138] | 0 |

All files: finite values, non-negative, outer-to-inner ordering, NOLD matches
MODEL_SPEC NF.

---

## 6. Shared Species DC Range Comparison

| Species | Proxy DC range | Reference DC range | Interpretation |
|---|---|---|---|
| HI | [1.0, 1.0] | [0.496, 35.65] | LTE vs converged NLTE |
| HeI | [1.0, 1.0] | [0.019, 45.41] | LTE vs converged NLTE |
| OI | [1.0, 1.0] | [0.347, 196.7] | LTE vs converged NLTE |
| CaI | [1.0, 1.0] | [0.139, 2.138] | LTE vs converged NLTE |
| Ca2 | [1.0, 1.0] | [0.006, 3.892] | LTE vs converged NLTE |

The proxy writes DC=1 (LTE); the reference contains converged NLTE departure
coefficients. These are fundamentally different physical states and cannot be
compared numerically.

---

## 7. Conclusion

Per AGENTS.md: "If the reference grid or atom differs, produce a structural
comparison but do not call it a validation of the current model."

The comparison is **STRUCTURAL ONLY**:
- File format is identical (verified by independent parsing)
- Column structure, depth ordering, and units are consistent
- NOLD matches MODEL_SPEC NF on both sides

The comparison is **NOT a numerical validation** of the proxy's physical
populations. A numerically comparable reference would require:
- Same model atom (same species, same level counts)
- Same or remappable grid (same ND or documented remapping)
- Same epoch or documented time difference
- Same CMFGEN version and atomic data

No such reference is currently available. The reference serves as a format
validation source only.

---

## 8. Recommendation

To enable numerical OUT-file validation, either:
1. Run the current model atom (SN_12C, 10 species) to convergence and compare
   the resulting _OUT files against a previous converged run of the same model.
2. Find or generate a reference with the same model atom and grid as the target.
3. Use the reference only for format validation (as done here).
