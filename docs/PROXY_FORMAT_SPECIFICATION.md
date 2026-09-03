# PROXY_FORMAT_SPECIFICATION.md

**Date:** 2026-08-21
**Scope:** Byte-level format specification for CMFGEN proxy population `_IN` files,
verified against `cmfgen_source/subs/regrid_log_dc_v1.f`.

---

## 1. File Format (VERIFIED from `regrid_log_dc_v1.f`)

Each species `_IN` file (e.g. `HI_IN`, `Ca2_IN`) has the following structure:

```
<blank line>
<format date>                        ! e.g. "21-Aug-2026" or "24-FEB-2004"
<blank line>
R(ND)  RLUM  NOLD  ND                ! header: 4 values, whitespace-separated
<blank line>
R  DI  ED  T  T1  T2  CLUMP_FAC  I   ! per-depth row (8 fields; I = depth index)
DC_1  DC_2  DC_3  DC_4  DC_5         ! DC block: NOLD values, 5 per line
DC_6  DC_7  ...                      ! (continuation lines as needed)
<blank line>
R  DI  ED  T  T1  T2  CLUMP_FAC  I   ! next depth
DC_1  DC_2  ...
...
```

### 1.1 Header

| Field | Type | Units | Description |
|---|---|---|---|
| R(ND) | float | 10^10 cm | Radius at innermost depth point (depth ND) |
| RLUM | float | Lsun | Luminosity (0.0 for SN models) |
| NOLD | int | — | Number of DC values per depth (= MODEL_SPEC NF) |
| ND | int | — | Number of depth points |

### 1.2 Per-Depth Row

| Column | Type | Units | Description | Used by CMFGEN? |
|---|---|---|---|---|
| R | float | 10^10 cm | Radius at this depth | Yes (interpolation coordinate) |
| DI | float | cm^-3 | Density of NEXT HIGHER ion stage | Yes (LTEPOP_WLD_V2) |
| ED | float | cm^-3 | Electron density | Yes (interpolation ref only; model ED from GAMMAS_IN) |
| T | float | 10^4 K | Temperature | Yes (interpolation ref only; model T from T_IN) |
| T1 | float | — | Unused (set to 0.0) | No |
| T2 | float | — | Unused (set to 0.0) | No |
| CLUMP_FAC | float | — | Clumping factor (1.0 = no clumping) | Yes |
| I | int | — | Depth index (1-based, outer-to-inner) | No (ignored by list-directed read) |

### 1.3 DC Block

- `NOLD` departure coefficient values per depth
- Written 5 per line, whitespace-separated
- Format: E17.7-ish (Fortran list-directed read accepts any reasonable float format)
- **DC are LINEAR** when format date is not one of the two legacy dates
  (`18-Nov-2010` = natural-log DC, `10-Dec-2010` = log10 DC)
- CMFGEN applies `LOG()` internally via `TAKE_LOGS=.TRUE.` for modern dates

### 1.4 Format Date

The proxy writes the current date (e.g. `21-Aug-2026`). Any date that is NOT
`18-Nov-2010` or `10-Dec-2010` triggers `TAKE_LOGS=.TRUE.`, meaning the file
contains LINEAR DC values. This is the proxy's default and is correct.

---

## 2. Depth Ordering (VERIFIED)

- **Outer-to-inner:** depth 1 = outermost (lowest density, highest velocity)
- Radius decreases with increasing depth index: R(1) > R(2) > ... > R(ND)
- This matches `SN_HYDRO_DATA` and CMFGEN's internal convention
- The proxy writes in this order throughout

---

## 3. DI Semantics (VERIFIED)

`DI` = density of the **next higher** ionization stage, NOT the current stage.

Evidence: `set_new_model_estimates.f:365-375` — the descending per-species loop
sets each lower stage's `DXzV_F` to the current stage's ground-state population.
`LTEPOP_WLD_V2` then uses `DXzV_F` as `DIC2` with `GION_C2` (next ion ground weight).

Example: `HI_IN`'s DI column = n(H II), not n(H I).

The proxy correctly computes DI via `compute_saha_ion_frac(..., ion_stage + 1, ...)`
and the `chi_to_next` Saha estimate for single-stage elements.

---

## 4. DC Semantics (VERIFIED)

- DC = departure coefficient = population / LTE_population
- DC = 1.0 means LTE (full local thermodynamic equilibrium)
- DC < 1.0 means sub-LTE (population below LTE)
- DC > 1.0 means super-LTE (population above LTE)
- In LTE mode (proxy default), all DC = 1.0
- In GREY mode, DC = W (dilution factor) for optically-thin regions
- In CONTINUATION mode, DC = regridded values from previous converged model

CMFGEN reconstructs populations as:
```
population(J) = EXP(ln(DC) + ln(LTE_pop)) = DC * LTE_pop
```

---

## 5. File Naming (VERIFIED)

Species _IN files are named `<SPECIES>_IN`:
- `HI_IN`, `HeI_IN`, `OI_IN`, `SkI_IN`, `Sk2_IN`, `CaI_IN`, `Ca2_IN`,
  `FeI_IN`, `Fe2_IN`, `FeIII_IN`

Non-species files (not population files):
- `GAMMAS_IN` — electron density input
- `T_IN` — temperature input
- `CFDAT__IN` — continuum frequency data

CMFGEN's `out2in.sh` renames `*OUT` → `*_IN` for warm restarts:
- `HIOUT` → `HI_IN`, `Ca2OUT` → `Ca2_IN`, etc.

---

## 6. Numerical Precision

- All values written in Fortran-compatible scientific notation (e.g. `1.2345678E+10`)
- Precision: ~8 significant digits (Python float repr)
- No blank lines within a depth record (would shift formatted reads)
- No trailing data after the last depth's DC block

---

## 7. Validation Checks (VERIFIED by independent parser)

Every generated _IN file must satisfy:
1. All values finite (no NaN, no Inf)
2. No negative DI, ED, or DC values
3. Radius monotonically decreasing (outer-to-inner)
4. ND consistent across all species files
5. R, ED, T arrays identical across all species files
6. NOLD = MODEL_SPEC NF for each species
7. CLUMP_FAC = 1.0
8. Header R(ND) matches innermost radius
9. No missing depth points
10. No duplicated depth indices

The independent parser (`tests/test_independent_in_parser.py`) verifies all of
these without importing the proxy generator.
