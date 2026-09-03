# Audit Report: `scripts/proxy_populations.py`

Generated **before** any code changes, as required by the task brief.
Every CMFGEN convention below was verified against the actual Fortran source
under `cmfgen_source/`; citations use `<ref_snippet .../>` with absolute paths.

## 0. Verified CMFGEN conventions (ground truth)

These are the non-negotiable facts the proxy must respect. Each was confirmed
by reading the source, not from comments.

### 0.1 `_IN` file reader = `REGRID_LOG_DC_V1`
<ref_snippet file="/viper/ptmp1/sbrennan/SN2024pba/cmfgen_source/subs/regrid_log_dc_v1.f" lines="6-7" />

Per-depth row read (clump-present branch):
<ref_snippet file="/viper/ptmp1/sbrennan/SN2024pba/cmfgen_source/subs/regrid_log_dc_v1.f" lines="162-167" />

**Column order:** `R, DI, ED, T, T1(unused), T2(unused), CLUMP_FAC` then a new
record with `NOLD` DC values. An 8th integer (depth index) on the data line is
ignored by list-directed read (next READ starts a new record). The proxy's
8-column data line is therefore format-compatible.

**Header:** `R(ND), RLUM, NOLD, ND` (line 139). `NOLD` = number of DC values
the file carries per depth.

### 0.2 DI = density of the NEXT HIGHER ion stage (CONFIRMED)
`REGRID_LOG_DC_V1` writes its interpolated `OLD_DI` into argument `DI`, which at
the call site is `ATM(ID)%DXzV_F`:
<ref_snippet file="/viper/ptmp1/sbrennan/SN2024pba/cmfgen_source/new_main/mod_subs/set_new_model_estimates.f" lines="166-170" />

`LTEPOP_WLD_V2` then uses `DXzV_F` as `DIC2` with `GION_C2` (next ion ground
weight) in the Saha expression `n_LTE ∝ DIC2 * ED * T^-1.5 / GION_C2 * exp(edge/kT)`:
<ref_snippet file="/viper/ptmp1/sbrennan/SN2024pba/cmfgen_source/new_main/mod_subs/set_new_model_estimates.f" lines="367-370" />

Decisive data-flow proof — the descending per-species loop sets each lower
stage's `DXzV_F` to the current stage's ground-state population:
<ref_snippet file="/viper/ptmp1/sbrennan/SN2024pba/cmfgen_source/new_main/mod_subs/set_new_model_estimates.f" lines="365-375" />

So `ATM(ID)%DXzV_F` holds the **next higher** stage's density. **The AGENTS.md
note (2026-08-02) is correct: `HI_IN`'s DI column = n(H II).** The proxy's
`compute_saha_ion_frac(..., ion_stage + 1, ...)` is correct. (One exploration
subagent inferred "current stage" from the variable name `OLD_DI` without
tracing the flow — that inference is wrong.)

### 0.3 DC are LINEAR when the format date is neither legacy date
<ref_snippet file="/viper/ptmp1/sbrennan/SN2024pba/cmfgen_source/subs/regrid_log_dc_v1.f" lines="127-135" />

<ref_snippet file="/viper/ptmp1/sbrennan/SN2024pba/cmfgen_source/subs/regrid_log_dc_v1.f" lines="176-181" />

- Format date `18-Nov-2010` → input is natural-log DC (LOG_TEN=1, no conversion).
- Format date `10-Dec-2010` → input is log10 DC (×log(10) to natural log).
- **Any other date (e.g. the proxy's `datetime.now()`, or `24-FEB-2004`) or no
  date → `TAKE_LOGS=.TRUE.` → input is LINEAR DC, `LOG()` applied.**
- `CHECK_DC` (the `|DC|<0.2 ⇒ treat as b-1` auto-detect) is ONLY active when
  there is **no** format date; the proxy always writes a date, so this never
  fires.

**The proxy writes a current date + linear DC ⇒ consistent.** Good.

### 0.4 Physical T and ED come from `T_IN` / `FIN_CAL_GRID`, NOT species `_IN`
<ref_snippet file="/viper/ptmp1/sbrennan/SN2024pba/cmfgen_source/new_main/mod_subs/set_new_model_estimates.f" lines="158-163" />

The species `_IN` file's T column is used **only** as an interpolation
coordinate for the `TR`/`RTX` DC-interpolation options (regrid_log_dc_v1.f
lines 560-591). It does **not** set the physical temperature structure.
**Safety rule #3:** writing a synthetic `T_ion` into the `_IN` T column cannot
alter the physical model T, but it *does* change DC interpolation in `TR`/`RTX`
mode. This must be opt-in and diagnosed.

### 0.5 `LTE_EST=T` completely bypasses the proxy `_IN` files (BLOCKING)
<ref_snippet file="/viper/ptmp1/sbrennan/SN2024pba/cmfgen_source/new_main/mod_subs/set_new_model_estimates.f" lines="261-270" />

When `LTE_EST=T`, `rd_control_variables.f` forces `DC_INTERP_METHOD='LTE'`:
<ref_snippet file="/viper/ptmp1/sbrennan/SN2024pba/cmfgen_source/new_main/mod_subs/rd_control_variables.f" lines="844-848" />

In that branch `REGRID_LOG_DC_V1` is **never called**: CMFGEN computes LTE DCs
via `SET_DC_LTE_V2` and sets `DXzV_F = 1.0E-200`. **The proxy `_IN` files are
not read at all.** Verified in real VADATs:

| model dir | `LTE_EST` | `DC_METH` | proxy `_IN` read? |
|---|---|---|---|
| `SN_12C/cmfgen/0.5days` | F | R | **yes** |
| `SN_12C_60days/cmfgen/30.0days` | **T** | R | **NO — ignored** |
| `SN_A/cmfgen/50.0days` | F | TR | yes |
| `X_v2/cmfgen/19.59days` | F | TR | yes |

### 0.6 `SET_DC_LTE_V2` stores log DC internally; DC=1 ⇒ LTE
<ref_snippet file="/viper/ptmp1/sbrennan/SN2024pba/cmfgen_source/subs/set_dc_lte_v2.f" lines="24-34" />

`XzVLTE` is the **natural log** of the departure coefficient; `0.0` means DC=1
(LTE). This is CMFGEN's *internal* representation after reading; the *file*
format is governed by §0.3.

### 0.7 `LTEPOP_WLD_V2` formula and constant
<ref_snippet file="/viper/ptmp1/sbrennan/SN2024pba/cmfgen_source/subs/ltepop_wld_v2.f" lines="66-76" />

Constant `2.07078E-22` (cm³ K^1.5), `GION_C2` = next-ion ground weight, returns
linear cm⁻³ plus `LOG_C2LTE`. The proxy's `compute_dc_nlte_reduced` uses the
same constant (`CMFGEN_C = 2.07078e-22`) — correct.

### 0.8 Abundance enforced by `SCALE_POPS`, not `SOLVEBA`
`SCALE_POPS` rescales populations so `SUM = POPHYD` (the VADAT abundance):
<ref_snippet file="/viper/ptmp1/sbrennan/SN2024pba/cmfgen_source/subs/scale_pops.f" lines="20-26" />

CMFGEN **re-derives** populations to match VADAT abundance each iteration. The
proxy need not match VADAT exactly, but a large mismatch produces large
first-iteration corrections (the reason the proxy swaps in VADAT abundances).

### 0.9 ED read linear from `GAMMAS_IN`, interpolated in log
<ref_snippet file="/viper/ptmp1/sbrennan/SN2024pba/cmfgen_source/subs/getelec_v2.f" lines="57-61" />

### 0.10 Depth ordering: outer→inner (depth 1 = outermost)
Confirmed by the AGENTS.md project convention and the `R(1) > R(ND)` grid. The
proxy uses outer→inner throughout — consistent.

### 0.11 `F_TO_S` and level count
`REGRID` reads `NOLD` DCs; `NZ = min(N, NOLD)` are directly interpolated
(`N` = model `N_F`). For levels `NZ+1..N`, only the `RTX` option extrapolates
via `T_EXCITE` (regrid_log_dc_v1.f 560-591); other options leave them to the
`F_TO_S` pass (609-619) which copies the lowest full level's DC across each
superlevel. In the real `MODEL` (SN_12C 0.5days) **N_S = N_F for every species**
(HI 12/12, HeI 10/10, Fe2 25/25, …), so `NOLD = min(N_S,N_F) = N_F` and the
concern is currently moot — but it is a latent bug for any future N_S < N_F
model with `DC_METH ≠ RTX`.

---

## 1. Issues affecting `DI`, `ED`, `T`, `DC`, level counts, output formatting

### BLOCKING

**B1. Proxy output is silently discarded when `LTE_EST=T`.**
The proxy runs its full grey-λ + Saha + DC pipeline and writes `_IN` files, but
if the target run uses `LTE_EST=T` (e.g. `SN_12C_60days/30.0days`), CMFGEN never
reads them (§0.5). The proxy gives no warning. A user can waste a full
optimization and then run CMFGEN from pure internal LTE. **Fix:** detect
`LTE_EST`/`DC_METH` in the target `VADAT`; if `LTE_EST=T`, either refuse to
write production `_IN` files or emit a prominent warning that they will be
ignored and that only `GAMMAS_IN`/`T_IN`/`FIN_CAL_GRID` matter.

**B2. No verification that the proxy's conventions match the *target* CMFGEN
binary/routine set.**
The proxy hardcodes `HDKT=4.7994145`, `SAHA_CONST=2.414e15`, `CMFGEN_C=2.07078e-22`.
`cmfgen.f` now computes `HDKT` from physical constants
(<ref_snippet file="/viper/ptmp1/sbrennan/SN2024pba/cmfgen_source/new_main/cmfgen.f" lines="117-122" />) and notes the old value `4.7994145`. The Saha constant in
`ltepop_wld_v2.f` is still the literal `2.07078E-22`, so the proxy matches
*that* routine — but this is unverified against the actually-built binary. The
brief requires: refuse production output if conventions cannot be verified
against the exact reader + a round-trip test.

### HIGH

**H1. `T_ion` written into the `_IN` T column in NLTE-aware mode (default when
`lte_floor=False`).**
`write_in_file(..., t_write=t_ion)` writes the hybrid radiation temperature into
every species `_IN` file. This does not change the physical T (§0.4) but, for
`DC_METH=TR`/`RTX` (used by `SN_A`, `X_v2`), it changes how DCs are interpolated
onto the new grid. This is a silent change to the interpolation model and
violates safety rule #3 unless explicitly opted in. **Fix:** default to writing
`T_local`; write `T_ion` only behind an explicit flag, with diagnostics showing
where `T_ion ≠ T_local`.

**H2. `smooth_depth_profile` (3-point median filter) is applied unconditionally
to every species' ion fraction.**
`generate_proxy_populations` always calls `smooth_depth_profile(ion_frac, ...)`
(lines ~1966). The brief (rule #6) says median filtering must be opt-in and
diagnostic because it can erase real ionization fronts. The code's own docstring
acknowledges ratio-clamping was removed for exactly this reason, but the median
filter that replaced it is still always-on. **Fix:** default off; add
`--smooth` opt-in; report which depths changed and by how much.

**H3. `min_ion_frac=1e-30` floor enters the ion population and hence DI.**
The floor is applied inside `compute_saha_ion_frac` (`f_ion = max(f_ion,
min_ion_frac)`) and again in the grey-λ loop (`stage_fracs = [max(f,
min_ion_frac)...]` then renormalized). The renormalization means the floor
*does* affect abundance closure and charge (rule #7). In cold outer zones where
the physical fraction is ~1e-100, flooring at 1e-30 inflates DI by 70 dex, which
directly inflates CMFGEN's LTE population estimate. This is intentional (per
AGENTS.md, to stay above FIXPOP) but is a numerical safeguard masquerading as a
physical value. **Fix:** separate physical fractions from floored output
fractions; report every depth where the floor dominates; exclude floored
contribution from charge/abundance accounting unless explicitly requested.

**H4. VADAT abundance swap silently overwrites hydro composition.**
Lines ~1711-1718 replace `hydro['mass_fracs']` with depth-constant VADAT
abundances, keeping hydro fracs only in `_hydro_mass_fracs`. This is physically
correct (CMFGEN enforces VADAT) but is a silent model change (rule #4): no
report of per-element, per-depth differences between hydro and VADAT. **Fix:**
emit a quantitative diff (max abs log difference per element, depths where they
differ by >N dex) to the diagnostics.

**H5. `compute_dc_dilute_lte` / `compute_dc` cap DC at 1e10–1e30.**
Caps are applied silently (lines ~789, ~849). A capped DC means the written
population is wrong by the cap factor. **Fix:** record capped depths and the
uncapped value; treat as a numerical safeguard, not a physical value.

**H6. `compute_dc_lte_cmfgen` can underflow to 0.0.**
`arg` is clipped to `[-700, 700]`; `exp(-700) ≈ 9.86e-305` (safe), but the
routine has no explicit nonzero floor, and any future change to the clip bound
would produce `0.0`, which `REGRID` would turn into `LOG(0)=-inf` (NaN). **Fix:**
explicit `max(val, 1e-300)` floor; add a round-trip finite-value check.

### MEDIUM

**M1. Level count `nlev = min(N_S, N_F)` is correct only because N_S=N_F today.**
If a future model has N_S < N_F and uses `DC_METH=R/ED/SPH_TAU`, levels
`N_S+1..N_F` are not extrapolated and rely on the `F_TO_S` copy pass with
uninitialized `DHEN` (§0.11). **Fix:** write `N_F` levels (filling
superlevels explicitly), or verify `N_S=N_F` and fail loudly otherwise.

**M2. Duplicate/contradictory DC formulas.**
Four DC routines exist: `compute_dc_lte_cmfgen`, `compute_dc`,
`compute_dc_dilute_lte`, `compute_dc_nlte_reduced`. Their docstrings disagree
on the sign of the Boltzmann exponent (the `compute_dc_dilute_lte` docstring
explicitly says the old `compute_dc` had "the wrong sign"). Three are reachable
from `generate_proxy_populations` via mode flags. This is exactly the
"duplicated or contradictory implementations" the brief flags. **Fix:**
collapse to one verified LTE formula + one named dilute-LTE strategy; delete the
others after tests pin behaviour.

**M3. `compute_saha_ion_frac` solves Saha in linear space with `exp` clip.**
The brief requires log-space Saha with `logsumexp`. The current code computes
`r_j = SAHA_CONST * T^1.5 / ne * exp(-chi/kt)` then a linear running product
`prod *= r[j]`. For long ion chains (Fe I→II→III→…) with extreme `r_j`, the
running product over/underflows even with the `±700` clip. **Fix:** rewrite in
log space (`logsumexp`), track log-fractions, exponentiate only at the end.

**M4. `parse_sn_hydro_data` fallback injects solar CARB=3e-3 when absent.**
Lines ~285-286 silently add a solar Carbon abundance if `CARB` is missing from
`SN_HYDRO_DATA`. This changes the composition model without VADAT cross-check
(rule #4). **Fix:** fail loudly or require explicit config; never invent
abundances.

**M5. `parse_f_oscdat` transition parser is fragile.**
The regex `(\d+)-\s*(\d+)\s+\d+\s*\|` and the "last 3 numbers before i-j"
heuristic misparse on non-standard layouts. Transitions are not used by the
population build, but a silent parse failure is undetected. **Fix:** validate or
drop transitions entirely (they are unused downstream).

**M6. `update_gammas_in_ne` patches `GAMMAS_IN` in place.**
It rewrites the ED section of the existing file. Rule #8 requires atomic writes
to a new output directory; the default `output_dir=model_dir` overwrites the
source. **Fix:** write to `output_dir`; only overwrite source with an explicit
flag.

**M7. `iterate_grey_lambda` mutates `hydro['ed']` and `hydro['kappa']` in place.**
Side effects on the input dict make the function non-reentrant and hide the
"input vs proxy" ED distinction the brief wants reported.

### LOW

**L1.** `HDKT`, `SAHA_CONST`, `M_H`, `SPEED_LIGHT` are hardcoded at module top
with no test citing their origin. `HDKT=4.7994145` matches the "old value"
comment in `cmfgen.f`; `SAHA_CONST=2.414e15` is `(2π m_e k/h²)^1.5 * 2` cgs —
undocumented in code.
**L2.** `detect_lte_floor_from_vadat` reads `USE_LAM_ES`/`FIX_T` and silently
chooses the proxy mode. This couples the proxy's physics mode to VADAT flags
without recording the decision in the diagnostics JSON.
**L3.** `write_in_file` writes `RLUM=0.0` (LSTAR=0 for SN). Correct for this
project, but undocumented; `DO_LAM_AUTO=T` with LSTAR=0 causes DTDR NaN (per
`prepare_IN_INTS.py`), so the zero must be paired with the right IN_ITS.
**L4.** Plotting code is interleaved with the build; no clean separation.
**L5.** `parse_cmfgen_convergence` regex-parses `OUTGEN`; brittle to log format
drift.

---

## 2. Depth-ordering assumptions

Every routine assumes **outer→inner** (depth 1 = outermost), matching
`SN_HYDRO_DATA` and CMFGEN. Verified consistent in: `parse_sn_hydro_data`,
`compute_tau` (integrates from depth 1 inward), `find_photosphere`,
`compute_w`/`compute_w_eff`, `write_in_file` (`r_nd = radius[nd-1]` = innermost),
`update_gammas_in_ne`. Plotting routines reverse to inner→outer for display
only. **No inner→outer assumption found in the build path.** Good.

## 3. Fallbacks that can silently change the physical model

| Location | Fallback | Risk |
|---|---|---|
| `parse_sn_hydro_data` ~285 | CARB=3e-3 solar default | M4 — invents abundance |
| `compute_saha_ion_frac` ~1165 | `n_stages==1` ⇒ even split | silent if no ionization data |
| `compute_tau` ~578 | kappa_es=0.2(1+X) if kappa=0 | changes optical depth / W / T_rad |
| `compute_t_rad` ~548 | `t_ph = t_local[nd-1]` if no τ=2/3 crossing | outer-edge T_rad guess |
| `get_atomic_info` ~215 | `gion_default=1.0` for unknown species | wrong Saha weight ratio |
| `build_saha_atomic_data` ~1085 | `g_ion=1.0` if next F_OSCDAT missing | wrong Saha weight ratio |
| `compute_dc_nlte_reduced` ~963 | `dc_val=1.0` if `n_lte_gs<1e-100` | hides a zero-LTE-population pathology |
| `iterate_grey_lambda` ~1445 | static W/T_ion if no elements with ionization data | silent NLTE→LTE switch |

## 4. NaN / overflow / underflow / zero-population / malformed-file risks

- `compute_dc*`: `exp` clipped to `±700` but no nonzero floor on the result
  (H6).
- `compute_saha_ion_frac`: `ne = max(ed[i], 1.0)` — if `ed` is NaN, `max`
  propagates NaN; no finite check.
- `compute_saha_ion_frac`: `kt_cm = 0.695 * T` with `T` already floored — safe,
  but the magic `0.695` (k_B/h*c in cm⁻¹/K) is undocumented.
- `get_excite_temp`: Newton step can oscillate; `t3==0` break leaves `t_exc`
  undefined-ish; no convergence-failure report.
- `parse_sn_hydro_data`: `int(re.search(r'\d+', lines[0]).group())` —
  `AttributeError` (crash) if header line has no digit; no graceful message.
- `parse_f_oscdat`: `float(parts[1])` with no guard — `IndexError`/`ValueError`
  on malformed level lines.
- `write_in_file`: if `ion_pop[i]` is NaN/inf, it is written verbatim → CMFGEN
  reads NaN → `LOG(NaN)=NaN` propagates.
- `update_gammas_in_ne`: if `nd` mismatches the file's actual ED-line count, it
  silently truncates or leaves stale lines.

## 5. Assumptions requiring confirmation from CMFGEN source / a real run

1. **Is the target run `LTE_EST=T` or `F`?** Determines whether the proxy is
   even read (B1). Must be read from the *actual* VADAT the run will use, not
   assumed.
2. **Which `DC_METH` will the run use?** Determines whether the `_IN` T column
   affects interpolation (H1) and whether N_S<N_F is safe (M1).
3. **Does the built CMFGEN binary use the literal `2.07078E-22` Saha constant
   or a constants-derived value?** Affects whether `CMFGEN_C` in
   `compute_dc_nlte_reduced` matches (B2). Resolve by inspecting the compiled
   `ltepop_wld_v2.f` / a round-trip test against the real binary.
4. **`GION_C2` source for each species:** is it the next stage's F_OSCDAT
   ground weight (proxy's `build_saha_atomic_data` assumption) or the MODEL-file
   `Gion` column? The proxy uses both in different branches — must pick one
   verified source.
5. **`W_C2(I,K)` (occupation probability / level dissolution):** `LTEPOP_WLD_V2`
   multiplies by `W_C2` (line 73). The proxy's LTE formula omits this factor.
   For SN ejecta with `DO_LEV_DISSOLUTION` enabled this could matter; needs a
   round-trip check.
6. **Round-trip tolerance:** no existing test reconstructs CMFGEN's LTE
   population from a written `_IN` file and compares to the target. Required by
   the brief (layer 8) and by acceptance criterion on first-iteration
   corrections.
7. **`FIXPOP` threshold:** the `min_ion_frac=1e-30` choice is justified against
   an undocumented FIXPOP value; confirm the exact threshold in the source.

---

## 6. Summary of priority counts

- **BLOCKING:** 2 (B1 silent discard when LTE_EST=T; B2 unverified constants
  vs. the real binary/reader)
- **HIGH:** 6 (H1 T_ion write, H2 always-on median filter, H3 floor enters
  charge/abundance, H4 silent abundance swap, H5 silent DC caps, H6 underflow
  to 0)
- **MEDIUM:** 7 (M1–M7)
- **LOW:** 5 (L1–L5)

The code's *file-format* conventions (DI=next stage, linear DC, outer→inner,
column order, units) are **correct and source-verified**. The serious problems
are behavioural: silent mode switches (B1, H1), always-on numerical filters
(H2), and floors/caps that leak into physical conservation (H3, H5). The
refactor will preserve the verified format layer and rebuild the strategy/
validation/reporting layers around it.

---

## 7. Repairs applied to V2 (2026-08-20 session)

### R1. DI for single-stage elements (HI, HeI, OI) — FIXED

**Root cause:** `build_saha_chains` only built Saha chains for consecutive
modelled stages. Single-stage elements (HI, HeI, OI — where the next ion
stage is not in MODEL_SPEC) had no chain, so `ion_density` returned `n_elem`
as DI. This is wrong: DI should be the next-stage density, which for hot
outer zones can be ~n_elem (H fully ionized) but for cooler zones or higher
ionization energies (He) is many orders of magnitude smaller.

**Fix:** `build_saha_chains` now also returns `chi_to_next`, a dict mapping
each modelled stage to its F_OSCDAT ionization energy and Saha gratio for
the transition to the next (possibly unmodelled) stage. `ion_density` uses
`chi_to_next` to estimate the next-stage density via a two-stage Saha
balance: `DI = n_elem × r/(1+r)`, where `r` is the Saha ratio. This is NOT
part of the closure chain — it is a physically-motivated DI estimate for
the _IN file, reported in diagnostics.

**Evidence:** At depth 1 (T=8279K, ne=1.717e9):
- HI: DI changed from 9.29e9 (n_elem) to 9.29e9 (n_elem × 0.998) — H is
  99.8% ionized, so DI ≈ n_elem is correct.
- HeI: DI changed from 9.36e9 (n_elem) to 2.14e6 — He is barely ionized
  (chi=198310 cm⁻¹), so DI is 5 orders of magnitude below n_elem.
- FeIII: DI changed from 9.36e9 (n_elem) to 0.72 — Fe IV is essentially
  absent at 8279K (chi=247220 cm⁻¹).

### R2. ED source: GAMMAS_IN instead of hydro ED — FIXED

**Root cause:** V2 used `hydro.ed` (from SN_HYDRO_DATA) for all physics
computations (Saha, LTEPOP, round-trip). CMFGEN's `GETELEC_V2` reads ED
from `GAMMAS_IN`, not from the species _IN ED column (audit §0.4/0.9).
The hydro ED and GAMMAS_IN ED differ by up to 295x at depth 1
(hydro: 9.449e9, GAMMAS_IN: 1.717e9).

**Fix:** `generate` now replaces `hydro.ed` with `inputs.gammas_in_ed`
when GAMMAS_IN ED is available. The ED source is reported in
`proxy_strategy.json` as `ed_source: "GAMMAS_IN"`.

### R3. Continuation mode (Mode C) — IMPLEMENTED

**Root cause:** `select_strategy` logged "continuation mode requested" but
`generate` never actually read or regridded the previous model's _IN files.

**Fix:** Added `read_in_file` (parses _IN format per `regrid_log_dc_v1.f`),
`regrid_log_dc` (log-space interpolation in R, matching REGRID_LOG_DC_V1),
and `build_continuation_species` (regrids DC/DI, computes LTE at new T/ED,
preserves regridded DCs without closure rescaling). `generate` now uses
these when `--strategy CONTINUATION --prev-model <dir>` is supplied.

### R4. Non-consecutive model stages — FIXED

**Root cause:** `build_saha_chains` iterated `while stage in stages` starting
from stage 1, breaking when a stage was missing. Non-consecutive stages
(e.g. Fe I and Fe III without Fe II) were never processed for `chi_to_next`.

**Fix:** Added a second pass in `build_saha_chains` that records
`chi_to_next` for ALL modelled stages, not just consecutive ones.

### R5. Round-trip metrics expanded — FIXED

**Root cause:** `roundtrip_error` reported only max/median/p95 relative
error, missing the log10 ratio, DC cap/floor counts, monotonicity, and
worst-location metrics required by the brief.

**Fix:** `roundtrip_error` now reports: max/median/p95 log10 ratio,
max/median/p95 relative error, DC min/max/p95/p99, n_dc_floored/capped
(counts and fractions), n_di_reversals, max_di_jump_dex, n_dc_jumps,
max_dc_jump_dex, worst_level/depth/log10.

### R6. Regression tests expanded — FIXED

Added 15 new tests (42 → 57 total):
- `TestExpandedPhysics`: 3-stage extreme ratios, extremely low next-ion
  density, DI estimation for single-stage, nonconsecutive model stages
- `TestEdgeCases`: reversed depth arrays, zero kappa, missing abundances,
  formatting invariance
- `TestContinuation`: read_in_file round-trip, regrid same/different grid,
  full continuation mode
- `TestSpeciesOrdering`: species order preserved
- `TestGammasInED`: GAMMAS_IN ED used and reported

### Test results

```
57 passed in 2.21s
```
