# Audit 2: Python Grid-Mapping Code (`convert_stella_swd_to_cmfgen.py`)

**Date:** 2026-08-20
**Auditor:** Devin (GLM-5.2 High), with exhaustive line-by-line subagent scan
**File:** `/viper/ptmp1/sbrennan/SN2024pba/scripts/convert_stella_swd_to_cmfgen.py` (6587 lines)
**Role:** Converts STELLA `.swd` output → CMFGEN `SN_HYDRO_DATA`. Invoked by `run_model_grid.py` via two call sites with ~75 CLI arguments.
**Cross-reference:** Audit 1 (`01_cmfgen_source_audit.md`) for the solver requirements.

---

## 0. CLI Contract (MUST be preserved in any refactor)

`run_model_grid.py` calls this script with **~75 arguments** across two call sites:

- **Batch path** (`build_cmfgen_conversion_cmd`, L1920-2130): shell command string, all args as `--key value` pairs.
- **Direct path** (`run_stella_conversion`, L3825-3934): `subprocess.run` list, same args.

The full argument list (from the subagent audit §16) defines the **non-negotiable CLI contract**. Any refactor MUST accept all ~75 arguments with the same names, types, and defaults, or `run_model_grid.py` will break. See the subagent audit for the complete table.

**Key positional args:** `swd_file`, `output_file` (SN_HYDRO_DATA path).
**Key keyword args:** `--target-days`, `--kappa`, `--tau-max`, `--tau-min`, `--model-name`, `--base-grid-target`, `--base-grid-method`, `--tau-breakpoints`, `--tau-zone-fractions`, `--tau-spacing`, `--composition-gradient-threshold`, `--min-zones-per-shell`, `--trace-threshold`, `--composition-interface-cap`, `--max-thin-zones`, `--max-thick-tau`, `--min-points`, `--max-points`, `--min-nd`, `--max-radius-jump`, `--max-temp-jump`, `--max-density-jump`, `--max-ne-jump`, `--max-logtau-jump`, `--max-velocity-jump`, `--max-sigma-jump`, `--max-cell-ratio`, `--thin-grid-tau`, `--thin-grid-priority`, `--outer-boundary-extra-points`, `--outer-boundary-ratio`, `--min-dr-r`, `--max-dr-r`, `--sigma-priority-threshold`, `--sigma-priority`, `--dsigma-priority-threshold`, `--dsigma-priority`, `--ne-priority-threshold`, `--ne-priority`, `--t-floor`, `--t-cap`, `--trace-zero-threshold`, `--velocity-max-kms`, `--v-min-km-s`, `--sigma-conformance-headroom`, `--max-abs-sigma`, `--enable-sigma-conformance`/`--no-sigma-conformance`, `--no-sigma-optimize`, `--plot-dir`, `--rho-file`, `--abn-file`, `--input-profile`, `--xni-file`, `--enable-adaptive-grid`, `--use-tau-amr`, `--smooth-velocity`, `--sigma-threshold`, `--smooth-width`, `--smooth-sigma`/`--no-smooth-sigma`, `--max-dsigma`, `--sigma-smooth-passes`, `--sigma-smooth-window`, `--no-smooth-wiggles`/`--smooth-wiggles`, `--wiggle-window`, `--wiggle-max-change`, `--velocity-smooth-window-kms`, `--smooth-sigma-scatter`/`--no-smooth-sigma-scatter`, `--max-dsigma-scatter`, `--sigma-scatter-window`, `--sigma-scatter-dv-threshold`, `--smooth-sigma-global`/`--no-smooth-sigma-global`, `--sigma-global-window`, `--sigma-global-max-rel-change`, `--sigma-global-preserve-threshold`, `--composition-aware`/`--no-composition-aware`, `--composition-elements`, `--no-composition-grid-insertion`.

---

## 1. Grid-Reduction / Downsampling Functions

| Function | Lines | Reduces | Criterion | Can delete physical features? |
|----------|-------|---------|-----------|------------------------------|
| `downsample_grid` | 1351-1379 | R, v, T, rho, kappa, mass_fracs | Log-spaced R (target_points) | **YES** — pure log spacing misses shocks, velocity extrema, composition interfaces |
| `downsample_grid_feature_aware` | 1382-1482 | Same | High-gradient regions (gradient_threshold=0.05) + log-spaced fill | **PARTIAL** — preserves high-gradient edges but may miss narrow features below threshold |
| `downsample_grid_tau_based` | 1732-2117 | Same + tau | Tau spacing (linear/log/auto) with composition-aware interface insertion | **PARTIAL** — composition-aware mode preserves interfaces; pure tau spacing misses velocity features |
| `adaptive_grid_refinement` | 2215-2800 | Same + n_e | Priority-based AMR: dr/r, dT/T, drho/rho, dNe/Ne, dlogtau, dv/v, dsigma, abs(sigma) | **NO** (insert-only), but degenerate merge (L2741-2796) can remove points |
| `thin_degenerate_grid` | 1329-1348 | Same | dr/r < min_dr_r (default 0.005) | **YES** — removes near-degenerate points, could be physical shock zones |
| `fill_grid_to_minimum` | 2927-3054 | Same | Inserts at largest cells until min_points | **NO** (insert-only) |
| `enforce_grid_spacing` | 3370-3519 | Same + kappa | Phase 1: removes dr/r < min_dr_r; Phase 2: inserts dr/r > max_dr_r; Phase 3: smooths cell ratio | **YES** — Phase 1 removes below-min points |
| `enforce_sigma_conformance` | 3057-3367 | Same + kappa | dsigma, cell ratio, dr/r violations | **NO** (insert-only) |
| `regularize_grid_spacing` | 1029-1112 | R only | Median filter on log(dR) | **NO** — nudges positions ≤1%, doesn't delete |

**Critical:** The degenerate merge in `adaptive_grid_refinement` (L2741-2796) removes zones where `dr/r < 0.2% AND dv/v < 2% AND dT/T < 5% AND dNe/Ne < 10%`. This can delete physical shock-front zones that happen to have small dr/r.

---

## 2. Interpolation Operations

**All interpolation is independent per-quantity.** There is **NO conservative shell-overlap remapping** anywhere in the code.

### 2a. Major interpolation sites

| Line | Method | What | Coordinate | Conservation |
|------|--------|------|------------|--------------|
| L873, L877 | `np.interp` | log_v (fine grid + smooth) | log_r | Independent (V only) |
| L2053-2055 | `np.interp` | v, T, rho, kappa | r_increasing | **VIOLATION** — each quantity independently |
| L2076 | `np.interp` | mass fractions | r_increasing | **VIOLATION** — elemental mass not conserved |
| L2288-2298 | Linear/log | r, v, T, rho, n_e, tau, mf | f (insertion fraction) | **VIOLATION** — each independently |
| L2905-2912 | Linear/log | Same | frac | **VIOLATION** |
| L2967-2975 | Linear/log | Same | f (optimal_insert_fraction) | **VIOLATION** |
| L3124-3133 | Linear/log | Same + kappa | f | **VIOLATION** |
| L3417-3426 | Linear/log | Same + kappa | f (optimal) | **VIOLATION** |
| L5720-5735 | Linear/log | Same | 0.5 (midpoint) | **VIOLATION** |
| L5874 | `np.interp` | kappa (STELLA) | r | Independent (kappa, not conserved) |
| L5978 | `np.interp` | Ni56 | r_tau_trim | **VIOLATION** — Ni56 mass fraction interpolated, NOT conservatively remapped |
| L6215-6216 | `np.interp` | Ni56 | r_old | **VIOLATION** |
| L6242 | `np.interp` | v | stella_r_full | **VIOLATION** |
| L6254-6255 | `np.interp` | rho, T | stella_r_full | **VIOLATION** |
| L6296-6309 | `np.interp` | mass fractions | stella_r_full | **VIOLATION** |
| L6374-6375 | `np.interp` | Ni56 | r_old | **VIOLATION** |
| L6387-6391, L6535-6536 | `np.interp` | ne (STELLA) | ne_stella_r | Independent (Ne) |
| L349-354 | `scipy.interpolate.interp1d` | Ca (from input_profile) | profile_r | Independent (Ca) |

### 2b. Critical conservation violations

1. **Ni56 mass fraction is interpolated** (L5978, L6215-6216, L6374-6375), not conservatively remapped. Total Ni56 mass can change after grid remapping. **BLOCKING for radioactive-isotope conservation.**

2. **Density is independently log-interpolated** (L2060, L2293, L2907, etc.). Shell mass = ∫ρ dV is NOT conserved. **HIGH for mass conservation.**

3. **Velocity is independently linearly interpolated** (L2058, L2290, L2905, etc.). Momentum = ∫ρV dV and kinetic energy = ∫½ρV² dV are NOT conserved. **HIGH for momentum/energy conservation.**

4. **Mass fractions are independently interpolated** (L2076, L2297-2298, L2911-2912, etc.). Elemental mass = ∫X_element ρ dV is NOT conserved. **HIGH for composition conservation.**

5. **Temperature is independently log-interpolated** (L2059, L2292, L2906, etc.). Thermal energy = ∫ρcvT dV is NOT conserved. **MEDIUM (T is re-derived from hydro in CMFGEN anyway).**

---

## 3. Smoothing Operations

| Function | Lines | Smooths | Default window | Default ON/OFF | Can erase shocks/reversals? |
|----------|-------|---------|----------------|----------------|---------------------------|
| `smooth_velocity_discontinuities` | 776-837 | velocity | 5 pts | **OFF** (`--smooth-velocity`) | **YES** — 3-pass moving average at \|sigma\|>2.0 |
| `smooth_velocity_wiggles` | 840-895 | velocity | 3 pts | **ON** (default True) | **NO** — capped at 2% rel change |
| `smooth_velocity_window_kms` | 898-946 | velocity | 100 km/s | OFF | **NO** — capped at window_kms |
| `smooth_sigma_transitions` | 949-1026 | velocity (sigma gradient) | 0.5% of log-r | OFF (default False) | **PARTIAL** — can smooth transitions at reversals |
| `smooth_sigma_scatter` | 1115-1238 | velocity (sigma median) | 7 pts | **ON** (default True) | **NO** — only in near-constant-V (dv/v<1%), protects jumps |
| `smooth_sigma_global` | 1241-1326 | velocity (savgol on log-v) | 7 pts | **ON** (default True) | **NO** — protects \|sigma\|>2.0 regions |
| `regularize_grid_spacing` | 1029-1112 | radius (dR sequence) | 5 pts | Conditional | **NO** — nudges ≤1% |
| `enforce_sigma_conformance._local_smooth_v` | 3180-3198 | velocity | 5% rel change | Conditional | **NO** — capped at 5% |

**Critical:** `smooth_velocity_discontinuities` (L776-837) applies 3 passes of 5-point moving average at regions where |sigma| > 2.0. This can erase physical ejecta interfaces. It is OFF by default but can be enabled via `--smooth-velocity`. The code itself warns at L801.

**Two smoothers are ON by default** (`smooth_velocity_wiggles` and `smooth_sigma_scatter` and `smooth_sigma_global`). While they have relative-change caps that protect large-scale structure, they are still smoothing operations that the prompt's "no smoothing in production mode" policy would forbid. **The refactor must make these opt-in, not opt-out, for the non-monotonic solver profile.**

---

## 4. Point Insertion and Deletion

### 4a. Insertion sites

| Function | Lines | Where inserted |
|----------|-------|----------------|
| `adaptive_grid_refinement._insert_midpoint` | 2287-2298 | Worst-scoring edge (dr/r, dT/T, drho/rho, dNe/Ne, dlogtau, dv/v, dsigma, abs(sigma)) |
| `adaptive_grid_refinement` (sigma pass) | 2598-2653 | dsigma > max_sigma_jump edges |
| `adaptive_grid_refinement` (grid pass) | 2660-2687 | Worst cell-ratio > max_cell_ratio |
| `insert_outer_boundary_points` | 2803-2924 | Outer tau<1 region, largest cells |
| `fill_grid_to_minimum` | 2927-3054 | Largest cells (by dr/r) |
| `enforce_grid_spacing` (Phase 2) | 3454-3473 | dr/r > max_dr_r |
| `enforce_grid_spacing` (Phase 3) | 3475-3507 | Worst cell-ratio pair |
| `enforce_sigma_conformance` | 3057-3367 | dsigma/dr_r/cell-ratio violations |
| `downsample_grid_tau_based` (composition) | 1903-2040 | Composition interfaces (Fe-core/Si/O/O-He/etc.) |
| `refine_outer_boundary` (main) | 5680-5743 | Outer tau<1, evenly spaced |

### 4b. Deletion sites

| Function | Lines | Where deleted |
|----------|-------|---------------|
| `thin_degenerate_grid` | 1329-1348 | dr/r < min_dr_r (0.5%) |
| `adaptive_grid_refinement` (coalescing) | 2718-2739 | Sliver zones (<10x smaller than both neighbours, non-original) |
| `adaptive_grid_refinement` (degenerate merge) | 2741-2796 | dr/r<0.2% AND dv/v<2% AND dT/T<5% AND dNe/Ne<10% |
| `enforce_grid_spacing` (Phase 1) | 3428-3452 | dr/r < min_dr_r |
| Tau trimming (main) | 5540-5562 | Outer: tau < tau_min; Inner: tau > tau_max |
| Velocity trimming (main) | 5567-5584 | Inner: V < v_min_km_s |
| Ghost zone trimming (main) | 5614-5637 | Outer: V > velocity_max_kms |

**Critical:** The tau/velocity/ghost trimming (L5540-5637) can delete zones that contain physical features. The tau trimming uses different kappa sources for outer vs inner (analytical vs constant), which is inconsistent. **HIGH.**

---

## 5. Ordering Assumptions

**Primary ordering:** All arrays are **outer→inner (decreasing R)** throughout.

### 5a. Array reversals (`[::-1]`)

22 reversal sites identified (see subagent audit §5 for full list). Key ones:
- L114: `read_swd` — sorts outer→inner after reading
- L301, L359, L392, L431: `read_abn`/`read_rho`/`read_input_profile`/`read_xni` — reverse inner→outer to outer→inner
- L1590-1591: `hybrid_interp_mass_fractions` — reverses for `np.interp` (needs increasing)
- L3541-3554: `write_sn_hydro_data` orientation guard — flips if inner→outer detected
- L6242, L6254-6255, L6296-6309: main re-interpolation — reverses for `np.interp`

### 5b. Silent flip risk

**L3541-3554**: The orientation guard in `write_sn_hydro_data` flips arrays if `r[0] < r[-1]` (inner→outer detected). This is defensive but could silently flip if upstream code accidentally reverses ordering. **MEDIUM.**

### 5c. Ordering assertions

- L95, L271, L311: docstrings state "outer→inner (decreasing r)"
- L717, L751: `compute_kappa_tau` and `sigma_velocity` assume outer→inner
- L3555-3561: Final check raises `ValueError` if non-monotonic

---

## 6. Centre-vs-Boundary Assumptions

**R is treated as cell BOUNDARIES** (zone edges), not cell centres.

Evidence:
- L717-724: `compute_kappa_tau` uses trapezoidal rule between adjacent points
- L755-772: `sigma_velocity` uses differences between adjacent points
- L1069: `dR = r[:-1] - r[1:]` — cell sizes as point differences

**No explicit shell mass or volume integration.** The code does not compute shell masses or volumes; it relies on STELLA's native Lagrangian mass coordinates (log10(M_ext_Msun) in .swd L5). **This is a fundamental gap for conservative remapping.**

---

## 7. Optical-Depth Integrations

**Tau is ALWAYS reintegrated from final density+opacity** using trapezoidal rule (`compute_kappa_tau`, L711-725). It is never carried from the input STELLA file. This is correct.

**Opacity sources used (not distinguished clearly):**
- STELLA NLTE Rosseland kappa (floored at kappa_es)
- Kramers analytical approximation
- Constant electron-scattering kappa (0.34 cm²/g)

The code does NOT distinguish between Rosseland, flux-mean, and electron-scattering for tau calculation — it uses a single kappa value per call. **MEDIUM** — the prompt requires distinguishing these.

---

## 8. Electron-Density Recalculations

- L528-595: `calculate_electron_density` — **Saha ionization** with fixed n_e_ref=1e8 cm⁻³. Per-element ionization potentials hardcoded. Two-step ionization for He.
- L3636-3654: `write_sn_hydro_data` — uses STELLA n_e if available (interpolated), else Saha fallback.
- Recomputed at L5832, L5919, L6196, L6360, L6394 after grid changes.

**Saha assumptions:** Fixed reference n_e, LTE ionization equilibrium, no self-consistency iteration, ground-state statistical weights hardcoded. **MEDIUM** — the Saha approximation may be inaccurate in non-LTE regions.

---

## 9. Radioactive-Isotope Remapping

**Ni56 mass fraction is interpolated, not conservatively remapped** (L5978, L6215-6216, L6374-6375). Total Ni56 mass is NOT conserved after grid remapping. **BLOCKING for radioactive-isotope conservation.**

Decay formula (L3676-3678):
```python
Ni56_t = Ni56_initial * exp(-time_days / 8.8)
Co56_t = Ni56_initial * (111.3/(111.3-8.8)) * (exp(-time_days/111.3) - exp(-time_days/8.8))
Fe56_t = Ni56_initial - Ni56_t - Co56_t
```

---

## 10. Velocity-Gradient and Sigma Calculations

### 10a. Primary sigma (`sigma_velocity`, L732-773)

**Formula matches CMFGEN exactly** (confirmed against `rv_sn_model_sniin.f` L262-269):
```python
T1 = (V[i-1]-V[i])/(R[i-1]-R[i])      # left dV/dR
T2 = (V[i]-V[i+1])/(R[i]-R[i+1])      # right dV/dR
dVdR = (T1*h_right + T2*h_left)/(h_left+h_right)  # distance-weighted
sigma = R*dVdR/V - 1                   # = dlnV/dlnR - 1
```

### 10b. Ghost-zone sigma (L5632)

**DIFFERENT formula:** `sigma = log(v[i]/v[i+1]) / log(r[i]/r[i+1]) - 1` (logarithmic central difference). Used only for counting zones with sigma < -1. **LOW** — not used for the actual SIGMA written to SN_HYDRO_DATA.

### 10c. Local sigma helpers

`_sigma_at_index` (L2120-2137) and `adaptive_grid_refinement._local_sigma` (L2300-2329) use the same formula as `sigma_velocity`. **OK.**

---

## 11. Array Length Inconsistency Risks

| Stage | Lines | What can mismatch | Risk |
|-------|-------|-------------------|------|
| After `regularize_grid_spacing` | 1029-1112 | **R only** — v, T, rho NOT updated | **HIGH** |
| After velocity re-interpolation | L6242 | **v only** — sigma NOT recomputed | **HIGH** |
| After density/T re-interpolation | L6254-6255 | **rho, T only** — n_e, tau NOT recomputed | **HIGH** |
| After mass fraction re-interpolation | L6286-6309 | **mass_fractions only** — n_e NOT recomputed | **HIGH** |
| After Ni56 resampling | L5978, L6216, L6375 | **ni56 only** — length mismatch | **HIGH** |

The code handles these by recomputing dependent quantities at later stages (L6196, L6202, L6360-6362, L6378), but the temporary inconsistencies are a fragility risk.

---

## 12. Negative/Zero/Duplicate/Non-Finite/Non-Monotonic Radius

| Issue | Where | Check? |
|-------|-------|--------|
| Negative R | L221 (`10**log_r`) | **NO check** |
| Zero R | — | **NO check** — would cause div-by-zero in sigma |
| Duplicate R | L1897, L1980 | Removed in tau-based downsampling |
| Non-finite R | — | **NO check anywhere** |
| Non-monotonic R | L3555-3561 | Final check before write raises ValueError |

**Missing checks:** No NaN/inf guard at any stage. No negative/zero R guard. Monotonicity checked only at final write. **HIGH.**

---

## 13. Hard-Coded Thresholds

See subagent audit §14 for the full table of ~50 hard-coded thresholds with values and units. Key ones:

| Threshold | Value | Units | Context |
|-----------|-------|-------|---------|
| `MIN_DR_OVER_R` | 1e-4 | fractional | AMR minimum cell |
| `DEGENERATE_DR_R` | 2e-3 | fractional (0.2%) | Degenerate merge |
| `COALESCE_RATIO` | 10.0 | ratio | Sliver coalescing |
| `ABUND_FLOOR` | 1e-30 | mass fraction | Prevent SETPOP crash |
| `t_Ni` | 8.8 | days | Ni56 e-folding |
| `t_Co` | 111.3 | days | Co56 e-folding |
| `DEVIATION_THRESHOLD` | 0.05 | fractional | Boundary homologization |

All AMR jump thresholds (`max_radius_jump`, `max_temp_jump`, etc.) are **fractional** (e.g. 0.10 = 10% change). `max_sigma_jump=1.0` is dimensionless. `max_cell_ratio=1.5` is a ratio. These are clearly documented in the argparse help text.

---

## 14. Pipeline Order (major stages)

1. L4917-5202: Argument parsing
2. L5208-5263: Read .swd, select epoch
3. L5263-5272: Extract snapshot (r, v, T, rho, kappa, ne, composition)
4. L5306-5482: Load composition (.swd > .rho > .abn > solar)
5. L5433-5462: Load Ni56 from .xni
6. L5469-5482: Normalize mass fractions
7. L5499-5521: Compute tau (3 opacity sources)
8. L5533-5565: Tau trimming (outer: analytical, inner: constant kappa)
9. L5567-5584: Velocity trimming (V < v_min)
10. L5614-5637: Ghost zone trimming (V > velocity_max)
11. L5657-5668: Save full STELLA data for re-interpolation
12. L5680-5743: Outer boundary refinement
13. L5746-5905: Adaptive grid workflow (downsample → AMR → recompute kappa)
14. L5907-5915: Legacy grid thinning (if not AMR)
15. L5923-5945: Outer boundary points
16. L5947-5971: Fill to minimum points
17. L5973-5978: Resample Ni56
18. L5987-6028: Smoothing passes (velocity, wiggles, window, sigma transitions, scatter)
19. L6178-6209: Enforce grid spacing
20. L6210-6216: Resample Ni56 again
21. L6218-6228: Regularize grid spacing
22. L6230-6337: Re-interpolate v, rho, T, mass_fracs from full STELLA data
23. L6338-6347: Global sigma smoothing
24. L6349-6375: Final sigma conformance
25. L6377-6394: Compute sigma, n_e, atom density
26. L6410-6444: Calculate luminosity
27. L6446-6457: Write SN_HYDRO_DATA
28. L6471-6583: Diagnostic plots

---

## 15. Non-Monotonic Velocity Handling

**The code PRESERVES non-monotonic velocity.** It does NOT force monotonicity.

- L5607-5637: Ghost zone trimming explicitly keeps zones with sigma < -1
- L5609: Comment: "The non-monotonic velocity profile (sigma < -1) from the STELLA hydro is PHYSICAL and must be conserved"
- No monotonic spline, no isotonic regression, no sorting by velocity, no abs-value replacement
- Smoothers have relative-change caps that preserve large-scale structure (except `smooth_velocity_discontinuities` which is OFF by default)

**This is correct for the modified CMFGEN solver.** The solver handles SIGMA < -1 via the Dessart+2015 guard + -0.999 cap (see Audit 1).

---

## 16. Conservation Checks

**NO explicit conservation checks** for shell mass, momentum, kinetic energy, or optical depth.

- Mass fractions are renormalized to sum=1 at each depth (L3810-3815, etc.) — ensures X_sum=1 but NOT total elemental mass conservation.
- No check that total mass = ∫ρ dV is conserved.
- No check that total momentum = ∫ρV dV is conserved.
- No check that total KE = ∫½ρV² dV is conserved.
- No check that total optical depth is conserved.

**BLOCKING for the prompt's conservation requirements.**

---

## 17. Output Files

| File | Format | Purpose |
|------|--------|---------|
| `SN_HYDRO_DATA` | Text (CMFGEN keyword-delimited) | Main output |
| 7 PNG plots | PNG | Diagnostic plots (comparison, adaptive grid, mass fraction, dr_vs_r, zone density, structure overview, sigma vs index) |

**No JSON or audit files** are written. The prompt requires `grid_audit.json`, `grid_features.json`, `grid_validation.json`, `grid_validation.txt`. **BLOCKING for the prompt's reporting requirements.**

---

## 18. Feature Detection

| Detector | Criterion | Line |
|----------|-----------|------|
| Composition interfaces | \|dX\|/X > gradient_threshold AND max(X_left, X_right) > min_abundance | L1497-1550 |
| High-gradient edges | dr/r, dv/v, dT/T, drho/rho > gradient_threshold | L1413-1424 |
| Sigma spikes | \|sigma\| > sigma_threshold (2.0) | L795 |
| Sigma scatter | \|d(sigma)\| > max_dsigma AND \|dv\|/v < dv_threshold | L1162-1167 |

**Missing detectors:**
- **No explicit shock detection** (density/pressure/velocity jump with Rankine-Hugoniot criteria)
- **No explicit velocity extrema detection** (local min/max of V, sign change of dV/dr)
- **No explicit optical-depth landmark detection** (tau=1, photosphere, ionization fronts)

**BLOCKING for the prompt's feature-detection requirements.**

---

## 19. Summary of Critical Issues

| # | Issue | Severity |
|---|-------|----------|
| 1 | No conservative shell-overlap remapping — all interpolation is independent per-quantity | **BLOCKING** |
| 2 | Ni56 mass fraction interpolated, not conservatively remapped | **BLOCKING** |
| 3 | No conservation checks (mass, momentum, KE, optical depth) | **BLOCKING** |
| 4 | No JSON/text audit or validation reports written | **BLOCKING** |
| 5 | No explicit shock or velocity-extremum detection | **HIGH** |
| 6 | `smooth_velocity_discontinuities` can erase physical shocks (OFF by default but enabled by `--smooth-velocity`) | **HIGH** |
| 7 | Three smoothers ON by default (wiggles, sigma_scatter, sigma_global) — violates "no smoothing in production" | **HIGH** |
| 8 | Degenerate merge can delete physical shock-front zones | **HIGH** |
| 9 | No NaN/inf/negative/zero R checks at most stages | **HIGH** |
| 10 | Array length inconsistencies after regularize_grid_spacing and re-interpolation | **HIGH** |
| 11 | Tau trimming uses inconsistent kappa sources (analytical vs constant) | **MEDIUM** |
| 12 | Opacity sources (Rosseland/flux-mean/electron-scattering) not clearly distinguished | **MEDIUM** |
| 13 | Saha Ne approximation with fixed n_e_ref, no self-consistency | **MEDIUM** |
| 14 | Silent array flip risk in orientation guard (L3541-3554) | **MEDIUM** |
| 15 | Ghost-zone sigma formula (L5632) differs from primary sigma | **LOW** |
| 16 | SIGMA formula matches CMFGEN exactly | **OK** |
| 17 | Non-monotonic velocity is preserved (not forced monotonic) | **OK** |
| 18 | R-grid strict monotonicity enforced | **OK** |
| 19 | Tau is reintegrated from final density+opacity | **OK** |
