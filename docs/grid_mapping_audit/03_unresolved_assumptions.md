# Audit 3: Unresolved Assumptions and Required Source/Tests

**Date:** 2026-08-20
**Status:** Phase 1 complete. Refactor is GATED on resolving the items below.

Per the user's "Phase-by-phase, confirm gates" instruction and the prompt's "When evidence is insufficient, do not guess" directive, the following must be resolved BEFORE any production refactor. Each item lists the exact source routine or test required.

---

## BLOCKING (must resolve before refactor)

### B1. `comp_j_blank.f` Eddington pre-computation path — is Fix 4 needed?

**The discrepancy:** `CMFGEN_IIn_Adaptation.md` "Fix 4" claims `comp_j_blank.f` was modified to wrap the `MOM_J_CMF_V11` call in the ACCURATE/Eddington pre-computation section (L371-376) in a `USE_LAM_ES` check. The actual `cmfgen_source/new_main/mod_subs/comp_j_blank.f` is **NOT modified** (`git status` clean; last upstream commit `3e80865`). The ACCURATE section calls `MOM_J_CMF_V11` unconditionally.

**Why it matters:** If `COMPUTE_EDDFAC=T` (first run, no EDDFACTOR file) and velocity is non-monotonic, the Eddington pre-computation calls `MOM_J_CMF_V11` (the stock moment solver, not `MOM_JREL_V9`). With `N_TYPE=N_ON_J` (the Fix 5 setting), this computes N/J (well-defined) rather than N/H (undefined at H=0). **It may be safe, but this is unverified.**

**What is required:**
1. Run a non-monotonic velocity model with `COMPUTE_EDDFAC=T` (delete any existing EDDFACTOR file) and `N_TYPE=N_ON_J`, `USE_J_REL=T`. Check whether the ACCURATE section crashes or converges.
2. If it crashes: apply Fix 4 (wrap the L371 `MOM_J_CMF_V11` call in `IF(USE_LAM_ES) ... ELSE ... END IF`, or better, in `IF(USE_J_REL) ... ELSE IF(USE_LAM_ES) ... ELSE ... END IF` to use `MOM_JREL_V9` in the ACCURATE section too).
3. If it converges: document that `N_ON_J` mitigates the issue, and no source change is needed.

**Exact source routine:** `cmfgen_source/new_main/mod_subs/comp_j_blank.f` L368-408 (ACCURATE/Eddington pre-computation loop).

### B2. Does `MOM_J_CMF_V11` with `N_ON_J` handle non-monotonic velocity, or must `MOM_JREL_V9` be used?

**The question:** The normal-iteration path (L624-717) has three branches: `MOM_J_CMF_V11` (default), `MOM_JREL_V9` (`USE_J_REL=T`), `USE_LAM_ES`. The Fix 5 workflow uses `USE_J_REL=T` → `MOM_JREL_V9`. But the ACCURATE section (B1) always calls `MOM_J_CMF_V11`. If `MOM_J_CMF_V11` with `N_ON_J` is safe for non-monotonic V, then B1 is resolved. If not, the ACCURATE section must also use `MOM_JREL_V9`.

**What is required:** Read `plane/mom_j_cmf_v11.f` and `plane/mom_jrel_v9.f` to compare how they handle `N_TYPE=N_ON_J` with non-monotonic velocity. Specifically:
- Does `MOM_J_CMF_V11` compute N/J or N/H for `N_ON_J`?
- Does it have any internal monotonic-V assumption (e.g. characteristic tracing that assumes V monotonic)?
- Does `MOM_JREL_V9` differ in its handling of velocity reversals?

**Exact source routines:** `cmfgen_source/plane/mom_j_cmf_v11.f`, `cmfgen_source/plane/mom_jrel_v9.f`, `cmfgen_source/plane/var_jrel_v5.f`.

### B3. Can STELLA produce V ≤ 0 or V < 0 in any snapshot?

**The question:** `rv_sn_model_sniin.f` L266 computes `SIGMA = R*dV/dR/V - 1`. If V=0, this is division by zero → NaN. If V < 0, SIGMA is computed but the sign is physically meaningful (infall). The `prepare_VADAT.py` L1287 takes `np.abs(vel_cms)` when writing `input_hydro.dat`, destroying any V < 0 sign.

**What is required:**
1. Check all STELLA `.swd` snapshots for the SN_A, SN_B, X_1, X_2, X_3 models: does V ever reach 0 or go negative?
2. If V can be 0: the Python grid code must trim or floor V > 0 before writing SN_HYDRO_DATA.
3. If V can be < 0: decide whether to preserve the sign (remove `np.abs()` in `prepare_VADAT.py` L1287) or document that the solver only supports outflowing V.

**Exact test:** `python3 -c "for each .swd: read V column, print min(V), count(V<=0)"` across all model `.swd` files.

### B4. Minimum points per velocity reversal — empirical determination

**The question:** The local solver has no explicit minimum-points-per-reversal requirement, but the SIGMA computation (3-point stencil) and `dGAMdR` (2-point stencil on SIGMA_RAY) become inaccurate if a reversal is under-resolved. The adaptation doc does not specify a minimum.

**What is required:** Build the synthetic single-reversal test model (prompt §"Required synthetic tests"), run the modified solver with varying points-per-reversal (2, 3, 5, 10, 20), and measure:
- Convergence (does the solver converge?)
- Spectrum accuracy (compare to high-resolution reference)
- SIGMA smoothness across the reversal

**Exact test:** Synthetic model + CMFGEN run with `VEL_LAW=12`, `USE_J_REL=T`, `N_ON_J`. Vary the grid resolution at the reversal.

---

## HIGH (should resolve before production, but refactor can proceed with configurable parameters)

### H1. Which opacity source should be used for tau-based grid construction?

**The question:** The Python code uses three opacity sources (STELLA NLTE Rosseland, Kramers, constant e-scattering) at different stages. `prepare_VADAT.py` L1295 uses constant e-scattering kappa (0.34) for the `input_hydro.dat` TAU column, with a detailed comment explaining that NLTE kappa can be 60x lower in cold gas. The Python converter uses different kappa for inner vs outer tau trimming.

**What is required:** Decide and document which opacity source is authoritative for:
- Grid construction (where to place points)
- Inner boundary selection (tau_max)
- Outer boundary selection (tau_min)
- Optical-depth reporting

**Recommendation:** Use constant e-scattering kappa (0.34) for grid construction and boundary selection (consistent with `input_hydro.dat`), and STELLA NLTE Rosseland for diagnostic/reporting tau. Document this clearly.

### H2. Should the three default-ON smoothers be disabled for the non-monotonic profile?

**The question:** `smooth_velocity_wiggles`, `smooth_sigma_scatter`, and `smooth_sigma_global` are ON by default. They have relative-change caps (2%, 1%, 5%) that protect large-scale structure, but they are still smoothing operations. The prompt's "no smoothing in production mode" policy would forbid them.

**What is required:** Decide whether the non-monotonic solver profile disables all three by default, or keeps them with their caps. The caps are designed to preserve reversals (e.g. `smooth_sigma_global` protects |sigma|>2.0 regions), so they may be safe. But the prompt is explicit: "no velocity smoothing; no global sigma smoothing; no median filtering of physical profiles."

**Recommendation:** Disable all three by default in the non-monotonic profile. Add a `--regularization-mode` flag that re-enables them for users who want smoothing.

### H3. Conservative remapping — is it required for solver compatibility, or only for physical accuracy?

**The question:** The solver reads R, V, SIGMA, T, RHO, Ne, and mass fractions as point values and interpolates them internally. It does NOT check mass/momentum/energy conservation. So conservative remapping is NOT required for the solver to run. However, the prompt requires it for physical accuracy.

**What is required:** Decide whether the refactor implements conservative shell-overlap remapping (significant effort, changes the entire interpolation pipeline) or keeps independent interpolation with conservation checks/reporting (less effort, identifies violations but doesn't fix them).

**Recommendation:** Implement conservative remapping for density and composition (mass conservation), keep independent interpolation for T and Ne (not conserved), and document the velocity mapping choice (shell-centre vs mass-weighted vs momentum-conserving).

### H4. Shell mass / volume — does STELLA .swd provide enough information for conservative remapping?

**The question:** Conservative remapping requires shell boundaries and shell masses. The Python code treats R as cell boundaries but does NOT compute shell masses or volumes — it relies on STELLA's Lagrangian mass coordinates (log10(M_ext_Msun) in .swd L5). To do conservative remapping, we need:
- Shell boundaries (R[i], R[i+1]) — available
- Shell mass (from STELLA mass coordinate) — available in .swd
- Shell volume (4π/3 * (R[i]^3 - R[i+1]^3)) — computable
- Density from M_shell / V_shell — computable

**What is required:** Verify that the STELLA .swd mass coordinate is read and available for conservative remapping. The `read_swd` function (L95-) reads `log_r` and other columns — check if mass coordinate is included.

**Exact source:** `convert_stella_swd_to_cmfgen.py` L95-230 (`read_swd` function), and the STELLA .swd format documentation.

---

## MEDIUM (can be resolved during refactor)

### M1. `VOL_EXP_FAC` uses SIGMA directly (no Dessart guard)

`rd_sn_data.f` L255-257: `VOL_EXP_FAC = 1/(1 - T1*V/R*(SIGMA+1)) / (1 - T1*V/R)^2`. The -0.999 SIGMA cap prevents `SIGMA+1` from reaching -1, but `VOL_EXP_FAC` accuracy at SIGMA≈-0.999 is questionable. **The Python code should report VOL_EXP_FAC as a diagnostic.**

### M2. `prepare_VADAT.py` `np.abs(vel_cms)` — is this intentional?

L1287: `v_kms = np.abs(hydro.vel_cms) / 1e5`. If STELLA always produces V > 0 (outflowing ejecta), this is harmless. If V < 0 is possible (infall), this destroys the sign. **See B3.**

### M3. Saha Ne approximation

The `calculate_electron_density` function uses a fixed n_e_ref=1e8 cm⁻³ and single-pass Saha. This may be inaccurate in non-LTE regions. **The refactor should keep this as a fallback but prefer STELLA Ne when available.**

---

## LOW (informational, no action required for refactor)

### L1. `MON_INTERP_NONMON` exists but is unused

`subs/mon_interp_nonmon.f` provides non-monotonic-capable interpolation, but the VEL_LAW=12 path uses `LIN_INTERP` (for V) and `MON_INTERP` (for density/T). This is fine because R is monotonic.

### L2. Ghost-zone sigma formula differs from primary

L5632 uses `log(v[i]/v[i+1])/log(r[i]/r[i+1]) - 1` (logarithmic) vs the primary `sigma_velocity` (distance-weighted linear). Used only for counting, not for the SIGMA written to SN_HYDRO_DATA. No action needed.

### L3. No characteristic/resonance handling in the solver

The modified solver uses the Sobolev approximation (local resonance) and short-characteristics formal solution. No characteristic-crossing or resonance-surface detection is added. The Python grid code does not need to add characteristic anchors.

---

## Required Tests Before Production Trust

| # | Test | Purpose | Status |
|---|------|---------|--------|
| T1 | Synthetic monotonic V model | Verify modifications don't break stock case | **NOT BUILT** |
| T2 | Synthetic single-reversal model | Verify Dessart guard + MOM_JREL_V9 handle one reversal | **NOT BUILT** |
| T3 | Synthetic multi-reversal model | Verify multiple reversals converge | **NOT BUILT** |
| T4 | Eddington pre-computation test (B1) | Verify ACCURATE section with N_ON_J | **NOT BUILT** |
| T5 | V=0 / V<0 test (B3) | Verify SIGMA behavior at V→0 | **NOT BUILT** |
| T6 | Grid resolution sweep (B4) | Determine minimum points-per-reversal | **NOT BUILT** |
| T7 | Conservation validation | Verify mass/momentum/KE/tau conservation | **NOT BUILT** |
| T8 | Reduced-grid vs full-grid spectrum comparison | Verify grid reduction doesn't degrade spectra | **NOT BUILT** |

---

## Decision Points for the User (before refactor begins)

1. **B1/B2:** Should I apply Fix 4 to `comp_j_blank.f` as a precaution, or wait for the test to determine if it's needed?
2. **B3:** Should the refactor assume V > 0 always (trim/floor V=0 zones), or handle V < 0 (remove `np.abs()` in `prepare_VADAT.py`)?
3. **H2:** Should the non-monotonic profile disable all three default-ON smoothers, or keep them with their caps?
4. **H3:** Should the refactor implement full conservative shell-overlap remapping (significant effort), or independent interpolation with conservation checks/reporting?
5. **H4:** Is the STELLA .swd mass coordinate available for conservative remapping, or should I check the .swd format first?
6. **Scope:** Should the refactor be a new module set (`source_audit`, `stella_io`, `shell_model`, etc. as the prompt specifies) that replaces `convert_stella_swd_to_cmfgen.py`, or an in-place refactor of the existing 6587-line file?
