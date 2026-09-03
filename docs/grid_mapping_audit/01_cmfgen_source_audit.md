# Audit 1: Local Modified CMFGEN Non-Monotonic Velocity Solver

**Date:** 2026-08-20
**Auditor:** Devin (GLM-5.2 High)
**Authoritative tree:** `/viper/ptmp1/sbrennan/SN2024pba/cmfgen_source/` (working-tree changes vs upstream `pitt-phyast/cmfsource` HEAD `dd1660a`)
**Reference (stock) tree:** `/viper/ptmp1/sbrennan/SN2024pba/cmfgen_reference/cur_cmf/` (extract of `cur_cmf_18jun25.tar.gz`)
**Cross-reference doc:** `/viper/ptmp1/sbrennan/jwst_cycle6/CMFGEN_IIn_Adaptation.md`

---

## 0. Method

The local modifications are **uncommitted working-tree changes** in `cmfgen_source/` (22 modified files, `git diff --stat HEAD`). The upstream repo (`pitt-phyast/cmfsource`) is the stock CMFGEN baseline. The audit below is derived from:

1. `git diff HEAD` for every modified Fortran file (the exact local changes).
2. Reading the full source of the key radiative-transfer routines.
3. Cross-referencing `CMFGEN_IIn_Adaptation.md`, which documents the intended modifications.
4. `grep` for remaining monotonicity assumptions across `subs/`, `plane/`, `new_main/`, `disp/`.

**A discrepancy was found between the adaptation doc and the actual source (see §3, BLOCKING).**

---

## 1. Relevant Local Source Files and Subroutines

### 1a. Files with uncommitted working-tree modifications (the local non-monotonic patch)

| File | Status | Role |
|------|--------|------|
| `new_main/subs/rv_sn_model_sniin.f` | **MODIFIED** | VEL_LAW=12 grid builder: reads `input_hydro.dat`, builds tau-based R grid, interpolates V (non-monotonic preserved), computes SIGMA, caps SIGMA at -0.999 |
| `new_main/cmfgen_sub.f` | MODIFIED | Passes `N_OB_INS` to `RV_SN_MODEL_SNIIN` (was hardcoded) |
| `new_main/mod_subs/rd_control_variables.f` | MODIFIED | VEL_LAW=12 now sets `SN_MODEL=.TRUE.`, `SN_HYDRO_MODEL=.TRUE.` (was missing) |
| `new_main/mod_subs/set_abund_clump.f` | MODIFIED | VELTYPE=12 branch: reads `input_hydro.dat` outer-to-inner (no reversal), interpolates T and estimates ED |
| `plane/fg_j_cmf_v13.f` | **MODIFIED** | Formal solution: Dessart+2015 GAM guard; SIGMA_RAY clamped to -0.999 for dGAMdR and GAMH derivatives |
| `plane/fg_j_cmf_v12.f` | MODIFIED | Same Dessart+2015 GAM guard as v13 (older version) |
| `subs/sobew.f` | MODIFIED | Sobolev optical depth: Dessart+2015 guard replaces `ABS(1+SIGMA*MU^2)` |
| `subs/sobew_grad.f` | MODIFIED | Same guard |
| `subs/sobew_grad_v2.f` | MODIFIED | Same guard |
| `subs/sobew_grad_v3.f` | MODIFIED | Same guard |
| `subs/sobjbar_sim.f` | MODIFIED | Same guard |
| `disp/subs/sobjbar.f` | MODIFIED | Same guard |
| `new_main/subs/rd_sn_data.f` | MODIFIED | Population consistency check: suppresses false-positive error when species abundance < 1e-15 |
| `new_main/subs/solveba_v13.f` | MODIFIED | Diagnostic: prints species/level index J of max corrections (first 5 iterations) |
| `new_main/subs/generate_full_matrix_v3.f` | MODIFIED | Guards `C_mat` multiplication against zero POPS (prevents NaN) |
| `subs/ltepop_wld_v2.f` | MODIFIED | Diagnostic print for Fe2 LTE population at depth 176 |
| `subs/regrid_log_dc_v1.f` | MODIFIED | Diagnostic print for Fe2_IN regridding at depth 176 |
| `misc/rewrite_old_model_data.f` | MODIFIED | `RETURN` → `STOP` (prevents fall-through) |
| `Makefile_definitions` | MODIFIED | DEBUG build settings |
| `disp/Makefile`, `new_main/Makefile`, `obs/Makefile` | MODIFIED | Build config |

### 1b. Key routines NOT modified but central to non-monotonic operation

| File | Role | Why it matters |
|------|------|----------------|
| `plane/mom_jrel_v9.f` | `MOM_JREL_V9` — moment solver designed for non-monotonic velocity | Invoked via `USE_J_REL=T` path in `comp_j_blank.f` L699-707. This is the Fix 5 path from the adaptation doc. **Not modified locally** — it is upstream code already present. |
| `new_main/mod_subs/comp_j_blank.f` | J computation dispatcher | Has three branches: `MOM_J_CMF_V11` (default, L624), `MOM_JREL_V9` (`USE_J_REL`, L699), `USE_LAM_ES` (L717). **See §3 for BLOCKING issue.** |
| `subs/lin_interp.f` | `LIN_INTERP` | Requires R strictly monotonic (increasing or decreasing). Used by `rv_sn_model_sniin.f` to interpolate V onto the R grid. **Never checks V monotonicity.** |
| `subs/mon_interp.f` / `mon_interp_fast.f` | `MON_INTERP` | Requires R strictly monotonic. Used by `set_abund_clump.f` for density/T interpolation. |
| `subs/mon_interp_nonmon.f` | `MON_INTERP_NONMON` | A non-monotonic-capable interpolator exists in the codebase but is **not used** by the VEL_LAW=12 path. |

---

## 2. Stock Monotonic-Flow Assumptions Removed or Relaxed

### 2a. SIGMA < -1 handling (the core non-monotonic change)

**Stock behavior:** `GAM = GAMH / ABS(1 + SIGMA*MU^2)` in the six Sobolev routines, and `GAM = T1*(1 + SIGMA*MU^2)` in `fg_j_cmf_v13/v12.f`. When `SIGMA < -1` and `MU^2 ≈ 1`, the denominator `1 + SIGMA*MU^2` crosses zero → division by zero / negative optical depth → NaN in the formal solution.

**Local modification (all 6 Sobolev routines + fg_j_cmf_v12/v13):** Replaced with the Dessart+2015 guard:
```fortran
T1 = 1.0_LDP + Z(I)*Z(I)/R(I)/R(I)*SIGMA(I)
T1 = MAX(SQRT(T1*T1 + 0.01_LDP), -0.2_LDP*(1.0_LDP + SIGMA(I)), 0.1_LDP)
GAM(I) = GAMH(I) / T1
```
This is algebraically smooth for **any** finite SIGMA (no singularity at SIGMA=-1). The `SQRT(T1*T1+0.01)` term regularizes the zero crossing; the `-0.2*(1+SIGMA)` and `0.1` floors prevent pathological values.

**In `fg_j_cmf_v13.f` additionally:** The derivative `dGAMdR` and half-point `GAMH` use SIGMA_RAY **directly** (not through the guard), so a separate clamp `IF(SIGMA_RAY < -0.999) SIGMA_RAY=-0.999` is applied **only** to the derivative computations (L962-980), not to GAM itself. The stock code had `IF(SIGMA_RAY .LE. -1.0) SIGMA_RAY=-0.999` for GAM (L939, now removed) — the local version removes that crude clamp from GAM (the guard handles it) but adds clamps to the derivative paths.

### 2b. VEL_LAW=12 grid construction (`rv_sn_model_sniin.f`)

**Stock behavior (from diff):** The stock `rv_sn_model_sniin.f` had a **bug**: `READ(LU,*) R_HYDRO(I),V_HYDRO,XD,...` — `V_HYDRO` was a scalar, not `V_HYDRO(I)`, so only the last row's velocity was kept. It also reversed `R_HYDRO` and `V_HYDRO` arrays (inner-to-outer → outer-to-inner) after reading.

**Local modification:**
- Fixed the read: `V_HYDRO(I)` (array, not scalar) — L147.
- Removed the array reversal: `R_HYDRO` and `V_HYDRO` are now read in the order they appear in `input_hydro.dat` (outer-to-inner, decreasing R), which matches `TA` (the tau-sorted R grid). No reversal needed. — L220-226.
- V is interpolated onto the R grid via `LIN_INTERP` preserving the non-monotonic profile. Boundary zones (1..NBND_INS and ND-NBND_INS..ND) are linearly extrapolated from the first/last hydro points. — L229-260.
- `NBND_INS` is now passed from VADAT (`N_OB_INS_IN`) instead of hardcoded 3. — L48-49.
- The upper bound check changed from `NBND_INS .GT. 3` to `NBND_INS .GT. ND/4`. — L58.
- Added the SIGMA cap: `IF(SIGMA(I) < -0.999) SIGMA(I)=-0.999` after the SIGMA computation loop. — (Documented in adaptation doc Fix 1; present in the working tree.)

### 2c. VEL_LAW=12 density/temperature reading (`set_abund_clump.f`)

**Stock behavior:** `READ(LUIN,*) R_HYDRO(NHYDRO-I+1), T1, RHO_HYDRO(NHYDRO-I+1), T_HYDRO(NHYDRO-I+1)` — reversed the array on read (assumed input was inner-to-outer).

**Local modification:** `READ(LUIN,*) R_HYDRO(I), T1, RHO_HYDRO(I), T_HYDRO(I)` — no reversal (input is outer-to-inner). Added T interpolation from hydro grid and ED estimation as singly-ionized POP_ATOM. — L322, L336-352.

### 2d. VEL_LAW=12 control variables (`rd_control_variables.f`)

**Stock behavior:** VEL_LAW=12 did **not** set `SN_MODEL` or `SN_HYDRO_MODEL` (a bug — these flags trigger SN-specific defaults like `JGREY_WITH_V_TERMS`, `RMDOT=1e-20`, etc.).

**Local modification:** Added `SN_MODEL=.TRUE.` and `SN_HYDRO_MODEL=.TRUE.` for VEL_LAW=12. — L178-179.

---

## 3. Monotonicity Assumptions That REMAIN Elsewhere in CMFGEN

### 3a. BLOCKING: `comp_j_blank.f` Eddington pre-computation path is NOT guarded

**The adaptation doc (Fix 4) claims** that `comp_j_blank.f` was modified to wrap the `MOM_J_CMF_V11` call in the ACCURATE/Eddington pre-computation section (L371-376) in a `USE_LAM_ES` check.

**The actual source:** `comp_j_blank.f` is **NOT** in the working-tree modifications (`git status` clean for this file; `git log` last touched by upstream commit `3e80865`). The ACCURATE section at L371-376 calls `MOM_J_CMF_V11` **unconditionally**:

```fortran
CALL TUNE(IONE,'MOM_J_CMF_ACC')
CALL MOM_J_CMF_V11(TA,CHIEXT,ESECEXT,VEXT,SIGMAEXT,REXT,
1               RJEXT,RSQHNU,VDOP_VEC_EXT,DELV_FRAC_MOM,
1              FL,dLOG_NU,INNER_BND_METH,DBB,IC,IB_STAB_FACTOR,
1              N_TYPE,H_CHK_OPTION,METHOD,COHERENT_ES,OUT_BC_TYPE,
1              FIRST_FREQ,NEW_FREQ,NCEXT,NPEXT,NDEXT)
```

There is **no** `IF(USE_LAM_ES)` or `IF(USE_J_REL)` guard around this call. The normal-iteration path (L624-717) has the three-way branch (`MOM_J_CMF_V11` / `MOM_JREL_V9` / `USE_LAM_ES`), but the Eddington pre-computation path does not.

**Consequence:** When `COMPUTE_EDDFAC=T` (first run, no EDDFACTOR file) and velocity is non-monotonic:
- If `N_TYPE=G_ONLY`: `MOM_J_CMF_V11` computes `G = N/H`. At a velocity reversal where `H→0`, G is undefined → J blows up. This is the original crash that Fix 4 was meant to prevent.
- If `N_TYPE=N_ON_J` (the Fix 5 setting): `MOM_J_CMF_V11` computes `N/J`, which is well-defined since `J>0`. **This path may be safe** even without the Fix 4 guard, because `N_ON_J` avoids the `N/H` division. However, `MOM_J_CMF_V11` (not `MOM_JREL_V9`) is still called — it is the stock moment solver, not the non-monotonic-aware one.

**Classification: BLOCKING for the audit, but may be mitigated by the Fix 5 `N_ON_J` setting.** This must be resolved before claiming the solver handles non-monotonic velocity in all code paths. The exact question: does `MOM_J_CMF_V11` with `N_TYPE=N_ON_J` produce correct, stable results for non-monotonic velocity, or does it require the `MOM_JREL_V9` solver? **This requires a test run (see §12).**

### 3b. R-grid monotonicity is REQUIRED (not relaxed)

`LIN_INTERP` (`subs/lin_interp.f` L4 comment: "R may be monotonically increasing, or decreasing") and `MON_INTERP` (`subs/mon_interp_fast.f` L75: "R values must be monotonic") both require R strictly monotonic. The VEL_LAW=12 path relies on these for V, density, and T interpolation. **Velocity is never checked for monotonicity by any interpolation routine.**

**Implication for the Python grid code:** The SN_HYDRO_DATA radius grid MUST be strictly monotonic (decreasing, outer-to-inner). Velocity may be non-monotonic. This is consistent with the Python code's current behavior.

### 3c. Frequency monotonicity (unchanged, irrelevant to grid mapping)

`set_cont_freq_v3.f`, `fg_j_cmf_v9.f`, etc. require frequency arrays monotonic. This is internal to CMFGEN and not affected by the grid mapping.

### 3d. `MON_INTERP_NONMON` exists but is unused by VEL_LAW=12

`subs/mon_interp_nonmon.f` provides a non-monotonic-capable interpolator, but the VEL_LAW=12 path uses `LIN_INTERP` (for V in `rv_sn_model_sniin.f`) and `MON_INTERP` (for density/T in `set_abund_clump.f`), both of which require monotonic R. This is fine because R is monotonic even when V is not.

---

## 4. Strictly Ordered Radial Grid Required? — YES

The modified solver requires R strictly monotonic (decreasing, outer-to-inner). This is enforced by:
- `LIN_INTERP` and `MON_INTERP` (error if non-monotonic R).
- `rv_sn_model_sniin.f` builds R from a log-tau grid (inherently monotonic).
- The Python `write_sn_hydro_data` orientation guard (L3555-3561) raises `ValueError` if R is non-monotonic.

**Velocity is NOT required to be monotonic.** This is the central design point of the local modifications.

---

## 5. Velocity May Be Zero or Change Sign?

### 5a. SIGMA computation requires V ≠ 0

`rv_sn_model_sniin.f` L266: `SIGMA(I) = R(I)*SIGMA(I)/V(I) - 1.0_LDP`. If `V(I)=0`, this is a division by zero. The SIGMA cap at -0.999 does NOT protect against V=0 (it only caps the result).

**The Python grid code** has a `--v-min-km-s` flag (default 0.0) that trims inner zones with V < this value. With the default (0.0), zones with V=0 are NOT trimmed. **If V can be exactly 0 in a STELLA snapshot, SIGMA will be NaN.** This is a HIGH-severity issue.

### 5b. Velocity sign change (V < 0)

The local solver does NOT check or forbid V < 0. The `input_hydro.dat` writer in `prepare_VADAT.py` (L1287) takes `np.abs(hydro.vel_cms)` when writing V — **this destroys the sign of V**. If the physical model has infalling material (V < 0), the `input_hydro.dat` will show |V|, and the SIGMA computation will treat it as outflowing.

**Classification: HIGH.** The `np.abs()` in `prepare_VADAT.py` L1287 is a sign-distorting operation that should be flagged. However, for the STELLA-to-CMFGEN path, STELLA ejecta are typically outflowing (V > 0), so this may not arise in practice. **Requires verification against actual STELLA snapshots.**

### 5c. Velocity reversal (dV/dr changes sign, V remains positive)

This is the intended use case (IIn CSM interaction). The local solver handles this via:
- The Dessart+2015 GAM guard (§2a).
- The SIGMA cap at -0.999 (§2b).
- The `MOM_JREL_V9` solver via `USE_J_REL=T` (normal iterations).
- `N_TYPE=N_ON_J` (avoids N/H division).

**This is supported in the normal-iteration path.** The Eddington pre-computation path (§3a) is the remaining question.

---

## 6. Extra Points Around Velocity Extrema Required?

**The local solver does NOT explicitly require extra points around velocity extrema.** There is no code in `rv_sn_model_sniin.f`, `fg_j_cmf_v13.f`, or the Sobolev routines that detects velocity extrema and demands extra grid points.

**However**, the SIGMA computation (L262-269) uses a 3-point stencil (I-1, I, I+1), and the Dessart+2015 guard's smoothness depends on SIGMA varying smoothly across the extremum. If the extremum is resolved by only 1-2 grid points, SIGMA will jump sharply, and the `dGAMdR` computation (which uses a 2-point stencil on SIGMA_RAY) will be inaccurate.

**Empirical requirement (from adaptation doc):** The X_2 model used `N_IB_INS=15`, `N_OB_INS=15` (inner/outer boundary inserted points) and a 150-point grid. The adaptation doc does not specify a minimum points-per-reversal count. **This is an unresolved empirical question (see §12).**

**Implication for the Python grid code:** The current code has no velocity-extremum-aware point insertion. The `adaptive_grid_refinement` function inserts points based on `dv/v` and `dsigma` thresholds, which indirectly resolves extrema, but there is no guaranteed minimum points-per-extremum. **This is a HIGH-severity gap for the refactor.**

---

## 7. Exact Definition and Use of Every Velocity-Gradient Quantity

### 7a. SIGMA (the primary velocity-gradient quantity)

**Definition** (from `rv_sn_model_sniin.f` L262-269):
```fortran
T1 = (V(I-1)-V(I))/(R(I-1)-R(I))      ! left-side dV/dR
T2 = (V(I)-V(I+1))/(R(I)-R(I+1))      ! right-side dV/dR
SIGMA(I) = (T1*(R(I)-R(I+1)) + T2*(R(I-1)-R(I))) / (R(I-1)-R(I))  ! distance-weighted dV/dR
SIGMA(I) = R(I)*SIGMA(I)/V(I) - 1.0   ! = dlnV/dlnR - 1
```

This is **SIGMA = dlnV/dlnR - 1** (dimensionless), computed via a distance-weighted average of the left and right finite-difference dV/dR, then non-dimensionalized as `R*dV/dR/V - 1`.

**Boundary:** `SIGMA(1) = R(1)*(V(1)-V(2))/(R(1)-R(2))/V(1) - 1` (one-sided), similarly `SIGMA(ND)`.

**Use in the solver:**
1. **Sobolev optical depth** (6 routines): `GAM = GAMH / GUARD(1 + SIGMA*MU^2)` where MU = Z/R (angle cosine). The guard is the Dessart+2015 formula. GAM is the angle-dependent Sobolev optical depth.
2. **Formal solution** (`fg_j_cmf_v13.f`): `GAM = T1 * GUARD(1 + SIGMA*MU^2)` where `T1 = 3.33564E-6 * V/R` (the CMF Doppler factor). `dGAMdR = GAM * SIGMA_clamped / R` (derivative for the formal solution Jacobian).
3. **Volume expansion factor** (`rd_sn_data.f` L255-257): `VOL_EXP_FAC = 1/(1 - T1*V/R*(SIGMA+1)) / (1 - T1*V/R)^2`. This uses SIGMA directly (no guard) for mass/volume computation. **If SIGMA < -1 and T1*V/R*(SIGMA+1) > 1, VOL_EXP_FAC diverges.** This is not guarded by the Dessart prescription.
4. **SIGMA cap**: `IF(SIGMA < -0.999) SIGMA = -0.999` in `rv_sn_model_sniin.f` (VEL_LAW=12 path). This cap is applied to the grid-point SIGMA before it is used by the solver. The ray-level `SIGMA_RAY` (interpolated onto ray points) is guarded separately in `fg_j_cmf_v13.f`.

**The Python `sigma_velocity` function (L732-773) uses the EXACT same formula.** This is confirmed: distance-weighted dV/dR, then `R*dV/dR/V - 1`. The Python code is correct here.

### 7b. dV/dR (not separately named)

Used internally in the SIGMA computation. Not exposed as a separate input or diagnostic.

### 7c. dlnV/dlnR (not separately computed)

Equivalent to `SIGMA + 1`. Not computed separately anywhere.

### 7d. GUARD(1 + SIGMA*MU^2) — the Dessart+2015 regularized denominator

```fortran
T1 = 1.0 + SIGMA*MU^2
T1 = MAX(SQRT(T1*T1 + 0.01), -0.2*(1.0+SIGMA), 0.1)
```
This is the quantity that replaces `ABS(1+SIGMA*MU^2)` in the stock code. It is smooth for all SIGMA, but its derivative with respect to SIGMA is not the same as the stock formula's. **The `dGAMdR` computation in `fg_j_cmf_v13.f` uses SIGMA_RAY directly (not through the guard), which is why the -0.999 clamp is applied there.**

---

## 8. Characteristic Crossings / Resonance Handling

**The local modifications do NOT add any characteristic-crossing or resonance-surface detection.** The Sobolev method inherently handles line resonances via the `GAM` (Sobolev optical depth) computation, which is local (no long-characteristic integration). The Dessart+2015 guard ensures `GAM` remains finite for SIGMA < -1.

**There is no explicit characteristic-tracing code** in the modified routines. The CMF formal solution (`fg_j_cmf_v13.f`) uses short characteristics along rays, with the Doppler factor `T1 = 3.33564E-6 * V/R` and the `GAM` term encoding the velocity-field geometry. Non-monotonic velocity changes the mapping between frequency and spatial location (a given frequency resonates at multiple radii), but the local solver does NOT add multi-resonance handling — it relies on the Sobolev approximation (local resonance) and the formal solution's ray-by-ray integration.

**Implication:** The Python grid code does not need to add characteristic-crossing anchors because the solver does not use them. The relevant quantity is SIGMA smoothness, not characteristic structure.

---

## 9. Restrictions on Cell Sizes, Velocity Jumps, or SIGMA Jumps

### 9a. No explicit cell-size restrictions from the solver

The solver does not enforce a minimum or maximum cell size. The `LIN_INTERP` and `MON_INTERP` routines only require monotonic R. However, the SIGMA computation (3-point stencil) becomes inaccurate if cells are too large relative to the velocity gradient scale.

### 9b. No explicit velocity-jump restrictions

The solver does not check for velocity jumps between adjacent cells. The SIGMA cap at -0.999 is the only restriction, and it acts on the computed SIGMA, not on the velocity jump directly.

### 9c. SIGMA jump restrictions (implicit)

The `dGAMdR` computation in `fg_j_cmf_v13.f` uses a 2-point stencil on `SIGMA_RAY`:
```fortran
dGAMdR(I,LS) = dGAMdR(I,LS) + T1*(MU^2)*(SIGMA_RAY(J)-SIGMA_RAY(K))/(R_RAY(J)-R_RAY(K))
```
If SIGMA jumps sharply between adjacent cells, this derivative is inaccurate. The Python code's `--max-sigma-jump` parameter (default 1.0) addresses this by inserting points where `|d(sigma)| > max_sigma_jump`. **This is a reasonable heuristic but is not derived from the solver's actual numerical requirements.**

### 9d. The -1.2 < SIGMA < -1.0 → -0.999 clamp (from AGENTS.md)

The AGENTS.md documents that this clamp exists in three places (`rv_sn_model_sniin.f`, `adjust_r_grid_v4.f`, `do_ng_accel_band_v4.f`) and that removing it causes catastrophic divergence. **Confirmed present in `rv_sn_model_sniin.f`** (the SIGMA cap). The other two files were not checked in this audit (they are STELLA files, not CMFGEN). **The Python grid code must NOT produce SIGMA in the range -1.2 to -1.0 without expecting it to be clamped to -0.999.**

---

## 10. Maximum Number of Turning Points / Velocity Reversals

**The local solver has NO explicit limit on the number of velocity reversals.** There is no counter, no array-size limit, and no check on the number of sign changes of dV/dr.

**However**, each reversal creates a region where SIGMA < -1 (capped to -0.999), and the Dessart+2015 guard regularizes the GAM computation. The more reversals, the more regions where the guard is active (rather than the physical formula), which may degrade accuracy. **There is no empirical test of the solver's behavior with many reversals (see §12).**

**Implication for the Python grid code:** No hard limit needs to be enforced, but the grid should resolve each reversal with enough points that SIGMA varies smoothly (see §6).

---

## 11. Required Input Flags and Their Defaults

### 11a. VADAT flags for non-monotonic velocity (from adaptation doc §"Key VADAT Settings")

| Flag | Value | Purpose | Set by |
|------|-------|---------|--------|
| `VEL_LAW` | 12 | Use `RV_SN_MODEL_SNIIN` (tau-based grid, non-monotonic V) | `prepare_VADAT.py` from `vel_law` YAML key |
| `USE_LAM_ES` | T or F | T: J=B(T), bypass moment solver. F: use moment solver | Fix 5 uses F (with `USE_J_REL=T`) |
| `USE_J_REL` | T | Use `MOM_JREL_V9` (non-monotonic moment solver) | `prepare_VADAT.py` for `static_snapshot` model type |
| `N_TYPE` | `N_ON_J` | Compute N/J (always defined). NOT `G_ONLY` (N/H undefined when H=0) | `prepare_VADAT.py` |
| `INCL_REL` | T | Include velocity terms in radiation transport | `prepare_VADAT.py` |
| `DO_DDT` | F | No time evolution (static snapshot) | YAML `model_type: static_snapshot` |
| `CHK_J` | `ABS_VAL` | Check J by absolute value (not relative) | `prepare_VADAT.py` |
| `METHOD` | `LOGMON` | Monotonic interpolation in log space | `prepare_VADAT.py` |
| `FIX_NE` | T (init/stable), F (final) | Ne fixed to hydro for bootstrap | `prepare_VADAT.py` |
| `FIX_T` | T | T fixed to SN_HYDRO_DATA | `prepare_VADAT.py` |
| `MAX_LAM` | 1.0D2 | 100x max population change per iteration | `prepare_VADAT.py` |
| `IBEG_NG` | 20 | Ng acceleration starts after 20 Lambda iterations | `prepare_VADAT.py` |
| `N_IB_INS` | 15 | Inner boundary inserted points | VADAT / YAML |
| `N_OB_INS` | 15 | Outer boundary inserted points | VADAT / YAML (now passed to `RV_SN_MODEL_SNIIN`) |
| `INC_RAD_DECAYS` | T | Include Ni56/Co56 decay heating | YAML |

### 11b. Default behavior when flags are NOT set

If `VEL_LAW != 12` (e.g., VEL_LAW=11 or default), the stock `set_rv_hydro_model_v3.f` is used, which has its own SIGMA cap at -0.999 but does NOT use the tau-based grid construction. **The non-monotonic modifications are ONLY active for VEL_LAW=12.**

If `USE_J_REL=F` and `USE_LAM_ES=F`, the stock `MOM_J_CMF_V11` is called with `N_TYPE` from VADAT. If `N_TYPE=G_ONLY` and velocity is non-monotonic, this crashes (N/H undefined). **The safe path requires `USE_J_REL=T` or `USE_LAM_ES=T` with `N_TYPE=N_ON_J`.**

---

## 12. Required Test Problems and Expected Convergence Behaviour

### 12a. Tests that DO NOT exist but are REQUIRED

**No test harness exists** in the local CMFGEN source or the Python grid code for verifying non-monotonic velocity handling. The adaptation doc describes empirical debugging on the X_2 model (5-45 days post-explosion) but not controlled synthetic tests.

**Required tests (BLOCKING for production trust):**

1. **Synthetic monotonic velocity model**: Homologous expansion (V ∝ R), no reversals. Must converge with all solver settings. Verifies the modifications don't break the stock case.

2. **Synthetic single-reversal model**: V increases to a maximum, then decreases (one dV/dr sign change). Must converge with `USE_J_REL=T`, `N_TYPE=N_ON_J`. Verifies the Dessart+2015 guard and `MOM_JREL_V9` handle one reversal.

3. **Synthetic multi-reversal model**: Several alternating acceleration/deceleration regions. Must converge. Determines the empirical minimum points-per-reversal.

4. **Eddington pre-computation test (§3a)**: Run with `COMPUTE_EDDFAC=T` (delete EDDFACTOR file) and non-monotonic velocity. Verify whether `MOM_J_CMF_V11` with `N_TYPE=N_ON_J` survives the ACCURATE section, or whether Fix 4 must be applied to `comp_j_blank.f`.

5. **V=0 test**: Verify SIGMA computation behavior when V approaches zero.

6. **Grid resolution sweep**: For a single-reversal model, vary the number of points around the reversal (2, 3, 5, 10, 20) and measure convergence rate and spectrum accuracy. Determines the minimum points-per-reversal.

### 12b. Expected convergence behaviour (from adaptation doc)

The X_2 model with the full fix chain (Fix 1-5) shows:
- INIT: corrections decrease ~6 orders of magnitude per iteration (1e112% → 1e72% over 23 iterations). Does NOT fully converge in 25 iterations but writes SCRTEMP and passes to STABLE.
- STABLE: continues convergence from INIT output.
- FINAL: full iterations with Ng acceleration.

**This is slow convergence** — the corrections start extremely high (1e112%) because of the proxy _IN initialization. The non-monotonic velocity itself is not the primary convergence bottleneck (the SIGMA cap and Dessart guard handle it); the bottleneck is the LTE initialization in cold outer layers.

---

## 13. Summary Classification

| # | Issue | Severity | Status |
|---|-------|----------|--------|
| 1 | `comp_j_blank.f` Eddington pre-computation path calls `MOM_J_CMF_V11` unconditionally (Fix 4 NOT applied) | **BLOCKING** | Must verify whether `N_ON_J` mitigates, or apply Fix 4 |
| 2 | V=0 causes SIGMA = NaN (division by zero in `rv_sn_model_sniin.f` L266) | **HIGH** | Python code must guarantee V > 0 (or handle V=0 zones) |
| 3 | `prepare_VADAT.py` L1287 takes `np.abs(vel_cms)` when writing `input_hydro.dat` — destroys V sign | **HIGH** | Must verify STELLA never produces V < 0, or remove abs() |
| 4 | No minimum-points-per-velocity-reversal requirement derived from solver tests | **HIGH** | Requires empirical test (§12a test 6) |
| 5 | `VOL_EXP_FAC` (`rd_sn_data.f` L255) uses SIGMA directly (no Dessart guard) — may diverge for SIGMA < -1 | **MEDIUM** | The -0.999 cap mitigates, but VOL_EXP_FAC accuracy at SIGMA≈-1 is questionable |
| 6 | No test harness for non-monotonic solver verification | **HIGH** | Must be built before production trust |
| 7 | No explicit characteristic/resonance handling in the modified solver | **LOW** | Sobolev approximation is local; not a grid-mapping concern |
| 8 | No limit on number of velocity reversals | **LOW** | No solver-imposed limit; accuracy degrades with more reversals |
| 9 | `MON_INTERP_NONMON` exists but is unused by VEL_LAW=12 | **LOW** | R is monotonic, so this is fine |
| 10 | SIGMA formula in Python matches CMFGEN exactly | **OK** | No action needed |
| 11 | R-grid strict monotonicity required (confirmed) | **OK** | Python code already enforces this |
| 12 | Velocity non-monotonicity is preserved (not forced monotonic) | **OK** | Python code already preserves it |

---

## 14. Key Conclusions for the Python Grid-Mapping Refactor

1. **The SIGMA formula is correct** in the Python code. Do not change it.

2. **R must be strictly monotonic** (decreasing, outer-to-inner). The Python code already enforces this. Do not relax it.

3. **Velocity must NOT be forced monotonic.** The Python code already preserves non-monotonic V. The solver handles SIGMA < -1 via the Dessart+2015 guard + -0.999 cap.

4. **V must be > 0** at all grid points (SIGMA divides by V). The Python code should guarantee this (trim V=0 zones, or floor V at a small positive value).

5. **The -0.999 SIGMA cap is the solver's safety net.** The Python code should NOT try to prevent SIGMA < -1 (that would require monotonicizing V). Instead, it should ensure SIGMA varies smoothly (no sharp jumps) so the cap is applied gradually, not as a step function.

6. **No characteristic/resonance anchors are needed.** The solver does not use them. Focus on SIGMA smoothness and velocity-extremum resolution instead.

7. **The Eddington pre-computation path (§3a) is a BLOCKING unknown.** Before production, verify whether `MOM_J_CMF_V11` with `N_ON_J` survives non-monotonic velocity in the ACCURATE section, or apply Fix 4 to `comp_j_blank.f`.

8. **Conservative remapping is NOT required by the solver** — the solver reads R, V, SIGMA, T, RHO, Ne, and mass fractions as point values and interpolates them itself. However, conservative remapping IS required for physical accuracy (mass, momentum, energy conservation), which is a separate concern from solver compatibility.
