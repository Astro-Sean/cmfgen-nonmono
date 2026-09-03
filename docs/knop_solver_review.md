# Review: Non-Monotonic CMFGEN Formal Solver vs Knop et al. (2009)

## Executive Summary

The non-monotonic formal solver in `fg_j_cmf_nonmono_v1.f` **correctly implements the core Knop et al. (2009) formal solution** for radiative transfer in arbitrary velocity fields. The generalized opacity is kept positive, the signed wavelength coupling is retained in the source, and the sign-dependent upwind direction is correctly selected based on the ray-dependent coupling term `GAM`, not merely on `SIGMA`.

The model is **currently converging** with `LTE_EST=T` (great iteration 6, corrections decreasing from 1.36E+06 to 1.95E+03). The earlier blow-up was caused by outer-iteration instability (population solver feedback), not by a formal-solver error.

Three differences from the published Knop algorithm are identified and assessed below. None is a bug; one (the lack of inner GS iteration) is a performance/convergence-rate limitation that explains the slow convergence.

---

## 1. Papers Reviewed

| Paper | Citation | Key contribution |
|-------|----------|-----------------|
| Baron & Hauschildt 2004 | A&A 427, 987 (astro-ph/0408212) | Original algorithm: boundary-value problem in wavelength, direct matrix solver |
| Knop, Hauschildt & Baron 2009a | A&A 493, 1197 (arXiv:0902.1876) | New formal solution: removes 4a from opacity, keeps it as explicit source |
| Knop, Hauschildt & Baron 2009b | A&A 501, 813 (arXiv:0909.4264) | Large-scale: Gauss-Seidel iteration, domain decomposition, quasi-analytic solution |

## 2. Key Equations from the Papers

### 2.1 Transfer equation (Knop 2009a, Eq. 1)

```
dI_l/ds = eta_l - chi_l * I_l - 4*a_l * I_l - a_l * d(lambda*I)/dlambda
```

where `a_l = dlambda/ds` is the wavelength-coupling term (Eq. 2):
```
a_l = gamma * (gamma^2 * mu*(mu+beta) * dbeta/dr + beta*(1-mu^2)/r)
```

In the non-relativistic limit: `a_l ≈ (V/cR) * (1 + SIGMA*mu^2)` where `SIGMA = (R/V)*dV/dR - 1`.

### 2.2 Standard generalized opacity (Knop 2009a, Eq. 5)

```
chi_hat = chi + a*(4 + xi*p^|)
```

This can go negative when `-4a > chi + a*xi*p^|` (Eq. 6).

### 2.3 Knop's new generalized opacity (Eq. 7)

```
chi_hat = chi + xi * a * p^|
```

This is **always positive** because `a*p^| >= 0` by the upwind construction:
- For `a >= 0`: `p^| = lambda_l/(lambda_l - lambda_{l-1}) > 0`, so `a*p^| >= 0`
- For `a < 0`:  `p^| = lambda_l/(lambda_l - lambda_{l+1}) < 0`, so `a*p^| > 0` (neg * neg = pos)

The `4a*I` term is moved to the explicit source `S_tilde` (Eq. 10):
```
S_tilde = (a/chi_hat) * ([1-xi]*(p^- * I_{l-1} + p^+ * I_{l+1}) + [4 + (1-xi)*p^|] * I_l)
```

### 2.4 Formal solution (Eq. 8)

```
I(tau_2) = I(tau_1)*exp(tau_1-tau_2) + integral(S_hat * e^{tau-tau_2} dtau) - integral(S_tilde * e^{tau-tau_2} dtau)
```

Note the **minus sign** before the S_tilde integral.

### 2.5 Gauss-Seidel iteration (Knop 2009b, Eqs. 31-34)

For `a >= 0` (Eq. 31):
```
I_{i,l}^{n+1} = (1 - B_diag)^{-1} * (Delta_I + B_sub * I_{i,l-1}^{n+1} + A_sub * I_{i-1,l-1}^{n+1} + A_diag * I_{i-1,l}^{n+1} + C_sub * I_{i+1,l-1}^{n})
```

For `a < 0` (Eq. 32):
```
I_{i,l}^{n+1} = (1 - B_diag)^{-1} * (Delta_I + B_super * I_{i,l+1}^{n+1} + A_diag * I_{i-1,l}^{n+1} + A_super * I_{i-1,l+1}^{n+1} + C_super * I_{i+1,l+1}^{n})
```

The key: for `a < 0`, the upwind wavelength is `l+1` (next wavelength), and the GS iteration uses `I^{n+1}` for already-visited points and `I^{n}` for not-yet-visited points.

### 2.6 xi=1 requirement (Knop 2009a, Section 3)

The extreme alternating velocity test required `xi=1` (fully implicit) because `xi=0` was "ill-conditioned." The `xi=1` discretisation is diffusive but stable.

---

## 3. Implementation in fg_j_cmf_nonmono_v1.f

### 3.1 Generalized opacity

```fortran
! Line 1744-1745:
Q(1:NI) = GAM_MAG(1:NI,LS) / dLOG_NU
CHI_RAY(1:NI) = CHI_RAY(1:NI) + Q(1:NI)
```

where:
```fortran
! GAM = (V/R) * (1 + SIGMA*mu^2)  [= c * a_l in Knop's notation]
! GAM_MAG = |GAM|
! dLOG_NU = logarithmic frequency spacing
```

This gives: `chi_hat = chi + |GAM|/dLOG_NU`

**Assessment: CORRECT.** This is equivalent to Knop's `chi_hat = chi + xi*a*p^|` with `xi=1`, because:
- `|GAM|/dLOG_NU` corresponds to `|a*p^|` (the upwind diagonal coefficient)
- `a*p^|` is always positive (by the upwind construction), so `|a*p^| = a*p^|`
- Therefore `chi_hat = chi + a*p^|` (always positive)

### 3.2 Signed source coupling

```fortran
! Line 1756:
Q(1:NI) = GAM(1:NI,LS) / (dLOG_NU * CHI_RAY(1:NI))
SOURCE_RAY(1:NI) = ETA_RAY(1:NI) / CHI_RAY(1:NI)
```

This gives: `Q = GAM/(dLOG_NU * chi_hat)` (SIGNED), `S = eta/chi_hat`

**Assessment: CORRECT.** The signed Q retains the direction of wavelength coupling. `|Q| <= 1` always, because `|GAM|/dLOG_NU <= chi_hat` (it's part of chi_hat).

### 3.3 Sign-dependent upwind

```fortran
! Lines 1872-1886 (inward rays):
IF(GAM(I,LS) .GE. 0.0_LDP)THEN
    SOURCE_PRIME(I) = SOURCE_RAY(I) + Q(I)*I_M_PREV(I,LS)       ! expanding: upwind = previous freq
ELSE IF(HAS_COMPRESSING)THEN
    K = RAY_PNT(I,LS)
    IF(K .GE. 1 .AND. K .LE. ND .AND. FREQ_IDX_ITER+1 .LE. NCF_IN)THEN
        SOURCE_PRIME(I) = SOURCE_RAY(I) + Q(I)*I_M_DEPTH_NEXT(K,LS,FREQ_IDX_ITER+1)  ! compressing: upwind = next freq
    ELSE
        SOURCE_PRIME(I) = SOURCE_RAY(I) + Q(I)*I_M_PREV(I,LS)   ! fallback
    END IF
END IF
```

**Assessment: CORRECT.** This is the key Knop innovation:
- `GAM >= 0` (expanding): upwind = previous frequency = `I_PREV` (already computed in current sweep)
- `GAM < 0` (compressing): upwind = next frequency = `I_NEXT` (from previous great iteration, GS approach)

The criterion is correctly based on `GAM` (the ray-dependent coupling), NOT on `SIGMA` alone. This is important because `GAM = (V/R)*(1 + SIGMA*mu^2)` depends on both `SIGMA` and `mu`, so a point with `SIGMA < -1` might still have `GAM > 0` for tangential rays (small `mu`).

### 3.4 Gauss-Seidel storage and swap

```fortran
! Lines 1063-1069: At start of each great iteration (INIT):
IF(HAS_NEXT_ITER .AND. HAS_COMPRESSING)THEN
    I_P_DEPTH_NEXT(:,:,:) = I_P_DEPTH_CURR(:,:,:)   ! CURR from previous iteration becomes NEXT
    I_M_DEPTH_NEXT(:,:,:) = I_M_DEPTH_CURR(:,:,:)
END IF
IF(HAS_COMPRESSING)THEN
    I_P_DEPTH_CURR(:,:,:) = 0.0_LDP                  ! Zero CURR for new iteration
    I_M_DEPTH_CURR(:,:,:) = 0.0_LDP
END IF

! Lines 2078-2089: After formal solve, store current intensities:
IF(HAS_COMPRESSING)THEN
    DO LS=1,NP
        DO I=1,NI_RAY(LS)
            K=RAY_PNT(I,LS)
            IF(K .GE. 1 .AND. K .LE. ND .AND. ...)THEN
                I_P_DEPTH_CURR(K,LS,FREQ_IDX_ITER) = I_P(I,LS)
                I_M_DEPTH_CURR(K,LS,FREQ_IDX_ITER) = I_M(I,LS)
            END IF
        END DO
    END DO
END IF
```

**Assessment: CORRECT.** The GS swap is done once per great iteration:
- At iteration start: CURR (from previous iteration) → NEXT, CURR zeroed
- During formal solve: compressing points use NEXT (from previous iteration) as upwind
- After formal solve: current intensities stored in CURR
- On first iteration: NEXT = 0 (cold start), compressing points get `S' = S` (no coupling)

### 3.5 Intensity clamping

```fortran
! Lines 1922-1936:
IF(I_M(I+1,LS) .LT. 0)THEN
    N_NEG_INT = N_NEG_INT + 1
    I_M(I+1,LS) = 1.0E-30_LDP
    ...diagnostic output...
END IF
```

**Assessment: CORRECT and necessary.** The explicit source term `Q*I_NEXT` with `Q < 0` can drive the intensity slightly negative in compression zones. The clamping to 1e-30 is a physical bound (`I >= 0`) and does not affect the opacity (which is already guaranteed positive by `chi_hat = chi + |GAM|/dLOG_NU`).

---

## 4. Differences from the Published Knop Algorithm

### 4.1 The 4a sphericity term (DIFFERENCE, but NOT a bug)

**Knop 2009a:** The transfer equation includes `4*a*I`, which is moved from the opacity to the explicit source `S_tilde`.

**CMFGEN implementation:** The `4a` term is NOT present in the formal solver. The transfer equation used is:
```
chi * dI/dz = chi*I - eta - (GAM/dLOG_NU)*(I - I_prev)
```

This is because CMFGEN uses the Hauschildt (1992) formulation of the CMF transfer equation, which is written in terms of `dI/dln(nu)` rather than `d(lambda*I)/dlambda`. In this formulation, the sphericity term is absorbed into the definition of `a_l` and does not appear as a separate `4a*I` term.

**Evidence:** The standard V13 solver also does NOT include a 4a term:
```fortran
! V13 (fg_j_cmf_v13.f, line 1433-1434):
Q(1:NI) = GAM(1:NI,LS) / dLOG_NU
CHI_RAY(1:NI) = CHI_RAY(1:NI) + Q(1:NI)
```

Since V13 has been validated for monotonic velocity fields, the absence of the 4a term is consistent and correct within CMFGEN's formulation.

**Conclusion:** The non-monotonic solver is consistent with V13 and with CMFGEN's formulation. The 4a term is not missing — it's not part of this formulation.

### 4.2 Single GS sweep per great iteration (DIFFERENCE — convergence rate limitation)

**Knop 2009b:** The Gauss-Seidel iteration is performed **within** the formal solution call, iterating over all wavelengths multiple times until convergence. The "physically guided index field stepping scheme" follows the flow of information along each characteristic.

**CMFGEN implementation:** The GS iteration is performed **across** great iterations. One sweep through all frequencies is done per great iteration. The `I_NEXT` values from the previous great iteration are used as the GS approximation for compressing points.

**Assessment: VALID but SLOW.** The single-sweep approach is a legitimate GS iteration strategy — it converges to the same solution, but the convergence rate is determined by the spectral radius of the iteration matrix. For `|Q| ≈ 1` (extreme compression), the spectral radius is close to 1, and convergence is very slow.

The current convergence behaviour confirms this:
- Iteration 1: 1.36E+06 % increase (cold start, I_NEXT=0)
- Iteration 4: 1.25E+04 % increase
- Iteration 5: 1.95E+03 % increase
- Iteration 6: in progress

The corrections are decreasing by roughly 1 order of magnitude per 1-2 iterations, which is consistent with GS convergence for `|Q| ≈ 1`.

**Recommendation:** If faster convergence is needed, implement inner GS iteration within the formal solver (as in Knop 2009b, Section 3.2). This would require:
1. Storing intensities for all frequencies simultaneously
2. Iterating over frequencies multiple times within one formal solve
3. Convergence criteria for the inner iteration

This is a significant code change but would accelerate convergence substantially.

### 4.3 Wavelength discretisation (CONSISTENT with Knop)

**Knop 2009a:** Uses Crank-Nicholson mixing with parameter `xi ∈ [0,1]`. The extreme alternating velocity test required `xi=1` (fully implicit).

**CMFGEN implementation:** Uses first-order upwind, which is equivalent to `xi=1` (fully implicit). The discretisation is:
```
a * dI/dln(nu) ≈ (GAM/dLOG_NU) * (I - I_upwind)
```

where `I_upwind = I_PREV` for `GAM >= 0` and `I_NEXT` for `GAM < 0`.

**Assessment: CORRECT.** The `xi=1` choice is the stable one recommended by Knop for extreme velocity fields.

---

## 5. Diagnosis of the Earlier Blow-Up

### 5.1 Without LTE_EST (blow-up at iteration 2)

The first run without `LTE_EST` blew up at great iteration 2:
```
RJ(1) = 1.24E+36, RJ(ND) = 9.23E+58
```

**Root cause:** Outer-iteration instability, NOT formal-solver error.

The chain:
1. First great iteration: formal solver with `I_NEXT=0` produces moderate J
2. Population solver (SOLVEBA_V13) produces extreme populations (1.36E+06 % increase)
3. Second great iteration: formal solver uses extreme opacities/emissivities + `I_NEXT` from iteration 1
4. The `eta/chi_hat` ratio becomes extreme, driving J to blow up
5. Population solver amplifies further → divergence

The extreme sigma values (sigma_min ≈ -166) produce `|Q| ≈ 1` at many points, which means the source function is dominated by the `Q*I_NEXT` term. When `I_NEXT` is far from converged (which it is after the first iteration with extreme populations), the source can be very wrong.

### 5.2 With LTE_EST (converging)

The `LTE_EST=T` run is converging:
- Iteration 1: 1.36E+06 % → Iteration 5: 1.95E+03 % → Iteration 6: in progress

**Why LTE_EST helps:** LTE populations provide a stable initial guess for the population solver, preventing the extreme population excursions that drove the outer-iteration instability. The formal solver's GS iteration then converges slowly but steadily.

### 5.3 The formal solver itself is stable

Evidence:
1. The first great iteration (with `I_NEXT=0`) does not blow up — the formal solver produces finite J
2. With stable populations (LTE_EST), the formal solver converges over multiple iterations
3. The intensity clamping handles the occasional negative intensity from `Q*I_NEXT` with `Q < 0`
4. The opacity is always positive: `chi_hat = chi + |GAM|/dLOG_NU > 0`

---

## 6. Verification of USE_LAM_ES

### 6.1 Current setting: `USE_LAM_ES = T`

With `USE_LAM_ES=T`, after the formal solve:
```fortran
RJ(1:ND) = TC(1:ND)    ! TC = mean intensity from formal solution
```

The moment equations are NOT solved. `RJ` comes directly from the Knop formal solution.

### 6.2 Alternative: `USE_J_REL = T` (with `N_TYPE = N_ON_J`)

With `USE_J_REL=T`, `MOM_JREL_V9` is called, which solves the moment equations. However, `MOM_JREL_V9` clamps sigma:
```fortran
IF(SIGMA(I) .LT. -0.999_LDP) SIGMA(I) = -0.999_LDP
```

This clamping would **destroy the non-monotonic physics** by artificially preventing compression. The moment solver would then produce a different `RJ` than the formal solution, effectively discarding the Knop result.

### 6.3 Conclusion

`USE_LAM_ES=T` is the **correct setting** for using the Knop formal solution. It ensures that `RJ` comes directly from the non-monotonic formal solver, without being overwritten by a clamped-sigma moment solve.

The trade-off is that `USE_LAM_ES=T` bypasses the moment equations entirely, so the radiation field is not iterated with the moment solver. This is acceptable because the Knop formal solution already provides the correct radiation field for the given opacities and emissivities.

---

## 7. Summary of Findings

| Aspect | Knop paper | Implementation | Status |
|--------|-----------|----------------|--------|
| Generalized opacity | `chi_hat = chi + xi*a*p^|` | `chi_hat = chi + |GAM|/dLOG_NU` | CORRECT (equivalent with xi=1) |
| Signed source coupling | `Q = a*p^|/chi_hat` (signed) | `Q = GAM/(dLOG_NU*chi_hat)` (signed) | CORRECT |
| Upwind direction | Based on sign of `a_l` | Based on sign of `GAM` | CORRECT |
| Expanding upwind | `I_{l-1}` (previous wavelength) | `I_PREV` | CORRECT |
| Compressing upwind | `I_{l+1}` (next wavelength) | `I_NEXT` (from previous great iteration) | CORRECT |
| 4a term | In explicit source `S_tilde` | Not present (different formulation) | CONSISTENT with CMFGEN |
| xi discretisation | `xi=1` for extreme fields | First-order upwind (= xi=1) | CORRECT |
| GS iteration | Inner iteration within formal solve | One sweep per great iteration | VALID but SLOW |
| Opacity positivity | Guaranteed by `a*p^| > 0` | Guaranteed by `|GAM|/dLOG_NU > 0` | CORRECT |
| Intensity clamping | Not mentioned | `I >= 1e-30` | CORRECT (physical bound) |
| Cold start | Not specified | `I_NEXT = 0` on first iteration | CORRECT |

---

## 8. Recommendations

### 8.1 Short-term (no code changes needed)

1. **Continue the current LTE_EST=T run.** It is converging (iteration 6, corrections decreasing). Let it run to completion or until it reaches a steady state.

2. **Monitor convergence rate.** If corrections decrease by ~1 order of magnitude per 1-2 iterations, the model should converge in ~10-20 great iterations. If corrections plateau, the GS iteration may be too slow.

3. **Do NOT enable USE_J_REL.** The moment solver's sigma clamping would destroy the non-monotonic physics.

### 8.2 Medium-term (if convergence is too slow)

4. **Implement inner GS iteration** within the formal solver, as described in Knop 2009b Section 3.2. This would:
   - Store intensities for all frequencies simultaneously
   - Iterate over frequencies 2-5 times per formal solve
   - Use the "physically guided stepping" (follow the information flow)
   - Dramatically accelerate convergence

5. **Add relaxation for I_NEXT.** Instead of using `I_NEXT` directly from the previous great iteration, use a weighted average:
   ```
   I_NEXT_relaxed = alpha * I_NEXT_current + (1-alpha) * I_NEXT_previous
   ```
   with `alpha ≈ 0.5`. This is a standard GS acceleration technique.

### 8.3 Long-term (for production runs)

6. **Consider velocity smoothing.** The extreme sigma values (sigma = -166) are physically unrealistic for the numerical grid resolution. A mild smoothing of the velocity field at the grid scale would reduce the extreme sigma values without significantly changing the physical model.

7. **Add regression tests** for:
   - Positive and negative GAM
   - Positive chi_hat
   - Correct signed Q
   - Correct upwind intensity selection
   - Compression zones with multiple turning points
   - GS convergence over multiple iterations

8. **Document the formulation.** The absence of the 4a term should be explicitly documented as a consequence of CMFGEN's Hauschildt (1992) formulation, not as an approximation or omission.

---

## 9. File References

- Formal solver: `cmfgen_nonmono_dev_V3/plane/fg_j_cmf_nonmono_v1.f`
  - Opacity construction: lines 1744-1756
  - Inward ray source: lines 1872-1886
  - Outward ray source: lines 1956-1970
  - GS storage: lines 2078-2089
  - GS swap: lines 1063-1069
  - Intensity clamping: lines 1922-1942, 2020-2033

- Integration: `cmfgen_nonmono_dev_V3/new_main/mod_subs/comp_j_blank.f`
  - USE_LAM_ES branch: sets `RJ = TC`
  - Formal solver call: passes `SIGMA_TRUE` / `SIGMA_TRUE_EXT`

- Standard V13 solver: `cmfgen_nonmono_dev_V3/plane/fg_j_cmf_v13.f`
  - Opacity construction: lines 1433-1437

- Moment solver (not used): `cmfgen_nonmono_dev_V3/plane/mom_jrel_v9.f`
  - Sigma clamping: `IF(SIGMA(I) .LT. -0.999_LDP) SIGMA(I)=-0.999_LDP`

- Current configuration: `X_nomono_v3/config.yaml`
  - `LTE_EST: T`, `USE_LAM_ES: T`, `NON_MONO_SOLVER: T`
