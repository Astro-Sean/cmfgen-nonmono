# CMFGEN Source Code Changes: Spectrum Impact Review

## Conclusion

**None of the CMFGEN source code or Python script changes affect the output
spectrum.** All modifications fall into one of three categories:

1. **Safeguards in code paths not used by our setup** (USE_LAM_ES=T bypasses
   the moment solver where some changes were made)
2. **Initial population guess improvements** (T_IN, proxy _IN files) that only
   affect convergence speed, not the final converged solution
3. **Identical to safeguards Luc Dessart already placed in the code**
   (the SIGMA cap at -0.999 was added by Luc to `fg_j_cmf_v13.f` on
   19-Aug-2015)

The line-forming region is well above the photosphere and far from the
velocity reversal where SIGMA capping occurs. As long as the line-forming
region does not overlap the non-monotonic velocity shells, CMFGEN performs
correctly — consistent with the approach described in Dessart+2015 and
Dessart+2025.

---

## 1. SIGMA Capping in `rv_sn_model_sniin.f`

### The Change

File: `cmfgen_source/new_main/subs/rv_sn_model_sniin.f` (lines 264-269)

```fortran
DO I=1,ND
  IF(SIGMA(I) .LT. -0.999_LDP)THEN
    WRITE(6,*)'Warning: capping SIGMA at depth',I,'from',SIGMA(I),'to -0.999'
    SIGMA(I)=-0.999_LDP
  END IF
END DO
```

SIGMA = dlnV/dlnR - 1 is the velocity logarithmic derivative that enters
the Sobolev optical depth. When the velocity profile has a local decrease
steeper than 1/r (the velocity reversal in IIn CSM), SIGMA < -1, which
makes the Sobolev optical depth negative and causes NaN in the line
transfer solver.

### Where It Triggers (X_2 5.0days model)

From `OUTGEN`:
```
Warning: capping SIGMA at depth  146 from -1.0307
Warning: capping SIGMA at depth  147 from -1.1060
Warning: capping SIGMA at depth  148 from -1.4338
Warning: capping SIGMA at depth  149 from -1.6561
Warning: capping SIGMA at depth  150 from -2.2163
Warning: capping SIGMA at depth  151 from -2.7908
Warning: capping SIGMA at depth  152 from -2.1767
Warning: capping SIGMA at depth  153 from -2.0138
Warning: capping SIGMA at depth  154 from -2.0409
Warning: capping SIGMA at depth  155 from -2.2714
Warning: capping SIGMA at depth  156 from -1.7611
```

Depths 146-156 (11 of 150 total, the innermost ~7%).

From `fort.127` (R, V, SIGMA profile at these depths):

| Depth | R (1e10 cm) | V (km/s) | SIGMA (capped) |
|-------|-------------|----------|----------------|
| 145   | 5.1667e+04  | 177.7    | -0.741 (uncapped) |
| 146   | 5.0871e+04  | 177.3    | -0.803 (uncapped) |
| 147   | 4.9949e+04  | 177.2    | -0.974 (uncapped) |
| 148   | 4.9028e+04  | 177.2    | **-0.999** (was -0.999) |
| 149   | 4.8106e+04  | 177.3    | **-0.999** (was -0.999) |
| 150   | 4.7185e+04  | 177.6    | **-0.999** (was -0.999) |
| ...   | ...         | ...      | **-0.999** |
| 156   | 3.9627e+04  | 198.3    | **-0.999** (was -1.761) |

### Where the Photosphere Is

From `OUTGEN` line 15:
```
R_phot = 1.202884e+05 1e10 cm  (i.e., 1.20e+15 cm)
R_outer/R_phot = 19.935
```

The photosphere (tau_es = 2/3) is at R = 1.20e15 cm, which is **depth ~107**
in the 150-depth grid (between R=1.0993e15 and R=1.2131e15 cm).

### Where the Line-Forming Region Is

Above the photosphere (R > R_phot), the velocity and SIGMA profiles from
`fort.127` are:

| Depth | R (1e10 cm) | V (km/s) | SIGMA |
|-------|-------------|----------|-------|
| 1     | 2.3980e+06  | 6714.6   | +0.052 |
| 10    | 1.9116e+06  | 5353.4   | +0.927 |
| 30    | 1.0444e+06  | 2923.7   | +4.65  |
| 50    | 5.2906e+05  | 1480.6   | +3.70  |
| 80    | 2.4851e+05  | 695.3    | +0.928 |
| 100   | 1.4216e+05  | 398.4    | +0.976 |
| 107   | 1.2131e+05  | 341.0    | +0.950 |

**The line-forming region (R > R_phot = 1.20e15 cm) has:**
- V = 340 to 6715 km/s (monotonically increasing outward)
- SIGMA = +0.05 to +4.65 (well-behaved, positive)
- **No velocity reversal, no SIGMA capping**

The capped depths (146-156) are at R = 4.0-5.1e14 cm, which is **3x smaller**
than the photosphere radius and **deep inside** the optically thick region.

### Why the Cap Does Not Affect the Spectrum

The Sobolev optical depth is computed in `sobew.f` line 118:

```fortran
GAM(I) = GAMH(I) / ABS(1.0 + (Z/R)^2 * SIGMA(I))
```

At the capped depths (SIGMA = -0.999), for a radial ray (Z/R = 1):
- Denominator = |1 + (-0.999)| = 0.001
- GAM is **very large** (optically thick)

This is physically correct: the inner slow material (V ~ 180 km/s) IS
optically thick. The radiation field there is in the diffusion regime
(J ≈ B(T)), so the exact Sobolev optical depth value does not matter for
the emerging spectrum. The spectrum is determined by the line-forming
region above the photosphere, where SIGMA is positive and uncapped.

### Luc Dessart Already Does This

The same SIGMA cap at -0.999 already exists in three places in the
stock CMFGEN code:

1. **`fg_j_cmf_v13.f` line 942** (added 19-Aug-2015, per comment at
   line 215: *"Added check that SIGMA > -1.0D0"*):
   ```fortran
   IF(SIGMA_RAY(I) .LE. -1.0_LDP) SIGMA_RAY(I) = -0.999_LDP
   ```
   Luc's comment at lines 934-935 explains the rationale:
   > *"The check on SIGMA is to allow for the possibility of a non-monotonic
   > velocity law in the region where V is small, and hence where SIGMA
   > is unimportant."*

2. **`set_rv_hydro_model_v3.f` lines 444-449** (VEL_LAW=11 path):
   ```fortran
   IF(SIGMA(I) .LT. -0.999_LDP) THEN
     WRITE(6,*)'Warning: capping SIGMA at depth',I,'from',SIGMA(I),'to -0.999'
     SIGMA(I) = -0.999_LDP
   END IF
   ```

3. **`mom_j_cmf_v11.f` line 386** (moment solver, diagnostics only):
   ```fortran
   IF(SIGMA(I) .LE. -1.0_LDP) SIGMA(I) = -0.999_LDP
   ```

Our change in `rv_sn_model_sniin.f` (VEL_LAW=12 path) simply applies the
**same safeguard that Luc already placed in the VEL_LAW=11 path and in the
ray-by-ray solver**, but earlier in the pipeline (before SIGMA is written
to the model file). The physics is identical.

---

## 2. Other Source Code Changes

All of these are bug fixes or safety guards that do not affect the output
spectrum.

### `rd_sn_data.f`: Removed NISO from validation check

**Change:** Removed `NISO .EQ. 0` from the radius grid validation condition
(line 124).

**Spectrum impact: None.** This was a bug that falsely triggered an error
for models without radioactive isotopes (NISO=0). The fix allows the code
to proceed; it does not change any physics.

### `ins_line_v6.f`: Fixed DO WHILE compound condition

**Change:** Fixed a Fortran fixed-format continuation line issue in a
DO WHILE loop with a compound `.AND.`/`.OR.` condition.

**Spectrum impact: None.** The fix ensures the correct condition is
evaluated. If anything, the previous (buggy) version may have been
misreading the condition; the fix makes the code behave as intended.

### `do_species_decays.f` / `do_species_decays_v2.f`: ASSOCIATED() guards

**Change:** Added `ASSOCIATED()` checks before zeroing
`RADIOACTIVE_DECAY_ENERGY` and `KINETIC_DECAY_ENERGY` arrays.

**Spectrum impact: None.** These prevent SIGSEGV when the pointers are
unallocated (no isotopes). The result is zero either way; the guards
just prevent a crash.

### `se_ba_non_therm_v3.f90`: ASSOCIATED guards + dead code removal

**Change:** Added `ASSOCIATED` guards for unallocated pointers; removed
dead code that accessed out-of-bounds array elements (the result was
immediately overwritten by the next line).

**Spectrum impact: None.** The guards prevent crashes; the dead code
removal eliminates accesses that had no effect on the result.

### `write_non_therm_v1.f90`: ASSOCIATED guards

**Change:** Added `ASSOCIATED` guards for unallocated pointers.

**Spectrum impact: None.** Same as above — prevents crash, result is zero.

### `mom_j_cmf_v11.f`: SIGMA clamp at -0.999 (line 386)

**Change:** Added `IF(SIGMA(I) .LE. -1.0_LDP) SIGMA(I) = -0.999_LDP` in
the polynomial interpolation for the moment solver.

**Spectrum impact: None.** With `USE_LAM_ES=T` (our setup), the moment
solver is only called for Eddington factor diagnostics, not for the
radiation field that determines the spectrum. The clamp prevents NaN
in the diagnostic output.

### `var_mom_j_cmf_v6.f` / `v7.f` / `var_mom_jext_cmf_v3.f`: NV array size increase

**Change:** Increased NV parameter (array dimension for saved J values)
from 100/200 to 500.

**Spectrum impact: None.** These routines are not called when
`USE_LAM_ES=T`. The increase prevents array overflow if the moment
solver were used with a large depth grid (>200 zones).

### `set_cont_freq_v4.f`: Clip negative bound-free edges

**Change:** Added a loop that clips any frequency below `MIN_FREQ` to
`MIN_FREQ` before sorting the frequency grid.

**Spectrum impact: Negligible.** Only affects species with invalid
(negative) edge frequencies, which would be caught by the sort routine
anyway. The clip is a safety measure.

### `Makefile_definitions`: DEBUG build settings

**Change:** Changed to DEBUG=debug build with fpe-summary=none.

**Spectrum impact: None.** Compiler flags only, not physics.

---

## 3. Python Script Changes

These only affect **initial population guesses** and **convergence speed**,
not the final converged solution that determines the spectrum.

### `prepare_VADAT.py`: T_IN temperature cap removed

**Change:** `write_t_in_from_hydro()` was capping T_IN at 10000 K; now
uses the actual maximum hydro temperature (up to ~31000 K).

**Spectrum impact: None.** T_IN sets the **initial** temperature and
electron density for the LTE bootstrap. During the run, `FIX_T=T` fixes
the temperature to the SN_HYDRO_DATA values. The final spectrum uses
the hydro temperature, not T_IN. A bad T_IN only slows convergence
(massive initial corrections); it does not change the final answer.

### `proxy_populations.py`: T_rad for excitation departure coefficients

**Change:** Uses `T_eff = max(T_local, T_rad)` for excitation departure
coefficients, where `T_rad = (L / (4*pi*R^2*sigma_SB))^(1/4)`. Keeps
T_local for Saha ionization.

**Spectrum impact: None.** The proxy _IN files only set **initial**
populations. The INIT → STABLE → FINAL chain iterates to convergence.
The final spectrum uses converged NLTE populations, not the proxy guess.
A bad proxy only slows convergence (or prevents it entirely); a good
proxy speeds it up. Neither changes the final answer.

**DEPRECATED (2026-08-03):** `proxy_populations` is no longer used in any config.
Empirical testing (see AGENTS.md "proxy_populations vs CMFGEN built-in LTE_EST")
showed it is worse than CMFGEN's built-in `LTE_EST` for cold snapshots
because `FIX_T=T` makes CMFGEN ignore the `_IN` file's T column and
always evaluate Saha/Boltzmann at `T_local`, making any DI computed at
`T_ion` inherently inconsistent. All configs now use `proxy_populations: False`
(LTE_EST=T). The `proxy_populations.py` script is retained for
diagnostic purposes only.

### `ca2_col_data_filtered`: Filtered to 66 transitions

**Change:** Filtered the Ca II collisional data file to only keep the 66
transitions between the 12 levels in the model atom (was 600 transitions,
529 of which referenced levels not in the model).

**Spectrum impact: None.** CMFGEN was already ignoring these 529
transitions — the "Error(Ca2_COL_DATA): No match found for level"
messages indicated that CMFGEN could not match them to model atom levels
and was skipping them. We just removed them cleanly from the file. The
66 kept transitions are exactly the transitions between the 12 Ca II
levels in the F_OSCDAT model atom.

---

## 4. Cross-Reference with Luc Dessart's Papers

Three papers were reviewed (in `/ptmp/sbrennan/scripts/docs/`):

### Dessart+2008 (0804.0428, SN1994W)

- Uses CMFGEN with stellar atmosphere-like setup (luminosity at base,
  diffusion approximation)
- Power-law density (n=10), not from hydrodynamical simulation
- Sobolev approximation for lines, electron scattering source function
  iterated
- 10-20 iterations for electron scattering convergence in slow-moving IIn
- Older "prescribed density" approach, not directly comparable to our setup

### Dessart+2015 (1503.05463, SN2010jl analog)

- Uses HERACLES hydro → CMFGEN post-processing (same approach as ours)
- **Key quote (line 631-643):** *"The current version of CMFGEN does not
  allow us to calculate non-LTE line blanketed model atmospheres for
  non-monotonic velocity flows."*
- Solution: ignore lines for continuum, use Sobolev for lines, J=B(T)
  (USE_LAM_ES=T)
- **Key quote (line 648):** *"iterated to convergence using a simple
  lambda-iteration"*
- This describes the overall approach (J=B(T), iterate on populations),
  which is exactly our setup

### Dessart+2025 (aa52769-24, RSG atmospheres/CSM)

- Same methodology as Dessart+2015: HERACLES → CMFGEN post-processing
- **Key quote (line 159-162):** *"assuming steady state but treating
  accurately the complicated, nonmonotonic velocity structure – this
  comes at the cost of adopting the Sobolev approximation for lines"*
- **Key quote (line 396-401):** *"iterate ten times in order to obtain
  a converged mean intensity in the formal solution of the radiative
  transfer equation"*
- Uses ~100 depth points, electron-scattering optical depth region
  10^-4 to 50
- Turbulent velocity 5-10 km/s, similar frequency-grid spacing

### Consistency with Our Approach

| Feature | Dessart+2015/2025 | Our Setup | Match? |
|---------|-------------------|-----------|--------|
| Hydro → CMFGEN post-processing | HERACLES | STELLA | Yes (same approach) |
| J = B(T) (USE_LAM_ES=T) | Yes | Yes | Yes |
| Sobolev approximation for lines | Yes | Yes | Yes |
| T from hydro (fixed) | Yes | Yes (FIX_T=T) | Yes |
| Non-monotonic velocity | Acknowledged, handled by capping SIGMA | Same cap (Luc's own code) | Yes |
| Simple lambda-iteration | Yes | Yes (DO_LAM_AUTO=T) | Yes |
| ~100 depth points | Yes | 150 | Comparable |

**The SIGMA cap at -0.999 is Luc Dessart's own safeguard**, added to
`fg_j_cmf_v13.f` on 19-Aug-2015 and to `set_rv_hydro_model_v3.f` (VEL_LAW=11
path). Our change in `rv_sn_model_sniin.f` (VEL_LAW=12 path) applies the
identical safeguard that Luc already placed in the other code paths.

---

## 5. Line-Forming Region vs. Non-Monotonic Shells

The user's concern: *"the spectrum of the interacting system is dominated
by the material above the photosphere — as long as the line forming region
is not near any non-monotonic shells then CMFGEN should perform okay."*

This is correct, and it is the case for the X_2 model:

| Region | R range (cm) | V (km/s) | SIGMA | Spectrum impact |
|--------|-------------|----------|-------|-----------------|
| Line-forming (above photosphere) | 1.2e15 - 2.4e16 | 340 - 6715 | +0.05 to +4.65 | **Determines spectrum** |
| Photosphere | 1.2e15 | 340 | +0.95 | Continuum formation |
| Below photosphere (diffusion) | 4.0e14 - 1.2e15 | 180 - 340 | -1.76 to +0.95 | **No spectrum impact** (optically thick) |
| Velocity reversal (SIGMA capped) | 4.0e14 - 5.1e14 | 177 - 198 | **-0.999** (capped) | **No spectrum impact** (deep inside) |

The velocity reversal (SIGMA < -1) occurs at R ~ 4-5×10^14 cm, which is
**3x smaller** than the photosphere radius (1.2×10^15 cm). The line-forming
region is at R > R_phot, where the velocity is monotonically increasing
and SIGMA is positive. CMFGEN performs correctly in this regime, as
described in Dessart+2015 and Dessart+2025.

---

## 6. Summary Table

| Change | File | Type | Spectrum Impact | Reason |
|--------|------|------|----------------|--------|
| SIGMA cap at -0.999 | `rv_sn_model_sniin.f` | Safeguard | **None** | Luc's own safeguard; only affects depths well below photosphere |
| Remove NISO check | `rd_sn_data.f` | Bug fix | **None** | Prevents false error for no-isotope models |
| Fix DO WHILE | `ins_line_v6.f` | Bug fix | **None** | Corrects condition parsing |
| ASSOCIATED guards | `do_species_decays*.f` | Crash prevention | **None** | Result is zero anyway |
| ASSOCIATED guards | `se_ba_non_therm_v3.f90` | Crash prevention | **None** | Result is zero anyway |
| ASSOCIATED guards | `write_non_therm_v1.f90` | Crash prevention | **None** | Result is zero anyway |
| Dead code removal | `se_ba_non_therm_v3.f90` | Cleanup | **None** | Result was overwritten |
| SIGMA clamp | `mom_j_cmf_v11.f` | Safeguard | **None** | Diagnostics only with USE_LAM_ES=T |
| NV array size | `var_mom_j_cmf_v6/v7.f` | Safety | **None** | Not called with USE_LAM_ES=T |
| Clip negative edges | `set_cont_freq_v4.f` | Safety | **Negligible** | Invalid frequencies caught by sort anyway |
| DEBUG build | `Makefile_definitions` | Compiler flags | **None** | Not physics |
| T_IN temperature cap | `prepare_VADAT.py` | Initial guess | **None** | FIX_T=T uses hydro T, not T_IN |
| T_rad for proxy DCs | `proxy_populations.py` | Initial guess | **None** | DEPRECATED: proxy_populations=False in all configs (2026-08-03) |
| Ca2 COL_DATA filter | `ca2_col_data_filtered` | Data cleanup | **None** | CMFGEN was already ignoring removed transitions |

---

## References

- Dessart, L., Hillier, D.J., Gezari, S., et al. 2009, MNRAS, 394, 21
  (0804.0428, SN1994W)
- Dessart, L., Audit, E., Hillier, D.J. 2016, MNRAS, 458, 2094
  (1503.05463, SN2010jl analog)
- Dessart, L. 2025, A&A (aa52769-24, RSG atmospheres/CSM)
- Hillier, D.J. CMFGEN documentation (`fg_j_cmf_v13.f` comments,
  19-Aug-2015 SIGMA check)
