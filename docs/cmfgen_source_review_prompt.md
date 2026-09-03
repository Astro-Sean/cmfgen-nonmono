# Prompt: Review of Local CMFGEN Source Modifications for Inaccuracies

You are reviewing a **forked copy of the CMFGEN radiative-transfer code**
(D. Hillier, University of Pittsburgh) for **physics and numerical
inaccuracies introduced by local modifications**. CMFGEN is a Sobolev-with-
continuum/CMF radiative transfer code for expanding stellar winds and
supernova ejecta. The upstream source is unmodified; **only the local
working-tree changes (uncommitted diff) are in scope.** Do not review the
upstream code itself — assume Hillier's original is correct unless a local
edit has touched it.

## How to see the changes

The code lives at `cmfgen_source/` (a git submodule). The local
modifications are **uncommitted** in the working tree. To see exactly what
was changed vs. upstream:

```bash
cd cmfgen_source
git diff                       # all local changes
git diff -- <file>             # one file
git log --oneline -5           # upstream history (unchanged)
```

The upstream `HEAD` is `dd1660a`. Everything in `git diff` is a local edit.
There are ~22 modified files; the physics-relevant ones are listed below
with a one-line description of what the edit does. **Read the diff, then
read enough surrounding context in each file to judge correctness** — many
of these edits interact (e.g. the Sobolev guard and the formal-solution
guard must be algebraically identical).

## Project context (essential for judging correctness)

1. **Application**: Type IIb / Type IIn supernova ejecta with
   **non-monotonic velocity profiles** (velocity inversions from shell
   interaction / double explosions). This is the *whole reason* most of
   these edits exist — upstream CMFGEN assumes monotonic Hubble-like flow
   and breaks for `SIGMA = dlnV/dlnR - 1 < -1`.

2. **Depth indexing**: CMFGEN indexes depth **outside-in**.
   `R(1) > R(2) > ... > R(ND)`. Depth 1 = outermost, depth ND = innermost.
   Keep this in mind when reading any loop that touches `R(I)` vs `R(I+1)`.

3. **SIGMA**: `SIGMA = dlnV/dlnR - 1`. For Hubble flow SIGMA=0. For a
   velocity *inversion* (deceleration then acceleration, or an inner shell
   overtaking outer material) SIGMA can go well below -1 (values down to
   **-22** occur in our models at the steepest part of the inversion).

4. **The Dessart+2015 guard** (the canonical fix for SIGMA<-1):
   ```
   T1 = 1 + SIGMA*mu^2
   T1 = MAX( SQRT(T1*T1 + 0.01), -0.2*(1+SIGMA), 0.1 )
   ```
   This replaces the raw `1 + SIGMA*mu^2` (which goes through zero / negative
   for SIGMA<-1) in the Sobolev optical-depth parameter GAM. The three terms
   handle: (a) `|sigma|>>1` → behaves like `|1+sigma*mu^2|` (transparent,
   correct); (b) `sigma << -1, mu~0` (tangential rays) → Dessart term
   `-0.2*(1+sigma)`; (c) `sigma ~ -1, mu ~ 1` → floor 0.1 prevents div-by-zero.

5. **CRITICAL — do NOT recommend removing the `-1.2 < SIGMA < -1.0 → -0.999`
   clamp** in `rv_sn_model_sniin.f`, `adjust_r_grid_v4.f`, or
   `do_ng_accel_band_v4.f`. This was tested empirically and caused
   catastrophic Ng-acceleration divergence (max % increase ~1e18 within
   ~100 iterations). The clamp is *required* for the Ng extrapolation
   polynomial to remain bounded, even though the Dessart guard makes the
   *local* Sobolev/GAM algebra safe. The clamp only fires in the narrow
   `-1.2..-1.0` band; SIGMA values below -1.2 pass through unclamped and
   are handled by the Dessart guard downstream. If you believe there is a
   problem in this region, propose alternatives (NG_BAND width, R-grid
   resolution around the SIGMA=-1 crossing, disabling Ng for cold-start
   iterations) — **do not propose removing the clamp.**

## Files to review (physics-relevant local edits)

### Sobolev optical-depth guard (6 files — must be mutually consistent)
The raw `GAM(I) = GAMH(I) / ABS(1 + (Z/R)^2 * SIGMA)` (or without ABS) was
replaced with the Dessart+2015 guard in all six Sobolev routines. **Check
that all six use the identical formula** (including the `+0.01` inside the
SQRT and the `-0.2*(1+SIGMA)` term), and that none accidentally dropped the
`ABS` semantics for the `sigma>0` regime where the original `ABS` mattered.

- `subs/sobew.f`
- `subs/sobew_grad.f`
- `subs/sobew_grad_v2.f`
- `subs/sobew_grad_v3.f`
- `subs/sobjbar_sim.f`
- `disp/subs/sobjbar.f`

### Formal solution: `plane/fg_j_cmf_v13.f` (and v12)
This is the most substantial and most important edit. Changes:
- Removed the `SIGMA_RAY <= -1.0 → -0.999` clamp (line ~942) so the real
  SIGMA feeds the Dessart guard.
- Added the Dessart guard to the GAM computation.
- **Re-clamped SIGMA to -0.999 *only* for the dGAMdR derivative and the
  half-point GAMH** (DIFFERENCE method), with a comment explaining that GAM
  itself is guarded but its derivative `dGAMdR = GAM*SIGMA/R` is not, and a
  raw SIGMA=-22 gives a pathological `dGAMdR` that corrupts DTAU and
  produces negative optical depths / negative mean intensities.

  **Specifically check**: is re-clamping the derivative while NOT clamping
  GAM internally consistent? The derivative of the *guarded* GAM is *not*
  `GAM*SIGMA/R` — it depends on which guard term is active (ABS / Dessart /
  floor). The clamp is a pragmatic substitute for the piecewise derivative.
  Is this substitute sound, or does it introduce a GAM/dGAMdR inconsistency
  worse than the one it fixes? Compare against the v12 file (which got the
  guard but NOT the derivative clamp) — are v12 and v13 now intentionally
  different, or is one of them wrong?

- `plane/fg_j_cmf_v12.f`: same guard added to GAM, clamp removed, but **no**
  derivative clamp. Confirm v12 is the older/alternate path and whether the
  asymmetry vs v13 is intended.

### Velocity-law / R-grid setup: `new_main/subs/rv_sn_model_sniin.f`
Substantial rewrite of the SN IIn R-grid and velocity interpolation. Changes:
- `NBND_INS` (number of inserted boundary zones) now passed in from VADAT
  via `N_OB_INS_IN` instead of hardcoded 3; upper bound changed from 3 to
  `ND/4`.
- **Bug fix**: `READ(LU,*) R_HYDRO(I),V_HYDRO,XD,...` → `V_HYDRO(I)`. The
  original read velocity into a scalar `V_HYDRO` (array element 0 /
  uninitialized) and left `V_HYDRO(I)` undefined. Confirm this is a genuine
  upstream bug and the fix is correct.
- **Removed the R_HYDRO / V_HYDRO reversal loop**, with a comment that the
  input `SN_HYDRO_DATA` is already written outer-to-inner (decreasing R),
  matching CMFGEN's depth convention. **Verify this**: does the input file
  actually come in decreasing-R order? If the input were increasing-R, this
  removal would silently reverse the velocity profile. Check
  `scripts/convert_stella_swd_to_cmfgen.py` (in the parent project) for the
  write order, or check the `RDINR` read path.
- New linear-extrapolation logic for `V` at the inserted boundary zones
  (inner and outer), using the local slope from `V_TEMP`. Check the
  indexing: `V(NBND_INS+I)=V_TEMP(I)` for `I=2..MND-1`, and the extrapolation
  loops `I=2..NBND_INS+1` (outer) and `I=ND-NBND_INS..ND-1` (inner). **Off-
  by-one errors here are the most likely bug** — verify the boundary zone
  indices line up with the R-grid insertion logic above (the
  `IF(NBND_INS .EQ. 3)THEN ... ELSE ... ENDIF` block that places
  `R(1+J), R(ND-J)`).
- `STOP` added after the invalid-`NBND_INS` error message (previously fell
  through). Confirm this is intentional.

### Statistical-equilibrium solver: `new_main/subs/solveba_v13.f`
- Tracks the species/level index `J` (JINC/JDEC) of the maximum STEQ
  correction, in addition to the depth `I`.
- Prints a diagnostic for the first 5 iterations only.
Check: the `INCREASE`/`DECREASE` sign convention. The code has
`IF(STEQ(J,I) .GT. DECREASE) DECREASE=STEQ(...)` and
`ELSE IF(STEQ .LT. INCREASE) INCREASE=STEQ(...)`. This means "increase" is
the most-negative STEQ and "decrease" is the most-positive — confirm this
matches CMFGEN's sign convention for STEQ (population correction
direction), and that the printed `INCREASE/100.0` scaling is meaningful.

### Population matrix: `new_main/subs/generate_full_matrix_v3.f`
- Wrapped the `C_MAT(I,J)=C_MAT(I,J)*POPS(J,L)` inner loop in
  `IF(ABS(POPS(J,L)) .GT. 0.0_LDP)`. This skips scaling when POPS is exactly
  zero. **Check**: is skipping (leaving C_MAT unscaled) the correct behavior
  when POPS=0, or should C_MAT be zeroed in that case? Mathematically
  `C_MAT * 0 = 0`, but skipping leaves the *old* C_MAT value. If C_MAT was
  pre-initialized to zero this is fine; if it carries prior values this is a
  bug. Read the surrounding code to determine C_MAT's initialization.

### Abundance/clumping init: `new_main/mod_subs/set_abund_clump.f`
- Changed hydro read from reversed `R_HYDRO(NHYDRO-I+1)` to in-order
  `R_HYDRO(I)` — same outer-to-inner assumption as `rv_sn_model_sniin.f`.
  **These two must agree on input ordering.**
- Added T interpolation from the hydro grid (log T in units of 10^4 K,
  clamped to [0.3, 1.0] i.e. 3000–10000 K) and an ED estimate
  `DENSITY/(mean atomic weight * amu)`, floored at 1.0.
  **Check units**: CMFGEN's internal ED is in cm^-3; is
  `DENSITY(K)/MEAN_ATOMIC_WEIGHT/ATOMIC_MASS_UNIT()` giving cm^-3? Confirm
  `DENSITY` is in g/cm^3 and `ATOMIC_MASS_UNIT()` is in g. Also check the
  `T` clamp range [0.3, 1.0] is in the right units (10^4 K → 3000–10000 K)
  and is physically sensible for SN ejecta initialization.

### Isotope population check: `new_main/subs/rd_sn_data.f`
- Added a guard `T3 = DENSITY/AT_MASS/1.66E-24 + T1` and only emits the
  "inconsistent total/isotope populations" error when the fractional
  difference exceeds 0.20 **and** either the species population or the
  isotope sum exceeds 1e-15 of the pure-element number density. This
  suppresses spurious warnings when both populations are negligible.
  **Check**: is `1.66E-24` consistent with the `ATOMIC_MASS_UNIT()` used
  elsewhere (the code mixes `1.66E-24_LDP` literals and
  `ATOMIC_MASS_UNIT()` calls)? Is `1.0E-15` a sensible threshold (absolute
  fraction of total element density)?

### Diagnostics (low-risk, but check for side effects)
- `subs/ltepop_wld_v2.f`: prints LTE population at depth 176 for Fe2.
- `subs/regrid_log_dc_v1.f`: prints DI at depth 176 for `Fe2_IN`.
- `misc/rewrite_old_model_data.f`: `RETURN` → `STOP` at end of program.
  Confirm this is a standalone tool program (where STOP is acceptable) and
  not a subroutine.
- `Makefile_definitions`: added `viper_intel_{nodebug,debug}` build blocks,
  changed default HOST/COMPILER to viper/intel, added `-lgfortran` to PGLIB.
  Check the `viper_intel_nodebug` LAPACK line uses `-lmkl_gf_lp64
  -lmkl_gnu_thread` (GNU MKL interface) — is that correct for ifort, or
  should it be `-lmkl_intel_lp64 -lmkl_intel_thread`? This is a likely
  build-system inaccuracy.

### Plumbing (low-risk, verify only)
- `new_main/cmfgen_sub.f`: passes `N_OB_INS` to `RV_SN_MODEL_SNIIN`.
- `new_main/mod_subs/rd_control_variables.f`: sets `SN_MODEL=.TRUE.` /
  `SN_HYDRO_MODEL=.TRUE.` for VELTYPE=12 (was missing).

## What to look for (in priority order)

1. **Algebraic/numerical errors** in the Dessart guard and dGAMdR: wrong
   sign, dropped `ABS`, missing `+0.01`, inconsistent formula between the 6
   Sobolev routines and `fg_j_cmf_v13.f`, or between v12 and v13.
2. **Off-by-one / index errors** in `rv_sn_model_sniin.f` boundary-zone
   insertion and velocity extrapolation — this is the highest-risk file.
3. **Input-ordering assumptions**: `rv_sn_model_sniin.f` and
   `set_abund_clump.f` both assume the hydro file is outer-to-inner
   (decreasing R). Verify against the writer
   (`scripts/convert_stella_swd_to_cmfgen.py`) and the reader (`RDINR`).
4. **Units consistency**: T in 10^4 K units, ED in cm^-3, the `1.66E-24`
   literal vs `ATOMIC_MASS_UNIT()`.
5. **The dGAMdR clamp vs guarded-GAM consistency** in `fg_j_cmf_v13.f` — is
   clamping the derivative to SIGMA=-0.999 while GAM uses the real SIGMA
   physically defensible, or does it create a GAM/dGAMdR mismatch?
6. **generate_full_matrix_v3.f**: is skipping the POPS=0 scaling correct
   given C_MAT's initialization?
7. **Build system**: MKL interface libs for ifort in
   `Makefile_definitions`.

## What NOT to flag

- Upstream CMFGEN code that the local diff did not touch.
- The `-1.2 < SIGMA < -1.0 → -0.999` clamp in `rv_sn_model_sniin.f`,
  `adjust_r_grid_v4.f`, `do_ng_accel_band_v4.f` (see context item 5).
- Style, formatting, or the diagnostic `WRITE` statements (unless they
  corrupt state).
- The mere presence of hardcoded depth 176 in the diagnostics.

## Output format

For each issue found, report:
- **File:line** (from the diff hunk or surrounding code)
- **Severity**: Critical (wrong physics / crash) / Major (wrong result,
  silent) / Minor (cosmetic, build-only)
- **Category**: algebra | indexing | units | input-ordering |
  build-system | other
- **What is wrong** (one sentence)
- **Why it is wrong** (the physics/numerics, 2-3 sentences)
- **Suggested fix** (concrete, minimal)

If a file's edits are correct, say so in one line — do not pad. End with a
short summary: number of Critical / Major / Minor issues, and whether the
Dessart guard is implemented consistently across all 7 files (6 Sobolev +
formal solution).
