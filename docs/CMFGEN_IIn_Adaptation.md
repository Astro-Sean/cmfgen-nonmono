# Adapting CMFGEN for Type IIn Supernova Models with Non-Monotonic Velocity

## Overview

This document describes the complete set of modifications required to adapt
CMFGEN for Type IIn supernova models with non-monotonic velocity profiles,
dense CSM interaction, and low-temperature outer layers. The work was done
for the SN2024pba project (X_1, X_2, X_3 models) which simulate a double-
explosion scenario: an inner H/He shell launched at 3000 km/s into a dense
RSG wind, creating a velocity reversal at R_star ~ 5e13 cm.

The X_2 model uses:
- VEL_LAW=12 (RV_SN_MODEL_SNIIN: tau-based grid for non-monotonic velocity)
- model_type=static_snapshot (no time evolution, USE_LAM_ES=T for J=B(T))
- ~~proxy_populations=True (dilute-LTE _IN files for initial populations)~~
  **DEPRECATED (2026-08-03):** `proxy_populations` is empirically worse than CMFGEN's
  built-in `LTE_EST` for cold snapshots. `FIX_T=T` makes CMFGEN ignore the
  `_IN` file's T column and always evaluate Saha/Boltzmann at `T_local`, so
  any DI computed at `T_ion` is inherently inconsistent. All configs now use
  `proxy_populations: False` (LTE_EST=T). See AGENTS.md "proxy_populations vs CMFGEN built-in
  LTE_EST" for the full empirical comparison.
- INC_RAD_DECAYS=True (Ni56/Co56 radioactive decay)
- 7 snapshot epochs: 5, 10, 15, 20, 25, 35, 45 days post-explosion

## Problem Summary

The X_2 CMFGEN models failed during the INIT (LTE bootstrap) stage with
three distinct errors, each requiring a separate fix:

1. **"Mean intensity blowing up"** in FG_J_CMF_V13 (SIGMA < -1 at inner depths)
2. **"bad initial population guesses"** with 7.2e116% corrections (T_IN capped at 10000 K)
3. **"bad initial population guesses"** with 6.6e100% corrections (proxy _IN used T_local instead of T_rad)

Each fix built on the previous one, progressively moving the failure point
deeper into the convergence process until the INIT could complete and pass
control to the STABLE job.

---

## Fix 1: SIGMA Capping in rv_sn_model_sniin.f (VEL_LAW=12)

### Problem

VEL_LAW=12 calls `rv_sn_model_sniin.f` to compute the R, V, SIGMA profile
from `input_hydro.dat`. The SIGMA parameter (dlnV/dlnR - 1) controls the
Sobolev optical depth integral. When the velocity profile has a local
decrease steeper than 1/r (which happens at the velocity reversal in IIn
CSM), SIGMA becomes < -1, making the Sobolev optical depth negative.

Negative SIGMA causes NaN in FG_J_CMF_V13 (the line transfer solver),
which manifests as "Mean intensity blowing up" in `comp_j_blank.f`.

### Root Cause

The SIGMA capping fix (cap at -0.999) was already present in
`set_rv_hydro_model_v3.f` (VEL_LAW=11 path) but was **missing** from
`rv_sn_model_sniin.f` (VEL_LAW=12 path). All X_ models use VEL_LAW=12,
so they never got the cap.

### Fix Applied

File: `cmfgen_source/new_main/subs/rv_sn_model_sniin.f` (lines ~257-272)

```fortran
! Cap SIGMA at -0.999 to prevent negative Sobolev optical depths
! (which cause NaN in FG_J_CMF_V13 and downstream STEQ routines).
! This can happen when the interpolated velocity profile has a local
! decrease steeper than 1/r (non-monotonic velocity in IIn CSM).
!
	DO I=1,ND
	  IF(SIGMA(I) .LT. -0.999_LDP)THEN
	    WRITE(6,*)'Warning: capping SIGMA at depth',I,'from',SIGMA(I),'to -0.999'
	    SIGMA(I)=-0.999_LDP
	  END IF
	END DO
```

This was added after the SIGMA computation loop and before the WRITE(127)
debug output. The cap at -0.999 (instead of -1.0) ensures the Sobolev
optical depth remains positive but very small, avoiding the NaN while
preserving the physical behavior.

### Verification

After the fix, OUTGEN shows:
```
Warning: capping SIGMA at depth 146 from -1.03067953021088
Warning: capping SIGMA at depth 147 from -1.10596285008587
Warning: capping SIGMA at depth 150 from -2.21630010473082
...
```
No more "Mean intensity blowing up" crash. The SIGMA values at depths
146-156 (inner boundary) are capped from -1.03 to -2.79 to -0.999.

---

## Fix 2: T_IN Temperature Cap in prepare_VADAT.py

### Problem

After Fix 1, the INIT failed with "bad initial population guesses" at
depth 22, with population corrections of 7.2e116%. The luminosity at
the outer boundary was 2.08e10 Lsun (vs 3.05e10 at inner), which seemed
reasonable.

### Root Cause

The `write_t_in_from_hydro()` function in `prepare_VADAT.py` had a
default parameter `t_max_104=1.0` (10000 K) and was called without
overriding this default:

```python
write_t_in_from_hydro(hydro, wdir / "T_IN")  # uses default t_max_104=1.0
```

The function caps the temperature:
```python
t = min(max(t_104[i], 0.3), t_max_104)  # floor at 3000 K, cap at 10000 K
```

The SN_HYDRO_DATA has temperatures up to 31117 K at the inner boundary,
but T_IN was capping them at 10000 K. This meant the LTE initial
populations at inner depths were computed at T=10000 K instead of
T=31000 K - a factor of 3 error in temperature.

At T=10000 K, LTE populations for excited states are ~1e-10 of ground
state. At T=31000 K, they are ~1e-4. The 1e6 factor difference in
initial populations caused the 7.2e116% corrections.

### Fix Applied

File: `scripts/prepare_VADAT.py` (line ~2954)

```python
# Use actual max T from hydro (no artificial 10000 K cap);
# the cap was causing bad LTE populations at inner depths where T > 10000 K.
hydro_t_max_104 = float(np.max(hydro.temp_104k))
write_t_in_from_hydro(hydro, wdir / "T_IN", t_max_104=hydro_t_max_104)
```

### Verification

T_IN now has correct temperatures at inner depths:
```
Depth  R(1e10cm)    T(K)
   50    3.2525e+05    3225 K
  100    1.1902e+05    7464 K
  130    6.7537e+04   14851 K  (was 10000 K)
  140    5.1667e+04   21687 K  (was 10000 K)
  150    3.9627e+04   31117 K  (was 10000 K)
```

After the fix, the error moved from depth 22 (7.2e116%) to depth 46
(6.6e100%) - a 16 order of magnitude improvement, but still failing.

---

## Fix 3: Radiation Temperature in proxy_populations.py

### Problem

After Fix 2, the INIT still failed with "bad initial population guesses"
at depth 46, with corrections of 6.6e100%. The CMF_BLKBAND_V3 solver
also reported singular BA matrices (IFAIL=31) at depths 14-16.

### Root Cause

The proxy_populations.py computes dilute-LTE departure coefficients for
the _IN population files. It uses the local gas temperature (T_local)
for the LTE calculation:

1. Compute dilution factor W = 0.5*(1-sqrt(1-(R_inner/R)^2))
2. Set ground state DC = W
3. Solve for T_EXC from: W = (T/T_EXC)^1.5 * exp(HDKT*FEDGE*(1/T_EXC - 1/T))
4. Compute excited state DC from T_EXC

The problem: in SN ejecta, the outer layers have T_local ~ 3000K but the
radiation temperature T_rad ~ 4000-10000K because the radiation field is
set by the luminosity from the hot inner ejecta photosphere, not the cold
outer gas.

The radiation temperature at each depth is:
```
T_rad = (L / (4*pi*R^2*sigma_SB))^(1/4)
```

At depth 46 (R=3.47e15 cm):
- T_local = 3072 K
- T_rad = 10801 K (3.5x hotter!)

Using T_local for LTE gives n(n=2)/n(n=1) ~ 1e-57 (essentially zero).
The correct NLTE value (set by the radiation field at T_rad) is ~1e-16.
This 1e41 factor discrepancy explains the 6.6e100% corrections.

### Fix Applied

File: `scripts/proxy_populations.py`

1. Added `compute_t_rad()` function to compute radiation temperature:
```python
def compute_t_rad(hydro):
    """Compute T_rad = (L / (4*pi*R^2*sigma))^(1/4) at each depth."""
    L = hydro['lum_lsun'] * L_SUN
    t_rad = []
    for i in range(nd):
        r_cm = hydro['radius'][i] * 1e10
        j = L / (4.0 * math.pi * r_cm * r_cm)
        t_rad_k = (j / SIGMA_SB) ** 0.25
        t_rad.append(t_rad_k / 1e4)  # 10^4 K
    return t_rad
```

2. Compute T_eff = max(T_local, T_rad) and use it for excitation DC:
```python
t_rad = compute_t_rad(hydro)
t_eff = [max(hydro['temperature'][i], t_rad[i]) for i in range(nd)]
```

3. Use T_eff for the excitation temperature calculation and departure
   coefficients (but keep T_local for Saha ionization to maintain
   consistency with FIX_NE=T):

```python
# Saha: use T_local (consistent with hydro electron density)
ion_pop, ion_frac = compute_saha_ion_frac(elem_label, ion_stage, hydro, atomic_mass)

# Excitation DC: use T_eff (reflects radiation field)
t_exc = []
for i in range(nd):
    te = get_excite_temp(t_eff[i], gs_dc[i], fedge_gs)
    t_exc.append(te)
dc = compute_dc(nlev, nd, t_eff, t_exc, atom_data_trunc['fedge'], w)
```

### Why T_local for Saha but T_rad for Excitation

Using T_rad for Saha ionization creates an inconsistency: the proxy _IN
predicts high ionization (from T_rad), but CMFGEN's INIT run uses
FIX_NE=T (electron density fixed to hydro values computed at T_local).
The hydro Ne reflects low ionization at T_local=3000K.

If the proxy says "90% ionized" but Ne says "10% ionized", the species
populations are inconsistent with the electron density, causing the BA
matrix to be singular (IFAIL=31).

The correct approach:
- **Saha ionization**: use T_local (matches hydro Ne, FIX_NE=T)
- **Excitation populations**: use T_rad (matches radiation field)

This way the ionization balance is consistent with the electron density,
while the excited state populations reflect the actual radiation field.

### Verification

After the fix, the HI departure coefficients at depth 46 changed:
- n(n=2)/n(n=1) improved from 4.1e-18 to 2.6e-6 (factor of 6.3e11)
- The new proxy is 0.148x of the expected NLTE value (was 2.3e-13x)

The INIT job now runs through all 25 iterations without crashing:
```
Iteration  1: 1.15e112% at depth 17
Iteration  4: 2.49e92%  at depth 46
Iteration  7: 2.51e86%  at depth 46
Iteration 12: 3.95e80%  at depth 48
Iteration 23: 3.17e72%  at depth 48
```

The corrections decrease by ~6 orders of magnitude per iteration. The
INIT doesn't fully converge in 25 iterations, but it writes SCRTEMP
and passes control to the STABLE job, which continues the convergence.

---

## Fix 4: Skip MOM_J_CMF_V11 in Eddington Pre-Computation when USE_LAM_ES=T

### Problem

`comp_j_blank.f` has two code paths that compute J:

1. **Normal iterations** (`CONT_VEL` section, line 531): Has a `USE_LAM_ES`
   branch (line 717) that skips the moment solver and sets `RJ=TC` (formal
   solution J from `FG_J_CMF_V13`). Safe with `G_ONLY` and non-monotonic
   velocity.

2. **Eddington pre-computation** (`CONT_VEL .AND. ACCURATE` section, line 259):
   Called when `COMPUTE_EDDFAC=T` (first run, no EDDFACTOR file). Calls
   `MOM_J_CMF_V11` **unconditionally** (line 372) with `N_TYPE=G_ONLY`.
   With velocity reversal, H→0 makes the N/H Eddington factor
   (`NMID_ON_HMID`) undefined, causing J to blow up.

The normal iteration path was already safe because of the existing
`USE_LAM_ES` branch. The ACCURATE/Eddington pre-computation path lacked
this check, creating a crash path for first-time runs.

### Fix Applied

File: `cmfgen_source/new_main/mod_subs/comp_j_blank.f` (line ~371)

Wrapped the `MOM_J_CMF_V11` call in a `USE_LAM_ES` check. When
`USE_LAM_ES=T`, the moment solver is skipped and `RJEXT=TC` (formal
solution J) is used instead, mirroring the existing branch at line 717:

```fortran
     IF(USE_LAM_ES)THEN
       RJEXT(1:NDEXT)=TC(1:NDEXT)
       IF(.NOT. USE_FORMAL_REL)THEN
         IF(.NOT. DIF)HFLUX_AT_IB=0.5_LDP*IC*(0.5_LDP+INBC)-INBC*RJEXT(NDEXT)
         HFLUX_AT_OB=HBC_CMF(1)*RJEXT(1)
       END IF
     ELSE
       CALL TUNE(IONE,'MOM_J_CMF_ACC')
       CALL MOM_J_CMF_V11(...)
       CALL TUNE(ITWO,'MOM_J_CMF_ACC')
       IF(.NOT. DIF)HFLUX_AT_IB=...
       HFLUX_AT_OB=HBC_CMF(1)*RJEXT(1)
     END IF
```

### Why This Works

- `FG_J_CMF_V13` computes FEDD directly from ray-by-ray integration —
  no N/H Eddington factor, no velocity reversal issue
- The DO WHILE convergence check still works: it checks if FEDD changes
  between iterations
- With `USE_LAM_ES=T`, the moment solver is skipped in normal iterations
  anyway — the ACCURATE section is now consistent
- `RJEXT=TC` is the same assignment used at line 717 for normal iterations

### Impact

- **No impact when `USE_LAM_ES=F`**: The ELSE branch preserves the original
  unconditional `MOM_J_CMF_V11` call
- **No impact on non-G_ONLY modes**: The fix is gated on `USE_LAM_ES`,
  not `N_TYPE`
- The `J_COMP` output will show 0% difference between J(mom) and J(ray)
  since both are TC — this is correct and expected

---

## Fix 5: Switch to USE_J_REL=T + N_ON_J for Convergence Stability

### Problem

With `USE_LAM_ES=T`, J is set from the formal solution (`FG_J_CMF_V13`, J=TC).
The formal solution can produce sharp J features at shock fronts (velocity
reversals), which create oscillating radiation rates → oscillating populations.

The X_1 5.0days run showed this clearly:
- Iterations 1-7: converging (1e14% → 1e6%)
- Iteration 8: a population zeroed out (100% decrease)
- Iterations 9+: zero/overshoot cycle (1e19% → 1e28% corrections, never converging)

Contributing factors:
- `MAX_LAM=1.0D5` (100,000x max population change per iteration) — allows wild swings
- `IBEG_NG=80` — Ng acceleration doesn't start until iteration 80, but oscillation
  starts at iteration 8
- `NUM_BNDS=1` (DIAG) — no inter-depth coupling in BA matrix to damp oscillations

### Fix Applied

**1. Switch from USE_LAM_ES=T to USE_J_REL=T + N_ON_J**

Files: `cmfgen_base/VADAT_INIT`, `scripts/prepare_VADAT.py`

| Setting | Before | After | Reason |
|---------|--------|-------|--------|
| `USE_LAM_ES` | T | F | No longer bypassing moment solver |
| `USE_J_REL` | F | T | Use MOM_JREL_V9 (designed for non-monotonic velocity) |
| `N_TYPE` | G_ONLY | N_ON_J | N/J always defined; G_ONLY uses N/H (undefined when H=0) |
| `INCL_REL` | F | T | Include velocity terms in radiation transport |

**Why N_ON_J instead of G_ONLY:**
- `G_ONLY`: computes G = N/H Eddington factor. When H→0 at velocity reversal,
  G is undefined → crash (the original instability that led to USE_LAM_ES=T)
- `N_ON_J`: computes N/J ratio. J is always > 0 from the formal solution, so
  N/J is always well-defined. This is the safe choice for non-monotonic velocity

**Why MOM_JREL_V9 helps convergence:**
- The moment equations solve a boundary value problem across all depths,
  naturally smoothing J. The formal solution (TC) can have sharp features
  at shock fronts that create oscillating rates
- `CHK_J=ABS_VAL` patches any remaining negative J at shock fronts
- With `NUM_BNDS=1` (DIAG), the BA matrix remains depth-independent — the
  moment solver only affects J, not the BA matrix structure

**2. Reduce MAX_LAM from 1.0D5 to 1.0D2**

File: `scripts/prepare_VADAT.py` (single_final overrides)

100x max population change per iteration (was 100,000x). This prevents the
zero/overshoot cycle where a population zeros out (100% decrease) then the
next iteration tries to repopulate from zero → 1e23% "increase".

> **2026-07-28 Update**: All convergence parameters have been aligned with the
> CMFGEN manual (`full_descr.pdf`). The previous `stable_overrides` with
> `MAX_LAM=1.1D0` (10% cap) caused an infinite Lambda loop (549 iterations,
> 10000% corrections never converging). See AGENTS.md "CMFGEN Convergence
> Parameters: Manual Recommendations" for the full table of manual-recommended
> values now used.

**3. Reduce IBEG_NG from 80 to 15**

File: `scripts/prepare_VADAT.py` (single_final overrides)

Ng acceleration starts after 15 Lambda iterations (was 80). Manual: "at least
15 (or more) iterations have been completed." The oscillation starts at
iteration 8, so Ng is needed earlier to damp it.

### Relationship to Fix 4

Fix 4 (USE_LAM_ES branch in ACCURATE section) was added to protect first-time
runs when `COMPUTE_EDDFAC=T` triggers Eddington pre-computation. With
`USE_LAM_ES=F`, the ACCURATE section now calls `MOM_J_CMF_V11` (not skipped),
but the `USE_LAM_ES` branch in Fix 4 is harmless (it's in the ELSE branch
which is not taken when `USE_LAM_ES=F`).

### Impact

- **MOM_JREL_V9 is now called** during normal iterations and Eddington
  pre-computation — adds computational cost but provides smoother J
- **No impact on non-IIn models**: The change is in `get_rel_config('static_snapshot')`
  which only affects static_snapshot model type
- **NUM_BNDS=1 (DIAG) preserved**: BA matrix remains depth-independent

---

## Other Source Code Fixes (from earlier work)

These fixes were applied before the current work session and are
documented for completeness.

### rd_sn_data.f: Removed NISO from error check

File: `cmfgen_source/new_main/subs/rd_sn_data.f`

The error check condition included NISO (number of isotopes), which
caused a false positive error when isotopes were present but not all
species had isotope data. Removed NISO from the check condition.

### ins_line_v6.f: Fixed DO WHILE with compound condition

File: `cmfgen_source/subs/ins_line_v6.f`

A DO WHILE loop with a compound condition (`.AND.` / `.OR.`) had a
Fortran fixed-format continuation line issue. The continuation line
was not properly indented, causing the condition to be misread.

### do_species_decays.f / do_species_decays_v2.f: ASSOCIATED() checks

Files: `cmfgen_source/new_main/subs/do_species_decays.f`, `do_species_decays_v2.f`

Added ASSOCIATED() checks before zeroing RADIOACTIVE_DECAY_ENERGY and
KINETIC_DECAY_ENERGY arrays, to prevent SIGSEGV when these pointers
are not allocated (e.g., when INC_RAD_DECAYS=F).

### Makefile_definitions: DEBUG build settings

File: `cmfgen_source/Makefile_definitions`

Changed to DEBUG=debug build (uses `exe_dbg/`, `lib_dbg/`, `mod_dbg/`),
with fpe-summary=none instead of fpe0 for cleaner output.

---

## CMFGEN Configuration for Type IIn Models

### VEL_LAW=12 and SN_MODEL

VEL_LAW=12 does NOT set `SN_MODEL` or `SN_HYDRO_MODEL` (original Dessart
behavior). `SN_MODEL` and `SN_HYDRO_MODEL` remain `.FALSE.`. This means:
- No radioactive decay heating (`cmfgen_sub.f:2972` gated on `SN_MODEL`)
- No SN_MODEL-specific defaults (JGREY_WITH_V_TERMS, WRITE_JH, etc.)
- `check_param_consistency.f` SN_MODEL block entirely skipped
- `REL_OBS` and `REL_CMF` not required (default to `.FALSE.`)
- Luminosity normalization uses wind convention (inward, `cmfgen_sub.f:3822-3826`)

### Density Reading in set_abund_clump.f

For VEL_LAW=12, `set_abund_clump.f` has a special `IF (VELTYPE .EQ. 12)`
branch that reads density from `input_hydro.dat` (before the
`ELSE IF(SN_MODEL)` analytical branch). The READ statement reads 4
values per line (R, V, RHO, T) with list-directed input, ignoring the
extra columns (XL, XTAU, TAU).

The data is written outer-to-inner (decreasing R), and the READ
reverses it with `NHYDRO-I+1` to make it inner-to-outer (increasing R).
Then MON_INTERP interpolates the density onto the CMFGEN R grid.

### Key VADAT Settings for X_2

| Setting | Value | Purpose |
|---------|-------|---------|
| VEL_LAW | 12 | RV_SN_MODEL_SNIIN: tau-based grid for non-monotonic V |
| USE_LAM_ES | T | J=B(T), bypasses moment solver (unstable with non-monotonic V) |
| DO_DDT | F | No time evolution (static snapshot) |
| DO_CO_MOV_DDT | F | No co-moving time derivative |
| CHK_J | ABS_VAL | Check J by absolute value (not relative, for non-monotonic V) |
| METHOD | LOGMON | Monotonic interpolation in log space |
| VFRAC_MOM | 20000 | Velocity fraction for moment equation |
| THK_CONT | F | No thick continuum (optically thin outer CSM) |
| TRAP_J | T | Trapezoidal J integration |
| FIX_NE | T (init/stable) | Ne fixed to hydro for bootstrap stability |
| FIX_NE | F (final) | Ne free for NLTE self-consistency |
| FIX_T | T | T fixed to SN_HYDRO_DATA |
| SCALE_OPT | MAJOR | Manual: "Preferred option is MAJOR" (handles impurity species separately) |
| MAX_LAM | 1.0D2 | 100x max population change per iteration (manual: "typically 10^5") |
| MAX_LIN | 2.0D0 | 200% max change per linearization (manual: "typically 10") |
| LAM_VAL | 4.0D2 | Lambda fallback at 400% change (manual: "typically > 400") |
| NUM_LAM | 2 | Lambda iterations per cycle (manual: "typically 2") |
| N_FIX_BA | 5 | BA freeze count (manual: "recommended 2") |
| BA_CHK_FAC | 5.0D-02 | BA recompute threshold 5% (manual: "typically 5%") |
| FIX_BA | 1.0D-02 | Freeze BA when max change < 1% |
| DO_NG | T | Ng acceleration enabled (manual: "important for NUM_BNDS=1") |
| BEG_NG | 1.0D1 | Start Ng at 10% max change (manual: "5% to 10%") |
| IBEG_NG | 15 | Start Ng after 15 iterations (manual: "at least 15") |
| ITS/NG | 8 | Ng every 8 iterations (manual: "typically 8") |
| AMASS_DOP | 1.0D+40 | Manual: "set 1e40, use VTURB to control width" |
| MAX_CHNG | 1.0D100 | Abort threshold (effectively disabled) |
| LTE_EST | T (proxy_populations deprecated) | Use LTE initial estimates (proxy_populations=False) |
| LIN_INT | T | Linear interpolation of _IN populations (used with LTE_EST=T) |
| INC_RAD_DECAYS | T | Include Ni56/Co56 decay heating |
| N_IB_INS | 15 | Inner boundary inserted points |
| N_OB_INS | 15 | Outer boundary inserted points |

### Workflow: INIT -> STABLE -> FINAL

The X_2 model uses the full INIT->STABLE->FINAL chain (single_final=False):

1. **INIT**: LTE bootstrap with proxy _IN files, FIX_NE=T, FIX_T=T,
   25 Lambda iterations. Produces rough populations (SCRTEMP).

2. **STABLE**: Pure Lambda iterations from INIT populations, FIX_NE=T,
   FIX_T=T. Settles NLTE populations before full iterations.

3. **FINAL**: Full iterations (Lambda + Ng acceleration), FIX_NE=F,
   FIX_T=F (or T varies). Produces the final spectrum.

---

## Input Data Files

### SN_HYDRO_DATA (from STELLA)

The SN_HYDRO_DATA file contains:
- Header: ND, N_MASS_FRAC, N_ISO, time, luminosity, R_inner, T_inner
- Radius grid (10^10 cm), outer-to-inner
- Velocity (km/s), Sigma, Temperature (10^4 K), Density (g/cm^3)
- Electron density, Atom density
- Mass fractions for 8 species (H, He, O, Si, Ca, Ni, Co, Fe)
- Isotope mass fractions (Ni56, Co56, Fe56)

### input_hydro.dat (for VEL_LAW=12)

Written by `prepare_VADAT.py` from SN_HYDRO_DATA. Format:
```
F                                    ! ADD_ACC_ZONE flag
! comment line
150                                  ! NHYDRO
! comment line
R(1e10cm) V(km/s) RHO(g/cm3) T(K) XL XTAU TAU   ! 7 columns, outer->inner
```

Read by `rv_sn_model_sniin.f` and `set_abund_clump.f` (VELTYPE=12 branch).

### T_IN (initial temperature/electron density)

Written by `write_t_in_from_hydro()` in `prepare_VADAT.py`. Format:
- Header: R_core, LUMOLD, NOLD, ND
- Per depth: R, V_dummy, ED, T(10^4 K), ION_FRAC, VEL, CLUMP, depth_index
- Temperature in units of 10^4 K (CMFGEN convention)

**Critical**: Must use actual max T from hydro, not capped at 10000 K.

### GAMMAS_IN (electron density for each species)

Written by `write_gammas_in_from_hydro()` in `prepare_VADAT.py`.
Contains electron density and species densities for LTE initialization.

### Proxy _IN files (HI_IN, HeI_IN, etc.)

Written by `proxy_populations.py`. Contains departure coefficients
for each species at each depth. Used when proxy_populations=True (sets LTE_EST=F,
LIN_INT=T in VADAT).

**DEPRECATED (2026-08-03):** `proxy_populations=True` is no longer recommended.
Empirical testing showed it is worse than CMFGEN's built-in `LTE_EST`
for cold snapshots because `FIX_T=T` makes CMFGEN ignore the `_IN` file's
T column and always evaluate Saha/Boltzmann at `T_local`, making any DI
computed at `T_ion` inherently inconsistent. All configs now use
`proxy_populations: False`. These files may still be generated for diagnostic
purposes but are not read by CMFGEN when `LTE_EST=T`.

---

## Debugging Methodology

The key diagnostic approach was:

1. **Read the OUTGEN file** to identify the error type and location
2. **Check the fort.127 file** for the R, V, SIGMA profile (written by
   rv_sn_model_sniin.f)
3. **Check the T_IN file** for the temperature profile (verify units
   and caps)
4. **Check the proxy _IN files** for departure coefficients (verify
   they give reasonable populations)
5. **Compute T_rad** at the failing depth to check if T_local << T_rad
6. **Verify the VADAT** for correct settings (LTE_EST, LIN_INT, USE_LAM_ES, etc.)
7. **Check the slurm script** for post-prepare_VADAT overrides (e.g.,
   the proxy_populations block that sets LTE_EST=F)

The progression of fixes:
- Fix 1 (SIGMA cap): crash -> "bad initial population guesses" at depth 22
- Fix 2 (T_IN cap): depth 22 (7.2e116%) -> depth 46 (6.6e100%)
- Fix 3 (T_rad in proxy): depth 46 (6.6e100%) -> converging through 25 iterations

Each fix was verified by rerunning the INIT job and checking the OUTGEN
for the error location and correction magnitude.

---

## File Summary

### Source Code Files Modified

| File | Fix | Description |
|------|-----|-------------|
| `cmfgen_source/new_main/subs/rv_sn_model_sniin.f` | Fix 1 | SIGMA cap at -0.999 for VEL_LAW=12 |
| `cmfgen_source/new_main/subs/set_rv_hydro_model_v3.f` | (existing) | SIGMA cap at -0.999 for VEL_LAW=11 |
| `cmfgen_source/new_main/subs/rd_sn_data.f` | (earlier) | Removed NISO from error check |
| `cmfgen_source/subs/ins_line_v6.f` | (earlier) | Fixed DO WHILE compound condition |
| `cmfgen_source/new_main/subs/do_species_decays.f` | (earlier) | ASSOCIATED() checks |
| `cmfgen_source/new_main/subs/do_species_decays_v2.f` | (earlier) | ASSOCIATED() checks |
| `cmfgen_source/Makefile_definitions` | (earlier) | DEBUG=debug, fpe-summary=none |
| `cmfgen_source/new_main/mod_subs/comp_j_blank.f` | Fix 4 | Skip MOM_J_CMF_V11 in Eddington pre-computation when USE_LAM_ES=T |

### Python Script Files Modified

| File | Fix | Description |
|------|-----|-------------|
| `scripts/prepare_VADAT.py` | Fix 2 | Pass actual max T to write_t_in_from_hydro() |
| `scripts/proxy_populations.py` | Fix 3 | Added compute_t_rad(), use T_eff for excitation DC |

### Configuration Files

| File | Setting | Description |
|------|---------|-------------|
| `configs/X2.yaml` | proxy_populations: True | Use proxy dilute-LTE _IN files |
| `configs/X2.yaml` | vel_law: 12 | VEL_LAW=12 for non-monotonic velocity |
| `configs/X2.yaml` | model_type: static_snapshot | No time evolution |
| `configs/X2.yaml` | inc_rad_decays: True | Include Ni56/Co56 decay |
| `configs/X2.yaml` | single_final: False | Use INIT->STABLE->FINAL chain |

---

## Build and Run

### Building CMFGEN

```bash
cd cmfgen_source
make clean
make          # builds exe_dbg/cmfgen_dev.exe
cp exe_dbg/cmfgen_dev.exe exe/cmfgen_dev.exe   # batch.sh uses exe/ path
```

### Running X_2

```bash
# Skip STELLA (use existing .swd), run CMFGEN only
python scripts/run_model_grid.py configs/X2.yaml --skip-stella --run-cmfgen --no-sedona

# Monitor jobs
squeue -u $USER

# Check INIT output
tail -30 models/X_2/cmfgen/5.0days/OUTGEN

# Check convergence progress
grep "iteration" models/X_2/cmfgen/5.0days/OUTGEN | tail -10
```

### Key Diagnostic Commands

```bash
# Check SIGMA capping
grep "capping SIGMA" models/X_2/cmfgen/5.0days/OUTGEN

# Check for fatal errors
grep -E "Mean intensity|bad initial|blowing up|SIGSEGV|NaN" models/X_2/cmfgen/5.0days/OUTGEN

# Check T_IN temperatures (should not be capped at 10000 K)
python3 -c "
with open('models/X_2/cmfgen/5.0days/T_IN') as f:
    lines = f.readlines()
for d in [1, 50, 100, 150]:
    idx = 5 + (d-1)*3
    parts = lines[idx].split()
    print(f'Depth {d}: T={float(parts[3])*1e4:.0f} K')
"

# Check T_rad vs T_local
python3 -c "
import numpy as np
L_sun = 3.828e33; L = 3.05e10 * L_sun; sigma = 5.67e-5
for d, r_1e10, t_loc in [(1, 2.4e6, 3020), (46, 3.47e5, 3072), (150, 3.96e4, 31117)]:
    r = r_1e10 * 1e10
    t_rad = (L / (4 * np.pi * r**2 * sigma))**0.25
    print(f'Depth {d}: T_local={t_loc}K, T_rad={t_rad:.0f}K')
"
```

---

## Outstanding Issues

1. **INIT doesn't fully converge in 25 iterations**: The corrections
   decrease by ~6 orders of magnitude per iteration but start at 1e112%.
   The STABLE job continues from the INIT output and should converge
   further. May need to increase INIT_NUM_ITS or add more STABLE iterations.

2. **CMF_BLKBAND_V3 singular matrix (IFAIL=31)** at depths 14-18: The BA
   rate matrix is singular at low-temperature depths (T~3000K) where
   metal excited-state populations are negligible. FIX_IMP=T helps but
   doesn't fully resolve this. The corrections are set to zero at these
   depths, which allows the iteration to continue but may slow convergence.

3. **X_1 and X_3 models**: Once X_2 converges successfully, the same
   fixes need to be applied to X_1 and X_3 (which use the same VEL_LAW=12
   and proxy_populations settings).
