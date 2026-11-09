# CMF_FLUX Spectral Correctness Audit

Scope: `cmfgen_jfnk/obs/` observed-frame pipeline for non-monotonic velocity
fields. Branch `main`, HEAD `8224a00bb2cc1c8890be6a7e1d5865a3eedea1e7`
(updated 09-Nov-2026 pass; the earlier `0cbcc3b0` reference was stale).
Baseline upstream tree `cmfgen_reference/cur_cmf` (cur_cmf_18jun25).
Executables: `exe/cmf_flux.exe`, `exe/obs_frame.exe` rebuilt 09-Nov 11:30
from this tree (`obs/` make, ifort `-O3`+AVX512); `nm` confirms
`mom_j_cmf_nonmono_v1`, `fg_j_cmf_nonmono_v1`, `obs_frame_sub_v9` linkage.

This document records what was inspected, what is confirmed, what is
suspected, and what remains untested. It distinguishes inspected fact from
hypothesis. Compilation is not evidence of correctness.

---

## 1. Executive summary

The observer-frame engine (`OBS_FRAME_SUB_V9`) is structurally capable of
non-monotonic velocity fields: it builds the ray in the Z coordinate,
inserts points until the projected-velocity change per step is below the
local Doppler scale, maps `nu_obs -> nu_cmf` point by point, and runs a
bounded two-direction frequency-bracket search that survives `T1`
reversals. The Sobolev line render (`DO_SOB_FLUX_OBS`) deposits Doppler
profiles into the stored CMF opacity/emissivity and re-solves the
electron-scattering field once. Flux is integrated over impact parameters
with the standard `mu dmu` weights.

Two confirmed defects were corrected in this pass (both small, isolated,
stock-preserving for monotonic outflow):

- `cmf_flux_sub_v5.f`: the ACCURATE extended grid still divided the
  turbulent velocity by `VEXT(1)`. On a non-monotonic or near-static base
  `VEXT(1)` can be zero, negative, or smaller than an interior `|VEXT|`,
  which drives the turbulent velocity and every Doppler width out of
  range. Now normalized by `max|VEXT|` exactly as the standard `LAW_V1`
  branch already does.
- `obs_frame_sub_v5`/`obs_frame_sub_v9.f`: the non-relativistic CMF-window
  bound used `MAXVAL(VMU_RAY)` (signed max). A strong negative projected
  velocity on a non-monotonic field makes `MINVAL(VMU)` more negative than
  `MAXVAL(VMU)` is positive, under-covering the blue/high-frequency CMF
  window; the per-frequency edge clamps then silently substitute boundary
  opacity. Now uses `MAXVAL(ABS(VMU_RAY))`, identical for monotonic
  outflow.

A second pass (09-Nov-2026) closed the three remaining code-level gaps:

- `COMP_JCONT_V4.INC`: under `DO_NONMONO_SOB` the `USE_LAM_ES=F` moment
  solve ran `MOM_J_CMF_V8` — clamped-sigma, bluer-frequency upwind only —
  while the population run used `MOM_J_CMF_NONMONO_V1`. The OBS `RJ`
  feeding `ESEC*RJ` and `RJ_CMF_ST` was therefore not the population
  field on reversing/`sigma<-1` regions. Now dispatched to
  `MOM_J_CMF_NONMONO_V1` with the population argument list, plus a
  one-sided `dLOG_NU` at `FIRST_FREQ` and the
  `DO_NONMONO_FREQ_MATCH` auto-set mirror.
- `cmf_flux_sub_v5.f`: `ACCURATE=T` (`INC_GRID`) is now rejected when the
  detected velocity field is non-monotonic (`N_NONMONO_OBS>0`), closing
  the `EXTEND_VTSIGMA` `log(V)`/`sigma+1`-floor hole. Monotonic
  `ACCURATE` unchanged; `DO_NONMONO_KNOP`'s `INC_GRID=F` rule unchanged.
- `obs_frame_sub_v9.f`: edge clamps at the three bounded hunts are now
  counted (worst overshoot vs window and vs the full `FREQ_CMF`, worst
  `(ML,LS,ray point)`), reported after the LS loop; out-of-grid clamps
  warn to widen `MIN_CF`/`MAX_CF`.
- `cmf_flux_sub_v5.f`: `OBS_FREQ` is checked finite and strictly
  decreasing after `INS_LINE_OBS_V4`; `CHK_OBS_FLUX` aborts on NaN/Inf
  and censuses negative bins before the `OBSFLUX` and both `OBSFRAME`
  writes.

Remaining items are listed as open risks in section 8. The largest open
question — whether the OBS `RJ` under the new dispatch numerically
matches the population field on a reversing field — is a *runtime*
question now that the code path is consistent; it still requires the
model A/B test in TODO T4/T6.

---

## 2. Files and routines investigated

| File | Role |
|---|---|
| `obs/cmf_flux_v5.f` | entry; reads RVTJ, populations, oscillator data, controls |
| `obs/cmf_flux_sub_v5.f` | driver: grids, opacities, J solves, Sobolev deposit, output |
| `obs/obs_frame_sub_v9.f` | observer-frame formal integral (the spectrum former) |
| `obs/ins_line_v5.f` | CMF frequency grid (continuum + BLANK resonance zones) |
| `obs/ins_line_obs_v4.f` | observer-frame frequency grid |
| `obs/COMP_JCONT_V4.INC` | CMF continuum J-solve dispatch (formal vs moment) |
| `obs/rd_cmf_flux_controls.f` | control parsing and fail-closed consistency checks |
| `obs/check_cmf_flux_param_consis.f` | parameter consistency |
| `obs/mod_cmf_obs.f`, `obs/mod_freq_obs.f` | shared storage (SIGMA_TRUE etc.) |
| `subs/extend_vtsigma.f` | extended-grid V/T/SIGMA interpolation |
| `plane/fg_j_cmf_v12.f`, `fg_j_cmf_v13.f`, `mom_j_cmf_v8.f`, `mom_j_cmf_nonmono_v1.f`, `fg_j_cmf_nonmono_v1.f` | J solvers (dispatch comparison) |
| `new_main/mod_subs/comp_j_blank.f` | population-solve J dispatch (for consistency comparison) |

Diffed `obs_frame_sub_v9.f` against `cmfgen_reference/cur_cmf` to separate
fork changes from upstream.

---

## 3. Execution path (Phase 1)

`CMF_FLUX_V5` -> reads RVTJ + POP + oscillator/line-link/photoionization
data -> `RD_OPTIONS_INTO_STORE`/`RD_CMF_FLUX_CONTROLS` ->
`CMF_FLUX_SUB_V5(ND,NC,NP,...)`.

Inside `CMF_FLUX_SUB_V5`:

1. `RD_CMF_FLUX_CONTROLS` parses flags; `CHECK_CMF_FLUX_PARAM_CONSIS`.
2. **Non-monotonic detector** (`cmf_flux_sub_v5.f:591-624`): counts
   `V(I+1) > V(I)*(1+1e-8)` and `SIGMA(I) < -1-1e-6`. If flagged and none
   of `DO_NONMONO_SOB`/`DO_NONMONO_KNOP`/`ALLOW_NONMONO_OBS` is set it
   `STOP`s. This is the fail-closed gate.
3. Under the nonmono flags, `SIGMA_TRUE_CORE` builds the unclamped
   `dlnV/dlnR-1`; undefined cells abort unless `ALLOW_SIGMA_UNDEF`.
4. `VINF = MAXVAL(ABS(V))` on nonmono paths (was `V(1)`); stock keeps
   `VINF=V(1)`.
5. `VTURB_VEC`/`VDOP_VEC` from `TURB_LAW` (`LAW_V1` fixed to `max|V|`
   normalization; `LAW_TAU1` cumulative tau_es).
6. Continuum frequencies -> `INS_LINE_V5` merges BLANK-line resonance
   zones -> `INS_LINE_OBS_V4` builds `OBS_FREQ` (uses `SOB_FREQ_IN_OBS`
   as `INCL_ALL_LINES` so SOB lines get observer-grid refinement).
7. `IMPAR` impact parameters + `NORDANGQW`/`GENANGQW` angular weights.
8. Optional `ACCURATE` extended grid: `REXT_COEF_V2`, `EXTEND_VTSIGMA`,
   `IMPAR`/`NORDANGQW` on the extended mesh.
9. Continuum loop over `NU`: `COMP_JCONT_V4.INC` solves CMF `RJ`/`FEDD`
   per frequency; `OPACITIES` + `ADD_LINE`/`SET_PROF_V6` build `ETA/CHI`;
   Sobolev EW estimators (`SOBEW_*`) fill line EWs; `DO_SOB_FLUX_OBS`
   deposits `ADD_SOB_LINE_PROF` into `ETA_CMF_ST`/`CHI_CMF_ST`.
10. `OBS_FRAME_SUB_V9` called (phase 1 continuum/blanketing), then under
    `DO_SOB_FLUX_OBS` a second electron-scattering solve and a second
    `OBS_FRAME_SUB_V9` (phase 2) update the resolved line profile.
11. Writes `OBSFLUX`, `OBSFRAME`, `IP_OBS`, `RTAU/ZTAU` diagnostics.

Non-monotonic treatment enters at: the detector (2), `SIGMA_TRUE` (3),
`VINF` (4), the `COMP_JCONT` dispatch (`DO_NONMONO_KNOP` ->
`FG_J_CMF_NONMONO_V1`; `DO_NONMONO_SOB` keeps clamped-sigma `V12`/`V8`),
the Sobolev EW estimators (9), the `ADD_SOB_LINE_PROF` deposit (9), and
the `OBS_FRAME_SUB_V9` ray/bracket logic (10).

Executable: `exe/cmf_flux.exe` (built from `obs/*.o` + `libsubs/libpgplt/
libdisp/libmod`). `obs_frame.exe` is a separate driver sharing
`obs_frame_sub_v9.o`.

---

## 4. Defects found

### 4.1 CONFIRMED + FIXED — extended-grid turbulent velocity (HIGH)

`cmf_flux_sub_v5.f` (was line ~1257):
`TA=VTURB_MIN+(VTURB_MAX-VTURB_MIN)*VEXT/VEXT(1)`.

The standard `LAW_V1` branch was already corrected to
`VTURB_MIN+(VTURB_MAX-VTURB_MIN)*MIN(1,|V|/max|V|)` with a `max|V|>0` floor
("bug fix, not flag-gated"). The extended-grid branch kept the old
`VEXT/VEXT(1)` form. `VEXT` comes from `EXTEND_VTSIGMA`, a **log-plane**
interpolation (`VEXT=exp(sum coef*log V)`), so on a non-monotonic field
`VEXT(1)` can be zero, negative, or smaller than an interior `|VEXT|`;
the ratio then leaves `[0,1]` and corrupts `VDOP_VEC_EXT`, i.e. every
Doppler width and line profile on the extended grid.

Even for monotonic outflow this is only safe because `VEXT(1)` happens to
be the largest velocity; for an IIn-like slow outer boundary it is not.

**Fix applied**: replaced with the identical `max|VEXT|` normalization and
`MIN(1,|VEXT|/max|VEXT|)` clamp plus the `max|VEXT|>1e-12` floor. For
monotonic outflow `VEXT(1)=max|VEXT|` so the result is identical to the
old expression; stock `ACCURATE` behaviour preserved.

Severity before fix: **HIGH** (silently wrong Doppler widths -> wrong
line profiles -> wrong spectrum, exactly the failure mode the audit was
asked to find).

### 4.2 CONFIRMED + FIXED — non-rel CMF window under-coverage (MEDIUM)

`obs_frame_sub_v9.f` (~line 665): `MAX_VMU=MAXVAL(VMU_RAY)` fed both
window bounds `1 +/- 1.2*MAX_VMU` in the `DO_REL_CORRECTIONS=F` branch.

The Doppler reach of `nu_obs` is set by `max|VMU|`. The relativistic
branch already uses the exact extremes `MAXVAL/MINVAL((1-VMU)*GAM)`. In
the non-rel branch the signed max misses a large negative projected
velocity (negative-V excursion on a non-monotonic field, or an
asymmetric ray), so the high-frequency CMF window `[STRT_INDX_CMF,
END_INDX_CMF]` is too narrow. `T1` then lands above `SM_FREQ_CMF(1)` and
the bounded clamp substitutes boundary opacity — a silent spectral
error, not an abort.

**Fix applied**: `MAX_VMU=MAXVAL(ABS(VMU_RAY))`. For monotonic outflow
the mirrored far-side ray makes `|VMU|` symmetric, so `MAXVAL(ABS)` ==
`MAXVAL(VMU)` and the window is unchanged. Only sign-changing fields get
the wider (correct) window.

Severity before fix: **MEDIUM** (only the `DO_REL_CORRECTIONS=F` path;
requires an asymmetric `VMU` excursion; produces a localized edge
artifact rather than a global failure).

### 4.3 OPEN — `EXTEND_VTSIGMA` log-V interpolation (MEDIUM, stock)

`subs/extend_vtsigma.f` interpolates `V`, `T`, `SIGMA` as
`exp(sum coef*log x)`. `LOG(V)` is undefined for `V<=0` (infall,
reversal through zero, or a zero base velocity) -> NaN `VEXT` ->
poisoned extended grid. For `V>0` non-monotonic fields the geometric
interpolation is an approximation that need not preserve a sharp
reversal. `SIGMAEXT` floors `SIGMA+1` at `1e-4`, clamping `sigma<-1`
cells to `~-0.9999` (consistent with the clamped-sigma SOB contract, but
it does erase the compression signature on the extended mesh).

This is **stock upstream** behaviour; changing the interpolation would
alter every monotonic `ACCURATE` result, so it was not changed. The
correct handling for non-monotonic fields is to reject `ACCURATE` for
that configuration — **now implemented**: `CMF_FLUX_SUB_V5` stops when
`ACCURATE .AND. N_NONMONO_OBS>0`, i.e. on the detected field (any
reversal or `sigma<-1` cell), under every flag that let the detector
pass. Monotonic `ACCURATE=T` is unaffected (see TODO T5).

### 4.4 RESOLVED — extended-grid solver assumes monotonic sigma (MEDIUM)

`COMP_JCONT_V4.INC` extended (`ACCURATE`) branch dispatches to
`FG_J_CMF_V12`/`MOM_J_CMF_V8` (or the relativistic variants) with
`SIGMAEXT`/`VEXT`. `DO_NONMONO_KNOP` is explicitly rejected here
(`INC_GRID=F` enforced in `RD_CMF_FLUX_CONTROLS`), but `DO_NONMONO_SOB`
is **not** rejected for `ACCURATE=T` at the key level. Under SOB the
clamped `sigma>=-1` is the intended contract (same as the population
continuum), so the solver choice is consistent — but the extended grid
still inherits the log-V interpolation limit of 4.3 and the extended
solvers carry no non-monotonic treatment. `DO_NONMONO_SOB+ACCURATE` on a
genuinely non-monotonic field is therefore an **unsafe configuration**,
now closed: the `N_NONMONO_OBS>0` rejection in `CMF_FLUX_SUB_V5` stops
the combination regardless of which nonmono flag is set, while leaving
`ACCURATE` free on monotonic fields (see TODO T5).

### 4.5 RESOLVED — OBS J-solve vs population J-solve consistency (MEDIUM)

Population dispatch (`comp_j_blank.f`): `DO_NONMONO_SOB+USE_LAM_ES=T` ->
`FG_J_CMF_V13`; `+USE_LAM_ES=F` -> `MOM_J_CMF_NONMONO_V1`;
`DO_NONMONO_KNOP` -> `FG_J_CMF_NONMONO_V1`.

OBS dispatch (`COMP_JCONT_V4.INC`): non-extended branch uses
`FG_J_CMF_V12` (`USE_LAM_ES=T`), `MOM_J_CMF_V8` (`USE_LAM_ES=F`), and
`FG_J_CMF_NONMONO_V1` only under `DO_NONMONO_KNOP`.

- `V13` vs `V12`: the formal solution is identical; `V13` only
  *additionally stores* the ALO diagonal for the BA solver, which OBS
  does not consume. **Not a defect.**
- `MOM_J_CMF_V8` vs `MOM_J_CMF_NONMONO_V1` under `DO_NONMONO_SOB+
  USE_LAM_ES=F`: **confirmed defect, now fixed.** `V8` computes midpoint
  coefficients `W=GAMH*(1+AV_SIGMA*NMID_ON_HMID)`,
  `EPS=GAMH*AV_SIGMA*NMID_ON_J/(1+W)`, `PSI=DTAUONQ*GAM*(1+SIGMA*K_ON_J)`
  and couples every cell to the *bluer* frequency only; at `sigma<-1`
  the `1+sigma*K` denominator degenerates and the upwind direction is
  wrong for compressing cells. Its revised-grid insertion
  `IT1=INT((V_SM(I)-V_SM(I+1))/T2)` also goes negative through a
  reversal, so no sub-grid points are added where they are most needed.
  The dispatch now calls `MOM_J_CMF_NONMONO_V1` under `DO_NONMONO_SOB`
  in the `USE_LAM_ES=F` leg, with the same argument list as the
  population call (`SIGMA,SIGMA_TRUE`, `INNER_BND_METH`,
  `IB_STAB_FACTOR`, `H_CHK_OPTION`, `FREQ_INDX`, `NCF`), preceded by
  `SET_MOM_NU_NEXT(NU(MIN(ML+1,NCF)))`. Ray moments still come from the
  `FG_J_CMF_V12` formal solve preceding the moment call, mirroring the
  population V13 pairing. Stock `V8` is retained when
  `DO_NONMONO_SOB=F`. Two companion pieces were required: a one-sided
  `dLOG_NU` at `FIRST_FREQ` (the nonmono solver consumes `dLOG_NU` in
  `chi_hat` even under `INIT`), and the `DO_NONMONO_FREQ_MATCH`
  auto-set/forwarding (`SET_MOM_FREQ_MATCH`, `SET_FG_FREQ_MATCH`)
  mirroring the VADAT enforcement — without it the OBS grid reads the
  lagged `J_NEXT` off the characteristic foot. Lagged-field state
  (`J_NEXT_FRQ`/`H_NEXT_FRQ`) cold-starts guarded on pass 1 and refills
  across the ES iterations; a `MOMJNXT_GS` checkpoint in the run
  directory is also restored when present. **Remaining: numerical A/B
  on a converged non-monotonic SOB model (TODO T4).**

### 4.6 OPEN — Sobolev line deposit coverage on the CMF grid (MEDIUM)

`INS_LINE_V5` only inserts resonance-zone refinement for
`TRANS_TYPE='BLA'` (BLANK) lines. Under `GLOBAL_LINE=SOB` **no** line gets
CMF resonance-zone refinement; the `NU` grid near a line is the
continuum grid. `ADD_SOB_LINE_PROF` deposits the Doppler/Voigt profile
onto whatever CMF cells bracket `F0` within `±W`, normalized by cell
width. If the continuum spacing near `F0` is coarser than the thermal
Doppler width, the deposited opacity is concentrated in 1-2 cells. The
observer-frame formal integral then has to resolve the line purely
through the `nu_obs -> nu_cmf` mapping and ray sampling. This is the
intended Sobolev-render design, but its accuracy depends on (a) the
observer grid (`INS_LINE_OBS_V4` does refine around SOB lines when
`SOB_FREQ_IN_OBS=T`) and (b) ray resolution. Not independently verified
here — flagged for the validation suite.

### 4.7 VERIFIED (structure) — two-phase `DO_SOB_FLUX_OBS` bookkeeping (LOW-MEDIUM)

Phase 2 re-solves electron scattering after line deposition and re-runs
`OBS_FRAME_SUB_V9`. Source trace confirms the bookkeeping:

- `NSOB_PHASE=2` only when `DO_SOB_FLUX_OBS=T`; the EW/deposit loop runs
  in phase 1 only (`IF(SOB_PHASE .EQ. 1)`), so `ADD_SOB_LINE_PROF`
  deposits each line exactly once into `ETA_CMF_ST`/`CHI_CMF_ST`.
- Phase 2 takes the `ES_COUNTER>1` reload branch per frequency
  (`SOB_PHASE>1` is an explicit term), so the deposited opacity is read
  back into `ETA`/`CHI` rather than recomputed; the `ES_COUNTER.EQ.1`
  store-back is idempotent (values reloaded, not regenerated).
- `RJ_CMF_ST` is updated unconditionally per frequency, so the phase-2
  `J` carries the line photons; the render fold-in applies `CLUMP_FAC`
  and the `ESEC*RJ`(+dust) convolution on the `ETA_RND_ST`/`CHI_RND_ST`
  copies, keeping the stored arrays as the phase-2 solve input.
- Phase 2 forces `COHERENT_ES=.TRUE.` and `COMPUTE_EDDFAC=.TRUE.` so the
  `ES_J_CONV`/`EDDFACTOR` record state built for a single pass cannot be
  misread on re-entry; `ACCESS_F` is restarted per phase.
- `OBS_FLUX` is zeroed inside `OBS_FRAME_SUB_V9`, so the phase-2 call
  rebuilds it with line profiles; `OBSFRAME` is overwritten by the
  phase-2 write (same `N_OBS` records, no stale tail). `OBSFLUX` is
  written in phase 1 only = continuum+BLANK. `OBSFRAME` (final) =
  continuum+BLANK+SOB profiles.
- Maser accounting preserved: `SOB_MASER_DROP` (cell count),
  `SOB_DROP_EMIS`/`SOB_INJ_EMIS` (emissivity sums); a >10% drop prints
  the lower-bound warning. `NEG_OPACITY` blend-group cells are dropped
  on the same footing as per-line masers.

Structure verified against the requirements; numerical confirmation
(no double counting, line flux in OBSFRAME) still needs the T6/T8
runtime tests.

### 4.8 LOW — ray-resolution heuristic (informational)

`NINS=|VMU(I-1)-VMU(I)|/MAX_DEL_V_RES_ZONE(I-1)+1` samples the projected
velocity at the `ND` depth points. It resolves reversals that occur
*between* depth points only to the extent the depth grid resolves `V(R)`.
Because `V` is only defined on the depth grid and interpolated, extrema
sit on grid points; `|dVMU|` is a reasonable local measure. `NI_MAX` now
uses a total-variation bound `max|V|+sum|dV|` with a runtime overflow
check that `STOP`s (fail-closed). No defect proven; convergence must be
demonstrated numerically (Test J).

### 4.9 LOW — flux integration (inspected, stock)

`OBS_FLUX(ML)+=HQW_AT_RMAX(LS)*PAR_FLUX` over rays, then
`*6.59934*R(1)^2` for Jy at 1 kpc. `HQW_AT_RMAX` is the `mu dmu`
quadrature weight; the `PdP ~ R^2 mu dmu` geometry is the standard
spherical flux. Inner-boundary diffusion injection
(`IB_J_CMF`/`IB_H_CMF` -> `I=J+3 mu H`) is applied when
`(DO_NONMONO_SOB .OR. DO_NONMONO_KNOP) .AND. DIF`, and `HOLLOW_CORE` is
forced false so core rays terminate at `R(ND)` and receive the injected
flux. Consistent with the stated design; luminosity conservation not yet
verified numerically.

---

## 5. Physical explanation of the corrected defects

The observer-frame transfer is a formal integral along each ray. At every
ray point the code needs the CMF opacity/emissivity at the *local*
comoving frequency `nu_cmf = nu_obs (1 - mu v/c) gamma`. Two things must
be right: (a) the `[STRT_INDX_CMF,END_INDX_CMF]` window must contain every
`nu_cmf` the ray will ask for, and (b) the Doppler widths used to build
those opacity/emissivity profiles must be correct.

Defect 4.2 shrank the window on the blue side for sign-changing `VMU`,
so `nu_cmf` fell off the grid and the code substituted edge opacity —
a smooth but wrong spectrum. Defect 4.1 corrupted the Doppler widths on
the extended grid, so line profiles had wrong widths/depths. Both are
exactly the "converged model, wrong spectrum" failure mode.

---

## 6. Supporting evidence

- Code inspection with line references above; fork-vs-upstream diff of
  `obs_frame_sub_v9.f` catalogued.
- `make` in `obs/` rebuilds `cmf_flux_sub_v5.o`, `obs_frame_sub_v9.o`,
  `exe/cmf_flux.exe`, `exe/obs_frame.exe` cleanly after both edits
  (ifort, only pre-existing format remarks).
- The `LAW_V1` standard path already contains the identical
  `max|V|` normalization, establishing the intended behaviour the
  extended path lacked.

No numerical spectrum comparison has been run yet; see section 8.

---

## 7. Proposed corrections (status)

| ID | Correction | Status |
|---|---|---|
| C1 | Extended-grid `VTURB` -> `max|VEXT|` normalization | DONE (prior pass) |
| C2 | Non-rel CMF window -> `MAXVAL(ABS(VMU_RAY))` | DONE (prior pass) |
| C3 | Reject `ACCURATE` for genuinely non-monotonic `V` (fail-closed on the detected field `N_NONMONO_OBS>0`), matching the `DO_NONMONO_KNOP` `INC_GRID=F` precedent | DONE (code; runtime matrix pending) |
| C4 | Route OBS `USE_LAM_ES=F` moment solve under `DO_NONMONO_SOB` to `MOM_J_CMF_NONMONO_V1` (with `SET_MOM_NU_NEXT`, one-sided `dLOG_NU` at `FIRST_FREQ`, `DO_NONMONO_FREQ_MATCH` mirror) | DONE (code; numerical A/B pending) |
| C5 | OBS_FRAME edge-clamp census (count, worst overshoot vs window and vs grid, location) | DONE (code) |
| C6 | Output sanity: `OBS_FREQ` strictly-decreasing abort; `CHK_OBS_FLUX` NaN abort + negative census at `OBSFLUX`/`OBSFRAME` writes | DONE (code) |
| C7 | Optional `CMF_FLUX_AUDIT` diagnostic output | TODO T7 |
| C8 | Validation suite A-J | TODO T6 |

All new code paths are pinned by
`unit_tests/python/test_obs_nonmono_dispatch.py` (37 source-level
checks) and `unit_tests/python/test_sob_prof.py`
(`ADD_SOB_LINE_PROF` driver tests). The production `obs/` build was
re-verified (`make all`, ifort `-O3` + AVX512); `exe/cmf_flux.exe`
relinked with `mom_j_cmf_nonmono_v1` linked and the new diagnostics
embedded.

---

## 8. Remaining uncertainties / open risks

1. **OBS `RJ` numerical consistency (4.5)** — the dispatch fix is
   source-verified, but the population-vs-OBS `J` comparison on a
   reversing field has not been run. Until that A/B executes, treat
   `DO_NONMONO_SOB+USE_LAM_ES=F` OBS spectra as corrected-in-code,
   unverified-in-numerics.
2. **Sobolev deposit fidelity (4.6)** — whether the deposited line
   opacity reproduces the correct Sobolev optical depth through the
   observer-frame formal integral on a coarse continuum grid.
3. **`ACCURATE`+nonmono validity (4.3/4.4)** — now fail-closed on the
   detected field; the extended-grid path remains untested for
   monotonic-code-limit equivalence and needs the T5 runtime matrix
   before `INC_GRID=T` is trusted anywhere.
4. **Ray convergence (4.8)** — `NINS`/`NI_MAX` adequacy on multi-reversal
   fields must be shown by grid-refinement tests (`N_INS_OBS`=3/5/7),
   not assumed.
5. **Frequency coverage** — the `VINF=MAXVAL(ABS(V))` fix widened the
   margins, but `NU_MAX_OBS/NU_MIN_OBS` still use `2*VINF`; extreme
   relativistic or slow-CSM cases may still under-cover. The new
   edge-clamp census will surface any residual gap in the log.
6. **Executable provenance** — `exe/cmf_flux.exe` rebuilt 09-Nov 11:30
   from this tree; run scripts invoke `exe/cmf_flux.exe`. Record
   `git rev-parse HEAD` in run logs when validating.
7. **Validation suite (Phase 7)** — none of tests A-J has been run on a
   live model; no spectrum has been compared against a reference. All
   correctness claims above are code-level, not numerical. In
   particular `MOM_J_CMF_NONMONO_V1` under `USE_LAM_ES=F` is flagged
   FAIL in the population run history (job 12316397, non-contracting
   lagged map on m15) — the OBS call inherits the same solver; whether
   it converges in a single-sweep OBS context is untested.

---

## 9. Regression-safety statement

The earlier corrections (C1/C2) reduce exactly to the previous
expressions for a monotonic outflow. The C3/C4/C5/C6 changes are
flag- or field-gated:

- C3 fires only when the detector already found a reversal/`sigma<-1`
  cell — a configuration that previously produced silent garbage.
- C4 dispatches only under `DO_NONMONO_SOB` in the `USE_LAM_ES=F` leg;
  `MOM_J_CMF_V8` remains the call for `DO_NONMONO_SOB=F` (stock
  unchanged, including the extended branch).
- `DO_NONMONO_FREQ_MATCH` defaults F in the obs controls and is only
  forced T under the same flags that force it in the population
  reader; explicit F is overridden with a printed banner, matching
  VADAT semantics.
- C5/C6 add counters and checks; they do not alter any computed
  quantity (the abort paths replace what were previously silent
  substitutions or crashes).
- One new control key was added (`DO_NONMONO_FREQ_MATCH`, obs-side
  mirror of the existing VADAT key); defaults preserve prior
  behaviour.
