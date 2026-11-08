# CMF_FLUX Spectral Correctness Audit

Scope: `cmfgen_jfnk/obs/` observed-frame pipeline for non-monotonic velocity
fields. Branch `main`, reference commit `0cbcc3b0`. Baseline upstream tree
`cmfgen_reference/cur_cmf` (cur_cmf_18jun25).

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

Remaining items are listed as open risks in section 8. The largest open
question is not in `OBS_FRAME_SUB_V9` itself but in whether the
pre-observer inputs (CMF `ETA/CHI`, `RJ`, the Sobolev-deposited line
opacity) are correct for a non-monotonic field; several of those inputs
are produced by solvers that assume `sigma >= -1`.

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
that configuration (see TODO T3) rather than silently interpolate.

### 4.4 OPEN — extended-grid solver assumes monotonic sigma (MEDIUM)

`COMP_JCONT_V4.INC` extended (`ACCURATE`) branch dispatches to
`FG_J_CMF_V12`/`MOM_J_CMF_V8` (or the relativistic variants) with
`SIGMAEXT`/`VEXT`. `DO_NONMONO_KNOP` is explicitly rejected here
(`INC_GRID=F` enforced in `RD_CMF_FLUX_CONTROLS`), but `DO_NONMONO_SOB`
is **not** rejected for `ACCURATE=T`. Under SOB the clamped `sigma>=-1`
is the intended contract (same as the population continuum), so the
solver choice is consistent — but the extended grid still inherits the
log-V interpolation limit of 4.3. Whether `DO_NONMONO_SOB+ACCURATE` is a
supported/tested cell is unresolved; treat as untested until a
convergence test exists.

### 4.5 OPEN — OBS J-solve vs population J-solve consistency (MEDIUM)

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
  USE_LAM_ES=F`: OBS uses `V8`, the population solve uses the sign-split
  non-monotonic moment solver. If `V8` assumes monotonic `sigma`, the OBS
  `RJ` (which feeds the electron-scattering source `ESEC*RJ` and the
  stored `RJ_CMF_ST`) can differ from the population `J`, making the
  emergent flux inconsistent with the populations. **This is the most
  important open consistency question.** Needs a direct comparison test
  (TODO T4).

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

### 4.7 OPEN — two-phase `DO_SOB_FLUX_OBS` bookkeeping (LOW-MEDIUM)

Phase 2 re-solves electron scattering after line deposition and re-runs
`OBS_FRAME_SUB_V9`. Comments state `OBSFLUX` remains the
continuum+BLANK product while `OBSFRAME`/the resolved-profile path is
updated. Need to confirm no double-counting of line emissivity and that
dropped masing cells (`SOB_MASER_DROP`, `SOB_DROP_EMIS`) are only a lower
bound, per the existing warnings. Diagnostics exist (`SOB_*` counters)
but were not exercised in this pass.

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
| C1 | Extended-grid `VTURB` -> `max|VEXT|` normalization | DONE (this pass) |
| C2 | Non-rel CMF window -> `MAXVAL(ABS(VMU_RAY))` | DONE (this pass) |
| C3 | Reject `ACCURATE` for genuinely non-monotonic `V` (fail-closed), matching the `DO_NONMONO_KNOP` `INC_GRID=F` precedent | TODO T3 |
| C4 | Resolve OBS-vs-population `RJ` solver mismatch under `DO_NONMONO_SOB+USE_LAM_ES=F` (route OBS to `MOM_J_CMF_NONMONO_V1` or document divergence) | TODO T4 |
| C5 | Optional `CMF_FLUX_AUDIT` diagnostic output | TODO T5 |
| C6 | Validation suite A-J | TODO T6 |

---

## 8. Remaining uncertainties / open risks

1. **OBS `RJ` consistency (4.5)** — the biggest open physics question.
   If `MOM_J_CMF_V8` gives a different `J` than `MOM_J_CMF_NONMONO_V1`
   on a reversing field, the electron-scattering source and stored `RJ`
   are inconsistent with the populations. Needs a numerical A/B test.
2. **Sobolev deposit fidelity (4.6)** — whether the deposited line
   opacity reproduces the correct Sobolev optical depth through the
   observer-frame formal integral on a coarse continuum grid.
3. **`ACCURATE`+nonmono validity (4.3/4.4)** — log-V interpolation and
   monotonic-sigma solvers; should be rejected explicitly until tested.
4. **Ray convergence (4.8)** — `NINS`/`NI_MAX` adequacy on multi-reversal
   fields must be shown by grid-refinement tests, not assumed.
5. **Frequency coverage** — the `VINF=MAXVAL(ABS(V))` fix widened the
   margins, but `NU_MAX_OBS/NU_MIN_OBS` still use `2*VINF`; extreme
   relativistic or slow-CSM cases may still under-cover.
6. **Executable provenance** — `exe/cmf_flux.exe` is the built artifact;
   confirm the run scripts invoke it (not a stale copy) before trusting
   any comparison.
7. **Validation suite (Phase 7)** — none of tests A-J has been run; no
   spectrum has been compared against a reference. All correctness
   claims above are code-level, not numerical.

---

## 9. Regression-safety statement

Both applied changes reduce exactly to the previous expressions for a
monotonic outflow (`VEXT(1)=max|VEXT|`; `|VMU|` symmetric), so the stock
`ACCURATE=F`/`ACCURATE=T` monotonic workflow is numerically unchanged.
No new controls were added; no defaults changed; no broad rewrite.
