# CMF_FLUX_OPTIMAL_INPUT_SETUP

Recommended `CMF_FLUX_PARAM` configurations for the three model classes:
Knop (`DO_NONMONO_KNOP`), Sobolev non-monotonic (`DO_NONMONO_SOB`), and
monotonic reference. Companion to `CMF_FLUX_SPECTRAL_AUDIT.md` and
`CMF_FLUX_TODO.md`. HEAD `8224a00bb2cc1c8890be6a7e1d5865a3eedea1e7`.

Status convention per setting: **[mandatory]** = required for the
configuration to be valid; **[physics]** = changes the physical model,
must mirror the population run; **[numerical]** = resolution/control only;
**[unverified]** = plausible default, no supporting test yet.

Nothing below is labelled optimal: these are *recommended starting
configurations*. Convergence tests (T6 suite) have not been run on a live
model at this revision.

---

## 0. Configuration traps found in this tree

- `NON_MONO_SOLVER` is a **dead key**. The existing
  `CMF_FLUX_PARAM_INIT` templates under `cmfgen_validation/` set
  `NON_MONO_SOLVER=T`, but `rd_cmf_flux_controls.f` reads only
  `DO_NONMONO_KNOP` / `DO_NONMONO_SOB`. An optional-key miss is silent,
  so a KNOP template carrying `NON_MONO_SOLVER` runs the *stock* V12/V8
  path — the non-mono detector would then abort the run on a non-mono
  RVTJ (fail-closed) or, with `ALLOW_NONMONO_OBS=T`, render a wrong
  spectrum. Always write `DO_NONMONO_KNOP` / `DO_NONMONO_SOB`.
- Cached `EDDFACTOR`, `ES_J_CONV`, and observer-grid files must be
  deleted when `DO_SOB_FLUX_OBS`, `SOB_FREQ_IN_OBS`, or the frequency
  controls change; the run script (`batobs.sh`) already deletes
  `ES_J_CONV`. A stale grid silently omits SOB line frequencies.
- `INC_GRID=T` now aborts whenever the *detected* velocity field is
  non-monotonic (`N_NONMONO_OBS>0`), under any flag — including the
  `ALLOW_NONMONO_OBS` escape route. This is by design (see audit 4.3/4.4).
- `SOB_SIGMA_MAX>0` clamps `SIGMA_TRUE` *for the OBS moment solve too*
  since the 09-Nov-2026 dispatch change (the solver now reads
  `SIGMA_TRUE` directly, not only the EW estimator). Leave it `<=0`
  (off) to keep OBS consistent with the unclamped population field.
- `ALLOW_SIGMA_UNDEF` and `ALLOW_NONMONO_OBS` downgrade fail-closed
  aborts to warnings. Keep both `F` in production; they exist for
  diagnosis, not spectra.

---

## 1. KNOP population models (`DO_NONMONO_KNOP=T` in VADAT)

```
T    [DO_NONMONO_KNOP]   ! [mandatory] selects FG_J_CMF_NONMONO_V1 for CMF J
F    [DO_NONMONO_SOB]    ! [mandatory] Knop populations carry no SOB treatment
T    [USE_LAM_ES]        ! [mandatory] RJ from the formal solve; USE_LAM_ES=F is rejected
F    [INC_GRID]          ! [mandatory] extended-grid Knop solve not implemented
BLANK [GLOBAL_LINE]      ! [mandatory] every line must be BLANK
F    [LAM_SET]           ! [mandatory] LAM_SET would give edge lines clamped-sigma SOB
F    [USE_FIXED_J]       ! [mandatory] rejected: bypasses the Knop solver
F    [DO_CMF_EW]         ! [mandatory] CMF EW estimator uses clamped SIGMA
F    [DO_SOB_LINES]      ! BLANK mode has no SOB line list to deposit
F    [SOB_FREQ_IN_OBS]   ! no SOB frequencies to add under BLANK
F    [DO_SOB_FLUX_OBS]   ! requires DO_NONMONO_SOB
```

All `TRANS_*` entries must be `BLANK` (or `GLOBAL_LINE=NONE` with every
present ion `BLANK`); the reader enforces this.

Also rejected under KNOP: `USE_FORMAL_REL`, `USE_J_REL`, `USE_DJDT_RTE`,
`PLANE_PARALLEL`, `PLANE_PARALLEL_NO_V` — the reader stops on each.

`DO_NONMONO_FREQ_MATCH` is auto-set `T` under `DO_NONMONO_KNOP`
(absent key -> T, explicit F overridden with a banner), because the OBS
frequency grid differs in `NCF` and spacing from the model grid and the
lagged `I_NEXT`/`J_NEXT` upwind foot lands between bins
(Ceulemans et al. 2025). **[physics-consistent solver state]**

---

## 2. SOB population models (`DO_NONMONO_SOB=T` in VADAT)

```
F    [DO_NONMONO_KNOP]   ! unless intentionally combining Knop J + SOB lines
T    [DO_NONMONO_SOB]    ! [mandatory] unclamped SIGMA_TRUE in EW + moment solve
SOB  [GLOBAL_LINE]       ! [mandatory] (or NONE with every TRANS_*=SOB)
T    [DO_SOB_LINES]      ! [mandatory for profiles] line list for EW + deposits
T    [SOB_FREQ_IN_OBS]   ! [mandatory for profiles] adds SOB freqs to obs grid
T    [DO_SOB_FLUX_OBS]   ! [mandatory for resolved SOB profiles]
F    [INC_GRID]          ! [mandatory] aborts on non-monotonic field anyway
F    [DO_CMF_EW]         ! [mandatory] clamped-sigma EW invalid for sigma<-1
F    [ALLOW_NONMONO_OBS] ! [mandatory] fail-closed policy
F    [ALLOW_SIGMA_UNDEF] ! [mandatory] abort on undefined SIGMA_TRUE cells
<=0  [SOB_SIGMA_MAX]     ! off: keep OBS on the population SIGMA_TRUE (see traps)
```

### 2a. Electron-scattering branch — two valid sub-modes

`USE_LAM_ES` mirrors whichever the population run used:

- `USE_LAM_ES=T`: `RJ` comes from `FG_J_CMF_V12` (the OBS partner of the
  population's `FG_J_CMF_V13`; identical formal solution, no ALO store —
  audit 4.5). This is the mode with a converged m15 population model
  (`USE_LAM_ES=F` is flagged FAIL in the population history: the lagged
  map did not contract, job 12316397).
- `USE_LAM_ES=F`: the OBS moment solve now calls
  `MOM_J_CMF_NONMONO_V1` — the same Knop sign-split solver the population
  run uses, with `SIGMA_TRUE`, `INNER_BND_METH`, `IB_STAB_FACTOR`,
  `H_CHK_OPTION`, `FREQ_INDX`, `NCF`, and `SET_MOM_NU_NEXT`. Ray moments
  come from `FG_J_CMF_V12` as before; `DO_NONMONO_FREQ_MATCH` auto-sets
  `T`. **[code-consistent; numerical A/B vs populations unverified — the
  same solver failed to contract on m15 in the population context]**

### 2b. Resolved Sobolev profiles

`DO_SOB_FLUX_OBS` triggers the two-phase render (Dessart, Audit &
Hillier 2015, Sec. 3.1): phase 1 solves continuum+BLANK and deposits
each `TRANS_TYPE='SOB'` line's Doppler profile into the stored CMF
opacity/emissivity; phase 2 re-solves the electron-scattering field on
the line-bearing arrays and re-renders `OBSFRAME` with resolved
profiles. Without it `OBSFRAME` is continuum+BLANK only. The
`SOB_MASER_DROP`/`SOB_DROP_EMIS`/`SOB_INJ_EMIS` diagnostics tell you how
much line emissivity survived; a >10% drop means the rendered line flux
is a lower bound, not a prediction.

`DO_SOB_DEP_WINGS`/`SOB_DEP_WING_EXT` (opt-in) add Voigt damping wings
to the deposit — off by default; the Gaussian core under-renders wings
of strong lines. **[unverified]**

---

## 3. Monotonic reference models

Stock settings are the validated configuration — do not set any
`DO_NONMONO_*` key. `GLOBAL_LINE` and `TRANS_*` follow the population
run; `INC_GRID=T` is allowed and unchanged. `DO_NONMONO_FREQ_MATCH`
auto-set does not fire without a nonmono flag.

---

## 4. Common settings

Frequency grids — three distinct controls, do not conflate:

- `MIN_CF`/`MAX_CF` **[physics-window]**: CMF frequency coverage. Must
  span the full Doppler reach of the requested observer range plus
  wings; the new edge-clamp census in `OBS_FRAME_SUB_V9` reports any
  residual gap and warns to widen these.
- `OBS_LAM_MIN`/`OBS_LAM_MAX` **[output range]**: observer-frame output
  window in Angstrom. Only restricts the *output*, not the physics.
- `F_LAM_BEG`/`F_LAM_END` **[line treatment]**: wavelength bracket for
  line/force-multiplier accounting; *not* a spectrum window.

Renderer and scattering:

- `DO_RELO=T`, `DO_CMF_RELO=T` **[physics]**: observer-frame and CMF
  relativistic corrections; match the population run for SNe.
- `DIF` **[physics]**: mirror the population run exactly; under
  `DO_NONMONO_SOB`/`DO_NONMONO_KNOP` + `DIF`, core rays get the
  diffusion inner-boundary injection `I(mu)=J+3*mu*H`.
- `NUM_ES` **[numerical]**: `2` for `COH_ES=F` (ES iterated); `1`
  suffices for `COH_ES=T` (single convolution).
- `COH_ES` **[physics]**: mirror the population run.
- `TAU_MAX=20`, `ES_DTAU=0.1`, `INT_METH=STAU` **[numerical]**: ray
  cutoff/optical-depth spacing — starting values from the template;
  convergence to be demonstrated (T6-J).
- `N_INS_OBS` **[numerical]**: ray insertion floor; 3 is the template
  default — run 3/5/7 for the ray-convergence test before trusting a
  narrow-feature spectrum. **[unverified]**
- `OBS_EXT_RAT=1.3`, `FRAC_DOP_OBS=0.5` **[numerical]**: observed-profile
  half-width (in Vinf) and Doppler-width spacing.
- `dV_OBS_PROF=300`, `dV_OBS_WING=1000`, `dV_OBS_BIG=4000` (km/s)
  **[numerical]**: observer-grid spacings. 300 km/s is the template
  default and a starting point only; narrow Type IIn features need the
  300/150/75 km/s refinement comparison (T6-M). **[unverified]**
- `N_TYPE=N_ON_J` **[solver option]**: used by both moment solvers;
  stable at velocity reversals per the template comment.
- `METHOD=LOGMON`, `FG_OPT=INT/INS` **[solver option]**: template
  defaults; `FG_OPT=INT/INS` is what the KNOP auto-set enforces on the
  population side, keep it consistent in OBS.
- `THK_CONT`, `VTURB_MIN/MAX`, `TURB_LAW`, `DIF`, `THK_LINE`
  **[physics]**: mirror the population model; do not retune for OBS.

## 5. Output semantics and distance scaling

- `OBSFRAME`: final observer-frame spectrum — continuum+BLANK(+SOB
  profiles when `DO_SOB_FLUX_OBS=T`), `Observed intensity (Janskys)` at
  **1 kpc**, frequencies monotonically decreasing (checked at run time).
- `OBSFLUX`: phase-1 product = continuum+BLANK only (plus `Luminosity`
  per depth). Unchanged by phase 2.
- `OBS_FREQ`, `IP_OBS`, `RTAU_DATA`, `ZTAU_DATA`: observer frequency
  vector and per-ray diagnostics from the last render.
- Distance scaling: `F(d) = F_1kpc * (1 kpc / d)^2`. The `6.59934*R(1)^2`
  factor is `1e23 * 2pi * 1e20 / (3.0856e21)^2` — Jy at 1 kpc.

## 6. What is *not* validated

No entry in this guide is backed by a converged numerical spectrum at
this revision. Specifically unverified: `N_INS_OBS`, `dV_OBS_PROF`,
`OBS_EXT_RAT`, `FRAC_DOP_OBS` convergence; `USE_LAM_ES=F` OBS `RJ` vs
population `J`; Sobolev deposit fidelity; luminosity conservation;
`DO_SOB_DEP_WINGS`. These are the T6/T8 runtime tests. Do not cite this
guide as a validation result.
