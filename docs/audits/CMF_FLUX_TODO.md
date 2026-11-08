# CMF_FLUX TODO — non-monotonic spectral correctness

Tracks the audit/fix work item by item. Each task lists files, the
correction, tests, results, and limits. A task is only done when the code
compiles **and** the stated test shows the behaviour, not merely when it
builds. See `CMF_FLUX_SPECTRAL_AUDIT.md` for findings and severity.

Status legend: `[x]` done+verified, `[~]` in progress, `[ ]` open.

---

## Completed

### T1 — Establish execution path and audit OBS_FRAME_SUB_V9
- Files read: `cmf_flux_v5.f`, `cmf_flux_sub_v5.f`, `obs_frame_sub_v9.f`,
  `ins_line_v5.f`, `ins_line_obs_v4.f`, `COMP_JCONT_V4.INC`,
  `rd_cmf_flux_controls.f`, `check_cmf_flux_param_consis.f`,
  `mod_cmf_obs.f`, `mod_freq_obs.f`, `subs/extend_vtsigma.f`, the J-solver
  set (`fg_j_cmf_v12/v13`, `mom_j_cmf_v8`, `*_nonmono_v1`), and
  `comp_j_blank.f` for the population dispatch.
- Result: full chain RVTJ/POP -> controls -> grids -> CMF J -> Sobolev
  deposit -> `OBS_FRAME_SUB_V9` -> `OBSFLUX`/`OBSFRAME` documented in the
  audit doc. Executable is `exe/cmf_flux.exe`.
- Limits: read-only audit; no runtime trace of a live model yet.

### T2 — Fix extended-grid turbulent-velocity scaling (audit 4.1)
- File: `obs/cmf_flux_sub_v5.f` (~line 1257, `IF(ACCURATE)` block).
- Correction: replaced `TA=VTURB_MIN+(VTURB_MAX-VTURB_MIN)*VEXT/VEXT(1)`
  with the `LAW_V1` normalization `MIN(1,|VEXT|/max|VEXT|)` plus the
  `max|VEXT|>1e-12` floor.
- Why: `VEXT(1)` can be 0/negative/non-max on a non-monotonic field;
  corrupts `VDOP_VEC_EXT` and all Doppler widths on the extended grid.
- Tests: `make cmf_flux_sub_v5.o` + full `make` -> clean build, both
  exes relinked. Monotonic equivalence argued analytically
  (`VEXT(1)=max|VEXT|`). **Numerical A/B (ACCURATE=F vs T, monotonic)
  still required — see T6.**
- Limits: does not address `EXTEND_VTSIGMA` log-V interpolation (T3) or
  the `LAW_TAU1`-vs-extended-grid inconsistency (noted, unchanged).

### T3 — Widen non-rel CMF frequency window for sign-changing VMU (audit 4.2)
- File: `obs/obs_frame_sub_v9.f` (~line 665).
- Correction: `MAX_VMU=MAXVAL(VMU_RAY)` -> `MAXVAL(ABS(VMU_RAY))`.
- Why: signed max under-covers the blue CMF window when `VMU` goes
  negative beyond the positive max; edge clamps then substitute boundary
  opacity silently.
- Tests: `make obs_frame_sub_v9.o` + full `make` -> clean build. Monotonic
  equivalence argued analytically (mirrored ray makes `|VMU|` symmetric).
  **Numerical confirmation on a sign-changing field still required.**
- Limits: only affects `DO_REL_CORRECTIONS=F`; rel path already used exact
  extremes.

---

## Open tasks (one at a time)

### T4 — Resolve OBS-vs-population `RJ` solver mismatch (audit 4.5) — HIGH
- Question: under `DO_NONMONO_SOB+USE_LAM_ES=F`, OBS uses `MOM_J_CMF_V8`
  while the population solve uses `MOM_J_CMF_NONMONO_V1`. Does `V8` assume
  monotonic `sigma` and give a different `J`?
- Steps: (a) read `mom_j_cmf_v8.f` sigma handling vs `mom_j_cmf_nonmono_v1.f`;
  (b) run a non-monotonic SOB model through OBS with `USE_LAM_ES=F` and
  compare `RJ_CMF_ST`/flux against the `USE_LAM_ES=T` (`V12`/`V13`, same
  formal solution) result; (c) if divergent, route the OBS moment solve to
  `MOM_J_CMF_NONMONO_V1` under `DO_NONMONO_SOB` (small dispatch change in
  `COMP_JCONT_V4.INC`, mirroring the `DO_NONMONO_KNOP` precedent).
- Test: monotonic regression (must be bit-identical when flags off) +
  non-monotonic A/B.
- Files: `obs/COMP_JCONT_V4.INC`, possibly `plane/mom_j_cmf_v8.f`.

### T5 — Reject `ACCURATE` for non-monotonic `V` (audit 4.3/4.4) — MEDIUM
- `EXTEND_VTSIGMA` uses `LOG(V)`; breaks for `V<=0` and geometrically
  smooths reversals; extended solvers assume monotonic `sigma`.
- Precedent: `DO_NONMONO_KNOP` already forces `INC_GRID=F`.
- Steps: in `CMF_FLUX_SUB_V5`, after the nonmono detector, if
  `ACCURATE .AND. N_NONMONO_OBS>0` (with a nonmono flag set) `STOP` with a
  clear message — OR scope it to `V<=0`/reversal in the interpolation
  stencil. Decide whether `DO_NONMONO_SOB` on a *monotonic* field may keep
  `ACCURATE` (probably yes; the detector only fires on real reversals).
- Test: monotonic `ACCURATE=T` must still run; non-monotonic `ACCURATE=T`
  must fail explicitly, not silently.
- Files: `obs/cmf_flux_sub_v5.f` or `obs/rd_cmf_flux_controls.f`.

### T6 — Validation suite A-J (Phase 7) — REQUIRED before claiming correctness
- A static atmosphere vs known RT solution.
- B monotonic vs stock/reference `cmf_flux.exe` (regression for T2/T3).
- C non-monotonic (>=2 extrema): finite spectrum, no discontinuities at
  reversals, coverage sufficient.
- D multiple resonance locations along one ray: all contribute.
- E continuum-only: integrated flux consistent with input luminosity.
- F isolated line vs independent reference.
- G overlapping lines: neither omitted nor double-counted.
- H coherent vs incoherent electron scattering wings.
- I `ACCURATE=F` vs `T` monotonic; nonmono extended-grid only after T5.
- J grid convergence: depth, `OBS_FREQ`, ray (`N_INS_OBS`) refinement.
- For each: record inputs, controls, pass criterion, result. Update this
  file and the audit doc.

### T7 — Optional `CMF_FLUX_AUDIT` diagnostic output (Phase 8)
- Separate file (not production output). Contents per the audit spec:
  model info, V range + reversal count, frequency coverage, interpolation
  stats, ray-resolution stats, transitions identified/rendered, continuum
  and line eta/chi totals, ES iteration info, inner-boundary flux,
  integrated flux/luminosity, warnings/failed checks.
- Default OFF; new control in `rd_cmf_flux_controls.f`.
- Files: `obs/cmf_flux_sub_v5.f`, `obs/obs_frame_sub_v9.f`,
  `obs/cmf_flux_cntrl_var_mod.f`, `obs/rd_cmf_flux_controls.f`.

### T8 — Sobolev deposit fidelity check (audit 4.6) — MEDIUM
- Verify `ADD_SOB_LINE_PROF` deposition reproduces the Sobolev optical
  depth through the observer-frame integral on a continuum-resolution CMF
  grid; check `SOB_MASER_DROP`/`SOB_DROP_EMIS`/`SOB_INJ_EMIS` accounting and
  the two-phase `DO_SOB_FLUX_OBS` for double-counting/stale arrays.
- Files: `obs/cmf_flux_sub_v5.f` (Sobolev block ~lines 1400-3100),
  `obs/ins_line_v5.f` (grid coverage for SOB lines).

### T9 — Update `cmfgen_obsidian/` notes for any behavior/control/diagnostic change
- Required by project rules for T5/T7 (new control) and the T2/T3 fixes.
- Files: affected notes under `cmfgen_obsidian/`.

---

## Standing constraints

- Preserve stock monotonic behaviour; new controls default to existing
  behaviour; fail closed on invalid configs; no silent clamps; keep
  OpenMP thread safety; check allocation status; small isolated changes.
- Do not re-audit the JFNK solver or unrelated CMFGEN components.
- Every modification needs a regression test; do not mark done on compile
  alone.
- Preserve the unrelated local edit in `new_main/subs/cmf_blkband_v3.f`
  and all untracked diagnostics; never reset the working tree.
