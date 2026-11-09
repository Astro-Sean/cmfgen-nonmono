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

### T4 — Resolve OBS-vs-population `RJ` solver mismatch (audit 4.5) — RESOLVED (code), runtime A/B open
- Finding (confirmed source defect): under `DO_NONMONO_SOB+USE_LAM_ES=F`
  the population solve uses `MOM_J_CMF_NONMONO_V1` (Knop sign-split,
  `SIGMA_TRUE`, lagged `J_NEXT`) while `COMP_JCONT_V4.INC` called
  `MOM_J_CMF_V8`. `V8` divides by `1+sigma*K_ON_J` and couples every cell
  to the bluer frequency only, so on reversing/`sigma<-1` regions the OBS
  continuum `J` feeding the emissivity was not the population field.
- Correction (`obs/COMP_JCONT_V4.INC`, non-extended branch):
  `IF(DO_NONMONO_SOB)` inside the `USE_LAM_ES=F` leg calls
  `SET_MOM_NU_NEXT(NU(MIN(ML+1,NCF)))` then `MOM_J_CMF_NONMONO_V1` with the
  `COMP_J_BLANK` argument order (`SIGMA,SIGMA_TRUE`, `INNER_BND_METH`,
  `IB_STAB_FACTOR`, `H_CHK_OPTION`, `FREQ_INDX`, `NCF`); `COHERENT=L_TRUE`
  and `TB` keep the OBS convention (incoherent ES folded into `ETA`, only
  coherent scatterers in `ESEC`). Stock `MOM_J_CMF_V8` retained in the
  `ELSE` for `DO_NONMONO_SOB=F`. Ray moments still come from the
  `FG_J_CMF_V12` formal solve that precedes the moment call, matching the
  population V13 pairing.
- Companion fix (same file, `FIRST_FREQ` block): one-sided
  `dLOG_NU=LOG(FL/NU(MIN(ML+1,NCF)))` under
  `DO_NONMONO_SOB .AND. .NOT. USE_LAM_ES`; the nonmono solver reads
  `dLOG_NU` in `chi_hat` even under `INIT` while `FL_OLD` is undefined at
  `ML=1`. Mirrors the `COMP_J_BLANK` fix. Other branches never read
  `dLOG_NU` at `FIRST_FREQ`.
- Companion change (`obs/rd_cmf_flux_controls.f`,
  `obs/cmf_flux_cntrl_var_mod.f`): `DO_NONMONO_FREQ_MATCH` read +
  forwarded to `SET_FG_FREQ_MATCH`/`SET_MOM_FREQ_MATCH`, with the same
  auto-set as `RD_CONTROL_VARIABLES`: absent key -> T under
  `DO_NONMONO_KNOP` or `DO_NONMONO_SOB+USE_LAM_ES=F`, explicit F
  overridden with a banner. Without it the OBS grid (different NCF and
  spacing) reads the lagged `J_NEXT` off the characteristic foot.
- Lagged-field state: `J_NEXT_FRQ`/`H_NEXT_FRQ` cold-start guarded
  (`J_NEXT_VALID=F`) on the first OBS pass and refill through the ES
  iterations; `RESTORE_MOM_JNEXT_GS` also picks up a `MOMJNXT_GS`
  checkpoint when present in the run directory.
- Verified: `make all` clean (ifort -O3+AVX512), `exe/cmf_flux.exe`
  relinked 09-Nov 11:30, `nm` shows `mom_j_cmf_nonmono_v1` symbols and
  the new `Calling MOM_J_CMF_NONMONO_V1` string. Source pins:
  `unit_tests/python/test_obs_nonmono_dispatch.py` (37 checks).
- **Remaining: numerical A/B on a converged non-monotonic SOB model —
  OBS `RJ` under the new dispatch vs `USE_LAM_ES=T` (V12/V13 formal J)
  and vs the population `RJ`; monotonic regression must show V8-path
  bit-identical behaviour with the flag off.**

### T5 — Reject `ACCURATE` for non-monotonic `V` (audit 4.3/4.4) — RESOLVED (code), runtime test open
- Implemented in `obs/cmf_flux_sub_v5.f` directly after the non-mono
  detector: `IF(ACCURATE .AND. N_NONMONO_OBS .GT. 0) STOP` naming
  `EXTEND_VTSIGMA` `log(V)` interpolation, the `SIGMAEXT` `sigma+1`
  floor, and the monotonic-sigma extended solvers.
- Check depends on the detected field (`N_NONMONO_OBS>0` = any reversal
  or `sigma<-1`), not on a flag, so it also closes the
  `ALLOW_NONMONO_OBS` escape; `DO_NONMONO_KNOP`'s `INC_GRID=F` key-level
  rejection in `rd_cmf_flux_controls.f` is unchanged.
- Monotonic `ACCURATE=T` unaffected under every flag.
- **Remaining: runtime matrix (monotonic+ACCURATE runs; reversing
  field+SOB+ACCURATE aborts; reversing+SOB+ACCURATE=F proceeds).**

### T5b — OBS_FRAME edge-clamp census (audit 4.2 follow-up) — DONE (code)
- `obs/obs_frame_sub_v9.f`: all three bounded clamp sites (loc 1, loc 2,
  loc 2b) now update `N_EDGE_CLAMP`, `CLAMP_FRAC_MAX` (vs SM window) and
  `CLAMP_GRID_MAX` (vs full `FREQ_CMF`), tracking worst `(ML,LS,I)` under
  `!$OMP CRITICAL`; report printed after the LS loop, with a
  `MIN_CF/MAX_CF` warning only when a clamp landed outside the real grid.
  Silent boundary-opacity substitution is now visible in the log.
- **Remaining: exercise on a non-monotonic model to confirm counts.**

### T5c — OBS output sanity checks — DONE (code)
- `obs/cmf_flux_sub_v5.f`: `OBS_FREQ` checked finite and strictly
  decreasing after `INS_LINE_OBS_V4` (abort on violation — the bracket
  hunts assume monotonic decreasing); new `CHK_OBS_FLUX` aborts on
  NaN/Inf and censuses negative bins before the `OBSFLUX` and both
  `OBSFRAME` writes. Negative bins are reported, not repaired.
- **Remaining: confirm on live output; add wavelength-range coverage
  report to T7 if needed.**

### T5d — COMP_OBS NU_STORE retention window — DONE (code), live OBS rerun open
- Finding (confirmed live defect, latent upstream): `COMP_OBS_V2`/`V3`
  store solved CMF intensities in `NU_STORE`/`IPLUS_STORE`, a sliding
  window of `NST_CMF` entries (20000 OBS / 6000 population).  When full,
  the shuffle kept only ~2 entries above `OBS_FREQ(OBS_INDX)`.  Rays
  with `FREQ_CONV_FAC>1` (mu below ~beta/2) need CMF frequencies up to
  `OBS_FREQ*MAX_FREQ_CONV_FAC` -- for VINF~27,000 km/s that is ~9% above
  `OBS_FREQ` -- so once the store filled the required entries had been
  discarded.  Before the T5c-era guard this interpolated a degenerate
  bracket into a silent NaN flux; after it, OBS jobs 12893194/12893195
  (`SN_A_p2p85r150_SOB_d39` pair) STOPped with `CMF_FREQ outside range`
  at `CMF_FREQ=1.47719 > NU_STORE(1)=1.47151`.
- Correction (`subs/comp_obs_v2.f`, `subs/comp_obs_v3.f`): save
  `MAX_FREQ_CONV_FAC=MAXVAL(FREQ_CONV_FAC)` at init next to
  `MIN_FREQ_CONV_FAC`, and raise the shuffle retention top from
  `OBS_FREQ(OBS_INDX)` to `OBS_FREQ(OBS_INDX)*MAX_FREQ_CONV_FAC`.  The
  `NU_STORE(NEXT_ST_LOC-2)` bottom gate is unchanged; the
  `CMF_FREQ outside range` STOP remains as backstop for genuine
  `MAX_CF` undercoverage or `NST_CMF` overflow (message now names both).
  `OBS_FREQ` is strictly decreasing, so a top retained for the current
  bin also covers all later bins.
- Tests: `unit_tests/python/test_comp_obs_retain.py` + new driver
  `unit_tests/f90/drv_comp_obs_retain.f90` -- NST=64, 250 CMF freqs, a
  mu=-0.5 ray (`conv~1.055`): old criterion reproduces the production
  STOP signature, fixed criterion completes all 40 bins with fluxes
  matching an independent evaluation to 2.8e-16.  `test_comp_obs.py`
  still 6/6.
- Limits: window width roughly doubles (log range
  `MAX_CONV/MIN_CONV` ~0.15-0.18 for VINF~3e4 km/s); `NST_CMF` headroom
  is ample at 20000 but only ~3x at 6000 for the population path --
  a too-small window fails loudly (`NST_CMF too small`), not silently.
  **Remaining: interactive-partition OBS rerun of the failed
  `SN_A_p2p85r150_SOB_d39` job to confirm end-to-end.**

### T10 — Recommended input guide + docs + regression run — DONE
- `docs/audits/CMF_FLUX_OPTIMAL_INPUT_SETUP.md` written: separate
  KNOP / SOB / monotonic configurations, trap inventory
  (`NON_MONO_SOLVER` dead key, `SOB_SIGMA_MAX` now clamps the OBS
  moment solve too, stale cached grids), common settings marked
  mandatory/physics/numerical/unverified, output semantics and the
  Jy-at-1-kpc distance scaling.
- Obsidian updates (09-Nov-2026): `reference/VADAT keys.md` (OBS
  freq-match key, INC_GRID rejection, stale "KNOP OBS removed" claim
  corrected), `physics/Model pipeline.md` (OBS solver consistency
  pass), `solvers/Sobolev non-mono solver.md` (dispatch-fix banner).
- Committed and pushed `8224a00b..4e5a64a6` (10 commits, one per
  file). Untracked pre-existing edits (`test_knop_lagged_field.py`
  and run products) left in place.
- Regression run: `unit_tests/run_all_tests.py` -> **31/31 files
  pass**, including `test_obs_nonmono_dispatch.py`.
- Limits: T6 model-based validation ladder still unrun; all source
  and build evidence only.

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
