# KNOP maturity gate — control-flow map (P8.3)

17-Nov-2026. Records where the lagged-field state lives, who can
write the spectrum, and where the gate landed. Implementation:
`obs/knop_spec_gate.f` (KNOP_SPEC_GATE_EVAL), wired before the
`DO SOB_PHASE` loop in `cmf_flux_sub_v5.f`; control key
`DO_KNOP_SPEC_GATE` (default T) in `rd_cmf_flux_controls.f` and
`cmf_flux_cntrl_var_mod.f`. Gate state lives in the module, not a
new counter.

Post-implementation corrections to the map below: `mom_refresh_conv`
stays F on the FG-refresh KNOP route so it is evidence-only, and
`mode=LAMBDA_GS` corroborates the frozen state. The predicate gates
on `n_frozen` / `frozen` / `mode` / `n_mature` plus the stall marker
and manifest termination.

## 1. State writers

| Record | Writer | Where |
|---|---|---|
| `GS_STATE` | `SOLVE_FOR_POPS` every call, KNOP only | `new_main/mod_subs/solve_for_pops.f:3114-3158` |
| maturity predicates | `GET_FG_GS_MATURITY` + `GET_MOM_GS_MATURITY` | `solve_for_pops.f:1088-1104` |
| freeze streak | `N_FROZEN` (consecutive forced-LAMBDA calls) | `solve_for_pops.f` locals, restored at 678-741 |
| refresh outcome | `mom_refresh_conv` flag in GS_STATE | `solve_for_pops.f:3141` |

GS_STATE keys (written this order; `n_frozen` is first because the
wrapper greps `^n_frozen=` -m1):

```
n_frozen  iter  frozen  mode  pops_held
maxch_app  maxch_attempt  se_res  se_res2
fg_bad_frac  fg_churn_frac
mom_bad5_frac  mom_bad50_frac  mom_negj_frac  mom_quar_frac
mom_nhuge  mom_maxchg  mom_refresh_conv
ba_rej_consec  ba_sat_cnt  lum_ratio
n_mature  n_rel  n_trial_acc  n_trials_rej  n_releases
n_raw_releases  n_precond  n_frozen_tot  n_precond_tot
n_rel_noprog  n_churn_streak
```

## 2. What "mature" means in code

- `GS_MATURE_TOL = 0.05` (`solve_for_pops.f:204`).
- Mature checkpoint = `GS_OK` AND last-sweep relative change
  `GS_CHG_FILE < 0.05` AND zero capped writes `GS_CAP_FILE = 0`,
  plus the MOM-side equivalent when `MOM_RAN`
  (`GS_RWARM_SKIP`, `solve_for_pops.f:1110-1117`).
- `N_MATURE` counts consecutive calls passing the raw maturity gate
  (`solve_for_pops.f:491-497`); published as `n_mature=` in GS_STATE.
- `GS_FG_IMMATURE` evaluated per call (`solve_for_pops.f:1276`).

## 3. Freeze/release cycle and stall detection

- `GS_FREEZE=T` forces `LAMBDA_ITERATION` while the lagged J_NEXT
  field matures; releases on maturity, re-freezes on drift
  (`GS_UNSAFE_TOL` path in `solve_for_pops.f`).
- Refresh calls: `FG_KNOP_REFRESH_V1` after each solve under
  `DO_NONMONO_KNOP` (`control_variable_mod.f:1145`);
  `MOM_KNOP_REFRESH_V1` under `DO_KNOP_REFRESH`
  (`rd_control_variables.f:3437`, auto-set on the KNOP/mom routes).
- Job-level watchdog (already exists): after `batch.sh` the generated
  slurm reads `GS_STATE`; `n_frozen >= KNOP_STALL_MAX` (60 default)
  touches `.CMFGEN_KNOP_STALLED` and exits 7
  (`cmfgen_base/submit_cmfgen.sh`, post-run block).
- Exit codes already allocated in slurm: 6 early-exit/stale RVTJ,
  7 knop_stall, 9 manifest verify reject.

## 4. Spectrum write sites (all inside cmf_flux.exe)

| File | Write site |
|---|---|
| `OBSFLUX` | `obs/cmf_flux_sub_v5.f:2455` |
| `OBSFRAME` | `obs/cmf_flux_sub_v5.f:2746` and `:3167` |
| `OBSFRAME` (alt drivers) | `obs/obs_frame.f:429`, `obs/obs_frame_v2.f:522` |

cmf_flux runs in `MODEL_DIR/obs/` via `batobs.sh`; `GS_STATE` is one
level up at `../GS_STATE`. `DO_NONMONO_KNOP` is already readable in
cmf_flux (`cmf_flux_cntrl_var_mod.f:200`, parsed by
`rd_cmf_flux_controls.f:518-528` and enforced consistent with the
model VADAT by `prepare_CMF_PARAM.py`).

## 5. Output-path coverage

| Path | Gate present? |
|---|---|
| Auto OBS after FINAL (`afterok` dependency, submit_cmfgen.sh:2823/2836) | YES - FINAL exits 7 on stall, OBS never starts |
| `skip_final` / OBS-only submission | NO - no model run to gate on; manifest finish only classifies |
| Manual `batobs.sh` / cmf_flux rerun in a run dir | NO - nothing reads `../GS_STATE` |
| In-model spectrum writes during iteration | N/A - OBSFLUX/OBSFRAME are not written by cmfgen_dev.exe |

The job-dependency gate covers the automatic chain only. Any
obs-invocation that does not route through `afterok` can render a
spectrum from stalled or immature populations.

## 6. Proposed gate location (for the deferred commit)

Single authoritative point: **inside cmf_flux, before the OBSFLUX /
OBSFRAME writes**, gated on `DO_NONMONO_KNOP` with a VADAT opt-out
(`DO_KNOP_MATURITY_GATE`, name TBD). It reads `../GS_STATE` (no new
convergence counter - the existing solver state, per the design
constraint):

| Evidence | Classification | Spectrum handling |
|---|---|---|
| `DO_NONMONO_KNOP=F`, or no readable `../GS_STATE` | not applicable | production allowed (SOB routes carry no lagged field) |
| `n_frozen >= KNOP_STALL_MAX`, or `.CMFGEN_KNOP_STALLED` present | knop_stall | block production output |
| `frozen=T`, `n_frozen > 0`, or `mode` = frozen-LAMBDA state | immature lagged field | block production output |
| `mom_refresh_conv=F` | refresh immature | block production output |
| population-runaway evidence (`ba_rej_consec`, `lum_ratio` wild) | population runaway | block production output |
| otherwise | converged | allowed |

"Blocked" should still permit a diagnostic render when an explicit
flag is set, but the output file must carry a machine-readable
maturity stamp in its header so a stalled result cannot pass as a
converged spectrum (ties into P10.1 self-describing output).

Open item for the commit: confirm the `mode=` vocabulary written by
MODE_NOW in GS_STATE (FULL/LAMBDA/...) so the frozen-state test reads
the right token, and decide the exact field set the predicate needs
beyond `n_frozen`, `frozen`, `mom_refresh_conv`.

## 7. Implemented gate - corrected precedence (17-Nov-2026)

KNOP_SPEC_GATE_EVAL (obs/knop_spec_gate.f) implements the gate at the
spectrum-writing boundary in cmf_flux.exe. The first version was
audited and corrected the same day:

- Manifest scan is segment-aware: the authoritative model-execution
  record is the last runs[] segment with "label":"FINAL"; open
  segments carry null exit_code/termination (the OBS render's own
  open segment must not poison the classification).
- GS_STATE now carries run_id= (SLURM_JOB_ID, 'local' for non-slurm
  runs; solve_for_pops.f). A converged-looking state whose run_id
  does not match the FINAL segment's slurm_job_id - or that has no
  run_id/manifest at all - classifies KNOP_PROVENANCE_UNVERIFIED (5),
  closing the warm-restart stale-state hole.
- .CMFGEN_KNOP_STALLED now records the stalled execution's job id;
  a stale marker cannot stall-classify a later run (kept on disk).
- Precedence: GS_STATE missing/malformed -> 4; FINAL open/abnormal
  -> 4; solver state -> 2/3/1; converged downgraded to 5 when
  provenance cannot be verified. States 2-4 block under
  DO_KNOP_SPEC_GATE=T; state 5 renders with marked output.
- 'knop_stall' termination defers to the solver-state stall
  classification; mom_refresh_conv, ba_rej_consec and lum_ratio
  remain evidence-only (written to KNOP_SPEC_STATUS, never gate).
