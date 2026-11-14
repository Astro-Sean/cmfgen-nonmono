# AGENTS_HISTORY.md - cmfgen_jfnk status log and reference inventory

> STALE-EVIDENCE NOTICE (20-Oct-2026): the solver stack was reworked
> through merge `bd23277c` (no-progress stops, precondition gate,
> KNOP bulk floor, churn arm).  Runtime observations and convergence
> verdicts in entries dated before that merge reflect older builds
> and are retained for provenance only - they are NOT reliable
> statements about current-run behaviour.  Treat them as history,
> not guidance.

This file holds the sections moved out of AGENTS.md on 10-Oct-2026
so AGENTS.md stays under the ~16 kB agent rule-injection limit.
Update conventions are unchanged: dated status entries (newest first
within the Status section) and the diff-vs-pristine inventory live
here. AGENTS.md section 7 points at this file.


## 7. Full diff vs pristine (as of 29-Sep-2026)

Modified: `.gitignore`, `Makefile`, `Makefile_definitions`,
`new_main/cmfgen.f`, `new_main/cmfgen_sub.f`,
`new_main/mod_subs/Makefile`,
`new_main/mod_subs/comp_j_blank.f`,
`new_main/mod_subs/control_variable_mod.f`,
`new_main/mod_subs/do_var_cont.f`,
`new_main/mod_subs/line_mod.f`,
`new_main/mod_subs/mod_cmfgen.f`,
`new_main/mod_subs/rd_control_variables.f`,
`new_main/mod_subs/solve_for_pops.f`,
`new_main/mod_subs/sub_sob_line_v3.f`,
`new_main/subs/Makefile`,
`new_main/subs/auto_add_ion.f`,
`new_main/subs/check_ions_present.f`,
`new_main/subs/do_ng_accel_band_v4.f`,
`new_main/subs/rd_sn_data.f` (NISO=-1 sentinel so
`Number of isotopes: 0` is accepted - same 2-line fix as V6;
upstream conflated "header missing" with "zero isotopes"),
`new_main/subs/rv_sn_model_sniin.f` (SIGMA finite-difference
denominator fix, 29-Sep-2026, see status note below),
`new_main/subs/solveba_v13.f`,
`obs/cmf_flux_cntrl_var_mod.f`, `obs/cmf_flux_sub_v5.f`,
`obs/cmf_flux_v5.f`, `obs/COMP_JCONT_V4.INC`, `obs/mod_cmf_obs.f`,
`obs/obs_frame_sub_v9.f`, `obs/rd_cmf_flux_controls.f`,
`plane/get_moms_non_rel.f`, `plane/Makefile`,
`spec_plt/Makefile`, `subs/Makefile`, `subs/gff.f`,
`subs/scr_read_v2.f`,
`subs/sobew.f`, `subs/sobew_grad_v2.f`, `subs/sobew_grad_v3.f`,
`subs/sobjbar_sim.f`.

OBS / spectrum code (18-Sep-2026, DO_NONMONO_SOB plumbing):
`obs/cmf_flux_cntrl_var_mod.f` (flag, init F),
`obs/mod_cmf_obs.f` (SIGMA_TRUE array),
`obs/rd_cmf_flux_controls.f` (reads [DO_NONMONO_SOB], echoes mode +
warns DO_CMF_EW keeps clamped sigma),
`obs/cmf_flux_v5.f` (allocates SIGMA_TRUE, re-derives it from R,V after
RD_RVTJ_VEC - RVTJ stores clamped sigma only),
`obs/cmf_flux_sub_v5.f` (SOBEW_GRAD_V2 gets SIGMA_TRUE when flag set;
COMP_EW_ABS_EW keeps clamped SIGMA - it divides by GAM and is not
sigma < -1 safe),
`obs_template/CMF_FLUX_PARAM_INIT` (new F [DO_NONMONO_SOB] key).

OBS / Knop formal-solver plug-in (19-Sep-2026, DO_NONMONO_KNOP):
`obs/cmf_flux_cntrl_var_mod.f` (flag, init F),
`obs/rd_cmf_flux_controls.f` (reads [DO_NONMONO_KNOP], echoes mode),
`obs/COMP_JCONT_V4.INC` (standard-grid branch calls
FG_J_CMF_NONMONO_V1 with SIGMA_TRUE+NCF when flag set; extended-grid
site stops with an explicit error -- SIGMA_TRUE_EXT not implemented),
`obs/obs_frame_sub_v9.f` (NI_MAX bound uses sum|dV|+max|V| instead of
V(1): non-monotonic rays traverse a velocity range multiple times;
also per-thread OpenMP work arrays, bounded band hunts, sqrt-radicand
guards -- see 30-Sep-2026 robustness entry),
`obs_template/CMF_FLUX_PARAM_INIT` (new F [DO_NONMONO_KNOP] key).
Pipeline: prepare_CMF_PARAM.py knop_mode_active() stamps
DO_NONMONO_KNOP=T only when non_mono_mode=KNOP.

OBS I_NEXT restart (Oct-2026): cmf_flux builds its own frequency grid,
so the FGINEXT_GS checkpoint cannot be restored by position.
FGINEXT_GS format VER=3 adds the writer's KN_FREQ grid ahead of the
intensity records; CMF_FLUX_SUB_V5 calls SET_FG_GS_OBS_GRID(NU,NCF)
before its frequency loop, and the solver init dispatches to
RESTORE_FG_INEXT_OBS, which log-log remaps the stored field onto the
OBS grid (slot s = field of freq s+1 on either grid).  (KC,LS) indices
carry over: OBS reads NP from RVTJ and builds rays with the same IMPAR
on the same R.  VER=2 files still restore positionally in the FINAL
path but are rejected in OBS (no grid to remap).  The OBS workdir is
MODEL/obs/, so the restore tries FGINEXT_GS then ../FGINEXT_GS.

Added: `new_main/subs/recompute_sigma_true.f`,
`plane/fg_j_cmf_nonmono_v1.f`, `plane/mom_j_cmf_nonmono_v1.f`,
`subs/sob_rh_couple.f`, `pgplot_stubs/`,
`obs_template/` (pipeline OBS-job template, copied from V6 with
`NON_MONO_SOLVER` stripped - key does not exist in dev7),
`build.sh`, `build_debug.sh` (29-Sep-2026; strict-exit wrappers --
nodebug now builds `make all`, see Sect. 5),
`docs/`, `cmfgen_obsidian/`, `knop_solver_files/`,
`knop_solver_files_full/` (review snapshots), this file plus
`README.md`. Project documentation lives in `cmfgen_obsidian/`
(audits, history, solvers, control, reference); keep new `.md`
notes inside the vault.

Also built by `./build.sh`: `exe/cmf_flux.exe`, `exe/obs_frame.exe`
(the `all` target orders `d_obs` after `d_dust`/`d_disp`; the narrower
`cmfgen` target does not build obs).

## 8. Status

### Banded TDMA memory pass: factorization verified, second-admission OOM (24-Oct-2026)

The low-memory `JFNK_PRECOND=JACOBI` + `JFNK_BANDED_TDMA=T` diagnostic
reached a usable admission on isolated `viper08` (job `12366624`,
iteration 743).  The block-Thomas factorization passed its
componentwise backward-error check at `5.08E-14` on `N=128730`, so the
compact algebraic factorization is valid.  GMRES still stagnated
immediately (`6.7069E11 -> 6.7065E11`, `info=3`), so the remaining
solver failure is direction quality, not memory or factorization.

Earlier attempts establish the memory sequence:

- job `12364717` (`JFNK_PRECOND=LU`, `JFNK_USE_BANDED=T`,
  `JFNK_BANDED_TDMA=T`) OOMed during setup at about 31 GB.
- compact `KL+KU+1` LUAB storage lowered job `12365268` to about
  27.3 GB, but viper07 was shared by four other CMFGEN jobs and the
  node-level OOM still fired.
- isolated job `12365656` ran the older LU-base/duplicate-storage
  binary and OOMed at 27.64 GB.
- job `12366624` survived the first low-memory admission, but retained
  activation-local `JFNK_CMAT_RAD_BANDS`, `JFNK_BANDED_LUAB`, and
  `JFNK_BTD_*` arrays through the `JFNK_FRZ_FAIL_LIM=2` /
  `JFNK_FRZ_FAIL_HOLD=8` hold.  The second admission at iteration 760
  allocated the new activation's temporary arrays on top of that stale
  state and was killed at 29.51 GB against the 30G cap.

The source now releases the activation-local banded operator and
TDMA/LUAB storage in `JFNK_FREE_BANDED_STATE_V1`, called from
`JFNK_INIT_V1` and at the end of `JFNK_SOLVE_CORR_V1`.  The factors are
always rebuilt for a new baseline, so retaining them was pure memory
risk.  #status/partial -- a fresh iteration-734 restart is needed to
verify a second admission under the 30G cap.

### FD-JVP + rank-4 secant diagnostic: operator verified, solve still FAIL (24-Oct-2026, job 12364305)

Config `configs/dessart_m15_2p82d_SOB_l14_lte_fds.yaml` (same warm
restart and controls as the secant run, but `JFNK_BA_JVP=F`) admitted
JFNK-FRZ at it735.  The FD-vs-analytic comparison on the first four
Arnoldi directions showed `|FD-BA*v|/|BA*v|` of only `1.4E-3` to
`1.7E-3` for both the diagonal-band and full-band actions.  That
rules out a first-order analytic-BA-JVP defect for this sampled
subspace.  #status/verified

The same four directions still exposed a weak Jacobi inverse:
preconditioner round-trip errors were `1.99E+1`, `2.33E-2`,
`9.60E+1`, and `1.35E+1`.  All four secant pairs were recorded and
the corrected apply reported `ready=T`, but GMRES moved the weighted
residual only `4.9366E10 -> 4.7654E10` before
`gmres_stagnation_abort`.  The returned Newton step had
`|s|_inf=4.0360E-8`, no Newton step committed, and the fallback BA
recovery made the full iteration cost about `766 s` although the
JFNK phase itself was about `65.5 s`.

Verdict: changing the JVP source does not generate the missing
descent direction.  The repeated failures are consistent with the
population-only frozen-J unknown space lacking a component that can
cancel the radiation-field residual, compounded by a poorly
conditioned Jacobi inverse on some Krylov directions.  A
block-Thomas banded-LU comparison was submitted as job 12364717
(`m15_2p82d_SOB_l14_lte_btd`, `JFNK_PRECOND=LU`,
`JFNK_USE_BANDED=T`, `JFNK_BANDED_TDMA=T`, secant off) to isolate
whether a direct inverse of the stored banded operator changes the
verdict; admission is pending at the time of this entry.

### Rank-4 secant JFNK diagnostic: mechanism live, fix FAIL (24-Oct-2026, job 12363588)

Commit `988340ad`, warm-restart config
`configs/dessart_m15_2p82d_SOB_l14_lte_sec.yaml` (npq_gd controls plus
`JFNK_SEC_R=4`; analytic BA-JVP, banded TRI path, Jacobi base PC,
LGMRES/GDEFL/recycling/cycle-stall retained).  Run resumed at it719,
left LAMBDA at it734, and admitted JFNK-FRZ at it735 on demand 1.19%
(vs `JFNK_FRZ_GATE=3%`).  #status/verified

The measured pair store was not contaminated by pre-Newton probes:
the first four eligible JVPs were GMRES Arnoldi directions, and the
preconditioner reported `secant update built r=1..4 ready=T`.  The
linear solve improved modestly relative to the `npq_gd` reference --
first-cycle residual ratio 0.916 vs 0.939 at it30; restart cosine
0.916 vs 0.940 -- but the cycle still stalled at
`hess_res=true_res=4.9125E10`, GDEFL retained four Ritz directions,
and Newton aborted with `info=3 gmres_stagnation_abort`.  The true
weighted residual stayed `5.3623E10`; no step was accepted.  The
fallback BA path then rejected at it738 (`MAXCH=123.46` >
`BA_MAXCH_FAIL=100`) and returned to LAMBDA.

Verdict: the rank-r inverse update is mathematically active and safe
as an opt-in diagnostic, but rank 4 does not repair the stored
banded-BA Krylov wall.  Do not use this run as evidence of model
convergence; it adds another bounded negative result consistent with
the cross-depth radiative-residual diagnosis.

Follow-up diagnostic submitted as job 12364305:
`configs/dessart_m15_2p82d_SOB_l14_lte_fds.yaml` keeps the secant rank
at 4 but switches `JFNK_BA_JVP` to F, so the measured pairs come from
finite-difference JVPs of the actual frozen residual.  The target
hydro file was pre-staged byte-identical before submission to avoid
the orchestrator's new-model inline-conversion miss.

### KNOP OBS formal-J path restored (06-Nov-2026, development 31a582ef)

Reversed the `be742860` (13-Oct) removal of the OBS-side Knop
wiring, adapted to the current tree. `cmf_flux` now dispatches the
non-extended CMF formal solve to `FG_J_CMF_NONMONO_V1` on
`SIGMA_TRUE` when `DO_NONMONO_KNOP=T`, with `SET_FG_NU_NEXT` per
call and `SET_FG_GS_OBS_GRID(NU,NCF)` before the frequency loop so
the version-3 `FGINEXT_GS` checkpoint remaps onto the OBS grid via
`RESTORE_FG_INEXT_OBS`. Controls restored in
`rd_cmf_flux_controls.f` with the same prerequisite checks as the
VADAT reader (BLANK lines, LAM_SET=F, USE_LAM_ES=T mandatory,
INC_GRID=F, DO_CMF_EW=F, no fixed-J/relativistic/plane-parallel
bypasses) plus `PRINT_KNOP_DEBUG` -> `SET_FG_KNOP_DEBUG`. The
fail-closed RVTJ detector now accepts `DO_NONMONO_KNOP` as a
confirmed mode. Ungated guards the path needs were restored:
`VINF=MAXVAL(ABS(V))`, LAW_V1/LAW_TAU1 VTURB bug fixes,
`MU_AT_RMAX` radicand floor, `CHI_FLR` division floor, and the
negative DBB/DTDR inner-boundary guards and diffusion-boundary
OBS_FRAME injection widened to KNOP. `prepare_CMF_PARAM.py`
(project `e3364fcc`) stamps the flag set again.
Verified: clean `./build.sh`, 17/17 tests, solver symbols resolved
in `exe/cmf_flux.exe`; dry-run param generation produces
`DO_NONMONO_KNOP=T` from both the yaml and VADAT-mirror paths.
Runtime smoke test on the converged `m15_2p82d_KNOP_l14_freet`
model (real RVTJ + FGINEXT_GS, 7-min bounded run): controls echoed
KNOP mode, the detector confirmed `DO_NONMONO_KNOP`, OUT_FLUX shows
`Calling FG_J_CMF_NONMONO_V1`, the FGINEXT_GS checkpoint was
accepted (`I_NEXT valid: YES`, tag 1145) with the OBS frequency
remap applied (file NCF=56377 -> OBS NCF=352489), all 70 depths
reported `SIGMA < -1` handled, and the run moved into the
ray-by-ray ES phase (`Using J from formal-solution`) before the
external timeout killed it mid-sweep. Queued npq job 12301378
resolves all of this at runtime (batobs.sh exe resolution +
runtime param regeneration), so no resubmit was needed - it is
now the live full-spectrum test. `#status/partial` -- machinery
runtime-verified on a real model; a complete spectrum file still
requires a full OBS run (npq 12301378).

### m15_2p82d_KNOP_np failure audit + resubmits (05-Nov-2026, later)

Job 12298049 (`_np`, cold KNOP, proxy pops) ran ~6h14m and stopped
at iteration 51 via the GS_PRECOND watchdog
(`50 no-progress pop-preconditioning iterations reach
GS_PRECOND_LIM`, demand pinned ~1.9E15).  The run loaded the
bd23277c-era binary at job start (15:19) -- it predates the
21-Oct stall fixes (`cc3024b9`, committed 20:03 while the job was
in flight) and the 05-Nov tranche.  The stop signature matches the
diagnosed mechanism exactly (per-depth MAJOR scale veto zeroing
corrections, LAMBDA-phase quarantine inert, progress-blind
precond counter); no NEW defect surfaced.  The STEQ residual
concentration at inner depths 66-70 is shared with the converging
SOB np run -- a snapshot property, not a KNOP code bug.  (Note
06-Nov: the "KNOP OBS remains unsupported" statement written here
is superseded - the path was restored in 31a582ef, see the entry
above.)

Resubmits on development `437e631e` (exes resolve at job start):
- KNOP: `dessart_m15_2p82d_KNOP_npq` job 12301377 (queued earlier
  with `DO_BA_QUAR=T`; gets cc3024b9 + the tranche at runtime).
  OBS follow-on 12301378 will now exercise the restored Knop
  formal-J path (31a582ef) instead of stopping at the rejection.
- SOB: `dessart_m15_2p82d_SOB_l14_lte_npq` job 12303015 (fresh
  clone of the np config, same VADAT overrides incl. explicit
  `DO_BA_RES_ACC: F`).  The prior np SOB (12298227) keeps running
  on the old binary.
`#status/partial` -- npq KNOP is the live test of the stall fix;
runtime outcome pending.

### Nonmono adaptation-gap closure: seven-item source-review tranche (05-Nov-2026, development 437e631e)

Commits 012db720..437e631e, one per item (items 6+7 bundled).
No solver physics touched; items 2-4 fill or clamp work arrays
consumed by non-rate diagnostics only.

- SNIIN hydro path (012db720): `RV_SN_MODEL_SNIIN` had two latent
  defects -- `LIN_INTERP` was called with the uninitialized output
  V as the target abscissa and clobbered TA (correct order is
  `LIN_INTERP(R,V,ND,R_HYDRO,V_HYDRO,NHYDRO)`), and the reversal
  loop ran `DO I=1,NHYDRO`, swapping each hydro pair twice
  (identity). Reversal now stops at `NHYDRO/2` so the file-order
  (increasing-R) table matches the decreasing CMFGEN grid that
  LIN_INTERP's strict-x check requires. Not flag-gated: both were
  dead code on every path. `test_rv_sniin.py` rewritten to assert
  correctness on monotone and repeated-V hydro tables (14/14).
- RDINR sigma (bd899027): the review premise (SIGMA unassigned)
  was not literal -- the shared tail assigns SIGMA from the
  hydro-table `MON_INT_FUNS_V2` cubic. The real gap was
  discretization mismatch with `SIGMA_TRUE` (centred differences
  of the final R,V). Under `RDINR .AND. (DO_NONMONO_SOB .OR.
  DO_NONMONO_KNOP)` the tail now calls `SIGMA_TRUE_CORE(SIGMA,...)`
  so the two arrays agree; stock RDINR keeps the cubic derivative.
  New `test_rv_hydro_rdin.py` + `drv_rv_hydro_rdin.f90` (9/9).
- Adiabatic coolant (fbebcac6): `EVAL_ADIABATIC_V3` fills a
  `SIG_AD` work array from `SIGMA`, switching to `SIGMA_TRUE` when
  allocated under either nonmono flag (SIGMA can be the
  forced-monotone-cubic derivative after a TAU_SPACE regrid).
  `B(I)` uses `SIG_AD`; ADIABAT_CHK prints the array actually used.
- Grey-T clamp (ab8f5c5f): `COMP_GREY_V2`/`COMP_GREY_V4` copy
  SIGMA into `SIGMA_GREY` and clamp `>= -0.999` under the nonmono
  flags before feeding `JGREY_WITH_FVT`, `MOM_JREL_GREY_V1/V2`,
  `FGREY_NOREL_V1`, `JGREY_HUB_DDT_V2/V5` -- the grey moment
  structures pole at sigma < -1. Grey is a T-revision diagnostic
  only (never STEQ/rates), so the local clamp is safe. Count logged
  on LUER under `CONV_DIAG_LVL>=1`. `cmfgen_sub.f` needed no
  call-site change: `GREY_T_ITERATE` -> V2 and `COMP_GREY_LST_IT`
  -> V4 funnel through the drivers.
- Dielectronic EW (e39df7c3): the `SOBEW` call gate widened from
  `DO_NONMONO_SOB` to `(DO_NONMONO_SOB .OR. DO_NONMONO_KNOP)`; the
  KNOP transfer solves J on the unclamped gradient, so the
  diagnostic now consumes the same `SIGMA_TRUE`.
- OBS boundary (437e631e, items 6+7): `CMF_FLUX_SUB_V5` now
  detects a non-monotonic RVTJ (non-strictly-decreasing V or
  sigma < -1, with print-precision tolerances) right after
  `RD_CMF_FLUX_CONTROLS`. `DO_NONMONO_SOB=T` confirms machinery and
  continues; `ALLOW_NONMONO_OBS=T` (new presence-probe obs key,
  default F, self-echoed by `RD_STORE_LOG`) downgrades to a
  warning; otherwise a loud LUER error + STOP. Makes silent
  monotonic processing of a KNOP model impossible.
  `mod_cmf_obs.f` SIGMA comment corrected to the actual contract.

Verification: clean `./build.sh` on the final tree, `Build
succeeded`, all three exes present, unit suite 17/17 (includes the
new RDINR test). The VELTYPE=12 and RDINR+KNOP smoke coverage is
synthetic (drv_sniin / drv_rv_hydro_rdin drivers) -- no
`unit_tests/systems/` model fixtures exist for these paths.
Caveat, same as the previous tranche: commits are file-disjoint
subsets of the one build+test-verified final state; intermediate
commits were not individually rebuilt.
#status/verified

### Lagged-field drift floor: diagnostics and gate fixes (22-Oct-2026, development 5a5b6bff)

Seven-item tranche against the lagged-field drift hanging at ~1e-1
vs the 1e-3 convergence target. Commits 352edcf1..5a5b6bff, one per
item. Solver physics untouched; all new output rides LUER under
`CONV_DIAG_LVL>=2` except the same-cell streak line (level 1).

- Per-pass drift split (352edcf1): `MOM_GS_CHG_P1/P2` record the max
  drift on the first (`NEW_FREQ=T`) vs repeat solver call per
  frequency; the `passes` row prints both plus the direction census
  (`MOM_SW_NCROW`/`MOM_SW_NROW`). Refresh re-calls are excluded --
  their residual is `MOM_RF_DCHG`.
- Argmax-jitter history (3da9b008): MOM and FG keep the last four
  sweeps' argmax tuples (`MOM_ARG_HI/HL/HV`, `FG_ARG_HK/HR/HL/HV`)
  and classify `cycle` (alternating values at a fixed cell =
  period-2 composite-map limit cycle), `margnl` (stable monotone
  argmax = contraction tail under-relaxation only slows), `bndry`
  (stable argmax at a boundary depth = HBC/NBC closure). Exported
  via `GET_MOM_GS_CLS`/`GET_FG_GS_CLS`. FG also gained true argmax
  coordinates (`FG_GS_ARGK/ARGR/ARGL`); the sweep store loops now
  record (depth,ray,freq) of the max drift instead of MAX-only.
- Accelerator echo (7f2530e7): the KNOP VADAT echo now reports
  `mom_rf_aa=`, `mom_om_v2=`, `mom_damp_min=` (via new
  `GET_MOM_DAMP_MIN`). Evaluation: `MOM_RF_AA` and
  `MOM_ADAPT_OMEGA_V2` already default T unconditionally, so no
  absent-key auto-set was needed.
- Boundary tag (dd4d4608): the drift row appends `bnd_row` when the
  argmax sits at depth 1/2 or ND-1/ND (H midpoint rows).
- `GS_SETTLE_TOL` (1464101c): new VADAT key, default -1 = off and
  the gate is then bit-identical. With a positive value, a
  `cycle`-classified floor at or below it exempts the store-churn
  term from refreezing (`L_SETTLE`); a cycle above it forces the
  LAMBDA phase via `settle_exceed`. `GS_UNSAFE` untouched.
- Bulk share (79e559f1): `MOM_GS_NLOW` counts cells above the 1e-3
  target; the `bulk` row prints the share. Gate-Bulk redirection of
  the gate's max-norm read was evaluated and NOT committed -- it
  stays a diagnostic until a run shows tail domination.
- Same-cell churn streak (5a5b6bff): `SOLVEBA_V13` exports the
  dominant demand argmax (`BA_ARG_J/I/PCT`); consecutive LAMBDA
  calls applying >10% to the same cell grow `N_SAMECELL`. At 2+ it
  bumps `N_CHURN_STREAK` early and waives the churn arm's
  `GS_NTOT>0` requirement (the bulk census can read 0 when the field
  never refresh-converges). JFNK, rejected, and non-LAMBDA calls
  reset it; exports are zeroed before each solve so stale argmax
  cannot fabricate a streak.

Verification: build clean; unit suite 16/16. `drv_mom_sweep.f90`
now runs the dump under `CONV_DIAG_LVL=2` with seeded pass,
bulk-share and argmax-history values plus three calls producing
`cycle` then `bndry`+`bnd_row`, and calls `DUMP_FG_INEXT_GS` with a
`cycle` seed -- runtime coverage of every new format.
#status/verified

### USE_LAM_ES un-forced; SOB defaults to the moment path (21-Oct-2026)

User directive: allow `USE_LAM_ES=F`, make it the SOB default, no
silent parameter changes, ALO and accelerators on by default.

`rd_control_variables.f` enforcement rewritten:

- **SOB**: `USE_LAM_ES` no longer forced to T. Absent or F runs the
  moment path (`MOM_J_CMF_NONMONO_V1` + fixed-J/`GET_MOM_ALO_DIAG`/
  `GET_MOM_CELL_QUAR` -- the KNOP+LAM_ES=F recipe); an explicit T
  keeps the V13 formal J. Both modes print a `====` banner to
  OUTGEN naming the active field solver; the `GLOBAL_LINE`
  override prints a banner before the `SOB auto-set:` line.
- **KNOP**: absent-key default stays T (verified pairing); an
  explicit `USE_LAM_ES=F` is now honoured with a banner (was
  silently flipped to T).
- **New SOB absent-key auto-sets** (explicit values always win):
  `DO_NONMONO_FREQ_MATCH=T` under `LAM_ES=F` (same lagged-map gain
  argument as the KNOP auto-set), `DO_FG13_ALO=T` under `LAM_ES=T`
  (V13 diagonal OAB ALO), `DO_BA_RES_ACC=T`, `DO_BA_RES_REUSE=T`
  (only when `DO_BA_RES_ACC` lands T), `NG_ABSTAIN=T`.

Pipeline: `_NON_MONO_LINE_SOB` (`prepare_VADAT.py`) and the
`run_model_grid.py` fallback now stamp `USE_LAM_ES=F` for SOB; the
selector forces T only on `cmfgen.use_lam_es: true` or an explicit
`vadat_overrides` key, and prints the pairing in all cases.

Audit notes carried forward: the 13-Oct frequency-grid analysis
(edge doublets at dLOG_NU~2E-11 dominating the advection
coefficients, second-order clamped-SIGMA K_ON_J mismatch) is now
quoted in the solver-mode banner and AGENTS.md; TRIDIAG cannot
help under non-mono (do_var_cont zeroes off-diagonal TX/TVX, so
TRI=DIAG) -- a tridiagonal linearization of MOM_J_CMF_NONMONO_V1
does not exist.

Tests: `test_vadat_keys.py` updated (explicit F honoured, absent
defaults F) -- 55/55; `test_knop_status_text.py` prose check
updated for the new banner text. Build clean on
viper_intel_nodebug.

Runtime status of SOB+LAM_ES=F is UNVERIFIED on this build; the
GS freeze/watchdog stack is KNOP-scoped, so a bare run's only
backstops are `MAX_CHNG` and the resacc rollback.

### No-progress stops and churn-arm gating (20-Oct-2026, development bd23277c)

Branch knop-no-progress-stops, merged --no-ff.  No-progress stops
and churn-arm gating for the m15
stall class.  N_FROZEN_TOT/N_PRECOND_TOT accumulate frozen and
preconditioned calls across release/refreeze refunds and reset only
on measured demand progress (progressive release, demand back under
GS_PRECOND_MAX, or LST_ITERATION); GS_FREEZE_LIM and GS_PRECOND_LIM
(KNOP auto-set 50) test the totals, and new key GS_REL_NOPROG_LIM
(KNOP auto-set 10) stops on counted no-progress releases -- the
earliest bound the 5-6-iteration release cycle reaches.  All three
totals ride GS_STATE (n_frozen_tot=, n_precond_tot=, n_rel_noprog=)
and are adopted on restart.  KN_RAY_BULK_TOL KNOP default raised to
6e-2 (above the measured 5.2e-2 refresh floor) and scoped to KNOP
(non-KNOP absent-key default 0, explicit values respected); BTOL<=0
now means 'no floor' in the commit gate as it already did in the
sweep exits.  The GS_PRECOND gate forces LAMBDA_ITERATION/FIXED_T
(the log claimed live LAMBDA but mode selection could take BA) and
prints a [PRECOND] sentence; GS_RAW_NUDGE no longer fires on
preconditioned calls.  The FG_TRIAL_CHURN arm is gated on its own
N_CHURN_STREAK >= FG_TRIAL_CHURN_MIN (default 3; 0 disables) -- the
freeze streak never reached FG_TRIAL_MIN because releases refund
it.  MOM sweep-table format fixed (15/10 descriptors for 14/9
args); status-text ES13.4 buffers; unsafe sentence says 'holding
populations' only when the hold is on.  Verified on
cmfgen_validation/knop_noprog_m15 (warm restart of the test10 it860
stall state, low limits): four pop_precond calls all ran LAMBDA
with the [PRECOND] sentence, N_PRECOND_TOT hit GS_PRECOND_LIM=4 and
the run stopped non-convergent at it865 instead of burning the
requeue; the refresh printed 'bulk rays: maxres=5.5e-2 floor=6e-2'
with 24 bulk exits (was 0 at the 2.5e-2 floor).  GS_STATE carries
the new counters.  Build clean, 14/15 unit-test files pass (the 2
pre-existing DO_NONMONO_SOB auto-set cases remain).

### FG_TRIAL_CHURN: churn-only arm on the FG marginal trial (20-Oct-2026)

New VADAT key `FG_TRIAL_CHURN` (optional, default 0.0 = off, stock
path unchanged; has effect only under `DO_NONMONO_KNOP` with the
`USE_LAM_ES=T` field path).  Commit `c917e695`.

Why: the m15_2p82d KNOP run (test10) reached a persistent
LAMBDA/GS hold - not the old release/refreeze limit cycle.  The
refresh floored at residual 0.06-0.15 (the documented lambda ~ -1
mode), the bad-cell share sat at 4-8e-4 (below the 1e-3 bulk
bound), but store churn asymptoted at 1.3-2.1e6 of ~6.5e8
significant cells, i.e. 2-3x the strict `MAX(50,NTOT/1000)` bound.
`GS_FG_IMMATURE` therefore never cleared, `FG_TRIAL` could never
arm (it requires a fully mature gate), and the run held
populations until the iteration budget exhausted at it860.

With `FG_TRIAL_CHURN > 0`, `solve_for_pops.f` snapshots
`GS_FG_PRECHURN` before the store-churn OR term and computes
`FG_CHURN_MARG`: the gate is immature ONLY through churn, and the
churn census sits within `FG_TRIAL_CHURN` times the bulk bound.
The FG-only marginal-release trial may then arm with every other
condition unchanged (no UNSAFE tier, resacc snapshot armed,
`N_FROZEN >= FG_TRIAL_MIN`, no luminosity block, trial budget
remaining).  The probe is the same damped LAMBDA step
(`T1=GS_LAM_CAP`), judged by the resacc accept/refund machinery -
a release that re-poisons the field rolls back like any other
rejected trial.  OUTGEN prints 'Ray field churn-marginal' with
the churn census, and the iteration reason is `fg_trial_churn`
(vs `fg_trial_lam`) so the audit scripts can distinguish the arm.

Suggested production value for the m15-class asymptote is ~3;
chosen through controlled testing, not assumed.

Tests: `test_knop_lagged_field.py` static checks (declaration,
0.0 default read, pre-churn snapshot, factor comparison, trial
admission, distinct reason); `test_vadat_keys.py` parse-and-land
case.  The same commit repairs the `test_vadat_keys` helper list
(`knop_status_mod.f90` was missing after the status-text merge
added `USE KNOP_STATUS_MOD` to `mom_j_cmf_nonmono_v1.f`, so the
driver build failed rc=1).  Known pre-existing failures left in
place: 2 `DO_NONMONO_SOB` auto-set cases in `test_vadat_keys` and
12 `test_knop_status_text` prose checks (verified on clean
development before this change).

#status/verified (clean ./build.sh, zero errors; exes rebuilt;
new-key parse test green; runtime arm on m15 pending)

### KN_NO_GAIN: chi floor at chi_scat under KNOP (18-Oct-2026)

New VADAT key `KN_NO_GAIN` (optional, default F, stock path
unchanged; has effect only under `DO_NONMONO_KNOP`).  When T,
`cmfgen_sub.f` floors `CHI` at `CHI_SCAT` immediately after the
stock `NEG_OPAC` fixups, at both sites (main iteration loop and the
post-iteration final pass), and marks the cell
`NEG_OPACITY`/`NEG_OPAC_FAC=0` exactly as the stock branch so the
BA linearization stays consistent.

Why: the `SRCE_CHK`/`ESEC_CHK` floors can leave total `CHI <
CHI_SCAT`, which makes the formal-solver source coefficient sum
`(CHI_SCAT+|GAM|/dLOG_NU)/(CHI+|GAM|/dLOG_NU) > 1` -- the cell is a
net photon amplifier (effective albedo `THETA` up to ~10).  The
bounded `FG_J_CMF_V13` path keeps nothing across great iterations,
but the KNOP lagged `I_P/I_M_DEPTH_NEXT` store recycles the
amplified intensity, so the gain compounds (~8.6x/iter measured).
Evidence: w19 run job 12268434 diverged with populations AND
temperature held; 100% of the runaway luminosity sat in one far-IR
line region (186.3-186.6 kA, OBSFLUX 6.4E16 Jy vs 1.5E5 Jy in the
SOB run).  Flooring `CHI` at `CHI_SCAT` caps every source
coefficient sum at 1, so the lagged map is non-expanding in the
sup norm.  Accepted cost: real maser (stimulated-emission)
amplification is suppressed at those frequencies; the inversions
come from unconverged trace-ion Rydberg levels, so this is a
diagnostic-grade loss, not a physics regression.

A per-iteration census prints to LUER after the luminosity line
(`KN_NO_GAIN: floored N cells (depths d1-d2); max chi_scat/chi= R
at nu= X (1e15 Hz) lambda= Y A`), accumulated in
`KN_NOGAIN_CNT/DMIN/DMAX/RMAX/NU` and reset each great iteration.
The key is echoed in the `KNOP effective` block in
`rd_control_variables.f` (`no_gain=`).  Armed in the production
w19 config `configs/dessart_m15_2p82d_KNOP_w19.yaml`.
#status/unverified (built; runtime arm pending)

### Non-mono/JFNK self-consistency enforcement (16-Oct-2026)

`RD_CONTROL_VARIABLES` now auto-configures dependent switches when
`DO_NONMONO_KNOP`, `DO_NONMONO_SOB` or `DO_JFNK` is set, replacing the
former hard STOPs for reconcilable mismatches.  Each override is
logged to LUER as `KNOP auto-set:` / `SOB auto-set:` /
`JFNK auto-set:` / `Trace auto-set:`.  The enforcement block sits
after the last dependent-key read (before the SOB/KNOP validation
section); `INC_GRID` is read later, so its KNOP override lives at the
INC_GRID read.

- KNOP: `GLOBAL_LINE=BLANK`, `LAM_SET=F`, `FG_OPT=INT/INS`,
  `USE_LAM_ES=T`, `DO_KNOP_REFRESH=T`, `N_LAM_RWARM=4` (when 0),
  `GS_REL_STABLE=2` (when 1), `INC_GRID=F` (auto-set, was STOP).
- SOB: `GLOBAL_LINE=SOB`, `USE_LAM_ES=T`.
- JFNK: `JFNK_MAX_SEC=3600` (when 0), `JFNK_LOGX/PTC/LS/CAUCHY=T`;
  a zero `JFNK_MXCH_TOL` warns instead of being raised.
- Trace (either non-mono flag, only when left at off defaults):
  `BA_DEM_PFLOOR=1e-30`, `BA_DEM_RELFRAC=0.01`, `BA_TR_MAXSWP=5`,
  `BA_TR_ABSMAX=1e-8`.
- Still STOP: `FG_OPT=DIFF` under KNOP, `JFNK_COUPLED=T` under SOB,
  both non-mono flags T, the V-coordinate reads (DDT/PURE_HUB/
  DO_HYDRO/GAMRAY/DO_CL) and the relativistic/plane-parallel modes.

The `USE_LAM_ES=F` moment-solve opt-in and `DO_KNOP_REFRESH=F` A/B
comparison are overridden by this pass (user decision 16-Oct-2026):
the moment path diverged on m15 under both pathways (jobs
12160690/12161004 warm, 12161660/12161668 cold) and is retained in
source only as reference code.

Verified: `./build.sh` clean (zero errors); a symlinked VADAT smoke
test of `exe/cmfgen_dev.exe` with deliberately inconsistent keys
logged all 15 expected auto-set lines in OUTGEN and ran on to the
iteration setup.  Not yet exercised through `run_model_grid.py`.

Follow-up (same day): SOB now emits a prominent warning at read that
SOB+USE_LAM_ES=T takes continuum J from the expansion-biased stock
formal solver with no within-sweep refresh - DO_NONMONO_KNOP is the
verified path on strongly non-monotonic fields.  JFNK admission in
`solve_for_pops.f` gained frozen-tier hysteresis relief: under
DO_NONMONO_SOB, once the LAMBDA-mode controller proves `MAXCH_CTRL <
VAL_EXIT_LAM` for `N_STABLE_LAM` consecutive calls, the frozen tier
admits with `JFNK_NREQ_EFF=1` (demand level and contraction tests
unchanged).  Previously every LAMBDA-designated call reset
`JFNK_FRZ_HITS`, so `JFNK_GATE_NREQ>1` was unreachable and JFNK
deferred on every eligible-era call (m15 SOB_l14_jfnk OUTGEN).

### JFNK review branch: restart cap, clipped reach, opt-in keys (12-Oct-2026)

Commits aa69daa..b4a6fa2 on master implement the JFNK_REVIEW.md
findings from probes 12155376/12155427 on
`m15_2p82d_SOB_jfnk_probe` (same deterministic it444 admission as the
failed production run 12153985).  Summary:

- GMRES now restarts on Arnoldi lucky breakdown/stagnation with a cap
  (JFNK_MAX_RESTART, probe used 2).  RESTART DIAG prints hess vs true
  residual: cycle 1's apparent 37% drop was illusory (hess 5.67e13 vs
  true 3.43e14, ratio 6.06) because the clipped FD JVP is not a fixed
  operator.  Second admission fell back to the inexact direction
  (pred_res/ref 0.87) instead of a blind Lambda abort.
- Every JFNK bail-out prints `jfnk_exit_reason`; nits line now carries
  maxch_est.
- New opt-in VADAT keys: JFNK_RSCALE=DIAG (row scaling by BA diag,
  Ne row by ED; EUCLIDEAN still default), JFNK_MXCH_TOL (max relative
  correction convergence test, default 0=off), JFNK_POP_FL (freeze
  near-empty levels via ROW_REG, default 0=off), JFNK_RIDGE_PIV
  (escalate ridge on measured pivot ratio, default 0=off).
- REACH/REACH-FD now price the commit-clipped step (Z_CLIP) and print
  max|z| before/after plus nclip; the |z|<1 restricted print is
  relabeled (|z|>1 zeroed).  cos(BAz,A_FDz) JJZ_DP corruption fixed;
  verify_cosine_bounds.py guards it (still flags the pre-fix OUTGEN).
- Per-depth LU pivot ratio + DGECON rcond printed per factorization;
  ridged=70 traced to YAML JFNK_RIDGE=1e-4 (config, not a trigger
  bug).
- The 16% ||Mz-v||/||v|| preconditioner residual is genuine LU
  backward error on blocks with equilibrated cond ~1e15 (pop rows
  carry 4.7e12; T/Ne rows only 4.7e-2 of |v|).  Frozen-row reduction
  is correct but changed nothing; the lever is row-space conditioning
  (POP_FL / ROW_TINY), not ridge.
- Issue 6: topF colprobe (HI J=23 K=70) FD-vs-BA |diff|/|BA|=14.6%
  with sign-flipped downward couplings; consistent with FD
  differentiating the dissolution factor w(n) that stored BA froze.
  topF h-scan added.
- Issue 10 report-only: V monotone-decreasing outward, sigma in
  [-120,-1.05], 18 jumps >20 km/s; S IV f=11 oscillator warning real;
  HeI 23/24 0.002 cm^-1 near-degeneracy benign.
- prepare_VADAT.py (outer repo 964c2f92) now prints the parsed-back
  final VADAT so the log cannot contradict the file.

Build `make DEBUG=nodebug all` clean (0 errors; three exes fresh).
Surrogate suite: only the four pre-existing T3/T4 failures.
unit_tests: 10/11; test_interp_guards flaky on uninitialized output
column jitter (routines untouched since baseline).  Opt-in probe with
the new keys not yet submitted.

### JFNK frozen-tier endgame review (11-Oct-2026)

OUTGEN for `m15_2p82d_SOB_jfnk_zlaw14_fresh` (job 12153985) confirmed
repeated frozen-J GMRES stagnation and a soft-column gap (`|diff|/|BA|`
about 0.83 at it384; stiff top-F column about 0.064). The initial
interpretation that fallback permanently prevented FULL-BA was too
strong: FULL-BA resumed at it399 after the cooldown and reduced demanded
MAXCH from 32.5% to 1.35% by it402, then plateaued around 1.4-1.5%.

The prior fallback forced Lambda and a fresh BA build for finite
stagnation/cooldown; observed cooldown calls took about 366 s versus
about 70 s for most FULL calls. A caller-mode pass-through build reached
the it408 cooldown path but passed an uninitialized `T1` to SOLVEBA and
stopped (job exit 6). The cap is now initialized from `COR_LIM_IN` at
entry. Job 12153986 was canceled by its failed dependency. The corrected
binary has not been exercised in a model run.

Source review found that `STEQ_MULTI_V10` includes direct ED and ion
density BA columns. The opt-in `JFNK_LAG_FACTORS` prototype therefore
freezes `POPION` and `MOD_LEV_DIS_BLK` factors only; it does not freeze
ED/T or direct ion-density columns. It is an experimental residual
change, not verified as the physical Jacobian. Surrogate T18 passes with
live/BA gap 2.11e4 and lagged FD-vs-BA error 5.0e-7; this does not
establish CMFGEN production behavior.

The prepared YAML uses `JFNK_FRZ_GATE=3%`, two contracting solves, and
`JFNK_LAG_FACTORS=T`. It also sets `FIX_T_AUTO=T` and `DO_T_AUTO=T` per
user choice. `DO_T_AUTO` releases global FIX_T at `MAXCH_CTRL<50%`, so
this is an early temperature release at current population demand;
`JFNK_TAUTO` only damps T corrections during JFNK, not ordinary FULL-BA.
Source snapshot committed as `b063fb7`; residual-count telemetry committed
as `8380f24`; production-file tests committed as `8e16f96`; VADAT-key
tests committed as `075e407`; the prepared YAML committed in the parent
repository as `a7696c68`. No redeployment has occurred.
Re-run via `scripts/run_model_grid.py` only after explicit confirmation.

Iteration 407 exposed a telemetry defect: `jfnk_total_f_evals` counts
Newton-core calls but omits residual callbacks inside `JFNK_JVP_V1`, so
OUTGEN said F evaluations=0 with 12 JVPs. The reports now print both
core F calls and `JFNK_N_RESIDUAL_EVALS` (all residual callbacks); T20
guards the log fields.

Verification: production `./build.sh` succeeded after the cap and
telemetry fixes; standalone solver tests 21/21, interface tests 9/9,
harness T15-T20 passed, and the separate production-file suite passed
11/11. Full surrogate suite still reports the existing T3-mild demand
threshold failure and three T4 LU-accuracy failures; these do not touch
the new factor path. Config resolution via `--print-resolved-config`
passed. Runtime behavior of the new fallback and factor mode remains
UNVERIFIED.

- Builds clean; `DO_NONMONO_SOB` live in binary: VERIFIED.
- Optimized and debug builds both complete from a clean tree
  (29-Sep-2026; `exe/` 13:25, `exe_dbg/` 13:22). See the
  29-Sep-2026 status entry for the SIGMA denominator fix and the
  `d_lte` link failure that the OBJ_FILES restriction works around.
- Runtime behaviour (flag off = unchanged; flag on = stable Sobolev on
  a non-monotonic model): UNVERIFIED until a controlled model run.
- VEL_LAW=12 (USE_HYDRO_GRID) runtime validation: BLOCKED -- the
  interactive partition rejected submission (AssocMaxJobsLimit, job
  12023231 cancelled). No `skip_stella` + `vel_law: 12` config exists
  yet; see the 29-Sep-2026 entry.
- `DO_NONMONO_KNOP` in cmf_flux (19-Sep-2026): VERIFIED at runtime on
  the TypeIIn_10d model (OBS job 11842074). FG_J_CMF_NONMONO_V1 ran the
  full formal J (~4 calls/freq, JNU_CLAMP_DIAG 0 events); OBSFLUX is a
  clean 93969-bin spectrum. LIMITATION (fixed same day, see below):
  obs_fin_10 (OBSFRAME ray product) carried ~13% negative/NaN bins.

### obs_frame_sub_v9 bidirectional frequency hunt (19-Sep-2026)

`obs_fin_10` corruption root cause: the per-ray frequency interpolation
hunts the comoving-frequency table with a warm-started index
(`K=INDX(I-1)`) but only walks K DOWN (`DO WHILE(T1.GT.SM_FREQ_CMF(K-1))
K=K-1`).  That assumes T1 = nu_obs*(1-VMU_RAY)*GAM is monotonic along
the ray -- true for monotonic V (VMU swept by Z), false for a
non-monotonic V law where VMU_RAY reverses and T1 jumps back up in
frequency.  The stale bracket then extrapolates ETA/CHI far outside the
table -> negative/NaN intensities (~13% of bins; OBSFLUX unaffected
because it does not do the per-ray transform).

Fix: added the upward hunt (`DO WHILE(T1.LT.SM_FREQ_CMF(K)) K=K+1`)
immediately after the existing down-hunt, so the bracket always
satisfies SM_FREQ_CMF(K-1) >= T1 >= SM_FREQ_CMF(K).  For monotonic V
the added loop never executes (unchanged stock behaviour).  Table-edge
cases still hit the existing STOP guards; the band window already
carries a 1.2*MAX_VMU margin + 10-index pad.  Verified the mirrored-ray
extension (R,V same at same R, Z/VMU negated) is exact for any
spherical V(R), and that the point-1 search is safe (ray point 1 is
always the outermost near-side point, VMU(1)>0).  Same routine serves
both cmf_flux.exe and obs_frame.exe (Makefile OBJ1/OBJ2).

Rebuilt cmf_flux.exe + obs_frame.exe; OBS resubmitted through
run_model_grid.py on TypeIIn_KNOP/config_resolved.yaml (skip_final=T
-> OBS only) as job 11842754.

Follow-up (same day): job 11842754 STOPped at 'Bad K indx (loc 2b)' --
T1 landed 3.3 ppm below the lowest tabulated CMF frequency (edge
discretization undershoot at the receding ray endpoint, ML=93969
I=214).  The up-hunt correctly detected it but aborted.  User added
0.1% edge-tolerance clamps at all three boundary hits (loc 1, loc 2
top -> K=2, loc 2b bottom -> K=SM_NCF): within-tolerance T1 clamps to
the edge bracket (tiny extrapolation of a smooth continuum); genuine
off-grid T1 still aborts with diagnostics.  cmf_flux.exe/obs_frame.exe
rebuilt clean (forced recompile after stale-.o suspicion -- on
NFS/lustre always `rm the .o` before make to rule out clock skew).

INFRASTRUCTURE BUG FOUND + FIXED (run_model_grid.py): the resubmit
wipe's restart_files preserve list had no POINT1/POINT2 -- the
SCR_READ_V2 checkpoint files (R,V,SIGMA,POPS + nonmono GS radiation
state).  They are the ONLY warm-restart state; losing them forces a
cold start.  Worse, POP<SPECIES> files (required by cmf_flux -- no
POP* => zero species => N_LINE_FREQ=0 => segfault in INDEXX) are only
rewritten at CMFGEN *completion*.  After repeated OBS-only
resubmissions wiped the TypeIIn_KNOP run dir, OBS jobs failed with
the INDEXX segfault (11843842).  Fix: POINT1/POINT2 + POP* (already
globbed) + VADAT (skip_final only) now all preserved; verified live
on an M15 OBS-only submission.

Recovery: no POP*/POINT* backups exist for TypeIIn_KNOP, so the only
path is a FINAL rerun to completion.  Per user request: resubmitted
configs/TypeIIn_10d_DEV7_KNOP.yaml with init_num_its=100 (FINAL
completes at it100 unconverged -> writes POP* -> chained OBS tests
spectral generation on the velocity-reversal state).  FINAL=11844120,
OBS=11844121; IN_ITS NUM_ITS set to 100 directly (the slurm only
re-instates FINAL_NUM_ITS on requeue/zero, the static IN_ITS=60 would
have capped at 60).

Outcome (19-Sep-2026, second half):

- FINAL ran cleanly through 3 requeue cycles (it1 cold start, ~4-5
  min/it, zero capped/NaN GS writes throughout, cells>0.5 collapsed
  706M -> ~10M).  User stopped the population solver manually at it62
  (IN_ITS NUM_ITS=0); CMFGEN finalised, wrote all 8 POP* + RVTJ +
  POINT1/2 and exited COMPLETED.  Chained OBS 11844121 ran but
  STOPped at 'Bad K indx (loc 2b)' with a NEW signature: T1=3.4954e-3
  vs SM_FREQ_CMF(SM_NCF)=3.4990e-3 at I=308, i.e. 0.101% below the
  edge -- just over the 0.1% tolerance.

- Root cause: not roundoff this time.  FREQ_CMF(NCF)=3.49897e-3 is the
  absolute grid bottom (END_INDX(93)=NCF, pad exhausted): the CMF
  frequency grid built for the SE solve under-covers the
  observer-frame Doppler reach at deep far-side ray points on the
  non-monotonic field.  dev6's equivalent up-hunt
  (K.LT.SM_NCF guard, no error check) silently extrapolated below the
  table -- one source of its corrupted bins.

- Fix (obs_frame_sub_v9.f): replaced edge extrapolation with a
  boundary clamp at all three edge sites.  T1 within 2% of the edge
  is clamped to the edge frequency so T2=0 (or 1 at the top) selects
  the boundary opacity exactly -- bounded, cannot go negative, and
  the deep-ray-point contribution is negligible anyway.  Grossly
  off-grid T1 (>2%) still aborts with diagnostics.  Comments updated
  to state the true mechanism (grid undercoverage, not roundoff).

- OBS resubmitted OBS-only via SKIP_FINAL=1 (job 11850776): all
  preserved files verified restored (POP*x8, RVTJ, VADAT, OUTGEN,
  POINT1/2).  COMPLETED in 15:47.

- RESULT: obs_fin_10 = 93,906 bins, ZERO negative / NaN / Inf / zero
  bins (previous run: 5,986 negative + 6,640 NaN ~= 13% corrupted).
  Zero spike bins >5x the 21-bin median.  8 dips <0.2x cluster at
  4923.5 A (He I 4922) and 5018 A (Fe II 5018) -- physical narrow
  absorption, not artifacts.  Balmer series strong: Ha peak 2.74e7 Jy
  (9.4x local continuum, dominant line), Hb 9.4x, Hd 4.8x.
  Qualitatively Type IIn-like: strong Balmer emission on a
  line-rich continuum, dominant Ha.

- Caveats: spectrum is from the it62 unconverged state (diagnostic
  only); MEANOPAC 'zero or negative flux-mean opacity' warning is
  non-fatal at this convergence level; the edge clamp used boundary
  opacity for the small set of off-grid deep ray points; line
  attribution not yet performed.

### obs_fin vs OBSFLUX discrepancy: broad blue emission is REAL; CV_BOUND under-converged (19-Sep-2026)

User report: obs_fin_10 spectrum shows an additional broad blue
emission component on many lines (e.g. Halpha wing at v ~ -2000 to
-5000 km/s), suspected gridding bug.  Investigation with the per-ray
diagnostics (WR_IP -> IP_DATA, WR_CMF_IP -> IP_CMF_DATA, IP_FG_DATA
file-existence trigger -> solver CV_BOUND per ray per freq):

- On monotonic M15, obs_fin/OBSFLUX agree at ratio ~1 everywhere.
  On the non-monotonic TypeIIn they differed by 55-94x in the Ha
  blue wing -- a non-monotonic-path effect only.

- Per-ray comparison (cont 5200A/7000A): IP_DATA (obs_frame
  independent ray integration) vs CV_BOUND*2 (=IPLUS -> IP_CMF_DATA
  -> OBSFLUX).  Agreement ~1.0 for rays with p > ~3e5 (units 1e10
  cm); interior rays p < R_shell ~ 1.4e5 differ by 300-500x.
  The discrepancy boundary is exactly the dense interaction shell
  (rho jumps to ~5e-13, V jumps 160->5800 km/s at R~1.4e15 cm,
  tau_es~1 to the shell top).

- DECISIVE test -- bolometric boundary flux.  Model luminosity
  L=9.85e8 Lsun => int_0^1 I mu dmu = 2H = L/(8*pi*R1^2) ~ 5.4e8
  cgs.  IP_DATA integrates to ~1.3e8 (order-consistent; residual
  factor ~2-4 from unconverged it62 state + quadrature).  CV_BOUND
  integrates to ~6.4e6 -- ~43-84x BELOW the luminosity the
  atmosphere must emit.  Energy conservation is violated by the
  solver's boundary intensity; obs_frame does not over-produce.

- Cross-check: int CV_BOUND dmu = RJ(1) = 7.54e-9 at 5200A exactly
  -> the solver is internally consistent (CV_BOUND is the same
  intensity field that fed J for SE).  RJ depth profile shows
  radiation trapped below the compressing shell: RJ(75)~8.5e-6
  vs RJ(1)~7.5e-9 (~1100x contrast).  obs_frame's interior
  intensity ~1e-5 ~ shell-depth J -- physically correct emergent
  intensity for a tau~1 line of sight to the hot shell.

- Mechanism: in compressing zones (GAM<0, 70/150 depths incl. the
  shell DeltaV~5600 km/s jump) the effective opacity is
  chi_hat = chi + |GAM|/dLOG_NU -- enormous across the jump, so
  DTAU_hat walls off the interior.  The compensating source
  |Q|*I_NEXT (lagged Gauss-Seidel upwind intensity, INEXT_RELAX=0.5)
  re-injects photons advected in from the adjacent frequency.  At
  it62 the lagged fields were still far from converged -> the shell
  acts as a wall -> CV_BOUND shows only outer-wind emission
  (~2e-8 ~ surface S) and OBSFLUX shows spurious deep P-Cygni-like
  absorption where obs_fin correctly shows shell emission.

- VERDICT: the broad blue emission component in obs_fin_10 is REAL
  -- it is line emission from the approaching shocked shell + fast
  inner ejecta (V ~ 4500-8000 km/s toward the observer), not a
  gridding artifact.  obs_fin_10 (OBS_FRAME_SUB_V9) is the
  trustworthy product.  OBSFLUX (CV_BOUND from FG_J_CMF_NONMONO_V1)
  is the defective product until the GS coupling converges.

- Pipeline bug found+fixed: FGINEXT_GS / MOMJNXT_GS (the GS lagged
  I_NEXT/J_NEXT checkpoints) were NOT in the run_model_grid.py
  preserve list -- every resubmission cold-started the coupling,
  preventing convergence through the workflow.  Added both to
  restart_files (POINT1/2 do NOT contain this state despite the
  comment claiming so).

- Follow-up: FINAL resubmitted (11851691, NUM_ITS=200, cold start as
  single_final always is) to let the compressing coupling converge;
  OBS 11851692 chained.  When it runs, re-create IP_FG_DATA (touch)
  in obs/ before the formal solve to re-check CV_BOUND bolometric
  recovery vs IP_DATA.

### SN-model parameter audit (18-Sep-2026)

Parameters checked for SN-snapshot correctness in the generated
VADAT/IN_ITS (M15_HeOCa sanity run):

- `IT_ON_T`, `DO_GT_AUTO`: grey T overwrite of hydro T - OFF (see
  note below).
- `FIX_NE`: was T (pinned hydro Saha estimate). Now F - hydro Ne is a
  rough estimate; CMFGEN solves Ne from charge conservation. Changed
  in the single_final + final profile dicts and the sanity config.
- `INIT_TEMP_V2`/`T_INIT_TAU`: dead path for `SN_T_OPT=USE_HYDRO`
  (T/ED set by SET_ABUND_CLUMP, untouched). Would re-interpolate T on
  a tau scale if SN_T_OPT=USE_T_IN - avoid that option.
- `FIX_T=F` + `FIX_T_AUTO=T` + `DO_T_AUTO=T` (IN_ITS): staged
  release - T held fixed until MAXCH<50%, then freed. Reasonable for
  a hydro-sourced T; OUTGEN "Temperature held fixed" is this
  mechanism, not a bug.
- `LTE_EST=T`: cold-start LTE guess only; populations are iterated.
  Fine.
- `ST_INT=40`/`END_INT=55`, `ND_QUAD=0`: extended-grid refinement
  range - dead because `N_INS=0`. If N_INS is ever >0 these restrict
  refinement to depths 40-55; for ejecta the whole grid should be
  refined instead.
- `GREY_TAU=10`: only consumed by the IT_ON_T grey loop - inert.
- `D2_EQ_D1=T` (IN_ITS): pins outermost BA depths - standard outer
  BC helper. Fine.
- `T_EXC=0.7` = T_EXCITE_MIN: floor for LTE departure coefficients.
  Fine.
- `ADJUST_R_DEFAULTS`: live tau-regridding control file (required
  when REV_RGRID=T, ignored when F). Now excluded from run-dir copies
  with ADJUST_CORRECTIONS/IT_SPECIFIER/table*.dat - a REV_RGRID=T run
  without it fails loudly instead of silently regridding.
- `SN_DEBUG=T`, `SN_DEBUG_LUM=T`: verbose SN diagnostics only -
  cosmetic.
- `DIF=F` + `IB_METH=ZERO_FLUX` + `LUM_METH=ZERO` + `LSTAR=0`:
  consistent zero-flux inner boundary for ejecta with no photosphere.
- `DO_NOMONO_KNOP=T`/`DO_NOMONO_SOB=F` stamped by the pipeline: V6
  keys, ignored by dev7 readers (dev7 flag is DO_NONMONO_SOB).

- OBS spectrum transfer for SOB populations (06-Oct-2026):
  `prepare_CMF_PARAM.py` stamps both DO_NONMONO_KNOP and DO_NONMONO_SOB,
  with GLOBAL_LINE=SOB and all TRANS_*=SOB. The OBS reader dispatches
  Knop formal J for the continuum/electron-scattering source and retains
  SIGMA_TRUE Sobolev EW/frequency handling. Main VADAT population modes
  remain mutually exclusive. `COH_ES=F` is used in the matched SOB and
  frozen-J configs; prepare_CMF_PARAM propagates it and ensures NUM_ES>=2
  plus IIn wing extents. Build verified; runtime OBS validation pending.

### Sanity run note (18-Sep-2026, jobs 11810955/11810956)

The config `configs/M15_HeOCa_25d_DEV7_sanity.yaml` (Stockholm
M15-He-O-Ca-evo snapshot, t=25.10 d, ND=150, monotonic homologous
ejecta) gives a **good starting setup for Stockholm-style models**:
CMFGEN passes VADAT checks, reads `SN_HYDRO_DATA` (incl. NISO=0 after
the rd_sn_data fix), starts Lambda iterations and does not immediately
blow up. Key settings that get it this far: `STOP_IF_BP=F` (pipeline's
REL_OBS=F/USE_LAM_ES=T trips the upstream SN check otherwise),
`LSTAR=0`, `IB_METH=ZERO_FLUX`, `VTURB=50`, `IT_ON_T=F`,
`DO_GT_AUTO=F`, `REV_RGRID=F`, `N_OB_INS=0`, `DIF=F`, 8 CPUs on the
interactive partition. Convergence to a finished spectrum is still
UNVERIFIED.

Note IT_ON_T = F works much better

**Why IT_ON_T must stay off for LSTAR=0 ejecta models (investigated
18-Sep-2026):** `IT_ON_T=T` runs a grey-atmosphere T re-estimate in
`SET_NEW_MODEL_ESTIMATES` before iteration 1. `COMP_GREY_V4` builds
TGREY from `HFLUX propto LUM`, i.e. the VADAT `LSTAR`. With `LSTAR=0`
(pure ejecta, no photospheric luminosity) the grey solution has
`HFLUX=0` and TGREY collapses to the `0.95*T_MIN` floor (~2850 K).
The correction `T = w*TGREY + (1-w)*T_hydro` applies with weight
`w = 1-exp(-TA/GREY_TAU) ~= 1` at every optically thick depth, so the
correct hydro T is replaced by the collapsed floor. Electron density,
LTE departure coefficients and initial populations are then computed
on the wrong T, and the first BA corrections explode. IT_ON_T is only
valid for photospheric models with a constant luminosity diffusing
through the atmosphere; for ejecta snapshots leave it F.

### TRIDIAG audit (18-Sep-2026, job 11812352)

Question: is SOL_METH=TRIDIAG (NUM_BNDS=3) stable with the current
static-snapshot setup? Answer: under this configuration the BA matrix
has EXACTLY ZERO off-diagonal depth blocks, so the block-Thomas solve
in CMF_BLKBAND_V3 reduces to independent per-depth solves - TRIDIAG
is numerically identical to DIAG.

Why all off-diagonals vanish here:

- VJ = 0: do_var_cont.f line ~331 `ELSE IF (USE_LAM_ES)THEN` is an
  EMPTY branch (upstream design, same in V6). USE_LAM_ES=T means TX/TVX
  (dJ/dchi,deta matrices) are never computed, so dJ_LOC=0 and VJ=0.
  Every J-coupling term in BA_UPDATE_V7 is VJ-multiplied.
- dZ (CMF resonance-zone line coupling, update_ba_for_line.f): off-
  diagonal entries are built only from dJ_LOC -> zero.
- sub_sob_line_v3 BA off-diagonals: VJ-multiplied -> zero.
- steq_advec_v4: returns immediately (INC_ADV=F).
- Diffusion T-coupling: dead (DIF=F; would need NUM_BNDS=ND anyway).
- All steq_* writes go to band M=DIAG_INDX only.
- CONT_VEL=T -> the F2DA/FC ray-by-ray path (which could fill VJ) is
  not taken.

Verified empirically: TRIDIAG job 11812352 and DIAG job 11811962 on
the same M15_HeOCa_25d snapshot produce bitwise-identical correction
trajectories (it5: 21892.4377234943 both; it6: 135.113165529951 vs
...52, a 1-ulp storage-path difference).

Implications:

- The earlier "TRIDIAG is unstable" observation was confounded by the
  IT_ON_T/DO_GT_AUTO/FIX_NE bugs fixed since - not by the solver.
- TRIDIAG adds no stability benefit in this config either; it is a
  no-op. Its only value here is exercising the banded-solver path.
- For TRIDIAG to carry REAL coupling you would need USE_LAM_ES=F
  (moment-equation J via MOM_J_CMF with clamped sigma) or INC_ADV=T -
  a different, riskier J solution. Not advisable for the non-mono SOB
  path where the moment solver is intentionally bypassed.
- The stale comment in prepare_VADAT.py claiming "TRIDIAG forces
  USE_SELECTIVE_JVAR -> VAR_MOM_J_CMF_NONMONO_V1" describes V6-only
  machinery that does not exist in dev7. Dev7 TRIDIAG uses the
  standard upstream CMF_BLKBAND_V3.
- prepare_VADAT SOB forced-dict still stamps SOL_METH=DIAG
  unconditionally - harmless under this config (equivalent), kept as
  conservative policy until real coupling exists (Knop plug-in).

Standalone test submitted: configs/M15_HeOCa_25d_DEV7_tridiag.yaml
(depth_independent_solver: false), model dir
dev7_sanity/M15_HeOCa_25d_DEV7_TRIDIAG. Generated files verified:
VADAT SOL_METH=TRI, MODEL_SPEC NUM_BNDS=3, N_TYPE=N_ON_J,
IT_ON_T=F, FIX_NE=F, DO_GT_AUTO=F. Jobs 11812352 (FINAL, running) /
11812354 (OBS, pending). Early iterations match the DIAG baseline.

### SN_HYDRO_DATA optional headers (18-Sep-2026)

Stockholm/ejecta-only SN_HYDRO_DATA files carry no "Energy Lum [lsun]"
or inner-boundary headers. The pipeline no longer requires them:

- proxy_populations_v2.parse_sn_hydro_data: derives r_inner_cm and
  t_inner_k from the innermost grid point, and estimates the ejecta
  luminosity as the blackbody at the tau=2/3 surface built from the
  kappa column (_ejecta_lum_lsun; falls back to the inner-boundary
  blackbody when kappa is unusable). Matches the header's established
  photospheric-flux meaning; the STELLA converter wrote |flux| at
  tau=2/3. For M15_HeOCa_25d this gives 3.2e9 Lsun (vs CMFGEN
  converged ~6.8e8 - right order for a proxy T_rad anchor).
- convert_sn_hydro_to_superlite.py: the bogus 1e42 Lsun default
  replaced by the same photospheric estimate.
- smooth_sigma_hydro_data.py: fixed-position header/section parsing
  replaced with keyword scans (no blank-separator assumption);
  energy_lum=None when absent (writer already skips it).
- run_model_grid.py SEDONA prep already falls back to the
  inner-boundary blackbody when the header is missing - unchanged.

### USE_LAM_ES=F / dev6 A/B test (18-Sep-2026)

Audit conclusion on dev6: its STANDARD solver path is bit-for-bit
upstream (do_var_cont.f, cmf_blkband_v3.f identical). USE_LAM_ES=F works
in both trees via the moment-equation J path; dev6's instability traces
to the parallel nonmono stack (SOLVE_FOR_POPS_NONMONO_V1 et al., gated
by NON_MONO_SOLVER = DO_NOMONO_KNOP .OR. DO_NOMONO_SOB), which this tree
deliberately does not import.

USE_LAM_ES=F requirements (upstream check_param_consistency.f): for an
SN_MODEL, one of USE_J_REL / USE_DJDT_RTE / USE_LAM_ES must be T. The
Stockholm reference (tmp/dev6_obs_test/M15_model, converged ~40 its)
ran USE_J_REL=T + REL_OBS=T + REL_CMF=T + TRIDIAG + NUM_BNDS=3 +
IT_ON_T=F - real tridiagonal coupling via VAR_JREL_V5 filling TX/TVX.
IT_ON_T=F was the decisive stability fix: with LSTAR=0, IT_ON_T=T runs
a grey-T re-estimate that collapses the temperature before iteration 1.

Pipeline changes enabling the test:

- prepare_VADAT.py enforce_sn_iin_physics() now takes user_overrides;
  keys the user sets explicitly in vadat_overrides are preserved even
  when not in user_overridable (previously DO_NOMONO_KNOP=T and
  REL_OBS/REL_CMF=F were stamped unconditionally).
- run_model_grid.py enforce_cmfgen_sn_iin_settings(): REL_OBS/REL_CMF
  stamps use _force_unless_user; an explicit YAML vadat_overrides value
  wins.
- plot_sn_hydro_cmfgen_inputs.py: keyword-scan headers; missing Energy
  Lum falls back to the inner-boundary blackbody for the title only.

A/B runs submitted (identical inputs, only the binary differs):

- dev6: configs/M15_HeOCa_25d_DEV6_jrel.yaml, install_dir
  cmfgen_nonmono_dev_V6, model dir dev6_sanity/M15_HeOCa_25d_DEV6_JREL.
  VADAT verified: USE_J_REL=T, USE_LAM_ES=F, REL_OBS=T, REL_CMF=T,
  SOL_METH=TRI, NUM_BNDS=3, DO_NOMONO_*=F (standard stack),
  IT_ON_T=F, FIX_NE=F. Jobs 11813665 (FINAL) / 11813666 (OBS).
  OUTGEN confirms "Calling MOM_JREL_V9".
- dev7: configs/M15_HeOCa_25d_DEV7_jrel.yaml, same solver config, model
  dir dev7_sanity/M15_HeOCa_25d_DEV7_JREL. Jobs 11813746 / 11813750.

Caveat: USE_LAM_ES=F puts J on the moment solver, whose sigma clamp is
inconsistent with non-monotonic fields - for the DO_NONMONO_SOB path
keep USE_LAM_ES=T + DIAG. The JREL+TRIDIAG config is the correct
reference/sanity configuration only for monotonic models.

### Frozen-Jacobian stall diagnosis + fix (18-Sep-2026)

The DIAG sanity run plateaued at MAXCH ~= 0.357 for ~175 iterations and
the job exited at it427 without converging. Mechanism:

- `COMP_BA=F` + `STORE_BA_INV=T` + `N_FIX_BA=0` + `FIX_BA=15`: once
  MAXCH < 15, `COMPUTE_BA` is permanently disabled and every iteration
  solves with the stored inverse of the last Jacobian (a stale-Jacobian
  fixed point). OUTGEN shows `read BAMATPNT` + `BA matrix **NOT**
  computed`.
- `DO_NG=T` could not rescue it: Ng acceleration requires
  MAXCH < BEG_NG (=0.1) and MAXCH sat at 0.357, so NEXT_NG stayed at
  1500 forever.
- Dominant residual was at depth 1 (~0.65% decrease) and depth 42
  (~0.16% increase); the stale inverse never refreshes the coupling
  that drives them.

Fix (YAML only, no source change): `FIX_BA=1.0D0` + `N_FIX_BA=10` in
vadat_overrides. N_FIX_BA=10 selects the periodic-recompute branch
(fresh BA when MAXCH > VAL_FIX_BA on the 10-iteration check, plus the
MAXCH_SUM escape) and FIX_BA=1 delays any freeze until below 1%.

Resubmitted (fresh start - the orchestrator refresh does not preserve
population `_IN` files, so these runs re-converge from scratch, a
cleaner test of the fix):

- DIAG: configs/M15_HeOCa_25d_DEV7_sanity.yaml, jobs 11816623 (FINAL) /
  11816624 (OBS). VADAT verified: N_FIX_BA=10, FIX_BA=1.0D0.
- TRIDIAG: configs/M15_HeOCa_25d_DEV7_tridiag.yaml, jobs 11816688 /
  11816689. VADAT SOL_METH=TRI, MODEL_SPEC NUM_BNDS=3 verified.
- Old BA stores were stashed under dev7_sanity/ba_stash/, not deleted.

### JREL+TRI segfault: upstream bug, not a dev-tree regression (18-Sep-2026)

The dev6/dev7 JREL A/B runs (11813665, 11813746) crashed identically:
~27-33 clean Lambda iterations, first real TRIDIAG BA solve returns
exactly 0.0000 corrections, segfault next iteration
(`forrtl: severe (174): SIGSEGV` + `array index out of bounds`).

Bounds-checked debug build (`make DEBUG=debug` -> exe_dbg/, ifort -O0
-check bounds -traceback; job 11815551 on tmp/jrel_dbg):

- Root cause (VERIFIED): `Error in LINPAC (or BLAS) routine DGETRF_2
  in CMF_BLKBAND, depth=1, IFAIL=11` -> `Error in CMF_BLKBAND_V3 -
  shutting code down`. The LU factorization of the banded matrix hits
  a zero pivot at depth 1 ON THE FIRST TRI SOLVE (debug run started
  in BA mode at it1; production runs reached it only after ~33
  Lambda its). The block-Thomas banded matrix is structurally
  singular at the outer boundary when JREL off-diagonal coupling
  populates NUM_BNDS=3.
- Production fallout (consistent): the failed TRI solve returns
  zero/garbage corrections ("0.0000" before each crash), corrupting
  later state -> `gff.f:101` bounds abort with `IGB=-2147483648`
  because `LOG10(T(I))` evaluated on `T(ND) <= 0` (NaN -> INT
  overflow) via genopaeta_v10.f:146 -> comp_opac -> cmfgen_sub.f:1485
  (DTDR inner-boundary opacity block).
- All TRI-path files (cmf_blkband_v3, solveba_v13, ba_update_v7,
  do_var_cont, var_jrel_v5, generate_full_matrix_v3) are
  bit-for-bit upstream in dev7 - this is an upstream bug in
  USE_J_REL=T + TRIDIAG, consistent with the JREL path being
  unexercised on current upstream. Fix is open; JREL runs should use
  SOL_METH=DIAG until the depth-1 band structure is understood.
- Latent bug found en route (kept, minimal): `check_ions_present.f`
  used uninitialized `DPTH_INDX` to index `POP_SPECIES` when all
  T(K) <= T_VAL on the startup path. Now initialized before the loop.
- `Error opening BAMATPNT in READBA IOSTAT=29` is benign (missing
  saved-inverse file on first start; handled).
- Temporary diagnostic prints in new_main/cmfgen_sub.f (iteration
  start + DTDR entry T(ND)/minT) were REMOVED 18-Sep-2026 during the
  Knop integration (they printed unconditionally on the default path).
  The subs/gff.f bad-T report remains in the working tree - it is in
  the -check bounds build only (exe_dbg); the optimized exe/ build is
  unaffected.

### DO_NONMONO_KNOP integration (18-Sep-2026)

Knop et al. (2009) non-monotonic formal solver ported from dev6 as a
minimal, opt-in path. Only the formal solver was ported - the V6
population solver fork (SOLVE_FOR_POPS_NONMONO_V1, JFNK, runaway
guards, moment solver MOM_J_CMF_NONMONO_V1) was deliberately NOT
imported.

Activation: `T [DO_NONMONO_KNOP]` in VADAT (orchestrator: set
`DO_NONMONO_KNOP: T` in cmfgen.vadat_overrides plus
`non_mono_mode: KNOP` - non_mono_mode alone does not activate the
enforcement block). Main VADAT population mode is mutually exclusive
with DO_NONMONO_SOB (Fortran STOP if both T). Requires GLOBAL_LINE=BLANK
(or NONE with every present ion TRANS_*=BLANK), LAM_SET=F, INC_GRID=F,
FG_OPT=INT/INS, no FIXED_J/rel/plane-parallel modes - all enforced in
rd_control_variables.f; GLOBAL_LINE and LAM_SET gates were added
03-Oct-2026. OBS uses the same prerequisites for Knop-only spectra, but
allows both OBS flags T for Knop formal J with SOB line frequencies.
DO_BA_RES_ACC=F warns.

J-update modes (both valid since 19-Sep-2026):
- USE_LAM_ES=T : J from the formal solution (Lambda iteration).
- USE_LAM_ES=F : J from MOM_J_CMF_NONMONO_V1 - the non-monotonic
  variable-Eddington-factor moment solve (ported from dev6):
  chi_hat = chi+|GAM|/dlnu with GAM=(V/R)c^-1(1+SIGMA_TRUE*K_ON_J),
  PSI/W diagonal uses |1+sigma*f| (always positive), signed coupling
  in the RHS via J_upwind = RSQJNU_PREV (Q>=0) or J_NEXT_FRQ(:,k+1)
  (Q<0, lagged Gauss-Seidel over great iterations).  J_NEXT_FRQ is
  checkpointed to MOMJNXT_GS by DUMP_MOM_JNEXT_GS (SCR_RITE_V2) and
  validated on restore (tag/dims/expanded-grid fingerprint) by
  RESTORE_MOM_JNEXT_GS inside the solver's first call; a failed
  restore leaves J_NEXT=0 and the N_LAM_RWARM warm-up covers it.
  The SOLVE_FOR_POPS warm-up skip now requires BOTH the FG I_NEXT
  and (when the mom solver ran) the MOM J_NEXT checkpoints to be
  mature (GET_MOM_GS_MATURITY accessor, same build-order workaround
  as GET_FG_GS_MATURITY).

  NOT ported (deliberate): VAR_MOM_J_CMF_NONMONO_V1 /
  EDD_J_VAR_NONMONO_V1 / J_NEXT_MD / STORE_J_NEXT_MD - the V6
  linearization solver for the selective-JVAR ALO stack.  Under
  DO_NONMONO_KNOP/DO_NONMONO_SOB the stock VAR_MOM_J_CMF_V12 is
  ALSO bypassed (19-Sep-2026): do_var_cont.f gained a fixed-J +
  diagonal-ALO branch that zeroes TX/TVX and fills only
  TX(K,K,1)=ALO_DJDCHI(K), TX(K,K,2)=ALO_DJDETA(K) from the formal
  solver's own per-frequency dJ/dchi, dJ/deta (GET_FG_ALO_DIAG).
  V12 linearizes the monotonic (clamped-sigma) moment equation and
  produced an all-NaN TX at the first post-warm-up BA iteration
  (M15 it9) - same failure V6's ~5000-line replacement stack was
  built to avoid.  NUM_BNDS>1 under KNOP gets no cross-depth
  J-coupling; port the V6 selective-JVAR stack properly if TRI is
  ever needed rather than re-enabling V12.

Files touched (diff vs cmfgen_reference/cur_cmf):
- plane/fg_j_cmf_nonmono_v1.f    NEW (copy of dev6 file + added
  GET_FG_GS_MATURITY accessor; new_main compiles before plane so a
  module USE is not possible there)
- plane/mom_j_cmf_nonmono_v1.f   NEW (dev6 file minus J_NEXT_MD /
  STORE_J_NEXT_MD which only serve VAR_MOM_J_CMF_NONMONO_V1; added
  GET_MOM_GS_MATURITY accessor)
- plane/Makefile                 fg_j_cmf_nonmono_v1.o +
  mom_j_cmf_nonmono_v1.o + rules
- new_main/mod_subs/control_variable_mod.f  DO_NONMONO_KNOP (def F),
  N_LAM_RWARM (def 4)
- new_main/mod_subs/mod_cmfgen.f NU_NEXT_FREQ scalar
- new_main/mod_subs/rd_control_variables.f  reads + guards;
  USE_LAM_ES=F under KNOP prints an informational note (non-mono
  moment solver) instead of STOPping since 19-Sep-2026
- new_main/mod_subs/comp_j_blank.f  Knop-vs-V13 dispatch,
  GET_MOMS_NON_REL under USE_LAM_ES, FIRST_FREQ one-sided dLOG_NU,
  MOM_J_CMF_NONMONO_V1-vs-V11 dispatch under USE_LAM_ES=F -
  all gated on DO_NONMONO_KNOP
- new_main/cmfgen_sub.f  SET_FG_GS_TAG(NITSF) + SET_MOM_GS_TAG(NITSF)
  after SCR_READ, NU_NEXT_FREQ in both freq loops,
  RECOMPUTE_SIGMA_TRUE gates extended to KNOP
- new_main/mod_subs/solve_for_pops.f  ITS_SINCE_START counter +
  N_LAM_RWARM LAMBDA warm-up (skipped on mature FGINEXT_GS +
  MOMJNXT_GS restores), FORCE_LAMBDA guards the mode controller,
  SIGMA_TRUE refresh gates
- subs/scr_read_v2.f  DUMP_FG_INEXT_GS + DUMP_MOM_JNEXT_GS at end of
  SCR_RITE_V2 (self-guard on unallocated arrays; zero cost stock)

Untested/known limitations:
- N_COMP=0 path only so far (monotonic sanity model): the
  compression-zone I_NEXT machinery is compiled but not exercised.
- USE_LAM_ES=F is now supported under KNOP via MOM_J_CMF_NONMONO_V1
  (ported 19-Sep-2026).  Before the port this config hard-STOPped -
  verified at runtime (18-Sep-2026, job 11818525 exit 6, message
  `DO_NONMONO_KNOP=T requires USE_LAM_ES=T`).  The pipeline does not
  force-correct USE_LAM_ES in KNOP mode (it stays forced in SOB
  mode), so the VADAT carries the user's choice explicitly.
  Validation status: M15_HeOCa_25d_DEV7_knop_MOM (monotonic, Q>=0
  branch only) runs and the MOM GS sweep max change decays;
  TypeIIn_10d_DEV7_KNOP_MOM exercises the Q<0 / J_NEXT branch at
  ~85% of depths - convergence there is the open question.
- FG_J_CMF_NONMONO_V1 warns `SIGMA_TRUE_IN is all zeros` at most once
  per process (SIGMA_ZERO_WARNED flag) when the field is identically
  zero - expected on a homologous model; retained as a diagnostic for
  genuinely non-monotonic inputs where a zero profile indicates a
  missing RECOMPUTE_SIGMA_TRUE.  The warning caught a real bug
  (19-Sep-2026): RECOMPUTE_SIGMA_TRUE's internal gate was
  `DO_NONMONO_SOB`-only, so KNOP-only runs kept the zero-initialised
  SIGMA_TRUE (masked on the homologous sanity model where sigma=0 is
  the true answer).  Gate fixed to
  `.NOT. DO_NONMONO_SOB .AND. .NOT. DO_NONMONO_KNOP`.

Validation (18-Sep-2026): `configs/M15_HeOCa_25d_DEV7_knop.yaml`
(monotonic M15 sanity snapshot + DO_NONMONO_KNOP=T) submitted via
run_model_grid.py, job 11817508. Confirmed in OUTGEN:
`Calling FG_J_CMF_NONMONO_V1`, chi_hat path, USE_LAM_ES branch,
`no compressing points - GS storage skipped`, GS restart fallback to
LAMBDA rebuild, `SOLVE_FOR_POPS: KNOP GS warm-up` forcing LAMBDA for
the first 4 its. Convergence of a full non-monotonic model is still
to be demonstrated (e.g. the TypeIIn 10d snapshot that diverged on
the SOB path).

### TypeIIn 10d non-monotonic SOB test (18-Sep-2026)

Config `configs/TypeIIn_10d_DEV7_SOB.yaml`, model dir
`TypeIIn/cmfgen/10.00days`, snapshot ~9.84 d post-explosion from the
existing TypeIIn STELLA .swd. ND=169 (adaptive grid), sigma_min=-98.06,
38.7% of zones sigma<-1 - genuinely non-monotonic, the intended
DO_NONMONO_SOB test case. Jobs 11817076 (FINAL) / 11817077 (OBS) on
partition `small` (24h): the interactive partition has
MaxMemPerNode=32000 and the estimator asked 35G (ND=169), so the job
sat PENDING(MaxMemPerLimit) until moved. Debug submissions on
interactive must keep the memory request <= ~31G.

Pipeline lesson (non_mono_mode activation): `cmfgen.non_mono_mode: SOB`
alone does NOT activate the orchestrator's non-monotonic block in
run_model_grid.py - activation requires a `DO_NOMONO_*` or
`NON_MONO_SOLVER` key set T in `cmfgen.vadat_overrides`
(run_model_grid.py ~8126). With only non_mono_mode set, the block is
skipped, `GLOBAL_LINE` is stamped BLANK (SN_IIN default) into the
resolved overrides, `n_ob_ins` is not forced to 0, and the forced SOB
key table never lands. The config therefore carries
`DO_NOMONO_SOB: T` + `DO_NONMONO_SOB: T` explicitly. Resolved overrides
verified: GLOBAL_LINE=SOB, all TRANS_*=SOB, DO_SOB_LINES=T,
SOB_FREQ_IN_OBS=T, NON_MONO_SOLVER=T, DO_NOMONO_KNOP=F, USE_LAM_ES=T,
SOL_METH=DIAG, VEL_LAW=11, N_OB_INS=0, DO_LAM_IT=T, DO_LAM_AUTO=T,
FIX_BA=1.0, N_FIX_BA=10. In-job INIT writes VADAT; verify
`T [DO_NONMONO_SOB]` once the FINAL starts.

Other settings: `FIX_T=T` + `FIX_T_AUTO=F` (hold shocked hydro T),
`DIF=T` + `IB_METH=DIFFUSION` (inner-boundary tau ~28 at conversion;
prepare_VADAT auto-selects ZERO_FLUX if tau<50 in the generated file -
check VADAT when written), `VTURB=20`, `LSTAR=0`, `n_ob_ins: 0`
(VEL_LAW=11/RDINR grid match), `proxy_populations: true`.

First submission (jobs 11817076/11817077) aborted at it9: LAMBDA
corrections grew ~1e33 -> ~1e109 %/iteration from it1
(depths 1-4 dominant, `FG_J_CMF_V13 - SIGMA .LT. -1.0` warning at
depth 1, sigma=-1.005) -> "bad initial population guesses". A
`LAM_SCALE_OPT=LIMIT` resubmit was prepared but reverted per user
directive: the clip only masks the divergence - if the system diverges
the cause is elsewhere (J at the reversal depths, the boundary
sigma=-1.005 cell, or the initial populations) and must surface.
Diagnose the divergence, do not damp it.

Non-local Sobolev coupling (19-Sep-2026) - Rybicki & Hummer 1978:

New file subs/sob_rh_couple.f implements the RH78 (ApJ 219, 654)
non-local coupling between distant Sobolev resonance points, in the
geometry-free form of Feldmeier & Nikutta (2006, A&A 446, 661;
astro-ph/0510806).  Motivation: the Dessart floor is an ad-hoc patch;
the physically correct treatment gives photons trapped at sigma ~ -1
a real escape/reprocessing channel at remote resonance points.

Method: for each line l and depth r,

  dJ_l(r) = 0.5 * INT dr' r' v' GAMH_l(r') S_l(r') K(r,r')
  K(r,r') = EXF(t)*EXF(t') * 1/sqrt[(r^2-r'^2)(r^2v'^2-r'^2v^2)]
  t = GAMH_l(r)/|1+mu^2 sigma(r)|, mu^2 from RH Eq.12

applied as dZNET_l = -dJ_l*CHIL_l/(ETAL_l*BB_COR), matching the
SOBJBAR_SIM convention ZNET_l = 1 - Jbar/(S_l*BB_COR) (same
conversion as BETAC_SIM).  Two-point kernel only - 3-point
directions would need the extra exp(-t'') extinction (bounded
overestimate without it).  Kernel uses UNCLAMPED |1+mu^2 sigma|;
the Dessart floor is only for the local beta integral.

VADAT flag: DO_NONLOC_SOB (default F).  Requires DO_NONMONO_SOB=T
(uses SIGMA_TRUE for the resonance geometry); rd_control_variables
warns and disables it otherwise.  Lagged (lambda-iteration)
correction - no Jacobian derivatives added.

Call site: sub_sob_line_v3.f after SOBJBAR_SIM, gated on
DO_NONMONO_SOB .AND. DO_NONLOC_SOB.  Geometry (partner pairs, mu^2,
Upsilon1, weights) is line-independent and cached in SOB_RH_MOD;
rebuilt only when the R/V signature changes.

Verified standalone (test_sob_rh driver): v~r gives 0 pairs
(monotonic LOS velocity along every ray); v~1/r couples all pairs
with dJ>0/ZNET down (correct sign); plateau+shock-rise profile
couples 3376 pairs with the largest correction in the dense rising
region.  Bookkeeping dZNET = dJ*chi/eta/BB_COR verified exactly.

First run: TypeIIn_10d_DEV7_SOB.yaml + DO_NONLOC_SOB=T, job 11820632.
Watch for: "SOB_RH: N non-local coupling pairs" in OUTGEN, whether
the STEQ explosion at sigma~ -1.005 depths is damped, and TUNE cost.

Result (job 11820632): SOB_RH fired - 9964 coupling pairs over 148
depths - but the run still aborted at it11 (~1e107) with the SAME
signature as floor-only: var 193 (Fe-group population eqn) dominant
at all depths, corrections growing ~10 orders/iter.  OBS dep 11820633
auto-cancelled.  Conclusion: two-point non-local coupling is a modest
(few-%) correction to Jbar and cannot de-stiffen net rates that are
catastrophically wrong in the sigma ~ -1 region.  The STEQ runaway
is not rescued by RH78 coupling alone.  Next options: (a) quantify
dJ/Jbar in production (add diagnostic), (b) 3-point coupling +
kink handling per Feldmeier-Nikutta, (c) damped/limited correction,
(d) accept Knop as the non-mono path for this model class.

RH78 revision 2 (19-Sep-2026, job 11821347) - Jacobian-clean split:

Diagnostic run 11820830 revealed the first implementation's real
defect: folding dJ into ZNET_SIM injected the remote term into the
BA Jacobian via ZNET_SIM/POPS(NUP) and (under LAMBDA_ITERATION) via
VB_SIM=-JBAR/ETAL.  The lagged remote term has no local-population
derivative, so both channels were wrong.  Fixed: SOB_RH_COUPLE now
returns DRC_SIM = CHIL*dJ/BB_COR (rate units, eta*ZNET scale) in a
separate line_mod array; sub_sob_line_v3 applies it only in the STEQ
residual T2 = ETAL*ZNET - DRC_SIM.  ZNET_SIM and JBAR untouched.
Verified against FN06 Eq.10-16 (kernel, mu^2, tau, pair cuts all
match); standalone driver re-verified all three velocity cases.

Result: STILL DIVERGES - it11 abort ~1e112, same ~10 orders/iter
growth.  Diagnostics: max dJ/J grows 62 -> 1e13 across iterations
(tracks diverging remote S), max dRATE 4e6 -> 1.7e8.  Conclusion:
the Jacobian artifact was real but not the driver.  The Sobolev
line-rate map itself is anti-contractive here: remote absorption
DRC demands local populations that cannot balance it, so each
lambda iterate amplifies.  Adding physically-correct nonlocal
coupling strengthens, not damps, the channel because dJ follows
the exploding source functions.  Consistent with the Knop+BLANK
control converging on the same hydro - the escape-probability
line formulation is the diverger for this sigma structure.

Files: subs/sob_rh_couple.f (DRC interface), line_mod.f DRC_SIM,
sub_sob_line_v3.f zeroes DRC_SIM per call + T2 subtraction only.
DO_NONLOC_SOB remains default-F.  NOTE: sob_rh_couple.f was twice
found reverted on disk mid-session (stale editor buffer?); if
editing it, verify `grep "SUBROUTINE SOB_RH_COUPLE"` shows the
13-arg (DRC,ZNET_SIM,...) signature before building.

### KNOP_MOM stability fixes + TypeIIn rerun (19-Sep-2026)

Three fixes applied after auditing the TypeIIn_10d_DEV7_KNOP_MOM
divergence (job 11832268, L_outer ~10x/iter to 6e14 by it11):

1. get_moms_non_rel.f - degenerate-J guards on every J-denominator
   ratio.  K/J -> 1/3 (isotropic limit, same threshold and philosophy
   as the FEDD guard in fg_j_cmf_nonmono_v1.f), H/J and N/J -> 0, and
   NMID_ON_J -> 0 where RSQ_JNU is sub-tiny (1.0E-30_LDP).  Unguarded
   quotients at EUV frequencies (JNU_STORE~1e-40) produced garbage
   Eddington factors that QFROMF integrated into grid-wide degenerate
   QF (bounded at 723/724 depths) -> near-singular tridiagonal.
   Verified: rerun reports ZERO "QF_NM bounded" lines (was hundreds).
2. Under-relaxation on lagged GS fields (omega=0.5):
   fg_j_cmf_nonmono_v1.f INEXT_RELAX applied to I_P/I_M_DEPTH_NEXT
   writes; mom_j_cmf_nonmono_v1.f JNEXT_RELAX applied to the
   J_NEXT_FRQ update.  Form: OLD + w*(NEW-OLD); cells with OLD<=0
   fill directly (cold start unaffected).  Fixed point unchanged;
   damps the lagged-feedback gain.  Knop-path only - stock solvers
   untouched.
3. comp_j_blank.f - upstream string bug fixed at both post-solve
   inner-flux blocks: `'ZERO FLUX'` (space) -> `'ZERO_FLUX'`
   (underscore), matching VADAT/solver spelling.  Previously the
   zero-flux branch never matched and HFLUX_AT_IB always used the
   diffusion formula.  Verified: inner luminosity now exactly 0.
   (Bug also exists in cmfgen_reference - flag upstream.)

Infra fix (cmfgen_base/submit_cmfgen.sh, both FINAL+OBS templates):
check_outgen_for_fatal's exclusion list covered the FG GS diagnostic
line but not the MOM one - `MOM_J_CMF_NONMONO_V1 ... capped/NaN
writes= 0 NCF=` matched the `NaN` fatal pattern, so every MOM-mode
job FAILED instead of requeuing at graceful stop (job 11834123,
exit 3:0).  Added exclusion `MOM_J_CMF_NONMONO_V1.*capped/NaN
writes= *0 NCF=` (only the all-clear line; capped>0 still trips).

Result (job 11834123, resubmitted fresh - diverged it11 checkpoints
moved aside first): STILL DIVERGES.  L_outer 2.9e8 -> 1.8e8 (it1-3,
healthy decay) then 6.2e8 -> 3.9e9 -> 3.9e10 -> 5.5e11 -> 5.6e12 ->
7.5e13 -> 1.2e15 -> ... -> 1.5e18 at it15, tracking the unguarded
run almost exactly.  Negative-J resets localize at expanded depths
~192-196 (sigma_true ~ -37, chi ~ 1.7e-5 EUV).  Conclusion: QF
degeneracy was a symptom amplifier; the driver is the compressing-
branch lagged feedback itself - |PSIPREV|*J_NEXT(k+1) injected once
per great iteration with loop gain >1 at sigma<<-1 depths.  omega=0.5
halves the increment but cannot change the sign.
GS max-change stayed ~1e5-1e6 throughout - the lagged field never
converged and warm-up LAMBDA remained forced, so the growth is in
the transfer solve itself, not the BA step.

### PSI_NM compressing-weight fix - divergence cured (19-Sep-2026)

Root cause pinned down by comparing with the formal solver's
compressing source: FG uses SOURCE_PRIME = S + |Q_k|*I_NEXT where the
explicit upwind weight is the SAME local |Q_k| that sits on the
implicit side.  The moment solver instead used
|PSIPREV| = DTAUONQ*|GAM*(1+sigma*K_ON_J_PREV)| - the explicit J_{k+1}
coupling weighted by the k-1 Eddington factor.  The mismatched weight
can exceed the diagonal's |Q_k|, so the lagged fixed-point map
J_next = F(J_next) is EXPANSIVE at strong-compression cells and J
grows ~10x per great iteration.

Fix: at every compressing injection site use PSI_NM
(= DTAUONQ*|Q_k|, already computed for the diagonal) as the J_{k+1}
coefficient instead of ABS(PSIPREV).  Equal explicit/implicit weights
make the RHS a convex combination - J_k <= max(S, J_{k+1}, neighbours)
- so the lagged map is non-expansive and cannot amplify.  Sites:
interior G_ONLY + else branches, both outer-boundary branches, and
the negative-solution rebuild path (5 places; expanding branch still
uses signed PSIPREV, unchanged).  This also mirrors the V11
convention's spirit (velocity factors at k, Eddington at the coupled
frequency) while giving the FG-consistent local-|Q| choice.

Result (job 11837622, fresh start): STABLE + CONTRACTING.
L_outer 2.86e8 -> 1.5e8 over its 1-5 and pinned at ~1.5-1.6e8 through
it12+ (prior runs were at 3e9..1e18 by then).  MAXCH bounded ~300-1100
at the shock-ramp depths - the same limit-cycle regime as the
USE_LAM_ES=T formal run.  GS max-change decaying (1e5 -> 5e4);
warm-up LAMBDA still forced while the GS field matures, so the BA
phase with the diagonal-ALO linearization is yet to be exercised on
this model.  QF_NM bounded count remains 0; inner L = 0 (ZERO_FLUX).

Also fixed earlier today (same session): get_moms_non_rel degenerate-J
guards (K/J->1/3, H/J,N/J->0, NMID_ON_J->0 at JNU<1e-30);
INEXT_RELAX/JNEXT_RELAX=0.5 under-relaxation on the lagged writes;
comp_j_blank 'ZERO FLUX'->'ZERO_FLUX' upstream typo; and the
submit_cmfgen.sh check_outgen_for_fatal exclusion for the
MOM_J_CMF_NONMONO_V1 'capped/NaN writes= 0 NCF=' line which had been
failing every MOM-mode requeue.

### In-run GS-immaturity gate ported from dev6 (19-Sep-2026, late)

Convergence audit of the TypeIIn KNOP_MOM run (11837622) vs Knop et al.
2009 Papers I+II (docs/Knop_et_al_2009_formal_solution.pdf,
docs/aa11793-09.pdf):

- Observed limit cycle: DO_LAM_AUTO releases RD_LAMBDA when MAXCH<50,
  then NUM_LAM=2 alternates 2 LAMBDA + 1 BA.  Every BA step produced
  1e4-1e8% corrections (it87 1.5e5, it102 6.4e5, it119 6.9e7) against
  a GS lagged field still moving by rel. change ~1e3-3e4 per sweep
  (~10-24M cells >50%/sweep).  MAXCH>1e5 then forces LAMBDA again -
  permanent oscillation, never converges.
- dev6 had the matching mechanism (dropped in the dev7 clean restart):
  an in-run gate in SOLVE_FOR_POPS_NONMONO_V1 that forced LAMBDA
  whenever the last GS sweep was still moving >50%, recorded cap
  events, or had significant negative-J resets (GS_UNSAFE_TOL=0.5).
- Ported minimally into solve_for_pops.f (DO_NONMONO_KNOP-gated only):
  GET_FG_GS_MATURITY now also returns the live FG_GS_MAXCHG and
  FG_INEXT_CAPPED (5 args); the MOM getter already returns
  MAXCHG/CAPPED/NEGJ.  When a BA iteration is pending and the GS
  field is immature the iteration is forced to LAMBDA (DIAG,
  MAX_LAM_COR) and FORCE_LAMBDA keeps the mode controller from
  re-requesting BA while the field stays immature.
- Also ported two dev6 audit fixes: NEWMOD resets the SAVE'd scalars
  MAXCH_SUM/COUNT_FULL_LAM_BA, and R_OLD is initialised before the
  REVISE_R_GRID branch (DO_NG_BAND_ACCEL_V4 reads it unconditionally
  when NG fires with REV_RGRID=F - previously uninitialised memory).

Paper-vs-dev7 structural difference (documented, not a bug): the paper
converges the wavelength boundary-value problem to ~1e-8 INSIDE each
formal solve for a fixed source function, then iterates the ALI with
the exact wavelength-tridiagonal Lambda* operator.  dev7 merges the GS
frequency coupling into the SE great-iteration loop (one sweep per
population update) and uses only the diagonal ALO.  The in-run gate
keeps the population solver in the contraction regime (LAMBDA) until
the merged loop matures; whether pure-LAMBDA reaches <0.1% on this
sigma~-70 model is the open question the restarted run must answer.

Build: both files recompiled clean, strings in exe/cmfgen_dev.exe
verified (23:04 build).  Job 11837622 requeued at 23:01:57 just before
the relink, so it picks up the gated binary at its next segment
boundary (~2h); no action needed.

### Frozen-population inner sweeps - Knop Paper-II nested iteration (19-Sep-2026, late)

Second surgical piece closing the Paper-II convergence setup.  The
paper converges the wavelength boundary-value problem internally for
a fixed source function before each ALI/source update; dev7 instead
merged the GS coupling into the SE loop (one sweep per population
update).  Rather than wrapping the fragile ML frequency loop (line-
storage bookkeeping, STEQ/BA accumulators are not re-entrant), the
nested structure is recovered at the great-iteration level:

- In SOLVE_FOR_POPS, when the just-completed sweep's lagged field is
  immature (FG or MOM max rel. change > GS_MATURE_TOL=5%, any capped
  write, or MOM significant-negative-J resets), SOLVEBA_V13 is not
  called - POPS are untouched, so the next great iteration is a pure
  transfer sweep that converges the coupling for FIXED populations.
  The SE update fires automatically on the first mature sweep, which
  is exactly the paper's ordering (converged BVP -> source update).
- Frozen iterations report MAXCH = 100*max(FG,MOM GS change) so
  CONVERGENCE_HISTORY tracks the inner iteration; they stay in
  LAMBDA mode via FORCE_LAMBDA and cannot fire NG/AV/smoothing/
  grid-revision (all gated on .NOT. GS_FREEZE).
- N_FROZEN_MAX=15 caps the freeze: if the lagged map cannot contract
  in 15 sweeps the update proceeds anyway (degrades to the previous
  merged-iteration behaviour, warning printed).
- Checkpoints continue to be written every frozen iteration
  (SCR_RITE_V2 -> POINT1/2 + FGINEXT_GS/MOMJNXT_GS), so the freeze
  state is requeue-safe; the cap counter N_FROZEN resets on NEWMOD.

Evaluation run: TypeIIn_10d_DEV7_KNOP_MOM resubmitted via
run_model_grid.py on the debug partition (FINAL 11853189 / OBS
11853190), cold start.  Watch for 'SOLVE_FOR_POPS: KNOP inner-GS
sweep' + 'frozen iter' lines and whether FG/MOM GS max-change
contracts across frozen sweeps; the first real SE update should
follow the first mature (<5%) sweep.

### Frozen-sweep cap removed (20-Sep-2026)

N_FROZEN_MAX=15 deleted from SOLVE_FOR_POPS; the freeze now persists
until the lagged field is genuinely mature.  Rationale: the first
FIX_T=F run (TypeIIn_10d_DEV7_KNOP_MOM_TFIT, job 11869230) reached
the cap at it16 with the field still moving ~1e6-1e9; the forced SE
update kicked FG/MOM max-change to ~1e16-1e17 and the next formal
solve died in comp_j_blank ("Mean intensity blowing up").  The cap
reintroduced the BA-on-immature-field failure on a timer.  N_FROZEN
is retained purely as the 'frozen iter N' diagnostic counter.
Consequence: a non-contracting GS map now stalls (frozen forever)
rather than degrading - watch consecutive segments for that failure
mode before considering a maturity-metric change (e.g. gating on the
cells>0.5 bulk fraction rather than single-cell max).

### Per-ray/per-frequency Knop refresh + asymmetric spacing (20-Sep-2026, late)

Structural fix for the one-bin-per-great-iteration lag that kept the
GS field immature forever (TypeIIn: FG maxchg ~3e4, MOM ~1e6, frozen
sweeps indefinitely).  Instead of waiting NCF great iterations for
red-to-blue information to cross a compressing zone, each sweep now
captures its transfer inputs and an internal alternating-direction
re-solve converges the lagged coupling before the populations are
updated.  Paper-II ordering is preserved: converged BVP field, then
SE update via the existing maturity gate (unchanged).

FG path (fg_j_cmf_nonmono_v1.f):
- Per-frequency capture during the ML sweep: KN_CHI_EXT/KN_ETA_EXT
  (single precision, extended depth grid) plus scalar/grid data.
- FG_KNOP_REFRESH_V1 (called from cmfgen_sub after 10000 CONTINUE)
  re-solves each ray's (space,frequency) intensity field in
  alternating direction passes on the captured inputs, using the same
  short-characteristics discretisation (DTAU, E-functions, A0-A4,
  monotonic slopes, boundary conditions).  Writes the converged field
  back into I_P/I_M_DEPTH_NEXT and updates the GS maturity metrics.
  Guards: DO_KNOP_REFRESH, HAS_COMPRESSING, complete sweep, INTEGRAL
  method only (TRIDIAG path is not refreshed - known limitation).
- VADAT switch DO_KNOP_REFRESH (default DO_NONMONO_KNOP) via
  control_variable_mod / rd_control_variables / SET_FG_KNOP_REFRESH.
- SET_FG_NU_NEXT(NU_NEXT_FREQ) passes nu(k+1) per call so compressing
  cells use the forward spacing ln(nu_k/nu_{k+1}) instead of the
  backward ln(nu_{k-1}/nu_k); falls back to dLOG_NU when the next
  frequency is unavailable.  MOM solver got the same fix via
  SET_MOM_NU_NEXT and per-depth spacing in mom_j_cmf_nonmono_v1.f
  (midpoint sign held in GAMH_SIGNED_NM).

MOM path (mom_j_cmf_nonmono_v1.f, USE_LAM_ES=F):
- KN_M_* capture store: per-frequency ETA/CHI/ESEC (SM grid), the
  post-GET_MOMS_NON_REL Eddington factors K_ON_J/NMID_ON_J/
  NMID_ON_HMID, boundary data (HBC/NBC/IN_HBC, JPLUS_IB/JMIN_IB,
  HNU_AT_IB), DBB/IC, dLOG_NU, nu(k+1); sweep-fixed args captured once.
- MOM_KNOP_REFRESH_V1 (cmfgen_sub, gated on DO_KNOP_REFRESH) re-drives
  MOM_J_CMF_NONMONO_V1 itself in alternating frequency order, so the
  tridiagonal assembly/THOMAS solve/safeguards are identical by
  construction.  MOM_KN_REFRESHING suppresses GET_MOMS_NON_REL (whose
  EDDFACTOR inputs exist only for the last sweep frequency) and the
  capture writes.  The lagged bluer field is seeded from
  J_NEXT_FRQ(:,L-1); lagged H from a scratch store.  Pass metrics are
  written back into MOM_GS_MAXCHG/MOM_INEXT_CAPPED/MOM_NEGJ_LRG so the
  frozen-iteration gate sees the refreshed field's self-consistency.
- Self-guards: no-op unless a complete KN_M store exists, some depth
  was compressing this sweep, and J_NEXT_VALID.

Pipeline (cmfgen_base/submit_cmfgen.sh):
- On requeue (file_requeue_count>0) with both EDDFACTOR and
  EDDFACTOR_INFO present, the leading COMP_F token in VADAT is flipped
  T->F after prepare_VADAT.py runs, so restarted segments reuse the
  persisted Eddington factors instead of re-Lambda-iterating them
  (previously each requeue reset the scattering J to zero at the
  first frequency).

Build: full `make cmfgen` clean; exe/cmfgen_dev.exe relinked.
Validation: TypeIIn_10d_DEV7_KNOP_MOM_TFIT (interactive, job
11877980 + requeue 11877981) was already running on the pre-refresh
binary; its next segment loads the new build.  Watch for
'FG_KNOP_REFRESH pass' / 'MOM_KNOP_REFRESH pass' lines and whether
FG/MOM GS max-change drops below the 5% maturity gate within the
refresh passes, then whether the first real SE update fires.

FREQ_IDX_ITER off-by-one (20-Sep-2026, later): during frequency k of
the sweep, FREQ_IDX_ITER is k-1 - the INIT branch (FIRST_FREQ=T at
ML=1) resets it to 0 and skips the NEW_FREQ increment.  So the lagged
I_P/I_M_DEPTH_NEXT store actually holds frequency s+1 in slot s
(slot NCF is never written = red boundary zeros; frequency 1 is
never stored).  The first refresh implementation assumed slot s =
freq s, so the capture armed at LK=FREQ_IDX_ITER never reached
KN_NCF, KN_STORE_READY stayed FALSE, and FG_KNOP_REFRESH_V1 returned
silently before its warning checks.  Fixed by: capture at
LK=FREQ_IDX_ITER+1 (KN_* slot = true frequency), SEED_FLD reading
STORE(:,:,l-1), WRITE_STORE writing STORE(:,:,l-1), and the driver
completeness check comparing FREQ_IDX_ITER to KN_NCF-1.  The
FGINEXT_GS checkpoint layout is unchanged (it was always shifted and
is self-consistent).  The MOM path uses FREQ_INDX (=ML) directly and
is unaffected.  Separately fixed the same day: COMP_JCONT_V4.INC
(the OBS-path call site) now calls SET_FG_NU_NEXT(NU(FREQ_INDX+1))
so OBS frequencies do not inherit the stale last-sweep nu(k+1).

### Refresh residual metric + convergence behaviour (20-Sep-2026, late night)

First live refresh runs (jobs 11880567/11881301/11881901, MOM_TFIT
11877980 requeue) exposed three issues the paper does not cover:

1. Faint-cell dominated residual.  |dV|/max(V,1e-25) gives res ~1e21
   when a cell leaves the numerical floor.  Residual is now normalized
   against the per-(ray,frequency) column peak (VMAX_L), so only cells
   that carry a significant fraction of the local intensity count.
   MOM refresh metric got the same treatment (KN_M store column peak).

2. Limit cycle / marginal contraction.  The undamped alternating GS
   converged fast (2e8 -> 0.66 in 6 passes) then either entered a
   period-4 cycle at ~14% residual (cold start) or crawled at ~1%/pass
   (warm restart, lambda ~ 0.99 eigenmode).  Refresh updates are now
   under/over-relaxed per ray: omega adapts between ~0.5 (residual
   growing/oscillating) and ~1.3 (smooth contraction); the fixed point
   is unchanged.  Pass history is printed (FG_KNOP_REFRESH: pass
   residuals).

3. Gate frozen by outliers.  GS_MATURE_TOL on the max-norm meant a
   handful of marginal cells (|Q|->1 at strong compressions) could
   keep the model in transfer-only mode forever.  FG_GS_NBAD/FG_GS_NTOT
   now count significant cells (>1e-4 of column peak) still above
   tolerance; solve_for_pops treats the field as immature only when
   the bad fraction exceeds ~0.1% (or a hard runaway threshold, or
   cap/NaN events).  Counters reset each sweep so a skipped refresh
   cannot open the gate.

Also: FREEZE_I_NEXT (JFNK freeze) now blocks the FG refresh writeback,
matching the in-sweep write guard.

Known open item: MOM_KNOP_REFRESH_V1 printed no output on the
MOM_TFIT job despite the MOM nonmono solver running - a silent guard
return.  A one-shot guard-state diagnostic line
('MOM_KNOP_REFRESH guard: ...') was added at driver entry to identify
which condition fails on the next segment.

### MOM refresh guard + SOLVEBA-on-immature-field fixes (20-Sep-2026, late night)

Guard diagnostic on job 11885512 showed WHY MOM_KNOP_REFRESH_V1 never
ran: KN_M_ANYCMP stayed FALSE all sweep.  HAS_COMPRESSING was only
tested on the midpoint factor (1 + AV_SIGMA_TRUE*NMID_ON_HMID), but
NMID_ON_HMID is floored at 0.01 in get_moms_rel.f, so the midpoint
factor stays ~0.02 > 0 even at sigma_true ~ -98.  The actual
compressing branch in the tridiagonal keys off the DEPTH-point
Q_NM(I)<0 (1 + SIGMA_TRUE*K_ON_J, K_ON_J~0.3), which does go negative
-- the coupling was active but invisible to the flag.  Fixed:
HAS_COMPRESSING now also set on depth-point GAM_SIGNED_NM<0.

Restored original-CMFGEN iteration semantics in solve_for_pops.f:
GS_FREEZE no longer skips SOLVEBA_V13.  Immature-field iterations now
run a forced-LAMBDA update (LAMBDA_ITERATION=T, FIXED_T=T,
COMPUTE_BA=T, TEMP_CHAR='DIAG', T1=MAX_LAM_COR), matching V6's
solve_for_pops_nonmono_v1 which always ran SOLVEBA_NONMONO_V1.  All
standard iteration diagnostics (LAMBDA iteration used, max %
increase/decrease, DEC_VEC/INC_VEC, MINSCALE) print every iteration.
The .NOT. GS_FREEZE gates on AUTO_SMOOTH_POPS / REVISE_R_GRID /
UNDO_IT / AVE_FLIPS remain -- conservative while the field matures.
N_FROZEN renamed in print to 'immature-field iter' (counter kept).

FG refresh adaptive relaxation revised after job 11881901 showed
omega=0.6 under-relaxation HURT (sawtooth envelopes converging ~10%
per direction pair vs undamped reaching 0.13): OMEGA_RAY now starts
at 1.0 and adapts on SAME-DIRECTION residual history (pass N vs N-2,
since alternating sweeps interleave expanding/compressing passes with
intrinsically different residual levels).  Same-direction growth ->
damp toward 0.45; slow same-direction contraction (ratio>0.85) ->
over-relax toward 1.45.  Also added bulk early-exit after pass>=6
when NBAD <= max(50,NTOT/1000), matching the maturity-gate level.

Both running jobs (11881901 timed out mid-it2 at the wall; MOM_TFIT
11877980 timed out at wall) resubmitted warm: FG=11885796,
MOM_TFIT=11886410 (cancelled 11885512/13 mid-sweep to pick up the
HAS_COMPRESSING fix immediately rather than waiting ~1.7h for the
requeue boundary).

### Single-precision capture-store NaN bug + int32 counter overflow (21-Sep-2026)

MOM_KNOP_REFRESH_V1 ran on job 11886410 (CMP=T after the
HAS_COMPRESSING fix) but every re-solve returned NaN rows
('non-finite XM_NM zeroed', SOURCE_NM=NaN at NM point 2) while the
in-sweep solve at the same frequencies was clean.  Root cause:
KN_M_ETASM/KN_M_CHISM/KN_M_ESECSM (and the K_ON_J/NMID stores) were
default REAL = single precision.  Perfectly valid double ETA values
below ~1.2E-38 (cold faint outer cells) underflow to exact 0.0 on
capture; the refresh's LOG(ETA_SM) interp then produces -Inf -> NaN
coefficients -> NaN source rows.  Fix: all six KN_M arrays promoted
to REAL(KIND=LDP); memory cost ~15 MB, irrelevant next to the 8 GB
FG store.  A one-shot store-integrity scan
('MOM_KNOP_REFRESH store scan: bad ETA=.. bad CHI=..') now runs at
first refresh entry to confirm zero non-positive entries.

FG side had the same single-precision store (KN_CHI_EXT/KN_ETA_EXT)
but its LOG(MAX(..,1E-300)) floor prevents NaN; the residual
corruption (values <1e-38 replayed as 1e-300) is physically
negligible at cells emitting ~1e-40, so the single-precision store
stays for its documented ~300 MB memory saving.

FG_GS_NTOT printed negative ('unconverged significant cells=
82210160 of -22911603'): int32 overflow.  With NCF=43475 and ~181
rays the last-pass significant-cell count is ~5e10.  Promoted
FG_GS_NBAD/FG_GS_NTOT/FG_GS_NCONV, the NBD_LS/NTT_LS ray arrays, the
worker NBAD/NTOT/NB_*/NT_*/NBAD_CUR/NTOT_CUR counters, the
GET_FG_GS_MATURITY NBAD/NTOT dummies, and the solve_for_pops
GS_NBAD/GS_NTOT receivers to INTEGER*8; print formats widened to
I12/I14 and MAX(50_8,...) used in the bulk comparisons.

Rebuilt cmfgen_dev.exe at 22:59; MOM_TFIT resubmitted warm as
11887664 (OBS 11887665) to exercise the LDP store + scan.

### MOM refresh GS cascade -> Jacobi snapshot (20-Sep-2026, 23:30)

Job 11887664 (LDP-store binary) ran the refresh cleanly -- no NaN
storm (capped counts ~3k/pass vs 724/call) -- but diverged: descending
passes spiked to 1E5-1E10 residuals with ~9E5 significant negative-J
resets per pass, while ascending passes sat at ~1E3.  Cause: on a
descending pass the routine's compressing read J_NEXT_FRQ(:,L+1) sees
a THIS-PASS write -- a true GS cascade along the coupling direction.
The in-sweep read is pure Jacobi (previous great iteration's field),
so the refresh's descending pass was strictly more aggressive than
the production semantics it is meant to converge.

Fix: JNEXT_RF_READ pass-start snapshot (module array).  While
MOM_KN_REFRESHING is set, all J_NEXT_FRQ(:,FREQ_INDX+1) read sites
(2 interior + 2 boundary + 1 negative-J retry) take JUP_NM from the
snapshot; writes still go through JNEXT_RELAX=0.5 under-relaxation.
This makes the refresh's compressing coupling identical to the
in-sweep lagged coupling.  The expanding side (RSQJNU_PREV chain)
is untouched.

MOM counters promoted to INTEGER*8 (MOM_INEXT_CAPPED, MOM_NEGJ_LRG,
MOM_GS_CAP_FILE, getter dummies, PASS_CAP/NEG, solve_for_pops
receivers): the NaN-storm path could overflow int32 (~3E10 events).
MOMJNXT_GS checkpoint VER bumped 2->3; older files are cleanly
rejected and J_NEXT rebuilt via Lambda warm-up.

Binary relinked 23:27; MOM_TFIT resubmitted warm as 11888720.

### MOM refresh containment + FG bad-cell metric fix (21-Sep-2026, 00:10)

Job 11888720 (Jacobi-snapshot binary) showed the snapshot works:
refresh contracted monotonically for 22 passes (8.1E4 -> ~44), then
an expansive mode blew the tail up to ~1E13.  negJ stayed ~9.7E5
cells/pass throughout -- about 3% of the ND x NCF field sign-flips
on EVERY pass, including in-sweep solves (the 200+ MB
MOM_J_NONMONO_ERRORS dump is mostly in-sweep negative-J rows at
sigma_true ~ -37 to -55 compressions).  The lagged fixed-point map
of the moment solve is genuinely non-contractive on that subspace;
it is a discretization-level instability, not a plumbing bug.

Containment added (mom_j_cmf_nonmono_v1.f):
- XM_NM magnitude cap at 1E40 alongside the non-finite zeroing, so a
  runaway-but-finite value cannot overflow to Inf inside the next
  solve's PSI*JUP/EPS*RSQJNU coupling terms.
- Refresh keeps JNEXT_BEST = field at the lowest-residual pass; on
  exit without convergence J_NEXT_FRQ is restored to it, and
  MOM_GS_MAXCHG reports BEST_CHG (the residual of the field actually
  stored).
- Abort once PASS_CHG > 1E3*BEST_CHG (post pass 3) -- the explosive
  tail has no chance of recovering and only wastes wall time.

FG side (fg_j_cmf_nonmono_v1.f): FG_GS_NBAD counted self-normalized
>5% changes on a ray's last pass.  Rays with no compressing cells
exit after pass 1, where the "change" is seed->exact-solve -- a
transient, not an unconverged state.  CLEAN_EXIT flag now zeroes a
ray's NB_O on RES<TOL or !HAS_CP exits, so FG_GS_NBAD reflects only
cells that are still unconverged.  (The 1.27% bad reported on job
11885796 was dominated by this false positive; the refresh itself
converged 8.4E6 -> 1.5E-2 in 9 passes.)

Binary relinked 00:08; MOM_TFIT resubmitted warm as 11889372.
FG job 11885796 picks this up at its next requeue.

### FG metric fix verified + gate opens (21-Sep-2026, 01:45)

Job 11885796 requeued onto the 00:32 binary and the CLEAN_EXIT/RES
fix confirmed: FG refresh max residual dropped from ~1.8E6 (clean-ray
seed-transient false positive) to 1.08E-1, and unconverged
significant cells from 1.1E8 (1.3%) to 1.98E6 (0.023%).  Pass
residuals converge 1.8E6 -> 0.053 in 10 passes under adaptive omega.
With NBAD/NTOT < 0.1% the bulk gate opens; population updates have
collapsed to physical scale (max Lambda change 4.5E+1 vs 9.5E+94
earlier).

MOM job 11889372: refresh containment stable across iterations --
contracts ~1-2 dex for 9-15 passes, expansive mode explodes, abort
fires, best-pass field restored (~1E4-5E2 residual), no NaN
propagation, ~60-80 s/call vs 350 s before.  The 1E40 XM cap is now
catching ~2.5E5 runaway cells/pass that previously became Inf.
The moment lagged-map instability is intrinsic to this model's
compression strengths (same negative-J rows appear in-sweep);
the MOM path stays in forced-Lambda (V6 gate semantics), which is
the correct conservative state.

### FIX_IMPURITY cleared on KNOP forced-LAMBDA sweeps (22-Sep-2026, 20:40)

Audit of FIX_IMP integration on the non-monotonic path (requested for
the TypeIIn_KNOP_TRIDIAG USE_LAM_ES=F run):

- Wiring verified end-to-end: VADAT FIX_IMP -> RD_FIX_IMP ->
  FIX_IMPURITY (BA mode only, via cmfgen_sub init or the mode
  controller) -> MOD_FIX_IMPURITY (solve_for_pops.f) ->
  FIXPOP_IN_BA_V3 inside GENERATE_FULL_MATRIX_V3 inside
  CMF_BLKBAND_V3 inside SOLVEBA_V13.  Semantics: freeze ions/species
  contributing <1E-15 of the population during the BA update.
  new_ba_up_tst.f and subs/fixpop.f are dead code on this path.
- Gap found: the three KNOP forced-LAMBDA gates (warm-up,
  immature-field, GS_FREEZE) do not clear FIX_IMPURITY, and
  MOD_FIX_IMPURITY is snapshotted before the gates run.  A model
  that reached BA (FIX_IMPURITY=T) and then hit an immature sweep
  would run its forced-LAMBDA relax with trace ions frozen --
  breaking the upstream "FIX_IMPURITY only in non-LAMBDA mode"
  invariant.  Same omission exists in the V6 reference
  (solve_for_pops_nonmono_v1.f).
- Fix (solve_for_pops.f, before CALL SOLVEBA_V13): single guard
  `IF(FORCE_LAMBDA) FIX_IMPURITY=.FALSE.; MOD_FIX_IMPURITY=.FALSE.`
  covering all three gates.  Default path unchanged: the reset only
  fires on forced-LAMBDA iterations where upstream semantics already
  require FIX_IMPURITY=F.
- Config note: with DO_LAM_IT=T + DO_LAM_AUTO=T the option is inert
  in practice -- the MAXCH<50 auto-switch does not restore
  FIX_IMPURITY (upstream soft-off near convergence); it can only
  engage via the CNT_LAM overflow path after the auto-switch.

make cmfgen clean; solve_for_pops.o rebuilt, exe/cmfgen_dev.exe
relinked 20:36.

### Iteration-cycle audit: MOD_FIXED_T + NG guard (23-Sep-2026, 08:45)

Audit of the KNOP convergence cycle in solve_for_pops.f against the
V6 reference (solve_for_pops_nonmono_v1.f) and stock upstream:

Fixed (same snapshot-ordering class as FIX_IMPURITY):
- MOD_FIXED_T/MOD_FIX_T_D_ST/MOD_FIX_T_D_END were computed from the
  incoming FIXED_T/RD_FIX_T before the KNOP gates force FIXED_T=T.
  A forced-LAMBDA iteration entered with FIXED_T=F would leave the
  R.E. row unprotected -- degenerate T equation, singular BA, hard
  STOP in SOLVEBA_V13 (or a garbage dT).  Masked so far by FIX_T=T.
  The FORCE_LAMBDA block now also forces MOD_FIXED_T=T with the
  full-depth range (identical to what the next call computes).
- NG acceleration guard had .NOT. GS_FREEZE but not
  .NOT. FORCE_LAMBDA: NEXT_NG is scheduled off LAST_NG, so a queued
  NG kick could land on a warm-up iteration with a mature restored
  field (GS_FREEZE=F) and extrapolate pops across the I_NEXT reinit.
  AV/UNDO are already protected implicitly via LAST_LAMBDA.

Checked, clean: ITS_SINCE_START/N_FROZEN/MAXCH_SUM/COUNT_FULL_LAM_BA
SAVE'd and NEWMOD-reset; LAST_LAMBDA updates on forced iters;
immature gate + freeze gate coverage has no gap; FIXPOP_IN_BA_V3
DIAG-band mechanics correct; LAPACK IFAIL tested in BLKBAND_V3;
IN_ITS re-read restores FIX_IMPURITY/FIXED_T on the DO_LAM_AUTO flip.

Noted divergence from V6 (not changed): a singular BA is a hard STOP
(SOLVEBA_V13) -- V6's per-depth singularity recovery and
N_BA_FAILURES Lambda-fallback were not ported.  SUCCESS is dead
post-call.  STOP surfaces the failure, which fits the
surface-divergence policy; a fallback would trade that for a silent
Lambda iteration on a bad solve.

Cosmetic: the "LAMBDA iteration used"/"BA matrix computed" headers
print before the KNOP gates, so a frozen BA iteration logs the BA
header without the LAMBDA line (the MAXCH tag still says LAMBDA).

make cmfgen clean; solve_for_pops.o rebuilt, exe/cmfgen_dev.exe
relinked 08:42.

### NG acceleration audit (23-Sep-2026, ~09:00)

Full audit of DO_NG_BAND_ACCEL_V4 + NG_MIT_BAND_OPT_V1 +
DO_TRANS_AND_NG_V1 + NGACCEL, and the call-site in solve_for_pops.f.

Verified correct:
- Iterate ordering: RDPOPS slot 1 = newest record (IREC); NGACCEL's
  D0/D1/D2 differences and the A,B least-squares solve are internally
  consistent with that ordering.  Singular-determinant fallback sets
  one coefficient to zero.
- Per-depth acceptance: LOCINC>10.1x or LOCDEC<0.09x or reversed T
  trend flags a bad depth; >3 bad depths cancels the whole NG and
  restores the newest iterate.  Caller commits POPS/R/V/SIGMA only on
  NG_DONE -- a failed NG leaves the state intact.
- Post-NG control: >50% inc AND dec forces a LAMBDA iteration with
  FIXED_T=T and FIX_IMPURITY=F -- the BA-only impurity invariant holds.
- SIGMA_TRUE is recomputed from the accelerated V after NG success
  (DO_NONMONO_SOB/DO_NONMONO_KNOP branch).
- R/V/SIGMA extrapolation is a no-op on this model: REV_RGRID=F and
  DO_HYDRO=F give WRITE_RVSIG=F, so all 4 history slots carry the
  current R/V/SIGMA and the commit writes back identical values.
  Hydro provenance is preserved.
- Scheduling: NG needs MAXCH < VAL_DO_NG (BEG_NG=10%), >= IBEG_NG=15
  iterations, >=4 iterations since the last NG (ITS_PER_NG=8 spacing),
  re-suspends (NEXT_NG=1500) if MAXCH climbs again, retries +8 after
  failure.  Call-site also excludes LST_ITERATION, GS_FREEZE and now
  FORCE_LAMBDA.

Fixed (upstream bugs in the CHK_NG validity check, dormant under
FIX_T=T because the block only runs when accelerated T differs):
- MEAN was never initialised before MEAN=MEAN+RATIO accumulation --
  the global T-consistency verdict was garbage whenever it ran.
- RATIO=(BIG_POPS-POPS)/T1 divided by T1=RDPOPS(NT,I,1)-RDPOPS(NT,I,2)
  with no guard: stationary T at any depth gives Inf/NaN contaminating
  MEAN/MIN_RATIO/MAX_RATIO.  Depths with T1=0 are now skipped and the
  mean is taken over the counted subset.

Findings (design-level, not patched):
- NG fits a single (A,B) pair over a 60-depth x all-levels band
  assuming one linearly-convergent map.  Under KNOP the per-iteration
  map includes the MOM/FG I_NEXT refresh; the 4-record history can
  span a LAMBDA->BA switch, and on a continuing model
  (NEXT_NG=2000 path) NG can fire on the first iteration of a
  restarted process.  The GS tag mechanism keeps I_NEXT consistent
  across restarts, so this is a narrow edge, not a violation.
- Per-depth bounds (10x/11x) are loose: the observed +87%/-413%
  kicks at iters 57/68 of the ES run passed legally.  With FIX_T=T
  the T-sign and CHK_NG checks are both dormant, so the magnitude
  bound is the only pops guard.  This is why DO_NG=F remains the
  right setting for this model: the kicks are "working as designed",
  not a bug -- the design just tolerates 5x population perturbations.

make cmfgen clean; do_ng_accel_band_v4.o rebuilt, exe/cmfgen_dev.exe
relinked 08:56.

### Full iteration-cycle audit (23-Sep-2026, ~09:30)

Traced the whole great-iteration in cmfgen_sub.f (loop label 20000)
plus solve_for_pops.f, solveba_v13.f, cmf_blkband_v3.f,
generate_full_matrix_v3.f, fixpop_in_ba_v3.f, fiddle_pop_corrections_v2.f,
specify_it_cycle_v3.f, comp_j_blank.f.

Cycle per iteration (verified):
  opacities -> per-freq inner f-loop [FG_J_CMF_NONMONO_V1 (rays->f,Jcmf)
  <-> MOM_J_CMF_NONMONO_V1 (J), or USE_LAM_ES=T: RJ=TC direct from rays]
  until |d f|<ACC_EDD_FAC or 20 iters -> EVALSE_QWVJ_V8 etc. accumulate
  SE residuals + BA Jacobian on the converged RJ -> STEQ_T (RE), STEQ_ED
  (charge cons.) -> FG/MOM_KNOP_REFRESH (mature lagged I_NEXT on this
  sweep's opacities) -> two-photon/chg-exch/advection/adiabatic terms ->
  SOLVE_FOR_POPS [KNOP gates -> FORCE_LAMBDA fixups -> SOLVEBA_V13 ->
  CMF_BLKBAND_V3 -> GENERATE_FULL_MATRIX_V3+FIXPOP -> DGEEQU/DGETRF/
  DGETRS -> POPS *= (1-STEQ*SCALE) MAJOR per-depth limiter] -> POPTOION,
  SUP_TO_FULL, hydro/clump (+SIGMA_TRUE refresh) -> COMP_J_CONV ->
  controller (DO_T_AUTO, MAXCH<EPS, DO_LAM_AUTO, IT_SPECIFIER, IN_ITS).

Maths verified correct:
- Newton step: A.x=F(residual), x converted absolute->fractional,
  multiplicative update POPS*(1-x*SCALE), positivity guaranteed by
  BIG_LIM<1. MAJOR limiter: scale set only by levels >1e-10*Ne,
  T capped at MAX_dT_COR=0.2 and floored at T_MIN=0.3.
- SE block: level SE + ionization eqs + species number conservation
  (EQ_SPECIES rows) + charge conservation (NT-1) + RE (NT). Correct.
- Fixed-point consistency: converged state satisfies SE given J(n),
  J satisfies transfer given n. SE rates use the inner-loop-converged
  RJ, not a stale field.
- KNOP maturity gates measure the self-consistency error of the field
  the SE step actually used (post-sweep refresh residual). Ordering
  sound: gate -> forced-Lambda or BA -> SOLVEBA always runs.
- Convergence (MAXCH<EPS=1%) requires LAST_LAMBDA != MAIN_COUNTER --
  cannot trigger during warm-up/freeze. NUM_ITS>=1 floor, +2 after NG.
- SIGMA_TRUE recomputed after every R/V change (grid revise, NG, hydro).
- SET_POPS_D2_EQ_D1 (T here): STEQ(:,2) zeroed pre-MAXCH, pops slaved
  to depth 1 post-update. Consistent.

Issues found:
1. STALE COMMENT, solve_for_pops.f:84-92 -- claims GS_FREEZE "skips the
   SOLVEBA update entirely ... transfer-only sweep for FIXED pops".
   Actual code forces a LAMBDA iteration and pops DO update (matching
   the later comment and V6 policy). On a divergent field the Lambda
   step still moves pops (TRIDIAG run's 1e24% corrections happened on
   forced-Lambda iters). A true skip-SOLVEBA freeze is NOT implemented.
   Comment is misleading; code matches V6.
2. MAX_CHNG_LIM=1.0D100 in VADAT -- the divergence guard (STOP when
   MAXCH>MAX_CHNG_LIM) is effectively disabled; nothing aborts a
   runaway. Config choice, but explains the TRIDIAG run continuing
   through 1e24% corrections.
3. MAXCH masking: when DO_LEVEL_CHK fires (all but top-10 inc + top-10
   dec < 10%), MAXCH reports the 10th-largest change. Up to 10 levels
   can still move >10% while MAXCH reads <1%; the EPS convergence test
   uses this masked value. Upstream design (ignore few bad levels) --
   flag when interpreting "converged".
4. FIX_NE=T: charge-conservation row zeroed (FIXPOP), Ne frozen at the
   initial value for the whole run (ED read from POPS(NT-1) each iter,
   never updated). Converged ionization is consistent with the frozen
   Ne, not necessarily with sum(z_i n_i). Config choice; matters if
   ionization moves far from the initial ED.
5. FIX_T=T + DO_T_AUTO=F: radiative equilibrium never enforced -- T is
   the input structure for the whole run. SE-converged, not
   energy-conserving. Deliberate per config.
6. BA under KNOP is block-diagonal quasi-Newton: TX/TVX radiation
   response zeroed by design, advection off-diagonals off (INC_ADV=F).
   The "full" iteration has no depth coupling -- convergence relies
   entirely on the outer J<->n map contracting (fails for the MOM
   path under strong compression, works for the ES path). TRIDIAG
   identical to DIAG in this configuration.
7. Inner f-iteration bails at 20 iters with only an OUTGEN warning --
   a non-converged f is used silently for the rest of that frequency.
8. R_GRID_REVISED forces LAMBDA_ITERATION but does not clear
   FIX_IMPURITY (pre-existing upstream edge; REV_RGRID=F here).
9. Jacobian reuse (N_FIX_BA=2, FIX_BA=1%): after MAXCH<1% the BA can
   freeze for >=2 iters; corrections then use the ORIG_POPS fraction
   conversion -- an approximation bounded by the MAXCH_SUM>3*FIX_BA
   recompute trigger.

No source changes made -- all findings are documentation-level or
config-level, not code bugs in the cycle itself.

### KNOP refresh performance pass (24-Sep-2026, commits 41083d9 + d0eb12f)

Audit of FG_KNOP_REFRESH_V1 / KNOP_RAY_SOLVE_V1 found the alternating-
direction fixed-point sweep rebuilt BUILD_FREQ outputs every pass even
though they are pass-invariant (depend only on captured opacities and
static ray geometry, never on the iterated field).  Implemented:

- 41083d9: pass-invariant coefficient pack (DTAU, A0-A4, SR, Q, inner-BC
  ETA/CHI at NI) built once per ray and replayed each pass; guarded by
  KN_RAY_PACK_MB (2 GB/thread), exact fallback to per-pass BUILD_FREQ.
- d0eb12f: DLN_N_SAV hoists the per-frequency LOG(FREQ/FG_NU_NEXT) out of
  three ray-point loops in the in-sweep INTEGRAL path.

Runtime validation (TypeIIn_10d_DEV7_KNOP_ES, FINAL job 11962582, 4 OMP
threads, cold start): refresh t= 749.5/711.3/660.9/646.2/633.2/677.1/
748.7 s over its 1-7 vs 690-773 s pre-change baseline; great-iteration
931-976 s vs 955-1112 s.  Numerics clean (0 diverged rays, 0 clamp
events, 0 cap/bad writes; per-pass residual history matches pre-change
behaviour).  Net refresh gain only ~3-8%, NOT the predicted ~3x:
profiling shows all 4 threads saturated and no compute hotspot -- the
sweeps are memory/latency-bound, so removing BUILD_FREQ flops helps
little.  Kept anyway: verified bit-identical values, small real gain.

Largest remaining lever (not exercised): the refresh is embarrassingly
parallel over KN_NP rays -- raising cpus_per_task (4 -> 16) should cut
it ~3-4x.  Secondary candidates, each needing its own validation: MOM-
path equivalent of the pack (mom_j_cmf_nonmono_v1.f rebuilds+refactors
the tridiagonal per pass); SEED_FLD/WRITE_STORE iterate over the slowest
store dimension (~5 MB strided walk per element group); KN_CHI_COEF/
KN_ETA_COEF (~1 GB double) could be single precision.

### USE_LAM_ES=F maturity fix: quarantine + sweep-drift gate (25-Sep-2026)

The TypeIIn 10d F path (TRIDIAG family) converged in Lambda but the
BA gate could never open: the original test (MOM_GS_MAXCHG<0.05 AND
MOM_NEGJ_LRG==0) is unreachable on this model because the MOM refresh
residual is pinned ~1.3-1.7 by marginally contractive compressing
cells and ~750-2000 significant negJ repairs occur per sweep even in
a settling field.  Three changes, all VADAT-gated (commits 28b0cdd,
7b1e570, 069ed05):

- MOM_RF_KEEPNEG (def F): repaired cells (negJ reset, non-finite,
  cap) keep their prior iterate in J_NEXT_FRQ instead of storing the
  artificial repair.  The repair values are not fixed-point values
  and were reseeding the lagged read with noise: sweep negJ fell
  ~9E5 -> ~2E3.  Quarantined cells are excluded from the residual
  metrics and counted separately (MOM_GS_NQUAR); the write loop
  skips them.
- MOM_RF_GS_DESC (def F): descending refresh passes read the live
  J_NEXT_FRQ(:,L+1) instead of the pass-start JNEXT_RF_READ snapshot
  (true GS along the red-to-blue compressing characteristic).  Tested
  on gsdesc (11966433): it1 best residual 40.8 vs baseline 13.9,
  later iterations track baseline -- neutral-to-worse, left off.
- MOM_KNOP_REFRESH stall exit (unconditional): same-direction
  residual ratio >0.90 twice after pass 8 exits early; ~420 s of
  oscillating passes cut to ~70-100 s/iter.
- MOM_GATE_BULK (def F) + MOM_SW_* snapshot: the bulk counters
  (nbad5/nbad50/ntot/negJ/quar) are now taken at refresh ENTRY, i.e.
  the sweep-level drift |XM_solve - stored J_NEXT| -- the change the
  populations actually see.  The refresh-pass residual measures map
  self-consistency and asymptotes; gating on it can never open BA.
  Freeze while >10% of cells drift >5% per sweep; unsafe on capped
  writes, >0.1% of cells drifting >50%, any cell >10x stored,
  negJ/quar >10%, or refresh residual >1e2.  MOM_SW_* is reset at
  each sweep start (guarded against the refresh's per-frequency
  resets) and snapshotted before the pass loop.

Result (gatebulk 11969828, small partition): freeze released at
it25 when nbad5_sw crossed 10%; first BA computed it28; thereafter
the same Lambda<->BA alternation the T path shows (BA kick ~1e4 ->
1-2 frozen Lambda its -> reopen).  Corrections pre-gate-open were
contracting 17->13%.  Containment throughout: abort guard + full
pre-refresh field restore on unconverged refresh.

Caveat for restarts: resubmitting an immature run replays the
4-iteration warm-up and rebuilds J_NEXT (sweep drift resets to
~2e5); populations carry over via *_IN files.  Only resubmit when
the binary/config change is worth replaying the maturation.

Winning config keys (now in TypeIIn_10d_DEV7_KNOP_TRIDIAG.yaml):
MOM_RF_KEEPNEG=T, MOM_GATE_BULK=T, MOM_RF_GS_DESC=F.

### DO_NONMONO_SOB working setup on Dessart m15 2.82d (27-Sep-2026)

Config `configs/dessart_m15_2p82d_SOB.yaml`, model dir
`cmfgen_validation/m15_2p82d_knop_SOB/cmfgen/2.82days`.
Job 12001805 (interactive, 8 cpu): stable Lambda contraction from a
cold LTE start -- max % increase 1.6e11 (it1) -> 64 (it5) -> 17 (it10)
-> 4.4 (it24), luminosity pinned at L=4.266e8 inner / ~1.04e10 outer
every iteration, DTDR=1.7e-3 nonzero.  Same binary and same
SN_HYDRO_DATA (Luc's 70-point m15 file, sigma=-1.037 outer, ~-157 at
the interface) that diverged 1e52 -> 3e113 -> NaN segfault at it11
under the old config.

The load-bearing fixes, taken from Luc's extracted converged VADAT
(`cmfgen_validation/jnk_extracted/VADAT`) and the Dessart+Audit+
Hillier 2015 setup (arXiv:1503.05463, sec 3.1: inner tau~30 diffusion
boundary, lines treated purely via Sobolev):

- LSTAR = SN_HYDRO_DATA header luminosity (4.265857e8 Lsun).  With
  LSTAR=0 the inner boundary is dead (DTDR=0 every iteration) and the
  statistical-equilibrium iteration has no flux anchor -- this alone
  turned divergence into contraction.
- DIF=T + IB_METH=DIFFUSION.  Inner tau_ross~37 sits in the paper's
  tau~30 diffusion regime.  Side effect: JBAR_W_EDD=T is legal again
  (PERTJFEAU_IBC hard-stops under DIF=F; the earlier JBAR_W_EDD=F /
  NEWPERTJD workaround is no longer needed).
- Solver details matching Luc: ALLOW_OL=F, DC_METH=LTE, THK_LINE=T,
  DO_DIS=T, INC_PEN=T, LAM_VAL=200, FIX_BA=10, N_FIX_BA=2,
  T_INIT_TAU=1.5, SIG_GAU_KMS=3000 (NEG_OPAC_OPT stays forced SRCE_CHK).
- INC_CHG must stay F: CHG_EXCH_DATA is not in the repo and
  prepare_VADAT refuses INC_CHG=T without it.
- Kept VEL_LAW=11 + RD_IN_R_GRID=T (stratified composition preserved).
  Luc used VEL_LAW=12 + RD_IN_R_GRID=F, which resamples onto a
  log-tau_es grid but switches to uniform VADAT abundances -- a
  different physical model; only worth it as a control experiment.
- DO_NONLOC_SOB=F (baseline paper comparison; RH78 coupling is a
  separate experiment and did not cure the TypeIIn divergence anyway).

Forced-by-orchestrator (unchanged): GLOBAL_LINE=SOB, TRANS_*=SOB,
DO_SOB_LINES=T, SOB_FREQ_IN_OBS=T, NON_MONO_SOLVER=T, USE_LAM_ES=T,
SOL_METH=DIAG, NUM_BNDS=1, VEL_LAW=11, RD_IN_R_GRID=T, N_OB_INS=0,
N_IB_INS=0, INC_GRID=F, FG_OPT=INT/INS, FIX_T=T, DO_LAM_IT=T,
DO_LAM_AUTO=T, DO_NG=T, COMP_BA=T, DO_NOMONO_KNOP=F.

Activation still requires BOTH DO_NOMONO_SOB=T and DO_NONMONO_SOB=T in
vadat_overrides; `non_mono_mode: SOB` alone does not open the
non-monotonic block.

Caveat: run still in progress at note time (it24 of up to ~140 inside
the 1:59:00 interactive wall; requeue configured).  The DO_LAM_AUTO BA
transition and final convergence are not yet verified on this run.

### Build provenance + VEL_LAW=12 SIGMA fix (29-Sep-2026)

External review of the dev7 source found a pre-existing numerical bug
and two clean-build failures; all three are fixed and committed:

- `f2f4784` -- `new_main/subs/rv_sn_model_sniin.f`: the interior
  SIGMA estimate divided the weighted average of the two one-sided
  dV/dR slopes by `R(I-1)-R(I)` instead of the full span
  `R(I-1)-R(I+1)`, overestimating `|dV/dR|` by up to 2x on a uniform
  grid and biasing `SIGMA = R/V*dV/dR - 1`.  Reached only under
  `VEL_LAW=12` (cmfgen_sub.f:1044 dispatch); the VEL_LAW=11 path
  (`SET_RV_HYDRO_MODEL_V3` + `adjust_r_grid_v4.f`) uses the correct
  log-difference form and was never affected.  NOT flag-gated --
  the formula is wrong for any smooth V(R).  Every validated run in
  this file used VEL_LAW=11, so no prior result is corrupted; any
  VEL_LAW=12 run predating f2f4784 carries the biased SIGMA.
- `12f1f1a` -- Makefile fixes: all recursive makes now run via
  `$(MAKE) -C dir` (non-zero exit on sub-make failure; the old
  `(cd dir; make)` chains let scripts report success after a failed
  compile).  `d_lte` restored to the `cmfgen` target (link needs
  lte_hydro objects; fc12e0f's removal only worked with stale
  archive members).  `mod_ray_mom_store.f` built as a prerequisite
  of `comp_j_blank.o` in new_main/mod_subs (KNOP dispatch imports
  HNU_AT_OB before plane/ is compiled).
- `22deecd` -- `build.sh` / `build_debug.sh` wrappers
  (set -euo pipefail) + `.gitignore` covers build*.log, mod_dbg/,
  lib_dbg/, exe_dbg/.
- `e3b2b84` -- `d_lte` builds `OBJ_FILES` only.  The debug build
  exposed a real missing-symbol bug in the LTE utility:
  `lte_hydro/lte_sub.f` calls `WR_ABUND_INFO_V2`, but only
  `WR_ABUND_INFO_V3` exists (`subs/wr_x_info.f`) -- an upstream
  V2->V3 rename that lte_sub.f never followed.  Under -O3 the
  optimizer removes the unreachable call and `main_lte.exe` links
  anyway; under -O0 the link fails.  The cmfgen link needs only the
  lte_hydro archive objects, so the utility exes are out of the
  `cmfgen` target (build them manually with `(cd lte_hydro && make)`
  under nodebug if ever needed).

Verified: `./build.sh` builds end to end with zero errors
(`exe/cmfgen_dev.exe` 13:25, contains the SIGMA fix) and
`./build_debug.sh` now completes too (`exe_dbg/cmfgen_dev.exe`
13:22, -O0 -check bounds -traceback).  Both binaries postdate all
source changes.

Reviewed, no change needed in this tree (items were V6/cmfgen_source
artifacts only):

- `set_abund_clump.f` hydro-T clamp: the 0.3..1.0 (3e3..1e4 K)
  clamp exists only in the cmfgen_source variant.  V7's file is
  pristine upstream -- T_HYDRO is read from input_hydro.dat but
  unused; initial T comes from the SN_HYDRO_DATA path
  (RD_SN_DATA/SET_NEW_MODEL_ESTIMATES).  No V7 edit.
- `generate_full_matrix_v3.f` zero-population guard: V7 keeps the
  upstream unconditional `C_MAT(:,J)*POPS(J,L)` column scaling.
  The `ABS(POPS).GT.0` guard that the review flagged as fragile
  exists only in the cmfgen_source variant.  Upstream semantics
  retained per minimal-diff policy.

Runtime validation: a controlled model run on the rebuilt binary is
still required before calling the tree usable for production.
`run_model_grid.py configs/TypeIIn_10d_DEV7_STD.yaml` resolves
install_dir to this tree and dry-runs clean (job graph: one CMFGEN
job).  An interactive submission was attempted but blocked by
AssocMaxJobsLimit and cancelled; resubmit when a slot frees.  A true
VEL_LAW=12 test additionally needs a config carrying `vel_law: 12`
with existing hydro input -- none currently sets both that and
skip_stella.

### Stall-point audit fixes (29-Sep-2026, commits a9ab6b7..5e7d5cd + ed48071)

A static review of the DO_NONMONO_KNOP stall paths flagged six
hazards; all were verified against this tree and fixed as follows.
Runtime validation on a real model is still pending -- these are
compile-verified numerical guards, not a demonstrated converged run.

- `a9ab6b7` -- exit-side hysteresis on the compressing-depth set
  (fg_j_cmf_nonmono_v1.f).  A depth near SIGMA_TRUE=-1 or V=0 flapped
  the set on every NG/grid-revision pass; each flap rebuilt the ray
  geometry and reallocated I_P/I_M_DEPTH_NEXT ZEROED, cold-starting
  the GS coupling and holding GS_FREEZE forever.  Entry test
  unchanged (sigma < -1 or V < 0); exit now requires sigma > -0.98
  or V > +1 km/s.  The stored mask still uses plain thresholds so
  the FGINEXT_GS fingerprint is unaffected.
- `be1213e` -- in-loop bound checks on all five ray-bracket DO WHILE
  searches (J_PNT, both H_PNT loops; J_INDX_NM, H_INDX_NM).  A depth
  that brackets no interval previously marched K past the array end
  -- a silent hang in unchecked builds.  Now STOPs with depth/ray/
  endpoint diagnostics.
- `72641f2` -- quarantine persistence counter (mom_j_cmf_nonmono_v1.f).
  XM_BAD_NM is re-evaluated per call so a structurally failing cell
  retries forever holding a stale iterate; MOM_QUAR_AGE counts
  consecutive quarantined sweeps per cell, MOM_GS_NPQUAR reports
  cells at age 3+, printed via GET_MOM_GS_MATURITY and
  DUMP_MOM_JNEXT_GS.
- `d57d105` -- BA_NREJ_LIM (VADAT, default 20, <=0 disables).  The
  trust radius floors at 1.0, so consecutive rejections could
  alternate LAMBDA/rejected-BA forever and finish un-converged.
  Now a loud STOP.  Also committed `5c198c3` (BA_SE_RES2 rank-based
  second residual channel) which was in the tree uncommitted.
- `5e7d5cd` + `ed48071` (main repo submit_cmfgen.sh) -- GS_STATE
  file written every SOLVE_FOR_POPS call under DO_NONMONO_KNOP;
  generated slurms refuse to requeue or report completion when the
  run ends a segment >= KNOP_STALL_MAX (default 60) consecutive
  frozen iterations; .CMFGEN_KNOP_STALLED marks the run dir.
  Optional VADAT key GS_FREEZE_LIM (default 0 = off) for an
  in-model stop.

Reviewed, no change needed:

- MOM_KNOP_REFRESH_V1 is bounded (KN_M_MAXPASS=60, divergence abort,
  same-direction stall exit, bulk-drift exit) and conditionally
  retains the best-pass field -- the "burns full budget then reverts"
  failure mode is already handled; the compound stall is covered by
  the watchdog above.
- FG_KNOP_REFRESH_V1 is per-ray bounded (KN_RAY_MAXPASS) via
  KNOP_RAY_SOLVE_V1.
- Checkpoint fingerprint (bitwise R/V/SIGMA + COMP_DEPTH_IDX match)
  is deliberately strict; the mask-flap trigger is removed by
  a9ab6b7, and remaining invalidations on genuine grid changes are
  correct behavior.
- Quarantine is not a permanently latched cell state -- XM_BAD_NM is
  rebuilt every call; the persisted-store risk is now visible via
  NPQUAR rather than silently degrading the solution.

### Marginal-release BA trial (30-Sep-2026, commit 7b07fdf)

The single-tier MOM maturity gate kept the FreeNe M15 run frozen
>150 iterations on ~0.12% drifting cells that asymptote but never
converge under the one-bin upwind.  The gate is now two-tier:
MOM_UNSAFE_B is the catastrophic backstop only (capped writes, >1%
of cells drifting >50%, a >10x cell drift, >10% negJ/quarantined);
the strict maturity terms moved to MOM_FREEZE_B.  VADAT key
MOM_TRIAL_MIN (default 30, <=0 off): after that many consecutive
frozen iterations, with DO_BA_RES_ACC armed and the FG side clean,
the iteration is released to a resacc-gated BA step that rolls back
on residual growth.  N_FROZEN keeps counting, so GS_FREEZE_LIM and
KNOP_STALL_MAX still bound a run whose field never matures.
Runtime-validated: NO -- watch for 'MOM field marginal' on the first
released iteration.
(Superseded 30-Sep-2026 PM: accepted trials now refund N_FROZEN to
MOM_TRIAL_MIN -- see "Accepted-trial refund" below.)

### Numerical-robustness audit patch set (30-Sep-2026, commits 5eecc10, 978cd0c, a5d1c20, 06ec8ff, f87e1bf)

External audit of the CMFGEN/CMF_FLUX plumbing around non-monotonic
velocity fields.  All sites verified against the source before
editing; the LAW_TAU1 loop bug is upstream, identical in
../cmfgen_reference/cur_cmf.  Non-flag-gated changes are bug fixes or
degenerate-input guards; each documents why the stock path changes.

- `5eecc10` -- obs/obs_frame_sub_v9.f.  PRIVATE allocatable work
  arrays were unallocated on every non-master thread (writes through
  unallocated/shared descriptors: UB); per-thread ALLOCATE inside the
  region.  NI_MAX now uses max|V| + sum|dV| with CEILING and requires
  MAX_DEL_V_RES_ZONE > 0.  Band hunts bounded at the FREQ_CMF edges;
  sqrt radicands (ray Z, relativistic gamma) guarded; P(LS) > R(NI)
  hunt bounded at NI=1.
- `978cd0c` -- obs/cmf_flux_sub_v5.f.  LAW_TAU1 out-of-bounds read at
  I=ND fixed; TA now the cumulative trapezoidal depth, R ordering
  validated.  LAW_V1 normalized by max|V| with |V| (was V(1)).  SRCE_CHK
  target guarded against zero/subnormal ETA_CONT; ZETA/THETA divide by
  floored CHI; Rosseland accumulator and mean-opacity normalization
  guarded; 1-exp(-x) series branch for x < 1E-4; MU_AT_RMAX radicands
  guarded.
- `a5d1c20` -- plane/fg_j_cmf_nonmono_v1.f.  chi_hat first-failure
  diagnostic made reachable (test now runs on the raw value, catches
  Inf as well as NaN, then clamps).  NEW_R_GRID uses relative
  tolerances for R and VDOP_FRAC instead of exact .NE.
- `06ec8ff` -- new_main/subs/recompute_sigma_true.f +
  new_main/cmfgen_sub.f.  SIGMA_TRUE_CORE returns cleanly for ND < 3;
  clipped cells report locations.  Mean-opacity denominators
  (RLUMST, J_INT, ROSS_MEAN) guarded; 1-exp(-x) and x*coth(x/2)-2
  given stable small-x branches.
- `f87e1bf` -- plane/mom_j_cmf_nonmono_v1.f.  K_ON_J/NMID_ON_J/
  NMID_ON_HMID sanitized after GET_MOMS_NON_REL, after refresh-store
  loads, and after the solve; HUP_NM checks H_NEXT_FRQ/HNEXT_RF_READ
  entries for finite values.  Pre-existing uncommitted 29-30-Sep WIP
  in this file remains uncommitted.

Build: make cmfgen + make -C obs clean (ifort -qopenmp, nodebug);
new diagnostic strings verified in exe/cmfgen_dev.exe, exe/cmf_flux.exe,
exe/obs_frame.exe.  NOTE: `make -C obs` hits a pre-existing missing
`f77` for disp/mod_disp.f; build it with ifort first or run
`(cd disp && make)` workaround as before.

Runtime validation: NOT DONE.  No model has run through
run_model_grid.py against these binaries; the OMP fix especially needs
the 1-thread vs multi-thread cmf_flux reproducibility check.  Audit
item 16 (full +3 a_nu I proportional term) deliberately not
implemented -- proposed as a future runtime flag
(KNOP_USE_FULL_NU_PROP_TERM), requires checkpoint-version gating so GS
state is never mixed between operator definitions.

### Accepted-trial refund + GS_TRIAL_LIM (30-Sep-2026, solve_for_pops.f)

Live evidence from job 12040589 (Dessart m15 dif config): the
marginal-release machinery pulsed correctly -- it30/it32 trials
released, resacc accepted both (nrej=0), lum_ratio moved
0.16 -> 1.03 -- but N_FROZEN counted every frozen iteration
including accepted-release ones, and the strict raw-mature gate
(MOM_SW_MAXCHG <= 10, observed >= ~1.3e3, asymptotes ~1e3 on
marginal cells under one-bin upwind) can never produce the 2
consecutive mature sweeps needed to clear it.  The run stopped at
GS_FREEZE_LIM=60 while making real progress; job 12040591 (dbg
baseline) hit the same wall.

Change: `TRIAL_FIRED` is set when the MOM_TRIAL release arm fires.
After SOLVEBA, an accepted trial (BA_REJ_THIS=.FALSE.) refunds
N_FROZEN to MOM_TRIAL_MIN -- not zero: the post-release field is
re-immature by construction, so the next trial must re-earn marginal
status -- and increments N_TRIAL_ACC.  Rejected trials keep counting
(no refund).  New VADAT key `GS_TRIAL_LIM` (default 15, <=0 restores
the old count-forever behavior) vetoes further trials once
N_TRIAL_ACC reaches it, so a pulsing-forever run is still bounded.
N_TRIAL_ACC is written to `GS_STATE` (`n_trial_acc=`, after the
refund, so the wrapper sees post-refund n_frozen first) and adopted
on warm restart; a genuine 2-sweep mature release resets it.
`MOM_TRIAL terms` prints a 9th `trlim` term plus `ntr=`.
Runtime-validated on itest job 12042808 (see below).

### Paper-fidelity audit fixes (30-Sep-2026, mom_j_cmf_nonmono_v1.f)

Term-by-term check of the discretization against Knop I (A&A 496,
295), Knop II (A&A 501, 813) and Ceulemans+ 2025 (arXiv:2505.08309).
Structure verified correct (chi_hat as |GAM|/dlnv, sign-of-GAM
upwinding, equal-weight compressing coupling, N_CAP=CEILING(T4)
frequency matching, boundary clamp).  Two defects fixed:

- Dead midpoint CHI_HAT loop read GAM_MAG_NM before it was populated
  for the current call (uninitialized on first call, stale after)
  and its output was unconditionally overwritten by the depth-point
  loop.  Removed.
- Expanding-row upwind coefficient PSIPREV used the LAGGED
  GAM_FAC_PREV factor, so a sign flip of 1+sigma*K_on_J across
  adjacent frequencies injected a NEGATIVE upwind source -- Knop I
  Eq. 9 makes that coefficient positive by construction.  Now
  |GAM*GAM_FAC_PREV| at all three assignment sites (the header
  comment already claimed this behavior).

### Diagnostics formatting + frozen= format bug (30-Sep-2026)

Bordered `----` tables added around ITER_SUM, MOM/FG refresh
pass+exit blocks, GS sweep diagnostics, MOM_TRIAL terms, quarantine,
retained-field delta and requeue reports (fix_diagnostics.py).
Runtime caught a real descriptor bug the static audit missed: the
GS/iter line `4(A,A)` was one pair too many, leaving the logical
GS_FREEZE on an A slot -- printed `frozen=` + raw bytes.  Fixed to
`3(A,A)`; all other writes re-verified arg-by-arg.

### misc/ build repairs (30-Sep-2026, unrelated to model exe)

`make` in misc/ failed on the generic link rule missing `-lmod`
(mod_cmfgen symbols), `RETURN` in a main program
(rewrite_old_model_data.f) and an `AES14.4` format typo
(step_init_est_v2.f).  Fixed those three; `plt_steq.f` has further
pre-existing bit-rot left alone -- it does not affect cmfgen_dev.exe.

### Runtime validation: 3-iteration itest PASSED (30-Sep-2026)

`dessart_m15_2p82d_KNOP_MOM_dif_itest.yaml` (submission: debug ->
interactive partition, init_num_its: 3, use_existing_sn_hydro, all
Knop/GS flags as the production config).  Job 12042808 COMPLETED
exit 0 in ~16 min.  Verified live: bordered tables render, 9-term
MOM_TRIAL line with trlim + ntr=, `frozen= T` (format fix),
GS_STATE n_frozen first key + n_trial_acc persisted, ITER_HISTORY
columns, clean RVTJ exit, no fatal/NaN.  (First itest 12042506
passed the same way on the pre-format-fix binary -- its only visible
defect was the frozen= bytes.)

Full model resubmitted through run_model_grid.py:
`dessart_m15_2p82d_KNOP_MOM_dbg_dif.yaml` (+ GS_TRIAL_LIM=15):
FINAL 12043642, OBS 12043643 on small (16 CPUs, 64G, NUM_ITS=500),
cold start on cmfgen_dev.exe rebuilt with the refund, PSIPREV fix,
dead-loop removal and all table formatting.  Production CMFGEN
default raised to 16 CPUs (run_model_grid.py CMFGEN_PROD_CPUS) per
the R1 thread-scaling analysis; debug/interactive stays at 8.
### Non-monotonic mode prerequisites hard-enforced (03-Oct-2026)

The DO_NONMONO_SOB / DO_NONMONO_KNOP line-treatment gates in
rd_control_variables.f were promoted from warnings to hard STOPs, and
rd_cmf_flux_controls.f now applies the same gates to CMF_FLUX_PARAM so
the OBS spectrum cannot be run with a line treatment inconsistent with
the populations:

- DO_NONMONO_SOB: GLOBAL_LINE=SOB, or NONE with every present ion
  TRANS_*=SOB; USE_LAM_ES=T required (the F branch runs the
  clamped-sigma MOM_J_CMF_V11, outside the validated setup);
  LAM_SET=T warns only.  DO_NONLOC_SOB unchanged (warns if SOB off).
- DO_NONMONO_KNOP: GLOBAL_LINE=BLANK, or NONE with every present ion
  TRANS_*=BLANK; LAM_SET=T now STOPs (it rewrites edge lines to SOB);
  USE_FIXED_J / USE_FORMAL_REL / USE_J_REL / USE_DJDT_RTE /
  PLANE_PARALLEL(_NO_V) STOP; INC_GRID/ACCURATE unchanged; both
  USE_LAM_ES branches remain valid (FG formal vs MOM_J_CMF_NONMONO_V1);
  DO_BA_RES_ACC=F warns.
- Both flags T STOP (pre-existing); both F leaves the upstream path.
- GLOBAL_LINE=NONE still requires a TRANS_* key for every species
  slot in VADAT - a missing key aborts in GET_KEY_STRING upstream of
  these checks (pre-existing GET_KEY_STRING behaviour).

Verified on a clone of the jfnk_smoke model dir (no run-model-grid
submission; gate-matrix exercised with the built exe): KNOP+SOB STOP,
KNOP+LAM_SET=T STOP, SOB+KNOP both T STOP, SOB+GLOBAL_LINE=BLANK STOP,
SOB+USE_LAM_ES=F STOP, SOB+SOB+LAM_ES=T proceeded past the reader,
KNOP+BLANK (reference config) proceeded, both flags F proceeded.
control_variable_mod.f flag comments updated; vault notes
[[VADAT keys]], [[Sobolev non-mono solver]], [[Knop formal solver]],
[[Model pipeline]] updated in the same change.

### JFNK order-of-operations audit fixes (03-Oct-2026, commits 9956744, e8abf97, 900d515)

An audit of the JFNK order of operations (admission -> driver ->
residual/JVP -> eval sweep -> fallback) found the architecture sound
but identified defects in state save/restore coverage and the
self-disable path:

- JFNK_DISABLED was a no-op (jfnk_solve_corr_v1.f): after
  JFNK_MAX_STAGNATE consecutive failures the disabled branch set the
  LAMBDA flags but left DO_JFNK_NOW true, so Newton ran anyway under
  mutated module flags.  Now clears DO_JFNK_NOW; dead FORCE_LAMBDA
  local removed.
- Shared helpers JFNK_SAVE_MOD_STEQ_V1 / JFNK_RESTORE_MOD_STEQ_V1
  (jfnk_cmfgen_interface_v1.f) now bracket every
  JFNK_COMPUTE_STEQ_V1 / eval-sweep site and cover SE%STEQ_ADV and
  SE%T_EHB in addition to SE%STEQ/STEQ_ED/STEQ_T/STEQ_T_EHB.  The
  advective and EHB arrays previously leaked perturbed values into
  GENERATE_FULL_MATRIX_V3 (SOLVEBA right-hand side) and the resacc
  metric.
- Frozen mode keeps the module STEQ at the full outer-iteration
  value for the whole activation (the radiative split is now
  bracketed, and the frozen baseline rebuild is bracketed too), so
  a same-activation SOLVEBA_V13 fallback solves the residual the
  outer iteration actually measured.
- Coupled eval admitted without a registered JFNK_SWEEP_CB now
  returns the budget-exhaustion sentinel instead of silently
  evaluating a collisional-only residual.
- JFNK_RAD_SWEEP_V1 guards the ACCESS_F park/restore on
  JFNK_EDD_START_REC > 0 (the capture block is skipped when
  EDD_CONT=F and ACCURATE=F, leaving the -1 sentinel).
- Stale FREEZE_I_NEXT comment in fg_j_cmf_nonmono_v1.f corrected.

Verified: ./build.sh clean (zero errors, only the pre-existing
executable-stack ld warning on cmfgen_sub.o); new helper symbols
and the no-callback sentinel string present in exe/cmfgen_dev.exe;
all changes inside JFNK-gated code paths so the default
(USE_JFNK=F) path is untouched.  Runtime validation through
run_model_grid.py still pending - do not treat as a working
feature until a coupled JFNK activation is exercised on a model.
Vault note [[JFNK solver]] updated in the same change.

Runs resubmitted 03-Oct on the fixed binary (`7221de9`; the stale-
binary jobs were scancelled and the old validation run dirs
deleted): SOB_mask 12073333 (+OBS 12073334,
cmfgen_validation_SOB/m15_2p82d_mask), KNOP_mask 12073335 (+OBS
12073336, cmfgen_validation/m15_2p82d_knop_mask), JFNK_opt
12073580 (+OBS 12073581), JFNK_hyb 12073582 (+OBS 12073583) under
cmfgen_validation_JFNK/.  SOB_v2 (12072559) kept running on the
pre-fix binary - its config exercises no JFNK-gated code.

### JFNK/NG follow-up robustness fixes (05-Oct-2026)

- Newton retry state is restored to the base point before another
  trust-region solve; the GMRES recompute pass is explicit; the
  stagnation test compares preconditioned norms; trust-region norm
  bookkeeping no longer combines coefficients from restarted bases.
- The PTC preconditioner now rescales the actual BA solve by the signed
  effective diagonal ratio. Null rows use the PTC diagonal; the cached
  effective preconditioner diagonal includes LU ridge and identity-row
  repairs (including the accepted banded matrix). The banded solve now
  shares log-x/PTC/mask post-processing. FD retry floors compare the
  perturbation norm h*||v||inf, not scalar h.
- JFNK admission requires USE_LAM_ES=T because the MOM J_NEXT_FRQ /
  H_NEXT_FRQ stores do not yet have an eval freeze. NG history resets on
  held populations; NG_DEFERRED preserves a +1 NG target across an
  auxiliary SCRTEMP counter increment.
- Coupled sweep evals cold-start then restore FOLD, restore EDDINGTON,
  ACCESS_F, and two-photon accumulators, stop the frequency loop after
  JFNK_EVAL_BROKEN, and convert the former eval-time STOP sites into a
  rejected sentinel residual. Outside JFNK_IN_EVAL, the original STOP
  behavior remains.

Verified 05-Oct: `./build.sh` completed with zero errors; all edited
Fortran units recompiled. `exe/cmfgen_dev.exe` contains the broken-eval
fallback and no-callback diagnostics. The build emits the existing
executable-stack linker warning and Makefile duplicate-recipe warning.
The convergence diagnostics tests passed 16/16 (OUTGEN parsing skipped
because no sample OUTGEN was present); targeted `check_convergence.py`
masked/zero/genuine-convergence fixtures passed 3/3. Changes to the
JFNK solver/interface and eval paths are gated; with USE_JFNK disabled
those solver/eval branches do not run. NG history bookkeeping only
changes NG decisions when DO_NG is enabled.

Runtime validation remains NOT DONE: `squeue -u sbrennan` was empty
and the JFNK validation run directories contained no OUTGEN or
ITER_HISTORY. Do not claim coupled JFNK production-ready until a run
through `scripts/run_model_grid.py` shows activation or deterministic
deferral/fallback.
### FG ray-grid V/SIGMA_TRUE fingerprint (05-Oct-2026, shock-robustness review D1)

`FG_J_CMF_NONMONO_V1`'s `NEW_R_GRID` test previously fired only on R
drift, VDOP_FRAC/SOLUTION_OPTIONS change, or compression-set
membership change (hysteresis).  A velocity update at fixed R -- the
DO_V_PARAMS wind revision is the known path; anything that perturbs V
or sigma_true without crossing the set boundaries -- left
V_COEF/SIGMA_COEF, GAM, chi_hat and the ray insertion spacing on the
old field.  `MOM_J_CMF_NONMONO_V1` already fingerprints V (`V_SM_SV`);
the FG solver now does the same, comparing V and SIGMA_TRUE_IN
against the stored main-grid segment of V_EXT/SIGMA_EXT at 1E-6
relative (floor 1.0 in natural units).  No new arrays; the extended
grid is itself the snapshot.

Verified 05-Oct: `./build.sh` clean (zero errors);
fg_j_cmf_nonmono_v1.f recompiled; `V drift`/`SIGMA_TRUE drift`
diagnostic strings present in exe/cmfgen_dev.exe.  The routine only
runs under DO_NONMONO_KNOP, so the default path is unchanged.
Runtime validation pending (run_model_grid.py).

Vault notes updated: [[Knop formal solver]] (rebuild triggers),
[[Failure modes and guards]] (guard table row).
### Shock-robustness review fixes batch 1 (05-Oct-2026)

Audit of non-monotonic velocity handling (shock fronts, velocity
inversions) produced five actionable findings; three landed here:

- D1 `266e06f`: NEW_R_GRID fingerprints V and SIGMA_TRUE (vs the
  V_EXT/SIGMA_EXT snapshot, 1E-6 relative) so a V update at fixed R
  rebuilds the ray geometry -- mirrors V_SM_SV in the MOM solver.
- D3/D4 `04157cc`: GRID_METH=VEL_SPACE rejected inside
  ADJUST_R_GRID_V4 under DO_NONMONO_SOB/DO_NONMONO_KNOP (DO_VEL_
  REGRID_V2 interpolates R(V), assumes monotone V);
  GAMRAY_TRANS=RAD_TRANS rejected under DO_NONMONO_KNOP
  (GAMRAY_SUB_V3/CMF_FORMAL_REL_V4_GAM is a signed-sigma solver fed
  the unclamped gradient).
- THK extension in FG_J_CMF_NONMONO_V1 floors the beta-law exponent
  source at sigma >= -0.999: SIGMA is not re-clamped under KNOP, so
  sigma(1) < -1 made BETA < 0 and V_EXT diverged.  Warns when the
  floor engages (THK_CONT=F is the pipeline default; latent defect).

Verified 05-Oct: `make cmfgen` clean (zero errors), edited files
recompiled, new diagnostic strings in exe/cmfgen_dev.exe.  Gates are
flag-gated; default path unchanged.  Runtime validation pending.
### Shock-robustness review: DO_NONLOC_SOB hardening (05-Oct-2026)

D5 `ff58a30`: sob_rh_couple.f now drops RH78 kernel pairs with
V(I)<=0 or V(J)<=0 (the mu^2 geometry uses V^2 and hides sign flips;
the pair weight carried signed V(J)) and reports the drop count.
GAMH_L uses |V| floored at 1 km/s -- the 1/V table could hold Inf at
a stagnation depth and leak into the net rates.  Infall regions now
get no non-local correction (documented limitation of the kernel
under DO_NONLOC_SOB).

Also note: a concurrent uncommitted solve_for_pops.f edit dropped a
`//` between two WRITE strings (line ~1035); restored so the build
passes.  That WIP (BA_LUM_GATE damping) is left uncommitted for its
owner.

### Shock-safe luminosity diagnostics (06-Oct-2026)

Inner/outer luminosity agreement and drift are no longer convergence
conditions: shocks permit physical energy deposition between the two
boundaries.  `BA_LUM_GATE` / `BA_LUM_DCHG` remain readable no-op keys
for restart compatibility.  `ITER_SUM` reports `lum diag` ratio and
relative drift with `[INFO]`; no luminosity quantity demotes, holds,
caps, or delays SOB or Knop/JFNK updates.  Build and runtime validation
are recorded with the associated submission.

### SOB-to-JFNK staged workflow (06-Oct-2026)

The production sequence keeps the hydro snapshot temperature fixed and
uses a converged SOB population checkpoint to initialize Knop/frozen
JFNK.  The target config sets `warm_restart_from` and
`warm_restart_population_only: true`.  The workflow requires identical
source/target `SN_HYDRO_DATA` and matching MODEL NT/ND, then stages only
POINT1/POINT2/SCRTEMP/MODEL.  It does not copy RVTJ, EDDFACTOR,
FGINEXT_GS, or MOMJNXT_GS across transfer operators.  The target rebuilds
its Knop radiation and lagged-frequency state while preserving the SOB
population basin.  SOB and JFNK configs use the same diffusion inner
boundary, thin outer boundary, dissolution/penetration physics, and
fixed-T controls.

### FIX_IMPURITY disabled in pipeline profiles (06-Oct-2026)

`FIX_IMP=T` was implicated in the SOB stall at job `12101701`: the
generated VADAT confirmed the switch, while CORRECTION_SUM reported
mostly zero corrections during full iterations.  The mechanism is
FIXPOP_IN_BA_V3 zeroing BA equations for low-abundance/low-fraction
stages at each depth; that suppresses exactly the population corrections
the full iteration must apply.  All production SOB/KNOP/JFNK YAML
profiles and the `NON_MONO_SOB_RECOMMENDED` recipe now stamp
`FIX_IMP=F`.  The VADAT key and FIXPOP_IN_BA_V3 remain for explicit
experiments, but are no longer part of the submitted solver recipe.
Fresh SOB submission: CMFGEN `12103504`, OBS `12103505`, convergence
monitor `12103506` on source commit `c13f017`; generated VADAT confirms
`F [FIX_IMP]`.

### Correction-cap defaults for the non-monotonic recipes (06-Oct-2026)

Live SOB job `12103504` is bounded, not divergent: the it21
`BA_PROBE_LAM` trial increased the SE residual and was rolled back;
subsequent Lambda calls reduced `SE_RES` from `6.33e13` to `4.06e13`
by it38.  `MAX_LAM=1e10` is already effectively unclamped and was not
limiting the run.  `MAX_LIN` is now `1e1` in all current SOB, Knop and
JFNK YAML profiles and in the non-monotonic/default recipes; saturated
or probe steps remain separately bounded by `BA_DAMP_CAP=1.3` and
residual acceptance.  Do not tighten `MAX_LAM` as a generic fix: small
Lambda caps made the LTE-to-NLTE cascade impractically slow.  Tune
`BA_DAMP_CAP`, `BA_ENTER_RAW`, or the residual gates only when a live
full-step failure shows that the ordinary `MAX_LIN` cap is binding.

### Fresh frozen-J JFNK submission (06-Oct-2026)

Submitted the frozen-J JFNK validation run
(`DO_NONMONO_KNOP=T`, `USE_JFNK=T`, `JFNK_COUPLED=F`) as a cold start
through `run_model_grid.py` on the standing resolved config
`cmfgen_validation_JFNK/m15_2p82d_jfnk_cpl/config_resolved.yaml`
(model `dessart_m15_2p82d_JFNK_frz_v1`).  The SOB-seeded
`warm_restart_from` keys are absent, so no checkpoint state is staged.
The emptied epoch directory was repopulated with the approved
`SN_HYDRO_DATA` (md5 `ea6c5d5bd2d1c0ae785d70d4c9a56053`, identical to
the SOB/KNOP validation copies) so `use_existing_sn_hydro` applies and
the inner-boundary fingerprint matches the aborted warm-restart design.
Jobs: CMFGEN `12104645` (apu, running), OBS `12104646` (afterok),
convergence monitor `12104647` (30 s cadence).  Source commit
`782c04c`; `cmfgen_dev.exe` (14:52) includes `9213576`.  Generated
VADAT echoes `T [DO_NONMONO_KNOP]`, `T [USE_JFNK]`, `F [JFNK_COUPLED]`.
NOTE: this resolved config predates the `FIX_IMP=F`/`MAX_LIN=1e1`
profile updates and still stamps `FIX_IMP=T`, `MAX_LIN=3.0D0` --
deliberately retained to reproduce the prior recipe; revisit if the
full-iteration corrections stall as in the SOB run.

### VEL_LAW 13/14 + JFNK consistency pass (08-Oct-2026, branch devin/m15-zlaw14-jfnk)

Work on the stalled `m15_2p82d_SOB` 2.82d convergence, all submitted
through `run_model_grid.py`. Commits `c958555`..`6296181`.

- `new_main/subs/rv_sn_model_snhyd.f` (new): `VEL_LAW=13` reads
  `SN_HYDRO_DATA` with the law-12 tau regridding; `VEL_LAW=14` builds
  a blended `dlog(tau_Ross)+1.5*dlog(R)` grid from the kappa column
  (thin layers get radial coverage, thick regions keep tau coverage,
  no 2-dex seam at the zone handoff). Dispatched in `cmfgen_sub.f`,
  flagged SN-model in `rd_control_variables.f`, built via
  `new_main/subs/Makefile`.
- `do_tau_regrid_v2.f`: TAU_OLD(1) uses floored TA, not FLUX_MEAN
  (zero outer opacity put -Inf into the tau grid).
- `cmfgen_sub.f`: DO_LAM_AUTO stable-count guard reverted -- it
  parked non-mono runs in LAMBDA at ~45% demand.
- JFNK consistency: `CMAT_RAD` cached across activations
  (`JFNK_CMAT_MAX_AGE=4`; fixed point is invariant to BA_rad);
  deterministic re-arm `JFNK_REARM_EVERY=10` admitted calls replaces
  the process-lifetime disable that silently reset on requeue;
  SOB zlaw14 config raised `JFNK_MAX_EVALS` 40 -> 200 (one GMRES
  cycle alone can spend ~60 JVPs).
- `solve_for_pops.f`: per-great-iteration `ITER_TYPE` banner naming
  the dispatched solve (LAMBDA/FULL/JFNK-FRZ/JFNK-CPL).
- `Makefile_definitions`: viper `INSTALL_DIR` is self-locating
  (this file's own directory) because ptmp1 and ptmp2 are two
  worktrees of the same repo and each must install into itself.
- `tests/f90` (30 checks: solver core vs stub int_mod, real
  interface routines) + `tests/python/jfnk_outgen_audit.py` (live
  OUTGEN validator + plots under `tests/plots/`). All tests pass.
- `build.sh` now runs the tests/f90 suite after linking; a failed
  check fails the build like a compile error. The suite also writes
  per-iterate |F|/dtau histories that
  `tests/python/plot_solver_tests.py` renders to
  `tests/plots/solver/newton_convergence.png` (plot step non-fatal).

Submitted: FINAL `12131626`, OBS `12131627`, monitor `12131628`
(config `dessart_m15_2p82d_SOB_jfnk_zlaw14_jfnk.yaml`,
`install_dir` ptmp2 worktree). VADAT verified `DO_JFNK=T`,
`JFNK_MAX_EVALS=200`, `VEL_LAW=14`; blended grid executed.
Runtime JFNK admission under SOB: UNVERIFIED -- still deferred
(`LAMBDA-designated iteration`) at it15, normal for the mode
controller while MAXCH is high.

### Status update (2026-10-08, later): null-row projection fix

- First live JFNK admissions (it31/33/35, FINAL 12131626) all
  aborted at the first GMRES matvec: `precond clip -- aborting`,
  info=3, zero useful evals. Root cause: identity-replaced near-null
  BA rows returned z_i=v_i with residual components ~1e12, so the
  |z|>1e8 nonlinear clip fired on every apply. Three consecutive
  aborts disabled JFNK for the rest of that process.
- Fix: `JFNK_ROW_REG(J,K)` mask persisted at factor time; masked
  components projected to zero at apply in all branches (LU, LU-fail
  Jacobi, DGETRS-fail, pure Jacobi, banded, log/PTC post-loop).
  Linear projection keeps M^-1 consistent; no clip event. Detection
  reads raw BA_A (ridge retries cannot un-mark). Counter
  `jfnk_n_reg_proj` in diagnostics.
- Commit 0aa5c4c on main + devin/m15-zlaw14-jfnk, pushed. Fresh
  cmfgen_dev.exe 17:26; 33/33 build-time tests pass.
- Job 12131626 graceful-stopped at it37, now PENDING
  (AssocMaxJobsLimit); requeue picks up the fixed binary and the
  JFNK disable resets on process start.
- Also fixed: .mod-shadowing hazard in tests/f90 (cwd .mod resolves
  ahead of -I; build_and_run.sh now cleans them).

### Status update (2026-10-09): cross-activation state disarm

Numerical audit of module state that survives between
`JFNK_SOLVE_CORR_V1` activations found two leaks; both fixed and
covered by the build-time suite (33/33 checks pass).

- Stale PTC scale/anchor (`jfnk_ptc_scale_ready`,
  `jfnk_ptc_anchor_ready`, `jfnk_ptc_x0`) persisted into the next
  activation's F(0) identity test, which runs before
  `JFNK_SET_PTC_SCALE_V1` re-arms.  Frozen mode would fail the test
  on every reactivation (permanent JFNK disable); coupled mode
  adopted a contaminated baseline and polluted every FD JVP.
  `jfnk_cmfgen_interface_v1.f` now clears both flags at INIT and
  clears `scale_ready` before the rebuild-baseline eval.
- `gmres_m_v1` kept its hookstep basis descriptor (`jj`, `beta`) in
  SAVE locals.  `JFNK_DIAGNOSTIC=T` runs the m-dimension experiment
  (m_test up to 40 > production m), leaving a `jj` that a later
  info==2 recompute indexes past the caller's m-sized h/v/y arrays.
  The pair is now module state (`jfnk_gmres_jj`,
  `jfnk_gmres_beta`), checkpointed per Arnoldi cycle, guarded by
  `jj <= m` on recompute, and saved/restored around the diagnostic
  experiment.
- Detail: `cmfgen_obsidian/solvers/JFNK solver.md`,
  "Cross-activation state disarm (09-Oct-2026)".

### Status update (2026-10-09): coupled JFNK under SOB + cooldown/re-probe

The SOB coupled tier is implemented and runtime-exercised.
`JFNK_SOB_LINE_PASS_V1` (cmfgen_sub.f 6241-6483) replays the
`SUB_SOB_LINE_V3` line loop at the perturbed populations at the end
of `JFNK_RAD_SWEEP_V1`; `JFNK_COUPLED=T` is no longer rejected
under `DO_NONMONO_SOB` (commit `0f6b9f6`).  First runtime test on
`m15_2p82d_SOB_jfnk_cpl` (jobs 12143251/12143252): eval machinery
verified (F_eval vs F_outer ~4-8% rel diff, deterministic evals,
clean state restore), but GMRES stalled on all three admissions --
the residual sits in near-null field-consistency directions of the
coupled Jacobian that the local BA preconditioner cannot reach.
No accepted coupled Newton step yet; admission gates and Lambda
contraction remain the working convergence path.

Lockout is now a cooldown, not permanent: `JFNK_MAX_STAGNATE`
(VADAT, default 3, clamp 1-20) consecutive failed solves park JFNK;
every `JFNK_REARM_EVERY` eligible calls (VADAT, default 10, clamp
1-500) a single-strike re-probe runs; a failed probe re-locks, a
successful one resets the strike count (commits `7e519e6`,
`8d7e068`).  Cooled calls print `JFNK DEFER` at admission.

JFNK OUTGEN prints are now one value per line with expanded labels
(`res` -> `residual`, `it` -> `iteration`, `dtau` ->
`pseudo-timestep`, etc.; commits `9131049`, `ffc18ce`).  Binary
relinked 16:51; 34/34 JFNK tests pass.

### Status update (2026-10-10): SCHUR preconditioner option + git re-init

Research-plan integration: `JFNK_PRECOND=SCHUR` adds a two-level
(n,T) Schur-complement preconditioner.  Routines
`jfnk_schur_form_v1` / `jfnk_schur_apply_v1` live in
`jfnk_solver_v1.f90` so the unit-test driver (`se_driver`) calls the
identical production code via `external` declarations - satisfies the
"exact Fortran in unit tests" requirement.  Per-depth pop blocks use
normal-equation LSQ (Tikhonov `lambda^2 I`, ridge rows, reg-row
normalization) WITHOUT row equilibration: measured bisection on the
shock surrogate shows the `rsf^2` reweighting corrupts `J_nn^-1`
(MZDIFF 0.25 vs 0.0075 raw).  Dense NDxND Schur complement via
DGETRF/DGETRS.

Surrogate verdict (se_driver, systems/shock): SCHUR INFO 3 vs LSQ
INFO 0 NITS 25 - the pop-only elimination drops the dominant T-row
regularization on masked-T specs.  LSQ remains the production
default; SCHUR is an opt-in VADAT value pending a cross-depth-T
surrogate (`tcoup` knob planned for gen_se_system.py).

Git: the worktree parent `cmfgen_nonmono_dev_V7_OLD/.git` is gone
(tree deleted; only a dangling `.git` pointer file remained).  The
tree was re-initialized in place; baseline commit `a15dba1` covers
all session solver work.  Prior history is lost with the parent.
Fresh run submitted against this state: job 12151424
(`m15_2p82d_SOB_jfnk_zlaw14_fresh`, cold start, JFNK_TAUTO=T,
JFNK_PRECOND=LSQ).  Binary relinked 23:07; unit tests pass except
the documented LU-on-shock rank-deficiency blowup (pre-existing).

### Status update (2026-10-12): DO_SOB_FLUX_OBS switched to profile injection

The EW deposition render (`DEP_SOB_OBS_FLUX`, commit `7cbc53b`) was
replaced before validation: it invented each line's profile shape
(flat/one-sided box over the resonance span), and job 12165745
(`m15_2p82d_SOB_jfnk_alo_probe`) showed the expected boxy emission
bumps and comb-like UV artifacts.

`ADD_SOB_LINE_PROF` (obs/cmf_flux_sub_v5.f, commit `eb3408a`) instead
writes each `TRANS_TYPE='SOB'` blend member's profile-integrated
opacity/emissivity (`CHIL_MAT`/`ETAL_MAT` x `CLUMP_FAC`) into the
stored `ETA_CMF_ST`/`CHI_CMF_ST` as a normalized Doppler core at the
line's own frequency, then re-runs `OBS_FRAME_SUB_V9` after the EW
loop so the observer's-frame formal integral draws resolved profiles
from the isovelocity geometry (Dessart et al. 2015, Sec. 3.1).  Core
weights normalized by local cell width preserve the integrated line
opacity on a coarse CMF grid; masering cells (`CHIL_MAT<0`)
contribute zero line opacity.  `OBSFRAME` is rewritten with the SOB
render; `OBSFLUX` keeps continuum+BLANK as a comparison product.
`ETA_DATA`/`CHI_DATA` are written before the EW block, so a
standalone `obs_frame.exe` render still lacks SOB lines.

Verified so far: incremental `make -C obs` clean, `add_sob_line_prof_`
and the `OBS_FRAME_SOB` TUNE label present in `exe/cmf_flux.exe`,
`dep_sob_obs_flux_` gone.  Runtime verification pending: OBS-only
resubmit of `dessart_m15_2p82d_SOB_jfnk_alo_probe.yaml`
(`skip_final: true`, `obs_line_treatment: sob`).

### Status update (2026-10-12): non-mono lagged-map fixes (a)+(c) probed

Two flag-gated stabilizers for `MOM_J_CMF_NONMONO_V1`, targeting the
observed LAM_ES=F gain channels (see `cmfgen_obsidian/audits/
MOM_KNOP_STIFFNESS_AUDIT.md` and the Sobolev non-mono solver note):

- `MOM_GAM_PIN` (commit `66756a3`): pins `GAM_FAC_PREV`/`GAMH_FAC_PREV`
  /`EPS_PREV`/`PSIPREV(1)` to their current-frequency evaluations so
  the explicit lagged weights can no longer exceed the implicit ones
  near the `1 + sigma*K = 0` resonance flip.
- `MOM_LAG_BND` (commit `6373327`): bounds the total lagged
  frequency-source weight per row by the assembled diagonal excess
  (per-row `LAG_SCL_NM`, applied to the interior assembly, both outer-
  boundary branches, and the negative-J rebuild path).  Sweeps report
  `lagbnd scaled=<n> mins=<min scale>` in OUTGEN.

Both default TRUE inside the non-mono solver only (same precedent as
`MOM_EPS_LIM`); the stock monotonic path is untouched.  VADAT keys
added in `rd_control_variables.f`.  Also committed: `b533eee`
(test_vadat_keys fix for module compile order + stale LAM_ES abort
expectation) and `bdf2f10` (docs).  A pre-existing user OBS maser
fix sits at `2db3338`, unrelated.

Build: `./build.sh` clean; F90 harness T1-T3 pass; python suite
10/11 files pass (`test_interp_guards.py` fails 3 checks, verified
pre-existing against the unmodified reference routines).

Probes (all frozen `exe_frozen` = sha 98cc0b41fccec88c, git bdf2f10):

| Job | Config | State |
|-----|--------|-------|
| 12167354 (+OBS 12167355) | SOB cold (`mom_cold`) | RUNNING - converging (pops ~50, L~4.7E6) |
| 12167426 (+OBS 12167427) | KNOP cold (`KNOP_MOM_cold`, resubmit of diverged 12165640 on stale 32b0d2a binary) | CANCELLED it19 - field never matured (mom floored ~3-8E3 at lagbnd-nulled pure-scattering rows) |
| 12167431 (+OBS 12167432) | SOB warm probe | RUNNING - converging |
| 12167439 (+OBS 12167440) | KNOP warm probe | RUNNING - freeze/release cycling productively (se_res 7.9E16 -> 1.2E16 on first live BA step) |
| 12168927 (+OBS 12168928) | KNOP cold with MOM_ABS_FLR=1D-3 (fix d) | submitted on binary 964638e (sha c50238ef) |

Old KNOP cold 12165640 (+dep 12165641) cancelled: its `exe_frozen`
was pinned at 19:35 to pre-fix commit `32b0d2a` and could never pick
up the fixes across requeues; it was already at `se_res=5E44`,
`mom_maxchg=6E6`, `n_frozen=14` at it14.

Early signal from SOB cold it1-6: luminosity bounded and slowly
falling (8.3E6 -> 7.6E6; the unfixed runs grew 1E7 -> 1E21+), zero
negJ/quarantine repairs, refresh bulk-converged by it3.  `lagbnd`
clips ~1.4M entries/sweep with min scale 0 (some rows have zero
diagonal excess, so their whole lagged source is nulled - watch
whether clipping contracts as the field settles).  NOT converged -
too early to classify; the prior SOB cold divergence was slow.

## 14-Oct-2026 (Devin session, cont.) - MOM_TRIAL_LAM Lambda-probe trial
KNOP warm (12167439) diagnosed in a freeze->release->reject->refreeze
limit cycle: every resacc-BA trial at demand ~1E7 applied a capped
400% step, next-call SE residual ballooned ~37x (9.4E13 -> 3.5E15),
resacc rolled it back.  LAMBDA-scheduled releases taking the damped
GS_LAM_CAP step worked instead (it524: residual 4.97E14 -> 8.8E13,
demand 2.2E8 -> 7.8E7 in one call).  New VADAT key MOM_TRIAL_LAM
(default 0 = off): MAXCH above the key fires the trial as a damped
GS_LAM_CAP LAMBDA step instead of the BA step.

Concurrency note: a parallel session added RLOC_SHOCK_POWER to the
same files (control_variable_mod.f, rd_control_variables.f,
drv_vadat.f90, test_vadat_keys.py); an earlier overwrite ate the
MOM_TRIAL_LAM edits once.  Re-applied on top; all 32 vadat tests pass.
Left UNCOMMITTED to avoid racing the parallel session's RLOC commit.

Probes on binary built 10:40 (pre-RLOC merge, has MOM_TRIAL_LAM):
- KNOP cold resubmit 12175438: full flag stack (H_IMPL, ABS_FLR=1e-3,
  RF_DIR=DESC, FAC_RLX=0.5) + TRIAL_LAM=1e3.  Fresh cold start.
- KNOP warm resubmit 12175441: base stack + TRIAL_LAM=1e3.  Warm
  restart lost MOMJNXT_GS (tag 532 vs NITSF 440, cancelled mid-it) ->
  cold reseed + pops from it440 checkpoint.

## 14-Oct-2026 (Devin session, cont.2) - SOB cold empty-ion fix + stall watchdog
SOB cold demand oscillation (bounded ~200% with spikes to 1.7E4 across
~400 its) traced to a single row: CSIX lev 1, row 246 in POPS, at the
outer compression depths (I~3, sigma_true~-1.04).  Population ~4.6E-26
vs depth max ~6.5E-4 (rel frac ~7E-23): above the 1E-30 absolute
PFLOOR but far below any physical threshold -- classic empty-ion
chatter, 65% of all 'largest' BA entries.  Runner-ups ~1-100% real
ejecta ions plus similar empty trace levels (CV, MgIV, CaV).
Enabled BA_DEM_RELFRAC=1.0D-6 in the cold config (drops rows below
1E-6 x depth-max from the demand metric; corrections still commit)
and resubmitted as 12176092 on binary 4f72887.

KNOP stall watchdog gap found+fixed: the n_frozen>=KNOP_STALL_MAX(60)
exit-7 block existed only in cmfgen_base/submit_cmfgen.sh, not the
cmfgen_base_dessart/cmfgen_base_fullspec templates these m15 runs
generate from -- a stalled KNOP run would have requeued forever.
Ported the watchdog block + GS_STATE cleanup lines to both templates
(bash -n clean).  Already-generated run slurms unaffected.

## 13-Oct-2026 (Devin session) - X_SOB 40d OBS render: boundary mirror + V12 floors
X_SOB OBS (DO_NONMONO_SOB populations, DO_NONMONO_KNOP=F ->
FG_J_CMF_V12 on clamped SIGMA) failed on "negative mean intensities"
at JNU_STORE(ND). Root cause: OBS template DIF=T while the population
VADAT ran DIF=F/IB_METH=ZERO_FLUX/THK_CONT=F -- the diffusion seed
I_CORE=(ETA+DBB*proj)/CHI went negative (DBB<0, optically thin inner
cells) and poisoned every core-ray sweep. Fixes: prepare_CMF_PARAM.py
mirrors DIF/IB_METH/THK_CONT from ../VADAT; fg_j_cmf_v12.f floors
I_CORE and swept I_M/I_P at 1E-30 (commits 4f72887..b35513e) with
large-value diagnostics capped at 3/ray. Job 12176536 (DIF=T + floor)
completed in ~2 min with zero sweep diagnostics and produced
obs_fin_10 with deposited Sobolev line structure. Resubmitted as
12176788 under the mirrored DIF=F/ZERO_FLUX boundary (population
consistency); compare spectra for boundary sensitivity.

## 13-Oct-2026 (Devin session, cont.) - X_SOB underluminosity + DIFFUSION repower
12176788 completed the full deposited render (clean, zero intensity
repairs) but the spectrum integrates to ~8.6e39 erg/s vs ~2e42
observed at +25.1d. Root cause is the population model, not the
render: DIF=F/IB_METH=ZERO_FLUX discarded the hydro inner-boundary
luminosity (SN_HYDRO_DATA header Lum=2.29e8 Lsun at tau_ross~99),
INC_SHOCK_POWER/INC_XRAYS/INC_RAD_DECAYS all F, and FIX_T=T -- the
ejecta only re-radiated stored heat (OUTGEN: L(d=1,ND)=0 for all
1692 its). X_SOB.yaml vadat_overrides changed to DIF=T +
IB_METH=DIFFUSION (the tau>=25 SOB auto-select the override had
blocked; m15 dif-probe and job-12001805 precedents). Resubmitted
through run_model_grid.py as FINAL 12177191 (warm_restart=T reuses
the converged starved checkpoint) + chained OBS 12177192. OBS will
mirror DIF=T/DIFFUSION; the V12 I_CORE floor protects the negative-
diffusion-seed case there.

Same session: pipeline luminosity guard (commit 4543df3d).
prepare_VADAT stamps LSTAR with precedence
cmfgen.inner_boundary_luminosity_lsun > SN_HYDRO_DATA header >
lum.xg, and check_inner_boundary_luminosity() refuses a
non-DIFFUSION boundary when LSTAR>0 and tau_inner>=thresh (25 SOB /
50 otherwise) unless cmfgen.allow_dead_inner_boundary: true.
run_model_grid.guard_inner_boundary_luminosity runs the same check
at submit when SN_HYDRO_DATA exists, so contradictory configs fail
before queueing. Sibling SN_B_X_grid yamls still carry
DIF=F/IB_METH=ZERO_FLUX and will now refuse until repowered or
acknowledged.

## 13-Oct-2026 (Devin session, cont.) - m15 2.82d MOM probes repowered to DIFFUSION

Boundary audit confirmed all four MOM probes ran
DIF=F/IB_METH=ZERO_FLUX on a luminous model (LSTAR=4.27e8 Lsun,
tau_ross=39.3) with no internal deposition (INC_SHOCK_POWER=F,
INC_RAD_DECAYS=F) -- the starved-boundary failure the new guard was
built for. Runtime evidence: RLUMST(ND)=0 exactly in every run; SOB
warm had converged to L(d=1)=3.7e6 (~115x under header), SOB cold
4.8e6; KNOP warm 2.9e8 was the FIX_T thermostat, not the true flow.
ZERO_FLUX had been inherited from the convergence-era recipe
("avoids negative DBB from non-monotonic T"), but T at the boundary
is smooth and monotonic (8.39->8.67e4 K inward) so DTDR>0 and DBB>0
-- DIFFUSION is safe, and dif_probe/dif2-4 proved it on this model
(715+ its, RLUMST(ND)=4.27e8).

Changes: all four configs now DIF=T + IB_METH=DIFFUSION; SOB pair
also JBAR_W_EDD=T (legal under DIF=T; dif_probe recipe), KNOP pair
keep JBAR_W_EDD=F (BLANK mode does not take the Jbar path). KNOP
BLANK-mode tau threshold is 50 vs SOB 25, so KNOP's tau=39.3 prints
the <50 warning but preserves the user override -- verified
"IB_METH: preserved (user override = DIFFUSION)" in OUTGEN.
warm_restart_from dropped from both probe configs so the preserve
list carries each run's own checkpoint; GS_STATE added to that
preserve list in run_model_grid.py so n_frozen/trial state survives
resubmission.

Resubmitted: SOB cold 12177450, SOB warm 12177452 (resumed it684,
own checkpoint), KNOP cold 12177454, KNOP warm 12177456 (resumed
it458, own checkpoint, GS_STATE intact). All report
RLUMST(ND)=4.265858E8 -- the boundary now carries the hydro
luminosity. No -ve-intensity stops; the KNOP nonmono DIFFUSION
branch (fg_j_cmf_nonmono_v1.f:2366, mom_j_cmf_nonmono_v1.f:2378) is
exercised for the first time at runtime. Surface L still starved-era
(SOB 3.7e6, KNOP warm 2.9e8); it must climb toward ~4.3e8 as the
repowered field equilibrates upward -- watch it as the convergence
metric alongside MAXCH.

## 13-Oct-2026 (Devin session, correction) - CMF_FLUX reset for SOB OBS

The earlier statement that OBS job 12177192 produced a completed
resolved Sobolev render is retracted. The earlier attribution of the
failure to the diffusion boundary and the claim that its floor fix made
the spectrum valid are also retracted. `obs/OUT_FLUX` ends with
`FQCOM_IBC_V2-- -ve intensities`; `SCRATCH/EWDATA` is empty, so the
Sobolev EW and profile-injection sections did not complete. The wrapper
still moved an `OBSFRAME` file to `obs_fin_10`; that file is not a valid
CMF_FLUX result for this job.

The CMF_FLUX sources were reset against
`cmfgen_reference/cur_cmf/obs`. The population-independent CMF_FLUX
path and `COMP_JCONT_V4.INC` are stock again. The only new spectral
behavior is gated by `DO_NONMONO_SOB` / `DO_SOB_FLUX_OBS`: use
`SIGMA_TRUE` for Sobolev EW estimates, set VINF from the full speed
range, and optionally inject Sobolev profiles for a second observer-frame
pass. Knop formal-J support was removed from CMF_FLUX; the config
preparer rejects Knop OBS rather than emitting an unverified spectrum.
`./build.sh` completed successfully and the changed CMF_FLUX objects
were rebuilt. `tests/test_prepare_cmf_param.py` passes (8 tests). The
broader `test_nonmono_vadat_enforcement.py` selection has five existing
failures in the already-dirty `run_model_grid.py` / `prepare_VADAT.py`
worktree, outside this CMF_FLUX change. The X_SOB `--dry-run` graph
contains only SEDONA because its YAML has `run_cmfgen: false`; no CMFGEN
job was submitted. The spectrum remains UNVERIFIED until an orchestrated
CMFGEN/OBS run passes the negative-intensity check and writes non-empty
EWDATA plus the final profile render.

## 13-Oct-2026 - KNOP held-population handoff guard (commit d1c51fa)

`new_main/cmfgen_sub.f` now defers the initial `DO_LAM_AUTO` transition
under `DO_NONMONO_KNOP` while `GS_POPS_HELD=T`. A held population step is
rolled back, so its `MAXCH_CTRL=0` must not signal convergence. The
monotonic path and SOB handoff are unchanged.

Verification: `./build.sh` passed, including 21 solver tests, 9
interface tests, and the analytical harness. Jobs 12179164 (cold) and
12179160 (resumed probe) were submitted through `run_model_grid.py`.
Cold it1 and probe it472-474 stayed `LAMBDA_GS -> LAMBDA_TFIX` with
`IN_ITS DO_LAM_IT=T`; the prior cold run had switched to full iterations
at it1 with `MAXCH_CTRL=0`. This verifies the handoff guard, not full
model convergence. The MOM path remains diagnostic-only on this field.

The frozen binaries report base HEAD 82b0e23. The source worktree was
dirty at build and submission, so that hash alone does not identify all
source included in those binaries. The new guard is recorded in commit
d1c51fa; other pre-existing worktree edits were left untouched.

## 13-Oct-2026 - KNOP refresh exit census + unmasked CORRECTION_SUM

Two diagnostic changes motivated by the m15 2.82d formal-J cold run
(interactive, job 12181576): the refresh residual plateaued at ~0.08
against KN_RAY_TOL=2e-3 with ~5.4% of significant cells still above
tolerance at it13, and CORRECTION_SUM read ~18 rows/depth demanding
corrections while the raw STEQ SOLUTION ARRAY showed ~1800/1839 rows
demanding changes.

`plane/fg_j_cmf_nonmono_v1.f`: per-sweep exit code (0 exact /
1 tol / 2 bulk / 3 same-direction stall / 4 pass cap / 5 diverged)
threaded out of SWEEP_FLD and SWEEP_FLD_ORD through
KNOP_RAY_SOLVE_V1; the refresh summary now prints a `ray exits`
count per class and, when stall/cap exits occur, their max residual
and bad-cell share. Output-only change; no solver behaviour touched.

`new_main/subs/solveba_v13.f`: CORRECTION_SUM now writes, per depth,
the all-row cumulative demand counts, then the floor-kept counts
(old masked view), then N_SKIP rows excluded by BA_DEM_PFLOOR /
BA_DEM_RELFRAC. The MAXCH control path still reads floor-passing
rows only; cmfgen_convergence.py reads columns 1-7 which now hold
the unmasked counts.

Verification: `make DEBUG=nodebug all` exit 0; `ray exits` strings
present in the rebuilt exe/cmfgen_dev.exe. Runtime check is the next
KNOP refresh block in any running job; the Aitken A/B arm
(dessart_m15_2p82d_KNOP_FORMAL_aitken) exercises KN_ORD_AITKEN=T
against the running cold arm.

### 13-Oct-2026 (cont.2) - KN_ORD_AA + KN_RAY_STALL + CORRECTION_SUM revert

`plane/fg_j_cmf_nonmono_v1.f`: three changes in the FG refresh.

1. `KN_ORD_AA` (VADAT, default F) + `KN_ORD_AA_M` (default 2):
   safeguarded type-II Anderson acceleration of the ordered pass
   iterate in SWEEP_FLD_ORD, mirroring the runtime-tested MOM
   refresh AA (weighted Gram, non-affine cell exclusion at AA_GCAP,
   Tikhonov ridge, DGESV, ||corr||>||f|| veto, 10|r| per-cell cap,
   positivity/cap fallback) plus a residual-acceptance rollback:
   if the pass after an applied mix does not improve, the pre-mix
   iterate is restored, history cleared and AA disabled for that
   solve.  History collection starts at pass 4 (lazy allocate);
   stored residuals are rescaled on OMEGA_RAY changes (exact for
   the affine map).  Per-refresh `AA: applied= vetoed= rejected=
   solve-fail=` counters print under the `ray exits` census.

2. `KN_RAY_STALL` (VADAT, default 0.90): the same-direction ratio
   threshold for the NSTALL early exit is now configurable in both
   sweeps; default preserves prior behaviour.

3. Aitken A/B verdict: KN_ORD_AITKEN=T produced 25-26 diverged
   rays/iteration vs 0 in control (residual 10.5->853 vs 1.0->0.17
   at matched iterations).  Failed arm; not for production.

`new_main/mod_subs/{control_variable_mod,rd_control_variables}.f`:
declarations + VADAT readers for KN_ORD_AA, KN_ORD_AA_M,
KN_RAY_STALL (setters SET_FG_ORD_AA, SET_FG_STALL_LIM).

`new_main/subs/solveba_v13.f`: CORRECTION_SUM reverted to the
original 8-column per-depth format but counts ALL rows (unmasked);
floor gating kept only on the DECREASE/INCREASE MAXCH scan.

Verification: ifort -syntax-only clean on all four files;
`make DEBUG=nodebug cmfgen` exit 0; KN_ORD_AA/KN_RAY_STALL/
`AA: applied` strings present in exe/cmfgen_dev.exe (relinked
20:31).  Runtime arms pending: freqmatch (DO_NONMONO_FREQ_MATCH=T,
job 12185982) plus AA/stall arms to be submitted.

## 13-Oct-2026 (evening, part 2) - KNOP stall default raised to 0.98

`plane/fg_j_cmf_nonmono_v1.f`,
`new_main/mod_subs/{control_variable_mod,rd_control_variables}.f`:
KN_RAY_STALL default 0.90 -> 0.98 in DATA init, VADAT reader
default, and all doc comments.  Evidence: the stall-relax arm cut
the unconverged fraction ~40x vs control at matched iterations
(0.2-0.3% vs 11-13%); the exit census had shown the ~12 stall
rays held nearly all bad cells while still contracting.

DO_AIT in SWEEP_FLD_ORD is now gated on .NOT. FG_ORD_AA so the
per-cell secant and the vector mix cannot compound on one iterate.

Runtime state (commit c288d94, exe relinked 21:30):
- stall arm 12186358 (KN_RAY_STALL=0.98): it4 residual 0.094,
  bad frac 0.68%, 9 stallers, cap=0 after it2 restart cleared.
- AA arm 12186580 (KN_ORD_AA=T, M=2): it4 residual 0.35,
  applied=69 vetoed=147 rejected=9; mixes largely vetoed, not
  beating stall-relax yet; opt-in stays.
- freqmatch arm 12185982: it7 residual 0.081, bad 7.5% -
  tracks control, mild benefit only.

Verification: three files compile clean under viper_intel_nodebug
flags; make DEBUG=nodebug cmfgen exit 0; 0.98 bit pattern present
in exe/cmfgen_dev.exe.  Parallel-session edits left unstaged.

## 13-Oct-2026 (evening, part 3) - KNOP field-maturity audit

Audit of `FG_KNOP_REFRESH_V1` -> `GET_FG_GS_MATURITY` ->
`SOLVE_FOR_POPS` -> `GS_STATE` / submit wrapper found that the ray
bulk census used `abs(new-old)/max(old,abs(new))`, while the residual
and tolerance use `abs(new-old)/max(old,1e-4*column_peak)`. The mismatch
can omit cells just above `KN_RAY_TOL` from both the per-ray bulk-exit
count and global `GS_NBAD`, allowing a slightly too-early maturity
pass. Both ordered and legacy sweeps now count with the same normalized
residual used by their residual metric.

Validation:
- Synthetic threshold case: for 101/100000 cells with normalized
  residual 0.002003 at tolerance 0.002, the old count was 0 although
  the exact count is 101 (> max(50,NTOT/1000)=100); the corrected
  predicate flags all 101.
- `plane/fg_j_cmf_nonmono_v1.f` compiles with production Intel flags.
- `./build.sh` exit 0, including all compile checks and JFNK T1-T3.
- `unit_tests/python/test_vadat_keys.py`: 32 passed, 0 failed.
- Existing formal-stall VADAT confirms DO_NONMONO_KNOP=T,
  KN_RAY_STALL=0.98, GS_REL_STABLE=2, GS_HOLD_POPS=T; its GS_STATE
  shows frozen=T, held=T, n_mature=0 at it8, consistent with its
  0.58% above-tolerance fraction still exceeding the 0.1% gate.

Integration caveats: FGINEXT_GS stores only max residual and cap count,
not the bulk census. This only controls skipping the process warm-up;
the new sweep's bulk gate runs before SOLVE_FOR_POPS can take a live
BA step. `GS_NCONV` counts stored compressing-depth cells but is
thresholded against `GS_NTOT` (the full significant ray census); confirm
whether that fraction is intentionally defined over the full field or
should use a separate stored-cell denominator. The wrapper also skips
its stall check when GS_STATE is absent; consider making that fail
closed for KNOP jobs.

No patched-binary model run was submitted; the run-level effect of the
corrected borderline bad-cell census remains unverified.

## 13-Oct-2026 (evening, part 4) - KNOP AA veto rescaled for marginal modes

Final stall-arm census (job 12186358, KN_RAY_STALL=0.98): bad-cell
fraction oscillated 0.21-0.85% across it1-8, never crossing the ~0.1%
maturity gate; ~8-10 stall rays held ~96% of the bad cells (worst ray
12, residual 0.076, exited at pass 16). Those rays contract at
lambda~0.98-0.99, so the correct type-II AA correction is
||f||/(1-lambda) ~ 50-100x ||f|| - the inherited MOM veto
(||corr||>||f||) rejected nearly every mix on arm 12186580
(vetoed ~3x applied) exactly where AA is designed to help.

SWEEP_FLD_ORD AA block: veto bound relaxed to non-finite or
||corr||>1E3||f|| (covers lambda<=0.999; the 10|r| per-cell cap and
the residual-rollback remain the real safeguards).  Opt-in behaviour
unchanged (KN_ORD_AA=F default); no other exits or caps touched.

Compile: production ifort flags, clean.  Runtime test pending on the
warm-restart arm off the it8 stall checkpoint.

## 14-Oct-2026 - X_SOB OBS_FRAME luminosity fix (verified)

Observer-frame spectrum (obs_fin_10) was ~170x below OBSFLUX
(7E6 vs 2.2E8 Lsun bolometric; FIR worst at ~400x). Two flag-gated
fixes in obs/, gated on DO_NONMONO_SOB.AND.DIF:

1. cmf_flux_sub_v5.f: local DTDR=(T(ND)-T(ND-1))/(R(ND-1)-R(ND)) is
   negative for the inward-decreasing base T; the DBB=1E-99 floor then
   zeroed DBC and starved the boundary. Under DO_NONMONO_SOB a negative
   DTDR is replaced by the population-solve value derived from
   STARS_LUM via the Rosseland integral (DTDR=4.36E-8 vs population
   5.36E-8).
2. obs_frame_sub_v9.f: new args IB_J_CMF/IB_H_CMF (J,H at R(ND) on the
   CMF grid, stored from RJ(ND) and HFLUX_AT_IB=DBB/(3*CHI(ND))) plus
   DO_IB_DIFF_FLUX. When set, HOLLOW_CORE is forced F (tau_ES(ND)~0.08
   was classifying the radiating interior as an empty cavity so core
   rays crossed it with zero seed), and a core ray reaching the
   boundary (SM_NRAY=NRAY) is seeded I(mu)=J+3*mu*H.

Verified on OBS job 12191101: IP_DATA pass-1 1.24E6 -> 1.99E8 Lsun;
obs_fin_10 bolometric 2.00E8 vs OBSFLUX 2.16E8 (93%); FIR ratio ~0.9,
optical ~2. The spectrum is still far too red: 91% of the luminosity
escapes longward of 1 um because FIX_T=T holds the gas at ~3000 K and
the injected luminosity free-streams through the FIR-transparent
envelope; the observed optical continuum needs a ~9e3 K photosphere.
INC_SHOCK_POWER + temperature release (IN_ITS pins DO_T_AUTO=F) is
the remaining lever; deferred pending user decision.

Temporary diagnostics used during the hunt (IB_DATA dump, RAY_DIAG
per-ray dump, IB_CNT seed counter) were removed after verification;
the diff is ~90 lines across cmf_flux_sub_v5.f, obs_frame_sub_v9.f,
obs_frame_v2.f.

## 14-Oct-2026 (part 2) - FIX_T=T enforced as pipeline policy

User directive: input SN_HYDRO_DATA temperatures come from the
hydrodynamic simulation and are authoritative; CMFGEN must not iterate
or auto-release T on any run. Enforcement points:

- prepare_VADAT.py apply_profile(): stamps FIX_T=T/FIX_T_AUTO=F after
  vadat_overrides, warns on any conflicting request.
- prepare_VADAT.py final override pass (same block that strips FIX_NE):
  strips FIX_T, FIX_T_AUTO, DO_T_AUTO, DO_GT_AUTO when set to a
  non-fixed value, with a WARNING.
- prepare_VADAT.py IT_SPECIFIER writer: bare FULL cycles rewritten to
  FULL(FIXT) -- SPECIFY_IT_CYCLE_V3 sets FIXED_T=.FALSE. for a bare
  FULL, which would bypass VADAT entirely.
- run_model_grid.py apply_fixed_j_workflow_overrides(): the J-release
  phase no longer injects FIX_T=F/FIX_T_AUTO=T, and the on-node VADAT
  rewrite in the generated j_release slurm no longer flips FIX_T.
- run_model_grid.py: warns when cmfgen.do_gt_auto=true is set (dead
  under FIX_T=T).
- prepare_IN_INTS.py: an explicit --gt-auto can no longer force
  DO_GT_AUTO=T when FIX_T=T.
- The stable convergence-mode profile tables no longer carry
  FIX_T=F/FIX_T_AUTO=T (dead under the policy).
- X_SOB.yaml corrected: FIX_T_AUTO=T/DO_T_AUTO=T -> F/F.

JFNK_FIXT auto-aligns to FIX_T in prepare_VADAT so it stays T with no
separate guard. Verified: apply_profile stamps T/F over a FIX_T=F
override for both profiles and both convergence modes; override strip
removes FIX_T/FIX_T_AUTO/DO_T_AUTO/DO_GT_AUTO=T requests while letting
benign FIX_T=T through; IT_SPECIFIER rewrite covers FULL and
FULL(FIXBA). Vault note added in 'VADAT keys.md'; policy documented in
top-level AGENTS.md section 23.

Addendum: two residual FIX_T_AUTO=T surfaces closed after the above:
- run_model_grid.py preserves VADAT verbatim on skip_final (OBS-only)
  resubmits; _pin_vadat_fixed_t() now rewrites FIX_T/FIX_T_AUTO/
  DO_T_AUTO/DO_GT_AUTO to T/F/F/F in any preserved VADAT, so a VADAT
  written before this policy cannot carry a stale auto-release flag.
  Verified on the live X_SOB 40d VADAT (FIX_T_AUTO T->F, DO_T_AUTO
  T->F, idempotent on second pass).
- scripts/submit_40d_staged.sh stage labels still claimed FIX_T=F for
  s2/s3; prepare_VADAT stamps FIX_T=T there anyway, so the labels were
  corrected rather than the behaviour.
- The dessart_m15_* validation configs still request FIX_T_AUTO=T/
  DO_T_AUTO=T in vadat_overrides; the strip pass removes them with a
  WARNING on every run (left in place, intent-vs-fact documented).

## 14-Oct-2026 - KNOP release-churn cycle and GS_FREEZE_LIM stop (arms done)

Two warm-restart arms off the it8 stall checkpoint ran to the designed
non-convergence stop (SOLVE_FOR_POPS GS_FREEZE_LIM=60 consecutive
frozen iterations -> STOP):

- stallpat (12188036 -> 12188938 -> 12190108; stall=0.999,
  MAXPASS=300, PRINT_KNOP_DEBUG): reached it121.  Both arms first
  RELEASED populations (stallpat it20, stallfm it13) after two
  consecutive mature sweeps -- first KNOP release on this model.
- stallfm (12188867 -> 12190107; freq_match=T, stall=0.98):
  reached it78.  Its bulk bad-frac declined monotonically
  3.1e-4 -> 1.8e-4 over the last ~30 its -- freq-match genuinely
  converges the bulk.

Mechanism (now quantified): release -> the second live-population
step (uncapped, demand ~1e21-1e98% on stiff levels CIII:2s2s / HI:17
at depth 69-70) displaces the field fixed point -> ~800-900K store
cells churn >50% -> GS_NCONV gate re-freezes -> store re-settles
~10x/iter under held pops -> re-release.  The cycle-2 spike was
~12-14% LARGER than cycle-1 in both arms -- a persistent limit
cycle, not a decaying transient.

Terminal blocker: after ~it60 (pat) / ~it18 (fm) the field locked
permanently frozen.  Bulk frac converged but GS_MAXCHG (the refresh
max residual, KN_RF_RESID at fg_j line ~5403) plateaued at 0.51-0.78
> GS_UNSAFE_TOL=0.5 -- a fixed-point inconsistency of the lagged
I_NEXT coupling at a handful of compressing hot-frequency cells.
Per-pass history shows NOT slow contraction but a fluctuating
plateau (~0.05, noisy), and the ordered sweep reads the live FLD
upwind (KN_ORD_SPV) -- so the intra-sweep GS spectral radius is ~1
on those cells; the lagged store is downstream of it, not the cause.

KN_DUMP/KNSTOR instrumentation (commit 4524716, cherry-picked into
cmfgen_jfnk_aa as 2d61d24, gated on PRINT_KNOP_DEBUG):
- KNDUMP: per-cell record in SWEEP_FLD_ORD (both march directions)
  for res>0.2 at pass>=8: pass, ray, depth row, freq idx, res,
  vold/vnew, GAM, JPR1/JPR2, NFM.
- KNSTOR: per-cell record in the in-sweep I_NEXT store update for
  |solved-stored|/stored>0.25: freq idx, ray, model depth, stored,
  solved, rel change, SIGMA_TRUE.
Caps: 2000 records per ray / per frequency.

Diagnostic arm dessart_m15_2p82d_KNOP_FORMAL_kndump warm-restarts
the it120 stallpat checkpoint for 4 its to classify the
non-contractive cells.

## 14-Oct-2026 (Devin session) - tree consolidation + per-cell P2 damping (commit 4754f8f)

Consolidation: `cmfgen_jfnk` is again the single active tree.
`cmfgen_jfnk_aa` (a cherry-pick mirror of this tree's solver commits
plus the OBS pipeline commits it lacked) and
`cmfgen_jfnk_stale_pre_1ec74d5` moved under `archive/`; the
install-dir `exe_frozen` cache was deleted (recreated per-run by
run_model_grid at submission).  All yaml `install_dir` entries now
point here.  Verified before the move: no `_aa`-unique files beyond
the exe cache; no pending/running job references the old path.

KN_ORD_P2DAMP (default F, commit 4754f8f): per-cell period-2
collapse in SWEEP_FLD_ORD.  Motivation: the kndump arm held
GS_MAXCHG at a stable 0.53 for it125-136 -- ~20 near-floor cells on
rays 3-6 flip-flopping ~2x/~0.5x per pass (lambda~-1), invisible to
the whole-ray KN_SGN_SUM detector (commit 635625c) because the
contracting bulk dominates its sign.  Significant cells track the
update-delta sign; KN_ORD_P2N (2) consecutive flips switch the
commit to the bracket midpoint (geometric when both endpoints are
positive - the observed flip is multiplicative so the midpoint is
the exact fixed point; arithmetic fallback otherwise); KN_ORD_P2REL
(3) flip-free passes release.  Fixed point unchanged; the nfm slot
needs no freeze since KN_ORD_SPV's upwind offset is static per
cell.  State: 3 int*1 arrays per ray (~3 B/cell), per-iteration
`p2damped writes=` counter prints only when on.

KN_GATE_DIVERGENCE (default F, dormant prototype): SOLVE_FOR_POPS
unsafe tier tests divergence (growth >KN_GATE_DIV_FAC over
KN_GATE_DIV_M sweeps, or non-finite) rather than the absolute
GS_UNSAFE_TOL.  No absolute backstop by design; owner opt-in only.

Validation arm `dessart_m15_2p82d_KNOP_FORMAL_p2damp` warm-restarts
the it136 kndump checkpoint with KN_ORD_P2DAMP=T.

## 14-Oct-2026 - Knop solver stability-audit fixes (code only)

Applied the outstanding items from the Knop-solver stability audit.
All changes sit inside DO_NONMONO_KNOP machinery; the stock path is
untouched.

- MOM refresh bulk exit (`MOM_RF_BULK_EXIT`, mom_j_cmf_nonmono_v1.f):
  per-frequency `MOM_GS_NHUGE` is now accumulated into `PASS_NHUGE`
  and the exit requires `PASS_NHUGE <= MAX(50,MOM_SW_NTOT/10000)` --
  the same bound the sweep gate's unsafe tier applies in
  SOLVE_FOR_POPS.  Previously a handful of cells drifting >1E6 could
  sit under a "bulk-converged" exit and have the locally exploded
  field retained.
- Retained-field delta loop: a non-finite entry on either side
  (JNEXT_SAVE still holds pre-heal poisoned slots) now sets
  `MOM_RF_DCHG=1E30` and counts as n5/n50, matching the sweep drift
  loop; before, a NaN silently slipped the significance floor or
  poisoned the MAX while the gate read the field as converged.
- New conditional pass-table rows: `huge-drift` (PASS_NHUGE) and
  `gam-floor` (cells at the `|1+sigma*K| <= MOM_GAM_FAC_MIN`
  resonance floor) so a refresh stall is attributable to resonance
  ill-conditioning vs marginal stiffness.  `MOM_GAM_FAC_MIN` moved
  to module scope for the count.
- Ray-side omega adaptation (SWEEP_FLD alternating, SWEEP_FLD_ORD
  ordered) no longer damps on a baseline pass still inside the
  warm-start seed transient: floors NPO>=5 (baseline NPO-2>=3) and
  NPO>=4 (baseline NPO-1>=3) respectively.  Recovery unfloored; the
  ORD sign-flip (NOSC) damping branch is floored too.
- SWEEP_FLD_ORD rescales the Aitken delta history `DLT_PRV` by
  OMEGA_RAY/T4 on an omega change, same as the existing AA_RH
  rescale -- the DEL_OLD/DEL_CUR pair no longer mixes two maps.

`reference/VADAT keys.md` gained the missing `MOM_RF_BULK_EXIT` row
and `reference/Diagnostics and OUTGEN.md` the two new grep targets.
Stock-path equivalence holds by construction: every changed line is
inside KNOP-gated refresh/ray-solve routines the stock path never
reaches.  `ifort -c` clean on both edited files and `./build.sh` succeeded
with the JFNK harness tests passing.  A `run_model_grid.py`
validation run is still required before this counts as validated
(per AGENTS.md section 6/7).

Submitted 14-Oct-2026 from committed state ef769d4e (frozen exe
hash 0d8fa5b4ecea4048): `dessart_m15_2p82d_KNOP_FORMAL_freqmatch`
(FG ordered sweep + AUDIT-4/5 + DO_NONMONO_FREQ_MATCH experiment)
as job 12204057 and `dessart_m15_2p82d_KNOP_MOM_probe`
(MOM_J_CMF_NONMONO_V1 warm-restart probe, covers the AUDIT-1/2/3
refresh path and MOM_TRIAL_LAM) as job 12204071, OBS 12204072.
Both on the interactive partition; results pending.

## 14-Oct-2026 (part 3) - Point-source irradiated 350 d SN A ejecta (PNT_SRCE)

New pipeline feature, no Fortran changes. `configs/point_source_sn_a.yaml`
(SN2024pba commit d958f034) submits `point_source_sn_a`: the
SN_A_g_T20k_Ni06 ejecta converted at 50.227 d and aged to 350 d by
`scripts/age_sn_hydro_data.py` (ballistic expansion R -> R + V*dt with
mass-conserving dilution, velocity monotonised at the 3950 km/s
reversal, Ni56 -> Co56 -> Fe56 decayed from the .xni initial profile).
An artificial inner point source heats the old ejecta:
`IB_METH=PNT_SRCE` + `PNT_SRCE=T` with `R_PNT=1.098D4 Rsun`
(7.6e14 cm, inside the 6.8e15 cm inner cavity) and `TEFF_PNT=0.7`
(L = 4 pi R^2 sigma T^4 = 1e42 erg/s = 2.61e8 Lsun, also stamped into
LSTAR for RLUMST).  Standard solver: SOL_METH=TRI (NUM_BNDS=3),
GLOBAL_LINE=BLANK, FIX_T=F + FIX_T_AUTO=T so the gas T responds to the
irradiation (new `cmfgen.solve_temperature` YAML opt-in gates the
hydro-T pin in prepare_VADAT.py/apply_profile and the skip_final
VADAT pin in run_model_grid.py).  `is_luminous_inner_boundary` now
counts PNT_SRCE so the dead-boundary guard does not reject it.

Frozen exe hash 0d8fa5b4ecea4048 (git f7562919) - stock PNT_SRCE path,
unchanged binary.  Jobs: CMFGEN 12204418, OBS 12204419 (afterok).
ND=140, tau_es(inner)=0.38 - optically thin ejecta, source flux mostly
transmitted/reprocessed.  Results pending.

## 14-Oct-2026 (part 4) - Non-mono audit v2: F4 fixes, F1/F2/F5 false positives

Commits 58a496dc..c3fe4c61 on development.  Re-verified all eight
findings of the audit-v2 brief against the source; four were false
positives or intentional design, documented in
cmfgen_obsidian/audits/NONMONO_AUDIT_V2.md.

Dispositions:

- F1 (refresh spacing): FALSE.  KN_DLN_P(LK)=dLOG_NU captures
  ln(nu_{k-1}/nu_k), the backward span expanding cells need --
  identical to in-sweep DLN_C=dLOG_NU.  The declaration comment was
  wrong (said forward); comment corrected, code untouched.
- F2 (missing 1/mu in Sobolev tau): FALSE, physics.  tau_S uses
  dv_los/ds=(V/R)(1+sigma*mu^2), which is finite for grazing rays;
  no 1/mu belongs.  Convention comment added to sobew/sobjbar_sim.
- F3 (SIGMA clamp destructive): BY DESIGN.  Flags are mutually
  exclusive; SIGMA_TRUE carries unclamped physics, SIGMA stays
  infrastructure-safe and is regenerated from R,V before re-clamp.
- F4 (feature-insertion gap): REAL, fixed (58a496dc).  U_GAP_FAC=5E-2
  floor rejects insertions leaving <5% of a uniform cell; post-build
  strict-decrease assert + tightest-cell report.
- F4b: N_UNDEF>0 now STOPs under a non-mono flag; opt-out key
  ALLOW_SIGMA_UNDEF (4385154a).
- F5 (dead J): FALSE.  J is the cubic-segment index in the V
  interpolation loop.
- F6 (K_ON_J 1E-4 clip): mechanism disproven (clip can only
  over-flag compression); added MOM_N_SIGBIG census printed as
  nsigT>1e4 to show whether |sigma|>1E4 ever occurs (9aa1fea5).
- F7 (hysteresis stale reads): FALSE.  Hysteresis-only cells have
  1+sigma*mu^2>0 for all mu -> expanding branch -> fresh reads.

Verification: ./build.sh clean (all three exes fresh, harness T1-T3
pass).  unit_tests: test_rv_snhyd 36/36 (new Steffen/U-coordinate
replica matches Fortran to ~5E-12; tuned adjacent-features case
proves the floor rejects a 0.036-cell gap), test_sigma_true 13/13,
test_vadat_keys 34/34.  Pre-existing failures unchanged:
test_interp_guards (3 checks, fails on reference routines too),
test_jfnk_solver (missing external binary).  Default path: all new
behaviour gated on DO_NONMONO_* or VEL_LAW=14; sobew/sobjbar_sim/fg
edits are comment-only.

## 14-Oct-2026 (part 5) - Baseline submission + Phase 0 solver-architecture audit

New mission brief (population-convergence modernization, Phases 0-5).
Directive: submit models with USE_LAM_ES=T; JFNK stays out of scope
and disabled.

Baseline run: configs/dessart_m15_2p82d_SOB_l14_base.yaml - cold
start, VEL_LAW=14, DO_NONMONO_SOB=T, USE_LAM_ES=T, JFNK off,
DIF=T/IB_METH=DIFFUSION (the ZERO_FLUX first attempt was correctly
rejected by the submission watchdog: inner tau ~37 with nonzero
LSTAR cannot discard flux).  SN_HYDRO_DATA staged from the proven
m15 2.82d file (Kappa column present for TAU_OPT=1).
Submitted from committed state c3fe4c61 (frozen exe): CMFGEN job
12205705, OBS job 12205706.

Phase 0 audit (static, complete): mapped SOLVE_FOR_POPS /
SOLVEBA_V13 / CMF_BLKBAND_V3 / DO_NG_BAND_ACCEL_V4 / NGACCEL and the
CMFGEN_SUB dispatch into a single note -
cmfgen_obsidian/audits/PHASE0_SOLVER_ARCHITECTURE.md.  Covers: the
per-iteration call graph with TUNE labels, the BA matrix lifecycle
(COMPUTE_BA / WRBAMAT / BA_MATRIX_STORED / N_ITS_TO_FIX_BA /
VAL_FIX_BA / MAXCH_SUM / CNT_FIX_BA / BCD factor cache), every POPS
mutation path (commit, masked/partial/full rollback, GS hold,
smoothing, AV averaging, NG, grid revision, LTE fill), the
Lambda/full state machine (hysteresis, BA_ENTER_RAW demand gate,
BA_DAMP_CAP, BA_SAT_LIM, T-release phasing), the resacc adjudication
(RES2 bulk channel, per-column mask, quarantine, soft rejects,
trust radius), NG window/cleanliness rules, and which metrics are
control inputs vs diagnostic-only.

Key structural findings for Phase 1-4 planning:

- BA reuse is already fully implemented and guarded (FIX_BA scheme
  + stored-flag + LU-cache presence checks).  Whether it wastes
  time is a measurement question, not an obvious gap.
- All known stale-matrix hazards are invalidated at the decision
  layer already (grid revision, rejection, NG kick, restart).
- Resacc already provides per-depth residual history (COL_RES /
  BAACC_COL_PREV / BAACC_COL_MIN) - the substrate a Phase 2
  depth-aware controller would have to reuse.
- BA_TRUST_CUR already zeroes over-trust columns; the documented
  reason is that a jointly-solved over-trust direction is not
  partially trustworthy - any future per-depth scaling must argue
  against that, not just add damping.
- NG history integrity is guarded by NG_CLEAN against every known
  contamination source; Anderson (Phase 3) inherits the same list
  plus grid revision.

Baseline metric extraction waits on 12205705 producing iterations
(~74 s/iter on viper07, SOB/LAM_ES path confirmed in OUTGEN).

=== 15-Oct-2026: Phase 1 implementation - DO_BA_RES_REUSE ===

Opt-in residual-driven BA-matrix recall (Phase 1 of the population
convergence brief).  New VADAT keys: DO_BA_RES_REUSE (F),
BA_REUSE_TOL (0.9), BA_REUSE_NSTALL (2).  Requires DO_BA_RES_ACC=T
(warns and is inert otherwise).  SOLVEBA_V13 exports BA_JUDGED /
BA_JUDGE_RATIO per armed-step adjudication (post/pre RES2 ratio;
nonfinite or saturated post-step residual counts as a stall).
SOLVE_FOR_POPS recalls the stored matrix while judged recall steps
contract (ratio <= TOL), past the FIX_BA demand/counter triggers;
BA_REUSE_NSTALL consecutive stalled recall judgments force one
rebuild.  Rejection (BA_REJ_THIS/BA_PART_WAS), MAXCH_SUM drift
(SUM_FIRED), grid revision (BA_MATRIX_STORED cleared), JFNK calls,
CON_SCL_T/=0, and LAMBDA next-call all still gate/veto.  Recall
assignment is a hard L_FALSE (COMP_BA=T in production VADAT would
make COMPUTE_BARDIN a no-op).  Diagnostics: ITER_REASON
'ba_reuse_stall', OUTGEN 'BA_RES_REUSE:' lines (CONV_DIAG_LVL).

Verified so far: all four edited sources compile clean; VADAT test
37/37 (key parse/land, inert-without-resacc warning, defaults).
Clean build running in an isolated worktree because a parallel
session is editing plane/{fg,mom}_j_cmf_nonmono_v1.f (MOMJNXT_GS
VER=4 reader, currently non-compiling) and building in the main
tree with /dev/shm/arwrap_sb (GPFS ar-mmap SIGBUS workaround -
root cause of the earlier phantom missing-.o and malformed-
archive build failures).

Baseline 12205705 at it 33: first FULL call at it 31 cost 378 s
and ended HELD (raw demand spiked 5e1 -> 1e7 on the fresh matrix);
LAMBDA warm-up otherwise 74.5 s/iter.

=== 15-Oct-2026 (audit session): refresh-exit coherence + MOM rollback ===

Commits aab11c44 (guards) and 9cf8eb78 (coherence/rollback) on
origin/development, built and verified in the main tree with the
/dev/shm/arwrap_sb ar wrapper (GPFS ar-mmap SIGBUS workaround).

aab11c44: staged every compound test that relied on .AND./.OR.
short-circuit evaluation with an allocation or bounds operand --
COMP_DEPTH_IDX hysteresis, KN_FREQ checkpoint validation, the two
NU_F(JF+1) remap marches, the KN_RES_PASS backward scan and record
guards, the KN_HOT_CNT run merge, MATCH_FREQ's KN_FREQ(L+NNF) span,
FACR_KONJ reshape (plus per-array DEALLOCATEs), MOM_DAMP_NM shape
check, MERGE(SIZE(K_ON_J)) in the MOM guard report, and
BA_QUAR_MASK(JK=0) in SOLVEBA's STEQ_ADV scan.

9cf8eb78: KNOP_RAY_SOLVE_V1 now returns RESID=1.0 on work-array
allocation failure (was 0.0 = converged); MOM_KNOP_REFRESH_V1
snapshots the pre-refresh J/H store after the poisoned-entry heal
so a retained pre-refresh field is finite; outer-AA history
allocation is all-or-nothing with a persistent MOM_AA_OUT_FAIL
latch keyed on (ND_NM,KN_M_NCF,MOM_AA_OUT_M) and a logged memory
estimate; new RESTORE_MOM_JNEXT_LATEST mirrors
RESTORE_FG_INEXT_LATEST so BA rejection rolls back the MOM
J_NEXT/H_NEXT store (refreshing the outer-AA pair and clearing
NHIST_OUT), wired into both SOLVEBA_V13 rollback sites under
(DO_NONMONO_KNOP .OR. DO_NONMONO_SOB).

Verified: make DEBUG=nodebug all clean (recompiled fg/mom/solveba),
symbols in exe/cmfgen_dev.exe, tests/f90/build_and_run.sh PASS.
The DO_BA_RES_REUSE hunks in solveba/control_variable_mod/
rd_control_variables/solve_for_pops/unit_tests belong to the
parallel session and were left uncommitted.

Runs submitted from committed state HEAD=a5469874 (Phase-1 commit
3bc23cde, pushed to origin/development by the sibling session's
sync; sibling commits aab11c44..a5469874 carry the MOM lagged-store
rollback + refresh telemetry work):

- dessart_m15_2p82d_SOB_l14_resacc   CMFGEN 12207306 / OBS 12207307
  Phase-1 CONTROL: DO_BA_RES_ACC=T, DO_BA_RES_REUSE=F.  Isolates the
  resacc overhead + rollback machinery vs the 12205705 baseline.
- dessart_m15_2p82d_SOB_l14_resreuse CMFGEN 12207308 / OBS 12207309
  Phase-1 EXPERIMENT: DO_BA_RES_ACC=T + DO_BA_RES_REUSE=T
  (BA_REUSE_TOL=0.9, BA_REUSE_NSTALL=2).

Success criterion per the brief: BA matrix work measured at ~5% of
per-iteration wall time (TIMING file: BA assembly 3.7% + LAPACK
factorization 1.5% of GIT), so the achievable ceiling of Phase-1
alone is a few percent of runtime.  The run pair measures builds vs
recalls, rejection counts, stall rebuilds, and physical equivalence
of the final pops/spectrum.  scripts/iter_phase0_metrics.py counts
the BA_RES_REUSE diagnostics.

## 14-Oct-2026 (Devin audit cont.) - JFNK applied-state verification (commit 6740dd5)

Audit finding (Phase 7): the Newton accept tests in JFNK_SOLVE_CORR_V1
all ran on unclipped iterates; the commit then clipped each component
to MAX_JFNK_COR / MAX_dT_COR and applied the result without ever
evaluating the residual at the state actually committed.  A heavily
clipped or leaked-state commit could therefore degrade the physical
state while still reporting SUCCESS=T.

Change (new_main/mod_subs/jfnk_solve_corr_v1.f): after the clipped
correction is written to POPS, the driver reconstructs
x_app = log(POPS/POPS_BASE) (frac form under [JFNK_LOGX]=F), runs one
JFNK_RESIDUAL_V1 eval at the applied point, masks fixed rows, strips
the PTC term, and compares |F(x_app)|_w with the entry |F|_w.
Degradation or a broken eval (jfnk_resid_nbad, non-finite x_app)
rolls POPS/ED/T back to POPS_PREV and falls through to the existing
caller-mode SOLVEBA fallback; the reject increments the stagnation
counter.  Frozen tier only - a coupled-tier check would cost a full
formal sweep, and coupled commits already carry
JFNK_UNCONVERGED_COMMIT.  No new VADAT key.

Verified: ifort syntax check clean; incremental make rebuilt
jfnk_solve_corr_v1.o into libsubs.a and relinked cmfgen_dev.exe;
new OUTGEN strings present in the binary;
tests/f90/build_and_run.sh: 36 PASS / 0 FAIL (solver 21, interface 9,
harness T1-T3 + 3 spot checks).  Prior surrogate suite (unit_tests/)
unchanged - it exercises jfnk_solver_v1.f90, not the driver.

## 14-Oct-2026 (Devin audit cont. 2) - TAUTO commit-range fix (commit 401eddd0)

Audit finding: under JFNK_TAUTO=T the residual, mask, rebuild, and
CMAT paths all treat T rows inside JFNK_FIX_T_D_ST:END as live
unknowns, but the commit block's POPS(NT,K) write and the REACH
clipped-step pricing skipped the range unconditionally - the damped
T corrections the solver computed for the released range were
silently dropped.  Added .NOT. JFNK_TAUTO_ACTIVE to both skip tests
(driver, jfnk_solve_corr_v1.f), matching the six guarded sites in
jfnk_cmfgen_interface_v1.f.  No behaviour change when JFNK_TAUTO=F.

Verified: archive member byte-identical to fresh compile; exe
relinked; tests/f90 unchanged (driver not unit-covered).
JFNK probe 12209146 (HEAD 6740dd5b frozen exe) queued to exercise
the applied-state verify eval on the m15 SOB model.

## 14-Oct-2026 (Devin audit cont. 3) - Live runtime results, both probes

Probe 12209146 (m15 SOB jfnk_probe, frozen exe at 6740dd5b,
JFNK_REACH_GATE=1D-3): resumed at it 461, ran to its 2h debug wall
at it 470 and finished cleanly; OBS 12209147 launched.  Six frozen
admissions, F(0) rel diff 7.68E-18, all six reach-skipped at demand
11-40% and fell back once to SOLVEBA - bounded admission cost
(~75-490 s/iter), no double correction, demand kept falling
(39% -> 7.8%) under BA.  Deterministic fallback under the intended
configuration demonstrated end-to-end.

Verify 12210267 (m15 SOB jfnk_verify, exe_frozen
655f2e4f = byte-identical to exe built at 401eddd0; resubmitted to
the small partition after AssocMaxJobsLimit on interactive):
warm-restarted from the probe checkpoint at it 469 with
JFNK_REACH_GATE=0 to force Newton.  Three forced admissions
(it 470-472) all hit gmres_stagnation_abort (info=3): Krylov res
flat (2.1706E2 -> 2.1687E2 in 4 its).  Stagnation counter reached
3/3, cooldown disarm printed "paused (cooldown 1 of 20)" and
caller-mode SOLVEBA ran once per iteration - no double correction,
clean fallback, re-arm scheduled after 20 eligible calls.
F(0) rel diff 2.8E-18 / 4.2E-18 on both admissions.

State-dependent operator finding (m15, it 470): preconditioner
audit ||Mz-v||/||v|| = 1.25E9 for v=F while unit-vector solves are
clean; FD-vs-BA JVP |diff|/|anal| = 46.4 at this state vs 6.6E-5
in the earlier probe and ~1E-2..1E-4 in zl14_fixt_sobes.  On
near-singular states (rcond ~1E-5 at depths 14-17, 42604
regularized near-null rows) the BA-LU preconditioner does not
represent the true FD Jacobian, so GMRES cannot progress; this is
the diagnosed operator limitation, not a plumbing defect.  On
better-conditioned states (zl14 run, commits at it 822) the same
machinery produces commits with PC accuracy ~1E-4.

The new applied-state verify eval (6740dd5b) is live in the verify
binary but has not yet fired - it only runs on a commit, and the
m15 admissions all stagnated before any step was accepted.

## 14-Oct-2026 (Devin audit cont. 4) - JFNK direction gates (commit d6d23ee6)

New opt-in VADAT keys JFNK_COS_MIN, JFNK_JZ_F_MAX,
JFNK_CLIP_FRAC_MAX (all default 0 = off) turn already-computed
reach-block metrics into admission defers.  Stage A tests the
stored-BA z/Jz metrics for free; stage B repeats them on the FD
image (one residual eval) - catches states like m15 it470 where the
stored operator reports reach while the true Jacobian maps z
orthogonally.  Defers take the info=-2 path (caller-mode BA once,
not stagnation); the -2 message now reads "admission gate".  The
precheck block entry condition extended to run when any of the
four gates is armed; the REACH_VAL test is now explicitly gated on
JFNK_REACH_GATE>0 (previously implicit via the entry condition).

Verified: ifort prod compile clean; incremental make rebuilt all
three objects (control_variable_mod, rd_control_variables,
jfnk_solve_corr_v1) into the archives; exe relinked, guard strings
present; tests/f90 36 PASS / 0 FAIL.

## 14-Oct-2026 (Devin audit cont. 5) - Direction gates live (d6d23ee6)

Guard probe dessart_m15_2p82d_SOB_jfnk_guard (CMFGEN job 12211053,
small partition, frozen exe at d6d23ee6 = byte-identical to current
exe) warm-restarted the same probe checkpoint as the verify run.
First FULL admission (it 473): stored-BA metrics at admission were
REACH clip max|z|=8.1E11 nclip=41101, |Jz|/|F|=5.0E3,
cos(F,Jz)=4.83E-3 - all three armed gates would fail; the stage-A
cos test fired first:

    JFNK GUARD DEFER: stored-BA non-descent cos(F,Jz)= 4.83E-03
                      < cos_min= 1.00E-02
    JFNK: admission gate skipped Newton; retaining caller mode
          (not counted as stagnation)

The admission cost one PC apply + one block matvec instead of a
doomed GMRES cycle; caller-mode BA ran once.  Compare the ungated
verify run: same checkpoint family, every admission burns a Newton
attempt, stagnates, and hits the 3-of-3 disarm.  All three keys
echoed in OUTGEN.

Verify job 12210267 update: re-armed at it ~492, stagnated 3 more
times, disarmed again (cycle 2, 40-call wait).  GMRES stagnation is
deterministic on this state - the preconditioner never covers the
HI outer-depth residual rows.

## 05-Nov-2026 (Devin audit cont. 6) - NG window fix (64ab8b82)

Convergence investigation of the slow verify/guard tails found the
outer-loop limiter: NG_CLEAN was reset by JFNK_THIS, the admission
flag, which is true on every eligible call while DO_JFNK is armed.
Deferred/cooldown/fallback calls write ordinary BA records yet still
zeroed the window, so NG acceleration could never reach its
4-clean-record requirement (NG done=F throughout; next NG deferred
indefinitely).  Population demand was contracting linearly at
~0.99/iteration -- slow fixed-point, not stagnation -- but the
endgame accelerator was suppressed by bookkeeping.

Fix: new module flag JFNK_COMMITTED (control_variable_mod), set from
APPLY_JFNK_UPDATE at the end of the JFNK branch in
jfnk_solve_corr_v1.f (verify-rejects, catastrophic, and all defer
paths leave it false); reset each call in solve_for_pops alongside
JFNK_UNCONVERGED_COMMIT.  NG_CLEAN now resets only when a Newton
correction was committed.  Deferred calls count toward the window
like the BA iterates they are.

Verified: incremental make rebuilt all three objects into the
archives; exe relinked; tests/f90 36 PASS / 0 FAIL.  No VADAT key;
behavioral change is internal bookkeeping only and matches the
comment's original intent.  Pushed to origin/development.

## 15-Oct-2026 (Devin audit cont. 7) - BA direction check DO_BA_DIR_CHK

Divergence audit follow-up: the residual-acceptance net (DO_BA_RES_ACC)
judges a committed BA step one call late and cannot distinguish a bad
direction from an oversized one - the it509 rebound (120x at +50% cap)
showed the direction was wrong, but only after the fact.  The JFNK
precheck already computes exactly the metrics that price the
linearization at the state the BA step uses (stored-BA cos(F,Jz),
|Jz|/|F|, REACH, then FD-vs-stored cos(F,A_FDz), |A_FDz|/|F|,
cos(BAz,A_FDz), row anatomy, COLPROBE), but only inside a Newton
admission that the demand gates usually defer.

New opt-in VADAT key DO_BA_DIR_CHK (default F, requires DO_JFNK):
SOLVE_FOR_POPS re-arms JFNK_THIS on every structurally eligible call
(past warmup, non-LAMBDA, not GS_FREEZE, USE_LAM_ES, not an rd_lam
probe), bypassing the demand gates; JFNK_COUPLED_THIS stays F (frozen
tier - the operator question needs no formal sweep per eval).  The
precheck entry condition gains .OR. DO_BA_DIR_CHK so the metrics run
with all four gates disarmed, and a forced defer at the block exit
(JFNK_REACH_OK=F -> info=-2) skips Newton unconditionally: the same
call's SOLVEBA_V13 fallback applies the real correction, unchanged.
Caller-mode BA once per dirchk call, not counted as stagnation;
JFNK_COMMITTED stays false so NG window bookkeeping is unaffected.

Audit context recorded while implementing: every m15 validation run
(fixt_noes/fixt_pop/fixt_sobes, verify, guard, probe, KNOP probes)
runs FIX_T=T - T is pinned to the hydro profile, FIX_T_AUTO/DO_T_AUTO
overrides are rejected by the pipeline guard.  The near-singular
operator state (it470: rcond~1E-5, FD-vs-BA=46.4, cos(F,Jz)=4.8E-3)
was therefore measured with all T rows masked - the inconsistent
linearization lives in the population rows, not the T/ED coupling.
Contrast zl14_fixt_pop (healthy state): REACH=1.0, cos(F,Jz)=1.0,
cos(F,A_FDz)=0.96 - the stored BA IS the true Jacobian there.
A free-T A/B on m15 cannot be pulled from existing logs because
free-T has never run on this model; the T-coupling question is moot
for the observed divergence.  Schur stays JFNK-only for now.

Vault: solvers/JFNK solver.md updated.

## 15-Oct-2026 - JFNK_TAUTO must not release a VADAT FIX_T pin

Found on the DO_BA_DIR_CHK probe (job 12220687, it473): the JFNK
CONTEXT block printed `T free to vary: T` and the `free-T JFNK` note
on a run with `FIX_T=T` in VADAT.  Two separate defects:

1. Print semantics.  `T free to vary` reported the IN_ITS schedule
   flag FIXED_T, which a `FULL` cycle token sets F
   (specify_it_cycle_v3.f) regardless of VADAT FIX_T.  The banner now
   reports the effective state `.NOT. FIXED_T .AND. .NOT. RD_FIX_T`,
   the free-T NOTE prints only under the same condition, and the
   `T held fixed over depths` line prints whenever MOD_FIX_T_D_ST>0
   (previously suppressed when FIXED_T=T, i.e. exactly when pinning
   was real).

2. Behaviour.  `JFNK_TAUTO_ACTIVE` keyed on `.NOT. JFNK_FIXED_T_LOC`
   (the schedule flag) only.  With `JFNK_TAUTO=T` plus VADAT
   `FIX_T=T`, TAUTO disengaged the depth-range T mask at all six mask
   sites - a committed Newton step would have applied damped T
   updates under a global FIX_T pin.  JFNK_TAUTO_ACTIVE now also
   requires `.NOT. RD_FIX_T`: a VADAT-level FIX_T is an absolute pin
   TAUTO may not release (TAUTO still governs VARFIXT ranges and
   schedule-freed T as before).

Harmless in the dirchk run itself (Newton is skipped after the
metric suite), but the invariant "FIX_T=T and FIX_T_AUTO=F implies T
always fixed" is now true inside the JFNK branch as well.

Build: ./build.sh exit 0, 36 tests PASS.  Commit pending.
Vault: solvers/JFNK solver.md TAUTO section updated.

## 15-Oct-2026 (2) - dirchk results on converging state; column defect found; lc probe submitted

DO_BA_DIR_CHK first data (job 12220687, m15 jfnk_probe warm restart,
it473-481, commit 6f8e5f29 binary).  On the converging state:

- F is 98% concentrated in HI packed-Rydberg rows lev 22/23 (n=16/17)
  at depths 65-70, collision-dominated (coll/F~1.0).  |F|=3.4e13.
- z=M^-1 F: 42670/128730 unknowns population-floor-frozen, ~41k more
  components at the +-2.3 commit clip -> ~65% of the vector dead or
  saturated.  |Jz|/|F|~5e3, cos(F,Jz)~+-0.1, REACH -1e6..-1e8 every
  iteration.  The frozen-J Newton problem is formulation-bound
  ill-posed here; direction gates correctly refuse.
- cos(BAz,A_FDz)=1.0000 and FD-vs-anal 1.2e-2 on the clipped step:
  the stored operator is faithful ALONG the committed direction.
- F(0) rel diff 1.06e-17 (evaluator == outer STEQ exactly); COLPROBE
  repeat-vs-pristine 0 (deterministic evals); PC healthy.

Localized defect (persistent, 6/6 iterations): COLPROBE topF idx
119558 = HI lev 23 (n=17 packed) at K66: |FDcol|=0.925 vs
|BAcol|=1.079, diff/BA=0.1696.  Element (row 21,col 23): FD -4.18e17
vs stored +5.1e15 coll (+3e9 rad).  No radiative lines exist between
packed levels in HI_F_OSCDAT, so the missing coupling is
non-radiative.  Column sums match (~-3.4e18): BA overweights the
diagonal outflow by ~8.6e17 while under-booking downward collisional
inflows into adjacent-n rows (steq_multi_v10.f:244 writes
BA(J,K)+=ED*CNM(K,J) for all pairs - investigate whether the packed
CNM bookkeeping misplaces downward n->n-k derivatives or whether an
n23-dependent outflow term in F_21 is omitted).  Latent: does not
load the committed direction (cos=1.0), but sits on the dominant
residual column.

JFNK-CPL tier experiment is unavailable on this model: JFNK_COUPLED
is KNOP-only (rejected under DO_NONMONO_SOB at VADAT read,
rd_control_variables.f:1252), and the KNOP/MOM pathway diverges
structurally on m15 (SIGMA vs SIGMA_TRUE).  The relax-J reachability
test is blocked by two independent gaps.

New run: dessart_m15_2p82d_SOB_jfnk_dirchk_lc (job 12221660, OBS
12221661) - dirchk on the diverging resacc limit-cycle checkpoint
(m15_2p82d_SOB_l14_resacc ~it169+, FULL->REJECT->LAMBDA cycle,
demand ~2 then 1.5e5 on rejected steps).  Frozen binary 4b7aaa24
(includes the TAUTO/RD_FIX_T banner fix).  Question: does
cos(BAz,A_FDz) break down on the diverging state (operator failure)
or stay ~1 (formulation-bound)?

Also fixed this session: JFNK_TAUTO_ACTIVE gated on RD_FIX_T +
banner prints effective T state (commit 4b7aaa24, pushed).

## 15-Oct-2026 (3) - dirchk_lc verdict: operator failure confirmed; BA_MAXCH_FAIL guard

dirchk on the diverging resacc limit-cycle state (jobs 12220687,
12221660) answered the motivating question: cos(BAz,A_FDz) breaks
down on the diverging state.  Clipped FD gives cos=0.87-0.9999
(stored BA faithful along the huge-|z| components), but with
|z|>1 rows zeroed cos collapses to 0.008-0.17 and |A_FD z|/|F| hits
1e2-4e5 - the true Jacobian amplifies the proposed step.  The
runaway rows are OIII/OIV excited levels at K60-69 demanding
z~1E9-3E10 on n~1E-5 populations.  Operator failure, not
formulation-bound.

The same logs exposed the commit-path hole: DO_LAM_AUTO's stock
handoff (cmfgen_sub.f, hardcoded MAXCH_CTRL<50) fired at it172 on a
47.7% demand, and the it173 full step committed a 400% update on a
raw 1E7 demand (largest decrease 1.0E+07% at CV:1s_2s_1Se d16) with
level-check=F and zero rejects - SCALE_OPT=MAJOR clamps but never
refuses, and DO_BA_RES_ACC=F in this probe means nothing vets the
committed state before it feeds the next iteration.  The mode
controller's demote is post-commit.

Changes (this commit):
- BA_MAXCH_FAIL (VADAT, default 1.0D30=off; the key existed in the
  template but was never registered).  SOLVEBA_V13 now rejects a
  committed non-LAMBDA step same-call when MAXCH_APP exceeds the
  limit or any population is non-finite: POPS restored to the
  pre-step POPS_OLD, BA_STEP_REJECTED raised, 1.5E5 sentinel -
  SOLVE_FOR_POPS forces a BA recompute and demotes via the existing
  BA_REJ_THIS path.  No lagged-field rollback needed: nothing
  consumes the rejected pops within the call.
- DO_LAM_AUTO handoff (cmfgen_sub.f) under DO_NONMONO_SOB/KNOP now
  additionally requires CNT_STABLE_LAM >= N_STABLE_LAM (already
  counted during RD_LAMBDA at MAXCH_CTRL < LAM_EXIT).  The
  cnt_lam_escape / BA_PROBE_LAM probes remain the bounded full-step
  channel, now protected by BA_MAXCH_FAIL.

Test run: dessart_m15_2p82d_SOB_jfnk_guard_lc - identical to
dirchk_lc (same resacc warm restart, DO_BA_DIR_CHK=T) plus
BA_MAXCH_FAIL=1.0D2.  Expected: the it173-equivalent catastrophic
step is rejected same-call and demoted instead of committed.

## 15-Oct-2026 (4) - JFNK_MAX_SEC wall-clock abort; converter fail-fast gate

Commit 5600205c (cmfgen_jfnk, pushed to origin/agent/jfnk-conv-2026-10-15):
JFNK_MAX_SEC (VADAT, default 0=off) bounds JFNK wall time per
activation.  The JFNK_MAX_EVALS counter alone could not fire on the
it191 failure mode - that admission spent ~615 s in precheck and
diagnostics while performing 0 residual/JVP evals.  JFNK_INIT_V1
stamps JFNK_T_START via OMP_GET_WTIME and clears JFNK_SEC_HIT;
JFNK_RESIDUAL_V1 returns the same 1.0D20 sentinel residual as the
eval budget once the elapsed limit passes, and JFNK_JVP_V1 treats a
time-exhausted eval identically so GMRES aborts instead of building
a Krylov vector on the sentinel.  Build clean (EXIT=0), all three
exes present.

Parent tree (uncommitted-in-git, scripts/prepare_SN_HYDRO_DATA
driver.py): fail-fast final-grid gate before write_sn_hydro_data.
Refuses to write when the final point count is outside
nzon_target*(1 +- grid_count_tol, default 0.5) or max|sigma| exceeds
--max-abs-sigma (default 60; the -70..-108 shocks in recent
snapshots are pathologies, not physics).  --allow-grid-mismatch
downgrades to a warning.  Synthetic-grid test still owed.

## 15-Oct-2026 (5) - guard_lc live evidence: BA_MAXCH_FAIL fires

Job 12223135 (dessart_m15_2p82d_SOB_jfnk_guard_lc, still RUNNING)
gave the first in-model evidence for the same-call reject:
- it203: "BA_MAXCH_FAIL: rejecting committed step -- committed max
  change 4.0000E+02 vs limit 1.0000E+02".  POPS restored, BA demoted
  to LAMBDA_TFIX (ba=REJECT), run continued.
- MAXCH trajectory after the reject: 4.056 (it200) -> 3.552 ->
  3.202 -> 2.979 (it205) - the guard prevented the catastrophic
  400% commit that dirchk_lc showed could fire unvetted, and
  Lambda contraction resumed.
- JFNK still defers each iteration on the demand-vs-freeze gate
  (jfnk demand ~3-4 vs frozen threshold 1E30): admission policy
  working as configured, not evidence of a JFNK failure.

## 15-Oct-2026 (6) - maxsec interactive run submitted (job 12224526)

dessart_m15_2p82d_SOB_jfnk_maxsec: clone of guard_lc plus
JFNK_MAX_SEC=1.2D2, submitted via run_model_grid.py with
submission: debug -> interactive partition (job 12224526, viper07,
1h59m wall).  Frozen exe git 5600205c.  Same resacc warm-restart
seed; gates disarmed and DO_BA_DIR_CHK=T so the full precheck suite
runs - a clean target for the 120 s wall budget to cut in.  VADAT
verified to carry 1.2D2 [JFNK_MAX_SEC].  OBS chained as 12224527.
Prior general-partition submissions 12224481/12224482 cancelled
before the partition switch.

Converter note: the fail-fast gate bound is --grid-sigma-max
(default 60), NOT --max-abs-sigma - that name is taken by the AMR
insertion threshold (default 5.0).

## 15-Oct-2026 (7) - JFNK_MAX_SEC reverted at user request

User judged a wall-clock cap inside JFNK too restrictive; commit
78e61516 reverts 5600205c in full (control var, VADAT key, module
state, residual/JVP checks).  grep confirms no JFNK_MAX_SEC /
JFNK_SEC_HIT / JFNK_T_START references remain.  Clean rebuild
(EXIT=0), all exes relinked; pushed to
origin/agent/jfnk-conv-2026-10-15.  maxsec config updated to drop
the key - it is now a plain guard_lc-equivalent probe of
BA_MAXCH_FAIL + dirchk cap + instrumentation.  Interactive test job
12224526 and chained OBS 12224527 were cancelled before the revert.

## 15-Oct-2026 (8) - replan runs submitted (Claude-suggested order)

Two new runs following the external review:

- dessart_m15_2p82d_SOB_l14_prod (FINAL 12225057, OBS 12225058):
  Lambda-only production lane seeded from the converged verify
  checkpoint (IREC 838, MAXCH ~5.4E-2 vs EPS_TERM=5.0D-02).  DO_JFNK=F,
  DO_NG=T, DIF=T + IB_METH=DIFFUSION restored (Luc ran DIF=T; the
  probe ZERO_FLUX was diagnostic-only).  Chained OBS retires the
  OBS_FIN_10 / OBS_FREQ coverage risk early.
- dessart_m15_2p82d_SOB_l14_purelam (FINAL 12225061, OBS 12225062):
  A/B arm from the same resacc limit-cycle seed as dirchk_lc -
  identical environment minus interventions (DO_JFNK=F, DO_NG=F,
  DO_BA_DIR_CHK=F).  Tests whether the period-3 ~1.6% demand cycle
  is self-inflicted by the interleaved BA/Ng pattern.

EPS_TERM semantics confirmed in new_main/cmfgen_sub.f (~line 5025):
termination fires on a single iteration with demanded MAXCH (%) <
EPS and schedules one final loop; units are percent, demand-side.
The verify run already satisfied it at it838-839.

## 15-Oct-2026 (9) - stiff-mode structure identified in dirchk_lc z

REACH-FD top-|z| tables (dirchk_lc OUTGEN): the dominant step
components are entire ion blocks moving together, not isolated
levels.  First admission: 10 OIV levels at K68 all carry
z ~ -3.3E12 (same sign, same magnitude) - a uniform ion-block
shift = near-null ionization-fraction sum mode vs the neighbouring
stage.  Second admission: NeIII/OIII at K64-65-69, z ~ +-1E10-11
alternating by depth, n ~ 1E-4..1E-3.  Supports the
sum-vs-difference hypothesis: the stiff direction is collective
(ion fraction), so a diagonal/LM ridge cannot fix it;
block elimination of the ionization-balance mode (or PE
elimination of the stiff pair) is the indicated fix.

Also observed: prod run (12225127, interactive) restarted from the
verify IREC-838 checkpoint with DIF=T now reports inner luminosity
4.27E8 (LSTAR) and outer 6.15E9 vs the ZERO_FLUX probes printing
0/0 - the boundary choice moves the luminosity diagnostic off the
floor immediately.

## 15-Oct-2026 (10) - live A/B evidence from the replan runs

prod (DIFFUSION boundary): warm restart at it850 re-measured demand
47.9 (from 3.46 pre-restart - checkpoint re-entry perturbs the SE
residual), then contracted 47.9 -> 2.32 over five Lambda iters.  At
it856 Ng fired on a 1.79 demand; it858 demanded 13.2 - a 7x Ng
rebound, caught live by the new demand instrumentation.  Scheduler
pushed the next Ng to it1500.  Residual depth summary (T2.3) shows
depths 66-70 carry ~95% of the residual column sum - inner-depth
dominated, consistent with a boundary-coupled slow mode.

purelam (no JFNK/Ng/dirchk, ZERO_FLUX): deterministic reproduction
of it170-178 demands (1180 -> 6.56), monotonic ~20%/iter
contraction, no limit cycle yet.  The decisive window is the
~1.6% plateau regime where dirchk_lc cycled; still ~15 iters away.

## 34. Trace-row MAXCH stall on m15 SOB l14 hydrot - RELFRAC fix (16-Oct-2026)

Run `cmfgen_validation/m15_2p82d_SOB_l14_hydrot` stalled at the
converged exit: MAXCH held at 0.3-0.9% (EPS_TERM=0.05) on dead
trace rows only.  SCRTEMP read at it381 confirmed every BA13-TRACE
offender is a dead stage at its depth (OI J=388 at d58 = 8.9e-39 of
the depth max; OSIX J=639-657 at d16 = 1e-23..1e-28).  Whole stages
are dead everywhere in this wind (OI max 2e-22, O2 6e-13, Sk I-III
<=5e-12, Ca I-III <=3e-13, Fe I/II <=2e-27).  Config fix only:
`BA_DEM_RELFRAC=1.0D-3` added to
`configs/dessart_m15_2p82d_SOB_l14_hydrot.yaml` (outer repo),
mirroring the luc_rfrac variant; corrections still commit on exempt
rows.  No source change.  Full details: convergence_issues.md
Sect. 34.

## 35. Stage-B trace polish made reachable - continuation routing (16-Oct-2026)

The converged-exit ELSE IF in `cmfgen_sub.f` armed BA_TRACE_PHASE
and then fell through to the CMFGEN_SUB return, so the Stage-B
polish was dead code.  New local `BA_STAGEB_ARM` (same tests as the
exit branch plus BA_TRACE_NEEDED) now excludes the exit and routes
through the ELSE continuation: BA_TRACE_PHASE is armed there,
BA_TR_NSWP incremented, and before `GOTO 20000` the sweep forces
COMPUTE_BA=T, LAMBDA_ITERATION=F, FIX_IMPURITY=F so BAMAT is
re-assembled at the held trace state (assembly is gated on
COMPUTE_BA .AND. .NOT. LAMBDA_ITERATION).  Default path unchanged.
Commit 545a3fe6; build clean.  Verification pending on queued job
12242816 (`dessart_m15_2p82d_SOB_l14_jfnk`, BA_TR_MAXSWP=3).

## 36. Convergence-exit hardening (16-Oct-2026)

Logic review of the convergence-control flow found four defects,
fixed in commit 23825fe1 (build clean, harness PASS):

- D: the converged-exit tests sat in the non-LST branch and required
  NUM_ITS_TO_DO .NE. 0, so the budgeted-final iteration (LST set when
  the counter hits 0) wrote out and STOPped without any convergence
  check.  BA_TRACE_NEEDED / BA_TRACE_DONE / BA_DEM_EXIT_OK and the
  new CONV_EXIT_MET are now evaluated at 9999 before the split; the
  LST branch reports whether the criteria were met, and the ELSE
  branch re-evaluates CONV_EXIT_MET after the MAXCH=100 hold-open
  ramps.  Iteration counts and write-out behaviour are unchanged.
- B: a one-time exit reason is now named in OUTGEN (CONV_DIAG_LVL>=1):
  clean / exempted demand saturated (>= BA_SAT_RAW) / Stage-B budget
  spent or still owed - the converged exit no longer looks identical
  whether trace demand was resolved or unresolvable.
- A: JFNK commits in JFNK_SOLVE_CORR_V1 already exported BA_TR_ABS
  but left BA_DEM_EXM=0, making exempted trace demand invisible to
  the exit hold on Newton iterations.  Both commit scans now also
  accumulate BA_DEM_EXM from the committed exempted change
  (best-effort proxy; the JFNK solve has no STEQ).  BA_DEM_EXM is
  initialized at VADAT read.
- C: verified the exit gates on demanded MAXCH in every mode by
  design - under DO_NONMONO_KNOP MAXCH_CTRL=MAXCH_APP is the applied
  step and would release on unsatisfied demand, so no code change;
  documented at the exit test and on the MAXCH_CTRL declaration.

## 37. Static-review bookkeeping fixes: JFNK relief reachability, floor sentinels (16-Oct-2026)

Static review of the 15/16-Oct hardening found the JFNK hysteresis
relief was unreachable and several flags/masks disagreed with their
docs.  Commits aab61f7d, 4a407476, 243cd3a1, 3b855cc9, e7d5d43b,
239b3262 (tip 239b3262 on development):

- The relief `JFNK_NREQ_EFF=1` keyed on `CNT_STABLE_LAM >=
  N_STABLE_LAM`, but the stable-exit release zeroes the count before
  the next eligible call.  New latch `LAM_STABLE_PROVEN`
  (CONTROL_VARIABLE_MOD) is set where the count earns the threshold
  and cleared on the instability/demotion resets; the relief keys on
  it and prints `JFNK: stable-LAMBDA relief active` on LUER.
- `BA_DEM_PFLOOR`/`BA_DEM_RELFRAC` defaults are now -1 = auto: under
  the non-mono flags they resolve to 1e-30/0.01, an explicit 0 is now
  honoured as off (previously rewritten), and the silent
  `BA_TR_MAXSWP=0 -> 5` coercion is removed so 0 disables Stage B.
- JFNK exempted-row scans in JFNK_SOLVE_CORR_V1 gained the
  `BA_DEM_PFLOOR > 0` guard so the mask matches SOLVEBA_V13 exactly.
- `MAXCH_TRACE_ABS` label/comments renamed to `BA_TR_ABS`.
- `JFNK_RAD_SWEEP_V1`: the alloc-fail tail now prints the diagnostic
  and poisons STEQ (the flag it toggled was only consumed by the
  skipped epilogue); the determinism bracket now also saves/restores
  MOM_GS_MAXCHG, MOM_RF_RESID_PRV, JNEXT_RELAX, INEXT_RELAX and
  KN_RF_RESID_PRV so eval order cannot perturb the production
  write-back weight.
- CONV_EXIT_WHY comment corrected: the "still owed" reason is only
  reachable on the LST_ITERATION path.
- README gains the VADAT table for BA_DEM_PFLOOR/RELFRAC,
  BA_TR_MAXSWP, BA_TR_ABSMAX, KN_ADAPT_RELAX, ALLOW_SIGMA_UNDEF,
  SOB_SIGMA_MAX.

./build.sh clean (ifort -O3, zero errors); all changed files
recompiled; `stable-LAMBDA relief` and `BA_TR_ABS` strings present in
exe/cmfgen_dev.exe.  Vault notes updated (JFNK admission, Convergence
criteria, VADAT keys).  Proof run to exercise the relief path: see
next entry / m15_2p82d_SOB_l14_jfnkrel.

## 38. 17-Oct-2026 -- jfnkrel validation runs: relief + coupled JFNK exercised (239b3262)

Two proof runs submitted against commit 239b3262 via
run_model_grid.py:

- m15_2p82d_SOB_l14_jfnkrel (job 12246975, obs 12246976): SOB clone of
  the run that deferred JFNK every call.  Outcome: the latch fires --
  `JFNK: stable-LAMBDA relief active (frozen-tier NREQ -> 1)` printed
  at it14/15 on the stable-exit-released call, and JFNK-FRZ then
  admitted repeatedly once demand fell below FRZ_GATE=3.0 (6
  admissions by ~it95).  Each solve short-circuits on the
  radiation-pinned gate (`residual unreachable by population step`)
  and falls back to BA -- correct: with J frozen the residual lives
  in rows a population step cannot move.  SE residual declining
  steadily; run healthy at write time.

- m15_2p82d_KNOP_jfnkrel (job 12247349, obs 12247350): KNOP coupled
  config (JFNK_COUPLED=T, ANYITER=T, COR_GATE=1e30, NREQ=1), warm
  restart from the m15_2p82d_SOB_l14_fixtest checkpoint (no KNOP
  I_NEXT state -> lagged store rebuilt over N_LAM_RWARM=4 + GS
  maturity gate).  GS field matured at it34-35; JFNK-CPL admitted at
  it35 and ran 30 coupled evals through JFNK_RAD_SWEEP_V1 (the new
  save/restore bracket): 18 JVPs, 1 LU factorisation, 0 NaN evals, 0
  failures.  GMRES stagnated (~7% residual drop, plateau 5.94e13);
  clean `gmres_stagnation_abort` + Lambda fallback.  HOWEVER: after
  the JFNK cycle the lagged field went runaway -- SE residual
  4.7e13 -> 2.2e31 by it40 and pinned there under GS_HOLD_POPS, FG
  churn [BAD] with capped writes.  The LAMBDA_GS rebuild did not
  recover it.  Interpretation: the coupled eval sequence (or the
  partial Newton step committed under near-null-row regularisation)
  left I_NEXT inconsistent with the moved field -- a solver/model
  behaviour finding, not a bookkeeping defect.  Watch for a
  watchdog/requeue; otherwise it burns to the wall limit.

Lesson: a clean JFNK abort restores the caller mode but NOT the
lagged radiation store the evals churned -- consider checkpointing /
restoring the I_NEXT store (or refusing admission) when the coupled
field is this marginal.

17-Oct-2026 (b7c2a38d): JFNK hardening pass (W1-W7) + tier-aware
reachability.  JFNK_FALLBACK_CAP is now the single source for the
LAMBDA->MAX_LAM / full->MAX_LIN cap mapping; every JFNK fallback is
capped by the caller's COR_LIM_IN (fixes the it35 +91% free-electron
step through a GS_LAM_CAP=3 window that ran the KNOP field to
se_res 2E31).  Reach gate is frozen-tier-only (stored BA is not the
coupled operator); reach skips count in JFNK_REACH_SKIP_STREAK and
the SOB stable-LAMBDA relief stands down at JFNK_REL_REACH_MAX=5.
JFNK_FD_COS_MIN=0.05 + JFNK_FD_AUTODEF=T now arm the FD descent gate
by default (stored-BA REACH 0.735 vs true-operator cos -0.23 on the
failing step).  LUM_BLOCKED revived via the LUM_DIV_* runaway
detector on |RLUMST(ND)| (KNOP only; demotes to damped GS_LAM_CAP
LAMBDA + holds pops; stops past LUM_DIV_CAP=1D15).  MAX_LAM/MAX_LIN
VADAT sanity checks; post-auto-set KNOP config re-print; distinct
catastrophic vs finite-stagnation messages.  Opt-in: GS_RAW_MATURE
(raw-max release criterion), JFNK_CPL_RELWAIT (release-window
coupled-tier veto).  Policy changes needing owner sign-off: coupled
tier no longer reach-gated, SOB relief stand-down, FD gate armed by
default.
Validation run: dessart_m15_2p82d_KNOP_jfnkrel_cap, CMFGEN job
12253518 (+ OBS 12253519), warm restart from SOB_l14_fixtest,
MAX_LAM=1E10 retained intentionally to exercise the new guards.

17-Oct-2026 (contd): KNOP FG_GS_MAXCHG denominator floor (commit
15797e32).  Live diagnosis of job 12255111 (b5b731ed binary): the
maturity gate was permanently pinned by the sweep-side max metric --
GS_SW_MAXCHG=1.56E263 from ABS(T1-T2)/T2 on a stored cell near the
1E-30 write floor, while every count-based gate passed (nbad frac
4E-4, nconv 0, commit rejects 84->3, div/heal/cap all zero).
GS_MAX_ALL exceeded GS_UNSAFE_TOL unconditionally, so
GS_FG_IMMATURE never cleared, FG_TRIAL could not arm, and
N_FROZEN was marching to the GS_FREEZE_LIM=60 non-convergent stop.
The NCONV counter already required MAX(T1,T2)>1E-4*VMAX_I but the
max metric did not -- floored at the ray column peak, matching the
MOM_GS_MAXCHG contract and the refresh-side MAX(VOLD,1E-4*VMAXF)
residual normalization.
Validation run: dessart_m15_2p82d_KNOP_w19 resubmitted as CMFGEN
job 12257428 (+ OBS 12257429), warm restart chained from w19's own
checkpoint (FGINEXT_GS carried, 16 immature-field iterations of
rebuilt I_NEXT preserved; GS_STATE preserved by the grid).

17-Oct-2026 (96df39ba): w19 runaway-containment audit fixes.  External
review of the w19 OUTGEN (L_out 9.5E9 -> 5.9E14 L_sun in 7 calls,
refresh residual 0.465 vs KN_RAY_TOL=2E-3, retained-field delta 21)
confirmed two code defects.  (1) KNOP_RAY_SOLVE_V1 committed the last
iterate of stall/cap/plateau ray exits (IXIT 3/4/6) to the lagged
store bounded only by FG_COMMIT_JUMP=50 -- an unconverged fixed point
could enter the store and poison the next refresh.  New VADAT key
FG_COMMIT_CONV (default F) skips WRITE_STORE on non-converged exits;
rejections print as 'nc rejects' (FG_NC_REJ).  (2) LUM_BLOCKED was
recomputed per call, so one sub-threshold iteration cleared the hold
at a still-runaway luminosity, and LUM_DIV_FAC=1E2 could not see a
~10x-per-call geometric runaway.  LUM_BLOCKED is now a latch released
only when L_out <= LUM_DIV_REC * pre-runaway level (default REC=1.0;
0 = never release); slow-streak tier LUM_DIV_FAC2=5 / LUM_DIV_N2=4
catches the geometric growth; LUM_DIV_CAP stop unchanged.  All new
keys default to the old behavior; build.sh clean.
Validation run: dessart_m15_2p82d_KNOP_w19 resubmitted with
FG_COMMIT_CONV=T (see config vadat_overrides).

17-Oct-2026 (d4d45893): KNOP acceptance-state audit -- residual-gated
  commits, capped provisional writes, bulk contraction.  External
  review of 96df39ba on origin/development confirmed the bulk exit
  (IXIT=2) still reached WRITE_STORE unconverged: it fires on the
  bad-cell fraction alone, and the w19 poison (residual 0.25-0.49,
  retained deltas 8-27x) committed through it.  Also confirmed the
  resubmitted w19 job 12260766 was built from 96df39ba only -- the
  FG_COMMIT_NC/FG_BULK_CONTR yaml keys were added post-submission and
  are absent from its VADAT.  This commit closes the gap:
  FG_COMMIT_CONV now tests MAX(RES_M,RES_P)>KN_RAY_TOL; FG_COMMIT_NC
  (>0) gives a capped provisional commit (WRITE_STORE_CAP, per-cell
  change <= NC*retained) instead of a skip; FG_BULK_CONTR requires
  same-direction contraction before the bulk exit may fire in either
  sweep; per-ray verdict census (committed/jump-rej/nc-rej/capped/
  div/unsolved) plus worst-cell persistence telemetry; AA rollback
  now recomputes the NBAD/NTOT census on the restored iterate
  (conditional, KN_ORD_AA path only).  All defaults preserve the old
  behavior; build.sh clean, unit harness PASS.
Validation run: dessart_m15_2p82d_KNOP_w19 resubmitted (kills
  12260766, which was running the weaker gate and healing -- se
  2.0E14 -> 8.0E13 over it81-86 with pops held, n_frozen=14/60).
17-Oct-2026 (fix): w19 job 12261599 died 11 min in on a format
  overrun in the new verdict-census WRITE ('(A,7(A,I5))' vs 8 item
  groups; forrtl severe 61 at fg_j_cmf_nonmono_v1.f:6062).  Fixed and
  rebuilt clean; resubmitted.
17-Oct-2026 (b398571a): BA_DEM_RELFRAC removed entirely.  The SOB
  l14_jfnk run (12258403, 242 its, budget-exhausted) showed why: the
  relative floor exempted ~700 rows/depth at >0.1% while MAXCH read
  ~0.08% and the converged-exit verdict fired -- an unconverged
  state certified as converged (user: "leaves the system
  unconverged and the output spectrum incorrect").  Exempted demand
  had asymptoted to a ~0.33% floor (FeII a6De[9/2] d4), never
  converging.  Removal touches CONTROL_VARIABLE_MOD, rd_control_variables
  (key + trace auto-set), SOLVEBA_V13 (5 sites), JFNK_SOLVE_CORR_V1
  (2 sites), CMFGEN_SUB and SOLVE_FOR_POPS arming conditions -- all
  now PFLOOR-only.  BA_DEM_PFLOOR (absolute 1E-30) retained; it only
  exempts physically-empty rows.  Outer repo configs stripped
  (commit 726fd106, 24 files) -- the key is now unrecognized and
  would warn.  Consequence: MAXCH on future SOB/KNOP runs will read
  the real trace demand (~0.3%+); converged exits now require the
  trace rows to actually settle.  Build clean, harness PASS.
18-Oct-2026 (fix): verdict-census crash root cause found -- on this
  ifort a widthless A followed by a parenthesized repeat group
  ('(A,N(A,I5))') raises forrtl 61 at runtime; both 12261599 and
  12266023 died on the same line.  Format now leads with the group
  '(6(A,I5),A,I5,A,I5,A)'.  Verified by standalone reproducer.

## 18-Oct-2026 - KNOP w19 (job 12267898): luminosity guard caught the seed-fill escape

The format-fixed run cleared the verdict census and completed two refreshes,
then failed at it90: L_out ran 0 -> 7.4E15 and the luminosity latch stopped
the run. Telemetry showed the residual-aware gates working (jump-rej=50,
nc gated), but a converged-classified ray seeded never-committed frequency
columns with unbounded values - FG_JUMP_MAX skipped any column whose store
peak sat at the 1E-30 write floor (the "first real write" cold-start
exemption). Committed delta reached ~3E7 while JUMP stayed <50.

Fix (commit bad4e470): fills into floor-level columns are now bounded by
FG_COMMIT_JUMP times the sibling-ray peak at the same field frequency;
only a column at the floor for EVERY ray keeps the cold-start exemption.
Telemetry now reports the max delta among committed rays separately
(committed-max=) so committed vs rejected deltas are distinguishable.

## 18-Oct-2026 - KN_NO_GAIN runtime-verified; BA_EX_RELFRAC exit-only floor added

Diagnostic job 12270690 (interactive, 8 cpu, warm from the poisoned
it99 FGINEXT_GS) verified KN_NO_GAIN: census holds at 5139 cells,
depths 1-59, max chi_scat/chi=15 at lambda 3.6E5 A; L_out decays
~1.5-2.2x per iteration (7.3E13 -> 1.6E12 over its 99-105) instead of
compounding 8.6x/iter.  The lagged-map photon-gain channel is closed.

Major/trace audit: the PFLOOR=1E-30 predicate is identical at all six
sites in both paths (SOLVEBA_V13 + JFNK).  Inverting rows are majors
(n_up 7E-7..2E3 cm-3, all >> floor) - no hidden-row leak.  Gaps found:
SRCE_CHK first branch has no CHI_SCAT bound (w19 switched to
ESEC_CHK); Stage-B is not a true lower-triangular solve; no hysteresis
on the classification; the 1e-6-relative demand class has no
mechanism to converge after RELFRAC removal.

Commits: dde0dcd3 (KN_NO_GAIN comment fix), 39b94dff (BA_EX_RELFRAC).
BA_EX_RELFRAC (VADAT, -1=off) exempts rows below relfrac*depth-max
from the EXIT demand BA_EX_MAXCH only; MAXCH, histogram, committed
scan and level-check keep every row.  BA_DEM_EXR reports the hidden
max every call and gates the exit through BA_DEM_EXIT_OK (small or
saturated).  JFNK commits mirror EXR on POPS_PREV.  Config: w19 +
l14 set BA_EX_RELFRAC=1E-6; w19 NEG_OPAC_OPT -> ESEC_CHK.

## 18-Oct-2026 - KN_NOSLOT sink census; monotonic A/B submitted

Commits: 5f0932d8 (sign-error comment fix in FG_J_CMF_NONMONO_V1),
2624d15d (KN_NOSLOT census).  A compressing ray-point (GAM<0) with no
bracketing compressing depth keeps chi_hat = chi + |GAM|/DLN but drops
|Q|*I_NEXT - absorb-only, a non-conservative sink.  Now counted per-ray
and reported once per great iteration; diagnostic only.

Validation: full ./build.sh clean (contended node made it ~2h);
tests/jfnk_core all PASS; KN_NOSLOT symbol in cmfgen_dev.exe.
Monotonic A/B submitted on a velocity-reversed SN_HYDRO_DATA (sigma
> -1 everywhere, zero compressing cells), 6 its, LTE_EST=T, all other
keys identical: MONO_KNOP job 12273306 vs MONO_V13 job 12273308.
Pass criterion: identical L_out / RLUMST / J to round-off, since the
KNOP expanding branch must reduce to the V13 sequence.

Live run: w19_lte (12271873, binary 39b94dff, no census) at it15 still
held (mature 0/2, bad frac ~1.5E-3 vs 1E-3 gate, GS_FREEZE_LIM=60).

## 19-Oct-2026 - KNOP maturity-gate repair + band-edge coupling

Commits df746176, af91495d, 27b9a4f5, f7ace371, 44c07359, 90f1104e,
dc41573a, cb19f307, 2696acfc (HEAD=2696acfc, pushed to
origin/development).

Correctness (group A): DO_KNOP_REFRESH is now read after USE_LAM_ES
is final and MOM_RF_ENABLED is re-synced on the post-auto-set value;
NOV_CONT=T is forced off under KNOP (GLOBAL_LINE!=BLANK would clear
CONT_VEL in COMP_J_BLANK and silently bypass the Knop continuum);
KN_ORD_P2DAMP default made FALSE in both rd_control_variables and
the FG module; KNOP knob echoes clamp before LOG10 and print
post-clamp values.

Gate bookkeeping (group B): sweep cap events, refresh cap events
and cap denominators are now separate censuses
(FG_SW_CAPPED/FG_SW_CAPTOT/FG_RF_CAPPED); SOLVE_FOR_POPS uses
same-census cap shares and falls back to a bulk criterion when the
refresh did not run (GS_NTOT=0) instead of letting a single
pathological cell pin GS_UNSAFE.

Stabilization (group C): DO_NONMONO_FREQ_MATCH auto-sets T under
KNOP when the VADAT key is absent (explicit F respected); the
bluest frequency now couples compressing cells to the lagged redder
slot instead of forcing Q=0, in both the sweep and BUILD_FREQ (the
matched NM_RAY/DLN_RAY stencil was hoisted above the INIT split).
The 3a*I nu-form cancellation and the red-end phantom-slot closure
were audited and documented in the source.

Build/tests: ./build.sh clean incl. the full unit_tests suite (now
run inside build.sh); 49/49 VADAT key tests PASS.  MONO_KNOP
12277527 (fix_temperature pinned) is queued; it will start on the
2696acfc binary and is the first real-model gate read on the new
code.

## 19-Oct-2026 - JFNK cache, PTC and globalization regressions

Commits: acfc885b (preconditioner cache reuse), 0960b767 (PTC/log
scaling, repaired-residual rejection and temperature-floor handling),
ffd698ff (restart trust bound, maxch convergence and GS_FREEZE gate).

`JFNK_FACTOR_PRECOND_V1` keeps cached BA diagonals valid and no longer
holds banded `REPLACE` state in a local allocatable. Active PTC factors
rebuild with the current `dtau`; log rows fold `sigma*exp(-x)` in
fractional units. JVPs abort on repaired perturbed residuals, and the
existing 1 K temperature floor is counted and included in the FD step
cap. Restarted GMRES retains earlier cycle contributions and bounds
the returned total step. Coupled solves require the residual tolerance
even when maxch passes. `GS_FREEZE` is asserted false outside the
`DO_NONMONO_KNOP` gate.

Validation: `./build.sh` passed, including the 13-file unit suite and
the `tests/f90` solver/interface tests. `tests/jfnk_core` passed T1-T3,
G1-G8; focused regressions cover PTC/log scaling, GMRES restart bounds,
failed JVP signaling, coupled maxch and source guards for cache reuse,
the 1 K floor and the KNOP freeze gate. No model run was submitted.
Adaptive multi-step FD selection and coupled evaluation repeat-check
remain open; complex-step remains unsupported and fails closed.

19-Oct-2026 — TRI hardening merged; JFNK+TRI probe submitted.
Commits 7166919b..db38161f (branch tridiag-hardening-t1-t4, merged
no-ff) plus 69c54e9e (INT-cast fix for the FG_NOSLOT_RAY sum; the
SUM(...,KIND=8) form did not compile under ifort).  TRI failures now
degrade: recall failure falls through to compute, compute failure
clears partial BCD cache and retries via the DIAG block solve
(CMF_DIAG_SOLVE); recall gates INQUIRE every depth; NUM_BNDS<3 forces
DIAG at read time; LAST_MATRIX flags the K=ND diagonal band build.

Validation: `make DEBUG=nodebug all` exit 0; tests/f90 suite PASS
(T2 Rosenbrock now converges in 8 iterates); unit_tests 13/13 files
pass.  `exe/cmfgen_dev.exe` at 69c54e9e.

Run: cmfgen_validation/m15_2p82d_SOB_l14_freet_jfnk_tri — warm restart
of the converged freet model (it=850 checkpoint) with DO_JFNK=T,
JFNK_USE_BANDED=T, SOL_METH=TRI, NUM_BNDS=3.  Jobs: CMFGEN 12285271,
OBS 12285272.  Note: JFNK ignores the FIX_X_DPTH T-pin while admitted
(mirrors MOD_FIX_T only) — solver robustness test, not a physics
variant; do not compare its inner T structure against freet.

19-Oct-2026 — resubmitted freet_jfnk_tri with DO_FG13_ALO=T.
Cancelled 12285271/12285272 (still pending) and resubmitted with the
FG_J_CMF_V13 diagonal ALO added: CMFGEN 12285416, OBS 12285417.
Config otherwise unchanged (DO_JFNK=T, JFNK_USE_BANDED=T, SOL_METH=TRI,
NUM_BNDS=3, warm restart from the converged freet checkpoint).

19-Oct-2026 — KNOP production candidate submitted: KN_l14_freet.
m15_2p82d_KNOP_l14_freet warm-restarts the converged freet checkpoint
(it=850) onto FG_J_CMF_NONMONO_V1 at commit 69c54e9e.  Physics deck
mirrors freet (thin-region T release FIX_X_DPTH=21, EXPO clumping,
DIFFUSION boundary, dielectronics, THK_CONT=T); guard stack from
w19_v2 (FG_COMMIT_*, KN_NO_GAIN, BA_RES_ACC, GS maturity gates,
N_LAM_RWARM=4).  JFNK admission retained (DO_JFNK=T, coupled tier).
Jobs: CMFGEN 12285736, OBS 12285737.  Prior real-field KNOP runs all
failed pre-guard-fix (w19 lum-divergence, w19_lte/v2 GS_FREEZE_LIM);
this is the first KNOP attempt on the real field with the bluest-
frequency finite guards in place.

19-Oct-2026 — KNOP field-drift audit fixes on knop-field-drift-fixes.
Static audit of the drifting KNOP_l14_freet run (job 12285736, refresh
residual rising 0.06->0.15 vs tol 2e-3 with pops held) found the map
algebra sound but eight defects keeping the field from a fixed point.
Committed individually on the branch (68bba35c..2bdedbff):

- FG_COMMIT_CONV/FG_COMMIT_NC/FG_BULK_CONTR defaults flipped to
  T/0.05/T — fail-closed commits were opt-in, so every generation
  committed unconverged stall/cap/plateau iterates by default.
- HEAL_STORE zeroes fully-poisoned depth rows (the NV.EQ.0 skip left
  NaN rows re-seeding at residual 1.0).
- In-sweep I_NEXT store write + churn metrics gated on NEW_FREQ —
  comp_j_blank re-enters the solver per J accuracy iteration, so the
  store relaxed ~1-(1-omega)^k per frequency and double-counted the
  maturity signal.
- RESTORE_FG_INEXT_LATEST runs the INIT-path value audit before
  accepting a rollback restore.
- FG_TRIAL excluded from GS_HOLD_POPS — held trial probes were
  discarded while consuming release budget.
- GS_RAW_RELEASES_TOT (integer*8) counts raw-cure frozen->live
  transitions, published as GS_STATE n_raw_releases=; GS_REL_N
  default 1->2.
- One-time OUTGEN note that the frequency-1 field has no lagged
  store slot (band edges carry zero memory; slot convention kept).
- Static checks: unit_tests/python/test_knop_lagged_field.py
  (6 checks, pytest -q passes, exit 0).

No compile per audit scope; runtime validation pending on the next
KNOP submission.  Consumer grep: FG_COMMIT_*/FG_BULK_CONTR only in
rd_control_variables + fg module; GS_REL_N only in the release/JFNK
damp windows; FG_TRIAL only in solve_for_pops — no unexpected sites.

19-Oct-2026 — KNOP sweep-march fixes on knop-gs-sweep-fixes.
Follow-up audit of the drift path (branch merged into development as
5a37d0c2, pushed): MARCH accumulated ABS(DCELL) into SGN_ACC so the
W6 sign-flip/period-2 damping in SWEEP_FLD could never fire — the
unordered sweep (KN_RAY_ORDERED=F, the production path on m15) ran
undamped while SWEEP_FLD_ORD's signed KN_SGN_SUM worked; SWEEP_FLD_ORD
boundary rows (I=1, I=NI) reprocessed out-of-window frequencies every
hot pass, polluting the pass census; march cap hits were folded into
NNEG instead of NCAP.  Commits 5b6b0488/d6627158/d98f1219 + tests
7eef81c6.  MARCH caller scan: single CALL MARCH inside SWEEP_FLD
(line ~6966), no other call sites — the signed accumulator cannot
leak into the in-sweep production store write.

21-Oct-2026 — SOB/moment-solver OUTGEN status text on sob-status-text.
Sentence-based narration for the lagged-field runtime output, same
style as the KNOP gate lines, shared by both drivers through the new
module new_main/mod_subs/knop_status_mod.f90 (KNOP_STATUS_MOD):
mom_state/sweep/refresh/pass/restart/direction sentence helpers.
solve_for_pops KNOP GS sweep report and MOM rows,
mom_j_cmf_nonmono_v1 GS sweep table, refresh pass table + exits,
direction-stall drops, guard line, restart/restore decisions,
accelerator and adaptive-relaxation summaries now print prose;
directions are words (compressing/expanding), booleans are yes/no,
and every number plus the [OK]/[WARN]/[BAD]/[HELD]/[SAT] tags are
preserved.  Machine records untouched: GS_STATE keys
(mom_refresh_conv=), ITER_SUM/ITER_HISTORY rows, MOMJNXT_GS
MAGIC/VER format, VADAT echo, audit-parsed strings.  'Current
great iteration count is' moved to (A,I0) so the audit regex
matches at count >= 1000 (was (A,I4), which glued).  Commits
55266e52..5c82f668, merged --no-ff as fc68c9a9 and pushed.
Static checks: unit_tests/python/test_knop_status_text.py
(7 checks; full suite 23 passed).  Syntax-compiled each edited
Fortran file against mod/ with ifort; no model run per audit
scope — runtime log review pending on the next SOB/KNOP
submission.

21-Oct-2026 — KNOP precondition stall fixes on development (runtime
audit of job 12298049, m15_2p82d_KNOP_np cold KNOP start: demand
pinned at 1E15..1E72 by degenerate rows, N_PRECOND_TOT heading to
the 50-call stop).  Three surgical corrections: (1)
fiddle_pop_corrections_v2.f MAJOR scale -- a row whose demanded
correction exceeds the per-element apply-clamp bound can no longer
veto the whole depth's update (nonmono-gated, stock path unchanged);
(2) solveba_v13.f -- DO_BA_QUAR extended to the LAMBDA phase: rows
demanding beyond BA_LQ_DEM% (default 1E12) for BA_LQ_PAT (default 8)
consecutive LAMBDA calls park into BA_QUAR_MASK and their STEQ is
suppressed on every call; new VADAT keys BA_LQ_DEM/BA_LQ_PAT read in
rd_control_variables.f, hit counter BA_LQ_HIT in
control_variable_mod.f; (3) solve_for_pops.f -- N_PRECOND_TOT counts
only non-improving calls vs the GS_PC_BEST tracker (10% margin, same
as the release audit).  Build clean; unit suite 16/16 after
widening the test_knop_lagged_field regex window for the moved
LAMBDA_ITERATION line.

22-Oct-2026 — external source review of the DO_NONMONO_SOB pathway:
four localized fixes on development.  (1) comp_j_blank.f -- the
FIRST_FREQ one-sided dLOG_NU fix was DO_NONMONO_KNOP-gated only;
widened to DO_NONMONO_SOB+USE_LAM_ES=F since MOM_J_CMF_NONMONO_V1 is
now the SOB default field solver and consumes dLOG_NU at the bluest
frequency (dLOG_NU is RADIATION_MOD state, stale there otherwise).
(2) solve_for_pops.f BA_SAT_CNT -- excluded JFNK_THIS calls: a JFNK
step returns the applied correction in MAXCH, not raw demand, so it
must not feed the saturated-demand demote counter.  (3)
LAM_STABLE_PROVEN now cleared at the cnt_lam_escape and
rgrid_revised sites -- a cadence escape is not a proven stable exit,
and a regrid invalidates the proof; a stale latch could arm the SOB
stable-LAMBDA JFNK relief (JFNK_NREQ_EFF -> 1) early.  (4) the
held-gate probe now requires BA_PROBE_LAM > 0: MAX(BA_PROBE_LAM,
RD_CNT_LAM) had let the "0 = off" key probe anyway on the RD_CNT_LAM
cadence.  Review's item 5 (SOB MAXCH_CTRL uses raw MAXCH, not
MAXCH_APP as KNOP does) left unchanged -- deliberate per the
comment block; flagged as the main SOB efficiency drag for a later
A/B.  Build clean, unit suite 16/16.

22-Oct-2026 — second tranche of review fixes on development, one
commit per item (67299053..63a13a44).  (2) rd_control_variables.f --
NG_ABSTAIN base default moved ahead of the enforcement block and the
read site no longer re-defaults it, so the SOB absent-key auto-set
survives; GET_KEY_STRING cannot re-probe a key right after
RD_STORE_LOG consumes it (I_DWN/I_UP cursor brackets the index), so
the presence probe must run at the enforcement block.  VADAT tests
added: absent -> T + logged auto-set; explicit F respected.  (5)
GS_FREEZE promoted to CONTROL_VARIABLE_MOD (re-derived every
SOLVE_FOR_POPS call, reset at entry); freeze-trial LAMBDA calls no
longer increment CNT_STABLE_LAM in either stability counter -- a
damped held-pops probe is neutral, not success/failure evidence --
and the DO_LAM_AUTO handoff in cmfgen_sub.f now also refuses to fire
while GS_FREEZE (KNOP-only flag, forced F outside).  (6)
check_param_consistency.f -- the SN_MODEL check requiring
USE_J_REL/USE_DJDT_RTE/USE_LAM_ES is skipped under DO_NONMONO_*;
explicit USE_LAM_ES=F is the honoured/announced pairing there.
(7) rd_control_variables.f -- BA_TR_ABSMAX auto-set trigger moved
from .EQ. 0.0 to .LT. 0.0 so explicit 0 = off survives (the .EQ. 0.0
form could only ever fire on an explicit user 0 and rewrote it);
solve_for_pops.f CONV_TRC_NEED and the cmfgen_sub.f BA_TRACE_NEEDED
replica now require BA_TR_ABSMAX > 0 so 0 disables the absolute
trace test.  (9) do_ng_accel_band_v4.f error strings corrected
V2 -> V4.  (10) rd_control_variables.f -- absent-key DO_BA_RES_ACC
auto-set to T under DO_NONMONO_KNOP ('KNOP auto-set:' log line);
every validated KNOP config already sets it explicitly, and
SOLVEBA_V13 covers the LAMBDA-phase calls plus the MOM/FG
lagged-field rollback, so baselines are unchanged.  Explicit F is
still honoured for controlled F baselines and keeps warning.
Build clean; unit suite 16/16 (57/57 VADAT).

21-Oct-2026 - JFNK deferral-path operator probe on development.
Motivation: m15_2p82d_SOB_l14_lte_np (~it600-710) admits JFNK-FRZ
every call and defers at the FD-cos gate; census shows |F|^2 entirely
in the top 1% of rows (HI lev 22/23, depths 66-70), where the stored
BA prices the clipped step at ~+9e9 per row while the FD JVP returns
~-1e6 - a per-row operator disagreement that global
cos(BAz,A_FDz)=1.0 masks.  The existing top-|F| column probe and
h-scan sit AFTER the defer, so repeated deferrals never produce the
evidence needed to separate stale-BA from FD-cancellation from a
truly unreachable direction.  Change is diagnostic-only: new VADAT
key JFNK_DEFER_PROBE (default 0 = off) runs the probe on the first N
FD-gate deferrals - the clipped-direction JVP at h x{1,1e2,1e4,1e6}
(scale 1 reuses the gate JVP), the unclipped direction once, and the
top-|F| Jacobian column at the same scales - logging per-row F_i,
BAz_i (clip + raw), FD(h), raw dF=F(x+hz)-F(x), actual ||hv||inf
(new JFNK_PROBE_H_LAST export from JFNK_JVP_V1), and BA source-iter
age (new JFNK_BA_SRC_IT stamped at 'BA matrix computed').  Per-row
verdict labels: h-unstable (cancellation) / h-dependent (nonlinear)
/ stable FD ~= BA / stable FD != BA (stale stored row).  The defer
decision, gates, thresholds and fallback are unchanged; no commit,
no field change, JFNK_PROBE_HSCALE reset after every eval, evals
respect NEV_CAP headroom.  Files: jfnk_solve_corr_v1.f,
jfnk_cmfgen_int_mod.f, jfnk_cmfgen_interface_v1.f,
control_variable_mod.f, rd_control_variables.f, solve_for_pops.f.

21-Oct-2026 (2): DEFER-PROBE runtime fix + first live read
(a18708c6). Bounded restart of the converged np model with
JFNK_FD_COS_MIN=1.0 forced real FD-gate deferrals; two probes ran
end-to-end in /ptmp/sbrennan/jfnk_deferprobe_test.  Fixed a format
mismatch that killed probe 1 mid-write ('(A,5(A,ES12.4))' pairs the
lead real arg with A, shifting strings onto ES -> forrtl(61); now
'(A,ES12.4,4(A,ES12.4))') and the src<0 BA-age print (disk-loaded BA
has no in-process stamp; prints 'src it=unknown').  Result: stored BA
and FD agree in magnitude AND direction on the Newton step
(|BAz|/|F| = |Az|/|F| = 5.12e10 over 6 decades of h, column diag
-1.18e19 vs -1.12e19, all 8 top rows 'FD = BA within 2x').  The
frozen-J deferral is not stale-BA: the Krylov direction itself is
unusable (|Az| ~ 5e10*|F|, cos(F,Az) = -0.21).  Gate veto correct.

21-Oct-2026 (3): OPTIONAL-arg segfault fix + full frozen-J experiment
(50307aa7).  A gates-off interactive run (npq_dbg, REACH/FD_COS/
FD_AUTODEF=0) let the real LU+GMRES solve run for the first time and
SIGSEGV'd in JFNK_FACTOR_PRECOND_V1's PRESENT(DO_FACTOR) -- five
no-arg calls through the implicit EXTERNAL interface left the
optional-arg protocol undefined.  All sites now pass .TRUE.
explicitly.  After the fix the solve ran three times: GMRES
stagnates at ~15% residual cut, the block-LU preconditioner is
~8e8 off on v=F (||Mz-v||/||v||), the inexact direction runs away
(|s|_w ~1e5-4e8) and every trust-region step is line-search
rejected.  Frozen-J Newton is confirmed unusable on this system -
not a gate artifact.  JVP diag caveat: on the GMRES basis direction
|FD|/|BA v| = 36, so BA fidelity is direction-dependent.

21-Oct-2026 (4): Review corrections to the full-solve interpretation
plus the JVP-diag extension.  (a) The PC audit's ||Mz-v||/||v||=8e8 on
v=F measures the POST-clip/projection z inside JFNK_PRECOND_V1, not raw
LU error -- v=e1 round-trips at 1e-4, so the factorization is fine and
the scalar mostly prices clip+projection loss on the direction F
produces (33026 ZMAX hits).  (b) The 36x JVP-diagnostic gap is a real
inf-norm worst-component result on the same masked v/h/state with PTC
on both sides; the residual is STEQ_coll + stored rad + CMAT_RAD*x and
BA_A is the analytic SE block plus PTC fold, so the gap reduces to how
faithful GENERATE_FULL_MATRIX_V3's analytic dSTEQ_coll is on the stiff
rows (vs a stale-matrix or missing-response question).  (c)
JFNK_ENRICH -- the anchor-column/capacitance machinery that repairs
exactly this stored-response gap -- is gated on JFNK_COUPLED_LOC and
runs plain LU in frozen-J mode.  Extended JVP_DIAG (one-shot) now
prints the top-8 |FD-BA*v| rows with K/J, species label, |F_i| and
row-max |BA|, then re-evaluates the SAME direction at h*0.1 and h*10
on the same snapshot: gap shrinking ~linearly with h = truncation,
h-invariant = structural analytic-BA misspecification.  Gates restored
on the npq_dbg yaml for routine runs; the gates-off config stays
debug-only.  User format fix in mom_j_cmf_nonmono_v1.f (trailing 'A'
descriptor for the 10th arg ' cells with sigT>1E4', previously hitting
format reversion) verified and committed.

21-Oct-2026 (5): JVP-diag h-scan first firing (852efafd, npq_dbg
gates-off interactive).  On the first post-LU GMRES basis vector:
|FD|=5.01e13, |BA*v|=1.56e12 (worst-component ratio 32).  All 8 worst
|FD-BA*v| components are the same stiff inner-boundary cluster (HI
21-23 at K65-70, one He2) -- a localized cluster, not a single-row
artifact.  h-scan on the SAME direction: h*0.1 -> |FD-BA*v|=5.97e14
(~1/h growth, NOT truncation), h*10 -> 6.19e12 (ratio ~4x residual).
A fixed-size jump model fits: FD_i ~ a + J0/(h*v_i) with post-jump
slope a~1.4e12 and jump J0~5e13*h0*v_i -- consistent with a
state-dependent switch inside the perturbed SE eval (dissolution /
REPLACE-class branch flipping on any nonzero perturbation along v)
rather than smooth curvature; the JVP clip-retry did NOT fire (no
new-constraint-crossing message), so the nonsmoothness is inside
JFNK_RESIDUAL_V1's own eval, not the monitored x-bounds.  Verdict:
the 36x is substantially an FD-eval artifact on this direction, not
proof BA is wrong -- and GMRES stagnation is consistent with a
contaminated (nonlinear, switch-polluted) secant as the Krylov
operator.  Frozen-J is partially rehabilitated: fixing the eval-side
nonsmoothness is now the prior task, before concluding on the
formulation.  Per-row FD_h print added (852efafd) to pin which row
carries the jump on the next firing.

21-Oct-2026 (6): switch hypothesis refuted; cancellation floor proven.
The three-mechanism probe (b65dd982) fired on npq_dbg: baseline
re-eval max|F(x0)-fx| = 0.0 exactly (no stale reference); repeat eval
diff = 1.7e-18 (deterministic); +/- direction test nearly
antisymmetric (mid << amplitude) with sign flips across h -- not a
one-sided switch.  Raw dF_i at the worst stiff rows is pinned near
~1e4 at small h regardless of step: the SE rows cancel ~1e19
inflow/outflow terms to ~1e12 residuals, so the deterministic eval
floor (~eps*1e19*sqrt(Nterms) ~ 1e4) exceeds the true linearized
signal h*Jv ~ 20 at h ~ sqrt(eps).  FD on those rows is noise/h --
meaningless at every small h.

21-Oct-2026 (7): resolved-h window found (7154f48b sweep).  Log sweep
on the live GMRES direction: |FD-BA*v|/|BA*v| = 2400 at heff=3e-12,
dropping to 0.044/0.018/0.018 at heff=3e-7/3e-6/3e-5 with raw dF
growing exactly 10x/decade (linear response, no truncation).  BA*v is
accurate to ~2% on the stiff rows once the FD resolves the signal.
Window in scalar h ~ 2e-7..2e-4; the noise floor |dF|~1e4 is measured
directly, confirming eps*term_scale.

21-Oct-2026 (8): JFNK_EPS honored in frozen mode (1d4ded6a, opt-in,
default unchanged); JFNK_CLIP_RETRY added (316e1e91, default 2) plus
fix: exhausted-retry path used to shrink H_STEP AFTER the last eval
then divide by it, inflating the kept JVP 10x -- part of the reported
FD-vs-BA mismatch on every prior admission was this artifact.
Directional clip cap (f8bbd7bd): the +-40/50 log-clip cap previously
used |v_i| symmetric, so a bound-saturated unknown moving back into
range strangled h to ~2e-9; now caps only outward motion (v_i*x_i>0),
matching the T-floor rule, and the binding idx is logged.  GMRES still
stagnated in all three experiment runs (12311079/12311336/12311542)
because the effective h never reached the window: vn of the Krylov
direction is ~2e3, so JFNK_EPS=1e-7 gives h ~1.9e-9 scalar.  Next:
JFNK_EPS=1e-4 puts scalar h ~2e-6 mid-window (12311716).

22-Oct-2026 (1): multi-direction JVP audit (JFNK_JVP_NDIR=4, cheap
norm line per direction + one-shot extended sweep).  Run 12311936:
|FD-BA*v|/|BA*v| = 1.77/1.78/1.67/1.28 % on the first four Krylov
directions -- BA is a uniformly accurate Jacobian on the explored
subspace, not just one direction.  GMRES still stagnated (~0.4% cut),
so the obstruction moved off the matvec entirely.

22-Oct-2026 (2): apply-side audit narrowed the failure.  PC audit
v=F gave ||Mz-v||/||v||=3.2E10.  Adding the reg-row split (257b1900)
and factored-M check (198ed8c3 precursor commit) showed the 6.6E22
residual lives on NON-reg pop rows even when M is reconstructed as
the factored matrix (BA+ridge, identity on reg rows): a genuine apply
failure, not an audit-vs-substitution artifact.  Weighted |F|^2 share
on non-reg rows = 1.0 -- the Krylov residual is entirely on rows the
preconditioned operator can reach, so unreachable-reg-row stagnation
is ruled out.  40079 rows (~31% of N, ~60%/block) are reg-substituted,
all on JFNK_POP_FL=1e-6.  Next probe: per-direction M^-1 J v
round-trip (198ed8c3) feeding JV and JV_ANAL back through the apply;
reports |M^-1Jv-v|/|v| on reachable rows plus reg share, with
pc_clip_event and ZMAX/reg counters saved/restored so the extra
applies cannot abort the live GMRES.  Submitted as 12312771.

22-Oct-2026 (3): Krylov machinery verified; FD shown inconsistent across
the basis.  Added JFNK_GMRES_STAG VADAT key (e49cad20, default 1D-3;
0 disables the 3-stall abort).  Run 12314393 (JACOBI, STAG=0,
MAXGMRES=300): RESTART DIAG printed hess_res vs true_res = 3.3e10 vs
5.1e12 at it=30 (156x) and 3.4e10 vs 1.5e13 at it=60 (455x) -- the
Arnoldi least-squares residual diverges from the true residual by two
orders of magnitude.  Root cause: the FD JVP picks an adaptive h (and
clip retries) per call, so the implied operator is mildly nonlinear in
v and the Krylov basis does not compose linearly.  Per-direction
FD-vs-BA agreement (~1.5%) never bounded this accumulated effect.
Also found: restart cap is hardcoded (n_restart>1 -> exit), so
MAXGMRES beyond ~2 cycles is moot.

22-Oct-2026 (4): analytic matvec experiment (f0734867, JFNK_BA_JVP,
default F).  Run 12314698 (JACOBI + BA diag-band matvec):
hess_res == true_res to 4 digits at both restarts (ratio 1.0000) --
with an exactly linear operator the GMRES machinery is verified
end-to-end.  But res still plateaus at 3.78e10 vs 4.13e10 initial
(~92% unreachable, flat over 60 its).  Conclusion: under the
diagonal-band operator, the weighted residual mass is genuinely
outside range(J_diag).  Open: whether the dropped cross-depth bands
(BA_ADV_TERM advection + ALO dJ/dopacity; only in SE(ID)%BA
off-diagonal bands, requires JFNK_USE_BANDED=T + NUM_BNDS=3/TRI)
carry the missing coupling.

22-Oct-2026 (5): banded-JVP run 12314940 (JACOBI + JFNK_USE_BANDED=T +
NUM_BNDS=3/SOL_METH=TRI, forced by prepare_VADAT) crashed with SIGBUS
21 min in, right after the admission census.  Root cause:
JFNK_FACTOR_PRECOND_V1 ran the full banded-LU assembly on
JFNK_USE_BANDED_LOC alone -- under JACOBI it allocated
JFNK_BANDED_LUAB (LDAB~10.9k x N~127k ~ 11 GB) plus the ARAW copy
(~11 GB) on top of CMAT_RAD_BANDS (~5.5 GB) and BA_A (~1.9 GB),
over the node memory ceiling (MaxRSS 29 GB at kill).  The apply path
prefers JFNK_BANDED_LU_READY over the requested preconditioner, so a
surviving run would also have silently applied banded LU instead of
Jacobi.  Fix: banded factorization now gated on PRECOND_TYPE='LU';
JACOBI/LSQ/SCHUR runs skip the assembly and apply what was asked.
Band storage for the analytic matvec (CMAT_RAD_BANDS) is armed
independently in JFNK_BUILD_BA_RAD_V1, so the banded-JVP experiment
under JACOBI is unaffected.  Regression added (64/64).

22-Oct-2026 (6): rerun as 12315454 on a57e78ce.  Full banded analytic
JVP under JACOBI now runs: band maxima printed on first JVP --
diag 4.46e18, lower/upper off-diag 9.7e12/6.7e12 (BA_ADV_TERM
advection + any ALO content; ~6 orders below the stiff-row diagonal
max).  GMRES with the complete banded operator: res 4.65e10 ->
4.43e10 (~5%), flat after it~30, hess_res == true_res = 1.0000 ratio
at both restarts -- the matvec, residual eval and Krylov bookkeeping
are consistent for the full operator.  Still capped by the hardcoded
restart limit (2 cycles => 60 its; MAXGMRES=300 never reached).
Step rejected safely (|s|_inf=2.6, |s|_w=3.4e7) -> Lambda fallback,
model continues.  Interpretation per the experiment table: operator
consistent, plateau persists -- cross-depth coupling does NOT explain
the diagonal-band plateau.  The missing reachability is not in the
dropped bands.  Next per plan: attainable-residual/conditioning
study (lift the restart cap to exercise the full budget) before any
residual-reformulation call.

22-Oct-2026 (7): configurable GMRES restart cap.  gmres_m_v1 hardcoded
'n_restart > 1' -> exit, so JFNK_MAXGMRES=300 never ran more than two
m=30 cycles.  New VADAT key JFNK_GMRES_RST (default 1 = legacy
behaviour) wires to module var jfnk_gmres_max_restart; imx still
bounds total iterations and the true-residual diag prints at every
restart.  Regression added (65/65).

22-Oct-2026 (8): convergence-bakeoff branch (off development,
jfnk-stability Phase-3A left on its own branch).  Findings 1-6
verified; DIAG_RDG shadowing fixed (bc7fb139).  Toy NLTE testbed in
tests/convergence_bakeoff/ links the production newton_hook_v1
unchanged; 32-instance ensemble split 11 TUNE / 21 HELDOUT before
tuning; DECISION_RULE.md committed (7c07dbe8) before the held-out
run.  Toy bug found+fixed: factor_ba_red read LU factors, not the
unfactored BA (ba_raw copy).  Held-out result: reduced-space trace
elimination (+a) is the decisive ingredient -- every +a variant
21/21 at median 12.0-14.3 sweep-equivalents vs stock 28.7
(2.0x-2.4x); no non-+a method beat stock.  Bare s2f went 16/21 at
equal cost (all stagnations on trace rows); s2c budgeted out 18/21
(~1 formal solve per JVP).  Recommendation (toy only): and2+a
(Anderson-2 + exact trace elimination), ranked ahead of
s2f+b+c1+a on fewer knobs at a 0.8-SE median gap.
Noise: the toy cannot reproduce the production |FD|/|BA v|=36x
discrepancy -- deterministic eval noise (Jacobi-count switching,
clips, cancellation at stiff=1e18) stays clean to h=1e-12 even at
the fixed point.  New VADAT JFNK_NOISE_DIAG (default F) logs
per-eval |dY| at |dX|=0 on the production residual to separate OMP
reduction-order from unbracketed sweep-state noise (runbook in
docs/audits/JFNK_CONVERGENCE_BAKEOFF.md).  Production port:
JFNK_TRACE_ELIM (VADAT, default F) -- GS form of +a: activation-
frozen BA_DEM_PFLOOR exempt-row mask applied through
JFNK_MASK_FIXED_UNKNOWNS so residuals/JVPs/preconditioned vectors
and the committed step stay in the major-row unknown set; Stage-B
BA_TRACE_PHASE polish resolves traces post-step.  UNVERIFIED on
real models; needs-human-validation runbook in the audit doc.

04-Nov-2026 (5db4cd5d, test fix 0fe26b1a): Phase-1 JVP diagnostics from the numerical-methods
review -- three opt-in VADAT keys, all default off, default path
unchanged.  JFNK_LINTEST runs the JVP operator-consistency probes
(repeat/homogeneity/additivity, 4 evals once per activation on the
first post-factorization direction, same jfnk_x and effective h).
JFNK_HSCAN forces the one-shot extended JVP block independent of the
10% FD-vs-BA trigger and adds the cancellation-floor SNR table per
worst row.  JFNK_TRES_FAC (>0) turns the RESTART DIAG true-residual
measurement into a GMRES acceptance criterion: a claimed res <= tol
must satisfy true_res <= fac*tol or the claim is rejected (restart
on the measured residual, else jfnk_exit_reason='gmres_tres_reject').
Solver-module name is jfnk_gmres_tres_fac (JFNK_TRES_FAC would
collide case-insensitively under the unqualified USEs in
solve_corr).  Unit tests T4a/T4b added to
tests/jfnk_core/test_newton_core.f90 (accept on consistent matvec;
reject on affine-offset matvec, hess 3e-11 vs true 4e-10).  Harness
verified T1-T4 + guards + ptc_sign; test_precond_blocks had a
pre-existing rand()/srand() failure under ifort (gfortran extensions)
-- fixed with '#ifdef __INTEL_COMPILER use ifport' in the five
routines that call them.  Full suite now green: tests/jfnk_core
4/4 binaries and unit_tests/python/test_jfnk_solver.py 78/78.

07-Nov-2026: external-audit fixes on the Phase-1 diagnostics
(follow-up to 5db4cd5d; review of the pushed branch).  Three defects
fixed: (1) LINTEST synthetic direction scaled to |u|_inf -- h was
chosen for the production direction's norm, so a unit-vs probe
measured clip branches, not operator nonlinearity.  (2) The
HSCAN/extended block now fences jfnk_resid_nbad, JFNK_N_NEW_CLIP and
JFNK_EVAL_BROKEN the same way LINTEST already did -- a broken 1E5*h
sweep trial no longer leaves status flags behind.  (3) The
JFNK_BA_JVP analytic path returned before the NaN/Inf guard; it now
applies the identical >10%-fails policy.  Also fail-closed on a
non-finite true_res in the tres gate (NaN > acc_tol is .false. and
would have been accepted).  T4 tests extended to T4a-T4f: exact
exit-reason asserts, restart-vs-no-headroom eval-count comparison,
left-preconditioned accept/reject with a non-cI diagonal M, and a
NaN-in-tres-matvec rejection test.

07-Nov-2026 (re-review follow-up): non-finite TRES no longer restarts.
A second review of 92f52a99 found that a NaN/Inf true_res took the
ordinary reject path -- restarting the GMRES cycle on an operator
that had just produced a non-finite action.  gmres_m_v1 now exits
immediately with jfnk_exit_reason='gmres_tres_nonfinite' (info=1,
non-converged) using ieee_is_finite, which classifies NaN and Inf
alike.  T4f updated to the new reason; T4g added: NaN on the tres
call WITH restart headroom asserts the same immediate exit and
matvec_ct == nan_at (zero matvecs after the measurement -> no
restart ran).

10-Nov-2026: external-audit round 3 (review of 92f52a99, applied on top
of 17b2759d).  Four changes, all opt-in or fail-closed:

1. Non-finite TRES generalized to ALL done-exits.  The 17b2759d guard
   only fired on claimed convergences (res <= tol); a stagnation exit
   that bypasses into the lb_hit restart could still re-enter the
   cycle on a contaminated measurement.  ieee_is_finite(true_res) now
   exits 'gmres_tres_nonfinite' at every done-exit.  T4h/T4i added:
   its==imx exit with injected NaN and +Inf, exact reason asserted.
2. JFNK_PRECOND_V1 final-output guard.  The per-depth apply only
   guarded its own intermediate result; the log/PTC rescale, T-ladder
   scaling and enrichment backsolve ran afterwards unguarded.  New
   pure routine jfnk_pc_final_guard_v1 (zeroes non-finite entries,
   counts) is called just before the apply counter, on the diag-only
   early return too, and raises jfnk_pc_clip_event so GMRES aborts.
   G10/G10b unit tests (test_guards.f90).
3. rcond gate fails closed: jfnk_pc_factor_block_v1 now requires
   ieee_is_finite(rcond_k) so a NaN DGECON estimate cannot pass the
   '< rcond_min' comparison vacuously.  T-RCOND-NAN added; caveat
   noted in the test header: real LAPACK lands poisoned inputs on
   rcond=0 or info_g=-5 anyway, so the clause covers BLAS variants.
4. JFNK_HADAPT (default F): adaptive FD step.  Reuses the one-shot
   h-sweep; jfnk_h_select_v1 picks the multiplier with the best
   adjacent-pair JVP agreement (30% tol) = plateau between truncation
   drift and cancellation noise.  Selected multiplier scales the
   post-cap H_STEP of later JVPs (JFNK_HADAPT_MULT in int_mod); the
   returned image is the sweep value at the selected step.  No
   agreeing pair -> jfnk_jvp_failed -> Newton defers.  G9a/b/c cover
   plateau pick, noise rejection, bad-eval skip.

Also: JFNK_REFINE and JFNK_PC_PAR_FACTOR confirmed dead (parsed,
unread).  VADAT text marked INACTIVE, JFNK_REFINE default T->F, and
the JFNK config banner warns once when either is set.

Harness: test_newton_core T1-T4i + test_guards G1-G10 +
test_precond_blocks (incl. T-RCOND-NAN) + test_ptc_sign all pass;
test_jfnk_solver.py 78/78.

22-Oct-2026 (late): plateau-state warm restart of the a397b326 audit
fixes submitted for live validation.  Run dir:
cmfgen_validation/jte_hadapt_p1 -- full clone of jte_AB_base (m15
2.82d, DO_JFNK=T, the bakeoff testbed that plateaued at 88% /
GMRES-stagnation+guard-defer at it506).  VADAT additions vs the
baseline clone: JFNK_LINTEST=T, JFNK_HSCAN=T, JFNK_HADAPT=T,
JFNK_TRES_FAC=2.0D0, JFNK_GMRES_RST=4.  All other JFNK gates,
budgets and admission controls identical to baseline.
Binary caveat: exe/cmfgen_dev.exe (17:27 build) = commit a397b326
PLUS uncommitted parallel-session changes in new_main/cmfgen_sub.f
(SOB dispatch-order hardening, foreign work) -- the tree was dirty
at build time, so the binary does not map exactly to the commit.
Jobs: interactive probe 12331116 (2h cap, viper07, started
19:47), chained small-partition continuation 12331117
(--dependency=afterany).  Success criterion set by user: JFNK makes
a measurable step towards convergence.  Notes: arming counters
reset on warm restart (hits frz=0/2 at it507); post-restart demand
transient 62->47->21 must fall under FRZ_GATE=3.0 for 2
consecutive its before JFNK re-admits.

22-Oct-2026 (cont.): jte_hadapt_p1 first live results.  Interactive
probe 12331432 (2h, gates relaxed JFNK_FD_AUTODEF=F/JFNK_FD_COS_MIN=0
-- required, else the admission gate defers upstream of JFNK_JVP_V1
and no diagnostic ever runs).  ~15 JFNK-FRZ activations it529-560:
HADAPT selected h x1E5 (heff=1.2920E-03, pair drift 8.4E-04) from
the plateau between truncation drift and the cancellation floor
(default heff=1.29e-8 sat at FD-vs-BA rel-diff 17.4, pure noise);
LINTEST at selected h: repeat exact, homog 1.9E-4, addit 7E-9.
TRES fired every activation: true_res/hess_res ratios 1.08x-3755x
(meaningful when the assembled step is sane, exploded when the
stagnated GMRES direction hits nonlinear walls).  Newton: 1 line
search accepted (alpha=1.2E-4, F unchanged); all other activations
abort (gmres_stagnation / too_many_rejected_steps).  Budget rung
pinned 4/4 (gmres<=15 evals<=12 per activation -- the rung ladder
caps budget below one m=30 cycle, so JFNK_GMRES_RST=4 was never
exercisable).  Verdict on FD step: SOLVED -- JVPs at the selected
h are linear and BA-consistent to ~1.5%.  Verdict on plateau:
preconditioner/GMRES is now the proven bottleneck -- residual
concentrated in ~10 HI lev 21-24 rows at depths 65-70, per-depth
LU cannot represent the coupling GMRES needs.

PTC_SIGNED=T live A/B (small job 12331117, started 22:12): per-depth
rcond at the compression-zone blocks improved ~10-1500x vs the
unsigned segment (d63 3.6E-7->4.1E-5, d64 1.6E-6->2.2E-4,
d70 1.2E-7->1.8E-4) -- direct evidence the unsigned PTC shift
degraded near-singular blocks.  GMRES outcome under signed PTC:
still stagnating in early activations; monitoring continues.

22-Oct-2026 (cont.2): jte_hadapt_p1 signed-PTC segment complete
picture (it563-592+, small job 12331117, commit a397b326):
38 JFNK-FRZ activations -> exits: 15 gmres_stagnation_abort,
21 too_many_rejected_steps, 2 trust_region_too_small, 3 line
search accepted.  All three acceptances negligible at print
precision (|F| unchanged: alpha=1.2E-4/|F|=3.3E10,
alpha=0.5/|F|=1.6E13 vacuous, alpha=3.9E-3/|F|=5.6E10 with
ared/pred=1.0005).  Signed PTC gains confirmed: rcond at the
compression depths 10-1500x better, step norms sane
(|s|_inf 0.03-0.08 vs 0.45-2E8 unsigned), one activation with
true_res = 1.004x hess_res (operator-consistent direction --
operator fine, GMRES is the limiter).  Pinned rows carry only
0.14% of |F| -- the residual lives in well-conditioned rows, so
the ~15-18% GMRES stall is NOT a regularization artifact: the
per-depth blocks cannot invert the cross-depth radiative
coupling concentrated at depths 63-70 (rad_i ~= coll_i on the
dominant HI lev 21-24 rows).

Parallel LSQ probe (interactive job 12333052, run jte_lsq_p1,
clone of jte_AB_base it506, JFNK_PRECOND=LSQ + signed PTC +
same diagnostics): first activation it520 -- the equilibrated
block LSQ fell back to Jacobi on 68/70 depth blocks
(ok=2 ridged=2 Jacobi=68), i.e. effectively unpreconditioned;
GMRES still reached res 3.33E10 -> 2.56E10 (23%) in 27 its,
similar to LU's ~15-18%.  |F| after solve unchanged ->
too_many_rejected_steps.  Interpretation: changing the block
SOLVE METHOD does not move the stall fraction -- the wall is
preconditioner STRUCTURE (missing cross-depth blocks), not
block conditioning or FD accuracy.  Remaining structural
options are banded LU (requires NUM_BNDS=3 + SOL_METH=TRI,
which changes the model population solve -- a model change,
not a VADAT flip) or the coupled SE+RTE tier.

Verdict vs user criterion "measurable step towards
convergence": JFNK machinery now verified correct and
instrumented (FD step resolved by HADAPT, GMRES overclaim
exposed by TRES, signed PTC fixes block conditioning, real
step attempts now occur vs zero acceptances in baseline
~500 its), but the outer plateau is NOT broken -- model
convergence oscillates 75-88% identical to baseline, demand
spikes persist, zero committed Newton progress in 38 signed
activations.  Report: PASS (diagnostics/mechanism) / FAIL on
plateau break -- blocker is now measured: cross-depth
coupling absent from the per-depth preconditioner.


## 23-Oct-2026 - Banded TDMA, BA-model falsification, secant update

Follow-on to the structural-wall finding above.  Three parallel
threads:

(1) JFNK_BANDED_TDMA (commit 19fa7a6b): block-Thomas factorization
of the NUM_BNDS=3 banded preconditioner replaces the infeasible
DGBTRF (~3.5E15 flops -> ~1E11).  NUM_BNDS=3 is a MODEL_SPEC key
(BA band storage) so METH_SOL=DIAG is preserved; warm restart
skips the NUM_BNDS=1 BAMAT and rebuilds.  Run jte_btd_p1 (small
job 12333785, clone jte_AB_base): factorization algebraically
exact (backward err 5E-14..1E-13 every activation) but GMRES
still stalls ~22% and every activation exits
too_many_rejected_steps -- identical to per-depth LU.  Offline
BAMAT band-norm parse: off-diagonal share is ~1E-6 of the
diagonal at the residual-dominant depths (63-70), i.e. nearest-
neighbour coupling is negligible exactly where it matters.
Superblock generalization added (band half-width MB depths =>
block-tridiag in MB*NT superblocks; NUM_BNDS=5/11/15 reachable,
ND=70 divisible), unexercised as of this entry.

(2) JFNK_BA_JVP probe (jte_bajvp_p1, small job 12334254): the
analytic full-band BA matvec as the JVP isolates the solver from
FD noise.  Result: with the TDMA preconditioner, GMRES drops the
Hessenberg residual ~1E4 in ONE iteration -- the preconditioner
is an excellent inverse of the banded BA operator (machinery
proven) -- but the true residual is unchanged (~3.5E10) and all
activations exit gmres_stagnation_abort.  Interpretation: the
banded BA is simply not the Jacobian of the residual on the
dominant rows.  JVP diagnostic rows show FD action ~1E14 vs
BA*v ~3E9 on H I lev 21-24 / He II rows at d66-69 where
max|BA_row|~1E19, and the dF response is non-smooth in h
(floor/clip crossings).  Three independent measurements agree:
BAMAT band shares ~1E-6, FD-vs-BA ~4 orders, BA-directions
don't reduce |F|.  The wall is the BA linearization model
itself, not the solver.

(3) JFNK_SEC_R (this commit): rank-r Broyden/secant preconditioner
update -- since the measured (u_i, J*u_i) pairs are produced free
by every JVP, the correction makes the applied inverse exact on
the recorded span: B'v = Bv + D*(G^-1 JV^T v), D_i = u_i - B(J*u_i),
G = JV^T JV (1E-12 ridge).  Frozen-tier analogue of the coupled-
tier JFNK_ENRICH idea, but pairs come from the GMRES directions
themselves (no extra evals).  Routines jfnk_secant_gram_v1 /
jfnk_secant_apply_v1 in jfnk_solver_v1.f90; recording in
JFNK_SEC_RECORD_V1; build+apply inside JFNK_PRECOND_V1 (recursive
base apply under JFNK_SEC_BASE); pairs cleared on baseline
rebuild, D rebuilt on factor rebuild.  Unit test T-SEC passes
(exact inversion on recorded span, err ~1.7E-10 at the Gram
ridge).  First run: jte_sec_p1 (SEC_R=15 + TDMA + signed PTC +
HADAPT).  Caveat to watch: the growing correction makes the
preconditioner non-stationary within a solve (no FGMRES); TRES
is the convergence authority.


## 23-Oct-2026 (cont.) - REG_ZERO_COL breakthrough; runaway-cap fix

Two decisive findings from the live probes:

(1) JFNK_REG_ZERO_COL=T (jte_rzc_p1, interactive job 12335859, LU
base + SEC_R=0): the factor-time zeroing of regularized columns --
present since audit commit e0ff8d53 but never run live -- collapses
the M^-1 round-trip error from ~1.4E12 to ~4E-3 on measured JVP
directions.  With the corrupted reg-column coupling removed, GMRES
converges the linear system for the first time in this model:
hess_res 4.7E10 -> 4.5E7 (~1E3 reduction) in 27 its, repeatable
across activations and budget rungs.  FD-vs-BA now agrees within
~28% on probe directions -- the BA model is a reasonable Jacobian
once the factorization stops contradicting itself.  Secant updates
(jte_secLU_p1 12335550, jte_sec_p1 12335911) built correctly to
rank 14-15 but cannot compensate for the corrupted base apply:
rt stays 1E2-1E5 and GMRES stalls ~5%.  Conclusion: REG_ZERO_COL
is a prerequisite, not an option, for any BA-derived
preconditioner.

(2) Runaway-cap miscalibration (commit f307b544): the RZC run
exposed that jfnk_runaway_nrm=1E2 is an ABSOLUTE weighted-step cap
calibrated for Euclidean RSCALE.  Under JFNK_RSCALE=DIAG,
|F|_w ~ 5E10 and every Krylov step (|s|_w ~ 3E10, ~0.6x ref_norm)
read as runaway.  The solver burned all 8 hookstep retries halving
del and aborted -- "Newton-level F evaluations: 0", i.e. the trust
region and line search never once evaluated a trial step in this
model's history.  Fixed by making the cap relative:
|s|_w > jfnk_runaway_fac * ||M^-1 F||_w (fac=1E2 preserves the
catastrophe semantics for the 1E31 singular-inverse class).  Same
fix on the Cauchy-direction check.  Validation run: jte_rzc2_p1
(job 12336336), same config as RZC + fixed exe.

(3) Residual anatomy (all runs): |F|^2 share in top 1% rows = 1.00;
dominant rows are H I lev22/23 at depths 67-70 + He II lev13
(F_i ~ -3.8E12, coll-dominated, NOT pinned).  The converged GMRES
step demands maxch ~ 2E10 on these rows; even |s|_inf~1% trial
steps were untested under the old cap.  The remaining question the
rzc2 run answers: does the now-functional trust region find ANY
accepted step, or are the dominant rows genuinely bound-limited
(population-floor optimum reachable only non-locally)?


## 23-Oct-2026 (cont.2) - First verified Newton commits; phantom-baseline defect

jte_rzc2_p1 (12336336, relative-runaway exe): first activation
reached the line search for the first time in this model's JFNK
history -- ared/pred = 1.0001, alpha = 3.8E-6 accepted, |F|_w
5.21e10 -> 5.00e10 -- but the second Newton iteration exhausted
the 200-eval budget and the commit-time applied-state eval then
read the budget sentinel (|F|_w ~ 1e46 = 1e20 * sqrt(sum RSCALE^2))
and rolled back the verified step.  Bypass fix (21b52cae):
JFNK_EVAL_BYPASS lets the single commit-verification eval through
the budget/wall-clock sentinels.

jte_rzc2_p1 resubmitted (12336793, bypass exe): THREE consecutive
verified commits -- it598 5.80e10 -> 5.38e10 (alpha 1.5e-5), it599
-> 5.15e10, it600 -> 5.11e10 (alpha 1.9e-6), each confirmed by the
applied-state check.  Entry-residual trajectory across admissions:
5.80 -> 5.38 -> 5.15 -> 5.11 -> 6.32 -> 4.67 -> 4.56e10 (net -21%
over 6 admissions; the 6.32 bump is BA-sweep interference between
activations, not JFNK divergence).  JFNK now admits on consecutive
great iterations.

New defect found at it602 (fixed in 1f56955a): after
max_newton_iterations, the attempt loop calls
JFNK_REBUILD_BASELINE_V1(BEST_X), which folds CMAT_RAD*X_CORR into
STEQ_RAD_FLAT and re-bases F_W_INIT.  The rebuilt |F(0)| = 1.49 is
the LINEAR model's prediction, not a measured residual -- a
phantom against the true ~6.04e10 the line search measured for
that state.  The phantom triggered converged_at_entry on the next
attempt, and the applied-state check compared moved-operator
applied (1.49) vs moved-baseline F_W_INIT (1.49): "no improvement"
-> rolled back a verified-improved state.  Fix: F_W_ENTRY
snapshots the activation-entry residual and is the sole commit
reference; F_W_INIT keeps its re-basing role in the attempt loop.
Residual risk accepted and documented: post-rebuild F_APP is
measured on the re-based operator, but the committed state is only
reachable via a measured BEST_F_W < entry gate.

Economics now the dominant question: each activation nets ~2-7%
|F|_w via alpha ~ 1e-5..1e-6 accepted steps, while ~70% of the
200-eval budget goes to TRES-rejected GMRES restart cycles whose
Hessenberg residuals explode (hess_res -> 1e12 while true_res
stays ~1e11-1e12).  JFNK admits consecutively, so progress
compounds; whether it outruns BA-sweep interference is what the
resumed run (12338500, entry-residual fix) measures.


## 12-Nov-2026 - True-residual gate + zero-progress rejection

Two confirmed defects from the pre-bypass OUTGEN fixed on
jfnk-precond-fixes (f8418465, 89d62588, tests 537b4da1):

(1) The Newton stagnation gate tested only the Hessenberg
least-squares residual, which is a Krylov-space MODEL residual.
Under the finite-difference matvec the basis need not span a fixed
operator; OUTGEN showed hess_res=2.85e10 vs true_res=4.30e11
(it=22) and hess_res=5.39e8 vs true_res=5.05e12 (it=90) on exits
the gate accepted as inexact directions (pred_res/ref=0.575).
gmres_m_v1 now measures ||M^-1(b - A x)||_w exactly once on every
info=1/2 exit (stagnation, max-iteration, restart-cap), reusing a
TRES or restart-diagnostic measurement when one already exists, and
exports it via jfnk_gmres_true_res/_ok.  newton_hook_v1 requires
BOTH pred_res < 0.9*ref_norm and a finite measured true_res <
0.9*ref_norm before accepting an info=1 direction; a missing or
non-finite measurement is unverifiable and aborts.  No extra
matvec on converged exits.

(2) An ared == 0 trial fell through to the rho < 0.75 retry branch,
which never touches n_reject: OUTGEN showed del halving ~6 -> 1e-8
over hundreds of zero-progress evals until trust_region_too_small,
each labeled 'ACCEPTED (rho<0.75)'.  ared <= 0 now takes the
rejection path with n_reject accounting, so no-progress ends the
solve in n_reject_max (=8) passes as too_many_rejected_steps; the
disposition label reads 'REJECTED (no decrease)'.

Regression tests (test_guards.f90): G11 drives newton_hook_v1 with
a deliberately inconsistent matvec (rank-34 projection on unit
basis calls, zero action on the accumulated-step call): hess_res
plateaus at 0.66*ref while true_res = ||b|| -- the pre-fix gate
accepted, the fix aborts as gmres_stagnation_abort.  G12 uses a
constant residual so every trial reports ared == 0; the solve must
abort via n_reject_max rather than loop rho<0.75 to the del floor.
All T1-T4/G1-G12/precond/PTC suites pass; interface suite 9/9.

23-Oct-2026 (sync/dev-stageb-hardening branch in cmfgen_jfnk_dev):
port of the four JFNK fixes from this branch (F_W_ENTRY commit
gate, measured true-residual stagnation gate, zero-progress
rejection, G13/G14 regressions) plus the Stage-B off-state rework
and production diagnostics, on top of origin/development 183f0dd.

Commits: 5a9fd24b F_W_ENTRY (trailers stripped), 4bcbadb2 true-res
gate (adds an early _ok/_res reset at gmres_m_v1 entry -- the
alloc-failure returns at :2009/:2045 precede the per-exit resets and
could have left a stale measurement for the Newton gate), ef9c9b8a
ared<=0 rejection, 4ff9a166 tests (this branch's G11/G12 land as
G13/G14; dev's own G11/G12 eval-cap tests kept), 3e41b1f4 Stage-B
state machine (shared BA_TRACE_EVAL in control_variable_mod; OFF
releases the exit, MAXSWP=0 master off, ABSMAX=0 kills only the abs
leg; test_trace_state BA1-BA8 links the production routine),
f64291e5 EXR export from committed state (was stale from the last
BA call on JFNK commits), 7bb565e6 arming line + final-status
branches (no more 'Convergence criteria met' + 'still owed'),
9af93895 exempt/trc + ALO field widths (ESn.nE3 forces the E on
three-digit exponents), 68a5b78d rebuild line reports info/iter.

Verification: make clean && make test in tests/jfnk_core -- all
suites green incl. new BA1-BA8 truth-table tests.  Live runs queued
on small: mom_stgb 12340500 (demand-arm), jte_stgb 12340501
(JFNK+arm), stgboff 12340849 (MAXSWP=0), stgbabs0 12340850
(ABSMAX=0), stgbswp1 12340851 (MAXSWP=1).

23-Oct-2026 (cont.): branch pushed fast-forward onto
origin/development (183f0dd4..ceaabab7) -- user directive was to
have the code on the development tree itself, not a PR.  Two
follow-ups: ceaabab7 stageB sweep counters I4; 26393245 comment
docs for res dual-use / jfnk_gits recompute semantics / dpotrf
fallback metric (GMRES audit reconciliation).

GMRES-breakdown follow-on (931a85eb): newton_hook_v1 now skips the
eta-tighten descent re-solve when the physical linear model already
shows amplification (||F - J*s||_w > ||F||_w) -- observed in OUTGEN
as hess_res 1.6e7 vs |F-Js| 5.0e15; the retried solve failed
identically and burned a full GMRES budget.  The step still goes
through the trial eval (ared <= 0 -> reject -> free hookstep
recompute).  G15 regression test (df_anti: definite basis action,
anti-aligned accumulated-step action, imx=2) asserts
jfnk_n_descent_retry == 0; pre-fix it burned 1+ re-solves.

Validation VADATs: all five pending stgb* jobs + jte_rzc2_p1's
slurm now point cmfdist at cmfgen_jfnk_dev (were erroneously at the
main worktree / jfnk branch).  rzc2 job 12338500 hit its 2h limit
mid-JFNK-activation.  jte_stgb_int submitted on interactive
(12341448) restarting from the stgb_p1 checkpoint with the new
exe; early OUTGEN shows clean E3 exponents and correct
below-threshold non-arming (demand 0.63% < 1%).

24-Oct-2026: JFNK restart-stall fix set on coupled-sob-jfnk.
Root causes found in job 12351111 (jte_lgmres_ab, m15 2.82d):
(a) `jfnk_lgmres=JFNK_LGMRES` in jfnk_solve_corr_v1.f was a
case-insensitive self-assignment -- JFNK_LGMRES=T never enabled
LGMRES in production (fixed to `jfnk_lgmres_on=JFNK_LGMRES`);
(b) restarted GMRES rebuilt a near-invariant Krylov space each
cycle (restart-residual cosine 0.99); (c) a late JVP failure
discarded the entire Newton direction (~1 hr of coupled evals).

Changes in jfnk_solver_v1.f90 / jfnk_solve_corr_v1.f /
control_variable_mod.f / rd_control_variables.f:
- JFNK_GDEFL/_K: GMRES-DR harmonic-Ritz deflation
  (jfnk_gdefl_extract_v1; DGESV H^-T e_j solve + DGEEV; smallest
  |theta| pairs; real spans for conjugate pairs; V*g left-PC /
  Z*g right-PC directions; images via Arnoldi relation -- no
  extra JVPs; shares the <=8 augmented-LS slots with LGMRES).
- JFNK_CYCSTALL/_FAC/_N: per-cycle stall early exit (info=1,
  reason gmres_cycle_stall, partial step kept).
- JFNK_RECYC: optional cross-solve augmentation retention (off).
- JVP soft-fail: later-column failure closes the cycle on the
  completed basis -> info=1 gmres_jvp_soft_exit, partial step
  rescaled to the trust radius and true-residual measured.
  First-residual / first-column failures stay info=3.
  jfnk_jvp_failed cleared before verification matvecs.
- Driver: LGMRES now actually enabled; combined aug capacity
  computed; new options echoed in the one-time JFNK CONFIG banner.

Tests: tests/jfnk_core all green incl. T7e (GDEFL), T7f/T7g
(mid-cycle + verification JVP soft exits), G7 (first-column hard
abort).  Synthetic stiff problem: plain 1.10 vs LGMRES+GDEFL 0.13
true_res at equal JVP count.  ./build.sh clean; exe 00:43.
test_jfnk_solver.py se:hard/T11 + se:shock/T12-13 failures are
pre-existing (fail identically at HEAD; rung/tauto replica drift,
unrelated to the GMRES changes).

Live: jte_lgmres_ab resubmitted as 12352827 (warm, it 676) with
LGMRES=T GDEFL=T(4) CYCSTALL=T(2E-2,2) RECYC=F.  Convergence and
the required cold-start retest are pending.

24-Oct-2026 (cont.): a2d501bb fixes the restart-cap residual audit.
After GMRES scales the accumulated step onto the trust radius, the
solver now invalidates the pre-scale cached norm and recomputes the
true residual for the returned step.  G5 reproduces the old mismatch
(1796.6 cached vs 2088.5 actual) and verifies the post-fix norm.

Harness repairs in the same commit restore the test contract rather
than changing physics: test_jfnk_solver.py accepts the five current
core verdicts; se_driver dispatches PC 0/1/3/4/5/6/7 to the selected
preconditioner instead of silently taking LU; hard/T11 uses the TSVD
committed-best-iterate path (info=2); shock/T12-T13 check the
adaptive-ladder invariants from the observed RELEASE/REVERT/PERT
transitions; T19 checks the budget-rung fallback path.  Clean
./build.sh: 17/17 Python test files passed.

Live: npq_dbg resubmitted through run_model_grid.py on interactive as
FINAL 12362078 (OBS dependency 12362079), warm_restart_from
m15_2p82d_SOB_l14_lte_np/cmfgen/2.82days.  The checkpoint supplies
POINT1/POINT2/SCRTEMP/MODEL with matching NT=1839, ND=70.  Model
convergence is pending.

24-Oct-2026 (later): compact-band diagonal-offset bug fixed.  The
[JFNK_BANDED_TDMA] allocation shrank LUAB to KL+KU+1 rows (commit
d32bf1ca) but every index still used the DGBTRF diagonal offset
KL+KU+1+i-j, so the guarded assembly dropped the entire lower
triangle and all lower superblocks.  Job 12366624's own
"superbands: max|L|= 0.00E+00" line was the fingerprint; the
5.08E-14 backward error was vacuous (i>j reads fell into the zeroed
headroom).  New local DIA_ROFF in JFNK_FACTOR_PRECOND_V1 selects
KL+KU (DGBTRF) or KU (compact); all LUAB/ARAW indexing now goes
through DIA_ROFF+1.  Non-TDMA path is bit-identical.  Warm
diagnostics npq_dbg/sec/npq_gd/fds (all flat 3.6-4.7 dex, every
admission info=3) and btdc cancelled; btdc resubmission pending on
a clean ./build.sh.  jte_lgmres_cold (12358056) left running.
See convergence_issues.md 24-Oct (cont).

24-Oct-2026 (evening): two more operator-layer defects found and
fixed; first committed JFNK Newton steps ever.

(a) 106f138e verified on btdc (12368192): correct compact layout
gives max|L|=1.00 and GMRES convergence in 4 its -- but the returned
step was still non-descent (cos(F,Js)=-0.08, |Js|~45|F|).  The
DIAG-weighted metric hid the stiff rows.

(b) 7b09326b extended JFNK_PRECOND_AUDIT_V1 to run under
JFNK_BTD_READY and added a banded J*M^-1 check (same analytic
matvec the JVP uses on z=M^-1 F).  bte/btf measured
|Jz-v|/|v|=3.6E+05, |Jz|=8.0E+17, cos=-0.14: the factored banded
operator was NOT the JVP operator even though the factor's own
backward error was 3.7E-14.

(c) 1f19bd81 root cause: LUAB kept the real energy/charge rows and
dF/dT (dF/dNe) columns while the JVP operator has them masked out
(BA_A frozen-row reduction + JFNK_MASK_FIXED_UNKNOWNS on both
sides).  The banded solve mixed masked directions into the
population step.  Mirrored the frozen-row reduction onto every
LUAB band block (zero T/Ne row+column, identity on diag block).

Verified on m15_2p82d_SOB_l14_lte_btg (12369040/41, EUCLIDEAN,
analytic JVPs, warm it735): J*M^-1 = I to 3.0E-13, cos=1.0, GMRES
converged in 1 iteration, cos(F,Js)=0.9998, ared/pred=0.77, five
consecutive ACCEPTED Newton steps -- the first committed JFNK
progress in the lineage.  |F|_w moved 5.86E+12->5.84E+12 (0.23%)
before the trust region starved: the Newton direction asks
|s|_inf~848 on dead rows and the +/-2.30 log clip (MAX_JFNK=10)
allows only ~5 ratchet steps per admission.  btj (12369227/28)
resubmitted with MAX_JFNK=1.0D6 to test whether larger dead-row
corrections drain the residual faster.  bte/btf cancelled after
their diagnostics were recorded; btdc kept as the pre-fix control.
(d) btj (12369227/28, MAX_JFNK=1.0D6 -> +/-13.82 log clip) answered
the clip question: the d26 dead row took its full 1e6-factor jump,
first-step ared rose 6x to 7.4E+10, but the admission still floored
at |F|_w 5.78E+12 (1.26% vs 0.23%).  After that move every trial
predicted negligible decrease -- the residual sits in the near-null
space of the frozen-J operator (dead rows that do not move F, counted
in full under BA_EX_RELFRAC=0).  The clip was a throttle, not the
wall.  Remaining question is outer-loop value: whether committed
JFNK corrections shorten the BA/Lambda path to margin over repeated
admissions.
    btj follow-through: the 11.99 dex spike never healed --
    OIV:2p3_2Do[5/2]@d16 pegged at the 1E7% demand ceiling (applied
    1E2), se residual inflated to 8.5E+13.  btj cancelled at it~749.
    Verdict: the +/-2.30 clip is protective, not just throttling;
    keep MAX_JFNK=10.  btg's second admission cut |F|_w 16.3%
    (5.11->4.28E+12, nine ACCEPTED steps, ared/pred~0.90) and the
    outer margin recovered to 4.04 -- the corrected path makes real
    per-admission progress.
(f) JFNK_CAP_STEP (commit 161c95a9, default F): per-component cap on
over-bound Newton/Cauchy step components replaces the uniform
inf-norm rescale, so a few dead-row demands no longer shrink the
whole direction.  T8 unit test (decoupled underpopulated row):
|F|=1.0e3 vs 2.8e3 uniform over 25 its.  Live test btk (12369984/85)
= btg + CAP_STEP.  Margin anatomy resolved: ~7 dex spikes are JFNK
commit MAXCH feeding BA_EX_MAXCH (decays in ~2 its); ~4.2 dex
baseline is the outer BA raw demand.  CORRECTION_SUM shows >10%
empty but ~1100-1400/1839 rows >0.1% -- 'mostly zeros' is outer-
loop work.
(g) btk first live result (job 12369984): JFNK_CAP_STEP cut |F|_w
62.3% in ONE step at it743 (btg's uniform-scaled best was 0.23%).
Post-commit kick: margin 8.2->11 dex, one transient dead row
(NeV 3p_3Se d26, -1E7%) that SELF-RESOLVED in ~6 calls (not the
btj permanent wall).  Baseline now trending down 8.17->4.07 dex vs
btg's flat ~4.1-4.9.  CORRECTION_SUM better than btg on every
column.  DIAG weighting disproven a third way: btdc (fixed
operator, |Jz-v|=6.5E-15) still rejects every DIAG admission -
the weak-row-dominated weighted problem produces non-descent
directions; use JFNK_RSCALE=EUCLIDEAN.  The ~6.4 dex 'full-call'
margins are the 1.5E5 BA_MAXCH_FAIL rejection sentinel, not a
residual - they clear once no row demands >=100% committed change.
CORRECTION_SUM histogram counts raw STEQ for ALL rows (PFLOOR and
RELFRAC do not remove entries; only quarantine or true convergence
does).  Cold run at 0.22 dex, nearly converged.

## 27-Oct-2026 - DO_ES_ALI on m15_2p82d_SOB_l14_lte_adp2_esali (commit 986a7b48)

Committed the three gated SOB/LAM_ES stall probes (all default F, stock
path unchanged): DO_SOB_CLAMPED_LINES (clamped-SIGMA line-rate
consistency probe), DO_SOB_FG_KNOP (Knop formal solver for
SOB+LAM_ES), DO_ES_ALI + ES_ALI_TOL/ES_ALI_MAX (Olson-Auer-Buchler
Jacobi update of the coherent-ES Lambda iterate on the stored
FG_J_CMF_V13 diagonal, inside the USE_LAM_ES J loop).

Live evidence (first 2h interactive segment, job 12412222, warm from
saw_restart_src/adp2): ES_ALI pass counts down freqs unconverged
2464 -> 0 with max final dJ/J reaching the 1D-4 target every great
iteration; run reached it706 at 75% model convergence before a manual
graceful stop (NUM_ITS=0 + READY=1 written directly to
IN_ITS/.CMFGEN_REQUEUE_STATE; the USR1 trap does not fire while
tcsh batch.sh runs in the foreground).

Resubmitted on the small partition for longer segments: job
12413522 (+ OBS 12413523), submission: normal, warm_restart: true,
warm_restart_from removed so the resubmit preserves this run's own
checkpoint instead of re-injecting adp2.  Resumed at it707 with
DO_ES_ALI echoed T.

Maintenance: the submit_cmfgen.sh wall cap did NOT fire at submit
time (no "Maintenance cap:" line in submit.log; timelimit stayed
23:59:00, so Slurm held the job "Reserved for maintenance" -- a 24h
wall does not fit before the 28-Oct 08:00 reservation).  The block
resolves correctly when re-run by hand (maint_start=1793170800, ~14h
to window); likely a transient scontrol failure at submit time since
stderr is suppressed.  Fixed by scontrol update JobId=12413522
TimeLimit=13:30:00 -- Slurm re-evaluates eligibility, USR1 still
fires 7140s before the shortened wall so the graceful checkpoint
works; requeues keep the 13:30:00 limit (bump back to 23:59:00 after
the window if longer chunks are wanted).

## 29-30-Oct-2026 - warmS restart NaN chain isolated and fixed

Diagnosis (KN_LUM_TRACE, f7c75c7c): DEP_RAD_EQ/RLUMST NaN traced to a
0*Inf at a zero-emissivity line resonance in the real sweep
(ZNET_SIM bracket 1 - RJ*CHIL/ETAL with ETAL_MAT=0, SIM=16 K=2 ML=1204)
and the matching dielectronic ZNET/BETAC divisions.  Fix 5600ff97 guards
the singular bracket only; LINE_QW_SUM/JBAR bookkeeping unchanged.

Second failure: DBG-T flagged T(2)=NaN on entry to it866 - a pre-fix
checkpoint had STEQ_T NaN at depths 1-2 which wrote NaN T into
SCRTEMP, so every restart (warm or cold) inherited it.  DBG-T's own
format bug (7 args vs 6 descriptors) had been masking the report;
fixed in 53eb0e88.  KN_REPAIR_T (8ad4d9c2) rebuilds bad T from the
nearest finite depth at iteration entry - T is re-solved the same
pass, so a poisoned checkpoint self-heals.

Open convergence issue (KNOP): lagged-field fixed point stalls on the
compression-shell cells (fg drift oscillates ~5-11, ~48M bad cells,
constant 56377 ALO fallbacks at sigma<-1 depths; exit gate exm
(BA_DEM_EXM) holds at ~60-120 dex).  ES_ALI reports dJ/J=0 because it
only measures the inner solve given the lagged field.  SOB+ESALI on
the same model converged cleanly to 1e-4 (it851).

## 03-Nov-2026 - architecture audit + BA null-pin and HADAPT lifecycle verified

Deliverables from the numerical-architecture review landed as two
tracked docs: `docs/audits/NONMONO_NUMERICAL_ARCHITECTURE_AUDIT_2026-11-03.md`
(commit 0e2940ae) records the solver ordering, active methods, and
classified findings; `docs/NUMERICAL_METHODS_DEVELOPMENT_PLAN.md`
(commit a47e6e9e) phases the follow-on investigations with acceptance
tests.

Two in-flight source changes were verified and committed:

- `CMF_BLKBAND_PIN_NULL_V1` (00f1c727, f235f325, e50a2b48, tests
  736fba2e/495f4420/58b27685): pins only an index whose row, column,
  and RHS are all exactly zero - an independent 0=0 equation whose
  free correction is set to zero rather than failing the whole depth
  block in DGETRF.  Live rows/columns and nonzero RHS keep the
  upstream solve.  Not flag-gated; ordinary blocks are unchanged.
- JFNK_HADAPT lifecycle (b5083bb3, b26176ab, tests a55a2668/
  b415dd7d/54dc9807): the selected FD-step multiplier is now scoped
  to the residual mode and linearization baseline that produced it.
  `JFNK_INIT_V1`, `JFNK_REBUILD_BASELINE_V1`, and a frozen/coupled
  mode change re-arm the h sweep and clear the one-shot diag gate, so
  a stale h cannot carry across a different linearization point.
  `JVP_DIAG_DONE` moved to module state as `JFNK_JVP_DIAG_DONE` for
  the same reason.  All opt-in under `JFNK_HADAPT` (default F).

Verification: `./build.sh` clean (exit 0, zero errors, all three
exes), `cmf_blkband_pin_null_v1_` in `cmfgen_dev.exe`, unit suite
20/20 files pass.  No model run yet - HADAPT JVP-level behaviour is
#status/unverified until a `run_model_grid.py` job exercises it.
Range pushed to origin/main: `e4c8f02b..24f4c970`.

Same-day continuation: replaced the ambiguous T4 LU XFAIL and added the
first JVP cycle-linearity measurement.

- `se_driver` `sepc` now emits per-depth block residuals (`MZDK`, with
  the accepted ridge-ladder index) and the unfactored-depth count
  (`LUBAD`).  Surrogate per-depth blocks reach cond ~1e23..1e39: every
  depth factors (LUBAD=0 on all specs) but the accepted solve is
  meaningless on the near-singular blocks - a real quality gap, not a
  metric artifact.
- Harness split (1c13a4d area, commits 91f59df5..22d63df9): T4 now
  gates the machinery contract (all depths factored, finite z, LSQ
  apply identity) and passes on all four specs; new T21 gates the
  worst per-depth backward error - PASS on std (2.4e-2), kept as named
  XFAIL on mild/hard/shock (3.9e8, 1.0e12, 4.8e13).  ZCORE stays
  printed, ungated.
- New T22 measures JVP cycle-linearity on the surrogate: JV(a v) is
  bitwise equal to a*JV(v) under the production h-rule, and the FD
  additivity defect stays within 7.5x the same-direction FD-vs-exact
  error over eps 1e-8..1e-2 on all specs - direction-dependent h does
  not detach the effective operator beyond the method's own error.
  Production-interface additivity (frozen/coupled residuals) remains
  open; it needs a residual evaluator bound to a real model state.
- Suite 20/20, 89 JFNK checks green.  Docs updated in place (audit
  finding 4 resolved, finding 1 progress note; dev plan Phase-1 marks).

## 09-Nov-2026 - Knop OBS render for SOB populations; FG_NOSLOT_RAY dealloc fix

First numerical A/B of the OBS render modes on the SN_A ERC models
(41-44d CMFGEN snapshots vs the SN2024pba DBSP spectrum at ~35.5d):

- `obs_line_treatment='knop'` added to `prepare_CMF_PARAM.py`:
  `DO_NONMONO_KNOP=T` + `GLOBAL_LINE`/`TRANS_*=BLANK` so the formal J
  runs through `FG_J_CMF_NONMONO_V1` on unclamped `SIGMA_TRUE` - real
  non-monotonic CMF line transfer instead of the deposited
  Doppler-core injection (`'sob'`, `ADD_SOB_LINE_PROF` + second
  `OBS_FRAME_SUB_V9` pass).
- Verified on `SN_A_p2p85_SOB_d44` 43.75d (OBS job 12883823, 51 min vs
  9.5 min deposited): the deposited render's broad 6000-7500 A
  pseudo-continuum hump disappeared, the strongest Ha-region structure
  moved to ~6250 A, and the absolute-flux offset vs DBSP went from
  x0.56 to x1.17.  `SN_A_E15p2p80_SOB_d42` finished identically.
- FG_NOSLOT_RAY STAT=151 fix (same defect class as NM_RAY/DLN_RAY,
  13-Oct): the INTEGRAL-block array was allocated at line ~2190 but
  missing from the NEW_R_GRID dealloc list, so the first
  compression-set rebuild on the amr fields re-ALLOCATEd it and
  `cmf_flux` died (jobs 12886069/70).  One-line dealloc added;
  resubmits 12889834/5 survived repeated rebuilds at depths 42/150.
- Pipeline enforcement: `run_model_grid.py` now forces
  `cmfgen.obs_line_treatment='knop'` under `non_mono_mode=SOB`
  ('sob'/'blank'/unset -> knop with a printed note; 'sob_ew' stays
  selectable since it writes EWDATA, not a spectrum), and
  `prepare_CMF_PARAM.obs_line_treatment()` defaults an unset or
  unreadable key to 'knop'.  Explicit 'sob'/'sob_ew'/'blank' still
  honored at the OBS layer for manual renders.
- Remaining four OBS jobs (d39, d39_amr, e15_amr, d44_amr) queued on
  the small partition (12893194-97) after the 1:30 interactive cap
  killed the amr runs - compression-set rebuild churn, not a solver
  fault.

Verification: `./build.sh` clean (31/31 test files); knop branch
echoed `DO_NONMONO_KNOP=T: CMF formal solution uses
FG_J_CMF_NONMONO_V1` + `DO_NONMONO_FREQ_MATCH=T` auto-set and the
signed-GAM/sigma<0 diagnostics (sigma_min=-0.73, 35 depths);
`obs_fin_10` regenerated for 2/6 models.  #status/partial until the
remaining four spectra land and the grid-level comparison is redone.

## 09-Nov-2026 (2) - COMP_OBS NU_STORE retention fix; p2p85r150 OBS verified

Root cause of the 12893194/95 `CMF_FREQ outside range` STOPs on
`SN_A_p2p85r150_SOB_d39` (+amr): the `COMP_OBS_V2`/`V3` sliding
`NU_STORE`/`IPLUS_STORE` window (NST_CMF=6000/20000) dropped CMF
frequencies only down to `OBS_FREQ`, but rays with
`FREQ_CONV_FAC>1` (transverse rays in the VINF~27k km/s field) need
`CMF_FREQ = OBS_FREQ*FREQ_CONV_FAC` up to ~9% higher.  Once the store
filled, the needed upper entries were discarded and the fail-closed
coverage guard fired (latent upstream defect; the needed frequencies
were solved - the CMF grid tops at MAX_CF=200 - the window just did
not retain them).

Fix (commits 415642cf, 74901201): save
`MAX_FREQ_CONV_FAC=MAXVAL(FREQ_CONV_FAC)` at init beside
`MIN_FREQ_CONV_FAC` and raise the shuffle retention top to
`OBS_FREQ*MAX_FREQ_CONV_FAC` in both routines, symmetric for the
population (V2) and OBS (V3) paths.  The `CMF_FREQ outside range`
STOP stays as backstop for genuine `MAX_CF` undercoverage or
`NST_CMF` overflow.

Regression: `unit_tests/f90/drv_comp_obs_retain.f90` +
`unit_tests/python/test_comp_obs_retain.py` - a ray geometry with
`FREQ_CONV_FAC>1` forces the shuffle path; on the pre-fix source it
reproduces the production signature (`CMF_FREQ=107.97 >
NU_STORE(1)=103.11`), on the fix all 40 observer bins are produced
with max rel flux diff 2.8e-16 vs independent evaluation.
`test_comp_obs.py` 6/6 still passes.

Live verification: interactive job 12901108 on
`SN_A_p2p85r150_SOB_d39` 38.75d - the exact model that STOPped at
17:53 - COMPLETED in 1h16m on the interactive partition (16 CPU,
30G).  `obs_fin_10`/`OBSFLUX` carry all 95408 continuum frequencies,
all-positive fluxes (2e-24..1.7e6 Jy), no NaN/Inf, plus
`hydro_fin_10`/`meanopac_fin`.  The earlier amr OBS pair
(12893196/97) is running on the fixed exe on the small partition.

Caveat: `test_mat5pen.py` fails on a `make -q` archive-member
freshness probe (`libsubs.a(x.o)` reports out-of-date for every
member after the post-build `rm -f *.o`) - a unit_tests
production-object infrastructure issue unrelated to this change;
the MAT5PEN kernel itself was not exercised.  #status/verified for
the triggering model; production resubmits of 12893194/95 (+amr
d39) can proceed on the fixed exe.

## 09-Nov-2026 - JFNK solver audit pass (jfnk_core, production .o linked)

Unit tests in `tests/jfnk_core/` now link `jfnk_solver_v1.o` extracted
verbatim from `lib/libmod.a` plus its production `.mod` (ifort, matching
flags) - previously the suite recompiled from source.  Audit findings
and fixes:

- `JFNK_PRECOND=LSQ` never passed its own acceptance test: `DSYRK`
  half-fills `gn_work`, the `mtst` copy read the uninitialized upper
  triangle.  Fixed (mirror + zero); LSQ now factors with expected
  ~1e-2 Tikhonov bias.
- Schur path: failed per-depth `DPOTRF` factors reached SCHUR_FORM /
  SCHUR_APPLY (DPOTRS cannot flag non-PD input); `JFNK_LU_OK` was set
  but never consulted.  Failed depths now get the identity factor -
  the LSQ solve degrades to a damped gradient step instead of NaN.
- Non-monotone line search read `jfnk_ls_nmhist(1:4)` while the write
  pointer cycled all 8 slots - stale window after wraparound (33/41
  accepts diverged on the test).  Reader now walks back from
  `jfnk_ls_nmidx` (true last-`nmwin_s` accepted norms).
- New opt-in VADAT keys `JFNK_PC_QUAL` (F) / `JFNK_PC_QUAL_TOL`
  (0.5): after a block factorization is accepted, probe
  `||M^-1 A v - v||` on a deterministic direction and feed a failure
  into the ITRY ridge ladder.  Closes the acceptance blind spot where
  neither the `x=1` micro-test nor `rcond(M)` measured `I - M^-1 A`.
- `pivrat` limitation noted in-code: post-DGEEQU equilibration
  self-normalizes the ratio, so `JFNK_RIDGE_PIV` barely triggers.
- S2 documented: rank-deficient Jnn under SCHUR is silently
  pseudo-inverted by lambda damping - acceptable, unflagged.

All jfnk_core suites + `tests/f90/build_and_run.sh` pass.  Vault:
`cmfgen_obsidian/solvers/JFNK solver.md` 09-Nov-2026 section.

## 10-Nov-2026 - Knop I_NEXT cold-start identified as the Luc-mismatch cause; OBS refresh wired

A/B isolation on `m15_2p82d_VAL_sob_jfnk` vs Luc Dessart's reference
`.fl` traced the spectral excess to `FG_J_CMF_NONMONO_V1`'s lagged
I_NEXT coupling, not the populations:

- Same converged population, three renders (band-normed 5500-6500 A):
  05-Nov sob render via stock `FG_J_CMF_V12` matched Luc in every
  continuum band (0.94-1.08) and Ha profile incl. e.s. wings; 10-Nov
  renders via `FG_J_CMF_NONMONO_V1` (both `obs_line_treatment=sob`
  and `knop`) showed a collapsed red/NIR continuum (x0.62-1.48 band
  ratios) and line cores 20-60x above it (Ha peak/cont 59 vs Luc 3).
- Mechanism: the compressing-branch upwind `I_NEXT` is lagged across
  great iterations; the OBS exe has no great-iteration loop, so every
  OBS sweep ran the cold start `I_NEXT=0`.  `chi_hat = chi +
  |GAM|/dLOG_NU` then keeps its |GAM| sink while the compensating
  source term vanishes - a pure photon sink that scales as 1/dLOG_NU
  on the fine OBS grid (~6 km/s bins vs ~35 km/s on the population
  grid) and as 1/CHI (worst in the optical continuum where CHI ~
  ESEC).  All 70 depths sit at sigma < -1, so the whole model is on
  the compressing branch; lines keep their emissivity while the
  continuum and the ESEC*J scattered emission collapse.
- Fix (this edit): `FG_KNOP_REFRESH_V1` is now called from
  `cmf_flux_sub_v5.f` inside the ES loop, after `10000 CONTINUE`
  (complete ML sweep), gated on `ES_COUNTER < NUM_ES_ITERATIONS` -
  it converges the lagged field within the sweep (Paper-II per-ray
  re-solve) and commits via WRITE_STORE to `I_*_DEPTH_NEXT`, which
  INIT leaves intact between sweeps (allocated once under
  FIRST_TIME; restore/zero only on the very first call).  Capture
  arms because the solver runs at every ML (NEW_FREQ=F on
  INACCURATE repeat passes), FREQ_IDX_ITER resets at each sweep's
  ML=1 (FIRST_FREQ->INIT), and `FG_OPT=INT/INS` selects INTEGRAL.
  `MOM_KNOP_REFRESH_V1` mirrored under `DO_NONMONO_SOB` for the
  USE_LAM_ES=F OBS route (inert under the knop recipe, which forces
  USE_LAM_ES=T).  `CMF_KN_RF_MB` env var feeds a new
  `SET_FG_RF_MB` setter - OBS grids ~6x finer than the population
  grid make the ~80 B/cell worker pack exceed the 32 GB
  `KN_RF_TOT_MB` default; unset keeps the default.  NUM_ES>=2 is
  required for any coupling maturity (refresh is skipped on the
  last sweep since nothing reads it); m15 resolved configs now
  carry `NUM_ES: '3'` for two refresh maturations.
- Secondary flags from the comparison (not the cause): OBS render
  VTURB=20 vs population/Luc 50 km/s; extra species (O I, Si I-II,
  Ca I-II, Fe I-II) vs Luc's MODEL_SPEC explain part of the red/NIR
  line forest; `adp2_tauto` warm-start population renders degenerate
  (separate issue).

Verification pending: job 12944101 (small, 150G, 8cpu, NUM_ES=3)
renders the refresh path against Luc; the cold-start knop products
are archived as `obs_fin_10_knop_cold_20261110`.  #status/unverified
until the spectrum lands.

## 10-Nov-2026 (later) - SOB pathway enforcement + OBS dispatch split, committed

Follow-up to the refresh fix above: the SOB pathway was reviewed
against Dessart, Audit & Hillier (2015, arXiv:1503.05463) Sec. 3.1
and the flags it raised were closed in source, pipeline and OBS.

- `DO_NONMONO_SOB` populations converge on the clamped-Sobolev
  approximation: Sobolev line optical depths use unclamped
  SIGMA_TRUE with the |1+sigma*mu^2| floor (bit-identical to
  upstream for sigma >= -1); the continuum formal solve runs
  FG_J_CMF_V13 on clamped SIGMA - upstream V12's own ray-level
  clamp, the same practice Hillier used in the Dessart runs.  This
  is intentionally less accurate than DO_NONMONO_KNOP and serves
  as its initial guess via warm_restart_from.
- `USE_LAM_ES` absent-key default flipped back to T under
  DO_NONMONO_SOB (cb5cce81): the 21-Oct F default routed populations
  to MOM_J_CMF_NONMONO_V1, which failed to contract on m15 (job
  12316397); the V13 route is the validated recipe (job 12001805 /
  VAL_sob_jfnk).  Explicit F stays honoured as an opt-in experiment
  and keeps the DO_KNOP_REFRESH lagged-map refresh.  The probe uses
  GET_KEY_STRING at the Rule-2 block, hundreds of reads past the
  consume at ~3286, because a probe adjacent to the read reports
  absent even when the key is present.
- OBS dispatch split (4b7f628e): COMP_JCONT_V4.INC routes to
  FG_J_CMF_NONMONO_V1 only under DO_NONMONO_KNOP.  SOB-alone OBS
  renders use clamped FG_J_CMF_V12, matching the solver the SOB
  populations were converged against; the 09-Nov routing of SOB
  through Knop produced the UV-heavy SED with collapsed red
  continuum (cold I_NEXT on the fine OBS grid).  The Knop render of
  SOB populations is now the warm-started DO_NONMONO_KNOP run, where
  flag and population solver agree.
- OBS nonmono hooks matched to the same routing (af0ed8f1):
  SET_FG_GS_OBS_GRID and FG_KNOP_REFRESH_V1 are KNOP-gated;
  MOM_KNOP_REFRESH_V1 covers the opt-in SOB+USE_LAM_ES=F route;
  SET_FG_RF_MB (cc6fc9da) lets CMF_KN_RF_MB raise the refresh
  worker budget past the 32 GB default on fine OBS grids.
- Regression: test_obs_nonmono_dispatch.py gained
  test_fg_dispatch_knop_only (a2e67a70) pinning the KNOP-only gate,
  the V12 fall-through and all three OBS hook gates.  All checks
  pass; ./build.sh clean; cmf_flux.exe, obs_frame.exe and
  cmfgen_dev.exe rebuilt.
- Pipeline mirrors: prepare_VADAT.py stamps USE_LAM_ES=T when the
  key is absent under SOB; prepare_CMF_PARAM.py mirrors the executed
  ../VADAT for line treatment and propagates COH_ES/NUM_ES/
  ES_WING_EXT/R_CMF_WING_EXT/THK_CONT/IB_METH/MIN_CF/MAX_CF into the
  OBS template; run_model_grid.py keeps the warm_restart_from
  plumbing for the SOB->KNOP handoff.  (scripts/ tree, not this repo)
- Docs: cmfgen_jfnk/AGENTS.md section 3 and section 4.3 updated;
  obsidian notes solvers/Sobolev non-mono solver.md and
  reference/VADAT keys.md updated (both gitignored, local only).

Commits pushed to origin/main: cc6fc9da, 4b7f628e, af0ed8f1,
cb5cce81, a2e67a70 (HEAD a2e67a70).  cmfgen_build for submitted
runs should pin this HEAD or later.

## 10-Nov-2026 (3) - FG field-drift audit defects A3, A1, C1 committed

First three confirmed defects from docs/audits/FG_FIELD_DRIFT_AUDIT.md
(the KNOP refresh drift investigation) landed in source, verified by
unit_tests/python/test_knop_ray.py (22/22) and test_knop_lagged_field.py
plus a full ./build.sh (all exes, all suites green; test_par_opac.py
fails only because the root make does not build opac/ - pre-existing,
passes 8/8 once the objects exist).

- A3 (a0366958): KN_RES_PASS updates now sit below the divergence check
  in both SWEEP_FLD and SWEEP_FLD_ORD, so a diverged or non-finite pass
  cannot poison the residual history the adaptive relaxation reads.
- A1 (a0366958): the maximum-envelope contraction estimate is replaced
  by a per-ray, same-direction, same-map estimator. Each field solve
  keeps a six-slot ratio ring (post-transient only), exports its median
  as RHO_M/RHO_P, the refresh pools the per-ray medians and drives
  INEXT_RELAX from the P90 with 50/50 smoothing and the [0.05,1] clamp;
  the old cross-refresh ratio stays as the low-sample fallback under
  FG_ADAPT_OMEGA_V2.
- C1 (a0366958, comment 6e82257c): FG_RF_CAPPED was all-pass cap/floor
  event counts plus the 0/1 diverged-ray flag divided by the last-pass
  cell census - a grinding ray reported ~NPASS times its true share.
  The numerator is now last-pass events per direction plus commit-phase
  caps (NCP_P/NNG_P per-pass counters, NCP_C in the store writers,
  NCL_LS/NNL_LS per ray). FG_INEXT_CAPPED keeps the all-pass
  cumulative for checkpoint comparability; diverged rays still tier
  separately via GS_NDIV. Not the active m15 drift cause (clamp
  events=0 live) but a latched false-UNSAFE risk removed.
- Tests (86069f0d, driver 542263e0): T16-T18 pin the contraction
  export and pooled estimator; T19 (capx, NONMONO_INEXT_CAP=1e-4)
  pins the census bound - NCAP=168 all-pass vs NCAP_LST=84 same-census
  under NTOT=90.
- Docs: docs/audits/FG_FIELD_DRIFT_AUDIT.md committed (f833919b);
  obsidian Knop formal solver.md and Diagnostics and OUTGEN.md updated
  (gitignored, local only).

Commits pushed to origin/main: a0366958, 542263e0, 86069f0d,
6e82257c, f833919b (HEAD f833919b).

## 11-Nov-2026 (2) - ES_ALI/BA audit landed; DO_ES_ALI auto-set under SOB; SOB over KNOP directive

The ES_ALI / ES_ALI_BA / BA_RES_ACC stack audit
(docs/audits/ES_ALI_BA_AUDIT.md) and the warm-arm experiment series
(docs/audits/SOB_ACCELERATION_RESULTS.md, X1-X4 plus cold X5)
completed.  Verdict: CASE A - DO_ES_ALI=T alone collapses the
ES-scattering-dominated tail (X2 converged ~30 new iterations from
the FRESH_sob checkpoint where unaccelerated X1 was still pinned);
DO_ES_ALI_BA added nothing (X3 == X2) and stays opt-in;
DO_BA_RES_ACC acted as a rejection/trust safety, not a speedup
(CASE C).  Slow-mode localisation in docs/audits/
SOB_SLOW_MODE_ANALYSIS.md; JFNK-CPL kept as a proposal only
(docs/audits/JFNK_CPL_PROPOSAL.md) since the continuum stack
sufficed on this model.

Guard work landed with the audit: CHECK_PARAM_CONSISTENCY now
hard-stops DO_ES_ALI_BA=T without DO_ES_ALI=T and warns when
DO_ES_ALI=T cannot engage (COH_ES off, no nonmono flag, or
DO_SOB_FG_KNOP); ES_ALI_NFLOOR / ES_ALI_BA_NFLOOR count denominator-
floor activations in the J loop and the BA fill (the floors differ
by design: 0.01 vs ES_ALI_BA_MIN=0.1, neither engaged on m15 - X3
showed 50 BA-floor hits with zero convergence change).

Policy change (user directive): DO_ES_ALI joins the SOB absent-key
auto-sets (T when COHERENT_ES and not DO_SOB_FG_KNOP, both USE_LAM_ES
routes carry a stored diagonal); explicit VADAT F still opts out.
New configs use DO_NONMONO_SOB + USE_LAM_ES=T; DO_NONMONO_KNOP is
retained for reference/diagnostic comparisons only.  The 'knop' OBS
line render produced negative OBSFLUX luminosities on the m15 SOB
models, so SOB-population configs now select obs_line_treatment:
blank (pipeline configs, not this repo).

Docs: AGENTS.md section 3 updated (auto-set list, solver
preference, DO_ES_ALI evidence); obsidian VADAT keys.md updated
(gitignored, local only).

## 11-Nov-2026 (3) - KNOP production-readiness audit: NOT READY

docs/audits/KNOP_PRODUCTION_AUDIT.md synthesises the full KNOP
audit chain (KNOP_AUDIT, KNOP_FORMAL_SOLVER_AUDIT,
KNOP_STALL_DIAGNOSIS, FIELD_DRIFT/FG_FIELD_DRIFT,
KNOP_ITEM1_RUNAWAY, ITEM2/POPULATION_SOLVER_FIXES,
NONMONO_COMPREHENSIVE, NONMONO_NUMERICAL_ARCHITECTURE,
CONVENTIONAL_SOLVER, FULL_ITERATION_CYCLE, JFNK_PRODUCTION_OPERATOR,
CMF_FLUX_KNOP_CONVERGENCE) plus live run evidence.

Verdict: NOT PRODUCTION READY.  No KNOP run has ever converged a
model end to end.  Live evidence this session: FRESH_knop
(job 13005611) STOPped non-convergent at it124 on
GS_REL_NOPROG_LIM (27 releases, demand oscillating 38-52, lagged
field never matured); VAL_knop_formal STOPped it1489 on luminosity
divergence (L_out=1.28E15, maxraw pinned at 1E7 cap, KN_NO_GAIN
flooring 8447 cells); FRESH_knop_d1a still running at ~50-75%;
VAL_knop_mom shows repeated formal-solve convergence failures;
KNOP_opt/VAL_knop_jfnk died at setup on a GENOSC_V9 data error.

Two kill mechanisms: field-map non-maturity (GS freeze/release
cycle stalls - the watchdog catches it but cannot cure it) and
the population-side runaway/null-space mode (H1 theorised,
POP_RUNAWAY census instrumented, remedies open).  OBS path has
verified defects D1-D4 (cold I_NEXT always, ~9.5x continuum
suppression) with fixes C1-C5 unapplied, and S6 is blocked on
blocker 1 anyway.  Accelerators are safe on KNOP (ES_ALI J-loop
clean on the live run) but unproven as a convergence fix - the
KNOP_esA/esB A/B never ran.

User directive recorded: DO_NONMONO_SOB + USE_LAM_ES=T is the
production route; KNOP is reference/diagnostic only.  Path-to-
production items ranked in the audit doc section 8.
